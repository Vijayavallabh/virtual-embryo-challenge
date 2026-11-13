"""Tests for tries K11 (experiments/t3_decomp/k11_factor_arms.py), on synthetic data only: no competition file is
opened, nothing is scored, and no prediction is compared with any target."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import eb_factor as eb  # noqa: E402
import k9_residual as k9  # noqa: E402
import k11_factor_arms as k11  # noqa: E402
import test_k9_residual as t9  # noqa: E402  (K9's synthetic world, reference and spy, reused)

SOURCE = HERE / "k11_factor_arms.py"
G = t9.G
CLIP = t9.CLIP
world = t9.world  # the module fixture: synthetic files in task x1's roles


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


# ----------------------------------------------------------------------------------------------- constants and tasks
def test_reused_files_are_pinned_as_they_were_before_k11():
    want = {"k9_residual.py": "1c58a0908b20cbb196bd3d1f64a522a7858e1d74b73fb1f73bf743c893af1f1c",
            "eb_factor.py": "36ec4c2daeff406c24731ca3814fda57a339f6192a66f08125d4d0dfbd2da9d8"}
    for name, sha in want.items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == sha, name


def test_the_fixed_constants_arms_and_blocked_lines():
    assert k11.DID_C == 1.0 and k11.C3_STRENGTH == 1.0
    assert k11.PERM_SEED == k9.PERM_SEED == 20261113 and k11.STREAM == k9.STREAM == harness.K4_STREAM
    assert (k11.MIN_GAIN, k11.MIN_AHEAD, k11.X2_FLOOR) == (0.3, 8, -0.3)
    assert k11.ARMS == ("raw", "raw_did", "raw_did_perm", "raw_c3", "raw_c3_perm") and k11.BASE == k9.BASE
    assert k11.CONTROLS == {"raw_did": "raw_did_perm", "raw_c3": "raw_c3_perm"}
    assert k11.NOT_BUILT_KEYS == ("source_mesp1_labels_in_contrast", "source_nonmesp1_labels")
    assert set(k11.BLOCKED) == {"raw_lvc", "raw_r8"}
    assert all(r.startswith("design not on record") for r in k11.BLOCKED.values())
    assert not set(k11.BLOCKED) & set(k11.ARMS)


def test_the_tasks_are_k9s_own():
    assert k11.TASKS is k9.TASKS and k11.KEYS == k9.KEYS
    assert all(c.applies(t) for c in k11.candidates()[1:] for t in k11.TASKS.values())
    assert not any(c.applies(t) for c in k11.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)
    assert [c.name for c in k11.candidates()] == [harness.FLOOR, *k11.ARMS]


def test_the_source_opens_no_file_beyond_the_tasks_inputs():
    src = SOURCE.read_text()
    src = src.replace(k11.__doc__, "")  # the docstring cites the report experiments/research_2026_11_runs/REPORT.md
    assert "research_2026_11_runs" not in src
    for call in ("read_h5ad", "load_stage", "read_stage", "h5py", "np.load", "open(", "boards/", "runs/", "data/"):
        assert call not in src, call
    assert re.findall(r'"[^"]*\.h5ad"', src) == []


# ----------------------------------------------------------------------------------------------- raw_did
def _mu(rng, L=5):
    mu_wt = rng.dirichlet(np.ones(G), size=L)
    w = rng.dirichlet(np.ones(L))
    return mu_wt, w


def test_the_did_contrast_is_l_minus_n_of_gene_factors_log_ratio():
    rng = np.random.default_rng(0)
    mu_wt, w = _mu(rng)
    mu_ko = mu_wt * np.exp(rng.normal(0, 0.5, mu_wt.shape))
    ml, mn = np.array([1, 1, 0, 0, 0], bool), np.array([0, 0, 1, 1, 0], bool)
    eps = 1e-5
    D = k11.did_contrast(mu_ko, mu_wt, w, ml, mn, eps)

    def lr(mask):
        ww = (w * mask)[:, None]
        return np.log(((ww * mu_ko).sum(0) + eps) / ((ww * mu_wt).sum(0) + eps))

    assert np.allclose(D, lr(ml) - lr(mn))
    # gene_factor's own log ratio before its rescale and clip, on N
    t3d = harness._t3_head().t3d
    phi, den = t3d.gene_factor(mu_ko, mu_wt, w, mn, eps, (1e-9, 1e9))
    assert np.allclose(np.log(phi) - np.log(phi).mean(), lr(mn) - lr(mn).mean())
    # a knockout effect shared by L and N cancels; nothing in L: no contrast
    shared = mu_wt * np.exp(rng.normal(0, 0.5, G))[None, :]
    assert np.allclose(k11.did_contrast(shared, mu_wt, w, ml, mn, 0.0), 0.0)
    assert np.array_equal(k11.did_contrast(mu_ko, mu_wt, w, np.zeros(5, bool), mn, eps), np.zeros(G))


def test_the_did_multiplier_and_its_strength():
    D = np.random.default_rng(1).normal(0, 0.4, G)
    assert np.allclose(k11.did_factors(D)["mesp1"], np.exp(D))
    assert np.all(k11.did_factors(D, 0.0)["mesp1"] == 1.0)
    assert set(k11.did_factors(D)) == {"mesp1"}


def test_lineage_classes_follow_the_lineage_maps_mesp1_flag():
    V = SimpleNamespace(mesp1=np.array(["yes", "partial", "no", "yes"]))
    lab = np.array([0, 1, 2, 3, 3, 1])
    assert list(k11.lineage_classes(V, lab)) == ["mesp1", "other", "other", "mesp1", "mesp1", "other"]


# ----------------------------------------------------------------------------------------------- raw_c3
def _sections(rng, K=6, L=3, n=60, outlier=None):
    lab = np.tile(np.arange(L), K * n // L)[: K * n]
    block = np.repeat(np.arange(K), n)
    base = rng.dirichlet(np.ones(G), size=L)
    C = np.vstack([rng.multinomial(200, base[lab[i]]) for i in range(K * n)]).astype(np.int32)
    if outlier is not None:
        C[block == outlier, 0] += 400  # one section with an artefact on gene 0
    mu_wt = rng.dirichlet(np.ones(G), size=L)
    return C, lab, block, mu_wt, np.full(L, 1.0 / L), np.ones(L, bool)


def test_section_estimates_are_eb_factors_jackknife():
    rng = np.random.default_rng(2)
    C, lab, block, mu_wt, w, mask = _sections(rng)
    y, ys, den, K = k11.section_estimates(C, lab, block, mu_wt, w, mask, 1e-5)
    jk = eb.jackknife(C, lab, block, mu_wt, w, mask, 1e-5)
    assert K == jk.n_blocks == 6
    assert np.array_equal(y, jk.y) and np.array_equal(den, jk.den)
    assert np.allclose((K - 1) / K * ((ys - ys.mean(0)) ** 2).sum(0), jk.v)
    with pytest.raises(ValueError, match="two sections"):
        k11.section_estimates(C, lab, np.zeros_like(block), mu_wt, w, mask, 1e-5)


def test_the_robust_factor_resists_one_artefact_section():
    rng = np.random.default_rng(3)
    clean = k11.section_estimates(*_sections(np.random.default_rng(3)), 1e-5)
    bad = k11.section_estimates(*_sections(rng, outlier=2), 1e-5)
    y_bad, rob_bad = bad[0][0], k11.robust_log_ratio(bad[0], bad[1])[0]
    assert abs(rob_bad - clean[0][0]) < abs(y_bad - clean[0][0])
    # identical sections: nothing to do, the robust value is the plain estimate
    one = _sections(np.random.default_rng(4), K=1, n=60)
    C = np.vstack([one[0]] * 5)
    lab, block = np.tile(one[1], 5), np.repeat(np.arange(5), 60)
    y, ys, _, _ = k11.section_estimates(C, lab, block, one[3], one[4], one[5], 1e-5)
    assert np.allclose(k11.robust_log_ratio(y, ys), y, atol=1e-10)


def test_c3_factor_at_strength_zero_is_the_raw_factor_and_otherwise_finalized():
    rng = np.random.default_rng(5)
    C, lab, block, mu_wt, w, mask = _sections(rng, outlier=1)
    t3d = harness._t3_head().t3d
    phi, den = t3d.gene_factor(t3d.props_type_means(C, lab, 3), mu_wt, w, mask, 1e-5, CLIP)
    frame = SimpleNamespace(phi=phi, C=C, lab=lab, block=block, mu_wt=mu_wt, weights=w, mask=mask, eps=1e-5,
                            clip=CLIP)
    assert np.array_equal(k11.c3_factor(frame, 0.0), phi)
    y, ys, d, _ = k11.section_estimates(C, lab, block, mu_wt, w, mask, 1e-5)
    got = k11.c3_factor(frame)
    assert np.array_equal(got, eb.finalize(k11.robust_log_ratio(y, ys), d, CLIP))
    assert got.min() >= CLIP[0] and got.max() <= CLIP[1] and not np.allclose(got, phi)
    assert np.allclose(eb.finalize(y, d, CLIP), phi)  # the plain estimate through finalize is the raw factor


def _c3_frame(outlier=1, seed=5):
    """A synthetic source frame whose raw factor is finalize of its own plain estimate y, so that a factor built as
    finalize(y + 0) is raw's factor bit for bit (gene_factor's ratio and exp of its log can differ in the last bit)."""
    C, lab, block, mu_wt, w, mask = _sections(np.random.default_rng(seed), outlier=outlier)
    y, _, den, _ = k11.section_estimates(C, lab, block, mu_wt, w, mask, 1e-5)
    return SimpleNamespace(phi=eb.finalize(y, den, CLIP), C=C, lab=lab, block=block, mu_wt=mu_wt, weights=w,
                           mask=mask, eps=1e-5, clip=CLIP)


def test_the_c3_control_permutes_the_shift_within_deciles_of_the_basis():
    frame = _c3_frame()
    basis = np.random.default_rng(11).gamma(2.0, size=G)
    y, ys, den, _ = k11.section_estimates(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask,
                                          frame.eps)
    shift = k11.robust_log_ratio(y, ys) - y
    perm = harness.t3_decile_permutation(shift, basis, k9.PERM_SEED)
    got = k11.c3_perm_factor(frame, basis)
    assert np.array_equal(got, eb.finalize(y + perm, den, CLIP))
    assert not np.allclose(got, k11.c3_factor(frame)) and not np.allclose(perm, shift)
    # the machinery raw_did_perm uses: the same permutation of the same vector
    assert np.array_equal(perm, k9.permuted_residuals({"g": shift}, basis, k11.PERM_SEED)["g"])
    for genes in np.array_split(np.argsort(basis, kind="stable"), 10):  # each decile keeps its values
        assert sorted(perm[genes]) == sorted(shift[genes])
    assert np.array_equal(k11.c3_perm_factor(frame, basis, strength=0.0), frame.phi)


def test_the_c3_control_is_raw_bit_for_bit_when_the_robust_value_is_the_plain_one(bs, monkeypatch):
    frame = _c3_frame()
    monkeypatch.setattr(k11, "robust_log_ratio", lambda y, ys: np.array(y, dtype=np.float64))  # m == y
    basis = np.random.default_rng(12).gamma(2.0, size=G)
    phi_perm = k11.c3_perm_factor(frame, basis)
    assert np.array_equal(phi_perm, frame.phi) and np.array_equal(k11.c3_factor(frame), frame.phi)
    rng = np.random.default_rng(13)
    ref = t9._reference(bs, rng)
    prep = dict(_prep(rng, ref), phi_raw=np.array(frame.phi), phi_c3=k11.c3_factor(frame), phi_c3_perm=phi_perm)
    raw = k11.arm_built(bs, ref, 200, 14, ref.V.p, "raw", prep)
    for arm in ("raw_c3", "raw_c3_perm"):
        b = k11.arm_built(bs, ref, 200, 14, ref.V.p, arm, prep)
        assert np.array_equal(b.X, raw.X) and np.array_equal(b.row, raw.row) and np.array_equal(b.lab, raw.lab), arm


# ----------------------------------------------------------------------------------------------- the arms
def _prep(rng, ref, D_scale=0.8):
    n_pool = ref.pool.C.shape[0]
    phi_raw = np.clip(np.exp(rng.normal(0, 0.6, G)), *CLIP)
    rho = k11.did_factors(rng.normal(0, D_scale, G))
    out = {"phi_raw": phi_raw, "rho_did": rho, "rho_did_perm": k9.permuted_residuals(rho, rng.gamma(2.0, size=G)),
           "phi_c3": np.clip(phi_raw * np.exp(rng.normal(0, 0.3, G)), *CLIP), "clip": CLIP,
           "classes": np.array(["mesp1", "other"])[rng.integers(0, 2, n_pool)]}
    out["phi_c3_perm"] = np.clip(phi_raw * np.exp(rng.normal(0, 0.3, G)), *CLIP)  # drawn last: the pins hold
    return out


def test_every_arm_draws_the_same_cells_and_did_moves_only_mesp1_cells_factor(bs, monkeypatch):
    rng = np.random.default_rng(6)
    ref = t9._reference(bs, rng)
    prep = _prep(rng, ref)
    calls = t9._spy(bs, monkeypatch)
    built = {arm: k11.arm_built(bs, ref, 200, 7, ref.V.p, arm, prep) for arm in k11.ARMS}
    rows, lab = built["raw"].row, built["raw"].lab
    assert all(np.array_equal(b.row, rows) and np.array_equal(b.lab, lab) for b in built.values())
    assert all(b.X.shape == built["raw"].X.shape for b in built.values())
    for a in k11.ARMS[1:]:
        assert not np.array_equal(built[a].X, built["raw"].X), a
    want = k9.arm_built(bs, ref, 200, 7, ref.V.p, k9.BASE, prep["phi_raw"], None, None, CLIP)
    assert np.array_equal(built["raw"].X, want.X)  # K9's raw arm exactly
    # raw: 1 call; raw_did and its control: one call per class; raw_c3 and its control: 1 call each; K9's raw: 1 call
    assert len(calls) == 1 + 2 + 2 + 1 + 1 + 1
    assert np.array_equal(calls[0][1], prep["phi_raw"]) and np.array_equal(calls[7][1], prep["phi_raw"])
    assert np.array_equal(calls[6][1], prep["phi_c3_perm"]) and np.array_equal(calls[6][0], calls[5][0])
    drawn = prep["classes"][rows]
    for k, key in ((1, "rho_did"), (3, "rho_did_perm")):
        (Cm, phim), (Co, phio) = calls[k], calls[k + 1]
        assert np.array_equal(Cm, ref.pool.C[rows][drawn == "mesp1"])
        assert np.array_equal(phim, np.clip(prep["phi_raw"] * prep[key]["mesp1"], *CLIP))
        assert np.array_equal(Co, ref.pool.C[rows][drawn == "other"]) and np.array_equal(phio, prep["phi_raw"])
    assert np.array_equal(calls[5][1], prep["phi_c3"])


def test_with_the_idea_off_every_arm_is_raw_bit_for_bit(bs):
    rng = np.random.default_rng(8)
    ref = t9._reference(bs, rng)
    prep = _prep(rng, ref)
    raw = k11.arm_built(bs, ref, 200, 9, ref.V.p, "raw", prep)
    off = dict(prep, rho_did=k11.did_factors(np.ones(G), 0.0), rho_did_perm=k11.did_factors(np.zeros(G)),
               phi_c3=prep["phi_raw"].copy(), phi_c3_perm=prep["phi_raw"].copy())
    for arm in k11.ARMS[1:]:
        b = k11.arm_built(bs, ref, 200, 9, ref.V.p, arm, off)
        assert np.array_equal(b.X, raw.X) and np.array_equal(b.row, raw.row) and np.array_equal(b.lab, raw.lab), arm
    # every cell outside the Mesp1 lineage: nothing for raw_did to do
    only_other = dict(prep, classes=np.full(ref.pool.C.shape[0], "other"))
    for arm in ("raw_did", "raw_did_perm"):
        assert np.array_equal(k11.arm_built(bs, ref, 200, 9, ref.V.p, arm, only_other).X, raw.X)


def _sha(*arrays) -> str:
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str((a.dtype.str, a.shape)).encode())
        h.update(a.tobytes())
    return h.hexdigest()


PINNED_UNIT = {  # taken at a21b484, before the raw_c3_perm fix
    "raw": "dd7ed70296c4a76943854a86b9e62b31cdb3c880bd7347de0b568594181163b9",
    "raw_did": "1f99a9caf45ba15ff340facf03eb9c164953c7159f51086abefc088de5efb4f3",
    "raw_did_perm": "6c7f10c21f75e8af74df67f256a22634c4c9494956e9fdd3bff9e66648d68ffc",
    "raw_c3": "2dcb46e460c538e717588e205d6f6643cb872bd6de9204e4f21eb43ded9ccc34"}


def test_the_existing_arms_are_pinned_to_a21b484_on_the_synthetic_cases(bs):
    """raw, raw_did, raw_did_perm and raw_c3 built as at a21b484 (checksums taken there, before the raw_c3_perm fix)."""
    rng = np.random.default_rng(6)
    ref = t9._reference(bs, rng)
    prep = _prep(rng, ref)
    got = {}
    for arm in ("raw", "raw_did", "raw_did_perm", "raw_c3"):
        b = k11.arm_built(bs, ref, 200, 7, ref.V.p, arm, prep)
        got[arm] = _sha(b.X, b.row, b.lab)
    assert got == PINNED_UNIT


PINNED_E2E = {  # taken at a21b484, before the raw_c3_perm fix
    "raw": "ab8b7e9adcda777150ed44be5767f94e1256f1c005fc653089c13e93afccb339",
    "raw_did": "87fceddbd0cf80d6275b840559fc9792ec85a106468ab665242586ca3a01e505",
    "raw_did_perm": "f2a89acfb41a1253264bd5ee963c364b2c7a0ae193df969087dbc0eec507bceb",
    "raw_c3": "44d2a7f66083a6c1b77274af9410212df75552cada311c026fa625c8719ed0ca",
    "factors": "e6d5b925cbd0ac7739fb2542ae4f11d3cd972fbea5df19706bcd9d9b9e6a3823",
    "record": ('{"cosine_log_phi_raw_did": -0.19853520182487364, "cosine_log_phi_raw_log_phi_c3": 0.9957836546265669, '
               '"fold": "T3/Mab21l2->Gata4rep2", "reference": "E8.75", "sections": 20, '
               '"source_mesp1_labels_in_contrast": 3, "source_nonmesp1_labels": 3}')}


def test_end_to_end_the_existing_arms_are_pinned_to_a21b484(world, monkeypatch):
    t = k11.TASKS["x1"]
    k11._PREP.clear()
    ctx = t9._ctx(t, world)
    with harness.guarded_reads():
        got = {}
        for arm in ("raw", "raw_did", "raw_did_perm", "raw_c3"):
            o = k11._build(ctx, arm)
            got[arm] = _sha(o.X, o.coords)
    prep = k11._PREP["fold"]
    got["factors"] = _sha(prep["phi_raw"], prep["rho_did"]["mesp1"], prep["rho_did_perm"]["mesp1"], prep["phi_c3"])
    got["record"] = json.dumps(prep["record"], sort_keys=True)
    k11._PREP.clear()
    assert got == PINNED_E2E


# ----------------------------------------------------------------------------------------------- end to end, synthetic
def test_end_to_end_raw_is_k9s_raw_and_the_arms_keep_its_cells(world, monkeypatch):
    t = k11.TASKS["x1"]
    k11._PREP.clear()
    k9._PREP.clear()
    ctx = t9._ctx(t, world)
    monkeypatch.setattr(t9.k7, "MIN_CELLS", 20)
    with harness.guarded_reads() as g:
        out = {arm: k11._build(ctx, arm) for arm in k11.ARMS}
        k9raw = k9._build(ctx, k9.BASE)
    assert not g["refused"]
    prep = k11._PREP["fold"]
    assert np.array_equal(prep["phi_raw"], k9._PREP["fold"]["phi_raw"])
    assert np.array_equal(out["raw"].X, k9raw.X) and np.array_equal(out["raw"].coords, k9raw.coords)
    for a in k11.ARMS[1:]:
        assert out[a].X.shape == out["raw"].X.shape and np.array_equal(out[a].coords, out["raw"].coords), a
    rec = prep["record"]
    assert rec["source_mesp1_labels_in_contrast"] >= 1 and rec["source_nonmesp1_labels"] >= 1
    assert json.loads(json.dumps(rec)) == rec
    assert (prep["classes"] == "mesp1").any() and (prep["classes"] == "other").any()
    assert not np.array_equal(out["raw_did"].X, out["raw"].X)
    assert not np.array_equal(out["raw_c3"].X, out["raw"].X)
    assert not np.array_equal(out["raw_c3_perm"].X, out["raw"].X)
    assert not np.array_equal(out["raw_c3_perm"].X, out["raw_c3"].X)
    fold = harness.K4_FOLDS[k11.FOLD_OF[t.key]]
    frame = harness.t3_source_frame(ctx.inputs["wt_train"], [ctx.inputs[r] for r in fold["ko_roles"]],
                                    harness._t3_head().load_rules(), fold["source_preset"])
    basis = harness.t3_mean_expression(ctx.inputs[ctx.task.reference_role])
    assert np.array_equal(prep["phi_c3_perm"], k11.c3_perm_factor(frame, basis, k9.PERM_SEED))
    assert k11.did_built(rec) and k11.not_built(rec) == []
    # the ideas switched off: bit for bit raw
    k11._PREP.clear()
    monkeypatch.setattr(k11, "DID_C", 0.0)
    monkeypatch.setattr(k11, "C3_STRENGTH", 0.0)
    with harness.guarded_reads():
        for a in (k11.DID, k11.C3, k11.C3_PERM):
            off = k11._build(ctx, a)
            assert np.array_equal(off.X, out["raw"].X) and np.array_equal(off.coords, out["raw"].coords), a
    k11._PREP.clear()
    k9._PREP.clear()


# ----------------------------------------------------------------------------------------------- the reading
def _result(key, scores, clean=True, commit="abc", record=None):
    rows = [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5, "severity": 0.6} for a in scores},
             "t3": {a: {"severity": {"state": "scored"}} for a in scores}} for i, s in enumerate(SEEDS)]
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
            "factor_record": record if record is not None else {"source_mesp1_labels_in_contrast": 3,
                                                                 "source_nonmesp1_labels": 3}}


def _scores(did=54.0, perm=50.0, c3=54.0, c3perm=50.0, raw=52.0):
    n = len(SEEDS)
    return {harness.FLOOR: [50.0] * n, "raw": [raw] * n, "raw_did": [did] * n, "raw_did_perm": [perm] * n,
            "raw_c3": [c3] * n, "raw_c3_perm": [c3perm] * n}


def _read(d: Path, k9_dir=None, records=None, **per_task):
    d.mkdir()
    for name, sc in per_task.items():
        (d / f"{name}.json").write_text(json.dumps(_result(k11.TASKS[name].key, sc,
                                                           record=(records or {}).get(name))))
    return k11.reading(d, k9_dir if k9_dir is not None else d / "no_k9")


def test_the_reading_per_arm(tmp_path):
    yes = _read(tmp_path / "yes", x1=_scores(), r12=_scores(), x2=_scores(did=51.8, c3=51.7))
    assert not yes["void"] and yes["arms"]["raw_did"]["supported"] and yes["arms"]["raw_c3"]["supported"]
    assert "K10" in yes["arms"]["raw_did"]["verdict"] and "first-choice" in yes["arms"]["raw_did"]["verdict"]
    assert yes["k9_reproduction"].startswith("not checked")
    assert set(yes["blocked"]) == {"raw_lvc", "raw_r8"}
    txt = k11.text(yes)
    assert "NOT the adoption rule" in txt and "never counts" in txt and "Blocked" in txt and "K10" in txt
    assert all(v["readable"] and v["not_built_on"] == [] for v in yes["arms"].values())
    # each arm's own control binds it, on x1 and on r12
    ctrl = _read(tmp_path / "ctrl", x1=_scores(perm=53.9), r12=_scores(), x2=_scores())
    assert not ctrl["arms"]["raw_did"]["supported"] and ctrl["arms"]["raw_c3"]["supported"]
    c3c = _read(tmp_path / "c3c", x1=_scores(), r12=_scores(c3perm=53.9), x2=_scores())
    assert c3c["arms"]["raw_did"]["supported"] and not c3c["arms"]["raw_c3"]["supported"]
    assert c3c["arms"]["raw_c3"]["verdict"].startswith("CLOSED")
    assert not c3c["arms"]["raw_c3"]["checks"]["r12: raw_c3 - raw_c3_perm"]
    assert c3c["tasks"]["r12"]["paired"]["raw_c3_perm - raw"]["mean"] == pytest.approx(1.9)
    c3x = _read(tmp_path / "c3x", x1=_scores(c3perm=53.9), r12=_scores(), x2=_scores())
    assert not c3x["arms"]["raw_c3"]["supported"]
    # x2 below -0.3 closes; r12 failing closes; x1 failing closes
    assert not _read(tmp_path / "x2", x1=_scores(), r12=_scores(), x2=_scores(c3=51.6))["arms"]["raw_c3"]["supported"]
    r12 = _read(tmp_path / "r12", x1=_scores(), r12=_scores(did=52.2, c3=52.2), x2=_scores())
    assert not any(v["supported"] for v in r12["arms"].values())
    assert r12["arms"]["raw_did"]["verdict"].startswith("CLOSED")
    x1 = _read(tmp_path / "x1", x1=_scores(c3=52.0), r12=_scores(), x2=_scores(c3=60.0))
    assert not x1["arms"]["raw_c3"]["supported"]


def test_did_built_needs_both_label_sets_in_the_factor_record():
    assert k11.did_built({"source_mesp1_labels_in_contrast": 2, "source_nonmesp1_labels": 1})
    for rec in ({"source_mesp1_labels_in_contrast": 0, "source_nonmesp1_labels": 4},
                {"source_mesp1_labels_in_contrast": 3, "source_nonmesp1_labels": 0},
                {"source_mesp1_labels_in_contrast": 3}, {}):
        assert not k11.did_built(rec) and k11.not_built(rec) == ["raw_did", "raw_did_perm"], rec


def test_did_contrast_without_mesp1_labels_is_zero_so_raw_did_would_be_raw(bs):
    rng = np.random.default_rng(15)
    mu_wt, w = _mu(rng)
    mu_ko = mu_wt * np.exp(rng.normal(0, 0.5, mu_wt.shape))
    D = k11.did_contrast(mu_ko, mu_wt, w, np.zeros(5, bool), np.array([0, 0, 1, 1, 0], bool), 1e-5)
    ref = t9._reference(bs, rng)
    prep = _prep(rng, ref)
    prep = dict(prep, rho_did=k11.did_factors(D))
    raw = k11.arm_built(bs, ref, 200, 16, ref.V.p, "raw", prep)
    assert np.array_equal(k11.arm_built(bs, ref, 200, 16, ref.V.p, "raw_did", prep).X, raw.X)


@pytest.mark.parametrize("fold,record", [
    ("x1", {"source_mesp1_labels_in_contrast": 0, "source_nonmesp1_labels": 3}),
    ("r12", {"source_mesp1_labels_in_contrast": 2, "source_nonmesp1_labels": 0})])
def test_raw_did_not_built_on_x1_or_r12_is_not_readable_not_closed(tmp_path, fold, record):
    # raw_did equals raw there (as did_contrast's zeros make it), which would otherwise read as CLOSED
    sc = {n: _scores() for n in k11.TASKS}
    sc[fold] = _scores(did=52.0, perm=52.0)
    res = _read(tmp_path / fold, records={fold: record}, **sc)
    did = res["arms"]["raw_did"]
    assert not res["void"] and not did["supported"] and not did["readable"] and did["not_built_on"] == [fold]
    assert did["verdict"].startswith("NOT READABLE") and not did["verdict"].startswith("CLOSED")
    assert res["tasks"][fold]["not_built"] == ["raw_did", "raw_did_perm"]
    assert res["arms"]["raw_c3"]["readable"] and res["arms"]["raw_c3"]["supported"]  # raw_c3 is read as before
    txt = k11.text(res)
    assert "NOT READABLE" in txt and "not built on this fold" in txt
    # passing checks do not make it readable either
    res2 = _read(tmp_path / f"{fold}_pass", records={fold: record}, x1=_scores(), r12=_scores(), x2=_scores())
    assert not res2["arms"]["raw_did"]["supported"] and res2["arms"]["raw_did"]["verdict"].startswith("NOT READABLE")


def test_raw_did_not_built_on_x2_only_is_labelled_and_read_as_before(tmp_path):
    rec = {"x2": {"source_mesp1_labels_in_contrast": 0, "source_nonmesp1_labels": 3}}
    res = _read(tmp_path / "x2", records=rec, x1=_scores(), r12=_scores(), x2=_scores(did=52.0, perm=52.0))
    assert res["arms"]["raw_did"]["readable"] and res["arms"]["raw_did"]["supported"]
    assert res["tasks"]["x2"]["not_built"] == ["raw_did", "raw_did_perm"]
    assert "not built on this fold: raw_did, raw_did_perm" in k11.text(res)


def test_ahead_in_eight_seeds_is_required(tmp_path):
    sc = _scores()
    sc["raw_c3"] = [70.0] * 7 + [40.0] * (len(SEEDS) - 7)
    res = _read(tmp_path / "a", x1=sc, r12=_scores(), x2=_scores())
    assert not res["arms"]["raw_c3"]["supported"] and res["tasks"]["x1"]["paired"]["raw_c3 - raw"]["ahead"] == 7


def test_the_reading_checks_k9s_committed_raw_scores(tmp_path):
    k9_dir = tmp_path / "k9"
    k9_dir.mkdir()
    for name in k11.TASKS:
        (k9_dir / f"{name}.json").write_text(json.dumps(_result(k9.TASKS[name].key, {"raw": [52.0] * len(SEEDS)})))
    ok = _read(tmp_path / "ok", k9_dir, x1=_scores(), r12=_scores(), x2=_scores())
    assert not ok["void"] and ok["k9_reproduction"] == "reproduced"
    near = _read(tmp_path / "near", k9_dir, x1=_scores(raw=52.0 + 5e-7), r12=_scores(), x2=_scores())
    assert not near["void"]
    off = _read(tmp_path / "off", k9_dir, x1=_scores(raw=52.01), r12=_scores(), x2=_scores())
    assert off["void"] and any("does not reproduce" in w for w in off["why"])
    (k9_dir / "x2.json").unlink()
    gone = _read(tmp_path / "gone", k9_dir, x1=_scores(), r12=_scores(), x2=_scores())
    assert gone["void"] and any("missing" in w for w in gone["why"])


def test_the_reading_is_void_on_missing_dirty_or_mixed_results(tmp_path):
    missing = _read(tmp_path / "m", x1=_scores(), r12=_scores())
    assert missing["void"] and any("missing task x2" in w for w in missing["why"])
    d = tmp_path / "dirty"
    d.mkdir()
    for name in k11.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k11.TASKS[name].key, _scores(), clean=name != "x2")))
    assert k11.reading(d, tmp_path / "none")["void"]
    d = tmp_path / "commits"
    d.mkdir()
    for i, name in enumerate(k11.TASKS):
        (d / f"{name}.json").write_text(json.dumps(_result(k11.TASKS[name].key, _scores(), commit=f"c{i % 2}")))
    assert k11.reading(d, tmp_path / "none")["void"]


def test_run_and_read_refuse_to_overwrite(tmp_path):
    (tmp_path / "x1.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k11.main(["run", "--task", "x1", "--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])
    out = tmp_path / "out"
    out.mkdir()
    (out / "reading.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k11.main(["read", "--runs", str(tmp_path), "--out", str(out)])
