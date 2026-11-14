"""Tests for tries K14 (experiments/t3_decomp/k14_form_arms.py), on synthetic data only: no competition file is opened,
nothing is scored, and no prediction is compared with any target."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

T3 = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import eb_factor as eb  # noqa: E402
import k7_wt_revert as k7  # noqa: E402
import k9_residual as k9  # noqa: E402
import k14_form_arms as k14  # noqa: E402
import test_k9_residual as t9  # noqa: E402  (K9's synthetic world and reference, reused)

CLIP = t9.CLIP
world = t9.world  # the module fixture: synthetic files in task x1's roles


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


def _prep(rng, phi_raw=None, y=None):
    y = rng.normal(0, 0.5, t9.G) if y is None else y
    den = rng.gamma(2.0, size=t9.G)
    phi_raw = eb.finalize(y, den, CLIP) if phi_raw is None else phi_raw
    return {"phi_raw": phi_raw, "phi_half_ctrl": k14.half_ctrl_factor(phi_raw), "clip": CLIP}


def _counts(bs, ref, n, seed):
    """The drawn cells' counts and smoothed proportions on K9's stream (construct without a factor)."""
    _, drawn, C, P = k14._drawn(bs, ref, n, seed, ref.V.p)
    return drawn, C, P


# ----------------------------------------------------------------------------------------------- constants and tasks
def test_the_fixed_constants_and_arms():
    assert k14.PERM_SEED == k9.PERM_SEED == 20261113 and k14.STREAM == k9.STREAM == harness.K4_STREAM
    assert (k14.MIN_GAIN, k14.MIN_AHEAD, k14.X2_FLOOR) == (0.3, 8, -0.3)
    assert k14.ARMS == ("raw", "raw_half", "raw_half_ctrl", "raw_xt", "raw_mn")
    assert k14.CONTROLS == {"raw_half": "raw_half_ctrl", "raw_xt": None, "raw_mn": None}
    assert len({k14.HALF_STREAM, k14.XT_STREAM, k14.MN_STREAM, k14.STREAM, f"{k14.STREAM}:assemble"}) == 5
    assert k14.TASKS is k9.TASKS and k14.RESULTS_K9.is_dir()
    assert [c.name for c in k14.candidates()] == [harness.FLOOR, *k14.ARMS]
    assert all(c.applies(t) for c in k14.candidates()[1:] for t in k14.TASKS.values())
    assert not any(c.applies(t) for c in k14.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)


# ----------------------------------------------------------------------------------------------- the half
@pytest.mark.parametrize("n", [200, 201, 1])
def test_the_half_is_floor_n_over_2_cells_and_another_seed_gives_another_half(n):
    m = k14.half_mask(n, 11)
    assert m.dtype == bool and m.shape == (n,) and m.sum() == n // 2
    assert np.array_equal(k14.half_mask(n, 11), m)  # reproducible from the build seed
    if n > 2:
        assert not np.array_equal(k14.half_mask(n, 12), m)


def test_raw_half_moves_the_half_and_keeps_the_others_as_drawn(bs):
    rng = np.random.default_rng(1)
    ref = t9._reference(bs, rng)
    prep = _prep(rng)
    drawn, C, _ = _counts(bs, ref, 200, 7)
    built = k14.arm_built(bs, ref, 200, 7, ref.V.p, "raw_half", prep)
    m = k14.half_mask(200, 7)
    assert np.array_equal(built.row, drawn.row) and np.array_equal(built.lab, drawn.lab)
    assert np.array_equal(built.X[~m], bs.t3d.lognorm(C[~m]))  # the other half: the reference's counts as drawn
    assert not (built.X[m] == bs.t3d.lognorm(C[m])).all(1).any()  # every moved cell changed
    # the half: apply_gene_factor on those cells, continuing K9's generator after the draws
    g = bs.t3d.stable_rng(7, k14.STREAM)
    bs.construct(ref, 200, g, f=ref.V.p)
    assert np.array_equal(built.X[m], bs.t3d.lognorm(bs.t3d.apply_gene_factor(C[m], prep["phi_raw"],
                                                                              ref.pool.P[drawn.row][m], g)))


def test_the_mask_of_all_cells_gives_raw_bit_for_bit_at_the_built_cells(bs):
    rng = np.random.default_rng(2)
    ref = t9._reference(bs, rng)
    prep = _prep(rng)
    raw = k14.arm_built(bs, ref, 200, 9, ref.V.p, "raw", prep)
    k9raw = k9.arm_built(bs, ref, 200, 9, ref.V.p, "raw", prep["phi_raw"], None, None, CLIP)
    assert np.array_equal(raw.X, k9raw.X) and np.array_equal(raw.row, k9raw.row)
    every = k14.arm_built(bs, ref, 200, 9, ref.V.p, "raw_half", prep, mask=np.ones(200, dtype=bool))
    assert np.array_equal(every.X, raw.X) and np.array_equal(every.row, raw.row) and np.array_equal(every.lab, raw.lab)
    half = k14.arm_built(bs, ref, 200, 9, ref.V.p, "raw_half", prep)
    assert not np.array_equal(half.X, raw.X)


# ----------------------------------------------------------------------------------------------- the control
def test_the_controls_factor_is_one_plus_phi_raw_over_two_within_the_clip_range():
    rng = np.random.default_rng(3)
    y = rng.normal(0, 1.5, 300)  # wide enough that genes reach both clip bounds
    den = rng.gamma(2.0, size=300)
    raw = eb.finalize(y, den, CLIP)
    raw[:5] = 1.0
    assert (raw == CLIP[0]).any() and (raw == CLIP[1]).any()
    ctrl = k14.half_ctrl_factor(raw)
    assert ctrl.dtype == np.float64 and np.array_equal(ctrl, (1.0 + raw) / 2.0)
    assert ((CLIP[0] <= ctrl) & (ctrl <= CLIP[1])).all()
    assert ctrl.min() == (1.0 + CLIP[0]) / 2.0 and ctrl.max() == (1.0 + CLIP[1]) / 2.0
    # the control differs from raw wherever phi_raw is not 1 and equals it where phi_raw is exactly 1
    assert np.all(ctrl[raw != 1.0] != raw[raw != 1.0]) and np.all(ctrl[raw == 1.0] == 1.0)


def _mean_case():
    """Genes at both clip bounds, moved down, moved up, and at exactly 1. den puts all its weight on the gene at 1, so
    finalize's rescale is exactly 1 and phi_raw = clip(exp(y)): the old control finalize(y / 2) = clip(sqrt(exp(y)))."""
    phi_in = np.array([0.1, 0.3, 0.7, 1.0, 1.0, 1.5, 3.0, 10.0])  # 0.1 and 10 land on the clip bounds
    y = np.log(phi_in)
    den = np.zeros(len(y))
    den[3] = 1.0
    phi_raw = eb.finalize(y, den, CLIP)
    old = eb.finalize(y / 2, den, CLIP)
    rng = np.random.default_rng(31)
    n = 200
    C = rng.poisson(6.0, size=(n, len(y))).astype(np.int32)
    C[:, 3] = np.maximum(C[:, 3], 1)  # every cell keeps a count at a gene with factor 1: no cell is ever emptied
    P = rng.dirichlet(np.ones(len(y)), size=n).astype(np.float32)
    return phi_raw, old, C, P


def _per_gene_means(C, phi, P, N, mask_of=None):
    """Per seed, the mean over the cells of each gene's built count: every cell through apply_gene_factor with phi,
    or (mask_of given) only the cells of mask_of(seed), the others as drawn."""
    t3d = harness._t3_head().t3d
    out = np.empty((N, C.shape[1]))
    for s in range(N):
        rng = np.random.default_rng([s, 3])
        if mask_of is None:
            b = t3d.apply_gene_factor(C, phi, P, rng)
        else:
            m = mask_of(s)
            b = C.copy()
            b[m] = t3d.apply_gene_factor(C[m], phi, P[m], rng)
        out[s] = b.mean(0)
    return out


def test_the_controls_mean_matches_raw_halfs_mean_per_gene_and_the_old_controls_does_not():
    phi_raw, old, C, P = _mean_case()
    assert phi_raw[0] == CLIP[0] and phi_raw[-1] == CLIP[1] and np.all(phi_raw[3:5] == 1.0)
    N = 1500
    half = _per_gene_means(C, phi_raw, P, N, mask_of=lambda s: k14.half_mask(len(C), s))
    ctrl = _per_gene_means(C, k14.half_ctrl_factor(phi_raw), P, N)
    prev = _per_gene_means(C, old, P, N)
    se = np.sqrt(half.var(0) / N + ctrl.var(0) / N)
    # tolerance: 5 Monte Carlo standard errors of the difference, plus 1e-9 for genes that never vary (factor 1)
    assert np.all(np.abs(half.mean(0) - ctrl.mean(0)) <= 5 * se + 1e-9)
    assert np.array_equal(half[:, 3:5], ctrl[:, 3:5])  # genes at exactly 1: as drawn in both
    # why the control changed: the old control finalize(y / 2) misses the mixture's mean on the genes moved more than
    # twofold, by more than 5 Monte Carlo standard errors of the difference
    big = (phi_raw > 2.0) | (phi_raw < 0.5)
    assert big.sum() == 4
    se_old = np.sqrt(half.var(0) / N + prev.var(0) / N)
    gap = prev.mean(0) - half.mean(0)
    assert np.all(np.abs(gap[big]) > 5 * se_old[big])


def test_the_controls_built_counts_differ_from_raws_where_phi_raw_is_not_1_and_equal_them_where_it_is(bs):
    phi_raw, _, C, P = _mean_case()
    t3d = bs.t3d
    g = np.random.default_rng(8)
    raw = t3d.apply_gene_factor(C, phi_raw, P, np.random.default_rng(8))
    ctrl = t3d.apply_gene_factor(C, k14.half_ctrl_factor(phi_raw), P, g)
    one = phi_raw == 1.0
    assert np.array_equal(ctrl[:, one], raw[:, one]) and np.array_equal(ctrl[:, one], C[:, one])
    assert np.all((ctrl[:, ~one] != raw[:, ~one]).any(0))


def test_the_control_is_built_through_k9s_raw_path_and_its_check_counts_the_changed_cells(bs):
    rng = np.random.default_rng(4)
    ref = t9._reference(bs, rng)
    prep = _prep(rng)
    chk: dict = {}
    built = k14.arm_built(bs, ref, 200, 5, ref.V.p, "raw_half_ctrl", prep, check=chk)
    want = k9.arm_built(bs, ref, 200, 5, ref.V.p, "raw", prep["phi_half_ctrl"], None, None, CLIP)
    assert np.array_equal(built.X, want.X)
    assert set(chk) == set(k14.CHECK_KEYS) and 0.0 < chk["share"] <= 1.0 and chk["mean_abs_off"] == 0.0


# ----------------------------------------------------------------------------------------------- the expected counts
def _small_case(seed=5):
    rng = np.random.default_rng(seed)
    G = 6
    C = rng.poisson(1.5, size=(8, G)).astype(np.int32)
    C[0] = [1, 0, 0, 0, 0, 0]  # a cell that can be left empty and restored
    C[1] = [0, 0, 0, 0, 0, 0]  # an empty cell
    C[2] = [2, 1, 0, 0, 0, 3]
    phi = np.array([0.25, 0.5, 1.0, 2.0, 3.5, 1.0])
    P = rng.dirichlet(np.ones(G), size=8).astype(np.float32)
    P[:, 4] = np.where(np.arange(8) % 2 == 0, 0.0, P[:, 4])  # some top-up entries with no smoothed mass
    return C, phi, P


def test_the_expected_counts_match_a_monte_carlo_mean_of_apply_gene_factor(bs):
    C, phi, P = _small_case()
    e = k14.expected_counts(C, phi, P)
    N = 20000
    rng = np.random.default_rng(0)
    draws = np.stack([bs.t3d.apply_gene_factor(C, phi, P, rng) for _ in range(N)]).astype(np.float64)
    mean, se = draws.mean(0), draws.std(0) / np.sqrt(N)
    # tolerance: 5 Monte Carlo standard errors plus 1e-9 for entries that never vary
    assert np.all(np.abs(mean - e) <= 5 * se + 1e-9)
    assert e[0, 0] > C[0, 0] * phi[0]  # the restore of the emptied cell is in the expectation
    assert np.array_equal(e[1], np.zeros(C.shape[1]))
    assert np.all(e[(C == 0) & (phi <= 1.0)[None, :]] == 0)  # a zero thinned or left: expected count 0


# ----------------------------------------------------------------------------------------------- the repair
def test_the_repair_brings_every_cells_moved_total_to_its_rounded_expectation(bs):
    rng = np.random.default_rng(6)
    ref = t9._reference(bs, rng, n=300)
    prep = _prep(rng)
    _, C, P = _counts(bs, ref, 300, 13)
    g = bs.t3d.stable_rng(13, k14.STREAM)
    bs.construct(ref, 300, g, f=ref.V.p)
    out, raw, chk = k14.kernel_counts(bs.t3d, C, prep["phi_raw"], P, g, "raw_xt", 13)
    e = k14.expected_counts(C, prep["phi_raw"], P)
    moved = k14.moved_entries(prep["phi_raw"], e)
    target = np.rint(np.where(moved, e, 0.0).sum(1))
    assert np.array_equal(np.where(moved, out, 0).sum(1), target)
    assert np.array_equal(out[~moved], raw[~moved])  # entries the step does not move: as drawn
    assert np.all(out[e == 0] == 0)
    have = np.where(moved, raw, 0).sum(1)
    at = have == target
    assert at.any() and (~at).any()
    assert np.array_equal(out[at], raw[at])  # a cell already at its total: untouched
    assert chk["share"] == pytest.approx(float((~at).mean()))
    assert chk["mean_abs_off"] == pytest.approx(float(np.abs(have - target).mean()))
    # raw's build of the same cells is K9's raw build
    k9raw = k9.arm_built(bs, ref, 300, 13, ref.V.p, "raw", prep["phi_raw"], None, None, CLIP)
    assert np.array_equal(bs.t3d.lognorm(raw), k9raw.X)


@pytest.mark.parametrize("arm", ["raw_xt", "raw_mn"])
def test_the_kernel_switched_off_is_raw_bit_for_bit(bs, arm):
    rng = np.random.default_rng(7)
    ref = t9._reference(bs, rng)
    prep = _prep(rng)
    raw = k14.arm_built(bs, ref, 200, 3, ref.V.p, "raw", prep)
    off = k14.arm_built(bs, ref, 200, 3, ref.V.p, arm, prep, repair=False)
    assert np.array_equal(off.X, raw.X) and np.array_equal(off.row, raw.row) and np.array_equal(off.lab, raw.lab)
    on = k14.arm_built(bs, ref, 200, 3, ref.V.p, arm, prep)
    assert not np.array_equal(on.X, raw.X)


# ----------------------------------------------------------------------------------------------- the redraw
def test_the_redraw_keeps_the_rounded_totals_unmoved_entries_and_zeros(bs):
    C, phi, P = _small_case()
    e = k14.expected_counts(C, phi, P)
    moved = k14.moved_entries(phi, e)
    assert not moved[:, phi == 1.0].any() and not moved[e == 0].any()
    S = np.where(moved, e, 0.0).sum(1)
    N = 20000
    outs = []
    for s in range(N):
        out, raw, _ = k14.kernel_counts(bs.t3d, C, phi, P, np.random.default_rng([s, 1]), "raw_mn", s)
        live = out.sum(1) > 0
        assert np.array_equal(np.where(moved, out, 0).sum(1)[live], np.rint(S)[live])
        assert np.array_equal(out[:, phi == 1.0], raw[:, phi == 1.0])  # a factor of exactly 1: as drawn
        assert np.all(out[e == 0] == 0)  # expected count 0 stays 0
        outs.append(out)
    restored = (np.rint(S) == 0) & (np.where(~moved, C, 0).sum(1) == 0)
    mean = np.stack(outs).astype(np.float64).mean(0)
    sd = np.stack(outs).astype(np.float64).std(0)
    ok = ~restored
    # each moved entry's mean is e rint(S) / S (H19's form); tolerance 5 Monte Carlo standard errors plus 1e-9
    want = np.where(moved, e * np.divide(np.rint(S), S, out=np.ones_like(S), where=S > 0)[:, None], C)
    assert np.all(np.abs(mean - want)[ok] <= 5 * sd[ok] / np.sqrt(N) + 1e-9)
    # and the expected counts themselves, up to that rounding of the total (|e rint(S) / S - e| per entry)
    mv = ok[:, None] & moved
    assert np.all(np.abs(mean - e)[mv] <= np.abs(want - e)[mv] + 5 * sd[mv] / np.sqrt(N) + 1e-9)


# ----------------------------------------------------------------------------------------------- the streams
def test_the_new_streams_differ_from_every_stream_the_build_uses_and_from_each_other(bs):
    seed = harness.BUILD_SEED_OFFSET + 3
    t3d = bs.t3d
    new = {name: t3d.stable_rng(seed, name).random(64) for name in (k14.HALF_STREAM, k14.XT_STREAM, k14.MN_STREAM)}
    old = [t3d.stable_rng(seed, k14.STREAM).random(64), t3d.stable_rng(seed, f"{k14.STREAM}:assemble").random(64),
           np.random.default_rng(k14.PERM_SEED).random(64), np.random.default_rng([k14.PERM_SEED, 2]).random(64),
           *(np.random.default_rng([k14.PERM_SEED, 1, p]).random(64) for p in range(5))]
    arrays = list(new.values()) + old
    for i in range(len(arrays)):
        for j in range(i + 1, len(arrays)):
            assert not np.array_equal(arrays[i], arrays[j]), (i, j)
    # the half's draw differs between build seeds and from the same draw on the build's own stream
    assert not np.array_equal(k14.half_mask(200, seed), k14.half_mask(200, seed + 1))
    m = np.zeros(200, dtype=bool)
    m[t3d.stable_rng(seed, k14.STREAM).choice(200, 100, replace=False)] = True
    assert not np.array_equal(m, k14.half_mask(200, seed))


# ----------------------------------------------------------------------------------------------- the check line
def test_the_check_line_fields_and_the_task_record():
    assert k14.check_line("x1", "raw_xt", {"share": 0.25, "mean_abs_off": 1.5}) == (
        "k14 check: task=x1 arm=raw_xt share=0.2500 mean_abs_off=1.5000")
    per = {"raw": {1: {"share": 0.0, "mean_abs_off": 0.0}},
           "raw_xt": {2: {"share": 0.5, "mean_abs_off": 2.0}, 1: {"share": 0.25, "mean_abs_off": 1.0}}}
    rec = k14.task_check(per)
    assert rec["raw_xt"]["share"] == pytest.approx(0.375) and rec["raw_xt"]["mean_abs_off"] == pytest.approx(1.5)
    assert list(rec["raw_xt"]["per_seed"]) == ["1", "2"] and "raw_half" not in rec


# ----------------------------------------------------------------------------------------------- end to end, synthetic
def test_end_to_end_raw_is_k9s_raw_and_every_arm_is_built(world, monkeypatch):
    monkeypatch.setattr(k7, "MIN_CELLS", 20)  # K9's own build of the same fold needs its groups (as k13's test)
    t = k14.TASKS["x1"]
    k14._PREP.clear()
    k14._CHECK.clear()
    k9._PREP.clear()
    ctx = t9._ctx(t, world)
    with harness.guarded_reads() as g:
        out = {arm: k14._build(ctx, arm) for arm in k14.ARMS}
        k9raw = k9._build(ctx, k9.BASE)
    assert not g["refused"]
    prep = k14._PREP["fold"]
    assert np.array_equal(prep["phi_raw"], k9._PREP["fold"]["phi_raw"])
    assert np.array_equal(prep["phi_half_ctrl"], (1.0 + prep["phi_raw"]) / 2.0)
    assert np.array_equal(out["raw"].X, k9raw.X) and np.array_equal(out["raw"].coords, k9raw.coords)
    for a in k14.ARMS[1:]:
        assert not np.array_equal(out[a].X, out["raw"].X), a
    rec = prep["record"]
    assert json.loads(json.dumps(rec)) == rec and rec["max_abs_log_finalize_y_over_phi_raw"] >= 0.0
    chk = k14.task_check(k14._CHECK)
    assert set(chk) == set(k14.ARMS) and chk["raw"]["share"] == 0.0
    assert chk["raw_xt"]["mean_abs_off"] > 0 and all(chk[a]["mean_abs_off"] == 0 for a in k14.ARMS if a != "raw_xt")
    assert json.loads(json.dumps(chk)) == chk
    k14._PREP.clear()
    k14._CHECK.clear()
    k9._PREP.clear()


# ----------------------------------------------------------------------------------------------- the reading
def _check(clean=True):
    return {a: {"share": 0.0 if a == "raw" else 0.5, "mean_abs_off": 1.0 if a == "raw_xt" else 0.0,
                "per_seed": {}} for a in k14.ARMS}


def _result(key, scores, clean=True, commit="abc", check=True):
    rows = [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5, "severity": 0.6} for a in scores},
             "t3": {a: {"severity": {"state": "scored"}} for a in scores}} for i, s in enumerate(SEEDS)]
    out = {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
           "factor_record": {}}
    if check:
        out["k14_check"] = _check()
    return out


def _scores(half=54.0, ctrl=50.0, xt=54.0, mn=54.0, raw=52.0):
    n = len(SEEDS)
    return {harness.FLOOR: [50.0] * n, "raw": [raw] * n, "raw_half": [half] * n, "raw_half_ctrl": [ctrl] * n,
            "raw_xt": [xt] * n, "raw_mn": [mn] * n}


def _k9_dir(d: Path, raw: float = 52.0) -> Path:
    d.mkdir()
    for name in k14.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k9.TASKS[name].key, {"raw": [raw] * len(SEEDS)})))
    return d


def _read(d: Path, k9_dir=None, **per_task):
    d.mkdir()
    for name, sc in per_task.items():
        (d / f"{name}.json").write_text(json.dumps(_result(k14.TASKS[name].key, sc)))
    return k14.reading(d, k9_dir if k9_dir is not None else _k9_dir(d / "k9"))


def test_the_reading_supported(tmp_path):
    res = _read(tmp_path / "yes", x1=_scores(), r12=_scores(), x2=_scores(half=51.8, xt=51.8, mn=51.8))
    assert not res["void"] and res["k9_reproduction"] == "reproduced"
    assert all(res["arms"][a]["supported"] and res["arms"][a]["verdict"].startswith("SUPPORTED") for a in k14.READ)
    txt = k14.text(res)
    assert "NOT the adoption rule" in txt and "k14 check: task=x1 arm=raw_xt share=0.5000 mean_abs_off=1.0000" in txt
    assert res["control"] == {"raw_half_ctrl": k14.HALF_CTRL_RULE} and "(1 + phi_raw) / 2" in txt
    assert "exactly as far as raw_half's mixture does" in txt


def test_the_reading_closed_by_x2(tmp_path):
    res = _read(tmp_path / "x2", x1=_scores(), r12=_scores(), x2=_scores(half=51.6, xt=51.6, mn=51.6))
    assert not res["void"]
    for a in k14.READ:
        assert not res["arms"][a]["supported"] and res["arms"][a]["verdict"] == "CLOSED: no request"
        assert not res["arms"][a]["checks"][f"x2: {a} - raw not below -0.3"]


@pytest.mark.parametrize("fold", ["x1", "r12"])
def test_the_reading_closed_by_the_control_for_raw_half_only(tmp_path, fold):
    sc = {n: _scores() for n in k14.TASKS}
    sc[fold] = _scores(ctrl=53.9)
    res = _read(tmp_path / fold, **sc)
    half = res["arms"]["raw_half"]
    assert not res["void"] and not half["supported"] and not half["checks"][f"{fold}: raw_half - raw_half_ctrl"]
    assert half["checks"][f"{fold}: raw_half - raw"]
    # the kernels have no control clause: a strong raw_half_ctrl changes nothing for them
    for a in ("raw_xt", "raw_mn"):
        assert res["arms"][a]["supported"] and not any("ctrl" in k for k in res["arms"][a]["checks"])
        assert len(res["arms"][a]["checks"]) == 3


def test_the_reading_needs_eight_seeds_ahead(tmp_path):
    sc = _scores()
    sc["raw_xt"] = [70.0] * 7 + [40.0] * (len(SEEDS) - 7)
    res = _read(tmp_path / "a", x1=sc, r12=_scores(), x2=_scores())
    assert not res["arms"]["raw_xt"]["supported"] and res["arms"]["raw_mn"]["supported"]


def test_the_reading_is_void_on_a_raw_mismatch(tmp_path):
    k9_dir = _k9_dir(tmp_path / "k9")
    ok = _read(tmp_path / "ok", k9_dir, x1=_scores(), r12=_scores(), x2=_scores())
    assert not ok["void"]
    off = _read(tmp_path / "off", k9_dir, x1=_scores(raw=52.01), r12=_scores(), x2=_scores())
    assert off["void"] and any("does not reproduce" in w for w in off["why"]) and k14.text(off).startswith("VOID")


def test_the_reading_is_void_without_k9s_committed_scores(tmp_path):
    res = _read(tmp_path / "r", tmp_path / "absent", x1=_scores(), r12=_scores(), x2=_scores())
    assert res["void"] and res["why"] == [f"K9's committed raw scores ({tmp_path / 'absent'}) are absent from this "
                                          f"checkout"]


def test_the_reading_is_void_on_missing_arm_seed_check_or_dirty_tree(tmp_path):
    k9_dir = _k9_dir(tmp_path / "k9")
    assert _read(tmp_path / "m", x1=_scores(), r12=_scores())["void"]
    sc = _scores()
    del sc["raw_mn"]
    assert _read(tmp_path / "arm", x1=sc, r12=_scores(), x2=_scores())["void"]
    for case in ("seed", "dirty", "check"):
        d = tmp_path / case
        d.mkdir()
        for name in k14.TASKS:
            r = _result(k14.TASKS[name].key, _scores(), clean=not (case == "dirty" and name == "x2"),
                        check=not (case == "check" and name == "r12"))
            if case == "seed" and name == "x1":
                r["rows"] = r["rows"][:-1]
            (d / f"{name}.json").write_text(json.dumps(r))
        assert k14.reading(d, k9_dir)["void"], case


# ----------------------------------------------------------------------------------------------- the commands
def test_run_writes_the_check_numbers_and_prints_the_lines(tmp_path, monkeypatch, capsys):
    def fake(*a, **kw):
        for arm in k14.ARMS:
            k14._CHECK.setdefault(arm, {})[1] = {"share": 0.5, "mean_abs_off": 0.0}
        return {"rows": []}

    monkeypatch.setattr(harness, "run_task", fake)
    monkeypatch.setitem(k14._PREP, "fold", {"record": {"fold": "f"}})
    runs = tmp_path / "runs"
    k14.main(["run", "--task", "r12", "--out", str(runs), "--data", str(tmp_path / "nowhere")])
    res = json.loads((runs / "r12.json").read_text())
    assert set(res["k14_check"]) == set(k14.ARMS) and res["factor_record"] == {"fold": "f"}
    printed = capsys.readouterr().out
    for arm in k14.ARMS:
        assert f"k14 check: task=r12 arm={arm} share=0.5000 mean_abs_off=0.0000" in printed
    k14._CHECK.clear()


def test_run_and_read_refuse_to_overwrite(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(harness, "run_task", lambda *a, **kw: called.append(a))
    runs = tmp_path / "runs"
    runs.mkdir()
    (runs / "x1.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k14.main(["run", "--task", "x1", "--out", str(runs), "--data", str(tmp_path / "nowhere")])
    assert not called and (runs / "x1.json").read_text() == "{}"
    out = tmp_path / "out"
    out.mkdir()
    for f in ("reading.json", "reading.txt"):
        (out / f).write_text("{}")
        with pytest.raises(SystemExit, match="exists"):
            k14.main(["read", "--runs", str(runs), "--out", str(out)])
        assert (out / f).read_text() == "{}"
        (out / f).unlink()
