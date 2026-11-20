"""Try K21 (experiments/t3_decomp/k21_factor_forms.py) and its shared transforms (factor_forms.py). Synthetic data,
hand-computed cases and hand-made score tables only: no competition file is loaded and no fold is scored."""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3_DIR = REPO / "experiments" / "t3_decomp"
if str(T3_DIR) not in sys.path:
    sys.path.insert(0, str(T3_DIR))
import factor_forms as ff  # noqa: E402
import k7_wt_revert as k7  # noqa: E402
import k8_s2_form as k8  # noqa: E402
import k16_revert_dose as k16  # noqa: E402
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402

G, N_POOL, N_DRAW = 14, 90, 150
GENES = [f"g{i}" for i in range(G)]
LN2 = math.log(2.0)


# ------------------------------------------------------------------------------------------------ constants


def test_arms_constants_and_the_shared_module():
    assert k21.ARMS == ("barred_form", "top100", "uni_top", "cap90", "uni_cap", "wtprod", "wtprod_perm")
    assert k21.BASE == k18.BARRED and k21.STREAM == k18.STREAM == k8.STREAM == "t3_s2_form"
    assert k21.NAMEABLE == ("top100", "cap90", "wtprod")
    assert k21.CONTROL_OF == {"top100": "uni_top", "cap90": "uni_cap", "wtprod": "wtprod_perm"}
    assert (k21.MIN_GAIN, k21.MIN_AHEAD, k21.MAX_EXTRA_GATED, k21.GUARD_MEAN, k21.GUARD_SEED) == (0.3, 8, 1, 0.1, 0.5)
    assert (ff.TOP_K, ff.CAP_Q, ff.MIN_CHANGED, ff.CLIP_TOL) == (100, 0.9, 10, 1e-9)
    assert ff.WT_EXPONENT == k7.EXPONENT == 0.5 and ff.CLIP == k7.CLIP == (0.25, 4.0)
    assert ff.LN_CLIP == pytest.approx(math.log(4.0), abs=0) and -math.log(ff.CLIP[0]) == pytest.approx(ff.LN_CLIP)
    assert k21.ff is ff  # the practice arms take their transforms from the module the builder applies
    assert k21.TASKS == k18.TASKS and k21.WORLDS == k18.WORLDS


def test_the_guards_reference_is_k18s_committed_files_at_their_pins():
    for f, sha in k21.K18_SHA256.items():
        assert hashlib.sha256((k21.K18_TASKS / f).read_bytes()).hexdigest() == sha
    rows, why = k21.k18_rows()
    assert not why and set(rows) == {(w, n) for w in k21.WORLDS for n in k21.TASKS}
    for v in rows.values():  # every seed holds K18's barred_form (values not read here)
        assert sorted(r["seed"] for r in v) == list(SEEDS) and all(k21.BASE in r["scores"] for r in v)


# ---------------------------------------------------------------------------------- the transforms, by hand


def test_top_k_by_hand():
    """|ln phi| = ln2, ln2, ln1.5, ln3, ln4, ln1.1; in descending order ln4 (gene 4), ln3 (gene 3), ln2 (genes 0 and
    1), ln1.5, ln1.1. K = 2 cuts between ln3 and ln2; K = 3 cuts between the two ln2 (a tie: refused); K = 4 keeps
    both ln2 genes."""
    phi = np.array([2.0, 0.5, 1.5, 3.0, 0.25, 1.1])
    out, rec = ff.top_k(phi, 2)
    assert np.array_equal(out, [1.0, 1.0, 1.0, 3.0, 0.25, 1.0])
    assert rec["genes_at_clip"] == 1 and rec["k"] == 2 and rec["cut"] == pytest.approx(math.log(3.0))
    with pytest.raises(ff.FormRefused, match="ties across the cut"):
        ff.top_k(phi, 3)
    out, _ = ff.top_k(phi, 4)
    assert np.array_equal(out, [2.0, 0.5, 1.0, 3.0, 0.25, 1.0])
    assert np.array_equal(ff.top_k(phi, 6)[0], phi) and np.array_equal(ff.top_k(phi, 50)[0], phi)  # K >= n: phi_M
    # the kept set does not depend on the gene order when no tie straddles the cut
    perm = np.array([5, 2, 0, 4, 1, 3])
    assert np.array_equal(ff.top_k(phi[perm], 2)[0], ff.top_k(phi, 2)[0][perm])


def test_top_k_refuses_when_the_clip_fills_the_kept_set():
    phi = np.array([4.0, 0.25, 4.0, 2.0, 1.0])
    with pytest.raises(ff.FormRefused, match="at the clip"):
        ff.top_k(phi, 3)
    with pytest.raises(ff.FormRefused, match="at the clip"):
        ff.top_k(phi, 2)
    out, rec = ff.top_k(phi, 4)
    assert np.array_equal(out, [4.0, 0.25, 4.0, 2.0, 1.0]) and rec["genes_at_clip"] == 3


def test_cap_by_hand():
    phi = np.array([4.0, 0.25, 2.0, 0.5, 1.5, 1.0])
    out, rec = ff.cap(phi, 0.5)  # off the clip |ln|: ln2, ln2, ln1.5, 0 -> median (ln2 + ln1.5) / 2
    t = (LN2 + math.log(1.5)) / 2.0
    assert rec["t"] == pytest.approx(t, rel=1e-15) and rec["genes_off_clip"] == 4 and rec["genes_at_clip"] == 2
    want = np.clip(phi, math.exp(-t), math.exp(t))
    assert np.array_equal(out, want)
    assert out[5] == 1.0 and out[4] == 1.5 and out[0] == pytest.approx(math.exp(t))
    q9 = float(np.quantile([LN2, LN2, math.log(1.5), 0.0], 0.9))  # numpy linear: 3 * 0.9 = 2.7 -> ln2
    assert ff.cap(phi)[1]["t"] == pytest.approx(q9, rel=1e-15) and q9 == pytest.approx(LN2)
    big = float(np.abs(np.log(phi)).max())
    assert np.array_equal(ff.cap_at(phi, big), phi) and np.array_equal(ff.cap_at(phi, big + 1.0), phi)
    with pytest.raises(ff.FormRefused, match="off the clip"):
        ff.cap(np.array([4.0, 0.25, 4.0]))


def test_uniform_and_the_controls_hold_their_tries_size():
    phi = np.array([4.0, 0.5, 2.0, 1.0])  # ln: a[2, -1, 1, 0]
    target = np.array([2.0, 0.5, 1.0, 1.0])  # ln: a[1, -1, 0, 0]
    out, rec = ff.uniform(phi, target)
    s = np.std([1.0, -1.0, 0.0, 0.0]) / np.std([2.0, -1.0, 1.0, 0.0])
    assert rec["s"] == pytest.approx(s, rel=1e-14) and np.allclose(out, phi ** s, rtol=1e-15)
    rng = np.random.default_rng(2)
    phi = np.clip(np.exp(rng.normal(0, 0.8, 300)), 0.25, 4.0)
    top, _ = ff.vector_form("top100", phi)
    uni_top, rt = ff.vector_form("uni_top", phi)
    cap, _ = ff.vector_form("cap90", phi)
    uni_cap, rc = ff.vector_form("uni_cap", phi)
    assert ff.sd(np.log(uni_top)) == pytest.approx(ff.sd(np.log(top)), rel=1e-12)
    assert ff.sd(np.log(uni_cap)) == pytest.approx(ff.sd(np.log(cap)), rel=1e-12)
    assert 0 < rt["s"] < 1 and 0 < rc["s"] < 1  # exponents below 1: controls only
    assert int((top == 1.0).sum()) == 200 and int((top != 1.0).sum()) == 100
    assert np.array_equal(ff.vector_form("identity", phi)[0], phi)
    with pytest.raises(ff.FormRefused, match="SD"):
        ff.uniform(np.ones(5), np.ones(5))
    with pytest.raises(ff.FormRefused, match="not a vector form"):
        ff.vector_form("soft50s", phi)
    with pytest.raises(ff.FormRefused, match="finite positive"):
        ff.vector_form("cap90", np.array([1.0, 0.0, 2.0]))


def test_at_exponent_and_wtprod_factors_by_hand():
    phi = np.array([2.0, 0.5, 4.0])
    assert ff.at_exponent(phi, 1.0) is phi
    assert np.array_equal(ff.at_exponent(phi, 2.0), [4.0, 0.25, 16.0])
    per = {"A": np.array([4.0, 1.0, 0.25]), "B": np.array([1.0, 4.0, 1.0])}
    pooled = np.array([1.0, 1.0, 4.0])
    p2, f2 = ff.wtprod_factors(phi, per, pooled)
    assert np.array_equal(p2["A"], [16.0, 0.25, 4.0]) and np.array_equal(p2["B"], [4.0, 1.0, 16.0])
    assert np.array_equal(f2, [4.0, 0.25, 64.0])
    # raised to the step's 0.5: phi_M * per_g ** 0.5, beyond the clip where both move up (4 * 2 = 8)
    assert np.allclose(p2["A"] ** 0.5, [4.0, 0.5, 2.0], rtol=1e-15) and (f2 ** 0.5)[2] == 8.0
    one, onef = ff.wtprod_factors(np.ones(3), per, pooled)  # phi_M = 1 leaves the wild-type-only factor exactly
    assert all(np.array_equal(one[g], per[g]) for g in per) and np.array_equal(onef, pooled)


def test_describe_and_the_pregate():
    rng = np.random.default_rng(4)
    phi = np.clip(np.exp(rng.normal(0, 0.6, 500)), 0.25, 4.0)
    pre = ff.pregate(phi)
    assert pre["dropped"] == {} and pre["n_genes"] == 500
    assert pre["arms"]["top100"]["genes_set_to_1"] == 400 and pre["arms"]["top100"]["genes_changed"] == 400
    assert pre["arms"]["cap90"]["t"] > 0 and pre["arms"]["uni_top"]["s"] < 1
    assert ff.active_arms(pre, k21.ARMS) == k21.ARMS
    d = ff.describe(phi, phi)
    assert d["genes_changed"] == 0 and d["gene_mean_ln"] == pytest.approx(np.log(phi).mean())
    # 100 or more genes at the clip: top100 and its control are dropped, the rest run
    many = phi.copy()
    many[:120] = 4.0
    pre = ff.pregate(many)
    assert set(pre["dropped"]) == {"top100", "uni_top"} and "at the clip" in pre["dropped"]["top100"]
    assert ff.active_arms(pre, k21.ARMS) == ("barred_form", "cap90", "uni_cap", "wtprod", "wtprod_perm")
    # every gene at the clip: the cap cannot be made
    allclip = np.where(np.arange(500) % 2 == 0, 4.0, 0.25)
    pre = ff.pregate(allclip)
    assert {"cap90", "uni_cap", "top100", "uni_top"} == set(pre["dropped"])
    # a factor of 1 but for 6 genes: top100 changes no gene and cap90 at most 6: both pairs are dropped
    few = np.ones(500)
    few[:6] = [2.0, 0.5, 1.5, 3.0, 0.7, 1.2]
    pre = ff.pregate(few)
    assert {"top100", "uni_top", "cap90", "uni_cap"} == set(pre["dropped"])
    assert "ties across the cut" in pre["dropped"]["top100"] and "changes 6 genes" in pre["dropped"]["cap90"]
    # wtprod's printed range, from the wild-type-only factors (never a drop)
    per = {"A": np.full(500, 4.0)}
    pre = ff.pregate(phi, per, np.ones(500))
    assert pre["arms"]["wtprod"]["max"] == pytest.approx(phi.max() * 2.0) and "wtprod" not in pre["dropped"]


# ------------------------------------------------------------------------------------------ the ledger row


def test_the_run_refuses_without_its_row(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text("| K18 (one row) | T3 | `experiments/t3_decomp/k18_built_files.py` |\n")
    with pytest.raises(SystemExit, match="no row"):
        k21.row_present(p)
    p.write_text(p.read_text() + "| K21 (one row) | T3 | `experiments/t3_decomp/k21_factor_forms.py` with ... |\n")
    assert k21.row_present(p) == "K21"
    p.write_text(p.read_text() + "| K22 (one row) | T3 | `experiments/t3_decomp/k21_factor_forms.py` again |\n")
    with pytest.raises(SystemExit, match="more than one"):
        k21.row_present(p)
    for argv in (["run", "--world", "plain", "--task", "g1", "--out", str(tmp_path / "o")],
                 ["pregate", "--out", str(tmp_path / "o")]):
        with pytest.MonkeyPatch.context() as m:
            m.setattr(k21, "TRIES_MD", tmp_path / "none.md")
            (tmp_path / "none.md").write_text("no rows\n")
            with pytest.raises(SystemExit, match="refused"):
                k21.main(argv)
    assert not (tmp_path / "o").exists()


# ------------------------------------------------------------------------------------------------ the arms


@pytest.fixture
def fold(monkeypatch):
    """K18's plain fold fixture: a synthetic pool and reference, the builder's construct, composition and assemble
    stubbed to log what they receive."""
    bs = harness._t3_head()
    rng = np.random.default_rng(5)
    C = rng.poisson(3.0, size=(N_POOL, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, 3, N_POOL)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    V = SimpleNamespace(p=np.array([0.5, 0.3, 0.2]))
    ref = bs.Reference(preset="gata4", rules={}, genes=GENES, V=V, pool=bs.t3d.Pool(C, lab, P, 3),
                       coords=rng.normal(size=(N_POOL, 3)))
    pool = ad.AnnData(X=bs.t3d.lognorm(C).astype(np.float32), var=pd.DataFrame(index=GENES))
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, N_POOL)]
    phi = np.exp(rng.normal(0, 0.5, G))
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {}}
    f_s2 = np.array([0.6, 0.25, 0.15])
    file_genes = GENES[::-1] + ["extra"]
    barred_file = np.clip(np.exp(rng.normal(0, 0.8, G + 1)), 0.25, 4.0)
    spec, _, _ = k18.barred_spec()
    log = {"construct": [], "assemble": [], "move": []}

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        assert r is ref and ratio is None and alphas is None
        rows = g.integers(0, N_POOL, n)  # the rows first, then the factor, as t3_decompose.build
        log["construct"].append({"phi": None if phi is None else np.asarray(phi).copy(), "w": w, "row": rows.copy()})
        Cr = C[rows]
        if phi is not None:
            Cr = bs.t3d.apply_gene_factor(Cr, np.asarray(phi) ** w if w != 1.0 else np.asarray(phi), P[rows], g)
        return bs.Built(X=bs.t3d.lognorm(Cr), lab=lab[rows], source=np.full(n, "wt", dtype=object), row=rows)

    def assemble(r, built, panel, g):
        log["assemble"].append({"row": built.row.copy(), "draw": g.random()})
        return built.X, None, {}

    monkeypatch.setattr(k7, "_prepared", lambda ctx: prep)
    monkeypatch.setattr(bs, "construct", construct)
    monkeypatch.setattr(bs, "composition", lambda v, s, params: (f_s2, {"states_missing": []}))
    monkeypatch.setattr(bs, "assemble", assemble)
    monkeypatch.setattr(harness, "count_for", lambda board, n: N_DRAW)
    monkeypatch.setattr(harness, "finish", lambda ctx, p, X, coords: X)
    monkeypatch.setattr(k18, "_BARRED", {"phi": barred_file, "genes": file_genes, "spec": spec, "record": {}})
    monkeypatch.setattr(ff, "TOP_K", 5)  # 14 genes here
    k18._FOLD.clear()

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    yield SimpleNamespace(bs=bs, ctx=ctx, prep=prep, log=log, pool=pool, C=C, P=P, ref=ref,
                          phi_barred=k18.in_gene_order(barred_file, file_genes, GENES))
    k18._FOLD.clear()


def test_barred_form_is_k18s_arm_value_for_value(fold):
    for seed in (0, 7):
        k18._FOLD.clear()
        want = k18._build_plain(fold.ctx("g1", seed), k18.BARRED)  # K18 as committed
        with k21.k18_arms(k21.ARMS):
            got = k18._build_plain(fold.ctx("g1", seed), k21.BASE)
        assert np.array_equal(got, want)
    # the patch is for the block only
    assert k18.ARMS == k21._K18_ORIGINAL["ARMS"] and k18.build_arm is k21._K18_ORIGINAL["build_arm"]
    assert k18.fold_factors is k21._K18_ORIGINAL["fold_factors"]


def test_every_arm_draws_the_barred_s1s_cells_and_applies_its_form(fold):
    X = {}
    with k21.k18_arms(k21.ARMS):
        for arm in k21.ARMS:
            fold.log["construct"].clear()
            fold.log["assemble"].clear()
            X[arm] = k18._build_plain(fold.ctx("g2", 3), arm)
            (c,) = fold.log["construct"]
            assert c["w"] == 1.0
            if arm in k21.VECTOR_OF:  # the shared transform of phi_M in the reference's order
                want, _ = ff.vector_form(k21.VECTOR_OF[arm], fold.phi_barred)
                assert np.array_equal(c["phi"], want)
            else:
                assert c["phi"] is None
            fold.log[arm] = (c["row"], fold.log["assemble"][0])
    base_rows, base_asm = fold.log[k21.BASE]
    for arm in k21.ARMS:  # one stream: the same rows and the same assemble draw in every arm
        rows, asm = fold.log[arm]
        assert np.array_equal(rows, base_rows) and np.array_equal(asm["row"], base_asm["row"])
        assert asm["draw"] == base_asm["draw"]
        if arm != k21.BASE:
            assert not np.array_equal(X[arm], X[k21.BASE])
    for t, c in k21.CONTROL_OF.items():
        assert not np.array_equal(X[t], X[c])


def test_wtprod_is_the_per_group_step_with_phi_m_squared(fold):
    bs = fold.bs
    with k21.k18_arms(k21.ARMS):
        got = k18._build_plain(fold.ctx("g1", 2), "wtprod")
        gotp = k18._build_plain(fold.ctx("g1", 2), "wtprod_perm")
        F = k18._FOLD["F"]
    built = bs.construct(fold.ref, N_DRAW, bs.t3d.stable_rng(2, k21.STREAM), f=None)
    for name, per, pooled, X in (("wtprod", fold.prep["per"], fold.prep["phi"], got),
                                 ("wtprod_perm", F["per_perm"], F["phi_perm"], gotp)):
        sq = fold.phi_barred ** 2
        want = k16.move_by_group(bs, fold.ref, built, fold.prep["group"], {g: sq * v for g, v in per.items()},
                                 sq * pooled, 0.5, bs.t3d.stable_rng(2, f"{k21.STREAM}:grouped"))
        assert np.array_equal(X, want.X), name
    # the control's factors are K7c's permutation (k16.permuted_factors' rule), not the originals
    basis = harness.t3_mean_expression(fold.pool)
    assert all(np.array_equal(F["per_perm"][g], harness.t3_decile_permutation(v, basis, k7.PERM_SEED))
               for g, v in fold.prep["per"].items())
    assert np.array_equal(F["phi_perm"], harness.t3_decile_permutation(fold.prep["phi"], basis, k7.PERM_SEED))
    assert not np.array_equal(F["per_perm"]["A"], fold.prep["per"]["A"])


def _ref(rng, n=90):
    t3d = harness._t3_head().t3d
    bs = harness._t3_head()
    C = rng.poisson(3.0, size=(n, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    return bs.Reference(preset="gata4", rules={}, genes=GENES, V=None, pool=t3d.Pool(C, rng.integers(0, 3, n), P, 3),
                        coords=rng.normal(size=(n, 3)))


@pytest.mark.parametrize("x", [1.0, 1.5, 2.0])
def test_the_practice_step_equals_the_builders_step(x):
    """K21's wtprod step (k16.move_by_group) and the builder's (the frozen apply_wt_revert) on the same rows, factors
    and generator give the same counts, at every base exponent; at phi_M = 1 both equal the frozen --wt-revert step."""
    bs = harness._t3_head()
    rng = np.random.default_rng(9)
    ref = _ref(rng)
    n = ref.pool.C.shape[0]
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, n)]
    row = rng.integers(0, n, 160)
    built = bs.Built(X=bs.t3d.lognorm(ref.pool.C[row]), lab=ref.pool.lab[row], source=np.full(160, "wt", dtype=object),
                     row=row)
    phi_m = np.clip(np.exp(rng.normal(0, 0.7, G)), 0.25, 4.0)
    per = {"A": np.exp(rng.normal(0, 0.4, G)), "B": np.exp(rng.normal(0, 0.4, G))}
    pooled = np.exp(rng.normal(0, 0.4, G))
    phi_mx = ff.at_exponent(phi_m, x)
    p2, f2 = ff.wtprod_factors(phi_mx, per, pooled)
    practice = k16.move_by_group(bs, ref, built, group, p2, f2, k7.EXPONENT, np.random.default_rng(31))
    build, _ = bs.apply_wt_revert(ref, built, f2, p2, group, np.random.default_rng(31))
    assert np.array_equal(practice.X, build.X) and not np.array_equal(practice.X, built.X)
    p1, f1 = ff.wtprod_factors(np.ones(G), per, pooled)
    a = k16.move_by_group(bs, ref, built, group, p1, f1, k7.EXPONENT, np.random.default_rng(31))
    b, _ = bs.apply_wt_revert(ref, built, pooled, per, group, np.random.default_rng(31))
    assert np.array_equal(a.X, b.X)


# --------------------------------------------------------------------------------------- the depth-matched row


def test_a_depth_matched_row_builds_this_trys_arms_with_k18s_setup(monkeypatch):
    """k18.dm_row inside k18_arms: every arm of this try and K18's check, the arms all on the one stream and the
    matched pool; the estimate, the scorer and de_genes stubbed."""
    from common import core_metrics

    from vec import evaluate

    bs = harness._t3_head()
    t3d = bs.t3d
    rng = np.random.default_rng(3)
    L, n_wt, n_ref, n_t = 3, 120, 80, 60
    C = rng.poisson(20.0, size=(n_wt, G)).astype(np.int32)
    lab = rng.integers(0, L, n_wt)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    V = SimpleNamespace(p=np.array([0.5, 0.3, 0.2]), L=L)
    rules = {"protocol": {"working_copy_fraction": 0.5, "split_seed": 0, "prediction_cells_multiplier": 2},
             "components": {"stage": {"min_cells": 5}}}
    T = SimpleNamespace(name="rep2", X_rel=rng.normal(size=(n_t, G)).astype(np.float32),
                        hash_rel=np.arange(n_t), hash_pool=np.arange(n_t + 20))
    ref_C = rng.poisson(20.0, size=(n_ref, G)).astype(np.int32)
    ctx = SimpleNamespace(targets=[T], rules=rules, V=V, ref_X=t3d.lognorm(ref_C), preset_name="gata4", genes=GENES,
                          wt=SimpleNamespace(pool=t3d.Pool(C, lab, P, L)), stages=[])
    group = np.array(["A", "B"])[rng.integers(0, 2, n_wt)]
    per = {"A": np.exp(rng.normal(0, 0.3, G)), "B": np.exp(rng.normal(0, 0.3, G))}
    phi = np.exp(rng.normal(0, 0.3, G))
    spec, _, _ = k18.barred_spec()
    monkeypatch.setattr(bs, "composition", lambda v, s, params: (np.array([0.6, 0.25, 0.15]), {"states_missing": []}))
    monkeypatch.setattr(ff, "TOP_K", 5)
    refw = bs.Reference(preset="gata4", rules=rules, genes=GENES, V=V, pool=ctx.wt.pool, coords=np.zeros((n_wt, 3)))
    phi_b = np.clip(np.exp(rng.normal(0, 0.8, G)), 0.25, 4.0)
    monkeypatch.setattr(t3d, "estimate", lambda c, t, rows, g: SimpleNamespace(ratio=np.full(L, 0.5)))
    monkeypatch.setattr(evaluate, "panel_metrics", lambda *a, **k: {})
    monkeypatch.setattr(core_metrics, "de_genes", lambda a, b: (np.arange(5), np.arange(5, 10), None))
    monkeypatch.setattr(k18, "_scored", lambda X, truth, ref, anchors, up, dn, seed: {
        "score": 50.0, "skills": {}, "severity": {"state": "scored", "beta": 0.1, "r2": 0.1}})
    real_build = t3d.build
    streams = []

    def build(c, n, g, f, src=None, **kw):
        streams.append((json.dumps(g.bit_generator.state, default=int), c.wt.pool if src is None else src))
        return real_build(c, n, g, f, src=src, **kw)

    monkeypatch.setattr(t3d, "build", build)
    with k21.k18_arms(k21.ARMS):
        F = k18.fold_factors(bs, refw, per, phi, group, ctx.ref_X.mean(0), phi_b, spec)
        st = {"ctx": ctx, "phi_b": phi_b, "ref_counts": ref_C, "ref_lab": rng.integers(0, L, n_ref),
              "ref_props": (ref_C / ref_C.sum(1, keepdims=True)), "F": F, "threads": 1}
        row = k18.dm_row(st, 0, 4)
    assert row["target"] == "rep2" and set(row["scores"]) == set(k21.ARMS) and "dm_B" in row["check"]
    arm_state = json.dumps(t3d.stable_rng(4, k21.STREAM).bit_generator.state, default=int)
    assert len(streams) == 1 + len(k21.ARMS) and all(s[0] == arm_state for s in streams[1:])
    assert {id(s[1]) for s in streams} == {id(streams[0][1])} and streams[0][1] is not ctx.wt.pool


# ----------------------------------------------------------------------------------------------- the rule


def _rows(scores: dict, gated: dict | None = None, noise: dict | None = None, skill_shift=None) -> list[dict]:
    gated = gated or {}
    out = []
    for s in SEEDS:
        sc = {a: v + 0.1 * s + (noise[a][s] if noise and a in noise else 0.0) for a, v in scores.items()}
        out.append({"seed": s, "scores": sc,
                    "skills": {a: {"de_score": 0.5 + (skill_shift or {}).get(a, 0.0), "severity_slope": 0.6}
                               for a in sc},
                    "t3": {a: {"severity": {"state": "r2_gate" if s < gated.get(a, 0) else "scored", "beta": 0.4,
                                            "r2": 0.05}} for a in sc}})
    return out


BASE_SCORE = 56.0


def _table(top=None, cap=None, wt=None, ctl=None) -> dict:
    t = {harness.FLOOR: 50.0, k21.BASE: BASE_SCORE, "top100": BASE_SCORE, "uni_top": BASE_SCORE,
         "cap90": BASE_SCORE, "uni_cap": BASE_SCORE, "wtprod": BASE_SCORE, "wtprod_perm": BASE_SCORE}
    for k, v in (("top100", top), ("cap90", cap), ("wtprod", wt)):
        if v is not None:
            t[k] = v
    t.update(ctl or {})
    return t


def _per(cells: dict, gated: dict | None = None, arms=k21.ARMS) -> dict:
    out: dict = {}
    for w in k21.WORLDS:
        a = (harness.FLOOR, *arms) if w == "plain" else arms
        out[w] = {n: k21.summarise(_rows({x: cells.get((w, n), _table())[x] for x in a},
                                         (gated or {}).get((w, n))), a) for n in k21.TASKS}
    return out


def test_qualified_needs_both_worlds_both_replicates_no_loss_the_gates_and_the_control():
    up = _table(cap=57.0)
    one_each = {("plain", "g1"): up, ("dm", "g2"): up}
    res = k21.decide(_per(one_each), k21.ARMS)
    assert res["status"]["cap90"]["qualified"] and res["names"] == ["cap90"]
    assert not res["status"]["top100"]["qualified"]
    assert res["status"]["cap90"]["ahead_cells"] == ["plain g1", "dm g2"]
    # one world only, or one replicate only: not qualified
    for cells in ({("plain", "g1"): up, ("plain", "g2"): up}, {("plain", "g1"): up, ("dm", "g1"): up}):
        assert not k21.decide(_per(cells), k21.ARMS)["status"]["cap90"]["qualified"]
    # behind in a cell
    res = k21.decide(_per(one_each | {("dm", "g1"): _table(cap=55.5)}), k21.ARMS)
    assert not res["status"]["cap90"]["qualified"] and res["status"]["cap90"]["behind_cells"] == ["dm g1"]
    # its gated seeds two above the base's in one cell
    res = k21.decide(_per(one_each, {("dm", "g1"): {"cap90": 2}}), k21.ARMS)
    assert not res["status"]["cap90"]["c"] and not res["names"]
    assert k21.decide(_per(one_each, {("dm", "g1"): {"cap90": 1}}), k21.ARMS)["names"] == ["cap90"]
    # its control as high where it is ahead: the size alone, not the shape
    ctl = _table(cap=57.0, ctl={"uni_cap": 56.9})
    res = k21.decide(_per({("plain", "g1"): ctl, ("dm", "g2"): up}), k21.ARMS)
    assert not res["status"]["cap90"]["d"] and not res["names"]
    # a control never qualifies or is named, however high
    res = k21.decide(_per({k: _table(ctl={"uni_top": 60.0, "uni_cap": 60.0, "wtprod_perm": 60.0})
                           for k in [(w, n) for w in k21.WORLDS for n in k21.TASKS]}), k21.ARMS)
    assert res["names"] == [] and set(res["status"]) == set(k21.NAMEABLE)


def test_the_ranking_and_at_most_two_names():
    all_up = {(w, n): _table(top=57.0, cap=58.0, wt=56.5) for w in k21.WORLDS for n in k21.TASKS}
    res = k21.decide(_per(all_up), k21.ARMS)
    assert res["ranked"] == ["cap90", "top100", "wtprod"] and res["names"] == ["cap90", "top100"]
    assert "file 1: cap90 on base B" in res["outcome"] and "file 2: top100 on base B" in res["outcome"]
    # more AHEAD cells rank first, whatever the mean
    mix = {(w, n): _table(top=56.5, cap=60.0) for w in k21.WORLDS for n in k21.TASKS}
    mix[("plain", "g2")] = _table(top=56.5, cap=56.1)
    res = k21.decide(_per(mix), k21.ARMS)
    assert res["ranked"][:2] == ["top100", "cap90"]
    # ties go to the fixed order
    tie = {(w, n): _table(top=57.0, cap=57.0) for w in k21.WORLDS for n in k21.TASKS}
    assert k21.decide(_per(tie), k21.ARMS)["names"] == ["top100", "cap90"]
    # a dropped arm is not run and cannot be named
    arms = ("barred_form", "cap90", "uni_cap", "wtprod", "wtprod_perm")
    res = k21.decide(_per(tie, arms=arms), arms)
    assert res["status"]["top100"] == {"ran": False, "qualified": False} and res["names"] == ["cap90"]


def test_on_pure_noise_a_name_is_rare():
    """Every arm = barred_form + noise of sd 1 per seed, no true gain; the plain world's g1 and g2 share their noise
    (one pool, one stream, one seed), the worst case. A false name costs one slot and no points."""
    rng = np.random.default_rng(12)
    named = 0
    for _ in range(300):
        per = {}
        for w in k21.WORLDS:
            per[w] = {}
            common = {a: rng.normal(0, 1.0, len(SEEDS)) for a in k21.ARMS}
            for n in k21.TASKS:
                noise = {a: common[a] if w == "plain" else rng.normal(0, 1.0, len(SEEDS)) for a in k21.ARMS
                         if a != k21.BASE}
                per[w][n] = k21.summarise(_rows({a: BASE_SCORE for a in k21.ARMS}, noise=noise), k21.ARMS)
        named += bool(k21.decide(per, k21.ARMS)["names"])
    assert named <= 6


def test_the_per_term_breakdown():
    rows = _rows({k21.BASE: 56.0, "top100": 57.0}, skill_shift={"top100": 0.02})
    s = k21.summarise(rows, (k21.BASE, "top100"))
    assert s["per_term"]["top100"]["de_score"] == pytest.approx(0.02) and s["per_term"]["top100"]["severity_slope"] == 0
    assert s["paired"]["top100 - barred_form"]["mean"] == pytest.approx(1.0)


# ------------------------------------------------------------------------------------ the reading of a folder


def _prov(commit="abc", clean=True, machine="m1"):
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "x"}


def _write(root: Path, cells=None, commit="abc", clean=True, machine="m1", shift=0.0, k18_dir=None, drop=None,
           arms=k21.ARMS, bad_arms_file=None):
    cells = cells or {}
    pre = ff.pregate(np.clip(np.exp(np.random.default_rng(0).normal(0, 0.6, 500)), 0.25, 4.0))
    if arms != k21.ARMS:
        pre["dropped"] = {a: "test" for a in k21.ARMS if a not in arms}
    root.mkdir(parents=True)
    rec = {"arms": list(arms), "dropped": pre["dropped"]}
    k18.write_task_file(root / "pregate.json", {"pregate": pre, "arms": list(arms), "provenance": _prov(commit, clean,
                                                                                                         machine)})
    refs = {}
    for w in k21.WORLDS:
        for n in k21.TASKS:
            a = (harness.FLOOR, *arms) if w == "plain" else arms
            rows = _rows({x: cells.get((w, n), _table())[x] for x in a})
            refs[(w, n)] = [{"seed": r["seed"], "scores": {k21.BASE: r["scores"][k21.BASE] - shift}} for r in rows]
            for r in rows:
                s = r["seed"]
                f = root / (f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json")
                f.parent.mkdir(parents=True, exist_ok=True)
                if drop == str(f.relative_to(root)):
                    continue
                r2 = {**r, "check": {k18.DM_CHECK: 60.0 + s}}
                body = ({"task": {"key": k21.TASKS[n].key}, "rows": [r], "seed": s, "pregate": rec,
                         "provenance": _prov(commit, clean, machine)} if w == "plain" else
                        r2 | {"target": k18.DM_TARGETS[n], "task": n, "pregate": rec, "arms": list(arms),
                              "provenance": _prov(commit, clean, machine)})
                if bad_arms_file == str(f.relative_to(root)):
                    body["pregate"] = {"arms": list(arms[:-1]), "dropped": {}}
                k18.write_task_file(f, body)
    k18_dir = k18_dir or root.parent / f"k18_{root.name}"
    k18_dir.mkdir(parents=True, exist_ok=True)
    for n in k21.TASKS:
        (k18_dir / f"plain_{n}.json").write_text(json.dumps({"rows": refs[("plain", n)]}))
    (k18_dir / "dm.json").write_text(json.dumps({"rows": {n: refs[("dm", n)] for n in k21.TASKS}}))
    pins = {f: hashlib.sha256((k18_dir / f).read_bytes()).hexdigest() for f in k21.K18_SHA256}
    return root, k18_dir, pins


@pytest.fixture
def u3d(tmp_path):
    rows = [{"target": t, "seed": s, "world": "depth_matched", "variant": "dm:transfer:mab21l2_all", "t3": 60.0 + s}
            for t in ("rep1", "rep2") for s in SEEDS]
    p = tmp_path / "orth_records.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_the_reading_of_a_run_folder(tmp_path, u3d, monkeypatch):
    up = _table(cap=57.0)
    root, refs, pins = _write(tmp_path / "ok", {("plain", "g1"): up, ("dm", "g2"): up})
    monkeypatch.setattr(k21, "K18_SHA256", pins)
    res = k21.reading(root, refs, u3d)
    assert not res["void"] and res["names"] == ["cap90"] and res["tries_counted"] == 6
    assert all(d["max_abs_diff"] == 0.0 and d["seeds_off"] == [] for d in res["guard"].values())
    assert len(res["files_sha256"]) == 41
    txt = k21.text(res)
    assert "WORLD plain" in txt and "WORLD dm" in txt and "file 1: cap90 on base B" in txt
    assert "per-term skill, cap90 - barred_form" in txt and "barred_form - floor" in txt and "pre-gate" in txt


@pytest.mark.parametrize("kw,why", [({"shift": 0.2}, "guard"), ({"clean": False}, "dirty tree"),
                                    ({"bad_arms_file": "dm/g1_seed_3.json"}, "arms or drops differ")])
def test_void_readings(tmp_path, u3d, monkeypatch, kw, why):
    root, refs, pins = _write(tmp_path / "v", **kw)
    monkeypatch.setattr(k21, "K18_SHA256", pins)
    res = k21.reading(root, refs, u3d)
    assert res["void"] and any(why in w for w in res["why"]) and "no file" in res["outcome"]
    assert "VOID" in k21.text(res)


def test_mixed_commits_or_machines_or_a_changed_k18_file_void(tmp_path, u3d, monkeypatch):
    root, refs, pins = _write(tmp_path / "m")
    monkeypatch.setattr(k21, "K18_SHA256", pins)
    p = root / "dm" / "g2_seed_5.json"
    body = json.loads(p.read_text())
    body["provenance"]["machine"] = "m2"
    p.write_text(json.dumps(body))
    res = k21.reading(root, refs, u3d)
    assert res["void"] and any("different machines" in w for w in res["why"])
    monkeypatch.setattr(k21, "K18_SHA256", {**pins, "dm.json": "0" * 64})
    root2, refs2, _ = _write(tmp_path / "k")
    res = k21.reading(root2, refs2, u3d)
    assert res["void"] and any("pinned sha256" in w for w in res["why"])


def test_a_dropped_arm_reads_with_fewer_tries(tmp_path, u3d, monkeypatch):
    arms = ("barred_form", "cap90", "uni_cap", "wtprod", "wtprod_perm")
    root, refs, pins = _write(tmp_path / "d", arms=arms)
    monkeypatch.setattr(k21, "K18_SHA256", pins)
    res = k21.reading(root, refs, u3d)
    assert not res["void"] and res["tries_counted"] == 4 and res["arms"] == list(arms)
    assert res["status"]["top100"]["ran"] is False


def test_a_missing_file_is_not_a_reading(tmp_path, u3d, monkeypatch, capsys):
    root, refs, pins = _write(tmp_path / "i", drop="plain_g2/seed_7.json")
    monkeypatch.setattr(k21, "K18_SHA256", pins)
    (root / "dm" / "g1_seed_0.json").write_text("{ not JSON: never opened")
    res = k21.reading(root, refs, u3d)
    assert res["incomplete"] and res["missing"] == ["plain_g2/seed_7.json"]
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k21.main(["read", "--runs", str(root), "--out", str(out)])
    assert e.value.code == 3 and not out.exists() and "not a reading" in capsys.readouterr().out


# ------------------------------------------------------------------------------------------- the plain run


@pytest.fixture
def fake_run(monkeypatch, tmp_path):
    calls = {"seeds": [], "log": None}
    lines = {}
    tries = tmp_path / "TRIES.md"
    tries.write_text("| K21 (one row) | T3 | `experiments/t3_decomp/k21_factor_forms.py` |\n")

    def run_task(t, cands, seeds, data_dir, log=print):
        calls["seeds"] += list(seeds)
        calls["log"] = log
        calls["names"] = [c.name for c in cands]
        for s in seeds:
            lines[s] = f"{t.key} seed {s}: floor 50.00, barred_form {60 - s / 9:.2f}, top100 {61 + s / 7:.2f}"
            log(lines[s])
        return {"task": {"key": t.key}, "rows": [{"seed": s} for s in seeds],
                "provenance": {"commit": "abc", "tree_clean": True}}

    phi = np.clip(np.exp(np.random.default_rng(1).normal(0, 0.6, 500)), 0.25, 4.0)

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=phi, genes=[f"g{i}" for i in range(500)], spec={}, record={})
        return {}

    monkeypatch.setattr(k21, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(k7, "_PREP", {"fold": {"record": {}}})
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"})
    return SimpleNamespace(calls=calls, lines=lines)


def test_the_plain_run_is_seed_by_seed_resumable_and_logs_progress_only(tmp_path, fake_run, capsys):
    out = tmp_path / "run"
    k21.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0-1", "--out", str(out)])
    assert fake_run.calls["seeds"] == [0, 1] and fake_run.calls["names"] == [harness.FLOOR, *k21.ARMS]
    printed = capsys.readouterr()
    for x in printed.out.splitlines():
        assert re.fullmatch(r"K21 plain g1 seed \d: scored|wrote \S+/seed_\d\.json: commit abc, clean True, "
                            r"machine m1", x), x
    for line in fake_run.lines.values():
        for score in re.findall(r"\d+\.\d\d", line):
            assert score not in printed.out and score not in printed.err
    assert sorted(p.name for p in (out / "plain_g1").iterdir()) == [
        "seed_0.json", "seed_0.json.prov.json", "seed_1.json", "seed_1.json.prov.json"]
    body = json.loads((out / "plain_g1" / "seed_0.json").read_text())
    assert body["pregate"]["arms"] == list(k21.ARMS) and body["provenance"]["machine"] == "m1"
    k21.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0-2", "--out", str(out)])  # resumes
    assert fake_run.calls["seeds"] == [0, 1, 2]
    assert "seed 0: present, skipped" in capsys.readouterr().out
    # back in place after the run
    assert k18.ARMS == k21._K18_ORIGINAL["ARMS"] and k18.build_arm is k21._K18_ORIGINAL["build_arm"]


def test_a_run_refuses_another_commit_a_dirty_tree_and_bad_seeds(tmp_path, fake_run, monkeypatch, capsys):
    out = tmp_path / "run"
    k21.main(["run", "--world", "plain", "--task", "g2", "--seeds", "3", "--out", str(out)])
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "def", "tree_clean": True, "machine": "m1", "utc": "x"})
    with pytest.raises(SystemExit):
        k21.main(["run", "--world", "plain", "--task", "g2", "--seeds", "4", "--out", str(out)])
    assert not (out / "plain_g2" / "seed_4.json").exists()
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": False, "machine": "m1", "utc": "x"})
    with pytest.raises(SystemExit, match="not clean"):
        k21.main(["run", "--world", "plain", "--task", "g2", "--seeds", "4", "--out", str(out)])
    for bad in ("10", "0-0,0", "a"):
        with pytest.raises((SystemExit, ValueError)):
            k21.parse_seeds(bad)
    assert k21.parse_seeds(None) == list(SEEDS) and k21.parse_seeds("2-4,7") == [2, 3, 4, 7]


def test_a_crash_is_masked_and_writes_nothing(tmp_path, fake_run, monkeypatch, capsys):
    def boom(*a, **k):
        raise RuntimeError("score 61.234 went wrong")

    monkeypatch.setattr(harness, "run_task", boom)
    with pytest.raises(SystemExit) as e:
        k21.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0", "--out", str(tmp_path / "r")])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "61.234" not in err and "<num>" in err
    assert not (tmp_path / "r" / "plain_g1" / "seed_0.json").exists()
    assert k18.ARMS == k21._K18_ORIGINAL["ARMS"]
