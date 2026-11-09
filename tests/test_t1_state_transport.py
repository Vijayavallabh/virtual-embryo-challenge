"""vec.t1_state_transport, O1's arms in experiments/t1_scope/pseudo_task.py, its readout and candidate g_mom (tries T8
to T11, experiments/harness/TRIES.md), on synthetic cells."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_t1_scope import E85_SIZES, E95_SIZES, GENES, HELD, _profiles

from vec import harness
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT, Board
from vec.evaluate import Stage
from vec.harness import CANDIDATES, PseudoTask
from vec.prep import log1p_cp10k
from vec.t1_scope import build_submission, classify, load_spec, reweight


def states_stage(seed, sizes, later=False, n_scale=1.0, target=False):
    """A T1-like stage where half of the V-CM cells sit in a second state (g30-g39 boosted 8x); with later=True that
    state has also moved within its type (g40-g44 boosted another 6x), so that a per-type mean sees half of it. With
    target=True every cell's g0-g7 are 4x up and g8-g15 4x down, so a scorer sees a real change from other stages."""
    rng = np.random.default_rng(seed)
    prof = _profiles()
    blocks, labels = [], []
    for lab, n0 in sizes.items():
        n = int(n0 * n_scale)
        mu = np.tile(prof[lab], (n, 1)) * rng.gamma(10.0, 0.1, size=(n, 1))
        if lab == "V-CM":
            second = np.arange(n) < n // 2
            mu[second, 30:40] *= 8.0
            if later:
                mu[second, 40:45] *= 6.0
        if target:
            mu[:, :8] *= 4.0
            mu[:, 8:16] *= 0.25
        C = rng.poisson(mu)
        C[np.arange(n), rng.integers(0, C.shape[1], n)] += 1
        blocks.append(C)
        labels += [lab] * n
    C = sparse.csr_matrix(np.vstack(blocks).astype(np.int32))
    a = ad.AnnData(X=log1p_cp10k(C))
    a.var_names = GENES
    a.obs_names = [f"s{seed}_{i}" for i in range(a.n_obs)]
    a.obs["celltype"] = labels
    a.layers["counts"] = C
    return a


def classified_pair(seed=1, n_scale=1.0):
    early = states_stage(seed, E85_SIZES, n_scale=n_scale)
    late = states_stage(seed + 1, E95_SIZES, later=True, n_scale=n_scale)
    classify(early, source="E8.5_RNA")
    classify(late, source="E9.5_RNA")
    return early, late


FALLBACKS = dict(earlier_fallback=load_spec().fallbacks("E8.5_RNA"), later_fallback=load_spec().fallbacks("E9.5_RNA"))


# ------------------------------------------------------------------------------------------------ building blocks
def test_hvg_columns_keep_detected_dispersed_genes_in_order():
    rng = np.random.default_rng(0)
    X = sparse.csr_matrix(rng.poisson(0.3, size=(400, 50)).astype(float))
    X = X.tolil()
    X[:, 5] = 0  # undetected, never chosen
    X = sparse.csr_matrix(X)
    cols = ts.hvg_columns(X, n_top=10)
    assert list(cols) == sorted(cols) and len(cols) == 10 and 5 not in cols
    assert len(ts.hvg_columns(X, n_top=500)) == X.shape[1] - 1  # fewer genes pass than asked for: only those


def test_joint_pcs_returns_each_stage_its_rows():
    early, late = classified_pair()
    pa, pb, rec = ts.joint_pcs(early, late, n_hvg=40, n_pcs=6)
    assert pa.shape == (early.n_obs, 6) and pb.shape == (late.n_obs, 6)
    assert rec == {"n_hvg": 40, "n_pcs": 6, "pca_seed": 0}


def test_sinkhorn_plan_has_the_uniform_marginals_at_the_tolerance():
    rng = np.random.default_rng(1)
    xa, xb = rng.normal(size=(120, 5)), rng.normal(size=(90, 5)) + 0.5
    plan, rec = ts.sinkhorn_plan(xa, xb)
    assert plan.shape == (120, 90) and abs(plan.sum() - 1) < 1e-6
    assert np.abs(plan.sum(1) - 1 / 120).max() < 1e-5 and np.abs(plan.sum(0) - 1 / 90).max() < 1e-5
    assert rec["converged"] and rec["median_cost"] > 0
    with pytest.raises(ValueError, match="no spread"):
        ts.sinkhorn_plan(np.zeros((3, 2)), np.zeros((4, 2)))


@pytest.mark.parametrize("n, k", [(30, 5), (300, 5), (760, 5), (1500, 10), (9000, 40)])
def test_metacell_count_is_the_row_rule(n, k):
    x = np.random.default_rng(3).normal(size=(n, 3))
    labels = ts.metacell_labels(x)
    assert labels.shape == (n,) and labels.min() == 0 and len(np.unique(labels)) == labels.max() + 1 == k


def test_metacells_handle_fewer_cells_than_the_minimum():
    x = np.random.default_rng(3).normal(size=(3, 2))
    assert ts.metacell_labels(x).max() == 2


# ------------------------------------------------------------------------------------------------ the factors
def test_fit_finds_the_within_type_change_and_skips_stand_ins():
    early, late = classified_pair()
    sf = ts.fit_state_factors(early, late, "earlier", ["cm", "vascular_endothelium", "epicardium", "surrounding"],
                              within={"surrounding": ["Foregut"]}, n_hvg=60, n_pcs=8, **FALLBACKS)
    assert set(sf.categories) == {"cm", "surrounding"}
    assert "stand-in" in sf.skipped["vascular_endothelium"] and "stand-in" in sf.skipped["epicardium"]
    cm = sf.categories["cm"]
    assert cm.n_metacells == int(np.clip(E85_SIZES["V-CM"] // 150, 5, 40)) and set(cm.index) <= set(early.obs_names)
    moved = cm.logf[:, 40:45].mean(1)
    assert moved.max() > 0.6 and moved.min() < 0.3  # the second state's metacells moved, the first state's did not
    assert cm.record["converged"] and cm.record["n_cells_a"] == E85_SIZES["V-CM"]
    assert cm.logf.shape == (cm.n_metacells, len(GENES))


def test_fit_on_the_later_stage_keeps_its_names_and_the_same_sign():
    early, late = classified_pair()
    sf = ts.fit_state_factors(early, late, "later", ["cm"], n_hvg=60, n_pcs=8, **FALLBACKS)
    cm = sf.categories["cm"]
    assert set(cm.index) <= set(late.obs_names) and cm.logf[:, 40:45].mean(1).max() > 0.6


def test_fit_refuses_an_unknown_side_and_mismatched_genes():
    early, late = classified_pair()
    with pytest.raises(ValueError, match="metacells_on"):
        ts.fit_state_factors(early, late, "middle", ["cm"])
    with pytest.raises(ValueError, match="different genes"):
        ts.fit_state_factors(early, late[:, GENES[::-1]].copy(), "later", ["cm"])


def test_fit_skips_a_category_with_too_few_cells():
    early, late = classified_pair()
    keep = np.flatnonzero(early.obs["category"].astype(str).to_numpy() == "cm")[:5]
    other = np.flatnonzero(early.obs["category"].astype(str).to_numpy() != "cm")
    small = early[np.concatenate([keep, other])].copy()
    sf = ts.fit_state_factors(small, late, "earlier", ["cm"], n_hvg=60, n_pcs=8, **FALLBACKS)
    assert not sf.categories and "fewer than 20" in sf.skipped["cm"]


def _factors(logf, sizes=(1, 3)):
    cf = ts.CategoryFactors("cm", {"a": 0, "b": 1, "c": 1, "d": 1}, np.array(sizes), np.array(logf, dtype=np.float64))
    return ts.StateFactors(genes=["g0", "g1"], metacells_on="earlier", categories={"cm": cf}, skipped={})


def test_finish_oracle_zeroes_below_the_bar_and_clips_without_damping():
    out = ts.finish_oracle(_factors([[0.04, -0.03], [2.0, -1.5]])).categories["cm"].logf
    assert out[0].tolist() == [0.0, 0.0]  # below log 1.05 in absolute value: exactly no change
    assert out[1] == pytest.approx([np.log(3.0), -np.log(3.0)])  # beyond log 3: clipped, not dropped


def test_finish_momentum_removes_the_weighted_mean_then_damps_then_zeroes_and_clips():
    sf = _factors([[0.3, 3.0], [-0.1, -1.0]])  # weighted means (sizes 1, 3) are 0 and 0
    out = ts.finish_momentum(sf, 0.5).categories["cm"].logf
    assert out[0] == pytest.approx([0.15, np.log(3.0)])  # damped; 1.5 is above log 3, so it is clipped
    assert out[1] == pytest.approx([-0.05, -0.5])
    small = ts.finish_momentum(_factors([[0.06, 0.0], [-0.02, 0.0]]), 0.5).categories["cm"].logf
    assert small.tolist() == [[0.0, 0.0], [0.0, 0.0]]  # after damping both are below log 1.05, so the zeroing sees them
    shifted = ts.finish_momentum(_factors([[1.3, 3.0], [0.9, -1.0]]), 0.5).categories["cm"].logf
    assert shifted == pytest.approx(ts.finish_momentum(sf, 0.5).categories["cm"].logf)  # a constant mean is removed


def test_r1_labels_match_the_builds_label_groups():
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_m1"))
    import feng_dirs

    labels = ["V-CM", "A-CM", "OFT-CM", "IFT-CM", "AVC-CM", "SV-CM", "Trabecular CM", "Endocardium", "BEC",
              "Endothelium", "Proepicardium", "Epicardium", "EPDC", "Blood", "Macrophage", "Early VIC", "aSHF", "JCF",
              "Foregut", "ST"]
    groups = feng_dirs.label_groups(labels)
    want = {lab for lab in labels if groups[lab] not in ts.R1_FENG_GROUPS}
    assert ts.r1_unbracketed_labels(labels) == want
    assert {"Blood", "Macrophage", "Early VIC"} <= want and not {"V-CM", "Endocardium", "BEC", "EPDC"} & want
    committed = json.loads((REPO_ROOT / "experiments" / "t1_hygiene" / "results" / "n12_hygiene_size" /
                            "n12_hygiene_size.json").read_text())
    assert set(ts.R1_FENG_GROUPS) == set(committed["groups"])  # the groups with a Feng direction, as N12 measured them


# ------------------------------------------------------------------------------------------------ applying them
def _toy_population(seed=0):
    """A small population with counts, categories and source names, and factor rows for two categories."""
    rng = np.random.default_rng(seed)
    n = 300
    C = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 2.0, size=(n, 30))).astype(np.int64))
    C[(C.sum(1).A1 == 0).nonzero()[0], 0] = 1
    C = sparse.csr_matrix(C)
    genes = [f"g{i}" for i in range(30)]
    cats = np.array(["cm"] * 200 + ["blood"] * 100, dtype=object)
    names = [f"cell{i}" for i in range(n)]
    pop = ad.AnnData(X=log1p_cp10k(C), obs=pd.DataFrame({"category": cats, "source_obs": names}, index=names),
                     var=pd.DataFrame(index=genes))
    pop.layers["counts"] = C
    meta_rng = np.random.default_rng(1)
    idx = {"cm": {names[i]: int(meta_rng.integers(4)) for i in range(200)},
           "blood": {names[i]: int(meta_rng.integers(2)) for i in range(200, 300)}}
    sizes = {"cm": np.bincount(list(idx["cm"].values()), minlength=4), "blood": np.bincount(list(idx["blood"].values()),
                                                                                             minlength=2)}
    logf = {"cm": meta_rng.normal(0.0, 0.5, size=(4, 30)), "blood": meta_rng.normal(0.0, 0.4, size=(2, 30))}
    cats_out = {c: ts.CategoryFactors(c, idx[c], sizes[c], logf[c]) for c in ("cm", "blood")}
    return pop, ts.StateFactors(genes=genes, metacells_on="earlier", categories=cats_out, skipped={})


def _cell_scale(C, F, row, held):
    C = sparse.csr_matrix(C, dtype=np.float64)
    entries = np.repeat(np.arange(C.shape[0]), np.diff(C.indptr))
    free = ~held[C.indices]
    f = F[row[entries], C.indices]
    num = np.bincount(entries, weights=C.data * free, minlength=C.shape[0])
    den = np.bincount(entries, weights=C.data * f * free, minlength=C.shape[0])
    return np.divide(num, den, out=np.ones(C.shape[0]), where=den > 0)


def test_rescale_holds_each_categorys_pseudobulk_in_expectation_with_the_depth_constant_included():
    pop, sf = _toy_population()
    C = sparse.csr_matrix(pop.layers["counts"], dtype=np.int64)
    held = np.array([i in (3, 7) for i in range(30)], dtype=bool)
    row, logf, blocks = ts._factor_rows(pop, sf, "category", "source_obs")
    F = np.exp(logf)
    F[:, held] = 1.0
    F0 = F.copy()
    rescaled, rec = ts.rescale_to_hold_pseudobulk(C, F.copy(), row, blocks, held)
    assert not np.allclose(rescaled, F0) and rec["cm"]["iterations"] >= 1
    s = _cell_scale(C, rescaled, row, held)
    depth = np.asarray(C.sum(1), dtype=np.float64).ravel()
    W = sparse.diags(1.0 / depth) @ sparse.csr_matrix(C, dtype=np.float64)
    for cat in ("cm", "blood"):
        offset, k, rows = blocks[cat]
        W_c = W[rows]
        target = np.asarray(W_c.sum(0)).ravel()
        entry_rows = np.repeat(np.arange(len(rows)), np.diff(W_c.indptr))
        value = W_c.data * rescaled[row[rows][entry_rows], W_c.indices] * s[rows][entry_rows]
        got = np.bincount(W_c.indices, weights=value, minlength=30)
        free = target > 0
        free[held] = False
        assert np.abs(got[free] / target[free] - 1.0).max() < 1e-6  # the expected mean CP10k of every free gene is held
        assert np.all(rescaled[offset:offset + k][:, held] == 1.0)


def test_apply_with_rescale_keeps_the_category_means_and_the_held_genes():
    early, late = classified_pair(seed=4, n_scale=1.5)
    sf = ts.finish_oracle(ts.fit_state_factors(early, late, "earlier", ["cm"], n_hvg=60, n_pcs=8, **FALLBACKS))
    shares = {"cm": 0.7, "surrounding": 0.3}
    pop, _ = reweight(early, shares, 800, 0, within={"surrounding": ["Foregut"]})
    out, rec = ts.apply_state_factors(pop, sf, HELD, [0, 11], rescale=True)
    before, after = sparse.csr_matrix(pop.layers["counts"]), sparse.csr_matrix(out.layers["counts"])
    held_cols = [GENES.index(g) for g in HELD]
    assert (before[:, held_cols] != after[:, held_cols]).nnz == 0  # the held-flat genes are the original counts
    cats = np.asarray(pop.obs["category"]).astype(str)
    assert np.array_equal(before[cats != "cm"].toarray(), after[cats != "cm"].toarray())  # no factors: unchanged
    m0, m1 = ts.category_mean_cp10k(before, cats, "cm"), ts.category_mean_cp10k(after, cats, "cm")
    well = m0 > 50
    assert well.sum() > 10 and np.abs(m1[well] / m0[well] - 1).max() < 0.05  # sampling noise only; the mean is held
    assert rec["n_edited"] == int((cats == "cm").sum()) and rec["rescale"]["cm"]["iterations"] >= 1
    assert after.sum() == pytest.approx(before.sum(), rel=0.02)  # the keep-depth constant keeps the total


def test_apply_without_rescale_changes_the_mean_as_the_factors_say():
    early, late = classified_pair(seed=4, n_scale=1.5)
    sf = ts.finish_oracle(ts.fit_state_factors(early, late, "earlier", ["cm"], n_hvg=60, n_pcs=8, **FALLBACKS))
    pop, _ = reweight(early, {"cm": 0.7, "surrounding": 0.3}, 800, 0, within={"surrounding": ["Foregut"]})
    plain, _ = ts.apply_state_factors(pop, sf, HELD, [0, 11], rescale=False)
    held, _ = ts.apply_state_factors(pop, sf, HELD, [0, 11], rescale=True)
    cats = np.asarray(pop.obs["category"]).astype(str)
    m0 = ts.category_mean_cp10k(pop.layers["counts"], cats, "cm")
    m_plain = ts.category_mean_cp10k(plain.layers["counts"], cats, "cm")
    m_held = ts.category_mean_cp10k(held.layers["counts"], cats, "cm")
    assert np.abs(m_held - m0).sum() < np.abs(m_plain - m0).sum()  # the rescale is what holds the mean


def test_cells_in_no_edit_and_cells_without_a_metacell_keep_their_counts():
    pop, sf = _toy_population(seed=2)
    mask = np.zeros(pop.n_obs, dtype=bool)
    mask[:20] = True  # R1's held types, say
    pop.obs.loc[pop.obs_names[250:], "source_obs"] = "not_a_source_cell"  # a cell drawn from outside the fit's pool
    out, _ = ts.apply_state_factors(pop, sf, [], 3, rescale=False, no_edit=mask)
    before, after = sparse.csr_matrix(pop.layers["counts"]), sparse.csr_matrix(out.layers["counts"])
    unmapped = np.flatnonzero(pop.obs["source_obs"].to_numpy() == "not_a_source_cell")
    assert np.array_equal(before[:20].toarray(), after[:20].toarray())
    assert np.array_equal(before[unmapped].toarray(), after[unmapped].toarray())
    assert (before[20:200] != after[20:200]).nnz > 0  # the other cells did move


def test_apply_refuses_other_genes():
    pop, sf = _toy_population()
    bad = replace(sf, genes=list(reversed(sf.genes)))
    with pytest.raises(ValueError, match="genes differ"):
        ts.apply_state_factors(pop, bad, [], 0, rescale=False)


# ------------------------------------------------------------------------------------------------ g_mom in the harness
def _t1_inputs(rng_seed=7):
    """E9.5-like and E8.5-like labels (the config's per-file rules apply) and a target, with within-type structure."""
    ref = states_stage(rng_seed, {"V-CM": 1200, "Endocardium": 350, "Blood": 100, "Foregut": 350}, later=True)
    start = states_stage(rng_seed + 1, {"V-CM": 700, "Endothelium": 250, "Blood": 80, "Foregut": 170})
    target = states_stage(rng_seed + 2, {"V-CM": 2200, "Endocardium": 800, "Blood": 200, "Foregut": 800}, later=True,
                          target=True)
    return ref, start, target


def _t1_task():
    return PseudoTask(key="t1gm", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      twin_board="T1:val")


def test_b0_population_is_the_population_b0_builds():
    ref, start, _ = _t1_inputs()
    board = Board(key="pseudo:x", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600)
    r0, s0 = harness._shallow(ref), harness._shallow(start)
    shares = harness.t1_chamber_shares(r0, "E9.5_RNA")
    sub, _ = build_submission(r0, board, shares, 400, 1000011, start=s0, damping=0.5, held_flat=HELD,
                              ref_source="E9.5_RNA", start_source="E8.5_RNA")
    r1, s1 = harness._shallow(ref), harness._shallow(start)
    pop = ts.b0_population(r1, s1, shares, 400, 1000011, HELD, ref_source="E9.5_RNA", start_source="E8.5_RNA",
                           damping=0.5)
    assert np.array_equal(pop[:, GENES].X.toarray(), sub.X)
    assert np.array_equal(np.asarray(pop.obs["celltype"]), np.asarray(pop.obs["celltype"]))


def _gmom_ctx(ref, start, seed=1000000):
    board = Board(key="pseudo:t1gm", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600)
    task = harness.TaskView.of(_t1_task())
    return harness.BuildContext(task=task, board=board, inputs={"ref": ref, "start": start}, seed=seed,
                                derived={"t1_held_flat": tuple(HELD)})


def test_g_mom_with_empty_factors_is_b0(monkeypatch):
    """With no within-type factor to apply, g_mom's cells are exactly B0's (the base is untouched and common)."""
    ref, start, _ = _t1_inputs()
    ctx = _gmom_ctx(ref, start)
    b0 = harness._b0_t1(harness.B0_DAMPING)(ctx)
    empty = ts.StateFactors(genes=list(ref.var_names), metacells_on="later", categories={}, skipped={})
    monkeypatch.setattr(harness, "_g_mom_factors", lambda *a, **k: empty)
    gm = harness._g_mom(ctx)
    assert np.array_equal(gm.X, b0.X)


def test_g_mom_edits_real_cells_within_a_type_and_leaves_r1s_held_types_alone(monkeypatch):
    ref, start, _ = _t1_inputs()
    harness._G_MOM_PREP.clear()
    ctx = _gmom_ctx(ref, start)
    b0 = harness._b0_t1(harness.B0_DAMPING)(ctx)
    gm = harness._g_mom(ctx)
    assert gm.X.shape == b0.X.shape and not np.array_equal(gm.X, b0.X)
    # the fold's factors are fitted once, on the reference's metacells, and removed of their per-category mean
    factors = harness._G_MOM_PREP["fold"]["factors"]
    assert factors.metacells_on == "later" and "cm" in factors.categories
    cm = factors.categories["cm"]
    assert set(cm.index) <= set(ref.obs_names)
    # the second run reuses the fold's factors: the fit is not repeated (the object stays the same)
    harness._g_mom(_gmom_ctx(ref, start, seed=1000001))
    assert harness._G_MOM_PREP["fold"]["factors"] is factors
    # within a chamber category the mean of the applied log factors is about zero per gene, up to zero and clip
    weighted = np.average(cm.logf, axis=0, weights=cm.sizes)
    assert np.abs(weighted).mean() < 0.25 * np.abs(cm.logf).mean()


def test_g_mom_keeps_r1_held_cells_at_their_b0_counts(monkeypatch):
    """Blood cells (no Feng group) get B0's momentum but no within-type factor: on the edit's own cells this reads as
    identical expression to B0's for those rows."""
    ref, start, _ = _t1_inputs()
    harness._G_MOM_PREP.clear()
    ctx = _gmom_ctx(ref, start)
    b0, gm = harness._b0_t1(harness.B0_DAMPING)(ctx), harness._g_mom(ctx)
    # the rows of Blood in B0's and g_mom's builds are the same rows: the same draw (b0_population's test)
    r0 = harness._shallow(ref)
    shares = harness.t1_chamber_shares(r0, "E9.5_RNA")
    pop = ts.b0_population(harness._shallow(ref), harness._shallow(start), shares, 600, ctx.seed, list(HELD),
                           ref_source="E9.5_RNA", start_source="E8.5_RNA", damping=harness.B0_DAMPING)
    labs = pop.obs["celltype"].astype(str).to_numpy()
    assert (labs == "Blood").any()
    assert np.array_equal(gm.X[labs == "Blood"], b0.X[labs == "Blood"])
    assert not np.array_equal(gm.X[labs == "V-CM"], b0.X[labs == "V-CM"])


def test_g_mom_builds_inside_the_read_guard_on_a_t1_pseudo_task(monkeypatch, tmp_path):
    ref, start, tgt = _t1_inputs()
    harness._G_MOM_PREP.clear()
    task = _t1_task()
    board = Board(key="pseudo:t1gm", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600,
                  truth_cells=200, ref_cells=300)
    calls = {"scored": [], "inputs": []}
    st = {"target": Stage(X=tgt.X.toarray()), "ref": Stage(X=ref.X.toarray())}

    def read(path, genes):
        calls["scored"].append(str(path).split(str(tmp_path) + "/")[1])
        s = st["target"] if str(path).endswith(task.target) else st["ref"]
        return ad.AnnData(X=np.asarray(s.X), var=pd.DataFrame(index=GENES))

    by_file = {"T1/E9.5_RNA.h5ad": ref, "T1/E8.5_RNA.h5ad": start}
    monkeypatch.setattr(harness, "read_stage", read)
    monkeypatch.setattr(harness, "load_stage", lambda path, **kw: by_file[str(path).split(str(tmp_path) + "/")[1]])
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    monkeypatch.setattr(harness, "derived_ok", lambda name, t: True)
    res = harness.run_task(task, [CANDIDATES["b0"], CANDIDATES["g_mom"]], seeds=(0,), data_dir=tmp_path,
                           boards={"T1:val": replace(board, key="T1:val")}, log=None)
    assert res["candidates"] == ["floor", "b0", "g_mom"] and res["rows"][0]["scores"]["floor"] == 50.0
    assert np.isfinite(res["rows"][0]["scores"]["g_mom"]) and res["params"]["g_mom"]["scope"] == "chamber"
    assert not any("results" in p for p in res["rows"][0]["reads"]["g_mom"])
    assert not any(p.startswith(("data/", "runs/")) for p in res["rows"][0]["reads"]["g_mom"])


def test_g_mom_applies_where_b0_does_and_needs_the_held_flat_list():
    tasks = {t.key: t for ts_ in harness.load_tasks().values() for t in ts_}
    gm = CANDIDATES["g_mom"]
    assert gm.applies(tasks["T1/E9.5->E10.5"]) and gm.applies(tasks["T1/E9.5->E10.5:chamber"])
    assert not gm.applies(tasks["T1/E8.5->E9.5"]) and not gm.applies(tasks["T2h/<=E9.5->E10.5"])
    assert gm.params["scope"] == CANDIDATES["b0"].params["scope"] and gm.derived == ("t1_held_flat",)


# ------------------------------------------------------------------------------------------------ O1's arms end to end
def _load_script(name):
    path = REPO_ROOT / "experiments" / "t1_scope" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"t1_scope_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


E2E_SCRIPT = """\
import importlib.util, json, sys
sys.path[:0] = [{src!r}, {tests!r}]
import anndata as ad
import numpy as np
from test_t1_scope import E85_SIZES, E95_SIZES, GENES, t1_board
from test_t1_state_transport import states_stage

spec = importlib.util.spec_from_file_location("pt", {pseudo!r})
pt = importlib.util.module_from_spec(spec); sys.modules["pt"] = pt; spec.loader.exec_module(pt)


def libraries(*stages):
    a = ad.concat(stages)
    a.obs_names = [f"c{{i}}_{{k + 1}}" for k, s in enumerate(stages) for i in range(s.n_obs)]
    return a


sizes85 = {{**E85_SIZES, "V-CM": 500, "Endothelium": 200, "JCF": 120, "Foregut": 220}}
sizes95 = {{"V-CM": 450, "Endocardium": 160, "Proepicardium": 90, "BEC": 90, "Blood": 80, "aSHF": 60, "Foregut": 180}}
files = {{"E8.5_RNA.h5ad": states_stage(60, sizes85),
         "E9.5_RNA.h5ad": libraries(states_stage(61, sizes95, later=True), states_stage(62, sizes95, later=True))}}
pt.load_stage = lambda path, lineage=True: files[__import__("pathlib").Path(path).name].copy()
pt.get_board = lambda key: t1_board(GENES, min_cells=10, max_cells=5000)
pt.git_state = lambda: {{"commit": "t", "uncommitted": [], "sha256": {{"x": "y"}}}}
pt.sha256_file = lambda p, chunk=0: "0" * 64
pt.provenance.git_state = lambda *a, **k: ("deadbeef", [])
pt.RESULTS = __import__("pathlib").Path({results!r})
run = pt.typer.main.get_command(pt.app).main
args = ["--ref", "E8.5_RNA.h5ad", "--truth-stage", "E9.5_RNA.h5ad", "--truth-library", "1",
        "--contamination-labels", "Foregut", "--seeds", "0", "--jobs", "1", "--no-smoke", "--n-cells", "300",
        "--frac", "0.5", "--dampings", "1.0", "--flip-dampings", "",
        "--only", "sibling_real,state_transport_1.0,state_transport_b0like_1.0", "--o1-arms", "--tag", "o1t"]
run(args, standalone_mode=False)
out = pt.RESULTS / "o1t"
scores = __import__("pandas").read_csv(out / "scores.csv")
arms = set(scores["candidate"])
want = {{"copy_last", "oracle", "shift_1", "sibling_real", "state_transport_1.0", "state_transport_b0like_1.0"}}
assert arms == want, arms
ran = json.loads((out / "run.json").read_text())
assert ran["o1"]["categories_fitted"] and ran["provenance"]["commit"] == "deadbeef" and ran["provenance"]["tree_clean"]
assert ran["provenance"]["inputs_sha256"] == ran["inputs"]
reuse = json.loads((out / "o1_reuse.json").read_text())["cells_reused_per_seed"]
assert set(reuse) == {{"oracle", "sibling_real", "state_transport_b0like_1.0"}} and "0" in reuse["sibling_real"]
factors = json.loads((out / "o1_factors.json").read_text())
assert "cm" in factors["categories"] and factors["categories"]["cm"]["n_cells_a"] == 500
# a second run with the same tag is refused before anything is read (the row's no-overwrite rule)
try:
    run(args, standalone_mode=False)
except pt.typer.BadParameter as e:
    assert "already holds" in str(e)
else:
    raise SystemExit("the second run was not refused")
print("E2E OK")
"""


def test_o1_arms_run_end_to_end_through_the_script(tmp_path):
    """The whole script on synthetic stages with --o1-arms, in a fresh process (a real run forks its fit child and its
    scoring workers; this process has already used OpenMP in other tests, so it must not be the parent)."""
    code = E2E_SCRIPT.format(src=str(REPO_ROOT / "src"), tests=str(REPO_ROOT / "tests"),
                             pseudo=str(REPO_ROOT / "experiments" / "t1_scope" / "pseudo_task.py"),
                             results=str(tmp_path / "results"))
    (tmp_path / "e2e.py").write_text(code)
    env = {**os.environ, "OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
           "NUMPY_MADVISE_HUGEPAGE": "0"}
    done = subprocess.run([sys.executable, "-I", "-B", str(tmp_path / "e2e.py")], capture_output=True, text=True,
                          timeout=900, env=env, cwd=tmp_path)
    assert done.returncode == 0, done.stderr[-2000:]
    assert "E2E OK" in done.stdout


# ------------------------------------------------------------------------------------------------ the readout
@pytest.fixture
def readout(monkeypatch, tmp_path):
    mod = _load_script("o1_readout")
    monkeypatch.setattr(mod, "RESULTS", tmp_path)
    monkeypatch.setattr(mod.provenance, "git_state", lambda *a, **k: ("abc123", []))
    return mod


ARMS = ("copy_last", "oracle", "shift_1", "sibling_real", "state_transport_1.0", "state_transport_b0like_1.0")


def _scores(gain_per_seed, prize=3.0, de_gap=0.0):
    """Scores of all arms over ten seeds: state_transport_1.0 ahead of shift_1 by gain_per_seed[s] points of
    mmd_u + variogram, sibling_real by `prize`, and de_score / de_direction gaps of `de_gap` skill."""
    rows = []
    for s in range(10):
        for arm in ARMS:
            base = {"copy_last": 0.60, "oracle": 0.70, "shift_1": 0.72}.get(arm, 0.72)
            mmd, var = base, base
            de_s, de_d = 0.80, 0.75
            if arm == "state_transport_1.0":
                mmd += gain_per_seed[s] / 100.0 / 0.5  # the gain in points split between mmd_u and the variogram at
                var += gain_per_seed[s] / 100.0 / 0.5  # weights 0.30 and 0.20: 100 x (0.3 + 0.2) x d = gain
                de_s, de_d = 0.80 + de_gap, 0.75 + de_gap
            if arm == "sibling_real":
                mmd += prize / 100.0 / 0.5
                var += prize / 100.0 / 0.5
            if arm == "state_transport_b0like_1.0":
                mmd, var = base - 0.05, base - 0.04
            rows.append({"seed": s, "candidate": arm, "score": 60.0, "skill:de_score": de_s, "skill:de_direction": de_d,
                         "skill:mmd_u": mmd, "skill:variogram": var})
    return pd.DataFrame(rows)


def _write_o1_run(mod, tag, lib, gains, prize=3.0, de_gap=0.0, commit="abc123", clean=True, n_cells=2500):
    d = mod.RESULTS / tag
    d.mkdir(parents=True)
    _scores(gains, prize, de_gap).to_csv(d / "scores.csv", index=False)
    (d / "o1_reuse.json").write_text(json.dumps({"cells_reused_per_seed": {
        "sibling_real": {str(s): {"cm": 0, "endocardium": 2} for s in range(10)}}}))
    (d / "run.json").write_text(json.dumps({"n_cells": n_cells, "seeds": list(range(10)), "truth": {"library": lib},
                                            "o1": {"arms": ["x"]},
                                            "provenance": {"commit": commit, "tree_clean": clean}}))


def _repro_pass(mod, tmp_path):
    (tmp_path / "o1_repro").mkdir()
    (tmp_path / "o1_repro" / "o1_repro.json").write_text(json.dumps({"verdict": "pass"}))


GOOD = [2.0, 2.5, 1.8, 2.2, 1.9, 2.1, 2.0, 1.7, 2.3, 2.1]  # mean 2.06, ten seeds ahead


def test_decide_clears_the_bar_only_with_the_mean_and_the_seeds(readout, tmp_path):
    _repro_pass(readout, tmp_path)
    _write_o1_run(readout, "o1a", "1", GOOD)
    _write_o1_run(readout, "o1b", "2", GOOD)
    readout.main(["decide", "--lib1", "o1a", "--lib2", "o1b"])
    res = json.loads((tmp_path / "o1_decision" / "o1_decision.json").read_text())
    dec = res["decision"]
    assert dec["verdict"] == "continue" and dec["bar_cleared"] and dec["seeds_ahead_of_10"] == 10
    assert dec["mean_transport_minus_shift_1"] == pytest.approx(np.mean(GOOD))
    assert dec["prize_P_mean_of_truths"] == pytest.approx(3.0) and dec["reach_R_ratio_of_means"] == pytest.approx(
        np.mean(GOOD) / 3.0)
    assert res["truths"]["1"]["seeds_ahead"] == 10
    t = (tmp_path / "o1_decision" / "o1_decision.txt").read_text()
    assert "verdict continue" in t and "Commit abc123, tree clean True" in t
    # the outputs are stamped and never overwritten
    assert res["stamp"]["commit"] == "abc123" and res["stamp"]["tree_clean"] is True and res["stamp"]["inputs"]
    with pytest.raises(SystemExit, match="already holds"):
        readout.main(["decide", "--lib1", "o1a", "--lib2", "o1b"])


@pytest.mark.parametrize("gains, verdict", [
    ([1.0] * 10, "stop"),  # the mean gain is under the bar
    ([3.0] * 7 + [-0.5] * 3, "stop"),  # the mean clears +1.5 (1.95) but only 7 of 10 seeds are ahead
])
def test_decide_stops_on_a_small_gain_or_too_few_seeds(readout, tmp_path, gains, verdict):
    _repro_pass(readout, tmp_path)
    _write_o1_run(readout, "o1a", "1", gains)
    _write_o1_run(readout, "o1b", "2", gains)
    readout.main(["decide", "--lib1", "o1a", "--lib2", "o1b"])
    dec = json.loads((tmp_path / "o1_decision" / "o1_decision.json").read_text())["decision"]
    assert dec["verdict"] == verdict and dec["bar_cleared"] is False


def test_decide_averages_seeds_across_the_two_truths_for_the_seeds_count(readout, tmp_path):
    """The seed count is on the mean of the two truths' differences of each seed; one bad truth seed can be carried."""
    _repro_pass(readout, tmp_path)
    a = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, -1.0, -2.0]  # 8 ahead here
    b = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]  # all ahead here, so the means are ahead on 10 seeds
    _write_o1_run(readout, "o1a", "1", a)
    _write_o1_run(readout, "o1b", "2", b)
    readout.main(["decide", "--lib1", "o1a", "--lib2", "o1b"])
    dec = json.loads((tmp_path / "o1_decision" / "o1_decision.json").read_text())["decision"]
    assert dec["verdict"] == "continue" and dec["seeds_ahead_of_10"] == 10


def test_decide_is_void_when_the_pseudobulk_was_not_held(readout, tmp_path):
    _repro_pass(readout, tmp_path)
    _write_o1_run(readout, "o1a", "1", GOOD, de_gap=0.03)
    _write_o1_run(readout, "o1b", "2", GOOD)
    readout.main(["decide", "--lib1", "o1a", "--lib2", "o1b"])
    dec = json.loads((tmp_path / "o1_decision" / "o1_decision.json").read_text())["decision"]
    assert dec["verdict"] == "void" and not dec["de_within_tolerance"] and dec["bar_cleared"]


def test_decide_refuses_without_a_passed_repro_or_from_a_dirty_or_mixed_run(readout, tmp_path):
    _write_o1_run(readout, "o1a", "1", GOOD)
    _write_o1_run(readout, "o1b", "2", GOOD)
    with pytest.raises(SystemExit, match="reproduction gate"):
        readout.main(["decide", "--lib1", "o1a", "--lib2", "o1b"])
    _repro_pass(readout, tmp_path)
    _write_o1_run(readout, "dirty", "2", GOOD, clean=False)
    with pytest.raises(SystemExit, match="clean tree"):
        readout.main(["decide", "--lib1", "o1a", "--lib2", "dirty"])
    _write_o1_run(readout, "other", "2", GOOD, commit="zzz")
    with pytest.raises(SystemExit, match="different commits"):
        readout.main(["decide", "--lib1", "o1a", "--lib2", "other"])
    _write_o1_run(readout, "big", "2", GOOD, n_cells=5000)
    with pytest.raises(SystemExit, match="n_cells"):
        readout.main(["decide", "--lib1", "o1a", "--lib2", "big"])


def _repro_scores(delta_score=0.0, delta_skill=0.0):
    rows = []
    for s in range(10):
        for arm in ("copy_last", "oracle", "shift_1"):
            rows.append({"seed": s, "candidate": arm, "score": 60.0 + delta_score, "skill:de_score": 0.8 + delta_skill,
                         "skill:de_direction": 0.7, "skill:mmd_u": 0.65, "skill:variogram": 0.9})
    return pd.DataFrame(rows)


def _write_repro(mod, tag, lib, committed, delta_score=0.0, delta_skill=0.0):
    for t, ds, dsk in ((tag, delta_score, delta_skill), (committed, 0.0, 0.0)):
        d = mod.RESULTS / t
        d.mkdir(parents=True, exist_ok=True)
        _repro_scores(ds, dsk).to_csv(d / "scores.csv", index=False)
        (d / "run.json").write_text(json.dumps({"n_cells": 5000, "seeds": list(range(10)), "frac": 0.1,
                                                "truth": {"library": lib},
                                                "provenance": {"commit": "abc", "tree_clean": True}}))


def test_repro_passes_within_its_tolerance_and_refuses_beyond_it(readout, tmp_path):
    _write_repro(readout, "new1", "1", "pseudo_E8.5_E9.5_lib1", delta_score=0.01, delta_skill=0.001)
    _write_repro(readout, "new2", "2", "pseudo_E8.5_E9.5_lib2")
    readout.main(["repro", "--lib1", "new1", "--lib2", "new2"])
    rec = json.loads((tmp_path / "o1_repro" / "o1_repro.json").read_text())
    assert rec["verdict"] == "pass" and rec["max_abs_differences"]["1"]["oracle"]["score"] == pytest.approx(0.01)
    with pytest.raises(SystemExit, match="already holds"):
        readout.main(["repro", "--lib1", "new1", "--lib2", "new2"])
    _write_repro(readout, "off1", "1", "pseudo_E8.5_E9.5_lib1", delta_score=0.05)
    with pytest.raises(SystemExit, match="reproduction gate failed"):
        readout.main(["repro", "--lib1", "off1", "--lib2", "new2", "--out", "repro2"])
    assert json.loads((tmp_path / "repro2" / "o1_repro.json").read_text())["verdict"] == "fail"
    with pytest.raises(SystemExit):
        _write_repro(readout, "off2", "2", "pseudo_E8.5_E9.5_lib2", delta_skill=0.01)
        readout.main(["repro", "--lib1", "new1", "--lib2", "off2", "--out", "repro3"])


def test_repro_refuses_a_run_at_the_wrong_cell_count(readout, tmp_path):
    _write_repro(readout, "new1", "1", "pseudo_E8.5_E9.5_lib1")
    _write_repro(readout, "new2", "2", "pseudo_E8.5_E9.5_lib2")
    run = json.loads((tmp_path / "new1" / "run.json").read_text())
    run["n_cells"] = 2500
    (tmp_path / "new1" / "run.json").write_text(json.dumps(run))
    with pytest.raises(SystemExit, match="n_cells"):
        readout.main(["repro", "--lib1", "new1", "--lib2", "new2"])


# ------------------------------------------------------------------------------------------------ T11's group margins
def _t11_result(key, margins, cands=("floor", "b0", "g_mom")):
    """A minimal `vec.harness.run_task` JSON: g_mom against b0 with the given mean margins (board points) per group."""
    weights = {"de_recovery": 0.25, "direction": 0.25, "distribution": 0.30, "covariance": 0.20}
    rows = []
    for s in range(10):
        base = {"de_recovery": 0.6, "direction": 0.6, "distribution": 0.7, "covariance": 0.7}
        change = {g: base[g] + margins.get(g, 0.0) / (100.0 * w) for g, w in weights.items()}
        rows.append({"seed": s, "groups": {"b0": base, "g_mom": change},
                     "scores": {"b0": 60.0, "g_mom": 60.0 + sum(margins.values())}})
    return {"task": {"key": key}, "candidates": list(cands), "rows": rows}


def _write_t11(tmp_path, a, b, folder="t11", **kw):
    run = tmp_path / folder
    run.mkdir()
    (run / "T1_E9.5_E10.5.json").write_text(json.dumps(_t11_result("T1/E9.5->E10.5", a, **kw)))
    (run / "T1_E9.5_E10.5_chamber.json").write_text(json.dumps(_t11_result("T1/E9.5->E10.5:chamber", b, **kw)))
    return run


def test_gmom_readout_passes_only_when_the_gain_sits_on_distribution_and_covariance(readout, tmp_path):
    good = {"distribution": 0.8, "covariance": 0.5, "de_recovery": 0.1, "direction": -0.1}
    run = _write_t11(tmp_path, good, good)
    readout.main(["gmom", "--runs", str(run), "--out", "groups"])
    rec = json.loads((tmp_path / "groups" / "gmom_groups.json").read_text())
    cut = rec["cuts"]["T1/E9.5->E10.5"]
    assert rec["rule2_gain_on_predicted_groups"] is True
    assert cut["margins"]["distribution"] == pytest.approx(0.8)
    assert cut["sum_of_margins"] == pytest.approx(1.3) and cut["score_margin"] == pytest.approx(1.3)
    assert cut["share_on_distribution_and_covariance"] == pytest.approx(1.3 / 1.3)
    assert rec["stamp"]["commit"] == "abc123" and len(rec["stamp"]["inputs"]) == 2
    with pytest.raises(SystemExit, match="already holds"):
        readout.main(["gmom", "--runs", str(run), "--out", "groups"])
    off = {"distribution": 0.2, "covariance": 0.1, "de_recovery": 0.5, "direction": 0.3}  # the gain sits elsewhere
    run2 = _write_t11(tmp_path, good, off, folder="t11b")
    readout.main(["gmom", "--runs", str(run2), "--out", "groups2"])
    rec2 = json.loads((tmp_path / "groups2" / "gmom_groups.json").read_text())
    assert rec2["rule2_gain_on_predicted_groups"] is False
    assert rec2["cuts"]["T1/E9.5->E10.5:chamber"]["share_on_distribution_and_covariance"] == pytest.approx(0.3 / 1.1)
    loss = {"distribution": -0.2, "covariance": -0.1}
    run3 = _write_t11(tmp_path, good, loss, folder="t11c")
    readout.main(["gmom", "--runs", str(run3), "--out", "groups3"])
    rec3 = json.loads((tmp_path / "groups3" / "gmom_groups.json").read_text())
    assert rec3["rule2_gain_on_predicted_groups"] is False
    assert rec3["cuts"]["T1/E9.5->E10.5:chamber"]["share_on_distribution_and_covariance"] is None


def test_gmom_readout_refuses_a_missing_cut_or_a_run_without_g_mom(readout, tmp_path):
    run = _write_t11(tmp_path, {"distribution": 1.0}, {"distribution": 1.0})
    (run / "T1_E9.5_E10.5.json").unlink()
    with pytest.raises(SystemExit, match="no result for T1/E9.5->E10.5"):
        readout.main(["gmom", "--runs", str(run), "--out", "g1"])
    run2 = _write_t11(tmp_path, {"distribution": 1.0}, {"distribution": 1.0}, folder="t11d", cands=("floor", "b0"))
    with pytest.raises(SystemExit, match="with b0 and g_mom"):
        readout.main(["gmom", "--runs", str(run2), "--out", "g2"])

