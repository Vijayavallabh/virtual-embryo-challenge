"""Heart tries H7 `b0_gmom` and H7c `b0_gmom_null` (experiments/harness/TRIES.md), on synthetic data only: T1 try
T11's within-type momentum (vec.t1_state_transport) on B0's T2 cells, with the harmonised lineage as the category."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec import t1_state_transport as ts
from vec.boards import Board
from vec.harness import BUILD_SEED_OFFSET, CANDIDATES, PseudoTask, load_tasks

GENES = [f"g{i}" for i in range(120)]
LINEAGES = np.array(["CM", "Endothelial", "Untracked", "Tiny"])


def _stage(n: int, seed: int, shift: float, tiny: int) -> ad.AnnData:
    """A T2 stage: CM and Endothelial lineages with a within-lineage state axis (a module of genes that a cell's state
    score moves), Untracked cells, and `tiny` cells of a lineage too small for T11's minimum at this stage. `shift`
    moves the CM states' module along the axis between stages, so the within-lineage change differs by state."""
    rng = np.random.default_rng(seed)
    types = rng.choice(3, size=n, p=[0.5, 0.35, 0.15])
    types[:tiny] = 3
    pi = np.random.default_rng(1).gamma(1.0, 1.0, len(GENES))
    pi /= pi.sum()
    a = np.zeros(len(GENES))
    a[:30] = np.random.default_rng(4).choice([-1.0, 1.0], 30) * 1.5
    s = rng.normal(size=n)
    lam = pi * np.exp(np.outer(s, a) + np.outer(shift * s * (types == 0), a) + 0.3 * (types == 1)[:, None])
    depth = np.maximum(np.round(np.exp(rng.normal(np.log(80.0), 0.4, n))), 5)[:, None]
    C = rng.poisson(lam / lam.sum(1, keepdims=True) * depth).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    X = np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32)
    centres = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [5.0, 5.0, 0.0]])
    out = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame(
        {"celltype": [f"type{t}" for t in types], "merged": LINEAGES[types]},
        index=[f"c{seed}_{i}" for i in range(n)]), var=pd.DataFrame(index=GENES))
    out.layers["counts"] = sparse.csr_matrix(C)
    out.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * 1.5).astype(np.float32)
    return out


def _heart_task(**over) -> PseudoTask:
    kw = {"key": "syn_h7", "test_board": "T2:heart:test_extrap", "task": "T2", "kind": "extrap",
          "target": "T2/E10.5.h5ad", "reference": "T2/E9.5.h5ad", "setting": "heart",
          "inputs": {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"}, "validation_board": "T2:heart:val_extrap"}
    return PseudoTask(**{**kw, **over})


BOARD = Board(key="pseudo:h7", task="T2", genes=list(GENES), needs_coords=True, setting="heart", min_cells=100,
              max_cells=700)


@pytest.fixture(scope="module")
def inputs():
    return {"ref": _stage(1600, seed=2, shift=0.0, tiny=60), "start": _stage(1100, seed=3, shift=0.8, tiny=10)}


def _ctx(inputs, seed: int = BUILD_SEED_OFFSET + 1) -> harness.BuildContext:
    return harness.BuildContext(task=harness.TaskView.of(_heart_task()), board=BOARD, inputs=inputs, seed=seed)


def _merged_of_draw(ctx) -> np.ndarray:
    pool, rows, _, _ = harness.draw_copy(ctx)
    return np.asarray(pool.obs["merged"]).astype(str)[rows]


def test_the_candidates_apply_where_b0s_t2_extrapolation_arms_do():
    tasks = {t.key: t for ts_ in load_tasks().values() for t in ts_}
    for name in ("b0_gmom", "b0_gmom_null"):
        c = CANDIDATES[name]
        assert c.applies(tasks["T2h/<=E9.5->E10.5"]) and c.applies(tasks["T2h/<=E8.75->E10.5"])
        assert not c.applies(tasks["T1/E9.5->E10.5"]) and not c.applies(tasks["T3/Gata4->Mab21l2"])
        assert not c.applies(_heart_task(inputs={"ref": "T2/E9.5.h5ad"}))
        assert all(c.applies(t) == CANDIDATES["b0_cfm"].applies(t) for t in tasks.values())
    assert CANDIDATES["b0_gmom_null"].params["split_seed"] == harness.G_MOM_NULL_SEED


def test_the_fit_uses_t11s_rules_on_the_lineages(inputs):
    """Lineages at both stages, Untracked never; a lineage below T11's minimum at a stage is skipped; metacells on the
    reference's cells by T11's rule; the PCA on the shared panel genes (capped components)."""
    finished, raw = harness.t2_gmom_pair_factors(inputs["start"], inputs["ref"])
    assert harness.t2_lineages(inputs["start"], inputs["ref"]) == ["CM", "Endothelial", "Tiny"]
    assert set(raw.categories) == {"CM", "Endothelial"} and "Untracked" not in raw.categories
    assert "Tiny" in raw.skipped and str(ts.MIN_CELLS) in raw.skipped["Tiny"]
    assert raw.genes == GENES and raw.metacells_on == "later"
    assert raw.record["n_pcs"] == ts.N_PCS and raw.record["n_hvg"] <= len(GENES)
    ref_names = np.asarray(inputs["ref"].obs_names)
    merged = np.asarray(inputs["ref"].obs["merged"]).astype(str)
    for c, cf in raw.categories.items():
        n = int((merged == c).sum())
        assert set(cf.index) == set(ref_names[merged == c])
        assert cf.n_metacells == min(n, int(np.clip(n // ts.METACELL_DIVISOR, ts.METACELL_MIN, ts.METACELL_MAX)))


def test_finished_factors_are_mean_free_per_lineage_before_the_zero_and_clip(inputs):
    finished, raw = harness.t2_gmom_pair_factors(inputs["start"], inputs["ref"])
    for c, cf in raw.categories.items():
        removed = (cf.logf - np.average(cf.logf, axis=0, weights=cf.sizes)) * harness.B0_DAMPING
        np.testing.assert_allclose(np.average(removed, axis=0, weights=cf.sizes), 0.0, atol=1e-12)
        np.testing.assert_array_equal(finished.categories[c].logf, ts._zero_and_clip(removed))
        out = finished.categories[c].logf
        assert np.abs(out).max() <= ts.CLIP_AT and not np.any((out != 0) & (np.abs(out) < ts.ZERO_BELOW))
    assert any(np.any(cf.logf != 0) for cf in finished.categories.values())  # the fixture moves something


def test_with_every_factor_one_b0_gmom_is_b0_to_the_last_bit(inputs):
    ctx = _ctx(inputs)
    b0 = CANDIDATES["b0"].build(ctx)

    def ones(sf):
        return ts.StateFactors(genes=sf.genes, metacells_on=sf.metacells_on, skipped=sf.skipped, record=sf.record,
                               categories={c: ts.CategoryFactors(c, cf.index, cf.sizes, np.zeros_like(cf.logf))
                                           for c, cf in sf.categories.items()})

    flat = harness._b0_gmom_arm(ctx, ones)
    np.testing.assert_array_equal(flat.X, b0.X)
    np.testing.assert_array_equal(flat.coords, b0.coords)


def test_b0_gmom_keeps_b0s_rows_and_coordinates_and_edits_only_fitted_lineages(inputs):
    ctx = _ctx(inputs)
    b0, gm = CANDIDATES["b0"].build(ctx), CANDIDATES["b0_gmom"].build(ctx)
    assert gm.n == b0.n == BOARD.max_cells
    np.testing.assert_array_equal(gm.coords, b0.coords)
    merged = _merged_of_draw(ctx)
    for kept in ("Untracked", "Tiny"):
        m = merged == kept
        assert m.any() and np.array_equal(gm.X[m], b0.X[m])
    cm = merged == "CM"
    assert not np.array_equal(gm.X[cm], b0.X[cm])
    # the edit's record: exactly the drawn cells of the fitted lineages are edited
    _, _, _, rec = harness.b0_gmom_population(ctx, harness._b0_gmom_factors(ctx))
    assert rec["n_edited"] == int(np.isin(merged, ["CM", "Endothelial"]).sum())
    assert rec["n_held_genes"] == 0


def test_the_fit_is_made_once_per_fold_and_other_seeds_draw_other_cells(inputs):
    ctx = _ctx(inputs)
    CANDIDATES["b0_gmom"].build(ctx)
    f = harness._B0_GMOM_PREP["fold"]["factors"]
    other = CANDIDATES["b0_gmom"].build(_ctx(inputs, seed=BUILD_SEED_OFFSET + 2))
    assert harness._B0_GMOM_PREP["fold"]["factors"] is f
    assert not np.array_equal(other.coords, CANDIDATES["b0"].build(ctx).coords)


def test_the_no_development_arm_has_b0_gmoms_draws_lineages_and_size(inputs):
    """H7c (T23's construction): the same base and draws, factors only for H7's lineages and cells of every such
    lineage covered by one of the two halves' metacells, at H7's root-mean-square size per lineage."""
    ctx = _ctx(inputs)
    b0, gm, null = (CANDIDATES[n].build(ctx) for n in ("b0", "b0_gmom", "b0_gmom_null"))
    np.testing.assert_array_equal(null.coords, b0.coords)
    assert not np.array_equal(null.X, gm.X)
    merged = _merged_of_draw(ctx)
    for kept in ("Untracked", "Tiny"):
        m = merged == kept
        assert np.array_equal(null.X[m], b0.X[m])
    finished = harness._B0_GMOM_PREP["fold"]["factors"]
    factors = harness._B0_GMOM_NULL_PREP["fold"]["factors"]
    assert harness._B0_GMOM_NULL_PREP["fold"]["ref"] is inputs["ref"]
    assert set(factors.categories) == set(finished.categories)
    free = np.ones(len(GENES), dtype=bool)
    ref_names = np.asarray(inputs["ref"].obs_names)
    ref_merged = np.asarray(inputs["ref"].obs["merged"]).astype(str)
    for c, cf in factors.categories.items():
        assert set(cf.index) == set(ref_names[ref_merged == c])  # every reference cell of the lineage has a factor
        goal = finished.categories[c]
        assert ts._factor_rms(cf.logf, cf.sizes, free) == pytest.approx(ts._factor_rms(goal.logf, goal.sizes, free),
                                                                        rel=0.1)
        assert cf.record["no_development"]["rms_target"] > 0
    _, _, _, rec_g = harness.b0_gmom_population(ctx, finished)
    _, _, _, rec_n = harness.b0_gmom_population(ctx, factors)
    assert rec_n["n_edited"] == rec_g["n_edited"]
    CANDIDATES["b0_gmom_null"].build(_ctx(inputs, seed=BUILD_SEED_OFFSET + 2))
    assert harness._B0_GMOM_NULL_PREP["fold"]["factors"] is factors


def test_interpolation_is_refused(inputs):
    from dataclasses import replace

    ctx = _ctx(inputs)
    with pytest.raises(ValueError, match="extrapolation"):
        harness.b0_gmom_population(replace(ctx, task=replace(ctx.task, kind="interp")), None)
