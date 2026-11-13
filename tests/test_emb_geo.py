"""The embryo geometry batch (branch campaign/emb-geo): fgw_bluenoise (an even-spacing relaxation of fgw's positions)
and fgw_depthdraw (which real cells of each bracket fgw draws, shaped by their depth towards the depth the mixture
should have), each fgw with one idea added. A partial or unbalanced coupling (the third line of the batch) is not built:
no record gives the share of mass it would leave unmatched. Synthetic data only, built here."""

from __future__ import annotations

import hashlib
from dataclasses import replace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import GENES
from scipy import sparse as sp
from scipy.spatial import cKDTree

from vec import harness
from vec import t2_build as tb
from vec.boards import Board
from vec.harness import CANDIDATES, PseudoTask, load_tasks

INTERP = dict(task="T2", setting="embryo", kind="interp", target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
              twin_board="T2:embryo:val_interp")


def _task(**over) -> PseudoTask:
    base = PseudoTask(key="x", test_board="T2:embryo:test_interp", task="T2", setting="embryo", kind="interp",
                      target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
                      inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"}, twin_board="T2:embryo:val_interp")
    return replace(base, **over)


def _counted(n, rng, base, depth, scale, flat_depth=False):
    """A counted stage: three types at three places, panel counts Poisson, the library the panel's own total (or, with
    flat_depth, every cell's library the same: counts drawn as a multinomial of one fixed total)."""
    types = rng.choice(3, size=n, p=[0.4, 0.3, 0.3])
    if flat_depth:
        p = base[types] / base[types].sum(1, keepdims=True)
        C = np.vstack([rng.multinomial(int(depth * 40), q) for q in p])
    else:
        C = rng.poisson(base[types] * depth * rng.lognormal(0.0, 0.5, size=(n, 1)))
        C[C.sum(1) == 0, 0] = 1
    X = np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32)
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float) * scale
    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"celltype": [f"type{t}" for t in types],
                                                         "group": [f"g{t}" for t in types]},
                                                        index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sp.csr_matrix(C.astype(np.int32))
    a.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * 1.5 * scale).astype(np.float32)
    return a


def _ctx(seed=7, flat_depth=False, n_before=900, n_after=1300, max_cells=700):
    rng = np.random.default_rng(29)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    pools = {"before": _counted(n_before, rng, base, 1.0, 1.0, flat_depth),
             "after": _counted(n_after, rng, base, 4.0, 1.5, flat_depth)}
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=max_cells)
    task = _task(key="syn_emb_geo")
    return harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs=pools, seed=seed), pools, board


def _digest(stage) -> str:
    """X exactly; coordinates to 1e-4 (a POT or BLAS build may move the last bits of a coupling's moves)."""
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(np.asarray(stage.X, np.float32)).tobytes())
    h.update(np.ascontiguousarray(np.round(np.asarray(stage.coords, np.float64), 4)).tobytes())
    return h.hexdigest()


# taken before this batch's edit (the parent commit's code), on _ctx(): the existing candidates stay as they were
PINNED = {
    "fgw": "24db4a5c7d50ee1d51652145dceba2f1f2ca65b12eab356153efd15aa8e12d48",
    "fgw_aligned": "2bc5769d0705ef756537c4ad052fd7fadf13aeafb8563a45ed257bd1440757df",
    "fgw_herd": "a6dd4d2377dd22419d92ef455b24fb3578d68e94d854dd5b807d0430c9775fd0",
    "fgw_midpoint": "927947b777bbec336e9cef091402c284db6cd5f79f4347878516758f27f5ec7d",
    "fgw_shape0": "d7df9b173a5b2bbea1bbc506cf0cb6b0a593f1fb110940994d6b6efd6e2f6aa9",
}


@pytest.fixture(scope="module")
def case():
    return _ctx()


@pytest.mark.parametrize("name", sorted(PINNED))
def test_existing_embryo_candidates_are_unchanged(case, name):
    ctx, _, _ = case
    assert _digest(CANDIDATES[name].build(ctx)) == PINNED[name]


# ------------------------------------------------------------------------------------------------ registration
def test_the_batch_registers_two_arms_beside_fgw_aligned_with_its_requirements():
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    for name in ("fgw_bluenoise", "fgw_depthdraw"):
        c = CANDIDATES[name]
        assert c.tasks == CANDIDATES["fgw_aligned"].tasks == frozenset({"T2"}) and not c.tuned_on
        for key in ("T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25"):
            assert c.applies(tasks[key]) == CANDIDATES["fgw_aligned"].applies(tasks[key]) is True
        for key in ("T2h/<=E9.5->E10.5", "T1/E9.5->E10.5"):
            assert not c.applies(tasks[key])
    assert CANDIDATES["fgw_bluenoise"].params["steps"] == tb.RELAX_STEPS == 1
    assert CANDIDATES["fgw_depthdraw"].params["rule"] == tb.DEPTH_DRAW_RULE
    assert not any("partial" in n or "unbalanced" in n for n in CANDIDATES)  # the third line is not built


def test_both_ideas_at_zero_strength_are_fgw_bit_for_bit(case):
    ctx, _, _ = case
    for kw in ({"relax_steps": 0}, {"depth_draw": False}, {"relax_steps": 0, "depth_draw": False}):
        assert _digest(harness._fgw(ctx, **kw)) == PINNED["fgw"]


# ------------------------------------------------------------------------------------------------ fgw_bluenoise
def _nn_cv(X) -> float:
    d = cKDTree(np.asarray(X, float)).query(np.asarray(X, float), k=2)[0][:, 1]
    return float(d.std() / d.mean())


def test_lloyd_relax_with_no_step_or_nothing_to_do_returns_its_input_exactly():
    rng = np.random.default_rng(3)
    C = rng.normal(size=(200, 3))
    P = np.vstack([C, rng.normal(size=(500, 3))])
    out, rec = tb.lloyd_relax(C, P, np.ones(len(P)), 0)
    assert np.array_equal(out, C) and rec["steps"] == 0
    out, rec = tb.lloyd_relax(C, C, np.ones(len(C)), 3)  # each point is its own pool point's centroid: nothing to do
    assert np.array_equal(out, C) and rec["static"] == [0, 0, 0]


def test_lloyd_relax_evens_the_spacing_and_keeps_count_centroid_and_rms():
    rng = np.random.default_rng(4)
    P = rng.uniform(-1.0, 1.0, size=(20000, 3))
    C = P[rng.choice(len(P), 400, replace=False)]
    out, rec = tb.lloyd_relax(C, P, np.ones(len(P)), 1)
    assert out.shape == C.shape and rec["exact_duplicates"] == 0
    np.testing.assert_allclose(out.mean(0), C.mean(0), atol=1e-12)
    assert tb.rms_radius(out) == pytest.approx(tb.rms_radius(C), rel=1e-12)
    assert _nn_cv(out) < 0.8 * _nn_cv(C)  # nearest-neighbour distances more even: the blue-noise property


def test_lloyd_relax_weights_pull_towards_the_heavier_pool_points():
    C = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    P = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [9.0, 0.0, 0.0], [11.0, 0.0, 0.0]])
    out, _ = tb.lloyd_relax(C, P, np.array([1.0, 3.0, 1.0, 1.0]), 1)
    # the steps give (0.5, 0, 0) and (10, 0, 0); the cloud then goes back to the input's centroid and RMS radius
    np.testing.assert_allclose(out[1] - out[0], [10.0, 0.0, 0.0])
    np.testing.assert_allclose(out.mean(0), C.mean(0))


def test_lloyd_relax_refuses_bad_inputs():
    C, P = np.zeros((3, 3)) + np.arange(3)[:, None], np.ones((5, 3))
    for w in (np.ones(4), -np.ones(5), np.zeros(5), np.full(5, np.nan)):
        with pytest.raises(ValueError, match="pool_weights"):
            tb.lloyd_relax(C, P, w, 1)
    with pytest.raises(ValueError):
        tb.lloyd_relax(C, P, np.ones(5), -1)


def test_fgw_bluenoise_keeps_fgws_cells_count_and_size_and_moves_only_positions(case):
    ctx, _, _ = case
    fgw, v = CANDIDATES["fgw"].build(ctx), CANDIDATES["fgw_bluenoise"].build(ctx)
    assert v.n == fgw.n == ctx.board.max_cells
    np.testing.assert_array_equal(v.X, fgw.X)  # the same drawn cells in the same order, the same expression
    np.testing.assert_allclose(np.asarray(v.coords, float).mean(0), np.asarray(fgw.coords, float).mean(0), atol=1e-3)
    assert tb.rms_radius(v.coords) == pytest.approx(tb.rms_radius(fgw.coords), rel=1e-5)
    assert not np.allclose(v.coords, fgw.coords)
    assert _nn_cv(v.coords) < _nn_cv(fgw.coords)


# ------------------------------------------------------------------------------------------------ fgw_depthdraw
def test_depth_tilt_puts_the_weighted_median_at_the_target():
    rng = np.random.default_rng(5)
    lib = rng.lognormal(np.log(200.0), 0.6, size=3000)
    for target in (90.0, 150.0, 260.0):
        lam, w = tb.depth_tilt(lib, target)
        s = np.sign(np.log(lib) - np.log(target))
        assert abs(float(s @ w)) / w.sum() < 1e-9
        assert (lam < 0) == (target < np.median(lib))  # a target below the median favours the shallower cells
        assert w.max() == pytest.approx(1.0) and w.min() >= 1e-12


def test_depth_tilt_has_nothing_to_do_when_the_cells_split_evenly_or_do_not_vary():
    lam, w = tb.depth_tilt(np.full(50, 30.0), 64.0)
    assert lam == 0.0 and np.array_equal(w, np.ones(50))
    lam, w = tb.depth_tilt(np.array([10.0, 20.0, 40.0, 80.0]), 30.0)  # two below, two above
    assert lam == 0.0 and np.array_equal(w, np.ones(4))
    with pytest.raises(ValueError, match="outside"):
        tb.depth_tilt(np.array([10.0, 20.0, 40.0]), 100.0)
    for bad in (np.array([1.0, 0.0]), np.array([1.0, np.nan]), np.array([])):
        with pytest.raises(ValueError, match="library sizes"):
            tb.depth_tilt(bad, 1.0)
    with pytest.raises(ValueError, match="target"):
        tb.depth_tilt(np.array([1.0, 2.0]), 0.0)


def _interp(ctx, **kw):
    before, after = ctx.inputs["before"], ctx.inputs["after"]
    g = ctx.board.genes
    Cb, Ca = (np.asarray(a.obsm["spatial_3D"], float)[:, :3] for a in (before, after))
    Eb, Ea = (harness.dense(a[:, g].X).astype(np.float32) for a in (before, after))
    return tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.25 / 0.75, harness.count_for(ctx.board, before.n_obs), ctx.seed,
                                 target_rms=harness.target_rms(ctx), **kw)


def _libs(ctx):
    return tuple(tb.cell_depths(ctx.inputs[r], ctx.board.genes) for r in ("before", "after"))


def test_depth_draw_changes_only_which_real_cells_each_bracket_draws(case):
    ctx, pools, _ = case
    libs = _libs(ctx)
    for a, lib in zip(pools.values(), libs, strict=True):
        assert np.median(lib) == tb.median_depth(a, ctx.board.genes)
    base, arm = _interp(ctx), _interp(ctx, depth_libs=libs)
    np.testing.assert_array_equal(arm.stage, base.stage)  # each bracket's count and the cells' order by bracket
    # the same coupling and moves: the pools agree up to the size step's similarity (it centres and scales the draw)
    P, Q = base.pool_coords - base.pool_coords.mean(0), arm.pool_coords - arm.pool_coords.mean(0)
    np.testing.assert_allclose(Q, P * (np.linalg.norm(Q) / np.linalg.norm(P)), atol=1e-9)
    assert arm.record["n_before"] == base.record["n_before"] and not np.array_equal(arm.rows, base.rows)
    for res in (base, arm):  # every cell at its own moved place (the pool), up to the free jitter of repeats
        own = res.pool_coords[np.where(res.stage == 0, res.rows, len(libs[0]) + res.rows)]
        assert np.mean(np.all(np.isclose(res.coords, own, atol=1e-9), axis=1)) > 0.95
    np.testing.assert_allclose(arm.pool_weights.sum(), 1.0)
    rec = arm.record["depth_draw"]
    target = tb.depth_target(*rec["median_depth"], 0.25 / 0.75)
    assert rec["target"] == target and rec["lam"][0] > 0 > rec["lam"][1]  # shallower bracket deeper, deeper shallower
    for s, lib in enumerate(libs):
        drawn_b, drawn_a = np.median(lib[base.rows[base.stage == s]]), np.median(lib[arm.rows[arm.stage == s]])
        assert abs(np.log(drawn_a / target)) < abs(np.log(drawn_b / target))  # each bracket's draw moves to the target


def test_fgw_depthdraw_keeps_fgws_count_and_real_cells_and_edits_no_counts(case):
    ctx, pools, _ = case
    fgw, v = CANDIDATES["fgw"].build(ctx), CANDIDATES["fgw_depthdraw"].build(ctx)
    assert v.n == fgw.n == ctx.board.max_cells and not np.array_equal(v.X, fgw.X)
    real = {row.tobytes() for a in pools.values() for row in np.asarray(a[:, GENES].X.todense(), np.float32)}
    assert all(row.tobytes() in real for row in np.asarray(v.X, np.float32))  # every row a released cell, unedited
    assert tb.rms_radius(v.coords) == pytest.approx(tb.rms_radius(fgw.coords), rel=0.01)  # the size rule's


def test_fgw_depthdraw_with_nothing_to_shape_is_fgw_bit_for_bit():
    ctx, _, _ = _ctx(flat_depth=True)  # every cell of a bracket at one library size: no depth to shape the draw by
    fgw, v = CANDIDATES["fgw"].build(ctx), CANDIDATES["fgw_depthdraw"].build(ctx)
    np.testing.assert_array_equal(v.X, fgw.X)
    np.testing.assert_array_equal(v.coords, fgw.coords)
    base, arm = _interp(ctx), _interp(ctx, depth_libs=_libs(ctx))
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
        np.testing.assert_array_equal(getattr(arm, f), getattr(base, f))
    assert arm.record["depth_draw"]["lam"] == [0.0, 0.0]


def test_depth_draw_refuses_the_options_that_also_set_the_rows_and_bad_libraries(case):
    ctx, _, _ = case
    libs = _libs(ctx)
    for kw in ({"herd": True}, {"t_mix": 1.0},
               {"draw": "midpoint", "log_depth": (np.log1p(libs[0]), np.log1p(libs[1]))}):
        with pytest.raises(ValueError, match="depth_libs"):
            _interp(ctx, depth_libs=libs, **kw)
    with pytest.raises(ValueError, match="one library size per cell"):
        _interp(ctx, depth_libs=(libs[0][:-1], libs[1]))
    with pytest.raises(ValueError, match="before, after"):
        _interp(ctx, depth_libs=(libs[0],))
