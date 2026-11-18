"""Embryo tries E47a and E47b of experiments/harness/TRIES.md, part A: the off-by-default position mode
`interpolate_clouds(..., positions="floor_exact")`, the harness candidates fgw_shape0x and fgw_shape0x_series and the
builder's --positions option, on synthetic arrays and stages only. With the default every output is what it was."""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.distance import cdist
from test_t2_build import _labelled_pair

from vec import harness
from vec import t2_build as tb
from vec.boards import Board
from vec.harness import CANDIDATES
from vec.io import COORD_KEY

REPO = Path(__file__).resolve().parents[1]
FLOOR = tb.POSITIONS[1]


def _two_types(rng, n_b=400, n_a=500):
    """Two brackets of two cell types each. Expression says the type (genes 0-4 for type 0, 5-9 for type 1) in both
    brackets, with a shared shift of the later one; in the earlier bracket type 0 sits at x < 0 and type 1 at x > 0, in
    the later one the clouds are of another shape, so a moved position differs from an earlier one."""
    tb_, ta = rng.integers(0, 2, n_b), rng.integers(0, 2, n_a)

    def expr(types, shift):
        E = rng.gamma(2.0, 0.2, size=(len(types), 10))
        E[types == 0, :5] += 3.0
        E[types == 1, 5:] += 3.0
        return E + shift

    Cb = rng.normal(size=(n_b, 3)) * [1.0, 1.0, 0.6]
    Cb[:, 0] = np.where(tb_ == 0, -1.0, 1.0) * (np.abs(Cb[:, 0]) + 1.5)
    Ca = rng.normal(size=(n_a, 3)) * [2.5, 1.0, 0.3]
    return Cb, Ca, expr(tb_, 0.0), expr(ta, 0.4), tb_, ta


KW = dict(seed=7, target_rms=5.0, n_anchor=150, k=4)


def test_default_positions_are_the_entry_bit_for_bit():
    rng = np.random.default_rng(1)
    Cb, Ca, Eb, Ea, _, _ = _two_types(rng)
    for tp in (None, 0.0, 1.0):
        base = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 300, t_pos=tp, **KW)
        same = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 300, t_pos=tp, positions=tb.POSITIONS[0], **KW)
        for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
            np.testing.assert_array_equal(getattr(base, f), getattr(same, f))
        assert base.record == same.record and "positions" not in base.record and "floor_exact" not in base.record


def test_floor_exact_keeps_the_rows_and_puts_every_cell_on_a_real_earlier_position():
    rng = np.random.default_rng(2)
    Cb, Ca, Eb, Ea, _, _ = _two_types(rng)
    n = 300
    s0 = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, n, t_pos=0.0, with_coupling=True, **KW)
    fx = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, n, t_pos=0.0, positions=FLOOR, with_coupling=True, **KW)
    # the same cells in the same order as fgw_shape0 (so every expression metric is its own), other positions
    np.testing.assert_array_equal(fx.rows, s0.rows)
    np.testing.assert_array_equal(fx.stage, s0.stage)
    assert not np.allclose(fx.coords, s0.coords)
    assert fx.record["positions"] == FLOOR and fx.record["floor_exact"]["fallback"] is False
    assert fx.record["t_pos"] == 0.0 and fx.record["n_jittered"] == 0
    # rebuild the cloud from its parts: before rows at their own positions, after rows at the assigned before cells
    Zb = fx.coupling["before_frame"]
    perm = np.random.default_rng([KW["seed"], 2]).permutation(n)  # the output order: output i is input row perm[i]
    order_b = np.flatnonzero(perm < (fx.stage == 0).sum())
    rows_b = fx.rows[order_b][np.argsort(perm[order_b])]  # the before rows in their drawn order
    at = fx.coupling["floor_exact_rows"]
    Z = np.vstack([Zb[rows_b], Zb[at]])[perm]
    want = Cb.mean(0) + (Z - Z.mean(0)) * (KW["target_rms"] / tb.rms_radius(Z))
    np.testing.assert_allclose(fx.coords, want, rtol=0, atol=1e-9)
    assert tb.rms_radius(fx.coords) == pytest.approx(KW["target_rms"])
    # the after rows' positions are before cells the before rows do not use, each used once
    assert len(np.unique(at)) == len(at) and not set(at) & set(rows_b)
    # the guard's source is the before stage alone, uniform
    assert fx.pool_coords.shape == (len(Cb), 3) and np.allclose(fx.pool_weights, 1.0 / len(Cb))
    # its before rows sit exactly where fgw_shape0 puts them, up to the cloud's own centring and scale
    sb = s0.coords[s0.stage == 0] - s0.coords[s0.stage == 0].mean(0)
    xb = fx.coords[fx.stage == 0] - fx.coords[fx.stage == 0].mean(0)
    np.testing.assert_allclose(xb / tb.rms_radius(xb), sb / tb.rms_radius(sb), atol=1e-9)


def test_the_assignment_is_the_exact_minimum_over_the_drawn_candidates():
    """Brute force on a tiny draw: no other one-to-one placement of the after rows on the same candidates costs less."""
    rng = np.random.default_rng(3)
    Cb, Ca, Eb, Ea, _, _ = _two_types(rng, n_b=60, n_a=60)
    n = 12
    fx = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, n, seed=11, target_rms=1.0, n_anchor=40, k=3, t_pos=0.0,
                               positions=FLOOR, with_coupling=True)
    c = fx.coupling
    ib, ia = c["anchors_before"], c["anchors_after"]
    centre, load = tb._joint_pc_loadings(Eb[ib], Ea[ia], tb.INTERP_N_PC)
    rows_b = np.flatnonzero(np.isin(np.arange(60), fx.rows[fx.stage == 0]))
    rows_a = tb.draw_balanced(60, n - len(rows_b), np.random.default_rng([11, 1]))
    at = c["floor_exact_rows"]
    cost = cdist((Ea[rows_a] - centre) @ load, (Eb[at] - centre) @ load, "sqeuclidean")
    best = cost.trace()
    for p in itertools.permutations(range(len(at))):
        assert cost[np.arange(len(at)), p].sum() >= best - 1e-9
    assert fx.record["floor_exact"]["mean_cost"] == pytest.approx(best / len(at))
    # the candidates are the uniform draw of the unused cells on the row's stream
    unused = np.setdiff1d(np.arange(60), rows_b)
    want = np.random.default_rng([11, tb.FLOOR_EXACT_STREAM]).choice(unused, len(rows_a), replace=False)
    assert sorted(at) == sorted(want)


def test_after_cells_land_among_earlier_cells_of_their_own_expression():
    rng = np.random.default_rng(4)
    Cb, Ca, Eb, Ea, type_b, type_a = _two_types(rng, n_b=600, n_a=600)
    fx = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 400, t_pos=0.0, positions=FLOOR, with_coupling=True, **{
        **KW, "target_rms": 1.0})
    perm = np.random.default_rng([KW["seed"], 2]).permutation(400)
    n_b = int((fx.stage == 0).sum())
    rows_a = fx.rows[np.flatnonzero(perm >= n_b)][np.argsort(perm[perm >= n_b])]
    at = fx.coupling["floor_exact_rows"]
    # an after cell of type k takes a type-k cell's position, but for the few that the candidates' type shares cannot
    # match one to one (a random draw against another random draw); a random placement would match about half
    assert (type_b[at] == type_a[rows_a]).mean() > 0.85
    # so every neighbourhood mixes the two brackets' cells of one type: x's sign says the type in the earlier geometry
    x = fx.coords[:, 0] - Cb.mean(0)[0]
    t_all = np.empty(400, dtype=int)
    t_all[fx.stage == 0], t_all[fx.stage == 1] = type_b[fx.rows[fx.stage == 0]], type_a[fx.rows[fx.stage == 1]]
    assert ((x > 0) == (t_all == 1)).mean() > 0.9


def test_the_drawn_positions_are_a_uniform_sample_of_the_earlier_bracket():
    """Before rows and candidates together are a uniform subset: every earlier cell is drawn about equally often."""
    rng = np.random.default_rng(5)
    Cb, Ca, Eb, Ea, _, _ = _two_types(rng, n_b=80, n_a=90)
    hits = np.zeros(80)
    for seed in range(200):
        fx = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 40, seed=seed, target_rms=1.0, n_anchor=30, k=3, t_pos=0.0,
                                   positions=FLOOR, with_coupling=True)
        used = np.r_[fx.rows[fx.stage == 0], fx.coupling["floor_exact_rows"]]
        assert len(np.unique(used)) == 40
        hits[used] += 1
    assert hits.min() > 0.35 * 200 * 0.5 and hits.max() < 1.65 * 200 * 0.5  # each cell about half the draws


def test_too_few_unused_cells_use_all_of_them_then_repeat_only_what_the_count_forces():
    """The fallback fires exactly when the count exceeds the earlier bracket's cells. Every unused cell is then a
    candidate and only the remainder repeats: with 100 earlier cells and a count of 130, 65 before rows and 65
    candidates, 35 unused cells, so 30 candidates repeat a position and 100 distinct positions are taken."""
    rng = np.random.default_rng(6)
    Cb, Ca, Eb, Ea, _, _ = _two_types(rng, n_b=100, n_a=300)
    for n, want in ((100, False), (101, True), (130, True), (260, True)):
        fx = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, n, seed=2, target_rms=1.0, n_anchor=60, k=3, t_pos=0.0,
                                   positions=FLOOR, with_coupling=True)
        assert fx.record["floor_exact"]["fallback"] is want, n
        assert np.unique(fx.coords.astype(np.float32), axis=0).shape[0] == n  # repeats are jittered apart
        used = np.r_[fx.rows[fx.stage == 0], fx.coupling["floor_exact_rows"]]
        assert len(np.unique(used)) == min(n, 100)  # every earlier cell is used before any repeats
        assert (fx.record["n_jittered"] > 0) is want
    assert len(fx.coupling["floor_exact_rows"]) == 130
    # each earlier cell's expected number of copies stays equal in the fallback (n = 150 over 100 earlier cells:
    # 1.5 copies per seed on average, 50 seeds)
    copies = np.zeros(100)
    for seed in range(50):
        fx = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 150, seed=seed, target_rms=1.0, n_anchor=60, k=3, t_pos=0.0,
                                   positions=FLOOR, with_coupling=True)
        np.add.at(copies, np.r_[fx.rows[fx.stage == 0], fx.coupling["floor_exact_rows"]], 1)
    assert copies.min() >= 50 and copies.max() <= 50 * 2  # 1 to 2 copies per seed, mean 1.5
    assert abs(copies.mean() - 75.0) < 1e-9 and copies.std() < 0.25 * 75


def test_floor_exact_refusals():
    rng = np.random.default_rng(7)
    Cb, Ca, Eb, Ea, _, _ = _two_types(rng, n_b=100, n_a=100)
    kw = dict(seed=1, target_rms=1.0, n_anchor=40, k=3)
    for bad in ({}, {"t_pos": 0.5}, {"t_pos": 1.0}):
        with pytest.raises(ValueError, match="t_pos=0.0"):
            tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 50, positions=FLOOR, **bad, **kw)
    with pytest.raises(ValueError, match="herd"):
        tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 50, t_pos=0.0, positions=FLOOR, herd=True, **kw)
    with pytest.raises(ValueError, match="draw='midpoint'"):
        tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 50, t_pos=0.0, positions=FLOOR, draw=tb.DRAWS[1], **kw)
    with pytest.raises(ValueError, match="positions must be"):
        tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 50, positions="exact", **kw)


# ------------------------------------------------------------------------------------------------ the harness
def _ctx(seed=1_000_003):
    before, after, panel = _labelled_pair(np.random.default_rng(61))
    for i, a in enumerate((before, after)):  # two clouds of other shapes, so moved and real positions differ
        a.obsm[COORD_KEY] = (np.random.default_rng(i).normal(size=(a.n_obs, 3)) * [[1.0, 1.0, 0.5], [3.0, 1.0, 0.4]][i]
                             ).astype(np.float32)
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(panel), needs_coords=True, setting="embryo",
                  min_cells=50, max_cells=500)  # most of the 600 earlier cells: a sample RMS near the pool's
    view = harness.TaskView(key="syn:E47", task="T2", kind="interp", setting="embryo",
                            test_board="T2:embryo:test_interp", roles=("before", "after"),
                            role_days={"before": 7.5, "after": 8.0}, role_stems={"before": "E7.5", "after": "E8.0"},
                            reference_role="before", target_day=7.75)
    return harness.BuildContext(task=view, board=board, inputs={"before": before, "after": after}, seed=seed)


def test_the_candidates_are_registered_like_their_controls():
    for new, ctl in (("fgw_shape0x", "fgw_shape0"), ("fgw_shape0x_series", "fgw_shape0_series")):
        c = CANDIDATES[new]
        assert c.tasks == CANDIDATES[ctl].tasks and c.params["scope"] == CANDIDATES[ctl].params["scope"]
        assert c.params["positions"] == FLOOR and c.params["stream"] == tb.FLOOR_EXACT_STREAM
        tasks = {t.key: t for ts in harness.load_tasks().values() for t in ts}
        for key in ("T2e/E6.75+E7.5->E7.25", "T2e/E7.25+E8.0->E7.5"):
            assert c.applies(tasks[key]) == CANDIDATES[ctl].applies(tasks[key])
    assert CANDIDATES["fgw_shape0x_series"].params["rule"] == tb.SERIES_RULE


def test_each_arm_has_its_controls_expression_and_can_differ_from_it_in_position():
    ctx = _ctx()
    for new, ctl in (("fgw_shape0x", "fgw_shape0"), ("fgw_shape0x_series", "fgw_shape0_series")):
        a, b = CANDIDATES[new].build(ctx), CANDIDATES[ctl].build(ctx)
        np.testing.assert_array_equal(a.X, b.X)  # the same cells, in the same order, with the same expression
        assert a.coords.shape == b.coords.shape and not np.allclose(a.coords, b.coords)  # the control can differ
    # the plain arm's expression is fgw's mixture as well
    np.testing.assert_array_equal(CANDIDATES["fgw_shape0x"].build(ctx).X, CANDIDATES["fgw"].build(ctx).X)
    # the file's positions are the position mode's, at the harness's count, seed and size rule
    st = CANDIDATES["fgw_shape0x"].build(ctx)
    g = list(ctx.board.genes)
    b, a = ctx.inputs["before"], ctx.inputs["after"]
    Eb, Ea = (np.asarray(x[:, g].X.toarray(), np.float32) for x in (b, a))
    res = tb.interpolate_clouds(np.asarray(b.obsm[COORD_KEY], float), np.asarray(a.obsm[COORD_KEY], float), Eb, Ea,
                                0.5, harness.count_for(ctx.board, b.n_obs), ctx.seed,
                                target_rms=harness.target_rms(ctx), t_pos=0.0, positions=FLOOR)
    np.testing.assert_allclose(st.coords, res.coords, rtol=1e-6, atol=1e-4)
    # the same seed builds the same file; another seed another draw
    np.testing.assert_array_equal(CANDIDATES["fgw_shape0x"].build(ctx).coords, st.coords)
    assert not np.array_equal(CANDIDATES["fgw_shape0x"].build(_ctx(1_000_004)).coords, st.coords)


# ------------------------------------------------------------------------------------------------ the builder
@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


def test_the_builder_option_is_off_by_default_and_refused_where_not_defined(bt):
    ap = bt.build_parser()
    base = ["embryo", "--method", "fgw_shape0", "--out", "runs/x/T2embryo_s2.h5ad"]
    assert ap.parse_args(base).positions == tb.POSITIONS[0]
    assert ap.parse_args([*base, "--positions", FLOOR]).positions == FLOOR
    src = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert "positions=a.positions)" in src and src.count('default="mixture"') == 1
    for argv in (["embryo", "--method", "fgw", "--positions", FLOOR],
                 ["embryo", "--method", "fgw_shape0", "--positions", FLOOR, "--herd"],
                 ["embryo", "--method", "fgw_shape0", "--positions", FLOOR, "--draw", "midpoint"]):
        a = ap.parse_args([*argv, "--out", "runs/x/T2embryo_x.h5ad"])
        with pytest.raises(SystemExit, match="--positions"):
            bt.build_embryo(a)  # refused before any board or stage is read
