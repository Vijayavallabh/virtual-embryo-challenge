"""Row E54 (the multi-stage flow, as specified): vec.t2_mmflow on synthetic arrays only. The numpy parts (thinning, PCA,
folds, the spline, chained OT, the training batch, retrieval, Kish pool use, the energy distance, the gate's draws and
rule) run everywhere; the torch parts (the net, the fit, dopri5) need the uv gpu group and are skipped without it (the
launcher's tests step refuses to run without torch, so on the DGX they run). No data is read."""

from __future__ import annotations

import itertools

import numpy as np
import pytest
from scipy import sparse

from vec import t2_mmflow as mf


def test_constants_are_the_specifications():
    assert (mf.THIN_DEPTH, mf.THIN_STREAM, mf.PCA_STREAM, mf.GATE_STREAM) == (30, 53, 54, 55)
    assert (mf.N_PC, mf.PCA_CELLS, mf.HIDDEN, mf.N_HIDDEN) == (30, 3000, 256, 3)
    assert (mf.STEPS, mf.BATCH, mf.LR, mf.TORCH_SEED, mf.RTOL) == (20_000, 1024, 1e-3, 20261122, 1e-5)
    assert (mf.GATE_DRAWS, mf.GATE_CELLS) == (5, 5000)
    assert mf.DIRECTIONS == ("snap", "anti")
    assert mf.stage_key(8.3) == 830 and mf.stage_key(6.75) == 675


# ------------------------------------------------------------------------------------------------ thinning and PCA
def _counts(rng, n=400, g=12, lam=6.0):
    C = rng.poisson(lam * rng.gamma(1.0, 1.0, size=(n, 1)), size=(n, g)).astype(np.int32)
    return sparse.csr_matrix(C)


def test_thinning_keeps_shallow_cells_and_thins_deep_cells_to_the_depth_on_average():
    rng = np.random.default_rng(1)
    C = _counts(rng)
    T = mf.thin_counts(C, 7.5)
    tot, after = np.asarray(C.sum(1)).ravel(), np.asarray(T.sum(1)).ravel()
    shallow = tot <= mf.THIN_DEPTH
    np.testing.assert_array_equal(T[shallow].toarray(), C[shallow].toarray())  # p = 1: every count kept
    assert np.all(T.toarray() <= C.toarray()) and np.all(T.data > 0)
    deep = tot > 3 * mf.THIN_DEPTH
    assert deep.sum() > 50 and abs(after[deep].mean() - mf.THIN_DEPTH) < 1.5


def test_thinning_draws_on_the_stages_own_stream():
    rng = np.random.default_rng(2)
    C = _counts(rng)
    a, b, c = mf.thin_counts(C, 8.0), mf.thin_counts(C, 8.0), mf.thin_counts(C, 7.5)
    assert (a != b).nnz == 0 and (a != c).nnz > 0
    # the stream [round(100 day), 53], every count drawn in canonical CSR order with p = min(1, 30 / total)
    Cs = sparse.csr_matrix(C, dtype=np.float64)
    Cs.sort_indices()
    tot = np.asarray(Cs.sum(1)).ravel()
    p = np.minimum(1.0, mf.THIN_DEPTH / np.maximum(tot, 1e-300))
    want = np.random.default_rng([800, 53]).binomial(Cs.data.astype(np.int64), np.repeat(p, np.diff(Cs.indptr)))
    np.testing.assert_array_equal(sparse.csr_matrix((want, Cs.indices, Cs.indptr), shape=Cs.shape).toarray(),
                                  a.toarray())


def test_the_pca_takes_the_same_number_of_cells_per_stage_on_stream_54():
    rng = np.random.default_rng(3)
    feats = {6.75: rng.normal(size=(50, 8)), 7.25: rng.normal(size=(80, 8)) + 1.0, 7.5: rng.normal(size=(60, 8))}
    centre, L = mf.fit_embedding(feats, n_pc=5, per_stage=1000)
    rows = {d: np.sort(np.random.default_rng([mf.stage_key(d), 54]).choice(len(F), 50, replace=False))
            for d, F in feats.items()}
    X = np.vstack([feats[d][rows[d]] for d in sorted(feats)])
    np.testing.assert_allclose(centre, X.mean(0))
    assert L.shape == (8, 5)
    np.testing.assert_allclose(L.T @ L, np.eye(5), atol=1e-10)
    _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    np.testing.assert_allclose(np.abs(L), np.abs(Vt[:5].T), atol=1e-10)
    np.testing.assert_allclose(mf.embed(X, centre, L), (X - centre) @ L)
    c2, _ = mf.fit_embedding(feats, n_pc=5, per_stage=20)  # fewer than every stage holds: 20 from each
    rows20 = {d: np.sort(np.random.default_rng([mf.stage_key(d), 54]).choice(len(F), 20, replace=False))
              for d, F in feats.items()}
    np.testing.assert_allclose(c2, np.vstack([feats[d][rows20[d]] for d in sorted(feats)]).mean(0))


# ------------------------------------------------------------------------------------------------ folds
RELEASED = (6.75, 7.25, 7.5, 8.0, 8.3, 8.5, 8.75)


@pytest.mark.parametrize("b, a, t, want", [
    (7.5, 8.5, 8.0, {8.0, 8.3}), (7.5, 8.3, 8.0, {8.0}), (7.5, 8.5, 8.3, {8.0, 8.3}), (7.5, 8.75, 8.5, {8.0, 8.3, 8.5}),
    (8.0, 8.75, 8.3, {8.3, 8.5}), (8.0, 8.75, 8.5, {8.3, 8.5}), (8.3, 8.75, 8.5, {8.5}), (7.25, 8.0, 7.5, {7.5}),
    (6.75, 7.5, 7.25, {7.25})])
def test_exclusion_is_the_target_and_every_stage_strictly_inside_the_bracket(b, a, t, want):
    assert mf.exclusion(b, a, t, RELEASED) == frozenset(want)


def test_neighbours():
    assert mf.neighbours(7.5, RELEASED) == (7.25, 8.0) and mf.neighbours(8.3, RELEASED) == (8.0, 8.5)
    with pytest.raises(ValueError):
        mf.neighbours(6.75, RELEASED)


# ------------------------------------------------------------------------------------------------ spline and batch
def _spline_at(y, knots, t):
    return mf.natural_spline3(y, knots, t)


def test_the_spline_passes_its_knots_and_is_natural_and_smooth():
    rng = np.random.default_rng(4)
    B, D = 7, 3
    y = rng.normal(size=(3, B, D))
    knots = np.stack([np.full(B, 7.25), np.full(B, 7.5), np.full(B, 8.0)]) + rng.uniform(0, 0.01, (1, B))
    for j in range(3):
        pos, _ = _spline_at(y, knots, knots[j])
        np.testing.assert_allclose(pos, y[j], atol=1e-12)
    h = 1e-4

    def second(t):
        return (_spline_at(y, knots, t + h)[0] - 2 * _spline_at(y, knots, t)[0] + _spline_at(y, knots, t - h)[0]) / h**2

    # natural ends: S'' is linear on each end segment, so its value at the end knot is extrapolated from two points
    np.testing.assert_allclose(2 * second(knots[0] + 20 * h) - second(knots[0] + 40 * h), 0.0, atol=1e-3)
    np.testing.assert_allclose(2 * second(knots[2] - 20 * h) - second(knots[2] - 40 * h), 0.0, atol=1e-3)
    assert np.abs(second(knots[1] - 20 * h)).max() > 0.1  # and not zero in the middle (the knots are not collinear)
    vl, vr = _spline_at(y, knots, knots[1] - 1e-9)[1], _spline_at(y, knots, knots[1] + 1e-9)[1]
    np.testing.assert_allclose(vl, vr, atol=1e-5)  # the velocity is continuous at the middle knot
    t = knots[0] + 0.3 * (knots[2] - knots[0])
    fd = (_spline_at(y, knots, t + 1e-6)[0] - _spline_at(y, knots, t - 1e-6)[0]) / 2e-6
    np.testing.assert_allclose(_spline_at(y, knots, t)[1], fd, atol=1e-5)  # the velocity is the derivative


def test_points_on_a_line_at_their_own_times_move_at_one_velocity():
    B = 5
    knots = np.stack([np.full(B, 7.0), np.full(B, 7.5), np.full(B, 8.25)])
    v = np.arange(B * 2, dtype=float).reshape(B, 2)
    y = np.stack([v * (k[0] - 7.0) for k in knots])
    t = np.linspace(7.0, 8.25, B)
    pos, vel = mf.natural_spline3(y, knots, t)
    np.testing.assert_allclose(vel, v, atol=1e-12)
    np.testing.assert_allclose(pos, v * (t - 7.0)[:, None], atol=1e-12)


def test_chained_ot_is_the_exact_assignment_between_consecutive_stages():
    rng = np.random.default_rng(5)
    a = rng.normal(size=(6, 2))
    perm = rng.permutation(6)
    b = a[perm] + 0.001
    c = b[rng.permutation(6)] + 0.001
    paths = mf.chain([a, b, c])
    np.testing.assert_allclose(paths[1], a + 0.001, atol=1e-12)  # each a finds its own shifted copy
    np.testing.assert_allclose(paths[2], a + 0.002, atol=1e-12)
    x, y = rng.normal(size=(5, 2)), rng.normal(size=(5, 2))  # against brute force
    best = min(itertools.permutations(range(5)), key=lambda p: ((x - y[list(p)]) ** 2).sum())
    np.testing.assert_allclose(mf.chain([x, y])[1], y[list(best)])
    with pytest.raises(ValueError):
        mf.chain([x, y[:4]])


def test_the_training_batch_takes_a_triple_and_a_time_inside_its_span():
    rng = np.random.default_rng(6)
    days = [6.75, 7.25, 7.5, 8.0]
    lat = [rng.normal(size=(40, 3)) + d for d in days]
    z, t, u = mf.training_batch(lat, days, 16, np.random.default_rng(9))
    assert z.shape == (16, 3) and t.shape == (16,) and u.shape == (16, 3)
    assert t.min() >= 6.75 and t.max() <= 8.0
    z2, t2, u2 = mf.training_batch(lat, days, 16, np.random.default_rng(9))
    np.testing.assert_array_equal(z, z2)
    np.testing.assert_array_equal(u, u2)
    with pytest.raises(ValueError):
        mf.training_batch(lat[:2], days[:2], 4, rng)
    with pytest.raises(ValueError):
        mf.training_batch(lat, [7.0, 6.0, 8.0, 9.0], 4, rng)


# ------------------------------------------------------------------------------------------------ retrieval
def test_retrieval_is_one_to_one_and_exact():
    rng = np.random.default_rng(7)
    pool = rng.normal(size=(9, 2))
    np.testing.assert_array_equal(mf.retrieve(pool[[4, 1, 7]], pool), [4, 1, 7])
    pred = rng.normal(size=(3, 2))
    rows = mf.retrieve(pred, pool)
    assert len(set(rows.tolist())) == 3

    def cost(r):
        return ((pred - pool[list(r)]) ** 2).sum()

    best = min(itertools.permutations(range(9), 3), key=cost)
    assert cost(rows) == pytest.approx(cost(best))
    two = np.vstack([pool[2], pool[2]])  # two states at one cell: still two distinct cells
    assert len(set(mf.retrieve(two, pool).tolist())) == 2
    with pytest.raises(mf.MmflowRefused):
        mf.retrieve(rng.normal(size=(10, 2)), pool)
    assert len(mf.retrieve(np.empty((0, 2)), pool)) == 0


def test_kish_pool_use():
    assert mf.kish_ess_frac(np.arange(5), 10) == pytest.approx(0.5)
    assert mf.kish_ess_frac(np.array([0, 0, 1, 1]), 10) == pytest.approx(0.2)
    assert mf.kish_ess_frac(np.array([0, 0, 0, 1]), 4) == pytest.approx(16 / 10 / 4)
    assert mf.kish_ess_frac(np.array([], dtype=int), 4) is None


# ------------------------------------------------------------------------------------------------ gate
def test_energy_distance():
    rng = np.random.default_rng(8)
    X = rng.normal(size=(300, 3))
    assert mf.energy_distance(X, X) == pytest.approx(0.0, abs=1e-12)
    Y = rng.normal(size=(250, 3)) + 1.0
    d = mf.energy_distance(X, Y, chunk=37)
    assert d > 0.5 and d == pytest.approx(mf.energy_distance(Y, X))
    from scipy.spatial.distance import cdist

    want = 2 * cdist(X, Y).mean() - cdist(X, X).mean() - cdist(Y, Y).mean()
    assert d == pytest.approx(want, rel=1e-10)


def test_gate_draws_follow_the_calendar_fraction_on_stream_55():
    lo, hi, tr = mf.gate_draws(1000, 2000, 3000, 0.25, 3, n=400)
    assert (len(lo), len(hi), len(tr)) == (300, 100, 400)
    rng = np.random.default_rng([3, 55])
    np.testing.assert_array_equal(lo, rng.choice(1000, 300, replace=False))
    np.testing.assert_array_equal(hi, rng.choice(2000, 100, replace=False))
    np.testing.assert_array_equal(tr, rng.choice(3000, 400, replace=False))
    lo, hi, tr = mf.gate_draws(10, 10, 10, 0.5, 0, n=40)  # small pools: with replacement
    assert (len(lo), len(hi), len(tr)) == (20, 20, 40)


def test_the_gate_stage_prefers_states_that_reach_the_held_stage():
    rng = np.random.default_rng(9)
    lo, hi, held = rng.normal(size=(800, 2)), rng.normal(size=(800, 2)) + 4.0, rng.normal(size=(800, 2)) + 2.0
    good = mf.gate_stage(lo, lo + 2.0, hi, hi - 2.0, held, 0.5, draws=2, n=200)
    assert good["beats"] and good["flow_mean"] < good["mixture_mean"] and len(good["flow"]) == 2
    bad = mf.gate_stage(lo, lo - 2.0, hi, hi + 2.0, held, 0.5, draws=2, n=200)
    assert not bad["beats"]


@pytest.mark.parametrize("beats, want", [
    ((1, 1, 1, 1, 1), True), ((1, 1, 1, 1, 0), True), ((0, 1, 1, 1, 1), True), ((1, 0, 1, 1, 1), False),
    ((1, 1, 1, 0, 0), False)])
def test_the_gate_rule_is_four_of_five_and_e75(beats, want):
    names = ("E7.25", "E7.5", "E8.0", "E8.25_late", "E8.5")
    stages = {n: {"beats": bool(b)} for n, b in zip(names, beats, strict=True)}
    assert mf.gate_pass(stages, "E7.5", 4) is want


# ------------------------------------------------------------------------------------------------ torch (gpu group)
def test_the_net_is_three_hidden_layers_of_256_silu_on_z_and_t():
    torch = pytest.importorskip("torch")
    net = mf.velocity_net(30)
    lin = [m for m in net if isinstance(m, torch.nn.Linear)]
    assert [(m.in_features, m.out_features) for m in lin] == [(31, 256), (256, 256), (256, 256), (256, 30)]
    assert sum(isinstance(m, torch.nn.SiLU) for m in net) == 3


def test_the_fit_is_deterministic_and_learns_a_drift():
    torch = pytest.importorskip("torch")
    pytest.importorskip("torchdiffeq")
    rng = np.random.default_rng(10)
    days = [7.0, 7.5, 8.0, 8.5]
    lat = [rng.normal(size=(200, 2)) * 0.1 + np.array([2.0 * (d - 7.0), 0.0]) for d in days]
    a, la = mf.fit_flow(lat, days, steps=300, batch=64, lr=1e-2)
    b, lb = mf.fit_flow(lat, days, steps=300, batch=64, lr=1e-2)
    np.testing.assert_array_equal(la, lb)
    for p, q in zip(a.state_dict().values(), b.state_dict().values(), strict=True):
        assert torch.equal(p, q)
    assert la[-50:].mean() < la[:50].mean()
    out = mf.integrate(a, lat[1], 7.5, 8.0)  # the drift is +2 per day in the first component
    assert abs(out[:, 0].mean() - lat[1][:, 0].mean() - 1.0) < 0.3
    back = mf.integrate(a, lat[1], 7.5, 7.0)
    assert out[:, 0].mean() > lat[1][:, 0].mean() > back[:, 0].mean()


def test_dopri5_carries_a_constant_field_exactly_and_leaves_t0_equal_t1_alone():
    torch = pytest.importorskip("torch")
    pytest.importorskip("torchdiffeq")
    net = mf.velocity_net(3)
    last = [m for m in net if isinstance(m, torch.nn.Linear)][-1]
    with torch.no_grad():
        last.weight.zero_()
        last.bias.copy_(torch.tensor([1.0, -2.0, 0.5]))
    z0 = np.random.default_rng(11).normal(size=(20, 3))
    np.testing.assert_allclose(mf.integrate(net, z0, 7.5, 7.75), z0 + 0.25 * np.array([1.0, -2.0, 0.5]), atol=1e-5)
    np.testing.assert_allclose(mf.integrate(net, z0, 8.0, 7.75), z0 - 0.25 * np.array([1.0, -2.0, 0.5]), atol=1e-5)
    np.testing.assert_array_equal(mf.integrate(net, z0, 7.5, 7.5), z0.astype(np.float32).astype(np.float64))
