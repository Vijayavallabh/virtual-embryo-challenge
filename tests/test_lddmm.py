"""vec.lddmm (row E53, as specified): geodesic shooting, the stored field, point transport, the Jacobian, the data term,
the k-means++ control points and the fit, on synthetic arrays only (CPU, float64). Nothing is read from data/.

Pinned: the row's constants; the geodesic equations match the Hamiltonian's autograd gradient; zero momenta give the
identity; carrying forward and back returns the points; the Jacobian determinant of the identity is 1 and of a
known stretch is its factor's product; the data term is zero between a cloud and itself and positive otherwise; the
fit at a fixed seed is deterministic, starts from the rigid data term and lowers it; the cache fits once per input."""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from vec import lddmm  # noqa: E402


def T(a):
    return torch.as_tensor(np.asarray(a, np.float64))


def _clouds(n=300, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(n, 3)) * [1.4, 1.0, 0.5]
    B = rng.normal(size=(n + 40, 3)) * [1.1, 1.1, 0.8] + [0.2, 0.0, 0.0]
    A, B = ((X - X.mean(0)) / np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()) for X in (A, B))
    FA = np.c_[A[:, :1], rng.normal(size=(n, 2))]
    FB = np.c_[B[:, :1], rng.normal(size=(n + 40, 2))]
    return A, B, FA, FB


def test_the_rows_constants():
    assert (lddmm.SUBSAMPLE, lddmm.N_CONTROL, lddmm.SIGMA_V) == (4000, 1500, 0.2)
    assert lddmm.DATA_SCALES == (0.05, 0.1, 0.2, 0.4) and lddmm.GAMMA == 10.0
    assert (lddmm.N_STEPS, lddmm.LBFGS_ITERS) == (10, 150)
    assert (lddmm.SUBSAMPLE_STREAM, lddmm.CONTROL_STREAM) == (51, 52)


def test_the_geodesic_equations_are_the_hamiltonians_gradient():
    rng = np.random.default_rng(1)
    q = T(rng.normal(size=(20, 3)) * 0.3).requires_grad_(True)
    p = T(rng.normal(size=(20, 3)) * 0.1).requires_grad_(True)
    H = 0.5 * (p * (lddmm._gauss(q, q, 0.2) @ p)).sum()
    gq, gp = torch.autograd.grad(H, (q, p))
    dq, dp = lddmm.hamiltonian_rhs(q.detach(), p.detach(), 0.2)
    np.testing.assert_allclose(dq.numpy(), gp.numpy(), atol=1e-12)
    np.testing.assert_allclose(dp.numpy(), -gq.numpy(), atol=1e-12)


def test_the_shooting_keeps_the_hamiltonian_and_stores_the_half_nodes():
    rng = np.random.default_rng(2)
    q, p = T(rng.normal(size=(15, 3)) * 0.4), T(rng.normal(size=(15, 3)) * 0.05)
    qs, ps = lddmm.shoot(q, p, 0.2, 10)
    assert len(qs) == len(ps) == 21

    def H(a, b):
        return float(0.5 * (b * (lddmm._gauss(a, a, 0.2) @ b)).sum())

    assert H(qs[-1], ps[-1]) == pytest.approx(H(q, p), rel=1e-6)
    s = lddmm._state(qs, ps, 0.35, 10)  # halfway between nodes 7 (0.35) exactly: the stored node
    assert torch.equal(s[0], qs[7])
    s = lddmm._state(qs, ps, 0.375, 10)
    np.testing.assert_allclose(s[0].numpy(), (0.5 * (qs[7] + qs[8])).numpy())


def test_zero_momenta_give_the_identity_and_its_jacobian_is_one():
    X = np.random.default_rng(3).normal(size=(50, 3))
    fit = lddmm.GeodesicFit(q0=X[:10].copy(), p0=np.zeros((10, 3)), record={})
    np.testing.assert_array_equal(fit.carry(X, 0.0, 0.5), X)
    np.testing.assert_allclose(fit.jacobian_det(X, 1.0, 0.3), 1.0, atol=1e-12)


def test_forward_then_back_returns_the_points_and_the_cuts_follow_the_grid():
    rng = np.random.default_rng(4)
    q0 = rng.normal(size=(30, 3)) * 0.5
    fit = lddmm.GeodesicFit(q0=q0, p0=rng.normal(size=(30, 3)) * 0.05, record={})
    X = rng.normal(size=(40, 3)) * 0.5
    Y = fit.carry(X, 0.0, 1.0)
    assert np.abs(Y - X).max() > 1e-3  # it moves
    np.testing.assert_allclose(fit.carry(Y, 1.0, 0.0), X, atol=1e-6)
    np.testing.assert_allclose(fit.carry(fit.carry(X, 0.0, 0.37), 0.37, 1.0), Y, atol=1e-6)
    np.testing.assert_allclose(fit.carry(Y, 1.0, 0.5), fit.carry(X, 0.0, 0.5), atol=1e-6)
    assert lddmm._cuts(0.0, 0.35, 10) == [0.0, pytest.approx(0.1), pytest.approx(0.2), pytest.approx(0.3), 0.35]
    assert lddmm._cuts(1.0, 0.5, 10)[:2] == [1.0, pytest.approx(0.9)] and lddmm._cuts(1.0, 0.5, 10)[-1] == 0.5


def test_the_jacobian_matches_finite_differences():
    rng = np.random.default_rng(5)
    fit = lddmm.GeodesicFit(q0=rng.normal(size=(25, 3)) * 0.4, p0=rng.normal(size=(25, 3)) * 0.1, record={})
    x = rng.normal(size=(1, 3)) * 0.3
    eps = 1e-6
    J = np.stack([(fit.carry(x + eps * e, 0.0, 0.5) - fit.carry(x - eps * e, 0.0, 0.5))[0] / (2 * eps)
                  for e in np.eye(3)], axis=1)
    assert fit.jacobian_det(x, 0.0, 0.5)[0] == pytest.approx(np.linalg.det(J), rel=1e-5)


def test_kmeanspp_is_seeded_distinct_and_spreads():
    X = np.r_[np.zeros((100, 3)), np.ones((100, 3)) * 10, np.full((1, 3), -10.0)]
    a = lddmm.kmeanspp(X, 3, np.random.default_rng([7, 52]))
    b = lddmm.kmeanspp(X, 3, np.random.default_rng([7, 52]))
    np.testing.assert_array_equal(a, b)
    assert len(set(a.tolist())) == 3
    assert {tuple(X[i]) for i in a} == {(0.0, 0.0, 0.0), (10.0, 10.0, 10.0), (-10.0, -10.0, -10.0)}
    np.testing.assert_array_equal(lddmm.kmeanspp(X[:5], 9, np.random.default_rng(0)), np.arange(5))
    assert len(set(lddmm.kmeanspp(np.zeros((6, 3)), 4, np.random.default_rng(0)).tolist())) == 4


def test_the_data_term_is_a_squared_mmd_with_the_median_bandwidth():
    A, B, FA, FB = _clouds(120)
    d = lddmm._Data(T(A), T(B), T(FA), T(FB))
    F = np.r_[FA, FB]
    D = np.sqrt(((F[:, None] - F[None]) ** 2).sum(-1))[np.triu_indices(len(F), 1)]
    assert d.bandwidth == pytest.approx(np.sort(D)[(len(D) - 1) // 2], rel=1e-9)  # torch's median: the lower one
    assert float(d(T(A))) > 0
    same = lddmm._Data(T(A), T(A), T(FA), T(FA))
    assert abs(float(same(T(A)))) < 1e-12
    with pytest.raises(ValueError, match="no spread"):
        lddmm._Data(T(A), T(B), T(np.zeros_like(FA)), T(np.zeros_like(FB)))


def test_the_fit_is_deterministic_lowers_the_data_term_and_records_flat_facts():
    A, B, FA, FB = _clouds(200)
    kw = {"n_sub": 150, "n_control": 40, "iters": 15, "device": "cpu"}
    f1 = lddmm.fit(A, B, FA, FB, 1_000_003, **kw)
    f2 = lddmm.fit(A, B, FA, FB, 1_000_003, **kw)
    np.testing.assert_array_equal(f1.p0, f2.p0)
    np.testing.assert_array_equal(f1.q0, f2.q0)
    r = f1.record
    assert r["data_ratio"] == pytest.approx(r["data_final"] / r["data_rigid"]) and r["data_ratio"] < 1.0
    assert r["device"] == "cpu" and 0 < r["lbfgs_n_iter"] <= 15
    assert all(not isinstance(v, (list, dict, tuple)) for v in r.values())  # flat values, no count
    # control points are rows of the earlier subsample (stream [seed, 51] then k-means++ on stream [seed, 52])
    rows = lddmm.subsample_rows(len(A), 150, np.random.default_rng([1_000_003, 51]))
    ic = lddmm.kmeanspp(A[rows], 40, np.random.default_rng([1_000_003, 52]))
    np.testing.assert_array_equal(f1.q0, A[rows][ic])
    f3 = lddmm.fit(A, B, FA, FB, 1_000_004, **kw)
    assert not np.array_equal(f3.q0, f1.q0)
    with pytest.raises(ValueError, match="one feature row per cell"):
        lddmm.fit(A, B, FA[:-1], FB, 1, **kw)


def test_a_repeated_request_answers_with_the_same_array(monkeypatch):
    """The two scorer builds of one build seed carry the same rows: the second request is the first's answer."""
    rng = np.random.default_rng(6)
    fit = lddmm.GeodesicFit(q0=rng.normal(size=(10, 3)), p0=rng.normal(size=(10, 3)) * 0.05, record={})
    X = rng.normal(size=(30, 3))
    first, jac = fit.carry(X, 0.0, 0.5), fit.jacobian_det(X, 1.0, 0.5)
    calls = []
    monkeypatch.setattr(lddmm, "carry", lambda *a, **k: calls.append(1))
    np.testing.assert_array_equal(fit.carry(X.copy(), 0.0, 0.5), first)
    np.testing.assert_array_equal(fit.jacobian_det(X, 1.0, 0.5), jac)
    assert calls == []
    fit.carry(X, 0.0, 0.5)[0, 0] = 99.0  # an answer is a copy
    np.testing.assert_array_equal(fit.carry(X, 0.0, 0.5), first)


def test_the_cache_fits_once_per_input_set(monkeypatch):
    calls = []

    def fake(Zb, Zt, Fb, Ft, seed, **kw):
        calls.append((seed, kw["iters"]))
        return lddmm.GeodesicFit(q0=np.zeros((1, 3)), p0=np.zeros((1, 3)), record={"seed": seed})

    monkeypatch.setattr(lddmm, "fit", fake)
    monkeypatch.setattr(lddmm, "_CACHE", {})
    A, B, FA, FB = _clouds(50)
    a = lddmm.fit_cached(A, B, FA, FB, 1)
    assert lddmm.fit_cached(A, B, FA, FB, 1) is a and calls == [(1, lddmm.LBFGS_ITERS)]
    lddmm.fit_cached(A, B + 1e-9, FA, FB, 1)
    lddmm.fit_cached(A, B, FA, FB, 2)
    assert len(calls) == 3
    for s in range(3, 3 + lddmm.CACHE_SIZE):
        lddmm.fit_cached(A, B, FA, FB, s)
    assert len(lddmm._CACHE) == lddmm.CACHE_SIZE
