"""Tests for experiments/t3_ext/context_map.py (row K26's features, map and stage 0's statistics). Synthetic arrays
only: no competition file, no outside table, no score."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
T3X = REPO / "experiments" / "t3_ext"
if str(T3X) not in sys.path:
    sys.path.insert(0, str(T3X))
import context_map as cmap  # noqa: E402


def _expr(rng, n_per_type=(300, 150, 40), n_genes=30, coupled=(("q", "p0", 0.9),)):
    """Synthetic pooled cells: genes p0..p{n-1} and q, three types (the last under the 100-cell floor), a library-size
    trend shared by every gene, and q coupled to one panel gene inside each type."""
    types, cols, lnN = [], [], []
    names = [f"p{i}" for i in range(n_genes)] + ["q", "rare"]
    for t, n in enumerate(n_per_type):
        ln_tot = rng.normal(9.0, 0.4, n)
        base = rng.normal(0, 1, (n, n_genes)) + 0.8 * ln_tot[:, None] + t
        q = 0.8 * ln_tot + rng.normal(0, 1, n)
        for _, p, w in coupled:
            j = int(p[1:])
            q = q + w * (base[:, j] - 0.8 * ln_tot)
        rare = np.zeros(n)
        X = np.column_stack([base, q, rare])
        cols.append(np.abs(X))
        types += [f"T{t}"] * n
        lnN.append(ln_tot)
    X = np.vstack(cols)
    X[:4, -1] = 1.0  # 'rare' detected in 4 cells only
    return cmap.Expression(X=sparse.csc_matrix(X), ln_total=np.concatenate(lnN), types=np.array(types), genes=names)


def test_type_correlations_remove_the_library_size_trend_and_mark_constant_columns():
    rng = np.random.default_rng(1)
    n = 400
    lt = rng.normal(9, 0.5, n)
    a = 2.0 * lt + rng.normal(0, 1, n)  # shares only the trend with b
    b = -3.0 * lt + rng.normal(0, 1, n)
    shared = rng.normal(0, 1, n)
    c = lt + shared
    d = 5 * lt + 2 * shared  # shares the residual with c
    C = cmap.type_correlations(np.column_stack([a, c, np.ones(n)]), np.column_stack([b, d]), lt)
    assert abs(C[0, 0]) < 0.15  # the trend alone gives no correlation after the regression
    assert C[1, 1] > 0.99
    assert np.all(np.isnan(C[2]))  # a constant query column has no residual: undefined


def test_raw_contexts_average_over_kept_types_weighted_by_cell_count():
    rng = np.random.default_rng(2)
    e = _expr(rng)
    panel = [f"p{i}" for i in range(30)]
    raw, rec = cmap.raw_contexts(e, ["q"], panel)
    assert rec["types"] == {"T0": 300, "T1": 150}  # T2 (40 cells) is under the floor
    want = np.zeros(30)
    for _t, rows in cmap.kept_types(e.types):
        D = e.columns([e.pos["q"]] + [e.pos[g] for g in panel])[rows]
        want += cmap.type_correlations(D[:, :1], D[:, 1:], e.ln_total[rows])[0] * rows.size
    np.testing.assert_allclose(raw[0], want / 450, atol=1e-12)
    assert raw[0, 0] > 0.4 and abs(raw[0, 5]) < 0.25


def test_a_panel_gene_absent_from_t1_is_nan_then_zero_after_the_finish():
    rng = np.random.default_rng(3)
    e = _expr(rng)
    panel = ["p0", "absent", *[f"p{i}" for i in range(1, 30)]]
    ctx, rec = cmap.contexts(e, ["q"], panel)
    v = ctx["q"]
    assert v[1] == 0.0 and rec["panel_absent_from_t1"] == 1
    ok = np.arange(len(panel)) != 1
    assert abs(v[ok].mean()) < 1e-9 and abs(v[ok].std() - 1.0) < 1e-9
    assert v[0] == v.max()  # the coupled panel gene stands out


def test_finish_context_residual_is_orthogonal_to_level_and_detection_and_z_scored():
    rng = np.random.default_rng(4)
    p = 60
    level, det = rng.uniform(0, 3, p), rng.uniform(0, 1, p)
    raw = 0.5 * level - 0.7 * det + rng.normal(0, 0.1, p)
    raw[7] = np.nan
    out = cmap.finish_context(raw, level, det)[0]
    ok = np.isfinite(raw)
    assert out[7] == 0.0
    assert abs(out[ok] @ (level[ok] - level[ok].mean())) < 1e-8 and abs(out[ok] @ (det[ok] - det[ok].mean())) < 1e-8
    assert abs(out[ok].mean()) < 1e-12 and abs(out[ok].std() - 1) < 1e-12


def test_a_gene_absent_or_under_the_detection_floor_is_refused():
    rng = np.random.default_rng(5)
    e = _expr(rng)
    panel = [f"p{i}" for i in range(30)]
    with pytest.raises(cmap.ContextRefused, match="not in the T1 gene list"):
        cmap.contexts(e, ["nope"], panel)
    det = cmap.detection(e)
    assert det[e.pos["rare"]] < cmap.DETECTION_FLOOR
    with pytest.raises(cmap.ContextRefused, match="fewer than 1%"):
        cmap.contexts(e, ["rare"], panel)


def test_detection_deciles_cover_only_genes_at_the_floor_in_ten_near_equal_runs():
    rng = np.random.default_rng(6)
    det = np.concatenate([rng.uniform(0.0, 0.009, 13), rng.uniform(0.01, 1.0, 97)])
    genes = [f"g{i}" for i in range(det.size)]
    dec = cmap.decile_of(det, genes)
    assert set(dec) == {genes[i] for i in range(13, 110)}
    sizes = np.bincount(list(dec.values()))
    assert sizes.size == 10 and sizes.max() - sizes.min() <= 1
    by = {b: [det[genes.index(g)] for g, x in dec.items() if x == b] for b in range(10)}
    assert all(max(by[b]) <= min(by[b + 1]) for b in range(9))


def test_draws_stay_in_the_decile_avoid_the_excluded_and_repeat_under_the_seed():
    det = np.linspace(0.02, 1.0, 200)
    genes = [f"g{i}" for i in range(200)]
    dec = cmap.decile_of(det, genes)
    excl = {"g3", "g5"}
    a = cmap.draw_same_decile(dec, "g0", 15, np.random.default_rng(20261122), excl)
    b = cmap.draw_same_decile(dec, "g0", 15, np.random.default_rng(20261122), excl)
    assert a == b and len(set(a)) == 15
    assert all(dec[g] == dec["g0"] for g in a) and not (set(a) & (excl | {"g0"}))
    with pytest.raises(cmap.ContextRefused, match="fewer than"):
        cmap.draw_same_decile(dec, "g0", 19, np.random.default_rng(0), excl)


def _latent(rng, n, p, d=4, A=None):
    """Contexts that live near a d-dimensional space shared by every gene (as co-expression profiles do), z-scored per
    row; returns (X, Z, A)."""
    A = rng.normal(0, 1, (d, p)) if A is None else A
    Z = rng.normal(0, 1, (n, d))
    X = Z @ A + 0.1 * rng.normal(0, 1, (n, p))
    return (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True), Z, A


def _linear_toy(rng, n, p=500, g=80, noise=0.05):
    """A response that is a known linear function of the context: Y = X W + 0.3 + noise, W of rank 2."""
    X, _, _ = _latent(rng, n, p)
    W = rng.normal(0, 1, (p, 2)) @ rng.normal(0, 1, (2, g)) / p
    Y = X @ W + 0.3 + noise * rng.normal(0, 1, (n, g)) * float(np.std(X @ W))
    return X, Y, W


def test_the_map_recovers_a_response_that_is_a_linear_function_of_the_context():
    rng = np.random.default_rng(7)
    X, Y, W = _linear_toy(rng, 41)
    m = cmap.fit_map(X[:40], Y[:40])
    assert m.rank == 2 and m.lam == 400.0 and m.n == 40
    assert np.linalg.matrix_rank(m.B, tol=1e-9) == 2
    pred = m.predict(X[40])
    truth = X[40] @ W + 0.3
    assert np.corrcoef(pred - m.y_mean, truth - m.y_mean)[0, 1] > 0.9


def test_the_map_returns_about_the_mean_when_the_context_is_pure_noise():
    rng = np.random.default_rng(8)
    n, p, g = 16, 500, 80
    X = rng.normal(0, 1, (n + 50, p))
    Y = rng.normal(0, 1, (n + 50, g)) + rng.normal(0, 1, g)  # a shared mean table, no link to X
    m = cmap.fit_map(X[:n], Y[:n])
    dev = np.array([np.linalg.norm(m.predict(x) - m.y_mean) for x in X[n:]])
    spread = np.array([np.linalg.norm(y - m.y_mean) for y in Y[n:]])
    assert np.median(dev / spread) < 0.3  # the prediction stays near the training mean table
    r = [np.corrcoef(m.predict(X[i]) - m.y_mean, Y[i] - m.y_mean)[0, 1] for i in range(n, n + 50)]
    assert abs(np.mean(r)) < 0.1


def test_the_intercept_is_the_training_mean_table_and_lambda_is_ten_n():
    rng = np.random.default_rng(9)
    X, Y, _ = _linear_toy(rng, 6, p=40, g=10)
    m = cmap.fit_map(X, Y)
    np.testing.assert_allclose(m.predict(np.zeros(X.shape[1])), Y.mean(0), atol=1e-12)
    assert m.lam == 60.0
    with pytest.raises(ValueError):
        cmap.fit_map(X[:1], Y[:1])


def _stage0_inputs(rng, linear: bool, n_pert=18):
    p, g = 500, 120
    names = [f"k{i}" for i in range(n_pert)]
    rnd = {h: [f"r{h}_{j}" for j in range(cmap.N_RANDOM)] for h in names}
    allg = names + [x for v in rnd.values() for x in v]
    if linear:
        X, _, _ = _latent(rng, len(allg), p)
    else:
        X = rng.normal(0, 1, (len(allg), p))
        X = (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True)
    ctx = {k: X[i] for i, k in enumerate(allg)}
    W = rng.normal(0, 1, (p, 2)) @ rng.normal(0, 1, (2, g)) / p
    mean = rng.normal(0, 0.05, g)
    scale = float(np.std(X @ W))
    if linear:
        labels = {k: mean + ctx[k] @ W + 0.05 * scale * rng.normal(0, 1, g) for k in names}
    else:
        labels = {k: mean + rng.normal(0, 1, g) for k in names}
    return ctx, labels, rnd


def test_stage_zero_passes_on_the_linear_toy_and_fails_on_the_noise_toy():
    rng = np.random.default_rng(10)
    ctx, labels, rnd = _stage0_inputs(rng, linear=True)
    res = cmap.leave_one_out(ctx, labels, rnd, ctx)
    assert res["passes"] and res["n_held"] == 18
    assert res["vs_mean"]["share_ok"] and res["vs_rand"]["margin_ok"]
    ctx, labels, rnd = _stage0_inputs(rng, linear=False)
    res = cmap.leave_one_out(ctx, labels, rnd, ctx)
    assert not res["passes"]


def _per(margins_mean, margins_rand):
    return {f"h{i}": {"rho_map": 0.5, "rho_mean": 0.5 - a, "rho_rand": 0.5 - b}
            for i, (a, b) in enumerate(zip(margins_mean, margins_rand, strict=True))}


def test_the_gate_needs_two_thirds_and_a_median_margin_above_five_hundredths_against_both():
    ok = [0.2, 0.2, -0.1]  # 2 of 3 ahead, median 0.2
    assert cmap.gate(_per(ok, ok))["passes"]
    assert not cmap.gate(_per(ok, [0.2, -0.1, -0.1]))["passes"]  # 1 of 3 against the random queries
    assert not cmap.gate(_per([0.05, 0.05, 0.3], ok))["passes"]  # median exactly 0.05: not above
    assert cmap.gate(_per([0.0501, 0.06, 0.3], ok))["passes"]
    assert not cmap.gate(_per(ok[:2], ok[:2]))["passes"]  # fewer than three held out
    per = _per(ok, ok)
    per["h0"]["rho_rand"] = None  # an undefined rho is not a win for the map
    g = cmap.gate(per)
    assert g["vs_rand"]["wins"] == 1 and not g["passes"]


def test_spearman_is_none_for_a_constant_vector():
    assert cmap.spearman(np.ones(5), np.arange(5.0)) is None
    assert cmap.spearman(np.arange(5.0), np.arange(5.0)) == pytest.approx(1.0)


def test_a_context_orthogonal_to_every_training_context_gets_exactly_the_mean_table():
    rng = np.random.default_rng(11)
    X, Y, _ = _linear_toy(rng, 8, p=40, g=10)
    m = cmap.fit_map(X, Y)
    q, _ = np.linalg.qr(X.T, mode="complete")
    x = q[:, -1] * 6.0  # outside the training contexts' span
    np.testing.assert_allclose(m.predict(x), Y.mean(0), atol=1e-10)


def test_load_t1_pools_both_files_through_the_repositorys_reader(tmp_path):
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(12)
    genes = [f"x{i}" for i in range(25)]
    for f, n in zip(cmap.T1_FILES, (30, 20), strict=True):
        C = rng.poisson(3.0, (n, 25)).astype(float)
        C[:, 0] = 1.0  # every cell holds a count of 1, so its counts recover exactly
        X = np.log1p(1e4 * C / C.sum(1, keepdims=True))
        a = ad.AnnData(X=sparse.csr_matrix(X.astype(np.float32)),
                       obs=pd.DataFrame({"celltype": [f"T{i % 3}" for i in range(n)]},
                                        index=[f"{f}_{i}" for i in range(n)]),
                       var=pd.DataFrame(index=genes))
        p = tmp_path / f
        p.parent.mkdir(parents=True, exist_ok=True)
        a.write_h5ad(p)
    e = cmap.load_t1(tmp_path)
    assert e.X.shape == (50, 25) and e.genes == genes and list(np.unique(e.types)) == ["T0", "T1", "T2"]
    assert np.all(np.isfinite(e.ln_total)) and e.ln_total.min() > np.log(25)
    det = cmap.detection(e)
    assert det[0] == 1.0 and 0 < det[1] <= 1.0
