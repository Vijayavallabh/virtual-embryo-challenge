"""Heart try H6 `b0_cfm` (experiments/harness/TRIES.md), on synthetic data only: the covariance filter of a lineage's
momentum (vec.t2_build.covariance_filter), its hook in vec.t1_scope.momentum_shift, the candidate and the build flag."""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness, t1_scope
from vec.boards import Board
from vec.harness import BUILD_SEED_OFFSET, CANDIDATES, PseudoTask, load_tasks
from vec.t2_build import CFM_MIN_CELLS, CFM_Z_CLIP, covariance_filter

BUILDS = Path(__file__).resolve().parents[1] / "experiments" / "builds"


def _counts(n: int, p: int, seed: int, programme: float = 0.0, module: int = 40, depth: float = 60.0,
            depth_sigma: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """A lineage's counts at a sparse panel's depth: per-cell depth (lognormal), gene proportions; with `programme` a
    coordinated module: cells' scores (standard normal) move the first `module` genes' proportions by exp(a s), a of
    both signs. Returns (counts, a)."""
    rng = np.random.default_rng(seed)
    pi = rng.gamma(0.5, 1.0, p)
    pi /= pi.sum()
    depth_n = np.maximum(np.round(np.exp(rng.normal(np.log(depth), depth_sigma, n))), 5)
    a = np.zeros(p)
    a[:module] = programme * rng.choice([-1.0, 1.0], module) * rng.uniform(0.5, 1.5, module)
    lam = pi * np.exp(np.outer(rng.normal(size=n), a))
    lam = lam / lam.sum(1, keepdims=True) * depth_n[:, None]
    C = rng.poisson(lam).astype(np.int64)
    return C[C.sum(1) > 0], a


def _cos(x: np.ndarray, y: np.ndarray) -> float:
    return float(x @ y / np.sqrt((x @ x) * (y @ y)))


def test_a_planted_programme_is_kept_and_orthogonal_noise_removed():
    C, a = _counts(2500, 300, seed=1, programme=2.5)
    rng = np.random.default_rng(2)
    noise = rng.normal(0.0, 0.05, a.size)  # gene-independent noise over all genes, most of the momentum's energy
    logf = 0.6 * a / np.linalg.norm(a) + noise
    out, rec = covariance_filter(logf, C)
    assert rec["eligible"] and rec["status"] == "filtered" and rec["directions"] >= 1
    assert _cos(out, a) > _cos(logf, a) + 0.2, (_cos(out, a), _cos(logf, a))  # the programme's share of it rises
    assert _cos(out, a) > 0.75, _cos(out, a)
    module = np.abs(a) > 0
    # the energy outside the programme's genes falls: the gene-independent part is removed
    assert (out[~module] @ out[~module]) / (out @ out) < (logf[~module] @ logf[~module]) / (logf @ logf)
    assert 0.0 < rec["retained_energy"] < 1.0
    assert rec["spearman"] is not None and rec["max_abs_ratio"] > 0.0


def test_the_norm_rule_keeps_the_unfiltered_factors_l2_norm_exactly():
    """Fix 3 of the row: renormalised to the unfiltered lineage's L2 norm, so the candidate is not a partial hold."""
    C, a = _counts(2500, 300, seed=3, programme=2.5)
    logf = 0.4 * a / np.linalg.norm(a) + np.random.default_rng(4).normal(0.0, 0.04, a.size)
    out, rec = covariance_filter(logf, C)
    assert rec["status"] == "filtered"
    assert np.linalg.norm(out) == pytest.approx(np.linalg.norm(logf), rel=1e-12)


def test_pure_noise_keeps_b0s_factor():
    """Fix 4: no eigenvalue above the edge, B0's factor unchanged, whatever the factor is."""
    C, _ = _counts(1500, 300, seed=5)
    logf = np.random.default_rng(6).normal(0.0, 0.3, 300)
    out, rec = covariance_filter(logf, C)
    assert rec["eligible"] and rec["directions"] == 0 and "no eigenvalue" in rec["status"]
    assert rec["retained_energy"] is None and rec["spearman"] is None
    np.testing.assert_array_equal(out, logf)


def test_the_residualisation_removes_the_depth_axis():
    """Fix 2: strongly varying depth with no programme. The unresidualised standardised covariance has a depth axis far
    above the edge; the filter sees none and keeps B0's factor."""
    C, _ = _counts(2000, 300, seed=7, depth=60.0, depth_sigma=1.0)
    lib = C.sum(1)
    X = np.log1p(1e4 * C / lib[:, None])
    Z = X - X.mean(0)
    Z = Z / Z.std(0, ddof=1).clip(1e-12)
    ev = np.linalg.eigvalsh(Z.T @ Z / (len(Z) - 1))
    edge = (1.0 + np.sqrt(300 / (len(Z) - 2))) ** 2
    assert ev[-1] > 1.5 * edge  # the fixture has a depth axis
    logf = np.random.default_rng(8).normal(0.0, 0.3, 300)
    out, rec = covariance_filter(logf, C)
    assert rec["directions"] == 0, rec
    np.testing.assert_array_equal(out, logf)


def test_a_depth_linked_trend_in_the_gene_mix_is_removed_with_the_depth_axis():
    """Fix 2: the gene mix itself trends with depth (per-gene log-share slopes), no programme: an axis the depth
    regression has to remove, or the filter would keep the depth-linked trend. A line in log depth leaves it."""
    rng = np.random.default_rng(15)
    n, p = 4000, 300
    pi = rng.gamma(0.5, 1.0, p)
    pi /= pi.sum()
    depth = np.maximum(np.round(np.exp(rng.normal(np.log(60.0), 0.9, n))), 5)
    beta = rng.normal(0.0, 0.35, p)
    lam = pi * np.exp(np.outer(np.log(depth / 60.0), beta))
    C = rng.poisson(lam / lam.sum(1, keepdims=True) * depth[:, None]).astype(np.int64)
    C = C[C.sum(1) > 0]
    X = np.log1p(1e4 * C / C.sum(1, keepdims=True))
    d = np.log(C.sum(1))
    B = np.column_stack([np.ones(len(d)), d - d.mean()])  # a straight line: the depth model this test rules out
    Q = np.linalg.qr(B)[0]
    res = X - Q @ (Q.T @ X)
    Z = np.clip(res / np.sqrt((res ** 2).sum(0) / (len(d) - 2)).clip(1e-12), -CFM_Z_CLIP, CFM_Z_CLIP)
    ev = np.linalg.eigvalsh(Z.T @ Z / (len(d) - 2))
    assert ev[-1] > 1.3 * (1.0 + np.sqrt(p / (len(d) - 2))) ** 2  # the fixture has a depth-linked axis
    logf = np.random.default_rng(16).normal(0.0, 0.3, p)
    out, rec = covariance_filter(logf, C)
    assert rec["directions"] == 0, rec
    np.testing.assert_array_equal(out, logf)


def test_a_lineage_with_too_few_cells_keeps_b0s_factor():
    C, a = _counts(CFM_MIN_CELLS - 1, 300, seed=9, programme=2.5)
    logf = 0.5 * a / np.linalg.norm(a) + 0.01
    out, rec = covariance_filter(logf, C)
    assert not rec["eligible"] and rec["directions"] == 0 and str(CFM_MIN_CELLS) in rec["status"]
    np.testing.assert_array_equal(out, logf)
    # counted on cells with counts: empty cells do not reach the threshold
    C2 = np.vstack([_counts(CFM_MIN_CELLS - 1, 300, seed=9, programme=2.5)[0], np.zeros((5, 300), np.int64)])
    assert not covariance_filter(logf, C2)[1]["eligible"]


def test_a_gene_that_does_not_vary_keeps_its_factor():
    C, a = _counts(2500, 300, seed=10, programme=2.5)
    C[:, 299] = 0  # never expressed in this lineage: its factor can act on no cell
    logf = 0.6 * a / np.linalg.norm(a) + np.random.default_rng(11).normal(0.0, 0.05, 300)
    logf[299] = -1.7
    out, rec = covariance_filter(logf, C)
    assert rec["status"] == "filtered" and out[299] == -1.7
    assert np.linalg.norm(out) == pytest.approx(np.linalg.norm(logf), rel=1e-12)


def test_held_genes_are_left_out_of_the_filter_and_the_norm_rule():
    C, a = _counts(2500, 300, seed=17, programme=2.5)
    logf = 0.6 * a / np.linalg.norm(a) + np.random.default_rng(18).normal(0.0, 0.05, 300)
    held = np.zeros(300, dtype=bool)
    held[[3, 150, 299]] = True
    logf[held] = 0.0
    out, rec = covariance_filter(logf, C, held)
    assert rec["status"] == "filtered"
    assert (out[held] == 0.0).all()  # held genes keep their factor exactly
    assert np.linalg.norm(out) == pytest.approx(np.linalg.norm(logf), rel=1e-12)  # the norm rule over the rest


def test_the_filter_is_deterministic_and_leaves_its_inputs_alone():
    C, a = _counts(1200, 200, seed=12, programme=2.5)
    logf = 0.6 * a / np.linalg.norm(a) + np.random.default_rng(13).normal(0.0, 0.05, 200)
    logf0, C0 = logf.copy(), C.copy()
    o1, r1 = covariance_filter(logf, sparse.csr_matrix(C))
    o2, r2 = covariance_filter(logf, C)
    np.testing.assert_array_equal(o1, o2)
    assert r1 == r2
    np.testing.assert_array_equal(logf, logf0)
    np.testing.assert_array_equal(C, C0)


def test_the_record_holds_no_count():
    C, a = _counts(1800, 200, seed=14, programme=2.5)
    _, rec = covariance_filter(0.6 * a / np.linalg.norm(a) + 0.01, C)
    assert rec["status"] == "filtered"
    assert set(rec) == {"eligible", "status", "directions", "retained_energy", "spearman", "max_abs_ratio"}


def test_factors_with_a_wrong_length_are_refused():
    with pytest.raises(ValueError, match="genes"):
        covariance_filter(np.zeros(4), np.ones((10, 5)))


# ---------------------------------------------------------------------------------------------------- momentum_shift


def _adata(C: np.ndarray, types: list[str], genes: list[str], prefix: str) -> ad.AnnData:
    a = ad.AnnData(X=np.log1p(1e4 * C / C.sum(1, keepdims=True).clip(1)).astype(np.float32),
                   obs=pd.DataFrame({"b0_type": types}, index=[f"{prefix}{i}" for i in range(len(C))]),
                   var=pd.DataFrame(index=genes))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    return a


@pytest.fixture(scope="module")
def shift_inputs():
    rng = np.random.default_rng(21)
    genes = [f"g{i}" for i in range(50)]
    mu = rng.gamma(2.0, 2.0, 50)
    start = _adata(rng.poisson(mu, (300, 50)), ["A"] * 200 + ["B"] * 100, genes, "s")
    ref = _adata(rng.poisson(mu * 1.5, (400, 50)), ["A"] * 250 + ["B"] * 120 + ["C"] * 30, genes, "r")
    return start, ref, ref[np.arange(0, 400, 3)].copy(), genes


def test_momentum_shift_without_a_filter_is_unchanged(shift_inputs):
    """The default None is the code path as before; an identity hook reproduces it, so the filter sits where it should
    and adds no draw of its own."""
    start, ref, like, _ = shift_inputs
    base, rep = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3])
    calls = []

    def identity(logf, counts, held):
        calls.append((counts.shape[0], float(logf.sum()), int(held.sum())))
        return logf, {"status": "same"}

    out, rep2 = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], logf_filter=identity)
    np.testing.assert_array_equal(base.layers["counts"].toarray(), out.layers["counts"].toarray())
    # the hook sees each shifted type's reference cells: type C has no momentum (no start cells) and no filter call
    assert sorted(c[0] for c in calls) == [120, 250] and all(c[2] == 0 for c in calls)
    by_type = {r["type"]: r for r in rep2.rows}
    assert by_type["A"]["filter"] == {"status": "same"} and "filter" not in by_type["C"]
    assert all("filter" not in r for r in rep.rows)


def test_a_held_gene_stays_flat_whatever_the_filter_returns(shift_inputs):
    start, ref, like, genes = shift_inputs
    seen = []

    def spread(logf, counts, held):  # a non-uniform factor: a uniform one is a no-op under keep_depth
        seen.append(held.copy())
        return np.linspace(-0.5, 0.5, logf.size), {"status": "spread"}

    out, _ = t1_scope.momentum_shift(start, ref, like, 0.5, held_flat=[genes[0]], key="b0_type", seed=[7, 3],
                                     logf_filter=spread)
    assert seen and all(h[0] and h.sum() == 1 for h in seen)
    old = sparse.csr_matrix(like.layers["counts"])
    new = sparse.csr_matrix(out.layers["counts"])
    assert (old[:, 0] != new[:, 0]).count_nonzero() == 0
    assert (old[:, 1:] != new[:, 1:]).count_nonzero() > 0


# ---------------------------------------------------------------------------------------------------- the candidate
GENES = [f"g{i}" for i in range(150)]


def _pool(n: int, seed: int, shift: float = 0.0) -> ad.AnnData:
    """A T2 stage: a CM lineage with a coordinated programme (strong enough to pass the edge), Endothelial cells, and
    Untracked cells (held flat by B0)."""
    rng = np.random.default_rng(seed)
    types = rng.choice(3, size=n, p=[0.5, 0.35, 0.15])
    pi = np.random.default_rng(1).gamma(1.0, 1.0, 150)
    pi /= pi.sum()
    a = np.zeros(150)
    a[:30] = 2.5 * rng.choice([-1.0, 1.0], 30) * rng.uniform(0.5, 1.5, 30)
    s = rng.normal(size=n) * (types == 0)
    lam = pi * np.exp(np.outer(s, a) + shift * (types == 0)[:, None])
    depth = np.maximum(np.round(np.exp(rng.normal(np.log(80.0), 0.4, n))), 5)[:, None]
    C = rng.poisson(lam / lam.sum(1, keepdims=True) * depth).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    X = np.log1p(1e4 * C / C.sum(1, keepdims=True)).astype(np.float32)
    centres = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 5.0, 0.0]])
    a2 = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame(
        {"celltype": [f"type{t}" for t in types], "merged": np.array(["CM", "Endothelial", "Untracked"])[types]},
        index=[f"c{seed}_{i}" for i in range(n)]), var=pd.DataFrame(index=GENES))
    a2.layers["counts"] = sparse.csr_matrix(C)
    a2.obsm["spatial_3D"] = (centres[types] + rng.normal(size=(n, 3)) * 1.5).astype(np.float32)
    return a2


def _heart_task(**over) -> PseudoTask:
    kw = {"key": "syn_h6", "test_board": "T2:heart:test_extrap", "task": "T2", "kind": "extrap",
          "target": "T2/E10.5.h5ad", "reference": "T2/E9.5.h5ad", "setting": "heart",
          "inputs": {"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"}, "validation_board": "T2:heart:val_extrap"}
    return PseudoTask(**{**kw, **over})


@pytest.fixture(scope="module")
def ctx():
    board = Board(key="pseudo:h6", task="T2", genes=list(GENES), needs_coords=True, setting="heart", min_cells=100,
                  max_cells=800)
    inputs = {"ref": _pool(1800, seed=2), "start": _pool(1200, seed=3, shift=0.6)}
    return harness.BuildContext(task=harness.TaskView.of(_heart_task()), board=board, inputs=inputs,
                                seed=BUILD_SEED_OFFSET + 1)


def test_b0_cfm_applies_on_the_heart_backtests_only():
    c = CANDIDATES["b0_cfm"]
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    assert c.applies(tasks["T2h/<=E9.5->E10.5"]) and c.applies(tasks["T2h/<=E8.75->E10.5"])
    assert not c.applies(tasks["T2e/E7.25+E8.0->E7.5"]) and not c.applies(tasks["T1/E9.5->E10.5"])
    assert not c.applies(tasks["T3/Gata4->Mab21l2"])
    assert not c.applies(_heart_task(inputs={"ref": "T2/E9.5.h5ad"}))  # no momentum start: B0's own rule


def test_b0_cfm_changes_only_the_moved_lineages_expression(ctx):
    """B0's draw and positions; Untracked cells stay flat; the expression of the filtered lineage changes."""
    b0, cfm = CANDIDATES["b0"].build(ctx), CANDIDATES["b0_cfm"].build(ctx)
    np.testing.assert_array_equal(b0.coords, cfm.coords)
    assert cfm.n == b0.n == ctx.board.max_cells
    pool, rows, _, _ = harness.draw_copy(ctx)
    merged = np.asarray(pool.obs["merged"]).astype(str)[rows]
    assert not np.array_equal(b0.X, cfm.X)
    untracked = merged == "Untracked"
    assert untracked.any() and np.array_equal(b0.X[untracked], cfm.X[untracked])  # held flat in both
    moved = merged == "CM"
    assert moved.any() and not np.array_equal(b0.X[moved], cfm.X[moved])


def test_cfm_goes_through_the_filter_and_only_on_extrapolation(ctx, monkeypatch):
    seen = []
    real = t1_scope.momentum_shift

    def spy(*a, **k):
        seen.append(k.get("logf_filter"))
        return real(*a, **k)

    monkeypatch.setattr(t1_scope, "momentum_shift", spy)
    harness._b0_t2(ctx)
    harness._b0_t2(ctx, cfm=True)
    from vec.t2_build import covariance_filter as cf
    assert seen == [None, cf]  # B0 passes none: unchanged
    interp = replace(ctx, task=replace(ctx.task, kind="interp"))
    with pytest.raises(ValueError, match="extrapolation"):
        harness._b0_t2(interp, cfm=True)


def test_build_test_cfm_flag_is_off_by_default():
    sys.path.insert(0, str(BUILDS))
    import build_test as bt

    ap = bt.build_parser()
    assert ap.parse_args(["heart", "--entry", "s1", "--out", "runs/x/h.h5ad"]).cfm is False
    assert ap.parse_args(["heart", "--entry", "s1", "--cfm", "--out", "runs/x/h.h5ad"]).cfm is True
