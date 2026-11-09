"""T12a and T12b (experiments/harness/TRIES.md), on synthetic data only: the strong-null down-tail prior's bins, gene
order and realisation (`vec.t1_scope`), its backtest candidates (`vec.harness`), the build option
(`experiments/builds/build_test.py`) and the stop-rule script (`experiments/t1_strong_null`)."""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import replace
from functools import partial
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_hidden_counts import FENG_GROUPS, T1_GENES, _reference_file, bt, t1_env  # noqa: F401 (the build's fixtures)

from vec import harness
from vec import t1_scope as ts
from vec.boards import Board
from vec.harness import BuildContext, PseudoTask, TaskView
from vec.prep import log1p_cp10k

EXPERIMENT = Path(__file__).resolve().parents[1] / "experiments" / "t1_strong_null"
sys.path.insert(0, str(EXPERIMENT))
import stop_rule  # noqa: E402
import strong_null  # noqa: E402

# ---------------------------------------------------------------------------------------------------------------------
# The rule's numbers (the rows T12a and T12b)


def test_the_constants_are_the_rows():
    assert ts.DP_LEVEL_EDGES == (100, 200, 400, 800, 1600, 3200) and ts.DP_PROMOTED_LEVELS == 4
    assert ts.DP_CHANGE_EDGES == (400, 1600) and ts.DP_CHANGE_TOP == 1600 and ts.DP_PROMOTED_CHANGE == (1, 2)
    assert ts.DP_BUDGET == 200 and ts.DP_BAND == -0.15
    assert ts.DP_F_MIN == 0.5 and ts.DP_MIN_PB == 0.25 and ts.DP_ITERATIONS == 3 and ts.DP_TOL == 0.02
    assert ts.DP_DEPTH_TOL == 0.01 and ts.DP_BOTTOM == 400
    assert ts.DP_MMD_BOUND == 0.0035 and ts.DP_MMD_CELLS == 5000
    assert stop_rule.LOSS_BOUND == -0.66 and stop_rule.APPROVED_LOSS == -0.33 and stop_rule.DIST_BOUND == -0.33
    assert stop_rule.RAW_MMD_BOUND == ts.DP_MMD_BOUND and stop_rule.RAW_VARIOGRAM_BOUND == 0.0004
    assert stop_rule.DF_MMD == (0.024, 0.05, 0.08, 0.11) and stop_rule.DF_VARIOGRAM == (0.0009, 0.003, 0.0055, 0.01)
    assert stop_rule.PAIRS == {"T1/E9.5->E10.5": ("b0", "b0_strong_null"),
                               "T1/E9.5->E10.5:chamber": ("b0", "b0_strong_null"),
                               "T1/E8.5->E9.5": ("b0_scope", "scope_strong_null")}
    assert harness.T12_GROUPS == ("CM", "Endocardium", "Endothelium", "Epicardium")


# ---------------------------------------------------------------------------------------------------------------------
# Bins and gene order


def test_bins_are_ranks_against_the_edges():
    n = 4000
    pb = np.linspace(9.0, 0.0, n)  # gene 0 is the highest
    dp = np.arange(n, dtype=float) - 2000.0  # gene 0 the most negative
    level, change = ts.down_prior_bins(pb, dp, ts.DownPriorTable(p=np.zeros((7, 4))))
    assert level[[0, 99, 100, 199, 200, 3199, 3200, n - 1]].tolist() == [0, 0, 1, 1, 2, 5, 6, 6]
    assert change[[0, 399, 400, 1599, 1600, 2399, 2400, n - 1]].tolist() == [0, 0, 1, 1, 2, 2, 3, 3]
    # ties in the level take the gene's position (a stable sort)
    level2, _ = ts.down_prior_bins(np.zeros(4), np.zeros(4), ts.DownPriorTable(p=np.zeros((3, 4)), level_edges=(1, 2),
                                                                               change_edges=(1,), change_top=1))
    assert level2.tolist() == [0, 1, 2, 2]


def test_level_only_table_falls_with_the_level_and_ignores_the_change():
    t = ts.DownPriorTable.level_only()
    assert t.p.shape == (7, 4) and (t.p[-1] == 0).all() and (np.diff(t.p, axis=0) < 0).all()
    assert all(len(set(row)) == 1 for row in t.p)


def test_order_takes_only_eligible_genes_best_first():
    table = ts.DownPriorTable(p=np.array([[0.0, 0.6, 0.9, 0.0], [0.0, 0.5, 0.5, 0.0], [0.0, 0.3, 0.3, 0.0]]),
                              level_edges=(3, 6), change_edges=(1, 2), change_top=1)
    genes = list("abcdefghij")
    level = np.array([0, 0, 1, 1, 1, 2, 0, 0, 1, 1])
    change = np.array([2, 1, 2, 1, 2, 2, 0, 3, 2, 2])
    #                  a  b  c  d  e  f  g  h  i  j   (f: the last level bin; g: the tail; h: predicted up)
    pb = np.ones(10)
    pb[9] = 0.2  # j is below MIN_LFC: it can never be a truth down-gene
    held = np.zeros(10, bool)
    held[0] = False
    held[8] = True  # i is held
    dp = np.array([0.0, 0.5, -0.2, 0.3, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0])
    order = ts.down_prior_order(table, level, change, pb, dp, genes, held)
    assert [genes[i] for i in order] == ["a", "b", "c", "e", "d"]  # 0.9 a; 0.6 b; 0.5 c, e, d by dp


def test_order_never_promotes_past_the_rank_800_cap_whatever_the_table_says():
    n = 4000
    pb = np.linspace(9.0, 0.5, n)  # gene 0 the highest level; every gene clears MIN_LFC
    dp = np.linspace(-1.0, 1.0, n)[np.random.default_rng(0).permutation(n)]  # the change bins independent of the level
    # a table that favours the low levels: the cap must still stop anything from level bin 4 on (level ranks 800+)
    p = np.tile(np.arange(7, dtype=np.float64)[:, None], (1, 4))
    table = ts.DownPriorTable(p=p)
    level, change = ts.down_prior_bins(pb, dp, table)
    order = ts.down_prior_order(table, level, change, pb, dp, [str(i) for i in range(n)], np.zeros(n, bool))
    assert len(order) and level[order].max() == 3 and (np.isin(change[order], ts.DP_PROMOTED_CHANGE)).all()
    assert order.max() < 800  # the genes' positions are their level ranks here, the highest level first


# ---------------------------------------------------------------------------------------------------------------------
# The realisation


def _population(n_cells=1500, seed=0):
    """A synthetic population: 100 genes in panel order, genes 0-29 high, 30-59 low, 60-99 off; half the cells
    covered."""
    rng = np.random.default_rng(seed)
    mu = np.concatenate([rng.uniform(40, 200, 30), rng.uniform(0.5, 4, 30), np.zeros(40)])
    C = rng.poisson(mu, (n_cells, 100)).astype(np.int32)
    C[:, 90] = 3  # a gene everywhere the same (a stand-in for a held gene)
    genes = [f"g{i}" for i in range(100)]
    a = ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), layers={"counts": sparse.csr_matrix(C)},
                   obs=pd.DataFrame({"celltype": ["V-CM" if i % 2 == 0 else "Blood" for i in range(n_cells)]},
                                    index=[f"c{i}" for i in range(n_cells)]))
    a.var_names = genes
    return a, genes


def _prior(table=None, **kw):
    table = table or ts.DownPriorTable(p=np.array([[0.0, 0.9, 0.8, 0.0], [0.0, 0.5, 0.4, 0.0], [0.0, 0.1, 0.1, 0.0]]),
                                       level_edges=(20, 30), change_edges=(10, 20), change_top=10)
    kw = {"budget": 5, "band": -0.2, "f_min": 0.3, "tol": 0.05, "mmd_bound": 0.5, **kw}  # the 100-gene world's
    # promoted genes are a twentieth of the panel, far more than the real panel's: its file-distance bound differs
    return ts.DownPrior(table=table, **kw)


def _run(pop, genes, covered=None, held=("g90",), prior=None, seed=11):
    pb_pop = np.asarray(pop.X.mean(axis=0)).ravel().astype(np.float64)
    offset = np.zeros(100)
    offset[30:40] = 0.3  # a block of genes the build calls down: they make the bottom of dp
    pb_ref = pb_pop + offset  # dp = -offset for the panel, near 0 elsewhere
    covered = np.array(pop.obs["celltype"] == "V-CM") if covered is None else covered
    return ts.apply_down_prior(pop, pb_ref, genes, covered, held, prior or _prior(), seed), pb_ref, covered


def test_realisation_moves_the_promoted_genes_to_the_band_and_nothing_else():
    pop, genes = _population()
    (out, rep), pb_ref, covered = _run(pop, genes)
    C0, C1 = sparse.csr_matrix(pop.layers["counts"]), sparse.csr_matrix(out.layers["counts"])
    # who was promoted: high genes with dp above the band, by P (level 0 near-tail/flat cells first, ties by the most
    # negative dp), never the held gene, never the block already at the bottom
    assert rep["n_promoted"] == 5 and set(rep["promoted"]) <= {f"g{i}" for i in range(30)}
    assert "g90" not in rep["promoted"]
    idx = [genes.index(g) for g in rep["promoted"]]
    dp_before = np.asarray(pop.X.mean(axis=0)).ravel() - pb_ref
    dp_after = np.asarray(out.X.mean(axis=0)).ravel() - pb_ref
    assert (dp_after[idx] <= rep["band_target_dp"] + 0.05).all() and (dp_after[idx] < dp_before[idx]).all()
    # the band is the rule's constant (here the test's -0.2), the promoted genes' dp after the edit
    assert rep["band_target_dp"] == -0.2
    # uncovered cells keep every count; the held gene is untouched everywhere
    assert (C1[~covered] != C0[~covered]).nnz == 0
    assert (C1[:, 90] != C0[:, 90]).nnz == 0
    # only the promoted genes fall in the covered cells' means; every other gene moves by the depth constant only
    other = np.ones(100, bool)
    other[idx] = False
    mean_cov0, mean_cov1 = np.asarray(C0[covered].mean(0)).ravel(), np.asarray(C1[covered].mean(0)).ravel()
    assert ((mean_cov1 - mean_cov0)[idx] < 0).all()
    assert (mean_cov1[other & (mean_cov0 > 1.0)] >= mean_cov0[other & (mean_cov0 > 1.0)] * 0.999).all()
    # each covered cell's expected library is kept; the realised libraries stay close to it on average
    lib0, lib1 = np.asarray(C0[covered].sum(1)).ravel(), np.asarray(C1[covered].sum(1)).ravel()
    assert abs(np.mean(lib1 / lib0) - 1.0) <= 0.01
    checks = rep["checks"]
    assert checks["passed"] and checks["band_ok"] and checks["expected_depth_ok"] and checks["n_promoted_ok"]
    assert checks["expected_depth_max_error"] < 1e-9
    # the file-distance check runs and is reported (an unbiased MMD can be slightly negative at a true distance of 0)
    assert checks["mmd_pre_post_ok"] and isinstance(checks["mmd_pre_post"], float)
    assert checks["pc_mean_shift_covered"] > 0.0
    assert rep["keep_depth"] is True and rep["budget"] == 5 and rep["f_min"] == 0.3
    assert 0 < min(rep["factors"].values()) and max(rep["factors"].values()) < 1  # one thinning factor per gene


def test_realisation_is_deterministic_and_records_no_cell_count():
    pop, genes = _population(seed=2)
    (a, rep_a), _, _ = _run(pop, genes)
    (b, rep_b), _, _ = _run(pop, genes)
    assert (sparse.csr_matrix(a.X) != sparse.csr_matrix(b.X)).nnz == 0 and rep_a == rep_b
    assert str(pop.n_obs) not in json.dumps({k: v for k, v in rep_a.items() if k != "factors"})


def test_a_gene_out_of_reach_is_skipped_for_the_next():
    pop, genes = _population()
    # a gene expressed only in the uncovered cells: thinning the covered cells cannot move it
    C = sparse.csr_matrix(pop.layers["counts"]).toarray()
    covered = np.array(pop.obs["celltype"] == "V-CM")
    C[covered, 0] = 0
    C[~covered, 0] = C[~covered, 0] + 200
    pop.layers["counts"], pop.X = sparse.csr_matrix(C), log1p_cp10k(sparse.csr_matrix(C))
    (_, rep), _, _ = _run(pop, genes)
    assert "g0" not in rep["promoted"] and rep["n_promoted"] == 5


def test_no_covered_cells_is_refused():
    pop, genes = _population()
    pb_ref = np.asarray(pop.X.mean(axis=0)).ravel()
    with pytest.raises(ValueError, match="covered cells are empty"):
        ts.apply_down_prior(pop, pb_ref, genes, np.zeros(pop.n_obs, bool), (), _prior(), 1)


def test_a_shortfall_or_an_empty_promotion_fails_the_checks():
    pop, genes = _population()
    (_, rep), _, _ = _run(pop, genes, prior=_prior(budget=31))  # the 30 high genes are all that can reach the band
    assert rep["n_promoted"] < 31 and not rep["checks"]["n_promoted_ok"] and not rep["checks"]["passed"]
    # nothing eligible at all (every gene held): the file would be M1's, but the checks must refuse it, not pass it
    (_, rep), _, _ = _run(pop, genes, held=tuple(genes), prior=_prior(budget=5))
    assert rep["n_promoted"] == 0 and not rep["checks"]["n_promoted_ok"] and not rep["checks"]["passed"]


def test_the_file_distance_above_its_bound_fails_the_checks():
    pop, genes = _population()
    (_, ok), _, _ = _run(pop, genes)
    (_, refused), _, _ = _run(pop, genes, prior=_prior(mmd_bound=ok["checks"]["mmd_pre_post"] - 1e-4))
    assert ok["checks"]["passed"] and not refused["checks"]["mmd_pre_post_ok"] and not refused["checks"]["passed"]
    # the numbers are the same whether or not the bound refuses them (reported either way)
    assert refused["checks"]["mmd_pre_post"] == ok["checks"]["mmd_pre_post"]


def test_the_prepost_distance_is_the_scorers_mmd_between_the_two_files():
    from common.core_metrics import mmd_unbiased

    pop, genes = _population()
    (out, rep), _, _ = _run(pop, genes)
    a = pop[:, genes].X
    b = out[:, genes].X
    a, b = a.toarray() if sparse.issparse(a) else a, b.toarray() if sparse.issparse(b) else b
    assert rep["checks"]["mmd_pre_post"] == pytest.approx(
        mmd_unbiased(b.astype(np.float32), a.astype(np.float32), n_pc=30), abs=1e-5)


def test_the_distance_check_reads_nothing_refused_under_the_harness_guard():
    pop, genes = _population()
    with harness.guarded_reads() as g:
        _run(pop, genes)
    assert g["refused"] == []


def test_sparse_genes_inside_the_cap_meet_the_band_at_t1_depth():
    # the row caps promotion at level rank 800 because the first-order solve undershoots the band for sparser genes;
    # at T1's depth (about 12K UMI per cell) the genes the cap admits meet it
    rng = np.random.default_rng(3)
    n_genes, n_cells = 3400, 1200
    weights = np.arange(1, n_genes + 1, dtype=np.float64) ** -1.1
    mu = 12000.0 * weights / weights.sum()
    C = sparse.csr_matrix(rng.poisson(mu, (n_cells, n_genes)).astype(np.int32))
    genes = [f"g{i}" for i in range(n_genes)]
    pop = ad.AnnData(X=log1p_cp10k(C), layers={"counts": C},
                     obs=pd.DataFrame({"celltype": ["V-CM"] * n_cells}, index=[f"c{i}" for i in range(n_cells)]))
    pop.var_names = genes
    pb_pop = np.asarray(pop.X.mean(axis=0)).ravel().astype(np.float64)
    pb_ref = pb_pop + rng.normal(0.0, 1e-4, n_genes)  # dp is a small noise, so the change bins do not follow the level
    prior = ts.DownPrior(table=ts.DownPriorTable.level_only(), mmd_bound=1.0)
    out, rep = ts.apply_down_prior(pop, pb_ref, genes, np.ones(n_cells, bool), (), prior, seed=[5, 4])
    idx = np.array([genes.index(g) for g in rep["promoted"]])
    level, _ = ts.down_prior_bins(pb_ref, pb_pop - pb_ref, prior.table)
    assert rep["n_promoted"] == ts.DP_BUDGET and level[idx].max() < ts.DP_PROMOTED_LEVELS
    assert rep["checks"]["band_ok"] and rep["checks"]["band_share"] >= 0.95 and rep["checks"]["n_promoted_ok"]


# ---------------------------------------------------------------------------------------------------------------------
# build_test.py t1


def test_strong_null_option_is_off_by_default_and_refused_until_approved(bt):  # noqa: F811
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "runs/x/T1_s1.h5ad"]).strong_null is False
    assert ap.parse_args(["t1", "--strong-null", "--out", "runs/x/T1_s1.h5ad"]).strong_null is True
    assert bt._t1_down_prior(False, True, True) is None  # off: nothing is read
    with pytest.raises(SystemExit, match="without R1 and R2"):
        bt._t1_down_prior(True, False, True)
    with pytest.raises(SystemExit, match="without R1 and R2"):
        bt._t1_down_prior(True, True, False)
    with pytest.raises(SystemExit, match="does not say approved"):  # the repository's config is not approved yet
        bt._t1_down_prior(True, True, True)


def test_the_repository_config_is_not_approved_yet():
    assert strong_null.load_config()["approved"] is False


def test_the_approval_switch_refuses_anything_but_true():
    for cfg in ({"approved": False}, {"approved": "true"}, {}, {"approved": 1}):
        with pytest.raises(strong_null.Refused, match="does not say approved"):
            strong_null.check_approved(cfg)
    assert strong_null.check_approved({"approved": True}) is None


def test_a_strong_null_build_refuses_an_existing_file(bt, tmp_path):  # noqa: F811
    out = tmp_path / "T1_s1.h5ad"
    out.write_text("x")
    with pytest.raises(SystemExit, match="never overwrites"):
        bt.main(["t1", "--strong-null", "--out", str(out)])


def test_option_off_writes_todays_file(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811 (the fixtures imported above)
    # the default build (R1 and R2 on, the union pinned), without --strong-null, is the file build_t1 writes today
    union = {"file": "f", "sha256": "ab" * 32, "genes": ["g5", "g7", "gone_x", "held_a"], "modules": {"D": 2},
             "sources": {}}
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: union)
    today = _reference_file(bt, t1_env, tmp_path / "today.h5ad",
                            held=tuple(sorted(set(T1_HELD_NAMES) | set(union["genes"]))))
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", "--out", str(out)])
    assert out.read_bytes() == today.read_bytes()
    assert "strong_null" not in json.loads(out.with_suffix(".build.json").read_text())


T1_HELD_NAMES = ["held_a", "held_b"]  # test_hidden_counts.T1_HELD, the synthetic pre-registered list


def _approved_world(bt, monkeypatch, **prior_kw):  # noqa: F811
    """The synthetic world's build with the prior approved: a small order on its own scale (62 genes) and the rule's
    constants made relative to it."""
    table = ts.DownPriorTable(p=np.array([[0.0, 0.9, 0.8, 0.0], [0.0, 0.7, 0.6, 0.0], [0.0, 0.1, 0.1, 0.0]]),
                              level_edges=(20, 40), change_edges=(12, 30), change_top=8)
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: {"file": "f", "sha256": "ab" * 32, "genes": ["g5"],
                                                                   "modules": {"D": 1}, "sources": {}})
    monkeypatch.setattr(bt.strong_null, "check_approved", lambda *a, **k: None)
    monkeypatch.setattr(bt.t1_scope.DownPriorTable, "level_only", classmethod(lambda cls: table))
    kw = {"budget": 4, "band": -0.3, "f_min": 0.2, "tol": 0.2, "mmd_bound": 5.0, **prior_kw}
    monkeypatch.setattr(bt.t1_scope, "DownPrior", partial(ts.DownPrior, **kw))


def test_option_on_moves_only_the_covered_cells(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    _approved_world(bt, monkeypatch)
    off, on = tmp_path / "off.h5ad", tmp_path / "on.h5ad"
    bt.main(["t1", "--out", str(off)])
    t1_env.reads.clear()
    bt.main(["t1", "--strong-null", "--out", str(on)])
    assert t1_env.reads == [("Feng", [""]), "E10.5", ("Feng", T1_GENES)]  # nothing read before the prior's refusals
    A, B = ad.read_h5ad(off).X, ad.read_h5ad(on).X
    A, B = A.toarray() if sparse.issparse(A) else A, B.toarray() if sparse.issparse(B) else B
    assert A.shape == B.shape
    rec = json.loads(on.with_suffix(".build.json").read_text())
    strong = rec["strong_null"]
    assert "level" in strong["selection"] and "table" not in strong
    checks = strong["checks"]
    assert checks["passed"] is True and checks["n_promoted_ok"] and checks["mmd_pre_post_ok"]
    assert isinstance(checks["mmd_pre_post"], float) and strong["constants"]["mmd_bound"] == 5.0
    assert strong["n_promoted"] == 4 and rec["n_cells"] == "max_cells"
    assert str(t1_env.board.max_cells) not in on.with_suffix(".build.json").read_text()
    # the same cells drawn and ordered: only rows that hold covered cells differ, and only some hold promoted genes
    changed = (A != B).any(axis=1)
    assert changed.any() and not changed.all()  # R1's cells (blood, macrophages, Early VIC) are the unchanged rows
    assert (A[~changed] == B[~changed]).all()
    # the promoted genes are thinner in the on file
    idx = [T1_GENES.index(g) for g in strong["promoted"]]
    assert (B[:, idx].mean(0) < A[:, idx].mean(0)).all()


def test_option_on_writes_no_file_when_a_check_refuses(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    _approved_world(bt, monkeypatch, mmd_bound=-1.0)  # a bound no estimate can meet: the file-distance check refuses
    out = tmp_path / "on.h5ad"
    with pytest.raises(SystemExit, match="target-free checks failed"):
        bt.main(["t1", "--strong-null", "--out", str(out)])
    assert not out.exists() and not out.with_suffix(".build.json").exists()


def test_option_on_refused_before_reading_without_r1(bt, t1_env, tmp_path):  # noqa: F811
    with pytest.raises(SystemExit, match="without R1 and R2"):
        bt.main(["t1", "--strong-null", "--no-hold-unbracketed", "--out", str(tmp_path / "x.h5ad")])
    assert t1_env.reads == []


# ---------------------------------------------------------------------------------------------------------------------
# The harness candidates


def _labelled(n, labels, seed, genes):
    rng = np.random.default_rng(seed)
    mu = rng.gamma(2.0, 8.0, len(genes))
    C = rng.poisson(mu, (n, len(genes))).astype(np.int32)
    C[np.arange(n), rng.integers(0, len(genes), n)] += 1
    a = ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), layers={"counts": sparse.csr_matrix(C)},
                   obs=pd.DataFrame({"celltype": np.resize(labels, n)}, index=[f"c{i}" for i in range(n)]))
    a.var_names = genes
    return a


GENES = [f"g{i}" for i in range(50)]
LABELS = ["V-CM", "Endocardium", "Blood", "Macrophage"]  # the last two are not among T12_GROUPS: R1's cells at test


def _ctx(task, inputs, derived=None, seed=harness.BUILD_SEED_OFFSET + 1):
    board = Board(key="pseudo", task="T1", genes=list(GENES), needs_coords=False, min_cells=50, max_cells=700)
    return BuildContext(task=TaskView.of(task), board=harness._candidate_board(board), inputs=inputs, seed=seed,
                        derived=derived or {})


def _small_prior(table):
    def prior(ctx):
        return ts.DownPrior(table=table, covered_labels=harness._t12_labels(ctx), budget=4, band=-0.3, f_min=0.2,
                            tol=0.2)
    return prior


def test_candidates_apply_where_their_inputs_allow():
    tasks = {t.key: t for t in harness.load_tasks()["T1:test"]}
    b0 = harness.CANDIDATES["b0_strong_null"]
    control = harness.CANDIDATES["scope_strong_null"]
    assert b0.applies(tasks["T1/E9.5->E10.5"]) and b0.applies(tasks["T1/E9.5->E10.5:chamber"])
    assert not b0.applies(tasks["T1/E8.5->E9.5"]) and control.applies(tasks["T1/E8.5->E9.5"])
    assert not control.applies(tasks["T1/E9.5->E10.5"])
    assert b0.params["scope"] == "chamber" == control.params["scope"]
    assert b0.derived == ("t1_held_flat",) and control.derived == ()
    assert "t1_down_prior_weak" not in harness.DERIVED


def test_both_candidates_order_genes_by_the_level_alone(monkeypatch):
    seen = []
    real = ts.DownPrior
    monkeypatch.setattr(ts, "DownPrior", lambda **kw: seen.append(kw["table"]) or real(**kw))
    monkeypatch.setattr(harness, "_t12_labels", lambda ctx: frozenset({"V-CM"}))
    prior = harness._t12_prior(None)
    assert prior.covered_labels == frozenset({"V-CM"})
    assert seen[0].p.tolist() == ts.DownPriorTable.level_only().p.tolist()


def test_candidate_shares_b0s_draw_and_moves_only_its_cells(monkeypatch):
    task = PseudoTask(key="T1/E9.5->E10.5", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      validation_board="T1:val")
    ref, start = _labelled(1500, LABELS, 1, GENES), _labelled(1200, LABELS, 2, GENES)
    table = ts.DownPriorTable(p=np.array([[0.0, 0.9, 0.8, 0.0], [0.0, 0.7, 0.6, 0.0], [0.0, 0.1, 0.1, 0.0]]),
                              level_edges=(15, 30), change_edges=(10, 25), change_top=5)
    derived = {"t1_held_flat": ("g49",)}
    base = harness.CANDIDATES["b0"].build(_ctx(task, {"ref": ref, "start": start}, derived))
    monkeypatch.setattr(harness, "_t12_prior", _small_prior(table))
    changed = harness.CANDIDATES["b0_strong_null"].build(_ctx(task, {"ref": ref, "start": start}, derived))
    assert base.X.shape == changed.X.shape
    row_changed = (base.X != changed.X).any(axis=1)
    assert row_changed.any() and not row_changed.all()  # the same cells; only the prior's cells differ
    # the chosen genes are down on average in the changed file, and the held gene never moves
    delta = changed.X.mean(0) - base.X.mean(0)
    assert delta.min() < -0.02 and np.isclose(delta[GENES.index("g49")], 0.0, atol=1e-3)


def test_control_candidate_needs_no_table_and_touches_only_its_cells(monkeypatch):
    task = PseudoTask(key="T1/E8.5->E9.5", test_board="T1:test", task="T1", kind="extrap", target="T1/E9.5_RNA.h5ad",
                      reference="T1/E8.5_RNA.h5ad", inputs={"ref": "T1/E8.5_RNA.h5ad"}, twin_board="T1:val")
    ref = _labelled(1500, LABELS, 3, GENES)
    table = ts.DownPriorTable.level_only(level_edges=(15, 30), change_edges=(10, 25), change_top=5)
    monkeypatch.setattr(harness, "_t12_prior", _small_prior(table))
    base = harness.CANDIDATES["b0_scope"].build(_ctx(task, {"ref": ref}))
    changed = harness.CANDIDATES["scope_strong_null"].build(_ctx(task, {"ref": ref}))
    assert (base.X != changed.X).any()
    assert base.X.shape == changed.X.shape


# ---------------------------------------------------------------------------------------------------------------------
# The stop rule


def _harness_result(key, base, change, delta, seeds=tuple(range(10)), dist=None, cov=None, raw=None):
    rows = []
    groups = {g: 0.5 for g in stop_rule.WEIGHTS}
    metrics = {m: 0.0 for m in ("de_score", "de_direction", "mmd_u", "variogram")}
    for s in seeds:
        dg = dict(groups)
        if dist is not None:
            dg["distribution"] += dist / (100 * stop_rule.WEIGHTS["distribution"])
        if cov is not None:
            dg["covariance"] += cov / (100 * stop_rule.WEIGHTS["covariance"])
        change_metrics = {**metrics, **(raw or {})}
        rows.append({"seed": s, "scores": {base: 60.0, change: 60.0 + delta}, "groups": {base: groups, change: dg},
                     "metrics": {base: metrics, change: change_metrics}})
    return {"task": {"key": key}, "rows": rows}


def _task(key, delta=0.1, **kw):
    return stop_rule.paired(_harness_result(key, *stop_rule.PAIRS[key], delta, **kw), *stop_rule.PAIRS[key])


def _clean():
    return {k: _task(k) for k in stop_rule.PAIRS}


def test_a_clean_run_proceeds():
    v = stop_rule.verdict(_clean())
    assert v["verdict"] == "proceed"
    assert not (v["S1_failures"] or v["S2_failures"] or v["captain_band"] or v["problems"])


def test_a_total_loss_below_the_stop_bound_stops():
    t = _clean()
    t["T1/E9.5->E10.5"] = _task("T1/E9.5->E10.5", -0.7)
    v = stop_rule.verdict(t)
    assert v["verdict"] == "STOP" and "T1/E9.5->E10.5" in v["S1_failures"][0] and "total" in v["S1_failures"][0]


def test_a_distribution_loss_cannot_hide_behind_a_de_score_gain():
    # a strong-null truth: the total is up (a de_score gain), the distribution and covariance groups together are down
    t = _clean()
    t["T1/E9.5->E10.5:chamber"] = _task("T1/E9.5->E10.5:chamber", 1.0, dist=-0.5, cov=-0.3)
    v = stop_rule.verdict(t)
    assert v["verdict"] == "STOP" and "distribution+covariance" in v["S1_failures"][0]
    assert v["S1_failures"][0].startswith("T1/E9.5->E10.5:chamber")


def test_a_loss_between_the_approved_price_and_the_stop_bound_goes_to_the_captain():
    for k, kw in (("T1/E9.5->E10.5", {"delta": -0.5}), ("T1/E9.5->E10.5:chamber", {"delta": 0.2, "dist": -0.4}),
                  (stop_rule.CONTROL_TASK, {"delta": -0.4})):
        t = _clean()
        t[k] = _task(k, **kw)
        v = stop_rule.verdict(t)
        assert v["verdict"] == "captain", (k, kw, v)
        assert v["captain_band"] and k in v["captain_band"][0] and not v["S1_failures"]


def test_edges_are_not_below_the_bound():
    t = _clean()
    t["T1/E9.5->E10.5:chamber"] = _task("T1/E9.5->E10.5:chamber", -0.66)  # at the stop bound: not below it, so the band
    v = stop_rule.verdict(t)
    assert v["verdict"] == "captain" and not v["S1_failures"]
    t["T1/E9.5->E10.5:chamber"] = _task("T1/E9.5->E10.5:chamber", -0.33)  # at the approved price: not below it
    assert stop_rule.verdict(t)["verdict"] == "proceed"


def test_the_controls_raw_bounds_stop_on_their_own():
    for raw, needle in (({"mmd_u": 0.004}, "raw mmd_u"), ({"variogram": 0.0005}, "raw variogram")):
        t = _clean()
        t[stop_rule.CONTROL_TASK] = _task(stop_rule.CONTROL_TASK, 0.0, raw=raw)
        v = stop_rule.verdict(t)
        assert v["verdict"] == "STOP" and needle in v["S2_failures"][0]
    t = _clean()  # just under the raw bounds passes, and the raw increases are reported
    t[stop_rule.CONTROL_TASK] = _task(stop_rule.CONTROL_TASK, 0.0,
                                      raw={"mmd_u": 0.99 * stop_rule.RAW_MMD_BOUND,
                                           "variogram": 0.99 * stop_rule.RAW_VARIOGRAM_BOUND})
    v = stop_rule.verdict(t)
    assert v["verdict"] == "proceed" and v["control_raw_increase"]["mmd_u"] == pytest.approx(0.99 * 0.0035)


def test_a_raw_bound_elsewhere_than_the_control_does_not_stop():
    t = _clean()
    t["T1/E9.5->E10.5"] = _task("T1/E9.5->E10.5", 0.1, raw={"mmd_u": 0.01, "variogram": 0.01})
    assert stop_rule.verdict(t)["verdict"] == "proceed"


def test_the_controls_points_bound_as_first_written_still_stops():
    t = _clean()
    t[stop_rule.CONTROL_TASK] = _task(stop_rule.CONTROL_TASK, 0.0, dist=-0.34)
    v = stop_rule.verdict(t)
    assert v["verdict"] == "STOP" and v["control_distributional_points"] == pytest.approx(-0.34)
    assert "distribution+covariance" in v["S2_failures"][0] and not v["S1_failures"]
    t[stop_rule.CONTROL_TASK] = _task(stop_rule.CONTROL_TASK, 0.0, dist=-0.32)
    assert stop_rule.verdict(t)["verdict"] == "proceed"


def test_an_incomplete_run_stops():
    t = _clean()
    t["T1/E9.5->E10.5"] = _task("T1/E9.5->E10.5", 0.1, seeds=(0, 1))
    v = stop_rule.verdict(t)
    assert v["verdict"] == "STOP" and v["problems"]
    t = _clean()
    del t[stop_rule.CONTROL_TASK]
    assert stop_rule.verdict(t)["verdict"] == "STOP"


def test_a_non_finite_score_stops_and_is_named():
    for mutate, expect in ((lambda r: r["scores"].__setitem__("b0_strong_null", float("nan")), "non-finite"),
                           (lambda r: r["groups"]["b0_strong_null"].__setitem__("covariance", float("nan")),
                            "non-finite"),
                           (lambda r: r["metrics"]["b0_strong_null"].__setitem__("mmd_u", float("inf")), "non-finite")):
        res = _harness_result("T1/E9.5->E10.5", "b0", "b0_strong_null", 0.1)
        res["rows"] = [{**r, "scores": dict(r["scores"]), "groups": {k: dict(v) for k, v in r["groups"].items()},
                        "metrics": {k: dict(v) for k, v in r["metrics"].items()}} for r in res["rows"]]
        mutate(res["rows"][0])
        t = _clean()
        t["T1/E9.5->E10.5"] = stop_rule.paired(res, "b0", "b0_strong_null")
        v = stop_rule.verdict(t)
        assert v["verdict"] == "STOP" and v["problems"] and expect in v["problems"][0], v


def test_a_non_finite_control_is_not_read_as_a_pass_of_its_raw_bounds():
    res = _harness_result(stop_rule.CONTROL_TASK, *stop_rule.PAIRS[stop_rule.CONTROL_TASK], 0.0)
    res["rows"] = [{**r, "metrics": {k: dict(v) for k, v in r["metrics"].items()}} for r in res["rows"]]
    res["rows"][0]["metrics"]["scope_strong_null"]["mmd_u"] = float("nan")
    t = _clean()
    t[stop_rule.CONTROL_TASK] = stop_rule.paired(res, *stop_rule.PAIRS[stop_rule.CONTROL_TASK])
    v = stop_rule.verdict(t)
    assert v["verdict"] == "STOP" and v["problems"] and v["control_raw_increase"] is None


def test_implied_points_are_reported_for_several_d_f():
    imp = stop_rule.implied_points(0.0035, 0.0004)
    assert imp["mmd_u"]["0.08"] == pytest.approx(-7.5 / 0.08 * 0.0035)
    assert imp["mmd_u"]["0.024"] == pytest.approx(-7.5 / 0.024 * 0.0035)
    assert imp["variogram"]["0.0055"] == pytest.approx(-5.0 / 0.0055 * 0.0004)
    assert set(imp["mmd_u"]) == {str(d) for d in stop_rule.DF_MMD}
    assert set(imp["variogram"]) == {str(d) for d in stop_rule.DF_VARIOGRAM}
    # the raw mmd_u bound is 0.33 points near the floor at d_f 0.08 (the row's statement), and more at a smaller d_f
    assert imp["mmd_u"]["0.08"] == pytest.approx(-0.328, abs=0.001) and imp["mmd_u"]["0.024"] < -1.0


def test_stop_rule_refuses_to_overwrite(tmp_path):
    (tmp_path / "x").write_text("y")
    with pytest.raises(SystemExit, match="already written"):
        stop_rule.main(["--runs", str(tmp_path), "--out", str(tmp_path)])


def test_the_prior_reads_nothing_refused_when_the_harness_imports_it_under_the_guard():
    code = ("import os, sys; sys.path.insert(0, os.path.abspath('src')); from vec import harness, t1_scope;"
            "import numpy as np, anndata as ad, pandas as pd;"
            "a = ad.AnnData(X=np.ones((3, 2), dtype=np.float32), "
            "obs=pd.DataFrame({'celltype': ['V-CM', 'Blood', 'BEC']}, index=['a', 'b', 'c']));"
            "from types import SimpleNamespace as N;"
            "ctx = N(inputs={'ref': a}, task=N(reference_role='ref'));"
            "\nwith harness.guarded_reads() as g:\n"
            "    labels = harness._t12_labels(ctx)\n"
            "print(sorted(labels), g['refused'])")
    out = subprocess.run([sys.executable, "-I", "-B", "-c", code], cwd=Path(__file__).resolve().parents[1],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip().splitlines()[-1] == "['BEC', 'V-CM'] []"


def _full_result(key, cands, deltas, commit="c0", clean=True, raw=None):
    """A harness result (the fields check_results reads) for stop_rule.main."""
    rows = []
    groups = {g: 0.5 for g in stop_rule.WEIGHTS}
    metrics = {m: 0.0 for m in ("de_score", "de_direction", "mmd_u", "variogram")}
    for s in range(10):
        scores = {c: 50.0 + deltas.get(c, 0.0) for c in cands}
        per_candidate = {c: ({**metrics, **(raw or {})} if deltas.get(c) is not None else metrics) for c in cands}
        rows.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "scores": scores,
                     "groups": {c: groups for c in cands}, "metrics": per_candidate})
    return {"task": {"key": key, "test_board": "T1:test"}, "rows": rows, "candidates": list(cands),
            "protocol": harness.PROTOCOL, "build": harness.BUILD, "anchors": "computed", "frac": harness.SUBSAMPLE_FRAC,
            "seeds": list(range(10)),
            "provenance": {"commit": commit, "tree_clean": clean, "config_sha256": "cf", "files_sha256": {},
                           "derived_sha256": {}}}


def test_stop_rule_main_reads_the_runs_and_writes_a_stamped_verdict(tmp_path, monkeypatch, own_repo):
    runs = tmp_path / "runs"
    runs.mkdir()
    for key, name in (("T1/E9.5->E10.5", "a"), ("T1/E9.5->E10.5:chamber", "b"), (stop_rule.CONTROL_TASK, "c")):
        base, change = stop_rule.PAIRS[key]
        res = _full_result(key, ("floor", base, change), {change: -0.1},
                           raw={"mmd_u": 0.001, "variogram": 0.0001} if key == stop_rule.CONTROL_TASK else None)
        (runs / f"{name}.json").write_text(json.dumps(res))
    monkeypatch.setattr(harness, "_sha256", lambda p: "cf")
    out = tmp_path / "stop"
    stop_rule.main(["--runs", str(runs), "--out", str(out)])
    res = json.loads((out / "stop_rule.json").read_text())
    assert res["verdict"]["verdict"] == "proceed" and res["stamp"]["commit"] and len(res["tasks"]) == 3
    assert set(res["inputs_sha256"]) == {"a.json", "b.json", "c.json"}
    assert res["verdict"]["control_raw_increase"]["mmd_u"] == pytest.approx(0.001)
    implied = res["implied_points_near_floor"][stop_rule.CONTROL_TASK]
    assert implied["mmd_u"]["0.08"] == pytest.approx(-7.5 / 0.08 * 0.001)
    text = (out / "stop_rule.txt").read_text()
    assert "T12a/T12b stop rule: proceed" in text and "implied points near the floor" in text and "d_f 0.024" in text
    # a dirty run is refused (check_results), nothing is written
    out2 = tmp_path / "stop2"
    res = _full_result("T1/E9.5->E10.5", ("floor", "b0", "b0_strong_null"), {}, clean=False)
    (runs / "a.json").write_text(json.dumps(res))
    with pytest.raises(harness.InconsistentResults, match="dirty"):
        stop_rule.main(["--runs", str(runs), "--out", str(out2)])
    assert not out2.exists()


def test_replace_constructs_candidate_board(monkeypatch):
    # guards the record: a board copied with replace keeps its genes (used by the candidates' synthetic context)
    b = Board(key="pseudo", task="T1", genes=list(GENES), needs_coords=False, min_cells=1, max_cells=2)
    assert replace(b, max_cells=3).genes == list(GENES)
