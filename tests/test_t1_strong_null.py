"""T12 (experiments/harness/TRIES.md), on synthetic data only: the strong-null down-tail prior's bins, table, gene
order and realisation (`vec.t1_scope`), its backtest candidates (`vec.harness`), the build option
(`experiments/builds/build_test.py`) and the calibration and stop-rule scripts (`experiments/t1_strong_null`)."""

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
from vec.evaluate import Stage
from vec.harness import BuildContext, PseudoTask, TaskView
from vec.prep import log1p_cp10k

EXPERIMENT = Path(__file__).resolve().parents[1] / "experiments" / "t1_strong_null"
sys.path.insert(0, str(EXPERIMENT))
import calibrate  # noqa: E402
import stop_rule  # noqa: E402
import strong_null  # noqa: E402

# ---------------------------------------------------------------------------------------------------------------------
# The rule's numbers (T12's row)


def test_the_constants_are_the_row():
    assert ts.DP_LEVEL_EDGES == (100, 200, 400, 800, 1600, 3200)
    assert ts.DP_CHANGE_EDGES == (400, 1600) and ts.DP_CHANGE_TOP == 1600 and ts.DP_PROMOTED_CHANGE == (1, 2)
    assert ts.DP_SHRINK == 100.0 and ts.DP_BUDGET == 200 and ts.DP_BAND == -0.15
    assert ts.DP_F_MIN == 0.5 and ts.DP_MIN_PB == 0.25 and ts.DP_ITERATIONS == 3 and ts.DP_TOL == 0.02
    assert ts.DP_DEPTH_TOL == 0.01 and ts.DP_BOTTOM == 400
    assert stop_rule.LOSS_BOUND == -0.66 and stop_rule.DIST_BOUND == -0.33
    assert stop_rule.PAIRS == {"T1/E9.5->E10.5": ("b0", "b0_strong_null"),
                               "T1/E9.5->E10.5:chamber": ("b0", "b0_strong_null"),
                               "T1/E8.5->E9.5": ("b0_scope", "scope_strong_null")}
    assert harness.T12_GROUPS == ("CM", "Endocardium", "Endothelium", "Epicardium")
    assert calibrate.FRACTION == {9.5: 0.2}


# ---------------------------------------------------------------------------------------------------------------------
# Bins, table, gene order


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


def test_table_shrinks_each_cell_towards_its_levels_marginal():
    h = np.array([[10.0, 0.0, 30.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
    n = np.array([[20.0, 0.0, 180.0, 0.0], [50.0, 0.0, 0.0, 0.0]])
    t = ts.DownPriorTable.from_counts(h, n, shrink=100.0, level_edges=(5,), change_edges=(1, 2), change_top=1)
    m = 40.0 / 200.0
    assert t.p[0, 0] == pytest.approx((10.0 + 100.0 * m) / 120.0)
    assert t.p[0, 1] == pytest.approx(m)  # an empty cell is its level's marginal
    assert t.p[0, 2] == pytest.approx((30.0 + 100.0 * m) / 280.0)
    assert t.p[1, 0] == pytest.approx(0.0)  # a level with no hits, including an empty cell
    assert ts.DownPriorTable.from_dict(json.loads(json.dumps(t.to_dict()))).to_dict() == t.to_dict()


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
    kw = {"budget": 5, "band": -0.2, "f_min": 0.3, "tol": 0.05, **kw}
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
    assert rep["checks"]["passed"] and rep["checks"]["band_ok"] and rep["checks"]["expected_depth_ok"]
    assert rep["checks"]["expected_depth_max_error"] < 1e-9
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


# ---------------------------------------------------------------------------------------------------------------------
# build_test.py t1


def test_strong_null_option_is_off_by_default_and_refused_without_its_pin(bt):  # noqa: F811
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "runs/x/T1_s1.h5ad"]).strong_null is False
    assert ap.parse_args(["t1", "--strong-null", "--out", "runs/x/T1_s1.h5ad"]).strong_null is True
    assert bt._t1_down_prior(False, True, True) is None  # off: nothing is read
    with pytest.raises(SystemExit, match="without R1 and R2"):
        bt._t1_down_prior(True, False, True)
    with pytest.raises(SystemExit, match="without R1 and R2"):
        bt._t1_down_prior(True, True, False)
    with pytest.raises(SystemExit, match="pins no table_sha256"):  # the repository's config pins none yet
        bt._t1_down_prior(True, True, True)


def test_the_repository_config_pins_no_table_yet():
    assert strong_null.load_config()["table_sha256"] is None
    assert strong_null.load_config()["table_file"].endswith("results/calibration/table_strong_chamber.json")


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


def test_option_on_moves_only_the_covered_cells(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    # a small table for the 62-gene synthetic world: levels and changes on its own scale
    table = ts.DownPriorTable(p=np.array([[0.0, 0.9, 0.8, 0.0], [0.0, 0.7, 0.6, 0.0], [0.0, 0.1, 0.1, 0.0]]),
                              level_edges=(20, 40), change_edges=(12, 30), change_top=8)
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: {"file": "f", "sha256": "ab" * 32, "genes": ["g5"],
                                                                   "modules": {"D": 1}, "sources": {}})
    monkeypatch.setattr(bt.strong_null, "load_table", lambda *a, **k: table)
    monkeypatch.setattr(bt.strong_null, "table_record", lambda *a, **k: {"file": "t.json", "sha256": "cd" * 32})
    monkeypatch.setattr(bt.t1_scope, "DownPrior", partial(ts.DownPrior, budget=4, band=-0.3, f_min=0.2, tol=0.2))
    off, on = tmp_path / "off.h5ad", tmp_path / "on.h5ad"
    bt.main(["t1", "--out", str(off)])
    t1_env.reads.clear()
    bt.main(["t1", "--strong-null", "--out", str(on)])
    assert t1_env.reads == [("Feng", [""]), "E10.5", ("Feng", T1_GENES)]  # nothing read before the prior's refusals
    A, B = ad.read_h5ad(off).X, ad.read_h5ad(on).X
    A, B = A.toarray() if sparse.issparse(A) else A, B.toarray() if sparse.issparse(B) else B
    assert A.shape == B.shape
    rec = json.loads(on.with_suffix(".build.json").read_text())
    assert rec["strong_null"]["table"] == {"file": "t.json", "sha256": "cd" * 32}
    assert rec["strong_null"]["checks"]["passed"] is True
    assert rec["strong_null"]["n_promoted"] == 4 and rec["n_cells"] == "max_cells"
    assert str(t1_env.board.max_cells) not in on.with_suffix(".build.json").read_text()
    # the same cells drawn and ordered: only rows that hold covered cells differ, and only some hold promoted genes
    changed = (A != B).any(axis=1)
    assert changed.any() and not changed.all()  # R1's cells (blood, macrophages, Early VIC) are the unchanged rows
    assert (A[~changed] == B[~changed]).all()
    # the promoted genes are thinner in the on file
    idx = [T1_GENES.index(g) for g in rec["strong_null"]["promoted"]]
    assert (B[:, idx].mean(0) < A[:, idx].mean(0)).all()


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


def test_candidates_apply_where_their_inputs_allow(monkeypatch):
    tasks = {t.key: t for t in harness.load_tasks()["T1:test"]}
    b0 = harness.CANDIDATES["b0_strong_null"]
    control = harness.CANDIDATES["scope_strong_null"]
    assert b0.applies(tasks["T1/E9.5->E10.5"]) and b0.applies(tasks["T1/E9.5->E10.5:chamber"])
    assert not b0.applies(tasks["T1/E8.5->E9.5"])  # its table's sources are not that task's inputs
    assert control.applies(tasks["T1/E8.5->E9.5"]) and not control.applies(tasks["T1/E9.5->E10.5"])
    d = harness.DERIVED["t1_down_prior_weak"]
    assert d.path.name == "table_weak.json" and d.sources() == ("T1/E8.5_RNA.h5ad", "T1/E9.5_RNA.h5ad")
    assert b0.params["scope"] == "chamber" == control.params["scope"]
    assert not harness.derived_ok("t1_down_prior_weak", tasks["T1/E8.5->E9.5"])
    assert harness.derived_ok("t1_down_prior_weak", tasks["T1/E9.5->E10.5"])


def test_down_prior_load_round_trips_the_committed_form(tmp_path):
    t = ts.DownPriorTable.from_counts(np.ones((7, 4)), 2 * np.ones((7, 4)))
    p = tmp_path / "t.json"
    p.write_text(json.dumps({**t.to_dict(), "extra": 1}))
    assert harness._down_prior_load(p).to_dict() == t.to_dict()


def test_candidate_shares_b0s_draw_and_moves_only_its_cells(monkeypatch):
    task = PseudoTask(key="T1/E9.5->E10.5", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      validation_board="T1:val")
    ref, start = _labelled(1500, LABELS, 1, GENES), _labelled(1200, LABELS, 2, GENES)
    table = ts.DownPriorTable(p=np.array([[0.0, 0.9, 0.8, 0.0], [0.0, 0.7, 0.6, 0.0], [0.0, 0.1, 0.1, 0.0]]),
                              level_edges=(15, 30), change_edges=(10, 25), change_top=5)
    derived = {"t1_held_flat": ("g49",), "t1_down_prior_weak": table}
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
    monkeypatch.setattr(harness, "_t12_control", _small_prior(table))
    base = harness.CANDIDATES["b0_scope"].build(_ctx(task, {"ref": ref}))
    changed = harness.CANDIDATES["scope_strong_null"].build(_ctx(task, {"ref": ref}))
    assert (base.X != changed.X).any()
    assert base.X.shape == changed.X.shape


def test_candidates_are_refused_a_data_read_under_the_guard(monkeypatch):
    # the table comes as a derived input (read by the harness); a candidate reading the file itself is refused
    task = PseudoTask(key="T1/E9.5->E10.5", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      validation_board="T1:val")
    with harness.guarded_reads(), pytest.raises(PermissionError, match="refused"):
        Path(harness.DERIVED["t1_down_prior_weak"].path).read_text()
    assert harness.CANDIDATES["b0_strong_null"].derived_for(task) == ("t1_held_flat", "t1_down_prior_weak")


# ---------------------------------------------------------------------------------------------------------------------
# The calibration script and the stop rule


def _stages(n_genes=3600, seed=0):
    rng = np.random.default_rng(seed)
    ref = rng.poisson(rng.gamma(1.0, 2.0, n_genes), (300, n_genes)).astype(np.float32)
    target = ref[:260].copy()
    # three highly expressed genes fall in the target (the truth's down list)
    mu = rng.gamma(1.0, 2.0, n_genes)
    base = rng.poisson(mu, (300, n_genes)).astype(np.float32)
    base[:, :3] += 60.0
    ref = np.log1p(base)
    tgt = base[:260].copy()
    tgt[:, :3] = np.round(tgt[:, :3] * 0.25)
    target = np.log1p(tgt)
    return Stage(X=target.astype(np.float32)), Stage(X=ref.astype(np.float32))


def test_calibrate_counts_the_truths_down_list_in_the_right_cells():
    target, ref = _stages()
    pb_ref = np.asarray(ref.X, dtype=np.float64).mean(0)
    noise = np.random.default_rng(1).normal(0.0, 1.0, pb_ref.size)
    noise[:3] = 0.0  # genes 0-2 sit in the middle of the change bins
    res = calibrate.calibrate_transition(target, ref, lambda seed: pb_ref + noise, seeds=(0, 1, 2))
    p = np.asarray(res["p"])
    assert p.shape == (7, 4) and res["seeds"] == [0, 1, 2] and res["shrink"] == ts.DP_SHRINK
    hits, genes = np.asarray(res["hits"]), np.asarray(res["genes"])
    assert genes.sum() == 3 * target.X.shape[1] and hits.sum() == sum(r["n_dn"] for r in res["per_seed"])
    # genes 0-2 are the highest level (level bin 0) and always in the down list; dp is 0 there, the middle of dp's
    # ranks, so change bin 2 (flat)
    assert hits[0, 2] >= 3 * 3 and p[0, 2] > p[3, 2]
    # the expected number of down-gene-seeds matches what the scorer returned for those seeds
    assert all(r["n_dn"] >= 3 for r in res["per_seed"]) and res["mean_n_dn"] >= 3


def test_calibrate_refuses_to_overwrite(tmp_path):
    (tmp_path / "x").write_text("y")
    with pytest.raises(SystemExit, match="already written"):
        calibrate.main(["--out", str(tmp_path)])


def _harness_result(key, base, change, delta, seeds=tuple(range(10)), dist=None):
    rows = []
    groups = {g: 0.5 for g in stop_rule.WEIGHTS}
    metrics = {m: 0.0 for m in ("de_score", "de_direction", "mmd_u", "variogram")}
    for s in seeds:
        dg = dict(groups)
        if dist is not None:
            dg["distribution"] += dist / (100 * stop_rule.WEIGHTS["distribution"])
        rows.append({"seed": s, "scores": {base: 60.0, change: 60.0 + delta}, "groups": {base: groups, change: dg},
                     "metrics": {base: metrics, change: metrics}})
    return {"task": {"key": key}, "rows": rows}


def test_stop_rule_decides_from_the_paired_numbers():
    ok = {k: stop_rule.paired(_harness_result(k, *stop_rule.PAIRS[k], 0.1), *stop_rule.PAIRS[k])
          for k in stop_rule.PAIRS}
    assert stop_rule.verdict(ok)["verdict"] == "proceed"
    loss = dict(ok)
    loss["T1/E9.5->E10.5"] = stop_rule.paired(_harness_result("T1/E9.5->E10.5", "b0", "b0_strong_null", -0.7),
                                              "b0", "b0_strong_null")
    v = stop_rule.verdict(loss)
    assert v["verdict"] == "STOP" and v["S1_failures"] and "T1/E9.5->E10.5" in v["S1_failures"][0]
    edge = dict(ok)  # exactly at the bound is not below it
    edge["T1/E9.5->E10.5:chamber"] = stop_rule.paired(_harness_result("c", "b0", "b0_strong_null", -0.66),
                                                      "b0", "b0_strong_null")
    assert stop_rule.verdict(edge)["verdict"] == "proceed" and not stop_rule.verdict(edge)["S1_failures"]
    dist = dict(ok)
    dist[stop_rule.CONTROL_TASK] = stop_rule.paired(
        _harness_result(stop_rule.CONTROL_TASK, *stop_rule.PAIRS[stop_rule.CONTROL_TASK], 0.0, dist=-0.34),
        *stop_rule.PAIRS[stop_rule.CONTROL_TASK])
    v = stop_rule.verdict(dist)
    assert v["verdict"] == "STOP" and v["control_distributional_points"] == pytest.approx(-0.34)
    dist[stop_rule.CONTROL_TASK] = stop_rule.paired(
        _harness_result(stop_rule.CONTROL_TASK, *stop_rule.PAIRS[stop_rule.CONTROL_TASK], 0.0, dist=-0.32),
        *stop_rule.PAIRS[stop_rule.CONTROL_TASK])
    assert stop_rule.verdict(dist)["verdict"] == "proceed"
    short = dict(ok)
    short["T1/E9.5->E10.5"] = stop_rule.paired(_harness_result("c", "b0", "b0_strong_null", 0.1, seeds=(0, 1)),
                                               "b0", "b0_strong_null")
    v = stop_rule.verdict(short)
    assert v["verdict"] == "STOP" and v["problems"]


def test_stop_rule_refuses_to_overwrite(tmp_path):
    (tmp_path / "x").write_text("y")
    with pytest.raises(SystemExit, match="already written"):
        stop_rule.main(["--runs", str(tmp_path), "--out", str(tmp_path)])


def test_strong_null_table_is_refused_unless_pinned_and_round_trips_when_it_is(tmp_path):
    t = ts.DownPriorTable.from_counts(np.ones((7, 4)), 2 * np.ones((7, 4)))
    f = tmp_path / "t.json"
    f.write_text(json.dumps(t.to_dict()))
    sha = strong_null.sha256_file(f)
    with pytest.raises(strong_null.Refused, match="pins no table_sha256"):
        strong_null.table_record({"table_file": "t.json", "table_sha256": None}, tmp_path)
    with pytest.raises(strong_null.Refused, match="no such file"):
        strong_null.table_record({"table_file": "none.json", "table_sha256": "ab"}, tmp_path)
    with pytest.raises(strong_null.Refused, match="is not the one"):
        strong_null.table_record({"table_file": "t.json", "table_sha256": "ab"}, tmp_path)
    assert strong_null.table_record({"table_file": "t.json", "table_sha256": sha}, tmp_path) == {"file": "t.json",
                                                                                                 "sha256": sha}
    assert strong_null.load_table({"table_file": "t.json", "table_sha256": sha}, tmp_path).to_dict() == t.to_dict()


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


def _mini_world(tmp_path, monkeypatch):
    """Everything calibrate.main reads, synthetic: three stages (E8.5, E9.5, E10.5) in panel order and the T1:val
    board."""
    genes = [f"g{i}" for i in range(50)]
    labels = ["V-CM", "Endocardium", "Blood", "BEC"]
    stages = {name: _labelled(800, labels, seed, genes) for seed, name in enumerate(
        ("E8.5_RNA.h5ad", "E9.5_RNA.h5ad", "E10.5_RNA.h5ad"))}
    board = Board(key="T1:val", task="T1", genes=genes, needs_coords=False, min_cells=10, max_cells=300)
    monkeypatch.setattr(calibrate, "load_boards", lambda: {"T1:val": board})
    monkeypatch.setattr(calibrate.harness, "task_board", lambda t, boards, n_t, n_r: board)
    monkeypatch.setattr(calibrate, "read_stage", lambda path, g: stages[Path(path).name].copy())
    monkeypatch.setattr(calibrate, "load_stage", lambda path, lineage: stages[Path(path).name].copy())
    monkeypatch.setattr(calibrate.harness, "filter_target",
                        lambda t, data, g: (Stage.from_anndata(stages[Path(t.target).name], g), {"kept": 800}))
    monkeypatch.setattr(calibrate.hygiene, "load_union", lambda: {"genes": ["g3"], "sha256": "ab" * 32})
    monkeypatch.setattr(calibrate.t1_scope, "held_flat_genes", lambda *a, **k: ["g49"])
    monkeypatch.setattr(calibrate.hygiene, "sha256_file", lambda p: "cd" * 32)
    dirs = {g: np.random.default_rng(0).normal(0.0, 0.5, len(genes)) for g in FENG_GROUPS}
    monkeypatch.setattr(calibrate.feng_dirs, "feng_directions",
                        lambda names, *a, **k: (dirs, {"sha256": "ef" * 32, "groups": {}}))
    return tmp_path / "calibration"


def test_calibrate_main_writes_three_stamped_tables_and_no_test_board_field(tmp_path, monkeypatch, own_repo):
    out = _mini_world(tmp_path, monkeypatch)
    # the real seed list is the adoption rule's; a smaller one keeps this synthetic check fast
    real_transition = calibrate.calibrate_transition
    monkeypatch.setattr(calibrate, "calibrate_transition",
                        lambda target, ref, predict: real_transition(target, ref, predict, seeds=(0, 1)))
    calibrate.main(["--out", str(out)])
    names = {p.name for p in out.iterdir()}
    assert {"table_strong_chamber.json", "table_strong_whole.json", "table_weak.json", "calibration.json",
            "calibration.txt"} <= names and "table_weak.json.stamp.json" in names
    weak = json.loads((out / "table_weak.json").read_text())
    assert np.asarray(weak["p"]).shape == (7, 4) and weak["seeds"] == [0, 1]
    assert weak["stamp"]["commit"] and "tree_clean" in weak["stamp"] and weak["inputs_sha256"]
    assert weak["transition"]["task"] == "T1/E8.5->E9.5" and weak["transition"]["truth_filter"] is None
    strong = json.loads((out / "table_strong_chamber.json").read_text())
    assert strong["transition"]["task"] == "T1/E9.5->E10.5:chamber" and strong["predictor"].startswith("m1")
    assert ts.DownPriorTable.from_dict(strong).p.shape == (7, 4)
    assert "T1:test" not in json.dumps(strong) and "max_cells" not in (out / "calibration.txt").read_text()
    with pytest.raises(SystemExit, match="already written"):  # it runs once
        calibrate.main(["--out", str(out)])


def _full_result(key, cands, deltas, commit="c0", clean=True):
    """A harness result (the fields check_results reads) for stop_rule.main."""
    rows = []
    groups = {g: 0.5 for g in stop_rule.WEIGHTS}
    metrics = {m: 0.0 for m in ("de_score", "de_direction", "mmd_u", "variogram")}
    for s in range(10):
        scores = {c: 50.0 + deltas.get(c, 0.0) for c in cands}
        rows.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "scores": scores,
                     "groups": {c: groups for c in cands}, "metrics": {c: metrics for c in cands}})
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
        res = _full_result(key, ("floor", base, change), {change: -0.1})
        (runs / f"{name}.json").write_text(json.dumps(res))
    monkeypatch.setattr(harness, "_sha256", lambda p: "cf")
    out = tmp_path / "stop"
    stop_rule.main(["--runs", str(runs), "--out", str(out)])
    res = json.loads((out / "stop_rule.json").read_text())
    assert res["verdict"]["verdict"] == "proceed" and res["stamp"]["commit"] and len(res["tasks"]) == 3
    assert set(res["inputs_sha256"]) == {"a.json", "b.json", "c.json"}
    assert "T12 stop rule: proceed" in (out / "stop_rule.txt").read_text()
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
