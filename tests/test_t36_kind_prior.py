"""Try T36 (experiments/t1_scope/t36_kind_prior.py), on synthetic data only: the base arm is kind_gap.py's own, the
prior arm applies b0_strong_null's prior to the same draw, the reading's rule and void conditions, and the refusals. No
test scores a prediction against a target: the scorer is replaced by a stand-in where a run is exercised."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import sparse
from test_hidden_counts import E105, T1_GENES, _t1_stage

from vec import harness
from vec import t1_scope as ts
from vec.adoption import MIN_AHEAD, MIN_GAIN, SEEDS
from vec.boards import RANKED_METRICS, Board

EXPERIMENT = Path(__file__).resolve().parents[1] / "experiments" / "t1_scope"
sys.path.insert(0, str(EXPERIMENT))
import kind_gap as kg  # noqa: E402
import t36_kind_prior as t36  # noqa: E402

KIND_GAP_SHA256 = "462d5ac635a6b8fa67effd7b5cb9362f72103efe2d4fc3a40985e3fbc2e5c23d"  # taken before T36's code
SMALL_TABLE = ts.DownPriorTable.level_only(level_edges=(15, 30), change_edges=(10, 25), change_top=5)


def test_kind_gap_is_t27s_as_committed():
    assert hashlib.sha256((EXPERIMENT / "kind_gap.py").read_bytes()).hexdigest() == KIND_GAP_SHA256


def test_constants():
    assert t36.ARMS == ("lib3_base", "lib3_base_sn") and t36.REPRO_TOL == 1e-6 and t36.PRIOR_SEED_STREAM == 4
    assert MIN_GAIN == 0.3 and MIN_AHEAD == 8 and len(SEEDS) == 10
    assert t36.VERDICTS == ("the prior adds in the kind geometry", "closed")
    assert "not adoption evidence either way" in t36.NOTE and "2026-11-10 (l)" in t36.NOTE
    assert t36.T27_RESULT == EXPERIMENT / "results" / "t27_kind_gap" / "kind_gap.json"


# ---------------------------------------------------------------------------------------------------------------------
# A synthetic E10.5 file with three libraries, prepared by kind_gap.prepare itself


def _file(seed=5):
    a = _t1_stage(E105, seed)
    lib = np.random.default_rng(seed).choice(["1", "2", "3"], a.n_obs, p=[0.25, 0.25, 0.5])
    a.obs_names = [f"cell{i}_{x}" for i, x in enumerate(lib)]
    return a


@pytest.fixture
def prepared(monkeypatch, tmp_path):
    a = _file()
    calls = []
    fake = lambda *args, **kw: calls.append(args) or a  # noqa: E731
    monkeypatch.setattr(kg, "load_stage", fake)
    board = Board(key="T1:val", task="T1", genes=list(T1_GENES), needs_coords=False, min_cells=100, max_cells=5000)
    monkeypatch.setattr(kg, "load_boards", lambda: {"T1:val": board})
    kg._D.clear()
    d = t36.prepare(tmp_path)
    assert len(calls) == 1  # one load, shared by kind_gap's geometry and the prior
    assert kg.load_stage is fake  # the capture is undone
    kg._D.update(d)
    yield SimpleNamespace(d=kg._D, a=a, n=harness.count_for(d["board"], len(d["ref_rows"])))
    kg._D.clear()


def _small(d, **kw):
    """b0_strong_null's prior on this world's scale: its covered labels kept, a table and constants for 62 genes."""
    return replace(d["t36_prior"], table=SMALL_TABLE, **{"budget": 4, "band": -0.3, "f_min": 0.2, "tol": 0.2, **kw})


def test_prepare_keeps_kind_gaps_geometry(prepared):
    d, a = prepared.d, prepared.a
    lib = kg.library_of(a.obs_names)
    assert set(d["half_a"]).isdisjoint(d["half_b"]) and (lib[d["half_a"]] == "3").all()
    assert set(d["ref_rows"]) == set(range(a.n_obs)) - set(d["half_a"])
    assert d["arms"] == kg.ARMS and d["truth_depth"] == 1.0
    assert d["t36_var_names"] == list(a.var_names) and d["t36_genes"] == list(T1_GENES)


def test_base_arm_is_kind_gaps_own_bit_for_bit(prepared):
    d = prepared.d
    for seed in (0, 3, 9):
        rows = kg.arm_rows("lib3_base", prepared.n, harness.BUILD_SEED_OFFSET + seed)
        X, info = t36.arm_x("lib3_base", prepared.n, seed)
        ref = d["X"][rows]
        assert info == {} and X.shape == ref.shape and X.dtype == ref.dtype
        assert (X != ref).nnz == 0 and np.array_equal(X.indices, ref.indices) and np.array_equal(X.data, ref.data)


def test_the_prior_is_b0_strong_nulls(prepared, monkeypatch):
    d = prepared.d
    p = d["t36_prior"]
    assert p.table.p.tolist() == ts.DownPriorTable.level_only().p.tolist()
    assert (p.table.level_edges, p.table.change_edges, p.table.change_top) == (
        ts.DP_LEVEL_EDGES, ts.DP_CHANGE_EDGES, ts.DP_CHANGE_TOP)
    assert (p.budget, p.band, p.f_min, p.tol, p.mmd_bound) == (
        ts.DP_BUDGET, ts.DP_BAND, ts.DP_F_MIN, ts.DP_TOL, ts.DP_MMD_BOUND)
    labels = np.asarray(prepared.a.obs["celltype"]).astype(str)
    ctx = SimpleNamespace(task=SimpleNamespace(reference_role="ref"),
                          inputs={"ref": prepared.a[d["ref_rows"]]})
    assert p.covered_labels == harness._t12_labels(ctx) and p.covered_labels  # the reference's labels, T12's groups
    assert p.covered_labels < set(labels)
    # pb_ref: the reference's mean log1p CP10k on the panel, as build_submission takes it
    ref = prepared.a[d["ref_rows"]]
    cols = ref.var_names.get_indexer(d["board"].genes)
    want = np.asarray(sparse.csr_matrix(ref.X[:, cols], dtype=np.float64).mean(axis=0)).ravel()
    assert np.array_equal(d["t36_pb_ref"], want)
    held = harness.DERIVED["t1_held_flat"]
    assert d["t36_held"] == tuple(held.load(held.path))
    # the step is called as build_submission calls it, with the build seed's stream 4
    seen = {}
    real = ts.apply_down_prior

    def spy(pop, pb_ref, genes, covered, held_genes, prior, seed):
        seen.update(pb_ref=pb_ref, genes=genes, covered=covered, held=held_genes, prior=prior, seed=seed)
        return real(pop, pb_ref, genes, covered, held_genes, prior, seed)

    monkeypatch.setattr(ts, "apply_down_prior", spy)
    t36.arm_x("lib3_base_sn", prepared.n, 2)
    rows = kg.arm_rows("lib3_base", prepared.n, harness.BUILD_SEED_OFFSET + 2)
    assert seen["seed"] == [harness.BUILD_SEED_OFFSET + 2, 4] and seen["prior"] is p
    assert seen["genes"] == list(T1_GENES) and seen["held"] == list(d["t36_held"])
    assert np.array_equal(seen["covered"], np.isin(labels[rows], sorted(p.covered_labels)))


def test_prior_arm_draws_the_same_cells_and_moves_only_covered_cells_and_promoted_genes(prepared):
    d = prepared.d
    d["t36_prior"] = _small(d)
    d["t36_held"] = ("held_a",)
    seed = 1
    rows = kg.arm_rows("lib3_base", prepared.n, harness.BUILD_SEED_OFFSET + seed)
    out, rep, cov = t36.prior_arm(d, rows, harness.BUILD_SEED_OFFSET + seed)
    assert rep["n_promoted"] == 4 and cov.any() and not cov.all()
    assert list(out.obs[ts.LABEL_KEY]) == list(d["t36_labels"][rows])  # the same cells, in the same order
    base_X, X = d["X"][rows], sparse.csr_matrix(out.X)
    changed = np.asarray((base_X != X).sum(1)).ravel() > 0
    assert changed.any() and not changed[~cov].any()  # uncovered cells exactly as drawn
    C0, C1 = sparse.csr_matrix(d["C"][rows]), sparse.csr_matrix(out.layers["counts"])
    assert (C1[~cov] != C0[~cov]).nnz == 0
    assert (C1[:, T1_GENES.index("held_a")] != C0[:, T1_GENES.index("held_a")]).nnz == 0
    promoted = [T1_GENES.index(g) for g in rep["promoted"]]
    other = np.setdiff1d(np.arange(len(T1_GENES)), promoted)
    m0, m1 = np.asarray(C0[cov].mean(0)).ravel(), np.asarray(C1[cov].mean(0)).ravel()
    assert (m1[promoted] < m0[promoted]).all()  # the promoted genes fall in the covered cells
    big = other[m0[other] > 5.0]
    assert (m1[big] >= 0.97 * m0[big]).all()  # every other gene moves by the depth constant only
    # the arm's prediction is that population's X; what it records is a number of genes, never of cells
    X2, info = t36.arm_x("lib3_base_sn", prepared.n, seed)
    assert (X2 != X).nnz == 0 and info["n_promoted_genes"] == 4 and 0 < info["covered_share"] < 1
    assert info["recompute_max_abs"] < 1e-6  # the synthetic file stores log1p CP10k of its counts exactly
    assert str(prepared.n) not in json.dumps({k: v for k, v in info.items() if k != "checks"})


def test_zero_budget_gives_lib3_base_bit_for_bit(prepared):
    d = prepared.d
    for prior in (replace(d["t36_prior"], budget=0), _small(d, budget=0)):
        d["t36_prior"] = prior
        for seed in (0, 4):
            base, _ = t36.arm_x("lib3_base", prepared.n, seed)
            X, info = t36.arm_x("lib3_base_sn", prepared.n, seed)
            assert info["n_promoted_genes"] == 0
            assert X.shape == base.shape and (X != base).nnz == 0
            assert np.array_equal(X.toarray(), base.toarray())


def test_held_flat_list_must_come_from_earlier_stages(monkeypatch):
    d = harness.DERIVED["t1_held_flat"]
    assert set(d.sources()) == {"T1/E8.5_RNA.h5ad", "T1/E9.5_RNA.h5ad"} and len(t36.held_flat()) > 0
    monkeypatch.setitem(harness.DERIVED, "t1_held_flat", replace(d, sources=lambda: ("T1/E10.5_RNA.h5ad",)))
    with pytest.raises(ValueError, match="not earlier than E10.5"):
        t36.held_flat()


def test_unknown_arm_is_refused(prepared):
    with pytest.raises(ValueError, match="unknown arm"):
        t36.arm_x("copy", prepared.n, 0)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


METRICS = ("mmd_u", "de_score")


def _run(diffs, base=None, **kw):
    base = [60.0 + s for s in SEEDS] if base is None else base
    rows = [{"seed": s, "scores": {"lib3_base": base[i], "lib3_base_sn": base[i] + diffs[i]},
             "skills": {"lib3_base": {m: 0.5 for m in METRICS}, "lib3_base_sn": {m: 0.5 + diffs[i] / 100 for m in
                                                                                 METRICS}},
             "prior": {"n_promoted_genes": 200, "checks": {"passed": True}}} for i, s in enumerate(SEEDS)]
    return {"try": "T36", "commit": "abc", "tree_clean": True, "split_seed": kg.SPLIT_SEED, "rows": rows, **kw}


def _t27(base=None):
    base = [60.0 + s for s in SEEDS] if base is None else base
    return {"split_seed": kg.SPLIT_SEED, "rows": [{"seed": s, "scores": {"lib3_base": base[i], "copy": 50.0}}
                                                  for i, s in enumerate(SEEDS)]}


def test_reading_passes_only_above_the_gain_with_eight_seeds_ahead():
    res = t36.reading(_run([0.5] * 10), _t27())
    assert not res["void"] and res["passes"] and res["verdict"] == "the prior adds in the kind geometry"
    assert math.isclose(res["paired"]["mean"], 0.5) and res["paired"]["ahead"] == 10
    assert math.isclose(res["skill_difference"]["mmd_u"], 0.005)
    assert res["genes_promoted_per_seed"] == {str(s): 200 for s in SEEDS}
    assert res["mean_score"]["lib3_base"] == pytest.approx(64.5)
    assert res["reproduction_largest_difference"] == 0.0
    # 8 ahead with a mean above the gain passes; 7 ahead does not, however large the mean
    assert t36.reading(_run([1.0] * 8 + [-0.1] * 2), _t27())["passes"]
    res7 = t36.reading(_run([3.0] * 7 + [-0.1] * 3), _t27())
    assert not res7["passes"] and res7["verdict"] == "closed"
    # a mean at the gain exactly is not above it
    at = t36.reading(_run([0.3] * 10), _t27())
    assert not at["passes"] and at["verdict"] == "closed"
    assert t36.reading(_run([-1.0] * 10), _t27())["verdict"] == "closed"
    for r in (res, res7):
        assert "not adoption evidence either way" in t36.text(r) and "2026-11-10 (l)" in t36.text(r)


def _mutated(fn):
    run = _run([0.5] * 10)
    fn(run)
    return run


@pytest.mark.parametrize("mutate, expect", [
    (lambda r: r["rows"].pop(), "seed is missing"),
    (lambda r: r["rows"].append(dict(r["rows"][0])), "seed is missing or repeated"),
    (lambda r: r["rows"][2]["scores"].pop("lib3_base_sn"), "arm is missing"),
    (lambda r: r["rows"][2].pop("prior"), "arm is missing"),
    (lambda r: r["rows"][3]["scores"].__setitem__("lib3_base_sn", float("nan")), "non-finite"),
    (lambda r: r["rows"][3]["skills"]["lib3_base"].__setitem__("mmd_u", float("inf")), "non-finite"),
    (lambda r: r.__setitem__("tree_clean", False), "dirty tree"),
    (lambda r: r.__setitem__("commit", "unknown"), "commit is unknown"),
    (lambda r: r.__setitem__("try", "T27"), "not T36's"),
    (lambda r: r.__setitem__("split_seed", 1), "split seed"),
])
def test_reading_voids(mutate, expect):
    res = t36.reading(_mutated(mutate), _t27())
    assert res["void"] and any(expect in w for w in res["why"]) and res["verdict"].startswith("none")
    assert t36.text(res).startswith("Try T36: VOID") and "not adoption evidence" in t36.text(res)


def test_reading_voids_unless_the_base_reproduces_t27():
    base = [60.0 + s for s in SEEDS]
    off = list(base)
    off[4] += 2e-6
    res = t36.reading(_run([0.5] * 10, base=off), _t27(base))
    assert res["void"] and "does not reproduce" in res["why"][0]
    near = list(base)
    near[4] += 5e-7
    assert not t36.reading(_run([0.5] * 10, base=near), _t27(base))["void"]
    t27 = _t27(base)
    t27["rows"].pop(0)
    assert t36.reading(_run([0.5] * 10), t27)["void"]
    assert t36.reading(None, _t27())["void"] and t36.reading(_run([0.5] * 10), None)["void"]


# ---------------------------------------------------------------------------------------------------------------------
# Run and read, with the scorer replaced (no prediction meets a target), and the refusals


def _fake_score(board, pred, target, ref, seed=0, anchors="computed", **kw):
    anchors_out = {m: {"floor": 0.0, "ceiling": 1.0} for m in RANKED_METRICS["T1"]}
    if anchors == "computed":
        return SimpleNamespace(score=50.0, skills={}, anchors=anchors_out)
    x = float(sparse.csr_matrix(pred.X).sum())  # a number of the prediction alone
    return SimpleNamespace(score=x / 1e3, skills={"mmd_u": x / 1e6}, anchors=anchors_out)


def test_run_then_read(prepared, monkeypatch, tmp_path):
    d = dict(prepared.d)
    d["t36_prior"] = _small(d)
    monkeypatch.setattr(t36, "prepare", lambda data: d)
    monkeypatch.setattr(t36, "score_board", _fake_score)
    data = tmp_path / "data"
    (data / "T1").mkdir(parents=True)
    (data / kg.FILE).write_bytes(b"synthetic")
    out = tmp_path / "run"
    t36.main(["run", "--out", str(out), "--data", str(data), "--jobs", "1"])
    res = json.loads((out / "t36.json").read_text())
    assert res["try"] == "T36" and res["seeds"] == list(SEEDS) and res["arms"] == list(t36.ARMS)
    assert res["commit"] and isinstance(res["tree_clean"], bool) and res["split_seed"] == kg.SPLIT_SEED
    assert set(res["inputs"]) == {kg.FILE, "experiments/t1_scope/kind_gap.py",
                                  "experiments/t1_direction/results/held_flat_genes.csv"}
    assert res["inputs"]["experiments/t1_scope/kind_gap.py"] == KIND_GAP_SHA256
    assert all(r["prior"]["n_promoted_genes"] == 4 for r in res["rows"])
    assert all(set(r["skills"]) == set(t36.ARMS) for r in res["rows"])
    # a T27 result that holds this run's base scores: the reading is not void when the tree was clean
    t27 = {"split_seed": kg.SPLIT_SEED, "rows": [{"seed": r["seed"], "scores": {"lib3_base": r["scores"]["lib3_base"]}}
                                                for r in res["rows"]]}
    res["tree_clean"] = True
    (out / "t36.json").write_text(json.dumps(res))
    (tmp_path / "t27.json").write_text(json.dumps(t27))
    t36.main(["read", "--runs", str(out), "--t27", str(tmp_path / "t27.json"), "--out", str(tmp_path / "read")])
    reading = json.loads((tmp_path / "read" / "reading.json").read_text())
    assert not reading["void"] and reading["verdict"] in t36.VERDICTS
    assert (tmp_path / "read" / "reading.txt").read_text().startswith("Try T36 (commit")
    # the refusals: neither command writes into a folder that holds anything
    with pytest.raises(SystemExit, match="never overwrites"):
        t36.main(["run", "--out", str(out), "--data", str(data), "--jobs", "1"])
    with pytest.raises(SystemExit, match="never overwrites"):
        t36.main(["read", "--runs", str(out), "--t27", str(tmp_path / "t27.json"), "--out", str(tmp_path / "read")])
    assert reading == json.loads((tmp_path / "read" / "reading.json").read_text())


def test_run_refuses_before_loading_anything(monkeypatch, tmp_path):
    out = tmp_path / "full"
    out.mkdir()
    (out / "x").write_text("x")
    monkeypatch.setattr(t36, "prepare", lambda data: pytest.fail("prepare ran"))
    with pytest.raises(SystemExit, match="never overwrites"):
        t36.main(["run", "--out", str(out)])
    afile = tmp_path / "afile"
    afile.write_text("x")
    with pytest.raises(SystemExit, match="never overwrites"):
        t36.main(["run", "--out", str(afile)])


def test_read_of_a_missing_run_is_void(tmp_path):
    t36.main(["read", "--runs", str(tmp_path / "none"), "--t27", str(tmp_path / "none.json"),
              "--out", str(tmp_path / "read")])
    res = json.loads((tmp_path / "read" / "reading.json").read_text())
    assert res["void"] and len(res["why"]) == 2
