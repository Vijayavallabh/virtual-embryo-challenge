"""experiments/calibration/calibrate.py (feat-005, A2): the fixed rules, on synthetic data."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from conftest import make_board

from vec.evaluate import Stage, subsample

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cal():
    spec = importlib.util.spec_from_file_location("calibration_calibrate",
                                                  REPO / "experiments" / "calibration" / "calibrate.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_anchor_match_uses_the_published_decimals(cal):
    assert cal._decimals(0.08359) == 5 and cal._decimals(0.0) == 1 and cal._decimals(1.5e-05) == 6
    assert cal.anchor_match(0.083594, 0.08359) == "exact"
    assert cal.anchor_match(0.0837, 0.08359) == "close"  # within 2%
    assert cal.anchor_match(0.09, 0.08359) == "off"
    assert cal.anchor_match(None, 0.08359) == "missing" and cal.anchor_match(float("nan"), 0.1) == "missing"


def test_counts_follow_the_ten_percent_working_copy(cal):
    board = make_board("T2", True)
    board.truth_cells, board.ref_cells = 2320, 2483
    assert cal.counts_ok(board, 46396, 24826)["ok"]  # T3:gata4 with replicate 2 as the truth
    assert not cal.counts_ok(board, 42257, 24826)["ok"]  # replicate 1 does not halve to truth_cells


def _stage(n, rng, types=("a", "b")):
    return Stage(X=rng.gamma(2.0, 1.0, size=(n, 6)).astype(np.float32), coords=None,
                 celltype=np.array([types[i % len(types)] for i in range(n)]))


def test_working_copy_baselines_put_copy_last_on_the_floor(cal):
    rng = np.random.default_rng(0)
    stages = {"target": _stage(400, rng), "ref": _stage(300, rng), "prev": _stage(200, rng)}
    preds = cal.predictions("T1", stages, {}, [f"g{i}" for i in range(6)], seed=3)
    floor = subsample(stages["ref"], cal.FRAC, 3)
    assert np.array_equal(preds["copy_last"]["wc"].X, floor.X)
    wc, full = preds["pseudobulk_shift"]["wc"], preds["pseudobulk_shift"]["full_delta"]
    assert wc.n == full.n == floor.n and not np.allclose(wc.X, full.X)
    ref, prev = stages["ref"], stages["prev"]
    d = ref.X[ref.celltype == "a"].mean(0) - prev.X[prev.celltype == "a"].mean(0)
    sel = floor.celltype == "a"
    assert np.allclose(full.X[sel], np.clip(floor.X[sel] + d, 0, None))


def test_gene_ko_zeroes_the_target_gene_on_both_constructions(cal):
    rng = np.random.default_rng(1)
    genes = [f"g{i}" for i in range(6)]
    stages = {"target": _stage(400, rng), "ref": _stage(300, rng), "ko_train": _stage(250, rng),
              "wt_train": _stage(260, rng)}
    preds = cal.predictions("T3", stages, {"target_gene": "g2"}, genes, seed=0)
    assert set(preds) == {"wt_identity", "shift_transfer", "gene_ko"}
    for c in cal.CONSTRUCTIONS:
        assert np.all(preds["gene_ko"][c].X[:, 2] == 0)
        assert np.array_equal(np.delete(preds["gene_ko"][c].X, 2, 1), np.delete(preds["shift_transfer"][c].X, 2, 1))


def test_summarise_applies_the_fixed_verdict(cal, tmp_path):
    runs = tmp_path / "runs"
    runs.mkdir()
    rows = []
    for c in cal.CONSTRUCTIONS:
        for s in cal.SEEDS:
            rows.append({"model": "pseudobulk_shift", "construction": c, "seed": s, "anchors": "computed",
                         "score": (48.6 if c == "wc" else 45.0) + 0.1 * (s % 2)})
        rows.append({"model": "pseudobulk_shift", "construction": c, "seed": 0, "anchors": "published",
                     "score": 48.58 if c == "wc" else 45.0})
    exact = {"m": {"floor": {"published": 0.1, "computed": 0.1, "match": "exact"},
                   "ceiling": {"published": 0.01, "computed": 0.01, "match": "exact"}}}
    for key in cal.VALIDATION_BOARDS:
        (runs / f"{key.replace(':', '__')}.json").write_text(json.dumps(
            {"board": key, "commit": "abc1234", "tree_clean": True, "counts": {"ok": True},
             "anchors_seed0": exact, "rows": rows if key == "T1:val" else []}))
    pub = tmp_path / "rows.json"
    pub.write_text(json.dumps({"rows": [{"board": "T1:val", "model": "pseudobulk_shift", "score": 48.58}]}))
    s = cal.summarise(pub, runs, tmp_path / "out")
    assert s["verdict_A2"] == "pass" and s["anchors_ok"]
    best = [c for c in s["comparisons"] if c.get("best")]
    assert len(best) == 1 and best[0]["construction"] == "wc" and best[0]["exact"]

    pub.write_text(json.dumps({"rows": [{"board": "T1:val", "model": "pseudobulk_shift", "score": 52.0}]}))
    assert cal.summarise(pub, runs, tmp_path / "out")["verdict_A2"] == "fail"  # 3.3 points off
    pub.write_text(json.dumps({"rows": [{"board": "T1:val", "model": "pseudobulk_shift", "score": 48.58},
                                        {"board": "T1:val", "model": "neural_ode", "score": 47.0},
                                        {"board": "T2:heart:val_interp", "model": "pseudobulk_shift", "score": 49.0}]}))
    s = cal.summarise(pub, runs, tmp_path / "out")  # neural_ode is not rebuilt; heart interp has no such row here
    assert s["rows_out_of_scope"] == [["T1:val", "neural_ode"]]
    assert s["rows_not_scored"] == [["T2:heart:val_interp", "pseudobulk_shift"]] and s["verdict_A2"] == "incomplete"
