"""Scorer builds in the emulator: "prefix" (veckit 0.1.1, pre-PR-#9 canonicalisation) and "main" (pinned cc3859a).

Synthetic clouds only. The helix below is strongly chiral with well-separated PCA axes, so its proper-mode and
mirror-mode Dice are far apart (about 0.25 at 2,000 cells per half), and the pre-fix SVD parity of its rows flips with
row order about half the time.
"""

from __future__ import annotations

import importlib.util
import inspect
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from conftest import GENES, make_board

import vec.evaluate as ev
from vec.boards import Board
from vec.evaluate import (
    SCORER_BUILDS,
    Stage,
    handedness,
    mirror_z,
    prefix_parity,
    score_board,
    score_board_seeds,
    scorer_build,
    split_half,
    subsample,
    svd_parity,
)

REPO = Path(__file__).resolve().parents[1]
N_GENES = len(GENES)
CANONICALISING = {"occupancy_dice", "sliced_wasserstein", "_sw_flip_spread"}


def _sm():
    import common.shape_metrics as sm

    return sm


def helix_coords(n, rng, scale=1.0, noise=0.35):
    s = rng.uniform(0.0, 2.5 * np.pi, n)
    C = np.c_[3.0 * np.cos(s), 1.5 * np.sin(s), 0.8 * s] * scale + rng.normal(scale=noise * scale, size=(n, 3))
    return C.astype(np.float32), s


def helix_halves(seed, n=2000):
    C, _ = helix_coords(2 * n, np.random.default_rng(seed))
    return C[:n], C[n:]  # truth, prediction


def helix_stage(n, rng, means_by_type, scale=1.0):
    coords, s = helix_coords(n, rng, scale)
    types = np.minimum((s / (2.5 * np.pi) * len(means_by_type)).astype(int), len(means_by_type) - 1)
    counts = rng.poisson(means_by_type[types])
    lib = counts.sum(1, keepdims=True).clip(min=1)
    X = np.log1p(1e4 * counts / lib).astype(np.float32)
    return Stage(X=X, coords=coords, celltype=np.array([f"type{t}" for t in types]))


@pytest.fixture(scope="module")
def helix_stages():
    rng = np.random.default_rng(20261019)
    base = rng.gamma(2.0, 2.0, size=(3, N_GENES))
    later = base.copy()
    later[:, :8] *= 4.0
    later[:, 8:16] *= 0.2
    ref = helix_stage(4000, np.random.default_rng(1), base, scale=0.8)
    target = helix_stage(4000, np.random.default_rng(2), later)
    pred = helix_stage(4000, np.random.default_rng(3), later)
    return ref, target, pred


def same(a, b) -> bool:
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


def frame_det(canonicalise, C) -> float:
    """Determinant of the rotation a canonicalisation applied (recovered from its output)."""
    C = np.asarray(C, float)
    X = C - C.mean(0)
    Z, rms = canonicalise(C)
    R, *_ = np.linalg.lstsq(X, Z * rms, rcond=None)
    return float(np.linalg.det(R))


# ------------------------------------------------------------------------------------------------ the two builds
def test_builds_differ_only_by_the_right_handed_flip():
    sm = _sm()
    pinned = ev._pinned_canonicalise()
    assert pinned.__module__ == "common.shape_metrics"
    assert "np.linalg.det(Vt) < 0" in inspect.getsource(pinned), "the pinned veckit is not the fixed build"
    with scorer_build("main"):
        assert sm._canonicalise is pinned
    with scorer_build("prefix"):
        assert sm._canonicalise is ev._canonicalise_prefix
        with scorer_build("main"):
            assert sm._canonicalise is pinned
        assert sm._canonicalise is ev._canonicalise_prefix
    assert sm._canonicalise is pinned

    rng = np.random.default_rng(0)
    parities = set()
    for _ in range(10):
        C, _ = helix_coords(500, rng)
        for cloud in (C, mirror_z(C)):
            s = svd_parity(cloud)
            parities.add(s)
            Zm, rms_m = pinned(cloud)
            Zp, rms_p = ev._canonicalise_prefix(cloud)
            assert rms_m == rms_p
            assert np.array_equal(Zp, Zm * np.array([1.0, 1.0, s]))  # the flip is the only difference, bit for bit
            assert frame_det(pinned, cloud) == pytest.approx(1.0)  # main: always right-handed
            assert frame_det(ev._canonicalise_prefix, cloud) == pytest.approx(s)  # prefix: as LAPACK returns it
    assert parities == {-1, 1}


def test_bad_build_is_rejected_and_the_pinned_function_restored(helix_stages):
    sm = _sm()
    pinned = ev._pinned_canonicalise()
    with pytest.raises(ValueError):
        with scorer_build("0.1.2"):
            pass
    with pytest.raises(RuntimeError, match="boom"):
        with scorer_build("prefix"):
            raise RuntimeError("boom")
    assert sm._canonicalise is pinned
    ref, target, pred = helix_stages
    with pytest.raises(ValueError):
        score_board(make_board("T2", True), pred, target, ref, build="pre")


@pytest.mark.parametrize("seed", range(4))
def test_mirror_twin_scores_the_same_under_prefix_and_differently_under_main(seed):
    truth, pred = helix_halves(seed)
    twin = mirror_z(pred)
    sm = _sm()
    with scorer_build("prefix"):
        assert sm.occupancy_dice(twin, truth)[0] == sm.occupancy_dice(pred, truth)[0]
        assert sm.sliced_wasserstein(twin, truth)[0] == sm.sliced_wasserstein(pred, truth)[0]
    with scorer_build("main"):
        gap = sm.occupancy_dice(pred, truth)[0] - sm.occupancy_dice(twin, truth)[0]
    assert gap > 0.1


def test_parity_flips_with_row_order_under_prefix():
    truth, pred = helix_halves(11)
    sm = _sm()
    with scorer_build("main"):
        proper = sm.occupancy_dice(pred, truth)[0]
        mirror = sm.occupancy_dice(mirror_z(pred), truth)[0]
    assert proper - mirror > 0.1
    rng = np.random.default_rng(5)
    seen: dict[int, set] = {}
    for _ in range(16):
        P = pred[rng.permutation(len(pred))]
        parity = prefix_parity(P, len(truth))
        with scorer_build("main"):
            assert sm.occupancy_dice(P, truth)[0] == proper  # main is blind to row order
        with scorer_build("prefix"):
            value = sm.occupancy_dice(P, truth)[0]
        h = handedness(P, truth)
        assert h["local_parity"] == h["parity_pred"] == parity
        assert value == (proper if h["prefix_mode"] == "proper" else mirror)
        assert round(value, 4) == h["dice_prefix"]
        seen.setdefault(parity, set()).add(value)
    assert set(seen) == {-1, 1}  # row order alone flips the file's parity
    # the truth's rows (and parity) are fixed here, so one file parity gives proper mode and the other mirror mode
    assert {frozenset(v) for v in seen.values()} == {frozenset({proper}), frozenset({mirror})}


@pytest.mark.parametrize("metric_seed", [0, 7])
def test_prefix_value_follows_the_parity_rule_and_d2_scale_ignore_the_build(metric_seed):
    truth, pred = helix_halves(100 + metric_seed, n=1500)
    pred = pred[:1200]  # unequal counts: _match_n subsamples the truth
    sm = _sm()
    h = handedness(pred, truth, seed=metric_seed, live_n_true=999)
    vals = {}
    for build in SCORER_BUILDS:
        with scorer_build(build):
            vals[build] = (sm.d2_distance(pred, truth, seed=metric_seed), sm.scale_log_ratio(pred, truth),
                           sm.d2_distance(mirror_z(pred), truth, seed=metric_seed),
                           round(sm.occupancy_dice(pred, truth, seed=metric_seed)[0], 4))
    assert vals["main"][:3] == vals["prefix"][:3]  # d2_shape and scale_log_ratio do not depend on the build
    assert vals["main"][0] == vals["main"][2]  # d2_shape is blind to mirroring
    assert vals["main"][3] == h["dice_proper"]
    assert vals["prefix"][3] == h["dice_prefix"]
    assert h["local_parity_n_true"] == 999 and h["local_parity"] == prefix_parity(pred, 999, 0)


# ------------------------------------------------------------------------------------------------ score_board
# Recorded on 2026-10-19 from the emulator as it was before the build switch (git HEAD c511e4c, src/vec/evaluate.py),
# on the helix stages above with the synthetic T2 board of conftest.
GOLDEN = {
    "default": {
        "kw": {},
        "score": 93.15207098813904,
        "skills": {"de_score": 1.0, "de_direction": 1.0, "mmd_u": 0.9881697301525224, "variogram": 1.0,
                   "d2_shape": 1.0, "occupancy_dice": 0.20731707317073225, "scale_log_ratio": 1.0,
                   "neighborhood_mmd": 0.9974086437104708},
        "metrics": {"de_score": 0.6667, "de_direction": 0.7706, "mmd_u": 0.0011, "variogram": 0.002529,
                    "d2_shape": 0.00969, "occupancy_dice": 0.6857, "scale_log_ratio": -0.009,
                    "neighborhood_mmd": 0.00736, "sliced_wasserstein": 0.0462},
        "anchors": {"d2_shape": {"floor": 0.01316, "ceiling": 0.01304},
                    "occupancy_dice": {"floor": 0.7241, "ceiling": 0.7377},
                    "scale_log_ratio": {"floor": -0.2213, "ceiling": -0.0231}},
    },
    "frac05_seed3": {
        "kw": {"frac": 0.5, "seed": 3},
        "score": 90.85265414899031,
        "skills": {"de_score": 1.0, "de_direction": 1.0, "mmd_u": 0.9990605035700864, "variogram": 1.0,
                   "d2_shape": 0.6271929824561404, "occupancy_dice": 0.27681660899654076, "scale_log_ratio": 1.0,
                   "neighborhood_mmd": 1.0},
        "metrics": {"de_score": 0.8235, "de_direction": 0.7997, "mmd_u": -0.00014, "variogram": 0.000556,
                    "d2_shape": 0.00537, "occupancy_dice": 0.8182, "scale_log_ratio": -0.0048,
                    "neighborhood_mmd": 0.00088, "sliced_wasserstein": 0.03208},
        "anchors": {"d2_shape": {"floor": 0.00653, "ceiling": 0.00367},
                    "occupancy_dice": {"floor": 0.8471, "ceiling": 0.8391},
                    "scale_log_ratio": {"floor": -0.2319, "ceiling": -0.0075}},
    },
}


@pytest.mark.parametrize("case", sorted(GOLDEN))
def test_main_build_reproduces_the_current_emulator(helix_stages, case):
    ref, target, pred = helix_stages
    board = make_board("T2", True)
    g = GOLDEN[case]
    runs = [score_board(board, pred, target, ref, **g["kw"]), score_board(board, pred, target, ref, **g["kw"],
                                                                          build="main")]
    for r in runs:
        assert r.score == pytest.approx(g["score"], rel=0, abs=1e-12)
        assert {k: r.metrics[k] for k in g["metrics"]} == g["metrics"]
        assert {k: r.skills[k] for k in g["skills"]} == pytest.approx(g["skills"], rel=0, abs=1e-12)
        assert {m: r.anchors[m] for m in g["anchors"]} == g["anchors"]
        assert r.extras["build"] == "main"
    assert runs[0].metrics.keys() == runs[1].metrics.keys()
    assert all(same(runs[0].metrics[k], runs[1].metrics[k]) for k in runs[0].metrics)


def test_score_board_reports_both_modes_and_the_local_parity(helix_stages):
    ref, target, pred = helix_stages
    board = make_board("T2", True)
    frac = 0.25
    for seed in range(12):  # a seed whose emulated truth half puts this candidate in mirror mode under prefix
        truth, _ = split_half(subsample(target, frac, seed), seed=0)
        if handedness(pred.coords, truth.coords, seed=seed)["prefix_mode"] == "mirror":
            break
    else:
        pytest.fail("no mirror-mode seed among 12")
    res = {b: score_board(board, pred, target, ref, seed=seed, frac=frac, build=b) for b in SCORER_BUILDS}
    hm, hp = res["main"].extras["handedness"], res["prefix"].extras["handedness"]
    assert {k: v for k, v in hm.items() if k != "panel_dice_matches_build"} == {
        k: v for k, v in hp.items() if k != "panel_dice_matches_build"}  # the diagnostics do not depend on the build
    assert hm["prefix_mode"] == "mirror" and hm["dice_proper"] - hm["dice_mirror"] > 0.05
    assert res["main"].metrics["occupancy_dice"] == hm["dice_proper"]
    assert res["prefix"].metrics["occupancy_dice"] == hp["dice_prefix"] == hp["dice_mirror"]
    assert hm["panel_dice_matches_build"] and hp["panel_dice_matches_build"]
    assert hm["local_parity_n_true"] == res["main"].extras["truth_cells"]  # synthetic board: no truth_cells
    assert hm["local_parity"] == prefix_parity(pred.coords, hm["local_parity_n_true"], 0)
    assert set(res["main"].extras["anchor_handedness"]) == {"floor", "ceiling"}
    for kind, h in res["prefix"].extras["anchor_handedness"].items():  # computed anchors follow the build too
        assert res["prefix"].anchors["occupancy_dice"][kind] == h["dice_prefix"]
        assert res["main"].anchors["occupancy_dice"][kind] == h["dice_proper"]
    for k, v in res["main"].metrics.items():
        if k not in CANONICALISING:
            assert same(v, res["prefix"].metrics[k]), k

    # with a board's truth_cells the local parity is computed for that count; published anchors skip floor/ceiling
    pub = Board(key="T2:synthetic", task="T2", genes=GENES, needs_coords=True, truth_cells=321,
                anchors=res["main"].anchors)
    summary = score_board_seeds(pub, pred, target, ref, seeds=(seed, seed + 1), frac=frac, anchors="published",
                                build="prefix")
    assert summary["build"] == "prefix"
    assert summary["handedness"]["local_parity_n_true"] == 321
    assert summary["handedness"]["local_parity"] == [prefix_parity(pred.coords, 321, 0)]
    assert all("anchor_handedness" not in r["extras"] for r in summary["per_seed"])


def test_score_board_mirror_twin_end_to_end(helix_stages):
    ref, target, pred = helix_stages
    twin = Stage(X=pred.X, coords=mirror_z(pred.coords), celltype=pred.celltype)
    board = make_board("T2", True)
    kw = {"seed": 1, "frac": 0.25}
    pre = [score_board(board, p, target, ref, build="prefix", **kw) for p in (pred, twin)]
    main = [score_board(board, p, target, ref, build="main", **kw) for p in (pred, twin)]
    assert pre[0].score == pre[1].score
    assert all(same(pre[0].metrics[k], pre[1].metrics[k]) for k in pre[0].metrics)
    assert main[0].metrics["occupancy_dice"] - main[1].metrics["occupancy_dice"] > 0.05
    assert main[0].score > main[1].score
    for k in main[0].metrics:
        if k not in CANONICALISING:
            assert same(main[0].metrics[k], main[1].metrics[k]), k


# ------------------------------------------------------------------------------------------------ anchor check
@pytest.fixture(scope="module")
def ac():
    spec = importlib.util.spec_from_file_location("scorer_build_anchor_check",
                                                  REPO / "experiments" / "scorer_build" / "anchor_check.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_anchor_check_exits_cleanly_without_answers(ac, tmp_path, capsys):
    (tmp_path / "T2").mkdir()
    for name in ("E7.25.h5ad", "E8.25_late.h5ad", "E9.5.h5ad"):  # released training stages are not answers
        (tmp_path / "T2" / name).touch()
    out = tmp_path / "out"
    rc = ac.main(["--data", str(tmp_path), "--index", str(tmp_path / "missing.json"), "--out", str(out)])
    assert rc == 0
    assert "No T2 validation answers" in capsys.readouterr().out
    assert not out.exists()


def test_anchor_check_finds_answer_files_by_stage(ac, tmp_path):
    t2 = tmp_path / "T2"
    t2.mkdir()
    for name in ("E7.25.h5ad", "E8.25_late.h5ad", "E9.5.h5ad", "E10.5.h5ad"):
        (t2 / name).touch()
    (tmp_path / "T1").mkdir()
    (tmp_path / "T1" / "E10.5.h5ad").touch()  # the T1 answer is not a T2 board's
    assert ac.find_answer(tmp_path, "E7.5") is None
    assert ac.find_answer(tmp_path, "E8.5") is None
    assert ac.find_answer(tmp_path, "E10.5") == t2 / "E10.5.h5ad"


def test_anchor_check_guards(ac, tmp_path):
    val = {"split": "val_extrap", "truth_cells": 10, "ref_cells": 10,
           "anchors": {m: {"floor": 0.0, "ceiling": 1.0} for m in ac.SHAPE_METRICS}}
    good = {"T2:heart:val_extrap": val, "T2:heart:test": {"split": "test", "truth_cells": 10}}
    path = tmp_path / "index.json"
    ac.screen_index(good, ["T2:heart:val_extrap"], path)  # a stripped repo index passes
    leaky = {**good, "T2:heart:test": {"split": "test", "anchors": {"d2_shape": {"floor": 0.0, "ceiling": 1.0}}}}
    with pytest.raises(SystemExit, match="non-validation"):
        ac.screen_index(leaky, ["T2:heart:val_extrap"], path)
    with pytest.raises(SystemExit, match="non-validation"):  # a validation key with a test split
        ac.screen_index({"T2:heart:val_extrap": {**val, "split": "test"}}, ["T2:heart:val_extrap"], path)
    for key in ("T2:heart:test", "T1:val", "T3:gata4"):
        with pytest.raises(SystemExit, match="not a T2 validation board"):
            ac.screen_index(good, [key], path)
    panels = tmp_path / "panels"
    panels.mkdir()
    (panels / "index.json").write_text(json.dumps(good))
    with pytest.raises(SystemExit, match="panels"):
        ac.load_index_entries(panels / "index.json", ["T2:heart:val_extrap"])
    path.write_text(json.dumps(good))
    entries = ac.load_index_entries(path, ["T2:heart:val_extrap"])
    assert set(entries) == {"T2:heart:val_extrap"} and set(entries["T2:heart:val_extrap"]["anchors"]) == set(
        ac.SHAPE_METRICS)


def test_anchor_check_classification_rules(ac):
    rows = [{"dice_main": 0.95, "dice_mirror": 0.82}]
    assert ac.classify_exact(0.95, rows) == "proper"
    assert ac.classify_exact(0.82, rows) == "mirror"
    assert ac.classify_exact(0.90, rows) == "neither"
    assert ac.classify_exact(0.95, [{"dice_main": 0.95, "dice_mirror": 0.95}]) == "tie"
    assert ac.classify_exact(0.95, rows + [{"dice_main": 0.80, "dice_mirror": 0.95}]) == "ambiguous"

    def stats(zp, zm, zd=0.0, zs=0.0):
        return {"dice_main": {"z_of_published": zp}, "dice_mirror": {"z_of_published": zm},
                "d2_shape": {"z_of_published": zd}, "scale_log_ratio": {"z_of_published": zs}}

    assert ac.classify_distribution(0.9, stats(-0.5, 12.0)) == "proper"
    assert ac.classify_distribution(0.9, stats(-9.0, 1.0)) == "mirror"
    assert ac.classify_distribution(0.9, stats(1.0, 1.5)) == "undetermined"
    assert ac.classify_distribution(0.9, stats(0.5, 12.0, zd=6.0)) == "undetermined (fingerprints off)"


def test_anchor_check_recovers_a_baked_protocol_and_its_dice_modes(ac):
    T, _ = helix_coords(20000, np.random.default_rng(1))
    R, _ = helix_coords(15000, np.random.default_rng(2), scale=0.8)
    res = ac.self_test(T, R, seeds=(7,), n_seeds=4, cases=ac.SELF_TEST_CASES[:1])
    (r,) = res.values()
    st = r["self_test"]
    assert st["passed"], st
    assert st["n_hits"] == 1
    assert set(st["true_modes"].values()) == {"proper", "mirror"}  # this bake has one anchor in each mode
    t = r["theory"]
    assert t["prefix_equals_parity_rule"] == t["mirror_twin_unchanged_under_prefix"] == t["rows"]
    assert t["d2_scale_build_independent"] == t["d2_scale_rows_under_both_builds"] == 2
