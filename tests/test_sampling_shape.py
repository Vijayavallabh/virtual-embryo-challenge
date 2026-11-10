"""experiments/harness/sampling_shape.py (tries E16 and E17 of experiments/harness/TRIES.md) on synthetic rows."""

from __future__ import annotations

import importlib.util
import sys

import pytest

from vec import harness
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "harness" / "sampling_shape.py"


@pytest.fixture(scope="module")
def ss():
    spec = importlib.util.spec_from_file_location("sampling_shape_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _skills(shape: float) -> dict:
    return {"de_score": 0.6, "de_direction": 0.7, "mmd_u": 0.65, "variogram": 0.6, "d2_shape": shape,
            "occupancy_dice": shape, "scale_log_ratio": 0.5, "neighborhood_mmd": 0.6}


def _result(ss, key, fgw, s0, s1, commit="abc", clean=True):
    rows = []
    for k, seed in enumerate(SEEDS):
        sc = {"floor": 50.0, "fgw": fgw[k], "fgw_shape0": s0[k], "fgw_shape1": s1[k]}
        rows.append({"seed": seed, "scores": sc, "skills": {a: _skills(0.5) for a in ss.ARMS}})
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean}}


def _pair(ss, main_s1, control_s1=None, **kw):
    n = len(SEEDS)
    return {ss.MAIN.key: _result(ss, ss.MAIN.key, [60.0] * n, [59.0] * n, main_s1, **kw),
            ss.CONTROL.key: _result(ss, ss.CONTROL.key, [60.0] * n, [59.0] * n, control_s1 or [59.0] * n, **kw)}


def test_the_tasks_pass_the_harness_checks_and_are_never_counted(ss):
    for t in ss.TASKS:
        assert harness.check_task(t) is t and not t.counted
        assert t.test_board == "T2:embryo:test_interp" and t.panel_board == "T2:embryo:val_interp"
    assert ss.MAIN.target == "T2/E8.25_late.h5ad" and ss.MAIN.inputs == {"before": "T2/E8.0.h5ad",
                                                                       "after": "T2/E8.5.h5ad"}
    for name in ss.ARMS[1:]:
        assert all(harness.CANDIDATES[name].applies(t) for t in ss.TASKS)


def test_the_rule_recommends_the_switch_only_above_the_bar_in_eight_seeds(ss):
    n = len(SEEDS)
    res = ss.reading(_pair(ss, [61.0] * n))  # D = 61 - 60 = +1 in every seed
    assert not res["void"] and res["passes"] and res["main"]["D_ahead"] == n
    assert res["main"]["D_mean"] == pytest.approx(1.0)
    assert ss.reading(_pair(ss, [60.2] * n))["passes"] is False  # +0.2: below the bar
    seven = [62.0] * 7 + [59.0] * (n - 7)  # mean +1.4 but ahead in 7 of 10
    res = ss.reading(_pair(ss, seven))
    assert res["main"]["D_mean"] > ss.BAR and res["main"]["D_ahead"] == 7 and not res["passes"]


def test_the_control_decides_nothing(ss):
    n = len(SEEDS)
    a = ss.reading(_pair(ss, [61.0] * n, control_s1=[59.0] * n))
    b = ss.reading(_pair(ss, [61.0] * n, control_s1=[70.0] * n))
    assert a["passes"] and b["passes"] and a["recommendation"] == b["recommendation"]


def test_the_reading_is_void_on_missing_or_mixed_results(ss):
    n = len(SEEDS)
    good = _pair(ss, [61.0] * n)
    assert ss.reading({ss.MAIN.key: good[ss.MAIN.key]})["void"]
    mixed = dict(good)
    mixed[ss.CONTROL.key] = _result(ss, ss.CONTROL.key, [60.0] * n, [59.0] * n, [59.0] * n, commit="other")
    assert ss.reading(mixed)["void"]
    assert ss.reading(_pair(ss, [61.0] * n, clean=False))["void"]


def test_the_reading_is_void_when_a_shape_variant_changes_expression(ss):
    n = len(SEEDS)
    res = _pair(ss, [61.0] * n)
    res[ss.MAIN.key]["rows"][3]["skills"]["fgw_shape1"]["mmd_u"] = 0.7
    out = ss.reading(res)
    assert out["void"] and not out["passes"] and any("mmd_u" in w for w in out["why"])
