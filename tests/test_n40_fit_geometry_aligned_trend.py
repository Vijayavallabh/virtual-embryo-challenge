"""Reading N40 (experiments/harness/n40_fit_geometry_aligned_trend.py): reading N39's fits for another arm and for two
bases, on stage days and synthetic score tables built here. Nothing is read from data/ or from a committed run, and
nothing is scored."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"
# the nine practice tasks by name: (before file, after file, target file)
NINE = {"half": ("T2/E7.5.h5ad", "T2/E8.5.h5ad", "T2/E8.0.h5ad"),
        "late": ("T2/E7.5.h5ad", "T2/E8.25_late.h5ad", "T2/E8.0.h5ad"),
        "new1": ("T2/E7.5.h5ad", "T2/E8.5.h5ad", "T2/E8.25_late.h5ad"),
        "new2": ("T2/E7.5.h5ad", "T2/E8.75.h5ad", "T2/E8.5.h5ad"),
        "x1": ("T2/E8.0.h5ad", "T2/E8.75.h5ad", "T2/E8.25_late.h5ad"),
        "x2": ("T2/E8.0.h5ad", "T2/E8.75.h5ad", "T2/E8.5.h5ad"),
        "x3": ("T2/E8.25_late.h5ad", "T2/E8.75.h5ad", "T2/E8.5.h5ad"),
        "e75": ("T2/E7.25.h5ad", "T2/E8.0.h5ad", "T2/E7.5.h5ad"),
        "e725": ("T2/E6.75.h5ad", "T2/E7.5.h5ad", "T2/E7.25.h5ad")}


@pytest.fixture(scope="module")
def n40():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n40_fit_geometry_aligned_trend",
                                                  HARNESS / "n40_fit_geometry_aligned_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tasks():
    return {n: SimpleNamespace(key=f"T2e/{n}", inputs={"before": b, "after": a}, target=t)
            for n, (b, a, t) in NINE.items()}


def geo(n40):
    return {n: n40.n39.geometry(*(n40.harness.STAGE_DAYS[f] for f in (b, a, t)), n40.n39.LIKE[n])
            for n, (b, a, t) in NINE.items()}


def tables(n40, own: dict[str, float], aligned: dict[str, float] | None = None):
    """Per-seed scores, constant over seeds: fgw at 60, fgw_aligned at 60 + aligned[task] (default 0), the arm at
    60 + own[task]."""
    aligned = aligned or {}
    return {f"T2e/{n}": {"fgw": np.full(10, 60.0), "fgw_aligned": np.full(10, 60.0 + aligned.get(n, 0.0)),
                         n40.ARM: np.full(10, 60.0 + own[n])} for n in NINE}


def test_the_constants_are_the_rows(n40):
    assert n40.ARM == "fgw_aligned_trend" and n40.RULE_QUANTITY == "own"
    assert n40.QUANTITIES == {"own": "fgw", "trend": "fgw_aligned"} and n40.PROMPTED == "e725"


def test_a_planted_line_is_recovered_and_the_rule_reads_the_gain_over_fgw(n40):
    g = geo(n40)
    own = {n: (1.0 + 10.0 * g[n]["g"] if n40.n39.LIKE[n] else 0.5) for n in NINE}  # +1 + 10 g on the classed tasks
    res = n40.check(tables(n40, own), tasks())
    v = res["quantities"]["own"]
    f = v["fits"]["all_seven"]
    assert f["n"] == 7 and f["a"] == pytest.approx(1.0) and f["b"] == pytest.approx(10.0)
    p = v["prediction_at_the_test"]["all_seven"]
    assert p["g_early"] == pytest.approx(0.3077 - 0.5, abs=1e-4) and p["y_early"] == pytest.approx(-0.923, abs=1e-3)
    assert p["y_late"] == pytest.approx(1.0 + 10.0 * (0.5385 - 0.5), abs=1e-3)
    assert res["not_put_forward"]
    assert v["tasks"]["e75"]["like"] is None and v["tasks"]["e725"]["ahead"] == 10


def test_a_gain_that_does_not_follow_the_geometry_is_not_withdrawn(n40):
    res = n40.check(tables(n40, {n: 2.0 for n in NINE}), tasks())
    v = res["quantities"]["own"]
    assert v["fits"]["all_seven"]["b"] == pytest.approx(0.0, abs=1e-9)
    assert v["prediction_at_the_test"]["without_e725"]["y_early"] == pytest.approx(2.0)
    assert not res["not_put_forward"]


def test_either_fit_at_or_below_zero_withdraws(n40):
    g = geo(n40)
    own = {n: (1.0 + 10.0 * g[n]["g"] if n40.n39.LIKE[n] else 0.0) for n in NINE}
    own["e725"] = -40.0  # the task with the largest g far below the line turns the full fit's slope around
    res = n40.check(tables(n40, own), tasks())
    p = res["quantities"]["own"]["prediction_at_the_test"]
    assert p["all_seven"]["y_early"] > 0 and p["without_e725"]["y_early"] < 0 and res["not_put_forward"]
    own = {n: (3.0 + 10.0 * g[n]["g"] if n40.n39.LIKE[n] else 0.0) for n in NINE}  # +1.08 at the earlier geometry
    own["e725"] = 60.0  # the same task far above the line steepens the full fit: it alone falls below zero there
    res = n40.check(tables(n40, own), tasks())
    p = res["quantities"]["own"]["prediction_at_the_test"]
    assert p["all_seven"]["y_early"] < 0 and p["without_e725"]["y_early"] > 0 and res["not_put_forward"]


def test_the_trend_part_is_measured_against_the_aligned_arm_and_decides_nothing(n40):
    g = geo(n40)
    aligned = {n: 3.0 for n in NINE}  # the aligned positions alone gain 3 everywhere
    own = {n: 3.0 + (10.0 * g[n]["g"] if n40.n39.LIKE[n] else 0.0) for n in NINE}  # the trend adds 10 g on top
    res = n40.check(tables(n40, own, aligned), tasks())
    trend = res["quantities"]["trend"]
    assert trend["base"] == "fgw_aligned"
    assert trend["fits"]["all_seven"]["a"] == pytest.approx(0.0, abs=1e-9)
    assert trend["fits"]["all_seven"]["b"] == pytest.approx(10.0)
    assert trend["prediction_at_the_test"]["all_seven"]["y_early"] < 0
    own_fit = res["quantities"]["own"]["prediction_at_the_test"]["all_seven"]
    assert own_fit["y_early"] == pytest.approx(3.0 - 1.923, abs=1e-3) and not res["not_put_forward"]


def test_the_reading_is_void_without_an_arm_and_the_text_states_the_outcome(n40, monkeypatch):
    tk = tasks()
    scores = tables(n40, {n: 1.0 for n in NINE})
    monkeypatch.setattr(n40.nt, "all_tasks", lambda: tk)
    monkeypatch.setattr(n40.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n40.reading([])
    assert not res["void"] and not res["not_put_forward"]
    txt = n40.text(res)
    assert "stays on record as an alternative" in txt and "quantity 'trend'" in txt and "truth like E7.5" in txt
    del scores["T2e/x3"]["fgw_aligned"]
    res = n40.reading([])
    assert res["void"] and any("no scores for fgw_aligned" in w for w in res["why"])
    assert n40.text(res).startswith("VOID")
    g = geo(n40)
    scores = tables(n40, {n: (1.0 + 10.0 * g[n]["g"] if n40.n39.LIKE[n] else 0.0) for n in NINE})
    monkeypatch.setattr(n40.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    assert "NOT put forward" in n40.text(n40.reading([]))
