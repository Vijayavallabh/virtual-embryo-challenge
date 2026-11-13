"""Reading N39 (experiments/harness/n39_fit_geometry_check.py): the trend line's weights against the mixture's
shares, on stage days and synthetic score tables built here. Nothing is read from data/ or from a committed run, and
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
def n39():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n39_fit_geometry_check", HARNESS / "n39_fit_geometry_check.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tasks():
    return {n: SimpleNamespace(key=f"T2e/{n}", inputs={"before": b, "after": a}, target=t)
            for n, (b, a, t) in NINE.items()}


def tables(n39, own: dict[str, float], second: float = -5.0):
    """Per-seed scores, constant over seeds: fgw at 60, the arm at 60 + own[task], the second upload at 60 + second."""
    return {f"T2e/{n}": {n39.BASE: np.full(10, 60.0), n39.ARM: np.full(10, 60.0 + own[n]),
                         n39.SECOND: np.full(10, 60.0 + second)} for n in NINE}


def test_the_lines_weights_sum_to_one_and_give_the_fitted_value(n39):
    days = [6.75, 7.25, 7.5, 8.0]
    w = n39.line_weights(days, 7.75)
    assert w.sum() == pytest.approx(1.0)
    y = np.array([0.3, -1.2, 2.0, 0.7])
    b, a = np.polyfit(days, y, 1)
    assert float(w @ y) == pytest.approx(a + b * 7.75)
    with pytest.raises(ValueError, match="one day"):
        n39.line_weights([7.5, 7.5], 7.5)


def test_the_sample_rule_is_row_e35s(n39):
    assert n39.sample_days(7.5, 8.5) == [6.75, 7.25, 7.5, 8.5]  # half and new1: E8.0 lies inside and is left out
    assert n39.sample_days(7.5, 8.3) == [6.75, 7.25, 7.5, 8.3]  # late
    assert n39.sample_days(7.5, 8.75) == [6.75, 7.25, 7.5, 8.75]  # new2
    assert n39.sample_days(8.0, 8.75) == [6.75, 7.25, 7.5, 8.0, 8.75]  # x1 and x2
    assert n39.sample_days(8.3, 8.75) == [6.75, 7.25, 7.5, 8.0, 8.3, 8.75]  # x3
    assert n39.sample_days(7.25, 8.0) == [6.75, 7.25, 8.0]  # the E7.5 backtest: its configured inputs
    assert n39.sample_days(6.75, 7.5) == [6.75, 7.5, 8.0]  # the E7.25 backtest: its configured inputs
    assert n39.sample_days(7.5, 8.0) == [6.75, 7.25, 7.5, 8.0]  # the test


def test_the_geometry_at_the_test_and_on_the_two_tasks_that_set_the_pattern(n39):
    early = n39.geometry(7.5, 8.0, 7.75, "before")
    late = n39.geometry(7.5, 8.0, 7.75, "after")
    assert early["t"] == pytest.approx(0.5)
    assert early["weight_before"] == pytest.approx(0.25 + 0.375 * 0.125 / 0.8125)  # 0.308 on E7.5
    assert early["weight_after"] == pytest.approx(0.25 + 0.375 * 0.625 / 0.8125)  # 0.538 on E8.0
    assert early["g"] == pytest.approx(0.3077 - 0.5, abs=1e-4) and late["g"] == pytest.approx(0.5385 - 0.5, abs=1e-4)
    e725 = n39.geometry(6.75, 7.5, 7.25, "before")  # the line gives the truth-like sample more than the mixture
    assert e725["share_like"] == pytest.approx(1 / 3) and e725["g"] == pytest.approx(0.4737 - 1 / 3, abs=1e-4)
    x1 = n39.geometry(8.0, 8.75, 8.3, "before")  # and here less
    assert x1["share_like"] == pytest.approx(0.6) and x1["g"] == pytest.approx(0.2978 - 0.6, abs=1e-3)
    own = n39.geometry(7.25, 8.0, 7.5, None)
    assert "g" not in own and own["like"] is None


def test_a_planted_line_is_recovered_and_the_advice_follows_the_prediction(n39):
    tk = tasks()
    geo = {n: n39.geometry(*(n39.harness.STAGE_DAYS[f] for f in (b, a, t)), n39.LIKE[n])
           for n, (b, a, t) in NINE.items()}
    own = {n: (1.0 + 10.0 * geo[n]["g"] if n39.LIKE[n] else 0.5) for n in NINE}  # y = 1 + 10 g on the classed tasks
    res = n39.check(tables(n39, own), tk)
    f = res["fit_all_seven"]
    assert f["n"] == 7 and f["a"] == pytest.approx(1.0) and f["b"] == pytest.approx(10.0)
    assert f["r"] == pytest.approx(1.0)
    p = res["prediction_at_the_test"]["all_seven"]
    assert p["y_early"] == pytest.approx(1.0 + 10.0 * (0.3077 - 0.5), abs=1e-3) and p["y_early"] < 0
    assert p["y_late"] == pytest.approx(1.0 + 10.0 * (0.5385 - 0.5), abs=1e-3)
    assert res["advice_changes_to_no"]
    assert res["tasks"]["e75"]["like"] is None and res["tasks"]["e725"]["ahead"] == 10
    flat = {n: 2.0 for n in NINE}  # the same gain whatever the geometry: the prediction at the test stays positive
    res = n39.check(tables(n39, flat), tk)
    assert res["fit_all_seven"]["b"] == pytest.approx(0.0, abs=1e-9) and not res["advice_changes_to_no"]
    assert np.isnan(res["fit_all_seven"]["r"])  # no spread in the gain: no correlation is reported
    assert res["prediction_at_the_test"]["without_e725"]["y_early"] == pytest.approx(2.0)


def test_either_fit_at_or_below_zero_changes_the_advice(n39):
    tk = tasks()
    geo = {n: n39.geometry(*(n39.harness.STAGE_DAYS[f] for f in (b, a, t)), n39.LIKE[n])
           for n, (b, a, t) in NINE.items()}
    own = {n: (3.0 + 10.0 * geo[n]["g"] if n39.LIKE[n] else 0.0) for n in NINE}  # +1.08 at the test's early world
    own["e725"] = 9.0  # one task far above the line pulls the full fit up; without it the line is the planted one
    res = n39.check(tables(n39, own), tk)
    p = res["prediction_at_the_test"]
    assert p["without_e725"]["y_early"] == pytest.approx(3.0 - 1.923, abs=1e-3) and p["all_seven"]["y_early"] > 0
    assert not res["advice_changes_to_no"]
    own = {n: (1.0 + 10.0 * geo[n]["g"] if n39.LIKE[n] else 0.0) for n in NINE}
    own["e725"] = 30.0  # the full fit is carried well above zero by one task; the fit without it is not
    res = n39.check(tables(n39, own), tk)
    p = res["prediction_at_the_test"]
    assert p["without_e725"]["y_early"] < 0 and res["advice_changes_to_no"]


def test_the_pair_value_and_the_d_zero_figures(n39):
    tk = tasks()
    own = {n: 0.0 for n in NINE}
    own.update(new2=0.4, e75=-0.2, e725=3.6)
    res = n39.check(tables(n39, own), tk)  # the second upload trails everywhere: the pair value is the own gain
    assert res["tasks"]["e725"]["pair_value"] == pytest.approx(3.6) and res["tasks"]["half"]["pair_value"] == 0.0
    z = res["n37_with_D_zero"]
    assert z["S0.45"] == pytest.approx(0.55 * (0.4 - 0.2 + 3.6) / 3)
    assert z["S0.45_without_e725"] == pytest.approx(0.55 * (0.4 - 0.2) / 2)


def test_the_reading_is_void_without_the_arm_and_the_text_states_the_advice(n39, monkeypatch):
    tk = tasks()
    scores = tables(n39, {n: 1.0 for n in NINE})
    monkeypatch.setattr(n39.nt, "all_tasks", lambda: tk)
    monkeypatch.setattr(n39.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n39.reading([])
    assert not res["void"] and not res["advice_changes_to_no"]
    txt = n39.text(res)
    assert "advice on request R-D: stands" in txt and "truth like E7.5" in txt and "own" in txt
    del scores["T2e/x3"][n39.ARM]
    res = n39.reading([])
    assert res["void"] and any("no scores for fgw_trend_render" in w for w in res["why"])
    assert n39.text(res).startswith("VOID")
