"""Reading N38 (experiments/harness/n38_second_beside_trend.py): the second-upload pair rule with the first upload
changed, on synthetic score tables built here. Nothing is read from data/ or from a committed run, and nothing is
scored."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"


@pytest.fixture(scope="module")
def n38():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n38_second_beside_trend", HARNESS / "n38_second_beside_trend.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n38, own: dict[str, dict[str, float | None]], first: dict[str, float] | None = None,
           base: dict[str, float] | None = None):
    """Per-seed scores on the nine tasks, constant over seeds: the built first upload fgw at 60; the new first upload
    at 60 + first[task] (default 0); the second upload as built at 60 + base[task] (default: 3 on the DEEP tasks, -2
    on the others); each candidate second upload at 60 + own[arm][task] (default 0; None = no scores there)."""
    nt = n38.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    base = {**{n: 3.0 for n in nt.DEEP}, **{n: -2.0 for n in nt.OTHER}, **(base or {})}
    first = {**{n: 0.0 for n in names}, **(first or {})}
    scores = {}
    for n in names:
        s = {n38.BUILT_FIRST: np.full(10, 60.0), n38.FIRST: np.full(10, 60.0 + first[n]),
             n38.BASE_SECOND: np.full(10, 60.0 + base[n])}
        for arm, by_task in own.items():
            if arm in s:
                continue
            v = by_task.get(n, 0.0)
            if v is not None:
                s[arm] = np.full(10, 60.0 + v)
        scores[tasks[n].key] = s
    return scores, tasks


def test_the_constants_are_the_rows_and_the_nine_tasks_own(n38):
    assert n38.FIRST == "fgw_trend_render" and n38.BASE_SECOND == n38.nt.RA == "fgw_aligned_series"
    assert n38.BUILT_FIRST == "fgw" and n38.LEAVE_OUT_P == 0.45 and n38.nt.BAR == 0.3
    assert "fgw" in n38.SECONDS and n38.BASE_SECOND not in n38.SECONDS and n38.FIRST not in n38.SECONDS


def test_a_second_upload_that_only_matches_the_one_as_built_is_not_advised(n38):
    same = {**{n: 3.0 for n in n38.nt.DEEP}, **{n: -2.0 for n in n38.nt.OTHER}}
    scores, tasks = tables(n38, {"fgw_series": same})
    v = n38.price(scores, tasks, candidates=("fgw_series",))
    r = v["second_uploads"]["fgw_series"]
    assert r["readable"] and r["S0.45"] == pytest.approx(0.0) and not r["passes"]
    assert v["advised_second"] == n38.BASE_SECOND


def test_a_second_upload_ahead_on_every_side_is_advised_and_the_largest_wins(n38):
    # DEEP: 1 above the second upload as built; OTHER: 1 above the first upload
    a = {**{n: 4.0 for n in n38.nt.DEEP}, **{n: 1.0 for n in n38.nt.OTHER}}
    b = {**{n: 5.0 for n in n38.nt.DEEP}, **{n: 1.0 for n in n38.nt.OTHER}}
    scores, tasks = tables(n38, {"fgw_series": a, "fgw_shape0_series": b})
    v = n38.price(scores, tasks, candidates=("fgw_series", "fgw_shape0_series"))
    ra, rb = v["second_uploads"]["fgw_series"], v["second_uploads"]["fgw_shape0_series"]
    assert ra["D"] == pytest.approx(1.0) and ra["O"] == pytest.approx(1.0) and ra["passes"]
    assert rb["D"] == pytest.approx(2.0) and rb["passes"] and v["advised_second"] == "fgw_shape0_series"


def test_giving_up_the_deeper_side_for_a_small_gain_elsewhere_fails(n38):
    # the arm loses the base's lead of 3 on every DEEP task and gains 1 on every OTHER task
    arm = {**{n: 0.0 for n in n38.nt.DEEP}, **{n: 1.0 for n in n38.nt.OTHER}}
    scores, tasks = tables(n38, {"fgw_early_specialist": arm})
    r = n38.price(scores, tasks, candidates=("fgw_early_specialist",))["second_uploads"]["fgw_early_specialist"]
    assert r["D"] == pytest.approx(-3.0) and r["O"] == pytest.approx(1.0)
    assert r["S0.45"] == pytest.approx(-0.8) and r["S0.3"] == pytest.approx(-0.2) and r["S0.6"] == pytest.approx(-1.4)
    assert not r["rule"]["a_S0.45_above_bar"] and not r["rule"]["b_S0.3_and_S0.6_not_below_0"] and not r["passes"]


def test_a_lower_worst_of_nine_fails_clause_c(n38):
    # the base second upload carries one weak task (x1) that the arm drops, while the arm gains everywhere else
    base = {"x1": -10.0}
    first = {"x1": -15.0}
    arm = {**{n: 5.0 for n in n38.nt.DEEP}, **{n: 2.0 for n in n38.nt.OTHER}, "x1": -14.0}
    scores, tasks = tables(n38, {"fgw_series": arm}, first=first, base=base)
    r = n38.price(scores, tasks, candidates=("fgw_series",))["second_uploads"]["fgw_series"]
    assert r["tasks"]["x1"] == pytest.approx(-4.0) and r["worst9_abs"] == pytest.approx(46.0)
    assert r["rule"]["a_S0.45_above_bar"] and not r["rule"]["c_worst9_not_lower"] and not r["passes"]


def test_the_value_is_measured_beside_the_new_first_upload_not_beside_fgw(n38):
    first = {n: 5.0 for n in n38.nt.OTHER}  # the new first upload already leads the OTHER tasks by 5
    arm = {**{n: 3.0 for n in n38.nt.DEEP}, **{n: 4.0 for n in n38.nt.OTHER}}  # +4 over fgw there, but under the first
    scores, tasks = tables(n38, {"fgw_series": arm}, first=first)
    v = n38.price(scores, tasks, candidates=("fgw_series",))
    r = v["second_uploads"]["fgw_series"]
    assert r["O"] == pytest.approx(0.0) and r["D"] == pytest.approx(0.0) and not r["passes"]
    assert v["base_pair_against_the_pair_as_built"]["e725"] == pytest.approx(5.0)
    assert v["base_pair_against_the_pair_as_built"]["half"] == pytest.approx(0.0)


def test_one_task_carrying_the_gain_fails_the_leave_one_out_clause(n38):
    arm = {**{n: 3.0 for n in n38.nt.DEEP}, **{n: -2.0 for n in n38.nt.OTHER}, "e725": 4.0}
    scores, tasks = tables(n38, {"fgw_series": arm})
    r = n38.price(scores, tasks, candidates=("fgw_series",))["second_uploads"]["fgw_series"]
    assert r["rule"]["a_S0.45_above_bar"] and not r["rule"]["d_leave_one_task_out_positive"] and not r["passes"]


def test_an_arm_missing_on_a_task_is_not_readable(n38):
    arm = {**{n: 5.0 for n in n38.nt.DEEP}, **{n: 5.0 for n in n38.nt.OTHER}, "x3": None}
    scores, tasks = tables(n38, {"fgw_aligned_hi": arm})
    v = n38.price(scores, tasks, candidates=("fgw_aligned_hi",))
    r = v["second_uploads"]["fgw_aligned_hi"]
    assert not r["readable"] and r["missing_on"] == ["x3"] and not r["passes"]
    assert v["advised_second"] == n38.BASE_SECOND


def test_the_reading_is_void_without_the_first_upload_and_the_text_names_the_advice(n38, monkeypatch):
    scores, tasks = tables(n38, {"fgw_series": {}})
    monkeypatch.setattr(n38.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n38.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n38.reading([])
    assert not res["void"] and res["advised_second"] == n38.BASE_SECOND
    txt = n38.text(res)
    assert "advised second upload: fgw_aligned_series" in txt and "not readable" in txt and "fgw_series" in txt
    del scores[tasks["half"].key][n38.FIRST]
    res = n38.reading([])
    assert res["void"] and any("no scores for fgw_trend_render" in w for w in res["why"])
    assert n38.text(res).startswith("VOID")
