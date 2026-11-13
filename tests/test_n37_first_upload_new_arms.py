"""Reading N37 (experiments/harness/n37_first_upload_new_arms.py): row N33's first-upload rule applied to the new
embryo arms, on synthetic score tables built here. Nothing is read from data/ or from a committed run, and nothing
is scored."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"


@pytest.fixture(scope="module")
def n37():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n37_first_upload_new_arms", HARNESS / "n37_first_upload_new_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n37, own: dict[str, dict[str, float]], second: dict[str, float] | None = None):
    """Per-seed scores on the nine tasks: fgw at 60 on every task and seed, the second upload at 60 + second[task]
    (default: 3 above on the DEEP tasks, 2 below on the others), and each arm at 60 + own[arm][task] (default 0)."""
    nt = n37.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    second = {**{n: 3.0 for n in nt.DEEP}, **{n: -2.0 for n in nt.OTHER}, **(second or {})}
    scores = {}
    for n in names:
        s = {n37.BASE: np.full(10, 60.0), n37.SECOND: np.full(10, 60.0 + second[n])}
        for arm, by_task in own.items():
            if by_task.get(n) is not None:
                s[arm] = np.full(10, 60.0 + by_task.get(n, 0.0))
        scores[tasks[n].key] = s
    return scores, tasks


def flat(n37, **over):
    return {**{n: 0.0 for n in [*n37.nt.DEEP, *n37.nt.OTHER]}, **over}


def test_the_rule_is_n33s_and_its_constants_are_the_nine_tasks_own(n37):
    assert n37.SECOND == n37.nt.RA == "fgw_aligned_series" and n37.BASE == "fgw"
    assert n37.LEAVE_OUT_P == 0.45 and n37.nt.BAR == 0.3 and n37.nt.PS == (0.3, 0.45, 0.6)
    assert n37.NEW == ("fgw_trend_render", "fgw_aligned_trend", "fgw_aligned_hi", "fgw_early_specialist")


def test_an_arm_that_gains_where_the_second_upload_does_not_lead_passes(n37):
    arm = flat(n37, new2=2.0, e75=2.0, e725=2.0)  # +2 on each OTHER task, where the second upload trails fgw
    scores, tasks = tables(n37, {"fgw_trend_render": arm})
    v = n37.price(scores, tasks, candidates=("fgw_trend_render",))
    r = v["first_uploads"]["fgw_trend_render"]
    assert r["readable"] and r["passes"] and v["requested_first_upload"] == "fgw_trend_render"
    assert r["D"] == pytest.approx(0.0) and r["O"] == pytest.approx(2.0) and r["S0.45"] == pytest.approx(1.1)
    assert r["own_minus_fgw"]["e725"] == pytest.approx(2.0) and r["tasks"]["half"] == pytest.approx(0.0)


def test_a_gain_under_the_second_uploads_lead_counts_for_nothing(n37):
    arm = flat(n37, **{n: 2.5 for n in n37.nt.DEEP})  # +2.5 on the DEEP tasks, where the second upload leads by 3
    scores, tasks = tables(n37, {"fgw_aligned_hi": arm})
    r = n37.price(scores, tasks, candidates=("fgw_aligned_hi",))["first_uploads"]["fgw_aligned_hi"]
    assert r["D"] == pytest.approx(0.0) and r["O"] == pytest.approx(0.0) and not r["passes"]
    assert not r["rule"]["a_S0.45_above_bar"]


def test_one_task_carrying_the_gain_fails_the_leave_one_out_clause(n37):
    arm = flat(n37, e725=4.0)  # S(0.45) = 0.55 x 4/3 = 0.73, all of it from one task
    scores, tasks = tables(n37, {"fgw_trend_render": arm})
    r = n37.price(scores, tasks, candidates=("fgw_trend_render",))["first_uploads"]["fgw_trend_render"]
    assert r["rule"]["a_S0.45_above_bar"] and r["rule"]["b_S0.3_and_S0.6_not_lower"]
    assert not r["rule"]["d_leave_one_task_out_positive"] and not r["passes"]
    assert r["S0.45_leave_one_task_out"]["e725"] == pytest.approx(0.0)


def test_a_loss_on_a_task_where_fgw_carries_the_pair_lowers_the_value(n37):
    arm = flat(n37, e725=4.0, new2=3.0, e75=-1.0)  # a loss on e75, where fgw leads the pair as built
    scores, tasks = tables(n37, {"fgw_trend_render": arm})
    r = n37.price(scores, tasks, candidates=("fgw_trend_render",))["first_uploads"]["fgw_trend_render"]
    assert r["tasks"]["e75"] == pytest.approx(-1.0) and r["O"] == pytest.approx(2.0)
    assert r["rule"]["a_S0.45_above_bar"] and r["rule"]["d_leave_one_task_out_positive"]
    # the pair as built is worst on the OTHER tasks (60, fgw's own), and the arm takes e75 down to 59: clause (c)
    assert r["worst9_abs"] == pytest.approx(59.0) and not r["rule"]["c_worst9_not_lower"] and not r["passes"]


def test_a_lower_worst_of_nine_fails_clause_c(n37):
    second = {n: -2.0 for n in n37.nt.DEEP}  # the second upload trails everywhere: fgw carries the pair on all nine
    arm = flat(n37, new2=5.0, e75=5.0, e725=5.0, half=-1.0)
    scores, tasks = tables(n37, {"fgw_trend_render": arm}, second=second)
    r = n37.price(scores, tasks, candidates=("fgw_trend_render",))["first_uploads"]["fgw_trend_render"]
    assert r["rule"]["a_S0.45_above_bar"] and not r["rule"]["c_worst9_not_lower"] and not r["passes"]


def test_an_arm_missing_on_a_task_is_not_readable_and_never_passes(n37):
    arm = flat(n37, new2=2.0, e75=2.0, e725=2.0)
    arm["x3"] = None  # no scores on one task
    other = flat(n37, new2=2.0, e75=2.0, e725=2.0)
    scores, tasks = tables(n37, {"fgw_early_specialist": arm, "fgw_trend_render": other})
    v = n37.price(scores, tasks, candidates=("fgw_early_specialist", "fgw_trend_render"))
    r = v["first_uploads"]["fgw_early_specialist"]
    assert not r["readable"] and r["missing_on"] == ["x3"] and not r["passes"]
    assert v["requested_first_upload"] == "fgw_trend_render"
    scores[tasks["half"].key]["fgw_trend_render"] = np.r_[np.full(9, 60.0), np.nan]
    v = n37.price(scores, tasks, candidates=("fgw_trend_render",))
    assert not v["first_uploads"]["fgw_trend_render"]["readable"]


def test_with_several_passing_the_largest_s045_is_requested_and_none_if_none(n37):
    a = flat(n37, new2=2.0, e75=2.0, e725=2.0)
    b = flat(n37, new2=3.0, e75=3.0, e725=3.0)
    scores, tasks = tables(n37, {"fgw_trend_render": a, "fgw_aligned_trend": b})
    v = n37.price(scores, tasks, candidates=("fgw_trend_render", "fgw_aligned_trend"))
    assert v["requested_first_upload"] == "fgw_aligned_trend"
    scores, tasks = tables(n37, {"fgw_trend_render": flat(n37)})
    assert n37.price(scores, tasks, candidates=("fgw_trend_render",))["requested_first_upload"] is None


def test_the_reading_is_void_without_the_pair_as_built_and_the_text_names_every_arm(n37, monkeypatch):
    arm = flat(n37, new2=2.0, e75=2.0, e725=2.0)
    scores, tasks = tables(n37, {"fgw_trend_render": arm})
    monkeypatch.setattr(n37.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n37.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n37.reading([])
    assert not res["void"] and res["requested_first_upload"] == "fgw_trend_render"
    txt = n37.text(res)
    assert "fgw_trend_render" in txt and "not readable" in txt and "fgw_early_specialist" in txt
    del scores[tasks["half"].key][n37.SECOND]
    res = n37.reading([])
    assert res["void"] and any("no scores for fgw_aligned_series" in w for w in res["why"])
    assert n37.text(res).startswith("VOID")
    monkeypatch.setattr(n37.nt, "_load", lambda dirs: (scores, ["x: a dirty tree"], []))
    assert n37.reading([])["void"]
