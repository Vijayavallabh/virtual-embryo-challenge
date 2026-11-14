"""Readings of tries E40 (experiments/harness/n44_embryo_e_readings.py): three existing rules applied to four new
arms and a check of the rescored arms, on synthetic score tables built here. Nothing is read from data/ or from a
committed run, and nothing is scored."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"
RUN = Path("RUN")


@pytest.fixture(scope="module")
def n44():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n44_embryo_e_readings", HARNESS / "n44_embryo_e_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n44, gains: dict[str, dict[str, float | None]] | None = None, spread: float = 0.0):
    """Per-seed scores on the nine tasks. Committed arms: the floor at 50; fgw and fgw_trend_render at 60; the second
    upload as built at 63 on the DEEP tasks and 58 on the others; fgw_aligned_trend at 61 and 57. Each new arm at its
    base + gains[arm][task] (default 0; None = no scores). spread adds a per-seed ramp to every arm alike, which
    leaves paired gains as set."""
    nt = n44.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    gains = gains or {}
    ramp = spread * np.arange(10)
    scores = {}
    for n in names:
        deep = n in nt.DEEP
        s = {"floor": 50.0 + 0 * ramp, n44.ENTRY: 60.0 + ramp, n44.LINE_FIRST: 60.0 + ramp,
             n44.SECOND: (63.0 if deep else 58.0) + ramp, n44.LINE_SECOND: (61.0 if deep else 57.0) + ramp}
        for arm in n44.ARMS:
            v = gains.get(arm, {}).get(n, 0.0)
            if v is not None:
                s[arm] = s[n44.BASE_OF[arm]] + v
        scores[tasks[n].key] = s
    return scores, tasks


def split(n44, scores: dict) -> tuple[dict, dict]:
    """The table as the reading meets it: the committed runs (no new arm) and the run (new arms and rescored arms)."""
    committed = {k: {a: v for a, v in s.items() if a not in n44.ARMS} for k, s in scores.items()}
    run = {k: {a: v.copy() for a, v in s.items() if a in n44.ARMS or a in n44.RESCORED} for k, s in scores.items()}
    return committed, run


def patch(n44, monkeypatch, committed: dict, run: dict, tasks: dict) -> None:
    monkeypatch.setattr(n44.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n44.nt, "_load", lambda dirs: ((run, [], ["run:1111111"]) if list(dirs) == [RUN]
                                                      else (committed, [], ["committed:0000000"])))


def test_the_constants_are_the_rows(n44):
    assert n44.ARMS == ("fgw_aligned_mix", "fgw_aligned_series_mn", "fgw_quad_render", "fgw_aligned_quad")
    assert n44.ADOPTION_ARMS == ("fgw_quad_render",)
    assert n44.BASE_OF == {"fgw_aligned_mix": "fgw_aligned_series", "fgw_aligned_series_mn": "fgw_aligned_series",
                           "fgw_quad_render": "fgw_trend_render", "fgw_aligned_quad": "fgw_aligned_trend"}
    assert n44.ALSO_AGAINST == {"fgw_aligned_mix": "fgw_aligned_trend", "fgw_quad_render": "fgw",
                                "fgw_aligned_quad": "fgw_aligned_series"}
    assert n44.RESCORED == ("floor", "fgw", "fgw_aligned_series", "fgw_trend_render", "fgw_aligned_trend")
    assert n44.ENTRY == "fgw" and n44.SECOND == "fgw_aligned_series" and n44.BACKTESTS == ("e75", "e725")
    assert n44.N36_RUNS.name == "tasks" and n44.N36_RUNS.parent.name == "n36_trend_second_upload"


def test_arms_equal_to_their_bases_close_the_line(n44):
    scores, tasks = tables(n44, spread=0.1)
    res = n44.check(scores, tasks)
    assert res["closed"] and res["adopted"] == []
    assert res["requested_first_upload"] is None and res["requested_second_upload"] is None
    assert res["adoption"]["fgw_quad_render"]["tasks"]["e725"]["mean_gain"] == pytest.approx(0.0)
    assert res["own_minus_base"]["fgw_aligned_series_mn"]["x1"] == pytest.approx(0.0)
    # the mixture equals the second upload here, an end of the interval between its halves
    assert res["mix_tasks_between_its_halves"] == 9


def test_a_gain_on_both_backtests_passes_the_adoption_rule(n44):
    scores, tasks = tables(n44, {"fgw_quad_render": {"e75": 0.5, "e725": 0.6}})
    res = n44.check(scores, tasks)
    a = res["adoption"]["fgw_quad_render"]
    assert a["keep"] and a["tasks"]["e75"]["seeds_ahead"] == 10
    assert a["tasks"]["e725"]["mean_gain"] == pytest.approx(0.6)
    assert res["adopted"] == ["fgw_quad_render"] and not res["closed"]


def test_the_adoption_rule_reads_the_quadratic_against_the_entry_not_against_the_line(n44):
    # the line's arm trails the entry by 1 on both backtests; the quadratic adds 0.5 to the line: -0.5 on the entry
    scores, tasks = tables(n44, {"fgw_quad_render": {"e75": 0.5, "e725": 0.5}})
    for n in n44.BACKTESTS:
        k = tasks[n].key
        scores[k][n44.LINE_FIRST] = scores[k][n44.LINE_FIRST] - 1.0
        scores[k]["fgw_quad_render"] = scores[k][n44.LINE_FIRST] + 0.5
    res = n44.check(scores, tasks)
    a = res["adoption"]["fgw_quad_render"]
    assert not a["keep"] and a["tasks"]["e75"]["mean_gain"] == pytest.approx(-0.5)
    assert res["own_minus_base"]["fgw_quad_render"]["e75"] == pytest.approx(0.5)
    assert res["own_minus_other"]["fgw_quad_render"]["e75"] == pytest.approx(-0.5)


def test_a_gain_on_one_backtest_and_a_loss_on_the_other_is_not_adopted(n44):
    scores, tasks = tables(n44, {"fgw_quad_render": {"e75": 2.0, "e725": -0.5}})
    res = n44.check(scores, tasks)
    a = res["adoption"]["fgw_quad_render"]
    assert not a["keep"] and a["tasks"]["e725"]["vetoes"] and res["adopted"] == []


def test_the_first_upload_rule_reads_the_pair_beside_the_second_upload(n44):
    scores, tasks = tables(n44, {"fgw_quad_render": {n: 1.0 for n in n44.nt.OTHER}})
    res = n44.check(scores, tasks)
    v = res["first_upload"]["first_uploads"]["fgw_quad_render"]
    assert v["O"] == pytest.approx(1.0) and v["D"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_first_upload"] == "fgw_quad_render" and not res["closed"]
    scores, tasks = tables(n44, {"fgw_quad_render": {n: 1.0 for n in n44.nt.DEEP}})  # hidden behind the second upload
    res = n44.check(scores, tasks)
    assert res["requested_first_upload"] is None


def test_the_second_upload_rule_reads_the_pair_beside_fgw_against_the_second_upload_as_built(n44):
    for arm in ("fgw_aligned_mix", "fgw_aligned_series_mn"):
        scores, tasks = tables(n44, {arm: {n: 1.0 for n in n44.nt.DEEP}})
        res = n44.check(scores, tasks)
        v = res["second_upload"]["second_uploads"][arm]
        assert v["D"] == pytest.approx(1.0) and v["O"] == pytest.approx(0.0) and v["passes"]
        assert res["requested_second_upload"] == arm and not res["closed"]
    # the aligned trend trails the second upload by 2 on the DEEP tasks: +1 there is still under it
    scores, tasks = tables(n44, {"fgw_aligned_quad": {n: 1.0 for n in n44.nt.DEEP}})
    res = n44.check(scores, tasks)
    assert res["requested_second_upload"] is None and res["closed"]
    assert res["own_minus_other"]["fgw_aligned_quad"]["half"] == pytest.approx(-1.0)


def test_an_arm_missing_on_a_task_is_not_readable_and_cannot_pass(n44):
    gains = {"fgw_aligned_mix": {**{n: 5.0 for n in n44.nt.DEEP}, "x3": None},
             "fgw_quad_render": {"e75": None, "e725": 5.0}}
    scores, tasks = tables(n44, gains)
    res = n44.check(scores, tasks)
    assert not res["readable"]["fgw_aligned_mix"] and "fgw_aligned_mix" not in res["own_minus_base"]
    assert res["mix_tasks_between_its_halves"] is None
    assert not res["second_upload"]["second_uploads"]["fgw_aligned_mix"]["readable"]
    assert res["adoption"]["fgw_quad_render"] is None and res["adopted"] == [] and res["closed"]


def test_the_mixture_between_its_halves_is_counted(n44):
    # half-way between the second upload and the aligned trend on the DEEP tasks (-1), above both on the others (+2)
    gains = {"fgw_aligned_mix": {**{n: -1.0 for n in n44.nt.DEEP}, **{n: 2.0 for n in n44.nt.OTHER}}}
    scores, tasks = tables(n44, gains)
    res = n44.check(scores, tasks)
    assert res["mix_tasks_between_its_halves"] == len(n44.nt.DEEP)
    assert res["own_minus_other"]["fgw_aligned_mix"]["half"] == pytest.approx(1.0)


def test_the_rescored_arms_must_match_the_record(n44):
    scores, tasks = tables(n44, spread=0.1)
    committed, run = split(n44, scores)
    worst, why = n44.rescored(committed, run, tasks)
    assert why == [] and set(worst) == set(n44.RESCORED) and max(worst.values()) == 0.0
    run[tasks["x2"].key]["fgw_aligned_series"][3] += 0.004  # inside the tolerance between runs
    worst, why = n44.rescored(committed, run, tasks)
    assert why == [] and worst["fgw_aligned_series"] == pytest.approx(0.004)
    run[tasks["e725"].key]["fgw"][0] += 0.02  # outside it
    worst, why = n44.rescored(committed, run, tasks)
    assert worst["fgw"] == pytest.approx(0.02)
    assert why == ["e725: fgw differs from its committed scores by 2.00e-02"]
    del run[tasks["half"].key]["fgw_trend_render"]
    _, why = n44.rescored(committed, run, tasks)
    assert "half: the run did not rescore fgw_trend_render" in why


def test_the_reading_uses_the_committed_bases_and_the_runs_new_arms(n44, monkeypatch):
    scores, tasks = tables(n44, {"fgw_aligned_series_mn": {n: 1.0 for n in n44.nt.DEEP}})
    committed, run = split(n44, scores)
    run[tasks["new1"].key]["fgw_aligned_series"] += 0.005  # a rescored arm inside the tolerance: the record's is used
    patch(n44, monkeypatch, committed, run, tasks)
    res = n44.reading([RUN])
    assert not res["void"] and res["runs"] == ["committed:0000000", "run:1111111"]
    assert res["rescored_max_abs_difference"]["fgw_aligned_series"] == pytest.approx(0.005)
    assert res["own_minus_base"]["fgw_aligned_series_mn"]["new1"] == pytest.approx(1.0)
    assert res["requested_second_upload"] == "fgw_aligned_series_mn"
    txt = n44.text(res)
    assert "second upload put forward: fgw_aligned_series_mn" in txt
    assert "rescored arms against their committed scores" in txt and "Readings of tries E40a to E40d" in txt
    assert "fgw_aligned_mix: its own mean lies between its two halves' means" in txt
    assert "fgw_quad_render minus fgw_trend_render, own means:" in txt and "    minus fgw: " in txt


def test_the_reading_is_void_on_a_rescored_mismatch_a_missing_base_or_a_missing_task(n44, monkeypatch):
    scores, tasks = tables(n44)
    committed, run = split(n44, scores)
    run[tasks["late"].key]["fgw_aligned_trend"][5] -= 0.3
    patch(n44, monkeypatch, committed, run, tasks)
    res = n44.reading([RUN])
    assert res["void"] and any("late: fgw_aligned_trend differs from its committed scores" in w for w in res["why"])
    assert n44.text(res).startswith("VOID")
    committed, run = split(n44, scores)
    del committed[tasks["half"].key]["fgw_trend_render"]
    patch(n44, monkeypatch, committed, run, tasks)
    res = n44.reading([RUN])
    assert res["void"] and any("half: no committed scores for fgw_trend_render" in w for w in res["why"])
    committed, run = split(n44, scores)
    del run[tasks["x1"].key]
    patch(n44, monkeypatch, committed, run, tasks)
    res = n44.reading([RUN])
    assert res["void"] and any("x1: no task file in the run" in w for w in res["why"])


def test_arms_equal_to_their_bases_read_as_closed_in_the_text(n44, monkeypatch):
    scores, tasks = tables(n44)
    committed, run = split(n44, scores)
    patch(n44, monkeypatch, committed, run, tasks)
    res = n44.reading([RUN])
    assert not res["void"] and res["closed"]
    txt = n44.text(res)
    assert "closed: no rule is passed" in txt and "adoption rule, fgw_quad_render against fgw" in txt
    scores2, _ = tables(n44, {"fgw_quad_render": {"e75": 0.5, "e725": 0.6}})
    committed, run = split(n44, scores2)
    patch(n44, monkeypatch, committed, run, tasks)
    assert "the adoption rule is passed by fgw_quad_render" in n44.text(n44.reading([RUN]))
