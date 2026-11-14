"""Readings of tries E41 (experiments/harness/n45_embryo_f_readings.py): three existing rules applied to two new arms
and a check of the rescored entry, on synthetic score tables built here. Nothing is read from data/ or from a
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
def n45():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n45_embryo_f_readings", HARNESS / "n45_embryo_f_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n45, gains: dict[str, dict[str, float | None]] | None = None, spread: float = 0.0):
    """Per-seed scores on the nine tasks. Committed arms: the floor at 50; fgw at 60; the second upload as built at
    63 on the DEEP tasks and 58 on the others. Each new arm at fgw + gains[arm][task] (default 0; None = no scores).
    spread adds a per-seed ramp to every arm alike, which leaves paired gains as set."""
    nt = n45.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    gains = gains or {}
    ramp = spread * np.arange(10)
    scores = {}
    for n in names:
        s = {"floor": 50.0 + 0 * ramp, n45.ENTRY: 60.0 + ramp, n45.SECOND: (63.0 if n in nt.DEEP else 58.0) + ramp}
        for arm in n45.ARMS:
            v = gains.get(arm, {}).get(n, 0.0)
            if v is not None:
                s[arm] = s[n45.BASE_OF[arm]] + v
        scores[tasks[n].key] = s
    return scores, tasks


def split(n45, scores: dict) -> tuple[dict, dict]:
    """The table as the reading meets it: the committed runs (no new arm) and the run (new arms and rescored arms)."""
    committed = {k: {a: v for a, v in s.items() if a not in n45.ARMS} for k, s in scores.items()}
    run = {k: {a: v.copy() for a, v in s.items() if a in n45.ARMS or a in n45.RESCORED} for k, s in scores.items()}
    return committed, run


def patch(n45, monkeypatch, committed: dict, run: dict, tasks: dict) -> None:
    monkeypatch.setattr(n45.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n45.nt, "_load", lambda dirs: ((run, [], ["run:1111111"]) if list(dirs) == [RUN]
                                                      else (committed, [], ["committed:0000000"])))


def test_the_constants_are_the_rows(n45):
    assert n45.ARMS == ("fgw_potts", "fgw_partial") and n45.ADOPTION_ARMS == n45.ARMS
    assert n45.BASE_OF == {"fgw_potts": "fgw", "fgw_partial": "fgw"}
    assert n45.RESCORED == ("floor", "fgw")
    assert n45.ENTRY == "fgw" and n45.SECOND == "fgw_aligned_series" and n45.BACKTESTS == ("e75", "e725")


def test_arms_equal_to_the_entry_close_the_line(n45):
    scores, tasks = tables(n45, spread=0.1)
    res = n45.check(scores, tasks)
    assert res["closed"] and res["adopted"] == []
    assert res["requested_first_upload"] is None and res["requested_second_upload"] is None
    assert res["adoption"]["fgw_potts"]["tasks"]["e725"]["mean_gain"] == pytest.approx(0.0)
    assert res["seeds_equal_to_base"]["fgw_partial"]["half"] == 10


def test_a_gain_on_both_backtests_passes_the_adoption_rule_for_either_arm(n45):
    for arm in n45.ARMS:
        scores, tasks = tables(n45, {arm: {"e75": 0.5, "e725": 0.6}})
        res = n45.check(scores, tasks)
        a = res["adoption"][arm]
        assert a["keep"] and a["tasks"]["e75"]["seeds_ahead"] == 10
        assert a["tasks"]["e725"]["mean_gain"] == pytest.approx(0.6)
        assert res["adopted"] == [arm] and not res["closed"]


def test_a_gain_on_one_backtest_and_a_loss_on_the_other_is_not_adopted(n45):
    scores, tasks = tables(n45, {"fgw_partial": {"e75": 2.0, "e725": -0.5}})
    res = n45.check(scores, tasks)
    a = res["adoption"]["fgw_partial"]
    assert not a["keep"] and a["tasks"]["e725"]["vetoes"] and res["adopted"] == []


def test_the_first_upload_rule_reads_the_pair_beside_the_second_upload(n45):
    scores, tasks = tables(n45, {"fgw_potts": {n: 1.0 for n in n45.nt.OTHER}})
    res = n45.check(scores, tasks)
    v = res["first_upload"]["first_uploads"]["fgw_potts"]
    assert v["O"] == pytest.approx(1.0) and v["D"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_first_upload"] == "fgw_potts" and not res["closed"]
    scores, tasks = tables(n45, {"fgw_potts": {n: 1.0 for n in n45.nt.DEEP}})  # hidden behind the second upload
    res = n45.check(scores, tasks)
    assert res["requested_first_upload"] is None


def test_the_second_upload_rule_reads_the_pair_beside_fgw_against_the_second_upload_as_built(n45):
    # the entry trails the second upload by 3 on the DEEP tasks: +1 there is still under it, +4.5 leads it by 1.5
    scores, tasks = tables(n45, {"fgw_partial": {n: 1.0 for n in n45.nt.DEEP}})
    res = n45.check(scores, tasks)
    assert res["requested_second_upload"] is None
    scores, tasks = tables(n45, {"fgw_partial": {n: 4.5 for n in n45.nt.DEEP}})
    res = n45.check(scores, tasks)
    v = res["second_upload"]["second_uploads"]["fgw_partial"]
    assert v["D"] == pytest.approx(1.5) and v["O"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_second_upload"] == "fgw_partial" and not res["closed"]


def test_an_arm_missing_on_a_task_is_not_readable_and_cannot_pass(n45):
    scores, tasks = tables(n45, {"fgw_partial": {"e75": None, "e725": 5.0}})
    res = n45.check(scores, tasks)
    assert not res["readable"]["fgw_partial"] and "fgw_partial" not in res["own_minus_base"]
    assert res["adoption"]["fgw_partial"] is None
    assert not res["first_upload"]["first_uploads"]["fgw_partial"]["readable"]
    assert res["readable"]["fgw_potts"] and res["closed"]


def test_the_count_of_seeds_equal_to_the_entry_is_reported(n45):
    scores, tasks = tables(n45)
    key = tasks["e725"].key
    arm = scores[key]["fgw_potts"].copy()
    arm[[1, 4, 7]] += 0.8
    arm[2] += 0.004  # a difference inside the tolerance between runs still counts as the entry
    scores[key]["fgw_potts"] = arm
    res = n45.check(scores, tasks)
    assert res["seeds_equal_to_base"]["fgw_potts"]["e725"] == 7


def test_the_rescored_entry_must_match_the_record(n45):
    scores, tasks = tables(n45, spread=0.1)
    committed, run = split(n45, scores)
    worst, why = n45.rescored(committed, run, tasks)
    assert why == [] and set(worst) == set(n45.RESCORED) and max(worst.values()) == 0.0
    run[tasks["e725"].key]["fgw"][0] += 0.02
    worst, why = n45.rescored(committed, run, tasks)
    assert worst["fgw"] == pytest.approx(0.02)
    assert why == ["e725: fgw differs from its committed scores by 2.00e-02"]
    del run[tasks["half"].key]["fgw"]
    _, why = n45.rescored(committed, run, tasks)
    assert "half: the run did not rescore fgw" in why


def test_the_reading_uses_the_committed_entry_and_the_runs_new_arms(n45, monkeypatch):
    scores, tasks = tables(n45, {"fgw_potts": {"e75": 0.5, "e725": 0.6}})
    committed, run = split(n45, scores)
    run[tasks["new1"].key]["fgw"] += 0.005  # the rescored entry inside the tolerance: the record's is used
    patch(n45, monkeypatch, committed, run, tasks)
    res = n45.reading([RUN])
    assert not res["void"] and res["runs"] == ["committed:0000000", "run:1111111"]
    assert res["rescored_max_abs_difference"]["fgw"] == pytest.approx(0.005)
    assert res["own_minus_base"]["fgw_potts"]["new1"] == pytest.approx(0.0)
    assert res["adopted"] == ["fgw_potts"]
    txt = n45.text(res)
    assert "the adoption rule is passed by fgw_potts" in txt and "Readings of tries E41a and E41b" in txt
    assert "rescored arms against their committed scores" in txt
    assert "seeds (of 10) in which fgw_partial scores as fgw:" in txt


def test_the_reading_is_void_on_a_rescored_mismatch_or_a_missing_task(n45, monkeypatch):
    scores, tasks = tables(n45)
    committed, run = split(n45, scores)
    run[tasks["late"].key]["fgw"][5] -= 0.3
    patch(n45, monkeypatch, committed, run, tasks)
    res = n45.reading([RUN])
    assert res["void"] and any("late: fgw differs from its committed scores" in w for w in res["why"])
    assert n45.text(res).startswith("VOID")
    committed, run = split(n45, scores)
    del run[tasks["x1"].key]
    patch(n45, monkeypatch, committed, run, tasks)
    res = n45.reading([RUN])
    assert res["void"] and any("x1: no task file in the run" in w for w in res["why"])
    committed, run = split(n45, scores)
    patch(n45, monkeypatch, committed, run, tasks)
    res = n45.reading([RUN])
    assert not res["void"] and res["closed"] and "closed: no rule is passed" in n45.text(res)
