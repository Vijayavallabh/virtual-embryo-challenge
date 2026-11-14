"""Readings of tries E42 (experiments/harness/n46_embryo_g_readings.py): the controls for fgw_partial, on synthetic
score tables built here. Nothing is read from data/ or from a committed run, and nothing is scored."""

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
def n46():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n46_embryo_g_readings", HARNESS / "n46_embryo_g_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n46, gains: dict[str, float | dict[str, float | None]] | None = None):
    """Per-seed scores on the nine tasks. The floor at 50; fgw at 60 plus a per-seed ramp; the second upload as built
    at 63 on the DEEP tasks and 58 on the others. Each of the four arms at fgw + its gain: one number for every task,
    or a dict per task (default 0; None = no scores)."""
    nt = n46.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    gains = gains or {}
    ramp = 0.1 * np.arange(10)
    scores = {}
    for n in names:
        s = {"floor": 50.0 + 0 * ramp, n46.ENTRY: 60.0 + ramp, n46.SECOND: (63.0 if n in nt.DEEP else 58.0) + ramp}
        for arm in n46.ORDER:
            g = gains.get(arm, 0.0)
            v = g.get(n, 0.0) if isinstance(g, dict) else g
            if v is not None:
                s[arm] = s[n46.ENTRY] + v
        scores[tasks[n].key] = s
    return scores, tasks


def split(n46, scores: dict) -> tuple[dict, dict]:
    """The table as the reading meets it: the committed runs (fgw_partial among them, no control) and the run (the
    controls and the rescored arms)."""
    committed = {k: {a: v for a, v in s.items() if a not in n46.CONTROLS} for k, s in scores.items()}
    run = {k: {a: v.copy() for a, v in s.items() if a in n46.CONTROLS or a in n46.RESCORED} for k, s in scores.items()}
    return committed, run


def patch(n46, monkeypatch, committed: dict, run: dict, tasks: dict) -> None:
    monkeypatch.setattr(n46.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n46.nt, "_load", lambda dirs: ((run, [], ["run:1111111"]) if list(dirs) == [RUN]
                                                      else (committed, [], ["committed:0000000"])))


def test_the_constants_are_the_rows(n46):
    assert n46.CONTROLS == ("fgw_still", "fgw_damp", "fgw_part_rand") and n46.PARTIAL == "fgw_partial"
    assert n46.LEAD_OVER == ("fgw_damp", "fgw_part_rand")
    assert n46.ORDER == ("fgw_still", "fgw_damp", "fgw_part_rand", "fgw_partial")
    assert n46.RESCORED == ("floor", "fgw", "fgw_partial")
    assert n46.ENTRY == "fgw" and n46.SECOND == "fgw_aligned_series" and n46.BACKTESTS == ("e75", "e725")
    assert n46.N45_RUNS.parent.name == "n45_embryo_f_arms"


def test_partial_ahead_of_both_moving_controls_is_matching_and_is_put_forward(n46):
    # the arm gains 2 everywhere; the controls that move as far gain 1: it leads them by 1 on both backtests
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": 1.0, "fgw_part_rand": 1.0, "fgw_still": 2.5})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "matching" and res["put_forward"] == "fgw_partial" and not res["closed"]
    assert res["lead"]["fgw_damp"]["keep"] and res["lead"]["fgw_part_rand"]["keep"]
    assert res["passes"]["fgw_still"], "a passing control does not displace the arm when the verdict is matching"


def test_partial_level_with_a_moving_control_is_moving_less_and_the_simplest_passing_arm_is_put_forward(n46):
    # a gain g on every task gives the pair rule S(0.45) = 0.55 g (the arm leads the pair on the three OTHER tasks
    # only), so 0.6 passes it and 0.5 would not
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": 1.9, "fgw_part_rand": 0.6, "fgw_still": 0.0})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less"
    assert not res["lead"]["fgw_damp"]["keep"] and res["lead"]["fgw_part_rand"]["keep"]
    assert res["passes"] == {"fgw_still": False, "fgw_damp": True, "fgw_part_rand": True, "fgw_partial": True}
    assert res["put_forward"] == "fgw_damp"
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": 1.9, "fgw_part_rand": 1.9, "fgw_still": 1.8})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less" and res["put_forward"] == "fgw_still"


def test_moving_less_with_no_passing_control_leaves_partial_as_the_arm(n46):
    # the controls gain 0.52 on every task: ahead of fgw on both backtests, but S(0.45) = 0.286 is under the pair
    # rule's bar; the arm gains 0.8 and passes both rules, and its lead of 0.28 over each control is under the bar
    scores, tasks = tables(n46, {"fgw_partial": 0.8, "fgw_damp": 0.52, "fgw_part_rand": 0.52, "fgw_still": 0.52})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less"
    assert res["passes"] == {"fgw_still": False, "fgw_damp": False, "fgw_part_rand": False, "fgw_partial": True}
    assert res["put_forward"] == "fgw_partial"


def test_nothing_passes_closes_the_line(n46):
    scores, tasks = tables(n46)
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less" and res["put_forward"] is None and res["closed"]
    assert res["own_minus_fgw"]["fgw_still"]["x2"] == pytest.approx(0.0)


def test_an_arm_needs_both_rules_to_pass(n46):
    # ahead of fgw on both backtests but hidden behind the second upload on the other tasks is still a pass of the
    # pair rule here (e75 and e725 are OTHER tasks); a gain on the DEEP tasks alone passes neither rule
    scores, tasks = tables(n46, {"fgw_still": {n: 1.0 for n in n46.nt.DEEP}})
    res = n46.check(scores, tasks)
    assert not res["adoption"]["fgw_still"]["keep"] and not res["passes"]["fgw_still"]
    scores, tasks = tables(n46, {"fgw_still": {"e75": 1.0, "e725": 1.0, "new2": -3.0}})
    res = n46.check(scores, tasks)
    assert res["adoption"]["fgw_still"]["keep"]
    assert not res["first_upload"]["first_uploads"]["fgw_still"]["passes"] and not res["passes"]["fgw_still"]


def test_a_control_missing_on_a_task_is_not_readable_and_cannot_lead_or_pass(n46):
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": {"e75": None}, "fgw_part_rand": 1.0})
    res = n46.check(scores, tasks)
    assert not res["readable"]["fgw_damp"] and res["lead"]["fgw_damp"] is None
    assert res["verdict"] == "moving less", "an unreadable control cannot make the verdict matching"
    assert res["adoption"]["fgw_damp"] is None and not res["passes"]["fgw_damp"]
    assert "fgw_damp" not in res["partial_minus_control"]


def test_the_rescored_arms_must_match_the_record(n46):
    scores, tasks = tables(n46, {"fgw_partial": 1.0})
    committed, run = split(n46, scores)
    worst, why = n46.rescored(committed, run, tasks)
    assert why == [] and set(worst) == set(n46.RESCORED) and max(worst.values()) == 0.0
    run[tasks["e725"].key]["fgw_partial"][0] += 0.02
    worst, why = n46.rescored(committed, run, tasks)
    assert worst["fgw_partial"] == pytest.approx(0.02)
    assert why == ["e725: fgw_partial differs from its committed scores by 2.00e-02"]
    del run[tasks["half"].key]["fgw"]
    _, why = n46.rescored(committed, run, tasks)
    assert "half: the run did not rescore fgw" in why


def test_the_reading_uses_the_committed_arms_and_the_runs_controls(n46, monkeypatch):
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": 1.9, "fgw_part_rand": 1.9, "fgw_still": 0.1})
    committed, run = split(n46, scores)
    run[tasks["new1"].key]["fgw_partial"] += 0.005  # a rescored arm inside the tolerance: the record's is used
    patch(n46, monkeypatch, committed, run, tasks)
    res = n46.reading([RUN])
    assert not res["void"] and res["runs"] == ["committed:0000000", "run:1111111"]
    assert res["rescored_max_abs_difference"]["fgw_partial"] == pytest.approx(0.005)
    assert res["own_minus_fgw"]["fgw_partial"]["new1"] == pytest.approx(2.0)
    assert res["verdict"] == "moving less" and res["put_forward"] == "fgw_damp"
    txt = n46.text(res)
    assert "verdict: moving less; arm put forward: fgw_damp (to a cold check before any request)" in txt
    assert "fgw_partial against fgw_damp:" in txt and "does not lead" in txt
    assert "Readings of tries E42a to E42c" in txt and "rescored arms against their committed scores" in txt


def test_the_reading_is_void_on_a_rescored_mismatch_a_missing_arm_or_a_missing_task(n46, monkeypatch):
    scores, tasks = tables(n46, {"fgw_partial": 1.0})
    committed, run = split(n46, scores)
    run[tasks["late"].key]["fgw"][5] -= 0.3
    patch(n46, monkeypatch, committed, run, tasks)
    res = n46.reading([RUN])
    assert res["void"] and any("late: fgw differs from its committed scores" in w for w in res["why"])
    assert n46.text(res).startswith("VOID")
    committed, run = split(n46, scores)
    del committed[tasks["x3"].key]["fgw_partial"]
    patch(n46, monkeypatch, committed, run, tasks)
    res = n46.reading([RUN])
    assert res["void"] and any("x3: no committed scores for fgw_partial" in w for w in res["why"])
    committed, run = split(n46, scores)
    del run[tasks["x1"].key]
    patch(n46, monkeypatch, committed, run, tasks)
    res = n46.reading([RUN])
    assert res["void"] and any("x1: no task file in the run" in w for w in res["why"])
    committed, run = split(n46, scores)
    patch(n46, monkeypatch, committed, run, tasks)
    res = n46.reading([RUN])
    assert not res["void"] and "arm put forward: fgw_partial" in n46.text(res)
