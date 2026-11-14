"""Reading of try E43 (experiments/harness/n47_embryo_h_readings.py): is the partial coupling's pass robust to the
label level, on synthetic score tables built here. Nothing is read from data/ or from a committed run, and nothing
is scored."""

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
def n47():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n47_embryo_h_readings", HARNESS / "n47_embryo_h_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n47, arm: float | dict | None = 0.0, partial: float = 1.0):
    """Per-seed scores on the nine tasks. The floor at 50; fgw at 60 plus a per-seed ramp; the second upload as built
    at 63 on the DEEP tasks and 58 on the others; fgw_partial at fgw + partial; the arm at fgw + its gain (one number
    for every task, a dict per task, or None for no scores at all)."""
    nt = n47.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    ramp = 0.1 * np.arange(10)
    scores = {}
    for n in names:
        s = {"floor": 50.0 + 0 * ramp, n47.ENTRY: 60.0 + ramp, n47.SECOND: (63.0 if n in nt.DEEP else 58.0) + ramp}
        s[n47.PARTIAL] = s[n47.ENTRY] + partial
        v = arm.get(n, 0.0) if isinstance(arm, dict) else arm
        if v is not None:
            s[n47.ARM] = s[n47.ENTRY] + v
        scores[tasks[n].key] = s
    return scores, tasks


def split(n47, scores: dict) -> tuple[dict, dict]:
    committed = {k: {a: v for a, v in s.items() if a != n47.ARM} for k, s in scores.items()}
    run = {k: {a: v.copy() for a, v in s.items() if a == n47.ARM or a in n47.RESCORED} for k, s in scores.items()}
    return committed, run


def patch(n47, monkeypatch, committed: dict, run: dict, tasks: dict) -> None:
    monkeypatch.setattr(n47.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n47.nt, "_load", lambda dirs: ((run, [], ["run:1111111"]) if list(dirs) == [RUN]
                                                      else (committed, [], ["committed:0000000"])))


def test_the_constants_are_the_rows(n47):
    assert n47.ARM == "fgw_partial_grp" and n47.PARTIAL == "fgw_partial" and n47.ENTRY == "fgw"
    assert n47.RESCORED == ("floor", "fgw", "fgw_partial") and n47.BACKTESTS == ("e75", "e725")
    assert n47.N45_RUNS.parent.name == "n45_embryo_f_arms"


def test_an_arm_that_passes_both_rules_makes_the_form_label_robust(n47):
    # a gain g on every task gives the pair rule S(0.45) = 0.55 g, so 1.0 passes both rules
    scores, tasks = tables(n47, arm=1.0)
    res = n47.check(scores, tasks)
    assert res["label_robust"] and res["adoption"]["keep"]
    assert res["first_upload"]["first_uploads"]["fgw_partial_grp"]["passes"]
    assert res["own_minus_partial"]["e75"] == pytest.approx(0.0)


def test_failing_either_rule_makes_it_not_label_robust(n47):
    scores, tasks = tables(n47, arm=0.0)  # the arm is fgw: neither rule
    assert not n47.check(scores, tasks)["label_robust"]
    scores, tasks = tables(n47, arm=0.5)  # ahead on both backtests, but S(0.45) = 0.275 is under the pair rule's bar
    res = n47.check(scores, tasks)
    assert res["adoption"]["keep"] and not res["first_upload"]["first_uploads"]["fgw_partial_grp"]["passes"]
    assert not res["label_robust"]
    scores, tasks = tables(n47, arm={"e75": 2.0, "e725": -0.5, "new2": 2.0})  # a veto on one backtest
    res = n47.check(scores, tasks)
    assert not res["adoption"]["keep"] and not res["label_robust"]


def test_an_arm_missing_on_a_task_is_not_readable_and_not_label_robust(n47):
    scores, tasks = tables(n47, arm={"x1": None, "e75": 1.0, "e725": 1.0})
    res = n47.check(scores, tasks)
    assert not res["readable"] and not res["label_robust"] and res["adoption"] is None
    assert "not readable" in n47.text({"void": False, "runs": ["x:0"], "rescored_max_abs_difference":
                                       {a: 0.0 for a in n47.RESCORED}, **res})


def test_the_reading_checks_the_rescored_arms_and_uses_the_committed_ones(n47, monkeypatch):
    scores, tasks = tables(n47, arm=1.2, partial=1.0)
    committed, run = split(n47, scores)
    run[tasks["new1"].key]["fgw_partial"] += 0.005
    patch(n47, monkeypatch, committed, run, tasks)
    res = n47.reading([RUN])
    assert not res["void"] and res["label_robust"]
    assert res["rescored_max_abs_difference"]["fgw_partial"] == pytest.approx(0.005)
    assert res["own_minus_partial"]["new1"] == pytest.approx(0.2)
    txt = n47.text(res)
    assert "label-robust: the arm passes both rules" in txt and "Reading of try E43" in txt
    run[tasks["late"].key]["fgw"][5] -= 0.3
    res = n47.reading([RUN])
    assert res["void"] and any("late: fgw differs from its committed scores" in w for w in res["why"])
    committed, run = split(n47, scores)
    del run[tasks["x1"].key]
    patch(n47, monkeypatch, committed, run, tasks)
    assert n47.reading([RUN])["void"]
    scores, tasks = tables(n47, arm=0.0)
    committed, run = split(n47, scores)
    patch(n47, monkeypatch, committed, run, tasks)
    assert "NOT label-robust" in n47.text(n47.reading([RUN]))
