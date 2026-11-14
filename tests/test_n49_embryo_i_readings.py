"""Reading of try E44 (experiments/harness/n49_embryo_i_readings.py): the reference-residual render on the embryo, on
synthetic score tables built here. Nothing is read from data/ or from a committed run, and nothing is scored."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"
RUN = Path("RUN")


@pytest.fixture(scope="module")
def n49():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n49_embryo_i_readings", HARNESS / "n49_embryo_i_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n49, arm: float | dict | None = 0.0, trend: float = -1.0):
    """Per-seed scores on the nine tasks. The floor at 50; fgw at 60 plus a per-seed ramp; the second upload as built
    at 63 on the DEEP tasks and 58 on the others; fgw_trend_render at fgw + trend; the arm at fgw + its gain (one
    number for every task, a dict per task, or None for no scores at all)."""
    nt = n49.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    ramp = 0.1 * np.arange(10)
    scores = {}
    for n in names:
        s = {"floor": 50.0 + 0 * ramp, n49.ENTRY: 60.0 + ramp, n49.SECOND: (63.0 if n in nt.DEEP else 58.0) + ramp}
        s[n49.TREND] = s[n49.ENTRY] + trend
        v = arm.get(n, 0.0) if isinstance(arm, dict) else arm
        if v is not None:
            s[n49.ARM] = s[n49.ENTRY] + v
        scores[tasks[n].key] = s
    return scores, tasks


def split(n49, scores: dict) -> tuple[dict, dict]:
    committed = {k: {a: v for a, v in s.items() if a != n49.ARM} for k, s in scores.items()}
    run = {k: {a: v.copy() for a, v in s.items() if a == n49.ARM or a in n49.RESCORED} for k, s in scores.items()}
    return committed, run


def patch(n49, monkeypatch, committed: dict, run: dict, tasks: dict, rec: dict | None = None) -> None:
    monkeypatch.setattr(n49.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n49.nt, "_load", lambda dirs: ((run, [], ["run:1111111"]) if list(dirs) == [RUN]
                                                      else (committed, [], ["committed:0000000"])))
    monkeypatch.setattr(n49, "corrections", lambda dirs: rec or {})


def test_the_constants_are_the_rows(n49):
    assert n49.ARM == "fgw_resid" and n49.TREND == "fgw_trend_render" and n49.ENTRY == "fgw"
    assert n49.SECOND == "fgw_aligned_series" and n49.RESCORED == ("floor", "fgw", "fgw_trend_render")
    assert n49.BACKTESTS == ("e75", "e725") and n49.N36_RUNS.parent.name == "n36_trend_second_upload"


def test_a_gain_everywhere_passes_as_a_first_upload(n49):
    # a gain g on every task gives the first-upload rule S(0.45) = 0.55 g (the arm leads the pair on the three OTHER
    # tasks only), so 1.0 passes it and both backtests; as a second upload beside fgw the pair (fgw, arm) is the arm
    # everywhere, which is below the pair as built on the six DEEP tasks (61 against 63)
    res = n49.check(*tables(n49, 1.0))
    assert res["adoption"]["keep"] and res["first_upload"]["first_uploads"]["fgw_resid"]["passes"]
    assert res["first_passes"] and not res["second_passes"] and not res["closed"]
    assert res["own_minus_fgw"]["x2"] == pytest.approx(1.0) and res["own_minus_trend"]["x2"] == pytest.approx(2.0)


def test_a_one_sided_arm_can_pass_as_a_second_upload_only(n49):
    # far ahead on the DEEP tasks and at fgw on the others: no gain on either counted backtest, so not a first upload;
    # as a second upload the pair (fgw, arm) beats the pair as built on every DEEP task
    scores, tasks = tables(n49, {n: 6.0 for n in n49.nt.DEEP})
    res = n49.check(scores, tasks)
    assert not res["adoption"]["keep"] and not res["first_passes"]
    assert res["second_upload"]["second_uploads"]["fgw_resid"]["passes"] and res["second_passes"] and not res["closed"]


def test_nothing_gained_closes_the_line(n49):
    res = n49.check(*tables(n49, 0.0))
    assert not res["first_passes"] and not res["second_passes"] and res["closed"]
    res = n49.check(*tables(n49, 0.3))
    assert not res["adoption"]["keep"] and res["closed"], "a gain of exactly 0.3 does not pass"


def test_the_first_upload_needs_both_rules(n49):
    # ahead on both backtests, 3 behind on the third OTHER task: the adoption rule keeps it and the pair rule does not
    scores, tasks = tables(n49, {"e75": 1.0, "e725": 1.0, "new2": -3.0})
    res = n49.check(scores, tasks)
    assert res["adoption"]["keep"] and not res["first_upload"]["first_uploads"]["fgw_resid"]["passes"]
    assert not res["first_passes"]


def test_an_arm_missing_on_a_task_is_not_readable_and_closes(n49):
    scores, tasks = tables(n49, {"x1": None})
    res = n49.check(scores, tasks)
    assert not res["readable"] and res["adoption"] is None and res["closed"]


def test_the_rescored_arms_must_match_the_record(n49):
    scores, tasks = tables(n49, 1.0)
    committed, run = split(n49, scores)
    worst, why = n49.rescored(committed, run, tasks)
    assert why == [] and set(worst) == set(n49.RESCORED) and max(worst.values()) == 0.0
    run[tasks["e725"].key]["fgw_trend_render"][0] += 0.02
    worst, why = n49.rescored(committed, run, tasks)
    assert worst["fgw_trend_render"] == pytest.approx(0.02)
    assert why == ["e725: fgw_trend_render differs from its committed scores by 2.00e-02"]
    del run[tasks["half"].key]["fgw"]
    _, why = n49.rescored(committed, run, tasks)
    assert "half: the run did not rescore fgw" in why


def test_the_corrections_are_read_from_the_task_files(n49, tmp_path):
    def rec(clip, size):
        return {"resid": {"rule": "r", "brackets": {
            "before": {"corrected": True, "neighbours": ["E7.25", "E8.0"], "w1": 2 / 3, "w2": 1 / 3, "L": 9 / 14,
                       "clipped_share": clip, "mean_abs_log_factor": size},
            "after": {"corrected": False}}}}

    one = {"task": {"key": "T2e/a"}, "build_records": {"fgw_resid": {"0": rec(0.0, 0.2), "1": rec(0.1, 0.4)}}}
    (tmp_path / "T2e_a.json").write_text(json.dumps(one))
    (tmp_path / "T2e_b.json").write_text(json.dumps({"task": {"key": "T2e/b"}}))
    got = n49.corrections([tmp_path])
    assert got["T2e/b"] == {} and set(got["T2e/a"]) == {"before"}
    b = got["T2e/a"]["before"]
    assert b["neighbours"] == ["E7.25", "E8.0"] and b["L"] == pytest.approx(9 / 14) and b["seeds"] == 2
    assert b["clipped_share"] == pytest.approx(0.05) and b["mean_abs_log_factor"] == pytest.approx(0.3)


def test_the_reading_uses_the_committed_arms_and_the_runs_arm(n49, monkeypatch):
    scores, tasks = tables(n49, 1.0)
    committed, run = split(n49, scores)
    run[tasks["new1"].key]["fgw"] += 0.005  # a rescored arm inside the tolerance: the record's is used
    rec = {tasks["e75"].key: {"before": {"neighbours": ["E6.75", "E8.0"], "L": 0.658, "seeds": 10,
                                         "clipped_share": 0.012, "mean_abs_log_factor": 0.2}}}
    patch(n49, monkeypatch, committed, run, tasks, rec)
    res = n49.reading([RUN])
    assert not res["void"] and res["runs"] == ["committed:0000000", "run:1111111"]
    assert res["rescored_max_abs_difference"]["fgw"] == pytest.approx(0.005)
    assert res["own_minus_fgw"]["new1"] == pytest.approx(1.0) and res["first_passes"]
    assert res["corrections"]["e75"]["before"]["L"] == 0.658 and res["corrections"]["half"] == {}
    txt = n49.text(res)
    assert "fgw_resid is put forward as a first upload (to a recomputation and a cold check" in txt
    assert "e75 before 0.658 0.012 0.200" in txt and "half none" in txt
    assert "first upload beside fgw_aligned_series (row N33's rule): S(0.45) +0.55" in txt
    json.dumps(res, default=float)


def test_the_reading_is_void_on_a_rescored_mismatch_a_missing_arm_or_a_missing_task(n49, monkeypatch):
    scores, tasks = tables(n49, 1.0)
    committed, run = split(n49, scores)
    run[tasks["late"].key]["fgw"][5] -= 0.3
    patch(n49, monkeypatch, committed, run, tasks)
    res = n49.reading([RUN])
    assert res["void"] and any("late: fgw differs from its committed scores" in w for w in res["why"])
    assert n49.text(res).startswith("VOID")
    committed, run = split(n49, scores)
    del committed[tasks["x3"].key]["fgw_trend_render"]
    patch(n49, monkeypatch, committed, run, tasks)
    res = n49.reading([RUN])
    assert res["void"] and any("x3: no committed scores for fgw_trend_render" in w for w in res["why"])
    committed, run = split(n49, scores)
    del run[tasks["x1"].key]
    patch(n49, monkeypatch, committed, run, tasks)
    res = n49.reading([RUN])
    assert res["void"] and any("x1: no task file in the run" in w for w in res["why"])
    scores, tasks = tables(n49, 0.0)
    committed, run = split(n49, scores)
    patch(n49, monkeypatch, committed, run, tasks)
    res = n49.reading([RUN])
    assert not res["void"] and "fgw_resid is not put forward: closed by this reading" in n49.text(res)
