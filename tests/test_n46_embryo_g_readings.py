"""Readings of tries E42 (experiments/harness/n46_embryo_g_readings.py): the controls for fgw_partial, on synthetic
score tables and movement records built here. Nothing is read from data/ or from a committed run, and nothing is
scored."""

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
def n46():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n46_embryo_g_readings", HARNESS / "n46_embryo_g_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n46, gains: dict[str, float | dict[str, float | None]] | None = None):
    """Per-seed scores on the nine tasks. The floor at 50; fgw at 60 plus a per-seed ramp; the second upload as built
    at 63 on the DEEP tasks and 58 on the others. Each of the five arms at fgw + its gain: one number for every task,
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


def records(n46, run: dict, partial=(0.5, 0.6), s=0.7) -> dict:
    """A movement record for every arm the run scored, every seed: fgw_partial's cells move `partial` of the entry's
    distance (earlier bracket, later bracket) and the matched control's exactly as far; the damped control s; the
    random control about s; the still control 0."""
    cell = {n46.PARTIAL: list(partial), n46.MATCHED: list(partial), n46.DAMP: [s, s], n46.RAND: [s + 0.01, s - 0.01],
            n46.STILL: [0.0, 0.0]}
    anchor = {**cell, n46.PARTIAL: [partial[0] - 0.1, partial[1] - 0.1]}
    return {k: {a: {seed: {"cell": list(cell[a]), "anchor": list(anchor[a])} for seed in range(10)}
                for a in n46.MOVED if a in arms} for k, arms in run.items()}


def patch(n46, monkeypatch, committed: dict, run: dict, tasks: dict, moves: dict | None = None) -> None:
    monkeypatch.setattr(n46.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n46.nt, "_load", lambda dirs: ((run, [], ["run:1111111"]) if list(dirs) == [RUN]
                                                      else (committed, [], ["committed:0000000"])))
    monkeypatch.setattr(n46, "movement", lambda dirs: records(n46, run) if moves is None else moves)


def test_the_constants_are_the_rows(n46):
    assert n46.CONTROLS == ("fgw_still", "fgw_damp", "fgw_part_rand", "fgw_damp_m") and n46.PARTIAL == "fgw_partial"
    assert n46.MATCHED == "fgw_damp_m"
    assert n46.ORDER == ("fgw_still", "fgw_damp", "fgw_part_rand", "fgw_damp_m", "fgw_partial")
    assert n46.MOVED == ("fgw_partial", "fgw_still", "fgw_damp", "fgw_part_rand", "fgw_damp_m")
    assert n46.RESCORED == ("floor", "fgw", "fgw_partial")
    assert n46.ENTRY == "fgw" and n46.SECOND == "fgw_aligned_series" and n46.BACKTESTS == ("e75", "e725")
    assert n46.N45_RUNS.parent.name == "n45_embryo_f_arms"
    assert n46.MATCH_TOL == 1e-6 and n46.MOVEMENT_KEY == "movement_ratio"


def test_partial_ahead_of_every_control_is_matching_and_is_put_forward(n46):
    # the arm gains 2 everywhere and every control gains 1: it leads each of the four by 1 on both backtests
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": 1.0, "fgw_part_rand": 1.0, "fgw_damp_m": 1.0,
                                 "fgw_still": 1.0})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "matching" and res["put_forward"] == "fgw_partial" and not res["closed"]
    assert all(res["lead"][c]["keep"] for c in n46.CONTROLS)
    assert res["passes"]["fgw_still"], "a passing control does not displace the arm when the verdict is matching"


@pytest.mark.parametrize("control", ["fgw_still", "fgw_damp", "fgw_part_rand", "fgw_damp_m"])
def test_one_control_that_is_not_led_makes_the_verdict_moving_less(n46, control):
    # every control gains 1 except one, which gains 1.9: the arm's lead of 0.1 over it is under the rule's bar
    gains = {"fgw_partial": 2.0, **{c: 1.0 for c in n46.CONTROLS}, control: 1.9}
    res = n46.check(*tables(n46, gains))
    assert res["verdict"] == "moving less" and not res["lead"][control]["keep"]
    assert all(res["lead"][c]["keep"] for c in n46.CONTROLS if c != control)
    assert res["put_forward"] == "fgw_still", "every arm passes here, so the first of the order is put forward"


def test_a_still_cloud_that_does_as_well_is_put_forward(n46):
    # a gain g on every task gives the pair rule S(0.45) = 0.55 g (an arm leads the pair on the three OTHER tasks
    # only), so 2.5 passes it; the arm does not lead a control that scores above it
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_still": 2.5, "fgw_damp": 1.0, "fgw_part_rand": 1.0,
                                 "fgw_damp_m": 1.0})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less" and not res["lead"]["fgw_still"]["keep"]
    assert res["passes"]["fgw_still"] and res["put_forward"] == "fgw_still"


def test_moving_less_puts_forward_the_first_passing_arm_of_the_order(n46):
    # 0.6 passes the pair rule (S(0.45) = 0.33) and 0.5 would not
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": 1.9, "fgw_part_rand": 0.6, "fgw_damp_m": 0.6,
                                 "fgw_still": 0.0})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less"
    assert not res["lead"]["fgw_damp"]["keep"] and res["lead"]["fgw_part_rand"]["keep"]
    assert res["passes"] == {"fgw_still": False, "fgw_damp": True, "fgw_part_rand": True, "fgw_damp_m": True,
                             "fgw_partial": True}
    assert res["put_forward"] == "fgw_damp"
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp_m": 1.9})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less" and res["put_forward"] == "fgw_damp_m"
    assert res["passes"] == {"fgw_still": False, "fgw_damp": False, "fgw_part_rand": False, "fgw_damp_m": True,
                             "fgw_partial": True}


def test_moving_less_with_no_passing_control_leaves_partial_as_the_arm(n46):
    # the controls gain 0.52 on every task: ahead of fgw on both backtests, but S(0.45) = 0.286 is under the pair
    # rule's bar; the arm gains 0.8 and passes both rules, and its lead of 0.28 over each control is under the bar
    scores, tasks = tables(n46, {"fgw_partial": 0.8, **{c: 0.52 for c in n46.CONTROLS}})
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less"
    assert res["passes"] == {**{c: False for c in n46.CONTROLS}, "fgw_partial": True}
    assert res["put_forward"] == "fgw_partial"


def test_nothing_passes_closes_the_line(n46):
    scores, tasks = tables(n46)
    res = n46.check(scores, tasks)
    assert res["verdict"] == "moving less" and res["put_forward"] is None and res["closed"]
    assert res["own_minus_fgw"]["fgw_still"]["x2"] == pytest.approx(0.0)


def test_an_arm_needs_both_rules_to_pass(n46):
    # a gain on the DEEP tasks alone passes neither rule; a gain on both backtests with a loss on the third OTHER
    # task passes the adoption rule and fails the pair rule
    scores, tasks = tables(n46, {"fgw_still": {n: 1.0 for n in n46.nt.DEEP}})
    res = n46.check(scores, tasks)
    assert not res["adoption"]["fgw_still"]["keep"] and not res["passes"]["fgw_still"]
    scores, tasks = tables(n46, {"fgw_damp_m": {"e75": 1.0, "e725": 1.0, "new2": -3.0}})
    res = n46.check(scores, tasks)
    assert res["adoption"]["fgw_damp_m"]["keep"]
    assert not res["first_upload"]["first_uploads"]["fgw_damp_m"]["passes"] and not res["passes"]["fgw_damp_m"]


def test_a_control_missing_on_a_task_is_not_readable_and_cannot_lead_or_pass(n46):
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": {"e75": None}, "fgw_part_rand": 1.0,
                                 "fgw_damp_m": 1.0, "fgw_still": 1.0})
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


def test_the_movement_record_is_read_from_the_task_files(n46, tmp_path):
    rec = {"s": 0.5, "movement_ratio": {"anchor": [0.4, 0.5], "cell": [0.45, 0.55]}}
    one = {"task": {"key": "T2e/a"}, "build_records": {"fgw_partial": {"0": rec, "1": {"s": 0.5}},
                                                       "fgw_potts": {"0": {"potts": {}}}}}
    (tmp_path / "T2e_a.json").write_text(json.dumps(one))
    (tmp_path / "T2e_b.json").write_text(json.dumps({"task": {"key": "T2e/b"}}))
    (tmp_path / "other.json").write_text("not read")
    got = n46.movement([tmp_path])
    assert set(got) == {"T2e/a", "T2e/b"} and set(got["T2e/a"]) == set(n46.MOVED)
    assert got["T2e/a"]["fgw_partial"] == {0: {"anchor": [0.4, 0.5], "cell": [0.45, 0.55]}, 1: None}
    assert got["T2e/a"]["fgw_damp_m"] == {} and got["T2e/b"]["fgw_still"] == {}


def test_the_matched_control_must_move_its_cells_as_far_as_the_arm(n46):
    scores, tasks = tables(n46, {"fgw_partial": 1.0})
    _, run = split(n46, scores)
    moves = records(n46, run, partial=(0.5, 0.6), s=0.7)
    table, worst, why = n46.moved(moves, run, tasks)
    assert why == [] and worst == 0.0
    assert table["fgw_partial"]["e75"] == {"cell": [pytest.approx(0.5), pytest.approx(0.6)],
                                           "anchor": [pytest.approx(0.4), pytest.approx(0.5)]}
    assert table["fgw_damp"]["x1"]["cell"] == [pytest.approx(0.7), pytest.approx(0.7)]
    assert table["fgw_still"]["half"]["cell"] == [0.0, 0.0]
    moves[tasks["late"].key]["fgw_damp_m"][7]["cell"][1] += 1e-3  # one seed, one bracket, off the arm's
    _, worst, why = n46.moved(moves, run, tasks)
    assert worst == pytest.approx(1e-3)
    assert why == ["late: fgw_damp_m moves its cells 1.00e-03 off fgw_partial (ratio to the entry's)"]
    moves[tasks["late"].key]["fgw_damp_m"][7]["cell"][1] -= 1e-3 - 1e-7  # inside the tolerance
    _, worst, why = n46.moved(moves, run, tasks)
    assert why == [] and worst == pytest.approx(1e-7, rel=1e-3)


def test_an_arm_the_run_scored_needs_its_movement_record_for_every_seed(n46):
    scores, tasks = tables(n46, {"fgw_partial": 1.0})
    _, run = split(n46, scores)
    moves = records(n46, run)
    del moves[tasks["x2"].key]["fgw_damp"][4]
    moves[tasks["x3"].key]["fgw_partial"][2] = None
    moves[tasks["new1"].key]["fgw_still"][0]["cell"] = [0.0]
    moves[tasks["new2"].key]["fgw_part_rand"][9]["anchor"] = [None, 0.7]
    table, _, why = n46.moved(moves, run, tasks)
    assert why == ["new1: no movement record for fgw_still (seeds [0])",
                   "x2: no movement record for fgw_damp (seeds [4])",
                   "x3: no movement record for fgw_partial (seeds [2])",
                   "new2: no movement record for fgw_part_rand (seeds [9])"]
    assert "x2" not in table["fgw_damp"] and "x2" in table["fgw_still"]
    # an arm the run did not score on a task (its build was refused) owes no record there
    _, run = split(n46, scores)
    del run[tasks["e75"].key]["fgw_damp_m"]
    moves = records(n46, run)
    table, worst, why = n46.moved(moves, run, tasks)
    assert why == [] and worst == 0.0 and "e75" not in table["fgw_damp_m"]


def test_the_reading_uses_the_committed_arms_and_the_runs_controls(n46, monkeypatch):
    scores, tasks = tables(n46, {"fgw_partial": 2.0, "fgw_damp": 1.9, "fgw_part_rand": 1.9, "fgw_damp_m": 1.9,
                                 "fgw_still": 0.1})
    committed, run = split(n46, scores)
    run[tasks["new1"].key]["fgw_partial"] += 0.005  # a rescored arm inside the tolerance: the record's is used
    patch(n46, monkeypatch, committed, run, tasks)
    res = n46.reading([RUN])
    assert not res["void"] and res["runs"] == ["committed:0000000", "run:1111111"]
    assert res["rescored_max_abs_difference"]["fgw_partial"] == pytest.approx(0.005)
    assert res["matched_max_abs_difference"] == 0.0
    assert res["own_minus_fgw"]["fgw_partial"]["new1"] == pytest.approx(2.0)
    assert res["verdict"] == "moving less" and res["put_forward"] == "fgw_damp"
    txt = n46.text(res)
    assert "verdict: moving less; arm put forward: fgw_damp (to a cold check before any request)" in txt
    assert "fgw_partial against fgw_damp:" in txt and "fgw_partial against fgw_damp_m:" in txt
    assert "does not lead" in txt
    assert "Readings of tries E42a to E42d" in txt and "rescored arms against their committed scores" in txt
    assert "fgw_partial 0.50 0.60, fgw_still 0.00 0.00, fgw_damp 0.70 0.70" in txt
    assert "fgw_damp_m against fgw_partial, largest difference of that share" in txt
    json.dumps(res, default=float)


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


def test_the_reading_is_void_without_the_movement_record_or_with_an_unmatched_control(n46, monkeypatch):
    scores, tasks = tables(n46, {"fgw_partial": 1.0})
    committed, run = split(n46, scores)
    moves = records(n46, run)
    del moves[tasks["half"].key]["fgw_partial"]
    patch(n46, monkeypatch, committed, run, tasks, moves)
    res = n46.reading([RUN])
    assert res["void"] and any("half: no movement record for fgw_partial" in w for w in res["why"])
    moves = records(n46, run)
    moves[tasks["e725"].key]["fgw_damp_m"][3]["cell"][0] = 0.9
    patch(n46, monkeypatch, committed, run, tasks, moves)
    res = n46.reading([RUN])
    assert res["void"] and any("e725: fgw_damp_m moves its cells 4.00e-01 off fgw_partial" in w for w in res["why"])
    assert res["matched_max_abs_difference"] == pytest.approx(0.4)
    assert n46.text(res).startswith("VOID")
