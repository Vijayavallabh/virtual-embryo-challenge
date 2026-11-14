"""Reading of try E45 (experiments/harness/n50_embryo_j_readings.py): the lineage-level twin of the movement-matched
control, on synthetic score tables and movement records built here. Nothing is read from data/ or from a committed
run, and nothing is scored."""

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
def n50():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n50_embryo_j_readings", HARNESS / "n50_embryo_j_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n50, arm: float | dict | None = 0.0, damp_m: float = 1.0, grp: float = 1.5):
    """Per-seed scores on the nine tasks. The floor at 50; fgw at 60 plus a per-seed ramp; the second upload as built
    at 63 on the DEEP tasks and 58 on the others; fgw_damp_m and fgw_partial_grp at fgw + their gains; the twin at
    fgw + its gain (one number for every task, a dict per task whose values may be per-seed arrays, or None for no
    scores at all)."""
    nt = n50.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    ramp = 0.1 * np.arange(10)
    scores = {}
    for n in names:
        s = {"floor": 50.0 + 0 * ramp, n50.ENTRY: 60.0 + ramp, n50.SECOND: (63.0 if n in nt.DEEP else 58.0) + ramp}
        s[n50.DAMP_M], s[n50.GRP] = s[n50.ENTRY] + damp_m, s[n50.ENTRY] + grp
        v = arm.get(n, 0.0) if isinstance(arm, dict) else arm
        if v is not None:
            s[n50.ARM] = s[n50.ENTRY] + v
        scores[tasks[n].key] = s
    return scores, tasks


def split(n50, scores: dict) -> tuple[dict, dict]:
    committed = {k: {a: v for a, v in s.items() if a != n50.ARM} for k, s in scores.items()}
    run = {k: {a: v.copy() for a, v in s.items() if a == n50.ARM or a in n50.RESCORED} for k, s in scores.items()}
    return committed, run


def records(n50, run: dict, grp=(0.8, 0.85), cell=(0.5, 0.6)) -> dict:
    """What the run records, every seed: fgw_partial_grp's cells move `grp` of the entry's distance and its anchors
    less; the twin's anchors and cells move by its two factors, which are those ratios; fgw_damp_m's by `cell`."""
    def rec(c, a, f):
        return {"ratio": {"cell": list(c), "anchor": list(a)}, "factors": None if f is None else list(f)}

    per = {n50.GRP: rec(grp, (grp[0] - 0.1, grp[1] - 0.1), None), n50.ARM: rec(grp, grp, grp),
           n50.DAMP_M: rec(cell, cell, cell)}
    return {k: {a: {seed: json.loads(json.dumps(per[a])) for seed in range(10)} for a in n50.MOVED if a in arms}
            for k, arms in run.items()}


def patch(n50, monkeypatch, committed: dict, run: dict, tasks: dict, moves: dict | None = None) -> None:
    monkeypatch.setattr(n50.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n50.nt, "_load", lambda dirs: ((run, [], ["run:1111111"]) if list(dirs) == [RUN]
                                                      else (committed, [], ["committed:0000000"])))
    monkeypatch.setattr(n50, "movement", lambda dirs: records(n50, run) if moves is None else moves)


def test_the_constants_are_the_rows(n50):
    assert n50.ARM == "fgw_damp_m_grp" and n50.DAMP_M == "fgw_damp_m" and n50.GRP == "fgw_partial_grp"
    assert n50.RESCORED == ("floor", "fgw", "fgw_partial_grp", "fgw_damp_m")
    assert n50.MOVED == ("fgw_partial_grp", "fgw_damp_m", "fgw_damp_m_grp")
    assert n50.BACKTESTS == ("e75", "e725") and n50.MIRRORED == ("e725", (3, 4, 5, 8, 9))
    assert n50.N46_RUNS.parent.name == "n46_embryo_g_arms" and n50.N47_RUNS.parent.name == "n47_embryo_h_arms"
    assert n50.MATCH_TOL == 1e-6


def test_a_twin_that_gains_everywhere_passes_both_rules(n50):
    # a gain g on every task gives the pair rule S(0.45) = 0.55 g, so 1.0 passes it and both backtests
    res = n50.check(*tables(n50, 1.0))
    assert res["adoption"]["keep"] and res["first_upload"]["first_uploads"]["fgw_damp_m_grp"]["passes"]
    assert res["passes"] and not res["closed"]
    assert res["own_minus_fgw"]["x2"] == pytest.approx(1.0) and res["own_minus_damp_m"]["x2"] == pytest.approx(0.0)
    assert res["own_minus_partial_grp"]["x2"] == pytest.approx(-0.5)


def test_a_twin_needs_both_rules(n50):
    res = n50.check(*tables(n50, 0.3))
    assert not res["adoption"]["keep"] and not res["passes"] and res["closed"], "a gain of exactly 0.3 does not pass"
    res = n50.check(*tables(n50, 0.52))
    assert res["adoption"]["keep"] and not res["passes"], "S(0.45) = 0.286 is under the pair rule's bar"
    res = n50.check(*tables(n50, {"e75": 1.0, "e725": 1.0, "new2": -3.0}))
    assert res["adoption"]["keep"] and not res["passes"]
    res = n50.check(*tables(n50, {n: 2.0 for n in n50.nt.DEEP}))
    assert not res["adoption"]["keep"] and res["closed"]


def test_the_mirrored_seeds_gain_is_reported_apart(n50):
    gain = np.zeros(10)
    gain[[3, 4, 5, 8, 9]] = 2.0  # the whole gain on the mirrored seeds: a mean of 1.0 with five seeds ahead
    res = n50.check(*tables(n50, {"e725": gain, "e75": 1.0, "new2": 1.0}))
    assert res["mirrored"] == {"task": "e725", "seeds": [3, 4, 5, 8, 9], "gain_on_them": pytest.approx(2.0),
                               "gain_on_the_others": pytest.approx(0.0)}
    assert not res["adoption"]["keep"], "five seeds ahead are not eight"


def test_the_rescored_arms_must_match_the_record(n50):
    scores, tasks = tables(n50, 1.0)
    committed, run = split(n50, scores)
    worst, why = n50.rescored(committed, run, tasks)
    assert why == [] and set(worst) == set(n50.RESCORED) and max(worst.values()) == 0.0
    run[tasks["e725"].key]["fgw_damp_m"][0] += 0.02
    worst, why = n50.rescored(committed, run, tasks)
    assert worst["fgw_damp_m"] == pytest.approx(0.02)
    assert why == ["e725: fgw_damp_m differs from its committed scores by 2.00e-02"]
    del run[tasks["half"].key]["fgw_partial_grp"]
    _, why = n50.rescored(committed, run, tasks)
    assert "half: the run did not rescore fgw_partial_grp" in why


def test_the_movement_record_is_read_from_the_task_files(n50, tmp_path):
    ratio = {"anchor": [0.7, 0.75], "cell": [0.8, 0.85]}
    one = {"task": {"key": "T2e/a"},
           "build_records": {"fgw_partial_grp": {"0": {"partial": {"s": 0.6}, "movement_ratio": ratio}},
                             "fgw_damp_m_grp": {"0": {"factors": [0.8, 0.85], "movement_ratio": ratio}, "1": {}}}}
    (tmp_path / "T2e_a.json").write_text(json.dumps(one))
    got = n50.movement([tmp_path])
    assert set(got["T2e/a"]) == set(n50.MOVED) and got["T2e/a"]["fgw_damp_m"] == {}
    assert got["T2e/a"]["fgw_partial_grp"] == {0: {"ratio": ratio, "factors": None}}
    assert got["T2e/a"]["fgw_damp_m_grp"] == {0: {"ratio": ratio, "factors": [0.8, 0.85]},
                                              1: {"ratio": None, "factors": None}}


def test_the_twin_must_move_its_cells_as_far_as_the_lineage_level_arm(n50):
    scores, tasks = tables(n50, 1.0)
    _, run = split(n50, scores)
    moves = records(n50, run)
    table, factors, worst, why = n50.moved(moves, run, tasks)
    assert why == [] and worst == 0.0
    assert table["fgw_damp_m_grp"]["e75"] == [pytest.approx(0.8), pytest.approx(0.85)]
    assert table["fgw_damp_m"]["e75"] == [pytest.approx(0.5), pytest.approx(0.6)]
    assert factors["fgw_damp_m_grp"]["late"] == [[0.8, 0.8], [0.85, 0.85]]
    assert factors["fgw_damp_m"]["late"] == [[0.5, 0.5], [0.6, 0.6]]
    moves[tasks["late"].key]["fgw_partial_grp"][7]["ratio"]["cell"][1] += 1e-3
    _, _, worst, why = n50.moved(moves, run, tasks)
    assert worst == pytest.approx(1e-3)
    assert why == ["late: fgw_damp_m_grp moves its cells 1.00e-03 off fgw_partial_grp (ratio to the entry's)"]


def test_a_damped_arms_record_must_equal_its_factors(n50):
    scores, tasks = tables(n50, 1.0)
    _, run = split(n50, scores)
    k = tasks["x2"].key
    moves = records(n50, run)
    moves[k]["fgw_damp_m_grp"][3]["ratio"]["anchor"] = [0.7, 0.75]  # the lineage-level arm's own record under its name
    moves[k]["fgw_damp_m"][0]["factors"] = None
    _, factors, _, why = n50.moved(moves, run, tasks)
    assert why == ["x2: fgw_damp_m's movement record holds no factors",
                   "x2: fgw_damp_m_grp's movement record is 1.00e-01 off its factors"]
    assert "x2" not in factors["fgw_damp_m"] and "x2" in factors["fgw_damp_m_grp"]
    moves = records(n50, run)
    del moves[k]["fgw_partial_grp"][4]
    moves[k]["fgw_damp_m_grp"][9]["ratio"]["cell"] = [None, 0.85]
    _, _, _, why = n50.moved(moves, run, tasks)
    assert why == ["x2: no movement record for fgw_partial_grp (seeds [4])",
                   "x2: no movement record for fgw_damp_m_grp (seeds [9])"]


def test_the_reading_uses_the_committed_arms_and_the_runs_twin(n50, monkeypatch):
    scores, tasks = tables(n50, 1.0)
    committed, run = split(n50, scores)
    run[tasks["new1"].key]["fgw_damp_m"] += 0.005  # a rescored arm inside the tolerance: the record's is used
    patch(n50, monkeypatch, committed, run, tasks)
    res = n50.reading([RUN])
    assert not res["void"] and res["runs"] == ["committed:0000000", "run:1111111"]
    assert res["rescored_max_abs_difference"]["fgw_damp_m"] == pytest.approx(0.005)
    assert res["matched_max_abs_difference"] == 0.0 and res["passes"]
    assert res["own_minus_damp_m"]["new1"] == pytest.approx(0.0)
    txt = n50.text(res)
    assert "fgw_damp_m_grp passes both rules: it goes to a recomputation and a cold check before any request" in txt
    assert "fgw_damp_m_grp's factors, earlier and later bracket, lowest to highest over seeds: half 0.80-0.80" in txt
    assert ("on e725, the twin's gain over fgw on the mirrored seeds [3, 4, 5, 8, 9]: +1.00; on the other five: "
            "+1.00") in txt
    json.dumps(res, default=float)
    scores, tasks = tables(n50, 0.0)
    committed, run = split(n50, scores)
    patch(n50, monkeypatch, committed, run, tasks)
    res = n50.reading([RUN])
    assert not res["void"] and res["closed"]
    assert "does not pass both rules: the line closes with no request (rows E42 and E43)" in n50.text(res)


def test_the_reading_is_void_on_a_mismatch_a_refused_twin_or_a_broken_record(n50, monkeypatch):
    scores, tasks = tables(n50, 1.0)
    committed, run = split(n50, scores)
    run[tasks["late"].key]["fgw"][5] -= 0.3
    patch(n50, monkeypatch, committed, run, tasks)
    res = n50.reading([RUN])
    assert res["void"] and any("late: fgw differs from its committed scores" in w for w in res["why"])
    assert n50.text(res).startswith("VOID")
    committed, run = split(n50, scores)
    del run[tasks["e75"].key]["fgw_damp_m_grp"]
    patch(n50, monkeypatch, committed, run, tasks)
    res = n50.reading([RUN])
    assert res["void"] and any(w.startswith("e75: the run did not score fgw_damp_m_grp") for w in res["why"])
    committed, run = split(n50, scores)
    del run[tasks["x1"].key]
    patch(n50, monkeypatch, committed, run, tasks)
    res = n50.reading([RUN])
    assert res["void"] and any("x1: no task file in the run" in w for w in res["why"])
    committed, run = split(n50, scores)
    del committed[tasks["x3"].key]["fgw_partial_grp"]
    patch(n50, monkeypatch, committed, run, tasks)
    res = n50.reading([RUN])
    assert res["void"] and any("x3: no committed scores for fgw_partial_grp" in w for w in res["why"])
    committed, run = split(n50, scores)
    moves = records(n50, run)
    moves[tasks["e725"].key]["fgw_partial_grp"][3]["ratio"]["cell"][0] = 0.4
    patch(n50, monkeypatch, committed, run, tasks, moves)
    res = n50.reading([RUN])
    assert res["void"] and res["matched_max_abs_difference"] == pytest.approx(0.4)
    patch(n50, monkeypatch, {}, {}, {})
    res = n50.reading([RUN])
    assert res["void"] and res["why"] == ["fgw_damp_m_grp was never compared with fgw_partial_grp"]
