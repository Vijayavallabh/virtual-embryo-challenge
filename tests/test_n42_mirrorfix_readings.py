"""Readings of tries E38 (experiments/harness/n42_mirrorfix_readings.py): three existing rules applied to four new
arms, on synthetic score tables built here. Nothing is read from data/ or from a committed run, and nothing is
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
def n42():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n42_mirrorfix_readings", HARNESS / "n42_mirrorfix_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n42, gains: dict[str, dict[str, float | None]] | None = None, spread: float = 0.0):
    """Per-seed scores on the nine tasks: fgw at 60; the second upload as built at 63 on the DEEP tasks and 58 on the
    others; fgw_series at 61 and 57; each new arm at its base + gains[arm][task] (default 0; None = no scores).
    spread adds a per-seed ramp to every arm alike, which leaves paired gains as set."""
    nt = n42.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    gains = gains or {}
    ramp = spread * np.arange(10)
    scores = {}
    for n in names:
        deep = n in nt.DEEP
        s = {n42.ENTRY: 60.0 + ramp, n42.SECOND: (63.0 if deep else 58.0) + ramp,
             "fgw_series": (61.0 if deep else 57.0) + ramp}
        for arm in n42.ARMS:
            v = gains.get(arm, {}).get(n, 0.0)
            if v is not None:
                s[arm] = s[n42.BASE_OF[arm]] + v
        scores[tasks[n].key] = s
    return scores, tasks


def test_the_constants_are_the_rows(n42):
    assert n42.ARMS == ("fgw_mirrorfix", "fgw_mirrorfix_series", "fgw_mirrorgate", "fgw_mirrorgate_series")
    assert n42.ADOPTION_ARMS == ("fgw_mirrorfix", "fgw_mirrorgate")
    assert n42.BASE_OF == {"fgw_mirrorfix": "fgw", "fgw_mirrorfix_series": "fgw_series", "fgw_mirrorgate": "fgw",
                           "fgw_mirrorgate_series": "fgw_series"}
    assert n42.ENTRY == "fgw" and n42.SECOND == "fgw_aligned_series" and n42.BACKTESTS == ("e75", "e725")


def test_arms_equal_to_their_bases_close_the_line(n42):
    scores, tasks = tables(n42, spread=0.1)
    res = n42.check(scores, tasks)
    assert res["closed"] and res["adopted"] == []
    assert res["requested_first_upload"] is None and res["requested_second_upload"] is None
    assert all(res["adoption"][a]["tasks"]["e725"]["mean_gain"] == pytest.approx(0.0) for a in n42.ADOPTION_ARMS)
    assert res["own_minus_base"]["fgw_mirrorfix_series"]["x1"] == pytest.approx(0.0)


def test_a_gain_on_both_backtests_passes_the_adoption_rule(n42):
    scores, tasks = tables(n42, {"fgw_mirrorgate": {"e75": 0.5, "e725": 0.6}})
    res = n42.check(scores, tasks)
    a = res["adoption"]["fgw_mirrorgate"]
    assert a["keep"] and a["tasks"]["e75"]["seeds_ahead"] == 10
    assert a["tasks"]["e725"]["mean_gain"] == pytest.approx(0.6)
    assert not res["adoption"]["fgw_mirrorfix"]["keep"]
    assert res["adopted"] == ["fgw_mirrorgate"] and not res["closed"]


def test_a_gain_on_one_backtest_and_a_loss_on_the_other_is_not_adopted(n42):
    scores, tasks = tables(n42, {"fgw_mirrorfix": {"e75": 2.0, "e725": -0.5}})
    res = n42.check(scores, tasks)
    a = res["adoption"]["fgw_mirrorfix"]
    assert not a["keep"] and a["tasks"]["e725"]["vetoes"] and res["adopted"] == []


def test_the_first_upload_rule_reads_the_pair_beside_the_second_upload(n42):
    scores, tasks = tables(n42, {"fgw_mirrorfix": {n: 1.0 for n in n42.nt.OTHER}})
    res = n42.check(scores, tasks)
    v = res["first_upload"]["first_uploads"]["fgw_mirrorfix"]
    assert v["O"] == pytest.approx(1.0) and v["D"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_first_upload"] == "fgw_mirrorfix" and not res["closed"]
    scores, tasks = tables(n42, {"fgw_mirrorfix": {n: 1.0 for n in n42.nt.DEEP}})  # hidden behind the second upload
    res = n42.check(scores, tasks)
    assert res["requested_first_upload"] is None


def test_the_second_upload_rule_reads_the_pair_beside_fgw_against_the_second_upload_as_built(n42):
    # fgw_series trails the second upload by 2 on the DEEP tasks; the correction adds 3 there: 1 above it
    scores, tasks = tables(n42, {"fgw_mirrorfix_series": {n: 3.0 for n in n42.nt.DEEP}})
    res = n42.check(scores, tasks)
    v = res["second_upload"]["second_uploads"]["fgw_mirrorfix_series"]
    assert v["D"] == pytest.approx(1.0) and v["O"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_second_upload"] == "fgw_mirrorfix_series" and not res["closed"]
    assert res["own_minus_base"]["fgw_mirrorfix_series"]["half"] == pytest.approx(3.0)
    scores, tasks = tables(n42, {"fgw_mirrorfix_series": {n: 1.0 for n in n42.nt.DEEP}})  # still under it
    res = n42.check(scores, tasks)
    assert res["requested_second_upload"] is None and res["closed"]


def test_an_arm_missing_on_a_task_is_not_readable_and_cannot_pass(n42):
    gains = {"fgw_mirrorfix_series": {**{n: 5.0 for n in n42.nt.DEEP}, "x3": None},
             "fgw_mirrorfix": {"e75": None, "e725": 5.0}}
    scores, tasks = tables(n42, gains)
    res = n42.check(scores, tasks)
    assert not res["readable"]["fgw_mirrorfix_series"] and "fgw_mirrorfix_series" not in res["own_minus_base"]
    assert not res["second_upload"]["second_uploads"]["fgw_mirrorfix_series"]["readable"]
    assert res["adoption"]["fgw_mirrorfix"] is None and res["adoption"]["fgw_mirrorgate"] is not None
    assert res["adopted"] == [] and res["closed"]


def test_the_reading_is_void_without_a_base_and_the_text_states_the_outcome(n42, monkeypatch):
    scores, tasks = tables(n42)
    monkeypatch.setattr(n42.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n42.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n42.reading([])
    assert not res["void"] and res["closed"]
    txt = n42.text(res)
    assert "closed: no rule is passed" in txt and "adoption rule, fgw_mirrorfix against fgw" in txt
    assert "fgw_mirrorfix_series minus fgw_series" in txt and "Readings of tries E38a to E38d" in txt
    assert "adoption rule, fgw_mirrorgate against fgw" in txt
    scores2, _ = tables(n42, {"fgw_mirrorfix": {"e75": 0.5, "e725": 0.6}})
    monkeypatch.setattr(n42.nt, "_load", lambda dirs: (scores2, [], ["synthetic:0000000"]))
    assert "the adoption rule is passed by fgw_mirrorfix" in n42.text(n42.reading([]))
    del scores["T2e/half"]["fgw_series"]
    monkeypatch.setattr(n42.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n42.reading([])
    assert res["void"] and any("no scores for fgw_series" in w for w in res["why"])
    assert n42.text(res).startswith("VOID")


def test_the_count_of_seeds_equal_to_the_base_shows_where_a_gate_fired(n42):
    scores, tasks = tables(n42)
    key = tasks["e725"].key
    gated = scores[key]["fgw_mirrorgate"].copy()
    gated[[1, 4, 7]] += 0.8  # the gate fires in three seeds: there the arm is no longer its base
    gated[2] += 0.004  # a difference inside the tolerance between runs still counts as the base
    scores[key]["fgw_mirrorgate"] = gated
    res = n42.check(scores, tasks)
    assert res["seeds_equal_to_base"]["fgw_mirrorgate"]["e725"] == 7
    assert res["seeds_equal_to_base"]["fgw_mirrorgate"]["half"] == 10
    assert "seeds (of 10) in which fgw_mirrorgate scores as fgw:" in n42.text({"void": False, "runs": ["x:0"], **res})
