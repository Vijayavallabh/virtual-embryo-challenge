"""Readings of tries E36 (experiments/harness/n41_polar_readings.py): three existing rules applied to two new arms, on
synthetic score tables built here. Nothing is read from data/ or from a committed run, and nothing is scored."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"


@pytest.fixture(scope="module")
def n41():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n41_polar_readings", HARNESS / "n41_polar_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n41, polar: dict[str, float] | None = None, as_polar: dict[str, float | None] | None = None,
           spread: float = 0.0):
    """Per-seed scores on the nine tasks: fgw at 60; the second upload as built at 63 on the DEEP tasks and 58 on the
    others; fgw_polar at fgw + polar[task] (default 0); fgw_aligned_polar_series at the second upload + as_polar[task]
    (default 0; None = no scores). spread adds a per-seed ramp to every arm alike, which leaves paired gains as set."""
    nt = n41.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    polar, as_polar = polar or {}, as_polar or {}
    ramp = spread * np.arange(10)
    scores = {}
    for n in names:
        second = 63.0 if n in nt.DEEP else 58.0
        s = {n41.ENTRY: 60.0 + ramp, n41.SECOND: second + ramp, "fgw_polar": 60.0 + polar.get(n, 0.0) + ramp}
        v = as_polar.get(n, 0.0)
        if v is not None:
            s["fgw_aligned_polar_series"] = second + v + ramp
        scores[tasks[n].key] = s
    return scores, tasks


def test_the_constants_are_the_rows(n41):
    assert n41.ARMS == ("fgw_polar", "fgw_aligned_polar_series") and n41.ENTRY == "fgw"
    assert n41.SECOND == "fgw_aligned_series" and n41.BACKTESTS == ("e75", "e725") and n41.ADOPTION_ARM == "fgw_polar"
    assert n41.BASE_OF == {"fgw_polar": "fgw", "fgw_aligned_polar_series": "fgw_aligned_series"}


def test_arms_equal_to_their_bases_close_the_line(n41):
    scores, tasks = tables(n41, spread=0.1)
    res = n41.check(scores, tasks)
    assert res["closed"] and not res["adopted"]
    assert res["requested_first_upload"] is None and res["requested_second_upload"] is None
    assert res["adoption"]["tasks"]["e75"]["mean_gain"] == pytest.approx(0.0)
    assert res["own_minus_base"]["fgw_aligned_polar_series"]["x1"] == pytest.approx(0.0)


def test_a_gain_on_both_backtests_passes_the_adoption_rule(n41):
    scores, tasks = tables(n41, polar={"e75": 0.5, "e725": 0.6})
    res = n41.check(scores, tasks)
    a = res["adoption"]
    assert a["keep"] and a["tasks"]["e75"]["seeds_ahead"] == 10
    assert a["tasks"]["e725"]["mean_gain"] == pytest.approx(0.6)
    assert res["adopted"] and not res["closed"]


def test_a_gain_on_one_backtest_and_a_loss_on_the_other_is_not_adopted(n41):
    scores, tasks = tables(n41, polar={"e75": 2.0, "e725": -0.5})
    res = n41.check(scores, tasks)
    assert not res["adoption"]["keep"] and res["adoption"]["tasks"]["e725"]["vetoes"] and not res["adopted"]


def test_the_first_upload_rule_reads_the_pair_beside_the_second_upload(n41):
    # ahead by 1 on every OTHER task, where the first upload leads its pair; nothing on the DEEP tasks
    scores, tasks = tables(n41, polar={n: 1.0 for n in n41.nt.OTHER})
    res = n41.check(scores, tasks)
    v = res["first_upload"]["first_uploads"]["fgw_polar"]
    assert v["O"] == pytest.approx(1.0) and v["D"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_first_upload"] == "fgw_polar" and not res["closed"]
    # a gain on the DEEP tasks only is hidden behind the second upload, which leads there by 3
    scores, tasks = tables(n41, polar={n: 1.0 for n in n41.nt.DEEP})
    res = n41.check(scores, tasks)
    assert res["first_upload"]["first_uploads"]["fgw_polar"]["D"] == pytest.approx(0.0)
    assert res["requested_first_upload"] is None


def test_the_second_upload_rule_reads_the_pair_beside_fgw(n41):
    # the polar second upload adds 1 where the second upload leads (DEEP); elsewhere fgw leads the pair either way
    scores, tasks = tables(n41, as_polar={n: 1.0 for n in n41.nt.DEEP})
    res = n41.check(scores, tasks)
    v = res["second_upload"]["second_uploads"]["fgw_aligned_polar_series"]
    assert v["D"] == pytest.approx(1.0) and v["O"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_second_upload"] == "fgw_aligned_polar_series" and not res["closed"]
    scores, tasks = tables(n41, as_polar={n: -1.0 for n in n41.nt.DEEP})
    res = n41.check(scores, tasks)
    assert res["requested_second_upload"] is None and res["closed"]


def test_an_arm_missing_on_a_task_is_not_readable_and_cannot_pass(n41):
    scores, tasks = tables(n41, as_polar={**{n: 5.0 for n in n41.nt.DEEP}, "x3": None})
    res = n41.check(scores, tasks)
    assert not res["readable"]["fgw_aligned_polar_series"] and "fgw_aligned_polar_series" not in res["own_minus_base"]
    assert not res["second_upload"]["second_uploads"]["fgw_aligned_polar_series"]["readable"]
    assert res["requested_second_upload"] is None


def test_the_reading_is_void_without_the_second_upload_and_the_text_states_the_outcome(n41, monkeypatch):
    scores, tasks = tables(n41)
    monkeypatch.setattr(n41.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n41.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n41.reading([])
    assert not res["void"] and res["closed"]
    txt = n41.text(res)
    assert "closed: no rule is passed" in txt and "adoption rule, fgw_polar against fgw" in txt
    assert "as a second upload beside fgw" in txt
    scores2, _ = tables(n41, polar={"e75": 0.5, "e725": 0.6})
    monkeypatch.setattr(n41.nt, "_load", lambda dirs: (scores2, [], ["synthetic:0000000"]))
    assert "the adoption rule is passed by fgw_polar" in n41.text(n41.reading([]))
    del scores["T2e/half"][n41.SECOND]
    monkeypatch.setattr(n41.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n41.reading([])
    assert res["void"] and any("no scores for fgw_aligned_series" in w for w in res["why"])
    assert n41.text(res).startswith("VOID")
