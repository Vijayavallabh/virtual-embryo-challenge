"""Readings of tries E39 (experiments/harness/n43_embryo_d_readings.py): three existing rules applied to three new
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
def n43():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n43_embryo_d_readings", HARNESS / "n43_embryo_d_readings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tables(n43, gains: dict[str, dict[str, float | None]] | None = None, spread: float = 0.0):
    """Per-seed scores on the nine tasks: fgw at 60; the second upload as built at 63 on the DEEP tasks and 58 on the
    others; fgw_series at 61 and 57; each new arm at its base + gains[arm][task] (default 0; None = no scores).
    spread adds a per-seed ramp to every arm alike, which leaves paired gains as set."""
    nt = n43.nt
    names = [*nt.DEEP, *nt.OTHER]
    tasks = {n: SimpleNamespace(key=f"T2e/{n}") for n in names}
    gains = gains or {}
    ramp = spread * np.arange(10)
    scores = {}
    for n in names:
        deep = n in nt.DEEP
        s = {n43.ENTRY: 60.0 + ramp, n43.SECOND: (63.0 if deep else 58.0) + ramp,
             "fgw_series": (61.0 if deep else 57.0) + ramp}
        for arm in n43.ARMS:
            v = gains.get(arm, {}).get(n, 0.0)
            if v is not None:
                s[arm] = s[n43.BASE_OF[arm]] + v
        scores[tasks[n].key] = s
    return scores, tasks


def test_the_constants_are_the_rows(n43):
    assert n43.ARMS == ("fgw_tissue", "fgw_tissue_series", "fgw_aligned_series_xt")
    assert n43.ADOPTION_ARMS == ("fgw_tissue",)
    assert n43.BASE_OF == {"fgw_tissue": "fgw", "fgw_tissue_series": "fgw_series",
                           "fgw_aligned_series_xt": "fgw_aligned_series"}
    assert n43.ENTRY == "fgw" and n43.SECOND == "fgw_aligned_series" and n43.BACKTESTS == ("e75", "e725")


def test_arms_equal_to_their_bases_close_the_line(n43):
    scores, tasks = tables(n43, spread=0.1)
    res = n43.check(scores, tasks)
    assert res["closed"] and res["adopted"] == []
    assert res["requested_first_upload"] is None and res["requested_second_upload"] is None
    assert res["adoption"]["fgw_tissue"]["tasks"]["e725"]["mean_gain"] == pytest.approx(0.0)
    assert res["own_minus_base"]["fgw_aligned_series_xt"]["x1"] == pytest.approx(0.0)
    assert res["seeds_equal_to_base"]["fgw_aligned_series_xt"]["half"] == 10


def test_a_gain_on_both_backtests_passes_the_adoption_rule(n43):
    scores, tasks = tables(n43, {"fgw_tissue": {"e75": 0.5, "e725": 0.6}})
    res = n43.check(scores, tasks)
    a = res["adoption"]["fgw_tissue"]
    assert a["keep"] and a["tasks"]["e75"]["seeds_ahead"] == 10
    assert a["tasks"]["e725"]["mean_gain"] == pytest.approx(0.6)
    assert res["adopted"] == ["fgw_tissue"] and not res["closed"]


def test_a_gain_on_one_backtest_and_a_loss_on_the_other_is_not_adopted(n43):
    scores, tasks = tables(n43, {"fgw_tissue": {"e75": 2.0, "e725": -0.5}})
    res = n43.check(scores, tasks)
    a = res["adoption"]["fgw_tissue"]
    assert not a["keep"] and a["tasks"]["e725"]["vetoes"] and res["adopted"] == []


def test_the_first_upload_rule_reads_the_pair_beside_the_second_upload(n43):
    scores, tasks = tables(n43, {"fgw_tissue": {n: 1.0 for n in n43.nt.OTHER}})
    res = n43.check(scores, tasks)
    v = res["first_upload"]["first_uploads"]["fgw_tissue"]
    assert v["O"] == pytest.approx(1.0) and v["D"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_first_upload"] == "fgw_tissue" and not res["closed"]
    scores, tasks = tables(n43, {"fgw_tissue": {n: 1.0 for n in n43.nt.DEEP}})  # hidden behind the second upload
    res = n43.check(scores, tasks)
    assert res["requested_first_upload"] is None


def test_the_second_upload_rule_reads_the_pair_beside_fgw_against_the_second_upload_as_built(n43):
    # the exact-total render adds 1 to the second upload where it leads (DEEP); elsewhere fgw leads the pair either way
    scores, tasks = tables(n43, {"fgw_aligned_series_xt": {n: 1.0 for n in n43.nt.DEEP}})
    res = n43.check(scores, tasks)
    v = res["second_upload"]["second_uploads"]["fgw_aligned_series_xt"]
    assert v["D"] == pytest.approx(1.0) and v["O"] == pytest.approx(0.0) and v["passes"]
    assert res["requested_second_upload"] == "fgw_aligned_series_xt" and not res["closed"]
    # fgw_series trails the second upload by 2 on the DEEP tasks: +1 there is still under it
    scores, tasks = tables(n43, {"fgw_tissue_series": {n: 1.0 for n in n43.nt.DEEP}})
    res = n43.check(scores, tasks)
    assert res["requested_second_upload"] is None and res["closed"]


def test_an_arm_missing_on_a_task_is_not_readable_and_cannot_pass(n43):
    gains = {"fgw_tissue_series": {**{n: 5.0 for n in n43.nt.DEEP}, "x3": None},
             "fgw_tissue": {"e75": None, "e725": 5.0}}
    scores, tasks = tables(n43, gains)
    res = n43.check(scores, tasks)
    assert not res["readable"]["fgw_tissue_series"] and "fgw_tissue_series" not in res["own_minus_base"]
    assert not res["second_upload"]["second_uploads"]["fgw_tissue_series"]["readable"]
    assert res["adoption"]["fgw_tissue"] is None and res["adopted"] == [] and res["closed"]


def test_the_count_of_seeds_equal_to_the_base_is_reported(n43):
    scores, tasks = tables(n43)
    key = tasks["e725"].key
    arm = scores[key]["fgw_aligned_series_xt"].copy()
    arm[[1, 4, 7]] += 0.8
    arm[2] += 0.004  # a difference inside the tolerance between runs still counts as the base
    scores[key]["fgw_aligned_series_xt"] = arm
    res = n43.check(scores, tasks)
    assert res["seeds_equal_to_base"]["fgw_aligned_series_xt"]["e725"] == 7
    txt = n43.text({"void": False, "runs": ["x:0"], **res})
    assert "seeds (of 10) in which fgw_aligned_series_xt scores as fgw_aligned_series:" in txt


def test_the_reading_is_void_without_a_base_and_the_text_states_the_outcome(n43, monkeypatch):
    scores, tasks = tables(n43)
    monkeypatch.setattr(n43.nt, "all_tasks", lambda: tasks)
    monkeypatch.setattr(n43.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n43.reading([])
    assert not res["void"] and res["closed"]
    txt = n43.text(res)
    assert "closed: no rule is passed" in txt and "adoption rule, fgw_tissue against fgw" in txt
    assert "fgw_tissue_series minus fgw_series" in txt and "Readings of tries E39a to E39c" in txt
    scores2, _ = tables(n43, {"fgw_tissue": {"e75": 0.5, "e725": 0.6}})
    monkeypatch.setattr(n43.nt, "_load", lambda dirs: (scores2, [], ["synthetic:0000000"]))
    assert "the adoption rule is passed by fgw_tissue" in n43.text(n43.reading([]))
    del scores["T2e/half"]["fgw_series"]
    monkeypatch.setattr(n43.nt, "_load", lambda dirs: (scores, [], ["synthetic:0000000"]))
    res = n43.reading([])
    assert res["void"] and any("no scores for fgw_series" in w for w in res["why"])
    assert n43.text(res).startswith("VOID")
