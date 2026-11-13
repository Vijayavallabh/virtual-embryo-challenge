"""experiments/harness/n35_second_upload_arms.py on synthetic per-seed scores written here: the rule as fixed in its
row, the refusal to rescore fgw or fgw_aligned_series, and the void cases. No committed result, data file or scorer is
read, and nothing is scored."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"


@pytest.fixture(scope="module")
def n35():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n35_second_upload_arms", HARNESS / "n35_second_upload_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SEEDS = list(range(10))
NAMES = ("half", "late", "new1", "x1", "x2", "x3", "new2", "e75", "e725")


def _write(d: Path, key: str, arms: dict[str, list[float]], clean=True, commit="a" * 40) -> None:
    d.mkdir(parents=True, exist_ok=True)
    rows = [{"seed": s, "scores": {a: v[s] for a, v in arms.items()}} for s in SEEDS]
    safe = "T2e_" + key.split("/", 1)[1].replace("+", "_").replace("->", "_")
    (d / f"{safe}.json").write_text(json.dumps({"task": {"key": key}, "rows": rows,
                                                "provenance": {"tree_clean": clean, "commit": commit}}))


def _setup(tmp_path, monkeypatch, n35, new: dict[str, dict[str, float]], base_second=None, extra_in_new=()):
    """Stored runs hold floor, fgw and fgw_aligned_series on every task; the new run holds floor and the two arms.
    new[arm][task name] is the arm's per-seed constant; fgw scores 50 and fgw_aligned_series base_second[task]."""
    tasks = n35.nt.all_tasks()
    stored, runs = tmp_path / "stored", tmp_path / "runs"
    base_second = base_second or dict.fromkeys(NAMES, 52.0)
    for n in NAMES:
        key = tasks[n].key
        _write(stored, key, {"floor": [50.0] * 10, "fgw": [50.0] * 10, "fgw_aligned_series": [base_second[n]] * 10})
        arms = {"floor": [50.0] * 10, **{a: [new[a][n]] * 10 for a in n35.NEW}}
        for a in extra_in_new:
            arms[a] = [50.0] * 10
        _write(runs, key, arms)
    monkeypatch.setattr(n35.nt, "STORED", (stored,))
    monkeypatch.setattr(n35, "E28_RUNS", tmp_path / "empty")
    (tmp_path / "empty").mkdir()
    return runs


def test_the_nine_tasks_and_arms_are_the_briefs(n35):
    assert set(n35.nt.all_tasks()) == set(NAMES) and set(n35.nt.DEEP) | set(n35.nt.OTHER) == set(NAMES)
    assert n35.NEW == ("fgw_aligned_hi", "fgw_early_specialist") and n35.BASE_SECOND == "fgw_aligned_series"
    assert n35.BASE_FIRST == "fgw" and n35.nt.PS == (0.3, 0.45, 0.6) and n35.nt.BAR == 0.3 and n35.nt.TOL == 0.01
    assert "fgw" not in n35.ARMS and "fgw_aligned_series" not in n35.ARMS  # never rescored
    for name in n35.NEW:
        assert all(n35.harness.CANDIDATES[name].applies(t) for t in n35.nt.all_tasks().values())


def test_an_arm_that_beats_the_pair_everywhere_is_put_forward(tmp_path, monkeypatch, n35):
    new = {"fgw_aligned_hi": dict.fromkeys(NAMES, 53.0), "fgw_early_specialist": dict.fromkeys(NAMES, 51.0)}
    res = n35.reading(_setup(tmp_path, monkeypatch, n35, new))
    assert not res["void"] and res["put_forward"] == "fgw_aligned_hi"
    v = res["second_uploads"]["fgw_aligned_hi"]
    assert v["S0.45"] == pytest.approx(1.0) and v["passes"] and v["worst9_abs"] == 53.0
    assert not res["second_uploads"]["fgw_early_specialist"]["passes"]  # -1 everywhere


def test_each_clause_of_the_rule_can_refuse(tmp_path, monkeypatch, n35):
    names = list(NAMES)
    # (a) a gain at the bar is not above it
    flat = dict.fromkeys(NAMES, 52.3)
    # (c) one task's pair mean drops more than 0.01 below the pair as built's worst
    worst = {**dict.fromkeys(NAMES, 56.0), "e725": 50.0}
    base = {**dict.fromkeys(NAMES, 52.0), "e725": 51.0}
    # (d) the gain sits on one OTHER task only
    one = {**dict.fromkeys(NAMES, 52.0), "e75": 70.0}
    res = n35.reading(_setup(tmp_path / "1", monkeypatch, n35, {"fgw_aligned_hi": flat, "fgw_early_specialist": one}))
    a, e = res["second_uploads"]["fgw_aligned_hi"], res["second_uploads"]["fgw_early_specialist"]
    assert a["S0.45"] == pytest.approx(0.3) and not a["rule"]["a_S0.45_above_bar"] and not a["passes"]
    assert e["rule"]["a_S0.45_above_bar"] and e["rule"]["c_worst9_not_lower"]
    assert not e["rule"]["d_leave_one_task_out_positive"] and e["S0.45_leave_one_task_out"]["e75"] == 0.0
    res = n35.reading(_setup(tmp_path / "2", monkeypatch, n35, {"fgw_aligned_hi": worst, "fgw_early_specialist": flat},
                             base_second=base))
    w = res["second_uploads"]["fgw_aligned_hi"]
    assert w["rule"]["a_S0.45_above_bar"] and not w["rule"]["c_worst9_not_lower"] and not w["passes"]
    # (b) a loss on the DEEP side large enough to turn S(0.6) negative while S(0.45) clears the bar
    split = {**{n: 50.0 for n in n35.nt.DEEP}, **{n: 58.0 for n in n35.nt.OTHER}}
    res = n35.reading(_setup(tmp_path / "3", monkeypatch, n35, {"fgw_aligned_hi": split, "fgw_early_specialist": flat}))
    b = res["second_uploads"]["fgw_aligned_hi"]
    assert b["D"] == pytest.approx(-2.0) and b["O"] == pytest.approx(6.0)
    assert b["S0.45"] > 0.3 and b["S0.6"] >= 0 and b["S0.3"] > 0  # here (b) holds
    split2 = {**{n: 48.0 for n in n35.nt.DEEP}, **{n: 54.5 for n in n35.nt.OTHER}}  # D = -2 (fgw's 50 floors it)
    res = n35.reading(_setup(tmp_path / "4", monkeypatch, n35, {"fgw_aligned_hi": split2,
                                                                 "fgw_early_specialist": flat}))
    b2 = res["second_uploads"]["fgw_aligned_hi"]
    assert b2["S0.45"] == pytest.approx(0.475) and b2["S0.6"] == pytest.approx(-0.2)
    assert b2["rule"]["a_S0.45_above_bar"] and not b2["rule"]["b_S0.3_and_S0.6_not_below_0"] and not b2["passes"]
    assert res["put_forward"] is None and len(names) == 9


def test_rescored_bases_or_missing_arms_void_the_reading(tmp_path, monkeypatch, n35):
    good = dict.fromkeys(NAMES, 53.0)
    runs = _setup(tmp_path / "r", monkeypatch, n35, {"fgw_aligned_hi": good, "fgw_early_specialist": good},
                  extra_in_new=("fgw",))
    res = n35.reading(runs)
    assert res["void"] and any("rescored ['fgw']" in w for w in res["why"])
    runs = _setup(tmp_path / "m", monkeypatch, n35, {"fgw_aligned_hi": good, "fgw_early_specialist": good})
    next(iter(sorted(runs.glob("*.json")))).unlink()
    res = n35.reading(runs)
    assert res["void"] and any("no scores for fgw_aligned_hi" in w for w in res["why"])


def test_read_refuses_to_overwrite_and_writes_both_files(tmp_path, monkeypatch, n35):
    good = dict.fromkeys(NAMES, 53.0)
    runs = _setup(tmp_path, monkeypatch, n35, {"fgw_aligned_hi": good, "fgw_early_specialist": good})
    out = tmp_path / "out"
    n35.main(["read", "--runs", str(runs), "--out", str(out)])
    assert {p.name for p in out.iterdir()} == {"n35.json", "n35.txt"}
    assert json.loads((out / "n35.json").read_text())["put_forward"] == "fgw_aligned_hi"
    with pytest.raises(SystemExit, match="not empty"):
        n35.main(["read", "--runs", str(runs), "--out", str(out)])
    assert np.isfinite(json.loads((out / "n35.json").read_text())["base_pair_worst9"])
