"""experiments/harness/n36_trend_second_upload.py (tries E35a and E35b) on stage names and on synthetic per-seed scores
written here: the sample rule, the extended tasks, the pair rule as fixed in its row, the void cases and the
not-buildable path. No committed result, data file or scorer is read, nothing is built from data, and nothing is
scored."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"


@pytest.fixture(scope="module")
def n36():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n36_trend_second_upload", HARNESS / "n36_trend_second_upload.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SEEDS = list(range(10))
NAMES = ("half", "late", "new1", "x1", "x2", "x3", "new2", "e75", "e725")
FIVE = {"E6.75", "E7.25", "E7.5", "E8.0", "E8.75"}
EXPECTED = {
    "half": {"E6.75", "E7.25", "E7.5", "E8.5"},
    "late": {"E6.75", "E7.25", "E7.5", "E8.25_late"},
    "new1": {"E6.75", "E7.25", "E7.5", "E8.5"},
    "new2": {"E6.75", "E7.25", "E7.5", "E8.75"},
    "x1": FIVE,
    "x2": FIVE,
    "x3": {"E6.75", "E7.25", "E7.5", "E8.0", "E8.25_late", "E8.75"},
    "e75": {"E6.75", "E7.25", "E8.0"},
    "e725": {"E6.75", "E7.5", "E8.0"},
}
KEYS = {"half": "T2e/E7.5+E8.5->E8.0", "late": "T2e/E7.5+E8.25_late->E8.0", "new1": "T2e/E7.5+E8.5->E8.25_late",
        "new2": "T2e/E7.5+E8.75->E8.5", "x1": "T2e/E8.0+E8.75->E8.25_late", "x2": "T2e/E8.0+E8.75->E8.5",
        "x3": "T2e/E8.25_late+E8.75->E8.5", "e75": "T2e/E7.25+E8.0->E7.5", "e725": "T2e/E6.75+E7.5->E7.25"}


def stems(files) -> set[str]:
    return {Path(f).name.removesuffix(".h5ad") for f in files}


# ------------------------------------------------------------------------------------------------ the sample rule
def test_the_sample_rule_gives_the_nine_expected_sets(n36):
    tasks = n36.all_tasks()
    assert set(tasks) == set(NAMES) and {n: t.key for n, t in tasks.items()} == KEYS
    for n, t in tasks.items():
        assert stems(t.inputs.values()) == EXPECTED[n], n
        assert stems(n36.trend_sample_files(t.inputs["before"], t.inputs["after"], t.target)) == EXPECTED[n]


def test_the_sample_rule_on_the_test_like_pair_and_its_refusals(n36):
    days = {"T2/E7.75.h5ad": 7.75}
    got = n36.trend_sample_files("T2/E7.5.h5ad", "T2/E8.0.h5ad", "T2/E7.75.h5ad", days=days)
    assert got == ["T2/E6.75.h5ad", "T2/E7.25.h5ad", "T2/E7.5.h5ad", "T2/E8.0.h5ad"]
    with pytest.raises(ValueError, match="no stage day"):  # E7.75 is not released: it has no day of its own
        n36.trend_sample_files("T2/E7.5.h5ad", "T2/E8.0.h5ad", "T2/E7.75.h5ad")
    with pytest.raises(ValueError, match="not inside its bracket"):
        n36.trend_sample_files("T2/E7.5.h5ad", "T2/E8.0.h5ad", "T2/E8.5.h5ad")
    # a stage strictly inside the bracket interval is never a sample
    assert "T2/E7.5.h5ad" not in n36.trend_sample_files("T2/E7.25.h5ad", "T2/E8.5.h5ad", "T2/E8.0.h5ad")
    assert n36.SERIES == ("T2/E6.75.h5ad", "T2/E7.25.h5ad", "T2/E7.5.h5ad", "T2/E8.0.h5ad")


def test_extended_tasks_keep_every_field_but_inputs(n36):
    base, ext = n36.nt.all_tasks(), n36.all_tasks()
    for n in NAMES:
        b, e = base[n], ext[n]
        assert dataclasses.replace(e, inputs=b.inputs) == b
        assert {r: e.inputs[r] for r in b.inputs} == dict(b.inputs)  # every configured role keeps its stage
        added = set(e.inputs) - set(b.inputs)
        assert all(r.startswith(n36.harness.TREND_SAMPLE_PREFIX) for r in added)
        assert {e.inputs[r] for r in added}.isdisjoint(set(b.inputs.values()))
        assert e.inputs["before"] == b.inputs["before"] and e.inputs["after"] == b.inputs["after"]
        assert n36.harness.TaskView.of(e).reference_role == "before"
        assert n36.harness.start_role(e) == n36.harness.start_role(b)  # "earlier" keeps its meaning
    # the two counted backtests get exactly the inputs configs/pseudo_tasks.yaml gives them
    for n in ("e75", "e725"):
        assert ext[n] == base[n]


def test_both_arms_apply_on_every_extended_task_and_the_trend_on_no_bracket_only_one(n36):
    base, ext = n36.nt.all_tasks(), n36.all_tasks()
    assert n36.NEW == ("fgw_trend_render", "fgw_aligned_trend") and n36.ARMS == ("floor", *n36.NEW)
    assert "fgw" not in n36.ARMS and "fgw_aligned_series" not in n36.ARMS  # never rescored
    for n in NAMES:
        for x in n36.NEW:
            assert n36.harness.CANDIDATES[x].applies(ext[n])
            assert n36.harness.CANDIDATES[x].applies(base[n]) == (n in ("e75", "e725"))
    assert n36.BASE_FIRST == "fgw" and n36.BASE_SECOND == "fgw_aligned_series" and n36.E33_CHECK == ("e75", "e725")


def test_extend_refuses_a_configured_input_outside_the_rule(n36):
    t = n36.nt.all_tasks()["half"]
    bad = dataclasses.replace(t, inputs={**t.inputs, "later": "T2/E8.75.h5ad"})
    with pytest.raises(ValueError, match="outside the sample rule"):
        n36.extend(bad)


# ------------------------------------------------------------------------------------------------ synthetic runs
def _write(d: Path, key: str, arms: dict[str, list[float]], inputs: dict, clean=True, commit="a" * 40,
           not_buildable=None, seeds=SEEDS) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    rows = [{"seed": s, "scores": {a: v[s] for a, v in arms.items()}} for s in seeds]
    r = {"task": {"key": key, "inputs": inputs}, "rows": rows, "provenance": {"tree_clean": clean, "commit": commit}}
    if not_buildable is not None:
        r["not_buildable"] = not_buildable
    p = d / f"{n36_safe(key)}.json"
    p.write_text(json.dumps(r))
    return p


def n36_safe(key: str) -> str:
    return "T2e_" + key.split("/", 1)[1].replace("+", "_").replace("->", "_")


def _setup(tmp_path, monkeypatch, n36, new: dict[str, dict[str, float]], base_second=None, extra_in_new=(),
           skip=(), commits=None, floor_new=50.0):
    """Stored runs hold floor, fgw and fgw_aligned_series on every task; the new run holds floor and the arms in `new`.
    new[arm][task name] is the arm's per-seed constant; fgw scores 50 and fgw_aligned_series base_second[task]."""
    tasks = n36.all_tasks()
    stored, runs = tmp_path / "stored", tmp_path / "runs"
    base_second = base_second or dict.fromkeys(NAMES, 52.0)
    for n in NAMES:
        t = tasks[n]
        _write(stored, t.key, {"floor": [50.0] * 10, "fgw": [50.0] * 10,
                               "fgw_aligned_series": [base_second[n]] * 10}, dict(n36.nt.all_tasks()[n].inputs))
        if n in skip:
            continue
        arms = {"floor": [floor_new] * 10, **{a: [v[n]] * 10 for a, v in new.items() if n in v}}
        for a in extra_in_new:
            arms[a] = [50.0] * 10
        _write(runs, t.key, arms, dict(t.inputs), commit=(commits or {}).get(n, "a" * 40))
    monkeypatch.setattr(n36.nt, "STORED", (stored,))
    monkeypatch.setattr(n36, "E28_RUNS", tmp_path / "empty")
    (tmp_path / "empty").mkdir(exist_ok=True)
    return runs


def both(a: dict, b: dict) -> dict:
    return {"fgw_trend_render": a, "fgw_aligned_trend": b}


def test_an_arm_that_beats_the_pair_everywhere_is_put_forward(tmp_path, monkeypatch, n36):
    res = n36.reading(_setup(tmp_path, monkeypatch, n36, both(dict.fromkeys(NAMES, 53.0), dict.fromkeys(NAMES, 51.0))))
    assert not res["void"] and res["put_forward"] == "fgw_trend_render"
    v = res["second_uploads"]["fgw_trend_render"]
    assert v["readable"] and v["passes"] and v["S0.45"] == pytest.approx(1.0) and v["worst9_abs"] == 53.0
    assert not res["second_uploads"]["fgw_aligned_trend"]["passes"]  # -1 everywhere
    d = res["described"]
    assert d["per_task"]["half"]["fgw_trend_render"] == {"mean": 53.0, "fgw_mean": 50.0, "paired_diff_from_fgw": 3.0,
                                                         "ahead_of_fgw": 10}
    assert d["fgw_trend_render_per_seed_for_the_E33_check"] == {"e75": [53.0] * 10, "e725": [53.0] * 10}


def test_when_both_pass_the_larger_S045_is_put_forward(tmp_path, monkeypatch, n36):
    res = n36.reading(_setup(tmp_path, monkeypatch, n36, both(dict.fromkeys(NAMES, 53.0), dict.fromkeys(NAMES, 54.0))))
    assert res["second_uploads"]["fgw_trend_render"]["passes"] and res["second_uploads"]["fgw_aligned_trend"]["passes"]
    assert res["put_forward"] == "fgw_aligned_trend"


def test_each_clause_of_the_rule_can_refuse(tmp_path, monkeypatch, n36):
    flat = dict.fromkeys(NAMES, 52.3)  # (a) a gain at the bar is not above it
    one = {**dict.fromkeys(NAMES, 52.0), "e75": 70.0}  # (d) the gain sits on one OTHER task only
    res = n36.reading(_setup(tmp_path / "1", monkeypatch, n36, both(flat, one)))
    a, e = res["second_uploads"]["fgw_trend_render"], res["second_uploads"]["fgw_aligned_trend"]
    assert a["S0.45"] == pytest.approx(0.3) and not a["rule"]["a_S0.45_above_bar"] and not a["passes"]
    assert e["rule"]["a_S0.45_above_bar"] and e["rule"]["b_S0.3_and_S0.6_not_below_0"]
    assert e["rule"]["c_worst9_not_lower"] and not e["rule"]["d_leave_one_task_out_positive"] and not e["passes"]
    assert e["S0.45_leave_one_task_out"]["e75"] == 0.0 and res["put_forward"] is None
    # (c) one task's pair mean drops more than 0.01 below the pair as built's worst
    worst = {**dict.fromkeys(NAMES, 56.0), "e725": 50.0}
    base = {**dict.fromkeys(NAMES, 52.0), "e725": 51.0}
    res = n36.reading(_setup(tmp_path / "2", monkeypatch, n36, both(worst, flat), base_second=base))
    w = res["second_uploads"]["fgw_trend_render"]
    assert w["rule"]["a_S0.45_above_bar"] and w["rule"]["d_leave_one_task_out_positive"]
    assert not w["rule"]["c_worst9_not_lower"] and not w["passes"]
    # (b) a DEEP loss turns S(0.6) negative while S(0.45) clears the bar
    split = {**{n: 48.0 for n in n36.nt.DEEP}, **{n: 54.5 for n in n36.nt.OTHER}}
    res = n36.reading(_setup(tmp_path / "3", monkeypatch, n36, both(split, flat)))
    b = res["second_uploads"]["fgw_trend_render"]
    assert b["S0.45"] == pytest.approx(0.475) and b["S0.6"] == pytest.approx(-0.2)
    assert b["rule"]["a_S0.45_above_bar"] and not b["rule"]["b_S0.3_and_S0.6_not_below_0"] and not b["passes"]
    assert res["put_forward"] is None


def test_void_conditions(tmp_path, monkeypatch, n36):
    good = both(dict.fromkeys(NAMES, 53.0), dict.fromkeys(NAMES, 53.0))

    def why(sub, **kw):
        res = n36.reading(_setup(tmp_path / sub, monkeypatch, n36, good, **kw) if "runs" not in kw else kw["runs"])
        assert res["void"]
        return " | ".join(res["why"])

    assert "rescored ['fgw']" in why("rescored", extra_in_new=("fgw",))
    assert "e725: no task file" in why("task", skip=("e725",))
    missing_arm = {"fgw_trend_render": dict.fromkeys(NAMES, 53.0),
                   "fgw_aligned_trend": {n: 53.0 for n in NAMES if n != "half"}}
    res = n36.reading(_setup(tmp_path / "arm", monkeypatch, n36, missing_arm))
    assert res["void"] and any("half: no scores for fgw_aligned_trend" in w for w in res["why"])
    assert "more than one commit" in why("commits", commits={"x3": "b" * 40})
    assert "differs between runs" in why("floor", floor_new=50.02)
    # a dirty tree, a missing seed, inputs that are not the rule's, a run without the floor
    runs = _setup(tmp_path / "dirty", monkeypatch, n36, good)
    t = n36.all_tasks()["new2"]
    _write(runs, t.key, {"floor": [50.0] * 10, **{a: [53.0] * 10 for a in n36.NEW}}, dict(t.inputs), clean=False)
    assert "a dirty tree" in why("dirty", runs=runs)
    runs = _setup(tmp_path / "seed", monkeypatch, n36, good)
    _write(runs, t.key, {"floor": [50.0] * 10, **{a: [53.0] * 10 for a in n36.NEW}}, dict(t.inputs), seeds=SEEDS[:9])
    assert "seeds are not" in why("seed", runs=runs)
    runs = _setup(tmp_path / "inputs", monkeypatch, n36, good)
    _write(runs, t.key, {"floor": [50.0] * 10, **{a: [53.0] * 10 for a in n36.NEW}},
           dict(n36.nt.all_tasks()["new2"].inputs))
    assert "new2: the run's inputs are not the sample rule's" in why("inputs", runs=runs)
    runs = _setup(tmp_path / "nofloor", monkeypatch, n36, good)
    _write(runs, t.key, {a: [53.0] * 10 for a in n36.NEW}, dict(t.inputs))
    assert "new2: the new run has no floor" in why("nofloor", runs=runs)


def test_the_stored_runs_must_not_hold_the_new_arms(tmp_path, monkeypatch, n36):
    runs = _setup(tmp_path, monkeypatch, n36, both(dict.fromkeys(NAMES, 53.0), dict.fromkeys(NAMES, 53.0)))
    t = n36.nt.all_tasks()["x1"]
    _write(tmp_path / "e28", t.key, {"floor": [50.0] * 10, "fgw_trend_render": [53.0] * 10}, dict(t.inputs))
    monkeypatch.setattr(n36, "E28_RUNS", tmp_path / "e28")
    res = n36.reading(runs)
    assert res["void"] and any("already holds ['fgw_trend_render']" in w for w in res["why"])


# ------------------------------------------------------------------------------------------------ not buildable
REFUSAL = "no extraembryonic label held by every sample: the trend factor is not defined"


def _fake_runner(score: dict[str, float]):
    """run_task's shape without data or scorer: each candidate is built (its build may raise) and given a constant."""
    calls = []

    def runner(t, cands, seeds, data_dir):
        calls.append([c.name for c in cands])
        rows = []
        for s in seeds:
            row = {"seed": s, "scores": {"floor": 50.0}, "skills": {"floor": {}}, "metrics": {"floor": {}},
                   "n_cells": {"floor": 1}, "groups": {"floor": {}}, "reads": {}}
            for c in cands:
                if c.name == "floor":
                    continue
                c.build(dataclasses.make_dataclass("Ctx", ["seed"])(seed=1_000_000 + s))
                row["scores"][c.name] = score[c.name]
                for k in ("skills", "metrics", "groups"):
                    row[k][c.name] = {}
                row["n_cells"][c.name], row["reads"][c.name] = 1, []
            rows.append(row)
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "candidates": [c.name for c in cands], "params": {c.name: {} for c in cands},
                "provenance": {"commit": "a" * 40, "tree_clean": True}}

    return runner, calls


def _patch_builds(monkeypatch, n36, fail_on: set[str], seed_fails=0):
    cands = n36.harness.CANDIDATES

    def ok(ctx):
        return "built"

    def bad(ctx):
        if ctx.seed - n36.harness.BUILD_SEED_OFFSET >= seed_fails:
            raise ValueError(REFUSAL + " (3 samples, 0.25 of cells)")
        return "built"

    for x in n36.NEW:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=bad if x in fail_on else ok))


def test_a_failing_build_stops_neither_the_other_arm_nor_the_task(tmp_path, monkeypatch, n36):
    _patch_builds(monkeypatch, n36, {"fgw_aligned_trend"}, seed_fails=3)
    runner, calls = _fake_runner({"fgw_trend_render": 53.0, "fgw_aligned_trend": 60.0})
    t = n36.all_tasks()["x2"]
    res = n36.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n36.NEW], ["floor", "fgw_trend_render"], ["floor", "fgw_aligned_trend"]]
    rec = res["not_buildable"]["fgw_aligned_trend"]
    assert rec["status"] == "not buildable on this task" and rec["seed"] == 3
    assert rec["message"] == "ValueError: " + REFUSAL + " (# samples, # of cells)"  # no data value
    assert all(set(r["scores"]) == {"floor", "fgw_trend_render"} for r in res["rows"])
    assert [r["seed"] for r in res["rows"]] == SEEDS and res["candidates"] == ["floor", "fgw_trend_render"]
    # both arms failing still leaves the floor scored
    _patch_builds(monkeypatch, n36, set(n36.NEW))
    runner, calls = _fake_runner({})
    res = n36.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n36.NEW) and calls[-1] == ["floor"]
    assert all(set(r["scores"]) == {"floor"} for r in res["rows"])
    # nothing fails: one run, nothing recorded
    _patch_builds(monkeypatch, n36, set())
    runner, calls = _fake_runner({"fgw_trend_render": 53.0, "fgw_aligned_trend": 54.0})
    res = n36.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n36.NEW]] and res["not_buildable"] == {}


def test_a_failure_outside_a_build_stops_the_task(tmp_path, monkeypatch, n36):
    _patch_builds(monkeypatch, n36, set())

    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError):
        n36.score_task(n36.all_tasks()["x2"], tmp_path, runner=runner)


def test_a_not_buildable_arm_is_not_readable_and_the_other_arm_is_still_read(tmp_path, monkeypatch, n36):
    runs = _setup(tmp_path, monkeypatch, n36, both(dict.fromkeys(NAMES, 53.0), dict.fromkeys(NAMES, 60.0)))
    _patch_builds(monkeypatch, n36, {"fgw_aligned_trend"})
    runner, _ = _fake_runner({"fgw_trend_render": 53.0, "fgw_aligned_trend": 60.0})
    t = n36.all_tasks()["e725"]
    res = n36.score_task(t, tmp_path, runner=runner)
    (runs / f"{n36_safe(t.key)}.json").write_text(json.dumps(res))
    out = n36.reading(runs)
    assert not out["void"]
    x = out["second_uploads"]["fgw_aligned_trend"]
    assert x["readable"] is False and x["status"] == "not readable" and x["passes"] is False
    assert set(x["not_buildable_on"]) == {"e725"} and x["not_buildable_on"]["e725"]["seed"] == 0
    assert out["second_uploads"]["fgw_trend_render"]["passes"] and out["put_forward"] == "fgw_trend_render"
    assert out["described"]["per_task"]["e725"]["fgw_aligned_trend"] is None
    assert "not readable: not buildable on this task on e725" in n36.text(out)


def test_read_refuses_to_overwrite_and_writes_both_files(tmp_path, monkeypatch, n36):
    runs = _setup(tmp_path, monkeypatch, n36, both(dict.fromkeys(NAMES, 53.0), dict.fromkeys(NAMES, 51.0)))
    out = tmp_path / "out"
    n36.main(["read", "--runs", str(runs), "--out", str(out)])
    assert {p.name for p in out.iterdir()} == {"n36.json", "n36.txt"}
    assert json.loads((out / "n36.json").read_text())["put_forward"] == "fgw_trend_render"
    with pytest.raises(SystemExit, match="not empty"):
        n36.main(["read", "--runs", str(runs), "--out", str(out)])
