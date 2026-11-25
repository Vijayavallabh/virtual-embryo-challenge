"""Row E58, gate G1 of row R59's file T2embryo_s2x_ovt23.h5ad, as the cold challenge of the step-change round rewrites
the design's G1: the read script experiments/coords/e58_mix_lean_read.py, the row check experiments/coords/e58_row.py
and the launcher experiments/coords/run_e58.sh. Synthetic run and gate records, synthetic ledgers and throwaway git
repositories only: nothing is read from data/ and nothing is scored. Two tests read committed practice records by code
(E55's reading and E47's task files: their sha256 pins and the tasks they hold) and print nothing of them.

Pinned: the arm is O's call with t_mix = min(1, t + 1/6) and nothing else changed (2/3 at t = 1/2, the builder's
--mix-lean later); the counted merge is two truths at equal weight (half and late are one truth, E8.0; new1 is
E8.25_late); new2 and x1 to x3 are reported and decide nothing; e75 and e725 are never run; each clause (P1 at least
+0.3 with 8 of 10 seeds ahead, P2 no group behind by E49's definition, P3 no refusal) and the idle case; the void
conditions; the J7 guard against E47's floor and E55's committed O figure, at 0.01, and its dry run on committed
records; the row's form and the fixed outcome words that the build launcher reads; the launcher's syntax and its
refusals (wrong commit 96, an unclean tree 98, a missing or malformed row E58 95), traced with bash -x in throwaway
clones."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from conftest import git, git_env

from vec import harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
LAUNCHER_REL = "experiments/coords/run_e58.sh"
ROW_REL = "experiments/coords/e58_row.py"
SEEDS = list(range(10))
COMMIT = "a" * 40
O = "fgw_shape0x_ov"  # noqa: E741  (the specification names the control O)
ARM = "fgw_shape0x_ov_t23"
GROUPS = ("shape_growth", "expression_change", "cell_state", "local_organisation")
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def e58():
    return _load("e58_mix_lean_read", "experiments/coords/e58_mix_lean_read.py")


@pytest.fixture(scope="module")
def er():
    return _load("e58_row", ROW_REL)


def _rows(er, outcome="not yet run (placed before its run)", first=None, cells=11):
    first = er.FIRST_CELL_OPENING + " (synthetic); diagnostic (J5, noise floor 6.54))" if first is None else first
    mid = [f"c{k}" for k in range(2, cells)]
    return "| " + " | ".join([first, *mid[:cells - 2], outcome]) + " |"


def _ledger(tmp_path: Path, *rows: str, name="TRIES.md") -> Path:
    p = tmp_path / name
    p.write_text("| ID | what |\n|---|---|\n| E57 (another row) | embryo |\n" + "".join(r + "\n" for r in rows))
    return p


def _pass_cell(er):
    return f"Run at the run commit {COMMIT}; in the read's fixed words: \"{er.PASS_TEXT}\" One try, counted."


# ------------------------------------------------------------------------------------------------ the arm
def test_the_constants(e58, er):
    assert e58.TASKS == ("half", "late", "new1", "new2", "x1", "x2", "x3")
    assert e58.COUNTED_TASKS == ("half", "late", "new1") and e58.REPORTED_TASKS == ("new2", "x1", "x2", "x3")
    assert e58.BARRED_TASKS == ("e75", "e725") and not set(e58.BARRED_TASKS) & set(e58.TASKS)
    assert (e58.BAR, e58.AHEAD, e58.GATE_TOL, e58.O_TOL, e58.BUILD) == (0.3, 8, 0.01, 0.01, "main")
    assert e58.O == O and e58.ARM == ARM and e58.OTHER_GROUPS == GROUPS[1:]
    assert e58.E55_READING_SHA256 == "e500c9fa4161720191562598af76f8fb28bc12b0c3b0dfa9fdd3e38a38153ab5"
    assert e58.E55_RUN_COMMIT == "4535f98c98ae1a2a8488185e8a99e4a2e8d92d86"
    assert list(e58.SEEDS) == SEEDS and e58.er is not None and e58.er.PASS_TEXT == er.PASS_TEXT
    for k, v in e58.CLAUSES.items():
        assert v in e58.PASS_RULE, k
    assert "e75 and e725 are not run" in e58.PASS_RULE and "idle case" in e58.PASS_RULE
    assert ARM not in harness.CANDIDATES  # the arm lives in the read script, not in vec.harness's registry


def test_the_arms_t_mix_is_the_builders_at_the_test_pair(e58):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    assert e58.arm_t_mix(0.5) == 2 / 3 == build_test.EMBRYO_OT23_T_MIX  # the same double: 1/2 + 1/6 rounds to it
    assert e58.arm_t_mix(0.8) == pytest.approx(0.8 + 1 / 6) and e58.arm_t_mix(0.9) == 1.0
    tasks = e58.nt.all_tasks()
    assert e58.frac_of(tasks["half"]) == 0.5 and e58.frac_of(tasks["late"]) == pytest.approx(0.625)
    assert e58.frac_of(tasks["new1"]) == pytest.approx(0.8)


def test_the_arm_is_os_call_with_the_leaned_t_mix_and_nothing_else(e58, monkeypatch):
    calls = []
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: calls.append(kw) or kw)
    for name in ("half", "late", "new1", "x2"):
        t = e58.nt.all_tasks()[name]
        ctx = SimpleNamespace(task=harness.TaskView.of(t), seed=harness.BUILD_SEED_OFFSET + 1)
        calls.clear()
        harness.CANDIDATES[O].build(ctx)
        e58.ARM_CANDIDATE.build(ctx)
        o_kw, arm_kw = calls
        assert "t_mix" not in o_kw or o_kw["t_mix"] is None
        assert arm_kw.pop("t_mix") == e58.arm_t_mix(e58.frac_of(t))
        assert arm_kw == {k: v for k, v in o_kw.items() if k != "t_mix"}  # t_pos 0, floor_exact, overlay, as O
        assert e58.ARM_CANDIDATE.applies(t) and harness.CANDIDATES[O].applies(t)


def test_interpolate_clouds_takes_the_lean_and_keeps_os_coupling():
    rng = np.random.default_rng(4)
    Cb, Ca = rng.normal(size=(300, 3)) * [3, 1, 0.5], rng.normal(size=(260, 3)) * [3.5, 1, 0.4]
    Eb, Ea = rng.gamma(2.0, 1.0, size=(300, 6)), rng.gamma(2.0, 1.2, size=(260, 6))
    kw = dict(target_rms=5.0, t_pos=0.0, positions="floor_exact", shape_form="overlay", n_anchor=80)
    o = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 120, 1_000_003, **kw)
    a = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 120, 1_000_003, t_mix=2 / 3, **kw)
    assert int((o.stage == 0).sum()) == 60 and int((a.stage == 0).sum()) == 40 and len(a.rows) == 120
    assert a.record["rotation_deg"] == o.record["rotation_deg"] and a.record["t_mix"] == 2 / 3
    rms = float(np.sqrt(((a.coords - a.coords.mean(0)) ** 2).sum(1).mean()))
    assert rms == pytest.approx(5.0, rel=1e-3)


# ------------------------------------------------------------------------------------------------ the pass rule
def _good():
    return {"board": np.full(10, 0.5), "groups": {g: np.zeros(10) for g in GROUPS[1:]}, "arm_refused": 0}


def _seven_ahead(mean_target=0.5):
    D = np.full(10, 1.5)
    D[:3] = -0.5
    return D * (mean_target / D.mean())


def test_decide_passes_when_every_clause_holds_and_at_the_bar_itself(e58, er, monkeypatch):
    r = e58.decide(_good())
    assert r["passes"] and r["outcome"] == "PASS" and r["outcome_text"] == er.PASS_TEXT and not r["idle"]
    monkeypatch.setattr(e58, "BAR", 0.375)  # a bar a float mean can equal exactly (ten 0.3s average below 0.3)
    d = _good()
    d["board"] = np.full(10, 0.375)  # "at least": the bar itself passes and is not idle
    assert e58.decide(d)["passes"] and not e58.decide(d)["idle"]
    d["board"] = np.full(10, np.nextafter(0.375, 0.0))
    assert not e58.decide(d)["passes"] and e58.decide(d)["idle"]


@pytest.mark.parametrize("case, clause", [("P1_bar", "P1"), ("P1_seeds", "P1"), ("P2", "P2"), ("P2_cell", "P2"),
                                          ("P3", "P3")])
def test_each_clause_can_stop(e58, er, case, clause):
    d = _good()
    if case == "P1_bar":
        d["board"] = np.full(10, 0.2999)
    elif case == "P1_seeds":
        d["board"] = _seven_ahead()
    elif case == "P2":
        D = np.full(10, -0.0006)
        D[:2] = 0.0001
        d["groups"]["local_organisation"] = D  # 8 of 10 below 0, mean below 0: behind, no tolerance
    elif case == "P2_cell":
        d["groups"]["cell_state"] = np.full(10, -1e-6)
    else:
        d["arm_refused"] = 1
    r = e58.decide(d)
    assert not r["passes"] and not r["clauses"][clause] and r["outcome"] == "STOP"
    assert all(v for k, v in r["clauses"].items() if k != clause)
    assert r["outcome_text"] in (er.STOP_TEXT, er.IDLE_TEXT)


def test_the_idle_case_and_a_group_with_seven_seeds_behind(e58, er):
    d = _good()
    d["board"] = np.full(10, -0.29)
    r = e58.decide(d)
    assert r["idle"] and r["outcome_text"] == er.IDLE_TEXT and not r["passes"]
    d["board"] = np.full(10, -0.5)  # a loss is not idle: a plain STOP
    r = e58.decide(d)
    assert not r["idle"] and r["outcome_text"] == er.STOP_TEXT
    d = _good()
    D = np.full(10, -0.001)
    D[:3] = 0.0001
    d["groups"]["expression_change"] = D  # 7 below: not behind by E49's definition
    assert e58.decide(d)["passes"]
    d["groups"]["shape_growth"] = np.full(10, -1.0)  # a fourth key is not a group the clause reads
    assert not e58.decide(d)["clauses"]["P2"]


# ------------------------------------------------------------------------------------------------ the reading
def _row(seed, shift, grp_shift):
    names = ("floor", O, ARM)
    scores = {n: 50.0 if n == "floor" else 47.0 + 0.01 * seed + shift.get(n, 0.0) for n in names}
    skills = {n: dict.fromkeys(METRICS, 0.5) for n in names}
    groups = {n: {g: 0.4 + grp_shift.get((n, g), 0.0) for g in GROUPS} for n in names}
    return {"seed": seed, "scores": scores, "skills": skills, "groups": groups}


def _facts(e58, task, refused=None):
    frac = e58.frac_of(e58.nt.all_tasks()[task])
    return {"board_split": "pseudo", "frac": frac, "n": {O: 100, ARM: 100},
            "records": {O: {"t_mix": frac, "t": frac, "n_jittered": 0, "target_rms": 1.5},
                        ARM: {"t_mix": e58.arm_t_mix(frac), "t": frac, "n_jittered": 0, "target_rms": 1.5}},
            "refused": {ARM: refused}}


def _recs(e58, shift=None, grp_shift=None, task_shift=None, refused_on=None):
    out = {}
    for name in e58.TASKS:
        rs = []
        for s in SEEDS:
            sh = dict(shift or {})
            for n, v in ((task_shift or {}).get(name) or {}).items():
                sh[n] = sh.get(n, 0.0) + v
            res = {"provenance": {"commit": COMMIT, "tree_clean": True}, "build": "main", "protocol": "bake",
                   "board": {"split": "pseudo"}, "rows": [_row(s, sh, grp_shift or {})]}
            rs.append({"seed": s, "run_commit": COMMIT, "builds": {"main": res},
                       "facts": _facts(e58, name, (refused_on or {}).get((name, s)))})
        out[name] = rs
    return out


GOOD = {ARM: 0.8}


def test_a_pass_writes_the_fixed_pass_words_on_their_own_line(e58, er, tmp_path):
    res = e58.reading(_recs(e58, GOOD), [], COMMIT)
    assert not res["void"] and res["passes"] and res["outcome"] == "PASS" and res["outcome_text"] == er.PASS_TEXT
    assert res["merged"]["board"]["mean"] == pytest.approx(0.8) and res["price"]["losses"] == {}
    assert set(res["price"]["per_truth"]) == {"T2/E8.0.h5ad", "T2/E8.25_late.h5ad"}
    text = e58.text(res)
    assert er.PASS_TEXT in text.splitlines() and not any(w in text for w in er.NOT_PASS)
    assert "--mix-lean later --captain-ruling R59" in text and "it moves the overlay's form" in text.splitlines()[0]
    assert all(line in text for line in e58.PRICE_LINES) and "e75 and e725 were not run" in text
    js, txt = tmp_path / "reading.json", tmp_path / "reading.txt"
    js.write_text(json.dumps({**res, "reading_code": {"head": COMMIT, "clean": True}}, default=float))
    txt.write_text(text)
    assert er.reading_commit(js, txt) == COMMIT


def test_half_and_late_are_one_truth_and_the_reported_tasks_decide_nothing(e58):
    # +1.4 on half alone: E8.0's truth gains 0.7, E8.25_late 0: the two-truth merge is 0.35 (by task it would be 0.467)
    res = e58.reading(_recs(e58, task_shift={"half": {ARM: 1.4}}), [], COMMIT)
    assert res["merged"]["board"]["mean"] == pytest.approx(0.35) and res["passes"]
    assert res["price"]["per_truth"] == {"T2/E8.0.h5ad": pytest.approx(0.7), "T2/E8.25_late.h5ad": pytest.approx(0.0)}
    res = e58.reading(_recs(e58, task_shift={"half": {ARM: 1.0}}), [], COMMIT)  # merged 0.25: idle, by task 0.33
    assert res["merged"]["board"]["mean"] == pytest.approx(0.25) and res["idle"] and not res["passes"]
    loss = {t: {ARM: -5.0} for t in ("new2", "x1", "x2", "x3")}
    res = e58.reading(_recs(e58, GOOD, task_shift=loss), [], COMMIT)
    assert res["passes"] and res["reported_deciding_nothing"]["reported_tasks_merge"]["mean"] == pytest.approx(-4.2)
    gain = {t: {ARM: 5.0} for t in ("new2", "x1", "x2", "x3")}
    res = e58.reading(_recs(e58, task_shift=gain), [], COMMIT)
    assert not res["passes"] and res["idle"] and res["merged"]["board"]["mean"] == pytest.approx(0.0)
    assert res["tasks"]["x1"]["counted"] is False and res["tasks"]["new1"]["counted"] is True
    assert res["barred_not_run"] == ["e75", "e725"] and "e75" not in res["tasks"]


def test_shape_growth_is_reported_and_decides_nothing(e58):
    res = e58.reading(_recs(e58, GOOD, grp_shift={(ARM, "shape_growth"): -0.2}), [], COMMIT)
    assert res["passes"]
    assert res["reported_deciding_nothing"]["shape_growth_counted_merge_recused"]["mean"] == pytest.approx(-0.2)
    res = e58.reading(_recs(e58, GOOD, grp_shift={(ARM, "cell_state"): -0.2}), [], COMMIT)
    assert not res["passes"] and not res["clauses"]["P2"]


def test_a_refusal_of_the_arm_is_a_stop_not_a_void(e58, er):
    res = e58.reading(_recs(e58, GOOD, refused_on={("x2", 5): "cloud guard"}), [], COMMIT)
    assert not res["void"] and not res["passes"] and not res["clauses"]["P3"] and res["arm_refused_task_seeds"] == 1
    text = e58.text(res)
    assert "STOP, not void" in text and er.STOP_TEXT in text.splitlines() and er.PASS_TEXT not in text


@pytest.mark.parametrize("breach", ["seed", "commit", "run_commit", "unclean", "protocol", "build", "test_board",
                                    "test_board_build", "missing_arm", "count", "o_t_mix", "arm_t_mix", "size",
                                    "gate", "missing_task"])
def test_the_void_conditions(e58, er, breach, tmp_path):
    recs = _recs(e58, GOOD)
    why = []
    r = recs["late"][3]
    res = r["builds"]["main"]
    f = r["facts"]
    if breach == "seed":
        recs["late"].pop(3)
    elif breach == "commit":
        res["provenance"]["commit"] = "b" * 40
    elif breach == "run_commit":
        r["run_commit"] = "b" * 40
    elif breach == "unclean":
        res["provenance"]["tree_clean"] = False
    elif breach == "protocol":
        res["protocol"] = "live"
    elif breach == "build":
        r["builds"] = {"prefix": res}
    elif breach == "test_board":
        f["board_split"] = "test_interp"
    elif breach == "test_board_build":
        res["board"]["split"] = "test_interp"
    elif breach == "missing_arm":
        res["rows"][0]["scores"].pop(ARM)
    elif breach == "count":
        f["n"][ARM] = 99
    elif breach == "o_t_mix":
        f["records"][O]["t_mix"] = 0.5  # late's t is 0.625
    elif breach == "arm_t_mix":
        f["records"][ARM]["t_mix"] = 2 / 3  # the test pair's value, not late's 0.625 + 1/6
    elif breach == "size":
        f["records"][ARM]["target_rms"] = 1.6
    elif breach == "gate":
        why = ["the gate did not pass on late at the run commit"]
    else:
        recs.pop("x3")
    out = e58.reading(recs, why, COMMIT)
    assert out["void"] and out["why"] and out["outcome"] == "VOID" and not out["passes"]
    text = e58.text(out)
    assert text.startswith(f"Row E58 (run commit {COMMIT}): VOID") and er.VOID_TEXT in text
    js, txt = tmp_path / "reading.json", tmp_path / "reading.txt"
    js.write_text(json.dumps({**out, "reading_code": {"head": COMMIT, "clean": True}}))
    txt.write_text(text)
    assert er.reading_commit(js, txt) is None


def test_a_refused_arm_is_not_voided_by_its_missing_record(e58):
    recs = _recs(e58, GOOD, refused_on={("half", 2): "shape: x"})
    recs["half"][2]["facts"]["records"][ARM] = {}
    out = e58.reading(recs, [], COMMIT)
    assert not out["void"] and not out["passes"] and not out["clauses"]["P3"]


def test_the_barred_tasks_are_refused_by_name(e58):
    for name in ("e75", "e725"):
        with pytest.raises(SystemExit, match="barred"):
            e58._task(name)
    with pytest.raises(SystemExit, match="unknown task"):
        e58._task("nope")
    assert e58._task("half").key == e58.nt.all_tasks()["half"].key


# ------------------------------------------------------------------------------------------------ the J7 guard
def _e47_dir(e58, tmp_path, shift_floor=0.0, seed_shift=None):
    d = tmp_path / "e47tasks"
    d.mkdir(exist_ok=True)
    for name in e58.TASKS:
        key = e58.nt.all_tasks()[name].key
        rows = [{"seed": s, "scores": {"floor": 50.0 + shift_floor + (seed_shift or {}).get(s, 0.0),
                                        "fgw_shape0x": 48.0 + 0.1 * s}, "skills": {}} for s in SEEDS]
        (d / f"T2e_{name}.json").write_text(json.dumps({
            "task": {"key": key}, "provenance": {"commit": "9967cf1" + "0" * 33, "tree_clean": True},
            "build": "main", "protocol": "bake", "rows": rows}))
    return d


def _e55_reading(tmp_path, per_task, commit="4535f98c98ae1a2a8488185e8a99e4a2e8d92d86", void=False):
    p = tmp_path / "e55.json"
    p.write_text(json.dumps({"void": void, "commit": commit,
                             "tasks": {n: {"board_minus_floor": {O: v}} for n, v in per_task.items()}}))
    return p


def test_the_o_guard_against_e55_and_the_floor_against_e47(e58, tmp_path):
    names = list(e58.TASKS)
    p = _e55_reading(tmp_path, dict.fromkeys(names, -2.0))
    o_ref = e58.e55_o_per_task(p, sha=None)
    assert o_ref[1] == [] and o_ref[0] == dict.fromkeys(names, -2.0)
    d = _e47_dir(e58, tmp_path)
    key = e58.nt.all_tasks()["late"].key
    fl = e58.e47_scores(key, d)
    assert fl[1] == []
    floor = np.full(10, 50.0)
    good = {"floor": floor, O: floor - 2.0}
    g = e58.gate_check("late", key, good, o_ref, fl)
    assert g["pass"] and g["o_off_e55"] == pytest.approx(0.0, abs=1e-12)
    wobble = good[O].copy()
    wobble[3] += 0.05  # O per seed moves; its ten-seed mean by 0.005: no per-seed check of O
    assert e58.gate_check("late", key, {**good, O: wobble}, o_ref, fl)["pass"]
    assert not e58.gate_check("late", key, {**good, O: good[O] + 0.011}, o_ref, fl)["pass"]
    f2 = floor.copy()
    f2[4] += 0.011  # the floor per seed against E47
    assert not e58.gate_check("late", key, {"floor": f2, O: f2 - 2.0}, o_ref, fl)["pass"]
    assert not e58.gate_check("late", key, {"floor": floor + 0.011, O: floor - 1.989}, o_ref, fl)["pass"]
    void = e58.e55_o_per_task(_e55_reading(tmp_path, dict.fromkeys(names, -2.0), void=True), sha=None)
    assert void[1] and not e58.gate_check("late", key, good, void, fl)["pass"]
    other = e58.e55_o_per_task(_e55_reading(tmp_path, dict.fromkeys(names, -2.0), commit="c" * 40), sha=None)
    assert other[1]
    assert e58.e55_o_per_task(p, sha="0" * 64)[1]  # another reading than the pinned one


def test_the_committed_e55_reading_and_e47_files_are_the_pinned_ones(e58):
    ref, why = e58.e55_o_per_task()
    assert why == [] and set(ref) == set(e58.TASKS)
    for name, sha in e58.E47_SHA256.items():
        assert hashlib.sha256((e58.E47_DIR / "tasks" / name).read_bytes()).hexdigest() == sha
    keys = {e58.nt.all_tasks()[n].key for n in e58.TASKS}
    for name in e58.TASKS:
        assert e58.e47_scores(e58.nt.all_tasks()[name].key)[1] == []
    assert len(keys) == len(e58.E47_SHA256) == 7


def test_the_j7_dry_run_on_synthetic_records(e58, tmp_path):
    names = list(e58.TASKS)
    d = _e47_dir(e58, tmp_path)  # C0 minus the floor: mean of 48 + 0.1 s - 50 = -1.55
    e49 = tmp_path / "e49.json"

    def write49(v):
        e49.write_text(json.dumps({"void": False, "commit": "d76ccd5" + "0" * 33,
                                   "tasks": {n: {O: {"board": {"mean": v}}} for n in names}}))

    write49(-0.45)  # O minus C0: O minus the floor = -0.45 - 1.55 = -2.0
    p55 = _e55_reading(tmp_path, dict.fromkeys(names, -2.0))
    r = e58.j7_dry_run(p55, e49, d, e55_sha=None, e49_sha=None)
    assert r["pass"] and r["worst"] == pytest.approx(0.0, abs=1e-12) and set(r["per_task"]) == set(names)
    write49(-0.439)
    r = e58.j7_dry_run(p55, e49, d, e55_sha=None, e49_sha=None)
    assert not r["pass"] and r["worst"] == pytest.approx(0.011)
    r = e58.j7_dry_run(p55, tmp_path / "missing.json", d, e55_sha=None, e49_sha=None)
    assert not r["pass"] and r["why"]


# ------------------------------------------------------------------------------------------------ the row (e58_row)
def test_the_row_check(er, tmp_path):
    good = _rows(er)
    cases = [((good,), True), ((), False), ((good, good), False), ((_rows(er, cells=10),), False),
             ((_rows(er, cells=12),), False),
             ((_rows(er, first="E58 (gate G1 of row R59's file: it changes the mix; one counted try)"),), False),
             ((_rows(er, first=er.FIRST_CELL_OPENING.replace("it moves the overlay's form", "it moves the mix")),),
              False),
             ((good.replace("| E58 ", "| E58a ", 1),), False), (("| E58 | x |",), False)]
    for k, (rows, want) in enumerate(cases):
        p = _ledger(tmp_path, *rows, name=f"T{k}.md")
        assert (not er.row_problems(p)) is want, rows
        assert (er.main(["row", str(p)]) == 0) is want
    assert er.row_problems(tmp_path / "missing.md") and er.main(["row", str(tmp_path / "missing.md")]) == 1
    escaped = _rows(er, first=er.FIRST_CELL_OPENING + " (blocks 1\\|2,3))")  # an escaped pipe stays in its cell
    assert not er.row_problems(_ledger(tmp_path, escaped, name="esc.md"))
    assert "it moves the overlay's form" in er.FIRST_CELL_OPENING and "`fgw_shape0x_ov_t23`" in er.FIRST_CELL_OPENING


def test_the_outcome_cell_must_hold_the_pass_words_and_nothing_else(er, tmp_path):
    cases = [(_pass_cell(er), True), ("not yet run", False), (er.PASS_TEXT.replace("PASS", "pass"), False),
             (_pass_cell(er) + " " + er.STOP_TEXT, False), (_pass_cell(er) + " G1: VOID.", False),
             (er.STOP_TEXT, False), (er.IDLE_TEXT, False), (er.VOID_TEXT, False),
             (er.PASS_TEXT[:-1], False)]
    for k, (cell, want) in enumerate(cases):
        p = _ledger(tmp_path, _rows(er, outcome=cell), name=f"P{k}.md")
        assert er.passed(p) is want, cell
        assert (er.main(["passed", str(p)]) == 0) is want
    p = _ledger(tmp_path, _rows(er, outcome=_pass_cell(er), cells=12), name="bad.md")
    assert not er.passed(p)  # a malformed row never passes
    assert er.main(["nonsense"]) == 2


def test_the_reading_check(er, tmp_path):
    def write(k, **over):
        r = {"void": False, "outcome": "PASS", "outcome_text": er.PASS_TEXT, "passes": True, "commit": COMMIT,
             "reading_code": {"head": COMMIT, "clean": True}}
        txt = over.pop("_txt", f"Row E58 (run commit {COMMIT}): ...\n{er.PASS_TEXT}\nNothing here shows ...\n")
        r.update(over)
        js, t = tmp_path / f"r{k}.json", tmp_path / f"r{k}.txt"
        js.write_text(json.dumps(r))
        t.write_text(txt)
        return js, t

    assert er.reading_commit(*write(0)) == COMMIT
    bad = [{"void": True}, {"outcome": "STOP"}, {"outcome_text": er.STOP_TEXT}, {"passes": False},
           {"commit": COMMIT[:12]}, {"reading_code": {"head": COMMIT, "clean": False}},
           {"reading_code": {"head": "b" * 40, "clean": True}}, {"_txt": f"x {er.PASS_TEXT}\n"},
           {"_txt": f"{er.PASS_TEXT}\n{er.STOP_TEXT}\n"}]
    for k, over in enumerate(bad, start=1):
        assert er.reading_commit(*write(k, **over)) is None, over
    js = tmp_path / "list.json"
    js.write_text("[]")
    assert er.reading_commit(js, write(99)[1]) is None
    out = subprocess.run([sys.executable, "-I", "-B", str(REPO / ROW_REL), "reading", *map(str, write(100))],
                         capture_output=True, text=True)
    assert out.returncode == 0 and out.stdout.strip() == COMMIT


def test_the_row_module_runs_on_the_standard_library_alone():
    src = (REPO / ROW_REL).read_text()
    imports = set(re.findall(r"^(?:from|import) ([a-z_]+)", src, flags=re.M))
    assert imports <= {"__future__", "json", "re", "sys", "pathlib"}


# ------------------------------------------------------------------------------------------------ the read's steps
def test_the_steps_refuse_without_a_well_formed_row_or_on_an_unclean_tree(e58, er, tmp_path, monkeypatch):
    monkeypatch.setattr(e58, "TRIES", _ledger(tmp_path, name="none.md"))
    with pytest.raises(SystemExit, match="E58 has no well-formed row"):
        e58._refuse_unless_ready(COMMIT)
    with pytest.raises(SystemExit, match="E58 has no well-formed row"):
        e58.main(["read", "--commit", COMMIT, "--runs", str(tmp_path), "--out", str(tmp_path / "o")])
    monkeypatch.setattr(e58, "TRIES", _ledger(tmp_path, _rows(er), name="ok.md"))
    monkeypatch.setattr(e58, "repo_state", lambda repo: SimpleNamespace(unclean=("x.py",), commit=COMMIT))
    with pytest.raises(SystemExit, match="not clean"):
        e58._refuse_unless_ready(COMMIT)


def test_the_run_step_refuses_without_the_seven_gate_files(e58, er, tmp_path, monkeypatch):
    monkeypatch.setattr(e58, "_refuse_unless_ready", lambda commit: None)
    built = []
    monkeypatch.setattr(e58, "one_seed", lambda *a, **k: built.append(a))
    runs = tmp_path / "runs"
    for name in e58.TASKS[:-1]:  # six of seven pass
        p = e58._gate_path(runs, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"pass": True, "commit": COMMIT}))
    with pytest.raises(SystemExit, match="refused"):
        e58.main(["run", "--task", "late", "--commit", COMMIT, "--runs", str(runs), "--data", str(tmp_path)])
    e58._gate_path(runs, "x3").write_text(json.dumps({"pass": True, "commit": "b" * 40}))
    with pytest.raises(SystemExit, match="refused"):
        e58.main(["run", "--task", "late", "--commit", COMMIT, "--runs", str(runs), "--data", str(tmp_path)])
    assert built == [] and not (runs / "run").exists()
    with pytest.raises(SystemExit, match="barred"):
        e58.main(["run", "--task", "e75", "--commit", COMMIT, "--runs", str(runs), "--data", str(tmp_path)])


def test_the_read_refuses_to_overwrite(e58, tmp_path):
    (tmp_path / "reading.txt").write_text("x")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        e58.main(["read", "--commit", COMMIT, "--runs", str(tmp_path), "--out", str(tmp_path)])


def test_the_wrappers_stand_o_in_for_a_refused_arm_and_never_meet_a_test_board(e58):
    from vec.boards import Board
    from vec.evaluate import Stage

    board = Board(key="pseudo:T2e/x", task="T2", genes=["a"], needs_coords=True, split="pseudo", min_cells=1,
                  max_cells=50, truth_cells=10)
    stage = Stage(X=np.ones((5, 1), np.float32), coords=np.zeros((5, 3)))
    ctx = SimpleNamespace(seed=harness.BUILD_SEED_OFFSET + 3)
    o = harness.Candidate(O, frozenset({"T2"}), (), lambda c: stage, {})

    def guard(c):
        raise ValueError("pseudo: cloud guard flags ['coverage']")

    with e58.captured() as seen:
        seen["board"] = board
        assert e58.observed(o, seen).build(ctx) is stage
        out = e58.guarded(harness.Candidate(ARM, frozenset({"T2"}), (), guard, {}), seen).build(ctx)
        assert out is stage and seen["refused"][(ARM, 3)] == "cloud guard"
        with pytest.raises(ValueError, match="else"):
            e58.guarded(harness.Candidate(ARM, frozenset({"T2"}), (), lambda c: (_ for _ in ()).throw(
                ValueError("something else")), {}), seen).build(ctx)
    test = Board(**{**board.__dict__, "key": "T2:embryo:test_interp", "split": "test_interp"})
    with e58.captured() as seen:
        seen["board"] = test
        with pytest.raises(RuntimeError, match="test board"):
            e58.observed(o, seen).build(ctx)
    cands = e58.candidates({"board": board, "stage": {}, "refused": {}})
    assert cands[0] is harness.CANDIDATES[harness.FLOOR] and [c.name for c in cands] == ["floor", O, ARM]


# ------------------------------------------------------------------------------------------------ the launcher
def test_the_launcher_syntax_and_contents():
    sh = (REPO / LAUNCHER_REL).read_text()
    assert subprocess.run(["bash", "-n", str(REPO / LAUNCHER_REL)]).returncode == 0
    code = [line.split("#")[0] for line in sh.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    assert "export PATH=$HOME/.local/bin:$PATH" in sh
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < sh.index("python3 -I -B experiments/coords/e58_row.py")
    tasks = re.search(r"^TASKS=\(([^)]*)\)$", sh, flags=re.M).group(1).split()
    assert tasks == ["half", "late", "new1", "new2", "x1", "x2", "x3"] and "e75" not in sh.split("set -u", 1)[1]
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh and 'TMPDIR="$S/tmp"' in sh
    for line in sh.splitlines():
        if "$PY " in line and not line.strip().startswith("#"):
            assert "/dev/null" in line or line.rstrip().endswith("\\"), line
    assert "upload" not in " ".join(code)
    assert "tests/test_e58_mix_lean_read.py" in sh and "tests/test_mix_lean_build.py" in sh
    i = {k: sh.index(k) for k in ("finish 96", "finish 98", "finish 95", "begin tests")}
    assert i["finish 96"] < i["finish 98"] < i["finish 95"] < i["begin tests"]


def _home(tmp_path: Path, name: str = "home") -> Path:
    """A home whose ~/.local/bin holds a `uv` that fails (the launcher puts that folder first on PATH)."""
    home = tmp_path / name
    fake = home / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text("#!/bin/sh\necho \"$*\" >> \"$HOME/uvcalls\"\nexit 1\n")
    (fake / "uv").chmod(0o755)
    return home


def _clone(tmp_path: Path, rows: tuple[str, ...]) -> Path:
    """A throwaway repository holding the launcher, the row check and a ledger with `rows`, committed."""
    repo = tmp_path / "checkout"
    for rel in (LAUNCHER_REL, ROW_REL):
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, repo / rel)
    tries = repo / "experiments/harness/TRIES.md"
    tries.parent.mkdir(parents=True)
    tries.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in rows))
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "rows")
    return repo


def _launch(repo: Path, home: Path, *argv: str, trace=False) -> subprocess.CompletedProcess:
    env = {**git_env(), "HOME": str(home)}
    return subprocess.run(["bash", *(["-x"] if trace else []), LAUNCHER_REL, *argv], cwd=repo, env=env,
                          capture_output=True, text=True, timeout=60)


def _scratch(home: Path, run: str) -> Path:
    return home / "v" / "scratch" / "e58" / run


def test_the_launcher_usage_odd_names_and_a_reused_folder(er, tmp_path):
    repo = _clone(tmp_path, (_rows(er),))
    home = _home(tmp_path)
    assert _launch(repo, home, "c" * 40).returncode == 2
    assert _launch(repo, home, "c" * 40, "a/b").returncode == 90 and not (home / "v").exists()
    assert _launch(repo, home, "c" * 40, ".x").returncode == 90
    _scratch(home, "r1").mkdir(parents=True)
    out = _launch(repo, home, "c" * 40, "r1")
    assert out.returncode == 91 and "never reused" in out.stderr and list(_scratch(home, "r1").iterdir()) == []


def test_a_wrong_commit_is_refused_96_and_its_trace_runs_nothing_but_git(er, tmp_path):
    repo = _clone(tmp_path, (_rows(er),))
    home = _home(tmp_path)
    out = _launch(repo, home, "0" * 40, "w1", trace=True)
    assert out.returncode == 96 and "HEAD is not" in out.stderr and "finish 96" in out.stderr
    assert "export PATH=" in out.stderr  # the PATH export ran first, in the trace
    assert not re.search(r"^\++ (uv|python3?) ", out.stderr, flags=re.M)
    s = _scratch(home, "w1")
    assert (s / "preflight.exit").read_text().strip() == "96" and (s / "done.exit").read_text().strip() == "96"
    assert not (home / "uvcalls").exists()


def test_an_unclean_tree_is_refused_98(er, tmp_path):
    repo = _clone(tmp_path, (_rows(er),))
    (repo / "experiments/harness/TRIES.md").write_text("changed\n")
    home = _home(tmp_path)
    out = _launch(repo, home, git(repo, "rev-parse", "HEAD"), "u1", trace=True)
    assert out.returncode == 98 and "not clean" in out.stderr and "finish 98" in out.stderr
    assert not re.search(r"^\++ (uv|python3?) ", out.stderr, flags=re.M)
    assert (_scratch(home, "u1") / "done.exit").read_text().strip() == "98" and not (home / "uvcalls").exists()


@pytest.mark.parametrize("rows", ["none", "two", "cells", "opening"])
def test_a_missing_or_malformed_row_e58_is_refused_95(er, tmp_path, rows):
    good = _rows(er)
    ledger = {"none": (), "two": (good, good), "cells": (_rows(er, cells=10),),
              "opening": (_rows(er, first="E58 (gate G1; one counted try)"),)}[rows]
    repo = _clone(tmp_path, ledger)
    home = _home(tmp_path)
    out = _launch(repo, home, git(repo, "rev-parse", "HEAD"), "m1", trace=True)
    assert out.returncode == 95 and "row E58" in out.stderr and "finish 95" in out.stderr
    assert re.search(r"^\++ python3 -I -B experiments/coords/e58_row.py row ", out.stderr, flags=re.M)
    s = _scratch(home, "m1")
    assert (s / "preflight.exit").read_text().strip() == "95" and not (s / "tests.exit").exists()
    assert not (home / "uvcalls").exists()


def test_a_well_formed_row_passes_preflight_and_a_failing_uv_stops_the_tests_step(er, tmp_path):
    repo = _clone(tmp_path, (_rows(er),))
    home = _home(tmp_path)
    out = _launch(repo, home, git(repo, "rev-parse", "HEAD"), "p1")
    s = _scratch(home, "p1")
    assert out.returncode == 1, out.stderr
    assert (s / "preflight.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "1"
    assert (s / "done.exit").read_text().strip() == "1" and not (s / "gate.exit").exists()
    calls = (home / "uvcalls").read_text().splitlines()
    assert len(calls) == 1 and "pytest" in calls[0] and "tests/test_e58_mix_lean_read.py" in calls[0]


def test_the_launcher_and_the_read_name_the_same_seven_tasks(e58):
    sh = (REPO / LAUNCHER_REL).read_text()
    tasks = re.search(r"^TASKS=\(([^)]*)\)$", sh, flags=re.M).group(1).split()
    assert tuple(tasks) == e58.TASKS
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    out = subprocess.run([sys.executable, "-I", "-B", str(REPO / "experiments/coords/e58_mix_lean_read.py"), "run",
                          "--help"], capture_output=True, text=True, env=env)
    assert out.returncode == 0 and "half" in out.stdout
