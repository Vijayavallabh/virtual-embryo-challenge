"""Row E59's read (experiments/coords/bf_read.py) and row check (experiments/coords/bf_row.py), on synthetic records:
the pass rule's clauses (a) to (f), the idle case, the void conditions, the guard, the fixed outcome words and the files
a PASS names, the G0 verdict, and the J7 dry run on committed records only. Nothing is built, read from data/ or
scored."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import bf_read as bf  # noqa: E402
import bf_row as br  # noqa: E402

C = "c" * 40
GROUPS = ("expression_change", "cell_state", "shape_growth", "local_organisation")
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


def _skills(base: float, shape: float = 0.0, nbhd: float = 0.0) -> dict:
    s = {m: base for m in METRICS}
    s["d2_shape"] += shape / 2
    s["occupancy_dice"] += shape / 2
    s["neighborhood_mmd"] += nbhd
    return s


def _record(task: str, s: int, gains: dict, prefix_gains: dict | None = None, **facts_over) -> dict:
    """One task-seed: every file's board, skills and groups, the arms at their base plus `gains[arm]` = (board, shape
    sum, neighbourhood) as functions of (task, seed)."""
    tasks = bf.nt.all_tasks()
    noise = 0.01 * np.sin(7 * s + len(task))
    scores = {"floor": 50.0, bf.O: 55.0 + noise, bf.SO: 56.0 + noise, bf.S2XS: 54.0 + noise}
    skills = {"floor": _skills(0.5), bf.O: _skills(0.6), bf.SO: _skills(0.62), bf.S2XS: _skills(0.58)}
    for arm in bf.ARMS:
        b, sh, nb = gains.get(arm, lambda t, s: (0.0, 0.0, 0.0))(task, s)
        base = bf.BASE_OF[arm]
        scores[arm] = scores[base] + b
        skills[arm] = {m: v for m, v in skills[base].items()}
        skills[arm]["d2_shape"] += sh / 2
        skills[arm]["occupancy_dice"] += sh / 2
        skills[arm]["neighborhood_mmd"] += nb
    groups = {n: dict.fromkeys(GROUPS, 0.6) for n in scores}

    def res(build, names, extra):
        row = {"seed": s, "scores": {n: scores[n] + extra.get(n, 0.0) for n in ("floor", *names)},
               "skills": {n: skills[n] for n in ("floor", *names)}, "groups": {n: groups[n] for n in ("floor", *names)}}
        return {"provenance": {"commit": C, "tree_clean": True}, "build": build, "protocol": "bake", "rows": [row],
                "board": {"split": "pseudo"}}

    pg = {} if prefix_gains is None else {a: f(task, s) for a, f in prefix_gains.items()}
    prior = "|".join(bf.prior_stages(tasks[task].target))
    recs = {n: {"target_rms": 6.0, "n_jittered": 0, "shape_interp": "overlay"} for n in (bf.O, bf.SO)}
    recs[bf.S2XS] = {"target_rms": 6.0, "n_jittered": 0}
    for arm in bf.ARMS:
        recs[arm] = {"target_rms": 6.0, "n_jittered": 0, "shape_interp": bf.FORM_OF[arm],
                     "shape_interp_prior": "before|after" if arm == bf.BF2_O else prior}
    facts = {"board_split": "pseudo", "prior": prior, "n_equal": dict.fromkeys(bf.ARMS, True),
             "x": {**dict.fromkeys((bf.O, bf.BF_O, bf.BF2_O, bf.BFNS_O), "xo"), bf.SO: "xso", bf.BF_SO: "xso",
                   bf.S2XS: "xs"},
             "refused": dict.fromkeys(bf.ARMS), "records": recs, "rms_off": dict.fromkeys(bf.ARMS, 0.0)}
    facts.update(facts_over)
    return {"task": task, "seed": s, "run_commit": C, "builds": {"main": res("main", bf.MAIN_NAMES, {}),
                                                                 "prefix": res("prefix", bf.PREFIX_NAMES, pg)},
            "facts": facts}


def _recs(gains: dict, prefix_gains=None, **over) -> dict:
    return {t: [_record(t, s, gains, prefix_gains, **over) for s in bf.SEEDS] for t in bf.TASKS}


def _ref(recs) -> tuple[dict, list]:
    """Committed guard figures that the synthetic run reproduces exactly."""
    tasks = bf.nt.all_tasks()
    ref = {"e47": {}, "o": {}, "so": {}}
    for name in bf.TASKS:
        rs = recs[name]
        ref["e47"][tasks[name].key] = {"floor": bf.series(rs, "floor", "board"),
                                       bf.S2XS: bf.series(rs, bf.S2XS, "board")}
        ref["o"][name] = float(np.mean(bf.series(rs, bf.O, "board") - bf.series(rs, "floor", "board")))
        ref["so"][name] = float(np.mean(bf.series(rs, bf.SO, "board") - bf.series(rs, bf.S2XS, "board")))
    return ref, []


G0_PASS = {"verdict": "PASS", "commit": C, "why": []}


def _good(board=1.0, shape=0.1, nbhd=0.0):
    return lambda t, s: (board + 0.05 * np.cos(s), shape, nbhd)


GOOD = {bf.BF_O: _good(), bf.BF_SO: _good(), bf.BF2_O: _good(0.5, 0.05), bf.BFNS_O: _good(0.5, 0.05)}


def _read(gains=GOOD, prefix_gains=None, g0=G0_PASS, ref=None, **over):
    recs = _recs(gains, prefix_gains, **over)
    return bf.reading(recs, [], C, g0, ref=_ref(recs) if ref is None else ref)


# ------------------------------------------------------------------------------------------------ the pass rule
def test_both_arms_pass_and_the_outcome_names_both_files():
    res = _read()
    assert not res["void"] and res["outcome"] == "PASS_BOTH" and res["outcome_text"] == br.PASS_TEXTS["PASS_BOTH"]
    assert res["files_passing"] == [br.FILE_O, br.FILE_SO] and res["tries_counted"] == 4
    assert all(all(v["clauses"].values()) for v in res["arms"].values())
    txt = bf.text(res)
    assert txt.splitlines()[-1] == br.PASS_TEXTS["PASS_BOTH"]


def test_a_failed_clause_e_stops_every_arm():
    gains = {**GOOD, bf.BF2_O: _good(1.0, 0.2)}  # the brackets-only control carries more shape than BF_O
    res = _read(gains)
    assert res["outcome"] == "STOP" and res["files_passing"] == [] and res["outcome_text"] == br.STOP_TEXT
    assert not res["clause_e"]["BF2_O"]["holds"] and res["clause_e"]["BFns_O"]["holds"]
    gains = {**GOOD, bf.BFNS_O: _good(1.0, 0.1)}  # a tie with the no-flip control is not ahead
    assert _read(gains)["outcome"] == "STOP"


@pytest.mark.parametrize("which, gain, outcome", [
    ("a", lambda t, s: (0.45 + 0.01 * np.cos(s), 0.1, 0.0), "PASS_SO"),  # under +0.5: (a) fails for BF_O
    ("b", lambda t, s: (1.0, 0.05, 0.0), "PASS_SO"),  # d2+Dice under +0.06
    ("b_truths", lambda t, s: (1.0, 0.2 if t in ("half", "late") else -0.01, 0.0), "PASS_SO"),  # 1 truth of 5
    ("c", lambda t, s: (1.0, 0.1, -0.03 if t == "e75" else 0.0), "PASS_SO"),  # one truth below -0.02
    ("f", lambda t, s: (-0.5 if t in ("half", "late") else 1.5, 0.1, 0.0), "PASS_SO"),  # E8.0's board not above 0
])
def test_each_clause_can_stop_bf_o_alone(which, gain, outcome):
    gains = {**GOOD, bf.BF_O: gain, bf.BF2_O: _good(0.3, 0.01), bf.BFNS_O: _good(0.3, 0.01)}
    res = _read(gains)
    assert res["outcome"] == outcome and res["files_passing"] == [br.FILE_SO], which
    clause = which[0]
    assert res["arms"]["BF_O"]["clauses"][clause] is False


def test_clause_d_reads_the_prefix_build():
    res = _read(prefix_gains={bf.BF_SO: lambda t, s: -0.5})
    assert res["outcome"] == "PASS_O" and res["arms"]["BF_SO"]["clauses"]["d"] is False
    assert res["outcome_text"] == br.PASS_TEXTS["PASS_O"]


def test_the_idle_case_reads_as_stop():
    idle = {**GOOD, bf.BF_O: _good(0.1, 0.1), bf.BF_SO: _good(-0.2, 0.1), bf.BF2_O: _good(0.0, 0.0),
            bf.BFNS_O: _good(0.0, 0.0)}
    res = _read(idle)
    assert res["outcome"] == "IDLE" and res["outcome_text"] == br.IDLE_TEXT and res["files_passing"] == []
    assert all(v["idle"] for v in res["arms"].values())
    nothing = {a: (lambda t, s: (0.0, 0.0, 0.0)) for a in bf.ARMS}  # the case where nothing changes
    res = _read(nothing)
    assert res["outcome"] == "IDLE" and "STOP" in res["outcome_text"]


# ------------------------------------------------------------------------------------------------ void
@pytest.mark.parametrize("over, words", [
    ({"refused": {**dict.fromkeys(bf.ARMS), bf.BF2_O: "cloud guard"}}, "refused"),
    ({"x": {**dict.fromkeys((bf.O, bf.BF2_O, bf.BFNS_O), "xo"), bf.BF_O: "other", bf.SO: "xso", bf.BF_SO: "xso",
            bf.S2XS: "xs"}}, "expression matrix"),
    ({"n_equal": {**dict.fromkeys(bf.ARMS, True), bf.BF_SO: False}}, "count"),
    ({"prior": "E6.75|E7.25|E7.5|E8.0|E8.25_late|E8.75|E9.5"}, "the truth must stay out"),
    ({"rms_off": {**dict.fromkeys(bf.ARMS, 0.0), bf.BF_O: 2e-6}}, "RMS radius"),
])
def test_the_void_conditions(over, words):
    res = _read(**over)
    assert res["void"] and res["outcome"] == "VOID" and res["outcome_text"] == br.VOID_TEXT
    assert any(words in w for w in res["why"]), res["why"][:3]
    assert bf.text(res).splitlines()[-1] == br.VOID_TEXT


def test_a_missing_seed_build_or_file_and_another_commit_are_void():
    recs = _recs(GOOD)
    ref = _ref(recs)
    recs["half"] = recs["half"][:9]
    assert bf.reading(recs, [], C, G0_PASS, ref=ref)["void"]
    recs = _recs(GOOD)
    del recs["x1"][3]["builds"]["prefix"]
    assert bf.reading(recs, [], C, G0_PASS, ref=ref)["void"]
    recs = _recs(GOOD)
    del recs["e75"][0]["builds"]["prefix"]["rows"][0]["scores"][bf.BF_SO]
    assert bf.reading(recs, [], C, G0_PASS, ref=ref)["void"]
    recs = _recs(GOOD)
    recs["new2"][5]["run_commit"] = "d" * 40
    assert bf.reading(recs, [], C, G0_PASS, ref=ref)["void"]
    recs = _recs(GOOD)
    recs["late"][2]["builds"]["main"]["rows"][0]["skills"][bf.BF_O]["mmd_u"] += 1e-6
    res = bf.reading(recs, [], C, G0_PASS, ref=ref)
    assert res["void"] and any("mmd_u" in w for w in res["why"])


def test_the_guard_voids_a_base_off_its_committed_figures():
    recs = _recs(GOOD)
    ref, _ = _ref(recs)
    for key, nudge in (("o", 0.02), ("so", -0.02)):
        bad = {**ref, key: {**ref[key], "new1": ref[key]["new1"] + nudge}}
        res = bf.reading(recs, [], C, G0_PASS, ref=(bad, []))
        assert res["void"] and any("new1" in w and "guard" in w for w in res["why"])
    k = bf.nt.all_tasks()["x2"].key
    bad = {**ref, "e47": {**ref["e47"], k: {**ref["e47"][k], bf.S2XS: ref["e47"][k][bf.S2XS] + 0.02}}}
    assert bf.reading(recs, [], C, G0_PASS, ref=(bad, []))["void"]
    ok = {**ref, "o": {**ref["o"], "new1": ref["o"]["new1"] + 0.009}}
    assert not bf.reading(recs, [], C, G0_PASS, ref=(ok, []))["void"]


def test_without_g0_passed_at_the_run_commit_the_read_is_void_and_a_g0_stop_counts_nothing():
    assert _read(g0=None)["void"] and _read(g0={**G0_PASS, "commit": "e" * 40})["void"]
    res = bf.reading({t: [] for t in bf.TASKS}, [], C, {"verdict": "STOP", "commit": C, "why": ["idle"]})
    assert res["outcome"] == "G0_STOP" and res["tries_counted"] == 0 and res["outcome_text"] == br.G0_STOP_TEXT
    assert bf.text(res).splitlines()[-1] == br.G0_STOP_TEXT


def test_g0s_verdict():
    def rec(share=0.85, ok=True, dice=0.6):
        return {"commit": C, "share": share, "guard_ok": ok, "guard_flags": [] if ok else ["block"],
                "idle_dice": dice, "z_extent_ratio": 1.0, "angle_thin_axis_z_deg": 10.0, "prior": "p"}

    good = {s: rec() for s in bf.G0_SEEDS}
    assert bf.g0_verdict(good, C)["verdict"] == "PASS"
    assert bf.g0_verdict({**good, bf.G0_SEEDS[3]: rec(share=0.69)}, C)["verdict"] == "STOP"
    assert bf.g0_verdict({**good, bf.G0_SEEDS[0]: rec(ok=False)}, C)["verdict"] == "STOP"
    idle = {s: rec(dice=0.95) for s in bf.G0_SEEDS}
    v = bf.g0_verdict(idle, C)
    assert v["verdict"] == "STOP" and v["idle"]
    assert bf.g0_verdict({s: r for s, r in good.items() if s != bf.G0_SEEDS[-1]}, C)["verdict"] == "STOP"
    assert bf.g0_verdict({**good, bf.G0_SEEDS[1]: {"commit": C, "error": "ShapeRefused: x"}}, C)["verdict"] == "STOP"
    assert bf.g0_verdict(good, "f" * 40)["verdict"] == "STOP"


def test_the_idle_dice_of_g0():
    import vec.shape_interp as si

    U = np.random.default_rng(0).normal(size=(3000, 3)) * [1.4, 1.0, 0.4]
    q = si.specimen_density(si.canonical_frame(U)[0])
    assert bf.idle_dice(U, q) > 0.8
    far = si.specimen_density(si.canonical_frame(np.random.default_rng(1).uniform(-1, 1, (3000, 3)) * [1, 1, 3])[0])
    assert bf.idle_dice(U, far) < bf.idle_dice(U, q)


# ------------------------------------------------------------------------------------------------ fixed choices
def test_the_merges_and_the_tasks():
    m = bf.merges()
    assert sorted(len(v) for v in m["five"].values()) == [1, 1, 2, 2, 3]
    assert m["e80"] == {bf.E80_TRUTH: ["half", "late"]} and m["five"][bf.E80_TRUTH] == ["half", "late"]
    assert set(bf.TASKS) == set(bf.nt.all_tasks()) and len(m["four"]) == 4
    assert bf.MAIN_NAMES == (bf.O, bf.SO, bf.S2XS, *bf.ARMS) and bf.PREFIX_NAMES == (bf.O, bf.SO, bf.BF_O, bf.BF_SO)
    assert (bf.BAR_A, bf.AHEAD, bf.BAR_B, bf.TRUTHS_B) == (0.5, 8, 0.06, 4)
    assert (bf.NB_MERGED, bf.NB_TRUTH, bf.PREFIX_TOL, bf.IDLE) == (-0.01, -0.02, 0.3, 0.3)
    assert (bf.G0_SHARE, bf.G0_IDLE_DICE, bf.GUARD_TOL) == (0.7, 0.9, 0.01)


def test_the_committed_reference_records_are_the_named_ones():
    assert bf.reference_problems() == []
    ref, why = bf.committed_figures()
    assert why == [] and set(ref["o"]) == set(ref["so"]) == set(bf.TASKS) and len(ref["e47"]) == 9


def test_the_j7_dry_run_on_committed_records(tmp_path):
    res = bf.j7_dry_run()
    assert res["pass"] and res["worst"] <= bf.GUARD_TOL and len(res["per_task"]) == 7
    moved = json.loads(bf.COMMITTED["e58"][0].read_text())
    moved["tasks"]["half"]["board_minus_floor"][bf.O] += 0.05
    p = tmp_path / "e58.json"
    p.write_text(json.dumps(moved))
    assert not bf.j7_dry_run(e58=p, check_sha=False)["pass"]  # the guard can fail there


# ------------------------------------------------------------------------------------------------ the row check
def _row(outcome="not yet run (placed before its run)", cells=11, first=None) -> str:
    first = br.FIRST_CELL_OPENING + " (synthetic))" if first is None else first
    return "| " + " | ".join([first, *[f"c{k}" for k in range(2, cells)], outcome]) + " |"


def _ledger(tmp_path, *rows) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in rows))
    return p


def test_the_row_check_and_the_files_a_pass_names(tmp_path):
    assert br.row_problems(_ledger(tmp_path, _row())) == []
    assert br.passed_files(_ledger(tmp_path, _row())) == ()
    for key, files in br.PASS_FILES.items():
        cell = f"Run at the run commit; in the read's fixed words: \"{br.PASS_TEXTS[key]}\""
        assert br.passed_files(_ledger(tmp_path, _row(cell))) == files
    for bad in (br.STOP_TEXT, br.IDLE_TEXT, br.G0_STOP_TEXT, br.VOID_TEXT, br.PASS_TEXTS["PASS_O"] + br.STOP_TEXT,
                br.PASS_TEXTS["PASS_O"] + br.PASS_TEXTS["PASS_SO"], br.PASS_TEXTS["PASS_BOTH"].replace("PASS", "Pass")):
        assert br.passed_files(_ledger(tmp_path, _row(bad))) == ()
    for rows in ((), (_row(), _row()), (_row(cells=10),), (_row(first="E59 (another gate)"),)):
        assert br.row_problems(_ledger(tmp_path, *rows))
    assert len(set(br.OUTCOME_TEXTS.values())) == len(br.OUTCOME_TEXTS)
    assert br.main(["row", str(_ledger(tmp_path, _row()))]) == 0
    assert br.main(["passed", str(_ledger(tmp_path, _row()))]) == 1


def test_the_committed_reading_check(tmp_path):
    def write(**over):
        r = {"void": False, "outcome": "PASS_O", "outcome_text": br.PASS_TEXTS["PASS_O"], "commit": C,
             "files_passing": [br.FILE_O], "reading_code": {"head": C, "clean": True}}
        txt = over.pop("_txt", f"synthetic\n{br.PASS_TEXTS['PASS_O']}\n")
        r.update(over)
        (tmp_path / "r.json").write_text(json.dumps(r))
        (tmp_path / "r.txt").write_text(txt)
        return br.reading_pass(tmp_path / "r.json", tmp_path / "r.txt")

    assert write() == (C, (br.FILE_O,))
    for over in ({"outcome": "STOP"}, {"void": True}, {"files_passing": [br.FILE_SO]}, {"commit": C[:7]},
                 {"reading_code": {"head": C, "clean": False}}, {"_txt": br.STOP_TEXT + "\n"},
                 {"outcome_text": br.PASS_TEXTS["PASS_SO"]}):
        assert write(**over) is None, over
    (tmp_path / "o.json").write_text(json.dumps({"outcome_text": br.IDLE_TEXT}))
    assert br.outcome_line(tmp_path / "o.json") == br.IDLE_TEXT
    (tmp_path / "o.json").write_text(json.dumps({"outcome_text": "GBF: PASS, trust me"}))
    assert br.outcome_line(tmp_path / "o.json") is None


def test_a_reading_written_by_the_read_passes_the_rows_reading_check(tmp_path, monkeypatch):
    res = {**_read(), "reading_code": {"head": C, "clean": True}}
    (tmp_path / "reading.json").write_text(json.dumps(res, default=float))
    (tmp_path / "reading.txt").write_text(bf.text(res))
    assert br.reading_pass(tmp_path / "reading.json", tmp_path / "reading.txt") == (C, (br.FILE_O, br.FILE_SO))
    assert br.outcome_line(tmp_path / "reading.json") == br.PASS_TEXTS["PASS_BOTH"]


def test_the_arms_hand_their_prior_to_their_one_interpolation_and_leave_the_harness_as_it_is(monkeypatch):
    seen = []

    def fake_ic(*a, **kw):
        seen.append(kw)
        return "res"

    def fake_fgw(ctx, **kw):
        return kw, bf.tb.interpolate_clouds(shape_form=kw["shape_form"])

    monkeypatch.setattr(bf.tb, "interpolate_clouds", fake_ic)
    monkeypatch.setattr(bf.harness, "_fgw", fake_fgw)
    holder = {"prior": "the task's prior"}
    for name in bf.ARMS:
        kw, out = bf.arm_candidate(name, holder).build(None)
        assert out == "res" and kw["expression"] == bf.EXPRESSION_OF[name] and kw["positions"] == "floor_exact"
        assert kw["t_pos"] == 0.0 and kw["shape_form"] == bf.FORM_OF[name]
    assert seen == [{"shape_form": "overlay_bayes", "shape_prior": "the task's prior"},
                    {"shape_form": "overlay_bayes", "shape_prior": "the task's prior"},
                    {"shape_form": "overlay_bayes_brackets"},
                    {"shape_form": "overlay_bayes_noflip", "shape_prior": "the task's prior"}]
    assert bf.tb.interpolate_clouds is fake_ic  # put back after each build
    assert "shape_prior" not in (REPO / "src/vec/harness.py").read_text()  # the harness's candidates are untouched


def test_the_launcher_lists_these_tests_and_the_tasks():
    sh = (REPO / "experiments/coords/run_bf.sh").read_text()
    assert "tests/test_bf_read.py" in sh and "tests/test_bf_form.py" in sh
    assert "TASKS=(half late new1 new2 e75 e725 x1 x2 x3)" in sh
    assert "SEEDS=(" + " ".join(str(s) for s in bf.G0_SEEDS) + ")" in sh
