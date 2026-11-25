"""Row H37's script (experiments/harness/h37_heart_world_blur.py) and its row module (experiments/harness/h37_row.py),
on hand-made records and temporary ledgers only: no data is read, nothing is built or scored.

The row's first cell and outcome words read back by code; the bars (B1, B2 labelling a PASS, IDLE and NOT PASSED as
STOP); G1 against row H35's parts; G0's verdict and words; every VOID condition the read holds; the refusals made
before anything loads."""

from __future__ import annotations

import copy
import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


h37 = _load("h37_heart_world_blur", REPO / "experiments" / "harness" / "h37_heart_world_blur.py")
hrow = h37.hrow
COMMIT = "a" * 40
SKILLS = ("de_score", "de_direction", "mmd_u", "variogram", "neighborhood_mmd", "d2_shape", "occupancy_dice",
          "scale_log_ratio")


def row_line(outcome="pending", first=hrow.FIRST_CELL) -> str:
    return f"| {first} | heart | x | x | x | x | x | x | x | x | {outcome} |\n"


# ------------------------------------------------------------------------------------------------ the row
def test_the_first_cell_is_fixed_and_holds_no_bar():
    assert "|" not in hrow.FIRST_CELL and hrow.FIRST_CELL.startswith("H37 (") and "six counted tries" in hrow.FIRST_CELL
    assert h37.FIRST_CELL == hrow.FIRST_CELL and h37.RULING_ROW == "R69" and h37.ROW == "H37"


def test_the_row_module_reads_the_row(tmp_path):
    t = tmp_path / "TRIES.md"
    t.write_text(row_line())
    assert hrow.main(["check", str(t)]) == 0 and hrow.main(["passed", str(t)]) == 4
    t.write_text(row_line() + row_line())
    assert hrow.main(["check", str(t)]) == 3
    t.write_text(row_line(first="H37 (another)"))
    assert hrow.main(["check", str(t)]) == 3
    t.write_text(row_line().replace(" x | x |", " x |", 1))
    assert hrow.main(["check", str(t)]) == 3
    for label in hrow.PASS_LABELS:
        t.write_text(row_line(f'In the read\'s words: "{h37.outcome_line(h37.PASS, label)}"'))
        assert hrow.main(["check", str(t)]) == 4 and hrow.main(["passed", str(t)]) == 0
        assert hrow.outcome(t.read_text()) == ("PASS", label)
    for line, code in ((h37.outcome_line(h37.STOP, h37.IDLE), 6), (h37.outcome_line(h37.STOP, h37.NOT_PASSED), 6),
                       (h37.outcome_line(h37.STOP, "G0", ["x"]), 6), (h37.outcome_line(h37.VOID, why=["v4: y"]), 5),
                       (h37.outcome_line(h37.PASS, "world pattern") + " " + h37.outcome_line(h37.STOP, h37.IDLE), 4),
                       ("OUTCOME H37: PASS (something else).", 4)):
        t.write_text(row_line(line))
        assert hrow.main(["passed", str(t)]) == code, line
    assert hrow.main(["passed", str(tmp_path / "none.md")]) == 3 and hrow.main(["other", "x"]) == 2


def test_the_ruling_row_is_read_by_its_first_cell(tmp_path):
    good = ("| R69 (the captain's order of 2026-11-25, ruling, no run: \"File WBH may use ...\"; built only after its "
            "gate row H37 reads PASS, only with `--captain-ruling R69`; not a try) | heart | x |\n")
    assert hrow.ruling_ok(good) and not hrow.ruling_ok(good + good)
    assert not hrow.ruling_ok(good.replace("ruling, no run", "ruling"))
    assert not hrow.ruling_ok(good.replace("`--captain-ruling R69`", "--captain-ruling R69"))
    assert not hrow.ruling_ok(good.replace("File WBH", "File W"))
    assert not hrow.ruling_ok(good.replace("R69 (the captain's order", "R69 (a note"))


def test_the_outcome_lines_are_fixed():
    assert h37.outcome_line(h37.PASS, "spread only") == "OUTCOME H37: PASS (spread only)."
    assert h37.outcome_line(h37.STOP, h37.IDLE) == "OUTCOME H37: STOP: IDLE."
    assert h37.outcome_line(h37.STOP, h37.NOT_PASSED) == "OUTCOME H37: STOP: NOT PASSED."
    assert h37.outcome_line(h37.STOP, "G0", ["a", "b"]) == "OUTCOME H37: STOP at G0: a; b."
    assert h37.outcome_line(h37.VOID, why=["v3: x"]) == "OUTCOME H37: VOID: v3: x; nothing is built."
    assert "heart_wbh_build.sh" in h37.follows_line(h37.PASS) and "nothing is built" in h37.follows_line(h37.STOP)


# ------------------------------------------------------------------------------------------------ the bars
def _tab(wc: dict, wp: dict, ahead: dict) -> dict:
    def st(m, a=10):
        return {"score": {"mean": m, "ahead": a, "behind": 10 - a, "lo": m, "hi": m, "se": 0.01}}
    return {w: {h37.FOLD: {h37.pname(h37.W, h37.CONTROL): st(wc[w], ahead[w]), h37.pname(h37.W, h37.P): st(wp[w])}}
            for w, _ in h37.WORLDS}


def test_the_bars():
    ws = [w for w, _ in h37.WORLDS]

    def both(v):
        return dict.fromkeys(ws, v)
    b = h37.bars(_tab(both(0.30), both(0.10), both(8)))
    assert (b["outcome"], b["label"]) == (h37.PASS, h37.WORLD_PATTERN)
    b = h37.bars(_tab(both(0.45), both(0.09), both(10)))
    assert (b["outcome"], b["label"]) == (h37.PASS, h37.SPREAD_ONLY)
    b = h37.bars(_tab(both(0.45), {ws[0]: 0.5, ws[1]: 0.0}, both(10)))
    assert b["label"] == h37.SPREAD_ONLY  # B2 in one world only
    b = h37.bars(_tab(both(0.45), both(0.5), {ws[0]: 7, ws[1]: 10}))
    assert (b["outcome"], b["label"]) == (h37.STOP, h37.NOT_PASSED)  # 7 of 10 seeds above 0 in one world
    b = h37.bars(_tab({ws[0]: 0.45, ws[1]: 0.29}, both(0.5), both(10)))
    assert b["label"] == h37.NOT_PASSED
    b = h37.bars(_tab(both(0.099), both(0.5), both(10)))
    assert (b["outcome"], b["label"]) == (h37.STOP, h37.IDLE)  # nothing changes: STOP
    b = h37.bars(_tab(both(-0.05), both(0.0), both(3)))
    assert b["label"] == h37.IDLE
    b = h37.bars(_tab({ws[0]: 0.10, ws[1]: 0.0}, both(0.0), both(5)))
    assert b["label"] == h37.NOT_PASSED  # |W - C| of exactly 0.10 is not under 0.10


# ------------------------------------------------------------------------------------------------ records
def blur_record(**over):
    return {"mode": "blocks", "recovery": {"exact": True}, "antithetic": True, "library_ratio_design": 1.001,
            "library_ratio": 1.05, "pb_max_abs_diff": 0.0004, "changed_cells": 0.9, **over}


def seed_row(seed: int, delta: dict, base: float = 60.0) -> dict:
    """One seed of a part: every arm's score in every truth, its skills, groups, records and checks."""
    scores, skills, groups = {}, {}, {}
    for q in h37.TRUTHS:
        scores[q] = {h37.FLOOR: 50.0, **{a: base + delta.get(a, 0.0) for a in h37.BUILT}}
        skills[q] = {a: dict.fromkeys(SKILLS, 0.5) for a in (h37.FLOOR, *h37.BUILT)}
        groups[q] = {a: {"g": 0.5} for a in (h37.FLOOR, *h37.BUILT)}
    return {"seed": seed, "scores": scores, "skills": skills, "groups": groups,
            "records": {h37.CONTROL: {"size": {}, "trend": {}}, **{a: blur_record() for a in h37.ARMS}},
            "checks": {"coords_equal": dict.fromkeys(h37.BUILT, True), "n_cells": dict.fromkeys(h37.BUILT, 100)}}


def parts_of(delta: dict, commit=COMMIT) -> dict:
    out = {}
    for key in h37.TASK_KEYS:
        out[key] = [{"axis": h37.AXIS, "task": key, "provenance": {"commit": commit, "tree_clean": True},
                     "rows": [seed_row(s, delta) for s in seeds]} for seeds in ((0, 1, 2, 3, 4), (5, 6, 7, 8, 9))]
    return out


def records(commit=COMMIT, parts=None):
    prov = {"commit": commit, "tree_clean": True}
    g0 = {"axis": h37.G0_AXIS, "commit": commit, "provenance": prov, "verdict": h37.PASS}
    parts = parts or parts_of({})
    rows = {key: [{"seed": r["seed"], "scores": r["scores"]} for p in ps for r in p["rows"]]
            for key, ps in parts.items()}
    gate = {"axis": h37.GATE_AXIS, "commit": commit, "provenance": prov, "verdict": h37.PASS, "rows": rows}
    return g0, gate


def test_a_clean_set_of_parts_reads_by_the_bars():
    parts = parts_of({h37.W: 0.5, h37.U: 0.2, h37.P: 0.3})
    g0, gate = records(parts=parts)
    assert h37.void_reasons(parts, g0, gate, COMMIT) == []
    res = h37.reading(parts, g0, gate, COMMIT)
    assert res["outcome"] == h37.PASS and res["label"] == h37.WORLD_PATTERN
    assert res["line"] == "OUTCOME H37: PASS (world pattern)." and "heart_wbh_build.sh" in res["follows"]
    assert "B1" in h37.text(res) and res["line"] in h37.text(res)
    json.dumps(res, default=float)
    parts = parts_of({h37.W: 0.02})
    res = h37.reading(parts, *records(parts=parts), COMMIT)
    assert res["line"] == "OUTCOME H37: STOP: IDLE." and "nothing is built" in res["follows"]


def _void(parts, g0=None, gate=None, commit=COMMIT):
    g, t = records(parts=parts)
    return h37.void_reasons(parts, g0 if g0 is not None else g, gate if gate is not None else t, commit)


@pytest.mark.parametrize("edit,word", [
    (lambda r: r["skills"]["own"][h37.W].update(d2_shape=0.6), "d2_shape"),
    (lambda r: r["records"][h37.U].update(library_ratio_design=1.03), "design library"),
    (lambda r: r["records"][h37.P].update(antithetic=False), "antithetic"),
    (lambda r: r["records"][h37.W]["recovery"].update(exact=False), "recovery"),
    (lambda r: r["records"][h37.W].update(pb_max_abs_diff=0.002), "mean X"),
    (lambda r: r["checks"]["coords_equal"].update(W=False), "coordinates"),
    (lambda r: r["checks"]["n_cells"].update(U=99), "cell numbers"),
    (lambda r: r["scores"]["own"].update(W=math.nan), "finite"),
    (lambda r: r["scores"]["own"].pop("P"), "missing from scores"),
])
def test_each_void_condition(edit, word):
    parts = parts_of({h37.W: 0.5})
    g0, gate = records(parts=copy.deepcopy(parts))
    edit(parts[h37.TASK_KEYS[0]][0]["rows"][1])
    why = h37.void_reasons(parts, g0, gate, COMMIT)
    assert any(word in w for w in why), why
    res = h37.reading(parts, g0, gate, COMMIT)
    assert res["outcome"] == h37.VOID and res["line"].startswith("OUTCOME H37: VOID: ")


def test_void_on_records_seeds_commits_and_the_tie():
    parts = parts_of({})
    g0, gate = records(parts=parts)
    assert any("no G0 record" in w for w in h37.void_reasons(parts, None, gate, COMMIT))
    assert any("gate said" in w for w in h37.void_reasons(parts, g0, {**gate, "verdict": h37.FAIL}, COMMIT))
    assert any("not of the run commit" in w for w in h37.void_reasons(parts, g0, gate, "b" * 40))
    short = copy.deepcopy(parts)
    short[h37.TASK_KEYS[1]][1]["rows"].pop()
    assert any("seeds" in w for w in _void(short))
    twice = copy.deepcopy(parts)
    twice[h37.TASK_KEYS[0]][1]["rows"].append(seed_row(0, {}))
    assert any("more than one part" in w for w in _void(twice))
    gone = {h37.TASK_KEYS[0]: parts[h37.TASK_KEYS[0]]}
    assert any("missing transition" in w for w in h37.void_reasons(gone, g0, gate, COMMIT))
    dirty = copy.deepcopy(parts)
    dirty[h37.TASK_KEYS[0]][0]["provenance"]["tree_clean"] = False
    assert any("unclean" in w for w in h37.void_reasons(dirty, g0, gate, COMMIT))
    tie = copy.deepcopy(gate)
    tie["rows"][h37.TASK_KEYS[0]][3]["scores"]["own"][h37.CONTROL] += 1e-12
    assert any("not the gate's" in w for w in h37.void_reasons(parts, g0, tie, COMMIT))
    assert h37.reading({}, g0, gate, COMMIT)["outcome"] == h37.VOID


# ------------------------------------------------------------------------------------------------ G1 and G0
def gate_rows(n_cells=100):
    rows = {}
    for key in h37.TASK_KEYS:
        rows[key] = []
        for s in h37.SEEDS:
            r = seed_row(s, {})
            rows[key].append({"seed": s, "n_cells": n_cells, "scores": r["scores"], "skills": r["skills"],
                              "groups": r["groups"], "outputs": {a: {"checks": dict.fromkeys(
                                  ("count recovery", "library"), True)} for a in h37.ARMS}})
    return rows


def committed_from(rows):
    out = {}
    for key, rs in rows.items():
        out[key] = [{"seed": r["seed"], "n_cells": {h37.H35_ARM: r["n_cells"]},
                     "scores": {q: {h37.H35_ARM: v[h37.CONTROL], h37.FLOOR: v[h37.FLOOR]}
                                for q, v in r["scores"].items()},
                     "skills": {q: {h37.H35_ARM: v[h37.CONTROL]} for q, v in r["skills"].items()},
                     "groups": {q: {h37.H35_ARM: v[h37.CONTROL]} for q, v in r["groups"].items()}} for r in rs]
    return json.loads(json.dumps(out))


def test_g1_holds_c_to_row_h35_exactly():
    rows = gate_rows()
    committed = committed_from(rows)
    assert h37.g1_reasons(rows, committed) == []
    assert h37.gate_verdict(rows, committed, [])[0] == h37.PASS
    nan = copy.deepcopy(rows)
    nan[h37.TASK_KEYS[0]][0]["skills"]["own"][h37.CONTROL]["de_score"] = math.nan
    assert h37.g1_reasons(nan, committed_from(nan)) == []  # a NaN equals a NaN
    off = copy.deepcopy(rows)
    off[h37.TASK_KEYS[1]][4]["scores"]["thin_0"][h37.CONTROL] += 1e-13
    why = h37.g1_reasons(off, committed)
    assert len(why) == 1 and "seed 4" in why[0] and h37.gate_verdict(off, committed, [])[0] == h37.FAIL
    fl = copy.deepcopy(rows)
    fl[h37.TASK_KEYS[0]][2]["scores"]["own"][h37.FLOOR] += 1
    assert "floor" in h37.g1_reasons(fl, committed)[0]
    assert "cell number" in h37.g1_reasons(gate_rows(101), committed)[0]
    assert h37.gate_verdict(rows, None, ["G1: x"])[0] == h37.FAIL
    bad = copy.deepcopy(rows)
    bad[h37.TASK_KEYS[0]][5]["outputs"][h37.U]["checks"]["library"] = False
    verdict, g1, practice = h37.gate_verdict(bad, committed, [])
    assert verdict == h37.STOP and g1 == [] and "seed 5 U: library" in practice[0]


def test_h35_committed_refuses_another_file(tmp_path):
    (tmp_path / "parts").mkdir()
    (tmp_path / "reading.json").write_text("{}")
    got, why = h37.h35_committed(tmp_path)
    assert got is None and any("sha256" in w for w in why)


def test_g0_reading_and_words():
    import numpy as np

    from vec import world_blur as wb

    prof = {n: np.array([0.3, -0.3]) for n in wb.PAIR_NAMES}
    prior = wb.Prior(("a", "b"), np.array([0.3, 0.3]), np.array([0.3, 0.3]), np.zeros(2, dtype=bool), prof,
                     {n: ["NCC", "Neural"] for n in wb.PAIR_NAMES})
    ok = dict.fromkeys(wb.OUTPUT_CHECKS, True)
    res = h37.g0_reading(prior, {"pass": True, "c": 1.0, "cut": 0.2}, ok)
    assert res["verdict"] == h37.PASS and res["line"] is None
    res = h37.g0_reading(prior, {"pass": False, "c": 0.1, "cut": 0.02}, {**ok, "dp rank": False})
    assert res["verdict"] == h37.STOP and res["failed"] == ["g00", "banked dp rank"]
    assert res["line"].startswith("OUTCOME H37: STOP at G0: G00's cut") and hrow.STOP_RE.search(res["line"])
    g = {"commit": COMMIT, "prior": prior.summary(), "g00": {"c": 0.1, "cut": 0.02}, **res}
    txt = h37.g0_text(g)
    assert "FAILS" in txt and res["line"] in txt
    assert set(h37.G0_CONDITION_WORDS) >= {*wb.prior_checks(prior), "g00", *(f"banked {k}" for k in ok)}


def test_the_arms_are_blurred_as_the_design_says():
    assert h37.arm_spec(h37.W) == {"mode": "blocks", "permute": False}
    assert h37.arm_spec(h37.U) == {"mode": "cells", "permute": False}
    assert h37.arm_spec(h37.P) == {"mode": "blocks", "permute": True}
    assert h37.FILE_SEED == 1000001 and h37.BANKED_NAME == "T2heart_organ_ts125_tk_half.h5ad"
    sums = (REPO / "experiments" / "builds" / "results" / "hearth_8a8bc3a_h2.SHA256SUMS").read_text()
    assert f"{h37.BANKED_SHA256}  tka/{h37.BANKED_NAME}" in sums
    assert h37.H35_RUN_COMMIT.startswith("53fc2674") and h37.H35_READING_SHA256.startswith("bf9746388ab5")


# ------------------------------------------------------------------------------------------------ refusals
def test_every_step_is_refused_without_an_open_row(tmp_path, monkeypatch):
    t = tmp_path / "TRIES.md"
    monkeypatch.setattr(h37, "TRIES", t)
    for text, word in (("", "exactly one row"), (row_line(first="H37 (x)"), "first cell"),
                       (row_line(h37.outcome_line(h37.STOP, h37.IDLE)), "already holds")):
        t.write_text(text)
        for call in (lambda: h37.refuse_g0(tmp_path / "g", COMMIT, None),
                     lambda: h37.refuse_gate(tmp_path / "g", COMMIT, tmp_path / "g0"),
                     lambda: h37.refuse_run(None, COMMIT, tmp_path / "g0", tmp_path / "gate"),
                     lambda: h37.refuse_read(tmp_path / "r", COMMIT)):
            with pytest.raises(SystemExit) as e:
                call()
            assert word in str(e.value)
    t.write_text(row_line())
    monkeypatch.setattr(h37, "git_state", lambda: {"commit": COMMIT, "tree_clean": False})
    with pytest.raises(SystemExit, match="tracked files"):
        h37.refuse_state(COMMIT)
    monkeypatch.setattr(h37, "git_state", lambda: {"commit": COMMIT, "tree_clean": True})
    with pytest.raises(SystemExit, match="HEAD"):
        h37.refuse_state("b" * 40)
    h37.refuse_state(COMMIT)  # nothing refused
    with pytest.raises(SystemExit, match="--banked"):
        h37.refuse_g0(tmp_path / "g", COMMIT, tmp_path / h37.BANKED_NAME)
    with pytest.raises(SystemExit, match="G0"):
        h37.refuse_gate(tmp_path / "g", COMMIT, tmp_path / "g0")
    (tmp_path / "g0").mkdir()
    g0, gate = records()
    (tmp_path / "g0" / h37.G0_JSON).write_text(json.dumps(g0))
    h37.refuse_gate(tmp_path / "g", COMMIT, tmp_path / "g0")
    with pytest.raises(SystemExit, match="gate"):
        h37.refuse_run(None, COMMIT, tmp_path / "g0", tmp_path / "gate")
    (tmp_path / "gate").mkdir()
    (tmp_path / "gate" / h37.GATE_JSON).write_text(json.dumps(gate))
    h37.refuse_run(None, COMMIT, tmp_path / "g0", tmp_path / "gate")
    (tmp_path / "r").mkdir()
    (tmp_path / "r" / "reading.json").write_text("{}")
    with pytest.raises(SystemExit, match="one read"):
        h37.refuse_read(tmp_path / "r", COMMIT)
