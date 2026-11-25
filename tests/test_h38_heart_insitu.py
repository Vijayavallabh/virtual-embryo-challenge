"""Row H38's read script (experiments/harness/h38_heart_insitu.py) and row module (experiments/harness/h38_row.py), on
hand-made parts, gate records, ledgers and chord files: no data file is read, nothing is built or scored. The bars on
fold E10.5 in both worlds (S1 to S5), IDLE read as STOP, VOID decided first (every condition code holds), the arms'
chords (ISC, anti-ISC, placebos deranged in their level bins with no fixed point), the J5 floor's arithmetic, and the
row's statements as code reads them."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "harness"))
import h38_row  # noqa: E402


@pytest.fixture(scope="module")
def h38():
    spec = importlib.util.spec_from_file_location("h38_heart_insitu_t", REPO / "experiments/harness/h38_heart_insitu.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


COMMIT = "c" * 40
SHAS = {"E10.5": "a" * 64, "E9.5": "b" * 64, "E8.75": "d" * 64}


def make(h38, delta=None, terms=None, commit=COMMIT, tweak=None):
    """Parts and a gate record: every arm's score is C's plus delta[arm] (per seed a little spread), every open term's
    skill C's plus terms[arm][term]; shape terms equal; checks clean."""
    delta = delta or {}
    terms = terms or {}
    parts = {}
    for key in h38.TASK_KEYS:
        rows = []
        for s in h38.SEEDS:
            scores, skills = {}, {}
            for q in h38.TRUTHS:
                scores[q] = {h38.FLOOR: 40.0}
                skills[q] = {h38.FLOOR: {m: 0.4 for m in (*h38.TERMS, *h38.SHAPE_TERMS)}}
                for a in h38.BUILT:
                    d = delta.get(a, 0.0)
                    scores[q][a] = 50.0 + 0.1 * s + d + (0.02 * (s % 3) if d else 0.0)
                    skills[q][a] = {m: 0.5 + terms.get(a, {}).get(m, 0.0) for m in h38.TERMS}
                    skills[q][a].update({m: 0.7 for m in h38.SHAPE_TERMS})
            checks = {"same_draw": True, "pre_size_coords_equal": True, "post_size_coords_equal": True,
                      "n_cells": {a: 100 for a in h38.BUILT}, "lib_mean": {a: 60.0 for a in h38.BUILT},
                      "faults": {a: [] for a in h38.COUNTED}, "fixed_points": {p: 0 for p in h38.PLACEBOS}}
            rows.append({"seed": s, "scores": scores, "skills": skills, "groups": {}, "checks": checks,
                         "records": {}})
        parts[key] = [{"axis": h38.AXIS, "task": key, "rows": rows, "chords_sha256": SHAS[h38.REF_OF[key]],
                       "provenance": {"commit": commit, "tree_clean": True}}]
    gate = {"axis": h38.GATE_AXIS, "verdict": "PASS", "commit": COMMIT, "j5_floor": 0.05,
            "provenance": {"commit": COMMIT, "tree_clean": True}, "bars": {"S1": 0.30, "S2": 0.20},
            "rows": {k: [{"seed": r["seed"], "scores": json.loads(json.dumps(r["scores"]))}
                         for r in parts[k][0]["rows"]] for k in h38.TASK_KEYS}}
    if tweak:
        tweak(parts, gate)
    return parts, gate


GOOD = {"ISC": 0.6, "anti-ISC": -0.4, "PL1": 0.1, "PL2": 0.0, "PL3": -0.1}
DE_UP = {"ISC": {"de_score": 0.01, "de_direction": 0.02}}


def test_pass_needs_every_bar_in_both_worlds(h38):
    res = h38.reading(*make(h38, GOOD, DE_UP), COMMIT, SHAS)
    assert res["outcome"] == "PASS" and res["line"] == "OUTCOME H38: PASS." and not res["void"]
    assert all(all(v for k, v in res["bars"][w].items() if k != "idle") for w, _ in h38.WORLDS)


@pytest.mark.parametrize("delta,terms,bar", [
    ({**GOOD, "ISC": 0.25}, DE_UP, "S1"),                               # under +0.30
    ({**GOOD, "PL2": 0.5}, DE_UP, "S2"),                                # a placebo within 0.20 of ISC
    ({**GOOD, "anti-ISC": -0.05}, DE_UP, "S3"),                         # the opposite sign not losing
    (GOOD, {"ISC": {"de_score": 0.01, "de_direction": -0.01}}, "S5"),   # a de term against
])
def test_a_failed_bar_reads_stop(h38, delta, terms, bar):
    res = h38.reading(*make(h38, delta, terms), COMMIT, SHAS)
    assert res["outcome"] == "STOP" and res["line"] == "OUTCOME H38: STOP."
    assert not all(res["bars"][w][bar] for w, _ in h38.WORLDS)


def test_s4_reads_each_transition(h38):
    def tweak(parts, gate):
        for r in parts[h38.TASK_KEYS[1]][0]["rows"]:
            for q in h38.TRUTHS:
                r["scores"][q]["ISC"] = r["scores"][q]["C"] - 0.2
                r["scores"][q]["anti-ISC"] = r["scores"][q]["C"] - 0.9
        for r in parts[h38.TASK_KEYS[0]][0]["rows"]:
            for q in h38.TRUTHS:
                r["scores"][q]["ISC"] = r["scores"][q]["C"] + 1.5
    res = h38.reading(*make(h38, GOOD, DE_UP, tweak=tweak), COMMIT, SHAS)
    assert res["outcome"] == "STOP" and not res["bars"]["own"]["S4"]


def test_nothing_changes_reads_stop_idle(h38):
    res = h38.reading(*make(h38, {"ISC": 0.01, "anti-ISC": -0.01}), COMMIT, SHAS)
    assert res["outcome"] == "STOP (IDLE)" and res["line"] == "OUTCOME H38: STOP (IDLE)."
    assert "nothing is built" in res["follows"]


def test_the_j5_floor_raises_the_bar(h38):
    parts, gate = make(h38, GOOD, DE_UP)
    gate["bars"] = {"S1": 0.9, "S2": 0.9}
    assert h38.reading(parts, gate, COMMIT, SHAS)["outcome"] == "STOP"


@pytest.mark.parametrize("what", ["shape", "lib", "fault", "fixed", "commit", "clean", "seed", "gate_c", "draw",
                                  "chords", "nan", "gate_missing", "gate_fail", "cells"])
def test_void_is_decided_first(h38, what):
    def tweak(parts, gate):
        r = parts[h38.TASK_KEYS[0]][0]["rows"][3]
        if what == "shape":
            r["skills"]["own"]["ISC"]["d2_shape"] += 1e-6
        elif what == "lib":
            r["checks"]["lib_mean"]["PL2"] = 60.0 * 1.02
        elif what == "fault":
            r["checks"]["faults"]["ISC"] = ["the P-carrying type CM is reported 'skipped'"]
        elif what == "fixed":
            r["checks"]["fixed_points"]["PL1"] = 1
        elif what == "commit":
            parts[h38.TASK_KEYS[1]][0]["provenance"]["commit"] = "e" * 40
        elif what == "clean":
            parts[h38.TASK_KEYS[1]][0]["provenance"]["tree_clean"] = False
        elif what == "seed":
            parts[h38.TASK_KEYS[0]][0]["rows"].pop()
        elif what == "gate_c":
            gate["rows"][h38.TASK_KEYS[0]][3]["scores"]["own"]["C"] += 1e-9
        elif what == "draw":
            r["checks"]["same_draw"] = False
        elif what == "chords":
            parts[h38.TASK_KEYS[0]][0]["chords_sha256"] = "f" * 64
        elif what == "nan":
            r["scores"]["own"]["anti-ISC"] = float("nan")
        elif what == "gate_missing":
            gate["axis"] = "other"
        elif what == "gate_fail":
            gate["verdict"] = "FAIL"
        elif what == "cells":
            r["checks"]["n_cells"]["PL3"] = 99
    res = h38.reading(*make(h38, GOOD, DE_UP, tweak=tweak), COMMIT, SHAS)
    assert res["outcome"] == "VOID" and res["why"] and res["line"] == "OUTCOME H38: VOID."
    assert "bars" not in res


def test_no_g0_statement_is_void(h38):
    assert h38.reading(*make(h38, GOOD, DE_UP), COMMIT, None)["outcome"] == "VOID"


def test_variants(h38):
    rng = np.random.default_rng(0)
    n = 80
    P = {"cm": np.where(rng.random(n) < 0.5, rng.normal(0, 0.3, n), np.nan),
         "endocardium": np.where(rng.random(n) < 0.3, rng.normal(0, 0.3, n), np.nan)}
    levels = {k: rng.normal(0, 1, n) for k in ("cm", "endocardium", "epicardium", "vascular_endothelium")}
    var, fixed = h38.variants(P, levels)
    assert set(var) == set(h38.COUNTED) and all(v == 0 for v in fixed.values())
    for k in P:
        np.testing.assert_array_equal(var["ISC"][k], P[k])
        np.testing.assert_array_equal(var["anti-ISC"][k], -P[k])
        held = np.isfinite(P[k])
        for p in h38.PLACEBOS:
            assert np.array_equal(np.isfinite(var[p][k]), held)
            assert sorted(var[p][k][held]) == sorted(P[k][held])
            assert not np.any(var[p][k][held] == P[k][held])
    a, _ = h38.variants(P, levels)
    np.testing.assert_array_equal(a["PL2"]["cm"], var["PL2"]["cm"])  # fixed by its seed


def test_j5_floor_arithmetic(h38):
    rows = {}
    for key in h38.TASK_KEYS:
        rs = []
        for s in h38.SEEDS:
            sc = {q: {"category": 50.0 + (0.3 if s % 2 else -0.1), "control": 50.0} for q in h38.TRUTHS}
            sk = {q: {a: {m: 0.5 for m in (*h38.TERMS, *h38.SHAPE_TERMS)} for a in ("category", "control")}
                  for q in h38.TRUTHS}
            rs.append({"seed": s, "scores": sc, "skills": sk})
        rows[key] = rs
    d = np.array([0.3 if s % 2 else -0.1 for s in h38.SEEDS])
    want = 2 * d.std(ddof=1) / np.sqrt(len(d))
    assert h38.j5_floor(rows) == pytest.approx(want)
    assert h38.j5_floor(None) is None


def test_chord_shas_need_the_row_statement_and_the_committed_files(h38, tmp_path, monkeypatch):
    chords = tmp_path / "chords"
    chords.mkdir()
    shas = {}
    for ref in ("E10.5", "E9.5", "E8.75"):
        p = chords / f"chords_{ref}.json"
        p.write_text(json.dumps({"ref": ref}))
        shas[ref] = hashlib.sha256(p.read_bytes()).hexdigest()
    monkeypatch.setattr(h38, "CHORDS", chords)
    monkeypatch.setattr(h38, "tracked_clean", lambda p: True)
    line = f"| {h38_row.FIRST_CELL} | heart | " + " | ".join(["x"] * 8) + " | {cell} |\n"
    g0 = h38_row.g0_statement(shas["E10.5"], shas["E9.5"], shas["E8.75"])
    got, why = h38.chord_shas(line.format(cell=g0))
    assert got == shas and not why
    assert h38.chord_shas(line.format(cell="pending"))[0] is None
    assert h38.chord_shas(line.format(cell=g0 + " OUTCOME H38: STOP."))[0] is None
    (chords / "chords_E9.5.json").write_text("{}")
    assert h38.chord_shas(line.format(cell=g0))[0] is None
    monkeypatch.setattr(h38, "tracked_clean", lambda p: False)
    (chords / "chords_E9.5.json").write_text(json.dumps({"ref": "E9.5"}))
    assert h38.chord_shas(line.format(cell=g0))[0] is None


def test_the_row_module(tmp_path):
    first = h38_row.FIRST_CELL
    line = f"| {first} | heart | " + " | ".join(["x"] * 8) + " | {cell} |\n"
    g0 = h38_row.g0_statement("a" * 64, "b" * 64, "d" * 64)
    assert h38_row.before_g0(line.format(cell="pending")) == []
    assert h38_row.before_g0(line.format(cell=g0))
    assert h38_row.g0_shas(line.format(cell=g0))[0] == ("a" * 64, "b" * 64, "d" * 64)
    assert h38_row.g0_shas(line.format(cell=g0 + " " + g0))[0] is None
    assert h38_row.passed_sha256(line.format(cell=g0 + " OUTCOME H38: PASS."))[0] == "a" * 64
    for bad in (g0 + " OUTCOME H38: STOP.", "OUTCOME H38: PASS. " + g0, g0 + " OUTCOME H38: PASS. OUTCOME H38: VOID.",
                g0 + " OUTCOME H38: STOP (IDLE).", "OUTCOME H38: PASS."):
        assert h38_row.passed_sha256(line.format(cell=bad))[0] is None, bad
    assert h38_row.row_reasons(line.replace("| heart |", "| embryo |").format(cell="p"))
    assert h38_row.row_reasons(line.format(cell="a | b"))
    assert h38_row.row_reasons(line.format(cell="p") * 2)
    assert h38_row.row_reasons("| H38 (another wording) | heart | x |\n")
    p = tmp_path / "TRIES.md"
    p.write_text(line.format(cell=g0 + " OUTCOME H38: PASS."))
    assert h38_row.main(["passed", str(p)]) == 0 and h38_row.main(["check", str(p)]) == 1
    assert h38_row.main(["g0", str(p)]) == 1 and h38_row.main(["nonsense", str(p)]) == 2
    with pytest.raises(ValueError):
        h38_row.outcome_statement("KEPT")


def test_the_outcome_lines_are_the_row_modules(h38):
    for w in ("PASS", "STOP", "STOP (IDLE)", "VOID"):
        assert h38_row.OUTCOME_FORM.fullmatch(h38_row.outcome_statement(w))
        assert w in h38.FOLLOWS
