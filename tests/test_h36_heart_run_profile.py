"""Row H36 (experiments/harness/h36_heart_run_profile.py, experiments/harness/h36_row.py): the heart's XPC read, on
synthetic arrays, hand-made tables and temporary ledgers only.

The row module: a well-formed row and every malformed one; the outcome cell read by code (the practice read's PASSED AT
PRACTICE, then T's PASS with its sha256, nothing that ends the line). The script: the placebos (eligible genes only,
within level bins, no fixed point); the bars B1 to B4 in both worlds, each failing alone; the void conditions, each
alone; the outcome statements; the arms on file D's synthetic context (C is OTS125's practice form, a zero profile is C
bit for bit, the arms keep C's draw, coordinates and cell number, X's profile is delta_p); the gate's G2 and its
refusal when a committed record is not at hand; the steps E, Q and T end to end on synthetic groups with their
refusals stubbed (T refused unless the read's outcome is committed and stated in the row); row H34's committed parts
on this transition are the files pinned (read by code; no figure printed). Nothing is read from data/."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_heart_trend_line_d import _ctx_d  # noqa: E402

COMMIT = "c" * 40
SHA = "d" * 64


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def h36():
    return _load("h36_for_tests", REPO / "experiments" / "harness" / "h36_heart_run_profile.py")


@pytest.fixture(scope="module")
def row():
    return _load("h36_row_for_tests", REPO / "experiments" / "harness" / "h36_row.py")


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


def ledger(row, outcome: str = "pending", first: str | None = None, cells: int | None = None) -> str:
    c = [first or row.FIRST_CELL, "heart", *[f"cell {i}" for i in range(3, 11)], outcome]
    if cells is not None:
        c = c[:cells] if cells < len(c) else c + ["x"] * (cells - len(c))
    return "| R56 (ruling) | heart | x |\n| " + " | ".join(c) + " |\n| H35 (other) | heart | x |\n"


# ------------------------------------------------------------------------------------------------ the row
def test_a_well_formed_row_and_every_malformed_one(row):
    assert row.row_reasons(ledger(row)) == []
    assert row.outcome_cell(ledger(row, "the outcome")) == "the outcome"
    assert row.row_reasons("| H34 (x) | heart |\n") == ["no row H36 in TRIES.md"]
    two = ledger(row) + "| " + row.FIRST_CELL + " | heart |\n"
    assert row.row_reasons(two) == ["2 rows open with H36, not one"]
    assert "10 cells" in row.row_reasons(ledger(row, cells=10))[0]
    assert "12 cells" in row.row_reasons(ledger(row, cells=12))[0]
    assert "first cell" in row.row_reasons(ledger(row, first=row.FIRST_CELL.replace("R56", "R57")))[0]
    assert row.row_reasons(ledger(row).replace("| heart | cell 3", "| T1 | cell 3")) == [
        "row H36's second cell is not 'heart'"]
    assert row.row_reasons(ledger(row).replace("pending |", "pending")) == ["row H36 does not close with |"]
    assert row.outcome_cell(ledger(row, cells=10)) is None
    assert row.FIRST_CELL.startswith("H36 (one row, one script; eleven counted heart arms on one backtest")


def test_the_outcome_cell_read_by_code(row, tmp_path):
    p = row.PASSED_AT_PRACTICE
    ok = f"Run at ... In the read's fixed words: \"{p}\" Then at the later commit: \"{row.pass_statement(SHA)}\""
    assert row.passed_sha256(ledger(row, ok)) == (SHA, [])
    cases = {
        "pending": "0 PASS statements",
        p: "0 PASS statements",
        row.pass_statement(SHA): "PASSED AT PRACTICE",
        f"{row.pass_statement(SHA)} {p}": "stands before",
        f"{p} {row.pass_statement(SHA)} {row.pass_statement('e' * 64)}": "2 PASS statements",
        f"{p} Try H36: T STOPPED (noise): OXP125 is not built. {row.pass_statement(SHA)}": "T STOPPED",
        f"{p} {row.pass_statement(SHA)} Try H36: VOID: x; nothing is built.": "VOID",
        f"Try H36: NOT PASSED (B1 in own): OXP125 is not built. {row.pass_statement(SHA)} {p}": "NOT PASSED",
        f"{p} {row.pass_statement(SHA[:-1])}": "0 PASS statements",
    }
    for cell, why in cases.items():
        sha, reasons = row.passed_sha256(ledger(row, cell))
        assert sha is None and any(why in r for r in reasons), (cell, reasons)
    assert row.passed_sha256(ledger(row, ok, cells=10))[0] is None
    f = tmp_path / "TRIES.md"
    f.write_text(ledger(row, ok))
    assert row.main(["check", str(f)]) == 0 and row.main(["passed", str(f)]) == 0
    f.write_text(ledger(row, p))
    assert row.main(["check", str(f)]) == 0 and row.main(["passed", str(f)]) == 1
    f.write_text(ledger(row, cells=9))
    assert row.main(["check", str(f)]) == 1 and row.main(["passed", str(f)]) == 1
    assert row.main(["check"]) == 2 and row.main(["check", str(tmp_path / "none.md")]) == 1


def test_the_launchers_read_the_same_module(row):
    run = (REPO / "experiments" / "harness" / "run_h36.sh").read_text()
    build = (REPO / "experiments" / "builds" / "heart_xpc_build.sh").read_text()
    assert "python3 -I -B experiments/harness/h36_row.py check experiments/harness/TRIES.md" in run
    assert "python3 -I -B experiments/harness/h36_row.py passed experiments/harness/TRIES.md" in build
    assert "export PATH=$HOME/.local/bin:$PATH" in run and "export PATH=$HOME/.local/bin:$PATH" in build
    assert "--run-profile $SHA --z-mode raw --captain-ruling R56" in build


# ------------------------------------------------------------------------------------------------ the placebos
def test_the_placebos_derange_eligible_genes_within_level_bins(h36):
    rng = np.random.default_rng(0)
    n = 120
    level = rng.normal(size=n)
    eligible = rng.random(n) < 0.8
    d = np.where(eligible, rng.normal(0, 0.5, n), 0.0)
    deltas, fixed = h36.arm_deltas(d, eligible, level)
    assert set(deltas) == {h36.X, h36.ANTI, *h36.PLACEBOS} and len(h36.PLACEBOS) == 9
    np.testing.assert_array_equal(deltas[h36.X], d)
    np.testing.assert_array_equal(deltas[h36.ANTI], -d)
    assert all(v == 0 for v in fixed.values())
    idx = np.flatnonzero(eligible)
    bins = np.array_split(idx[np.argsort(level[idx], kind="stable")], h36.PERM_BINS)
    for name, seed in h36.PLACEBO_SEEDS.items():
        p = h36.placebo_perm(level, eligible, seed)
        assert sorted(p) == list(range(n)) and np.all(p[~eligible] == np.flatnonzero(~eligible))
        assert np.all(p[eligible] != idx)
        for b in bins:
            assert set(p[b]) == set(b)  # within its level bin
        np.testing.assert_array_equal(deltas[name], d[p])
        np.testing.assert_array_equal(p, h36.placebo_perm(level, eligible, seed))
    assert len({tuple(h36.placebo_perm(level, eligible, s)) for s in h36.PLACEBO_SEEDS.values()}) == 9
    assert list(h36.PLACEBO_SEEDS.values()) == list(range(20261251, 20261260))


# ------------------------------------------------------------------------------------------------ the bars
def rows_of(h36, x=1.0, anti=-1.0, p3=0.0, x_de=0.6, worlds=("own", "CLEAN"), shape_off=0.0):
    """Hand-made rows: in the worlds named, X is C + x, anti-X C + anti, P3 moved by p3 and X's de_score x_de; in the
    other world every bar holds."""
    out = []
    for s in h36.SEEDS:
        base = 50.0 + 0.1 * s
        scores, skills = {}, {}
        for q in h36.TRUTHS:
            in_w = (q == "own" and "own" in worlds) or (q != "own" and "CLEAN" in worlds)
            v = {h36.FLOOR: 50.0, h36.CONTROL: base, h36.X: base + (x if in_w else 1.0),
                 h36.ANTI: base + (anti if in_w else -1.0)}
            v.update({p: base + 0.01 * (i - 4) + (p3 if p == "P3" and in_w else 0.0)
                      for i, p in enumerate(h36.PLACEBOS)})
            scores[q] = v
            skills[q] = {a: {"de_score": (x_de if a == h36.X and in_w else 0.6 if a == h36.X else 0.5),
                             "de_direction": 0.6 if a == h36.X else 0.5,
                             **{m: 0.7 + (shape_off if a == h36.X else 0.0) for m in h36.SHAPE_TERMS}}
                         for a in v}
        out.append({"seed": s, "scores": scores, "skills": skills,
                    "checks": {"zero_profile_is_c": True, "same_draw": True, "pre_size_coords_equal": True,
                               "post_size_coords_equal": True, "n_cells": {a: 100 for a in h36.BUILT},
                               "mean_library_ratio": {a: 1.0 for a in h36.BUILT}}})
    return out


def test_the_bars_pass_together_and_each_fails_alone(h36):
    b = h36.bars(rows_of(h36))
    assert all(b[w][k]["holds"] for w in ("own", "CLEAN") for k in ("B1", "B2", "B3", "B4"))
    assert b["own"]["B2"]["need"] == pytest.approx(0.3) and b["own"]["B1"]["above"] == 10
    for kw, fails in (({"x": 0.29}, {"B1", "B2"}), ({"anti": -0.29}, {"B3"}), ({"x_de": 0.4}, {"B4"}),
                      ({"p3": 1.5}, {"B2"})):
        for w in ("own", "CLEAN"):
            b = h36.bars(rows_of(h36, worlds=(w,), **kw))
            assert {k for k in ("B1", "B2", "B3", "B4") if not b[w][k]["holds"]} == fails, (kw, w)
            other = "CLEAN" if w == "own" else "own"
            assert all(b[other][k]["holds"] for k in ("B1", "B2", "B3", "B4"))
    # B2's bar rises with the placebos' spread: two SDs of their ten-seed means
    rows = rows_of(h36)
    for r in rows:
        for q in h36.TRUTHS:
            for i, p in enumerate(h36.PLACEBOS):
                r["scores"][q][p] += 0.4 * (i - 4)
    b = h36.bars(rows)
    assert b["own"]["B2"]["need"] == pytest.approx(2 * 0.41 * np.std(np.arange(9) - 4, ddof=1))
    assert not b["own"]["B2"]["holds"]


def parts_of(h36, rows, commit=COMMIT, sha=SHA, **over):
    pc = {"x_is_delta_p": True, "anti_is_minus_delta_p": True,
          "placebo_fixed_points": {p: 0 for p in h36.PLACEBOS}, **over.pop("pc", {})}
    return [{"axis": h36.AXIS, "task": h36.TASK, "rows": rows[i:i + 2], "practice_profile_sha256": sha,
             "profile_checks": pc, "provenance": {"commit": commit, "tree_clean": True}, **over}
            for i in range(0, 10, 2)]


def rec_of(h36, step, verdict="PASS", commit=COMMIT, **kw):
    return {"row": "H36", "step": step, "verdict": verdict, "commit": commit,
            "provenance": {"commit": commit, "tree_clean": True}, **kw}


def test_the_reading_says_passed_not_passed_or_void(h36, row):
    gate, prof = rec_of(h36, "gate"), rec_of(h36, "profile", practice_profile_sha256=SHA)
    res = h36.reading(parts_of(h36, rows_of(h36)), gate, prof, COMMIT)
    assert res["outcome"] == h36.PASSED and res["line"] == row.PASSED_AT_PRACTICE and not res["void"]
    assert "priced question" in res["follows"] and "(T) runs" in res["follows"]
    res = h36.reading(parts_of(h36, rows_of(h36, worlds=("CLEAN",), anti=0.0)), gate, prof, COMMIT)
    assert res["line"] == "Try H36: NOT PASSED (B3 in CLEAN): OXP125 is not built."
    res = h36.reading(parts_of(h36, rows_of(h36)), rec_of(h36, "gate", "FAIL"), prof, COMMIT)
    assert res["outcome"] == h36.VOID and res["line"].startswith("Try H36: VOID: v: the gate step said 'FAIL'")
    assert res["line"].endswith("; nothing is built.") and "bars" not in res
    assert h36.text(res).startswith(res["line"])


def test_every_void_condition_alone(h36):
    gate, prof = rec_of(h36, "gate"), rec_of(h36, "profile", practice_profile_sha256=SHA)
    assert h36.void_reasons(parts_of(h36, rows_of(h36)), gate, prof, COMMIT) == []

    def why(parts=None, g=gate, p=prof, commit=COMMIT, unreadable=()):
        return h36.void_reasons(parts if parts is not None else parts_of(h36, rows_of(h36)), g, p, commit,
                                unreadable)

    rows = rows_of(h36)
    cases = {
        "no gate record": why(g=None), "the gate step said": why(g=rec_of(h36, "gate", "FAIL")),
        "not of the run commit": why(g=rec_of(h36, "gate", commit="e" * 40)),
        "the profile step said": why(p=rec_of(h36, "profile", "STOP", practice_profile_sha256=SHA)),
        "cannot be read": why(unreadable=["x/h36.json"]), "no part": why(parts=[]),
        "more than one part": why(parts=parts_of(h36, rows) + parts_of(h36, rows)[:1]),
        "not [0, 1, 2": why(parts=parts_of(h36, rows)[1:]),
        "more than one commit": why(parts=parts_of(h36, rows)[:4] + parts_of(h36, rows, commit="e" * 40)[4:]),
        "other than the run commit": why(commit="e" * 40),
        "unclean tree": why(parts=[{**p, "provenance": {"commit": COMMIT, "tree_clean": False}}
                                   for p in parts_of(h36, rows)]),
        "sha256 the profile step recorded": why(parts=parts_of(h36, rows, sha="e" * 64)),
        "X's profile is not delta_p": why(parts=parts_of(h36, rows, pc={"x_is_delta_p": False})),
        "leaves a gene on its own value": why(parts=parts_of(
            h36, rows, pc={"placebo_fixed_points": {**{p: 0 for p in h36.PLACEBOS}, "P2": 1}})),
        "d2_shape is not C's": why(parts=parts_of(h36, rows_of(h36, shape_off=2e-9))),
    }
    for k, v in (("zero_profile_is_c", False), ("same_draw", False), ("post_size_coords_equal", False)):
        bad = rows_of(h36)
        bad[3]["checks"][k] = v
        cases[k.replace("_", " ") + " fails"] = why(parts=parts_of(h36, bad))
    bad = rows_of(h36)
    bad[0]["checks"]["n_cells"][h36.X] = 99
    cases["cell numbers are not C's"] = why(parts=parts_of(h36, bad))
    bad = rows_of(h36)
    bad[0]["checks"]["mean_library_ratio"]["P4"] = 1.011
    cases["more than 1% from C's"] = why(parts=parts_of(h36, bad))
    bad = rows_of(h36)
    del bad[2]["scores"]["thin_0"]["P9"]
    cases["missing from scores"] = why(parts=parts_of(h36, bad))
    bad = rows_of(h36)
    bad[2]["scores"]["own"][h36.X] = float("nan")
    cases["not a finite number"] = why(parts=parts_of(h36, bad))
    for want, got in cases.items():
        assert got and any(want in g for g in got), (want, got)
    assert len(cases["d2_shape is not C's"]) == len(h36.SEEDS) * len(h36.TRUTHS) * len(h36.SHAPE_TERMS)


# ------------------------------------------------------------------------------------------------ the arms
def test_the_arms_on_file_ds_context(h36, bt):
    """C is OTS125's practice form (h28.trend_expression at 1.25, then D_f's size and half); a zero profile is C bit for
    bit; every arm keeps C's draw, coordinates and cell number and only the expression moves; X's and anti-X's
    expression differ from C's in opposite directions on the genes the profile moves most."""
    ctx = _ctx_d(bt)
    G = len(ctx.board.genes)
    rng = np.random.default_rng(3)
    eligible = np.ones(G, dtype=bool)
    d = np.clip(rng.normal(0, 0.6, G), -1.0, 1.0)
    deltas, _ = h36.arm_deltas(d, eligible, rng.normal(size=G))
    target = {"target": 7.0, "base": 5.6, "ratio": 1.25, "source": "synthetic"}
    out, extra = h36.build_arms(ctx, target, deltas)
    c = extra["checks"]
    assert c["zero_profile_is_c"] and c["same_draw"] and c["pre_size_coords_equal"] and c["post_size_coords_equal"]
    assert set(c["n_cells"]) == set(h36.BUILT) and len(set(c["n_cells"].values())) == 1
    assert all(abs(v - 1) <= h36.DEPTH_TOL for v in c["mean_library_ratio"].values())
    st, _, _, _ = h36.h28.trend_expression(ctx, 1.25)
    ref, _ = h36.h34.size_and_half(ctx, st, target, h36.CONTROL)
    np.testing.assert_array_equal(out[h36.CONTROL].X, ref.X)
    np.testing.assert_array_equal(out[h36.CONTROL].coords, ref.coords)
    assert extra["trend"]["strength"] == 1.25 and extra["trend"]["placebo"] is False
    for name in h36.COUNTED:
        assert not np.array_equal(out[name].X, out[h36.CONTROL].X), name
    up = d > 0.5  # over-detected in the run: X moves them down, anti-X up
    mean = {k: np.asarray(out[k].X)[:, up].mean() for k in (h36.CONTROL, h36.X, h36.ANTI)}
    assert mean[h36.X] < mean[h36.CONTROL] < mean[h36.ANTI]


# ------------------------------------------------------------------------------------------------ the gate
def test_the_gate_g2_and_its_records_not_at_hand(h36):
    rows = [{"seed": s, "scores": {q: {h36.COPY: 50.0, h36.T10: 51.0 + s, h36.T15: 52.0, h36.FLOOR: 50.0}
                                   for q in h36.TRUTHS}, "draw": {}} for s in h36.SEEDS]
    why, found = h36.gate_reasons(rows, {}, None, None, [dict(r) for r in rows])
    assert any("row H29b's reading.json is not at hand" in w for w in why)
    assert found["G2"]["unequal"] == [] and not any(w.startswith("G2") for w in why)
    theirs = [json.loads(json.dumps(r)) for r in rows]
    theirs[4]["scores"]["thin_0"][h36.T15] += 1e-12
    why, found = h36.gate_reasons(rows, {}, None, None, theirs)
    assert found["G2"]["unequal"] == ["trend_s1.5 CLEAN"] and any(w.startswith("G2:") for w in why)
    why, _ = h36.gate_reasons(rows, {}, None, None, None)
    assert "G2: row H34's committed parts are not at hand as the files pinned here" in why
    why, _ = h36.gate_reasons(rows[:9], {}, None, None, rows)
    assert why == ["G1: the gate did not build every seed"]


def test_row_h34s_committed_parts_on_this_transition_are_the_files_pinned(h36, tmp_path):
    rows = h36.h34_committed_rows()
    assert rows is not None and [r["seed"] for r in rows] == list(h36.SEEDS)
    assert all({h36.T10, h36.T15} <= set(r["scores"]["own"]) for r in rows)
    means = {n: h36.world_means(rows, n) for n in (h36.T10, h36.T15)}  # formed, never printed
    assert all(np.isfinite(v) for m in means.values() for v in m.values())
    for seeds in h36.H34_PART_SHA256:
        (tmp_path / seeds).mkdir()
        src = (h36.H34_PARTS / seeds / "h34.json").read_bytes()
        (tmp_path / seeds / "h34.json").write_bytes(src if seeds != "4-5" else src + b" ")
    assert h36.h34_committed_rows(tmp_path) is None


# ------------------------------------------------------------------------------------------------ the refusals
def test_every_step_is_refused_without_a_well_formed_row(h36, row, monkeypatch, tmp_path):
    f = tmp_path / "TRIES.md"
    monkeypatch.setattr(h36, "TRIES", f)
    for text in ("| H34 (x) | heart |\n", ledger(row, cells=10), ledger(row, first="H36 (another wording)")):
        f.write_text(text)
        for step in ("entry", "profile", "gate", "run", "read"):
            with pytest.raises(SystemExit, match="H36 is refused"):
                h36.refuse(step, COMMIT, tmp_path / "out", tmp_path / "before")
        with pytest.raises(SystemExit, match="H36 is refused"):
            h36.refuse_test(COMMIT, tmp_path)


def test_each_step_wants_the_step_before_it_at_the_run_commit(h36, monkeypatch, tmp_path):
    monkeypatch.setattr(h36, "refuse_state", lambda commit, what: None)
    for step, name, verdict in (("profile", "entry", "ENTER"), ("gate", "profile", "PASS"), ("run", "gate", "PASS"),
                                ("read", "gate", "PASS")):
        before = tmp_path / f"before_{step}"
        before.mkdir()
        with pytest.raises(SystemExit, match=f"no {name} record"):
            h36.refuse(step, COMMIT, None, before)
        (before / f"{name}.json").write_text(json.dumps(rec_of(h36, name, "STOP")))
        with pytest.raises(SystemExit, match="said 'STOP'"):
            h36.refuse(step, COMMIT, None, before)
        (before / f"{name}.json").write_text(json.dumps(rec_of(h36, name, verdict, commit="e" * 40)))
        with pytest.raises(SystemExit, match="not of the run commit"):
            h36.refuse(step, COMMIT, None, before)
        (before / f"{name}.json").write_text(json.dumps(rec_of(h36, name, verdict)))
        h36.refuse(step, COMMIT, None, before)
    out = tmp_path / "out"
    out.mkdir()
    (out / "h36.json").write_text("{}")
    with pytest.raises(SystemExit, match="no step twice"):
        h36.refuse("run", COMMIT, out, tmp_path / "before_run")


# ------------------------------------------------------------------------------------------------ E, Q and T
def raw_of(h36, profile=None, seed=0, e95=None):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_xpc_table import synth

    gs, _ = synth(h36.xpc, seed=seed, profile=profile, e95=e95)
    raw = {h36.xpc.MERFISH: {}, h36.xpc.RNA: {}}
    for (s, inst, k), C in gs.items():
        libs = None if inst == h36.xpc.MERFISH else np.array(
            ["3"] * C.shape[0] if s == "E10.5" else ["1", "2"] * (C.shape[0] // 2) + ["1"] * (C.shape[0] % 2))
        raw[inst].setdefault(s, {})[k] = (C, libs)
    return raw


def test_e_q_and_t_end_to_end(h36, row, monkeypatch, tmp_path, capsys):
    from test_xpc_table import GENES

    monkeypatch.setattr(h36, "refuse_state", lambda commit, what: None)
    monkeypatch.setattr(h36, "git_state", lambda: {"commit": COMMIT, "tree_clean": True})
    monkeypatch.setattr(h36.xpc, "panel_genes", lambda: list(GENES))
    seen = []
    p = np.random.default_rng(5).normal(0, 0.5, len(GENES))
    p9 = np.random.default_rng(6).normal(0, 0.3, len(GENES))  # E9.5's own run profile: the practice profile's
    full = raw_of(h36, profile=p, e95=p9)

    def read(data_dir, stages, genes):
        seen.append(tuple(stages))
        return {inst: {s: full[inst][s] for s in stages} for inst in full}

    out = tmp_path / "res"
    assert h36.run_entry(out / "entry", tmp_path, COMMIT, read=read) == 0
    assert json.loads((out / "entry" / "entry.json").read_text())["verdict"] == "ENTER"
    assert h36.run_profile(out / "profile", tmp_path, COMMIT, out / "entry", read=read) == 0
    assert seen == [("E8.5", "E9.5", "E10.5"), ("E8.5", "E9.5")]  # Q reads no E10.5 file
    prof = json.loads((out / "profile" / "profile.json").read_text())
    assert prof["verdict"] == "PASS" and prof["failed"] == []
    # the practice record read back as the run step reads it (by gene name, in the transition's order)
    with pytest.raises(ValueError):
        h36.practice_delta(out / "profile", list(GENES)[::-1][:-1])
    d, e, level, sha = h36.practice_delta(out / "profile", list(GENES)[::-1])
    rec = json.loads((out / "profile" / "practice_profile.json").read_text())
    np.testing.assert_array_equal(d, np.asarray(rec["delta"])[::-1])
    assert sha == prof["practice_profile_sha256"] and level.shape == e.shape == (len(GENES),)
    assert h36.xpc.spearman(np.asarray(rec["delta"])[np.asarray(rec["eligible"])],
                            p9[np.asarray(rec["eligible"])]) > 0.9
    # T: refused until the read's PASSED AT PRACTICE is committed and in the row's outcome cell
    reading_dir = tmp_path / "res"
    (reading_dir / "reading.json").write_text(json.dumps({"row": "H36", "void": False, "outcome": h36.PASSED,
                                                          "line": row.PASSED_AT_PRACTICE}))
    f = tmp_path / "TRIES.md"
    monkeypatch.setattr(h36, "TRIES", f)
    monkeypatch.setattr(h36, "HERE", tmp_path)
    tracked = {"yes": True}
    monkeypatch.setattr(h36, "git", lambda *a: "res/reading.json" if tracked["yes"] else "")
    prof_path = tmp_path / "xpc_profile.json"
    f.write_text(ledger(row, "pending"))
    with pytest.raises(SystemExit, match="outcome cell does not hold"):
        h36.run_test(tmp_path, COMMIT, reading_dir, profile_path=prof_path, read=read)
    f.write_text(ledger(row, f"\"{row.PASSED_AT_PRACTICE}\""))
    tracked["yes"] = False
    with pytest.raises(SystemExit, match="not committed"):
        h36.run_test(tmp_path, COMMIT, reading_dir, profile_path=prof_path, read=read)
    tracked["yes"] = True
    capsys.readouterr()
    assert h36.run_test(tmp_path, COMMIT, reading_dir, profile_path=prof_path, read=read) == 0
    line = capsys.readouterr().out.strip()
    text = prof_path.read_text()
    assert line == row.pass_statement(h36.xpc.sha256_text(text))
    t = json.loads((reading_dir / "test_profile" / "test_profile.json").read_text())
    assert t["verdict"] == "PASS" and t["line"] == line
    assert h36.xpc.spearman(np.asarray(json.loads(text)["delta"]), p - p9 / 2) > 0.9  # E10.5's less half E9.5's
    np.testing.assert_array_equal(h36.xpc.load_profile(prof_path, h36.xpc.sha256_text(text), list(GENES)),
                                  np.asarray(json.loads(text)["delta"]))
    with pytest.raises(SystemExit, match="exists"):
        h36.run_test(tmp_path, COMMIT, reading_dir, profile_path=prof_path, read=read)
    f.write_text(ledger(row, f"{row.PASSED_AT_PRACTICE} {line}"))
    assert row.passed_sha256(f.read_text())[0] == h36.xpc.sha256_text(text)


def test_t_stops_on_a_profile_that_fails_and_writes_no_file(h36, row, monkeypatch, tmp_path, capsys):
    from test_xpc_table import GENES

    monkeypatch.setattr(h36, "refuse_state", lambda commit, what: None)
    monkeypatch.setattr(h36.xpc, "panel_genes", lambda: list(GENES))
    full = raw_of(h36, profile=None, seed=7)

    def read(data_dir, stages, genes):
        return {inst: {s: full[inst][s] for s in stages} for inst in full}

    (tmp_path / "reading.json").write_text(json.dumps({"row": "H36", "void": False, "outcome": h36.PASSED,
                                                       "line": row.PASSED_AT_PRACTICE}))
    f = tmp_path / "TRIES.md"
    f.write_text(ledger(row, row.PASSED_AT_PRACTICE))
    monkeypatch.setattr(h36, "TRIES", f)
    monkeypatch.setattr(h36, "HERE", tmp_path)
    monkeypatch.setattr(h36, "git", lambda *a: "reading.json")
    prof_path = tmp_path / "xpc_profile.json"
    capsys.readouterr()
    assert h36.run_test(tmp_path, COMMIT, tmp_path, profile_path=prof_path, read=read) == 3
    line = capsys.readouterr().out.strip()
    assert line.startswith("Try H36: T STOPPED (") and line.endswith("): OXP125 is not built.")
    assert not prof_path.exists() and (tmp_path / "test_profile" / "xpc_profile_not_passed.json").exists()
    f.write_text(ledger(row, f"{row.PASSED_AT_PRACTICE} {line}"))
    assert row.passed_sha256(f.read_text())[0] is None
