#!/usr/bin/env python
"""Row H35 of experiments/harness/TRIES.md: the counted screen of row R55's four trend forms (the heart's trend-form
line, the HEART PORTFOLIO design of the step-change round as its cold challenge fixes it; OKS125 struck). A screen for
a set of files chosen by headline only under row R22's rule: not an adoption read and not evidence of development. It
builds no test file and changes no file, rule or default.

The forms are vec.trend_forms', the heart builder's own code for row R55's files OTK125, OTX125, OTW125 and OTP125:
OTS125's step (row H28's strength-scaled b0_trend step at 1.25) with one change each, the direction of each lineage's
log factor:
  category  the merged Endothelial lineage split by chamber category (Endothelial/endocardium and
            Endothelial/vascular_endothelium at every stage); every other cell keeps its merged type, so CM is untouched
  contrast  each lineage's factor less s-bar, the unweighted per-gene mean of the trend factors fitted on the whole
            reference stage's context; each moved lineage's L2 norm held at the banked factor's
  whole     the reference's point on the whole reference stage; the draw and momentum_shift stay on D's pool
  portion   the slope over the stages before the reference; each moved lineage's L2 norm held at the banked factor's
            (on the 1.75-day backtest a single quarter-day step, E8.25_late to E8.5)
The control is OTS125's practice form: the banked step at 1.25 on row H29b's D_f (counted); the copy D_f is a reference
(not counted). Row H34's script is imported, not copied (and through it rows H28's, H29b's, H27's and H26's): the
transitions, D's pool and D_f's size and half, the scorer's calls and the G1 route. No row H34 figure is read.

  g0    target-free and not a try: on the test inputs as the heart builder loads them for row R55's files (its own
        context, taken at its call to the form step: nothing is built or written), each form's applied log factors
        against OTS125's on the drawn cells. VOID, a code fault, if a form's draw or coordinates are not OTS125's. Per
        form: DROP if its mean |log f| exceeds 3 times OTS125's (wild); IDLE if the changed share of the step, the sum
        over the drawn cells of ||lf_form - lf_OTS125||^2 over the sum of ||lf_OTS125||^2, is under 0.02, or (category)
        the split collapses to the merged lineage at every stage and in the draw (fall-through); a form that cannot be
        formed is DROP; PASS otherwise. It writes each form's word and the name of the condition that set it, and
        nothing else: the figures are sums over the test board's drawn cells, which would give back their number, so
        they stay in memory. Exit 0 when a form passes, 3 when none does (STOP: nothing further runs), 5 VOID.
  gate  first on the practice tasks, counting nothing: the copy D_f built and scored, and G1: its record, draw, means
        and paired statistics equal row H29b's committed D_f exactly (row H34's G1 route); and in every seed's context
        the banked step at strength 1.0 is vec.harness._b0_t2(ctx, trend=True) bit for bit (row H28's V5) and each form
        with its change out is the banked step at 1.25 bit for bit. Exit 0 PASS, 3 FAIL (VOID).
  run   one transition, a seed range, one JSON under --out: the control and the forms G0 passed (counted) and the copy
        D_f once more (the gate's tie, not counted), each through D_f's size and half and the checked builder, scored in
        worlds own and CLEAN with the floor. Refused unless G0 and the gate of this commit passed. Prints only
        "<task> seed <s> (<sec> s) done" and "wrote <path>".
  read  the void conditions first, then each form's word on fold E10.5, form minus control paired by seed: DROP if
        below -0.30 in either world with at least 8 of 10 seeds below 0; IDLE if within +-0.10 in both worlds; KEPT
        otherwise; a form G0 did not pass keeps G0's word. STOP when no form is KEPT. Writes reading.json and
        reading.txt (never overwritten); prints the outcome line and what follows, and nothing else.
  refusals  the refusals of one step and nothing else (the launcher asks before it marks a step as started).
Every step is refused before anything loads unless TRIES.md holds exactly one row H35 of eleven cells whose first cell
is FIRST_CELL and whose outcome cell holds no outcome yet, the tracked tree is clean, HEAD is the commit given and vec
is this checkout's own. That no step is started twice is the launcher's to hold (experiments/harness/run_h35.sh). The
figures go to the captain or fresh analysts, never beside a test headline; the main session acts on the outcome words
alone, and the control's level is never set beside OTS125's headline (J8).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2]  # this checkout
sys.path.insert(0, str(HERE / "src"))  # this checkout's vec first: a shared environment may name another checkout

from vec import harness, t2_build, trend_forms  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.core import counts_of  # noqa: E402
from vec.prep import load_stage  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h34 = _load("h34_heart_organ_trend", HERE / "experiments" / "harness" / "h34_heart_organ_trend.py")
h28, h29b = h34.h28, h34.h29b
h27, h26, h25, payoff = h34.h27, h34.h26, h34.h25, h34.payoff

ROW = "H35"
RULING_ROW = "R55"
# Row H35's first cell, exactly: every step and both launchers (run_h35.sh, experiments/builds/heart_trend_line_h.sh)
# refuse a row H35 whose first cell is not this text
FIRST_CELL = "H35 (the counted screen of row R55's trend forms: ten tries; G0 first, target-free, not a try)"
N_CELLS = 11  # a ledger row's cells; the eleventh is the outcome cell
OUTCOME_CELL = N_CELLS - 1
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h35-heart-trend-forms"
G0_AXIS, GATE_AXIS = AXIS + ":g0", AXIS + ":gate"
PART = "h35.json"
G0_JSON, G0_TXT = "g0.json", "g0.txt"
GATE_JSON, GATE_TXT = "gate.json", "gate.txt"
SOURCE = "experiments/harness/h35_heart_trend_forms.py (row H35)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h35_heart_trend_forms"
BUILDER = HERE / "experiments" / "builds" / "build_test.py"

FORMS = trend_forms.FORMS  # category, contrast, whole, portion: the design's order
FILES = trend_forms.FILES
STRENGTH = 1.25  # OTS125's
TASK_KEYS = h34.TASK_KEYS  # the two heart backtests on E10.5
UNITS = h34.UNITS
FOLD = h34.FOLD
SEEDS = h34.SEEDS  # 0 to 9
OWN, CLEAN, TRUTHS, WORLDS = h34.OWN, h34.CLEAN, h34.TRUTHS, h34.WORLDS
POOL, D_F = h34.POOL, h34.D_F
PRACTICE_BUDGET = h34.PRACTICE_BUDGET
TOL, FACTOR_RTOL = h34.TOL, h34.FACTOR_RTOL  # 1e-9 on a shape term's skill; relative 1e-12 on a budget
DE_TERMS, TERMS, SHAPE_TERMS = h34.DE_TERMS, h34.TERMS, h34.SHAPE_TERMS
TAG = h34.TAG
FLOOR = harness.FLOOR
CONTROL = "control"  # OTS125's practice form, counted
COPY = "copy"  # row H29b's D_f, a reference: G1's subject and the gate's tie
BANKED = trend_forms.BANKED
DROP_BAR, MIN_SEEDS = 0.30, 8  # DROP: a paired mean below -0.30 with at least 8 of 10 seeds below 0, in either world
IDLE_BAR = 0.10  # IDLE: a paired mean within +-0.10 in both worlds
G0_IDLE_SHARE = 0.02  # G0's idle: the changed share of the step under this
G0_WILD_RATIO = 3.0  # G0's wild: the mean |log f| above this times OTS125's
PASS, FAIL = "PASS", "FAIL"
KEPT, DROP, IDLE, STOP, VOID = "KEPT", "DROP", "IDLE", "STOP", "VOID"
WORDS = (KEPT, DROP, IDLE)
GATE_EXIT = {PASS: 0, FAIL: 3}
G0_EXIT = {PASS: 0, STOP: 3, VOID: 5}
G0_ARGV = ("heart", "--entry", "organ_scope", "--trend-strength", "1.25", "--trend-form", "category", "--z-mode", "raw",
           "--captain-ruling", RULING_ROW)  # the builder's command for OTK125: its context is every file's
# The outcome as the read and G0 print it, and as the build launcher reads it from row H35's outcome cell
OUTCOME_RE = re.compile(rf"OUTCOME {ROW}: (STOP: )?"
                        + "; ".join(rf"{f} ({'|'.join(WORDS)})" for f in FORMS) + r"\.")
OUTCOME_MARK = f"OUTCOME {ROW}:"


def pname(a: str, b: str) -> str:
    return f"{a} - {b}"


RULE_PAIRS = tuple((f, CONTROL) for f in FORMS)
DIAGNOSTIC_PAIRS = ((CONTROL, COPY), *((f, COPY) for f in FORMS))


# ------------------------------------------------------------------------------------------------ the row
def ledger_row(text: str) -> list[str] | None:
    """Row H35's cells when the ledger text holds exactly one row H35 of eleven cells; None otherwise."""
    rows = [ln for ln in text.splitlines() if ln.startswith(f"| {ROW} ")]
    if len(rows) != 1:
        return None
    ln = rows[0].rstrip()
    if not ln.endswith("|"):
        return None
    cells = [c.strip() for c in ln[1:-1].split("|")]
    return cells if len(cells) == N_CELLS else None


def row_formed(text: str) -> bool:
    """Exactly one row H35, eleven cells, the first cell exactly FIRST_CELL."""
    cells = ledger_row(text)
    return cells is not None and cells[0] == FIRST_CELL


def screen_outcome(text: str) -> dict | None:
    """Row H35's outcome as its outcome cell holds it, read by code: None unless the row is formed (row_formed) and its
    outcome cell holds exactly one outcome statement, either VOID or the words of the four forms in the read's form
    (OUTCOME_RE), with STOP exactly when no form is KEPT. Returns {"outcome": VOID | STOP | KEPT, "words": form -> word,
    "allowed": the KEPT forms in the design's order}."""
    if not row_formed(text):
        return None
    cell = ledger_row(text)[OUTCOME_CELL]
    if cell.count(OUTCOME_MARK) != 1:
        return None
    if f"{OUTCOME_MARK} {VOID}" in cell:
        return {"outcome": VOID, "words": {}, "allowed": []}
    m = OUTCOME_RE.search(cell)
    if m is None:
        return None
    words = dict(zip(FORMS, m.groups()[1:], strict=True))
    allowed = [f for f in FORMS if words[f] == KEPT]
    if bool(m.group(1)) == bool(allowed):  # STOP with a KEPT form, or no STOP and none KEPT
        return None
    return {"outcome": STOP if m.group(1) else KEPT, "words": words, "allowed": allowed}


def row_open(text: str) -> bool:
    """The row is formed and its outcome cell holds no outcome yet (a screen is read once)."""
    return row_formed(text) and OUTCOME_MARK not in ledger_row(text)[OUTCOME_CELL]


# ------------------------------------------------------------------------------------------------ refusals
def git_state() -> dict:
    """This checkout's HEAD and whether its tracked files are as committed (from the script's own location)."""
    return h29b.run_provenance()


def vec_home() -> Path:
    return Path(harness.__file__).resolve().parents[2]


def refuse_state(commit: str | None, what: str = "run") -> None:
    """The refusals of every step, before anything loads."""
    try:
        text = TRIES.read_text(encoding="utf-8")
    except OSError:
        text = ""
    if ledger_row(text) is None:
        sys.exit(f"{ROW} is refused: experiments/harness/TRIES.md must hold exactly one row {ROW} of {N_CELLS} cells "
                 f"(a row comes before its {what})")
    if not row_formed(text):
        sys.exit(f"{ROW} is refused: row {ROW}'s first cell is not the text this script fixes: {FIRST_CELL}")
    if not row_open(text):
        sys.exit(f"{ROW} is refused: row {ROW}'s outcome cell already holds an outcome (the screen is read once)")
    if vec_home() != HERE:
        sys.exit(f"{ROW} is refused: vec is imported from {vec_home()}, not from this checkout")
    state = git_state()
    if state.get("tree_clean") is not True:
        sys.exit(f"{ROW} is refused: tracked files differ from the commit")
    if not commit or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit {commit!r}")


def _load_json(path: Path) -> dict | None:
    try:
        rec = json.loads(Path(path).read_text())
        return rec if isinstance(rec, dict) else None
    except (OSError, ValueError):
        return None


def _of_commit(rec: dict, commit: str | None) -> bool:
    prov = rec.get("provenance") or {}
    return bool(commit and rec.get("commit") == commit and prov.get("commit") == commit
                and prov.get("tree_clean") is True)


def g0_record_reasons(g0: dict | None, commit: str | None) -> list[str]:
    """v3: no G0 record of the run commit, or one that passed no form."""
    if not isinstance(g0, dict) or g0.get("axis") != G0_AXIS:
        return ["v3: no G0 record"]
    if not _of_commit(g0, commit):
        return ["v3: the G0 record is not of the run commit"]
    if g0.get("verdict") != PASS or not passed_forms(g0):
        return [f"v3: G0 said {g0.get('verdict')!r} and passed no form"]
    return []


def passed_forms(g0: dict) -> list[str]:
    """The forms G0 passed, in the design's order."""
    forms = g0.get("forms") if isinstance(g0, dict) else None
    return [f for f in FORMS if isinstance(forms, dict) and (forms.get(f) or {}).get("verdict") == PASS]


def gate_record_reasons(gate: dict | None, commit: str | None) -> list[str]:
    """v3: the gate failed, or its record is missing or of another commit than the run commit."""
    if not isinstance(gate, dict) or gate.get("axis") != GATE_AXIS:
        return ["v3: no gate record"]
    if gate.get("verdict") != PASS:
        return [f"v3: the gate said {gate.get('verdict')!r}, not {PASS}"]
    if not _of_commit(gate, commit):
        return ["v3: the gate record is not of the run commit"]
    return []


def refuse_g0(out: Path, commit: str | None) -> None:
    refuse_state(commit, "G0")
    if (Path(out) / G0_JSON).exists() or (Path(out) / G0_TXT).exists():
        sys.exit(f"{ROW} is refused: {out} already holds a G0 record")


def refuse_gate(out: Path, commit: str | None, g0_dir: Path) -> None:
    refuse_state(commit, "gate")
    why = g0_record_reasons(_load_json(Path(g0_dir) / G0_JSON), commit)
    if why:
        sys.exit(f"{ROW} is refused: G0 of {commit!r} has passed no form ({why[0]})")
    if (Path(out) / GATE_JSON).exists() or (Path(out) / GATE_TXT).exists():
        sys.exit(f"{ROW} is refused: {out} already holds a gate record")


def refuse_run(path: Path | None, commit: str | None, g0_dir: Path, gate_dir: Path) -> None:
    refuse_state(commit)
    why = g0_record_reasons(_load_json(Path(g0_dir) / G0_JSON), commit)
    if why:
        sys.exit(f"{ROW} is refused: G0 of {commit!r} has passed no form ({why[0]})")
    if gate_record_reasons(_load_json(Path(gate_dir) / GATE_JSON), commit):
        sys.exit(f"{ROW} is refused: the gate of {commit!r} has not said {PASS} ({Path(gate_dir) / GATE_JSON})")
    if path is not None and Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


def refuse_read(out: Path, commit: str | None) -> None:
    refuse_state(commit, "read")
    if (Path(out) / "reading.json").exists() or (Path(out) / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")


# ------------------------------------------------------------------------------------------------ the fixed words
SCREEN = "A screen, not evidence of development."
NOTHING = (f"What follows: nothing is built under row {RULING_ROW}; the line ends there, with no second read of these "
           "forms under this or another row.")
BUILD = (f"What follows: each KEPT form's file is built twice by experiments/builds/heart_trend_line_h.sh under ruling "
         f"row {RULING_ROW}, before any score of its own; the figures go to the captain or to fresh analysts, never "
         "beside a test headline, and a file goes in only on the captain's word at the time.")


def outcome_line(outcome: str, words: dict | None = None, why=()) -> str:
    """The outcome statement, wording fixed: VOID with every condition found; else each form's word in the design's
    order, after "STOP: " when no form is KEPT (OUTCOME_RE reads it back)."""
    if outcome == VOID:
        return f"{OUTCOME_MARK} {VOID}: {'; '.join(why)}; nothing is built."
    body = "; ".join(f"{f} {words[f]}" for f in FORMS) + "."
    return f"{OUTCOME_MARK} {'STOP: ' if outcome == STOP else ''}{body}"


def follows_line(outcome: str) -> str:
    return f"{SCREEN} {BUILD if outcome == KEPT else NOTHING}"


def verdict_of(words: dict) -> str:
    return KEPT if any(w == KEPT for w in words.values()) else STOP


# ------------------------------------------------------------------------------------------------ G0 (test inputs)
class _Taken(Exception):
    """Raised in place of the builder's form step once its arguments are taken: nothing is built or written."""


def builder():
    """The heart builder (experiments/builds/build_test.py), loaded as a module of this checkout."""
    return _load("build_test_for_h35", BUILDER)


def builder_context(data_dir: Path | None = None, bt=None) -> tuple:
    """The heart builder's own context for row R55's files and the whole E10.5, taken at its call to the form step
    (build_test._trend_form_step, stopped there): the builder's refusals, its stages, its organ cut, its view, its
    test board and its seed. Nothing is built or written, and the builder module is as it was afterwards. `bt` is the
    builder module (default: loaded afresh from this checkout)."""
    bt = bt or builder()
    if data_dir is not None:
        bt.DATA = Path(data_dir)
    seen: list = []
    own = bt._trend_form_step

    def take(ctx, strength, form, whole):
        seen.append((ctx, float(strength), whole))
        raise _Taken

    bt._trend_form_step = take
    try:
        bt.build_heart(bt.build_parser().parse_args([*G0_ARGV, "--out", str(RESULTS / "g0" / "never_written.h5ad")]))
    except _Taken:
        pass
    finally:
        bt._trend_form_step = own
    if len(seen) != 1:
        raise RuntimeError("the heart builder did not reach its form step exactly once")
    ctx, strength, whole = seen[0]
    if strength != STRENGTH:
        raise RuntimeError(f"the builder's strength is {strength}, not {STRENGTH}")
    return ctx, whole


def _same(a, b) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return bool(np.array_equal(np.asarray(a), np.asarray(b)))


def count_free_trend(rec: dict) -> dict:
    """b0_trend_factors' record without a cell number (the builder's _trend_record's fields)."""
    keys = ("lineage", "stages", "days", "own_factor", "weights", "cosine", "norm_ratio", "mean_abs_logf_trend",
            "mean_abs_logf_b0")
    return {"stages": list(rec["stages"]), "days": [float(d) for d in rec["days"]],
            "budget_days": float(rec["budget_days"]),
            "lineages": [{k: r[k] for k in keys} for r in rec["lineages"]],
            "median_cosine": float(rec["median_cosine"]), "median_norm_ratio": float(rec["median_norm_ratio"])}


FORM_FIELDS = {"category": ("split_types", "collapsed"), "contrast": ("shared_lineages", "shared_norm", "lineages"),
               "whole": ("lineages",), "portion": ("lineages", "banked_factor_kept")}
# G0 writes each form's word and the name of the condition that set it, and nothing else (the review of aeb070e): its
# figures are sums over the test board's drawn cells divided by their number, from per-cell factors released data
# reproduce, so a written figure (rounded or not) gives back the hidden test count. They stay in memory. Condition ->
# the word it sets, and the fixed words g0.txt gives it (no numeral, nothing computed).
G0_CONDITIONS = {"none": PASS, "wild": DROP, "not formed": DROP, "fall-through": IDLE, "idle": IDLE,
                 "code fault": VOID}
G0_CONDITION_WORDS = {
    "none": "no failing condition",
    "wild": "its mean absolute log factor over the drawn cells is above the bar's multiple of OTS125's",
    "not formed": "the form cannot be formed on the test inputs",
    "fall-through": "its split collapses to the merged lineage at every stage and in the draw",
    "idle": "its changed share of the step over the drawn cells is under the bar",
    "code fault": "its draw or coordinates are not OTS125's",
}


def _g0_word(condition: str) -> dict:
    return {"verdict": G0_CONDITIONS[condition], "condition": condition}


def g0_assess(ctx, whole, form: str, banked, banked_applied: dict) -> tuple[dict, dict]:
    """One form against the banked step: (what G0 writes, the figures behind it). What G0 writes is the form's word
    and the name of the condition that set it; the figures (the changed share, the two mean absolute log factors and
    their ratio, the fall-through flag, each type's status) are computed over the test draw and are returned for use
    in memory only. VOID on a code fault (another draw or other coordinates), DROP if the form cannot be formed or is
    wild, IDLE if it falls through (category) or is idle, PASS otherwise."""
    try:
        args = trend_forms.form_arguments(ctx, STRENGTH, form, whole)
    except (ValueError, t2_build.UnmappedT2LabelError):
        return _g0_word("not formed"), {}
    if not (_same(args.rows, banked.rows) and _same(args.coords, banked.coords)):
        return _g0_word("code fault"), {}
    app = trend_forms.applied(ctx, args)
    st = trend_forms.g0_statistics(app, banked_applied)
    fell = form == "category" and all(args.record["collapsed"].values())
    figures = {**st, "fall_through": bool(fell), "statuses": app["statuses"]}
    if st["wild_ratio"] > G0_WILD_RATIO:
        return _g0_word("wild"), figures
    if fell:
        return _g0_word("fall-through"), figures
    if st["changed_share"] < G0_IDLE_SHARE:
        return _g0_word("idle"), figures
    return _g0_word("none"), figures


def g0_form(ctx, whole, form: str, banked, banked_applied: dict) -> dict:
    """What G0 writes for one form: its word and the name of the condition that set it (g0_assess)."""
    return g0_assess(ctx, whole, form, banked, banked_applied)[0]


def g0_reading(ctx, whole) -> dict:
    """Every form's word and condition, and the verdict: VOID on any code fault, PASS when a form passes, STOP
    otherwise; on STOP or VOID the outcome line and what follows. Nothing in it is computed from the test draw."""
    banked = trend_forms.form_arguments(ctx, STRENGTH, None)
    banked_applied = trend_forms.applied(ctx, banked)
    forms = {f: g0_form(ctx, whole, f, banked, banked_applied) for f in FORMS}
    faults = [f"G0: {f}: {G0_CONDITION_WORDS['code fault']} (a code fault)" for f in FORMS
              if forms[f]["verdict"] == VOID]
    verdict = VOID if faults else PASS if any(forms[f]["verdict"] == PASS for f in FORMS) else STOP
    words = {f: forms[f]["verdict"] for f in FORMS}
    line = (outcome_line(VOID, why=faults) if verdict == VOID
            else outcome_line(STOP, words) if verdict == STOP else None)
    return {"verdict": verdict, "forms": forms, "why": faults, "line": line,
            "follows": None if line is None else follows_line(STOP)}


G0_TEXT_HEAD = ("Row {row}'s G0 (commit {commit}): target-free and not a try. On the test inputs as the heart builder "
                "loads them for row {ruling}'s files, each form's applied log factors against OTS125's on the drawn "
                "cells. G0 writes each form's word and the condition that set it, and nothing else: its figures are "
                "sums over the test board's drawn cells and would give back their number, so none is written.")


def g0_text(g: dict) -> str:
    lines = [G0_TEXT_HEAD.format(row=ROW, commit=g["commit"], ruling=RULING_ROW)]
    for f in FORMS:
        r = g["forms"][f]
        lines.append(f"G0 {f} ({FILES[f]}): {r['verdict']} ({r['condition']}: {G0_CONDITION_WORDS[r['condition']]}).")
    lines.append(f"G0 {g['verdict']}.")
    if g.get("line"):
        lines += [g["line"], g["follows"]]
    return "\n".join(lines) + "\n"


def run_g0(out: Path, data_dir: Path | None, commit: str | None, bt=None) -> int:
    """The G0 step: the refusals first, then the builder's own context on the test inputs, every form's word and
    condition; the record holds those and the run's provenance, and no figure."""
    out = Path(out)
    refuse_g0(out, commit)
    ctx, whole = builder_context(data_dir, bt)
    res = g0_reading(ctx, whole)
    state = git_state()
    if state.get("tree_clean") is not True or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: tracked files or HEAD changed while G0 ran; no G0 record is written")
    g = {"row": ROW, "axis": G0_AXIS, "commit": commit, "provenance": state,
         "label": "G0 of row H35: target-free, not a try; the test inputs as the builder loads them, no target, no "
                  "score; each form's word and the condition that set it, and no figure", **res, "source": SOURCE}
    out.mkdir(parents=True, exist_ok=True)
    (out / G0_JSON).write_text(json.dumps(g, indent=1, default=float))
    (out / G0_TXT).write_text(g0_text(g))
    for f in FORMS:
        print(f"G0 {f} {res['forms'][f]['verdict']}", flush=True)
    print(f"G0 {res['verdict']}", flush=True)
    if res["line"]:
        print(res["line"], flush=True)
        print(res["follows"], flush=True)
    return G0_EXIT[res["verdict"]]


# ------------------------------------------------------------------------------------------------ the transitions
@dataclass
class Transition:
    """Row H34's transition (D's pool the context's reference input, the earlier stages whole) and the whole
    reference stage, which the contrast and whole forms read."""

    base: object  # h34.Transition
    whole: object  # the reference stage as released, loaded as the builder loads its inputs

    def __getattr__(self, name):
        return getattr(self.base, name)


def load_transition(key: str, data_dir: Path) -> Transition:
    """Row H34's transition, and the whole reference stage loaded as row H34 loads it before its organ cut."""
    tr = h34.load_transition(key, data_dir)
    t = payoff.heart_task(key)
    genes = load_boards()[t.panel_board].genes
    whole = load_stage(Path(data_dir) / t.inputs[tr.view.reference_role], genes=genes, drop_duplicates=True)
    pool = tr.inputs[tr.view.reference_role]
    if not (pool.n_obs < whole.n_obs and set(pool.obs_names) <= set(whole.obs_names)):
        raise ValueError(f"{key}: D's pool is not a part of the whole reference stage")
    return Transition(base=tr, whole=whole)


def context(tr: Transition, s: int) -> harness.BuildContext:
    return h34.context(tr.base, s)


# ------------------------------------------------------------------------------------------------ the gate
def identity(ctx, whole) -> dict:
    """In one seed's context, exactly: the banked step at strength 1.0 is harness._b0_t2(ctx, trend=True) (row H28's
    V5), and each form with its change out is the banked step at 1.25."""
    one = trend_forms.form_step(ctx, 1.0, None)[0]
    b0 = harness._b0_t2(ctx, trend=True)
    banked = trend_forms.form_step(ctx, STRENGTH, None)[0]
    zeroed = {}
    for f in FORMS:
        st = trend_forms.form_step(ctx, STRENGTH, f, whole, zero=True)[0]
        zeroed[f] = bool(_same(st.X, banked.X) and _same(st.coords, banked.coords))
    s1 = bool(_same(one.X, b0.X) and _same(one.coords, b0.coords))
    return {"s1_equals_harness_b0_trend": s1, "zeroed_equals_banked": zeroed, "ok": s1 and all(zeroed.values())}


def g1_reasons(rows: dict, facts: dict, committed: dict | None, h27_rows: dict | None) -> tuple[list[str], dict]:
    """G1, row H34's route with the copy D_f in its control's place: the reading is the file row H29b records, and
    (R), (D), (M) and (P) each equal their entries there exactly."""
    found: dict = {"R": None, "D": None, "M": None, "P": None}
    why = []
    complete = all(key in rows and key in facts and [r["seed"] for r in rows[key]] == list(SEEDS)
                   for key in TASK_KEYS)
    if not complete:
        why.append("G1: the copy D_f was not built and scored on both transitions and every seed")
    if not h34.committed_ok(committed):
        why.append("G1: row H29b's reading.json is not at hand as the file its row records (missing, another sha256, "
                   "or void)")
    elif complete:
        found["R"] = {key: h34.record_differs(key, facts[key], committed) for key in TASK_KEYS}
        why += [f"G1 (R): {key}: {', '.join(off)} not row H29b's committed" for key, off in found["R"].items() if off]
        found["D"] = {key: h34.pool_draw(rows[key]) == committed["pools_short"][key][POOL] for key in TASK_KEYS}
        why += [f"G1 (D): {key}: the copy's draw is not row H29b's committed pools_short.{POOL}"
                for key, same in found["D"].items() if not same]
        mine = h34.means(rows, COPY)
        off = [[w, u, mine[w][u] - committed["tables"][w][u]["mean"][D_F]] for w, _ in h34.G1_WORLDS
               for u in h34.G1_UNITS if mine[w][u] != committed["tables"][w][u]["mean"][D_F]]
        n = len(h34.G1_WORLDS) * len(h34.G1_UNITS)
        found["M"] = {"means": n, "unequal": off, "largest_abs_diff": max((abs(d) for _, _, d in off), default=0.0)}
        if off:
            why.append(f"G1 (M): {len(off)} of {n} means are not row H29b's committed D_f means (first: {off[0][0]} "
                       f"{off[0][1]})")
        if h27_rows is None:
            why.append("G1 (P): row H27's committed parts are not at hand")
        else:
            off = []
            for other in h34.G1_ARMS:
                got = h34.paired(rows, COPY, h27_rows, other)
                for w, _ in h34.G1_WORLDS:
                    for u in h34.G1_UNITS:
                        want = committed["tables"][w][u]["stated"][h29b.pname(D_F, other)]
                        off += [[w, u, h29b.pname(D_F, other), k] for k in h34.STAT_FIELDS if got[w][u][k] != want[k]]
            n = len(h34.G1_ARMS) * len(h34.G1_WORLDS) * len(h34.G1_UNITS) * len(h34.STAT_FIELDS)
            found["P"] = {"numbers": n, "unequal": off}
            if off:
                why.append(f"G1 (P): {len(off)} of {n} paired statistics are not row H29b's committed (first: "
                           f"{' '.join(off[0])})")
    return why, found


def gate_reasons(rows: dict, facts: dict, committed: dict | None, h27_rows: dict | None) -> tuple[list[str], dict]:
    """G1 (g1_reasons) and the identities of every seed's context (identity)."""
    why, found = g1_reasons(rows, facts, committed, h27_rows)
    seen = [(key, r) for key in TASK_KEYS for r in rows.get(key, [])]
    bad = [(key, r["seed"]) for key, r in seen if (r.get("identity") or {}).get("ok") is not True]
    found["I"] = {"contexts": len(seen), "failed": len(bad)}
    if bad:
        why.append(f"I: {len(bad)} of {len(seen)} contexts fail an identity (first: {bad[0][0]} seed {bad[0][1]})")
    return why, found


GATE_LINES = {
    PASS: ("GATE PASS: G1 and the identities hold (the copy D_f reproduces row H29b's committed D_f bit for bit; the "
           "banked step at strength 1 is the harness's b0_trend and each form with its change out is the banked step, "
           "bit for bit); the run may follow at this commit."),
    FAIL: ("GATE FAIL: G1 or an identity does not hold; the run step refuses, no try is counted and the outcome is "
           "VOID."),
}


def gate_text(g: dict) -> str:
    lines = [f"Row {ROW}'s gate (commit {g['commit']}): first, counting nothing. It builds and scores the copy D_f and "
             "nothing else; no counted arm is built or scored in it. G1 and the identities by exact equality, with "
             "no tolerance.",
             f"G1: the copy D_f against row H29b's committed D_f, its reading.json of sha256 "
             f"{h34.H29B_READING_SHA256}: (R), (D), (M) and (P) as row H34's gate reads them.",
             "Identities, in every seed's context: the banked step at strength 1 against vec.harness._b0_t2(ctx, "
             "trend=True); each form with its change out against the banked step at 1.25."]
    f = g.get("found") or {}
    if f.get("M") is not None:
        lines.append(f"(M) found: {f['M']['means']} means compared, {len(f['M']['unequal'])} unequal.")
    if f.get("P") is not None:
        lines.append(f"(P) found: {f['P']['numbers']} paired statistics compared, {len(f['P']['unequal'])} unequal.")
    if f.get("I") is not None:
        lines.append(f"Identities found: {f['I']['contexts']} contexts, {f['I']['failed']} failed.")
    lines.append(GATE_LINES[g["verdict"]])
    if g["verdict"] == FAIL:
        why = [f"v3: the gate fails ({'; '.join(g['why'])})"]
        lines += ["WHY: " + "; ".join(g["why"]), outcome_line(VOID, why=why), follows_line(VOID)]
    return "\n".join(lines) + "\n"


def gate_seed(tr: Transition, s: int, data_dir: Path) -> dict:
    """One seed of the gate: the copy D_f built and scored; the identities checked, nothing else built or scored."""
    t0 = time.time()
    ctx = context(tr, s)
    with harness.guarded_reads(harness.REPO_ROOT, data_dir):
        stage, rec, draw = h34.build_control(ctx, tr.target)
        ident = identity(ctx, tr.whole)
    scores, skills, groups = h34.score_arms(tr.base, s, {COPY: stage})
    return {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": {COPY: int(stage.n)},
            "size_record": rec, "draw": draw, "identity": ident, "scores": scores, "skills": skills,
            "groups": groups, "seconds": round(time.time() - t0, 1)}


def run_gate(out: Path, data_dir: Path, commit: str | None, g0_dir: Path, h29b_reading: Path | None = None,
             h27_parts: Path | None = None) -> int:
    out = Path(out)
    refuse_gate(out, commit, g0_dir)
    committed, h27_rows = h34.load_h29b_reading(h29b_reading), h34.load_h27_rows(h27_parts)
    rows, facts, why = {}, {}, []
    for key in TASK_KEYS:
        try:
            tr = load_transition(key, data_dir)
        except t2_build.UnmappedT2LabelError as e:
            why.append(f"G1: {key}: the copy cannot be built: labels without a T2 chamber category on the reference "
                       f"stage: {sorted(e.labels)}")
            continue
        facts[key], rows[key] = tr.facts, []
        for s in SEEDS:
            rows[key].append(gate_seed(tr, s, data_dir))
            print(f"{key} seed {s} ({rows[key][-1]['seconds']} s) done", flush=True)
        del tr
    state = git_state()
    if state.get("tree_clean") is not True or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: tracked files or HEAD changed while the gate ran; no gate record is written")
    found: dict = {}
    try:
        more, found = gate_reasons(rows, facts, committed, h27_rows)
        why += [w for w in more if w not in why]
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:  # a malformed committed record
        why.append(f"G1: a committed record cannot be read ({type(e).__name__}: {e})")
    g = {"row": ROW, "axis": GATE_AXIS, "commit": commit, "provenance": state,
         "label": "the gate of row H35: counting nothing; the copy D_f alone is built and scored, and its figures are "
                  "row H29b's own", "verdict": FAIL if why else PASS, "why": why, "found": found,
         "h29b_reading_sha256": h34.H29B_READING_SHA256, "transitions": facts, "rows": rows, "source": SOURCE}
    out.mkdir(parents=True, exist_ok=True)
    (out / GATE_JSON).write_text(json.dumps(g, indent=1, default=float))
    (out / GATE_TXT).write_text(gate_text(g))
    print(f"GATE {g['verdict']}", flush=True)
    if why:
        print(outcome_line(VOID, why=[f"v3: the gate fails ({'; '.join(why)})"]), flush=True)
        print(follows_line(VOID), flush=True)
    return GATE_EXIT[g["verdict"]]


# ------------------------------------------------------------------------------------------------ the arms
def arm_record(rec: dict | None, ctx) -> dict:
    """An arm's step record for the read's checks: its form and strength, the budget before the strength, D's pool as
    the reference, the shifted types and each type's status, the form's own fields."""
    if rec is None:
        return {"form": None, "strength": None}
    form = rec["form"]
    out = {"form": form, "strength": float(rec["strength"]), "budget_days": float(rec["budget_days"]),
           "effective_budget_days": float(rec["strength"]) * float(rec["budget_days"]),
           "ref_cells": int(ctx.inputs[ctx.task.reference_role].n_obs), "shifted": list(rec["shifted"]),
           "types": list(rec["types"]), "trend": count_free_trend(rec["trend"])}
    for k in FORM_FIELDS.get(form, ()):
        out[k] = rec[k]
    return out


def build_arms(ctx, whole, target: dict, forms=FORMS) -> tuple[dict, dict]:
    """Every arm of one seed on D's context: the copy D_f (harness._copy), the control (the banked step at 1.25) and
    each form of `forms` (vec.trend_forms.form_step), each then through D_f's size and half and the checked builder
    (h34.size_and_half); and the identity checks the read's v5 holds them to."""
    pre = {COPY: (harness._copy(ctx), None, None, None), CONTROL: trend_forms.form_step(ctx, STRENGTH, None)}
    for f in forms:
        pre[f] = trend_forms.form_step(ctx, STRENGTH, f, whole)
    pool, rows, _, _ = harness.draw_copy(ctx)
    drawn = counts_of(pool[rows], list(ctx.board.genes))
    copy = pre[COPY][0]
    checks = {"same_draw": bool(all(_same(v[2], rows) for k, v in pre.items() if k != COPY)),
              "pre_size_coords_equal": bool(all(_same(v[0].coords, copy.coords) for v in pre.values())),
              "changed_counts": {k: int((v[1] != drawn).nnz) for k, v in pre.items() if k != COPY}}
    out, recs = {}, {}
    for name, (st, _, _, rec) in pre.items():
        out[name], srec = h34.size_and_half(ctx, st, target, name)
        recs[name] = {"size": srec, "trend": arm_record(rec, ctx)}
    checks["post_size_coords_equal"] = bool(all(_same(v.coords, out[COPY].coords) for v in out.values()))
    checks["n_cells"] = {k: int(v.n) for k, v in out.items()}
    return out, {"records": recs, "checks": checks, "draw": h34.draw_record(ctx, copy)}


def run_transition(key: str, seeds, out: Path, data_dir: Path, commit: str | None, g0_dir: Path,
                   gate_dir: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse_run(path, commit, g0_dir, gate_dir)
    if key not in TASK_KEYS:
        sys.exit(f"{key}: not a transition of row {ROW}")
    forms = passed_forms(_load_json(Path(g0_dir) / G0_JSON))
    tr = load_transition(key, data_dir)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        ctx = context(tr, s)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            preds, extra = build_arms(ctx, tr.whole, tr.target, forms)
        scores, skills, groups = h34.score_arms(tr.base, s, preds)
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": extra["checks"]["n_cells"],
               "checks": extra["checks"], "records": extra["records"], "draw": extra["draw"], "scores": scores,
               "skills": skills, "groups": groups, "seconds": round(time.time() - t0, 1)}
        rows_out.append(row)
        print(f"{key} seed {s} ({row['seconds']} s) done", flush=True)
    res = {"label": "row H35: the control (OTS125's practice form) and the trend forms G0 passed on row H29b's D_f "
                    "form; a screen on one truth, not an adoption result", "axis": AXIS, "task": tr.key,
           "forms": forms, "arms": [COPY, CONTROL, *forms], "truth_values": list(TRUTHS), **tr.facts,
           "strength": STRENGTH, "truth_seed": payoff.TRUTH_SEED, "half_seed_base": h25.HALF_SEED_BASE,
           "departures": h27.DEPARTURES, "seeds": list(seeds), "rows": rows_out, "source": SOURCE,
           "provenance": git_state()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


def load_parts(root: Path) -> tuple[dict[str, list[dict]], list[str]]:
    """Every part under root, by transition, and the files named as parts that cannot be read as one."""
    out: dict[str, list[dict]] = {}
    unreadable = []
    for p in sorted(Path(root).rglob(PART)):
        try:
            r = json.loads(p.read_text())
            if not (isinstance(r, dict) and r.get("axis") == AXIS and isinstance(r.get("task"), str)
                    and isinstance(r.get("rows"), list)):
                raise ValueError("not a part of this row")
        except (OSError, ValueError):
            unreadable.append(str(p))
            continue
        out.setdefault(r["task"], []).append(r)
    return out, unreadable


# ------------------------------------------------------------------------------------------------ the void
def _finite(v) -> bool:
    return isinstance(v, int | float) and not isinstance(v, bool) and math.isfinite(v)


def identities(key: str, rows: list[dict], built: tuple) -> list[str]:
    """v5: every arm's draw, coordinates (before and after the size step and the half) and cell number are the copy's,
    and no arm's shape term is more than 1e-9 from the control's in any seed and truth."""
    why = []
    for r in rows:
        c = r["checks"]
        tag = f"v5: {key} seed {r['seed']}"
        for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
            if c.get(k) is not True:
                why.append(f"{tag}: {k.replace('_', ' ')} fails")
        if set(r["n_cells"]) != set(built) or len(set(r["n_cells"].values())) != 1:
            why.append(f"{tag}: the arms' cell numbers are not the copy's")
        for q in TRUTHS:
            for name in built:
                if name == CONTROL:
                    continue
                for m in SHAPE_TERMS:
                    if not abs(r["skills"][q][name][m] - r["skills"][q][CONTROL][m]) <= TOL:
                        why.append(f"{tag} {q}: {name}'s {m} is not the control's")
    return why


def step_records(key: str, rows: list[dict], forms: list[str], pool_cells: int) -> list[str]:
    """v5, the arms' records: each counted arm holds its own form (banked for the control) at strength 1.25, its
    transition's budget before the strength (relative 1e-12) and D's pool as the step's reference."""
    why = []
    want = PRACTICE_BUDGET[key]
    for r in rows:
        tag = f"v5: {key} seed {r['seed']}"
        for name, form in ((CONTROL, BANKED), *((f, f) for f in forms)):
            rec = r["records"][name]["trend"]
            if rec.get("form") != form or rec.get("strength") != STRENGTH:
                why.append(f"{tag}: {name}'s record does not hold form {form} at strength {STRENGTH:g}")
            if not abs(rec.get("budget_days", float("nan")) - want) <= FACTOR_RTOL * want:
                why.append(f"{tag}: {name}'s budget is {rec.get('budget_days')}, not {want}")
            if rec.get("ref_cells") != pool_cells:
                why.append(f"{tag}: {name}'s reference is not D's pool")
    return why


def gate_scores(gate: dict) -> dict[str, dict[int, dict]]:
    return {key: {r["seed"]: r["scores"] for r in gate["rows"][key]} for key in TASK_KEYS}


def gate_tie(rows: dict[str, list[dict]], gate: dict) -> list[str]:
    """v5: the copy and the floor as the run step scores them are exactly the gate's, in every transition, seed and
    truth."""
    try:
        theirs = gate_scores(gate)
        return [f"v5: {key} seed {r['seed']} {q}: the run step's {name} score is not the gate's"
                for key in TASK_KEYS for r in rows[key] for q in TRUTHS for name in (COPY, FLOOR)
                if r["scores"][q][name] != theirs[key][r["seed"]][q][name]]
    except (KeyError, TypeError, IndexError, AttributeError) as e:
        return [f"v3: the gate record holds no score to hold the run step's to ({type(e).__name__}: {e})"]


def void_reasons(parts: dict[str, list[dict]], g0: dict | None, gate: dict | None, commit: str | None,
                 unreadable=()) -> list[str]:
    """Why the read is void (empty when it may be read), by the conditions code can hold.
      v3  G0 or the gate failed, or a record is missing or of another commit than the run commit
      v4  a transition, seed, truth or arm is missing (the arms: the copy, the control and the forms G0 passed), a seed
          occurs twice, a score is not a finite number, a part cannot be read or is of another transition, or the
          parts come from more than one commit, another commit than the run commit or an unclean tree
      v5  an arm is not the arm the row defines (identities, step_records), or the copy or the floor as the run step
          scores them is not exactly the gate's (gate_tie)"""
    v3 = g0_record_reasons(g0, commit) + gate_record_reasons(gate, commit)
    why = v3 + [f"v4: a part cannot be read: {p}" for p in unreadable]
    if g0_record_reasons(g0, commit):
        return why
    forms = passed_forms(g0)
    built = (COPY, CONTROL, *forms)
    every = [p for ps in parts.values() for p in ps]
    complete = True
    for key in TASK_KEYS:
        ps = parts.get(key)
        if not ps:
            why.append(f"v4: missing transition {key}")
            complete = False
            continue
        rows, twice = h26.merge_rows(ps)
        if twice:
            why.append(f"v4: {key}: seeds {twice} occur in more than one part")
            complete = False
        if [r["seed"] for r in rows] != list(SEEDS):
            why.append(f"v4: {key}: seeds {[r['seed'] for r in rows]}, not {list(SEEDS)}")
            complete = False
        for r in rows:
            for field in ("scores", "skills", "groups"):
                gone = sorted({f"{q}/{a}" for q in TRUTHS for a in (FLOOR, *built)
                               if a not in r.get(field, {}).get(q, {})})
                if gone:
                    why.append(f"v4: {key} seed {r['seed']}: missing from {field}: {gone[:6]}")
                    complete = False
            if set(r.get("records", {})) != set(built) or "checks" not in r or "n_cells" not in r:
                why.append(f"v4: {key} seed {r['seed']}: records, checks or cell numbers missing, or other arms")
                complete = False
            bad = sorted({f"{q}/{a}" for q in TRUTHS for a in (FLOOR, *built)
                          if a in r.get("scores", {}).get(q, {}) and not _finite(r["scores"][q][a])})
            if bad:
                why.append(f"v4: {key} seed {r['seed']}: a score is not a finite number: {bad[:6]}")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"v4: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not every:
        return why
    if len({p["provenance"]["commit"] for p in every}) > 1:
        why.append("v4: the parts come from more than one commit")
    if not all(p["provenance"]["tree_clean"] is True for p in every):
        why.append("v4: a part comes from an unclean tree")
    if not commit or any(str(p["provenance"]["commit"]) != commit for p in every):
        why.append(f"v4: a part comes from a commit other than the run commit {commit!r}")
    if any(p.get("forms") != forms for p in every):
        why.append("v4: a part's forms are not the ones G0 passed")
    if not complete:
        return why
    merged = {}
    for key in TASK_KEYS:
        rows = merged[key] = h26.merge_rows(parts[key])[0]
        why += identities(key, rows, built) + step_records(key, rows, forms,
                                                           parts[key][0]["pred_info"][POOL]["pool_cells"])
    return why + ([] if v3 else gate_tie(merged, gate))


# ------------------------------------------------------------------------------------------------ the tables and words
def rule_rows(parts: dict[str, list[dict]], gate: dict) -> dict[str, list[dict]]:
    """The rows the tables read: the run step's, with the copy's scores the gate's (v5 holds them equal)."""
    theirs = gate_scores(gate)
    return {key: [{**r, "scores": {q: {**r["scores"][q], COPY: theirs[key][r["seed"]][q][COPY]}
                                   for q in r["scores"]}} for r in h26.merge_rows(parts[key])[0]]
            for key in TASK_KEYS}


def compare(rows, keys, truths, a: str, b: str) -> dict:
    return h34.compare(rows, keys, truths, a, b)


def pairs_of(forms) -> tuple:
    return (*((f, CONTROL) for f in forms), (CONTROL, COPY), *((f, COPY) for f in forms),
            *((a, FLOOR) for a in (COPY, CONTROL, *forms)))


def tables(rows: dict[str, list[dict]], forms) -> dict:
    return {w: {u: {pname(a, b): compare(rows, ks, ts, a, b) for a, b in pairs_of(forms)} for u, ks in UNITS.items()}
            for w, ts in WORLDS}


def form_word(tab: dict, form: str) -> dict:
    """One read form's word on fold E10.5, form minus control: DROP if below -DROP_BAR in either world with at least
    MIN_SEEDS of the 10 seeds below 0; IDLE if within +-IDLE_BAR in both worlds; KEPT otherwise. A difference of
    exactly 0 is neither above nor below 0; a mean of exactly -0.30 is not below -0.30; a mean of exactly +-0.10 is
    within +-0.10."""
    st = {w: tab[w][FOLD][pname(form, CONTROL)]["score"] for w, _ in WORLDS}
    below = {w: bool(s["mean"] < -DROP_BAR and s["behind"] >= MIN_SEEDS) for w, s in st.items()}
    within = {w: bool(abs(s["mean"]) <= IDLE_BAR) for w, s in st.items()}
    word = DROP if any(below.values()) else IDLE if all(within.values()) else KEPT
    return {"word": word, "below": below, "within": within}


def reading(parts: dict[str, list[dict]], g0: dict | None, gate: dict | None, commit: str | None,
            unreadable=()) -> dict:
    try:
        why = void_reasons(parts, g0, gate, commit, unreadable)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:  # a malformed part
        why = [f"v4: a part cannot be read ({type(e).__name__}: {e})"]
    res = {"row": ROW, "commit": commit, "void": bool(why), "why": why, "screen": SCREEN,
           "g0": {k: (g0 or {}).get(k) for k in ("verdict", "commit", "why")},
           "gate": {k: (gate or {}).get(k) for k in ("verdict", "commit", "why", "found")}}
    if why:
        res.update(outcome=VOID, words={}, line=outcome_line(VOID, why=why))
    else:
        forms = passed_forms(g0)
        rows = rule_rows(parts, gate)
        tab = tables(rows, forms)
        read = {f: form_word(tab, f) for f in forms}
        words = {f: read[f]["word"] if f in read else g0["forms"][f]["verdict"] for f in FORMS}
        outcome = verdict_of(words)
        res.update(outcome=outcome, words=words, read=read, forms=forms, tables=tab,
                   g0_words={f: g0["forms"][f]["verdict"] for f in FORMS if f not in forms},
                   budgets={"test": h28.budget_at_test(STRENGTH),
                            **{k: STRENGTH * b for k, b in PRACTICE_BUDGET.items()}},
                   records={k: rows[k][0]["records"] for k in TASK_KEYS},
                   changed_counts={k: rows[k][0]["checks"].get("changed_counts") for k in TASK_KEYS},
                   d_targets={k: parts[k][0]["d_targets"] for k in TASK_KEYS},
                   pred_info={k: parts[k][0]["pred_info"] for k in TASK_KEYS},
                   truth_info={k: parts[k][0]["truth_info"] for k in TASK_KEYS},
                   line=outcome_line(outcome, words))
    res["follows"] = follows_line(res["outcome"])
    return res


RULE = (f"RULE, fixed in the row before its run, per form on fold E10.5, form minus control paired by seed (the "
        f"control is OTS125's practice form, the banked step at {STRENGTH:g} on row H29b's D_f): DROP if the paired "
        f"mean is below -{DROP_BAR:.2f} in own or in CLEAN with at least {MIN_SEEDS} of the 10 seeds below 0 there; "
        f"IDLE if the paired mean is within +-{IDLE_BAR:.2f} in both worlds; KEPT otherwise. A form G0 did not pass "
        "keeps G0's word (DROP: wild or not formed; IDLE: idle or fallen through). STOP if no form is KEPT. VOID is "
        "decided first, and under VOID nothing is read.")
UNPRICED = ("Unpriced (J1): the test's size at the two-day horizon and its halving, which the practice arms take "
            "before the checked builder with row H27's size analogue through row H29b's ratio; E12.5's sampling, run "
            "profile and endothelium; the test slope's last step across the portion-to-dissected change, which the "
            "backtests do not have; every moved cell's counts drawn afresh in each form (a count draw rides on each "
            "step). OTP125's practice form on the 1.75-day backtest is a single quarter-day step (E8.25_late to E8.5). "
            "J5: one truth (E10.5), so every figure is a labelled diagnostic. J6: the test budget, 0.625 day, is 1.33 "
            "and 4 times the practice budgets (0.469 and 0.156 day): a practice difference is a range from zero, never "
            "a point. No form is ranked by these figures and no upload order comes from them.")


def _st(s: dict) -> str:
    return (f"{s['mean']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10 above, {s['behind']}/10 below "
            f"(se {s['se']:.3f})")


def _terms(c: dict) -> str:
    return ", ".join(f"{m} {c[m]:+.4f}" for m in (*TERMS, *SHAPE_TERMS))


def text(res: dict) -> str:
    lines = [res["line"], res["follows"],
             f"Row {ROW} (commit {res['commit']}): the counted screen of row {RULING_ROW}'s trend forms on row H29b's "
             "D_f form, one truth (E10.5; the two backtests merged seed by seed into fold E10.5), worlds own and "
             "CLEAN, seeds 0 to 9, released stages only. No figure here is to be read beside a test score, and the "
             "control's level is never set beside OTS125's headline. The read builds nothing.",
             f"G0: {res['g0'].get('verdict')} at commit {res['g0'].get('commit')}. GATE: "
             f"{res['gate'].get('verdict')} at commit {res['gate'].get('commit')}."]
    if res["void"]:
        return "\n".join([*lines, RULE, "VOID, the conditions found, one by one:", *res["why"]]) + "\n"
    lines.append(RULE)
    for f in FORMS:
        if f in res["read"]:
            r = res["read"][f]
            lines.append(f"{f} ({FILES[f]}): {r['word']}; " + "; ".join(
                f"{w} {_st(res['tables'][w][FOLD][pname(f, CONTROL)]['score'])}" for w, _ in WORLDS))
        else:
            lines.append(f"{f} ({FILES[f]}): {res['words'][f]} at G0, not read.")
    lines.append(UNPRICED)
    lines.append("TABLES (paired by seed: mean [min, max] seeds above and below 0 of 10 (se); per-term mean skill "
                 "differences): beyond the rule's pairs on fold E10.5, every entry is a practice diagnostic and gates "
                 "nothing:")
    for w, _ in WORLDS:
        for u in UNITS:
            for p, c in res["tables"][w][u].items():
                lines.append(f"{TAG} {w} {u} {p}: {_st(c['score'])}; {_terms(c)}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def _write_new(path: Path, content: str) -> None:
    with open(path, "x", encoding="utf-8") as f:
        f.write(content)


def cmd_g0(a) -> None:
    sys.exit(run_g0(Path(a.out), Path(a.data) if a.data else None, a.commit))


def cmd_gate(a) -> None:
    sys.exit(run_gate(Path(a.out), Path(a.data), a.commit, Path(a.g0), Path(a.h29b), Path(a.h27)))


def cmd_run(a) -> None:
    refuse_run(Path(a.out) / PART, a.commit, Path(a.g0), Path(a.gate))
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS}")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data), a.commit, Path(a.g0),
                   Path(a.gate))


def cmd_read(a) -> None:
    out = Path(a.out)
    refuse_read(out, a.commit)
    parts, unreadable = load_parts(Path(a.runs))
    res = reading(parts, _load_json(Path(a.g0) / G0_JSON), _load_json(Path(a.gate) / GATE_JSON), a.commit,
                  unreadable)
    out.mkdir(parents=True, exist_ok=True)
    _write_new(out / "reading.json", json.dumps(res, indent=1, default=float))
    _write_new(out / "reading.txt", text(res))
    print(res["line"])
    print(res["follows"])


def cmd_refusals(a) -> None:
    """The refusals of one step, made before anything loads, and nothing else: exit 0 when none holds."""
    if a.step == "run":
        if a.g0 is None or a.gate is None:
            sys.exit("refusals --step run needs --g0 and --gate")
        refuse_run(None, a.commit, Path(a.g0), Path(a.gate))
    elif a.out is None:
        sys.exit(f"refusals --step {a.step} needs --out")
    elif a.step == "g0":
        refuse_g0(Path(a.out), a.commit)
    elif a.step == "gate":
        if a.g0 is None:
            sys.exit("refusals --step gate needs --g0")
        refuse_gate(Path(a.out), a.commit, Path(a.g0))
    else:
        refuse_read(Path(a.out), a.commit)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    z = sub.add_parser("g0")
    z.add_argument("--out", required=True)
    z.add_argument("--commit", required=True, help="row H35's run commit (full hash; HEAD must be it)")
    z.add_argument("--data", default=None, help="the data folder (default: the builder's own, data/raw)")
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--g0", required=True, help="G0's folder (its g0.json must have passed a form at --commit)")
    g.add_argument("--commit", required=True)
    g.add_argument("--data", default=str(harness.DATA_DIR))
    g.add_argument("--h29b", default=str(h34.H29B_READING), help="row H29b's committed reading.json")
    g.add_argument("--h27", default=str(h34.H27_PARTS), help="row H27's committed parts")
    r = sub.add_parser("run")
    r.add_argument("--task", required=True, help=f"one of {TASK_KEYS}")
    r.add_argument("--seeds", required=True, help="a range such as 0-1")
    r.add_argument("--out", required=True)
    r.add_argument("--g0", required=True)
    r.add_argument("--gate", required=True)
    r.add_argument("--commit", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--g0", required=True)
    d.add_argument("--gate", required=True)
    d.add_argument("--out", default=str(RESULTS))
    d.add_argument("--commit", required=True)
    f = sub.add_parser("refusals")
    f.add_argument("--step", required=True, choices=("g0", "gate", "run", "read"))
    f.add_argument("--commit", required=True)
    f.add_argument("--out", help="the step's output folder (g0, gate, read)")
    f.add_argument("--g0", help="G0's folder (gate, run)")
    f.add_argument("--gate", help="the gate's folder (run)")
    a = ap.parse_args(argv)
    {"g0": cmd_g0, "gate": cmd_gate, "run": cmd_run, "read": cmd_read, "refusals": cmd_refusals}[a.cmd](a)


if __name__ == "__main__":
    main()
