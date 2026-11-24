#!/usr/bin/env python
"""Try K29 (experiments/harness/TRIES.md; every step refuses without row K29, on an unclean tracked tree or at a HEAD
other than the commit given): the practice read that T3's file cD7 on C2 is bound to (ruling row R36). A labelled
analogue that can only strike, never select; not the adoption rule. It asks, on the one released knockout at the test's
stage and reference, with the barred factor: what does deepening the composition cut gain when the states are that
knockout's own pre-registered list, and what does it cost when they are another list (the scenario's)?

Built on k18_built_files.py, imported unchanged: row K18's four cells (plain g1, plain g2, dm g1, dm g2), its stream
(`t3_s2_form`), seeds, barred spec and factor, plain build, depth-matched row, stamp, task-file writer, folder rule,
masking and guard distance. For a run, three of k18's module names are replaced inside this process only (`k18_arms`):
ARMS (this try's arms), build_arm (this try's) and fold_factors (K18's, plus this try's compositions), restored on exit,
as rows K21 to K24 do it. Also imported unchanged: k21_factor_forms.py (the seed parser, K18's committed rows at their
pins, _read_json, _gated), k22_outside_forms.py (the frozen-file check, the stamp check) and k7_wt_revert.py (_paired).

The arms. Each is K18's barred_form call (bs.construct with the barred spec's keywords through prediction_inputs: the
Mab21l2 factor at exponent 1, no depth or stage step), every arm on stable_rng(seed, t3_s2_form), with only the spec's
composition block changed; the builder's composition routine (build_submission.composition, t3_decompose.rule_comp)
cuts the listed states by 2^-m before renormalising:
  own_h, own_f, own_d   rules.yaml presets.gata4 cm_states (nine labels; seven exist, V-CM and IFT-CM are absent from
                        the vocabulary, which the gata4 run's meta lists) at m = 0.358906935493353, 0.717813870986706
                        and 1.3826521560990648
  oth_h, oth_d          the scenario's five depleted states (rules.yaml presets.gata4 scenario) at the first and third m
  oth_f                 the scenario at the second m: K18's barred_form value for value (its composition is checked
                        equal to the barred spec's, array for array). THE GUARD, not counted while inside it
The three m are fixed here and checked against the committed decision file (the second is the barred s1's, the first
half of it, the third and its depletion factor 0.3835131214774091 the allowed s1's); the two lists are read from the
frozen rules.yaml. Every arm's composition is checked on each world's vocabulary before anything is drawn (the own list
missing exactly V-CM and IFT-CM, the scenario missing nothing, every cut state at 2^-m before renormalising); a failed
check raises, so no task file is written.

Steps (each with --commit, the full hash of the run commit, which must be HEAD; outputs in the run folder R):
  gate --out R            THE TARGET-FREE GATE, before the guard: the frozen builder's composition, given the decision
                          file's submissions.allowed.s1 composition on WT E8.75's vocabulary (the builder's own
                          load_reference), finds all five scenario states (a); each state's (f_i / p_i) / (f_j / p_j),
                          for every uncut label j, equals the decision file's depletion factor 0.3835131214774091 within
                          a relative 1e-12 (b); rule_comp with presets.gata4 cm_states finds exactly seven labels, V-CM
                          and IFT-CM missing (c). R/gate.json; it prints counts and ratios only. Exit 5 when a check
                          fails (nothing is scored, no try counts; the reading is VOID naming the check); a crash exits
                          1 (INCOMPLETE)
  run --part guard --world plain --task g1|g2 [--seeds 0-4] --out R   the guard part: oth_f only, into R/guard/;
  run --part guard --world dm [--workers 8 --threads 2] --out R       refused without a passing gate.json of this commit
                          and machine, or with K18's committed files off their pins
  guard --runs R          the guard check (R/guard.json, made once): oth_f within K17's tolerance (ten-seed mean 0.1,
                          every seed 0.5 board points) of K18's committed barred_form in all four cells. Exit 6 when it
                          fails (VOID with one counted try, oth_f; the run step never starts)
  run [--part score] --world plain|dm ... --out R   a scoring part: the six arms; refused without a passing gate.json
                          and guard.json of this commit and machine. THE FIVE TRIES COUNT FROM THE START OF THE RUN
                          STEP: the first scoring part to pass its refusals writes R/run_started.json
  read --runs R [--out DIR]   the fixed reading: DIR/reading.json and DIR/reading.txt (DIR defaults to
                          experiments/t3_decomp/results/k29_cd7_read), never over an existing one; it prints the outcome
                          line, the per-cell table of G, Gd and C and the labels

Layout of R: gate.json, guard/plain_<task>/seed_<s>.json, guard/dm/<task>_seed_<s>.json, guard.json, run_started.json,
plain_<task>/seed_<s>.json, dm/<task>_seed_<s>.json: one task file per seed and cell, each with its <file>.prov.json
(k18.write_task_file); dm/setup.json and guard/dm/setup.json. A rerun at the same commit on the same machine skips the
seeds already written; a folder holding a file of another commit or machine is refused (k18.check_folder). No log holds
a score: progress lines only, K18's timing line per dm row, the fixed self-check line, the gate's counts and ratios; a
crash prints its traceback with every decimal number masked (k18.masked).

The reading (fixed in row K29, cells 8a and 8b; `reading`, `decide`, `labels`):
  INCOMPLETE  gate.json, guard.json or a seed file missing (the scoring part's only after a passing guard): not a
              reading; read opens no seed file, writes nothing and exits 3. A read made at or after 2026-11-25 06:00
              UTC with a file still missing writes instead the VOID reading "INCOMPLETE past its cutoff".
  VOID        checked first, and it wins (no strike clause is applied): the gate failed ("the target-free gate failed
              ([check])", no try counts); the guard failed at the guard step (one counted try, oth_f); a seed file that
              misses an arm (the plain ones the floor too) or is not its part's, task's or seed's; mixed commits or
              machines; a dirty tree; K18's committed task files off their pins; the guard step's oth_f, read again,
              off K18's rows; the run step's oth_f off the guard step's beyond K17's tolerance. Outcome "Try K29: VOID:
              [reasons]; no file." Once the run step has started, the five tries count, and oth_f also when it fails
              its guard in the run step.
  Per cell, paired by seed: G = own_f - own_h; Gd = own_d - own_h; C = oth_d - oth_h. AHEAD: mean above 0.3 with at
              least 8 of 10 seeds positive. BEHIND: mean below -0.3.
  STRIKE      S1 Gd is AHEAD in no plain cell, or in no dm cell. S2 G is BEHIND in two or more cells. S3 C's four-cell
              mean is below -0.3, or C is BEHIND in two or more cells. S4 in any cell own_d or oth_d has more gated
              seeds (r2_gate, beta_gate, no_response) than its _h arm plus 1. Outcome "Try K29: STRUCK ([clauses]): no
              file; C2 stays."
  OTHERWISE   "Try K29: NOT STRUCK: the file may be built and uploaded under the captain's bound yes (ruling row
              R36)." with the per-cell table of G, Gd and C beside it.
  LABELLED, deciding nothing: "a fit of this fold, not an effect" if (L1) C is AHEAD in at least as many of the four
              cells as Gd; (L2) Gd is AHEAD in a dm cell, and in each dm cell where it is AHEAD its mean over the seeds
              without a gate flip is 0.3 or less (a gate flip: a seed where own_d and own_h have different severity
              states; a cell with no flip-free seed counts as 0.3 or less); (L3) counting a dm cell for Gd only where
              its flip-free mean is above 0.3, Gd's AHEAD cells (one at least) lie in one world, or G is AHEAD somewhere
              and all its AHEAD cells lie in one world.
  REPORTED, deciding nothing: each arm against oth_f; own_d - own_f; per-term skills; gated counts; the realised
              cut-state shares. Both files' gate bars belong to the file build before upload, not to this read.
Nothing changes after any score: the three magnitudes, 0.3, 8 of 10, -0.3, two cells, plus 1, the seeds, 0.1, 0.5, the
relative 1e-12 and the cutoff of 2026-11-25 06:00 UTC.

Exit codes: 0 done; 1 a crash (masked traceback; seeds written so far stay; the gate's: INCOMPLETE); 2 refused (no
row K29, an unclean tree, another HEAD, the frozen T3 files off their pins, K18's committed files off their pins, a
decision file or rules whose magnitudes, depletion factor or lists are not the row's, a part without the passing
gate.json or guard.json it needs); 3 INCOMPLETE; 4 a reading, gate.json or guard.json is there already (never
overwritten); 5 the gate failed; 6 the guard failed.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import sys
import tempfile
import time
import traceback
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (_paired, MIN_GAIN, MIN_AHEAD, FULL_PANEL)
import k18_built_files as k18  # noqa: E402  (worlds, tasks, stream, barred spec and factor, builds, scoring; unchanged)
import k21_factor_forms as k21  # noqa: E402  (seed parser, K18's committed rows, _read_json, _gated; unchanged)
import k22_outside_forms as k22  # noqa: E402  (the frozen-file check, the stamp check; unchanged)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.evaluate import SEVERITY_GATED  # noqa: E402

SCRIPT = "experiments/t3_decomp/k29_cd7_read.py"
ROW_SCRIPT = f"`{SCRIPT}`"  # row K29 names this script
ROW_LABEL = "K29"
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"  # looked up at call time (row_present)
STREAM = k18.STREAM  # "t3_s2_form"
TASKS = k18.TASKS  # g1, g2
WORLDS = k18.WORLDS  # plain, dm
CELLS = ("plain g1", "plain g2", "dm g1", "dm g2")
DM_CELLS = ("dm g1", "dm g2")

# The three magnitudes (log2 of the cut; the states go to 2^-m of their wild-type share before renormalising)
L2_H = 0.358906935493353
L2_F = 0.717813870986706
L2_D = 1.3826521560990648
DEPLETION_D = 0.3835131214774091  # the decision file's submissions.allowed.s1 depletion_factor, 2^-L2_D
OWN, OTH = "own", "oth"  # rules.yaml presets.gata4 cm_states; its scenario's depleted_states
OWN_ABSENT = ("V-CM", "IFT-CM")  # the own list's labels the vocabulary lacks (the gata4 run's meta lists them)
OWN_FOUND = 7
ARM_OF = {"own_h": (OWN, L2_H), "own_f": (OWN, L2_F), "own_d": (OWN, L2_D),
          "oth_h": (OTH, L2_H), "oth_f": (OTH, L2_F), "oth_d": (OTH, L2_D)}
GUARD_ARM = "oth_f"
COUNTED = ("own_h", "own_f", "own_d", "oth_h", "oth_d")  # the five counted tries
ARMS = (GUARD_ARM, *COUNTED)
STEPS = {"G": ("own_f", "own_h"), "Gd": ("own_d", "own_h"), "C": ("oth_d", "oth_h")}
GATED_PAIRS = (("own_d", "own_h"), ("oth_d", "oth_h"))  # S4
FLIP_PAIR = ("own_d", "own_h")  # L2, L3: a gate flip compares these two arms' severity states seed by seed
REPORTED_PAIRS = (*((x, GUARD_ARM) for x in COUNTED), ("own_d", "own_f"))
MIN_GAIN, MIN_AHEAD = k7.MIN_GAIN, k7.MIN_AHEAD  # AHEAD: mean above 0.3 with at least 8 of 10 seeds positive
MAX_BEHIND = -k7.MIN_GAIN  # BEHIND: mean below -0.3; S3's four-cell mean below -0.3
MIN_CELLS = 2  # S2, S3: "two or more cells"
MAX_EXTRA_GATED = 1  # S4: "plus 1"
GUARD_MEAN, GUARD_SEED = k18.GUARD_MEAN, k18.GUARD_SEED  # 0.1 and 0.5 (K17's tolerance)
GATE_REL = 1e-12  # the gate: each cut over the uncut labels against the depletion factor, relative
REL_TOL = 1e-9  # the arms' checks: a cut state's share before renormalising against 2^-m, relative
GATE_CHECKS = {"a": "a: the five scenario states found, states_missing empty",
               "b": "b: each state's (f_i / p_i) / (f_j / p_j) over the uncut labels j equals the depletion factor "
                    "within a relative 1e-12",
               "c": "c: rule_comp with presets.gata4 cm_states finds exactly seven labels, V-CM and IFT-CM missing"}
CUTOFF = datetime(2026, 11, 25, 6, 0, tzinfo=UTC)
CUTOFF_TEXT = "2026-11-25 06:00 UTC"
K18_TASKS = k21.K18_TASKS  # experiments/t3_decomp/results/k18_built_files/tasks
K18_SHA256 = {"plain_g1.json": "5ec09ef583862ab47ab078e72c7aed807c89b4d9c234c57e3b0727ff9e460ed9",
              "plain_g2.json": "cfc5ede7e7185ba5e16d4774d975b9e140e60bbd29736fb85c0fdb8e99f4634f",
              "dm.json": "1ef3e26347c65f915c1630130da02a655c155b912392180001cba315fe135b60"}
RESULTS = HERE / "results" / "k29_cd7_read"
GATE_FILE, GUARD_FILE, GUARD_DIR, RUN_MARK = "gate.json", "guard.json", "guard", "run_started.json"
EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_GATE, EXIT_GUARD = 2, 3, 4, 5, 6
VOID_LINE = "Try K29: VOID: {}; no file."
STRUCK_LINE = "Try K29: STRUCK ({}): no file; C2 stays."
NOT_STRUCK_LINE = ("Try K29: NOT STRUCK: the file may be built and uploaded under the captain's bound yes (ruling row "
                   "R36).")
GATE_VOID = "the target-free gate failed ({})"
LATE = "INCOMPLETE past its cutoff"
FIT_LABEL = "a fit of this fold, not an effect"
_K29: dict = {}  # the two state lists, loaded before any arm builds


class Refused(Exception):
    """A refusal of the whole step (exit status 2)."""


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


def _now() -> datetime:
    return datetime.now(UTC)


def past_cutoff() -> bool:
    return _now() >= CUTOFF


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row and the preconditions


def row_present(path: Path | None = None) -> str:
    """Row K29: exactly one line of TRIES.md that opens with '| K29 ' and names this script (exit 2 otherwise). Other
    rows may cite the script; they neither stand in for the row nor block it."""
    path = TRIES_MD if path is None else Path(path)
    try:
        lines = Path(path).read_text().splitlines()
    except OSError:
        lines = []
    rows = [ln for ln in lines if ln.startswith(f"| {ROW_LABEL} ") and ROW_SCRIPT in ln]
    if len(rows) != 1:
        refuse(f"{'no' if not rows else 'more than one'} row {ROW_LABEL} of {Path(path).name} names {ROW_SCRIPT}; the "
               "row is committed before any step")
    return ROW_LABEL


DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing written; INCOMPLETE, not VOID")


def head_ok(st: dict, commit: str) -> None:
    if not commit or st["commit"] != commit:
        refuse(f"HEAD is {st['commit']}, not the run commit {commit}")


def preconditions(commit: str) -> tuple[str, dict]:
    """Before anything loads: row K29, a clean tracked tree, HEAD the commit given, the frozen T3 files at their
    pins. Returns (row id, stamp)."""
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    head_ok(st, commit)
    k22.frozen_ok()
    return rid, st


def step_passed(out: Path, name: str, st: dict, what: str) -> None:
    """A part starts only after a passing gate.json (guard.json) of this commit and machine (exit 2 otherwise)."""
    p = Path(out) / name
    if not p.is_file():
        refuse(f"{p} is missing: {what} runs and passes first")
    try:
        g = json.loads(p.read_text())
    except ValueError:
        refuse(f"{p} is not JSON")
    if not isinstance(g, dict) or g.get("passes") is not True:
        refuse(f"{what} did not pass: no part runs after it")
    prov = g.get("provenance") or {}
    if prov.get("commit") != st["commit"] or prov.get("machine") != st["machine"]:
        refuse(f"{name} is from another commit or machine: refused")


# ---------------------------------------------------------------------------------------------------------------------
# The compositions


def state_lists(rules: dict) -> dict[str, list[str]]:
    """The frozen rules' own list (presets.gata4 cm_states) and other list (its scenario's depleted states)."""
    pre = rules["presets"]["gata4"]
    return {OWN: [str(s) for s in pre["cm_states"]], OTH: [str(s) for s in pre["scenario"]["depleted_states"]]}


def allowed_s1(dec: dict) -> dict:
    return (((dec.get("submissions") or {}).get("allowed") or {}).get("s1") or {}).get("composition") or {}


def check_spec(spec: dict, dec: dict, lists: dict) -> None:
    """The row's magnitudes and lists against the committed files: the barred s1's composition is the scenario's five
    states at the second magnitude, the first is half of it, and the decision file's allowed s1 composition is the
    scenario at the third with its depletion factor 0.3835131214774091. Refused otherwise."""
    comp = spec.get("composition") or {}
    if comp.get("kind") in (None, "none") or list(comp.get("depleted_states") or []) != lists[OTH]:
        raise Refused("the barred s1's composition is not the scenario's five depleted states")
    if comp.get("log2_magnitude") != L2_F or L2_H != L2_F * 0.5:
        raise Refused(f"the barred s1's magnitude {comp.get('log2_magnitude')!r} is not the row's second, {L2_F!r} "
                      "(and the first half of it)")
    allowed = allowed_s1(dec)
    if allowed.get("kind") in (None, "none") or list(allowed.get("depleted_states") or []) != lists[OTH] \
            or allowed.get("log2_magnitude") != L2_D or allowed.get("depletion_factor") != DEPLETION_D:
        raise Refused("the decision file's allowed s1 composition is not the scenario at the row's third magnitude "
                      f"{L2_D!r} with the depletion factor {DEPLETION_D!r}")
    if len(lists[OTH]) != 5 or len(lists[OWN]) != 9 or any(s not in lists[OWN] for s in OWN_ABSENT):
        raise Refused("the frozen rules' lists are not the row's (five scenario states; nine own labels holding V-CM "
                      "and IFT-CM)")


def load_lists(decision: Path) -> tuple[dict, dict, dict]:
    """The barred spec, the decision file and the frozen rules as k18.barred_spec resolves them; the row's checks; the
    two lists held in _K29 for fold_factors. Returns (spec, decision, rules)."""
    spec, dec, rules = k18.barred_spec(Path(decision))
    lists = state_lists(rules)
    check_spec(spec, dec, lists)
    _K29.clear()
    _K29.update(lists=lists)
    return spec, dec, rules


def gate_eval(bs, V, dec: dict, rules: dict) -> dict:
    """The target-free gate, on a vocabulary (WT E8.75's, as the builder loads it), the decision file and the frozen
    rules: (a) the builder's composition with the decision file's submissions.allowed.s1 composition finds all five
    scenario states; (b) for each of them and every uncut label j, (f_i / p_i) / (f_j / p_j) equals the decision
    file's depletion factor within a relative GATE_REL; (c) rule_comp with presets.gata4 cm_states finds exactly seven
    labels, V-CM and IFT-CM missing. Counts and ratios only; no share is kept."""
    comp = allowed_s1(dec)
    lists = state_lists(rules)
    f, info = bs.composition(V, {"composition": comp}, {})
    p = np.asarray(V.p, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    found = [str(s) for s in info.get("states_found") or []]
    missing = [str(s) for s in info.get("states_missing") or []]
    cut = V.index(found)
    uncut = np.setdiff1d(np.arange(p.size), cut)
    D = float(comp.get("depletion_factor"))
    ratios = {}
    for t, i in zip(found, cut, strict=True):
        r = (f[i] / p[i]) / (f[uncut] / p[uncut]) if uncut.size else np.array([np.nan])
        ratios[t] = {"min": float(np.min(r)), "max": float(np.max(r))}
    checks = {"a": not missing and sorted(found) == sorted(lists[OTH]) and len(found) == 5}
    checks["b"] = bool(checks["a"] and uncut.size and D == DEPLETION_D and all(
        abs(v["min"] / D - 1.0) <= GATE_REL and abs(v["max"] / D - 1.0) <= GATE_REL for v in ratios.values()))
    _, own_missing = bs.t3d.rule_comp(V, {"deplete": dict.fromkeys(lists[OWN], 0.0)}, {}, L2_D)
    own_missing = [str(s) for s in own_missing]
    own_found = len(lists[OWN]) - len(own_missing)
    checks["c"] = own_missing == list(OWN_ABSENT) and own_found == OWN_FOUND
    failed = [GATE_CHECKS[k] for k in ("a", "b", "c") if not checks[k]]
    return {"passes": not failed, "failed": failed, "checks": checks, "scenario_found": len(found),
            "scenario_missing": missing, "ratios": ratios, "depletion_factor": D, "uncut_labels": int(uncut.size),
            "own_found": own_found, "own_listed": len(lists[OWN]), "own_missing": own_missing}


def gate_lines(res: dict) -> list[str]:
    """The gate's printed lines: counts and ratios only."""
    head = "passes" if res["passes"] else "FAILS (" + "; ".join(res["failed"]) + ")"
    return [f"K29 target-free gate: {head}",
            f"  scenario states found {res['scenario_found']} of 5, missing {res['scenario_missing'] or 'none'}",
            "  (f_i / p_i) / (f_j / p_j) over the uncut labels j, least and most: " + ", ".join(
                f"{t} {v['min']!r} {v['max']!r}" for t, v in res["ratios"].items())
            + f"; the depletion factor {res['depletion_factor']!r}",
            f"  rule_comp with cm_states: {res['own_found']} of {res['own_listed']} labels found, missing "
            f"{', '.join(res['own_missing']) or 'none'}"]


def arm_spec(spec: dict, arm: str, lists: dict) -> dict:
    """The barred s1's spec with only its composition block changed: the arm's list at the arm's magnitude (the kind
    stays the barred s1's). The builder's composition routine reads the kind, the states and the magnitude."""
    which, l2 = ARM_OF[arm]
    out = copy.deepcopy(spec)
    out["composition"] = {"kind": spec["composition"]["kind"], "depleted_states": list(lists[which]),
                          "log2_magnitude": float(l2), "depletion_factor": float(2.0 ** -l2),
                          "magnitude_source": f"row K29, arm {arm}"}
    return out


def composition_figures(V, f, info: dict, l2: float) -> dict:
    """Target-free: each cut state's wild-type share, its share after renormalising, and its share before renormalising
    relative to its wild-type share (the renormalising divisor recovered from the uncut labels)."""
    p = np.asarray(V.p, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    found = [str(s) for s in info.get("states_found") or []]
    idx = V.index(found)
    rest = np.setdiff1d(np.arange(p.size), idx)
    div = float(p[rest].sum() / f[rest].sum()) if rest.size and f[rest].sum() > 0 else None
    states = {}
    for t, i in zip(found, idx, strict=True):
        after = float(f[i] / p[i]) if p[i] > 0 else None
        states[t] = {"wt_share": float(p[i]), "share": float(f[i]), "ratio_after": after,
                     "ratio_before_renormalising": None if after is None or div is None else after * div}
    return {"log2_magnitude": float(l2), "cut": float(2.0 ** -l2), "states_found": found,
            "states_missing": [str(s) for s in info.get("states_missing") or []], "cut_share_wt": float(p[idx].sum()),
            "cut_share": float(f[idx].sum()), "renormalising_divisor": div, "states": states}


def composition_checks(figs: dict) -> list[str]:
    """Why the arms' compositions are not the row's (empty: they are): the own list missing exactly V-CM and IFT-CM,
    the scenario missing nothing, every cut state at 2^-m of its wild-type share before renormalising."""
    why = []
    for a in ARMS:
        v = figs.get(a) or {}
        which = ARM_OF[a][0]
        if which == OWN and list(v.get("states_missing") or []) != list(OWN_ABSENT):
            why.append(f"{a}: the own list misses {v.get('states_missing')}, not exactly {list(OWN_ABSENT)}")
        if which == OTH and (v.get("states_missing") or len(v.get("states") or {}) != 5):
            why.append(f"{a}: the scenario misses {v.get('states_missing')}")
        cut = 2.0 ** -ARM_OF[a][1]
        if not v.get("states") or any(s.get("ratio_before_renormalising") is None
                                      or abs(s["ratio_before_renormalising"] / cut - 1.0) > REL_TOL
                                      for s in v["states"].values()):
            why.append(f"{a}: a cut state is not at 2^-m of its wild-type share before renormalising")
    return why


def arm_compositions(bs, V, spec: dict, lists: dict) -> dict:
    """Each arm's spec and composition through the builder's routine, their target-free figures and the checks; oth_f's
    composition must equal the barred spec's array for array. A failed check raises (no task file is written)."""
    specs, f, figs = {}, {}, {}
    for a in ARMS:
        specs[a] = arm_spec(spec, a, lists)
        f[a], info = bs.composition(V, specs[a], {})
        figs[a] = composition_figures(V, f[a], info, ARM_OF[a][1])
    why = composition_checks(figs)
    f_barred, _ = bs.composition(V, spec, {})
    if not np.array_equal(np.asarray(f[GUARD_ARM]), np.asarray(f_barred)):
        why.append("oth_f's composition is not the barred s1's value for value")
    if why:
        raise RuntimeError("the arms' compositions fail the row's checks: " + "; ".join(why))
    return {"specs": specs, "f": f, "figures": figs}


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K18's fold_factors (the original, called unchanged) plus the six arms' specs and compositions on this world's
    reference vocabulary."""
    if not _K29.get("lists"):
        raise RuntimeError("the state lists are not loaded (k29_cd7_read.py loads them before any arm)")
    F = _K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    F["k29"] = arm_compositions(bs, ref.V, spec, _K29["lists"])
    return F


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """One arm's drawn cells: K18's barred_form call (prediction_inputs' keywords into bs.construct, on
    stable_rng(seed, STREAM)) with the arm's spec, whose composition block alone differs from the barred s1's."""
    if arm not in ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    rng = bs.t3d.stable_rng(seed, STREAM)
    kw, _ = bs.prediction_inputs(ref, F["k29"]["specs"][arm], {}, F["phi_barred"])
    if kw["ratio"] is not None or kw["alphas"] is not None or kw["w"] != 1.0:
        raise RuntimeError("the barred s1 is the Mab21l2 factor at exponent 1 with no depth ratio or stage mixing")
    return bs.construct(ref, n, rng, **kw)


_K18_ORIGINAL = {"ARMS": k18.ARMS, "build_arm": k18.build_arm, "fold_factors": k18.fold_factors}


@contextmanager
def k18_arms(arms: tuple[str, ...]):
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's, restored on exit and on an
    exception. K18's file is not changed; its plain build and depth-matched row then build this try's arms."""
    if not set(arms) <= set(ARMS) or GUARD_ARM not in arms:
        raise RuntimeError(f"arms {arms} are not a set of this try's arms with the guard arm")
    k18.ARMS, k18.build_arm, k18.fold_factors = tuple(arms), build_arm, fold_factors
    k18._FOLD.clear()
    try:
        yield
    finally:
        k18.ARMS = _K18_ORIGINAL["ARMS"]
        k18.build_arm = _K18_ORIGINAL["build_arm"]
        k18.fold_factors = _K18_ORIGINAL["fold_factors"]
        k18._FOLD.clear()


def part_arms(part: str) -> tuple[str, ...]:
    return (GUARD_ARM,) if part == "guard" else ARMS


def candidates(arms: tuple[str, ...]) -> list:
    def one(arm):
        which, l2 = ARM_OF[arm]
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: k18._build_plain(ctx, arm),
                                 params={"form": arm, "stream": STREAM, "states": which, "log2_magnitude": l2,
                                         "factor": "the barred s1's (build_submission.py --branch barred --submission "
                                                   "s1)", "guard": arm == GUARD_ARM},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in arms)]


# ---------------------------------------------------------------------------------------------------------------------
# The gate


def cmd_gate(a) -> None:
    """The target-free gate before the guard: R/gate.json (made once); exit 0 when it passes, 5 when a check fails
    (nothing is scored; the reading is VOID naming the check), 1 on a crash (INCOMPLETE: no gate.json)."""
    rid, st = preconditions(a.commit)
    out = Path(a.out)
    p = out / GATE_FILE
    if p.exists():
        print(f"refused: {p} exists: the gate is made once", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    try:
        spec, dec, rules = load_lists(Path(a.decision))
        bs = harness._t3_head()
        wt = Path(a.data) / "T2" / "E8.75.h5ad"
        ref = bs.load_reference(wt, "gata4", rules, threads=a.threads)
        res = gate_eval(bs, ref.V, dec, rules)
        out.mkdir(parents=True, exist_ok=True)
        k18.write_task_file(p, {**res, "row": rid, "composition": allowed_s1(dec), "wt_sha256": bs.sha256(wt),
                                "decision_sha256": bs.sha256(Path(a.decision)), "provenance": st})
    except Refused as e:
        refuse(str(e))
    except Exception:  # noqa: BLE001  (a crash of the gate: masked traceback, exit 1, INCOMPLETE)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print("\n".join(gate_lines(res)), flush=True)
    if not res["passes"]:
        sys.exit(EXIT_GATE)


# ---------------------------------------------------------------------------------------------------------------------
# Running


def progress_log(part: str, world: str, task: str):
    """The harness's per-seed line lists every arm's score; this prints the part (guard only), world, task and seed."""
    lead = "K29 guard" if part == "guard" else "K29"

    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"{lead} {world} {task} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


def k18_rows(refs: Path = K18_TASKS) -> tuple[dict, list[str]]:
    """K18's committed barred_form rows per (world, task), from its task files at this script's pins (K21's reader)."""
    saved = dict(k21.K18_SHA256)
    k21.K18_SHA256.clear()
    k21.K18_SHA256.update(K18_SHA256)
    try:
        return k21.k18_rows(refs)
    finally:
        k21.K18_SHA256.clear()
        k21.K18_SHA256.update(saved)


def guard_refs(refs: Path = K18_TASKS) -> dict:
    """What the guard compares with: K18's committed rows, refused when a file is missing or off its pin."""
    rows, why = k18_rows(refs)
    if why:
        raise Refused("; ".join(why) + ": refused")
    return rows


def mark_run_started(out: Path, st: dict) -> None:
    """The run step's start, which counts the five tries: R/run_started.json, made once by the first scoring part
    that passes its refusals (an exclusive create; the other parts find it there)."""
    try:
        with open(Path(out) / RUN_MARK, "x") as fh:
            json.dump({"row": ROW_LABEL, "counted": list(COUNTED), "commit": st["commit"], "machine": st["machine"],
                       "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}, fh)
    except FileExistsError:
        pass


def run_plain(a, st: dict, out: Path, part: str) -> None:
    """One plain task, seed by seed through harness.run_task (K8's harness), written as
    <root>/plain_<task>/seed_<s>.json with the world's compositions; seeds already written at this commit on this
    machine are skipped."""
    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / f"plain_{a.task}"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    arms = part_arms(part)
    lead = "K29 guard" if part == "guard" else "K29"
    with k18_arms(arms):
        for s in k21.parse_seeds(a.seeds):
            path = folder_ / f"seed_{s}.json"
            if path.exists():
                print(f"{lead} plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            res = harness.run_task(TASKS[a.task], candidates(arms), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log(part, "plain", a.task))
            got = res["provenance"]
            if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
                raise k22._Unclean(f"the harness's stamp says commit {got.get('commit')}, clean "
                                   f"{got.get('tree_clean')}: no task file written (INCOMPLETE)")
            k22._stamp_ok(st)
            fold = k18._FOLD.get("F")
            if not fold or "k29" not in fold:
                raise RuntimeError("the fold's compositions were not made: no task file written")
            res["provenance"]["machine"] = st["machine"]
            res |= {"world": "plain", "part": part, "task_source": SCRIPT, "seed": s, "arms": list(arms),
                    "compositions": fold["k29"]["figures"], "barred_factor": rec}
            prov = k18.write_task_file(path, res)
            print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}, machine {prov['machine']}",
                  flush=True)


def run_dm(a, st: dict, out: Path, part: str) -> None:
    """The depth-matched world, set up as k18.run_dm and rows K21 to K24 set it up (the same functions in the same
    order), then each (replicate, seed) row through k18.dm_row with this part's arms; rows already written at this
    commit on this machine are skipped."""
    import multiprocessing as mp

    import orth_diag as od

    root = out / GUARD_DIR if part == "guard" else out
    folder_ = root / "dm"
    folder_.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder_, st)
    seeds = k21.parse_seeds(a.seeds)
    bs = harness._t3_head()
    t3d = bs.t3d
    spec, dec, rules = k18.barred_spec(Path(a.decision))
    if list(rules["protocol"]["seeds"]) != list(SEEDS):
        raise RuntimeError("the frozen script's seeds are not vec.adoption.SEEDS")
    paths = od.context_paths(Path(a.data), with_targets=True)
    n_par = a.workers * a.threads
    by_name = {v: k for k, v in k18.DM_TARGETS.items()}
    arms = part_arms(part)
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=Path(tmp))
        if sorted(T.name for T in ctx.targets) != sorted(k18.DM_TARGETS.values()):
            raise RuntimeError(f"targets {[T.name for T in ctx.targets]}, expected {sorted(k18.DM_TARGETS.values())}")
        k18.quiet_selfcheck(t3d, ctx)
        genes = list(ctx.genes)
        phi_b, frec = k18.barred_factor(Path(a.run_dir), dec, spec, genes)
        ref_C, ref_lab, ref_P, gap = od.reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        wt_path = Path(paths["wt"])
        series = bs.wt_revert_series(wt_path, {n: Path(a.data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                               if n != wt_path.stem})
        phi, per, group, rrec = bs.wt_revert_factors(wt_path, series, genes, ctx.wt.pool.C)
        refw = bs.Reference(preset=ctx.preset_name, rules=rules, genes=genes, V=ctx.V, pool=ctx.wt.pool,
                            coords=np.zeros((ctx.wt.pool.C.shape[0], 3)))
        files = [*(Path(a.data) / f for f in od.WT_FILES.values()), *(Path(a.data) / f for f in od.KO_FILES.values()),
                 *series.values(), Path(a.decision), Path(a.run_dir) / "factors" / f"{od.BARRED}.csv"]
        with k18_arms(arms):
            F = k18.fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b,
                                 spec)
            comps = F["k29"]["figures"]
            setup = {"world": "dm", "part": part, "barred_factor": frec,
                     "wt_revert": {k: v for k, v in rrec.items() if k != "series"},
                     "reference_counts_largest_gap": gap, "compositions": comps,
                     "files_sha256": {str(f): t3d.sha256(f) for f in dict.fromkeys(files)}}
            sp = folder_ / "setup.json"
            if sp.exists():
                old = json.loads(sp.read_text())
                if {k: old.get(k) for k in setup} != json.loads(json.dumps(setup, default=float)):
                    sys.exit(f"{sp} records another set-up (inputs, factor or compositions): not run")
            else:
                k18.write_task_file(sp, setup | {"provenance": st})
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder_ / f"{by_name[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                k22._stamp_ok(st)
                path = folder_ / f"{by_name[row['target']]}_seed_{row['seed']}.json"
                k18.write_task_file(path, row | {"world": "dm", "part": part, "task": by_name[row["target"]],
                                                 "arms": list(arms), "compositions": comps, "provenance": st,
                                                 "task_source": SCRIPT})
                print(f"wrote {path}", flush=True)

            if a.workers > 1 and len(todo) > 1:
                if not bs.openmp_limited():
                    raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
                with mp.get_context("fork").Pool(min(a.workers, len(todo))) as p:
                    for row in p.imap_unordered(k18._dm_task, todo):
                        write(row)
            else:
                for t in todo:
                    write(k18._dm_task(t))


def cmd_run(a) -> None:
    if a.world == "plain" and a.task not in TASKS:
        refuse("--world plain needs --task g1 or g2")
    if a.world != "plain" and a.task is not None:
        refuse("--task is for --world plain only")
    _, st = preconditions(a.commit)
    out = Path(a.out)
    step_passed(out, GATE_FILE, st, "the target-free gate")
    if a.part == "guard":
        try:
            guard_refs()
        except Refused as e:
            refuse(str(e))
    else:
        step_passed(out, GUARD_FILE, st, "the guard part")
    try:
        load_lists(Path(a.decision))
        out.mkdir(parents=True, exist_ok=True)
        if a.part == "score":
            mark_run_started(out, st)
        (run_plain if a.world == "plain" else run_dm)(a, st, out, a.part)
    except Refused as e:
        refuse(str(e))
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seeds written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)


# ---------------------------------------------------------------------------------------------------------------------
# The guard


def guard_files() -> list[str]:
    return ([f"{GUARD_DIR}/plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"{GUARD_DIR}/dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def score_files() -> list[str]:
    return ([f"plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS])


def cell_of(rel: str) -> tuple[str, int]:
    """(cell, seed) of a task file's path relative to its part's root."""
    seed = int(rel.rsplit("_", 1)[1][:-len(".json")])
    if rel.startswith("plain_"):
        return f"plain {rel.split('/')[0][len('plain_'):]}", seed
    if rel.startswith("dm/"):
        return f"dm {rel[3:5]}", seed
    raise ValueError(f"unknown task file {rel}")


def rows_of(runs: Path, files: list[str], prefix: str, part: str, want: tuple[str, ...], why: list[str], sha: dict,
            prov: dict, comps: dict | None = None) -> dict:
    """Rows per cell from task files (plain: rows[0]; dm: the row), each checked for its part, task, seed and arms (the
    plain ones the floor too); a file that fails is a reason in `why` and gives no row. Each file's composition record
    goes into `comps` when given."""
    rows: dict = {c: [] for c in CELLS}
    for f in files:
        r, sha[f] = k21._read_json(Path(runs) / f)
        prov[f] = r.get("provenance") or {}
        if comps is not None:
            comps[f] = r.get("compositions")
        cell, seed = cell_of(f[len(prefix):])
        world, task = cell.split()
        if r.get("part") != part:
            why.append(f"{f}: not a file of the {part} part")
            continue
        if world == "plain":
            if (r.get("task") or {}).get("key") != TASKS[task].key or [q.get("seed") for q in r.get("rows") or []] \
                    != [seed]:
                why.append(f"{f}: another task or seed")
                continue
            row, need = r["rows"][0], (harness.FLOOR, *want)
        else:
            if r.get("task") != task or r.get("seed") != seed:
                why.append(f"{f}: another task or seed")
                continue
            row, need = r, want
        miss = [x for x in need if x not in (row.get("scores") or {}) or x not in (row.get("skills") or {})
                or x not in (row.get("t3") or {})]
        if miss:
            why.append(f"{f}: missing arms {miss}")
            continue
        rows[cell].append(row)
    return {c: sorted(v, key=lambda q: q["seed"]) for c, v in rows.items()}


def guard_eval(rows: dict, refs: dict) -> dict:
    """The guard: oth_f within K17's tolerance (k18.guard_passes: ten-seed mean within 0.1, every seed within 0.5) of
    K18's committed barred_form in all four cells (the ten seeds needed in each)."""
    cells = {}
    for c in CELLS:
        rr = rows.get(c) or []
        ok = [q["seed"] for q in rr] == list(SEEDS)
        d = k18.distance(rr, GUARD_ARM, refs.get(tuple(c.split())) or [], k18.BARRED) if ok else None
        cells[c] = {"distance": d, "passes": bool(ok and k18.guard_passes(d))}
    passes = all(v["passes"] for v in cells.values())
    return {"cells": cells, "failed": [] if passes else [GUARD_ARM], "passes": passes}


def match_eval(rows: dict, grows: dict) -> dict:
    """The run step's oth_f against the guard step's, seed by seed, within the same tolerance in all four cells."""
    cells = {}
    for c in CELLS:
        rr, gg = rows.get(c) or [], grows.get(c) or []
        ok = [q["seed"] for q in rr] == list(SEEDS) and [q["seed"] for q in gg] == list(SEEDS)
        d = k18.distance(rr, GUARD_ARM, gg, GUARD_ARM) if ok else None
        cells[c] = {"distance": d, "passes": bool(ok and k18.guard_passes(d))}
    return {"cells": cells, "passes": all(v["passes"] for v in cells.values())}


def cmd_guard(a) -> None:
    """The guard check over R/guard/: R/guard.json (made once); exit 0 when the guard holds, 6 when it fails (VOID: no
    scoring part starts; oth_f counts one try), 3 when a guard file is missing, 2 without a passing gate.json, when the
    files are not this commit's and machine's or K18's files are off their pins."""
    rid, st = preconditions(a.commit)
    runs = Path(a.runs)
    step_passed(runs, GATE_FILE, st, "the target-free gate")
    p = runs / GUARD_FILE
    if p.exists():
        print(f"refused: {p} exists: the guard check is made once", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    missing = [f for f in guard_files() if not (runs / f).exists()]
    if missing:
        print(f"INCOMPLETE: {len(missing)} guard file(s) missing", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    try:
        refs = guard_refs()
    except Refused as e:
        refuse(str(e))
    try:
        why: list[str] = []
        sha, prov = {}, {}
        rows = rows_of(runs, guard_files(), f"{GUARD_DIR}/", "guard", (GUARD_ARM,), why, sha, prov)
        bad = [f for f, q in prov.items() if q.get("tree_clean") is not True or q.get("commit") != st["commit"]
               or q.get("machine") != st["machine"]]
        if bad or why:
            refuse(f"the guard files are not this commit's and machine's, or not well formed ({(bad + why)[:4]})")
        res = guard_eval(rows, refs)
        k18.write_task_file(p, {**res, "row": rid, "files_sha256": sha, "provenance": st})
    except SystemExit:
        raise
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print(f"guard: {'passes' if res['passes'] else 'FAILS'}", flush=True)
    if not res["passes"]:
        sys.exit(EXIT_GUARD)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def summarise(rows: list[dict], arms: tuple[str, ...]) -> dict:
    """Per arm: the score's mean and sd, each metric's mean skill, the severity states, seeds scored and gated, the
    median beta and r2. Paired (k7._paired): G, Gd and C; each counted arm against oth_f; own_d - own_f; oth_f against
    the floor (plain). Per-term: each pair's mean per-seed skill difference. Gate flips: the seeds where own_d and own_h
    have different severity states, and Gd's mean over the other (flip-free) seeds, None when there is none."""
    out: dict = {"arms": {}, "paired": {}, "per_term": {}}
    for x in arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "severity_states": {v: sum(z["state"] == v for z in sev) for v in sorted({z["state"] for z in sev})},
            "severity_scored": int(sum(v["state"] == "scored" for v in sev)), "severity_gated": k21._gated(sev),
            "median_beta": k18._median([v.get("beta") for v in sev]),
            "median_r2": k18._median([v.get("r2") for v in sev])}
    pairs = dict(STEPS) | {f"{x} - {y}": (x, y) for x, y in REPORTED_PAIRS}
    if harness.FLOOR in arms:
        pairs[f"{GUARD_ARM} - {harness.FLOOR}"] = (GUARD_ARM, harness.FLOOR)
    out["paired"] = {k: k7._paired(rows, x, y) for k, (x, y) in pairs.items()}
    for k, (x, y) in pairs.items():
        if y == harness.FLOOR:
            continue
        out["per_term"][k] = {m: float(np.mean([q["skills"][x][m] - q["skills"][y][m] for q in rows]))
                              for m in rows[0]["skills"][y]}
    x, h = FLIP_PAIR
    same = [q["t3"][x]["severity"]["state"] == q["t3"][h]["severity"]["state"] for q in rows]
    free = [q["scores"][x] - q["scores"][h] for q, s in zip(rows, same, strict=True) if s]
    out["flips"] = {"seeds": [q["seed"] for q, s in zip(rows, same, strict=True) if not s],
                    "flip_free_seeds": len(free), "flip_free_mean_gd": float(np.mean(free)) if free else None}
    out["gated_flips"] = {f"{a} - {b}": [q["seed"] for q in rows if (q["t3"][a]["severity"]["state"] in SEVERITY_GATED)
                                         != (q["t3"][b]["severity"]["state"] in SEVERITY_GATED)]
                          for a, b in GATED_PAIRS}
    return out


def _ahead(v: dict) -> bool:
    return bool(v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD)


def _behind(v: dict) -> bool:
    return bool(v["mean"] < MAX_BEHIND)


def _worlds(cells: list[str]) -> list[str]:
    return sorted({c.split()[0] for c in cells})


def labels(per: dict, ahead: dict) -> dict:
    """The labelled statement of cell 8b, deciding nothing: FIT_LABEL when any of L1 to L3 holds.
    L1  C is AHEAD in at least as many of the four cells as Gd.
    L2  Gd is AHEAD in a dm cell, and in each dm cell where it is AHEAD its mean over the flip-free seeds is 0.3 or
        less (a cell with no flip-free seed counts as 0.3 or less).
    L3  counting a dm cell for Gd only where its flip-free mean is above 0.3, Gd's AHEAD cells (one at least) lie in
        one world; or G is AHEAD somewhere and all its AHEAD cells lie in one world."""
    n_gd, n_c = sum(ahead["Gd"].values()), sum(ahead["C"].values())
    l1 = n_c >= n_gd
    ff = {c: per[c]["flips"]["flip_free_mean_gd"] for c in DM_CELLS}
    up_dm = [c for c in DM_CELLS if ahead["Gd"][c]]
    l2 = bool(up_dm) and all(ff[c] is None or ff[c] <= MIN_GAIN for c in up_dm)
    gd_cells = [c for c in CELLS if ahead["Gd"][c] and (c.startswith("plain") or (ff[c] is not None
                                                                                  and ff[c] > MIN_GAIN))]
    g_cells = [c for c in CELLS if ahead["G"][c]]
    l3 = len(_worlds(gd_cells)) == 1 or (bool(g_cells) and len(_worlds(g_cells)) == 1)
    return {"L1": bool(l1), "L2": bool(l2), "L3": bool(l3), "statement": FIT_LABEL if (l1 or l2 or l3) else None,
            "ahead_cells": {"Gd": n_gd, "C": n_c}, "flip_seeds": {c: per[c]["flips"]["seeds"] for c in DM_CELLS},
            "flip_free_mean_gd": ff, "gd_cells_counted": gd_cells, "g_ahead_cells": g_cells}


def decide(per: dict) -> dict:
    """The strike rule fixed in row K29 (S1 to S4), on summarise()'s output per cell, then the labels L1 to L3."""
    d = {s: {c: per[c]["paired"][s] for c in CELLS} for s in STEPS}
    ahead = {s: {c: _ahead(v) for c, v in d[s].items()} for s in STEPS}
    behind = {s: {c: _behind(v) for c, v in d[s].items()} for s in STEPS}
    plain = [c for c in CELLS if c.startswith("plain")]
    s1 = not any(ahead["Gd"][c] for c in plain) or not any(ahead["Gd"][c] for c in DM_CELLS)
    s2 = sum(behind["G"].values()) >= MIN_CELLS
    c_mean = float(np.mean([d["C"][c]["mean"] for c in CELLS]))
    s3 = c_mean < MAX_BEHIND or sum(behind["C"].values()) >= MIN_CELLS
    gated = {c: {a: per[c]["arms"][a]["severity_gated"] for a in ARMS} for c in CELLS}
    s4_where = [f"{c}: {x} {gated[c][x]} against {h} {gated[c][h]}" for c in CELLS for x, h in GATED_PAIRS
                if gated[c][x] > gated[c][h] + MAX_EXTRA_GATED]
    s4 = bool(s4_where)
    clauses = [k for k, v in (("S1", s1), ("S2", s2), ("S3", s3), ("S4", s4)) if v]
    return {"strike": {"S1": bool(s1), "S2": bool(s2), "S3": bool(s3), "S4": s4}, "clauses": clauses,
            "struck": bool(clauses), "c_four_cell_mean": c_mean, "s4_cells": s4_where,
            "table": {s: {c: {"mean": float(v["mean"]), "sd": float(v["sd"]), "ahead_seeds": int(v["ahead"]),
                              "n": int(v["n"]), "AHEAD": ahead[s][c], "BEHIND": behind[s][c]}
                          for c, v in d[s].items()} for s in STEPS},
            "gated": gated, "label": labels(per, ahead),
            "outcome": STRUCK_LINE.format(", ".join(clauses)) if clauses else NOT_STRUCK_LINE}


INCOMPLETE = ("Try K29: INCOMPLETE: not a reading; gate.json, guard.json or a seed file is missing (a crash, a refusal "
              "such as a dirty tree, out of memory or a kill); nothing was written. If the cause is outside the code, "
              "rerun the missing part into this folder at the run's commit on the same machine (the runs skip what is "
              "there); after a committed code fix, rerun every step at the new commit into a fresh folder and never "
              f"read this one. At or after {CUTOFF_TEXT} a read writes the VOID reading '{LATE}' instead")


def run_started(runs: Path) -> bool:
    """The run step has started: R/run_started.json, or a scoring file, is there."""
    return (Path(runs) / RUN_MARK).exists() or any((Path(runs) / f).exists() for f in score_files())


def _missing(missing: list[str], counted: list[str], started: bool, base: dict) -> dict:
    """Files missing: INCOMPLETE (not a reading) before the cutoff; at or after it, the VOID reading 'INCOMPLETE past
    its cutoff', with the tries already counted."""
    if not past_cutoff():
        return {"incomplete": True, "void": False, "missing": missing, "outcome": INCOMPLETE}
    why = [f"{LATE} ({CUTOFF_TEXT}): {len(missing)} file(s) missing, e.g. {', '.join(missing[:4])}"]
    return base | {"void": True, "late": True, "missing": missing, "why": why, "tries_counted": len(counted),
                   "counted_arms": counted, "run_step_started": started, "outcome": VOID_LINE.format(why[0])}


def gate_report(g: dict) -> dict:
    return {k: g.get(k) for k in ("passes", "failed", "checks", "scenario_found", "scenario_missing", "ratios",
                                  "depletion_factor", "uncut_labels", "own_found", "own_listed", "own_missing")}


def reading(runs: Path, refs: Path = K18_TASKS) -> dict:
    """Apply the reading of row K29 to a run folder: the gate, the guard step, then the run step. VOID is checked first
    and wins; the strike rule is applied only to a run that is not void. No seed file is opened while one the reading
    needs is missing."""
    runs = Path(runs)
    started = run_started(runs)
    sha, prov = {}, {}
    base: dict = {"incomplete": False, "files_sha256": sha, "provenance": prov}
    five = list(COUNTED) if started else []
    if not (runs / GATE_FILE).exists():
        return _missing([GATE_FILE], five, started, base)
    gate, sha[GATE_FILE] = k21._read_json(runs / GATE_FILE)
    prov[GATE_FILE] = gate.get("provenance") or {}
    base["gate"] = gate_report(gate)
    if gate.get("passes") is not True:
        why = [GATE_VOID.format("; ".join(gate.get("failed") or ["unrecorded"]))]
        return base | {"void": True, "why": why, "tries_counted": len(five), "counted_arms": five,
                       "run_step_started": started, "outcome": VOID_LINE.format(why[0])}
    missing = [f for f in (GUARD_FILE, *guard_files()) if not (runs / f).exists()]
    if missing:  # guard.json is opened only once every guard file is there
        return _missing(missing, five, started, base)
    g, sha[GUARD_FILE] = k21._read_json(runs / GUARD_FILE)
    prov[GUARD_FILE] = g.get("provenance") or {}
    base["guard_part"] = {k: g.get(k) for k in ("cells", "failed", "passes")}
    guard_step = g.get("passes") is True
    missing = [f for f in score_files() if not (runs / f).exists()] if guard_step else []
    if missing:
        return _missing(missing, five, started, base)
    why: list[str] = []
    gsha: dict = {}
    grows = rows_of(runs, guard_files(), f"{GUARD_DIR}/", "guard", (GUARD_ARM,), why, gsha, prov)
    if gsha != g.get("files_sha256"):
        why.append("guard.json: its files' sha256 are not the guard files'")
    sha |= gsha
    if not guard_step:  # the run step never starts after a failed guard: oth_f counts one try
        counted = [*five, GUARD_ARM]
        why = ["the guard failed in the guard step: oth_f is off K18's committed barred_form", *why]
        return base | {"void": True, "why": why, "tries_counted": len(counted), "counted_arms": counted,
                       "run_step_started": started, "outcome": VOID_LINE.format("; ".join(why))}
    comps: dict = {}
    rows = rows_of(runs, score_files(), "", "score", ARMS, why, sha, prov, comps)
    for f, p in prov.items():
        if p.get("tree_clean") is not True:
            why.append(f"{f}: a dirty tree")
    commits = {p.get("commit") for p in prov.values()}
    machines = {p.get("machine") for p in prov.values()}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(map(str, commits))}")
    if len(machines) > 1:
        why.append(f"results from different machines {sorted(map(str, machines))}")
    k18r, kwhy = k18_rows(refs)
    why += kwhy
    again = None if kwhy else guard_eval(grows, k18r)
    if again is not None and not again["passes"]:
        why.append("the guard: the guard step's oth_f, read again, is off K18's committed barred_form")
    match = match_eval(rows, grows)
    for c, v in match["cells"].items():
        if not v["passes"]:
            dd = v["distance"]
            why.append(f"the run step's oth_f is off the guard step's in {c}" + (
                f" by {dd['mean_diff']:+.3f} in the mean and {dd['max_abs_diff']:.3f} in one seed" if dd else
                " (a seed is missing)") + f" (limits {GUARD_MEAN} and {GUARD_SEED})")
    out_of_guard = (again is not None and not again["passes"]) or not match["passes"]
    per: dict = {}
    for c, rr in rows.items():
        if [q["seed"] for q in rr] != list(SEEDS):
            why.append(f"{c}: not the ten seeds")
            continue
        per[c] = summarise(rr, (harness.FLOOR, *ARMS) if c.startswith("plain") else ARMS)
    counted = [*COUNTED, *((GUARD_ARM,) if out_of_guard else ())]
    mark = runs / RUN_MARK
    base |= {"guard": again, "guard_match": match, "tries_counted": len(counted), "counted_arms": counted,
             "run_step_started": True, "run_started": json.loads(mark.read_text()) if mark.exists() else None,
             "cells": per, "compositions": {w: next((c for f, c in comps.items() if f.split("/")[0].startswith(w)),
                                                   None) for w in WORLDS}}
    if why:
        return base | {"void": True, "why": why, "outcome": VOID_LINE.format("; ".join(why))}
    dec = decide(per)
    return base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines)),
                   "decision": dec, "outcome": dec["outcome"]}


def reading_code() -> dict:
    st = k18.stamp()
    h = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731
    return {"script": SCRIPT, "script_sha256": h(__file__), "k18_sha256": h(HERE / "k18_built_files.py"),
            "k21_sha256": h(HERE / "k21_factor_forms.py"), "k22_sha256": h(HERE / "k22_outside_forms.py"),
            "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}


def _f(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def _dist(d: dict | None) -> str:
    if d is None:
        return "n/a"
    return (f"mean {d['mean_diff']:+.4f}, largest seed {d['max_abs_diff']:.4f}, seeds differing "
            f"{d['seeds_off'] or 'none'}")


def _comp_lines(comps: dict | None) -> list[str]:
    lines = []
    for a, v in (comps or {}).items():
        lines.append(f"    {a}: cut 2^-{v['log2_magnitude']!r} = {v['cut']:.4f}; cut states' share "
                     f"{v['cut_share_wt']:.4f} (wild type) -> {v['cut_share']:.4f}; missing "
                     f"{v['states_missing'] or 'none'}; "
                     + ", ".join(f"{t} {s['wt_share']:.4f} -> {s['share']:.4f}" for t, s in v["states"].items()))
    return lines


def table_lines(dec: dict) -> list[str]:
    """The per-cell table of G, Gd and C."""
    lines = ["per cell, paired by seed (AHEAD: mean > 0.3 with >= 8 of 10 positive; BEHIND: mean < -0.3):"]
    for s, cells in dec["table"].items():
        a, b = STEPS[s]
        lines.append(f"  {s} = {a} - {b}: " + "; ".join(
            f"{c} {v['mean']:+.2f} (sd {v['sd']:.2f}, {v['ahead_seeds']} of {v['n']} positive)"
            + (" AHEAD" if v["AHEAD"] else "") + (" BEHIND" if v["BEHIND"] else "") for c, v in cells.items()))
    return lines


def label_lines(dec: dict) -> list[str]:
    lab = dec["label"]
    ff = ", ".join(f"{c} {_f(v)}" for c, v in lab["flip_free_mean_gd"].items())
    return [f"labelled, deciding nothing: {lab['statement'] or 'no label'}; L1 {lab['L1']} (AHEAD cells: C "
            f"{lab['ahead_cells']['C']}, Gd {lab['ahead_cells']['Gd']}); L2 {lab['L2']} (gate flips of own_d and own_h "
            f"{lab['flip_seeds']}; Gd's flip-free means {ff}); L3 {lab['L3']} (Gd's AHEAD cells counted "
            f"{lab['gd_cells_counted'] or 'none'}; G's AHEAD cells {lab['g_ahead_cells'] or 'none'})"]


def summary_lines(res: dict) -> list[str]:
    """What `read` prints: the outcome line, then for a decided reading the table of G, Gd and C and the labels."""
    lines = [res["outcome"]]
    if not res.get("incomplete") and not res.get("void"):
        lines += table_lines(res["decision"]) + label_lines(res["decision"])
    return lines


def text(res: dict) -> str:
    lines = [res["outcome"]]
    if res.get("incomplete"):
        miss = res["missing"]
        return "\n".join([*lines, f"missing: {len(miss)} file(s), e.g. " + ", ".join(miss[:6])]) + "\n"
    lines.append(f"tries counted: {res['tries_counted']} ({', '.join(res['counted_arms']) or 'none'}); the run step "
                 f"started: {res['run_step_started']}")
    gt = res.get("gate")
    if gt:
        lines.append(f"target-free gate: passes {gt.get('passes')}; failed {gt.get('failed') or 'none'}; scenario "
                     f"states found {gt.get('scenario_found')} of 5; the own list {gt.get('own_found')} of "
                     f"{gt.get('own_listed')} (missing {gt.get('own_missing')}); depletion factor "
                     f"{gt.get('depletion_factor')!r}; " + ", ".join(
                         f"{t} {v['min']!r} to {v['max']!r}" for t, v in (gt.get("ratios") or {}).items()))
    gp = res.get("guard_part") or {}
    if gp:
        lines.append(f"guard step: passes {gp.get('passes')}, failed {gp.get('failed') or 'none'}")
        for c, v in (gp.get("cells") or {}).items():
            lines.append(f"  guard step {c}: passes {v['passes']}; {_dist(v.get('distance'))}")
    for c, v in ((res.get("guard") or {}).get("cells") or {}).items():
        lines.append(f"  the guard step's oth_f read again, {c}: passes {v['passes']}; {_dist(v['distance'])}")
    for c, v in ((res.get("guard_match") or {}).get("cells") or {}).items():
        lines.append(f"  the run step's oth_f against the guard step's, {c}: passes {v['passes']}; "
                     f"{_dist(v['distance'])}")
    if res["void"]:
        lines.append("VOID (checked first; no strike clause is applied): " + "; ".join(res["why"]))
    else:
        dec = res["decision"]
        lines.append(f"commit {res['commit']}, machine {res['machine']}.")
        lines.append("strike: " + ", ".join(f"{k} {'yes' if v else 'no'}" for k, v in dec["strike"].items())
                     + f"; C's four-cell mean {dec['c_four_cell_mean']:+.2f}; S4 cells {dec['s4_cells'] or 'none'}")
        lines += table_lines(dec) + label_lines(dec)
    lines.append("reported, deciding nothing (each arm against oth_f, own_d - own_f, per-term skills, gated counts, "
                 "the realised cut-state shares; both files' gate bars belong to the file build before upload):")
    for c, v in (res.get("cells") or {}).items():
        lines.append(f"{c}:")
        for a, s in v["arms"].items():
            lines.append(f"    {a:<6} {s['mean']:.2f} (sd {s['sd']:.2f}); severity scored {s['severity_scored']}, "
                         f"gated {s['severity_gated']} {s['severity_states']}, median beta {_f(s['median_beta'])}, "
                         f"r2 {_f(s['median_r2'])}; skills " + ", ".join(f"{m} {x:.3f}" for m, x in
                                                                           s["skills"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
        for k, terms in v["per_term"].items():
            lines.append(f"    per-term skill, {k}: " + ", ".join(f"{m} {x:+.3f}" for m, x in terms.items()))
        lines.append(f"    gate flips of own_d and own_h (seeds): {v['flips']['seeds'] or 'none'}; gated-state "
                     "changes: " + "; ".join(f"{k} {s or 'none'}" for k, s in v["gated_flips"].items()))
    for w, comps in (res.get("compositions") or {}).items():
        lines.append(f"realised cut-state shares, {w} world (target-free):")
        lines += _comp_lines(comps)
    lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    lines.append("g1 and g2 are two replicates of one knockout run; oth_f is the guard (K18's barred_form).")
    if res.get("reading_code"):
        rc = res["reading_code"]
        lines.append(f"read by {rc['script']} {rc['script_sha256'][:12]} (k18 {rc['k18_sha256'][:12]}, k21 "
                     f"{rc['k21_sha256'][:12]}, k22 {rc['k22_sha256'][:12]}) at commit {rc['commit']}, clean "
                     f"{rc['tree_clean']}, machine {rc['machine']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    rid = row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse(DIRTY)
    head_ok(st, a.commit)
    out = Path(a.out)
    targets = [out / "reading.json", out / "reading.txt"]
    if any(p.exists() for p in targets):
        print(f"refused: {out} holds a reading already; it is never overwritten", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    try:
        res = reading(Path(a.runs))
        if not res.get("incomplete"):
            res["reading_code"] = reading_code()
            res["row"] = rid
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
        shown = summary_lines(res)
    except Exception:  # noqa: BLE001  (a crash after the files were opened: masked traceback, nothing written)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if res.get("incomplete"):
        print("INCOMPLETE: not a reading; nothing written", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    out.mkdir(parents=True, exist_ok=True)
    for p, body in zip(targets, (js, t), strict=True):
        with open(p, "x") as fh:  # exclusive: an existing reading is never overwritten
            fh.write(body)
    print("\n".join(shown), flush=True)
    print(f"wrote {targets[0]} and {targets[1]}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    gt = sub.add_parser("gate")
    gt.add_argument("--out", required=True)
    gt.add_argument("--threads", type=int, default=8)
    r = sub.add_parser("run")
    r.add_argument("--part", choices=("guard", "score"), default="score")
    r.add_argument("--world", choices=WORLDS, required=True)
    r.add_argument("--task", choices=sorted(TASKS))
    r.add_argument("--seeds", default=None)
    r.add_argument("--out", required=True)
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--threads", type=int, default=8)
    r.add_argument("--run-dir", default=str(k18.RUN_DIR))
    for q in (gt, r):
        q.add_argument("--data", default=str(harness.DATA_DIR))
        q.add_argument("--decision", default=str(k18.DECISION))
    gd = sub.add_parser("guard")
    gd.add_argument("--runs", required=True)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", default=str(RESULTS))
    for q in (gt, r, gd, d):
        q.add_argument("--commit", required=True)
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "guard": cmd_guard, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
