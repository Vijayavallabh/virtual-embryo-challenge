#!/usr/bin/env python
"""Row H37 of experiments/harness/TRIES.md: the gate of row R69's file WBH, the heart's world-blur hedge
(HEART_MOON2_design.md's P1 with every fix of its cold challenge, HEART_MOON2_challenge.md section 5; the blur is
vec.world_blur). Six counted heart tries: arms W, U and P on the two heart backtests on E10.5. A labelled diagnostic on
one truth; it builds no test file and changes no file, rule or default.

  C   the control, row H35's category arm on row H29b's D_f (OTK125's practice form): vec.trend_forms.form_step at
      strength 1.25 under the category form, then D_f's size and half (h34.size_and_half), as row H35 built it
  W   C plus the blocked blur (32 KMeans blocks on the coordinates, paired by size rank, antithetic worlds), WBH's form
  U   C plus the same blur with a world per cell (random antithetic cell pairs, stream [seed, 65])
  P   W with the spread permuted across genes (stream [seed, 64]): the same sizes, the wrong pattern
The prior reads row R69's five files only (T2 E8.75 and E9.5, T3 Gata4 KO replicates 1 and 2 and the Mab21l2 KO): no
E10.5 cell enters it, so the practice truth stays out.

  g0    target-free and not a try: the prior and its checks (at least 3 pairs with 2 or more groups; median s in
        [0.10, 1.00]; the KO-against-WT pairs' median spread at most twice the replicate pair's, fix 4c); G00 on the
        released E9.5 in the scorer's geometry (the cut (1 - e^(-c/2)) / 2 at least 0.10, fix 4a); and the blur of the
        banked half (OTK125's half, refused unless it is the file of its recorded sha256) at the file's seed with
        G0's output checks (exact count recovery; pseudobulk Spearman >= 0.99 and max gap <= 0.15; dp Spearman >=
        0.99 against E10.5 and against E9.5, fix 4b; the design library within 2%; X changed in at least half the
        cells; antithetic pairs). On the banked half it writes words only (fix 4d). Exit 0 PASS, 3 STOP.
  gate  counting nothing. G1: C rebuilt and scored on both transitions at seeds 0 to 9 equals row H35's committed
        parts exactly (row H35's reading.json of sha256 H35_READING_SHA256 at run commit H35_RUN_COMMIT: the category
        arm's scores, skills, groups and cell number and the floor's scores, in every truth). G0 on the practice
        outputs: on each seed's C, W, U and P are built (never scored) and hold G0's output checks (dp against the
        transition's reference stage and E9.5). Exit 0 PASS, 3 FAIL (VOID), 4 STOP (G0 on the practice outputs).
  run   one transition, a seed range, one JSON under --out: C, W, U and P built, through the checked builder, and
        scored in worlds own and CLEAN with the floor. The six tries count from the run step's start. Refused unless
        G0 and the gate of this commit passed. Prints only "<task> seed <s> (<sec> s) done" and "wrote <path>".
  read  VOID first (VOID_LIST), then on fold E10.5, paired by seed, in worlds own and CLEAN: B1 W - C >= +0.30 with
        at least 8 of 10 seeds above 0, in both worlds; B2 W - P >= +0.10 in both worlds. PASS = B1, labelled "world
        pattern" when B2 holds and "spread only" when it does not. IDLE (|W - C| < 0.10 in both worlds) and NOT
        PASSED are STOP. U - W and the per-term lines are printed and decide nothing. Writes reading.json and
        reading.txt (never overwritten); prints the outcome line and what follows, and nothing else.
  refusals  the refusals of one step and nothing else (the launcher asks before it marks a step as started).
Every step is refused before anything loads unless TRIES.md holds exactly one row H37 of eleven cells whose first cell
is h37_row.FIRST_CELL and whose outcome cell holds no outcome yet, the tracked tree is clean, HEAD is the commit given
and vec is this checkout's own. That no step is started twice is the launcher's to hold
(experiments/harness/run_h37.sh).
The figures go to the captain or fresh analysts, never beside a test headline; the main session acts on the outcome
words alone.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import multiprocessing
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2]  # this checkout
sys.path.insert(0, str(HERE / "src"))  # this checkout's vec first: a shared environment may name another checkout

from vec import harness, trend_forms  # noqa: E402
from vec import world_blur as wb  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.core import counts_of  # noqa: E402
from vec.evaluate import Stage  # noqa: E402
from vec.prep import load_stage  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


hrow = _load("h37_row", HERE / "experiments" / "harness" / "h37_row.py")
h35 = _load("h35_heart_trend_forms", HERE / "experiments" / "harness" / "h35_heart_trend_forms.py")
h34, h29b, h26, payoff = h35.h34, h35.h29b, h35.h26, h35.payoff

ROW, RULING_ROW, FIRST_CELL = hrow.ROW, hrow.RULING_ROW, hrow.FIRST_CELL
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h37-heart-world-blur"
G0_AXIS, GATE_AXIS = AXIS + ":g0", AXIS + ":gate"
PART = "h37.json"
G0_JSON, G0_TXT = "g0.json", "g0.txt"
GATE_JSON, GATE_TXT = "gate.json", "gate.txt"
SOURCE = "experiments/harness/h37_heart_world_blur.py (row H37)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h37_heart_world_blur"
BOARD_KEY = "T2:heart:test_extrap"
# the banked half (OTK125's: build_test.py heart --entry organ_scope --trend-strength 1.25 --trend-form category
# --z-mode raw --captain-ruling R55, then halve.py at 1000001), as experiments/builds/results/
# hearth_8a8bc3a_h2.SHA256SUMS records it
BANKED_NAME = "T2heart_organ_ts125_tk_half.h5ad"
BANKED_SHA256 = "9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54"
FILE_SEED = harness.BUILD_SEED_OFFSET + 1  # WBH's seed, 1000001
# row H35's committed reading and parts, G1's subject
H35_RESULTS = h35.RESULTS
H35_READING_SHA256 = "bf9746388ab5079f6d35350e7a459989c9e3215a2337d34e278e8f0c5ca7637d"
H35_RUN_COMMIT = "53fc2674479c2fe18fec91831523ba8d8dc7aab4"
H35_ARM = "category"  # row H35's name for C

STRENGTH = h35.STRENGTH  # 1.25
FORM = "category"
TASK_KEYS, UNITS, FOLD, SEEDS = h35.TASK_KEYS, h35.UNITS, h35.FOLD, h35.SEEDS
TRUTHS, WORLDS = h35.TRUTHS, h35.WORLDS
TERMS, SHAPE_TERMS, TOL, TAG = h35.TERMS, h35.SHAPE_TERMS, h35.TOL, h35.TAG
FLOOR = harness.FLOOR
CONTROL, W, U, P = "C", "W", "U", "P"
ARMS = (W, U, P)  # counted: six tries, three arms on two transitions
BUILT = (CONTROL, *ARMS)
B1_BAR, B1_SEEDS, B2_BAR, IDLE_BAR = 0.30, 8, 0.10, 0.10
PASS, FAIL, STOP, VOID, IDLE, NOT_PASSED = "PASS", "FAIL", "STOP", "VOID", "IDLE", "NOT PASSED"
WORLD_PATTERN, SPREAD_ONLY = hrow.PASS_LABELS
G0_EXIT = {PASS: 0, STOP: 3}
GATE_EXIT = {PASS: 0, FAIL: 3, STOP: 4}
MARK = hrow.MARK
WORKERS = 10  # the gate's processes (one per transition and seed pair)


def pname(a: str, b: str) -> str:
    return f"{a} - {b}"


def arm_spec(name: str) -> dict:
    """How an arm is blurred: W blocks, U cells, P blocks with the spread permuted."""
    return {W: {"mode": wb.BLOCKS, "permute": False}, U: {"mode": wb.CELLS, "permute": False},
            P: {"mode": wb.BLOCKS, "permute": True}}[name]


# ------------------------------------------------------------------------------------------------ refusals
def git_state() -> dict:
    return h29b.run_provenance()


def vec_home() -> Path:
    return Path(harness.__file__).resolve().parents[2]


def refuse_state(commit: str | None, what: str = "run") -> None:
    """The refusals of every step, before anything loads."""
    try:
        text = TRIES.read_text(encoding="utf-8")
    except OSError:
        text = ""
    if hrow.cells(text) is None:
        sys.exit(f"{ROW} is refused: experiments/harness/TRIES.md must hold exactly one row {ROW} of "
                 f"{hrow.N_CELLS} cells (a row comes before its {what})")
    if not hrow.formed(text):
        sys.exit(f"{ROW} is refused: row {ROW}'s first cell is not the text this script fixes: {FIRST_CELL}")
    if not hrow.is_open(text):
        sys.exit(f"{ROW} is refused: row {ROW}'s outcome cell already holds an outcome (the gate is read once)")
    if vec_home() != HERE:
        sys.exit(f"{ROW} is refused: vec is imported from {vec_home()}, not from this checkout")
    state = git_state()
    if state.get("tree_clean") is not True:
        sys.exit(f"{ROW} is refused: tracked files differ from the commit")
    if not commit or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit {commit!r}")


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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
    if not isinstance(g0, dict) or g0.get("axis") != G0_AXIS:
        return ["v3: no G0 record"]
    if not _of_commit(g0, commit):
        return ["v3: the G0 record is not of the run commit"]
    if g0.get("verdict") != PASS:
        return [f"v3: G0 said {g0.get('verdict')!r}, not {PASS}"]
    return []


def gate_record_reasons(gate: dict | None, commit: str | None) -> list[str]:
    if not isinstance(gate, dict) or gate.get("axis") != GATE_AXIS:
        return ["v3: no gate record"]
    if gate.get("verdict") != PASS:
        return [f"v3: the gate said {gate.get('verdict')!r}, not {PASS}"]
    if not _of_commit(gate, commit):
        return ["v3: the gate record is not of the run commit"]
    return []


def refuse_g0(out: Path, commit: str | None, banked: Path | None) -> None:
    refuse_state(commit, "G0")
    if (Path(out) / G0_JSON).exists() or (Path(out) / G0_TXT).exists():
        sys.exit(f"{ROW} is refused: {out} already holds a G0 record")
    if banked is not None:
        b = Path(banked)
        if b.name != BANKED_NAME or not b.is_file() or sha256_of(b) != BANKED_SHA256:
            sys.exit(f"{ROW} is refused: --banked must be {BANKED_NAME} of sha256 {BANKED_SHA256}")


def refuse_gate(out: Path, commit: str | None, g0_dir: Path) -> None:
    refuse_state(commit, "gate")
    why = g0_record_reasons(_load_json(Path(g0_dir) / G0_JSON), commit)
    if why:
        sys.exit(f"{ROW} is refused: G0 of {commit!r} has not passed ({why[0]})")
    if (Path(out) / GATE_JSON).exists() or (Path(out) / GATE_TXT).exists():
        sys.exit(f"{ROW} is refused: {out} already holds a gate record")


def refuse_run(path: Path | None, commit: str | None, g0_dir: Path, gate_dir: Path) -> None:
    refuse_state(commit)
    why = g0_record_reasons(_load_json(Path(g0_dir) / G0_JSON), commit)
    if why:
        sys.exit(f"{ROW} is refused: G0 of {commit!r} has not passed ({why[0]})")
    if gate_record_reasons(_load_json(Path(gate_dir) / GATE_JSON), commit):
        sys.exit(f"{ROW} is refused: the gate of {commit!r} has not said {PASS} ({Path(gate_dir) / GATE_JSON})")
    if path is not None and Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


def refuse_read(out: Path, commit: str | None) -> None:
    refuse_state(commit, "read")
    if (Path(out) / "reading.json").exists() or (Path(out) / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")


# ------------------------------------------------------------------------------------------------ the fixed words
HEDGE = "A hedge's gate on one truth, not evidence of development."
BUILD = (f"What follows: file WBH is built twice by experiments/builds/heart_wbh_build.sh under ruling row "
         f"{RULING_ROW}, before any score of its own; the figures go to the captain or to fresh analysts, never beside "
         "a test headline, and the file goes in only on the captain's word at the time.")
NOTHING = (f"What follows: nothing is built under row {RULING_ROW}; the line ends there, with no second read of this "
           "blur under this or another row.")


def outcome_line(outcome: str, label: str | None = None, why=()) -> str:
    """The outcome statement, wording fixed (h37_row reads it back): PASS with its label; STOP with IDLE, NOT PASSED
    or the G0 conditions; VOID with every condition found."""
    if outcome == VOID:
        return f"{MARK} {VOID}: {'; '.join(why)}; nothing is built."
    if outcome == PASS:
        return f"{MARK} {PASS} ({label})."
    if label in (IDLE, NOT_PASSED):
        return f"{MARK} {STOP}: {label}."
    return f"{MARK} {STOP} at {label}: {'; '.join(why)}."


def follows_line(outcome: str) -> str:
    return f"{HEDGE} {BUILD if outcome == PASS else NOTHING}"


# ------------------------------------------------------------------------------------------------ shared inputs
def board_genes() -> list[str]:
    return list(load_boards()[BOARD_KEY].genes)


def log_cp10k(counts) -> np.ndarray:
    C = np.asarray(counts.toarray() if hasattr(counts, "toarray") else counts, dtype=np.float64)
    tot = C.sum(1, keepdims=True)
    return np.log1p(np.divide(C * wb.LIBRARY, tot, out=np.zeros_like(C), where=tot > 0))


def stage_x(data_dir: Path, name: str, genes: list[str]) -> np.ndarray:
    """A released stage's log1p CP10k on the board's genes, loaded as the heart builder loads its inputs."""
    a = load_stage(Path(data_dir) / name, genes=genes, drop_duplicates=True)
    return log_cp10k(counts_of(a, genes))


# ------------------------------------------------------------------------------------------------ G0
G0_CONDITION_WORDS = {
    "pairs": "fewer than 3 prior pairs have 2 or more groups",
    "median s": "the prior's median spread is outside [0.10, 1.00]",
    "ko against replicate": "the KO-against-WT pairs' median spread exceeds twice the replicate pair's",
    "g00": "G00's cut (1 - e^(-c/2)) / 2 on E9.5 is under 0.10",
    **{f"banked {k}": f"on the banked half, {k} fails" for k in wb.OUTPUT_CHECKS},
}


def banked_checks(banked: Path, s: np.ndarray, refs: dict) -> dict:
    """G0's output checks of the blur of the banked half at the file's seed: words only (no figure leaves here)."""
    import anndata as ad

    a = ad.read_h5ad(banked)
    try:
        out, err = wb.blur(a.X, np.asarray(a.obsm["spatial_3D"])[:, :3], s, FILE_SEED), None
    except (wb.CountRecoveryError, wb.HoldError, ValueError) as e:
        out, err = None, type(e).__name__
    return wb.output_checks(a.X, out, refs, err)


def g0_reading(prior: wb.Prior, g00: dict | None, banked: dict) -> dict:
    """Every G0 condition and the verdict: PASS when all hold, STOP otherwise (with the outcome line)."""
    checks = {**wb.prior_checks(prior), "g00": bool(g00 and g00.get("pass")),
              **{f"banked {k}": bool(v) for k, v in banked.items()}}
    failed = [k for k, ok in checks.items() if not ok]
    verdict = STOP if failed else PASS
    line = outcome_line(STOP, "G0", [G0_CONDITION_WORDS[k] for k in failed]) if failed else None
    return {"verdict": verdict, "checks": checks, "failed": failed, "line": line,
            "follows": follows_line(STOP) if failed else None}


def g0_text(g: dict) -> str:
    sm = g["prior"]
    lines = [f"Row {ROW}'s G0 (commit {g['commit']}): target-free and not a try. The world prior from row "
             f"{RULING_ROW}'s five files (no E10.5 cell), G00 on the released E9.5, and the blur of the banked half "
             "at the file's seed, of which only words are written.",
             f"Prior: pairs used {sm['pairs_used']}; groups per pair {sm['groups']}; flagged genes (unexpressed in "
             f"every used group, given the median spread) {sm['flagged']} of {sm['n_genes']}; median s "
             f"{sm['median_s']:.4f}; KO-against-WT median spread {sm['ko_wt_median_spread']:.4f} against the "
             f"replicate pair's {sm['replicate_median_spread']:.4f}.",
             ("G00 (E9.5, the scorer's geometry): " + (f"c {g['g00']['c']:.4f}, cut {g['g00']['cut']:.4f} (bar "
                                                       f"{wb.G00_MIN_CUT:.2f})." if g.get("g00") else "not formed."))]
    for k, ok in g["checks"].items():
        lines.append(f"G0 {k}: {'holds' if ok else 'FAILS: ' + G0_CONDITION_WORDS[k]}.")
    lines.append(f"G0 {g['verdict']}.")
    if g.get("line"):
        lines += [g["line"], g["follows"]]
    return "\n".join(lines) + "\n"


def run_g0(out: Path, data_dir: Path, commit: str | None, banked: Path) -> int:
    out = Path(out)
    refuse_g0(out, commit, banked)
    genes = board_genes()
    prior = wb.load_prior(data_dir, genes)
    x95 = stage_x(data_dir, "T2/E9.5.h5ad", genes)
    x105 = stage_x(data_dir, "T2/E10.5.h5ad", genes)
    refs = {"E10.5": x105.mean(0), "E9.5": x95.mean(0)}
    usable = bool(np.isfinite(prior.s).all())
    g00 = wb.g00(x95, prior.s, FILE_SEED) if usable else None
    banked_words = banked_checks(banked, prior.s, refs) if usable else {k: False for k in wb.OUTPUT_CHECKS}
    res = g0_reading(prior, g00, banked_words)
    state = git_state()
    if state.get("tree_clean") is not True or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: tracked files or HEAD changed while G0 ran; no G0 record is written")
    g = {"row": ROW, "axis": G0_AXIS, "commit": commit, "provenance": state,
         "label": "G0 of row H37: target-free, not a try; the prior's and G00's figures come from released files, "
                  "the banked half's checks are words only", "prior": prior.summary(), "g00": g00,
         "s": [float(v) for v in prior.s], "flagged": [bool(v) for v in prior.flagged], **res, "source": SOURCE}
    out.mkdir(parents=True, exist_ok=True)
    (out / G0_JSON).write_text(json.dumps(g, indent=1, default=float))
    (out / G0_TXT).write_text(g0_text(g))
    print(f"G0 {res['verdict']}", flush=True)
    if res["line"]:
        print(res["line"], flush=True)
        print(res["follows"], flush=True)
    return G0_EXIT[res["verdict"]]


def prior_of_record(g0: dict, genes) -> np.ndarray:
    """The spread G0 recorded (the gate and the run take it from there, so every arm uses G0's prior)."""
    s = np.asarray(g0["s"], dtype=np.float64)
    if s.shape != (len(genes),) or not np.isfinite(s).all():
        raise ValueError("G0's record holds no spread for the board's genes")
    return s


# ------------------------------------------------------------------------------------------------ the arms
def build_control(ctx, tr) -> tuple[Stage, dict]:
    """C: row H35's category arm on D's context, then D_f's size and half (h34.size_and_half), as row H35 built it."""
    st, _, _, rec = trend_forms.form_step(ctx, STRENGTH, FORM, tr.whole)
    stage, srec = h34.size_and_half(ctx, st, tr.target, FORM)
    return stage, {"size": srec, "trend": h35.arm_record(rec, ctx)}


def blurred(ctx, c: Stage, s: np.ndarray, name: str) -> tuple[Stage, dict]:
    """An arm from C: the blur of C's X at C's coordinates, then the checked builder at those coordinates."""
    spec = arm_spec(name)
    out = wb.blur(c.X, c.coords, s, ctx.seed, mode=spec["mode"], permute=spec["permute"])
    stage = h29b._checked(ctx, np.asarray(out.X), np.asarray(c.coords), name)
    return stage, out.record


def arms_of(ctx, tr, s: np.ndarray) -> tuple[dict, dict]:
    c, crec = build_control(ctx, tr)
    preds, recs = {CONTROL: c}, {CONTROL: crec}
    for name in ARMS:
        preds[name], recs[name] = blurred(ctx, c, s, name)
    checks = {"coords_equal": {k: bool(np.array_equal(v.coords, c.coords)) for k, v in preds.items()},
              "n_cells": {k: int(v.n) for k, v in preds.items()}}
    return preds, {"records": recs, "checks": checks}


# ------------------------------------------------------------------------------------------------ the gate
_GATE: dict = {}  # the gate's shared inputs, set before its processes fork


def _transition(key: str, data_dir: Path):
    tr = _GATE.setdefault("transitions", {}).get(key)
    if tr is None:
        tr = _GATE["transitions"][key] = h35.load_transition(key, data_dir)
    return tr


def gate_task(job: tuple) -> dict:
    """One transition and seed of the gate: C built and scored; W, U and P built from it and checked, never scored."""
    key, s, data_dir = job
    t0 = time.time()
    tr = _transition(key, Path(data_dir))
    ctx = h35.context(tr, s)
    with harness.guarded_reads(harness.REPO_ROOT, Path(data_dir)):
        c, crec = build_control(ctx, tr)
        refs = {"reference": tr.base.ref_stage.X.astype(np.float64).mean(0), "E9.5": _GATE["pb_e95"]}
        outputs = {}
        for name in ARMS:
            spec = arm_spec(name)
            try:
                b, err = wb.blur(c.X, c.coords, _GATE["s"], ctx.seed, mode=spec["mode"],
                                 permute=spec["permute"]), None
            except (wb.CountRecoveryError, wb.HoldError, ValueError) as e:
                b, err = None, f"{type(e).__name__}: {e}"
            outputs[name] = {"checks": wb.output_checks(c.X, b, refs, err), "record": None if b is None else b.record,
                             "error": err}
    scores, skills, groups = h34.score_arms(tr.base, s, {CONTROL: c})
    return {"task": key, "seed": s, "n_cells": int(c.n), "records": crec, "outputs": outputs, "scores": scores,
            "skills": skills, "groups": groups, "seconds": round(time.time() - t0, 1)}


def h35_committed(root: Path) -> tuple[dict | None, list[str]]:
    """Row H35's committed parts by transition, when its reading.json is the file of H35_READING_SHA256 and every
    part is of row H35's run commit from a clean tree; else None and why."""
    root = Path(root)
    why = []
    reading = root / "reading.json"
    if not reading.is_file() or sha256_of(reading) != H35_READING_SHA256:
        why.append(f"G1: row H35's reading.json is not the file of sha256 {H35_READING_SHA256}")
    parts, unreadable = h35.load_parts(root / "parts")
    if unreadable:
        why.append(f"G1: row H35's parts cannot be read: {unreadable[:2]}")
    every = [p for ps in parts.values() for p in ps]
    if any((p.get("provenance") or {}).get("commit") != H35_RUN_COMMIT
           or (p.get("provenance") or {}).get("tree_clean") is not True for p in every):
        why.append(f"G1: a part of row H35 is not of its run commit {H35_RUN_COMMIT} from a clean tree")
    if any(H35_ARM not in (p.get("arms") or []) for p in every) or set(parts) != set(TASK_KEYS):
        why.append(f"G1: row H35's parts do not hold the {H35_ARM} arm on both transitions")
    if why:
        return None, why
    return {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}, []


def _same_json(x, y) -> bool:
    """Exact equality as the records hold it (a NaN equals a NaN; floats by their shortest repr)."""
    return json.dumps(x, sort_keys=True, default=float) == json.dumps(y, sort_keys=True, default=float)


def g1_reasons(rows: dict[str, list[dict]], committed: dict | None) -> list[str]:
    """G1: per transition, seed and truth, C's scores, skills and groups and cell number are row H35's category arm's,
    and the floor's scores are row H35's, exactly."""
    if committed is None:
        return ["G1: row H35's committed parts are not at hand"]
    why = []
    for key in TASK_KEYS:
        mine = {r["seed"]: r for r in rows.get(key, [])}
        theirs = {r["seed"]: r for r in committed.get(key, [])}
        if sorted(mine) != list(SEEDS) or sorted(theirs) != list(SEEDS):
            why.append(f"G1: {key}: seeds {sorted(mine)} built, {sorted(theirs)} committed, not {list(SEEDS)}")
            continue
        for s in SEEDS:
            a, b = mine[s], theirs[s]
            off = [f"{q} {what}" for q in TRUTHS for what, x, y in (
                ("score", a["scores"][q][CONTROL], b["scores"][q][H35_ARM]),
                ("skills", a["skills"][q][CONTROL], b["skills"][q][H35_ARM]),
                ("groups", a["groups"][q][CONTROL], b["groups"][q][H35_ARM]),
                ("floor", a["scores"][q][FLOOR], b["scores"][q][FLOOR])) if not _same_json(x, y)]
            if a["n_cells"] != b["n_cells"][H35_ARM]:
                off.append("cell number")
            if off:
                why.append(f"G1: {key} seed {s}: C is not row H35's {H35_ARM} arm ({', '.join(off[:3])})")
    return why


def practice_reasons(rows: dict[str, list[dict]]) -> list[str]:
    """G0 on the practice outputs: every arm of every seed holds every output check."""
    why = []
    for key in TASK_KEYS:
        for r in rows.get(key, []):
            for name in ARMS:
                bad = [k for k, ok in r["outputs"][name]["checks"].items() if not ok]
                if bad:
                    why.append(f"{key} seed {r['seed']} {name}: {', '.join(bad)}")
    return why


GATE_LINES = {
    PASS: ("GATE PASS: G1 holds (C rebuilt is row H35's category arm bit for bit, with the floor) and every practice "
           "output holds G0's checks; the run may follow at this commit."),
    FAIL: "GATE FAIL: G1 does not hold; the run step refuses, no try is counted and the outcome is VOID.",
    STOP: ("GATE STOP: a practice output fails G0's checks; the run step refuses, no try is counted and the outcome "
           "is STOP at G0."),
}


def gate_verdict(rows: dict, committed: dict | None, committed_why: list[str]) -> tuple[str, list[str], list[str]]:
    g1 = committed_why + (g1_reasons(rows, committed) if committed is not None else [])
    if not committed_why and committed is None:
        g1.append("G1: row H35's committed parts are not at hand")
    practice = practice_reasons(rows)
    return (FAIL if g1 else STOP if practice else PASS), g1, practice


def gate_text(g: dict) -> str:
    lines = [f"Row {ROW}'s gate (commit {g['commit']}): counting nothing. C is rebuilt and scored on both transitions "
             f"at seeds 0 to 9 and held to row H35's committed parts (reading.json of sha256 {H35_READING_SHA256}); "
             "on each seed's C the arms W, U and P are built, never scored, and held to G0's output checks.",
             f"G1: {len(g['g1'])} conditions found.", f"G0 on the practice outputs: {len(g['practice'])} failures.",
             GATE_LINES[g["verdict"]]]
    if g["verdict"] == FAIL:
        lines += ["WHY: " + "; ".join(g["g1"]), g["line"], g["follows"]]
    elif g["verdict"] == STOP:
        lines += ["WHY: " + "; ".join(g["practice"][:20]), g["line"], g["follows"]]
    return "\n".join(lines) + "\n"


def run_gate(out: Path, data_dir: Path, commit: str | None, g0_dir: Path, h35_root: Path | None = None,
             workers: int = WORKERS) -> int:
    out = Path(out)
    refuse_gate(out, commit, g0_dir)
    committed, committed_why = h35_committed(h35_root or H35_RESULTS)
    genes = board_genes()
    _GATE.update(s=prior_of_record(_load_json(Path(g0_dir) / G0_JSON), genes),
                 pb_e95=stage_x(data_dir, "T2/E9.5.h5ad", genes).mean(0), transitions={})
    jobs = [(key, s, str(data_dir)) for key in TASK_KEYS for s in SEEDS]
    rows: dict[str, list[dict]] = {}
    with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("fork")) as ex:
        for r in ex.map(gate_task, jobs):
            rows.setdefault(r["task"], []).append(r)
            print(f"{r['task']} seed {r['seed']} ({r['seconds']} s) done", flush=True)
    state = git_state()
    if state.get("tree_clean") is not True or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: tracked files or HEAD changed while the gate ran; no gate record is written")
    try:
        verdict, g1, practice = gate_verdict(rows, committed, committed_why)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:  # a malformed committed record
        verdict, g1, practice = FAIL, [f"G1: a committed record cannot be read ({type(e).__name__}: {e})"], []
    line = (outcome_line(VOID, why=[f"v3: the gate fails ({'; '.join(g1)})"]) if verdict == FAIL
            else outcome_line(STOP, "G0 (practice outputs)", practice[:20]) if verdict == STOP else None)
    g = {"row": ROW, "axis": GATE_AXIS, "commit": commit, "provenance": state,
         "label": "the gate of row H37: counting nothing; C alone is scored, and its figures are row H35's own",
         "verdict": verdict, "g1": g1, "practice": practice, "line": line,
         "follows": None if line is None else follows_line(STOP), "h35_reading_sha256": H35_READING_SHA256,
         "rows": rows, "source": SOURCE}
    out.mkdir(parents=True, exist_ok=True)
    (out / GATE_JSON).write_text(json.dumps(g, indent=1, default=float))
    (out / GATE_TXT).write_text(gate_text(g))
    print(f"GATE {verdict}", flush=True)
    if line:
        print(line, flush=True)
        print(g["follows"], flush=True)
    return GATE_EXIT[verdict]


# ------------------------------------------------------------------------------------------------ the run
def run_transition(key: str, seeds, out: Path, data_dir: Path, commit: str | None, g0_dir: Path,
                   gate_dir: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse_run(path, commit, g0_dir, gate_dir)
    if key not in TASK_KEYS:
        sys.exit(f"{key}: not a transition of row {ROW}")
    s = prior_of_record(_load_json(Path(g0_dir) / G0_JSON), board_genes())
    tr = h35.load_transition(key, data_dir)
    rows_out = []
    for seed in seeds:
        t0 = time.time()
        ctx = h35.context(tr, seed)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            preds, extra = arms_of(ctx, tr, s)
        scores, skills, groups = h34.score_arms(tr.base, seed, preds)
        row = {"seed": seed, "build_seed": harness.BUILD_SEED_OFFSET + seed, **extra, "scores": scores,
               "skills": skills, "groups": groups, "seconds": round(time.time() - t0, 1)}
        rows_out.append(row)
        print(f"{key} seed {seed} ({row['seconds']} s) done", flush=True)
    res = {"label": "row H37: C (row H35's category arm on D_f) and the arms W, U and P of the world blur; a gate on "
                    "one truth, not an adoption result", "axis": AXIS, "task": tr.key, "arms": list(BUILT),
           "truth_values": list(TRUTHS), **tr.facts, "strength": STRENGTH, "seeds": list(seeds), "rows": rows_out,
           "source": SOURCE, "provenance": git_state()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


def load_parts(root: Path) -> tuple[dict[str, list[dict]], list[str]]:
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


VOID_LIST = ("v3 G0 or the gate failed, or a record is missing or of another commit than the run commit; v4 a "
             "transition, seed, truth or arm missing or twice, a score not finite, a part unreadable or of another "
             "transition, parts from more than one commit, another commit than the run commit or an unclean tree; v5 "
             "an arm's d2_shape, occupancy_dice or scale_log_ratio more than 1e-9 from C's in any seed and truth, "
             "other coordinates or another cell number than C's, an inexact count recovery, a broken antithetic "
             "pair, a design library more than 2% from C's, a gene's mean X more than 1e-3 from C's, or C or the "
             "floor as the run scores them not the gate's")


def arm_reasons(key: str, rows: list[dict]) -> list[str]:
    why = []
    for r in rows:
        tag = f"v5: {key} seed {r['seed']}"
        ch = r["checks"]
        if not all(ch["coords_equal"].get(k) is True for k in BUILT):
            why.append(f"{tag}: an arm's coordinates are not C's")
        if set(ch["n_cells"]) != set(BUILT) or len(set(ch["n_cells"].values())) != 1:
            why.append(f"{tag}: the arms' cell numbers are not C's")
        for name in ARMS:
            rec = r["records"][name]
            if rec["recovery"]["exact"] is not True:
                why.append(f"{tag}: {name}'s count recovery is not exact")
            if rec["antithetic"] is not True:
                why.append(f"{tag}: {name} breaks an antithetic pair")
            if not abs(rec["library_ratio_design"] - 1) <= wb.LIBRARY_TOL:
                why.append(f"{tag}: {name}'s design library is more than 2% from C's")
            if not rec["pb_max_abs_diff"] <= wb.PB_TOL:
                why.append(f"{tag}: a gene's mean X in {name} is more than {wb.PB_TOL} from C's")
            for q in TRUTHS:
                for m in SHAPE_TERMS:
                    if not abs(r["skills"][q][name][m] - r["skills"][q][CONTROL][m]) <= TOL:
                        why.append(f"{tag} {q}: {name}'s {m} is not C's")
    return why


def gate_tie(rows: dict[str, list[dict]], gate: dict) -> list[str]:
    """v5: C and the floor as the run scores them are exactly the gate's, in every transition, seed and truth."""
    try:
        theirs = {key: {r["seed"]: r["scores"] for r in gate["rows"][key]} for key in TASK_KEYS}
        return [f"v5: {key} seed {r['seed']} {q}: the run's {name} score is not the gate's"
                for key in TASK_KEYS for r in rows[key] for q in TRUTHS for name in (CONTROL, FLOOR)
                if r["scores"][q][name] != theirs[key][r["seed"]][q][name]]
    except (KeyError, TypeError, IndexError, AttributeError) as e:
        return [f"v3: the gate record holds no score to hold the run's to ({type(e).__name__}: {e})"]


def void_reasons(parts: dict[str, list[dict]], g0: dict | None, gate: dict | None, commit: str | None,
                 unreadable=()) -> list[str]:
    v3 = g0_record_reasons(g0, commit) + gate_record_reasons(gate, commit)
    why = v3 + [f"v4: a part cannot be read: {p}" for p in unreadable]
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
            for fld in ("scores", "skills", "groups"):
                gone = sorted({f"{q}/{a}" for q in TRUTHS for a in (FLOOR, *BUILT)
                               if a not in (r.get(fld) or {}).get(q, {})})
                if gone:
                    why.append(f"v4: {key} seed {r['seed']}: missing from {fld}: {gone[:6]}")
                    complete = False
            if set(r.get("records") or {}) != set(BUILT) or "checks" not in r:
                why.append(f"v4: {key} seed {r['seed']}: records or checks missing, or other arms")
                complete = False
            bad = sorted({f"{q}/{a}" for q in TRUTHS for a in (FLOOR, *BUILT)
                          if a in (r.get("scores") or {}).get(q, {}) and not _finite(r["scores"][q][a])})
            if bad:
                why.append(f"v4: {key} seed {r['seed']}: a score is not a finite number: {bad[:6]}")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"v4: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not every:
        return why
    if len({(p.get("provenance") or {}).get("commit") for p in every}) > 1:
        why.append("v4: the parts come from more than one commit")
    if not all((p.get("provenance") or {}).get("tree_clean") is True for p in every):
        why.append("v4: a part comes from an unclean tree")
    if not commit or any(str((p.get("provenance") or {}).get("commit")) != commit for p in every):
        why.append(f"v4: a part comes from a commit other than the run commit {commit!r}")
    if not complete:
        return why
    merged = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
    for key in TASK_KEYS:
        why += arm_reasons(key, merged[key])
    return why + ([] if v3 else gate_tie(merged, gate))


# ------------------------------------------------------------------------------------------------ the bars
PAIRS = ((W, CONTROL), (U, CONTROL), (P, CONTROL), (W, P), (U, W), *((a, FLOOR) for a in BUILT))


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: {pname(a, b): h34.compare(rows, ks, ts, a, b) for a, b in PAIRS} for u, ks in UNITS.items()}
            for w, ts in WORLDS}


def bars(tab: dict) -> dict:
    """The rule on fold E10.5. B1 in a world: W - C at least +0.30 with at least 8 of 10 seeds above 0. B2 in a world:
    W - P at least +0.10. IDLE in a world: |W - C| under 0.10. PASS when B1 holds in both worlds (labelled world
    pattern when B2 holds in both, spread only otherwise); IDLE when IDLE holds in both; NOT PASSED otherwise."""
    def fold(w, a, b):
        return tab[w][FOLD][pname(a, b)]["score"]
    b1 = {w: bool(fold(w, W, CONTROL)["mean"] >= B1_BAR and fold(w, W, CONTROL)["ahead"] >= B1_SEEDS)
          for w, _ in WORLDS}
    b2 = {w: bool(fold(w, W, P)["mean"] >= B2_BAR) for w, _ in WORLDS}
    idle = {w: bool(abs(fold(w, W, CONTROL)["mean"]) < IDLE_BAR) for w, _ in WORLDS}
    if all(b1.values()):
        word, label = PASS, WORLD_PATTERN if all(b2.values()) else SPREAD_ONLY
    elif all(idle.values()):
        word, label = STOP, IDLE
    else:
        word, label = STOP, NOT_PASSED
    return {"b1": b1, "b2": b2, "idle": idle, "outcome": word, "label": label}


def reading(parts: dict[str, list[dict]], g0: dict | None, gate: dict | None, commit: str | None,
            unreadable=()) -> dict:
    try:
        why = void_reasons(parts, g0, gate, commit, unreadable)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:  # a malformed part
        why = [f"v4: a part cannot be read ({type(e).__name__}: {e})"]
    res = {"row": ROW, "commit": commit, "void": bool(why), "why": why, "hedge": HEDGE,
           "g0": {k: (g0 or {}).get(k) for k in ("verdict", "commit", "failed")},
           "gate": {k: (gate or {}).get(k) for k in ("verdict", "commit", "g1", "practice")}}
    if why:
        res.update(outcome=VOID, label=None, line=outcome_line(VOID, why=why))
    else:
        rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
        tab = tables(rows)
        b = bars(tab)
        res.update(outcome=b["outcome"], label=b["label"], bars=b, tables=tab,
                   blur={key: {r["seed"]: {n: r["records"][n] for n in ARMS} for r in rows[key]} for key in TASK_KEYS},
                   line=outcome_line(b["outcome"], b["label"]))
    res["follows"] = follows_line(res["outcome"])
    return res


RULE = (f"RULE, fixed in row {ROW} before its run, on fold E10.5 (the two backtests merged seed by seed), paired by "
        f"seed, in worlds own and CLEAN: B1 W - C at least +{B1_BAR:.2f} with at least {B1_SEEDS} of 10 seeds above "
        f"0, in both worlds; B2 W - P at least +{B2_BAR:.2f} in both worlds. PASS = B1, labelled world pattern when "
        f"B2 holds and spread only when it does not. IDLE (|W - C| under {IDLE_BAR:.2f} in both worlds) and NOT PASSED "
        "are STOP. VOID is decided first, and under VOID nothing is read.")
UNPRICED = ("Cannot show: one truth, so one world draw; whether E12.5 shares E10.5's world (if so the file loses at "
            "test whatever practice says); the practice references are portion samples, the test's a dissected "
            "heart; practice offsets mix development with world. P carries the same count noise as W, so B2 guards "
            "spread and count noise alike. Seed SE 0.010-0.030 (row H28). No figure here is read beside a test score.")


def _st(s: dict) -> str:
    return (f"{s['mean']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10 above, {s['behind']}/10 below "
            f"(se {s['se']:.3f})")


def _terms(c: dict) -> str:
    return ", ".join(f"{m} {c[m]:+.4f}" for m in (*TERMS, *SHAPE_TERMS))


def text(res: dict) -> str:
    lines = [res["line"], res["follows"],
             f"Row {ROW} (commit {res['commit']}): the gate of row {RULING_ROW}'s file WBH on row H29b's D_f form, one "
             "truth (E10.5; the two backtests merged seed by seed into fold E10.5), worlds own and CLEAN, seeds 0 to "
             "9, released stages only. The read builds nothing.",
             f"G0: {res['g0'].get('verdict')} at commit {res['g0'].get('commit')}. GATE: "
             f"{res['gate'].get('verdict')} at commit {res['gate'].get('commit')}.", RULE]
    if res["void"]:
        return "\n".join([*lines, "VOID, the conditions found, one by one:", *res["why"]]) + "\n"
    b = res["bars"]
    lines.append(f"B1 {b['b1']}; B2 {b['b2']}; IDLE {b['idle']}.")
    for w, _ in WORLDS:
        for p in (pname(W, CONTROL), pname(W, P), pname(U, W)):
            lines.append(f"fold {w} {p}: {_st(res['tables'][w][FOLD][p]['score'])}")
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
    sys.exit(run_g0(Path(a.out), Path(a.data), a.commit, Path(a.banked)))


def cmd_gate(a) -> None:
    sys.exit(run_gate(Path(a.out), Path(a.data), a.commit, Path(a.g0), Path(a.h35), a.workers))


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
        if a.banked is None:
            sys.exit("refusals --step g0 needs --banked")
        refuse_g0(Path(a.out), a.commit, Path(a.banked))
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
    z.add_argument("--commit", required=True, help="row H37's run commit (full hash; HEAD must be it)")
    z.add_argument("--banked", required=True, help=f"the banked half, {BANKED_NAME} of its recorded sha256")
    z.add_argument("--data", default=str(harness.DATA_DIR))
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--g0", required=True, help="G0's folder (its g0.json must say PASS at --commit)")
    g.add_argument("--commit", required=True)
    g.add_argument("--data", default=str(harness.DATA_DIR))
    g.add_argument("--h35", default=str(H35_RESULTS), help="row H35's committed results folder")
    g.add_argument("--workers", type=int, default=WORKERS)
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
    f.add_argument("--banked", help="the banked half (g0)")
    a = ap.parse_args(argv)
    {"g0": cmd_g0, "gate": cmd_gate, "run": cmd_run, "read": cmd_read, "refusals": cmd_refusals}[a.cmd](a)


if __name__ == "__main__":
    main()
