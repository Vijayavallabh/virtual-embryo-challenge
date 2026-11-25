#!/usr/bin/env python
"""Row T73 (experiments/harness/TRIES.md): gate G0 of row R58's file T1_s2rcMSBm20 (T1's line 11, the step-change
round's T1 PORTFOLIO design as its cold challenge words the gate). A labelled in-sample diagnostic, never an adoption
result: every task scores one truth (E10.5) over one day, the test spans two (J6).

The question: on row T62's five practice tasks, does MS's momentum at the per-day pace T1_s2rcMSBm20 applies at the test
(1.0 of the released day's change per day: damping 2.0 over two days) beat the paces T1_s2rcMSBm10 (0.5 a day) and
T1_s2rcMSB (0.25 a day) apply, each with the atlas step at MSB's per-day size (1.5 x the practice fraction, 0.3)?

Tasks, draws, seeds, protocol and scoring are row T62's and row T65's (experiments/t1_m1/ms_read.py, imported
unchanged): the harness's T1/E9.5->E10.5:chamber (primary) and T1/E9.5->E10.5, start E8.5, bake protocol, build main,
build seed 1,000,000 + seed, seeds 0-9; row T60's three library tasks with copy_last scored first (a reference, not
counted). One cell draw per task and seed; four arms on it, each MS's step (B0's per-category damped mean momentum
through t1_scope.momentum_shift on the covered cells, stream [seed, 1], eps 0.5, min_cells 20, keep_depth; then the
atlas step on the same cells, stream [seed, 2]; then g_mom as row T62's candidate's, stream [seed, 5]):
  m1rc_gmom_mom  damping 0.5, fraction 0.2: row T65's counted arm, by ms_read.ms_population itself (not counted; the
                 gate: within 0.02 board points of row T65's committed m1rc_gmom_mom on every task and seed)
  ms025_x15      damping 0.25, fraction 0.3: T1_s2rcMSB per day (counted)
  ms05_x15       damping 0.5,  fraction 0.3: T1_s2rcMSBm10 per day (counted)
  ms1_x15        damping 1.0,  fraction 0.3: T1_s2rcMSBm20 per day (counted)
Fifteen tries (three counted arms on five tasks). Every cell outside the covered set keeps the gate arm's bytes in every
counted arm: each run checks this per seed and arm before scoring and stops otherwise.

The statement, fixed in row T73 (d1 = ms1_x15 - ms05_x15 and d2 = ms1_x15 - ms025_x15, paired per seed):
PASS ("m20's step is AHEAD") only if all hold: on the chamber cut the mean of d1 at least +0.3 with at least 8 of 10
seeds above 0; on the whole cut the mean of d1 above 0; no valid library task with a mean of d1 below -0.3; on the
chamber cut the mean of d2 at least +0.3 with at least 8 of 10 seeds above 0. STOP otherwise, and also when two counted
arms are one file on every seed of a task. VOID (nothing follows; the tries count): a task, seed, arm or sha256
missing; task files of more than one commit, not at the row's commit, or from an unclean tree; the gate failing on any
task and seed; a recorded fraction, damping or stream not the row's; the relabelling's statistics, settings, route or
authors_labels_read not the row's; a counted arm changing a cell outside the covered set, or a covered set not the
builder's. A library with inverted anchors on any seed is invalid, not void.

  run-cut  --cut chamber|whole --commit C   one E10.5 cut. Refused before anything loads unless row T73 is in TRIES.md
           naming this script, HEAD is C, the task file does not exist, the tree is clean by vec.provenance's rule,
           the source guards pass and the relabelling's committed checks pass. Prints no score.
  run-lib  --library 1|2|3 --commit C        one library task; the same refusals. Prints no score.
  read     --commit C   the statement, after every void condition; writes t73.json and t73.txt; prints the outcome
           words and what follows only.
  j7       reads committed results only (no row needed, no try, nothing written): the gate's 0.02 tolerance dry-run on
           committed reproduction data other than the pair the row reads, row T65's m1rc_gmom against row T62's
           committed m1rc_gmom, per task (row T73 quotes it before it is committed).

Usage: through experiments/t1_m1/t73_m20.sh on a checkout holding data/derived/feng_relabel_c/ (one step per call).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec.adoption import SEEDS
from vec.boards import get_board, load_boards
from vec.evaluate import Stage
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
TRIES = ROOT / "experiments" / "harness" / "TRIES.md"
sys.path.insert(0, str(HERE))
import ms_read as ms  # noqa: E402  (T65: its draws, MS population, guards and checks)

rc, br = ms.rc, ms.br
gv, pt, rr = rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs

ROW = "T73"
SCRIPT_PATH = "experiments/t1_m1/m20_read.py"
OUT = "t73_m20"
OUT_FILES = ("t73.json", "t73.txt")
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T65_TASKS = RESULTS / "t65_ms" / "tasks"
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"

GATE = ms.MOM  # m1rc_gmom_mom: row T65's counted arm, rebuilt here as the gate (not counted)
MS025, MS05, MS1 = "ms025_x15", "ms05_x15", "ms1_x15"
COUNTED = (MS025, MS05, MS1)
ARMS = (GATE, *COUNTED)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first (it sets the anchors)
LIB_ARMS = (LIB_REFERENCE, *ARMS)
PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2: one day of Feng's E9.5 -> E14.5 bracket
MSB_MULT = 1.5  # T1_s2rcMSB's atlas multiplier (build_test.T1_MSB_MULT)
# arm -> (MS's damping, the atlas step's multiple of the practice fraction)
SETTINGS = {GATE: (0.5, 1.0), MS025: (0.25, MSB_MULT), MS05: (0.5, MSB_MULT), MS1: (1.0, MSB_MULT)}
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 15
REPRO_TOL = ms.REPRO_TOL  # 0.02 board points
FRACTION_TOL = 1e-12
BAR = 0.3
MIN_AHEAD = 8
WHOLE_ABOVE = 0.0
LIB_FLOOR = -0.3
LABEL = ("try T73 (gate G0 of row R58's file T1_s2rcMSBm20), a labelled in-sample diagnostic, never an adoption "
         "result: MS's momentum at 0.25, 0.5 and 1.0 of the released day's change per day with the atlas step at 1.5 x "
         "the practice fraction, on row T62's five tasks, after the gate arm m1rc_gmom_mom reproduces row T65's; every "
         "task scores one truth (E10.5) over one day, the test spans two")
BESIDE = ("J6: the momentum's fit spans one day on every practice task and the test's two (horizon ratio 2), so any "
          "gain is reported as 0 to the practice value.",)
FOLLOWS = {"PASS": ("What follows: row R58's file T1_s2rcMSBm20 is built twice with one sha256 by "
                    "experiments/builds/t1_lm_m20_line.sh. Nothing is built or changed by this read."),
           "STOP": "What follows: no T1_s2rcMSBm20 and no rerun under another row.",
           "VOID": "What follows: nothing; the tries count, no file is built and nothing is rerun under this row."}


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def ms_population_at(pop, start, ref, directions: Mapping, fraction: float, groups: Mapping[str, str],
                     held_m1: Sequence[str], seed: int, shift_kw: Mapping, damping: float):
    """ms_read.ms_population with the momentum's damping as an argument: on the covered cells B0's per-category
    momentum at `damping` (stream [seed, 1]), then the atlas step at `fraction` (stream [seed, 2]); every other cell
    as drawn. At damping 0.5 it is ms_population's arm byte for byte (tested)."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cov = ms.covered_mask(labels, groups, directions)
    parts, rep = [], None
    if cov.any():
        idx = np.flatnonzero(cov)
        mom, rep = t1_scope.momentum_shift(start, ref, pop[idx].copy(), float(damping), list(held_m1),
                                           label_key=t1_scope.LABEL_KEY, seed=[seed, ms.MOMENTUM_STREAM], **shift_kw)
        moved, _ = t1_scope.external_shift(mom, directions, fraction, groups, list(held_m1),
                                           label_key=t1_scope.LABEL_KEY, seed=[seed, ms.DIRECTION_STREAM])
        parts.append((idx, moved))
    if not cov.all():
        idx = np.flatnonzero(~cov)
        parts.append((idx, pop[idx].copy()))
    return t1_scope._stack_rows(parts, pop.n_obs), rep


def arm_population(arm: str, pop, start, ref, directions, fraction, groups, held_m1, seed, shift_kw):
    """The arm's population before g_mom: the gate by ms_read.ms_population itself, the counted arms at their
    settings."""
    damping, mult = SETTINGS[arm]
    if arm == GATE:
        return ms.ms_population(pop, start, ref, directions, fraction, groups, held_m1, seed, shift_kw)
    return ms_population_at(pop, start, ref, directions, mult * fraction, groups, held_m1, seed, shift_kw, damping)


def covered_record(pop, groups: Mapping[str, str], directions: Mapping, gate_X, arm_X, arm: str,
                   rep: Mapping | None) -> dict:
    """One counted arm's covered cells against the gate arm on the same draw: every row whose label is not covered
    has the gate's bytes, or the run stops."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cov = ms.covered_mask(labels, groups, directions)
    a = np.ascontiguousarray(dense(gate_X), dtype=np.float32)
    b = np.ascontiguousarray(dense(arm_X), dtype=np.float32)
    if a.shape != b.shape or a.shape[0] != len(labels):
        raise SystemExit(f"{ROW}: {GATE} and {arm} are not the same cells")
    equal = a[~cov].tobytes() == b[~cov].tobytes()
    if not equal:
        raise SystemExit(f"{ROW}: a cell outside the covered set differs between {GATE} and {arm}; the run stops")
    present = sorted(set(labels))
    return {"covered_set": sorted(lab for lab in present if groups.get(lab) in directions),
            "builder_covered_set": ms.builder_covered(present, directions), "covered_cells": int(cov.sum()),
            "covered_cells_changed": int((a[cov] != b[cov]).any(axis=1).sum()), "non_covered_equal": bool(equal),
            "shift": None if rep is None else list(rep.get("rows") or [])}


def cut_arms(ctx: harness.BuildContext, factors, directions: Mapping, fraction: float,
             held_m1: Sequence[str]) -> tuple[dict[str, Stage], dict]:
    """The four arms on one harness cut and seed, each one build_submission call with ms_read.cut_arms's draw keywords
    and its arm's population then recomputed_read.rc_hook as post_cells, finished as the harness finishes a
    candidate; and the counted arms' covered records against the gate."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    start, start_src = harness._shallow(ctx.inputs["start"]), ctx.task.role_stems["start"]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    t1_scope.classify(start, source=start_src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, directions)
    hook = rc.rc_hook(factors, ctx.seed, held, with_direction, covered)
    shift_kw = {"start_fallback": spec.fallbacks(start_src), "ref_fallback": spec.fallbacks(src)}
    out, seen = {}, {}
    for arm in ARMS:
        got = []

        def post(pop, arm=arm, got=got):
            moved, rep = arm_population(arm, pop, start, ref, directions, fraction, groups, held_m1, ctx.seed, shift_kw)
            got.append((pop, rep))
            return hook(moved)

        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src, post_cells=post)
        out[arm] = harness.finish(ctx, ref, dense(sub.X), None)
        seen[arm] = got[0]
    recs = {}
    for arm in COUNTED:
        pop, rep = seen[arm]
        recs[arm] = covered_record(pop, groups, directions, out[GATE].X, out[arm].X, arm,
                                   None if rep is None else rep.as_dict())
    return out, recs


def library_arms(e85, e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                 n_cells: int = rr.N_CELLS) -> tuple[dict, dict]:
    """The four arms on row T60's oracle draw (ms_read.library_arms's draw and settings), as AnnData, and the counted
    arms' covered records against the gate."""
    base0, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    groups, covered, with_direction = rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), directions)
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    shift_kw = {"within": within, "start_fallback": fb["E8.5"], "ref_fallback": fb["E9.5"], "min_cells": rr.MIN_CELLS,
                "eps": rr.EPS}
    out, reps = {}, {}
    for arm in ARMS:
        pop, rep = arm_population(arm, base0, e85, e95, directions, fraction, groups, held_m1, seed, shift_kw)
        out[arm] = hook(pop)
        reps[arm] = rep
    recs = {arm: covered_record(base0, groups, directions, out[GATE].X, out[arm].X, arm,
                                None if reps[arm] is None else reps[arm].as_dict()) for arm in COUNTED}
    return out, recs


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def row_names_script(tries: Path | None = None) -> bool:
    """TRIES.md holds exactly one row T73, whose first cell starts 'T73 (' and names this script as a whole path."""
    tries = TRIES if tries is None else Path(tries)
    lines = [ln for ln in tries.read_text(encoding="utf-8").splitlines() if ln.startswith(f"| {ROW} ")]
    if len(lines) != 1:
        return False
    first = lines[0][2:].split("|", 1)[0]
    whole = re.search(r"(^|[^A-Za-z0-9_./-])" + re.escape(SCRIPT_PATH) + r"([^A-Za-z0-9_./-]|$)", first)
    return first.startswith(f"{ROW} (") and whole is not None


def refuse(out: Path, name: str, commit: str | None) -> tuple[dict, dict]:
    """The refusals before anything loads (module docstring). Returns the start state and the relabelling checks."""
    if not row_names_script():
        sys.exit(f"{ROW} is refused: TRIES.md has no single row {ROW} whose first cell starts '{ROW} (' and names "
                 f"{SCRIPT_PATH}")
    gv.refuse_existing(out, (name,))
    started = gv.tree_state()
    if not commit or started.get("commit") != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit given by --commit")
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    bad = ms.source_guard(Path(__file__)) + [f"ms_read.py: {x}" for x in ms.source_guard()]
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    checks = rc.relabel_checks()
    if not checks["passes"]:
        sys.exit(f"{ROW} is refused: {checks['reasons']}")
    if ms.DAMPING != 0.5 or SETTINGS[GATE] != (ms.DAMPING, 1.0):
        sys.exit(f"{ROW} is refused: the gate arm is not row T65's (damping 0.5, the practice fraction)")
    return started, checks


def check_fraction(fraction: float, where: str) -> None:
    if abs(fraction - PRACTICE_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")


def _task_fields(rows: Sequence[Mapping], fraction: float) -> dict:
    first = rows[0]["covered"] if rows else {}
    return {"fraction": fraction, "settings": {a: {"damping": SETTINGS[a][0], "fraction": SETTINGS[a][1] * fraction}
                                               for a in ARMS},
            "momentum_stream": ms.MOMENTUM_STREAM, "direction_stream": ms.DIRECTION_STREAM,
            "shift_rows": {a: (first.get(a) or {}).get("shift") for a in COUNTED} if rows else None}


def cmd_run_cut(a) -> None:
    out = RESULTS / a.out / "tasks"
    started, checks = refuse(out, CUT_FILES[a.cut], a.commit)
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[CUTS[a.cut]]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, _ = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                 (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref_stage.n)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    preregistered = t1_scope.held_flat_genes()
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    ref_role = harness.base_role(t)
    directions, ext = rc.recomputed_inputs(list(inputs[ref_role].var_names))
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    check_fraction(fraction, t.key)
    held_m1, union = rc._held_m1(preregistered)
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms, recs = cut_arms(ctx, factors, directions, fraction, held_m1)
        shas = {k: ms.prediction_sha(v.X) for k, v in arms.items()}
        rows.append({"seed": s, "floor": float(f.score),
                     "arms": gv.score_arms(board, target, ref_stage, s, f.anchors, arms), "sha256": shas,
                     "covered": recs, "seconds": round(time.time() - t0, 1)})
        del arms
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows,
           "relabel": rc._relabel_fields(ext, checks), **_task_fields(rows, fraction),
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / CUT_FILES[a.cut]).write_text(json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_lib(a) -> None:
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started, checks = refuse(out, name, a.commit)
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    board = get_board(rr.BOARD)
    spec = t1_scope.load_spec()
    chamber = spec.in_scope("chamber")
    stage = {k: load_stage(data / f, lineage=False) for k, f in rr.FILES.items()}
    for k, s in stage.items():
        if list(s.var_names) != board.genes:
            sys.exit(f"{k}: genes differ from {rr.BOARD}'s gene order")
        t1_scope.classify(s, source=f"{k}_RNA", spec=spec)
    e85, e95, e105 = stage["E8.5"], stage["E9.5"], stage["E10.5"]
    within = {"surrounding": list(rr.CONTAMINATION_LABELS)}
    fb = {k: spec.fallbacks(f"{k}_RNA") for k in rr.FILES}
    held = t1_scope.held_flat_genes()
    ref_stage = pt.to_stage(e95)
    half_a, _ = rr.split_library(e105.obs_names, a.library)
    truth = e105[rr.chamber_rows(e105, half_a, chamber)].copy()
    shares = t1_scope.load_shares(t1_scope.composition(truth), spec)
    truth_stage = pt.to_stage(truth)
    factors = br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    directions, ext = rc.recomputed_inputs(list(e95.var_names))
    fraction = feng_dirs.bracket_fraction(9.5, 10.5)
    check_fraction(fraction, f"library {a.library}")
    held_m1, union = rc._held_m1(held)
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms, recs = library_arms(e85, e95, shares, factors, held, held_m1, s, within, fb, directions, fraction)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        shas = {k: ms.prediction_sha(v.X) for k, v in stages.items()}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s), "sha256": shas,
                     "covered": recs, "seconds": round(time.time() - t0, 1)})
        del arms, stages
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows,
           "relabel": rc._relabel_fields(ext, checks), **_task_fields(rows, fraction),
           "split": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], feng_dirs.RELABEL_CONFIG,
                                       rc.RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def committed_scores(name: str, arm: str, tasks: Path) -> dict[int, float]:
    """seed -> the arm's score in a committed task file."""
    rows = json.loads((tasks / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][arm]["score"]) for r in rows}


def settings_ok(r: Mapping) -> list[str]:
    bad = []
    if not abs(float(r.get("fraction", float("nan"))) - PRACTICE_FRACTION) <= FRACTION_TOL:
        bad.append(f"the practice fraction is {r.get('fraction')}, not {PRACTICE_FRACTION}")
    want = {a: {"damping": SETTINGS[a][0], "fraction": SETTINGS[a][1] * PRACTICE_FRACTION} for a in ARMS}
    got = r.get("settings") or {}
    for arm in ARMS:
        g = got.get(arm) or {}
        if g.get("damping") != want[arm]["damping"] or not abs(
                float(g.get("fraction", float("nan"))) - want[arm]["fraction"]) <= FRACTION_TOL:
            bad.append(f"{arm}'s damping or fraction is not the row's")
    if r.get("momentum_stream") != ms.MOMENTUM_STREAM or r.get("direction_stream") != ms.DIRECTION_STREAM:
        bad.append("a stream is not the row's")
    return bad


def guards(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
           source_bad: Sequence[str], row_commit: str | None = None, seeds: Sequence[int] = SEEDS) -> list[str]:
    """Every void condition of row T73; an empty list means none fired."""
    reasons, commits = [], set()
    for k, r in loaded.items():
        if r is None:
            reasons.append(f"{k}: no task file")
            continue
        rows = r.get("rows") or []
        if sorted(int(x["seed"]) for x in rows) != list(seeds):
            reasons.append(f"{k}: seeds {[x['seed'] for x in rows]}, not {list(seeds)}")
        want = LIB_ARMS if k in LIBRARIES else ARMS
        absent = sorted({a for x in rows for a in want if a not in x["arms"] or a not in (x.get("sha256") or {})})
        if absent or not rows:
            reasons.append(f"{k}: arms {absent} missing")
            continue
        stamp = r.get("stamp") or {}
        commits.add(stamp.get("commit"))
        if stamp.get("tree_clean") is not True:
            reasons.append(f"{k}: its tree was not clean at one commit")
        if row_commit is not None and stamp.get("commit") != row_commit:
            reasons.append(f"{k}: it did not run at the row's commit")
        ref = committed.get(k, {})
        worst = max((abs(br.score_of(x, GATE) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                     for x in rows), default=float("inf"))
        if not worst <= REPRO_TOL:
            reasons.append(f"{k}: the gate failed: {GATE} differs from row T65's committed {GATE} by up to "
                           f"{worst:.3g} (tolerance {REPRO_TOL:g})")
        reasons += [f"{k}: {x}" for x in settings_ok(r)]
        rel = r.get("relabel") or {}
        if rel.get("sha256") != rc.STATS_SHA256 or rel.get("stats_sha256_hashed") != rc.STATS_SHA256:
            reasons.append(f"{k}: ran on another recomputed statistics file")
        if rel.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: ran under another settings file")
        if rel.get("route") != rc.GROUPS_OPTION or rel.get("record_route") != rc.GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        if rel.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the record did not say authors_labels_read false")
        for x in rows:
            for arm in COUNTED:
                c = (x.get("covered") or {}).get(arm) or {}
                if c.get("non_covered_equal") is not True:
                    reasons.append(f"{k}: seed {x['seed']}: {arm}: non_covered_equal is not true")
                if "covered_set" not in c or sorted(c.get("covered_set") or []) != sorted(
                        c.get("builder_covered_set") or []):
                    reasons.append(f"{k}: seed {x['seed']}: {arm}: the covered set is not the builder's")
    if len(commits) > 1:
        reasons.append(f"the task files ran at more than one commit: {sorted(map(str, commits))}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def identical_pairs(r: Mapping) -> list[str]:
    """Pairs of counted arms that are one file on every seed of the task (by prediction sha256)."""
    out = []
    for i, a in enumerate(COUNTED):
        for b in COUNTED[i + 1:]:
            if all(x["sha256"][a] == x["sha256"][b] for x in r["rows"]):
                out.append(f"{a} and {b}")
    return out


def statement(tasks1: Mapping[str, Mapping], tasks2: Mapping[str, Mapping], valid: Sequence[str],
              same: Mapping[str, Sequence[str]]) -> tuple[bool, list[str]]:
    """Row T73's statement: PASS, and the parts that failed (words only). tasks1: d1 = ms1 - ms05 per task; tasks2:
    d2 = ms1 - ms025 per task; same: per task, the pairs of counted arms that are one file on every seed."""
    c1, w1, c2 = tasks1["chamber"], tasks1["whole"], tasks2["chamber"]
    failed = []
    if not (c1["mean"] >= BAR and c1["ahead"] >= MIN_AHEAD):
        failed.append(f"on the chamber cut {MS1} - {MS05} is not +{BAR:g} or more with {MIN_AHEAD} of 10 seeds ahead")
    if not w1["mean"] > WHOLE_ABOVE:
        failed.append(f"on the whole cut the mean of {MS1} - {MS05} is not above 0")
    low = [lib for lib in valid if tasks1[lib]["mean"] < LIB_FLOOR]
    if low:
        failed.append(f"a valid library task's mean of {MS1} - {MS05} is below {LIB_FLOOR:g} (library "
                      f"{', '.join(low)})")
    if not (c2["mean"] >= BAR and c2["ahead"] >= MIN_AHEAD):
        failed.append(f"on the chamber cut {MS1} - {MS025} is not +{BAR:g} or more with {MIN_AHEAD} of 10 seeds "
                      "ahead")
    for k, pairs in same.items():
        for p in pairs:
            failed.append(f"{p} are one file on every seed of {k}")
    return not failed, failed


def t73_reading(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
                source_bad: Sequence[str] = (), row_commit: str | None = None) -> dict:
    """Row T73's reading: the void conditions (the gate among them) first, then the statement."""
    reasons = guards(loaded, committed, checks, source_bad, row_commit)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons, "follows": FOLLOWS["VOID"]}
    d1 = {k: ms.diff(r["rows"], MS1, MS05) for k, r in loaded.items()}
    d2 = {k: ms.diff(r["rows"], MS1, MS025) for k, r in loaded.items()}
    d3 = {k: ms.diff(r["rows"], MS05, MS025) for k, r in loaded.items()}
    valid = [lib for lib in LIBRARIES if ms.library_valid(loaded[lib])]
    same = {k: identical_pairs(r) for k, r in loaded.items()}
    ahead, failed = statement(d1, d2, valid, {k: v for k, v in same.items() if v})
    outcome = "PASS" if ahead else "STOP"
    return {"outcome": outcome, "ahead": ahead, "failed": failed, "valid_libraries": valid, "ms1_minus_ms05": d1,
            "ms1_minus_ms025": d2, "ms05_minus_ms025": d3, "identical": same,
            "j6_range": sorted([0.0, d1["chamber"]["mean"]]), "beside": list(BESIDE), "follows": FOLLOWS[outcome]}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return [f"Outcome: {ROW} VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]]
    first = (f"Outcome: {ROW} PASS (m20's step is AHEAD)" if r["ahead"] else
             f"Outcome: {ROW} STOP (" + "; ".join(r["failed"]) + ")")
    return [first, *r["beside"], r["follows"]]


def t73_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        for title, key in ((f"d1 = {MS1} - {MS05} (the statement's), paired per seed on one draw:", "ms1_minus_ms05"),
                           (f"d2 = {MS1} - {MS025} (the statement reads the chamber cut only):", "ms1_minus_ms025"),
                           (f"{MS05} - {MS025} (reported, deciding nothing):", "ms05_minus_ms025")):
            lines.append(title)
            lines += [ms._line(k, p) for k, p in r[key].items()]
        lines.append(f"Valid libraries: {', '.join(r['valid_libraries']) or 'none'}")
        lo, hi = r["j6_range"]
        lines.append(f"J6 range (0 to the chamber cut's mean of d1): {lo:+.2f} to {hi:+.2f}")
    lines.append(f"Tries counted: {res['tries']}")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, row_commit: str | None) -> dict:
    if not row_names_script():
        sys.exit(f"{ROW} is refused: TRIES.md has no single row {ROW} naming {SCRIPT_PATH}")
    gv.refuse_existing(out, OUT_FILES)
    if not row_commit:
        sys.exit(f"{ROW} is refused: the read needs the row's commit (--commit)")
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    paths = {k: out / "tasks" / f for k, f in rc.task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    committed = {k: committed_scores(f, GATE, T65_TASKS) for k, f in rc.task_files().items()}
    checks = rc.relabel_checks(need_stats=False)
    bad = ms.source_guard(Path(__file__)) + [f"ms_read.py: {x}" for x in ms.source_guard()]
    reading = t73_reading(loaded, committed, checks, bad, row_commit)
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT, "row_commit": row_commit}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / OUT_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / OUT_FILES[1]).write_text(t73_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def j7(t65: Path = T65_TASKS, t62: Path = T62_TASKS) -> dict:
    """The gate's tolerance dry-run on committed reproduction data other than the pair row T73 reads: row T65's
    m1rc_gmom against row T62's committed m1rc_gmom, per task: the largest per-seed difference and the seeds on both
    sides. Committed results only; nothing is written."""
    out = {}
    for k, f in rc.task_files().items():
        a, b = committed_scores(f, ms.BASE, t65), committed_scores(f, ms.BASE, t62)
        both = sorted(set(a) & set(b))
        out[k] = {"seeds_compared": len(both), "largest_difference": max((abs(a[s] - b[s]) for s in both),
                                                                         default=float("nan"))}
    worst = max(v["largest_difference"] for v in out.values())
    return {"tasks": out, "largest_difference": worst, "tolerance": REPRO_TOL,
            "holds": bool(all(v["seeds_compared"] == len(SEEDS) for v in out.values()) and worst <= REPRO_TOL)}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T73")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", default=OUT)
    p.add_argument("--commit", required=True, help="the run commit (full hash): HEAD must be it")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T73")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", default=OUT)
    p.add_argument("--commit", required=True, help="the run commit (full hash): HEAD must be it")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T73's reading")
    p.add_argument("--out", default=OUT)
    p.add_argument("--commit", required=True, help="the row's commit: every task file's stamp must be at it")
    sub.add_parser("j7", help="the gate's tolerance dry-run on committed results (rows T65 and T62)")
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out, a.commit)
    elif a.cmd == "j7":
        res = j7()
        for k, v in res["tasks"].items():
            print(f"{k}: {v['seeds_compared']} seeds, largest difference {v['largest_difference']:.4g}")
        print(f"largest {res['largest_difference']:.4g}, tolerance {REPRO_TOL:g}: holds {res['holds']}")
    else:
        {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)


if __name__ == "__main__":
    main()
