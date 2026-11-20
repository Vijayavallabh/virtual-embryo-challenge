#!/usr/bin/env python
"""Try T64 (experiments/harness/TRIES.md): candidate D of decision log 2026-11-20 (r), read on row T62's five practice
tasks. D is the recomputed-groups file (row R13) with the team's released per-category momentum on the cells that
ruling R1 holds (the drawn cells whose label has no recomputed direction).

The practice analogue of R13's file is row T62's `m1rc_gmom`, taken unchanged from `recomputed_read.cut_pair` and
`recomputed_read.library_pair` (their CANDIDATE; their other arms are built and not scored). It is the gate and the
base, not counted. The counted arm `m1rc_gmom_heldmom` is that base with one change: the drawn cells whose label has no
recomputed direction get B0's per-category damped momentum at harness.B0_DAMPING (0.5),
    log f = 0.5 x log((m_later + 0.5) / (m_earlier + 0.5)) per category, over the category's own cells,
f = 1 on R2's held genes (the direction step's list), keep_depth, on build_submission's momentum stream [seed, 1]; a
category with fewer than 20 own cells, or drawn from stand-ins at either stage, keeps f = 1 (`t1_scope.momentum_shift`).
Then g_mom as the base's (`recomputed_read.rc_hook`: the same factors, stream [seed, 5], the same mask, which holds the
held labels, so they draw nothing there). Every other cell keeps the base's bytes; each run checks that per seed before
scoring and stops otherwise.
  Cuts: one `t1_scope.build_submission` call with recomputed_read.cut_pair's keywords for its candidate plus start
  (the task's start stage, E8.5), start_source and damping (T61's b0_gmom route). The fit runs E8.5 -> E9.5.
  Libraries: T60's oracle draw, `gmom_combo.moved_by_feng` on the covered rows (stream [seed, 2]), `momentum_shift` from
  released E8.5 to released E9.5 with T61's library settings on the uncovered rows, stacked by `t1_scope._stack_rows`
  (tested equal to build_submission's path on the same draw). E10.5 enters no fit.
It imports recomputed_read.py (T62) and, through it, bets_read.py (T61) and gmom_combo.py (T26), and copies none of
their arithmetic. No builder option exists for D; none is written here.

  run-cut  One E10.5 cut. Refused before anything loads unless row T64 is in TRIES.md, the output file does not exist,
           the tree is clean by vec.provenance's rule, this script's and recomputed_read.py's source guards pass, the
           relabelling's committed checks pass (`recomputed_read.relabel_checks`) and the damping is 0.5. Prints no
           score.
  run-lib  One library task; the same refusals. Prints no score.
  read     T64's statement, applied before anything is printed; writes t64.json and t64.txt; prints the outcome words
           and what follows only.

Usage: through experiments/t1_m1/t64_heldmom.sh on the DGX campaign checkout ~/v/vec-campaign-d (it alone holds
data/derived/feng_relabel_c/), at the row's commit, committing nothing until the five task files exist.
"""

from __future__ import annotations

import argparse
import hashlib
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
import recomputed_read as rc  # noqa: E402  (T62: its candidate's arms, guards and checks)

br, gc = rc.br, rc.gc
bt, gv, pt, rr = rc.bt, rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs

ROW = "T64"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"
OUT_FILES = ("t64.json", "t64.txt")

BASE = rc.CANDIDATE  # m1rc_gmom: the gate and the base, not counted
ARM = "m1rc_gmom_heldmom"
ARMS = (BASE, ARM)
GATED = (BASE,)
COUNTED = (ARM,)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first (it sets the anchors)
LIB_ARMS = (LIB_REFERENCE, BASE, ARM)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 5
FIRST_TRY = 413  # T63 ended at T1 413
TRY_RANGE = (FIRST_TRY, FIRST_TRY + TRIES_COUNT)  # T1 413 to 418

DAMPING = harness.B0_DAMPING  # 0.5, checked before anything loads
MOMENTUM_STREAM = 1  # build_submission's momentum stream [seed, 1]
REPRO_TOL = br.REPRO_TOL  # 0.02 board points: m1rc_gmom against row T62's committed m1rc_gmom, per task and seed
FRACTION_TOL = 1e-12
AHEAD_BAR = br.AHEAD_BAR  # 0.3
MIN_AHEAD = br.MIN_AHEAD  # 8 of 10
WHOLE_ABOVE = 0.0  # T61's own clause: the whole cut's mean strictly above 0
LIB_FLOOR = -0.3
MIN_LIBRARIES = br.MIN_LIBRARIES  # 2
PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2
GROUP_WEIGHTS = rc.GROUP_WEIGHTS
EXPECTED_HELD = "Blood"  # on practice the only held label E9.5 holds; any other drawn held label is named
BANNED_IMPORTS = ("perm_read", "eb_read", "step_read", "deep_read", "over_read", "lib_step_read", "ctl_range_read",
                  "n48_m1_direction_shrink", "n12_hygiene_size", "double_step_read")
BANNED_PREFIX = "n5" + "4"  # no script of row N54 (spelled in pieces so that recomputed_read's guard passes this file)

AHEAD_WORDS = {True: "D is AHEAD", False: "D is NOT AHEAD"}
BESIDE = ("The arm reverses R1 (the captain's ruling of 2026-11-05) for the held cells: only his ruling naming R1 puts "
          "it in a file.",
          "On practice the arm moved Blood cells alone; Macrophage and Early VIC, which D also moves at the test, are "
          "unpriced.")
FOLLOWS_AHEAD = ("What follows: one question to the captain, QD: reverse R1 for this one file, or keep it. It is "
                 "priced both ways, its default is keep R1, and only a reply naming R1 reverses it; then the builder "
                 "option is written, reviewed and rehearsed before any build. Nothing is built or changed by this "
                 "read.")
FOLLOWS_END = ("What follows: candidate D ends. No question goes to the captain, nothing is built, and nothing is "
               "rerun under another row.")
QD_FIXED = (
    "Beside QD (fixed in row T64 before the run): the library figures are bounds, never a centre (they sit on a mix "
    "chosen with the truth's proportions, J3); the practice priced Blood alone; R1's own readout on the same truth "
    "found no sign skill on Blood (sign agreement 0.50 full and 0.55 cross-fitted, E8.5 -> E9.5 momentum against the "
    "E10.5 change); row T63's placebo of the atlas step's size gained +0.30 on the chamber cut.",
    "An earlier read, voided because its reproduction gate failed, had the form before your ruling ahead of R1 "
    "together with a gene hold by only 0.16 / 0.13 points. The two changes cannot be separated there. The gain is "
    "small, and it points towards this file.")
LABEL = ("try T64, a labelled in-sample diagnostic, never an adoption result: candidate D of decision log 2026-11-20 "
         "(r), R13's file with the team's released per-category momentum on the cells ruling R1 holds, by its practice "
         "analogue m1rc_gmom_heldmom against m1rc_gmom (row T62's) on T62's five tasks; every task scores one truth "
         "(E10.5) over one day, the test spans two; on practice the held cells drawn are Blood alone")


# ---------------------------------------------------------------------------------------------------------------------
# Guards on this source


def source_guard(path: Path | None = None) -> list[str]:
    """recomputed_read's guard on this file (no other route to the atlas groups, no token of row N54's listing), and no
    import of a script that forms the atlas direction on another route (BANNED_IMPORTS), of any row-N54 script, or any
    dynamic import. Returns the failures."""
    path = Path(__file__ if path is None else path)
    bad = rc.source_guard(path)
    text = path.read_text(encoding="utf-8")
    for m in re.finditer(r"^\s*(?:from\s+(\S+)\s+import\s+([^\n#]+)|import\s+([^\n#]+))", text, re.M):
        names = [x for x in re.split(r"[\s,.()]+", " ".join(x for x in m.groups() if x)) if x]
        hit = sorted(set(names) & set(BANNED_IMPORTS)) + sorted(x for x in names if x.startswith(BANNED_PREFIX))
        if hit:
            bad.append(f"the source imports {hit}")
    for token in ("import" + "lib", "__imp" + "ort__", "spec_from" + "_file_location", "T1_s2" + ".h5ad"):
        if token in text:
            bad.append(f"the source names {token}")
    return bad


# ---------------------------------------------------------------------------------------------------------------------
# The held cells


def held_set(labels, directions: Mapping) -> set[str]:
    """The labels R1 holds under --feng-groups recomputed: no recomputed direction for their group (the builder's
    `_labels_without_direction`)."""
    return set(bt._labels_without_direction(labels, sorted(directions)))


def uncovered_set(labels, groups: Mapping[str, str], directions: Mapping) -> set[str]:
    """The labels build_submission leaves uncovered (its `direction_groups.get(label) in directions` test)."""
    return {str(lab) for lab in labels if groups.get(str(lab)) not in directions}


def prediction_sha(X) -> str:
    """The sha256 of an arm's scored matrix (cells x genes, float32, C order)."""
    return hashlib.sha256(np.ascontiguousarray(dense(X), dtype=np.float32).tobytes()).hexdigest()


def held_record(pop, groups: Mapping[str, str], directions: Mapping, base_X, arm_X, shift: Mapping | None) -> dict:
    """The per-seed record of the held cells, after the byte check: every row whose label (pop's obs, in row order) is
    outside the held set has the base's bytes in the arm, or the run stops."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    present = sorted(set(labels))
    hs, us = held_set(present, directions), uncovered_set(present, groups, directions)
    held = np.isin(labels, sorted(hs))
    a = np.ascontiguousarray(dense(base_X), dtype=np.float32)
    b = np.ascontiguousarray(dense(arm_X), dtype=np.float32)
    if a.shape != b.shape or a.shape[0] != len(labels):
        raise SystemExit(f"{ROW}: the two arms are not the same cells ({a.shape}, {b.shape}, {len(labels)} labels)")
    equal = a[~held].tobytes() == b[~held].tobytes()
    if not equal:
        raise SystemExit(f"{ROW}: a cell outside the held set differs between {BASE} and {ARM}; the run stops")
    drawn = {lab: int((labels == lab).sum()) for lab in sorted(hs)}
    return {"held_set": sorted(hs), "uncovered_set": sorted(us), "held_labels_drawn": drawn,
            "held_cells": int(held.sum()), "held_cells_changed": int((a[held] != b[held]).any(axis=1).sum()),
            "non_held_equal": bool(equal), "shift": None if shift is None else list(shift.get("rows") or [])}


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def cut_arms(ctx: harness.BuildContext, factors, directions: Mapping, fraction: float,
             held_m1: Sequence[str]) -> tuple[dict[str, Stage], dict]:
    """m1rc_gmom (recomputed_read.cut_pair's CANDIDATE, as called there) and m1rc_gmom_heldmom (one build_submission
    call with cut_pair's keywords for its candidate plus start, start_source and damping) on one harness cut and seed,
    each finished as the harness finishes a candidate, and the held cells' record."""
    base = rc.cut_pair(ctx, factors, directions, fraction, held_m1)[BASE]
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, directions)
    hook, seen = rc.rc_hook(factors, ctx.seed, held, with_direction, covered), []

    def post(pop):
        seen.append(pop)  # the drawn and moved population before g_mom, with its obs labels in row order
        return hook(pop)

    sub, report = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src,
                                            held_flat=list(held_m1), directions=directions, direction_fraction=fraction,
                                            direction_groups=groups, start=harness._shallow(ctx.inputs["start"]),
                                            start_source=ctx.task.role_stems["start"], damping=DAMPING,
                                            post_cells=post)
    arm = harness.finish(ctx, ref, dense(sub.X), None)
    return {BASE: base, ARM: arm}, held_record(seen[0], groups, directions, base.X, arm.X, report.get("shift"))


def library_population(e85, e95, half_shares: Mapping[str, float], held_m1: Sequence[str], seed: int, within: Mapping,
                       fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                       n_cells: int = rr.N_CELLS) -> tuple:
    """The arm's population before g_mom on a library task: T60's oracle draw, the direction on the covered rows
    (`gmom_combo.moved_by_feng`, stream [seed, 2]) and B0's momentum from released E8.5 to released E9.5 on the
    uncovered rows (T61's library settings, stream [seed, MOMENTUM_STREAM]). Returns (population, momentum report or
    None, label -> group map, covered labels, groups with a direction)."""
    base0, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    groups, covered, with_direction = rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), directions)
    moved = gc.moved_by_feng(base0, directions, fraction, groups, list(held_m1), [seed, 2])
    labels = base0.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cov = np.array([groups.get(lab) in directions for lab in labels], dtype=bool)
    idx, cov_idx = np.flatnonzero(~cov), np.flatnonzero(cov)
    if not idx.size:
        return moved, None, groups, covered, with_direction
    rest, rep = t1_scope.momentum_shift(e85, e95, base0[idx].copy(), DAMPING, list(held_m1),
                                        seed=[seed, MOMENTUM_STREAM], within=within, start_fallback=fb["E8.5"],
                                        ref_fallback=fb["E9.5"], min_cells=rr.MIN_CELLS, eps=rr.EPS)
    parts = ([(cov_idx, moved[cov_idx].copy())] if cov_idx.size else []) + [(idx, rest)]
    return t1_scope._stack_rows(parts, base0.n_obs), rep, groups, covered, with_direction


def library_arms(e85, e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                 n_cells: int = rr.N_CELLS) -> tuple[dict, dict]:
    """m1rc_gmom (recomputed_read.library_pair's CANDIDATE, as called there) and m1rc_gmom_heldmom on the same draw,
    as AnnData, and the held cells' record."""
    base = rc.library_pair(e95, half_shares, factors, held, held_m1, seed, within, fb, directions, fraction,
                           n_cells=n_cells)[BASE]
    pop, rep, groups, covered, with_direction = library_population(e85, e95, half_shares, held_m1, seed, within, fb,
                                                                   directions, fraction, n_cells=n_cells)
    arm = rc.rc_hook(factors, seed, held, with_direction, covered)(pop)
    rec = held_record(pop, groups, directions, base.X, arm.X, None if rep is None else rep.as_dict())
    return {BASE: base, ARM: arm}, rec


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def refuse(out: Path, name: str) -> tuple[dict, dict]:
    """The refusals before anything loads: the row, an existing output file, an unclean tree, both source guards, the
    relabelling's committed checks and the damping. Returns the start state and the relabelling's checks."""
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, (name,))
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    checks = rc.relabel_checks()
    if not checks["passes"]:
        sys.exit(f"{ROW} is refused: {checks['reasons']}")
    if DAMPING != 0.5:
        sys.exit(f"{ROW} is refused: the damping is {DAMPING}, not the row's 0.5 (harness.B0_DAMPING)")
    return started, checks


def check_fraction(fraction: float, where: str) -> None:
    if abs(fraction - PRACTICE_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")


def _constants(held_m1: Sequence[str]) -> dict:
    return {"edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1),
            "momentum_eps": rr.EPS, "momentum_min_cells": rr.MIN_CELLS}


def cmd_run_cut(a) -> None:
    out = RESULTS / a.out / "tasks"
    started, checks = refuse(out, CUT_FILES[a.cut])
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
        arms, rec = cut_arms(ctx, factors, directions, fraction, held_m1)
        shas = {k: prediction_sha(v.X) for k, v in arms.items()}
        rows.append({"seed": s, "floor": float(f.score), "arms": gv.score_arms(board, target, ref_stage, s, f.anchors,
                                                                               arms),
                     "sha256": shas, "held": rec, "seconds": round(time.time() - t0, 1)})
        del arms
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "damping": DAMPING, "momentum_stream": MOMENTUM_STREAM, "relabel": rc._relabel_fields(ext, checks),
           "shift_rows": rows[0]["held"]["shift"] if rows else None, "constants": _constants(held_m1),
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / CUT_FILES[a.cut]).write_text(json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_lib(a) -> None:
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started, checks = refuse(out, name)
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
        arms, rec = library_arms(e85, e95, shares, factors, held, held_m1, s, within, fb, directions, fraction)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        shas = {k: prediction_sha(v.X) for k, v in stages.items()}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s), "sha256": shas,
                     "held": rec, "seconds": round(time.time() - t0, 1)})
        del arms, stages
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "damping": DAMPING, "momentum_stream": MOMENTUM_STREAM, "relabel": rc._relabel_fields(ext, checks),
           "shift_rows": rows[0]["held"]["shift"] if rows else None,
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD, **_constants(held_m1)},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], feng_dirs.RELABEL_CONFIG,
                                       rc.RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def task_files() -> dict[str, str]:
    return rc.task_files()


def committed_t62(name: str, tasks: Path | None = None) -> dict[int, float]:
    """seed -> score of m1rc_gmom in row T62's committed task file."""
    rows = json.loads(((T62_TASKS if tasks is None else tasks) / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][BASE]["score"]) for r in rows}


def library_valid(r: Mapping) -> bool:
    return not any(x.get("inverted_anchors") for x in r["rows"])


def idle(r: Mapping) -> bool:
    """The two arms are one file on every seed (by prediction sha256)."""
    return all(x["sha256"][ARM] == x["sha256"][BASE] for x in r["rows"])


def guards(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
           source_bad: Sequence[str], seeds: Sequence[int] = SEEDS) -> list[str]:
    """Every void condition of row T64; an empty list means none fired."""
    reasons = []
    commits = set()
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
        ref = committed.get(k, {})
        worst = max((abs(br.score_of(x, BASE) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                     for x in rows), default=float("inf"))
        if not worst <= REPRO_TOL:
            reasons.append(f"{k}: the gate failed: {BASE} differs from row T62's committed {BASE} by up to "
                           f"{worst:.3g} (tolerance {REPRO_TOL:g})")
        if not abs(float(r.get("fraction", float("nan"))) - PRACTICE_FRACTION) <= FRACTION_TOL:
            reasons.append(f"{k}: the practice fraction is {r.get('fraction')}, not {PRACTICE_FRACTION}")
        if r.get("damping") != DAMPING:
            reasons.append(f"{k}: the damping is {r.get('damping')}, not {DAMPING}")
        if r.get("momentum_stream") != MOMENTUM_STREAM:
            reasons.append(f"{k}: the momentum stream is {r.get('momentum_stream')}, not {MOMENTUM_STREAM}")
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
            h = x.get("held") or {}
            if h.get("non_held_equal") is not True:
                reasons.append(f"{k}: seed {x['seed']}: non_held_equal is not true")
            if sorted(h.get("held_set") or []) != sorted(h.get("uncovered_set") or []) or "held_set" not in h:
                reasons.append(f"{k}: seed {x['seed']}: the held set is not build_submission's uncovered set")
    if len(commits) > 1:
        reasons.append(f"the task files ran at more than one commit: {sorted(map(str, commits))}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def diff(rows: Sequence[Mapping]) -> dict:
    """d = m1rc_gmom_heldmom - m1rc_gmom per seed (`bets_read.paired`: mean, sd, seeds above 0, range, points by group),
    with the seeds below 0 and both arms' mean scores."""
    p = br.paired(rows, ARM, BASE)
    p["behind"] = int(sum(br.score_of(r, ARM) - br.score_of(r, BASE) < 0 for r in rows))
    p["mean_scores"] = {a: float(np.mean([br.score_of(r, a) for r in rows])) for a in ARMS}
    return p


def statement(tasks: Mapping[str, Mapping], valid: Sequence[str]) -> tuple[bool, list[str]]:
    """Row T64's statement (T61's rule): AHEAD, and the parts that failed (words only, no figure)."""
    c, w = tasks["chamber"], tasks["whole"]
    failed = []
    if not c["mean"] > AHEAD_BAR:
        failed.append(f"the chamber cut's mean is not above +{AHEAD_BAR:g}")
    if not c["ahead"] >= MIN_AHEAD:
        failed.append(f"fewer than {MIN_AHEAD} of 10 seeds are above 0 on the chamber cut")
    if not w["mean"] > WHOLE_ABOVE:
        failed.append("the whole cut's mean is not above 0")
    if len(valid) < MIN_LIBRARIES:
        failed.append(f"fewer than {MIN_LIBRARIES} library tasks are valid")
    low = [lib for lib in valid if tasks[lib]["mean"] < LIB_FLOOR]
    if low:
        failed.append(f"a valid library task's mean is below {LIB_FLOOR:g} (library {', '.join(low)})")
    return not failed, failed


def held_summary(r: Mapping) -> dict:
    """A task's held cells: labels drawn summed over seeds, cells changed per seed, and the momentum's rows per
    category (the first seed's report)."""
    drawn: dict[str, int] = {}
    for x in r["rows"]:
        for lab, n in (x["held"].get("held_labels_drawn") or {}).items():
            drawn[lab] = drawn.get(lab, 0) + int(n)
    keys = ("type", "status", "n_target", "n_start", "n_ref", "start_cells", "ref_cells", "n_up", "n_down",
            "mean_abs_log_f", "max_abs_log_f", "median_depth_scale")
    return {"held_labels_drawn": drawn, "held_cells_changed": [int(x["held"]["held_cells_changed"]) for x in r["rows"]],
            "momentum_rows": [{k: v for k, v in row.items() if k in keys} for row in (r.get("shift_rows") or [])]}


def t64_reading(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
                source_bad: Sequence[str] = ()) -> dict:
    """Row T64's reading: the void conditions (the gate among them) first, then the statement."""
    reasons = guards(loaded, committed, checks, source_bad)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons, "follows": FOLLOWS_END}
    tasks = {k: diff(r["rows"]) for k, r in loaded.items()}
    valid = [lib for lib in LIBRARIES if library_valid(loaded[lib])]
    ahead, failed = statement(tasks, valid)
    held = {k: held_summary(r) for k, r in loaded.items()}
    other = {k: sorted(lab for lab in h["held_labels_drawn"] if lab != EXPECTED_HELD) for k, h in held.items()}
    return {"outcome": AHEAD_WORDS[ahead], "ahead": ahead, "failed": failed, "valid_libraries": valid,
            "tasks": tasks, "held": held, "idle_tasks": [k for k, r in loaded.items() if idle(r)],
            "other_held_labels": {k: v for k, v in other.items() if v},
            "j6_range": sorted([0.0, tasks["chamber"]["mean"]]), "beside": list(BESIDE),
            "follows": FOLLOWS_AHEAD if ahead else FOLLOWS_END}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]]
    first = f"Outcome: {r['outcome']}" + ("" if r["ahead"] else " (" + "; ".join(r["failed"]) + ")")
    lines = [first, *r["beside"]]
    lines += [f"the arm moved no cell on {k}: nothing of D is priced there" for k in r["idle_tasks"]]
    lines += [f"held labels other than Blood drawn on {k}: {', '.join(v)}" for k, v in r["other_held_labels"].items()]
    return lines + [r["follows"]]


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} behind {p['behind']} of "
            f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}; mean scores "
            f"{BASE} {p['mean_scores'][BASE]:.2f}, {ARM} {p['mean_scores'][ARM]:.2f}")


def t64_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        lines.append(f"d = {ARM} - {BASE}, paired per seed on the same draw:")
        lines += [_line(k, p) for k, p in r["tasks"].items()]
        for k, h in r["held"].items():
            lines.append(f"  {k}: held labels drawn over the ten seeds {h['held_labels_drawn']}; held cells changed "
                         f"per seed {h['held_cells_changed']}")
            for m in h["momentum_rows"]:
                lines.append(f"    momentum {m.get('type')}: {m.get('status')}; own cells {m.get('n_start')} at the "
                             f"start, {m.get('n_ref')} at the reference ({m.get('start_cells')}, {m.get('ref_cells')});"
                             f" mean |log f| {m.get('mean_abs_log_f', 0.0):.3f}")
        lines.append(f"Valid libraries: {', '.join(r['valid_libraries']) or 'none'}; tasks where the arm moved no "
                     f"cell: {', '.join(r['idle_tasks']) or 'none'}")
        lines.append(f"J6 range (0 to the chamber cut's mean; a negative mean is reported as [mean, 0]): "
                     f"{r['j6_range'][0]:+.2f} to {r['j6_range'][1]:+.2f}")
        if r["ahead"]:
            means = ", ".join(f"{k} {p['mean']:+.2f}" for k, p in r["tasks"].items())
            lines.append(f"QD, priced both ways: the five means {means}; J6 {r['j6_range'][0]:+.2f} to "
                         f"{r['j6_range'][1]:+.2f}")
            lines += list(QD_FIXED)
    lines.append(f"Tries counted: {res['tries']} (T1 {TRY_RANGE[0]} to {TRY_RANGE[1]})")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path) -> dict:
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md")
    gv.refuse_existing(out, OUT_FILES)
    started = gv.tree_state()
    paths = {k: out / "tasks" / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    committed = {k: committed_t62(f) for k, f in task_files().items()}
    checks = rc.relabel_checks(need_stats=False)  # the stats file's hash is in each task file, checked by `guards`
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    reading = t64_reading(loaded, committed, checks, bad)
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / OUT_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / OUT_FILES[1]).write_text(t64_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T64")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T64")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T64's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out)
    else:
        {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)


if __name__ == "__main__":
    main()
