#!/usr/bin/env python
"""Try T65 (experiments/harness/TRIES.md): candidate MS of the T1 order list after decision log 2026-11-20 (r), read on
row T62's five practice tasks. MS is R13's file (set t1rc, `t1 --g-mom --base lib3 --feng-groups recomputed`) with one
change: the drawn cells whose label has a recomputed group (CM, Endocardium, Endothelium, Epicardium) also get B0's
per-category damped mean momentum, applied before the atlas step and g_mom (arm 'all'). Cells without a group keep R1's
hold. The facet is expression only.

The practice analogue of R13's file is row T62's `m1rc_gmom`, taken unchanged from `recomputed_read.cut_pair` and
`recomputed_read.library_pair` (their CANDIDATE; their other arms are built and not scored). It is the gate and the
base, not counted. Two counted arms, each that base with one change on the covered cells (the cells whose label's
group has a recomputed direction, build_submission's `direction_groups.get(label) in directions` test):
  m1rc_gmom_mom      B0's per-category damped momentum at harness.B0_DAMPING (0.5) through `t1_scope.momentum_shift`,
                     log f = 0.5 x log((m_later + 0.5) / (m_earlier + 0.5)) per category over the category's own cells,
                     f = 1 on R2's held genes (the direction step's list), keep_depth, stream [seed, 1] over the covered
                     cells alone; a category with fewer than 20 own cells, or drawn from stand-ins at either stage,
                     keeps f = 1. Then the atlas step on the same cells exactly as the base takes it
                     (`t1_scope.external_shift`, the same groups, fraction, held genes and stream [seed, 2]), then g_mom
                     as the base's (`recomputed_read.rc_hook`: the same factors, stream [seed, 5], the same mask).
  m1rc_gmom_momperm  the same, with each category's log factor permuted across its non-held genes within
                     harness.LP_BINS (10) bins of the category's mean CP10k at the fit's later end (its own cells there,
                     the rows `momentum_shift` hands its `logf_filter`), by `harness.t3_decile_permutation` with
                     default_rng(PERM_SEED = 20261110), one generator per category: the size-preserving placebo (each
                     category keeps its non-held values, the held genes stay at f = 1).
The covered cells' draw, the uncovered cells and every stream are the base's: each run checks per seed and arm, before
scoring, that every cell whose label is not covered has the base's bytes, and stops otherwise.
  Cuts: one `t1_scope.build_submission` call per arm with the base's draw keywords and a post_cells hook that takes the
  momentum, the atlas step and g_mom (`ms_population`, then `recomputed_read.rc_hook`); the momentum is fitted from the
  task's start (E8.5) to its reference (E9.5) as T61's b0_gmom fits it (start and reference fallbacks of the chamber
  config for their file stems, eps 0.5, min_cells 20).
  Libraries: T60's oracle draw, the same `ms_population` with the momentum fitted from released E8.5 to released E9.5
  with T61's library settings (within surrounding, the E8.5 and E9.5 fallbacks, min_cells 20, eps 0.5). E10.5 enters
  no fit.
It imports recomputed_read.py (T62) and, through it, bets_read.py (T61) and gmom_combo.py (T26). No builder option for
MS is written here; the builder option is written later and only if this read says AHEAD.

  run-cut  One E10.5 cut. Refused before anything loads unless row T65 is in TRIES.md, HEAD is --commit, the output
           file does not exist, the tree is clean by vec.provenance's rule, this script's and recomputed_read.py's
           source guards pass, the relabelling's committed checks pass and the damping is 0.5. Prints no score.
  run-lib  One library task; the same refusals. Prints no score.
  read     T65's statement, applied before anything is printed; writes t65.json and t65.txt; prints the outcome words
           and what follows only.

Usage: through experiments/t1_m1/t65_ms.sh on the DGX campaign checkout ~/v/vec-campaign-d (it alone holds
data/derived/feng_relabel_c/), at the row's commit, committing nothing until the five task files exist.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec.adoption import SEEDS
from vec.boards import get_board, load_boards
from vec.core import mean_cp10k
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

ROW = "T65"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"
OUT_FILES = ("t65.json", "t65.txt")

BASE = rc.CANDIDATE  # m1rc_gmom: the gate and the base, not counted
MOM = "m1rc_gmom_mom"
MOMPERM = "m1rc_gmom_momperm"
ARMS = (BASE, MOM, MOMPERM)
GATED = (BASE,)
COUNTED = (MOM, MOMPERM)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first (it sets the anchors)
LIB_ARMS = (LIB_REFERENCE, *ARMS)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 10
FIRST_TRY = 418  # T64 ended at T1 418
TRY_RANGE = (FIRST_TRY, FIRST_TRY + TRIES_COUNT)  # T1 418 to 428

DAMPING = harness.B0_DAMPING  # 0.5, checked before anything loads
MOMENTUM_STREAM = 1  # build_submission's momentum stream [seed, 1]
DIRECTION_STREAM = 2  # build_submission's direction stream [seed, 2]
PERM_SEED = 20261110  # the placebo's one generator per category
PERM_BINS = harness.LP_BINS  # 10
REPRO_TOL = br.REPRO_TOL  # 0.02 board points: m1rc_gmom against row T62's committed m1rc_gmom, per task and seed
FRACTION_TOL = 1e-12
AHEAD_BAR = br.AHEAD_BAR  # 0.3
MIN_AHEAD = br.MIN_AHEAD  # 8 of 10
WHOLE_ABOVE = 0.0  # T61's clause: the whole cut's mean strictly above 0
LIB_FLOOR = -0.3
MIN_LIBRARIES = br.MIN_LIBRARIES  # 2
MIN_BEATS = 8  # mom above momperm on the chamber cut in at least 8 of 10 seeds
PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2
BANNED_IMPORTS = ("perm_read", "eb_read", "step_read", "deep_read", "over_read", "lib_step_read", "ctl_range_read",
                  "n48_m1_direction_shrink", "n12_hygiene_size", "double_step_read", "heldmom_read")
BANNED_PREFIX = "n5" + "4"  # no script of row N54 (spelled in pieces so that recomputed_read's guard passes this file)

AHEAD_WORDS = {True: "MS is AHEAD", False: "MS is NOT AHEAD"}
BESIDE = ("MS changes M1's design as ruled (the direction in place of momentum on the direction-group cells, decision "
          "log 2026-10-19 (c)): only the captain's ruling naming the arm puts it in a file.",
          "J6: the momentum's fit spans one day on every practice task and the test's two (horizon ratio 2), so any "
          "gain is reported as 0 to the practice value.")
FOLLOWS_AHEAD = ("What follows: one priced question to the captain. Only on his yes is the builder option written, "
                 "reviewed and dry-run traced, and the file built twice with one sha256 in a campaign checkout where "
                 "the default command first rebuilds ad8f24c7. Nothing is built or changed by this read.")
FOLLOWS_END = ("What follows: candidate MS ends. No question goes to the captain, nothing is built, and nothing is "
               "rerun under another row.")
LABEL = ("try T65, a labelled in-sample diagnostic, never an adoption result: candidate MS of the T1 order list after "
         "decision log 2026-11-20 (r), R13's file with B0's released per-category momentum added on the "
         "direction-group cells before the atlas step (arm 'all'), by its practice analogue m1rc_gmom_mom and the "
         "permuted placebo m1rc_gmom_momperm against m1rc_gmom (row T62's) on T62's five tasks; every task scores "
         "one truth (E10.5) over one day, the test spans two")


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
    for token in ("import" + "lib", "__imp" + "ort__", "spec_from" + "_file_location", "T1_s2" + ".h5ad",
                  "index" + ".json"):
        if token in text:
            bad.append(f"the source names {token}")
    return bad


# ---------------------------------------------------------------------------------------------------------------------
# The placebo's factor


def perm_filter(seed: int = PERM_SEED, bins: int = PERM_BINS) -> tuple[Callable, list[dict]]:
    """momentum_shift's `logf_filter` for the placebo, and the list its records go to (one per category whose factor
    is computed, in momentum_shift's order). For one category: logf (damped, held genes already 0), counts (the
    category's own cells at the later end, in gene order) and held (bool by gene). The non-held genes are ranked by the
    category's mean CP10k there and their log factors permuted within `bins` runs of that ranking by
    `harness.t3_decile_permutation` with default_rng(seed); the held genes keep 0."""
    records: list[dict] = []

    def filt(logf, counts, held):
        logf = np.asarray(logf, dtype=np.float64)
        free = np.flatnonzero(~np.asarray(held, dtype=bool))
        basis = mean_cp10k(counts)[free]
        out = np.zeros_like(logf)
        out[free] = harness.t3_decile_permutation(logf[free], basis, seed=seed, bins=bins)
        rec = {"seed": int(seed), "bins": int(bins), "free_genes": int(free.size),
               "sum_abs_before": float(np.abs(logf[free]).sum()), "sum_abs_after": float(np.abs(out[free]).sum()),
               "same_values": bool(np.array_equal(np.sort(logf[free]), np.sort(out[free]))),
               "moved_genes": int((out[free] != logf[free]).sum())}
        records.append(rec)
        return out, rec

    return filt, records


# ---------------------------------------------------------------------------------------------------------------------
# The covered cells


def covered_mask(labels, groups: Mapping[str, str], directions: Mapping) -> np.ndarray:
    """build_submission's covered test: a cell whose label's group has a direction."""
    return np.array([groups.get(str(lab)) in directions for lab in labels], dtype=bool)


def ms_population(pop, start, ref, directions: Mapping, fraction: float, groups: Mapping[str, str],
                  held_m1: Sequence[str], seed: int, shift_kw: Mapping, logf_filter: Callable | None = None):
    """The arm's population before g_mom: on the covered cells B0's momentum (`momentum_shift` from `start` to `ref`,
    DAMPING, the held genes `held_m1`, stream [seed, MOMENTUM_STREAM], `shift_kw`, `logf_filter`), then the atlas step
    on the same cells (`external_shift`, stream [seed, DIRECTION_STREAM], as build_submission and moved_by_feng take
    it); every other cell as drawn; stacked by `t1_scope._stack_rows`. Returns (population, momentum report or None)."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cov = covered_mask(labels, groups, directions)
    parts, rep = [], None
    if cov.any():
        idx = np.flatnonzero(cov)
        mom, rep = t1_scope.momentum_shift(start, ref, pop[idx].copy(), DAMPING, list(held_m1),
                                           label_key=t1_scope.LABEL_KEY, seed=[seed, MOMENTUM_STREAM],
                                           logf_filter=logf_filter, **shift_kw)
        moved, _ = t1_scope.external_shift(mom, directions, fraction, groups, list(held_m1),
                                           label_key=t1_scope.LABEL_KEY, seed=[seed, DIRECTION_STREAM])
        parts.append((idx, moved))
    if not cov.all():
        idx = np.flatnonzero(~cov)
        parts.append((idx, pop[idx].copy()))
    return t1_scope._stack_rows(parts, pop.n_obs), rep


def builder_covered(labels, directions: Mapping) -> list[str]:
    """The labels the builder moves under --feng-groups recomputed: every present label but
    `build_test._labels_without_direction` on the groups with a direction (the set R1 holds)."""
    present = sorted({str(x) for x in labels})
    return sorted(set(present) - set(bt._labels_without_direction(present, sorted(directions))))


def prediction_sha(X) -> str:
    """The sha256 of an arm's scored matrix (cells x genes, float32, C order)."""
    return hashlib.sha256(np.ascontiguousarray(dense(X), dtype=np.float32).tobytes()).hexdigest()


def covered_record(pop, groups: Mapping[str, str], directions: Mapping, base_X, arm_X, arm: str,
                   rep: Mapping | None, perm: Sequence[Mapping] | None = None) -> dict:
    """The per-seed record of one arm's covered cells, after the byte check: every row whose label (pop's obs, in row
    order) is not covered has the base's bytes in the arm, or the run stops."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cov = covered_mask(labels, groups, directions)
    a = np.ascontiguousarray(dense(base_X), dtype=np.float32)
    b = np.ascontiguousarray(dense(arm_X), dtype=np.float32)
    if a.shape != b.shape or a.shape[0] != len(labels):
        raise SystemExit(f"{ROW}: {BASE} and {arm} are not the same cells ({a.shape}, {b.shape}, {len(labels)} labels)")
    equal = a[~cov].tobytes() == b[~cov].tobytes()
    if not equal:
        raise SystemExit(f"{ROW}: a cell outside the covered set differs between {BASE} and {arm}; the run stops")
    present = sorted(set(labels))
    return {"covered_set": sorted(lab for lab in present if groups.get(lab) in directions),
            "builder_covered_set": builder_covered(present, directions),
            "covered_cells": int(cov.sum()), "covered_cells_changed": int((a[cov] != b[cov]).any(axis=1).sum()),
            "non_covered_equal": bool(equal), "shift": None if rep is None else list(rep.get("rows") or []),
            "perm": None if perm is None else [dict(x) for x in perm]}


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def cut_arms(ctx: harness.BuildContext, factors, directions: Mapping, fraction: float,
             held_m1: Sequence[str]) -> tuple[dict[str, Stage], dict]:
    """m1rc_gmom (recomputed_read.cut_pair's CANDIDATE, as called there), m1rc_gmom_mom and m1rc_gmom_momperm (one
    build_submission call each with the base's draw keywords and `ms_population` then rc_hook as post_cells) on one
    harness cut and seed, each finished as the harness finishes a candidate, and the covered cells' records."""
    base = rc.cut_pair(ctx, factors, directions, fraction, held_m1)[BASE]
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
    out, recs = {BASE: base}, {}
    for arm in COUNTED:
        filt, perm = perm_filter() if arm == MOMPERM else (None, None)
        seen = []

        def post(pop, filt=filt, seen=seen):
            moved, rep = ms_population(pop, start, ref, directions, fraction, groups, held_m1, ctx.seed, shift_kw,
                                       filt)
            seen.append((pop, rep))
            return hook(moved)

        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src, post_cells=post)
        out[arm] = harness.finish(ctx, ref, dense(sub.X), None)
        pop, rep = seen[0]
        recs[arm] = covered_record(pop, groups, directions, base.X, out[arm].X, arm,
                                   None if rep is None else rep.as_dict(), perm)
    return out, recs


def library_arms(e85, e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                 n_cells: int = rr.N_CELLS) -> tuple[dict, dict]:
    """m1rc_gmom (recomputed_read.library_pair's CANDIDATE, as called there), m1rc_gmom_mom and m1rc_gmom_momperm on
    the same oracle draw, as AnnData, and the covered cells' records."""
    base = rc.library_pair(e95, half_shares, factors, held, held_m1, seed, within, fb, directions, fraction,
                           n_cells=n_cells)[BASE]
    base0, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    groups, covered, with_direction = rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), directions)
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    shift_kw = {"within": within, "start_fallback": fb["E8.5"], "ref_fallback": fb["E9.5"], "min_cells": rr.MIN_CELLS,
                "eps": rr.EPS}
    out, recs = {BASE: base}, {}
    for arm in COUNTED:
        filt, perm = perm_filter() if arm == MOMPERM else (None, None)
        pop, rep = ms_population(base0, e85, e95, directions, fraction, groups, held_m1, seed, shift_kw, filt)
        out[arm] = hook(pop)
        recs[arm] = covered_record(base0, groups, directions, base.X, out[arm].X, arm,
                                   None if rep is None else rep.as_dict(), perm)
    return out, recs


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def refuse(out: Path, name: str | Sequence[str], commit: str | None) -> tuple[dict, dict]:
    """The refusals before anything loads: the row, HEAD not the run commit, an existing output file, an unclean tree,
    both source guards, the relabelling's committed checks and the damping. Returns the start state and the checks."""
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, (name,) if isinstance(name, str) else tuple(name))
    started = gv.tree_state()
    if not commit or started.get("commit") != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit given by --commit")
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


def _task_fields(rows: Sequence[Mapping], held_m1: Sequence[str]) -> dict:
    first = rows[0]["covered"] if rows else {}
    return {"damping": DAMPING, "momentum_stream": MOMENTUM_STREAM, "direction_stream": DIRECTION_STREAM,
            "perm_seed": PERM_SEED, "perm_bins": PERM_BINS,
            "shift_rows": {a: (first.get(a) or {}).get("shift") for a in COUNTED} if rows else None,
            "constants": _constants(held_m1)}


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
        shas = {k: prediction_sha(v.X) for k, v in arms.items()}
        rows.append({"seed": s, "floor": float(f.score), "arms": gv.score_arms(board, target, ref_stage, s, f.anchors,
                                                                               arms),
                     "sha256": shas, "covered": recs, "seconds": round(time.time() - t0, 1)})
        del arms
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "relabel": rc._relabel_fields(ext, checks), **_task_fields(rows, held_m1),
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
        shas = {k: prediction_sha(v.X) for k, v in stages.items()}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s), "sha256": shas,
                     "covered": recs, "seconds": round(time.time() - t0, 1)})
        del arms, stages
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "relabel": rc._relabel_fields(ext, checks), **_task_fields(rows, held_m1),
           "split": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD},
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


def idle(r: Mapping, arm: str) -> bool:
    """The arm and the base are one file on every seed (by prediction sha256)."""
    return all(x["sha256"][arm] == x["sha256"][BASE] for x in r["rows"])


def guards(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
           source_bad: Sequence[str], row_commit: str | None = None, seeds: Sequence[int] = SEEDS) -> list[str]:
    """Every void condition of row T65; an empty list means none fired."""
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
        if row_commit is not None and stamp.get("commit") != row_commit:
            reasons.append(f"{k}: it did not run at the row's commit")
        ref = committed.get(k, {})
        worst = max((abs(br.score_of(x, BASE) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                     for x in rows), default=float("inf"))
        if not worst <= REPRO_TOL:
            reasons.append(f"{k}: the gate failed: {BASE} differs from row T62's committed {BASE} by up to "
                           f"{worst:.3g} (tolerance {REPRO_TOL:g})")
        if not abs(float(r.get("fraction", float("nan"))) - PRACTICE_FRACTION) <= FRACTION_TOL:
            reasons.append(f"{k}: the practice fraction is {r.get('fraction')}, not {PRACTICE_FRACTION}")
        for field, want_v in (("damping", DAMPING), ("momentum_stream", MOMENTUM_STREAM),
                              ("direction_stream", DIRECTION_STREAM), ("perm_seed", PERM_SEED),
                              ("perm_bins", PERM_BINS)):
            if r.get(field) != want_v:
                reasons.append(f"{k}: the {field.replace('_', ' ')} is {r.get(field)}, not {want_v}")
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
                if arm == MOMPERM and not all(p.get("seed") == PERM_SEED and p.get("bins") == PERM_BINS
                                              and p.get("same_values") is True for p in (c.get("perm") or [])):
                    reasons.append(f"{k}: seed {x['seed']}: the placebo's permutation is not the row's")
    if len(commits) > 1:
        reasons.append(f"the task files ran at more than one commit: {sorted(map(str, commits))}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def diff(rows: Sequence[Mapping], change: str, base: str) -> dict:
    """d = change - base per seed (`bets_read.paired`: mean, sd, seeds above 0, range, points by group), with the
    seeds below 0 and both arms' mean scores."""
    p = br.paired(rows, change, base)
    p["behind"] = int(sum(br.score_of(r, change) - br.score_of(r, base) < 0 for r in rows))
    p["mean_scores"] = {a: float(np.mean([br.score_of(r, a) for r in rows])) for a in (base, change)}
    return p


def statement(tasks: Mapping[str, Mapping], valid: Sequence[str], beats: Mapping) -> tuple[bool, list[str]]:
    """Row T65's statement: AHEAD, and the parts that failed (words only, no figure). tasks: d = mom - base per task;
    beats: d = mom - momperm on the chamber cut."""
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
    if not beats["ahead"] >= MIN_BEATS:
        failed.append(f"{MOM} beats {MOMPERM} on the chamber cut in fewer than {MIN_BEATS} of 10 seeds")
    if not beats["mean"] > 0:
        failed.append(f"the chamber cut's mean of {MOM} minus {MOMPERM} is not positive")
    return not failed, failed


def covered_summary(r: Mapping) -> dict:
    """A task's covered cells per arm: cells changed per seed, the momentum's rows per category (the first seed's
    report) and the placebo's permutation records (the first seed's)."""
    keys = ("type", "status", "n_target", "n_start", "n_ref", "start_cells", "ref_cells", "n_up", "n_down",
            "mean_abs_log_f", "max_abs_log_f", "median_depth_scale")
    out = {}
    for arm in COUNTED:
        first = r["rows"][0]["covered"][arm]
        out[arm] = {"covered_set": first.get("covered_set"),
                    "covered_cells_changed": [int(x["covered"][arm]["covered_cells_changed"]) for x in r["rows"]],
                    "momentum_rows": [{k: v for k, v in row.items() if k in keys}
                                      for row in (first.get("shift") or [])],
                    "perm": first.get("perm")}
    return out


def t65_reading(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]], checks: Mapping,
                source_bad: Sequence[str] = (), row_commit: str | None = None) -> dict:
    """Row T65's reading: the void conditions (the gate among them) first, then the statement."""
    reasons = guards(loaded, committed, checks, source_bad, row_commit)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons, "follows": FOLLOWS_END}
    tasks = {k: diff(r["rows"], MOM, BASE) for k, r in loaded.items()}
    placebo = {k: diff(r["rows"], MOMPERM, BASE) for k, r in loaded.items()}
    versus = {k: diff(r["rows"], MOM, MOMPERM) for k, r in loaded.items()}
    valid = [lib for lib in LIBRARIES if library_valid(loaded[lib])]
    ahead, failed = statement(tasks, valid, versus["chamber"])
    return {"outcome": AHEAD_WORDS[ahead], "ahead": ahead, "failed": failed, "valid_libraries": valid,
            "tasks": tasks, "placebo": placebo, "versus": versus,
            "covered": {k: covered_summary(r) for k, r in loaded.items()},
            "idle_tasks": {arm: [k for k, r in loaded.items() if idle(r, arm)] for arm in COUNTED},
            "j6_range": sorted([0.0, tasks["chamber"]["mean"]]), "beside": list(BESIDE),
            "follows": FOLLOWS_AHEAD if ahead else FOLLOWS_END}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]]
    first = f"Outcome: {r['outcome']}" + ("" if r["ahead"] else " (" + "; ".join(r["failed"]) + ")")
    lines = [first, *r["beside"]]
    lines += [f"{arm} moved no cell on {k}: nothing of it is priced there"
              for arm, ks in r["idle_tasks"].items() for k in ks]
    return lines + [r["follows"]]


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    b, c = p["base"], p["change"]
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} behind {p['behind']} of "
            f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}; mean scores "
            f"{b} {p['mean_scores'][b]:.2f}, {c} {p['mean_scores'][c]:.2f}")


def t65_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        for title, key in ((f"d = {MOM} - {BASE} (the statement's d), paired per seed on the same draw:", "tasks"),
                           (f"d = {MOMPERM} - {BASE} (reported, deciding nothing):", "placebo"),
                           (f"d = {MOM} - {MOMPERM} (the statement reads the chamber cut only):", "versus")):
            lines.append(title)
            lines += [_line(k, p) for k, p in r[key].items()]
        lines.append("Covered cells (reported, deciding nothing):")
        for k, per in r["covered"].items():
            for arm, h in per.items():
                lines.append(f"  {k} {arm}: covered labels {h['covered_set']}; covered cells changed per seed "
                             f"{h['covered_cells_changed']}")
                for m in h["momentum_rows"]:
                    lines.append(f"    momentum {m.get('type')}: {m.get('status')}; own cells {m.get('n_start')} at "
                                 f"the start, {m.get('n_ref')} at the reference ({m.get('start_cells')}, "
                                 f"{m.get('ref_cells')}); mean |log f| {m.get('mean_abs_log_f', 0.0):.3f}")
                for p in h["perm"] or []:
                    lines.append(f"    permutation: seed {p['seed']}, {p['bins']} bins, {p['moved_genes']} of "
                                 f"{p['free_genes']} genes moved, sum |log f| {p['sum_abs_before']:.3f} before and "
                                 f"{p['sum_abs_after']:.3f} after")
        lines.append(f"Valid libraries: {', '.join(r['valid_libraries']) or 'none'}")
        lines.append(f"J6 range (0 to the chamber cut's mean of {MOM} - {BASE}; a negative mean is reported as "
                     f"[mean, 0]): {r['j6_range'][0]:+.2f} to {r['j6_range'][1]:+.2f}")
    lines.append(f"Tries counted: {res['tries']} (T1 {TRY_RANGE[0]} to {TRY_RANGE[1]})")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, row_commit: str | None) -> dict:
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md")
    gv.refuse_existing(out, OUT_FILES)
    if not row_commit:
        sys.exit(f"{ROW} is refused: the read needs the row's commit (--commit)")
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    paths = {k: out / "tasks" / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    committed = {k: committed_t62(f) for k, f in task_files().items()}
    checks = rc.relabel_checks(need_stats=False)  # the stats file's hash is in each task file, checked by `guards`
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    reading = t65_reading(loaded, committed, checks, bad, row_commit)
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT, "row_commit": row_commit}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / OUT_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / OUT_FILES[1]).write_text(t65_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T65")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True, help="the run commit (full hash): HEAD must be it")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T65")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True, help="the run commit (full hash): HEAD must be it")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T65's reading")
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True, help="the row's commit: every task file's stamp must be at it")
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out, a.commit)
    else:
        {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)


if __name__ == "__main__":
    main()
