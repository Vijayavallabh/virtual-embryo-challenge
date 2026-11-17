#!/usr/bin/env python
"""Try T52 (experiments/harness/TRIES.md): M1's Feng step with a smaller pseudocount in the direction, against the same
step inflated uniformly. A labelled in-sample read on tries T44's and T46's routes (experiments/t1_m1/step_read.py,
experiments/t1_m1/lib_step_read.py).

M1 multiplies each covered cell's counts by exp(fraction x d), d per group the log ratio of mean CP10k between Feng's
E9.5 and E14.5 hearts with a pseudocount of 0.5 CP10k (dirlib.lfc_ratio; experiments/t1_direction/config.yaml,
metrics.pseudocount_cp10k). Tries T44 to T46 found the step inflated uniformly (1.5 to 4 times) scores better on every
released truth of the E10.5 cuts. This read puts the other way of delivering more to the same cuts: the direction made
with a smaller pseudocount (0.1 and 0.02), at the entry's fraction and at 1.5 and 2 times it, beside the uniformly
inflated steps already read.

  run    One truth, 10 seeds. Truths: `pooled`, `lib1`, `lib2`, `lib3` (the chamber cut T1/E9.5->E10.5:chamber exactly
         as try T46's run makes them, lib_step_read.truth_and_contamination) and `whole` (the cut T1/E9.5->E10.5 with
         its target exactly as try T44's run loads it). Arms per seed, all from one cell draw: b0, m1_r_gmom (fraction
         0.2) and the uniform arms m1_r_gmom_x15, _x2, _x3 (step_read.step_arms), with _x4 and try T46's placebo
         m1_r_gmom_perm on the four chamber truths (m1_r_gmom_x3 on `whole` has no committed score and is not gated);
         six new arms, the build of m1_r_gmom with Feng's directions made with another pseudocount and nothing else
         changed (m1_r_gmom_e01, _e01_x15, _e01_x2 at pseudocount 0.1 and fractions 0.2, 0.3, 0.4; m1_r_gmom_e002,
         _e002_x15, _e002_x2 at 0.02); on the four chamber truths two placebos, m1_r_gmom_e01_x15_perm and
         m1_r_gmom_e002_x15_perm, the builds of _e01_x15 and _e002_x15 along try T46's permuted directions made from
         that arm's own directions (lib_step_read.placebo_directions: the same two strata, the stream [build seed, 46]).
         The directions are made by feng_dirs.feng_directions with a deep copy of the loaded config whose
         metrics.pseudocount_cp10k is the arm's; no file and no shared object changes. The task file holds what try
         T46's holds (every arm's score, sha256 of the predicted matrix, scored zero share and detected scored genes per
         cell, per group the mean absolute log factor on the scored genes, the truth's row and contamination counts as
         checks never printed), every arm's fraction and pseudocount, each direction set's per-group sha256, per seed
         each placebo's permuted direction's sha256 per group and whether a new placebo's permutation of positions is
         m1_r_gmom_perm's, and per direction set and group the squared size (sum of d squared over the moved genes) over
         the entry's and the share of it on genes whose Feng E9.5 mean is under 0.5 CP10k (ratios). Refused on an
         existing task file, unless the committed N22 passed, unless experiments/harness/TRIES.md holds row T52, and
         unless the loaded config's pseudocount is the entry's 0.5.
  read   T52's reading of the five truths' files (rule below, fixed here and in the row before the run).

The reading. Void unless the five truths hold the ten seeds at one commit on a clean tree (tree_clean exactly true); no
file records an anchor failure; every arm the truth should hold is on every row with its fraction and pseudocount; on
the four chamber truths b0, m1_r_gmom, _x15, _x2, _x3, _x4 and m1_r_gmom_perm equal, seed for seed, their committed
scores in experiments/t1_m1/results/t46_m1_lib_step/tasks/, and on `whole` b0, m1_r_gmom, _x15 and _x2 theirs in
experiments/t1_m1/results/t44_m1_step/tasks/ (the whole cut's file), each to 1e-9; no new arm scores exactly as
m1_r_gmom on every seed of a truth; no placebo's recorded prediction sha256 equals its real arm's on every seed; the new
arms' recorded directions differ from the entry's in every group. Also void, so that the checks above cannot pass
vacuously: a row without a sha256 string for every arm, or without the sparsity and log-factor records.
Tests. On one truth: try T46's single-task test (vec.adoption.keep_change on that truth alone: "ahead" for a mean
paired gain above 0.3 with at least 8 of 10 seeds ahead, "behind" the same reversed, else "level"). On the two cuts
(`pooled` and `whole`): the adoption rule's test (vec.adoption.keep_change, both cuts: "passes"), with the reversed sign
for "loses", else "neither".
Statement 1, per new arm: against m1_r_gmom on the two cuts (passes, loses, neither) and its label on each library
truth.
Statement 2, decompression against uniform inflation at about the same delivered size (read back on library 3's cells
with the entry's ruler: _e01 0.173 against _x15 0.175; _e01_x2 0.347 against _x3 0.349): m1_r_gmom_e01 against
m1_r_gmom_x15 and m1_r_gmom_e01_x2 against m1_r_gmom_x3, each on the two cuts and on each library truth. "the smaller
pseudocount scores better at the same delivered size" if both comparisons pass on the two cuts; "uniform inflation
scores better" if both lose; otherwise "no difference at the bar", with both comparisons' numbers. The 0.02 arms have
no uniform arm of the same delivered size: each is reported against the two uniform arms that bracket it and carries no
statement (_x4 is not built on `whole`, so _e002_x2 against _x4 is reported on the chamber truths only).
Statement 3, per chamber truth and NN in 01, 002: D = (_eNN_x15 - m1_r_gmom) - (_eNN_x15_perm - m1_r_gmom_perm) per
seed with the single-task test: "the gain needs the direction's genes" where D is ahead and _eNN_x15 is ahead of
m1_r_gmom; "no gain on this truth" where _eNN_x15 is not ahead of m1_r_gmom; otherwise "not separated".
Reported: per truth every arm's mean board score and the arm with the highest mean among m1_r_gmom, the uniform arms and
the new arms (placebos and b0 left out; the first in the file's arm order on a tie); each new arm's paired difference
from m1_r_gmom and from the uniform arm of the nearest delivered size present on that truth, with seeds ahead, per-seed
range and sd; every metric's mean skill difference for the comparisons of statements 1 and 2 and the part on the two de
terms and on mmd_u and the variogram; each arm's scored zero share and detected scored genes per cell. The read puts
nothing forward by itself.

Limits, repeated in the text (and in the row):
- One released day, in-sample, on cuts earlier choices were made on.
- A smaller pseudocount lets weakly expressed genes' ratios, which rest on few counts, move cells further, so a gain can
  be the larger factors on such genes and a loss their noise: the placebo and the uniform arms bound that and do not
  remove it.
- The delivered sizes quoted are read-backs on E10.5 library 3's cells, not on the cuts' own reference cells, and equal
  delivered size along the entry's direction does not make two arms equal in anything else (the smaller pseudocount's
  direction is 3 to 7 times larger in squared size).
- The pseudocounts 0.1 and 0.02 were chosen before the run from the read-back table and no other value is read.
- Nothing here says how far any hidden truth moves.

NOT an adoption result: Feng's E9.5 cells lie inside these pseudo-tasks' external windows and the cuts are ones earlier
T1 choices were made on. The result can only go to the team captain.

Usage: on the DGX campaign checkout, repo root, all five truths at one commit, committing nothing until the five task
files exist, and on the machine T44 and T46 ran on, so that the rescored arms reproduce. The environment variables are
docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eps_read.py run --truth pooled --out t52_m1_eps
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eps_read.py run --truth lib1 --out t52_m1_eps
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eps_read.py run --truth lib2 --out t52_m1_eps
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eps_read.py run --truth lib3 --out t52_m1_eps
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eps_read.py run --truth whole --out t52_m1_eps
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/eps_read.py read --out t52_m1_eps
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS, keep_change
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import eb_read as er  # noqa: E402
import feng_dirs  # noqa: E402
import gmom_combo as gc  # noqa: E402
import lib_step_read as lsr  # noqa: E402
import step_read as sr  # noqa: E402

gv, hygiene, n12 = gc.gv, gc.hygiene, er.n12
dl = feng_dirs.dl
ROW = "T52"
CHAMBER_TRUTHS = lsr.TRUTHS  # pooled, lib1, lib2, lib3: the chamber cut as try T46 makes them
LIBRARIES = tuple(lsr.LIBRARIES)
WHOLE = "whole"  # the cut T1/E9.5->E10.5 as try T44 loads it
TRUTHS = (*CHAMBER_TRUTHS, WHOLE)
CUTS = (lsr.POOLED, WHOLE)  # the two cuts of the adoption rule's test
TASK_OF = {**{k: gc.TASK_CHAMBER for k in CHAMBER_TRUTHS}, WHOLE: gc.TASK_WHOLE}
BASE = sr.BASE  # m1_r_gmom
FRACTION, FRACTION_TOL, REPRO_TOL = sr.FRACTION, sr.FRACTION_TOL, sr.REPRO_TOL  # 0.2, 1e-12, 1e-9
ENTRY_EPS = 0.5  # the entry's pseudocount (config metrics.pseudocount_cp10k); the run refuses another
EPS = {"e01": 0.1, "e002": 0.02}  # chosen before the run from the read-back table; no other value is read
EPS_TOL = 1e-12
UNIFORM = {"m1_r_gmom_x15": 1.5, "m1_r_gmom_x2": 2.0, "m1_r_gmom_x3": 3.0, "m1_r_gmom_x4": 4.0}
UNIFORM_WHOLE = ("m1_r_gmom_x15", "m1_r_gmom_x2", "m1_r_gmom_x3")  # no _x4 on the whole cut
STEPS = {"": 1.0, "_x15": 1.5, "_x2": 2.0}  # a new arm's multiplier of the entry's fraction
PERM_ENTRY = lsr.PLACEBO_BASE  # m1_r_gmom_perm, try T46's placebo at the base fraction
PLACEBO_STEP = "_x15"  # the new placebos are the _x15 builds along permuted directions
PLACEBO_STREAM = lsr.PLACEBO_STREAM  # 46


def new_arm(tag: str, step: str) -> str:
    return f"{BASE}_{tag}{step}"


def new_placebo(tag: str) -> str:
    return f"{new_arm(tag, PLACEBO_STEP)}_perm"


NEW_ARMS = {new_arm(t, s): (t, m) for t in EPS for s, m in STEPS.items()}  # arm -> (pseudocount tag, multiplier)
NEW_PLACEBOS = {new_placebo(t): new_arm(t, PLACEBO_STEP) for t in EPS}  # placebo -> its real arm
PLACEBO_OF = {PERM_ENTRY: BASE, **NEW_PLACEBOS}
REPRO_CHAMBER = ("b0", BASE, *UNIFORM, PERM_ENTRY)  # against try T46's committed files
REPRO_WHOLE = ("b0", BASE, "m1_r_gmom_x15", "m1_r_gmom_x2")  # against try T44's committed whole-cut file
T46_TASKS = RESULTS / "t46_m1_lib_step" / "tasks"
T44_WHOLE = RESULTS / "t44_m1_step" / "tasks" / gc.CUT_FILES[gc.TASK_WHOLE]
# delivered size read back on E10.5 library 3's cells with the entry's ruler (the 2026-11-14 and 2026-11-17 target-free
# checks, as the brief quotes them); _e01_x15 is not quoted and is 1.5 x _e01's (the read-back is linear in the
# fraction)
DELIVERED = {"m1_r_gmom_x15": 0.175, "m1_r_gmom_x2": 0.233, "m1_r_gmom_x3": 0.349, "m1_r_gmom_x4": 0.463,
             "m1_r_gmom_e01": 0.173, "m1_r_gmom_e01_x15": 1.5 * 0.173, "m1_r_gmom_e01_x2": 0.347,
             "m1_r_gmom_e002": 0.207, "m1_r_gmom_e002_x15": 0.312, "m1_r_gmom_e002_x2": 0.416}
MATCHED = (("m1_r_gmom_e01", "m1_r_gmom_x15"), ("m1_r_gmom_e01_x2", "m1_r_gmom_x3"))  # statement 2's two comparisons
BRACKETS = {"m1_r_gmom_e002": ("m1_r_gmom_x15", "m1_r_gmom_x2"), "m1_r_gmom_e002_x15": ("m1_r_gmom_x2", "m1_r_gmom_x3"),
            "m1_r_gmom_e002_x2": ("m1_r_gmom_x3", "m1_r_gmom_x4")}
LOW_CP10K = 0.5  # the share of the squared size on genes whose Feng E9.5 mean is under this (the entry's pseudocount)
PASSES, LOSES, NEITHER = "passes", "loses", "neither"
LEVEL = {"ahead": "ahead", "behind": "behind", "neither": "level"}
SMALLER_BETTER = "the smaller pseudocount scores better at the same delivered size"
UNIFORM_BETTER = "uniform inflation scores better"
NO_DIFFERENCE = "no difference at the bar"
NEEDS_GENES = "the gain needs the direction's genes"
NO_GAIN = "no gain on this truth"
NOT_SEPARATED = "not separated"
TASK_FILE = "{truth}.json"
LABEL_T52 = ("try T52, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with Feng's directions made with a pseudocount of 0.1 and 0.02 CP10k instead of 0.5, "
             "at 1, 1.5 and 2 times the entry's fraction, beside the step inflated uniformly, on tries T44's and T46's "
             "routes (the chamber cut pooled and library by library, and the whole cut); Feng's E9.5 cells lie inside "
             "these pseudo-tasks' external windows and the cuts are ones earlier T1 choices were made on; the result "
             "can only go to the team captain")
LIMITS = (
    "One released day, in-sample, on cuts earlier choices were made on.",
    "A smaller pseudocount lets weakly expressed genes' ratios, which rest on few counts, move cells further, so a "
    "gain can be the larger factors on such genes and a loss their noise: the placebo and the uniform arms bound that "
    "and do not remove it.",
    "The delivered sizes quoted are read-backs on E10.5 library 3's cells, not on the cuts' own reference cells, and "
    "equal delivered size along the entry's direction does not make two arms equal in anything else (the smaller "
    "pseudocount's direction is 3 to 7 times larger in squared size).",
    "The pseudocounts 0.1 and 0.02 were chosen before the run from the read-back table and no other value is read.",
    "Nothing here says how far any hidden truth moves.",
)


def uniform_of(truth: str) -> dict[str, float]:
    return {a: m for a, m in UNIFORM.items() if truth != WHOLE or a in UNIFORM_WHOLE}


def arms_of(truth: str) -> tuple[str, ...]:
    """The arms a truth's file holds, in the order the run builds them."""
    chamber = truth != WHOLE
    return ("b0", BASE, *uniform_of(truth), *((PERM_ENTRY,) if chamber else ()), *NEW_ARMS,
            *(NEW_PLACEBOS if chamber else ()))


def arm_spec(truth: str, fraction: float = FRACTION) -> dict[str, dict[str, float]]:
    """Every directional arm's fraction and pseudocount on this truth (b0 has neither)."""
    out = {BASE: {"fraction": float(fraction), "pseudocount": ENTRY_EPS}}
    out.update({a: {"fraction": float(fraction) * m, "pseudocount": ENTRY_EPS} for a, m in uniform_of(truth).items()})
    if truth != WHOLE:
        out[PERM_ENTRY] = {"fraction": float(fraction), "pseudocount": ENTRY_EPS}
    out.update({a: {"fraction": float(fraction) * m, "pseudocount": EPS[t]} for a, (t, m) in NEW_ARMS.items()})
    if truth != WHOLE:
        out.update({p: {"fraction": float(fraction) * STEPS[PLACEBO_STEP], "pseudocount": EPS[NEW_ARMS[a][0]]}
                    for p, a in NEW_PLACEBOS.items()})
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The directions


def eps_config(eps: float, cfg: Mapping | None = None) -> dict:
    """A deep copy of the direction config (the loaded file unless given) with metrics.pseudocount_cp10k = eps; the
    original object and the file are untouched."""
    out = copy.deepcopy(dict(feng_dirs.load_config() if cfg is None else cfg))
    out["metrics"]["pseudocount_cp10k"] = float(eps)
    return out


def eps_directions(genes: Sequence[str], eps: float, stats_path: Path | None = None) -> tuple[dict, dict]:
    """feng_dirs.feng_directions on `genes` with a copy of the config whose pseudocount is eps."""
    return feng_dirs.feng_directions(list(genes), cfg=eps_config(eps), stats_path=stats_path)


def vector_sha256(v) -> str:
    return hashlib.sha256(np.ascontiguousarray(np.asarray(v, dtype=np.float64)).tobytes()).hexdigest()


def directions_sha256(directions: Mapping[str, np.ndarray]) -> dict[str, str]:
    return {g: vector_sha256(directions[g]) for g in sorted(directions)}


def start_means(genes: Sequence[str], record: Mapping, stats_path: Path | None = None) -> dict[str, np.ndarray]:
    """Per group of a direction record, the group's Feng mean CP10k at the bracket's start (E9.5) on `genes` (NaN where
    unmapped), from the same statistics file and labels the directions used."""
    path = Path(stats_path) if stats_path is not None else feng_dirs.clean_stats_path()
    feng_dirs.fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    idx = dl.sym_index(st["genes"], list(genes))
    start = record["bracket"][0]
    return {g: dl.on_t1(dl.feng_sel(st, start, c["feng_labels"])["cp"], idx) for g, c in record["groups"].items()}


def moved_mask(v, genes: Sequence[str], held: Sequence[str]) -> np.ndarray:
    """The genes the step moves: a finite direction value and not held (external_shift's factor is 1 elsewhere)."""
    held_set = set(map(str, held))
    return np.isfinite(np.asarray(v, dtype=np.float64)) & np.array([str(g) not in held_set for g in genes])


def direction_stats(sets: Mapping[str, Mapping[str, np.ndarray]], means: Mapping[str, np.ndarray],
                    genes: Sequence[str], held: Sequence[str]) -> dict:
    """Per direction set (`entry` and the new ones) and group: the squared size (sum of d squared over the moved genes)
    over the entry's, and the share of it on genes whose Feng start mean is under LOW_CP10K. Ratios only."""
    entry = sets["entry"]
    out = {}
    for tag, dirs in sets.items():
        rec = {}
        for g in sorted(dirs):
            v = np.asarray(dirs[g], dtype=np.float64)
            m = moved_mask(v, genes, held)
            sq = float((v[m] ** 2).sum())
            e = np.asarray(entry[g], dtype=np.float64)
            sq_entry = float((e[moved_mask(e, genes, held)] ** 2).sum())
            low = m & (np.asarray(means[g]) < LOW_CP10K)
            rec[g] = {"squared_size_over_entry": sq / sq_entry if sq_entry > 0 else None,
                      "share_on_low_start_genes": float((v[low] ** 2).sum()) / sq if sq > 0 else None}
        out[tag] = rec
    return out


def positions_of(directions: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Each group's vector of its own positions, NaN where the direction has none: permuted as a direction is, it
    shows the permutation of positions itself."""
    return {g: np.where(np.isfinite(np.asarray(v, dtype=np.float64)), np.arange(len(v), dtype=np.float64), np.nan)
            for g, v in directions.items()}


def same_positions(ctx: harness.BuildContext, a: Mapping, b: Mapping, held: Sequence[str]) -> bool:
    """Whether try T46's placebo permutes the positions of directions a and b the same way for this seed."""
    pa = lsr.placebo_directions(ctx, positions_of(a), held)
    pb = lsr.placebo_directions(ctx, positions_of(b), held)
    return set(pa) == set(pb) and all(np.array_equal(pa[g], pb[g], equal_nan=True) for g in pa)


def mean_abs_log_factor(v, genes: Sequence[str], held: Sequence[str], scored: Sequence[str], fraction: float) -> float:
    """The mean absolute log factor over the scored genes at `fraction` (lib_step_read.scored_log_factors' measure)."""
    return lsr.scored_log_factors({"g": v}, {"g": v}, genes, held, scored, fraction)["g"][BASE]


# ---------------------------------------------------------------------------------------------------------------------
# The run


def build_seed(ctx: harness.BuildContext, b0, truth: str, directions: Mapping, eps_dirs: Mapping[str, Mapping],
               groups: Mapping, fraction: float, held_m1: Sequence[str],
               held_factors: Sequence[str]) -> tuple[dict[str, Stage], dict, dict]:
    """One seed's arms in arms_of(truth)'s order; per group and directional arm the mean absolute log factor on the
    scored genes; and the placebos' record (each permuted direction's sha256 per group, and whether a new placebo's
    permutation of positions is m1_r_gmom_perm's)."""
    chamber = truth != WHOLE
    preds = {"b0": b0(ctx), **sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors,
                                           multipliers=uniform_of(truth))}
    used = {a: (directions, f) for a, f in sr.arm_fractions(fraction, uniform_of(truth)).items()}
    perms = {}
    if chamber:
        perms[PERM_ENTRY] = lsr.placebo_directions(ctx, directions, held_m1)
        preds[PERM_ENTRY] = gc.m1_r_arms(ctx, perms[PERM_ENTRY], groups, fraction, held_m1, held_factors)[BASE]
        used[PERM_ENTRY] = (perms[PERM_ENTRY], float(fraction))
    for tag, dirs in eps_dirs.items():
        steps = {new_arm(tag, s): m for s, m in STEPS.items() if s}
        built = sr.step_arms(ctx, dirs, groups, fraction, held_m1, held_factors, multipliers=steps)
        fr = sr.arm_fractions(fraction, steps)
        for a, p in built.items():
            name = new_arm(tag, "") if a == BASE else a
            preds[name] = p
            used[name] = (dirs, fr[a])
    placebo = {"sha256": {}, "same_positions_as_m1_r_gmom_perm": {}}
    if chamber:
        placebo["sha256"][PERM_ENTRY] = directions_sha256(perms[PERM_ENTRY])
        for tag, dirs in eps_dirs.items():
            p = new_placebo(tag)
            perms[p] = lsr.placebo_directions(ctx, dirs, held_m1)
            f = float(fraction) * STEPS[PLACEBO_STEP]
            preds[p] = gc.m1_r_arms(ctx, perms[p], groups, f, held_m1, held_factors)[BASE]
            used[p] = (perms[p], f)
            placebo["sha256"][p] = directions_sha256(perms[p])
            placebo["same_positions_as_m1_r_gmom_perm"][p] = same_positions(ctx, directions, dirs, held_m1)
    genes = list(ctx.inputs[ctx.task.reference_role].var_names)
    scored = list(ctx.board.genes)
    log_factors = {g: {a: mean_abs_log_factor(d[g], genes, held_m1, scored, f) for a, (d, f) in used.items()}
                   for g in sorted(directions)}
    return {a: preds[a] for a in arms_of(truth)}, log_factors, placebo


def truth_and_contamination(t: harness.PseudoTask, data: Path, genes, truth: str) -> tuple[Stage, int | None]:
    """A chamber truth as try T46's run makes it; `whole` as try T44's run loads its target."""
    if truth != WHOLE:
        return lsr.truth_and_contamination(t, data, genes, truth)
    if t.target_filter is not None:
        target, rec = harness.filter_target(t, data, genes)
        return target, int(rec["contamination"]) if "contamination" in rec else None
    return Stage.from_anndata(read_stage(data / t.target, genes), genes), None


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T52 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T52 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    eps_cfg = feng_dirs.load_config()["metrics"]["pseudocount_cp10k"]
    if float(eps_cfg) != ENTRY_EPS:
        sys.exit(f"T52 is refused: the direction config's pseudocount is {eps_cfg}, not the entry's {ENTRY_EPS}")
    out = RESULTS / a.out
    name = TASK_FILE.format(truth=a.truth)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK_OF[a.truth]]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    truth, contamination = truth_and_contamination(t, data, genes, a.truth)
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T52, "task": t.key, "truth": a.truth, "truth_rows": int(truth.n),
           "truth_contamination": contamination, "seeds": list(seeds)}
    floors, failures = lsr.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        (out / "tasks").mkdir(parents=True, exist_ok=True)
        (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")
        sys.exit(f"T52 is refused on {a.truth}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T52 is void)")
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    pool_genes = list(pool.var_names)
    directions, ext = feng_dirs.feng_directions(pool_genes)
    eps_sets = {tag: eps_directions(pool_genes, eps) for tag, eps in EPS.items()}
    eps_dirs = {tag: d for tag, (d, _) in eps_sets.items()}
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    means = start_means(pool_genes, ext)  # on the directions' own gene order
    sets = {"entry": directions, **eps_dirs}
    dir_record = {"entry": {"pseudocount": ENTRY_EPS, "sha256": directions_sha256(directions)},
                  **{tag: {"pseudocount": EPS[tag], "sha256": directions_sha256(d)} for tag, d in eps_dirs.items()}}
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds, log_factors, placebo = build_seed(ctx, b0, a.truth, directions, eps_dirs, groups, fraction, held_m1,
                                                 preregistered)
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms,
                     "sha256": {k: lsr.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: lsr.sparsity(p) for k, p in preds.items()}, "scored_log_factors": log_factors,
                     "placebo": placebo, "seconds": round(time.time() - t0, 1)})
        print(f"{a.truth} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    spec = arm_spec(a.truth, fraction)
    res.update({"fraction": fraction, "fractions": {k: v["fraction"] for k, v in spec.items()},
                "pseudocounts": {k: v["pseudocount"] for k, v in spec.items()}, "directions": dir_record,
                "direction_stats": direction_stats(sets, means, pool_genes, held_m1),
                "placebo": {"arms": {p: r for p, r in PLACEBO_OF.items() if p in spec}, "stream": PLACEBO_STREAM,
                            "strata": ["scored genes", "other moved genes"]},
                "external": {"entry": ext, **{tag: rec for tag, (_, rec) in eps_sets.items()}}, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), feng_dirs.CONFIG,
                                            n22])})
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(truth: str, rows: Sequence[Mapping], committed_rows: Sequence[Mapping]) -> dict:
    """Per seed, each rescored arm against its committed score; holds when every seed and arm is there and each is
    within REPRO_TOL."""
    arms = REPRO_WHOLE if truth == WHOLE else REPRO_CHAMBER
    by_seed = {r["seed"]: r for r in committed_rows}
    worst, absent = {a: 0.0 for a in arms}, []
    for r in rows:
        old = by_seed.get(r["seed"])
        for arm in arms:
            if old is None or arm not in old["arms"] or arm not in r["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - old["arms"][arm]["score"]))
    holds = bool(rows) and bool(committed_rows) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "absent": absent, "holds": holds}


def _off(got, want: float, tol: float) -> bool:
    return not isinstance(got, int | float) or isinstance(got, bool) or abs(got - want) > tol


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    """The void conditions. results: truth -> task file (a missing truth is absent); committed: truth -> the committed
    rows its rescored arms are checked against (try T46's file per chamber truth, try T44's whole-cut file)."""
    reasons = []
    missing = [k for k in TRUTHS if k not in results]
    if missing:
        reasons.append(f"missing truths {missing}")
    for k, r in results.items():
        if r.get("anchor_failures"):
            reasons.append(f"{k}: the floor's anchors invert or go missing on seeds "
                           f"{[f['seed'] for f in r['anchor_failures']]}")
        seeds = [row["seed"] for row in r["rows"]]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        if r.get("task") != TASK_OF.get(k):
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("truth") != k:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}")
        absent = sorted({a for row in r["rows"] for a in arms_of(k) if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if _off(r.get("fraction"), FRACTION, FRACTION_TOL):
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        spec = arm_spec(k)
        fr, ps = r.get("fractions") or {}, r.get("pseudocounts") or {}
        bad_fr = [a for a, v in spec.items() if _off(fr.get(a), v["fraction"], FRACTION_TOL)]
        bad_ps = [a for a, v in spec.items() if _off(ps.get(a), v["pseudocount"], EPS_TOL)]
        if bad_fr:
            reasons.append(f"{k}: wrong fraction on {bad_fr}")
        if bad_ps:
            reasons.append(f"{k}: wrong pseudocount on {bad_ps}")
        dirs = r.get("directions") or {}
        want_ps = {"entry": ENTRY_EPS, **EPS}
        if set(dirs) != set(want_ps) or any(_off((dirs[t] or {}).get("pseudocount"), e, EPS_TOL)
                                            for t, e in want_ps.items()):
            reasons.append(f"{k}: the direction record does not hold the entry's and each new pseudocount's set")
        else:
            entry = dirs["entry"].get("sha256") or {}
            for tag in EPS:
                new = dirs[tag].get("sha256") or {}
                same = sorted(g for g in entry if new.get(g) == entry[g])
                if not entry or set(new) != set(entry) or same:
                    reasons.append(f"{k}: the {tag} directions do not differ from the entry's in every group "
                                   f"(groups {sorted(new)} against {sorted(entry)}; equal in {same})")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the truths ran at {len(commits)} commits")
    repro = {}
    for k, r in results.items():
        repro[k] = reproduction(k, r["rows"], committed.get(k) or [])
        if not repro[k]["holds"]:
            src = "T44's whole-cut file" if k == WHOLE else "T46's file"
            reasons.append(f"{k}: {', '.join(REPRO_WHOLE if k == WHOLE else REPRO_CHAMBER)} do not reproduce their "
                           f"committed scores ({src}) within {REPRO_TOL:g} per seed")
    for k, r in results.items():
        arms = arms_of(k)
        if not r["rows"] or not all(a in row["arms"] for row in r["rows"] for a in arms):
            continue
        same = [a for a in NEW_ARMS if all(row["arms"][a]["score"] == row["arms"][BASE]["score"] for row in r["rows"])]
        if same:
            reasons.append(f"{k}: new arms score exactly as {BASE} on every seed (directions unchanged?): {same}")
        if any(not isinstance((row.get("sha256") or {}).get(a), str) for row in r["rows"] for a in arms):
            reasons.append(f"{k}: sha256 not recorded for every arm")
        for p, real in PLACEBO_OF.items():
            if p in arms and all((row.get("sha256") or {}).get(p) is not None
                                 and row["sha256"][p] == (row.get("sha256") or {}).get(real) for row in r["rows"]):
                reasons.append(f"{k}: the placebo {p}'s predictions have {real}'s sha256 on every seed (permutation "
                               "idle?)")
        unrecorded = sorted({a for row in r["rows"] for a in arms
                             if not {"zero_share", "detected_per_cell"} <= set((row.get("sparsity") or {}).get(a, {}))})
        lf_bad = any(not isinstance(row.get("scored_log_factors"), Mapping) or not row["scored_log_factors"]
                     or any(not isinstance(v, Mapping) or not set(arms) - {"b0"} <= set(v)
                            for v in row["scored_log_factors"].values()) for row in r["rows"])
        if unrecorded or lf_bad:
            reasons.append(f"{k}: the sparsity or scored log-factor records are missing (arms {unrecorded})")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def one_truth(r: Mapping, truth: str, arm: str, base: str) -> dict:
    """lib_step_read.truth_comparison (try T46's single-task test and its reported numbers) with the label ahead,
    behind or level."""
    c = lsr.truth_comparison(r, truth, arm, base)
    return {**c, "label": LEVEL[c["test"]]}


def two_cuts(results: Mapping[str, Mapping], arm: str, base: str) -> dict:
    """The adoption rule's test on the two cuts (pooled and whole): passes; with the sign reversed, loses; else
    neither; and each cut's single-truth numbers."""
    pair = {c: results[c] for c in CUTS}
    fwd, rev = sr.compare(pair, base, arm), sr.compare(pair, arm, base)
    label = PASSES if fwd["passes"] else LOSES if rev["passes"] else NEITHER
    return {"base": base, "change": arm, "label": label, "rule": fwd["rule"]["reason"],
            "reversed_rule": rev["rule"]["reason"], "cuts": {c: one_truth(results[c], c, arm, base) for c in CUTS}}


def did(r: Mapping, truth: str, real: str, base: str, placebo: str, placebo_base: str) -> dict:
    """Per seed (real - base) - (placebo - placebo_base) in board points, with the single-task test and the reversed
    test, the mean, seeds ahead, per-seed range and sd, and the parts as the difference of the two pairs' points."""
    def step(row, a, b):
        return row["arms"][a]["score"] - row["arms"][b]["score"]

    paired = {row["seed"]: (step(row, placebo, placebo_base), step(row, real, base)) for row in r["rows"]}
    fwd = keep_change({truth: paired})["tasks"][truth]
    rev = keep_change({truth: {s: (b, a) for s, (a, b) in paired.items()}})["tasks"][truth]
    d = np.array([b - a for a, b in paired.values()])
    sd = float(d.std(ddof=1)) if len(d) > 1 else float("nan")
    ahead, behind = bool(fwd.get("passes")), bool(rev.get("passes"))
    pr, pp = lsr.points_by_part(r["rows"], base, real), lsr.points_by_part(r["rows"], placebo_base, placebo)
    return {"quantity": f"({real} - {base}) - ({placebo} - {placebo_base})", "mean": float(d.mean()),
            "seeds_ahead": int((d > 0).sum()), "of": len(d), "range": [float(d.min()), float(d.max())], "sd": sd,
            "ahead": ahead, "behind": behind, "test": "ahead" if ahead else "behind" if behind else "neither",
            "points": {"de": pr["de"] - pp["de"], "mmd_u_variogram": pr["mmd_u_variogram"] - pp["mmd_u_variogram"]}}


def gene_identity(r: Mapping, truth: str, tag: str) -> dict:
    """Statement 3 on one chamber truth for one pseudocount."""
    arm = new_arm(tag, PLACEBO_STEP)
    gate = one_truth(r, truth, arm, BASE)
    d = did(r, truth, arm, BASE, new_placebo(tag), PERM_ENTRY)
    label = NO_GAIN if not gate["ahead"] else NEEDS_GENES if d["ahead"] else NOT_SEPARATED
    return {"label": label, "gate": gate, "did": d}


def matched_statement(comps: Sequence[Mapping]) -> str:
    """Statement 2 from its two comparisons on the two cuts."""
    labels = [c["label"] for c in comps]
    if all(x == PASSES for x in labels):
        return SMALLER_BETTER
    if all(x == LOSES for x in labels):
        return UNIFORM_BETTER
    return NO_DIFFERENCE


def nearest_uniform(arm: str, truth: str) -> str:
    """The uniform arm built on this truth whose quoted delivered size is nearest the arm's (the smaller on a tie)."""
    cands = list(uniform_of(truth))
    return min(cands, key=lambda u: (abs(DELIVERED[u] - DELIVERED[arm]), DELIVERED[u]))


def candidates(truth: str) -> tuple[str, ...]:
    """The arms the highest-mean report ranks: m1_r_gmom, the uniform arms and the new arms."""
    return (BASE, *uniform_of(truth), *NEW_ARMS)


def t52_reading(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    out["statement_1"] = {a: {"two_cuts": two_cuts(results, a, BASE),
                              "libraries": {k: one_truth(results[k], k, a, BASE) for k in LIBRARIES}}
                          for a in NEW_ARMS}
    s2 = [{"two_cuts": two_cuts(results, a, u), "libraries": {k: one_truth(results[k], k, a, u) for k in LIBRARIES}}
          for a, u in MATCHED]
    out["statement_2"] = {"comparisons": s2, "statement": matched_statement([c["two_cuts"] for c in s2])}
    out["outcome"] = out["statement_2"]["statement"]
    brackets = {}
    for a, us in BRACKETS.items():
        brackets[a] = {}
        for u in us:
            on_cuts = all(u in uniform_of(c) for c in CUTS)
            brackets[a][u] = {"two_cuts": two_cuts(results, a, u) if on_cuts else None,
                              "libraries": {k: one_truth(results[k], k, a, u) for k in LIBRARIES}}
    out["brackets_0.02"] = brackets
    out["statement_3"] = {k: {tag: gene_identity(results[k], k, tag) for tag in EPS} for k in CHAMBER_TRUTHS}
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in results[k]["rows"]])) for a in arms_of(k)}
             for k in TRUTHS}
    out["mean_scores"] = means
    out["highest_mean"] = {k: max(candidates(k), key=lambda a, k=k: means[k][a]) for k in TRUTHS}
    out["new_arm_differences"] = {
        k: {a: {"against_m1_r_gmom": one_truth(results[k], k, a, BASE),
                "nearest_uniform": nearest_uniform(a, k),
                "against_nearest_uniform": one_truth(results[k], k, a, nearest_uniform(a, k))} for a in NEW_ARMS}
        for k in TRUTHS}
    out["sparsity"] = {k: {a: {q: float(np.mean([row["sparsity"][a][q] for row in results[k]["rows"]]))
                               for q in ("zero_share", "detected_per_cell")} for a in arms_of(k)} for k in TRUTHS}
    out["scored_log_factors"] = {
        k: {g: {a: float(np.mean([row["scored_log_factors"][g][a] for row in results[k]["rows"]]))
                for a in results[k]["rows"][0]["scored_log_factors"][g]}
            for g in results[k]["rows"][0]["scored_log_factors"]} for k in TRUTHS}
    out["direction_stats"] = {k: results[k].get("direction_stats") for k in TRUTHS}
    out["same_positions_as_m1_r_gmom_perm"] = {
        k: {p: all(row.get("placebo", {}).get("same_positions_as_m1_r_gmom_perm", {}).get(p) is True
                   for row in results[k]["rows"]) for p in NEW_PLACEBOS} for k in CHAMBER_TRUTHS}
    return out


def _sd(x) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:.2f}"


def _line(c: Mapping) -> str:
    sk = ", ".join(f"{m} {x:+.4f}" for m, x in c["skills"].items())
    return (f"{c['change']} against {c['base']}: {c['mean']:+.2f} points, per-seed sd {_sd(c['sd'])}, ahead "
            f"{c['seeds_ahead']} of {c['of']}, per seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}, "
            f"{c['label']}; de terms {c['points']['de']:+.2f}, mmd_u and variogram "
            f"{c['points']['mmd_u_variogram']:+.2f}; skills: {sk}")


def _two_cut_lines(tc: Mapping, indent: str = "    ") -> list[str]:
    lines = [f"{indent}{tc['change']} against {tc['base']} on the two cuts: {tc['label']} ({tc['rule']})"]
    lines += [f"{indent}    on {c}: " + _line(x) for c, x in tc["cuts"].items()]
    return lines


def t52_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"Statement 2, decompression against uniform inflation at about the same delivered size: "
             f"{r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, g in v["reproduction"].items():
        lines.append(f"{k}: largest per-seed difference from the committed scores: "
                     + ", ".join(f"{a} {x:.3g}" for a, x in g["largest_difference"].items())
                     + f" (gate {REPRO_TOL:g})")
    if r["valid"]:
        for i, c in enumerate(r["statement_2"]["comparisons"]):
            lines.append(f"Statement 2, comparison {i + 1}:")
            lines += _two_cut_lines(c["two_cuts"])
            lines += [f"    on {k}: " + _line(x) for k, x in c["libraries"].items()]
        lines.append("Statement 1, each new arm against m1_r_gmom:")
        for s in r["statement_1"].values():
            lines += _two_cut_lines(s["two_cuts"])
            lines += [f"        on {k}: " + _line(x) for k, x in s["libraries"].items()]
        lines.append("The 0.02 arms against the two uniform arms that bracket each (no statement):")
        for a, us in r["brackets_0.02"].items():
            for u, x in us.items():
                if x["two_cuts"] is None:
                    lines.append(f"    {a} against {u} on the two cuts: not made ({u} is not built on {WHOLE})")
                else:
                    lines += _two_cut_lines(x["two_cuts"])
                lines += [f"        on {k}: " + _line(y) for k, y in x["libraries"].items()]
        lines.append("Statement 3, gene identity of the new directions (chamber truths):")
        for k, by in r["statement_3"].items():
            for tag, s in by.items():
                d = s["did"]
                lines.append(f"    on {k}, {tag}: {s['label']}; gate: " + _line(s["gate"]))
                lines.append(f"        on {k}, {tag}: {d['quantity']}: {d['mean']:+.2f} points, per-seed sd "
                             f"{_sd(d['sd'])}, ahead {d['seeds_ahead']} of {d['of']}, per seed {d['range'][0]:+.2f} "
                             f"to {d['range'][1]:+.2f}, test: {d['test']}; de terms {d['points']['de']:+.2f}, mmd_u "
                             f"and variogram {d['points']['mmd_u_variogram']:+.2f}")
        for k in TRUTHS:
            lines.append(f"truth {k}: highest mean {r['highest_mean'][k]}; mean scores "
                         + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
            for x in r["new_arm_differences"][k].values():
                lines.append(f"    {k}: " + _line(x["against_m1_r_gmom"]))
                lines.append(f"    {k}: nearest delivered size {x['nearest_uniform']}: "
                             + _line(x["against_nearest_uniform"]))
            sp = r["sparsity"][k]
            lines.append("    scored genes, mean over seeds: " + "; ".join(
                f"{a} zero share {x['zero_share']:.4f}, detected per cell {x['detected_per_cell']:.2f}"
                for a, x in sp.items()))
            lines.append("    mean |log factor| on the scored genes by group: " + "; ".join(
                f"{g}: " + ", ".join(f"{a} {x:.4f}" for a, x in fs.items())
                for g, fs in r["scored_log_factors"][k].items()))
            ds = r["direction_stats"][k] or {}
            lines.append("    directions (squared size over the entry's; share on genes with a Feng E9.5 mean under "
                         f"{LOW_CP10K:g} CP10k): " + "; ".join(
                             f"{tag} {g} " + ", ".join("n/a" if x is None else f"{x:.3f}" for x in rec.values())
                             for tag, by in ds.items() for g, rec in by.items()))
            if k in r["same_positions_as_m1_r_gmom_perm"]:
                lines.append("    permutation of positions the same as m1_r_gmom_perm's on every seed: " + ", ".join(
                    f"{p} {x}" for p, x in r["same_positions_as_m1_r_gmom_perm"][k].items()))
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t46_tasks: Path = T46_TASKS, t44_whole: Path = T44_WHOLE) -> dict:
    """Read the truths present in out/tasks against try T46's committed per-truth files and try T44's whole-cut file;
    write t52.{json,txt}."""
    gv.refuse_existing(out, ("t52.json", "t52.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / TASK_FILE.format(truth=k) for k in TRUTHS
               if (out / "tasks" / TASK_FILE.format(truth=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    files = {k: Path(t46_tasks) / TASK_FILE.format(truth=k) for k in CHAMBER_TRUTHS}
    files[WHOLE] = Path(t44_whole)
    files = {k: p for k, p in files.items() if p.exists()}
    committed = {k: json.loads(p.read_text())["rows"] for k, p in files.items()}
    res = {"label": LABEL_T52, "reading": t52_reading(results, committed),
           "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *files.values()])
    gv.write_outputs(out, "t52", res, t52_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T52's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one truth of try T52")
    p.add_argument("--truth", required=True, choices=list(TRUTHS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T52's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
