#!/usr/bin/env python
"""Try T58 (experiments/harness/TRIES.md): M1's step beyond twice its size, with the controls that show whether a gain
needs Feng's own genes, at 3, 4 and 5 times the practice step. A labelled in-sample read on the routes of tries T45
(experiments/t1_m1/perm_read.py), T46 (experiments/t1_m1/lib_step_read.py), T50 (experiments/t1_m1/over_read.py) and T54
(experiments/t1_m1/deep_read.py). Built to brief AM.

Tries T44 to T54 found that a larger M1 step scores better on the released day. The controls that show a gain to need
Feng's own genes (try T45's placebo direction and its there-and-back arms) exist at one and two times the practice step
only. This read runs them at 3, 4 and 5 times. Its use is fixed in docs/PLAN.md decision log 2026-11-18 (c): it is
quoted in the captain's question; it sets no step and is no evidence on the hidden truth's depth, which stays closed.

  run    One world, 10 seeds. Worlds (six):
           pooled   the chamber cut T1/E9.5->E10.5:chamber with the board's own truth (harness.filter_target, as try
                    T45's chamber run and try T46's pooled truth load it: over_read.truth_and_contamination);
           whole    the whole cut T1/E9.5->E10.5 with its truth as try T45's (and T44's) whole run loads it;
           lib3     the chamber cut with the truth library 3 as try T46 makes it (over_read.truth_and_contamination);
           lib2     the same with library 2;
           ref67    try T54's worlds of those names (deep_read.world_reference: the reference thinned by try T54's own
           ref50    code, p and stream, before anything else is made), the truth library 3.
         Arms per seed, all from one cell draw:
           b0; m1_r_gmom (fraction 0.2); m1_r_gmom_x2, _x3, _x4, _x5 (try T50's names and builds:
             step_read.step_arms on feng_dirs.feng_directions at over_read.MULTIPLIERS' values);
           m1_r_gmom_p2 (try T45's name and build), m1_r_gmom_p3, _p4, _p5: try T45's placebo build
             (gmom_combo.m1_r_arms on perm_read.placebo_directions, try T15's permutation, one per group for every
             fraction and seed) at fractions 0.4, 0.6, 0.8 and 1.0;
           m1_r_gmom_x3b2, _x4b2, _x5b2: try T45's there-and-back code (perm_read.there_and_back) from twice the step,
             fractions [0.6, -0.2], [0.8, -0.4] and [1.0, -0.6]: its first edit is the x3, x4 or x5 arm's (the same
             float, bit for bit), its second edit returns to twice the step on that code's own second stream [seed,
             perm_read.SECOND_STREAM]; x_kb2 returns to x_2's expected counts approximately (exactly only for a depth
             constant taken from expectations) and keeps x_k's zeroing (vec.core.apply_factors never refills a zero).
         Each placebo arm's record is the sha256 of the directions its own build received, per seed (as try T45's
         perm_read.placebo_record); a record that differs between seeds says so and voids.
         Not run here: try T46's placebo (a random permutation on stream 46) and try T45's second permutation (q1,
         q2); the placebo is try T45's first within-bin permutation only.
         Refused on an existing task file, unless the committed N22 passed, and unless experiments/harness/TRIES.md
         holds row T58.
  read   T58's reading of the six worlds' files (below, fixed here and repeated in the ledger row before the run).
         Writes t58.json and t58.txt and refuses to overwrite them.

The reading. Void unless: the six worlds hold the ten seeds with every arm and its fraction(s) at one commit on a clean
tree; no file records an anchor failure; every arm for which a committed score exists reproduces it seed for seed to
1e-9, namely (each try's arms that this read also builds, from the try's own ARMS):
  pooled  try T45's chamber file (experiments/t1_m1/results/t45_m1_perm/tasks/T1_E9.5_E10.5_chamber.json): b0,
          m1_r_gmom, m1_r_gmom_x2, m1_r_gmom_p2; try T46's pooled.json: b0, m1_r_gmom, _x2, _x3, _x4; try T50's
          pooled.json: b0, m1_r_gmom, _x2, _x3, _x4, _x5;
  whole   try T45's whole file (T1_E9.5_E10.5.json): b0, m1_r_gmom, m1_r_gmom_x2, m1_r_gmom_p2; try T44's whole file
          (experiments/t1_m1/results/t44_m1_step/tasks/T1_E9.5_E10.5.json): b0, m1_r_gmom, m1_r_gmom_x2;
  lib3    try T46's lib3.json: b0, m1_r_gmom, _x2, _x3, _x4; try T50's lib3.json: b0, m1_r_gmom, _x2 to _x5;
  lib2    the same in try T46's and try T50's lib2.json;
  ref67   try T54's ref67.json: b0, m1_r_gmom, _x2 to _x5;
  ref50   try T54's ref50.json: the same;
the thinned worlds pass try T54's checks (their p and stream, the depth ratio within 0.01 of p, no cell emptied, the
rebuild check); each placebo arm records the directions its build received on every seed as one set (no "differs
between seeds"), the four arms the same set, with Feng's groups and differing from Feng's in every group; that set is
the same in the six files and equals try T45's committed record (m1_r_gmom_p1 and m1_r_gmom_p2 in both of try T45's
files); no placebo or there-and-back arm scores exactly as its twin (the step arm of its first fraction: p_k's x_k,
x_kb2's x_k) on every seed of a world; no placebo arm p_k (k = 3, 4, 5) scores exactly as p_2, and no step arm x_k as
x_2, on every seed of a world (an idle fraction; in whole x_3 to x_5 have no committed score); every file's stamp
records a commit. The prediction sha256 is not compared across pooled, lib3 and lib2: harness.task_board sizes the
board from the truth, so the drawn cells differ between those worlds.
The test is try T46's single-task test inside one world ("ahead": a mean paired gain above 0.3 with at least 8 of 10
seeds ahead; "behind" the same reversed; otherwise "level"), each quantity paired by seed.
Statement 1, per world and multiple k of 3, 4 and 5, each with its mean, seeds ahead, per-seed range and label: the
step's increment R_k = x_k - x_2; the placebo's increment P_k = p_k - p_2; the specific increment S_k = R_k - P_k.
Statement 2, per world and k, the same way: V_k = x_kb2 - x_2, x_k's extra zeroing and redraw at about twice the
step's expected counts (a reading of the closed depth lever, never a reason to add redraw or sparsity), and N_k = x_k -
x_kb2, the step's expected change at about x_k's zero pattern.
Statement 3, the summary the use reads, on the six readings (k = 3 and 4 in pooled, whole and lib3), gated as try
T45's statements are gated on its (h):
  "the step's gain beyond twice its size is not in all six readings" if any of the six readings of R_k is level or
    behind, naming each such reading with R_k's and P_k's means, and where S_k reads ahead there "the placebo loses
    more" (the readings not named may still gain clearly);
  otherwise "the gain beyond twice the step needs Feng's own genes" if S_3 and S_4 read ahead in all six and P_3 and
    P_4 are ahead in none, quoting each of the six readings with R_k's and P_k's means ("not ahead" includes a placebo
    that gains up to 0.3 points on the mean, or more on fewer than 8 seeds);
  "a placebo gains as well; Feng's direction more" if S_3 and S_4 read ahead in all six and P_k is ahead in any,
    naming each such reading with R_k's and P_k's means;
  otherwise "not shown", naming each of the six where S_k is level or behind.
The three that read S_k carry "one placebo permutation, try T45's first" in the outcome itself. This is stricter than
brief AM's statement 3 (S_3 and S_4 ahead alone): S_k ahead is the expected result when a random direction of this
size loses at 3 to 5 times, whether or not the step gains, so S_k alone cannot show a gain that needs Feng's genes.
Decision log 2026-11-18 (c) fixes the use on S_3 and S_4 alone and names try T46's placebo; this read runs try T45's
first within-bin permutation in its place. The log's condition (S_3 and S_4 ahead in all six) is returned beside the
outcome as s_all_six_ahead with the six labels and decides nothing. Also beside the outcome, deciding nothing: each of
the six readings' R_k, P_k and S_k means and seeds ahead (an increment near 0.3 points gives an outcome that noise
decides), and each placebo arm's groups whose applied size is more than 10% from Feng's (perm_read.size_flags).
Reported and read by nothing: every arm's mean per world; S_5; each quantity's parts on the two de terms and on mmd_u
with the variogram; the placebo's applied size on the non-held genes beside Feng's, per fraction (try T45's record).

NOT an adoption result: a labelled in-sample diagnostic, as tries T44 to T57 are; the result can only go to the team
captain, and no outcome is evidence on Task 1's read depth, which stays closed. The limits are LIMITS below, repeated
in the text.

Usage: on the DGX campaign checkout, repo root, all six worlds at one commit, committing nothing until the six task
files exist, and on the machine tries T44, T45, T46, T50 and T54 ran on, so that the rescored arms reproduce. The
environment variables are docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/ctl_range_read.py run --world pooled --out t58_ctl_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/ctl_range_read.py run --world whole --out t58_ctl_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/ctl_range_read.py run --world lib3 --out t58_ctl_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/ctl_range_read.py run --world lib2 --out t58_ctl_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/ctl_range_read.py run --world ref67 --out t58_ctl_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/ctl_range_read.py run --world ref50 --out t58_ctl_range
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/ctl_range_read.py read --out t58_ctl_range
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import anndata as ad
import numpy as np

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import deep_read as dr  # noqa: E402
import over_read as ob  # noqa: E402
import part_depth_read as pdr  # noqa: E402
import part_range_read as prr  # noqa: E402
import perm_read as pm  # noqa: E402

ls, sr, gc, er = ob.ls, ob.sr, ob.gc, ob.er
gv, hygiene, n12, feng_dirs, zr = ob.gv, ob.hygiene, ob.n12, ob.feng_dirs, prr.zr
CHAMBER, WHOLE_CUT = gc.TASK_CHAMBER, gc.TASK_WHOLE
ROW = "T58"
STEM = "t58"
POOLED, WHOLE, LIB3, LIB2 = ls.POOLED, "whole", ob.HEART, "lib2"
WORLDS = (POOLED, WHOLE, LIB3, LIB2, "ref67", "ref50")
TASK_OF = {k: (WHOLE_CUT if k == WHOLE else CHAMBER) for k in WORLDS}
REFERENCE_OF = {"ref67": "ref67", "ref50": "ref50"}  # try T54's worlds, under their own names
TRUTH_OF = {POOLED: POOLED, WHOLE: WHOLE, LIB3: LIB3, LIB2: LIB2, "ref67": LIB3, "ref50": LIB3}
STATEMENT_WORLDS = (POOLED, WHOLE, LIB3)  # statement 3's worlds
BASE = sr.BASE  # m1_r_gmom, fraction 0.2
X2, X3, X4, X5 = ob.X2, ob.X3, ob.X4, ob.X5
STEP_AT = {2: X2, 3: X3, 4: X4, 5: X5}  # multiple -> try T50's arm
STEP_MULTIPLIERS = {a: float(ob.MULTIPLIERS[a]) for a in STEP_AT.values()}  # try T50's values
MULTIPLES = (3, 4, 5)  # the multiples the statements read
STATEMENT_MULTIPLES = (3, 4)  # statement 3's
PLACEBO_AT = {2: pm.P2, **{k: f"{BASE}_p{k}" for k in MULTIPLES}}  # multiple -> placebo arm (p2 is try T45's)
BACK_AT = {k: f"{BASE}_x{k}b2" for k in MULTIPLES}  # multiple of the first edit -> there-and-back arm
TWIN = {**{PLACEBO_AT[k]: STEP_AT[k] for k in PLACEBO_AT}, **{BACK_AT[k]: STEP_AT[k] for k in BACK_AT}}
IDLE_TWIN = {PLACEBO_AT[k]: PLACEBO_AT[2] for k in MULTIPLES}  # a placebo arm left at p2's fraction scores as p2
ARMS = ("b0", BASE, *STEP_AT.values(), *PLACEBO_AT.values(), *BACK_AT.values())
IDLE_STEP = {STEP_AT[k]: X2 for k in MULTIPLES}  # a step arm left at x2's fraction scores as x2
FRACTION, FRACTION_TOL, REPRO_TOL = sr.FRACTION, sr.FRACTION_TOL, sr.REPRO_TOL  # 0.2, 1e-12, 1e-9
T44_TASKS = RESULTS / "t44_m1_step" / "tasks"
T45_TASKS = RESULTS / "t45_m1_perm" / "tasks"
T46_TASKS = RESULTS / "t46_m1_lib_step" / "tasks"
T50_TASKS = RESULTS / "t50_m1_over_step" / "tasks"
T54_TASKS = RESULTS / "t54_m1_deep_step" / "tasks"
# The committed files per world (try -> how the try names that world's file), and each try's arms (its own ARMS)
SOURCES = {POOLED: ("t45", "t46", "t50"), WHOLE: ("t45", "t44"), LIB3: ("t46", "t50"), LIB2: ("t46", "t50"),
           "ref67": ("t54",), "ref50": ("t54",)}
SOURCE_ARMS = {"t44": sr.ARMS, "t45": pm.ARMS, "t46": ls.ARMS, "t50": ob.ARMS, "t54": dr.ARMS}
TASK_FILE = "{world}.json"


def repro_arms(src: str) -> tuple[str, ...]:
    """The arms of this read that try `src` committed scores for (from the try's own ARMS), in this read's order."""
    return tuple(a for a in ARMS if a in SOURCE_ARMS[src])


def unheld_arms() -> dict[str, tuple[str, ...]]:
    """Per world, the arms that no committed file holds (from the gates the read applies)."""
    return {k: tuple(a for a in ARMS if not any(a in repro_arms(s) for s in SOURCES[k])) for k in WORLDS}


def _increment(hi: str, lo: str) -> str:
    return f"{hi} - {lo}"


R_AT = {k: _increment(STEP_AT[k], X2) for k in MULTIPLES}  # R_k = x_k - x_2
P_AT = {k: _increment(PLACEBO_AT[k], PLACEBO_AT[2]) for k in MULTIPLES}  # P_k = p_k - p_2

NEEDS_GENES = "the gain beyond twice the step needs Feng's own genes"
PLACEBO_TOO = "a placebo gains as well; Feng's direction more"
NO_GAIN = "the step's gain beyond twice its size is not in all six readings"
NOT_SHOWN = "not shown"
ONE_PERMUTATION = "one placebo permutation, try T45's first"
FALLS_FASTER = "the placebo loses more"
LABEL_T58 = ("try T58, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng step at 1, 2, 3, 4 and 5 times, try T45's placebo direction at 2, 3, 4 "
             "and 5 times and try T45's there-and-back arms from 3, 4 and 5 times back to twice the step, on the "
             "chamber cut's pooled truth, library 3's and library 2's, the reference thinned to two thirds and a "
             "half (try T54's draws) and the whole cut; the result can only go to the team captain, it sets no step, "
             "and no outcome is evidence on Task 1's read depth, which stays closed")
LIMITS = (
    "One released day, in-sample.",
    "One placebo permutation (try T45's first).",
    "A placebo that gains as much as the step says the gain goes with factors of this size profile or with redrawn "
    "counts, not that Feng's direction is wrong.",
    "Where R_k is negative (the step has turned) a positive S_k only says the placebo falls faster.",
    "A there-and-back arm redraws more than either end.",
    "Each difference is taken inside one world.",
    "A thinned reference is one draw and is not a deeper truth.",
    "Nothing here says how deep any hidden truth is, and no statement may be used with a figure for it.",
    "x_kb2 returns to x_2's expected counts approximately (exactly only for a depth constant taken from expectations) "
    "and keeps x_k's zeroing.",
    "V_k reads x_k's extra zeroing and redraw at twice the step's expected counts; it is a reading of the closed depth "
    "lever and never a reason to add redraw or sparsity; N_k compares x_k with an arm that shares its zeros, not with "
    "x_2.",
    "The six readings of statement 3 are two cuts (chamber, whole); pooled and lib3 share the reference, the builds "
    "and part of the truth, and S_3 and S_4 share x_2 and p_2.",
    "Try T46's placebo (random permutation on stream 46) is not run here; only try T45's first within-bin "
    "permutation, and not its second (q1, q2).",
    "Decision log 2026-11-18 (c) names try T46's placebo and reads S_3 and S_4 alone; this read runs try T45's first "
    "within-bin permutation in its place, gates on R_k and splits on P_k. The log's own condition (S_3 and S_4 ahead "
    "in all six readings) is reported beside the outcome and decides nothing: with R_k level and a placebo that loses "
    "as k grows it holds whether or not the step gains.",
    "P_k 'not ahead' by try T46's test includes a placebo that gains up to 0.3 points on the mean, or more on fewer "
    "than 8 of 10 seeds; the outcome that the gain needs Feng's own genes quotes each of the six readings' R_k and P_k "
    "means.",
    "An increment near 0.3 points gives an outcome that noise decides between the four; each of the six readings' "
    "R_k, P_k and S_k means and seeds ahead are printed beside the outcome.",
    "The placebo's applied size on the non-held genes can differ from Feng's (the permutation moves values onto and "
    "off the held genes, more as k grows); a larger placebo loses more and raises S_k. Groups more than 10% from "
    "Feng's are flagged beside statement 3 (perm_read.size_flags); the flag decides nothing.",
)
STATEMENTS = (
    "Statement 1, per world and multiple k of 3, 4 and 5, each paired by seed with its mean, seeds ahead, per-seed "
    "range and its label by try T46's single-task test: the step's increment R_k = x_k - x_2; the placebo's increment "
    "P_k = p_k - p_2; the specific increment S_k = R_k - P_k.",
    "Statement 2, per world and k, the same way: V_k = x_kb2 - x_2, x_k's extra zeroing and redraw at about twice the "
    "step's expected counts (a reading of the closed depth lever, never a reason to add redraw or sparsity), and N_k = "
    "x_k - x_kb2, the step's expected change at about x_k's zero pattern.",
    f"Statement 3, on the six readings (k = 3 and 4 in pooled, whole and lib3), gated as try T45's statements are on "
    f"its (h): '{NO_GAIN}' if any of the six readings of R_k is level or behind, naming each with R_k's and P_k's "
    f"means and, where S_k reads ahead there, '{FALLS_FASTER}' (the other readings may still gain); otherwise "
    f"'{NEEDS_GENES}' if S_3 and S_4 read ahead by try T46's test in all six and P_3 and P_4 in none, quoting each of "
    f"the six with R_k's and P_k's means ('not ahead' includes a placebo gain below 0.3 points); '{PLACEBO_TOO}' if "
    "S_3 and S_4 read ahead in all six and P_k is ahead in any, naming each such reading with R_k's and P_k's means; "
    f"otherwise '{NOT_SHOWN}', naming each of the six where S_k is level or behind. The three that read S_k carry "
    f"'{ONE_PERMUTATION}' in the outcome. Beside the outcome, deciding nothing: whether S_3 and S_4 read ahead in "
    "all six (decision log 2026-11-18 (c)'s condition, with the six labels), each of the six readings' R_k, P_k and "
    "S_k means and seeds ahead, and the placebo arms' applied sizes more than 10% from Feng's.",
)
CHOICES = (
    "The there-and-back arms call perm_read.there_and_back with the fraction of twice the step (m1_r_gmom_x2's) and "
    "the multiplier k/2: its first edit is then (2 x 0.2) x (k/2), the same float as the x_k arm's 0.2 x k (both are "
    "the one rounding of the exact product), and its second edit (1 - k/2) x (2 x 0.2) the same float as (2 - k) x "
    "0.2; no copy of try T45's code.",
    "The three there-and-back arms share try T45's second stream (index perm_read.SECOND_STREAM after the build seed), "
    "as try T45's two there-and-back arms do (each is its own build and uses it once); the brief names that code's own "
    "stream, so no new stream index is drawn.",
    "A placebo or there-and-back arm's twin is the step arm of its first fraction (p_k's and x_kb2's twin is x_k): an "
    "idle permutation or an idle second edit scores exactly as x_k.",
    "The placebo record is checked against try T45's committed m1_r_gmom_p1 and m1_r_gmom_p2 records in both of try "
    "T45's files; a missing try T45 file voids.",
    "S_k is the paired per-seed difference (x_k - x_2) - (p_k - p_2), tested as try T45's (a) is (the increments as "
    "arms, perm_read._arm_difference, with try T50's comparison); R_k, P_k, V_k and N_k are try T50's comparison of "
    "the two arms.",
    "Statement 3 names each of the six readings that is not ahead, with its label ('level' or 'behind'); the gate on "
    "R_k and the split on P_k follow try T45's statements ((h) gates; (b) passing gives 'a placebo gains as well').",
    "Each placebo arm's record is the direction sha256 of the object its own build received, collected per seed and "
    "folded by perm_read.placebo_record; the read also voids a record equal to Feng's in any group (try T45's "
    "condition).",
    "A placebo arm p_k (k = 3, 4, 5) scoring exactly as p_2 on every seed of a world voids (an idle fraction); the "
    "prediction sha256 is not compared across worlds, because the board's cell count follows each world's truth.",
    "The whole cut's truth is loaded as tries T44 and T45 load it (harness.filter_target when the task has a target "
    "filter, otherwise the released target on the panel); its file records truth 'whole'.",
    "The thinned worlds' checks are try T57's statement of try T54's (part_range_read._thinning_reasons).",
    "The run's --out folder for the command lines is t58_ctl_range.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def arm_fractions(fraction: float) -> dict[str, float | list[float]]:
    """Every arm's fraction; the there-and-back arms' two, in the order applied, as perm_read.there_and_back computes
    them from twice the step and the multiplier k/2."""
    fr = sr.arm_fractions(fraction, STEP_MULTIPLIERS)
    twice = fr[X2]
    out = {**fr, **{PLACEBO_AT[k]: fr[STEP_AT[k]] for k in PLACEBO_AT}}
    for k, a in BACK_AT.items():
        m = float(k) / 2.0
        out[a] = [twice * m, (1.0 - m) * twice]
    return out


def build_seed(ctx: harness.BuildContext, b0, directions: Mapping, placebo: Mapping, groups: Mapping, fraction: float,
               held_m1: Sequence[str], held_factors: Sequence[str], used: dict | None = None) -> dict[str, Stage]:
    """One seed's arms: b0, m1_r_gmom and x2 to x5 as try T50 builds them (step_read.step_arms); the placebo arms as
    try T45 builds p2 (gmom_combo.m1_r_arms on the placebo directions) at each fraction; the there-and-back arms by try
    T45's perm_read.there_and_back from twice the step with the multiplier k/2. `used`, when given, receives each
    placebo arm's directions object as its build received it (as perm_read.perm_arms' "directions")."""
    preds = {"b0": b0(ctx), **sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors,
                                           multipliers=STEP_MULTIPLIERS)}
    fr = arm_fractions(fraction)
    for a in PLACEBO_AT.values():
        preds[a] = gc.m1_r_arms(ctx, placebo, groups, fr[a], held_m1, held_factors)[BASE]
        if used is not None:
            used[a] = placebo
    twice = fr[X2]
    for k, a in BACK_AT.items():
        preds[a] = pm.there_and_back(ctx, directions, groups, twice, held_m1, held_factors, multiplier=float(k) / 2.0)
    return preds


# ---------------------------------------------------------------------------------------------------------------------
# The worlds and the run


def world_route(world: str, t, data: Path, genes: Sequence[str], inputs: Mapping[str, ad.AnnData], role: str,
                ref_stage: Stage, read_full: Callable[[], ad.AnnData]
                ) -> tuple[Mapping[str, ad.AnnData], Stage, Stage, int | None, dict | None, str]:
    """The world's inputs, scorer's reference stage, truth, contamination count, thinning record and truth name. A
    thinned world: the reference thinned first by deep_read.world_reference under its own name; then the truth: the
    whole cut's as tries T44 and T45 load it, otherwise over_read.truth_and_contamination (try T46's for these
    names)."""
    if world not in WORLDS:
        raise ValueError(f"unknown world {world!r}; one of {list(WORLDS)}")
    thinning = None
    if world in REFERENCE_OF:
        inputs, ref_stage, thinning = dr.world_reference(REFERENCE_OF[world], inputs, role, ref_stage, genes,
                                                         read_full)
    if world == WHOLE:
        truth, rec = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                      (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
        contamination = None if rec is None else int(rec["contamination"])
    else:
        truth, contamination, _ = ob.truth_and_contamination(t, data, genes, TRUTH_OF[world])
    return inputs, ref_stage, truth, contamination, thinning, TRUTH_OF[world]


def world_placebo(pool_genes: Sequence[str]) -> tuple[dict, dict, dict, str]:
    """Feng's directions as tries T45 and T50 load them, try T45's placebo on them (perm_read.placebo_directions on
    eb_read.feng_halves), their external record and the statistics' sha256. Raises ValueError where the statistics the
    halves came from are not the directions'."""
    directions, ext = feng_dirs.feng_directions(list(pool_genes))
    halves, stats_sha = er.feng_halves(list(pool_genes), directions)
    if stats_sha != ext["sha256"]:
        raise ValueError("the statistics the halves came from are not the ones M1's direction came from")
    return directions, pm.placebo_directions(directions, halves), ext, stats_sha


def applied_sizes(directions: Mapping, placebo: Mapping, held: Sequence[str], genes: Sequence[str], fraction: float
                  ) -> dict:
    """Per placebo arm, try T45's applied-size record at its fraction (Feng and the placebo, non-held genes)."""
    fr = arm_fractions(fraction)
    return {a: {"fraction": fr[a], "feng": pm.applied_size(directions, held, genes, fr[a]),
                "placebo": pm.applied_size(placebo, held, genes, fr[a])} for a in PLACEBO_AT.values()}


def seed_rows(view, cboard, board, truth: Stage, ref_stage: Stage, inputs: Mapping[str, ad.AnnData], floors: Mapping,
              seeds: Sequence[int], directions: Mapping, placebo: Mapping, groups: Mapping, fraction: float,
              held_m1: Sequence[str], preregistered: Sequence[str], label: str = "",
              built_sha: dict[str, set[str]] | None = None) -> list[dict]:
    """Per seed: the arms built from the world's inputs and scored against the world's truth with the world's reference
    stage and that seed's anchors, recorded as tries T50 and T54 record them. `built_sha`, when given, collects per
    placebo arm the direction sha256 of what its build received on each seed (try T45's built_sha)."""
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        used = {}
        preds = build_seed(ctx, b0, directions, placebo, groups, fraction, held_m1, preregistered, used=used)
        if built_sha is not None:
            for a in PLACEBO_AT.values():
                built_sha.setdefault(a, set()).add(json.dumps(pm.direction_sha(used[a]), sort_keys=True))
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "anchors": ob.anchor_record(f.anchors), "arms": arms,
                     "sha256": {k: ls.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: ls.sparsity(p) for k, p in preds.items()}, "seconds": round(time.time() - t0, 1)})
        print(f"{label} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    return rows


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T58 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T58 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = TASK_FILE.format(world=a.world)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK_OF[a.world]]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    role = harness.base_role(t)
    try:  # a thinned reference comes before anything else is made (try T54's order)
        inputs, ref_stage, truth, contamination, thinning, truth_name = world_route(
            a.world, t, data, genes, inputs, role, ref_stage, lambda: read_stage(data / t.reference, None))
    except ValueError as e:
        sys.exit(f"T58 is refused on {a.world}: {e}")
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T58, "task": t.key, "world": a.world, "truth": truth_name, "truth_rows": int(truth.n),
           "truth_contamination": contamination, "seeds": list(seeds), "thinning": thinning}
    floors, failures = ls.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        zr.write_task_file(out, name, res)
        sys.exit(f"T58 is refused on {a.world}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T58 is void)")
    pool = inputs[role]
    try:
        directions, placebo, ext, stats_sha = world_placebo(list(pool.var_names))
    except ValueError as e:
        sys.exit(f"T58 is refused on {a.world}: {e}")
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    built_sha = {p: set() for p in PLACEBO_AT.values()}
    rows = seed_rows(view, cboard, board, truth, ref_stage, inputs, floors, seeds, directions, placebo, groups,
                     fraction, held_m1, preregistered, label=a.world, built_sha=built_sha)
    res.update({"fraction": fraction, "fractions": arm_fractions(fraction), "multipliers": dict(STEP_MULTIPLIERS),
                "second_stream": pm.SECOND_STREAM, "feng_sha256": pm.direction_sha(directions),
                "placebo_sha256": pm.placebo_record(built_sha),
                "applied_size": applied_sizes(directions, placebo, held_m1, list(pool.var_names), fraction),
                "feng_stats_sha256": stats_sha, "external": ext, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    zr.write_task_file(out, name, res)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _rows(r: Mapping) -> list:
    return ob._rows(r)


def committed_files(t44_tasks: Path = T44_TASKS, t45_tasks: Path = T45_TASKS, t46_tasks: Path = T46_TASKS,
                    t50_tasks: Path = T50_TASKS, t54_tasks: Path = T54_TASKS) -> dict[str, dict[str, Path]]:
    """try -> world -> the committed task file the read reproduces (SOURCES): tries T44 and T45 by cut file, tries T46
    and T50 by truth, try T54 by world."""
    folder = {"t44": Path(t44_tasks), "t45": Path(t45_tasks), "t46": Path(t46_tasks), "t50": Path(t50_tasks),
              "t54": Path(t54_tasks)}

    def file(src: str, k: str) -> Path:
        if src in ("t44", "t45"):
            return folder[src] / gc.CUT_FILES[TASK_OF[k]]
        return folder[src] / TASK_FILE.format(world=TRUTH_OF[k] if src in ("t46", "t50") else k)

    return {src: {k: file(src, k) for k in WORLDS if src in SOURCES[k]} for src in folder}


def t45_placebo_records(committed: Mapping[str, Mapping[str, Mapping]]) -> dict[str, dict | None]:
    """Try T45's committed placebo record per cut file (m1_r_gmom_p1 and m1_r_gmom_p2), None where a file is
    missing."""
    out = {}
    for k in (POOLED, WHOLE):
        r = (committed.get("t45") or {}).get(k)
        out[gc.CUT_FILES[TASK_OF[k]]] = None if r is None else (r.get("placebo_sha256") or {})
    return out


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Mapping]]) -> dict:
    """The void conditions. results: world -> task file; committed: try -> world -> its committed file (a missing file
    absent)."""
    reasons = []
    missing = [k for k in WORLDS if k not in results]
    if missing:
        reasons.append(f"missing worlds {missing}")
    want_fr = arm_fractions(FRACTION)
    repro, placebos = {}, {}
    for k, r in results.items():
        rows = _rows(r)
        if r.get("anchor_failures"):
            reasons.append(f"{k}: the floor's anchors invert or go missing on seeds "
                           f"{[f['seed'] for f in r['anchor_failures']]}")
        seeds = [row["seed"] for row in rows]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in rows for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if k in TASK_OF and r.get("task") != TASK_OF[k]:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("world") != k:
            reasons.append(f"{k}: the file holds world {r.get('world')!r}")
        if k in TRUTH_OF and r.get("truth") != TRUTH_OF[k]:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}, not {TRUTH_OF[k]}")
        if not pm._fraction_ok(r.get("fraction"), FRACTION):
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items() if not pm._fraction_ok(got_fr.get(a), f)]
        if bad:
            reasons.append(f"{k}: arm fractions {got_fr}, not {want_fr} (wrong fraction on {bad})")
        repro[k] = {}
        for src in SOURCES.get(k, ()):
            arms = repro_arms(src)
            rep = pdr.reproduction(rows, committed.get(src, {}).get(k), arms)
            repro[k][src] = rep
            if not rep["holds"]:
                reasons.append(f"{k}: {', '.join(arms)} do not reproduce try {src.upper()}'s committed scores within "
                               f"{REPRO_TOL:g} per seed")
        th = r.get("thinning")
        if k in REFERENCE_OF:
            reasons.extend(prr._thinning_reasons(k, th, results))
        elif th is not None:
            reasons.append(f"{k}: the file records a thinning")
        rec = r.get("placebo_sha256") or {}
        recs = [rec.get(a) for a in PLACEBO_AT.values()]
        feng = r.get("feng_sha256") or {}
        split = [a for a, x in zip(PLACEBO_AT.values(), recs, strict=True)
                 if isinstance(x, Mapping) and pm.DIFFERS in x]
        if split:
            reasons.append(f"{k}: placebo arms {split} were built with different directions on different seeds")
        elif not all(isinstance(x, Mapping) and x for x in recs) or any(x != recs[0] for x in recs):
            reasons.append(f"{k}: the placebo arms {list(PLACEBO_AT.values())} do not record one set of permuted "
                           "directions")
        elif not feng or set(recs[0]) != set(feng):
            reasons.append(f"{k}: the placebo's groups {sorted(recs[0])} are not Feng's {sorted(feng)}")
        else:
            same_g = sorted(g for g in feng if recs[0][g] == feng[g])
            if same_g:
                reasons.append(f"{k}: the placebo arms record Feng's direction in groups {same_g}")
            placebos[k] = json.dumps(recs[0], sort_keys=True)
        if rows and not absent:
            same = [a for a, twin in TWIN.items()
                    if all(row["arms"][a]["score"] == row["arms"][twin]["score"] for row in rows)]
            if same:
                reasons.append(f"{k}: placebo or there-and-back arms score exactly as their twin on every seed: {same}")
            idle = [a for a, twin in IDLE_TWIN.items()
                    if all(row["arms"][a]["score"] == row["arms"][twin]["score"] for row in rows)]
            if idle:
                reasons.append(f"{k}: placebo arms score exactly as {PLACEBO_AT[2]} on every seed (an idle fraction): "
                               f"{idle}")
            idle_x = [a for a, twin in IDLE_STEP.items()
                      if all(row["arms"][a]["score"] == row["arms"][twin]["score"] for row in rows)]
            if idle_x:
                reasons.append(f"{k}: step arms score exactly as {X2} on every seed (an idle fraction): {idle_x}")
    if len(set(placebos.values())) > 1:
        reasons.append(f"the placebo's recorded sha256 differ between the worlds' files {sorted(placebos)}")
    for cut, t45 in t45_placebo_records(committed).items():
        if t45 is None:
            reasons.append(f"try T45's committed file {cut} is missing (its placebo record cannot be checked)")
            continue
        want = {json.dumps(t45.get(a), sort_keys=True) for a in pm.PLACEBO_ARMS}
        off = sorted(k for k, v in placebos.items() if want != {v})
        if off:
            reasons.append(f"the placebo's recorded sha256 differ from try T45's committed record ({cut}) in {off}")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    no_commit = sorted(k for k, r in results.items() if not (r.get("stamp") or {}).get("commit"))
    if no_commit:
        reasons.append(f"no commit recorded on {no_commit}")
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the worlds ran at {len(commits)} commits")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def comparison(r: Mapping, world: str, change: str, base: str) -> dict:
    """Try T50's comparison (try T46's single-task test on one world, the reversed test, the points by part) with the
    label ahead / behind / level."""
    return ob.comparison(r, world, change, base)


def with_increments(r: Mapping) -> dict:
    """A world's file with R_k and P_k as arms of each row (perm_read._arm_difference: score, skills, groups and
    metrics), seed by seed, so that S_k is try T50's comparison of the pair (P_k, R_k)."""
    def incs(arms: Mapping) -> dict:
        return {**{R_AT[k]: pm._arm_difference(arms[STEP_AT[k]], arms[X2]) for k in MULTIPLES},
                **{P_AT[k]: pm._arm_difference(arms[PLACEBO_AT[k]], arms[PLACEBO_AT[2]]) for k in MULTIPLES}}

    return {**r, "rows": [{**row, "arms": {**row["arms"], **incs(row["arms"])}} for row in _rows(r)]}


def world_quantities(r: Mapping, world: str) -> dict:
    """Per k of 3, 4, 5: R_k, P_k, S_k (statement 1) and V_k, N_k (statement 2), each try T50's comparison."""
    inc = with_increments(r)
    out = {}
    for k in MULTIPLES:
        out[k] = {"R": comparison(r, world, STEP_AT[k], X2), "P": comparison(r, world, PLACEBO_AT[k], PLACEBO_AT[2]),
                  "S": comparison(inc, world, R_AT[k], P_AT[k]),
                  "V": comparison(r, world, BACK_AT[k], X2), "N": comparison(r, world, STEP_AT[k], BACK_AT[k])}
    return out


def _beside(k: str, m: int, q: Mapping) -> str:
    return f"{k} R_{m} {q['R']['label']} {q['R']['mean']:+.2f} beside P_{m} {q['P']['label']} {q['P']['mean']:+.2f}"


def statement_3(quantities: Mapping[str, Mapping[int, Mapping]]) -> dict:
    """Statement 3 on the six readings (k = 3 and 4 in pooled, whole and lib3), gated on R_k as try T45's statements are
    on its (h):
      NO_GAIN if any reading of R_k is not ahead, naming each with R_k's and P_k's means (and FALLS_FASTER where S_k
        reads ahead there);
      otherwise NEEDS_GENES if every S_k reads ahead and no P_k does;
      PLACEBO_TOO if every S_k reads ahead and some P_k does, naming those with R_k's and P_k's means;
      otherwise NOT_SHOWN, naming each S_k reading that is level or behind.
    The three that read S_k carry ONE_PERMUTATION in the outcome; NEEDS_GENES quotes the six readings' R_k and P_k
    means ('P not ahead' includes a placebo gain below 0.3 points). Returned beside the outcome and deciding nothing:
    s_all_six_ahead, decision log 2026-11-18 (c)'s own condition (S_3 and S_4 ahead in all six), with the six labels
    in `readings`."""
    six = [(k, m) for k in STATEMENT_WORLDS for m in STATEMENT_MULTIPLES]
    lab = {x: {q: quantities[x[0]][x[1]][q]["label"] for q in ("R", "P", "S")} for x in six}
    readings = {f"{k} S_{m}": lab[(k, m)]["S"] for k, m in six}
    r_readings = {f"{k} R_{m}": lab[(k, m)]["R"] for k, m in six}
    p_readings = {f"{k} P_{m}": lab[(k, m)]["P"] for k, m in six}
    off = {f"{k} S_{m}": lab[(k, m)]["S"] for k, m in six if lab[(k, m)]["S"] != "ahead"}
    no_gain = [(k, m) for k, m in six if lab[(k, m)]["R"] != "ahead"]
    placebo_up = [(k, m) for k, m in six if lab[(k, m)]["P"] == "ahead"]
    if no_gain:
        kind = NO_GAIN
        said = f"{NO_GAIN} (" + "; ".join(
            _beside(k, m, quantities[k][m]) + (f", S_{m} ahead: {FALLS_FASTER}" if lab[(k, m)]["S"] == "ahead" else "")
            for k, m in no_gain) + ")"
    elif not off and not placebo_up:
        kind = NEEDS_GENES
        said = f"{NEEDS_GENES} ({ONE_PERMUTATION}; " + "; ".join(_beside(k, m, quantities[k][m]) for k, m in six) + ")"
    elif not off:
        kind = PLACEBO_TOO
        said = f"{PLACEBO_TOO} ({ONE_PERMUTATION}; placebo ahead: " + "; ".join(
            _beside(k, m, quantities[k][m]) for k, m in placebo_up) + ")"
    else:
        kind = NOT_SHOWN
        said = f"{NOT_SHOWN} ({ONE_PERMUTATION}; " + ", ".join(f"{n} {v}" for n, v in off.items()) + ")"
    return {"statement": said, "kind": kind, "readings": readings, "s_all_six_ahead": not off, "r_readings": r_readings,
            "p_readings": p_readings, "not_ahead": off, "r_not_ahead": [f"{k} R_{m}" for k, m in no_gain],
            "p_ahead": [f"{k} P_{m}" for k, m in placebo_up]}


def t58_reading(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Mapping]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    q = {k: world_quantities(results[k], k) for k in WORLDS}
    s3 = statement_3(q)
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in _rows(results[k])])) for a in ARMS}
             for k in WORLDS}
    out.update({"outcome": s3["statement"], "statement_3": s3,
                "quantities": {k: {str(m): x for m, x in by.items()} for k, by in q.items()},
                "mean_scores": means,
                "applied_size": {k: results[k].get("applied_size") for k in WORLDS},
                "size_flags": {k: {a: pm.size_flags((results[k].get("applied_size") or {}).get(a))
                                   for a in PLACEBO_AT.values()} for k in WORLDS},
                "thinning": {k: {x: (results[k].get("thinning") or {}).get(x) for x in
                                 ("p", "stream", "depth_ratio", "rebuild_largest_difference", "any_cell_emptied")}
                             for k in REFERENCE_OF}})
    return out


def _line(name: str, c: Mapping) -> str:
    return (f"{name} ({c['change']} minus {c['base']}): {c['mean']:+.2f} points, ahead {c['seeds_ahead']} of "
            f"{c['of']}, per seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}, test: {c['label']}; de terms "
            f"{c['points']['de']:+.2f}, mmd_u and variogram {c['points']['mmd_u_variogram']:+.2f} (read by nothing)")


def _size_line(arm: str, sizes: Mapping | None) -> str:
    if not sizes:
        return f"{arm}: not recorded"
    ratios = pm.size_ratios(sizes)
    return (f"{arm} at fraction {sizes.get('fraction')} (mean |log f|, sum log f^2, non-held genes): " + "; ".join(
        f"{g} Feng {fs['mean_abs_logf']:.4f}, {fs['sum_sq_logf']:.3f}, placebo "
        f"{sizes['placebo'][g]['mean_abs_logf']:.4f}, {sizes['placebo'][g]['sum_sq_logf']:.3f}, ratio "
        + ", ".join("n/a" if x is None else f"{x:.3f}" for x in ratios.get(g, {}).values())
        for g, fs in sorted(sizes.get("feng", {}).items()) if g in sizes.get("placebo", {})))


def t58_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"Outcome: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, by in v["reproduction"].items():
        for src, x in by.items():
            lines.append(f"{k}: largest per-seed difference from try {src.upper()}'s committed scores: "
                         + ", ".join(f"{a} {d:.3g}" for a, d in x["largest_difference"].items())
                         + f" (gate {REPRO_TOL:g})")
    if r["valid"]:
        s3 = r["statement_3"]
        lines.append(f"Statement 3 (the summary the use reads): {s3['statement']}; the six readings: "
                     + ", ".join(f"{n} {lab}" for n, lab in s3["readings"].items())
                     + "; beside them (the gate and the split): "
                     + ", ".join(f"{n} {lab}" for n, lab in {**s3["r_readings"], **s3["p_readings"]}.items()))
        lines.append("Beside statement 3, deciding nothing: decision log 2026-11-18 (c)'s condition (S_3 and S_4 ahead "
                     f"in pooled, whole and lib3), all six ahead: {s3['s_all_six_ahead']} ("
                     + ", ".join(f"{n} {lab}" for n, lab in s3["readings"].items()) + ")")
        lines.append("Beside statement 3, the six readings' means and seeds ahead: " + "; ".join(
            f"{k} k={m}: " + ", ".join(
                f"{name}_{m} {c['mean']:+.2f} ({c['seeds_ahead']} of {c['of']} ahead)"
                for name, c in ((n, r["quantities"][k][str(m)][n]) for n in ("R", "P", "S")))
            for k in STATEMENT_WORLDS for m in STATEMENT_MULTIPLES))
        flags = {(k, a): f for k, by in r["size_flags"].items() for a, f in by.items() if f}
        lines.append("Beside statement 3, the placebo's applied size more than "
                     f"{pm.SIZE_FLAG:.0%} from Feng's (reported, decides nothing): " + (
                         "; ".join(f"{k} {a}: " + ", ".join(f) for (k, a), f in flags.items()) if flags else "none"))
        lines.append("Statement 1 (R_k = x_k - x_2, P_k = p_k - p_2, S_k = R_k - P_k; S_5 reported and read by "
                     "nothing) and statement 2 (V_k = x_kb2 - x_2, N_k = x_k - x_kb2), per world and k, by try T46's "
                     "test inside the world:")
        for k in WORLDS:
            for m in MULTIPLES:
                x = r["quantities"][k][str(m)]
                for name in ("R", "P", "S", "V", "N"):
                    lines.append(f"    {k}, {name}_{m}: " + _line(f"{name}_{m}", x[name]))
        lines.append("Reported and read by nothing:")
        for k in WORLDS:
            lines.append(f"    world {k}: mean scores "
                         + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
        for k in WORLDS:
            sizes = r["applied_size"].get(k) or {}
            for a in PLACEBO_AT.values():
                lines.append(f"    world {k}: the placebo's applied size beside Feng's, " + _size_line(a, sizes.get(a)))
        for k, th in r["thinning"].items():
            lines.append(f"    world {k}: realised depth ratio {th['depth_ratio']:.4f} (p {dr.THINNED[k]['p']:.4f}, "
                         f"stream {th['stream']}), rebuild check's largest difference "
                         f"{th['rebuild_largest_difference']:.3g}, any cell left with no count: "
                         f"{th['any_cell_emptied']}")
    lines.append("The statements as the code reads them (fixed before the run):")
    lines.extend(f"  - {x}" for x in STATEMENTS)
    lines.append("Arms no committed file holds, per world: " + "; ".join(
        f"{k} {len(a)} ({', '.join(a)})" for k, a in unheld_arms().items()))
    lines.append("Choices the brief left open:")
    lines.extend(f"  - {x}" for x in CHOICES)
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain; it sets no step; no outcome is evidence on Task 1's read depth, "
                 "which stays closed")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t44_tasks: Path = T44_TASKS, t45_tasks: Path = T45_TASKS, t46_tasks: Path = T46_TASKS,
         t50_tasks: Path = T50_TASKS, t54_tasks: Path = T54_TASKS) -> dict:
    """Read the worlds present in out/tasks against the committed files of tries T44, T45, T46, T50 and T54; write
    t58.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / TASK_FILE.format(world=k) for k in WORLDS
               if (out / "tasks" / TASK_FILE.format(world=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    files = committed_files(t44_tasks, t45_tasks, t46_tasks, t50_tasks, t54_tasks)
    found = {src: {k: p for k, p in by.items() if p.exists()} for src, by in files.items()}
    committed = {src: {k: json.loads(p.read_text()) for k, p in by.items()} for src, by in found.items()}
    res = {"label": LABEL_T58, "statements": list(STATEMENTS), "choices": list(CHOICES), "limits": list(LIMITS),
           "unheld_arms": unheld_arms(), "reading": t58_reading(results, committed),
           "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *sorted({p for by in found.values() for p in by.values()})])
    gv.write_outputs(out, STEM, res, t58_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T58's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one world of try T58")
    p.add_argument("--world", required=True, choices=list(WORLDS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T58's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
