#!/usr/bin/env python
"""Try T48 (experiments/harness/TRIES.md): M1's step set by its own calibration, on the two E10.5 cuts. A labelled read
of one rule on try T44's route (experiments/t1_m1/step_read.py).

M1 applies Feng's E9.5 -> E14.5 change per group as the count factor exp(fraction x d), d the log ratio of mean CP10k
made with the config's pseudocount. Read back with try T47's share code (experiments/t1_m1/share_read.py), M1's step
delivers only part of its nominal size in the units in which d is made; in those units even a true fraction of Feng's
change reads under its nominal size, so the rule's target is not the reading of the change M1 aims at (see the limits).
The rule, with no free number: set the nominal fraction so that M1's own step, read back on the task's reference cells,
delivers the bracket's fraction times the scale.

  calibrated_fraction(pool, directions, groups, held, fraction, target_share)
         cells: every cell of the pool whose label maps to a direction group, all of them (not the cells a build
         draws from the pool at the chamber shares, which may repeat or leave cells out). delivered(f):
         t1_scope.external_shift at f on those cells (the build's held list, stream [CAL_ROOT, CAL_STREAM], which no
         build step uses, its defaults otherwise), read before -> after with share_read's code: per group share_stats
         of dirlib.lfc_ratio of the mean CP10k (the config's pseudocount) on the genes M1 moves (a finite direction,
         not held), a group under the config's min_cells_t1 giving no value (share_read.calibration's rule), combined
         with share_read.combine weighted by the groups' cell numbers among those cells (the pool's group counts). The
         answer is the f in [fraction, 5 x fraction] with delivered(f) = target_share, found by bisection to a relative
         tolerance of 1e-3 on the share, on the same stream at every step. It refuses (ValueError) a target outside
         [delivered(fraction), delivered(5 x fraction)], a delivered share that does not rise across the bracket, and a
         bisection that does not reach the tolerance in MAX_STEPS steps; it never extrapolates. With a target within
         the tolerance of delivered(fraction) it returns the fraction itself.
  calibrate  Both answers and, reported only: the stream noise at each answer (sd of delivered over the ten further
         streams [CAL_ROOT, CAL_STREAM, k], k = 1..10, and that sd divided by the slope of delivered between the base
         fraction and the answer: the fraction's sd); and per arm's fraction (0.2, 0.3, 0.4 and both answers), on the
         rule's stream, the delivered share and the share over f in four bins of the group's mean CP10k before the step
         ([0, 0.1), [0.1, 1), [1, 10), >= 10), and the covered cells' mean detected genes and total counts per cell
         before and after the step, and, given the board's scored genes (the panel), the combined share over the
         panel genes M1 moves and the share over f.
  run    One E10.5 cut, 10 seeds, on try T44's route. Refused unless experiments/harness/TRIES.md holds row T48, unless
         the committed N22 passed, on an existing task file, and unless T47's committed reading
         (results/t47_share_read/t47.json) is there, not void (its stamp's tree clean) and its raw L finite and in (0.5,
         1.5). Before the seed loop the calibration runs once on the task's reference pool (the build's directions,
         groups and held list; the board's panel genes for the reported panel share only): T48a's fraction for
         target_share = 0.2 x L and T48b's for target_share = 0.2, written to the task file with delivered(0.2), the
         delivered shares at the answers, the steps, L and the sha256 of t47.json. Arms per seed, all from one cell
         draw: b0, m1_r_gmom (0.2), m1_r_gmom_x15 (0.3) and m1_r_gmom_x2 (0.4) as T44 builds them (step_read.step_arms);
         m1_r_gmom_calL (T48a) and m1_r_gmom_cal (T48b), gmom_combo.m1_r_arms' m1_r_gmom at the calibrated fraction with
         every other argument the base's. Per seed and arm the file also records the built cells' mean detected panel
         genes and mean summed panel CP10k per cell (reported only; a build's total counts are not in its Stage).
  read   T48's reading of both cuts' files (fixed in the row before the run):
         Void unless both cuts hold the ten seeds at one commit on a clean tree; b0 and m1_r_gmom equal, seed for seed,
         their committed scores in results/t42_m1_eb/tasks/ and m1_r_gmom_x15 and m1_r_gmom_x2 theirs in
         results/t44_m1_step/tasks/, each to 1e-9; each file's task is its cut, its base fraction 0.2 and every arm's
         fraction its own; both calibrated fractions lie in (0.2, 1.0) and the scaled one is the smaller; the delivered
         shares at the answers equal their targets (0.2 x L and 0.2) to 1e-3 relative; both cuts record the same
         t47_sha256, the same L and the same two calibrated fractions (one reference pool, one calibration); neither
         calibrated arm scores exactly as m1_r_gmom on every seed of both cuts.
         The test: the adoption rule (vec.adoption.keep_change, paired by seed, both cuts) for each calibrated arm
         against m1_r_gmom. Statement, on T48a alone: "the calibrated step passes on both cuts" if it passes; "the
         calibrated step loses" if its mean difference is below zero on both cuts with at most 2 of 10 seeds ahead on
         each; otherwise "not separated". T48b's test is reported and decides nothing.
         Reported: per cut every arm's mean board score; each calibrated arm's mean paired difference from m1_r_gmom,
         m1_r_gmom_x15 and m1_r_gmom_x2 with seeds ahead and per-seed range; every metric's mean skill difference and
         the part of each gain on the two de terms, mmd_u and the variogram (100 x the group's weight x the skill
         difference); per calibrated arm and cut, the gain over m1_r_gmom that a straight line through T44's committed
         m1_r_gmom_x15 and m1_r_gmom_x2 per-seed scores predicts at the calibrated fraction (marked extrapolated
         outside [0.3, 0.4]), and the residual (observed minus predicted) with its per-seed range; per cut every arm's
         mean detected panel genes and summed panel CP10k per cell; the calibration record; and, when the statement is
         "the calibrated step loses" or "not separated", the line that the reproduction and the calibration record are
         checked for a fault first and that, if none is found, the outcome stands as recorded and bears on the rule's
         target. The read puts nothing forward by itself.

NOT an adoption result: a labelled in-sample diagnostic, as tries T26, T42 and T44 are; Feng's E9.5 cells lie inside
these pseudo-tasks' external windows, the cuts are ones earlier T1 choices were made on, and L comes from T47's read of
the released E9.5 cells. Its result can only go to the team captain.

Limits (each sentence is in LIMITS and in the ledger row):
  In the units d is made in, M1's step of 0.2 reads 0.107 on E9.5's cells and 0.117 on E10.5 library 3's cells, that
  is 0.53 to 0.58 of its nominal size; a true fraction of Feng's change also reads under its nominal size in these
  units, so that ratio is not a shortfall against 1.0.
  At the test the rule would run once per T1 base (lib3 and mid) on that base's E10.5 cells, with fraction 0.4 and
  target 0.4 x L; by T47's numbers it would give about 1.43 times the nominal fraction, against about 1.56 times on
  these cuts, whose pool is E9.5's; the multiples assume a proportional delivered; the rule's multiple at the test is
  not known from these cuts.
  The weights are the pool's group counts; the build draws at the reference's own chamber shares, so its expected mix
  is the same, and one seed's draw moves the fraction by under 0.001 on synthetic pools (by up to 2% only if the mix
  were far from the pool's).
  The calibration draws on the fixed stream [48, 48]; a port to the test build needs its own row and a stream from the
  build seed.
  delivered(f) on one stream carries the draw's noise: the calibration record gives the sd of delivered at each answer
  over ten further streams, and the fraction is known only to about that sd divided by the slope of delivered; on a
  small pool that noise can make the rule refuse, so a test build that calls it unchanged is not assured an answer.
  T48a's fraction is 0.31 if delivered is proportional and lower if it is convex; where it falls below 0.3, T44's line
  is an extrapolation.
  T44's x1.5 arm passed on every seed; a pass is expected from T44 and is not evidence that the calibration target is
  right, nor about the fraction at the test.
  A loss or a not-separated outcome contradicts T44 at a neighbouring fraction; the reproduction and the calibration
  record are checked for a fault first; if none is found, the outcome stands as recorded and bears on the rule's
  target, not on T44's arms.
  The rule matches one d^2-weighted share; at the calibrated fraction genes above 10 CP10k move about 1.5 to 2 times
  the target along d while weakly expressed genes stay several times under; a gain may come from over-moving
  well-expressed genes, not from delivering the direction.
  delivered(f) rises with f; how far from proportional, and in which direction, depends on the pool (synthetic: about
  linear on one pool, delivered/f from 0.33 at 0.2 to 0.76 at 1.0 on another); the record's steps_by_arm share_over_f
  at 0.2, 0.3 and 0.4 gives the real shape.
  A larger step leaves total counts per cell unchanged but lowers detected genes per cell (synthetic: 19.5 genes of
  1,033 at the T48a answer against 5.0 at 0.2; a true fifth lowers them by 19.5 too); any gain on mmd_u or the
  variogram may travel through that sparsity, which the closed read-depth lever also moved, so a T48 result reopens
  nothing there and adoption for that reason needs the captain's named answer on that lever.
  L is T47's raw between-type slope at E9.5, a fixed input from another read; it is not computable from the test's
  base cells, the rule's fraction moves with it (L_corr 0.955 would give about 10 to 15% more, depending on the
  curvature of delivered), and the choice of raw over corrected is not tested here.
  Recomputing L on the test's base cells (E10.5) would need Feng's contrasts at a matching stage; against E9.5
  contrasts it would mix a day's change into the scale.
  The cuts share one reference pool and one calibration; their day holds more than 0.2L along d (T47's A 0.254), so
  the read favours any larger step; it says nothing about the pace from E10.5 to E12.5 or about the rule on the test's
  bases.
  T48b's fraction is about 0.37 to 0.38 if delivered is proportional and lower if it is convex, near T44's x2 arm
  (0.40); it reports nothing T44 did not.
  In these units even a true fraction of Feng's change reads under its nominal size (synthetic: a true fifth reads
  0.68 of nominal), so the target 0.2 x L is not the reading of the change M1 aims at; the rule's answer lies above
  the true-change-equivalent fraction (synthetic 0.39 and 0.44 against 0.34), and the excess falls on the
  well-expressed genes, which the nominal step already moves about as a true fifth would (synthetic: 0.105 against
  0.110 above 10 CP10k).
  The calibrated fraction carries a stream noise of about sd / slope (recorded); it is negligible beside the
  uncertainty of the target itself.
  A pass says that on the one E9.5 -> E10.5 pair a step at the calibrated fraction beats 0.2, as T44 already showed at
  0.3 and 0.4; it is no evidence that 0.2 x L is the right target, nor about the multiple the rule would give on E10.5
  bases with fraction 0.4.
  The cuts' day holds more than 0.2 x L along d, and the rule's multiple shrinks or grows with the base fraction
  according to the pool's curvature; the read cannot transfer its fraction or its gain to the test.
  L is a fixed input from T47's read of released E9.5 cells; it is target-free here but not computable from a test
  build's base cells; a port of the rule carries L as a frozen constant with its own row.
  passes: T44's direction holds at the calibrated fraction, nothing about the target or the test; loses or not
  separated: after a fault check, evidence against the rule's target; T48b reports nothing beyond T44's x2.
  T48 opens no closed lever; any use of a T48 gain that runs through sparsity waits for the captain's named answer on
  the read-depth lever.

Usage: on the DGX campaign checkout, repo root, both cuts at one commit, committing nothing until both task files exist,
on the machine T42 and T44 ran on, so that their arms reproduce.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/cal_read.py run --out t48_m1_cal \
      --task "T1/E9.5->E10.5:chamber"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/cal_read.py run --out t48_m1_cal \
      --task "T1/E9.5->E10.5"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/cal_read.py read --out t48_m1_cal
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
from scipy import sparse

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import GROUPS, REPO_ROOT, load_boards
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import share_read as shr  # noqa: E402  (sets the paths of dirlib, hygiene and n12)
import step_read as sr  # noqa: E402

er, gc, gv, feng_dirs, hygiene = sr.er, sr.gc, sr.gv, sr.feng_dirs, sr.hygiene
dl, n12 = shr.dl, sr.n12

ROW = "T48"
STEM = "t48"
CUT_FILES = gc.CUT_FILES
T42_TASKS = sr.T42_TASKS
T44_TASKS = RESULTS / "t44_m1_step" / "tasks"
T47_JSON = RESULTS / "t47_share_read" / "t47.json"
BASE = sr.BASE  # m1_r_gmom
FRACTION = sr.FRACTION  # 0.2
FRACTION_TOL = sr.FRACTION_TOL
T44_ARMS = tuple(sr.MULTIPLIERS)  # m1_r_gmom_x15, m1_r_gmom_x2
CAL_SCALED, CAL_PLAIN = "m1_r_gmom_calL", "m1_r_gmom_cal"  # T48a (target 0.2 x L), T48b (target 0.2)
CAL_ARMS = (CAL_SCALED, CAL_PLAIN)
ARMS = ("b0", BASE, *T44_ARMS, *CAL_ARMS)
CAL_ROOT = 48  # row T48's root: the calibration runs once per run, before any build seed
CAL_STREAM = 48  # [root, 48]: an index no build step uses
UPPER = 5.0  # the bracket's upper end, as a multiple of the fraction (the brief's [fraction, 5 x fraction])
REL_TOL = 1e-3  # on the share, relative to the target
MAX_STEPS = 60  # bisection halvings before the rule refuses
L_RANGE = (0.5, 1.5)  # T47's raw L must lie in this open interval
CAL_RANGE = (FRACTION, UPPER * FRACTION)  # (0.2, 1.0), open: a calibrated fraction at either end voids the read
REPRO_TOL = sr.REPRO_TOL  # 1e-9
REPRO_T42 = ("b0", BASE)
REPRO_T44 = T44_ARMS
LOSES_MAX_AHEAD = 2
PASSES = "the calibrated step passes on both cuts"
LOSES = "the calibrated step loses"
NOT_SEPARATED = "not separated"
T1_GROUP_WEIGHT = {members[0][0]: w for _, w, members in GROUPS["T1"]}  # metric -> its group's weight (one each)
DE_TERMS = ("de_score", "de_direction")
SPREAD_STREAMS = 10  # further streams [CAL_ROOT, CAL_STREAM, k], k = 1..10, for the noise of delivered (reported only)
CP_BINS = ((0.0, 0.1), (0.1, 1.0), (1.0, 10.0), (10.0, float("inf")))  # mean CP10k before the step (reported only)
CP_BIN_NAMES = ("[0, 0.1)", "[0.1, 1)", "[1, 10)", ">= 10")
FAULT_LINE = ("a loss or a not-separated outcome contradicts T44 at a neighbouring fraction; the reproduction and the "
              "calibration record are checked for a fault first; if none is found, the outcome stands as recorded and "
              "bears on the rule's target, not on T44's arms")
LIMITS = (
    "In the units d is made in, M1's step of 0.2 reads 0.107 on E9.5's cells and 0.117 on E10.5 library 3's cells, "
    "that is 0.53 to 0.58 of its nominal size; a true fraction of Feng's change also reads under its nominal size "
    "in these units, so that ratio is not a shortfall against 1.0.",
    "At the test the rule would run once per T1 base (lib3 and mid) on that base's E10.5 cells, with fraction 0.4 "
    "and target 0.4 x L; by T47's numbers it would give about 1.43 times the nominal fraction, against about 1.56 "
    "times on these cuts, whose pool is E9.5's; the multiples assume a proportional delivered; the rule's multiple "
    "at the test is not known from these cuts.",
    "The weights are the pool's group counts; the build draws at the reference's own chamber shares, so its "
    "expected mix is the same, and one seed's draw moves the fraction by under 0.001 on synthetic pools (by up to "
    "2% only if the mix were far from the pool's).",
    "The calibration draws on the fixed stream [48, 48]; a port to the test build needs its own row and a stream "
    "from the build seed.",
    "delivered(f) on one stream carries the draw's noise: the calibration record gives the sd of delivered at each "
    "answer over ten further streams, and the fraction is known only to about that sd divided by the slope of "
    "delivered; on a small pool that noise can make the rule refuse, so a test build that calls it unchanged is not"
    " assured an answer.",
    "T48a's fraction is 0.31 if delivered is proportional and lower if it is convex; where it falls below 0.3, "
    "T44's line is an extrapolation.",
    "T44's x1.5 arm passed on every seed; a pass is expected from T44 and is not evidence that the calibration "
    "target is right, nor about the fraction at the test.",
    "A loss or a not-separated outcome contradicts T44 at a neighbouring fraction; the reproduction and the "
    "calibration record are checked for a fault first; if none is found, the outcome stands as recorded and bears "
    "on the rule's target, not on T44's arms.",
    "The rule matches one d^2-weighted share; at the calibrated fraction genes above 10 CP10k move about 1.5 to 2 "
    "times the target along d while weakly expressed genes stay several times under; a gain may come from "
    "over-moving well-expressed genes, not from delivering the direction.",
    "delivered(f) rises with f; how far from proportional, and in which direction, depends on the pool (synthetic: "
    "about linear on one pool, delivered/f from 0.33 at 0.2 to 0.76 at 1.0 on another); the record's steps_by_arm "
    "share_over_f at 0.2, 0.3 and 0.4 gives the real shape.",
    "A larger step leaves total counts per cell unchanged but lowers detected genes per cell (synthetic: 19.5 genes"
    " of 1,033 at the T48a answer against 5.0 at 0.2; a true fifth lowers them by 19.5 too); any gain on mmd_u or "
    "the variogram may travel through that sparsity, which the closed read-depth lever also moved, so a T48 result "
    "reopens nothing there and adoption for that reason needs the captain's named answer on that lever.",
    "L is T47's raw between-type slope at E9.5, a fixed input from another read; it is not computable from the "
    "test's base cells, the rule's fraction moves with it (L_corr 0.955 would give about 10 to 15% more, depending "
    "on the curvature of delivered), and the choice of raw over corrected is not tested here.",
    "Recomputing L on the test's base cells (E10.5) would need Feng's contrasts at a matching stage; against E9.5 "
    "contrasts it would mix a day's change into the scale.",
    "The cuts share one reference pool and one calibration; their day holds more than 0.2L along d (T47's A 0.254),"
    " so the read favours any larger step; it says nothing about the pace from E10.5 to E12.5 or about the rule on "
    "the test's bases.",
    "T48b's fraction is about 0.37 to 0.38 if delivered is proportional and lower if it is convex, near T44's x2 "
    "arm (0.40); it reports nothing T44 did not.",
    "In these units even a true fraction of Feng's change reads under its nominal size (synthetic: a true fifth "
    "reads 0.68 of nominal), so the target 0.2 x L is not the reading of the change M1 aims at; the rule's answer "
    "lies above the true-change-equivalent fraction (synthetic 0.39 and 0.44 against 0.34), and the excess falls on"
    " the well-expressed genes, which the nominal step already moves about as a true fifth would (synthetic: 0.105 "
    "against 0.110 above 10 CP10k).",
    "The calibrated fraction carries a stream noise of about sd / slope (recorded); it is negligible beside the "
    "uncertainty of the target itself.",
    "A pass says that on the one E9.5 -> E10.5 pair a step at the calibrated fraction beats 0.2, as T44 already "
    "showed at 0.3 and 0.4; it is no evidence that 0.2 x L is the right target, nor about the multiple the rule "
    "would give on E10.5 bases with fraction 0.4.",
    "The cuts' day holds more than 0.2 x L along d, and the rule's multiple shrinks or grows with the base fraction"
    " according to the pool's curvature; the read cannot transfer its fraction or its gain to the test.",
    "L is a fixed input from T47's read of released E9.5 cells; it is target-free here but not computable from a "
    "test build's base cells; a port of the rule carries L as a frozen constant with its own row.",
    "passes: T44's direction holds at the calibrated fraction, nothing about the target or the test; loses or not "
    "separated: after a fault check, evidence against the rule's target; T48b reports nothing beyond T44's x2.",
    "T48 opens no closed lever; any use of a T48 gain that runs through sparsity waits for the captain's named "
    "answer on the read-depth lever.",
)
LABEL_T48 = ("try T48, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) at the fraction its own step needs, read back on the task's reference cells with T47's "
             "share code, to deliver 0.2 x L (T48a; L from T47's committed reading) and 0.2 (T48b, reported), on T44's "
             "route; Feng's E9.5 cells lie inside these pseudo-tasks' external windows, the cuts are ones earlier T1 "
             "choices were made on; the result can only go to the team captain")


# ---------------------------------------------------------------------------------------------------------------------
# The rule


def _mean_cp10k(X, rows: np.ndarray) -> np.ndarray:
    """share_read.calibration's mean CP10k (expm1 of .X in float64, averaged over the rows), without densifying."""
    S = X[np.flatnonzero(rows)]
    if sparse.issparse(S):
        S = sparse.csr_matrix(S, dtype=np.float64, copy=True)
        S.data = np.expm1(S.data)
        return np.asarray(S.mean(axis=0), dtype=np.float64).ravel()
    return np.expm1(np.asarray(S, np.float64)).mean(0)


def _settings(eps: float | None, min_cells: float | None) -> tuple[float, float]:
    if eps is None or min_cells is None:
        cfg = feng_dirs.load_config()
        eps = cfg["metrics"]["pseudocount_cp10k"] if eps is None else eps
        min_cells = cfg["metrics"]["min_cells_t1"] if min_cells is None else min_cells
    return float(eps), float(min_cells)


def prepare(pool, directions: Mapping[str, np.ndarray], groups: Mapping[str, str], held: Sequence[str],
            eps: float | None = None, min_cells: float | None = None) -> dict:
    """The rule's fixed part: the covered cells (every pool cell whose label maps to a direction group, as
    gmom_combo.moved_by_feng selects them; not the cells a build draws at the chamber shares), each group's cells, gene
    mask (a finite direction, not held), weight (its cells among the covered cells: the pool's group counts) and mean
    CP10k before the step. A group under min_cells keeps its weight but gives no share (combine then skips it)."""
    eps, min_cells = _settings(eps, min_cells)
    labels = pool.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    covered = np.array([groups.get(lab) in directions for lab in labels], dtype=bool)
    sub = pool[np.flatnonzero(covered)].copy()
    sl = labels[covered]
    held_set = set(map(str, held))
    not_held = np.array([g not in held_set for g in pool.var_names], dtype=bool)
    cell_group = np.array([groups.get(lab, "") for lab in sl], dtype=object)
    members, before, weights, masks = {}, {}, {}, {}
    for g in sorted(directions):
        m = cell_group == g
        masks[g] = np.isfinite(np.asarray(directions[g], np.float64)) & not_held
        weights[g] = float(m.sum())
        if m.sum() == 0 or m.sum() < min_cells:
            continue
        members[g] = m
        before[g] = {"cp": _mean_cp10k(sub.X, m)}
    return {"sub": sub, "directions": directions, "groups": groups, "held": list(held), "eps": eps,
            "min_cells": min_cells, "members": members, "before": before, "weights": weights, "masks": masks}


def delivered_on(prep: Mapping, f: float, seed: Sequence[int] = (CAL_ROOT, CAL_STREAM)) -> dict:
    """delivered(f) on a prepared pool: external_shift at f on the covered cells, read before -> after with
    share_read's share_stats and combine (weighted by the groups' cells)."""
    moved = _moved(prep, f, seed)
    shares = _shares(prep, moved, {g: prep["masks"][g] for g in prep["members"]})
    return {"share": shr.combine(shares, prep["weights"])[0], "groups": shares}


def _moved(prep: Mapping, f: float, seed: Sequence[int]):
    moved, _ = t1_scope.external_shift(prep["sub"], prep["directions"], f, prep["groups"], prep["held"],
                                       label_key=t1_scope.LABEL_KEY, seed=list(seed))
    return moved


def _shares(prep: Mapping, moved, masks: Mapping[str, np.ndarray]) -> dict[str, float]:
    """Per group share_stats' share of before -> after over the group's gene mask (NaN for a group without cells)."""
    shares = {g: float("nan") for g in prep["directions"]}
    for g, m in prep["members"].items():
        after = {"cp": _mean_cp10k(moved.X, m)}
        st = shr.share_stats(dl.lfc_ratio(prep["before"][g], after, prep["eps"]), prep["directions"][g], masks[g],
                             None)
        shares[g] = st["share"]
    return shares


def _cell_depth(C) -> dict[str, float]:
    C = sparse.csr_matrix(C)
    return {"detected_genes": float(np.diff(C.indptr).mean()) if C.shape[0] else float("nan"),
            "total_counts": float(np.asarray(C.sum(axis=1)).mean()) if C.shape[0] else float("nan")}


def step_record(prep: Mapping, f: float, seed: Sequence[int] = (CAL_ROOT, CAL_STREAM),
                panel_genes: Sequence[str] | None = None) -> dict:
    """Reported only, one external_shift at f on the rule's stream: the delivered share, and per bin of the group's mean
    CP10k before the step (CP_BINS) the combined share and the share over f; the covered cells' mean detected genes
    and total counts per cell before and after the step; given the board's scored genes, the combined share and the
    share over f over the panel genes M1 moves (the same step, no further draw)."""
    moved = _moved(prep, f, seed)
    rec = _step_record(prep, f, moved)
    if panel_genes is not None:
        in_panel = np.isin(np.asarray(prep["sub"].var_names, dtype=str), np.asarray(list(panel_genes), dtype=str))
        masks = {g: prep["masks"][g] & in_panel for g in prep["members"]}
        s = shr.combine(_shares(prep, moved, masks), prep["weights"])[0]
        rec["panel"] = {"share": s, "share_over_f": s / f, "genes": {g: int(m.sum()) for g, m in masks.items()}}
    return rec


def _step_record(prep: Mapping, f: float, moved) -> dict:
    share = shr.combine(_shares(prep, moved, {g: prep["masks"][g] for g in prep["members"]}), prep["weights"])[0]
    bins = {}
    for name, (lo, hi) in zip(CP_BIN_NAMES, CP_BINS, strict=True):
        masks = {g: prep["masks"][g] & (prep["before"][g]["cp"] >= lo) & (prep["before"][g]["cp"] < hi)
                 for g in prep["members"]}
        s = shr.combine(_shares(prep, moved, masks), prep["weights"])[0]
        bins[name] = {"share": s, "share_over_f": s / f,
                      "genes": {g: int(m.sum()) for g, m in masks.items()}}
    return {"fraction": float(f), "share": share, "share_over_f": share / f, "bins": bins,
            "before": _cell_depth(t1_scope.counts_of(prep["sub"])), "after": _cell_depth(moved.layers["counts"])}


def stream_spread(prep: Mapping, f: float, base_fraction: float, delivered_base: float, delivered_f: float,
                  seed: Sequence[int] = (CAL_ROOT, CAL_STREAM), n: int = SPREAD_STREAMS) -> dict:
    """Reported only: delivered at f on n further streams [*seed, k], k = 1..n, their sd, the slope of delivered
    between the base fraction and f on the rule's stream, and the fraction's sd (the sd over the slope)."""
    vals = [delivered_on(prep, f, (*seed, k))["share"] for k in range(1, n + 1)]
    slope = (delivered_f - delivered_base) / (f - base_fraction) if f != base_fraction else float("nan")
    sd = float(np.std(vals, ddof=1))
    return {"streams": [[*map(int, seed), k] for k in range(1, n + 1)], "values": [float(v) for v in vals],
            "sd": sd, "slope": float(slope), "fraction_sd": float(sd / slope) if np.isfinite(slope) and slope > 0
            else float("nan")}


def delivered(pool, directions, groups, held, f: float, *, eps: float | None = None, min_cells: float | None = None,
              seed: Sequence[int] = (CAL_ROOT, CAL_STREAM)) -> float:
    """delivered(f): the combined share M1's own step at f delivers on the pool's covered cells."""
    return delivered_on(prepare(pool, directions, groups, held, eps, min_cells), f, seed)["share"]


def _close(x: float, target: float) -> bool:
    return bool(np.isfinite(x) and abs(x - target) <= REL_TOL * abs(target))


def calibrated_fraction(pool, directions: Mapping[str, np.ndarray], groups: Mapping[str, str], held: Sequence[str],
                        fraction: float, target_share: float, *, eps: float | None = None,
                        min_cells: float | None = None, seed: Sequence[int] = (CAL_ROOT, CAL_STREAM)) -> dict:
    """The fraction f in [fraction, UPPER x fraction] at which delivered(f) = target_share (relative tolerance REL_TOL
    on the share), by bisection on one stream. Returns the fraction, delivered(fraction), the delivered share at the
    answer and the number of bisection steps (0 when the fraction itself delivers the target). Refuses (ValueError) a
    target outside the bracket's delivered range, a delivered share that does not rise across the bracket, and a
    bisection that does not reach the tolerance; it never extrapolates. Reads no file once eps and min_cells are
    given (else the config, for those two only)."""
    _check_args(fraction, target_share)
    return _calibrated(prepare(pool, directions, groups, held, eps, min_cells), fraction, target_share, seed)


def _check_args(fraction: float, target_share: float) -> None:
    if not (np.isfinite(fraction) and fraction > 0 and np.isfinite(target_share) and target_share != 0):
        raise ValueError(f"the rule needs a positive finite fraction and a finite non-zero target, got {fraction}, "
                         f"{target_share}")


def _calibrated(prep: Mapping, fraction: float, target_share: float, seed: Sequence[int]) -> dict:
    """calibrated_fraction on a prepared pool."""
    _check_args(fraction, target_share)
    lo, hi = float(fraction), float(UPPER * fraction)
    d_lo = delivered_on(prep, lo, seed)["share"]
    out = {"base_fraction": lo, "target": float(target_share), "bracket": [lo, hi], "stream": [int(s) for s in seed],
           "delivered_at_base": d_lo}
    if not np.isfinite(d_lo):
        raise ValueError(f"delivered({lo:g}) is not finite: no group gives a share")
    if _close(d_lo, target_share):
        return {**out, "fraction": lo, "delivered": d_lo, "steps": 0}
    d_hi = delivered_on(prep, hi, seed)["share"]
    if not (np.isfinite(d_hi) and d_hi > d_lo):
        raise ValueError(f"delivered does not rise across the bracket: {d_lo:.6g} at {lo:g}, {d_hi:.6g} at {hi:g}")
    if _close(d_hi, target_share):
        return {**out, "fraction": hi, "delivered": d_hi, "steps": 0}
    if not d_lo < target_share < d_hi:
        raise ValueError(f"the target {target_share:.6g} lies outside the bracket's delivered range [{d_lo:.6g}, "
                         f"{d_hi:.6g}] on [{lo:g}, {hi:g}]: the rule refuses, it does not extrapolate")
    for step in range(1, MAX_STEPS + 1):
        mid = (lo + hi) / 2.0
        d = delivered_on(prep, mid, seed)["share"]
        if _close(d, target_share):
            return {**out, "fraction": mid, "delivered": d, "steps": step}
        if d < target_share:
            lo = mid
        else:
            hi = mid
    raise ValueError(f"the bisection did not reach {REL_TOL:g} relative on the share in {MAX_STEPS} steps")


def calibrate(pool, directions, groups, held, fraction: float, L: float, t47_sha256: str, *,
              eps: float | None = None, min_cells: float | None = None,
              seed: Sequence[int] = (CAL_ROOT, CAL_STREAM), panel_genes: Sequence[str] | None = None) -> dict:
    """Both calibrations (T48a at fraction x L, T48b at fraction) on one prepared pool, and their record; reported
    beside them: each answer's stream noise (stream_spread) and each arm's step (step_record, with the share over the
    board's panel genes when they are given). The panel genes enter no fraction."""
    t0 = time.time()
    prep = prepare(pool, directions, groups, held, eps, min_cells)
    targets = {CAL_SCALED: float(fraction) * float(L), CAL_PLAIN: float(fraction)}
    res = {a: _calibrated(prep, fraction, t, seed) for a, t in targets.items()}
    fractions = {a: r["fraction"] for a, r in res.items()}
    rec = {"L": float(L), "t47": gv._rel(T47_JSON), "t47_sha256": t47_sha256, "base_fraction": float(fraction),
           "targets": targets, "fractions": fractions, "delivered_at_base": res[CAL_PLAIN]["delivered_at_base"],
           "delivered": {a: r["delivered"] for a, r in res.items()}, "steps": {a: r["steps"] for a, r in res.items()},
           "stream": res[CAL_PLAIN]["stream"], "bracket": res[CAL_PLAIN]["bracket"], "rel_tol": REL_TOL,
           "weights": dict(prep["weights"]), "seconds_rule": round(time.time() - t0, 1)}
    rec["noise"] = {a: stream_spread(prep, r["fraction"], float(fraction), r["delivered_at_base"], r["delivered"], seed)
                    for a, r in res.items()}
    rec["steps_by_arm"] = {a: step_record(prep, f, seed, panel_genes)
                           for a, f in arm_fractions(fraction, fractions).items()}
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


# ---------------------------------------------------------------------------------------------------------------------
# The run


def _num(x) -> bool:
    return isinstance(x, int | float) and not isinstance(x, bool) and bool(np.isfinite(x))


def t47_scale(path: Path = T47_JSON) -> tuple[float, str]:
    """T47's raw L and the file's sha256; the run is refused if the file is missing, void (unreadable, or its stamp's
    tree not clean) or its L is not finite and in L_RANGE."""
    if not Path(path).exists():
        sys.exit(f"{ROW} is refused: {gv._rel(path)} is missing (L comes from T47's committed reading)")
    try:
        rec = json.loads(Path(path).read_text())
        clean = rec["stamp"]["tree_clean"]
        L = rec["read"]["reading"]["L"]["value"]
    except (ValueError, KeyError, TypeError) as e:
        sys.exit(f"{ROW} is refused: {gv._rel(path)} is void (unreadable: {e!r})")
    if clean is not True:
        sys.exit(f"{ROW} is refused: {gv._rel(path)} is void (its tree was not clean)")
    if not _num(L) or not L_RANGE[0] < L < L_RANGE[1]:
        sys.exit(f"{ROW} is refused: T47's L {L!r} is not finite and in ({L_RANGE[0]:g}, {L_RANGE[1]:g})")
    return float(L), gv.sha256_file(Path(path))


def arm_fractions(fraction: float, cal_fractions: Mapping[str, float]) -> dict[str, float]:
    return {**sr.arm_fractions(fraction), **{a: float(cal_fractions[a]) for a in CAL_ARMS}}


def cal_arms(ctx: harness.BuildContext, directions: Mapping, groups: Mapping, fraction: float,
             held_m1: Sequence[str], held_factors: Sequence[str],
             cal_fractions: Mapping[str, float]) -> dict[str, Stage]:
    """T44's arms (step_read.step_arms) and, per calibrated arm, gmom_combo.m1_r_arms' m1_r_gmom at its fraction with
    every other argument the same object."""
    out = sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors)
    for a in CAL_ARMS:
        out[a] = gc.m1_r_arms(ctx, directions, groups, cal_fractions[a], held_m1, held_factors)[BASE]
    return out


def cell_stats(stage: Stage) -> dict[str, float]:
    """Reported only: a built prediction's mean detected panel genes per cell and mean summed panel CP10k per cell (a
    Stage holds log1p CP10k on the panel, not counts)."""
    X = np.asarray(stage.X, dtype=np.float64)
    return {"detected_genes": float((X > 0).sum(1).mean()), "panel_cp10k": float(np.expm1(X).sum(1).mean())}


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"{ROW} is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"{ROW} is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = re.sub(r"[^A-Za-z0-9.]+", "_", a.task).strip("_") + ".json"
    gv.refuse_existing(out / "tasks", (name,))
    L, t47_sha = t47_scale(T47_JSON)
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[a.task]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, _ = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                 (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref_stage.n)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    directions, ext = feng_dirs.feng_directions(list(pool.var_names))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    try:  # once per run, on a shallow copy of the reference pool (the inputs' obs untouched), before any build
        cal = calibrate(harness._shallow(pool), directions, groups, held_m1, fraction, L, t47_sha,
                        panel_genes=list(genes))
    except ValueError as e:
        sys.exit(f"{ROW} is refused: the calibration refuses ({e})")
    print(f"{t.key} calibration: " + ", ".join(f"{k} {v:.4f}" for k, v in cal["fractions"].items())
          + f" (delivered(0.2) {cal['delivered_at_base']:.4f}, {cal['seconds']} s)", flush=True)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = {"b0": b0(ctx), **cal_arms(ctx, directions, groups, fraction, held_m1, preregistered,
                                           cal["fractions"])}
        arms = er.score_arms(board, target, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "cells": {a: cell_stats(p) for a, p in
                                                                                 preds.items()},
                     "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res = {"label": LABEL_T48, "task": t.key, "fraction": fraction, "fractions": arm_fractions(fraction,
                                                                                              cal["fractions"]),
           "multipliers": dict(sr.MULTIPLIERS), "calibration": cal, "external": ext, "seeds": list(seeds),
           "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, REPO_ROOT / union["file"],
                                       feng_dirs.clean_stats_path(), n22, T47_JSON])}
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(rows: Sequence[Mapping], committed: Sequence[Mapping], arms: Sequence[str]) -> dict:
    """Per seed, `arms` against a committed run's scores; holds when every seed is there and each is within
    REPRO_TOL."""
    by_seed = {r["seed"]: r for r in committed}
    missing = [r["seed"] for r in rows if r["seed"] not in by_seed]
    worst = {a: 0.0 for a in arms}
    absent = []
    for r in rows:
        if r["seed"] not in by_seed:
            continue
        for arm in arms:
            if arm not in r["arms"] or arm not in by_seed[r["seed"]]["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - by_seed[r["seed"]]["arms"][arm]["score"]))
    holds = bool(committed) and not missing and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "missing_seeds": missing, "absent": absent, "holds": holds}


def calibration_checks(k: str, cal: Mapping | None) -> list[str]:
    """One file's calibration record: both fractions inside CAL_RANGE (open), the scaled one the smaller, and each
    delivered share at its answer equal to its target (FRACTION x L and FRACTION) within REL_TOL relative."""
    if not isinstance(cal, Mapping):
        return [f"{k}: no calibration record"]
    reasons = []
    fr, dv, L = cal.get("fractions") or {}, cal.get("delivered") or {}, cal.get("L")
    if not all(_num(fr.get(a)) for a in CAL_ARMS) or not _num(L):
        return [f"{k}: the calibration record lacks a fraction or L"]
    out = [a for a in CAL_ARMS if not CAL_RANGE[0] < fr[a] < CAL_RANGE[1]]
    if out:
        reasons.append(f"{k}: calibrated fraction outside ({CAL_RANGE[0]:g}, {CAL_RANGE[1]:g}) on {out}: "
                       f"{ {a: fr[a] for a in out} }")
    if not fr[CAL_SCALED] < fr[CAL_PLAIN]:
        reasons.append(f"{k}: the scaled calibrated fraction {fr[CAL_SCALED]:.6g} is not the smaller "
                       f"({fr[CAL_PLAIN]:.6g})")
    targets = {CAL_SCALED: FRACTION * L, CAL_PLAIN: FRACTION}
    off = [a for a in CAL_ARMS if not (_num(dv.get(a)) and abs(dv[a] - targets[a]) <= REL_TOL * abs(targets[a]))]
    if off:
        reasons.append(f"{k}: delivered share off its target by more than {REL_TOL:g} relative on {off} "
                       f"(delivered { {a: dv.get(a) for a in off} }, targets { {a: targets[a] for a in off} })")
    return reasons


def calibration_identity(results: Mapping[str, Mapping]) -> list[str]:
    """The cuts share one reference pool, so one calibration: every cut's record must give the same t47_sha256, the
    same L and the same two calibrated fractions (exactly)."""
    cals = {k: r.get("calibration") for k, r in results.items()}
    if len(cals) < 2 or not all(isinstance(c, Mapping) for c in cals.values()):
        return []

    def key(c):
        return {"t47_sha256": c.get("t47_sha256"), "L": c.get("L"),
                "fractions": {a: (c.get("fractions") or {}).get(a) for a in CAL_ARMS}}

    keys = {k: key(c) for k, c in cals.items()}
    first = next(iter(keys.values()))
    if all(v == first for v in keys.values()):
        return []
    return [f"the calibration differs between the cuts (t47_sha256, L or the calibrated fractions): {keys}"]


def validity(results: Mapping[str, Mapping], t42: Mapping[str, Mapping], t44: Mapping[str, Mapping]) -> dict:
    """The void conditions. results, t42, t44: cut -> task file (a missing cut is absent)."""
    reasons = []
    missing = [k for k in CUT_FILES if k not in results]
    if missing:
        reasons.append(f"missing cuts {missing}")
    for k, r in results.items():
        seeds = [row["seed"] for row in r["rows"]]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in r["rows"] for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if r.get("task") != k:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if not _num(r.get("fraction")) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        cal = r.get("calibration")
        cal_reasons = calibration_checks(k, cal)
        reasons += cal_reasons
        if isinstance(cal, Mapping) and all(_num((cal.get("fractions") or {}).get(a)) for a in CAL_ARMS):
            want = arm_fractions(FRACTION, cal["fractions"])
            got = r.get("fractions") or {}
            bad = [a for a, f in want.items() if not _num(got.get(a)) or abs(got[a] - f) > FRACTION_TOL]
            if bad:
                reasons.append(f"{k}: arm fractions {got}, not {want} (wrong fraction on {bad})")
    reasons += calibration_identity(results)
    commits = {r["stamp"]["commit"] for r in results.values()}
    unclean = sorted(k for k, r in results.items() if r["stamp"]["tree_clean"] is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the cuts ran at {len(commits)} commits")
    repro = {}
    for k, r in results.items():
        repro[k] = {"t42": reproduction(r["rows"], t42.get(k, {}).get("rows", []), REPRO_T42),
                    "t44": reproduction(r["rows"], t44.get(k, {}).get("rows", []), REPRO_T44)}
        if not repro[k]["t42"]["holds"]:
            reasons.append(f"{k}: b0 and m1_r_gmom do not reproduce T42's committed scores within {REPRO_TOL:g} per "
                           "seed")
        if not repro[k]["t44"]["holds"]:
            reasons.append(f"{k}: m1_r_gmom_x15 and m1_r_gmom_x2 do not reproduce T44's committed scores within "
                           f"{REPRO_TOL:g} per seed")
    complete = all(a in row["arms"] for r in results.values() for row in r["rows"] for a in ARMS)
    if results and complete:
        same = [a for a in CAL_ARMS
                if all(row["arms"][a] == row["arms"][BASE] for r in results.values() for row in r["rows"])]
        if same:
            reasons.append(f"calibrated arms score exactly as {BASE} on every seed of both cuts: {same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def gain_parts(skills: Mapping[str, float]) -> dict:
    """The part of a mean gain on each T1 metric (100 x its group's weight x the mean skill difference) and on the two
    de terms together."""
    parts = {m: 100.0 * T1_GROUP_WEIGHT[m] * float(v) for m, v in skills.items() if m in T1_GROUP_WEIGHT}
    parts["de_terms"] = sum(parts.get(m, 0.0) for m in DE_TERMS)
    return parts


def compare(results: Mapping[str, Mapping], base: str, change: str) -> dict:
    """step_read.compare (the adoption rule, paired by seed, both cuts; mean differences, seeds ahead, per-seed range,
    metric, group and skill differences) with each cut's gain parts added."""
    out = sr.compare(results, base, change)
    for cut in out["cuts"].values():
        cut["points_by_metric"] = gain_parts(cut["skills"])
    return out


def statement(comp: Mapping) -> str:
    if comp["passes"]:
        return PASSES
    if all(c["mean"] < 0 and c["ahead"] <= LOSES_MAX_AHEAD for c in comp["cuts"].values()):
        return LOSES
    return NOT_SEPARATED


def t44_line(results: Mapping[str, Mapping], t44: Mapping[str, Mapping]) -> dict:
    """Reported only: per calibrated arm and cut, the gain over m1_r_gmom that a straight line through T44's committed
    m1_r_gmom_x15 and m1_r_gmom_x2 per-seed scores predicts at the cut's calibrated fraction, the observed gain, and
    the residual (observed minus predicted) with its per-seed range; extrapolated outside [0.3, 0.4]."""
    fr44 = sr.arm_fractions(FRACTION)
    x15, x2 = T44_ARMS
    f15, f2 = fr44[x15], fr44[x2]
    out = {}
    for a in CAL_ARMS:
        out[a] = {}
        for k, r in results.items():
            f = float(r["calibration"]["fractions"][a])
            w = (f - f15) / (f2 - f15)
            committed = {row["seed"]: row["arms"] for row in t44[k]["rows"]}
            pred, obs = [], []
            for row in r["rows"]:
                c = committed[row["seed"]]
                base = row["arms"][BASE]["score"]
                pred.append(c[x15]["score"] + w * (c[x2]["score"] - c[x15]["score"]) - base)
                obs.append(row["arms"][a]["score"] - base)
            res = np.array(obs) - np.array(pred)
            out[a][k] = {"fraction": f, "weight_on_x2": float(w), "extrapolated": not f15 <= f <= f2,
                         "predicted": float(np.mean(pred)), "observed": float(np.mean(obs)),
                         "residual": float(res.mean()), "residual_range": [float(res.min()), float(res.max())]}
    return out


def cell_report(results: Mapping[str, Mapping]) -> dict:
    """Reported only: per cut and arm, the mean over seeds of the built cells' detected panel genes and summed panel
    CP10k per cell (rows without the record are skipped)."""
    out = {}
    for k, r in results.items():
        rows = [row["cells"] for row in r["rows"] if isinstance(row.get("cells"), Mapping)]
        if rows:
            out[k] = {a: {m: float(np.mean([c[a][m] for c in rows])) for m in ("detected_genes", "panel_cp10k")}
                      for a in ARMS if all(a in c for c in rows)}
    return out


def _cal_key(cal: Mapping) -> dict:
    return {k: v for k, v in cal.items() if k not in ("seconds", "seconds_rule")}


def t48_reading(results: Mapping[str, Mapping], t42: Mapping[str, Mapping], t44: Mapping[str, Mapping]) -> dict:
    v = validity(results, t42, t44)
    out = {"valid": v["valid"], "validity": v,
           "calibration": {k: r.get("calibration") for k, r in results.items()}}
    cals = [r.get("calibration") for r in results.values() if isinstance(r.get("calibration"), Mapping)]
    out["calibration_same_on_both_cuts"] = (len(cals) == len(CUT_FILES)
                                            and all(_cal_key(c) == _cal_key(cals[0]) for c in cals))
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    comps = {a: {other: compare(results, other, a) for other in (BASE, *T44_ARMS)} for a in CAL_ARMS}
    out["comparisons"] = comps
    out["outcome"] = statement(comps[CAL_SCALED][BASE])
    out["t48b_test"] = {"passes": comps[CAL_PLAIN][BASE]["passes"], "would_read": statement(comps[CAL_PLAIN][BASE]),
                        "decides": False}
    out["mean_scores"] = {k: {a: float(np.mean([row["arms"][a]["score"] for row in r["rows"]])) for a in ARMS}
                          for k, r in results.items()}
    out["t44_line"] = t44_line(results, t44)
    out["cells"] = cell_report(results)
    out["fault_line"] = FAULT_LINE if out["outcome"] in (LOSES, NOT_SEPARATED) else None
    return out


def t48_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"Statement (T48a, {CAL_SCALED} against {BASE}): {r['outcome']}"]
    if r.get("fault_line"):
        lines.append(f"Read as: {r['fault_line']}")
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, g in v["reproduction"].items():
        w42, w44 = g["t42"]["largest_difference"], g["t44"]["largest_difference"]
        lines.append(f"{k}: largest per-seed difference from committed scores: "
                     + ", ".join(f"{a} {x:.3g}" for a, x in {**w42, **w44}.items()) + f" (gate {REPRO_TOL:g})")
    for k, cal in r["calibration"].items():
        if not isinstance(cal, Mapping):
            lines.append(f"calibration {k}: none")
            continue
        lines.append(f"calibration {k}: L {cal.get('L')} (t47.json sha256 {cal.get('t47_sha256')}); delivered(0.2) "
                     f"{cal.get('delivered_at_base')}; fractions {cal.get('fractions')}; delivered at the answers "
                     f"{cal.get('delivered')} against targets {cal.get('targets')}; steps {cal.get('steps')}; stream "
                     f"{cal.get('stream')}; {cal.get('seconds')} s")
        for a, nz in (cal.get("noise") or {}).items():
            lines.append(f"    noise {a} (reported): sd of delivered over {len(nz.get('values', []))} further streams "
                         f"{nz.get('sd')}, slope {nz.get('slope')}, fraction sd {nz.get('fraction_sd')}")
        for a, st in (cal.get("steps_by_arm") or {}).items():
            bins = ", ".join(f"{b} {x['share_over_f']:.3f}" for b, x in st["bins"].items())
            lines.append(f"    step {a} at {st['fraction']:.4f} (reported): share {st['share']:.4f}, over f "
                         f"{st['share_over_f']:.3f}; over f by mean CP10k before: {bins}; detected genes per cell "
                         f"{st['before']['detected_genes']:.1f} -> {st['after']['detected_genes']:.1f}, total counts "
                         f"{st['before']['total_counts']:.0f} -> {st['after']['total_counts']:.0f}"
                         + (f"; panel genes: share {st['panel']['share']:.4f}, over f {st['panel']['share_over_f']:.3f}"
                            if isinstance(st.get("panel"), Mapping) else ""))
    lines.append(f"calibration identical on both cuts (one reference stage): {r['calibration_same_on_both_cuts']}")
    for k, ms in r.get("mean_scores", {}).items():
        lines.append(f"mean scores {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in ms.items()))
    if "t48b_test" in r:
        t = r["t48b_test"]
        lines.append(f"T48b ({CAL_PLAIN} against {BASE}), reported, decides nothing: "
                     f"{'passes' if t['passes'] else 'does not pass'} (would read: {t['would_read']})")
    for a, by in r.get("comparisons", {}).items():
        for other, comp in by.items():
            lines.append(f"{a} against {other}: {'PASSES' if comp['passes'] else 'does not pass'} the rule "
                         f"({comp['rule']['reason']}){'' if other == BASE else ' (reported)'}")
            for k, cut in comp["cuts"].items():
                sk = ", ".join(f"{m} {x:+.4f}" for m, x in cut["skills"].items())
                pp = ", ".join(f"{m} {x:+.2f}" for m, x in cut["points_by_metric"].items())
                lines.append(f"    {k}: {cut['mean']:+.2f} points, ahead {cut['ahead']} of {cut['of']}, per seed "
                             f"{cut['range'][0]:+.2f} to {cut['range'][1]:+.2f}; skills: {sk}; points: {pp}")
    for a, by in r.get("t44_line", {}).items():
        for k, x in by.items():
            lines.append(f"{a} {k}, T44's line at fraction {x['fraction']:.4f}"
                         f"{' (extrapolated)' if x['extrapolated'] else ''}: predicted {x['predicted']:+.2f}, observed "
                         f"{x['observed']:+.2f}, residual {x['residual']:+.2f} (per seed {x['residual_range'][0]:+.2f} "
                         f"to {x['residual_range'][1]:+.2f}) (reported)")
    for k, by in r.get("cells", {}).items():
        lines.append(f"cells {k} (reported): " + ", ".join(
            f"{a} {c['detected_genes']:.1f} genes, {c['panel_cp10k']:.0f} panel CP10k" for a, c in by.items()))
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines += [f"Limit: {s}" for s in LIMITS]
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t42_tasks: Path = T42_TASKS, t44_tasks: Path = T44_TASKS) -> dict:
    """Read the cuts present in out/tasks against T42's and T44's committed task files; write t48.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / f for k, f in CUT_FILES.items() if (out / "tasks" / f).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    committed = {}
    for name, folder in (("t42", t42_tasks), ("t44", t44_tasks)):
        committed[name] = {k: folder / f for k, f in CUT_FILES.items() if (folder / f).exists()}
    t42 = {k: json.loads(p.read_text()) for k, p in committed["t42"].items()}
    t44 = {k: json.loads(p.read_text()) for k, p in committed["t44"].items()}
    res = {"label": LABEL_T48, "reading": t48_reading(results, t42, t44),
           "stamps": {k: v["stamp"] for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *committed["t42"].values(), *committed["t44"].values()])
    gv.write_outputs(out, STEM, res, t48_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T48's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one cut of try T48")
    p.add_argument("--task", required=True, choices=list(CUT_FILES))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T48's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
