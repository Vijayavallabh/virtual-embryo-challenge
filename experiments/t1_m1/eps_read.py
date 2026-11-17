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
         that arm's own directions (lib_step_read.placebo_directions: the same two strata, the stream [build seed, 46]),
         reported only; and three expression-matched placebos, m1_r_gmom_sperm (fraction 0.2), m1_r_gmom_e01_x15_sperm
         and m1_r_gmom_e002_x15_sperm (fraction 0.3), the builds along each direction permuted within the same two
         strata and, inside each, within deciles of the gene's mean CP10k on the cut's own reference cells of the
         direction's group (all reference cells for a group with none), stream [build seed, 52]: statement 3 reads
         these (see there). The directions are made by feng_dirs.feng_directions with a deep copy of the loaded config
         whose metrics.pseudocount_cp10k is the arm's; no file and no shared object changes. The task file holds what
         try T46's holds (every arm's score, sha256 of the predicted matrix, scored zero share and detected scored genes
         per cell, per group the mean absolute log factor on the scored genes, the truth's row and contamination counts
         as checks never printed), every arm's fraction and pseudocount, each direction set's per-group sha256, per seed
         each placebo's permuted direction's sha256 per group and whether a new placebo's permutation of positions is
         its base placebo's (m1_r_gmom_perm's, m1_r_gmom_sperm's), and per direction set and group the squared size
         (sum of d squared over the moved genes) over the entry's and the shares of it on genes whose Feng mean is under
         0.5 CP10k at E9.5, at E14.5 and at either end (ratios). Count-space records, all target-free: per seed, arm
         and group the largest absolute log factor on any gene (external_shift's ShiftReport max_abs_log_f; no clip is
         applied on this route), the number of scored genes with |log f| above log 10, and the share of the group's
         reference-cell CP10k the factors move before the redraw and depth renormalisation (sum over genes of the
         reference mean CP10k x |f - 1| over the sum of the means; all genes and scored genes); per seed and arm the
         mean absolute difference of the scored genes' pseudobulk log1p CP10k from m1_r_gmom's prediction; per
         direction set and group the agreement of the half-A and half-B Feng directions (correlation over the moved
         genes and feng_dirs.eb_shrink's per-decile factors). Refused on an existing task file, unless the committed
         N22 passed, unless experiments/harness/TRIES.md holds row T52, and unless the loaded config's pseudocount is
         the entry's 0.5. A seed makes about 17 calls of gmom_combo.m1_r_arms on a chamber truth against try T46's 7:
         budget about 2.5 times T46's time per truth.
  read   T52's reading of the five truths' files (rule below, fixed here and in the row before the run).

The reading. Void unless the five truths hold the ten seeds at one commit on a clean tree (tree_clean exactly true); no
file records an anchor failure; every arm the truth should hold is on every row with its fraction and pseudocount; on
the four chamber truths b0, m1_r_gmom, _x15, _x2, _x3, _x4 and m1_r_gmom_perm equal, seed for seed, their committed
scores in experiments/t1_m1/results/t46_m1_lib_step/tasks/, and on `whole` b0, m1_r_gmom, _x15 and _x2 theirs in
experiments/t1_m1/results/t44_m1_step/tasks/ (the whole cut's file), each to 1e-9, and on the four chamber truths
the same arms' prediction sha256 equal T46's recorded ones seed for seed (a scorer-free check; T44's file holds none);
no new arm scores exactly as m1_r_gmom on every seed of a truth; no placebo's recorded prediction sha256 equals its real
arm's on every seed; the new arms' recorded directions differ from the entry's in every group. Also void, so that the
checks above cannot pass vacuously: a row without a sha256 string for every arm, or without the sparsity, log-factor,
count-space factor, pseudobulk and per-cut ruler read-back records. Also void (added after review): a direction set
whose own record from feng_dirs.feng_directions (the file's `external`, the config the direction code was actually
given) holds another pseudocount than the set's (0.5 for the entry's, 0.1 and 0.02 for the new ones): the per-arm
`pseudocounts` record is written from the script's constants and cannot catch a wrong config on the route.
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
scores better" if both lose; "no difference at the bar" if both are neither; otherwise "mixed: <comparison 1's label>
at the smaller size, <comparison 2's label> at the larger" (changed from the brief, which said "no difference at the
bar" for every other case: that wording was false when one comparison passes and the other loses, and misleading when
one passes and the other is neither), always with both comparisons' numbers. "The same delivered size" is the entry's
ruler along the entry's direction on library 3's cells, not a count-matched size: beside each comparison the read
prints both arms' count-space sizes on each cut (the share of reference-cell CP10k moved and the pseudobulk change
from m1_r_gmom; the scored genes' moved share split into the group's weakly and well expressed halves, at the median
reference level of the scored genes), so a reader sees how far the pair differs in counts moved and where. The
matching itself is checked where the pair is scored (added after review): per seed the run applies each real
directional arm's direction step alone (external_shift, M1's held genes, keep_depth, no clip, one fixed stream
[build seed, 53] for every arm) to the cut's own reference cells and records, per group with reference cells, the
slope (least squares with an intercept) of the realised log ratio of the group's mean CP10k, after over before, with
the entry's pseudocount 0.5, on the entry's direction over its moved genes. Beside each statement-2 pair the read
prints both arms' slopes on each cut and flags the pair on a cut where, in any group, the two arms differ by more than
half the smaller gap between the uniform arm and its neighbours on the ladder m1_r_gmom, _x15, _x2, _x3 (, _x4). A
flag on either cut adds "(not ruler-matched on the cuts' reference cells: ...)" to the statement, whose words then
claim only what was compared. The 0.02 arms have no uniform arm of the
same delivered size: each is reported against the two uniform arms that bracket it and carries no statement (_x4 is not
built on `whole`, so _e002_x2 against _x4 is reported on the chamber truths only).
Statement 3, per chamber truth and NN in 01, 002: D = (_eNN_x15 - m1_r_gmom) - (_eNN_x15_sperm - m1_r_gmom_sperm) per
seed with the single-task test: "the gain needs gene identity (not separated from the entry's own direction-specific
gain)" where D is ahead and _eNN_x15 is ahead of m1_r_gmom (changed from the brief's "the gain needs the direction's
genes", which read as if the decompressed genes carry the gain: D mixes the step from fraction 0.2 to 0.3 with the
change of direction, and a permuted placebo does more damage the larger the direction it permutes, so D can come out
ahead whenever the gate passes even when the truth moves along the entry's direction only; the text prints that
caveat on the same line as every statement-3 label: statement 3 checks that the gain is not generic movement of that
size, never which genes carry it); "no gain on this truth" where _eNN_x15 is not ahead of m1_r_gmom; otherwise "not
separated"; and "not read (permutation differs)" on a truth where any seed's _eNN_x15_sperm permutes positions
otherwise than m1_r_gmom_sperm (D's design rests on the two pairs sharing one permutation). Changed from the brief,
which read D with try T46's placebos (_eNN_x15_perm, m1_r_gmom_perm): T46's permutation keeps only the scored/other
strata, so it moves a decompressed direction's large factors, which sit on genes weak at one end, onto well-expressed
genes; the placebo pair then moves far more counts than the real pair and D is mostly the placebo's damage, ahead by
construction once the gate passes. The expression-matched placebo keeps each factor within its decile of
reference-cell expression, so the placebo moves about the counts the real direction moves (to within a decile). D with
T46's placebos is still computed and printed, marked reported and not read. Statement 3 says something about gene
identity only: it does not bear on decompression.
Statements 1 and 3 can both come out positive when decompression does no better than its delivered size (statement 1
inherits T44's uniform-inflation gain, since the new _x15 and _x2 arms deliver at least as much along the entry's
direction as m1_r_gmom; statement 3 inherits the entry's own direction-specific gain, T46); only statement 2 compares
decompression with size.
Reported: per truth every arm's mean board score and the arm with the highest mean among m1_r_gmom, the uniform arms and
the new arms (placebos and b0 left out; the first in the file's arm order on a tie); each new arm's paired difference
from m1_r_gmom and from the uniform arm of the nearest delivered size present on that truth (the pairing can differ
between the truths: _x4 is not built on `whole`; the text names each case), with seeds ahead, per-seed range and sd;
every metric's mean skill difference for the comparisons of statements 1 and 2 and the part on the two de terms and on
mmd_u and the variogram; each arm's scored zero share and detected scored genes per cell; the count-space records and
the half-split agreement above. _e01_x15's delivered size (0.2595) is not a quoted read-back: it is 1.5 x _e01's 0.173
on the read-back's linearity in the fraction, and is marked extrapolated wherever it is used. The read puts nothing
forward by itself.

Limits, repeated in the text (and in the row):
- One released day, in-sample, on cuts earlier choices were made on.
- A smaller pseudocount lets the ratios of genes weakly expressed at either end, which rest on few counts, move cells
  further, so a gain can be the larger factors on such genes and a loss their noise; up-moves on genes with zero counts
  are idle, and much of the count change can be the redraw and the depth renormalisation on well-expressed genes: the
  placebos and the uniform arms bound that and do not remove it.
- The delivered sizes quoted are read-backs on E10.5 library 3's cells, not on the cuts' own reference cells, and equal
  delivered size along the entry's direction does not make two arms equal in anything else (the smaller pseudocount's
  direction is 3 to 7 times larger in squared size): ruler-matched does not mean count-matched, and the count-space
  sizes printed beside statement 2 show how far each pair differs.
- Statement 1 against m1_r_gmom does not separate decompression from step size; only statement 2 compares them.
- Statements 1 and 3 can both be positive when decompression does no better than its delivered size; only statement 2
  bears on decompression.
- The expression-matched placebo matches the counts the real direction moves only to within a decile of reference-cell
  expression.
- All arms of a seed share one cell draw (the same cells, composition and count); they differ in the count redraw of
  the direction step and of g_mom's factors (the same streams with other factors), so the cell-draw noise is shared
  within a pair and the redraw noise is not.
- The pooled chamber truth and the whole cut come from the same released E10.5 file and the same reference, so the
  adoption rule's two cuts are not two independent replications (as in try T44).
- The pseudocounts 0.1 and 0.02 were chosen before the run from the read-back table and no other value is read.
- The pseudocounts and statement 2's matched pairs come from a target-free read-back on E10.5 library 3's cells, which
  are one of the four chamber truths: the match is in-sample for lib3.
- Statement 2's "neither" and "no difference at the bar" are not equivalence, and "mixed" supports neither side.
- Statement 3's D mixes the step from fraction 0.2 to 0.3 with the change of direction, and a permuted placebo does
  more damage the larger the direction it permutes, so D can come out ahead whenever the gate passes even when the
  truth moves along the entry's direction only: it checks that the gain is not generic movement of that size, never
  which genes carry it.
- The per-cut read-back applies the direction step alone (no R1, R2 or g_mom) to the cut's reference cells on one
  fixed stream, and its slope, like the quoted sizes, is the entry's ruler along the entry's direction: it is least
  sensitive on weakly expressed genes, where the smaller pseudocount's direction differs.
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
from scipy import sparse

from vec import core, harness, t1_scope
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
SPERM_ENTRY = "m1_r_gmom_sperm"  # the expression-matched placebo at the base fraction: statement 3's placebo base
SPERM_STREAM = 52  # [build seed, 52]: the expression-matched permutation's stream, which no other step uses
READBACK_STREAM = 53  # [build seed, 53]: the per-cut ruler read-back's redraw, which no other step uses
LEVEL_BINS = 10  # deciles of reference-cell expression within each stratum (the reviewer's proposal; not tuned)


def new_arm(tag: str, step: str) -> str:
    return f"{BASE}_{tag}{step}"


def new_placebo(tag: str) -> str:
    return f"{new_arm(tag, PLACEBO_STEP)}_perm"


def new_splacebo(tag: str) -> str:
    return f"{new_arm(tag, PLACEBO_STEP)}_sperm"


NEW_ARMS = {new_arm(t, s): (t, m) for t in EPS for s, m in STEPS.items()}  # arm -> (pseudocount tag, multiplier)
NEW_PLACEBOS = {new_placebo(t): new_arm(t, PLACEBO_STEP) for t in EPS}  # T46-style placebo -> its real arm (reported)
NEW_SPLACEBOS = {new_splacebo(t): new_arm(t, PLACEBO_STEP) for t in EPS}  # expression-matched placebo -> its real arm
PLACEBO_OF = {PERM_ENTRY: BASE, **NEW_PLACEBOS, SPERM_ENTRY: BASE, **NEW_SPLACEBOS}
BASE_PLACEBO = {**{p: PERM_ENTRY for p in NEW_PLACEBOS}, **{p: SPERM_ENTRY for p in NEW_SPLACEBOS}}
SAME_KEY = {PERM_ENTRY: "same_positions_as_m1_r_gmom_perm", SPERM_ENTRY: "same_positions_as_m1_r_gmom_sperm"}
REPRO_CHAMBER = ("b0", BASE, *UNIFORM, PERM_ENTRY)  # against try T46's committed files (score and sha256)
REPRO_WHOLE = ("b0", BASE, "m1_r_gmom_x15", "m1_r_gmom_x2")  # against try T44's committed whole-cut file
T46_TASKS = RESULTS / "t46_m1_lib_step" / "tasks"
T44_WHOLE = RESULTS / "t44_m1_step" / "tasks" / gc.CUT_FILES[gc.TASK_WHOLE]
# delivered size read back on E10.5 library 3's cells with the entry's ruler (the 2026-11-14 and 2026-11-17 target-free
# checks, as the brief quotes them); _e01_x15 is not quoted: its 0.2595 is extrapolated as 1.5 x _e01's 0.173 (the
# read-back is linear in the fraction) and is marked so wherever it is used
EXTRAPOLATED = ("m1_r_gmom_e01_x15",)
DELIVERED = {"m1_r_gmom_x15": 0.175, "m1_r_gmom_x2": 0.233, "m1_r_gmom_x3": 0.349, "m1_r_gmom_x4": 0.463,
             "m1_r_gmom_e01": 0.173, "m1_r_gmom_e01_x15": 1.5 * 0.173, "m1_r_gmom_e01_x2": 0.347,
             "m1_r_gmom_e002": 0.207, "m1_r_gmom_e002_x15": 0.312, "m1_r_gmom_e002_x2": 0.416}
MATCHED = (("m1_r_gmom_e01", "m1_r_gmom_x15"), ("m1_r_gmom_e01_x2", "m1_r_gmom_x3"))  # statement 2's two comparisons
BRACKETS = {"m1_r_gmom_e002": ("m1_r_gmom_x15", "m1_r_gmom_x2"), "m1_r_gmom_e002_x15": ("m1_r_gmom_x2", "m1_r_gmom_x3"),
            "m1_r_gmom_e002_x2": ("m1_r_gmom_x3", "m1_r_gmom_x4")}
LOW_CP10K = 0.5  # the shares of the squared size on genes whose Feng mean is under this (the entry's pseudocount)
BIG_LOG_F = float(np.log(10.0))  # scored genes whose |log factor| exceeds this are counted (no clip on the route)
PASSES, LOSES, NEITHER = "passes", "loses", "neither"
LEVEL = {"ahead": "ahead", "behind": "behind", "neither": "level"}
SMALLER_BETTER = "the smaller pseudocount scores better at the same delivered size"
UNIFORM_BETTER = "uniform inflation scores better"
NO_DIFFERENCE = "no difference at the bar"
MIXED = "mixed: {} at the smaller size, {} at the larger"
UNMATCHED = " (not ruler-matched on the cuts' reference cells: {})"
NEEDS_GENES = "the gain needs gene identity (not separated from the entry's own direction-specific gain)"
S3_CAVEAT = ("statement 3's D mixes the step from fraction 0.2 to 0.3 with the change of direction, and a permuted "
             "placebo does more damage the larger the direction it permutes, so D can come out ahead whenever the gate "
             "passes even when the truth moves along the entry's direction only: it checks that the gain is not "
             "generic movement of that size, never which genes carry it")
NO_GAIN = "no gain on this truth"
NOT_SEPARATED = "not separated"
NOT_READ_PERM = "not read (permutation differs)"
TASK_FILE = "{truth}.json"
LABEL_T52 = ("try T52, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with Feng's directions made with a pseudocount of 0.1 and 0.02 CP10k instead of 0.5, "
             "at 1, 1.5 and 2 times the entry's fraction, beside the step inflated uniformly, on tries T44's and T46's "
             "routes (the chamber cut pooled and library by library, and the whole cut); Feng's E9.5 cells lie inside "
             "these pseudo-tasks' external windows and the cuts are ones earlier T1 choices were made on; the result "
             "can only go to the team captain")
LIMITS = (
    "One released day, in-sample, on cuts earlier choices were made on.",
    "A smaller pseudocount lets the ratios of genes weakly expressed at either end, which rest on few counts, move "
    "cells further, so a gain can be the larger factors on such genes and a loss their noise; up-moves on genes with "
    "zero counts are idle, and much of the count change can be the redraw and the depth renormalisation on "
    "well-expressed genes: the placebos and the uniform arms bound that and do not remove it.",
    "The delivered sizes quoted are read-backs on E10.5 library 3's cells, not on the cuts' own reference cells, and "
    "equal delivered size along the entry's direction does not make two arms equal in anything else (the smaller "
    "pseudocount's direction is 3 to 7 times larger in squared size): ruler-matched does not mean count-matched, and "
    "the count-space sizes printed beside statement 2 show how far each pair differs.",
    "Statement 1 against m1_r_gmom does not separate decompression from step size; only statement 2 compares them.",
    "Statements 1 and 3 can both be positive when decompression does no better than its delivered size; only "
    "statement 2 bears on decompression.",
    "The expression-matched placebo matches the counts the real direction moves only to within a decile of "
    "reference-cell expression.",
    "All arms of a seed share one cell draw (the same cells, composition and count); they differ in the count redraw "
    "of the direction step and of g_mom's factors (the same streams with other factors), so the cell-draw noise is "
    "shared within a pair and the redraw noise is not.",
    "The pooled chamber truth and the whole cut come from the same released E10.5 file and the same reference, so the "
    "adoption rule's two cuts are not two independent replications (as in try T44).",
    "The pseudocounts 0.1 and 0.02 were chosen before the run from the read-back table and no other value is read.",
    "The pseudocounts and statement 2's matched pairs come from a target-free read-back on E10.5 library 3's cells, "
    "which are one of the four chamber truths: the match is in-sample for lib3.",
    "Statement 2's \"neither\" and \"no difference at the bar\" are not equivalence, and \"mixed\" supports neither "
    "side.",
    S3_CAVEAT[0].upper() + S3_CAVEAT[1:] + ".",
    "The per-cut read-back applies the direction step alone (no R1, R2 or g_mom) to the cut's reference cells on one "
    "fixed stream, and its slope, like the quoted sizes, is the entry's ruler along the entry's direction: it is least "
    "sensitive on weakly expressed genes, where the smaller pseudocount's direction differs.",
    "Nothing here says how far any hidden truth moves.",
)


def uniform_of(truth: str) -> dict[str, float]:
    return {a: m for a, m in UNIFORM.items() if truth != WHOLE or a in UNIFORM_WHOLE}


def arms_of(truth: str) -> tuple[str, ...]:
    """The arms a truth's file holds, in the order the run builds them."""
    chamber = truth != WHOLE
    return ("b0", BASE, *uniform_of(truth), *((PERM_ENTRY,) if chamber else ()), *NEW_ARMS,
            *((*NEW_PLACEBOS, SPERM_ENTRY, *NEW_SPLACEBOS) if chamber else ()))


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
        out[SPERM_ENTRY] = {"fraction": float(fraction), "pseudocount": ENTRY_EPS}
        out.update({p: {"fraction": float(fraction) * STEPS[PLACEBO_STEP], "pseudocount": EPS[NEW_ARMS[a][0]]}
                    for p, a in NEW_SPLACEBOS.items()})
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


def feng_means(genes: Sequence[str], record: Mapping, end: int,
               stats_path: Path | None = None) -> dict[str, np.ndarray]:
    """Per group of a direction record, the group's Feng mean CP10k at the bracket's start (end=0, E9.5) or end (end=1,
    E14.5) on `genes` (NaN where unmapped), from the same statistics file and labels the directions used."""
    path = Path(stats_path) if stats_path is not None else feng_dirs.clean_stats_path()
    feng_dirs.fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    idx = dl.sym_index(st["genes"], list(genes))
    stage = record["bracket"][end]
    return {g: dl.on_t1(dl.feng_sel(st, stage, c["feng_labels"])["cp"], idx) for g, c in record["groups"].items()}


def start_means(genes: Sequence[str], record: Mapping, stats_path: Path | None = None) -> dict[str, np.ndarray]:
    return feng_means(genes, record, 0, stats_path)


def end_means(genes: Sequence[str], record: Mapping, stats_path: Path | None = None) -> dict[str, np.ndarray]:
    return feng_means(genes, record, 1, stats_path)


def moved_mask(v, genes: Sequence[str], held: Sequence[str]) -> np.ndarray:
    """The genes the step moves: a finite direction value and not held (external_shift's factor is 1 elsewhere)."""
    held_set = set(map(str, held))
    return np.isfinite(np.asarray(v, dtype=np.float64)) & np.array([str(g) not in held_set for g in genes])


def direction_stats(sets: Mapping[str, Mapping[str, np.ndarray]], means: Mapping[str, np.ndarray],
                    genes: Sequence[str], held: Sequence[str], end: Mapping[str, np.ndarray] | None = None) -> dict:
    """Per direction set (`entry` and the new ones) and group: the squared size (sum of d squared over the moved genes)
    over the entry's, and the share of it on genes whose Feng start mean (`means`, E9.5) is under LOW_CP10K; with `end`
    (the E14.5 means) also the shares on genes under LOW_CP10K at the end and at either end. Ratios only. These are
    direction-space numbers: they do not say where the predictions' counts change (the count-space records do)."""
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

            def share(low, v=v, sq=sq):
                return float((v[low] ** 2).sum()) / sq if sq > 0 else None

            low_start = m & (np.asarray(means[g]) < LOW_CP10K)
            rec[g] = {"squared_size_over_entry": sq / sq_entry if sq_entry > 0 else None,
                      "share_on_low_start_genes": share(low_start)}
            if end is not None:
                low_end = m & (np.asarray(end[g]) < LOW_CP10K)
                rec[g]["share_on_low_end_genes"] = share(low_end)
                rec[g]["share_on_low_either_end_genes"] = share(low_start | low_end)
        out[tag] = rec
    return out


def positions_of(directions: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Each group's vector of its own positions, NaN where the direction has none: permuted as a direction is, it
    shows the permutation of positions itself."""
    return {g: np.where(np.isfinite(np.asarray(v, dtype=np.float64)), np.arange(len(v), dtype=np.float64), np.nan)
            for g, v in directions.items()}


def same_positions(ctx: harness.BuildContext, a: Mapping, b: Mapping, held: Sequence[str], permute=None) -> bool:
    """Whether a placebo's permutation (try T46's, lib_step_read.placebo_directions, unless `permute` is given as
    permute(ctx, directions, held)) moves the positions of directions a and b the same way for this seed."""
    permute = lsr.placebo_directions if permute is None else permute
    pa = permute(ctx, positions_of(a), held)
    pb = permute(ctx, positions_of(b), held)
    return set(pa) == set(pb) and all(np.array_equal(pa[g], pb[g], equal_nan=True) for g in pa)


def reference_levels(ctx: harness.BuildContext, groups: Mapping[str, str],
                     directions: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Per direction group, each gene's mean CP10k over the reference role's cells whose label maps to that group (the
    cells the step moves), on the reference role's gene order; over all the reference role's cells for a group with
    none there (the step then moves nothing). Target-free: the cut's own reference cells only."""
    ref = ctx.inputs[ctx.task.reference_role]
    C = sparse.csr_matrix(core.counts_of(ref), dtype=np.float64)
    depth = np.asarray(C.sum(axis=1)).ravel()
    cp = sparse.diags(np.where(depth > 0, 1e4 / np.where(depth > 0, depth, 1.0), 0.0)) @ C
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cell_group = np.array([groups.get(lab, "") for lab in labels], dtype=object)
    every = np.asarray(cp.mean(axis=0)).ravel()
    out = {}
    for g in sorted(directions):
        rows = np.flatnonzero(cell_group == g)
        out[g] = np.asarray(cp[rows].mean(axis=0)).ravel() if rows.size else every.copy()
    return out


def stratified_permuted_directions(directions: Mapping[str, np.ndarray], genes: Sequence[str], held: Sequence[str],
                                   scored: Sequence[str], seed, levels: Mapping[str, np.ndarray]) -> dict:
    """Try T46's placebo with expression strata added: each group's direction permuted within the moved genes (finite,
    not held) of each of try T46's two strata (scored, other) and, inside each, within LEVEL_BINS equal-count bins of
    the genes ordered by `levels[group]` (stable order; reference-cell mean CP10k). NaN and held positions stay. The
    bins depend on the moved-gene set and the levels only, not on the direction's values, so two directions with the
    same moved genes are permuted the same way. Groups in sorted order, scored before other, bins lowest first, one
    generator default_rng([*seed, SPERM_STREAM])."""
    held_set, scored_set = set(map(str, held)), set(map(str, scored))
    movable = np.array([str(g) not in held_set for g in genes], dtype=bool)
    on_board = np.array([str(g) in scored_set for g in genes], dtype=bool)
    rng = np.random.default_rng([*np.atleast_1d(seed).astype(int).tolist(), SPERM_STREAM])
    out = {}
    for g in sorted(directions):
        v = np.array(directions[g], dtype=np.float64, copy=True)
        lev = np.asarray(levels[g], dtype=np.float64)
        moved = np.isfinite(v) & movable
        for stratum in (moved & on_board, moved & ~on_board):
            idx = np.flatnonzero(stratum)
            for part in np.array_split(idx[np.argsort(lev[idx], kind="stable")], LEVEL_BINS):
                v[part] = v[part][rng.permutation(len(part))]
        out[g] = v
    return out


def splacebo_directions(ctx: harness.BuildContext, directions: Mapping, held_m1: Sequence[str],
                        levels: Mapping[str, np.ndarray]) -> dict:
    """The expression-matched permuted directions on the reference role's gene order, M1's held genes, the board's
    scored genes and the build seed (as lib_step_read.placebo_directions)."""
    genes = list(ctx.inputs[ctx.task.reference_role].var_names)
    return stratified_permuted_directions(directions, genes, held_m1, list(ctx.board.genes), ctx.seed, levels)


def mean_abs_log_factor(v, genes: Sequence[str], held: Sequence[str], scored: Sequence[str], fraction: float) -> float:
    """The mean absolute log factor over the scored genes at `fraction` (lib_step_read.scored_log_factors' measure)."""
    return lsr.scored_log_factors({"g": v}, {"g": v}, genes, held, scored, fraction)["g"][BASE]


def factor_stats(v, genes: Sequence[str], held: Sequence[str], scored: Sequence[str], fraction: float,
                 level) -> dict:
    """Count-space records of one group's step, target-free: log f = fraction x d (0 where NaN and on held genes, no
    clip: external_shift's own formula). max_abs_log_f: the largest |log f| on any gene (ShiftReport's);
    scored_over_log10: scored genes with |log f| above log 10; cp10k_moved_share (and _scored): the sum over genes
    (scored genes) of the group's reference mean CP10k x |f - 1|, over the sum of the means, i.e. the share of the
    group's mean profile the factors move before the redraw and the depth renormalisation; cp10k_moved_share_scored_weak
    and _strong: the scored part split at the median reference level of the scored genes (weak: below it; strong: the
    rest), over the same sum, so a gain confined to weakly expressed genes is visible in count space."""
    held_set, scored_set = set(map(str, held)), set(map(str, scored))
    on_board = np.array([str(g) in scored_set for g in genes], dtype=bool)
    held_mask = np.array([str(g) in held_set for g in genes], dtype=bool)
    logf = float(fraction) * np.nan_to_num(np.asarray(v, dtype=np.float64), nan=0.0)
    logf[held_mask] = 0.0
    lev = np.asarray(level, dtype=np.float64)
    moved = lev * np.abs(np.expm1(logf))
    total = float(lev.sum())
    weak = on_board & (lev < np.median(lev[on_board])) if on_board.any() else on_board
    strong = on_board & ~weak

    def share(mask):
        return float(moved[mask].sum()) / total if total > 0 else None

    return {"max_abs_log_f": float(np.abs(logf).max()) if logf.size else 0.0,
            "scored_over_log10": int(((np.abs(logf) > BIG_LOG_F) & on_board).sum()),
            "cp10k_moved_share": share(np.ones(len(lev), dtype=bool)), "cp10k_moved_share_scored": share(on_board),
            "cp10k_moved_share_scored_weak": share(weak), "cp10k_moved_share_scored_strong": share(strong)}


def pseudobulk(pred: Stage) -> np.ndarray:
    """The prediction's pseudobulk on the scored genes: log1p of the mean CP10k over its cells (X is log1p CP10k)."""
    X = np.asarray(pred.X.toarray() if hasattr(pred.X, "toarray") else pred.X, dtype=np.float64)
    return np.log1p(np.expm1(X).mean(axis=0))


def pseudobulk_changes(preds: Mapping[str, Stage], base: str = BASE) -> dict[str, float]:
    """Per arm, the mean absolute difference of its scored genes' pseudobulk log1p CP10k from `base`'s prediction."""
    ref = pseudobulk(preds[base])
    return {a: float(np.abs(pseudobulk(p) - ref).mean()) for a, p in preds.items()}


def half_agreement(genes: Sequence[str], held: Sequence[str], pseudocounts: Mapping[str, float],
                   stats_path: Path | None = None) -> dict:
    """Per direction set (tag -> pseudocount) and group, how far Feng's own half-A and half-B directions agree at that
    pseudocount (feng_dirs.direction_halves with the config copy of eps_config): their correlation over the moved genes
    with a value in both halves, and feng_dirs.eb_shrink's per-decile factors (signal over signal plus noise, lowest
    expression first). Target-free. A set whose halves cannot be read records the reason."""
    path = Path(stats_path) if stats_path is not None else feng_dirs.clean_stats_path()
    feng_dirs.fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    held_set = set(map(str, held))
    movable = np.array([str(g) not in held_set for g in genes], dtype=bool)
    out = {}
    for tag, eps in pseudocounts.items():
        cfg = eps_config(eps)
        dirs, _, _ = feng_dirs.directions_from_stats(st, list(genes), cfg)
        try:
            halves = feng_dirs.direction_halves(st, list(genes), cfg)
        except ValueError as e:
            out[tag] = {"error": str(e)}
            continue
        rec = {}
        for g in sorted(set(dirs) & set(halves)):
            h = halves[g]
            ok = np.isfinite(h["d_a"]) & np.isfinite(h["d_b"]) & movable
            corr = (float(np.corrcoef(h["d_a"][ok], h["d_b"][ok])[0, 1])
                    if ok.sum() > 1 and np.std(h["d_a"][ok]) > 0 and np.std(h["d_b"][ok]) > 0 else None)
            _, factors = feng_dirs.eb_shrink(dirs[g], h["d_a"], h["d_b"], h["mean_cp"])
            rec[g] = {"correlation": corr, "eb_shrink_factors": [float(x) for x in factors]}
        out[tag] = rec
    return out


def ruler_arms(truth: str, directions: Mapping, eps_dirs: Mapping[str, Mapping],
               fraction: float) -> dict[str, tuple[Mapping, float]]:
    """The real directional arms of a truth (m1_r_gmom, the uniform arms, the new arms; no placebo) with the directions
    and the very fraction each is built with (step_read.arm_fractions)."""
    out = {a: (directions, f) for a, f in sr.arm_fractions(fraction, uniform_of(truth)).items()}
    for tag, dirs in eps_dirs.items():
        fr = sr.arm_fractions(fraction, {new_arm(tag, s): m for s, m in STEPS.items() if s})
        out.update({(new_arm(tag, "") if a == BASE else a): (dirs, f) for a, f in fr.items()})
    return out


def ruler_readback(ctx: harness.BuildContext, groups: Mapping[str, str], entry: Mapping[str, np.ndarray],
                   arms: Mapping[str, tuple[Mapping, float]], held: Sequence[str]) -> dict[str, dict[str, float]]:
    """Target-free, on the cut's own reference cells: per arm, its direction step alone (t1_scope.external_shift, the
    held genes, keep_depth, no clip, stream [build seed, READBACK_STREAM] for every arm); per direction group with
    reference cells, the least-squares slope (with an intercept) of the realised log ratio of the group's mean CP10k,
    after over before, with the entry's pseudocount, on the entry's direction over its moved genes: the entry's ruler
    along the entry's direction, read where the pair is scored. group -> arm -> slope (None where the moved genes'
    direction has no spread)."""
    ref = ctx.inputs[ctx.task.reference_role]
    genes = list(ref.var_names)
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    cell_group = np.array([groups.get(lab, "") for lab in labels], dtype=object)
    C0 = sparse.csr_matrix(core.counts_of(ref))
    seed = [*np.atleast_1d(ctx.seed).astype(int).tolist(), READBACK_STREAM]
    present = [g for g in sorted(entry) if (cell_group == g).any()]
    before = {g: core.mean_cp10k(C0[np.flatnonzero(cell_group == g)]) for g in present}
    out = {g: {} for g in present}
    for a, (dirs, f) in arms.items():
        moved_ref, _ = t1_scope.external_shift(ref, dirs, f, groups, held, seed=seed)
        C1 = sparse.csr_matrix(core.counts_of(moved_ref))
        for g in present:
            d = np.asarray(entry[g], dtype=np.float64)
            m = moved_mask(d, genes, held)
            after = core.mean_cp10k(C1[np.flatnonzero(cell_group == g)])
            r = np.log((after + ENTRY_EPS) / (before[g] + ENTRY_EPS))
            x, y = d[m] - d[m].mean(), r[m] - r[m].mean()
            den = float((x * x).sum())
            out[g][a] = float((x * y).sum()) / den if m.sum() > 1 and den > 0 else None
    return out


def readback_means(r: Mapping) -> dict[str, dict[str, float | None]]:
    """Per group and arm, the mean over seeds of the per-cut read-back (None where a seed has none)."""
    rows = r.get("rows") or []
    first = (rows[0].get("ruler_readback") or {}) if rows else {}
    out = {}
    for g in sorted(first):
        out[g] = {}
        for a in first[g] or {}:
            xs = [((row.get("ruler_readback") or {}).get(g) or {}).get(a) for row in rows]
            ok = all(isinstance(x, int | float) and not isinstance(x, bool) for x in xs)
            out[g][a] = float(np.mean(xs)) if ok else None
    return out


def ruler_match(r: Mapping, truth: str, arm: str, uniform: str) -> dict:
    """Whether a statement-2 pair is ruler-matched on a cut's own reference cells: per group the two arms' mean
    read-backs and half the smaller gap between the uniform arm and its neighbours on the ladder m1_r_gmom, then the
    uniform arms in multiplier order (those built on this truth); the pair is flagged when, in any group, the two arms
    differ by more than that half gap, or a value is missing, or no group has a read-back."""
    means = readback_means(r)
    ladder = [BASE, *sorted(uniform_of(truth), key=UNIFORM.get)]
    i = ladder.index(uniform)
    neighbours = [ladder[j] for j in (i - 1, i + 1) if 0 <= j < len(ladder)]
    groups, flagged = {}, []
    for g, by in means.items():
        vals = [by.get(x) for x in (arm, uniform, *neighbours)]
        if any(v is None for v in vals):
            groups[g] = {"arm": vals[0], "uniform": vals[1], "half_gap": None, "flagged": True}
            flagged.append(g)
            continue
        half = 0.5 * min(abs(vals[1] - v) for v in vals[2:])
        bad = abs(vals[0] - vals[1]) > half
        groups[g] = {"arm": vals[0], "uniform": vals[1], "half_gap": half, "flagged": bool(bad)}
        if bad:
            flagged.append(g)
    return {"arm": arm, "uniform": uniform, "neighbours": neighbours, "groups": groups, "flagged_groups": flagged,
            "matched": bool(groups) and not flagged}


# ---------------------------------------------------------------------------------------------------------------------
# The run


def build_seed(ctx: harness.BuildContext, b0, truth: str, directions: Mapping, eps_dirs: Mapping[str, Mapping],
               groups: Mapping, fraction: float, held_m1: Sequence[str],
               held_factors: Sequence[str]) -> tuple[dict[str, Stage], dict, dict, dict]:
    """One seed's arms in arms_of(truth)'s order; per group and directional arm the mean absolute log factor on the
    scored genes; the placebos' record (each permuted direction's sha256 per group, and whether a new placebo's
    permutation of positions is its base placebo's: m1_r_gmom_perm's for the _perm placebos, m1_r_gmom_sperm's for the
    _sperm ones); and per group and directional arm the count-space factor records (factor_stats)."""
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
    placebo = {"sha256": {}, SAME_KEY[PERM_ENTRY]: {}, SAME_KEY[SPERM_ENTRY]: {}}
    levels = reference_levels(ctx, groups, directions)
    if chamber:
        placebo["sha256"][PERM_ENTRY] = directions_sha256(perms[PERM_ENTRY])
        f = float(fraction) * STEPS[PLACEBO_STEP]
        for tag, dirs in eps_dirs.items():
            p = new_placebo(tag)
            perms[p] = lsr.placebo_directions(ctx, dirs, held_m1)
            preds[p] = gc.m1_r_arms(ctx, perms[p], groups, f, held_m1, held_factors)[BASE]
            used[p] = (perms[p], f)
            placebo["sha256"][p] = directions_sha256(perms[p])
            placebo[SAME_KEY[PERM_ENTRY]][p] = same_positions(ctx, directions, dirs, held_m1)

        def spermute(c, d, h):
            return splacebo_directions(c, d, h, levels)

        perms[SPERM_ENTRY] = spermute(ctx, directions, held_m1)
        preds[SPERM_ENTRY] = gc.m1_r_arms(ctx, perms[SPERM_ENTRY], groups, fraction, held_m1, held_factors)[BASE]
        used[SPERM_ENTRY] = (perms[SPERM_ENTRY], float(fraction))
        placebo["sha256"][SPERM_ENTRY] = directions_sha256(perms[SPERM_ENTRY])
        for tag, dirs in eps_dirs.items():
            p = new_splacebo(tag)
            perms[p] = spermute(ctx, dirs, held_m1)
            preds[p] = gc.m1_r_arms(ctx, perms[p], groups, f, held_m1, held_factors)[BASE]
            used[p] = (perms[p], f)
            placebo["sha256"][p] = directions_sha256(perms[p])
            placebo[SAME_KEY[SPERM_ENTRY]][p] = same_positions(ctx, directions, dirs, held_m1, permute=spermute)
    genes = list(ctx.inputs[ctx.task.reference_role].var_names)
    scored = list(ctx.board.genes)
    log_factors = {g: {a: mean_abs_log_factor(d[g], genes, held_m1, scored, f) for a, (d, f) in used.items()}
                   for g in sorted(directions)}
    factors = {g: {a: factor_stats(d[g], genes, held_m1, scored, f, levels[g]) for a, (d, f) in used.items()}
               for g in sorted(directions)}
    return {a: preds[a] for a in arms_of(truth)}, log_factors, placebo, factors


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
    ends = end_means(pool_genes, ext)
    halves = half_agreement(pool_genes, held_m1, {"entry": ENTRY_EPS, **EPS})
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
        preds, log_factors, placebo, factors = build_seed(ctx, b0, a.truth, directions, eps_dirs, groups, fraction,
                                                          held_m1, preregistered)
        readback = ruler_readback(ctx, groups, directions, ruler_arms(a.truth, directions, eps_dirs, fraction), held_m1)
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms,
                     "sha256": {k: lsr.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: lsr.sparsity(p) for k, p in preds.items()}, "scored_log_factors": log_factors,
                     "factor_stats": factors, "pseudobulk_change": pseudobulk_changes(preds),
                     "ruler_readback": readback, "placebo": placebo, "seconds": round(time.time() - t0, 1)})
        print(f"{a.truth} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    spec = arm_spec(a.truth, fraction)
    res.update({"fraction": fraction, "fractions": {k: v["fraction"] for k, v in spec.items()},
                "pseudocounts": {k: v["pseudocount"] for k, v in spec.items()}, "directions": dir_record,
                "direction_stats": direction_stats(sets, means, pool_genes, held_m1, ends), "half_agreement": halves,
                "placebo": {"arms": {p: r for p, r in PLACEBO_OF.items() if p in spec},
                            "base_placebo": {p: b for p, b in BASE_PLACEBO.items() if p in spec},
                            "stream": PLACEBO_STREAM, "strata": ["scored genes", "other moved genes"],
                            "expression_matched": {"arms": [p for p in (SPERM_ENTRY, *NEW_SPLACEBOS) if p in spec],
                                                   "stream": SPERM_STREAM, "bins": LEVEL_BINS,
                                                   "levels": "mean CP10k on the reference role's cells of the group"}},
                "external": {"entry": ext, **{tag: rec for tag, (_, rec) in eps_sets.items()}}, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), feng_dirs.CONFIG,
                                            n22])})
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(truth: str, rows: Sequence[Mapping], committed_rows: Sequence[Mapping]) -> dict:
    """Per seed, each rescored arm against its committed score; on the chamber truths also its prediction sha256
    against try T46's recorded one (T44's whole-cut file holds none). Holds when every seed and arm is there, each score
    is within REPRO_TOL and, on a chamber truth, each sha256 is equal."""
    arms = REPRO_WHOLE if truth == WHOLE else REPRO_CHAMBER
    sha_arms = () if truth == WHOLE else REPRO_CHAMBER
    by_seed = {r["seed"]: r for r in committed_rows}
    worst, absent, sha_off = {a: 0.0 for a in arms}, [], []
    for r in rows:
        old = by_seed.get(r.get("seed"))
        for arm in arms:
            if old is None or arm not in old["arms"] or arm not in (r.get("arms") or {}):
                absent.append(f"{r.get('seed')}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - old["arms"][arm]["score"]))
        for arm in sha_arms:
            want = ((old or {}).get("sha256") or {}).get(arm)
            if not isinstance(want, str) or (r.get("sha256") or {}).get(arm) != want:
                sha_off.append(f"{r.get('seed')}:{arm}")
    scores_hold = bool(rows) and bool(committed_rows) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "absent": absent, "sha256_off": sha_off, "scores_hold": scores_hold,
            "sha256_holds": not sha_off, "holds": scores_hold and not sha_off}


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
        rows = r.get("rows") or []
        seeds = [row.get("seed") for row in rows]
        if (any(not isinstance(x, int) for x in seeds) or sorted(seeds) != list(SEEDS)
                or list(r.get("seeds", seeds)) != list(SEEDS)):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        if r.get("task") != TASK_OF.get(k):
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("truth") != k:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}")
        absent = sorted({a for row in rows for a in arms_of(k) if a not in (row.get("arms") or {})})
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
        # the pseudocount each direction set was actually made with, as feng_directions recorded the config it was
        # given (not arm_spec's constants, which the pseudocounts record above copies)
        ext = r.get("external") if isinstance(r.get("external"), Mapping) else {}
        bad_ext = [t for t, e in want_ps.items()
                   if _off((ext.get(t) if isinstance(ext.get(t), Mapping) else {}).get("pseudocount_cp10k"), e,
                           EPS_TOL)]
        if bad_ext:
            reasons.append(f"{k}: wrong pseudocount in the direction code's own record (external) for {bad_ext}")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the truths ran at {len(commits)} commits")
    repro = {}
    for k, r in results.items():
        repro[k] = reproduction(k, r.get("rows") or [], committed.get(k) or [])
        if not repro[k]["scores_hold"]:
            src = "T44's whole-cut file" if k == WHOLE else "T46's file"
            reasons.append(f"{k}: {', '.join(REPRO_WHOLE if k == WHOLE else REPRO_CHAMBER)} do not reproduce their "
                           f"committed scores ({src}) within {REPRO_TOL:g} per seed")
        if not repro[k]["sha256_holds"]:
            reasons.append(f"{k}: the rescored arms' predictions do not have T46's recorded sha256 on "
                           f"{repro[k]['sha256_off'][:6]}{' ...' if len(repro[k]['sha256_off']) > 6 else ''}")
    for k, r in results.items():
        arms = arms_of(k)
        if not r.get("rows") or not all(a in (row.get("arms") or {}) for row in r["rows"] for a in arms):
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
        directional = set(arms) - {"b0"}
        fs_bad = any(not isinstance(row.get("factor_stats"), Mapping) or not row["factor_stats"]
                     or any(not isinstance(v, Mapping) or not directional <= set(v)
                            for v in row["factor_stats"].values()) for row in r["rows"])
        pb_bad = any(not isinstance(row.get("pseudobulk_change"), Mapping)
                     or not set(arms) <= set(row["pseudobulk_change"]) for row in r["rows"])
        if fs_bad or pb_bad:
            reasons.append(f"{k}: the count-space factor or pseudobulk records are missing")
        real = set(candidates(k))  # the real directional arms: the read-back covers these
        rb_bad = any(not isinstance(row.get("ruler_readback"), Mapping) or not row["ruler_readback"]
                     or any(not isinstance(v, Mapping) or not real <= set(v)
                            for v in row["ruler_readback"].values()) for row in r["rows"])
        if rb_bad:
            reasons.append(f"{k}: the per-cut ruler read-back records are missing")
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


def same_on_every_seed(r: Mapping, placebo: str) -> bool:
    """Whether every row records that `placebo` permuted positions as its base placebo did (BASE_PLACEBO)."""
    key = SAME_KEY[BASE_PLACEBO[placebo]]
    return bool(r["rows"]) and all(((row.get("placebo") or {}).get(key) or {}).get(placebo) is True
                                   for row in r["rows"])


def gene_identity(r: Mapping, truth: str, tag: str) -> dict:
    """Statement 3 on one chamber truth for one pseudocount: D with the expression-matched placebos (read), and D with
    try T46's placebos (reported, not read: its placebo moves the decompressed factors onto well-expressed genes)."""
    arm = new_arm(tag, PLACEBO_STEP)
    sp = new_splacebo(tag)
    gate = one_truth(r, truth, arm, BASE)
    d = did(r, truth, arm, BASE, sp, SPERM_ENTRY)
    d_t46 = did(r, truth, arm, BASE, new_placebo(tag), PERM_ENTRY)
    same = same_on_every_seed(r, sp)
    if not same:
        label = NOT_READ_PERM
    else:
        label = NO_GAIN if not gate["ahead"] else NEEDS_GENES if d["ahead"] else NOT_SEPARATED
    return {"label": label, "gate": gate, "did": d, "same_positions": same,
            "did_t46_placebo_not_read": {**d_t46, "same_positions": same_on_every_seed(r, new_placebo(tag))}}


def matched_statement(comps: Sequence[Mapping]) -> str:
    """Statement 2 from its two comparisons on the two cuts: both pass, both lose, both neither, or mixed (each
    comparison's label, the smaller delivered size first)."""
    labels = [c["label"] for c in comps]
    if all(x == PASSES for x in labels):
        return SMALLER_BETTER
    if all(x == LOSES for x in labels):
        return UNIFORM_BETTER
    if all(x == NEITHER for x in labels):
        return NO_DIFFERENCE
    return MIXED.format(*labels)


def unmatched(comps: Sequence[Mapping]) -> list[str]:
    """The statement-2 pairs flagged by the per-cut read-back, as 'arm against uniform on cut (groups)'."""
    return [f"{m['arm']} against {m['uniform']} on {c} ({', '.join(m['flagged_groups']) or 'no group read'})"
            for x in comps for c, m in x["ruler_match"].items() if not m["matched"]]


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
    s2 = [{"two_cuts": two_cuts(results, a, u), "libraries": {k: one_truth(results[k], k, a, u) for k in LIBRARIES},
           "ruler_match": {c: ruler_match(results[c], c, a, u) for c in CUTS}} for a, u in MATCHED]
    out["statement_2"] = {"comparisons": s2, "statement": matched_statement([c["two_cuts"] for c in s2]),
                          "not_ruler_matched": unmatched(s2)}
    if out["statement_2"]["not_ruler_matched"]:
        out["statement_2"]["statement"] += UNMATCHED.format("; ".join(out["statement_2"]["not_ruler_matched"]))
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
    out["half_agreement"] = {k: results[k].get("half_agreement") for k in TRUTHS}
    out["same_positions"] = {k: {p: same_on_every_seed(results[k], p) for p in BASE_PLACEBO} for k in CHAMBER_TRUTHS}
    out["count_space"] = {k: count_space(results[k], arms_of(k)) for k in TRUTHS}
    pairing = {a: {k: nearest_uniform(a, k) for k in TRUTHS} for a in NEW_ARMS}
    out["nearest_uniform_differs"] = {a: p for a, p in pairing.items() if len(set(p.values())) > 1}
    out["extrapolated_delivered"] = {a: DELIVERED[a] for a in EXTRAPOLATED}
    return out


def count_space(r: Mapping, arms: Sequence[str]) -> dict:
    """Per directional arm, over seeds: per group the largest |log f| and the most scored genes beyond log 10 on any
    seed, and the mean shares of reference CP10k moved (all genes, scored genes); and the mean pseudobulk change from
    m1_r_gmom. Target-free."""
    rows = r["rows"]
    groups = sorted(rows[0]["factor_stats"])
    out = {}
    for a in arms:
        if a == "b0":
            continue
        rec = {"pseudobulk_change": float(np.mean([row["pseudobulk_change"][a] for row in rows])), "groups": {}}
        for g in groups:
            fs = [row["factor_stats"][g][a] for row in rows]

            def mean_of(q, fs=fs):
                xs = [x.get(q) for x in fs if x.get(q) is not None]
                return float(np.mean(xs)) if xs else None

            rec["groups"][g] = {"max_abs_log_f": float(max(x["max_abs_log_f"] for x in fs)),
                                "scored_over_log10": int(max(x["scored_over_log10"] for x in fs)),
                                **{q: mean_of(q) for q in ("cp10k_moved_share", "cp10k_moved_share_scored",
                                                           "cp10k_moved_share_scored_weak",
                                                           "cp10k_moved_share_scored_strong")}}
        out[a] = rec
    return out


def _sd(x) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:.2f}"


def _line(c: Mapping) -> str:
    sk = ", ".join(f"{m} {x:+.4f}" for m, x in c["skills"].items())
    return (f"{c['change']} against {c['base']}: {c['mean']:+.2f} points, per-seed sd {_sd(c['sd'])}, ahead "
            f"{c['seeds_ahead']} of {c['of']}, per seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}, "
            f"{c['label']}; de terms {c['points']['de']:+.2f}, mmd_u and variogram "
            f"{c['points']['mmd_u_variogram']:+.2f}; skills: {sk}")


def _f3(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def _did_line(d: Mapping) -> str:
    return (f"{d['quantity']}: {d['mean']:+.2f} points, per-seed sd {_sd(d['sd'])}, ahead {d['seeds_ahead']} of "
            f"{d['of']}, per seed {d['range'][0]:+.2f} to {d['range'][1]:+.2f}, test: {d['test']}; de terms "
            f"{d['points']['de']:+.2f}, mmd_u and variogram {d['points']['mmd_u_variogram']:+.2f}")


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
        cs = r["count_space"]

        def size_of(a: str, k: str) -> str:
            c = cs[k][a]
            shares = ", ".join(f"{g} {_f3(x['cp10k_moved_share'])} (scored weak "
                               f"{_f3(x['cp10k_moved_share_scored_weak'])}, strong "
                               f"{_f3(x['cp10k_moved_share_scored_strong'])})" for g, x in c["groups"].items())
            mark = " (extrapolated)" if a in EXTRAPOLATED else ""
            return (f"{a}: ruler size {DELIVERED.get(a, float('nan')):.4g}{mark}, share of reference CP10k moved by "
                    f"group {shares}, pseudobulk change from m1_r_gmom {c['pseudobulk_change']:.4f}")

        def match_line(m: Mapping, k: str) -> str:
            groups = "; ".join(f"{g} {_f3(x['arm'])} against {_f3(x['uniform'])}, half gap {_f3(x['half_gap'])}"
                               + (" FLAGGED" if x["flagged"] else "") for g, x in m["groups"].items())
            verdict = "ruler-matched" if m["matched"] else "NOT ruler-matched"
            return (f"    read-back on {k}'s reference cells (entry's ruler, slope on the entry's direction; "
                    f"{m['arm']} against {m['uniform']}, half the smaller gap to {', '.join(m['neighbours'])}): "
                    f"{verdict}; {groups or 'no group read'}")

        for i, c in enumerate(r["statement_2"]["comparisons"]):
            lines.append(f"Statement 2, comparison {i + 1}:")
            lines += _two_cut_lines(c["two_cuts"])
            lines += [f"    on {k}: " + _line(x) for k, x in c["libraries"].items()]
            for k in CUTS:
                lines.append(f"    count-space sizes on {k} (ruler-matched is not count-matched): "
                             + "; ".join(size_of(a, k) for a in (c["two_cuts"]["change"], c["two_cuts"]["base"])))
                lines.append(match_line(c["ruler_match"][k], k))
        lines.append("Statement 1, each new arm against m1_r_gmom (does not separate decompression from step size):")
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
        lines.append("Statement 3, gene identity of the new directions (chamber truths; the expression-matched "
                     "placebos; bears on gene identity, not on decompression):")
        for k, by in r["statement_3"].items():
            for tag, s in by.items():
                lines.append(f"    on {k}, {tag}: {s['label']} ({S3_CAVEAT}); gate: " + _line(s["gate"]))
                lines.append(f"        on {k}, {tag}: " + _did_line(s["did"]))
                lines.append(f"        on {k}, {tag}: reported, not read (try T46's placebo, which moves the "
                             "decompressed factors onto well-expressed genes): "
                             + _did_line(s["did_t46_placebo_not_read"]))
        if r["nearest_uniform_differs"]:
            lines.append("The nearest-delivered-size uniform arm differs between the truths for: " + "; ".join(
                f"{a} (" + ", ".join(f"{k} {u}" for k, u in p.items()) + ")"
                for a, p in r["nearest_uniform_differs"].items()))
        lines.append("Extrapolated delivered sizes (not quoted read-backs): " + ", ".join(
            f"{a} {x:.4g}" for a, x in r["extrapolated_delivered"].items()))
        for k in TRUTHS:
            lines.append(f"truth {k}: highest mean {r['highest_mean'][k]}; mean scores "
                         + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
            for a, x in r["new_arm_differences"][k].items():
                lines.append(f"    {k}: " + _line(x["against_m1_r_gmom"]))
                mark = " (the arm's size is extrapolated)" if a in EXTRAPOLATED else ""
                lines.append(f"    {k}: nearest delivered size {x['nearest_uniform']}{mark}: "
                             + _line(x["against_nearest_uniform"]))
            sp = r["sparsity"][k]
            lines.append("    scored genes, mean over seeds: " + "; ".join(
                f"{a} zero share {x['zero_share']:.4f}, detected per cell {x['detected_per_cell']:.2f}"
                for a, x in sp.items()))
            lines.append("    mean |log factor| on the scored genes by group: " + "; ".join(
                f"{g}: " + ", ".join(f"{a} {x:.4f}" for a, x in fs.items())
                for g, fs in r["scored_log_factors"][k].items()))
            lines.append("    count space by arm (largest |log f| on any gene and most scored genes beyond log 10 on "
                         "any seed, no clip on the route; mean shares of reference CP10k moved, all and scored genes; "
                         "pseudobulk change from m1_r_gmom): " + "; ".join(
                             f"{a} pseudobulk {c['pseudobulk_change']:.4f} " + ", ".join(
                                 f"{g} max {x['max_abs_log_f']:.3f} over-log10 {x['scored_over_log10']} moved "
                                 f"{_f3(x['cp10k_moved_share'])}/{_f3(x['cp10k_moved_share_scored'])}"
                                 for g, x in c["groups"].items()) for a, c in r["count_space"][k].items()))
            ds = r["direction_stats"][k] or {}
            lines.append("    directions (squared size over the entry's; shares on genes with a Feng mean under "
                         f"{LOW_CP10K:g} CP10k at E9.5, at E14.5, at either end): " + "; ".join(
                             f"{tag} {g} " + ", ".join("n/a" if x is None else f"{x:.3f}" for x in rec.values())
                             for tag, by in ds.items() for g, rec in by.items()))
            ha = r["half_agreement"][k] or {}
            lines.append("    Feng half-A against half-B directions (correlation over the moved genes; eb_shrink "
                         "factor per decile, lowest expression first): " + "; ".join(
                             f"{tag}: {by['error']}" if "error" in by else f"{tag} " + ", ".join(
                                 f"{g} {_f3(x['correlation'])} [" + " ".join(f"{f:.2f}" for f in x["eb_shrink_factors"])
                                 + "]" for g, x in by.items()) for tag, by in ha.items()))
            if k in r["same_positions"]:
                lines.append("    permutation of positions the same as the base placebo's on every seed: " + ", ".join(
                    f"{p} (against {BASE_PLACEBO[p]}) {x}" for p, x in r["same_positions"][k].items()))
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
