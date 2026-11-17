#!/usr/bin/env python
"""Try T50 (experiments/harness/TRIES.md): how far M1's step can over-deliver before the score turns, and whether the
flat point follows the truth's depth. A labelled read on try T46's route (experiments/t1_m1/lib_step_read.py).

Tries T44 to T47 found larger steps better on every released truth of the E10.5 chamber cut, the gain needing Feng's
own genes and not count noise, and library 3's curve (the dissected heart, the deepest library) flattening between 3
and 4 times the step while libraries 1 and 2 (embryo portions) still rise. This read extends the curve to 5, 6 and 8
times the step, and puts the same steps to library 3's truth thinned to about the depth of libraries 1 and 2.

  run    One truth of the chamber cut T1/E9.5->E10.5:chamber, 10 seeds, on T46's route. `pooled`, `lib1`, `lib2`,
         `lib3`: T46's truths, built by T46's own functions (lib_step_read.truth_and_contamination). `lib3_t33`,
         `lib3_t67`, `lib3_t33b`: lib3's rows exactly as T46 makes them (the same rows in the same order), with every
         integer count c of those cells, on every gene of the released file, replaced by a draw from Binomial(c, p):
         p = 1/3, 2/3 and 1/3, one generator per name on the stream [50, 1] (lib3_t33), [50, 2] (lib3_t67) and
         [50, 3] (lib3_t33b, added after review: a second draw at p = 1/3, everything else as lib3_t33's, statement
         3's self-check), cells in the truth's order and within a cell the file's gene order (CSR order). The stream
         does not depend on the build seed: one thinned truth per name, the same for all ten seeds. How the matrix
         reaches the scorer:
         harness.filter_target hands the scorer Stage.from_anndata of the released file's rows, i.e. the file's .X
         (log1p CP10k of the cell's counts over every gene of the file) on the panel's columns; T1's scRNA file holds
         no counts layer. The integer counts are recovered by the repository's own routine, vec.core.counts_of
         (vec.prep.recover_counts, the one vec.prep.load_stage uses, refused above a residual of 1e-2), from the full
         file's rows; the thinned counts are turned back by vec.prep.log1p_cp10k (the released convention, which the
         core and the T1 head use for every count edit) and Stage.from_anndata on the panel. Before thinning, the run
         checks that the recovered counts, through the same two routines, give back lib3's matrix as the scorer reads
         it (largest absolute difference at most REBUILD_TOL, recorded); otherwise it refuses: the thinning would not
         be of the matrix the scorer reads. The thinned truth keeps lib3's rows, so its board, its limits, its cell
         count and every prediction are lib3's; the floor's anchors and the scoring are made for the thinned truth as
         the harness makes them for any target. Arms per seed, all from one cell draw: b0 and m1_r_gmom as T42 builds
         them, and step_read.step_arms at multipliers 1.5, 2, 3, 4 (T46's m1_r_gmom_x15, _x2, _x3, _x4, rescored for
         the gate) and 5, 6, 8 (m1_r_gmom_x5, _x6, _x8: fractions 1.0, 1.2 and 1.6), each the same call with the
         fraction multiplied and every other argument the same object. No placebo arm. The route takes a fraction
         above 1 with no clip, guard or refusal: external_shift only refuses a non-finite fraction, M1's build passes
         no max_abs_log_factor (so no per-gene factor is clipped), and apply_factors refuses only non-finite or
         non-positive factors (exp of a finite log factor is neither). The task file holds every arm's fraction, the
         truth's name, the truth's row count and contamination count (checks, compared by code and never printed),
         per seed and arm the sha256 of the predicted matrix as the scorer reads it and the scored genes' zero share
         and detected genes per cell, per seed the floor's and ceiling's metric values (this pseudo-task's computed
         anchors), and, for a thinned truth, p, the stream, the realised depth ratio (the sum of
         the truth's counts after thinning over the sum before: a ratio), the rebuild check's largest difference and
         whether any cell was left with no count. Refused on an existing task file, unless the committed N22 passed,
         and unless experiments/harness/TRIES.md holds row T50.
  read   T50's reading of the six truths' files (rule below, fixed in the row before the run).

The reading. Void unless the seven truths hold the ten seeds at one commit on a clean tree (tree_clean exactly true); no
file records an anchor failure; every arm is on every row, and each file's task is the chamber cut, its truth its own
name, its base fraction 0.2 and every arm's fraction its multiplier times 0.2; on pooled, lib1, lib2 and lib3 the six
rescored arms (b0, m1_r_gmom, _x15, _x2, _x3, _x4) equal, seed for seed, their committed scores in
experiments/t1_m1/results/t46_m1_lib_step/tasks/, each to 1e-9; on lib3_t33, lib3_t67 and lib3_t33b every arm's
recorded prediction sha256 equals, seed for seed, the one in this run's lib3 file; each thinned truth records its own p
and stream, its realised depth ratio is within 0.01 of its p, it records that no truth cell was left with no count
(any_cell_emptied exactly false; added after review) and its rebuild check's largest difference is at most REBUILD_TOL
(added after review); on every truth no new arm (_x5, _x6, _x8) scores exactly as m1_r_gmom_x4 on every seed. Compared
by code and never printed: the three library truths' row counts and contamination counts add up to the pooled truth's,
and each thinned truth's equal lib3's. A file with no rows is a void reason, never an error.
The test on one truth is T46's single-task test (vec.adoption.keep_change on that truth alone): "ahead" is a mean
paired gain above 0.3 with at least 8 of 10 seeds ahead; "behind" the same with the sign reversed; otherwise "level".
Statement 1, over-delivery, per truth (all six): each of x5, x6 and x8 against x4: ahead, behind or level; and the
multiplier with the highest mean board score among 1, 1.5, 2, 3, 4, 5, 6 and 8 (the smallest on a tie; "still rising
at 8" if that is 8).
Statement 2, the worst cases a request would quote, on lib3 (and reported the same way on the other five truths): x4
against x3, x5 against x3, x6 against x3, x6 against x4, x8 against x4 and x8 against x5, each with its mean paired
difference, seeds ahead, per-seed range and its label by the single-task test.
Statement 3, depth. Per seed D33 = (x4 - x2 on lib3_t33) - (x4 - x2 on lib3), and D67 the same with lib3_t67 (the same
predictions scored against two truths, paired by seed). "The larger steps gain more against the shallower truth" if D33
is ahead by the single-task test; "The larger steps gain less against the shallower truth" if it is behind; otherwise
"no difference at the bar". D67 gets the same label and is reported, as are (x8 - x4) on each thinned truth minus the
same on lib3, and, beside them and read by nothing, lib1's and lib2's (x4 - x2) minus lib3's (different truths, not
paired in kind).
Statement 3's self-check (added after review; each thinned truth is one binomial draw shared by all ten seeds, so 8 of
10 seeds counts agreement within one draw, not across draws): the draw null, per seed (x4 - x2 on lib3_t33b) - (x4 -
x2 on lib3_t33), the same predictions scored against two draws at p = 1/3. If it is not "level" by the single-task
test, statement 3 is void and only reported (D33's and D67's statements read DEPTH_VOID, their labels are kept as
label_if_valid); otherwise statement 3 stands, with the null's mean and per-seed range printed beside D33. The void
of statement 3 does not void statements 1 and 2. (x8 - x4) between the two draws is reported beside it. A "more" or
"less" reading is a property of scoring E9.5-depth predictions against a shallower truth of the same composition (in
log space); it shows no biological change and no change in Feng's share (see the limits).
Reported beside statement 3 and read by nothing: D33 and D67 on each metric's raw values (free of each truth's
anchors), each of lib3's and the thinned truths' mean ceiling minus floor per metric, and the curve's shape on lib3 and
the thinned truths side by side (statement 1's labels and argmax), the comparison that does not depend on scale.
The outcome line gives per truth the three labels against x4 first, then the argmax, marked "(not separated)" when no
label is ahead or behind (added after review).
Reported: per truth every arm's mean board score, each neighbour gain per 0.1 of fraction, every metric's mean skill and
raw difference for the comparisons of statement 2, the part of each on the two de terms and on mmd_u and the variogram
(as T46 computes them), each thinned truth's realised depth ratio, rebuild check's largest difference and whether any
cell was left with no count, and whether the predictions' sha256 agree across the truths. The read puts nothing forward
by itself.

NOT an adoption result: Feng's E9.5 cells lie inside this pseudo-task's external window, the cut is one earlier T1
choices were made on, and every truth is a part of a released target read in-sample. The result can only go to the
team captain.

Limits, repeated in the text and in the row:
- One released day, read in-sample on the chamber cut only.
- Multiples of 5 and more apply Feng's whole five-day change or more to a one-day truth; they are read to find where
  the score turns, not as candidates.
- Thinning lowers depth and raises count noise together and cannot tell them apart.
- A thinned dissected-heart truth is not an embryo-portion truth, and libraries 1 and 2 also differ from it in
  contamination share and capture, so agreement with libraries 1 and 2 is shown beside it and proves nothing about
  kind.
- Nothing here says how deep any hidden truth is, and no statement may be used with a figure for one.
- No multiplier read from T50 is a candidate setting; the test's step stays the bracket fraction by rule.
- No outcome of T50 is evidence on T1 read depth, a lever closed after an exposure that this read does not reopen.
- Board points are not comparable across truths: each truth, the thinned ones included, has its own floor and ceiling;
  D33 and D67 mix a rescaled anchor with any change in the curve's shape.
- Where adjacent multipliers differ by less than their per-seed spread the argmax is not separated; the argmax alone
  is not a finding.
- 'Level' in statement 2 is not 'no loss'.
- Each thinned truth is one draw shared by the ten seeds; statement 3 stands only when the draw null reads level.
- A 'more' or 'less' reading of statement 3 is a property of depth mismatch in log space at fixed composition.
- The depth-ratio check catches a wrong p or a wrong stream only.
(LIMITS below holds the full sentences, printed in the text.)

Usage: on the DGX campaign checkout, repo root, all seven truths at one commit, committing nothing until the seven
task files exist, and on the machine T46 ran on, so that the rescored arms reproduce. The environment variables are
docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py run --truth pooled --out t50_m1_over_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py run --truth lib1 --out t50_m1_over_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py run --truth lib2 --out t50_m1_over_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py run --truth lib3 --out t50_m1_over_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py run --truth lib3_t33 --out t50_m1_over_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py run --truth lib3_t67 --out t50_m1_over_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py run --truth lib3_t33b --out t50_m1_over_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/over_read.py read --out t50_m1_over_step
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from vec import core, harness, prep, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS, keep_change
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import lib_step_read as ls  # noqa: E402

er, gc, sr, feng_dirs = ls.er, ls.gc, ls.sr, ls.feng_dirs
gv, hygiene, n12 = ls.gv, ls.hygiene, ls.n12
TASK = ls.TASK
T46_TASKS = RESULTS / "t46_m1_lib_step" / "tasks"
POOLED, LIBRARIES, HEART = ls.POOLED, ls.LIBRARIES, ls.HEART
THINNED = {"lib3_t33": {"p": 1.0 / 3.0, "stream": [50, 1]}, "lib3_t67": {"p": 2.0 / 3.0, "stream": [50, 2]},
           "lib3_t33b": {"p": 1.0 / 3.0, "stream": [50, 3]}}
DRAW_NULL = ("lib3_t33b", "lib3_t33")  # a second draw at p = 1/3 against the first: statement 3's self-check
TRUTHS = (*ls.TRUTHS, *THINNED)
BASE = sr.BASE  # m1_r_gmom, fraction 0.2
OLD = dict(ls.MULTIPLIERS)  # T46's four scaled arms
NEW = {"m1_r_gmom_x5": 5.0, "m1_r_gmom_x6": 6.0, "m1_r_gmom_x8": 8.0}
MULTIPLIERS = {**OLD, **NEW}
X2, X3, X4, X5 = "m1_r_gmom_x2", "m1_r_gmom_x3", "m1_r_gmom_x4", "m1_r_gmom_x5"
X6, X8 = "m1_r_gmom_x6", "m1_r_gmom_x8"
ARMS = ("b0", BASE, *MULTIPLIERS)
RESCORED = ("b0", BASE, *OLD)  # against T46's committed files, on T46's four truths
ROW = "T50"
FRACTION = sr.FRACTION  # 0.2
FRACTION_TOL = sr.FRACTION_TOL
REPRO_TOL = sr.REPRO_TOL  # 1e-9
DEPTH_TOL = 0.01  # the realised depth ratio against p
REBUILD_TOL = 1e-5  # recovered counts through log1p_cp10k against lib3's matrix as the scorer reads it (float32)
WORST_CASES = ((X4, X3), (X5, X3), (X6, X3), (X6, X4), (X8, X4), (X8, X5))  # (change, base), statement 2
STILL_RISING = "still rising at 8"
MORE = "The larger steps gain more against the shallower truth"
LESS = "The larger steps gain less against the shallower truth"
NO_DIFFERENCE = "no difference at the bar"
DEPTH_VOID = ("void: the second thinning draw (lib3_t33b) against the first (lib3_t33) is not level by the single-task "
              "test, so one draw alone can pass the bar; statement 3 is only reported")
NOT_SEPARATED = "(not separated)"
TASK_FILE = "{truth}.json"
LABEL_T50 = ("try T50, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng step multiplied by 1, 1.5, 2, 3, 4, 5, 6 and 8, on T46's route, with the "
             "chamber cut's truth taken pooled, library by library and as library 3 thinned to a third (two draws) and "
             "two thirds of its depth; Feng's E9.5 cells lie inside this pseudo-task's external window and the cut "
             "is one earlier T1 choices were made on; the result can only go to the team captain")
LIMITS = (
    "One released day, read in-sample on the chamber cut only.",
    "Multiples of 5 and more apply Feng's whole five-day change or more to a one-day truth; they are read to find "
    "where the score turns, not as candidates.",
    "Thinning lowers depth and raises count noise together and cannot tell them apart.",
    "A thinned dissected-heart truth is not an embryo-portion truth, and libraries 1 and 2 also differ from it in "
    "contamination share and capture, so agreement with libraries 1 and 2 is shown beside it and proves nothing about "
    "kind.",
    "Nothing here says how deep any hidden truth is, and no statement may be used with a figure for one.",
    "No multiplier read from T50 is a candidate setting; the test's step stays the bracket fraction by rule.",
    "No outcome of T50 is evidence on T1 read depth, a lever closed after an exposure that this read does not reopen.",
    "Board points are not comparable across truths: each truth, the thinned ones included, has its own floor and "
    "ceiling; D33 and D67 mix a rescaled anchor with any change in the curve's shape, so the raw-metric differences "
    "and the thinned truths' labels beside lib3's are reported with them.",
    "Where adjacent multipliers differ by less than their per-seed spread the argmax is not separated; the argmax "
    "alone is not a finding.",
    "A label of 'level' in statement 2 means the difference is not separated at a 0.3 mean with 8 of 10 seeds on one "
    "in-sample day; 'level' is not 'no loss'.",
    "Each thinned truth is one binomial draw shared by all ten seeds, so its seeds agree within one draw, not across "
    "draws; statement 3 stands only when a second draw at p = 1/3 (lib3_t33b) reads level against the first.",
    "A 'more' or 'less' reading of statement 3 is a property of scoring E9.5-depth predictions against a shallower "
    "truth of the same composition (in log space); it shows no biological change and no change in Feng's share, and "
    "says nothing about which turn applies at the test without a figure for the hidden truth's depth, which this read "
    "may not use.",
    "The depth-ratio check catches a wrong p or a wrong stream only; it says nothing about whether the thinned depth "
    "matches libraries 1 and 2.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The thinned truths


def thin_counts(C, name: str) -> sparse.csr_matrix:
    """Every integer count c replaced by a draw from Binomial(c, p), p and the generator the name's (THINNED), cells in
    row order and within a cell the columns' order (CSR order). No dependence on any build seed."""
    spec = THINNED[name]
    C = sparse.csr_matrix(C, dtype=np.int64, copy=True)
    C.sort_indices()
    rng = np.random.default_rng(list(spec["stream"]))
    out = sparse.csr_matrix((rng.binomial(C.data, spec["p"]), C.indices.copy(), C.indptr.copy()), shape=C.shape)
    out.eliminate_zeros()
    return out


def as_scored(counts, obs: pd.DataFrame, var_names: Sequence[str], genes) -> Stage:
    """Counts as the scorer reads a target: log1p CP10k (vec.prep's) over every gene given, then Stage.from_anndata
    on the panel (the route harness.filter_target hands the scorer)."""
    a = ad.AnnData(X=prep.log1p_cp10k(counts), obs=obs.copy(), var=pd.DataFrame(index=list(var_names)))
    return Stage.from_anndata(a, genes)


def lib3_counts(t: harness.PseudoTask, data: Path, genes) -> tuple[Stage, int, sparse.csr_matrix, pd.DataFrame, list]:
    """lib3 exactly as T46 makes it (its Stage and contamination count), and the integer counts of those rows, in that
    order, on every gene of the released file (vec.core.counts_of), with their obs and the file's gene names."""
    target, names, contaminated = ls.filtered_target(t, data, genes)
    rows = ls.library_rows(names, HEART)
    truth = target.take(rows)
    full = read_stage(data / t.target, None)
    if not full.obs_names.is_unique:
        raise ValueError(f"{t.target}: cell names are not unique; lib3's rows cannot be found in the full file")
    pos = pd.Index(full.obs_names).get_indexer(np.asarray(names)[rows])
    if (pos < 0).any():
        raise ValueError(f"{t.target}: lib3's cells are not all in the full file")
    sub = full[pos]
    C = core.counts_of(sub)
    return truth, int(contaminated[rows].sum()), C, sub.obs[[c for c in ("celltype",) if c in sub.obs]].copy(), \
        list(full.var_names)


def thinned_truth(t: harness.PseudoTask, data: Path, genes, name: str) -> tuple[Stage, int, dict]:
    """lib3's rows with their counts thinned (thin_counts), as the scorer reads them; lib3's contamination count; the
    record (p, stream, the realised depth ratio, the rebuild check's largest difference, whether a cell was left
    empty). Refused when the recovered counts do not give back lib3's matrix within REBUILD_TOL."""
    lib3, contamination, C, obs, var_names = lib3_counts(t, data, genes)
    rebuilt = as_scored(C, obs, var_names, genes)
    worst = float(np.abs(rebuilt.X.astype(np.float64) - lib3.X.astype(np.float64)).max()) if lib3.X.size else 0.0
    if rebuilt.X.shape != lib3.X.shape or not worst <= REBUILD_TOL:
        raise ValueError(f"{name}: the counts recovered from the released file do not give back lib3's matrix as the "
                         f"scorer reads it (largest difference {worst:.3g} > {REBUILD_TOL:g}); T50 does not thin "
                         "normalised values")
    if lib3.celltype is not None and not np.array_equal(rebuilt.celltype, lib3.celltype):
        raise ValueError(f"{name}: the rebuilt rows' cell types differ from lib3's")
    thinned = thin_counts(C, name)
    before = float(C.sum())
    rec = {"p": THINNED[name]["p"], "stream": list(THINNED[name]["stream"]), "of": HEART,
           "depth_ratio": float(thinned.sum()) / before if before > 0 else float("nan"),
           "rebuild_largest_difference": worst,
           "any_cell_emptied": bool(((np.asarray(thinned.sum(1)).ravel() == 0)
                                     & (np.asarray(C.sum(1)).ravel() > 0)).any())}
    return as_scored(thinned, obs, var_names, genes), contamination, rec


def truth_and_contamination(t: harness.PseudoTask, data: Path, genes, truth: str) -> tuple[Stage, int, dict | None]:
    """T46's truths by T46's own function; a thinned truth by thinned_truth."""
    if truth in THINNED:
        return thinned_truth(t, data, genes, truth)
    if truth not in ls.TRUTHS:
        raise ValueError(f"unknown truth {truth!r}; one of {list(TRUTHS)}")
    return (*ls.truth_and_contamination(t, data, genes, truth), None)


# ---------------------------------------------------------------------------------------------------------------------
# The run


def arm_fractions(fraction: float) -> dict[str, float]:
    return sr.arm_fractions(fraction, MULTIPLIERS)


def scaled_arms(ctx: harness.BuildContext, directions: Mapping, groups: Mapping, fraction: float,
                held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """m1_r_gmom and the seven scaled arms: step_read.step_arms at this read's multipliers."""
    return sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors, multipliers=MULTIPLIERS)


def build_seed(ctx: harness.BuildContext, b0, directions: Mapping, groups: Mapping, fraction: float,
               held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """One seed's arms in the order the run builds them: b0, m1_r_gmom, the scaled arms."""
    return {"b0": b0(ctx), **scaled_arms(ctx, directions, groups, fraction, held_m1, held_factors)}


def anchor_record(anchors: Mapping) -> dict:
    """The floor's and ceiling's metric values per ranked metric (this pseudo-task's computed anchors, from released
    data), so that a reader can tell a rescaled anchor from a curve that moved."""
    return {m: {"floor": float(v["floor"]), "ceiling": float(v["ceiling"])} for m, v in anchors.items()}


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T50 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T50 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = TASK_FILE.format(truth=a.truth)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    try:
        truth, contamination, thinning = truth_and_contamination(t, data, genes, a.truth)
    except ValueError as e:
        sys.exit(f"T50 is refused on {a.truth}: {e}")
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T50, "task": t.key, "truth": a.truth, "truth_rows": int(truth.n),
           "truth_contamination": contamination, "seeds": list(seeds)}
    if thinning is not None:
        res["thinning"] = thinning
    floors, failures = ls.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        (out / "tasks").mkdir(parents=True, exist_ok=True)
        (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")
        sys.exit(f"T50 is refused on {a.truth}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T50 is void)")
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    directions, ext = feng_dirs.feng_directions(list(pool.var_names))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = build_seed(ctx, b0, directions, groups, fraction, held_m1, preregistered)
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "anchors": anchor_record(f.anchors), "arms": arms,
                     "sha256": {k: ls.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: ls.sparsity(p) for k, p in preds.items()}, "seconds": round(time.time() - t0, 1)})
        print(f"{a.truth} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res.update({"fraction": fraction, "fractions": arm_fractions(fraction), "multipliers": dict(MULTIPLIERS),
                "external": ext, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(results: Mapping[str, Mapping], committed: Mapping[str, Mapping]) -> dict:
    """On T46's four truths, per seed, each rescored arm against T46's committed file for that truth; holds when every
    truth, seed and arm is there and each is within REPRO_TOL."""
    worst = {k: {a: 0.0 for a in RESCORED} for k in ls.TRUTHS}
    absent = []
    for k in ls.TRUTHS:
        if k not in results:
            continue
        old = {r["seed"]: r for r in (committed.get(k) or {}).get("rows", [])}
        for r in _rows(results[k]):
            for arm in RESCORED:
                o = old.get(r["seed"])
                if o is None or arm not in o["arms"] or arm not in r["arms"]:
                    absent.append(f"{k}:{r['seed']}:{arm}")
                    continue
                worst[k][arm] = max(worst[k][arm], abs(r["arms"][arm]["score"] - o["arms"][arm]["score"]))
    holds = (all(k in results and _rows(results[k]) for k in ls.TRUTHS) and not absent
             and all(w <= REPRO_TOL for v in worst.values() for w in v.values()))
    return {"largest_difference": worst, "absent": absent, "holds": holds}


def _rows(r: Mapping) -> list:
    """A task file's rows; a file with no rows key holds none (a void reason, never a KeyError)."""
    return list(r.get("rows") or [])


def _count(x) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


def _sizes_agree(results: Mapping[str, Mapping], key: str) -> list[str]:
    """The library truths' sizes add up to the pooled truth's, and each thinned truth's equals lib3's (compared by
    code, never printed)."""
    bad = []
    if all(k in results for k in ls.TRUTHS):
        sizes = [results[k].get(key) for k in ls.TRUTHS]
        if not all(_count(x) for x in sizes) or sum(sizes[1:]) != sizes[0]:
            bad.append("the library truths do not add up to the pooled truth")
    if HEART in results:
        for k in THINNED:
            if k in results and not (_count(results[k].get(key)) and results[k].get(key) == results[HEART].get(key)):
                bad.append(f"{k} does not equal lib3")
    return bad


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Mapping]) -> dict:
    """The void conditions. results: truth -> task file (a missing truth is absent); committed: T46's truth -> its
    committed task file."""
    reasons = []
    missing = [k for k in TRUTHS if k not in results]
    if missing:
        reasons.append(f"missing truths {missing}")
    want_fr = arm_fractions(FRACTION)
    for k, r in results.items():
        if r.get("anchor_failures"):
            reasons.append(f"{k}: the floor's anchors invert or go missing on seeds "
                           f"{[f['seed'] for f in r['anchor_failures']]}")
        seeds = [row["seed"] for row in _rows(r)]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in _rows(r) for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if r.get("task") != TASK:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("truth") != k:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items()
               if not isinstance(got_fr.get(a), int | float) or abs(got_fr[a] - f) > FRACTION_TOL]
        if bad:
            reasons.append(f"{k}: arm fractions {got_fr}, not {want_fr} (wrong fraction on {bad})")
        if any(not isinstance((row.get("sha256") or {}).get(a), str) for row in _rows(r) for a in ARMS):
            reasons.append(f"{k}: sha256 not recorded for every arm")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the truths ran at {len(commits)} commits")
    repro = reproduction(results, committed)
    if not repro["holds"]:
        reasons.append(f"{', '.join(RESCORED)} do not reproduce T46's committed scores on {', '.join(ls.TRUTHS)} "
                       f"within {REPRO_TOL:g} per seed")
    for key, what in (("truth_rows", "cells"), ("truth_contamination", "contamination cells")):
        for b in _sizes_agree(results, key):
            reasons.append(f"{what}: {b}")
    for k, spec in THINNED.items():
        if k not in results:
            continue
        th = results[k].get("thinning") or {}
        if th.get("p") != spec["p"] or list(th.get("stream") or []) != list(spec["stream"]):
            reasons.append(f"{k}: the file records p {th.get('p')!r} and stream {th.get('stream')!r}, not "
                           f"{spec['p']!r} and {spec['stream']}")
        ratio = th.get("depth_ratio")
        if not isinstance(ratio, int | float) or not abs(ratio - spec["p"]) <= DEPTH_TOL:
            reasons.append(f"{k}: realised depth ratio {ratio!r} is not within {DEPTH_TOL:g} of p {spec['p']:.4f}")
        if th.get("any_cell_emptied") is not False:
            reasons.append(f"{k}: the thinning left a truth cell with no count (or did not record whether it did): "
                           f"any_cell_emptied {th.get('any_cell_emptied')!r}; such a cell reaches the scorer as an "
                           "all-zero row")
        rebuild = th.get("rebuild_largest_difference")
        if not isinstance(rebuild, int | float) or not rebuild <= REBUILD_TOL:
            reasons.append(f"{k}: the rebuild check's largest difference {rebuild!r} is not at most {REBUILD_TOL:g}")
        lib3 = {row["seed"]: row for row in _rows(results.get(HEART, {}))}
        differ = sorted({f"{row['seed']}:{a}" for row in _rows(results[k]) for a in ARMS
                         if (row.get("sha256") or {}).get(a) is None
                         or (lib3.get(row["seed"], {}).get("sha256") or {}).get(a) != row["sha256"][a]})
        if HEART not in results or differ:
            reasons.append(f"{k}: the predictions' sha256 do not equal lib3's seed for seed (seed:arm {differ})")
    for k, r in results.items():
        if not _rows(r) or not all(a in row["arms"] for row in _rows(r) for a in ARMS):
            continue
        same = [a for a in NEW if all(row["arms"][a]["score"] == row["arms"][X4]["score"] for row in _rows(r))]
        if same:
            reasons.append(f"{k}: new arms score exactly as {X4} on every seed (built at x4's fraction?): {same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def level(c: Mapping) -> str:
    """The single-task test's label: ahead, behind or level."""
    return "ahead" if c["ahead"] else "behind" if c["behind"] else "level"


def comparison(r: Mapping, truth: str, change: str, base: str) -> dict:
    """lib_step_read.truth_comparison (step_read.compare on this truth alone, its sd, the reversed test, the points
    by part) with the label ahead / behind / level."""
    c = ls.truth_comparison(r, truth, change, base=base)
    return {**c, "label": level(c)}


def turn(mean_scores: Mapping[str, float]) -> dict:
    """The multiplier with the highest mean board score among 1 (m1_r_gmom) and the seven scaled arms (the smallest on
    a tie); STILL_RISING when that is 8."""
    by_m = {1.0: mean_scores[BASE], **{m: mean_scores[a] for a, m in MULTIPLIERS.items()}}
    best = max(by_m, key=lambda m: (by_m[m], -m))
    return {"by_multiplier": {f"{m:g}": v for m, v in by_m.items()}, "best": best,
            "statement": STILL_RISING if best == max(MULTIPLIERS.values()) else f"turns at {best:g}",
            "tie_rule": "smallest multiplier"}


def neighbour_of(arm: str) -> str:
    """The arm at the next smaller multiplier (m1_r_gmom for x1.5)."""
    ms = sorted(MULTIPLIERS.items(), key=lambda x: x[1])
    i = [a for a, _ in ms].index(arm)
    return BASE if i == 0 else ms[i - 1][0]


def _multiplier(arm: str) -> float:
    return 1.0 if arm == BASE else MULTIPLIERS[arm]


def neighbours(r: Mapping, truth: str) -> dict:
    """Each scaled arm against its next smaller neighbour, with the gain per 0.1 of fraction."""
    out = {}
    for a in MULTIPLIERS:
        b = neighbour_of(a)
        c = comparison(r, truth, a, b)
        tenths = (MULTIPLIERS[a] - _multiplier(b)) * FRACTION / 0.1
        out[a] = {"base": b, "mean": c["mean"], "label": c["label"], "per_tenth": c["mean"] / tenths}
    return out


def _increment(row: Mapping, hi: str, lo: str) -> float:
    return row["arms"][hi]["score"] - row["arms"][lo]["score"]


def paired_difference(other: Mapping, ref: Mapping, hi: str, lo: str) -> dict:
    """Per seed, (hi - lo on `other`) - (hi - lo on `ref`), with the single-task test (keep_change on the pair (ref's
    increment, other's increment)) and the reversed test; mean, seeds ahead, per-seed range and sd."""
    a = {row["seed"]: _increment(row, hi, lo) for row in _rows(ref)}
    b = {row["seed"]: _increment(row, hi, lo) for row in _rows(other)}
    seeds = sorted(set(a) & set(b))
    paired = {s: (a[s], b[s]) for s in seeds}
    fwd = keep_change({"d": paired})["tasks"]["d"]
    rev = keep_change({"d": {s: (y, x) for s, (x, y) in paired.items()}})["tasks"]["d"]
    d = np.array([b[s] - a[s] for s in seeds])
    sd = float(d.std(ddof=1)) if len(d) > 1 else float("nan")
    ahead, behind = bool(fwd.get("passes")), bool(rev.get("passes"))
    return {"quantity": f"({hi} - {lo}) on {other['truth']} - ({hi} - {lo}) on {ref['truth']}",
            "mean": float(d.mean()) if len(d) else float("nan"), "seeds_ahead": int((d > 0).sum()), "of": len(d),
            "range": [float(d.min()), float(d.max())] if len(d) else [float("nan")] * 2, "sd": sd,
            "ahead": ahead, "behind": behind, "label": "ahead" if ahead else "behind" if behind else "level",
            "raw_by_metric": raw_paired_difference(other, ref, hi, lo)}


def raw_paired_difference(other: Mapping, ref: Mapping, hi: str, lo: str) -> dict:
    """The same quantity per metric on the raw metric values (each arm's recorded 'metrics'), mean over seeds: free of
    each truth's anchors, so a reader can tell a rescaled anchor from a curve that moved. Read by nothing."""
    a = {row["seed"]: row for row in _rows(ref)}
    b = {row["seed"]: row for row in _rows(other)}
    seeds = sorted(set(a) & set(b))
    out = {}
    for m in ls.WEIGHTS:
        def inc(row, m=m):
            return row["arms"][hi]["metrics"].get(m, np.nan) - row["arms"][lo]["metrics"].get(m, np.nan)
        d = np.array([inc(b[s]) - inc(a[s]) for s in seeds], dtype=float)
        if len(d) and np.isfinite(d).all():
            out[m] = {"mean": float(d.mean()), "range": [float(d.min()), float(d.max())]}
    return out


def anchor_gap(r: Mapping) -> dict | None:
    """Per metric, the mean over seeds of the ceiling's minus the floor's value (each seed's computed anchors), when the
    task file records them; None otherwise. Read by nothing."""
    rows = _rows(r)
    if not rows or not all(isinstance(row.get("anchors"), Mapping) for row in rows):
        return None
    out = {}
    for m in ls.WEIGHTS:
        g = [row["anchors"][m]["ceiling"] - row["anchors"][m]["floor"] for row in rows if m in row["anchors"]]
        if len(g) == len(rows):
            out[m] = float(np.mean(g))
    return out


def depth_statement(d33: Mapping) -> str:
    return MORE if d33["ahead"] else LESS if d33["behind"] else NO_DIFFERENCE


def depth_null(results: Mapping[str, Mapping]) -> dict:
    """Statement 3's self-check: per seed, (x4 - x2 on lib3_t33b) - (x4 - x2 on lib3_t33), two thinning draws at
    p = 1/3 scored with the same predictions; it must read level by the single-task test, or statement 3 is void.
    (x8 - x4) between the two draws is reported beside it."""
    b, a = (results[k] for k in DRAW_NULL)
    d = paired_difference(b, a, X4, X2)
    return {"x4_minus_x2": d, "x8_minus_x4": paired_difference(b, a, X8, X4), "holds": d["label"] == "level"}


def outcome_part(s1: Mapping) -> str:
    """One truth's statement 1 for the outcome line: the three labels against x4 first, then the argmax, marked not
    separated when no label is ahead or behind."""
    labels = ", ".join(f"x{a.rsplit('_x', 1)[1]} {lab}" for a, lab in s1["against_x4"].items())
    best = f"argmax {s1['argmax']:g}"
    if not any(x in ("ahead", "behind") for x in s1["against_x4"].values()):
        best += f" {NOT_SEPARATED}"
    elif s1["statement"] == STILL_RISING:
        best += f" ({STILL_RISING})"
    return f"{labels}; {best}"


def sha_agreement(results: Mapping[str, Mapping]) -> dict:
    """Per arm, whether every seed's predicted-matrix sha256 is the same on every truth present."""
    out = {}
    for a in ARMS:
        differ, unrecorded = [], []
        for s in SEEDS:
            shas = []
            for r in results.values():
                row = next((x for x in r["rows"] if x["seed"] == s), None)
                shas.append(None if row is None else (row.get("sha256") or {}).get(a))
            if any(x is None for x in shas):
                unrecorded.append(s)
            elif len(set(shas)) > 1:
                differ.append(s)
        out[a] = {"agree": not differ and not unrecorded, "seeds_differing": differ, "seeds_unrecorded": unrecorded}
    return {"arms": out, "all_agree": all(v["agree"] for v in out.values())}


def t50_reading(results: Mapping[str, Mapping], committed: Mapping[str, Mapping]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in results[k]["rows"]])) for a in ARMS}
             for k in TRUTHS}
    over = {k: {a: comparison(results[k], k, a, X4) for a in NEW} for k in TRUTHS}
    turns = {k: turn(means[k]) for k in TRUTHS}
    out["statement_1"] = {k: {"against_x4": {a: over[k][a]["label"] for a in NEW}, "argmax": turns[k]["best"],
                              "statement": turns[k]["statement"]} for k in TRUTHS}
    out["statement_2"] = {k: {f"{c} against {b}": comparison(results[k], k, c, b) for c, b in WORST_CASES}
                          for k in TRUTHS}
    lib3 = results[HEART]
    d33 = paired_difference(results["lib3_t33"], lib3, X4, X2)
    d67 = paired_difference(results["lib3_t67"], lib3, X4, X2)
    null = depth_null(results)

    def said(d):
        return depth_statement(d) if null["holds"] else DEPTH_VOID

    out["statement_3"] = {
        "valid": null["holds"], "draw_null": null,
        "D33": {**d33, "statement": said(d33), "label_if_valid": depth_statement(d33)},
        "D67": {**d67, "statement": said(d67), "label_if_valid": depth_statement(d67)},
        "x8_minus_x4": {k: paired_difference(results[k], lib3, X8, X4) for k in THINNED},
        "reported_not_read": {k: paired_difference(results[k], lib3, X4, X2) for k in ("lib1", "lib2")},
        "shapes": {k: out["statement_1"][k] for k in (HEART, *THINNED)},
        "anchor_gap": {k: anchor_gap(results[k]) for k in (HEART, *THINNED)}}
    out["outcome"] = ("over-delivery: " + "; ".join(f"{k} {outcome_part(out['statement_1'][k])}" for k in TRUTHS)
                      + f"; depth: {out['statement_3']['D33']['statement']}")
    out["over_comparisons"] = over
    out["mean_scores"] = means
    out["turn"] = turns
    out["neighbours"] = {k: neighbours(results[k], k) for k in TRUTHS}
    out["depth_ratio"] = {k: results[k]["thinning"]["depth_ratio"] for k in THINNED}
    out["thinning"] = {k: {x: results[k]["thinning"].get(x) for x in
                           ("p", "stream", "depth_ratio", "rebuild_largest_difference", "any_cell_emptied")}
                       for k in THINNED}
    out["sha256_across_truths"] = sha_agreement(results)
    return out


def _parts(c: Mapping) -> str:
    raw = ", ".join(f"{m} skill {c['skills'].get(m, float('nan')):+.4f} raw {c['metrics'].get(m, float('nan')):+.4g}"
                    for m in ls.WEIGHTS)
    return (f"de terms {c['points']['de']:+.2f}, mmd_u and variogram {c['points']['mmd_u_variogram']:+.2f}; {raw}")


def _comparison_line(c: Mapping) -> str:
    return (f"{c['change']} against {c['base']}: {c['mean']:+.2f} points, ahead {c['seeds_ahead']} of {c['of']}, per "
            f"seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}, per-seed sd {c['sd']:.2f}, test: {c['label']}; "
            + _parts(c))


def _paired_line(d: Mapping) -> str:
    return (f"{d['quantity']}: {d['mean']:+.2f} points, ahead {d['seeds_ahead']} of {d['of']}, per seed "
            f"{d['range'][0]:+.2f} to {d['range'][1]:+.2f}, per-seed sd {d['sd']:.2f}")


def _raw_line(raw: Mapping) -> str:
    if not raw:
        return "not recorded"
    return ", ".join(f"{m} {v['mean']:+.4g} (per seed {v['range'][0]:+.4g} to {v['range'][1]:+.4g})"
                     for m, v in raw.items())


def t50_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"Outcome: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, w in v["reproduction"]["largest_difference"].items():
        lines.append(f"{k}: largest per-seed difference from T46's committed scores: "
                     + ", ".join(f"{a} {x:.3g}" for a, x in w.items()) + f" (gate {REPRO_TOL:g})")
    if r["valid"]:
        lines.append("Statement 1, over-delivery (each new arm against m1_r_gmom_x4 by the single-task test; the "
                     "argmax among 1, 1.5, 2, 3, 4, 5, 6 and 8, smallest on a tie):")
        for k in TRUTHS:
            s1 = r["statement_1"][k]
            lines.append(f"    on {k}: " + ", ".join(f"{a} {lab}" for a, lab in s1["against_x4"].items())
                         + f"; argmax {s1['statement']}"
                         + ("" if any(x in ("ahead", "behind") for x in s1["against_x4"].values())
                            else f" {NOT_SEPARATED}"))
        lines.append("Statement 2, the worst cases a request would quote (read on lib3; reported on the others):")
        for k in (HEART, *(x for x in TRUTHS if x != HEART)):
            for c in r["statement_2"][k].values():
                lines.append(f"    on {k}: " + _comparison_line(c))
        s3 = r["statement_3"]
        null = s3["draw_null"]
        lines.append("Statement 3, self-check (two thinning draws at p = 1/3, the same predictions): "
                     f"{_paired_line(null['x4_minus_x2'])}, test: {null['x4_minus_x2']['label']}; statement 3 "
                     + ("stands" if s3["valid"] else "is void and only reported"))
        lines.append(f"    reported: {_paired_line(null['x8_minus_x4'])}, test: {null['x8_minus_x4']['label']}")
        lines.append(f"Statement 3, depth: D33, {_paired_line(s3['D33'])}, test: {s3['D33']['label']}: "
                     f"{s3['D33']['statement']}; beside it the draw null's mean {null['x4_minus_x2']['mean']:+.2f}, "
                     f"per seed {null['x4_minus_x2']['range'][0]:+.2f} to {null['x4_minus_x2']['range'][1]:+.2f}")
        lines.append(f"    D67, reported: {_paired_line(s3['D67'])}, test: {s3['D67']['label']}: "
                     f"{s3['D67']['statement']}")
        for key in ("D33", "D67"):
            lines.append(f"    {key} on the raw metrics (free of each truth's anchors; read by nothing): "
                         + _raw_line(s3[key]["raw_by_metric"]))
        for d in s3["x8_minus_x4"].values():
            lines.append(f"    reported: {_paired_line(d)}, test: {d['label']}")
        for d in s3["reported_not_read"].values():
            lines.append(f"    beside them, read by nothing (different truths, not paired in kind): {_paired_line(d)}")
        lines.append("    the curve's shape side by side (labels against x4; argmax), independent of each truth's "
                     "scale: " + "; ".join(f"{k} {outcome_part(s)}" for k, s in s3["shapes"].items()))
        for k, g in s3["anchor_gap"].items():
            lines.append(f"    {k}: mean ceiling minus floor per metric (read by nothing): "
                         + ("not recorded" if g is None else ", ".join(f"{m} {x:+.4g}" for m, x in g.items())))
        for k in THINNED:
            th = r["thinning"][k]
            emptied = th["any_cell_emptied"]
            lines.append(f"truth {k}: realised depth ratio {r['depth_ratio'][k]:.4f} (p {THINNED[k]['p']:.4f}, "
                         f"stream {th['stream']}), rebuild check's largest difference "
                         f"{th['rebuild_largest_difference']:.3g}, any cell left with no count: {emptied}")
        for k in TRUTHS:
            lines.append(f"truth {k}: mean scores " + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
            lines.append("    each against its next smaller: " + ", ".join(
                f"{a} against {n['base']} {n['mean']:+.2f} ({n['label']}), {n['per_tenth']:+.2f} per 0.1 of fraction"
                for a, n in r["neighbours"][k].items()))
        sa = r["sha256_across_truths"]
        lines.append("predictions' sha256 across the truths: "
                     + ("the same on every seed and arm" if sa["all_agree"] else
                        "differ: " + "; ".join(f"{a} on seeds {x['seeds_differing'] or x['seeds_unrecorded']}"
                                               for a, x in sa["arms"].items() if not x["agree"])))
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t46_tasks: Path = T46_TASKS) -> dict:
    """Read the truths present in out/tasks against T46's committed task files; write t50.{json,txt}."""
    gv.refuse_existing(out, ("t50.json", "t50.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / TASK_FILE.format(truth=k) for k in TRUTHS
               if (out / "tasks" / TASK_FILE.format(truth=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    files = {k: Path(t46_tasks) / TASK_FILE.format(truth=k) for k in ls.TRUTHS
             if (Path(t46_tasks) / TASK_FILE.format(truth=k)).exists()}
    committed = {k: json.loads(p.read_text()) for k, p in files.items()}
    res = {"label": LABEL_T50, "reading": t50_reading(results, committed),
           "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *files.values()])
    gv.write_outputs(out, "t50", res, t50_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T50's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one truth of try T50")
    p.add_argument("--truth", required=True, choices=list(TRUTHS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T50's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
