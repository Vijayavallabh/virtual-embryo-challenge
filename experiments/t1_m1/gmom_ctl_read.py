#!/usr/bin/env python
"""Try T59 (experiments/harness/TRIES.md): g_mom's larger damping, read for development and for depth. A labelled
in-sample read on the routes of tries T21, T23, T49, T50 and T54. Built to brief AN; its use is fixed in decision log
2026-11-18 (c).

Try T49 found g_mom's within-type factors finished at dampings 0.75 and 1.0 (in place of B0's 0.5) above g_mom on both
counted Task 1 cuts. This read asks two things T49 did not: whether the larger damping's gain needs development and
the assignment of vectors to states (try T23's no-development twin and try T21's permuted control, finished at the
larger dampings), and whether its sign holds where the truth is shallower or deeper than the reference (the worlds of
tries T50 and T54). A labelled in-sample diagnostic, never an adoption result: the damping is a size chosen on the
cuts that kept g_mom, and no outcome is evidence on Task 1's read depth, which stays closed.

  run --cut chamber|whole   Part 1: one counted E10.5 cut, form A (try T11's), 10 seeds, one cell draw per seed. Arms:
         b0, g_mom (0.5), g_mom_d075, g_mom_d1 exactly as try T49 builds them (gmom_size_read.form_a_arms on the fold's
         fit, finished by finish_momentum at 0.5, 0.75 and 1.0);
         g_mom_null (0.5), g_mom_null_d075, g_mom_null_d1: try T23's no-development factors. The fold's cross-fit is
         made once by the registered function itself (harness._g_mom_null_factors, its cache cleared, called with
         g_mom's own finished factors, finish_no_development watched); the two raw directions it hands to
         finish_no_development are then finished again by finish_no_development with only the damping changed: T23's
         finish matches the twin's size per category to g_mom's (the cell-weighted root mean square over the metacells
         and the genes not held flat, by bisection on the damping), so at damping d the target is g_mom's factors
         finished at d, and the record's reference damping is d. At 0.5 the factors are the registered function's own
         object, and the run stops unless a fresh finish at 0.5 equals it bit for bit;
         g_mom_perm (0.5), g_mom_perm_d075, g_mom_perm_d1: try T21's derangement (permute_factors, seed
         harness.G_MOM_PERM_SEED) applied to g_mom's factors finished at each damping.
         Every arm goes through harness._g_mom_arm (g_mom's draws, B0 base, fit and stream [seed, 5]), with a
         transform that checks it was handed g_mom's own factors and returns the arm's.
  run --world lib3|lib3_t67|lib3_t33|ref67|ref50   Part 2: one world of the chamber cut, form B (the entry's), 10
         seeds: b0, m1_r_gmom (0.5), m1_r_gmom_g075, m1_r_gmom_g1 as try T49's form B builds them
         (gmom_size_read.form_b_arms), inside the world made by try T50's or try T54's own code (part_depth_read's
         world_route: over_read.truth_and_contamination for lib3 and the thinned truths, deep_read.world_reference
         for the thinned references, with their own p and streams; in a thinned-reference world the fit is made on the
         thinned reference, as try T54's arms are).
  Refused on an existing task file, unless the committed N22 passed, and unless experiments/harness/TRIES.md holds row
  T59 (try T57's refusals).
  read   T59's reading (below, fixed here and repeated in the ledger row before the run). Writes t59.json and t59.txt
         and refuses to overwrite them.

The reading. Void unless: both cuts and the five worlds hold the ten seeds with every arm and its damping at one
commit on a clean tree; no anchor failure; every arm with a committed score reproduces it to 1e-9 (on the cuts b0,
g_mom, g_mom_d075 and g_mom_d1 seed for seed against try T49's task files; g_mom_null and g_mom_perm against try T23's
and T21's committed record, experiments/harness/results/t20_t24_g_mom_controls/decisions_T1_test.json, which holds
per cut each arm's mean and per-seed sd and its gain over b0 with seeds ahead, and seed for seed as well when the run's
uncommitted task files are given with --t20-tasks; in the worlds b0 and m1_r_gmom seed for seed against try T50's and
T54's task files); the thinned worlds pass the checks of the try they come from (their own p and stream, a realised
depth ratio within 0.01 of p, no cell emptied, the rebuild check; in try T50's worlds every prediction's sha256 equal
to lib3's seed for seed and the truth's row and contamination counts equal to lib3's; in try T54's worlds no arm
scoring exactly as in lib3 on every seed; in every world the recorded M1 fraction is 0.2 within its tolerance, as
tries T50 and T54 check it); the twin's split seed (the registered function's seed default) and the derangement's seed
(the permuted factors' own record) are try T23's and T21's (20261110), the derangement is the same at every damping,
and the factors used at 0.5 for g_mom, the twin and the permuted set have the sha256 of the registered builds' own
factors (harness._G_MOM_PREP's, harness._G_MOM_NULL_PREP's and permute_factors of the former at G_MOM_PERM_SEED); the
twin's size match succeeded in every category at every damping (g_mom's size there lies between the twin's sizes at
the two ends of SIZE_MATCH_DAMPING, so the bisection could reach it, or g_mom has no size there, so both are zero; a
size the clip caps below g_mom's leaves the damping found short of an end, so the damping alone does not show it);
and no control scores exactly as its twin on every seed of a cut (no larger arm of g_mom, the twin or the permuted set
as the same family's arm at 0.5, and no twin or permuted arm as g_mom at the same damping).
Reading 1, per cut and for the dampings 0.75 and 1.0, each paired by seed, with its mean, seeds ahead, per-seed range
and its label by try T46's single-task test ("ahead": a mean paired gain above 0.3 with at least 8 of 10 seeds ahead;
"behind" the same reversed; otherwise "level"): the larger damping's gain G = g_mom_d - g_mom; the twin's Z =
g_mom_null_d - g_mom_null; the development-specific increment D = G - Z; the permuted control's Q = g_mom_perm_d -
g_mom_perm; the assignment-specific increment A = G - Q.
Statement 1 (changed from brief AN after review: the brief's words, "the larger damping's gain needs development and
the assignment", say more than D and A ahead can show, since a twin that carries most of the gain still leaves D
ahead, and its condition never asked that G itself is ahead): "the larger damping's gain is ahead and exceeds the
no-development twin's and the permuted control's gains at g_mom's size by more than 0.3 on both cuts" if at damping
1.0 G, D and A all read ahead on both cuts (six readings); otherwise "not shown", naming each of the six that is level
or behind. Printed beside it, read by nothing: per cut at damping 1.0 the shares Z/G and Q/G of the means (how much of
g_mom's gain the controls reproduce), and per damping the twin's size match's largest relative error.
Reading 2, per world: m1_r_gmom_g1 - m1_r_gmom and m1_r_gmom_g075 - m1_r_gmom, the same way. Statement 2, on the
damping 1.0: "stays positive at every depth read" if in all five worlds the mean is positive with at least 8 of 10
seeds ahead; "it turns" if in some world it is behind by try T46's test, naming it; otherwise "not separated", naming
every world that meets neither. Every world a statement names is printed with its mean, its seeds ahead and the sign
of its mean (the bars are not symmetric; see the limits).
Reported and read by nothing: every arm's mean; each difference's parts on the two de terms and on mmd_u with the
variogram; per damping the target-free records of the finished factors for g_mom, the twin and the permuted set (try
T49's record: share not zeroed, share at the clip, mean absolute log factor), with each family's cell-weighted root
mean square per category and the twin's size match.

Limits (written before the run; each goes into the ledger row and into t59.txt):
- One E8.5 to E9.5 fit serves every seed and both cuts: the seeds measure draw and scoring noise only, n = 1 for the
  fit, and '8 of 10 ahead' is not replication.
- Twice the damping is not twice g_mom: it brings back entries zeroed at 0.5 and moves the clip's share. The twin and
  the permuted set are finished the same way, so D and A compare like with like only as far as the records show their
  sizes agree. The twin is matched to g_mom's size on the cell-weighted root mean square only, not on the share not
  zeroed or the share at the clip, which are printed beside it; the permuted set keeps g_mom's values exactly.
- 'Not shown' is not 'does not need development': it mixes controls that gain as much as g_mom with noise (D and A
  are differences of four scored arms, so their per-seed noise is about sqrt(2) times G's), and a D or A read 'behind'
  means the twin or the permuted set gains more than g_mom. Statement 1 never says the gain needs development: a twin
  that reproduces most of the gain can still leave D ahead, so Z/G and Q/G are printed beside it.
- Statement 1 is read only in try T11's form A on the two cuts and Statement 2 only in the entry's form B in the
  worlds: nothing reads the development-dependence in the entry's form or the depth sign in form A, and Statement 1 is
  not carried into the entry's form. The two cuts share the fit, the reference and the truth's stage, so its six
  readings are closer to three.
- Statement 2's bars are not symmetric: 'stays positive' needs a mean above 0 with 8 of 10 seeds ahead, 'it turns' a
  mean below -0.3 with 8 of 10 behind (try T46's test), so a steady reversal smaller than 0.3 reads 'not separated';
  every world a statement names is printed with its mean, seeds ahead and sign.
- Each thinned world is one draw (the second draws of tries T50 and T54 and their draw-null checks are not used), and
  try T50's three worlds share their predictions bit for bit and differ only in the truth, so they are not independent
  of one another.
- A larger arm redraws more counts: apply_factors draws every count of a call from one stream, so a changed factor
  misaligns every later draw, and larger factors add resampling variance; per-seed differences hold that noise.
- 'Stays positive' has no size threshold: a gain that falls to a few hundredths still stays positive, and each mean is
  printed.
- A thinned reference is one draw, changes the floor, the de terms' reference, B0's shift and g_mom's fit, and is not
  a deeper truth; a thinned truth is one draw shared by all ten seeds. Board points are not comparable across worlds.
- One released day, in-sample: the cuts are the ones that kept g_mom and the worlds are made from the chamber cut's
  library 3; Feng's E9.5 cells lie inside these pseudo-tasks' external windows.
- Nothing here says how deep any hidden truth is, and no statement may be used with a figure for one. No outcome is
  evidence on Task 1's read depth, a lever closed after an exposure that this read does not reopen.
- Whatever it reads, the damping is a size chosen on the cuts that kept g_mom: no outcome changes g_mom's damping in
  any entry or at test, and the result can only go to the team captain.
- The twin's and the permuted control's scores at 0.5 are committed only as per-cut means, per-seed sds and gains over
  b0 (their per-seed task files are not in the repository), so without --t20-tasks their reproduction is checked on
  those aggregates, not seed for seed.

Usage: on the DGX campaign checkout, repo root, every cut and world at one commit, committing nothing until the seven
task files exist, and on the machine tries T49, T50 and T54 ran on, so that the rescored arms reproduce. The
environment variables are docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py run --cut chamber --out t59_gmom_ctl
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py run --cut whole --out t59_gmom_ctl
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py run --world lib3 --out t59_gmom_ctl
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py run --world lib3_t67 --out t59_gmom_ctl
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py run --world lib3_t33 --out t59_gmom_ctl
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py run --world ref67 --out t59_gmom_ctl
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py run --world ref50 --out t59_gmom_ctl
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_ctl_read.py read --out t59_gmom_ctl
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import sys
import time
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
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
import deep_read as dr  # noqa: E402
import gmom_size_read as gs  # noqa: E402
import over_read as ob  # noqa: E402
import part_depth_read as pdr  # noqa: E402

er, gc, ls, feng_dirs = gs.er, gs.gc, ob.ls, gs.feng_dirs
gv, hygiene, n12 = gs.gv, gs.hygiene, gs.n12
ROW = "T59"
STEM = "t59"
CUTS = {"chamber": gc.TASK_CHAMBER, "whole": gc.TASK_WHOLE}
CUT_FILES = gc.CUT_FILES
TASK = ob.TASK  # the chamber cut, every world's
HEART = ob.HEART  # lib3
T50_WORLDS = (HEART, "lib3_t67", "lib3_t33")  # try T50's truths; the inputs as they are
T54_WORLDS = ("ref67", "ref50")  # try T54's thinned references; the truth lib3 as released
WORLDS = (*T50_WORLDS, *T54_WORLDS)
T50_THINNED = tuple(k for k in T50_WORLDS if k != HEART)
THINNING_OF = {**{k: ob.THINNED[k] for k in T50_THINNED}, **{k: dr.THINNED[k] for k in T54_WORLDS}}
SHALLOWER, DEEPER = pdr.SHALLOWER, pdr.DEEPER
SIDE = {**dict.fromkeys(T50_WORLDS, SHALLOWER), **dict.fromkeys(T54_WORLDS, DEEPER)}
TASK_FILE = "{world}.json"

BASE_DAMPING = gs.BASE_DAMPING  # 0.5, g_mom's
LEVELS = gs.LEVELS  # 0.5, 0.75, 1.0
LARGER = (0.75, 1.0)
SUFFIX = {0.5: "", 0.75: "_d075", 1.0: "_d1"}
FAMILIES = ("g_mom", "g_mom_null", "g_mom_perm")
CUT_ARMS = {f"{f}{SUFFIX[d]}": d for f in FAMILIES for d in LEVELS}  # arm -> damping, form A
CUT_ARM_NAMES = ("b0", *CUT_ARMS)
WORLD_ARMS = dict(gs.FORM_B)  # m1_r_gmom, m1_r_gmom_g075, m1_r_gmom_g1 -> multiple of B0's damping
WORLD_DAMPINGS = {a: BASE_DAMPING * m for a, m in WORLD_ARMS.items()}
WORLD_ARM_NAMES = ("b0", *WORLD_ARMS)
WORLD_SUFFIX = {0.5: "", 0.75: "_g075", 1.0: "_g1"}
M1 = gs.BASE_B  # m1_r_gmom
REPRO_TOL = gs.REPRO_TOL  # 1e-9
DAMPING_TOL = gs.DAMPING_TOL
DEPTH_TOL, REBUILD_TOL = ob.DEPTH_TOL, ob.REBUILD_TOL
NULL_SEED, PERM_SEED = harness.G_MOM_NULL_SEED, harness.G_MOM_PERM_SEED  # 20261110, tries T23's and T21's rows
T49_TASKS = RESULTS / "t49_gmom_size" / "tasks"
T49_ARMS = ("b0", "g_mom", "g_mom_d075", "g_mom_d1")  # seed for seed against try T49's files
T20_DECISIONS = REPO_ROOT / "experiments" / "harness" / "results" / "t20_t24_g_mom_controls" / "decisions_T1_test.json"
T20_ARMS = ("g_mom_null", "g_mom_perm")  # tries T23's and T21's candidates
T50_TASKS = ob.RESULTS / "t50_m1_over_step" / "tasks"
T54_TASKS = dr.RESULTS / "t54_m1_deep_step" / "tasks"
WORLD_REPRO_ARMS = ("b0", M1)
FRACTION, FRACTION_TOL = pdr.FRACTION, pdr.FRACTION_TOL  # M1's 0.2 and its tolerance, tries T50's and T54's check
HOLD_AHEAD = 8
EXCEEDS = ("the larger damping's gain is ahead and exceeds the no-development twin's and the permuted control's gains "
           "at g_mom's size by more than 0.3 on both cuts")
NOT_SHOWN = "not shown"
STAYS = "stays positive at every depth read"
TURNS = "it turns"
NOT_SEPARATED = "not separated"
LABEL_T59 = ("try T59, a labelled in-sample diagnostic, never an adoption result: g_mom's within-type factors finished "
             "at dampings 0.75 and 1.0 against B0's 0.5, read against try T23's no-development twin and try T21's "
             "permuted control finished the same way on the two counted E10.5 cuts (try T11's form), and in the "
             "entry's form in the worlds of tries T50 and T54 (library 3's truth as released and thinned to two thirds "
             "and a third; the reference thinned to two thirds and a half); the damping is a size chosen on the cuts "
             "that kept g_mom; the result can only go to the team captain, and no outcome is evidence on Task 1's "
             "read depth, which stays closed")
LIMITS = tuple(" ".join(x.split()) for x in __doc__.split("\nLimits")[1].split("\n\n")[0].split("\n- ")[1:])
CHOICES = (
    "The twin at a larger damping: the registered cross-fit's two raw directions (captured from the one call of "
    "finish_no_development that harness._g_mom_null_factors makes) finished by finish_no_development with the target "
    "g_mom's factors finished at that damping (T23's size match, per category, to g_mom at the same damping) and the "
    "record's reference damping set to that damping (it enters the record only).",
    "The permuted set at a larger damping: permute_factors(g_mom's factors finished at that damping, G_MOM_PERM_SEED); "
    "the derangement depends only on the seed and the metacell counts, so it is the same at every damping (checked).",
    "Tries T21's and T23's per-seed scores are not committed; their reproduction uses the committed aggregates (each "
    "arm's mean and per-seed sd, and its gain over b0 with seeds ahead) to 1e-9, and seed for seed only when the "
    "uncommitted task files are given (--t20-tasks).",
    "'No control scores exactly as its twin' is read two ways, both void reasons: a larger arm of any family scoring "
    "as that family's arm at 0.5 on every seed of a cut, and a twin or permuted arm scoring as g_mom at its damping on "
    "every seed of a cut.",
    "Try T54's 'thinning idle' check is made against lib3 (the inputs as they are and the truth lib3, try T54's ref100 "
    "world), since ref100 is not one of this read's worlds.",
    "The label of D and A is the single-task test on the pair (control's increment, g_mom's increment) per seed, as "
    "try T46's difference-in-differences is tested; the parts of D and A are the differences of the parts.",
    "Besides try T49's records, each family's cell-weighted root mean square per category and damping (the size T23 "
    "matches) and the twin's size-match record are reported, read by nothing.",
    "Statement 1's words and condition (after review): 'exceeds ... at g_mom's size by more than 0.3' in place of "
    "'needs', and G ahead on both cuts at damping 1.0 added to D and A ahead (a tighter bar, never a looser one); Z/G "
    "and Q/G are the ratios of the paired means at damping 1.0 (nan where G's mean is 0).",
    "A failed size match voids: a category whose twin was finished at no damping while g_mom has a size there, or "
    "whose g_mom size lies outside the twin's sizes at the two ends of SIZE_MATCH_DAMPING (the bisection cannot reach "
    "it; where the clip caps the twin's size it keeps the first damping at the cap, short of an end, so a test of the "
    "damping found against the ends would miss it), recorded at run time by watching ts._finish_to_size; no tolerance "
    "is set on the relative error, whose largest value per damping is printed.",
    "The factors used at 0.5 are recorded by sha256 (genes, categories, metacell indexes, sizes, log factors) beside "
    "the registered builds' own (harness._G_MOM_PREP's and _G_MOM_NULL_PREP's at the end of the run, and "
    "permute_factors of the former at harness.G_MOM_PERM_SEED); the twin's seed is the registered function's default "
    "(its one call passes none) and the derangement's is the permuted factors' own record.",
    "Statement 2's bars are kept as the brief sets them (a symmetric bar would loosen 'it turns' or change what 'stays "
    "positive' says); every world a statement names is printed with its mean, seeds ahead and sign instead.",
    "Each world's recorded M1 fraction is checked against 0.2 within part_depth_read's FRACTION_TOL, as tries T50 "
    "and T54 check it.",
)


# ---------------------------------------------------------------------------------------------------------------------
# Part 1: the controls' factors at each damping

_CTL: dict = {}


@contextmanager
def _capture_no_development(captured: list) -> Iterator[None]:
    real = ts.finish_no_development

    def spy(first, second, target, held_flat, reference_damping):
        captured.append({"first": first, "second": second, "target": target, "held_flat": list(held_flat),
                         "reference_damping": reference_damping})
        return real(first, second, target, held_flat, reference_damping)

    ts.finish_no_development = spy
    try:
        yield
    finally:
        ts.finish_no_development = real


@contextmanager
def _capture_size_match(captured: list) -> Iterator[None]:
    """Watch ts._finish_to_size (the twin's size match, called once per category by finish_no_development): per call,
    g_mom's size there and the twin's sizes at the two ends of SIZE_MATCH_DAMPING, computed with the module's own
    _zero_and_clip and _factor_rms. Target-free."""
    real = ts._finish_to_size

    def spy(mean_removed, sizes, free, target_rms):
        ends = [ts._factor_rms(ts._zero_and_clip(mean_removed * float(d)), sizes, free) for d in ts.SIZE_MATCH_DAMPING]
        out = real(mean_removed, sizes, free, target_rms)
        captured.append({"rms_target": float(target_rms), "rms_at_smallest_damping": ends[0],
                         "rms_at_largest_damping": ends[1],
                         "bracketed": bool(ends[0] <= float(target_rms) <= ends[1])})
        return out

    ts._finish_to_size = spy
    try:
        yield
    finally:
        ts._finish_to_size = real


def _by_category(sf: ts.StateFactors, calls: Sequence[Mapping]) -> dict:
    """The size-match calls of one finish_no_development, keyed by the category each finished (they run in the order
    of the finished categories)."""
    if len(calls) != len(sf.categories):
        raise RuntimeError(f"{len(calls)} size matches for {len(sf.categories)} finished categories of the twin")
    return dict(zip(sf.categories, calls, strict=True))


def factors_sha256(sf: ts.StateFactors) -> str:
    """The sha256 of a set of finished factors: genes, categories, metacell indexes, sizes, log factors (what
    gs.same_factors compares)."""
    h = hashlib.sha256(json.dumps(list(map(str, sf.genes))).encode())
    for c, cf in sf.categories.items():
        h.update(json.dumps([str(c), sorted((str(k), int(v)) for k, v in cf.index.items())]).encode())
        for x in (np.asarray(cf.sizes), np.asarray(cf.logf)):
            h.update(str(x.dtype).encode() + str(x.shape).encode() + np.ascontiguousarray(x).tobytes())
    return h.hexdigest()


def control_factors(ctx: harness.BuildContext) -> dict:
    """The fold's factors (try T49's fold_factors) and, at every damping, try T23's twin and try T21's permuted set;
    made once per fold. The twin's cross-fit is made by the registered harness._g_mom_null_factors itself (its cache
    cleared, called with g_mom's own finished factors, finish_no_development watched); its own factors stand at 0.5,
    and the run stops unless exactly one finish was made, at 0.5 and on g_mom's own factors, and a fresh finish at 0.5
    equals it bit for bit. Every size match of the twin is watched (_capture_size_match)."""
    fold = gs.fold_factors(ctx)
    hit = _CTL.get("fold")
    if hit is not None and hit["fold"] is fold:
        return hit
    finished = fold["finished"]
    harness._G_MOM_NULL_PREP.clear()
    captured: list = []
    matches: dict = {}
    calls: list = []
    with _capture_no_development(captured), _capture_size_match(calls):
        own = harness._g_mom_null_factors(ctx, finished[BASE_DAMPING])
    if (len(captured) != 1 or captured[0]["target"] is not finished[BASE_DAMPING]
            or captured[0]["reference_damping"] != BASE_DAMPING):
        raise RuntimeError(f"_g_mom_null_factors finished {len(captured)} times, not once at {BASE_DAMPING:g} on "
                           "g_mom's own factors")
    matches[BASE_DAMPING] = _by_category(own, calls)
    c = captured[0]
    if not gs.same_factors(ts.finish_no_development(c["first"], c["second"], finished[BASE_DAMPING], c["held_flat"],
                                                    BASE_DAMPING), own):
        raise RuntimeError("a fresh finish of the captured cross-fit at 0.5 is not _g_mom_null_factors' own factors")
    null = {BASE_DAMPING: own}
    for d in LEVELS:
        if d == BASE_DAMPING:
            continue
        calls = []
        with _capture_size_match(calls):
            null[d] = ts.finish_no_development(c["first"], c["second"], finished[d], c["held_flat"], d)
        matches[d] = _by_category(null[d], calls)
    null = {d: null[d] for d in LEVELS}
    perm = {d: ts.permute_factors(finished[d], PERM_SEED) for d in LEVELS}
    _CTL.clear()
    hit = _CTL["fold"] = {"fold": fold, "g_mom": finished, "g_mom_null": null, "g_mom_perm": perm,
                          "held_flat": c["held_flat"], "size_match": matches}
    return hit


def at_half_record(ctl: Mapping) -> dict:
    """The identity facts at 0.5: per family the sha256 of the factors the arm used and of the registered build's own
    factors (harness._G_MOM_PREP's and harness._G_MOM_NULL_PREP's as they stand, and permute_factors of the former at
    harness.G_MOM_PERM_SEED, as harness._g_mom_perm makes it); the twin's split seed is the registered function's
    default (its one call passes none) and the derangement's the permuted factors' own record."""
    reg_g = (harness._G_MOM_PREP.get("fold") or {}).get("factors")
    reg_null = (harness._G_MOM_NULL_PREP.get("fold") or {}).get("factors")
    registered = {"g_mom": reg_g, "g_mom_null": reg_null,
                  "g_mom_perm": None if reg_g is None else ts.permute_factors(reg_g, harness.G_MOM_PERM_SEED)}
    sha = {f: {"used": factors_sha256(ctl[f][BASE_DAMPING]),
               "registered": None if registered[f] is None else factors_sha256(registered[f])} for f in FAMILIES}
    null_seed = inspect.signature(harness.g_mom_no_development_factors).parameters["seed"].default
    perm_seed = (ctl["g_mom_perm"][BASE_DAMPING].record.get("permutation") or {}).get("seed")
    return {"sha256": sha, "null_seed": null_seed, "perm_seed": perm_seed}


def control_arms(ctx: harness.BuildContext) -> dict[str, Stage]:
    """g_mom_null and g_mom_perm at 0.5, 0.75 and 1.0: `harness._g_mom_arm` with a transform that checks it was
    handed g_mom's own factors and returns the control's factors at the arm's damping."""
    ctl = control_factors(ctx)
    fold = ctl["fold"]

    def at(family: str, d: float) -> Callable[[ts.StateFactors], ts.StateFactors]:
        def transform(got: ts.StateFactors) -> ts.StateFactors:
            gs._base_factors(fold, got)
            return ctl[family][d]
        return transform

    return {f"{f}{SUFFIX[d]}": harness._g_mom_arm(ctx, at(f, d)) for f in ("g_mom_null", "g_mom_perm") for d in LEVELS}


def cut_arms(ctx: harness.BuildContext, b0) -> dict[str, Stage]:
    """One seed's arms of part 1: b0, try T49's form A arms, the controls."""
    a = gs.form_a_arms(ctx)
    return {"b0": b0(ctx), **a, **control_arms(ctx)}


def size_records(by_d: Mapping[float, ts.StateFactors], held_flat: Sequence[str]) -> dict:
    """Per damping and category, the cell-weighted root mean square over the metacells and the free genes (the size
    try T23 matches). Target-free; no count."""
    out = {}
    for d, sf in sorted(by_d.items()):
        free = np.array([g not in set(map(str, held_flat)) for g in sf.genes], dtype=bool)
        out[f"{d:g}"] = {c: ts._factor_rms(cf.logf, cf.sizes, free) for c, cf in sf.categories.items()}
    return out


def control_records(ctl: Mapping, held_flat: Sequence[str]) -> dict:
    """The target-free records of each family's finished factors (try T49's, per damping), their sizes, the twin's
    size match per category and damping (finish_no_development's record and, where watched, g_mom's size against the
    twin's sizes at the two ends of SIZE_MATCH_DAMPING), and the derangement per damping."""
    watched = ctl.get("size_match") or {}
    match = {f"{d:g}": {c: {**cf.record.get("no_development", {}), **((watched.get(d) or {}).get(c) or {})}
                        for c, cf in ctl["g_mom_null"][d].categories.items()} for d in LEVELS}
    perms = {f"{d:g}": (ctl["g_mom_perm"][d].record.get("permutation") or {}) for d in LEVELS}
    return {"factors": {f: gs.factor_records(ctl[f], held_flat) for f in FAMILIES},
            "sizes": {f: size_records(ctl[f], held_flat) for f in FAMILIES},
            "null_match": match,
            "null_categories": {f"{d:g}": sorted(ctl["g_mom_null"][d].categories) for d in LEVELS},
            "permutation": perms,
            "permutation_same_at_every_damping": all(perms[f"{d:g}"] == perms[f"{BASE_DAMPING:g}"] for d in LEVELS)}


# ---------------------------------------------------------------------------------------------------------------------
# The run


def _refusals(a) -> Path:
    if ROW not in n12.tries_rows():
        sys.exit(f"T59 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T59 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    return n22


def _write(out: Path, name: str, res: Mapping) -> None:
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_cut(a) -> None:
    n22 = _refusals(a)
    key = CUTS[a.cut]
    out = RESULTS / a.out
    name = CUT_FILES[key]
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[key]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, _ = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                 (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T59, "part": "cut", "task": t.key, "seeds": list(seeds), "dampings": dict(CUT_ARMS),
           "reproduction_tolerance": REPRO_TOL}
    floors, failures = ls.all_anchors(t, board, target, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        _write(out, name, res)
        sys.exit(f"T59 is refused on {t.key}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T59 is void)")
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    preregistered = t1_scope.held_flat_genes()
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    gs._FOLD.clear()
    _CTL.clear()
    harness._G_MOM_PREP.clear()
    harness._G_MOM_NULL_PREP.clear()
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = cut_arms(ctx, b0)
        arms = er.score_arms(board, target, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    ctl = _CTL["fold"]
    half = at_half_record(ctl)
    res.update({"null_seed": half["null_seed"], "perm_seed": half["perm_seed"], "at_half": half["sha256"],
                "records": control_records(ctl, ctl["held_flat"]),
                "limits": list(LIMITS), "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, n22])})
    _write(out, name, res)


def cmd_run_world(a) -> None:
    n22 = _refusals(a)
    out = RESULTS / a.out
    name = TASK_FILE.format(world=a.world)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    role = harness.base_role(t)
    try:  # a thinned reference comes before anything else is made (try T54's order)
        inputs, ref_stage, truth, contamination, thinning = pdr.world_route(
            a.world, t, data, genes, inputs, role, ref_stage, lambda: read_stage(data / t.reference, None))
    except ValueError as e:
        sys.exit(f"T59 is refused on {a.world}: {e}")
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T59, "part": "world", "task": t.key, "world": a.world,
           "truth": a.world if a.world in T50_WORLDS else HEART, "truth_rows": int(truth.n),
           "truth_contamination": contamination, "seeds": list(seeds), "thinning": thinning,
           "dampings": dict(WORLD_DAMPINGS), "reproduction_tolerance": REPRO_TOL}
    floors, failures = ls.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        _write(out, name, res)
        sys.exit(f"T59 is refused on {a.world}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T59 is void)")
    pool = inputs[role]
    directions, ext = feng_dirs.feng_directions(list(pool.var_names))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    gs._FOLD.clear()
    harness._G_MOM_PREP.clear()
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = {"b0": b0(ctx), **gs.form_b_arms(ctx, directions, groups, fraction, held_m1, preregistered)}
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "anchors": ob.anchor_record(f.anchors), "arms": arms,
                     "sha256": {k: ls.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: ls.sparsity(p) for k, p in preds.items()}, "seconds": round(time.time() - t0, 1)})
        print(f"{a.world} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res.update({"fraction": fraction, "external": ext,
                "factor_records": gs.factor_records(gs._FOLD["fold"]["finished"], preregistered),
                "limits": list(LIMITS), "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    _write(out, name, res)


# ---------------------------------------------------------------------------------------------------------------------
# The reading: validity


def _rows(r: Mapping) -> list:
    return list(r.get("rows") or [])


def _per_seed(rows: Sequence[Mapping], committed: Mapping[int, float], arm: str, absent: list, tag: str) -> float:
    worst = 0.0
    for r in rows:
        old = committed.get(r["seed"])
        if old is None or arm not in r["arms"]:
            absent.append(f"{tag}:{r['seed']}:{arm}")
            continue
        worst = max(worst, abs(float(r["arms"][arm]["score"]) - float(old)))
    return worst


def cut_reproduction(key: str, rows: Sequence[Mapping], t49: Mapping | None, t20: Mapping | None,
                     t20_rows: Mapping | None) -> dict:
    """b0 and form A against try T49's file seed for seed; g_mom_null and g_mom_perm against try T20-T24's committed
    aggregates (mean, per-seed sd, gain over b0 and seeds ahead) and, when given, its task file seed for seed."""
    absent: list = []
    worst: dict = {}
    t49_by = {r["seed"]: r for r in _rows(t49 or {})}
    for arm in T49_ARMS:
        old = {s: (r.get("arms", {}).get(arm) or {}).get("score") for s, r in t49_by.items()}
        worst[arm] = _per_seed(rows, {s: v for s, v in old.items() if v is not None}, arm, absent, "t49")
    aggregate: dict = {}
    summaries = ((t20 or {}).get("summaries") or {}).get(key) or {}
    decisions = {(d.get("base"), d.get("change")): d for d in (t20 or {}).get("decisions", [])}
    for arm in T20_ARMS:
        sc = [r["arms"][arm]["score"] for r in rows if arm in r["arms"]]
        gain = [r["arms"][arm]["score"] - r["arms"]["b0"]["score"] for r in rows if arm in r["arms"]
                and "b0" in r["arms"]]
        summ = summaries.get(arm)
        dec = (decisions.get(("b0", arm)) or {}).get("tasks", {}).get(key)
        if len(sc) != len(SEEDS) or summ is None or dec is None:
            absent.append(f"t20:{arm}")
            aggregate[arm] = None
            continue
        rec = {"mean": abs(float(np.mean(sc)) - float(summ[0])), "sd": abs(float(np.std(sc, ddof=1)) - float(summ[1])),
               "gain_over_b0": abs(float(np.mean(gain)) - float(dec["mean_gain"])),
               "seeds_ahead": int(sum(g > 0 for g in gain)), "seeds_ahead_committed": dec.get("seeds_ahead")}
        aggregate[arm] = rec
        worst[arm] = max(rec["mean"], rec["sd"], rec["gain_over_b0"])
        if rec["seeds_ahead"] != rec["seeds_ahead_committed"]:
            absent.append(f"t20:{arm}: seeds ahead of b0 {rec['seeds_ahead']}, committed "
                          f"{rec['seeds_ahead_committed']}")
        if t20_rows is not None:
            old = {r["seed"]: (r.get("scores") or {}).get(arm) for r in _rows(t20_rows)}
            worst[arm] = max(worst[arm], _per_seed(rows, {s: v for s, v in old.items() if v is not None}, arm,
                                                   absent, "t20 rows"))
    holds = bool(rows) and bool(t49_by) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "aggregates": aggregate, "absent": absent, "holds": holds,
            "seed_for_seed_t21_t23": t20_rows is not None}


def world_reproduction(rows: Sequence[Mapping], committed: Mapping | None) -> dict:
    absent: list = []
    old_by = {r["seed"]: r for r in _rows(committed or {})}
    worst = {}
    for arm in WORLD_REPRO_ARMS:
        old = {s: (r.get("arms", {}).get(arm) or {}).get("score") for s, r in old_by.items()}
        worst[arm] = _per_seed(rows, {s: v for s, v in old.items() if v is not None}, arm, absent, "committed")
    holds = bool(rows) and bool(old_by) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "absent": absent, "holds": holds}


def _dampings_off(got: Mapping | None, want: Mapping[str, float]) -> list[str]:
    got = got or {}
    return [a for a, d in want.items() if not isinstance(got.get(a), int | float) or abs(got[a] - d) > DAMPING_TOL]


def _common_reasons(k: str, r: Mapping, arms: Sequence[str], dampings: Mapping[str, float]) -> list[str]:
    reasons = []
    if r.get("anchor_failures"):
        reasons.append(f"{k}: the floor's anchors invert or go missing on seeds "
                       f"{[f['seed'] for f in r['anchor_failures']]}")
    seeds = [row["seed"] for row in _rows(r)]
    if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
        reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
    absent = sorted({a for row in _rows(r) for a in arms if a not in row["arms"]})
    if absent:
        reasons.append(f"{k}: arms {absent} missing")
    off = _dampings_off(r.get("dampings"), dampings)
    if off:
        reasons.append(f"{k}: dampings {r.get('dampings')}, not {dict(dampings)} (wrong damping on {off})")
    return reasons


def _complete(r: Mapping, arms: Sequence[str]) -> bool:
    return bool(_rows(r)) and all(a in row["arms"] for row in _rows(r) for a in arms)


def _same_on_every_seed(r: Mapping, a: str, b: str) -> bool:
    return all(row["arms"][a] == row["arms"][b] for row in _rows(r))


def cut_identity_reasons(k: str, r: Mapping) -> list[str]:
    """No control scores exactly as its twin on every seed of a cut: no larger arm as its family's arm at 0.5, and no
    twin or permuted arm as g_mom at its damping."""
    if not _complete(r, CUT_ARM_NAMES):
        return []
    same = [f"{f}{SUFFIX[d]} = {f}" for f in FAMILIES for d in LARGER
            if _same_on_every_seed(r, f"{f}{SUFFIX[d]}", f)]
    same += [f"{f}{SUFFIX[d]} = g_mom{SUFFIX[d]}" for f in ("g_mom_null", "g_mom_perm") for d in LEVELS
             if _same_on_every_seed(r, f"{f}{SUFFIX[d]}", f"g_mom{SUFFIX[d]}")]
    if not same:
        return []
    return [f"{k}: an arm scores exactly as its twin on every seed (a damping or a control idle?): {same}"]


def record_reasons(k: str, r: Mapping) -> list[str]:
    """The twin's and the derangement's seeds are tries T23's and T21's, the derangement is the same at every damping,
    the factors used at 0.5 have the registered builds' sha256, and the twin's size match succeeded in every category
    at every damping."""
    reasons = []
    if r.get("null_seed") != NULL_SEED or r.get("perm_seed") != PERM_SEED:
        reasons.append(f"{k}: the twin's and the derangement's seeds are {r.get('null_seed')!r} and "
                       f"{r.get('perm_seed')!r}, not tries T23's and T21's ({NULL_SEED} and {PERM_SEED})")
    half = r.get("at_half") or {}
    off = [f for f in FAMILIES if not (isinstance((half.get(f) or {}).get("used"), str)
                                       and half[f]["used"] == half[f].get("registered"))]
    if off:
        reasons.append(f"{k}: the factors used at 0.5 are not recorded with the registered candidates' sha256 on {off}")
    rec = r.get("records") or {}
    reasons += size_match_reasons(k, rec.get("null_match"))
    if rec.get("permutation_same_at_every_damping") is not True:
        reasons.append(f"{k}: the derangement is not recorded as the same at every damping")
    perm = (rec.get("permutation") or {}).get(f"{BASE_DAMPING:g}") or {}
    if perm.get("seed") != PERM_SEED:
        reasons.append(f"{k}: the derangement's record holds seed {perm.get('seed')!r}, not {PERM_SEED}")
    cats = rec.get("null_categories") or {}
    if not cats or len({json.dumps(v) for v in cats.values()}) != 1:
        reasons.append(f"{k}: the twin's categories are not recorded as the same at every damping")
    return reasons


def size_match_reasons(k: str, match: Mapping | None) -> list[str]:
    """The twin's size match per damping and category: void where it is not recorded, where the twin was finished at
    no damping while g_mom has a size there, or where g_mom's size lies outside the twin's sizes at the two ends of
    SIZE_MATCH_DAMPING (the bisection cannot reach it). Where g_mom has no size (0) both are zero; that is a match."""
    match = match or {}
    if sorted(match) != sorted(f"{d:g}" for d in LEVELS) or not all(match.values()):
        return [f"{k}: the twin's size match is not recorded at every damping"]
    failed = []
    for d, cats in match.items():
        for c, m in cats.items():
            target = m.get("rms_target")
            if isinstance(target, int | float) and target == 0:
                continue
            if m.get("damping") is None or m.get("bracketed") is not True:
                failed.append(f"{c} at {d} (damping found {m.get('damping')!r}, g_mom's size within the twin's "
                              f"range {m.get('bracketed')!r})")
    if not failed:
        return []
    return [f"{k}: the twin's size match failed (no damping, or g_mom's size outside the twin's range over "
            f"SIZE_MATCH_DAMPING): {failed}"]


def thinning_reasons(k: str, th: Mapping | None) -> list[str]:
    """Tries T50's and T54's own checks of a thinned world's record."""
    spec, th = THINNING_OF[k], th or {}
    reasons = []
    if th.get("p") != spec["p"] or list(th.get("stream") or []) != list(spec["stream"]):
        reasons.append(f"{k}: the file records p {th.get('p')!r} and stream {th.get('stream')!r}, not {spec['p']!r} "
                       f"and {spec['stream']}")
    ratio = th.get("depth_ratio")
    if not isinstance(ratio, int | float) or not abs(ratio - spec["p"]) <= DEPTH_TOL:
        reasons.append(f"{k}: realised depth ratio {ratio!r} is not within {DEPTH_TOL:g} of p {spec['p']:.4f}")
    if th.get("any_cell_emptied") is not False:
        reasons.append(f"{k}: the thinning left a cell with no count (or did not record whether it did): "
                       f"any_cell_emptied {th.get('any_cell_emptied')!r}")
    rebuild = th.get("rebuild_largest_difference")
    if not isinstance(rebuild, int | float) or not rebuild <= REBUILD_TOL:
        reasons.append(f"{k}: the rebuild check's largest difference {rebuild!r} is not at most {REBUILD_TOL:g}")
    return reasons


def world_reasons(worlds: Mapping[str, Mapping]) -> list[str]:
    reasons = []
    lib3 = worlds.get(HEART)
    for k, r in worlds.items():
        if r.get("task") != TASK:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("world") != k:
            reasons.append(f"{k}: the file holds world {r.get('world')!r}")
        want_truth = k if k in T50_WORLDS else HEART
        if r.get("truth") != want_truth:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}, not {want_truth}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g} (tries T50's and T54's check)")
        if k == HEART:
            if r.get("thinning") is not None:
                reasons.append(f"{k}: the file records a thinning")
            continue
        reasons += thinning_reasons(k, r.get("thinning"))
        rows = _rows(r)
        if k in T50_THINNED:  # try T50's checks: the truth changes only
            lib3_rows = {row["seed"]: row for row in _rows(lib3 or {})}
            differ = sorted({f"{row['seed']}:{a}" for row in rows for a in WORLD_ARM_NAMES
                             if (row.get("sha256") or {}).get(a) is None
                             or (lib3_rows.get(row["seed"], {}).get("sha256") or {}).get(a) != row["sha256"][a]})
            if lib3 is None or differ:
                reasons.append(f"{k}: the predictions' sha256 do not equal this run's lib3 file seed for seed "
                               f"(seed:arm {differ})")
            for key in ("truth_rows", "truth_contamination"):
                if lib3 is None or not (ob._count(r.get(key)) and r.get(key) == lib3.get(key)):
                    reasons.append(f"{k}: {key.replace('_', ' ')} do not equal lib3's (try T50's check)")
        else:  # try T54's check: the thinning is not idle (against lib3, this read's unthinned world)
            lib3_rows = {row["seed"]: row for row in _rows(lib3 or {})}
            if lib3_rows and _complete(r, WORLD_ARM_NAMES):
                same = [a for a in WORLD_ARM_NAMES
                        if all(row["seed"] in lib3_rows and a in lib3_rows[row["seed"]]["arms"]
                               and row["arms"][a]["score"] == lib3_rows[row["seed"]]["arms"][a]["score"]
                               for row in rows)]
                if same:
                    reasons.append(f"{k}: arms score exactly as in {HEART} on every seed (thinning idle?): {same}")
    return reasons


def validity(cuts: Mapping[str, Mapping], worlds: Mapping[str, Mapping], committed: Mapping) -> dict:
    """The void conditions. cuts: cut key -> task file; worlds: world -> task file; committed: {"t49": cut -> file,
    "t20": the decisions record, "t20_rows": cut -> uncommitted task file (optional), "worlds": world -> try T50's or
    T54's file}."""
    reasons = []
    missing = [k for k in CUT_FILES if k not in cuts] + [k for k in WORLDS if k not in worlds]
    if missing:
        reasons.append(f"missing cuts or worlds {missing}")
    files = {**cuts, **worlds}
    for k, r in cuts.items():
        reasons += _common_reasons(k, r, CUT_ARM_NAMES, CUT_ARMS)
        if r.get("task") != k:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        reasons += record_reasons(k, r)
        reasons += cut_identity_reasons(k, r)
    for k, r in worlds.items():
        reasons += _common_reasons(k, r, WORLD_ARM_NAMES, WORLD_DAMPINGS)
    reasons += world_reasons(worlds)
    commits = {(r.get("stamp") or {}).get("commit") for r in files.values()}
    unclean = sorted(k for k, r in files.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the cuts and worlds ran at {len(commits)} commits")
    repro = {}
    for k, r in cuts.items():
        repro[k] = cut_reproduction(k, _rows(r), (committed.get("t49") or {}).get(k), committed.get("t20"),
                                    (committed.get("t20_rows") or {}).get(k))
        if not repro[k]["holds"]:
            reasons.append(f"{k}: {', '.join(T49_ARMS)} (try T49's) or {', '.join(T20_ARMS)} (tries T23's and T21's) "
                           f"do not reproduce their committed scores within {REPRO_TOL:g}")
    for k, r in worlds.items():
        repro[k] = world_reproduction(_rows(r), (committed.get("worlds") or {}).get(k))
        if not repro[k]["holds"]:
            src = "try T50's" if k in T50_WORLDS else "try T54's"
            reasons.append(f"{k}: {', '.join(WORLD_REPRO_ARMS)} do not reproduce {src} committed scores within "
                           f"{REPRO_TOL:g} per seed")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


# ---------------------------------------------------------------------------------------------------------------------
# The reading: the differences and the statements


def paired(by_seed: Mapping[int, tuple[float, float]]) -> dict:
    """Per seed (base, change): the difference change - base with its mean, seeds ahead, per-seed range and sd, and
    try T46's single-task test (vec.adoption.keep_change on this one comparison, and reversed)."""
    fwd = keep_change({"d": dict(by_seed)})["tasks"]["d"]
    rev = keep_change({"d": {s: (y, x) for s, (x, y) in by_seed.items()}})["tasks"]["d"]
    d = np.array([y - x for _, (x, y) in sorted(by_seed.items())], dtype=float)
    ahead, behind = bool(fwd.get("passes")), bool(rev.get("passes"))
    return {"mean": float(d.mean()) if len(d) else float("nan"), "seeds_ahead": int((d > 0).sum()), "of": len(d),
            "range": [float(d.min()), float(d.max())] if len(d) else [float("nan")] * 2,
            "sd": float(d.std(ddof=1)) if len(d) > 1 else float("nan"),
            "ahead": ahead, "behind": behind, "label": "ahead" if ahead else "behind" if behind else "level"}


def _score(row: Mapping, arm: str) -> float:
    return float(row["arms"][arm]["score"])


def increment(rows: Sequence[Mapping], base: str, change: str) -> dict:
    """change - base, paired by seed, with its parts (try T46's points by part)."""
    out = paired({row["seed"]: (_score(row, base), _score(row, change)) for row in rows})
    return {**out, "quantity": f"{change} - {base}", "points": ls.points_by_part(rows, base, change)}


def _parts_minus(a: Mapping, b: Mapping) -> dict:
    return {"by_metric": {m: a["by_metric"][m] - b["by_metric"][m] for m in a["by_metric"]},
            "de": a["de"] - b["de"], "mmd_u_variogram": a["mmd_u_variogram"] - b["mmd_u_variogram"]}


def difference_of_increments(rows: Sequence[Mapping], real: tuple[str, str], control: tuple[str, str]) -> dict:
    """(real change - real base) - (control change - control base), paired by seed; the test on the pair (control's
    increment, real's increment); the parts are the difference of the two increments' parts."""
    def inc(row, pair):
        return _score(row, pair[1]) - _score(row, pair[0])

    out = paired({row["seed"]: (inc(row, control), inc(row, real)) for row in rows})
    parts = _parts_minus(ls.points_by_part(rows, *real), ls.points_by_part(rows, *control))
    return {**out, "quantity": f"({real[1]} - {real[0]}) - ({control[1]} - {control[0]})", "points": parts}


def reading_one(cuts: Mapping[str, Mapping]) -> dict:
    """Per cut and larger damping: G, Z, D = G - Z, Q, A = G - Q."""
    out = {}
    for k, r in cuts.items():
        rows = _rows(r)
        out[k] = {}
        for d in LARGER:
            s = SUFFIX[d]
            g, z, q = ("g_mom", f"g_mom{s}"), ("g_mom_null", f"g_mom_null{s}"), ("g_mom_perm", f"g_mom_perm{s}")
            out[k][f"{d:g}"] = {"G": increment(rows, *g), "Z": increment(rows, *z),
                                "D": difference_of_increments(rows, g, z), "Q": increment(rows, *q),
                                "A": difference_of_increments(rows, g, q)}
    return out


def _share(part: float, whole: float) -> float:
    return float(part) / float(whole) if whole != 0 else float("nan")


def statement_one(r1: Mapping) -> dict:
    """At damping 1.0: G, D and A ahead on both cuts, or 'not shown' naming each of the six that is level or behind.
    Beside it, read by nothing: per cut the shares Z/G and Q/G of the paired means (how much of g_mom's gain the
    controls reproduce), so that 'exceeds' is never read as 'needs'."""
    six = {f"{q} on {k}": r1[k]["1"][q]["label"] for k in CUT_FILES for q in ("G", "D", "A")}
    short = [f"{name} {lab}" for name, lab in six.items() if lab != "ahead"]
    text = EXCEEDS if not short else f"{NOT_SHOWN} ({'; '.join(short)})"
    shares = {k: {"Z/G": _share(r1[k]["1"]["Z"]["mean"], r1[k]["1"]["G"]["mean"]),
                  "Q/G": _share(r1[k]["1"]["Q"]["mean"], r1[k]["1"]["G"]["mean"])} for k in CUT_FILES}
    return {"statement": text, "readings": six, "holds": not short, "control_share_of_gain": shares}


def size_match_errors(cuts: Mapping[str, Mapping]) -> dict:
    """Per cut and damping the twin's size match's largest relative error over the categories (None where none is
    recorded). Read by nothing; printed beside statement 1."""
    out = {}
    for k, r in cuts.items():
        match = (r.get("records") or {}).get("null_match") or {}
        out[k] = {d: max((m["relative_error"] for m in cats.values()
                          if isinstance(m.get("relative_error"), int | float)), default=None)
                  for d, cats in match.items()}
    return out


def reading_two(worlds: Mapping[str, Mapping]) -> dict:
    return {k: {f"{d:g}": increment(_rows(worlds[k]), M1, f"{M1}{WORLD_SUFFIX[d]}") for d in LARGER} for k in WORLDS}


def _stays(c: Mapping) -> bool:
    return c["mean"] > 0 and c["seeds_ahead"] >= HOLD_AHEAD


def _sign(x: float) -> str:
    return "positive" if x > 0 else "negative" if x < 0 else "zero"


def _named(k: str, c: Mapping) -> str:
    return (f"{k}, {SIDE[k]}: mean {c['mean']:+.2f}, {c['seeds_ahead']} of {c['of']} seeds ahead, "
            f"{_sign(c['mean'])}")


def statement_two(r2: Mapping) -> dict:
    """On the damping 1.0: stays positive in all five worlds; it turns where a world is behind; otherwise not
    separated, naming every world that meets neither. Each world named is printed with its mean, seeds ahead and the
    sign of its mean, since the bars are not symmetric (a mean above 0 with 8 of 10 ahead stays; it turns only below
    -0.3 with 8 of 10 behind)."""
    at1 = {k: r2[k]["1"] for k in WORLDS}
    behind = [k for k, c in at1.items() if c["behind"]]
    neither = [k for k, c in at1.items() if not _stays(c) and not c["behind"]]
    if all(_stays(c) for c in at1.values()):
        text = STAYS
    elif behind:
        text = f"{TURNS} ({'; '.join(_named(k, at1[k]) for k in behind)})"
    else:
        text = f"{NOT_SEPARATED} ({'; '.join(_named(k, at1[k]) for k in neither)})"
    return {"statement": text, "behind": behind, "neither": neither,
            "signs": {k: _sign(c["mean"]) for k, c in at1.items()}}


def t59_reading(cuts: Mapping[str, Mapping], worlds: Mapping[str, Mapping], committed: Mapping) -> dict:
    v = validity(cuts, worlds, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    r1, r2 = reading_one(cuts), reading_two(worlds)
    s1, s2 = statement_one(r1), statement_two(r2)
    out.update({"outcome": f"statement 1: {s1['statement']}; statement 2: {s2['statement']}",
                "reading_1": r1, "statement_1": s1, "size_match_largest_relative_error": size_match_errors(cuts),
                "reading_2": r2, "statement_2": s2})
    out["mean_scores"] = {k: {a: float(np.mean([_score(row, a) for row in _rows(r)])) for a in arms}
                          for k, r, arms in [*((k, cuts[k], CUT_ARM_NAMES) for k in CUT_FILES),
                                             *((k, worlds[k], WORLD_ARM_NAMES) for k in WORLDS)]}
    out["records"] = {k: cuts[k].get("records") for k in CUT_FILES}
    out["world_factor_records"] = {k: worlds[k].get("factor_records") for k in WORLDS}
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The text and the read


def _line(c: Mapping) -> str:
    p = c["points"]
    return (f"{c['quantity']}: {c['mean']:+.2f} points, ahead {c['seeds_ahead']} of {c['of']}, per seed "
            f"{c['range'][0]:+.2f} to {c['range'][1]:+.2f}, test: {c['label']}; de terms {p['de']:+.2f}, mmd_u and "
            f"variogram {p['mmd_u_variogram']:+.2f}")


def t59_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"Outcome: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, g in v["reproduction"].items():
        lines.append(f"{k}: largest difference from the committed scores: "
                     + ", ".join(f"{a} {x:.3g}" for a, x in g["largest_difference"].items()) + f" (gate {REPRO_TOL:g}"
                     + (")" if "seed_for_seed_t21_t23" not in g else
                        "; g_mom_null and g_mom_perm " + ("seed for seed)" if g["seed_for_seed_t21_t23"]
                                                          else "on the committed aggregates)")))
    if r["valid"]:
        lines.append(f"Statement 1: {r['statement_1']['statement']}")
        lines.append("    beside it, read by nothing ('exceeds' is not 'needs'): at damping 1 the controls' gains as "
                     "shares of g_mom's: " + "; ".join(
                         f"{k} Z/G {x['Z/G']:.2f}, Q/G {x['Q/G']:.2f}"
                         for k, x in r["statement_1"]["control_share_of_gain"].items()))
        lines.append("    the twin's size match, largest relative error per damping: " + "; ".join(
            f"{k} " + ", ".join(f"{d} {e if e is None else format(e, '.3g')}" for d, e in by_d.items())
            for k, by_d in r["size_match_largest_relative_error"].items()))
        for k, by_d in r["reading_1"].items():
            for d, qs in by_d.items():
                for q, c in qs.items():
                    lines.append(f"    {k}, damping {d}, {q}: " + _line(c))
        lines.append(f"Statement 2 (damping 1): {r['statement_2']['statement']}")
        for k, by_d in r["reading_2"].items():
            for d, c in by_d.items():
                lines.append(f"    {k} ({SIDE[k]}), damping {d}: " + _line(c))
        lines.append("Reported and read by nothing:")
        for k, ms in r["mean_scores"].items():
            lines.append(f"    mean scores {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in ms.items()))
        for k, rec in r["records"].items():
            for fam, by_d in (rec or {}).get("factors", {}).items():
                for d, cats in by_d.items():
                    size = ((rec.get("sizes") or {}).get(fam) or {}).get(d) or {}
                    lines.append(f"    factors {k} {fam} damping {d}: " + "; ".join(
                        f"{c} not zeroed {x['share_not_zeroed']:.3f}, at clip {x['share_at_clip']:.4f}, mean |log f| "
                        f"{x['mean_abs_log_not_zeroed']:.4f} not zeroed, {x['mean_abs_log_all_entries']:.4f} all "
                        f"entries, rms {size.get(c, float('nan')):.4f}" for c, x in cats.items()))
            for d, cats in ((rec or {}).get("null_match") or {}).items():
                lines.append(f"    twin's size match {k} damping {d}: " + "; ".join(
                    f"{c} damping found {m.get('damping')}, relative error {m.get('relative_error')}, g_mom's size "
                    f"within the twin's range {m.get('bracketed')}"
                    for c, m in cats.items()))
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.extend(f"Choice: {x}" for x in CHOICES)
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t49_tasks: Path = T49_TASKS, t20_decisions: Path = T20_DECISIONS, t50_tasks: Path = T50_TASKS,
         t54_tasks: Path = T54_TASKS, t20_tasks: Path | None = None) -> dict:
    """Read the cuts and worlds present in out/tasks against the committed files of tries T49, T20 to T24, T50 and
    T54 (and, when given, the T20 to T24 run's task files); write t59.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    cut_paths = {k: out / "tasks" / f for k, f in CUT_FILES.items() if (out / "tasks" / f).exists()}
    world_paths = {k: out / "tasks" / TASK_FILE.format(world=k) for k in WORLDS
                   if (out / "tasks" / TASK_FILE.format(world=k)).exists()}
    cuts = {k: json.loads(p.read_text()) for k, p in cut_paths.items()}
    worlds = {k: json.loads(p.read_text()) for k, p in world_paths.items()}
    paths = []

    def load(p: Path) -> dict | None:
        if not Path(p).exists():
            return None
        paths.append(Path(p))
        return json.loads(Path(p).read_text())

    committed = {"t49": {k: load(Path(t49_tasks) / f) for k, f in CUT_FILES.items()},
                 "t20": load(t20_decisions),
                 "t20_rows": None if t20_tasks is None else {k: load(Path(t20_tasks) / f)
                                                             for k, f in CUT_FILES.items()},
                 "worlds": {k: load(Path(t50_tasks if k in T50_WORLDS else t54_tasks) / TASK_FILE.format(world=k))
                            for k in WORLDS}}
    res = {"label": LABEL_T59, "limits": list(LIMITS), "choices": list(CHOICES),
           "reading": t59_reading(cuts, worlds, committed),
           "stamps": {k: v.get("stamp") for k, v in {**cuts, **worlds}.items()}}
    res["stamp"] = gv.stamp(started, [*cut_paths.values(), *world_paths.values(), *paths])
    gv.write_outputs(out, STEM, res, t59_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out, t20_tasks=None if a.t20_tasks is None else Path(a.t20_tasks))
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T59's reading is void and goes to the captain")


def cmd_run(a) -> None:
    if (a.cut is None) == (a.world is None):
        sys.exit("T59's run takes exactly one of --cut and --world")
    (cmd_run_cut if a.cut is not None else cmd_run_world)(a)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one cut (part 1) or one world (part 2) of try T59")
    p.add_argument("--cut", choices=list(CUTS))
    p.add_argument("--world", choices=list(WORLDS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T59's reading")
    p.add_argument("--out", required=True)
    p.add_argument("--t20-tasks", default=None,
                   help="the T20 to T24 run's task files (uncommitted), for a seed-for-seed check of g_mom_null and "
                        "g_mom_perm")
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
