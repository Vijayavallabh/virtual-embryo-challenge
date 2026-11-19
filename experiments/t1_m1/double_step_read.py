#!/usr/bin/env python
"""Try T63 (experiments/harness/TRIES.md): the ruled recomputed-groups file for T1's last upload, read at the test's
fraction of the outside atlas's bracket (2/5) on row T62's five practice tasks, with a placebo direction and the
epicardium-held arm at that size.

The ruled file is `build_test.py t1 --g-mom --base lib3 --feng-groups recomputed` (decision log 2026-11-19 (ii)). Row
T62 read its practice analogue `m1rc_gmom` at 1/5 of the atlas's E9.5 -> E14.5 bracket, the one-day practice fraction;
the builder applies 2/5 at the test. This script builds the same analogue at 2/5 on T62's tasks, draws, seeds,
protocol and scoring. It imports recomputed_read.py (T62) and, through it, bets_read.py (T61) and gmom_combo.py (T26),
and copies none of their arithmetic. It does not import perm_read.py, eb_read.py, step_read.py or any other script that
forms the atlas direction on another route (`source_guard`, tested); the placebo is written here and tested byte-equal
to perm_read.placebo_directions on synthetic directions.

What 2/5 on a one-day task is. The step's size is the bracket fraction, the time gap over the bracket's five days. At
the test, 2/5 spans two days; on a one-day practice task it is twice the step per elapsed day. So the doubled arm
matches the test's absolute step, not its step per day, and the ratio of the test's gap to every task's is 2 for every
arm here (J6). A pass shows that the test's absolute step does no harm on one-day truths at the released depth ratio;
it is not a check of the two-day horizon and says nothing about how deep the hidden truth is.

Tasks: row T62's five (the two E9.5 -> E10.5 cuts of configs/pseudo_tasks.yaml, chamber primary, start E8.5; row T60's
three library tasks as T61 and T62 ran them). Every task scores one truth, E10.5. Seeds 0 to 9.

Arms per seed, all from one cell draw (the draw none_gmom's post_cells hook receives):
  none_gmom              gate, not counted: T62's none_gmom (`bets_read.file_hook`), rebuilt
  m1rc_gmom              gate, not counted: T62's candidate at 1/5, rebuilt the same way
  m1rc_gmom_x2           counted: m1rc_gmom with the direction at 2/5; nothing else changes
  m1rc_gmom_x2_epi_held  counted: m1rc_gmom_x2 with the Epicardium-mapped cells given back their drawn counts after
                         the direction step (`recomputed_read.hold_group`), as T62's epi-held arm at 1/5
  m1rc_gmom_x2_epi_eb    counted, reported only: m1rc_gmom_x2 with the Epicardium-mapped cells taken from the same draw
                         moved by the Epicardium direction shrunk by `feng_dirs.eb_shrink` (the recomputed statistics'
                         own halves, ten expression bins), on the same stream [seed, 2]; every other cell is
                         m1rc_gmom_x2's byte for byte
  placebo_x1, placebo_x2 counted: the candidate on the placebo direction (offset 0) at 1/5 and 2/5
  placebo_x2_q1 .. _q4   counted: the candidate on placebo offsets 1 to 4 at 2/5
On each library task copy_last is scored first, as T61 and T62 scored it (a reference, not counted).

The placebo (row T45's: try T15's permuted control). Per recomputed group g in harness.LP_GROUPS, the finite values of
the recomputed direction permuted across genes by harness.t3_decile_permutation within harness.LP_BINS = 10 runs of
genes ranked by the group's mean CP10k over the two bracket ends (`feng_dirs.direction_halves` on the recomputed
statistics with the recomputed group map, 'mean_cp'), seed [harness.LP_PERM_SEED + offset, LP_GROUPS.index(g)]; NaN
stays NaN. A group outside LP_GROUPS or without its own halves is refused. Computed once per task, the same for every
seed. It keeps each expression bin's values (the between-bin part of the direction is kept, printed) and moves values
across held and non-held genes alike (R2's held genes then drop some; the applied size on the non-held genes is
printed beside the candidate's, flagged beyond 10%). The g_mom step and its mask are the candidate's.

  run-cut  One E10.5 cut. Refused before anything loads unless row T63 is in TRIES.md, the output file does not exist,
           the tree is clean by vec.provenance's rule, this script's and recomputed_read.py's source guards pass and
           the relabelling's committed checks pass (`recomputed_read.relabel_checks`). Prints no score.
  run-lib  One library task; the same refusals. Prints no score.
  read     T63's statements, applied before anything is printed; writes t63.json and t63.txt; prints the outcome words
           and what follows only.

Usage: on the DGX campaign checkout ~/v/vec-campaign-d (it alone holds data/derived/feng_relabel_c/), repo root, at the
row's commit, committing nothing until the five task files exist:
  export NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
  uv run --frozen --no-sync python -I -B experiments/t1_m1/double_step_read.py run-cut --out t63_double_step \
      --cut chamber
  (and --cut whole, run-lib --library 1, 2, 3), then, at the commit that holds the five task files,
  uv run --frozen --no-sync python -I -B experiments/t1_m1/double_step_read.py read --out t63_double_step
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
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
import recomputed_read as rc  # noqa: E402  (T62: its arms' pieces, guards and checks)

br, gc = rc.br, rc.gc
bt, gv, pt, rr = rc.bt, rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs

ROW = "T63"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"
T62_RELIABILITY = T62_TASKS / "reliability.json"
OUT_FILES = ("t63.json", "t63.txt")

CONTROL, X1 = rc.CONTROL, rc.CANDIDATE  # none_gmom, m1rc_gmom: the gates
X2 = "m1rc_gmom_x2"
X2_EPI_HELD = "m1rc_gmom_x2_epi_held"
X2_EPI_EB = "m1rc_gmom_x2_epi_eb"
P1, P2 = "placebo_x1", "placebo_x2"
QS = tuple(f"placebo_x2_q{i}" for i in range(1, 5))
PLACEBO_OFFSET = {P1: 0, P2: 0, **{q: i for i, q in enumerate(QS, start=1)}}
X2_PLACEBOS = (P2, *QS)  # S4 at 2/5 is decided over all five permutations
X1_PLACEBOS = (P1,)  # S4 at 1/5: one permutation, reported only
GATED = (CONTROL, X1)
ARMS = (CONTROL, X1, X2, X2_EPI_HELD, X2_EPI_EB, P1, *X2_PLACEBOS)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first
LIB_ARMS = (LIB_REFERENCE, *ARMS)
COUNTED = (X2, X2_EPI_HELD, X2_EPI_EB, P1, *X2_PLACEBOS)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 45
FIRST_TRY = 368  # T62 ended at T1 368
TRY_RANGE = (FIRST_TRY, FIRST_TRY + TRIES_COUNT)  # T1 368 to 413

PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2: one day of the E9.5 -> E14.5 bracket
TEST_FRACTION = 0.4  # the builder's at the test: E10.5 -> E12.5
FRACTION_TOL = 1e-12
FRACTIONS = {X1: PRACTICE_FRACTION, X2: TEST_FRACTION, X2_EPI_HELD: TEST_FRACTION, X2_EPI_EB: TEST_FRACTION,
             P1: PRACTICE_FRACTION, **{q: TEST_FRACTION for q in X2_PLACEBOS}}
MULTIPLIER = {X1: 1.0, X2: 2.0, X2_EPI_HELD: 2.0, X2_EPI_EB: 2.0, P1: 1.0, **{q: 2.0 for q in X2_PLACEBOS}}
EPI_GROUP, EPI_TASK = rc.EPI_GROUP, rc.EPI_TASK  # "Epicardium", library 3
REPRO_TOL = rc.REPRO_TOL  # 1e-9
BAR = 0.3
MIN_SEEDS = 8
WHOLE_FLOOR = -0.3
LIB_FLOOR = -0.3
MIN_LIBRARIES = 2
SIZE_FLAG = 0.10  # a placebo whose applied size on the non-held genes is this far from the candidate's is named
BINS = harness.LP_BINS
GROUP_WEIGHTS = rc.GROUP_WEIGHTS
BANNED_IMPORTS = ("perm_read", "eb_read", "step_read", "deep_read", "over_read", "lib_step_read", "ctl_range_read",
                  "n48_m1_direction_shrink", "n12_hygiene_size", "n54" + "_relabel_listing",
                  "n54" + "_authors_agreement")  # spelled in pieces so that recomputed_read's guard passes this file

S1_WORDS = {True: "the doubled step is AHEAD of the no-outside-data analogue",
            False: "the doubled step is NOT AHEAD of the no-outside-data analogue"}
S1_QUALIFIER = ("the test's absolute step on one-day truths at the released depth ratio; not a check of the two-day "
                "horizon and nothing about the hidden truth's depth; the record predicted a pass")
S2_WORDS = ("GAINS", "LOSES", "level")
S2_QUALIFIER = ("decides nothing and is never evidence for any fraction or step size; the record predicted GAINS")
S3_WORDS = ("PAYS", "COSTS", "within noise", "unread")
NEEDS = "the gain needs the direction's own genes"
MORE = "a placebo gains as well; the direction more"
PLACEBO_GAINS = "a placebo gains as well"
NOT_SEPARATED = "not separated"
S4_WORDS = (NEEDS, MORE, PLACEBO_GAINS, NOT_SEPARATED)

BUILD = "build"
ASK = "ask"
VOID = "void"
RULED = "the ruled recomputed file (`t1 --g-mom --base lib3 --feng-groups recomputed`)"
THIRD = ("a third file, the recomputed groups with the epicardium group held (it needs an off-by-default builder "
         "option, a check build and his ruling)")
NONE_FILE = "the no-outside-data file (`t1 --g-mom --base lib3 --feng-groups none`, row R8's, built)"

LABEL = ("try T63, a labelled in-sample diagnostic, never an adoption result: T1's ruled recomputed-groups file at the "
         "test's fraction of the atlas bracket (2/5) against its single step and the no-outside-data analogue, with a "
         "placebo direction (five within-bin permutations at 2/5, one at 1/5) and the epicardium step at 2/5, on row "
         "T62's five tasks; every task scores one truth (E10.5) over one day, so 2/5 is twice the test's step per day")


# ---------------------------------------------------------------------------------------------------------------------
# Guards on this source


def source_guard(path: Path | None = None) -> list[str]:
    """recomputed_read's guard on this file (every feng_directions call names the recomputed groups; no token of the
    other groups, their statistics or row N54's listing), and no import of a script that forms the atlas direction on
    another route (BANNED_IMPORTS) or any dynamic import. Returns the failures."""
    path = Path(__file__ if path is None else path)
    bad = rc.source_guard(path)
    text = path.read_text(encoding="utf-8")
    for m in re.finditer(r"^\s*(?:from\s+(\S+)\s+import\s+([^\n#]+)|import\s+([^\n#]+))", text, re.M):
        names = re.split(r"[\s,.()]+", " ".join(x for x in m.groups() if x))
        hit = sorted(set(names) & set(BANNED_IMPORTS))
        if hit:
            bad.append(f"the source imports {hit}")
    for token in ("import" + "lib", "__imp" + "ort__", "spec_from" + "_file_location"):
        if token in text:
            bad.append(f"the source names {token}")
    return bad


# ---------------------------------------------------------------------------------------------------------------------
# The placebo and the shrunk epicardium direction (target-free)


def recomputed_halves(genes: Sequence[str], directions: Mapping) -> dict:
    """The recomputed statistics' own halves for every group with a direction (`feng_dirs.direction_halves`, the
    recomputed group map, the bracket E9.5 -> E14.5), through feng_dirs' checks. Refused unless every group with a
    direction has halves of its own."""
    path = feng_dirs.recomputed_stats_path()
    feng_dirs.check_recomputed(path)
    st = feng_dirs.dl.load_stats(path)
    feng_dirs.check_recomputed_keys(st)
    cfg = feng_dirs.groups_config(feng_dirs.load_config(), feng_dirs.load_relabel_config())
    halves = feng_dirs.direction_halves(st, list(genes), cfg, feng_dirs.BRACKET, groups=set(directions))
    missing = sorted(set(directions) - set(halves))
    if missing:
        raise ValueError(f"groups without halves of their own: {missing}: no placebo or shrinkage can be formed")
    return halves


def placebo_directions(directions: Mapping[str, np.ndarray], halves: Mapping[str, Mapping[str, np.ndarray]],
                       offset: int = 0) -> dict[str, np.ndarray]:
    """Row T45's placebo (try T15's permuted control) on the recomputed directions: per group, the finite values
    permuted across genes within LP_BINS bins of the group's mean CP10k on the bracket, seed [LP_PERM_SEED + offset,
    LP_GROUPS.index(group)]. NaN stays NaN. A group outside LP_GROUPS or without halves is refused."""
    out = {}
    for g, d in directions.items():
        if g not in harness.LP_GROUPS:
            raise ValueError(f"{g}: not one of try T15's groups {list(harness.LP_GROUPS)}: it has no permutation seed")
        if g not in halves:
            raise ValueError(f"{g}: no halves give a mean CP10k for this direction: it cannot be permuted")
        d = np.asarray(d, dtype=np.float64).copy()
        ok = np.isfinite(d)
        d[ok] = harness.t3_decile_permutation(d[ok], np.asarray(halves[g]["mean_cp"], dtype=np.float64)[ok],
                                              [harness.LP_PERM_SEED + offset, harness.LP_GROUPS.index(g)],
                                              harness.LP_BINS)
        out[g] = d
    return out


def placebo_set(directions: Mapping, halves: Mapping) -> dict[str, dict[str, np.ndarray]]:
    """arm -> its placebo directions (placebo_x1 and placebo_x2 share offset 0)."""
    cache = {o: placebo_directions(directions, halves, o) for o in sorted(set(PLACEBO_OFFSET.values()))}
    return {arm: cache[o] for arm, o in PLACEBO_OFFSET.items()}


def epi_shrunk(directions: Mapping, halves: Mapping) -> tuple[dict, list[float] | None]:
    """The directions with the Epicardium group's shrunk by feng_dirs.eb_shrink (its halves, LP_BINS bins of mean
    CP10k); the other groups unchanged (the same arrays). Without an Epicardium direction: unchanged, None."""
    if EPI_GROUP not in directions:
        return dict(directions), None
    h = halves[EPI_GROUP]
    shrunk, factors = feng_dirs.eb_shrink(directions[EPI_GROUP], h["d_a"], h["d_b"], h["mean_cp"], BINS)
    return {**directions, EPI_GROUP: shrunk}, [float(f) for f in factors]


def direction_sha(directions: Mapping[str, np.ndarray]) -> dict[str, str]:
    return {g: hashlib.sha256(np.ascontiguousarray(d, dtype=np.float64).tobytes()).hexdigest()
            for g, d in sorted(directions.items())}


def free_genes(held: Sequence[str], genes: Sequence[str]) -> np.ndarray:
    return ~np.isin(np.asarray(list(genes), dtype=object), list(held))


def applied_size(directions: Mapping[str, np.ndarray], held: Sequence[str], genes: Sequence[str], fraction: float
                 ) -> dict[str, dict[str, float]]:
    """Per group, the log factors at `fraction` (fraction x d, NaN as 0) on the non-held genes: mean absolute value
    and root mean square. Target-free."""
    free = free_genes(held, genes)
    out = {}
    for g, d in sorted(directions.items()):
        logf = float(fraction) * np.nan_to_num(np.asarray(d, dtype=np.float64), nan=0.0)[free]
        out[g] = {"mean_abs_logf": float(np.abs(logf).mean()) if logf.size else 0.0,
                  "rms_logf": float(np.sqrt((logf * logf).mean())) if logf.size else 0.0}
    return out


def size_ratios(placebo: Mapping, candidate: Mapping, held: Sequence[str], genes: Sequence[str]) -> dict:
    """Per group, the placebo's applied size over the candidate's (mean |log f| and RMS on the non-held genes; the
    fraction cancels) and whether either is more than SIZE_FLAG from 1."""
    p, c = applied_size(placebo, held, genes, 1.0), applied_size(candidate, held, genes, 1.0)
    out = {}
    for g in sorted(candidate):
        r = {k: (p[g][k] / c[g][k] if c[g][k] > 0 else float("nan")) for k in ("mean_abs_logf", "rms_logf")}
        r["flagged"] = bool(any(not abs(v - 1.0) <= SIZE_FLAG for v in r.values()))
        out[g] = r
    return out


def between_bin_share(directions: Mapping[str, np.ndarray], halves: Mapping[str, Mapping[str, np.ndarray]],
                      free: np.ndarray | None = None) -> dict[str, float]:
    """Per group, the share of the direction's variance between the placebo's bins (each finite gene given its bin's
    mean; with `free`, the bins stay over all finite genes and means and variances are over the free ones): the part
    the placebo keeps. NaN when d has no variance."""
    out = {}
    for g, d in sorted(directions.items()):
        d = np.asarray(d, dtype=np.float64)
        ok = np.isfinite(d)
        v, basis = d[ok], np.asarray(halves[g]["mean_cp"], dtype=np.float64)[ok]
        keep = np.ones(v.size, dtype=bool) if free is None else np.asarray(free, dtype=bool)[ok]
        means = np.empty_like(v)
        for genes in np.array_split(np.argsort(basis, kind="stable"), BINS):
            genes = genes[keep[genes]]
            if genes.size:
                means[genes] = v[genes].mean()
        total = float(v[keep].var()) if keep.any() else 0.0
        out[g] = float(means[keep].var()) / total if total > 0 else float("nan")
    return out


def direction_record(directions: Mapping, placebos: Mapping, shrunk: Mapping, eb_factors, halves: Mapping,
                     held_m1: Sequence[str], genes: Sequence[str]) -> dict:
    """What the task file keeps of the directions: the sha256 of each group's candidate, placebo and shrunk direction,
    the placebos' applied sizes over the candidate's, the between-bin shares and the Epicardium shrinkage factors."""
    free = free_genes(held_m1, genes)
    return {"candidate_sha256": direction_sha(directions),
            "placebo_sha256": {a: direction_sha(p) for a, p in placebos.items()},
            "placebo_offsets": dict(PLACEBO_OFFSET),
            "epi_shrunk_sha256": direction_sha(shrunk),
            "epi_eb_factors": eb_factors,
            "epi_eb_factor_mean": None if eb_factors is None else float(np.mean(eb_factors)),
            "size_ratios": {a: size_ratios(placebos[a], directions, held_m1, genes) for a in placebos},
            "between_bin_share_free": between_bin_share(directions, halves, free)}


def prediction_sha(stage) -> str:
    """The sha256 of an arm's scored matrix (cells x board genes, float32, C order)."""
    return hashlib.sha256(np.ascontiguousarray(dense(stage.X), dtype=np.float32).tobytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def cut_arms(ctx: harness.BuildContext, factors, directions: Mapping, placebos: Mapping, shrunk: Mapping,
             fraction: float, held_m1: Sequence[str]) -> dict[str, Stage]:
    """Every arm of one seed on a harness cut, from one draw (none_gmom's post_cells input), each finished as the
    harness finishes a candidate. none_gmom and m1rc_gmom are recomputed_read.cut_pair's (tested byte-equal); the
    doubled arms run the same build at 2 x fraction; the placebo arms the same build on the placebo directions, with
    the candidate's g_mom hook and mask."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, directions)
    _, covered_h, with_direction_h = rc.direction_kit(labels, rc.without_epicardium(directions))
    control_hook, drawn = br.file_hook(factors, ctx.seed, held), []
    hook = rc.rc_hook(factors, ctx.seed, held, with_direction, covered)
    epi_hook = rc.rc_hook(factors, ctx.seed, held, with_direction_h, covered_h)
    f2 = 2.0 * float(fraction)

    def control(pop):
        drawn.append(pop)  # the draw itself: no direction and no momentum come before post_cells on this arm
        return control_hook(pop)

    def build(dirs=None, fr=None, post=None) -> Stage:
        kw = ({"held_flat": held} if dirs is None else
              {"held_flat": list(held_m1), "directions": dirs, "direction_fraction": fr, "direction_groups": groups})
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src, post_cells=post, **kw)
        return harness.finish(ctx, ref, dense(sub.X), None)

    def eb_cells(moved):
        other = gc.moved_by_feng(drawn[0], shrunk, f2, groups, list(held_m1), [ctx.seed, 2])
        return hook(rc.hold_group(moved, other, groups))

    arms = {CONTROL: build(post=control)}  # first: it hands the draw to the arms that splice cells
    arms[X1] = build(directions, fraction, hook)
    arms[X2] = build(directions, f2, hook)
    arms[X2_EPI_HELD] = build(directions, f2, lambda moved: epi_hook(rc.hold_group(moved, drawn[0], groups)))
    arms[X2_EPI_EB] = build(directions, f2, eb_cells)
    for arm in (P1, *X2_PLACEBOS):
        arms[arm] = build(placebos[arm], MULTIPLIER[arm] * float(fraction), hook)
    return arms


def library_arms(e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping, placebos: Mapping,
                 shrunk: Mapping, fraction: float, n_cells: int = rr.N_CELLS) -> dict:
    """Every arm of one seed on a library task, from T60's oracle draw (recomputed_read.library_pair's route:
    `gmom_combo.moved_by_feng` on stream [seed, 2], then g_mom), as AnnData."""
    base, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    labels = e95.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, directions)
    _, covered_h, with_direction_h = rc.direction_kit(labels, rc.without_epicardium(directions))
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    epi_hook = rc.rc_hook(factors, seed, held, with_direction_h, covered_h)
    f2 = 2.0 * float(fraction)

    def moved(dirs, fr):
        return gc.moved_by_feng(base, dirs, fr, groups, list(held_m1), [seed, 2])

    x2 = moved(directions, f2)
    arms = {CONTROL: br.file_hook(factors, seed, held)(base), X1: hook(moved(directions, fraction)), X2: hook(x2),
            X2_EPI_HELD: epi_hook(rc.hold_group(x2, base, groups)),
            X2_EPI_EB: hook(rc.hold_group(x2, moved(shrunk, f2), groups))}
    for arm in (P1, *X2_PLACEBOS):
        arms[arm] = hook(moved(placebos[arm], MULTIPLIER[arm] * float(fraction)))
    return arms


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def refuse(out: Path, name: str) -> tuple[dict, dict]:
    """The refusals before anything loads: the row, an existing output file, an unclean tree, both source guards and
    the relabelling's committed checks. Returns the start state and the relabelling's checks."""
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
    return started, checks


def arm_fractions(fraction: float) -> dict[str, float]:
    """Every direction arm's applied fraction (the task file records it; the read checks it against FRACTIONS)."""
    return {arm: m * float(fraction) for arm, m in MULTIPLIER.items()}


def check_fractions(fraction: float, where: str) -> None:
    if abs(fraction - PRACTICE_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")
    test = feng_dirs.bracket_fraction(10.5, 12.5)
    if abs(test - TEST_FRACTION) > FRACTION_TOL or abs(2.0 * fraction - test) > FRACTION_TOL:
        sys.exit(f"{where}: the test's fraction is {test}, not {TEST_FRACTION} = twice the practice fraction")


def _inputs_for(genes: Sequence[str], held_m1: Sequence[str]) -> tuple[dict, dict, dict, dict, list | None, dict]:
    directions, ext = rc.recomputed_inputs(list(genes))
    halves = recomputed_halves(genes, directions)
    placebos = placebo_set(directions, halves)
    shrunk, eb = epi_shrunk(directions, halves)
    record = direction_record(directions, placebos, shrunk, eb, halves, held_m1, genes)
    return directions, ext, placebos, shrunk, eb, record


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
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    check_fractions(fraction, t.key)
    held_m1, union = rc._held_m1(preregistered)
    dir_genes = list(inputs[ref_role].var_names)
    directions, ext, placebos, shrunk, _, record = _inputs_for(dir_genes, held_m1)
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms = cut_arms(ctx, factors, directions, placebos, shrunk, fraction, held_m1)
        shas = {k: prediction_sha(v) for k, v in arms.items()}
        rows.append({"seed": s, "floor": float(f.score), "arms": gv.score_arms(board, target, ref_stage, s, f.anchors,
                                                                               arms),
                     "sha256": shas, "seconds": round(time.time() - t0, 1)})
        del arms
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "fractions": arm_fractions(fraction), "relabel": rc._relabel_fields(ext, checks), "directions": record,
           "constants": {"edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1),
                         "bins": BINS, "perm_seed": harness.LP_PERM_SEED},
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
    fraction = feng_dirs.bracket_fraction(9.5, 10.5)
    check_fractions(fraction, f"library {a.library}")
    held_m1, union = rc._held_m1(held)
    directions, ext, placebos, shrunk, _, record = _inputs_for(list(e95.var_names), held_m1)
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms = library_arms(e95, shares, factors, held, held_m1, s, within, fb, directions, placebos, shrunk,
                            fraction)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        shas = {k: prediction_sha(v) for k, v in stages.items()}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s), "sha256": shas,
                     "seconds": round(time.time() - t0, 1)})
        del arms, stages
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "fractions": arm_fractions(fraction), "relabel": rc._relabel_fields(ext, checks), "directions": record,
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD,
                         "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1),
                         "bins": BINS, "perm_seed": harness.LP_PERM_SEED},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], feng_dirs.RELABEL_CONFIG,
                                       rc.RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def task_files() -> dict[str, str]:
    return rc.task_files()


def committed_t62(name: str, tasks: Path | None = None) -> dict[str, dict[int, float]]:
    """arm -> seed -> score of the gated arms in row T62's committed task file."""
    rows = json.loads(((T62_TASKS if tasks is None else tasks) / name).read_text())["rows"]
    return {arm: {int(r["seed"]): float(r["arms"][arm]["score"]) for r in rows} for arm in GATED}


def diff(rows: Sequence[Mapping], plus: Sequence[str], minus: Sequence[str]) -> dict:
    """Per seed, the sum of `plus` arms' scores less the sum of `minus` arms' (board points): mean, sd, seeds above and
    below 0, range, and each scored group's points (100 x weight x the mean group-skill difference)."""
    def val(r, key):
        return sum(r["arms"][a][key] for a in plus) - sum(r["arms"][a][key] for a in minus)

    d = np.array([val(r, "score") for r in rows], dtype=np.float64)
    groups = {g: float(100.0 * w * np.mean([sum(r["arms"][a]["groups"][g] for a in plus)
                                             - sum(r["arms"][a]["groups"][g] for a in minus) for r in rows]))
              for g, w in GROUP_WEIGHTS.items()}
    return {"plus": list(plus), "minus": list(minus), "mean": float(d.mean()),
            "sd": float(d.std(ddof=1)) if d.size > 1 else 0.0, "ahead": int((d > 0).sum()),
            "behind": int((d < 0).sum()), "of": int(d.size), "range": [float(d.min()), float(d.max())],
            "groups": groups}


def ahead(p: Mapping) -> bool:
    """The single-task test, ahead: a mean above +0.3 with at least 8 of 10 seeds above 0."""
    return bool(p["mean"] > BAR and p["ahead"] >= MIN_SEEDS)


def behind(p: Mapping) -> bool:
    return bool(p["mean"] < -BAR and p["behind"] >= MIN_SEEDS)


def passes(chamber: Mapping, whole: Mapping) -> bool:
    """The bar of rows T61 and T62: ahead on the chamber cut and the whole cut's mean at least -0.3."""
    return bool(ahead(chamber) and whole["mean"] >= WHOLE_FLOOR)


def library_valid(r: Mapping) -> bool:
    return not any(x.get("inverted_anchors") for x in r["rows"])


def guards(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[str, Mapping[int, float]]],
           checks: Mapping, source_bad: Sequence[str], seeds: Sequence[int] = SEEDS) -> list[str]:
    """Every void condition of row T63; an empty list means none fired."""
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
        for arm in GATED:
            ref = committed.get(k, {}).get(arm, {})
            worst = max((abs(br.score_of(x, arm) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                         for x in rows), default=float("inf"))
            if not worst <= REPRO_TOL:
                reasons.append(f"{k}: {arm} differs from row T62's committed {arm} by {worst:.3g} (tolerance "
                               f"{REPRO_TOL:g})")
        got = r.get("fractions") or {}
        if not abs(float(r.get("fraction", float("nan"))) - PRACTICE_FRACTION) <= FRACTION_TOL:
            reasons.append(f"{k}: the practice fraction is {r.get('fraction')}, not {PRACTICE_FRACTION}")
        for arm, want_f in FRACTIONS.items():
            if not abs(float(got.get(arm, float("nan"))) - want_f) <= FRACTION_TOL:
                reasons.append(f"{k}: {arm}'s fraction is {got.get(arm)}, not {want_f}")
        rel = r.get("relabel") or {}
        if rel.get("sha256") != rc.STATS_SHA256 or rel.get("stats_sha256_hashed") != rc.STATS_SHA256:
            reasons.append(f"{k}: ran on another recomputed statistics file")
        if rel.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: ran under another settings file")
        if rel.get("route") != rc.GROUPS_OPTION or rel.get("record_route") != rc.GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        if rel.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the record did not say authors_labels_read false")
        reasons += [f"{k}: {x}" for x in placebo_reasons(r.get("directions") or {})]
        for change, base in ((X2, X1), (P1, X1), (P2, X2)):
            if all(x["sha256"][change] == x["sha256"][base] for x in rows):
                reasons.append(f"{k}: {change} built the same prediction as {base} on every seed")
    if len(commits) > 1:
        reasons.append(f"the task files ran at more than one commit: {sorted(map(str, commits))}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def placebo_reasons(d: Mapping) -> list[str]:
    """The placebo's record: every placebo arm's directions differ from the candidate's in every group; placebo_x1 and
    placebo_x2 are one permutation (offset 0); the five permutations differ from each other in every group."""
    cand, pl = d.get("candidate_sha256") or {}, d.get("placebo_sha256") or {}
    if not cand or set(pl) != set(PLACEBO_OFFSET):
        return ["the placebo record is missing"]
    out = []
    for arm, shas in pl.items():
        if set(shas) != set(cand):
            out.append(f"{arm}: its groups {sorted(shas)} are not the candidate's {sorted(cand)}")
        elif any(shas[g] == cand[g] for g in cand):
            out.append(f"{arm}: a group's placebo is the candidate's direction")
    if pl.get(P1) != pl.get(P2):
        out.append("placebo_x1 and placebo_x2 are not one permutation")
    for a, b in itertools.combinations(X2_PLACEBOS, 2):
        if any(pl.get(a, {}).get(g) == pl.get(b, {}).get(g) for g in cand):
            out.append(f"{a} and {b} share a permutation in some group")
    return out


def s1_statement(tasks: Mapping[str, Mapping], valid: Sequence[str]) -> dict:
    c, w = tasks["chamber"], tasks["whole"]
    cuts = passes(c, w)
    low = [lib for lib in valid if tasks[lib]["mean"] < LIB_FLOOR]
    libraries = bool(len(valid) >= MIN_LIBRARIES and not low)
    holds = bool(cuts and libraries)
    return {"holds": holds, "words": S1_WORDS[holds], "cuts_ahead": cuts, "library_check": libraries,
            "libraries_below_floor": low}


def s2_statement(p: Mapping) -> str:
    return "GAINS" if ahead(p) else "LOSES" if behind(p) else "level"


def s3_statement(p: Mapping | None, valid: bool, directed: bool) -> str:
    if not (valid and directed) or p is None:
        return "unread"
    return "PAYS" if ahead(p) else "COSTS" if behind(p) else "within noise"


def s4_category(beats: Sequence[bool], placebo_passes: Sequence[bool]) -> str:
    """Rows T45 and T58's categories over the permutations at one size: the direction beats the placebo only if it
    beats every permutation; a placebo gains if any permutation passes against none_gmom."""
    beat, gains = bool(beats) and all(beats), any(placebo_passes)
    if beat and not gains:
        return NEEDS
    if beat and gains:
        return MORE
    if gains:
        return PLACEBO_GAINS
    return NOT_SEPARATED


def s4_statement(loaded: Mapping[str, Mapping], candidate: str, placebos: Sequence[str]) -> dict:
    rows = {k: loaded[k]["rows"] for k in ("chamber", "whole")}
    per = {}
    for p in placebos:
        cand = {k: diff(v, [candidate], [p]) for k, v in rows.items()}
        plac = {k: diff(v, [p], [CONTROL]) for k, v in rows.items()}
        per[p] = {"candidate_minus_placebo": cand, "placebo_minus_none": plac,
                  "beats": passes(cand["chamber"], cand["whole"]), "placebo_passes": passes(plac["chamber"],
                                                                                            plac["whole"])}
    word = s4_category([v["beats"] for v in per.values()], [v["placebo_passes"] for v in per.values()])
    worst = min(per, key=lambda p: per[p]["candidate_minus_placebo"]["chamber"]["mean"])
    return {"words": word, "permutations": list(placebos), "per_permutation": per, "worst": worst}


def follows(s1_holds: bool | None, s3: str | None, s4_x2: str | None, void: bool = False) -> dict:
    """The one 'what follows' line for a combination of outcomes (fixed in row T63 before the run). The captain's
    ruling (the recomputed file) stands on every branch unless he changes it; nothing is built or changed by the
    read."""
    if void:
        return {"key": VOID, "files": [], "line": (
            "What follows: the check is VOID; no rerun under row T63 and its tries count. One message to the captain: "
            "the check is void and his ruling (the recomputed file) stands; the ruled file is not built until his "
            "word.")}
    if s1_holds not in (True, False) or s3 not in S3_WORDS or s4_x2 not in S4_WORDS:
        raise ValueError(f"not an outcome: {s1_holds!r}, {s3!r}, {s4_x2!r}")
    if s1_holds and s3 != "COSTS" and s4_x2 in (NEEDS, MORE):
        line = (f"What follows: {RULED} is built as ruled; the result goes to the captain beside it; his ruling "
                f"stands unless he changes it.")
        if s4_x2 == MORE:
            line += " The message says that a part of the gain comes with any direction of this size profile."
        return {"key": BUILD, "files": [RULED], "line": line}
    files = [RULED]
    if s3 == "COSTS":
        files.append(THIRD)
    if not s1_holds or s4_x2 in (PLACEBO_GAINS, NOT_SEPARATED):
        files.append(NONE_FILE)
    key = ASK + ("+third" if THIRD in files else "") + ("+none" if NONE_FILE in files else "")
    return {"key": key, "files": files, "line": (
        "What follows: one message to the captain before any build, with the result, naming: " + "; ".join(files)
        + ". His ruling (the recomputed file) stands unless he changes it; nothing is built until he answers.")}


FOLLOW_KEYS = (BUILD, ASK + "+third", ASK + "+none", ASK + "+third+none", VOID)


def reliability_beside(path: Path | None = None) -> dict | None:
    """Row T62's committed, target-free reliability figures for the Epicardium group (beside S3)."""
    path = T62_RELIABILITY if path is None else Path(path)
    if not path.exists():
        return None
    g = json.loads(path.read_text())["groups"].get(EPI_GROUP) or {}
    if not g.get("own_direction"):
        return {"own_direction": False}
    return {"own_direction": True, "cells": g["cells"], "spearman_all": g["spearman_all"],
            "eb_factor_mean": g["eb_factor_mean"], "decile_spearman": [x["spearman"] for x in g["deciles"]]}


def t63_reading(loaded: Mapping[str, Mapping | None], committed: Mapping, checks: Mapping,
                source_bad: Sequence[str] = ()) -> dict:
    """Row T63's statements and what follows, applied to the task files."""
    reasons = guards(loaded, committed, checks, source_bad)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons, "follows": follows(None, None, None, void=True)}
    rows = {k: r["rows"] for k, r in loaded.items()}
    valid = [lib for lib in LIBRARIES if library_valid(loaded[lib])]
    vs_none = {k: diff(v, [X2], [CONTROL]) for k, v in rows.items()}
    s1 = s1_statement(vs_none, valid)
    vs_single = {k: diff(v, [X2], [X1]) for k, v in rows.items()}
    s2 = s2_statement(vs_single["chamber"])
    epi = {k: diff(v, [X2], [X2_EPI_HELD]) for k, v in rows.items()}
    directed = (loaded[EPI_TASK].get("relabel") or {}).get("epicardium_directed") is True
    s3 = s3_statement(epi[EPI_TASK], EPI_TASK in valid, directed)
    s4_x2 = s4_statement(loaded, X2, X2_PLACEBOS)
    s4_x1 = s4_statement(loaded, X1, X1_PLACEBOS)
    reported = {
        "x1_minus_none": {k: diff(v, [X1], [CONTROL]) for k, v in rows.items()},
        "eb_minus_x2": {k: diff(v, [X2_EPI_EB], [X2]) for k, v in rows.items()},
        "eb_minus_epi_held": {k: diff(v, [X2_EPI_EB], [X2_EPI_HELD]) for k, v in rows.items()},
        "placebo_minus_none": {p: {k: diff(v, [p], [CONTROL]) for k, v in rows.items()} for p in (P1, *X2_PLACEBOS)},
        "candidate_minus_placebo": {p: {k: diff(v, [X2 if p in X2_PLACEBOS else X1], [p]) for k, v in rows.items()}
                                    for p in (P1, *X2_PLACEBOS)},
        "increment_vs_placebo_increment": {k: diff(v, [X2, P1], [X1, P2]) for k, v in rows.items()},
        "mean_scores": {k: {a: float(np.mean([br.score_of(x, a) for x in v]))
                            for a in (LIB_ARMS if k in LIBRARIES else ARMS)} for k, v in rows.items()}}
    sizes = {k: (r.get("directions") or {}) for k, r in loaded.items()}
    flags = sorted({f"{k}: {a} {g}" for k, d in sizes.items() for a, gs in (d.get("size_ratios") or {}).items()
                    for g, x in gs.items() if x.get("flagged")})
    f = follows(s1["holds"], s3, s4_x2["words"])
    return {"outcome": "read", "S1": s1, "S1_tasks": vs_none, "S2": s2, "S2_tasks": vs_single, "S3": s3,
            "S3_tasks": epi, "S3_task_valid": EPI_TASK in valid, "S3_epicardium_directed": directed,
            "S4_x2": s4_x2, "S4_x1_reported": s4_x1, "valid_libraries": valid, "reported": reported,
            "j6_range": sorted([0.0, vs_none["chamber"]["mean"]]), "size_flags": flags,
            "directions": sizes, "follows": f}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]["line"]]
    s4 = r["S4_x2"]
    return [f"S1: \"{r['S1']['words']}\" ({S1_QUALIFIER}; cuts: {'AHEAD' if r['S1']['cuts_ahead'] else 'NOT AHEAD'}; "
            f"library check: {'passed' if r['S1']['library_check'] else 'failed'})",
            f"S2 (m1rc_gmom_x2 - m1rc_gmom on the chamber cut): {r['S2']} ({S2_QUALIFIER})",
            f"S3 (the epicardium step at 2/5, m1rc_gmom_x2 - m1rc_gmom_x2_epi_held on library {EPI_TASK}; a valid "
            f"library: {'yes' if r['S3_task_valid'] else 'no'}): {r['S3']}",
            f"S4 at 2/5: \"{s4['words']}\" (decided over {len(s4['permutations'])} within-bin permutations, offsets "
            f"0 to {len(s4['permutations']) - 1}; the worst for the direction: {s4['worst']}); at 1/5, reported only, "
            f"one permutation: \"{r['S4_x1_reported']['words']}\"",
            r["follows"]["line"]]


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} behind {p['behind']} of "
            f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}")


def t63_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        rep = r["reported"]
        blocks = [("S1: m1rc_gmom_x2 - none_gmom", r["S1_tasks"]), ("S2: m1rc_gmom_x2 - m1rc_gmom", r["S2_tasks"]),
                  ("S3: m1rc_gmom_x2 - m1rc_gmom_x2_epi_held", r["S3_tasks"]),
                  ("reported: m1rc_gmom - none_gmom (T62's pair, rebuilt)", rep["x1_minus_none"]),
                  ("reported beside S3: m1rc_gmom_x2_epi_eb - m1rc_gmom_x2 (the shrunk epicardium step against the "
                   "full one)", rep["eb_minus_x2"]),
                  ("reported beside S3: m1rc_gmom_x2_epi_eb - m1rc_gmom_x2_epi_held", rep["eb_minus_epi_held"]),
                  ("reported beside S4: (m1rc_gmom_x2 - m1rc_gmom) - (placebo_x2 - placebo_x1), offset 0 (row T45's "
                   "increments)", rep["increment_vs_placebo_increment"])]
        for p in (P1, *X2_PLACEBOS):
            blocks.append((f"S4: candidate - {p} (offset {PLACEBO_OFFSET[p]})", rep["candidate_minus_placebo"][p]))
            blocks.append((f"S4: {p} - none_gmom", rep["placebo_minus_none"][p]))
        for title, per_task in blocks:
            lines.append(title)
            lines += [_line(k, p) for k, p in per_task.items()]
        lines.append("Mean scores: " + "; ".join(f"{k}: " + ", ".join(f"{a} {v:.2f}" for a, v in m.items())
                                                 for k, m in rep["mean_scores"].items()))
        for k, d in r["directions"].items():
            sr_ = d.get("size_ratios") or {}
            sizes = "; ".join(f"{a} " + ", ".join(f"{g} {x['mean_abs_logf']:.3f}/{x['rms_logf']:.3f}"
                                                  for g, x in gs.items()) for a, gs in sr_.items())
            lines.append(f"{k}: placebo applied size over the candidate's (mean |log f|, RMS; non-held genes): {sizes}")
            lines.append(f"{k}: between-bin share of each direction's variance (kept by the placebo): " + ", ".join(
                f"{g} {v:.3f}" for g, v in (d.get("between_bin_share_free") or {}).items())
                + f"; Epicardium eb_shrink factor mean {d.get('epi_eb_factor_mean')}")
        lines.append(f"Placebo sizes more than {SIZE_FLAG:.0%} from the candidate's: "
                     f"{', '.join(r['size_flags']) or 'none'}")
        rel = res.get("epicardium_reliability_t62")
        if rel:
            lines.append(f"Beside S3, row T62's committed reliability listing (target-free) for Epicardium: {rel}")
        lines.append(f"Valid libraries: {', '.join(r['valid_libraries']) or 'none'}; S1's libraries below "
                     f"{LIB_FLOOR}: {', '.join(r['S1']['libraries_below_floor']) or 'none'}")
        lines.append(f"J6 range (0 to the chamber cut's m1rc_gmom_x2 - none_gmom; the ratio of the test's gap to every "
                     f"task's is 2, so the gain counts zero at the centre): {r['j6_range'][0]:+.2f} to "
                     f"{r['j6_range'][1]:+.2f}")
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
    reading = t63_reading(loaded, committed, checks, bad)
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT,
           "epicardium_reliability_t62": reliability_beside()}
    res["stamp"] = gv.stamp(started, [p for p in paths.values() if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / OUT_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / OUT_FILES[1]).write_text(t63_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T63")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T63")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T63's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out)
    else:
        {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)


if __name__ == "__main__":
    main()
