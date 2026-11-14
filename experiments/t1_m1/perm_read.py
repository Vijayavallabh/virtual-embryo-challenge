#!/usr/bin/env python
"""Try T45 (experiments/harness/TRIES.md): what a larger M1 step's gain on the E10.5 cuts is made of. A labelled read on
try T44's route (experiments/t1_m1/step_read.py) with three controls.

If try T44 reads "a larger step scores better on these cuts", three explanations are open: (1) the direction: on day
one the real change along Feng's direction is larger than a fifth of the bracket; (2) count noise: M1 applies its
factors in count space (vec.core.apply_factors: Binomial thinning below 1, a Poisson top-up above 1), so a larger step
redraws more of every covered cell's counts (try T43c: a count redraw alone gains 3.8 points on these cuts; T1's depth
lever is closed); (3) any set of factors with this size profile, whatever the genes. The controls separate them.

  run    One E10.5 cut, 10 seeds, on the harness's draws and scoring exactly as T44's run; refused unless the committed
         N22 passed, and unless experiments/harness/TRIES.md holds row T45. Arms per seed, all from one cell draw:
         b0, m1_r_gmom (fraction 0.2, the entry's form) and m1_r_gmom_x2 (fraction 0.4) as T44 builds them; T45a
         m1_r_gmom_x0, the same build at fraction 0 (external_shift accepts 0: every factor is exp(0) = 1, so no count
         is redrawn); T45b m1_r_gmom_p1 and T45c m1_r_gmom_p2, the builds at 0.2 and 0.4 with try T15's permuted control
         (harness.t3_decile_permutation within LP_BINS bins of feng_dirs.direction_halves' mean_cp on M1's bracket,
         seed [LP_PERM_SEED, LP_GROUPS.index(group)]) in place of Feng's directions, one permutation per group for both
         fractions and every seed; T45d m1_r_gmom_x2b, "there and back": the build at 0.4 with, before g_mom's
         factors, a second application of Feng's directions at -0.2 to the covered cells (gmom_combo.moved_by_feng:
         external_shift's own rule, R2's held genes) on stream [seed, SECOND_STREAM], through build_submission's
         post_cells hook. Each arm is gmom_combo.m1_r_arms' m1_r_gmom with only its fraction, its directions or (T45d)
         the hook's first edit changed. The task file holds every arm's fraction(s), each group's sha256 of Feng's
         direction and of the directions object each placebo arm was built with, each group's applied size (mean
         |log f| and sum of log f^2 over the non-held genes at the base fraction) for Feng and the placebo, and per
         seed and arm the target-free redraw measures of the covered cells before g_mom against the drawn cells: the
         share of non-zero counts changed, the share of non-zero counts set to zero, and the mean detected genes per
         cell (drawn and edited).
  read   T45's reading of both cuts' files (rule below, fixed in the row before the run).

What a control shares with its twin and does not: every arm reaches external_shift with the same cells, groups, held
genes and stream ([seed, 2]); downstream each covered cell's keep-depth constant follows the step, and with it which
counts are thinned or topped up, and g_mom's factors (stream [seed, 5]) then act on different counts. A per-seed
difference includes resampling noise the arms do not share. T45d's first edit is bit for bit the x2 arm's; its second
edit redraws the covered cells' counts again: held genes keep their counts; every other gene of a covered cell, NaN-
direction genes included, is redrawn, because the depth constant scales its factor of 1. x2b's expected counts are the
base's, but on a synthetic build it carries about 1.6 times x2's per-entry count variance and zeroing, and 2.4 times
x2's excess over the base. (d) passing measures the gain from that much redraw, not from the larger step's, and (c) is
shifted by the difference, so either count-noise value overstates the side the noise favours; the read prints the
redraw measures beside (c) and (d). The placebo permutes held and non-held genes together, as T15's did; the read
prints each group's applied size for Feng and the placebo, and a difference of more than 10% is named beside (a) and
(b).

NOT an adoption result: a labelled in-sample diagnostic, as tries T26, T42 and T44 are. Feng's E9.5 cells lie inside
these pseudo-tasks' external windows, the cuts are ones earlier T1 choices were made on, and the truth of a one-day cut
differs from its reference in the kind of sample as well as in age. The result can only go to the team captain.

The reading. Void unless both cuts hold the ten seeds at one commit on a clean tree; b0 and m1_r_gmom equal, seed for
seed, their committed scores in experiments/t1_m1/results/t42_m1_eb/tasks/ and m1_r_gmom_x2 its committed scores in
experiments/t1_m1/results/t44_m1_step/tasks/, each to 1e-9; each file's task is its cut, its base fraction 0.2 and each
arm's fraction(s) as above; no control scores exactly as m1_r_gmom on every seed of both cuts; m1_r_gmom_p1 and
m1_r_gmom_p2 record the same permuted directions, and these differ from Feng's in every group. The test, used eight
times, is the adoption rule's (vec.adoption.keep_change, paired by seed, both cuts); "A over B" is A as the change, B as
the base: (a) compares the larger step's increment with the placebo's increment, (x2 - m1_r_gmom) against (p2 - p1)
paired by seed; the level comparison x2 over p2 is reported as (a0) and decides nothing; (b) p2 over p1; (c) x2 over
x2b; (d) x2b over m1_r_gmom; (e) m1_r_gmom over p1; (f) p1 over x0; (g) m1_r_gmom over x0; (h) x2 over m1_r_gmom, the
adoption rule on these files (it reproduces T44's comparison). The two statements are read only if (h) passes;
otherwise both read "no larger-step gain to explain" and the comparisons are reported only. Gene identity: "a placebo
gains as well; the direction more" if (a) and (b) pass; "a placebo gains as well" if (b) passes and (a) does not;
"needs the direction's genes" if (a) passes and (b) does not; otherwise "not separated". Count noise: "count noise of
x2b's size gains" if (d) passes and (c) does not; "is not count noise" if (c) passes and (d) does not; "part of each"
if both pass (the text gives the two sizes); otherwise "not separated". (a0), (e), (f) and (g) decide nothing.
Reported: per cut every arm's mean board score, each comparison's mean paired difference with seeds ahead and per-seed
range, every metric's mean skill difference, the redraw measures beside (c) and (d), the placebo's applied size beside
(a) and (b), and beside (b) and (d) try T43c's committed difference (b0_mn against b0 in
experiments/harness/results/t43_t1_mix/tasks/), what a full count redraw does on the same cuts. The read puts nothing
forward by itself.

Usage: on the DGX campaign checkout, repo root, both cuts at one commit, committing nothing until both task files
exist, and on the machine T42 and T44 ran on, so that b0, m1_r_gmom and m1_r_gmom_x2 reproduce. The environment
variables are docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/perm_read.py run --out t45_m1_perm \
      --task "T1/E9.5->E10.5:chamber"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/perm_read.py run --out t45_m1_perm \
      --task "T1/E9.5->E10.5"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/perm_read.py read --out t45_m1_perm
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import re
import sys
import time
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT, load_boards
from vec.core import counts_of
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
T42_TASKS = RESULTS / "t42_m1_eb" / "tasks"
T44_TASKS = RESULTS / "t44_m1_step" / "tasks"
T43_TASKS = REPO_ROOT / "experiments" / "harness" / "results" / "t43_t1_mix" / "tasks"
sys.path.insert(0, str(HERE))
import eb_read as er  # noqa: E402
import feng_dirs  # noqa: E402
import gmom_combo as gc  # noqa: E402
import step_read as sr  # noqa: E402

gv, hygiene, n12 = gc.gv, gc.hygiene, er.n12
CUT_FILES = gc.CUT_FILES
ROW = "T45"
BASE = sr.BASE  # m1_r_gmom
FRACTION = er.FRACTION  # 0.2, Feng's step on an E9.5 -> E10.5 cut
X2, X0, P1, P2, X2B = "m1_r_gmom_x2", "m1_r_gmom_x0", "m1_r_gmom_p1", "m1_r_gmom_p2", "m1_r_gmom_x2b"
CONTROLS = (X0, P1, P2, X2B)
ARMS = ("b0", BASE, X2, *CONTROLS)
PLACEBO_ARMS = (P1, P2)
SECOND_STREAM = 15  # T45d's second edit: an index no existing step uses (build: 1, 2, 4; g_mom: 5; others 0-14, 20+)
FRACTION_TOL = sr.FRACTION_TOL
REPRO_TOL = sr.REPRO_TOL  # 1e-9: same machine and code path
REPRO = {"t42": ("b0", BASE), "t44": (X2,)}
M1_ARMS = (BASE, X2, *CONTROLS)  # the arms built through gmom_combo.m1_r_arms, in build order
X2_INC, P2_INC = f"{X2} - {BASE}", f"{P2} - {P1}"  # (a)'s paired increments: the larger step's and the placebo's
COMPARISONS = {"a": (P2_INC, X2_INC), "a0": (P2, X2), "b": (P1, P2), "c": (X2B, X2), "d": (BASE, X2B),
               "e": (P1, BASE), "f": (X0, P1), "g": (X0, BASE), "h": (BASE, X2)}  # letter -> (base, change)
REPORTED_ONLY = ("a0", "e", "f", "g")  # reported, decide nothing; (h) gates the two statements
T43C = ("b0", "b0_mn")  # (base, change): try T43c's full count redraw, reported beside (b) and (d)
NEEDS_GENES, PLACEBO_GAINS, PLACEBO_LESS, NOT_SEPARATED = (
    "needs the direction's genes", "a placebo gains as well", "a placebo gains as well; the direction more",
    "not separated")
IS_NOISE, NOT_NOISE, PART_OF_EACH = "count noise of x2b's size gains", "is not count noise", "part of each"
NO_GAIN = "no larger-step gain to explain"
SIZE_FLAG = 0.10  # the reviewer's 10%: a placebo applied size this far from Feng's is named beside (a) and (b)
REDRAW_KEYS = ("share_nonzero_changed", "share_nonzero_zeroed", "detected_per_cell", "detected_per_cell_drawn")
LABEL_T45 = ("try T45, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) at fractions 0, 0.2 and 0.4 of Feng's bracket, with try T15's permuted placebo at 0.2 "
             "and 0.4, and a there-and-back step (0.4 then -0.2), on T44's route; Feng's E9.5 cells lie inside these "
             "pseudo-tasks' external windows, the cuts are ones earlier T1 choices were made on, and the truth of a "
             "one-day cut differs from its reference in the kind of sample as well as in age; the result can only go "
             "to the team captain")


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def arm_fractions(fraction: float) -> dict[str, float | list[float]]:
    """Every M1 arm's fraction; T45d's two, in the order applied."""
    f = float(fraction)
    return {BASE: f, X2: f * 2.0, X0: 0.0, P1: f, P2: f * 2.0, X2B: [f * 2.0, -f]}


def placebo_directions(directions: Mapping[str, np.ndarray], halves: Mapping[str, Mapping[str, np.ndarray]]
                       ) -> dict[str, np.ndarray]:
    """Try T15's permuted control (harness._late_proxy with permute=True), unchanged, on M1's directions: per group, the
    finite values permuted across genes within LP_BINS bins of the group's mean CP10k on M1's bracket, seed
    [LP_PERM_SEED, LP_GROUPS.index(group)]. NaN stays NaN. A group outside LP_GROUPS or without halves is refused."""
    out = {}
    for g, d in directions.items():
        if g not in harness.LP_GROUPS:
            raise ValueError(f"{g}: not one of try T15's groups {list(harness.LP_GROUPS)}: it has no permutation seed")
        if g not in halves:
            raise ValueError(f"{g}: Feng's halves give no mean CP10k for this direction: it cannot be permuted")
        d = np.asarray(d, dtype=np.float64).copy()
        ok = np.isfinite(d)
        d[ok] = harness.t3_decile_permutation(d[ok], np.asarray(halves[g]["mean_cp"], dtype=np.float64)[ok],
                                              [harness.LP_PERM_SEED, harness.LP_GROUPS.index(g)], harness.LP_BINS)
        out[g] = d
    return out


def direction_sha(directions: Mapping[str, np.ndarray]) -> dict[str, str]:
    return {g: hashlib.sha256(np.ascontiguousarray(d, dtype=np.float64).tobytes()).hexdigest()
            for g, d in sorted(directions.items())}


def applied_size(directions: Mapping[str, np.ndarray], held: Sequence[str], genes: Sequence[str], fraction: float
                 ) -> dict[str, dict[str, float]]:
    """Per group, the log factors external_shift applies at `fraction` (NaN as 0, held genes left out): their mean
    absolute value and their sum of squares over the non-held genes. Target-free."""
    free = ~np.isin(np.asarray(list(genes), dtype=object), list(held))
    out = {}
    for g, d in sorted(directions.items()):
        logf = float(fraction) * np.nan_to_num(np.asarray(d, dtype=np.float64), nan=0.0)[free]
        out[g] = {"mean_abs_logf": float(np.abs(logf).mean()) if logf.size else 0.0,
                  "sum_sq_logf": float((logf * logf).sum())}
    return out


def _dense(M) -> np.ndarray:
    return M.toarray() if hasattr(M, "toarray") else np.asarray(M)


def redraw_measures(drawn, edited) -> dict[str, float]:
    """Target-free measures of how much of the covered cells' counts an arm redrew: among the drawn non-zero counts,
    the share changed and the share set to zero; the mean detected genes per cell, edited and drawn."""
    B, A = _dense(drawn).astype(np.int64), _dense(edited).astype(np.int64)
    if A.shape != B.shape:
        raise ValueError(f"the edited covered cells {A.shape} are not the drawn ones {B.shape}")
    nz = B > 0
    n = int(nz.sum())
    return {"share_nonzero_changed": float((A[nz] != B[nz]).sum() / n) if n else 0.0,
            "share_nonzero_zeroed": float((A[nz] == 0).sum() / n) if n else 0.0,
            "detected_per_cell": float((A > 0).sum(axis=1).mean()) if A.shape[0] else 0.0,
            "detected_per_cell_drawn": float(nz.sum(axis=1).mean()) if B.shape[0] else 0.0}


@contextlib.contextmanager
def _redraw_probe(directions: Mapping, groups: Mapping, out: list) -> Iterator[None]:
    """For the duration, append to `out`, at every g_mom hook call, redraw_measures of the covered cells as drawn (the
    input of the build's own external_shift call, stream [seed, 2]) against the same cells as they enter g_mom. Every
    value passes through unchanged. Restored on exit."""
    real_shift, real_hook = t1_scope.external_shift, gc.gmom_hook
    drawn = []

    def shift(target_like, *args, **kwargs):
        if [int(x) for x in np.atleast_1d(kwargs.get("seed", 0))][-1:] == [2]:
            drawn.append(counts_of(target_like).copy())
        return real_shift(target_like, *args, **kwargs)

    def hook(*args, **kwargs):
        gmom = real_hook(*args, **kwargs)

        def probed(pop):
            labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
            covered = np.flatnonzero([groups.get(lab) in directions for lab in labels])
            out.append(redraw_measures(drawn[-1], counts_of(pop)[covered]))
            return gmom(pop)

        return probed

    t1_scope.external_shift, gc.gmom_hook = shift, hook
    try:
        yield
    finally:
        t1_scope.external_shift, gc.gmom_hook = real_shift, real_hook


def second_edit(directions: Mapping, groups: Mapping, fraction: float, held: Sequence[str], seed: int):
    """T45d's second edit: Feng's directions at `fraction` on the covered cells (gmom_combo.moved_by_feng, through
    external_shift with its own rule), R2's held genes held, stream [seed, SECOND_STREAM]; other cells as they are."""
    def edit(pop):
        return gc.moved_by_feng(pop, directions, fraction, groups, held, [seed, SECOND_STREAM])

    return edit


def placebo_record(built_sha: Mapping[str, set[str]]) -> dict:
    """Per placebo arm, the group -> sha256 of the directions it was built with when every seed built it with the same
    ones; otherwise {DIFFERS: [each seed's record]}, which the reading voids."""
    return {a_: (json.loads(next(iter(v))) if len(v) == 1 else {DIFFERS: sorted(v)}) for a_, v in built_sha.items()}


DIFFERS = "differs between seeds"


@contextlib.contextmanager
def _before_gmom(edit) -> Iterator[None]:
    """For the duration, gmom_combo.m1_r_arms' g_mom hook runs `edit` on the population first: the post_cells hook
    becomes g_mom(edit(pop)). Restored on exit."""
    real = gc.gmom_hook

    def hook(*args, **kwargs):
        gmom = real(*args, **kwargs)
        return lambda pop: gmom(edit(pop))

    gc.gmom_hook = hook
    try:
        yield
    finally:
        gc.gmom_hook = real


def there_and_back(ctx: harness.BuildContext, directions: Mapping, groups: Mapping, fraction: float,
                   held_m1: Sequence[str], held_factors: Sequence[str]) -> Stage:
    """T45d: gmom_combo.m1_r_arms at 2 * fraction (its first edit is the x2 arm's) with the second edit at -fraction
    put in build_submission's post_cells hook before g_mom's factors."""
    with _before_gmom(second_edit(directions, groups, -float(fraction), held_m1, ctx.seed)):
        return gc.m1_r_arms(ctx, directions, groups, 2.0 * float(fraction), held_m1, held_factors)[BASE]


def perm_arms(ctx: harness.BuildContext, directions: Mapping, placebo: Mapping, groups: Mapping, fraction: float,
              held_m1: Sequence[str], held_factors: Sequence[str]) -> dict:
    """The six M1 arms of one seed, each gmom_combo.m1_r_arms' m1_r_gmom as step_read.step_arms builds it: on Feng's
    directions (base, x2, x0) and on the placebo (p1, p2), and T45d. The same context, groups, held genes and streams
    for each. Returns {"arms": arm -> Stage, "directions": arm -> the directions object the arm was built with,
    "redraw": arm -> redraw_measures of its covered cells before g_mom}."""
    f = float(fraction)
    spec = {BASE: (directions, f), X2: (directions, 2.0 * f), X0: (directions, 0.0), P1: (placebo, f),
            P2: (placebo, 2.0 * f), X2B: (directions, None)}
    arms, used, redraw = {}, {}, {}
    for arm in M1_ARMS:
        dirs, fr = spec[arm]
        rec = []
        with _redraw_probe(dirs, groups, rec):
            if arm == X2B:
                arms[arm] = there_and_back(ctx, dirs, groups, f, held_m1, held_factors)
            else:
                arms[arm] = gc.m1_r_arms(ctx, dirs, groups, fr, held_m1, held_factors)[BASE]
        used[arm], redraw[arm] = dirs, (rec[-1] if rec else None)
    return {"arms": arms, "directions": used, "redraw": redraw}


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T45 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T45 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = re.sub(r"[^A-Za-z0-9.]+", "_", a.task).strip("_") + ".json"
    gv.refuse_existing(out / "tasks", (name,))
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
    halves, stats_sha = er.feng_halves(list(pool.var_names), directions)
    if stats_sha != ext["sha256"]:
        sys.exit("the statistics the halves came from are not the ones M1's direction came from")
    placebo = placebo_directions(directions, halves)  # once: the same for both fractions and every build seed
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows, built_sha = [], {a_: set() for a_ in PLACEBO_ARMS}
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        built = perm_arms(ctx, directions, placebo, groups, fraction, held_m1, preregistered)
        for a_ in PLACEBO_ARMS:  # the directions object each placebo arm was built with, hashed on every seed
            built_sha[a_].add(json.dumps(direction_sha(built["directions"][a_]), sort_keys=True))
        arms = er.score_arms(board, target, ref_stage, s, f.anchors, {"b0": b0(ctx), **built["arms"]})
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "redraw": built["redraw"],
                     "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res = {"label": LABEL_T45, "task": t.key, "fraction": fraction, "fractions": arm_fractions(fraction),
           "second_stream": SECOND_STREAM, "feng_sha256": direction_sha(directions),
           "placebo_sha256": placebo_record(built_sha),
           "applied_size": {"fraction": fraction, "feng": applied_size(directions, held_m1, pool.var_names, fraction),
                            "placebo": applied_size(placebo, held_m1, pool.var_names, fraction)},
           "feng_stats_sha256": stats_sha, "external": ext, "seeds": list(seeds), "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, REPO_ROOT / union["file"],
                                       feng_dirs.clean_stats_path(), n22])}
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(rows: Sequence[Mapping], ref_rows: Sequence[Mapping], arms: Sequence[str]) -> dict:
    """Per seed, `arms` against a committed file's scores; holds when every seed is there and each is within
    REPRO_TOL."""
    by_seed = {r["seed"]: r for r in ref_rows}
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
    holds = bool(ref_rows) and not missing and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "missing_seeds": missing, "absent": absent, "holds": holds}


def _fraction_ok(got, want) -> bool:
    if isinstance(want, list):
        return (isinstance(got, list) and len(got) == len(want)
                and all(_fraction_ok(g, w) for g, w in zip(got, want, strict=True)))
    return isinstance(got, int | float) and not isinstance(got, bool) and abs(got - want) <= FRACTION_TOL


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Mapping]]) -> dict:
    """The void conditions. results: cut -> task file (a missing cut is absent); committed: "t42" and "t44" -> cut ->
    their committed task files."""
    reasons = []
    missing = [k for k in CUT_FILES if k not in results]
    if missing:
        reasons.append(f"missing cuts {missing}")
    want_fr = arm_fractions(FRACTION)
    for k, r in results.items():
        seeds = [row["seed"] for row in r["rows"]]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in r["rows"] for a in ARMS if a not in row["arms"]})
        if absent:
            reasons.append(f"{k}: arms {absent} missing")
        if r.get("task") != k:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if not _fraction_ok(r.get("fraction"), FRACTION):
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items() if not _fraction_ok(got_fr.get(a), f)]
        if bad:
            reasons.append(f"{k}: arm fractions {got_fr}, not {want_fr} (wrong fraction on {bad})")
        feng, placebo = r.get("feng_sha256") or {}, r.get("placebo_sha256") or {}
        recorded = [placebo.get(a) for a in PLACEBO_ARMS]
        if any(isinstance(p, Mapping) and DIFFERS in p for p in recorded):
            reasons.append(f"{k}: a placebo arm was built with different directions on different seeds")
        elif any(not p for p in recorded) or any(p != recorded[0] for p in recorded):
            reasons.append(f"{k}: {P1} and {P2} do not record the same permuted directions")
        elif not feng or set(recorded[0]) != set(feng):
            reasons.append(f"{k}: the placebo's groups {sorted(recorded[0])} are not Feng's {sorted(feng)}")
        else:
            same = sorted(g for g in feng if recorded[0][g] == feng[g])
            if same:
                reasons.append(f"{k}: the placebo equals Feng's direction in groups {same}")
    commits = {r["stamp"]["commit"] for r in results.values()}
    unclean = sorted(k for k, r in results.items() if r["stamp"]["tree_clean"] is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the cuts ran at {len(commits)} commits")
    repro = {}
    for k, r in results.items():
        repro[k] = {}
        for src, arms in REPRO.items():
            repro[k][src] = reproduction(r["rows"], committed.get(src, {}).get(k, {}).get("rows", []), arms)
            if not repro[k][src]["holds"]:
                reasons.append(f"{k}: {', '.join(arms)} do not reproduce {src.upper()}'s committed scores within "
                               f"{REPRO_TOL:g} per seed")
    complete = all(a in row["arms"] for r in results.values() for row in r["rows"] for a in ARMS)
    if results and complete:
        same = [a for a in CONTROLS
                if all(row["arms"][a] == row["arms"][BASE] for r in results.values() for row in r["rows"])]
        if same:
            reasons.append(f"controls score exactly as {BASE} on every seed of both cuts: {same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def statements(p: Mapping[str, bool]) -> dict[str, str]:
    """The two statements about the larger step, from which comparisons pass; both read NO_GAIN unless (h), the larger
    step over the base on these files, passes."""
    if not p["h"]:
        return {"gene identity": NO_GAIN, "count noise": NO_GAIN}
    if p["b"]:
        genes = PLACEBO_LESS if p["a"] else PLACEBO_GAINS
    elif p["a"]:
        genes = NEEDS_GENES
    else:
        genes = NOT_SEPARATED
    if p["d"] and p["c"]:
        noise = PART_OF_EACH
    elif p["d"]:
        noise = IS_NOISE
    elif p["c"]:
        noise = NOT_NOISE
    else:
        noise = NOT_SEPARATED
    return {"gene identity": genes, "count noise": noise}


def _arm_difference(change: Mapping, base: Mapping) -> dict:
    """change - base for an arm record's score and every skill, group and metric both hold."""
    out = {"score": float(change["score"]) - float(base["score"])}
    for f in ("skills", "groups", "metrics"):
        keys = sorted(set(change.get(f, {})) & set(base.get(f, {})))
        if keys:
            out[f] = {m: float(change[f][m]) - float(base[f][m]) for m in keys}
    return out


def with_increments(results: Mapping[str, Mapping]) -> dict:
    """The task files with (a)'s two paired increments as arms of each row: X2_INC = x2 - m1_r_gmom and P2_INC =
    p2 - p1, seed by seed, so that (a) is the adoption rule on (P2_INC, X2_INC) as the (base, change) pair."""
    return {k: {**r, "rows": [{**row, "arms": {**row["arms"],
                                               X2_INC: _arm_difference(row["arms"][X2], row["arms"][BASE]),
                                               P2_INC: _arm_difference(row["arms"][P2], row["arms"][P1])}}
                              for row in r["rows"]]}
            for k, r in results.items()}


def size_flags(sizes: Mapping | None) -> list[str]:
    """Groups whose placebo applied size (mean |log f| or sum of log f^2, non-held genes) differs from Feng's by more
    than SIZE_FLAG, relative to Feng's."""
    if not sizes:
        return []
    out = []
    for g, fs in sorted(sizes.get("feng", {}).items()):
        ps = sizes.get("placebo", {}).get(g, {})
        for m, x in fs.items():
            y = ps.get(m)
            if y is None or (x == 0 and y != 0) or (x != 0 and abs(y - x) / abs(x) > SIZE_FLAG):
                out.append(f"{g} {m}: Feng {x:.4g}, placebo {y if y is None else f'{y:.4g}'}")
    return out


def mean_redraw(results: Mapping[str, Mapping]) -> dict:
    """Per cut and M1 arm, the mean over seeds of each redraw measure; an arm or cut without the record is left out."""
    out = {}
    for k, r in results.items():
        per = {}
        for a in M1_ARMS:
            recs = [row.get("redraw", {}).get(a) for row in r["rows"]]
            if recs and all(isinstance(x, Mapping) for x in recs):
                per[a] = {m: float(np.mean([x[m] for x in recs])) for m in REDRAW_KEYS if all(m in x for x in recs)}
        if per:
            out[k] = per
    return out


def harness_rows(task_file: Mapping, arms: Sequence[str]) -> dict:
    """A harness task file (rows of per-candidate scores and skills) in the arms form eb_read.compare reads."""
    rows = [{"seed": r["seed"], "arms": {a: {"score": float(r["scores"][a]),
                                             **({"skills": dict(r["skills"][a])}
                                                if isinstance(r.get("skills", {}).get(a), Mapping) else {})}
                                         for a in arms}}
            for r in task_file["rows"]]
    return {"rows": rows}


def t43c(t43: Mapping[str, Mapping]) -> dict | None:
    """Try T43c's committed difference, b0_mn against b0, paired by seed on both cuts; None if a cut or an arm is
    missing."""
    if not all(k in t43 for k in CUT_FILES):
        return None
    if not all(a in r["scores"] for k in CUT_FILES for r in t43[k]["rows"] for a in T43C):
        return None
    return sr.compare({k: harness_rows(t43[k], T43C) for k in CUT_FILES}, *T43C)


def t45_reading(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Mapping]],
                t43: Mapping[str, Mapping] | None = None) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    inc = with_increments(results)
    comps = {c: sr.compare(inc, base, change) for c, (base, change) in COMPARISONS.items()}
    out["comparisons"] = comps
    out["statements"] = statements({c: comps[c]["passes"] for c in comps})
    out["outcome"] = "; ".join(f"{k}: {s}" for k, s in out["statements"].items())
    if out["statements"]["count noise"] == PART_OF_EACH:
        out["count_noise_note"] = ("both (c) and (d) pass: part of each; sizes per cut, (c) "
                                   + ", ".join(f"{c_['mean']:+.2f}" for c_ in comps["c"]["cuts"].values()) + "; (d) "
                                   + ", ".join(f"{c_['mean']:+.2f}" for c_ in comps["d"]["cuts"].values()))
    if out["statements"]["gene identity"] == PLACEBO_LESS:
        out["gene_identity_note"] = ("both (a) and (b) pass; sizes per cut, (a) "
                                     + ", ".join(f"{c_['mean']:+.2f}" for c_ in comps["a"]["cuts"].values())
                                     + "; (b) "
                                     + ", ".join(f"{c_['mean']:+.2f}" for c_ in comps["b"]["cuts"].values()))
    out["mean_scores"] = {k: {a: float(np.mean([row["arms"][a]["score"] for row in r["rows"]])) for a in ARMS}
                          for k, r in results.items()}
    out["redraw"] = mean_redraw(results)
    out["applied_size"] = {k: r.get("applied_size") for k, r in results.items()}
    out["size_flags"] = {k: size_flags(r.get("applied_size")) for k, r in results.items()}
    t43c_comp = t43c(t43 or {})
    if t43c_comp is not None:
        out["t43c"] = t43c_comp
    return out


def t45_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"Outcome: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, by_src in v["reproduction"].items():
        lines.append(f"{k}: largest per-seed difference from the committed scores: " + ", ".join(
            f"{a} {w:.3g} ({src.upper()})" for src, g in by_src.items() for a, w in g["largest_difference"].items())
            + f" (gate {REPRO_TOL:g})")
    for k, ms in r.get("mean_scores", {}).items():
        lines.append(f"mean scores {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in ms.items()))
    if "gene_identity_note" in r:
        lines.append(f"Gene identity: {r['gene_identity_note']}")
    if "count_noise_note" in r:
        lines.append(f"Count noise: {r['count_noise_note']}")
    if r.get("statements", {}).get("count noise") == NO_GAIN:
        lines.append("(h), the larger step over m1_r_gmom on these files, does not pass: the two statements are not "
                     "read and the comparisons are reported only")

    def redraw_lines(arms: Sequence[str]) -> None:
        for k, per in r.get("redraw", {}).items():
            lines.append(f"    redraw of the covered cells before g_mom, {k}: " + "; ".join(
                f"{a} changed {per[a]['share_nonzero_changed']:.4f}, zeroed {per[a]['share_nonzero_zeroed']:.4f}, "
                f"detected per cell {per[a]['detected_per_cell']:.1f} (drawn {per[a]['detected_per_cell_drawn']:.1f})"
                for a in arms if a in per and len(per[a]) == len(REDRAW_KEYS)))
        if not r.get("redraw"):
            lines.append("    redraw measures: not recorded in the task files")

    def size_lines() -> None:
        for k, flags in r.get("size_flags", {}).items():
            sizes = (r.get("applied_size") or {}).get(k)
            if not sizes:
                lines.append(f"    placebo applied size, {k}: not recorded")
                continue
            lines.append(f"    applied size at fraction {sizes.get('fraction')}, {k} (mean |log f|, sum log f^2, "
                         "non-held genes): " + "; ".join(
                             f"{g} Feng {fs['mean_abs_logf']:.4f}, {fs['sum_sq_logf']:.3f}, placebo "
                             f"{sizes['placebo'][g]['mean_abs_logf']:.4f}, {sizes['placebo'][g]['sum_sq_logf']:.3f}"
                             for g, fs in sorted(sizes.get("feng", {}).items()) if g in sizes.get("placebo", {})))
            if flags:
                lines.append(f"    the placebo's applied size differs from Feng's by more than {SIZE_FLAG:.0%}: "
                             + "; ".join(flags))

    def block(name: str, comp: Mapping) -> None:
        lines.append(f"{name}: {comp['change']} over {comp['base']}: "
                     f"{'PASSES' if comp['passes'] else 'does not pass'} ({comp['rule']['reason']})")
        for k, cut in comp["cuts"].items():
            sk = ", ".join(f"{m} {x:+.4f}" for m, x in cut["skills"].items())
            lines.append(f"    {k}: {cut['mean']:+.2f} points, ahead {cut['ahead']} of {cut['of']}, per seed "
                         f"{cut['range'][0]:+.2f} to {cut['range'][1]:+.2f}; skills: {sk}")

    for c, comp in r.get("comparisons", {}).items():
        note = {"h": " gates the two statements"}.get(c, " reported, decides nothing" if c in REPORTED_ONLY else "")
        block(f"({c}){note}", comp)
        if c in ("a", "b"):
            size_lines()
        if c in ("c", "d"):
            redraw_lines((X2, X2B) if c == "c" else (X2B, BASE))
        if c in ("b", "d") and "t43c" in r:
            block("    beside it, try T43c's full count redraw (committed files, reported)", r["t43c"])
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def _load(folder: Path) -> tuple[dict, list[Path]]:
    paths = {k: folder / f for k, f in CUT_FILES.items() if (folder / f).exists()}
    return {k: json.loads(p.read_text()) for k, p in paths.items()}, list(paths.values())


def read(out: Path, t42_tasks: Path = T42_TASKS, t44_tasks: Path = T44_TASKS, t43_tasks: Path = T43_TASKS) -> dict:
    """Read the cuts present in out/tasks against T42's and T44's committed task files, with T43c's reported beside;
    write t45.{json,txt}."""
    gv.refuse_existing(out, ("t45.json", "t45.txt"))
    started = gv.tree_state()
    results, present = _load(out / "tasks")
    t42, t42_paths = _load(t42_tasks)
    t44, t44_paths = _load(t44_tasks)
    t43, t43_paths = _load(t43_tasks)
    res = {"label": LABEL_T45, "reading": t45_reading(results, {"t42": t42, "t44": t44}, t43),
           "stamps": {k: v["stamp"] for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present, *t42_paths, *t44_paths, *t43_paths])
    gv.write_outputs(out, "t45", res, t45_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T45's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one cut of try T45")
    p.add_argument("--task", required=True, choices=list(CUT_FILES))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T45's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
