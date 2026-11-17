#!/usr/bin/env python
"""Try T49 (experiments/harness/TRIES.md): the size of g_mom's within-type factors on the two counted T1 cuts.

g_mom (try T11) takes, per reference-stage metacell, the log ratio of its mean CP10k to its transported start-stage
partners' mean (pseudocount 0.5 CP10k), removes the category's mean, damps by B0's 0.5, zeroes below log 1.05 and clips
at log 3 (`vec.t1_state_transport.finish_momentum`). Its controls (tries T20 to T25) read its sign, its assignment, a
no-development twin and its composition, never its size. This read puts two larger dampings, 0.75 and 1.0 (1.5 and 2
times B0's), to the same cuts, in try T11's own form and in the entry's form.

  run    One E10.5 cut, 10 seeds, the harness's draws and scoring exactly as try T42's run. The fold's raw fit is made
         once per run by `harness._g_mom_factors` itself (`fold_factors` clears its cache, calls it, and keeps the fit
         it hands to `finish_momentum`), then finished three times with `finish_momentum` at 0.5, 0.75 and 1.0 (the same
         zeroing and clip, acting on the larger values); the one at 0.5 is `_g_mom_factors`' own object, and the run
         stops unless a fresh finish at 0.5 equals it bit for bit. Arms per seed, all from one cell draw:
           form A, try T11's (B0's cells and shift with the factors on top, `harness._g_mom_arm` with a transform that
           returns the factors finished at the arm's damping): b0, g_mom (0.5), g_mom_d075, g_mom_d1;
           form B, the entry's (try T42's m1_r_gmom, `gmom_combo.m1_r_arms`, with the factors its `gmom_hook` applies
           replaced by the ones finished at the arm's damping and nothing else changed): m1_r_gmom (0.5),
           m1_r_gmom_g075, m1_r_gmom_g1.
         The task file holds every arm's damping, the reproduction tolerance and, per damping, target-free records of
         the finished factors over the free genes (every fitted gene not on the pre-registered held-flat list, which
         is held at factor 1 when applied; genes off the board are kept): per category the share of entries not
         zeroed, the share at the clip (|log f| at or above log 3), the mean absolute log factor over the entries not
         zeroed and over all entries (zeros included), and the share of entries zeroed at 0.5 and not zeroed at that
         damping. Refused on an existing task file, unless the committed N22 passed, and unless
         experiments/harness/TRIES.md holds row T49. No arm holds T23's no-development factors at the larger sizes:
         the part of a larger gain that needs no development is not controlled (see the limits).
  read  T49's reading of both cuts' files (rule below, fixed in the row before the run).

What a larger arm shares with its base and what it does not: the same fit object, counts, cells, held genes, R1's
unedited labels, rescale and stream index ([seed, 5] for the factors; in form B also M1's [seed, 2]); the factors
differ. But `apply_factors` draws all of a call's Binomial and Poisson values from that one stream, so a changed factor
anywhere misaligns every later draw: on the tests' synthetic fold 71 to 74% of a larger arm's output entries differ
from its base's, not only the resized entries. Form B builds m1_r as a by-product of each `m1_r_arms` call; it is
discarded.

In both forms the read is a labelled in-sample diagnostic: form A selects a constant (B0's damping applied to g_mom's
factors) on the cuts that kept g_mom; form B, the entry's form, is in-sample as tries T26, T42 and T44 are (R1, R2 and
the held genes were chosen on these cuts, and Feng's E9.5 cells lie inside these pseudo-tasks' external windows).

The reading. Void unless both cuts hold the ten seeds at one commit on a clean tree; b0 and g_mom equal, seed for seed,
their committed scores in experiments/harness/results/t11_g_mom/tasks/ and m1_r_gmom its committed scores in
experiments/t1_m1/results/t42_m1_eb/tasks/, to 1e-9; each file's task is its cut and its dampings 0.5, 0.75 and 1.0;
and no larger arm scores exactly as its base (the same record) on every seed of both cuts. The arithmetic is the
adoption rule's (vec.adoption.keep_change, paired by seed, both cuts: a mean gain above 0.3 with at least 8 of 10 seeds
ahead on each, no mean loss beyond 0.3); applied to a constant chosen on these cuts it keeps nothing. Form A:
g_mom_d075 against g_mom (A1), g_mom_d1 against g_mom (A2), g_mom_d1 against g_mom_d075 (A3). The arm named: g_mom_d1
if it passes the arithmetic against g_mom and against g_mom_d075; otherwise g_mom_d075 if it passes against g_mom;
otherwise g_mom_d1 if it alone passes against g_mom. Form B, the same three comparisons on the entry's form (B1, B2,
B3), a labelled diagnostic with one use: the arm of form A is named only if its twin in form B does not lose more than
the adoption rule's 0.3 (vec.adoption.MAX_LOSS) to m1_r_gmom on either cut (its mean difference at or above -0.3 on
both cuts). Statement: "a larger g_mom scores above g_mom by the adoption rule's arithmetic on the cuts that kept g_mom
(a size chosen on these cuts; not a keep): <arm>"; or, when the twin loses beyond 0.3, "above g_mom by the adoption
rule's arithmetic in try T11's form and loses more than 0.3 in the entry's form: not put forward (<arm>)"; with no arm,
"a larger g_mom loses (A1 <chamber mean> / <whole mean>, A2 <chamber mean> / <whole mean>)" if both larger arms' mean
differences from g_mom are below zero on both cuts with at most 2 of 10 seeds ahead on each; otherwise "not separated
(largest mean gain <x>)", x the largest of A1's and A2's mean differences over the two cuts. Reported: per cut every
arm's mean board score; each comparison's mean paired difference with seeds ahead and per-seed range; every metric's
mean skill difference and its part in board points (100 x the group weight x the skill difference: the two de terms,
mmd_u and the variogram); the factors' records; the limits below.

Limits (written before the run; each goes into the ledger row and into t49.txt):
- T49 chooses a value of one constant (B0's damping applied to g_mom) on the two cuts that kept g_mom and have read it
  seven times. Whatever it reads, it keeps nothing under rule 6, it is counted as tries, and its statement is a
  labelled in-sample diagnostic in both forms.
- T47's 0.58 is a calibration of Feng's direction over M1's genes. It is not measured for g_mom's factors, whose
  zeroing removes most low-mean entries. g_mom's raw fit holds sampling noise, transport mismatch and the amplification
  of each state's own deviation as well as development (T23: about 40% of the gain with nothing developing; N21: R^2
  0.35 to 0.64 from amplification). So a delivered share below 1 is expected of a noisy estimate and is not evidence
  that the size is too small. No outcome is read as 'g_mom was under-sized by a factor of about 3'.
- Form B's veto reads the sign of an in-sample twin, and M1's step there is itself too small (T44). A veto or its
  absence says nothing about g_mom's size in the entry at test.
- A larger arm redraws every count of every edited cell on a stream it does not share with its base, and its larger
  factors add resampling variance (about three times g_mom's at damping 1 on synthetic counts). Per-seed differences
  hold that noise, and a small part of any gain may be added spread rather than size.
- 'loses' is a sign on these cuts, not a size; 'not separated' includes steady gains below 0.3.
- Run in the DGX checkout where T42 and T44 ran; a 1e-9 miss on b0, g_mom or m1_r_gmom voids the read uncounted and
  goes to the captain.
- Twice the damping is not twice g_mom. It brings back the entries zeroed at 0.5, which are mostly sampling noise of
  low-mean genes, and grows the mean factor about fourfold over all entries; the clip does not bind. A larger arm's
  result is the result of that combined change and is not read as a size of the development g_mom measured.
- One E8.5 -> E9.5 fit serves every seed and both cuts. The seeds measure draw and scoring noise only, n = 1 for the
  fit, and '8 of 10 ahead' is not replication.
- The two cuts share their truth stage, reference and fit and are not independent evidence. Their truth differs from
  the reference in kind and depth (T27), and a gain from a larger spread may come from that difference rather than
  from development.
- The order rule only makes the choice coarser. Any damping it names is chosen on these cuts and is a tuned constant,
  not a pre-registered one.
- Form B's sign is in-sample and mixes g_mom's size with M1's step, which is itself too small on these cuts (T44). The
  veto can only hold back a form A result; it can never support one.
- No outcome of T49 sets a damping at the test pair. There the fit crosses the change of sample kind and depth (T27,
  N21) that the E8.5 -> E9.5 fit does not, its leading direction does not carry from the backtest's fit (N21: |cos|
  0.20 and 0.29), and a larger damping would enlarge a library contrast. A pass reads the size of the E8.5 -> E9.5 fit
  against an E10.5 truth only.
- No value of T49's statement changes g_mom's damping in any entry or at test. A pass is a tuned in-sample size on
  cuts where 40% of g_mom's gain needs no development (T23). A veto or 'loses' is a sign on these cuts. 'Not
  separated' includes steady gains below 0.3. Each value goes to the captain with these sentences, and any change is
  the captain's ruling.
- A gain at a larger damping fits three readings the read cannot tell apart: a faster pace on E9.5 -> E10.5 than on
  E8.5 -> E9.5, a scorer that rewards state-aligned spread or count noise, and the change of sample kind in the truth.

Usage: on the DGX campaign checkout, repo root, both cuts at one commit, committing nothing until both task files
exist, and on the machine T11 and T42 ran on, so that b0, g_mom and m1_r_gmom reproduce. The environment variables are
docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_size_read.py run --out t49_gmom_size \
      --task "T1/E9.5->E10.5:chamber"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_size_read.py run --out t49_gmom_size \
      --task "T1/E9.5->E10.5"
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/gmom_size_read.py read --out t49_gmom_size
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import MAX_LOSS, SEEDS
from vec.boards import GROUPS, REPO_ROOT, load_boards
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
T11_TASKS = REPO_ROOT / "experiments" / "harness" / "results" / "t11_g_mom" / "tasks"
T42_TASKS = RESULTS / "t42_m1_eb" / "tasks"
sys.path.insert(0, str(HERE))
import eb_read as er  # noqa: E402
import feng_dirs  # noqa: E402
import gmom_combo as gc  # noqa: E402

gv, hygiene, n12 = gc.gv, gc.hygiene, er.n12
CUT_FILES = gc.CUT_FILES
ROW = "T49"
BASE_DAMPING = harness.B0_DAMPING  # 0.5, g_mom's
FORM_A = {"g_mom": 1.0, "g_mom_d075": 1.5, "g_mom_d1": 2.0}  # arm -> multiple of B0_DAMPING
FORM_B = {"m1_r_gmom": 1.0, "m1_r_gmom_g075": 1.5, "m1_r_gmom_g1": 2.0}
DAMPINGS = {a: BASE_DAMPING * m for a, m in {**FORM_A, **FORM_B}.items()}  # 0.5, 0.75, 1.0
LEVELS = tuple(sorted(set(DAMPINGS.values())))
BASE_A, BASE_B = "g_mom", "m1_r_gmom"
LARGER_A, LARGER_B = ("g_mom_d075", "g_mom_d1"), ("m1_r_gmom_g075", "m1_r_gmom_g1")
TWIN = dict(zip(LARGER_A, LARGER_B, strict=True))  # form A arm -> its form B twin
ARMS = ("b0", *FORM_A, *FORM_B)
COMPARISONS_A = {"A1": (BASE_A, "g_mom_d075"), "A2": (BASE_A, "g_mom_d1"), "A3": ("g_mom_d075", "g_mom_d1")}
COMPARISONS_B = {"B1": (BASE_B, "m1_r_gmom_g075"), "B2": (BASE_B, "m1_r_gmom_g1"),
                 "B3": ("m1_r_gmom_g075", "m1_r_gmom_g1")}  # (base, change)
REPRO_TOL = 1e-9  # b0 and g_mom against T11's committed scores, m1_r_gmom against T42's: same machine and code path
REPRO = {"b0": "t11", "g_mom": "t11", "m1_r_gmom": "t42"}
DAMPING_TOL = 1e-12
LOSES_MAX_AHEAD = 2
VETO_LOSS = MAX_LOSS  # form B's veto: the twin's mean loss beyond the adoption rule's 0.3 on either cut
PASSES = ("a larger g_mom scores above g_mom by the adoption rule's arithmetic on the cuts that kept g_mom (a size "
          "chosen on these cuts; not a keep)")
VETOED = ("above g_mom by the adoption rule's arithmetic in try T11's form and loses more than 0.3 in the entry's "
          "form: not put forward")
LOSES = "a larger g_mom loses"
NOT_SEPARATED = "not separated"
POINTS = {m: 100.0 * w * mw for _, w, ms in GROUPS["T1"] for m, mw in ms}  # board points per unit of a metric's skill
LABEL_T49 = ("try T49: the size of g_mom's within-type factors (dampings 0.75 and 1.0 against B0's 0.5) on the two "
             "counted E10.5 cuts; in both forms a labelled in-sample diagnostic: form A selects a constant on the cuts "
             "that kept g_mom, form B (M1 with R1, R2 and g_mom's factors) is the entry's form, whose one use is the "
             "veto; Feng's E9.5 cells lie inside these pseudo-tasks' external windows")
LIMITS = tuple(" ".join(x.split()) for x in __doc__.split("\nLimits")[1].split("\n\n")[0].split("\n- ")[1:])


# ---------------------------------------------------------------------------------------------------------------------
# The fold's fit and its three finishes

_FOLD: dict = {}


@contextmanager
def _capture_finish(captured: list) -> Iterator[None]:
    real = ts.finish_momentum

    def spy(fit, damping):
        captured.append((fit, damping))
        return real(fit, damping)

    ts.finish_momentum = spy
    try:
        yield
    finally:
        ts.finish_momentum = real


def same_factors(a: ts.StateFactors, b: ts.StateFactors) -> bool:
    """Bit for bit: the same genes, categories, metacell indexes, sizes and log factors."""
    if a.genes != b.genes or list(a.categories) != list(b.categories):
        return False
    return all(a.categories[c].index == b.categories[c].index
               and np.array_equal(a.categories[c].sizes, b.categories[c].sizes)
               and a.categories[c].logf.dtype == b.categories[c].logf.dtype
               and a.categories[c].logf.tobytes() == b.categories[c].logf.tobytes() for c in a.categories)


def fold_factors(ctx: harness.BuildContext) -> dict:
    """The fold's raw fit and its finishes at LEVELS, made once per pair of loaded inputs (the harness hands every seed
    the same objects). The fit is the very object `harness._g_mom_factors` passes to `finish_momentum`: its cache is
    cleared, it is called on classified shallow copies of the inputs (as `gmom_combo.m1_r_arms` hands them to it) with
    `finish_momentum` watched, and its own finished factors stand at 0.5. 0.75 and 1.0 are `finish_momentum` of that fit
    at those dampings. Raises unless exactly one finish at 0.5 was made and a fresh finish at 0.5 equals it bit for
    bit."""
    objs = (ctx.inputs[ctx.task.reference_role], ctx.inputs["start"])
    hit = _FOLD.get("fold")
    if hit is not None and all(a is b for a, b in zip(hit["inputs"], objs, strict=True)):
        return hit
    role = ctx.task.reference_role
    ref, src = harness._shallow(objs[0]), ctx.task.role_stems[role]
    start, start_src = harness._shallow(objs[1]), ctx.task.role_stems["start"]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    t1_scope.classify(start, source=start_src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    harness._G_MOM_PREP.clear()
    captured: list = []
    with _capture_finish(captured):
        own = harness._g_mom_factors(ctx, ref, src, start, start_src, shares)
    if len(captured) != 1 or captured[0][1] != BASE_DAMPING:
        raise RuntimeError(f"_g_mom_factors finished {len(captured)} times (dampings {[d for _, d in captured]}), "
                           f"not once at {BASE_DAMPING:g}")
    raw = captured[0][0]
    if not same_factors(ts.finish_momentum(raw, BASE_DAMPING), own):
        raise RuntimeError("a fresh finish of the captured fit at 0.5 is not _g_mom_factors' own factors")
    finished = {d: own if d == BASE_DAMPING else ts.finish_momentum(raw, d) for d in LEVELS}
    _FOLD.clear()
    hit = _FOLD["fold"] = {"inputs": objs, "raw": raw, "finished": finished}
    return hit


def _base_factors(fold: Mapping, got: ts.StateFactors) -> None:
    if got is not fold["finished"][BASE_DAMPING]:
        raise RuntimeError("the arm's build reached factors other than the fold's own finish at 0.5 (another fit?)")


def factor_records(finished: Mapping[float, ts.StateFactors], held_flat: Sequence[str]) -> dict:
    """Per damping and category, over the (metacell, gene) entries of the free genes (the fitted genes not in
    held_flat, which are held at factor 1 when applied; genes off the board are kept): the share not zeroed, the share
    at the clip (|log f| at or above CLIP_AT), the mean absolute log factor over the entries not zeroed and over all
    entries (zeros included), and the share of entries zeroed at BASE_DAMPING and not zeroed at this damping (0 at
    BASE_DAMPING itself). finished must hold BASE_DAMPING. Target-free; no count."""
    base = finished[BASE_DAMPING]
    out = {}
    for d, sf in sorted(finished.items()):
        if sf.genes != base.genes:
            raise ValueError("the finishes are on different gene orders")
        free = ~np.isin(np.asarray(sf.genes, dtype=object), np.asarray(list(held_flat), dtype=object))
        cats = {}
        for c, cf in sf.categories.items():
            a = np.abs(cf.logf[:, free])
            a0 = np.abs(base.categories[c].logf[:, free])
            nz = a > 0
            cats[c] = {"share_not_zeroed": float(nz.mean()) if a.size else 0.0,
                       "share_at_clip": float((a >= ts.CLIP_AT).mean()) if a.size else 0.0,
                       "mean_abs_log_not_zeroed": float(a[nz].mean()) if nz.any() else 0.0,
                       "mean_abs_log_all_entries": float(a.mean()) if a.size else 0.0,
                       "share_unzeroed_vs_0.5": float(((a0 == 0) & nz).mean()) if a.size else 0.0}
        out[f"{d:g}"] = cats
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def form_a_arms(ctx: harness.BuildContext, arms: Mapping[str, float] = FORM_A) -> dict[str, Stage]:
    """Try T11's form: `harness._g_mom_arm` with a transform that returns the fold's factors finished at the arm's
    damping (g_mom's own object at 0.5), after checking it was handed g_mom's own factors."""
    fold = fold_factors(ctx)

    def at(d: float) -> Callable[[ts.StateFactors], ts.StateFactors]:
        def transform(got: ts.StateFactors) -> ts.StateFactors:
            _base_factors(fold, got)
            return fold["finished"][d]
        return transform

    return {a: harness._g_mom_arm(ctx, at(BASE_DAMPING * m)) for a, m in arms.items()}


@contextmanager
def hook_factors(fold: Mapping, damping: float) -> Iterator[None]:
    """`harness._g_mom_factors` as `gmom_combo.gmom_hook` reaches it, returning the fold's factors finished at
    `damping` after the real call (whose cache hands back g_mom's own factors, checked)."""
    real = harness._g_mom_factors

    def at(ctx, ref, src, start, start_src, shares):
        _base_factors(fold, real(ctx, ref, src, start, start_src, shares))
        return fold["finished"][damping]

    harness._g_mom_factors = at
    try:
        yield
    finally:
        harness._g_mom_factors = real


def form_b_arms(ctx: harness.BuildContext, directions: Mapping, groups: Mapping, fraction: float,
                held_m1: Sequence[str], held_factors: Sequence[str],
                arms: Mapping[str, float] = FORM_B) -> dict[str, Stage]:
    """The entry's form: try T42's m1_r_gmom (`gmom_combo.m1_r_arms`, every argument the same object) with the factors
    its hook applies finished at the arm's damping; m1_r, built by the same call, is discarded."""
    fold = fold_factors(ctx)
    out = {}
    for a, m in arms.items():
        with hook_factors(fold, BASE_DAMPING * m):
            out[a] = gc.m1_r_arms(ctx, directions, groups, fraction, held_m1, held_factors)[BASE_B]
    return out


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T49 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T49 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
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
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    _FOLD.clear()
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = {"b0": b0(ctx), **form_a_arms(ctx),
                 **form_b_arms(ctx, directions, groups, fraction, held_m1, preregistered)}
        arms = er.score_arms(board, target, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms, "seconds": round(time.time() - t0, 1)})
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res = {"label": LABEL_T49, "task": t.key, "dampings": dict(DAMPINGS), "reproduction_tolerance": REPRO_TOL,
           "fraction": fraction, "external": ext, "seeds": list(seeds),
           "factor_records": factor_records(_FOLD["fold"]["finished"], preregistered), "limits": list(LIMITS),
           "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, REPO_ROOT / union["file"],
                                       feng_dirs.clean_stats_path(), n22])}
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _committed_score(kind: str, row: Mapping, arm: str) -> float | None:
    """A committed row's score of `arm`: T11's harness rows hold `scores`, T42's rows `arms`."""
    if kind == "t11":
        return row.get("scores", {}).get(arm)
    return row.get("arms", {}).get(arm, {}).get("score")


def reproduction(rows: Sequence[Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    """Per seed, b0 and g_mom against T11's committed scores and m1_r_gmom against T42's; holds when every seed and
    arm is there and each difference is within REPRO_TOL. committed: "t11" / "t42" -> that cut's committed rows."""
    by = {k: {r["seed"]: r for r in v} for k, v in committed.items()}
    worst = {a: 0.0 for a in REPRO}
    missing, absent = [], []
    for r in rows:
        for arm, kind in REPRO.items():
            old_row = by.get(kind, {}).get(r["seed"])
            if old_row is None:
                missing.append(f"{kind}:{r['seed']}")
                continue
            old = _committed_score(kind, old_row, arm)
            if old is None or arm not in r["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(float(r["arms"][arm]["score"]) - float(old)))
    holds = (all(committed.get(k) for k in set(REPRO.values())) and not missing and not absent
             and all(w <= REPRO_TOL for w in worst.values()))
    return {"largest_difference": worst, "missing_seeds": sorted(set(missing)), "absent": absent, "holds": holds,
            "tolerance": REPRO_TOL}


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Sequence[Mapping]]]) -> dict:
    """The void conditions. results: cut -> task file (a missing cut is absent); committed: cut -> {"t11": rows,
    "t42": rows} (a missing file is absent)."""
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
        got = r.get("dampings") or {}
        bad = [a for a, d in DAMPINGS.items()
               if not isinstance(got.get(a), int | float) or abs(got[a] - d) > DAMPING_TOL]
        if bad:
            reasons.append(f"{k}: dampings {got}, not {DAMPINGS} (wrong damping on {bad})")
    commits = {r["stamp"]["commit"] for r in results.values()}
    unclean = sorted(k for k, r in results.items() if r["stamp"]["tree_clean"] is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the cuts ran at {len(commits)} commits")
    repro = {}
    for k, r in results.items():
        repro[k] = reproduction(r["rows"], committed.get(k, {}))
        if not repro[k]["holds"]:
            reasons.append(f"{k}: b0 and g_mom (T11's) and m1_r_gmom (T42's) do not reproduce their committed scores "
                           f"within {REPRO_TOL:g} per seed")
    complete = all(a in row["arms"] for r in results.values() for row in r["rows"] for a in ARMS)
    if results and complete:
        same = [f"{c} = {b}" for b, cs in ((BASE_A, LARGER_A), (BASE_B, LARGER_B)) for c in cs
                if all(row["arms"][c] == row["arms"][b] for r in results.values() for row in r["rows"])]
        if same:
            reasons.append(f"larger arms score exactly as their base on every seed (built at 0.5?): {same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro, "commits": sorted(c for c in commits if c)}


def compare(results: Mapping[str, Mapping], base: str, change: str) -> dict:
    """eb_read.compare (the adoption rule, paired by seed, both cuts; the mean paired differences of skills, groups and
    metrics) with each cut's per-seed range and each metric's part of the difference in board points."""
    out = er.compare(results, base, change)
    for k, r in results.items():
        d = [row["arms"][change]["score"] - row["arms"][base]["score"] for row in r["rows"]]
        cut = out["cuts"][k]
        cut["range"] = [float(min(d)), float(max(d))]
        cut["points"] = {m: POINTS[m] * x for m, x in cut["skills"].items() if m in POINTS}
    return out


def _loses(comps: Mapping[str, Mapping]) -> bool:
    return all(c["mean"] < 0 and c["ahead"] <= LOSES_MAX_AHEAD
               for name in ("A1", "A2") for c in comps[name]["cuts"].values())


def statement(comps_a: Mapping[str, Mapping], comps_b: Mapping[str, Mapping]) -> dict:
    """The statement from form A's three comparisons, with form B's veto on the arm put forward."""
    p = {c: comps_a[c]["passes"] for c in COMPARISONS_A}
    if p["A2"] and p["A3"]:
        arm, basis = "g_mom_d1", "passes against g_mom and against g_mom_d075"
    elif p["A1"]:
        arm, basis = "g_mom_d075", "passes against g_mom"
    elif p["A2"]:
        arm, basis = "g_mom_d1", "alone passes against g_mom"
    else:
        arm, basis = None, None
    if arm is None:
        means = {n: [comps_a[n]["cuts"][k]["mean"] for k in CUT_FILES if k in comps_a[n]["cuts"]] for n in ("A1", "A2")}
        if _loses(comps_a):
            parts = [f"{n} " + " / ".join(f"{x:+.2f}" for x in m) for n, m in means.items()]
            text = f"{LOSES} (" + ", ".join(parts) + ")"
        else:
            text = f"{NOT_SEPARATED} (largest mean gain {max(x for m in means.values() for x in m):+.2f})"
        return {"statement": text, "arm": None, "basis": None, "veto": None, "put_forward": False}
    twin = TWIN[arm]
    twin_cmp = comps_b["B1" if twin == "m1_r_gmom_g075" else "B2"]
    loses = any(c["mean"] < -VETO_LOSS for c in twin_cmp["cuts"].values())
    veto = {"twin": twin, "means": {k: c["mean"] for k, c in twin_cmp["cuts"].items()}, "loses": loses}
    text = f"{VETOED} ({arm})" if loses else f"{PASSES}: {arm}"
    return {"statement": text, "arm": arm, "basis": basis, "veto": veto, "put_forward": not loses}


def t49_reading(results: Mapping[str, Mapping], committed: Mapping[str, Mapping[str, Sequence[Mapping]]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    comps_a = {c: compare(results, *pair) for c, pair in COMPARISONS_A.items()}
    comps_b = {c: compare(results, *pair) for c, pair in COMPARISONS_B.items()}
    st = statement(comps_a, comps_b)
    out.update({"outcome": st["statement"], "statement": st, "form_a": comps_a, "form_b_diagnostic": comps_b})
    out["mean_scores"] = {k: {a: float(np.mean([row["arms"][a]["score"] for row in r["rows"]])) for a in ARMS}
                          for k, r in results.items()}
    out["factor_records"] = {k: r.get("factor_records") for k, r in results.items()}
    return out


def t49_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"Statement: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, g in v["reproduction"].items():
        w = g["largest_difference"]
        lines.append(f"{k}: largest per-seed difference from the committed scores: b0 {w['b0']:.3g}, g_mom "
                     f"{w['g_mom']:.3g} (T11's), m1_r_gmom {w['m1_r_gmom']:.3g} (T42's) (gate {REPRO_TOL:g})")
    st = r.get("statement")
    if st and st["arm"]:
        vt = st["veto"]
        lines.append(f"Arm: {st['arm']} ({st['basis']}); its twin {vt['twin']} against m1_r_gmom, mean "
                     + ", ".join(f"{k} {x:+.2f}" for k, x in vt["means"].items())
                     + (f": loses more than {VETO_LOSS:g} (veto)" if vt["loses"]
                        else f": does not lose more than {VETO_LOSS:g}"))
    for k, ms in r.get("mean_scores", {}).items():
        lines.append(f"mean scores {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in ms.items()))

    def block(name: str, comp: Mapping) -> None:
        lines.append(f"({name}) {comp['change']} against {comp['base']}: the adoption rule's arithmetic "
                     f"{'passes' if comp['passes'] else 'does not pass'}, not a keep ({comp['rule']['reason']})")
        for k, cut in comp["cuts"].items():
            sk = ", ".join(f"{m} {x:+.4f} ({cut['points'].get(m, float('nan')):+.2f} pts)"
                           for m, x in cut["skills"].items())
            lines.append(f"    {k}: {cut['mean']:+.2f} points, ahead {cut['ahead']} of {cut['of']}, per seed "
                         f"{cut['range'][0]:+.2f} to {cut['range'][1]:+.2f}; skills: {sk}")

    for c, comp in r.get("form_a", {}).items():
        block(c, comp)
    if r.get("form_b_diagnostic"):
        lines.append("Form B, the entry's form: a labelled in-sample diagnostic; its one use is the veto")
        for c, comp in r["form_b_diagnostic"].items():
            block(c, comp)
    for k, rec in r.get("factor_records", {}).items():
        for d, cats in (rec or {}).items():
            lines.append(f"factors {k} damping {d}: " + "; ".join(
                f"{c} not zeroed {x['share_not_zeroed']:.3f}, at clip {x['share_at_clip']:.4f}, mean |log f| "
                f"{x['mean_abs_log_not_zeroed']:.4f} not zeroed, {x['mean_abs_log_all_entries']:.4f} all entries, "
                f"zeroed at 0.5 and not here {x['share_unzeroed_vs_0.5']:.3f}" for c, x in cats.items()))
    lines.append("Limits (written before the run; the result goes to the team captain):")
    lines.extend(f"- {x}" for x in LIMITS)
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t11_tasks: Path = T11_TASKS, t42_tasks: Path = T42_TASKS) -> dict:
    """Read the cuts present in out/tasks against T11's and T42's committed task files; write t49.{json,txt}."""
    gv.refuse_existing(out, ("t49.json", "t49.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / f for k, f in CUT_FILES.items() if (out / "tasks" / f).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    paths, committed = [], {}
    for k, f in CUT_FILES.items():
        committed[k] = {}
        for kind, folder in (("t11", t11_tasks), ("t42", t42_tasks)):
            if (folder / f).exists():
                paths.append(folder / f)
                committed[k][kind] = json.loads((folder / f).read_text())["rows"]
    res = {"label": LABEL_T49, "reading": t49_reading(results, committed), "limits": list(LIMITS),
           "stamps": {k: v["stamp"] for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *paths])
    gv.write_outputs(out, "t49", res, t49_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T49's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one cut of try T49")
    p.add_argument("--task", required=True, choices=list(CUT_FILES))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T49's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
