#!/usr/bin/env python
"""Try T54 (experiments/harness/TRIES.md): M1's step when the truth is deeper than the cells the prediction is built
from. A labelled in-sample read on try T50's route (experiments/t1_m1/over_read.py), with the reference thinned instead
of the truth.

Try T50 found that the multiple of M1's step with the highest score moves up as the truth gets shallower against the
cells the prediction is built from. Every released truth of the E10.5 cuts is shallower than its reference, so the other
side, a truth deeper than the prediction's cells, has never been read. This read makes the reference shallower.

  run    One world of the chamber cut T1/E9.5->E10.5:chamber, 10 seeds, with the truth `lib3` exactly as try T46's run
         makes it (unthinned; over_read.truth_and_contamination). Worlds:
           ref100   the cut's inputs as they are (the gate: try T50's route bit for bit);
           ref67    every cell of the reference stage (E9.5) has every integer count c replaced by a draw from
                    Binomial(c, 2/3), one generator on the stream [54, 1];
           ref67b   the same on the stream [54, 3] (a second draw: statement 3's self-check);
           ref50    Binomial(c, 1/2) on the stream [54, 2].
         The thinning is try T50's own (over_read.thin_counts: cells in the pool's order, within a cell the file's gene
         order), on the integer counts of the reference pool as the route loads it (vec.prep.load_stage's counts layer,
         recovered from the file's .X by vec.prep.recover_counts and read by vec.core.counts_of), before anything else
         is made; it does not depend on the build seed: one thinned reference per world for all ten seeds. The
         thinned counts reach every form of the reference the route reads:
           - the input pool (ctx.inputs["ref"], the reference role): counts layer = the thinned counts (int32, as
             load_stage holds them), .X = vec.prep.log1p_cp10k of them over every gene of the file, obs, var and obsm
             the pool's. Through it: B0's draw, its chamber shares (labels only) and its momentum shift (start ->
             reference type means); every M1 arm's draw, Feng step (external_shift) and g_mom hook; g_mom's fit
             (harness._g_mom_factors, whose module cache is keyed by the input objects' identity and so refits on the
             thinned pool). The start stage (E8.5, ctx.inputs["start"]) is the released one.
           - the scorer's reference stage (Stage.from_anndata on the panel): over_read.as_scored of the thinned counts
             (log1p CP10k over every gene of the file, then the panel), as try T50 turns thinned counts into what the
             scorer reads. Through it: the floor's anchors (the floor is a subsample of it, and it is the de terms'
             reference) and every arm's scoring (the de terms' reference); the board's limits read its cell count
             only, which thinning keeps.
         Before thinning, the run checks that the recovered counts give back both forms as the route reads them (the
         pool's .X through log1p_cp10k, and the scorer's matrix through as_scored; largest absolute difference at most
         over_read.REBUILD_TOL, both recorded), that the full file's cells and genes are the pool's in the same order,
         and refuses otherwise. Arms per seed, all from one cell draw: b0, m1_r_gmom (fraction 0.2) and
         step_read.step_arms at 1.5, 2, 3, 4, 5 and 6 times the fraction (m1_r_gmom_x15, _x2, _x3, _x4, _x5, _x6), built
         as try T50 builds them. The task file holds what try T50's holds (every arm's fraction, the truth's name, its
         row and contamination counts as checks compared by code and never printed, per seed and arm the prediction's
         sha256, the scored genes' zero share and detected scored genes per cell, per seed the anchors), the world's
         name and, for a thinned world, p, the stream, the realised depth ratio (thinned counts over released counts: a
         ratio), the rebuild checks' largest differences and whether any cell was left with no count. Refused on an
         existing task file, unless the committed N22 passed, and unless experiments/harness/TRIES.md holds row T54.
  read   T54's reading of the four worlds' files (rule below, fixed in the row before the run). Writes t54.json and
         t54.txt and refuses to overwrite them.

The reading. Void unless the four worlds hold the ten seeds with every arm and its fraction at one commit on a clean
tree (tree_clean exactly true), each file's task the chamber cut, its truth lib3 and its world its own name; no file
records an anchor failure; in ref100 nothing is thinned and every arm equals, seed for seed, its committed score in
experiments/t1_m1/results/t50_m1_over_step/tasks/lib3.json to 1e-9; each thinned world records its own p and stream,
its realised depth ratio is within 0.01 of its p, no cell is emptied and its rebuild checks hold; in no thinned world
does an arm score exactly as in ref100 on every seed.
The test on one world's paired differences is try T46's single-task test: "ahead" is a mean paired gain above 0.3 with
at least 8 of 10 seeds ahead; "behind" the same reversed; otherwise "level".
Statement 1, the turn, per world: the multiplier with the highest mean score among 1, 1.5, 2, 3, 4, 5 and 6 (the
smallest on a tie), and the labels of x4, x5 and x6 against x3, of x3 against x2 and of x2 against x1.
Statement 2, the cases a request would quote, per world: x4 against x3, x5 against x3, x4 against x2 and x5 against x2,
each with its mean paired difference, seeds ahead, per-seed range and label.
Statement 3, depth, deep side. Per seed D67 = (x4 - x2 in ref67) - (x4 - x2 in ref100), and D50 the same with ref50.
"The larger steps gain less when the truth is deeper than the reference" if D67 and D50 are both behind; "more" if both
are ahead; "no difference at the bar" if both are level; otherwise "mixed" with both labels. Statement 3 is void, and
only reported, unless the draw null, (x4 - x2 in ref67b) - (x4 - x2 in ref67), reads level. Beside it, read by nothing:
the same differences on each scored metric's raw value, each world's floor-to-ceiling gap per metric, and the curve's
shape per world side by side.
Reported: per world every arm's mean score, neighbour gains per 0.1 of fraction, the parts on the two de terms and on
mmd_u and the variogram for statement 2's comparisons, each arm's scored zero share and detected genes per cell. The
read puts nothing forward by itself.

NOT an adoption result: Feng's E9.5 cells lie inside this pseudo-task's external window, the cut is one earlier T1
choices were made on, and the truth is a part of a released target read in-sample. The result can only go to the team
captain, and no outcome is evidence on Task 1's read depth, a lever that stays closed.

Usage: on the DGX campaign checkout, repo root, all four worlds at one commit, committing nothing until the four task
files exist, and on the machine T50 ran on, so that ref100 reproduces. The environment variables are docs/COMPUTE.md's
for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/deep_read.py run --world ref100 --out t54_m1_deep_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/deep_read.py run --world ref67 --out t54_m1_deep_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/deep_read.py run --world ref67b --out t54_m1_deep_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/deep_read.py run --world ref50 --out t54_m1_deep_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/deep_read.py read --out t54_m1_deep_step
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from vec import core, harness, prep, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage
from vec.io import read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import over_read as ob  # noqa: E402

ls, sr, gc, er = ob.ls, ob.sr, ob.gc, ob.er
gv, hygiene, n12, feng_dirs = ob.gv, ob.hygiene, ob.n12, ob.feng_dirs
TASK = ob.TASK
HEART = ob.HEART  # lib3, unthinned, the truth of every world
T50_LIB3 = RESULTS / "t50_m1_over_step" / "tasks" / "lib3.json"
GATE = "ref100"
THINNED = {"ref67": {"p": 2.0 / 3.0, "stream": [54, 1]}, "ref50": {"p": 1.0 / 2.0, "stream": [54, 2]},
           "ref67b": {"p": 2.0 / 3.0, "stream": [54, 3]}}
WORLDS = (GATE, "ref67", "ref67b", "ref50")
DRAW_NULL = ("ref67b", "ref67")  # a second draw at p = 2/3 against the first: statement 3's self-check
BASE = sr.BASE  # m1_r_gmom, fraction 0.2
X15, X2, X3, X4, X5, X6 = "m1_r_gmom_x15", "m1_r_gmom_x2", ob.X3, ob.X4, ob.X5, ob.X6
MULTIPLIERS = {a: m for a, m in ob.MULTIPLIERS.items() if a != ob.X8}  # x1.5 to x6, T50's names and builds
ARMS = ("b0", BASE, *MULTIPLIERS)
ROW = "T54"
FRACTION = sr.FRACTION  # 0.2
FRACTION_TOL = sr.FRACTION_TOL
REPRO_TOL = sr.REPRO_TOL  # 1e-9
DEPTH_TOL = ob.DEPTH_TOL  # 0.01
REBUILD_TOL = ob.REBUILD_TOL  # 1e-5
TURN_AGAINST_X3 = (X4, X5, X6)
QUOTED = ((X4, X3), (X5, X3), (X4, X2), (X5, X2))  # (change, base), statement 2
STILL_RISING = "still rising at 6"
MORE = "The larger steps gain more when the truth is deeper than the reference"
LESS = "The larger steps gain less when the truth is deeper than the reference"
NO_DIFFERENCE = "no difference at the bar"
MIXED = "mixed"
DEPTH_VOID = ("void: the second thinning draw (ref67b) against the first (ref67) is not level by the single-task test, "
              "so one draw alone can pass the bar; statement 3 is only reported")
TASK_FILE = "{world}.json"
LABEL_T54 = ("try T54, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng step multiplied by 1, 1.5, 2, 3, 4, 5 and 6, on T50's route, against "
             "library 3's truth of the chamber cut as released, with the reference stage (E9.5) as released and "
             "thinned to two thirds (two draws) and a half of its depth; Feng's E9.5 cells lie inside this "
             "pseudo-task's external window and the cut is one earlier T1 choices were made on; the result can only go "
             "to the team captain, and no outcome is evidence on Task 1's read depth, which stays closed")
LIMITS = (
    "A thinned reference is not a deeper truth: thinning lowers the reference's depth and raises its count noise "
    "together, changes the floor (the copy of the reference) and with it every skill's scale, and changes g_mom's fit, "
    "so board points are not comparable across worlds and statement 3 is read with the raw metrics and the shapes.",
    "Thinning the reference also changes B0's momentum shift and the M1 arms' momentum-free cells (they are drawn from "
    "the thinned pool), and the de terms' reference in the scorer, which is the same reference stage; the start stage "
    "(E8.5) is not thinned, so g_mom's fit and B0's momentum compare a released start with a thinned reference.",
    "One released day, in-sample, the chamber cut and one library only.",
    "One thinning draw per world, with a single second draw as self-check; a level null does not certify the draw.",
    "Nothing here says how deep any hidden truth is or how deep it is against the cells an upload is built from, and "
    "no statement may be used with a figure for either.",
    "No outcome is evidence on Task 1's read depth, which stays closed.",
    "Multiples of 5 and 6 apply Feng's whole five-day change or more to a one-day truth; they are read to find where "
    "the score turns, not as candidates. No multiplier read here is a candidate setting; the test's step stays the "
    "bracket fraction by rule.",
    "Where adjacent multipliers differ by less than their per-seed spread the argmax is not separated; the argmax "
    "alone is not a finding.",
    "A label of 'level' in statement 2 means the difference is not separated at a 0.3 mean with 8 of 10 seeds on one "
    "in-sample day; 'level' is not 'no loss'.",
    "The depth-ratio check catches a wrong p or a wrong stream only.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The thinned reference


@contextmanager
def _t50_name(world: str) -> Iterator[None]:
    """over_read.thin_counts reads p and the stream from over_read.THINNED by name: the world's spec under its own
    name for the call, restored after it."""
    missing = object()
    old = ob.THINNED.get(world, missing)
    ob.THINNED[world] = THINNED[world]
    try:
        yield
    finally:
        if old is missing:
            ob.THINNED.pop(world, None)
        else:
            ob.THINNED[world] = old


def thin_counts(C, world: str) -> sparse.csr_matrix:
    """Try T50's thinning (over_read.thin_counts) at the world's p on the world's stream: every integer count c a draw
    from Binomial(c, p), cells in row order and within a cell the columns' order. No build seed enters."""
    with _t50_name(world):
        return ob.thin_counts(C, world)


def _largest_difference(a, b) -> float:
    """The largest absolute entry of a - b (in float64; two sparse matrices are compared without densifying)."""
    if sparse.issparse(a) and sparse.issparse(b):
        d = sparse.csr_matrix(a, dtype=np.float64) - sparse.csr_matrix(b, dtype=np.float64)
        return float(np.abs(d.data).max()) if d.nnz else 0.0
    a = a.toarray() if sparse.issparse(a) else np.asarray(a)
    b = b.toarray() if sparse.issparse(b) else np.asarray(b)
    return float(np.abs(a.astype(np.float64) - b.astype(np.float64)).max()) if a.size else 0.0


def thinned_world(world: str, inputs: Mapping[str, ad.AnnData], role: str, full: ad.AnnData, ref_stage: Stage,
                  genes: Sequence[str]) -> tuple[dict[str, ad.AnnData], Stage, dict]:
    """The world's inputs (the reference role's pool thinned, every other role as loaded) and the scorer's reference
    stage, both from one thinning of the pool's integer counts; and the record. full: the reference file as read (every
    gene, its own .X); ref_stage: the scorer's reference stage as the route reads it. Refused (ValueError) when the
    full file's cells, labels or genes are not the pool's in the same order, or when the recovered counts do not give
    back the pool's .X or the scorer's matrix within REBUILD_TOL."""
    pool = inputs[role]
    if full.n_obs != pool.n_obs or list(full.var_names) != list(pool.var_names):
        raise ValueError(f"{world}: the reference file's cells or genes are not the input pool's")
    obs = full.obs[[c for c in ("celltype",) if c in full.obs]].copy()
    if "celltype" in obs and not np.array_equal(np.asarray(obs["celltype"]).astype(str),
                                                np.asarray(pool.obs["celltype"]).astype(str)):
        raise ValueError(f"{world}: the reference file's cell types are not the input pool's in its order")
    C = core.counts_of(pool)
    var_names = list(full.var_names)
    pool_diff = _largest_difference(prep.log1p_cp10k(C), pool.X)
    rebuilt = ob.as_scored(C, obs, var_names, genes)
    scored_diff = (_largest_difference(rebuilt.X, ref_stage.X) if rebuilt.X.shape == ref_stage.X.shape
                   else float("inf"))
    if not (pool_diff <= REBUILD_TOL and scored_diff <= REBUILD_TOL):
        raise ValueError(f"{world}: the counts recovered from the reference do not give back the matrices the route "
                         f"reads (the pool's .X: largest difference {pool_diff:.3g}; the scorer's: {scored_diff:.3g}; "
                         f"bound {REBUILD_TOL:g}); T54 does not thin normalised values")
    if ref_stage.celltype is not None and not np.array_equal(rebuilt.celltype, ref_stage.celltype):
        raise ValueError(f"{world}: the rebuilt cells' types differ from the scorer's reference stage's")
    thinned = thin_counts(C, world)
    new_pool = ad.AnnData(X=prep.log1p_cp10k(thinned), obs=pool.obs.copy(), var=pool.var.copy(),
                          obsm=dict(pool.obsm), uns=dict(pool.uns))
    new_pool.layers["counts"] = sparse.csr_matrix(thinned, dtype=np.int32)
    out = {r: (new_pool if a is pool else a) for r, a in inputs.items()}
    before = float(C.sum())
    rec = {"world": world, "p": THINNED[world]["p"], "stream": list(THINNED[world]["stream"]),
           "of": "the reference stage (E9.5), every cell and gene of the file",
           "depth_ratio": float(thinned.sum()) / before if before > 0 else float("nan"),
           "rebuild_largest_difference": max(pool_diff, scored_diff),
           "rebuild_pool_largest_difference": pool_diff, "rebuild_scored_largest_difference": scored_diff,
           "any_cell_emptied": bool(((np.asarray(thinned.sum(1)).ravel() == 0)
                                     & (np.asarray(C.sum(1)).ravel() > 0)).any())}
    return out, ob.as_scored(thinned, obs, var_names, genes), rec


def world_reference(world: str, inputs: Mapping[str, ad.AnnData], role: str, ref_stage: Stage, genes: Sequence[str],
                    read_full: Callable[[], ad.AnnData]) -> tuple[Mapping[str, ad.AnnData], Stage, dict | None]:
    """ref100: the inputs and the scorer's reference stage as loaded, untouched (read_full is not called); a thinned
    world: thinned_world on the reference file read by read_full."""
    if world == GATE:
        return inputs, ref_stage, None
    if world not in THINNED:
        raise ValueError(f"unknown world {world!r}; one of {list(WORLDS)}")
    return thinned_world(world, inputs, role, read_full(), ref_stage, genes)


# ---------------------------------------------------------------------------------------------------------------------
# The run


def arm_fractions(fraction: float) -> dict[str, float]:
    return sr.arm_fractions(fraction, MULTIPLIERS)


def build_seed(ctx: harness.BuildContext, b0, directions: Mapping, groups: Mapping, fraction: float,
               held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """One seed's arms in the order try T50's run builds them: b0, m1_r_gmom, the scaled arms up to x6."""
    return {"b0": b0(ctx), **sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors,
                                          multipliers=MULTIPLIERS)}


def seed_rows(view, cboard, board, truth: Stage, ref_stage: Stage, inputs: Mapping[str, ad.AnnData], floors: Mapping,
              seeds: Sequence[int], directions: Mapping, groups: Mapping, fraction: float, held_m1: Sequence[str],
              preregistered: Sequence[str], label: str = "") -> list[dict]:
    """Per seed: the arms built from the world's inputs and scored against lib3 with the world's reference stage and
    that seed's anchors, as try T50's run records them."""
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = build_seed(ctx, b0, directions, groups, fraction, held_m1, preregistered)
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "anchors": ob.anchor_record(f.anchors), "arms": arms,
                     "sha256": {k: ls.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: ls.sparsity(p) for k, p in preds.items()}, "seconds": round(time.time() - t0, 1)})
        print(f"{label} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    return rows


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T54 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T54 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
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
    try:  # the thinning comes before anything else is made
        inputs, ref_stage, thinning = world_reference(a.world, inputs, role, ref_stage, genes,
                                                      lambda: read_stage(data / t.reference, None))
    except ValueError as e:
        sys.exit(f"T54 is refused on {a.world}: {e}")
    try:
        truth, contamination, _ = ob.truth_and_contamination(t, data, genes, HEART)
    except ValueError as e:
        sys.exit(f"T54 is refused on {a.world}: {e}")
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T54, "task": t.key, "world": a.world, "truth": HEART, "truth_rows": int(truth.n),
           "truth_contamination": contamination, "seeds": list(seeds), "thinning": thinning}
    floors, failures = ls.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        (out / "tasks").mkdir(parents=True, exist_ok=True)
        (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")
        sys.exit(f"T54 is refused on {a.world}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T54 is void)")
    pool = inputs[role]
    directions, ext = feng_dirs.feng_directions(list(pool.var_names))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    rows = seed_rows(view, cboard, board, truth, ref_stage, inputs, floors, seeds, directions, groups, fraction,
                     held_m1, preregistered, label=a.world)
    res.update({"fraction": fraction, "fractions": arm_fractions(fraction), "multipliers": dict(MULTIPLIERS),
                "external": ext, "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _rows(r: Mapping) -> list:
    return ob._rows(r)


def _named(r: Mapping, world: str) -> dict:
    """A world's file with its world's name as the 'truth' the comparisons print (every world's truth is lib3)."""
    return {**r, "truth": world}


def reproduction(gate: Mapping, committed: Mapping) -> dict:
    """ref100 against try T50's committed lib3 file, per seed and arm; holds when every seed and arm is there and each
    is within REPRO_TOL."""
    worst = {a: 0.0 for a in ARMS}
    absent = []
    old = {r["seed"]: r for r in _rows(committed)}
    for r in _rows(gate):
        o = old.get(r["seed"])
        for arm in ARMS:
            if o is None or arm not in o.get("arms", {}) or arm not in r["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - o["arms"][arm]["score"]))
    holds = bool(_rows(gate)) and bool(old) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "absent": absent, "holds": holds}


def validity(results: Mapping[str, Mapping], committed: Mapping) -> dict:
    """The void conditions. results: world -> task file (a missing world is absent); committed: try T50's lib3 file."""
    reasons = []
    missing = [k for k in WORLDS if k not in results]
    if missing:
        reasons.append(f"missing worlds {missing}")
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
        if r.get("truth") != HEART:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}, not {HEART}")
        if r.get("world") != k:
            reasons.append(f"{k}: the file holds world {r.get('world')!r}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items()
               if not isinstance(got_fr.get(a), int | float) or abs(got_fr[a] - f) > FRACTION_TOL]
        if bad:
            reasons.append(f"{k}: arm fractions {got_fr}, not {want_fr} (wrong fraction on {bad})")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the worlds ran at {len(commits)} commits")
    repro = None
    if GATE in results:
        if results[GATE].get("thinning") is not None:
            reasons.append(f"{GATE}: the file records a thinning")
        repro = reproduction(results[GATE], committed)
        if not repro["holds"]:
            reasons.append(f"{GATE}: {', '.join(ARMS)} do not reproduce try T50's committed lib3 scores within "
                           f"{REPRO_TOL:g} per seed")
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
            reasons.append(f"{k}: the thinning left a reference cell with no count (or did not record whether it "
                           f"did): any_cell_emptied {th.get('any_cell_emptied')!r}")
        rebuild = th.get("rebuild_largest_difference")
        if not isinstance(rebuild, int | float) or not rebuild <= REBUILD_TOL:
            reasons.append(f"{k}: the rebuild check's largest difference {rebuild!r} is not at most {REBUILD_TOL:g}")
        gate_rows = {row["seed"]: row for row in _rows(results.get(GATE, {}))}
        rows = _rows(results[k])
        if gate_rows and rows and all(a in row["arms"] for row in rows for a in ARMS):
            same = [a for a in ARMS
                    if all(row["seed"] in gate_rows and a in gate_rows[row["seed"]]["arms"]
                           and row["arms"][a]["score"] == gate_rows[row["seed"]]["arms"][a]["score"] for row in rows)]
            if same:
                reasons.append(f"{k}: arms score exactly as in {GATE} on every seed (thinning idle?): {same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def comparison(r: Mapping, world: str, change: str, base: str) -> dict:
    """try T50's comparison (T46's single-task test on one world, its sd, the reversed test, the points by part) with
    the label ahead / behind / level."""
    return ob.comparison(r, world, change, base)


def _multiplier(arm: str) -> float:
    return 1.0 if arm == BASE else MULTIPLIERS[arm]


def turn(mean_scores: Mapping[str, float]) -> dict:
    """The multiplier with the highest mean board score among 1 (m1_r_gmom) and the six scaled arms (the smallest on a
    tie); STILL_RISING when that is 6."""
    by_m = {1.0: mean_scores[BASE], **{m: mean_scores[a] for a, m in MULTIPLIERS.items()}}
    best = max(by_m, key=lambda m: (by_m[m], -m))
    return {"by_multiplier": {f"{m:g}": v for m, v in by_m.items()}, "best": best,
            "statement": STILL_RISING if best == max(MULTIPLIERS.values()) else f"turns at {best:g}",
            "tie_rule": "smallest multiplier"}


def neighbour_of(arm: str) -> str:
    """The arm at the next smaller multiplier (m1_r_gmom for x1.5)."""
    ms = [a for a, _ in sorted(MULTIPLIERS.items(), key=lambda x: x[1])]
    i = ms.index(arm)
    return BASE if i == 0 else ms[i - 1]


def neighbours(r: Mapping, world: str) -> dict:
    """Each scaled arm against its next smaller neighbour, with the gain per 0.1 of fraction."""
    out = {}
    for a in MULTIPLIERS:
        b = neighbour_of(a)
        c = comparison(r, world, a, b)
        tenths = (MULTIPLIERS[a] - _multiplier(b)) * FRACTION / 0.1
        out[a] = {"base": b, "mean": c["mean"], "label": c["label"], "per_tenth": c["mean"] / tenths}
    return out


def sparsity_means(r: Mapping) -> dict:
    """Per arm, the mean over seeds of the scored zero share and of the detected scored genes per cell; None where a
    row does not record it."""
    out = {}
    for a in ARMS:
        recs = [(row.get("sparsity") or {}).get(a) for row in _rows(r)]
        if not recs or any(not isinstance(x, Mapping) or not {"zero_share", "detected_per_cell"} <= set(x)
                           for x in recs):
            out[a] = None
            continue
        out[a] = {q: float(np.mean([x[q] for x in recs])) for q in ("zero_share", "detected_per_cell")}
    return out


def statement_one(r: Mapping, world: str, means: Mapping[str, float]) -> dict:
    t = turn(means)
    return {"argmax": t["best"], "statement": t["statement"], "by_multiplier": t["by_multiplier"],
            "against_x3": {a: comparison(r, world, a, X3)["label"] for a in TURN_AGAINST_X3},
            "x3_against_x2": comparison(r, world, X3, X2)["label"],
            "x2_against_x1": comparison(r, world, X2, BASE)["label"]}


def depth_statement(d67: Mapping, d50: Mapping) -> str:
    if d67["behind"] and d50["behind"]:
        return LESS
    if d67["ahead"] and d50["ahead"]:
        return MORE
    if d67["label"] == "level" and d50["label"] == "level":
        return NO_DIFFERENCE
    return f"{MIXED} (D67 {d67['label']}, D50 {d50['label']})"


def statement_three(results: Mapping[str, Mapping], s1: Mapping) -> dict:
    gate = _named(results[GATE], GATE)
    d67 = ob.paired_difference(_named(results["ref67"], "ref67"), gate, X4, X2)
    d50 = ob.paired_difference(_named(results["ref50"], "ref50"), gate, X4, X2)
    b, a = (_named(results[k], k) for k in DRAW_NULL)
    null = ob.paired_difference(b, a, X4, X2)
    holds = null["label"] == "level"
    label = depth_statement(d67, d50)
    return {"valid": holds, "draw_null": {**null, "holds": holds},
            "D67": d67, "D50": d50, "statement": label if holds else DEPTH_VOID, "label_if_valid": label,
            "shapes": {k: s1[k] for k in WORLDS},
            "anchor_gap": {k: ob.anchor_gap(results[k]) for k in WORLDS}}


def outcome_part(s1: Mapping) -> str:
    labels = ", ".join(f"x{a.rsplit('_x', 1)[1]} {lab}" for a, lab in s1["against_x3"].items())
    return (f"against x3 {labels}; x3 against x2 {s1['x3_against_x2']}; x2 against x1 {s1['x2_against_x1']}; "
            f"argmax {s1['argmax']:g}" + (f" ({STILL_RISING})" if s1["statement"] == STILL_RISING else ""))


def t54_reading(results: Mapping[str, Mapping], committed: Mapping) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in results[k]["rows"]])) for a in ARMS}
             for k in WORLDS}
    s1 = {k: statement_one(results[k], k, means[k]) for k in WORLDS}
    out["statement_1"] = s1
    out["statement_2"] = {k: {f"{c} against {b}": comparison(results[k], k, c, b) for c, b in QUOTED} for k in WORLDS}
    out["statement_3"] = statement_three(results, s1)
    out["outcome"] = ("the turn: " + "; ".join(f"{k} {outcome_part(s1[k])}" for k in WORLDS)
                      + f"; depth: {out['statement_3']['statement']}")
    out["mean_scores"] = means
    out["neighbours"] = {k: neighbours(results[k], k) for k in WORLDS}
    out["sparsity"] = {k: sparsity_means(results[k]) for k in WORLDS}
    out["thinning"] = {k: {x: results[k]["thinning"].get(x) for x in
                           ("p", "stream", "depth_ratio", "rebuild_largest_difference", "any_cell_emptied")}
                       for k in THINNED}
    return out


def _comparison_line(c: Mapping) -> str:
    return ob._comparison_line(c)


def _paired_line(d: Mapping) -> str:
    return ob._paired_line(d)


def t54_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"Outcome: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    if v["reproduction"] is not None:
        lines.append(f"{GATE}: largest per-seed difference from try T50's committed lib3 scores: "
                     + ", ".join(f"{a} {x:.3g}" for a, x in v["reproduction"]["largest_difference"].items())
                     + f" (gate {REPRO_TOL:g})")
    if r["valid"]:
        lines.append("Statement 1, the turn (each by the single-task test; the argmax among 1, 1.5, 2, 3, 4, 5 and 6, "
                     "smallest on a tie):")
        for k in WORLDS:
            lines.append(f"    in {k}: {outcome_part(r['statement_1'][k])}")
        lines.append("Statement 2, the cases a request would quote:")
        for k in WORLDS:
            for c in r["statement_2"][k].values():
                lines.append(f"    in {k}: " + _comparison_line(c))
        s3 = r["statement_3"]
        null = s3["draw_null"]
        lines.append("Statement 3, self-check (two thinning draws at p = 2/3): "
                     f"{_paired_line(null)}, test: {null['label']}; statement 3 "
                     + ("stands" if s3["valid"] else "is void and only reported"))
        lines.append(f"Statement 3, depth (deep side): {s3['statement']}")
        for key in ("D67", "D50"):
            d = s3[key]
            lines.append(f"    {key}: {_paired_line(d)}, test: {d['label']}")
            lines.append(f"    {key} on the raw metrics (read by nothing): " + ob._raw_line(d["raw_by_metric"]))
        lines.append(f"    beside it the draw null's mean {null['mean']:+.2f}, per seed {null['range'][0]:+.2f} to "
                     f"{null['range'][1]:+.2f}")
        lines.append("    the curve's shape per world side by side (read by nothing): "
                     + "; ".join(f"{k} {outcome_part(s)}" for k, s in s3["shapes"].items()))
        for k, g in s3["anchor_gap"].items():
            lines.append(f"    {k}: mean ceiling minus floor per metric (read by nothing): "
                         + ("not recorded" if g is None else ", ".join(f"{m} {x:+.4g}" for m, x in g.items())))
        for k in THINNED:
            th = r["thinning"][k]
            lines.append(f"world {k}: realised depth ratio {th['depth_ratio']:.4f} (p {THINNED[k]['p']:.4f}, stream "
                         f"{th['stream']}), rebuild checks' largest difference {th['rebuild_largest_difference']:.3g}, "
                         f"any cell left with no count: {th['any_cell_emptied']}")
        for k in WORLDS:
            lines.append(f"world {k}: mean scores " + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
            lines.append("    each against its next smaller: " + ", ".join(
                f"{a} against {n['base']} {n['mean']:+.2f} ({n['label']}), {n['per_tenth']:+.2f} per 0.1 of fraction"
                for a, n in r["neighbours"][k].items()))
            lines.append("    scored zero share and detected scored genes per cell: " + ", ".join(
                f"{a} " + ("not recorded" if sp is None else
                           f"{sp['zero_share']:.4f} / {sp['detected_per_cell']:.2f}")
                for a, sp in r["sparsity"][k].items()))
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t50_file: Path = T50_LIB3) -> dict:
    """Read the worlds present in out/tasks against try T50's committed lib3 file; write t54.{json,txt}."""
    gv.refuse_existing(out, ("t54.json", "t54.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / TASK_FILE.format(world=k) for k in WORLDS
               if (out / "tasks" / TASK_FILE.format(world=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    t50_file = Path(t50_file)
    committed = json.loads(t50_file.read_text()) if t50_file.exists() else {}
    res = {"label": LABEL_T54, "reading": t54_reading(results, committed),
           "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *([t50_file] if t50_file.exists() else [])])
    gv.write_outputs(out, "t54", res, t54_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T54's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one world of try T54")
    p.add_argument("--world", required=True, choices=list(WORLDS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T54's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
