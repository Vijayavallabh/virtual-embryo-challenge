#!/usr/bin/env python
"""Try T56 (experiments/harness/TRIES.md): the common part of Feng's direction at equal movement, where the truth is
shallower and deeper than the reference. A labelled in-sample read on the routes of tries T50
(experiments/t1_m1/over_read.py), T54 (experiments/t1_m1/deep_read.py) and T55 (experiments/t1_m1/part_read.py).
Built to brief AJ.

Try T55 found that M1's step along the common part of Feng's direction alone, at the entry's squared size
(m1_r_gmom_cn), scores above the entry's step on both E10.5 cuts. Tries T50 and T54 found that for a larger step the
gain turns with how deep the truth is against the cells the prediction is built from. The common-part form does not
enlarge the step; row T55 fixed that it is read in the worlds of tries T50 and T54 before anything is said about it.

  run    One world of the chamber cut T1/E9.5->E10.5:chamber, 10 seeds. Worlds (seven):
           lib3       the truth library 3 as try T50 makes it, the inputs as they are;
           lib3_t67   try T50's thinned truths of those names, made by try T50's own code with its own p and streams
           lib3_t33   (over_read.truth_and_contamination; lib3_t33b is try T50's second draw at p = 1/3); the
           lib3_t33b  predictions are lib3's (the inputs as they are);
           ref67      try T54's worlds of those names, made by try T54's own code with its own p and streams
           ref67b     (deep_read.world_reference: the reference thinned in the input pool, the floor's anchors, the
           ref50      scorer's reference, every build and g_mom's fit; ref67b is try T54's second draw at p = 2/3), the
                      truth library 3 as released.
         Arms per seed, all from one cell draw: b0, m1_r_gmom (fraction 0.2) and m1_r_gmom_x2 (0.4) as tries T50 and
         T54 build them (step_read.step_arms on feng_dirs.feng_directions), and try T55's forms built by try T55's own
         functions (part_read.all_directions on the entry's directions as try T55 makes them, and each arm the build of
         the pooled arm of its fraction, gmom_combo.m1_r_arms, with the form's directions in place of the entry's):
         m1_r_gmom_c1, m1_r_gmom_cn (0.2) and m1_r_gmom_cn_x2 (0.4). The forms' directions depend on Feng's statistics,
         the pool's gene names and the held list only, so they are the same in every world; the task file records
         them as try T55 does (sha256, a_n and a_p per group), and per group whether the entry's directions as tries
         T50 and T54 load them and as try T55 loads them have the same sha256. The task file also holds what tries T50
         and T54 hold for the world (the thinning record with ratios only, per seed and arm the prediction's sha256,
         the scored genes' zero share and detected genes per cell, per seed the anchors; the truth's row and
         contamination counts as checks, never printed). Refused on an existing task file, unless the committed N22
         passed, unless experiments/harness/TRIES.md holds row T56, and where try T55's directions refuse.
  read   T56's reading of the seven worlds' files (below, fixed here and repeated in the ledger row before the run).
         Writes t56.json and t56.txt and refuses to overwrite them.

The reading (brief AJ and its amendment of 2026-11-17). Void unless: the seven worlds hold the ten seeds with every arm
and its fraction at one commit on a clean tree; no file records an anchor failure; b0, m1_r_gmom and m1_r_gmom_x2
equal, seed for seed and to 1e-9, their committed scores in experiments/t1_m1/results/t50_m1_over_step/tasks/ for
lib3, lib3_t67, lib3_t33 and lib3_t33b, and in experiments/t1_m1/results/t54_m1_deep_step/tasks/ for ref67, ref67b and
ref50; in lib3 the three forms' arms equal, seed for seed and to 1e-9, their committed scores in
experiments/t1_m1/results/t55_m1_part/tasks/lib3.json; the thinned worlds pass the checks of the try they come from
(their own p and stream, depth ratio within 0.01 of p, no cell emptied, the rebuild check); in lib3_t67, lib3_t33 and
lib3_t33b every arm's prediction sha256 equals, seed for seed, the one in this run's lib3 file; the truth's row and
contamination counts of try T50's thinned worlds equal lib3's (try T50's check, compared by code, never printed); per
group the entry's directions as tries T50 and T54 load them and as try T55 loads them have the same sha256 in every
world; the recorded directions' sha256 and multiples are the same in the seven files and equal those in try T55's lib3
file; no form's arm scores exactly as the entry's arm of its fraction on every seed of a world.
Per world, two differences: m1_r_gmom_cn minus m1_r_gmom, and m1_r_gmom_cn_x2 minus m1_r_gmom_x2, each with its mean,
seeds ahead, per-seed range and its label by try T46's single-task test ("ahead": a mean paired gain above 0.3 with at
least 8 of 10 seeds ahead; "behind" the same reversed; otherwise "level").
The statement: "the common part's gain at equal movement stays positive at every depth read" if in each of the seven
worlds both differences have a positive mean with at least 8 of 10 seeds ahead; "it turns" if in some world a
difference is behind by try T46's test, where that world has a second draw (lib3_t33 and lib3_t33b; ref67 and ref67b)
only if the same difference is behind in both draws, naming each such world and its side ("truth shallower than the
reference" for try T50's worlds, "truth deeper than the reference" for try T54's); otherwise "not separated", naming
every world and difference that meets neither condition, with its side.
Printed beside it and read by nothing: each world's two means; the two draws' differences side by side; each
difference on the raw metrics and its parts on the two de terms and on mmd_u and the variogram; from try T55's
committed files, the same two differences on libraries 1 and 2 as released; m1_r_gmom_c1 minus m1_r_gmom per world;
every arm's mean score per world.

NOT an adoption result: a labelled in-sample diagnostic, as tries T44 to T55 are; the result can only go to the team
captain, and no outcome is evidence on Task 1's read depth, which stays closed. The limits are LIMITS below, repeated in
the text.

Usage: on the DGX campaign checkout, repo root, all seven worlds at one commit, committing nothing until the seven task
files exist, and on the machine tries T50, T54 and T55 ran on, so that the rescored arms reproduce. The environment
variables are docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py run --world lib3 --out t56_cn_depth
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py run --world lib3_t67 --out t56_cn_depth
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py run --world lib3_t33 --out t56_cn_depth
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py run --world lib3_t33b --out t56_cn_depth
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py run --world ref67 --out t56_cn_depth
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py run --world ref67b --out t56_cn_depth
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py run --world ref50 --out t56_cn_depth
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \\
    uv run --frozen --no-sync python -I -B experiments/t1_m1/part_depth_read.py read --out t56_cn_depth
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
import part_read as pr  # noqa: E402

ls, sr, gc, er, zr = ob.ls, ob.sr, ob.gc, ob.er, pr.zr
gv, hygiene, n12, feng_dirs = ob.gv, ob.hygiene, ob.n12, ob.feng_dirs
TASK = ob.TASK  # T1/E9.5->E10.5:chamber
HEART = ob.HEART  # lib3
ROW = "T56"
STEM = "t56"
T50_WORLDS = (HEART, "lib3_t67", "lib3_t33", "lib3_t33b")  # try T50's truths; the inputs as they are
T54_WORLDS = ("ref67", "ref67b", "ref50")  # try T54's thinned references; the truth lib3 as released
WORLDS = (*T50_WORLDS, *T54_WORLDS)
T50_THINNED = tuple(k for k in T50_WORLDS if k != HEART)  # these worlds change the truth only
DRAWS = (("lib3_t33", "lib3_t33b"), ("ref67", "ref67b"))  # a world and its second draw (tries T50's and T54's)
SECOND_DRAW_OF = {k: pair for pair in DRAWS for k in pair}
SHALLOWER, DEEPER = "truth shallower than the reference", "truth deeper than the reference"
SIDE = {**dict.fromkeys(T50_WORLDS, SHALLOWER), **dict.fromkeys(T54_WORLDS, DEEPER)}
THINNING_OF = {**{k: ob.THINNED[k] for k in T50_WORLDS if k != HEART}, **{k: dr.THINNED[k] for k in T54_WORLDS}}
BASE, X2 = sr.BASE, ob.X2  # m1_r_gmom (0.2), m1_r_gmom_x2 (0.4)
UNIFORM_MULTIPLIERS = {X2: 2.0}
C1, CN, CN_X2 = f"{BASE}_c1", f"{BASE}_cn", f"{BASE}_cn_x2"
FORM_ARMS = {a: pr.NEW_ARMS[a] for a in (C1, CN, CN_X2)}  # arm -> (try T55's form, the entry's arm of its fraction)
ENTRY_ARMS = ("b0", BASE, X2)
ARMS = (*ENTRY_ARMS, *FORM_ARMS)
DIFFERENCES = ((CN, BASE), (CN_X2, X2))  # (change, base): the read
REPORTED = ((C1, BASE),)  # read by nothing
FRACTION, FRACTION_TOL, REPRO_TOL = sr.FRACTION, sr.FRACTION_TOL, sr.REPRO_TOL  # 0.2, its tolerance, 1e-9
DEPTH_TOL, REBUILD_TOL = ob.DEPTH_TOL, ob.REBUILD_TOL  # 0.01, 1e-5 (the same in tries T50 and T54)
HOLD_AHEAD = 8  # seeds ahead of 10 for "stays positive"
T50_TASKS = RESULTS / "t50_m1_over_step" / "tasks"
T54_TASKS = RESULTS / "t54_m1_deep_step" / "tasks"
T55_TASKS = RESULTS / "t55_m1_part" / "tasks"
T55_LIB3 = T55_TASKS / "lib3.json"
T55_LIBRARIES = ("lib1", "lib2")  # try T55's committed libraries 1 and 2 as released: printed, read by nothing
TASK_FILE = "{world}.json"
HOLDS = "the common part's gain at equal movement stays positive at every depth read"
TURNS = "it turns"
NOT_SEPARATED = "not separated"
LABEL_T56 = ("try T56, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng step along the common part of Feng's direction alone at the entry's "
             "squared size (try T55's cn form), at the step and twice the step, against the entry's step of the same "
             "fraction, in the worlds of tries T50 (library 3's truth as released and thinned to two thirds and a "
             "third, two draws at a third) and T54 (the reference thinned to two thirds, two draws, and a half) on the "
             "chamber cut; the result can only "
             "go to the team captain, and no outcome is evidence on Task 1's read depth, which stays closed")
LIMITS = (
    "Each difference is taken inside one world, in that world's own points, and is not compared across worlds.",
    "A thinned truth or reference is one draw for all seeds, which the paired difference inside a world shares between "
    "its two arms.",
    "A thinned reference also changes the floor, the de terms' reference, B0's shift and g_mom's fit.",
    "'Stays positive at every depth read' covers truths 0.3 to 1.8 times as deep as the reference and no other.",
    "One released day, in-sample, the chamber cut and one library only.",
    "Equal movement is equal squared size of the direction, not equal counts moved.",
    "The released day's low share along the group-specific parts may be that day's pace and says nothing about the two "
    "days of the test.",
    "Nothing here says how deep any hidden truth is, and no statement may be used with a figure for it.",
    "lib3 reproduces try T55's committed scores by rule and adds no evidence: the statement rests on the six thinned "
    "worlds, all made from the one library where the form was ahead at both fractions.",
    "'Stays positive' says nothing about size: a gain that falls to a few hundredths still stays positive, and each "
    "mean is printed.",
    "'It turns' is a difference behind by 0.3 in a world's own points, which a rescaled floor and ceiling can produce "
    "with no raw metric changing sign, so it is read beside the raw-metric differences.",
    "Try T50's worlds change the truth only and try T54's change the inputs, the floor, the de terms' reference, B0's "
    "shift and g_mom's fit: two manipulations, not one axis.",
)
CHOICES = (
    "b0, m1_r_gmom and m1_r_gmom_x2 are built on the entry's directions as tries T50 and T54 load them "
    "(feng_dirs.feng_directions); the forms on the entry's directions as try T55 loads them (share_read.feng_inputs, "
    "feng_dirs.directions_from_stats); the run records per group, over every group either route gives, whether the "
    "two agree by sha256, and the read voids a world whose record is missing, empty or not all true (amendment).",
    "The forms are built only for the three arms the brief names (c1, cn, cn_x2); try T55's direction record "
    "(part_read.all_directions) covers all five of its forms and is recorded whole, so that it can be compared with "
    "try T55's lib3 file by part_read.direction_values (every form's sha256, a_n and a_p per group).",
    "'Stays positive' uses the difference's mean (above zero) and seeds ahead (strictly above zero, at least 8) as "
    "tries T46 to T55 count them; 'behind' is try T46's label (the reversed single-task test: a mean below -0.3 with "
    "at least 8 of 10 seeds behind), which the amendment's 'at most 2 of 10 seeds ahead' equals only when no seed "
    "ties; the amendment's statement drops the brief's second turn condition (a negative mean with at most 2 seeds "
    "ahead, level by the test).",
    "A difference behind in one draw of a world with a second draw (and not in the other) is not a turn and is named "
    "under 'not separated'; a draw where that difference stays positive is not named. 'It turns' wins over 'not "
    "separated' when both apply, and the differences meeting neither are still recorded.",
    "The side of lib3 is try T50's ('truth shallower than the reference'), as for every try T50 world.",
    "Each file must hold its own world's name, the chamber cut and its truth (the world's name for try T50's worlds, "
    "lib3 for try T54's), and lib3 must record no thinning: the same file checks tries T50 and T54 make.",
    "The truth's row and contamination checks are try T50's own (over_read._sizes_agree: each thinned truth's counts "
    "equal lib3's), run on the T50 worlds only; try T54's worlds, whose truth is lib3 as released, are not checked "
    "against lib3's counts, since neither try makes that check.",
    "The two draws side by side and libraries 1 and 2 from try T55's committed files are computed with try T50's "
    "comparison (over_read.comparison); a library file that is missing or lacks a seed or an arm is printed as not "
    "available and voids nothing.",
)


# ---------------------------------------------------------------------------------------------------------------------
# The worlds and the directions


def world_route(world: str, t, data: Path, genes: Sequence[str], inputs: Mapping[str, ad.AnnData], role: str,
                ref_stage: Stage, read_full: Callable[[], ad.AnnData]
                ) -> tuple[Mapping[str, ad.AnnData], Stage, Stage, int, dict | None]:
    """The world's inputs, scorer's reference stage, truth, contamination count and thinning record. Try T50's worlds:
    the truth by over_read.truth_and_contamination under the world's name, the inputs and the reference as loaded.
    Try T54's worlds: the reference thinned first by deep_read.world_reference, then lib3 as released."""
    if world in T50_WORLDS:
        truth, contamination, thinning = ob.truth_and_contamination(t, data, genes, world)
        return inputs, ref_stage, truth, contamination, thinning
    if world in T54_WORLDS:
        inputs, ref_stage, thinning = dr.world_reference(world, inputs, role, ref_stage, genes, read_full)
        truth, contamination, _ = ob.truth_and_contamination(t, data, genes, HEART)
        return inputs, ref_stage, truth, contamination, thinning
    raise ValueError(f"unknown world {world!r}; one of {list(WORLDS)}")


def world_directions(pool_genes: Sequence[str], held_m1: Sequence[str]) -> tuple[dict, dict, dict, dict]:
    """The entry's directions as tries T50 and T54 load them, try T55's forms built by try T55's own functions and
    their record, and the external record (try T50's, try T55's source, and per group, over every group either route
    gives, whether the two routes' entry directions agree by sha256; the read voids where they do not). Depends on
    the pool's gene names and the held list only. Raises ValueError where try T55's directions refuse."""
    pool_genes = list(pool_genes)
    directions, ext = feng_dirs.feng_directions(pool_genes)
    cfg = feng_dirs.load_config()
    st, src = pr.share_read.feng_inputs()
    entry, _, _ = feng_dirs.directions_from_stats(st, pool_genes, cfg, pr.BRACKET)
    forms, rec = pr.all_directions(entry, pool_genes, held_m1)
    agree = {g: g in directions and g in entry
             and zr.direction_sha256(directions[g]) == zr.direction_sha256(entry[g])
             for g in sorted({*directions, *entry})}
    return directions, forms, rec, {**ext, "t55_source": src, "entry_agrees_with_t55s": agree}


# ---------------------------------------------------------------------------------------------------------------------
# The run


def arm_fractions(fraction: float) -> dict[str, float]:
    fr = sr.arm_fractions(fraction, UNIFORM_MULTIPLIERS)
    return {**fr, **{a: fr[p] for a, (_, p) in FORM_ARMS.items()}}


def build_seed(ctx: harness.BuildContext, b0, directions: Mapping, forms: Mapping[str, Mapping], groups: Mapping,
               fraction: float, held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """One seed's arms: b0, m1_r_gmom and m1_r_gmom_x2 as try T50 builds them (step_read.step_arms), then each form's
    arm as try T55 builds it (the build of the entry's arm of its fraction with the form's directions in place of the
    entry's and every other argument the same)."""
    preds = {"b0": b0(ctx), **sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors,
                                           multipliers=UNIFORM_MULTIPLIERS)}
    fr = arm_fractions(fraction)
    for arm, (key, _) in FORM_ARMS.items():
        preds[arm] = gc.m1_r_arms(ctx, forms[key], groups, fr[arm], held_m1, held_factors)[BASE]
    return preds


def seed_rows(view, cboard, board, truth: Stage, ref_stage: Stage, inputs: Mapping[str, ad.AnnData], floors: Mapping,
              seeds: Sequence[int], directions: Mapping, forms: Mapping, groups: Mapping, fraction: float,
              held_m1: Sequence[str], preregistered: Sequence[str], label: str = "") -> list[dict]:
    """Per seed: the arms built from the world's inputs and scored against the world's truth with the world's reference
    stage and that seed's anchors, recorded as tries T50 and T54 record them."""
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows = []
    for s in seeds:
        t0 = time.time()
        f = floors[s]
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = build_seed(ctx, b0, directions, forms, groups, fraction, held_m1, preregistered)
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "anchors": ob.anchor_record(f.anchors), "arms": arms,
                     "sha256": {k: ls.prediction_sha256(p) for k, p in preds.items()},
                     "sparsity": {k: ls.sparsity(p) for k, p in preds.items()}, "seconds": round(time.time() - t0, 1)})
        print(f"{label} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    return rows


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T56 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T56 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
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
        inputs, ref_stage, truth, contamination, thinning = world_route(
            a.world, t, data, genes, inputs, role, ref_stage, lambda: read_stage(data / t.reference, None))
    except ValueError as e:
        sys.exit(f"T56 is refused on {a.world}: {e}")
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
    stamp_inputs = [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()]]
    res = {"label": LABEL_T56, "task": t.key, "world": a.world, "truth": a.world if a.world in T50_WORLDS else HEART,
           "truth_rows": int(truth.n), "truth_contamination": contamination, "seeds": list(seeds),
           "thinning": thinning}
    floors, failures = ls.all_anchors(t, board, truth, ref_stage, seeds)
    if failures:
        res.update({"anchor_failures": failures, "rows": [], "stamp": gv.stamp(started, stamp_inputs)})
        zr.write_task_file(out, name, res)
        sys.exit(f"T56 is refused on {a.world}: the floor's anchors fail on seeds {[f['seed'] for f in failures]} "
                 "(recorded in the task file; T56 is void)")
    pool = inputs[role]
    preregistered = t1_scope.held_flat_genes()
    union = hygiene.load_union()
    held_m1 = gc.held_feng(union, preregistered)
    try:
        directions, forms, rec, ext = world_directions(list(pool.var_names), held_m1)
    except ValueError as e:
        sys.exit(str(e))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    rows = seed_rows(view, cboard, board, truth, ref_stage, inputs, floors, seeds, directions, forms, groups, fraction,
                     held_m1, preregistered, label=a.world)
    res.update({"fraction": fraction, "fractions": arm_fractions(fraction), "multipliers": dict(UNIFORM_MULTIPLIERS),
                "form_arms": {k: list(v) for k, v in FORM_ARMS.items()}, "external": ext, "directions": rec,
                "rows": rows,
                "stamp": gv.stamp(started, [*stamp_inputs, t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG,
                                            REPO_ROOT / union["file"], feng_dirs.clean_stats_path(), n22])})
    zr.write_task_file(out, name, res)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _rows(r: Mapping) -> list:
    return ob._rows(r)


def committed_file(world: str, t50_tasks: Path = T50_TASKS, t54_tasks: Path = T54_TASKS) -> Path:
    """The committed task file the world's entry arms reproduce: try T50's for its worlds, try T54's for its."""
    return Path(t50_tasks if world in T50_WORLDS else t54_tasks) / TASK_FILE.format(world=world)


def reproduction(rows: Sequence[Mapping], committed: Mapping | None, arms: Sequence[str]) -> dict:
    """Per seed, each arm against its committed score; holds when every seed and arm is there and each is within
    REPRO_TOL."""
    worst = dict.fromkeys(arms, 0.0)
    absent = []
    old = {r["seed"]: r for r in _rows(committed or {})}
    for r in rows:
        o = old.get(r["seed"])
        for arm in arms:
            if o is None or arm not in o.get("arms", {}) or arm not in r["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - o["arms"][arm]["score"]))
    holds = bool(rows) and bool(old) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "absent": absent, "holds": holds}


def _values_key(v: Mapping) -> str:
    return json.dumps(sorted((f"{a}|{b}", x) for (a, b), x in v.items()))


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Mapping], t55_lib3: Mapping | None) -> dict:
    """The void conditions. results: world -> task file; committed: world -> its committed try T50 or T54 file;
    t55_lib3: try T55's committed lib3 file (None when absent)."""
    reasons = []
    missing = [k for k in WORLDS if k not in results]
    if missing:
        reasons.append(f"missing worlds {missing}")
    want_fr = arm_fractions(FRACTION)
    repro = {}
    values = {}
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
        if r.get("task") != TASK:
            reasons.append(f"{k}: the file holds task {r.get('task')!r}")
        if r.get("world") != k:
            reasons.append(f"{k}: the file holds world {r.get('world')!r}")
        want_truth = k if k in T50_WORLDS else HEART
        if r.get("truth") != want_truth:
            reasons.append(f"{k}: the file holds truth {r.get('truth')!r}, not {want_truth}")
        if not isinstance(r.get("fraction"), int | float) or abs(r["fraction"] - FRACTION) > FRACTION_TOL:
            reasons.append(f"{k}: fraction {r.get('fraction')!r}, not {FRACTION:g}")
        got_fr = r.get("fractions") or {}
        bad = [a for a, f in want_fr.items()
               if not isinstance(got_fr.get(a), int | float) or abs(got_fr[a] - f) > FRACTION_TOL]
        if bad:
            reasons.append(f"{k}: arm fractions {got_fr}, not {want_fr} (wrong fraction on {bad})")
        if any(not isinstance((row.get("sha256") or {}).get(a), str) for row in rows for a in ARMS):
            reasons.append(f"{k}: sha256 not recorded for every arm")
        # the entry's arms against try T50's or T54's committed file
        rep = reproduction(rows, committed.get(k), ENTRY_ARMS)
        repro[k] = {"entry_arms": rep}
        if not rep["holds"]:
            src = "try T50's" if k in T50_WORLDS else "try T54's"
            reasons.append(f"{k}: {', '.join(ENTRY_ARMS)} do not reproduce {src} committed scores within "
                           f"{REPRO_TOL:g} per seed")
        if k == HEART:
            rep55 = reproduction(rows, t55_lib3, tuple(FORM_ARMS))
            repro[k]["form_arms"] = rep55
            if not rep55["holds"]:
                reasons.append(f"{k}: {', '.join(FORM_ARMS)} do not reproduce try T55's committed lib3 scores within "
                               f"{REPRO_TOL:g} per seed")
        # the thinning, as the try it comes from checks it
        th = r.get("thinning")
        if k == HEART:
            if th is not None:
                reasons.append(f"{k}: the file records a thinning")
        elif k in THINNING_OF:
            spec, th = THINNING_OF[k], th or {}
            if th.get("p") != spec["p"] or list(th.get("stream") or []) != list(spec["stream"]):
                reasons.append(f"{k}: the file records p {th.get('p')!r} and stream {th.get('stream')!r}, not "
                               f"{spec['p']!r} and {spec['stream']}")
            ratio = th.get("depth_ratio")
            if not isinstance(ratio, int | float) or not abs(ratio - spec["p"]) <= DEPTH_TOL:
                reasons.append(f"{k}: realised depth ratio {ratio!r} is not within {DEPTH_TOL:g} of p "
                               f"{spec['p']:.4f}")
            if th.get("any_cell_emptied") is not False:
                reasons.append(f"{k}: the thinning left a cell with no count (or did not record whether it did): "
                               f"any_cell_emptied {th.get('any_cell_emptied')!r}")
            rebuild = th.get("rebuild_largest_difference")
            if not isinstance(rebuild, int | float) or not rebuild <= REBUILD_TOL:
                reasons.append(f"{k}: the rebuild check's largest difference {rebuild!r} is not at most "
                               f"{REBUILD_TOL:g}")
        # these worlds change the truth only: every arm's prediction is lib3's, seed for seed
        if k in T50_THINNED:
            lib3 = {row["seed"]: row for row in _rows(results.get(HEART, {}))}
            differ = sorted({f"{row['seed']}:{a}" for row in rows for a in ARMS
                             if (row.get("sha256") or {}).get(a) is None
                             or (lib3.get(row["seed"], {}).get("sha256") or {}).get(a) != row["sha256"][a]})
            if HEART not in results or differ:
                reasons.append(f"{k}: the predictions' sha256 do not equal this run's lib3 file seed for seed "
                               f"(seed:arm {differ})")
        # the entry's directions as tries T50 and T54 load them and as try T55 loads them
        agree = (r.get("external") or {}).get("entry_agrees_with_t55s")
        if not isinstance(agree, Mapping) or not agree or not all(x is True for x in agree.values()):
            off = sorted(g for g, x in agree.items() if x is not True) if isinstance(agree, Mapping) else agree
            reasons.append(f"{k}: the entry's directions as tries T50 and T54 load them and as try T55 loads them do "
                           f"not have the same sha256 in every group (or it is not recorded): {off!r}")
        # a form's arm built at the entry's arm (idle form, wrong directions or wrong fraction)
        if rows and not absent:
            same = [a for a, (_, p) in FORM_ARMS.items()
                    if all(row["arms"][a]["score"] == row["arms"][p]["score"] for row in rows)]
            if same:
                reasons.append(f"{k}: form arms score exactly as the entry's arm of their fraction on every seed: "
                               f"{same}")
        values[k] = pr.direction_values(r.get("directions") or {})
        if values[k] is None:
            reasons.append(f"{k}: no full direction record")
    for key, what in (("truth_rows", "cells"), ("truth_contamination", "contamination cells")):
        for b in ob._sizes_agree(results, key):  # try T50's check between its truths (compared, never printed)
            reasons.append(f"the truth's {what}: {b}")
    known = {k: _values_key(v) for k, v in values.items() if v is not None}
    if len(set(known.values())) > 1:
        reasons.append(f"the recorded directions' sha256 or multiples differ between the worlds' files {sorted(known)}")
    t55_values = pr.direction_values((t55_lib3 or {}).get("directions") or {})
    if t55_values is None:
        reasons.append("try T55's lib3 file holds no full direction record")
    else:
        off = sorted(k for k, v in known.items() if v != _values_key(t55_values))
        if off:
            reasons.append(f"the recorded directions' sha256 or multiples differ from try T55's lib3 file in {off}")
    commits = {(r.get("stamp") or {}).get("commit") for r in results.values()}
    unclean = sorted(k for k, r in results.items() if (r.get("stamp") or {}).get("tree_clean") is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the worlds ran at {len(commits)} commits")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def comparison(r: Mapping, world: str, change: str, base: str) -> dict:
    """Try T50's comparison (try T46's single-task test on one world, its sd, the reversed test, the points by part
    and the raw metric differences) with the label ahead / behind / level."""
    return ob.comparison(r, world, change, base)


def holds(c: Mapping) -> bool:
    """'Stays positive': a positive mean with at least HOLD_AHEAD of 10 seeds ahead."""
    return c["mean"] > 0 and c["seeds_ahead"] >= HOLD_AHEAD


def behind(c: Mapping) -> bool:
    """Behind by try T46's single-task test (the reversed test)."""
    return c["label"] == "behind"


def turns(diffs: Mapping[str, Mapping[str, Mapping]], world: str, name: str) -> bool:
    """The difference counts as a turn in this world: behind, and where the world has a second draw, behind in both
    draws."""
    return all(behind(diffs[k][name]) for k in SECOND_DRAW_OF.get(world, (world,)))


def _item(world: str, name: str) -> dict:
    return {"world": world, "difference": name, "side": SIDE[world]}


def _items_text(items: Sequence[Mapping]) -> str:
    return "; ".join(f"{x['world']} ({x['side']}): {x['difference']}" for x in items)


def statement(diffs: Mapping[str, Mapping[str, Mapping]]) -> dict:
    """diffs: world -> difference name -> comparison. HOLDS when every difference in every world stays positive; TURNS
    when a difference is behind in some world (in both draws where the world has two), naming each such world and its
    side; otherwise NOT_SEPARATED, naming every world and difference that meets neither, with its side."""
    turned = [_item(k, n) for k in WORLDS for n in diffs[k] if turns(diffs, k, n)]
    neither = [_item(k, n) for k in WORLDS for n, c in diffs[k].items() if not holds(c) and not turns(diffs, k, n)]
    if all(holds(c) for k in WORLDS for c in diffs[k].values()):
        said = HOLDS
    elif turned:
        said = f"{TURNS} ({_items_text(turned)})"
    else:
        said = f"{NOT_SEPARATED} ({_items_text(neither)})"
    return {"statement": said, "turned": turned, "neither": neither,
            "turned_in": sorted({x["world"] for x in turned}, key=WORLDS.index),
            "neither_in": sorted({x["world"] for x in neither}, key=WORLDS.index)}


def _name(a: str, b: str) -> str:
    return f"{a} minus {b}"


def library_differences(t55_libs: Mapping[str, Mapping | None]) -> dict:
    """From try T55's committed files, the two differences on libraries 1 and 2 as released (read by nothing); None
    where a file is missing or lacks a seed or an arm."""
    out = {}
    need = {a for pair in DIFFERENCES for a in pair}
    for k in T55_LIBRARIES:
        r = t55_libs.get(k)
        rows = _rows(r or {})
        complete = sorted(row["seed"] for row in rows) == list(SEEDS) and all(
            a in row.get("arms", {}) for row in rows for a in need)
        out[k] = {_name(a, b): comparison(r, k, a, b) for a, b in DIFFERENCES} if complete else None
    return out


def draws_side_by_side(diffs: Mapping[str, Mapping[str, Mapping]]) -> dict:
    """Per world with a second draw and per difference, the two draws' comparisons."""
    return {f"{a} and {b}": {n: {a: diffs[a][n], b: diffs[b][n]} for n in diffs[a]} for a, b in DRAWS}


def t56_reading(results: Mapping[str, Mapping], committed: Mapping[str, Mapping], t55_lib3: Mapping | None,
                t55_libs: Mapping[str, Mapping | None] | None = None) -> dict:
    v = validity(results, committed, t55_lib3)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    diffs = {k: {_name(a, b): comparison(results[k], k, a, b) for a, b in DIFFERENCES} for k in WORLDS}
    st = statement(diffs)
    out.update({"outcome": st["statement"], "statement": st, "differences": diffs,
                "means": {k: {n: c["mean"] for n, c in diffs[k].items()} for k in WORLDS},
                "draws": draws_side_by_side(diffs), "libraries_1_and_2": library_differences(t55_libs or {}),
                "reported": {k: {_name(a, b): comparison(results[k], k, a, b) for a, b in REPORTED} for k in WORLDS},
                "mean_scores": {k: {a: float(np.mean([row["arms"][a]["score"] for row in _rows(results[k])]))
                                    for a in ARMS} for k in WORLDS},
                "thinning": {k: {x: (results[k].get("thinning") or {}).get(x) for x in
                                 ("p", "stream", "depth_ratio", "rebuild_largest_difference", "any_cell_emptied")}
                             for k in THINNING_OF},
                "directions": results[HEART]["directions"]})
    return out


def _raw_line(c: Mapping) -> str:
    return ", ".join(f"{m} {c['metrics'].get(m, float('nan')):+.4g}" for m in ls.WEIGHTS)


def _difference_line(c: Mapping) -> str:
    return (f"{c['change']} minus {c['base']}: {c['mean']:+.2f} points, ahead {c['seeds_ahead']} of {c['of']}, per "
            f"seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}, test: {c['label']}; de terms "
            f"{c['points']['de']:+.2f}, mmd_u and variogram {c['points']['mmd_u_variogram']:+.2f} (read by nothing); "
            f"raw metrics {_raw_line(c)} (read by nothing)")


def t56_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"Outcome: {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    for k, rep in v["reproduction"].items():
        for what, x in rep.items():
            lines.append(f"{k}: {what.replace('_', ' ')}, largest per-seed difference from the committed scores: "
                         + ", ".join(f"{a} {d:.3g}" for a, d in x["largest_difference"].items())
                         + f" (gate {REPRO_TOL:g})")
    if r["valid"]:
        st = r["statement"]
        lines.append(f"The statement: {st['statement']}")
        if st["turned"] and st["neither"]:
            lines.append("    also meeting neither condition: " + _items_text(st["neither"]))
        lines.append("    each world's two means: " + "; ".join(
            f"{k} ({SIDE[k]}): " + ", ".join(f"{n} {m:+.2f}" for n, m in r["means"][k].items()) for k in WORLDS))
        lines.append("The two draws side by side:")
        for pair, by_name in r["draws"].items():
            for n, cs in by_name.items():
                lines.append(f"    {pair}, {n}: " + " | ".join(
                    f"{k} {c['mean']:+.2f}, ahead {c['seeds_ahead']} of {c['of']}, per seed {c['range'][0]:+.2f} to "
                    f"{c['range'][1]:+.2f}, test: {c['label']}" for k, c in cs.items()))
        lines.append("Per world, the two differences at equal movement (try T46's single-task test inside the world):")
        for k in WORLDS:
            for c in r["differences"][k].values():
                lines.append(f"    in {k}: " + _difference_line(c))
        lines.append("Reported and read by nothing:")
        for k, cs in r["libraries_1_and_2"].items():
            if cs is None:
                lines.append(f"    {k} as released (try T55's committed file): not available")
                continue
            for c in cs.values():
                lines.append(f"    {k} as released (try T55's committed file): " + _difference_line(c))
        for k in WORLDS:
            for c in r["reported"][k].values():
                lines.append(f"    in {k}: " + _difference_line(c))
        for k in WORLDS:
            lines.append(f"    world {k}: mean scores "
                         + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
        for k, th in r["thinning"].items():
            lines.append(f"    world {k}: realised depth ratio {th['depth_ratio']:.4f} (p {THINNING_OF[k]['p']:.4f}, "
                         f"stream {th['stream']}), rebuild check's largest difference "
                         f"{th['rebuild_largest_difference']:.3g}, any cell left with no count: "
                         f"{th['any_cell_emptied']}")
        d = r["directions"]
        lines.append("    the forms' directions (the same in the seven files and in try T55's lib3 file): " + "; ".join(
            f"{g}: a_n {e[pr.N]:.4f}, a_p {e[pr.P]:.4f}" for g, e in d["entry"].items()))
    lines.append("Choices the brief left open:")
    lines.extend(f"  - {x}" for x in CHOICES)
    lines.extend(f"Limit: {x}" for x in LIMITS)
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain; no outcome is evidence on Task 1's read depth, which stays closed")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t50_tasks: Path = T50_TASKS, t54_tasks: Path = T54_TASKS, t55_tasks: Path = T55_TASKS) -> dict:
    """Read the worlds present in out/tasks against the committed files of tries T50, T54 and T55 (t55_tasks holds
    try T55's lib3, lib1 and lib2 files); write t56.{json,txt}."""
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / TASK_FILE.format(world=k) for k in WORLDS
               if (out / "tasks" / TASK_FILE.format(world=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    files = {k: committed_file(k, t50_tasks, t54_tasks) for k in WORLDS}
    found = {k: p for k, p in files.items() if p.exists()}
    committed = {k: json.loads(p.read_text()) for k, p in found.items()}
    t55_files = {k: Path(t55_tasks) / TASK_FILE.format(world=k) for k in (HEART, *T55_LIBRARIES)}
    t55_found = {k: p for k, p in t55_files.items() if p.exists()}
    t55 = {k: json.loads(p.read_text()) for k, p in t55_found.items()}
    res = {"label": LABEL_T56, "choices": list(CHOICES), "limits": list(LIMITS),
           "reading": t56_reading(results, committed, t55.get(HEART), {k: t55.get(k) for k in T55_LIBRARIES}),
           "stamps": {k: v.get("stamp") for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *found.values(), *t55_found.values()])
    gv.write_outputs(out, STEM, res, t56_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T56's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one world of try T56")
    p.add_argument("--world", required=True, choices=list(WORLDS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T56's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
