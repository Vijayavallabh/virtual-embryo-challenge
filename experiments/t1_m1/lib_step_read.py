#!/usr/bin/env python
"""Try T46 (experiments/harness/TRIES.md): M1's Feng step by the kind of truth, and beyond twice the step. A labelled
read on try T44's route (experiments/t1_m1/step_read.py), with the chamber cut's truth taken pooled and library by
library.

Try T44 read M1's step at 1.5 and 2 times its fraction on the two E10.5 cuts and found both better on both cuts, the
response still rising at twice the step. Its truth differs from its reference in the kind of sample as well as in age.
The released E10.5 file's three libraries are of two kinds: libraries 1 and 2 are embryo portions (the kind E9.5's two
libraries are, shallower than E9.5), library 3 a dissected-heart sample. The test pair is dissected heart to dissected
heart, so what a larger step gains through the change of kind does not carry to the test; what it gains through a day
of age may. This read puts the same steps, and two larger ones, to each library as the truth.

  run    One truth of the chamber cut T1/E9.5->E10.5:chamber, 10 seeds, on T44's draws and scoring. `pooled`: the cut's
         target as T44's run loads it (harness.filter_target). `lib1`, `lib2`, `lib3`: the same filtered target with
         only the rows whose cell name's library (pseudo_task.library_of: the suffix after the last '_') is 1, 2 or 3,
         in their order; the restriction comes after the task's own filter, so a library keeps its share of the
         contamination cells the filter draws. The board, the floor's anchors and the scoring are then made for that
         truth as the harness makes them for any target (harness.task_board on the truth's row count); the reference
         stage, the inputs and every build argument are the cut's. Arms per seed, all from one cell draw: b0 and
         m1_r_gmom (fraction 0.2) as T42 and T44 build them, and step_read.step_arms at multipliers 1.5, 2, 3 and 4
         (m1_r_gmom_x15, _x2, _x3, _x4: fractions 0.3, 0.4, 0.6, 0.8), nothing else changed. The task file holds every
         arm's fraction, the truth's name, the truth's row count (a check, compared by code and never printed) and, per
         seed and arm, the sha256 of the predicted matrix as the scorer reads it. Refused on an existing task file,
         unless the committed N22 passed, and unless experiments/harness/TRIES.md holds row T46.
  read   T46's reading of the four truths' files (rule below, fixed in the row before the run).

Where a prediction depends on the truth: only through the board. harness.task_board gives the cut's synthetic board
max_cells = max(5000, 3 x max(truth_cells, ref_cells)) and min_cells = min(1000, truth_cells), with truth_cells half of
10% of the truth's rows; the candidate's board keeps both limits, and count_for draws max_cells cells. A library truth
small enough that 3 x truth_cells is at most max(5000, 3 x ref_cells) leaves the count, and so every prediction, as on
the pooled truth; otherwise the count, every draw and every prediction differ. The sha256 agreement across truths, which
the read reports, says which held.

NOT an adoption result: Feng's E9.5 cells lie inside this pseudo-task's external window, the cut is one earlier T1
choices were made on, and a library truth is a part of a released target read in-sample. The result can only go to
the team captain.

The reading. Void unless the four truths hold the ten seeds at one commit on a clean tree; each file's task is the
chamber cut, its truth its own name, its base fraction 0.2 and each arm's fraction its multiplier times 0.2; on
`pooled`, b0 and m1_r_gmom equal, seed for seed, their committed scores in
experiments/t1_m1/results/t42_m1_eb/tasks/ and m1_r_gmom_x15 and m1_r_gmom_x2 theirs in
experiments/t1_m1/results/t44_m1_step/tasks/ (the chamber cut's files), each to 1e-9; the three library truths' row
counts add up to the pooled truth's; and no scaled arm scores exactly as m1_r_gmom (the same record) on every seed of a
truth. The test on one truth is the adoption rule's single-task test (vec.adoption.keep_change on that truth alone): a
mean paired gain above 0.3 with at least 8 of 10 seeds ahead, scaled arm against m1_r_gmom ("ahead"); "behind" is the
same with the sign reversed. On the kind, from twice the step (m1_r_gmom_x2): "the larger step's gain is the kind's" if
it is ahead on lib3 and on neither lib1 nor lib2; "the larger step's gain is not the kind's" if it is ahead on both lib1
and lib2; otherwise "mixed". The text always gives the three library gains of twice the step and, for each, the part on
the two de terms and the part on mmd_u and the variogram, in board points (100 x the board's weight x the mean skill
difference, a NaN skill counting 0 as on the server), because the portion libraries are shallower than E9.5 and count
noise is rewarded more there. Where the response turns, per truth: the multiplier with the highest mean board score
among 1, 1.5, 2, 3 and 4 (the smallest on a tie), and "still rising at 4" if that is 4. Reported: per truth every arm's
mean board score, each scaled arm's mean paired difference from m1_r_gmom with seeds ahead and per-seed range, every
metric's mean skill difference, and whether the predictions' sha256 agree across the truths. The read puts nothing
forward by itself.

Usage: on the DGX campaign checkout, repo root, all four truths at one commit, committing nothing until the four task
files exist, and on the machine T42 and T44 ran on, so that the pooled truth reproduces. The environment variables are
docs/COMPUTE.md's for repository Python on the DGX.
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/lib_step_read.py run --truth pooled --out t46_m1_lib_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/lib_step_read.py run --truth lib1 --out t46_m1_lib_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/lib_step_read.py run --truth lib2 --out t46_m1_lib_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/lib_step_read.py run --truth lib3 --out t46_m1_lib_step
  NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/lib_step_read.py read --out t46_m1_lib_step
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np

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
sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_scope"))
import eb_read as er  # noqa: E402
import feng_dirs  # noqa: E402
import gmom_combo as gc  # noqa: E402
import pseudo_task as pt  # noqa: E402
import step_read as sr  # noqa: E402

gv, hygiene, n12 = gc.gv, gc.hygiene, er.n12
TASK = gc.TASK_CHAMBER
CUT_FILE = gc.CUT_FILES[TASK]
T42_FILE = sr.T42_TASKS / CUT_FILE
T44_FILE = RESULTS / "t44_m1_step" / "tasks" / CUT_FILE
POOLED = "pooled"
LIBRARIES = {"lib1": "1", "lib2": "2", "lib3": "3"}  # truth -> library tag (pseudo_task.library_of)
PORTIONS, HEART = ("lib1", "lib2"), "lib3"
TRUTHS = (POOLED, *LIBRARIES)
BASE = sr.BASE  # m1_r_gmom, fraction 0.2
MULTIPLIERS = {"m1_r_gmom_x15": 1.5, "m1_r_gmom_x2": 2.0, "m1_r_gmom_x3": 3.0, "m1_r_gmom_x4": 4.0}
ARMS = ("b0", BASE, *MULTIPLIERS)
TWICE = "m1_r_gmom_x2"
ROW = "T46"
FRACTION = sr.FRACTION  # 0.2
FRACTION_TOL = sr.FRACTION_TOL
REPRO_TOL = sr.REPRO_TOL  # 1e-9
REPRO = {"b0": "t42", BASE: "t42", "m1_r_gmom_x15": "t44", "m1_r_gmom_x2": "t44"}  # arm -> committed file, on pooled
KIND = "the larger step's gain is the kind's"
NOT_KIND = "the larger step's gain is not the kind's"
MIXED = "mixed"
STILL_RISING = "still rising at 4"
DE_TERMS, DIST_TERMS = ("de_score", "de_direction"), ("mmd_u", "variogram")
WEIGHTS = {m: gw * mw for _, gw, ms in GROUPS["T1"] for m, mw in ms}  # the board's weight of each metric
TASK_FILE = "{truth}.json"
LABEL_T46 = ("try T46, a labelled in-sample diagnostic, never an adoption result: M1 (with R1, R2 and g_mom's factors, "
             "the entry's form) with its Feng step multiplied by 1, 1.5, 2, 3 and 4, on T44's route, with the chamber "
             "cut's truth taken pooled and library by library (libraries 1 and 2 embryo portions, library 3 a "
             "dissected heart); Feng's E9.5 cells lie inside this pseudo-task's external window and the cut is one "
             "earlier T1 choices were made on; the result can only go to the team captain")


# ---------------------------------------------------------------------------------------------------------------------
# The truths


@dataclass
class _NamedStage(Stage):
    """A Stage that also keeps the cell names of the AnnData it was made from (used only while filter_target runs)."""

    names: np.ndarray | None = None

    @classmethod
    def from_anndata(cls, a, genes, need_coords=False):
        s = Stage.from_anndata(a, genes, need_coords=need_coords)
        return cls(X=s.X, coords=s.coords, celltype=s.celltype, names=np.asarray(a.obs_names).astype(str))


@contextmanager
def _named_stages() -> Iterator[None]:
    old = harness.Stage
    harness.Stage = _NamedStage
    try:
        yield
    finally:
        harness.Stage = old


def filtered_target(t: harness.PseudoTask, data: Path, genes) -> tuple[Stage, np.ndarray]:
    """harness.filter_target's Stage (the very rows, in order) and the cell names of those rows."""
    if t.target_filter is None:
        raise ValueError(f"{t.key} has no target filter; T46 reads the chamber cut only")
    with _named_stages():
        named, _ = harness.filter_target(t, data, genes)
    return Stage(X=named.X, coords=named.coords, celltype=named.celltype), named.names


def library_rows(names: Sequence[str], truth: str) -> np.ndarray:
    """Rows, in order, whose cell name's library is the truth's."""
    return np.flatnonzero(pt.library_of(names) == LIBRARIES[truth])


def restrict(target: Stage, names: Sequence[str], truth: str) -> Stage:
    return target.take(library_rows(names, truth))


def load_truth(t: harness.PseudoTask, data: Path, genes, truth: str) -> Stage:
    """`pooled`: the cut's target exactly as T44's run loads it; a library: that target's rows of that library."""
    if truth == POOLED:
        return harness.filter_target(t, data, genes)[0]
    if truth not in LIBRARIES:
        raise ValueError(f"unknown truth {truth!r}; one of {list(TRUTHS)}")
    target, names = filtered_target(t, data, genes)
    return restrict(target, names, truth)


def prediction_sha256(pred: Stage) -> str:
    """sha256 of the predicted matrix as the scorer reads it (dtype, shape and bytes)."""
    X = np.ascontiguousarray(pred.X)
    h = hashlib.sha256()
    h.update(f"{X.dtype.str}{X.shape}".encode())
    h.update(X.tobytes())
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# The run


def arm_fractions(fraction: float) -> dict[str, float]:
    return sr.arm_fractions(fraction, MULTIPLIERS)


def scaled_arms(ctx: harness.BuildContext, directions: Mapping, groups: Mapping, fraction: float,
                held_m1: Sequence[str], held_factors: Sequence[str]) -> dict[str, Stage]:
    """m1_r_gmom and the four scaled arms: step_read.step_arms at this read's multipliers."""
    return sr.step_arms(ctx, directions, groups, fraction, held_m1, held_factors, multipliers=MULTIPLIERS)


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows():
        sys.exit(f"T46 is refused: {ROW} has no row in experiments/harness/TRIES.md (a row comes before its run)")
    n22 = RESULTS / a.n22 / "n22.json"
    if not gc.n22_passed(n22):
        sys.exit(f"T46 is refused: {gv._rel(n22)} does not hold a passing N22 (T26's route needs N22 passed)")
    out = RESULTS / a.out
    name = TASK_FILE.format(truth=a.truth)
    gv.refuse_existing(out / "tasks", (name,))
    started = gv.tree_state()
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    truth = load_truth(t, data, genes, a.truth)
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, truth.n, ref_stage.n)
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
        f = gv.floor_anchors(t, board, truth, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        preds = {"b0": b0(ctx), **scaled_arms(ctx, directions, groups, fraction, held_m1, preregistered)}
        arms = er.score_arms(board, truth, ref_stage, s, f.anchors, preds)
        rows.append({"seed": s, "floor": float(f.score), "arms": arms,
                     "sha256": {k: prediction_sha256(p) for k, p in preds.items()},
                     "seconds": round(time.time() - t0, 1)})
        print(f"{a.truth} seed {s}: " + ", ".join(f"{k} {v['score']:.2f}" for k, v in arms.items()), flush=True)
    res = {"label": LABEL_T46, "task": t.key, "truth": a.truth, "truth_rows": int(truth.n), "fraction": fraction,
           "fractions": arm_fractions(fraction), "multipliers": dict(MULTIPLIERS), "external": ext,
           "seeds": list(seeds), "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, ts.R1_CONFIG, REPO_ROOT / union["file"],
                                       feng_dirs.clean_stats_path(), n22])}
    (out / "tasks").mkdir(parents=True, exist_ok=True)
    (out / "tasks" / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def reproduction(rows: Sequence[Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    """On the pooled truth, per seed, each arm of REPRO against its committed file's score (committed: 't42' / 't44' ->
    that file's rows); holds when every seed and arm is there and each is within REPRO_TOL."""
    by_seed = {src: {r["seed"]: r for r in rs} for src, rs in committed.items()}
    worst = {a: 0.0 for a in REPRO}
    absent = []
    for r in rows:
        for arm, src in REPRO.items():
            old = by_seed.get(src, {}).get(r["seed"])
            if old is None or arm not in old["arms"] or arm not in r["arms"]:
                absent.append(f"{r['seed']}:{arm}")
                continue
            worst[arm] = max(worst[arm], abs(r["arms"][arm]["score"] - old["arms"][arm]["score"]))
    holds = bool(rows) and not absent and all(w <= REPRO_TOL for w in worst.values())
    return {"largest_difference": worst, "absent": absent, "holds": holds}


def validity(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    """The void conditions. results: truth -> task file (a missing truth is absent); committed: as `reproduction`."""
    reasons = []
    missing = [k for k in TRUTHS if k not in results]
    if missing:
        reasons.append(f"missing truths {missing}")
    want_fr = arm_fractions(FRACTION)
    for k, r in results.items():
        seeds = [row["seed"] for row in r["rows"]]
        if sorted(seeds) != list(SEEDS) or list(r.get("seeds", seeds)) != list(SEEDS):
            reasons.append(f"{k}: seeds {seeds}, not {list(SEEDS)}")
        absent = sorted({a for row in r["rows"] for a in ARMS if a not in row["arms"]})
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
    commits = {r["stamp"]["commit"] for r in results.values()}
    unclean = sorted(k for k, r in results.items() if r["stamp"]["tree_clean"] is not True)
    if unclean:
        reasons.append(f"tree not clean on {unclean}")
    if len(commits) > 1:
        reasons.append(f"the truths ran at {len(commits)} commits")
    repro = None
    if POOLED in results:
        repro = reproduction(results[POOLED]["rows"], committed)
        if not repro["holds"]:
            reasons.append(f"{POOLED}: {', '.join(REPRO)} do not reproduce their committed scores (T42's, T44's) "
                           f"within {REPRO_TOL:g} per seed")
    if all(k in results for k in TRUTHS):
        sizes = [results[k].get("truth_rows") for k in TRUTHS]
        if not all(isinstance(x, int) for x in sizes) or sum(sizes[1:]) != sizes[0]:
            reasons.append("the library truths' cells do not add up to the pooled truth's")
    for k, r in results.items():
        if not all(a in row["arms"] for row in r["rows"] for a in ARMS):
            continue
        same = [a for a in MULTIPLIERS if all(row["arms"][a] == row["arms"][BASE] for row in r["rows"])]
        if same:
            reasons.append(f"{k}: scaled arms score exactly as {BASE} on every seed (built at the base fraction?): "
                           f"{same}")
    return {"valid": not reasons, "reasons": reasons, "reproduction": repro,
            "commits": sorted(c for c in commits if c)}


def _nz(x) -> float:
    x = float(x)
    return x if np.isfinite(x) else 0.0


def points_by_part(rows: Sequence[Mapping], base: str, change: str) -> dict:
    """The mean paired gain in board points split by metric (100 x weight x skill difference, NaN as 0), and its de
    part (de_score, de_direction) and distribution part (mmd_u, variogram)."""
    per = {m: float(np.mean([100.0 * w * (_nz(row["arms"][change]["skills"].get(m, np.nan))
                                          - _nz(row["arms"][base]["skills"].get(m, np.nan))) for row in rows]))
           for m, w in WEIGHTS.items()}
    return {"by_metric": per, "de": sum(per[m] for m in DE_TERMS), "mmd_u_variogram": sum(per[m] for m in DIST_TERMS)}


def truth_comparison(r: Mapping, truth: str, arm: str) -> dict:
    """A scaled arm against m1_r_gmom on one truth: step_read.compare (keep_change on this truth alone, the mean paired
    difference, seeds ahead, per-seed range, metric differences), the reversed test, and the points by part."""
    fwd = sr.compare({truth: r}, BASE, arm)
    rev = sr.compare({truth: r}, arm, BASE)
    cut = fwd["cuts"][truth]
    ahead, behind = bool(fwd["rule"]["tasks"][truth]["passes"]), bool(rev["rule"]["tasks"][truth]["passes"])
    return {"base": BASE, "change": arm, "mean": cut["mean"], "seeds_ahead": cut["ahead"], "of": cut["of"],
            "range": cut["range"], "skills": cut["skills"], "groups": cut["groups"], "metrics": cut["metrics"],
            "ahead": ahead, "behind": behind, "test": "ahead" if ahead else "behind" if behind else "neither",
            "points": points_by_part(r["rows"], BASE, arm)}


def kind_statement(comps: Mapping[str, Mapping[str, Mapping]]) -> str:
    """comps: truth -> arm -> truth_comparison. From twice the step on the three libraries."""
    ahead = {k: comps[k][TWICE]["ahead"] for k in LIBRARIES}
    if ahead[HEART] and not any(ahead[k] for k in PORTIONS):
        return KIND
    if all(ahead[k] for k in PORTIONS):
        return NOT_KIND
    return MIXED


def turn(mean_scores: Mapping[str, float]) -> dict:
    """The multiplier with the highest mean board score among 1 (m1_r_gmom), 1.5, 2, 3 and 4; the smallest on a tie."""
    by_m = {1.0: mean_scores[BASE], **{m: mean_scores[a] for a, m in MULTIPLIERS.items()}}
    best = max(by_m, key=lambda m: (by_m[m], -m))
    top = max(MULTIPLIERS.values())
    return {"by_multiplier": {f"{m:g}": v for m, v in by_m.items()}, "best": best,
            "statement": STILL_RISING if best == top else f"turns at {best:g}"}


def sha_agreement(results: Mapping[str, Mapping]) -> dict:
    """Per arm, whether every seed's predicted-matrix sha256 is the same on every truth present; seeds that differ."""
    out = {}
    for a in ARMS:
        differ, unrecorded = [], []
        for s in SEEDS:
            shas = []
            for r in results.values():
                row = next((x for x in r["rows"] if x["seed"] == s), None)
                shas.append(None if row is None else row.get("sha256", {}).get(a))
            if any(x is None for x in shas):
                unrecorded.append(s)
            elif len(set(shas)) > 1:
                differ.append(s)
        out[a] = {"agree": not differ and not unrecorded, "seeds_differing": differ, "seeds_unrecorded": unrecorded}
    return {"arms": out, "all_agree": all(v["agree"] for v in out.values())}


def t46_reading(results: Mapping[str, Mapping], committed: Mapping[str, Sequence[Mapping]]) -> dict:
    v = validity(results, committed)
    out = {"valid": v["valid"], "validity": v}
    if not v["valid"]:
        out["outcome"] = "VOID"
        return out
    comps = {k: {a: truth_comparison(results[k], k, a) for a in MULTIPLIERS} for k in TRUTHS}
    means = {k: {a: float(np.mean([row["arms"][a]["score"] for row in results[k]["rows"]])) for a in ARMS}
             for k in TRUTHS}
    out["outcome"] = kind_statement(comps)
    out["kind"] = {k: {"gain": comps[k][TWICE]["mean"], "de": comps[k][TWICE]["points"]["de"],
                       "mmd_u_variogram": comps[k][TWICE]["points"]["mmd_u_variogram"], "test": comps[k][TWICE]["test"]}
                   for k in LIBRARIES}
    out["comparisons"] = comps
    out["mean_scores"] = means
    out["turn"] = {k: turn(means[k]) for k in TRUTHS}
    out["sha256_across_truths"] = sha_agreement(results)
    return out


def t46_text(res: Mapping) -> str:
    """The reading as text; no cell count is printed."""
    r = res["reading"]
    lines = [res["label"], f"On the kind (twice the step): {r['outcome']}"]
    v = r["validity"]
    if not r["valid"]:
        lines.append("Void: " + "; ".join(v["reasons"]))
    if v["reproduction"] is not None:
        w = v["reproduction"]["largest_difference"]
        lines.append(f"{POOLED}: largest per-seed difference from the committed scores: "
                     + ", ".join(f"{a} {x:.3g}" for a, x in w.items()) + f" (gate {REPRO_TOL:g})")
    for k, g in r.get("kind", {}).items():
        lines.append(f"    {TWICE} on {k}: {g['gain']:+.2f} points ({g['test']}); de terms {g['de']:+.2f}, mmd_u and "
                     f"variogram {g['mmd_u_variogram']:+.2f}")
    for k in r.get("mean_scores", {}):
        lines.append(f"truth {k}: mean scores " + ", ".join(f"{a} {x:.2f}" for a, x in r["mean_scores"][k].items()))
        tu = r["turn"][k]
        lines.append(f"    response: {tu['statement']} (mean score by multiplier "
                     + ", ".join(f"{m} {x:.2f}" for m, x in tu["by_multiplier"].items()) + ")")
        for a, c in r["comparisons"][k].items():
            sk = ", ".join(f"{m} {x:+.4f}" for m, x in c["skills"].items())
            lines.append(f"    {a} against {BASE}: {c['mean']:+.2f} points, ahead {c['seeds_ahead']} of {c['of']}, per "
                         f"seed {c['range'][0]:+.2f} to {c['range'][1]:+.2f}, test: {c['test']}; de terms "
                         f"{c['points']['de']:+.2f}, mmd_u and variogram {c['points']['mmd_u_variogram']:+.2f}; "
                         f"skills: {sk}")
    if "sha256_across_truths" in r:
        sa = r["sha256_across_truths"]
        lines.append("predictions' sha256 across the truths: "
                     + ("the same on every seed and arm" if sa["all_agree"] else
                        "differ: " + "; ".join(f"{a} on seeds {x['seeds_differing'] or x['seeds_unrecorded']}"
                                               for a, x in sa["arms"].items() if not x["agree"])))
    lines.append("Limit: in-sample labelled diagnostic; never adoption evidence; puts nothing forward by itself; it "
                 "can only go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, t42_file: Path = T42_FILE, t44_file: Path = T44_FILE) -> dict:
    """Read the truths present in out/tasks against T42's and T44's committed chamber files; write t46.{json,txt}."""
    gv.refuse_existing(out, ("t46.json", "t46.txt"))
    started = gv.tree_state()
    present = {k: out / "tasks" / TASK_FILE.format(truth=k) for k in TRUTHS
               if (out / "tasks" / TASK_FILE.format(truth=k)).exists()}
    results = {k: json.loads(p.read_text()) for k, p in present.items()}
    files = {src: p for src, p in (("t42", Path(t42_file)), ("t44", Path(t44_file))) if p.exists()}
    committed = {src: json.loads(p.read_text())["rows"] for src, p in files.items()}
    res = {"label": LABEL_T46, "reading": t46_reading(results, committed),
           "stamps": {k: v["stamp"] for k, v in results.items()}}
    res["stamp"] = gv.stamp(started, [*present.values(), *files.values()])
    gv.write_outputs(out, "t46", res, t46_text(res))
    return res


def cmd_read(a) -> None:
    res = read(RESULTS / a.out)
    if res["reading"]["outcome"] == "VOID":
        sys.exit("T46's reading is void and goes to the captain")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="one truth of try T46")
    p.add_argument("--truth", required=True, choices=list(TRUTHS))
    p.add_argument("--out", required=True)
    p.add_argument("--n22", default="n22_additivity")
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T46's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
