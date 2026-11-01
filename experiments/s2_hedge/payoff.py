#!/usr/bin/env python
"""s2 hedge payoffs (docs/PLAN.md section 8): a board's score when the prediction's tissue mix and the truth's differ.

Sensitivity analysis for the s2 pre-registration (17-23 Nov). Not an adoption result and not a choice: nothing here
picks an entry, and every setting below was fixed before the first run. On a released pseudo-task the truth is
re-sampled without replacement, stratified by category, to a grid of tissue mixes along one axis; the prediction (B0,
as the harness builds it) is built at the same grid of mixes; every prediction is scored against every truth with the
harness's scoring (bake protocol, scorer build main, 10% working copies, each truth's floor and ceiling computed per
seed and shared by the predictions, build seed 1,000,000 + seed). The payoff matrix score[prediction mix][truth mix]
shows what a hedge gains where the truth's mix moves away from the one s1 assumes, and what it costs where it does not.

Axes (T1/E9.5->E10.5:chamber, the truth the harness's chamber cut: E10.5's chamber categories plus 5% surrounding
cells, configs/pseudo_tasks.yaml):
  t1-blood  blood's share of the chamber cells, the other chamber categories in proportion and the truth's 5%
            surrounding cells holding their share (vec.t1_scope.with_share(..., "blood", x, within=<chamber
            categories>), the T1 build's --blood-share)
  t1-cm     the cardiomyocyte share of heart cells, every category outside the heart holding its share
            (with_share(..., "cm", x), the T1 build's --cm-share)
  heart-scope  T2h/<=E9.5->E10.5: how much tissue surrounds the heart core (vec.prep core: the CM and Endocardium
            lineage groups), as the in-plane scope log(rms_xy full / rms_xy core) the size procedure's scope
            hypotheses use. Truths: E10.5 as released ("own"), and E10.5 with surrounding cells removed until its scope
            is a half or none of its own, either farthest from the core centroid first ("radial", a tighter dissection)
            or at random ("thin"). Predictions: B0 from E9.5 ("hold", the s1 analogue: scope held), B0 shrunk in plane
            by half its scope ("half_size", today's size-only H2 s2), and B0 from E9.5 trimmed to half its scope
            ("half_radial", "half_thin") or to none ("core_radial") before the draw. Every heart prediction is scaled
            in plane about its centroid by the truth's core in-plane RMS over E9.5's (z held, as the pre-registered
            xy_grow_z_hold mode does), an oracle core size, so the rows differ in scope only; half_size then shrinks.
  t1-grid   T1/E9.5->E10.5:chamber, blood x cardiomyocyte shares (GRID2) for truths and predictions alike. Unlike the
            two axes above, every prediction starts from the truth's own chamber mix (the chamber cut's categories
            without its surrounding cells, renormalised: an oracle base) with the two shares set, so a row differs
            from s1's analogue ("own/own") only on the hedged shares, as s1 and s2 do at test. Run in parts
            (--truth-part K/N), merged by summarise.
The T1 axes' prediction mix is the reference's own chamber mix (vec.harness.t1_chamber_shares, B0's) with the axis
share set to each grid value; "own" is B0 itself. The truth's mix is the chamber cut's own with the axis share set to
each grid value; "own" is the cut unchanged. Truth re-sampling uses one fixed seed per grid value (TRUTH_SEED), so a
grid value is one fixed truth across the scoring seeds.

Usage (DGX, repo root, BLAS threads capped; one process per seed range):
  uv run python experiments/s2_hedge/payoff.py t1-blood --seeds 0-4 --out runs/s2_hedge/t1_blood_a
  uv run python experiments/s2_hedge/payoff.py heart-scope --seeds 0-9 --out runs/s2_hedge/heart-scope/0-9
  uv run python experiments/s2_hedge/payoff.py summarise --runs runs/s2_hedge/t1-blood/0-1 runs/s2_hedge/t1-blood/2-3 \
      --out experiments/s2_hedge/results/t1_blood
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from vec import harness, prep, t1_scope
from vec.boards import load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, dense, read_stage
from vec.prep import load_stage

TRUTH_SEED = 20261102  # truth re-sampling (fixed before any run)
TASK = "T1/E9.5->E10.5:chamber"
GRIDS = {  # axis -> (category, grid of shares); "own" is added as the unmodified mix
    "t1-blood": ("blood", (0.03, 0.08, 0.15, 0.35)),
    "t1-cm": ("cm", (0.15, 0.40, 0.58, 0.75)),
}
LABEL = ("sensitivity analysis for the s2 pre-registration, not an adoption result: truths re-sampled to set tissue "
         "mixes, predictions built at the same mixes")


def stratified(categories: np.ndarray, target: dict[str, float], rng: np.random.Generator) -> np.ndarray:
    """Rows re-sampled without replacement to the category shares `target` (summing to 1) at the largest size the
    available cells allow; a category with share 0 is dropped, a category with a share but no cells is refused."""
    cats = {c: np.flatnonzero(categories == c) for c, w in target.items() if w > 0}
    empty = sorted(c for c, ix in cats.items() if len(ix) == 0)
    if empty:
        raise ValueError(f"no cells of {empty} to re-sample")
    total = min(len(ix) / target[c] for c, ix in cats.items())
    rows = [rng.choice(ix, int(round(target[c] * total)), replace=False) for c, ix in cats.items()]
    return np.sort(np.concatenate(rows))


def axis_mix(base: dict[str, float], axis: str, value) -> dict[str, float]:
    """`base` with the axis share set to `value` ("own": unchanged). t1-blood: blood's share of the chamber categories
    present in `base` (as the T1 build's --blood-share), so a truth's surrounding cells keep their share; t1-cm: the
    cardiomyocyte share of the config's heart categories."""
    if value == "own":
        return dict(base)
    cat, _ = GRIDS[axis]
    chamber = set(t1_scope.load_spec().in_scope("chamber"))
    within = [c for c in base if c in chamber] if axis == "t1-blood" else None
    return t1_scope.with_share(base, cat, float(value), within=within)


def chamber_cut(t, data_dir: Path, genes):
    """The harness's chamber cut (vec.harness.filter_target), with each kept row's category."""
    spec = t1_scope.load_spec()
    a = read_stage(data_dir / t.target, genes)
    t1_scope.classify(a, source=harness.stem(t.target))
    cat = np.asarray(a.obs[t1_scope.CATEGORY_KEY]).astype(str)
    f = t.target_filter
    kept = np.flatnonzero(np.isin(cat, sorted(set(spec.in_scope(f["keep"])))))
    pool = np.flatnonzero(np.isin(cat, list(f.get("contamination_from", ()))))
    c = float(f.get("contamination", 0.0))
    k = min(len(pool), int(round(c / (1.0 - c) * len(kept)))) if c > 0 else 0
    extra = np.random.default_rng(int(f["seed"])).choice(pool, k, replace=False) if k else np.array([], int)
    rows = np.sort(np.concatenate([kept, extra]))
    return a[rows], cat[rows]


def run_t1(axis: str, seeds, out: Path, data_dir: Path) -> dict:
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    cut, cut_cat = chamber_cut(t, data_dir, genes)
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes)
    inputs = {r: load_stage(data_dir / f, lineage=False) for r, f in t.inputs.items()}
    role = harness.base_role(t)
    pool, src = harness._shallow(inputs[role]), harness.stem(t.inputs[role])
    start, start_src = harness._shallow(inputs["start"]), harness.stem(t.inputs["start"])
    held = list(t1_scope.held_flat_genes())
    cut_mix = t1_scope.composition(list(cut_cat))
    ref_mix = harness.t1_chamber_shares(pool, src)
    values = ["own", *GRIDS[axis][1]]
    truths = {}
    for q in values:
        mix = axis_mix(cut_mix, axis, q)
        rows = (np.arange(cut.n_obs) if q == "own" else
                stratified(cut_cat, mix, np.random.default_rng([TRUTH_SEED, values.index(q)])))
        truths[str(q)] = (Stage.from_anndata(cut[rows], genes), t1_scope.composition(list(cut_cat[rows])))
    base_board = harness.task_board(t, boards, truths["own"][0].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    n = harness.count_for(cboard, pool.n_obs)  # as the harness's B0 (one count for every prediction)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": held})
        preds = {}
        for p in values:
            mix = axis_mix(ref_mix, axis, p)
            sub, _ = t1_scope.build_submission(pool, cboard, mix, n, ctx.seed, start=start, damping=harness.B0_DAMPING,
                                               held_flat=held, ref_source=src, start_source=start_src)
            preds[str(p)] = harness.finish(ctx, pool, dense(sub.X), None)
        row = {"seed": s, "scores": {}, "skills": {}}
        for q, (truth, _) in truths.items():
            board = harness.task_board(t, boards, truth.n, ref_stage.n)
            f = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, s), truth, ref_stage, seed=s,
                            anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
            row["scores"][q], row["skills"][q] = {"floor": f.score}, {"floor": f.skills}
            for p, pred in preds.items():
                r = score_board(board, pred, truth, ref_stage, seed=s, anchors=f.anchors, build=harness.BUILD,
                                protocol=harness.PROTOCOL)
                row["scores"][q][p], row["skills"][q][p] = r.score, r.skills
        row["seconds"] = round(time.time() - t0, 1)
        rows_out.append(row)
        print(f"{axis} seed {s} ({row['seconds']} s): " + "; ".join(
            f"truth {q}: " + ", ".join(f"{p} {v:.2f}" for p, v in sc.items()) for q, sc in row["scores"].items()),
            flush=True)
    res = {"label": LABEL, "axis": axis, "category": GRIDS[axis][0], "task": t.key, "values": [str(v) for v in values],
           "truth_mix": {q: m for q, (_, m) in truths.items()}, "truth_cells": {q: s.n for q, (s, _) in truths.items()},
           "pred_mix": {str(p): axis_mix(ref_mix, axis, p) for p in values}, "n_pred": n, "truth_seed": TRUTH_SEED,
           "seeds": list(seeds), "rows": rows_out, "provenance": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{axis}.json").write_text(json.dumps(res, indent=1, default=float))
    return res


GRID2 = {"blood": ("own", 0.05, 0.40), "cm": (0.15, "own", 0.45, 0.60, 0.75)}


def grid_points() -> list[tuple[object, object]]:
    return [(b, c) for b in GRID2["blood"] for c in GRID2["cm"]]


def point_label(b, c) -> str:
    return f"b{b}_c{c}"


def grid_mix(base: dict[str, float], b, c) -> dict[str, float]:
    """`base` with blood's share of the chamber categories set to b and the cardiomyocyte share of heart cells to c
    ("own": unchanged); the two settings commute (each keeps the other's ratio)."""
    mix = axis_mix(base, "t1-blood", b)
    return axis_mix(mix, "t1-cm", c)


def run_t1_grid(seeds, part: tuple[int, int], out: Path, data_dir: Path) -> dict:
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[TASK]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    cut, cut_cat = chamber_cut(t, data_dir, genes)
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes)
    inputs = {r: load_stage(data_dir / f, lineage=False) for r, f in t.inputs.items()}
    role = harness.base_role(t)
    pool, src = harness._shallow(inputs[role]), harness.stem(t.inputs[role])
    start, start_src = harness._shallow(inputs["start"]), harness.stem(t.inputs["start"])
    held = list(t1_scope.held_flat_genes())
    cut_mix = t1_scope.composition(list(cut_cat))
    chamber = set(t1_scope.load_spec().in_scope("chamber"))
    oracle = {k: v for k, v in cut_mix.items() if k in chamber and v > 0}
    oracle = {k: v / sum(oracle.values()) for k, v in oracle.items()}
    points = grid_points()
    mine = [pt for i, pt in enumerate(points) if i % part[1] == part[0]]
    truths = {}
    for b, c in mine:
        mix = grid_mix(cut_mix, b, c)
        i = points.index((b, c))
        rows = (np.arange(cut.n_obs) if (b, c) == ("own", "own") else
                stratified(cut_cat, mix, np.random.default_rng([TRUTH_SEED, 1000 + i])))
        truths[point_label(b, c)] = (Stage.from_anndata(cut[rows], genes), t1_scope.composition(list(cut_cat[rows])))
    base_board = harness.task_board(t, boards, cut.n_obs, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    n = harness.count_for(cboard, pool.n_obs)
    rows_out = []
    for s_ in seeds:
        t0 = time.time()
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s_,
                                   derived={"t1_held_flat": held})
        row = {"seed": s_, "scores": {}, "skills": {}}
        floors = {}  # truth -> (board, floor anchors), computed before any prediction
        for q, (truth, _) in truths.items():
            board = harness.task_board(t, boards, truth.n, ref_stage.n)
            fl = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, s_), truth, ref_stage, seed=s_,
                             anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
            row["scores"][q], row["skills"][q] = {"floor": fl.score}, {"floor": fl.skills}
            floors[q] = (board, fl.anchors)
        # One prediction in memory at a time: holding all 15 took about 20 GB per process (2026-11-01). Every builder
        # and scorer call is seeded on its own (no global RNG), so the order does not change any score.
        for b, c in points:
            sub, _ = t1_scope.build_submission(pool, cboard, grid_mix(oracle, b, c), n, ctx.seed, start=start,
                                               damping=harness.B0_DAMPING, held_flat=held, ref_source=src,
                                               start_source=start_src)
            pred, p_ = harness.finish(ctx, pool, dense(sub.X), None), point_label(b, c)
            del sub
            for q, (truth, _) in truths.items():
                board, anchors = floors[q]
                r = score_board(board, pred, truth, ref_stage, seed=s_, anchors=anchors, build=harness.BUILD,
                                protocol=harness.PROTOCOL)
                row["scores"][q][p_], row["skills"][q][p_] = r.score, r.skills
            del pred
        row["seconds"] = round(time.time() - t0, 1)
        rows_out.append(row)
        print(f"t1-grid part {part[0]}/{part[1]} seed {s_} ({row['seconds']} s)", flush=True)
    labels = [point_label(b, c) for b, c in points]
    res = {"label": LABEL, "axis": "t1-grid", "category": "blood x cm", "task": t.key,
           "pred_values": labels, "truth_values": labels, "part": list(part),
           "truth_mix": {q: m for q, (_, m) in truths.items()},
           "truth_cells": {q: st.n for q, (st, _) in truths.items()},
           "pred_mix": {point_label(b, c): grid_mix(oracle, b, c) for b, c in points}, "n_pred": n,
           "truth_seed": TRUTH_SEED, "seeds": list(seeds), "rows": rows_out, "provenance": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    (out / "t1-grid.json").write_text(json.dumps(res, indent=1, default=float))
    return res


HEART_TASK = "T2h/<=E9.5->E10.5"
HEART_PREDS = ("hold", "half_size", "half_radial", "half_thin", "core_radial")
HEART_TRUTHS = ("own", "radial_0.5", "radial_0", "thin_0.5", "thin_0")
# PLAN 4.3 step 5's second released transition, defined here as a diagnostic (adoption false; not in the harness
# config). E8.25_late predicts E8.5 on a synthetic board that borrows the heart interpolation board's gene panel. No
# earlier stage is an input, so the predictions are copies (no momentum); step 5 reads only d2 and Dice, which the
# coordinates set.
HEART_EXTRA_TASKS = {
    "T2h/E8.25_late->E8.5": {"key": "T2h/E8.25_late->E8.5", "test_board": "T2:heart:test_extrap", "task": "T2",
                             "kind": "extrap", "setting": "heart", "target": "T2/E8.5.h5ad",
                             "reference": "T2/E8.25_late.h5ad", "inputs": {"ref": "T2/E8.25_late.h5ad"},
                             "twin_board": "T2:heart:val_interp", "adoption": False},
}


def heart_task(key: str):
    """A heart pseudo-task from the harness config, or one of HEART_EXTRA_TASKS (checked by the harness's rules)."""
    tasks = {x.key: x for x in harness.load_tasks()["T2:heart:test_extrap"]}
    if key in tasks:
        return tasks[key]
    if key in HEART_EXTRA_TASKS:
        return harness.check_task(harness.PseudoTask(**HEART_EXTRA_TASKS[key]))
    raise SystemExit(f"unknown heart task {key!r}; one of {sorted(tasks) + sorted(HEART_EXTRA_TASKS)}")


def xy_rms(C: np.ndarray) -> float:
    X = C[:, :2] - C[:, :2].mean(0)
    return float(np.sqrt((X ** 2).sum(1).mean()))


def xy_scope(C: np.ndarray, core: np.ndarray) -> float:
    """log(in-plane RMS of the whole cloud / in-plane RMS of its core), the size procedure's in-plane scope."""
    return float(np.log(xy_rms(C) / xy_rms(C[core])))


def trim_mask(C: np.ndarray, core: np.ndarray, frac: float, mode: str, rng: np.random.Generator) -> np.ndarray:
    """Cells kept when surrounding (non-core) cells are removed, farthest in plane from the core centroid first
    ("radial") or in random order ("thin"), until the in-plane scope is at most `frac` of its own (the fewest removed
    that reach it, by bisection on the removal order); core cells are always kept."""
    target = frac * xy_scope(C, core)
    shell = np.flatnonzero(~core)
    if mode == "radial":
        d = np.linalg.norm(C[shell, :2] - C[core, :2].mean(0), axis=1)
        order = shell[np.argsort(-d, kind="stable")]
    elif mode == "thin":
        order = rng.permutation(shell)
    else:
        raise ValueError(mode)

    def keep_after(k: int) -> np.ndarray:
        keep = np.ones(len(C), dtype=bool)
        keep[order[:k]] = False
        return keep

    lo, hi = 0, len(order)
    while lo < hi:
        mid = (lo + hi) // 2
        m = keep_after(mid)
        if xy_scope(C[m], core[m]) <= target:
            hi = mid
        else:
            lo = mid + 1
    return keep_after(lo)


def run_heart(seeds, out: Path, data_dir: Path, task: str = HEART_TASK) -> dict:
    t = heart_task(task)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    truths = {}
    for i, name in enumerate(HEART_TRUTHS):
        if name == "own":
            keep = np.ones(ta.n_obs, dtype=bool)
        else:
            mode, frac = name.split("_")
            keep = trim_mask(TC, tcore, float(frac), mode, np.random.default_rng([TRUTH_SEED, i]))
        truths[name] = (Stage.from_anndata(ta[keep], genes, need_coords=True),
                        {"cells": int(keep.sum()), "core_share": float(tcore[keep].mean()),
                         "xy_scope": xy_scope(TC[keep], tcore[keep])})
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    pool = inputs["ref"]
    pcore = np.asarray(pool.obs["core"], dtype=bool)
    PC = np.asarray(pool.obsm[COORD_KEY], dtype=float)[:, :3]
    s_pool = xy_scope(PC, pcore)
    f_core = xy_rms(TC[tcore]) / xy_rms(PC[pcore])  # oracle core size, in plane (z held)
    trims = {"hold": None, "half_size": None, "half_radial": ("radial", 0.5), "half_thin": ("thin", 0.5),
             "core_radial": ("radial", 0.0)}
    pools, pred_info = {}, {}
    for i, name in enumerate(HEART_PREDS):
        tr = trims[name]
        keep = (np.ones(pool.n_obs, dtype=bool) if tr is None else
                trim_mask(PC, pcore, tr[1], tr[0], np.random.default_rng([TRUTH_SEED, 100 + i])))
        pools[name] = pool if tr is None else pool[keep].copy()
        pred_info[name] = {"pool_cells": int(keep.sum()), "core_share": float(pcore[keep].mean()),
                           "xy_scope": xy_scope(PC[keep], pcore[keep])}
    base_board = harness.task_board(t, boards, truths["own"][0].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    builder = "b0" if "start" in inputs else "copy"  # B0's momentum needs an earlier stage
    rows_out = []
    for s in seeds:
        t0 = time.time()
        preds = {}
        for name in HEART_PREDS:
            ctx = harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pools[name]},
                                       seed=harness.BUILD_SEED_OFFSET + s)
            st = harness._b0_t2(ctx) if builder == "b0" else harness._copy(ctx)
            C = np.asarray(st.coords, dtype=float).copy()
            c = C[:, :2].mean(0)
            f = f_core * (np.exp(-0.5 * s_pool) if name == "half_size" else 1.0)
            C[:, :2] = c + (C[:, :2] - c) * f
            preds[name] = Stage(X=st.X, coords=C.astype(np.float32), celltype=st.celltype)
        row = {"seed": s, "scores": {}, "skills": {}}
        for q, (truth, _) in truths.items():
            # the validation geometry's board (its published counts are the untrimmed files'), with the trimmed
            # truth's own count: scoring reads the board's task and, for diagnostics only, truth_cells
            board = replace(base_board, truth_cells=truth.n)
            fl = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, s), truth, ref_stage, seed=s,
                             anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
            row["scores"][q], row["skills"][q] = {"floor": fl.score}, {"floor": fl.skills}
            for p, pred in preds.items():
                r = score_board(board, pred, truth, ref_stage, seed=s, anchors=fl.anchors, build=harness.BUILD,
                                protocol=harness.PROTOCOL)
                row["scores"][q][p], row["skills"][q][p] = r.score, r.skills
        row["seconds"] = round(time.time() - t0, 1)
        rows_out.append(row)
        print(f"heart-scope seed {s} ({row['seconds']} s): " + "; ".join(
            f"truth {q}: " + ", ".join(f"{p} {v:.2f}" for p, v in sc.items()) for q, sc in row["scores"].items()),
            flush=True)
    res = {"label": LABEL, "axis": "heart-scope", "category": "scope", "task": t.key, "builder": builder,
           "pred_values": list(HEART_PREDS), "truth_values": list(HEART_TRUTHS),
           "truth_info": {q: i for q, (_, i) in truths.items()}, "pred_info": pred_info,
           "pool_xy_scope": s_pool, "oracle_core_xy_factor": f_core, "truth_seed": TRUTH_SEED,
           "seeds": list(seeds), "rows": rows_out, "provenance": provenance()}
    out.mkdir(parents=True, exist_ok=True)
    (out / "heart-scope.json").write_text(json.dumps(res, indent=1, default=float))
    return res


def provenance() -> dict:
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                           text=True).stdout.strip()
    return {"commit": commit, "tree_clean": dirty == "", "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}


def summarise(runs: list[Path], out: Path) -> dict:
    parts = [json.loads(p.read_text()) for d in runs for p in sorted(Path(d).glob("*.json"))]
    axes, tasks = {p["axis"] for p in parts}, {p["task"] for p in parts}
    if len(axes) != 1 or len(tasks) != 1:
        raise SystemExit(f"runs mix axes {sorted(axes)} or tasks {sorted(tasks)}")
    commits = {p["provenance"]["commit"] for p in parts}
    if len(commits) != 1 or not all(p["provenance"]["tree_clean"] for p in parts):
        raise SystemExit(f"runs come from different commits or a dirty tree: {sorted(commits)}")
    merged: dict[int, dict] = {}
    for part_ in parts:
        for r in part_["rows"]:
            m = merged.setdefault(r["seed"], {"seed": r["seed"], "scores": {}, "skills": {}})
            twice = set(m["scores"]) & set(r["scores"])
            if twice:
                raise SystemExit(f"seed {r['seed']}: truth(s) {sorted(twice)} appear twice")
            m["scores"].update(r["scores"])
            m["skills"].update(r["skills"])
    rows = [merged[k] for k in sorted(merged)]
    seeds = [r["seed"] for r in rows]
    head = dict(parts[0])
    for k in ("truth_mix", "truth_cells", "truth_info"):
        if k in head:
            head[k] = {q: v for part_ in parts for q, v in part_[k].items()}
    want = head.get("truth_values", head.get("values"))
    short = [r["seed"] for r in rows if set(r["scores"]) != set(want)]
    if short:
        raise SystemExit(f"seeds {short} miss some truths")
    pv = head.get("pred_values", head.get("values"))
    tv = head.get("truth_values", head.get("values"))
    mean = {q: {p: float(np.mean([r["scores"][q][p] for r in rows])) for p in ["floor", *pv]} for q in tv}
    sd = {q: {p: float(np.std([r["scores"][q][p] for r in rows], ddof=1)) if len(rows) > 1 else float("nan")
              for p in pv} for q in tv}
    # each prediction's gain over the s1 analogue (the first row), paired over seeds
    s1 = pv[0]
    gain = {q: {p: {"mean": float(np.mean(d)), "ahead": int((d > 0).sum()), "of": len(d)}
                for p in pv[1:] for d in [np.array([r["scores"][q][p] - r["scores"][q][s1] for r in rows])]}
            for q in tv}
    keep = [k for k in ("label", "axis", "category", "task", "values", "pred_values", "truth_values", "truth_mix",
                        "truth_cells", "truth_info", "pred_mix", "pred_info", "n_pred", "pool_xy_scope",
                        "oracle_core_xy_factor", "truth_seed") if k in head]
    res = {k: head[k] for k in keep}
    res.update(seeds=seeds, mean=mean, sd=sd, gain_over_first=gain, provenance=head["provenance"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "payoff.json").write_text(json.dumps(res, indent=1, default=float))
    cat = head["category"]
    lines = [f"{head['axis']} on {head['task']} ({LABEL}); {len(seeds)} seeds; commit {head['provenance']['commit']}",
             f"rows: prediction ({cat}); columns: truth ({cat}); cells: mean board score (floor row last); "
             f"[gain over {s1}, seeds ahead]",
             "pred \\ truth " + " ".join(f"{q:>17s}" for q in tv)]
    for p in pv:
        cells = []
        for q in tv:
            g = gain[q].get(p)
            cells.append(f"{mean[q][p]:6.2f}" + (f" [{g['mean']:+5.2f} {g['ahead']:>2d}]" if g else " " * 11))
        lines.append(f"{p:>12s} " + " ".join(f"{c:>17s}" for c in cells))
    lines.append(f"{'floor':>12s} " + " ".join(f"{mean[q]['floor']:17.2f}" for q in tv))
    if "truth_mix" in head:
        heart = t1_scope.load_spec().heart
        lines.append("truth mixes: " + "; ".join(
            f"{q}: blood {m.get('blood', 0.0):.3f}, cm/heart {t1_scope.share_of(m, 'cm', heart):.3f} "
            f"({head['truth_cells'][q]} cells)" for q in tv for m in [head["truth_mix"][q]]))
    if "truth_info" in head:
        lines.append("truths: " + "; ".join(
            f"{q}: {i['cells']} cells, core {i['core_share']:.2f}, xy scope {i['xy_scope']:.3f}"
            for q, i in head["truth_info"].items()))
        lines.append("predictions (pool): " + "; ".join(
            f"{p}: {i['pool_cells']} cells, core {i['core_share']:.2f}, xy scope {i['xy_scope']:.3f}"
            for p, i in head["pred_info"].items()))
    (out / "payoff.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return res


def hedge(runs: list[Path], s1: str, out: Path | None = None) -> dict:
    """For each candidate s2: the expected better of (s1, s2) over the truth grid, paired by seed (the captain marks
    the better-scoring upload), under two stated priors: uniform over the truths, and one third on the truth that
    equals s1's mix ("own" or "own/own") with the rest uniform. Also each candidate's worst truth. Informs the s2
    pre-registration; it chooses nothing."""
    parts = [json.loads(q.read_text()) for d in runs for q in sorted(Path(d).glob("*.json"))]
    merged: dict[int, dict] = {}
    for part_ in parts:
        for r in part_["rows"]:
            merged.setdefault(r["seed"], {}).update(r["scores"])
    head = parts[0]
    pv = head.get("pred_values", head.get("values"))
    tv = head.get("truth_values", head.get("values"))
    if s1 not in pv:
        raise SystemExit(f"{s1!r} is not a prediction; one of {pv}")
    seeds = sorted(merged)
    if any(set(merged[k]) != set(tv) for k in seeds):
        raise SystemExit("some seeds miss some truths")
    own = next((q for q in tv if q in ("own", "bown_cown")), None)
    priors = {"uniform": {q: 1 / len(tv) for q in tv}}
    if own is not None:
        rest = (2 / 3) / (len(tv) - 1)
        priors["own_third"] = {q: (1 / 3 if q == own else rest) for q in tv}
    res = {"s1": s1, "seeds": seeds, "priors": priors, "candidates": {}}
    for p2 in pv:
        per_seed = {name: [sum(w[q] * max(merged[k][q][s1], merged[k][q][p2]) for q in tv) for k in seeds]
                    for name, w in priors.items()}
        worst = min(tv, key=lambda q: np.mean([max(merged[k][q][s1], merged[k][q][p2]) for k in seeds]))
        res["candidates"][p2] = {
            **{name: {"mean": float(np.mean(v)), "sd": float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")}
               for name, v in per_seed.items()},
            "worst_truth": worst,
            "worst_mean": float(np.mean([max(merged[k][worst][s1], merged[k][worst][p2]) for k in seeds]))}
    order = sorted(pv, key=lambda c: -res["candidates"][c]["uniform"]["mean"])
    lines = [f"{head['axis']}: expected better of s1={s1} and s2 over {len(tv)} truths, {len(seeds)} seeds "
             "(s2 = s1 is s1 alone)", f"{'s2':>14s} " + " ".join(f"{n:>16s}" for n in priors) + "   worst truth"]
    for c in order:
        r = res["candidates"][c]
        lines.append(f"{c:>14s} " + " ".join(f"{r[n]['mean']:8.2f} +/- {r[n]['sd']:4.2f}" for n in priors)
                     + f"   {r['worst_truth']} ({r['worst_mean']:.2f})")
    print("\n".join(lines))
    if out is not None:
        out.mkdir(parents=True, exist_ok=True)
        (out / f"hedge_{s1.replace('/', '_')}.json").write_text(json.dumps(res, indent=1, default=float))
        (out / f"hedge_{s1.replace('/', '_')}.txt").write_text("\n".join(lines) + "\n")
    return res


STEP5_D2_GAIN = 0.15  # docs/PLAN.md 4.3 step 5: >= 0.15 d2 skill, Dice not lowered, on both released transitions


def step5(runs: list[Path], out: Path | None = None) -> dict:
    """PLAN 4.3 step 5's bar for each scope prediction against hold on the released target (truth "own"): the mean
    gain in d2_shape skill (at least 0.15) and in occupancy_dice skill (not below 0), paired by seed."""
    parts = [json.loads(q.read_text()) for d in runs for q in sorted(Path(d).glob("heart-scope.json"))]
    tasks = {p["task"] for p in parts}
    if len(tasks) != 1:
        raise SystemExit(f"runs mix tasks {sorted(tasks)}")
    rows = sorted((r for p in parts for r in p["rows"]), key=lambda r: r["seed"])
    sk = [r["skills"]["own"] for r in rows]
    res = {"task": tasks.pop(), "seeds": [r["seed"] for r in rows], "d2_gain_bar": STEP5_D2_GAIN, "preds": {}}
    lines = [f"PLAN 4.3 step 5 on {res['task']} (the target as released; {len(rows)} seeds): mean gain over hold in "
             "skill (seeds ahead)"]
    for p in parts[0]["pred_values"][1:]:
        d2 = np.array([s[p]["d2_shape"] - s["hold"]["d2_shape"] for s in sk])
        dice = np.array([s[p]["occupancy_dice"] - s["hold"]["occupancy_dice"] for s in sk])
        ok = bool(d2.mean() >= STEP5_D2_GAIN and dice.mean() >= 0)
        res["preds"][p] = {"d2_gain": float(d2.mean()), "d2_ahead": int((d2 > 0).sum()),
                           "dice_gain": float(dice.mean()), "dice_ahead": int((dice > 0).sum()), "passes": ok}
        lines.append(f"{p:>12s}  d2 {d2.mean():+.3f} ({(d2 > 0).sum()}/{len(d2)})  dice {dice.mean():+.3f} "
                     f"({(dice > 0).sum()}/{len(dice)})  {'passes' if ok else 'fails'}")
    print("\n".join(lines))
    if out is not None:
        out.mkdir(parents=True, exist_ok=True)
        (out / "step5.json").write_text(json.dumps(res, indent=1, default=float))
        (out / "step5.txt").write_text("\n".join(lines) + "\n")
    return res


def parse_seeds(text: str) -> tuple[int, ...]:
    lo, _, hi = text.partition("-")
    return tuple(range(int(lo), int(hi) + 1)) if hi else tuple(int(x) for x in text.split(","))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for axis in GRIDS:
        p = sub.add_parser(axis)
        p.add_argument("--seeds", default="0-9")
        p.add_argument("--out", required=True)
        p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("t1-grid")
    p.add_argument("--seeds", default="0-9")
    p.add_argument("--truth-part", default="0/1", help="K/N: this process scores truths K, K+N, K+2N, ...")
    p.add_argument("--out", required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("heart-scope")
    p.add_argument("--seeds", default="0-9")
    p.add_argument("--task", default=HEART_TASK, help="a heart pseudo-task key, or one of HEART_EXTRA_TASKS")
    p.add_argument("--out", required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("step5")
    p.add_argument("--runs", nargs="+", required=True)
    p.add_argument("--out")
    p = sub.add_parser("hedge")
    p.add_argument("--runs", nargs="+", required=True)
    p.add_argument("--s1", required=True, help="the prediction label that stands for s1 (e.g. own, bown_cown, hold)")
    p.add_argument("--out")
    p = sub.add_parser("summarise")
    p.add_argument("--runs", nargs="+", required=True)
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "summarise":
        summarise([Path(r) for r in a.runs], Path(a.out))
        return
    if a.cmd == "hedge":
        hedge([Path(r) for r in a.runs], a.s1, Path(a.out) if a.out else None)
        return
    if a.cmd == "t1-grid":
        k, _, n_ = a.truth_part.partition("/")
        run_t1_grid(parse_seeds(a.seeds), (int(k), int(n_)), Path(a.out), Path(a.data))
        return
    if a.cmd == "heart-scope":
        run_heart(parse_seeds(a.seeds), Path(a.out), Path(a.data), a.task)
        return
    if a.cmd == "step5":
        step5([Path(r) for r in a.runs], Path(a.out) if a.out else None)
        return
    run_t1(a.cmd, parse_seeds(a.seeds), Path(a.out), Path(a.data))


if __name__ == "__main__":
    sys.exit(main())
