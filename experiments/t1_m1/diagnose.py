#!/usr/bin/env python
"""M1 diagnostic on the T1 E9.5 -> E10.5 pseudo-tasks: B0 against B0 with Feng's direction (M1's direction step).

NOT an adoption result. The harness refuses Feng on every T1 pseudo-task: Feng's E9.5 cells lie inside this task's
version of the test's external window, (E8.5, E11.5] (vec.harness.external_window), while at test they lie outside
(E9.5, E13.5]. The geometry is the pre-registered readout's (experiments/t1_direction, Feng E9.5 -> E14.5 on the
E9.5 -> E10.5 transition), which set M1's direction source. This run checks the build end to end on the board
emulator (count-space edits on real cells, B0's composition, the held-flat list), with every setting fixed before it
ran: the direction is the readout's primary predictor, its size the time fraction of Feng's bracket the step covers
(1/5 here, 2/5 at test), and cells without a Feng group keep B0's damped momentum. Nothing is chosen on its result.

Scoring is the harness's: bake protocol, scorer build main, 10% working copies, each seed's floor and ceiling shared by
the candidates, build seed = 1,000,000 + seed, so B0 here is the harness's B0 cell for cell.

Usage (DGX, repo root, BLAS threads capped; ~15 min per pseudo-task):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 uv run python experiments/t1_m1/diagnose.py \
      --task "T1/E9.5->E10.5:chamber" --out experiments/t1_m1/results/diag_v1
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec.adoption import SEEDS
from vec.boards import load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import feng_dirs  # noqa: E402

LABEL = ("diagnostic, not an adoption result: Feng's E9.5 cells lie inside this pseudo-task's external window "
         "(the harness refuses them here); geometry of the pre-registered T1 direction readout")


def m1_build(ctx: harness.BuildContext, directions: dict, groups: dict, fraction: float):
    """B0 (the harness's `_b0_t1(B0_DAMPING)`) with Feng's direction on the cells of Feng's groups."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    shares = harness.t1_chamber_shares(ref, src)
    start = harness._shallow(ctx.inputs["start"])
    held = list(ctx.derived["t1_held_flat"])
    sub, rep = t1_scope.build_submission(ref, ctx.board, shares, harness.count_for(ctx.board, ref.n_obs), ctx.seed,
                                         start=start, damping=harness.B0_DAMPING, held_flat=held, ref_source=src,
                                         start_source=ctx.task.role_stems["start"], directions=directions,
                                         direction_fraction=fraction, direction_groups=groups)
    return harness.finish(ctx, ref, dense(sub.X), None), rep


def run(task_key: str, out: Path, seeds, data_dir: Path) -> dict:
    tasks = {t.key: t for t in harness.load_tasks()["T1:test"]}
    t = tasks[task_key]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, filt = (harness.filter_target(t, data_dir, genes) if t.target_filter is not None else
                    (Stage.from_anndata(read_stage(data_dir / t.target, genes), genes), None))
    ref = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref.n)
    inputs = {r: load_stage(data_dir / f, lineage=False) for r, f in t.inputs.items()}
    pool = inputs[harness.base_role(t)]
    directions, ext_record = feng_dirs.feng_directions(list(pool.var_names))
    groups = feng_dirs.label_groups(pool.obs["celltype"].astype(str).unique())
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    held = t1_scope.held_flat_genes()
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    b0 = harness._b0_t1(harness.B0_DAMPING)
    rows, ext_rows = [], None
    for s in seeds:
        t0 = time.time()
        f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, s), target, ref, seed=s, anchors="computed",
                        build=harness.BUILD, protocol=harness.PROTOCOL)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": held})
        preds = {"b0": b0(ctx)}
        preds["m1_feng"], rep = m1_build(ctx, directions, groups, fraction)
        ext_rows = ext_rows or rep["external"]["rows"]
        row = {"seed": s, "scores": {"floor": f.score}, "skills": {"floor": f.skills}, "metrics": {"floor": f.metrics}}
        for name, p in preds.items():
            r = score_board(board, p, target, ref, seed=s, anchors=f.anchors, build=harness.BUILD,
                            protocol=harness.PROTOCOL)
            row["scores"][name], row["skills"][name], row["metrics"][name] = r.score, r.skills, r.metrics
        row["seconds"] = round(time.time() - t0, 1)
        rows.append(row)
        print(f"{t.key} seed {s}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["scores"].items()), flush=True)
    d = np.array([r["scores"]["m1_feng"] - r["scores"]["b0"] for r in rows])
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    clean = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                           text=True).stdout.strip() == ""
    res = {"label": LABEL, "task": t.key, "board": board.key, "target_filter_record": filt, "fraction": fraction,
           "groups": {lab: g for lab, g in groups.items() if g in directions}, "external": ext_record,
           "external_rows": ext_rows, "seeds": list(seeds), "rows": rows,
           "summary": {"m1_feng_minus_b0": {"mean": float(d.mean()), "ahead": int((d > 0).sum()), "of": len(d)},
                       "mean_scores": {k: float(np.mean([r["scores"][k] for r in rows])) for k in rows[0]["scores"]}},
           "provenance": {"commit": commit, "tree_clean": clean,
                          "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    path = out / (re.sub(r"[^A-Za-z0-9.]+", "_", t.key).strip("_") + ".json")
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}: m1_feng - b0 = {d.mean():+.2f}, ahead in {(d > 0).sum()} of {len(d)}", flush=True)
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--task", required=True, choices=["T1/E9.5->E10.5", "T1/E9.5->E10.5:chamber"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", default="0-9")
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    lo, _, hi = a.seeds.partition("-")
    seeds = tuple(range(int(lo), int(hi) + 1)) if hi else tuple(int(x) for x in a.seeds.split(","))
    if seeds != tuple(SEEDS)[: len(seeds)]:
        print(f"note: seeds {seeds} differ from the harness's {tuple(SEEDS)}", flush=True)
    run(a.task, Path(a.out), seeds, Path(a.data))


if __name__ == "__main__":
    main()
