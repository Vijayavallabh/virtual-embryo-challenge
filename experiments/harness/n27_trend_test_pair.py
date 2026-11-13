#!/usr/bin/env python
"""Diagnostic N27 of experiments/harness/TRIES.md (not a try): the trend factor of heart candidate b0_trend at the
heart's TEST pair, against B0's own factor there. Target-free: it reads the five released heart stages the test task
has as inputs (E10.5 only as the test pair's reference stage, under A20); nothing is compared with a truth, nothing is
built and nothing is scored. It prints and writes lineage names, stage roles, weights, cosines, size ratios and mean
absolute log factors only: no cell number of any stage, no board limit, no anchor.

For each pool (the first upload's: all of E10.5; the second upload's: E10.5 trimmed as build_test.py heart --entry
s2_trim trims it), per lineage B0 shifts at the test pair: the stages usable for the trend, the transition weights of
the least-squares slope, the cosine and L2-norm ratio of the trend's log factor against B0's own (half the E9.5 ->
E10.5 step), and the mean absolute log factor of each (vec.harness.b0_trend_factors, B0_TREND_RULE).

Run (branch campaign/trend): python -I -B experiments/harness/n27_trend_test_pair.py --out <folder>
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np

from vec import harness, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.io import COORD_KEY
from vec.prep import load_stage

BOARD = "T2:heart:test_extrap"
ROLES = {"earliest": "T2/E8.25_late.h5ad", "earlier2": "T2/E8.5.h5ad", "earlier": "T2/E8.75.h5ad",
         "start": "T2/E9.5.h5ad", "ref": "T2/E10.5.h5ad"}
DROP = ("n_cells",)  # per-lineage cell numbers of the released stages: not written (the row's promise)


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=False).stdout.strip()


def _view(board_key: str) -> harness.TaskView:
    days = {r: float(harness.STAGE_DAYS[f]) for r, f in ROLES.items()}
    return harness.TaskView(key="test:T2h", task="T2", kind="extrap", setting="heart", test_board=board_key,
                            roles=tuple(ROLES), role_days=days,
                            role_stems={r: harness.stem(f) for r, f in ROLES.items()}, reference_role="ref",
                            target_day=float("nan"))  # never read by the trend factor


def _pool_record(view, board, inputs: dict, pool) -> dict:
    ctx = harness.BuildContext(task=view, board=harness._candidate_board(board), inputs={**inputs, "ref": pool},
                               seed=harness.BUILD_SEED_OFFSET)
    genes = list(pool.var_names)
    types = sorted(set(np.asarray(pool.obs["merged"]).astype(str)))
    _, rec = harness.b0_trend_factors(ctx, genes, types)
    lin = []
    for row in rec["lineages"]:
        row = {k: v for k, v in row.items() if k not in DROP}
        row["usable"] = {r: (r in row["stages"]) for r in rec["stages"]}
        lin.append(row)
    out = {k: v for k, v in rec.items() if k != "lineages"}
    out["lineages"] = lin
    out["n_genes"] = len(genes)
    out["median_mean_abs_logf_b0"] = float(np.median([r["mean_abs_logf_b0"] for r in lin])) if lin else float("nan")
    out["median_mean_abs_logf_trend"] = (float(np.median([r["mean_abs_logf_trend"] for r in lin])) if lin
                                         else float("nan"))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    out.mkdir(parents=True, exist_ok=True)
    board = load_boards()[BOARD]
    genes = board.genes
    view = _view(board.key)
    inputs = {r: load_stage(harness.DATA_DIR / f, genes=genes, drop_duplicates=True) for r, f in ROLES.items()}
    ref = inputs["ref"]
    core = np.asarray(ref.obs["core"], dtype=bool)
    C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
    keep = t2_build.radial_trim_mask(C, core, 0.5)  # as build_test.py heart --entry s2_trim
    res = {"label": "diagnostic N27 (not a try): b0_trend's factor at the heart's test pair against B0's own; "
                    "target-free, nothing built or scored, no cell number written",
           "board": board.key, "stages": dict(ROLES),
           "rules": {"trend": harness.B0_TREND_RULE},
           "pools": {"s1": _pool_record(view, board, inputs, ref),
                     "s2_trim": _pool_record(view, board, inputs, ref[keep].copy())},
           "provenance": {"commit": _git("rev-parse", "HEAD"),
                          "tree_clean": _git("status", "--porcelain", "--untracked-files=no") == "",
                          "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                          "files_sha256": {f: harness._sha256(harness.DATA_DIR / f) for f in ROLES.values()}}}
    (out / "n27.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    lines = [res["label"], f"commit {res['provenance']['commit'][:7]}, tree clean {res['provenance']['tree_clean']}"]
    for name, p in res["pools"].items():
        lines.append(f"pool {name}: stages {p['stages']} at days {p['days']}, budget {p['budget_days']} days, "
                     f"{p['n_trend']} lineage(s) with a trend factor; median cosine {p['median_cosine']:.3f}, median "
                     f"size against B0's {p['median_norm_ratio']:.3f}; median mean |log f|: B0 "
                     f"{p['median_mean_abs_logf_b0']:.3f}, trend {p['median_mean_abs_logf_trend']:.3f}")
        for r in p["lineages"]:
            w = ", ".join(f"{x:.3f}" for x in r["weights"])
            lines.append(f"    {r['lineage']:<20} stages {'/'.join(r['stages'])}; weights {w}; cosine "
                         f"{r['cosine']:+.3f}; size {r['norm_ratio']:.3f}; mean |log f| B0 "
                         f"{r['mean_abs_logf_b0']:.3f}, trend {r['mean_abs_logf_trend']:.3f}"
                         + ("; keeps B0's own factor" if r["own_factor"] else ""))
    (out / "n27.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
