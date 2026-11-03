#!/usr/bin/env python
"""Embryo diagnostic (2026-11-03; not an adoption result, released stages only): each embryo stage's detection depth,
and mmd_u / variogram of plain copies and bracket mixtures against each embryo backtest target, as the harness scores
them (the bake protocol's truth half). It found the depth gap behind tries E7 and E8 of experiments/harness/TRIES.md.

Usage (DGX, repo root, BLAS threads capped as docs/COMPUTE.md says):
  uv run python experiments/harness/diag_emb_mix.py --out experiments/harness/results/diag_emb_mix
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage, panel_metrics, subsample, veckit_halves
from vec.prep import load_stage
from vec.t2_build import library_sizes

DATA = REPO_ROOT / "data" / "raw" / "T2"
STAGES = ("E6.75", "E7.25", "E7.5", "E8.0")
TASKS = (("E7.25", "E6.75", "E7.5", 2 / 3), ("E7.5", "E7.25", "E8.0", 1 / 3))
N_PRED = 5000


def _mix(pb: Stage, pa: Stage, t: float, n: int, seed: int) -> Stage:
    rng = np.random.default_rng(seed)
    nb = int(round((1 - t) * n))
    ib = rng.choice(pb.n, nb, replace=pb.n < nb)
    ia = rng.choice(pa.n, n - nb, replace=pa.n < n - nb)
    return Stage(X=np.vstack([pb.X[ib], pa.X[ia]]), coords=np.vstack([pb.coords[ib], pa.coords[ia]]))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    genes = list(load_boards()["T2:embryo:val_interp"].genes)
    loaded = {s: load_stage(DATA / f"{s}.h5ad", genes=genes, lineage=False, drop_duplicates=True) for s in STAGES}
    depth = {}
    for s, ad_ in loaded.items():
        lib = library_sizes(ad_.X, ad_.layers["counts"])
        det = np.asarray((ad_.layers["counts"] > 0).sum(1)).ravel()
        depth[s] = {"cells": int(ad_.n_obs), "median_counts_per_cell": float(np.median(lib)),
                    "q10_q90_counts": [float(np.percentile(lib, 10)), float(np.percentile(lib, 90))],
                    "median_panel_genes_detected": float(np.median(det))}
    mixes = {}
    for target, before, after, t in TASKS:
        full = Stage.from_anndata(load_stage(DATA / f"{target}.h5ad", genes=genes, lineage=False), genes,
                                  need_coords=True)
        truth = veckit_halves(subsample(full, 0.1, 0), 0)[0]
        pb, pa = (Stage.from_anndata(loaded[s], genes, need_coords=True) for s in (before, after))
        ref = subsample(pb, 0.1, 0)
        preds = {"copy_before": subsample(pb, min(1.0, N_PRED / pb.n), 1),
                 "copy_after": subsample(pa, min(1.0, N_PRED / pa.n), 1),
                 "mix_nominal_t": _mix(pb, pa, t, N_PRED, 1), "mix_half": _mix(pb, pa, 0.5, N_PRED, 1)}
        mixes[f"{target} from {before}+{after} (t = {t:.3f})"] = {
            k: {m: float(v) for m, v in panel_metrics("T2", p, truth, ref, seed=0).items()
                if m in ("mmd_u", "variogram", "de_score", "de_direction", "neighborhood_mmd")}
            for k, p in preds.items()}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=REPO_ROOT).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                                text=True, cwd=REPO_ROOT).stdout.strip())
    res = {"label": "diagnostic, not an adoption result", "commit": commit, "tree_clean": not dirty,
           "depth": depth, "mixtures": mixes}
    (out / "diag_emb_mix.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
