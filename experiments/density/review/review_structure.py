"""Reviewer's structure-only checks (no expression values): row order vs z for every released T2 file (incl. E6.75),
duplicate coordinates, and a recomputation of a few diag_order numbers with independent code.

Usage (DGX, repo root): uv run python experiments/density/review/review_structure.py

Copied into the repo from the independent review of 2026-10-19 (originally run from runs/density_review/);
only script paths changed, plus lint-only edits (line wraps; variable O renamed). Outputs stay on the DGX in
runs/density_review/.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import kendalltau, spearmanr
from sklearn.neighbors import NearestNeighbors

from vec.evaluate import SUBSAMPLE_FRAC
from vec.prep import load_stage

STAGES = ["E6.75", "E7.25", "E8.0", "E8.25_late", "E8.75", "E9.5"]


def canon(C):
    X = C - C.mean(0)
    rms = float(np.sqrt((X ** 2).sum(1).mean()))
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    if np.linalg.det(Vt) < 0:
        Vt[-1] *= -1
    return X @ Vt.T / rms, rms


def main():
    for st in STAGES:
        a = load_stage(f"data/raw/T2/{st}.h5ad", lineage=False)
        C = np.asarray(a.obsm["spatial_3D"], dtype=float)
        n = len(C)
        r = np.arange(n)
        z = C[:, 2]
        # fraction of adjacent row pairs where z increases (0 = perfectly non-increasing)
        inc = float((np.diff(z) > 0).mean())
        n_dup = n - len(np.unique(C[:, :3], axis=0))
        # how many z-runs (maximal blocks of equal z in file order) vs distinct z values
        runs = int(1 + (np.diff(z) != 0).sum())
        Z, rms = canon(C[:, :3])
        d, _ = NearestNeighbors(n_neighbors=2).fit(C[:, :3]).kneighbors(C[:, :3])
        print(f"{st}: n={n} dims={C.shape[1]} z_levels={len(np.unique(z))} z_runs={runs} "
              f"spearman(row,z)={spearmanr(r, z)[0]:.4f} kendall={kendalltau(r, z)[0]:.4f} "
              f"frac_adjacent_z_increase={inc:.5f} dup_coords={n_dup} rms={rms:.1f} nn_full={np.median(d[:, 1]):.2f} "
              f"truth_cells={int(round(n * SUBSAMPLE_FRAC)) // 2}", flush=True)


if __name__ == "__main__":
    main()
