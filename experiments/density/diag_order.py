"""Structure-only diagnostics for the density experiment (no expression values): row order vs coordinates, how much
of the tissue a contiguous row block covers, and how the median nearest-neighbour spacing scales with count.

Usage (DGX): uv run python experiments/density/diag_order.py --out runs/density/diag_order.csv
Stages whose files are missing (the release stages E7.5 and E10.5 before 2026-10-20) are skipped with a message.
n_z_runs counts maximal runs of equal z in file order: n_z_runs == n_z_levels means strictly section by section.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from run_density import resolve  # same folder
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

from vec.evaluate import SUBSAMPLE_FRAC, Stage, split_half, subsample
from vec.prep import load_stage

STAGES = ["E9.5", "E8.75", "E8.25_late", "E8.0", "E7.25", "E6.75", "E10.5", "E7.5"]
TARGETS = {"E9.5", "E8.75", "E8.0"}
RATIOS = [0.25, 0.5, 1, 2, 4, 8]


def canon(C):
    X = C - C.mean(0)
    rms = float(np.sqrt((X ** 2).sum(1).mean()))
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    if np.linalg.det(Vt) < 0:
        Vt[-1] *= -1
    return X @ Vt.T / rms, rms


def voxels(Z, n_grid=16, extent=3.0):
    v = np.floor((np.clip(Z, -extent, extent - 1e-9) + extent) / (2 * extent) * n_grid).astype(int)
    return set(map(tuple, v))


def median_nn(C):
    d, _ = NearestNeighbors(n_neighbors=2).fit(C).kneighbors(C)
    return float(np.median(d[:, 1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/density/diag_order.csv")
    args = ap.parse_args()
    rows = []
    for st in STAGES:
        path, note = resolve(f"data/raw/T2/{st}.h5ad")
        if path is None:
            print(f"{st}: skipped ({note})")
            continue
        a = load_stage(str(path), lineage=False)
        C = np.asarray(a.obsm["spatial_3D"], dtype=float)[:, :3]
        n = len(C)
        Z, rms = canon(C)
        r = np.arange(n)
        base = dict(stage=st, n_cells=n, n_z_levels=len(np.unique(C[:, 2])),
                    n_z_runs=int(1 + (np.diff(C[:, 2]) != 0).sum()), rms_radius=rms,
                    rho_row_z=spearmanr(r, C[:, 2])[0],
                    rho_row_pc1=spearmanr(r, Z[:, 0])[0], rho_row_pc2=spearmanr(r, Z[:, 1])[0],
                    rho_row_pc3=spearmanr(r, Z[:, 2])[0], nn_full=median_nn(C))
        if st not in TARGETS:
            rows.append(base | dict(kind="stage"))
            continue
        # Same bookkeeping as run_density.py, seed 0: working copy, truth half, pool outside the working copy.
        stage = Stage(X=np.arange(n, dtype=np.int64)[:, None], coords=C)
        wc = subsample(stage, SUBSAMPLE_FRAC, 0)
        truth = split_half(wc, seed=0)[0]
        pool = np.setdiff1d(r, wc.X[:, 0])
        tc = truth.n
        vox_all = voxels(Z)
        z_all = np.unique(C[:, 2])
        rows.append(base | dict(kind="stage", truth_cells=tc, nn_truth=median_nn(truth.coords)))
        rng = np.random.default_rng(0)
        for ratio in RATIOS:
            m = min(int(round(ratio * tc)), len(pool))
            rand = pool[np.sort(rng.permutation(len(pool))[:m])]
            centre = len(pool) // 2
            block = pool[centre - m // 2: centre - m // 2 + m]
            head = pool[:m]
            out = dict(stage=st, kind=f"ratio_{ratio}x", truth_cells=tc, n_cells=m)
            for name, idx in (("random", rand), ("block", block), ("head", head)):
                out[f"vox_cover_{name}"] = len(voxels(Z[idx]) & vox_all) / len(vox_all)
                out[f"z_cover_{name}"] = len(np.unique(C[idx, 2])) / len(z_all)
                out[f"rms_ratio_{name}"] = canon(C[idx])[1] / rms
            out["nn_random"] = median_nn(C[rand])
            rows.append(out)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in row.items()})
    for row in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})


if __name__ == "__main__":
    main()
