"""Shared helpers for the coordinate experiments (experiments/coords). Geometry only; no expression values leave
the DGX. All functions take raw obsm["spatial_3D"] arrays in each file's own convention."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]

# Released MERFISH files (training stages only). Nominal developmental times in E-days; E8.25_late is "slightly
# later than E8.25" per the organisers, so its nominal time is uncertain (see REPORT.md).
FILES = {
    "E6.75": "data/raw/T2/E6.75.h5ad",
    "E7.25": "data/raw/T2/E7.25.h5ad",
    "E8.0": "data/raw/T2/E8.0.h5ad",
    "E8.25_late": "data/raw/T2/E8.25_late.h5ad",
    "E8.75": "data/raw/T2/E8.75.h5ad",
    "E9.5": "data/raw/T2/E9.5.h5ad",
    "E9.5_mab21l2_ko": "data/raw/T3/E9.5_mab21l2_ko.h5ad",
}
SETTING = {"E6.75": "embryo", "E7.25": "embryo", "E8.0": "embryo", "E8.25_late": "heart", "E8.75": "heart",
           "E9.5": "heart", "E9.5_mab21l2_ko": "T3"}


def coords_of(a) -> np.ndarray:
    return np.asarray(a.obsm["spatial_3D"], dtype=np.float64)[:, :3]


def rms_radius(C: np.ndarray) -> float:
    """The scorer's RMS radius: sqrt(mean squared distance to the centroid), all three axes, raw units."""
    X = C - C.mean(0)
    return float(np.sqrt((X ** 2).sum(1).mean()))


def rms_xy(C: np.ndarray) -> float:
    X = C[:, :2] - C[:, :2].mean(0)
    return float(np.sqrt((X ** 2).sum(1).mean()))


def principal_sds(C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Principal SDs (descending) and the unit principal axes as columns; population covariance."""
    X = C - C.mean(0)
    w, V = np.linalg.eigh(X.T @ X / X.shape[0])
    order = np.argsort(w)[::-1]
    return np.sqrt(np.clip(w[order], 0, None)), V[:, order]


def snap_z(z: np.ndarray, tol: float = 1e-3) -> tuple[np.ndarray, np.ndarray, float]:
    """Group z values into section levels: values closer than `tol` are one level (floating-point duplicates).

    Returns (level index per cell, level z value (mean of the raw values in the level), largest within-level spread).
    """
    u = np.unique(z)
    gid = np.concatenate([[0], np.cumsum(np.diff(u) > tol)])
    n_lev = int(gid[-1]) + 1
    levels = np.array([u[gid == g].mean() for g in range(n_lev)])
    spread = max(float(np.ptp(u[gid == g])) for g in range(n_lev))
    lev = gid[np.searchsorted(u, z)]
    return lev.astype(int), levels, spread


def inplane_nn(C: np.ndarray, lev: np.ndarray, k: int = 1, min_cells: int = 3) -> np.ndarray:
    """Distance (xy only) from every cell to its k-th nearest neighbour in the same section level (NaN if the level
    has too few cells). Zero distances (duplicate coordinates) are kept."""
    out = np.full(C.shape[0], np.nan)
    for g in np.unique(lev):
        idx = np.flatnonzero(lev == g)
        if idx.size < max(min_cells, k + 1):
            continue
        d, _ = cKDTree(C[idx, :2]).query(C[idx, :2], k=k + 1)
        out[idx] = d[:, k]
    return out


def local_density_2d(C: np.ndarray, lev: np.ndarray, k: int = 10) -> np.ndarray:
    """In-plane local density k / (pi r_k^2) per cell (cells per unit^2), from the k-th in-plane neighbour."""
    r = inplane_nn(C, lev, k=k, min_cells=k + 5)
    with np.errstate(divide="ignore", invalid="ignore"):
        return k / (np.pi * r ** 2)


def vertical_offset(C: np.ndarray, lev: np.ndarray, levels: np.ndarray, dk: int = 1) -> np.ndarray:
    """xy distance from each cell to the nearest cell dk levels above (NaN when that level is empty)."""
    out = np.full(C.shape[0], np.nan)
    trees = {g: cKDTree(C[lev == g, :2]) for g in np.unique(lev) if (lev == g).sum() > 0}
    for g in np.unique(lev):
        if g + dk not in trees:
            continue
        idx = np.flatnonzero(lev == g)
        d, _ = trees[g + dk].query(C[idx, :2], k=1)
        out[idx] = d
    return out


def occupied_area(xy: np.ndarray, pixel: float) -> float:
    """Area covered by a set of points on a square grid of the given pixel size (occupied pixels x pixel^2)."""
    if xy.shape[0] == 0:
        return 0.0
    v = np.floor(xy / pixel).astype(np.int64)
    return float(np.unique(v, axis=0).shape[0] * pixel ** 2)


def summarize(x) -> dict:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"median": np.nan, "p10": np.nan, "p90": np.nan, "mean": np.nan}
    return {"median": float(np.median(x)), "p10": float(np.percentile(x, 10)), "p90": float(np.percentile(x, 90)),
            "mean": float(x.mean())}
