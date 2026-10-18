"""Shared helpers for the coordinate experiments (experiments/coords). Geometry only; no expression values leave
the DGX. All functions take raw obsm["spatial_3D"] arrays in each file's own convention.

Conventions used by every script here (2026-10-19 revision, after the independent review):
- Geometry is measured after dropping duplicated records (`load`, i.e. vec.prep.load_stage(..., drop_duplicates=True));
  scorer emulation (ceilings, mirror Dice, scale skill) uses the files as released (`raw_coords`), because the scorer
  draws its truth halves from them.
- Section levels are paired by their z difference (`level_partners`), never by index, so a stack gap is not counted
  as one step.
- A level that holds two physical sections (several obs["sample_id"] values at one z; only the Mab21l2 KO records
  them) is left out of pair statistics (`multi_section_levels`) and counted as separate sections (`section_ids`).
- Random draws use a fresh generator per stage, seeded by the stage name (`stage_rng`), so a stage's result does not
  depend on which other stages were processed or in what order.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
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
# Validation answers released on 2026-10-20 (training input from then on). The paths follow the naming of the released
# T2 files; the data manifest was not consulted (it needs the team token, which stays on the local machine), so every
# script takes paths on the command line and these are only defaults.
P3_FILES = {
    "E7.5": "data/raw/T2/E7.5.h5ad",
    "E8.5": "data/raw/T2/E8.5.h5ad",
    "E10.5": "data/raw/T2/E10.5.h5ad",
}
SETTING = {"E6.75": "embryo", "E7.25": "embryo", "E7.5": "embryo", "E8.0": "embryo", "E8.25_late": "heart",
           "E8.5": "heart", "E8.75": "heart", "E9.5": "heart", "E10.5": "heart", "E9.5_mab21l2_ko": "T3"}
T_NOMINAL = {"E6.75": 6.75, "E7.25": 7.25, "E7.5": 7.5, "E8.0": 8.0, "E8.25_late": 8.25, "E8.5": 8.5, "E8.75": 8.75,
             "E9.5": 9.5, "E10.5": 10.5, "E9.5_mab21l2_ko": 9.5}


def nominal_t(name: str) -> float:
    """Nominal E-day of a stage name (E8.25_late -> 8.25)."""
    if name in T_NOMINAL:
        return T_NOMINAL[name]
    m = re.match(r"^E(\d+(?:\.\d+)?)", name)
    if not m:
        raise ValueError(f"cannot read an E-day from stage name {name!r}")
    return float(m.group(1))


def stage_specs(items: list[str] | None, default: dict[str, str]) -> dict[str, str]:
    """Stage name -> path from command-line items `NAME=PATH` (the default mapping when none are given)."""
    if not items:
        return dict(default)
    out = {}
    for it in items:
        name, sep, path = it.partition("=")
        if not sep or not name or not path:
            raise SystemExit(f"--stage expects NAME=PATH, got {it!r}")
        out[name] = path
    return out


def stage_rng(name: str, seed: int = 0) -> np.random.Generator:
    """A fresh generator for one stage, seeded by its name and a base seed (stable across runs and processes)."""
    digest = hashlib.sha256(f"{seed}:{name}".encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


def load(path: str | Path, drop_duplicates: bool = True):
    """vec.prep.load_stage with duplicated records dropped (the default for geometry). Returns (adata, summary)."""
    from vec.prep import load_stage

    a = load_stage(path, drop_duplicates=drop_duplicates)
    return a, dedup_summary(a)


def dedup_summary(a) -> dict:
    """What load_stage(..., drop_duplicates=True) removed, from uns["qc"] (zeros when nothing was checked)."""
    qc = a.uns.get("qc") or {}
    return {"n_input": int(qc.get("n_input", a.n_obs)), "n_dropped": int(qc.get("n_dropped", 0)),
            "dropped_by_gap": {str(k): int(v) for k, v in dict(qc.get("by_gap", {})).items()},
            "dropped_max_xy_offset": float(qc.get("max_xy_offset", 0.0)),
            "below_min_library": int(qc.get("below_min_library", 0))}


def dropped_per_level(a, levels: np.ndarray, tol: float = 1e-3) -> np.ndarray:
    """Number of dropped duplicate records per (de-duplicated) level, from uns["qc"]["by_level_pair"]."""
    out = np.zeros(levels.size, dtype=int)
    qc = a.uns.get("qc") or {}
    pairs = qc.get("by_level_pair")
    if pairs is None or len(pairs) == 0:
        return out
    pairs = pd.DataFrame(pairs)
    for z, n in pairs.groupby("dropped_z")["n"].sum().items():
        i = int(np.argmin(np.abs(levels - float(z))))
        if abs(levels[i] - float(z)) <= tol:
            out[i] += int(n)
    return out


def raw_coords(path: str | Path) -> np.ndarray:
    """obsm["spatial_3D"] of a file as released (duplicates included), for scorer emulation."""
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    try:
        return np.asarray(a.obsm["spatial_3D"], dtype=np.float64)[:, :3].copy()
    finally:
        a.file.close()


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


def level_step(levels: np.ndarray) -> float:
    """Typical spacing of sorted levels (median of the positive differences; a few gaps do not move it)."""
    d = np.diff(np.asarray(levels, float))
    d = d[d > 0]
    return float(np.median(d)) if d.size else float("nan")


def level_partners(levels: np.ndarray, k: int, step: float | None = None, tol_frac: float = 0.25) -> np.ndarray:
    """Index of the level whose z lies k steps above each level (|z_h - z_g - k step| <= tol_frac step), else -1.

    Levels are paired by their z difference, not by index: the Mab21l2 KO jumps 30 where its step is 10, so the next
    index there sits three steps away, and a level next to a gap has no partner on that side.
    """
    levels = np.asarray(levels, float)
    if levels.size == 0:
        return np.zeros(0, dtype=int)
    step = level_step(levels) if step is None else float(step)
    target = levels + k * step
    j = np.clip(np.searchsorted(levels, target), 0, levels.size - 1)
    jm = np.clip(j - 1, 0, levels.size - 1)
    best = np.where(np.abs(levels[j] - target) <= np.abs(levels[jm] - target), j, jm)
    ok = np.abs(levels[best] - target) <= tol_frac * step
    return np.where(ok, best, -1).astype(int)


def multi_section_levels(a, lev: np.ndarray) -> list[int]:
    """Levels holding more than one physical section (several obs["sample_id"] values); [] without sample ids."""
    if "sample_id" not in a.obs:
        return []
    n = pd.Series(np.asarray(a.obs["sample_id"]).astype(str)).groupby(np.asarray(lev)).nunique()
    return sorted(int(g) for g, v in n.items() if v > 1)


def section_ids(a, lev: np.ndarray) -> np.ndarray:
    """Physical section per cell: (level, sample_id) when the file records sample ids, else the level."""
    lev = np.asarray(lev)
    if "sample_id" not in a.obs:
        return lev.astype(str)
    return np.char.add(np.char.add(lev.astype(str), ":"), np.asarray(a.obs["sample_id"]).astype(str))


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


def vertical_offset(C: np.ndarray, lev: np.ndarray, levels: np.ndarray, dk: int = 1,
                    step: float | None = None) -> np.ndarray:
    """xy distance from each cell to the nearest cell in the level dk steps above, paired by z difference (NaN when
    no level sits there or it is empty)."""
    out = np.full(C.shape[0], np.nan)
    partner = level_partners(levels, dk, step)
    present = set(np.unique(lev).tolist())
    trees = {g: cKDTree(C[lev == g, :2]) for g in present}
    for g in present:
        h = int(partner[g])
        if h < 0 or h not in trees:
            continue
        idx = np.flatnonzero(lev == g)
        d, _ = trees[h].query(C[idx, :2], k=1)
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
