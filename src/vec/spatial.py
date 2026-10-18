"""Spatial data QC for the MERFISH stages: z-level snapping, duplicated records and a per-file intake summary.

`obsm["spatial_3D"][:, 2]` is a section index times a per-stage z unit. The raw values carry float noise (E9.5 stores 98
distinct values on 40 levels, differing by up to 1e-13), a stack can skip levels (the Mab21l2 KO jumps 30 where its step
is 10) and one level can hold two sections. Some records are stored more than once: identical count rows at the same xy
in an adjacent level (copied sections) or in the same level (a section stored twice with ~0.001 xy offsets). The
released heart and embryo files also hold copies whose xy was moved (identical rows one level away, 0.06-63 units off;
`xy_tol=inf` catches them) and interpolated records (half-integer counts equal to the mean of two records in the levels
above and below; `qc_stage` counts them as non-integer rows). See experiments/qc/ for the numbers per file.

- `snap_z`: raw z to level index per cell, level values and step.
- `find_duplicate_records`: mask of redundant copies (the first occurrence in file order is kept) plus a small report.
- `qc_stage`: aggregate summary of one file for intake checks (`vec data qc FILE...`).

The scorer draws its truth halves from the organisers' files as they are, duplicates included, so the board emulator
loads targets unchanged; only prediction pipelines drop duplicates (`vec.prep.load_stage(..., drop_duplicates=True)`).
"""

from __future__ import annotations

from typing import NamedTuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from vec.io import COORD_KEY

NOISE_REL = 1e-9  # z values closer than this fraction of max|z| are one level whatever the step
STEP_FRAC = 0.25  # default snapping tolerance, as a fraction of the step
SECTION_KEYS = ("sample_id", "section")  # obs columns naming the physical section, checked against the z levels
MAX_LISTED_LEVELS = 500  # per-level lists are left out of qc_stage beyond this (z is then not sectioned)


class ZLevels(NamedTuple):
    """Snapped z: `index` per cell (0..n_levels-1 in ascending z), ascending level values, and the grid step."""

    index: np.ndarray
    levels: np.ndarray
    step: float

    @property
    def grid(self) -> np.ndarray:
        """Level positions in steps above the lowest level: whole numbers on a regular stack, jumps over gaps."""
        if self.levels.size < 2 or not np.isfinite(self.step):
            return np.zeros(self.levels.size)
        return (self.levels - self.levels[0]) / self.step

    def gaps(self, min_steps: float = 1.5) -> list[dict]:
        """Consecutive levels more than `min_steps` steps apart, with the number of whole levels missing between."""
        d = np.diff(self.grid)
        return [{"from": round(float(self.levels[i]), 6), "to": round(float(self.levels[i + 1]), 6),
                 "steps": round(float(d[i]), 3),
                 "missing": int(round(float(d[i]))) - 1} for i in np.flatnonzero(d > min_steps)]


def _chain(u: np.ndarray, tol: float) -> np.ndarray:
    """Group label per sorted value: a new group starts wherever the gap to the previous value exceeds `tol`."""
    return np.concatenate([[0], np.cumsum(np.diff(u) > tol)]).astype(np.int64)


def _level_values(u: np.ndarray, cnt: np.ndarray, group: np.ndarray) -> np.ndarray:
    """Count-weighted mean of each group of sorted values (exactly the value itself for a single-valued group)."""
    base = u[np.flatnonzero(np.r_[True, np.diff(group) > 0])]
    off = np.bincount(group, weights=(u - base[group]) * cnt) / np.bincount(group, weights=cnt)
    return base + off


def _grid_step(levels: np.ndarray) -> float:
    """Typical spacing of sorted levels, robust to missing levels (whole multiples) and a minority of sub-step gaps."""
    d = np.diff(levels)
    d = d[d > 0]
    if d.size == 0:
        return float("nan")
    s = float(np.median(d[d >= 0.4 * np.quantile(d, 0.75)]))
    k = np.rint(d / s)
    on = (k >= 1) & (np.abs(d / s - k) <= 0.2)
    return float(np.median(d[on] / k[on])) if on.any() else s


def snap_z(z, tol: float | None = None) -> ZLevels:
    """Snap raw z values to section levels.

    Sorted distinct values are merged wherever consecutive ones are at most `tol` apart, so float noise and a section
    stored twice at a small z offset collapse onto one level. Nothing is rounded against a fixed origin, so levels on
    x.5 values (or at half-step grid positions) stay distinct. With `tol=None` the tolerance is a quarter of the step,
    which is first estimated after merging float noise only. `step` is the typical spacing of consecutive levels,
    robust to missing levels and to a minority of sub-step spacings; it is nan with fewer than two levels. `grid` and
    `gaps()` on the result give level positions in steps and the missing stretches.
    """
    z = np.asarray(z, dtype=np.float64).ravel()
    if not np.isfinite(z).all():
        raise ValueError("z contains NaN or infinite values")
    if z.size == 0:
        return ZLevels(np.zeros(0, dtype=np.int64), np.zeros(0), float("nan"))
    u, inv, cnt = np.unique(z, return_inverse=True, return_counts=True)
    if tol is None:
        noise = NOISE_REL * max(1.0, float(np.abs(u).max()))
        step0 = _grid_step(_level_values(u, cnt, _chain(u, noise)))
        tol = STEP_FRAC * step0 if np.isfinite(step0) else noise
    elif not tol >= 0:
        raise ValueError(f"tol must be >= 0, got {tol}")
    group = _chain(u, tol)
    levels = _level_values(u, cnt, group)
    return ZLevels(index=group[inv.ravel()], levels=levels, step=_grid_step(levels))


# --- row identity -------------------------------------------------------------------------------------------------

_M1, _M2, _GOLD = np.uint64(0xBF58476D1CE4E5B9), np.uint64(0x94D049BB133111EB), np.uint64(0x9E3779B97F4A7C15)


def _mix64(x: np.ndarray) -> np.ndarray:
    """splitmix64 finaliser on a uint64 array (wraps modulo 2**64)."""
    x = x ^ (x >> np.uint64(30))
    x = x * _M1
    x = x ^ (x >> np.uint64(27))
    x = x * _M2
    return x ^ (x >> np.uint64(31))


def _canonical_csr(M) -> sparse.csr_matrix:
    """float64 CSR with summed duplicates, sorted indices and no explicit zeros: equal rows have equal storage."""
    C = sparse.csr_matrix(M, dtype=np.float64, copy=True)
    C.sum_duplicates()
    C.eliminate_zeros()
    return C


def row_hashes(C: sparse.csr_matrix) -> np.ndarray:
    """64-bit content hash per row of a canonical CSR matrix: equal rows always hash equal (collisions ~n^2/2^64)."""
    nnz = np.diff(C.indptr)
    h = _mix64(C.indices.astype(np.uint64) * _GOLD + _mix64(np.ascontiguousarray(C.data).view(np.uint64)))
    keys = _mix64(nnz.astype(np.uint64) + _GOLD)
    rows = np.flatnonzero(nnz)
    if rows.size:  # reduceat over the starts of non-empty rows sums exactly each row's entries
        keys[rows] += np.add.reduceat(h, C.indptr[rows])
    return keys


def _groups(keys: np.ndarray):
    """Row indices (in file order) of each hash value shared by two or more rows."""
    order = np.lexsort((np.arange(keys.size), keys))
    k = keys[order]
    starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    ends = np.r_[starts[1:], k.size]
    for s, e in zip(starts, ends, strict=True):
        if e - s > 1:
            yield order[s:e]


def _rows_equal(C: sparse.csr_matrix, i: int, j: int) -> bool:
    a, b = slice(C.indptr[i], C.indptr[i + 1]), slice(C.indptr[j], C.indptr[j + 1])
    return bool(np.array_equal(C.indices[a], C.indices[b]) and np.array_equal(C.data[a], C.data[b]))


def _identical_rows(C: sparse.csr_matrix, min_library: float | None) -> dict:
    keys = row_hashes(C)
    _, first, counts = np.unique(keys, return_index=True, return_counts=True)
    out = {"repeats": int(C.shape[0] - first.size), "groups": int((counts > 1).sum()),
           "largest_group": int(counts.max()) if counts.size else 0,
           "empty_rows": int((np.diff(C.indptr) == 0).sum())}
    if min_library is not None:
        repeat = np.ones(C.shape[0], dtype=bool)
        repeat[first] = False
        out["repeats_min_library"] = int((repeat & (np.asarray(C.sum(1)).ravel() >= min_library)).sum())
    return out


def identical_rows(M, min_library: float | None = None) -> dict:
    """Rows repeating the content of an earlier row anywhere in the matrix (position ignored), by 64-bit row hash.

    With `min_library`, also counts the repeats whose row sum is at least that (low-count rows coincide by chance).
    """
    return _identical_rows(_canonical_csr(M), min_library)


# --- duplicated records ---------------------------------------------------------------------------------------------

def _close_pairs(xy: np.ndarray, tol: float) -> np.ndarray:
    """Index pairs (i < j) whose xy agree within `tol` on both axes."""
    if xy.shape[0] <= 64:
        d = np.abs(xy[:, None, :] - xy[None, :, :]).max(-1)
        return np.argwhere(np.triu(d <= tol, k=1))
    return cKDTree(xy).query_pairs(tol, p=np.inf, output_type="ndarray")


def _candidate_pairs(xy: np.ndarray, zc: np.ndarray, xy_tol: float, step: float, max_level_gap: int) -> np.ndarray:
    """Local index pairs within one group of identical rows that meet the xy and level criteria.

    With an infinite `xy_tol`, members sorted by level are paired with their successor when within `max_level_gap`
    levels: linear in the group size, and the same components as all qualifying pairs (a 1-D threshold graph).
    """
    if np.isinf(xy_tol):
        order = np.argsort(zc, kind="stable")
        ok = _level_steps(np.diff(zc[order]), step) <= max_level_gap
        return np.c_[order[:-1][ok], order[1:][ok]]
    pairs = _close_pairs(xy, xy_tol)
    if pairs.size:
        pairs = pairs[_level_steps(zc[pairs[:, 0]] - zc[pairs[:, 1]], step) <= max_level_gap]
    return pairs.reshape(-1, 2)


def _find(parent: np.ndarray, i: int) -> int:
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


def _roots(parent: np.ndarray) -> np.ndarray:
    """Root of every node by pointer jumping (parent[i] <= i throughout, so roots are first occurrences)."""
    root = parent.copy()
    while True:
        nxt = root[root]
        if np.array_equal(nxt, root):
            return root
        root = nxt


def _level_steps(dz: np.ndarray, step: float) -> np.ndarray:
    """|dz| in whole steps (0 for the same level); without a step only the same level counts as 0."""
    dz = np.abs(dz)
    if not np.isfinite(step):
        return np.where(dz == 0, 0, np.iinfo(np.int64).max)
    return np.rint(dz / step).astype(np.int64)


def _duplicates(C: sparse.csr_matrix, coords: np.ndarray, xy_tol: float, min_library: float,
                max_level_gap: int) -> tuple[np.ndarray, dict]:
    n = C.shape[0]
    xy = coords[:, :2]
    zl = snap_z(coords[:, 2])
    zc = zl.levels[zl.index]
    lib = np.asarray(C.sum(1)).ravel()
    parent = np.arange(n)
    for members in _groups(row_hashes(C)):
        for a, b in members[_candidate_pairs(xy[members], zc[members], xy_tol, zl.step, max_level_gap)]:
            if _rows_equal(C, a, b):
                ra, rb = _find(parent, a), _find(parent, b)
                parent[max(ra, rb)] = min(ra, rb)  # the root is the first occurrence in file order
    root = _roots(parent)
    linked = root != np.arange(n)
    mask = linked & (lib[root] >= min_library)  # copies share their row, hence their library
    dropped = np.flatnonzero(mask)
    kept_z, drop_z = zc[root[dropped]], zc[dropped]
    gap = _level_steps(drop_z - kept_z, zl.step)
    by_pair = (pd.DataFrame({"kept_z": np.round(kept_z, 6), "dropped_z": np.round(drop_z, 6), "gap": gap})
               .value_counts().rename("n").reset_index().sort_values(["kept_z", "dropped_z"], ascending=False))
    offset = np.abs(xy[dropped] - xy[root[dropped]]).max(1) if dropped.size else np.zeros(0)
    report = {
        "n_records": int(n),
        "n_redundant": int(dropped.size),
        "n_with_copies": int(np.unique(root[dropped]).size),
        "by_gap": {str(g): int(c) for g, c in zip(*np.unique(gap, return_counts=True), strict=True)},
        "by_level_pair": [{"kept_z": float(r.kept_z), "dropped_z": float(r.dropped_z), "gap": int(r.gap),
                           "n": int(r.n)} for r in by_pair.itertuples()],
        "max_xy_offset": float(offset.max()) if offset.size else 0.0,
        "below_min_library": int((linked & ~mask).sum()),
        "xy_tol": float(xy_tol) if np.isfinite(xy_tol) else "any", "min_library": float(min_library),
        "max_level_gap": int(max_level_gap),
        "z_step": float(zl.step),
    }
    return mask, report


def find_duplicate_records(adata: ad.AnnData, xy_tol: float = 0.01, min_library: float = 20, max_level_gap: int = 1,
                           layer: str = "counts", coord_key: str = COORD_KEY) -> tuple[np.ndarray, dict]:
    """Flag records stored more than once.

    A record is a copy when another record has an identical row in `layers[layer]` (compared exactly), xy within
    `xy_tol` on both axes (`xy_tol=inf` ignores xy) and a z level at most `max_level_gap` levels away (same level
    included; levels from `snap_z`, distance in steps, so a level across a gap is not adjacent). Copies chain (a record
    copied into two neighbouring levels is one record) and every copy after the first occurrence in file order is
    flagged. Rows with fewer than `min_library` total counts are never flagged, since low-count rows can coincide by
    chance; such matches are counted under `below_min_library`. Rows are hashed first, so only identical rows are
    compared.

    Returns (mask, report): mask is True for redundant records; the report counts them by level gap (`by_gap`) and by
    (kept level, dropped level) pair (`by_level_pair`), with the largest xy offset between a copy and its kept record.
    """
    if layer not in adata.layers:
        raise KeyError(f"layers[{layer!r}] is required (identical count rows define a copy)")
    if coord_key not in adata.obsm:
        raise KeyError(f"obsm[{coord_key!r}] is required (copies sit at the same xy in the same or an adjacent level)")
    coords = np.asarray(adata.obsm[coord_key], dtype=np.float64)[:, :3]
    return _duplicates(_canonical_csr(adata.layers[layer]), coords, xy_tol, min_library, max_level_gap)


# --- per-file QC ------------------------------------------------------------------------------------------------------

def _zkey(v: float) -> str:
    return str(round(float(v), 6))


def _dense_levels(per_level: np.ndarray, ratio: float = 1.6) -> dict[int, float]:
    """Levels holding more than `ratio` times the mean cell count of their neighbouring levels, with that ratio.

    A hint for a level holding two sections (about 2x); natural section-to-section variation can reach ~1.8x.
    """
    c = per_level.astype(float)
    if c.size < 2:
        return {}
    left, right = np.r_[np.nan, c[:-1]], np.r_[c[1:], np.nan]
    r = c / np.nanmean(np.vstack([left, right]), axis=0)
    return {int(i): round(float(r[i]), 2) for i in np.flatnonzero(r > ratio)}


def _integrality(C: sparse.csr_matrix) -> dict:
    frac = np.abs(C.data - np.rint(C.data))
    bad = frac > 1e-6
    rows = np.repeat(np.arange(C.shape[0]), np.diff(C.indptr))[bad]
    return {"non_integer_entries": int(bad.sum()), "non_integer_rows": int(np.unique(rows).size),
            "half_integer_entries": int((np.abs(frac - 0.5) <= 1e-6).sum()),
            "negative_entries": int((C.data < 0).sum())}


def _agreement(labels: np.ndarray, level: np.ndarray, levels: np.ndarray) -> dict:
    """Does each section label sit on one z level, and each level hold one label?"""
    df = pd.DataFrame({"s": labels, "l": level})
    per_label = df.groupby("s")["l"].unique()
    per_level = df.groupby("l")["s"].unique()
    return {
        "n": int(per_label.size),
        "one_level_per_label": bool((per_label.map(len) == 1).all()),
        "one_label_per_level": bool((per_level.map(len) == 1).all()),
        "labels_on_several_levels": {str(k): sorted(round(float(levels[i]), 6) for i in v)
                                     for k, v in per_label.items() if len(v) > 1},
        "levels_with_several_labels": {_zkey(levels[k]): sorted(map(str, v))
                                       for k, v in per_level.items() if len(v) > 1},
    }


def _axis_summary(C: np.ndarray) -> dict:
    lo, hi, mean = C.min(0), C.max(0), C.mean(0)
    extent, mid = hi - lo, (hi + lo) / 2
    tol = 1e-3 * np.maximum(extent, 1.0)
    centred_on = "mean" if (np.abs(mean) <= tol).all() else "bbox" if (np.abs(mid) <= tol).all() else "no"
    return {"min": lo.tolist(), "max": hi.tolist(), "extent": extent.tolist(), "mean": mean.tolist(),
            "bbox_mid": mid.tolist(), "centred_on": centred_on}


def qc_stage(adata: ad.AnnData, xy_tol: float = 0.01, min_library: float = 20, max_level_gap: int = 1,
             layer: str = "counts", coord_key: str = COORD_KEY) -> dict:
    """Intake QC of one file as read by `anndata.read_h5ad` (not `load_stage` output); every value is an aggregate.

    Reports cells, genes, duplicated obs names, count-layer integrality, rows repeating an earlier row anywhere, and,
    when `obsm[coord_key]` exists: raw z values, levels, step, spacings, gaps, off-grid and unusually dense levels
    (a level holding two sections), section labels against levels (`sample_id`, `section`), row order (Spearman of row
    index vs z level and the number of z runs), duplicated records (`find_duplicate_records` at `xy_tol`, and again
    with xy ignored as `duplicates_any_xy`) and coordinate ranges and centring. Files without coordinates (T1 scRNA)
    get the row-identity part only.
    """
    has_counts = layer in adata.layers
    C = _canonical_csr(adata.layers[layer] if has_counts else adata.X)
    out: dict = {"cells": int(adata.n_obs), "genes": int(adata.n_vars),
                 "obs_names_duplicated": int(adata.obs_names.duplicated().sum()),
                 "matrix": f"layers[{layer!r}]" if has_counts else "X"}
    if has_counts:
        out["counts"] = _integrality(C)
    out["identical_rows"] = _identical_rows(C, min_library if has_counts else None)
    if coord_key not in adata.obsm:
        return _json_safe(out)

    coords = np.asarray(adata.obsm[coord_key], dtype=np.float64)[:, :3]
    z = coords[:, 2]
    zl = snap_z(z)
    n_lev = zl.levels.size
    per_level = np.bincount(zl.index, minlength=n_lev)
    spread = np.zeros(n_lev)
    np.maximum.at(spread, zl.index, z - zl.levels[zl.index])
    low = np.zeros(n_lev)
    np.minimum.at(low, zl.index, z - zl.levels[zl.index])
    grid = zl.grid
    off_grid = np.abs(grid - np.rint(grid)) > 0.1
    spacing = pd.Series(np.round(np.diff(zl.levels), 6)).value_counts().sort_index()
    out["z"] = {
        "raw_values": int(np.unique(z).size),
        "levels": int(n_lev),
        "step": zl.step,
        "spacings": {f"{k:g}": int(v) for k, v in spacing.items()},
        "max_level_spread": float((spread - low).max()) if n_lev else 0.0,
        "gaps": zl.gaps(),
        "off_grid_levels": [round(float(v), 6) for v in zl.levels[off_grid]],
        "dense_levels": {_zkey(zl.levels[i]): {"cells": int(per_level[i]), "neighbour_ratio": r}
                         for i, r in _dense_levels(per_level).items()},
        "continuous": bool(n_lev > MAX_LISTED_LEVELS and n_lev > 0.2 * adata.n_obs),
    }
    if n_lev <= MAX_LISTED_LEVELS:
        out["z"]["level_values"] = np.round(zl.levels, 6).tolist()
        out["z"]["cells_per_level"] = per_level.tolist()
    out["sections"] = {"levels": int(n_lev)}
    for key in SECTION_KEYS:
        if key in adata.obs:
            out["sections"][key] = _agreement(adata.obs[key].astype(str).to_numpy(), zl.index, zl.levels)
    rho = float(spearmanr(np.arange(adata.n_obs), zl.index).statistic) if n_lev > 1 else float("nan")
    runs = int(1 + np.count_nonzero(np.diff(zl.index))) if adata.n_obs else 0
    out["row_order"] = {"spearman_row_vs_z": rho, "z_runs": runs, "grouped_by_level": runs == n_lev,
                        "direction": "descending" if rho < 0 else "ascending" if rho > 0 else "none"}
    if has_counts:
        out["duplicates"] = _duplicates(C, coords, xy_tol, min_library, max_level_gap)[1]
        out["duplicates_any_xy"] = _duplicates(C, coords, np.inf, min_library, max_level_gap)[1]
    out["coords"] = _axis_summary(coords)
    out["coords"]["duplicate_positions"] = int(coords.shape[0] - np.unique(coords, axis=0).shape[0])
    return _json_safe(out)


def _json_safe(o):
    """Plain Python types for json.dumps: numpy scalars and arrays converted, nan/inf to None."""
    if isinstance(o, dict):
        return {str(k): _json_safe(v) for k, v in o.items()}
    if isinstance(o, list | tuple | np.ndarray):
        return [_json_safe(v) for v in (o.tolist() if isinstance(o, np.ndarray) else o)]
    if isinstance(o, bool | np.bool_):
        return bool(o)
    if isinstance(o, int | np.integer):
        return int(o)
    if isinstance(o, float | np.floating):
        return float(o) if np.isfinite(o) else None
    return o


def summary_row(file: str, qc: dict) -> dict:
    """One flat row per file from a `qc_stage` dict, for the CLI table and CSV."""
    z, dup, order = qc.get("z", {}), qc.get("duplicates", {}), qc.get("row_order", {})
    anyxy = qc.get("duplicates_any_xy", {})
    sec = qc.get("sections", {})
    labels = [f"{k}:{'ok' if v['one_level_per_label'] and v['one_label_per_level'] else 'MISMATCH'}"
              for k, v in sec.items() if isinstance(v, dict)]
    return {
        "file": file,
        "cells": qc["cells"],
        "genes": qc["genes"],
        "obs_names_duplicated": qc["obs_names_duplicated"],
        "non_integer_count_rows": qc.get("counts", {}).get("non_integer_rows", ""),
        "identical_rows": qc["identical_rows"]["repeats"],
        "identical_rows_min_library": qc["identical_rows"].get("repeats_min_library", ""),
        "z_raw_values": z.get("raw_values", ""),
        "z_levels": z.get("levels", ""),
        "z_step": "" if z.get("step") is None else round(z["step"], 6),
        "z_spacings": ";".join(f"{k}x{v}" for k, v in z.get("spacings", {}).items()),
        "z_gaps": ";".join(f"{g['from']:g}..{g['to']:g} ({g['missing']} missing)" for g in z.get("gaps", [])),
        "z_max_level_spread": z.get("max_level_spread", ""),
        "z_off_grid_levels": len(z.get("off_grid_levels", [])) if z else "",
        "z_dense_levels": ";".join(f"{k}:{v['cells']} ({v['neighbour_ratio']}x)"
                                   for k, v in z.get("dense_levels", {}).items()),
        "section_labels": ";".join(labels),
        "spearman_row_vs_z": "" if order.get("spearman_row_vs_z") is None else round(order["spearman_row_vs_z"], 4),
        "z_runs": order.get("z_runs", ""),
        "dup_redundant": dup.get("n_redundant", ""),
        "dup_by_gap": ";".join(f"gap{k}:{v}" for k, v in dup.get("by_gap", {}).items()),
        "dup_max_xy_offset": dup.get("max_xy_offset", ""),
        "dup_below_min_library": dup.get("below_min_library", ""),
        "dup_any_xy_redundant": anyxy.get("n_redundant", ""),
        "dup_any_xy_by_gap": ";".join(f"gap{k}:{v}" for k, v in anyxy.get("by_gap", {}).items()),
        "dup_any_xy_max_xy_offset": anyxy.get("max_xy_offset", ""),
        "duplicate_positions": qc.get("coords", {}).get("duplicate_positions", ""),
        "centred_on": qc.get("coords", {}).get("centred_on", ""),
    }
