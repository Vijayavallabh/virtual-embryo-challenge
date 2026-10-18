"""Characterise records that `find_duplicate_records` does not flag. Aggregates only (run on the DGX):
    uv run python experiments/qc/probe_records.py data/raw/T2/*.h5ad data/raw/T3/*.h5ad \
        --out experiments/qc/results/probe_records_2026-10-19.json

A. Identical count rows (library >= 20) anywhere in a file: level gap and xy offset between each repeat and the first
   occurrence, per level pair, with an affine fit of repeat xy on first-occurrence xy (a re-registered copy fits an
   affine map exactly; unrelated cells do not).
B. Rows with non-integer counts: fractional parts, levels, and whether 2 * row = A + B for a record A near it in xy
   (within --radius, levels within 2) and any record B (hash lookup).
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from vec.spatial import row_hashes, snap_z


def canonical(M) -> sparse.csr_matrix:
    C = sparse.csr_matrix(M, dtype=np.float64, copy=True)
    C.sum_duplicates()
    C.eliminate_zeros()
    return C


def same_row(C, i, j) -> bool:
    a, b = slice(C.indptr[i], C.indptr[i + 1]), slice(C.indptr[j], C.indptr[j + 1])
    return bool(np.array_equal(C.indices[a], C.indices[b]) and np.array_equal(C.data[a], C.data[b]))


def affine_fit(src: np.ndarray, dst: np.ndarray) -> dict:
    A = np.c_[src, np.ones(len(src))]
    coef, *_ = np.linalg.lstsq(A, dst, rcond=None)
    resid = dst - A @ coef
    lin = coef[:2].T
    return {"rms_residual": float(np.sqrt((resid ** 2).sum(1).mean())), "det": float(np.linalg.det(lin)),
            "rotation_deg": float(np.degrees(np.arctan2(lin[1, 0], lin[0, 0]))),
            "shift": [float(v) for v in coef[2]]}


def part_a(C, xy, zc, step, lib, min_library) -> dict:
    keys = row_hashes(C)
    order = np.lexsort((np.arange(keys.size), keys))
    k = keys[order]
    starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    ends = np.r_[starts[1:], k.size]
    rec = []
    for s, e in zip(starts, ends, strict=True):
        if e - s < 2:
            continue
        m = order[s:e]
        first = m[0]
        if lib[first] < min_library:
            continue
        for j in m[1:]:
            if same_row(C, first, j):
                d = xy[j] - xy[first]
                rec.append((round(float(zc[first]), 6), round(float(zc[j]), 6),
                            int(np.rint((zc[j] - zc[first]) / step)), d[0], d[1], first, j))
    if not rec:
        return {"repeats": 0}
    df = pd.DataFrame(rec, columns=["z_first", "z_repeat", "gap", "dx", "dy", "i", "j"])
    df["cheb"] = np.maximum(df.dx.abs(), df.dy.abs())
    df["colocated"] = df.cheb <= 0.01
    out = {"repeats": int(len(df)), "colocated": int(df.colocated.sum()),
           "by_gap": {str(g): int(c) for g, c in df.gap.value_counts().sort_index().items()},
           "by_gap_not_colocated": {str(g): int(c)
                                    for g, c in df[~df.colocated].gap.value_counts().sort_index().items()},
           "xy_offset_quantiles_not_colocated": ([float(v) for v in df[~df.colocated].cheb.quantile([0, .1, .5, .9, 1])]
                                                 if (~df.colocated).any() else [])}
    pairs = []
    for (z1, z2), g in df.groupby(["z_first", "z_repeat"]):
        row = {"z_first": z1, "z_repeat": z2, "gap": int(g.gap.iloc[0]), "n": int(len(g)),
               "colocated": int(g.colocated.sum()), "median_dx": float(g.dx.median()),
               "median_dy": float(g.dy.median()), "sd_dx": float(g.dx.std()) if len(g) > 1 else 0.0,
               "sd_dy": float(g.dy.std()) if len(g) > 1 else 0.0}
        if len(g) >= 10 and not g.colocated.all():
            row["affine"] = affine_fit(xy[g.i.to_numpy()], xy[g.j.to_numpy()])
        pairs.append(row)
    out["by_level_pair"] = sorted(pairs, key=lambda r: -r["n"])[:25]
    out["level_pairs"] = len(pairs)
    return out


def part_b(C, xy, zc, zl, radius) -> dict:
    frac = np.abs(C.data - np.rint(C.data))
    bad = frac > 1e-6
    row_of = np.repeat(np.arange(C.shape[0]), np.diff(C.indptr))
    rows = np.unique(row_of[bad])
    out = {"rows": int(rows.size)}
    if not rows.size:
        return out
    out["fraction_values"] = {f"{v:.4f}": int(c) for v, c in
                              pd.Series(np.round(frac[bad], 4)).value_counts().head(8).items()}
    per_level = Counter(np.round(zc[rows], 6).tolist())
    out["levels"] = {str(k): int(v) for k, v in sorted(per_level.items())}
    other = bad & (np.abs(frac - 0.5) > 1e-6)
    out["half_only_rows"] = int(rows.size - np.unique(row_of[other]).size)
    keys = row_hashes(C)
    lookup: dict[int, list[int]] = {}
    for idx, h in enumerate(keys.tolist()):
        lookup.setdefault(h, []).append(idx)
    tree = cKDTree(xy)
    lev = zl.index
    found = Counter()
    n_found = 0
    for r in rows:
        cand = [a for a in tree.query_ball_point(xy[r], r=radius, p=np.inf) if a != r and abs(lev[a] - lev[r]) <= 2]
        if not cand:
            continue
        Dd = 2 * C[r].toarray() - C[cand].toarray()  # 2 * row - A for each nearby candidate A (dense, 500 genes)
        ok = (Dd >= 0).all(1) & (np.abs(Dd - np.rint(Dd)) <= 1e-6).all(1)
        if not ok.any():
            continue
        hit = False
        for t, h in zip(np.flatnonzero(ok), row_hashes(canonical(Dd[ok])).tolist(), strict=True):
            for b in lookup.get(h, []):
                if b != r and np.array_equal(C[b].toarray().ravel(), Dd[t]):
                    a = cand[t]
                    found[(int(lev[a] - lev[r]), int(lev[b] - lev[r]),
                           bool(np.abs(xy[a] - xy[r]).max() <= 0.01), bool(np.abs(xy[b] - xy[r]).max() <= 0.01))] += 1
                    hit = True
        n_found += hit
    out["rows_with_2r_eq_A_plus_B"] = int(n_found)
    out["decompositions"] = [{"gap_A": k[0], "gap_B": k[1], "A_colocated": k[2], "B_colocated": k[3], "n": int(v)}
                             for k, v in found.most_common(12)]
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("files", nargs="+", type=Path)
    p.add_argument("--min-library", type=float, default=20)
    p.add_argument("--radius", type=float, default=30.0)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    res = []
    for f in args.files:
        a = ad.read_h5ad(f)
        C = canonical(a.layers["counts"])
        xyz = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)[:, :3]
        zl = snap_z(xyz[:, 2])
        zc = zl.levels[zl.index]
        lib = np.asarray(C.sum(1)).ravel()
        r = {"file": str(f), "identical": part_a(C, xyz[:, :2], zc, zl.step, lib, args.min_library),
             "non_integer": part_b(C, xyz[:, :2], zc, zl, args.radius)}
        print(json.dumps(r))
        res.append(r)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
