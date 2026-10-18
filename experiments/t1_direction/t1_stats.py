"""Per-group sufficient statistics of a released T1 stage (E8.5, E9.5; E10.5 once released).

For every (cell type, library, random half) group: the number of cells, and per gene the sum of log1p-CP10k (.X as
released), the sum of CP10k (expm1 of .X) and the number of cells with a non-zero value; plus the summed per-cell depth
(UMIs recovered from the smallest non-zero value, N = 1e4 / expm1(min)). Any grouping (harmonised types, halves, the
whole population) is then a sum of rows. Halves are drawn per cell type with numpy's default_rng(seed) over the sorted
type names, as in the 2026-10-18 research run (same seed, same assignment).

Usage (DGX, repo root): uv run python experiments/t1_direction/t1_stats.py <h5ad> <stage> <out.npz> [--label-key K]
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

MAX_LIBRARIES = 8  # more distinct obs-name suffixes than this means the suffix is not a library tag


def sha256_file(path: str | Path, chunk: int = 1 << 23) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def library_tags(obs_names) -> tuple[np.ndarray, str]:
    lib = pd.Series(np.asarray(obs_names, dtype=str)).str.split("_").str[-1].to_numpy(dtype=object)
    if len(set(lib)) > MAX_LIBRARIES:
        return np.full(len(lib), "0", dtype=object), "obs-name suffix is not a library tag; one library assumed"
    return lib, "library = obs-name suffix after the last '_'"


def compute(path: str | Path, stage: str, label_key: str = "celltype", seed: int = 0) -> dict:
    a = ad.read_h5ad(path)
    if label_key not in a.obs:
        raise KeyError(f"{path}: obs has no {label_key!r} column (columns: {list(a.obs.columns)}); pass --label-key")
    X = a.X.tocsr().astype(np.float32) if sparse.issparse(a.X) else sparse.csr_matrix(a.X, dtype=np.float32)
    X.eliminate_zeros()
    n, _ = X.shape
    if (np.diff(X.indptr) == 0).any():
        raise ValueError(f"{path}: {int((np.diff(X.indptr) == 0).sum())} cells have no expressed gene")
    ct = a.obs[label_key].astype(str).values
    lib, lib_note = library_tags(a.obs_names)
    rng = np.random.default_rng(seed)
    half = np.empty(n, dtype=object)
    for c in np.unique(ct):  # stratified random halves within each type
        idx = np.where(ct == c)[0]
        perm = rng.permutation(idx)
        half[perm[: len(idx) // 2]] = "A"
        half[perm[len(idx) // 2:]] = "B"
    keys = pd.Series([f"{c}|{li}|{h}" for c, li, h in zip(ct, lib, half, strict=True)])
    uk, inv = np.unique(keys.values.astype(str), return_inverse=True)
    K = len(uk)
    Ind = sparse.csr_matrix((np.ones(n, np.float32), (inv, np.arange(n))), shape=(K, n))
    S_log = np.asarray((Ind @ X).todense(), dtype=np.float64)
    Xc = X.copy()
    Xc.data = np.expm1(Xc.data)
    S_cp = np.asarray((Ind @ Xc).todense(), dtype=np.float64)
    Xb = X.copy()
    Xb.data = np.ones_like(Xb.data)
    S_nz = np.asarray((Ind @ Xb).todense(), dtype=np.float64)
    ncell = np.asarray(Ind.sum(1)).ravel()
    mins = np.minimum.reduceat(X.data, X.indptr[:-1])
    depth = 1e4 / np.expm1(mins)
    return {"genes": np.asarray(a.var_names, dtype=str), "keys": uk.astype(str), "n": ncell, "S_log": S_log,
            "S_cp": S_cp, "S_nz": S_nz, "depth_sum": np.bincount(inv, weights=depth, minlength=K),
            "stage": np.array(stage), "label_key": np.array(label_key), "seed": np.array(seed),
            "source_sha256": np.array(sha256_file(path)), "library_note": np.array(lib_note),
            "median_depth": np.array(float(np.median(depth)))}


def save(stats: dict, out: str | Path) -> None:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f".{out.stem}.partial.npz")
    try:
        np.savez_compressed(tmp, **stats)
        tmp.replace(out)
    finally:
        if tmp.exists():
            tmp.unlink()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("h5ad")
    ap.add_argument("stage")
    ap.add_argument("out")
    ap.add_argument("--label-key", default="celltype")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    st = compute(a.h5ad, a.stage, a.label_key, a.seed)
    save(st, a.out)
    print(a.stage, "cells", int(st["n"].sum()), "genes", len(st["genes"]), "groups", len(st["keys"]),
          "median depth", float(st["median_depth"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
