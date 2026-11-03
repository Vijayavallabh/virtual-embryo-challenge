#!/usr/bin/env python
"""Detection depth of every released stage (2026-11-03; a diagnostic, not an adoption result): median counts per cell
(each cell's library size, the N its released values were normalised by) and median genes detected, duplicates
dropped on the MERFISH files. Behind tries E7, E8 and H3 of experiments/harness/TRIES.md.

Usage (DGX, repo root, BLAS threads capped as docs/COMPUTE.md says):
  uv run python experiments/harness/diag_depth.py --out experiments/harness/results/diag_depth
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from vec.boards import REPO_ROOT
from vec.prep import load_stage
from vec.t2_build import library_sizes

DATA = REPO_ROOT / "data" / "raw"
MERFISH = ("T2/E6.75", "T2/E7.25", "T2/E7.5", "T2/E8.0", "T2/E8.25_late", "T2/E8.5", "T2/E8.75", "T2/E9.5",
           "T2/E10.5", "T3/E8.75_gata4_ko_1", "T3/E8.75_gata4_ko_2", "T3/E9.5_mab21l2_ko")
SCRNA = ("T1/E8.5_RNA", "T1/E9.5_RNA", "T1/E10.5_RNA")


def _summary(lib: np.ndarray, det: np.ndarray, n: int) -> dict:
    return {"cells": int(n), "median_counts_per_cell": float(np.median(lib)),
            "q10_q90_counts": [float(np.percentile(lib, 10)), float(np.percentile(lib, 90))],
            "median_genes_detected": float(np.median(det))}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    res = {}
    for f in MERFISH:
        st = load_stage(DATA / f"{f}.h5ad", lineage=False, drop_duplicates=True)
        lib = library_sizes(st.X, st.layers["counts"])
        res[f] = _summary(lib, np.asarray((st.layers["counts"] > 0).sum(1)).ravel(), st.n_obs)
    for f in SCRNA:  # released as log(1 + 1e4 c / N); a nucleus's smallest value is a single count
        X = sparse.csr_matrix(ad.read_h5ad(DATA / f"{f}.h5ad").X, dtype=np.float64)
        X.eliminate_zeros()
        lib = 1e4 / np.expm1(np.minimum.reduceat(X.data, X.indptr[:-1]))
        res[f] = _summary(lib, np.diff(X.indptr), X.shape[0])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=REPO_ROOT).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                                text=True, cwd=REPO_ROOT).stdout.strip())
    doc = {"label": "diagnostic, not an adoption result", "commit": commit, "tree_clean": not dirty, "stages": res}
    (out / "diag_depth.json").write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()
