#!/usr/bin/env python
"""T1's pool depth record for row R37's pooled base (docs/DATA_POLICY.md 2026-11-04 (f)(3)), from released E10.5 alone.

Row R37 of experiments/harness/TRIES.md (the captain's named exception of 2026-11-24) lets one file, R13's command on
the pooled base (build_test.py t1 --g-mom --base pooled --feng-groups recomputed --captain-ruling R37, file
T1_s2rc_pooled), change T1's depth. This script prints the target-free record the policy's entry (f)(3) asks for: the
median per-cell library size (total counts per cell, vec.t1_bets.median_library) of the pool each base draws from, as a
ratio to the median over E10.5's whole chamber pool (every E10.5 cell whose label configs/t1_chamber_classes.yaml puts
in a chamber category):
  pooled  every E10.5 cell in a chamber category (the pooled base; T1_s2rc_pooled)
  lib3    library 3's E10.5 cells in a chamber category (the dissected-heart library; R13's file and its lines)
The draw of every T1 file takes chamber-category cells only (its shares are vec.harness.t1_chamber_shares), so these are
the cells each base can draw; the pooled pool is the whole chamber pool, so its ratio is 1 by construction and is
printed as a check. A cell's library is the suffix after the last underscore of its name, as build_test.py reads it.

Nothing else is read: no board file (so no test board's count), no built file and no build record. Only ratios of
medians of released cells are printed, no cell number. Run from the repository root (DGX, BLAS threads capped):
  uv run --frozen --no-sync python -I -B experiments/builds/t1_pool_depth_record.py
"""

from __future__ import annotations

import json
import sys

import numpy as np

from vec import harness, t1_bets, t1_scope
from vec.prep import load_stage

E105 = harness.DATA_DIR / "T1" / "E10.5_RNA.h5ad"
SOURCE = "E10.5_RNA"
DISSECTED_LIBRARY = "3"  # build_test.T1_DISSECTED_LIBRARY: E10.5's dissected-heart library
POOLS = {
    "pooled": "every E10.5 cell in a chamber category: the pooled base (row R37's file, T1_s2rc_pooled)",
    "lib3": "library 3's E10.5 cells in a chamber category: the lib3 base (R13's file, T1_s2rc, and its lines)",
}


def libraries(names) -> np.ndarray:
    """Each cell's library: the suffix after the last underscore of its name, as build_test.build_t1 reads it."""
    return np.array([x.rsplit("_", 1)[-1] if "_" in x else "" for x in map(str, names)])


def pool_rows(ref) -> dict[str, np.ndarray]:
    """The rows of `ref` (E10.5) in each pool, and in E10.5's whole chamber pool ("chamber"). Refused unless the cell
    names separate library 3 from at least one other library, as build_test.py refuses its lib3 base."""
    spec = t1_scope.load_spec()
    cats = np.asarray(t1_scope.classify(harness._shallow(ref), source=SOURCE, spec=spec, key_added=None)).astype(str)
    chamber = np.isin(cats, spec.in_scope("chamber"))
    library = libraries(ref.obs_names)
    kinds = set(library.tolist())
    if not all(k.isdigit() for k in kinds) or DISSECTED_LIBRARY not in kinds or len(kinds) < 2:
        sys.exit("the reference's cell names do not separate library 3 from the others (every name must end in "
                 "_<library number>, and library 3 must be one of at least two)")
    rows = {"chamber": np.flatnonzero(chamber), "pooled": np.flatnonzero(chamber),
            "lib3": np.flatnonzero(chamber & (library == DISSECTED_LIBRARY))}
    if any(len(r) == 0 for r in rows.values()):
        sys.exit("a pool has no chamber cell")
    return rows


def record(ref) -> dict:
    """The record: each pool's median per-cell library size over the whole chamber pool's, and the pooled pool's over
    library 3's. Ratios only."""
    rows = pool_rows(ref)
    med = {k: t1_bets.median_library(ref[r]) for k, r in rows.items()}
    if not med["chamber"] > 0:
        sys.exit("E10.5's whole chamber pool has a median library size of zero")
    ratios = {k: med[k] / med["chamber"] for k in POOLS}
    return {
        "record": "docs/DATA_POLICY.md 2026-11-04 (f)(3) for row R37's pooled base: the median per-cell library size "
                  "of the pool each T1 base draws from, as a ratio to E10.5's whole chamber pool",
        "source": f"{SOURCE} as released (data/raw/T1/E10.5_RNA.h5ad); no board, test count, built file or build "
                  "record is read",
        "statistic": "median over the pool's cells of the total counts per cell (vec.t1_bets.median_library)",
        "denominator": "E10.5's whole chamber pool: every E10.5 cell whose label configs/t1_chamber_classes.yaml puts "
                       "in a chamber category",
        "pools": {k: {"cells": POOLS[k], "median_library_over_whole_chamber_pool": ratios[k]} for k in POOLS},
        "pooled_over_lib3": ratios["pooled"] / ratios["lib3"],
    }


def main(argv=None) -> None:
    if argv:
        sys.exit("t1_pool_depth_record.py takes no arguments: it reads released E10.5 alone")
    print(json.dumps(record(load_stage(E105, lineage=False)), indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
