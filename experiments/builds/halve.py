"""Keep a uniform random half of a built file's cells: the operation that made the heart file uploaded on 2026-11-18
(docs/PLAN.md decision log 2026-11-18 (j), departure 3), as a committed tool.

  halve.py IN OUT [--seed 1000001]

The rows kept are numpy.sort(numpy.random.default_rng(seed).choice(n_obs, size=n_obs // 2, replace=False)): half of
the rows without replacement, in their order in IN. Nothing else changes: same variables, obs columns, X dtype and X
type; the file is written with anndata's write_h5ad(compression="gzip") and no stamp, so the same input and seed give
the same bytes. The seed is never 0 and never below vec.harness.BUILD_SEED_OFFSET. The tool prints "ok" and the
output's sha256 only, never a cell count (a test board's count is never shown). OUT must not be IN.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys

import anndata as ad
import numpy as np

from vec.harness import BUILD_SEED_OFFSET

DEFAULT_SEED = BUILD_SEED_OFFSET + 1


def half_rows(n_obs: int, seed: int) -> np.ndarray:
    """The rows kept: a uniform random half without replacement, sorted (the file's order kept)."""
    return np.sort(np.random.default_rng(seed).choice(n_obs, size=n_obs // 2, replace=False))


def check_seed(seed: int) -> int:
    if seed < BUILD_SEED_OFFSET:  # BUILD_SEED_OFFSET > 0, so this refuses 0 too
        raise SystemExit(f"--seed must be at least vec.harness.BUILD_SEED_OFFSET ({BUILD_SEED_OFFSET}) and never 0")
    return seed


def halve(src: str, out: str, seed: int = DEFAULT_SEED) -> str:
    """Write the halved file and return its sha256."""
    check_seed(seed)
    if os.path.realpath(src) == os.path.realpath(out):
        raise SystemExit("OUT must not be IN: the input is read whole into memory and would be replaced")
    a = ad.read_h5ad(src)
    idx = half_rows(a.n_obs, seed)
    b = a[idx].copy()
    assert list(b.var_names) == list(a.var_names), "variables changed"
    assert list(b.obs.columns) == list(a.obs.columns), "obs columns changed"
    assert b.X.dtype == a.X.dtype, "X dtype changed"
    assert type(b.X) is type(a.X), "X type changed"
    b.write_h5ad(out, compression="gzip")
    h = hashlib.sha256()
    with open(out, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src", metavar="IN")
    ap.add_argument("out", metavar="OUT")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    a = ap.parse_args(argv)
    check_seed(a.seed)
    digest = halve(a.src, a.out, a.seed)
    print("ok")
    print(digest)


if __name__ == "__main__":
    main(sys.argv[1:])
