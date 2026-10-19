#!/usr/bin/env python
"""Replicate independence check before the Gata4 run (rules.yaml note 4, lead, 2026-10-19).

Note 4: before the gata4 run, count identical count rows shared by the two replicate files; if they share more than 1%
of either file's rows, the borrowed-replicate rules (D2, D2u, D4-D6, D8, W1) are recorded as not independent.

How (fixed 2026-10-19, before any Gata4 file was downloaded): both files as released (no records dropped: a record
present in both files is what this looks for), layers["counts"] on the genes both files hold, in one order, compared
exactly (64-bit row hashes from vec.spatial.row_hashes, then an exact comparison of every candidate pair). Empty rows
are never counted as shared. The verdict uses every non-empty row, as the note says, against each file's total row
count; the same count among rows with at least MIN_LIBRARY total counts (vec.spatial's threshold for chance
coincidences) is reported next to it, so the record shows whether low-count rows drive the verdict.

Writes a JSON record (aggregates only: counts, fractions, file hashes, the commit) and prints one line. Run on the DGX
from the repository root before run_gata4.sh, then bring the record back and commit it next to the decision file:

    uv run python experiments/t3_decomp/replicate_overlap.py REP1.h5ad REP2.h5ad \\
        --out runs/t3_decomp/replicate_overlap.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from vec.spatial import row_hashes

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
THRESHOLD = 0.01  # note 4: more than 1% of either file's rows
MIN_LIBRARY = 20  # vec.spatial: rows with fewer total counts can coincide by chance
BORROWED = "D2, D2u, D4-D6, D8 and W1"


def sha256(path, chunk=1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def git_state() -> dict:
    def run(*cmd):
        try:
            return subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception as e:  # noqa: BLE001
            return f"unavailable: {e!r}"
    return {"commit": run("git", "rev-parse", "HEAD"),
            "status": run("git", "status", "--porcelain", "--", "src/vec", "experiments/t3_decomp")}


def counts_on(a: ad.AnnData, genes: list[str]) -> sparse.csr_matrix:
    """layers["counts"] on `genes`, in that order, as canonical float64 CSR (equal rows have equal storage)."""
    if "counts" not in a.layers:
        raise KeyError("layers['counts'] is required (note 4 compares count rows)")
    C = sparse.csr_matrix(sparse.csr_matrix(a.layers["counts"], dtype=np.float64)[:, a.var_names.get_indexer(genes)])
    C.sum_duplicates()
    C.eliminate_zeros()
    return C


def _equal(A: sparse.csr_matrix, i: int, B: sparse.csr_matrix, j: int) -> bool:
    a, b = slice(A.indptr[i], A.indptr[i + 1]), slice(B.indptr[j], B.indptr[j + 1])
    return bool(np.array_equal(A.indices[a], B.indices[b]) and np.array_equal(A.data[a], B.data[b]))


def shared_rows(A: sparse.csr_matrix, B: sparse.csr_matrix) -> tuple[np.ndarray, np.ndarray]:
    """Masks of the non-empty rows of A equal to some row of B, and of the non-empty rows of B equal to one of A."""
    ha, hb = row_hashes(A), row_hashes(B)
    full_a, full_b = np.diff(A.indptr) > 0, np.diff(B.indptr) > 0
    in_a, in_b = np.zeros(A.shape[0], dtype=bool), np.zeros(B.shape[0], dtype=bool)
    common = np.intersect1d(ha[full_a], hb[full_b])
    by_hash: dict[int, list[int]] = {}
    for j in np.flatnonzero(full_b & np.isin(hb, common)):
        by_hash.setdefault(int(hb[j]), []).append(int(j))
    for i in np.flatnonzero(full_a & np.isin(ha, common)):
        for j in by_hash.get(int(ha[i]), ()):
            if _equal(A, int(i), B, j):  # the hash only proposes a pair; equal storage decides
                in_a[i] = in_b[j] = True
    return in_a, in_b


def overlap(A: sparse.csr_matrix, B: sparse.csr_matrix) -> dict:
    """Note 4's counts for two count matrices over the same genes in the same order."""
    n = [int(A.shape[0]), int(B.shape[0])]
    if min(n) == 0:
        raise ValueError(f"a replicate file has no rows ({n})")
    in_a, in_b = shared_rows(A, B)
    lib = [np.asarray(M.sum(axis=1)).ravel() for M in (A, B)]
    shared = [int(in_a.sum()), int(in_b.sum())]
    shared_lib = [int((in_a & (lib[0] >= MIN_LIBRARY)).sum()), int((in_b & (lib[1] >= MIN_LIBRARY)).sum())]
    frac = [shared[k] / n[k] for k in (0, 1)]
    return {"rows": n, "empty_rows": [int((np.diff(M.indptr) == 0).sum()) for M in (A, B)],
            "shared_rows": shared, "shared_fraction": [round(f, 6) for f in frac],
            "shared_rows_min_library": shared_lib,
            "shared_fraction_min_library": [round(shared_lib[k] / n[k], 6) for k in (0, 1)],
            "independent": bool(max(frac) <= THRESHOLD)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rep1", type=Path)
    ap.add_argument("rep2", type=Path)
    ap.add_argument("--out", type=Path, default=Path("runs/t3_decomp/replicate_overlap.json"))
    args = ap.parse_args(argv)
    for p in (args.rep1, args.rep2):
        if not p.is_file():
            print(f"replicate_overlap.py: {p} not found; nothing was run", file=sys.stderr)
            return 2
    if args.rep1.resolve() == args.rep2.resolve():
        print("replicate_overlap.py: both paths name the same file", file=sys.stderr)
        return 2
    a, b = ad.read_h5ad(args.rep1), ad.read_h5ad(args.rep2)
    in_b = set(b.var_names)
    genes = [g for g in a.var_names if g in in_b]
    if not genes:
        print("replicate_overlap.py: the files share no genes", file=sys.stderr)
        return 2
    rec = {"rule": "rules.yaml note 4 (lead, 2026-10-19)", "date": datetime.now(UTC).isoformat(timespec="seconds"),
           "git": git_state(), "script_sha256": sha256(Path(__file__)), "threshold": THRESHOLD,
           "min_library": MIN_LIBRARY,
           "files": [{"path": str(p), "sha256": sha256(p), "genes": int(x.n_vars)}
                     for p, x in ((args.rep1, a), (args.rep2, b))],
           "genes_compared": len(genes), **overlap(counts_on(a, genes), counts_on(b, genes))}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rec, indent=2) + "\n")
    s, n, f, sl = rec["shared_rows"], rec["rows"], rec["shared_fraction"], rec["shared_rows_min_library"]
    verdict = ("independent" if rec["independent"]
               else f"NOT independent: record {BORROWED} as not independent")
    print(f"note 4: {s[0]} of {n[0]} rep1 rows ({f[0]:.2%}) and {s[1]} of {n[1]} rep2 rows ({f[1]:.2%}) also occur in "
          f"the other file ({sl[0]} and {sl[1]} with >= {MIN_LIBRARY} counts): {verdict}; wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
