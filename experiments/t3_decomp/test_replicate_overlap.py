"""Tests for the note-4 replicate overlap check (experiments/t3_decomp/replicate_overlap.py), on synthetic files."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("replicate_overlap", HERE / "replicate_overlap.py")
ro = importlib.util.module_from_spec(spec)
sys.modules["replicate_overlap"] = ro
spec.loader.exec_module(ro)

GENES = [f"g{i}" for i in range(20)]


def _adata(rows, genes=GENES) -> ad.AnnData:
    C = sparse.csr_matrix(np.asarray(rows, dtype=np.float32))
    a = ad.AnnData(X=C.copy(), layers={"counts": C})
    a.var_names = list(genes)
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    return a


def _pair(n: int = 200, copied: int = 3, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Two independent Poisson count tables (20 genes, about 100 counts a row), the first `copied` rows of the first
    also placed in the second."""
    rng = np.random.default_rng(seed)
    r1, r2 = rng.poisson(5, size=(n, 20)).astype(float), rng.poisson(5, size=(n, 20)).astype(float)
    r2[n - copied:] = r1[:copied]
    return r1, r2


def _counts(a: ad.AnnData, genes=GENES) -> sparse.csr_matrix:
    return ro.counts_on(a, list(genes))


def test_shared_rows_are_found_in_both_directions():
    r1, r2 = _pair()
    rec = ro.overlap(_counts(_adata(r1)), _counts(_adata(r2)))
    assert rec["rows"] == [200, 200]
    assert rec["shared_rows"] == [3, 3]
    assert rec["shared_fraction"] == [0.015, 0.015]
    assert rec["independent"] is False  # 1.5% > 1%


def test_exactly_one_percent_is_still_independent():
    r1, r2 = _pair(copied=2)
    rec = ro.overlap(_counts(_adata(r1)), _counts(_adata(r2)))
    assert rec["shared_fraction"] == [0.01, 0.01]
    assert rec["independent"] is True  # the note says more than 1%


def test_empty_rows_never_count_and_low_count_rows_are_reported_apart():
    r1, r2 = _pair(copied=0)
    low = np.zeros(20)
    low[3] = 2.0  # a two-count row, present in both files
    r1 = np.vstack([r1, np.zeros(20), low])
    r2 = np.vstack([r2, np.zeros(20), low])
    rec = ro.overlap(_counts(_adata(r1)), _counts(_adata(r2)))
    assert rec["empty_rows"] == [1, 1]
    assert rec["shared_rows"] == [1, 1]  # the low-count row only
    assert rec["shared_rows_min_library"] == [0, 0]


def test_gene_order_is_aligned_before_comparing():
    r1, r2 = _pair()
    perm = np.random.default_rng(1).permutation(20)
    b = _adata(r2[:, perm], genes=[GENES[i] for i in perm])
    rec = ro.overlap(_counts(_adata(r1)), _counts(b))
    assert rec["shared_rows"] == [3, 3]


def test_a_hash_match_alone_is_not_a_shared_row(monkeypatch):
    r1, r2 = _pair()
    monkeypatch.setattr(ro, "row_hashes", lambda M: np.zeros(M.shape[0], dtype=np.uint64))  # every row collides
    rec = ro.overlap(_counts(_adata(r1)), _counts(_adata(r2)))
    assert rec["shared_rows"] == [3, 3]


def test_a_row_repeated_in_one_file_counts_once_per_row():
    r1, r2 = _pair(copied=1)
    r2[0] = r1[0]  # the same rep1 row now sits twice in rep2
    rec = ro.overlap(_counts(_adata(r1)), _counts(_adata(r2)))
    assert rec["shared_rows"] == [1, 2]


def test_main_writes_the_record_and_refuses_bad_paths(tmp_path, capsys):
    r1, r2 = _pair()
    extra = np.hstack([r2, np.ones((200, 1))])  # rep2 holds one more gene; only shared genes are compared
    p1, p2 = tmp_path / "rep1.h5ad", tmp_path / "rep2.h5ad"
    _adata(r1).write_h5ad(p1)
    _adata(extra, genes=GENES + ["extra"]).write_h5ad(p2)
    out = tmp_path / "overlap.json"
    assert ro.main([str(p1), str(p2), "--out", str(out)]) == 0
    rec = json.loads(out.read_text())
    assert rec["genes_compared"] == 20
    assert rec["shared_rows"] == [3, 3] and rec["independent"] is False
    assert [f["genes"] for f in rec["files"]] == [20, 21]
    assert "NOT independent" in capsys.readouterr().out
    missing = tmp_path / "missing.h5ad"
    assert ro.main([str(p1), str(missing), "--out", str(tmp_path / "x.json")]) == 2
    assert ro.main([str(p1), str(p1), "--out", str(tmp_path / "y.json")]) == 2
    assert not (tmp_path / "x.json").exists() and not (tmp_path / "y.json").exists()
