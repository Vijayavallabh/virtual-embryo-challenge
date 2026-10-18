"""Load released stages into one uniform form.

Every loaded stage has: CSR float32 `.X` in log1p CP10k with explicit zeros removed, unique obs names, an integer
`layers["counts"]` (taken from the MERFISH files, or recovered exactly from `.X` for the T1 scRNA), and coarse lineage
columns from `configs/lineage_map.yaml` (`obs["group"]`, `obs["core"]`, `obs["mesp1"]`).
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import anndata as ad
import numpy as np
import yaml
from scipy import sparse

from vec.boards import REPO_ROOT

LINEAGE_MAP = REPO_ROOT / "configs" / "lineage_map.yaml"
TARGET_SUM = 1e4


@lru_cache
def lineage_map(path: Path = LINEAGE_MAP) -> dict:
    spec = yaml.safe_load(Path(path).read_text())
    undefined = sorted(set(spec["labels"].values()) - set(spec["groups"]))
    if undefined:
        raise ValueError(f"labels map to undefined groups: {undefined}")
    return spec


def recover_counts(X) -> tuple[sparse.csr_matrix, np.ndarray, float]:
    """Integer counts from `.X = log1p(1e4 * c / N)`.

    Each cell's smallest non-zero value corresponds to a count of 1, so N = 1e4 / expm1(min non-zero) and
    c = expm1(x) * N / 1e4. Returns (counts, library sizes N, largest distance of any recovered count from an integer),
    which should be ~0 when every cell has at least one gene with a single count.
    """
    X = sparse.csr_matrix(X, dtype=np.float64)
    X.eliminate_zeros()
    n = X.shape[0]
    nnz_per_row = np.diff(X.indptr)
    if (nnz_per_row == 0).any():
        raise ValueError(f"{int((nnz_per_row == 0).sum())} cells have no expressed genes")
    row_min = np.minimum.reduceat(X.data, X.indptr[:-1])
    lib = TARGET_SUM / np.expm1(row_min)
    vals = np.expm1(X.data) * np.repeat(lib, nnz_per_row) / TARGET_SUM
    rounded = np.rint(vals)
    resid = float(np.abs(vals - rounded).max()) if vals.size else 0.0
    counts = sparse.csr_matrix((rounded.astype(np.int32), X.indices.copy(), X.indptr.copy()), shape=X.shape)
    assert counts.shape[0] == n
    return counts, lib, resid


def log1p_cp10k(counts) -> sparse.csr_matrix:
    """log(1 + 1e4 * c / sum c) per cell, as float32 CSR (cells with zero counts stay zero)."""
    C = sparse.csr_matrix(counts, dtype=np.float64)
    lib = np.asarray(C.sum(1)).ravel()
    scale = np.divide(TARGET_SUM, lib, out=np.zeros_like(lib), where=lib > 0)
    C = sparse.diags(scale) @ C
    C.data = np.log1p(C.data)
    return C.astype(np.float32).tocsr()


def add_lineage(a: ad.AnnData, label_key: str = "celltype") -> list[str]:
    """Annotate obs with group/core/mesp1 from the lineage map; returns labels that were not in the map."""
    spec = lineage_map()
    labels = a.obs[label_key].astype(str) if label_key in a.obs else np.array(["Unknown"] * a.n_obs)
    groups = [spec["labels"].get(lab, "Unknown") for lab in labels]
    a.obs["group"] = np.asarray(groups)
    a.obs["core"] = np.array([spec["groups"][g]["core"] for g in groups], dtype=bool)
    a.obs["mesp1"] = np.array([spec["groups"][g]["mesp1"] for g in groups])
    return sorted({lab for lab in set(labels) if lab not in spec["labels"]})


def load_stage(path: str | Path, genes: list[str] | None = None, lineage: bool = True) -> ad.AnnData:
    """Read a released h5ad into the uniform form described in the module docstring."""
    a = ad.read_h5ad(path)
    X = sparse.csr_matrix(a.X, dtype=np.float32)
    X.eliminate_zeros()
    if "counts" in a.layers:
        counts = sparse.csr_matrix(a.layers["counts"]).astype(np.int32)
        counts.eliminate_zeros()
    else:
        counts, _, resid = recover_counts(X)
        if resid > 1e-2:
            raise ValueError(f"{path}: counts not recoverable from X (max residual {resid:.3g})")
    keep_obs = [c for c in a.obs.columns if c in ("celltype", "cm_celltype", "genotype", "sample_id", "condition")]
    out = ad.AnnData(X=X, obs=a.obs[keep_obs].copy(), var=a.var[[]].copy(), obsm=dict(a.obsm))
    out.layers["counts"] = counts
    out.obs_names_make_unique()
    if genes is not None:
        missing = [g for g in genes if g not in out.var_names]
        if missing:
            raise KeyError(f"{path}: {len(missing)} panel genes missing, e.g. {missing[:5]}")
        out = out[:, genes].copy()
    if lineage:
        unmapped = add_lineage(out)
        out.uns["unmapped_labels"] = unmapped
    return out
