"""Load released stages into one uniform form.

Every loaded stage has: CSR float32 `.X` in log1p CP10k with explicit zeros removed, unique obs names, an integer
`layers["counts"]` (taken from the MERFISH files, or recovered exactly from `.X` for the T1 scRNA), and coarse lineage
columns from `configs/lineage_map.yaml` (`obs["group"]`, `obs["core"]`, `obs["mesp1"]`).

Prediction pipelines load with `drop_duplicates=True` (copied and interpolated records removed); the board emulator
keeps the default False, because the scorer's truth halves come from the files as released (see `load_stage`).
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from vec.boards import REPO_ROOT
from vec.io import COORD_KEY
from vec.spatial import find_duplicate_records

LINEAGE_MAP = REPO_ROOT / "configs" / "lineage_map.yaml"
TARGET_SUM = 1e4


@lru_cache
def lineage_map(path: Path = LINEAGE_MAP) -> dict:
    spec = yaml.safe_load(Path(path).read_text())
    spec.setdefault("file_overrides", {})
    targets = set(spec["labels"].values()).union(*(set(o.values()) for o in spec["file_overrides"].values()))
    undefined = sorted(targets - set(spec["groups"]))
    if undefined:
        raise ValueError(f"labels map to undefined groups: {undefined}")
    for name, g in spec["groups"].items():
        if "merged" not in g:
            raise ValueError(f"group {name} has no `merged` entry (use null for no cross-stage meaning)")
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


def add_lineage(a: ad.AnnData, label_key: str = "celltype", source: str | None = None) -> list[str]:
    """Annotate obs with group/core/mesp1/merged from the lineage map; returns labels that were not in the map.

    `source` is the file stem; labels listed under the map's `file_overrides` for that stem are remapped first.
    `merged` is "Untracked" for groups without a stable cross-stage meaning.
    """
    spec = lineage_map()
    mapping = {**spec["labels"], **spec["file_overrides"].get(source, {})}
    labels = a.obs[label_key].astype(str) if label_key in a.obs else np.array(["Unknown"] * a.n_obs)
    groups = [mapping.get(lab, "Unknown") for lab in labels]
    a.obs["group"] = np.asarray(groups)
    a.obs["core"] = np.array([spec["groups"][g]["core"] for g in groups], dtype=bool)
    a.obs["mesp1"] = np.array([spec["groups"][g]["mesp1"] for g in groups])
    a.obs["merged"] = np.array([spec["groups"][g]["merged"] or "Untracked" for g in groups])
    return sorted({lab for lab in set(labels) if lab not in mapping})


def _duplicate_filter(a: ad.AnnData) -> tuple[np.ndarray, dict]:
    """Rows to keep and the `uns["qc"]` record for `load_stage(..., drop_duplicates=True)`.

    Two kinds of synthetic record are dropped (experiments/qc): copies, i.e. identical count rows in the same or an
    adjacent z level, whether at the same xy or moved (`find_duplicate_records` with `xy_tol=inf`; identical rows of
    >= 20 counts never occur by chance in the released files), and interpolated records, i.e. rows with non-integer
    counts, which equal the average of two cells in neighbouring sections.
    """
    qc: dict = {"drop_duplicates": True, "n_input": int(a.n_obs)}
    if "counts" not in a.layers or COORD_KEY not in a.obsm:
        return np.ones(a.n_obs, dtype=bool), {
            **qc, "n_dropped": 0, "note": f"no counts layer or obsm['{COORD_KEY}']: duplicate records need both"}
    copies, rep = find_duplicate_records(a, xy_tol=np.inf)
    C = sparse.csr_matrix(a.layers["counts"])
    fractional = np.abs(C.data - np.rint(C.data)) > 1e-6
    interpolated = np.zeros(a.n_obs, dtype=bool)
    interpolated[np.repeat(np.arange(a.n_obs), np.diff(C.indptr))[fractional]] = True
    interpolated &= ~copies
    drop = copies | interpolated
    pairs = pd.DataFrame(rep["by_level_pair"], columns=["kept_z", "dropped_z", "gap", "n"]).astype(
        {"kept_z": np.float64, "dropped_z": np.float64, "gap": np.int64, "n": np.int64})
    qc.update(n_dropped=int(drop.sum()), n_copies=int(copies.sum()), n_interpolated=int(interpolated.sum()),
              dropped_rows=np.flatnonzero(drop).astype(np.int64), by_gap=rep["by_gap"], by_level_pair=pairs,
              below_min_library=rep["below_min_library"], max_xy_offset=rep["max_xy_offset"], xy_tol="any",
              min_library=rep["min_library"], max_level_gap=rep["max_level_gap"])
    return ~drop, qc


def load_stage(path: str | Path, genes: list[str] | None = None, lineage: bool = True,
               drop_duplicates: bool = False) -> ad.AnnData:
    """Read a released h5ad into the uniform form described in the module docstring.

    `drop_duplicates=True` first removes synthetic records (copies at the same or a moved xy in the same or an adjacent
    z level, and interpolated records with non-integer counts; see `_duplicate_filter`), keeping the first occurrence in
    file order, and records what was dropped in `uns["qc"]` (`dropped_rows` are positions in the file). Prediction
    pipelines pass True. The default stays False because the scorer draws its truth halves from the organisers' files
    as released, so the board emulator must load targets as they are. Counts are never truncated: they are int32 when
    every value is an integer and float32 otherwise. Files without coordinates or a counts layer (T1 scRNA) are
    returned whole, with a note in `uns["qc"]`.
    """
    a = ad.read_h5ad(path)
    qc = None
    if drop_duplicates:
        keep, qc = _duplicate_filter(a)
        if not keep.all():
            a = a[keep].copy()
    X = sparse.csr_matrix(a.X, dtype=np.float32)
    X.eliminate_zeros()
    if "counts" in a.layers:
        counts = sparse.csr_matrix(a.layers["counts"], dtype=np.float64)
        if np.all(np.abs(counts.data - np.rint(counts.data)) <= 1e-6):
            counts.data = np.rint(counts.data)
            counts = counts.astype(np.int32)
        else:
            counts = counts.astype(np.float32)
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
        unmapped = add_lineage(out, source=Path(path).name.removesuffix(".gz").removesuffix(".h5ad"))
        out.uns["unmapped_labels"] = unmapped
    if qc is not None:
        out.uns["qc"] = qc
    return out
