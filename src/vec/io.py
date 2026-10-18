"""Loading stages, aligning to a board's gene panel, writing and checking submission files.

The checks mirror what the portal's free format check reports (genes, order, cell count, coordinates, missing values,
structure) plus two convention warnings the scorer only reports: non-normalised scale and broadcast populations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from vec.boards import Board

COORD_KEY = "spatial_3D"
MAX_FILE_BYTES = 1200 * 1024 * 1024
TARGET_LIBRARY_SIZE = 1e4  # released data is log(1 + 1e4 * c / sum c)


def dense(X) -> np.ndarray:
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def read_stage(path: str | Path, genes: list[str] | None = None) -> ad.AnnData:
    """Read an h5ad; if `genes` is given, subset and reorder columns to that panel (all genes must be present)."""
    a = ad.read_h5ad(path)
    if genes is not None:
        missing = [g for g in genes if g not in a.var_names]
        if missing:
            raise KeyError(f"{path}: {len(missing)} panel genes missing, e.g. {missing[:5]}")
        a = a[:, genes].copy()
    return a


def make_submission(X, genes: list[str], coords: np.ndarray | None = None) -> ad.AnnData:
    """Build a submission AnnData: float32 non-negative X in panel order, optional (n, 3) coordinates."""
    X = dense(X).astype(np.float32, copy=False)
    if X.shape[1] != len(genes):
        raise ValueError(f"X has {X.shape[1]} genes, panel has {len(genes)}")
    X = np.clip(X, 0.0, None)
    a = ad.AnnData(X=X)
    a.var_names = list(genes)
    a.obs_names = [f"cell_{i}" for i in range(a.n_obs)]
    if coords is not None:
        coords = np.asarray(coords, dtype=np.float32)
        if coords.shape != (a.n_obs, 3):
            raise ValueError(f"coords must be ({a.n_obs}, 3), got {coords.shape}")
        a.obsm[COORD_KEY] = coords
    return a


def write_submission(a: ad.AnnData, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path, compression="gzip")
    return path


@dataclass
class CheckReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    n_cells: int = 0

    @property
    def ok(self) -> bool:
        return not self.errors


def check_submission(path: str | Path, board: Board) -> CheckReport:
    """Validate a prediction file against a board spec without scoring it."""
    path = Path(path)
    rep = CheckReport()
    if path.stat().st_size > MAX_FILE_BYTES:
        rep.errors.append(f"file is {path.stat().st_size / 2**20:.0f} MB, limit is 1200 MB")
    a = ad.read_h5ad(path)
    rep.n_cells = a.n_obs

    genes = [str(g) for g in a.var_names]
    if genes != board.genes:
        if set(genes) == set(board.genes) and len(genes) == len(board.genes):
            rep.errors.append("genes match the panel but are in a different order")
        else:
            missing = [g for g in board.genes if g not in set(genes)]
            extra = [g for g in genes if g not in set(board.genes)]
            rep.errors.append(f"gene panel mismatch: {len(missing)} missing (e.g. {missing[:5]}), "
                              f"{len(extra)} extra (e.g. {extra[:5]})")

    if a.n_obs < board.min_cells or (board.max_cells is not None and a.n_obs > board.max_cells):
        rep.errors.append(f"{a.n_obs} cells outside [{board.min_cells}, {board.max_cells}]")

    X = dense(a.X)
    if not np.isfinite(X).all():
        rep.errors.append("X contains NaN or infinite values")
    elif (X < 0).any():
        rep.errors.append("X contains negative values (data centred upstream?)")
    else:
        lib = np.median(np.expm1(np.clip(X, 0, 50)).sum(1))
        if board.task == "T1" and not 0.5 * TARGET_LIBRARY_SIZE <= lib <= 2 * TARGET_LIBRARY_SIZE:
            rep.warnings.append(f"median implied library size {lib:.0f}, expected about 1e4")
        if X.shape[0] > 1 and float(X.var(0).mean()) < 1e-8:
            rep.warnings.append("all cells are identical: not a population")

    if board.needs_coords:
        C = a.obsm.get(COORD_KEY)
        if C is None:
            rep.errors.append(f"obsm['{COORD_KEY}'] is required for {board.task}")
        else:
            C = np.asarray(C)
            if C.ndim != 2 or C.shape[0] != a.n_obs or C.shape[1] < 3:
                rep.errors.append(f"obsm['{COORD_KEY}'] must be (n_cells, 3), got {C.shape}")
            elif not np.isfinite(C[:, :3]).all():
                rep.errors.append(f"obsm['{COORD_KEY}'] contains NaN or infinite values")
            else:
                n_dup = C.shape[0] - np.unique(C[:, :3], axis=0).shape[0]
                if n_dup:
                    rep.warnings.append(f"{n_dup} cells repeat the exact coordinates of another cell (copied "
                                        "positions? kNN and shape metrics see them as zero-distance neighbours)")
    if "celltype" in a.obs:
        rep.warnings.append("obs['celltype'] is ignored by the scorer; it is safe to drop")
    return rep
