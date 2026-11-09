#!/usr/bin/env python
"""Row N14 of experiments/harness/TRIES.md (all boards; a diagnostic, not a try): the integrity pass of the ten frozen
upload files. Insurance against a silent build-path error that vec check (the schema), the byte-identical rehearsals
(determinism) and N5 (a tail count) would not see: X and coordinates paired wrongly, a column-order slip, a wrong source
stage, a size step not applied as pre-registered.

Target-free (the N1/N2/N5 precedent): it reads the ten frozen files in runs/archive/builds_02abf39, their released
source stages (T1 E10.5; heart E10.5; embryo E7.5 and E8.0; T3 WT E8.75, each loaded as its builder loads it), the
committed size JSONs, build records and stamps, and the board gene panels through vec.harness._candidate_board. It reads
no target, computes no score, reads no anchor, and writes and prints no cell count: every number in its output is a
correlation, a share of a fixed-size sample, a ratio or a pass/fail, and the size step is asserted, never shown.

The checks, per file, with thresholds fixed in the row before the run (RULES below):
  (a) per-gene pseudobulk (mean of X) of the file against its source stage, Pearson and Spearman over the panel, a
      column-order check: both correlations at the floor, 0.95 (T1 and T3 0.90); the embryo's source is the equal-weight
      mean of its two brackets.
  (b) T2 and T3: Moran's I profile of the file against its source on the source's 200 most variable genes (the recipe
      of veckit's morans_I_agreement at the pinned git, k = 15), which is what catches X/coordinate misalignment:
      Pearson of the two profiles at least 0.5, and the file's mean Moran's I at least 0.35 of the source's (a partial
      misalignment shrinks the profile but barely changes its correlation); the embryo takes the bracket that agrees
      best. A shuffled-rows control of the file is reported, deciding nothing.
  (c) X finite and non-negative everywhere; on T1 the implied library sum(expm1(X)) of the 0.1th and 99.9th percentile
      cells within 1e4 x [0.99, 1.01] (T1's X covers all genes, and released and moved cells are log1p CP10k of whole
      counts, so every cell sits at 1e4).
  (d) heart: the file's 3-D RMS radius equals the pre-registered size JSON's forecast for its hypothesis (s1 H1, s2 H2),
      assert-only; the embryo likewise against its raw-mode forecast. The value is never printed or stored.
  (e) provenance: the file's hash against the archive's SHA256SUMS, its stamp and its build record (no build record
      keeps a cell's source rows, so the row-to-cell check of the sharpener has nothing to read and this pass says so);
      and row-level provenance by content on a fixed-size sample (K_SAMPLE draws with replacement, so the sample size
      never depends on a file's size): T1 the share of sampled rows that are exact E10.5 rows, the embryo the share that
      are exact rows of E7.5 or E8.0 and how they split between the two, T3 the share whose coordinates are exactly a
      WT E8.75 cell's, the heart the share whose z is a z level of E10.5.
Any flag means: stop, trace the build path, take it to the captain, and no file is rebuilt or uploaded until the cause
is written down; no threshold moves after the run; the row records the outcome. It runs once on a clean tree.

Usage (DGX, repo root; after the row is committed):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
      uv run --frozen --no-sync python -I -B experiments/research_2026_11/diagnostics.py integrity \
      --out experiments/research_2026_11/results/integrity_02abf39
(the same as `python -I -B experiments/research_2026_11/integrity.py`, whose --archive, --data and --out default to the
registered paths). It refuses to run while its outputs exist and on a tree that is not clean by vec.provenance's rule.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import warnings
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse, stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))

import hygiene  # noqa: E402

from vec import harness, provenance, t2_build  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.io import COORD_KEY  # noqa: E402
from vec.prep import load_stage  # noqa: E402

ARCHIVE = REPO / "runs" / "archive" / "builds_02abf39"
OUT = Path(__file__).resolve().parent / "results" / "integrity_02abf39"
SIZES = {"heart": REPO / "experiments/coords/results/p3/size_heart_E12.5.json",
         "embryo": REPO / "experiments/coords/results/p3/size_embryo_E7.75.json"}
DECISIONS = REPO / "experiments/t3_decomp/results/gata4/gata4_decisions.json"
SOURCES = {"t1": ("T1/E10.5_RNA.h5ad",), "heart": ("T2/E10.5.h5ad",), "embryo": ("T2/E7.5.h5ad", "T2/E8.0.h5ad"),
           "t3": ("T2/E8.75.h5ad",)}

# The row's thresholds, fixed before the run (experiments/harness/TRIES.md N14); no threshold moves after it.
K_SAMPLE = 2000       # rows drawn with replacement per file for the row-level provenance: never a function of the file
SAMPLE_SEED = 20261109
PB_FLOOR = {"t1": 0.90, "heart": 0.95, "embryo": 0.95, "t3": 0.90}
MORAN_K = 15          # the recipe of veckit's morans_I_agreement (T2/metrics.py at the pinned git)
MORAN_GENES = 200
MORAN_FLOOR = 0.5
MORAN_LEVEL_FLOOR = 0.35
LIBRARY_TARGET = 1e4
LIBRARY_BAND = (0.99, 1.01)       # 0.1th and 99.9th percentile of the implied library, over 1e4
LIBRARY_QUANTILES = (0.001, 0.999)
SIZE_TOL = {"heart": 1e-4, "embryo": 1e-3}  # abs(log(RMS radius / forecast)); the embryo's jitter moves a few cells
EXACT_FLOOR = {"t1": 0.05, "embryo": 0.99, "t3": 0.05, "heart": 0.05}  # (e); see the module docstring
EMBRYO_SHARE = (0.40, 0.60)       # share of the embryo's matched rows that come from E7.5 (t = 1/2)
BUILD_SEED = harness.BUILD_SEED_OFFSET + 1
LABEL = ("N14 (experiments/harness/TRIES.md; all boards, a diagnostic, not a try): the integrity pass of the ten "
         "frozen upload files against their released source stages (target-free: no score, no anchor, no cell count)")
RULES = {"pb_floor": PB_FLOOR, "moran": {"k": MORAN_K, "genes": MORAN_GENES, "floor": MORAN_FLOOR,
                                         "level_floor": MORAN_LEVEL_FLOOR},
         "library": {"target": LIBRARY_TARGET, "band": list(LIBRARY_BAND), "quantiles": list(LIBRARY_QUANTILES)},
         "size_tol": SIZE_TOL, "k_sample": K_SAMPLE, "sample_seed": SAMPLE_SEED, "exact_floor": EXACT_FLOOR,
         "embryo_share": list(EMBRYO_SHARE)}
ROW_KEYS = ("rows", "source_rows", "source_row", "row_ids", "source_obs", "source_index", "partners")
CHUNK = 5000
_MISSING = object()


@dataclass(frozen=True)
class Job:
    """One frozen upload file and how it was built (experiments/builds/build_all.sh)."""

    file: str           # path in the archive
    board: str          # board key
    kind: str           # "t1", "heart", "embryo" or "t3"
    entry: str | None = None        # the build record's "entry" (T1, heart, embryo)
    cm_share: float | None = None   # T1's s2: --cm-share
    hypothesis: str | None = None   # heart: the size hypothesis
    branch: str | None = None       # T3
    submission: str | None = None   # T3

    @property
    def record(self) -> str:
        """The build record's path in the archive (build_test.py's <name>.build.json, the T3 builder's sidecar)."""
        return self.file[:-len(".h5ad")] + (".json" if self.kind == "t3" else ".build.json")


# The ten frozen files of experiments/builds/build_all.sh (sets s1, s2, t3alt)
JOBS = (Job("builds/T1_s1.h5ad", "T1:test", "t1", entry="M1"),
        Job("builds/T1_s2.h5ad", "T1:test", "t1", entry="M1", cm_share=0.60),
        Job("builds/T2heart_s1.h5ad", "T2:heart:test_extrap", "heart", entry="B0 s1 (H1)", hypothesis="H1"),
        Job("builds/T2heart_s2.h5ad", "T2:heart:test_extrap", "heart", entry="B0 s2_trim (H2)", hypothesis="H2"),
        Job("builds/T2embryo_s1.h5ad", "T2:embryo:test_interp", "embryo", entry="fgw"),
        Job("builds/T2embryo_s2.h5ad", "T2:embryo:test_interp", "embryo", entry="fgw_shape0"),
        Job("t3_build/allowed_s1.h5ad", "T3:bcatenin", "t3", branch="allowed", submission="s1"),
        Job("t3_build/refused_s1.h5ad", "T3:bcatenin", "t3", branch="refused", submission="s1"),
        Job("t3_build/barred_s1.h5ad", "T3:bcatenin", "t3", branch="barred", submission="s1"),
        Job("t3_build/barred_s2.h5ad", "T3:bcatenin", "t3", branch="barred", submission="s2"))


# ---------------------------------------------------------------------------------------------------------------------
# The checks: pure functions of given arrays
def _rows(X, chunk: int = CHUNK) -> Iterator[np.ndarray]:
    """X (dense or sparse) as dense blocks of at most `chunk` rows, in row order."""
    for i in range(0, X.shape[0], chunk):
        part = X[i:i + chunk]
        yield part.toarray() if sparse.issparse(part) else np.asarray(part)


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(1 << 22):
            h.update(block)
    return h.hexdigest()


def pseudobulk(X) -> np.ndarray:
    """The per-gene mean of X over cells, in float64 (the scorer's pb), accumulated in blocks."""
    total = np.zeros(X.shape[1], dtype=np.float64)
    for part in _rows(X):
        total += part.astype(np.float64).sum(0)
    return total / X.shape[0]


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson's r (nan for a constant input)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return float(np.corrcoef(a, b)[0, 1])


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman's rho (nan for a constant input)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # scipy warns about a constant input; nan fails the check that reads it
        return float(stats.spearmanr(a, b).statistic)


def pb_agreement(pb_file: np.ndarray, pb_source: np.ndarray, floor: float) -> dict:
    """(a): Pearson and Spearman of the two pseudobulks; passes only when both reach `floor` (nan fails)."""
    r, rho = _corr(pb_file, pb_source), _spearman(pb_file, pb_source)
    return {"pearson": r, "spearman": rho, "floor": float(floor), "pass": bool(r >= floor and rho >= floor)}


def knn_graph(coords: np.ndarray, k: int = MORAN_K) -> np.ndarray:
    """The k nearest spatial neighbours of every cell, self dropped (the recipe of veckit's morans_I_profile)."""
    from sklearn.neighbors import NearestNeighbors

    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(np.asarray(coords)).kneighbors(np.asarray(coords))
    return idx[:, 1:]


def top_variance_genes(X, n_genes: int = MORAN_GENES) -> np.ndarray:
    """veckit's morans_I_agreement's gene set (chosen on the truth, here the source stage), as column indices."""
    D = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return np.argsort(-D.var(0))[:n_genes]


def moran_profile(X, idx: np.ndarray, genes: np.ndarray) -> np.ndarray:
    """Moran's I per gene on a kNN weight graph (veckit's morans_I_profile at the pinned git), in blocks of rows so that
    a large cloud stays within memory: I_g = sum_i z_ig lag_ig / sum_i z_ig^2, with lag_ig the mean of z over i's k
    neighbours (row-normalised weights)."""
    M = X[:, genes]
    Z = (M.toarray() if sparse.issparse(M) else np.asarray(M)).astype(np.float64)
    Z = Z - Z.mean(0)
    lag = np.empty_like(Z)
    for i in range(0, Z.shape[0], CHUNK):
        lag[i:i + CHUNK] = Z[idx[i:i + CHUNK]].mean(1)
    num, den = (Z * lag).sum(0), (Z**2).sum(0) + 1e-12
    return num / den


def moran_check(profile_file: np.ndarray, profile_source: np.ndarray, profile_shuffled: np.ndarray) -> dict:
    """(b)'s readings for one source: the Pearson of the profiles, the level ratio and the shuffled-rows control."""
    level_file, level_source = float(profile_file.mean()), float(profile_source.mean())
    ratio = level_file / level_source if level_source > 0 else float("nan")
    return {"agreement": _corr(profile_file, profile_source), "level_ratio": ratio,
            "shuffled_control_agreement": _corr(profile_shuffled, profile_source),
            "file_mean_I": level_file, "source_mean_I": level_source}


def x_sanity(X) -> dict:
    """(c): whether every value of X is finite and non-negative, in blocks."""
    finite, nonneg = True, True
    for part in _rows(X):
        finite = finite and bool(np.isfinite(part).all())
        nonneg = nonneg and not bool((part < 0).any())
    return {"finite": finite, "non_negative": nonneg, "pass": finite and nonneg}


def implied_library(X) -> dict:
    """(c), T1: the per-cell implied library sum(expm1(X)) over 1e4, its median and its registered quantiles; passes
    when those quantiles lie inside LIBRARY_BAND. Cells that are not finite fail the check above and pass nothing here:
    a non-finite library gives nan quantiles, so this reports failure."""
    lib = np.concatenate([np.expm1(part.astype(np.float64)).sum(1) for part in _rows(X)]) / LIBRARY_TARGET
    lo, hi = np.quantile(lib, list(LIBRARY_QUANTILES))
    return {"median_over_1e4": float(np.median(lib)), "p0.1_over_1e4": float(lo), "p99.9_over_1e4": float(hi),
            "band": list(LIBRARY_BAND), "pass": bool(LIBRARY_BAND[0] <= lo and hi <= LIBRARY_BAND[1])}


def row_digests(X) -> set[bytes]:
    """16-byte digests of X's rows as float32 bytes (any sparsity): a row is matched by digest, never stored."""
    out: set[bytes] = set()
    for part in _rows(X, 2000):
        for row in np.asarray(part, dtype=np.float32):
            out.add(hashlib.blake2b(row.tobytes(), digest_size=16).digest())
    return out


def exact_row_mask(rows, digests: set[bytes]) -> np.ndarray:
    """Which of the sampled rows (dense float32) are byte-identical to a row digested in `digests`."""
    return np.array([hashlib.blake2b(np.asarray(r, dtype=np.float32).tobytes(), digest_size=16).digest() in digests
                     for r in rows], dtype=bool)


def coord_keys(coords: np.ndarray) -> set[bytes]:
    """Exact-match keys of (n, 3) float32 positions."""
    return {row.tobytes() for row in np.asarray(coords, dtype=np.float32)}


def in_keys(coords: np.ndarray, keys: set[bytes]) -> np.ndarray:
    return np.array([row.tobytes() in keys for row in np.asarray(coords, dtype=np.float32)], dtype=bool)


def size_matches(coords, kind: str, size_json, hypothesis: str | None, tol: float) -> bool:
    """(d), assert-only: whether the file's 3-D RMS radius equals the pre-registered size JSON's forecast (the heart:
    the scope hypothesis's full-cloud forecast in the pre-registered xy_grow_z_hold mode; the embryo: its raw-mode
    forecast) within `tol` on the log scale. Returns only True or False, so no size can reach an output."""
    size, _ = t2_build.load_size_result(size_json)
    if size.get("valid") is not True:
        return False
    pre = size.get("preregistered") or {}
    if kind == "heart":
        if pre.get("z_mode") != "xy_grow_z_hold":
            return False
        target = ((pre.get("scope") or {}).get(hypothesis) or {}).get("pred_rms_full")
    else:
        if pre.get("z_mode") != "raw":
            return False
        target = pre.get("pred_rms_full")
    if isinstance(target, bool) or not isinstance(target, int | float) or not np.isfinite(target) or target <= 0:
        return False
    rms = t2_build.rms_radius(np.asarray(coords, dtype=np.float64))
    return bool(abs(np.log(rms / float(target))) <= tol)


def source_row_fields(obj, path: str = "") -> list[str]:
    """JSON paths in a build record that hold per-cell source rows: a key named as one (ROW_KEYS) whose value is a list
    of at least 10 whole numbers or strings. (The builders write none: the sharpener's row-to-cell check needs such a
    field, so a build record without one is said to have nothing to read.)"""
    out: list[str] = []
    if isinstance(obj, Mapping):
        for k, v in obj.items():
            here = f"{path}.{k}" if path else str(k)
            if str(k) in ROW_KEYS and isinstance(v, list) and len(v) >= 10 and all(isinstance(x, int | str) for x in v):
                out.append(here)
            out += source_row_fields(v, here)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += source_row_fields(v, f"{path}[{i}]")
    return out


def _lookup(obj, dotted: str):
    """obj's value at a dotted path, or _MISSING."""
    for part in dotted.split("."):
        if not isinstance(obj, Mapping) or part not in obj:
            return _MISSING
        obj = obj[part]
    return obj


def record_expectations(job: Job, file_sha: str, source_sha: Mapping[str, str], decisions_sha: str,
                        union_sha: str | None) -> dict:
    """What the build record of `job` should say, dotted path -> value: the builders' own fields, by name only."""
    if job.kind == "t3":
        return {"branch": job.branch, "submission": job.submission, "inputs.board.key": job.board,
                "inputs.n_cells": "max_cells", "inputs.seed": BUILD_SEED,
                "inputs.wt.sha256": source_sha["T2/E8.75.h5ad"],
                "output.sha256": file_sha, "output.check.ok": True, "output.n_cells": "max_cells",
                "decision_file.sha256": decisions_sha, "hashes.mismatch": []}
    expected = {"entry": job.entry, "board": job.board, "dry_run": False, "seed": BUILD_SEED, "n_cells": "max_cells"}
    if job.kind == "t1":
        expected |= {"hold_unbracketed": True, "herd": False, "feng_hygiene.sha256": union_sha}
    if job.cm_share is not None:
        expected["cm_share_of_heart"] = job.cm_share
    if job.hypothesis is not None:  # the heart's size step: the scope hypothesis the build applied
        expected["size.hypothesis"] = job.hypothesis
    return expected


def check_record(job: Job, archive: Path, expected: Mapping[str, object]) -> dict:
    """(e): the build record's fields against `expected`. A field the record does not hold is listed under missing (an
    older record), one it holds with another value is a problem; only the field's name is reported, never its value."""
    path = Path(archive) / job.record
    out = {"path": job.record, "present": path.is_file(), "problems": [], "missing": [], "source_row_fields": []}
    if not out["present"]:
        return out
    rec = json.loads(path.read_text())
    out["source_row_fields"] = source_row_fields(rec)
    for field, want in expected.items():
        have = _lookup(rec, field)
        if have is _MISSING:
            out["missing"].append(field)
        elif not _equal(have, want):
            out["problems"].append(field)
    return out


def _equal(have, want) -> bool:
    """Equality of a record's field and its expected value; floats within 1e-9."""
    if isinstance(want, float) and isinstance(have, int | float) and not isinstance(have, bool):
        return bool(abs(float(have) - want) <= 1e-9)
    return bool(have == want)


def check_stamp(path: Path, file_sha: str) -> dict:
    """(e): the file's build stamp (vec.provenance): its sha256 against the file, its commit and clean flag."""
    sp = provenance.stamp_path(path)
    if not sp.is_file():
        return {"present": False, "sha256_matches": None, "commit": None, "clean": None}
    try:
        s = provenance.read_stamp(sp)
    except (OSError, ValueError) as e:
        return {"present": True, "error": str(e), "sha256_matches": False, "commit": None, "clean": None}
    return {"present": True, "sha256_matches": s["sha256"] == file_sha, "commit": s["commit"], "clean": s["clean"]}


# ---------------------------------------------------------------------------------------------------------------------
# Sources and files
@dataclass
class Source:
    """A source stage on a board's panel, as its builder loads it, and what this pass caches from it."""

    rel: str
    X: sparse.csr_matrix                    # float32, cells x panel genes, panel order
    coords: np.ndarray | None               # float32 (n, 3)
    sha256: str
    _pb: np.ndarray | None = None
    _moran: tuple[np.ndarray, np.ndarray] | None = None
    _digests: set[bytes] | None = None
    _coord_keys: set[bytes] | None = None

    @property
    def pb(self) -> np.ndarray:
        if self._pb is None:
            self._pb = pseudobulk(self.X)
        return self._pb

    def moran(self) -> tuple[np.ndarray, np.ndarray]:
        """(the source's 200 most variable genes, its Moran's I profile on them)."""
        if self._moran is None:
            genes = top_variance_genes(self.X)
            self._moran = (genes, moran_profile(self.X, knn_graph(self.coords), genes))
        return self._moran

    @property
    def digests(self) -> set[bytes]:
        if self._digests is None:
            self._digests = row_digests(self.X)
        return self._digests

    @property
    def keys(self) -> set[bytes]:
        if self._coord_keys is None:
            self._coord_keys = coord_keys(self.coords)
        return self._coord_keys


class Sources:
    """The source stages, loaded once per (kind, panel) and kept across the files that use them."""

    def __init__(self, data: Path):
        self.data, self.loaded = Path(data), {}

    def get(self, kind: str, rel: str, panel: Sequence[str]) -> Source:
        key = (rel, hashlib.sha256("\n".join(panel).encode()).hexdigest())
        if key not in self.loaded:
            path = self.data / rel
            if kind == "t1":
                a = load_stage(path, lineage=False)   # T1's builder: the released cells, no duplicate records
                pos = {g: i for i, g in enumerate(a.var_names)}
                X = sparse.csr_matrix(a.X, dtype=np.float32)[:, [pos[g] for g in panel]]
                coords = None
            else:
                a = load_stage(path, genes=list(panel), drop_duplicates=True)   # the T2 and T3 builders
                X = sparse.csr_matrix(a.X, dtype=np.float32)
                coords = np.asarray(a.obsm[COORD_KEY], dtype=np.float32)[:, :3]
            self.loaded[key] = Source(rel=rel, X=X, coords=coords, sha256=sha256_file(path))
        return self.loaded[key]


def load_submission(path: Path, panel: Sequence[str], coords: bool) -> tuple[object, np.ndarray | None, bool]:
    """A file's X, its coordinates when the board needs them, and whether its gene names are the panel in order."""
    a = ad.read_h5ad(path)
    same = [str(g) for g in a.var_names] == list(panel)
    return a.X, (np.asarray(a.obsm[COORD_KEY], dtype=np.float32)[:, :3] if coords else None), same


def read_sums(path: Path) -> dict[str, str]:
    """An archive's SHA256SUMS: relative path -> sha256."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if line.strip():
            sha, rel = line.split(maxsplit=1)
            out[rel.strip().lstrip("*")] = sha
    return out


# ---------------------------------------------------------------------------------------------------------------------
# One file
def _flag(flags: list[str], job: Job, what: str) -> None:
    flags.append(f"{job.file}: {what}")


def check_file(index: int, job: Job, archive: Path, sources: Sources, panel: Sequence[str], sums: Mapping[str, str],
               sizes: Mapping[str, object], decisions_sha: str, union_sha: str | None) -> tuple[dict, list[str]]:
    """Every check on one file: its result (no cell count) and its flags."""
    flags: list[str] = []
    path = Path(archive) / job.file
    file_sha = sha256_file(path)
    if sums.get(job.file) != file_sha:
        _flag(flags, job, "(e) sha256 differs from the archive's SHA256SUMS")
    panel = list(panel)
    needs_coords = job.kind != "t1"
    X, coords, same_order = load_submission(path, panel, needs_coords)
    if not same_order:
        _flag(flags, job, "(a) the file's gene names are not the board panel in order")
    srcs = [sources.get(job.kind, rel, panel) for rel in SOURCES[job.kind]]
    rng = np.random.default_rng([SAMPLE_SEED, index])
    sample = rng.integers(0, X.shape[0], K_SAMPLE)          # with replacement: K_SAMPLE rows whatever the file's size
    res: dict = {"board": job.board, "sha256": file_sha, "sources": [s.rel for s in srcs],
                 "panel_names_in_order": bool(same_order)}

    # (a) pseudobulk agreement, the embryo against the equal-weight mean of its two brackets
    pb_file = pseudobulk(X)
    pb_source = srcs[0].pb if len(srcs) == 1 else np.mean([s.pb for s in srcs], axis=0)
    res["pb"] = pb_agreement(pb_file, pb_source, PB_FLOOR[job.kind])
    if len(srcs) > 1:
        res["pb"]["per_source"] = {s.rel: {"pearson": _corr(pb_file, s.pb), "spearman": _spearman(pb_file, s.pb)}
                                   for s in srcs}
    if not res["pb"]["pass"]:
        _flag(flags, job, f"(a) pseudobulk Pearson {res['pb']['pearson']:.4f} / Spearman {res['pb']['spearman']:.4f} "
                          f"below {PB_FLOOR[job.kind]}")

    # (b) Moran's I: the file's cells on the file's own coordinates, against the source's profile
    if needs_coords:
        idx = knn_graph(coords)
        perm = np.random.default_rng([SAMPLE_SEED, index, 1]).permutation(X.shape[0])
        per = {}
        for s in srcs:
            genes, prof_src = s.moran()
            per[s.rel] = moran_check(moran_profile(X, idx, genes), prof_src, moran_profile(X[perm], idx, genes))
        best_rel = max(per, key=lambda r: (-np.inf if np.isnan(per[r]["agreement"]) else per[r]["agreement"]))
        best = per[best_rel]
        res["moran"] = {**best, "source": best_rel, "per_source": per, "floor": MORAN_FLOOR,
                        "level_floor": MORAN_LEVEL_FLOOR,
                        "pass": bool(best["agreement"] >= MORAN_FLOOR and best["level_ratio"] >= MORAN_LEVEL_FLOOR)}
        if not res["moran"]["pass"]:
            _flag(flags, job, f"(b) Moran's I agreement {best['agreement']:.4f} (floor {MORAN_FLOOR}), level ratio "
                              f"{best['level_ratio']:.4f} (floor {MORAN_LEVEL_FLOOR})")

    # (c) finite, non-negative X; T1's implied library
    res["x"] = x_sanity(X)
    if not res["x"]["pass"]:
        _flag(flags, job, f"(c) X finite {res['x']['finite']}, non-negative {res['x']['non_negative']}")
    if job.kind == "t1":
        res["library"] = implied_library(X)
        if not res["library"]["pass"]:
            _flag(flags, job, f"(c) implied library outside {list(LIBRARY_BAND)} x 1e4 (p0.1 "
                              f"{res['library']['p0.1_over_1e4']:.4f}, p99.9 {res['library']['p99.9_over_1e4']:.4f})")

    # (d) the size step, assert-only
    if job.kind in SIZES and job.kind in sizes:
        ok = size_matches(coords, job.kind, sizes[job.kind], job.hypothesis, SIZE_TOL[job.kind])
        res["size_as_preregistered"] = ok
        if not ok:
            _flag(flags, job, "(d) the size step does not match the pre-registered size JSON")

    # (e) provenance: the stamp, the build record, row-level content on a fixed-size sample
    res["stamp"] = check_stamp(path, file_sha)
    if res["stamp"]["present"] and res["stamp"]["sha256_matches"] is not True:
        _flag(flags, job, "(e) the build stamp's sha256 differs from the file")
    res["record"] = check_record(job, archive, record_expectations(job, file_sha, {s.rel: s.sha256 for s in srcs},
                                                                   decisions_sha, union_sha))
    if res["record"]["problems"]:
        _flag(flags, job, f"(e) build record fields differ from the build's: {', '.join(res['record']['problems'])}")
    if res["record"]["source_row_fields"]:
        _flag(flags, job, f"(e) a build record keeps source rows ({', '.join(res['record']['source_row_fields'])}): "
                          "the row-to-cell check needs writing before this pass can read them")
    prov: dict = {"k_sample": K_SAMPLE, "floor": EXACT_FLOOR[job.kind]}
    if job.kind == "t1":
        hit = exact_row_mask(_as_rows(X, sample), srcs[0].digests)
        prov["exact_row_fraction"] = float(hit.mean())
        prov["pass"] = bool(prov["exact_row_fraction"] >= EXACT_FLOOR["t1"])
    elif job.kind == "embryo":
        rows = _as_rows(X, sample)
        before, after = exact_row_mask(rows, srcs[0].digests), exact_row_mask(rows, srcs[1].digests)
        either = before | after
        prov["exact_row_fraction"] = float(either.mean())
        prov["share_from_first_bracket"] = float(before.sum() / either.sum()) if either.any() else None
        prov["share_band"] = list(EMBRYO_SHARE)
        share_ok = prov["share_from_first_bracket"] is not None and (
            EMBRYO_SHARE[0] <= prov["share_from_first_bracket"] <= EMBRYO_SHARE[1])
        prov["pass"] = bool(prov["exact_row_fraction"] >= EXACT_FLOOR["embryo"] and share_ok)
    elif job.kind == "t3":
        prov["exact_coord_fraction"] = float(in_keys(coords[sample], srcs[0].keys).mean())
        prov["pass"] = bool(prov["exact_coord_fraction"] >= EXACT_FLOOR["t3"])
    else:  # heart: z held by the size step (the jittered copies of repeated rows aside), so z is a z level of E10.5
        prov["z_level_fraction"] = float(np.isin(coords[sample, 2], np.unique(srcs[0].coords[:, 2])).mean())
        prov["pass"] = bool(prov["z_level_fraction"] >= EXACT_FLOOR["heart"])
    res["provenance"] = prov
    if not prov["pass"]:
        _flag(flags, job, "(e) row-level provenance " + ", ".join(f"{k} {v}" for k, v in prov.items()
                                                                  if isinstance(v, float)))
    res["flags"] = flags
    return res, flags


def _as_rows(X, rows: np.ndarray) -> np.ndarray:
    part = X[rows]
    return (part.toarray() if sparse.issparse(part) else np.asarray(part)).astype(np.float32)


# ---------------------------------------------------------------------------------------------------------------------
# The pass and its outputs
def run_pass(jobs: Sequence[Job], archive: Path, data: Path, panels: Mapping[str, Sequence[str]],
             sizes: Mapping[str, object] = SIZES, decisions: Path = DECISIONS,
             union_sha: str | None = None, sums: Mapping[str, str] | None = None) -> dict:
    """Every job's checks and the inputs' sha256 (no cell count). panels: board key -> its gene names in panel order
    (from vec.harness._candidate_board); sizes: board kind -> pre-registered size JSON; union_sha: T1's pinned R2
    union."""
    archive = Path(archive)
    sums = read_sums(archive / "SHA256SUMS") if sums is None else sums
    sources = Sources(Path(data))
    decisions_sha = sha256_file(decisions)
    files, flags = {}, []
    for i, job in enumerate(jobs):
        res, f = check_file(i, job, archive, sources, panels[job.board], sums, sizes, decisions_sha, union_sha)
        files[job.file] = res
        flags += f
    inputs = {"archive_sha256sums": sha256_file(archive / "SHA256SUMS"),
              "files_sha256": {j.file: files[j.file]["sha256"] for j in jobs},
              "sources_sha256": {s.rel: s.sha256 for s in sources.loaded.values()},
              "panels_sha256": {b: hashlib.sha256("\n".join(p).encode()).hexdigest() for b, p in panels.items()},
              "size_json_sha256": {k: sha256_file(v) for k, v in sizes.items()},
              "decisions_sha256": decisions_sha}
    return {"label": LABEL, "verdict": "flag" if flags else "no flag", "flags": flags, "rules": RULES, "files": files,
            "source_rows_in_build_records": any(r["record"]["source_row_fields"] for r in files.values()),
            "inputs": inputs}


def _fmt(x: float | None, nd: int = 4) -> str:
    return "n/a" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def report_text(res: dict, stamp: Mapping, archive: str) -> str:
    """The pass as text: the verdict, the flags and one line per file (no cell count, no size, no anchor)."""
    lines = [f"{LABEL} (commit {stamp['commit']}, tree clean {stamp['tree_clean']}, UTC {stamp['utc']})",
             f"Archive: {archive}",
             f"VERDICT: {res['verdict']}" + ("" if not res["flags"] else " -- stop, trace, and take it to the captain; "
                                             "no file is rebuilt or uploaded until the cause is written down"),
             "Build records keeping per-cell source rows: "
             + ("some (the sharpener's row-to-cell check applies)" if res["source_rows_in_build_records"]
                else "none (the row-to-cell check has nothing to read)")]
    lines += [f"FLAG {f}" for f in res["flags"]]
    for rel, r in res["files"].items():
        parts = [f"{rel} ({r['board']}; sources {', '.join(r['sources'])})",
                 f"(a) pb r {_fmt(r['pb']['pearson'])}, rho {_fmt(r['pb']['spearman'])}, floor {r['pb']['floor']}"]
        if "moran" in r:
            m = r["moran"]
            parts.append(f"(b) Moran r {_fmt(m['agreement'])} vs {m['source']} (floor {m['floor']}), level ratio "
                         f"{_fmt(m['level_ratio'])} (floor {m['level_floor']}), shuffled control r "
                         f"{_fmt(m['shuffled_control_agreement'])}")
        parts.append(f"(c) finite {r['x']['finite']}, non-negative {r['x']['non_negative']}")
        if "library" in r:
            lb = r["library"]
            parts.append(f"implied library / 1e4: median {_fmt(lb['median_over_1e4'], 5)}, p0.1 "
                         f"{_fmt(lb['p0.1_over_1e4'], 5)}, p99.9 {_fmt(lb['p99.9_over_1e4'], 5)}")
        if "size_as_preregistered" in r:
            parts.append(f"(d) size as pre-registered: {r['size_as_preregistered']}")
        p, rc = r["provenance"], r["record"]
        provs = ", ".join(f"{k} {_fmt(v)}" for k, v in p.items() if isinstance(v, float) and k != "floor")
        parts.append(f"(e) {provs} (floor {p['floor']}); stamp sha256 ok {r['stamp']['sha256_matches']}, clean "
                     f"{r['stamp']['clean']}; build record present {rc['present']}, fields differing "
                     f"{rc['problems'] or 'none'}, missing {rc['missing'] or 'none'}")
        lines.append("  " + "; ".join(parts))
    lines.append(f"Inputs (sha256): {json.dumps(res['inputs'], sort_keys=True)}")
    return "\n".join(lines) + "\n"


def tree_state(exempt: Sequence[Path]) -> dict:
    """HEAD and whether the tree is clean by vec.provenance's rule, the outputs exempt (None when git cannot tell)."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    return {"commit": commit, "clean": None if unclean is None else not unclean}


def main(argv=None, jobs: Sequence[Job] = JOBS, panels: Mapping[str, Sequence[str]] | None = None,
         sizes: Mapping[str, object] = SIZES, decisions: Path = DECISIONS) -> None:
    ap = argparse.ArgumentParser(description="N14: the integrity pass of the ten frozen upload files")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--archive", default=str(ARCHIVE))
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    out = Path(a.out)
    names = {"json": out / "integrity.json", "text": out / "integrity.txt"}
    exempt = [*names.values(), *(provenance.stamp_path(p) for p in names.values())]
    if out.exists() and any(out.iterdir()):
        sys.exit(f"N14 has already written its outputs in {out}; a rerun needs a new experiments/harness/TRIES.md row")
    started = tree_state(exempt)
    if started["clean"] is not True:
        sys.exit("N14 runs once on a clean tree (vec.provenance's rule): commit or remove what keeps it from clean")
    if panels is None:
        panels = {j.board: list(harness._candidate_board(load_boards()[j.board]).genes) for j in jobs}
    union_sha = hygiene.load_config().get("union_sha256")
    res = run_pass(jobs, Path(a.archive), Path(a.data), panels, sizes=sizes, decisions=decisions, union_sha=union_sha)
    ended = tree_state(exempt)
    clean = (None if started["clean"] is None or ended["clean"] is None
             else started["clean"] and ended["clean"] and started["commit"] == ended["commit"])
    stamp = {"commit": ended["commit"], "tree_clean": clean, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out.mkdir(parents=True, exist_ok=True)
    body = {**res, "provenance": {**stamp, "archive": _relative(a.archive), "data": _relative(a.data)}}
    names["json"].write_text(json.dumps(body, indent=1, sort_keys=True) + "\n")
    text = report_text(res, stamp, _relative(a.archive))
    names["text"].write_text(text)
    for p in names.values():
        provenance.write_stamp(p)
    print(text, end="")


def _relative(path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO).as_posix()
    except ValueError:
        return str(path)


if __name__ == "__main__":
    main()
