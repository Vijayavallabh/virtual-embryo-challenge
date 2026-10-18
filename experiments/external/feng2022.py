"""Feng et al. 2022 (GEO GSE193346), CD1 mouse heart atlas: window filter on the per-cell stage column.

Keep rule (docs/DATA_POLICY.md; docs/PLAN.md section 4.1): stage E9.5, or stage E14.5 or later; strain CD1. Every
other stage (E10.5-E13.5 in this atlas, and anything else strictly between E9.5 and E14.5) is dropped before anything
is written. A stage label that does not parse aborts the run (fail closed). Kept rows are re-checked against the
protected windows of every task named with --tasks (default T1), and any disagreement aborts.

Inputs: the legal subsets already on the DGX (paths in sources.yaml). The 2026-10-18 research run produced them by
filtering the GEO object on this same stage column; this script re-applies the rule as committed code:
  cd1_E9.5_E14.5_obs.csv     per-cell metadata (stage, strain, zone, batch, cell_type_spec) of the E9.5 and E14.5 cells
  cd1_E9.5_E14.5_counts.npz  raw UMI counts of the same cells (scipy sparse), rows aligned with the obs table
  cd1_stats.npz              per-group sufficient statistics; group keys are 'stage|zone|cell_type_spec|half', built
                             from the per-cell stage column, so a key's stage token is that column's value
Not taken: legal/counts_legal.mtx (all 14 legal stages). Its per-cell stage labels exist only in the GEO metadata
table, deleted from the DGX on 2026-10-19, and its barcodes carry no stage, so the rule cannot be re-applied to it.
To bring those cells in later, register the source, fetch the metadata table and run `filter_obs` on it (it streams
the table and never writes an in-window row), then subset the counts by the kept barcodes.

Usage (DGX, repo root; seconds):
  uv run python experiments/external/feng2022.py [--src DIR] [--dst data/external/feng2022_gse193346] [--tasks T1]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import filterio as fio  # noqa: E402
import windows  # noqa: E402

SOURCE = "feng2022"
STRAINS = ("CD1",)
START_EDAY, LATE_MIN_EDAY = 9.5, 14.5
CODE = [Path(__file__), Path(windows.__file__), Path(fio.__file__)]


def keep_stage(stage: str, tasks=("T1",)) -> tuple[bool, str]:
    """The stage half of the keep rule: E9.5 exactly, or E14.5 and later. Unknown labels raise."""
    e = windows.parse_stage(stage)
    if not (e == START_EDAY or e >= LATE_MIN_EDAY):
        return False, f"dropped: stage {stage} (not E9.5 and before E14.5)"
    if not windows.allowed(e, tasks):  # the rule above already excludes every window; this guards edits to it
        raise AssertionError(f"keep rule would keep {stage}, which lies in a protected window of {tasks}")
    return True, "kept"


def decide(values: tuple[str, ...], tasks=("T1",)) -> tuple[bool, str]:
    """Keep rule for one row, from its (stage, strain) values."""
    stage, strain = values
    keep, reason = keep_stage(stage, tasks)
    if keep and strain not in STRAINS:
        return False, f"dropped: strain {strain}"
    return keep, reason


def filter_obs(src: str | Path, dst: str | Path, tasks=("T1",)) -> dict:
    """Per-cell filter on the stage and strain columns of a metadata table (legal subset or the full GEO table)."""
    return fio.filter_csv(src, dst, ["stage", "strain"], lambda v: decide(v, tasks), keep_positions=True)


def filter_counts(src: str | Path, dst: str | Path, positions: np.ndarray, n_rows: int) -> dict:
    """Subset a cells x genes count matrix to the rows the obs filter kept (rows must align with the obs table)."""
    X = sparse.load_npz(src).tocsr()
    if X.shape[0] != n_rows:
        raise ValueError(f"{src}: {X.shape[0]} rows but the obs table has {n_rows}; rows must align")
    Y = X[positions]
    dst = Path(dst)
    tmp = dst.with_name(f".{dst.stem}.partial.npz")
    try:
        sparse.save_npz(tmp, Y, compressed=True)
        tmp.replace(dst)
    finally:
        if tmp.exists():
            tmp.unlink()
    return {"rows_in": int(X.shape[0]), "rows_kept": int(Y.shape[0]), "genes": int(Y.shape[1])}


def filter_stats(src: str | Path, dst: str | Path, tasks=("T1",)) -> dict:
    """Group-level filter of the sufficient statistics: keep a group if its stage token passes the stage rule."""
    z = np.load(src, allow_pickle=True)  # written by the 2026-10-18 research run (object arrays)
    keys = np.array([str(k) for k in z["keys"]])
    stage = np.array([k.split("|")[0] for k in keys])
    verdict = {s: keep_stage(s, tasks) for s in np.unique(stage)}
    keep = np.array([verdict[s][0] for s in stage])
    n = np.asarray(z["n"], dtype=np.float64)
    out = {"keys": keys[keep].astype(str), "genes": np.array([str(g) for g in z["genes"]]).astype(str)}
    for k in z.files:
        if k in ("keys", "genes"):
            continue
        arr = np.asarray(z[k])
        if arr.ndim == 0 or arr.shape[0] != len(keys):
            raise ValueError(f"{src}: field {k!r} is not per-group (shape {arr.shape}); refusing to guess")
        out[k] = arr[keep]
    fio.save_npz_atomic(dst, **out)
    cells = Counter()
    for s, ni, kp in zip(stage, n, keep, strict=True):
        cells[(s, bool(kp))] += ni
    return {"groups_in": int(len(keys)), "groups_kept": int(keep.sum()),
            "cells_by_stage": [{"stage": s, "kept": kp, "cells": int(c)} for (s, kp), c in sorted(cells.items())],
            "cells_kept": int(n[keep].sum()), "cells_dropped": int(n[~keep].sum())}


def run(src: Path, dst: Path, inputs: dict, tasks=("T1",)) -> dict:
    dst.mkdir(parents=True, exist_ok=True)
    paths = {k: src / v for k, v in inputs.items()}
    for k, p in paths.items():
        if not p.exists():
            raise FileNotFoundError(f"{k}: {p}")
    outputs = {"obs": "cd1_E9.5_E14.5_obs.csv", "counts": "cd1_E9.5_E14.5_counts.npz", "stats": "cd1_stats.npz"}
    rep = {"obs": filter_obs(paths["obs"], dst / outputs["obs"], tasks)}
    pos = rep["obs"].pop("positions")
    rep["counts"] = filter_counts(paths["counts"], dst / outputs["counts"], pos, rep["obs"]["rows_in"])
    rep["stats"] = filter_stats(paths["stats"], dst / outputs["stats"], tasks)
    return {"inputs": {k: {"path": str(p), "sha256": fio.sha256_file(p)} for k, p in paths.items()},
            "outputs": {v: {"from": k, "sha256": fio.sha256_file(dst / v), **rep[k]} for k, v in outputs.items()}}


def main(argv=None) -> int:
    cfg = yaml.safe_load((fio.HERE / "sources.yaml").read_text())[SOURCE]
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src", default=cfg["src_dir"], help="directory holding the legal subsets")
    ap.add_argument("--dst", default=cfg["dst"], help="output directory (must be under data/, which is gitignored)")
    ap.add_argument("--tasks", nargs="+", default=list(cfg["tasks"]), help="tasks the data may serve (union)")
    a = ap.parse_args(argv)
    tasks = windows.check_tasks(a.tasks)
    src = Path(a.src).expanduser()
    dst = (fio.REPO / a.dst).resolve() if not Path(a.dst).is_absolute() else Path(a.dst)
    if fio.REPO / "data" not in dst.parents:
        raise SystemExit(f"--dst must be inside {fio.REPO / 'data'} (gitignored); got {dst}")
    res = run(src, dst, cfg["inputs"], tasks)
    manifest = {"source": SOURCE, "citation": cfg["citation"], "accession": cfg["accession"], "tasks": list(tasks),
                "keep_rule": cfg["keep"], "not_taken": cfg.get("not_taken", []), **res,
                "provenance": fio.provenance(CODE)}
    path = fio.write_manifest(dst, manifest)
    for name, o in res["outputs"].items():
        print(name, {k: v for k, v in o.items() if k not in ("by_key", "sha256", "cells_by_stage")})
    print("stages kept (cells):", {r["stage"]: r["cells"] for r in res["outputs"]["cd1_stats.npz"]["cells_by_stage"]
                                   if r["kept"]})
    print("manifest:", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
