"""I/O for the external-data window filters: row filters that never write a dropped row, whole-file gates, manifests.

Every filter here decides from metadata alone (a stage or embryo column, or a file's embryo number) and writes only
kept rows or kept files. Row filters stream the input and write kept rows to a temporary file in the destination
directory, renamed into place only when the whole input passed (a failed run leaves no partial output). Whole files
that pass a gate are hard-linked into the destination when it shares a filesystem with the source (no extra disk; the
copy survives deletion of the source), else copied.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from collections.abc import Callable, Sequence
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
# Field size: some metadata tables carry long free-text columns.
csv.field_size_limit(sys.maxsize)


def sha256_file(path: str | Path, chunk: int = 1 << 23) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def _git(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=REPO, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def provenance(code_files: Sequence[str | Path]) -> dict:
    """git commit, uncommitted changes to the filter code, and the code's sha256 (recorded in every manifest)."""
    rel = [str(Path(p).resolve().relative_to(REPO)) for p in code_files]
    return {"git_commit": _git("rev-parse", "HEAD"),
            "git_uncommitted": [ln for ln in _git("status", "--porcelain", "--", *rel).splitlines() if ln],
            "code_sha256": {p: sha256_file(REPO / p) for p in rel},
            "run_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")}


def _tmp(dst: Path) -> Path:
    return dst.with_name(f".{dst.name}.partial")


def filter_csv(src: str | Path, dst: str | Path, key_columns: Sequence[str],
               decide: Callable[[tuple[str, ...]], tuple[bool, str]], keep_positions: bool = False) -> dict:
    """Stream `src` row by row and write only the rows that `decide` keeps.

    `decide(values)` gets the row's values in `key_columns` and returns (keep, reason); it is evaluated once per
    distinct key and cached. Dropped rows are never written: kept rows go to a temporary file that is renamed to `dst`
    only after the last row passed. Malformed rows (wrong field count) abort the run. Returns counts per reason and
    per key, and with `keep_positions` the 0-based data-row positions kept (to subset aligned matrices).
    """
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp(dst)
    cache: dict[tuple[str, ...], tuple[bool, str]] = {}
    by_reason: Counter = Counter()
    by_key: Counter = Counter()
    kept_pos: list[int] = []
    n_in = n_kept = 0
    try:
        with open(src, newline="") as fin, open(tmp, "w", newline="") as fout:
            reader, writer = csv.reader(fin), csv.writer(fout, lineterminator="\n")
            header = next(reader)
            missing = [c for c in key_columns if c not in header]
            if missing:
                raise KeyError(f"{src}: key columns {missing} not in header")
            idx = [header.index(c) for c in key_columns]
            writer.writerow(header)
            for row in reader:
                if len(row) != len(header):
                    raise ValueError(f"{src}: data row {n_in} has {len(row)} fields, header has {len(header)}")
                key = tuple(row[i] for i in idx)
                if key not in cache:
                    cache[key] = decide(key)
                keep, reason = cache[key]
                by_reason[reason] += 1
                by_key[(key, keep)] += 1
                if keep:
                    writer.writerow(row)
                    if keep_positions:
                        kept_pos.append(n_in)
                    n_kept += 1
                n_in += 1
        os.replace(tmp, dst)
    finally:
        if tmp.exists():
            tmp.unlink()
    out = {"rows_in": n_in, "rows_kept": n_kept, "rows_dropped": n_in - n_kept, "by_reason": dict(by_reason),
           "by_key": [{"key": list(k), "kept": kept, "rows": n} for (k, kept), n in sorted(by_key.items())]}
    if keep_positions:
        out["positions"] = np.asarray(kept_pos, dtype=np.int64)
    return out


def link_or_copy(src: str | Path, dst: str | Path, overwrite: bool = False) -> str:
    """Place a file that passed a whole-file gate: hard link if possible, else copy. Returns the mode used."""
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if os.path.samefile(src, dst):
            return "exists (same file)"
        if not overwrite:
            raise FileExistsError(f"{dst} exists and is a different file (pass --overwrite to replace it)")
        dst.unlink()
    try:
        os.link(src, dst)
        return "hardlink"
    except OSError:
        tmp = _tmp(dst)
        shutil.copy2(src, tmp)
        os.replace(tmp, dst)
        return "copy"


def save_npz_atomic(dst: str | Path, **arrays) -> None:
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.stem}.partial.npz")
    try:
        np.savez_compressed(tmp, **arrays)
        os.replace(tmp, dst)
    finally:
        if tmp.exists():
            tmp.unlink()


def write_manifest(dst_dir: str | Path, manifest: dict) -> Path:
    path = Path(dst_dir) / "MANIFEST.json"
    tmp = _tmp(path)
    tmp.write_text(json.dumps(manifest, indent=1, default=_json_default) + "\n")
    os.replace(tmp, path)
    return path


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def verify_manifest_output(path: str | Path) -> dict:
    """Check that a file under data/external was written by a committed filter: the MANIFEST.json beside it (or in a
    parent directory) lists it with a matching sha256. Returns that manifest entry; raises otherwise."""
    path = Path(path).resolve()
    for d in path.parents:
        m = d / "MANIFEST.json"
        if m.exists():
            man = json.loads(m.read_text())
            rel = str(path.relative_to(d))
            entry = man.get("outputs", {}).get(rel)
            if entry is None:
                raise ValueError(f"{path} is not listed in {m}")
            if entry.get("sha256") and entry["sha256"] != sha256_file(path):
                raise ValueError(f"{path} does not match the sha256 recorded in {m}")
            return entry
    raise FileNotFoundError(f"no MANIFEST.json above {path}: external data must come through experiments/external")
