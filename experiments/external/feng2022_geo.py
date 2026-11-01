"""Feng et al. 2022 (GEO GSE193346), CD1 atlas: the clean reduction from the raw GEO files (feat-007).

The legal subsets in data/external/feng2022_gse193346/ came from screened research code (2026-10-18). This script
redoes the reduction from GEO as committed code, as EXTERNAL_DATA.md (X002) requires before any submission uses Feng.
The captain's go-ahead (2026-11-01): the download contains E10.5-E13.5 cells; the filter drops them before anything is
written, and the raw files are deleted afterwards.

Steps (DGX, repo root; `uv run --group external python experiments/external/feng2022_geo.py`; ~30 min, mostly the
download and the parse):
1. Download GEO's CD1 metadata table and CD1 Seurat object into a private folder under data/ (mode 700), recording
   their sha256.
2. Stream the metadata table and write only the rows feng2022.py's keep rule keeps (stage E9.5, or E14.5 and later, on
   the per-cell stage column; strain CD1). Dropped rows are never written and are reported as one total, not by stage.
3. Parse the Seurat object in memory with the pure-Python `rdata` reader (no R) and take the raw counts of the kept
   cells only (matched by cell name).
4. Write the outputs (below) and a MANIFEST.json, compare them with the 2026-10-18 subsets, and delete the raw files.

Outputs (data/external/feng2022_gse193346_geo/):
  cd1_legal_obs.csv        GEO metadata rows of the kept cells (every legal stage), columns unchanged
  cd1_legal_counts.npz     their raw UMI counts, cells x genes (scipy CSR int32), rows aligned with the obs table
  cd1_legal_genes.txt      gene symbols, the count matrix's columns
  cd1_stats.npz            per-group sufficient statistics, keys 'stage|zone|cell_type_spec|half', fields n, S_log
                           (sum of log1p CP10k), S_cp (sum of CP10k), S_nz (cells detected), depth_sum, genes, keys;
                           the halves are a seeded random split within each group (both halves summed are what the T1
                           direction code uses)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng2022  # noqa: E402
import filterio as fio  # noqa: E402
import windows  # noqa: E402

SOURCE = "feng2022_geo"
CODE = [Path(__file__), Path(feng2022.__file__), Path(windows.__file__), Path(fio.__file__)]
HALF_SEED = 20261101
CP10K = 1e4


# ------------------------------------------------------------------------------------------------ download


def download(url: str, dst: Path, retries: int = 3) -> dict:
    """Stream `url` to `dst` (via a partial file), returning its size and sha256."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".partial")
    for attempt in range(1, retries + 1):
        h, size = hashlib.sha256(), 0
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
                while block := r.read(1 << 23):
                    f.write(block)
                    h.update(block)
                    size += len(block)
            tmp.replace(dst)
            return {"url": url, "bytes": size, "sha256": h.hexdigest()}
        except OSError as e:
            if tmp.exists():
                tmp.unlink()
            if attempt == retries:
                raise
            print(f"download attempt {attempt} failed ({e}); retrying", flush=True)
            time.sleep(10 * attempt)
    raise AssertionError("unreachable")


# ------------------------------------------------------------------------------------------------ metadata


def filter_geo_obs(src_gz: Path, dst: Path, tasks=("T1",)) -> dict:
    """Stream the gzipped GEO metadata table and write only kept rows (feng2022.decide on stage and strain). The first
    column holds the cell names. Returns the kept cell names in row order, kept cells per stage, and one total of
    dropped rows (nothing about the dropped rows beyond that number)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.name}.partial")
    kept_names: list[str] = []
    kept_by_stage: Counter = Counter()
    n_in = n_dropped = 0
    cache: dict[tuple[str, str], bool] = {}
    try:
        with gzip.open(src_gz, "rt", newline="") as fin, open(tmp, "w", newline="") as fout:
            reader, writer = csv.reader(fin), csv.writer(fout, lineterminator="\n")
            header = next(reader)
            missing = [c for c in ("stage", "strain") if c not in header]
            if missing:
                raise KeyError(f"{src_gz.name}: columns {missing} not in the header {header[:12]}...")
            i_stage, i_strain = header.index("stage"), header.index("strain")
            writer.writerow(header)
            for row in reader:
                if len(row) != len(header):
                    raise ValueError(f"{src_gz.name}: data row {n_in} has {len(row)} fields, header {len(header)}")
                key = (row[i_stage], row[i_strain])
                if key not in cache:
                    cache[key] = feng2022.decide(key, tasks)[0]
                n_in += 1
                if cache[key]:
                    writer.writerow(row)
                    kept_names.append(row[0])
                    kept_by_stage[row[i_stage]] += 1
                else:
                    n_dropped += 1
        tmp.replace(dst)
    finally:
        if tmp.exists():
            tmp.unlink()
    if len(set(kept_names)) != len(kept_names):
        raise ValueError(f"{src_gz.name}: duplicate cell names among the kept rows")
    return {"rows_in": n_in, "rows_kept": len(kept_names), "rows_dropped": n_dropped,
            "kept_by_stage": dict(sorted(kept_by_stage.items(), key=lambda kv: windows.parse_stage(kv[0]))),
            "names": kept_names}


# ------------------------------------------------------------------------------------------------ R objects (rdata)


def _types():
    from rdata.parser import RObjectType

    return RObjectType


def deref(o):
    """Follow rdata's reference objects to the object they stand for."""
    T = _types()
    while o is not None and o.info.type == T.REF:
        o = o.referenced_object
    return o


def _name(tag) -> str:
    T = _types()
    t = deref(tag)
    if t.info.type == T.SYM:
        t = deref(t.value)
    if t.info.type != T.CHAR or not isinstance(t.value, bytes):
        raise ValueError(f"unexpected tag of type {t.info.type}")
    return t.value.decode("utf-8")


def pairlist(o):
    """(name, value) pairs of an R pairlist (attributes and S4 slots are pairlists)."""
    T = _types()
    o = deref(o)
    while o is not None and o.info.type in (T.LIST, T.ATTRLIST):
        car, cdr = o.value
        yield (_name(o.tag) if o.tag is not None else None), deref(car)
        o = deref(cdr)


def attributes(o) -> dict:
    o = deref(o)
    return dict(pairlist(o.attributes)) if o.attributes is not None else {}


def strings(o) -> list[str | None]:
    """An R character vector as Python strings (NA as None)."""
    o = deref(o)
    out = []
    for c in o.value:
        v = deref(c).value
        out.append(None if v is None else v.decode("utf-8", "replace"))
    return out


def named_list(o) -> dict:
    o = deref(o)
    names = strings(attributes(o)["names"])
    return dict(zip(names, (deref(v) for v in o.value), strict=True))


def slot(o, name: str):
    a = attributes(o)
    if name not in a:
        raise KeyError(f"no slot {name!r}; slots: {sorted(k for k in a if k)}")
    return a[name]


def r_class(o) -> str:
    return strings(attributes(o)["class"])[0]


def dgc_matrix(m) -> tuple[sparse.csc_matrix, list[str] | None, list[str] | None]:
    """An R dgCMatrix (genes x cells) as scipy CSC, with its row and column names."""
    if r_class(m) != "dgCMatrix":
        raise ValueError(f"expected a dgCMatrix, got {r_class(m)}")
    dim = np.asarray(slot(m, "Dim").value, dtype=np.int64)
    X = sparse.csc_matrix((np.asarray(slot(m, "x").value), np.asarray(slot(m, "i").value, dtype=np.int64),
                           np.asarray(slot(m, "p").value, dtype=np.int64)), shape=(int(dim[0]), int(dim[1])))
    dn = slot(m, "Dimnames")
    rows = cols = None
    if deref(dn).value:
        r, c = (deref(v) for v in deref(dn).value)
        rows = strings(r) if r.value is not None else None
        cols = strings(c) if c.value is not None else None
    return X, rows, cols


def seurat_counts(obj) -> tuple[sparse.csc_matrix, list[str], list[str], dict]:
    """The RNA assay's raw counts of a Seurat object (Assay v3/v4 `counts` slot; Assay5 `layers$counts`)."""
    if r_class(obj) != "Seurat":
        raise ValueError(f"expected a Seurat object, got {r_class(obj)}")
    assays = named_list(slot(obj, "assays"))
    if "RNA" not in assays:
        raise KeyError(f"no RNA assay; assays: {sorted(assays)}")
    rna = assays["RNA"]
    kind = r_class(rna)
    if kind == "Assay5":
        X, genes, cells = dgc_matrix(named_list(slot(rna, "layers"))["counts"])
    else:
        X, genes, cells = dgc_matrix(slot(rna, "counts"))
    if genes is None or cells is None:
        raise ValueError("the counts matrix has no row or column names")
    return X, genes, cells, {"assay_class": kind, "genes": X.shape[0], "cells": X.shape[1], "nnz": int(X.nnz)}


def take_cells(X: sparse.csc_matrix, cells: list[str], names: list[str]) -> sparse.csr_matrix:
    """The columns of `X` (genes x cells) for the cell `names`, in that order, as integer counts, cells x genes."""
    pos = {c: j for j, c in enumerate(cells)}
    missing = [n for n in names if n not in pos]
    if missing:
        raise KeyError(f"{len(missing)} kept cells are not in the count matrix, e.g. {missing[:3]}")
    Y = X[:, [pos[n] for n in names]].T.tocsr()
    if Y.nnz and not np.allclose(Y.data, np.round(Y.data)):
        raise ValueError("the counts are not integers")
    Y.data = np.round(Y.data).astype(np.int32)
    return Y


# ------------------------------------------------------------------------------------------------ statistics


def group_stats(C: sparse.csr_matrix, stage, zone, cell_type, names, genes, seed: int = HALF_SEED) -> dict:
    """Per-group sufficient statistics over cells (rows of integer counts C). Groups are stage|zone|cell_type; within
    each, the cells sorted by name are split at random into halves A (the larger) and B, seeded per group."""
    stage, zone, cell_type, names = (np.asarray(v).astype(str) for v in (stage, zone, cell_type, names))
    base = np.char.add(np.char.add(np.char.add(np.char.add(stage, "|"), zone), "|"), cell_type)
    lib = np.asarray(C.sum(1)).ravel().astype(np.float64)
    if (lib <= 0).any():
        raise ValueError("a kept cell has no counts")
    cp = sparse.diags(CP10K / lib) @ C.astype(np.float64)
    lg = cp.copy()
    lg.data = np.log1p(lg.data)
    nz = C.copy()
    nz.data = np.ones_like(nz.data)
    keys, n, S_log, S_cp, S_nz, depth = [], [], [], [], [], []
    for k, g in enumerate(sorted(set(base))):
        rows = np.flatnonzero(base == g)
        rows = rows[np.argsort(names[rows], kind="stable")]
        perm = np.random.default_rng([seed, k]).permutation(len(rows))
        half_a = len(rows) - len(rows) // 2
        for h, sel in (("A", rows[perm[:half_a]]), ("B", rows[perm[half_a:]])):
            if len(sel) == 0:
                continue
            keys.append(f"{g}|{h}")
            n.append(len(sel))
            S_log.append(np.asarray(lg[sel].sum(0)).ravel())
            S_cp.append(np.asarray(cp[sel].sum(0)).ravel())
            S_nz.append(np.asarray(nz[sel].sum(0)).ravel())
            depth.append(lib[sel].sum())
    return {"keys": np.array(keys), "genes": np.asarray(genes).astype(str), "n": np.asarray(n, dtype=np.float64),
            "S_log": np.vstack(S_log), "S_cp": np.vstack(S_cp), "S_nz": np.vstack(S_nz),
            "depth_sum": np.asarray(depth, dtype=np.float64)}


def totals_by_group(st: dict) -> dict[str, dict[str, np.ndarray]]:
    """Both halves summed per stage|zone|cell_type group."""
    out: dict[str, dict] = {}
    for i, k in enumerate(np.asarray(st["keys"]).astype(str)):
        g = k.rsplit("|", 1)[0]
        d = out.setdefault(g, {"n": 0.0, "S_log": 0.0, "S_cp": 0.0, "S_nz": 0.0, "depth_sum": 0.0})
        for f in d:
            d[f] = d[f] + np.asarray(st[f][i], dtype=np.float64)
    return out


def compare_with_old(new: dict, new_genes: list[str], new_obs_names: list[str], new_counts: sparse.csr_matrix,
                     new_stage: np.ndarray, old_dir: Path) -> dict:
    """The new outputs against the 2026-10-18 subsets: same E9.5 and E14.5 cells and counts; same group totals."""
    rep: dict = {}
    old_stats = np.load(old_dir / "cd1_stats.npz", allow_pickle=True)
    old = {k: old_stats[k] for k in old_stats.files}
    old_genes = [str(g) for g in old["genes"]]
    rep["genes_same_order"] = old_genes == list(new_genes)
    gpos = {g: j for j, g in enumerate(new_genes)}
    common = [g for g in old_genes if g in gpos]
    rep["genes_old_new_common"] = [len(old_genes), len(new_genes), len(common)]
    if rep["genes_same_order"]:
        oi = ni = np.arange(len(old_genes))
    else:
        opos = {g: j for j, g in enumerate(old_genes)}
        oi, ni = np.array([opos[g] for g in common]), np.array([gpos[g] for g in common])
    to, tn = totals_by_group(old), totals_by_group(new)
    shared = sorted(set(to) & set(tn))
    rep["groups_old_new_shared"] = [len(to), len(tn), len(shared)]
    worst = {f: 0.0 for f in ("n", "S_log", "S_cp", "S_nz", "depth_sum")}
    for g in shared:
        for f in worst:
            a, b = np.atleast_1d(to[g][f]), np.atleast_1d(tn[g][f])
            if a.ndim and a.size > 1:
                a, b = a[oi], b[ni]
            scale = np.maximum(np.abs(a), 1.0)
            worst[f] = max(worst[f], float(np.max(np.abs(a - b) / scale)))
    rep["group_totals_max_rel_diff"] = worst
    obs_old = (old_dir / "cd1_E9.5_E14.5_obs.csv").read_text().splitlines()
    old_names = [ln.split(",", 1)[0] for ln in obs_old[1:]]
    sel = np.flatnonzero(np.isin(new_stage, ["E9.5", "E14.5"]))
    names_new = [new_obs_names[i] for i in sel]
    rep["e9.5_e14.5_cells_old_new"] = [len(old_names), len(names_new), len(set(old_names) & set(names_new))]
    if set(old_names) == set(names_new) and rep["genes_same_order"]:
        Xo = sparse.load_npz(old_dir / "cd1_E9.5_E14.5_counts.npz").tocsr()
        pos = {n: i for i, n in enumerate(names_new)}
        Xn = new_counts[sel][[pos[n] for n in old_names]]
        rep["e9.5_e14.5_counts_identical"] = bool(Xo.shape == Xn.shape and (Xo != Xn).nnz == 0)
    else:
        rep["e9.5_e14.5_counts_identical"] = None
    return rep


# ------------------------------------------------------------------------------------------------ run


def run(cfg: dict, tasks, raw_dir: Path, dst: Path, old_dir: Path | None, keep_raw: bool) -> dict:
    raw_dir.mkdir(parents=True, exist_ok=True)
    os.chmod(raw_dir, 0o700)
    dst.mkdir(parents=True, exist_ok=True)
    raw = {k: raw_dir / Path(u).name for k, u in cfg["urls"].items()}
    t0 = time.time()
    got = {}
    for k, u in cfg["urls"].items():
        got[k] = download(u, raw[k])
        print(f"downloaded {raw[k].name}: {got[k]['bytes'] / 1e9:.2f} GB ({time.time() - t0:.0f} s)", flush=True)
    obs = filter_geo_obs(raw["metadata"], dst / "cd1_legal_obs.csv", tasks)
    names = obs.pop("names")
    print(f"metadata: kept {obs['rows_kept']} of {obs['rows_in']} cells; by stage {obs['kept_by_stage']}", flush=True)
    from rdata.parser import parse_file

    t1 = time.time()
    obj = deref(parse_file(raw["seurat"]).object)
    X, genes, cells, info = seurat_counts(obj)
    del obj
    print(f"parsed the Seurat object in {time.time() - t1:.0f} s: {info}", flush=True)
    C = take_cells(X, cells, names)
    del X
    with open(dst / "cd1_legal_obs.csv", newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        cols = {c: header.index(c) for c in ("stage", "zone", "cell_type_spec")}
        rows = list(rd)
    if [r[0] for r in rows] != names:
        raise AssertionError("the written obs table and the kept cell names disagree")
    stage = np.array([r[cols["stage"]] for r in rows])
    st = group_stats(C, stage, [r[cols["zone"]] for r in rows], [r[cols["cell_type_spec"]] for r in rows], names,
                     genes)
    tmp = dst / ".cd1_legal_counts.partial.npz"
    sparse.save_npz(tmp, C, compressed=True)
    tmp.replace(dst / "cd1_legal_counts.npz")
    (dst / "cd1_legal_genes.txt").write_text("\n".join(genes) + "\n")
    fio.save_npz_atomic(dst / "cd1_stats.npz", **st)
    compared = compare_with_old(st, genes, names, C, stage, old_dir) if old_dir and old_dir.exists() else None
    outputs = {n: {"sha256": fio.sha256_file(dst / n)} for n in
               ("cd1_legal_obs.csv", "cd1_legal_counts.npz", "cd1_legal_genes.txt", "cd1_stats.npz")}
    outputs["cd1_legal_obs.csv"].update(obs)
    outputs["cd1_legal_counts.npz"].update({"cells": int(C.shape[0]), "genes": int(C.shape[1])})
    outputs["cd1_stats.npz"].update({"groups": int(len(st["keys"])), "half_seed": HALF_SEED})
    manifest = {"source": SOURCE, "citation": cfg["citation"], "accession": cfg["accession"], "tasks": list(tasks),
                "keep_rule": cfg["keep"], "raw_inputs": got, "seurat": info, "outputs": outputs,
                "compared_with_2026_10_18_subsets": compared, "raw_deleted": not keep_raw,
                "provenance": fio.provenance(CODE)}
    if not keep_raw:
        shutil.rmtree(raw_dir)
    path = fio.write_manifest(dst, manifest)
    print(f"compared with the 2026-10-18 subsets: {compared}", flush=True)
    print(f"manifest: {path}; raw files {'kept in ' + str(raw_dir) if keep_raw else 'deleted'}", flush=True)
    return manifest


def main(argv=None) -> int:
    cfg = yaml.safe_load((fio.HERE / "sources.yaml").read_text())[SOURCE]
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raw", default=cfg["raw_dir"], help="private folder for the raw GEO files (under data/)")
    ap.add_argument("--dst", default=cfg["dst"], help="output folder (under data/, which is gitignored)")
    ap.add_argument("--old", default=cfg["compare_with"], help="the 2026-10-18 subsets to compare with")
    ap.add_argument("--keep-raw", action="store_true", help="keep the raw files (default: delete them at the end)")
    a = ap.parse_args(argv)
    tasks = windows.check_tasks(cfg["tasks"])
    paths = {k: (fio.REPO / v).resolve() for k, v in (("raw", a.raw), ("dst", a.dst), ("old", a.old))}
    for k in ("raw", "dst"):
        if fio.REPO / "data" not in paths[k].parents:
            raise SystemExit(f"--{k} must be inside {fio.REPO / 'data'} (gitignored); got {paths[k]}")
    try:
        run(cfg, tasks, paths["raw"], paths["dst"], paths["old"], a.keep_raw)
    except BaseException:
        if paths["raw"].exists():
            print(f"FAILED; the raw files are still in {paths['raw']}: rerun, or delete that folder", flush=True)
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())

