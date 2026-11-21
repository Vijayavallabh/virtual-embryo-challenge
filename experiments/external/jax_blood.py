"""Row T66 (experiments/harness/TRIES.md): the erythroid direction of JAX/Qiu 2024 (register row X003), recomputed from
raw counts.

The kill test of row T66 (experiments/t1_m1/jaxblood_read.py) moves the drawn Blood cells by the per-gene change of an
erythroid group of JAX nuclei between an early end (embryos 1-12, E8.5) and a late end (embryos 59-63, E15.5-E16.5).
This script makes that change and nothing else, with the fixed settings of experiments/external/jax_blood.yaml:

  1. The group is recomputed from the nuclei's raw counts by row N54's marker-score machinery, imported from
     experiments/external/feng_relabel.py (log1p_cp10k, unique_index, control_sets, set_score, assign; none of it is
     copied), with N54's released-data gene sets. The controls' gene means are taken per end over the end's nuclei and
     averaged over the two ends with equal weight (feng_relabel.gene_means's rule, with the ends as its stages), summed
     embryo by embryo so that no end is held in memory as log values at once (`end_means`, tested equal to gene_means).
  2. Each assigned erythroid nucleus with at least umi_target UMIs is downsampled to exactly umi_target UMIs (without
     replacement); a nucleus with fewer is left out. Per embryo with at least min_embryo_nuclei such nuclei, the mean
     CP10k (vec.core.mean_cp10k); per end, the mean over its embryos with equal weight.
  3. direction[g] = log((m_late[g] + eps) / (m_early[g] + eps)) per JAX gene (dirlib.lfc_ratio's form).

What it reads, and nothing else: the per-embryo count matrices of the row's embryos (each checked against the JAX
folder's MANIFEST.json), from meta/df_cell.csv the columns of CELL_WHITELIST alone (a reader that refuses any other
column, and refuses a whitelist naming an annotation, cluster or embedding column; celltype_update and
major_trajectory are refused by name), and from meta/df_gene.csv gene_id and gene_short_name. JAX's own labels and the
stats/*.npz files keyed by them are never read. An embryo outside ROW_EMBRYOS refuses before any file is opened.

It runs once (`preflight`, all checked before any JAX file is opened): ledger row T66 in TRIES.md; X003's row in
EXTERNAL_DATA.md amended to name this use (AMENDMENT_PHRASE); the committed yaml with sha256 YAML_SHA256 and status
frozen, equal to its version at HEAD; a clean tree at the commit given by --commit; the output folder absent and inside
the checkout's own data/ directory. It writes STARTED.json before the first JAX file is opened, then the direction file,
the per-embryo table and MANIFEST.json (inputs' and outputs' sha256, the yaml's sha256, the columns read, the embryos
read, the gate, per-embryo counts, the sets' genes found and controls). It never removes anything; a folder left by a
stopped run blocks a second run (a rerun is a new row).

Run (DGX campaign checkout, repo root, through experiments/t1_m1/jaxblood.sh's `reduce` step; rdata from the external
dependency group):
  uv run --frozen --no-sync python -I -B experiments/external/jax_blood.py --commit <row commit>
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import feng2022_geo  # noqa: E402  (the rdata reader's dgCMatrix helper)
import feng_relabel as fr  # noqa: E402  (row N54's marker-score machinery and gate helpers)
import filterio as fio  # noqa: E402

from vec import provenance  # noqa: E402
from vec.core import mean_cp10k  # noqa: E402

SOURCE = "jax_blood"
ROW = "T66"
CONFIG = fio.HERE / "jax_blood.yaml"
YAML_SHA256 = "750cd3c3c15e0d6cf20e02b219a99652bfa5e36dada9376078db1878be7f6064"
TRIES = fio.REPO / "experiments" / "harness" / "TRIES.md"
REGISTER = fio.REPO / "EXTERNAL_DATA.md"
REGISTER_ROW = "X003"
AMENDMENT_PHRASE = "model input for the held Blood direction"
CODE = [Path(__file__), Path(fr.__file__), Path(feng2022_geo.__file__), Path(fio.__file__), fr.PANEL_SOURCE]
STARTED = "STARTED.json"
ROW_EMBRYOS = {"early": tuple(range(1, 13)), "late": (59, 60, 61, 62, 63)}
ENDS = ("early", "late")
# The only df_cell columns this module may read, fixed in code; the cell id is read by name, not by position.
CELL_WHITELIST = ("cell_id", "embryo_id", "day", "experimental_batch")
GENE_COLUMNS = ("gene_id", "gene_short_name")
# Named author columns refused whatever else the pattern finds (the row's whitelist check)
REFUSED_COLUMNS = ("celltype_update", "major_trajectory")
_EMBRYO = re.compile(r"^embryo_(\d+)$")


class Refused(RuntimeError):
    """A rule of jax_blood.yaml or of this module stops the run."""


# ------------------------------------------------------------------------------------------------ config


def load_config(path: str | Path = CONFIG) -> dict:
    cfg = yaml.safe_load(Path(path).read_text())
    check_config(cfg)
    return cfg


def check_config(cfg: Mapping) -> None:
    """Structural checks that need no data: the row, the embryo lists, the sets' sources and the settings."""
    if cfg.get("row") != ROW:
        raise Refused(f"the yaml is not row {ROW}'s")
    check_embryos(cfg["embryos"])
    for name, s in cfg["sets"].items():
        if s.get("role") not in fr.ROLES:
            raise Refused(f"set {name}: role {s.get('role')!r} is not one of {fr.ROLES}")
        bad = [g for g, src in s["genes"].items() if src != "released_e95"]
        if bad:
            raise Refused(f"set {name}: genes {bad} are not sourced to the released E9.5 label means")
        if len(s["genes"]) < int(cfg["score"]["min_genes_found"]):
            raise Refused(f"set {name} has fewer genes than min_genes_found")
    if cfg["sets"].get("blood", {}).get("label") != cfg["erythroid_label"]:
        raise Refused("the erythroid label is not the blood set's label")
    for k in ("umi_target", "min_embryo_nuclei"):
        if int(cfg[k]) < 1:
            raise Refused(f"{k} must be at least 1")
    if float(cfg["eps"]) <= 0:
        raise Refused("eps must be positive")


def check_embryos(embryos: Mapping[str, Iterable[int]]) -> dict[str, tuple[int, ...]]:
    """The yaml's lists must be the row's lists exactly (ROW_EMBRYOS): an embryo outside them refuses."""
    got = {e: tuple(sorted(int(x) for x in embryos.get(e, ()))) for e in ENDS}
    if set(embryos) != set(ENDS) or got != ROW_EMBRYOS:
        raise Refused(f"the embryo lists {dict(embryos)} are not row {ROW}'s {ROW_EMBRYOS}")
    return got


def embryo_number(embryo_id: str) -> int:
    m = _EMBRYO.match(str(embryo_id).strip())
    if not m:
        raise Refused(f"malformed embryo_id {embryo_id!r}")
    return int(m.group(1))


def end_of(n: int) -> str:
    for e in ENDS:
        if n in ROW_EMBRYOS[e]:
            return e
    raise Refused(f"embryo {n} is outside row {ROW}'s lists")


def scoring_config(cfg: Mapping) -> dict:
    """The parts of the yaml feng_relabel's assign, classes and named_genes read."""
    return {"sets": cfg["sets"], "assign": cfg["assign"], "dropped": list(cfg.get("dropped") or [])}


# ------------------------------------------------------------------------------------------------ the gate


def register_amended(path: Path | None = None) -> bool:
    """X003's row in EXTERNAL_DATA.md names this use (AMENDMENT_PHRASE)."""
    for line in Path(path or REGISTER).read_text(encoding="utf-8").splitlines():
        if line.startswith(f"| {REGISTER_ROW} |"):
            return AMENDMENT_PHRASE in " ".join(line.split())
    return False


def preflight(cfg_path: Path, commit: str, dst: Path) -> dict:
    """Everything the run needs, checked before any JAX file is opened. Returns the gate the manifest records."""
    if not fr.row_in_tries(ROW, TRIES):
        raise Refused(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if not register_amended():
        raise Refused(f"{REGISTER_ROW}'s row in EXTERNAL_DATA.md does not name this use ({AMENDMENT_PHRASE!r}): the "
                      "register is amended before any read")
    if Path(cfg_path).resolve() != CONFIG.resolve():
        raise Refused(f"only the committed {CONFIG} is run, not {cfg_path}")
    raw = Path(cfg_path).read_bytes()
    if fio.sha256_file(cfg_path) != YAML_SHA256:
        raise Refused(f"the yaml's sha256 is not the row's {YAML_SHA256}")
    cfg = yaml.safe_load(raw)
    if cfg.get("status") != "frozen":
        raise Refused(f"the yaml's status is {cfg.get('status')!r}: it is run only when 'frozen'")
    check_config(cfg)
    if fr.committed_at_head(CONFIG) != raw:
        raise Refused("the yaml differs from its committed version at HEAD")
    head, unclean = provenance.git_state(fio.REPO)
    if unclean is None or unclean:
        raise Refused(f"the tree is not clean ({unclean}): the reducer runs once, on a committed tree")
    if head != commit:
        raise Refused(f"HEAD is {head}, not the run commit {commit}")
    if Path(dst).exists():
        raise Refused(f"{dst} exists: the reducer never writes over a folder (a rerun is a new row)")
    return {"row": ROW, "commit": head, "clean": True, "yaml_sha256": YAML_SHA256}


# ------------------------------------------------------------------------------------------------ the readers


def check_columns(columns: Iterable[str]) -> tuple[str, ...]:
    """The df_cell columns to read must lie inside CELL_WHITELIST; none may be a refused author column or look like
    an annotation, cluster or embedding column (feng_relabel.ANNOTATION_PATTERN)."""
    cols = tuple(columns)
    bad = [c for c in cols if c not in CELL_WHITELIST or c in REFUSED_COLUMNS or fr.ANNOTATION_PATTERN.search(c)]
    if bad:
        raise Refused(f"columns {bad} may not be read: only {list(CELL_WHITELIST)} are read")
    return cols


def read_cells(path: str | Path, embryos: Iterable[int], columns: Sequence[str] = CELL_WHITELIST) -> dict:
    """The whitelisted columns of df_cell, for the rows of `embryos` alone: each row's other fields are dropped as it
    is parsed, rows of other embryos are dropped unread past their embryo_id, and a requested column outside the
    whitelist refuses before the file is opened. Returns {column: values} and the rows dropped as other embryos."""
    cols = check_columns(columns)
    want = {int(e) for e in embryos}
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        missing = [c for c in cols if c not in header]
        if missing:
            raise Refused(f"{Path(path).name}: columns {missing} are not in its header")
        idx = [header.index(c) for c in cols]
        e_idx = header.index("embryo_id")
        rows, other = [], 0
        for i, row in enumerate(reader):
            if len(row) != len(header):
                raise Refused(f"{Path(path).name}: data row {i} has {len(row)} fields, header {len(header)}")
            if embryo_number(row[e_idx]) not in want:
                other += 1
                continue
            rows.append([row[j] for j in idx])
    arr = np.array(rows, dtype=str).reshape(len(rows), len(cols))
    return {"columns": {c: arr[:, k] for k, c in enumerate(cols)}, "rows_other_embryos": other}


def read_genes(path: str | Path) -> tuple[list[str], list[str]]:
    """gene_id and gene_short_name of the JAX gene table (no cells), by name."""
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        missing = [c for c in GENE_COLUMNS if c not in header]
        if missing:
            raise Refused(f"{Path(path).name}: columns {missing} are not in its header")
        i, j = (header.index(c) for c in GENE_COLUMNS)
        rows = [(r[i], r[j]) for r in reader]
    return [a for a, _ in rows], [b for _, b in rows]


def load_rds(path: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """One embryo's raw counts from its .rds (a dgCMatrix, genes x nuclei), as integer CSR nuclei x genes, with the
    gene ids (row names) and nucleus names (column names). Anything else refuses."""
    from rdata.parser import parse_file

    obj = feng2022_geo.deref(parse_file(path).object)
    X, genes, cells = feng2022_geo.dgc_matrix(obj)
    if genes is None or cells is None:
        raise Refused(f"{path.name}: the count matrix has no row or column names")
    C = X.T.tocsr()
    if C.nnz and not np.allclose(C.data, np.round(C.data)):
        raise Refused(f"{path.name}: the counts are not integers")
    C.data = np.round(C.data).astype(np.int32)
    return C, list(genes), list(cells)


# ------------------------------------------------------------------------------------------------ the arithmetic


def end_means(sums: Mapping[str, np.ndarray], counts: Mapping[str, int]) -> np.ndarray:
    """feng_relabel.gene_means from per-end sums: each gene's mean log1p CP10k per end, averaged over the ends with
    equal weight. Refuses an end with no nucleus."""
    per = []
    for e in ENDS:
        if not counts.get(e):
            raise Refused(f"no nucleus at the {e} end")
        per.append(np.asarray(sums[e], dtype=np.float64) / counts[e])
    return np.mean(per, axis=0)


def set_indices(cfg: Mapping, symbols: Sequence[str]) -> tuple[dict, dict, dict]:
    """Per set, its genes found among JAX's unique symbols (feng_relabel.unique_index) and their columns; refuses a set
    with fewer than min_genes_found."""
    uidx = fr.unique_index(list(symbols))
    found, set_idx = {}, {}
    for name, s in cfg["sets"].items():
        f = [g for g in s["genes"] if g in uidx]
        if len(f) < int(cfg["score"]["min_genes_found"]):
            raise Refused(f"set {name}: {len(f)} of its genes are in JAX's unique symbols; the yaml needs "
                          f"{cfg['score']['min_genes_found']}")
        found[name], set_idx[name] = f, [uidx[g] for g in f]
    return found, set_idx, uidx


def controls(cfg: Mapping, symbols: Sequence[str], means: np.ndarray, set_idx: Mapping, uidx: Mapping) -> dict:
    """feng_relabel.control_sets with N54's exclusions: every set gene, every duplicated symbol, every gene the yaml
    names (feng_relabel.named_genes) and every panel gene (feng_relabel.panel_genes)."""
    sc = cfg["score"]
    dup = set(range(len(symbols))) - set(uidx.values())
    kept_out = fr.named_genes(scoring_config(cfg)) | fr.panel_genes()
    exclude = {i for idx in set_idx.values() for i in idx} | dup | {uidx[g] for g in kept_out if g in uidx}
    ctrl = fr.control_sets(means, list(symbols), set_idx, exclude, int(sc["control_bins"]),
                           int(sc["control_per_gene"]), int(sc["control_seed"]))
    empty = [n for n, c in ctrl.items() if len(c) == 0]
    if empty:
        raise Refused(f"sets {empty} drew no control gene")
    return ctrl


def label_nuclei(L: sparse.csr_matrix, cfg: Mapping, set_idx: Mapping, ctrl: Mapping) -> np.ndarray:
    """Every nucleus's recomputed label (feng_relabel.set_score per set, feng_relabel.assign)."""
    scores = {name: fr.set_score(L, set_idx[name], ctrl[name]) for name in cfg["sets"]}
    return fr.assign(scores, scoring_config(cfg))


def downsample(C: sparse.csr_matrix, target: int, rng: np.random.Generator) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Each row of integer counts C with at least `target` UMIs, downsampled to exactly `target` UMIs without
    replacement (multivariate hypergeometric over its non-zero entries, in row order from one generator); rows with
    fewer are left out. Returns the downsampled rows and the boolean mask of the rows kept."""
    C = sparse.csr_matrix(C, dtype=np.int64)
    C.sort_indices()
    lib = np.asarray(C.sum(1)).ravel()
    keep = lib >= target
    data = []
    for i in np.flatnonzero(keep):
        a, b = C.indptr[i], C.indptr[i + 1]
        data.append(rng.multivariate_hypergeometric(C.data[a:b], int(target), method="count"))
    rows = np.flatnonzero(keep)
    indptr = np.concatenate([[0], np.cumsum(np.diff(C.indptr)[rows])])
    indices = np.concatenate([C.indices[C.indptr[i]:C.indptr[i + 1]] for i in rows]) if rows.size else np.zeros(
        0, dtype=np.int32)
    vals = np.concatenate(data) if data else np.zeros(0, dtype=np.int64)
    out = sparse.csr_matrix((vals, indices, indptr), shape=(rows.size, C.shape[1]))
    out.eliminate_zeros()
    return out, keep


def embryo_mean(C_ery: sparse.csr_matrix, target: int, seed: Sequence[int]) -> tuple[np.ndarray | None, int]:
    """The mean CP10k (vec.core.mean_cp10k) over one embryo's erythroid nuclei after downsampling to `target` UMIs, and
    the number of nuclei that entered (None when none did)."""
    ds, keep = downsample(C_ery, target, np.random.default_rng(list(seed)))
    if ds.shape[0] == 0:
        return None, 0
    return mean_cp10k(ds), int(keep.sum())


def end_direction(per_embryo: Mapping[int, tuple[np.ndarray | None, int]], min_nuclei: int, eps: float) -> dict:
    """Per end, the equal-weight mean over its embryos with at least min_nuclei entering nuclei; the direction
    log((m_late + eps) / (m_early + eps)). An end with no such embryo gives no direction (None)."""
    m, used = {}, {}
    for e in ENDS:
        rows = [v for n, (v, k) in sorted(per_embryo.items()) if end_of(n) == e and v is not None and k >= min_nuclei]
        used[e] = [n for n, (v, k) in sorted(per_embryo.items())
                   if end_of(n) == e and v is not None and k >= min_nuclei]
        m[e] = np.mean(rows, axis=0) if rows else None
    if m["early"] is None or m["late"] is None:
        return {"m_early": m["early"], "m_late": m["late"], "direction": None, "embryos_used": used}
    d = np.log((m["late"] + eps) / (m["early"] + eps))
    return {"m_early": m["early"], "m_late": m["late"], "direction": d, "embryos_used": used}


# ------------------------------------------------------------------------------------------------ the run


def _rel(p: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(fio.REPO))
    except ValueError:
        return str(p)


def reduce(cfg: Mapping, src: Path, gate: Mapping, loader=load_rds, verify=fio.verify_manifest_output) -> dict:
    """The two passes over the row's embryos (module docstring): returns the arrays and the record. `loader` reads one
    embryo's counts (load_rds; the tests hand a synthetic one); `verify` checks an input against its manifest."""
    if not gate or gate.get("row") != ROW or gate.get("clean") is not True:
        raise Refused(f"no gate of row {ROW} on a clean tree: the reducer runs only through main(), after preflight")
    embryos = check_embryos(cfg["embryos"])
    inp = cfg["inputs"]
    wanted = sorted(n for e in ENDS for n in embryos[e])
    paths = {"cells": src / inp["cells"], "genes": src / inp["genes"],
             **{f"embryo_{n}": src / inp["counts"].format(n=n) for n in wanted}}
    inputs = {k: {"file": _rel(p), "sha256": verify(p).get("sha256")} for k, p in paths.items()}
    cells = read_cells(paths["cells"], wanted)
    col = cells["columns"]
    gene_ids, symbols = read_genes(paths["genes"])
    sym_of = dict(zip(gene_ids, symbols, strict=True))
    cfg_s = scoring_config(cfg)
    # pass 1: counts in memory, per-end sums of log1p CP10k
    store, info = {}, {}
    sums, counts = {}, {e: 0 for e in ENDS}
    genes_first = None
    for n in wanted:
        C, genes, names = loader(paths[f"embryo_{n}"])
        if genes_first is None:
            genes_first = list(genes)
        elif list(genes) != genes_first:
            raise Refused(f"embryo {n}: its genes differ from the first embryo's")
        mine = np.flatnonzero(np.char.equal(col["embryo_id"], f"embryo_{n}"))
        listed = {str(c): i for c, i in zip(col["cell_id"][mine], mine, strict=True)}
        pos = [j for j, c in enumerate(names) if c in listed]
        C = sparse.csr_matrix(C)[pos]
        nm = [names[j] for j in pos]
        order = np.argsort(np.asarray(nm, dtype=str), kind="stable")
        C, nm = C[order], [nm[j] for j in order]
        lib = np.asarray(C.sum(1)).ravel()
        empty = lib <= 0
        C, nm = C[~empty], [x for x, z in zip(nm, empty, strict=True) if not z]
        rows = [listed[x] for x in nm]
        e = end_of(n)
        L = fr.log1p_cp10k(C)
        sums[e] = sums.get(e, 0.0) + np.asarray(L.sum(0)).ravel()
        counts[e] += C.shape[0]
        store[n] = C
        info[n] = {"end": e, "nuclei_in_matrix": len(names), "nuclei_listed": int(mine.size),
                   "nuclei_read": int(C.shape[0]), "empty_dropped": int(empty.sum()),
                   "days": sorted(set(col["day"][rows].tolist())),
                   "experimental_batch": sorted(set(col["experimental_batch"][rows].tolist()))}
    missing = [g for g in genes_first if g not in sym_of]
    if missing:
        raise Refused(f"{len(missing)} genes of the count matrices are not in the gene table, e.g. {missing[:3]}")
    jsym = [sym_of[g] for g in genes_first]
    found, set_idx, uidx = set_indices(cfg, jsym)
    means = end_means(sums, counts)
    ctrl = controls(cfg, jsym, means, set_idx, uidx)
    # pass 2: labels, erythroid nuclei, downsampled means
    per_embryo, labels_count = {}, {}
    for n in wanted:
        C = store.pop(n)
        lab = label_nuclei(fr.log1p_cp10k(C), cfg, set_idx, ctrl)
        ery = lab == cfg["erythroid_label"]
        labels_count[n] = {lab_: int(np.sum(lab == lab_)) for lab_ in fr.all_labels(cfg_s)}
        v, k = embryo_mean(C[np.flatnonzero(ery)], int(cfg["umi_target"]), [int(cfg["downsample_seed"]), n])
        per_embryo[n] = (v, k)
        info[n].update(erythroid_assigned=int(ery.sum()), erythroid_entering=k)
    res = end_direction(per_embryo, int(cfg["min_embryo_nuclei"]), float(cfg["eps"]))
    gnames = np.asarray(jsym, dtype=str)
    record = {"embryos_read": wanted, "embryo_lists": {e: list(embryos[e]) for e in ENDS},
              "columns_read": list(CELL_WHITELIST), "gene_columns_read": list(GENE_COLUMNS),
              "authors_labels_read": False, "stats_npz_read": False, "inputs": inputs,
              "rows_other_embryos_dropped": cells["rows_other_embryos"], "per_embryo": info,
              "labels_per_embryo": labels_count, "embryos_used": res["embryos_used"],
              "sets": {k: {"found": found[k], "controls": [str(gnames[i]) for i in ctrl[k]]} for k in cfg["sets"]},
              "has_direction": res["direction"] is not None}
    arrays = {"gene_id": np.asarray(genes_first, dtype=str), "symbol": gnames,
              "m_early": res["m_early"] if res["m_early"] is not None else np.full(len(gnames), np.nan),
              "m_late": res["m_late"] if res["m_late"] is not None else np.full(len(gnames), np.nan),
              "direction": res["direction"] if res["direction"] is not None else np.full(len(gnames), np.nan),
              "embryo": np.array(sorted(per_embryo), dtype=np.int64),
              "embryo_mean_cp10k": np.vstack([per_embryo[n][0] if per_embryo[n][0] is not None
                                              else np.full(len(gnames), np.nan) for n in sorted(per_embryo)]),
              "embryo_nuclei": np.array([per_embryo[n][1] for n in sorted(per_embryo)], dtype=np.int64)}
    return {"arrays": arrays, "record": record}


def write_outputs(dst: Path, cfg: Mapping, cfg_path: Path, gate: Mapping, res: Mapping) -> dict:
    out = cfg["output"]
    fio.save_npz_atomic(dst / out["direction"], **res["arrays"])
    tmp = dst / f".{out['embryos']}.partial"
    with open(tmp, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["embryo", "end", "days", "experimental_batch", "nuclei_read", "erythroid_assigned",
                    "erythroid_entering"])
        for n, r in sorted(res["record"]["per_embryo"].items()):
            w.writerow([n, r["end"], ";".join(r["days"]), ";".join(r["experimental_batch"]), r["nuclei_read"],
                        r["erythroid_assigned"], r["erythroid_entering"]])
    tmp.replace(dst / out["embryos"])
    manifest = {"source": SOURCE, "row": ROW, "register_row": REGISTER_ROW,
                "config": {"file": _rel(cfg_path), "sha256": fio.sha256_file(cfg_path), "status": cfg.get("status")},
                "gate": dict(gate), **res["record"],
                "settings": {k: cfg[k] for k in ("umi_target", "downsample_seed", "min_embryo_nuclei", "eps",
                                                 "erythroid_label")},
                "outputs": {out["direction"]: {"sha256": fio.sha256_file(dst / out["direction"])},
                            out["embryos"]: {"sha256": fio.sha256_file(dst / out["embryos"])}},
                "provenance": fio.provenance(CODE)}
    fio.write_manifest(dst, manifest)
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--commit", required=True, help="the run commit (full hash): HEAD must be it")
    a = ap.parse_args(argv)
    try:  # nothing of JAX's is opened in this block
        cfg = load_config(CONFIG)
        dst = fr.output_dir_ok(fio.REPO / cfg["output"]["dir"])
        gate = preflight(CONFIG, a.commit, dst)
    except (Refused, fr.Refused) as e:
        raise SystemExit(f"refused: {e}") from None
    dst.mkdir(parents=True)
    (dst / STARTED).write_text(json.dumps({"gate": gate, "config_sha256": YAML_SHA256}, indent=1) + "\n")
    res = reduce(cfg, (fio.REPO / cfg["inputs"]["dir"]), gate)
    write_outputs(dst, cfg, CONFIG, gate, res)
    print(f"written {_rel(dst)}")  # no count
    return 0


if __name__ == "__main__":
    sys.exit(main())
