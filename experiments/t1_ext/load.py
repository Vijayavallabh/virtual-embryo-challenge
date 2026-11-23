#!/usr/bin/env python
"""Load one fetched Task 1 outside atlas: window removal FIRST, raw counts only, a fixed cell call, a QC record.

Order, per atlas (D/raw/<ACC>/MANIFEST.json, written by fetch.py):
  1. The atlas entry is the one the plan was made from, the atlas is fetchable, and atlases.yaml is frozen (main()).
  2. WINDOW FIRST: every sample's stage is read again from its own full record stored in the manifest (the same
     record the plan read and the readers pinned; t1x_common.record_stage). A sample that reads inside (E9.5, E13.5]
     is removed before any of its files is opened; a sample whose stage cannot be established refuses the atlas.
     PINS: every sample of the manifest must be pinned in atlases.yaml, its stored record hashing to its pin and its
     end and stage its pin's; every series record must hash to its pin; the pinned samples must cover both ends.
     The canonical processing text, rebuilt from the manifest's records, must hash to the atlas's processing_cleared.
     Only then is a file opened. Every cell carries its sample's stage, and drop_window_cells
     removes (and counts) any cell in the window again before any expression value is used. The QC record holds the
     cells per stage before and after, and the cells in the window recomputed per cell from its own sample's record;
     read_reduced repeats that check on every row of cells.csv.
  3. Gates per file, before any expression value is used: the file's sha256 equals the fetch manifest's; every count is
     a non-negative integer; for 10x files every barcode ends in "-1" (cellranger aggr renumbers libraries); a
     Multiome features file must have the feature_type column, and only 'Gene Expression' rows are kept; in a dense
     dge table no row name may be a summary or annotation name (nUMI, Cluster, ...) and no row may equal the column
     sums or non-zero counts of the others. The suffix and integer gates do not detect a per-sample re-split of a
     downsampled aggr matrix: the processing clearance (by text hash) is the guard.
  4. Columns read, and nothing else (a whitelist in code; a name matching the labeller's annotation pattern can never
     be read): 10x barcodes (first column, the only one split out), features (id, symbol, feature_type), the matrix;
     Microwell-seq dge: the gene column (row names checked first), the barcode header and the counts. No author
     label, cluster, embedding, filter or normalisation is read; the QC record lists the columns read and the checks.
  5. The atlas's fixed cell call (atlases.yaml qc: UMI, genes and mitochondrial fraction), the same for every sample.
  6. Genes by symbol: duplicated symbols in a sample are dropped (as dirlib.sym_index), the atlas keeps the symbols
     every sample has.
Outputs (D/reduced/<ACC>/, outside the repository): counts.npz (CSR int32, cells x genes), genes.txt, cells.csv (cell,
gsm, stage label, stage day), qc.json, MANIFEST.json (each output's sha256, the inputs' sha256, tasks ["T1"]).

  uv run --frozen --no-sync python -I -B experiments/t1_ext/load.py GSE230531 --dest ~/v/vec-ext/t1 [--delete-raw]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import io as spio
from scipy import sparse

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import t1x_common as cm  # noqa: E402

Refused = cm.Refused
# The only columns any loader may read, per table kind. Fixed in code; no config widens it.
WHITELIST = {"barcodes": ("barcode",), "features": ("id", "symbol", "feature_type"), "dge": ("gene",)}
GENE_EXPRESSION = "Gene Expression"


def annotation_pattern() -> re.Pattern:
    import feng_relabel

    return feng_relabel.ANNOTATION_PATTERN


# ------------------------------------------------------------------------------------------------ tables


class TsvTable:
    """A headerless tab-separated file read by position: names[i] is column i. Only a whitelisted name is ever read,
    and only that column is split out of each line (split at most i times; the rest of the line is never parsed), so an
    author column beside it is never parsed (the first review's point 5)."""

    def __init__(self, path: Path, names: Sequence[str]):
        self.path, self._names = Path(path), tuple(names)

    def _open(self):
        op = gzip.open if self.path.name.endswith(".gz") else open
        return op(self.path, "rt", newline="")

    def names(self) -> tuple[str, ...]:
        with self._open() as f:
            first = next((ln for ln in f if ln.strip()), "")
        width = first.rstrip("\n").count("\t") + 1 if first else 0
        return self._names[:width]

    def column(self, name: str) -> np.ndarray:
        i = self._names.index(name)
        out = []
        with self._open() as f:
            for ln in f:
                if not ln.strip():
                    continue
                parts = ln.rstrip("\n").split("\t", i + 1)
                if len(parts) <= i:
                    raise Refused(f"{self.path.name}: a line has no column {i} ({name})")
                out.append(parts[i])
        return np.array(out, dtype=object)


class Reader:
    """Reads columns through the whitelist and records what was read."""

    def __init__(self):
        self.read: list[str] = []

    def column(self, table, kind: str, name: str) -> np.ndarray:
        if name not in WHITELIST[kind] or annotation_pattern().search(name):
            raise Refused(f"column {name!r} of a {kind} table is not on the loader's whitelist")
        self.read.append(f"{kind}:{name}")
        return np.asarray(table.column(name))


def check_counts(M) -> None:
    data = np.asarray(M.data)
    if data.size and (data.min() < 0 or not np.all(np.equal(np.mod(data, 1), 0))):
        raise Refused("a count is negative or not an integer: authors' corrected or normalised values are refused")


def check_suffix(barcodes: np.ndarray, suffix: str | None) -> dict:
    if suffix is None:
        return {"barcode_suffix": "not applicable (no library suffix in this form)"}
    bad = int(sum(not str(b).endswith(suffix) for b in barcodes))
    if bad:
        raise Refused(f"{bad} barcodes do not end in {suffix!r}: an aggregated or merged matrix is refused")
    return {"barcode_suffix": f"all {len(barcodes)} end in {suffix!r}"}


def tenx_from_tables(bt, ft, M, multiome: bool, suffix: str | None, rd: Reader) -> dict:
    """A 10x sample from its barcodes table, features table and matrix (genes x barcodes, Matrix Market)."""
    barcodes = rd.column(bt, "barcodes", "barcode").astype(str)
    gates = check_suffix(barcodes, suffix)
    symbols = rd.column(ft, "features", "symbol").astype(str)
    rd.column(ft, "features", "id")
    M = sparse.csc_matrix(M)
    if M.shape != (len(symbols), len(barcodes)):
        raise Refused(f"matrix shape {M.shape} does not match {len(symbols)} features x {len(barcodes)} barcodes")
    check_counts(M)
    keep = np.ones(len(symbols), dtype=bool)
    if multiome:
        if "feature_type" not in ft.names():
            raise Refused("a Multiome features file without the feature_type column is refused")
        keep = rd.column(ft, "features", "feature_type").astype(str) == GENE_EXPRESSION
        gates["gene_expression_rows"] = int(keep.sum())
        gates["other_rows_dropped"] = int((~keep).sum())
    C = sparse.csr_matrix(M[np.flatnonzero(keep)].T.astype(np.int64))
    gates["integer_counts"] = True
    return {"C": C, "barcodes": barcodes, "genes": symbols[keep], "gates": gates}


# Row names of a dense gene-by-cell table that are not genes (the reviews' blocker 5): per-cell summaries, author
# annotations and metadata rows, matched against the WHOLE row name, case ignored. Not checked against a full mouse
# gene list: a real gene named like one of these refuses the file, which fails closed. "time" is listed only as itself
# or "time point" (Timeless is a gene); "pca" only as itself or with a number (Pcare is a gene).
SUMMARY_ROW = re.compile(
    r"n_?umis?|n_?genes?|n_?counts?|n_?features?|n_?reads?|n_?cells?|ncount(?:[._]\S*)?|nfeature(?:[._]\S*)?|"
    r"total(?:[._ ]\S*)?|umis?(?:[._ ]?counts?)?|genes?(?:[._ ]?counts?)?|library[._ ]?size|size[._ ]?factors?|"
    r"percent(?:[._ ]?\S*)?|pct(?:[._ ]\S*)?|mito\S*|cluster\S*|cell[._ ]?types?\S*|louvain\S*|leiden\S*|"
    r"seurat\S*|batch\S*|samples?\S*|stages?\S*|orig\.ident|ident\S*|barcodes?|cell(?:[._ ]?(?:ids?|names?|"
    r"barcodes?))?|doublet\S*|scrublet\S*|replicates?\S*|time|time[._ ]?points?|age|embryos?\S*|genotypes?\S*|"
    r"sex|conditions?\S*|groups?\S*|annot\S*|labels?\S*|class\S*|umap\S*|tsne\S*|t_sne\S*|pca|pca[._ ]?\d+|"
    r"pc[._ ]?\d+|harmony\S*|embed\S*|lineage\S*|subtype\S*|scores?(?:[._ ]\S*)?|phase|\S*[._]snn\S*|snn[._]\S*|"
    r"res\.\S*|s\.score|g2m\.score", re.IGNORECASE)
# The one real mouse gene symbol the labeller's pattern matches whole (stannin); any other spelling of it refuses.
ANNOTATION_GENE_EXCEPTIONS = frozenset({"Snn"})


def check_row_names(genes: Sequence[str]) -> dict:
    """Refuse a dense table any of whose row names fully matches SUMMARY_ROW or the labeller's annotation pattern, or
    splits into two or more tokens one of which fully matches that pattern (RNA_snn_res.0.8). A substring is not
    enough: Epcam holds 'pca'."""
    pat = annotation_pattern()
    bad = []
    for g in genes:
        g = str(g).strip()
        toks = [t for t in re.split(r"[^A-Za-z0-9]+", g) if t]
        if SUMMARY_ROW.fullmatch(g) or (pat.fullmatch(g) and g not in ANNOTATION_GENE_EXCEPTIONS) or \
                (len(toks) >= 2 and any(pat.fullmatch(t) for t in toks)):
            bad.append(g)
    if bad:
        raise Refused(f"dense-table rows {bad[:5]} are summaries or annotations, not genes: the file is refused")
    return {"row_names": f"checked: none of {len(genes)} row names is a summary or annotation name"}


def check_summary_rows(V: np.ndarray, chunk: int = 2048) -> dict:
    """Refuse a dense table with a row equal to the column sums of the other rows (a UMI total) or, when not all
    zero, to the number of non-zero other rows per column (a gene count)."""
    V = np.asarray(V, dtype=np.float64)
    if V.shape[0] < 3:
        return {"summary_rows": "not checked (fewer than three rows)"}
    tot, nz = V.sum(0), (V > 0).sum(0)
    for i in range(0, V.shape[0], chunk):
        X = V[i:i + chunk]
        eq_sum = np.all(2.0 * X == tot[None, :], axis=1)
        eq_nz = np.all(X == (nz[None, :] - (X > 0)), axis=1) & np.any(X != 0, axis=1)
        hit = np.flatnonzero(eq_sum | eq_nz)
        if hit.size:
            raise Refused(f"dense-table row {i + int(hit[0])} equals the column sums or non-zero counts of the other "
                          f"rows: a summary row, the file is refused")
    return {"summary_rows": "checked: no row equals the column sums or the non-zero counts of the other rows"}


def dge_from_table(dt, rd: Reader) -> dict:
    """A Microwell-seq dge sample: dt.names() is the header (the gene column, then one name per barcode); every
    non-gene column must be a barcode, so a column named like an annotation refuses the file before it is read. Every
    row name is checked before any count is parsed, and no row may be a per-cell summary."""
    names = list(dt.names())
    if not names or names[0] != "gene":
        raise Refused("a dge table's first column must be the gene column")
    pat = annotation_pattern()
    bad = [n for n in names[1:] if pat.search(n) or not re.fullmatch(r"[A-Za-z0-9._:-]+", n)]
    if bad:
        raise Refused(f"dge columns {bad[:3]} are not barcodes: the file is refused")
    genes = rd.column(dt, "dge", "gene").astype(str)
    gates = check_row_names(genes)
    rd.read.append("dge:row names (checked against the summary and annotation names before any count)")
    V = np.asarray(dt.counts(), dtype=np.float64)
    if V.shape != (len(genes), len(names) - 1):
        raise Refused("dge counts do not match its gene column and header")
    gates.update(check_summary_rows(V))
    rd.read.append("dge:counts")
    rd.read.append("dge:header barcodes")
    M = sparse.csr_matrix(V)  # genes x barcodes
    del V
    check_counts(M)
    gates.update({"barcode_suffix": "not applicable (no library suffix in this form)", "integer_counts": True})
    return {"C": sparse.csr_matrix(M.T.astype(np.int64)), "barcodes": np.array(names[1:], dtype=str),
            "genes": genes, "gates": gates}


class DgeFile:
    """A Microwell-seq dge.txt(.gz): a header of barcodes (with or without a leading gene-column name), then one row
    per gene. Separator: tab when the header holds one, else whitespace. Counts are parsed only when counts() is
    called, after the row names were checked; a value that is not a number refuses the file."""

    def __init__(self, path: Path):
        self.path = Path(path)
        op = gzip.open if self.path.name.endswith(".gz") else open
        with op(self.path, "rt") as f:
            text = f.read()
        lines = [ln for ln in text.split("\n") if ln.strip()]
        sep = "\t" if "\t" in lines[0] else None
        head = lines[0].split(sep)
        first = lines[1].split(sep) if len(lines) > 1 else []
        if len(head) == len(first):
            head = head[1:]
        elif len(head) != len(first) - 1:
            raise Refused(f"{self.path.name}: header and rows do not line up")
        self._names = ["gene", *[h.strip().strip('"') for h in head]]
        self._rows = [ln.split(sep) for ln in lines[1:]]
        self._genes = np.array([r[0].strip('"') for r in self._rows], dtype=object)
        self._width = len(head)
        del text, lines

    def names(self):
        return tuple(self._names)

    def column(self, name):
        if name != "gene":
            raise Refused(f"column {name!r}")
        return self._genes

    def counts(self):
        if not self._rows:
            return np.zeros((0, self._width))
        try:
            return np.array([[float(x) for x in r[1:]] for r in self._rows], dtype=np.float64)
        except ValueError:
            raise Refused(f"{self.path.name}: a value that is not a number: the file is refused") from None


# ------------------------------------------------------------------------------------------------ window first


STAGE_RECORD_KEYS = ("record", "series_protocol", "label")
# The first review's replacement text for the QC record, extended by what this code does.
IN_WINDOW_HOW = ("recomputed per cell from its sample's stored full record (the whole record re-read by the stage "
                 "rule, the plug-day widening applied, the point compared with the cell's day) and from the day "
                 "itself; equals 0 by construction of window_filter; it does not test the stage rule itself. "
                 "read_reduced repeats it on every row of cells.csv.")


def sample_stage(gsm: str, s: Mapping) -> dict:
    """The stage rule on one sample's stored full record (the record the plan read and the readers pinned) and its
    series' summary and design (for the plug-day convention)."""
    missing = [k for k in STAGE_RECORD_KEYS if k not in s]
    if missing:
        raise Refused(f"{gsm}: the manifest lacks {missing} of its own record: the atlas is refused")
    # the series-wide fields (ruling row R27's exception) as the plan computed them; absent, no field is series text
    return cm.record_stage(s["record"], s["series_protocol"], s.get("series_text"))


def stage_record(s: Mapping) -> dict:
    """What the reduction keeps of a sample to re-read its stage: its record, its series' text and its label (and the
    series-wide fields when the fetch stored them)."""
    out = {k: s[k] for k in STAGE_RECORD_KEYS}
    if "series_text" in s:
        out["series_text"] = s["series_text"]
    return out


def window_filter(samples: Mapping[str, Mapping]) -> tuple[dict, dict]:
    """Stage FIRST: each sample's stage read again from its own stored record. A WINDOW sample is removed here (none of
    its files is opened); a sample that is not PERMITTED, whose readings are not one point, or whose point is not its
    label refuses the atlas. Returns the kept samples and the record."""
    kept, removed = {}, []
    for gsm, s in samples.items():
        r = sample_stage(gsm, s)
        if r["mark"] == "WINDOW":
            removed.append({"gsm": gsm, "low": r["low"], "high": r["high"]})
            continue
        if r["mark"] != "PERMITTED" or r["window_reading"] or r["nominal_low"] != r["nominal_high"]:
            raise Refused(f"{gsm}: its stage cannot be established from its own record ({r['reason'] or 'a range'})")
        if cm.label_day(s["label"]) != r["nominal_low"]:
            raise Refused(f"{gsm}: its record reads E{r['nominal_low']:g}, the manifest says {s['label']}")
        kept[gsm] = {**s, "day": r["nominal_low"]}
    return kept, {"samples_removed_window": removed, "samples_kept": sorted(kept)}


def check_pinned(acc: str, a: Mapping, man: Mapping) -> dict:
    """The first line, again at load: every sample of the manifest pinned, its stored record hashing to its pin (and to
    the manifest's record_sha256), its end and stage its pin's; the samples cover both ends; every series record of
    the atlas hashing to its pin. A sample that is not pinned is never loaded: the atlas is refused."""
    pins = cm.check_pins(acc, a)
    if not cm.pins_cover_both_ends(pins):
        raise Refused(f"{acc}: the pinned samples do not cover both ends: nothing is loaded")
    for g, s in man["samples"].items():
        p = pins["samples"].get(g)
        if p is None:
            raise Refused(f"{acc}: {g} is not pinned: the atlas is refused")
        h = cm.record_sha256("SAMPLE", g, s.get("record", []))
        if not (h == s.get("record_sha256") == p["record_sha256"]) or p["end"] != s.get("end") or \
                p["stage"] != s.get("label"):
            raise Refused(f"{acc}: {g}: its stored record, end or stage is not its pin's: the atlas is refused")
    if {s.get("end") for s in man["samples"].values()} != {"early", "late"}:
        raise Refused(f"{acc}: the manifest's samples do not cover both ends")
    recs = man.get("series_records") or {}
    for sr in a["series"]:
        if sr not in recs or sr not in pins["series"] or \
                cm.record_sha256("SERIES", sr, recs[sr]) != pins["series"][sr]["record_sha256"]:
            raise Refused(f"{acc}: series {sr}: its stored record does not hash to its pin")
    return pins


def check_processing(acc: str, a: Mapping, man: Mapping) -> str:
    """The canonical processing text rebuilt from the fetch manifest must hash to the atlas's processing_cleared."""
    if "series_records" not in man or any("record" not in s for s in man["samples"].values()):
        raise Refused(f"{acc}: the fetch manifest holds no records: the atlas is refused")
    st = man.get("series_text")
    for g, s in man["samples"].items():
        if "series_text" in s and (st is None or s["series_text"] != st.get(s["series"], {})):
            raise Refused(f"{acc}: {g}: its series-wide fields are not its series' in the manifest")
    text = cm.processing_text(acc, [(s, man["series_records"].get(s, [])) for s in a["series"]],
                              {g: {"series": s["series"], "text": cm.sample_text(s["record"])}
                               for g, s in man["samples"].items()}, st)
    sha = cm.text_sha256(text)
    ok, why = cm.clearance(a, sha)
    if not ok:
        raise Refused(f"{acc}: {why}")
    return sha


def check_screen(acc: str, a: Mapping, man: Mapping) -> dict | None:
    """For a screened atlas: atlases.yaml records its screen kept, the manifest was fetched under that screen, and
    every fetched barcodes and features file the manifest records has the sha256 of the screened file (row T72, cell 3;
    read_sample re-hashes each file against the manifest before it is opened)."""
    e = cm.screen_kept(acc, a)
    if e is None:
        return None
    if man.get("screen") != e:
        raise Refused(f"{acc}: the files were not fetched under the screen atlases.yaml records")
    for g, s in man["samples"].items():
        for role in cm.SCREEN_ROLES:
            f = s["files"].get(role) or {}
            if g not in e or f.get("name") != e[g][role]["name"] or f.get("sha256") != e[g][role]["sha256"]:
                raise Refused(f"{acc}: {g} {role}: not the file the readers screened")
    return e


def cells_in_window(gsm: np.ndarray, day: np.ndarray, records: Mapping[str, Mapping]) -> int:
    """Cells in the window, recomputed per cell from its own sample's stored full record (re-read by the stage rule,
    the convention widening applied) and from its day: a cell counts when its sample is not PERMITTED, has a
    reading that touches the window, does not read as one point equal to the cell's day, or the day itself is in the
    window."""
    gsm, day = np.asarray(gsm, dtype=str), np.asarray(day, dtype=np.float64)
    bad = np.zeros(len(gsm), dtype=bool)
    for g in np.unique(gsm):
        if g not in records:
            raise Refused(f"{g}: a cell's sample has no stored record")
        r = sample_stage(g, records[g])
        m = gsm == g
        ok_sample = r["mark"] == "PERMITTED" and not r["window_reading"] and r["nominal_low"] == r["nominal_high"]
        bad[m] = True if not ok_sample else (day[m] != r["nominal_low"])
    bad |= ~drop_window_cells(day)
    return int(bad.sum())


def drop_window_cells(day: np.ndarray) -> np.ndarray:
    """Per cell: True where the cell's stage is outside the T1 window (windows.allowed with tasks T1)."""
    day = np.asarray(day, dtype=np.float64)
    if not np.isfinite(day).all():
        raise Refused("a cell has no stage")
    return np.array([not cm.touches_window(d, d) for d in day], dtype=bool)


def cell_call(C: sparse.csr_matrix, genes: np.ndarray, qc: Mapping) -> np.ndarray:
    """The atlas's fixed cell call: UMI, detected genes and mitochondrial fraction, per cell of one sample."""
    umi = np.asarray(C.sum(1)).ravel().astype(np.float64)
    ng = np.diff(C.indptr)
    mt = np.array([str(g).lower().startswith("mt-") for g in genes])
    mito = np.asarray(C[:, np.flatnonzero(mt)].sum(1)).ravel() / np.maximum(umi, 1.0) if mt.any() else \
        np.zeros(C.shape[0])
    return (umi >= float(qc["min_umi"])) & (ng >= int(qc["min_genes"])) & (mito <= float(qc["max_mito"]))


def unique_symbols(C: sparse.csr_matrix, genes: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray, int]:
    g = np.asarray(genes, dtype=str)
    vals, counts = np.unique(g, return_counts=True)
    dup = set(vals[counts > 1])
    keep = np.array([x not in dup for x in g])
    return C[:, np.flatnonzero(keep)], g[keep], int((~keep).sum())


# ------------------------------------------------------------------------------------------------ one atlas


def read_sample(acc: str, form: Mapping, a: Mapping, s: Mapping, sdir: Path, rd: Reader) -> dict:
    files = {}
    for role, f in s["files"].items():
        p = sdir / f["name"]
        if not p.exists() or cm.sha256_file(p) != f["sha256"]:
            raise Refused(f"{acc} {p.name}: missing or not the file the fetch manifest records")
        files[role] = p
    if a["form"] == "tenx_mtx":
        M = spio.mmread(str(files["matrix"]), spmatrix=False)
        return tenx_from_tables(TsvTable(files["barcodes"], WHITELIST["barcodes"]),
                                TsvTable(files["features"], WHITELIST["features"]), M,
                                bool(a.get("multiome")), form.get("barcode_suffix"), rd)
    if a["form"] == "dge_txt":
        return dge_from_table(DgeFile(files["dge"]), rd)
    raise Refused(f"unknown form {a['form']!r}")


def load_atlas(acc: str, dest: Path, cfg: Mapping | None = None, check_frozen: bool = True,
               reader=None) -> dict:
    """The reduction of one atlas (see the module docstring). `reader(acc, form, atlas, sample, folder, Reader)` is
    read_sample unless a test passes a synthetic one."""
    cfg = cfg or cm.load_atlases()
    if check_frozen:
        cm.require_frozen(cfg)
    a = cm.atlas(acc, cfg)
    if a["status"] not in cm.FETCHABLE:
        raise Refused(f"{acc} has status {a['status']!r}: it is never loaded")
    rdir = cm.writable(Path(dest) / "raw" / acc)
    man_path = rdir / "MANIFEST.json"
    man = json.loads(man_path.read_text())
    if man.get("atlas_entry_sha256") != cm.entry_sha256(acc, cfg):
        raise Refused(f"{acc}: fetched under another atlas entry")
    kept, wrec = window_filter(man["samples"])  # FIRST: before any file is opened
    check_pinned(acc, a, man)  # the first line: every sample pinned, its record its pin's, both ends
    processing_sha = check_processing(acc, a, man)  # the cleared text, rebuilt from the manifest, before any file
    screen = check_screen(acc, a, man)  # a screened member is loaded only while atlases.yaml records it kept
    records = {g: stage_record(s) for g, s in kept.items()}
    form = cfg["forms"][a["form"]]
    rd = Reader()
    parts, qc = [], {"per_sample": {}}
    for gsm, s in sorted(kept.items()):
        x = (reader or read_sample)(acc, form, a, s, rdir / gsm, rd)
        C, genes = sparse.csr_matrix(x["C"]), np.asarray(x["genes"], dtype=str)
        # the cell call and each cell's total on ALL its Gene Expression rows, before duplicated symbols are dropped
        # (row T72: CP10k is counts x 10^4 over the cell's total on Gene Expression rows)
        called = cell_call(C, genes, a["qc"])
        total = np.asarray(C.sum(1)).ravel().astype(np.int64)
        C, genes, ndup = unique_symbols(C, genes)
        day = np.full(int(called.sum()), s["day"])
        ok = drop_window_cells(day)
        qc["per_sample"][gsm] = {"stage": s["label"], "day": s["day"], "barcodes": int(C.shape[0]),
                                 "cells_called": int(called.sum()), "cells_kept": int(ok.sum()),
                                 "duplicated_symbols_dropped": ndup, "gates": x["gates"]}
        Ck = C[np.flatnonzero(called)][np.flatnonzero(ok)]
        bc = np.asarray(x["barcodes"], dtype=str)[np.flatnonzero(called)][np.flatnonzero(ok)]
        tot = total[np.flatnonzero(called)][np.flatnonzero(ok)]
        parts.append((gsm, s["label"], s["day"], Ck, genes, bc, tot))
    if not parts:
        raise Refused(f"{acc}: no sample left after the window filter")
    common = sorted(set.intersection(*[set(p[4]) for p in parts]))
    mats, cells = [], []
    for gsm, lab, day, Ck, genes, bc, tot in parts:
        idx = {g: i for i, g in enumerate(genes)}
        mats.append(Ck[:, [idx[g] for g in common]])
        cells += [(f"{gsm}:{b}", gsm, lab, day, int(t)) for b, t in zip(bc, tot, strict=True)]
    C = sparse.vstack(mats).tocsr().astype(np.int32)
    days = np.array([c[3] for c in cells], dtype=np.float64)
    before = {lab: int(sum(1 for c in cells if c[2] == lab)) for lab in sorted({c[2] for c in cells})}
    keep = drop_window_cells(days)
    if not keep.all():  # cannot happen after window_filter; the count is the proof either way
        C, cells = C[np.flatnonzero(keep)], [c for c, k in zip(cells, keep, strict=True) if k]
    after = {lab: int(sum(1 for c in cells if c[2] == lab)) for lab in before}
    in_window = cells_in_window(np.array([c[1] for c in cells], dtype=str),
                                np.array([c[3] for c in cells], dtype=np.float64), records)
    qc.update({"accession": acc, "window": {**wrec, "cells_per_stage_before": before, "cells_per_stage_after": after,
                                            "cells_in_window_after": in_window,
                                            "cells_in_window_after_how": IN_WINDOW_HOW},
               "processing_sha256": processing_sha,
               "columns_read": sorted(set(rd.read)), "genes_common": len(common),
               "genes_per_sample": {p[0]: len(p[4]) for p in parts}, "qc_rule": a["qc"], "tasks": list(cm.TASKS),
               "authors_objects_read": "none: no label, cluster, embedding, filter or normalisation column exists in "
                                       "what is read, and the whitelist refuses any other column"})
    if in_window:
        raise Refused("a cell inside the window survived: nothing is written")
    return {"C": C, "genes": common, "cells": cells, "qc": qc, "fetch_manifest_sha256": cm.sha256_file(man_path),
            "chosen": man["chosen"], "stage_records": records, "screen": screen}


def write_reduced(acc: str, dest: Path, red: Mapping, cfg: Mapping | None = None) -> Path:
    import filterio as fio

    out = cm.writable(Path(dest) / "reduced" / acc)
    out.mkdir(parents=True, exist_ok=True)
    C = red["C"]
    np.savez_compressed(out / "counts.npz", data=C.data, indices=C.indices, indptr=C.indptr,
                        shape=np.array(C.shape))
    (out / "genes.txt").write_text("\n".join(red["genes"]) + "\n")
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["cell", "gsm", "stage", "day", "total_ge"])
    w.writerows(red["cells"])
    (out / "cells.csv").write_text(buf.getvalue())
    cm.dump_json(red["qc"], out / "qc.json")
    outputs = {n: {"sha256": cm.sha256_file(out / n)} for n in ("counts.npz", "genes.txt", "cells.csv", "qc.json")}
    man = {"source": "t1_ext_load", "accession": acc, "tasks": list(cm.TASKS), "outputs": outputs,
           "inputs": {"fetch_manifest_sha256": red["fetch_manifest_sha256"]}, "chosen": red["chosen"],
           "atlas_entry_sha256": cm.entry_sha256(acc, cfg), "window": red["qc"]["window"],
           "columns_read": red["qc"]["columns_read"], "stage_records": red["stage_records"],
           "processing_sha256": red["qc"]["processing_sha256"], "screen": red.get("screen"),
           "provenance": fio.provenance([Path(__file__), HERE / "t1x_common.py"])}
    return fio.write_manifest(out, man)


def check_reduction_pins(acc: str, a: Mapping, man: Mapping, cfg: Mapping) -> None:
    """THE PARKED BRANCH'S KNOWN BLOCKER, closed: a reduction is read only when it was made under the atlas entry in
    force now (pins included) and every sample it holds is pinned now, its stored record hashing to its pin and its
    label its pin's stage, the pinned samples covering both ends; its processing text is cleared now; and a screened
    member is still recorded kept. A reduction left on disk from an earlier entry, other pins, a lapsed clearance or a
    later drop is refused, not read."""
    if a["status"] not in cm.FETCHABLE:
        raise Refused(f"{acc} has status {a['status']!r}: its reduction is not read")
    if man.get("atlas_entry_sha256") != cm.entry_sha256(acc, cfg):
        raise Refused(f"{acc}: the reduction was made under another atlas entry (pins included): load again")
    pins = cm.check_pins(acc, a)
    recs = man.get("stage_records") or {}
    if not recs:
        raise Refused(f"{acc}: the reduction holds no stage records: it cannot be checked against the pins")
    for g, r in recs.items():
        p = pins["samples"].get(g)
        if p is None:
            raise Refused(f"{acc}: the reduction holds {g}, which is not pinned now")
        if cm.record_sha256("SAMPLE", g, r.get("record", [])) != p["record_sha256"] or r.get("label") != p["stage"]:
            raise Refused(f"{acc}: {g}: the reduction's record or stage is not its pin's")
    if {pins["samples"][g]["end"] for g in recs} != {"early", "late"}:
        raise Refused(f"{acc}: the reduction's samples do not cover both pinned ends")
    bad = cm.pins_complete(a, pins)
    if bad:
        raise Refused(f"{acc}: {'; '.join(bad)}")
    ok, why = cm.clearance(a, man.get("processing_sha256") or "")
    if not ok:
        raise Refused(f"{acc}: {why}")
    if cm.screen_kept(acc, a) != man.get("screen"):  # a screened member is still recorded kept, with the same files
        raise Refused(f"{acc}: the reduction was made under another screen record")


def read_reduced(acc: str, dest: Path, cfg: Mapping | None = None) -> dict:
    """The reduction, each file checked against its manifest and against the pins in force (check_reduction_pins);
    refused unless it records zero cells in the window AND every row of cells.csv, re-checked on read against its own
    sample's stored record (cells_in_window) and by drop_window_cells, is outside the window."""
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    out = cm.writable(Path(dest) / "reduced" / acc)
    man = json.loads((out / "MANIFEST.json").read_text())
    for n, e in man["outputs"].items():
        if cm.sha256_file(out / n) != e["sha256"]:
            raise Refused(f"{acc} {n}: does not match its manifest")
    if man.get("tasks") != list(cm.TASKS) or man["window"].get("cells_in_window_after") != 0:
        raise Refused(f"{acc}: the reduction does not record T1-only use with zero cells in the window")
    check_reduction_pins(acc, a, man, cfg)
    z = np.load(out / "counts.npz")
    C = sparse.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    genes = (out / "genes.txt").read_text().split()
    with open(out / "cells.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    gsm = np.array([r["gsm"] for r in rows], dtype=str)
    day = np.array([float(r["day"]) for r in rows], dtype=np.float64)
    if "stage_records" not in man:
        raise Refused(f"{acc}: the reduction holds no stage records: it cannot be re-checked")
    if not drop_window_cells(day).all() or cells_in_window(gsm, day, man["stage_records"]):
        raise Refused(f"{acc}: a row of cells.csv is in the window or does not match its sample's record")
    if not set(gsm) <= set(man["stage_records"]):
        raise Refused(f"{acc}: a row of cells.csv names a sample the reduction holds no pinned record for")
    return {"C": C, "genes": genes, "cell": np.array([r["cell"] for r in rows]),
            "gsm": gsm, "stage": np.array([r["stage"] for r in rows]), "day": day, "manifest": man,
            "total_ge": np.array([float(r.get("total_ge") or "nan") for r in rows], dtype=np.float64),
            "manifest_sha256": cm.sha256_file(out / "MANIFEST.json")}


RECHECK = ("suffix", "own_peaks", "whitelist_sized")  # (a), (b) and (d) again on the fetched files (row T72, LOAD)


def recheck_screens(acc: str, a: Mapping, cfg: Mapping, man: Mapping, rdir: Path) -> dict:
    """Checks (a), (b) and (d) of cell 3, as the atlas lists them, on the FETCHED barcodes and features files (each
    re-hashed against the manifest first); refused when one says no (a failed (p1) at the row's load)."""
    names = [n for n in a.get("screens") or [] if n in RECHECK]
    if not names:
        return {}
    files, sha = {}, {}
    for g, s in man["samples"].items():
        files[g], sha[g] = {}, {}
        for role in cm.SCREEN_ROLES:
            f = s["files"][role]
            p = rdir / g / f["name"]
            if not p.is_file() or cm.sha256_file(p) != f["sha256"]:
                raise Refused(f"{acc} {g}: {f['name']} is missing or not the file the fetch manifest records")
            files[g][role], sha[g][role] = p, f["sha256"]
    res = cm.sibling("screens").run_checks({**a, "screens": names}, cfg["screen_rules"], files, sha)
    failed = [n for n in names if not res[n]["passed"]]
    if failed:
        raise Refused(f"{acc}: checks {failed} of the screens fail on the fetched files")
    return {n: True for n in names}


def row_load(commit: str, repo: Path | None = None, guard=None, cfg: Mapping | None = None, reader=None) -> dict:
    """ROW T72'S LOAD STEP, at the row's commit: every guard of cell 8, data/external/t1_ext_t72/FETCH.json present,
    then data/external/t1_ext_t72/reduced/ created exclusively; per fetched member, the load of load_atlas (window
    removal first, the pins, the clearance, the screened copies' sha256, Gene Expression rows only, the cell call)
    with checks (a), (b) and (d) again on the fetched files, and every pinned sample loaded. A failure of any of these
    is a failed (p1) for that member: recorded in reduced/LOAD.json (the pre-gate reads it), its reduction not
    written, and the step goes on (exit 0). Returns LOAD.json's content."""
    rt = cm.sibling("row_t72")
    repo = Path(repo or cm.REPO)
    g = (guard or rt.guards)("load", commit, repo)
    cfg = cfg or cm.load_atlases(repo / rt.ATLASES_REL)
    ext = rt.row_external(repo)
    fp = ext / "FETCH.json"
    if not fp.is_file():
        raise Refused(f"{fp} is missing: the fetch step did not finish")
    fx = json.loads(fp.read_text())
    if fx.get("commit") != commit:
        raise Refused(f"{fp} was written at another commit")
    rt.make_exclusive(ext / "reduced")
    out = {**g, "fetch_sha256": cm.sha256_file(fp), "p1": {}, "dropped_by_screens": fx.get("dropped") or {},
           "load": {}}
    for acc in sorted(fx["atlases"]):
        try:
            a = cm.atlas(acc, cfg)
            rdir = cm.writable(ext / "raw" / acc)
            man = json.loads((rdir / "MANIFEST.json").read_text())
            if cm.sha256_file(rdir / "MANIFEST.json") != fx["atlases"][acc]["manifest_sha256"]:
                raise Refused(f"{acc}: the fetch manifest is not the one FETCH.json records")
            pins = cm.check_pins(acc, a)
            if set(man["samples"]) != set(pins["samples"]):
                raise Refused(f"{acc}: the fetched samples are not exactly the pinned ones")
            checks = recheck_screens(acc, a, cfg, man, rdir)
            red = load_atlas(acc, ext, cfg, check_frozen=False, reader=reader)  # the guards checked the freeze
            if red["qc"]["window"]["samples_removed_window"] or \
                    set(red["qc"]["window"]["samples_kept"]) != set(pins["samples"]):
                raise Refused(f"{acc}: not every pinned sample loaded")
            p = write_reduced(acc, ext, red, cfg)
            out["p1"][acc] = {"passed": True, "rechecked": checks}
            out["load"][acc] = {"manifest_sha256": cm.sha256_file(p), "qc": red["qc"]}
        except Refused as e:
            out["p1"][acc] = {"passed": False, "reason": str(e)}
    cm.dump_json(out, ext / "reduced" / "LOAD.json")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--commit", required=True, help="the row's commit, full sha (HEAD must be it)")
    a = ap.parse_args(argv)
    try:
        out = row_load(a.commit)
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    for acc, p1 in sorted(out["p1"].items()):
        if p1["passed"]:
            q = out["load"][acc]["qc"]
            per = {g: {k: x[k] for k in ("stage", "barcodes", "cells_called", "cells_kept")}
                   for g, x in q["per_sample"].items()}
            print(f"loaded\t{acc}\t(p1) passed\t{json.dumps(per)}\tgenes {q['genes_common']}")
        else:
            print(f"loaded\t{acc}\t(p1) FAILED\t{p1['reason']}")
    for acc, failed in out["dropped_by_screens"].items():
        print(f"loaded\t{acc}\tdropped by the screens ({failed}): not loaded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
