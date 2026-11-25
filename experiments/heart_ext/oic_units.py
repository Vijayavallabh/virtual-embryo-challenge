"""OIC125's units: the outside files of rows X021 and X022, binned on one 12 um grid (fetch_oic.py's constants).

- MOSTA (X022), one bin1 GEM file per section: the gene column, the DNB coordinates x and y and the raw UMI-count
  column (MIDCounts) only; DNB (x, y) goes to unit (x // 24, y // 24), 24 x 24 DNB = 12 um. No other column is read.
- Xenium (X021), the decoded-transcript table: gene name, x_location and y_location (um), the decoding quality
  (qv >= 20, the vendor's own threshold for a decoded transcript) and the gene flag (is_gene, or codeword_category,
  or else the control-name prefixes); (x, y) goes to unit (floor(x / 12), floor(y / 12)). Never cell_id,
  overlaps_nucleus, nucleus_distance or any other column: the vendor's segmentation enters nowhere.

A section's units: the integer grid indices (bx, by), the total count of every gene (the reader's plot), and the
counts of the board's panel genes (CSR, units x panel). The panel genes the source holds are those with any count.
Every file is checked against the fetch's MANIFEST.json (name, size, sha256; the sizes against fetch_oic.py's pins)
and its sha256 is recomputed from the bytes as they are read; a mismatch refuses the section.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import fetch_oic  # noqa: E402

GENE_COLS = ("geneID", "geneName", "gene", "gene_id")
COUNT_COLS = ("MIDCounts", "MIDCount", "UMICount", "UMICounts", "MIDcounts")
XENIUM_QV_MIN = 20.0
XENIUM_GENE_CATEGORIES = ("predesigned_gene", "custom_gene")
XENIUM_CONTROL_PREFIXES = ("NegControl", "BLANK", "Blank", "UnassignedCodeword", "DeprecatedCodeword", "Intergenic",
                           "antisense", "GenomicControl")
XENIUM_COLUMNS = ("feature_name", "x_location", "y_location", "qv")


class Refused(Exception):
    pass


# ---- the manifest -------------------------------------------------------------------------------------------------
def manifest_files(data_dir: Path) -> dict[str, dict]:
    """The fetch's MANIFEST.json, every record checked against fetch_oic.py's pinned list: name and URL pinned, a whole
    file's size the pinned size, the file on disk of the recorded size. Refused otherwise."""
    data_dir = Path(data_dir)
    path = data_dir / "MANIFEST.json"
    try:
        doc = json.loads(path.read_text())
    except (OSError, ValueError) as e:
        raise Refused(f"no readable MANIFEST.json in {data_dir}: {e}") from e
    pinned = {e["name"]: e for e in fetch_oic.PINNED}
    out = {}
    for r in doc.get("files", []):
        e = pinned.get(r.get("name"))
        if e is None or r.get("url") != e["url"]:
            raise Refused(f"MANIFEST.json records {r.get('name')!r} from {r.get('url')!r}: not a pinned file")
        if e["bytes"] is not None and r.get("size") != e["bytes"]:
            raise Refused(f"MANIFEST.json records {r['name']} at {r.get('size')} bytes, pinned {e['bytes']}")
        if e["kind"] == "zip_member" and r.get("archive_bytes") != e["archive_bytes"]:
            raise Refused(f"MANIFEST.json records {r['name']} from a bundle of another size")
        f = data_dir / r["name"]
        if not f.is_file() or f.stat().st_size != r.get("size"):
            raise Refused(f"{f} is missing or not of its recorded size")
        if not isinstance(r.get("sha256"), str) or len(r["sha256"]) != 64:
            raise Refused(f"MANIFEST.json gives no sha256 for {r['name']}")
        out[r["name"]] = {**r, "study": e["study"], "section": e["section"], "stage": e["stage"]}
    return out


class _Hashing(io.RawIOBase):
    """A read-through file that hashes every byte read (the section's sha256 is checked once it is read whole)."""

    def __init__(self, fh):
        super().__init__()
        self.fh, self.h, self.n = fh, hashlib.sha256(), 0

    def readable(self):
        return True

    def readinto(self, b):
        data = self.fh.read(len(b))
        n = len(data)
        b[:n] = data
        self.h.update(data)
        self.n += n
        return n


# ---- accumulation -------------------------------------------------------------------------------------------------
UNIT_OFF, UNIT_BITS = 2 ** 20, 22  # a unit's key: (bx + 2^20) * 2^22 + (by + 2^20)
GENE_BITS = 12  # a (unit, panel gene) key: unit key * 2^12 + gene index (a panel of at most 4096 genes)


def _reduce(keys: np.ndarray, vals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    uk, inv = np.unique(keys, return_inverse=True)
    return uk, np.bincount(inv.ravel(), weights=vals).astype(np.int64)


class Grid:
    """Counts summed into units: totals over every gene and panel counts, keyed by the unit's grid indices. Chunks are
    reduced as they come and merged every few chunks, so memory stays at the size of the section's units."""

    def __init__(self, panel: list[str]):
        self.panel = [str(g) for g in panel]
        if len(self.panel) >= 2 ** GENE_BITS:
            raise ValueError("a panel of more than 4095 genes")
        self.pidx = {g: i for i, g in enumerate(self.panel)}
        self.tot = [(np.zeros(0, np.int64), np.zeros(0, np.int64))]
        self.pan = [(np.zeros(0, np.int64), np.zeros(0, np.int64))]
        self.held = np.zeros(len(self.panel), bool)
        self.n_rows = 0

    @staticmethod
    def key(bx: np.ndarray, by: np.ndarray) -> np.ndarray:
        return (np.asarray(bx, np.int64) + UNIT_OFF) * 2 ** UNIT_BITS + (np.asarray(by, np.int64) + UNIT_OFF)

    def panel_index(self, names) -> np.ndarray:
        """Each name's panel index, -1 off the panel (through the names' categories: one lookup per distinct name)."""
        import pandas as pd

        cat = pd.Categorical(np.asarray(names).astype(str))
        look = np.array([self.pidx.get(str(c), -1) for c in cat.categories], np.int64)
        codes = np.asarray(cat.codes, np.int64)
        return np.where(codes >= 0, look[np.maximum(codes, 0)] if len(look) else -1, -1)

    def add_indexed(self, gi: np.ndarray, bx: np.ndarray, by: np.ndarray, counts: np.ndarray) -> None:
        counts = np.asarray(counts, np.int64)
        k = self.key(bx, by)
        self.n_rows += len(k)
        self.tot.append(_reduce(k, counts))
        gi = np.asarray(gi, np.int64)
        m = gi >= 0
        if m.any():
            self.held[np.unique(gi[m])] = True
            self.pan.append(_reduce(k[m] * 2 ** GENE_BITS + gi[m], counts[m]))
        if len(self.tot) > 16:
            self.tot = [_reduce(np.concatenate([a for a, _ in self.tot]), np.concatenate([b for _, b in self.tot]))]
        if len(self.pan) > 16:
            self.pan = [_reduce(np.concatenate([a for a, _ in self.pan]), np.concatenate([b for _, b in self.pan]))]

    def add(self, genes, bx: np.ndarray, by: np.ndarray, counts: np.ndarray) -> None:
        self.add_indexed(self.panel_index(genes), bx, by, counts)

    def finish(self) -> dict:
        tk, tv = _reduce(np.concatenate([a for a, _ in self.tot]), np.concatenate([b for _, b in self.tot]))
        pk, pv = _reduce(np.concatenate([a for a, _ in self.pan]), np.concatenate([b for _, b in self.pan]))
        unit, gene = pk // 2 ** GENE_BITS, pk % 2 ** GENE_BITS
        rows = np.searchsorted(tk, unit)
        if len(unit) and not np.array_equal(tk[np.minimum(rows, len(tk) - 1)], unit):
            raise RuntimeError("a panel count outside the units")
        M = sparse.csr_matrix((pv, (rows, gene)), shape=(len(tk), len(self.panel)), dtype=np.int64)
        M.sum_duplicates()
        bx = tk // 2 ** UNIT_BITS - UNIT_OFF
        by = tk % 2 ** UNIT_BITS - UNIT_OFF
        return {"bx": bx.astype(np.int64), "by": by.astype(np.int64), "total": tv, "panel": M,
                "genes": self.panel, "held": [g for g, h in zip(self.panel, self.held, strict=True) if h],
                "n_rows": self.n_rows}


# ---- MOSTA bin1 GEM -----------------------------------------------------------------------------------------------
def _gem_header(text) -> list[str]:
    while True:
        line = text.readline()
        if not line:
            raise Refused("a GEM file with no header line")
        if not line.startswith("#"):
            return line.rstrip("\r\n").split("\t")


def mosta_units(path: Path, panel: list[str], sha256: str | None = None, chunk_rows: int = 5_000_000,
                bin_dnb: int = fetch_oic.MOSTA_BIN_DNB) -> dict:
    """One MOSTA section's units at bin `bin_dnb` (the gene, x, y and count columns only)."""
    import pandas as pd

    path = Path(path)
    grid = Grid(panel)
    with open(path, "rb") as raw:
        hr = _Hashing(raw)
        with gzip.GzipFile(fileobj=io.BufferedReader(hr, 1 << 22)) as gz:
            text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
            cols = _gem_header(text)
            gene = next((c for c in GENE_COLS if c in cols), None)
            cnt = next((c for c in COUNT_COLS if c in cols), None)
            if gene is None or cnt is None or "x" not in cols or "y" not in cols:
                raise Refused(f"{path.name}: columns {cols} lack a gene, x, y or UMI-count column")
            reader = pd.read_csv(text, sep="\t", header=None, names=cols, usecols=[gene, "x", "y", cnt],
                                 chunksize=chunk_rows, dtype={gene: "category", "x": np.int64, "y": np.int64,
                                                              cnt: np.int64}, engine="c")
            for df in reader:
                cat = df[gene].cat
                look = np.array([grid.pidx.get(str(c), -1) for c in cat.categories], np.int64)
                codes = np.asarray(cat.codes, np.int64)
                gi = np.where(codes >= 0, look[np.maximum(codes, 0)], -1) if len(look) else np.full(len(df), -1)
                grid.add_indexed(gi, df["x"].to_numpy() // bin_dnb, df["y"].to_numpy() // bin_dnb,
                                 df[cnt].to_numpy())
            while gz.read(1 << 22):  # to the end, so the hash covers every byte
                pass
        buf = bytearray(1 << 22)
        while hr.readinto(buf):
            pass
    got = hr.h.hexdigest()
    if sha256 is not None and got != sha256:
        raise Refused(f"{path.name}: sha256 {got} is not the manifest's {sha256}")
    out = grid.finish()
    out.update(unit_um=fetch_oic.MOSTA_BIN_UM, sha256=got, source="X022")
    return out


# ---- Xenium decoded transcripts -----------------------------------------------------------------------------------
def xenium_gene_mask(names: np.ndarray, extra: dict) -> np.ndarray:
    """The gene flag: is_gene when the table has it, else codeword_category, else the control-name prefixes."""
    if "is_gene" in extra:
        return np.asarray(extra["is_gene"], bool)
    if "codeword_category" in extra:
        return np.isin(np.asarray(extra["codeword_category"]).astype(str), XENIUM_GENE_CATEGORIES)
    names = np.asarray(names).astype(str)
    return np.array([not n.startswith(XENIUM_CONTROL_PREFIXES) for n in names], bool)


def xenium_batches(path: Path, batch_rows: int = 5_000_000):
    """The transcript table's allowed columns, batch by batch (pyarrow; only these columns are ever read)."""
    import pyarrow.parquet as pq

    pf = pq.ParquetFile(path)
    names = set(pf.schema_arrow.names)
    missing = [c for c in XENIUM_COLUMNS if c not in names]
    if missing:
        raise Refused(f"{Path(path).name}: the transcript table lacks {missing}")
    extra = [c for c in ("is_gene", "codeword_category") if c in names][:1]
    for b in pf.iter_batches(batch_size=batch_rows, columns=[*XENIUM_COLUMNS, *extra]):
        yield {c: b.column(c).to_numpy(zero_copy_only=False) for c in b.schema.names}


def xenium_units(batches, panel: list[str], grid_um: float = fetch_oic.XENIUM_GRID_UM) -> dict:
    """The Xenium section's units on the 12 um grid, from batches of the allowed columns."""
    grid = Grid(panel)
    for d in batches:
        names = np.asarray(d["feature_name"], dtype=object)
        if len(names) and isinstance(names[0], bytes | bytearray):
            names = np.array([n.decode() for n in names], dtype=object)
        names = names.astype(str)
        extra = {k: d[k] for k in ("is_gene", "codeword_category") if k in d}
        keep = (np.asarray(d["qv"], np.float64) >= XENIUM_QV_MIN) & xenium_gene_mask(names, extra)
        x = np.floor(np.asarray(d["x_location"], np.float64)[keep] / grid_um).astype(np.int64)
        y = np.floor(np.asarray(d["y_location"], np.float64)[keep] / grid_um).astype(np.int64)
        grid.add(names[keep], x, y, np.ones(int(keep.sum()), np.int64))
    out = grid.finish()
    out.update(unit_um=grid_um, source="X021")
    return out


def xenium_file_units(path: Path, panel: list[str], sha256: str | None = None) -> dict:
    got = file_sha256(path)
    if sha256 is not None and got != sha256:
        raise Refused(f"{Path(path).name}: sha256 {got} is not the manifest's {sha256}")
    out = xenium_units(xenium_batches(path), panel)
    out["sha256"] = got
    return out


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(1 << 23)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


# ---- the units' own file ------------------------------------------------------------------------------------------
def save_units(u: dict, path: Path) -> str:
    """Writes the units as an .npz (written whole, then renamed) and returns its sha256."""
    path = Path(path)
    tmp = path.with_name(path.name + ".writing.npz")
    M = u["panel"].tocsr()
    np.savez_compressed(tmp, bx=u["bx"], by=u["by"], total=u["total"], data=M.data, indices=M.indices,
                        indptr=M.indptr, shape=np.array(M.shape), genes=np.array(u["genes"], dtype=object).astype(str),
                        held=np.array(u["held"], dtype=object).astype(str), unit_um=np.array(u["unit_um"]),
                        source_sha256=np.array(u.get("sha256", "")), source=np.array(u.get("source", "")))
    tmp.replace(path)
    return file_sha256(path)


def load_units(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    M = sparse.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    return {"bx": z["bx"], "by": z["by"], "total": z["total"], "panel": M, "genes": [str(g) for g in z["genes"]],
            "held": [str(g) for g in z["held"]], "unit_um": float(z["unit_um"]), "sha256": str(z["source_sha256"]),
            "source": str(z["source"])}
