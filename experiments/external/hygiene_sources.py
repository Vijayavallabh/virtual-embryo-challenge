#!/usr/bin/env python
"""R2's external gene lists (experiments/harness/TRIES.md R2; EXTERNAL_DATA.md X007 and X008), extracted by committed
code from the two hashed source files, so that no gene symbol is typed by hand (docs/PLAN.md 4.1, R2's provenance rule).

Each source file's size and sha256 are checked against configs/t1_hygiene.yaml (`sources`) before anything is read.

- van den Brink et al. 2017, Supplementary Information PDF (X007), turned into text by Poppler's `pdftotext -layout`
  (its version goes into the extracts' headers):
  - Supplementary Table 1, column 'Gene symbol' (module D): the lines holding a symbol, a log2 fold change and a
    p-value between the table's caption and the next table's caption. Their number must equal the count the caption
    states ("Only the N genes ..."). The rows are read whole to find and count them; only the symbols are written, the
    one column X007 registers.
  - Supplementary Table 5, column 'Gene Symbol' (read by N12 in place of D, deciding nothing): the one-word lines
    between the table's 'Gene Symbol' header and its footnote.
- Bakken et al. 2018, S6 Table (X008), an xlsx read with the standard library (zipfile and xml.etree). Module N2 is the
  rows with type_of_gene ncRNA, nuclear.prop >= 0.9 and expression in all 11 matched clusters in nuclei and in cells
  (inh.clusters_nuc + exc.clusters_nuc = 11, inh.clusters_cell + exc.clusters_cell = 11), thresholds from the config's
  `n2_filter`. 'NA' never passes a threshold.

Before anything is written, every list must have the config's registered number of rows (`extracts`: 32, 140, 19), no
symbol may repeat within a list, and Table 1's symbols must all be in Table 5. Each extract is a TSV under
experiments/external/extracts/ whose '#' header lines name the source (citation, URL, sha256, size), the extractor (this
file's sha256, and the tool with its version) and the date (UTC); the rows hold only the registered columns, as printed
or stored (D and Table 5: the gene symbol, X007; N2: the X008 columns). The script prints each extract's sha256, which
configs/t1_hygiene.yaml pins (`extracts`).

Requirements and reproducibility. Re-extracting needs Poppler's `pdftotext` (22.02.0 wrote the committed extracts), a
system program outside uv.lock that only this extraction uses; the xlsx needs only the standard library. The test of
record is `--check`: it re-extracts and compares the columns and rows with the committed extracts, writing nothing. An
extract's bytes also hold this file's sha256, the tool's version (Python's for N2) and the date, so a rerun writes the
same bytes, and so the pinned sha256, only with `--date` set to the extracts' date and the same extractor, Poppler and
Python versions.

Usage (repo root; the downloads stay outside the repository):
  .venv/bin/python -I -B experiments/external/hygiene_sources.py --vdb SI.pdf --bakken S6.xlsx [--date YYYY-MM-DD]
  .venv/bin/python -I -B experiments/external/hygiene_sources.py --vdb SI.pdf --bakken S6.xlsx --check
"""

from __future__ import annotations

import argparse
import hashlib
import math
import platform
import re
import subprocess
import sys
import zipfile
from collections.abc import Mapping, Sequence
from datetime import UTC, date, datetime
from pathlib import Path
from xml.etree import ElementTree as ET

import yaml

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "configs" / "t1_hygiene.yaml"
THIS = Path(__file__).resolve()

TABLE1 = (r"\s*Supplementary table 1\. ", r"\s*Supplementary table 2\. ")  # captions that open and close each table
TABLE5 = (r"\s*Supplementary table 5\. ", r"\s*Supplementary table 6\. ")
TABLE1_ROW = re.compile(r"\s*([A-Za-z][A-Za-z0-9.\-]*)\s+(-?\d+\.\d+)\s+(\d(?:\.\d+)?E[-+]\d+)\s*")
TABLE1_STATED = re.compile(r"Only the (\d+) genes")
TABLE5_HEADER = "Gene Symbol"
SYMBOL = re.compile(r"[A-Za-z][A-Za-z0-9.\-]*")
N2_COLUMNS = ("gene", "type_of_gene", "inh.clusters_nuc", "exc.clusters_nuc", "inh.clusters_cell",
              "exc.clusters_cell", "nuclear.prop")
MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DOC_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(1 << 22):
            h.update(block)
    return h.hexdigest()


def check_source(path: Path, src: Mapping) -> str:
    """The file's sha256, after its size and sha256 are checked against the config; refuses any other file."""
    size, sha = path.stat().st_size, sha256_file(path)
    if size != int(src["bytes"]) or sha != src["sha256"]:
        sys.exit(f"{path}: {size} bytes, sha256 {sha}; configs/t1_hygiene.yaml pins {src['bytes']} bytes, sha256 "
                 f"{src['sha256']} (refused: not the registered file)")
    return sha


# ---------------------------------------------------------------------------------------------------------------------
# van den Brink 2017 (X007): the supplementary PDF through pdftotext


def poppler_version() -> str:
    """pdftotext's version line, e.g. 'pdftotext version 22.02.0' (Poppler)."""
    r = subprocess.run(["pdftotext", "-v"], capture_output=True, text=True)
    line = next((ln.strip() for ln in (r.stderr + r.stdout).splitlines() if ln.startswith("pdftotext version")), None)
    if line is None:
        sys.exit("pdftotext -v printed no version line")
    return line


def pdf_text(path: Path) -> str:
    """The PDF's text in its physical layout (pdftotext -layout, UTF-8)."""
    r = subprocess.run(["pdftotext", "-layout", "-enc", "UTF-8", str(path), "-"], capture_output=True, check=True)
    return r.stdout.decode("utf-8")


def _block(lines: Sequence[str], opens: str, closes: str) -> list[str]:
    """The lines from the one that matches `opens` (at its start) up to, not including, the next one matching
    `closes`; each caption must occur exactly once."""
    starts = [i for i, ln in enumerate(lines) if re.match(opens, ln)]
    ends = [i for i, ln in enumerate(lines) if re.match(closes, ln)]
    if len(starts) != 1 or len(ends) != 1 or ends[0] <= starts[0]:
        sys.exit(f"the captions {opens!r} and {closes!r} are not found once each, in order")
    return list(lines[starts[0]:ends[0]])


def vdb_table1(text: str) -> list[tuple[str, str, str]]:
    """Supplementary Table 1: (gene symbol, log2 fold change, p-value) per row, in the table's order."""
    block = _block(text.splitlines(), *TABLE1)
    rows = [m.groups() for ln in block if (m := TABLE1_ROW.fullmatch(ln))]
    stated = TABLE1_STATED.search(" ".join(block))
    if stated is None or int(stated.group(1)) != len(rows):
        sys.exit(f"Supplementary Table 1: {len(rows)} rows read, its caption states "
                 f"{stated.group(1) if stated else 'no count'}")
    return rows


def vdb_table5(text: str) -> list[str]:
    """Supplementary Table 5: the gene symbols, in the table's order (one word per line under 'Gene Symbol' until the
    footnote; page footers, which hold several words, are skipped)."""
    block = _block(text.splitlines(), *TABLE5)
    heads = [i for i, ln in enumerate(block) if ln.strip() == TABLE5_HEADER]
    if len(heads) != 1:
        sys.exit(f"Supplementary Table 5: {len(heads)} '{TABLE5_HEADER}' headers")
    out = []
    for ln in block[heads[0] + 1:]:
        s = ln.strip()
        if s.startswith("*"):  # the footnote ends the table
            return out
        if SYMBOL.fullmatch(s):
            out.append(s)
    sys.exit("Supplementary Table 5: no footnote after the symbols")


# ---------------------------------------------------------------------------------------------------------------------
# Bakken 2018 (X008): the S6 xlsx with the standard library


def _column(ref: str | None) -> str:
    m = re.match(r"[A-Z]+", ref or "")
    if m is None:
        sys.exit(f"an xlsx cell without a column reference ({ref!r})")
    return m.group(0)


def _shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    out = []
    with z.open("xl/sharedStrings.xml") as f:
        for _, el in ET.iterparse(f):
            if el.tag == f"{MAIN}si":  # plain text, or rich-text runs (phonetic runs are left out)
                parts = [el.findtext(f"{MAIN}t") or ""]
                parts += [r.findtext(f"{MAIN}t") or "" for r in el.findall(f"{MAIN}r")]
                out.append("".join(parts))
                el.clear()
    return out


def _cell_text(c: ET.Element, strings: Sequence[str]) -> str:
    kind = c.get("t")
    if kind == "inlineStr":
        return "".join(t.text or "" for t in c.iter(f"{MAIN}t"))
    v = c.findtext(f"{MAIN}v")
    if v is None:
        return ""
    return strings[int(v)] if kind == "s" else v  # numbers as stored


def xlsx_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """The workbook's first sheet: the first row's names, and the other rows keyed by them (cell text, shared strings
    resolved; an empty cell is absent)."""
    with zipfile.ZipFile(path) as z:
        sheet = ET.fromstring(z.read("xl/workbook.xml")).find(f"{MAIN}sheets/{MAIN}sheet")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in rels if r.get("Id") == sheet.get(f"{DOC_REL}id"))
        part = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
        strings = _shared_strings(z)
        rows = []
        with z.open(part) as f:
            for _, el in ET.iterparse(f):
                if el.tag == f"{MAIN}row":
                    rows.append({_column(c.get("r")): _cell_text(c, strings) for c in el.findall(f"{MAIN}c")})
                    el.clear()
    if not rows:
        sys.exit(f"{path}: the first sheet is empty")
    names = rows[0]
    if len(set(names.values())) != len(names):
        sys.exit(f"{path}: repeated column names in the first row")
    return list(names.values()), [{names[col]: v for col, v in r.items() if col in names} for r in rows[1:]]


def _number(text: str) -> float:
    try:
        return float(text)
    except ValueError:  # 'NA' or an empty cell: never passes a threshold
        return math.nan


def bakken_n2(columns: Sequence[str], rows: Sequence[Mapping[str, str]], f: Mapping) -> list[list[str]]:
    """Module N2: the rows passing `f` (the config's n2_filter), with N2_COLUMNS as stored, in the sheet's order."""
    missing = [c for c in N2_COLUMNS if c not in columns]
    if missing:
        sys.exit(f"S6 Table: columns {missing} are missing")

    def num(r, col):
        return _number(r.get(col, ""))

    keep = [r for r in rows
            if r.get("type_of_gene") == f["type_of_gene"]
            and num(r, "nuclear.prop") >= float(f["min_nuclear_prop"])
            and num(r, "inh.clusters_nuc") + num(r, "exc.clusters_nuc") == float(f["clusters_nuc"])
            and num(r, "inh.clusters_cell") + num(r, "exc.clusters_cell") == float(f["clusters_cell"])]
    return [[r.get(c, "") for c in N2_COLUMNS] for r in keep]


# ---------------------------------------------------------------------------------------------------------------------
# Extract files


def extract_text(meta: Sequence[tuple[str, str]], columns: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    """'# key: value' header lines, then a tab-separated table."""
    for v in [*columns, *(x for r in rows for x in r)]:
        if any(ch in v for ch in "\t\r\n"):
            sys.exit(f"a value holds a tab or a line break: {v!r}")
    lines = [f"# {k}: {v}" for k, v in meta] + ["\t".join(columns)] + ["\t".join(r) for r in rows]
    return "\n".join(lines) + "\n"


def table_rows(text: str) -> tuple[list[str], list[list[str]]]:
    """An extract's column names and rows (its '#' lines left out)."""
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]


def extracts(vdb: Path, bakken: Path, cfg: Mapping) -> dict[str, dict]:
    """The three lists with their columns, rows and header, after every check; nothing is written here."""
    srcs = cfg["sources"]
    vdb_sha, bak_sha = check_source(vdb, srcs["vdb2017_si"]), check_source(bakken, srcs["bakken2018_s6"])
    text, poppler = pdf_text(vdb), poppler_version()
    t1, t5 = vdb_table1(text), vdb_table5(text)
    n2 = bakken_n2(*xlsx_rows(bakken), cfg["n2_filter"])
    lists = {
        "D": {"what": "R2 module D: van den Brink 2017 Supplementary Table 1, column 'Gene symbol'",
              "columns": ["gene_symbol"], "rows": [[r[0]] for r in t1],  # X007 registers the symbols only
              "tool": f"pdftotext -layout -enc UTF-8 ({poppler}, Poppler)", "sha": vdb_sha},
        "D_table5": {"what": "van den Brink 2017 Supplementary Table 5, column 'Gene Symbol' (N12 reads it in place "
                             "of D; it decides nothing)",
                     "columns": ["gene_symbol"], "rows": [[s] for s in t5],
                     "tool": f"pdftotext -layout -enc UTF-8 ({poppler}, Poppler)", "sha": vdb_sha},
        "N2": {"what": "R2 module N2: Bakken 2018 S6 Table rows with type_of_gene {type_of_gene}, nuclear.prop >= "
                       "{min_nuclear_prop}, inh.clusters_nuc + exc.clusters_nuc = {clusters_nuc} and inh.clusters_cell "
                       "+ exc.clusters_cell = {clusters_cell}".format(**cfg["n2_filter"]),
               "columns": ["gene_symbol", *N2_COLUMNS[1:]], "rows": n2,
               "tool": f"Python {platform.python_version()} zipfile and xml.etree.ElementTree (standard library)",
               "sha": bak_sha},
    }
    for name, ent in lists.items():
        symbols = [r[0] for r in ent["rows"]]
        want = int(cfg["extracts"][name]["n"])
        if len(symbols) != want or len(set(symbols)) != len(symbols):
            sys.exit(f"{name}: {len(symbols)} rows ({len(set(symbols))} distinct symbols); the config registers {want}")
    missing = {r[0] for r in lists["D"]["rows"]} - {r[0] for r in lists["D_table5"]["rows"]}
    if missing:
        sys.exit(f"{len(missing)} of Table 1's symbols are not in Table 5")
    return lists


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vdb", type=Path, required=True, help="van den Brink 2017 Supplementary Information PDF (X007)")
    ap.add_argument("--bakken", type=Path, required=True, help="Bakken 2018 S6 Table xlsx (X008)")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--check", action="store_true", help="compare with the committed extracts; write nothing")
    ap.add_argument("--date", type=date.fromisoformat,
                    help="the date the headers record, YYYY-MM-DD (default: today, UTC); the extracts' date writes "
                         "their bytes again")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(a.config.read_text())
    lists = extracts(a.vdb, a.bakken, cfg)
    if a.check:
        bad = []
        for name, ent in lists.items():
            path = REPO / cfg["extracts"][name]["path"]
            same = path.exists() and table_rows(path.read_text(encoding="utf-8")) == (ent["columns"], ent["rows"])
            print(f"{name}: {'same rows as' if same else 'DIFFERS from'} {cfg['extracts'][name]['path']}")
            bad += [] if same else [name]
        sys.exit(1 if bad else 0)
    day = (a.date or datetime.now(UTC).date()).isoformat()
    extractor = f"{THIS.relative_to(REPO)} (sha256 {sha256_file(THIS)})"
    for name, ent in lists.items():
        dst = cfg["extracts"][name]
        src = cfg["sources"][dst["source"]]
        meta = [("extract", ent["what"]),
                ("register", f"EXTERNAL_DATA.md {src['register']}; experiments/harness/TRIES.md R2"),
                ("source", " ".join(str(src["citation"]).split())), ("source_url", src["url"]),
                ("source_sha256", ent["sha"]), ("source_bytes", str(src["bytes"])),
                ("extractor", extractor), ("tool", ent["tool"]), ("date", day), ("rows", str(len(ent["rows"])))]
        data = extract_text(meta, ent["columns"], ent["rows"]).encode("utf-8")
        path = REPO / dst["path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        print(f"{name}: {len(ent['rows'])} rows -> {dst['path']} sha256 {hashlib.sha256(data).hexdigest()}")
    t1 = {r[0] for r in lists["D"]["rows"]}
    t5 = {r[0] for r in lists["D_table5"]["rows"]}
    print(f"Table 1's {len(t1)} symbols are all in Table 5's {len(t5)}: {t1 <= t5}")


if __name__ == "__main__":
    main()
