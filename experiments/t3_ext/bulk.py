#!/usr/bin/env python
"""Bulk whole-embryo knockout units (stage 9, campaign/t3-ext-bulk; README "Bulk whole-embryo units").

A bulk unit is one knockout gene of one study at one pinned stage: its homozygous or null embryos against wild-type
embryos of the same series, line and stage, one library (one embryo, or one pool of one genotype) per replicate. The
counts come from GEO through two forms, or (phase 2, ae.py) from one ArrayExpress resource:

  per_sample    one plain text count table in each pinned sample's OWN GEO folder (the single-cell own-folder rule,
                DENY and sample_file_patterns apply; .txt, .tsv, .csv or .tab, optionally gzipped; XLSX, RDS and
                binary forms are refused). It holds the declared gene_columns and exactly one count column; with more
                columns, count_column names the header, which must tie to this sample by the exact tie rule. Its
                first line is screened before any value is parsed (a held-out hit refuses with the text withheld; a
                review term needs a review_cleared entry naming the sample), and a single count column that ties by
                the exact rule to another sample is refused (round 2).
  series_table  ONE text count table of a screened series' OWN folder. Before any pin, the header step
                (`fetch.py header`) streams the file, stops at the first newline (at most HEADER_MAX_BYTES of
                decompressed text), writes only that line and its sha256, screens every header and the whole line,
                and prints the tie table: every header is a declared gene column or ties to exactly one sample of the
                screened series by EXACT equality with its accession, its whole title, a whole description value or
                the value part of one characteristics entry (the text after the first ': '), and no sample ties to
                two headers; otherwise the whole table is refused. Each pinned sample names its column; the table is
                pinned in pinned_series_files {series, name, header_sha256} and fetched only after the pins; process
                refuses it unless its first line hashes to header_sha256 and reads it with usecols (the gene columns
                and the pinned columns only; the other columns' bytes pass through memory, never parsed).
  arrayexpress  ae.py (phase 2).

The math (bulk_log_ratio, pure numpy): per library i, x_gi its counts (runs or lanes of one library summed) and T_i its
total over ALL gene rows of its file or column ('__*' and 'N_*' rows left out); y_gi = log((x_gi + 10) / T_i);
L_g = mean over knockout libraries of y_gi minus the mean over control libraries, clipped to log [0.25, 4]. Not
measured (NaN): the median of x_gi under 20 in either arm; a gene absent from any pinned library; the knocked-out
gene's own row; Xist, Tsix, the Y genes of the sex call and every chrY panel gene (always); the other chrX panel genes
when the sex check finds an imbalance or cannot run. k = (m/(m+1))^2 * d (t3x_common.shrink_k_bulk), m the smaller
number of libraries left in an arm after the guards; stages are averaged with equal weight. Counts only: finite,
non-negative and integral within 1e-6; a table whose columns all sum to within 1% of 1e6, a second count column
without count_column, and a column or file named as FPKM, RPKM, TPM, CPM, RSEM, normalised or expected counts are
refused (a series or matrix table so named at the offer, the header step and the read). A GEO sample whose own title,
source name or characteristics carry heterozygous wording is never pinned, in either arm (round 2).

The guards (fail closed; every number fixed by the reviewed design before any count; outside libraries and the WT
E8.75 reference statistics only; verdicts to QC.json; no gene-level value is printed): G1 library size, G2 counts
only, G3 sex, G4 outlier and grouping, G5 whole embryo, G6 somites, G7 stage and design (the pins), G8 duplicates,
G9 the genotype in the counts, G10 every single-cell refusal kept, G12 composition.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import sys
import time
import urllib.parse
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import fetch as fe  # noqa: E402
import t3x_common as cm  # noqa: E402
import t3x_screen as sc  # noqa: E402
import t3x_tables as tb  # noqa: E402

Refused = cm.Refused
WHOLE = "__whole__"
NOT_TEXT = "not a plain text table (XLSX, RDS, h5 and other binary forms are refused this round)"
NORMALISED = re.compile(r"fpkm|rpkm|tpm|(?<![a-z])cpm|rsem|normali[sz]|expected[ _.-]?counts?", re.IGNORECASE)
HEADER_MAX_COMPRESSED = 16 * 1024 ** 2   # the header step reads at most this many bytes off the wire
READ_BLOCK = 4096   # round 2 (the first reviewer's lesser point 1): before the pins, at most one block of this size
#                     passes through memory past the first newline; it is never parsed or written
# round 2 (the first reviewer's blocker 5): heterozygous wording in a GEO sample's own title, source name or
# characteristics refuses a pin in either arm (a false hit costs one sample)
HET_WORDING = re.compile(r"(?<![A-Za-z])het(?:erozygo\w*)?(?![A-Za-z])|\+/-|-/\+", re.IGNORECASE)
VERSION = re.compile(r"\.\d+$")
SPECIAL_ROW = ("__", "N_")


def refuse_normalised_name(name: str) -> None:
    """Round 2 (the first reviewer's blocker 4): a count table NAMED as normalised or as expected counts is refused at
    the offer, at the header step and at the read, for the series form and the accession-keyed form alike."""
    if NORMALISED.search(name):
        raise Refused(f"{name}: named as normalised or expected counts: refused (counts only)")


def het_wording(rec: dict) -> str | None:
    """The reason a GEO sample's own record carries heterozygous wording (its title, source name or a characteristics
    value), or None."""
    for k in sorted(rec):
        if k == "title" or k.startswith(("source_name", "characteristics")):
            if any(HET_WORDING.search(v) for v in rec[k]):
                return f"heterozygous wording in the sample's own record ({k}): never pinned this round"
    return None


# ---------------------------------------------------------------------------------------------------------------------
# GEO: the files a bulk dataset may take (derived from the hashed fields; verify_pins derives them again)


def _links(rec: dict) -> list[str]:
    return [v for k in sorted(rec) if k.startswith("supplementary_file") for v in rec[k]]


def sample_files(ds: dict, rec: dict, gsm: str) -> tuple[list[dict], list[dict]]:
    """(the per-sample count table offered, as [{name, url}]; every other link as {name, reason}). per_sample only:
    a file of the sample's own folder (fetch.own_file), naming no held accession, not denied, a plain text table
    (BULK_TEXT) matched by sample_file_patterns. A sample with more than one such table offers none of them."""
    pats = [re.compile(p, re.IGNORECASE) for p in ds.get("sample_file_patterns") or []]
    cf = ds.get("count_form")
    ok, other = [], []
    for f in _links(rec):
        if f.upper() == "NONE":
            continue
        name = Path(urllib.parse.urlparse(f).path).name
        named = fe._named_held(ds, f)
        u, why = fe.own_file(gsm, f)
        if named:
            why = f"url names a held or excluded accession {sorted(named)}"
        elif u and cf != "per_sample":
            why = f"count_form {cf}: the counts come from the pinned count table, never a sample file"
        elif u and fe.denied(name):
            why = "a denied file type"
        elif u and not cm.BULK_TEXT.search(name):
            why = NOT_TEXT
        elif u and not any(p.search(name) for p in pats):
            why = "not matched by the dataset's sample_file_patterns"
        if why:
            other.append({"name": name, "reason": why})
        else:
            ok.append({"name": name, "url": u})
    if len(ok) > 1:
        names = [f["name"] for f in ok]
        other += [{"name": n, "reason": f"the sample lists {len(ok)} text tables of its own ({', '.join(names)}): one "
                   "text table per sample is taken, so none is"} for n in names]
        ok = []
    return ok, other


def series_files(ds: dict, gse: str, fields: dict) -> tuple[list[dict], list[dict]]:
    """(the series count tables offered, as [{name, url}]; every other link as {name, reason}). series_table only: a
    file of the series' own folder (fetch.own_series_file), naming no held accession, not denied, a plain text table."""
    cf = ds.get("count_form")
    ok, other = [], []
    for f in [v for k in fields if k.startswith("supplementary_file") for v in fields[k]]:
        if f.upper() == "NONE":
            continue
        name = Path(urllib.parse.urlparse(f).path).name
        named = fe._named_held(ds, f)
        u, why = fe.own_series_file(gse, f)
        if named:
            why = f"url names a held or excluded accession {sorted(named)}"
        elif u and cf != "series_table":
            why = f"count_form {cf}: no series file is taken"
        elif u and fe.DENY.search(name):
            why = "a denied file type"
        elif u and not cm.BULK_TEXT.search(name):
            why = NOT_TEXT
        elif u and NORMALISED.search(name):
            why = "named as normalised or expected counts"
        if why:
            other.append({"name": name, "reason": why})
        else:
            ok.append({"name": name, "url": u})
    return ok, other


# ---------------------------------------------------------------------------------------------------------------------
# Stages: stage_from_series, and the superseries rule


def series_stage(fields: dict, field: str) -> tuple[float | None, str]:
    """(E-day, how) from one series field: every value of it read by the stage parser; one value exactly, no range,
    at or before E9.5; else None with the reason."""
    vals = fields.get(field, [])
    stated = set()
    for v in vals:
        if fe._STAGE_RANGE.search(v):
            return None, f"range in the series field {field}"
        stated |= {float(x) for x in fe._STAGE_VAL.findall(v)}
    if len(stated) != 1:
        return None, f"the series field {field} states {len(stated)} stage values"
    e = stated.pop()
    try:
        cm.stage_ok(f"E{e:g}")
    except Refused as err:
        return None, f"{err} (series field {field})"
    return e, f"series field {field}"


def row_stage(ds: dict, rec: dict, series_fields: dict | None) -> tuple[float | None, str]:
    """A bulk sample's stage: under stage_from_series the own series' field gives it, and a sample whose own fields
    state another stage (or a range, or two values) is refused; else the single-cell rule (fetch.sample_stage)."""
    f = ds.get("stage_from_series")
    if not f:
        return fe.sample_stage(rec, cm.stage_fields(ds))
    e, how = series_stage(series_fields or {}, f)
    if e is None:
        return None, how
    es, hows = fe.sample_stage(rec, cm.stage_fields(ds))
    if es is None and not (hows == "0 stage values" or hows.startswith("no listed stage field")):
        return None, f"stage_from_series gives E{e:g} but the sample's own fields give {hows}"
    if es is not None and es != e:
        return None, f"stage_from_series gives E{e:g} but the sample's own fields give E{es:g}"
    return e, how


def held_relations(ds: dict, fields: dict) -> list[str]:
    """The held series a series record's relation values name ('SuperSeries of: X', 'SubSeries of: X', any held
    accession named in a relation)."""
    held = fe.held_set(ds)
    out = set()
    for v in fields.get("relation", []):
        out |= {m.group(0).upper() for m in fe._ACCESSION.finditer(v)} & {h.upper() for h in held}
    return sorted(out)


def adjust_plan(ds: dict, out: dict, metas: dict) -> None:
    """The bulk part of a GEO plan (fetch.plan_from_meta calls it for a bulk dataset): the declaration the plan is
    bound to, and every sample's stage under stage_from_series."""
    acc = ds["accession"]
    out["bulk"] = declaration(ds)
    sf = metas.get(acc, {}).get("series", {})
    for rows in out["tables"].values():
        for r in rows:
            r["stage"], r["stage_from"] = row_stage(ds, r["fields"], sf)
            het = het_wording(r["fields"])
            if het:
                r["not_pinnable"] = "; ".join(x for x in (r.get("not_pinnable"), het) if x)


def declaration(ds: dict) -> dict:
    keys = ("modality", "count_form", "id_type", "gene_columns", "count_column", "stage_from_series", "somite_field",
            "sex_field", "withheld_lines", "source_name_pattern")
    return {k: ds[k] for k in keys if k in ds}


def plan_lines(pl: dict) -> list[str]:
    b = pl.get("bulk")
    if not b:
        return []
    return [f"# BULK (modality bulk): count_form {b.get('count_form')}, id_type {b.get('id_type')}, gene_columns "
            f"{b.get('gene_columns')}, count_column {b.get('count_column')}, stage_from_series "
            f"{b.get('stage_from_series')}; sex field {b.get('sex_field')}, somite field {b.get('somite_field')}",
            "#   a bulk unit: homozygous or null embryos against wild type of the same series, line and stage; hets, "
            "point, ENU, catalytic, hypomorph and isoform alleles are never pinned (README check B1)",
            f"#   at least {cm.BULK_MIN_LIBS} libraries per arm per stage; a series count table is read only after the "
            "header step, its tie table and the pins"]


# ---------------------------------------------------------------------------------------------------------------------
# The header step and the exact tie rule


def first_line(read, gz: bool, label: str) -> bytes:
    """The bytes before the first newline of a stream (read(n) -> bytes), decompressed incrementally when gz, one
    trailing CR removed. Refused past HEADER_MAX_BYTES of decompressed text, past HEADER_MAX_COMPRESSED read, or at
    the end of the stream without a newline."""
    d = zlib.decompressobj(16 + zlib.MAX_WBITS) if gz else None
    buf = bytearray()
    n_in = 0
    eof = False
    while True:
        if d is not None and d.unconsumed_tail:
            chunk = d.decompress(d.unconsumed_tail, READ_BLOCK)
        else:
            block = read(READ_BLOCK)
            if not block:
                eof = True
                chunk = d.flush() if d is not None else b""
            else:
                n_in += len(block)
                if n_in > HEADER_MAX_COMPRESSED:
                    raise Refused(f"{label}: no newline within {HEADER_MAX_COMPRESSED} bytes read: refused")
                try:
                    chunk = d.decompress(block, READ_BLOCK) if d is not None else block
                except zlib.error:
                    raise Refused(f"{label}: not a gzip stream: refused") from None
        buf += chunk
        i = buf.find(b"\n")
        if i >= 0:
            if i > cm.HEADER_MAX_BYTES:
                break
            line = bytes(buf[:i])
            return line[:-1] if line.endswith(b"\r") else line
        if len(buf) > cm.HEADER_MAX_BYTES:
            break
        if eof:
            raise Refused(f"{label}: the file ends before its first newline: refused")
    raise Refused(f"{label}: a first line longer than {cm.HEADER_MAX_BYTES} bytes: refused")


def header_sha256(line: bytes) -> str:
    return hashlib.sha256(line).hexdigest()


def sep_of(name: str) -> str:
    low = name.lower().removesuffix(".gz")
    return "," if low.endswith(".csv") else "\t"


def split_header(line: bytes, name: str) -> list[str]:
    """The header fields of a first line: UTF-8 (strict), split on the file's separator (comma for .csv, else tab);
    a field wrapped in one pair of double quotes is unquoted; any other quote character refuses."""
    try:
        text = line.decode("utf-8")
    except UnicodeDecodeError:
        raise Refused(f"{name}: the first line is not UTF-8: refused") from None
    out = []
    for f in text.split(sep_of(name)):
        if len(f) >= 2 and f[0] == '"' and f[-1] == '"' and '"' not in f[1:-1]:
            f = f[1:-1]
        if '"' in f:
            raise Refused(f"{name}: a header field with a quote character inside: refused")
        out.append(f)
    return out


def tie_candidates(gsm: str, rec: dict) -> dict[str, str]:
    """The exact strings a header may equal to tie to this sample: its accession, a whole title value, a whole
    description value, the value part (after the first ': ') of one characteristics entry."""
    c = {gsm: "accession"}
    for v in rec.get("title", []):
        c.setdefault(v, "title")
    for v in rec.get("description", []):
        c.setdefault(v, "description")
    for v in rec.get("characteristics_ch1", []):
        if ": " in v:
            c.setdefault(v.split(": ", 1)[1], "characteristics value")
    return c


def tie_table(fields: list[str], rows: list[dict], gene_columns: list[str]) -> tuple[list[dict], list[str]]:
    """The exact tie rule (no tokens): ([{header, sample, how} per header in file order], problems). A header that is
    a declared gene column ties to none; every other header must equal exactly one sample's candidate (tie_candidates);
    no sample may tie to two headers; every gene column appears exactly once; no header appears twice. Readers
    confirm the table; they do not decide ties."""
    samples: dict[str, dict] = {}
    for r in rows:
        samples.setdefault(r["accession"], r["fields"])
    cands = {g: tie_candidates(g, rec) for g, rec in samples.items()}
    problems, out, tied = [], [], {}
    dup = sorted({h for h in fields if fields.count(h) > 1})
    if dup:
        problems.append(f"headers appear twice: {dup}")
    for gcol in gene_columns:
        if fields.count(gcol) != 1:
            problems.append(f"the gene column {gcol!r} appears {fields.count(gcol)} times (headerless, or not the "
                            "declared column)")
    for i, h in enumerate(fields):
        if h in gene_columns:
            out.append({"index": i, "header": h, "sample": None, "how": "gene column"})
            continue
        hits = sorted(g for g, c in cands.items() if h in c)
        if len(hits) != 1:
            problems.append(f"header {i} ties to {len(hits)} samples ({hits})")
            out.append({"index": i, "header": h, "sample": None, "how": f"ties to {len(hits)} samples"})
            continue
        g = hits[0]
        out.append({"index": i, "header": h, "sample": g, "how": cands[g][h]})
        tied.setdefault(g, []).append(h)
    for g, hs in tied.items():
        if len(hs) > 1:
            problems.append(f"{g} ties to {len(hs)} headers")
    return out, problems


def screen_header(ds: dict, series: str, fields: list[str], text: str) -> tuple[bool, list[str]]:
    """(held, uncleared hits) of the header: every field and the whole line through t3x_screen, review terms cleared
    only by a review_cleared entry naming the series."""
    held, bad = False, []
    extra, cleared = ds.get("screen_extra_refuse"), ds.get("review_cleared")
    void = any(sc.extended(t) for t in [text, *fields])
    for t in [text, *fields]:
        for kind, term in sc.screen_text(series, t, extra, cleared, void=void):
            held |= kind == "held_out"
            bad.append(f"{kind}" + (f" term {term!r}" if kind == "review" else ""))
    return held, sorted(set(bad))


def header_path(dest: Path, acc: str, name: str) -> Path:
    return cm.outside_repo(dest) / "meta" / acc / f"header_{name}.json"


def record_header(ds: dict, dest: Path, series: str, name: str, url: str, line: bytes, rows: list[dict],
                  series_sha: str) -> dict:
    """The header record: written to meta/<ACC>/header_<name>.json (and the line to header_<name>.txt). A held-out hit
    withholds the line: only the number of columns is written. A review term not cleared, or a tie problem, refuses
    the table (status refused)."""
    acc = ds["accession"]
    refuse_normalised_name(name)
    fields = split_header(line, name)
    text = line.decode("utf-8")
    rec = {"accession": acc, "series": series, "name": name, "url": url, "n_columns": len(fields),
           "header_sha256": header_sha256(line), "series_metadata_sha256": series_sha,
           "gene_columns": ds["gene_columns"], "created": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
    held, bad = screen_header(ds, series, fields, text)
    hp = header_path(dest, acc, name)
    tp = hp.with_suffix(".txt")
    tp.unlink(missing_ok=True)
    if held:
        rec.update(status="refused", withheld=True, reason="a held-out hit in the header (text withheld)")
        cm.dump_json(rec, hp)
        return rec
    rec.update(withheld=False, text=text, fields=fields, screen=bad)
    ties, problems = tie_table(fields, rows, ds["gene_columns"])
    rec["ties"] = ties
    rec["tie_problems"] = problems
    if bad:
        rec.update(status="refused", reason=f"header screen: {bad}")
    elif problems:
        rec.update(status="refused", reason="tie rule: " + "; ".join(problems[:10]))
    else:
        rec["status"] = "ok"
    cm.dump_json(rec, hp)
    tp.write_text(text + "\n", encoding="utf-8")
    return rec


def header(acc: str, dest: Path, series: str | None = None, name: str | None = None,
           datasets: Path | None = None) -> dict:
    """The header step of a series_table dataset (an arrayexpress one goes to ae.header): from the stored plan (status
    ok, the committed datasets.yaml), the named file must be one the series' own hashed fields offer; its first line is
    streamed (fetch._open), screened and tied, and only the line and its record are written."""
    ds = cm.dataset(acc, datasets)
    if not cm.is_bulk(ds) or ds.get("count_form") not in ("series_table", "arrayexpress"):
        raise Refused(f"{acc}: the header step is for a bulk series_table or arrayexpress dataset")
    if ds["count_form"] == "arrayexpress":
        import ae

        return ae.header(acc, dest, datasets)
    dest = cm.outside_repo(dest)
    pl = fe._load_plan(acc, dest, datasets)
    if pl.get("withheld") or series not in (pl.get("series") or {}):
        raise Refused(f"{acc}: {series} is not a screened series of the stored plan")
    b = pl["series"][series]
    if cm.series_metadata_sha256(series, b["fields"]) != b["series_metadata_sha256"]:
        raise Refused(f"{series}: the stored series fields no longer hash to the plan's value: plan again")
    refuse_normalised_name(name or "")
    derived, _ = series_files(ds, series, b["fields"])
    if b.get("files") != derived:
        raise Refused(f"{series}: the stored series file list differs from the one its fields give: plan again")
    offered = {f["name"]: f for f in derived}
    if name not in offered or fe.duplicate_names(derived):
        raise Refused(f"{series}: {name!r} is not a count table the series' own folder offers ({sorted(offered)})")
    url = offered[name]["url"]
    with fe._open(url) as r:
        line = first_line(r.read, name.lower().endswith(".gz"), name)
    rows = [r for rs in pl["tables"].values() for r in rs]
    return record_header(ds, dest, series, name, url, line, rows, b["series_metadata_sha256"])


def header_lines(rec: dict) -> list[str]:
    if isinstance(rec.get("ties"), dict):   # the ArrayExpress record (ae.record_header)
        import ae

        return ae.header_lines(rec)
    out = [f"header {rec['name']} of {rec['series']}: status {rec['status']}; {rec['n_columns']} columns; "
           f"header_sha256 {rec['header_sha256']}" + (f"; {rec.get('reason')}" if rec.get("reason") else "")]
    if rec.get("withheld"):
        return out + ["  the header is withheld (a held-out hit): do not pin anything of this table"]
    out.append(f"  the line: {rec['text']}")
    if not fe.plain(rec["text"]):
        out.append(f"  escaped: {ascii(rec['text'])}")
    out.append("  tie table (exact equality only; confirm it, do not decide ties):")
    for t in rec["ties"]:
        h = t["header"] if fe.plain(t["header"]) else ascii(t["header"])
        out.append(f"    [{t['index']}] {h} -> {t['sample'] or '-'} ({t['how']})")
    return out


def stored_headers(dest: Path, acc: str) -> dict:
    """Every header record stored in meta/<ACC>/ (read by the plan, which carries them for verify_pins)."""
    mdir = cm.outside_repo(dest) / "meta" / acc
    out = {}
    for p in sorted(mdir.glob("header_*.json")) if mdir.is_dir() else []:
        if p.is_symlink():
            raise Refused(f"{p}: a symlinked header record: refused")
        rec = json.loads(p.read_text(encoding="utf-8"))
        out[rec["name"]] = rec
    return out


def _header_check(ds: dict, pl: dict, rows: list[dict]) -> dict:
    """The pinned count table's stored header record, checked: in the plan, status ok, the pinned header_sha256 (and
    the hash of its stored line), the url its series' hashed fields give, and its ties derived again from its line and
    the plan's rows equal to the stored ones."""
    t = ds["pinned_series_files"][0]
    refuse_normalised_name(t["name"])
    rec = (pl.get("count_table_headers") or {}).get(t["name"])
    if not rec or rec.get("status") != "ok" or rec.get("withheld"):
        raise Refused(f"{t['name']}: no passing header record in the plan: run the header step, then plan again")
    line = rec["text"].encode("utf-8")
    if header_sha256(line) != t["header_sha256"] or rec["header_sha256"] != t["header_sha256"]:
        raise Refused(f"{t['name']}: the stored header does not hash to the pinned header_sha256: refused")
    if ds["count_form"] == "series_table":
        b = pl["series"][t["series"]]
        derived, _ = series_files(ds, t["series"], b["fields"])
        url = {f["name"]: f["url"] for f in derived}.get(t["name"])
        if url is None or url != rec["url"] or fe.duplicate_names(derived):
            raise Refused(f"{t['name']}: not a count table its series' own hashed fields offer: refused")
        ties, problems = tie_table(split_header(line, t["name"]), rows, ds["gene_columns"])
        if problems or ties != rec["ties"]:
            raise Refused(f"{t['name']}: the tie table derived again differs from the stored one: refused")
    return rec


# ---------------------------------------------------------------------------------------------------------------------
# Pins checked against the plan (fetch.verify_pins dispatches here for a bulk dataset)


def verify_pins(ds: dict, pl: dict, live: bool = False) -> list[dict]:
    """Each pinned series and sample of a bulk GEO dataset checked against the plan, as fetch.verify_pins does for a
    single-cell dataset (series hashes, own-series membership, held links, the sample hash, the screen again, the
    stage under the bulk rule, design_refuse), plus the count source: per_sample its one file derived again from the
    hashed fields; series_table its column, which must be the header the stored header's tie table ties to it."""
    if ds.get("count_form") == "arrayexpress":
        import ae

        return ae.verify_pins(ds, pl, live)
    acc = ds["accession"]
    pins = ds.get("pinned") or []
    if not pins:
        raise Refused(f"{acc}: the pinned block is empty: nothing is fetched or processed until two readers pin "
                      "samples")
    if pl.get("withheld") or "tables" not in pl or acc not in pl["tables"]:
        raise Refused(f"{acc}: the plan holds no readable table for its own series")
    if (pl.get("bulk") or None) != declaration(ds):
        raise Refused(f"{acc}: the plan was made under another bulk declaration: plan again")
    what = "live" if live else "stored"
    sp = cm.check_bulk_series_pins(acc, ds, True)
    blocks = pl.get("series") or {}
    if set(blocks) != set(cm.all_screened(ds)):
        raise Refused(f"{acc}: the plan's series blocks are not exactly the screened series")
    for e in sp:
        b = blocks[e["accession"]]
        h = cm.series_metadata_sha256(e["accession"], b["fields"])
        if h != e["series_metadata_sha256"] or h != b["series_metadata_sha256"]:
            raise Refused(f"series {e['accession']}: its {what} text hashes to {h}, not to the pinned "
                          f"{e['series_metadata_sha256']}: read and pin again")
        bad = held_relations(ds, b["fields"])
        if bad:
            raise Refused(f"series {e['accession']}: its relation names the held series {bad}: refused")
    rows_all = [r for rs in pl["tables"].values() for r in rs]
    hdr = _header_check(ds, pl, rows_all) if ds["count_form"] == "series_table" else None
    tie_of = {t["sample"]: t["header"] for t in (hdr or {}).get("ties", []) if t["sample"]}
    own = {r["accession"]: r for r in pl["tables"][acc]}
    elsewhere = {r["accession"]: s for s, rows in pl["tables"].items() if s != acc for r in rows}
    sfields = blocks[acc]["fields"]
    per: dict[str, dict] = {}
    for p in pins:
        g = p["accession"]
        if g not in own:
            if g in elsewhere:
                raise Refused(f"{g} belongs to another series ({elsewhere[g]}), not to {acc}")
            raise Refused(f"{g} is not in the plan table of {acc}")
        row = own[g]
        if not all(fe.plain_key(k) for k in row["fields"]):
            raise Refused(f"{g}: a field name that is not a plain GEO key: refused")
        why = fe.pin_blocker(ds, row["fields"]) or het_wording(row["fields"])
        if why:
            raise Refused(f"{g}: {why}: refused")
        h = cm.metadata_sha256(g, row["fields"])
        if h != p["metadata_sha256"] or h != row["metadata_sha256"]:
            raise Refused(f"{g}: its {what} metadata hashes to {h}, not to the pinned {p['metadata_sha256']}: read "
                          "and pin again")
        r = sc.screen_record(g, row["fields"], ds.get("screen_extra_refuse"), ds.get("review_cleared"))
        if r["verdict"] != "pass":
            raise Refused(f"{g}: screen verdict {r['verdict']} on the pinned sample: refused")
        e, how = row_stage(ds, row["fields"], sfields)
        if e is None or e != cm.stage_ok(p["stage"]):
            raise Refused(f"{g}: pinned stage {p['stage']} but the record gives {how if e is None else f'E{e:g}'}")
        for pat in ds.get("design_refuse") or []:
            if re.search(pat, fe.sample_text(row["fields"]), re.IGNORECASE):
                raise Refused(f"{g}: design_refuse {pat!r} matches a pinned sample")
        if ds["count_form"] == "per_sample":
            derived, _ = sample_files(ds, row["fields"], g)
            if row.get("files") != derived:
                raise Refused(f"{g}: the stored file list differs from the one its pinned fields give: plan again")
            allowed = {f["name"]: f for f in derived}
            if p["files"] != [f["name"] for f in derived] or fe.denied(p["files"][0]):
                raise Refused(f"{g}: pinned files {p['files']} are not the one count table of its plan entry")
            src = {"files": [allowed[p["files"][0]]]}
        else:
            if tie_of.get(g) != p["column"]:
                raise Refused(f"{g}: pinned column {p['column']!r} is not the header the tie table ties to it "
                              f"({tie_of.get(g)!r}): refused")
            src = {"column": p["column"]}
        s = per.setdefault(g, {"accession": g, "arm": p["arm"], "genes": [], "stage": e, **src,
                               "metadata_sha256": h, "row": row})
        s["genes"].append(p["gene"])
    return list(per.values())


# ---------------------------------------------------------------------------------------------------------------------
# Fetch (fetch.fetch dispatches here)


def _clean_folder(sdir: Path, names: set[str], man: dict, prefix: str) -> None:
    if sdir.is_dir():
        for q in sorted(sdir.iterdir()):
            if q.name in names:
                continue
            if not q.is_file() or q.is_symlink():
                raise Refused(f"{q}: not a plain file in a pinned folder: remove it by hand")
            q.unlink()
            man["files"].pop(f"{prefix}/{q.name}", None)


def _local_first_line(p: Path, name: str) -> bytes:
    with open(p, "rb") as f:
        return first_line(f.read, name.lower().endswith(".gz"), name)


def fetch(acc: str, dest: Path, datasets: Path | None = None, max_bytes: int = fe.MAX_BYTES_DEFAULT) -> dict:
    """The pinned count files of a bulk dataset only: per_sample each pinned sample's one table into
    raw/<ACC>/<GSM>/; series_table and arrayexpress the one pinned count table into raw/<ACC>/_series/, its first
    line checked against header_sha256 right after the download (a file whose first line differs is deleted and
    refused). Each manifest record holds the url, size, sha256, the sample's or series' pinned text hash, the pins in
    force and (for a table) header_sha256; a file is reused only when all equal the current ones. Nothing else is
    requested; nothing is unpacked."""
    ds = cm.dataset(acc, datasets)
    if not ds.get("pinned"):
        raise Refused(f"{acc}: the pinned block is empty: nothing is fetched")
    if ds.get("id_type") == "ensembl":
        cm.load_panel_ids()   # refused while PANEL_IDS_SHA256 is None or the file does not hash to it
    dest = cm.outside_repo(dest)
    pl = fe._load_plan(acc, dest, datasets)
    pinned = fe.verify_pins(ds, pl)
    pp = dest / "meta" / acc / "plan.json"
    rdir = dest / "raw" / acc
    table = (ds.get("pinned_series_files") or [None])[0]
    folders = [rdir / s["accession"] for s in pinned if s.get("files")] + ([rdir / "_series"] if table else [])
    for d in (dest / "raw", rdir, *folders):
        if d.is_symlink():
            raise Refused(f"{d}: a symlinked folder: remove it by hand")
    mf = rdir / "MANIFEST.json"
    if mf.is_symlink():
        raise Refused(f"{mf}: a symlinked manifest: remove it by hand")
    man = json.loads(mf.read_text()) if mf.is_file() else {"accession": acc, "files": {}}
    cur = cm.pins_sha256(ds)
    keep = {s["accession"] for s in pinned if s.get("files")} | ({"_series"} if table else set())
    man["files"] = {k: v for k, v in man["files"].items() if k.split("/")[0] in keep}
    total = 0
    for s in pinned:
        if not s.get("files"):
            continue
        gsm = s["accession"]
        sdir = cm.outside_repo(rdir / gsm)
        f = s["files"][0]
        _clean_folder(sdir, {f["name"]}, man, gsm)
        key, p = f"{gsm}/{f['name']}", sdir / f["name"]
        rec = man["files"].get(key)
        if (rec and p.is_file() and not p.is_symlink() and rec.get("metadata_sha256") == s["metadata_sha256"]
                and rec.get("url") == f["url"] and rec.get("pins_sha256") == cur
                and cm.sha256_file(p) == rec["sha256"]):
            total += rec["size"]
            continue
        man["files"].pop(key, None)
        p.unlink(missing_ok=True)
        size, digest = fe._get(f["url"], p, max_bytes - total)
        total += size
        man["files"][key] = {"url": f["url"], "size": size, "sha256": digest, "gsm": gsm,
                             "metadata_sha256": s["metadata_sha256"], "pins_sha256": cur,
                             "fetched": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
        cm.dump_json(man, mf)
    if table:
        rec_h = pl["count_table_headers"][table["name"]]
        url = rec_h["url"]
        shash = pl["series"][table["series"]]["series_metadata_sha256"]
        sdir = cm.outside_repo(rdir / "_series")
        _clean_folder(sdir, {table["name"]}, man, "_series")
        key, p = f"_series/{table['name']}", sdir / table["name"]
        rec = man["files"].get(key)
        if not (rec and p.is_file() and not p.is_symlink() and rec.get("series_metadata_sha256") == shash
                and rec.get("url") == url and rec.get("pins_sha256") == cur
                and rec.get("header_sha256") == table["header_sha256"] and cm.sha256_file(p) == rec["sha256"]):
            man["files"].pop(key, None)
            p.unlink(missing_ok=True)
            getter = fe._get_ae if ds["count_form"] == "arrayexpress" else fe._get
            size, digest = getter(url, p, max_bytes - total)
            if header_sha256(_local_first_line(p, table["name"])) != table["header_sha256"]:
                p.unlink()
                raise Refused(f"{table['name']}: the downloaded file's first line does not hash to the pinned "
                              "header_sha256: deleted, refused")
            man["files"][key] = {"url": url, "size": size, "sha256": digest, "series": table["series"],
                                 "series_metadata_sha256": shash, "pins_sha256": cur,
                                 "header_sha256": table["header_sha256"],
                                 "fetched": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
            total += size
            cm.dump_json(man, mf)
        else:
            total += rec["size"]
    man.update(plan_sha256=cm.sha256_file(pp), d=ds["d"], total_bytes=total, pinned=ds["pinned"],
               pinned_series=ds["pinned_series"], pinned_series_files=ds.get("pinned_series_files") or [],
               pins_sha256=cur, modality="bulk")
    cm.dump_json(man, mf)
    return man


# ---------------------------------------------------------------------------------------------------------------------
# Pin templates (never write datasets.yaml)


def pin_template(acc: str, dest: Path, arm: str, gene: str, gsms: list[str], datasets: Path | None = None,
                 column: str | None = None) -> str:
    """Ready-to-paste pins of NAMED samples of a bulk GEO dataset: per_sample with its one count table; series_table
    with the column the stored header's tie table ties to it (a --column given must equal it). Refused for a held
    or not pinnable sample, a stage that does not parse under the bulk rule, a sample with no table or no tie."""
    if arm not in cm.PIN_ARMS:
        raise Refused(f"arm must be exactly one of {cm.PIN_ARMS}")
    if not gsms or not all(isinstance(g, str) and re.fullmatch(r"GSM\d+", g) for g in gsms) or len(set(gsms)) != len(
            gsms):
        raise Refused("name every sample to pin, each once (--gsm or --gsm-file; there is no bulk selector)")
    if column is not None and len(gsms) != 1:
        raise Refused("--column names the header of ONE sample")
    pp, pl = fe._stored_plan(acc, dest, datasets)
    ds = cm.dataset(acc, datasets)
    rows = {r["accession"]: r for r in pl["tables"][acc]}
    tie_of = {}
    if ds["count_form"] == "series_table":
        hdrs = [h for h in stored_headers(dest, acc).values() if h.get("status") == "ok"]
        if len(hdrs) != 1:
            raise Refused(f"{acc}: run the header step on the one count table first ({len(hdrs)} passing headers)")
        rows_all = [r for rs in pl["tables"].values() for r in rs]
        ties, problems = tie_table(split_header(hdrs[0]["text"].encode("utf-8"), hdrs[0]["name"]), rows_all,
                                   ds["gene_columns"])
        if problems or ties != hdrs[0]["ties"]:
            raise Refused(f"{acc}: the stored header's tie table does not derive again: run the header step again")
        tie_of = {t["sample"]: t["header"] for t in ties if t["sample"]}
    out = [f"# pin-template {acc} (bulk, {ds['count_form']}): arm {arm}, gene {gene}; {len(gsms)} sample(s) named; "
           f"from {pp}", "# READ each sample's entry in its plan table before pasting (README checks B1 to B10).",
           f"# A bulk unit needs at least {cm.BULK_MIN_LIBS} libraries per arm per stage."]
    sfields = pl["series"][acc]["fields"]
    for g in gsms:
        r = rows.get(g)
        if r is None:
            raise Refused(f"{g} is not in the plan table of {acc}")
        if r["verdict"] != "pass":
            raise Refused(f"{g}: screen verdict {r['verdict']}: not pinned here")
        why = r.get("not_pinnable") or fe.pin_blocker(ds, r["fields"]) or het_wording(r["fields"])
        if why:
            raise Refused(f"{g}: not pinnable: {why}")
        h = cm.metadata_sha256(g, r["fields"])
        if h != r["metadata_sha256"]:
            raise Refused(f"{g}: the stored fields hash to {h}, not to the table's value: plan again")
        e, how = row_stage(ds, r["fields"], sfields)
        if e is None or e != r["stage"]:
            raise Refused(f"{g}: no single stage under the dataset's stage rule ({how}): it cannot be pinned")
        if ds["count_form"] == "per_sample":
            derived, _ = sample_files(ds, r["fields"], g)
            if r["files"] != derived or len(derived) != 1:
                raise Refused(f"{g}: not exactly one count table of its own folder offered: it cannot be pinned")
            src = f"        files: [{json.dumps(derived[0]['name'])}]"
        else:
            col = tie_of.get(g)
            if col is None or (column is not None and column != col):
                raise Refused(f"{g}: no header ties to it exactly" if col is None else
                              f"{g}: --column {column!r} is not the header the tie table ties to it ({col!r})")
            src = f"        column: {json.dumps(col)}"
        out += [f"      # ---- {g}: its fields as read in the plan table", *fe._comment_lines(r["fields"], "      "),
                f"      - accession: {g}", f"        arm: {arm}", f"        gene: {json.dumps(gene)}",
                f"        stage: E{e:g}", src, f'        metadata_sha256: "{h}"']
    return "\n".join(out) + "\n"


def table_template(acc: str, dest: Path, datasets: Path | None = None) -> str:
    """The ready-to-paste pinned_series_files entry of the one count table whose header step passed."""
    ds = cm.dataset(acc, datasets)
    if not cm.is_bulk(ds) or ds.get("count_form") not in ("series_table", "arrayexpress"):
        raise Refused(f"{acc}: a count table is pinned only for a bulk series_table or arrayexpress dataset")
    hdrs = [h for h in stored_headers(dest, acc).values() if h.get("status") == "ok"]
    if len(hdrs) != 1:
        raise Refused(f"{acc}: {len(hdrs)} passing header records: run the header step on the one count table")
    h = hdrs[0]
    return "\n".join([f"# table-template {acc}: {h['name']} of {h['series']}, {h['n_columns']} columns; READ the tie "
                      "table the header step printed before pasting",
                      "# Paste under the dataset's `pinned_series_files:` key.",
                      f"      - series: {h['series']}", f"        name: {json.dumps(h['name'])}",
                      f'        header_sha256: "{h["header_sha256"]}"']) + "\n"


# ---------------------------------------------------------------------------------------------------------------------
# Reading counts (pure, no network)


def _read_table(src, name: str, usecols=None) -> pd.DataFrame:
    """A text table as strings (no NA conversion), the separator from the name, gzip from the name."""
    try:
        tab = sep_of(name) == "\t"
        # round 2 (the first reviewer's lesser point 3): a tab-separated table is read without quote handling, so a
        # quote in an unpinned text column cannot merge rows
        return pd.read_csv(src, sep=sep_of(name), header=0, dtype=str, keep_default_na=False, na_filter=False,
                           usecols=usecols, compression="gzip" if name.lower().endswith(".gz") else None,
                           engine="c", **({"quoting": csv.QUOTE_NONE} if tab else {}))
    except (ValueError, pd.errors.ParserError) as e:
        raise Refused(f"{name}: not a table this intake reads ({type(e).__name__}): refused") from None


def counts_column(values, label: str) -> np.ndarray:
    """One count column as float64: every value a finite, non-negative number integral within INTEGRAL_TOL."""
    try:
        x = np.array([float(v) for v in values], dtype=np.float64)
    except ValueError:
        raise Refused(f"{label}: a value that is not a number: refused (counts only)") from None
    if not np.isfinite(x).all() or (x < 0).any():
        raise Refused(f"{label}: a negative or non-finite value: refused (counts only)")
    if (np.abs(x - np.round(x)) > cm.INTEGRAL_TOL).any():
        raise Refused(f"{label}: a non-integral value (expected counts, FPKM or TPM): refused (counts only)")
    return np.round(x)


def screen_sample_header(ds: dict, gsm: str, fields: list[str], text: str) -> None:
    """Round 2 (the first reviewer's blocker 2): a per-sample count file's first line, screened before any value is
    parsed. A held-out hit refuses with the text withheld; a review term refuses unless a review_cleared entry names
    this sample."""
    held, bad = screen_header(ds, gsm, fields, text)
    if held:
        raise Refused(f"{gsm}: a held-out hit in a per-sample header (text withheld): refused")
    if bad:
        raise Refused(f"{gsm}: per-sample header screen: {bad}: refused")


def read_count_file(path: Path, name: str, gene_columns: list[str], count_column: str | None,
                    tie_hits, gsm: str, screen) -> tuple[list[str], list[str] | None, np.ndarray, str]:
    """One per-sample count table: (ids, symbols or None, counts, the count column's header). Its first line is
    screened first (screen(fields, text), before any value is parsed). The first line must hold the declared
    gene_columns (a headerless file is refused); exactly one other column is the count column, unless count_column
    names one, which must tie to this sample by the exact rule (tie_hits(header) -> the samples it ties to, exactly
    [gsm]); more than one other column without count_column is refused (a Length column, a second sample). Round 2
    (the first reviewer's blocker 3): a single count column that ties by the exact rule to ANOTHER sample is refused
    (its hits must be [] or exactly [gsm])."""
    with open(path, "rb") as f:
        line = first_line(f.read, name.lower().endswith(".gz"), name)
    hdr = split_header(line, name)
    screen(hdr, line.decode("utf-8"))
    missing = [c for c in gene_columns if hdr.count(c) != 1]
    if missing:
        raise Refused(f"{name}: its first line does not hold the declared gene column(s) {missing} once: a headerless "
                      "or other table, refused")
    rest = [h for h in hdr if h not in gene_columns]
    if count_column is None:
        if len(rest) != 1:
            raise Refused(f"{name}: {len(rest)} columns besides the gene columns and no count_column: refused (a "
                          "second integral column such as Length, or a table of several samples)")
        col = rest[0]
    else:
        if hdr.count(count_column) != 1 or count_column in gene_columns:
            raise Refused(f"{name}: count_column {count_column!r} is not one header of the file: refused")
        col = count_column
    hits = tie_hits(col)
    if count_column is not None or len(rest) > 1:
        if hits != [gsm]:
            raise Refused(f"{name}: the count column {col!r} does not tie to this sample by the exact tie rule: "
                          "refused")
    elif hits not in ([], [gsm]):
        raise Refused(f"{name}: the count column ties by the exact rule to another sample: refused")
    if NORMALISED.search(col) or NORMALISED.search(name):
        raise Refused(f"{name}: a column or file named as normalised or expected counts: refused (counts only)")
    df = _read_table(path, name, usecols=[*gene_columns, col])
    ids = df[gene_columns[0]].tolist()
    sym = df[gene_columns[1]].tolist() if len(gene_columns) > 1 else None
    return ids, sym, counts_column(df[col].tolist(), f"{name} {col}"), col


def read_series_table(path: Path, name: str, gene_columns: list[str], columns: list[str],
                      header_sha: str) -> tuple[list[str], list[str] | None, dict[str, np.ndarray]]:
    """The pinned count table: refused when it is named as normalised or expected counts, or unless its first line
    hashes to header_sha; read with usecols (the gene columns and the pinned columns only), so no unpinned column is
    parsed."""
    refuse_normalised_name(name)
    line = _local_first_line(path, name)
    if header_sha256(line) != header_sha:
        raise Refused(f"{name}: its first line does not hash to the pinned header_sha256: refused")
    hdr = split_header(line, name)
    for c in [*gene_columns, *columns]:
        if hdr.count(c) != 1:
            raise Refused(f"{name}: column {c!r} is not exactly one header of the table: refused")
    for c in columns:
        if NORMALISED.search(c):
            raise Refused(f"{name}: a column named as normalised or expected counts: refused (counts only)")
    df = _read_table(path, name, usecols=[*gene_columns, *columns])
    ids = df[gene_columns[0]].tolist()
    sym = df[gene_columns[1]].tolist() if len(gene_columns) > 1 else None
    return ids, sym, {c: counts_column(df[c].tolist(), f"{name} {c}") for c in columns}


def counts_only(totals: list[float], label: str) -> None:
    """The CPM or TPM signature: every library total within TPM_SIGNATURE_TOL of 1e6 refuses the dataset."""
    t = np.asarray(totals, dtype=np.float64)
    if len(t) and (np.abs(t - 1e6) <= cm.TPM_SIGNATURE_TOL * 1e6).all():
        raise Refused(f"{label}: every library sums to within 1% of 1e6 (a CPM or TPM table): refused (counts only)")


def special(i: str) -> bool:
    return i.startswith(SPECIAL_ROW)


def library_total(ids: list[str], x: np.ndarray) -> float:
    """T_i: the library's total over all gene rows, HTSeq '__*' and STAR 'N_*' rows left out."""
    keep = np.array([not special(i) for i in ids], dtype=bool)
    return float(x[keep].sum()) if len(x) else 0.0


class IdMap:
    """The id route of one dataset: ensembl (version stripped, through panel_ids.tsv) or symbol (process.to_panel's
    rule: exact symbol, then the committed ALIASES; duplicates summed)."""

    def __init__(self, id_type: str, genes: list[str], panel: dict | None):
        self.id_type, self.genes, self.panel = id_type, list(genes), panel
        self._cache: dict[int, tuple[list, list[str], dict[str, int]]] = {}
        if id_type == "ensembl":
            if panel is None:
                raise Refused("the Ensembl route needs the pinned gene-id map (panel_ids.tsv)")
            self.id_to_symbol = {i: s for s, i in zip(panel["symbols"], panel["ids"], strict=True) if i}
            self.extra = dict(panel.get("extra") or {})

    def _stripped(self, ids: list[str]) -> tuple[list[str], dict[str, int]]:
        """(the ids, version stripped, and their row positions), cached per id list; two rows of one stripped id
        refuse the dataset."""
        c = self._cache.get(id(ids))
        if c is not None and c[0] is ids:
            return c[1], c[2]
        out = [VERSION.sub("", i) for i in ids]
        pos: dict[str, int] = {}
        dup = 0
        for k, i in enumerate(out):
            if special(i):
                continue
            if i in pos:
                dup += 1
            pos.setdefault(i, k)
        if dup:
            raise Refused(f"{dup} table rows repeat a stripped Ensembl id: refused")
        self._cache[id(ids)] = (ids, out, pos)
        return out, pos

    def strip(self, ids: list[str]) -> list[str]:
        if self.id_type != "ensembl":
            return list(ids)
        return self._stripped(ids)[0]

    def to_panel(self, ids: list[str], x: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
        """(counts on the panel genes in the reference order, presence, record)."""
        if self.id_type == "symbol":
            import process as pr

            keep = [k for k, i in enumerate(ids) if not special(i)]
            M = sparse.csr_matrix(np.asarray(x, dtype=np.float64)[keep][None, :])
            C, present, rec = pr.to_panel(M, [ids[k] for k in keep], ["Gene Expression"] * len(keep), self.genes)
            return np.asarray(C.toarray()).ravel().astype(np.float64), present, rec
        _, pos = self._stripped(ids)
        out = np.zeros(len(self.genes))
        present = np.zeros(len(self.genes), dtype=bool)
        for j, gid in enumerate(self.panel["ids"]):
            if gid is not None and gid in pos:
                out[j] = x[pos[gid]]
                present[j] = True
        return out, present, {"route": "ensembl", "matched": int(present.sum()), "of": len(self.genes),
                              "na_in_map": int(sum(i is None for i in self.panel["ids"]))}

    def lookup(self, ids: list[str], x: np.ndarray, symbol: str) -> float | None:
        """The count of one gene by symbol (a sex gene or the knocked-out gene); None when it cannot be found."""
        if self.id_type == "symbol":
            hits = [k for k, i in enumerate(ids) if i == symbol]
            if not hits and symbol in cm.ALIASES:
                hits = [k for k, i in enumerate(ids) if i == cm.ALIASES[symbol]]
            return float(x[hits].sum()) if hits else None
        gid = self.extra.get(symbol)
        if gid is None and symbol in self.panel["symbols"]:
            gid = self.panel["ids"][self.panel["symbols"].index(symbol)]
        if gid is None:
            return None
        k = self._stripped(ids)[1].get(gid)
        return float(x[k]) if k is not None else None

    def symbol_disagreements(self, ids: list[str], symbols: list[str] | None) -> int | None:
        if symbols is None or self.id_type != "ensembl":
            return None
        s = self.strip(ids)
        return int(sum(1 for i, y in zip(s, symbols, strict=True) if i in self.id_to_symbol
                       and self.id_to_symbol[i] != y))


# ---------------------------------------------------------------------------------------------------------------------
# The math and the guards (pure numpy; synthetic tests run these)


def log_profile(x: np.ndarray, T: float) -> np.ndarray:
    return np.log((np.asarray(x, dtype=np.float64) + cm.PSEUDOCOUNT) / T)


def bulk_log_ratio(X_ko: np.ndarray, X_ctl: np.ndarray, T_ko, T_ctl, present, nan_mask=None) -> np.ndarray:
    """L (P,) of one unit at one stage: the mean over knockout libraries of y_gi = log((x_gi + 10) / T_i) minus the
    mean over control libraries, clipped to log [0.25, 4]; NaN where a gene is absent, where the median count over an
    arm's libraries is under MIN_UMI in either arm, and where nan_mask is True. X_*: libraries x panel genes."""
    X_ko, X_ctl = np.atleast_2d(np.asarray(X_ko, float)), np.atleast_2d(np.asarray(X_ctl, float))
    yk = np.vstack([log_profile(x, t) for x, t in zip(X_ko, T_ko, strict=True)]).mean(0)
    yc = np.vstack([log_profile(x, t) for x, t in zip(X_ctl, T_ctl, strict=True)]).mean(0)
    L = np.clip(yk - yc, *cm.LOG_CLIP)
    low = (np.median(X_ko, 0) < cm.MIN_UMI) | (np.median(X_ctl, 0) < cm.MIN_UMI)
    bad = low | ~np.asarray(present, bool)
    if nan_mask is not None:
        bad |= np.asarray(nan_mask, bool)
    return np.where(bad, np.nan, L)


def own_gene_ratio(x_ko: list[float], x_ctl: list[float], T_ko, T_ctl) -> float:
    """G9: the knocked-out gene's own clipped log ratio (the same mean of logs)."""
    yk = np.mean([math.log((x + cm.PSEUDOCOUNT) / t) for x, t in zip(x_ko, T_ko, strict=True)])
    yc = np.mean([math.log((x + cm.PSEUDOCOUNT) / t) for x, t in zip(x_ctl, T_ctl, strict=True)])
    return float(np.clip(yk - yc, *cm.LOG_CLIP))


def own_gene_guard(L_own: float | None) -> tuple[bool, str]:
    """G9: refused when the own gene cannot be found or its clipped log ratio is >= 0 (the knockout arm not lower).
    Only the verdict is recorded, never the value."""
    if L_own is None:
        return False, "the knocked-out gene is not found in the counts: the genotype cannot be checked"
    if L_own >= 0:
        return False, "genotype not supported in the counts (the knockout arm is not lower on its own gene)"
    return True, "the knockout arm is lower on its own gene"


def sex_call(xist: float | None, y: list[float | None], T: float) -> str:
    """G3: 'male' when the Y genes' CPM sum is at least SEX_CPM and Xist's CPM under it, 'female' the reverse, else
    'unclear'; 'not possible' when a sex gene is missing."""
    if xist is None or any(v is None for v in y) or T <= 0:
        return "not possible"
    yc = sum(y) / T * 1e6
    xc = xist / T * 1e6
    if yc >= cm.SEX_CPM and xc < cm.SEX_CPM:
        return "male"
    if yc < cm.SEX_CPM and xc >= cm.SEX_CPM:
        return "female"
    return "unclear"


def sex_balance(calls_ko: list[str], calls_ctl: list[str]) -> tuple[str, dict]:
    """G3: 'balanced', 'imbalanced' (the male share differs by more than SEX_SHARE_GAP) or 'not possible' (a sex gene
    absent, or an arm with no clear call)."""
    rec = {}
    if "not possible" in calls_ko + calls_ctl:
        return "not possible", rec
    shares = {}
    for arm, calls in (("ko", calls_ko), ("control", calls_ctl)):
        clear = [c for c in calls if c in ("male", "female")]
        if not clear:
            return "not possible", rec
        shares[arm] = clear.count("male") / len(clear)
    rec["male_share"] = shares
    return ("imbalanced" if abs(shares["ko"] - shares["control"]) > cm.SEX_SHARE_GAP else "balanced"), rec


def _spearman(a, b) -> float:
    ra, rb = rankdata(a), rankdata(b)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def outlier_drop(Y: np.ndarray) -> tuple[list[bool], list[dict]]:
    """G4: in an arm of 3 or more libraries (rows of log CPM on the present panel genes), a library is dropped when
    its Spearman with the leave-one-out arm mean has a robust z under -OUTLIER_Z AND lies more than OUTLIER_RHO_GAP
    under the arm's median Spearman. Returns (keep per library, records)."""
    n = len(Y)
    if n < 3:
        return [True] * n, [{"rho": None, "note": "arm under 3 libraries: no outlier test"} for _ in range(n)]
    rho = np.array([_spearman(Y[i], np.delete(Y, i, 0).mean(0)) for i in range(n)])
    med = float(np.nanmedian(rho))
    mad = float(np.nanmedian(np.abs(rho - med)))
    keep, recs = [], []
    for r in rho:
        if mad > 0:
            z = (r - med) / (1.4826 * mad)
        else:
            z = -math.inf if r < med else 0.0
        drop = bool(np.isfinite(r) and z < -cm.OUTLIER_Z and r < med - cm.OUTLIER_RHO_GAP)
        keep.append(not drop)
        recs.append({"rho": float(r), "z": (None if not np.isfinite(z) else float(z)), "dropped": drop})
    return keep, recs


def control_grouping(Y_ko: np.ndarray, Y_ctl: np.ndarray) -> tuple[bool, dict]:
    """G4: the unit is refused when the mean Spearman within the controls is under the mean Spearman across the arms
    minus OUTLIER_RHO_GAP (the controls do not group: a label swap or a mixed design)."""
    within = [_spearman(Y_ctl[i], Y_ctl[j]) for i in range(len(Y_ctl)) for j in range(i + 1, len(Y_ctl))]
    cross = [_spearman(a, b) for a in Y_ko for b in Y_ctl]
    w, c = float(np.nanmean(within)), float(np.nanmean(cross))
    return not (w < c - cm.OUTLIER_RHO_GAP), {"within_control": w, "cross_arm": c}


def whole_embryo_guard(profile: np.ndarray, ref: dict, present) -> tuple[bool, dict]:
    """G5: the profile's Spearman with pb must exceed its Spearman with every single reference group (margin >
    BULK_WHOLE_MARGIN, strict; present genes). A check against single- or few-lineage samples, not proof of a whole
    embryo. pb is the reference's whole-embryo pseudo-bulk (tb.whole_den)."""
    pr_ = np.asarray(present, bool)
    p = np.asarray(profile, float)[pr_]
    rho_pb = _spearman(p, tb.whole_den(ref)[pr_])
    best = max((_spearman(p, np.asarray(ref["mu"][g])[pr_]), g) for g in ref["groups"])
    margin = rho_pb - best[0]
    ok = bool(np.isfinite(rho_pb) and np.isfinite(best[0]) and margin > cm.BULK_WHOLE_MARGIN)
    return ok, {"rho_whole": rho_pb, "rho_best_group": best[0], "best_group": best[1], "margin": float(margin),
                "threshold": cm.BULK_WHOLE_MARGIN, "rule": "margin > threshold"}


def composition_shares(profile: np.ndarray, ref: dict, present) -> dict[str, float]:
    """G12: the NNLS fit of a mean panel proportion profile onto the reference groups' mean profiles mu[g] (present
    genes), the weights normalised to shares."""
    pr_ = np.asarray(present, bool)
    A = np.vstack([np.asarray(ref["mu"][g], float)[pr_] for g in ref["groups"]]).T
    w, _ = nnls(A, np.asarray(profile, float)[pr_])
    s = w.sum()
    return {g: (float(v / s) if s > 0 else 0.0) for g, v in zip(ref["groups"], w, strict=True)}


def composition_guard(P_ko: np.ndarray, P_ctl: np.ndarray, ref: dict, present) -> tuple[bool, dict]:
    """G12: refused when any group with a control share of at least COMPOSITION_MIN_SHARE has a knockout-to-control
    share ratio above COMPOSITION_RATIO_MAX or under its inverse. Group-level shares are recorded."""
    sk = composition_shares(np.mean(P_ko, 0), ref, present)
    scn = composition_shares(np.mean(P_ctl, 0), ref, present)
    bad = []
    for g in ref["groups"]:
        if scn[g] >= cm.COMPOSITION_MIN_SHARE:
            r = sk[g] / scn[g]
            if r > cm.COMPOSITION_RATIO_MAX or r < 1.0 / cm.COMPOSITION_RATIO_MAX:
                bad.append(g)
    return not bad, {"shares_ko": sk, "shares_control": scn, "groups_moved": bad}


def somite_record(ko: list[int | None], ctl: list[int | None]) -> tuple[bool, dict]:
    """G6: where both arms record somites, non-overlapping ranges refuse; the median difference (knockout minus
    control), the ranges and the overlap count are recorded. Nothing is removed."""
    k = [v for v in ko if v is not None]
    c = [v for v in ctl if v is not None]
    if not k or not c:
        return True, {"recorded": False}
    rec = {"recorded": True, "range_ko": [min(k), max(k)], "range_control": [min(c), max(c)],
           "median_difference": float(np.median(k) - np.median(c)),
           "ko_within_control_range": int(sum(min(c) <= v <= max(c) for v in k)),
           "complete": len(k) == len(ko) and len(c) == len(ctl)}
    ok = not (max(k) < min(c) or max(c) < min(k))
    return ok, rec


def duplicate_vectors(libs: list[dict]) -> list[list[str]]:
    """G8: groups of pinned libraries with identical count vectors (ids and counts)."""
    seen: dict[str, list[str]] = {}
    for lb in libs:
        seen.setdefault(lb["vector_sha256"], []).append(lb["label"])
    return [v for v in seen.values() if len(v) > 1]


# ---------------------------------------------------------------------------------------------------------------------
# Sex and somite fields


SOMITE = re.compile(r"^(\d+)(?:-somite stage| somites)?$")


def field_value(rec: dict, spec: str | None) -> str | None:
    """The one value of a sample field: 'characteristics_ch1:<key>' (the value part of the entry whose key equals
    <key>) or a plain field name; None when absent or not exactly one value."""
    if not spec:
        return None
    if spec.startswith("characteristics_ch1:"):
        key = spec.split(":", 1)[1]
        vals = [v.split(": ", 1)[1] for v in rec.get("characteristics_ch1", []) if ": " in v
                and v.split(": ", 1)[0] == key]
    else:
        vals = list(rec.get(spec, []))
    return vals[0] if len(vals) == 1 else None


def recorded_sex(v: str | None) -> str | None:
    return v if v in ("male", "female") else None


def somites(v: str | None) -> int | None:
    m = SOMITE.fullmatch(v) if isinstance(v, str) else None
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------------------------------------------------------------
# Orchestration (process.run dispatches here)


def _check_manifest(ds, pl, man, rdir, pinned, cur, table) -> dict[str, str]:
    """Every pinned file bound to the pin it was fetched under; no other file in a pinned folder; no symlink."""
    count_sha = {}
    for s in pinned:
        if not s.get("files"):
            continue
        gsm = s["accession"]
        sdir = rdir / gsm
        names = {f["name"] for f in s["files"]}
        if sdir.is_symlink():
            raise Refused(f"{gsm}: a symlinked sample folder: refused")
        for q in sorted(sdir.iterdir()) if sdir.is_dir() else []:
            if q.is_symlink():
                raise Refused(f"{gsm}/{q.name}: a symlink in a pinned sample's folder: refused")
            if q.name not in names or f"{gsm}/{q.name}" not in man["files"]:
                raise Refused(f"{gsm}/{q.name}: a file in a pinned sample's folder that is not its pinned file with a "
                              "manifest record: refused (fetch again)")
        for f in s["files"]:
            key = f"{gsm}/{f['name']}"
            rec = man["files"].get(key)
            if rec is None:
                raise Refused(f"{key}: no manifest record: refused (fetch again)")
            if rec.get("metadata_sha256") != s["metadata_sha256"] or rec.get("url") != f["url"] or rec.get(
                    "pins_sha256") != cur:
                raise Refused(f"{key}: fetched under another metadata_sha256, url or pins: refused (fetch again)")
            p = sdir / f["name"]
            if p.is_symlink() or not p.is_file() or cm.sha256_file(p) != rec["sha256"]:
                raise Refused(f"{key}: missing or its sha256 differs from the manifest")
            count_sha[key] = rec["sha256"]
    if table:
        sdir = rdir / "_series"
        if sdir.is_symlink():
            raise Refused("_series: a symlinked folder: refused")
        for q in sorted(sdir.iterdir()) if sdir.is_dir() else []:
            if q.is_symlink() or q.name != table["name"]:
                raise Refused(f"_series/{q.name}: not the pinned count table: refused (fetch again)")
        key = f"_series/{table['name']}"
        rec = man["files"].get(key)
        hrec = pl["count_table_headers"][table["name"]]
        shash = pl["series"][table["series"]]["series_metadata_sha256"]
        if (rec is None or rec.get("url") != hrec["url"] or rec.get("pins_sha256") != cur
                or rec.get("series_metadata_sha256") != shash or rec.get("header_sha256") != table["header_sha256"]):
            raise Refused(f"{key}: no manifest record under these pins, url, series text and header: refused (fetch "
                          "again)")
        p = sdir / table["name"]
        if p.is_symlink() or not p.is_file() or cm.sha256_file(p) != rec["sha256"]:
            raise Refused(f"{key}: missing or its sha256 differs from the manifest")
        count_sha[key] = rec["sha256"]
    return count_sha


def _tie_hits(pl: dict):
    """hits(header) -> the samples of every plan table the header ties to by the exact rule (sorted)."""
    rows = [r for rs in pl["tables"].values() for r in rs]
    samples = {}
    for r in rows:
        samples.setdefault(r["accession"], r["fields"])
    cands = {g: tie_candidates(g, rec) for g, rec in samples.items()}

    def hits(h: str) -> list[str]:
        return sorted(g for g, c in cands.items() if h in c)
    return hits


def read_libraries(ds: dict, pl: dict, pinned: list[dict], rdir: Path) -> tuple[list[dict], dict]:
    """Each pinned library's (ids, symbols, counts): per_sample from its one file, series_table from its pinned
    column, arrayexpress from the sum of its pinned runs' columns."""
    cf = ds["count_form"]
    libs, rec = [], {}
    if cf == "per_sample":
        for s in pinned:
            f = s["files"][0]
            g = s["accession"]
            ids, sym, x, col = read_count_file(rdir / g / f["name"], f["name"], ds["gene_columns"],
                                               ds.get("count_column"), _tie_hits(pl), g,
                                               lambda fields, text, g=g: screen_sample_header(ds, g, fields, text))
            libs.append({"label": s["accession"], "pin": s, "ids": ids, "symbols": sym, "x": x, "column": col})
        return libs, rec
    table = ds["pinned_series_files"][0]
    if cf == "series_table":
        cols = [s["column"] for s in pinned]
    else:
        cols = [r for s in pinned for r in s["runs"]]
    ids, sym, X = read_series_table(rdir / "_series" / table["name"], table["name"], ds["gene_columns"], cols,
                                    table["header_sha256"])
    for s in pinned:
        if cf == "series_table":
            x, col = X[s["column"]], s["column"]
        else:
            x, col = np.sum([X[r] for r in s["runs"]], axis=0), "+".join(s["runs"])
        libs.append({"label": s["accession"], "pin": s, "ids": ids, "symbols": sym, "x": x, "column": col})
    rec["columns_read"] = len(cols)
    return libs, rec


def _sex_somite(ds: dict, s: dict) -> tuple[str | None, int | None]:
    if ds["count_form"] == "arrayexpress":
        return s.get("sex"), somites(s.get("somite"))
    return recorded_sex(field_value(s["row"]["fields"], ds.get("sex_field"))), somites(
        field_value(s["row"]["fields"], ds.get("somite_field")))


def run_bulk(acc: str, dest: Path, wt: Path, refstats: Path, out: Path, datasets: Path | None = None,
             panel_path: Path | None = None, panel_sha: str | None = None) -> dict:
    """One fetched bulk dataset to one unit table per knockout gene and QC.json, keeping process.run's checks in
    order: the stale output deleted first, the plan (passing, the committed datasets.yaml), the manifest under these
    pins, verify_pins (the screen again before any count is read), every file bound to its pin, no symlink, each count
    file's sha256; then the counts (G2, G8), the id route, and per unit the guards and the math."""
    import process as pr

    ds = cm.dataset(acc, datasets)
    if not cm.is_bulk(ds):
        raise Refused(f"{acc}: not a bulk dataset")
    if not ds.get("pinned"):
        raise Refused(f"{acc}: the pinned block is empty: nothing is processed")
    dest, out = cm.outside_repo(dest), cm.outside_repo(out)
    adir = out / acc
    if adir.is_symlink():
        raise Refused(f"{adir}: a symlinked output folder: remove it by hand")
    if adir.exists():
        shutil.rmtree(adir)
    pl = fe._load_plan(acc, dest, datasets)
    plan_sha = cm.sha256_file(dest / "meta" / acc / "plan.json")
    rdir = dest / "raw" / acc
    mp = rdir / "MANIFEST.json"
    for d in (dest / "raw", rdir, mp):
        if d.is_symlink():
            raise Refused(f"{d}: a symlink in the raw folder: refused")
    if not mp.is_file():
        raise Refused(f"{mp}: no manifest; fetch first")
    man = json.loads(mp.read_text())
    if man.get("plan_sha256") != plan_sha:
        raise Refused(f"{acc}: the manifest was not fetched from this plan")
    cur = cm.pins_sha256(ds)
    if (man.get("pins_sha256") != cur or man.get("pinned") != ds["pinned"]
            or man.get("pinned_series") != ds["pinned_series"]
            or (man.get("pinned_series_files") or []) != (ds.get("pinned_series_files") or [])):
        raise Refused(f"{acc}: the samples, arms, genes, stages or count table are not exactly as pinned when the "
                      "files were fetched: fetch again")
    pinned = fe.verify_pins(ds, pl)
    table = (ds.get("pinned_series_files") or [None])[0]
    count_sha = _check_manifest(ds, pl, man, rdir, pinned, cur, table)
    ref = tb.load_ref_stats(refstats)
    if ref["wt_sha256"] != cm.sha256_file(wt):
        raise Refused("the WT file is not the one the reference statistics were computed from")
    genes = ref["genes"]
    panel = None
    if ds["id_type"] == "ensembl" or panel_path is not None:
        panel = cm.load_panel_ids(genes, panel_path, panel_sha)
    else:
        try:
            panel = cm.load_panel_ids(genes)
        except Refused:
            panel = None   # the symbol route runs without the map; chromosomes are then unknown (G3)
    idm = IdMap(ds["id_type"], genes, panel)
    # G8 (files): two pinned per-sample files of one sha256 refuse the dataset
    shas = [v for k, v in count_sha.items() if not k.startswith("_series/")]
    if len(shas) != len(set(shas)):
        raise Refused(f"{acc}: two pinned libraries share one count file (sha256): refused (G8)")
    try:
        libs, read_rec = read_libraries(ds, pl, pinned, rdir)
    except Refused as e:   # G2 and the readers' refusals fail the dataset's QC (combine then refuses it)
        libs, read_rec = None, {"refused": str(e)}
    qc = {"accession": acc, "modality": "bulk", "count_form": ds["count_form"], "id_type": ds["id_type"],
          "tier": ds["tier"], "study": ds["study"], "d": ds["d"], "screened_samples": [s["accession"] for s in pinned],
          "pinned": ds["pinned"], "pinned_series": ds["pinned_series"], "pins_sha256": cur, "plan_sha256": plan_sha,
          "manifest_sha256": cm.sha256_file(mp), "count_files_sha256": count_sha, "read": read_rec,
          "panel_genes_present": 0, "panel_genes": len(genes),
          "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "git": pr._git_head()}
    if table:
        qc["pinned_series_files"] = ds["pinned_series_files"]
    if libs is None:
        qc.update(status="failed", reason=read_rec["refused"])
        cm.dump_json(qc, adir / "QC.json")
        return qc
    # G2: counts only (each column was read through counts_column); the CPM/TPM signature over every library
    for lb in libs:
        lb["T"] = library_total(lb["ids"], lb["x"])
        h = hashlib.sha256()
        h.update("\n".join(lb["ids"]).encode("utf-8"))
        h.update(np.asarray(lb["x"], dtype=np.float64).tobytes())
        lb["vector_sha256"] = h.hexdigest()
    try:
        counts_only([lb["T"] for lb in libs], acc)
    except Refused as e:
        qc.update(status="failed", reason=str(e))
        cm.dump_json(qc, adir / "QC.json")
        return qc
    dups = duplicate_vectors(libs)
    if dups:
        raise Refused(f"{acc}: pinned libraries with identical count vectors {dups}: refused (G8)")
    presence, panel_rec = [], {}
    for lb in libs:
        lb["xp"], pres, panel_rec[lb["label"]] = idm.to_panel(lb["ids"], lb["x"])
        presence.append(pres)
        lb["sex_counts"] = {g: idm.lookup(lb["ids"], lb["x"], g) for g in cm.SEX_GENES}
        lb["symbol_disagreements"] = idm.symbol_disagreements(lb["ids"], lb["symbols"])
    present = np.logical_and.reduce(presence)
    qc.update(panel_match=panel_rec, panel_genes_present=int(present.sum()), panel_genes=len(genes))
    if present.mean() < cm.MIN_PANEL_MATCH:
        qc.update(status="failed", reason=f"only {present.sum()} of {len(genes)} panel genes matched")
        cm.dump_json(qc, adir / "QC.json")
        return qc
    chrom = panel["chromosome"] if panel is not None else [None] * len(genes)
    always = np.array([g in cm.ALWAYS_NAN_BULK or g in cm.SEX_Y_GENES or c == "Y" for g, c in
                       zip(genes, chrom, strict=True)])
    chrx = np.array([c == "X" for c in chrom])
    chrom_known = all(c is not None for c in chrom)
    lib_rec = {}
    for lb in libs:
        s = lb["pin"]
        sex, som = _sex_somite(ds, s)
        lb["sex_recorded"], lb["somites"] = sex, som
        lb["sex_call"] = sex_call(lb["sex_counts"]["Xist"], [lb["sex_counts"][g] for g in cm.SEX_Y_GENES], lb["T"])
        lb["panel_total"] = float(lb["xp"][present].sum())
        lib_rec[lb["label"]] = {"total": lb["T"], "panel_total": lb["panel_total"], "sex_call": lb["sex_call"],
                                "sex_recorded": sex, "somites": som, "column": lb["column"],
                                "symbol_disagreements": lb["symbol_disagreements"]}
    qc["libraries"] = lib_rec
    by = {lb["label"]: lb for lb in libs}
    qc["units"] = {}
    for u in cm.pinned_units(ds["pinned"]):
        uqc = bulk_unit(ds, u, by, ref, present, always, chrx, chrom_known, idm)
        if uqc.get("status") == "ok":
            uqc.update(write_bulk_unit(adir, acc, u["gene"], uqc.pop("_L"), ref, {
                "tier": ds["tier"], "study": ds["study"], "d": ds["d"], "ko_samples": u["ko"],
                "control_samples": u["control"], "present": present.tolist(), "git": qc["git"],
                "k": uqc["k_by_stage"], "m": uqc["m_by_stage"]}))
        else:
            uqc.pop("_L", None)
        qc["units"][u["gene"]] = uqc
    qc["status"] = "ok" if any(v["status"] == "ok" for v in qc["units"].values()) else "left_out"
    cm.dump_json(qc, adir / "QC.json")
    return qc


def bulk_unit(ds, u, by, ref, present, always, chrx, chrom_known, idm) -> dict:
    """One knockout gene's bulk unit: per stage the guards in order (G1, G3 drop, G4 drop, G5 drop, the floor, G4
    grouping, G5 arm mean, G6, G12, G9) and the math; stages averaged with equal weight. Returns the unit QC with the
    log table under '_L' when it passes."""
    genes = ref["genes"]
    gene = u["gene"]
    own = np.array([g == gene or cm.ALIASES.get(g) == gene for g in genes])
    uqc = {"gene": gene, "g5_rule": "bulk margin > 0 (BULK_WHOLE_MARGIN), decision 2026-11-21", "stages": {},
           "k_by_stage": {}, "m_by_stage": {}}
    Ls = []
    for e in u["stages"]:
        st = {"dropped": {}}
        uqc["stages"][f"E{e:g}"] = st
        arms = {a: [by[g] for g in u[a] if cm.stage_ok(f"E{by[g]['pin']['stage']:g}") == e]
                for a in ("ko", "control")}
        tot = [lb["T"] for a in arms for lb in arms[a]]
        st["library_total_ratio"] = (max(tot) / min(tot)) if tot and min(tot) > 0 else None
        # G1 library size
        for a in arms:
            kept = []
            for lb in arms[a]:
                if lb["T"] < cm.BULK_MIN_LIB_TOTAL or lb["panel_total"] < cm.BULK_MIN_PANEL_TOTAL:
                    st["dropped"][lb["label"]] = "G1 library size"
                else:
                    kept.append(lb)
            arms[a] = kept
        # G3 a clear sex call contradicting the recorded sex drops the library
        for a in arms:
            kept = []
            for lb in arms[a]:
                if lb["sex_call"] in ("male", "female") and lb["sex_recorded"] and lb["sex_call"] != lb["sex_recorded"]:
                    st["dropped"][lb["label"]] = "G3 sex call contradicts the recorded sex"
                else:
                    kept.append(lb)
            arms[a] = kept

        def logcpm(lb):
            return np.log1p(lb["xp"][present] / lb["T"] * 1e6)

        # G4 outliers
        st["outliers"] = {}
        for a in arms:
            if not arms[a]:
                continue
            keep, recs = outlier_drop(np.vstack([logcpm(lb) for lb in arms[a]]))
            st["outliers"][a] = dict(zip([lb["label"] for lb in arms[a]], recs, strict=True))
            for lb, k in zip(arms[a], keep, strict=True):
                if not k:
                    st["dropped"][lb["label"]] = "G4 outlier"
            arms[a] = [lb for lb, k in zip(arms[a], keep, strict=True) if k]

        def prop(lb):
            v = np.where(present, lb["xp"], 0.0)
            return v / v.sum() if v.sum() > 0 else v

        # G5 single-lineage check: a control that fails is dropped; knockouts are recorded only
        st["whole_embryo"] = {}
        for a in arms:
            kept = []
            for lb in arms[a]:
                ok, rec = whole_embryo_guard(prop(lb), ref, present)
                st["whole_embryo"][lb["label"]] = {**rec, "passed": ok, "arm": a}
                if a == "control" and not ok:
                    st["dropped"][lb["label"]] = "G5 dominated by one lineage (not more like pb than any one group)"
                else:
                    kept.append(lb)
            arms[a] = kept
        m = min(len(arms["ko"]), len(arms["control"]))
        st["libraries"] = {a: [lb["label"] for lb in arms[a]] for a in arms}
        st["m"] = m
        if m < cm.BULK_MIN_LIBS:
            st["skipped"] = f"fewer than {cm.BULK_MIN_LIBS} libraries left in an arm after the guards"
            continue
        Yk = np.vstack([logcpm(lb) for lb in arms["ko"]])
        Yc = np.vstack([logcpm(lb) for lb in arms["control"]])
        ok, rec = control_grouping(Yk, Yc)
        st["grouping"] = rec
        if not ok:
            return {**uqc, "status": "refused_grouping", "reason": f"E{e:g}: the controls do not group (G4)"}
        ok, rec = whole_embryo_guard(np.mean([prop(lb) for lb in arms["control"]], 0), ref, present)
        st["whole_embryo_control_mean"] = rec
        if not ok:
            return {**uqc, "status": "refused_whole_embryo",
                    "reason": f"E{e:g}: the control arm's mean is dominated by one lineage (not more like pb than any "
                              "one group) (G5)"}
        ok, rec = somite_record([lb["somites"] for lb in arms["ko"]], [lb["somites"] for lb in arms["control"]])
        st["somites"] = rec
        if not ok:
            return {**uqc, "status": "refused_somites",
                    "reason": f"E{e:g}: the arms' somite ranges do not overlap (G6)"}
        ok, rec = composition_guard(np.vstack([prop(lb) for lb in arms["ko"]]),
                                    np.vstack([prop(lb) for lb in arms["control"]]), ref, present)
        st["composition"] = rec
        if not ok:
            return {**uqc, "status": "refused_composition",
                    "reason": f"E{e:g}: a group's share moves more than {cm.COMPOSITION_RATIO_MAX:g}x (G12)"}
        xo_k = [idm.lookup(lb["ids"], lb["x"], gene) for lb in arms["ko"]]
        xo_c = [idm.lookup(lb["ids"], lb["x"], gene) for lb in arms["control"]]
        L_own = None if any(v is None for v in xo_k + xo_c) else own_gene_ratio(
            xo_k, xo_c, [lb["T"] for lb in arms["ko"]], [lb["T"] for lb in arms["control"]])
        ok, verdict = own_gene_guard(L_own)
        st["genotype_in_counts"] = verdict
        if not ok:
            return {**uqc, "status": "refused_genotype", "reason": f"E{e:g}: {verdict} (G9)"}
        balance, srec = sex_balance([lb["sex_call"] for lb in arms["ko"]], [lb["sex_call"] for lb in arms["control"]])
        st["sex"] = {"balance": balance, **srec}
        nan_mask = always | own
        if balance != "balanced":
            if not chrom_known:
                return {**uqc, "status": "refused_sex_chromosomes",
                        "reason": f"E{e:g}: the sex check is {balance} and the gene-id map carries no chromosome, so "
                                  "the chrX panel genes cannot be set not measured (G3)"}
            nan_mask = nan_mask | chrx
        L = bulk_log_ratio(np.vstack([lb["xp"] for lb in arms["ko"]]), np.vstack([lb["xp"] for lb in arms["control"]]),
                           [lb["T"] for lb in arms["ko"]], [lb["T"] for lb in arms["control"]], present, nan_mask)
        k = cm.shrink_k_bulk(m, ds["d"])
        uqc["k_by_stage"][f"E{e:g}"] = k
        uqc["m_by_stage"][f"E{e:g}"] = m
        st.update(k=k, passed=True)
        Ls.append(k * L)
    if not Ls:
        return {**uqc, "status": "left_out", "reason": "no stage passes the guards (fail closed: left out)"}
    import process as pr

    return {**uqc, "status": "ok", "_L": pr._nanmean_rows(np.vstack(Ls))}


def write_bulk_unit(out_dir: Path, acc: str, gene: str, L: np.ndarray, ref: dict, meta: dict) -> dict:
    """unit.json (L under '__whole__', modality bulk) and the finished per-dataset table (pooled form only, against
    the reference's whole-embryo pseudo-bulk)."""
    name = f"{acc}__{gene}".replace("+", "_and_").replace("/", "_")
    out_dir = cm.outside_repo(out_dir)
    unit = {"accession": acc, "gene": gene, "genes": ref["genes"], "refstats_sha256": ref["file_sha256"],
            "modality": "bulk", "L": {WHOLE: [None if not np.isfinite(v) else float(v) for v in L]}, **meta}
    up = out_dir / f"{name}.unit.json"
    cm.dump_json(unit, up)
    phi = tb.finish(L, tb.whole_den(ref))
    tables = tb.write_table(out_dir / "tables" / name, name, ref["genes"], phi, None,
                            {**meta, "modality": "bulk", "per_group": "none (bulk)"})
    return {"unit_file": str(up), "unit_sha256": cm.sha256_file(up), "tables": tables}
