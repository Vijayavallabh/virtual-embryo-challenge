#!/usr/bin/env python
"""The many-line ArrayExpress resource of the bulk units (stage 9, phase 2; README "The many-line resource").

E-MTAB-7322 only (AE_ACCESSIONS). Every url is derived from the accession, on one host (fetch.AE_HOSTS), with no query
and no port:
  IDF and SDRF   https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/<last 3 digits>/<ACC>/Files/<ACC>.idf.txt, .sdrf.txt
  count matrix   https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/<ACC>/<ACC>-raw-counts.tsv.undecorated
  methods text   .../atlas/experiments/<ACC>/<ACC>-analysis-methods.tsv

plan downloads the IDF, the SDRF and the methods text only (metadata; no count). Text that is not UTF-8, a line or
paragraph separator inside a line, an IDF key or SDRF header outside the closed grammar, or a row of another width
refuses the dataset with its text withheld. Then, BEFORE the screen, the views, the hashes and the tables, every SDRF
row whose line column ('Characteristics[mutant line]') equals a withheld_lines value exactly is dropped, with every
other row of a Source Name, ENA sample or individual that a withheld row names; only '<gene>: n rows withheld
(reason)' is printed. A line value that differs from a withheld value but equals it after case folding and removing
spaces, '-', '_' and '.' refuses the dataset. A sample is one Source Name: its fields map '<3-digit column number>|
<header>' to the list of its rows' values in file order, hashed with t3x_common.canonical_metadata(Source Name, ...).
The series block is the IDF (each key under 'idf|<key>'), the SDRF header row and the methods text, hashed alike and
pinned as pinned_series. The screen runs on every series field and name and on every Source Name's record; a held-out
hit anywhere withholds the text.

The mechanical comparisons, each printing its number of disagreements (a disagreeing Source Name is NOT PINNABLE):
  (a) all rows of a Source Name agree on every Characteristics and Factor Value column;
  (b) the Source Name's genotype token (source_name_pattern's group 'genotype') maps through
      t3x_common.GENOTYPE_MAP to Characteristics[genotype], which is one of the map's values;
  (c) the Source Name's gene token equals the line equals the allele's prefix (Comment[allele] before its first '<');
  (d) the stage is 'embryonic day X' exactly, the somite value fullmatches '^(\\d+)-somite stage$', and no
      Characteristics or Factor Value has leading or trailing whitespace;
  (e) exactly one ENA sample per Source Name, one Source Name per individual and per ENA sample, and no Source Name
      or ENA sample under two lines;
  (f) every run is an ERR accession of its own record (its FASTQ_URI names it) and of no other Source Name;
  (g) the pin hash equals the plan hash (verify_pins).
Only 'embryonic day 9.5' is pinnable (days 10 and 10.5 fall under the stage cap); hets and the line 'baseline' are
never pinned. The committed rule fixes the line set: every line not withheld with at least BULK_MIN_LIBS homozygous and
BULK_MIN_LIBS wild-type pinnable E9.5 Source Names; verify_pins refuses pins that are not exactly that set.
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import bulk  # noqa: E402
import fetch as fe  # noqa: E402
import t3x_common as cm  # noqa: E402
import t3x_screen as sc  # noqa: E402

Refused = cm.Refused
AE_ACCESSIONS = ("E-MTAB-7322",)
COL = {"source": "Source Name", "ers": "Comment[ENA_SAMPLE]", "line": "Characteristics[mutant line]",
       "genotype": "Characteristics[genotype]", "stage": "Characteristics[developmental stage]",
       "somite": "Characteristics[somite stage]", "sex": "Characteristics[sex]",
       "individual": "Characteristics[individual]", "allele": "Comment[allele]", "run": "Comment[ENA_RUN]"}
FASTQ = "Comment[FASTQ_URI]"
PINNABLE_STAGE = "embryonic day 9.5"
BASELINE = "baseline"
_BR = r"\[[A-Za-z][A-Za-z0-9_ ]*\]"
SDRF_HEADER = re.compile(
    r"(?:Source Name|Sample Name|Extract Name|Labeled Extract Name|Label|Assay Name|Scan Name|Hybridization Name|"
    r"Material Type|Protocol REF|Term Source REF|Term Accession Number|Performer|Technology Type|Array Design REF|"
    r"Array Data File|Derived Array Data File|Normalization Name|"
    rf"Comment{_BR}|Characteristics{_BR}|Factor Value{_BR}|Unit{_BR}|Parameter Value{_BR})")
IDF_KEY = re.compile(rf"[A-Za-z][A-Za-z0-9 ]*(?:{_BR})?")
STAGE_VALUE = re.compile(r"embryonic day (\d+(?:\.\d+)?)")
SOMITE_VALUE = re.compile(r"^(\d+)-somite stage$")
ERR_RX = re.compile(r"ERR\d{6,7}")
VIEW_MAX_DISTINCT = 60
WITHHELD = "<withheld>"


# ---------------------------------------------------------------------------------------------------------------------
# Urls (derived from the accession only)


def _check_acc(acc: str) -> None:
    if acc not in AE_ACCESSIONS:
        raise Refused(f"{acc}: the arrayexpress form takes {AE_ACCESSIONS} only")


def biostudies_folder(acc: str) -> str:
    _check_acc(acc)
    m = re.fullmatch(r"(E-[A-Z]+-)(\d+)", acc)
    return f"https://ftp.ebi.ac.uk/biostudies/fire/{m.group(1)}/{m.group(2)[-3:].zfill(3)}/{acc}/Files"


def atlas_folder(acc: str) -> str:
    _check_acc(acc)
    return f"https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/{acc}"


def urls(acc: str) -> dict[str, str]:
    return {"idf": f"{biostudies_folder(acc)}/{acc}.idf.txt", "sdrf": f"{biostudies_folder(acc)}/{acc}.sdrf.txt",
            "methods": f"{atlas_folder(acc)}/{acc}-analysis-methods.tsv",
            "matrix": f"{atlas_folder(acc)}/{acc}-raw-counts.tsv.undecorated"}


def matrix_name(acc: str) -> str:
    return f"{acc}-raw-counts.tsv.undecorated"


# ---------------------------------------------------------------------------------------------------------------------
# Parsing (no network)


class Withhold(Refused):
    """A refusal whose text must not be written."""


def text_lines(raw: bytes, label: str) -> list[str]:
    """Strict UTF-8, split on '\\n' only, one trailing CR stripped per line, a final empty line dropped; any other
    line or paragraph separator inside a line refuses (text withheld)."""
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        raise Withhold(f"{label}: undecodable bytes: refused, text withheld") from None
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines = lines[:-1]
    out = []
    for i, ln in enumerate(lines, 1):
        if ln.endswith("\r"):
            ln = ln[:-1]
        if any(c in fe.LINE_SEPARATORS for c in ln):
            raise Withhold(f"{label}: a line or paragraph separator inside line {i}: refused, text withheld")
        out.append(ln)
    return out


def parse_idf(raw: bytes) -> dict[str, list[str]]:
    """{'idf|<key>': [values]} in file order (a repeated key's values appended), '#' lines and other lines under
    '_other_lines'; a key outside the closed grammar refuses (text withheld)."""
    out: dict[str, list[str]] = {}
    for i, ln in enumerate(text_lines(raw, "IDF"), 1):
        if not ln.strip():
            continue
        if ln.startswith("#"):
            out.setdefault(fe.OTHER, []).append(ln)
            continue
        f = ln.split("\t")
        if not IDF_KEY.fullmatch(f[0]):
            raise Withhold(f"IDF line {i}: a key outside the closed grammar: refused, text withheld")
        out.setdefault(f"idf|{f[0]}", []).extend(f[1:])
    return out


def parse_sdrf(raw: bytes) -> tuple[list[str], list[list[str]]]:
    """(header, rows): every header in the closed grammar, every row of the header's width; values verbatim (no
    whitespace stripped: comparison (d) reads it)."""
    lines = text_lines(raw, "SDRF")
    if not lines:
        raise Withhold("SDRF: empty: refused")
    header = lines[0].split("\t")
    bad = [i for i, h in enumerate(header) if not SDRF_HEADER.fullmatch(h)]
    if bad:
        raise Withhold(f"SDRF: {len(bad)} header(s) outside the closed grammar: refused, text withheld")
    rows = []
    for i, ln in enumerate(lines[1:], 2):
        if not ln.strip():
            continue
        r = ln.split("\t")
        if len(r) != len(header):
            raise Withhold(f"SDRF line {i}: {len(r)} fields, not {len(header)}: refused, text withheld")
        rows.append(r)
    for k, name in COL.items():
        if header.count(name) != 1:
            raise Withhold(f"SDRF: the column {name} ({k}) appears {header.count(name)} times: refused")
    return header, rows


def _norm(v: str) -> str:
    return re.sub(r"[\s\-_.]", "", v).casefold()


def withhold_rows(ds: dict, header: list[str], rows: list[list[str]]) -> tuple[list[list[str]], dict, int]:
    """(the rows kept, {gene: (n rows withheld, reason)}, n rows withheld because their Source Name, ENA sample or
    individual is named by a withheld row). A near-match of a withheld value refuses."""
    wl = {e["line"]: e["reason"] for e in ds.get("withheld_lines") or []}
    li = header.index(COL["line"])
    near = sorted({w for v in {r[li] for r in rows} for w in wl if v != w and _norm(v) == _norm(w)})
    if near:
        raise Withhold(f"a line value near the withheld line(s) {near} (a case, space or separator variant): the "
                       "dataset is refused")
    counts = {w: [0, wl[w]] for w in wl}
    linked_keys = set()
    keep = []
    for r in rows:
        if r[li] in wl:
            counts[r[li]][0] += 1
            for c in ("source", "ers", "individual"):
                linked_keys.add((c, r[header.index(COL[c])]))
        else:
            keep.append(r)
    kept, linked = [], 0
    for r in keep:
        if any((c, r[header.index(COL[c])]) in linked_keys for c in ("source", "ers", "individual")):
            linked += 1
        else:
            kept.append(r)
    return kept, {w: tuple(v) for w, v in counts.items()}, linked


def records(header: list[str], rows: list[list[str]]) -> dict[str, dict]:
    """{Source Name: fields} with fields '<3-digit column number>|<header>' -> [values of its rows], file order."""
    si = header.index(COL["source"])
    out: dict[str, dict] = {}
    for r in rows:
        f = out.setdefault(r[si], {})
        for i, (h, v) in enumerate(zip(header, r, strict=True)):
            f.setdefault(f"{i + 1:03d}|{h}", []).append(v)
    return out


def _vals(fields: dict, header: list[str], name: str) -> list[str]:
    return fields.get(f"{header.index(name) + 1:03d}|{name}", [])


def _one(vals: list[str]) -> str | None:
    s = list(dict.fromkeys(vals))
    return s[0] if len(s) == 1 else None


def comparisons(ds: dict, header: list[str], recs: dict[str, dict]) -> dict[str, dict[str, str]]:
    """{Source Name: {letter: reason}} of the mechanical comparisons (a) to (f)."""
    pat = re.compile(ds["source_name_pattern"])
    out: dict[str, dict[str, str]] = {sn: {} for sn in recs}
    ers_of, ind_of, run_of = {}, {}, {}
    for sn, f in recs.items():
        for k, vals in f.items():
            h = k.split("|", 1)[1]
            if h.startswith(("Characteristics[", "Factor Value[")):
                if len(set(vals)) > 1:
                    out[sn].setdefault("a", f"its rows disagree on {h}")
                if any(v != v.strip() for v in vals):
                    out[sn].setdefault("d", f"a value of {h} with leading or trailing whitespace")
        line = _one(_vals(f, header, COL["line"]))
        geno = _one(_vals(f, header, COL["genotype"]))
        if line != BASELINE:
            m = pat.fullmatch(sn)
            if geno not in cm.GENOTYPE_MAP.values():
                out[sn].setdefault("b", "its genotype value is not one of the committed map's values")
            elif m is None:
                out[sn].setdefault("b", "the Source Name does not match source_name_pattern")
            elif cm.GENOTYPE_MAP.get(m.group("genotype")) != geno:
                out[sn].setdefault("b", "its genotype token does not map to its genotype value")
            allele = _one(_vals(f, header, COL["allele"]))
            prefix = allele.split("<", 1)[0] if allele else None
            if m is None or line is None or not prefix or not (m.group("gene") == line == prefix):
                out[sn].setdefault("c", "the Source Name's gene token, the line and the allele prefix differ")
        stage = _one(_vals(f, header, COL["stage"]))
        if stage is None or not STAGE_VALUE.fullmatch(stage):
            out[sn].setdefault("d", "its stage is not one value 'embryonic day X'")
        som = _one(_vals(f, header, COL["somite"]))
        if line != BASELINE and (som is None or not SOMITE_VALUE.fullmatch(som)):
            out[sn].setdefault("d", "its somite value is not one value '<n>-somite stage'")
        ers = list(dict.fromkeys(_vals(f, header, COL["ers"])))
        if len(ers) != 1:
            out[sn].setdefault("e", f"{len(ers)} ENA samples")
        for e in ers:
            ers_of.setdefault(e, set()).add(sn)
        for i in set(_vals(f, header, COL["individual"])):
            ind_of.setdefault(i, set()).add(sn)
        runs = _vals(f, header, COL["run"])
        fq = _vals(f, header, FASTQ) if FASTQ in header else None
        for k, r in enumerate(runs):
            if not ERR_RX.fullmatch(r):
                out[sn].setdefault("f", "a run that is not an ERR accession")
            elif fq is not None and r not in fq[k]:
                out[sn].setdefault("f", "a row's FASTQ_URI does not name its run")
            run_of.setdefault(r, set()).add(sn)
    lines_of_sn = {sn: set(_vals(f, header, COL["line"])) for sn, f in recs.items()}
    for sns in ers_of.values():
        lines = set().union(*(lines_of_sn[s] for s in sns))
        if len(sns) > 1 or len(lines) > 1:
            for s in sns:
                out[s].setdefault("e", "an ENA sample under two Source Names or two lines")
    for sns in ind_of.values():
        if len(sns) > 1:
            for s in sns:
                out[s].setdefault("e", "an individual under two Source Names")
    for sn, ls in lines_of_sn.items():
        if len(ls) > 1:
            out[sn].setdefault("e", "a Source Name under two lines")
    for sns in run_of.values():
        if len(sns) > 1:
            for s in sns:
                out[s].setdefault("f", "a run under two Source Names")
    return out


def stage_of(value: str | None) -> tuple[float | None, str]:
    m = STAGE_VALUE.fullmatch(value) if isinstance(value, str) else None
    if not m:
        return None, "no single stage value 'embryonic day X'"
    e = float(m.group(1))
    try:
        cm.stage_ok(f"E{e:g}")
    except Refused as err:
        return None, str(err)
    if value != PINNABLE_STAGE:
        return e, f"only '{PINNABLE_STAGE}' is pinnable"
    return e, "ok"


def line_set(rows: list[dict]) -> dict[str, dict[str, list[str]]]:
    """The committed rule: every line (withheld lines were dropped before) other than 'baseline' with at least
    BULK_MIN_LIBS pinnable homozygous and BULK_MIN_LIBS pinnable wild-type Source Names at E9.5."""
    by: dict[str, dict[str, list[str]]] = {}
    for r in rows:
        if r["not_pinnable"] or r["line"] in (None, BASELINE):
            continue
        arm = {"homozygous knockout": "ko", "wild type genotype": "control"}.get(r["genotype"])
        if arm:
            by.setdefault(r["line"], {"ko": [], "control": []})[arm].append(r["source_name"])
    return {k: v for k, v in sorted(by.items()) if len(v["ko"]) >= cm.BULK_MIN_LIBS
            and len(v["control"]) >= cm.BULK_MIN_LIBS}


def plan_from_texts(ds: dict, idf: bytes, sdrf: bytes, methods: bytes, headers: dict | None = None) -> dict:
    """The plan of the resource from its three texts (no network)."""
    acc = ds["accession"]
    _check_acc(acc)
    out = {"accession": acc, "own_series": acc, "repository": "arrayexpress", "bulk": bulk.declaration(ds),
           "screen": {"n_samples": 0}, "tables": {}, "withheld": True, "count_table_headers": headers or {}}
    try:
        sfields = parse_idf(idf)
        header, rows = parse_sdrf(sdrf)
        mlines = text_lines(methods, "methods")
        rows, wcounts, linked = withhold_rows(ds, header, rows)
    except Withhold as e:
        return {**out, "status": "refused", "reason": str(e)}
    out["withheld_rows"] = {g: {"rows": n, "reason": why} for g, (n, why) in wcounts.items()}
    out["withheld_linked_rows"] = linked
    sfields["sdrf|header"] = list(header)
    sfields["methods|lines"] = mlines
    recs = records(header, rows)
    held = {h.upper() for h in fe.held_set(ds)}
    named = sorted({m.group(0).upper() for t in [*(v for vs in sfields.values() for v in vs),
                                                   *(v for r in rows for v in r)]
                    for m in fe._ACCESSION.finditer(t)} & held)
    # the screen: every series field and name, every Source Name's record
    extra, cleared = ds.get("screen_extra_refuse"), ds.get("review_cleared")
    texts = [(f"{k} (field name)", k) for k in sfields] + [(k, v) for k, vs in sfields.items() for v in vs]
    void = any(sc.extended(t) for _, t in texts)
    refusals, held_hit = [], False
    for key, t in texts:
        for kind, term in sc.screen_text(acc, t, extra, cleared, void=void):
            held_hit |= kind == "held_out"
            refusals.append(f"{acc} {key}: {kind}" + (f" term {term!r}" if kind == "review" else ""))
    scr = {}
    for sn, f in recs.items():
        r = sc.screen_record(sn, f, extra, cleared)
        scr[sn] = r
        if r["verdict"] != "pass":
            held_hit |= r["verdict"] == "held_out"
            terms = sorted({h[2] for h in r["hits"] if h[1] == "review"})
            refusals.append(f"{sn}: {r['verdict']}" + (f" (review terms {terms})" if terms else ""))
    out["screen"]["n_samples"] = len(recs)
    if held_hit or named:
        out["tables"] = {acc: [{"source_name": WITHHELD, "verdict": "held_out" if held_hit else "held accession"}]}
        out["screen"]["refusals"] = [fe.WITHHELD]
        return {**out, "status": "refused",
                "reason": "screen: a held-out hit or a held accession in the resource (text withheld)"}
    out["withheld"] = False
    out["series"] = {acc: {"accession": acc, "fields": sfields, "role": "screened series",
                           "series_metadata_sha256": cm.series_metadata_sha256(acc, sfields)}}
    comp = comparisons(ds, header, recs)
    rows_out = []
    for sn, f in recs.items():
        why = dict(comp[sn])
        line = _one(_vals(f, header, COL["line"]))
        geno = _one(_vals(f, header, COL["genotype"]))
        stage_v = _one(_vals(f, header, COL["stage"]))
        e, how = stage_of(stage_v)
        reasons = [f"({k}) {v}" for k, v in sorted(why.items())]
        if scr[sn]["verdict"] != "pass":
            reasons.append(f"screen verdict {scr[sn]['verdict']}")
        if line == BASELINE:
            reasons.append("the line 'baseline': never a knockout line or a control")
        if geno == cm.GENOTYPE_MAP["het"]:
            reasons.append("heterozygous: never pinned in this round")
        if how != "ok":
            reasons.append(f"stage: {how}")
        ers = list(dict.fromkeys(_vals(f, header, COL["ers"])))
        rows_out.append({"source_name": sn, "accession": ers[0] if len(ers) == 1 else None, "fields": f,
                         "metadata_sha256": cm.metadata_sha256(sn, f), "verdict": scr[sn]["verdict"],
                         "hits": [list(h) for h in scr[sn]["hits"]], "clearance_void": scr[sn]["clearance_void"],
                         "comparisons": why, "not_pinnable": "; ".join(reasons) or None, "line": line,
                         "genotype": geno, "stage": e, "stage_value": stage_v,
                         "somite": _one(_vals(f, header, COL["somite"])), "sex": _one(_vals(f, header, COL["sex"])),
                         "individual": _one(_vals(f, header, COL["individual"])),
                         "runs": list(dict.fromkeys(_vals(f, header, COL["run"])))})
    out["tables"] = {acc: rows_out}
    out["comparisons"] = {k: sorted(sn for sn, c in comp.items() if k in c) for k in "abcdef"}
    out["line_set"] = line_set(rows_out)
    out["views"] = views(header, recs, rows_out)
    out["screen"]["refusals"] = refusals
    if refusals:
        return {**out, "status": "refused", "reason": "screen: " + "; ".join(refusals[:10])}
    if ds.get("pinned"):
        try:
            out["pinned"] = [{k: v for k, v in s.items() if k != "row"} for s in verify_pins(ds, out, live=True)]
            out["units"] = cm.pinned_units(ds["pinned"])
        except Refused as e:
            return {**out, "status": "refused", "reason": f"pins: {e}"}
    return {**out, "status": "ok"}


def views(header: list[str], recs: dict, rows: list[dict]) -> dict:
    """The distinct-value view (per column, value -> number of Source Names; past VIEW_MAX_DISTINCT values only their
    number and the pattern view, digits as '#'), and the per-line view (genotype by stage, sex, somite range)."""
    cols = {}
    for i, h in enumerate(header):
        key = f"{i + 1:03d}|{h}"
        vals: dict[str, int] = {}
        for f in recs.values():
            for v in dict.fromkeys(f.get(key, [])):
                vals[v] = vals.get(v, 0) + 1
        pats: dict[str, int] = {}
        for v, n in vals.items():
            p = re.sub(r"[0-9]+", "#", v)
            pats[p] = pats.get(p, 0) + n
        cols[key] = {"distinct": len(vals), "patterns": dict(sorted(pats.items())),
                     **({"values": dict(sorted(vals.items()))} if len(vals) <= VIEW_MAX_DISTINCT else {})}
    per_line: dict[str, dict] = {}
    for r in rows:
        v = per_line.setdefault(str(r["line"]), {"genotype_by_stage": {}, "sex": {}, "somites": [],
                                                 "pinnable": {"ko": 0, "control": 0}})
        k = f"{r['genotype']} | {r['stage_value']}"
        v["genotype_by_stage"][k] = v["genotype_by_stage"].get(k, 0) + 1
        v["sex"][str(r["sex"])] = v["sex"].get(str(r["sex"]), 0) + 1
        s = bulk.somites(r["somite"]) if r["somite"] else None
        if s is not None:
            v["somites"].append(s)
        arm = {"homozygous knockout": "ko", "wild type genotype": "control"}.get(r["genotype"])
        if arm and not r["not_pinnable"]:
            v["pinnable"][arm] += 1
    for v in per_line.values():
        v["somite_range"] = [min(v["somites"]), max(v["somites"])] if v["somites"] else None
        del v["somites"]
    return {"columns": cols, "per_line": dict(sorted(per_line.items()))}


# ---------------------------------------------------------------------------------------------------------------------
# The plan step (network: three small text files)


def plan(acc: str, dest: Path, datasets: Path | None = None, texts: dict | None = None) -> dict:
    """Download the IDF, the SDRF and the methods text (or take texts={'idf','sdrf','methods'} bytes: tests), plan,
    write plan.json and the plan table. plan.json says in_progress from the start."""
    ds = cm.dataset(acc, datasets)
    _check_acc(acc)
    dest = cm.outside_repo(dest)
    mdir = dest / "meta" / acc
    cm.dump_json({"accession": acc, "status": "in_progress", "reason": "a plan run started and did not finish"},
                 mdir / "plan.json")
    files = {}
    try:
        if texts is None:
            texts = {}
            for k in ("idf", "sdrf", "methods"):
                u = urls(acc)[k]
                p = mdir / Path(u).name
                size, digest = fe._get_ae(u, p, 256 * 1024 ** 2)
                files[k] = {"url": u, "file": str(p), "size": size, "sha256": digest}
                texts[k] = p.read_bytes()
    except Refused as e:
        pl = {"accession": acc, "own_series": acc, "status": "refused", "reason": str(e), "withheld": True,
              "tables": {}, "screen": {"n_samples": 0, "refusals": [str(e)]}, "files": files,
              "datasets_sha256": cm.datasets_sha256(datasets), "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
              "table_files": []}
        cm.dump_json(pl, mdir / "plan.json")
        return pl
    pl = plan_from_texts(ds, texts["idf"], texts["sdrf"], texts["methods"], bulk.stored_headers(dest, acc))
    pl.update(datasets_sha256=cm.datasets_sha256(datasets), files=files, d=ds["d"], tier=ds["tier"],
              study=ds["study"], created=time.strftime("%Y-%m-%d %H:%M:%S %Z"))
    pl["table_files"] = write_tables(pl, mdir)
    cm.dump_json(pl, mdir / "plan.json")
    return pl


def write_tables(pl: dict, mdir: Path) -> list[str]:
    """The one plan table of the resource: the withheld counts, the series block, the views, the comparisons, the
    line set and every Source Name's record verbatim."""
    acc = pl["accession"]
    L = [f"# plan table {acc} (ArrayExpress, bulk): plan status {pl['status']}; written {pl.get('created', '')}",
         "# Values are verbatim from the IDF, the SDRF and the methods text (nothing stripped).",
         "# A value that is not plain printable ASCII is followed by an 'escaped:' line, its Python ascii()."]
    for g, v in (pl.get("withheld_rows") or {}).items():
        L.append(f"# {g}: {v['rows']} rows withheld ({v['reason']})")
    if "withheld_linked_rows" in pl:
        L.append(f"# {pl['withheld_linked_rows']} more rows withheld: their Source Name, ENA sample or individual is "
                 "named by a withheld row")
    if pl.get("withheld"):
        L.append(f"# WITHHELD ({pl.get('reason')}): no text of this resource is written. Do not pin.")
    else:
        L += bulk.plan_lines(pl)
        L.append(f"# status reason: {pl.get('reason', 'none')}")
        b = pl["series"][acc]
        L += ["", f"==== SERIES {acc}  series_metadata_sha256 {b['series_metadata_sha256']}  (IDF, SDRF header row, "
                  "methods text)", *fe.field_lines(b["fields"])]
        L += ["", "# VIEWS: distinct values per column (value -> Source Names), patterns (digits as #)"]
        for k, v in pl["views"]["columns"].items():
            L.append(f"  {k}: {v['distinct']} distinct")
            for val, n in (v.get("values") or {}).items():
                shown = val if fe.plain(val) else ascii(val)
                L.append(f"    {shown} -> {n}")
            if "values" not in v:
                for p, n in v["patterns"].items():
                    L.append(f"    pattern {p if fe.plain(p) else ascii(p)} -> {n}")
        L += ["", "# PER LINE: genotype | stage counts, sex, somite range, pinnable E9.5 Source Names"]
        for line, v in pl["views"]["per_line"].items():
            L.append(f"  {line}: {json.dumps(v, sort_keys=True)}")
        L += ["", "# MECHANICAL COMPARISONS (a disagreeing Source Name is NOT PINNABLE)"]
        for k, sns in pl["comparisons"].items():
            L.append(f"  ({k}) {len(sns)} disagreements" + (f": {sns}" if sns else ""))
        L += ["", f"# LINE SET (the committed rule: lines not withheld with at least {cm.BULK_MIN_LIBS} homozygous and "
                  f"{cm.BULK_MIN_LIBS} wild-type pinnable Source Names at '{PINNABLE_STAGE}'; no hand-picking)"]
        for line, v in pl["line_set"].items():
            L.append(f"  {line}: ko {v['ko']}; control {v['control']}")
        L += ["", "# SOURCE NAMES"]
        for r in pl["tables"][acc]:
            L += ["", f"==== {r['source_name']}  ENA sample {r['accession']}  verdict {r['verdict']}  metadata_sha256 "
                      f"{r['metadata_sha256']}",
                  f"  line {r['line']}; genotype {r['genotype']}; stage {r['stage_value']}; runs {r['runs']}",
                  f"  NOT PINNABLE: {r['not_pinnable']}" if r["not_pinnable"] else "  pinnable",
                  f"  screen hits: {r['hits'] or 'none'}", *fe.field_lines(r["fields"])]
    p = cm.outside_repo(mdir / f"plan_table_{acc}.txt")
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(L) + "\n", encoding="utf-8")
    return [str(p)]


# ---------------------------------------------------------------------------------------------------------------------
# The header step and the run ties


def run_ties(fields: list[str], rows: list[dict], gene_columns: list[str]) -> tuple[dict, list[str]]:
    """The ties of the matrix header (exact equality with an ERR run, an ENA sample or a Source Name of a record
    kept): ({'by_header': {header: Source Name}, 'untied': n, 'complete': [...], 'incomplete': [...],
    'other_tie': [...]}, problems). A header tying to two records, a header appearing twice and a gene column not
    appearing once are problems (the table is refused). A header tying to no kept record is never read (counted
    only). A record whose runs are not all headers is incomplete, one tied through a non-run header is 'other_tie':
    neither is pinnable."""
    keys: dict[str, set] = {}
    for r in rows:
        for k in [*r["runs"], r["accession"], r["source_name"]]:
            if k:
                keys.setdefault(k, set()).add(r["source_name"])
    problems = []
    dup = sorted({h for h in fields if fields.count(h) > 1})
    if dup:
        problems.append(f"{len(dup)} headers appear twice")
    for g in gene_columns:
        if fields.count(g) != 1:
            problems.append(f"the gene column {g!r} appears {fields.count(g)} times")
    by, untied = {}, 0
    for h in fields:
        if h in gene_columns:
            continue
        sns = keys.get(h, set())
        if len(sns) > 1:
            problems.append(f"a header ties to {len(sns)} Source Names")
        elif sns:
            by[h] = next(iter(sns))
        else:
            untied += 1
    complete, incomplete, other = [], [], []
    run_set = {h for h in by if ERR_RX.fullmatch(h)}
    for r in rows:
        tied = [h for h, sn in by.items() if sn == r["source_name"]]
        if any(h not in r["runs"] for h in tied):
            other.append(r["source_name"])
        elif r["runs"] and all(x in run_set for x in r["runs"]):
            complete.append(r["source_name"])
        else:
            incomplete.append(r["source_name"])
    return {"by_header": by, "untied": untied, "complete": complete, "incomplete": incomplete,
            "other_tie": other}, problems


def header(acc: str, dest: Path, datasets: Path | None = None) -> dict:
    """The header step of the resource: from the stored plan (ok), the matrix url derived from the accession, its
    first line streamed (fetch._open), screened and tied by run; only the record is written (the line is kept in the
    record for the re-derivation, never printed; untied headers are counted, never listed)."""
    ds = cm.dataset(acc, datasets)
    _check_acc(acc)
    dest = cm.outside_repo(dest)
    pl = fe._load_plan(acc, dest, datasets)
    if pl.get("withheld") or acc not in (pl.get("series") or {}):
        raise Refused(f"{acc}: the stored plan holds no readable record")
    name, url = matrix_name(acc), urls(acc)["matrix"]
    with fe._open(url) as r:
        line = bulk.first_line(r.read, False, name)
    return record_header(ds, dest, pl, name, url, line)


def record_header(ds: dict, dest: Path, pl: dict, name: str, url: str, line: bytes) -> dict:
    acc = ds["accession"]
    fields = bulk.split_header(line, name)
    text = line.decode("utf-8")
    rec = {"accession": acc, "series": acc, "name": name, "url": url, "n_columns": len(fields),
           "header_sha256": bulk.header_sha256(line), "gene_columns": ds["gene_columns"],
           "series_metadata_sha256": pl["series"][acc]["series_metadata_sha256"],
           "created": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
    held, bad = bulk.screen_header(ds, acc, fields, text)
    hp = bulk.header_path(dest, acc, name)
    if held:
        rec.update(status="refused", withheld=True, reason="a held-out hit in the header (text withheld)")
        cm.dump_json(rec, hp)
        return rec
    ties, problems = run_ties(fields, pl["tables"][acc], ds["gene_columns"])
    rec.update(withheld=False, text=text, screen=bad, ties=ties, tie_problems=problems)
    if bad:
        rec.update(status="refused", reason=f"header screen: {bad}")
    elif problems:
        rec.update(status="refused", reason="tie rule: " + "; ".join(problems[:10]))
    else:
        rec["status"] = "ok"
    cm.dump_json(rec, hp)
    return rec


def header_lines(rec: dict) -> list[str]:
    out = [f"header {rec['name']}: status {rec['status']}; {rec['n_columns']} columns; header_sha256 "
           f"{rec['header_sha256']}" + (f"; {rec.get('reason')}" if rec.get("reason") else "")]
    if rec.get("withheld"):
        return out
    t = rec["ties"]
    out += [f"  headers tied to a kept Source Name: {len(t['by_header'])}; headers tied to none (never read): "
            f"{t['untied']}",
            f"  Source Names with every run a header: {len(t['complete'])}; with a run missing (not pinnable): "
            f"{len(t['incomplete'])} {t['incomplete']}; tied through a non-run header (not pinnable): {t['other_tie']}"]
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Pins


def verify_pins(ds: dict, pl: dict, live: bool = False) -> list[dict]:
    """The pins of the resource checked against the plan: the plan made under this declaration; the series pin; the
    count table's stored header (its pinned hash, the url derived from the accession, its run ties derived again);
    each pin's Source Name kept, pinnable, its ENA sample, its runs (all of them, every one a header tied to it), its
    stage E9.5, its arm from its genotype, its gene its line, its hash, the screen again (comparison (g)); and the
    pins are exactly the committed rule's line set, every eligible Source Name of each line and arm."""
    acc = ds["accession"]
    pins = ds.get("pinned") or []
    if not pins:
        raise Refused(f"{acc}: the pinned block is empty: nothing is fetched or processed")
    if pl.get("withheld") or acc not in (pl.get("tables") or {}):
        raise Refused(f"{acc}: the plan holds no readable table")
    if (pl.get("bulk") or None) != bulk.declaration(ds):
        raise Refused(f"{acc}: the plan was made under another declaration (withheld_lines, source_name_pattern ...): "
                      "plan again")
    what = "live" if live else "stored"
    sp = cm.check_bulk_series_pins(acc, ds, True)
    b = pl["series"][acc]
    h = cm.series_metadata_sha256(acc, b["fields"])
    if h != sp[0]["series_metadata_sha256"] or h != b["series_metadata_sha256"]:
        raise Refused(f"{acc}: its {what} IDF, SDRF header and methods text hash to {h}, not to the pinned value: read "
                      "and pin again")
    t = ds["pinned_series_files"][0]
    rec = (pl.get("count_table_headers") or {}).get(t["name"])
    if not rec or rec.get("status") != "ok" or rec.get("withheld"):
        raise Refused(f"{t['name']}: no passing header record in the plan: run the header step, then plan again")
    line = rec["text"].encode("utf-8")
    if (bulk.header_sha256(line) != t["header_sha256"] or rec["header_sha256"] != t["header_sha256"]
            or rec["url"] != urls(acc)["matrix"]):
        raise Refused(f"{t['name']}: the stored header is not the pinned one (hash or url): refused")
    rows = pl["tables"][acc]
    ties, problems = run_ties(bulk.split_header(line, t["name"]), rows, ds["gene_columns"])
    if problems or ties != rec["ties"]:
        raise Refused(f"{t['name']}: the run ties derived again differ from the stored ones: refused")
    by_sn = {r["source_name"]: r for r in rows}
    ls = line_set(rows)
    per: dict[str, dict] = {}
    for p in pins:
        sn = p["source_name"]
        r = by_sn.get(sn)
        if r is None:
            raise Refused(f"{sn}: not a kept Source Name of the plan (withheld, or unknown)")
        if r["not_pinnable"]:
            raise Refused(f"{sn}: not pinnable: {r['not_pinnable']}")
        if r["accession"] != p["accession"]:
            raise Refused(f"{sn}: pinned with ENA sample {p['accession']}, its record holds {r['accession']}")
        hh = cm.metadata_sha256(sn, r["fields"])
        if hh != p["metadata_sha256"] or hh != r["metadata_sha256"]:
            raise Refused(f"{sn}: its {what} record hashes to {hh}, not to the pinned {p['metadata_sha256']} "
                          "(comparison (g)): read and pin again")
        s = sc.screen_record(sn, r["fields"], ds.get("screen_extra_refuse"), ds.get("review_cleared"))
        if s["verdict"] != "pass":
            raise Refused(f"{sn}: screen verdict {s['verdict']}: refused")
        if p["runs"] != r["runs"] or sn not in ties["complete"]:
            raise Refused(f"{sn}: pinned runs {p['runs']} are not exactly its record's runs, each a header tied to it")
        if r["stage"] is None or r["stage_value"] != PINNABLE_STAGE or cm.stage_ok(p["stage"]) != r["stage"]:
            raise Refused(f"{sn}: pinned stage {p['stage']}, its record says {r['stage_value']!r}")
        want = {"ko": cm.GENOTYPE_MAP["hom"], "control": cm.GENOTYPE_MAP["wt"]}[p["arm"]]
        if r["genotype"] != want:
            raise Refused(f"{sn}: pinned as {p['arm']}, its genotype is {r['genotype']!r}")
        if p["gene"] != r["line"]:
            raise Refused(f"{sn}: pinned for gene {p['gene']}, its line is {r['line']}")
        x = per.setdefault(sn, {"accession": p["accession"], "source_name": sn, "arm": p["arm"], "genes": [],
                                "stage": r["stage"], "runs": list(r["runs"]), "metadata_sha256": hh, "row": r,
                                "sex": bulk.recorded_sex(r["sex"]), "somite": r["somite"]})
        x["genes"].append(p["gene"])
    pinned_sets = {}
    for p in pins:
        pinned_sets.setdefault(p["gene"], {"ko": set(), "control": set()})[p["arm"]].add(p["source_name"])
    if set(pinned_sets) != set(ls):
        raise Refused(f"the pinned lines {sorted(pinned_sets)} are not the committed rule's line set {sorted(ls)}: no "
                      "hand-picking")
    for g, v in ls.items():
        for arm in ("ko", "control"):
            if pinned_sets[g][arm] != set(v[arm]):
                raise Refused(f"{g} {arm}: the pinned Source Names are not exactly the rule's eligible ones")
    return list(per.values())


def _stored(acc: str, dest: Path, datasets: Path | None) -> dict:
    pp = cm.outside_repo(dest) / "meta" / acc / "plan.json"
    if not pp.is_file():
        raise Refused(f"{pp}: no plan; run plan first")
    pl = json.loads(pp.read_text(encoding="utf-8"))
    if pl.get("withheld") or acc not in (pl.get("tables") or {}) or "line_set" not in pl:
        raise Refused(f"{acc}: the plan is withheld or holds no readable table ({pl.get('reason')})")
    return pl


def pin_template(acc: str, dest: Path, arm: str, gene: str, sources: list[str], datasets: Path | None = None) -> str:
    """Ready-to-paste pins of the NAMED Source Names of one line and one arm. Refused unless the names are exactly the
    committed rule's eligible Source Names of that line and arm (no hand-picking), each pinnable, of that genotype,
    with every run a header of the stored header record (run the header step first)."""
    cm.dataset(acc, datasets)
    if arm not in cm.PIN_ARMS:
        raise Refused(f"arm must be exactly one of {cm.PIN_ARMS}")
    if not sources or len(set(sources)) != len(sources):
        raise Refused("name every Source Name to pin, each once (--source or --source-file)")
    pl = _stored(acc, dest, datasets)
    rows = {r["source_name"]: r for r in pl["tables"][acc]}
    ls = line_set(pl["tables"][acc])
    if gene not in ls:
        raise Refused(f"{gene}: not a line of the committed rule's line set {sorted(ls)}")
    if set(sources) != set(ls[gene][arm]):
        raise Refused(f"{gene} {arm}: name exactly the rule's eligible Source Names ({len(ls[gene][arm])}), not "
                      f"{len(sources)}")
    hdrs = [h for h in bulk.stored_headers(dest, acc).values() if h.get("status") == "ok"]
    if len(hdrs) != 1:
        raise Refused(f"{acc}: run the header step first ({len(hdrs)} passing header records)")
    complete = set(hdrs[0]["ties"]["complete"])
    out = [f"# pin-template {acc} (ArrayExpress): line {gene}, arm {arm}; {len(sources)} Source Names, exactly the "
           "committed rule's eligible ones", "# READ each record in the plan table before pasting (README checks B1 "
           "to B10)."]
    for sn in sources:
        r = rows[sn]
        if r["not_pinnable"] or sn not in complete:
            raise Refused(f"{sn}: not pinnable ({r['not_pinnable'] or 'a run missing from the header'})")
        h = cm.metadata_sha256(sn, r["fields"])
        if h != r["metadata_sha256"]:
            raise Refused(f"{sn}: the stored record no longer hashes to the table's value: plan again")
        out += [f"      - accession: {r['accession']}", f"        source_name: {json.dumps(sn)}",
                f"        arm: {arm}", f"        gene: {json.dumps(gene)}", f"        stage: E{r['stage']:g}",
                f"        runs: [{', '.join(r['runs'])}]", f'        metadata_sha256: "{h}"']
    return "\n".join(out) + "\n"


def series_pin_template(acc: str, dest: Path, datasets: Path | None = None) -> str:
    cm.dataset(acc, datasets)
    pl = _stored(acc, dest, datasets)
    b = pl["series"][acc]
    h = cm.series_metadata_sha256(acc, b["fields"])
    if h != b["series_metadata_sha256"]:
        raise Refused(f"{acc}: the stored series fields no longer hash to the table's value: plan again")
    return "\n".join([f"# pin-template {acc} --series: the IDF, the SDRF header row and the methods text, as read",
                      f"      - accession: {acc}", f'        series_metadata_sha256: "{h}"']) + "\n"
