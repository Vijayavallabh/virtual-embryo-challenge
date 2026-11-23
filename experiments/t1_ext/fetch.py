#!/usr/bin/env python
"""Plan, pin, then fetch, one Task 1 outside atlas from GEO into a folder OUTSIDE the repository (--dest).

Two lines of defence keep window samples out (README, "Pinned samples"):
  FIRST, THE PINS. A sample is used only when atlases.yaml pins it: its accession, the end it serves, its stage as two
  readers read it, the sha256 of its full verbatim record (t1x_common.record_text) and who read it and when; each
  series of the atlas is pinned by the sha256 of its series record. The readers read the full record line by line
  before they pin. A changed record no longer hashes to its pin and is refused.
  SECOND, THE STAGE RULE (t1x_common.record_stage). It reads the whole record and can only refuse: a pinned sample it
  reads as WINDOW or UNKNOWN, or at a stage other than its pin's, refuses the atlas.

  plan ACC --dest D     downloads only the family SOFT metadata of each series the atlas lists (D/meta/<ACC>/), and
                        writes D/meta/<ACC>/plan.json, plan_table.tsv, records.txt (with records/<ACC>.txt, one file
                        per record) and processing_text.txt. Every sample of every listed series gets one mark: USED
                        (pinned, its record hashing to the pin, read by the stage rule as permitted at the pin's stage,
                        the tissue and assay of its end, one file per role in its own GEO folder), WINDOW (a reading,
                        widened by the plug-day convention, meets (E9.5, E13.5]: never fetched; only the fields that
                        gave a reading are kept), NOT_USED (permitted but not pinned, or not needed) or UNKNOWN (no
                        reading, an unparsed stage-like token, conflicting readings or conventions). The full record of
                        every sample the stage rule lets through is printed (records.txt), with its sha256, so that
                        the readers can read it and pin it. The plan refuses when: a pin does not match its record, its
                        end or its stage; a pinned sample is not usable; the pinned samples do not cover both ends; a
                        series is not pinned or its record changed; an UNKNOWN sample the atlas would use exists; the
                        processing text is not cleared (processing_cleared); an uncleared protocol term; the protocol
                        gate (chemistry, enrichment, library strategy, molecule the same at both ends). A held atlas is
                        planned and never passes.
  pin-template ACC --dest D --gsm G [--gsm G ...]
                        prints, from the stored plan, the ready-to-paste `pinned:` block for the NAMED samples (no bulk
                        selector) and every series of the atlas, and the processing_cleared line for those samples (it
                        writes the processing text they give to meta/<ACC>/processing_text_pinned.txt). No hash is
                        typed by hand. read_by and date are placeholders the atlas check refuses until the readers
                        replace them. It never writes atlases.yaml.
  fetch-meta ACC --dest ~/v/scratch/t72_screen --plan-sha256 S
                        THE SCREENS' FETCH (row T72, cell 3), before the row's commit and before any matrix: S is the
                        sha256 of the stored plan.json (a plan nobody read cannot be fetched); the named samples
                        (cm.SCREEN_SAMPLES, the six of the row) must be candidates of that plan at their ends; the
                        series' metadata is downloaded and planned again and their records and files must be unchanged.
                        Only then are their barcodes and features files, and nothing else, downloaded from each sample's
                        own folder to ~/v/scratch/t72_screen/screen/<ACC>/<GSM>/ with FETCH_META.json; every
                        never-fetched accession (cm.NEVER_FETCHED), every matrix and every RAW.tar is refused by name.
                        It prints each file's name, bytes and sha256, then runs screens.py (the facts of checks (a) to
                        (d) and the yaml screen block).
  fetch ACC --dest D --plan-sha256 S
                        Every check of the stored and the live plan (pins included); for an atlas with screens only
                        while atlases.yaml records its screen kept, and every fetched barcodes and features file must
                        have the sha256 of the screened file the yaml records. Every data file comes from the sample's
                        own GEO folder only (never a series RAW.tar, a bigwig, a never-fetch sample's or file's name);
                        D/raw/<ACC>/MANIFEST.json holds url, bytes and sha256 per file, each sample's full record, its
                        record_sha256 and pin, the series records, the series-wide fields and the screen.

Requests go only to ftp.ncbi.nlm.nih.gov over https with a generic user agent (geo.get). Exit status 2 on a
refusal.

  uv run --frozen --no-sync python -I -B experiments/t1_ext/fetch.py plan GSE230531 --dest ~/v/scratch/t72_screen
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections.abc import Callable, Mapping
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import geo  # noqa: E402
import t1x_common as cm  # noqa: E402

Refused = cm.Refused
TISSUE_KEY = re.compile(r"tissue|organ|cell type|source|strain|genotype", re.IGNORECASE)
SERIES_PROTOCOL = cm.SERIES_PROTOCOL
SAMPLE_PROTOCOL = cm.SAMPLE_TEXT
_CHEM = re.compile(r"(?<![a-z0-9])(?:v|version\s?)([123])(?![0-9])")
MAX_FILE_BYTES = 3 * 1024 ** 3
WITHHELD = "<not kept: window sample>"
file_name = cm.file_name
supplementary = cm.supplementary
stage_fields = cm.stage_fields


# ------------------------------------------------------------------------------------------------ record fields


def tissue_fields(rec: Mapping) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for ch in rec.get("characteristics_ch1", []):
        k, _, v = ch.partition(":")
        if TISSUE_KEY.search(k):
            out.setdefault(f"characteristics: {k.strip()}", []).append(v.strip())
    for k in ("source_name_ch1", "title", "library_strategy", "molecule_ch1", "instrument_model", "platform_id"):
        if rec.get(k):
            out[k] = list(rec[k])
    return out


def protocol_text(rec: Mapping, keys) -> list[str]:
    return [v for k in keys for v in rec.get(k, [])]


def role_files(gsm: str, rec: Mapping, roles: Mapping[str, str]) -> tuple[dict | None, str]:
    """The sample's one file per role, from its own GEO folder (geo.own_file), never a denied type or an archive;
    (None, reason) when a role has no file or more than one."""
    out: dict[str, dict] = {}
    for role, pat in roles.items():
        hits = []
        for url in supplementary(rec):
            name = file_name(url)
            if not re.search(pat, name, re.IGNORECASE):
                continue
            u, why = geo.own_file(gsm, url)
            if u is None or geo.DENY.search(name):
                return None, f"{role} file {name!r} not fetchable ({why or 'a denied or archive type'})"
            hits.append({"name": name, "url": u})
        if len(hits) != 1:
            return None, f"{len(hits)} files for role {role}"
        out[role] = hits[0]
    return out, ""


def chemistry(rec: Mapping) -> list[str]:
    txt = " ".join(" ".join(rec.get(k, [])) for k in ("extract_protocol_ch1", "description", "data_processing",
                                                       "library_strategy", "characteristics_ch1", "title")).lower()
    return sorted({f"v{v}" for v in _CHEM.findall(txt)})


def enrichment(rec: Mapping, rx: re.Pattern) -> list[str]:
    txt = " ".join(" ".join(rec.get(k, [])) for k in ("extract_protocol_ch1", "treatment_protocol_ch1",
                                                       "source_name_ch1", "characteristics_ch1", "title",
                                                       "description")).lower()
    return sorted({m.group(0) for m in rx.finditer(txt)})


def protocol_hits(texts, rx: re.Pattern) -> list[dict]:
    out = []
    for t in texts:
        for sent in re.split(r"(?<=[.;])\s+", t):
            for m in rx.finditer(sent.lower()):
                out.append({"term": m.group(0), "text": sent.strip()})
    return out


# ------------------------------------------------------------------------------------------------ the plan


def would_use(a: Mapping, rec: Mapping) -> tuple[bool, str]:
    """Whether the atlas would use this sample whatever its stage: its library strategy is the atlas's assay and its
    title, source name or characteristics match an end's tissue. An UNKNOWN sample that would be used refuses the
    atlas; one that would not (another tissue or assay) is listed with the reason."""
    strat = " ".join(rec.get("library_strategy", []))
    if not re.search(a.get("assay_rule", {}).get("library_strategy", "."), strat, re.IGNORECASE):
        return False, f"library strategy {strat!r} is not the atlas's assay"
    tissue_txt = " ".join([*rec.get("title", []), *rec.get("source_name_ch1", []), *rec.get("characteristics_ch1", [])])
    ends = [e for e in ("early", "late") if re.search(a[e]["tissue"], tissue_txt, re.IGNORECASE)]
    if not ends:
        return False, "neither end's tissue"
    return True, f"the assay and the tissue of the {' and '.join(ends)} end"


def classify(a: Mapping, cfg: Mapping, gsm: str, rec: Mapping, series_protocol: list[str],
             series_text: Mapping[str, list[str]] | None = None) -> dict:
    """One sample's row before the pins are applied. Its stage is read from the whole of its own record
    (cm.record_stage), the plug-day convention from its own and its series' protocol text. A sample with any reading
    that touches the window keeps only the fields and file names that gave a reading or a leftover. A sample the rule
    lets through keeps its full record and its record_sha256. A sample at a bracket stage with the end's tissue, the
    assay and one file per role is a candidate (mark USED here; plan_from_meta keeps it USED only if it is pinned)."""
    names = [file_name(u) for u in supplementary(rec)]
    r = cm.record_stage(rec, series_protocol, series_text)
    use, why_use = would_use(a, rec)
    base = {"gsm": gsm, "mark": r["mark"], "low": r["low"], "high": r["high"], "nominal_low": r["nominal_low"],
            "nominal_high": r["nominal_high"], "reason": r["reason"], "readings": r["readings"],
            "leftovers": r["leftovers"], "would_use": use, "would_use_reason": why_use,
            "series_text_fields": r["series_text_fields"],
            "convention": {k: r["convention"][k] for k in ("name", "shift", "reason")}}
    if r["window_reading"] or r["mark"] == "WINDOW":
        read = {k: v for k, v in rec.items() if k in r["fields_read"]}
        return {**base, "stage_fields": read, "files": [n for n in names if n in r["files_read"]],
                "tissue_fields": WITHHELD}
    row = {**base, "stage_fields": cm.stage_fields(rec), "protocol_fields": cm.protocol_fields(rec), "files": names,
           "convention_text": cm.convention_text(rec, series_protocol), "convention": r["convention"],
           "tissue_fields": tissue_fields(rec)}
    if r["mark"] == "UNKNOWN":
        return row
    row.update({"mark": "NOT_USED", "record": cm.stored_record(rec),
                "record_sha256": cm.record_sha256("SAMPLE", gsm, rec)})
    point = r["nominal_low"] if r["nominal_low"] == r["nominal_high"] else None
    end = label = None
    for e in ("early", "late"):
        for lab in a[e]["labels"]:
            if point is not None and point == cm.label_day(lab):
                end, label = e, lab
    if end is None:
        row["reason"] = "not a bracket stage" if point is not None else "not one stage (a span of readings)"
        return row
    strat = " ".join(rec.get("library_strategy", []))
    if not re.search(a.get("assay_rule", {}).get("library_strategy", "."), strat, re.IGNORECASE):
        row["reason"] = f"library strategy {strat!r} is not the atlas's assay"
        return row
    tissue_txt = " ".join([*rec.get("title", []), *rec.get("source_name_ch1", []), *rec.get("characteristics_ch1", [])])
    if not re.search(a[end]["tissue"], tissue_txt, re.IGNORECASE):
        row["reason"] = f"not the {end} end's tissue ({a[end]['tissue']!r})"
        return row
    nf = cm.never_fetch_refusal(a, gsm=gsm)
    if nf:
        row["reason"] = nf
        return row
    files, why = role_files(gsm, rec, cfg["forms"][a["form"]]["roles"])
    if files is None:
        row["reason"] = why
        return row
    rx_enr = re.compile(cfg["enrichment_terms"], re.IGNORECASE)
    rx_pro = re.compile(cfg["protocol_terms"], re.IGNORECASE)
    row.update({"mark": "USED", "candidate": {"end": end, "label": label}, "end": end, "label": label,
                "role_files": files, "reason": "",
                "chemistry": chemistry(rec), "enrichment": enrichment(rec, rx_enr),
                "library_strategy": rec.get("library_strategy", []), "molecule": rec.get("molecule_ch1", []),
                "instrument": rec.get("instrument_model", []),
                "protocol_hits": protocol_hits(protocol_text(rec, SAMPLE_PROTOCOL), rx_pro),
                "processing": cm.sample_text(rec)})
    return row


def apply_pins(acc: str, rows: list[dict], pins: Mapping, a: Mapping | None = None) -> list[str]:
    """The first line: a candidate stays USED only when it is pinned and its record, end and stage match its pin;
    every other candidate is NOT_USED ('not pinned'). Returns the refusals (a pin that does not match, a pinned sample
    that is not a usable candidate, a pin whose sample is not in the atlas's series)."""
    refusals = []
    by = {r["gsm"]: r for r in rows}
    for g, p in pins["samples"].items():
        r = by.get(g)
        if r is None:
            refusals.append(f"pinned sample {g} is not in the atlas's series")
            continue
        r["pinned"] = {k: p[k] for k in ("end", "stage", "record_sha256", "read_by", "date")}
        if r["mark"] != "USED":
            refusals.append(f"pinned sample {g} is not usable: {r['mark']} ({r['reason']})")
            continue
        bad = []
        if r["record_sha256"] != p["record_sha256"]:
            bad.append(f"its record hashes to {r['record_sha256']}, not to its pin {p['record_sha256']} (the record "
                       f"changed since the readers read it: read and pin again)")
        if r["end"] != p["end"]:
            bad.append(f"pinned for the {p['end']} end, the rule places it at the {r['end']} end")
        if r["label"] != p["stage"]:
            bad.append(f"pinned at {p['stage']}, but its record states {r['label']}")
        if bad:
            refusals.append(f"pinned sample {g}: " + "; ".join(bad))
            r.update({"mark": "NOT_USED", "reason": "its pin does not match: " + "; ".join(bad)})
    for r in rows:
        if r["mark"] == "USED" and r["gsm"] not in pins["samples"]:
            r.update({"mark": "NOT_USED", "reason": "not pinned: a candidate at the "
                      f"{r['end']} end ({r['label']}); the readers read its full record (records.txt) before a pin"})
    if not cm.pins_cover_both_ends(pins):
        refusals.append(f"the pinned samples do not cover both ends (pinned ends: "
                        f"{sorted({p['end'] for p in pins['samples'].values()})}): nothing is fetched or loaded")
    refusals += cm.pins_complete(a or {}, pins)
    return refusals


def series_pin_refusals(a: Mapping, series_sha: Mapping[str, str], pins: Mapping) -> list[str]:
    out = []
    for s in a["series"]:
        p = pins["series"].get(s)
        if p is None:
            out.append(f"series {s} is not pinned (record_sha256 {series_sha[s]})")
        elif p["record_sha256"] != series_sha[s]:
            out.append(f"series {s}: its record hashes to {series_sha[s]}, not to its pin {p['record_sha256']} "
                       f"(read and pin again)")
    return out


def plan_from_meta(acc: str, metas: Mapping[str, Mapping], cfg: Mapping | None = None) -> dict:
    """The plan of one atlas from its series' parsed SOFT (geo.parse_soft output, per series)."""
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    if a["status"] == "in_use":
        raise Refused(f"{acc} is the atlas in use: it is read only through its recomputed file, never planned")
    if a["status"] == "dropped":
        raise Refused(f"{acc} was dropped by the readers' screens: it is never planned or fetched again")
    if set(metas) != set(a["series"]):
        raise Refused(f"{acc}: metadata given for {sorted(metas)}, the atlas lists {a['series']}")
    pins = cm.check_pins(acc, a)
    rows, refusals = [], []
    seen: dict[str, str] = {}
    series_records = {s: metas[s].get("series", {}) for s in a["series"]}
    series_sha = {s: cm.record_sha256("SERIES", s, series_records[s]) for s in a["series"]}
    series_protocol = {s: protocol_text(series_records[s], SERIES_PROTOCOL) for s in a["series"]}
    # the stage-parser exception (ruling row R27): per series, the protocol and processing fields byte-identical
    # across EVERY sample record of that series (all of them, in-window and other-assay samples included)
    series_text = {s: cm.series_wide_fields(metas[s]["samples"]) for s in a["series"]}
    # row T72, cell 3: the sha256 of each protocol and processing field of EVERY sample record of the series (window
    # samples included, and nothing else of a window sample), for the readers' byte-identity comparison
    field_sha = {s: cm.protocol_field_sha256(metas[s]["samples"]) for s in a["series"]}
    for series in a["series"]:
        for gsm, rec in metas[series]["samples"].items():
            if gsm in seen:
                continue  # a sample listed by two of the atlas's series is classified once
            seen[gsm] = series
            rows.append({"series": series, **classify(a, cfg, gsm, rec, series_protocol[series],
                                                      series_text[series])})
    unknown = [r["gsm"] for r in rows if r["mark"] == "UNKNOWN" and r["would_use"]]
    if unknown:
        refusals.append(f"UNKNOWN stage in samples the atlas would use {unknown}: the atlas is refused")
    refusals += apply_pins(acc, rows, pins, a)
    refusals += series_pin_refusals(a, series_sha, pins)
    chosen = {}
    for e in ("early", "late"):
        for lab in a[e]["labels"]:
            if any(r["mark"] == "USED" and r["label"] == lab for r in rows):
                chosen[e] = lab
                break
        if e not in chosen:
            refusals.append(f"no USED sample at the {e} end ({a[e]['labels']})")
    for r in rows:  # a fallback label not chosen is not used (pins name one stage per end, so this is a guard)
        if r["mark"] == "USED" and chosen.get(r["end"]) != r["label"]:
            r.update({"mark": "NOT_USED", "reason": f"{r['label']} is not the chosen {r['end']} end"})
    used = [r for r in rows if r["mark"] == "USED"]
    rx_pro = re.compile(cfg["protocol_terms"], re.IGNORECASE)
    series_hits = []
    for series in a["series"]:
        series_hits += protocol_hits(protocol_text(series_records[series], SERIES_PROTOCOL), rx_pro)
    cleared = a.get("protocol_terms_cleared") or {}
    hits = [h for h in series_hits + [h for r in used for h in r["protocol_hits"]] if h["term"] not in cleared]
    if hits:
        refusals.append(f"protocol terms not cleared: {sorted({h['term'] for h in hits})}")
    gate = {}
    gate_cleared = a.get("protocol_gate_cleared") or {}
    for key in ("chemistry", "enrichment", "library_strategy", "molecule"):
        sides = {e: sorted({v for r in used if r["end"] == e for v in r[key]}) for e in ("early", "late")}
        gate[key] = {**sides, "same": sides["early"] == sides["late"], "cleared": key in gate_cleared}
        if not gate[key]["same"] and key not in gate_cleared:
            refusals.append(f"protocol gate: {key} differs between the ends ({sides})")
    if chosen.get("early") and chosen.get("late"):
        try:
            cm.bracket_ok(cm.label_day(chosen["early"]), cm.label_day(chosen["late"]))
        except Refused as e:
            refusals.append(str(e))
    text = processing_for(acc, a, series_records, {r["gsm"]: (r["series"], r["record"]) for r in used},
                          series_text)
    sha = cm.text_sha256(text)
    ok, why = cm.clearance(a, sha)
    if not ok:
        refusals.append(why)
    if a["status"] == "held":
        refusals.append("the atlas is held: planned only, never fetched")
    marks = {m: sum(r["mark"] == m for r in rows) for m in cm.MARKS}
    return {"accession": acc, "status": "passed" if not refusals else "refused", "refusals": refusals,
            "atlas_entry_sha256": cm.entry_sha256(acc, cfg), "series": a["series"], "chosen": chosen,
            "marks": marks, "protocol_gate": gate, "series_protocol_hits": series_hits, "samples": rows,
            "series_records": {s: cm.stored_record(series_records[s]) for s in a["series"]},
            "series_record_sha256": series_sha, "series_protocol": series_protocol, "series_text": series_text,
            "protocol_field_sha256": field_sha, "processing_text": text, "processing_sha256": sha,
            "processing_clearance": why if ok else "not cleared", "tasks": list(cm.TASKS)}


def processing_for(acc: str, a: Mapping, series_records: Mapping, samples: Mapping[str, tuple[str, object]],
                   series_text: Mapping[str, Mapping[str, list[str]]]) -> str:
    """The canonical processing text of these samples ({GSM: (GSE, record)}), with each series' series-wide fields."""
    return cm.processing_text(acc, [(s, series_records[s]) for s in a["series"]],
                              {g: {"series": s, "text": cm.sample_text(rec)} for g, (s, rec) in samples.items()},
                              series_text)


TABLE_COLUMNS = ("gsm", "series", "mark", "would_use", "end", "label", "pinned", "record_sha256", "nominal_low",
                 "nominal_high", "low", "high", "reason", "would_use_reason", "convention", "readings", "leftovers",
                 "stage_fields", "protocol_fields", "files", "tissue_fields", "chemistry", "enrichment",
                 "protocol_hits")


def records_text(pl: Mapping) -> str:
    """Every series record and the full record of every sample the stage rule lets through (mark USED or NOT_USED),
    each after a header line naming its sha256; the bytes between a header and the next are exactly the canonical
    record text (records/<ACC>.txt holds the same bytes, so `sha256sum records/<ACC>.txt` gives the hash)."""
    out = [f"# RECORDS of {pl['accession']}: canonical form {cm.RECORD_FORMAT} (README, 'Canonical record text').",
           "# Every series record, then the full record of every sample the stage rule lets through. Read each line,",
           "# the protocol and processing text included, and every 'escaped:' line, before a pin."]
    for s in pl["series"]:
        out.append(f"==== SERIES {s}  record_sha256 {pl['series_record_sha256'][s]}")
        out.append(cm.record_text("SERIES", s, pl["series_records"][s]).rstrip("\n"))
        st = sorted((pl.get("series_text") or {}).get(s) or {})
        out.append(f"==== SERIES-WIDE FIELDS of {s} (byte-identical in every sample record of the series, read as "
                   f"series text under ruling row R27: no stage reading, pooling words refuse, the stage comes from "
                   f"each sample's own title and characteristics): {', '.join(st) if st else '(none)'}")
    for r in pl["samples"]:
        if "record" not in r:
            continue
        cand = r.get("candidate")
        out.append(f"==== SAMPLE {r['gsm']} ({r['series']})  mark {r['mark']}  "
                   f"candidate {cand['end'] + ' ' + cand['label'] if cand else 'no'}  "
                   f"pinned {'yes' if r.get('pinned') else 'no'}  record_sha256 {r['record_sha256']}")
        out.append(cm.record_text("SAMPLE", r["gsm"], r["record"]).rstrip("\n"))
    return "\n".join(out) + "\n"


def write_table(pl: Mapping, mdir: Path) -> Path:
    """plan_table.tsv: one row per sample (its stored fields verbatim, JSON-escaped, so a hidden character shows), the
    plan's verdict, the records (records_text), then the canonical processing text. records.txt holds the records,
    records/<ACC>.txt each record's exact bytes, processing_text.txt exactly the processing text."""
    lines = ["\t".join(TABLE_COLUMNS)]
    for r in pl["samples"]:
        lines.append("\t".join(json.dumps(r.get(c), default=str) if not isinstance(r.get(c), str) else
                               json.dumps(r[c])[1:-1] for c in TABLE_COLUMNS))
    lines.append("")
    lines.append(f"# status {pl['status']}; refusals {json.dumps(pl['refusals'])}; chosen {json.dumps(pl['chosen'])}; "
                 f"marks {json.dumps(pl['marks'])}")
    recs = records_text(pl)
    lines.append(recs.rstrip("\n"))
    lines.append(f"# PROCESSING TEXT, canonical form {cm.PROCESSING_FORMAT}, sha256 {pl['processing_sha256']} "
                 f"(the same bytes are processing_text.txt; clearance: {pl['processing_clearance']})")
    path = cm.outside_repo(Path(mdir) / "plan_table.tsv")
    path.write_text("\n".join(lines) + "\n" + pl["processing_text"], encoding="utf-8")
    cm.outside_repo(Path(mdir) / "records.txt").write_text(recs, encoding="utf-8")
    rdir = cm.outside_repo(Path(mdir) / "records")
    rdir.mkdir(parents=True, exist_ok=True)
    for old in rdir.glob("*.txt"):
        old.unlink()
    for s in pl["series"]:
        (rdir / f"{s}.txt").write_bytes(cm.record_text("SERIES", s, pl["series_records"][s]).encode("utf-8"))
    for r in pl["samples"]:
        if "record" in r:
            (rdir / f"{r['gsm']}.txt").write_bytes(cm.record_text("SAMPLE", r["gsm"], r["record"]).encode("utf-8"))
    tp = cm.outside_repo(Path(mdir) / "processing_text.txt")
    tp.write_bytes(pl["processing_text"].encode("utf-8"))
    if cm.sha256_file(tp) != pl["processing_sha256"]:
        raise Refused("processing_text.txt does not hash to the plan's processing_sha256")
    cm.outside_repo(Path(mdir) / "protocol_field_sha256.tsv").write_text(
        "".join(ln + "\n" for ln in field_sha_lines(pl)), encoding="utf-8")
    return path


def field_sha_lines(pl: Mapping) -> list[str]:
    """One line per sample record of every series and per protocol or processing field: 'field_sha256 <GSE> <GSM>
    <field> <sha256 or (absent)>' (row T72, cell 3; window samples included, nothing else of them)."""
    out = []
    for s in pl["series"]:
        for g, fields in sorted((pl.get("protocol_field_sha256") or {}).get(s, {}).items()):
            out += [f"field_sha256\t{s}\t{g}\t{k}\t{h}" for k, h in sorted(fields.items())]
    return out


def default_getter(url: str, dst: Path) -> None:
    geo.get(url, dst, MAX_FILE_BYTES)


def read_series(series: str, mdir: Path, getter: Callable[[str, Path], None] | None = None) -> tuple[dict, str]:
    dst = cm.outside_repo(mdir / f"{series}_family.soft.gz")
    (getter or default_getter)(geo.soft_url(series), dst)
    raw = gzip.decompress(dst.read_bytes())
    return geo.parse_soft(geo.decode(raw, series)), cm.sha256_file(dst)


def plan(acc: str, dest: Path, getter=None, cfg: Mapping | None = None) -> dict:
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    if a["status"] == "in_use":
        raise Refused(f"{acc} is the atlas in use: nothing is planned or fetched for it")
    mdir = cm.outside_repo(Path(dest) / "meta" / acc)
    mdir.mkdir(parents=True, exist_ok=True)
    cm.dump_json({"accession": acc, "status": "in_progress"}, mdir / "plan.json")
    metas, shas = {}, {}
    for s in a["series"]:
        metas[s], shas[s] = read_series(s, mdir, getter)
    pl = plan_from_meta(acc, metas, cfg)
    pl["soft_sha256"] = shas
    write_table(pl, mdir)
    cm.dump_json(pl, mdir / "plan.json")
    return pl


# ------------------------------------------------------------------------------------------------ pin template


def _comment(key: str, values) -> list[str]:
    return [f"        # {key}[{i}]: {v if cm._plain(v) else 'escaped ' + v.encode('unicode_escape').decode('ascii')}"
            for i, v in enumerate(values, 1)]


def pin_template(acc: str, dest: Path, gsms: list[str], cfg: Mapping | None = None) -> str:
    """The ready-to-paste pinned block for the NAMED samples and every series of the atlas, and the processing_cleared
    line for those samples, from the stored plan. Refused for an unnamed or unknown sample, a sample the stage rule
    does not let through or that is not a candidate at a bracket end, a stored record that no longer hashes to the
    plan's value, and a set of samples that does not cover both ends."""
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    if a["status"] not in cm.FETCHABLE:
        raise Refused(f"{acc} has status {a['status']!r}: it is never pinned")
    if not gsms or not all(isinstance(g, str) and re.fullmatch(r"GSM\d+", g) for g in gsms):
        raise Refused("name every sample to pin with --gsm GSMnnn (one per sample; there is no bulk selector)")
    if len(set(gsms)) != len(gsms):
        raise Refused("a sample is named twice")
    mdir = cm.outside_repo(Path(dest) / "meta" / acc)
    pp = mdir / "plan.json"
    if not pp.is_file():
        raise Refused(f"{pp}: no plan; run plan first")
    pl = json.loads(pp.read_text(encoding="utf-8"))
    if "samples" not in pl or "series_records" not in pl:
        raise Refused(f"{acc}: the stored plan holds no records: plan again")
    rows = {r["gsm"]: r for r in pl["samples"]}
    series_records = {s: pl["series_records"][s] for s in a["series"]}
    series_text = pl.get("series_text") or {}
    named = cm.named_gsms(a)
    if named and set(gsms) != set(named):
        raise Refused(f"the strict form pins exactly the samples the atlas names ({sorted(named)}), not "
                      f"{sorted(gsms)}")
    picked = {}
    for g in gsms:
        r = rows.get(g)
        if r is None:
            raise Refused(f"{g} is not in the plan of {acc}")
        if "record" not in r or not r.get("candidate"):
            raise Refused(f"{g}: not a candidate ({r['mark']}: {r['reason']}): it cannot be pinned")
        h = cm.record_sha256("SAMPLE", g, r["record"])
        if h != r["record_sha256"]:
            raise Refused(f"{g}: the stored record hashes to {h}, not to the plan's value: plan again")
        st = cm.record_stage(r["record"], pl["series_protocol"][r["series"]], series_text.get(r["series"]))
        if st["mark"] != "PERMITTED" or st["nominal_low"] != st["nominal_high"] or \
                st["nominal_low"] != cm.label_day(r["candidate"]["label"]):
            raise Refused(f"{g}: the stage rule no longer reads it as {r['candidate']['label']} ({st['reason']})")
        if named and named[g] != r["candidate"]["end"]:
            raise Refused(f"{g}: the atlas names it at the {named[g]} end, the rule places it at the "
                          f"{r['candidate']['end']} end")
        picked[g] = r
    ends = {r["candidate"]["end"] for r in picked.values()}
    if ends != {"early", "late"}:
        raise Refused(f"the named samples cover the ends {sorted(ends)}, not both: nothing would be fetched")
    for e in ("early", "late"):
        st = {r["candidate"]["label"] for r in picked.values() if r["candidate"]["end"] == e}
        if len(st) > 1:
            raise Refused(f"the named samples of the {e} end are at two stages {sorted(st)}: name one stage")
    text = processing_for(acc, a, series_records, {g: (r["series"], r["record"]) for g, r in picked.items()},
                          series_text)
    tp = cm.outside_repo(mdir / "processing_text_pinned.txt")
    tp.write_bytes(text.encode("utf-8"))
    out = [f"# pin-template {acc}: {len(gsms)} sample(s) named; from {pp} (sha256 {cm.sha256_file(pp)})",
           "# READ each named sample's full record (records/<GSM>.txt) and every series record (records/<GSE>.txt)",
           f"# line by line, and {tp.name}, before pasting. Two readers each confirm the checklist (README).",
           "# Replace READER_1, READER_2 and the date; a reader who doubts does not pin. Paste under the atlas.",
           "    pinned:", "      samples:"]
    for g in gsms:
        r = picked[g]
        d = cm.record_dict(r["record"])
        out += [f"        # ---- {g} ({r['series']}): title and characteristics as read in its record",
                *_comment("title", d.get("title", [])),
                *_comment("characteristics_ch1", d.get("characteristics_ch1", [])),
                f"        - accession: {g}",
                f"          end: {r['candidate']['end']}",
                f"          stage: {r['candidate']['label']}      # the rule's reading; the readers confirm it",
                f'          record_sha256: "{r["record_sha256"]}"',
                "          read_by: [READER_1, READER_2]",
                '          date: "YYYY-MM-DD"']
    out.append("      series:")
    for s in a["series"]:
        h = cm.record_sha256("SERIES", s, series_records[s])
        if h != pl["series_record_sha256"][s]:
            raise Refused(f"{s}: the stored series record hashes to {h}, not to the plan's value: plan again")
        out += [*_comment("title", cm.record_dict(series_records[s]).get("title", [])),
                f"        - accession: {s}",
                f'          record_sha256: "{h}"',
                "          read_by: [READER_1, READER_2]",
                '          date: "YYYY-MM-DD"']
    out.append(f'    processing_cleared: {{sha256: "{cm.text_sha256(text)}", read_by: [READER_1, READER_2], '
               f'date: "YYYY-MM-DD"}}')
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------------------------------------ the fetch


def stored_plan(acc: str, dest: Path) -> tuple[Path, dict]:
    p = cm.outside_repo(Path(dest) / "meta" / acc / "plan.json")
    if not p.exists():
        raise Refused(f"{acc}: no plan at {p}")
    return p, json.loads(p.read_text())


def check_plan(acc: str, dest: Path, plan_sha256: str, cfg: Mapping | None = None) -> dict:
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    if a["status"] not in cm.FETCHABLE:
        raise Refused(f"{acc} has status {a['status']!r}: it is never fetched")
    pins = cm.check_pins(acc, a)
    if not cm.pins_cover_both_ends(pins):
        raise Refused(f"{acc}: the pinned samples do not cover both ends: nothing is fetched")
    p, pl = stored_plan(acc, dest)
    if cm.sha256_file(p) != plan_sha256:
        raise Refused(f"{acc}: the plan's sha256 is not the one given: read the plan table, then pass its sha256")
    if pl.get("status") != "passed":
        raise Refused(f"{acc}: the plan did not pass ({pl.get('refusals')})")
    if pl.get("atlas_entry_sha256") != cm.entry_sha256(acc, cfg):
        raise Refused(f"{acc}: the plan was made from another atlas entry (pins included): plan again")
    if cm.text_sha256(pl.get("processing_text", "")) != pl.get("processing_sha256"):
        raise Refused(f"{acc}: the stored processing text does not hash to the plan's processing_sha256")
    ok, why = cm.clearance(a, pl["processing_sha256"])
    if not ok:
        raise Refused(f"{acc}: {why}")
    return pl


def _used_signature(pl: Mapping) -> dict:
    return {r["gsm"]: {"role_files": r["role_files"], "record": r["record"], "record_sha256": r["record_sha256"],
                       "label": r["label"], "end": r["end"]}
            for r in pl["samples"] if r["mark"] == "USED"}


def check_used_against_pins(acc: str, a: Mapping, pl: Mapping) -> dict:
    """Each USED sample of a plan: pinned, its stored record hashing to its pin, its end and stage its pin's, and the
    stage rule reading its stored record as permitted at that stage; the USED samples cover both ends; every series
    record hashes to its pin. Returns the pins."""
    pins = cm.check_pins(acc, a)
    used = [r for r in pl["samples"] if r["mark"] == "USED"]
    for r in used:
        p = pins["samples"].get(r["gsm"])
        if p is None:
            raise Refused(f"{r['gsm']}: USED but not pinned: refused")
        h = cm.record_sha256("SAMPLE", r["gsm"], r["record"])
        if not (h == r["record_sha256"] == p["record_sha256"]) or p["end"] != r["end"] or p["stage"] != r["label"]:
            raise Refused(f"{r['gsm']}: its record, end or stage is not its pin's: refused")
        st = cm.record_stage(r["record"], pl["series_protocol"][r["series"]],
                             (pl.get("series_text") or {}).get(r["series"]))
        if st["mark"] != "PERMITTED" or st["nominal_low"] != st["nominal_high"] or \
                st["nominal_low"] != cm.label_day(r["label"]):
            raise Refused(f"{r['gsm']}: its stored record no longer reads as {r['label']} ({st['reason']})")
        why = cm.never_fetch_refusal(a, gsm=r["gsm"])
        if why:
            raise Refused(f"{r['gsm']}: {why}: refused")
    if {r["end"] for r in used} != {"early", "late"}:
        raise Refused(f"{acc}: the USED samples do not cover both ends")
    bad = cm.pins_complete(a, pins)
    if bad:
        raise Refused(f"{acc}: {'; '.join(bad)}")
    for s in a["series"]:
        h = cm.record_sha256("SERIES", s, pl["series_records"][s])
        if s not in pins["series"] or pins["series"][s]["record_sha256"] != h:
            raise Refused(f"series {s}: its record does not hash to its pin: refused")
    return pins


def _live_checked(acc: str, dest: Path, plan_sha256: str, getter, cfg: Mapping) -> tuple[dict, dict, dict]:
    """Before any data file: the stored plan is checked (its sha256, its atlas entry with the pins, the processing
    clearance), then the series' metadata is downloaded again and planned again: the live plan must pass, its
    processing text must still hash to processing_cleared, its USED samples, their full records (each hashing to its
    pin), their files and the series-wide fields must be the stored plan's. Returns (plan, atlas, pins)."""
    pl = check_plan(acc, dest, plan_sha256, cfg)
    a = cm.atlas(acc, cfg)
    check_used_against_pins(acc, a, pl)
    live_dir = cm.outside_repo(Path(dest) / "meta" / acc / "live")
    live_dir.mkdir(parents=True, exist_ok=True)
    live = plan_from_meta(acc, {s: read_series(s, live_dir, getter)[0] for s in a["series"]}, cfg)
    if live["status"] != "passed":
        raise Refused(f"{acc}: the live metadata no longer plans as passed ({live['refusals']})")
    if live["processing_sha256"] != pl["processing_sha256"] or _used_signature(live) != _used_signature(pl) or \
            live["series_record_sha256"] != pl["series_record_sha256"] or \
            live["series_text"] != (pl.get("series_text") or {}):
        raise Refused(f"{acc}: the live metadata differs from the plan that was read (processing text, USED samples, "
                      f"their records or files, the series records, the series-wide fields): plan again and read it")
    pins = check_used_against_pins(acc, a, live)
    return pl, a, pins


def _get_data(a: Mapping, gsm: str, f: Mapping, dst_dir: Path, getter) -> dict:
    """One data file of one USED sample: only from the sample's own GEO folder, never a denied type, an archive or a
    file the never-fetch list names."""
    u, why = geo.own_file(gsm, f["url"])
    if u is None or geo.DENY.search(f["name"]):
        raise Refused(f"{gsm}: {f['name']} is not fetchable ({why or 'a denied or archive type'})")
    nf = cm.never_fetch_refusal(a, gsm=gsm, url=u)
    if nf:
        raise Refused(f"{gsm}: {nf}")
    dst = cm.outside_repo(dst_dir / gsm / f["name"])
    (getter or default_getter)(u, dst)
    return {"name": f["name"], "url": u, "bytes": dst.stat().st_size, "sha256": cm.sha256_file(dst)}


SCREEN_ROLES = cm.SCREEN_ROLES


def screen_dir_ok(dest: Path, screen_root: Path | None = None) -> Path:
    """The screens' folder: exactly ~/v/scratch/t72_screen (cm.SCREEN_ROOT; row T72, cell 3), outside every checkout."""
    root = Path(screen_root or cm.SCREEN_ROOT).expanduser().resolve()
    d = cm.outside_repo(dest)
    if d != root:
        raise Refused(f"the screens write only under {root}, not {d}")
    return d


def _requested(a: Mapping, gsm: str, role: str, f: Mapping, matrix_rx: str) -> str:
    """One barcodes or features file of a named sample, refused by name when it is a matrix, a series archive, a
    never-fetched accession's (row T72's literal list or the yaml's), outside the sample's own folder or denied."""
    if role not in SCREEN_ROLES:
        raise Refused(f"{gsm}: role {role!r}: the screens fetch barcodes and features only")
    name = f["name"]
    if re.search(matrix_rx, name, re.IGNORECASE) or re.search(r"matrix|\.mtx|\.h5$", name, re.IGNORECASE):
        raise Refused(f"{gsm}: {name} is a matrix file: no matrix is fetched before the row's commit")
    u, why = geo.own_file(gsm, f["url"])
    if u is None or geo.DENY.search(name):
        raise Refused(f"{gsm}: {name} is not fetchable ({why or 'a denied or archive type'})")
    nf = cm.never_fetched_refusal(gsm, name, u) or cm.never_fetch_refusal(a, gsm=gsm, url=u)
    if nf:
        raise Refused(f"{gsm}: {nf}")
    return u


def fetch_meta(acc: str, dest: Path, plan_sha256: str, getter=None, cfg: Mapping | None = None,
               screen_root: Path | None = None, samples: Mapping[str, Mapping[str, str]] | None = None) -> dict:
    """THE SCREENS' FETCH (row T72, cell 3), before the row's commit and before any matrix: only the barcodes and
    features files of the six named samples (cm.SCREEN_SAMPLES), each from the sample's own folder, into
    ~/v/scratch/t72_screen/screen/<ACC>/<GSM>/ (never the row's data folders). Before any request: the stored plan
    (meta/<ACC>/plan.json under the same folder) hashes to `plan_sha256` and was made under the atlas entry in force;
    every named sample is a candidate of the plan at its named end, its stored record hashing to the plan's value and
    read by the stage rule as permitted at that stage; the series' metadata is downloaded again and planned again and
    the named samples' records and files are unchanged. Every never-fetched accession, every matrix and every RAW.tar
    is refused by name. Writes screen/<ACC>/FETCH_META.json (url, bytes, sha256 per file) and returns it."""
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    want = dict((samples or cm.SCREEN_SAMPLES).get(acc) or {})
    if a["status"] not in cm.FETCHABLE or a.get("form") != "tenx_mtx" or not a.get("screens") or not want:
        raise Refused(f"{acc}: the screens' fetch is for row T72's two 10x members only")
    if dict(cm.named_gsms(a)) != want:
        raise Refused(f"{acc}: the atlas names {cm.named_gsms(a)}, not row T72's samples {want}")
    root = screen_dir_ok(dest, screen_root)
    pp = root / "meta" / acc / "plan.json"
    if not pp.is_file() or cm.sha256_file(pp) != plan_sha256:
        raise Refused(f"{acc}: {pp} is missing or does not hash to the plan sha256 given (read the plan, then pass "
                      f"its sha256)")
    pl = json.loads(pp.read_text())
    if pl.get("atlas_entry_sha256") != cm.entry_sha256(acc, cfg):
        raise Refused(f"{acc}: the plan was made under another atlas entry: plan again")
    rows = {r["gsm"]: r for r in pl.get("samples", [])}
    for g, end in want.items():
        r = rows.get(g)
        if r is None or "record" not in r or not r.get("candidate") or r["candidate"]["end"] != end:
            raise Refused(f"{g}: not a candidate of the plan at the {end} end: nothing is fetched")
        if cm.record_sha256("SAMPLE", g, r["record"]) != r["record_sha256"]:
            raise Refused(f"{g}: the stored record does not hash to the plan's value: plan again")
        st = cm.record_stage(r["record"], pl["series_protocol"][r["series"]],
                             (pl.get("series_text") or {}).get(r["series"]))
        if st["mark"] != "PERMITTED" or st["nominal_low"] != st["nominal_high"] or \
                st["nominal_low"] != cm.label_day(r["candidate"]["label"]):
            raise Refused(f"{g}: the stage rule no longer reads it as {r['candidate']['label']} ({st['reason']})")
    matrix_rx = cfg["forms"][a["form"]]["roles"]["matrix"]
    todo = []
    for g, end in want.items():
        for role in SCREEN_ROLES:
            f = rows[g]["role_files"][role]
            todo.append((g, end, role, f, _requested(a, g, role, f, matrix_rx)))  # every refusal before any request
    live_dir = cm.outside_repo(root / "meta" / acc / "live")
    live_dir.mkdir(parents=True, exist_ok=True)
    live = plan_from_meta(acc, {s: read_series(s, live_dir, getter)[0] for s in a["series"]}, cfg)
    lrows = {r["gsm"]: r for r in live["samples"]}
    for g in want:
        lr = lrows.get(g) or {}
        if lr.get("record_sha256") != rows[g]["record_sha256"] or lr.get("role_files") != rows[g].get("role_files") \
                or lr.get("candidate") != rows[g]["candidate"]:
            raise Refused(f"{g}: the live record or its files differ from the plan that was read: plan again")
    out_dir = cm.outside_repo(root / "screen" / acc)
    out = {"accession": acc, "plan_sha256": plan_sha256, "atlas_entry_sha256": pl["atlas_entry_sha256"],
           "tasks": list(cm.TASKS), "roles": list(SCREEN_ROLES), "samples": {}}
    for g, end, role, f, u in todo:
        dst = cm.outside_repo(out_dir / g / f["name"])
        (getter or default_getter)(u, dst)
        out["samples"].setdefault(g, {"end": end, "files": {}})["files"][role] = {
            "name": f["name"], "url": u, "bytes": dst.stat().st_size, "sha256": cm.sha256_file(dst)}
    cm.dump_json(out, out_dir / "FETCH_META.json")
    return out


def fetch(acc: str, dest: Path, plan_sha256: str, getter=None, cfg: Mapping | None = None) -> dict:
    """The USED samples' files, after every check of _live_checked. For an atlas with screens, only while atlases.yaml
    records the screen kept, and each fetched barcodes and features file must have the sha256 of the screened file
    the yaml records (row T72, cell 3: the fetched files are tied to the screened copies). Writes
    D/raw/<ACC>/MANIFEST.json."""
    cfg = cfg or cm.load_atlases()
    pl, a, pins = _live_checked(acc, dest, plan_sha256, getter, cfg)
    screen = cm.screen_kept(acc, a)
    rdir = cm.outside_repo(Path(dest) / "raw" / acc)
    rdir.mkdir(parents=True, exist_ok=True)
    man = {"accession": acc, "plan_sha256": plan_sha256, "atlas_entry_sha256": pl["atlas_entry_sha256"],
           "chosen": pl["chosen"], "form": a["form"], "tasks": list(cm.TASKS),
           "series_records": pl["series_records"], "series_record_sha256": pl["series_record_sha256"],
           "series_text": pl.get("series_text") or {},
           "processing_sha256": pl["processing_sha256"], "pins": {"samples": pins["samples"],
                                                                  "series": pins["series"]},
           "screen": screen, "samples": {}}
    for r in pl["samples"]:
        if r["mark"] != "USED":
            continue
        files = {}
        for role, f in r["role_files"].items():
            files[role] = _get_data(a, r["gsm"], f, rdir, getter)
            if screen is not None and role in SCREEN_ROLES:
                want = screen[r["gsm"]][role]
                if files[role]["name"] != want["name"] or files[role]["sha256"] != want["sha256"]:
                    raise Refused(f"{r['gsm']}: {f['name']} is not the file the readers screened (sha256 "
                                  f"{files[role]['sha256']}, the yaml records {want['sha256']})")
        man["samples"][r["gsm"]] = {"series": r["series"], "end": r["end"], "label": r["label"],
                                    "record": r["record"], "record_sha256": r["record_sha256"],
                                    "series_protocol": pl["series_protocol"][r["series"]],
                                    "series_text": (pl.get("series_text") or {}).get(r["series"], {}),
                                    "files": files}
    cm.dump_json(man, rdir / "MANIFEST.json")
    return man


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("plan")
    p1.add_argument("accession")
    p1.add_argument("--dest", required=True, type=Path)
    for name in ("fetch-meta", "fetch"):
        p2 = sub.add_parser(name)
        p2.add_argument("accession")
        p2.add_argument("--dest", required=True, type=Path)
        p2.add_argument("--plan-sha256", required=True)
    p3 = sub.add_parser("pin-template")
    p3.add_argument("accession")
    p3.add_argument("--dest", required=True, type=Path)
    p3.add_argument("--gsm", action="append", required=True)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "plan":
            pl = plan(a.accession, a.dest)
            print(json.dumps({k: pl[k] for k in ("accession", "status", "refusals", "chosen", "marks")}, indent=1))
            print("\n".join(field_sha_lines(pl)))
            return 0 if pl["status"] == "passed" else 2
        if a.cmd == "pin-template":
            print(pin_template(a.accession, a.dest, a.gsm), end="")
            return 0
        if a.cmd == "fetch-meta":
            fm = fetch_meta(a.accession, a.dest, a.plan_sha256)
            for g, s in sorted(fm["samples"].items()):
                for role, f in sorted(s["files"].items()):
                    print(f"fetched\t{a.accession}\t{g}\t{role}\t{f['name']}\tbytes {f['bytes']}\tsha256 {f['sha256']}")
            return cm.sibling("screens").main([a.accession, "--dest", str(a.dest)])
        man = fetch(a.accession, a.dest, a.plan_sha256)
        print(json.dumps({"accession": a.accession, "samples": sorted(man["samples"])}, indent=1))
        return 0
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
