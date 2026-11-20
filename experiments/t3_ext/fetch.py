#!/usr/bin/env python
"""Plan, pin and fetch one cleared outside-knockout dataset from GEO into a folder OUTSIDE the repository.

The first line of defence is the PINNED SAMPLE LIST in datasets.yaml (README, "The pinned sample list"); the text
screen (t3x_screen) is the second. Every step reads the committed datasets.yaml only (its sha256 must equal
t3x_common.DATASETS_SHA256); no step takes another dataset list, and the NEVER_FETCH accessions are refused whatever
the list says.

  plan ACC --dest D    downloads only the family SOFT metadata of every series in screen_series
                       (D/meta/<ACC>/<series>_family.soft.gz, sha256 recorded) and writes D/meta/<ACC>/plan.json plus
                       one plan table per series, D/meta/<ACC>/plan_table_<SERIES>.txt: EVERY sample of the series, its
                       accession, the verbatim text of every field (all are screened), the sha256 of their canonical
                       serialisation (t3x_common.canonical_metadata), the screen's verdict and hits, the per-sample
                       files, the parsed stage and chemistry, and the selectors' SUGGESTED arm (a guess, never used).
                       A held-out hit anywhere refuses the dataset and withholds every text of it from the plan (only
                       accessions, verdicts and the fields hit are written). A review term not cleared refuses it too.
                       If the dataset already has pins, each pinned sample's LIVE metadata must hash to its pin.
  pin-template ACC --dest D --arm ko|control --gene G (--gsm GSM ... | --title REGEX | --suggested)
                       prints a ready-to-paste `pinned` block for the samples the selector matches, from the stored
                       plan (hashes recomputed from the stored fields and checked against the table). It never writes
                       datasets.yaml: the readers paste it after reading each sample's table entry.
  check-pins ACC ...   the gate before any fetch: each accession is cleared, has a non-empty pinned block of the
                       right form, and the committed datasets.yaml hashes to its pin. Exit 2 otherwise.
  datasets-sha         prints datasets.yaml's sha256 and t3x_common.DATASETS_SHA256 (exit 2 if they differ).
  fetch ACC --dest D   downloads exactly the pinned samples' pinned files to D/raw/<ACC>/<GSM>/ (.part, then renamed)
                       and records each file's url, size and sha256, the pins and their sha256 in
                       D/raw/<ACC>/MANIFEST.json. It refuses a dataset with no pins, a plan that did not pass or was
                       made from another datasets.yaml, a pinned accession that is not in its own series' plan table (or
                       belongs to another series), a pinned sample whose stored metadata no longer hashes to its pin or
                       no longer passes the screen, a file not in that sample's plan entry. Nothing is unpacked, no
                       count is read, and no file of an unpinned sample is requested.

Requests go only to ftp.ncbi.nlm.nih.gov over https, with a generic user agent and nothing that names a team, a
person or a challenge. Exit status 2 on a refusal.

  PATH=$HOME/.local/bin:$PATH uv run --frozen --no-sync python -I -B experiments/t3_ext/fetch.py plan GSE168768 \\
      --dest ~/v/vec-ext/t3
"""

from __future__ import annotations

import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import t3x_common as cm  # noqa: E402
import t3x_screen as sc  # noqa: E402

Refused = cm.Refused
HOSTS = ("ftp.ncbi.nlm.nih.gov",)
USER_AGENT = "python-urllib (public metadata and data download)"
# never fetched, whatever a dataset's patterns say
DENY = re.compile(r"fastq|\.bam$|\.bai$|fragments|atac|peaks?[._]|\.bw$|bigwig|\.bed(\.gz)?$|raw_feature_bc_matrix|"
                  r"cloupe|\.sra$|_raw\.tar$|\.h5seurat$", re.IGNORECASE)
_STAGE_KEY = re.compile(r"stage|age|time|day", re.IGNORECASE)
_STAGE_VAL = re.compile(r"(?<![a-z0-9])(?:e\s?|embryonic day\s*)(\d+(?:\.\d+)?)", re.IGNORECASE)
_CHEM = re.compile(r"(?<![a-z0-9])v\s?([23])(?![0-9])", re.IGNORECASE)
MAX_BYTES_DEFAULT = 40 * 1024 ** 3
WITHHELD = "<withheld: held-out hit in this dataset>"


# ---------------------------------------------------------------------------------------------------------------------
# GEO metadata


def soft_url(series: str) -> str:
    if not re.fullmatch(r"GSE\d+", series):
        raise Refused(f"{series}: only GEO series are fetched")
    stub = series[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{series}/soft/{series}_family.soft.gz"


def https_url(url: str) -> str:
    """A GEO file URL as https on an allowed host, else Refused."""
    u = url.strip()
    if u.startswith("ftp://"):
        u = "https://" + u[len("ftp://"):]
    p = urllib.parse.urlparse(u)
    if p.scheme != "https" or p.hostname not in HOSTS:
        raise Refused(f"{url}: not an allowed public repository URL")
    return u


def parse_soft(text: str) -> dict:
    """{'series': {key: [values]}, 'samples': {GSM: {key: [values]}}} from a family SOFT text. Values are stripped of
    surrounding whitespace and otherwise kept verbatim, in file order."""
    series: dict[str, list[str]] = {}
    samples: dict[str, dict[str, list[str]]] = {}
    cur, cur_kind = None, None
    for line in text.splitlines():
        if line.startswith("^"):
            kind, _, ident = line[1:].partition("=")
            cur_kind, cur = kind.strip().upper(), ident.strip()
            if cur_kind == "SAMPLE":
                samples.setdefault(cur, {})
            continue
        if not line.startswith("!") or "=" not in line:
            continue
        key, _, val = line[1:].partition("=")
        key, val = key.strip(), val.strip()
        if cur_kind == "SAMPLE" and key.startswith("Sample_"):
            samples[cur].setdefault(key[len("Sample_"):], []).append(val)
        elif cur_kind == "SERIES" and key.startswith("Series_"):
            series.setdefault(key[len("Series_"):], []).append(val)
    return {"series": series, "samples": samples}


def sample_text(rec: dict) -> str:
    keys = ("title", "source_name_ch1", "characteristics_ch1", "description")
    return " | ".join(f"{k}: " + "; ".join(rec.get(k, [])) for k in keys).lower()


def sample_stage(rec: dict) -> tuple[float | None, str]:
    """(E-day, how) from the characteristics (keys naming stage, age, time or day), else the title. Two different
    values or a range give None (fail closed); so does a value after E9.5."""
    found = []
    for ch in rec.get("characteristics_ch1", []):
        k, _, v = ch.partition(":")
        if _STAGE_KEY.search(k):
            found.append(v.strip())
    if not found:
        found = [t for t in rec.get("title", [])]
    vals = set()
    for f in found:
        if re.search(r"\d\s*(?:-|–|to)\s*e?\s*\d", f, re.IGNORECASE):
            return None, f"range in {f!r}"
        vals |= {float(x) for x in _STAGE_VAL.findall(f)}
    if len(vals) != 1:
        return None, f"{len(vals)} stage values"
    e = vals.pop()
    try:
        cm.stage_ok(f"E{e:g}")
    except Refused as err:
        return None, str(err)
    return e, "characteristics" if rec.get("characteristics_ch1") else "title"


def sample_chemistry(rec: dict) -> str | None:
    txt = " ".join(" ".join(rec.get(k, [])) for k in ("extract_protocol_ch1", "description", "data_processing",
                                                       "library_strategy", "characteristics_ch1", "title"))
    vals = set(_CHEM.findall(txt))
    return f"v{vals.pop()}" if len(vals) == 1 else None


def _rx_ok(spec: dict, title: str, text: str) -> bool:
    s = title if spec.get("on", "text") == "title" else text
    return all(re.search(p, s, re.IGNORECASE) for p in spec.get("all", [])) and \
        not any(re.search(p, s, re.IGNORECASE) for p in spec.get("none", []))


def resolve_arm(ds: dict, rec: dict) -> tuple[str | None, str | None]:
    """The selectors' SUGGESTION for one sample: ('control', gene_hint or None), ('ko', gene) or (None, None). Since
    stage 3 it chooses nothing: only the pinned block does."""
    title = " ".join(rec.get("title", [])).lower()
    text = sample_text(rec)
    where = title if ds.get("exclude_match_on", "text") == "title" else text
    if any(re.search(p, where, re.IGNORECASE) for p in ds.get("exclude_match") or []):
        return None, None
    arms = ds.get("arms") or {}
    ko = arms.get("ko") or []
    if arms.get("control") and _rx_ok(arms["control"], title, text):
        hint = next((k["gene"] for k in ko if re.search(k["all"][0], title, re.IGNORECASE)), None)
        return "control", hint
    for k in ko:
        if _rx_ok(k, title, text):
            return "ko", k["gene"]
    return None, None


def sample_files(ds: dict, rec: dict) -> tuple[list[dict], list[str]]:
    """(per-sample files the dataset's patterns allow, as {name, url}; every other supplementary file name)."""
    pats = [re.compile(p, re.IGNORECASE) for p in ds.get("sample_file_patterns") or []]
    ok, other = [], []
    for f in [v for k in sorted(rec) if k.startswith("supplementary_file") for v in rec[k]]:
        if f.upper() == "NONE":
            continue
        name = Path(urllib.parse.urlparse(f).path).name
        if DENY.search(name) or not any(p.search(name) for p in pats):
            other.append(name)
            continue
        try:
            ok.append({"name": name, "url": https_url(f)})
        except Refused:
            other.append(name)
    return ok, other


def screen_all(ds: dict, metas: dict[str, dict]) -> dict:
    """The held-out screen over every sample and every series-level field of every screened series."""
    records, refusals, held = [], [], False
    extra, cleared = ds.get("screen_extra_refuse"), ds.get("review_cleared")
    pinned = {p["accession"] for p in ds.get("pinned") or []}
    for series, meta in metas.items():
        # every series-level text field (title, summary, overall_design, type, relations ...) is screened alike: a
        # held-out term refuses, and so does a review term not cleared for this series; contact fields are skipped
        texts = []
        for key, vals in meta["series"].items():
            if key.startswith("contact_"):
                continue
            for v in vals:
                texts.append((key, Path(urllib.parse.urlparse(v).path).name if key == "supplementary_file" else v))
        void = any(sc.extended(t) for _, t in texts)
        for key, txt in texts:
            for kind, term in sc.screen_text(series, txt, extra, cleared, void=void):
                held |= kind == "held_out"
                refusals.append(f"{series} {key}: {kind}" + (f" term {term!r}" if kind == "review" else ""))
        for gsm, rec in meta["samples"].items():
            r = sc.screen_record(gsm, rec, extra, cleared)
            records.append({"series": series, **r})
            if r["verdict"] != "pass":
                held |= r["verdict"] == "held_out"
                terms = sorted({h[2] for h in r["hits"] if h[1] == "review"})
                refusals.append(f"{gsm}: {r['verdict']}" + (f" (review terms {terms})" if terms else "")
                                + (f" (clearances void: {r['clearance_void']})" if r["clearance_void"] else ""))
        for gsm, rec in meta["samples"].items():
            if not resolve_arm(ds, rec)[0] and gsm not in pinned:
                continue
            for pat in ds.get("design_refuse") or []:
                if re.search(pat, sample_text(rec), re.IGNORECASE):
                    refusals.append(f"{gsm}: design_refuse {pat!r} matched a suggested or pinned sample")
    return {"n_samples": len(records), "records": records, "refusals": refusals, "held_out": held}


def _row(ds: dict, series: str, gsm: str, rec: dict, scr: dict, in_series: list[str]) -> dict:
    files, other = sample_files(ds, rec)
    arm, gene = resolve_arm(ds, rec)
    e, how = sample_stage(rec)
    return {"accession": gsm, "series": series, "in_series": in_series, "fields": rec,
            "metadata_sha256": cm.metadata_sha256(gsm, rec), "verdict": scr["verdict"],
            "hits": [list(h) for h in scr["hits"]], "clearance_void": scr["clearance_void"],
            "suggested": {"arm": arm, "gene": gene}, "stage": e, "stage_from": how,
            "chemistry": sample_chemistry(rec), "files": files, "other_files": other}


def plan_from_meta(ds: dict, metas: dict[str, dict]) -> dict:
    """The plan for one dataset from parsed SOFT metadata (no network): the screen, one table per series, and, when
    the dataset has pins, their check against this (live) metadata."""
    acc = ds["accession"]
    if acc not in metas:
        raise Refused(f"{acc}: its own series' metadata is missing")
    scr = screen_all(ds, metas)
    by_id = {(r["series"], r["sample"]): r for r in scr["records"]}
    where: dict[str, list[str]] = {}
    for series, meta in metas.items():
        for gsm in meta["samples"]:
            where.setdefault(gsm, []).append(series)
    out = {"accession": acc, "own_series": acc, "screen": {"n_samples": scr["n_samples"]}}
    if scr["held_out"]:
        # nothing of a dataset with a held-out hit is written: accessions, verdicts and the fields hit only
        out["tables"] = {s: [{"accession": g, "verdict": by_id[(s, g)]["verdict"],
                              "fields_hit": sorted({h[0] for h in by_id[(s, g)]["hits"]})}
                             for g in meta["samples"]] for s, meta in metas.items()}
        out["withheld"] = True
        out["screen"]["refusals"] = [WITHHELD]
        return {**out, "status": "refused", "reason": "screen: a held-out hit in the dataset (text withheld)"}
    out["withheld"] = False
    out["tables"] = {s: [_row(ds, s, g, rec, by_id[(s, g)], where[g]) for g, rec in meta["samples"].items()]
                     for s, meta in metas.items()}
    out["screen"]["refusals"] = scr["refusals"]
    out["selector_notes"] = _selector_notes(ds, out["tables"][acc])
    if scr["refusals"]:
        return {**out, "status": "refused", "reason": "screen: " + "; ".join(scr["refusals"][:10])}
    pins = ds.get("pinned") or []
    if pins:
        try:
            out["pinned"] = [{k: v for k, v in s.items() if k != "row"} for s in verify_pins(ds, out, live=True)]
            out["units"] = cm.pinned_units(pins)
        except Refused as e:
            return {**out, "status": "refused", "reason": f"pins: {e}"}
    return {**out, "status": "ok"}


def _selector_notes(ds: dict, rows: list[dict]) -> list[str]:
    notes = []
    for st, want in (ds.get("expect") or {}).items():
        e = cm.stage_ok(st)
        got = {a: sum(1 for r in rows if r["stage"] == e and r["suggested"]["arm"] == a) for a in ("ko", "control")}
        if got != want:
            notes.append(f"E{e:g}: the scouts expected {want}, the selectors suggest {got}")
    return notes


def verify_pins(ds: dict, pl: dict, live: bool = False) -> list[dict]:
    """Each pinned sample checked against the plan's own-series table: present there (not only in another series),
    its stored (or, with live, just-fetched) metadata hashing to the pin, screened clean again now, its parsed stage
    equal to the pin's, its files in its table entry, its chemistry equal to its controls' where the dataset reads
    chemistry. Returns [{accession, arm, genes, stage, files: [{name, url}], metadata_sha256, row}]. Refused on any
    failure, with the reason."""
    acc = ds["accession"]
    pins = ds.get("pinned") or []
    if not pins:
        raise Refused(f"{acc}: the pinned block is empty: nothing is fetched or processed until two readers pin "
                      "samples (README, 'The pinned sample list')")
    if pl.get("withheld") or "tables" not in pl or acc not in pl["tables"]:
        raise Refused(f"{acc}: the plan holds no readable table for its own series")
    own = {r["accession"]: r for r in pl["tables"][acc]}
    elsewhere = {r["accession"]: s for s, rows in pl["tables"].items() if s != acc for r in rows}
    what = "live" if live else "stored"
    per: dict[str, dict] = {}
    for p in pins:
        g = p["accession"]
        if g not in own:
            if g in elsewhere:
                raise Refused(f"{g} belongs to another series ({elsewhere[g]}), not to {acc}")
            raise Refused(f"{g} is not in the plan table of {acc}")
        row = own[g]
        h = cm.metadata_sha256(g, row["fields"])
        if h != p["metadata_sha256"] or h != row["metadata_sha256"]:
            raise Refused(f"{g}: its {what} metadata hashes to {h}, not to the pinned {p['metadata_sha256']}: the "
                          "text the readers read has changed; read and pin again")
        r = sc.screen_record(g, row["fields"], ds.get("screen_extra_refuse"), ds.get("review_cleared"))
        if r["verdict"] != "pass":
            raise Refused(f"{g}: screen verdict {r['verdict']} on the pinned sample: refused")
        e, how = sample_stage(row["fields"])
        if e is None or e != cm.stage_ok(p["stage"]):
            raise Refused(f"{g}: pinned stage {p['stage']} but the record gives {how if e is None else f'E{e:g}'}")
        allowed = {f["name"]: f for f in row["files"]}
        bad = [f for f in p["files"] if f not in allowed or DENY.search(f)]
        if bad:
            raise Refused(f"{g}: pinned files {bad} are not per-sample files of its plan entry")
        for pat in ds.get("design_refuse") or []:
            if re.search(pat, sample_text(row["fields"]), re.IGNORECASE):
                raise Refused(f"{g}: design_refuse {pat!r} matches a pinned sample")
        s = per.setdefault(g, {"accession": g, "arm": p["arm"], "genes": [], "stage": e,
                               "files": [allowed[f] for f in p["files"]], "metadata_sha256": h, "row": row})
        s["genes"].append(p["gene"])
    if ds.get("chemistry"):
        for u in cm.pinned_units(pins):
            for k in u["ko"]:
                ck = sample_chemistry(per[k]["row"]["fields"])
                for c in u["control"]:
                    if per[c]["stage"] != per[k]["stage"]:
                        continue
                    cc = sample_chemistry(per[c]["row"]["fields"])
                    if not ck or ck != cc:
                        raise Refused(f"{u['gene']}: {k} chemistry {ck} and its control {c} chemistry {cc} differ "
                                      "or are unreadable")
    return list(per.values())


def write_tables(pl: dict, mdir: Path) -> list[str]:
    """One human-readable plan table per series (the readers read these)."""
    written = []
    for series, rows in pl["tables"].items():
        lines = [f"# plan table {series} (dataset {pl['accession']}, own series {pl['own_series']}): {len(rows)} "
                 f"samples; plan status {pl['status']}; written {pl.get('created', '')}",
                 "# Field values are verbatim from GEO's family SOFT (surrounding whitespace stripped by the parser).",
                 "# metadata_sha256 is the sha256 of t3x_common.canonical_metadata(accession, fields): README, 'The "
                 "pinned sample list'.",
                 "# 'suggested' is the selectors' guess and chooses nothing; only datasets.yaml's pinned block does."]
        if pl.get("withheld"):
            lines.append("# WITHHELD: a held-out hit in this dataset; no text of it is written. Do not pin.")
            for r in rows:
                lines.append(f"{r['accession']}  verdict {r['verdict']}  fields hit {r['fields_hit']}")
        else:
            for r in rows:
                lines.append("")
                lines.append(f"==== {r['accession']}  verdict {r['verdict']}  metadata_sha256 {r['metadata_sha256']}")
                lines.append(f"  in series: {', '.join(r['in_series'])}")
                stage = "none" if r["stage"] is None else f"E{r['stage']:g}"
                lines.append(f"  suggested (selectors, a guess): arm {r['suggested']['arm']}, gene "
                             f"{r['suggested']['gene']}; parsed stage {stage} ({r['stage_from']}); "
                             f"parsed chemistry {r['chemistry']}")
                lines.append(f"  per-sample files allowed: {[f['name'] for f in r['files']] or 'none'}")
                lines.append(f"  other supplementary files (never fetched): {r['other_files'] or 'none'}")
                lines.append(f"  screen hits: {r['hits'] or 'none'}"
                             + (f"; clearances void: {r['clearance_void']}" if r["clearance_void"] else ""))
                for k, vals in r["fields"].items():
                    for v in vals:
                        lines.append(f"  {k}: {v}")
        p = cm.outside_repo(mdir / f"plan_table_{series}.txt")
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        written.append(str(p))
    return written


# ---------------------------------------------------------------------------------------------------------------------
# Network (only here)


def _get(url: str, dst: Path, max_bytes: int) -> tuple[int, str]:
    url = https_url(url)
    dst = cm.outside_repo(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(dst.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    n = 0
    import hashlib

    h = hashlib.sha256()
    with urllib.request.urlopen(req, timeout=120) as r, open(part, "wb") as f:
        while block := r.read(1 << 20):
            n += len(block)
            if n > max_bytes:
                f.close()
                part.unlink(missing_ok=True)
                raise Refused(f"{url}: over the byte limit ({max_bytes}); nothing kept")
            h.update(block)
            f.write(block)
    part.replace(dst)
    return n, h.hexdigest()


def plan(acc: str, dest: Path, datasets: Path | None = None, metas: dict | None = None) -> dict:
    ds = cm.dataset(acc, datasets)
    if ds.get("repository") != "geo":
        raise Refused(f"{acc}: only GEO is fetched")
    dest = cm.outside_repo(dest)
    mdir = dest / "meta" / acc
    soft_files = {}
    if metas is None:
        metas = {}
        for series in ds["screen_series"]:
            p = mdir / f"{series}_family.soft.gz"
            size, digest = _get(soft_url(series), p, 2 * 1024 ** 3)
            soft_files[series] = {"url": soft_url(series), "file": str(p), "size": size, "sha256": digest}
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
                metas[series] = parse_soft(f.read())
    pl = plan_from_meta(ds, metas)
    pl.update(datasets_sha256=cm.datasets_sha256(datasets), soft_files=soft_files,
              created=time.strftime("%Y-%m-%d %H:%M:%S %Z"), d=ds["d"], tier=ds["tier"], study=ds["study"])
    pl["table_files"] = write_tables(pl, mdir)
    cm.dump_json(pl, mdir / "plan.json")
    return pl


def _load_plan(acc: str, dest: Path, datasets: Path | None) -> dict:
    pp = cm.outside_repo(dest) / "meta" / acc / "plan.json"
    if not pp.is_file():
        raise Refused(f"{pp}: no plan; run plan first")
    pl = json.loads(pp.read_text(encoding="utf-8"))
    if pl.get("status") != "ok":
        raise Refused(f"{acc}: the plan did not pass ({pl.get('reason')})")
    if pl.get("datasets_sha256") != cm.datasets_sha256(datasets):
        raise Refused(f"{acc}: datasets.yaml changed since the plan: plan again")
    return pl


def fetch(acc: str, dest: Path, datasets: Path | None = None, max_bytes: int = MAX_BYTES_DEFAULT) -> dict:
    ds = cm.dataset(acc, datasets)
    if not ds.get("pinned"):
        raise Refused(f"{acc}: the pinned block is empty: nothing is fetched")
    dest = cm.outside_repo(dest)
    pl = _load_plan(acc, dest, datasets)
    pinned = verify_pins(ds, pl)
    pp = dest / "meta" / acc / "plan.json"
    rdir = dest / "raw" / acc
    mf = rdir / "MANIFEST.json"
    man = json.loads(mf.read_text()) if mf.is_file() else {"accession": acc, "files": {}}
    man["files"] = {k: v for k, v in man["files"].items() if k.split("/")[0] in {s["accession"] for s in pinned}}
    total = 0
    for s in pinned:
        for f in s["files"]:
            if DENY.search(f["name"]):
                raise Refused(f"{f['name']}: a denied file type")
            key = f"{s['accession']}/{f['name']}"
            p = rdir / s["accession"] / f["name"]
            rec = man["files"].get(key)
            if rec and p.is_file() and cm.sha256_file(p) == rec["sha256"]:
                total += rec["size"]
                continue
            size, digest = _get(f["url"], p, max_bytes - total)
            total += size
            man["files"][key] = {"url": f["url"], "size": size, "sha256": digest, "gsm": s["accession"],
                                 "fetched": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
            cm.dump_json(man, mf)
    man.update(plan_sha256=cm.sha256_file(pp), d=ds["d"], total_bytes=total, pinned=ds["pinned"],
               pins_sha256=cm.pins_sha256(ds["pinned"]))
    cm.dump_json(man, mf)
    return man


# ---------------------------------------------------------------------------------------------------------------------
# Pinning helpers (never write datasets.yaml)


def pin_template(acc: str, dest: Path, arm: str, gene: str, gsms: list[str] | None = None, title: str | None = None,
                 suggested: bool = False, datasets: Path | None = None) -> str:
    """A ready-to-paste YAML pinned block for the samples one selector matches (exactly one of gsms, title, suggested),
    from the stored plan. Refused for a withheld plan, a held-out sample, a sample whose stage does not parse or
    that has no allowed per-sample file, and when the stored fields no longer hash to the table's value."""
    if arm not in cm.PIN_ARMS:
        raise Refused(f"arm must be exactly one of {cm.PIN_ARMS}")
    if sum(bool(x) for x in (gsms, title, suggested)) != 1:
        raise Refused("give exactly one selector: --gsm, --title or --suggested")
    cm.dataset(acc, datasets)
    pp = cm.outside_repo(dest) / "meta" / acc / "plan.json"
    if not pp.is_file():
        raise Refused(f"{pp}: no plan; run plan first")
    pl = json.loads(pp.read_text(encoding="utf-8"))
    if pl.get("withheld"):
        raise Refused(f"{acc}: a held-out hit in the dataset: nothing of it is pinned")
    rows = pl["tables"][acc]
    if gsms:
        missing = sorted(set(gsms) - {r["accession"] for r in rows})
        if missing:
            raise Refused(f"{missing} are not in the plan table of {acc}")
        pick = [r for r in rows if r["accession"] in gsms]
    elif title:
        rx = re.compile(title, re.IGNORECASE)
        pick = [r for r in rows if any(rx.search(t) for t in r["fields"].get("title", []))]
    else:
        pick = [r for r in rows if r["suggested"]["arm"] == arm
                and (arm == "control" or r["suggested"]["gene"] == gene)]
    if not pick:
        raise Refused("the selector matches no sample")
    out = [f"# pin-template {acc}: arm {arm}, gene {gene}; {len(pick)} sample(s); from {pp}",
           "# READ each sample's entry in its plan table before pasting. The arm and gene are what you asserted.",
           "# A control names the knockout gene it controls; pin it once per gene it serves."]
    for r in pick:
        if r["verdict"] == "held_out":
            raise Refused(f"{r['accession']}: held out: never pinned")
        h = cm.metadata_sha256(r["accession"], r["fields"])
        if h != r["metadata_sha256"]:
            raise Refused(f"{r['accession']}: the stored fields hash to {h}, not to the table's value: plan again")
        if r["stage"] is None:
            raise Refused(f"{r['accession']}: no single stage parses ({r['stage_from']}): it cannot be pinned")
        if not r["files"]:
            raise Refused(f"{r['accession']}: no per-sample file the dataset's patterns allow")
        if r["verdict"] == "review":
            out.append(f"# {r['accession']} has review terms {sorted({h[2] for h in r['hits'] if h[1] == 'review'})}"
                       ": pin it only with review_cleared entries the readers agree")
        files = ", ".join(json.dumps(f["name"]) for f in r["files"])
        out += [f"      - accession: {r['accession']}",
                f"        arm: {arm}",
                f"        gene: {json.dumps(gene)}",
                f"        stage: E{r['stage']:g}",
                f"        files: [{files}]",
                f'        metadata_sha256: "{h}"']
    return "\n".join(out) + "\n"


def check_gate(accs: list[str], datasets: Path | None = None) -> list[str]:
    """The launcher's 'pinned' gate: [] if every accession is cleared and pinned and the committed dataset list
    hashes to its pin; else the reasons."""
    bad = []
    for a in accs:
        try:
            ds = cm.dataset(a, datasets)
            if not ds.get("pinned"):
                bad.append(f"{a}: the pinned block is empty")
        except Refused as e:
            bad.append(f"{a}: {e}")
    return bad


def main(argv: list[str] | None = None) -> int:
    import argparse
    import hashlib

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="step", required=True)
    for name in ("plan", "fetch"):
        s = sub.add_parser(name)
        s.add_argument("accession")
        s.add_argument("--dest", required=True, type=Path, help="folder outside the repository")
        if name == "fetch":
            s.add_argument("--max-gb", type=float, default=40.0)
    t = sub.add_parser("pin-template")
    t.add_argument("accession")
    t.add_argument("--dest", required=True, type=Path)
    t.add_argument("--arm", required=True)
    t.add_argument("--gene", required=True)
    t.add_argument("--gsm", action="append")
    t.add_argument("--title")
    t.add_argument("--suggested", action="store_true")
    g = sub.add_parser("check-pins")
    g.add_argument("accessions", nargs="+")
    sub.add_parser("datasets-sha")
    a = ap.parse_args(argv)
    try:
        if a.step == "plan":
            pl = plan(a.accession, a.dest)
            print(f"plan {pl['status']}: {a.accession}: {pl['screen']['n_samples']} samples screened; tables "
                  f"{pl['table_files']}")
            if pl["status"] != "ok":
                print(f"refused: {a.accession}: {pl['reason']}")
                return 2
            for note in pl.get("selector_notes") or []:
                print(f"  selector note: {note}")
            if pl.get("pinned"):
                print(f"  pins verified against the live metadata: {[s['accession'] for s in pl['pinned']]}; units "
                      f"{[u['gene'] for u in pl['units']]}")
            else:
                print("  no pins yet: read the plan tables, then pin (README, 'The pinned sample list')")
        elif a.step == "fetch":
            man = fetch(a.accession, a.dest, max_bytes=int(a.max_gb * 1024 ** 3))
            print(f"fetched: {a.accession}: {len(man['files'])} files of the pinned samples, {man['total_bytes']} "
                  "bytes")
        elif a.step == "pin-template":
            print(pin_template(a.accession, a.dest, a.arm, a.gene, a.gsm, a.title, a.suggested), end="")
        elif a.step == "check-pins":
            bad = check_gate(a.accessions)
            for b in bad:
                print(f"refused: {b}")
            if bad:
                return 2
            print(f"pinned gate passed: {a.accessions}; datasets.yaml sha256 {cm.DATASETS_SHA256}")
        else:
            got = hashlib.sha256(cm.DATASETS.read_bytes()).hexdigest()
            print(f"datasets.yaml sha256 {got}\nt3x_common.DATASETS_SHA256 {cm.DATASETS_SHA256}")
            if got != cm.DATASETS_SHA256:
                print("refused: they differ: put the first value into t3x_common.DATASETS_SHA256 in the same commit")
                return 2
    except Refused as e:
        print(f"refused: {e}")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
