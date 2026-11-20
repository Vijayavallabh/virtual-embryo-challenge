#!/usr/bin/env python
"""Fetch one cleared outside-knockout dataset from GEO into a folder OUTSIDE the repository.

Two steps, both refused for an accession that is not in datasets.yaml with status cleared:

  plan  ACC --dest D   downloads only the series' family SOFT metadata (D/meta/<ACC>/<series>_family.soft.gz, sha256
                       recorded), screens EVERY sample of every series in screen_series with t3x_screen (any held-out
                       hit, an unscreenable record or an unresolved review term refuses the dataset), applies
                       design_refuse, resolves each sample's arm, gene, stage (capped at E9.5) and chemistry, checks the
                       expected counts, and writes D/meta/<ACC>/plan.json: the chosen samples and their per-sample
                       supplementary files (never series-level files, never FASTQ, ATAC, fragments or BAM).
  fetch ACC --dest D   downloads exactly the files plan.json lists to D/raw/<ACC>/<GSM>/, each written to a .part file
                       and renamed when complete, and records each file's url, size and sha256 in
                       D/raw/<ACC>/MANIFEST.json.
                       Nothing is unpacked. It refuses unless plan.json passed and datasets.yaml still has the sha256
                       the plan recorded.

Requests go only to ftp.ncbi.nlm.nih.gov over https, with a generic user agent and nothing that names a team, a
person or a challenge. No count is read here. Exit status 2 on a refusal.

  PATH=$HOME/.local/bin:$PATH uv run --frozen --no-sync python -I -B experiments/t3_ext/fetch.py plan GSE168768 \\
      --dest ~/v/vec-ext/t3
"""

from __future__ import annotations

import gzip
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


# ---------------------------------------------------------------------------------------------------------------------
# GEO metadata


def soft_url(series: str) -> str:
    if not re.fullmatch(r"GSE\d+", series):
        raise Refused(f"{series}: only GEO series are fetched in stage 1")
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
    """{'series': {key: [values]}, 'samples': {GSM: {key: [values]}}} from a family SOFT text."""
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
    """(arm, gene) for one sample: ('control', gene_hint or None), ('ko', gene) or (None, None)."""
    title = " ".join(rec.get("title", [])).lower()
    text = sample_text(rec)
    if any(re.search(p, text, re.IGNORECASE) for p in ds.get("exclude_match") or []):
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


def screen_all(ds: dict, metas: dict[str, dict]) -> dict:
    """The held-out screen over every sample (and the series title and design) of every screened series."""
    records, refusals = [], []
    for series, meta in metas.items():
        for key in ("title", "overall_design"):
            for v in meta["series"].get(key, []):
                got = sc.match(v)
                if got and got[0] == "held_out":
                    refusals.append(f"{series} {key}: held-out term {got[1]!r}")
        for v in meta["series"].get("summary", []):
            got = sc.match(v)
            if got:
                refusals.append(f"{series} summary: {got[0]} term {got[1]!r} (refuse for review)")
        for v in meta["series"].get("supplementary_file", []):
            got = sc.match(Path(v).name)
            if got:
                refusals.append(f"{series} series file name: {got[0]} term {got[1]!r}")
        for gsm, rec in meta["samples"].items():
            r = sc.screen_record(gsm, rec, ds.get("screen_extra_refuse"), ds.get("review_cleared"))
            records.append({"series": series, **r})
            if r["verdict"] != "pass":
                refusals.append(f"{gsm}: {r['verdict']} ({sorted({h[2] for h in r['hits']})})")
        for gsm, rec in meta["samples"].items():
            if not resolve_arm(ds, rec)[0]:
                continue
            for pat in ds.get("design_refuse") or []:
                if re.search(pat, sample_text(rec), re.IGNORECASE):
                    refusals.append(f"{gsm}: design_refuse {pat!r} matched a sample of an arm")
    return {"n_samples": len(records), "records": records, "refusals": refusals}


def plan_from_meta(ds: dict, metas: dict[str, dict]) -> dict:
    """The plan for one dataset from parsed SOFT metadata (no network): screen, resolve arms, choose files."""
    acc = ds["accession"]
    scr = screen_all(ds, metas)
    out = {"accession": acc, "screen": {"n_samples": scr["n_samples"], "refusals": scr["refusals"],
                                        "hits": [r for r in scr["records"] if r["hits"]]}}
    if scr["refusals"]:
        return {**out, "status": "refused", "reason": "screen: " + "; ".join(scr["refusals"][:10])}
    used = {cm.stage_ok(s) for s in ds["stages_used"]}
    pats = [re.compile(p, re.IGNORECASE) for p in ds.get("sample_file_patterns") or []]
    own = metas.get(acc) or next(iter(metas.values()))
    chosen, problems = [], []
    for gsm, rec in own["samples"].items():
        arm, gene = resolve_arm(ds, rec)
        if arm is None:
            continue
        e, how = sample_stage(rec)
        if e is None:
            problems.append(f"{gsm} ({arm} {gene}): stage unresolved ({how})")
            continue
        if e not in used:
            continue
        files = []
        for f in [v for k in sorted(rec) if k.startswith("supplementary_file") for v in rec[k]]:
            if f.upper() == "NONE":
                continue
            name = Path(urllib.parse.urlparse(f).path).name
            if DENY.search(name) or not any(p.search(name) for p in pats):
                continue
            files.append({"name": name, "url": https_url(f)})
        if not files:
            problems.append(f"{gsm} ({arm} {gene}): no per-sample file matches the dataset's patterns")
            continue
        chosen.append({"gsm": gsm, "title": "; ".join(rec.get("title", [])), "arm": arm, "gene": gene,
                       "stage": e, "stage_from": how, "chemistry": sample_chemistry(rec) if ds.get("chemistry")
                       else None, "files": files, "record": rec})
    # pair knockouts with controls of the same stage (and chemistry where required)
    kos = [c for c in chosen if c["arm"] == "ko"]
    ctls = [c for c in chosen if c["arm"] == "control"]
    units = {}
    for c in kos:
        same = [k for k in ctls if k["stage"] == c["stage"]
                and (not ds.get("chemistry") or (c["chemistry"] and k["chemistry"] == c["chemistry"]))]
        if ds.get("chemistry") and not c["chemistry"]:
            problems.append(f"{c['gsm']} ({c['gene']}): chemistry unreadable")
            continue
        own_ctl = [k for k in same if k["gene"] == c["gene"]]
        pick = own_ctl or [k for k in same if k["gene"] is None]
        if not pick:
            problems.append(f"{c['gsm']} ({c['gene']}, E{c['stage']:g}): no matched control")
            continue
        u = units.setdefault(c["gene"], {"gene": c["gene"], "ko": [], "control": []})
        u["ko"].append(c["gsm"])
        u["control"] += [k["gsm"] for k in pick if k["gsm"] not in u["control"]]
    for st, want in (ds.get("expect") or {}).items():
        e = cm.stage_ok(st)
        got = {a: sum(1 for c in chosen if c["stage"] == e and c["arm"] == a) for a in ("ko", "control")}
        if got != want:
            problems.append(f"E{e:g}: expected {want}, found {got}")
    if problems:
        return {**out, "status": "refused", "reason": "selection: " + "; ".join(problems[:10]),
                "chosen": [{k: v for k, v in c.items() if k != "record"} for c in chosen]}
    if not units:
        return {**out, "status": "refused", "reason": "no knockout arm with a matched control"}
    keep = {g for u in units.values() for g in u["ko"] + u["control"]}
    return {**out, "status": "ok", "units": list(units.values()),
            "samples": [c for c in chosen if c["gsm"] in keep]}


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
        raise Refused(f"{acc}: only GEO is fetched in stage 1")
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
    pl.update(datasets_sha256=cm.sha256_file(Path(datasets or cm.DATASETS)), soft_files=soft_files,
              created=time.strftime("%Y-%m-%d %H:%M:%S %Z"), d=ds["d"], tier=ds["tier"], study=ds["study"])
    cm.dump_json(pl, mdir / "plan.json")
    return pl


def fetch(acc: str, dest: Path, datasets: Path | None = None, max_bytes: int = MAX_BYTES_DEFAULT) -> dict:
    ds = cm.dataset(acc, datasets)
    dest = cm.outside_repo(dest)
    import json

    pp = dest / "meta" / acc / "plan.json"
    if not pp.is_file():
        raise Refused(f"{pp}: no plan; run plan first")
    pl = json.loads(pp.read_text())
    if pl.get("status") != "ok":
        raise Refused(f"{acc}: the plan did not pass ({pl.get('reason')})")
    if pl.get("datasets_sha256") != cm.sha256_file(Path(datasets or cm.DATASETS)):
        raise Refused(f"{acc}: datasets.yaml changed since the plan: plan again")
    rdir = dest / "raw" / acc
    mf = rdir / "MANIFEST.json"
    man = json.loads(mf.read_text()) if mf.is_file() else {"accession": acc, "files": {}}
    total = 0
    for s in pl["samples"]:
        for f in s["files"]:
            if DENY.search(f["name"]):
                raise Refused(f"{f['name']}: a denied file type")
            key = f"{s['gsm']}/{f['name']}"
            p = rdir / s["gsm"] / f["name"]
            rec = man["files"].get(key)
            if rec and p.is_file() and cm.sha256_file(p) == rec["sha256"]:
                total += rec["size"]
                continue
            size, digest = _get(f["url"], p, max_bytes - total)
            total += size
            man["files"][key] = {"url": f["url"], "size": size, "sha256": digest, "gsm": s["gsm"],
                                 "fetched": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
            cm.dump_json(man, mf)
    man.update(plan_sha256=cm.sha256_file(pp), d=ds["d"], total_bytes=total)
    cm.dump_json(man, mf)
    return man


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=("plan", "fetch"))
    ap.add_argument("accession")
    ap.add_argument("--dest", required=True, type=Path, help="folder outside the repository")
    ap.add_argument("--datasets", type=Path, default=None)
    ap.add_argument("--max-gb", type=float, default=40.0)
    a = ap.parse_args(argv)
    try:
        if a.step == "plan":
            pl = plan(a.accession, a.dest, a.datasets)
            if pl["status"] != "ok":
                print(f"refused: {a.accession}: {pl['reason']}")
                return 2
            n = {arm: sum(1 for s in pl["samples"] if s["arm"] == arm) for arm in ("ko", "control")}
            print(f"plan ok: {a.accession}: {pl['screen']['n_samples']} samples screened; chosen {n}; "
                  f"units {[u['gene'] for u in pl['units']]}")
            for s in pl["samples"]:
                print(f"  {s['gsm']} {s['arm']:7s} {s['gene'] or '-':14s} E{s['stage']:g} {s['title']}")
        else:
            man = fetch(a.accession, a.dest, a.datasets, int(a.max_gb * 1024 ** 3))
            print(f"fetched: {a.accession}: {len(man['files'])} files, {man['total_bytes']} bytes")
    except Refused as e:
        print(f"refused: {e}")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
