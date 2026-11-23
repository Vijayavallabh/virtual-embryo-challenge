#!/usr/bin/env python
"""What the two readers need for checks (a) to (d) of row T72's cell 3, computed on the barcodes and features files
`fetch.py fetch-meta` fetched into ~/v/scratch/t72_screen/screen/<ACC>/ (candidate MX2s; Task 1 only), BEFORE any
matrix is fetched or opened. The readers decide; this prints counts and yes/no facts only, never a cell, a barcode or a
gene, and the yaml block the main session fills from both readers' hand-backs.

  (a) every barcode ends in -1 (no aggr): per barcodes file, its line count (as `zcat | wc -l`) and the count of
      barcodes not ending in -1;
  (b) X019 (multiome): the two features files have different sha256; at least one Peaks interval of each file overlaps
      an interval of the other file without being identical to it (per-library peak calls; subsets of one merged peak
      set never partly overlap); the two files' Gene Expression rows are the same features in the same order;
  (c) no hashing or multiplexing feature: the feature_type counts per features file, the rows of a type outside the
      atlas's assay, hashtag-like names, the column counts, barcode lines with more than one column;
  (d) X020: each barcodes file has the line count of the 10x 3' v2 barcode whitelist (737,280).
The same functions run again on the fetched files at the row's load (load.py), where a failure is a failed (p1).

  uv run --frozen --no-sync python -I -B experiments/t1_ext/screens.py GSE244783 --dest ~/v/scratch/t72_screen
Prints the facts and the paste block; exit 0 when every check of the atlas says yes, 2 otherwise (or refused).
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import itertools
import json
import re
import sys
from collections.abc import Iterator, Mapping
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import t1x_common as cm  # noqa: E402

Refused = cm.Refused
SUFFIX = "-1"
WHITELIST_V2 = 737280  # cell 3 (d): the 10x 3' v2 barcode whitelist's line count
GE = "Gene Expression"


def _open(path: Path):
    return gzip.open(path, "rt", newline="") if Path(path).name.endswith(".gz") else open(path, newline="")


def _lines(path: Path) -> Iterator[list[str]]:
    with _open(path) as f:
        for ln in f:
            ln = ln.rstrip("\n").rstrip("\r")
            if ln.strip():
                yield ln.split("\t")


def newline_count(path: Path) -> int:
    """The line count `zcat | wc -l` prints: the number of newline characters."""
    n = 0
    with _open(path) as f:
        while block := f.read(1 << 20):
            n += block.count("\n")
    return n


def barcode_facts(path: Path) -> dict:
    """Counts only: newline count, barcodes, barcodes not ending in -1, lines with more than one column."""
    n = bad = extra = 0
    for parts in _lines(path):
        n += 1
        bad += not parts[0].endswith(SUFFIX)
        extra += len(parts) > 1
    return {"wc_l": newline_count(path), "barcodes": n, "not_suffix": bad, "lines_with_extra_columns": extra}


def _interval(parts: list[str]) -> tuple[str, int, int] | None:
    if len(parts) >= 6:
        try:
            return parts[3], int(parts[4]), int(parts[5])
        except ValueError:
            pass
    m = re.fullmatch(r"([^:\s]+)[:_-](\d+)[-_](\d+)", parts[0])
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


def feature_facts(path: Path, rules: Mapping) -> dict:
    """Counts per feature type, the column counts seen, hashtag-like names; kept in memory only (never written): the
    Peaks intervals and the Gene Expression rows (id, name) in order, for check (b)."""
    rx = re.compile(rules["hashing_names"], re.IGNORECASE)
    types: dict[str, int] = {}
    widths: dict[int, int] = {}
    hashing = unparsed = 0
    peaks: list[tuple[str, int, int]] = []
    ge: list[tuple[str, str]] = []
    for parts in _lines(path):
        widths[len(parts)] = widths.get(len(parts), 0) + 1
        t = parts[2] if len(parts) >= 3 else GE
        types[t] = types.get(t, 0) + 1
        hashing += any(rx.search(x.strip()) for x in parts[:2])
        if t == "Peaks":
            iv = _interval(parts)
            if iv is None:
                unparsed += 1
            else:
                peaks.append(iv)
        elif t == GE:
            ge.append((parts[0], parts[1] if len(parts) > 1 else ""))
    return {"types": types, "widths": {str(k): v for k, v in sorted(widths.items())}, "hashing_names": hashing,
            "peaks_rows": len(peaks), "peaks_unparsed": unparsed, "_peaks": peaks, "_ge": ge}


def partial_overlap(a: list[tuple[str, int, int]], b: list[tuple[str, int, int]]) -> bool:
    """Whether some interval of `a` overlaps an interval of `b` (same chromosome, start < other's end and other's start
    < end) without being identical to it. The relation is symmetric, so one answer serves both files."""
    by: dict[str, list[tuple[int, int]]] = {}
    for c, s, e in b:
        by.setdefault(c, []).append((s, e))
    for c in by:
        by[c].sort()
    span = {c: max(e - s for s, e in v) for c, v in by.items()}
    starts = {c: [s for s, _ in v] for c, v in by.items()}
    for c, s, e in a:
        if c not in by:
            continue
        lo = bisect.bisect_left(starts[c], s - span[c])
        hi = bisect.bisect_left(starts[c], e)
        for s2, e2 in by[c][lo:hi]:
            if s < e2 and s2 < e and (s, e) != (s2, e2):
                return True
    return False


def check_a(bc: Mapping[str, Mapping]) -> dict:
    n = sum(x["barcodes"] for x in bc.values())
    bad = sum(x["not_suffix"] for x in bc.values())
    return {"passed": n > 0 and bad == 0 and all(x["barcodes"] > 0 for x in bc.values()),
            "barcodes_not_ending_in_-1": bad}


def check_b(ft: Mapping[str, Mapping], fsha: Mapping[str, str]) -> dict:
    """(b) on two or more features files (each pair must hold every clause)."""
    gs = sorted(ft)
    pairs = list(itertools.combinations(gs, 2))
    out = {"pairs": {}}
    ok = len(gs) >= 2
    for x, y in pairs:
        p = {"different_sha256": fsha[x] != fsha[y],
             "peaks_partly_overlap": partial_overlap(ft[x]["_peaks"], ft[y]["_peaks"]),
             "gene_expression_rows_same_in_order": ft[x]["_ge"] == ft[y]["_ge"] and bool(ft[x]["_ge"]),
             "peaks_rows": [ft[x]["peaks_rows"], ft[y]["peaks_rows"]],
             "peaks_rows_unparsed": [ft[x]["peaks_unparsed"], ft[y]["peaks_unparsed"]]}
        p["passed"] = bool(p["different_sha256"] and p["peaks_partly_overlap"]
                           and p["gene_expression_rows_same_in_order"] and not any(p["peaks_rows_unparsed"])
                           and all(p["peaks_rows"]))
        ok &= p["passed"]
        out["pairs"][f"{x}|{y}"] = p
    out["passed"] = bool(ok)
    return out


def check_c(bc: Mapping, ft: Mapping, rules: Mapping, multiome: bool) -> dict:
    allowed_types = set(rules["feature_types"]["multiome" if multiome else "tenx"])
    allowed_widths = {"3", "6"} if multiome else {"2", "3"}
    extra_bc = sum(x["lines_with_extra_columns"] for x in bc.values())
    other_types = {g: sum(v for t, v in x["types"].items() if t not in allowed_types) for g, x in ft.items()}
    bad_widths = {g: sum(v for w, v in x["widths"].items() if w not in allowed_widths) for g, x in ft.items()}
    hashing = {g: x["hashing_names"] for g, x in ft.items()}
    ok = extra_bc == 0 and not any(other_types.values()) and not any(bad_widths.values()) and \
        not any(hashing.values())
    return {"passed": bool(ok), "barcode_lines_with_extra_columns": extra_bc,
            "feature_rows_of_other_types": other_types, "feature_rows_of_other_widths": bad_widths,
            "hashtag_like_feature_names": hashing}


def check_d(bc: Mapping[str, Mapping]) -> dict:
    n = {g: x["wc_l"] for g, x in bc.items()}
    return {"passed": bool(n) and all(v == WHITELIST_V2 for v in n.values()), "wc_l": n, "whitelist": WHITELIST_V2}


CHECKS = {"suffix": "a", "own_peaks": "b", "no_hashing": "c", "whitelist_sized": "d"}


def run_checks(a: Mapping, rules: Mapping, files: Mapping[str, Mapping[str, Path]],
               sha: Mapping[str, Mapping[str, str]]) -> dict:
    """The atlas's checks on {GSM: {barcodes: path, features: path}}; sha: {GSM: {role: sha256}}."""
    bc = {g: barcode_facts(f["barcodes"]) for g, f in sorted(files.items())}
    ft = {g: feature_facts(f["features"], rules) for g, f in sorted(files.items())}
    out = {"per_file": {g: {"barcodes": {k: bc[g][k] for k in ("wc_l", "not_suffix")},
                            "features": {"types": ft[g]["types"]}} for g in sorted(files)}}
    for name in a["screens"]:
        if name == "suffix":
            out[name] = check_a(bc)
        elif name == "own_peaks":
            out[name] = check_b(ft, {g: sha[g]["features"] for g in ft})
        elif name == "no_hashing":
            out[name] = check_c(bc, ft, rules, bool(a.get("multiome")))
        elif name == "whitelist_sized":
            out[name] = check_d(bc)
        else:
            raise Refused(f"unknown screen {name!r}")
    out["passed"] = all(out[n]["passed"] for n in a["screens"])
    return out


def screen(acc: str, dest: Path, cfg: Mapping | None = None, screen_root: Path | None = None) -> dict:
    """The checks on the screens' fetch of one atlas: every file re-hashed against FETCH_META.json first; no matrix."""
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    if a["status"] not in cm.FETCHABLE or not a.get("screens"):
        raise Refused(f"{acc}: no screens for this atlas")
    fetch = cm.sibling("fetch")
    root = fetch.screen_dir_ok(dest, screen_root)
    sdir = root / "screen" / acc
    mp = sdir / "FETCH_META.json"
    if not mp.is_file():
        raise Refused(f"{acc}: no FETCH_META.json: run `fetch.py fetch-meta` first")
    fm = json.loads(mp.read_text())
    if fm.get("accession") != acc or set(fm.get("roles") or []) != set(cm.SCREEN_ROLES):
        raise Refused(f"{acc}: FETCH_META.json is not the screens' fetch of this atlas (barcodes and features only)")
    files, sha = {}, {}
    for g, s in sorted(fm["samples"].items()):
        if set(s["files"]) != set(cm.SCREEN_ROLES):
            raise Refused(f"{g}: the screens' fetch holds roles {sorted(s['files'])}")
        files[g], sha[g] = {}, {}
        for role, f in s["files"].items():
            p = sdir / g / f["name"]
            if not p.is_file() or cm.sha256_file(p) != f["sha256"]:
                raise Refused(f"{g}: {f['name']} is missing or not the file FETCH_META.json records")
            files[g][role], sha[g][role] = p, f["sha256"]
    res = run_checks(a, cfg["screen_rules"], files, sha)
    return {"accession": acc, "fetch_meta_sha256": cm.sha256_file(mp), "checks": res,
            "files": {g: {r: {"name": fm["samples"][g]["files"][r]["name"], "sha256": sha[g][r]} for r in sha[g]}
                      for g in sha}}


def paste_block(res: Mapping) -> str:
    """The yaml block the main session fills from BOTH readers' hand-backs (never from this output alone)."""
    out = ["    screen:",
           "      verdict: READERS_VERDICT   # kept only if both readers said yes to every point; otherwise set the",
           "                                 # atlas's status: dropped, verdict: dropped and failed: [points]",
           "      files:"]
    for g, f in sorted(res["files"].items()):
        out.append(f"        {g}:")
        for r in cm.SCREEN_ROLES:
            out.append(f'          {r}: {{name: {f[r]["name"]}, sha256: "{f[r]["sha256"]}"}}')
    out += ["      read_by: [READER_1, READER_2]", '      date: "YYYY-MM-DD"']
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accession")
    ap.add_argument("--dest", required=True, type=Path)
    a = ap.parse_args(argv)
    try:
        res = screen(a.accession, a.dest)
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr)
        return 2
    ch = res["checks"]
    summary = {CHECKS[k]: ("yes" if ch[k]["passed"] else "no") for k in ch if k in CHECKS}
    print(json.dumps({"accession": res["accession"], "screens": summary, "checks": ch,
                      "fetch_meta_sha256": res["fetch_meta_sha256"]}, indent=1))
    print("# the screen block for experiments/t1_ext/atlases.yaml (fill from both readers' hand-backs):")
    print(paste_block(res), end="")
    return 0 if ch["passed"] else 2


if __name__ == "__main__":
    sys.exit(main())
