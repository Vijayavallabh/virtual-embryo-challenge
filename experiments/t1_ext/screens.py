#!/usr/bin/env python
"""The specification's SCREENS (a) to (d) on one atlas's fetched barcode and feature files, BEFORE any matrix is
fetched or opened (candidate MX2s, ruling row R27; Task 1 only).

Input: D/raw/<ACC>/HEADS.json and the files it lists (written by `fetch.py heads`; every file is re-hashed first). No
matrix is read: the matrix role is never in the heads. Each check returns yes or no with counts, never a cell, a
barcode or a gene; the checks an atlas runs are its `screens` list in atlases.yaml, their constants are
atlases.yaml `screen_rules`.
  suffix          (a) every barcode of every sample ends in -1 (cellranger aggr renumbers libraries -2, -3, ...);
  own_peaks       (b) Multiome only: the samples' features files have different sha256 and different sets of Peaks
                  rows (cellranger-arc count per library gives each library its own peak set; aggr gives one shared
                  set), every sample has Peaks rows and Gene Expression rows. The Peaks rows' ids are hashed in memory
                  to compare the sets and never written; the counts of Gene Expression rows are reported (only those
                  rows are ever read as expression, by load.py);
  no_hashing      (c) no hashing or pooling field: every barcode line has one column (no per-cell sample, stage or tag
                  column); every features line has 2 or 3 columns (Multiome: 6, or 3); every feature type is one of
                  the atlas's assay (screen_rules.feature_types); no feature id or name looks like a hashtag, CMO or
                  antibody tag (screen_rules.hashing_names);
  whitelist_sized (d) every sample's barcode file has one of the 10x whitelist sizes (screen_rules.whitelist_sizes):
                  the matrices are unfiltered and the team's own cell call applies.
Output: D/screens/<ACC>/screen.json (the checks, the heads' sha256, passed yes or no, tasks ["T1"]), and on stdout the
summary and the line the main session pastes into experiments/t1_ext/screens.yaml (kept when every check said yes,
dropped otherwise; a dropped member's matrices are never fetched). Exit 0 when every check passed, 2 otherwise.

  uv run --frozen --no-sync python -I -B experiments/t1_ext/screens.py GSE244783 --dest ~/v/vec-ext/t1
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
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


def _lines(path: Path) -> Iterator[list[str]]:
    op = gzip.open if path.name.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        for ln in f:
            ln = ln.rstrip("\n").rstrip("\r")
            if ln.strip():
                yield ln.split("\t")


def barcode_facts(path: Path) -> dict:
    """Counts only: barcodes, barcodes not ending in -1, lines with more than one column."""
    n = bad = extra = 0
    for parts in _lines(path):
        n += 1
        bad += not parts[0].endswith(SUFFIX)
        extra += len(parts) > 1
    return {"barcodes": n, "not_suffix": bad, "lines_with_extra_columns": extra}


def feature_facts(path: Path, rules: Mapping) -> dict:
    """Counts per feature type, the column counts seen, hashtag-like names, and (Multiome) a hash of the set of
    Peaks ids, with the ids kept in memory only for the cross-sample comparison."""
    rx = re.compile(rules["hashing_names"], re.IGNORECASE)
    types: dict[str, int] = {}
    widths: dict[int, int] = {}
    hashing = 0
    peaks: set[str] = set()
    for parts in _lines(path):
        widths[len(parts)] = widths.get(len(parts), 0) + 1
        t = parts[2] if len(parts) >= 3 else "Gene Expression"
        types[t] = types.get(t, 0) + 1
        hashing += any(rx.search(x.strip()) for x in parts[:2])
        if t == "Peaks":
            peaks.add(parts[0])
    h = hashlib.sha256("\n".join(sorted(peaks)).encode()).hexdigest() if peaks else None
    return {"types": types, "widths": {str(k): v for k, v in sorted(widths.items())}, "hashing_names": hashing,
            "peaks_rows": len(peaks), "peaks_set_sha256": h, "_peaks": peaks}


def run_screens(a: Mapping, rules: Mapping, heads: Mapping, rdir: Path) -> dict:
    """The checks the atlas lists, on the heads' files (each re-hashed against HEADS.json first)."""
    multiome = bool(a.get("multiome"))
    bc, ft, fsha = {}, {}, {}
    for gsm, s in sorted(heads["samples"].items()):
        for role in ("barcodes", "features"):
            f = s["files"][role]
            p = rdir / gsm / f["name"]
            if not p.is_file() or cm.sha256_file(p) != f["sha256"]:
                raise Refused(f"{gsm}: {f['name']} is missing or not the file HEADS.json records")
        bc[gsm] = barcode_facts(rdir / gsm / s["files"]["barcodes"]["name"])
        ft[gsm] = feature_facts(rdir / gsm / s["files"]["features"]["name"], rules)
        fsha[gsm] = s["files"]["features"]["sha256"]
    checks: dict[str, dict] = {}
    for name in a["screens"]:
        if name == "suffix":
            n = sum(x["barcodes"] for x in bc.values())
            bad = sum(x["not_suffix"] for x in bc.values())
            checks[name] = {"passed": n > 0 and bad == 0 and all(x["barcodes"] > 0 for x in bc.values()),
                            "barcodes": n, "not_ending_in_-1": bad}
        elif name == "own_peaks":
            gs = sorted(ft)
            pairs = list(itertools.combinations(gs, 2))
            same_file = sum(fsha[x] == fsha[y] for x, y in pairs)
            same_peaks = sum(ft[x]["peaks_set_sha256"] == ft[y]["peaks_set_sha256"] for x, y in pairs)
            shared = {f"{x}|{y}": len(ft[x]["_peaks"] & ft[y]["_peaks"]) for x, y in pairs}
            ok = (len(gs) >= 2 and same_file == 0 and same_peaks == 0
                  and all(ft[g]["peaks_rows"] > 0 and ft[g]["types"].get("Gene Expression", 0) > 0 for g in gs))
            checks[name] = {"passed": bool(ok), "samples": len(gs), "pairs_with_identical_features_file": same_file,
                            "pairs_with_identical_peaks_rows": same_peaks,
                            "peaks_rows": {g: ft[g]["peaks_rows"] for g in gs},
                            "peaks_rows_shared_per_pair": shared,
                            "gene_expression_rows": {g: ft[g]["types"].get("Gene Expression", 0) for g in gs}}
        elif name == "no_hashing":
            allowed_types = set(rules["feature_types"]["multiome" if multiome else "tenx"])
            allowed_widths = {"3", "6"} if multiome else {"2", "3"}
            extra_bc = sum(x["lines_with_extra_columns"] for x in bc.values())
            other_types = {g: sum(v for t, v in x["types"].items() if t not in allowed_types) for g, x in ft.items()}
            bad_widths = {g: sum(v for w, v in x["widths"].items() if w not in allowed_widths) for g, x in ft.items()}
            hashing = {g: x["hashing_names"] for g, x in ft.items()}
            ok = extra_bc == 0 and not any(other_types.values()) and not any(bad_widths.values()) and \
                not any(hashing.values())
            checks[name] = {"passed": bool(ok), "barcode_lines_with_extra_columns": extra_bc,
                            "feature_rows_of_other_types": other_types, "feature_rows_of_other_widths": bad_widths,
                            "hashtag_like_feature_names": hashing}
        elif name == "whitelist_sized":
            sizes = {int(x) for x in rules["whitelist_sizes"]}
            n = {g: x["barcodes"] for g, x in bc.items()}
            checks[name] = {"passed": bool(n) and all(v in sizes for v in n.values()), "barcodes_per_sample": n,
                            "whitelist_sizes": sorted(sizes)}
        else:
            raise Refused(f"unknown screen {name!r}")
    return checks


def screen(acc: str, dest: Path, cfg: Mapping | None = None) -> dict:
    cfg = cfg or cm.load_atlases()
    a = cm.atlas(acc, cfg)
    if a["status"] not in cm.FETCHABLE or not a.get("screens"):
        raise Refused(f"{acc}: no screens for this atlas")
    rdir = cm.outside_repo(Path(dest) / "raw" / acc)
    hp = rdir / "HEADS.json"
    if not hp.is_file():
        raise Refused(f"{acc}: no HEADS.json: run `fetch.py heads` first")
    heads = json.loads(hp.read_text())
    if heads.get("atlas_entry_sha256") != cm.entry_sha256(acc, cfg):
        raise Refused(f"{acc}: the heads were fetched under another atlas entry")
    if "matrix" in heads.get("roles", []):
        raise Refused(f"{acc}: the heads list a matrix: refused")
    checks = run_screens(a, cfg["screen_rules"], heads, rdir)
    out = {"accession": acc, "tasks": list(cm.TASKS), "heads_sha256": cm.sha256_file(hp),
           "plan_sha256": heads["plan_sha256"], "atlas_entry_sha256": heads["atlas_entry_sha256"],
           "samples": sorted(heads["samples"]), "checks": checks,
           "passed": all(c["passed"] for c in checks.values())}
    p = cm.dump_json(out, Path(dest) / "screens" / acc / "screen.json")
    out["screen_sha256"] = cm.sha256_file(p)
    return out


def paste_line(res: Mapping) -> str:
    failed = [k for k, c in res["checks"].items() if not c["passed"]]
    if not failed:
        return f'  {res["accession"]}: {{verdict: kept, screen_sha256: "{res["screen_sha256"]}", date: "YYYY-MM-DD"}}'
    return (f'  {res["accession"]}: {{verdict: dropped, screen_sha256: "{res["screen_sha256"]}", date: "YYYY-MM-DD", '
            f'failed: [{", ".join(failed)}]}}')


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
    summary = {k: ("yes" if c["passed"] else "no") for k, c in res["checks"].items()}
    print(json.dumps({"accession": res["accession"], "screens": summary, "checks": res["checks"],
                      "screen_sha256": res["screen_sha256"]}, indent=1))
    print("# paste under members: in experiments/t1_ext/screens.yaml (replace the date):")
    print(paste_line(res))
    return 0 if res["passed"] else 2


if __name__ == "__main__":
    sys.exit(main())
