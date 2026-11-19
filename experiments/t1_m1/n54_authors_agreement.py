"""FOR THE METHOD REPORT ONLY. Row N54's one read of the authors' labels: per direction group, how far the recomputed
groups of Feng's permitted cells (experiments/external/feng_relabel.py) agree with the authors' cell_type_spec.

This is the only code in the repository that reads cell_type_spec next to the recomputed labels, and it is kept apart
on purpose (feng_relabel.yaml, `authors_labels`): agreement with the authors' labels selects nothing. No gene,
threshold, margin, set or fallback is chosen, edited or kept by what it shows; nothing imports this module or reads
its output (tests/test_feng_relabel.py checks the labeller, feng_dirs, the builder, the harness and the listing).

It runs once, only after the captain's ruling on T1's last upload (recomputed file or none file) is committed, so it
cannot inform that ruling (the final review of the port, lesser point 1): it refuses unless --upload-ruling names a
docs/PLAN.md decision-log entry naming N54 and containing the fixed phrase "rules on T1's last upload"
(feng_relabel.UPLOAD_RULING_PHRASE) that is not the relabelling's own ruling, is dated on or after the listing's
commit and was not yet in docs/PLAN.md at that commit, read there (feng_relabel.later_ruling); the listing's n54.json
exists, records the same recomputed sha256 as --frozen-sha256 and is committed at HEAD as written; the file and its
manifest carry that sha256; the tree is clean; and its own output does not exist.
It writes aggregates only, per stage and group: the cells under each labelling, those in both, the share of the
authors' group recovered, the share of the recomputed group the authors placed elsewhere, and the Jaccard index. No
per-cell, per-gene or per-cluster detail. Output: experiments/t1_m1/results/n54_authors_agreement/agreement.json.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in ("t1_direction", "external"):
    sys.path.insert(0, str(REPO / "experiments" / _p))
sys.path.insert(0, str(HERE))
import feng_dirs  # noqa: E402
import filterio as fio  # noqa: E402

AUTHORS_COLUMN = "cell_type_spec"  # read here and nowhere else
LISTING = HERE / "results" / "n54_relabel_listing" / "n54.json"
OUT = HERE / "results" / "n54_authors_agreement"


def agreement(authors: Mapping[str, np.ndarray], recomputed: Mapping[str, np.ndarray]) -> dict:
    """Per group: masks of the same cells under the authors' and the recomputed labelling -> aggregate agreement."""
    out = {}
    for g in recomputed:
        a, r = np.asarray(authors[g], dtype=bool), np.asarray(recomputed[g], dtype=bool)
        both, union = int((a & r).sum()), int((a | r).sum())
        out[g] = {"authors": int(a.sum()), "recomputed": int(r.sum()), "both": both,
                  "recovered_share": both / int(a.sum()) if a.any() else None,
                  "elsewhere_share": 1.0 - both / int(r.sum()) if r.any() else None,
                  "jaccard": both / union if union else None}
    return out


def group_masks(labels: np.ndarray, groups: Mapping[str, Sequence[str]]) -> dict[str, np.ndarray]:
    return {g: np.isin(labels, list(labs)) for g, labs in groups.items()}


def read_authors(obs_path: Path, cells: Sequence[str]) -> np.ndarray:
    """The authors' label of each of `cells`, by name, from the clean reduction's obs table."""
    with open(obs_path, newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        j = header.index(AUTHORS_COLUMN)
        lab = {row[0]: row[j] for row in rd}
    return np.array([lab[c] for c in cells], dtype=str)


def main(argv=None) -> None:
    from vec import provenance

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--frozen-sha256", required=True, help="the recomputed statistics file's sha256, as N54 read it")
    ap.add_argument("--method-report", action="store_true", required=True,
                    help="acknowledges that this output is for the method report only and decides nothing")
    ap.add_argument("--upload-ruling", required=True,
                    help="the captain's docs/PLAN.md decision-log entry on T1's last upload (naming N54 and "
                         "containing \"rules on T1's last upload\"), made after the listing's commit")
    a = ap.parse_args(argv)
    if (OUT / "agreement.json").exists():
        sys.exit("the agreement has been computed once already; it is not recomputed")
    if not LISTING.exists():
        sys.exit("N54's listing has not run: the agreement comes after it")
    import feng_relabel

    if feng_relabel.committed_at_head(LISTING) != LISTING.read_bytes():
        sys.exit("N54's listing is not committed at HEAD as written: the agreement comes after the listing's commit")
    raw = feng_relabel.committed_at_head(feng_relabel.RECORD)
    run_ruling = ((json.loads(raw).get("gate") or {}).get("ruling")) if raw is not None else None
    why = feng_relabel.later_ruling(a.upload_ruling, LISTING, feng_relabel.ROW, excluded=[run_ruling],
                                    phrase=feng_relabel.UPLOAD_RULING_PHRASE)
    if why is not None:
        sys.exit(f"the agreement comes only after the captain's ruling on T1's last upload: {why}")
    listed = json.loads(LISTING.read_text())["inputs"]["recomputed"]["sha256"]
    rcfg = feng_dirs.load_relabel_config()
    rpath = feng_dirs.recomputed_stats_path(rcfg)
    entry = fio.verify_manifest_output(rpath)
    if not (entry["sha256"] == listed == a.frozen_sha256):
        sys.exit("the recomputed file's sha256, the listing's and --frozen-sha256 differ: the file is not the one "
                 "N54 read")
    start = provenance.git_state(provenance.REPO_ROOT)
    if start[1] is None or start[1]:
        sys.exit(f"the tree is not clean ({start[1]})")
    with open(rpath.parent / rcfg["output"]["cells"], newline="") as f:
        rd = csv.reader(f)
        header = next(rd)
        rows = [(r[0], r[header.index("stage")], r[header.index("label")]) for r in rd]
    cells = [r[0] for r in rows]
    stage = np.array([r[1] for r in rows], dtype=str)
    rlab = np.array([r[2] for r in rows], dtype=str)
    alab = read_authors(REPO / rcfg["inputs"]["dir"] / rcfg["inputs"]["obs"], cells)
    cfg = feng_dirs.load_config()
    authors_groups = {g: list(cfg["external_labels"][g][0] or []) for g in rcfg["groups"]}
    report = {"label": "FOR THE METHOD REPORT ONLY: agreement of the recomputed groups with the authors' labels; it "
                       "selects nothing",
              "upload_ruling": a.upload_ruling, "recomputed_sha256": entry["sha256"],
              "relabel_config_sha256": fio.sha256_file(feng_dirs.RELABEL_CONFIG),
              "authors_groups": authors_groups, "recomputed_groups": rcfg["groups"],
              "by_stage": {s: agreement(group_masks(alab[stage == s], authors_groups),
                                        group_masks(rlab[stage == s], rcfg["groups"])) for s in rcfg["stages"]},
              "stamp": {"commit": start[0], "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "agreement.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report["by_stage"], indent=1))


if __name__ == "__main__":
    main()
