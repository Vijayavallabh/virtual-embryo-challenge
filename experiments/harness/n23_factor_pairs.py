#!/usr/bin/env python
"""Diagnostic N23 of experiments/harness/TRIES.md (target-free, not a try): is the per-gene factor of the series
rendering (tries E24a to E24c) a property of the run, or development?

A cold reviewer's check of 2026-11-12 on E24's recommendation. For pairs of RELEASED embryo stages it computes the
factor exactly as the candidate does (vec.t2_build.series_labels and series_factor: the matched extraembryonic labels,
eps, rescale and clip of SERIES_RULE) and reports its size and its rank correlation with the test pair's factor
(E7.5 -> E8.0). A pair across the depth step between E7.5 and E8.0 holds a change of run plus development; a pair of
neighbouring stages on one side of it holds development, plus whatever differs between their two samples.

It reads released stages only: no target stage of any board, no score, no candidate build.

Usage (repo root, BLAS threads capped as docs/COMPUTE.md says):
  uv run --frozen --no-sync python -I -B experiments/harness/n23_factor_pairs.py \
    --out experiments/harness/results/n23_factor_pairs
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from vec import harness, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.prep import load_stage

DATA = REPO_ROOT / "data" / "raw" / "T2"
BOARD = "T2:embryo:test_interp"
DAY = {"E6.75": 6.75, "E7.25": 7.25, "E7.5": 7.5, "E8.0": 8.0, "E8.25_late": 8.3, "E8.5": 8.5, "E8.75": 8.75,
       "E9.5": 9.5, "E10.5": 10.5}
STEP = 7.75  # the depth step lies between E7.5 and E8.0 (docs/DATA_NOTES.md)
TEST_PAIR = ("E7.5", "E8.0")
PAIRS = (TEST_PAIR, ("E7.5", "E8.25_late"), ("E7.5", "E8.5"), ("E7.5", "E8.75"), ("E7.25", "E8.0"),
         ("E6.75", "E7.25"), ("E7.25", "E7.5"), ("E6.75", "E7.5"),
         ("E8.0", "E8.25_late"), ("E8.25_late", "E8.5"), ("E8.5", "E8.75"), ("E8.0", "E8.5"), ("E8.0", "E8.75"),
         ("E8.75", "E9.5"), ("E9.5", "E10.5"))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    commit = harness._git("rev-parse", "HEAD")
    clean = harness._git("status", "--porcelain", "--untracked-files=no") == ""
    genes = list(load_boards()[BOARD].genes)
    stages = {}
    for s in DAY:
        st = load_stage(DATA / f"{s}.h5ad", genes=genes, drop_duplicates=True)
        _, C, lib, _ = t2_build._stage_counts(st, genes)
        stages[s] = {"C": C, "names": np.asarray(st.obs["celltype"]).astype(str),
                     "groups": np.asarray(st.obs["group"]).astype(str), "depth": float(np.median(lib))}
    logf, rec = {}, {}
    for frm, to in PAIRS:
        labs = t2_build.series_labels(stages[frm]["names"], stages[frm]["groups"], stages[to]["names"])
        key = f"{frm}->{to}"
        kind = "across the step" if (DAY[frm] < STEP) != (DAY[to] < STEP) else "one side"
        if not labs:
            rec[key] = {"kind": kind, "days": DAY[to] - DAY[frm], "labels": [], "note": "no matched label"}
            continue
        phi = t2_build.series_factor(stages[frm]["C"], stages[frm]["names"], stages[to]["C"], stages[to]["names"], labs)
        logf[(frm, to)] = np.log(phi)
        lo, hi = t2_build.SERIES_CLIP
        rec[key] = {"kind": kind, "days": DAY[to] - DAY[frm], "labels": labs,
                    "sd_log_factor": float(np.log(phi).std()),
                    "genes_beyond_2x": int((np.abs(np.log(phi)) > np.log(2.0)).sum()),
                    "genes_at_the_clip": int(((phi <= lo * 1.0001) | (phi >= hi * 0.9999)).sum()),
                    "depth_ratio": stages[to]["depth"] / stages[frm]["depth"]}
    ref = logf[TEST_PAIR]
    for (frm, to), v in logf.items():
        if (frm, to) != TEST_PAIR:
            rec[f"{frm}->{to}"]["spearman_with_test_pair"] = float(spearmanr(ref, v).statistic)
    res = {"what": __doc__.split("\n")[0], "commit": commit, "tree_clean": clean, "board_panel": BOARD,
           "n_genes": len(genes), "rule": t2_build.SERIES_RULE, "median_counts_per_cell":
               {s: stages[s]["depth"] for s in DAY}, "pairs": rec}
    lines = [f"N23 (commit {commit}, clean {clean}; target-free): the series factor between released stages, "
             f"{len(genes)} panel genes",
             "  pair                      kind             days  labels  SD log factor  beyond 2x  at clip  depth ratio"
             "  Spearman with E7.5->E8.0"]
    for key, r in rec.items():
        if not r["labels"]:
            lines.append(f"  {key:<25} {r['kind']:<16} {r['days']:.2f}  none")
            continue
        rho = r.get("spearman_with_test_pair")
        lines.append(f"  {key:<25} {r['kind']:<16} {r['days']:.2f}  {len(r['labels']):6d}  {r['sd_log_factor']:13.3f}  "
                     f"{r['genes_beyond_2x']:9d}  {r['genes_at_the_clip']:7d}  {r['depth_ratio']:11.2f}  "
                     + ("(itself)" if rho is None else f"{rho:+.2f}"))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "factor_pairs.json").write_text(json.dumps(res, indent=1))
    (out / "factor_pairs.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
