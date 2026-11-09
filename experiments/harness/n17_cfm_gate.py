#!/usr/bin/env python
"""N17 (experiments/harness/TRIES.md; a diagnostic, not a try): the target-free gate of heart try H6, candidate b0_cfm.

H6 filters each lineage's momentum by the lineage's own within-stage covariance (vec.t2_build.covariance_filter,
CFM_RULE): the momentum is projected onto the directions of the depth-residualised, gene-standardised log1p CP10k of the
lineage's reference-stage cells whose eigenvalues clear the Marchenko-Pastur edge, then renormalised to its unfiltered
L2 norm.

For each momentum pair B0 uses, this script runs that filter on the candidate's own code path (t1_scope.momentum_shift
with logf_filter=covariance_filter, on the harness's types: harmonised `merged` lineages, Untracked held flat) and
records per lineage the number of directions above the edge, the retained energy and the Spearman correlation between
filtered and unfiltered factors. The pairs: the test pair (E9.5 -> E10.5, on the test board's panel; E10.5 is unlocked,
A20) and the two counting backtests' momentum pairs (their start and reference stages, on their panel).

The gate, fixed in row N17 before this runs: a pair passes if at least one lineage with at least 200 reference cells has
a direction above the edge (otherwise the filter is a no-op there and the try would test nothing), and every such
lineage that has a direction has its Spearman inside [0.30, 0.95] (below: the momentum destroyed; above: a no-op). A
lineage without a direction keeps B0's factor and is reported; a lineage under 200 cells, or one B0 does not move
(Untracked; fewer than 20 cells at a stage), is reported and not judged. The gate passes only if all three pairs pass;
then H6's row is committed and goes to the captain with its plain summary before the counted run. If it fails, H6 is
pruned by its gate, not run and not counted (as E10).

Target-free: it reads the reference and start stages of the pairs (E10.5 only as the test pair's reference, A20) and no
score, board limit or built file. It writes and prints no cell count: no cell numbers, no eigenvalues, no
Marchenko-Pastur edge (from which a lineage's cell number could be solved), only directions, shares, correlations and a
boolean per lineage.

Usage (DGX, repo root, BLAS threads capped):
    OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
        uv run --frozen --no-sync python -I -B experiments/harness/n17_cfm_gate.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import textwrap
import time
from pathlib import Path

import numpy as np

from vec import harness, provenance, t1_scope, t2_build
from vec.boards import load_boards
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "n17_cfm_gate"
BOARD = "T2:heart:test_extrap"
TEST_PAIR = ("test/E9.5->E10.5", "T2/E9.5.h5ad", "T2/E10.5.h5ad")  # build_test.py heart: start, reference; A20
BAND = (0.30, 0.95)  # the row's Spearman band for every lineage with at least 200 reference cells and a direction
DATA = harness.DATA_DIR
OUT_NAME = "n17_cfm_gate"


def momentum_pairs() -> list[dict]:
    """The momentum pairs (start, reference, panel board) of the test build and the board's counting backtests, from the
    pseudo-task config; the counting backtests' own panel is their panel board's."""
    tasks = harness.load_tasks()[BOARD]
    pairs = [{"key": TEST_PAIR[0], "start": TEST_PAIR[1], "ref": TEST_PAIR[2], "panel_board": BOARD}]
    pairs += [{"key": t.key, "start": t.inputs["start"], "ref": t.reference, "panel_board": t.panel_board}
              for t in tasks if t.counted]
    return pairs


def typed(a, tag: str):
    """vec.harness._b0_t2's types: the harmonised lineage (obs["merged"]), Untracked held flat by a stage-tagged one."""
    m = np.asarray(a.obs["merged"]).astype(str)
    return harness._shallow(a, b0_type=np.where(m == "Untracked", f"Untracked@{tag}", m))


def lineage_rows(start, ref) -> list[dict]:
    """Per lineage of the reference stage: B0's momentum status, and the filter's record for lineages with momentum. The
    candidate's own code path: momentum_shift with its filter (the applied cells are discarded); the filter uses no
    random numbers, so these values do not depend on the seed."""
    _, rep = t1_scope.momentum_shift(typed(start, "start"), typed(ref, "end"), typed(ref, "end"), harness.B0_DAMPING,
                                     (), key="b0_type", seed=[harness.BUILD_SEED_OFFSET, 3],
                                     logf_filter=t2_build.covariance_filter)
    rows = []
    for r in rep.rows:
        rec = r.get("filter")
        rows.append({"lineage": r["type"], "momentum": r["status"], "eligible": bool(rec and rec["eligible"]),
                     "filter": None if rec is None else rec["status"],
                     "directions": 0 if rec is None else rec["directions"],
                     "retained_energy": None if rec is None else rec["retained_energy"],
                     "spearman": None if rec is None else rec["spearman"],
                     "max_abs_ratio": None if rec is None else rec["max_abs_ratio"]})
    return rows


def judge_pair(rows: list[dict], band: tuple[float, float] = BAND) -> dict:
    """The row's rule on one pair: at least one eligible lineage (200 reference cells or more) with a direction, and
    every eligible lineage with a direction inside the band (a Spearman that is not finite counts as outside)."""
    low, high = band
    eligible = [r for r in rows if r["eligible"]]
    acting = [r for r in eligible if r["directions"] >= 1]
    outside = [r["lineage"] for r in acting if r["spearman"] is None or not low <= r["spearman"] <= high]
    return {"eligible": [r["lineage"] for r in eligible], "acting": [r["lineage"] for r in acting],
            "outside_band": outside, "passes": bool(acting) and not outside}


def verdict(pairs: dict[str, dict]) -> str:
    bad = [k for k, p in pairs.items() if not p["judged"]["passes"]]
    if not bad:
        return "pass: H6 may be registered and run (after its row and the captain's sign-off)"
    return "prune H6 by its gate, not run and not counted: fails on " + "; ".join(bad)


def _file_sha256(path: Path) -> str:
    return harness._sha256(path)


def _panel_sha256(genes) -> str:
    return hashlib.sha256("\n".join(map(str, genes)).encode()).hexdigest()


def report_text(res: dict) -> str:
    lines = [f"N17: target-free gate of heart try H6 `b0_cfm` (reads reference and start stages only; no score, no "
             f"cell count). Band {BAND[0]:.2f} to {BAND[1]:.2f}, lineages with at least {t2_build.CFM_MIN_CELLS} "
             "reference cells and a direction above the edge. Directions are eigenvalues above (1 + sqrt(p / n))^2."]
    for key, p in res["pairs"].items():
        lines += ["", f"{key} (start {p['start']}, reference {p['ref']}, panel {p['panel_board']}): "
                      f"{'passes' if p['judged']['passes'] else 'does NOT pass'}"]
        lines.append(f"  eligible lineages: {', '.join(p['judged']['eligible']) or 'none'}; acting: "
                     f"{', '.join(p['judged']['acting']) or 'none'}; outside the band: "
                     f"{', '.join(p['judged']['outside_band']) or 'none'}")
        for r in p["lineages"]:
            if r["filter"] is None:
                what = r["momentum"]
            elif r["directions"] == 0:
                what = r["filter"]
            else:
                what = (f"{r['directions']} direction(s), retained energy {r['retained_energy']:.3f}, Spearman "
                        + ("n/a" if r["spearman"] is None else f"{r['spearman']:.3f}")
                        + f", max abs factor ratio {r['max_abs_ratio']:.2f}")
            lines.append(f"    {r['lineage']:<22s} {what}")
    s = res["stamp"]
    lines += ["", f"Inputs (sha256): {json.dumps(res['inputs']['files'], sort_keys=True)}",
              f"Panels (sha256 of the gene list): {json.dumps(res['inputs']['panels'], sort_keys=True)}",
              f"Config: {res['inputs']['config_sha256']}",
              f"Commit {s['commit']}, clean tree {s['tree_clean']}, UTC {s['utc']}",
              "", f"verdict: {res['verdict']}"]
    return "\n".join(textwrap.fill(ln, 118, subsequent_indent="      ") if len(ln) > 118 else ln for ln in lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="N17: the target-free gate of heart try H6")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N17 has already written its outputs; a rerun needs a new TRIES row")
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    boards = load_boards()
    pairs = momentum_pairs()
    panels = {p["panel_board"]: list(boards[p["panel_board"]].genes) for p in pairs}  # the panels, never the limits
    res = {"what": "row N17 of experiments/harness/TRIES.md (a diagnostic, not a try): the target-free gate of H6",
           "rule": {"band": list(BAND), "lineage_min_cells": t2_build.CFM_MIN_CELLS, "filter": t2_build.CFM_RULE},
           "pairs": {}}
    files = {}
    for p in pairs:
        genes = panels[p["panel_board"]]
        start = load_stage(DATA / p["start"], genes=genes, drop_duplicates=True)
        ref = load_stage(DATA / p["ref"], genes=genes, drop_duplicates=True)
        for f in (p["start"], p["ref"]):
            files[f] = _file_sha256(DATA / f)
        rows = lineage_rows(start, ref)
        res["pairs"][p["key"]] = {"start": p["start"], "ref": p["ref"], "panel_board": p["panel_board"],
                                  "lineages": rows, "judged": judge_pair(rows)}
        print(f"{p['key']}: {'passes' if res['pairs'][p['key']]['judged']['passes'] else 'does not pass'}", flush=True)
    res["verdict"] = verdict(res["pairs"])
    res["inputs"] = {"files": files, "panels": {b: _panel_sha256(g) for b, g in panels.items()},
                     "config_sha256": _file_sha256(harness.CONFIG)}
    res["stamp"] = {"commit": commit, "tree_clean": unclean == [],
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{OUT_NAME}.json").write_text(json.dumps(res, indent=1) + "\n")
    text = report_text(res)
    (out / f"{OUT_NAME}.txt").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
