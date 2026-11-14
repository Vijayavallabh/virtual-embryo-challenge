#!/usr/bin/env python
"""Row N48 of experiments/harness/TRIES.md (T1; a diagnostic, not a try): is there anything to shrink in M1's real
direction? Target-free.

Why. The ledger prunes "a noise-aware, subtype-resolved Feng direction" and "a James-Stein scale of Feng's step" with
the sentence "a split-half check of Feng's direction finds nothing to shrink". The seventh independent sweep of the
record (2026-11-14) called that closure a hard call: the check's result is not committed, and a counted route that
scores M1 arms on the E10.5 cuts exists (row T26). This script makes the check and commits it.

It reads Feng's clean statistics (X002, manifest-checked, through feng_dirs) and the gene names of the released
T1/E9.5_RNA.h5ad (names only, as N12 and N20: no matrix, no cell); no E10.5 file, target, score, board or built file,
and it counts, prints and stores no T1 cell.

Per direction group of M1 (vec.harness.LP_GROUPS: CM, Endocardium, Endothelium, Epicardium), with d = Feng's
E9.5 -> E14.5 log ratio of mean CP10k on the T1 gene order (M1's real direction), d_a and d_b the same from each half
of Feng's cells, and the decile basis (the mean CP10k over the two ends; feng_dirs.direction_halves):
  the ten decile factors of feng_dirs.eb_shrink, try T18's rule and constants unchanged, lowest expression first;
  the scalar James-Stein scale over the genes with every value, max(0, mean(d^2) - noise) / mean(d^2) with
  noise = mean(((d_a - d_b) / 2)^2): the one factor a James-Stein scale would put on the whole step;
  both again without the pre-registered held-flat genes, and the correlation of the two halves' directions (reported;
  they decide nothing).

Rule, fixed in the row before anything is read: try T18's stop (n20_proxy_profile.t18_stop) applied to the real
direction with its constants. The decile shrinkage is a no-op unless, in each of CM and Endocardium (the groups the
T1 backtests read), some factor of expression deciles 4-10 is below 0.9. The scalar scale is a no-op unless it is
below 0.9 in each of the two groups. A no-op closes that form by this stop. A form that is not a no-op is open: its
counted read (M1 with the shrunk direction against M1 on the two E10.5 cuts, T26's route) gets its own row before it
runs, and its result could only go to the captain, since no adoption-rule route exists for M1's direction.

  uv run --frozen --no-sync python -I -B experiments/t1_m1/n48_m1_direction_shrink.py

Writes experiments/t1_m1/results/n48_m1_direction_shrink/n48.{json,txt}, stamped (commit, clean-tree flag, UTC, the
inputs' sha256); refused unless the row is committed, the tree is clean and no output exists (a rerun is a new row).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in ("t1_hygiene", "t1_direction", "external"):
    sys.path.insert(0, str(REPO / "experiments" / _p))
sys.path.insert(0, str(HERE))
import dirlib as dl  # noqa: E402
import feng_dirs  # noqa: E402
import filterio as fio  # noqa: E402
import hygiene  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402

from vec import harness, provenance, t1_scope  # noqa: E402

ROW = "N48"
T1_FILE = REPO / "data" / "raw" / "T1" / "E9.5_RNA.h5ad"
OUT = HERE / "results" / "n48_m1_direction_shrink"
REAL = feng_dirs.BRACKET  # E9.5 -> E14.5, M1's direction
STOP_GROUPS = ("CM", "Endocardium")  # the groups the T1 backtests read (n20_proxy_profile.READ_GROUPS)
NO_OP = 0.9  # T18's stop: a factor at or above this barely changes the direction
EXPRESSED_FROM = 3  # index of decile 4 (factors are lowest expression first), as T18's stop


def js_scale(d, d_a, d_b) -> dict:
    """The scalar James-Stein scale of a direction d from its two half-sample versions: over the genes with all three
    values, noise = mean(((d_a - d_b) / 2)^2) (eb_shrink's estimate of d's noise variance), and the scale is
    max(0, mean(d^2) - noise) / mean(d^2); 1.0 where d is all zero, None where no gene has every value."""
    d, d_a, d_b = (np.asarray(x, dtype=np.float64) for x in (d, d_a, d_b))
    ok = np.isfinite(d) & np.isfinite(d_a) & np.isfinite(d_b)
    if not ok.any():
        return {"genes": 0, "noise": None, "mean_sq": None, "scale": None}
    noise = float(np.mean(((d_a[ok] - d_b[ok]) / 2.0) ** 2))
    mean_sq = float(np.mean(d[ok] ** 2))
    scale = max(mean_sq - noise, 0.0) / mean_sq if mean_sq > 0 else 1.0
    return {"genes": int(ok.sum()), "noise": noise, "mean_sq": mean_sq, "scale": float(scale)}


def half_correlation(d_a, d_b) -> float | None:
    """Pearson correlation of the two halves' directions over the genes both have (None under three genes or where a
    half is constant)."""
    a, b = np.asarray(d_a, dtype=np.float64), np.asarray(d_b, dtype=np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return None
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def group_shrink(d, halves: dict, held) -> dict:
    """One group's record: the decile factors and the scalar scale, with and without the held-flat genes."""
    d = np.asarray(d, dtype=np.float64)
    held = np.asarray(held, dtype=bool)
    d_a, d_b, basis = halves["d_a"], halves["d_b"], halves["mean_cp"]
    _, factors = feng_dirs.eb_shrink(d, d_a, d_b, basis, harness.LP_BINS)
    free = np.where(held, np.nan, d)
    _, factors_free = feng_dirs.eb_shrink(free, d_a, d_b, basis, harness.LP_BINS)
    return {"genes_with_value": int(np.isfinite(d).sum()), "eb_factors": [float(f) for f in factors],
            "js": js_scale(d, d_a, d_b),
            "eb_factors_without_held": [float(f) for f in factors_free], "js_without_held": js_scale(free, d_a, d_b),
            "half_correlation": half_correlation(d_a, d_b)}


def stop(groups: dict) -> dict:
    """The row's rule. A form is open only if it would change the direction in each of STOP_GROUPS: the decile form if
    some factor of deciles 4-10 is below NO_OP there, the scalar form if the scale is. A group without a record counts
    as no."""
    eb = {g: (bool(min(groups[g]["eb_factors"][EXPRESSED_FROM:]) < NO_OP) if g in groups else None)
          for g in STOP_GROUPS}
    js = {g: (None if g not in groups or groups[g]["js"]["scale"] is None else bool(groups[g]["js"]["scale"] < NO_OP))
          for g in STOP_GROUPS}
    return {"no_op_at": NO_OP, "deciles": f"{EXPRESSED_FROM + 1}-{harness.LP_BINS}", "eb_below_no_op": eb,
            "js_below_no_op": js, "eb_open": all(v is True for v in eb.values()),
            "js_open": all(v is True for v in js.values())}


def text(report: dict) -> str:
    s = report["stop"]
    lines = ["N48: the reliability shrinkage of M1's real direction (Feng E9.5 -> E14.5), target-free",
             f"decile form (T18's rule): {'OPEN: it would change the direction' if s['eb_open'] else 'no-op: closed'} "
             f"(some factor of deciles {s['deciles']} below {s['no_op_at']}: {s['eb_below_no_op']})",
             f"scalar James-Stein scale: {'OPEN: it would change the step' if s['js_open'] else 'no-op: closed'} "
             f"(scale below {s['no_op_at']}: {s['js_below_no_op']})",
             "  group: genes; scalar scale (without held genes); halves' correlation; decile factors, lowest "
             "expression first"]
    for g, r in report["groups"].items():
        js, jf, hc = r["js"]["scale"], r["js_without_held"]["scale"], r["half_correlation"]
        lines.append(f"  {g}: {r['genes_with_value']}; {'none' if js is None else f'{js:.4f}'} "
                     f"({'none' if jf is None else f'{jf:.4f}'}); {'none' if hc is None else f'{hc:.3f}'}; "
                     + " ".join(f"{f:.3f}" for f in r["eb_factors"]))
    lines += [f"Cells per group at the two ends (Feng): {report['real_cells']}",
              f"Commit {report['stamp']['commit']}, clean at end {report['stamp']['tree_clean_at_end']}, "
              f"UTC {report['stamp']['utc']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--t1-file", default=str(T1_FILE))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if ROW not in n12.tries_rows():
        sys.exit(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if not feng_dirs.x002_amended():
        sys.exit("EXTERNAL_DATA.md's X002 row does not carry its amendment: a use of an external dataset is logged "
                 "before it enters anything")
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{ROW} has already written its outputs; a rerun is a new row (experiments/harness/TRIES.md)")
    start = provenance.git_state(provenance.REPO_ROOT)
    if start[1] is None or start[1]:
        sys.exit(f"the tree is not clean ({start[1]}): {ROW} runs once, on a committed tree")
    cfg = feng_dirs.load_config()
    genes = n12.t1_gene_names(Path(a.t1_file))
    path = feng_dirs.clean_stats_path()
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    dirs, cells, idx = feng_dirs.directions_from_stats(st, genes, cfg, REAL)
    halves = feng_dirs.direction_halves(st, genes, cfg, REAL)
    held_set = set(t1_scope.held_flat_genes())
    held = np.array([g in held_set for g in genes], dtype=bool)
    groups = {g: group_shrink(dirs[g], halves[g], held) for g in harness.LP_GROUPS if g in dirs and g in halves}
    report = {"label": "N48: target-free reliability shrinkage of M1's real direction (no E10.5 file, no target, no "
                       "score, no board, no built file, no T1 cell count)",
              "bracket": list(REAL), "stop": stop(groups), "groups": groups,
              "real_cells": {g: cells[g] for g in sorted(cells)},
              "n_t1_genes": len(genes), "mapped_genes": int((idx >= 0).sum()), "held_flat_genes": int(held.sum()),
              "inputs": {"t1_file": str(Path(a.t1_file).resolve().relative_to(REPO)),
                         "t1_file_sha256": hygiene.sha256_file(a.t1_file),
                         "feng_stats": str(path.relative_to(REPO)), "feng_stats_sha256": entry["sha256"],
                         "held_flat_sha256": hygiene.sha256_file(t1_scope.HELD_FLAT_GENES)},
              "stamp": {"commit": start[0], "tree_clean_at_start": True, "tree_clean_at_end": None,
                        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    end_commit, end_unclean = provenance.git_state(provenance.REPO_ROOT)
    report["stamp"]["tree_clean_at_end"] = None if end_unclean is None else (not end_unclean and end_commit == start[0])
    out.mkdir(parents=True, exist_ok=True)
    (out / "n48.json").write_text(json.dumps(report, indent=1) + "\n")
    (out / "n48.txt").write_text(text(report))
    print(text(report), end="")


if __name__ == "__main__":
    main()
