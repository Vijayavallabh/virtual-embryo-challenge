#!/usr/bin/env python
"""Rows N20 and N20b of experiments/harness/TRIES.md (T1; a diagnostic, not a try): the late-proxy direction's profile,
read
before tries T13-T18 run. Target-free: it reads Feng's clean statistics (X002, manifest-checked, through feng_dirs), the
committed extracts of X007 and X008 and the gene names of the released T1/E9.5_RNA.h5ad (names only, as N12: no
matrix, no cell); no E10.5 file, target, score, board or built file, and it counts, prints and stores no T1 cell.

Per direction group of M1 (vec.harness.LP_GROUPS: CM, Endocardium, Endothelium, Epicardium), with log f = 2/7 x Feng's
E14.5 -> E18 log ratio of mean CP10k (a one-day step's share of the 3.5-day bracket; the T1 gene order; the
pre-registered held-flat list at 0) and d_real = Feng's E9.5 -> E14.5 log ratio (M1's real direction):
  Feng's cells at each end (E14.5, E18; E9.5, E14.5) and the T1 genes with a value;
  the proxy's size: the mean abs(log f) and the genes it moves by more than 10% (abs(log f) > log 1.1);
  the proxy's similarity to the real direction, over the genes with a value in both: Spearman correlation, cosine,
  the share of the proxy's 400 largest abs(d_proxy) genes whose d_real has the same sign, and the overlap of the 400
  largest abs(d) of each;
  in the proxy's top 400 up and down by log f: how many genes of each of R2's modules (D, A1, N1, N2), as N12 does for
  M1's real factor;
  the reliability shrinkage (feng_dirs.eb_shrink) per expression decile: its factors, one per decile, which every gene
  of the decile shares (it cannot reorder a decile).

Its gate (N20's row, fixed before anything is read): the tries T13-T18 run only if Feng's statistics hold a stage
E14.5 and a stage E18, the direction group of each within-type read (READ_GROUPS: CM, and Endocardium for the
endocardium read) has at least min_cells_ext (10) cells at both ends of the proxy bracket, and at least 20,000 T1 genes
have a Feng value. Otherwise the tries are not run and do not count.

Two more numbers, fixed in the row before anything is read, that late_proxy.py reads and that stop none of the gate's
tries:
  the landing floor: CM's Spearman correlation of the proxy with the real direction, each from the other half of Feng's
  cells (proxy half A against real half B, and the reverse, their mean; over the genes with all four values), must be
  at least 0.2. Below it, a pass of T16, T17 or T18 is platform evidence only, with no recommendation for M1. The
  directions share Feng's E14.5 with opposite signs, so the plain correlation of the pooled directions is pulled down
  by E14.5's sampling noise; disjoint halves remove that. The plain correlation is reported and decides nothing.
  T18's stop: T18 does not run (and does not count) unless, in each of CM and Endocardium, some factor of expression
  deciles 4-10 is below 0.9; otherwise its direction barely changes on those genes, and the variants' bar, which asks
  for a gain on reads of both groups, could not be met.

Usage (DGX, repo root; refused unless the row is committed, X002's amendment is in EXTERNAL_DATA.md and the tree is
clean):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/t1_m1/n20_proxy_profile.py

Writes experiments/t1_m1/results/n20b_proxy_profile/n20_proxy_profile.{json,txt}, stamped (commit, vec.provenance's
clean-tree flag, UTC, the inputs' sha256); an existing output is refused (a rerun is a new row). The tries refuse to
run on any commit but N20's, so nothing is committed between N20 and the tournament's last decision.
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

T1_FILE = REPO / "data" / "raw" / "T1" / "E9.5_RNA.h5ad"
OUT = HERE / "results" / "n20b_proxy_profile"
PROXY = harness.LP_BRACKET
REAL = feng_dirs.BRACKET
FRACTION = 2 / 7  # a one-day step's share of the 3.5-day proxy bracket
TOP = 400
MIN_MAPPED = 20000
GATE_GROUP = "CM"
READ_GROUPS = {"T1/E8.5->E9.5:cm": "CM", "T1/E9.5->E10.5:cm": "CM",
               "T1/E9.5->E10.5:endocardium": "Endocardium"}  # late_proxy.WITHIN: each read's direction group
LANDING_FLOOR = 0.2  # CM's cross-half Spearman of the proxy with the real direction; below it no recommendation for M1
EB_NO_OP = 0.9  # T18's stop: a decile factor at or above this barely changes the direction
EB_EXPRESSED_FROM = 3  # index of decile 4 (factors are lowest expression first): low deciles' changes are already small


def _rel(path: Path) -> str:
    return str(Path(path).resolve().relative_to(REPO)) if Path(path).resolve().is_relative_to(REPO) else str(path)


def cosine(a: np.ndarray, b: np.ndarray) -> float | None:
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    return None if na == 0 or nb == 0 else float(a @ b / (na * nb))


def group_profile(group: str, d_proxy: np.ndarray, d_real: np.ndarray, halves: dict, genes: list[str],
                  held: np.ndarray, modules: dict[str, set[str]]) -> dict:
    """The row's numbers for one group (proxy log f = FRACTION x d_proxy, the held-flat genes at 0, NaN as no value)."""
    logf = n12.m1_log_factor(d_proxy, FRACTION, held)
    both = np.isfinite(d_proxy) & np.isfinite(d_real)
    top_proxy = np.argsort(-np.abs(np.nan_to_num(d_proxy)), kind="stable")[:TOP]
    top_real = np.argsort(-np.abs(np.nan_to_num(d_real)), kind="stable")[:TOP]
    in_both = top_proxy[both[top_proxy]]
    agree = float(np.mean(np.sign(d_proxy[in_both]) == np.sign(d_real[in_both]))) if in_both.size else None
    rank = n12.rankings(logf, genes)
    lists = {}
    for side in ("up", "down"):
        idx = rank[side][:TOP]
        names = {genes[i] for i in idx if logf[i] != 0}
        lists[f"top_{TOP}_{side}"] = {m: len(names & syms) for m, syms in sorted(modules.items())}
    eb, factors = feng_dirs.eb_shrink(d_proxy, halves["d_a"], halves["d_b"], halves["mean_cp"], harness.LP_BINS)
    return {"genes_with_value": int(np.isfinite(d_proxy).sum()),
            "proxy_mean_abs_log_f": float(np.abs(logf).mean()), "proxy_genes_moved_gt_10pct":
            int((np.abs(logf) > np.log(1.1)).sum()),
            "spearman_with_real": n12.spearman(d_proxy[both], d_real[both]),
            "cosine_with_real": cosine(d_proxy[both], d_real[both]),
            "top400_sign_agreement_with_real": agree,
            "top400_abs_overlap_with_real": len(set(top_proxy.tolist()) & set(top_real.tolist())),
            "modules": lists, "eb_factors": [round(f, 4) for f in factors],
            "eb_spearman": n12.spearman(np.nan_to_num(d_proxy), np.nan_to_num(eb))}


def landing_floor(st: dict, bracket: dict | None, genes: list[str], cfg: dict) -> dict:
    """The landing floor: CM's Spearman of the proxy with the real direction, each from the other half of Feng's cells
    (proxy half A against real half B, and the reverse, the mean of the two), over the genes with all four values, so
    that the sampling noise of E14.5, which the two directions share with opposite signs, cannot pull it down."""
    record = {"group": GATE_GROUP, "spearman_cross_half": None, "floor": LANDING_FLOOR, "met": False}
    if bracket is None:
        return record
    proxy = feng_dirs.direction_halves(bracket, genes, cfg, PROXY, groups={GATE_GROUP}).get(GATE_GROUP)
    real = feng_dirs.direction_halves(st, genes, cfg, REAL, groups={GATE_GROUP}).get(GATE_GROUP)
    if proxy is None or real is None:
        return record
    ok = np.isfinite(proxy["d_a"]) & np.isfinite(proxy["d_b"]) & np.isfinite(real["d_a"]) & np.isfinite(real["d_b"])
    ab, ba = n12.spearman(proxy["d_a"][ok], real["d_b"][ok]), n12.spearman(proxy["d_b"][ok], real["d_a"][ok])
    value = None if ab is None or ba is None else (ab + ba) / 2.0
    record.update(spearman_cross_half=value, met=value is not None and value >= LANDING_FLOOR)
    return record


def t18_stop(groups: dict) -> dict:
    """T18's target-free stop: in each group of the within-type reads, does any factor of expression deciles 4-10 fall
    below EB_NO_OP? T18 runs only if so in every group (a group with no factors counts as no)."""
    factors = {g: groups[g]["eb_factors"] if g in groups else None for g in sorted(set(READ_GROUPS.values()))}
    below = {g: None if f is None else bool(min(f[EB_EXPRESSED_FROM:]) < EB_NO_OP) for g, f in factors.items()}
    return {"no_op_at": EB_NO_OP, "deciles": f"{EB_EXPRESSED_FROM + 1}-{harness.LP_BINS}", "below_no_op": below,
            "runs": all(v is True for v in below.values())}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--t1-file", default=str(T1_FILE))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if "N20b" not in n12.tries_rows():
        sys.exit("N20b has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if not feng_dirs.x002_amended():
        sys.exit(f"EXTERNAL_DATA.md's X002 row does not carry the amendment (its phrase {feng_dirs.X002_AMENDMENT!r}): "
                 "a use of an external dataset is logged before it enters anything")
    if out.exists() and any(out.iterdir()):
        sys.exit("N20 has already written its outputs; a rerun is a new row (experiments/harness/TRIES.md)")
    start = provenance.git_state(provenance.REPO_ROOT)
    if start[1] is None or start[1]:
        sys.exit(f"the tree is not clean ({start[1]}): N20 runs once, on the commit the tries will run at")
    cfg = feng_dirs.load_config()
    genes = n12.t1_gene_names(Path(a.t1_file))
    path = feng_dirs.clean_stats_path()
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    gate = {"stages_present": all(s in {k.split("|")[0] for k in st["keys"]} for s in PROXY)}
    bracket = feng_dirs.bracket_stats(st, PROXY) if gate["stages_present"] else None
    proxy_dirs, proxy_cells, idx = (feng_dirs.directions_from_stats(bracket, genes, cfg, PROXY) if bracket is not None
                                    else ({}, {}, np.full(len(genes), -1)))
    real_dirs, real_cells, _ = feng_dirs.directions_from_stats(st, genes, cfg, REAL)
    gate["cm_at_both_ends"] = GATE_GROUP in proxy_dirs
    gate["direction_available"] = {read: g in proxy_dirs for read, g in READ_GROUPS.items()}
    gate["every_read_has_a_direction"] = all(gate["direction_available"].values())
    gate["mapped_genes"] = int((idx >= 0).sum())
    gate["mapped_genes_at_least"] = gate["mapped_genes"] >= MIN_MAPPED
    proceed = all(gate[k] for k in ("stages_present", "cm_at_both_ends", "every_read_has_a_direction",
                                    "mapped_genes_at_least"))
    modules = {m: set(hygiene.map_symbols(syms, genes)["mapped"]) for m, syms in
               hygiene.build_modules([str(g) for g in st["genes"]]).items()}
    held_set = set(t1_scope.held_flat_genes())
    held = np.array([g in held_set for g in genes], dtype=bool)
    groups = {}
    halves = feng_dirs.direction_halves(bracket, genes, cfg, PROXY) if bracket is not None else {}
    for g in harness.LP_GROUPS:
        if g in proxy_dirs and g in real_dirs:
            groups[g] = group_profile(g, proxy_dirs[g], real_dirs[g], halves[g], genes, held, modules)
    landing = landing_floor(st, bracket, genes, cfg)
    stop = t18_stop(groups)
    report = {"label": "N20: target-free profile of the late-proxy direction (no E10.5 file, no target, no score, no "
                       "board, no built file, no T1 cell count); it gates tries T13-T18 by its row's rule",
              "gate": gate, "proceed_to_tries": proceed,
              "landing_floor": landing, "proxy_close_to_real": landing["met"], "t18_stop": stop,
              "t18_runs": stop["runs"],
              "proxy_cells": {g: proxy_cells[g] for g in sorted(proxy_cells)},
              "real_cells": {g: real_cells[g] for g in sorted(real_cells)},
              "groups": groups, "n_t1_genes": len(genes),
              "inputs": {"t1_file": _rel(a.t1_file), "t1_file_sha256": hygiene.sha256_file(a.t1_file),
                         "feng_stats": _rel(path), "feng_stats_sha256": entry["sha256"],
                         "t1_hygiene_config_sha256": hygiene.sha256_file(hygiene.CONFIG),
                         "held_flat_sha256": hygiene.sha256_file(t1_scope.HELD_FLAT_GENES)},
              "stamp": {"commit": start[0], "tree_clean_at_start": None if start[1] is None else not start[1],
                        "tree_clean_at_end": None,
                        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    end_commit, end_unclean = provenance.git_state(provenance.REPO_ROOT)
    report["stamp"]["tree_clean_at_end"] = None if end_unclean is None else (not end_unclean and end_commit == start[0])
    out.mkdir(parents=True, exist_ok=True)
    (out / "n20_proxy_profile.json").write_text(json.dumps(report, indent=1) + "\n")
    lines = [f"N20: proceed to the tries: {'YES' if proceed else 'NO'} (gate {gate})",
             f"Landing floor (T16-T18): {landing}: {'met' if landing['met'] else 'not met (no recommendation for M1)'}",
             f"T18's stop: {stop}: {'T18 runs if the gate passes' if stop['runs'] else 'T18 does not run or count'}",
             "  group  genes  mean|logf|  moved>10%  spearman(real)  cosine  top400 same-sign  top400 overlap"]
    for g, r in groups.items():
        sp, co, ag = r["spearman_with_real"], r["cosine_with_real"], r["top400_sign_agreement_with_real"]
        lines.append(f"  {g}: {r['genes_with_value']}  {r['proxy_mean_abs_log_f']:.4f}  "
                     f"{r['proxy_genes_moved_gt_10pct']}  {sp if sp is None else round(sp, 3)}  "
                     f"{co if co is None else round(co, 3)}  {ag if ag is None else round(ag, 3)}  "
                     f"{r['top400_abs_overlap_with_real']}")
    lines += [f"Cells per group (Feng): proxy {report['proxy_cells']}; real {report['real_cells']}",
              f"Commit {report['stamp']['commit']}, clean at end {report['stamp']['tree_clean_at_end']}, "
              f"UTC {report['stamp']['utc']}"]
    text = "\n".join(lines) + "\n"
    (out / "n20_proxy_profile.txt").write_text(text)
    print(text, end="")
    if not proceed:
        sys.exit("N20's gate says no: tries T13-T18 are not run (they do not count)")


if __name__ == "__main__":
    main()
