#!/usr/bin/env python
"""Diagnostic N29 of experiments/harness/TRIES.md (T3; target-free, a diagnostic and not a try): is there a
per-lineage-group part of the wild-type-only revert factor of tries K7 beyond its pooled part?

Why. T3's first-choice file (the allowed branch's s1) applies the Gata4 factor, measured in the non-Mesp1 compartments
against WT E8.75, so it already holds the pooled part of WT E8.75's own deviation. A reviewer proposes adding to s1
only each lineage group's deviation beyond the pooled part. If that residual is small next to the per-group factor
itself, the proposal is a no-op and is pruned without a try.

What it reads. The seven released wild-type MERFISH stage files with the full panel (k7_wt_revert.FULL_PANEL: WT E8.75,
the reference, and the six others), loaded exactly as the harness loads K7's inputs (vec.prep.load_stage on the panel of
K7's tasks, drop_duplicates=True). It computes the factors with K7's own functions, unchanged
(k7_wt_revert.group_means, trend_at and revert_factors). It opens no knockout file and holds no path to one, compares
nothing with a fold or backtest target, builds no candidate and no test-board file, and computes no score.

The statistic (every constant below is a module constant):
  r_g     the per-gene log factor K7b applies to lineage group g, before any exponent (log of revert_factors' per-group
          factor), for every group K7b covers (the groups revert_factors returns).
  r_pool  the per-gene log factor K7a applies (log of revert_factors' pooled factor).
  e_g     r_g - r_pool.
  w_g     group g's share of WT E8.75's cells among the covered groups (the cells K7 weights by: those with counts).
  F       sqrt(sum_g w_g |e_g|^2) / sqrt(sum_g w_g |r_g|^2), norms over the panel's genes (K7b's factor is defined on
          every panel gene).
GATE, fixed in N29's row before the run: verdict "prune" if F < F_PRUNE (0.20: the per-group factor is the pooled one
within a fifth of its size, so the residual is a no-op); otherwise "open" (a counted try on the two cross-knockout folds
may be registered; its row would also fix a second criterion, the residual's correlation with s1's own factor, which
N29 does not read because that factor comes from a knockout file).

Reported, deciding nothing: F over the Mesp1-lineage groups (configs/lineage_map.yaml flag mesp1 "yes"), over the
strict non-Mesp1 groups (rules.yaml lineage.nonmesp1_groups) and over the rest (flag "partial", or "no" outside that
list); per group |e_g| / |r_g| and the cosine of r_g with r_pool; the number of groups; the SD of r_pool over the
genes; and F_nm, the same F with the pooled part taken over the strict non-Mesp1 groups only (t3_decompose.gene_factor,
the function K7 calls, with that group mask), the compartments s1's factor is measured in. Group names are written; no
per-gene value, no gene name, no group share and no cell number.

Outputs, in experiments/t3_decomp/results/n29_residual_gate/: n29.json and n29.txt, each with the commit, the
clean-tree flag, the UTC time and the inputs' sha256, and a vec.provenance build stamp (<file>.stamp.json) each. N29
runs once: it refuses when its output folder already holds anything (a rerun needs a new TRIES row).

Usage (repo root, after N29's row is committed):
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/t3_decomp/n29_residual_gate.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (the factor of tries K7, unchanged)

OUT = HERE / "results" / "n29_residual_gate"

F_PRUNE = 0.20  # the gate: "prune" if F < F_PRUNE, otherwise "open"
REFERENCE = "E8.75"  # the reference of K7's Gata4 tasks and of the test (WT E8.75)
PANEL_BOARD = "T3:gata4"  # the panel K7's tasks g1 and g2 load their inputs on (their twin and validation board)
MESP1_FLAG = "yes"  # configs/lineage_map.yaml: a group is Mesp1-lineage if its mesp1 flag is this
VERDICT_PRUNE, VERDICT_OPEN = "prune", "open"


def stage_file(stage: str) -> str:
    """A released wild-type MERFISH stage file, relative to the data folder."""
    return f"T2/{stage}.h5ad"


def weights(sizes: dict[str, int], groups: list[str]) -> dict[str, float]:
    """Each group's share of the reference's cells among `groups`."""
    n = np.array([float(sizes[g]) for g in groups])
    if n.sum() <= 0:
        raise ValueError("the covered groups hold no cells")
    return {g: float(x) for g, x in zip(groups, n / n.sum(), strict=True)}


def fraction(r: dict[str, np.ndarray], r_pool: np.ndarray, w: dict[str, float], groups=None) -> float | None:
    """F = sqrt(sum_g w_g |r_g - r_pool|^2) / sqrt(sum_g w_g |r_g|^2) over `groups` (default: every group of `w`);
    None when no group is given or the denominator is zero. The shares need not sum to one over a subset: they
    cancel."""
    groups = list(w) if groups is None else [g for g in groups if g in w]
    if not groups:
        return None
    num = sum(w[g] * float(np.sum((r[g] - r_pool) ** 2)) for g in groups)
    den = sum(w[g] * float(np.sum(r[g] ** 2)) for g in groups)
    return None if den <= 0 else float(np.sqrt(num / den))


def cosine(a: np.ndarray, b: np.ndarray) -> float | None:
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    return None if na == 0 or nb == 0 else float(a @ b / (na * nb))


def verdict(f: float | None) -> str:
    """The gate fixed in N29's row: "prune" if F < F_PRUNE, otherwise "open" (an undefined F does not prune)."""
    return VERDICT_PRUNE if f is not None and f < F_PRUNE else VERDICT_OPEN


def lineage_sets(groups: list[str], flags: dict[str, str], nonmesp1: list[str]) -> dict[str, list[str]]:
    """The covered groups split by the lineage flag: Mesp1-lineage (flag "yes"), strict non-Mesp1 (rules.yaml's list),
    and the rest."""
    mesp1 = [g for g in groups if flags.get(g) == MESP1_FLAG]
    strict = [g for g in groups if g in set(nonmesp1)]
    return {"mesp1": mesp1, "nonmesp1": strict, "other": [g for g in groups if g not in set(mesp1) | set(strict)]}


def factors(ref, others: list, ref_day: float, other_days: list[float]) -> dict:
    """K7's factors from the reference's and the other samples' AnnData, through K7's own functions: the covered
    groups, r_g and r_pool (logs, before any exponent), the reference's group sizes, and the stacks K7 builds the
    pooled factor from (for the reported non-Mesp1-pooled form)."""
    ref_means, ref_sizes = k7.group_means(ref)
    trend = k7.trend_at([k7.group_means(a)[0] for a in others], other_days, ref_day)
    phi, per, rec = k7.revert_factors(ref_means, ref_sizes, trend)
    groups = list(rec["groups"])
    return {"groups": groups, "r": {g: np.log(per[g]) for g in groups}, "r_pool": np.log(phi),
            "sizes": {g: ref_sizes[g] for g in groups},
            "mu_trend": np.stack([trend[g] for g in groups]), "mu_ref": np.stack([ref_means[g] for g in groups])}


def nonmesp1_pool(fx: dict, strict: list[str]) -> np.ndarray | None:
    """log of t3_decompose.gene_factor (K7's call, same eps and clip) of the trend against the reference over the
    strict non-Mesp1 groups only, weighted by the reference's group sizes; None if no such group is covered."""
    from vec import harness

    mask = np.array([g in set(strict) for g in fx["groups"]])
    if not mask.any():
        return None
    w = np.array([float(fx["sizes"][g]) for g in fx["groups"]])
    phi, _ = harness._t3_head().t3d.gene_factor(fx["mu_trend"], fx["mu_ref"], w / w.sum(), mask, k7.EPS, k7.CLIP)
    return np.log(phi)


def summarise(fx: dict, flags: dict[str, str], nonmesp1: list[str]) -> dict:
    """N29's statistic, its gate and the reported quantities: summaries only (no per-gene value, no share, no count)."""
    groups, r, r_pool = fx["groups"], fx["r"], fx["r_pool"]
    w = weights(fx["sizes"], groups)
    f = fraction(r, r_pool, w)
    sets = lineage_sets(groups, flags, nonmesp1)
    r_nm = nonmesp1_pool(fx, sets["nonmesp1"])
    per = {}
    for g in groups:
        nr = float(np.linalg.norm(r[g]))
        per[g] = {"lineage": next(k for k, v in sets.items() if g in v),
                  "residual_ratio": None if nr == 0 else float(np.linalg.norm(r[g] - r_pool) / nr),
                  "cosine_with_pool": cosine(r[g], r_pool)}
    return {"F": f, "gate": {"F_prune": F_PRUNE, "rule": f'"prune" if F < {F_PRUNE}, otherwise "open"'},
            "verdict": verdict(f),
            "reported": {"F_by_lineage": {k: fraction(r, r_pool, w, v) for k, v in sets.items()},
                         "groups_by_lineage": sets,
                         "F_nonmesp1_pool": None if r_nm is None else fraction(r, r_nm, w),
                         "n_groups": len(groups), "n_genes": int(r_pool.size),
                         "sd_r_pool": float(r_pool.std()), "per_group": per}}


def _fmt(x: float | None) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def report_text(res: dict) -> str:
    rep = res["reported"]
    lines = ["N29 (T3; target-free diagnostic, not a try): the per-group part of K7's wild-type revert factor beyond "
             "its pooled part",
             f"Verdict: {res['verdict']} (F = {_fmt(res['F'])}; gate: {res['gate']['rule']})", "",
             f"Reported, deciding nothing: {rep['n_groups']} groups covered, {rep['n_genes']} panel genes, SD of the "
             f"pooled log factor {rep['sd_r_pool']:.4f}",
             "  F by lineage: " + ", ".join(f"{k} {_fmt(v)} ({', '.join(rep['groups_by_lineage'][k]) or 'none'})"
                                           for k, v in rep["F_by_lineage"].items()),
             f"  F with the pooled part over the strict non-Mesp1 groups only: {_fmt(rep['F_nonmesp1_pool'])}",
             "  per group: |e_g|/|r_g|, cosine of r_g with r_pool"]
    lines += [f"    {g:<20} {v['lineage']:<9} {_fmt(v['residual_ratio'])}  {_fmt(v['cosine_with_pool'])}"
              for g, v in rep["per_group"].items()]
    s = res["stamp"]
    lines += ["", f"Inputs (sha256): {res['inputs']}",
              f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="N29: the per-group residual of K7's wild-type revert factor")
    ap.add_argument("--data", default=None, help="the released data (default: data/raw)")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N29 has already written its outputs; a rerun needs a new TRIES row")

    from vec import harness, provenance
    from vec.boards import load_boards
    from vec.prep import lineage_map, load_stage

    started = provenance.git_state(provenance.REPO_ROOT)
    data = Path(a.data) if a.data else harness.DATA_DIR
    boards = load_boards()
    harness.check_index(boards)
    genes = boards[PANEL_BOARD].genes
    files = {s: stage_file(s) for s in k7.FULL_PANEL}
    others = [s for s in k7.FULL_PANEL if s != REFERENCE]
    loaded = {s: load_stage(data / f, genes=genes, drop_duplicates=True) for s, f in files.items()}
    for s, st in loaded.items():
        if [g for g in genes if g not in st.var_names]:
            raise ValueError(f"{files[s]} lacks panel genes")
    fx = factors(loaded[REFERENCE], [loaded[s] for s in others], harness.STAGE_DAYS[files[REFERENCE]],
                 [harness.STAGE_DAYS[files[s]] for s in others])
    flags = {g: str(v["mesp1"]) for g, v in lineage_map()["groups"].items()}
    nonmesp1 = list(harness._t3_head().load_rules()["lineage"]["nonmesp1_groups"])
    res = summarise(fx, flags, nonmesp1)
    ended = provenance.git_state(provenance.REPO_ROOT)
    clean = None if started[1] is None or ended[1] is None else (not started[1] and not ended[1]
                                                                 and started[0] == ended[0])
    res = {"label": "N29 (experiments/harness/TRIES.md): the per-group residual of K7's wild-type revert factor; "
                    "released wild-type stages only", "verdict": res.pop("verdict"), **res,
           "factor_rule": k7.WT_REVERT_RULE, "reference": REFERENCE, "panel_board": PANEL_BOARD,
           "inputs": {f: harness._sha256(data / f) for f in files.values()},
           "stamp": {"commit": ended[0], "tree_clean": clean,
                     "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    text = report_text(res)
    for p, body in {out / "n29.json": json.dumps(res, indent=1) + "\n", out / "n29.txt": text}.items():
        p.write_text(body)
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
