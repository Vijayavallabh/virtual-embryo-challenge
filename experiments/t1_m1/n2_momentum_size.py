#!/usr/bin/env python
"""Row N2 of experiments/harness/TRIES.md (T1; a diagnostic, not a try): how far M1's momentum moves the cells that
have no Feng group, the size of what ruling R1 would remove. Target-free: it reads the test build's two inputs, the
released E9.5 and the unlocked E10.5, as `build_test.py t1` loads them, and Feng's statistics only to learn which
groups carry a direction. It reads no target, no score, no board and no built file, and it prints and writes no cell
count.

What M1 does to these cells at test (build_t1 -> vec.t1_scope.build_submission):
- Which cells. A drawn cell whose label maps (feng_dirs.label_groups) to a group with a Feng direction
  (feng_dirs.feng_directions) moves along it at 2/5 of the bracket; every other drawn cell gets `momentum_shift` with
  its chamber category's factor. The momentum categories and their labels are derived here from that code.
- The factor of category t on gene g: log f = B0_DAMPING x [log(m_E10.5 + eps) - log(m_E9.5 + eps)], m the mean CP10k
  over all of t's cells at that stage (its labels with a Feng direction included), with momentum_shift's defaults,
  which build_t1 leaves alone (eps 0.5, no clip); f = 1 on the held-flat genes, and for a category with fewer than
  min_cells cells at a stage or drawn from stand-ins.
- vec.core.apply_factors applies f in count space: a count c becomes Binomial(c, f s) or c + Poisson((f s - 1) c), s
  the cell's keep-depth constant (sum c / sum f c over its genes not held flat). Zeros stay zero; held-flat counts stay.

Per momentum category it reports: the genes whose factor moves by more than 10% and 25% (|f - 1| > 0.10, > 0.25)
among the genes with at least one count in the E10.5 cells the momentum moves (the denominator; held-flat genes
included, at f = 1), and the mean |log f| over them; for context, those cells' CP10k share on the moving genes, the
keep-depth constant and the same counts for f s at its median; and the share of s1's and s2's cells that the category
holds and that the momentum moves, from E10.5's chamber mix (harness.t1_chamber_shares, as build_t1) and the frozen
s2 rule (docs/S2_RULES.md, `--cm-share`, applied by t1_scope.with_share), to two decimals. Before anything is written,
the factors are checked against momentum_shift's own report on the same cells.

Usage (DGX, repo root):
  OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
      uv run --frozen --no-sync python -I -B experiments/t1_m1/n2_momentum_size.py \
      --out experiments/t1_m1/results/n2_momentum_size
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec.boards import REPO_ROOT
from vec.core import counts_of, mean_cp10k
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import feng_dirs  # noqa: E402

REF, START = "E10.5_RNA", "E9.5_RNA"  # build_t1's reference and start files
KEY, LABEL_KEY = t1_scope.CATEGORY_KEY, t1_scope.LABEL_KEY
THRESHOLDS = (0.10, 0.25)  # |f - 1| above these (the row)
S2_RULES, S2_PIN = REPO_ROOT / "docs" / "S2_RULES.md", REPO_ROOT / "docs" / "S2_RULES.sha256"
# momentum_shift's settings that build_t1 leaves at their defaults (it passes no shift_kw)
SHIFT = {k: p.default for k, p in inspect.signature(t1_scope.momentum_shift).parameters.items()
         if k in ("eps", "min_cells", "max_abs_log_factor", "keep_depth", "standin_momentum")}


def s2_cm_share(rules: Path = S2_RULES, pin: Path = S2_PIN) -> float:
    """T1's s2 rule, the cardiomyocyte share of heart cells (`--cm-share`), read from the frozen rules after their
    sha256 check (docs/S2_RULES.md, "Before each s2 build")."""
    data = rules.read_bytes()
    if hashlib.sha256(data).hexdigest() != pin.read_text().split()[0]:
        sys.exit(f"{rules} differs from its frozen sha256 in {pin}")
    found = set(re.findall(r"--cm-share (\d+(?:\.\d+)?)", data.decode()))
    if len(found) != 1:
        sys.exit(f"{rules}: expected one --cm-share value for T1's s2, found {sorted(found)}")
    return float(found.pop())


def category_logf(S, R, start, ref, t: str, held: np.ndarray, spec) -> tuple[np.ndarray, str]:
    """momentum_shift's log factor for category t (its rows, formula and settings) and its status. S, R: the counts of
    start and ref on ref's genes."""
    mc = SHIFT["min_cells"]
    s_rows, s_how = t1_scope._type_rows(start, t, KEY, LABEL_KEY, None, spec.fallbacks(START), mc)
    r_rows, r_how = t1_scope._type_rows(ref, t, KEY, LABEL_KEY, None, spec.fallbacks(REF), mc)
    if len(s_rows) < mc or len(r_rows) < mc:
        return np.zeros(R.shape[1]), "skipped: fewer than min_cells cells at a stage"
    if not SHIFT["standin_momentum"] and (s_how.startswith("stand-in") or r_how.startswith("stand-in")):
        return np.zeros(R.shape[1]), "skipped: stand-in cells at a stage"
    eps = SHIFT["eps"]
    logf = harness.B0_DAMPING * (np.log(mean_cp10k(R[r_rows]) + eps) - np.log(mean_cp10k(S[s_rows]) + eps))
    if SHIFT["max_abs_log_factor"] is not None:
        logf = np.clip(logf, -SHIFT["max_abs_log_factor"], SHIFT["max_abs_log_factor"])
    logf[held] = 0.0
    return logf, "shifted"


def keep_depth_constant(C, logf: np.ndarray, held: np.ndarray) -> np.ndarray:
    """apply_factors' keep-depth constant of each cell (rows of counts C): sum c / sum f c over the genes not held."""
    if not SHIFT["keep_depth"]:
        return np.ones(C.shape[0])
    free = (~held).astype(np.float64)
    num, den = C @ free, C @ (np.exp(logf) * free)
    return np.divide(num, den, out=np.ones(C.shape[0]), where=den > 0)


def factor_size(logf: np.ndarray, expressed: np.ndarray, held: np.ndarray, cp10k: np.ndarray,
                scale: float = 1.0) -> dict:
    """The factor exp(logf), times `scale` on the genes not held flat, among the `expressed` genes (the denominator):
    genes with |f - 1| above each threshold (up and down), the mean |log f| over the denominator and, with cp10k the
    moved cells' mean CP10k per gene, the share of their CP10k on those genes."""
    lf = np.where(held, 0.0, logf + np.log(scale))[expressed]
    f, w = np.exp(lf), cp10k[expressed]
    n = int(expressed.sum())
    out = {"denominator_genes": n, "held_flat_in_denominator": int((held & expressed).sum()),
           "mean_abs_log_f": float(np.abs(lf).mean()) if n else 0.0,
           "max_abs_log_f": float(np.abs(lf).max()) if n else 0.0}
    for thr in THRESHOLDS:
        m = np.abs(f - 1.0) > thr
        out[f"gt_{thr:.2f}"] = {"genes": int(m.sum()), "up": int((m & (f > 1)).sum()), "down": int((m & (f < 1)).sum()),
                                "share_of_denominator": float(m.sum() / n) if n else 0.0,
                                "cp10k_share": float(w[m].sum() / w.sum()) if w.sum() > 0 else 0.0}
    return out


def cell_shares(shares: dict, moved_fraction: dict) -> dict:
    """Per momentum category, the share of all cells it holds and the share the momentum moves (its share times the
    fraction of its E10.5 pool whose labels have no Feng direction; a category's cells are drawn uniformly from its
    pool), with their totals, to two decimals."""
    rows = {c: {"category": float(shares.get(c, 0.0)), "moved": float(shares.get(c, 0.0)) * fr}
            for c, fr in moved_fraction.items()}
    rows["total"] = {k: sum(r[k] for r in rows.values()) for k in ("category", "moved")}
    return {c: {k: round(v, 2) for k, v in r.items()} for c, r in rows.items()}


def label_mix(labels: np.ndarray) -> dict:
    """Each label's share of a category's cells at a stage, to two decimals."""
    u, n = np.unique(labels, return_counts=True)
    return {str(lab): round(float(k) / len(labels), 2) for lab, k in zip(u, n, strict=True)}


def expressed(C, rows: np.ndarray) -> np.ndarray:
    """Genes with at least one count in the given rows of counts C."""
    return np.diff(C[np.flatnonzero(rows)].tocsc().indptr) > 0


def check_against_momentum_shift(start, ref, rows: np.ndarray, held_genes: list[str], spec, mine: dict) -> None:
    """momentum_shift itself on the moved E10.5 cells, with build_t1's arguments: its report must match the factors and
    keep-depth constants computed here, or the script stops (it would describe another factor than the build's)."""
    _, rep = t1_scope.momentum_shift(start, ref, ref[rows].copy(), harness.B0_DAMPING, held_genes, label_key=LABEL_KEY,
                                     start_fallback=spec.fallbacks(START), ref_fallback=spec.fallbacks(REF), seed=0)
    if {r["type"] for r in rep.rows} != set(mine):
        sys.exit("momentum_shift moved other categories than the ones derived here")
    for r in rep.rows:
        logf, status, s = mine[r["type"]]
        a = np.abs(logf)
        same = ((r["status"] == "shifted") == (status == "shifted")
                and (r["n_up"], r["n_down"]) == (int((logf > np.log(1.1)).sum()), int((logf < -np.log(1.1)).sum()))
                and np.allclose([r["mean_abs_log_f"], r["max_abs_log_f"], r["median_depth_scale"]],
                                [a.mean(), a.max(), np.median(s)], rtol=1e-9, atol=1e-12))
        if not same:
            sys.exit(f"{r['type']}: the factor computed here differs from momentum_shift's report")


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], capture_output=True, text=True, cwd=REPO_ROOT,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"})


def _sizes(d: dict, what: str) -> str:
    """One line of factor sizes (genes are counted, never cells)."""
    a, b = d["gt_0.10"], d["gt_0.25"]
    return (f"  {what}: {d['denominator_genes']:,} genes ({d['held_flat_in_denominator']} held flat); "
            f"|f - 1| > 0.10: {a['genes']:,} ({a['share_of_denominator']:.1%}; {a['up']:,} up, {a['down']:,} down); "
            f"> 0.25: {b['genes']:,} ({b['share_of_denominator']:.1%}; {b['up']:,} up, {b['down']:,} down); "
            f"mean |log f| {d['mean_abs_log_f']:.3f}; the cells' CP10k on those genes {a['cp10k_share']:.2f} and "
            f"{b['cp10k_share']:.2f}")


def _mix(m: dict) -> str:
    return ", ".join(f"{lab} {v:.2f}" for lab, v in m.items()) or "none"


def report_text(res: dict) -> str:
    st = res["settings"]
    clip = "none" if st["max_abs_log_factor"] is None else st["max_abs_log_factor"]
    lines = [f"N2 (experiments/harness/TRIES.md; T1, a diagnostic, not a try): {res['label']}",
             f"Feng groups with a direction: {', '.join(res['feng_groups_with_direction'])}; momentum categories: "
             f"{', '.join(res['categories'])}",
             f"Factor per category: log f = {st['damping']} x [log(m_E10.5 + {st['eps']}) - log(m_E9.5 + {st['eps']})]"
             f", clip {clip}, f = 1 on {st['held_flat_genes']} held-flat genes, keep_depth {st['keep_depth']}; s2: "
             f"cardiomyocytes at {st['s2_cm_share']:.2f} of heart cells (docs/S2_RULES.md)"]
    for c, r in res["categories"].items():
        sh, s, fs = r["share_of_cells"], r["keep_depth_constant"], r["factor_at_median_s"]
        a, b = fs["gt_0.10"], fs["gt_0.25"]
        lines += [f"{c} ({r['status']}): moved labels {', '.join(r['moved_labels'])} "
                  f"({r['moved_fraction_of_category']:.2f} of the category at E10.5); on a Feng direction: "
                  f"{', '.join(r['feng_labels']) or 'none'}; label mix E9.5 {_mix(r['label_mix']['E9.5'])}; E10.5 "
                  f"{_mix(r['label_mix']['E10.5'])}",
                  _sizes(r["factor"], "genes expressed in the moved E10.5 cells"),
                  f"  keep-depth constant s: median {s['median']:.2f} (10% {s['q10']:.2f}, 90% {s['q90']:.2f}); at "
                  f"the median, |f s - 1| > 0.10: {a['genes']:,} genes ({a['up']:,} up, {a['down']:,} down), > 0.25: "
                  f"{b['genes']:,} ({b['up']:,} up, {b['down']:,} down)"]
        if "factor_whole_category" in r:
            lines.append(_sizes(r["factor_whole_category"], "genes expressed in any of the category's E10.5 cells"))
        lines.append(f"  share of cells: s1 {sh['s1']['category']:.2f} in the category, {sh['s1']['moved']:.2f} moved "
                     f"by the momentum; s2 {sh['s2']['category']:.2f} and {sh['s2']['moved']:.2f}")
    tot = res["share_of_cells_total"]
    lines.append(f"Moved by the momentum, all categories: s1 {tot['s1']['moved']:.2f} of its cells, s2 "
                 f"{tot['s2']['moved']:.2f}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    ref = load_stage(Path(a.data) / "T1" / f"{REF}.h5ad", lineage=False)  # as build_t1
    start = load_stage(Path(a.data) / "T1" / f"{START}.h5ad", lineage=False)
    spec = t1_scope.load_spec()
    s2_cm = s2_cm_share()
    s1 = harness.t1_chamber_shares(harness._shallow(ref), REF)
    shares = {"s1": s1, "s2": t1_scope.with_share(s1, "cm", s2_cm)}
    t1_scope.classify(ref, LABEL_KEY, REF, spec)
    t1_scope.classify(start, LABEL_KEY, START, spec)
    genes = list(ref.var_names)
    directions, _ = feng_dirs.feng_directions(genes)  # its record holds Feng's cell counts: not kept
    groups = feng_dirs.label_groups(ref.obs[LABEL_KEY].astype(str).unique())
    labels, cats = ref.obs[LABEL_KEY].astype(str).to_numpy(), ref.obs[KEY].astype(str).to_numpy()
    s_labels, s_cats = start.obs[LABEL_KEY].astype(str).to_numpy(), start.obs[KEY].astype(str).to_numpy()
    covered = np.array([groups.get(lab) in directions for lab in labels])  # build_submission's split
    held_genes = t1_scope.held_flat_genes()
    held = np.isin(genes, held_genes)
    R, S = counts_of(ref), counts_of(start, genes)
    moving = [c for c in s1 if (~covered & (cats == c)).any()]  # s1: the chamber categories with cells
    mine, out, fraction = {}, {}, {}
    for c in moving:
        pool, moved = cats == c, (cats == c) & ~covered
        logf, status = category_logf(S, R, start, ref, c, held, spec)
        s = keep_depth_constant(R[np.flatnonzero(moved)], logf, held)
        mine[c], fraction[c] = (logf, status, s), moved.sum() / pool.sum()
        expr, cp = expressed(R, moved), mean_cp10k(R[np.flatnonzero(moved)])
        out[c] = {"status": status, "moved_labels": sorted(set(labels[moved])),
                  "feng_labels": sorted(set(labels[pool & covered])),
                  "label_mix": {"E9.5": label_mix(s_labels[s_cats == c]), "E10.5": label_mix(labels[pool])},
                  "moved_fraction_of_category": round(float(fraction[c]), 2),
                  "factor": factor_size(logf, expr, held, cp),
                  "keep_depth_constant": {k: float(np.percentile(s, q)) for k, q in
                                          (("median", 50), ("q10", 10), ("q90", 90))},
                  "factor_at_median_s": factor_size(logf, expr, held, cp, float(np.median(s)))}
        if (pool & covered).any():  # the momentum moves only part of the category
            out[c]["factor_whole_category"] = factor_size(logf, expressed(R, pool), held,
                                                          mean_cp10k(R[np.flatnonzero(pool)]))
    check_against_momentum_shift(start, ref, np.flatnonzero(np.isin(cats, moving) & ~covered), held_genes, spec, mine)
    by_entry = {e: cell_shares(sh, fraction) for e, sh in shares.items()}
    for c in moving:
        out[c]["share_of_cells"] = {e: by_entry[e][c] for e in by_entry}
    head = _git("rev-parse", "HEAD")
    status = _git("status", "--porcelain", "--untracked-files=no")
    res = {"label": "the size of M1's momentum on cells without a Feng group (target-free: no target, no score, no "
                    "board, no built file, no cell count)",
           "inputs": {"start": f"T1/{START}.h5ad", "reference": f"T1/{REF}.h5ad"},
           "settings": {"damping": harness.B0_DAMPING, **{k: v for k, v in SHIFT.items() if k != "min_cells"},
                        "held_flat_genes": int(held.sum()),
                        "thresholds": list(THRESHOLDS), "s2_cm_share": s2_cm,
                        "s2_rules_sha256": S2_PIN.read_text().split()[0]},
           "feng_groups_with_direction": sorted(directions),
           "label_groups": {lab: groups[lab] for lab in sorted(set(labels[np.isin(cats, list(s1))]))},
           "categories": out,
           "share_of_cells_total": {e: by_entry[e]["total"] for e in by_entry},
           "provenance": {"commit": head.stdout.strip() if head.returncode == 0 else None,
                          "tree_clean": (status.stdout.strip() == "") if status.returncode == 0 else None,
                          "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}}
    d = Path(a.out)
    d.mkdir(parents=True, exist_ok=True)
    (d / "n2_momentum_size.json").write_text(json.dumps(res, indent=1) + "\n")
    text = report_text(res)
    (d / "n2_momentum_size.txt").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
