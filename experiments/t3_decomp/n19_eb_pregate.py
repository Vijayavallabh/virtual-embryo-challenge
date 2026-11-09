#!/usr/bin/env python
"""Row N19 of experiments/harness/TRIES.md (T3; a diagnostic, not a try): the target-free pre-gate of tries K4 and K4c
(the empirical-Bayes shrinkage of the T3 source factor, experiments/t3_decomp/eb_factor.py). It asks, on each counted
fold's source, whether the shrinkage would change anything: if it is a no-op on both sources, K4 and K4c are pruned
without running.

Target-free and source-files only. For each of the two counted folds (the entries of vec.harness.K4_FOLDS with
counted=True) it reads only the fold's source: the knockout file(s) and their matched wild type, loaded as the harness
loads them (the fold's panel, `load_stage(drop_duplicates=True)`), and computes K2a's raw factor, K4's factor and K4c's
control factor with the code the candidates use (`vec.harness.t3_source_frame`, `t3_k4_fit`). Every file is read as its
own fold's source only: each knockout file is the other fold's target, and N19 compares nothing against a target and
computes no score.

The gate, per source (fixed in N19's row before it runs):
  predicted change: for the raw and the EB factor, dp_g = log1p(1e4 p_g phi_g) - log1p(1e4 p_g), with p the factor's
    own denominator read as proportions (den / sum(den), the w-weighted WT mean proportion over the masked labels): the
    factor's effect at the source WT's expression level.
  the top 100 genes by the raw factor's |dp| (fewer if the panel is smaller), and among them
    median_weight: the median shrinkage weight B (eb_factor.eb_weights: 1 = no shrinkage),
    n_moved: how many move more than 20 ranks between the raw and the EB |dp| rankings over the panel (ties by panel
    order).
  a source is a no-op if median_weight >= 0.95 and n_moved < 10.
  Verdict: "prune" if both sources are no-ops (K4's kill criterion 1: not run, not counted), otherwise "run".

Reported, deciding nothing: the sections used, the amplitude constant c, each arm's robust SD of log phi against the
raw factor's, the per-decile tau2, the quantiles of B over the panel, and the Spearman correlations of log phi with the
raw factor's. Per-gene values are not written (summaries only), and no cell count is printed or stored.

Outputs, in experiments/t3_decomp/results/n19_eb_pregate/ (no target-derived content): n19_eb_pregate.json and
n19_eb_pregate.txt, each with the commit, the clean-tree flag (vec.provenance), the UTC time and the inputs' sha256, and
a vec.provenance build stamp (<file>.stamp.json) each. N19 runs once: it refuses when its output folder already holds
anything (a rerun needs a new TRIES row).

Usage (DGX, repo root; after N19's row is committed):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/t3_decomp/n19_eb_pregate.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "n19_eb_pregate"

TOP = 100  # the top genes by the raw factor's |predicted dp|
RANK_MOVE = 20  # "moved": more than this many ranks between the raw and EB |dp| rankings
WEIGHT_NOOP = 0.95  # a source is a no-op if the median weight among the top genes is at least this
MOVED_NOOP = 10  # and fewer than this many of them moved


def predicted_change(phi, den) -> np.ndarray:
    """The factor's predicted change per gene at the source WT's expression level: log1p(1e4 p phi) - log1p(1e4 p),
    with p = den / sum(den), the factor's own denominator read as proportions."""
    den = np.asarray(den, dtype=np.float64)
    p = den / max(den.sum(), 1e-12)
    return np.log1p(1e4 * p * np.asarray(phi, dtype=np.float64)) - np.log1p(1e4 * p)


def ranks(values) -> np.ndarray:
    """Each gene's rank when `values` are sorted from the largest down (ties in panel order)."""
    order = np.argsort(-np.asarray(values, dtype=np.float64), kind="stable")
    out = np.empty(order.size, dtype=np.int64)
    out[order] = np.arange(order.size)
    return out


def source_gate(raw, eb, weights, den, top: int = TOP, rank_move: int = RANK_MOVE) -> dict:
    """N19's gate on one source (the module docstring): the top genes by the raw factor's |predicted dp|, their median
    shrinkage weight, how many of them move more than `rank_move` ranks, and whether the source is a no-op."""
    weights = np.asarray(weights, dtype=np.float64)
    dp_raw, dp_eb = predicted_change(raw, den), predicted_change(eb, den)
    rank_raw, rank_eb = ranks(np.abs(dp_raw)), ranks(np.abs(dp_eb))
    top_genes = np.flatnonzero(rank_raw < min(top, rank_raw.size))
    median_weight = float(np.median(weights[top_genes]))
    n_moved = int((np.abs(rank_raw[top_genes] - rank_eb[top_genes]) > rank_move).sum())
    return {"top": int(top_genes.size), "median_weight": median_weight, "n_moved": n_moved,
            "noop": bool(median_weight >= WEIGHT_NOOP and n_moved < MOVED_NOOP)}


def _spearman(a, b) -> float | None:
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.std() == 0 or b.std() == 0:
        return None
    return float(spearmanr(a, b).correlation)


def source_record(frame, fit) -> dict:
    """One source's gate and its reported, non-deciding quantities (summaries only: no per-gene value)."""
    from vec import harness

    eb = harness._t3_eb()
    gate = source_gate(fit.raw, fit.eb, fit.weights, fit.jack.den)
    log_raw, log_eb, log_ctrl = np.log(fit.raw), np.log(fit.eb), np.log(fit.ctrl)
    sd_raw = eb.robust_sd(log_raw)

    def amp_ratio(x) -> float | None:
        return None if sd_raw <= 0 else eb.robust_sd(x) / sd_raw

    quantiles = ("0.1", "0.25", "0.5", "0.75", "0.9")
    return {"gate": gate, "sections_used": int(fit.jack.n_blocks), "levels_over_files": int(frame.n_levels),
            "amplitude_constant": float(fit.constant),
            "robust_sd": {"raw": sd_raw, "eb": eb.robust_sd(log_eb), "ctrl": eb.robust_sd(log_ctrl)},
            "robust_sd_ratio_to_raw": {"eb": amp_ratio(log_eb), "ctrl": amp_ratio(log_ctrl)},
            "tau2_by_decile": [float(x) for x in fit.tau2],
            "weight_quantiles": {q: float(np.quantile(fit.weights, float(q))) for q in quantiles},
            "spearman_log_phi": {"raw_eb": _spearman(log_raw, log_eb), "raw_ctrl": _spearman(log_raw, log_ctrl)}}


def _fmt(x: float | None) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def verdict(sources: dict[str, dict]) -> str:
    """K4's kill criterion 1: "prune" if both sources are no-ops, otherwise "run"."""
    return "prune" if all(r["gate"]["noop"] for r in sources.values()) else "run"


def report_text(res: dict) -> str:
    lines = ["N19: the target-free pre-gate of tries K4 and K4c (experiments/t3_decomp/eb_factor.py)",
             f"Verdict: {res['verdict']} (prune only if both sources are no-ops: median weight >= {WEIGHT_NOOP} and "
             f"fewer than {MOVED_NOOP} of the top {TOP} genes move more than {RANK_MOVE} ranks)", ""]
    for key, r in res["sources"].items():
        g = r["gate"]
        lines += [f"{key} source: median weight among the top {g['top']} genes {g['median_weight']:.3f}, moved "
                  f">{RANK_MOVE} ranks: {g['n_moved']}, no-op: {g['noop']}",
                  f"  sections used {r['sections_used']}, amplitude constant {r['amplitude_constant']:.3f}, robust SD "
                  f"of log phi raw {r['robust_sd']['raw']:.4f}, EB {r['robust_sd']['eb']:.4f} "
                  f"({_fmt(r['robust_sd_ratio_to_raw']['eb'])} x raw), control {r['robust_sd']['ctrl']:.4f} "
                  f"({_fmt(r['robust_sd_ratio_to_raw']['ctrl'])} x raw)",
                  "  tau2 by decile of the denominator (lowest first): "
                  + ", ".join(f"{x:.4g}" for x in r["tau2_by_decile"]),
                  "  weight quantiles over the panel: "
                  + ", ".join(f"{q}: {v:.3f}" for q, v in r["weight_quantiles"].items()),
                  f"  Spearman of log phi with the raw factor's: EB {r['spearman_log_phi']['raw_eb']}, control "
                  f"{r['spearman_log_phi']['raw_ctrl']}", ""]
    s = res["stamp"]
    lines += [f"Inputs (sha256): {res['inputs']}",
              f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="N19: the target-free pre-gate of K4 and K4c")
    ap.add_argument("--data", default=None, help="the released data (default: data/raw)")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N19 has already written its outputs; a rerun needs a new TRIES row")

    from vec import harness, provenance
    from vec.boards import load_boards
    from vec.prep import load_stage

    started = provenance.git_state(provenance.REPO_ROOT)
    data = Path(a.data) if a.data else harness.DATA_DIR
    boards = load_boards()
    harness.check_index(boards)
    tasks = {t.key: t for t in harness.load_tasks()["T3:bcatenin"]}
    rules = harness._t3_head().load_rules()
    sources, inputs = {}, {}
    for key, fold in harness.K4_FOLDS.items():
        if not fold["counted"]:
            continue
        t = tasks[key]
        genes = boards[t.panel_board].genes
        roles = (*fold["ko_roles"], "wt_train")
        loaded = {r: load_stage(data / t.inputs[r], genes=genes, drop_duplicates=True) for r in roles}
        frame = harness.t3_source_frame(loaded["wt_train"], [loaded[r] for r in fold["ko_roles"]], rules,
                                        fold["source_preset"])
        sources[key] = source_record(frame, harness.t3_k4_fit(frame))
        inputs.update({t.inputs[r]: harness._sha256(data / t.inputs[r]) for r in roles})
    ended = provenance.git_state(provenance.REPO_ROOT)
    clean = None if started[1] is None or ended[1] is None else (not started[1] and not ended[1]
                                                                 and started[0] == ended[0])
    res = {"label": "N19 (experiments/harness/TRIES.md): the target-free pre-gate of K4 and K4c; source files only",
           "verdict": verdict(sources), "sources": sources, "inputs": inputs,
           "settings": {"top": TOP, "rank_move": RANK_MOVE, "weight_noop": WEIGHT_NOOP, "moved_noop": MOVED_NOOP,
                        "bins": harness.K4_BINS, "control_seed": harness.K4_CTRL_SEED},
           "stamp": {"commit": ended[0], "tree_clean": clean,
                     "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    text = report_text(res)
    paths = {out / "n19_eb_pregate.json": json.dumps(res, indent=1) + "\n", out / "n19_eb_pregate.txt": text}
    for p, body in paths.items():
        p.write_text(body)
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
