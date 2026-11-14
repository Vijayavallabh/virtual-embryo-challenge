#!/usr/bin/env python
"""Diagnostic N26 (experiments/harness/TRIES.md), the target-free gate of heart tries H7 and H7c: does the
within-lineage factor H7 would apply to a stage's cells agree with what those cells' lineage did next?

It reads the released heart stages E8.25_late, E8.5, E8.75 and E9.5 only. E10.5, the one target of both heart
backtests, is never opened, so nothing here compares a predictor with a counted target: a diagnostic, not a try.

For a triple of consecutive stages (A, B, C):
  predicted  H7's own finished factors for B's cells, `harness.t2_gmom_pair_factors(A, B)`: metacells on B, the log
             ratio of each metacell's mean to its transported A partners' mean, the lineage's cell-weighted mean
             removed per gene, damped, zeroed and clipped. This is what H7 adds to B's cells when B is the reference.
  realised   the next step measured by the same machinery, `ts.fit_state_factors(B, C, "earlier", ...)`: metacells on
             B, the log ratio of the transported C partners' mean to the metacell's own mean, raw. It is carried onto
             the predicted fit's metacells cell by cell (each B cell takes the vector of its realised metacell; a
             predicted metacell's realised vector is the mean over its cells), and the lineage's cell-weighted mean is
             removed per gene.
  agreement  per lineage, the Pearson correlation over (metacell, gene) entries, each metacell weighted by its cells,
             both sides centred per gene. Per triple, S = the cell-weighted mean of the lineage correlations.
  null       the predicted metacell rows permuted at random within each lineage, N_PERM times (one generator per
             triple): the share of permutations whose S is at least the observed S is p.
Triples: (E8.25_late, E8.5, E8.75) and (E8.5, E8.75, E9.5), the two that end at or before E9.5.

Gate, fixed in the row before the run: H7 and H7c run as counted tries only if S > 0 with p <= 0.05 on BOTH triples.
Anything else: both are dropped, not run and not counted. Reported, deciding nothing: each lineage's correlation, the
same statistic with the predicted factors before damping, zeroing and clipping, and N21's two concordance numbers
between the two consecutive fits (leading-axis |cosine|, per-gene size Spearman).

  python experiments/harness/n26_heart_persistence.py --out experiments/harness/results/n26_heart_persistence
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from vec import harness
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT, load_boards
from vec.prep import load_stage

sys.path.insert(0, str(REPO_ROOT / "experiments" / "t1_scope"))
import g_mom_verify as gv  # noqa: E402  (N21's concordance of two fits, unchanged)

STAGES = {"E8.25_late": "T2/E8.25_late.h5ad", "E8.5": "T2/E8.5.h5ad", "E8.75": "T2/E8.75.h5ad", "E9.5": "T2/E9.5.h5ad"}
TRIPLES = (("E8.25_late", "E8.5", "E8.75"), ("E8.5", "E8.75", "E9.5"))
PANEL_BOARD = "T2:heart:val_extrap"  # the heart backtests' panel
N_PERM = 2000
PERM_SEED = 20261113
ALPHA = 0.05
NEVER = "E10.5"  # the heart backtests' target: this script must not open it


def realised_on(pred_cf: ts.CategoryFactors, real_cf: ts.CategoryFactors) -> np.ndarray:
    """The realised next-step vectors carried onto the predicted fit's metacells: per predicted metacell, the mean over
    its cells of the vector of each cell's realised metacell. Both fits hold the same cells of stage B."""
    if set(pred_cf.index) != set(real_cf.index):
        raise ValueError(f"{pred_cf.category}: the two fits hold different cells of the middle stage")
    out = np.zeros_like(pred_cf.logf)
    n = np.zeros(pred_cf.n_metacells)
    for name, m in pred_cf.index.items():
        out[m] += real_cf.logf[real_cf.index[name]]
        n[m] += 1
    if not np.array_equal(n.astype(np.int64), pred_cf.sizes):
        raise ValueError(f"{pred_cf.category}: metacell sizes do not match the index")
    return out / n[:, None]


def centred(F: np.ndarray, w: np.ndarray) -> np.ndarray:
    """F with its metacell-weighted mean removed per gene (w sums to 1)."""
    return F - (w[:, None] * F).sum(0)


def weighted_corr(P: np.ndarray, R: np.ndarray, w: np.ndarray) -> float:
    """Pearson correlation over (metacell, gene) entries, metacell weights w (sum 1), both centred per gene; NaN when
    either side has no variance."""
    P, R = centred(P, w), centred(R, w)
    cov = float((w[:, None] * P * R).sum())
    vp, vr = float((w[:, None] * P * P).sum()), float((w[:, None] * R * R).sum())
    return cov / np.sqrt(vp * vr) if vp > 0 and vr > 0 else float("nan")


def triple_reading(a, b, c, names: tuple[str, str, str], pair_index: int) -> dict:
    finished, raw = harness.t2_gmom_pair_factors(a, b)
    s, e = harness._shared_genes(b, c)
    cats = sorted(finished.categories)
    real = ts.fit_state_factors(s, e, "earlier", cats, key=harness.T2_LINEAGE_KEY, n_hvg=e.n_vars)
    if list(real.genes) != list(finished.genes):
        raise ValueError("the two fits hold different genes")
    unfinished = ts.finish_momentum(raw, 1.0)  # reported only: its zeroing and clipping still apply at damping 1
    rows, mats = {}, {}
    for cat in cats:
        if cat not in real.categories:
            rows[cat] = {"skipped": real.skipped.get(cat, "no realised fit")}
            continue
        pc = finished.categories[cat]
        w = pc.sizes / pc.sizes.sum()
        R = realised_on(pc, real.categories[cat])
        mats[cat] = (pc.logf, R, w, int(pc.sizes.sum()))
        raw_cf = raw.categories[cat]
        rows[cat] = {"n_cells": int(pc.sizes.sum()), "n_metacells": int(pc.n_metacells),
                     "corr_finished": weighted_corr(pc.logf, R, w),
                     "corr_before_damping": weighted_corr(unfinished.categories[cat].logf, R, w),
                     "corr_raw_mean_free": weighted_corr(raw_cf.logf, R, raw_cf.sizes / raw_cf.sizes.sum()),
                     "share_entries_nonzero": float((pc.logf != 0).mean())}
    use = [cat for cat in mats if np.isfinite(rows[cat]["corr_finished"])]
    cells = np.array([mats[cat][3] for cat in use], dtype=float)

    def stat(perm_rng=None) -> float:
        vals = []
        for cat in use:
            P, R, w, _ = mats[cat]
            if perm_rng is not None:
                P = P[perm_rng.permutation(P.shape[0])]
            vals.append(weighted_corr(P, R, w))
        vals = np.array(vals)
        ok = np.isfinite(vals)
        return float((cells[ok] * vals[ok]).sum() / cells[ok].sum()) if ok.any() else float("nan")

    s_obs = stat()
    rng = np.random.default_rng([PERM_SEED, pair_index])
    null = np.array([stat(rng) for _ in range(N_PERM)])
    null = null[np.isfinite(null)]
    p = float((np.sum(null >= s_obs) + 1) / (len(null) + 1)) if np.isfinite(s_obs) and len(null) else float("nan")
    return {"triple": list(names), "lineages": rows, "lineages_used": use, "skipped_predicted": dict(finished.skipped),
            "S": s_obs, "p_upper": p, "null_mean": float(null.mean()) if len(null) else None,
            "null_sd": float(null.std(ddof=1)) if len(null) > 1 else None, "null_q95": float(np.quantile(null, 0.95))
            if len(null) else None, "passes": bool(np.isfinite(s_obs) and s_obs > 0 and p <= ALPHA)}, finished


def run(data: Path) -> dict:
    if any(NEVER in f for f in STAGES.values()):
        raise ValueError("N26 must not open the heart backtests' target")
    genes = load_boards()[PANEL_BOARD].genes
    stage = {k: load_stage(data / f, genes=genes, drop_duplicates=True) for k, f in STAGES.items()}
    out, fits = {"triples": []}, []
    for i, (a, b, c) in enumerate(TRIPLES):
        res, finished = triple_reading(stage[a], stage[b], stage[c], (a, b, c), i)
        out["triples"].append(res)
        fits.append(finished)
    # reported only: N21's concordance of the two consecutive fits (E8.25_late -> E8.5 against E8.5 -> E8.75), and of
    # the second with the third (E8.5 -> E8.75 against E8.75 -> E9.5)
    third, _ = harness.t2_gmom_pair_factors(stage["E8.75"], stage["E9.5"])
    conc = {}
    for label, (x, y) in {"E8.25_late->E8.5 vs E8.5->E8.75": (fits[0], fits[1]),
                          "E8.5->E8.75 vs E8.75->E9.5": (fits[1], third)}.items():
        both = sorted(set(x.categories) & set(y.categories))
        conc[label] = {cat: gv.concordance(x, y, cat, set()) for cat in both}
    out["consecutive_fit_concordance"] = conc
    out["gate"] = {"rule": f"S > 0 with p <= {ALPHA} on both triples", "passes": bool(all(t["passes"] for t in
                                                                                         out["triples"]))}
    out["settings"] = {"n_perm": N_PERM, "perm_seed": PERM_SEED, "alpha": ALPHA, "damping": harness.B0_DAMPING,
                       "panel_board": PANEL_BOARD, "stages": STAGES, "zero_below": ts.ZERO_BELOW,
                       "clip_at": ts.CLIP_AT}
    return out


def text(res: dict) -> str:
    lines = [f"N26: gate {'PASSES' if res['gate']['passes'] else 'FAILS'} ({res['gate']['rule']}); commit "
             f"{res['commit'][:7]}, tree clean {res['tree_clean']}"]
    for t in res["triples"]:
        a, b, c = t["triple"]
        lines.append(f"  {a} -> {b} predicts {b} -> {c}: S {t['S']:+.3f}, permutation null mean {t['null_mean']:+.3f} "
                     f"sd {t['null_sd']:.3f} (95th {t['null_q95']:+.3f}), p {t['p_upper']:.4f}: "
                     f"{'passes' if t['passes'] else 'fails'}")
        for cat, r in t["lineages"].items():
            if "skipped" in r:
                lines.append(f"      {cat:<28} skipped: {r['skipped']}")
            else:
                lines.append(f"      {cat:<28} cells {r['n_cells']:>6}, metacells {r['n_metacells']:>3}, corr "
                             f"{r['corr_finished']:+.3f} (before damping, zeroing and clipping "
                             f"{r['corr_before_damping']:+.3f}; raw {r['corr_raw_mean_free']:+.3f}), entries moved "
                             f"{r['share_entries_nonzero']:.2f}")
    for label, rows in res["consecutive_fit_concordance"].items():
        ax = [r["axis_abs_cosine"] for r in rows.values() if np.isfinite(r["axis_abs_cosine"])]
        sp = [r["size_spearman"] for r in rows.values() if np.isfinite(r["size_spearman"])]
        lines.append(f"  reported only, N21's concordance, {label}: leading-axis |cosine| median "
                     f"{np.median(ax):.2f} over {len(ax)} lineages; per-gene size Spearman median {np.median(sp):.2f}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    commit = harness._git("rev-parse", "HEAD")
    clean = harness._git("status", "--porcelain", "--untracked-files=no") == ""
    res = run(Path(a.data))
    res.update(commit=commit, tree_clean=clean)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "n26.json").write_text(json.dumps(res, indent=1))
    (out / "n26.txt").write_text(text(res))
    print(text(res), end="")


if __name__ == "__main__":
    main()
