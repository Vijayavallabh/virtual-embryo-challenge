"""Reading N25 (experiments/harness/TRIES.md): a scratch reading the AI assistant ran on 2026-11-13 BEFORE writing any
row, kept here as it was run (only this paragraph is new). It compares predictors with fold targets, so rule 4 counts
it: one try each on T3, the heart and the embryo. Its output is in results/n25_revert_reading/. The registered
follow-up is try K7 (k7_wt_revert.py).

Usage (repo root): uv run --frozen --no-sync python -I -B experiments/t3_decomp/n25_revert_reading.py OUT_DIR

Released and unlocked files only; no scorer, no build, no test board: does a reference
sample's own deviation from the other released samples predict the per-gene change from it to another sample?

Per sample and lineage group (>= MIN cells): y = log(mean within-cell panel proportion + EPS), centred by its median
over genes. For an ordered pair (r -> t) and an allowed set A (WT samples, never r or t):
  centroid   c_L = mean over A of y_{s,L}                       pred_L = c_L - y_{r,L}   ("revert to the centroid")
  trend      a_L + b_L * day fitted on A; value at r's day       pred_L = trend_L(day_r) - y_{r,L} (+ b_L * gap)
  global     d_r = median over groups of (y_{r,L} - c_L)         pred_L = -d_r            (one profile for all groups)
  last step  0.5 * (y_{r,L} - y_{prev(r),L})                     (B0's damped momentum, in this space)
Against the true change y_{t,L} - y_{r,L}: Spearman, least-squares multiplier, and squared error relative to "no
change" at multiplier 1 and 0.5. Printed per pair as medians over groups.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

from vec import t2_build
from vec.boards import REPO_ROOT
from vec.prep import load_stage

D = REPO_ROOT / "data" / "raw"
WT = ("E6.75", "E7.25", "E7.5", "E8.0", "E8.25_late", "E8.5", "E8.75", "E9.5", "E10.5")
DAY = {"E6.75": 6.75, "E7.25": 7.25, "E7.5": 7.5, "E8.0": 8.0, "E8.25_late": 8.3, "E8.5": 8.5, "E8.75": 8.75,
       "E9.5": 9.5, "E10.5": 10.5}
KO = {"E8.75_gata4_ko_1": "E8.75", "E8.75_gata4_ko_2": "E8.75", "E9.5_mab21l2_ko": "E9.5"}
NON_MESP1 = ("Endoderm", "Ectoderm", "Neural", "NCC")
EPS, MIN = 1e-5, 200
OUT = Path(sys.argv[1])

Y, N = {}, {}
genes = None
raw = {}
for sub, names in (("T2", WT), ("T3", tuple(KO))):
    for s in names:
        a = load_stage(D / sub / f"{s}.h5ad", drop_duplicates=True)
        raw[s] = a
        genes = set(a.var_names) if genes is None else genes & set(a.var_names)
genes = sorted(genes)
for s, a in raw.items():
    a = a[:, genes]
    _, C, lib, _ = t2_build._stage_counts(a, genes)
    C = sparse.csr_matrix(C, dtype=np.float64)
    grp = np.asarray(a.obs["group"]).astype(str)
    tot = np.asarray(C.sum(1)).ravel()
    P = sparse.diags(1.0 / np.maximum(tot, 1.0)) @ C
    Y[s], N[s] = {}, {}
    for g in sorted(set(grp)):
        rows = np.flatnonzero(grp == g)
        if g == "Unknown" or len(rows) < MIN:
            continue
        y = np.log(np.asarray(P[rows].mean(0)).ravel() + EPS)
        Y[s][g] = y - np.median(y)
        N[s][g] = int(len(rows))
    print(f"{s}: {a.n_obs} cells, median panel counts {np.median(tot):.0f}, groups {len(Y[s])}", flush=True)
del raw


def stats(pred, truth):
    rho = float(spearmanr(pred, truth).statistic)
    k = float(pred @ truth / max(pred @ pred, 1e-12))
    base = float(truth @ truth)
    return {"rho": rho, "mult": k, "rel1": float(((truth - pred) ** 2).sum() / base),
            "rel05": float(((truth - 0.5 * pred) ** 2).sum() / base), "sd_truth": float(truth.std()),
            "sd_pred": float(pred.std())}


def predictors(r, t, allowed):
    """Per group shared by r and t: the predictions of the change r -> t from samples in `allowed` (never r or t)."""
    A = [s for s in allowed if s not in (r, t)]
    out = {}
    day_r = DAY[r]
    gap = (DAY[t] if t in DAY else DAY[KO[t]]) - day_r
    cent = {}
    for g in Y[r]:
        have = [s for s in A if g in Y[s]]
        if len(have) >= 2:
            cent[g] = np.mean([Y[s][g] for s in have], axis=0)
    d_glob = np.median([Y[r][g] - c for g, c in cent.items()], axis=0) if len(cent) >= 3 else None
    prev = [s for s in WT if DAY[s] < day_r]
    prev = prev[-1] if prev else None
    for g in Y[r]:
        if g not in Y[t]:
            continue
        truth = Y[t][g] - Y[r][g]
        rec = {"n_r": N[r][g], "n_t": N[t][g]}
        if g in cent:
            rec["centroid"] = stats(cent[g] - Y[r][g], truth)
            have = [s for s in A if g in Y[s]]
            rec["k"] = len(have)
            if len(have) >= 3:
                x = np.array([DAY[s] for s in have])
                M = np.stack([Y[s][g] for s in have])
                xm = x.mean()
                b = ((x - xm)[:, None] * (M - M.mean(0))).sum(0) / max(((x - xm) ** 2).sum(), 1e-12)
                a0 = M.mean(0) + b * (day_r - xm)
                rec["trend"] = stats(a0 - Y[r][g], truth)
                rec["trend+slope"] = stats(a0 - Y[r][g] + b * gap, truth)
                rec["slope_only"] = stats(b * gap, truth) if abs(gap) > 0 else None
        if d_glob is not None:
            rec["global"] = stats(-d_glob, truth)
        if prev is not None and prev != t and g in Y[prev]:
            rec["last_step"] = stats(0.5 * (Y[r][g] - Y[prev][g]), truth)
        out[g] = rec
    return out


def show(title, r, t, allowed, groups=None):
    res = predictors(r, t, allowed)
    if groups is not None:
        res = {g: v for g, v in res.items() if g in groups}
    line = f"{title:<46} groups {len(res):2d}"
    summ = {}
    for name in ("centroid", "trend", "trend+slope", "slope_only", "global", "last_step"):
        vals = [v[name] for v in res.values() if v.get(name)]
        if not vals:
            continue
        summ[name] = {k: float(np.median([x[k] for x in vals])) for k in ("rho", "mult", "rel1", "rel05")}
        line += f" | {name}: rho {summ[name]['rho']:+.2f} mult {summ[name]['mult']:+.2f} rel {summ[name]['rel1']:.2f}/{summ[name]['rel05']:.2f}"
    print(line, flush=True)
    return {"pair": [r, t], "allowed": list(allowed), "groups": res, "summary": summ}


ALL = {}
print("\n== T3 folds: the matched wild type as reference, a knockout as the target; allowed = every other WT sample ==")
for ko, r in KO.items():
    ALL[f"{r}->{ko}|nonmesp1"] = show(f"{r} -> {ko} (non-Mesp1 groups)", r, ko, WT, NON_MESP1)
    ALL[f"{r}->{ko}|all"] = show(f"{r} -> {ko} (all groups)", r, ko, WT)
    late = tuple(s for s in WT if DAY[s] >= 8.0)
    ALL[f"{r}->{ko}|nonmesp1|late"] = show(f"{r} -> {ko} (non-Mesp1; WT from E8.0 on)", r, ko, late, NON_MESP1)
print("\n== knockout against knockout: replicate and cross-knockout agreement of the per-gene change from the WT ==")
for a_, b_ in (("E8.75_gata4_ko_1", "E8.75_gata4_ko_2"),):
    for g in NON_MESP1:
        if g in Y[a_] and g in Y[b_] and g in Y["E8.75"]:
            fa, fb = Y[a_][g] - Y["E8.75"][g], Y[b_][g] - Y["E8.75"][g]
            print(f"  {g}: Spearman of the two Gata4 replicates' changes {spearmanr(fa, fb).statistic:+.2f}; SD {fa.std():.2f}, {fb.std():.2f}")
print("\n== WT pairs, next stage as target (extrapolation-like); allowed = earlier WT samples only, then all other WT ==")
for i in range(1, len(WT) - 1):
    r, t = WT[i], WT[i + 1]
    ALL[f"{r}->{t}|earlier"] = show(f"{r} -> {t} (allowed: earlier)", r, t, WT[:i])
    ALL[f"{r}->{t}|all"] = show(f"{r} -> {t} (allowed: all other)", r, t, WT)
print("\n== WT pairs two stages ahead ==")
for i in range(1, len(WT) - 2):
    r, t = WT[i], WT[i + 2]
    ALL[f"{r}->{t}|earlier"] = show(f"{r} -> {t} (allowed: earlier)", r, t, WT[:i])
print("\n== the test pair's reference: E10.5's deviation from the earlier samples (no target; sizes only) ==")
r = "E10.5"
for g in Y[r]:
    have = [s for s in WT[:-1] if g in Y[s]]
    if len(have) >= 3:
        c = np.mean([Y[s][g] for s in have], axis=0)
        print(f"  {g:<20} samples {len(have)}  SD of (E10.5 - centroid) {float((Y[r][g] - c).std()):.2f}")
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "revert_diag.json").write_text(json.dumps(ALL, indent=1))
