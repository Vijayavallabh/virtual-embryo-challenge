"""Part 2b: fit the physical section pitch p (in xy units) by an isotropy collapse (DGX).

For stacked pairs k sections apart with in-plane offset dxy, an isotropic tissue predicts a dissimilarity equal to the
in-plane curve at r = sqrt(dxy^2 + (k p)^2). We bin stacked pairs by (k, dxy) for dxy < 40, read the in-plane curve
(pairs in the same section, d >= 4, isotonic fit) at r(p), and choose p on a grid by weighted least squares.
s = step / p is the stretch of the stored z relative to x/y (s = 1: z already in xy units). The lead's decision of
2026-10-19: the z-unit factor is unresolved and never used by any procedure; this is a characterisation.

Revision 2026-10-19 (independent review):
- duplicated records are dropped first, levels k steps apart are paired by z difference (the KO's 30-unit gap is three
  steps, not one), levels holding two physical sections are left out, and each stage has its own generator;
- the flat range is reported as open-ended where it reaches the grid floor or the part of the grid where the
  smallest-offset k = 1 pairs read the in-plane curve below its first distance (np.interp clamps there, so the fit
  cannot tell those pitches apart): p_lo_open / p_hi_open, p_clamp, and the text columns range_p / range_s;
- the k set changes the answer (review item 5), so three sets are fitted: k = 1-3 (the main fit), k = 1, k = 1-6.

Also writes the stacked dissimilarity profile against dxy for every k, which separates a small physical pitch (steep
profile at small dxy) from a large one (flat profile while dxy << pitch).

Writes runs/coords/zpitch_fit.csv, zpitch_profile.csv.
Usage: uv run python experiments/coords/zpitch_fit.py [--out runs/coords] [--stage NAME=PATH ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.isotonic import IsotonicRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coordlib import (  # noqa: E402
    FILES,
    coords_of,
    level_partners,
    level_step,
    load,
    multi_section_levels,
    snap_z,
    stage_rng,
    stage_specs,
)

R_IN = 150.0
BIN_IN = 2.0
DXY_MAX = 40.0
BIN_Z = 2.0
D_MIN_IN = 4.0
KS = (1, 2, 3, 4, 5, 6)
K_SETS = {"1-3": (1, 2, 3), "1": (1,), "1-6": (1, 2, 3, 4, 5, 6)}
P_GRID = np.round(np.arange(3.0, 40.01, 0.25), 2)
MEASURES = ("gamma_pc", "mism_group", "mism_label")


def measures(P, lab_c, lab_f, i, j):
    return {"gamma_pc": 0.5 * ((P[i] - P[j]) ** 2).sum(1), "mism_group": (lab_c[i] != lab_c[j]).astype(float),
            "mism_label": (lab_f[i] != lab_f[j]).astype(float)}


def fit_pitch(dgrid: np.ndarray, gfit: np.ndarray, obs: np.ndarray, w: np.ndarray, dxy2: np.ndarray, kk: np.ndarray,
              grid: np.ndarray = P_GRID) -> dict:
    """Weighted least-squares pitch on `grid`, with the flat range and whether its ends are identified.

    The flat range is the set of grid points within 10% of the way from the minimum SSE to the median SSE. Its lower
    end is open when it reaches the grid floor or p_clamp, the pitch below which the smallest-offset pairs of the
    smallest k read the in-plane curve below its first distance (clamped). Its upper end is open at the grid top."""
    sse = np.array([float((w * (obs - np.interp(np.sqrt(dxy2 + (kk * p) ** 2), dgrid, gfit)) ** 2).sum() / w.sum())
                    for p in grid])
    best = int(np.argmin(sse))
    tol = sse[best] + 0.1 * (np.median(sse) - sse[best])
    near = grid[sse <= tol]
    kmin = kk == kk.min()
    p_clamp = float(np.sqrt(max(dgrid[0] ** 2 - float(dxy2[kmin].min()), 0.0)) / kk.min())
    lo, hi = float(near.min()), float(near.max())
    lo_open = bool(lo <= grid[0] or lo < p_clamp)
    hi_open = bool(hi >= grid[-1])
    p_id = max(float(grid[0]), p_clamp)
    lo_txt = f"(open, <={p_id:.2f}" if lo_open else f"[{lo:g}"
    hi_txt = f">={grid[-1]:g}, open)" if hi_open else f"{hi:g}]"
    return {"p_best": float(grid[best]), "p_lo": lo, "p_hi": hi, "p_lo_open": lo_open, "p_hi_open": hi_open,
            "p_clamp": p_clamp, "p_identified_min": p_id, "p_best_in_clamped_zone": bool(grid[best] < p_clamp),
            "range_p": f"{lo_txt}, {hi_txt}", "rmse_best": float(np.sqrt(sse[best])), "sse": sse}


def main(out: Path, files: dict[str, str], seed: int = 0) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fit_rows, prof_rows = [], []
    for stage, path in files.items():
        rng = stage_rng(stage, seed)
        a, dd = load(path, drop_duplicates=True)
        C = coords_of(a)
        lev, levels, _ = snap_z(C[:, 2])
        step = level_step(levels)
        skip = set(multi_section_levels(a, lev))
        P = PCA(n_components=20, random_state=0).fit_transform(a.X.toarray())
        lab_c = np.asarray(a.obs["group"]).astype(str)
        lab_f = np.asarray(a.obs["celltype"]).astype(str)
        trees = {g: (np.flatnonzero(lev == g), cKDTree(C[lev == g, :2])) for g in np.unique(lev) if g not in skip}
        partners = {k: level_partners(levels, k, step) for k in KS}

        nb_in = int(R_IN / BIN_IN)
        acc_in = {m: np.zeros(nb_in) for m in MEASURES}
        cnt_in, dsum_in = np.zeros(nb_in), np.zeros(nb_in)
        nb_z = int(DXY_MAX / BIN_Z)
        acc_z = {(k, m): np.zeros(nb_z) for k in KS for m in MEASURES}
        cnt_z = {k: np.zeros(nb_z) for k in KS}
        d2sum_z = {k: np.zeros(nb_z) for k in KS}
        n_level_pairs = dict.fromkeys(KS, 0)
        for g, (idx, tree) in trees.items():
            if idx.size < 20:
                continue
            anc = idx[rng.choice(idx.size, min(500, idx.size), replace=False)]
            nbrs = tree.query_ball_point(C[anc, :2], r=R_IN)
            ii = np.repeat(anc, [len(n) for n in nbrs])
            jj = idx[np.concatenate([np.asarray(n, dtype=int) for n in nbrs])]
            keep = ii != jj
            ii, jj = ii[keep], jj[keep]
            d = np.linalg.norm(C[ii, :2] - C[jj, :2], axis=1)
            b = np.minimum((d / BIN_IN).astype(int), nb_in - 1)
            cnt_in += np.bincount(b, minlength=nb_in)
            dsum_in += np.bincount(b, weights=d, minlength=nb_in)
            for m, v in measures(P, lab_c, lab_f, ii, jj).items():
                acc_in[m] += np.bincount(b, weights=v, minlength=nb_in)
            for k in KS:
                h = int(partners[k][g])
                if h < 0 or h not in trees:
                    continue
                idx2, tree2 = trees[h]
                n_level_pairs[k] += 1
                nb2 = tree2.query_ball_point(C[idx, :2], r=DXY_MAX)
                ii = np.repeat(idx, [len(n) for n in nb2])
                if ii.size == 0:
                    continue
                jj = idx2[np.concatenate([np.asarray(n, dtype=int) for n in nb2])]
                d = np.linalg.norm(C[ii, :2] - C[jj, :2], axis=1)
                b = np.minimum((d / BIN_Z).astype(int), nb_z - 1)
                cnt_z[k] += np.bincount(b, minlength=nb_z)
                d2sum_z[k] += np.bincount(b, weights=d ** 2, minlength=nb_z)
                for m, v in measures(P, lab_c, lab_f, ii, jj).items():
                    acc_z[(k, m)] += np.bincount(b, weights=v, minlength=nb_z)

        dmean_in = dsum_in / np.maximum(cnt_in, 1)
        ok_in = (cnt_in >= 200) & (dmean_in >= D_MIN_IN)
        for m in MEASURES:
            g_in = acc_in[m] / np.maximum(cnt_in, 1)
            dgrid = dmean_in[ok_in]
            gfit = IsotonicRegression(increasing=True).fit_transform(dgrid, g_in[ok_in])
            per_k = {}
            for k in KS:
                okz = cnt_z[k] >= 50
                per_k[k] = (acc_z[(k, m)][okz] / cnt_z[k][okz], cnt_z[k][okz], d2sum_z[k][okz] / cnt_z[k][okz],
                            np.full(okz.sum(), k))
                for bi in np.flatnonzero(okz):
                    prof_rows.append({"stage": stage, "measure": m, "k": k, "dxy_lo": bi * BIN_Z,
                                      "dxy_rms": float(np.sqrt(d2sum_z[k][bi] / cnt_z[k][bi])),
                                      "n_pairs": int(cnt_z[k][bi]), "value": float(acc_z[(k, m)][bi] / cnt_z[k][bi])})
            for kname, kset in K_SETS.items():
                obs, w, dxy2, kk = (np.concatenate([per_k[k][i] for k in kset]) for i in range(4))
                if obs.size == 0:
                    continue
                f = fit_pitch(dgrid, gfit, obs, w, dxy2, kk)
                sse = f.pop("sse")
                s_lo = np.nan if f["p_hi_open"] else step / f["p_hi"]
                s_hi = np.nan if f["p_lo_open"] else step / f["p_lo"]
                s_lo_txt = f"(open, <={step / P_GRID[-1]:.2f}" if f["p_hi_open"] else f"[{s_lo:.2f}"
                s_hi_txt = f">={step / f['p_identified_min']:.2f}, open)" if f["p_lo_open"] else f"{s_hi:.2f}]"
                fit_rows.append({"stage": stage, "measure": m, "k_set": kname, "step": step,
                                 "n_dropped_duplicates": dd["n_dropped"],
                                 "levels_left_out": ";".join(f"{levels[g]:g}" for g in sorted(skip)),
                                 "n_level_pairs": ";".join(f"k{k}:{n_level_pairs[k]}" for k in kset), **f,
                                 "s_best": step / f["p_best"], "s_lo": s_lo, "s_hi": s_hi,
                                 "range_s": f"{s_lo_txt}, {s_hi_txt}",
                                 "rmse_at_p_eq_step": float(np.sqrt(sse[int(np.argmin(np.abs(P_GRID - step)))])),
                                 "inplane_d_min": float(dgrid[0]), "inplane_floor": float(gfit[0]),
                                 "inplane_sill": float(gfit[-1])})
        f = pd.DataFrame([r for r in fit_rows if r["stage"] == stage and r["measure"] == "mism_group"])
        print(stage, "step", step, "left out", sorted(skip), f[["k_set", "p_best", "range_p", "range_s", "rmse_best",
                                                                "rmse_at_p_eq_step"]].round(3).to_dict("records"),
              flush=True)
    pd.DataFrame(fit_rows).to_csv(out / "zpitch_fit.csv", index=False)
    pd.DataFrame(prof_rows).to_csv(out / "zpitch_profile.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    ap.add_argument("--stage", action="append", help="NAME=PATH (repeatable; default: every released file)")
    args = ap.parse_args()
    main(Path(args.out), stage_specs(args.stage, FILES))
