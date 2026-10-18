"""Part 2b: fit the physical section pitch p (in xy units) by an isotropy collapse (DGX).

For stacked pairs k sections apart with in-plane offset dxy, an isotropic tissue predicts a dissimilarity equal to the
in-plane curve at r = sqrt(dxy^2 + (k p)^2). We bin stacked pairs by (k, dxy) for k = 1..3 and dxy < 40, read the
in-plane curve (pairs in the same section, d >= 4, isotonic fit) at r(p), and choose p on a grid by weighted least
squares. s = step / p is the stretch of the stored z relative to x/y (s = 1: z already in xy units).

Also writes the k = 1 dissimilarity profile against dxy, which separates a small physical pitch (steep profile at
small dxy) from a large one (flat profile while dxy << pitch).

Writes runs/coords/zpitch_fit.csv, zpitch_profile.csv. Usage: uv run python experiments/coords/zpitch_fit.py
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
from coordlib import FILES, coords_of, snap_z  # noqa: E402

from vec.prep import load_stage  # noqa: E402

R_IN = 150.0
BIN_IN = 2.0
DXY_MAX = 40.0
BIN_Z = 2.0
KS = (1, 2, 3)
P_GRID = np.round(np.arange(3.0, 40.01, 0.25), 2)


def measures(P, lab_c, lab_f, i, j):
    return {"gamma_pc": 0.5 * ((P[i] - P[j]) ** 2).sum(1), "mism_group": (lab_c[i] != lab_c[j]).astype(float),
            "mism_label": (lab_f[i] != lab_f[j]).astype(float)}


def main(out: Path, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    fit_rows, prof_rows = [], []
    for stage, path in FILES.items():
        a = load_stage(path)
        C = coords_of(a)
        lev, levels, _ = snap_z(C[:, 2])
        step = float(np.median(np.diff(levels)))
        P = PCA(n_components=20, random_state=0).fit_transform(a.X.toarray())
        lab_c = np.asarray(a.obs["group"]).astype(str)
        lab_f = np.asarray(a.obs["celltype"]).astype(str)
        trees = {g: (np.flatnonzero(lev == g), cKDTree(C[lev == g, :2])) for g in np.unique(lev)}

        nb_in = int(R_IN / BIN_IN)
        acc_in = {m: np.zeros(nb_in) for m in ("gamma_pc", "mism_group", "mism_label")}
        cnt_in, dsum_in = np.zeros(nb_in), np.zeros(nb_in)
        nb_z = int(DXY_MAX / BIN_Z)
        acc_z = {(k, m): np.zeros(nb_z) for k in KS for m in acc_in}
        cnt_z = {k: np.zeros(nb_z) for k in KS}
        d2sum_z = {k: np.zeros(nb_z) for k in KS}
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
                if g + k not in trees:
                    continue
                idx2, tree2 = trees[g + k]
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
        ok_in = (cnt_in >= 200) & (dmean_in >= 4.0)
        for m in acc_in:
            g_in = acc_in[m] / np.maximum(cnt_in, 1)
            dgrid = dmean_in[ok_in]
            gfit = IsotonicRegression(increasing=True).fit_transform(dgrid, g_in[ok_in])
            obs, w, dxy2, kk = [], [], [], []
            for k in KS:
                okz = cnt_z[k] >= 50
                obs.append(acc_z[(k, m)][okz] / cnt_z[k][okz])
                w.append(cnt_z[k][okz])
                dxy2.append(d2sum_z[k][okz] / cnt_z[k][okz])
                kk.append(np.full(okz.sum(), k))
                for bi in np.flatnonzero(okz):
                    prof_rows.append({"stage": stage, "measure": m, "k": k, "dxy_lo": bi * BIN_Z,
                                      "dxy_rms": float(np.sqrt(d2sum_z[k][bi] / cnt_z[k][bi])),
                                      "n_pairs": int(cnt_z[k][bi]), "value": float(acc_z[(k, m)][bi] / cnt_z[k][bi])})
            obs, w, dxy2, kk = map(np.concatenate, (obs, w, dxy2, kk))
            sse = []
            for p in P_GRID:
                r = np.sqrt(dxy2 + (kk * p) ** 2)
                pred = np.interp(r, dgrid, gfit)  # clamps outside the fitted range
                sse.append(float((w * (obs - pred) ** 2).sum() / w.sum()))
            sse = np.array(sse)
            best = int(np.argmin(sse))
            # rough 'flat region': grid points within 10% of the minimum excess error over the sill-variance scale
            tol = sse[best] + 0.1 * (np.median(sse) - sse[best])
            near = P_GRID[sse <= tol]
            pb = float(P_GRID[best])
            fit_rows.append({"stage": stage, "measure": m, "step": step, "p_best": pb,
                             "p_lo": float(near.min()), "p_hi": float(near.max()), "s_best": step / pb,
                             "s_lo": step / float(near.max()), "s_hi": step / float(near.min()),
                             "rmse_best": float(np.sqrt(sse[best])), "rmse_at_p_eq_step":
                             float(np.sqrt(sse[int(np.argmin(np.abs(P_GRID - step)))])),
                             "inplane_floor": float(gfit[0]), "inplane_sill": float(gfit[-1])})
        f = pd.DataFrame([r for r in fit_rows if r["stage"] == stage])
        print(stage, "step", step, f[["measure", "p_best", "p_lo", "p_hi", "s_best", "rmse_best",
                                       "rmse_at_p_eq_step"]].round(3).to_dict("records"), flush=True)
    pd.DataFrame(fit_rows).to_csv(out / "zpitch_fit.csv", index=False)
    pd.DataFrame(prof_rows).to_csv(out / "zpitch_profile.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    main(Path(ap.parse_args().out))
