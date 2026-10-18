"""Part 2: is z in the same units as x/y? Variogram matching between in-plane and across-section pairs (DGX).

Idea: if tissue is statistically isotropic at small scales, two cells a physical distance r apart are equally
dissimilar whichever direction separates them. In-plane pairs give the dissimilarity curve gamma_xy(d) against a
known xy distance d. Pairs stacked k sections apart (xy offset < eps) give gamma_z(k). Reading gamma_z(k) off
gamma_xy gives the physical vertical separation h_k; pitch_k = sqrt(h_k^2 - mean dxy^2) / k is the section pitch
in xy units, and s = step / pitch is how much the file's z is stretched relative to x/y (s = 1: same units).
The lead's decision of 2026-10-19: the z-unit factor is unresolved and is never used by any procedure; this script is
a characterisation.

Dissimilarity measures: (a) semivariance of 20 expression PCs (log1p CP10k), (b) coarse-group mismatch rate,
(c) fine-label mismatch rate. Known biases: section misregistration inflates gamma_z (pitch too large, s < 1);
anatomy elongated along z deflates gamma_z (pitch too small, s > 1).

Revision 2026-10-19: duplicated records are dropped first; levels k steps apart are paired by z difference (not index),
so a stack gap is not one step; levels holding two physical sections are left out; each stage draws from its own
generator seeded by its name.

Writes runs/coords/variogram_xy.csv, variogram_z.csv, zscale_estimates.csv.
Usage: uv run python experiments/coords/zscale.py [--out runs/coords] [--stage NAME=PATH ...]
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

R_MAX = 120.0
BIN = 4.0
EPS = 6.0
K_MAX = 10
N_ANCHOR_PER_LEVEL = 400


def pair_stats(P, lab_c, lab_f, i, j):
    dp = 0.5 * ((P[i] - P[j]) ** 2).sum(1)
    return dp, (lab_c[i] != lab_c[j]).astype(float), (lab_f[i] != lab_f[j]).astype(float)


def main(out_dir: Path, files: dict[str, str], seed: int = 0) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    xy_rows, z_rows, est_rows = [], [], []
    for stage, path in files.items():
        rng = stage_rng(stage, seed)
        a, dd = load(path, drop_duplicates=True)
        C = coords_of(a)
        lev, levels, _ = snap_z(C[:, 2])
        step = level_step(levels)
        skip = set(multi_section_levels(a, lev))
        X = a.X.toarray()
        P = PCA(n_components=20, random_state=0).fit_transform(X)
        lab_c = np.asarray(a.obs["group"]).astype(str)
        lab_f = np.asarray(a.obs["celltype"]).astype(str)
        trees = {g: (np.flatnonzero(lev == g), cKDTree(C[lev == g, :2])) for g in np.unique(lev) if g not in skip}

        # in-plane pairs
        nb = int(R_MAX / BIN)
        acc = np.zeros((nb, 5))  # count, sum d, sum gamma_pc, sum mism_c, sum mism_f
        for idx, tree in trees.values():
            if idx.size < 20:
                continue
            anc = rng.choice(idx.size, min(N_ANCHOR_PER_LEVEL, idx.size), replace=False)
            nbrs = tree.query_ball_point(C[idx[anc], :2], r=R_MAX)
            ii = np.repeat(idx[anc], [len(n) for n in nbrs])
            jj = idx[np.concatenate([np.asarray(n, dtype=int) for n in nbrs])]
            keep = ii != jj
            ii, jj = ii[keep], jj[keep]
            d = np.linalg.norm(C[ii, :2] - C[jj, :2], axis=1)
            gp, mc, mf = pair_stats(P, lab_c, lab_f, ii, jj)
            b = np.minimum((d / BIN).astype(int), nb - 1)
            for col, v in enumerate([np.ones_like(d), d, gp, mc, mf]):
                acc[:, col] += np.bincount(b, weights=v, minlength=nb)
        cnt = np.maximum(acc[:, 0], 1)
        xy = pd.DataFrame({"stage": stage, "bin_lo": np.arange(nb) * BIN, "n_pairs": acc[:, 0],
                           "d_mean": acc[:, 1] / cnt, "gamma_pc": acc[:, 2] / cnt, "mism_group": acc[:, 3] / cnt,
                           "mism_label": acc[:, 4] / cnt})
        xy_rows.append(xy)

        # stacked pairs k steps apart, levels paired by z difference
        for k in range(1, K_MAX + 1):
            partner = level_partners(levels, k, step)
            tot = np.zeros(5)
            n_level_pairs = 0
            for g, (idx, _) in trees.items():
                h = int(partner[g])
                if h < 0 or h not in trees:
                    continue
                idx2, tree2 = trees[h]
                n_level_pairs += 1
                nbrs = tree2.query_ball_point(C[idx, :2], r=EPS)
                ii = np.repeat(idx, [len(n) for n in nbrs])
                if ii.size == 0:
                    continue
                jj = idx2[np.concatenate([np.asarray(n, dtype=int) for n in nbrs])]
                d2 = ((C[ii, :2] - C[jj, :2]) ** 2).sum(1)
                gp, mc, mf = pair_stats(P, lab_c, lab_f, ii, jj)
                tot += [ii.size, d2.sum(), gp.sum(), mc.sum(), mf.sum()]
            n = max(tot[0], 1)
            z_rows.append({"stage": stage, "k": k, "dz_raw": k * step, "n_level_pairs": n_level_pairs,
                           "n_pairs": tot[0], "mean_dxy2": tot[1] / n, "gamma_pc": tot[2] / n,
                           "mism_group": tot[3] / n, "mism_label": tot[4] / n})

        # invert gamma_xy at each gamma_z(k)
        zdf = pd.DataFrame([r for r in z_rows if r["stage"] == stage])
        ok = xy["n_pairs"] >= 200
        for meas in ["gamma_pc", "mism_group", "mism_label"]:
            dgrid = xy.loc[ok, "d_mean"].to_numpy()
            gfit = IsotonicRegression(increasing=True).fit_transform(dgrid, xy.loc[ok, meas].to_numpy())
            for _, r in zdf.iterrows():
                gz = r[meas]
                if r["n_pairs"] <= 0 or gz <= gfit[0] or gz >= gfit[-1] or np.ptp(gfit) <= 0:
                    h = np.nan
                else:
                    # first crossing of the monotone curve
                    dstar = float(np.interp(gz, gfit, dgrid))
                    h = np.sqrt(max(dstar ** 2 - r["mean_dxy2"], 0.0))
                pitch = h / r["k"] if np.isfinite(h) else np.nan
                est_rows.append({"stage": stage, "measure": meas, "k": int(r["k"]), "step": step,
                                 "gamma_z": gz, "h_xy_units": h, "pitch_xy_units": pitch,
                                 "z_stretch_s": step / pitch if pitch and np.isfinite(pitch) and pitch > 0
                                 else np.nan,
                                 "sill_fraction": (gz - gfit[0]) / max(gfit[-1] - gfit[0], 1e-12)})
        e = pd.DataFrame([r for r in est_rows if r["stage"] == stage])
        summ = e[(e["k"] <= 4)].groupby("measure")["z_stretch_s"].median().round(3).to_dict()
        print(f"{stage}: step {step:g}; dropped {dd['n_dropped']} duplicates; levels left out "
              f"{sorted(float(levels[g]) for g in skip)}; median s over k<=4: {summ}", flush=True)

    pd.concat(xy_rows).to_csv(out_dir / "variogram_xy.csv", index=False)
    pd.DataFrame(z_rows).to_csv(out_dir / "variogram_z.csv", index=False)
    pd.DataFrame(est_rows).to_csv(out_dir / "zscale_estimates.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    ap.add_argument("--stage", action="append", help="NAME=PATH (repeatable; default: every released file)")
    args = ap.parse_args()
    main(Path(args.out), stage_specs(args.stage, FILES))
