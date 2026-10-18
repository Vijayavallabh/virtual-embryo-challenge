"""Part 3: handedness of each file's frame, from anatomical landmarks (DGX).

For each file we build anatomical direction vectors in the raw frame:
  AP  cranial direction: (heart) arterial-pole minus venous-pole cells; (body) anterior minus posterior tissue
  DV  dorsal direction: dorsal tissue (neural tube, dorsal foregut, dorsal paraxial mesoderm) minus core cardiac cells
  LR  left direction: expression-weighted centroid offset of a left-side marker (Pitx2, Nodal, Lefty1, Lefty2)
      within the cells where it is lateralised (mesoderm and cardiac groups), relative to the same cells' centroid
Handedness index h = det[AP, DV, LR] on unit vectors after Gram-Schmidt (AP, DV orthogonalised). In a right-handed
physical frame a normal (situs solitus) embryo gives h = -1: x = left, y = cranial, z = ventral is right-handed, and
det[cranial, dorsal, left] = -1. Opposite signs in two files = mirror images. Bootstrap over cells gives the sign
stability; a permutation test gives a z-score for each marker's offset.

Also: anatomical axes in the scorer's canonical frame (centre, PCA, det +1; full cloud and core), the best proper vs
improper sign flip aligning consecutive stages, and veckit occupancy_dice for (A, B) vs (mirror A, B).

Writes runs/coords/handedness_vectors.csv, handedness_index.csv, canonical_axes.csv, frame_pairs.csv.
Usage: uv run python experiments/coords/handedness.py
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coordlib import FILES, coords_of  # noqa: E402

from vec.prep import load_stage  # noqa: E402

LR_GROUPS = {"CM", "CardiacProgenitor", "PharyngealMesoderm", "LPM", "Pericardium", "EpicardiumST", "Endocardium"}

# Landmark definitions per file: (obs column, labels) for each pole. Chosen from label names before computing any
# geometry; see REPORT.md for the label tables.
LANDMARKS = {
    "E8.0": {
        "heart_arterial": ("cm_celltype", ["aSHF1", "aSHF2"]),
        "heart_venous": ("cm_celltype", ["pSHF"]),
        "body_anterior": ("celltype", ["Forebrain"]),
        "body_posterior": ("celltype", ["Primitive Streak", "Caudal Epiblast", "Allantois"]),
        "dorsal": ("celltype", ["D-FG"]),
        "ventral": ("celltype", ["V-FG"]),
    },
    "E8.25_late": {
        "heart_arterial": ("cm_celltype", ["aSHF1", "aSHF2"]),
        "heart_venous": ("cm_celltype", ["pSHF"]),
        "body_anterior": ("celltype", ["Forebrain"]),
        "body_posterior": ("celltype", ["Caudal Epiblast", "Hindgut"]),
        "dorsal": ("celltype", ["Neural Tube"]),
        "ventral": ("core", [True]),
    },
    "E8.75": {
        "heart_arterial": ("cm_celltype", ["OFT/RV-CM"]),
        "heart_venous": ("cm_celltype", ["aCM1", "aCM2", "JCF/SV-CM"]),
        "body_anterior": ("celltype", ["Forebrain", "aPHM"]),
        "body_posterior": ("celltype", ["Caudal Epiblast", "Hindgut", "SOM"]),
        "dorsal": ("celltype", ["Neural Tube"]),
        "ventral": ("core", [True]),
    },
    "E9.5": {
        "heart_arterial": ("celltype", ["OFT-CM"]),
        "heart_venous": ("celltype", ["SV-CM"]),
        "body_anterior": ("celltype", ["Branch Arch"]),
        "body_posterior": ("celltype", ["Hepatocytes", "ST"]),
        "dorsal": ("celltype", ["Dorsal-PAM"]),
        "ventral": ("core", [True]),
    },
    "E9.5_mab21l2_ko": {
        "heart_arterial": ("celltype", ["OFT/RV-CM"]),
        "heart_venous": ("celltype", ["SV-CM"]),
        "body_anterior": ("celltype", ["Branchial Arch"]),
        "body_posterior": ("celltype", ["Hepatocyte"]),
        "dorsal": ("celltype", ["Paraxial Mesoderm", "Neural Tube"]),
        "ventral": ("core", [True]),
    },
}
LR_GENES = ["Pitx2", "Nodal", "Lefty2", "Lefty1"]


def mask(a, spec):
    col, labs = spec
    v = np.asarray(a.obs[col])
    return np.isin(v, labs)


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v * np.nan


def handed_index(ap, dv, lr):
    e1 = unit(ap)
    e2 = unit(dv - (dv @ e1) * e1)
    e3 = unit(lr)
    return float(np.linalg.det(np.stack([e1, e2, e3])))


def weighted_offset(C, w):
    if w.sum() <= 0:
        return np.full(3, np.nan)
    return (w[:, None] * C).sum(0) / w.sum() - C.mean(0)


def canonical_basis(C):
    """The scorer's canonical frame: rows of Vt (PCA axes), last row flipped so det = +1."""
    X = C - C.mean(0)
    _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    if np.linalg.det(Vt) < 0:
        Vt[-1] *= -1
    return Vt


PROPER = [np.diag(f) for f in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]


def main(out: Path, n_boot: int = 200, n_perm: int = 200, seed: int = 0) -> None:
    from common.shape_metrics import occupancy_dice

    rng = np.random.default_rng(seed)
    vec_rows, h_rows, can_rows = [], [], []
    store = {}
    for stage, spec in LANDMARKS.items():
        a = load_stage(FILES[stage])
        C = coords_of(a)
        core = np.asarray(a.obs["core"], dtype=bool)
        groups = np.asarray(a.obs["group"]).astype(str)
        m = {k: mask(a, s) for k, s in spec.items()}
        lrcells = np.isin(groups, list(LR_GROUPS))
        vecs = {
            "AP_heart": C[m["heart_arterial"]].mean(0) - C[m["heart_venous"]].mean(0),
            "AP_body": C[m["body_anterior"]].mean(0) - C[m["body_posterior"]].mean(0),
            "DV": C[m["dorsal"]].mean(0) - C[m["ventral"]].mean(0),
        }
        counts = {k: int(v.sum()) for k, v in m.items()}
        expr = {}
        for g in LR_GENES:
            if g not in a.var_names:
                continue
            w = np.asarray(a[:, g].X.toarray()).ravel()
            expr[g] = w
            sub = lrcells & np.isfinite(w)
            off = weighted_offset(C[sub], w[sub])
            # permutation z-score of the offset length
            null = np.array([np.linalg.norm(weighted_offset(C[sub], rng.permutation(w[sub])))
                             for _ in range(n_perm)])
            vecs[f"LR_{g}"] = off
            vec_rows.append({"stage": stage, "vector": f"LR_{g}", "x": off[0], "y": off[1], "z": off[2],
                             "length": float(np.linalg.norm(off)), "n_cells": int(sub.sum()),
                             "frac_expressing": float((w[sub] > 0).mean()),
                             "perm_z": float((np.linalg.norm(off) - null.mean()) / max(null.std(), 1e-12))})
        for k in ("AP_heart", "AP_body", "DV"):
            v = vecs[k]
            vec_rows.append({"stage": stage, "vector": k, "x": v[0], "y": v[1], "z": v[2],
                             "length": float(np.linalg.norm(v)), "n_cells": counts, "frac_expressing": np.nan,
                             "perm_z": np.nan})
        # handedness index for every AP x LR definition, with bootstrap sign stability
        for ap_key, lr_key in itertools.product(["AP_heart", "AP_body"], [k for k in vecs if k.startswith("LR_")]):
            h = handed_index(vecs[ap_key], vecs["DV"], vecs[lr_key])
            g = lr_key[3:]
            boots = []
            for _ in range(n_boot):
                idx = rng.integers(0, C.shape[0], C.shape[0])
                Cb, mb = C[idx], {k: v[idx] for k, v in m.items()}
                lb, wb = lrcells[idx], expr[g][idx]
                if ap_key == "AP_heart":
                    ap = Cb[mb["heart_arterial"]].mean(0) - Cb[mb["heart_venous"]].mean(0)
                else:
                    ap = Cb[mb["body_anterior"]].mean(0) - Cb[mb["body_posterior"]].mean(0)
                dv = Cb[mb["dorsal"]].mean(0) - Cb[mb["ventral"]].mean(0)
                lr = weighted_offset(Cb[lb], wb[lb])
                boots.append(np.sign(handed_index(ap, dv, lr)))
            boots = np.array(boots)
            h_rows.append({"stage": stage, "AP": ap_key, "LR": lr_key, "h": h,
                           "boot_frac_same_sign": float((boots == np.sign(h)).mean()),
                           "cos_AP_DV": float(unit(vecs[ap_key]) @ unit(vecs["DV"])),
                           "cos_AP_LR": float(unit(vecs[ap_key]) @ unit(vecs[lr_key])),
                           "cos_DV_LR": float(unit(vecs["DV"]) @ unit(vecs[lr_key]))})
        # canonical frames (full and core)
        for cloud, sel in [("full", np.ones(C.shape[0], bool)), ("core", core)]:
            Vt = canonical_basis(C[sel])
            for k, v in vecs.items():
                u = Vt @ unit(v)
                can_rows.append({"stage": stage, "cloud": cloud, "vector": k, "pc1": u[0], "pc2": u[1], "pc3": u[2]})
            store[(stage, cloud)] = (C[sel], Vt, {k: Vt @ unit(v) for k, v in vecs.items()})
        hh = pd.DataFrame([r for r in h_rows if r["stage"] == stage])
        print(stage, "landmark counts", counts, flush=True)
        print(hh[["AP", "LR", "h", "boot_frac_same_sign", "cos_AP_DV"]].round(3).to_string(index=False), flush=True)

    # consecutive / paired frames
    pairs = [("E8.0", "E8.25_late"), ("E8.25_late", "E8.75"), ("E8.75", "E9.5"), ("E9.5", "E9.5_mab21l2_ko"),
             ("E8.75", "E9.5_mab21l2_ko")]
    pair_rows = []
    for (sa, sb), cloud in itertools.product(pairs, ["full", "core"]):
        CA, _, ua = store[(sa, cloud)]
        CB, _, ub = store[(sb, cloud)]
        keys = [k for k in ("AP_heart", "AP_body", "DV", "LR_Pitx2") if k in ua and k in ub]
        def score(F, ua=ua, ub=ub, keys=keys):
            return float(np.mean([(F @ ua[k]) @ ub[k] for k in keys]))
        prop = max(score(F) for F in PROPER)
        impr = max(score(-F) for F in PROPER)
        mir = CA * np.array([-1.0, 1.0, 1.0])
        d_p = occupancy_dice(CA, CB, seed=0)[0]
        d_m = occupancy_dice(mir, CB, seed=0)[0]
        pair_rows.append({"A": sa, "B": sb, "cloud": cloud, "axes_used": ",".join(keys),
                          "best_proper_mean_cos": prop, "best_improper_mean_cos": impr,
                          "dice_proper": d_p, "dice_mirrored": d_m})
    pd.DataFrame(vec_rows).to_csv(out / "handedness_vectors.csv", index=False)
    pd.DataFrame(h_rows).to_csv(out / "handedness_index.csv", index=False)
    pd.DataFrame(can_rows).to_csv(out / "canonical_axes.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(out / "frame_pairs.csv", index=False)
    print(pd.DataFrame(pair_rows).round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    main(Path(ap.parse_args().out))
