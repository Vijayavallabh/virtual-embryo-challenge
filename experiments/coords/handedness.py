"""Part 3: handedness of each file's frame (DGX).

A. Landmark index. Anatomical direction vectors in the raw frame (duplicated records dropped first):
     AP  cranial: (heart) arterial pole minus venous pole; (body) anterior minus posterior tissue
     DV  dorsal: dorsal tissue minus ventral tissue (core-cardiac cells; ventral foregut when the dorsal landmark is
         the dorsal foregut)
     LR  left: weighted centroid offset of a left-side marker (Pitx2 decides; Lefty2, Nodal, Lefty1 for context)
         within a cell set, relative to the same cells' centroid
   h = det[AP, DV orthogonalised to AP, LR] on unit vectors, i.e. the signed part of LR normal to the AP-DV plane.
   x = left, y = cranial, z = ventral is right-handed and det[cranial, dorsal, left] = -1, so a right-handed frame of a
   situs-solitus embryo gives h < 0; opposite signs in two files mean mirror images.
   Landmarks come from two sources, both applied to the released files: LANDMARKS, hand-picked from label names on
   2026-10-18 (released files only), and LANDMARK_RULE, a label-pattern rule fixed on 2026-10-19 before the release,
   which also picks landmarks for E8.5 and E10.5.
   Robustness (the 2026-10-19 review): LR cell sets (default = mesoderm and cardiac groups, non-cardiac mesoderm, CM
   only, all cells) x weights (expression, binary) x both AP definitions. Rule verdict per file and landmark source:
   "resolved" only when h has the same sign for Pitx2 in every variant, otherwise "unresolved"; the final verdict
   needs every available source resolved with one sign. A secondary verdict, added after the first run of
   2026-10-19 and labelled as such, ignores variants with |h| < H_INFORMATIVE (LR within about 6 degrees of the AP-DV
   plane, which carries no handedness information) and needs at least half the variants to remain. The bootstrap
   same-sign fraction (default variant) is kept, but it cannot see systematic errors from ill-conditioned axes; the
   variants can. E8.0 is unresolved (lead decision 2026-10-19). E6.75, E7.25 and E7.5 have no heart poles and no
   lateralised marker: not testable this way.
B. Mirror cost on matched shapes: veckit occupancy_dice of half vs half and of mirror(half) vs half, split halves of
   the 10% working copy of the file as released (as the scorer sees it), MIRROR_SEEDS seeds.
C. Shape alignment, proper vs improper: the best rigid fit of stage A onto stage B over all proper rotations and over
   all improper ones (mirror, then rotate). Random-rotation search, then symmetric ICP from the best starts; scored by
   the symmetric mean nearest-neighbour distance on RMS-normalised clouds (and occupancy Dice on the scorer's grid).
   Calibration: split halves of one stage (how chiral the shape is at this resolution) and consecutive released
   stages; at the release, E7.5 against E7.25 and E8.0 (run_p3.sh). Verdict per pair: "proper" or "improper" when
   every seed agrees, else "undecided". usable_verdict (rule added 2026-10-19 after the first run, before any
   release file was seen): a pair counts only if both stages are chiral, i.e. their own split halves give "proper"
   in every seed; otherwise a mirror image cannot be told apart and the pair is "uninformative".
D. Canonical frames (the scorer's PCA frame) and the diagonal-flip comparison of consecutive stages, kept for
   continuity; the review showed it is uninformative when the PCA axis order changes between stages.

Writes to --out: handedness_vectors.csv, handedness_index.csv, handedness_robustness.csv, handedness_verdict.csv,
landmarks_used.csv, mirror_dice.csv, mirror_dice_summary.csv, shape_alignment.csv, shape_alignment_summary.csv,
canonical_axes.csv, frame_pairs.csv.
Usage: uv run python experiments/coords/handedness.py [--out runs/coords] [--stage NAME=PATH ...] [--pairs A:B ...]
"""

from __future__ import annotations

import argparse
import itertools
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coordlib import FILES, coords_of, load, raw_coords, stage_rng, stage_specs  # noqa: E402

LR_SETS = {
    "default": {"CM", "CardiacProgenitor", "PharyngealMesoderm", "LPM", "Pericardium", "EpicardiumST", "Endocardium"},
    "noncardiac_meso": {"CardiacProgenitor", "PharyngealMesoderm", "LPM", "Pericardium", "EpicardiumST"},
    "CM_only": {"CM"},
    "all_cells": None,
}
WEIGHTS = ("expression", "binary")
LR_GENES = ["Pitx2", "Nodal", "Lefty2", "Lefty1"]
VERDICT_GENE = "Pitx2"
LEAD_VERDICTS = {"E8.0": "unresolved (lead decision 2026-10-19)"}
H_INFORMATIVE = 0.1  # secondary verdict only: |h| below this is treated as no information
LABEL_COLUMNS = ("celltype", "cm_celltype")
MIN_LANDMARK_CELLS = 50
MIRROR = np.array([-1.0, 1.0, 1.0])
MIRROR_SEEDS = 10
ALIGN_SEEDS = 5
ALIGN_N = 2000

# Hand-picked landmarks for the released files: (obs column, labels) per pole, chosen from label names before any
# geometry was computed (2026-10-18).
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

# Label-pattern rule (fixed 2026-10-19, before the release). Per landmark, tiers of regular expressions tried in
# order; within a tier the obs columns are tried in LABEL_COLUMNS order, and the first column whose matching labels
# hold at least MIN_LANDMARK_CELLS cells gives the landmark. Ventral is the core-cardiac cells, except that a dorsal
# landmark from the dorsal foregut is paired with the ventral foregut.
LANDMARK_RULE = {
    "heart_arterial": [[r"OFT"], [r"^aSHF"]],
    "heart_venous": [[r"SV"], [r"^IFT", r"^A-CM$", r"^aCM\d*$"], [r"^pSHF$"]],
    "body_anterior": [[r"^Forebrain$"], [r"^Branch(ial)? Arch$"], [r"^aPHM$"]],
    "body_posterior": [[r"^Caudal Epiblast$", r"^Primitive Streak$", r"^Hindgut$"], [r"^SOM$", r"^Allantois$"],
                       [r"^Hepatocytes?$", r"^ST$"]],
    "dorsal": [[r"^Neural Tube$"], [r"^Dorsal-PAM$", r"^Paraxial Mesoderm$"], [r"^D-FG$"]],
}
VENTRAL_FOR_DORSAL_FOREGUT = ("celltype", [r"^V-FG$"])

DEFAULT_PAIRS = [("E6.75", "E7.25"), ("E7.25", "E8.0"), ("E6.75", "E8.0"), ("E8.25_late", "E8.75"),
                 ("E8.75", "E9.5"), ("E9.5", "E9.5_mab21l2_ko")]


@dataclass
class StageData:
    name: str
    C: np.ndarray
    core: np.ndarray
    groups: np.ndarray
    labels: dict[str, np.ndarray]
    expr: dict[str, np.ndarray] = field(default_factory=dict)
    raw_C: np.ndarray | None = None
    n_dropped: int = 0
    set_masks: dict[str, np.ndarray] = field(default_factory=dict)


def read_stage(name: str, path: str, need_raw: bool) -> StageData:
    a, dd = load(path, drop_duplicates=True)
    labels = {c: np.asarray(a.obs[c]).astype(str) for c in LABEL_COLUMNS if c in a.obs}
    expr = {g: np.asarray(a[:, g].X.toarray()).ravel().astype(np.float64) for g in LR_GENES if g in a.var_names}
    groups = np.asarray(a.obs["group"]).astype(str)
    set_masks = {k: np.ones(len(groups), bool) if v is None else np.isin(groups, list(v)) for k, v in LR_SETS.items()}
    return StageData(name=name, C=coords_of(a).copy(), core=np.asarray(a.obs["core"], dtype=bool), groups=groups,
                     labels=labels, expr=expr, raw_C=raw_coords(path) if need_raw else None,
                     n_dropped=dd["n_dropped"], set_masks=set_masks)


# ----------------------------------------------------------------------------------------------- A. landmarks
def rule_landmarks(sd: StageData) -> dict[str, tuple[str, list[str]]]:
    """Landmarks picked by LANDMARK_RULE from the labels present in the file (missing landmarks are left out)."""
    out = {}
    for lm, tiers in LANDMARK_RULE.items():
        for tier in tiers:
            hit = None
            for col in LABEL_COLUMNS:
                if col not in sd.labels:
                    continue
                vc = pd.Series(sd.labels[col]).value_counts()
                labs = [lab for lab in vc.index if any(re.search(p, lab) for p in tier)]
                if labs and int(vc[labs].sum()) >= MIN_LANDMARK_CELLS:
                    hit = (col, sorted(labs))
                    break
            if hit:
                out[lm] = hit
                break
    if "dorsal" in out:
        if any(re.search(r"^D-FG$", lab) for lab in out["dorsal"][1]):
            col, pats = VENTRAL_FOR_DORSAL_FOREGUT
            present = np.unique(sd.labels[col]) if col in sd.labels else []
            labs = sorted(lab for lab in present if any(re.search(p, lab) for p in pats))
            if labs:
                out["ventral"] = (col, labs)
        else:
            out["ventral"] = ("core", [True])
    return out


def landmark_masks(sd: StageData, spec: dict[str, tuple[str, list]]) -> dict[str, np.ndarray]:
    masks = {}
    for lm, (col, labs) in spec.items():
        v = sd.core if col == "core" else sd.labels.get(col)
        if v is not None:
            masks[lm] = np.isin(v, labs)
    return masks


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v * np.nan


def handed_index(ap, dv, lr) -> float:
    e1 = unit(ap)
    e2 = unit(dv - (dv @ e1) * e1)
    return float(np.linalg.det(np.stack([e1, e2, unit(lr)])))


def axis_vectors(C: np.ndarray, masks: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    def diff(a, b):
        if a in masks and b in masks and masks[a].any() and masks[b].any():
            return C[masks[a]].mean(0) - C[masks[b]].mean(0)
        return None

    out = {"AP_heart": diff("heart_arterial", "heart_venous"), "AP_body": diff("body_anterior", "body_posterior"),
           "DV": diff("dorsal", "ventral")}
    return {k: v for k, v in out.items() if v is not None}


def weighted_offset(C, w):
    if w.sum() <= 0:
        return np.full(3, np.nan)
    return (w[:, None] * C).sum(0) / w.sum() - C.mean(0)


def lr_vector(sd: StageData, gene: str, lr_set: str, weight: str, idx: np.ndarray | None = None) -> np.ndarray:
    C, sel, w = sd.C, sd.set_masks[lr_set], sd.expr[gene]
    if idx is not None:
        C, sel, w = C[idx], sel[idx], w[idx]
    ww = w[sel] if weight == "expression" else (w[sel] > 0).astype(float)
    return weighted_offset(C[sel], ww)


def landmark_handedness(sd: StageData, source: str, spec: dict, rng: np.random.Generator, n_boot: int, n_perm: int):
    masks = landmark_masks(sd, spec)
    vecs = axis_vectors(sd.C, masks)
    used = [{"stage": sd.name, "source": source, "landmark": lm, "column": col, "labels": "|".join(map(str, labs)),
             "n_cells": int(masks[lm].sum()) if lm in masks else 0} for lm, (col, labs) in spec.items()]
    vec_rows, rob_rows, idx_rows = [], [], []
    for k, v in vecs.items():
        vec_rows.append({"stage": sd.name, "source": source, "vector": k, "x": v[0], "y": v[1], "z": v[2],
                         "length": float(np.linalg.norm(v)), "perm_z": np.nan})
    aps = [k for k in ("AP_heart", "AP_body") if k in vecs]
    if "DV" not in vecs or not aps:
        return used, vec_rows, rob_rows, idx_rows
    for gene in sd.expr:
        for lr_set, weight in itertools.product(LR_SETS, WEIGHTS):
            lr = lr_vector(sd, gene, lr_set, weight)
            if not np.isfinite(lr).all():
                continue
            for ap in aps:
                rob_rows.append({"stage": sd.name, "source": source, "AP": ap, "gene": gene, "lr_set": lr_set,
                                 "weight": weight, "h": handed_index(vecs[ap], vecs["DV"], lr),
                                 "cos_AP_DV": float(unit(vecs[ap]) @ unit(vecs["DV"])),
                                 "cos_DV_LR": float(unit(vecs["DV"]) @ unit(lr))})
        # default variant: permutation z of the offset length and bootstrap sign stability
        sel = sd.set_masks["default"]
        w = sd.expr[gene][sel]
        off = weighted_offset(sd.C[sel], w)
        null = np.array([np.linalg.norm(weighted_offset(sd.C[sel], rng.permutation(w))) for _ in range(n_perm)])
        vec_rows.append({"stage": sd.name, "source": source, "vector": f"LR_{gene}", "x": off[0], "y": off[1],
                         "z": off[2], "length": float(np.linalg.norm(off)),
                         "perm_z": float((np.linalg.norm(off) - null.mean()) / max(null.std(), 1e-12))})
    boots = {(ap, g): [] for ap in aps for g in sd.expr}
    n = len(sd.C)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        bm = {k: v[idx] for k, v in masks.items()}
        bv = axis_vectors(sd.C[idx], bm)
        if "DV" not in bv:
            continue
        for g in sd.expr:
            lr = lr_vector(sd, g, "default", "expression", idx)
            for ap in aps:
                if ap in bv:
                    boots[(ap, g)].append(np.sign(handed_index(bv[ap], bv["DV"], lr)))
    for (ap, g), b in boots.items():
        h = [r["h"] for r in rob_rows if r["AP"] == ap and r["gene"] == g and r["lr_set"] == "default"
             and r["weight"] == "expression"]
        if not h:
            continue
        idx_rows.append({"stage": sd.name, "source": source, "AP": ap, "LR": f"LR_{g}", "h": h[0],
                         "boot_frac_same_sign": float(np.mean(np.array(b) == np.sign(h[0]))) if b else np.nan,
                         "cos_AP_DV": float(unit(vecs[ap]) @ unit(vecs["DV"]))})
    return used, vec_rows, rob_rows, idx_rows


def rule_verdict(rob: pd.DataFrame, h_min: float = 0.0) -> tuple[str, int, int]:
    """Verdict from the Pitx2 variants. h_min > 0 (secondary verdict) drops variants with |h| < h_min first and needs
    at least half of the variants to remain."""
    r = rob[rob.gene == VERDICT_GENE] if len(rob) else rob
    if len(r) == 0:
        return "not testable (landmarks or the marker missing)", 0, 0
    n_all = len(r)
    r = r[r.h.abs() >= h_min]
    neg, pos = int((r.h < 0).sum()), int((r.h > 0).sum())
    if 2 * len(r) < n_all:
        return f"not enough informative variants ({len(r)} of {n_all})", neg, pos
    if pos == 0:
        return "resolved: h < 0 in every variant (right-handed frame, situs solitus)", neg, pos
    if neg == 0:
        return "resolved: h > 0 in every variant (mirror image of a right-handed situs-solitus frame)", neg, pos
    return "unresolved: the sign of h changes across variants", neg, pos


def combine_verdicts(texts: list[str]) -> str:
    """Final verdict over landmark sources: every testable source resolved with one sign."""
    texts = [t for t in texts if t.startswith(("resolved", "unresolved"))]
    if not texts:
        return "not testable by landmarks"
    if all(t.startswith("resolved: h < 0") for t in texts):
        return "resolved: h < 0 (right-handed frame, situs solitus)"
    if all(t.startswith("resolved: h > 0") for t in texts):
        return "resolved: h > 0 (mirror image)"
    return "unresolved"


# ----------------------------------------------------------------------------------------------- B. mirror Dice
def mirror_dice(sd: StageData, seeds: int) -> list[dict]:
    from common.shape_metrics import occupancy_dice

    from vec.evaluate import Stage, split_half, subsample

    S = Stage(X=np.zeros((len(sd.raw_C), 1), np.float32), coords=sd.raw_C)
    rows = []
    for s in range(seeds):
        h1, h2 = split_half(subsample(S, 0.1, s), seed=0)
        rows.append({"stage": sd.name, "seed": s, "n_half": int(h1.n),
                     "dice_proper": float(occupancy_dice(h1.coords, h2.coords, seed=0)[0]),
                     "dice_mirrored": float(occupancy_dice(h1.coords * MIRROR, h2.coords, seed=0)[0])})
    return rows


# ----------------------------------------------------------------------------------------------- C. shape alignment
def random_rotations(n: int, rng: np.random.Generator) -> np.ndarray:
    """n uniformly distributed proper rotations (unit quaternions), plus the identity first."""
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q.T
    R = np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
                  np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
                  np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], 1)
    return np.concatenate([np.eye(3)[None], R])


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Proper rotation R and translation t minimising sum |R p + t - q|^2."""
    mp, mq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - mp).T @ (Q - mq))
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, mq - R @ mp


def chamfer(A: np.ndarray, B: np.ndarray, tree_b: cKDTree | None = None) -> float:
    """Symmetric mean nearest-neighbour distance."""
    tree_b = tree_b or cKDTree(B)
    return 0.5 * float(tree_b.query(A)[0].mean() + cKDTree(A).query(B)[0].mean())


def normalise(C: np.ndarray) -> np.ndarray:
    X = C - C.mean(0)
    return X / np.sqrt((X ** 2).sum(1).mean())


def dice_aligned(A: np.ndarray, B: np.ndarray, n_grid: int = 16, extent: float = 3.0) -> float:
    """Occupancy Dice of two already aligned, RMS-normalised clouds on the scorer's grid (16^3 over [-3, 3])."""
    def occ(Z):
        v = np.floor((np.clip(Z, -extent, extent - 1e-9) + extent) / (2 * extent) * n_grid).astype(int)
        return set(map(tuple, v))

    oa, ob = occ(A), occ(B)
    return 2 * len(oa & ob) / max(len(oa) + len(ob), 1)


def best_fit(A: np.ndarray, B: np.ndarray, improper: bool, rng: np.random.Generator, n_rot: int = 600,
             n_coarse: int = 400, n_refine: int = 8, iters: int = 25) -> tuple[float, float]:
    """Best rigid fit of A onto B (both RMS-normalised) over proper rotations, or over improper ones (A mirrored,
    then rotated). Returns (symmetric mean NN distance, occupancy Dice) at the best fit."""
    A = A * MIRROR if improper else A
    tree_b = cKDTree(B)
    Rs = random_rotations(n_rot, rng)
    Ac = A[rng.choice(len(A), min(n_coarse, len(A)), replace=False)]
    Bc = B[rng.choice(len(B), min(n_coarse, len(B)), replace=False)]
    tree_bc = cKDTree(Bc)
    coarse = np.array([chamfer(Ac @ R.T, Bc, tree_bc) for R in Rs])
    best = (np.inf, None)
    for i in np.argsort(coarse)[:n_refine]:
        R, t = Rs[i], np.zeros(3)
        for _ in range(iters):
            AT = A @ R.T + t
            _, jb = tree_b.query(AT)
            _, ja = cKDTree(AT).query(B)
            R, t = kabsch(np.vstack([A, A[ja]]), np.vstack([B[jb], B]))
        AT = A @ R.T + t
        d = chamfer(AT, B, tree_b)
        if d < best[0]:
            best = (d, AT)
    return float(best[0]), dice_aligned(best[1], B)


def align_task(name_a: str, name_b: str, kind: str, CA: np.ndarray, CB: np.ndarray | None, seed: int) -> dict:
    """One seed of the proper-vs-improper test. CB None: split halves of CA (same-stage calibration)."""
    rng = stage_rng(f"{name_a}|{name_b}|{kind}", seed)
    if CB is None:
        perm = rng.permutation(len(CA))
        CA, CB = CA[perm[: len(CA) // 2]], CA[perm[len(CA) // 2:]]
    A = normalise(CA[rng.choice(len(CA), min(ALIGN_N, len(CA)), replace=False)])
    B = normalise(CB[rng.choice(len(CB), min(ALIGN_N, len(CB)), replace=False)])
    dp, dice_p = best_fit(A, B, False, rng)
    di, dice_i = best_fit(A, B, True, rng)
    return {"A": name_a, "B": name_b, "kind": kind, "seed": seed, "n": int(min(len(A), len(B))),
            "nn_dist_proper": dp, "nn_dist_improper": di, "rel_margin": di / dp - 1.0,
            "dice_proper": dice_p, "dice_improper": dice_i}


def summarise_alignment(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (a, b, kind), g in df.groupby(["A", "B", "kind"], sort=False):
        prop = (g.nn_dist_proper < g.nn_dist_improper).all()
        impr = (g.nn_dist_improper < g.nn_dist_proper).all()
        rows.append({"A": a, "B": b, "kind": kind, "seeds": len(g),
                     "nn_dist_proper_mean": g.nn_dist_proper.mean(), "nn_dist_improper_mean": g.nn_dist_improper.mean(),
                     "rel_margin_mean": g.rel_margin.mean(), "rel_margin_sd": g.rel_margin.std(ddof=1),
                     "dice_proper_mean": g.dice_proper.mean(), "dice_improper_mean": g.dice_improper.mean(),
                     "verdict": "proper" if prop else "improper" if impr else "undecided"})
    out = pd.DataFrame(rows)
    chiral = {r.A: r.verdict == "proper" for r in out.itertuples() if r.kind == "same_stage_halves"}

    def usable(r) -> str:
        if r.kind == "same_stage_halves":
            return "chiral" if r.verdict == "proper" else "not chiral: split halves do not tell a mirror image apart"
        missing = [s for s in (r.A, r.B) if s not in chiral]
        if missing:
            return f"uninformative: no split-half calibration for {', '.join(missing)}"
        if not (chiral[r.A] and chiral[r.B]):
            return "uninformative: " + ", ".join(s for s in (r.A, r.B) if not chiral[s]) + " not chiral"
        return r.verdict

    out["usable_verdict"] = [usable(r) for r in out.itertuples()]
    return out


# ----------------------------------------------------------------------------------------------- D. canonical frames
def canonical_basis(C):
    """The scorer's canonical frame: rows of Vt (PCA axes), last row flipped so det = +1."""
    X = C - C.mean(0)
    _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    if np.linalg.det(Vt) < 0:
        Vt[-1] *= -1
    return Vt


PROPER = [np.diag(f) for f in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]


def main(out: Path, files: dict[str, str], landmark_stages: list[str], pairs: list[tuple[str, str]],
         mirror_stages: list[str], same_stages: list[str], n_boot: int = 200, n_perm: int = 200, seed: int = 0,
         jobs: int = 8) -> None:
    from joblib import Parallel, delayed

    out.mkdir(parents=True, exist_ok=True)
    needed = sorted(set(landmark_stages) | {s for p in pairs for s in p} | set(mirror_stages) | set(same_stages),
                    key=list(files).index)
    data = {s: read_stage(s, files[s], need_raw=s in mirror_stages) for s in needed}
    for s in needed:
        print(f"loaded {s}: {len(data[s].C)} cells after dropping {data[s].n_dropped} duplicates", flush=True)

    # A. landmarks
    used, vec_rows, rob_rows, idx_rows, verdicts, store = [], [], [], [], [], {}
    for s in landmark_stages:
        sd = data[s]
        sources = {"rule": rule_landmarks(sd)} | ({"hand": LANDMARKS[s]} if s in LANDMARKS else {})
        per_source = {}
        for src, spec in sources.items():
            u, v, r, i = landmark_handedness(sd, src, spec, stage_rng(f"{s}|{src}", seed), n_boot, n_perm)
            used += u
            vec_rows += v
            rob_rows += r
            idx_rows += i
            rob = pd.DataFrame(r)
            text, neg, pos = rule_verdict(rob)
            text_i, neg_i, pos_i = rule_verdict(rob, H_INFORMATIVE)
            per_source[src] = (text, text_i)
            cos = rob[rob.gene == VERDICT_GENE].cos_AP_DV.abs().max() if len(rob) else np.nan
            verdicts.append({"stage": s, "source": src, "rule_verdict": text, "pitx2_variants_neg": neg,
                             "pitx2_variants_pos": pos, "max_abs_cos_AP_DV": cos,
                             "secondary_verdict_informative_only": text_i, "pitx2_informative_neg": neg_i,
                             "pitx2_informative_pos": pos_i})
            if src == ("hand" if s in LANDMARKS else "rule"):  # axes for the canonical-frame comparison (D)
                vecs = axis_vectors(sd.C, landmark_masks(sd, spec))
                if VERDICT_GENE in sd.expr:
                    vecs[f"LR_{VERDICT_GENE}"] = lr_vector(sd, VERDICT_GENE, "default", "expression")
                store[s] = vecs
        final = combine_verdicts([t for t, _ in per_source.values()])
        final_i = combine_verdicts([t for _, t in per_source.values()])
        if s in LEAD_VERDICTS:
            final = LEAD_VERDICTS[s] + ("" if final == "unresolved" else f"; rule gave: {final}")
            final_i = LEAD_VERDICTS[s] + ("" if final_i == "unresolved" else f"; secondary rule gave: {final_i}")
        for v in verdicts:
            if v["stage"] == s:
                v["final_verdict"] = final
                v["final_secondary_verdict"] = final_i
        print(s, "verdicts:", per_source, "->", final, "| secondary:", final_i, flush=True)

    # B. mirror cost on matched shapes
    mir = pd.DataFrame([r for s in mirror_stages for r in mirror_dice(data[s], MIRROR_SEEDS)])
    if len(mir):
        mir["diff"] = mir.dice_proper - mir.dice_mirrored
        mir_sum = mir.groupby("stage", sort=False).agg(dice_proper_mean=("dice_proper", "mean"),
                                                       dice_proper_sd=("dice_proper", "std"),
                                                       dice_mirrored_mean=("dice_mirrored", "mean"),
                                                       dice_mirrored_sd=("dice_mirrored", "std"),
                                                       diff_mean=("diff", "mean"), diff_min=("diff", "min"),
                                                       n_half=("n_half", "first")).reset_index()
    else:
        mir_sum = pd.DataFrame()

    # C. shape alignment, proper vs improper
    tasks = [(s, s, "same_stage_halves", data[s].C, None, k) for s in same_stages for k in range(ALIGN_SEEDS)]
    tasks += [(a, b, "pair", data[a].C, data[b].C, k) for a, b in pairs for k in range(ALIGN_SEEDS)]
    al = pd.DataFrame(Parallel(n_jobs=jobs)(delayed(align_task)(*t) for t in tasks)) if tasks else pd.DataFrame()
    al_sum = summarise_alignment(al) if len(al) else pd.DataFrame()

    # D. canonical frames and diagonal flips of consecutive landmark stages
    can_rows, pair_rows = [], []
    for s, vecs in store.items():
        Vt = canonical_basis(data[s].C)
        for k, v in vecs.items():
            u = Vt @ unit(v)
            can_rows.append({"stage": s, "cloud": "full", "vector": k, "pc1": u[0], "pc2": u[1], "pc3": u[2]})
    order = [s for s in landmark_stages if s in store]
    for sa, sb in zip(order[:-1], order[1:], strict=True):
        ua = {k: canonical_basis(data[sa].C) @ unit(v) for k, v in store[sa].items()}
        ub = {k: canonical_basis(data[sb].C) @ unit(v) for k, v in store[sb].items()}
        keys = [k for k in ("AP_heart", "AP_body", "DV", f"LR_{VERDICT_GENE}") if k in ua and k in ub]
        if "DV" not in keys or not {"AP_heart", "AP_body"} & set(keys):
            continue  # no anatomical frame to compare (stages without landmarks)

        def score(F, ua=ua, ub=ub, keys=keys):
            return float(np.mean([(F @ ua[k]) @ ub[k] for k in keys]))

        pair_rows.append({"A": sa, "B": sb, "axes_used": ",".join(keys),
                          "best_proper_mean_cos": max(score(F) for F in PROPER),
                          "best_improper_mean_cos": max(score(-F) for F in PROPER)})

    for name, rows in (("landmarks_used", used), ("handedness_vectors", vec_rows), ("handedness_index", idx_rows),
                       ("handedness_robustness", rob_rows), ("handedness_verdict", verdicts),
                       ("canonical_axes", can_rows), ("frame_pairs", pair_rows)):
        pd.DataFrame(rows).to_csv(out / f"{name}.csv", index=False)
    mir.to_csv(out / "mirror_dice.csv", index=False)
    mir_sum.to_csv(out / "mirror_dice_summary.csv", index=False)
    al.to_csv(out / "shape_alignment.csv", index=False)
    al_sum.to_csv(out / "shape_alignment_summary.csv", index=False)
    pd.set_option("display.width", 220)
    if verdicts:
        print(pd.DataFrame(verdicts).to_string(index=False))
    if rob_rows:
        rob = pd.DataFrame(rob_rows)
        print(rob[rob.gene == VERDICT_GENE].pivot_table(index=["stage", "source", "AP"], columns=["lr_set", "weight"],
                                                       values="h").round(2).to_string())
    if len(mir_sum):
        print(mir_sum.round(3).to_string(index=False))
    if len(al_sum):
        print(al_sum.round(3).to_string(index=False))


def _pairs(items: list[str] | None, default) -> list[tuple[str, str]]:
    if items is None:
        return list(default)
    out = []
    for it in items:
        a, sep, b = it.partition(":")
        if not sep:
            raise SystemExit(f"--pairs expects A:B, got {it!r}")
        out.append((a, b))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    ap.add_argument("--stage", action="append", help="NAME=PATH (repeatable; default: every released file)")
    ap.add_argument("--landmark-stages", nargs="*", default=None, help="default: every stage given")
    ap.add_argument("--pairs", nargs="*", default=None, help="A:B pairs for the shape-alignment test")
    ap.add_argument("--mirror-stages", nargs="*", default=None, help="default: every stage given")
    ap.add_argument("--same-stages", nargs="*", default=None, help="split-half calibration; default: every stage")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    files = stage_specs(args.stage, FILES)
    stages = list(files)
    main(Path(args.out), files,
         landmark_stages=stages if args.landmark_stages is None else args.landmark_stages,
         pairs=_pairs(args.pairs, [p for p in DEFAULT_PAIRS if p[0] in files and p[1] in files]),
         mirror_stages=stages if args.mirror_stages is None else args.mirror_stages,
         same_stages=stages if args.same_stages is None else args.same_stages, jobs=args.jobs)
