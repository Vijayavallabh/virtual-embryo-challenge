"""OIC125's pure computations (row R70; HEART_MOON1_design.md P1 with every fix text of HEART_MOON1_challenge.md
section 5 and its section-6 points where cheap): marker ranks, unit typing, the heart region, contrasts,
residual chords, the consensus and G0's checks. Arrays in, arrays out: nothing here reads a file, a target or a score.

Fixed here, before any outside matrix is read (the row states the same numbers):
- CATS, the four chamber categories (configs/t2_chamber_classes.yaml): cm, endocardium, epicardium and
  vascular_endothelium.
- Marker rule (design step 1): per category, every panel gene ranked by its one-vs-rest AUROC on the reference stage's
  organ-cut cells typed by their labels (log1p CP10k over the panel; ties by gene name). The ranks are computed and
  hashed from the reference alone; a source's markers are the first N_MARKERS ranked genes that the source holds.
- Unit typing (design step 2, FIX 4): a unit's CP10k over G_s (the panel genes the source holds), each marker scaled by
  its mean over the section's units (the per-gene platform efficiency cancels, up to each unit's
  renormalisation); a_k = mean over k's markers, score
  s_k = a_k / sum_j a_j. Units with fewer than MIN_UNIT_COUNTS counts over G_s are not typed.
- Heart region (FIX 4): the largest 8-connected component of units with argmax cm at margin >= MARGIN_REGION, dilated
  by DILATE_UM (25 um); a reader confirms it on a plot of total counts; inside it a unit takes its argmax at margin
  >= MARGIN_TYPE, else it is dropped.
- Second typing (G0(a)): a KNN-vote of the reference's labelled cells (15 neighbours) in a PCA of log1p(e / gene mean)
  over G_s; agreement overall >= AGREE_ALL, Jaccard per category >= AGREE_JACCARD.
- Contrasts (step 3, FIX 5): m_kg the category mean of unit (cell) CP10k over G_s; C_kg = log(m_kg + PSEUDO) minus its
  mean over the centre set E (categories with >= MIN_REF_CELLS reference cells and G0(a)'s counts in the source, the
  same set at both ends), then each category's vector less its median over G_s.
- Chord (FIX 1): per category k and source stage s, C_skg = alpha + b C_ref,kg + r_skg by least squares over the guarded
  genes (m_kg >= GUARD_CP10K in the source and in the reference: section 6 asks 5 CP10k over the design's 1); P_skg =
  a_s r_skg; b is printed, never applied; a study's stages are averaged where held.
- Consensus (FIX 2): P_kg = the median over studies, kept where at least two studies hold it and all share its sign,
  never on the typing genes; clipped at |P| <= log 2 (the clip count is reported).
- G0(c) (FIX 3): per lineage, Spearman of two independent studies' residuals >= SPEARMAN_MIN and above the 95th
  percentile of N_PERM permutations within PERM_BINS level bins (seed PERM_SEED); a lineage that fails gets no P.
- G0(d): at least N_P_MIN entries with |P| >= P_MIN (under FIX 1's "add", |P| is the change).
"""

from __future__ import annotations

import hashlib
import json
import warnings
import zlib

import numpy as np
from scipy import ndimage, sparse, stats

CATS = ("cm", "endocardium", "epicardium", "vascular_endothelium")
N_MARKERS = 20
MIN_UNIT_COUNTS = 20
MARGIN_REGION = 0.2
MARGIN_TYPE = 0.1
DILATE_UM = 25.0
KNN = 15
KNN_PCS = 30
AGREE_ALL = 0.85
AGREE_JACCARD = 0.70
MIN_TYPED = {"cm": 200, "endocardium": 200, "vascular_endothelium": 200, "epicardium": 100}
MIN_SHARED_GENES = 150
MIN_REF_CELLS = 30
MIN_START_CELLS = 20  # momentum_shift's min_cells: a category it would skip at the start stage carries no P
PSEUDO = 0.5
GUARD_CP10K = 5.0
CLIP = float(np.log(2.0))
P_MIN = 0.10
N_P_MIN = 50
SPEARMAN_MIN = 0.25
N_PERM = 1000
PERM_SEED = 20261126
PERM_BINS = 10
TARGET_DAY = 12.5  # E12.5, the test's target
P1_DAY = 20.5  # P1: birth at E19.5 is P0

RULES = {"cats": CATS, "n_markers": N_MARKERS, "min_unit_counts": MIN_UNIT_COUNTS, "margin_region": MARGIN_REGION,
         "margin_type": MARGIN_TYPE, "dilate_um": DILATE_UM, "knn": KNN, "knn_pcs": KNN_PCS, "agree_all": AGREE_ALL,
         "agree_jaccard": AGREE_JACCARD, "min_typed": MIN_TYPED, "min_shared_genes": MIN_SHARED_GENES,
         "min_ref_cells": MIN_REF_CELLS, "min_start_cells": MIN_START_CELLS, "pseudo": PSEUDO,
         "guard_cp10k": GUARD_CP10K, "clip": CLIP, "p_min": P_MIN, "n_p_min": N_P_MIN, "spearman_min": SPEARMAN_MIN,
         "n_perm": N_PERM, "perm_seed": PERM_SEED, "perm_bins": PERM_BINS, "target_day": TARGET_DAY, "p1_day": P1_DAY}


def sha256_json(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def cp10k(C) -> np.ndarray:
    """Dense CP10k of a counts matrix (rows renormalised over its columns); a row with no count stays 0."""
    C = sparse.csr_matrix(C, dtype=np.float64)
    lib = np.asarray(C.sum(1)).ravel()
    scale = np.divide(1e4, lib, out=np.zeros_like(lib), where=lib > 0)
    return np.asarray((sparse.diags(scale) @ C).todense())


# ---- the marker rule ----------------------------------------------------------------------------------------------
def auroc_columns(X: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """One-vs-rest AUROC of every column of X for the rows in `pos` (Mann-Whitney U over n_pos n_neg; ties averaged)."""
    pos = np.asarray(pos, bool)
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return np.full(X.shape[1], np.nan)
    R = stats.rankdata(X, axis=0)
    return (R[pos].sum(0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def marker_ranks(C_ref, cats_ref, genes) -> dict:
    """Per category, every gene in descending one-vs-rest AUROC (ties by gene name) on the reference's cells, from
    log1p CP10k over all `genes`. The reference alone: no outside value enters."""
    genes = [str(g) for g in genes]
    X = np.log1p(cp10k(C_ref))
    cats_ref = np.asarray(cats_ref).astype(str)
    out = {}
    for k in CATS:
        au = auroc_columns(X, cats_ref == k)
        order = sorted(range(len(genes)), key=lambda i: (-(au[i] if np.isfinite(au[i]) else -1.0), genes[i]))
        out[k] = [genes[i] for i in order]
    return {"rule": "one-vs-rest AUROC on log1p CP10k over the panel, reference organ-cut cells typed by label",
            "genes": genes, "ranks": out}


def markers_for(ranks: dict, held, n: int = N_MARKERS) -> dict[str, list[str]]:
    """The first n ranked genes per category among the genes the source holds."""
    held = set(map(str, held))
    return {k: [g for g in ranks["ranks"][k] if g in held][:n] for k in CATS}


def typing_genes(markers: dict[str, list[str]]) -> set[str]:
    return {g for k in CATS for g in markers[k]}


# ---- unit typing and the heart region -----------------------------------------------------------------------------
def unit_scores(C_units, genes_gs, markers: dict, min_counts: int = MIN_UNIT_COUNTS) -> dict:
    """Marker scores of every unit: s (units x 4, rows summing to 1; NaN for an untyped unit), argmax and margin."""
    C = sparse.csr_matrix(C_units, dtype=np.float64)
    tot = np.asarray(C.sum(1)).ravel()
    E = cp10k(C)
    mu = E.mean(0)
    idx = {g: i for i, g in enumerate(map(str, genes_gs))}
    A = np.zeros((E.shape[0], len(CATS)))
    for j, k in enumerate(CATS):
        cols = [idx[g] for g in markers[k] if g in idx and mu[idx[g]] > 0]
        if cols:
            A[:, j] = (E[:, cols] / mu[cols]).mean(1)
    den = A.sum(1)
    ok = (tot >= min_counts) & (den > 0)
    S = np.full_like(A, np.nan)
    S[ok] = A[ok] / den[ok, None]
    srt = np.sort(np.where(ok[:, None], S, 0.0), axis=1)
    margin = np.where(ok, srt[:, -1] - srt[:, -2], np.nan)
    arg = np.where(ok, np.argmax(np.where(ok[:, None], S, -1.0), axis=1), -1)
    return {"S": S, "argmax": arg, "margin": margin, "typed": ok}


def disk(radius_units: float) -> np.ndarray:
    r = int(np.floor(radius_units))
    yy, xx = np.mgrid[-r: r + 1, -r: r + 1]
    return (xx ** 2 + yy ** 2) <= radius_units ** 2 + 1e-9


def heart_region(bx: np.ndarray, by: np.ndarray, cm_seed: np.ndarray, unit_um: float,
                 dilate_um: float = DILATE_UM) -> np.ndarray:
    """The largest 8-connected component of the seed units on the unit grid, dilated by dilate_um, restricted to the
    units present. Returns a boolean mask over the units (all False when no seed unit exists)."""
    bx, by = np.asarray(bx, np.int64), np.asarray(by, np.int64)
    seed = np.asarray(cm_seed, bool)
    if not seed.any():
        return np.zeros(len(bx), bool)
    x0, y0 = bx.min(), by.min()
    ix, iy = bx - x0, by - y0
    img = np.zeros((int(iy.max()) + 1, int(ix.max()) + 1), bool)
    img[iy[seed], ix[seed]] = True
    lab, n = ndimage.label(img, structure=np.ones((3, 3), bool))
    sizes = ndimage.sum(img, lab, index=np.arange(1, n + 1))
    best = int(np.argmax(sizes)) + 1  # ties: the lowest label (scan order)
    comp = lab == best
    grown = ndimage.binary_dilation(comp, structure=disk(dilate_um / unit_um))
    return grown[iy, ix]


def region_typing(sc: dict, region: np.ndarray, margin: float = MARGIN_TYPE) -> np.ndarray:
    """Inside the region, a unit's category index (argmax at margin >= `margin`), else -1."""
    ok = region & sc["typed"] & (np.nan_to_num(sc["margin"], nan=-1.0) >= margin)
    return np.where(ok, sc["argmax"], -1)


def knn_vote(C_units, C_ref, cats_ref, genes_gs, k: int = KNN, n_pcs: int = KNN_PCS) -> np.ndarray:
    """The second typing: each unit's majority category among its k nearest reference cells, in a PCA (fitted on the
    reference) of log1p(CP10k / the dataset's gene mean) over G_s. Index into CATS, or -1 for another category."""
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    def feats(C):
        E = cp10k(C)
        mu = E.mean(0)
        return np.log1p(np.divide(E, mu, out=np.zeros_like(E), where=mu > 0))

    Fr, Fu = feats(C_ref), feats(C_units)
    n = int(min(n_pcs, Fr.shape[1], max(1, Fr.shape[0] - 1)))
    pca = PCA(n_components=n, random_state=0, svd_solver="full").fit(Fr)
    Zr, Zu = pca.transform(Fr), pca.transform(Fu)
    nn = NearestNeighbors(n_neighbors=min(k, len(Zr))).fit(Zr)
    _, ind = nn.kneighbors(Zu)
    cats_ref = np.asarray(cats_ref).astype(str)
    code = np.array([CATS.index(c) if c in CATS else len(CATS) for c in cats_ref])
    votes = np.zeros((len(Zu), len(CATS) + 1), int)
    for j in range(ind.shape[1]):
        np.add.at(votes, (np.arange(len(Zu)), code[ind[:, j]]), 1)
    win = np.argmax(votes, axis=1)  # ties: the first in CATS order, "other" last
    return np.where(win < len(CATS), win, -1)


def agreement(t1: np.ndarray, t2: np.ndarray) -> dict:
    """Over the units typed by the marker rule (t1 >= 0): the share whose vote agrees, and per category the Jaccard of
    the two typings' unit sets."""
    t1, t2 = np.asarray(t1), np.asarray(t2)
    typed = t1 >= 0
    share = float((t1[typed] == t2[typed]).mean()) if typed.any() else float("nan")
    jac = {}
    for j, k in enumerate(CATS):
        a, b = typed & (t1 == j), typed & (t2 == j)
        u = int((a | b).sum())
        jac[k] = float((a & b).sum() / u) if u else float("nan")
    return {"share": share, "jaccard": jac, "typed": {k: int((t1 == j).sum()) for j, k in enumerate(CATS)}}


# ---- contrasts, residual chords and the consensus -----------------------------------------------------------------
def category_means(C, cat_idx: np.ndarray) -> dict[str, np.ndarray]:
    """Per category, the mean over its units (cells) of CP10k over the matrix's columns."""
    E = cp10k(C)
    cat_idx = np.asarray(cat_idx)
    return {k: E[cat_idx == j].mean(0) for j, k in enumerate(CATS) if (cat_idx == j).any()}


def contrasts(means: dict[str, np.ndarray], centre: list[str], pseudo: float = PSEUDO) -> dict[str, np.ndarray]:
    """C_kg = log(m_kg + pseudo) minus its mean over the centre set, then less its median over the genes (FIX 5:
    only categories of the centre set get a contrast)."""
    L = {k: np.log(means[k] + pseudo) for k in centre}
    mid = np.mean([L[k] for k in centre], axis=0)
    out = {}
    for k in centre:
        c = L[k] - mid
        out[k] = c - np.median(c)
    return out


def residual(c_src: np.ndarray, c_ref: np.ndarray, guard: np.ndarray) -> tuple[np.ndarray, dict]:
    """Least squares of the source's contrast on the reference's over the guarded genes: c_src = alpha + b c_ref + r.
    Returns r (NaN off the guard) and the fit (alpha, b, n): b is printed, never applied."""
    g = np.asarray(guard, bool)
    r = np.full(len(c_src), np.nan)
    if g.sum() < 3:
        return r, {"alpha": float("nan"), "b": float("nan"), "n": int(g.sum())}
    X = np.vstack([np.ones(int(g.sum())), c_ref[g]]).T
    coef, *_ = np.linalg.lstsq(X, c_src[g], rcond=None)
    r[g] = c_src[g] - X @ coef
    return r, {"alpha": float(coef[0]), "b": float(coef[1]), "n": int(g.sum())}


def a_of(t_target: float, t_ref: float, t_src: float) -> float:
    """The chord's share of the way from the reference to the source, in days."""
    return (t_target - t_ref) / (t_src - t_ref)


def study_mean(vals: list[np.ndarray]) -> np.ndarray:
    """The mean over a study's stages where held (NaN where no stage holds the entry)."""
    V = np.vstack(vals)
    held = np.isfinite(V)
    s = np.where(held, V, 0.0).sum(0)
    n = held.sum(0)
    return np.divide(s, n, out=np.full(V.shape[1], np.nan), where=n > 0)


def consensus(P_by_study: dict[str, np.ndarray], excluded: np.ndarray, clip: float = CLIP) -> tuple[np.ndarray, dict]:
    """FIX 2: the median over studies, kept where at least two studies hold it and all share its sign, never on the
    excluded (typing) genes; clipped at |P| <= clip. NaN = none."""
    V = np.vstack(list(P_by_study.values()))
    held = np.isfinite(V)
    n = held.sum(0)
    sgn = np.sign(np.where(held, V, 0.0))
    pos, neg = (sgn > 0).sum(0), (sgn < 0).sum(0)
    agree = ((pos == n) | (neg == n)) & (n >= 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # an entry no study holds: its median is NaN, kept as none
        med = np.nanmedian(np.where(held, V, np.nan), axis=0) if V.shape[0] else np.full(V.shape[1], np.nan)
    keep = agree & ~np.asarray(excluded, bool) & np.isfinite(med)
    P = np.where(keep, med, np.nan)
    clipped = int((np.abs(np.nan_to_num(P)) > clip).sum())
    P = np.where(keep, np.clip(P, -clip, clip), np.nan)
    return P, {"kept": int(keep.sum()), "clipped": clipped, "held_by_two": int((n >= 2).sum())}


def level_bins(level: np.ndarray, bins: int = PERM_BINS) -> list[np.ndarray]:
    return [b for b in np.array_split(np.argsort(level, kind="stable"), bins) if len(b)]


def spearman_test(ra: np.ndarray, rb: np.ndarray, level: np.ndarray, key: str, n_perm: int = N_PERM,
                  seed: int = PERM_SEED, bins: int = PERM_BINS) -> dict:
    """G0(c) for one lineage and one pair of studies: Spearman of the residuals over the genes both hold, and the 95th
    percentile of n_perm permutations of rb within level bins (seeded by seed and the lineage's name)."""
    both = np.isfinite(ra) & np.isfinite(rb) & np.isfinite(level)
    n = int(both.sum())
    if n < 10:
        return {"n": n, "rho": float("nan"), "p95": float("nan"), "pass": False}
    a, b, lv = ra[both], rb[both], level[both]
    rho = float(stats.spearmanr(a, b).statistic)
    rng = np.random.default_rng([seed, zlib.crc32(key.encode())])
    groups = level_bins(lv, bins)
    ar = stats.rankdata(a)
    perm = np.empty(n_perm)
    for i in range(n_perm):
        bp = b.copy()
        for gidx in groups:
            bp[gidx] = b[rng.permutation(gidx)]
        perm[i] = stats.pearsonr(ar, stats.rankdata(bp)).statistic
    p95 = float(np.percentile(perm, 95))
    return {"n": n, "rho": rho, "p95": p95, "pass": bool(rho >= SPEARMAN_MIN and rho > p95)}


def derange_within_bins(level: np.ndarray, key: str, seed: int, bins: int = PERM_BINS) -> np.ndarray:
    """A permutation of the entries that moves each to another entry of its level bin (equal-count bins by rank, at
    most `bins`, each of at least two entries), never itself: a cyclic shift of a seeded shuffle within each bin."""
    level = np.asarray(level, np.float64)
    n = len(level)
    nb = int(max(1, min(bins, n // 2)))
    rng = np.random.default_rng([seed, zlib.crc32(key.encode())])
    perm = np.arange(n)
    for b in np.array_split(np.argsort(level, kind="stable"), nb):
        if len(b) < 2:
            continue
        o = rng.permutation(b)
        perm[o] = np.roll(o, -1)
    return perm
