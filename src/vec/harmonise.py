"""Map cells from any MERFISH file into a reference file's cell-type vocabulary, and compare compositions across files.

Each released file has its own labels (the Mab21l2 KO says vCM, LPM, Paraxial Mesoderm, Branchial Arch; WT E9.5 says
V-CM, A-CM, Peri, Dorsal-PAM, aPHM...), so comparing compositions across files needs relabelling on the shared panel
genes. The transfer model is

    counts on the shared genes -> log1p CP10k -> per-gene scaling -> PCA -> classifier

with the scaling and the PCA fitted on the reference only. The classifier is multinomial logistic regression
("logreg") or a k-nearest-neighbour vote ("knn"); `classifier="auto"` picks the one with the higher cross-validated
macro-F1 on the reference. It predicts the reference's fine labels. Three levels never contradict each other: lineage
groups (configs/lineage_map.yaml, including the per-file overrides `vec.prep.load_stage` writes to `obs["group"]`)
are sums of fine-label probabilities, and the coarser merged level (the map's `merged` entries, "Untracked" for
groups without a stable cross-stage meaning) sums group probabilities. Cells whose top probability is below
`threshold` are "Unassigned".

Composition across files: quote top call and EM together (`composition_bracket`, `bracket_difference`). Counting top
calls is the best estimator only when the classifier has seen the query's distribution (within a file, estimated vs
true change slope about 1.0). Trained on one z-half of a stack and applied to the other half, the analogue of
wild-type-trained calls on a knockout, top call understated label-level changes by about a quarter while EM was
close to unbiased (experiments/harmonise, stage `shifted`). So the top-call/EM bracket is the default for any query
from another file (a knockout, another stage, another embryo); top call alone is for in-distribution, within-file
use. Cluster-bootstrap intervals resample sections within each population only: with one embryo per condition they
measure section sampling, not embryo-to-embryo or genotype variation.

Sections: `section_ids` and `z_levels` snap z to section levels with `vec.spatial.snap_z` (float noise, levels off
integer values such as E8.75's -240.473939 + 14k, missing levels). `trim_mask` and `trim_sensitivity` drop end
sections of a stack; `z_profile` and `z_location` show where each category sits along it, to check that two stacks
cover the same anatomy before comparing them.

`novelty` is the larger of two ratios, each scaled by the reference's own 99th percentile: the mean distance to the k
nearest reference cells in PCA space and the distance from the reference PCA subspace. Above 1, a cell is more unusual
than 99% of reference cells on at least one of the two. It does not detect cell types missing from the reference on
real data (experiments/harmonise). For cells that are in the model's own training set the kNN part is biased low (the
cell's own record is one of its neighbours) unless `predict(..., exclude_self=True)` removes it.

Depth: median counts per cell are 75 at E9.5, 95 in the Mab21l2 KO and 130 at E8.75. When the reference is deeper,
thin it to the query's depth with `fit_transfer(depth=median_depth(query))` (E8.75 -> E9.5 group accuracy +2
points, experiments/harmonise). When the query is deeper, `predict` thins it binomially to the reference depth and
averages `depth_draws=5` independent draws (one draw lost about 2 points of group accuracy; five matched no thinning
while removing the technical shift seen by `novelty`).

    model = fit_transfer(wt, "celltype")
    calls = predict(model, ko)                                   # one row per query cell
    d = bracket_difference(oof_wt, calls, model, "group", section_ids(wt_kept), section_ids(ko))
"""

from __future__ import annotations

import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, log_loss
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler

from vec.io import COORD_KEY
from vec.prep import lineage_map, log1p_cp10k
from vec.spatial import ZLevels, snap_z

UNASSIGNED = "Unassigned"
UNTRACKED = "Untracked"
METHODS = ("logreg", "knn")
LEVELS = ("label", "group", "merged")
BRACKET = ("top_call", "em")  # the default pair for cross-file compositions
ESTIMATORS = ("top_call", "em", "mean_prob")


# ---------------------------------------------------------------------------------------------------------------------
# Expression on a gene panel


def shared_genes(*adatas: ad.AnnData) -> list[str]:
    """Genes present in every AnnData, in the order of the first."""
    if not adatas:
        raise ValueError("need at least one AnnData")
    common = set(adatas[0].var_names)
    for a in adatas[1:]:
        common &= set(a.var_names)
    return [g for g in adatas[0].var_names if g in common]


def _column_index(a: ad.AnnData, genes: list[str]) -> np.ndarray:
    idx = a.var_names.get_indexer(list(genes))
    if (idx < 0).any():
        missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
        raise KeyError(f"{len(missing)} model genes are missing from the AnnData (e.g. {missing[:5]}); "
                       "fit the model on shared_genes(reference, query)")
    return idx


def panel_counts(a: ad.AnnData, genes: list[str]) -> sparse.csr_matrix | None:
    """Integer counts on `genes` (CSR, in that order), or None when the AnnData has no counts layer."""
    if "counts" not in a.layers:
        return None
    return sparse.csr_matrix(a.layers["counts"])[:, _column_index(a, genes)]


def expression(a: ad.AnnData, genes: list[str], counts=None) -> np.ndarray:
    """Dense float32 log1p CP10k on `genes`.

    Renormalised from counts on `genes` when counts are given or the AnnData has a counts layer (identical to `.X`
    when `genes` is the file's whole panel); otherwise the matching columns of `.X` as stored.
    """
    if counts is None:
        counts = panel_counts(a, genes)
    if counts is not None:
        return log1p_cp10k(counts).toarray()
    X = a.X[:, _column_index(a, genes)]
    return (X.toarray() if sparse.issparse(X) else np.asarray(X)).astype(np.float32)


def thin_counts(counts, frac: float, seed: int = 0) -> sparse.csr_matrix:
    """Keep each molecule independently with probability `frac` (binomial thinning); int32 CSR."""
    if not 0.0 < frac <= 1.0:
        raise ValueError(f"frac must be in (0, 1], got {frac}")
    C = sparse.csr_matrix(counts, dtype=np.int64, copy=True)
    if frac < 1.0:
        C.data = np.random.default_rng(seed).binomial(C.data, frac)
        C.eliminate_zeros()
    return C.astype(np.int32)


def match_depth(counts, target: float, seed: int = 0, tol: float = 0.05) -> tuple[sparse.csr_matrix, float]:
    """Thin `counts` so that its median counts per cell is about `target`.

    No-op unless the median exceeds the target by more than `tol` (counts cannot be added to a shallower file).
    Returns (counts, fraction of molecules kept).
    """
    C = sparse.csr_matrix(counts)
    med = float(np.median(np.asarray(C.sum(1)).ravel()))
    if med == 0 or med <= target * (1.0 + tol):
        return C, 1.0
    frac = float(target / med)
    return thin_counts(C, frac, seed), frac


def median_depth(a: ad.AnnData, genes: list[str] | None = None) -> float | None:
    """Median counts per cell on `genes` (default: all genes); None without a counts layer."""
    counts = panel_counts(a, list(a.var_names) if genes is None else list(genes))
    return None if counts is None else float(np.median(np.asarray(counts.sum(1)).ravel()))


# ---------------------------------------------------------------------------------------------------------------------
# Sections along z


def z_levels(a: ad.AnnData, tol: float | None = None) -> ZLevels:
    """Section levels of `obsm['spatial_3D'][:, 2]` from `vec.spatial.snap_z`: per-cell level index (0 = lowest z),
    ascending level values and the stack's step."""
    return snap_z(np.asarray(a.obsm[COORD_KEY])[:, 2], tol=tol)


def section_ids(a: ad.AnnData, tol: float | None = None) -> np.ndarray:
    """Section label per cell, "z<level>", with z snapped to section levels by `vec.spatial.snap_z`.

    z is a section index times a per-stage unit and carries float noise (WT E9.5 stores 98 raw values on 40 levels);
    levels need not sit on integers (E8.75's are -240.473939 + 14k), so nothing is rounded against a fixed origin.
    The label shows the level value rounded to 6 decimals; `z_levels(a).index` gives the rank along z.
    """
    zl = z_levels(a, tol)
    names = np.array([f"z{v:.12g}" for v in np.round(zl.levels, 6) + 0.0])
    return names[zl.index] if names.size else np.array([], dtype=str)


def subset_levels(zl: ZLevels, rows) -> ZLevels:
    """The same stack restricted to some cells (`rows`: positions or a boolean mask); levels and step are kept."""
    return ZLevels(index=np.asarray(zl.index)[rows], levels=zl.levels, step=zl.step)


def z_position(zl: ZLevels) -> np.ndarray:
    """Relative position of each cell's level along the stack by z value: 0 at the lowest level, 1 at the highest
    (missing levels count as distance)."""
    if zl.levels.size < 2:
        return np.zeros(np.asarray(zl.index).size)
    return (zl.levels[zl.index] - zl.levels[0]) / (zl.levels[-1] - zl.levels[0])


def trim_mask(zl: ZLevels, low: int = 0, high: int = 0) -> np.ndarray:
    """Cells kept after dropping the `low` lowest and `high` highest section levels of the stack."""
    n = zl.levels.size
    idx = np.asarray(zl.index)
    return (idx >= low) & (idx <= n - 1 - high)


def default_trims(a: str = "a", b: str = "b", k: int = 6, both=(2, 4, 6)) -> dict[str, tuple[int, int, int, int]]:
    """The 8 z-trims of the 2026-10-19 review, as {name: (a low, a high, b low, b high)} numbers of end sections
    dropped: none, both ends of both stacks by each of `both`, and each end of each stack alone by `k`."""
    trims = {"none": (0, 0, 0, 0)}
    trims.update({f"both ends -{j}": (j, j, j, j) for j in both})
    trims.update({f"{a} low -{k}": (k, 0, 0, 0), f"{a} high -{k}": (0, k, 0, 0),
                  f"{b} low -{k}": (0, 0, k, 0), f"{b} high -{k}": (0, 0, 0, k)})
    return trims


# ---------------------------------------------------------------------------------------------------------------------
# Lineage vocabularies


def group_membership(classes, groups, group_of: dict[str, str]) -> np.ndarray:
    """(n_classes, n_groups) 0/1 matrix: column j marks the labels whose group is groups[j]."""
    own = np.array([group_of[c] for c in classes])
    return (own[:, None] == np.asarray(groups)[None, :]).astype(float)


def group_lookup(labels, group_map: dict[str, str] | None = None) -> dict[str, str]:
    """Lineage group of each label (configs/lineage_map.yaml unless `group_map` is given); unmapped -> "Unknown"."""
    gm = dict(lineage_map()["labels"]) if group_map is None else dict(group_map)
    return {str(lab): gm.get(str(lab), "Unknown") for lab in labels}


def lineage_groups_of(a: ad.AnnData, label_key: str = "celltype", group_key: str = "group") -> dict[str, str] | None:
    """label -> group as annotated in `a.obs[group_key]` (vec.prep.load_stage applies the lineage map's per-file
    overrides there), or None when a column is missing or a label carries more than one group."""
    if label_key not in a.obs or group_key not in a.obs:
        return None
    pairs = pd.DataFrame({"label": a.obs[label_key].astype(str).to_numpy(),
                          "group": a.obs[group_key].astype(str).to_numpy()}).drop_duplicates()
    if pairs["label"].duplicated().any():
        return None
    return dict(zip(pairs["label"], pairs["group"], strict=True))


def group_order(groups) -> np.ndarray:
    """Groups in lineage-map order, then any others alphabetically."""
    present = {str(g) for g in groups}
    known = [g for g in lineage_map()["groups"] if g in present]
    return np.array(known + sorted(present - set(known)), dtype=object)


def merged_lookup(groups) -> dict[str, str]:
    """Merged category of each lineage group (the map's `merged` entry); "Untracked" for null or unknown groups."""
    spec = lineage_map()["groups"]
    return {str(g): (spec.get(str(g)) or {}).get("merged") or UNTRACKED for g in groups}


def merged_order(merged) -> np.ndarray:
    """Merged categories in lineage-map order, others alphabetically, "Untracked" last."""
    present = {str(m) for m in merged}
    known: list[str] = []
    for g in lineage_map()["groups"].values():
        m = g.get("merged")
        if m and m in present and m not in known:
            known.append(m)
    rest = sorted(present - set(known) - {UNTRACKED})
    return np.array(known + rest + ([UNTRACKED] if UNTRACKED in present else []), dtype=object)


# ---------------------------------------------------------------------------------------------------------------------
# Model


def _clip(X, bound):
    return np.clip(X, -bound, bound)


def make_pipeline(method: str = "logreg", n_pcs: int = 50, k: int = 15, C: float = 1.0, clip: float = 10.0,
                  seed: int = 0, n_jobs: int = 1) -> Pipeline:
    """Per-gene scaling (clipped at +-`clip` SD) -> PCA -> classifier. The first three steps are the embedding."""
    steps = [
        ("scale", StandardScaler()),
        ("clip", FunctionTransformer(_clip, kw_args={"bound": clip})),
        ("pca", PCA(n_components=n_pcs, random_state=seed)),
    ]
    if method == "logreg":
        steps += [("pcscale", StandardScaler()), ("clf", LogisticRegression(C=C, max_iter=2000))]
    elif method == "knn":
        steps += [("clf", KNeighborsClassifier(n_neighbors=k, n_jobs=n_jobs))]
    else:
        raise ValueError(f"unknown classifier {method!r}; choose from {METHODS}")
    return Pipeline(steps)


@dataclass
class TransferModel:
    """A fitted label-transfer model. Probabilities are over `classes`; `groups` are their lineage groups and `merged`
    the merged categories of those groups. `train_obs_names` are the reference cells it was fitted on."""

    genes: list[str]
    label_key: str
    method: str
    classes: np.ndarray
    groups: np.ndarray
    group_of: dict[str, str]
    pipeline: Pipeline
    priors: np.ndarray
    depth: float | None
    neighbours: NearestNeighbors | None
    novelty_scale: float
    residual_scale: float
    threshold: float = 0.5
    cv: pd.DataFrame | None = None
    dropped: dict[str, str] = field(default_factory=dict)
    params: dict = field(default_factory=dict)
    merged: np.ndarray = field(default_factory=lambda: np.array([], dtype=object))
    merged_of: dict[str, str] = field(default_factory=dict)
    train_obs_names: np.ndarray | None = None

    @property
    def membership(self) -> np.ndarray:
        """(n_classes, n_groups) indicator: which group each fine label belongs to."""
        return group_membership(self.classes, self.groups, self.group_of)

    @property
    def merged_membership(self) -> np.ndarray:
        """(n_groups, n_merged) indicator: which merged category each group belongs to."""
        return group_membership(self.groups, self.merged, self.merged_of)

    def categories(self, level: str) -> np.ndarray:
        """Categories at `level` ("label", "group" or "merged")."""
        _check_level(level)
        return {"label": self.classes, "group": self.groups, "merged": self.merged}[level]

    def projection(self, level: str) -> np.ndarray:
        """(n_classes, n_categories) 0/1 matrix summing fine-label probabilities into `level` categories."""
        _check_level(level)
        if level == "label":
            return np.eye(len(self.classes))
        if level == "group":
            return self.membership
        return self.membership @ self.merged_membership


def _check_level(level: str) -> None:
    if level not in LEVELS:
        raise ValueError(f"level must be one of {LEVELS}, got {level!r}")


def _vocabulary(classes, group_map):
    """(group_of, group order, merged_of, merged order) for fine labels `classes`."""
    group_of = group_lookup(classes, group_map)
    groups = group_order(group_of.values())
    merged_of = merged_lookup(groups)
    return group_of, groups, merged_of, merged_order(merged_of.values())


def _prepare(reference, label_key, genes, exclude_labels, min_class_cells, depth, seed, depth_tol):
    """Training matrix for the kept reference cells: (X, y, keep mask, median depth, dropped labels)."""
    labels = reference.obs[label_key].astype(str).to_numpy()
    keep = ~np.isin(labels, list(exclude_labels))
    dropped = {str(lab): "excluded" for lab in np.unique(labels[~keep])}
    vals, cnt = np.unique(labels[keep], return_counts=True)
    rare = vals[cnt < min_class_cells]
    dropped.update({str(lab): f"fewer than {min_class_cells} cells" for lab in rare})
    keep &= ~np.isin(labels, rare)
    if not keep.any():
        raise ValueError("no reference cells left after excluding labels")
    counts = panel_counts(reference, genes)
    ref_depth = None
    if counts is not None:
        counts = counts[np.flatnonzero(keep)]
        if depth is not None:
            counts, _ = match_depth(counts, depth, seed=seed, tol=depth_tol)
        ref_depth = float(np.median(np.asarray(counts.sum(1)).ravel()))
        X = log1p_cp10k(counts).toarray()
    else:
        X = expression(reference, genes)[keep]
    return X, labels[keep], keep, ref_depth, dropped


def _embed(pipe, X):
    """(scaled and clipped expression, PCA scores, distance from the PCA subspace) under a fitted pipeline."""
    Xs = pipe[:2].transform(X)
    Z = pipe["pca"].transform(Xs)
    return Xs, Z, np.linalg.norm(Xs - pipe["pca"].inverse_transform(Z), axis=1)


def _fit_one(X, y, method, n_pcs, k, C, clip, seed, n_jobs, novelty, k_novelty, novelty_quantile):
    """Fitted pipeline and novelty reference (neighbour index, kNN-distance scale, residual scale) or None."""
    n_pcs = int(max(1, min(n_pcs, X.shape[1] - 1, X.shape[0] - 1)))
    pipe = make_pipeline(method, n_pcs=n_pcs, k=k, C=C, clip=clip, seed=seed, n_jobs=n_jobs).fit(X, y)
    nov = None
    if novelty:
        _, Z, res = _embed(pipe, X)
        kk = int(min(k_novelty, len(X) - 1))
        nn = NearestNeighbors(n_neighbors=kk, n_jobs=n_jobs).fit(Z)
        d, _ = nn.kneighbors(Z, n_neighbors=kk + 1)  # first neighbour is the cell itself
        nov = (nn, float(np.quantile(d[:, 1:].mean(1), novelty_quantile)), float(np.quantile(res, novelty_quantile)))
    return pipe, nov


def _apply(pipe, classes, X):
    """Class probabilities over `classes` (zero for classes the fitted model never saw)."""
    P = pipe.predict_proba(X)
    fitted = np.asarray(pipe.classes_)
    if not np.array_equal(fitted, classes):
        full = np.zeros((X.shape[0], len(classes)))
        full[:, np.searchsorted(classes, fitted)] = P
        P = full
    return P


def _novelty(pipe, nov, X, self_rows=None):
    """(kNN component, residual component) of the novelty score, each over its reference 99th percentile.

    self_rows: per query cell, the row of its own record in the neighbour index (-1 if none); that record is left
    out of the cell's neighbours, as when the scale was set on the reference.
    """
    if nov is None:
        return np.full(X.shape[0], np.nan), np.full(X.shape[0], np.nan)
    nn, knn_scale, res_scale = nov
    _, Z, res = _embed(pipe, X)
    if self_rows is None or not (np.asarray(self_rows) >= 0).any():
        d, _ = nn.kneighbors(Z)
        return d.mean(1) / knn_scale, res / res_scale
    k = nn.n_neighbors
    d, idx = nn.kneighbors(Z, n_neighbors=k + 1)
    own = idx == np.asarray(self_rows)[:, None]
    drop = np.where(own.any(1), own.argmax(1), k)  # the cell's own record, else the farthest of the k + 1
    keep = np.ones(d.shape, dtype=bool)
    keep[np.arange(len(d)), drop] = False
    return d[keep].reshape(len(d), k).mean(1) / knn_scale, res / res_scale


def _frame(P, nov_parts, voc, index, threshold, group_threshold, novelty_threshold, proba):
    classes, groups, M, merged, MM = voc
    nov_knn, nov_res = nov_parts
    nov = np.fmax(nov_knn, nov_res)
    rows = np.arange(len(P))
    i = P.argmax(1)
    conf = P[rows, i]
    PG = P @ M
    j = PG.argmax(1)
    gconf = PG[rows, j]
    PM = PG @ MM
    m = PM.argmax(1)
    mconf = PM[rows, m]
    novel = nov > novelty_threshold if novelty_threshold is not None else np.zeros(len(P), dtype=bool)
    gthr = threshold if group_threshold is None else group_threshold
    label_raw = np.asarray(classes, dtype=object)[i]
    group_raw = np.asarray(groups, dtype=object)[j]
    merged_raw = np.asarray(merged, dtype=object)[m]
    df = pd.DataFrame({
        "label": np.where((conf >= threshold) & ~novel, label_raw, UNASSIGNED),
        "label_raw": label_raw,
        "label_conf": conf,
        "group": np.where((gconf >= gthr) & ~novel, group_raw, UNASSIGNED),
        "group_raw": group_raw,
        "group_conf": gconf,
        "merged": np.where((mconf >= gthr) & ~novel, merged_raw, UNASSIGNED),
        "merged_raw": merged_raw,
        "merged_conf": mconf,
        "novelty": nov,
        "novelty_knn": nov_knn,
        "novelty_residual": nov_res,
    }, index=index)
    if proba:
        df = pd.concat([df, pd.DataFrame(P, index=index, columns=[f"p:{c}" for c in classes]),
                        pd.DataFrame(PG, index=index, columns=[f"pg:{g}" for g in groups]),
                        pd.DataFrame(PM, index=index, columns=[f"pm:{c}" for c in merged])], axis=1)
    return df


def _stratified_subsample(y, n, min_per_class, seed):
    rng = np.random.default_rng(seed)
    frac = n / len(y)
    idx = []
    for c in np.unique(y):
        members = np.flatnonzero(y == c)
        take = min(len(members), max(min_per_class, int(round(len(members) * frac))))
        idx.append(rng.choice(members, take, replace=False))
    return np.sort(np.concatenate(idx))


def _cv_arrays(X, y, index, *, method, n_folds, groups, seed, n_pcs, k, C, clip, n_jobs, novelty, k_novelty,
               novelty_quantile, threshold, group_map):
    classes = np.unique(y)
    group_of, gorder, merged_of, morder = _vocabulary(classes, group_map)
    voc = (classes, gorder, group_membership(classes, gorder, group_of), morder,
           group_membership(gorder, morder, merged_of))
    if groups is None:
        splits = StratifiedKFold(n_folds, shuffle=True, random_state=seed).split(X, y)
    else:
        splits = StratifiedGroupKFold(n_folds, shuffle=True, random_state=seed).split(X, y, groups)
    P = np.zeros((len(y), len(classes)))
    nov_knn, nov_res = np.full(len(y), np.nan), np.full(len(y), np.nan)
    fold = np.full(len(y), -1)
    for f, (tr, te) in enumerate(splits):
        pipe, nov = _fit_one(X[tr], y[tr], method, n_pcs, k, C, clip, seed, n_jobs, novelty, k_novelty,
                             novelty_quantile)
        P[te] = _apply(pipe, classes, X[te])
        nov_knn[te], nov_res[te] = _novelty(pipe, nov, X[te])
        fold[te] = f
    df = _frame(P, (nov_knn, nov_res), voc, index, threshold, None, None, proba=True)
    df.insert(0, "fold", fold)
    true_group = [group_of[lab] for lab in y]
    df.insert(0, "true_merged", [merged_of[g] for g in true_group])
    df.insert(0, "true_group", true_group)
    df.insert(0, "true", y)
    return df


def cross_validate(reference: ad.AnnData, label_key: str = "celltype", *, method: str = "logreg",
                   genes: list[str] | None = None, n_folds: int = 5, groups=None, n_pcs: int = 50, k: int = 15,
                   C: float = 1.0, clip: float = 10.0, threshold: float = 0.5,
                   exclude_labels=("Unknown",), min_class_cells: int = 10, group_map: dict[str, str] | None = None,
                   depth: float | None = None, depth_tol: float = 0.05, max_cells: int | None = None,
                   novelty: bool = True, k_novelty: int = 15, novelty_quantile: float = 0.99, seed: int = 0,
                   n_jobs: int = 1) -> pd.DataFrame:
    """Out-of-fold predictions on the reference itself, one row per kept cell.

    Same columns as `predict` plus `true`, `true_group`, `true_merged` and `fold`. Folds are stratified by label; with
    `groups` (one id per reference cell, e.g. `section_ids`) whole groups are held out together
    (StratifiedGroupKFold). `max_cells` evaluates on a stratified subsample. Scaling, PCA and classifier are refitted
    inside every fold. group_map defaults to the reference's own `obs["group"]` pairs (see `lineage_groups_of`), else
    the lineage map. `depth` thins the reference first, as in `fit_transfer`.
    """
    genes = list(reference.var_names) if genes is None else list(genes)
    X, y, keep, _, dropped = _prepare(reference, label_key, genes, exclude_labels, max(min_class_cells, n_folds),
                                      depth, seed, depth_tol)
    index = reference.obs_names[keep]
    grp = None if groups is None else np.asarray(groups)[keep]
    if max_cells is not None and max_cells < len(y):
        sub = _stratified_subsample(y, max_cells, n_folds, seed)
        X, y, index = X[sub], y[sub], index[sub]
        grp = None if grp is None else grp[sub]
    gm = group_map if group_map is not None else lineage_groups_of(reference, label_key)
    df = _cv_arrays(X, y, index, method=method, n_folds=n_folds, groups=grp, seed=seed, n_pcs=n_pcs, k=k, C=C,
                    clip=clip, n_jobs=n_jobs, novelty=novelty, k_novelty=k_novelty,
                    novelty_quantile=novelty_quantile, threshold=threshold, group_map=gm)
    df.attrs.update(method=method, dropped=dropped, grouped=groups is not None)
    return df


def fit_transfer(reference: ad.AnnData, label_key: str = "celltype", *, genes: list[str] | None = None,
                 classifier: str = "auto", n_pcs: int = 50, k: int = 15, C: float = 1.0, clip: float = 10.0,
                 threshold: float = 0.5, cv_folds: int = 5, cv_groups=None, cv_max_cells: int | None = None,
                 exclude_labels=("Unknown",), min_class_cells: int = 10, group_map: dict[str, str] | None = None,
                 depth: float | None = None, depth_tol: float = 0.05, k_novelty: int = 15,
                 novelty_quantile: float = 0.99, seed: int = 0, n_jobs: int = 1) -> TransferModel:
    """Fit a label-transfer model on a reference AnnData.

    reference: log1p CP10k `.X` and, ideally, an integer `layers["counts"]` (as from `vec.prep.load_stage`).
    label_key: obs column with the fine labels. genes: panel to use (default: all reference genes; pass
    `shared_genes(reference, query)` when the query lacks some). classifier: "logreg", "knn" or "auto" (stratified
    `cv_folds`-fold CV on the reference, higher macro-F1 wins; the table is kept in `model.cv`). exclude_labels:
    labels left out of training ("Unknown" is not a cell type); labels with fewer than `min_class_cells` cells are
    dropped too (both recorded in `model.dropped`). group_map: label -> coarse group (default: the reference's own
    `obs["group"]`, which carries the lineage map's per-file overrides, else the lineage map). depth: thin the
    reference to this median counts per cell first (use the query's depth when it is shallower).
    """
    if classifier not in (*METHODS, "auto"):
        raise ValueError(f"classifier must be one of {(*METHODS, 'auto')}, got {classifier!r}")
    genes = list(reference.var_names) if genes is None else list(genes)
    X, y, keep, ref_depth, dropped = _prepare(reference, label_key, genes, exclude_labels,
                                              max(min_class_cells, cv_folds), depth, seed, depth_tol)
    gm = group_map if group_map is not None else lineage_groups_of(reference, label_key)
    params = dict(n_pcs=n_pcs, k=k, C=C, clip=clip, k_novelty=k_novelty, novelty_quantile=novelty_quantile,
                  seed=seed, exclude_labels=list(exclude_labels), min_class_cells=min_class_cells,
                  depth_target=depth, n_train=int(len(y)))
    cv, method = None, classifier
    if classifier == "auto":
        index = reference.obs_names[keep]
        grp = None if cv_groups is None else np.asarray(cv_groups)[keep]
        Xc, yc = X, y
        if cv_max_cells is not None and cv_max_cells < len(y):
            sub = _stratified_subsample(y, cv_max_cells, cv_folds, seed)
            Xc, yc, index = X[sub], y[sub], index[sub]
            grp = None if grp is None else grp[sub]
        rows = []
        for m in METHODS:
            t0 = time.perf_counter()
            oof = _cv_arrays(Xc, yc, index, method=m, n_folds=cv_folds, groups=grp, seed=seed, n_pcs=n_pcs, k=k,
                             C=C, clip=clip, n_jobs=n_jobs, novelty=False, k_novelty=k_novelty,
                             novelty_quantile=novelty_quantile, threshold=threshold, group_map=gm)
            rows.append({"method": m, **transfer_metrics(oof), "seconds": time.perf_counter() - t0})
        cv = pd.DataFrame(rows)
        method = str(cv.loc[cv["macro_f1"].idxmax(), "method"])
    pipe, (nn, knn_scale, res_scale) = _fit_one(X, y, method, n_pcs, k, C, clip, seed, n_jobs, True, k_novelty,
                                                novelty_quantile)
    classes = np.asarray(pipe.classes_)
    group_of, groups, merged_of, merged = _vocabulary(classes, gm)
    return TransferModel(
        genes=genes, label_key=label_key, method=method, classes=classes, groups=groups, group_of=group_of,
        pipeline=pipe, priors=np.array([(y == c).mean() for c in classes]), depth=ref_depth, neighbours=nn,
        novelty_scale=knn_scale, residual_scale=res_scale, threshold=threshold, cv=cv, dropped=dropped,
        params=params, merged=merged, merged_of=merged_of,
        train_obs_names=np.asarray(reference.obs_names[keep]).astype(str),
    )


def predict(model: TransferModel, query: ad.AnnData, *, threshold: float | None = None,
            group_threshold: float | None = None, novelty_threshold: float | None = None, depth_match: bool = True,
            depth_draws: int = 5, depth_tol: float = 0.05, seed: int = 0, proba: bool = True,
            exclude_self: bool = False) -> pd.DataFrame:
    """Map query cells into the model's vocabulary; one row per query cell (index = query obs names).

    Columns: `label` (fine label, or "Unassigned" when its probability is below `threshold`, default the model's),
    `label_raw` (top label regardless), `label_conf`; `group`, `group_raw`, `group_conf` for the lineage groups and
    `merged`, `merged_raw`, `merged_conf` for the merged level (threshold `group_threshold`, default `threshold`);
    `novelty` (max of `novelty_knn` and `novelty_residual`, see the module docstring); and, with `proba`,
    probabilities `p:<label>`, `pg:<group>` and `pm:<merged>`. `novelty_threshold` also unassigns cells whose novelty
    exceeds it. With `depth_match` and a counts layer, a query deeper than the reference (median counts per cell) is
    thinned to the reference depth first; `df.attrs["depth_frac"]` records the fraction of molecules kept.
    Probabilities and novelty are averaged over `depth_draws` independent thinning draws.

    Novelty of cells the model was trained on is biased low, because each one's own record sits among its nearest
    reference neighbours. `exclude_self=True` drops that record for query cells whose obs name is a training cell
    (use it when the query is the model's own reference). It stays opt-in because obs names can repeat across files;
    a warning is raised when most query names are training names.
    """
    self_rows = None
    names = None if model.train_obs_names is None else pd.Index(model.train_obs_names)
    if names is not None and names.is_unique:
        pos = names.get_indexer(np.asarray(query.obs_names).astype(str))
        if exclude_self:
            self_rows = pos
        elif query.n_obs and (pos >= 0).mean() > 0.5:
            warnings.warn(f"{(pos >= 0).mean():.0%} of query cells have the obs names of the model's training cells; "
                          "their novelty is biased low unless predict(..., exclude_self=True)", stacklevel=2)
    elif exclude_self:
        raise ValueError("exclude_self needs unique model.train_obs_names (refit the model with this version)")
    counts = panel_counts(query, model.genes)
    frac = 1.0
    draws = [counts]
    if depth_match and counts is not None and model.depth is not None:
        thinned, frac = match_depth(counts, model.depth, seed=seed, tol=depth_tol)
        draws = [thinned]
        if frac < 1.0:
            draws += [thin_counts(counts, frac, seed + i) for i in range(1, max(1, depth_draws))]
    P, nov_knn, nov_res = 0.0, 0.0, 0.0
    novelty_ref = (model.neighbours, model.novelty_scale, model.residual_scale)
    for c in draws:
        X = expression(query, model.genes, counts=c)
        P = P + _apply(model.pipeline, model.classes, X) / len(draws)
        nk, nr = _novelty(model.pipeline, novelty_ref, X, self_rows)
        nov_knn, nov_res = nov_knn + nk / len(draws), nov_res + nr / len(draws)
    thr = model.threshold if threshold is None else threshold
    voc = (model.classes, model.groups, model.membership, model.merged, model.merged_membership)
    df = _frame(P, (nov_knn, nov_res), voc, query.obs_names, thr, group_threshold, novelty_threshold, proba)
    df.attrs.update(depth_frac=frac, depth_draws=len(draws), method=model.method, exclude_self=bool(exclude_self))
    return df


# ---------------------------------------------------------------------------------------------------------------------
# Evaluation


def transfer_metrics(pred: pd.DataFrame, true_col: str = "true", true_group_col: str = "true_group",
                     true_merged_col: str = "true_merged") -> dict:
    """Accuracy and macro-F1 of the top calls (fine, group and, when a truth column exists, merged level), coverage
    and accuracy of the thresholded calls, and log-loss when probabilities are present. Macro-F1 averages over the
    categories present in the truth."""
    y = pred[true_col].astype(str).to_numpy()
    yhat = pred["label_raw"].astype(str).to_numpy()
    assigned = pred["label"].to_numpy() != UNASSIGNED
    out = {
        "n": int(len(pred)),
        "accuracy": float(accuracy_score(y, yhat)),
        "macro_f1": float(f1_score(y, yhat, labels=np.unique(y), average="macro", zero_division=0)),
        "coverage": float(assigned.mean()),
        "assigned_accuracy": float((yhat[assigned] == y[assigned]).mean()) if assigned.any() else float("nan"),
    }
    pcols = [c for c in pred.columns if c.startswith("p:")]
    if pcols:
        classes = np.array([c[2:] for c in pcols])
        inside = np.isin(y, classes)
        if inside.any():
            P = pred[pcols].to_numpy()[inside]
            out["log_loss"] = float(log_loss(y[inside], P / P.sum(1, keepdims=True), labels=classes))
    if true_group_col in pred:
        g = pred[true_group_col].astype(str).to_numpy()
        ghat = pred["group_raw"].astype(str).to_numpy()
        gassigned = pred["group"].to_numpy() != UNASSIGNED
        out.update(
            group_accuracy=float(accuracy_score(g, ghat)),
            group_macro_f1=float(f1_score(g, ghat, labels=np.unique(g), average="macro", zero_division=0)),
            group_coverage=float(gassigned.mean()),
            group_assigned_accuracy=float((ghat[gassigned] == g[gassigned]).mean()) if gassigned.any()
            else float("nan"),
        )
    if true_merged_col in pred and "merged_raw" in pred:
        mt = pred[true_merged_col].astype(str).to_numpy()
        mhat = pred["merged_raw"].astype(str).to_numpy()
        out.update(merged_accuracy=float(accuracy_score(mt, mhat)),
                   merged_macro_f1=float(f1_score(mt, mhat, labels=np.unique(mt), average="macro",
                                                  zero_division=0)))
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Composition estimators


def em_priors(P, train_priors, weights=None, init=None, max_iter: int = 1000,
              tol: float = 1e-8) -> tuple[np.ndarray, np.ndarray]:
    """Class priors of a new population from posteriors computed under the training priors (Saerens, Latinne and
    Decaestecker 2002). Corrects classify-and-count's pull toward the training composition under prior shift.

    weights: optional per-cell weights (e.g. bootstrap multiplicities); init: starting priors (warm start).
    Returns (priors, adjusted posteriors).
    """
    P = np.asarray(P, dtype=float)
    w0 = np.asarray(train_priors, dtype=float)
    if (w0 <= 0).any():
        raise ValueError("training priors must be positive")
    cw = np.full(len(P), 1.0 / len(P)) if weights is None else np.asarray(weights, dtype=float)
    cw = cw / cw.sum()
    q = w0.copy() if init is None else np.asarray(init, dtype=float).copy()
    W = P
    for _ in range(max_iter):
        W = P * (q / w0)
        W /= np.clip(W.sum(1, keepdims=True), 1e-300, None)
        q_new = cw @ W
        done = np.abs(q_new - q).max() < tol
        q = q_new
        if done:
            break
    return q, W


def _em_block(P, w0, Wb, Q, max_iter, tol):
    """EM for a block of populations at once. Wb: (n_cells, b) column-normalised weights; Q: (K, b) starting priors.

    Same update and stopping rule as `em_priors`, per column: q <- (q / w0) * P' (w / (P (q / w0))).
    """
    Q = np.array(Q, dtype=float)
    active = np.arange(Wb.shape[1])
    for _ in range(max_iter):
        R = Q[:, active] / w0[:, None]
        S = np.clip(P @ R, 1e-300, None)
        Qn = R * (P.T @ (Wb[:, active] / S))
        done = np.abs(Qn - Q[:, active]).max(0) < tol
        Q[:, active] = Qn
        active = active[~done]
        if not active.size:
            break
    return Q


def em_priors_many(P, train_priors, weights, init=None, max_iter: int = 1000, tol: float = 1e-8,
                   block: int = 256) -> np.ndarray:
    """`em_priors` for many weightings of the same cells at once (e.g. bootstrap resamples), vectorised.

    weights: (n_cells, B), one column per population (normalised internally). init: (K,) or (B, K) starting priors.
    Returns (B, K) priors; column b matches `em_priors(P, train_priors, weights[:, b], init)` to within `tol`.
    """
    P = np.asarray(P, dtype=float)
    w0 = np.asarray(train_priors, dtype=float)
    if (w0 <= 0).any():
        raise ValueError("training priors must be positive")
    W = np.asarray(weights, dtype=float)
    if W.ndim != 2 or W.shape[0] != P.shape[0]:
        raise ValueError(f"weights must be (n_cells, B) with n_cells = {P.shape[0]}, got {W.shape}")
    B, K = W.shape[1], P.shape[1]
    Q0 = np.broadcast_to(w0 if init is None else np.asarray(init, dtype=float), (B, K))
    out = np.empty((B, K))
    for s in range(0, B, block):
        Wb = W[:, s:s + block]
        out[s:s + block] = _em_block(P, w0, Wb / Wb.sum(0, keepdims=True), Q0[s:s + block].T, max_iter, tol).T
    return out


def _resample(clusters, n_boot: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """(cluster index per cell, (n_boot, n_clusters) multiplicities): whole clusters drawn with replacement."""
    ids, inv = np.unique(np.asarray(clusters), return_inverse=True)
    mult = np.random.default_rng(seed).multinomial(len(ids), np.full(len(ids), 1.0 / len(ids)), size=n_boot)
    return inv.ravel(), mult


def bootstrap_em(P, train_priors, clusters, n_boot: int = 2000, seed: int = 0, tol: float = 1e-7,
                 max_iter: int = 1000, block: int = 256, n_jobs: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """EM priors on all cells and on `n_boot` cluster-bootstrap resamples (warm-started from the full estimate).

    The resamples are those of `bootstrap_draws` with the same `clusters`, `n_boot` and `seed`, so EM and top-call
    draws pair up. Blocks of `block` resamples run vectorised, `n_jobs` blocks at a time in threads (the result does
    not depend on `n_jobs`). Returns (priors on all cells, (n_boot, K) draws).
    """
    P = np.asarray(P, dtype=float)
    w0 = np.asarray(train_priors, dtype=float)
    q0, _ = em_priors(P, w0, max_iter=max_iter)
    inv, mult = _resample(clusters, n_boot, seed)
    D = np.empty((n_boot, P.shape[1]))

    def run(s):
        Wb = mult[s:s + block][:, inv].T.astype(float)
        Wb /= Wb.sum(0, keepdims=True)
        D[s:s + Wb.shape[1]] = _em_block(P, w0, Wb, np.repeat(q0[:, None], Wb.shape[1], 1), max_iter, tol).T

    starts = range(0, n_boot, block)
    if n_jobs > 1 and len(starts) > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(starts))) as ex:
            list(ex.map(run, starts))
    else:
        for s in starts:
            run(s)
    return q0, D


def composition(pred: pd.DataFrame, level: str = "group", method: str = "hard",
                model: TransferModel | None = None) -> pd.Series:
    """Fraction of cells per category from `predict` output.

    method: "hard" counts the thresholded calls (Unassigned included), "raw" the top calls, "soft" averages the
    probabilities, "em" re-estimates the priors by EM (needs `model`; group and merged fractions sum fine-label
    priors). level: "label", "group" or "merged". For a query from another file use `composition_bracket`.
    """
    _check_level(level)
    if method in ("hard", "raw"):
        col = level if method == "hard" else f"{level}_raw"
        return pred[col].value_counts(normalize=True)
    if method == "soft":
        prefix = {"group": "pg:", "label": "p:", "merged": "pm:"}[level]
        cols = [c for c in pred.columns if c.startswith(prefix)]
        return pd.Series(pred[cols].to_numpy().mean(0), index=[c[len(prefix):] for c in cols])
    if method == "em":
        if model is None:
            raise ValueError("method='em' needs the model (training priors)")
        q, _ = em_priors(pred[[f"p:{c}" for c in model.classes]].to_numpy(), model.priors)
        return pd.Series(q @ model.projection(level), index=model.categories(level))
    raise ValueError(f"unknown method {method!r}")


def one_hot(labels, categories=None) -> pd.DataFrame:
    """Cells x categories 0/1 frame (labels outside `categories` get an all-zero row)."""
    labels = np.asarray(labels).astype(str)
    cats = list(pd.unique(labels)) if categories is None else [str(c) for c in categories]
    return pd.DataFrame((labels[:, None] == np.asarray(cats)[None, :]).astype(float), columns=cats)


def bootstrap_draws(weights, clusters, n_boot: int = 2000, seed: int = 0) -> np.ndarray:
    """(n_boot, n_categories) compositions from resampling whole clusters (e.g. sections) with replacement.

    weights: cells x categories (one-hot calls or probabilities). Each draw is sum(weights) / number of cells over
    the resampled clusters, so fractions are of all cells.
    """
    W = np.asarray(weights, dtype=float)
    inv, mult = _resample(clusters, n_boot, seed)
    S = np.zeros((mult.shape[1], W.shape[1]))
    np.add.at(S, inv, W)
    n = np.bincount(inv, minlength=mult.shape[1]).astype(float)
    return (mult @ S) / (mult @ n)[:, None]


def bootstrap_composition(weights: pd.DataFrame, clusters, n_boot: int = 2000, seed: int = 0,
                          ci: float = 0.95) -> pd.DataFrame:
    """Composition with cluster-bootstrap percentile intervals: columns estimate, lo, hi, se."""
    W = weights.to_numpy(dtype=float)
    D = bootstrap_draws(W, clusters, n_boot, seed)
    a = (1.0 - ci) / 2.0
    return pd.DataFrame({"estimate": W.mean(0), "lo": np.quantile(D, a, axis=0),
                         "hi": np.quantile(D, 1.0 - a, axis=0), "se": D.std(0, ddof=1)}, index=weights.columns)


def bootstrap_difference(weights_a: pd.DataFrame, clusters_a, weights_b: pd.DataFrame, clusters_b,
                         n_boot: int = 2000, seed: int = 0, ci: float = 0.95,
                         pseudo: float = 1e-3) -> pd.DataFrame:
    """Composition of b minus a, with independent cluster bootstraps of each population.

    Columns: a, b, diff (b - a), diff_lo, diff_hi, log2_ratio (log2((b + pseudo) / (a + pseudo))) with its interval,
    and ci_excludes_0. Categories are aligned by name (missing ones count as 0). The intervals cover resampling of
    clusters within each population only (with sections of one embryo per condition: section sampling, not embryo or
    genotype variation).
    """
    cols = list(weights_a.columns) + [c for c in weights_b.columns if c not in set(weights_a.columns)]
    Wa = weights_a.reindex(columns=cols, fill_value=0.0).to_numpy(dtype=float)
    Wb = weights_b.reindex(columns=cols, fill_value=0.0).to_numpy(dtype=float)
    Da = bootstrap_draws(Wa, clusters_a, n_boot, seed)
    Db = bootstrap_draws(Wb, clusters_b, n_boot, seed + 1)
    ea, eb = Wa.mean(0), Wb.mean(0)
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    diff = Db - Da
    lr = np.log2((Db + pseudo) / (Da + pseudo))
    out = pd.DataFrame({
        "a": ea, "b": eb, "diff": eb - ea,
        "diff_lo": np.quantile(diff, lo_q, axis=0), "diff_hi": np.quantile(diff, hi_q, axis=0),
        "log2_ratio": np.log2((eb + pseudo) / (ea + pseudo)),
        "log2_lo": np.quantile(lr, lo_q, axis=0), "log2_hi": np.quantile(lr, hi_q, axis=0),
    }, index=cols)
    out["ci_excludes_0"] = (out["diff_lo"] > 0) | (out["diff_hi"] < 0)
    return out


def _probabilities(pred: pd.DataFrame, model: TransferModel) -> np.ndarray:
    return pred[[f"p:{c}" for c in model.classes]].to_numpy(dtype=float)


def _estimates(pred, model, level, subset, estimators, clusters, n_boot, seed, em):
    """Categories, point estimates and cluster-bootstrap draws (None without clusters) per estimator."""
    _check_level(level)
    bad = sorted(set(estimators) - set(ESTIMATORS))
    if bad:
        raise ValueError(f"unknown estimators {bad}; choose from {ESTIMATORS}")
    P = _probabilities(pred, model)
    A = model.projection(level)
    cats = [str(c) for c in model.categories(level)]
    cols = None
    if subset is not None:
        if level != "label":
            raise ValueError("subset is a set of fine labels; use it with level='label'")
        cats = [str(c) for c in subset]
        where = {str(c): i for i, c in enumerate(model.classes)}
        missing = [c for c in cats if c not in where]
        if missing:
            raise KeyError(f"subset labels not in the model: {missing}")
        cols = np.array([where[c] for c in cats])
    res = None
    if clusters is not None:
        inv, mult = _resample(clusters, n_boot, seed)
        G = sparse.csr_matrix((np.ones(inv.size), (inv, np.arange(inv.size))), shape=(mult.shape[1], inv.size))
        res = (G, mult)

    def ratio(num, den):
        est = num.sum(0) / den.sum()
        if res is None:
            return est, None
        G, mult = res
        tot = mult @ (G @ den)
        return est, (mult @ (G @ num)) / np.where(tot > 0, tot, np.nan)[:, None]

    if subset is None:
        top = one_hot(pred[f"{level}_raw"].astype(str).to_numpy(), cats).to_numpy()
        den = np.ones(len(P))
        mean = P @ A
    else:
        inside = np.isin(pred["label_raw"].astype(str).to_numpy(), cats)
        top = one_hot(pred["label_raw"].astype(str).to_numpy(), cats).to_numpy()
        den = inside.astype(float)
        Ps = P[:, cols]
        mean = Ps / np.clip(Ps.sum(1, keepdims=True), 1e-300, None) * inside[:, None]
    point, draws = {}, {}
    if "top_call" in estimators:
        point["top_call"], draws["top_call"] = ratio(top, den)
    if "mean_prob" in estimators:
        point["mean_prob"], draws["mean_prob"] = ratio(mean, den)
    if "em" in estimators:
        if em is None:
            em = (bootstrap_em(P, model.priors, clusters, n_boot, seed) if clusters is not None
                  else (em_priors(P, model.priors)[0], None))
        q0, D = em

        def proj(v):
            if cols is None:
                return v @ A
            s = v[..., cols]
            return s / s.sum(-1, keepdims=True)

        point["em"] = proj(np.asarray(q0))
        draws["em"] = proj(np.asarray(D)) if (D is not None and res is not None) else None
    return cats, point, draws


def composition_bracket(pred: pd.DataFrame, model: TransferModel, level: str = "group", clusters=None, *,
                        subset=None, estimators=BRACKET, n_boot: int = 2000, seed: int = 0, ci: float = 0.95,
                        em=None) -> pd.DataFrame:
    """Composition of a population by top call and by EM together: the default for a query from another file.

    Top call counts each cell's most probable category (at `level`: argmax of the summed probabilities); EM
    re-estimates the fine-label priors (`em_priors`) and sums them to `level`. The two bracket the composition: top
    call is pulled toward the reference's composition when the query is out of distribution, EM can over-correct
    when cells of one type differ between conditions. `subset` (fine labels, level "label") gives the composition
    within those labels: top call over cells whose top label is in the subset, EM priors renormalised over it.
    "mean_prob" (mean probability) can be added to `estimators` as a sensitivity column.

    Rows: categories. Columns: one per estimator, `bracket_low`/`bracket_high` (min and max of the top-call and EM
    estimates) and, with `clusters` (e.g. `section_ids`), `<estimator>_lo`/`_hi` percentile intervals from resampling
    whole clusters and their envelope `envelope_lo`/`envelope_hi`. `em`: precomputed `bootstrap_em(P, model.priors,
    clusters, n_boot, seed)` (or `(priors, None)` without clusters) to avoid recomputing it.
    """
    cats, point, draws = _estimates(pred, model, level, subset, tuple(estimators), clusters, n_boot, seed, em)
    a = (1.0 - ci) / 2.0
    out = pd.DataFrame(index=pd.Index(cats, name=level))
    for e in estimators:
        out[e] = point[e]
        if draws[e] is not None:
            out[f"{e}_lo"] = np.nanquantile(draws[e], a, axis=0)
            out[f"{e}_hi"] = np.nanquantile(draws[e], 1.0 - a, axis=0)
    br = [e for e in BRACKET if e in estimators]
    if br:
        out["bracket_low"], out["bracket_high"] = out[br].min(1), out[br].max(1)
        if all(f"{e}_lo" in out for e in br):
            out["envelope_lo"] = out[[f"{e}_lo" for e in br]].min(1)
            out["envelope_hi"] = out[[f"{e}_hi" for e in br]].max(1)
    return out


def bracket_difference(pred_a: pd.DataFrame, pred_b: pd.DataFrame, model: TransferModel, level: str = "group",
                       clusters_a=None, clusters_b=None, *, subset=None, estimators=BRACKET, n_boot: int = 2000,
                       seed: int = 0, seeds: tuple[int, int] | None = None, ci: float = 0.95, em_a=None, em_b=None,
                       pseudo: float = 1e-3) -> pd.DataFrame:
    """Composition of b minus a by top call and EM (the cross-file default), with independent cluster bootstraps.

    Both populations must be predictions in the same model's vocabulary (e.g. the reference's out-of-fold calls from
    `cross_validate` and a query's from `predict`). Columns: `a_<est>`, `b_<est>`, `<est>_diff`,
    `<est>_log2_ratio` (with `pseudo`), `bracket_low`/`bracket_high` (min and max of the top-call and EM
    differences) and, with clusters for both, `<est>_diff_lo`/`_hi`, `<est>_excludes_0`, `envelope_lo`/`_hi` and
    `bracket_excludes_0` (top call and EM both exclude 0 on the same side). Population a is resampled with
    `seeds[0]` (default `seed`) and b with `seeds[1]` (default `seed + 1`); `em_a`/`em_b` are precomputed
    `bootstrap_em` results with those seeds.

    The intervals cover resampling of clusters within each population only. With sections of one embryo per
    condition they measure section sampling, not embryo-to-embryo, coverage or genotype variation.
    """
    sa, sb = (seed, seed + 1) if seeds is None else seeds
    est = tuple(estimators)
    with_ci = clusters_a is not None and clusters_b is not None
    cats, pa, da = _estimates(pred_a, model, level, subset, est, clusters_a if with_ci else None, n_boot, sa, em_a)
    _, pb, db = _estimates(pred_b, model, level, subset, est, clusters_b if with_ci else None, n_boot, sb, em_b)
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    out = pd.DataFrame(index=pd.Index(cats, name=level))
    for e in est:
        out[f"a_{e}"], out[f"b_{e}"] = pa[e], pb[e]
        out[f"{e}_diff"] = pb[e] - pa[e]
        out[f"{e}_log2_ratio"] = np.log2((pb[e] + pseudo) / (pa[e] + pseudo))
        if da[e] is not None and db[e] is not None:
            diff = db[e] - da[e]
            out[f"{e}_diff_lo"] = np.nanquantile(diff, lo_q, axis=0)
            out[f"{e}_diff_hi"] = np.nanquantile(diff, hi_q, axis=0)
            out[f"{e}_excludes_0"] = (out[f"{e}_diff_lo"] > 0) | (out[f"{e}_diff_hi"] < 0)
    br = [e for e in BRACKET if e in est]
    if br:
        d = out[[f"{e}_diff" for e in br]]
        out["bracket_low"], out["bracket_high"] = d.min(1), d.max(1)
        if all(f"{e}_diff_lo" in out for e in br):
            lo = out[[f"{e}_diff_lo" for e in br]]
            hi = out[[f"{e}_diff_hi" for e in br]]
            out["envelope_lo"], out["envelope_hi"] = lo.min(1), hi.max(1)
            out["bracket_excludes_0"] = (lo.min(1) > 0) | (hi.max(1) < 0)
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Stack coverage and z-trim sensitivity


def z_profile(pred: pd.DataFrame, zl: ZLevels, level: str = "group", n_bins: int = 10,
              categories=None) -> pd.DataFrame:
    """Top-call composition along the stack: cells binned by relative z position (`z_position`, 0 = lowest level)
    into `n_bins` equal-width bins. Rows: bins; columns bin_lo, bin_hi, n_cells, n_levels, then one fraction per
    category (`categories`, default those called)."""
    pos = z_position(zl)
    b = np.minimum((pos * n_bins).astype(int), n_bins - 1)
    calls = pred[f"{level}_raw"].astype(str).to_numpy()
    cats = list(pd.unique(calls)) if categories is None else [str(c) for c in categories]
    frac = pd.crosstab(b, calls, normalize="index").reindex(index=range(n_bins), columns=cats, fill_value=0.0)
    out = pd.DataFrame({"bin_lo": np.arange(n_bins) / n_bins, "bin_hi": (np.arange(n_bins) + 1) / n_bins,
                        "n_cells": np.bincount(b, minlength=n_bins),
                        "n_levels": [int(np.unique(np.asarray(zl.index)[b == i]).size) for i in range(n_bins)]})
    return pd.concat([out, frac.fillna(0.0).reset_index(drop=True)], axis=1)


def z_location(pred: pd.DataFrame, zl: ZLevels, level: str = "group", ends: float = 0.1) -> pd.DataFrame:
    """Where each category's cells (top calls) sit along the stack: n, share of all cells, mean and quartiles of the
    relative z position, and the share of the category's cells in the lowest and highest `ends` of the stack."""
    pos = z_position(zl)
    calls = pd.Series(pred[f"{level}_raw"].astype(str).to_numpy())
    rows = []
    for c in pd.unique(calls):
        p = pos[(calls == c).to_numpy()]
        rows.append({"category": c, "n": int(p.size), "share": float(p.size / pos.size), "mean_pos": float(p.mean()),
                     "q25": float(np.quantile(p, 0.25)), "median_pos": float(np.median(p)),
                     "q75": float(np.quantile(p, 0.75)), f"low{ends:g}_share": float((p <= ends).mean()),
                     f"high{ends:g}_share": float((p >= 1 - ends).mean())})
    return pd.DataFrame(rows).set_index("category")


def trim_sensitivity(pred_a: pd.DataFrame, zl_a: ZLevels, pred_b: pd.DataFrame, zl_b: ZLevels,
                     model: TransferModel, *, trims: dict | None = None, specs: dict | None = None,
                     estimators=BRACKET, n_boot: int = 1000, seed: int = 0, ci: float = 0.95,
                     n_jobs: int = 1) -> pd.DataFrame:
    """b - a composition differences when end sections of each stack are dropped (stack-coverage sensitivity).

    trims: {name: (a low, a high, b low, b high)} end levels dropped (default `default_trims()`, 8 trims). specs:
    {name: (level, subset)}, default {"group": ("group", None)}; e.g. ("label", cm_labels) for within-CM fractions.
    Sections (snapped z levels) are the bootstrap clusters; EM draws are computed once per trimmed population
    (`n_jobs` threads, see `bootstrap_em`).
    Long table: trim, spec, category, estimator, a, b, diff, diff_lo, diff_hi, excludes_0, cells and levels kept.
    Trims leaving either stack with fewer than two levels are skipped and listed in `attrs["skipped_trims"]`.
    """
    trims = default_trims() if trims is None else trims
    specs = {"group": ("group", None)} if specs is None else specs
    est = tuple(estimators)
    Pa, Pb = _probabilities(pred_a, model), _probabilities(pred_b, model)
    cache: dict = {}

    def em_for(tag, P, mask, clusters, s):
        key = (tag, *mask)
        if "em" in est and key not in cache:
            cache[key] = bootstrap_em(P, model.priors, clusters, n_boot, s, n_jobs=n_jobs)
        return cache.get(key)

    rows, skipped = [], []
    for tname, (la, ha, lb, hb) in trims.items():
        ma, mb = trim_mask(zl_a, la, ha), trim_mask(zl_b, lb, hb)
        ca, cb = np.asarray(zl_a.index)[ma], np.asarray(zl_b.index)[mb]
        if np.unique(ca).size < 2 or np.unique(cb).size < 2:
            skipped.append(tname)
            continue
        em_a = em_for("a", Pa[ma], (la, ha), ca, seed)
        em_b = em_for("b", Pb[mb], (lb, hb), cb, seed + 1)
        for sname, (level, subset) in specs.items():
            d = bracket_difference(pred_a[ma], pred_b[mb], model, level, ca, cb, subset=subset, estimators=est,
                                   n_boot=n_boot, seeds=(seed, seed + 1), ci=ci, em_a=em_a, em_b=em_b)
            for cat in d.index:
                for e in est:
                    rows.append({"trim": tname, "spec": sname, "category": cat, "estimator": e,
                                 "a": d.at[cat, f"a_{e}"], "b": d.at[cat, f"b_{e}"], "diff": d.at[cat, f"{e}_diff"],
                                 "diff_lo": d.at[cat, f"{e}_diff_lo"], "diff_hi": d.at[cat, f"{e}_diff_hi"],
                                 "excludes_0": bool(d.at[cat, f"{e}_excludes_0"]), "cells_a": int(ma.sum()),
                                 "cells_b": int(mb.sum()), "levels_a": int(np.unique(ca).size),
                                 "levels_b": int(np.unique(cb).size)})
    if skipped:
        warnings.warn(f"trims skipped (fewer than two levels left in a stack): {skipped}", stacklevel=2)
    out = pd.DataFrame(rows)
    out.attrs["skipped_trims"] = skipped
    return out


def trim_summary(table: pd.DataFrame, untrimmed: str = "none") -> pd.DataFrame:
    """Per spec, category and estimator from `trim_sensitivity`: the untrimmed difference, its range over trims,
    whether every trim's interval excludes 0 on the side of the untrimmed difference, and whether the sign flips."""
    rows = []
    for (spec, cat, e), g in table.groupby(["spec", "category", "estimator"], sort=False):
        base = g.loc[g["trim"] == untrimmed, "diff"]
        d0 = float(base.iloc[0]) if len(base) else float("nan")
        side = np.sign(d0) if np.isfinite(d0) else 0.0
        same_side = (g["diff_lo"] > 0) if side > 0 else (g["diff_hi"] < 0)
        rows.append({"spec": spec, "category": cat, "estimator": e, "diff_untrimmed": d0,
                     "diff_min": float(g["diff"].min()), "diff_max": float(g["diff"].max()),
                     "lo_min": float(g["diff_lo"].min()), "hi_max": float(g["diff_hi"].max()),
                     "excludes_0_all_trims": bool(side != 0 and same_side.all()),
                     "n_trims_excluding_0": int(same_side.sum()) if side != 0 else 0,
                     "sign_flips": bool((np.sign(g["diff"]) != side).any()), "n_trims": int(len(g))})
    return pd.DataFrame(rows)
