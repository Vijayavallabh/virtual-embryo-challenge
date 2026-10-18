"""Map cells from any MERFISH file into a reference file's cell-type vocabulary.

Each released file has its own labels (the Mab21l2 KO says vCM, LPM, Paraxial Mesoderm, Branchial Arch; WT E9.5 says
V-CM, A-CM, Peri, Dorsal-PAM, aPHM...), so comparing compositions across files needs relabelling on the shared panel
genes. The transfer model is

    counts on the shared genes -> log1p CP10k -> per-gene scaling -> PCA -> classifier

with the scaling and the PCA fitted on the reference only. The classifier is multinomial logistic regression
("logreg") or a k-nearest-neighbour vote ("knn"); `classifier="auto"` picks the one with the higher cross-validated
macro-F1 on the reference. It predicts the reference's fine labels. Coarse lineage groups (configs/lineage_map.yaml)
are sums of fine-label probabilities, so fine and coarse calls never contradict each other. Cells whose top
probability is below `threshold` are "Unassigned". `novelty` is the larger of two ratios, each scaled by the
reference's own 99th percentile: the mean distance to the k nearest reference cells in PCA space (unlike any reference
cell along the reference's main axes) and the distance from the reference PCA subspace (expression the reference
axes do not explain, e.g. marker genes that never vary in the reference). Above 1, a cell is more unusual than 99% of
reference cells on at least one of the two.

Depth: median counts per cell are 75 at E9.5, 95 in the Mab21l2 KO and 130 at E8.75. When the reference is deeper,
thin it to the query's depth with `fit_transfer(depth=median_depth(query))` (E8.75 -> E9.5 group accuracy 0.655 ->
0.675, experiments/harmonise). When the query is deeper, `predict` thins it binomially to the reference depth and
averages `depth_draws=5` independent draws; one draw lost about 2 points of group accuracy (E9.5 -> E8.75: 0.769 vs
0.788 unthinned), five draws matched no thinning (0.790) while removing the technical shift seen by `novelty`.

Also here: out-of-fold predictions on the reference (`cross_validate`), classification metrics, composition estimates
(hard calls, mean probabilities, or EM re-estimation of the priors under prior shift) and cluster-bootstrap intervals
over sections.

    model = fit_transfer(wt, "celltype")
    calls = predict(model, ko)          # one row per query cell
"""

from __future__ import annotations

import time
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

UNASSIGNED = "Unassigned"
METHODS = ("logreg", "knn")


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


def section_ids(a: ad.AnnData, decimals: int = 0) -> np.ndarray:
    """Section label per cell from the z coordinate of `obsm['spatial_3D']`.

    In the MERFISH files z is a section index times a per-stage factor; rounding removes floating-point noise
    (WT E9.5 stores 98 distinct z values for 40 sections).
    """
    z = np.round(np.asarray(a.obsm[COORD_KEY])[:, 2].astype(float), decimals) + 0.0
    return np.array([f"z{v:g}" for v in z])


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
    """A fitted label-transfer model. Probabilities are over `classes`; `groups` are their lineage groups."""

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

    @property
    def membership(self) -> np.ndarray:
        """(n_classes, n_groups) indicator: which group each fine label belongs to."""
        return group_membership(self.classes, self.groups, self.group_of)


def group_membership(classes, groups, group_of: dict[str, str]) -> np.ndarray:
    """(n_classes, n_groups) 0/1 matrix: column j marks the labels whose group is groups[j]."""
    own = np.array([group_of[c] for c in classes])
    return (own[:, None] == np.asarray(groups)[None, :]).astype(float)


def group_lookup(labels, group_map: dict[str, str] | None = None) -> dict[str, str]:
    """Lineage group of each label (configs/lineage_map.yaml unless `group_map` is given); unmapped -> "Unknown"."""
    gm = dict(lineage_map()["labels"]) if group_map is None else dict(group_map)
    return {str(lab): gm.get(str(lab), "Unknown") for lab in labels}


def group_order(groups) -> np.ndarray:
    """Groups in lineage-map order, then any others alphabetically."""
    present = {str(g) for g in groups}
    known = [g for g in lineage_map()["groups"] if g in present]
    return np.array(known + sorted(present - set(known)), dtype=object)


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


def _novelty(pipe, nov, X):
    """(kNN component, residual component) of the novelty score, each over its reference 99th percentile."""
    if nov is None:
        return np.full(X.shape[0], np.nan), np.full(X.shape[0], np.nan)
    nn, knn_scale, res_scale = nov
    _, Z, res = _embed(pipe, X)
    d, _ = nn.kneighbors(Z)
    return d.mean(1) / knn_scale, res / res_scale


def _frame(P, nov_parts, classes, groups, M, index, threshold, group_threshold, novelty_threshold, proba):
    nov_knn, nov_res = nov_parts
    nov = np.fmax(nov_knn, nov_res)
    rows = np.arange(len(P))
    i = P.argmax(1)
    conf = P[rows, i]
    PG = P @ M
    j = PG.argmax(1)
    gconf = PG[rows, j]
    novel = nov > novelty_threshold if novelty_threshold is not None else np.zeros(len(P), dtype=bool)
    gthr = threshold if group_threshold is None else group_threshold
    label_raw = np.asarray(classes, dtype=object)[i]
    group_raw = np.asarray(groups, dtype=object)[j]
    df = pd.DataFrame({
        "label": np.where((conf >= threshold) & ~novel, label_raw, UNASSIGNED),
        "label_raw": label_raw,
        "label_conf": conf,
        "group": np.where((gconf >= gthr) & ~novel, group_raw, UNASSIGNED),
        "group_raw": group_raw,
        "group_conf": gconf,
        "novelty": nov,
        "novelty_knn": nov_knn,
        "novelty_residual": nov_res,
    }, index=index)
    if proba:
        df = pd.concat([df, pd.DataFrame(P, index=index, columns=[f"p:{c}" for c in classes]),
                        pd.DataFrame(PG, index=index, columns=[f"pg:{g}" for g in groups])], axis=1)
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
    group_of = group_lookup(classes, group_map)
    gorder = group_order(group_of.values())
    M = group_membership(classes, gorder, group_of)
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
    df = _frame(P, (nov_knn, nov_res), classes, gorder, M, index, threshold, None, None, proba=True)
    df.insert(0, "fold", fold)
    df.insert(0, "true_group", [group_of[lab] for lab in y])
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

    Same columns as `predict` plus `true`, `true_group` and `fold`. Folds are stratified by label; with `groups`
    (one id per reference cell, e.g. `section_ids`) whole groups are held out together (StratifiedGroupKFold).
    `max_cells` evaluates on a stratified subsample. Scaling, PCA and classifier are refitted inside every fold.
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
    df = _cv_arrays(X, y, index, method=method, n_folds=n_folds, groups=grp, seed=seed, n_pcs=n_pcs, k=k, C=C,
                    clip=clip, n_jobs=n_jobs, novelty=novelty, k_novelty=k_novelty,
                    novelty_quantile=novelty_quantile, threshold=threshold, group_map=group_map)
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
    dropped too (both recorded in `model.dropped`). group_map: label -> coarse group (default: the lineage map).
    depth: thin the reference to this median counts per cell first (use the query's depth when it is shallower).
    """
    if classifier not in (*METHODS, "auto"):
        raise ValueError(f"classifier must be one of {(*METHODS, 'auto')}, got {classifier!r}")
    genes = list(reference.var_names) if genes is None else list(genes)
    X, y, keep, ref_depth, dropped = _prepare(reference, label_key, genes, exclude_labels,
                                              max(min_class_cells, cv_folds), depth, seed, depth_tol)
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
                             novelty_quantile=novelty_quantile, threshold=threshold, group_map=group_map)
            rows.append({"method": m, **transfer_metrics(oof), "seconds": time.perf_counter() - t0})
        cv = pd.DataFrame(rows)
        method = str(cv.loc[cv["macro_f1"].idxmax(), "method"])
    pipe, (nn, knn_scale, res_scale) = _fit_one(X, y, method, n_pcs, k, C, clip, seed, n_jobs, True, k_novelty,
                                                novelty_quantile)
    classes = np.asarray(pipe.classes_)
    group_of = group_lookup(classes, group_map)
    return TransferModel(
        genes=genes, label_key=label_key, method=method, classes=classes,
        groups=group_order(group_of.values()), group_of=group_of, pipeline=pipe,
        priors=np.array([(y == c).mean() for c in classes]), depth=ref_depth, neighbours=nn,
        novelty_scale=knn_scale, residual_scale=res_scale, threshold=threshold, cv=cv, dropped=dropped,
        params=params,
    )


def predict(model: TransferModel, query: ad.AnnData, *, threshold: float | None = None,
            group_threshold: float | None = None, novelty_threshold: float | None = None, depth_match: bool = True,
            depth_draws: int = 5, depth_tol: float = 0.05, seed: int = 0, proba: bool = True) -> pd.DataFrame:
    """Map query cells into the model's vocabulary; one row per query cell (index = query obs names).

    Columns: `label` (fine label, or "Unassigned" when its probability is below `threshold`, default the model's),
    `label_raw` (top label regardless), `label_conf`; `group`, `group_raw`, `group_conf` for the lineage groups
    (threshold `group_threshold`, default `threshold`); `novelty` (max of `novelty_knn` and `novelty_residual`, see
    the module docstring); and, with `proba`, probabilities `p:<label>` and `pg:<group>`. `novelty_threshold` also
    unassigns cells whose novelty exceeds it. With `depth_match` and a counts layer, a query deeper than the
    reference (median counts per cell) is thinned to the reference depth first; `df.attrs["depth_frac"]` records the
    fraction of molecules kept. Probabilities and novelty are averaged over `depth_draws` independent thinning draws
    (a single draw throws information away; see the module docstring).
    """
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
        nk, nr = _novelty(model.pipeline, novelty_ref, X)
        nov_knn, nov_res = nov_knn + nk / len(draws), nov_res + nr / len(draws)
    nov = (nov_knn, nov_res)
    thr = model.threshold if threshold is None else threshold
    df = _frame(P, nov, model.classes, model.groups, model.membership, query.obs_names, thr, group_threshold,
                novelty_threshold, proba)
    df.attrs.update(depth_frac=frac, depth_draws=len(draws), method=model.method)
    return df


# ---------------------------------------------------------------------------------------------------------------------
# Evaluation and composition


def transfer_metrics(pred: pd.DataFrame, true_col: str = "true", true_group_col: str = "true_group") -> dict:
    """Accuracy and macro-F1 of the top calls (fine and group level), coverage and accuracy of the thresholded calls,
    and log-loss when probabilities are present. Macro-F1 averages over the labels present in the truth."""
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
    return out


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


def composition(pred: pd.DataFrame, level: str = "group", method: str = "hard",
                model: TransferModel | None = None) -> pd.Series:
    """Fraction of cells per category from `predict` output.

    method: "hard" counts the thresholded calls (Unassigned included), "raw" the top calls, "soft" averages the
    probabilities, "em" re-estimates the priors by EM (needs `model`; group fractions sum fine-label priors).
    level: "group" or "label".
    """
    if level not in ("group", "label"):
        raise ValueError("level must be 'group' or 'label'")
    if method in ("hard", "raw"):
        col = level if method == "hard" else f"{level}_raw"
        return pred[col].value_counts(normalize=True)
    if method == "soft":
        prefix = "pg:" if level == "group" else "p:"
        cols = [c for c in pred.columns if c.startswith(prefix)]
        return pd.Series(pred[cols].to_numpy().mean(0), index=[c[len(prefix):] for c in cols])
    if method == "em":
        if model is None:
            raise ValueError("method='em' needs the model (training priors)")
        q, _ = em_priors(pred[[f"p:{c}" for c in model.classes]].to_numpy(), model.priors)
        if level == "group":
            return pd.Series(q @ model.membership, index=model.groups)
        return pd.Series(q, index=model.classes)
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
    ids, inv = np.unique(np.asarray(clusters), return_inverse=True)
    S = np.zeros((len(ids), W.shape[1]))
    np.add.at(S, inv, W)
    n = np.bincount(inv, minlength=len(ids)).astype(float)
    mult = np.random.default_rng(seed).multinomial(len(ids), np.full(len(ids), 1.0 / len(ids)), size=n_boot)
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
    and ci_excludes_0. Categories are aligned by name (missing ones count as 0).
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
