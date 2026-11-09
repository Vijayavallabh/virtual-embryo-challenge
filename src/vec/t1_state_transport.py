"""T1 head: within-type state transport, edited into real cells' counts per metacell.

The machinery of the O1 diagnostic arms (tries T8 to T10, `experiments/t1_scope/pseudo_task.py`) and of candidate g_mom
(try T11, `vec.harness`), fixed in their experiments/harness/TRIES.md rows before their runs. Nothing here changes
B0, M1 or any other entry: the module is new and its functions are called only under the new candidate name and the
new script option.

1. A transport within each cell category. `fit_state_factors` couples two stages' cells of one category (POT's entropic
   Sinkhorn on squared Euclidean distances in a joint PCA; epsilon a fixed fraction of the median cost, in log space for
   numerical stability), groups one stage's cells into k-means metacells (K = clip(floor(n / 150), 5, 40), random_state
   0) and takes, per metacell, the log ratio of mean CP10k between the later and the earlier stage's cells, the
   metacell's own mean against its transported partners' coupling-weighted mean (pseudocount 0.5, the shift's). Which
   stage carries the metacells is the caller's: `metacells_on="earlier"` for O1 (the edited cells are the E8.5 cells
   that move towards the sibling library), `"later"` for g_mom (the edited cells are the reference cells, which keep
   moving as they did since the stage before). A category with too few cells at a stage, or drawn from stand-in labels
   there, gets no factors, as in `vec.t1_scope.momentum_shift` (min_cells 20, standin off).
2. Two finishes, each fixed in its row. `finish_oracle` (O1): the factors as they are, zeroed below |log f| < log 1.05
   and clipped to +-log 3. `finish_momentum` (g_mom): the per-category cell-weighted mean of the log factor removed
   per gene (a within-type edit, on top of the per-type momentum of B0), damped, then zeroed and clipped the same way.
3. `apply_state_factors` edits a re-weighted population's counts by the factor of each cell's metacell (`vec.core.
   apply_factors`: binomial thinning and Poisson top-up of its own counts, zeros kept, held-flat genes and the cell's
   expected depth kept). With `rescale=True` (O1), each category's per-gene factors are first multiplied by the constant
   that keeps its mean CP10k unchanged, the per-cell keep-depth constant of apply_factors included, so the per-type
   pseudobulk is held and the edit changes only the within-type shape (the row's required order: zero and clip, then
   the rescale).

`r1_unbracketed_labels` names the labels of R1 (experiments/harness/TRIES.md): cells whose label maps to no group with a
Feng direction, which g_mom leaves alone at test and on its backtests alike (its row).
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path

import anndata as ad
import numpy as np
import yaml
from scipy import sparse

from vec.boards import REPO_ROOT
from vec.core import apply_factors, counts_of, log1p_cp10k, mean_cp10k
from vec.t1_scope import CATEGORY_KEY, LABEL_KEY, ChamberSpec, _type_rows, classify, load_spec, momentum_shift, reweight

N_HVG = 2000
N_PCS = 30
HVG_MIN_DETECTION = 0.01  # the rule of experiments/t1_scope/label_evidence.py: dispersion among genes in >= 1% of cells
OT_EPSILON_FRACTION = 0.05  # of the median cost
OT_MAX_ITER = 5000
OT_TOL = 1e-5  # the marginals' L1 error at which Sinkhorn stops (recorded per category, with its value at the stop)
METACELL_DIVISOR = 150
METACELL_MIN = 5
METACELL_MAX = 40
KMEANS_SEED = 0
PSEUDOCOUNT = 0.5
MIN_CELLS = 20  # as vec.t1_scope.momentum_shift
ZERO_BELOW = float(np.log(1.05))
CLIP_AT = float(np.log(3.0))
RESCALE_ITERS = 100
RESCALE_TOL = 1e-6
R1_FENG_GROUPS = ("CM", "Endocardium", "Endothelium", "Epicardium")  # the groups with a Feng direction (build_t1, N12)
R1_CONFIG = REPO_ROOT / "experiments" / "t1_direction" / "config.yaml"

Seed = int | Sequence[int]


@dataclass
class CategoryFactors:
    """One category's per-metacell log factors: `logf[m, g]` for the metacell cells `index` names (cell name -> m)."""

    category: str
    index: dict[str, int]
    sizes: np.ndarray
    logf: np.ndarray
    record: dict = field(default_factory=dict)

    @property
    def n_metacells(self) -> int:
        return len(self.sizes)


@dataclass
class StateFactors:
    """Per-category metacell factors on one gene order, and what `fit_state_factors` skipped (category -> why)."""

    genes: list[str]
    metacells_on: str
    categories: dict[str, CategoryFactors]
    skipped: dict[str, str]
    record: dict = field(default_factory=dict)


def hvg_columns(X, n_top: int = N_HVG, min_detection: float = HVG_MIN_DETECTION) -> np.ndarray:
    """Columns of the n_top highest dispersion genes (variance / mean of the entries of X, log1p CP10k) among genes
    detected in at least min_detection of the rows (experiments/t1_scope/label_evidence.py's rule, on these rows);
    ties in column order. Fewer columns when fewer genes pass."""
    X = sparse.csr_matrix(X, dtype=np.float64)
    mean = np.asarray(X.mean(0)).ravel()
    sq = np.asarray(X.multiply(X).mean(0)).ravel()
    det = np.asarray((X > 0).mean(0)).ravel()
    ok = (det >= min_detection) & (mean > 0)
    disp = np.where(ok, (sq - mean**2) / np.maximum(mean, 1e-12), -np.inf)
    top = np.argsort(-disp, kind="stable")[:n_top]
    return np.sort(top[ok[top]])


def joint_pcs(a: ad.AnnData, b: ad.AnnData, n_hvg: int = N_HVG, n_pcs: int = N_PCS,
              seed: int = 0) -> tuple[np.ndarray, np.ndarray, dict]:
    """The two stages' cells in one PCA: scikit-learn's randomized PCA of the pooled, centred log1p CP10k on the top
    n_hvg dispersion genes of the pooled cells (rows of a first). Returns (a's PCs, b's PCs, a record)."""
    from sklearn.decomposition import PCA

    X = sparse.vstack([sparse.csr_matrix(a.X), sparse.csr_matrix(b.X)], format="csr")
    cols = hvg_columns(X, n_hvg)
    k = min(n_pcs, min(X.shape[0], len(cols)) - 1)
    if k < 1:
        raise ValueError("too few cells or variable genes for a PCA")
    pcs = PCA(n_components=k, svd_solver="randomized", random_state=seed).fit_transform(
        X[:, cols].toarray().astype(np.float32))
    pcs = np.asarray(pcs, dtype=np.float64)
    return pcs[:a.n_obs], pcs[a.n_obs:], {"n_hvg": int(len(cols)), "n_pcs": int(k), "pca_seed": int(seed)}


def sinkhorn_plan(xa, xb, epsilon_fraction: float = OT_EPSILON_FRACTION, max_iter: int = OT_MAX_ITER,
                  tol: float = OT_TOL) -> tuple[np.ndarray, dict]:
    """Balanced entropic transport between two point clouds with uniform masses: squared Euclidean cost divided by its
    median, Sinkhorn in log space with epsilon = epsilon_fraction (so epsilon_fraction x the median cost in the
    original units). Returns (plan with rows a and columns b, a record with the last marginal error)."""
    import ot

    xa, xb = np.asarray(xa, np.float64), np.asarray(xb, np.float64)
    M = ot.dist(xa, xb, metric="sqeuclidean")
    median_cost = float(np.median(M))
    if not median_cost > 0:
        raise ValueError("the cost matrix has no spread")
    a, b = np.full(len(xa), 1.0 / len(xa)), np.full(len(xb), 1.0 / len(xb))
    T, log = ot.sinkhorn(a, b, M / median_cost, epsilon_fraction, method="sinkhorn_log", numItermax=max_iter,
                         stopThr=tol, log=True, warn=False)
    errs = [float(e) for e in log.get("err", [])]
    return np.asarray(T, dtype=np.float64), {"median_cost": median_cost, "last_marginal_error": errs[-1] if errs else
                                              float("nan"), "converged": bool(errs and errs[-1] <= tol)}


def metacell_labels(x, seed: int = KMEANS_SEED) -> np.ndarray:
    """k-means metacells of the rows of x (K = clip(floor(n / 150), 5, 40), at most n; random_state fixed, ten
    restarts), labels 0 ... K - 1, relabelled to run without gaps."""
    from sklearn.cluster import KMeans

    n = np.asarray(x).shape[0]
    K = min(n, int(np.clip(n // METACELL_DIVISOR, METACELL_MIN, METACELL_MAX)))
    if K <= 1:
        return np.zeros(n, dtype=np.int64)
    labels = KMeans(n_clusters=K, n_init=10, random_state=seed).fit_predict(np.asarray(x, np.float64))
    return np.unique(labels, return_inverse=True)[1].astype(np.int64)


def fit_state_factors(earlier: ad.AnnData, later: ad.AnnData, metacells_on: str, categories: Iterable[str], *,
                      key: str = CATEGORY_KEY, label_key: str = LABEL_KEY,
                      within: Mapping[str, Sequence[str]] | None = None,
                      earlier_fallback: Mapping[str, Sequence[str]] | None = None,
                      later_fallback: Mapping[str, Sequence[str]] | None = None,
                      min_cells: int = MIN_CELLS, standin: bool = False, pseudocount: float = PSEUDOCOUNT,
                      n_hvg: int = N_HVG, n_pcs: int = N_PCS) -> StateFactors:
    """Per category, the per-metacell log ratio later / earlier of mean CP10k, the metacell's own mean against its
    transported partners' coupling-weighted mean (module docstring, step 1). earlier, later: classified stages on one
    gene order. metacells_on: "earlier" or "later" (the stage whose cells the edit will move; its row names are the
    index later cells are matched on). A category's cells at a stage are those of `vec.t1_scope._type_rows` (`within`
    restricts them to labels; stand-ins only when it has fewer than min_cells own cells); a category short of min_cells
    at a stage or drawn from stand-ins there gets no factors (with `standin` it does). Returns the raw log factors
    (`finish_oracle` or `finish_momentum` make them usable)."""
    if metacells_on not in ("earlier", "later"):
        raise ValueError("metacells_on must be 'earlier' or 'later'")
    genes = list(later.var_names)
    if list(earlier.var_names) != genes:
        raise ValueError("the two stages hold different genes")
    Ce, Cl = counts_of(earlier, genes), counts_of(later, genes)
    pe, pl, record = joint_pcs(earlier, later, n_hvg, n_pcs)
    fits: dict[str, CategoryFactors] = {}
    skipped: dict[str, str] = {}
    for cat in categories:
        rows_e, how_e = _type_rows(earlier, cat, key, label_key, within, earlier_fallback, min_cells)
        rows_l, how_l = _type_rows(later, cat, key, label_key, within, later_fallback, min_cells)
        if len(rows_e) < min_cells or len(rows_l) < min_cells:
            skipped[cat] = f"fewer than {min_cells} cells at a stage ({len(rows_e)} earlier, {len(rows_l)} later)"
            continue
        if not standin and (how_e.startswith("stand-in") or how_l.startswith("stand-in")):
            skipped[cat] = "stand-in cells at a stage (identity, not development)"
            continue
        on_earlier = metacells_on == "earlier"
        rows_a, rows_b = (rows_e, rows_l) if on_earlier else (rows_l, rows_e)
        ca, cb = (Ce, Cl) if on_earlier else (Cl, Ce)
        xa, xb = (pe, pl) if on_earlier else (pl, pe)
        names = np.asarray((earlier if on_earlier else later).obs_names)[rows_a].astype(str)
        if len(set(names)) != len(names):
            raise ValueError(f"{cat}: cell names are not unique, so drawn cells cannot be matched to their metacells")
        plan, sink = sinkhorn_plan(xa[rows_a], xb[rows_b])
        labels = metacell_labels(xa[rows_a])
        logf = _metacell_log_ratio(ca[rows_a], cb[rows_b], plan, labels, on_earlier, pseudocount)
        sizes = np.bincount(labels).astype(np.int64)
        fits[cat] = CategoryFactors(category=cat, index=dict(zip(names.tolist(), labels.tolist(), strict=True)),
                                    sizes=sizes, logf=logf, record={
                                        "category": cat, "n_cells_a": int(len(rows_a)), "n_cells_b": int(len(rows_b)),
                                        "cells_a": how_e if on_earlier else how_l,
                                        "cells_b": how_l if on_earlier else how_e, "n_metacells": int(len(sizes)),
                                        **sink, "median_abs_log_f": float(np.median(np.abs(logf)))})
    record.update({"metacells_on": metacells_on, "categories": sorted(fits), "skipped": dict(skipped)})
    return StateFactors(genes=genes, metacells_on=metacells_on, categories=fits, skipped=skipped, record=record)


def _metacell_log_ratio(ca: sparse.csr_matrix, cb: sparse.csr_matrix, plan: np.ndarray, labels: np.ndarray,
                        on_earlier: bool, pseudocount: float) -> np.ndarray:
    """Per metacell (labels on the rows of ca) and gene: log((later mean CP10k + pc) / (earlier mean + pc)), the
    metacell's own mean (over its cells) against the mean of cb weighted by the coupling mass of its cells."""
    n_meta = int(labels.max()) + 1
    lib_b = np.asarray(cb.sum(1), dtype=np.float64).ravel()
    scale_b = np.divide(1e4, lib_b, out=np.zeros_like(lib_b), where=lib_b > 0)
    cp_b = sparse.diags(scale_b) @ sparse.csr_matrix(cb, dtype=np.float64)
    out = np.empty((n_meta, ca.shape[1]), dtype=np.float64)
    for m in range(n_meta):
        rows = labels == m
        own = mean_cp10k(ca[rows])
        mass = plan[rows].sum(0)
        partner = np.asarray(cp_b.T @ (mass / mass.sum())).ravel()
        later_m, earlier_m = (partner, own) if on_earlier else (own, partner)
        out[m] = np.log(later_m + pseudocount) - np.log(earlier_m + pseudocount)
    return out


def _finished(sf: StateFactors, logf: Mapping[str, np.ndarray], note: str) -> StateFactors:
    cats = {c: replace(cf, logf=logf[c], record={**cf.record, "median_abs_log_f_finished":
                                                 float(np.median(np.abs(logf[c]))), "finish": note})
            for c, cf in sf.categories.items()}
    return replace(sf, categories=cats, record={**sf.record, "finish": note})


def _zero_and_clip(logf: np.ndarray, zero_below: float = ZERO_BELOW, clip_at: float = CLIP_AT) -> np.ndarray:
    out = np.clip(logf, -clip_at, clip_at)
    out[np.abs(out) < zero_below] = 0.0
    return out


def finish_oracle(sf: StateFactors) -> StateFactors:
    """O1's factors: |log f| below log 1.05 set to 0 and the rest clipped to +-log 3 (no damping)."""
    return _finished(sf, {c: _zero_and_clip(cf.logf) for c, cf in sf.categories.items()},
                     "zero below log 1.05, clip at log 3")


def finish_momentum(sf: StateFactors, damping: float) -> StateFactors:
    """g_mom's factors: per gene the cell-weighted mean over the category's metacells removed, damped, then zeroed
    (|log f| below log 1.05) and clipped (+-log 3)."""
    logf = {}
    for c, cf in sf.categories.items():
        mean = np.average(cf.logf, axis=0, weights=cf.sizes)
        logf[c] = _zero_and_clip((cf.logf - mean) * float(damping))
    return _finished(sf, logf, f"per-category mean removed, damped by {float(damping):g}, zero below log 1.05, clip "
                               "at log 3")


def r1_unbracketed_labels(labels: Iterable[str], config=R1_CONFIG) -> set[str]:
    """The labels among `labels` that map to no group with a Feng direction (`experiments/t1_direction/config.yaml`'s
    harmonisation regexes, ventricular CM labels counted as CM, as `experiments/t1_m1/feng_dirs.py` `label_groups`): the
    cells R1 keeps at their own expression, which g_mom leaves alone."""
    cfg = yaml.safe_load(Path(config).read_text())
    rules, include = cfg["harmonise"], set(cfg.get("cm_includes", ()))
    out = set()
    for lab in sorted(set(map(str, labels))):
        group = next((g for g, rx in rules if re.search(rx, lab)), f"other:{lab}")
        group = "CM" if group in include else group
        if group not in R1_FENG_GROUPS:
            out.add(lab)
    return out


def b0_population(ref: ad.AnnData, start: ad.AnnData, shares: Mapping[str, float], n_cells: int, seed: int,
                  held_flat: Sequence[str], *, ref_source: str, start_source: str, damping: float,
                  label_key: str = LABEL_KEY, spec: ChamberSpec | None = None) -> ad.AnnData:
    """B0's cells with their counts: `vec.t1_scope.build_submission`'s steps without the panel (the same `reweight` and
    `momentum_shift` calls, seeds and defaults, the same stand-ins), so a candidate that edits them on top of B0 starts
    from the very cells B0's own build draws (tests/test_t1_state_transport.py pins the two against each other). ref,
    start: as `_b0_t1` hands them over (`vec.harness._shallow`); both are classified in place, as there."""
    spec = spec or load_spec()
    classify(ref, label_key, ref_source, spec)
    fallback = spec.fallbacks(ref_source)
    pop, _ = reweight(ref, shares, n_cells, seed, label_key=label_key, within=None, fallback=fallback, spec=spec)
    classify(start, label_key, start_source, spec)
    shifted, _ = momentum_shift(start, ref, pop, damping, list(held_flat), label_key=label_key, within=None,
                                start_fallback=spec.fallbacks(start_source), ref_fallback=fallback, seed=[seed, 1])
    return shifted


def _factor_rows(pop: ad.AnnData, sf: StateFactors, key: str, source_key: str) -> tuple[np.ndarray, np.ndarray, dict]:
    """Per-cell row of the factor matrix (the metacells of every category in sf.categories, in that order, then one row
    of ones for cells without a factor: the row index), the matrix (rows x genes, log factors) and the blocks."""
    cats = pop.obs[key].astype(str).to_numpy()
    src = pop.obs[source_key].astype(str).to_numpy()
    n_rows = sum(cf.n_metacells for cf in sf.categories.values()) + 1
    logf = np.zeros((n_rows, len(sf.genes)), dtype=np.float64)
    row = np.full(pop.n_obs, n_rows - 1, dtype=np.int64)
    blocks = {}
    offset = 0
    for cat, cf in sf.categories.items():
        logf[offset:offset + cf.n_metacells] = cf.logf
        idx = np.flatnonzero(cats == cat)
        m = np.fromiter((cf.index.get(s, -1) for s in src[idx]), dtype=np.int64, count=len(idx))
        hit = m >= 0
        row[idx[hit]] = offset + m[hit]
        blocks[cat] = (offset, cf.n_metacells, idx[hit])
        offset += cf.n_metacells
    return row, logf, blocks


def rescale_to_hold_pseudobulk(C, factors: np.ndarray, row: np.ndarray,
                               blocks: Mapping[str, tuple[int, int, np.ndarray]], held: np.ndarray,
                               iters: int = RESCALE_ITERS, tol: float = RESCALE_TOL) -> tuple[np.ndarray, dict]:
    """Per category and gene, the multiplicative constant on every metacell's factor that keeps the category's mean
    CP10k in expectation under `apply_factors` (the sum of c * f * s / depth over the category's cells), with each
    cell's keep-depth constant s included and recomputed as the factors change (a fixed-point iteration, at most
    `iters`, stopping when the largest |log| of a step's constants is below tol). Held genes stay at 1. Returns
    (factors, per-category record); factors is modified in place."""
    C = sparse.csr_matrix(C, dtype=np.int64)
    n = C.shape[0]
    depth = np.asarray(C.sum(1), dtype=np.float64).ravel()
    depth = np.where(depth > 0, depth, 1.0)
    weights = sparse.diags(1.0 / depth) @ sparse.csr_matrix(C, dtype=np.float64)
    row_of = np.repeat(np.arange(n), np.diff(C.indptr))
    free_entry = ~held[C.indices]
    counts = C.data.astype(np.float64)
    numer = np.bincount(row_of, weights=counts * free_entry, minlength=n)
    per = {cat: weights[rows].tocsr() for cat, (_, _, rows) in blocks.items() if len(rows)}
    targets = {cat: np.asarray(w.sum(0)).ravel() for cat, w in per.items()}
    rec = {cat: {"iterations": 0, "max_abs_log_rho": 0.0} for cat in per}
    active = list(per)
    for it in range(1, iters + 1):
        f_entry = factors[row[row_of], C.indices]
        den = np.bincount(row_of, weights=counts * f_entry * free_entry, minlength=n)
        scale = np.divide(numer, den, out=np.ones(n), where=den > 0)
        step = 0.0
        for cat in list(active):
            offset, k, rows = blocks[cat]
            w = per[cat]
            cell_scale = scale[rows]
            metacell = row[rows]
            entry_row = np.repeat(np.arange(len(rows)), np.diff(w.indptr))
            value = w.data * factors[metacell[entry_row], w.indices] * cell_scale[entry_row]
            denominator = np.bincount(w.indices, weights=value, minlength=w.shape[1])
            ok = (denominator > 0) & ~held & (targets[cat] > 0)
            rho = np.ones(w.shape[1])
            rho[ok] = targets[cat][ok] / denominator[ok]
            factors[offset:offset + k] *= rho
            this = float(np.abs(np.log(rho[ok])).max()) if ok.any() else 0.0
            rec[cat] = {"iterations": it, "max_abs_log_rho": this}
            step = max(step, this)
            if this < tol:
                active.remove(cat)
        if step < tol:
            break
    return factors, rec


def apply_state_factors(pop: ad.AnnData, sf: StateFactors, held_flat: Iterable[str], seed: Seed, *, rescale: bool,
                        no_edit: np.ndarray | None = None, key: str = CATEGORY_KEY,
                        source_key: str = "source_obs") -> tuple[ad.AnnData, dict]:
    """Edit a re-weighted population (`vec.t1_scope.reweight`'s obs[key] and obs[source_key] on real cells, with counts
    in layers["counts"]) by the factor of each cell's metacell, through `apply_factors` (keep_depth). A cell without a
    metacell (its category has no factors, or its source cell is not among them) and a cell in `no_edit` (a boolean
    mask, R1's held types for g_mom) keep their counts. rescale: hold each category's mean CP10k (module docstring,
    step 3). Returns a copy with the edited counts and X = log1p CP10k, and a record."""
    if list(pop.var_names) != sf.genes:
        raise ValueError("the population's genes differ from the factors'")
    row, logf, blocks = _factor_rows(pop, sf, key, source_key)
    if no_edit is not None:
        mask = np.asarray(no_edit, dtype=bool)
        row[mask] = logf.shape[0] - 1
        blocks = {c: (o, k, r[~mask[r]]) for c, (o, k, r) in blocks.items()}
    held_set = set(map(str, held_flat))
    held = np.array([g in held_set for g in sf.genes], dtype=bool)
    C = counts_of(pop)
    factors = np.exp(logf)
    factors[:, held] = 1.0
    record = {"n_cells": int(pop.n_obs), "n_edited": int((row != logf.shape[0] - 1).sum()),
              "n_held_genes": int(held.sum()),
              "by_category": {cat: {"n_edited": int(len(rows))} for cat, (_, _, rows) in blocks.items()}}
    if rescale:
        factors, res = rescale_to_hold_pseudobulk(C, factors, row, blocks, held)
        record["rescale"] = res
    new, scale, dead = apply_factors(C, factors, row, held, keep_depth=True, seed=seed)
    out = pop.copy()
    out.layers["counts"] = new.astype(np.int32)
    out.X = log1p_cp10k(new)
    record["median_depth_scale"] = float(np.median(scale))
    record["restored_cells"] = int(dead.sum())
    return out, record


def category_mean_cp10k(C, cats: np.ndarray, category: str) -> np.ndarray:
    """A category's mean CP10k per gene over the given counts (for checks and tests)."""
    return mean_cp10k(sparse.csr_matrix(C)[np.asarray(cats) == category])


__all__ = ["CLIP_AT", "METACELL_DIVISOR", "METACELL_MAX", "METACELL_MIN", "OT_EPSILON_FRACTION", "PSEUDOCOUNT",
           "R1_FENG_GROUPS", "ZERO_BELOW", "CategoryFactors", "StateFactors", "apply_state_factors", "b0_population",
           "category_mean_cp10k", "finish_momentum", "finish_oracle", "fit_state_factors", "hvg_columns", "joint_pcs",
           "metacell_labels", "r1_unbracketed_labels", "rescale_to_hold_pseudobulk", "sinkhorn_plan"]
