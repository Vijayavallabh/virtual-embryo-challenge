"""Released-data bets on T1's last file (try T61 of experiments/harness/TRIES.md; every step here is off by default).

T1's last file as built (row R8) is the first upload's design with the outside atlas's step taken out: library 3's cells
at their own chamber mix, every cell at its E10.5 expression, then g_mom's within-type factors on cardiomyocytes and
endocardium. Try T61 reads released-data variants of it on the practice tasks. Two of its arm families change which
cells are drawn, and both touch the post-exposure block (docs/DATA_POLICY.md 2026-11-04 (d): T1's composition, except
by its legal route, a rule from released and legal data committed before any run). They are written here so that the
practice read and the test builder run the same code:

  trend_shares     share_trend: each chamber category's share moved along the released stages' own like-for-like
                   trend, own x (later / earlier) ** lam, renormalised. Only categories with a share in `own` take part
                   (no births); a category with fewer than `min_cells` own cells at either end of the trend keeps its
                   own share's weight (ratio 1). lam = 0 gives `own` back exactly.
  momentum_axes    each category's released momentum axis: damping-1 log ratio of mean CP10k between two stages of
                   the same kind, pseudocount 0.5, the held genes at 0, skipped (None) where momentum_shift skips the
                   type (fewer than min_cells cells, or stand-in cells, at either stage).
  advance_scores   z per cell: the cell's log1p CP10k projected on its category's unit axis, centred and scaled to unit
                   standard deviation within the category's cells of the pool (0 where the category has no axis).
  advance_reweight `vec.t1_scope.reweight` with each own pool cut into ADVANCE_BINS z-quantile bins and the category's
                   quota split over them in proportion to each bin's sum of exp(tau z), by `vec.core.quotas` at the
                   count; each bin drawn uniformly. Shares are set by quotas, never by importance weights
                   (docs/DATA_POLICY.md 2026-11-04 (f)(2)). It changes the shares of states inside a category
                   (composition under the block, (f)(1)) and never a category's share. tau = 0 is `reweight` itself,
                   bit for bit (the same call).
  median_library   the drawn population's median library size, (f)(3)'s report.
  assemble         `vec.t1_scope.build_submission`'s last lines (the post_cells hook, the board's gene order), for a
                   population drawn outside build_submission; tests pin it to build_submission on the no-shift path.

Nothing here reads a file, a target or a score. Every existing output is unchanged: no existing function calls this
module.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from vec.boards import Board
from vec.core import counts_of, draw_least_used, mean_cp10k, quotas
from vec.io import make_submission
from vec.t1_scope import (
    CATEGORY_KEY,
    LABEL_KEY,
    ChamberSpec,
    ReweightReport,
    _seed_record,
    _take,
    _type_rows,
    classify,
    load_shares,
    load_spec,
    reweight,
)

TREND_MIN_CELLS = 20  # as momentum_shift's min_cells: fewer own cells at either end of the trend, no trend
AXIS_MIN_CELLS = 20
AXIS_EPS = 0.5  # momentum_shift's pseudocount
LAMBDAS = (0.5, 1.0)  # share_trend's two strengths, fixed in row T61 before any run
TAUS = (0.5, 1.0)  # advance's two strengths, likewise
ADVANCE_BINS = 5  # advance's state bins per category: quintiles of z in the category's own pool, likewise


# ---------------------------------------------------------------------------------------------------------------------
# share_trend


def chamber_counts(stage: ad.AnnData, source: str, spec: ChamberSpec | None = None) -> dict[str, int]:
    """Own cells per chamber category (configs/t1_chamber_classes.yaml, the file's own overrides), from a classified
    shallow copy: `stage` is not modified."""
    spec = spec or load_spec()
    a = ad.AnnData(obs=stage.obs[[LABEL_KEY]].copy())
    cats = classify(a, source=source, spec=spec, key_added=None)
    chamber = spec.in_scope("chamber")
    return {c: int((cats == c).sum()) for c in chamber}


def trend_shares(own: Mapping[str, float], later: Mapping[str, int], earlier: Mapping[str, int], lam: float,
                 min_cells: int = TREND_MIN_CELLS) -> tuple[dict[str, float], dict]:
    """share_trend's shares: own[c] x (later share of c / earlier share of c) ** lam over the categories with a positive
    own share, renormalised. later, earlier: own-cell counts per chamber category at the two ends of the trend (the
    released stages read like for like), each turned into shares over its chamber categories. A category with fewer
    than min_cells own cells at either end gets ratio 1. Returns the shares and a record (shares and ratios, no
    count)."""
    if not (np.isfinite(lam) and lam >= 0):
        raise ValueError(f"lam must be finite and >= 0, got {lam}")
    keep = {c: float(s) for c, s in own.items() if float(s) > 0}
    if not keep:
        raise ValueError("own shares are all zero")
    tl, te = sum(later.values()), sum(earlier.values())
    if tl <= 0 or te <= 0:
        raise ValueError("a trend end has no chamber cell")
    if lam == 0:
        return dict(own), {"lam": 0.0, "ratios": {c: 1.0 for c in keep}, "trended": []}
    ratios, trended = {}, []
    for c in keep:
        nl, ne = int(later.get(c, 0)), int(earlier.get(c, 0))
        if nl >= min_cells and ne >= min_cells:
            ratios[c] = (nl / tl) / (ne / te)
            trended.append(c)
        else:
            ratios[c] = 1.0
    raw = {c: keep[c] * ratios[c] ** lam for c in keep}
    total = sum(raw.values())
    out = {c: (raw[c] / total if c in raw else 0.0) for c in own}
    return out, {"lam": float(lam), "ratios": ratios, "trended": sorted(trended), "min_cells": int(min_cells)}


# ---------------------------------------------------------------------------------------------------------------------
# advance


def momentum_axes(start: ad.AnnData, later: ad.AnnData, categories: Sequence[str], held_flat: Sequence[str], *,
                  key: str = CATEGORY_KEY, label_key: str = LABEL_KEY,
                  within: Mapping[str, Sequence[str]] | None = None,
                  start_fallback: Mapping[str, Sequence[str]] | None = None,
                  later_fallback: Mapping[str, Sequence[str]] | None = None,
                  min_cells: int = AXIS_MIN_CELLS, eps: float = AXIS_EPS) -> tuple[dict[str, np.ndarray | None], dict]:
    """Per category, momentum_shift's factor at damping 1 from `start` to `later` (log((m_later + eps) / (m_start +
    eps)), m the mean CP10k over the type's rows, `t1_scope._type_rows` as momentum_shift picks them), on `later`'s
    genes, the held genes at 0; None where momentum_shift would skip the type (fewer than min_cells rows or stand-in
    rows at either stage). Both stages classified (obs[key])."""
    genes = list(later.var_names)
    S, L = counts_of(start, genes), counts_of(later, genes)
    held = np.isin(genes, list(held_flat))
    axes, rec = {}, {}
    for t in categories:
        s_rows, s_how = _type_rows(start, t, key, label_key, within, start_fallback, min_cells)
        l_rows, l_how = _type_rows(later, t, key, label_key, within, later_fallback, min_cells)
        if len(s_rows) < min_cells or len(l_rows) < min_cells:
            axes[t], rec[t] = None, "skipped: fewer than min_cells cells at a stage"
            continue
        if s_how.startswith("stand-in") or l_how.startswith("stand-in"):
            axes[t], rec[t] = None, "skipped: stand-in cells at a stage"
            continue
        logf = np.log(mean_cp10k(L[l_rows]) + eps) - np.log(mean_cp10k(S[s_rows]) + eps)
        logf = np.asarray(logf, dtype=np.float64).ravel()
        logf[held] = 0.0
        axes[t], rec[t] = logf, "axis"
    return axes, rec


def advance_scores(pool: ad.AnnData, axes: Mapping[str, np.ndarray | None], key: str = CATEGORY_KEY) -> np.ndarray:
    """z per row of `pool` (classified; its X is log1p CP10k on the axes' genes): the projection of the cell's
    expression on its category's unit axis, centred and divided by its standard deviation (ddof 0) over the category's
    rows of the pool; 0 for a category without an axis or with no spread."""
    X = sparse.csr_matrix(pool.X, dtype=np.float64)
    cats = pool.obs[key].astype(str).to_numpy()
    z = np.zeros(pool.n_obs, dtype=np.float64)
    for c, axis in axes.items():
        rows = np.flatnonzero(cats == c)
        if axis is None or not len(rows):
            continue
        axis = np.asarray(axis, dtype=np.float64)
        if axis.shape != (X.shape[1],):
            raise ValueError(f"axis of {c!r} has shape {axis.shape}, expected ({X.shape[1]},)")
        norm = float(np.linalg.norm(axis))
        if norm == 0:
            continue
        p = np.asarray(X[rows] @ (axis / norm)).ravel()
        sd = float(p.std())
        if sd > 0:
            z[rows] = (p - p.mean()) / sd
    return z


def z_bins(zp: np.ndarray, n_bins: int = ADVANCE_BINS) -> np.ndarray:
    """Each pool row's bin, 0 to n_bins - 1, by the quantiles of its z within the pool (np.quantile's default method;
    a row on an edge goes to the upper bin). A pool whose z has no spread is one bin (0)."""
    zp = np.asarray(zp, dtype=np.float64)
    if not len(zp) or float(zp.max() - zp.min()) == 0.0:
        return np.zeros(len(zp), dtype=np.int64)
    edges = np.quantile(zp, np.arange(1, n_bins) / n_bins)
    return np.searchsorted(edges, zp, side="right").astype(np.int64)


def bin_quotas(zp: np.ndarray, bins: np.ndarray, tau: float, k: int) -> dict[int, int]:
    """The cells each bin gets of a category's quota k: shares proportional to the sum of exp(tau z) over the bin's
    pool rows (what a draw with those weights would give in expectation), turned into whole cells by `vec.core.quotas`
    at k. Shares set through quotas at the count, never through importance weights (docs/DATA_POLICY.md 2026-11-04
    (f)(2)); tau = 0 gives each bin its pool share."""
    w = np.exp(tau * np.asarray(zp, dtype=np.float64))
    mass = {int(b): float(w[bins == b].sum()) for b in np.unique(bins)}
    total = sum(mass.values())
    q = quotas({str(b): v / total for b, v in mass.items()}, k)
    return {int(b): n for b, n in q.items()}


def advance_reweight(adata: ad.AnnData, target_shares, n_cells: int, seed, z: np.ndarray, tau: float, *,
                     key: str = CATEGORY_KEY, label_key: str = LABEL_KEY,
                     within: Mapping[str, Sequence[str]] | None = None,
                     fallback: Mapping[str, Sequence[str]] | None = None,
                     spec: ChamberSpec | None = None) -> tuple[ad.AnnData, ReweightReport]:
    """`reweight` with each category's own pool split into ADVANCE_BINS bins by its z quantiles (`z_bins`) and its
    quota split over the bins by `bin_quotas` (the sum of exp(tau z) per bin, through `vec.core.quotas` at the count);
    each bin is drawn uniformly by `vec.core.draw_least_used`, bins in ascending order, one generator. A pool of
    stand-in labels is drawn as reweight draws it (its cells belong to another category's states). The category
    quotas, the plan's order, the stand-ins, the final shuffle and the obs columns are reweight's; tau = 0 calls
    reweight itself. z: one value per row of adata (`advance_scores`)."""
    if not np.isfinite(tau):
        raise ValueError(f"tau must be finite, got {tau}")
    if tau == 0:
        return reweight(adata, target_shares, n_cells, seed, key=key, label_key=label_key, within=within,
                        fallback=fallback, spec=spec)
    spec = spec or load_spec()
    shares = load_shares(target_shares, spec)
    z = np.asarray(z, dtype=np.float64)
    if z.shape != (adata.n_obs,) or not np.isfinite(z).all():
        raise ValueError("z must be finite, one value per row")
    if key not in adata.obs:
        raise KeyError(f"obs has no {key!r} column: run classify() first")
    cats = adata.obs[key].astype(str).to_numpy()
    labels = adata.obs[label_key].astype(str).to_numpy() if label_key in adata.obs else None
    for name, m in (("within", within), ("fallback", fallback)):
        unknown = sorted(set(m or {}) - set(spec.categories))
        if unknown:
            raise ValueError(f"{name} names undefined categories {unknown}")
    q = quotas(shares, n_cells)
    rng = np.random.default_rng(seed)
    plan = []
    for c in spec.categories:
        if q[c] == 0:
            continue
        own = np.flatnonzero(cats == c)
        pool, how, labs = own, "own", sorted(set(labels[own])) if labels is not None else []
        if within and c in within:
            pool = own[np.isin(labels[own], list(within[c]))]
            how, labs = "own, restricted", sorted(set(labels[pool]))
        if len(pool) == 0 and fallback and c in fallback:
            pool = np.flatnonzero(np.isin(labels, list(fallback[c])))
            how, labs = "stand-in", sorted(set(labels[pool]))
        if len(pool) == 0:
            raise ValueError(f"category {c!r} has a quota of {q[c]} cells but no cells to draw from")
        plan.append((how == "stand-in", c, pool, how, labs, len(own)))
    plan.sort(key=lambda p: p[0])
    used = np.zeros(adata.n_obs, dtype=np.int64)
    picks, roles, rows = [], [], []
    for standin, c, pool, how, labs, n_own in plan:
        if standin:
            idx, n_reused = draw_least_used(pool, q[c], used, rng)
            bin_rec = bin_reuse = None
        else:
            bins = z_bins(z[pool])
            bq = bin_quotas(z[pool], bins, tau, q[c])
            parts, n_reused, bin_reuse = [], 0, {}
            for b in sorted(bq):
                got, rep = draw_least_used(pool[bins == b], bq[b], used, rng)
                parts.append(got)
                n_reused += rep
                bin_reuse[str(b)] = rep / bq[b] if bq[b] else 0.0  # a fraction of the bin's quota, no count
            idx = np.concatenate(parts) if parts else np.array([], dtype=np.int64)
            bin_rec = {str(b): bq[b] / q[c] for b in sorted(bq)}  # shares within the category, no count
        picks.append(idx)
        roles.append(np.full(len(idx), c, dtype=object))
        rows.append({"category": c, "target_share": shares[c], "quota": q[c], "own_cells": n_own,
                     "pool": int(len(pool)), "drawn_from": how, "labels": ";".join(labs),
                     "unique": int(len(np.unique(idx))), "reused": n_reused,
                     "realised_share": q[c] / n_cells if n_cells else float("nan"),
                     "advance": {"tau": float(tau), "bin_shares": bin_rec, "bin_reused": bin_reuse}})
    idx = np.concatenate(picks) if picks else np.array([], dtype=np.int64)
    role = np.concatenate(roles) if roles else np.array([], dtype=object)
    perm = rng.permutation(len(idx))
    idx, role = idx[perm], role[perm]
    out = _take(adata, idx)
    out.obs["source_category"] = cats[idx]
    out.obs["source_obs"] = np.asarray(adata.obs_names)[idx]
    out.obs[key] = role
    order = {c: i for i, c in enumerate(spec.categories)}
    rows.sort(key=lambda r: order[r["category"]])
    return out, ReweightReport(n_cells=int(n_cells), seed=_seed_record(seed), rows=rows)


# ---------------------------------------------------------------------------------------------------------------------
# The panel


def median_library(pop: ad.AnnData) -> float:
    """The median library size (total counts per cell) of a drawn population, from its counts (docs/DATA_POLICY.md
    2026-11-04 (f)(3): a share rule's row reports it against E10.5's)."""
    return float(np.median(np.asarray(counts_of(pop).sum(axis=1)).ravel()))


def assemble(pop: ad.AnnData, board: Board, post_cells: Callable[[ad.AnnData], ad.AnnData] | None = None) -> ad.AnnData:
    """The last lines of `t1_scope.build_submission`: the hook, then the board's genes as float32 log1p CP10k."""
    if post_cells is not None:
        pop = post_cells(pop)
    cols = pd.Index(pop.var_names).get_indexer(board.genes)
    if (cols < 0).any():
        missing = [g for g, i in zip(board.genes, cols, strict=True) if i < 0]
        raise KeyError(f"{len(missing)} board genes are missing from the reference, e.g. {missing[:5]}")
    X = pop.X[:, cols]
    return make_submission(X.toarray() if sparse.issparse(X) else np.asarray(X), board.genes)
