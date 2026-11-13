"""The shared transport-then-adjust core (docs/PLAN.md section 5): the T1 and T2 heads build their predictions
through it, and the pre-registered T3 script runs identical copies of its T3 primitives.

Task heads hold what differs between boards and pass it in: the pools of real cells, the target composition, the
expression factors and, on the coordinate boards, the rules for positions. The core turns them into cells. It knows no
stage, cell type or board: no size, composition or class list of any stage is written here (tests/test_core.py
checks the module's numeric literals against an allow-list and its strings against the released class lists).

Steps implemented now:

1. Allocate a composition to whole cell counts: `quotas` (shares by name; floors, then the largest remainders, ties by
   name, and a share of 0 never gets a cell) and `alloc` (an array of shares; floors, then the largest remainders).
2. Draw real cells from pools, under one of three draw policies:
   - `draw_least_used`: every unused row before any repeat, then whole passes over the pool, each without replacement;
     the use counts carry over between calls, so a pool that overlaps rows drawn before takes its unused rows first;
   - `draw_balanced`: random rows in random order, without replacement up to the pool size and above it each pool row
     n // N or n // N + 1 times;
   - `draw_then_repeat`: without replacement up to the pool size; above it the whole pool once, then the rest drawn
     with replacement.
3. Edit expression in count space, so counts stay non-negative integers (`counts_of` recovers them, `mean_cp10k` gives a
   population's mean CP10k per gene), and a cell the edit would empty keeps its counts:
   - `apply_factors`: per-(type, gene) factors, binomial thinning below 1 and a Poisson top-up of the cell's own counts
     above 1 (exact zeros stay zero), each cell's expected depth kept on request;
   - `apply_gene_factor`: a per-gene factor, binomial thinning below 1 and a Poisson top-up from smoothed proportions
     above 1;
   - `rescale_depth`: a per-cell depth change, `thin` (binomial thinning) below 1 and a Poisson top-up from smoothed
     proportions above 1.
4. log1p CP10k: `log1p_cp10k` (vec.prep's, the released convention).
5. Couple two populations with positions and move cells between them (feat-010, 2026-11-01): `fgw_coupling` (a
   balanced fused Gromov-Wasserstein coupling, POT's exact solver, each cost scaled to a maximum of 1),
   `coupling_rotation` (the proper rotation that best maps one centred cloud onto the other under a coupling),
   `barycentric` (each row's coupling-weighted partner) and `nearest_anchors` (k nearest anchors of each point). The
   T2 head turns them into displacement interpolation (`vec.t2_build.interpolate_clouds`). `aligned_coupling`
   (2026-11-03, embryo try E10 of experiments/harness/TRIES.md): an exact transport of two clouds already in one frame
   on a fused cost of features and positions, which tells apart the mirror-image sides that the within-cloud distances
   and left-right symmetric features cannot.
6. Herd a draw (2026-11-03, T1 try T5 of experiments/harness/TRIES.md): `herd_rows` swaps drawn rows for undrawn rows
   of the same pool when the swap brings the drawn rows' mean expression closer to the pool's, a variance reduction of
   the sample's pseudobulk with no change to what is drawn from.

Steps PLAN section 5 plans, which later features add here: unbalanced coupling of populations without positions
(feat-015); births from the best legal source (feat-014); leave-one-transition-out calibration of a few shrinkage
parameters (feat-008).

The task heads and the boundary (a head decides what goes in; the core does the work on cells):
- T1, `vec.t1_scope`: chamber categories, target shares, stand-in labels and the momentum factors (damped log ratios of
  mean CP10k between two stages). It calls `quotas`, `draw_least_used`, `counts_of`, `mean_cp10k`, `apply_factors` and
  `log1p_cp10k`.
- T2, `vec.t2_build`: the count rule, the coordinate head (duplicate jitter, proper rotations, the pre-registered size
  procedure) and the cloud guard. It calls `draw_balanced`, for the submitted rows and for the guard's reference draws.
- T3, `experiments/t3_decomp`: the genotype head (compositions, the per-gene factor, depth ratios, stage mixing). Its
  script is pre-registered and hash-pinned (build_submission.py refuses a decision file whose script sha256 differs
  from the file on disk), so it keeps its own primitives, written there first, and the core holds identical copies:
  `alloc`, `draw_rows` (here `draw_then_repeat`), `thin`, `rescale_depth` and `apply_gene_factor`. tests/test_core.py
  pins the two together: the same bodies, and the same outputs and generator states for the same inputs and seeds.
  The script's `lognorm` is not a copy: it computes log1p CP10k in another order and agrees with `log1p_cp10k` within
  one float32 ulp.

The reference baselines (`vec.baselines`, `vec baseline`) sit outside the core: they make floor rows and practice
files, never a prediction submitted for a score.

The functions that moved here keep their old names in the heads, as the same objects: `vec.t1_scope.quotas`, `_draw`
(`draw_least_used`), `counts_of`, `mean_cp10k` and `apply_factors`, and `vec.t2_build._balanced_rows`
(`draw_balanced`).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import anndata as ad
import numpy as np
from scipy import sparse

from vec.prep import TARGET_SUM, log1p_cp10k, recover_counts

__all__ = ["aligned_coupling", "alloc", "apply_factors", "apply_gene_factor", "barycentric", "counts_of",
           "coupling_rotation", "draw_balanced", "draw_least_used", "fgw_coupling", "nearest_anchors",
           "draw_then_repeat", "log1p_cp10k", "mean_cp10k", "quotas", "rescale_depth", "thin"]

Seed = int | Sequence[int]


# ---------------------------------------------------------------------------------------------------------------------
# 1. Allocation


def quotas(shares: Mapping[str, float], n: int) -> dict[str, int]:
    """Whole cell counts per category summing to n: floors of share x n, then the largest remainders get one more
    (ties by category name). A category with share 0 never gets a cell."""
    if n < 0:
        raise ValueError(f"n must be >= 0, got {n}")
    raw = {c: float(s) * n for c, s in shares.items()}
    q = {c: int(np.floor(v + 1e-9)) for c, v in raw.items()}
    left = n - sum(q.values())
    for c in sorted(raw, key=lambda c: (-(raw[c] - q[c]), c))[:max(left, 0)]:
        q[c] += 1
    return q


def alloc(n, f) -> np.ndarray:
    """Largest-remainder allocation of n cells to shares f (floors, then one more for the largest remainders, ties in
    np.argsort's order). Frozen copy: experiments/t3_decomp/t3_decompose.py `alloc`."""
    raw = n * np.asarray(f, dtype=np.float64)
    k = np.floor(raw).astype(int)
    rem = int(n - k.sum())
    if rem > 0:
        k[np.argsort(-(raw - k))[:rem]] += 1
    return k


# ---------------------------------------------------------------------------------------------------------------------
# 2. Draws from pools of real cells


def draw_least_used(pool: np.ndarray, k: int, used: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, int]:
    """k rows from `pool` (unique row indices), least-used first: every unused row before any repeat, then whole
    passes over the pool, each without replacement. Updates `used`; returns (rows, picks that repeat a cell)."""
    pool = np.asarray(pool)
    before = used[pool].copy()
    taken = np.zeros(len(pool), dtype=np.int64)
    parts, rest = [], k
    while rest > 0:
        level = before + taken
        cand = np.flatnonzero(level == level.min())
        sel = rng.choice(cand, min(rest, len(cand)), replace=False)
        taken[sel] += 1
        parts.append(pool[sel])
        rest -= len(sel)
    used[pool] += taken
    picks = np.concatenate(parts) if parts else np.array([], dtype=np.int64)
    return picks, int(k - np.count_nonzero((before == 0) & (taken > 0)))


def draw_balanced(n_pool: int, n: int, rng: np.random.Generator) -> np.ndarray:
    """n row indices into a pool of n_pool, in random order: without replacement when n <= n_pool, otherwise every row
    n // n_pool times plus n % n_pool distinct rows once more (a plain draw with replacement would leave about 1/e of
    the pool out at n = n_pool and repeat other rows several times)."""
    k, r = divmod(n, n_pool)
    rows = np.concatenate([np.tile(np.arange(n_pool), k), rng.choice(n_pool, r, replace=False)])
    return rng.permutation(rows).astype(np.int64)


def draw_then_repeat(pool_idx, k, rng) -> np.ndarray:
    """k rows of `pool_idx`: without replacement when k <= len(pool_idx); otherwise the whole pool once, in pool order,
    then k - len(pool_idx) rows drawn with replacement. Empty for k <= 0 or an empty pool. Frozen copy:
    experiments/t3_decomp/t3_decompose.py `draw_rows`."""
    m = len(pool_idx)
    if k <= 0 or m == 0:
        return np.array([], dtype=int)
    if k <= m:
        return rng.choice(pool_idx, k, replace=False)
    return np.concatenate([pool_idx, rng.choice(pool_idx, k - m, replace=True)])


# ---------------------------------------------------------------------------------------------------------------------
# 3. Count space


def counts_of(adata: ad.AnnData, genes: Sequence[str] | None = None) -> sparse.csr_matrix:
    """Integer counts (CSR int64) in the order of `genes` (default: the AnnData's own), from layers["counts"] or
    recovered exactly from the log1p CP10k `.X` (`vec.prep.recover_counts`)."""
    if "counts" in adata.layers:
        C = sparse.csr_matrix(adata.layers["counts"])
        if C.dtype.kind == "f":
            if C.nnz and float(np.abs(C.data - np.rint(C.data)).max()) > 1e-6:
                raise ValueError("layers['counts'] holds non-integer values")
            C = sparse.csr_matrix((np.rint(C.data), C.indices, C.indptr), shape=C.shape)
        C = C.astype(np.int64)
    else:
        C, _, resid = recover_counts(adata.X)
        if resid > 1e-2:
            raise ValueError(f"counts are not recoverable from X (largest residual {resid:.3g})")
        C = C.astype(np.int64)
    if genes is not None and list(adata.var_names) != list(genes):
        idx = adata.var_names.get_indexer(list(genes))
        if (idx < 0).any():
            missing = [g for g, i in zip(genes, idx, strict=True) if i < 0]
            raise KeyError(f"{len(missing)} genes missing, e.g. {missing[:5]}")
        C = C[:, idx]
    C.eliminate_zeros()
    return C


def mean_cp10k(C) -> np.ndarray:
    """Mean over cells of 1e4 * c / (cell total), per gene."""
    C = sparse.csr_matrix(C, dtype=np.float64)
    lib = np.asarray(C.sum(1)).ravel()
    scale = np.divide(TARGET_SUM, lib, out=np.zeros_like(lib), where=lib > 0)
    return np.asarray((sparse.diags(scale) @ C).mean(0)).ravel()


def apply_factors(C, factors: np.ndarray, type_index: np.ndarray, held: np.ndarray | None = None,
                  keep_depth: bool = True, seed: Seed = 0,
                  kernel: str = "independent") -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Edit integer counts by per-(type, gene) factors in count space.

    C: (cells x genes) counts; factors: (types x genes), > 0; type_index: each cell's row of `factors`; held: genes
    left exactly as they are. With keep_depth, each cell's non-held factors are multiplied by the constant that keeps
    its expected total count unchanged. A count c with factor f < 1 becomes Binomial(c, f), with f > 1 c + Poisson
    ((f - 1) c), with f == 1 stays c. Zeros stay zero. A cell left with no counts gets its original counts back.
    Returns (counts as CSR int64, per-cell depth constant, restored-cell mask).

    kernel (default "independent": the per-count draws above, every existing output unchanged). A cell's moved entries
    M are its non-held entries with a count above zero before the step and f != 1 after the depth constant; S is the
    sum of c f over M (the step's expected total there). Entries outside M stay as they are under every kernel.
    "exact_total" (heart arm b0_trend_xt, vec.harness.B0_TREND_XT_RULE; the repair form): the independent draws above,
    unchanged and on the same generator (output Y); then per cell, with d = (sum of Y over M) - rint(S): if d > 0, d
    units are removed from Y on M by one multivariate hypergeometric draw over Y's units there (every unit equally
    likely, without replacement); if d < 0, -d units are added on M by one multinomial draw with probabilities c f / S;
    if d == 0 or M is empty, nothing is drawn. The repair's draws follow all the independent draws on the same stream,
    cells in row order, so a step where every cell has d == 0 gives the independent kernel's output bit for bit. Each
    cell's total over M is then rint(S) (unless the restore rule below brings its original counts back).
    "multinomial" (heart arm b0_trend_mn, vec.harness.B0_TREND_MN_RULE; a whole-cell resampling, not a repair): per
    cell the counts on M are redrawn together as Multinomial(rint(S), c f / S); each count's expectation is c f x
    rint(S) / S. With keep_depth every free nonzero count with f != 1 is in M, so this resamples every such count of a
    shifted cell and adds a variance of about c per count even where f is about 1. One stream, cells in row order.
    Under every kernel zeros stay zero and the restore rule is the same.
    """
    C = sparse.csr_matrix(C, dtype=np.int64, copy=True)
    C.sort_indices()
    n = C.shape[0]
    factors = np.asarray(factors, dtype=np.float64)
    if not np.all(np.isfinite(factors)) or (factors <= 0).any():
        raise ValueError("factors must be finite and positive")
    type_index = np.asarray(type_index)
    held = np.zeros(C.shape[1], dtype=bool) if held is None else np.asarray(held, dtype=bool)
    row_of = np.repeat(np.arange(n), np.diff(C.indptr))
    c = C.data.astype(np.float64)
    f = factors[type_index[row_of], C.indices]
    free = ~held[C.indices]
    scale = np.ones(n)
    if keep_depth:
        num = np.bincount(row_of, weights=c * free, minlength=n)
        den = np.bincount(row_of, weights=c * f * free, minlength=n)
        scale = np.divide(num, den, out=np.ones(n), where=den > 0)
        f = np.where(free, f * scale[row_of], 1.0)
    else:
        f = np.where(free, f, 1.0)
    if kernel not in ("independent", "exact_total", "multinomial"):
        raise ValueError(f"unknown kernel {kernel!r}; known: independent, exact_total, multinomial")
    rng = np.random.default_rng(seed)
    new = C.data.copy()
    if kernel == "multinomial":
        moved = (f != 1.0) & (C.data > 0)
        w = np.where(moved, c * f, 0.0)
        for i in np.unique(row_of[moved]):
            lo, hi = C.indptr[i], C.indptr[i + 1]
            m = moved[lo:hi]
            wi = w[lo:hi][m]
            s = float(wi.sum())
            new[lo:hi][m] = rng.multinomial(int(np.rint(s)), wi / s)
    else:
        down, up = f < 1.0, f > 1.0
        new[down] = rng.binomial(C.data[down], f[down])
        new[up] = C.data[up] + rng.poisson((f[up] - 1.0) * C.data[up])
        if kernel == "exact_total":
            _repair_totals(new, C, c * f, (f != 1.0) & (C.data > 0), row_of, rng)
    dead = np.bincount(row_of, weights=new, minlength=n) == 0
    if dead.any():
        back = dead[row_of]
        new[back] = C.data[back]
    out = sparse.csr_matrix((new, C.indices.copy(), C.indptr.copy()), shape=C.shape)
    out.eliminate_zeros()
    return out, scale, dead


def _repair_totals(new: np.ndarray, C: sparse.csr_matrix, w: np.ndarray, moved: np.ndarray, row_of: np.ndarray,
                   rng: np.random.Generator) -> np.ndarray:
    """apply_factors' "exact_total" repair, in place on `new` (C's data order): per cell, d = (sum of new over the
    moved entries) - rint(sum of w over them); d > 0 removes d units by one multivariate hypergeometric draw over the
    moved entries' units, d < 0 adds -d units by one multinomial draw with probabilities w / sum(w) over them. Cells
    with d == 0 draw nothing; cells in row order. Returns the per-cell d (int64)."""
    n = C.shape[0]
    target = np.rint(np.bincount(row_of[moved], weights=w[moved], minlength=n)).astype(np.int64)
    have = np.bincount(row_of[moved], weights=new[moved], minlength=n).astype(np.int64)
    d = have - target
    for i in np.flatnonzero(d):
        lo, hi = C.indptr[i], C.indptr[i + 1]
        idx = lo + np.flatnonzero(moved[lo:hi])
        if d[i] > 0:
            new[idx] -= rng.multivariate_hypergeometric(new[idx], int(d[i]))
        else:
            wi = w[idx]
            new[idx] += rng.multinomial(int(-d[i]), wi / wi.sum())
    return d


def thin(C, p, rng) -> np.ndarray:
    """Binomial thinning with per-cell keep probability p; a cell that would become empty keeps its counts (dense
    counts, cells x genes). Frozen copy: experiments/t3_decomp/t3_decompose.py `thin`."""
    p = np.clip(np.broadcast_to(np.asarray(p, dtype=np.float64).reshape(-1, 1), (C.shape[0], 1)), 0.0, 1.0)
    out = rng.binomial(C.astype(np.int64), p).astype(np.int32)
    empty = out.sum(1) == 0
    out[empty] = C[empty]
    return out


def rescale_depth(C, ratio, P, rng) -> np.ndarray:
    """Per-cell depth change by `ratio`: thinning below 1, Poisson top-up from smoothed proportions P above 1. Frozen
    copy: experiments/t3_decomp/t3_decompose.py `rescale_depth`."""
    ratio = np.asarray(ratio, dtype=np.float64)
    out = thin(C, np.minimum(ratio, 1.0), rng)
    up = ratio > 1.0
    if up.any():
        n = C[up].sum(1).astype(np.float64)
        out[up] = out[up] + rng.poisson(P[up].astype(np.float64) * (n * (ratio[up] - 1.0))[:, None]).astype(np.int32)
    return out


def apply_gene_factor(C, phi, P, rng) -> np.ndarray:
    """Per-gene factor in count space: binomial thinning where phi < 1, Poisson top-up from P where phi > 1; a cell left
    empty keeps its counts. Frozen copy: experiments/t3_decomp/t3_decompose.py `apply_gene_factor`."""
    out = C.copy()
    down = phi < 1.0
    if down.any():
        out[:, down] = rng.binomial(C[:, down].astype(np.int64), phi[down]).astype(np.int32)
    up = phi > 1.0
    if up.any():
        n = C.sum(1, keepdims=True).astype(np.float64)
        out[:, up] += rng.poisson(P[:, up].astype(np.float64) * n * (phi[up] - 1.0)).astype(np.int32)
    empty = out.sum(1) == 0
    out[empty] = C[empty]
    return out


# ---------------------------------------------------------------------------------------------------------------------
# 5. Coupling and displacement


def fgw_coupling(features_a, features_b, coords_a, coords_b, alpha: float, max_iter: int,
                 emd_max_iter: int) -> np.ndarray:
    """A balanced fused Gromov-Wasserstein coupling of two populations with positions (POT's exact conditional
    gradient, square loss, uniform masses).

    The feature cost is the squared Euclidean distance between rows of features_a and features_b; the structure costs
    are the Euclidean distances within each cloud, so the coupling ignores how the clouds are placed or rotated. Each
    cost is divided by its maximum, which makes `alpha` (the structure term's weight: 0 features only, 1 structure
    only) free of units. max_iter bounds the conditional-gradient steps, emd_max_iter the network simplex of each
    step's exact transport. Returns the coupling (rows a, columns b; its entries sum to 1)."""
    import ot

    FA, FB = np.asarray(features_a, dtype=np.float64), np.asarray(features_b, dtype=np.float64)
    A, B = np.asarray(coords_a, dtype=np.float64), np.asarray(coords_b, dtype=np.float64)
    if FA.shape[0] != A.shape[0] or FB.shape[0] != B.shape[0]:
        raise ValueError("each population needs one feature row per position")
    if not 0.0 <= alpha <= 1.0:
        raise ValueError(f"alpha must lie in [0, 1], got {alpha}")
    M = ot.dist(FA, FB)
    CA = ot.dist(A, A, metric="euclidean")
    CB = ot.dist(B, B, metric="euclidean")
    for C in (M, CA, CB):
        top = C.max()
        if top > 0:
            C /= top
    p = np.full(A.shape[0], 1.0 / A.shape[0])
    q = np.full(B.shape[0], 1.0 / B.shape[0])
    T = ot.gromov.fused_gromov_wasserstein(M, CA, CB, p, q, loss_fun="square_loss", alpha=alpha, max_iter=max_iter,
                                           numItermaxEmd=emd_max_iter)
    return np.asarray(T, dtype=np.float64)


def aligned_coupling(features_a, features_b, coords_a, coords_b, alpha: float, emd_max_iter: int) -> np.ndarray:
    """An exact optimal transport between two populations whose positions are already in one frame, on a fused cost:
    (1 - alpha) times the squared Euclidean distance between rows of features_a and features_b plus alpha times the
    squared Euclidean distance between positions, each divided by its maximum as in `fgw_coupling` (uniform masses, so
    two populations of the same size are paired one to one). Unlike `fgw_coupling` it reads where the cells sit in
    the common frame, so it tells apart cells that within-cloud distances and features cannot: the two mirror-image
    sides of a body whose expression is left-right symmetric. emd_max_iter bounds the network simplex. Returns the
    coupling (rows a, columns b; its entries sum to 1)."""
    import ot

    FA, FB = np.asarray(features_a, dtype=np.float64), np.asarray(features_b, dtype=np.float64)
    A, B = np.asarray(coords_a, dtype=np.float64), np.asarray(coords_b, dtype=np.float64)
    if FA.shape[0] != A.shape[0] or FB.shape[0] != B.shape[0]:
        raise ValueError("each population needs one feature row per position")
    if A.shape[1] != B.shape[1]:
        raise ValueError("the two populations' positions need the same dimensions")
    if not 0.0 <= alpha <= 1.0:
        raise ValueError(f"alpha must lie in [0, 1], got {alpha}")
    M = ot.dist(FA, FB)
    D = ot.dist(A, B)
    for C in (M, D):
        top = C.max()
        if top > 0:
            C /= top
    p = np.full(A.shape[0], 1.0 / A.shape[0])
    q = np.full(B.shape[0], 1.0 / B.shape[0])
    return np.asarray(ot.emd(p, q, (1.0 - alpha) * M + alpha * D, numItermax=emd_max_iter), dtype=np.float64)


def coupling_rotation(coords_a, coords_b, coupling) -> np.ndarray:
    """The proper rotation R (det +1) minimising sum_ij T_ij |R a_i - b_j|^2 for centred clouds a and b under the
    coupling T (weighted Kabsch): R maps a's frame onto b's, so b's points in a's frame are `coords_b @ R`."""
    A, B = np.asarray(coords_a, dtype=np.float64), np.asarray(coords_b, dtype=np.float64)
    T = np.asarray(coupling, dtype=np.float64)
    if T.shape != (A.shape[0], B.shape[0]):
        raise ValueError(f"coupling shape {T.shape} does not match the clouds ({A.shape[0]}, {B.shape[0]})")
    H = A.T @ T @ B
    U, _, Vt = np.linalg.svd(H)
    D = np.eye(H.shape[0])
    if np.linalg.det(Vt.T @ U.T) < 0:
        D[-1, -1] = -1.0
    return Vt.T @ D @ U.T


def barycentric(coupling, values) -> np.ndarray:
    """Each row's coupling-weighted mean of `values` (one row per column of the coupling)."""
    T = np.asarray(coupling, dtype=np.float64)
    w = T.sum(1)
    if (w <= 0).any():
        raise ValueError("a row of the coupling carries no mass")
    return (T @ np.asarray(values, dtype=np.float64)) / w[:, None]


def nearest_anchors(points, anchors, k: int) -> tuple[np.ndarray, np.ndarray]:
    """The k nearest anchors of each point (Euclidean): (distances, indices), each points x k."""
    from scipy.spatial import cKDTree

    k = min(int(k), np.asarray(anchors).shape[0])
    d, i = cKDTree(np.asarray(anchors, dtype=np.float64)).query(np.asarray(points, dtype=np.float64), k=k)
    return np.asarray(d).reshape(len(points), k), np.asarray(i).reshape(len(points), k)


HERD_PROPOSALS = 20  # random swap proposals per drawn row (a size of the computation; more only matches closer)
HERD_BATCH = 256  # proposals evaluated together; a batch's accepted swaps are applied at once


def herd_rows(X, pool, picked, rng: np.random.Generator, proposals: int | None = None,
              target: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    """Swap rows of `picked` (distinct rows drawn from `pool`) for undrawn rows of `pool` whenever a swap brings the
    mean of the picked rows of X closer, in squared distance over all columns, to `target` (default: the pool's mean).

    Random proposals (HERD_PROPOSALS per picked row unless given), drawn HERD_BATCH at a time on distinct positions
    and distinct undrawn rows, are taken in order and each is accepted only if it improves on the mean as it stands
    after the swaps accepted before it (exact greedy: the distance falls with every accepted swap). The picked rows keep
    their positions in the draw (a swapped-in row takes the place of the row it replaces), so a draw's order and its
    size are unchanged. With repeats, or nothing left to swap in, the draw comes back as it is.

    Returns (rows, record: squared distance before and after, proposals, swaps accepted)."""
    pool, picked = np.asarray(pool, dtype=np.int64), np.asarray(picked, dtype=np.int64).copy()
    m = len(picked)
    out = np.setdiff1d(pool, picked)
    rec = {"rows": int(m), "pool": int(len(pool)), "proposals": 0, "accepted": 0}
    if m == 0 or len(out) == 0 or len(np.unique(picked)) < m:
        return picked, {**rec, "skipped": "repeats or nothing to swap in"}
    X = sparse.csr_matrix(X, dtype=np.float64)
    goal = np.asarray(X[pool].mean(0)).ravel() if target is None else np.asarray(target, dtype=np.float64)
    norms = np.asarray(X.multiply(X).sum(1)).ravel()
    e = np.asarray(X[picked].sum(0)).ravel() / m - goal
    rec["sq_distance_before"] = float(e @ e)
    total = HERD_PROPOSALS * m if proposals is None else int(proposals)
    done = 0
    while done < total:
        b = int(min(HERD_BATCH, m, len(out), total - done))
        pos = rng.choice(m, b, replace=False)
        cand = rng.choice(len(out), b, replace=False)
        r_out, r_in = picked[pos].copy(), out[cand].copy()
        Xo, Xi = X[r_out], X[r_in]
        a, o = Xi @ e, Xo @ e  # each proposal's rows against the current error, kept current below
        cross = np.asarray(Xo.multiply(Xi).sum(1)).ravel()
        q = norms[r_in] + norms[r_out] - cross - cross  # |x_in - x_out|^2
        for k in range(b):
            lin = (a[k] - o[k]) / m  # e.d for d = (x_in - x_out) / m
            if lin + lin + q[k] / (m * m) < 0:
                d = (np.asarray(Xi[k].todense()).ravel() - np.asarray(Xo[k].todense()).ravel()) / m
                e = e + d
                a, o = a + Xi @ d, o + Xo @ d
                picked[pos[k]], out[cand[k]] = r_in[k], r_out[k]
                rec["accepted"] += 1
        done += b
    rec["proposals"], rec["sq_distance_after"] = int(done), float(e @ e)
    return picked, rec
