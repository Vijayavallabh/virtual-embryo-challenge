"""Row K26's features and map (experiments/harness/TRIES.md; the fold read is experiments/t3_decomp/k26_query_map.py).

A knocked-out gene's CONTEXT is a vector over the T3 panel genes built from the released wild-type T1 files alone
(E8.5 and E9.5, whole transcriptome): the gene's within-cell-type co-expression with each panel gene. Nothing here
reads a knockout file, a fold target or a score; the map is fitted on log tables handed to it.

Features (fixed in the row):
  * T1 E8.5 and E9.5 pooled, cell types as annotated (obs `celltype`, pooled by name), types with >= MIN_TYPE_CELLS;
  * per type: the Pearson correlation of log1p CP10k of gene G with each panel gene, after regressing both on
    ln(total counts) (OLS with an intercept, inside the type);
  * the average over types weighted by type cell count (types where G or the panel gene does not vary inside the type
    are left out of that gene's average, the weights renormalised);
  * the vector residualised across panel genes on each panel gene's pooled mean level (mean log1p CP10k) and pooled
    detection rate (share of cells with a count) by OLS with an intercept; then z-scored across panel genes. Panel
    genes absent from T1, or with no defined correlation, are left out of both fits and set to 0 after the z-score;
  * a gene detected in fewer than DETECTION_FLOOR of the pooled cells has no context (refused).

The map (fixed in the row): reduced-rank ridge from context (n x p) to log table (n x g), rank RANK, lambda =
LAMBDA_PER x n, intercept the training mean table: with Yc = Y - mean(Y) (the contexts are not centred, so the
intercept is exactly the mean table), B_ridge = X' (X X' + lambda I)^-1 Yc, B = B_ridge V_k V_k' with V_k the top-k
right singular vectors of X B_ridge, and the prediction at x is mean(Y) + x B; at a context orthogonal to every
training context it is the mean table. No tuning.

Stage 0's statistics (leave_one_out) and the detection deciles used for random queries live here as well, as pure
functions on arrays.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

MIN_TYPE_CELLS = 100
DETECTION_FLOOR = 0.01
RANK = 2
LAMBDA_PER = 10.0
N_RANDOM = 20
DECILES = 10
T1_FILES = ("T1/E8.5_RNA.h5ad", "T1/E9.5_RNA.h5ad")
TYPE_KEY = "celltype"


class ContextRefused(Exception):
    """A gene that has no context (absent from T1, or under the detection floor)."""


# ---------------------------------------------------------------------------------------------------------------------
# The expression data


@dataclass
class Expression:
    """Pooled wild-type cells: X (cells x genes, log1p CP10k, CSC or dense), ln total counts per cell, each cell's
    type, the gene names."""

    X: object
    ln_total: np.ndarray
    types: np.ndarray
    genes: list[str]

    def __post_init__(self):
        self.ln_total = np.asarray(self.ln_total, dtype=np.float64)
        self.types = np.asarray(self.types).astype(str)
        self.genes = [str(g) for g in self.genes]
        self.pos = {g: i for i, g in enumerate(self.genes)}
        if len(self.pos) != len(self.genes):
            raise ValueError("duplicated gene names in the T1 expression")
        if self.X.shape != (self.ln_total.size, len(self.genes)) or self.types.size != self.ln_total.size:
            raise ValueError("the expression's shapes disagree")

    def columns(self, idx) -> np.ndarray:
        """Dense float64 columns (cells x len(idx))."""
        sub = self.X[:, list(idx)]
        return np.asarray(sub.toarray() if hasattr(sub, "toarray") else sub, dtype=np.float64)


def load_t1(data: Path) -> Expression:
    """T1 E8.5 and E9.5 through vec.prep.load_stage (counts recovered from .X; no lineage columns), pooled on the
    genes of the first file (both files must hold the same genes in the same order)."""
    from scipy import sparse

    from vec.prep import load_stage

    parts = [load_stage(Path(data) / f, lineage=False) for f in T1_FILES]
    genes = list(parts[0].var_names)
    for a in parts[1:]:
        if list(a.var_names) != genes:
            raise ValueError("the T1 files hold different genes or another order")
    X = sparse.vstack([sparse.csr_matrix(a.X) for a in parts]).tocsc()
    tot = np.concatenate([np.asarray(a.layers["counts"].sum(1)).ravel() for a in parts]).astype(np.float64)
    if np.any(tot <= 0):
        raise ValueError("a T1 cell has no counts")
    types = np.concatenate([a.obs[TYPE_KEY].astype(str).to_numpy() for a in parts])
    return Expression(X=X, ln_total=np.log(tot), types=types, genes=genes)


def detection(expr: Expression) -> np.ndarray:
    """Per gene, the share of pooled cells with a non-zero value."""
    X = expr.X
    if hasattr(X, "getnnz"):
        nz = np.asarray((X != 0).sum(0)).ravel()
    else:
        nz = (np.asarray(X) != 0).sum(0)
    return nz.astype(np.float64) / X.shape[0]


def mean_level(expr: Expression) -> np.ndarray:
    return np.asarray(expr.X.mean(0)).ravel().astype(np.float64)


def kept_types(types: np.ndarray, min_cells: int = MIN_TYPE_CELLS) -> list[tuple[str, np.ndarray]]:
    """(type, row indices) for the types with at least min_cells cells, sorted by name."""
    types = np.asarray(types).astype(str)
    out = []
    for t in sorted(set(types)):
        rows = np.flatnonzero(types == t)
        if rows.size >= min_cells:
            out.append((t, rows))
    return out


def residual_on(Z: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Y minus its OLS fit on the columns of Z (Z holds the intercept column)."""
    coef, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    return Y - Z @ coef


def type_correlations(xq: np.ndarray, XP: np.ndarray, ln_total: np.ndarray) -> np.ndarray:
    """Within one type: corr(residual of each query column, residual of each panel column), both regressed on
    [1, ln total]; (q x p) with NaN where a residual has no variance."""
    Z = np.column_stack([np.ones_like(ln_total), ln_total])
    RQ = residual_on(Z, np.asarray(xq, dtype=np.float64))
    RP = residual_on(Z, np.asarray(XP, dtype=np.float64))
    nq = np.sqrt((RQ ** 2).sum(0))
    npn = np.sqrt((RP ** 2).sum(0))
    tol = 1e-12 * max(1.0, float(np.sqrt(len(ln_total))))
    with np.errstate(invalid="ignore", divide="ignore"):
        C = (RQ.T @ RP) / np.outer(nq, npn)
    C[nq <= tol, :] = np.nan
    C[:, npn <= tol] = np.nan
    return C


def raw_contexts(expr: Expression, query: list[str], panel: list[str], min_cells: int = MIN_TYPE_CELLS
                 ) -> tuple[np.ndarray, dict]:
    """The type-count-weighted average correlation (len(query) x len(panel)), NaN where undefined or the panel gene is
    absent from T1; and a record (types kept, their cell counts)."""
    qi = [expr.pos[g] for g in query]
    present = [g in expr.pos for g in panel]
    pi = [expr.pos[g] for g in panel if g in expr.pos]
    num = np.zeros((len(query), len(pi)))
    wsum = np.zeros((len(query), len(pi)))
    rec = {"types": {}}
    D = expr.columns(qi + pi)  # the needed columns once, dense (cells x (q + p))
    for t, rows in kept_types(expr.types, min_cells):
        C = type_correlations(D[rows, :len(qi)], D[rows, len(qi):], expr.ln_total[rows])
        ok = np.isfinite(C)
        num += np.where(ok, C, 0.0) * rows.size
        wsum += ok * rows.size
        rec["types"][t] = int(rows.size)
    with np.errstate(invalid="ignore", divide="ignore"):
        avg = np.where(wsum > 0, num / np.where(wsum > 0, wsum, 1.0), np.nan)
    out = np.full((len(query), len(panel)), np.nan)
    out[:, np.flatnonzero(present)] = avg
    return out, rec


def finish_context(raw: np.ndarray, level: np.ndarray, det: np.ndarray) -> np.ndarray:
    """Residualise each row across panel genes on [1, level, det] (OLS over the defined entries), z-score, and set the
    undefined entries to 0."""
    raw = np.atleast_2d(np.asarray(raw, dtype=np.float64))
    level = np.asarray(level, dtype=np.float64)
    det = np.asarray(det, dtype=np.float64)
    out = np.zeros_like(raw)
    for i, r in enumerate(raw):
        ok = np.isfinite(r) & np.isfinite(level) & np.isfinite(det)
        if ok.sum() < 4:
            raise ContextRefused(f"row {i}: fewer than 4 panel genes with a defined correlation")
        Z = np.column_stack([np.ones(ok.sum()), level[ok], det[ok]])
        e = residual_on(Z, r[ok])
        sd = float(np.std(e))
        if not sd > 0:
            raise ContextRefused(f"row {i}: the residualised context has no variance")
        out[i, ok] = (e - e.mean()) / sd
    return out


def contexts(expr: Expression, query: list[str], panel: list[str], det: np.ndarray | None = None,
             level: np.ndarray | None = None, floor: float = DETECTION_FLOOR,
             min_cells: int = MIN_TYPE_CELLS) -> tuple[dict[str, np.ndarray], dict]:
    """{gene: context (len(panel),)} for each query gene; refuses (ContextRefused) a gene absent from T1 or detected in
    fewer than `floor` of the cells. The panel genes' level and detection are the pooled ones."""
    det = detection(expr) if det is None else np.asarray(det, dtype=np.float64)
    level = mean_level(expr) if level is None else np.asarray(level, dtype=np.float64)
    for g in query:
        if g not in expr.pos:
            raise ContextRefused(f"{g} is not in the T1 gene list: no context")
        if det[expr.pos[g]] < floor:
            raise ContextRefused(f"{g} is detected in fewer than {floor:.0%} of the T1 cells: refused")
    raw, rec = raw_contexts(expr, list(query), list(panel), min_cells)
    lv = np.array([level[expr.pos[g]] if g in expr.pos else np.nan for g in panel])
    dv = np.array([det[expr.pos[g]] if g in expr.pos else np.nan for g in panel])
    fin = finish_context(raw, lv, dv)
    rec["panel_absent_from_t1"] = int(sum(g not in expr.pos for g in panel))
    rec["breadth"] = {g: {"detection": float(det[expr.pos[g]]),
                          "mean_abs_raw_corr": float(np.nanmean(np.abs(raw[i]))) if np.isfinite(raw[i]).any()
                          else None}
                      for i, g in enumerate(query)}
    return {g: fin[i] for i, g in enumerate(query)}, rec


# ---------------------------------------------------------------------------------------------------------------------
# Detection deciles and random queries


def decile_of(det: np.ndarray, genes: list[str], floor: float = DETECTION_FLOOR, bins: int = DECILES
              ) -> dict[str, int]:
    """Each gene at or above the floor gets its detection decile among those genes (a stable sort by detection, cut
    into `bins` runs of consecutive ranks whose sizes differ by at most one, as harness.t3_decile_permutation cuts)."""
    det = np.asarray(det, dtype=np.float64)
    elig = np.flatnonzero(det >= floor)
    order = elig[np.argsort(det[elig], kind="stable")]
    out = {}
    for b, chunk in enumerate(np.array_split(order, bins)):
        for i in chunk:
            out[genes[i]] = b
    return out


def draw_same_decile(deciles: dict[str, int], gene: str, n: int, rng: np.random.Generator,
                     exclude: set[str]) -> list[str]:
    """n distinct genes of `gene`'s decile, none in `exclude` nor `gene` itself, drawn without replacement from the
    decile's genes in sorted order."""
    if gene not in deciles:
        raise ContextRefused(f"{gene} has no detection decile (absent or under the floor)")
    pool = sorted(g for g, b in deciles.items() if b == deciles[gene] and g != gene and g not in exclude)
    if len(pool) < n:
        raise ContextRefused(f"{gene}'s decile holds {len(pool)} eligible genes, fewer than {n}")
    return [pool[i] for i in rng.choice(len(pool), size=n, replace=False)]


# ---------------------------------------------------------------------------------------------------------------------
# The map


@dataclass
class Map:
    y_mean: np.ndarray
    B: np.ndarray
    rank: int
    lam: float
    n: int

    def predict(self, x) -> np.ndarray:
        return self.y_mean + np.asarray(x, dtype=np.float64) @ self.B


def fit_map(X, Y, rank: int = RANK, lambda_per: float = LAMBDA_PER) -> Map:
    """The reduced-rank ridge map of the row (see the module): lambda = lambda_per x n, rank min(rank, numerical rank
    of the ridge fit)."""
    X = np.atleast_2d(np.asarray(X, dtype=np.float64))
    Y = np.atleast_2d(np.asarray(Y, dtype=np.float64))
    n = X.shape[0]
    if n < 2 or Y.shape[0] != n:
        raise ValueError("the map needs at least two training perturbations, one table each")
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(Y))):
        raise ValueError("contexts and tables must be finite")
    lam = float(lambda_per) * n
    ym = Y.mean(0)
    Yc = Y - ym
    B_ridge = X.T @ np.linalg.solve(X @ X.T + lam * np.eye(n), Yc)
    F = X @ B_ridge
    _, s, Vt = np.linalg.svd(F, full_matrices=False)
    k = int(min(rank, int((s > 1e-12 * max(1.0, float(s[0]) if s.size else 0.0)).sum())))
    B = B_ridge @ Vt[:k].T @ Vt[:k] if k > 0 else np.zeros_like(B_ridge)
    return Map(y_mean=ym, B=B, rank=k, lam=lam, n=n)


# ---------------------------------------------------------------------------------------------------------------------
# Stage 0


def spearman(a, b) -> float | None:
    from scipy.stats import spearmanr

    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.size < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return None
    return float(spearmanr(a, b).statistic)


GATE_SHARE = (2, 3)  # at least two thirds of the held-out perturbations
GATE_MARGIN = 0.05  # median Spearman margin above this
MIN_HELD = 3  # fewer held-out perturbations: stage 0 cannot pass


def leave_one_out(ctx: dict[str, np.ndarray], labels: dict[str, np.ndarray], random_queries: dict[str, list[str]],
                  random_ctx: dict[str, np.ndarray]) -> dict:
    """Stage 0 over the perturbations in `labels` (the outside genes only): per held-out gene h, the map fitted on the
    others and evaluated at h's context, against (i) the others' mean table and (ii) the same map at each of h's random
    queries (the mean of their Spearman correlations). Returns per-gene figures and the gate."""
    genes = sorted(labels)
    per = {}
    for h in genes:
        train = [g for g in genes if g != h]
        y = labels[h]
        row = {"n_train": len(train), "rho_map": None, "rho_mean": None, "rho_rand": None}
        if len(train) >= 2:
            m = fit_map(np.vstack([ctx[g] for g in train]), np.vstack([labels[g] for g in train]))
            row["rho_map"] = spearman(m.predict(ctx[h]), y)
            row["rho_mean"] = spearman(m.y_mean, y)
            rr = [spearman(m.predict(random_ctx[q]), y) for q in random_queries[h]]
            rr = [v for v in rr if v is not None]
            row["rho_rand"] = float(np.mean(rr)) if rr else None
            row["rank"] = m.rank
        per[h] = row
    return {"per_gene": per, **gate(per)}


def gate(per: dict) -> dict:
    """The map beats BOTH controls, each on at least two thirds of the held-out perturbations (strictly higher rho) and
    by a median margin above GATE_MARGIN; an undefined rho counts as not beaten, an undefined margin as -inf."""
    n = len(per)
    out = {"n_held": n, "share_needed": f"{GATE_SHARE[0]}/{GATE_SHARE[1]}", "margin_needed": GATE_MARGIN}
    ok = n >= MIN_HELD
    for c in ("mean", "rand"):
        marg = []
        for r in per.values():
            a, b = r.get("rho_map"), r.get(f"rho_{c}")
            marg.append(a - b if a is not None and b is not None else -np.inf)
        marg = np.array(marg, dtype=np.float64)
        wins = int((marg > 0).sum())
        med = float(np.median(marg)) if n else None
        share_ok = wins * GATE_SHARE[1] >= GATE_SHARE[0] * n and n > 0
        marg_ok = med is not None and np.isfinite(med) and med > GATE_MARGIN
        out[f"vs_{c}"] = {"wins": wins, "median_margin": (med if med is None or np.isfinite(med) else None),
                          "share_ok": bool(share_ok), "margin_ok": bool(marg_ok)}
        ok = ok and share_ok and marg_ok
    out["passes"] = bool(ok)
    if n < MIN_HELD:
        out["why"] = f"{n} held-out perturbations (< {MIN_HELD}): stage 0 cannot pass"
    return out
