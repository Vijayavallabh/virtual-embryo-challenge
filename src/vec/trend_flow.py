"""The heart's learned flow kept to what is left once each lineage's mean change is the entry's: `b0_trend` with one
per-cell, per-gene factor D inside each moved lineage (candidates `b0_trend_flow`, `b0_trend_flow_perm`,
`b0_trend_flow_pool`, `b0_trend_flow_pool_perm` of `vec.harness`; heart only).

Cells, coordinates, number, size, the per-lineage trend factor F and b0_trend's streams are b0_trend's. The final
factor on drawn cell i, gene g is F(lineage, g) x D(i, g), applied in b0_trend's one count-space call
(`vec.t1_scope.momentum_shift(..., cell_factors=D)`, the same `apply_factors`, depth rule and stream). D:

- Latent move. The count model and its fit are H10's (`vec.latent_extrap.fold_model`, one fit per task shared by every
  arm and seed, and by `b0_lat_ratio` and `b0_lat_flow`). Each drawn cell of a moved lineage (B0's moved lineages,
  `vec.latent_extrap.b0_moved_types`, which are b0_trend's) has its own latent z (the fold's encoder mean) and a move
  dz over H10's time budget tau = B0_DAMPING x the last gap. The move is centred within the lineage's drawn cells:
  du = dz - mean(dz over the lineage's drawn cells); a lineage whose drawn cells all have the same dz gets du = 0
  exactly.
- Ratio. r(i, g) = decoder mean at z + du over decoder mean at z (both at library 1), the log clipped at +-CLIP as the
  ratio arm clips it (`ratio_log_factors`, the lines of `vec.latent_extrap.ratio_counts` that make the factor; that
  function also applies them, so it cannot be used unchanged). A cell with du = 0 gets r = 1 undecoded.
- Recentring. D(i, g) = r(i, g) / G(lineage, g), with G set so that each gene's expected total over the lineage's
  drawn cells after b0_trend's depth rule is what it is without D: sum_i c F D s'_i = sum_i c F s_i, with c the drawn
  cells' counts, F the trend factor and s_i (s'_i) the per-cell constant `vec.core.apply_factors` uses to keep each
  cell's expected total without (with) D. s' changes with D, so the count-weighted mean of r, G0 = sum_i c r / sum_i c,
  does not give that identity; G starts at G0 and is iterated, G <- G x (expected total with D / without D) with s'
  recomputed each time, until the largest relative departure over the lineage's genes is at most RECENTRE_TOL
  (`recentre_after_depth`). A gene with no count in the lineage, a held-flat gene and every cell of a lineage
  b0_trend does not move keep D = 1.

The four arms differ only in where dz comes from:
1. flow: H10's field exactly (`vec.latent_extrap.fold_flow` and `flow_latent`, every FLOW_* constant unchanged).
2. flow_perm, its placebo: the same, with each lineage's mean development kept and only the within-lineage residuals
   permuted across stages before the pairs are formed (`permute_residuals`): u = z - mean z(lineage, true stage), the
   u's permuted across the stages within each lineage (`stage_permutation`: each stage keeps its number of cells per
   lineage; one permutation drawn from default_rng([FLOW_SEED, PERM_STREAM])), each stage's own lineage mean added
   back; the drawn reference cells are still the cells that are moved, from their own z.
3. flow_pool: every cell's latent centred on the mean of its own lineage at its own stage, u = z - mean z(lineage,
   stage); H10's pairs, weights, network and constants on u in place of z over every transition between consecutive
   input stages; the drawn reference cells' u (their z less their lineage's reference-stage mean) integrated over tau,
   dz that move.
4. flow_pool_perm: 3 with the residuals centred on their true (lineage, stage) means first and then permuted across
   stages within lineage as in 2 (no mean added back, no re-centring); the drawn cells' u is their z less the true
   reference stage's lineage mean.
Every field is fitted once per task with H10's FLOW_SEED (one field per backtest, the same on every build seed).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import anndata as ad
import numpy as np
from scipy import sparse

from vec import latent_extrap as le
from vec.core import counts_of

ARMS = ("flow", "flow_perm", "flow_pool", "flow_pool_perm")
NAMES = {"flow": "b0_trend_flow", "flow_perm": "b0_trend_flow_perm", "flow_pool": "b0_trend_flow_pool",
         "flow_pool_perm": "b0_trend_flow_pool_perm"}
PERM_STREAM = 3  # the placebos' stage permutation: default_rng([FLOW_SEED, 3]); flow_pairs uses streams 1 and 2
RECENTRE_TOL = 1e-12  # float tolerance of the recentring's fixed point (largest relative departure; not a model knob)
RECENTRE_MAX_ITER = 10_000  # the recentring raises if it has not reached RECENTRE_TOL by then
RULE = ("b0_trend x D: D(i, g) = r(i, g) / G(lineage, g), r = decoded mean at z + du over decoded mean at z (log "
        "clipped at +-log 3), du = dz - mean dz over the lineage's drawn cells, G set so that each gene's expected "
        "total over the lineage's drawn cells after b0_trend's depth rule is unchanged (iterated to 1e-12 relative); "
        "dz the move over B0's budget from H10's field (flow), its placebo with the within-lineage residuals permuted "
        "across stages and each stage's lineage mean kept (flow_perm), the field on latents centred per lineage and "
        "stage (flow_pool) or that on the centred residuals permuted across stages (flow_pool_perm)")
_FIELDS: dict = {}  # arm -> {"fold", "key", "objects", "field", "ref_means"}: one field per arm and task


# ------------------------------------------------------------------------------------------------ the placebo and pool
def stage_permutation(labs: Sequence[np.ndarray], seed: Sequence[int]) -> np.ndarray:
    """The stage index of every cell of the concatenated stages (stage 0's cells first) after the stages are permuted
    among the cells of each lineage (sorted lineage order, one generator default_rng(seed)): each stage keeps its
    number of cells per lineage."""
    lab = np.concatenate([np.asarray(x).astype(str) for x in labs])
    stage = np.concatenate([np.full(len(x), k, dtype=np.int64) for k, x in enumerate(labs)])
    rng = np.random.default_rng([int(s) for s in seed])
    out = stage.copy()
    for lin in sorted(set(lab)):
        idx = np.flatnonzero(lab == lin)
        out[idx] = stage[idx][rng.permutation(len(idx))]
    return out


def permute_stages(zs: Sequence[np.ndarray], labs: Sequence[np.ndarray],
                   seed: Sequence[int]) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """The stages' latents and lineages with the stage of every cell permuted within its lineage (`stage_permutation`);
    within a stage the cells keep their concatenation order."""
    Z = np.concatenate([np.asarray(z, dtype=np.float64) for z in zs])
    lab = np.concatenate([np.asarray(x).astype(str) for x in labs])
    new = stage_permutation(labs, seed)
    return [Z[new == k] for k in range(len(zs))], [lab[new == k] for k in range(len(zs))]


def permute_residuals(zs: Sequence[np.ndarray], labs: Sequence[np.ndarray], seed: Sequence[int],
                      add_means: bool) -> tuple[list[np.ndarray], list[np.ndarray], list[dict]]:
    """The placebos' stages: every cell's residual u = z - mean z(lineage, its true stage) goes to the stage
    `stage_permutation` gives it (within its lineage; each stage keeps its number of cells per lineage). With
    add_means the destination stage's true lineage mean is added back (flow_perm: each lineage's mean development is
    kept, only the within-lineage residuals are permuted); without, the stages hold the permuted residuals as they are
    (flow_pool_perm). Within a stage the cells keep their concatenation order. Returns the stages' latents, their
    lineages and each true stage's lineage means."""
    lab = np.concatenate([np.asarray(x).astype(str) for x in labs])
    us, means = centre_by_stage(zs, labs)
    U = np.concatenate(us)
    new = stage_permutation(labs, seed)
    if add_means:
        U = U + np.array([means[k][lin] for k, lin in zip(new, lab, strict=True)]).reshape(U.shape)
    return [U[new == k] for k in range(len(zs))], [lab[new == k] for k in range(len(zs))], means


def centre_by_stage(zs: Sequence[np.ndarray], labs: Sequence[np.ndarray]) -> tuple[list[np.ndarray], list[dict]]:
    """u = z - mean z(lineage, stage) for every cell of every stage; returns the u's and each stage's lineage means."""
    us, means = [], []
    for z, lab in zip(zs, labs, strict=True):
        z = np.asarray(z, dtype=np.float64)
        lab = np.asarray(lab).astype(str)
        u = np.empty_like(z)
        m = {}
        for lin in sorted(set(lab)):
            sel = lab == lin
            m[lin] = z[sel].mean(0)
            u[sel] = z[sel] - m[lin]
        us.append(u)
        means.append(m)
    return us, means


def arm_field(arm: str, fold: le.FoldModel, stages: Sequence[tuple[str, float, ad.AnnData]], pool: ad.AnnData,
              start: ad.AnnData, lineages: Sequence[str]) -> tuple[le.FlowField, dict | None]:
    """The arm's field (fitted once per process per fold model object, stages and lineages) and, for the pooled arms,
    the true reference stage's lineage means the drawn cells are centred on. Arm "flow" is
    `vec.latent_extrap.fold_flow`; the placebos' stages are `permute_residuals`'."""
    if arm == "flow":
        return le.fold_flow(fold, stages, pool, start, lineages), None
    if arm not in ARMS:
        raise ValueError(f"arm must be one of {ARMS}, got {arm!r}")
    key = (tuple((r, float(d)) for r, d, _ in stages), tuple(lineages))
    objects = [s[2] for s in stages]
    hit = _FIELDS.get(arm)
    if (hit is None or hit["fold"] is not fold or hit["key"] != key
            or any(a is not b for a, b in zip(hit["objects"], objects, strict=True))):
        zs, labs, left_out = le._stage_latents(fold, stages, pool, start)
        pooled = arm.startswith("flow_pool")
        ref_means = None
        if arm.endswith("_perm"):
            zs, labs, means = permute_residuals(zs, labs, [le.FLOW_SEED, PERM_STREAM], add_means=not pooled)
            if pooled:
                ref_means = means[-1]
        elif pooled:
            zs, means = centre_by_stage(zs, labs)
            ref_means = means[-1]
        pairs, prec = le.flow_pairs(zs, labs, [d for _, d, _ in stages], lineages)
        ff = le.fit_flow(pairs, lineages)
        ff.record = {**prec, **ff.record, "cells_left_out": left_out, "stages": [r for r, _, _ in stages],
                     "lineages": list(lineages), "permuted": arm.endswith("_perm"),
                     "pooled": arm.startswith("flow_pool")}
        for a in [a for a, h in _FIELDS.items() if h["fold"] is not fold]:
            del _FIELDS[a]
        hit = _FIELDS[arm] = {"fold": fold, "key": key, "objects": objects, "field": ff, "ref_means": ref_means}
    return hit["field"], hit["ref_means"]


def latent_moves(arm: str, fold: le.FoldModel, stages: Sequence[tuple[str, float, ad.AnnData]], pool: ad.AnnData,
                 start: ad.AnnData, types: np.ndarray, moved: Mapping, z_own: np.ndarray,
                 damping: float) -> tuple[np.ndarray, le.FlowField | dict]:
    """dz of every drawn cell (zero outside the moved lineages) and the field's record. Arm "flow": z_new of
    `vec.latent_extrap.flow_latent` less z_own; the others integrate their own field (`vec.latent_extrap.integrate`)
    from z_own (flow_perm) or from u_own = z_own - the lineage's reference-stage mean (the pooled arms)."""
    types = np.asarray(types).astype(str)
    z_own = np.asarray(z_own, dtype=np.float64)
    if arm == "flow":
        z_new, rec = le.flow_latent(fold, stages, pool, start, types, moved, z_own, damping)
        return z_new - z_own, rec
    if stages[-1][2] is not pool or stages[-2][2] is not start:
        raise ValueError("the stages must end with `start` and the reference")
    lineages = sorted(moved)
    ff, ref_means = arm_field(arm, fold, stages, pool, start, lineages)
    tau = float(damping) * (float(stages[-1][1]) - float(stages[-2][1]))
    dz = np.zeros_like(z_own)
    m = np.isin(types, lineages)
    if m.any():
        x0 = z_own[m].copy()
        if ref_means is not None:
            x0 = x0 - np.array([ref_means[t] for t in types[m]])
        li = np.searchsorted(np.asarray(ff.lineages), types[m])
        dz[m] = le.integrate(ff.velocity, x0, li, tau) - x0
    return dz, {**ff.record, "tau_days": tau, "rk4_steps": le.FLOW_RK4_STEPS}


# ------------------------------------------------------------------------------------------------ D
def centre_moves(dz: np.ndarray, types: np.ndarray, lineages: Sequence[str]) -> np.ndarray:
    """du = dz - mean dz over each lineage's drawn cells; a lineage whose cells all have the same dz (and every cell
    outside `lineages`) gets du = 0 exactly."""
    dz = np.asarray(dz, dtype=np.float64)
    types = np.asarray(types).astype(str)
    du = np.zeros_like(dz)
    for lin in lineages:
        sel = types == lin
        d = dz[sel]
        if not len(d) or (d == d[0]).all():
            continue
        du[sel] = d - d.mean(0)
    return du


def ratio_log_factors(model, z: np.ndarray, du: np.ndarray, chunk: int = le.CELL_CHUNK) -> tuple[np.ndarray, int]:
    """log r = log decoded mean(z + du) - log decoded mean(z), both at library 1, a gene whose two means are zero 0,
    clipped at +-CLIP (the lines of `vec.latent_extrap.ratio_counts` that make its factor); a cell with du = 0 gets 0
    undecoded. Returns (cells x genes, how many entries had both means zero)."""
    z, du = np.asarray(z, dtype=np.float64), np.asarray(du, dtype=np.float64)
    moved = np.flatnonzero(np.any(du != 0, axis=1))
    out = None
    both_zero = 0
    for lo in range(0, len(moved), chunk):
        sel = moved[lo:lo + chunk]
        ones = np.ones(len(sel))
        mu_new, _ = model.decode(z[sel] + du[sel], ones)
        mu_own, _ = model.decode(z[sel], ones)
        mu_new, mu_own = np.asarray(mu_new, dtype=np.float64), np.asarray(mu_own, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            lf = np.log(mu_new) - np.log(mu_own)
        zero = (mu_new == 0) & (mu_own == 0)
        both_zero += int(zero.sum())
        lf[zero] = 0.0
        lf = np.clip(lf, -le.CLIP, le.CLIP)
        if out is None:
            out = np.zeros((len(z), lf.shape[1]))
        out[sel] = lf
    return out, both_zero


def recentre(r: np.ndarray, C, types: np.ndarray, lineages: Sequence[str]) -> np.ndarray:
    """D = r / G per lineage, G(g) = sum_i c(i, g) r(i, g) / sum_i c(i, g) over the lineage's cells; D = 1 for a gene
    without a count in the lineage and for every cell outside `lineages`."""
    r = np.asarray(r, dtype=np.float64)
    Cd = sparse.csr_matrix(C).toarray().astype(np.float64) if sparse.issparse(C) else np.asarray(C, dtype=np.float64)
    types = np.asarray(types).astype(str)
    D = np.ones_like(r)
    for lin in lineages:
        sel = types == lin
        if not sel.any():
            continue
        c, rl = Cd[sel], r[sel]
        den = c.sum(0)
        num = (c * rl).sum(0)
        has = den > 0
        G = np.where(has, num / np.where(has, den, 1.0), 1.0)
        D[sel] = np.where(has[None, :], rl / G[None, :], 1.0)
    return D


def depth_scale(C: np.ndarray, Fc: np.ndarray, held: np.ndarray) -> np.ndarray:
    """`vec.core.apply_factors`' per-cell depth constant (keep_depth) for dense counts C and per-cell factors Fc: the
    sum of the free genes' counts over the sum of their counts x factors (1 where that is 0)."""
    free = ~np.asarray(held, dtype=bool)
    num = (C * free).sum(1)
    den = (C * Fc * free).sum(1)
    return np.divide(num, den, out=np.ones(len(C)), where=den > 0)


def depth_departure(C: np.ndarray, Fc: np.ndarray, D: np.ndarray, held: np.ndarray) -> np.ndarray:
    """Per gene, the relative departure of the expected total after the depth rule with D, sum_i c F D s'_i, from the
    one without D, sum_i c F s_i (NaN for a held gene or a gene without a count)."""
    held = np.asarray(held, dtype=bool)
    FD = Fc * np.where(held[None, :], 1.0, D)
    with_d = (C * FD * depth_scale(C, FD, held)[:, None]).sum(0)
    without = (C * Fc * depth_scale(C, Fc, held)[:, None]).sum(0)
    out = np.full(C.shape[1], np.nan)
    ok = (without > 0) & ~held
    out[ok] = with_d[ok] / without[ok] - 1.0
    return out


def recentre_after_depth(r: np.ndarray, C, Fc: np.ndarray, held: np.ndarray, types: np.ndarray,
                         lineages: Sequence[str]) -> tuple[np.ndarray, dict]:
    """D = r / G per lineage with G such that sum_i c F D s'_i = sum_i c F s_i for every free gene with a count (s, s'
    the depth rule's per-cell constants without and with D, `depth_scale`): G starts at `recentre`'s count-weighted
    mean and is multiplied per gene by (1 + departure), s' recomputed, until the largest departure is at most
    RECENTRE_TOL (raises after RECENTRE_MAX_ITER). D = 1 for a gene without a count in the lineage, a held gene and
    every cell outside `lineages`. Returns D and, over lineages and genes, the departure after the depth rule of the
    one-step G (max and p90) and of the final G (max), and the largest number of iterations."""
    r = np.asarray(r, dtype=np.float64)
    Cd = sparse.csr_matrix(C).toarray().astype(np.float64) if sparse.issparse(C) else np.asarray(C, dtype=np.float64)
    Fc = np.asarray(Fc, dtype=np.float64)
    held = np.asarray(held, dtype=bool)
    types = np.asarray(types).astype(str)
    D = recentre(r, Cd, types, lineages)
    D[:, held] = 1.0
    one_step, final, iters = [], [], 0
    for lin in lineages:
        sel = types == lin
        if not sel.any():
            continue
        c, f, Dl = Cd[sel], Fc[sel], D[sel]
        dep = depth_departure(c, f, Dl, held)
        ok = np.isfinite(dep)
        one_step.extend(np.abs(dep[ok]).tolist())
        k = 0
        while ok.any() and np.abs(dep[ok]).max() > RECENTRE_TOL:
            if k == RECENTRE_MAX_ITER:
                raise RuntimeError(f"the recentring of {lin!r} did not reach {RECENTRE_TOL} in {k} iterations")
            Dl[:, ok] = Dl[:, ok] / (1.0 + dep[ok])[None, :]
            dep = depth_departure(c, f, Dl, held)
            k += 1
        iters = max(iters, k)
        final.extend(np.abs(dep[ok]).tolist())
        D[sel] = Dl
    one = np.asarray(one_step)
    return D, {"max_rel_departure_one_step": float(one.max()) if one.size else 0.0,
               "p90_rel_departure_one_step": float(np.quantile(one, 0.9)) if one.size else 0.0,
               "max_rel_departure": float(max(final)) if final else 0.0, "recentre_iterations": int(iters)}


def _cos_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    na, nb = np.linalg.norm(a, axis=1), float(np.linalg.norm(b))
    ok = (na > 0) & (nb > 0)
    out = np.full(len(a), np.nan)
    out[ok] = (a[ok] @ b) / (na[ok] * nb)
    return out


def cell_factors(arm: str, ctx, pool: ad.AnnData, rows: np.ndarray, like: ad.AnnData,
                 damping: float) -> tuple[object, dict]:
    """For b0_trend's drawn cells `like` (typed as `_b0_t2_shift` types them, obs["b0_type"]) of `pool` (rows `rows`):
    the callable `vec.t1_scope.momentum_shift` calls with the trend's per-cell factors and the held mask, which returns
    D (`recentre_after_depth`) and completes the build's record under "check", and that record (shares and ratios
    only)."""
    t = ctx.task
    start = ctx.inputs["start"]
    genes = list(map(str, like.var_names))
    if list(map(str, pool.var_names)) != genes:
        raise ValueError("the drawn cells' genes must be the pool's")
    types = np.asarray(like.obs[le.T2_TYPE_KEY]).astype(str)
    moved = le.b0_moved_types(le._typed(start, "start"), le._typed(pool, "end"), types, le.T2_TYPE_KEY)
    lineages = sorted(moved)
    fold = le.fold_model(pool, start, genes)
    z_own = fold.z_ref[np.asarray(rows, dtype=np.int64)]
    stages = [(r, d, ctx.inputs[r]) for r, d in le.flow_stages(t)]
    dz, frec = latent_moves(arm, fold, stages, pool, start, types, moved, z_own, damping)
    du = centre_moves(dz, types, lineages)
    logr, both_zero = ratio_log_factors(fold.model, z_own, du)
    C = counts_of(like)
    out: dict = {"flow": frec, "fold": fold.record}

    def factors(Fc: np.ndarray, held: np.ndarray) -> np.ndarray:
        if logr is None:
            D, L = np.ones(C.shape), np.zeros(C.shape)
            dep = {"max_rel_departure_one_step": 0.0, "p90_rel_departure_one_step": 0.0, "max_rel_departure": 0.0,
                   "recentre_iterations": 0}
        else:
            D, dep = recentre_after_depth(np.exp(logr), C, Fc, held, types, lineages)
            L = logr
        rec = factor_record(D, L, C, types, lineages, du, dz)
        rec.update(dep)
        rec.update({"arm": arm, "candidate": NAMES[arm], "rule": "vec.trend_flow.RULE",
                    "n_both_zero_means": both_zero, "latent_sha256": fold.record.get("latent_sha256"),
                    "params_sha256": frec.get("params_sha256"), "tau_days": frec.get("tau_days")})
        out["check"] = rec
        return D

    return factors, out


def factor_record(D, logr, C, types, lineages, du, dz) -> dict:
    """Shares and ratios only: the lineages moved, the SD of log D over the entries it acts on (a count in a drawn cell
    of a moved lineage), the share of those entries at the clip, the mean cosine between du and the lineage's
    uncentred mean move, and per lineage the within share of the move (sum ||du||^2 / sum ||dz||^2 over its drawn
    cells, as `vec.latent_extrap.lineage_diagnostics` defines it) with its median. The departure after the depth rule
    is `recentre_after_depth`'s."""
    Cd = sparse.csr_matrix(C).toarray().astype(np.float64)
    types = np.asarray(types).astype(str)
    in_moved = np.isin(types, list(lineages))
    acting = (Cd > 0) & in_moved[:, None]
    L = np.log(D[acting]) if acting.any() else np.zeros(0)
    at_clip = np.abs(logr[acting]) >= le.CLIP - 1e-12 if acting.any() else np.zeros(0, dtype=bool)
    cos, within = [], {}
    for lin in lineages:
        sel = types == lin
        if not sel.any():
            continue
        cos.extend(_cos_rows(du[sel], dz[sel].mean(0)).tolist())
        ss = float((np.asarray(dz[sel]) ** 2).sum())
        within[lin] = float((np.asarray(du[sel]) ** 2).sum()) / ss if ss > 0 else float("nan")
    cos = np.asarray(cos, dtype=np.float64)
    return {"lineages_moved": list(lineages), "sd_log_d": float(L.std()) if L.size else 0.0,
            "share_at_clip": float(at_clip.mean()) if at_clip.size else 0.0,
            "mean_cosine_du_mean_move": float(np.nanmean(cos)) if np.isfinite(cos).any() else float("nan"),
            "within_share": within,
            "median_within_share": float(np.median(list(within.values()))) if within else float("nan")}


def check_line(task_key: str, seed: int, rec: dict) -> str:
    c = rec["check"]
    return (f"check {c['candidate']} task={task_key} seed={seed} lineages={','.join(c['lineages_moved'])} "
            f"sd_log_d={c['sd_log_d']:.4g} share_at_clip={c['share_at_clip']:.4g} "
            f"max_rel_departure={c['max_rel_departure']:.3g} "
            f"max_rel_departure_one_step={c['max_rel_departure_one_step']:.3g} "
            f"p90_rel_departure_one_step={c['p90_rel_departure_one_step']:.3g} "
            f"recentre_iterations={c['recentre_iterations']} mean_cos_du_move={c['mean_cosine_du_mean_move']:.4g} "
            f"median_within_share={c['median_within_share']:.4g} "
            f"latent_sha256={c['latent_sha256']} params_sha256={c['params_sha256']}")


__all__ = ["ARMS", "NAMES", "PERM_STREAM", "RECENTRE_MAX_ITER", "RECENTRE_TOL", "RULE", "arm_field", "centre_by_stage",
           "centre_moves", "cell_factors", "check_line", "depth_departure", "depth_scale", "factor_record",
           "latent_moves", "permute_residuals", "permute_stages", "ratio_log_factors", "recentre",
           "recentre_after_depth", "stage_permutation"]
