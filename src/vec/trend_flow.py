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
- Recentring. D(i, g) = r(i, g) / G(lineage, g), G = sum_i c(i, g) r(i, g) / sum_i c(i, g) over the lineage's drawn
  cells (c the drawn cells' counts), so sum_i c(i, g) D(i, g) = sum_i c(i, g). A gene with no count in the lineage keeps
  D = 1, and so does every cell of a lineage b0_trend does not move.

The four arms differ only in where dz comes from:
1. flow: H10's field exactly (`vec.latent_extrap.fold_flow` and `flow_latent`, every FLOW_* constant unchanged).
2. flow_perm, its placebo: the same, with the stage of every cell permuted within its lineage before the pairs are
   formed (`stage_permutation`: each stage keeps its number of cells per lineage; one permutation drawn from
   default_rng([FLOW_SEED, PERM_STREAM])); the drawn reference cells are still the cells that are moved, from their own
   z.
3. flow_pool: every cell's latent centred on the mean of its own lineage at its own stage, u = z - mean z(lineage,
   stage); H10's pairs, weights, network and constants on u in place of z over every transition between consecutive
   input stages; the drawn reference cells' u (their z less their lineage's reference-stage mean) integrated over tau,
   dz that move.
4. flow_pool_perm: 3 on stages permuted within lineage as in 2 (the centring then uses the permuted stages, and the
   drawn cells' u uses the permuted reference stage's lineage mean).
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
RULE = ("b0_trend x D: D(i, g) = r(i, g) / G(lineage, g), r = decoded mean at z + du over decoded mean at z (log "
        "clipped at +-log 3), du = dz - mean dz over the lineage's drawn cells, G the count-weighted mean of r over "
        "them; dz the move over B0's budget from H10's field (flow), its stage-permuted placebo (flow_perm), the field "
        "on latents centred per lineage and stage (flow_pool) or that on permuted stages (flow_pool_perm)")
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
    the reference stage's lineage means the drawn cells are centred on. Arm "flow" is `vec.latent_extrap.fold_flow`."""
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
        if arm.endswith("_perm"):
            zs, labs = permute_stages(zs, labs, [le.FLOW_SEED, PERM_STREAM])
        ref_means = None
        if arm.startswith("flow_pool"):
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


def _cos_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    na, nb = np.linalg.norm(a, axis=1), float(np.linalg.norm(b))
    ok = (na > 0) & (nb > 0)
    out = np.full(len(a), np.nan)
    out[ok] = (a[ok] @ b) / (na[ok] * nb)
    return out


def cell_factors(arm: str, ctx, pool: ad.AnnData, rows: np.ndarray, like: ad.AnnData,
                 damping: float) -> tuple[np.ndarray, dict]:
    """D for b0_trend's drawn cells `like` (typed as `_b0_t2_shift` types them, obs["b0_type"]) of `pool` (rows
    `rows`), and the build's record (shares and ratios only)."""
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
    if logr is None:
        D = np.ones(C.shape)
        logr = np.zeros(C.shape)
    else:
        D = recentre(np.exp(logr), C, types, lineages)
    rec = factor_record(D, logr, C, types, lineages, du, dz)
    rec.update({"arm": arm, "candidate": NAMES[arm], "rule": "vec.trend_flow.RULE", "n_both_zero_means": both_zero,
                "latent_sha256": fold.record.get("latent_sha256"), "params_sha256": frec.get("params_sha256"),
                "tau_days": frec.get("tau_days")})
    return D, {"check": rec, "flow": frec, "fold": fold.record}


def factor_record(D, logr, C, types, lineages, du, dz) -> dict:
    """Shares and ratios only: the lineages moved, the SD of log D over the entries it acts on (a count in a drawn cell
    of a moved lineage), the share of those entries at the clip, the largest relative departure of sum_i c D from
    sum_i c over lineages and genes, the mean cosine between du and the lineage's uncentred mean move."""
    Cd = sparse.csr_matrix(C).toarray().astype(np.float64)
    types = np.asarray(types).astype(str)
    in_moved = np.isin(types, list(lineages))
    acting = (Cd > 0) & in_moved[:, None]
    L = np.log(D[acting]) if acting.any() else np.zeros(0)
    at_clip = np.abs(logr[acting]) >= le.CLIP - 1e-12 if acting.any() else np.zeros(0, dtype=bool)
    dep, cos = 0.0, []
    for lin in lineages:
        sel = types == lin
        if not sel.any():
            continue
        tot = Cd[sel].sum(0)
        has = tot > 0
        if has.any():
            dep = max(dep, float(np.max(np.abs((Cd[sel] * D[sel]).sum(0)[has] - tot[has]) / tot[has])))
        cos.extend(_cos_rows(du[sel], dz[sel].mean(0)).tolist())
    cos = np.asarray(cos, dtype=np.float64)
    return {"lineages_moved": list(lineages), "sd_log_d": float(L.std()) if L.size else 0.0,
            "share_at_clip": float(at_clip.mean()) if at_clip.size else 0.0,
            "max_rel_departure": dep,
            "mean_cosine_du_mean_move": float(np.nanmean(cos)) if np.isfinite(cos).any() else float("nan")}


def check_line(task_key: str, seed: int, rec: dict) -> str:
    c = rec["check"]
    return (f"check {c['candidate']} task={task_key} seed={seed} lineages={','.join(c['lineages_moved'])} "
            f"sd_log_d={c['sd_log_d']:.4g} share_at_clip={c['share_at_clip']:.4g} "
            f"max_rel_departure={c['max_rel_departure']:.3g} mean_cos_du_move={c['mean_cosine_du_mean_move']:.4g} "
            f"latent_sha256={c['latent_sha256']} params_sha256={c['params_sha256']}")


__all__ = ["ARMS", "NAMES", "PERM_STREAM", "RULE", "arm_field", "centre_by_stage", "centre_moves", "cell_factors",
           "check_line", "factor_record", "latent_moves", "permute_stages", "ratio_log_factors", "recentre",
           "stage_permutation"]
