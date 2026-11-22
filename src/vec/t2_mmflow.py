"""Row E54 of experiments/harness/TRIES.md (the moonshot "multi-stage flow, then retrieval of real cells", as specified
by a score-blind judge; off by default, used only by experiments/coords/e54_mmflow_read.py and the harness arms
fgw_mmflow_snap / fgw_mmflow_anti).

A cell-state latent and a continuous-time velocity field fitted through every released embryo-series stage at once:

1. Thinning (the embedding only): each cell's panel counts binomially thinned to min(its own panel total, THIN_DEPTH)
   (each count kept with chance min(1, THIN_DEPTH / total)), stream [stage key, THIN_STREAM], the stage key being
   round(100 x the stage's day); features log1p of the thinned counts.
2. PCA: N_PC components of the pooled training stages, PCA_CELLS cells per stage (the same number from each, fewer
   only if a stage holds fewer), stream [stage key, PCA_STREAM]; every cell of every training stage is projected.
3. The velocity net: an MLP of three hidden layers of HIDDEN SiLU units taking (z, t), t the stage's day.
4. Training: multi-marginal conditional flow matching. Each step draws BATCH cells of every training stage, couples
   consecutive stages by minibatch exact OT (squared Euclidean in the latent; with equal uniform marginals the exact
   plan is a permutation, found by the assignment problem), chains the couplings into BATCH paths through all stages,
   and on each path takes one chained triple of consecutive stages (uniformly) and a time uniform in its span; the
   target is the velocity of the natural cubic spline through that triple. STEPS steps of Adam (LR), torch seed
   TORCH_SEED, torch.use_deterministic_algorithms(True); the batch draws come from numpy's default_rng(TORCH_SEED).
5. Integration: dopri5 (torchdiffeq), rtol RTOL, the library's default atol.
6. Retrieval: predicted states against a bracket's whole pool, one to one, by the exact minimum-cost assignment
   (scipy.optimize.linear_sum_assignment, rectangular) on squared Euclidean distance in the latent.

torch and torchdiffeq are imported only by `fit_flow` and `integrate` (the uv gpu group); everything else is numpy.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
from scipy import sparse

THIN_DEPTH = 30  # counts per cell after thinning (a cell with fewer keeps all of its own)
THIN_STREAM = 53  # [stage key, 53]: the thinning
PCA_STREAM = 54  # [stage key, 54]: the PCA's cells per stage
GATE_STREAM = 55  # [draw, 55]: the gate's draws
N_PC = 30
PCA_CELLS = 3000
HIDDEN = 256
N_HIDDEN = 3
STEPS = 20_000
BATCH = 1024
LR = 1e-3
TORCH_SEED = 20261122
RTOL = 1e-5
GATE_DRAWS = 5
GATE_CELLS = 5000
DIRECTIONS = ("snap", "anti")  # the arm (towards the target day) and its placebo (away from it by the same span)
RULE = ("mmflow: panel counts binomially thinned, each count kept with p = min(1, 30 / total) (stream [stage key, "
        "53]); 30 PCs of log1p on the pooled training stages, 3,000 cells per stage (stream [stage key, 54]); MLP 3 x "
        "256 SiLU on (z, t); multi-marginal CFM with minibatch exact OT between consecutive stages and a natural "
        "cubic spline through each chained triple; 20,000 steps, batch 1,024 per stage, Adam 1e-3, torch seed "
        "20261122, deterministic algorithms; dopri5 rtol 1e-5; exact one-to-one retrieval against the bracket's "
        "whole pool")


class MmflowRefused(ValueError):
    """The arm cannot be built as specified (one-to-one retrieval needs at least as many pool cells as states)."""


def stage_key(day: float) -> int:
    """A stage's integer key for its random streams: round(100 x day)."""
    return int(round(100.0 * float(day)))


# ------------------------------------------------------------------------------------------------ embedding
def thin_counts(counts, day: float, depth: int = THIN_DEPTH) -> sparse.csr_matrix:
    """Each cell's counts binomially thinned with p = min(1, depth / its total), stream [stage key, THIN_STREAM].
    Counts are rounded to integers first; the entries are visited in canonical CSR order (sorted, no zeros)."""
    C = sparse.csr_matrix(counts, dtype=np.float64)
    C.data = np.rint(C.data)
    C.eliminate_zeros()
    C.sort_indices()
    if (C.data < 0).any():
        raise ValueError("counts must not be negative")
    total = np.asarray(C.sum(1)).ravel()
    p = np.ones_like(total)
    np.divide(float(depth), total, out=p, where=total > depth)
    rng = np.random.default_rng([stage_key(day), THIN_STREAM])
    data = rng.binomial(C.data.astype(np.int64), np.repeat(p, np.diff(C.indptr))).astype(np.float64)
    out = sparse.csr_matrix((data, C.indices.copy(), C.indptr.copy()), shape=C.shape)
    out.eliminate_zeros()
    return out


def features(thinned) -> np.ndarray:
    """log1p of the thinned counts, dense float64."""
    return np.log1p(np.asarray(sparse.csr_matrix(thinned, dtype=np.float64).toarray()))


def fit_embedding(feats: Mapping[float, np.ndarray], n_pc: int = N_PC,
                  per_stage: int = PCA_CELLS) -> tuple[np.ndarray, np.ndarray]:
    """The PCA of the pooled training stages (day -> cells x genes features): the same number of cells from each,
    min(per_stage, the smallest stage), each stage's drawn on [stage key, PCA_STREAM] without replacement. Returns
    (centre, loadings: genes x components)."""
    m = min(per_stage, *(len(F) for F in feats.values()))
    parts = []
    for day in sorted(feats):
        rows = np.sort(np.random.default_rng([stage_key(day), PCA_STREAM]).choice(len(feats[day]), m, replace=False))
        parts.append(np.asarray(feats[day], dtype=np.float64)[rows])
    X = np.vstack(parts)
    centre = X.mean(0)
    _, _, Vt = np.linalg.svd(X - centre, full_matrices=False)
    return centre, Vt[: min(n_pc, Vt.shape[0])].T


def embed(F: np.ndarray, centre: np.ndarray, loadings: np.ndarray) -> np.ndarray:
    return (np.asarray(F, dtype=np.float64) - centre) @ loadings


# ------------------------------------------------------------------------------------------------ folds
def exclusion(before_day: float, after_day: float, target_day: float, released: Sequence[float]) -> frozenset:
    """The stages a task's fit leaves out: its target and every released stage strictly inside its bracket."""
    return frozenset({float(target_day), *(float(d) for d in released if before_day < d < after_day)})


def neighbours(day: float, released: Sequence[float]) -> tuple[float, float]:
    """The released stages just before and just after `day` (the gate's brackets)."""
    lo = [d for d in released if d < day]
    hi = [d for d in released if d > day]
    if not lo or not hi:
        raise ValueError(f"stage {day} has no released stage on both sides")
    return max(lo), min(hi)


# ------------------------------------------------------------------------------------------------ training batch
def natural_spline3(y: np.ndarray, knots: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Position and velocity at t of the natural cubic spline through three knots, path by path.
    y: (3, B, D) values, knots: (3, B) increasing times, t: (B,) inside [knots[0], knots[2]]."""
    y0, y1, y2 = y
    t0, t1, t2 = knots
    h0, h1 = (t1 - t0)[:, None], (t2 - t1)[:, None]
    if (h0 <= 0).any() or (h1 <= 0).any():
        raise ValueError("knots must increase")
    m1 = 3.0 * ((y2 - y1) / h1 - (y1 - y0) / h0) / (h0 + h1)  # the middle second derivative; zero at both ends
    first = (t <= t1)[:, None]
    ya, yb = np.where(first, y0, y1), np.where(first, y1, y2)
    ma, mb = np.where(first, 0.0, m1), np.where(first, m1, 0.0)
    ta, h = np.where(first[:, 0], t0, t1)[:, None], np.where(first, h0, h1)
    B = (t[:, None] - ta) / h
    A = 1.0 - B
    pos = A * ya + B * yb + ((A ** 3 - A) * ma + (B ** 3 - B) * mb) * h * h / 6.0
    vel = (yb - ya) / h - (3.0 * A * A - 1.0) / 6.0 * h * ma + (3.0 * B * B - 1.0) / 6.0 * h * mb
    return pos, vel


def chain(batches: Sequence[np.ndarray]) -> np.ndarray:
    """Chained minibatch exact OT: paths (K, B, D) through K equal-size batches in stage order, each consecutive pair
    coupled by the exact assignment on squared Euclidean distance (the exact OT plan between equal uniform
    marginals)."""
    from concurrent.futures import ThreadPoolExecutor

    from scipy.optimize import linear_sum_assignment
    from scipy.spatial.distance import cdist

    n = len(batches[0])
    if any(len(b) != n for b in batches):
        raise ValueError("every stage's batch must be of the same size")

    def plan(k: int) -> np.ndarray:  # the pairs are independent; the solver releases the GIL, so threads run them
        r, c = linear_sum_assignment(cdist(batches[k - 1], batches[k], "sqeuclidean"))
        step = np.empty(n, dtype=np.int64)
        step[r] = c
        return step

    with ThreadPoolExecutor(max(1, len(batches) - 1)) as ex:
        steps = list(ex.map(plan, range(1, len(batches))))
    idx = [np.arange(n)]
    for step in steps:
        idx.append(step[idx[-1]])
    return np.stack([np.asarray(b)[i] for b, i in zip(batches, idx, strict=True)])


def training_batch(latents: Sequence[np.ndarray], days: Sequence[float], batch: int,
                   rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One step's (z_t, t, u_t): `batch` cells of each stage (stages in increasing day), chained paths, one chained
    triple per path chosen uniformly, t uniform in its span, the natural spline's position and velocity there."""
    days = np.asarray(days, dtype=np.float64)
    if len(days) < 3 or not np.all(np.diff(days) > 0):
        raise ValueError("training needs at least three stages in increasing day")
    draws = [np.asarray(L, dtype=np.float64)[rng.choice(len(L), batch, replace=len(L) < batch)] for L in latents]
    paths = chain(draws)
    k = rng.integers(0, len(days) - 2, batch)
    lo, hi = days[k], days[k + 2]
    t = lo + rng.random(batch) * (hi - lo)
    cols = np.arange(batch)
    y = np.stack([paths[k + j, cols] for j in range(3)])
    knots = np.stack([days[k + j] for j in range(3)])
    pos, vel = natural_spline3(y, knots, t)
    return pos, t, vel


# ------------------------------------------------------------------------------------------------ torch
def velocity_net(dim: int, hidden: int = HIDDEN, n_hidden: int = N_HIDDEN):
    import torch

    layers, width = [], dim + 1
    for _ in range(n_hidden):
        layers += [torch.nn.Linear(width, hidden), torch.nn.SiLU()]
        width = hidden
    layers.append(torch.nn.Linear(width, dim))
    return torch.nn.Sequential(*layers)


class _Field:
    """The net as an ODE right-hand side f(t, z)."""

    def __init__(self, net):
        self.net = net

    def __call__(self, t, z):
        import torch

        return self.net(torch.cat([z, t.expand(z.shape[0], 1).to(z.dtype)], 1))


def fit_flow(latents: Sequence[np.ndarray], days: Sequence[float], *, steps: int = STEPS, batch: int = BATCH,
             lr: float = LR, seed: int = TORCH_SEED, device: str = "cpu"):
    """The velocity net fitted as in the module docstring; returns (net in eval mode, per-step losses)."""
    import torch

    torch.use_deterministic_algorithms(True)
    torch.manual_seed(seed)
    dim = int(np.asarray(latents[0]).shape[1])
    net = velocity_net(dim).to(device)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    rng = np.random.default_rng(seed)
    order = np.argsort(np.asarray(days, dtype=np.float64))
    lat = [latents[i] for i in order]
    ds = [float(days[i]) for i in order]
    losses = np.empty(steps, dtype=np.float64)
    for s in range(steps):
        z, t, u = training_batch(lat, ds, batch, rng)
        z, t, u = (torch.as_tensor(a, dtype=torch.float32, device=device) for a in (z, t[:, None], u))
        loss = ((net(torch.cat([z, t], 1)) - u) ** 2).mean()
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        losses[s] = float(loss.detach())
    return net.eval(), losses


def integrate(net, z0: np.ndarray, t0: float, t1: float, device: str = "cpu") -> np.ndarray:
    """z0 (cells x components) carried from day t0 to day t1 by the field (dopri5, rtol RTOL), all cells together."""
    import torch
    from torchdiffeq import odeint

    z0 = np.asarray(z0, dtype=np.float32)
    if float(t0) == float(t1) or len(z0) == 0:
        return z0.astype(np.float64).copy()
    with torch.no_grad():
        z = torch.as_tensor(z0, device=device)
        ts = torch.tensor([float(t0), float(t1)], dtype=torch.float32, device=device)
        out = odeint(_Field(net), z, ts, method="dopri5", rtol=RTOL)
    return out[-1].cpu().numpy().astype(np.float64)


# ------------------------------------------------------------------------------------------------ retrieval
def retrieve(pred: np.ndarray, pool: np.ndarray) -> np.ndarray:
    """For each predicted state a distinct pool row, by the exact minimum total squared distance (rectangular
    assignment, predicted states against the whole pool). Refused when the pool holds fewer rows than the states."""
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial.distance import cdist

    pred, pool = np.asarray(pred, dtype=np.float64), np.asarray(pool, dtype=np.float64)
    if len(pred) > len(pool):
        raise MmflowRefused("one-to-one retrieval needs at least as many pool cells as predicted states")
    if len(pred) == 0:
        return np.array([], dtype=np.int64)
    r, c = linear_sum_assignment(cdist(pred, pool, "sqeuclidean"))
    out = np.empty(len(pred), dtype=np.int64)
    out[r] = c
    return out


def kish_ess_frac(rows: np.ndarray, pool_size: int) -> float | None:
    """Kish's effective sample size of the pool's use (each pool row weighted by how often it is used), as a fraction
    of the pool."""
    w = np.bincount(np.asarray(rows, dtype=np.int64), minlength=int(pool_size)).astype(np.float64)
    s2 = float(w @ w)
    return None if s2 <= 0.0 else float(w.sum() ** 2 / s2) / float(pool_size)


def energy_distance(X: np.ndarray, Y: np.ndarray, chunk: int = 512) -> float:
    """The energy distance 2 E|X - Y| - E|X - X'| - E|Y - Y'| (V-statistics over all pairs, Euclidean)."""
    from scipy.spatial.distance import cdist

    def mean_dist(A, B):
        tot = 0.0
        for lo in range(0, len(A), chunk):
            tot += float(cdist(A[lo:lo + chunk], B).sum())
        return tot / (len(A) * len(B))

    X, Y = np.asarray(X, dtype=np.float64), np.asarray(Y, dtype=np.float64)
    return 2.0 * mean_dist(X, Y) - mean_dist(X, X) - mean_dist(Y, Y)


def gate_draws(lo_pool: int, hi_pool: int, truth_pool: int, f: float, s: int,
               n: int = GATE_CELLS) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Draw s of the gate (stream [s, GATE_STREAM]): round(f n) cells of the later neighbour, the rest of the earlier
    one (the calendar fraction f), and n cells of the held stage; without replacement where the pool allows."""
    rng = np.random.default_rng([s, GATE_STREAM])
    n_hi = int(round(f * n))
    n_lo = n - n_hi
    lo = rng.choice(lo_pool, n_lo, replace=n_lo > lo_pool)
    hi = rng.choice(hi_pool, n_hi, replace=n_hi > hi_pool)
    tr = rng.choice(truth_pool, n, replace=n > truth_pool)
    return lo, hi, tr


def gate_stage(lo_lat, lo_fwd, hi_lat, hi_bwd, truth_lat, f: float, draws: int = GATE_DRAWS,
               n: int = GATE_CELLS) -> dict:
    """One held stage of the gate: per draw, the energy distance to the held stage's cells of the flow-predicted states
    (the earlier neighbour's cells carried forward, the later's backward) and of the same cells as measured (the
    calendar-fraction mixture). The flow beats the mixture when its mean over the draws is lower."""
    flow, mix = [], []
    for s in range(draws):
        lo, hi, tr = gate_draws(len(lo_lat), len(hi_lat), len(truth_lat), f, s, n)
        T = np.asarray(truth_lat)[tr]
        flow.append(energy_distance(np.vstack([np.asarray(lo_fwd)[lo], np.asarray(hi_bwd)[hi]]), T))
        mix.append(energy_distance(np.vstack([np.asarray(lo_lat)[lo], np.asarray(hi_lat)[hi]]), T))
    return {"flow": flow, "mixture": mix, "flow_mean": float(np.mean(flow)), "mixture_mean": float(np.mean(mix)),
            "beats": bool(np.mean(flow) < np.mean(mix))}


def gate_pass(stages: Mapping[str, Mapping], must: str, need: int = 4) -> bool:
    """The gate: the flow beats the mixture on at least `need` held stages and on `must`."""
    beats = [bool(v["beats"]) for v in stages.values()]
    return bool(sum(beats) >= need and stages.get(must, {}).get("beats", False))
