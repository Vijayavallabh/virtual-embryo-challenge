# ruff: noqa: E501  (long report lines kept as the decision panel ran this script on 2026-10-21)
"""Monte Carlo for the floor case of approved rule 7: reading G (next best gain) vs reading C (next count up).

Rule 7 (src/vec/t2_build.py COUNT_RULE): on the validation target, 10 paired seeds score shuffle50 predictions at
f = 0.25, 0.5, 0.75, 1 x max_cells. A lower level passes if its paired mean gain over f=1 exceeds 0.3 points and it is
ahead (difference > 0) in at least 8 of 10 seeds; the largest mean gain wins (ties to the higher count). The level
carries over to the test board as the same fraction of the test board's max_cells.

Floor case: the winner is 0.25 and the test board has max_T / truth_T = r_T in [3, 4), so 0.25 x max_T < truth_T is
disallowed while 0.5 and 0.75 x max_T (>= 1.5 x truth_T) are allowed.
  G: among the passing allowed levels, the largest measured mean gain (ties to the higher count); none -> max_cells.
  C: the passing allowed level with the fewest cells; none -> max_cells.
They differ only when 0.5 and 0.75 both pass and gain(0.75) >= gain(0.5): G -> 0.75 x max_T, C -> 0.5 x max_T.

A world is a true score curve g(position) (points) plus geometry. The position at which a level sits:
  mapping "x":   log2(n / truth)        validation f*r_V, test f*r_T   (the scorer's occupancy/neighbourhood effects)
  mapping "f":   log2(n / max_cells)    same on both boards             (stress test requested by the task)
  mapping "abs": log2(n / max_V)        test shifted by Delta = log2(max_T / max_V) (absolute-count stress test)
Per-seed noise of the validation rerun: e[seed, level] ~ N(0, Sigma) with Var(e_i - e_j) = A + B (1 - rho^h),
h = |log2 n_i - log2 n_j| in doublings, (A, B, rho) fitted to the shuffle50 paired-difference SDs of the 20-seed
review (calibrate.py), then scaled so the 1-doubling paired SD equals sd1. Optional multivariate-t tails.
True test score of each choice: g_T(position of the chosen test level), g_T = g + optional linear tilt in position.
D = g_T(C's level) - g_T(G's level) in points; D > 0 means C is better.

Run from the repo root:  uv run --no-sync python -I <this file> [--quick]
Reads only calibration.json in this directory (written by calibrate.py). Writes results.json, results_table.md.
"""

from __future__ import annotations

import json
import math
import sys
import time
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
from scipy.special import ndtr

OUT = Path(__file__).resolve().parent
F = np.array([0.25, 0.5, 0.75, 1.0])
LF = np.log2(F)
MIN_GAIN, MIN_AHEAD, NSEED = 0.3, 8, 10
QUICK = "--quick" in sys.argv

# ------------------------------------------------------------------------------------------------ noise calibration
CAL = json.loads((OUT / "calibration.json").read_text())
S1, S2, S3 = (CAL["shuffle50_sd_by_sep"][k] for k in ("1.0", "2.0", "3.0"))


def fit_noise(s1, s2, s3):
    """Var(h) = A + B (1 - rho^h) through the three (separation, variance) points."""
    v1, v2, v3 = s1 ** 2, s2 ** 2, s3 ** 2
    rho = (v3 - v2) / (v2 - v1)
    B = (v2 - v1) / (rho * (1 - rho))
    A = v1 - B * (1 - rho)
    return A, B, rho


A_FIT, B_FIT, RHO_FIT = fit_noise(S1, S2, S3)
SD1_MEDIAN = S1  # median 1-doubling paired SD over the four review targets (shuffle50, 20 seeds)


def noise_cov(model: str, sd1: float) -> np.ndarray:
    """Covariance of per-level noise e at the four validation levels (log2 n separations = LF differences)."""
    H = np.abs(LF[:, None] - LF[None, :])
    if model == "fit":
        k = (sd1 / S1) ** 2
        A, B, rho = A_FIT * k, B_FIT * k, RHO_FIT
    elif model == "indep":  # equal paired SD at every separation
        A, B, rho = sd1 ** 2, 0.0, 0.5
    elif model == "strong":  # strongly correlated nested draws: rho 0.6 per doubling, small nugget
        rho = 0.6
        A = 0.1 * sd1 ** 2
        B = 0.9 * sd1 ** 2 / (1 - rho)
    else:
        raise ValueError(model)
    return 0.5 * A * np.eye(4) + 0.5 * B * rho ** H


def paired_sd(model: str, sd1: float, h: float) -> float:
    if model == "fit":
        k = (sd1 / S1) ** 2
        return math.sqrt(A_FIT * k + B_FIT * k * (1 - RHO_FIT ** h))
    if model == "indep":
        return sd1
    rho = 0.6
    return math.sqrt(0.1 * sd1 ** 2 + 0.9 * sd1 ** 2 / (1 - rho) * (1 - rho ** h))


# ------------------------------------------------------------------------------------------------ curve families
def U(rng, n, lo, hi):
    return rng.uniform(lo, hi, n)


def logU(rng, n, lo, hi):
    return np.exp(rng.uniform(math.log(lo), math.log(hi), n))


def c_(a):
    return a[:, None]


def cliff(p, pos):
    """Below 1x truth (x mapping only): occupancy_dice thins the truth to the submitted count."""
    return -c_(p["P"]) * np.maximum(0.0, -pos)


def g_inc(p, pos):
    z = pos - c_(p["pos0"])
    return c_(p["a"]) * np.sign(z) * np.abs(z) ** c_(p["c"])


def g_dec(p, pos):
    z = np.maximum(pos - c_(p["pos0"]), 0.0)
    return -c_(p["b"]) * z ** c_(p["q"])


def g_decplat(p, pos):
    z = np.maximum(pos - c_(p["pos0"]), 0.0)
    return -c_(p["A"]) * (1.0 - np.exp(-z / c_(p["w"])))


def g_quad(p, pos):
    z = pos - c_(p["pp"])
    return np.where(z < 0, -c_(p["DL"]) * z ** 2, -c_(p["DR"]) * z ** 2)


def g_bump(p, pos):
    z = pos - c_(p["pp"])
    w = np.where(z < 0, c_(p["wL"]), c_(p["wR"]))
    return c_(p["H"]) * np.exp(-z ** 2 / (2 * w ** 2))


def g_knee(p, pos):
    z = np.maximum(pos, 0.0)
    return -c_(p["K"]) * np.exp(-z / c_(p["w"]))


def g_rugged(p, pos):
    """Independent level effects, defined only at the four level positions (f mapping)."""
    out = np.full(pos.shape, np.nan)
    for j, lf in enumerate(LF):
        out = np.where(np.isclose(pos, lf), c_(p["vals"][:, j]), out)
    if np.isnan(out).any():
        raise ValueError("rugged family needs positions at the levels (f mapping)")
    return out


@dataclass
class Family:
    name: str
    mapping: str  # "x", "f" or "abs"
    shape: str
    desc: str
    ranges: dict = field(default_factory=dict)


def sample_params(fam: Family, rng, n):
    r = fam.ranges
    p = {}
    if fam.shape == "inc":
        p["a"] = U(rng, n, *r.get("a", (0.3, 1.2)))
        p["c"] = U(rng, n, *r.get("c", (0.6, 1.0)))
        p["pos0"] = np.full(n, r.get("pos0", 0.0))
    elif fam.shape == "dec":
        p["b"] = U(rng, n, *r.get("b", (0.2, 1.5)))
        p["q"] = logU(rng, n, *r.get("q", (0.5, 2.0)))
        p["pos0"] = np.full(n, r.get("pos0", 0.0))
    elif fam.shape == "decplat":
        p["A"] = U(rng, n, *r.get("A", (0.3, 2.5)))
        p["w"] = U(rng, n, *r.get("w", (0.2, 1.5)))
        lo, hi = r.get("pos0", (0.0, 0.0))
        p["pos0"] = U(rng, n, lo, hi) if hi > lo else np.full(n, lo)
    elif fam.shape == "quad":
        p["pp"] = U(rng, n, *r.get("pp", (0.0, 2.0)))
        p["DL"] = U(rng, n, *r.get("DL", (0.1, 1.0)))
        p["DR"] = U(rng, n, *r.get("DR", (0.1, 1.0)))
    elif fam.shape == "bump":
        p["pp"] = U(rng, n, *r.get("pp", (0.0, 2.0)))
        p["H"] = U(rng, n, *r.get("H", (0.3, 2.5)))
        p["wL"] = U(rng, n, *r.get("wL", (0.25, 1.5)))
        p["wR"] = U(rng, n, *r.get("wR", (0.25, 1.5)))
    elif fam.shape == "knee":
        p["K"] = U(rng, n, *r.get("K", (0.0, 1.5)))
        p["w"] = U(rng, n, *r.get("w", (0.1, 0.8)))
    elif fam.shape == "gp":
        p["sig"] = U(rng, n, *r.get("sig", (0.2, 1.0)))
        p["ell"] = U(rng, n, *r.get("ell", (0.25, 1.0)))
        blo, bhi = r.get("b", (0.0, 0.0))
        p["b"] = U(rng, n, blo, bhi) if bhi > blo else np.full(n, blo)
        p["pos0"] = np.full(n, r.get("pos0", 0.0))
    elif fam.shape == "rugged":
        mu = U(rng, n, *r.get("mu", (0.0, 1.0)))
        tau = U(rng, n, *r.get("tau", (0.2, 1.5)))
        vals = mu[:, None] + tau[:, None] * rng.standard_normal((n, 4))
        vals[:, 3] = 0.0
        p["vals"] = vals
    else:
        raise ValueError(fam.shape)
    if fam.mapping == "x":
        p["P"] = U(rng, n, *r.get("P", (1.0, 6.0)))
    return p


def gp_curves(fam: Family, p, posV, posT, rng):
    """Smooth trend -b*max(pos - pos0, 0) plus a Gaussian-process wiggle (squared-exponential kernel in log2 units,
    amplitude sig, length ell), drawn jointly at the validation and test positions so that nearby positions share it."""
    P = np.concatenate([posV, posT], axis=1)
    d2 = (P[:, :, None] - P[:, None, :]) ** 2
    s2 = (p["sig"] ** 2)[:, None, None]
    K = s2 * np.exp(-d2 / (2 * (p["ell"] ** 2)[:, None, None])) + 1e-6 * s2 * np.eye(P.shape[1])
    w = (np.linalg.cholesky(K) @ rng.standard_normal((P.shape[0], P.shape[1], 1)))[..., 0]

    def trend(pos):
        g = -c_(p["b"]) * np.maximum(pos - c_(p["pos0"]), 0.0)
        return g + cliff(p, pos) if fam.mapping == "x" else g

    return trend(posV) + w[:, :4], trend(posT) + w[:, 4:]


SHAPES = dict(inc=g_inc, dec=g_dec, decplat=g_decplat, quad=g_quad, bump=g_bump, knee=g_knee, rugged=g_rugged)


def curve(fam: Family, p, pos):
    g = SHAPES[fam.shape](p, pos)
    if fam.mapping == "x":
        g = g + cliff(p, pos)
    return g


# ------------------------------------------------------------------------------------------------ scenario
@dataclass
class Scenario:
    rV: tuple = (4.0, 8.0)  # validation max/truth, uniform
    rT: tuple = (3.0, 4.0)  # test max/truth, uniform on [lo, hi)
    noise: str = "fit"
    sd1: float = SD1_MEDIAN
    tails: float | None = None  # multivariate-t df, None = Gaussian
    tilt_mu: float = 0.0  # test-curve tilt, points per doubling of count (positive favours more cells)
    tilt_sd: float = 0.0
    delta: tuple = (-1.5, 1.5)  # abs mapping: log2(max_T / max_V), uniform
    sd1_range: tuple | None = None  # draw the 1-doubling paired SD per world, log-uniform (overrides sd1)


def positions(fam: Family, sc: Scenario, rng, n):
    if fam.mapping == "x":
        rV = U(rng, n, *sc.rV) if sc.rV[1] > sc.rV[0] else np.full(n, sc.rV[0])
        rT = U(rng, n, *sc.rT) if sc.rT[1] > sc.rT[0] else np.full(n, sc.rT[0])
        posV = np.log2(rV)[:, None] + LF[None, :]
        posT = np.log2(rT)[:, None] + LF[None, :]
        extra = dict(rV=rV, rT=rT)
    elif fam.mapping == "f":
        posV = np.broadcast_to(LF, (n, 4)).copy()
        posT = posV.copy()
        extra = {}
    elif fam.mapping == "abs":
        dlt = U(rng, n, *sc.delta) if sc.delta[1] > sc.delta[0] else np.full(n, sc.delta[0])
        posV = np.broadcast_to(LF, (n, 4)).copy()
        posT = posV + dlt[:, None]
        extra = dict(delta=dlt)
    else:
        raise ValueError(fam.mapping)
    return posV, posT, extra


def simulate_chunk(fam: Family, sc: Scenario, rng, n):
    p = sample_params(fam, rng, n)
    posV, posT, extra = positions(fam, sc, rng, n)
    if fam.shape == "gp":
        gV, gT = gp_curves(fam, p, posV, posT, rng)
    else:
        gV = curve(fam, p, posV)  # (n, 4) true validation scores
        gT = curve(fam, p, posT)  # (n, 4) true test scores at f = 0.25 .. 1
    if sc.tilt_sd > 0 or sc.tilt_mu != 0:
        eta = sc.tilt_mu + sc.tilt_sd * rng.standard_normal(n)
        gT = gT + eta[:, None] * (posT - posT[:, 3:4])
    true_gain = gV[:, :3] - gV[:, 3:4]
    L = np.linalg.cholesky(noise_cov(sc.noise, sc.sd1))
    Z = rng.standard_normal((n, NSEED, 4)) @ L.T
    if sc.sd1_range is not None:
        Z = Z * (logU(rng, n, *sc.sd1_range) / sc.sd1)[:, None, None]
    if sc.tails is not None:
        nu = sc.tails
        Z = Z * np.sqrt((nu - 2) / rng.chisquare(nu, (n, NSEED, 1)))
    d = true_gain[:, None, :] + Z[:, :, :3] - Z[:, :, 3:4]  # (n, seeds, 3) paired differences vs 1max
    gain = d.mean(axis=1)
    ahead = (d > 0).sum(axis=1)
    passing = (gain > MIN_GAIN) & (ahead >= MIN_AHEAD)
    sc_w = np.where(passing, gain, -np.inf)
    idx = 2 - np.argmax(sc_w[:, ::-1], axis=1)  # largest gain, exact ties to the higher count
    winner = np.where(passing.any(axis=1), idx, 3)
    floor = winner == 0
    both = passing[:, 1] & passing[:, 2]
    disagree = floor & both & (gain[:, 2] >= gain[:, 1])
    # choices in the floor case (index into F)
    G = np.where(both, np.where(gain[:, 2] >= gain[:, 1], 2, 1), np.where(passing[:, 1], 1,
                 np.where(passing[:, 2], 2, 3)))
    C = np.where(passing[:, 1], 1, np.where(passing[:, 2], 2, 3))
    D = gT[:, 1] - gT[:, 2]  # C's level minus G's level when they disagree
    rows = np.arange(n)
    out = dict(
        n=n, n_floor=int(floor.sum()), n_dis=int(disagree.sum()),
        n_win=np.bincount(winner, minlength=4),
        n_floor_G=np.bincount(G[floor], minlength=4), n_floor_C=np.bincount(C[floor], minlength=4),
        D=D[disagree],
        # in disagreement: each choice vs keeping max_cells, and vs the best allowed level (0.5, 0.75, 1)
        C_vs_max=(gT[:, 1] - gT[:, 3])[disagree], G_vs_max=(gT[:, 2] - gT[:, 3])[disagree],
        oracle=(gT[:, 1:4].max(axis=1))[disagree] - gT[:, 2][disagree],
        floor_D=np.where(G[floor] == C[floor], 0.0, (gT[rows, C] - gT[rows, G])[floor]),
        true_gap_V=(gV[:, 1] - gV[:, 2])[disagree],  # true validation g(0.5) - g(0.75): is the 0.75>=0.5 read noise?
        meas_gap_V=(gain[:, 1] - gain[:, 2])[disagree],
    )
    for k, v in extra.items():
        out["x_" + k] = v[disagree]
    if fam.shape in ("quad", "bump"):
        out["x_pp"] = p["pp"][disagree]
    return out


def run(fam: Family, sc: Scenario, n_total: int, seed: int, chunk: int = 100_000):
    rng = np.random.default_rng(seed)
    acc = None
    left = n_total
    while left > 0:
        n = min(chunk, left)
        o = simulate_chunk(fam, sc, rng, n)
        if acc is None:
            acc = {k: ([v] if isinstance(v, np.ndarray) and k not in ("n_win", "n_floor_G", "n_floor_C") else v)
                   for k, v in o.items()}
        else:
            for k, v in o.items():
                if k in ("n_win", "n_floor_G", "n_floor_C"):
                    acc[k] = acc[k] + v
                elif isinstance(v, np.ndarray):
                    acc[k].append(v)
                else:
                    acc[k] += v
        left -= n
    for k, v in list(acc.items()):
        if isinstance(v, list):
            acc[k] = np.concatenate(v)
    return acc


def summarise(acc, sc: Scenario):
    n, nf, nd = acc["n"], acc["n_floor"], acc["n_dis"]
    D = acc["D"]
    s = dict(n=n, p_floor=nf / n, p_dis=nd / n, p_dis_given_floor=(nd / nf) if nf else float("nan"), n_dis=nd,
             p_win=(acc["n_win"] / n).tolist())
    h = abs(LF[2] - LF[1])
    sd_real = paired_sd(sc.noise, sc.sd1, h)
    if nd:
        s.update(
            p_C_better=float((D > 1e-9).mean()), p_G_better=float((D < -1e-9).mean()),
            p_tie=float((np.abs(D) <= 1e-9).mean()),
            mean_D=float(D.mean()), sd_D=float(D.std(ddof=1)) if nd > 1 else float("nan"),
            se_mean_D=float(D.std(ddof=1) / math.sqrt(nd)) if nd > 1 else float("nan"),
            q05=float(np.quantile(D, 0.05)), q50=float(np.quantile(D, 0.5)), q95=float(np.quantile(D, 0.95)),
            # one realised test scoring with noise of a single-seed paired difference at 0.585 doublings
            p_C_realised_better=float(np.mean(ndtr(D / sd_real))),
            mean_C_vs_max=float(acc["C_vs_max"].mean()), mean_G_vs_max=float(acc["G_vs_max"].mean()),
            mean_G_regret=float(acc["oracle"].mean()),
            mean_C_regret=float((acc["oracle"] - D).mean()),
            p_trueV_05_ge_075=float((acc["true_gap_V"] > 0).mean()),
            unconditional_gain_C=float(D.sum() / n),
        )
    else:
        s.update(p_C_better=float("nan"), mean_D=float("nan"), unconditional_gain_C=0.0)
    return s


# ------------------------------------------------------------------------------------------------ families
X_FAMILIES = [
    Family("X0 more-is-better (offline)", "x", "inc",
           "g = a*u^c above 1x, a 0.3-1.2 pts/doubling, c 0.6-1 (offline shuffle50: +0.35..1.33/doubling)"),
    Family("X1 decreasing above 1x", "x", "dec", "g = -b*u^q, b 0.2-1.5 pts/doubling, q 0.5-2 (log-uniform)"),
    Family("X1s decrease then plateau", "x", "decplat", "g = -A(1-exp(-u/w)), A 0.3-2.5 pts, w 0.2-1.5 doublings"),
    Family("X2q unimodal, log-quadratic", "x", "quad",
           "peak at 1x-4x, drop 0.1-1 pts at one doubling either side (independent sides)"),
    Family("X2b unimodal bump to plateau", "x", "bump",
           "peak at 1x-4x, height 0.3-2.5 pts, widths 0.25-1.5 doublings per side, flat far away"),
    Family("X3 flat plus penalty below 1x", "x", "knee",
           "g = -K exp(-u/w) above 1x (K 0-1.5, w 0.1-0.8 doublings): flat, knee just above 1x"),
    Family("X4 GP wiggle, flat trend", "x", "gp",
           "flat above 1x plus a smooth random wiggle in log2 x (amplitude 0.2-1 pts, length 0.25-1 doublings)"),
    Family("X4d GP wiggle, decreasing trend", "x", "gp",
           "-b*u (b 0.2-1) plus the same wiggle", dict(b=(0.2, 1.0))),
]
F_FAMILIES = [
    Family("F0 increasing in f", "f", "inc", "g = a*(log2 f + 2)^c", dict(pos0=-2.0)),
    Family("F1 decreasing in f", "f", "dec", "g = -b*(log2 f + 3)^q", dict(pos0=-3.0)),
    Family("F1s decrease-plateau in f", "f", "decplat", "decline starts at f 1/8-0.35", dict(pos0=(-3.0, -1.5))),
    Family("F2q unimodal in f, quad", "f", "quad", "peak at f 0.18-1", dict(pp=(-2.5, 0.0))),
    Family("F2b unimodal bump in f", "f", "bump", "peak at f 0.18-1", dict(pp=(-2.5, 0.0))),
    Family("F3 rugged level effects", "f", "rugged",
           "g(0.25), g(0.5), g(0.75) iid N(mu, tau^2) vs g(1)=0; mu 0-1, tau 0.2-1.5"),
    Family("F4 GP wiggle in f", "f", "gp", "flat plus smooth random wiggle in log2 f, carried over by f"),
]
A_FAMILIES = [
    Family("A1 decreasing in n (abs)", "abs", "dec", "as F1, test shifted by Delta", dict(pos0=-3.0)),
    Family("A2b unimodal bump in n (abs)", "abs", "bump", "as F2b, test shifted by Delta", dict(pp=(-2.5, 0.0))),
    Family("A2q unimodal quad in n (abs)", "abs", "quad", "as F2q, test shifted by Delta", dict(pp=(-2.5, 0.0))),
    Family("A0 increasing in n (abs)", "abs", "inc", "as F0, test shifted by Delta", dict(pos0=-2.0)),
]
ALL = X_FAMILIES + F_FAMILIES + A_FAMILIES


def pooled(rows, weights):
    """Mix families by prior weight; conditioning on disagreement reweights each by P(disagree | family)."""
    num = {k: 0.0 for k in ("p_C_better", "mean_D", "p_G_better", "p_C_realised_better")}
    den = 0.0
    second = 0.0
    for name, w in weights.items():
        s = rows[name]
        if not s["n_dis"]:
            continue
        pw = w * s["p_dis"]
        den += pw
        for k in num:
            num[k] += pw * s[k]
        second += pw * (s["sd_D"] ** 2 + s["mean_D"] ** 2) if s["n_dis"] > 1 else pw * s["mean_D"] ** 2
    if den == 0:
        return None
    out = {k: v / den for k, v in num.items()}
    out["sd_D"] = math.sqrt(max(second / den - out["mean_D"] ** 2, 0.0))
    out["p_dis"] = den / sum(weights.values())
    out["post_weights"] = {name: (w * rows[name]["p_dis"]) / den for name, w in weights.items()
                           if rows[name]["n_dis"]}
    out["unconditional_gain_C"] = sum(w * rows[name]["unconditional_gain_C"] for name, w in weights.items()) / sum(
        weights.values())
    return out


def x_bins(parts):
    """P(C better) and mean D across x families (each case weighted by its family's P(disagree)/n_dis, i.e. an
    equal prior), binned by k = log2(r_V / r_T) and, for the unimodal families, by the peak position in x."""
    D = np.concatenate([p[0] for p in parts])
    k = np.concatenate([p[1] for p in parts])
    pp = np.concatenate([p[2] for p in parts])
    w = np.concatenate([np.full(len(p[0]), p[3]) for p in parts])
    out = {}
    for lo, hi in ((0.0, 0.5), (0.5, 0.9), (0.9, 1.29), (1.29, 1.5)):
        m = (k >= lo) & (k < hi)
        if m.any():
            out[f"k {lo:.2f}-{hi:.2f}"] = dict(share=float(w[m].sum() / w.sum()),
                                               p_C_better=float((w[m] * (D[m] > 1e-9)).sum() / w[m].sum()),
                                               mean_D=float((w[m] * D[m]).sum() / w[m].sum()))
    for lo, hi in ((1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 4.01)):
        m = (2 ** pp >= lo) & (2 ** pp < hi)
        if m.any():
            out[f"peak {lo:g}x-{hi:g}x (unimodal families)"] = dict(
                share_of_unimodal=float(w[m].sum() / w[~np.isnan(pp)].sum()),
                p_C_better=float((w[m] * (D[m] > 1e-9)).sum() / w[m].sum()),
                mean_D=float((w[m] * D[m]).sum() / w[m].sum()))
    return out


def fmt(s, keys):
    return " | ".join(s.get(k, "") if isinstance(s.get(k), str) else
                      (f"{s[k]:.3g}" if isinstance(s.get(k), float) and not math.isnan(s[k]) else
                       (str(s[k]) if k in s else "")) for k in keys)


def main():
    t0 = time.time()
    N_BASE = 200_000 if QUICK else 300_000
    N_RARE = 400_000 if QUICK else 1_200_000  # families where the floor case needs a lot of luck
    N_SENS = 50_000 if QUICK else 60_000
    rare = {"X0 more-is-better (offline)", "X3 flat plus penalty below 1x", "F0 increasing in f",
            "A0 increasing in n (abs)"}
    base = Scenario()
    results = dict(calibration=dict(A=A_FIT, B=B_FIT, rho=RHO_FIT, sd1_median=SD1_MEDIAN,
                                    paired_sd_075_vs_1=paired_sd("fit", SD1_MEDIAN, abs(LF[2])),
                                    paired_sd_05_vs_1=paired_sd("fit", SD1_MEDIAN, 1.0),
                                    paired_sd_025_vs_1=paired_sd("fit", SD1_MEDIAN, 2.0),
                                    paired_sd_05_vs_075=paired_sd("fit", SD1_MEDIAN, abs(LF[2] - LF[1]))),
                   base={}, sens={}, pooled={}, diag={})
    seed = 20261021
    print(f"noise fit: A={A_FIT:.4f} B={B_FIT:.4f} rho={RHO_FIT:.4f}; paired SDs at base sd1={SD1_MEDIAN:.3f}: "
          f"0.25v1 {results['calibration']['paired_sd_025_vs_1']:.3f}, 0.5v1 {results['calibration']['paired_sd_05_vs_1']:.3f}, "
          f"0.75v1 {results['calibration']['paired_sd_075_vs_1']:.3f}, 0.5v0.75 {results['calibration']['paired_sd_05_vs_075']:.3f}")

    # ---- base runs, every family
    xdiag_parts = []
    for i, fam in enumerate(ALL):
        n = N_RARE if fam.name in rare else N_BASE
        acc = run(fam, base, n, seed + i)
        s = summarise(acc, base)
        results["base"][fam.name] = s
        # diagnostics for x families: when does G win?
        if fam.mapping == "x" and acc["n_dis"]:
            k = np.log2(acc["x_rV"] / acc["x_rT"])
            xdiag_parts.append((acc["D"], k, acc.get("x_pp", np.full(acc["n_dis"], np.nan)), s["p_dis"] / acc["n_dis"]))
            Dd = acc["D"]
            diag = dict(p_C_better_k_lt_1_29=float((Dd[k < 1.29] > 1e-9).mean()) if (k < 1.29).any() else None,
                        p_C_better_k_ge_1_29=float((Dd[k >= 1.29] > 1e-9).mean()) if (k >= 1.29).any() else None,
                        share_k_ge_1_29=float((k >= 1.29).mean()))
            if "x_pp" in acc:
                pp = acc["x_pp"]
                diag["peak_x_quartiles_in_disagreement"] = [float(2 ** np.quantile(pp, q)) for q in (0.25, 0.5, 0.75)]
            for lo, hi in ((0.0, 0.5), (0.5, 0.9), (0.9, 1.29), (1.29, 1.5)):
                m = (k >= lo) & (k < hi)
                if m.sum() >= 20:
                    diag[f"k {lo:.2f}-{hi:.2f}"] = dict(n=int(m.sum()), p_C_better=float((Dd[m] > 1e-9).mean()),
                                                        mean_D=float(Dd[m].mean()))
            results["diag"][fam.name] = diag
        if fam.mapping == "abs" and acc["n_dis"]:
            dl = acc["x_delta"]
            Dd = acc["D"]
            results["diag"][fam.name] = {f"p_C_better_delta_{lo:+.1f}_{hi:+.1f}": float((Dd[(dl >= lo) & (dl < hi)] > 1e-9).mean())
                                         for lo, hi in ((-1.5, -0.75), (-0.75, 0.0), (0.0, 0.75), (0.75, 1.5))
                                         if ((dl >= lo) & (dl < hi)).any()}
        print(f"[{time.time() - t0:6.1f}s] {fam.name:34s} n={n:>9,d} P(floor)={s['p_floor']:.4g} "
              f"P(dis)={s['p_dis']:.3g} P(dis|floor)={s['p_dis_given_floor']:.3g} n_dis={s['n_dis']} "
              f"P(C better)={s['p_C_better']:.3f} meanD={s['mean_D']:+.3f}")

    xdiag = x_bins(xdiag_parts)
    names_x = [f.name for f in X_FAMILIES]
    results["diag"]["x families incl. GP, pooled (equal prior), by k = log2(r_V/r_T) and by peak"] = xdiag
    fewer = ["X1 decreasing above 1x", "X1s decrease then plateau", "X2q unimodal, log-quadratic",
             "X2b unimodal bump to plateau"]
    smooth_x = [n for n in names_x if not n.startswith("X4")]
    priors = {
        "x families incl. GP wiggle, equal prior": {n: 1.0 for n in names_x},
        "x smooth families, equal prior": {n: 1.0 for n in smooth_x},
        "x families, offline-heavy prior (X0 0.6, X3 0.1, fewer-is-better 0.075 each)":
            {**{n: 0.075 for n in fewer}, "X0 more-is-better (offline)": 0.6, "X3 flat plus penalty below 1x": 0.1},
        "fewer-is-better x families only": {n: 1.0 for n in fewer},
        "f-mapping stress families, equal prior": {f.name: 1.0 for f in F_FAMILIES},
        "abs-mapping stress families, equal prior": {f.name: 1.0 for f in A_FAMILIES},
        "all families: mapping x 1/2, f 1/4, abs 1/4":
            {**{f.name: 0.5 / len(X_FAMILIES) for f in X_FAMILIES}, **{f.name: 0.25 / len(F_FAMILIES) for f in F_FAMILIES},
             **{f.name: 0.25 / len(A_FAMILIES) for f in A_FAMILIES}},
    }
    for pname, w in priors.items():
        results["pooled"][pname] = pooled(results["base"], w)

    # ---- one-at-a-time sensitivities, x families (pooled under the equal prior) plus the two stress mappings
    sens = {
        "noise sd1 0.45 (large heart truth half)": replace(base, sd1=0.45),
        "noise sd1 1.35 (noisiest target)": replace(base, sd1=1.35),
        "noise independent across levels": replace(base, noise="indep"),
        "noise strongly correlated (rho 0.6/doubling)": replace(base, noise="strong"),
        "noise t-tails (df 4)": replace(base, tails=4.0),
        "noise sd1 per world, log-uniform 0.4-1.35": replace(base, sd1_range=(0.4, 1.35)),
        "r_V in [4, 5]": replace(base, rV=(4.0, 5.0)),
        "r_V in [6, 8]": replace(base, rV=(6.0, 8.0)),
        "r_V = 8 (worst case for C)": replace(base, rV=(8.0, 8.0)),
        "r_V in [3, 4) (0.25x below truth on validation)": replace(base, rV=(3.0, 3.999)),
        "r_T = 3 exactly": replace(base, rT=(3.0, 3.0)),
        "r_T = 3.9": replace(base, rT=(3.9, 3.9)),
        "transfer tilt sd 0.3 pts/doubling": replace(base, tilt_sd=0.3),
        "transfer tilt sd 0.6 pts/doubling": replace(base, tilt_sd=0.6),
        "test favours more cells: tilt mean +0.3, sd 0.3": replace(base, tilt_mu=0.3, tilt_sd=0.3),
        "test favours more cells: tilt mean +0.6, sd 0.3": replace(base, tilt_mu=0.6, tilt_sd=0.3),
    }
    fam_sets = {"x": X_FAMILIES, "f": F_FAMILIES, "abs": A_FAMILIES}
    for j, (sname, sc) in enumerate(sens.items()):
        results["sens"][sname] = {}
        for mname, fams in fam_sets.items():
            rows = {}
            for i, fam in enumerate(fams):
                if mname != "x" and sname.startswith(("r_V", "r_T")):
                    continue  # geometry of truth only matters on the x mapping
                n = N_SENS * (2 if fam.name in rare else 1)
                rows[fam.name] = summarise(run(fam, sc, n, seed + 1000 * (j + 1) + i), sc)
            if not rows:
                continue
            pw = pooled(rows, {k: 1.0 for k in rows})
            results["sens"][sname][mname] = dict(pooled=pw, per_family={k: {kk: v[kk] for kk in (
                "p_floor", "p_dis", "n_dis", "p_C_better", "mean_D")} for k, v in rows.items()})
        xs = results["sens"][sname].get("x", {}).get("pooled")
        print(f"[{time.time() - t0:6.1f}s] sens {sname:52s} x-pooled P(C better)="
              f"{xs['p_C_better'] if xs else float('nan'):.3f} meanD={xs['mean_D'] if xs else float('nan'):+.3f}")

    # ---- parameter-range sensitivities for the fewer-is-better families
    ranges = {
        "X1 shallow (b 0.2-0.6)": Family("X1 decreasing above 1x", "x", "dec", "", dict(b=(0.2, 0.6))),
        "X1 steep (b 1-3)": Family("X1 decreasing above 1x", "x", "dec", "", dict(b=(1.0, 3.0))),
        "X2q peak 1x-2x": Family("X2q", "x", "quad", "", dict(pp=(0.0, 1.0))),
        "X2q peak 2x-4x": Family("X2q", "x", "quad", "", dict(pp=(1.0, 2.0))),
        "X2q peak 2x-4x, shallow right side (DR 0.05-0.3)": Family("X2q", "x", "quad", "", dict(pp=(1.0, 2.0),
                                                                                         DR=(0.05, 0.3))),
        "X2b peak 2x-4x": Family("X2b", "x", "bump", "", dict(pp=(1.0, 2.0))),
        "X2b peak 1.5x-3x, narrow (w 0.25-0.6)": Family("X2b", "x", "bump", "", dict(pp=(0.585, 1.585),
                                                                                    wL=(0.25, 0.6), wR=(0.25, 0.6))),
        "X2b peak 1x-8x": Family("X2b", "x", "bump", "", dict(pp=(0.0, 3.0))),
        "X4 GP, small wiggle (sig 0.05-0.3)": Family("X4", "x", "gp", "", dict(sig=(0.05, 0.3))),
        "X4 GP, large wiggle (sig 0.5-1.5)": Family("X4", "x", "gp", "", dict(sig=(0.5, 1.5))),
        "X4 GP, short length (ell 0.1-0.3)": Family("X4", "x", "gp", "", dict(ell=(0.1, 0.3))),
        "X4 GP, long length (ell 1-2)": Family("X4", "x", "gp", "", dict(ell=(1.0, 2.0))),
        "X4d GP on decreasing trend, small wiggle (sig 0.05-0.3)": Family("X4d", "x", "gp", "",
                                                                          dict(b=(0.2, 1.0), sig=(0.05, 0.3))),
        "F4 GP in f, small wiggle (sig 0.05-0.3)": Family("F4", "f", "gp", "", dict(sig=(0.05, 0.3))),
        "F3 rugged, small level effects (tau 0.1-0.4)": Family("F3", "f", "rugged", "", dict(tau=(0.1, 0.4))),
        "F3 rugged, large level effects (tau 1-3)": Family("F3", "f", "rugged", "", dict(tau=(1.0, 3.0))),
    }
    results["ranges"] = {}
    for j, (rname, fam) in enumerate(ranges.items()):
        s = summarise(run(fam, base, 200_000 if not QUICK else N_SENS * 2, seed + 50_000 + j), base)
        results["ranges"][rname] = {k: s[k] for k in ("p_floor", "p_dis", "n_dis", "p_C_better", "mean_D", "sd_D",
                                                      "q05", "q95") if k in s}
        print(f"[{time.time() - t0:6.1f}s] range {rname:50s} P(dis)={s['p_dis']:.3g} n={s['n_dis']} "
              f"P(C better)={s['p_C_better']:.3f} meanD={s['mean_D']:+.3f}")

    # ---- shift scan: every mapping is a shift s between validation and test level positions on the curve
    # (x: s = log2(r_T/r_V) in [-1.42, 0); f: s = 0; abs: s = log2(max_T/max_V)). Smooth shapes from the abs set.
    results["shift_scan"] = {}
    for j, sh in enumerate((-1.4, -1.0, -0.8, -0.5, 0.0, 0.5, 1.0)):
        sc = replace(base, delta=(sh, sh))
        rows = {}
        for i, fam in enumerate([f for f in A_FAMILIES if f.shape != "inc"]):
            rows[fam.name] = summarise(run(fam, sc, 120_000 if not QUICK else 40_000, seed + 70_000 + 10 * j + i), sc)
        pw = pooled(rows, {k: 1.0 for k in rows})
        results["shift_scan"][f"{sh:+.1f}"] = dict(pooled=pw, per_family={k: {kk: v[kk] for kk in (
            "p_dis", "n_dis", "p_C_better", "mean_D")} for k, v in rows.items()})
        print(f"[{time.time() - t0:6.1f}s] shift {sh:+.1f}: pooled P(C better)={pw['p_C_better']:.3f} "
              f"meanD={pw['mean_D']:+.3f}")

    results["elapsed_s"] = time.time() - t0
    (OUT / "results.json").write_text(json.dumps(results, indent=1, default=float))
    write_table(results)
    print(f"done in {time.time() - t0:.1f}s")


def write_table(R):
    L = []
    w = L.append
    c = R["calibration"]
    w("# Floor-case reading: G (next best gain) vs C (next count up), Monte Carlo results\n")
    w(f"Noise: per-seed paired SD fitted to the 20-seed shuffle50 review (A={c['A']:.3f}, B={c['B']:.3f}, "
      f"rho={c['rho']:.3f}); base 1-doubling paired SD {c['sd1_median']:.2f} points. Paired SDs at the rule's "
      f"levels: 0.25 vs 1 {c['paired_sd_025_vs_1']:.2f}, 0.5 vs 1 {c['paired_sd_05_vs_1']:.2f}, 0.75 vs 1 "
      f"{c['paired_sd_075_vs_1']:.2f}, 0.5 vs 0.75 {c['paired_sd_05_vs_075']:.2f}.\n")
    w("D = true test score of C's choice (0.5 x max_T) minus G's (0.75 x max_T), points, in the disagreement case. "
      "Base geometry r_V = max_V/truth_V ~ U[4, 8], r_T ~ U[3, 4); no transfer tilt.\n")
    px, ps = R["pooled"]["x families incl. GP wiggle, equal prior"], R["pooled"]["x smooth families, equal prior"]
    pa, pf = R["pooled"]["abs-mapping stress families, equal prior"], R["pooled"]["f-mapping stress families, equal prior"]
    xs = [d["x"]["pooled"] for d in R["sens"].values() if d.get("x") and d["x"].get("pooled")]
    w("## Headline\n")
    w(f"- x mapping (score depends on n/truth), all x families: disagreement in {px['p_dis']:.2%} of releases; "
      f"P(C better) {px['p_C_better']:.2f}; mean D {px['mean_D']:+.2f} (sd {px['sd_D']:.2f}). Smooth families only: "
      f"{ps['p_dis']:.2%}, {ps['p_C_better']:.2f}, {ps['mean_D']:+.2f} (sd {ps['sd_D']:.2f}).")
    w(f"- abs mapping: P(C better) {pa['p_C_better']:.2f}, mean D {pa['mean_D']:+.2f}. f-mapping stress test: "
      f"{pf['p_C_better']:.2f}, {pf['mean_D']:+.2f} (G wins only through rugged or wiggly effects carried over by f).")
    w(f"- One-at-a-time sensitivities, x mapping: P(C better) {min(d['p_C_better'] for d in xs):.2f} to "
      f"{max(d['p_C_better'] for d in xs):.2f}, mean D {min(d['mean_D'] for d in xs):+.2f} to "
      f"{max(d['mean_D'] for d in xs):+.2f}. Break-even test tilt toward more cells: {ps['mean_D'] / 0.585:+.2f} "
      f"(smooth) to {px['mean_D'] / 0.585:+.2f} (all x) points per doubling.")
    w(f"- Expected gain of C over G per release (floor geometry assumed): {px['unconditional_gain_C']:+.4f} points "
      f"(all x), {ps['unconditional_gain_C']:+.4f} (smooth x).\n")
    w("## Base, per family\n")
    w("| family | P(floor) | P(disagree) | P(dis given floor) | n disagree | P(C better) | mean D | sd D | 5% | 95% | "
      "P(C better, one realised scoring) | P(true val g0.5 > g0.75) |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for name, s in R["base"].items():
        if s["n_dis"]:
            w(f"| {name} | {s['p_floor']:.3g} | {s['p_dis']:.2g} | {s['p_dis_given_floor']:.2g} | {s['n_dis']} | "
              f"{s['p_C_better']:.3f} | {s['mean_D']:+.3f} | {s['sd_D']:.3f} | {s['q05']:+.3f} | {s['q95']:+.3f} | "
              f"{s['p_C_realised_better']:.3f} | {s['p_trueV_05_ge_075']:.3f} |")
        else:
            w(f"| {name} | {s['p_floor']:.3g} | 0 (of {s['n']:,}) | - | 0 | - | - | - | - | - | - | - |")
    w("\n## Pooled over families (posterior weight = prior x P(disagree | family))\n")
    w("| prior | P(disagree) per release | P(C better) | mean D | sd D | E[gain of C per release] | posterior weights |")
    w("|---|---|---|---|---|---|---|")
    for pname, s in R["pooled"].items():
        if s is None:
            w(f"| {pname} | 0 | - | - | - | - | - |")
            continue
        pw = "; ".join(f"{k.split(' ')[0]} {v:.2f}" for k, v in s["post_weights"].items() if v >= 0.005)
        w(f"| {pname} | {s['p_dis']:.2g} | {s['p_C_better']:.3f} | {s['mean_D']:+.3f} | {s['sd_D']:.3f} | "
          f"{s['unconditional_gain_C']:+.2g} | {pw} |")
    w("\n## One-at-a-time sensitivity (pooled, equal prior within each mapping)\n")
    w("| assumption changed | x: P(dis) | x: P(C better) | x: mean D | f: P(C better) | f: mean D | abs: P(C better) | "
      "abs: mean D |")
    w("|---|---|---|---|---|---|---|---|")

    def cell(d, k):
        if not d or d.get("pooled") is None:
            return "-"
        v = d["pooled"][k]
        return f"{v:+.3f}" if k == "mean_D" else f"{v:.3g}" if k == "p_dis" else f"{v:.3f}"

    bx = R["pooled"]["x families incl. GP wiggle, equal prior"]
    bf = R["pooled"]["f-mapping stress families, equal prior"]
    ba = R["pooled"]["abs-mapping stress families, equal prior"]
    w(f"| (base) | {bx['p_dis']:.3g} | {bx['p_C_better']:.3f} | {bx['mean_D']:+.3f} | {bf['p_C_better']:.3f} | "
      f"{bf['mean_D']:+.3f} | {ba['p_C_better']:.3f} | {ba['mean_D']:+.3f} |")
    for sname, d in R["sens"].items():
        w(f"| {sname} | {cell(d.get('x'), 'p_dis')} | {cell(d.get('x'), 'p_C_better')} | {cell(d.get('x'), 'mean_D')} | "
          f"{cell(d.get('f'), 'p_C_better')} | {cell(d.get('f'), 'mean_D')} | {cell(d.get('abs'), 'p_C_better')} | "
          f"{cell(d.get('abs'), 'mean_D')} |")
    w("\n## Family parameter ranges (base scenario)\n")
    w("| family variant | P(disagree) | n disagree | P(C better) | mean D | 5% | 95% |")
    w("|---|---|---|---|---|---|---|")
    for rname, s in R["ranges"].items():
        if s.get("n_dis"):
            w(f"| {rname} | {s['p_dis']:.2g} | {s['n_dis']} | {s['p_C_better']:.3f} | {s['mean_D']:+.3f} | "
              f"{s['q05']:+.3f} | {s['q95']:+.3f} |")
        else:
            w(f"| {rname} | 0 | 0 | - | - | - | - |")
    w("\n## Shift scan (smooth abs-mapping shapes A1, A2b, A2q, equal prior)\n")
    w("s = test level positions minus validation level positions on the true curve, in doublings of count. "
      "x mapping: s = log2(r_T/r_V), -1.42 to 0; f mapping: s = 0; abs mapping: s = log2(max_T/max_V).\n")
    w("| s | P(disagree) | P(C better) | mean D |")
    w("|---|---|---|---|")
    for sh, d in R.get("shift_scan", {}).items():
        pw = d["pooled"]
        w(f"| {sh} | {pw['p_dis']:.2g} | {pw['p_C_better']:.3f} | {pw['mean_D']:+.3f} |")
    ext = OUT / "shift_scan_extended.txt"
    if ext.exists():
        w("\nExtended to larger negative shifts (test levels far lower on the curve than validation levels; "
          "per-family P(C better)/mean D; A1 at -3.0 is a tie on its flat tail):\n")
        w(ext.read_text().strip())
    w("\n## Diagnostics\n")
    for k, v in R["diag"].items():
        w(f"- {k}: {json.dumps(v)}")
    comp = OUT / "floor_composition.txt"
    if comp.exists():
        w("\n## What each reading picks across all floor cases (floor_composition.py, base, x families)\n")
        rows = comp.read_text().strip().splitlines()
        w("| " + rows[0] + " |")
        w("|" + "---|" * (rows[0].count("|") + 1))
        for r in rows[1:]:
            w("| " + r + " |")
    ver = OUT / "verify_rule.txt"
    if ver.exists():
        w("\n## Rule cross-check (verify_rule.py)\n")
        w(ver.read_text().strip())
    (OUT / "results_table.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
