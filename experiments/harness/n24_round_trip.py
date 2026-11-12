#!/usr/bin/env python
"""Diagnostic N24 of experiments/harness/TRIES.md (target-free, not a try): how faithful is the depth top-up of the
renderings (tries E19c, E24, E25), read on a round trip inside one deep released stage?

The second research pass (experiments/research_2026_11_runs/REPORT2.md, recipe 4) proposed it: thin a deep stage's
cells to the shallow bracket's depth, render them back up, and compare the result with held-out deep cells of the
same stage by the scorer's own distances. Nothing but that one released stage is read: no practice target, no test
file. It bounds what a sharper top-up could recover and can choose its settings without any target.

Per stage (E8.0, E8.25_late, E8.5): the cells are split once into halves A and B. A's panel counts are thinned
binomially to 1 / RATIO of their depth (RATIO = 3.42, the test pair's depth ratio as diagnostic N23 read it). Each
variant renders the thinned cells back to their original depth; the scorer's mmd_unbiased and variogram_score
(common.core_metrics, as the emulator calls them) are then taken against B on log(1 + 1e4 c / N), over five scorer
seeds.

Variants, fixed before the run:
  thin          the thinned cells as they are (the floor of this trip)
  true          A's own deep counts (the ceiling: two halves of one sample)
  poisson_k15   today's top-up: Poisson draws from the mean proportions of the cell's 15 nearest thinned cells (30 PCs)
  poisson_k30, poisson_k50    the same with larger neighbourhoods
  post_k15, post_k30, post_k50    a posterior top-up: per gene a Gamma prior with the neighbourhood rate as mean and a
                dispersion measured on deep cells (half B, method of moments around their own kNN means), updated by
                the cell's thinned counts; the top-up is a Poisson draw from a sample of that posterior
  self          a Poisson top-up from the cell's own thinned proportions (no neighbours; a reference, not a candidate)

Reading, fixed before the run. A variant is taken to a confirmation try only if on all three stages it lowers BOTH
distances to B below poisson_k15's by more than the seed spread (the SD of poisson_k15's value over the five scorer
seeds); among those the one with the lowest mean ratio to poisson_k15 over stages and distances. If none qualifies the
top-up stays as built. Reported either way: how much of the gap between thin and true each variant closes
(1 - (d - d_true) / (d_thin - d_true)), the ceiling for any sharper top-up.

Usage (repo root):
  uv run --frozen --no-sync python -I -B experiments/harness/n24_round_trip.py \
    --out experiments/harness/results/n24_round_trip
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from vec import harness, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.prep import load_stage

DATA = REPO_ROOT / "data" / "raw" / "T2"
BOARD = "T2:embryo:test_interp"
STAGES = ("E8.0", "E8.25_late", "E8.5")
RATIO = 3.42
SPLIT_SEED = 20261116
THIN_SEED = 20261117
RENDER_SEED = 20261118
N_RENDER = 6000  # cells of half A rendered and scored
N_TRUTH = 6000  # cells of half B scored against
SCORER_SEEDS = (0, 1, 2, 3, 4)
KS = (15, 30, 50)
N_PC = t2_build.RENDER_N_PC
ALPHA_CAP = 1e6
VARIANTS = ("thin", "true", *[f"poisson_k{k}" for k in KS], *[f"post_k{k}" for k in KS], "self")
BASE = "poisson_k15"


def lognorm(C) -> np.ndarray:
    C = sparse.csr_matrix(C, dtype=np.float64)
    tot = np.asarray(C.sum(1)).ravel()
    X = sparse.diags(1e4 / np.maximum(tot, 1.0)) @ C
    X.data = np.log1p(X.data)
    return np.asarray(X.toarray(), dtype=np.float32)


def knn_rates(C: sparse.csr_matrix, rows: np.ndarray, k: int) -> np.ndarray:
    """Mean over each row's k nearest cells of the pool (itself included; Euclidean in the first N_PC principal
    components of the pool's log-normalised counts) of counts / total: t2_build._smoothed_props with k free."""
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    E = lognorm(C)
    Z = PCA(n_components=int(min(N_PC, E.shape[0] - 1, E.shape[1])), svd_solver="full").fit_transform(
        np.asarray(E, dtype=np.float64))
    nn = NearestNeighbors(n_neighbors=k).fit(Z).kneighbors(Z[rows], return_distance=False)
    tot = np.asarray(C.sum(1)).ravel().astype(np.float64)
    P = sparse.diags(1.0 / np.maximum(tot, 1.0)) @ sparse.csr_matrix(C, dtype=np.float64)
    A = sparse.csr_matrix((np.full(nn.size, 1.0 / k), nn.ravel(), np.arange(0, nn.size + 1, k)),
                          shape=(len(rows), C.shape[0]))
    return np.asarray((A @ P).toarray(), dtype=np.float64)


def dispersions(C_deep: sparse.csr_matrix, k: int, rng) -> np.ndarray:
    """Per-gene overdispersion of deep counts around their kNN means: phi = sum((c - mu)^2 - mu) / sum(mu^2), floored
    at 0, on a sample of the deep cells."""
    rows = rng.choice(C_deep.shape[0], min(4000, C_deep.shape[0]), replace=False)
    S = knn_rates(C_deep, rows, k)
    c = np.asarray(C_deep[rows].toarray(), dtype=np.float64)
    mu = S * c.sum(1, keepdims=True)
    num = ((c - mu) ** 2 - mu).sum(0)
    den = (mu ** 2).sum(0)
    return np.where(den > 0, np.maximum(num / np.maximum(den, 1e-300), 0.0), 0.0)


def render(variant: str, C_thin: sparse.csr_matrix, rows: np.ndarray, phi_by_k: dict, rng) -> np.ndarray:
    c = np.asarray(C_thin[rows].toarray(), dtype=np.float64)
    lib = c.sum(1, keepdims=True)
    extra = (RATIO - 1.0) * lib
    if variant == "self":
        return c + rng.poisson(c / np.maximum(lib, 1.0) * extra)
    kind, k = variant.split("_k")
    k = int(k)
    S = knn_rates(C_thin, rows, k)
    if kind == "poisson":
        return c + rng.poisson(S * extra)
    alpha = np.minimum(1.0 / np.maximum(phi_by_k[k], 1.0 / ALPHA_CAP), ALPHA_CAP)[None, :]
    rate = alpha / np.maximum(S, 1e-300) + lib  # Gamma prior with mean S, updated by the cell's own counts
    lam = rng.gamma(alpha + c, 1.0) / rate
    lam = np.where(S > 0, lam, 0.0)
    return c + rng.poisson(lam * extra)


def main(argv=None) -> None:
    from common.core_metrics import mmd_unbiased, variogram_score

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    commit = harness._git("rev-parse", "HEAD")
    clean = harness._git("status", "--porcelain", "--untracked-files=no") == ""
    genes = list(load_boards()[BOARD].genes)
    res = {"what": __doc__.split("\n")[0], "commit": commit, "tree_clean": clean, "ratio": RATIO, "stages": {}}
    lines = [f"N24 (commit {commit}, clean {clean}; target-free): the top-up's round trip, thinning to 1/{RATIO}"]
    for si, stage in enumerate(STAGES):
        st = load_stage(DATA / f"{stage}.h5ad", genes=genes, drop_duplicates=True)
        _, C, _, _ = t2_build._stage_counts(st, genes)
        perm = np.random.default_rng([SPLIT_SEED, si]).permutation(C.shape[0])
        ia, ib = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2:])
        Ca, Cb = sparse.csr_matrix(C[ia]), sparse.csr_matrix(C[ib])
        trng = np.random.default_rng([THIN_SEED, si])
        Ct = sparse.csr_matrix((trng.binomial(Ca.data, 1.0 / RATIO), Ca.indices.copy(), Ca.indptr.copy()),
                               shape=Ca.shape)
        Ct.eliminate_zeros()
        ok = np.flatnonzero(np.asarray(Ct.sum(1)).ravel() > 0)
        rrng = np.random.default_rng([RENDER_SEED, si])
        rows = np.sort(rrng.choice(ok, min(N_RENDER, len(ok)), replace=False))
        truth = lognorm(Cb[np.sort(rrng.choice(Cb.shape[0], min(N_TRUTH, Cb.shape[0]), replace=False))])
        phi = {k: dispersions(Cb, k, np.random.default_rng([RENDER_SEED, si, k])) for k in KS}
        out = {}
        for vi, v in enumerate(VARIANTS):
            vr = np.random.default_rng([RENDER_SEED, si, 100 + vi])
            if v == "thin":
                counts = np.asarray(Ct[rows].toarray(), dtype=np.float64)
            elif v == "true":
                counts = np.asarray(Ca[rows].toarray(), dtype=np.float64)
            else:
                counts = render(v, Ct, rows, phi, vr)
            X = lognorm(counts)
            m = [float(mmd_unbiased(X, truth, seed=s)) for s in SCORER_SEEDS]
            g = [float(variogram_score(X, truth, seed=s)) for s in SCORER_SEEDS]
            out[v] = {"mmd_u": float(np.mean(m)), "mmd_u_sd": float(np.std(m, ddof=1)),
                      "variogram": float(np.mean(g)), "variogram_sd": float(np.std(g, ddof=1)),
                      "median_counts": float(np.median(counts.sum(1)))}
        for v in VARIANTS:
            for met in ("mmd_u", "variogram"):
                gap = out["thin"][met] - out["true"][met]
                out[v][f"{met}_gap_closed"] = float(1.0 - (out[v][met] - out["true"][met]) / gap) if gap > 0 else None
        res["stages"][stage] = {"cells": int(C.shape[0]), "median_counts_deep": float(np.median(np.asarray(
            Ca.sum(1)).ravel())), "median_phi": {str(k): float(np.median(phi[k])) for k in KS}, "variants": out}
        deep = res["stages"][stage]["median_counts_deep"]
        lines.append(f"{stage} ({C.shape[0]} cells; median panel counts {deep:.0f}; median dispersion at k = 15: "
                     f"{np.median(phi[15]):.3f}):")
        for v in VARIANTS:
            o = out[v]
            lines.append(f"  {v:<12} mmd_u {o['mmd_u']:.5f} (sd {o['mmd_u_sd']:.5f}; gap closed "
                         f"{o['mmd_u_gap_closed']:.3f})   variogram {o['variogram']:.5f} (sd {o['variogram_sd']:.5f}; "
                         f"gap closed {o['variogram_gap_closed']:.3f})   median counts {o['median_counts']:.0f}")
    cands = [v for v in VARIANTS if v not in ("thin", "true", "self", BASE)]
    qualifies = {}
    for v in cands:
        ok_all, ratios = True, []
        for stage in STAGES:
            o, b = res["stages"][stage]["variants"][v], res["stages"][stage]["variants"][BASE]
            for met in ("mmd_u", "variogram"):
                ok_all &= bool(b[met] - o[met] > b[f"{met}_sd"])
                ratios.append(o[met] / b[met] if b[met] > 0 else float("nan"))
        qualifies[v] = {"passes": ok_all, "mean_ratio_to_base": float(np.mean(ratios))}
    winners = [v for v in cands if qualifies[v]["passes"]]
    pick = min(winners, key=lambda v: qualifies[v]["mean_ratio_to_base"]) if winners else None
    res["reading"] = {"qualifies": qualifies, "pick": pick}
    lines.append("reading: " + (f"{pick} goes to a confirmation try" if pick else "no variant qualifies; the top-up "
                                "stays as built") + "; mean ratio to poisson_k15: "
                 + ", ".join(f"{v} {q['mean_ratio_to_base']:.3f}{' (passes)' if q['passes'] else ''}"
                             for v, q in qualifies.items()))
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "round_trip.json").write_text(json.dumps(res, indent=1))
    (out_dir / "round_trip.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
