"""Reviewer's large-truth-half variant of review_run.py: score_board(frac=F) with F=0.3 emulates a board whose truth
half is ~8k cells (the heart val_extrap board has truth_cells 8,393 and max_cells 3x that).

Original docstring follows.

Index bookkeeping is re-derived directly with numpy (not via Stage-of-row-indices) and checked against
vec.evaluate.subsample/split_half on the real stage. Conditions follow the RNG recipe described in the experiment so
that its stored rows can be reproduced, plus new conditions:
  shuffle50_<r>x   random subset at r x truth_cells with 50% of the cells' expression rows permuted among themselves
                   (keeps pseudobulk, cell multiset and coordinates; damages only local organisation), i.e. an
                   intermediate-quality prediction for the count question
and a fourth, small-truth target embryo_E7.25 (ref E6.75; truth half ~665 cells).

Usage (DGX, repo root):
  uv run python experiments/density/review/review_big.py --target heart_E9.5 --seed 0 --frac 0.3 \
    --out runs/density_review/raw_big

Copied into the repo from the independent review of 2026-10-19 (originally run from runs/density_review/);
only script paths changed, plus lint-only edits (line wraps; variable O renamed). Outputs stay on the DGX in
runs/density_review/.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from sklearn.neighbors import NearestNeighbors

from vec.boards import RANKED_METRICS, get_board
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, split_half, subsample
from vec.prep import load_stage

TARGETS = {
    "heart_E9.5": ("T2:heart:val_extrap", "data/raw/T2/E9.5.h5ad", "data/raw/T2/E8.75.h5ad"),
    "heart_E8.75": ("T2:heart:val_extrap", "data/raw/T2/E8.75.h5ad", "data/raw/T2/E8.25_late.h5ad"),
    "embryo_E8.0": ("T2:embryo:val_interp", "data/raw/T2/E8.0.h5ad", "data/raw/T2/E7.25.h5ad"),
    "embryo_E7.25": ("T2:embryo:val_interp", "data/raw/T2/E7.25.h5ad", "data/raw/T2/E6.75.h5ad"),
}
RATIOS = [0.5, 1, 2, 3]
JIT = []

ALL_SEEDS = ([f"random_{r}x" for r in (0.5, 1, 2, 3)] + [f"refstage_{r}x" for r in (0.5, 1, 2)]
             + [f"shuffle50_{r}x" for r in (1, 2, 3)])
REPRO_EXTRA = []
EXTRA = ["count_log_ratio", "_dice_voxel_over_nn", "morans_I_agreement"]


FRAC = SUBSAMPLE_FRAC


def wc_indices(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    m = max(2, int(round(n * FRAC)))
    return np.sort(rng.choice(n, m, replace=False))


def half_indices(m: int):
    p = np.random.default_rng(0).permutation(m)
    return np.sort(p[: m // 2]), np.sort(p[m // 2:])


def mnn(c: np.ndarray) -> float:
    d, _ = NearestNeighbors(n_neighbors=2).fit(c).kneighbors(c)
    return float(np.median(d[:, 1]))


def build(pool, ref_pool, tc, seed):
    n_pool = len(pool)
    sizes = {r: min(int(round(r * tc)), n_pool) for r in RATIOS}
    perm = np.random.default_rng([seed, 1]).permutation(n_pool)
    rand = {r: pool[np.sort(perm[:n])] for r, n in sizes.items()}
    n_max = max(sizes.values())
    centre = int(np.random.default_rng([seed, 2]).integers(n_max // 2, n_pool - (n_max - n_max // 2) + 1))
    out = {}
    for i, (r, n) in enumerate(sizes.items()):
        out[f"random_{r}x"] = ("target", rand[r], None, None)
        out[f"block_{r}x"] = ("target", pool[centre - n // 2: centre - n // 2 + n], None, None)
        out[f"head_{r}x"] = ("target", pool[:n], None, None)
        out[f"shuffle50_{r}x"] = ("target", rand[r], None, (0.5, i))
    for k, (r, f, u) in enumerate(JIT):
        out[f"jitter_{r}x_{f}{u}"] = ("target", rand[r], (f, u, k), None)
    rs = {r: min(int(round(r * tc)), len(ref_pool)) for r in RATIOS}
    rperm = np.random.default_rng([seed, 4]).permutation(len(ref_pool))
    nm = max(rs.values())
    rc = int(np.random.default_rng([seed, 5]).integers(nm // 2, len(ref_pool) - (nm - nm // 2) + 1))
    for r, n in rs.items():
        out[f"refstage_{r}x"] = ("ref", ref_pool[np.sort(rperm[:n])], None, None)
        out[f"refblock_{r}x"] = ("ref", ref_pool[rc - n // 2: rc - n // 2 + n], None, None)
        out[f"refhead_{r}x"] = ("ref", ref_pool[:n], None, None)
    return out


def fin(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=sorted(TARGETS))
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", default="runs/density_review/raw_big")
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--frac", type=float, default=0.3)
    args = ap.parse_args()
    global FRAC
    FRAC = args.frac
    bkey, tpath, rpath = TARGETS[args.target]
    board = get_board(bkey)
    t0 = time.time()
    target = Stage.from_anndata(load_stage(tpath, genes=board.genes, lineage=False), board.genes, need_coords=True)
    ref = Stage.from_anndata(load_stage(rpath, genes=board.genes, lineage=False), board.genes, need_coords=True)
    seed = args.seed

    wc = wc_indices(target.n, seed)
    h0, h1 = half_indices(len(wc))
    truth_idx, other_idx = wc[h0], wc[h1]
    W = subsample(target, FRAC, seed)
    T, other = split_half(W, seed=0)
    assert np.array_equal(target.coords[truth_idx], T.coords) and np.array_equal(target.X[truth_idx], T.X)
    assert np.array_equal(target.coords[other_idx], other.coords) and np.array_equal(target.X[other_idx], other.X)
    ref_wc = wc_indices(ref.n, seed + 1)
    assert np.array_equal(ref.X[ref_wc], subsample(ref, FRAC, seed + 1).X)
    pool = np.setdiff1d(np.arange(target.n), wc)
    ref_pool = np.setdiff1d(np.arange(ref.n), ref_wc)
    # content-level leakage check: pool cells whose coordinates coincide exactly with a working-copy cell
    wc_set = {tuple(c) for c in target.coords[wc].tolist()}
    clash = int(sum(tuple(c) in wc_set for c in target.coords[pool].tolist()))
    rwc_set = {tuple(c) for c in ref.coords[ref_wc].tolist()}
    rclash = int(sum(tuple(c) in rwc_set for c in ref.coords[ref_pool].tolist()))
    tc = len(truth_idx)
    nn_truth = mnn(T.coords)
    print(f"[{args.target} s{seed}] target={target.n} ref={ref.n} truth={tc} pool={len(pool)} "
          f"ref_pool={len(ref_pool)} coord_clash={clash} ref_coord_clash={rclash} ({time.time() - t0:.1f}s)",
          flush=True)

    want = ALL_SEEDS + (REPRO_EXTRA if seed < 5 else [])
    if args.only:
        want = [w for w in want if w in args.only]
    conds = build(pool, ref_pool, tc, seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{args.target}_s{seed}.jsonl"
    done = set()
    if path.exists():
        done = {json.loads(x)["condition"] for x in path.read_text().splitlines() if x.strip()}
    for label in want:
        if label in done:
            continue
        source, idx, jit, shuf = conds[label]
        assert len(np.intersect1d(idx, wc if source == "target" else ref_wc)) == 0
        pred = (target if source == "target" else ref).take(idx)
        nn_own = mnn(pred.coords)
        jr = 0.0
        if jit is not None:
            f, u, k = jit
            jr = f * (nn_own if u == "own" else nn_truth)
            noise = np.random.default_rng([seed, 3, k]).normal(size=pred.coords.shape) * (jr / np.sqrt(3))
            pred = Stage(X=pred.X, coords=(pred.coords + noise).astype(np.float32), celltype=pred.celltype)
        if shuf is not None:
            p, ri = shuf
            rng = np.random.default_rng([seed, 7, ri])
            n = pred.n
            sel = rng.choice(n, int(round(p * n)), replace=False)
            X = pred.X.copy()
            X[sel] = X[sel[rng.permutation(len(sel))]]
            pred = Stage(X=X, coords=pred.coords, celltype=None)
        t1 = time.time()
        res = score_board(board, pred, target, ref, seed=seed, anchors="computed", frac=FRAC)
        row = dict(target=args.target, frac=FRAC, seed=seed, condition=label, source=source, n_pred=pred.n,
                   truth_cells=tc,
                   pool=len(pool), ref_pool=len(ref_pool), coord_clash=clash, nn_own=nn_own, nn_truth=nn_truth,
                   jitter_rms=jr, score=res.score, group_skill=res.group_skill,
                   skills={m: fin(res.skills.get(m)) for m in RANKED_METRICS["T2"]},
                   metrics={m: fin(res.metrics.get(m)) for m in RANKED_METRICS["T2"] + EXTRA},
                   floor={m: fin(res.anchors.get(m, {}).get("floor")) for m in RANKED_METRICS["T2"]},
                   ceiling={m: fin(res.anchors.get(m, {}).get("ceiling")) for m in RANKED_METRICS["T2"]},
                   inverted=res.extras["inverted_anchors"], seconds=time.time() - t1)
        with path.open("a") as fh:
            fh.write(json.dumps(row) + "\n")
        print(f"[{args.target} s{seed}] {label:22s} n={pred.n:6d} score={res.score:7.3f} ({row['seconds']:.1f}s)",
              flush=True)
    print(f"[{args.target} s{seed}] done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
