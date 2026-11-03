#!/usr/bin/env python
"""Embryo coupling diagnostic (2026-11-03; target-free: it reads no target stage and no score). The gate of embryo try
E10 in experiments/harness/TRIES.md, fixed there before this ran.

For each bracket pair of the embryo board (the two backtests' and the test's) and each build seed, it measures, for
the entry's moves (the FGW coupling) and for try E10's (vec.core.aligned_coupling in the rotated common frame):
- whether the anchors' partners agree along every direction: the smallest eigenvalue of the symmetrised whitened
  cross-correlation of anchor and partner positions in the before stage's frame (about 1 when partners sit where
  their anchors do, about 0 along a direction where they sit at random, below 0 where they are mirrored);
- whether the moved cloud keeps its brackets' thickness: the thinnest principal SD over the RMS radius of the moved
  pool at the pair's fraction, against the thinner bracket's;
- the partners' agreement on the lineage map's merged labels, the rotation, the coupling's support, the anchors' mean
  shift, the shape at t_pos = 0 (the s2's geometry) and how many of each cell's 15 nearest neighbours stay among them.

The interpolation is the entry's own function (vec.t2_build.interpolate_clouds) on the stages loaded as the harness
loads them, so the couplings are the ones the harness and the s1 build use at these seeds. The test board is not in
boards/ until the board intake; every pair uses T2:embryo:val_interp's panel, the backtests' own, and the report says
so. The test pair's verdict at the s1 seed (1,000,001) is reported for the captain whatever the gate says.

Usage (DGX, repo root, BLAS threads capped as docs/COMPUTE.md says):
  uv run python experiments/harness/diag_emb_coupling.py --out experiments/harness/results/diag_emb_coupling
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from vec.boards import REPO_ROOT, load_boards
from vec.harness import BUILD_SEED_OFFSET
from vec.io import COORD_KEY, dense
from vec.prep import load_stage
from vec.t2_build import MOVE_COUPLINGS, interpolate_clouds

DATA = REPO_ROOT / "data" / "raw" / "T2"
GENES_BOARD = "T2:embryo:val_interp"
# (key, before, after, fraction): the two counted embryo backtests and the test pair
PAIRS = (("T2e/E6.75+E7.5->E7.25", "E6.75", "E7.5", 2 / 3),
         ("T2e/E7.25+E8.0->E7.5", "E7.25", "E8.0", 1 / 3),
         ("test/E7.5+E8.0->E7.75", "E7.5", "E8.0", 1 / 2))
TEST_PAIR = PAIRS[2][0]
SEEDS = tuple(range(10))  # build seeds BUILD_SEED_OFFSET + s, as the harness's; s = 1 is the s1 build's seed
S1_SEED = 1
N_DRAW = 5000  # the draw's size changes neither the coupling nor the moved pool, which is what this measures
NN_K = 15  # neighbours per cell, as the scorer's neighbourhoods
NN_SAMPLE = 2000  # cells per stage in the neighbour check
NN_STREAM = 11  # its random stream, apart from the interpolation's 0-7
TRACKED_OUT = ("Untracked", "Unknown")

# The gate, fixed in TRIES.md (row E10) before this ran; decided by the AI assistant under the captain's delegation of
# 2026-11-02 (the ideas workflow's proposal, unchanged).
DEFECT_XCORR = 0.3  # the entry's partners: smallest whitened cross-correlation below this is a defect ...
DEFECT_THIN = 0.9  # ... or a moved cloud thinner than this fraction of the thinner bracket
FIX_XCORR = 0.6  # the fix works: the aligned partners' smallest cross-correlation at least this ...
FIX_LABEL_DROP = 0.05  # ... and their label agreement at most this far below the entry's
BACKTEST_MIN = 6  # draws (of 10) with the defect and a working fix, on each backtest pair
TEST_MIN = 3  # ... on the test pair (or the s1 seed's draw)


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True, cwd=REPO_ROOT).stdout.strip()


def _shape(C: np.ndarray, w: np.ndarray | None = None) -> list[float]:
    """Principal SDs over the RMS radius, largest first (weights optional)."""
    w = np.full(len(C), 1.0 / len(C)) if w is None else np.asarray(w, float) / np.sum(w)
    X = C - w @ C
    ev = np.clip(np.linalg.eigvalsh((X * w[:, None]).T @ X), 0.0, None)[::-1]
    return [float(v) for v in np.sqrt(ev / ev.sum())]


def _xcorr(X: np.ndarray, Y: np.ndarray) -> list[float]:
    """Eigenvalues (smallest first) of the symmetrised whitened cross-correlation of paired rows of X and Y."""
    X, Y = X - X.mean(0), Y - Y.mean(0)

    def isqrt(S):
        ev, V = np.linalg.eigh(S)
        return V @ np.diag(ev ** -0.5) @ V.T

    K = isqrt(X.T @ X / len(X)) @ (X.T @ Y / len(X)) @ isqrt(Y.T @ Y / len(Y))
    return [float(v) for v in np.linalg.eigvalsh((K + K.T) / 2)]


def _nn_kept(before: np.ndarray, after: np.ndarray, rows: np.ndarray) -> float:
    """Mean share of each sampled cell's NN_K nearest sampled neighbours that are still among them after the move."""
    a, b = before[rows], after[rows]
    ia = cKDTree(a).query(a, k=NN_K + 1)[1][:, 1:]
    ib = cKDTree(b).query(b, k=NN_K + 1)[1][:, 1:]
    return float(np.mean([len(np.intersect1d(x, y)) / NN_K for x, y in zip(ia, ib, strict=True)]))


def _measure(res, t: float, lab_b: np.ndarray, lab_a: np.ndarray, nn_rows: tuple[np.ndarray, np.ndarray]) -> dict:
    c = res.coupling
    ib, ia, Tm = c["anchors_before"], c["anchors_after"], c["moves"]
    Zb, Za = c["before_frame"], c["after_frame"]
    partners = (Tm @ Za[ia]) / Tm.sum(1, keepdims=True)
    xc = _xcorr(Zb[ib], partners)
    lb, la = lab_b[ib], lab_a[ia]
    tracked = ~np.isin(lb, TRACKED_OUT)[:, None] & ~np.isin(la, TRACKED_OUT)[None, :]
    same = (lb[:, None] == la[None, :]) & tracked
    mass = float((Tm * tracked).sum())
    mb, ma = c["shift_before"], c["shift_after"]
    w = np.r_[np.full(len(Zb), (1.0 - t) / len(Zb)), np.full(len(Za), t / len(Za))]
    at_t = np.vstack([Zb + t * mb, Za + (1.0 - t) * ma])
    at_0 = np.vstack([Zb, Za + ma])
    sb, sa, st = _shape(Zb), _shape(Za), _shape(at_t, w)
    kept = (1.0 - t) * _nn_kept(Zb, Zb + t * mb, nn_rows[0]) + t * _nn_kept(Za, Za + (1.0 - t) * ma, nn_rows[1])
    label = float((Tm * same).sum() / mass) if mass > 0 else None
    return {"xcorr_eigs": xc, "xcorr_min": xc[0], "label_agreement": label,
            "coupling_nonzero": int((Tm > 0).sum()), "rotation_deg": res.record["rotation_deg"],
            "mean_anchor_shift": res.record["mean_anchor_shift"],
            "shape_before": sb, "shape_after": sa, "shape_at_t": st, "shape_at_t_pos0": _shape(at_0, w),
            "shape_mapped_after_half_t_pos0": _shape(Za + ma),
            "thin_ratio": st[-1] / min(sb[-1], sa[-1]), "nn15_kept": kept}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, nargs="*", default=list(SEEDS))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    commit, clean = _git("rev-parse", "HEAD"), _git("status", "--porcelain", "--untracked-files=no") == ""
    genes = list(load_boards()[GENES_BOARD].genes)
    stages = sorted({s for _, b, af, _ in PAIRS for s in (b, af)})
    loaded = {s: load_stage(DATA / f"{s}.h5ad", genes=genes, drop_duplicates=True) for s in stages}
    report = {"what": __doc__.split("\n")[0], "commit": commit, "clean": clean, "genes_board": GENES_BOARD,
              "note": "the test board is not in boards/ before the board intake; every pair uses the backtests' "
                      "panel (T2:embryo:val_interp)",
              "gate": {"DEFECT_XCORR": DEFECT_XCORR, "DEFECT_THIN": DEFECT_THIN, "FIX_XCORR": FIX_XCORR,
                       "FIX_LABEL_DROP": FIX_LABEL_DROP, "BACKTEST_MIN": BACKTEST_MIN, "TEST_MIN": TEST_MIN,
                       "S1_SEED": S1_SEED},
              "pairs": {}}
    for key, sb_, sa_, t in PAIRS:
        before, after = loaded[sb_], loaded[sa_]
        Cb = np.asarray(before.obsm[COORD_KEY], float)[:, :3]
        Ca = np.asarray(after.obsm[COORD_KEY], float)[:, :3]
        Eb, Ea = dense(before[:, genes].X).astype(np.float32), dense(after[:, genes].X).astype(np.float32)
        lab_b, lab_a = before.obs["merged"].astype(str).to_numpy(), after.obs["merged"].astype(str).to_numpy()
        draws = []
        for s in a.seeds:
            seed = BUILD_SEED_OFFSET + s
            t0 = time.time()
            rng = np.random.default_rng([seed, NN_STREAM])
            nn_rows = (rng.choice(len(Cb), min(NN_SAMPLE, len(Cb)), replace=False),
                       rng.choice(len(Ca), min(NN_SAMPLE, len(Ca)), replace=False))
            row = {"seed": s, "build_seed": seed}
            fgw_T = None
            for mc in MOVE_COUPLINGS:
                res = interpolate_clouds(Cb, Ca, Eb, Ea, t, N_DRAW, seed, target_rms=1.0, move_coupling=mc,
                                         with_coupling=True)
                if fgw_T is None:
                    fgw_T = res.coupling["fgw"]
                elif not np.array_equal(fgw_T, res.coupling["fgw"]):
                    raise RuntimeError("the FGW coupling is not reproducible within a seed")
                row[mc] = _measure(res, t, lab_b, lab_a, nn_rows)
            f, al = row["fgw"], row["aligned"]
            row["defect"] = bool(f["xcorr_min"] < DEFECT_XCORR or f["thin_ratio"] < DEFECT_THIN)
            row["fixed"] = bool(al["xcorr_min"] >= FIX_XCORR and (f["label_agreement"] is None or (
                al["label_agreement"] is not None and al["label_agreement"] >= f["label_agreement"] - FIX_LABEL_DROP)))
            row["counted"] = row["defect"] and row["fixed"]
            row["seconds"] = round(time.time() - t0, 1)
            draws.append(row)
            print(f"{key} seed {s}: fgw xcorr {f['xcorr_min']:.2f} thin {f['thin_ratio']:.2f} label "
                  f"{f['label_agreement']} | aligned xcorr {al['xcorr_min']:.2f} thin {al['thin_ratio']:.2f} "
                  f"label {al['label_agreement']} | defect {row['defect']} fixed {row['fixed']} "
                  f"({row['seconds']} s)", flush=True)
        n_counted = sum(r["counted"] for r in draws)
        s1 = next((r for r in draws if r["seed"] == S1_SEED), None)
        need = TEST_MIN if key == TEST_PAIR else BACKTEST_MIN
        passes = n_counted >= need or (key == TEST_PAIR and s1 is not None and s1["counted"])
        report["pairs"][key] = {"before": sb_, "after": sa_, "t": t, "cells": [int(before.n_obs), int(after.n_obs)],
                                "n_defect": sum(r["defect"] for r in draws), "n_counted": n_counted, "need": need,
                                "passes": bool(passes), "draws": draws}
    p = report["pairs"]
    backtests = [k for k, *_ in PAIRS if k != TEST_PAIR]
    if len(a.seeds) != len(SEEDS):
        verdict = "partial run: no verdict"
    elif not p[TEST_PAIR]["passes"]:
        verdict = "prune E10: no defect the fix removes at test"
    elif all(p[k]["passes"] for k in backtests):
        verdict = "run E10: the defect shows and the fix works on the test pair and both backtest pairs"
    else:
        verdict = "prune E10 and report a test-only risk: the backtests cannot judge the fix"
    report["verdict"] = verdict
    s1 = next((r for r in p[TEST_PAIR]["draws"] if r["seed"] == S1_SEED), None)
    if s1 is not None:
        report["s1_draw"] = {"defect": s1["defect"], "fgw_thin_ratio": s1["fgw"]["thin_ratio"],
                             "fgw_xcorr_min": s1["fgw"]["xcorr_min"], "aligned_thin_ratio": s1["aligned"]["thin_ratio"],
                             "aligned_xcorr_min": s1["aligned"]["xcorr_min"]}
    (out / "diag_emb_coupling.json").write_text(json.dumps(report, indent=1) + "\n")
    print(verdict)


if __name__ == "__main__":
    main()
