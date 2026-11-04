#!/usr/bin/env python
"""Embryo coupling diagnostic (2026-11-03; target-free: it reads no target stage and no score). The gate of embryo try
E10 in experiments/harness/TRIES.md, fixed there before this ran; with --test-panel, row N1 of the same file.

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
loads them, so the couplings are the ones the harness and the s1 build use at these seeds. Without --test-panel every
pair uses T2:embryo:val_interp's panel, the backtests' own (the test board was not in boards/ before the board intake),
and the report says so. The test pair's verdict at the s1 seed (1,000,001) is reported for the captain whatever the
gate says.

--test-panel (row N1, registered at 945de7d; a diagnostic, not a try: no target, no score): the test pair alone (E7.5,
E8.0 -> E7.75, t = 1/2), so only E7.5 and E8.0 are loaded, on the test board's own panel (T2:embryo:test_interp, read
through vec.boards.load_boards), with the entry's moves (fgw), at the same build seeds; at the s1 seed its coupling is
the one the s1 and s2 builds use. The same measurements, the s1 seed's draw reported separately, written to
diag_emb_coupling_testpanel.json and a short .txt that ends with N1's verdict against the values recorded on the
validation panel (RECORDED). Nothing it prints or writes is a cell count, a stage's or the draw's.

Usage (DGX, repo root, BLAS threads capped as docs/COMPUTE.md says):
  uv run --frozen --no-sync python -I -B experiments/harness/diag_emb_coupling.py \
    --out experiments/harness/results/diag_emb_coupling
  uv run --frozen --no-sync python -I -B experiments/harness/diag_emb_coupling.py --test-panel \
    --out experiments/harness/results/diag_emb_coupling_testpanel
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import textwrap
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

# Row N1 (--test-panel), fixed in TRIES.md before it runs. The values it is read against, as (low, high, what) at their
# two decimals: TRIES.md, "Diagnostics behind these choices", 2026-11-03 (the test pair's row of the table and the
# paragraph under it; 5baba46, T2:embryo:val_interp's panel).
TEST_BOARD = "T2:embryo:test_interp"
RECORDED = {
    "xcorr_min": (-0.50, -0.26, "fgw's smallest whitened anchor-partner cross-correlation, every draw"),
    "thin": (0.21, 0.25, "the moved cloud's thinnest SD over its RMS radius, every draw"),
    "s1_xcorr_min": (-0.31, -0.31, "the s1 seed's draw (1,000,001): the cross-correlation"),
    "s1_thin_ratio": (0.67, 0.67, "the s1 seed's draw: its thinnest axis over the thinner bracket's"),
    "s2_pool_thin": (0.47, 0.51, "the s2's pool (t_pos = 0): its thinnest SD over its RMS radius, every draw"),
}
ALSO_RECORDED = {"nn15_kept": (0.44, 0.47), "thin_ratio": (0.58, 0.70)}  # the table's 15-NN; row E11: shown, not judged
HOLD = 0.03  # "hold": every value within this of its recorded range (the row's expectation)
NOTE = 0.05  # "risk note": the s1 seed's thinnest axis or the s2 pool's leaves its range by more (the row's fixed rule)
NOTE_KEYS = ("s1_thin_ratio", "s2_pool_thin")
COUNT_FIELDS = ("coupling_nonzero",)  # an integer tied to the anchor count: left out of N1's outputs


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True, cwd=REPO_ROOT,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"}).stdout.strip()


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


def _inputs(before, after, genes: list[str]) -> tuple[np.ndarray, ...]:
    """A pair's coordinates, panel expression and merged labels, as the harness reads them."""
    Cb = np.asarray(before.obsm[COORD_KEY], float)[:, :3]
    Ca = np.asarray(after.obsm[COORD_KEY], float)[:, :3]
    Eb, Ea = dense(before[:, genes].X).astype(np.float32), dense(after[:, genes].X).astype(np.float32)
    lab_b, lab_a = before.obs["merged"].astype(str).to_numpy(), after.obs["merged"].astype(str).to_numpy()
    return Cb, Ca, Eb, Ea, lab_b, lab_a


def _nn_rows(seed: int, n_b: int, n_a: int) -> tuple[np.ndarray, np.ndarray]:
    """The cells of each stage in the neighbour check, from its own random stream."""
    rng = np.random.default_rng([seed, NN_STREAM])
    return rng.choice(n_b, min(NN_SAMPLE, n_b), replace=False), rng.choice(n_a, min(NN_SAMPLE, n_a), replace=False)


def _hundredths(v: float) -> int:
    """A value at the recorded precision, as it prints at two decimals."""
    return round(float(f"{v:.2f}") * 100)


def n1_verdict(values: dict[str, list[float]]) -> tuple[str, dict[str, float]]:
    """Row N1's verdict from the values per RECORDED key, compared in whole hundredths, and each key's departure from
    its range (0 inside it): "risk note" if the s1 seed's thinnest axis or the s2 pool's leaves its range by more than
    NOTE, "hold" if every value stays within HOLD of its range, otherwise "moved, within the note threshold"."""
    dep = {k: max(max(_hundredths(lo) - _hundredths(v), _hundredths(v) - _hundredths(hi), 0) for v in values[k])
           for k, (lo, hi, _) in RECORDED.items()}
    if any(dep[k] > _hundredths(NOTE) for k in NOTE_KEYS):
        verdict = "risk note"
    elif all(d <= _hundredths(HOLD) for d in dep.values()):
        verdict = "hold"
    else:
        verdict = "moved, within the note threshold"
    return verdict, {k: d / 100 for k, d in dep.items()}


def _span(v) -> str:
    lo, hi = f"{min(v):.2f}", f"{max(v):.2f}"
    return lo if lo == hi else f"{lo} to {hi}"


def run_n1(a, out: Path, commit: str, clean: bool) -> str:
    """Row N1 (--test-panel): the test pair on the test board's panel with the entry's moves. Writes the json and the
    txt; returns the txt, whose last line is the verdict."""
    key, sb_, sa_, t = PAIRS[2]
    boards = load_boards()
    if TEST_BOARD not in boards:
        sys.exit(f"{TEST_BOARD} is not in boards/: run the board intake first")
    if not a.seeds:
        sys.exit("--seeds: give at least one seed")
    genes = list(boards[TEST_BOARD].genes)
    stages = (load_stage(DATA / f"{s}.h5ad", genes=genes, drop_duplicates=True) for s in (sb_, sa_))
    Cb, Ca, Eb, Ea, lab_b, lab_a = _inputs(*stages, genes)
    draws = []
    for s in a.seeds:
        seed, t0 = BUILD_SEED_OFFSET + s, time.time()
        res = interpolate_clouds(Cb, Ca, Eb, Ea, t, N_DRAW, seed, target_rms=1.0, with_coupling=True)
        m = _measure(res, t, lab_b, lab_a, _nn_rows(seed, len(Cb), len(Ca)))
        f = {k: v for k, v in m.items() if k not in COUNT_FIELDS}
        draws.append({"seed": s, "build_seed": seed, "fgw": f, "seconds": round(time.time() - t0, 1)})
        print(f"{key} seed {s}, test panel: fgw xcorr {f['xcorr_min']:.2f} thin {f['shape_at_t'][-1]:.2f} "
              f"({f['thin_ratio']:.2f} of the thinner bracket's) s2 pool {f['shape_at_t_pos0'][-1]:.2f} "
              f"15-NN kept {f['nn15_kept']:.2f} ({draws[-1]['seconds']} s)", flush=True)
    fs = [d["fgw"] for d in draws]
    s1 = [d["fgw"] for d in draws if d["seed"] == S1_SEED]
    values = {"xcorr_min": [f["xcorr_min"] for f in fs], "thin": [f["shape_at_t"][-1] for f in fs],
              "s1_xcorr_min": [f["xcorr_min"] for f in s1], "s1_thin_ratio": [f["thin_ratio"] for f in s1],
              "s2_pool_thin": [f["shape_at_t_pos0"][-1] for f in fs]}
    full = sorted(a.seeds) == list(SEEDS)
    verdict, dep = n1_verdict(values) if full else ("partial run: no verdict", None)
    seeds = (f"{BUILD_SEED_OFFSET + SEEDS[0]:,}-{BUILD_SEED_OFFSET + SEEDS[-1]:,}" if full
             else ", ".join(f"{BUILD_SEED_OFFSET + s:,}" for s in a.seeds))
    labels = [f["label_agreement"] for f in fs if f["label_agreement"] is not None]
    mirrored = sum(v < 0 for v in values["xcorr_min"])
    lines = textwrap.wrap(
        f"N1 (experiments/harness/TRIES.md; a diagnostic, not a try): the embryo coupling on the test pair ({key}, "
        f"t = 1/2) on {TEST_BOARD}'s panel, the entry's moves (fgw), build seeds {seeds}. Target-free: reads {sb_} and "
        f"{sa_} only, no target and no score. Commit {commit}, clean tree: {clean}. Recorded: TRIES.md, diagnostics "
        f"of 2026-11-03 (5baba46, {GENES_BOARD}'s panel), compared at two decimals.", 118)
    lines += ["", f"  {'recorded':<16}{'now':<16}{'departure':<11}measure"]
    lines += [f"  {_span((lo, hi)):<16}{_span(values[k]) if values[k] else '-':<16}"
              f"{f'{dep[k]:.2f}' if dep else '-':<11}{what}" for k, (lo, hi, what) in RECORDED.items()]
    lines += [""] + textwrap.wrap(
        f"Reported, not judged: partners mirrored (cross-correlation below 0) in {mirrored} of {len(fs)} draws; "
        f"15-NN kept {_span([f['nn15_kept'] for f in fs])} (recorded "
        f"{_span(ALSO_RECORDED['nn15_kept'])}); thinnest axis over the thinner bracket's "
        f"{_span([f['thin_ratio'] for f in fs])} (recorded {_span(ALSO_RECORDED['thin_ratio'])}); merged-label "
        f"agreement {_span(labels) if labels else '-'}; the s1 seed's s2 pool "
        f"{_span([f['shape_at_t_pos0'][-1] for f in s1]) if s1 else '-'}.", 118)
    lines += textwrap.wrap(
        f"Rule (row N1, fixed before the run): hold if every value stays within {HOLD:.2f} of its recorded range; a "
        f"risk note to the captain if the s1 seed's thinnest axis or the s2 pool's leaves its recorded range by more "
        f"than {NOTE:.2f}; otherwise moved, within the note threshold.", 118)
    lines.append(f"verdict: {verdict}")
    txt = "\n".join(lines) + "\n"
    report = {"what": "row N1 of experiments/harness/TRIES.md (a diagnostic, not a try): the test pair on the test "
                      "board's panel, the entry's moves (fgw)",
              "commit": commit, "clean": clean, "genes_board": TEST_BOARD, "n_genes": len(genes), "pair": key,
              "before": sb_, "after": sa_, "t": t, "seeds": list(a.seeds), "s1_seed": S1_SEED,
              "recorded": {k: list(v) for k, v in RECORDED.items()}, "also_recorded": ALSO_RECORDED,
              "recorded_at": f"5baba46, {GENES_BOARD}'s panel (TRIES.md, diagnostics of 2026-11-03)",
              "rule": {"hold": HOLD, "note": NOTE, "note_keys": list(NOTE_KEYS)},
              "values": values, "departure": dep, "verdict": verdict,
              "s1_draw": next((d for d in draws if d["seed"] == S1_SEED), None), "draws": draws}
    (out / "diag_emb_coupling_testpanel.json").write_text(json.dumps(report, indent=1) + "\n")
    (out / "diag_emb_coupling_testpanel.txt").write_text(txt)
    return txt


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, nargs="*", default=list(SEEDS))
    ap.add_argument("--test-panel", action="store_true",
                    help="row N1 of experiments/harness/TRIES.md: the test pair alone, on T2:embryo:test_interp's "
                         "panel, the entry's moves only")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    commit, clean = _git("rev-parse", "HEAD"), _git("status", "--porcelain", "--untracked-files=no") == ""
    if a.test_panel:
        print(run_n1(a, out, commit, clean), end="")
        return
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
        Cb, Ca, Eb, Ea, lab_b, lab_a = _inputs(before, after, genes)
        draws = []
        for s in a.seeds:
            seed = BUILD_SEED_OFFSET + s
            t0 = time.time()
            nn_rows = _nn_rows(seed, len(Cb), len(Ca))
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
