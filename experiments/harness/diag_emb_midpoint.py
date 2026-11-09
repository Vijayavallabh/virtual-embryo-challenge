#!/usr/bin/env python
"""N16 (experiments/harness/TRIES.md; 2026-11-09): the midpoint draw's target-free gate. Row E13 (fgw_midpoint) is
registered and run only if this gate passes on every pair; a no-op is not run and not counted.

Target-free: it reads the bracket stages of the three pairs only (E6.75 + E7.5, E7.25 + E8.0 on the backtests' panel;
E7.5 + E8.0, the unlocked E7.5 under A20, on the test board's panel), no target stage and no score, so it is not a try
(N1's precedent). It prints and writes ratios, not a cell count, a board limit or an anchor.

For each pair and build seed (vec.harness.BUILD_SEED_OFFSET + 0-9) it runs the entry's interpolation with
draw="midpoint" (the same call the candidate fgw_midpoint makes, at a fixed draw of N_DRAW cells that is not any board's
count) and measures:

- ESS/pool per bracket: Kish's effective size of the draw's expected cell shares over the bracket's pool
  (`record["midpoint"]["ess_frac"]`).
- The library check, on the draw's own expected shares rather than a finite sample, so sampling noise cannot decide it:
  the median library size N (vec.t2_build.library_sizes, the N the scorer normalises by) of the whole mixture and of
  each bracket, tilted against uniform. Reported, deciding nothing: the same ratio for the actual finite draw.
- The effect-size floor (the sharpener's fix 2): the MMD between the tilted draw and the same seed's uniform draw,
  against the MMD between that uniform draw and a second uniform draw of the same size (its own random streams). The
  MMD is the biased V-statistic with an RBF kernel on the build's joint PCs (the anchors' PCA of step 2, every cell
  projected), bandwidth the median pairwise distance of the pooled uniform draws (a fixed-size sample). Reported,
  deciding nothing: the same ratio for a no-op (the tilted draw's strata and quotas with flat weights), the level the
  gate's ratio takes when the tilt does nothing; it sits below 1, since such a draw is stratified and the uniform draws
  are not.
- Integrity: the drawn stratum shares of each bracket equal the pool's up to the rounding of one cell; the cloud guard
  flags nothing on the tilted draw.

The rule (row N16, fixed before the run; thresholds below), per pair: ESS/pool >= ESS_MIN in every bracket and seed; the
library median within LIBRARY_TOL in every seed (mixture and each bracket); the MMD ratio >= MMD_RATIO_MIN in at least
MMD_SEEDS_MIN of the 10 seeds; no cloud-guard flag and the strata exact in every seed. The gate passes only if every
pair passes. A partial run (--seeds) gives no verdict. Reported, deciding nothing: each bracket's mean signed offset s
along the local axis, the pool's and the tilted draw's (the size of the tilt and the fit signature of row E13).

Results go to experiments/harness/results/diag_emb_midpoint/ (a folder that holds no target-derived content), stamped
with the commit, the clean-tree flag (vec.provenance), the UTC time and the input files' hashes; it refuses to run if
outputs exist there (a rerun needs a new row).

Usage (DGX, repo root):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
    uv run --frozen --no-sync python -I -B experiments/harness/diag_emb_midpoint.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import textwrap
import time
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist, pdist

from vec import provenance
from vec.boards import REPO_ROOT, load_boards
from vec.core import draw_balanced
from vec.harness import BUILD_SEED_OFFSET
from vec.io import COORD_KEY, dense
from vec.prep import load_stage
from vec.t2_build import _draw_weighted, cloud_guard, interpolate_clouds, library_sizes, log_panel_counts

HERE = Path(__file__).resolve().parent
DATA = REPO_ROOT / "data" / "raw" / "T2"
OUT = HERE / "results" / "diag_emb_midpoint"
VAL_BOARD = "T2:embryo:val_interp"  # the backtests' panel (the harness's panel_board for both embryo pseudo-tasks)
TEST_BOARD = "T2:embryo:test_interp"  # the test build's panel; only its genes are read, never a count, limit or anchor
# (key, before, after, the mixing fraction t, the panel's board)
PAIRS = (("T2e/E6.75+E7.5->E7.25", "E6.75", "E7.5", 2 / 3, VAL_BOARD),
         ("T2e/E7.25+E8.0->E7.5", "E7.25", "E8.0", 1 / 3, VAL_BOARD),
         ("test/E7.5+E8.0->E7.75", "E7.5", "E8.0", 1 / 2, TEST_BOARD))
SEEDS = tuple(range(10))  # build seeds BUILD_SEED_OFFSET + s, the harness's
N_DRAW = 5000  # the gate's fixed draw (the coupling diagnostics' N_DRAW): not a board's count, never reported
TARGET_RMS = 1.0  # the scale is irrelevant to this gate; the size rule needs the target's day

# The gate, fixed in row N16 before this ran (the designer's and the sharpener's numbers, unchanged)
ESS_MIN = 0.5  # effective share of each bracket's pool, in every seed
LIBRARY_TOL = 0.05  # the median library size of the mixture and of each bracket, tilted against uniform, in every seed
MMD_RATIO_MIN = 3.0  # MMD(tilted, uniform) >= 3 x MMD(uniform, second uniform) ...
MMD_SEEDS_MIN = 8  # ... in at least 8 of the 10 seeds, on every pair

UNIFORM2_STREAMS = (10, 11)  # the second uniform draw (before, after), apart from the build's 0-9
BANDWIDTH_STREAM = 12  # the sample that sets the kernel bandwidth
FLAT_STREAMS = (13, 14)  # the no-op reference: the tilted draw's strata and quotas with flat weights (reported only)
BANDWIDTH_SAMPLE = 2000  # rows of the pooled uniform draws (a size of the computation)
BLOCK = 2048  # kernel-mean block size (memory only; the result does not depend on it)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def genes_sha256(genes) -> str:
    return hashlib.sha256("\n".join(map(str, genes)).encode()).hexdigest()


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """The median of `values` under non-negative `weights` (the smallest value whose cumulative share reaches 1/2)."""
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="stable")
    cum = np.cumsum(np.asarray(weights, dtype=np.float64)[order])
    if cum[-1] <= 0:
        raise ValueError("the weights are all zero")
    return float(values[order][np.searchsorted(cum, 0.5 * cum[-1])])


def rbf_gamma(X: np.ndarray, seed: int) -> float:
    """The RBF rate 1 / (2 sigma^2) with sigma the median pairwise distance of a BANDWIDTH_SAMPLE-row sample of X (its
    own stream, so numerator and denominator of the ratio use the same bandwidth for the seed)."""
    rng = np.random.default_rng([seed, BANDWIDTH_STREAM])
    sample = X[rng.choice(len(X), min(BANDWIDTH_SAMPLE, len(X)), replace=False)]
    sigma = float(np.median(pdist(sample)))
    return 1.0 / (2.0 * max(sigma, 1e-12) ** 2)


def kernel_mean(A: np.ndarray, B: np.ndarray, gamma: float) -> float:
    """The mean of exp(-gamma |a - b|^2) over all pairs of rows of A and B."""
    total = 0.0
    for lo in range(0, len(A), BLOCK):
        total += float(np.exp(-gamma * cdist(A[lo:lo + BLOCK], B, "sqeuclidean")).sum())
    return total / (len(A) * len(B))


def mmd_v(xx: float, yy: float, xy: float) -> float:
    """The biased (V-statistic) MMD from the three kernel means; never negative, so it can be a ratio's denominator."""
    return float(np.sqrt(max(xx + yy - 2.0 * xy, 0.0)))


def effect_size(tilted: np.ndarray, uni: np.ndarray, uni2: np.ndarray, seed: int,
                flat: np.ndarray | None = None) -> dict[str, float]:
    """MMD(tilted, uniform), MMD(uniform, second uniform) and their ratio: the same bandwidth for both (from the pooled
    uniform draws), the biased V-statistic with an RBF kernel (the row's measure, on the build's joint PCs). With `flat`
    (a draw with the tilted draw's stratum quotas but flat weights: a no-op) also the same ratio for it, the level the
    gate's ratio takes when the tilt does nothing (below 1: such a draw is stratified, the uniform draws are not);
    reported, deciding nothing."""
    gamma = rbf_gamma(np.vstack([uni, uni2]), seed)
    tt, uu, vv = kernel_mean(tilted, tilted, gamma), kernel_mean(uni, uni, gamma), kernel_mean(uni2, uni2, gamma)
    tu, uv = kernel_mean(tilted, uni, gamma), kernel_mean(uni, uni2, gamma)
    d_tilt, d_null = mmd_v(tt, uu, tu), mmd_v(uu, vv, uv)
    out = {"d_tilt": d_tilt, "d_null": d_null, "ratio": d_tilt / max(d_null, 1e-12)}
    if flat is not None:
        d_flat = mmd_v(kernel_mean(flat, flat, gamma), uu, kernel_mean(flat, uni, gamma))
        out["ratio_flat"] = d_flat / max(d_null, 1e-12)
    return out


def library_ratios(res, m: dict, t: float, lib_b: np.ndarray, lib_a: np.ndarray) -> dict[str, float]:
    """The tilted draw's median library size over the uniform draw's, on the draw's expected shares, for the mixture
    (the pool of both brackets) and for each bracket; and the same for the actual finite draw (the whole mixture)."""
    n_b, n_a = len(lib_b), len(lib_a)
    w = np.asarray(res.pool_weights, dtype=np.float64)
    w0 = np.r_[np.full(n_b, (1.0 - t) / n_b), np.full(n_a, t / n_a)]  # the uniform draw's expected shares
    lib = np.r_[lib_b, lib_a]

    def ratio(x: np.ndarray, tilted: np.ndarray, uniform: np.ndarray) -> float:
        return weighted_median(x, tilted) / weighted_median(x, uniform)

    out = {"mixture": ratio(lib, w, w0), "before": ratio(lib_b, w[:n_b], w0[:n_b]),
           "after": ratio(lib_a, w[n_b:], w0[n_b:])}
    drawn = np.r_[lib_b[res.rows[res.stage == 0]], lib_a[res.rows[res.stage == 1]]]
    uni = np.r_[lib_b[m["uniform_before"]], lib_a[m["uniform_after"]]]
    out["mixture_finite_draw"] = float(np.median(drawn) / np.median(uni))
    return out


def strata_exact(rows: np.ndarray, strata: np.ndarray) -> bool:
    """The drawn cells' stratum shares equal the pool's up to the rounding of one cell (the draw's quotas)."""
    k = int(strata.max()) + 1
    drawn = np.bincount(strata[rows], minlength=k) / len(rows)
    pool = np.bincount(strata, minlength=k) / len(strata)
    return bool(np.abs(drawn - pool).max() <= 1.0 / len(rows) + 1e-9)


def measure_seed(res, Eb: np.ndarray, Ea: np.ndarray, lib_b: np.ndarray, lib_a: np.ndarray, t: float, s: int,
                 seed: int) -> dict:
    """One build seed of one pair: every measure row N16 uses or reports, from the midpoint interpolation `res` (asked
    with with_coupling=True). Nothing here is a cell count."""
    m = res.coupling["midpoint"]
    rows_b, rows_a = res.rows[res.stage == 0], res.rows[res.stage == 1]
    u_b, u_a = np.asarray(m["uniform_before"]), np.asarray(m["uniform_after"])
    rec = res.record["midpoint"]
    ratios = library_ratios(res, m, t, lib_b, lib_a)
    # the build's joint PCs: the anchors' PCA of step 2 with every cell projected (the tilt's own, before residualising)
    Pb = (np.asarray(Eb, np.float64) - m["centre"]) @ m["loadings"]
    Pa = (np.asarray(Ea, np.float64) - m["centre"]) @ m["loadings"]
    tilted = np.vstack([Pb[rows_b], Pa[rows_a]])
    uni = np.vstack([Pb[u_b], Pa[u_a]])
    uni2 = np.vstack([Pb[draw_balanced(len(Pb), len(u_b), np.random.default_rng([seed, UNIFORM2_STREAMS[0]]))],
                      Pa[draw_balanced(len(Pa), len(u_a), np.random.default_rng([seed, UNIFORM2_STREAMS[1]]))]])
    # the no-op reference (reported, deciding nothing): the same strata and quotas with flat weights
    flat_b = _draw_weighted(m["strata_before"], np.ones(len(Pb)), len(u_b),
                            np.random.default_rng([seed, FLAT_STREAMS[0]]))[0]
    flat_a = _draw_weighted(m["strata_after"], np.ones(len(Pa)), len(u_a),
                            np.random.default_rng([seed, FLAT_STREAMS[1]]))[0]
    flat = np.vstack([Pb[flat_b], Pa[flat_a]])
    eff = effect_size(tilted, uni, uni2, seed, flat=flat)
    guard = cloud_guard(res.coords, res.pool_coords, seed=seed, source_weights=res.pool_weights)
    exact = strata_exact(rows_b, m["strata_before"]) and strata_exact(rows_a, m["strata_after"])

    def tilt_means(s_arr: np.ndarray, expected: np.ndarray) -> dict[str, float]:
        return {"pool": float(s_arr.mean()), "tilted": float((expected * s_arr).sum() / expected.sum())}

    return {"seed": s, "build_seed": seed, "ess_before": rec["ess_frac"]["before"],
            "ess_after": rec["ess_frac"]["after"],
            "library_ratio": {k: ratios[k] for k in ("mixture", "before", "after")},
            "library_max_shift": max(abs(ratios[k] - 1.0) for k in ("mixture", "before", "after")),
            "library_ratio_finite_draw": ratios["mixture_finite_draw"],
            "mmd_tilted_uniform": eff["d_tilt"], "mmd_uniform_uniform": eff["d_null"], "mmd_ratio": eff["ratio"],
            "mmd_ratio_flat": eff["ratio_flat"],
            "guard_flags": len(guard.flags), "guard_kinds": sorted(f.split(":")[0] for f in guard.flags),
            "strata_exact": exact,
            "s_mean": {"before": tilt_means(m["s_before"], m["expected_before"]),
                       "after": tilt_means(m["s_after"], m["expected_after"])}}


def pair_verdict(draws: list[dict]) -> dict:
    """Row N16's rule on one pair's per-seed measures: pass or fail, and what decided it."""
    ess = [min((d["ess_before"] or 0.0), (d["ess_after"] or 0.0)) for d in draws]
    ess_ok = all(v >= ESS_MIN for v in ess)
    library_ok = all(d["library_max_shift"] <= LIBRARY_TOL for d in draws)
    ratios = [d["mmd_ratio"] for d in draws]
    mmd_seeds = sum(r >= MMD_RATIO_MIN for r in ratios)
    guard_ok = all(d["guard_flags"] == 0 for d in draws)
    exact = all(d["strata_exact"] for d in draws)
    passes = ess_ok and library_ok and mmd_seeds >= MMD_SEEDS_MIN and guard_ok and exact
    return {"passes": bool(passes), "ess_min": float(min(ess)), "ess_ok": bool(ess_ok),
            "library_max_shift": float(max(d["library_max_shift"] for d in draws)), "library_ok": bool(library_ok),
            "mmd_ratio_min": float(min(ratios)), "mmd_ratio_median": float(np.median(ratios)),
            "mmd_ratio_max": float(max(ratios)), "mmd_seeds_at_floor": int(mmd_seeds),
            "mmd_ratio_flat_median": float(np.median([d.get("mmd_ratio_flat", np.nan) for d in draws])),
            "mmd_ok": bool(mmd_seeds >= MMD_SEEDS_MIN), "guard_clean": bool(guard_ok), "strata_exact": bool(exact)}


def gate_verdict(pairs: dict[str, dict], full: bool) -> str:
    """Pass only if every pair passes; a partial run gives no verdict."""
    if not full:
        return "partial run: no verdict"
    return "pass" if all(p["verdict"]["passes"] for p in pairs.values()) else "fail"


def _load_pair(sb_: str, sa_: str, genes: list[str]) -> tuple:
    """A pair's coordinates, panel expression, log panel depths and library sizes, as the harness loads them."""
    before, after = (load_stage(DATA / f"{s}.h5ad", genes=genes, drop_duplicates=True) for s in (sb_, sa_))
    Cb, Ca = (np.asarray(a.obsm[COORD_KEY], float)[:, :3] for a in (before, after))
    Eb, Ea = dense(before.X).astype(np.float32), dense(after.X).astype(np.float32)
    log_depth = (log_panel_counts(before, genes), log_panel_counts(after, genes))
    lib_b, lib_a = library_sizes(before.X, before.layers["counts"]), library_sizes(after.X, after.layers["counts"])
    return Cb, Ca, Eb, Ea, log_depth, np.asarray(lib_b, np.float64), np.asarray(lib_a, np.float64)


def _wrap(text: str) -> list[str]:
    return textwrap.wrap(" ".join(text.split()), 118)


def report_text(report: dict) -> str:
    """The gate's report: per pair the measures and its verdict, the rule, the verdict last (ratios only, no counts)."""
    s = report["stamp"]
    lines = _wrap(f"N16 (experiments/harness/TRIES.md; a diagnostic, not a try): the midpoint draw's target-free gate "
                  f"(build seeds {report['seed_span']}). It reads the pairs' bracket stages only, no target and no "
                  f"score. Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}.")
    for key, p in report["pairs"].items():
        v, n = p["verdict"], len(p["draws"])
        lines += ["", f"  {key} (t = {p['t']:.3g}; panel {p['panel_board']}): {'pass' if v['passes'] else 'FAIL'}"]
        lines.append(f"    ESS/pool per bracket: smallest {v['ess_min']:.3f} (needs {ESS_MIN} in every seed: "
                     f"{v['ess_ok']}); library shift, tilted over uniform, largest {v['library_max_shift']:.3f}")
        lines.append(f"    (needs at most {LIBRARY_TOL} in every seed: {v['library_ok']})")
        lines.append(f"    MMD(tilted, uniform) over MMD(uniform, second uniform): {v['mmd_ratio_min']:.2f} to "
                     f"{v['mmd_ratio_max']:.2f}, median {v['mmd_ratio_median']:.2f}")
        lines.append(f"    (needs {MMD_RATIO_MIN:g}x in at least {MMD_SEEDS_MIN} of {n} seeds: "
                     f"{v['mmd_seeds_at_floor']} of {n} at the floor: {v['mmd_ok']}); a no-op's level, reported only: "
                     f"median {v['mmd_ratio_flat_median']:.2f}")
        lines.append(f"    cloud guard clean in every seed: {v['guard_clean']}; stratum shares exact in every seed: "
                     f"{v['strata_exact']}")
    lines += [""] + _wrap(
        f"Rule (row N16, fixed before the run): ESS/pool at least {ESS_MIN} in every bracket and seed; the median "
        f"library size of the mixture and of each bracket within {LIBRARY_TOL:.0%} of the uniform draw's in every "
        f"seed (on the draw's expected shares); MMD(tilted, uniform) at least {MMD_RATIO_MIN:g}x MMD(uniform, second "
        f"uniform) on the build's joint PCs in at least {MMD_SEEDS_MIN} of 10 seeds; no cloud-guard flag; exact "
        "stratum shares; on every pair. A pass allows row E13; a fail means no run and no count. Reported only: the "
        "finite-draw library ratio, each bracket's mean signed offset, and the ratio of a no-op (the same strata and "
        "quotas with flat weights) in the JSON.")
    lines.append(f"verdict: {report['verdict']}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--seeds", type=int, nargs="*", default=list(SEEDS))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N16 has already written its outputs; a rerun needs a new TRIES row")
    if not a.seeds:
        sys.exit("--seeds: give at least one seed")
    boards = load_boards()
    for board in (VAL_BOARD, TEST_BOARD):
        if board not in boards:
            sys.exit(f"{board} is not in boards/: run the board intake first")
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    full = sorted(a.seeds) == list(SEEDS)
    report = {"what": __doc__.split("\n")[0], "pairs": {}, "seeds": list(a.seeds),
              "seed_span": (f"{BUILD_SEED_OFFSET + SEEDS[0]:,}-{BUILD_SEED_OFFSET + SEEDS[-1]:,}" if full
                            else ", ".join(f"{BUILD_SEED_OFFSET + s:,}" for s in a.seeds)),
              "rule": {"ESS_MIN": ESS_MIN, "LIBRARY_TOL": LIBRARY_TOL, "MMD_RATIO_MIN": MMD_RATIO_MIN,
                       "MMD_SEEDS_MIN": MMD_SEEDS_MIN, "streams_second_uniform": list(UNIFORM2_STREAMS),
                       "stream_bandwidth": BANDWIDTH_STREAM}}
    inputs: dict[str, dict] = {"stages": {}, "panels": {}}
    for key, sb_, sa_, t, board in PAIRS:
        genes = list(boards[board].genes)
        inputs["panels"][board] = {"genes_sha256": genes_sha256(genes)}
        for s_ in (sb_, sa_):
            if s_ not in inputs["stages"]:
                inputs["stages"][s_] = {"path": f"data/raw/T2/{s_}.h5ad", "sha256": sha256_file(DATA / f"{s_}.h5ad")}
        Cb, Ca, Eb, Ea, log_depth, lib_b, lib_a = _load_pair(sb_, sa_, genes)
        draws = []
        for s in a.seeds:
            t0 = time.time()
            seed = BUILD_SEED_OFFSET + s
            res = interpolate_clouds(Cb, Ca, Eb, Ea, t, N_DRAW, seed, target_rms=TARGET_RMS, draw="midpoint",
                                     log_depth=log_depth, with_coupling=True)
            row = measure_seed(res, Eb, Ea, lib_b, lib_a, t, s, seed)
            row["seconds"] = round(time.time() - t0, 1)
            draws.append(row)
            print(f"{key} seed {s}: ESS/pool {row['ess_before']:.3f} / {row['ess_after']:.3f}, library shift "
                  f"{row['library_max_shift']:.3f}, MMD ratio {row['mmd_ratio']:.2f}, guard flags "
                  f"{row['guard_flags']}, strata exact {row['strata_exact']} ({row['seconds']} s)", flush=True)
        report["pairs"][key] = {"before": sb_, "after": sa_, "t": float(t), "panel_board": board,
                                "verdict": pair_verdict(draws), "draws": draws}
    report["verdict"] = gate_verdict(report["pairs"], full)
    report["inputs"] = inputs
    report["stamp"] = {"commit": commit, "tree_clean": unclean == [],
                       "unclean_count": None if unclean is None else len(unclean),
                       "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out.mkdir(parents=True, exist_ok=True)
    (out / "diag_emb_midpoint.json").write_text(json.dumps(report, indent=1) + "\n")
    text = report_text(report)
    (out / "diag_emb_midpoint.txt").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
