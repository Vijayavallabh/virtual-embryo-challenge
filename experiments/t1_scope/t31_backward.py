#!/usr/bin/env python
"""Try T31 (experiments/harness/TRIES.md, registered before the run): which of E10.5's libraries best stands in for
another stage's sample across a REAL one-day gap? A labelled diagnostic; no upload changes and no cell's depth changes.

Why. Tries T27, T29 and T30 priced the base of T1's first upload (the pooled E10.5 chamber cells, or library 3's
alone) at zero time gap, and T30's deeper-truth scenario was made by thinning. The third cold review of 2026-11-12
said the open question (request R-B) lacks evidence at a real gap, and proposed predicting E10.5's library 3 from
E9.5 with two kinds of base. That cannot be built: E9.5's two libraries are one kind and one depth (the census of
docs/DATA_NOTES.md), and no released stage follows E10.5. The mirror image can: predict E9.5 from E10.5's libraries.
The scorer does not know the direction of time.

Geometry. The truth is every chamber cell of E9.5 (both libraries: embryo portions, about 14,000 counts per chamber
cell), as released. The reference, hence the floor, is every E10.5 cell. Every arm draws real E10.5 chamber cells at
the TRUTH's own chamber mix, so that only the base differs:
  copy          a uniform draw of the reference (no scope, E10.5's own mix)
  pooled_base   libraries in the proportions the E10.5 file holds, per category (the s1 base rule)
  lib3_base     library 3's cells only: the deepest E10.5 library (about 10,500), but a dissected-heart sample
  lib12_base    libraries 1 and 2 only: embryo portions like the truth, but at about a third of its depth
So the truth is deeper than every base, as in T30's deeper-truth scenario, with no thinning anywhere; and kind and
depth pull in opposite directions (the kind-matched base is the shallowest), which the test will not have if E12.5 is
a dissected-heart sample.

Reading, fixed before the run. T30's deeper-truth scenario "holds at a real gap" if lib3_base beats pooled_base by
more than 1.0 point with at least 8 of 10 seeds ahead; it is "contradicted" if pooled_base beats lib3_base by more than
1.0 with at least 8 of 10; otherwise "unclear". lib12_base against pooled_base is reported as kind against depth. No
recommendation follows from any outcome: request R-B keeps its condition (ruled on released data alone, or declined).

Usage (repo root; a few minutes):
  OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 uv run --frozen --no-sync python -I -B \
      experiments/t1_scope/t31_backward.py --out runs/t1_pricing/t31_<stamp>
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kind_gap as kg  # noqa: E402  (the library split and the balanced draw, unchanged)

from vec import harness, t1_scope  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import RANKED_METRICS, load_boards  # noqa: E402
from vec.core import counts_of, quotas  # noqa: E402
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample  # noqa: E402
from vec.prep import load_stage  # noqa: E402

REF_FILE, REF_SOURCE = "T1/E10.5_RNA.h5ad", "E10.5_RNA"
TRUTH_FILE, TRUTH_SOURCE = "T1/E9.5_RNA.h5ad", "E9.5_RNA"
ARMS = ("copy", "pooled_base", "lib3_base", "lib12_base")
READINGS = (("lib3_base", "pooled_base"), ("lib12_base", "pooled_base"), ("lib3_base", "lib12_base"))
MIN_GAIN = 1.0
MIN_AHEAD = 8
_D: dict = {}


def prepare(data_dir: Path) -> dict:
    spec = t1_scope.load_spec()
    chamber = [c for c in spec.categories if c in set(spec.in_scope("chamber"))]

    def cats(a, source):
        b = harness._shallow(a)
        t1_scope.classify(b, source=source, spec=spec)
        return np.asarray(b.obs["category"]).astype(str)

    ref = load_stage(data_dir / REF_FILE, lineage=False)
    truth = load_stage(data_dir / TRUTH_FILE, lineage=False)
    if list(ref.var_names) != list(truth.var_names):
        raise ValueError("the two stages do not share one gene order")
    cat_r, cat_t = cats(ref, REF_SOURCE), cats(truth, TRUTH_SOURCE)
    lib = kg.library_of(ref.obs_names)
    t_rows = np.flatnonzero(np.isin(cat_t, chamber))
    c = pd.Series(cat_t[t_rows]).value_counts()
    mix = {k: float(c.get(k, 0)) / float(c.sum()) for k in chamber}
    X = sparse.csr_matrix(ref.X)
    tot_r = np.asarray(counts_of(ref).sum(1)).ravel()
    tot_t = np.asarray(counts_of(truth).sum(1)).ravel()
    is3 = lib == kg.DISSECTED
    pools = {"3": {k: np.flatnonzero((cat_r == k) & is3) for k in chamber},
             "12": {k: np.flatnonzero((cat_r == k) & ~is3) for k in chamber}}
    lib3_share = {k: float(((cat_r == k) & is3).sum()) / float(max((cat_r == k).sum(), 1)) for k in chamber}
    boards = load_boards()
    harness.check_index(boards)
    board = replace(boards[kg.TWIN], key="pseudo:T1/E10.5->E9.5 backward", split="pseudo", anchors={},
                    **harness.pseudo_counts(len(t_rows), ref.n_obs))
    return {"X": X, "x_truth": sparse.csr_matrix(truth.X)[t_rows], "chamber": chamber, "mix": mix, "pools": pools,
            "lib3_share": lib3_share, "tot": tot_r, "board": board, "n_ref": int(ref.n_obs),
            "truth_cells": int(len(t_rows)), "truth_median_library": float(np.median(tot_t[t_rows])),
            "pool_sizes": {a: {k: int(len(v)) for k, v in p.items()} for a, p in pools.items()}}


def arm_rows(arm: str, n: int, seed: int) -> np.ndarray:
    d = _D
    rng = np.random.default_rng([seed, ARMS.index(arm)])
    if arm == "copy":
        return np.sort(kg.draw_balanced(d["n_ref"], n, rng))
    q = quotas(d["mix"], n)
    out = []
    for c in d["chamber"]:
        k = int(q.get(c, 0))
        if arm == "pooled_base":
            k3 = int(round(k * d["lib3_share"][c]))
            out += [_draw(d["pools"]["3"][c], k3, rng), _draw(d["pools"]["12"][c], k - k3, rng)]
        else:
            out.append(_draw(d["pools"]["3" if arm == "lib3_base" else "12"][c], k, rng))
    return np.sort(np.concatenate(out))


def _draw(pool: np.ndarray, k: int, rng) -> np.ndarray:
    if k <= 0:
        return np.array([], dtype=np.int64)
    if len(pool) == 0:
        raise ValueError("a category has no cell in a library it must be drawn from")
    return pool[kg.draw_balanced(len(pool), k, rng)]


def score_seed(seed: int) -> dict:
    d = _D
    t0 = time.time()
    X, board = d["X"], d["board"]
    target, ref = Stage(X=d["x_truth"]), Stage(X=X)
    f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, seed), target, ref, seed=seed, anchors="computed",
                    build=harness.BUILD, protocol=harness.PROTOCOL)
    if abs(f.score - 50.0) > 1e-9 or set(RANKED_METRICS["T1"]) - set(f.anchors):
        raise ValueError(f"seed {seed}: the floor scores {f.score}")
    n = harness.count_for(board, d["n_ref"])
    row = {"seed": seed, "scores": {}, "skills": {}, "median_library": {}, "max_repeat": {}}
    for arm in ARMS:
        rows = arm_rows(arm, n, harness.BUILD_SEED_OFFSET + seed)
        r = score_board(board, Stage(X=X[rows]), target, ref, seed=seed, anchors=f.anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        row["scores"][arm], row["skills"][arm] = r.score, r.skills
        row["median_library"][arm] = float(np.median(d["tot"][rows]))
        row["max_repeat"][arm] = int(np.unique(rows, return_counts=True)[1].max())
    row["seconds"] = round(time.time() - t0, 1)
    print(f"seed {seed}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["scores"].items()), flush=True)
    return row


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    ap.add_argument("--jobs", type=int, default=5)
    a = ap.parse_args(argv)
    commit = harness._git("rev-parse", "HEAD")
    clean = harness._git("status", "--porcelain", "--untracked-files=no") == ""
    _D.update(prepare(Path(a.data)))
    d = _D
    if a.jobs > 1:
        with mp.get_context("fork").Pool(a.jobs) as pool:
            rows = pool.map(score_seed, list(SEEDS))
    else:
        rows = [score_seed(s) for s in SEEDS]
    rows = sorted(rows, key=lambda r: r["seed"])
    readings = {}
    for x, y in READINGS:
        diff = np.array([r["scores"][x] - r["scores"][y] for r in rows])
        readings[f"{x} - {y}"] = {"mean": float(diff.mean()), "sd": float(diff.std(ddof=1)),
                                  "ahead": int((diff > 0).sum()),
                                  "skills": {m: float(np.mean([r["skills"][x][m] - r["skills"][y][m] for r in rows]))
                                             for m in rows[0]["skills"][x]}}
    lead = readings["lib3_base - pooled_base"]
    if lead["mean"] > MIN_GAIN and lead["ahead"] >= MIN_AHEAD:
        verdict = "holds at a real gap"
    elif lead["mean"] < -MIN_GAIN and lead["ahead"] <= len(rows) - MIN_AHEAD:
        verdict = "contradicted"
    else:
        verdict = "unclear"
    res = {"try": "T31", "commit": commit, "tree_clean": clean, "seeds": list(SEEDS), "arms": list(ARMS),
           "protocol": harness.PROTOCOL, "build": harness.BUILD, "frac": SUBSAMPLE_FRAC, "truth_mix": d["mix"],
           "truth_cells": d["truth_cells"], "truth_median_library": d["truth_median_library"],
           "library3_share_of_each_category": d["lib3_share"], "pool_sizes": d["pool_sizes"], "rows": rows,
           "mean_score": {x: float(np.mean([r["scores"][x] for r in rows])) for x in ARMS},
           "mean_skills": {x: {m: float(np.mean([r["skills"][x][m] for r in rows])) for m in rows[0]["skills"][x]}
                           for x in ARMS}, "readings": readings, "verdict": verdict}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "backward.json").write_text(json.dumps(res, indent=1))
    lines = [f"Try T31 (commit {commit}, clean {clean}): T30's deeper-truth scenario at a real one-day gap: {verdict}",
             "mean board score: " + ", ".join(f"{k} {v:.2f}" for k, v in res["mean_score"].items()),
             f"truth: {d['truth_cells']} E9.5 chamber cells, median library {d['truth_median_library']:.0f}; median "
             "library of each arm, first seed: "
             + ", ".join(f"{k} {v:.0f}" for k, v in rows[0]["median_library"].items())
             + "; largest number of times one cell is drawn, first seed: "
             + ", ".join(f"{k} {v}" for k, v in rows[0]["max_repeat"].items())]
    for k, v in readings.items():
        lines.append(f"  {k}: {v['mean']:+.2f} (sd {v['sd']:.2f}, ahead {v['ahead']} of {len(rows)}); "
                     + ", ".join(f"{m} {x:+.3f}" for m, x in v["skills"].items()))
    (out / "backward.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
