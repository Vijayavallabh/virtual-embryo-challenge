#!/usr/bin/env python
"""Try T27 (experiments/harness/TRIES.md, registered before the run): what the kind of E10.5's libraries is worth to a
dissected-heart truth, read inside E10.5 at zero time gap. A labelled diagnostic for a ruling; no upload changes.

Why. The unlocked E10.5 file holds three libraries (the obs-name suffix). A target-free reading of that file
(2026-11-12) found two kinds: libraries 1 and 2 are embryo portions like E9.5 (43-47% of their cells in the chamber
categories, the rest pharyngeal, surrounding, progenitor and pericardial), and library 3 is a dissected-heart sample
(94% chamber cells, about twice as deep). Inside one chamber category, library 3 differs from library 1 or 2 in several
hundred genes at the scorer's 0.25 threshold, against about 100 between libraries 1 and 2. The E12.5 test hearts were
dissected (docs/DATA_NOTES.md), so library 3 is the like-for-like E10.5 sample. T1's s1 draws its chamber cells from
all three libraries at their pooled chamber mix; about two thirds of them are library 3's.

Geometry (E10.5 only, no time gap). Library 3 is split once into halves A and B (SPLIT_SEED). The truth is half A,
every cell as released (a real dissected-heart sample with its own non-chamber cells). The reference, as at test, is
the pooled stage: every E10.5 cell outside half A. Scored with the T1 board emulator (bake protocol, build main, 10%
working copies, computed anchors) on a synthetic board with the T1 validation board's panel.

Arms, n cells each by the count rule, real cells only, none from half A:
  copy              a uniform draw of the reference pool
  s1_base           chamber categories at the pooled E10.5 chamber mix, each category drawn from library 3 (half B)
                    and from libraries 1-2 in the proportion the full file holds: the entry's base as built
  s1_base_lib3mix   the same library proportions at library 3's own chamber mix (the mix alone)
  lib3_base         chamber categories at library 3's mix, every cell from library 3 (half B): the like-for-like base
  lib12_base        chamber categories at library 3's mix, every cell from libraries 1-2 (the other kind)
  lib3_whole        a uniform draw of half B with its non-chamber cells: library 3's sibling half, the upper reference

Readings (paired by seed): lib3_base - s1_base (the base's kind and mix together), s1_base_lib3mix - s1_base (the mix),
lib3_base - s1_base_lib3mix (the kind), lib3_whole - lib3_base (what dropping library 3's non-chamber cells costs).
Rule, fixed in the row: a ruling request for a library-3 base is prepared only if lib3_base beats s1_base by more than
1.0 board point, ahead in at least 8 of 10 seeds. The reading is an upper bound for the test, where two days of
development lie between the base and the truth.

Usage (repo root; about 15 minutes):
  OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 uv run --frozen --no-sync python -I -B \
      experiments/t1_scope/kind_gap.py --out runs/t1_kind_gap/<stamp>
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from vec import harness, t1_scope
from vec.adoption import SEEDS
from vec.boards import RANKED_METRICS, REPO_ROOT, load_boards
from vec.core import counts_of, draw_balanced, quotas
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.prep import load_stage

FILE = "T1/E10.5_RNA.h5ad"
SOURCE = "E10.5_RNA"
TWIN = "T1:val"
DISSECTED = "3"  # the library of dissected hearts (94% chamber cells)
SPLIT_SEED = 20261112
ARMS = ("copy", "s1_base", "s1_base_lib3mix", "lib3_base", "lib12_base", "lib3_whole")
READINGS = (("lib3_base", "s1_base"), ("s1_base_lib3mix", "s1_base"), ("lib3_base", "s1_base_lib3mix"),
            ("lib3_whole", "lib3_base"), ("lib12_base", "s1_base"))
MIN_GAIN = 1.0
MIN_AHEAD = 8
_D: dict = {}


def library_of(obs_names) -> np.ndarray:
    return pd.Series(np.asarray(obs_names, dtype=str)).str.split("_").str[-1].to_numpy()


def prepare(data_dir: Path) -> dict:
    a = load_stage(data_dir / FILE, lineage=False)
    spec = t1_scope.load_spec()
    b = harness._shallow(a)
    t1_scope.classify(b, source=SOURCE, spec=spec)
    cat = np.asarray(b.obs["category"]).astype(str)
    lib = library_of(a.obs_names)
    chamber = [c for c in spec.categories if c in set(spec.in_scope("chamber"))]
    rows3 = np.flatnonzero(lib == DISSECTED)
    perm = np.random.default_rng(SPLIT_SEED).permutation(rows3)
    half_a, half_b = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2:])
    in_a = np.zeros(a.n_obs, dtype=bool)
    in_a[half_a] = True
    ref_rows = np.flatnonzero(~in_a)
    other = np.flatnonzero(lib != DISSECTED)
    is_ch = np.isin(cat, chamber)

    def shares(rows) -> dict[str, float]:
        c = pd.Series(cat[rows][is_ch[rows]]).value_counts()
        return {k: float(c.get(k, 0)) / float(c.sum()) for k in chamber}

    pooled_mix = shares(np.arange(a.n_obs))  # the entry's: every E10.5 cell's chamber mix (harness.t1_chamber_shares)
    lib3_mix = shares(half_b)
    lib3_share = {c: float(((cat == c) & (lib == DISSECTED)).sum()) / float(max((cat == c).sum(), 1)) for c in chamber}
    X = sparse.csr_matrix(a.X)
    tot = np.asarray(counts_of(a).sum(1)).ravel()
    boards = load_boards()
    harness.check_index(boards)
    board = replace(boards[TWIN], key="pseudo:T1/E10.5 kind gap", split="pseudo", anchors={},
                    **harness.pseudo_counts(len(half_a), len(ref_rows)))
    pools = {"b": {c: half_b[cat[half_b] == c] for c in chamber},
             "o": {c: other[cat[other] == c] for c in chamber}}
    return {"X": X, "cat": cat, "lib": lib, "chamber": chamber, "half_a": half_a, "half_b": half_b,
            "ref_rows": ref_rows, "pooled_mix": pooled_mix, "lib3_mix": lib3_mix, "lib3_share": lib3_share,
            "tot": tot, "board": board, "pools": pools, "n_obs": int(a.n_obs)}


def _draw(pool: np.ndarray, k: int, rng) -> np.ndarray:
    if k <= 0:
        return np.array([], dtype=np.int64)
    if len(pool) == 0:
        raise ValueError("a category has no cell in a library it must be drawn from")
    return pool[draw_balanced(len(pool), k, rng)]


def arm_rows(arm: str, n: int, seed: int) -> np.ndarray:
    d = _D
    rng = np.random.default_rng([seed, ARMS.index(arm)])
    if arm == "copy":
        return np.sort(d["ref_rows"][draw_balanced(len(d["ref_rows"]), n, rng)])
    if arm == "lib3_whole":
        return np.sort(d["half_b"][draw_balanced(len(d["half_b"]), n, rng)])
    mix = d["pooled_mix"] if arm == "s1_base" else d["lib3_mix"]
    q = quotas(mix, n)
    out = []
    for c in d["chamber"]:
        k = int(q.get(c, 0))
        if arm in ("s1_base", "s1_base_lib3mix"):
            k3 = int(round(k * d["lib3_share"][c]))
            out += [_draw(d["pools"]["b"][c], k3, rng), _draw(d["pools"]["o"][c], k - k3, rng)]
        elif arm == "lib3_base":
            out.append(_draw(d["pools"]["b"][c], k, rng))
        else:  # lib12_base
            out.append(_draw(d["pools"]["o"][c], k, rng))
    return np.sort(np.concatenate(out))


def score_seed(seed: int) -> dict:
    d = _D
    t0 = time.time()
    X, board = d["X"], d["board"]
    target, ref = Stage(X=X[d["half_a"]]), Stage(X=X[d["ref_rows"]])
    f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, seed), target, ref, seed=seed, anchors="computed",
                    build=harness.BUILD, protocol=harness.PROTOCOL)
    if abs(f.score - 50.0) > 1e-9 or set(RANKED_METRICS["T1"]) - set(f.anchors):
        raise ValueError(f"seed {seed}: the floor scores {f.score}")
    n = harness.count_for(board, len(d["ref_rows"]))
    row = {"seed": seed, "scores": {}, "skills": {}, "median_library": {}}
    for arm in ARMS:
        rows = arm_rows(arm, n, harness.BUILD_SEED_OFFSET + seed)
        r = score_board(board, Stage(X=X[rows]), target, ref, seed=seed, anchors=f.anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        row["scores"][arm], row["skills"][arm] = r.score, r.skills
        row["median_library"][arm] = float(np.median(d["tot"][rows]))
    row["seconds"] = round(time.time() - t0, 1)
    print(f"seed {seed}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["scores"].items()), flush=True)
    return row


def summarise(rows: list[dict]) -> dict:
    out = {"mean_score": {a: float(np.mean([r["scores"][a] for r in rows])) for a in ARMS},
           "mean_skills": {a: {m: float(np.mean([r["skills"][a][m] for r in rows])) for m in rows[0]["skills"][a]}
                           for a in ARMS},
           "median_library_first_seed": rows[0]["median_library"], "readings": {}}
    for x, y in READINGS:
        diff = np.array([r["scores"][x] - r["scores"][y] for r in rows])
        out["readings"][f"{x} - {y}"] = {
            "mean": float(diff.mean()), "sd": float(diff.std(ddof=1)), "ahead": int((diff > 0).sum()),
            "seeds": int(len(diff)),
            "skills": {m: float(np.mean([r["skills"][x][m] - r["skills"][y][m] for r in rows]))
                       for m in rows[0]["skills"][x]}}
    lead = out["readings"]["lib3_base - s1_base"]
    out["rule"] = {"min_gain": MIN_GAIN, "min_ahead": MIN_AHEAD,
                   "ruling_request": bool(lead["mean"] > MIN_GAIN and lead["ahead"] >= MIN_AHEAD)}
    return out


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
    res = {"commit": commit, "tree_clean": clean, "file": FILE, "split_seed": SPLIT_SEED, "seeds": list(SEEDS),
           "protocol": harness.PROTOCOL, "build": harness.BUILD, "frac": SUBSAMPLE_FRAC,
           "pooled_chamber_mix": d["pooled_mix"], "library3_chamber_mix": d["lib3_mix"],
           "library3_share_of_each_category": d["lib3_share"],
           "median_library_whole_stage": float(np.median(d["tot"])), "rows": rows, "summary": summarise(rows)}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "kind_gap.json").write_text(json.dumps(res, indent=1))
    s = res["summary"]
    lines = [f"Try T27 (commit {commit}, clean {clean}): ruling request by the rule: {s['rule']['ruling_request']}",
             "mean board score: " + ", ".join(f"{k} {v:.2f}" for k, v in s["mean_score"].items())]
    for k, v in s["readings"].items():
        lines.append(f"  {k}: {v['mean']:+.2f} (sd {v['sd']:.2f}, ahead {v['ahead']} of {v['seeds']}); "
                     + ", ".join(f"{m} {x:+.3f}" for m, x in v["skills"].items()))
    lines.append("median library size, first seed: "
                 + ", ".join(f"{k} {v:.0f}" for k, v in s["median_library_first_seed"].items())
                 + f"; whole stage {res['median_library_whole_stage']:.0f}")
    (out / "kind_gap.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
