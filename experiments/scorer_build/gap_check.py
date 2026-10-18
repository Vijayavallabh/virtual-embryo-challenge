#!/usr/bin/env python
"""Proper-minus-mirror occupancy_dice gap on a released stage's split halves, under both scorer builds.

Released training stages only (default T2 E9.5, the heart stage nearest the E10.5 validation board). Two parts:

1. Halves. For each half size and seed, draw 2 x half cells, split them into a truth half A and a prediction half B
   (as the board ceiling does) and record: Dice in proper mode (main build), Dice in mirror mode (main build on the
   z-mirrored prediction) and their gap; the prefix build's Dice and the mode it took (SVD parities); the theory
   checks (prefix Dice = the parity rule; a z-mirrored prediction scores exactly like the original under prefix;
   main Dice = proper); and the prediction's local parity (metric seed 0, n_true = |A|). Half size "wc" is the
   emulator's own working copy (10% of the stage, vec.evaluate.split_half seed 0); numeric sizes draw 2 x half cells
   with default_rng(seed) and split with veckit's split_half(seed).
2. Emulator. vec.evaluate.score_board on a pseudo-board (target = --stage, reference = --ref, gene panel of
   --genes-board, computed anchors) under both builds. The candidate is the other half of the working copy with its
   rows shuffled and a random proper rotation, so its parity coin is independent of the ceiling's. The panel's Dice
   must equal the proper value under main and the parity-rule value under prefix; every metric that does not
   canonicalise must be identical across builds.

Outputs aggregate numbers only: gap_check.json, gap_summary.csv and gap_check.txt in --out.
Run on the DGX (repo root), in tmux with a log under runs/logs/:
  OMP_NUM_THREADS=8 uv run python experiments/scorer_build/gap_check.py --out runs/scorer_build/gap_check
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from vec.boards import BOARDS_DIR, REPO_ROOT, Board, read_genes
from vec.evaluate import (
    SUBSAMPLE_FRAC,
    Stage,
    handedness,
    mirror_z,
    score_board,
    scorer_build,
    split_half,
    subsample,
)

BUILD_DEPENDENT = {"occupancy_dice", "sliced_wasserstein", "_sw_flip_spread"}


def load_coords(path) -> np.ndarray:
    import h5py

    with h5py.File(path, "r") as f:
        o = f["obsm"]["spatial_3D"]
        if isinstance(o, h5py.Dataset):
            C = o[()]
        else:
            C = np.stack([o[c][()] for c in list(o.attrs.get("column-order"))], 1)
    return np.asarray(C, dtype=np.float32)[:, :3]


def rand_rot(rng) -> np.ndarray:
    """A uniformly random proper rotation (unit quaternion)."""
    q = rng.normal(size=4)
    a, b, c, d = q / np.linalg.norm(q)
    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]])


def half_row(B, A) -> dict:
    """Dice of prediction B against truth A under both builds, with the theory checks."""
    import common.shape_metrics as sm

    h = handedness(B, A, seed=0)
    with scorer_build("main"):
        main = sm.occupancy_dice(B, A, seed=0)[0]
    with scorer_build("prefix"):
        pre = sm.occupancy_dice(B, A, seed=0)[0]
        pre_twin = sm.occupancy_dice(mirror_z(B), A, seed=0)[0]
    return {**{k: h[k] for k in ("dice_proper", "dice_mirror", "dice_gap", "parity_pred", "parity_truth",
                                 "prefix_mode", "local_parity")},
            "dice_prefix": round(pre, 4), "n_half": int(len(A)),
            "main_equals_proper": round(main, 4) == h["dice_proper"],
            "prefix_equals_parity_rule": round(pre, 4) == h["dice_prefix"],
            "mirror_twin_unchanged_under_prefix": abs(pre_twin - pre) < 1e-12}


def halves_part(C, half_sizes, seeds) -> list[dict]:
    from common.core_metrics import split_half as veckit_split_half

    stage = Stage(X=np.zeros((len(C), 1), np.float32), coords=C)
    rows = []
    for hs in half_sizes:
        for s in range(seeds):
            if hs == "wc":
                A, B = (x.coords for x in split_half(subsample(stage, SUBSAMPLE_FRAC, s), seed=0))
            else:
                S = C[np.random.default_rng(s).choice(len(C), 2 * int(hs), replace=False)]
                ia, ib = veckit_split_half(len(S), seed=s)
                A, B = S[ia], S[ib]
            rows.append({"half_size": hs, "seed": s, **half_row(B, A)})
        print(f"half size {hs} done", flush=True)
    return rows


def summarise(rows) -> list[dict]:
    out = []
    for hs in dict.fromkeys(r["half_size"] for r in rows):
        R = [r for r in rows if r["half_size"] == hs]
        gap = np.array([r["dice_gap"] for r in R])
        out.append({"half_size": hs, "n_half": R[0]["n_half"], "seeds": len(R),
                    "dice_proper_mean": round(float(np.mean([r["dice_proper"] for r in R])), 4),
                    "dice_mirror_mean": round(float(np.mean([r["dice_mirror"] for r in R])), 4),
                    "gap_mean": round(float(gap.mean()), 4),
                    "gap_sd": round(float(gap.std(ddof=1)), 4) if len(gap) > 1 else 0.0,
                    "gap_min": round(float(gap.min()), 4), "gap_max": round(float(gap.max()), 4),
                    "prefix_mirror_mode_frac": round(float(np.mean([r["prefix_mode"] == "mirror" for r in R])), 3),
                    "main_equals_proper": sum(r["main_equals_proper"] for r in R),
                    "prefix_equals_parity_rule": sum(r["prefix_equals_parity_rule"] for r in R),
                    "mirror_twin_unchanged_under_prefix": sum(r["mirror_twin_unchanged_under_prefix"] for r in R)})
    return out


def load_stage(path, genes) -> Stage:
    import anndata as ad

    return Stage.from_anndata(ad.read_h5ad(path), genes, need_coords=True)


def _diagnostics(res) -> dict:
    """A candidate's handedness diagnostics without the build-specific consistency flag."""
    return {k: v for k, v in res.extras["handedness"].items() if k != "panel_dice_matches_build"}


def emulator_part(stage_path, ref_path, genes_board: str, seeds) -> list[dict]:
    genes = read_genes(BOARDS_DIR / f"{genes_board.replace(':', '__')}.genes.txt")
    # a pseudo-board: only the gene panel is borrowed; no index.json field (anchors, counts) is read
    board = Board(key=f"pseudo:{Path(stage_path).stem}<-{Path(ref_path).stem}", task="T2", genes=genes,
                  needs_coords=True)
    target, ref = load_stage(stage_path, genes), load_stage(ref_path, genes)
    rows = []
    for s in range(seeds):
        _, other = split_half(subsample(target, SUBSAMPLE_FRAC, s), seed=0)
        rng = np.random.default_rng(500 + s)
        cand = other.take(rng.permutation(other.n))
        cand = Stage(X=cand.X, coords=(cand.coords.astype(float) @ rand_rot(rng).T).astype(np.float32),
                     celltype=cand.celltype)
        res = {b: score_board(board, cand, target, ref, seed=s, build=b) for b in ("main", "prefix")}
        m, p = res["main"], res["prefix"]
        same_other = all(m.metrics[k] == p.metrics[k] or (m.metrics[k] is None and p.metrics[k] is None)
                         for k in m.metrics if k not in BUILD_DEPENDENT)
        hand = m.extras["handedness"]
        rows.append({"seed": s, "score_main": round(m.score, 3), "score_prefix": round(p.score, 3),
                     "score_diff_prefix_minus_main": round(p.score - m.score, 3),
                     "dice_main": m.metrics["occupancy_dice"], "dice_prefix": p.metrics["occupancy_dice"],
                     "skill_dice_main": round(m.skills["occupancy_dice"], 4),
                     "skill_dice_prefix": round(p.skills["occupancy_dice"], 4),
                     "dice_anchors_main": m.anchors["occupancy_dice"],
                     "dice_anchors_prefix": p.anchors["occupancy_dice"],
                     "handedness": hand, "anchor_handedness": m.extras.get("anchor_handedness"),
                     "panel_dice_matches_build": [m.extras["handedness"]["panel_dice_matches_build"],
                                                  p.extras["handedness"]["panel_dice_matches_build"]],
                     "handedness_same_under_both_builds": _diagnostics(m) == _diagnostics(p),
                     "other_metrics_identical": same_other})
        print(f"emulator seed {s} done: {rows[-1]['score_main']} (main) {rows[-1]['score_prefix']} (prefix)",
              flush=True)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", default=str(REPO_ROOT / "data" / "raw" / "T2" / "E9.5.h5ad"))
    ap.add_argument("--ref", default=str(REPO_ROOT / "data" / "raw" / "T2" / "E8.75.h5ad"))
    ap.add_argument("--genes-board", default="T2:heart:val_extrap", help="board whose gene list the pseudo-board uses")
    ap.add_argument("--half-sizes", default="wc,1185,2687,5374,8393,13000")
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--emulator-seeds", type=int, default=3, help="0 skips the score_board part")
    ap.add_argument("--out", default=str(REPO_ROOT / "runs" / "scorer_build" / "gap_check"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    C = load_coords(a.stage)
    sizes = [x if x == "wc" else int(x) for x in a.half_sizes.split(",") if x.strip()]
    rows = halves_part(C, sizes, a.seeds)
    summary = summarise(rows)
    report = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "stage": Path(a.stage).name, "n_cells": int(len(C)),
              "halves_summary": summary, "halves_rows": rows}
    if a.emulator_seeds:
        report["emulator"] = emulator_part(a.stage, a.ref, a.genes_board, a.emulator_seeds)
        report["emulator_reference"] = Path(a.ref).name
    report["sec"] = round(time.time() - t0, 1)
    (out / "gap_check.json").write_text(json.dumps(report, indent=1, default=float) + "\n")
    with open(out / "gap_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    lines = [f"gap check {report['generated']}: {report['stage']} ({report['n_cells']} cells), {a.seeds} seeds"]
    for r in summary:
        lines.append(f"half {r['half_size']} ({r['n_half']} cells): proper {r['dice_proper_mean']} mirror "
                     f"{r['dice_mirror_mean']} gap {r['gap_mean']} (sd {r['gap_sd']}, range {r['gap_min']} to "
                     f"{r['gap_max']}); prefix mirror mode {r['prefix_mirror_mode_frac']}; checks main=proper "
                     f"{r['main_equals_proper']}/{r['seeds']}, prefix=rule {r['prefix_equals_parity_rule']}/"
                     f"{r['seeds']}, twin {r['mirror_twin_unchanged_under_prefix']}/{r['seeds']}")
    for e in report.get("emulator", []):
        lines.append(f"emulator seed {e['seed']}: score main {e['score_main']} prefix {e['score_prefix']} "
                     f"(diff {e['score_diff_prefix_minus_main']}); Dice main {e['dice_main']} prefix "
                     f"{e['dice_prefix']} (proper {e['handedness']['dice_proper']}, mirror "
                     f"{e['handedness']['dice_mirror']}, mode {e['handedness']['prefix_mode']}, local parity "
                     f"{e['handedness']['local_parity']:+d}); panel Dice matches build {e['panel_dice_matches_build']}"
                     f"; other metrics identical {e['other_metrics_identical']}")
    lines.append(f"{report['sec']} s")
    (out / "gap_check.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
