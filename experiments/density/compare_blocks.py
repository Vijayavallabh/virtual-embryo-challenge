"""Compare the contiguous-block conditions of the first run (v1: every block nested around one centre per seed,
chosen so that the 8x block fits) with the rerun (v2: one centre per count level, uniform over the positions where
that block fits). head and refhead do not depend on the centre and must reproduce exactly.

Usage (DGX, repo root):
  uv run python experiments/density/compare_blocks.py --old runs/density/raw --new runs/density/raw_v2 \
    --out runs/density/summary_v2
Writes blocks_v1_vs_v2.csv (per target x pattern x ratio) and blocks_calls.csv (one row per block call: position
and paired delta, no cell-level data), and prints the tables used in the report.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PATTERNS = ["block", "head", "refblock", "refhead"]
TARGETS = ["heart_E9.5", "heart_E8.75", "embryo_E8.0"]


def load(raw: Path) -> pd.DataFrame:
    rows = []
    for p in sorted(raw.glob("*.jsonl")):
        for line in p.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows.append(dict(target=r["target"], seed=r["seed"], condition=r["condition"], pattern=r["pattern"],
                                 ratio=float(r["ratio"]), n=r["n_pred"], pool=r["pool_cells"],
                                 block_start=r.get("block_start"), score=r["score"],
                                 d2=r["skills"]["d2_shape"], dice=r["skills"]["occupancy_dice"],
                                 nbhd=r["skills"]["neighborhood_mmd"], shape=r["group_skill"]["shape_growth"]))
    return pd.DataFrame(rows).drop_duplicates(["target", "seed", "condition"], keep="last")


def paired(df: pd.DataFrame, base: pd.DataFrame) -> pd.DataFrame:
    """Blocks with the paired delta against random sampling of the same source at the same count and seed."""
    b = df[df["pattern"].isin(PATTERNS)].copy()
    b["base_pattern"] = np.where(b["pattern"].str.startswith("ref"), "refstage", "random")
    ref = base[base["pattern"].isin(["random", "refstage"])].set_index(["target", "seed", "pattern", "ratio"])
    b = b.join(ref[["score"]], on=["target", "seed", "base_pattern", "ratio"], rsuffix="_samen")
    b["delta"] = b["score"] - b["score_samen"]
    span = (b["pool"] - b["n"]).astype(float)
    b["rel_pos"] = np.where(span > 0, b["block_start"] / span.where(span > 0, 1.0), np.nan)
    return b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", default="runs/density/raw")
    ap.add_argument("--new", default="runs/density/raw_v2")
    ap.add_argument("--out", default="runs/density/summary_v2")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    old, new = load(Path(args.old)), load(Path(args.new))
    # random/refstage rows are identical in both directories (seed_v2.py copies them); take them from v2.
    v1, v2 = paired(old, new).assign(version="v1"), paired(new, new).assign(version="v2")
    calls = pd.concat([v1, v2], ignore_index=True)
    calls[["version", "target", "seed", "pattern", "ratio", "n", "pool", "block_start", "rel_pos", "score",
           "score_samen", "delta", "d2", "dice", "nbhd", "shape"]].to_csv(out / "blocks_calls.csv", index=False,
                                                                           float_format="%.6g")

    for pat in ("head", "refhead"):
        j = v1[v1["pattern"] == pat].merge(v2[v2["pattern"] == pat], on=["target", "seed", "condition"])
        print(f"{pat}: {len(j)} calls in both runs; max |score v2 - v1| = {np.abs(j.score_x - j.score_y).max():.3g}")

    g = calls.groupby(["version", "target", "pattern", "ratio"])
    summ = pd.DataFrame({
        "n_seeds": g.size(), "score_mean": g["score"].mean(), "delta_mean": g["delta"].mean(),
        "delta_sd": g["delta"].std(ddof=1), "seeds_lose": g["delta"].apply(lambda s: int((s < 0).sum())),
        "delta_min": g["delta"].min(), "delta_max": g["delta"].max(), "rel_pos_min": g["rel_pos"].min(),
        "rel_pos_max": g["rel_pos"].max(), "block_start_min": g["block_start"].min(),
        "block_start_max": g["block_start"].max(), "d2_mean": g["d2"].mean(), "dice_mean": g["dice"].mean(),
        "shape_mean": g["shape"].mean()}).reset_index()
    wide = summ.pivot_table(index=["target", "pattern", "ratio"], columns="version",
                            values=["delta_mean", "delta_sd", "seeds_lose", "delta_min", "delta_max", "rel_pos_min",
                                    "rel_pos_max", "block_start_min", "block_start_max", "n_seeds"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    wide.to_csv(out / "blocks_v1_vs_v2.csv", index=False, float_format="%.6g")

    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    print("\n== paired delta vs random sampling of the same source at the same count: mean ± sd [seeds < 0] "
          "(position = block start / (pool - n), 0 = first rows)")
    for pat in ("block", "refblock"):
        print(f"-- {pat}")
        for t in TARGETS:
            for r in sorted(calls["ratio"].unique()):
                cells = []
                for v in ("v1", "v2"):
                    z = summ[(summ.version == v) & (summ.target == t) & (summ.pattern == pat) & (summ.ratio == r)]
                    if len(z):
                        z = z.iloc[0]
                        cells.append(f"{v}: {z.delta_mean:+6.1f} ± {z.delta_sd:4.1f} [{z.seeds_lose}/{z.n_seeds}] "
                                     f"range {z.delta_min:+6.1f}..{z.delta_max:+6.1f} pos {z.rel_pos_min:.2f}-"
                                     f"{z.rel_pos_max:.2f}")
                print(f"{t:12s} {r:5g}x  " + " | ".join(cells))

    print("\n== pooled over the three targets (mean of target means; seeds < 0 over all calls)")
    for r in (1.0, 8.0):
        for pat in PATTERNS:
            line = []
            for v in ("v1", "v2"):
                z = calls[(calls.version == v) & (calls.pattern == pat) & (calls.ratio == r)]
                tm = z.groupby("target")["delta"].mean()
                line.append(f"{v} {tm.mean():+6.1f} (targets {tm.min():+.1f} to {tm.max():+.1f}; "
                            f"{int((z.delta < 0).sum())}/{len(z)} lose)")
            print(f"{r:g}x {pat:9s} " + " | ".join(line))
    for v in ("v1", "v2"):
        z = calls[(calls.version == v) & calls.pattern.isin(["block", "refblock"])]
        print(f"{v}: block/refblock calls that beat random at the same count: {int((z.delta > 0).sum())}/{len(z)}; "
              + ", ".join(f"{q.target} {q.pattern} {q.ratio:g}x s{q.seed} {q.delta:+.1f}"
                          for q in z[z.delta > 0].itertuples()))

    print("\n== v2: does the position of the block matter? Spearman(distance from mid-file, delta), within "
          "(target, ratio) after removing their means; and mean delta of blocks in the outer 10% of positions")
    for pat in ("block", "refblock"):
        z = v2[(v2.pattern == pat) & v2.rel_pos.notna()].copy()
        z["resid"] = z["delta"] - z.groupby(["target", "ratio"])["delta"].transform("mean")
        z["edge_dist"] = (z["rel_pos"] - 0.5).abs()
        rho, p = spearmanr(z["edge_dist"], z["resid"])
        edge = z[(z.rel_pos < 0.1) | (z.rel_pos > 0.9)]
        mid = z[(z.rel_pos >= 0.1) & (z.rel_pos <= 0.9)]
        print(f"{pat}: n={len(z)} rho={rho:+.2f} (p={p:.2g}); residual delta outer 10%: {edge.resid.mean():+.1f} "
              f"(n={len(edge)}), middle 80%: {mid.resid.mean():+.1f} (n={len(mid)})")


if __name__ == "__main__":
    main()
