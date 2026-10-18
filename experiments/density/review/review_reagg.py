"""Reviewer's independent re-aggregation of the experiment's raw JSONL (runs/density/raw), without summarise.py.

Usage (DGX, repo root): uv run python experiments/density/review/review_reagg.py [raw_dir]

Copied into the repo from the independent review of 2026-10-19 (originally run from runs/density_review/);
only script paths changed, plus lint-only edits (line wraps; variable O renamed). Outputs stay on the DGX in
runs/density_review/.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("runs/density/raw")
T = ["heart_E9.5", "heart_E8.75", "embryo_E8.0"]
M = ["de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
     "neighborhood_mmd"]

rows = []
for p in sorted(RAW.glob("*.jsonl")):
    for line in p.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            o = dict(target=r["target"], seed=r["seed"], cond=r["condition"], n=r["n_pred"], tc=r["truth_cells"],
                     score=r["score"], nn_own=r["nn_own"], nn_truth=r["nn_truth"], jr=r["jitter_rms"])
            for m in M:
                o["s_" + m] = r["skills"][m]
                o["v_" + m] = r["metrics"][m]
                o["f_" + m] = r["floor"][m]
                o["c_" + m] = r["ceiling"][m]
            o["vox_nn"] = r["metrics"].get("_dice_voxel_over_nn")
            o["vox_nn_ceil"] = r["ceiling_extra"].get("_dice_voxel_over_nn")
            o.update({"g_" + k: v for k, v in r["group_skill"].items()})
            rows.append(o)
df = pd.DataFrame(rows)
print("rows", len(df), "unique (target,seed,cond)", len(df.drop_duplicates(["target", "seed", "cond"])))
S = df.pivot_table(index=["target", "seed"], columns="cond", values="score")


def paired(c, base):
    d = S[c] - S[base]
    return d


def show(title, pairs):
    print(f"\n== {title}")
    for c, base in pairs:
        line = []
        for t in T:
            d = paired(c, base).loc[t]
            line.append(f"{t}: {d.mean():+.2f} (sd {d.std(ddof=1):.2f}) [{int((d < 0).sum())}/5] "
                        f"min {d.min():+.2f} max {d.max():+.2f}")
        allv = paired(c, base)
        print(f"{c} - {base}: " + " | ".join(line) + f" || pooled mean of target means "
              f"{np.mean([paired(c, base).loc[t].mean() for t in T]):+.2f}, seeds<0 {int((allv < 0).sum())}/15")


show("near-perfect count sweep vs random_1x", [(f"random_{r}x", "random_1x") for r in ("0.25", "0.5", "2", "4", "8")])
show("floor-quality count sweep vs refstage_1x",
     [(f"refstage_{r}x", "refstage_1x") for r in ("0.25", "0.5", "2", "4", "8")])
show("near-perfect 2x vs 4x and 8x", [("random_2x", "random_4x"), ("random_2x", "random_8x")])
show("floor-quality 2x vs 4x and 8x", [("refstage_2x", "refstage_4x"), ("refstage_2x", "refstage_8x")])
show("blocks vs random same n", [(f"{p}_{r}x", f"random_{r}x") for p in ("block", "head") for r in ("1", "8")])
show("ref blocks vs refstage same n", [(f"{p}_{r}x", f"refstage_{r}x") for p in ("refblock", "refhead")
                                       for r in ("1", "4")])
show("jitter vs same cells", [(f"jitter_1x_{f}own", "random_1x") for f in ("0.25", "0.5", "1.0", "2.0")]
     + [(f"jitter_4x_{f}truth", "random_4x") for f in ("0.25", "0.5", "1.0", "2.0")]
     + [(f"jitter_4x_{f}own", "random_4x") for f in ("0.25", "0.5", "1.0", "2.0")])

print("\n== score mean (sd) per condition (key ones)")
for c in ["random_1x", "refstage_1x", "random_2x", "random_8x"]:
    print(c, " | ".join(f"{t}: {S[c].loc[t].mean():.2f} ({S[c].loc[t].std(ddof=1):.2f})" for t in T))

print("\n== nbhd raw / ceiling (mean over seeds), random")
for t in T:
    q = df[df.target == t]
    vals = []
    for r in ("0.25", "0.5", "1", "2", "4", "8"):
        z = q[q.cond == f"random_{r}x"]
        vals.append(f"{r}x {np.mean(z.v_neighborhood_mmd / z.c_neighborhood_mmd):.2f}")
    z = q[q.cond == "random_1x"]
    print(t, "floor/ceil", f"{np.mean(z.f_neighborhood_mmd / z.c_neighborhood_mmd):.0f}",
          f"(ratio of means {z.f_neighborhood_mmd.mean() / z.c_neighborhood_mmd.mean():.0f})", " ".join(vals))
print("per-seed argmin of raw nbhd over random counts:")
for t in T:
    q = df[(df.target == t) & df.cond.str.startswith("random_")]
    am = q.loc[q.groupby("seed").v_neighborhood_mmd.idxmin(), ["seed", "cond"]]
    print(t, dict(zip(am.seed, am.cond, strict=True)))

print("\n== dice voxel/NN ratio (truth side after count matching) random")
for t in T:
    q = df[df.target == t]
    print(t, " ".join(f"{r}x {q[q.cond == f'random_{r}x'].vox_nn.mean():.1f}" for r in ("0.25", "0.5", "1", "2")),
          f"ceil {q[q.cond == 'random_1x'].vox_nn_ceil.mean():.1f}")

print("\n== loss decomposition below 1x (board points = 25 x group-skill change; shape group gated)")
G = ["g_expression_change", "g_cell_state", "g_shape_growth", "g_local_organisation"]
for t in T:
    q = df[df.target == t].set_index(["seed", "cond"])
    for r in ("0.25", "0.5"):
        d = {g: 25 * (q.xs(f"random_{r}x", level="cond")[g] - q.xs("random_1x", level="cond")[g]).mean() for g in G}
        dd = {m: (q.xs(f"random_{r}x", level="cond")["s_" + m] - q.xs("random_1x", level="cond")["s_" + m]).mean()
              for m in M}
        print(t, r, {k[2:]: round(v, 2) for k, v in d.items()}, "skill deltas:",
              {k: round(v, 3) for k, v in dd.items()})

print("\n== anchors (mean over seeds of random_1x rows)")
for t in T:
    z = df[(df.target == t) & (df.cond == "random_1x")]
    print(t, " ".join(f"{m}: {z['f_' + m].mean():.4g}/{z['c_' + m].mean():.4g}" for m in M))
print("\n== jitter rms at f=1 and nn values")
for t in T:
    q = df[df.target == t]
    print(t, "nn_truth", round(q.nn_truth.mean(), 2), "1x own jr", round(q[q.cond == "jitter_1x_1.0own"].jr.mean(), 2),
          "4x own", round(q[q.cond == "jitter_4x_1.0own"].jr.mean(), 2),
          "4x truth", round(q[q.cond == "jitter_4x_1.0truth"].jr.mean(), 2))
