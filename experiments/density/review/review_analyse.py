"""Reviewer's analysis of runs/density_review/raw (and reproduction check against runs/density/raw).

Usage (DGX, repo root): uv run python experiments/density/review/review_analyse.py > runs/density_review/analysis.txt

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
from scipy import stats

MINE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("runs/density_review/raw")
THEIRS = Path("runs/density/raw")
M = ["de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
     "neighborhood_mmd"]


def load(d: Path) -> pd.DataFrame:
    rows = []
    for p in sorted(d.glob("*.jsonl")):
        for line in p.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                o = dict(target=r["target"], seed=r["seed"], cond=r["condition"], n=r["n_pred"], score=r["score"],
                         tc=r["truth_cells"], frac=r.get("frac", 0.1))
                for m in M:
                    o["s_" + m] = r["skills"][m]
                    o["v_" + m] = r["metrics"][m]
                    o["c_" + m] = r["ceiling"][m]
                    o["f_" + m] = r["floor"][m]
                o["vox_nn"] = r["metrics"].get("_dice_voxel_over_nn")
                o.update({"g_" + k: v for k, v in r["group_skill"].items()})
                rows.append(o)
    return pd.DataFrame(rows)


mine = load(MINE)
print(f"reviewer rows: {len(mine)}; unique: {len(mine.drop_duplicates(['target', 'seed', 'cond']))}")
if THEIRS.exists() and MINE == Path("runs/density_review/raw"):
    theirs = load(THEIRS)
    j = mine.merge(theirs, on=["target", "seed", "cond"], suffixes=("_m", "_t"))
    sk = np.max(np.abs(np.stack([j[f"s_{m}_m"].fillna(-9) - j[f"s_{m}_t"].fillna(-9) for m in M])), axis=0)
    print(f"\n== reproduction: {len(j)} overlapping (target, seed, condition) rows; "
          f"max |score diff| = {np.abs(j.score_m - j.score_t).max():.3g}; max |skill diff| = {sk.max():.3g}; "
          f"n_pred equal: {bool((j.n_m == j.n_t).all())}")
    print("conditions compared:", sorted(j.cond.unique()))
    bad = j[np.abs(j.score_m - j.score_t) > 1e-9]
    if len(bad):
        print(bad[["target", "seed", "cond", "score_m", "score_t"]].to_string())

S = mine.pivot_table(index=["target", "seed"], columns="cond", values="score")
TG = [t for t in ["heart_E9.5", "heart_E8.75", "embryo_E8.0", "embryo_E7.25"] if t in set(mine.target)]


def ci(d):
    d = d.dropna()
    n = len(d)
    m = d.mean()
    se = d.std(ddof=1) / np.sqrt(n)
    h = stats.t.ppf(0.975, n - 1) * se
    return m, h, int((d < 0).sum()), n


def table(pairs, title):
    print(f"\n== {title}: paired mean, 95% CI half-width [seeds<0 / n]")
    for a, b in pairs:
        if a not in S or b not in S:
            continue
        cells = []
        for t in TG:
            m, h, k, n = ci(S[a].loc[t] - S[b].loc[t])
            cells.append(f"{t}: {m:+.2f} ±{h:.2f} [{k}/{n}]")
        d_all = pd.concat([S[a].loc[t] - S[b].loc[t] for t in TG if t != "embryo_E7.25"])
        m, h, k, n = ci(d_all)
        print(f"{a} - {b}: " + " | ".join(cells) + f" || 3 orig targets pooled: {m:+.2f} ±{h:.2f} [{k}/{n}]")


rs = ["0.5", "2", "4", "8", "3"]
table([(f"random_{r}x", "random_1x") for r in rs], "near-perfect (random) vs 1x")
table([(f"refstage_{r}x", "refstage_1x") for r in rs], "floor-quality (refstage) vs 1x")
table([(f"shuffle50_{r}x", "shuffle50_1x") for r in rs], "intermediate (50% expression shuffled) vs 1x")
table([("random_2x", "random_4x"), ("random_2x", "random_8x"), ("refstage_2x", "refstage_4x"),
       ("refstage_2x", "refstage_8x"), ("shuffle50_2x", "shuffle50_4x"), ("shuffle50_2x", "shuffle50_8x"),
       ("random_2x", "random_3x"), ("refstage_2x", "refstage_3x"), ("shuffle50_2x", "shuffle50_3x")],
      "2x vs higher counts")
table([("jitter_1x_2.0own", "random_1x"), ("jitter_4x_2.0truth", "random_4x")], "jitter at 2 NN (same cells)")

print("\n== does 4x tolerate the same physical jitter better? loss(1x/own f=2) - loss(4x/truth f=2), paired by seed")
for t in TG:
    if "jitter_1x_2.0own" not in S:
        break
    l1 = S["jitter_1x_2.0own"].loc[t] - S["random_1x"].loc[t]
    l4 = S["jitter_4x_2.0truth"].loc[t] - S["random_4x"].loc[t]
    m, h, k, n = ci(l4 - l1)
    print(f"{t}: loss1x {l1.mean():+.2f}, loss4x {l4.mean():+.2f}, (loss4x - loss1x) {m:+.2f} ±{h:.2f} "
          f"[4x worse in {int((l4 < l1).sum())}/{n}]")

print("\n== mean score (sd) and mean skills per condition")
for t in TG:
    q = mine[mine.target == t]
    print(f"-- {t} (truth_cells {int(q.tc.iloc[0])}, frac {q.frac.iloc[0]})")
    for c in sorted(q.cond.unique(), key=lambda c: (c.split("_")[0], float(c.split("_")[1].rstrip("x"))
                                                    if c.split("_")[1].rstrip("x").replace(".", "").isdigit()
                                                    else 0, c)):
        z = q[q.cond == c]
        sk = " ".join(f"{m[:5]}={z['s_' + m].mean():.3f}" for m in M)
        print(f"{c:22s} n={int(z.n.mean()):6d} seeds={len(z):2d} score={z.score.mean():6.2f} "
              f"({z.score.std(ddof=1):.2f}) "
              f"nbhd/ceil={np.mean(z.v_neighborhood_mmd / z.c_neighborhood_mmd):5.2f} "
              f"vox/nn={z.vox_nn.mean():.1f} | {sk}")

print("\n== per-target argmin of mean score loss and of raw nbhd over counts (random / refstage / shuffle50)")
for t in TG:
    for pat in ("random", "refstage", "shuffle50"):
        cs = [c for c in S.columns if c.startswith(pat + "_") and c.endswith("x")]
        if not cs:
            continue
        means = {c: S[c].loc[t].mean() for c in cs}
        best = max(means, key=means.get)
        q = mine[(mine.target == t) & (mine.cond.isin(cs))]
        nb = q.groupby("cond").apply(lambda z: np.mean(z.v_neighborhood_mmd / z.c_neighborhood_mmd))
        print(f"{t} {pat}: best mean score {best} ({means[best]:.2f}); raw-nbhd argmin {nb.idxmin()}; "
              + ", ".join(f"{c.split('_')[1]} {means[c]:.2f}" for c in sorted(cs, key=lambda c: float(
                  c.split('_')[1].rstrip('x')))))
