"""Summarise the density experiment: per-condition mean/sd over seeds, paired deltas, plots.

Usage (DGX): uv run python experiments/density/summarise.py --raw runs/density/raw_v2 --out runs/density/summary_v2
Outputs: results_long.csv (one row per target x seed x condition), summary_by_condition.csv (mean/sd per target x
condition plus paired deltas), summary_pooled.csv (target means averaged), PNG plots.

Paired deltas (same target and seed, so the same working copy and truth half):
  d_*   against the 1x condition of the same family: random_1x for target cells (random, block, head, jitter),
        shuffle50_1x for shuffle50, refstage_1x for reference cells (refstage, refblock, refhead)
  dr_*  against the random subset of the same source at the same count level: random for target cells (for jitter
        and shuffle50 these are the very same cells, unperturbed), refstage for reference cells
  dm_*  against the same pattern at max_cells (the 1max level), for random, shuffle50 and refstage
Rows written before the max_cells levels existed have no basis/level fields; they are read as basis "truth" with
level = ratio.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vec.boards import ABS_METRICS, LOWER_IS_BETTER, RANKED_METRICS  # noqa: E402

METRICS = RANKED_METRICS["T2"]
SHORT = {"de_score": "de", "de_direction": "dir", "mmd_u": "mmd", "variogram": "vario", "d2_shape": "d2",
         "occupancy_dice": "dice", "scale_log_ratio": "scale", "neighborhood_mmd": "nbhd"}
TARGET_ORDER = ["heart_E9.5", "heart_E8.75", "embryo_E8.0", "heart_E10.5", "embryo_E7.5"]
PATTERN_RANK = {"random": 0, "shuffle50": 1, "block": 2, "head": 3, "jitter_own": 4, "jitter_truth": 5,
                "refstage": 6, "refblock": 7, "refhead": 8}
ROW_KEYS = ("target", "board", "seed", "condition", "source", "pattern", "basis", "level", "ratio", "jitter_f",
            "jitter_unit", "block_start", "block_draw", "n_pred", "truth_cells", "ref_cells", "max_cells",
            "max_cells_source", "pool_cells", "n_ratio", "nn_own", "nn_truth", "jitter_rms", "score", "seconds")


def gap_fraction(m, v, f, c):
    """x = signed distance from the ceiling in floor-ceiling units (skill = 1/(1+x), clipped at 1 for x < 0)."""
    if v is None or f is None or c is None:
        return np.nan
    if m in ABS_METRICS:
        v, f, c = abs(v), abs(f), abs(c)
    d_floor = abs(c - f)
    if d_floor < 1e-12:
        return np.nan
    return ((v - c) if m in LOWER_IS_BETTER else (c - v)) / d_floor


def load(raw: Path) -> pd.DataFrame:
    rows = []
    for p in sorted(raw.glob("*.jsonl")):
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            out = {k: r.get(k) for k in ROW_KEYS}
            out.update({f"group_{k}": v for k, v in r["group_skill"].items()})
            for m in METRICS:
                v, f, c = r["metrics"].get(m), r["floor"].get(m), r["ceiling"].get(m)
                out[f"skill_{m}"] = r["skills"].get(m)
                out[f"val_{m}"], out[f"floor_{m}"], out[f"ceil_{m}"] = v, f, c
                out[f"x_{m}"] = gap_fraction(m, v, f, c)
            for m, v in r["metrics"].items():
                if m not in METRICS:
                    out[f"val_{m}"] = v
            for m, v in r.get("ceiling_extra", {}).items():
                out[f"ceil_{m}"] = v
            out["nbhd_over_ceiling"] = (r["metrics"]["neighborhood_mmd"] / r["ceiling"]["neighborhood_mmd"]
                                        if r["ceiling"].get("neighborhood_mmd") else np.nan)
            out["inverted_anchors"] = ",".join(r.get("inverted_anchors", []))
            rows.append(out)
    df = pd.DataFrame(rows)
    df["jitter_f"] = df["jitter_f"].fillna(0.0)
    df["jitter_unit"] = df["jitter_unit"].fillna("")
    df["basis"] = df["basis"].fillna("truth")
    df["level"] = pd.to_numeric(df["level"]).fillna(pd.to_numeric(df["ratio"]))
    for c in ("ratio", "block_start", "max_cells"):
        df[c] = pd.to_numeric(df[c])
    n0 = len(df)
    df = df.drop_duplicates(["target", "seed", "condition"], keep="last").reset_index(drop=True)
    if len(df) < n0:
        print(f"WARNING: dropped {n0 - len(df)} duplicate (target, seed, condition) rows, keeping the last")
    return df


def condition_order(df: pd.DataFrame) -> list[str]:
    keys = df[["condition", "pattern", "basis", "level", "jitter_f"]].drop_duplicates("condition")
    keys = keys.assign(pr=keys["pattern"].map(PATTERN_RANK).fillna(99), br=(keys["basis"] == "max").astype(int))
    return list(keys.sort_values(["pr", "br", "level", "jitter_f"])["condition"])


def summarise(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    value_cols = (["score", "n_pred", "n_ratio", "pool_cells", "nn_own", "jitter_rms", "nbhd_over_ceiling"]
                  + [f"group_{g}" for g in ("expression_change", "cell_state", "shape_growth", "local_organisation")]
                  + [f"skill_{m}" for m in METRICS] + [f"x_{m}" for m in METRICS] + [f"val_{m}" for m in METRICS]
                  + [f"ceil_{m}" for m in METRICS] + [f"floor_{m}" for m in METRICS]
                  + ["val__dice_voxel_over_nn", "val_count_log_ratio", "val_energy_distance",
                     "val_morans_I_agreement"])
    value_cols = [c for c in value_cols if c in df]
    delta_cols = ["score"] + [f"skill_{m}" for m in METRICS] + [f"x_{m}" for m in METRICS]
    by_cond = df.set_index(["target", "seed", "condition"])[delta_cols]
    d = df.copy()
    d["base_1x"] = np.select([d["pattern"] == "shuffle50", d["source"] == "ref"], ["shuffle50_1x", "refstage_1x"],
                             "random_1x")
    d = d.join(by_cond, on=["target", "seed", "base_1x"], rsuffix="_ref1x")
    d["base_pattern"] = np.where(d["source"] == "ref", "refstage", "random")
    rnd = df[df["pattern"].isin(["random", "refstage"])].set_index(["target", "seed", "pattern", "basis", "level"])
    d = d.join(rnd[delta_cols], on=["target", "seed", "base_pattern", "basis", "level"], rsuffix="_samen")
    d["base_max"] = np.where(d["pattern"].isin(["random", "shuffle50", "refstage"]), d["pattern"] + "_1max", "")
    d = d.join(by_cond, on=["target", "seed", "base_max"], rsuffix="_max")
    for c in delta_cols:
        d[f"d_{c}"] = d[c] - d[f"{c}_ref1x"]
        d[f"dr_{c}"] = d[c] - d[f"{c}_samen"]
        d[f"dm_{c}"] = d[c] - d[f"{c}_max"]
    dcols = [f"{p}_{c}" for p in ("d", "dr", "dm") for c in delta_cols]
    g = d.groupby(["target", "condition"])
    mean = g[value_cols + dcols].mean()
    sd = g[value_cols + dcols].std(ddof=1)
    n = g.size().rename("n_seeds")
    neg = [g[f"{p}_score"].apply(lambda s: int((s < 0).sum())).rename(f"seeds_{p}score_neg") for p in ("d", "dr", "dm")]
    meta = g[["pattern", "basis", "level", "ratio", "jitter_f", "jitter_unit", "truth_cells", "max_cells",
              "max_cells_source"]].first()
    bstart = g["block_start"].agg(["min", "max"]).add_prefix("block_start_")
    summ = pd.concat([meta, n, *neg, bstart, mean.add_suffix("_mean"), sd.add_suffix("_sd")], axis=1).reset_index()
    summ["condition"] = pd.Categorical(summ["condition"], condition_order(df), ordered=True)
    targets = [t for t in TARGET_ORDER if t in set(df["target"])]
    targets += sorted(set(df["target"]) - set(targets))
    summ["target"] = pd.Categorical(summ["target"], targets, ordered=True)
    summ = summ.sort_values(["target", "condition"])
    pooled_cols = [c for c in summ.columns if c.endswith("_mean")]
    first = summ.groupby("condition", observed=True)[["pattern", "basis", "level", "ratio", "jitter_f",
                                                      "jitter_unit"]].first()
    pooled = summ.groupby("condition", observed=True)[pooled_cols].mean()
    pooled_sd = summ.groupby("condition", observed=True)[pooled_cols].std(ddof=1).add_suffix("_across_targets_sd")
    keep_sd = [c for c in pooled_sd.columns if c.startswith(("score", "d_score", "dr_score", "dm_score"))]
    n_t = summ.groupby("condition", observed=True).size().rename("n_targets")
    pooled = pd.concat([first, n_t, pooled, pooled_sd[keep_sd]], axis=1)
    return summ, pooled.reset_index()


def print_tables(summ: pd.DataFrame):
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.max_rows", 500)
    for t, s in summ.groupby("target", observed=True):
        mx = s["max_cells"].dropna()
        mx_txt = f", max_cells={int(mx.iloc[0])} ({s['max_cells_source'].dropna().iloc[0]})" if len(mx) else ""
        print(f"\n=== {t} (truth_cells={int(s['truth_cells'].iloc[0])}{mx_txt}) : mean (sd) over seeds ===")
        tab = pd.DataFrame({"n": s["n_pred_mean"].round(0).astype(int).values,
                            "seeds": s["n_seeds"].values}, index=s["condition"].values)
        tab["score"] = [f"{m:.2f} ({v:.2f})" for m, v in zip(s["score_mean"], s["score_sd"], strict=True)]
        for p, lab in (("d", "dFam1x"), ("dr", "dSameN"), ("dm", "dMax")):
            if s[f"{p}_score_mean"].notna().any():
                tab[lab] = [f"{m:+.2f} ({v:.2f})" if np.isfinite(m) else "" for m, v in
                            zip(s[f"{p}_score_mean"], s[f"{p}_score_sd"], strict=True)]
        for m in METRICS:
            tab[SHORT[m]] = [f"{a:.3f}({b:.3f})" for a, b in zip(s[f"skill_{m}_mean"], s[f"skill_{m}_sd"],
                                                                strict=True)]
        print(tab.to_string())
        print(f"--- {t}: gap fraction x (mean; <0 = better than the split-half ceiling) and raw nbhd/ceiling ---")
        tab2 = pd.DataFrame(index=s["condition"].values)
        for m in METRICS:
            tab2[SHORT[m]] = s[f"x_{m}_mean"].round(3).values
        tab2["nbhd/ceil"] = s["nbhd_over_ceiling_mean"].round(3).values
        tab2["nn_own"] = s["nn_own_mean"].round(2).values
        tab2["jit_rms"] = s["jitter_rms_mean"].round(2).values
        print(tab2.to_string())


def plots(df: pd.DataFrame, summ: pd.DataFrame, out: Path):
    targets = list(summ["target"].cat.categories)
    colors = {"random": "#1f77b4", "shuffle50": "#2ca02c", "block": "#d62728", "head": "#ff7f0e",
              "refstage": "#7f7f7f", "refblock": "#9467bd", "refhead": "#8c564b"}
    truth_ratios = sorted(df.loc[df["basis"] == "truth", "ratio"].dropna().unique())
    fig, axes = plt.subplots(1, len(targets), figsize=(5 * len(targets), 4.2), sharey=True, squeeze=False)
    for ax, t in zip(axes[0], targets, strict=True):
        s = summ[summ["target"] == t]
        for pat, col in colors.items():
            q = s[(s["pattern"] == pat) & (s["basis"] == "truth")].sort_values("ratio")
            if len(q):
                ax.errorbar(q["ratio"], q["score_mean"], yerr=q["score_sd"], marker="o", color=col, label=pat,
                            capsize=3, ls="--" if pat.startswith("ref") else "-")
            q = s[(s["pattern"] == pat) & (s["basis"] == "max")].sort_values("ratio")
            if len(q):
                ax.errorbar(q["ratio"], q["score_mean"], yerr=q["score_sd"], marker="s", mfc="none", color=col,
                            capsize=2, ls=":")
        mx = s["max_cells"].dropna()
        if len(mx):
            ax.axvline(float(mx.iloc[0]) / float(s["truth_cells"].iloc[0]), color="k", lw=0.8, ls="-.")
        ax.set_xscale("log", base=2)
        ax.set_xticks(truth_ratios)
        ax.set_xticklabels([f"{r:g}" for r in truth_ratios])
        ax.set_title(t)
        ax.set_xlabel("submitted cells / truth_cells")
        ax.axhline(50, color="k", lw=0.5, ls=":")
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel(f"emulator board score (mean ± sd, {int(summ['n_seeds'].max())} seeds)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if (summ["basis"] == "max").any():
        handles.append(plt.Line2D([], [], color="k", marker="s", mfc="none", ls=":"))
        labels.append("levels at fractions of max_cells (dash-dot line = max_cells)")
    fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 5), frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(out / "score_vs_count.png", dpi=130)
    plt.close(fig)

    rnd = summ[(summ["pattern"] == "random") & (summ["basis"] == "truth")]
    fig, axes = plt.subplots(2, 4, figsize=(18, 7.5))
    for ax, m in zip(axes.ravel(), METRICS, strict=True):
        for t, mk in zip(targets, ["o", "s", "^", "D", "v"], strict=False):
            q = rnd[rnd["target"] == t].sort_values("ratio")
            ax.errorbar(q["ratio"], q[f"x_{m}_mean"], yerr=q[f"x_{m}_sd"], marker=mk, capsize=3, label=t)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log", base=2)
        ax.set_title(f"{m}: gap fraction x (random)")
        ax.set_xlabel("cells / truth_cells")
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "gap_fraction_vs_count_random.png", dpi=130)
    plt.close(fig)

    jit = summ[summ["pattern"].isin(["jitter_own", "jitter_truth"])]
    if not len(jit):
        return
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
    for ax, (col, lab) in zip(axes, [("d_score", "Δ board score vs random 1x"),
                                     ("d_skill_neighborhood_mmd", "Δ neighborhood_mmd skill"),
                                     ("d_skill_occupancy_dice", "Δ occupancy_dice skill")], strict=True):
        for t, mk in zip(targets, ["o", "s", "^", "D", "v"], strict=False):
            for (r, unit), ls in zip([(1, "own"), (4, "own"), (4, "truth")], ["-", "--", ":"], strict=True):
                q = jit[(jit["target"] == t) & (jit["ratio"] == r) & (jit["jitter_unit"] == unit)]
                q = q.sort_values("jitter_f")
                base = summ[(summ["target"] == t) & (summ["condition"] == f"random_{r}x")]
                xs = [0.0] + list(q["jitter_f"])
                ys = list(base[f"{col}_mean"]) + list(q[f"{col}_mean"])
                ax.plot(xs, ys, marker=mk, ls=ls, label=f"{t} {r}x/{unit}")
        ax.set_xlabel("jitter RMS / median NN distance")
        ax.set_title(lab)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "jitter.png", dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="runs/density/raw_v2")
    ap.add_argument("--out", default="runs/density/summary_v2")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df = load(Path(args.raw))
    df.to_csv(out / "results_long.csv", index=False, float_format="%.6g")
    summ, pooled = summarise(df)
    summ.to_csv(out / "summary_by_condition.csv", index=False, float_format="%.6g")
    pooled.to_csv(out / "summary_pooled.csv", index=False, float_format="%.6g")
    print(f"{len(df)} rows; conditions per target: {df.groupby('target').size().to_dict()}")
    bad = df[df["inverted_anchors"] != ""]
    print(f"rows with inverted anchors: {len(bad)}")
    print_tables(summ)
    plots(df, summ, out)


if __name__ == "__main__":
    main()
