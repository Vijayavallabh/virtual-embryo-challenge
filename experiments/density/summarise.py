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

Decide mode (--decide), for the release rerun (docs/P3_RUNBOOK.md step 4.5): applies approved rule 7 of
docs/RELEASE_RULES_SUMMARY.md (2026-10-21) by code, vec.t2_build.count_override, to the raw rows' intermediate-quality
condition (shuffle50) at the max_cells levels: heart_E10.5 decides the heart test board, embryo_E7.5 the embryo test
board, each separately (vec.t2_build.COUNT_RULE_TARGETS). Rows that do not measure the rule's levels on the validation
board (no single max_cells from the board's published limit, a level capped at the pool) keep max_cells, as an
incomplete rerun does. It writes only <out>/count_decision.json and prints each board's decision; the summary above is
left to the plain run. --raw and --out are required here (the defaults are the pre-release run's). Refused, with
nothing written: uncommitted changes to the rule code (src/vec/t2_build.py, this file) unless --allow-uncommitted
(recorded); without --overwrite, an existing decision in --out or the release record
experiments/density/results/release/count_decision.json, whatever --out says (another --out cannot sidestep the first
decision; --overwrite keeps a replaced decision next to the new one, leaves the release record outside --out where it
is, and the new decision lists them under earlier_decisions); a raw line that is not a JSON row; raw rows with no rule
row for either target. vec.t2_build.n_cells_from_decision reads the file and passes a board's "override" to
vec.t2_build.n_cells on that test board (null keeps max_cells); a level below that board's lower limit gives way to the
passing level ranked next in rule 7's order (largest mean gain, ties to the higher count) that the board allows, else
max_cells (vec.t2_build.COUNT_FLOOR_RULE: the user approved "the next passing count the board allows, or max_cells if
none does" on 2026-10-21; reading "next" in rule 7's order is the AI assistant's and awaits the user's confirmation).
  uv run python experiments/density/summarise.py --raw runs/density/raw_release --out \
    experiments/density/results/release --decide
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import time
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


# --- decide mode -----------------------------------------------------------------------------------------------------
DECISION_FILE = "count_decision.json"
RELEASE_OUT = "experiments/density/results/release"  # the release decision's --out (the command above), from the repo
RULE_FILES = ("src/vec/t2_build.py", "experiments/density/summarise.py")  # the code that applies approved rule 7


def release_record() -> Path:
    """The release decision's canonical record, <repo>/experiments/density/results/release/count_decision.json (the
    --out of the release command in the module docstring). `decide` refuses while it exists whatever --out says,
    unless overwrite: another --out cannot sidestep the first decision (the guard of experiments/t1_direction/
    readout.py)."""
    from vec.boards import REPO_ROOT

    return REPO_ROOT / RELEASE_OUT / DECISION_FILE


def _sorted(values) -> list:
    """Numbers in numeric order, then anything else (None, text) in repr order: for records and messages."""
    def key(v):
        number = isinstance(v, int | float) and not isinstance(v, bool)
        return (0, v, "") if number else (1, 0, repr(v))
    return sorted(values, key=key)


def rule_rows(raw: Path, condition: str, targets) -> tuple[dict, dict, int]:
    """The raw rows approved rule 7 reads: per target, {(seed, condition label): row} for `condition` at the max_cells
    levels (basis "max"), from every *.jsonl in `raw` in name order; a later row replaces an earlier one with the same
    target, seed and condition, as `load` keeps the last. Also the sha256 of each file read and the rows replaced. A
    line that is not a JSON object, such as the partial last line a killed job leaves, is refused with its file and
    line number."""
    rows: dict[str, dict] = {t: {} for t in targets}
    sha, replaced = {}, 0
    for p in sorted(raw.glob("*.jsonl")):
        data = p.read_bytes()
        sha[p.name] = hashlib.sha256(data).hexdigest()
        for i, line in enumerate(data.decode(errors="replace").splitlines(), 1):
            if not line.strip():
                continue
            try:
                r, why = json.loads(line), "not a JSON object"
            except json.JSONDecodeError as e:
                r, why = None, str(e)
            if not isinstance(r, dict):
                raise SystemExit(f"REFUSED: {p} line {i} is not a raw row ({why}). A killed job leaves a partial last "
                                 "line: delete it, rerun that target and seed (finished conditions are skipped), then "
                                 "decide. Nothing written")
            if r.get("target") in rows and r.get("pattern") == condition and r.get("basis") == "max":
                key = (r.get("seed"), r.get("condition"))
                replaced += key in rows[r["target"]]
                rows[r["target"]][key] = r
    return rows, sha, replaced


def row_problems(by_level: dict[str, list[dict]], levels) -> list[str]:
    """Why a validation target's rule rows (by level label) do not measure approved rule 7's levels on the validation
    board; empty when they do. The rows must agree on one max_cells, the validation board's published limit
    (max_cells_source "board", not the documented cap rule or --max-cells), and each of `levels` must have drawn
    round(f x max_cells) cells. run_density caps a level at the pool, and a capped level is not the rule's: a 1max
    short of max_cells is not max_cells, and a lower level capped like 1max is not a lower count."""
    rows = [r for rs in by_level.values() for r in rs]
    if not rows:
        return []
    out = []
    caps = _sorted({r.get("max_cells") for r in rows})
    max_cells = caps[0] if len(caps) == 1 else None
    if isinstance(max_cells, bool) or not isinstance(max_cells, int) or max_cells < 1:
        out.append(f"the rows disagree on max_cells or lack it ({caps}): their levels are not one board's")
        max_cells = None
    sources = sorted({str(r.get("max_cells_source")) for r in rows})
    if sources != ["board"]:
        out.append(f"max_cells came from {sources}, not only from the validation board's published limit")
    if max_cells is not None:
        for level in levels:
            f = float(level.removesuffix("max"))
            want = int(round(f * max_cells))
            got = _sorted({r.get("n_pred") for r in by_level.get(level, [])})
            if got and got != [want]:
                pools = _sorted({r.get("pool_cells") for r in by_level[level]})
                out.append(f"{level} drew {got} cells, not round({f:g} x {max_cells}) = {want} (pool {pools}): not "
                           "the rule's level")
    return out


def rule_code_state(repo: Path) -> dict:
    """The state of the code that applies approved rule 7: the git HEAD, the paths that keep the tree from being clean
    (vec.provenance.git_state), the uncommitted changes to RULE_FILES among them (None when git cannot tell) and the
    sha256 of each of RULE_FILES."""
    from vec.provenance import git_state

    commit, unclean = git_state(repo)
    return {"git_commit": commit, "git_unclean": unclean,
            "rule_files_uncommitted": None if unclean is None else [p for p in unclean if p in RULE_FILES],
            "rule_files_sha256": {p: hashlib.sha256((repo / p).read_bytes()).hexdigest() for p in RULE_FILES
                                  if (repo / p).is_file()}}


def _replaced_name(path: Path) -> Path:
    """A free name next to `path` for the record that --overwrite replaces: <stem>.replaced-<UTC time>[-k]<suffix>."""
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    for k in itertools.count():
        p = path.with_name(f"{path.stem}.replaced-{stamp}{f'-{k}' if k else ''}{path.suffix}")
        if not p.exists():
            return p
    raise AssertionError("unreachable")


def decide(raw: Path, out: Path, overwrite: bool = False, allow_uncommitted: bool = False) -> dict:
    """Decide mode: approved rule 7, applied by vec.t2_build.count_override, to each test board's validation target
    (vec.t2_build.COUNT_RULE_TARGETS) in the rerun's raw rows. Writes <out>/count_decision.json (strict JSON, for
    committing; vec.t2_build.n_cells_from_decision reads it) and returns its content: per board the override for
    vec.t2_build.n_cells, the decision and its record, what the rows say about the counts (max_cells and its source,
    cells drawn per level, row problems), plus the inputs' sha256, the git state and the rule code's sha256.

    The rule keeps max_cells "if the rerun cannot run": a missing level, seed or 1max row keeps it (count_override), and
    so do rows that do not measure the rule's levels on the validation board (`row_problems`). Refused, with nothing
    written: a missing raw directory; an existing decision without `overwrite`, in `out` or at the release record
    (`release_record`, whatever `out` says: the first decision is the record, and another --out cannot sidestep it);
    uncommitted changes to RULE_FILES, unless `allow_uncommitted` (recorded); a raw line that is not a JSON row; and no
    rule row for any validation target (the wrong directory, or the rerun has not run: max_cells then stands with no
    decision file). `overwrite` keeps a replaced decision in `out` next to the new one and leaves the release record
    outside `out` where it is; the new decision lists every decision it replaces or joins (`earlier_decisions`: path
    and sha256, and where a replaced one was kept). A target without rows is decided as an incomplete rerun, and the
    printout flags it."""
    from vec import t2_build as tb
    from vec.boards import REPO_ROOT

    path = out / DECISION_FILE
    record = release_record()
    if not raw.is_dir():
        raise SystemExit(f"REFUSED: no raw directory {raw}")
    if path.exists() and not overwrite:
        raise SystemExit(f"REFUSED: {path} exists (the first decision is the record); pass --overwrite to replace it "
                         "(the replaced record is kept next to it)")
    joined = record if record.exists() and record.resolve() != path.resolve() else None  # the release record elsewhere
    if joined is not None and not overwrite:
        raise SystemExit(f"REFUSED: the release decision {record} exists (the first decision is the record, and "
                         "another --out does not sidestep it). Pass --overwrite to decide again: with --out "
                         f"{record.parent} the replaced record is kept next to the new one; with another --out it "
                         "stays where it is and the new decision names it")
    code = rule_code_state(REPO_ROOT)
    dirty = code["rule_files_uncommitted"]
    if (dirty is None or dirty) and not allow_uncommitted:
        what = "git cannot tell whether they are committed" if dirty is None else "\n  ".join(dirty)
        raise SystemExit(f"REFUSED: uncommitted changes to the files that apply approved rule 7 (pre-registration):\n"
                         f"  {what}\nCommit them first, or pass --allow-uncommitted (recorded)")
    targets = tb.COUNT_RULE_TARGETS
    rows, sha, replaced = rule_rows(raw, tb.COUNT_RULE_CONDITION, targets.values())
    if not any(rows.values()):
        raise SystemExit(f"REFUSED: no {tb.COUNT_RULE_CONDITION} rows at the max_cells levels for "
                         f"{' or '.join(targets.values())} in {raw}: the wrong directory, or the rerun has not run. "
                         "Nothing written; if the rerun cannot run, approved rule 7 keeps max_cells (n_cells with no "
                         "override)")
    levels = (tb.COUNT_RULE_BASE, *tb.COUNT_RULE_LEVELS)
    boards = {}
    for board, target in targets.items():
        scores: dict[str, dict] = {}
        by_level: dict[str, list[dict]] = {}
        for (seed, cond), r in rows[target].items():
            level = str(cond).removeprefix(f"{tb.COUNT_RULE_CONDITION}_")
            scores.setdefault(level, {})[seed] = r.get("score")
            by_level.setdefault(level, []).append(r)
        override, record = tb.count_override(scores)
        problems = row_problems(by_level, levels)
        decision = record["decision"]
        if problems:
            decision = (f"keep max_cells: the rows do not measure approved rule 7's levels on the validation board "
                        f"({'; '.join(problems)}), and the rule keeps max_cells when the rerun cannot run")
            if override is not None:
                decision += f" (count_override on these rows alone gives {override!r})"
        sel = list(rows[target].values())
        boards[board] = {
            "validation_target": target, "override": None if problems else override, "decision": decision,
            "rows": len(sel), "row_problems": problems, "max_cells": _sorted({r.get("max_cells") for r in sel}),
            "max_cells_source": sorted({str(r.get("max_cells_source")) for r in sel}),
            "n_pred": {lv: _sorted({r.get("n_pred") for r in rs}) for lv, rs in sorted(by_level.items())},
            "record": record}
    earlier, kept = [], None  # the decisions this one replaces (overwrite in `out`) or joins (the release record)
    if path.exists():
        kept = _replaced_name(path)
        earlier.append({"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "kept_as": str(kept)})
    if joined is not None:
        earlier.append({"path": str(joined), "sha256": hashlib.sha256(joined.read_bytes()).hexdigest(),
                        "kept_as": str(joined)})
    result = {
        "what": ("the cell count of each T2 test board under approved rule 7 (docs/RELEASE_RULES_SUMMARY.md), applied "
                 "by vec.t2_build.count_override to the density rerun's raw rows; vec.t2_build.n_cells_from_decision "
                 "passes a board's override to vec.t2_build.n_cells on that test board (null keeps max_cells; a "
                 "level below that board's lower limit gives way to the passing level ranked next in rule 7's order "
                 "that the board allows, else max_cells: vec.t2_build.COUNT_FLOOR_RULE)"),
        "rule": tb.COUNT_RULE, "approved": tb.RULES_APPROVED, "raw": str(raw), "inputs_sha256": sha,
        "rows_replaced": replaced, "earlier_decisions": earlier, "allow_uncommitted": bool(allow_uncommitted), **code,
        "boards": boards,
    }
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    out.mkdir(parents=True, exist_ok=True)
    if kept is not None:  # --overwrite keeps the replaced record next to the new one
        path.rename(kept)
    path.write_text(text)
    print(f"count rule (approved {tb.RULES_APPROVED}) on {raw}:")
    for board, b in boards.items():
        print(f"  {board} board ({b['validation_target']}): override {b['override']!r}; {b['decision']}")
        if not b["rows"]:
            print(f"    WARNING: no rows for {b['validation_target']} in {raw}: check the directory; if its rerun "
                  "did not run, the rule keeps max_cells")
    if dirty is None or dirty:
        print(f"  WARNING (allowed, recorded): uncommitted rule code: {dirty if dirty else 'git cannot tell'}")
    if kept is not None:
        print(f"  replaced decision kept as {kept}")
    if joined is not None:
        print(f"  WARNING (--overwrite): the release decision {joined} stays where it is, and this one is a second "
              "decision beside it (named in earlier_decisions); vec.t2_build.n_cells_from_decision reads whichever "
              "file it is given")
    print(f"wrote {path}")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=None, help="raw rows directory [runs/density/raw_v2]; required with --decide")
    ap.add_argument("--out", default=None, help="output directory [runs/density/summary_v2]; required with --decide")
    ap.add_argument("--decide", action="store_true",
                    help=f"decide mode: apply approved rule 7 to the raw rows and write only <out>/{DECISION_FILE}")
    ap.add_argument("--overwrite", action="store_true",
                    help=f"with --decide: decide again although a {DECISION_FILE} exists in --out (the replaced one is "
                         f"kept next to it) or in {RELEASE_OUT} (left where it is when --out is elsewhere)")
    ap.add_argument("--allow-uncommitted", action="store_true",
                    help="with --decide: decide although the rule code has uncommitted changes (recorded)")
    args = ap.parse_args()
    if args.decide:
        missing = [opt for opt, v in (("--raw", args.raw), ("--out", args.out)) if v is None]
        if missing:
            ap.error(f"--decide needs {' and '.join(missing)}, e.g. --raw runs/density/raw_release --out "
                     f"{RELEASE_OUT} (the defaults are the pre-release run's)")
        decide(Path(args.raw), Path(args.out), overwrite=args.overwrite, allow_uncommitted=args.allow_uncommitted)
        return
    out = Path(args.out or "runs/density/summary_v2")
    out.mkdir(parents=True, exist_ok=True)
    df = load(Path(args.raw or "runs/density/raw_v2"))
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
