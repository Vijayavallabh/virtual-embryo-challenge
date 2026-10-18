"""Print the markdown tables used in the density report from summary_by_condition.csv (stdlib only, runs anywhere).
Conditions missing from the CSV are skipped, so the same script serves the first run, v2 and the release rerun.

Usage: python3 experiments/density/report_tables.py [results_dir]
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

RES = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results"
ORDER = ["heart_E9.5", "heart_E8.75", "embryo_E8.0", "heart_E10.5", "embryo_E7.5"]
REFS = {"heart_E9.5": "E8.75", "heart_E8.75": "E8.25_late", "embryo_E8.0": "E7.25", "heart_E10.5": "E9.5",
        "embryo_E7.5": "E7.25"}
RATIOS = ["0.25", "0.5", "1", "2", "4", "8"]
MAX_LEVELS = ["0.25", "0.5", "0.75", "1"]
METRICS = ["de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd"]
SHORT = ["de", "dir", "mmd_u", "vario", "d2", "dice", "scale", "nbhd"]

with open(RES / "summary_by_condition.csv") as fh:
    rows = {(r["target"], r["condition"]): r for r in csv.DictReader(fh)}
present = {t for t, _ in rows}
TARGETS = [t for t in ORDER if t in present] + sorted(present - set(ORDER))


def has(t, c):
    return (t, c) in rows


def g(t, c, k):
    v = rows.get((t, c), {}).get(k, "")
    return float(v) if v not in ("", "nan") else float("nan")


def ms(t, c, k, nd=2, sign=False):
    m, s = g(t, c, f"{k}_mean"), g(t, c, f"{k}_sd")
    return f"{m:+.{nd}f} ({s:.{nd}f})" if sign else f"{m:.{nd}f} ({s:.{nd}f})"


def seeds(t, c, key):
    r = rows[(t, c)]
    return f"[{r[key]}/{r['n_seeds']}]"


def table_random():
    print("| target | cells/truth | n | score | Δ vs 1x [seeds<0] | " + " | ".join(SHORT) + " |")
    print("|---|---|---|---|---|" + "---|" * len(SHORT))
    for t in TARGETS:
        for r in RATIOS:
            c = f"random_{r}x"
            if not has(t, c):
                continue
            skills = " | ".join(ms(t, c, f"skill_{m}", 3) for m in METRICS)
            print(f"| {t} | {r}x | {g(t, c, 'n_pred_mean'):.0f} | {ms(t, c, 'score')} | "
                  f"{ms(t, c, 'd_score', sign=True)} {seeds(t, c, 'seeds_dscore_neg')} | {skills} |")


def table_nbhd_raw():
    print("| target | nbhd floor/ceiling | " + " | ".join(f"{r}x" for r in RATIOS) + " |")
    print("|---|---|" + "---|" * len(RATIOS))
    for t in TARGETS:
        if not has(t, "random_1x"):
            continue
        fc = g(t, "random_1x", "floor_neighborhood_mmd_mean") / g(t, "random_1x", "ceil_neighborhood_mmd_mean")
        vals = " | ".join(f"{g(t, f'random_{r}x', 'nbhd_over_ceiling_mean'):.2f}" for r in RATIOS)
        print(f"| {t} | {fc:.0f} | {vals} |")


def table_refstage():
    print("| target | cells/truth | score | Δ vs refstage 1x [seeds<0] | nbhd skill | dice skill | de skill |")
    print("|---|---|---|---|---|---|---|")
    for t in TARGETS:
        for r in RATIOS:
            c = f"refstage_{r}x"
            if not has(t, c):
                continue
            print(f"| {t} | {r}x | {ms(t, c, 'score')} | {ms(t, c, 'd_score', sign=True)} "
                  f"{seeds(t, c, 'seeds_dscore_neg')} | {ms(t, c, 'skill_neighborhood_mmd', 3)} | "
                  f"{ms(t, c, 'skill_occupancy_dice', 3)} | {ms(t, c, 'skill_de_score', 3)} |")


def table_shuffle50():
    if not any(has(t, "shuffle50_1x") for t in TARGETS):
        return
    print("| target | cells/truth | score | Δ vs shuffle50 1x [seeds<0] | cost of the shuffle vs random same n | "
          "nbhd skill |")
    print("|---|---|---|---|---|---|")
    for t in TARGETS:
        for r in RATIOS:
            c = f"shuffle50_{r}x"
            if not has(t, c):
                continue
            print(f"| {t} | {r}x | {ms(t, c, 'score')} | {ms(t, c, 'd_score', sign=True)} "
                  f"{seeds(t, c, 'seeds_dscore_neg')} | {ms(t, c, 'dr_score', sign=True)} | "
                  f"{ms(t, c, 'skill_neighborhood_mmd', 3)} |")


def table_max():
    """Count levels at fractions of max_cells, with the paired change against submitting max_cells itself."""
    if not any(has(t, "random_1max") for t in TARGETS):
        return
    print("| target | max_cells (source; x truth) | pattern | level | n | score | Δ vs 1max [seeds<0] |")
    print("|---|---|---|---|---|---|---|")
    for t in TARGETS:
        if not has(t, "random_1max"):
            continue
        r1 = rows[(t, "random_1max")]
        mx = float(r1["max_cells"])
        head = f"{mx:.0f} ({r1['max_cells_source']}; {mx / float(r1['truth_cells']):.2f}x)"
        for pat in ("random", "shuffle50", "refstage"):
            for lv in MAX_LEVELS:
                c = f"{pat}_{lv}max"
                if not has(t, c):
                    continue
                print(f"| {t} | {head} | {pat} | {lv} max | {g(t, c, 'n_pred_mean'):.0f} | {ms(t, c, 'score')} | "
                      f"{ms(t, c, 'dm_score', sign=True)} {seeds(t, c, 'seeds_dmscore_neg')} |")
                head = ""


def table_blocks():
    print("| target | cells/truth | random | block (Δ same n) | head (Δ same n) | refstage | refblock (Δ) | "
          "refhead (Δ) |")
    print("|---|---|---|---|---|---|---|---|")
    for t in TARGETS:
        for r in RATIOS:
            if not has(t, f"block_{r}x"):
                continue
            cells = [f"{g(t, f'random_{r}x', 'score_mean'):.1f}"]
            for p in ("block", "head"):
                c = f"{p}_{r}x"
                cells.append(f"{g(t, c, 'score_mean'):.1f} ({g(t, c, 'dr_score_mean'):+.1f} ± "
                             f"{g(t, c, 'dr_score_sd'):.1f})")
            cells.append(f"{g(t, f'refstage_{r}x', 'score_mean'):.1f}")
            for p in ("refblock", "refhead"):
                c = f"{p}_{r}x"
                cells.append(f"{g(t, c, 'score_mean'):.1f} ({g(t, c, 'dr_score_mean'):+.1f} ± "
                             f"{g(t, c, 'dr_score_sd'):.1f})")
            print(f"| {t} | {r}x | " + " | ".join(cells) + " |")


def table_block_positions():
    """Where the blocks started (row of the pool), as a check that the centre draw varies across seeds."""
    print("| target | cells/truth | block start min-max (pool) | refblock start min-max (ref pool) |")
    print("|---|---|---|---|")
    for t in TARGETS:
        for r in RATIOS:
            cells = []
            for p in ("block", "refblock"):
                c = f"{p}_{r}x"
                lo, hi, pool = g(t, c, "block_start_min"), g(t, c, "block_start_max"), g(t, c, "pool_cells_mean")
                cells.append("" if math.isnan(lo) else f"{lo:.0f}-{hi:.0f} ({pool:.0f})")
            if any(cells):
                print(f"| {t} | {r}x | " + " | ".join(cells) + " |")


def table_block_metrics():
    print("| target | condition | score | " + " | ".join(SHORT) + " | shape group |")
    print("|---|---|---|" + "---|" * len(SHORT) + "---|")
    for t in TARGETS:
        for c in ("random_1x", "block_1x", "head_1x", "block_8x", "head_8x"):
            if not has(t, c):
                continue
            skills = " | ".join(f"{g(t, c, f'skill_{m}_mean'):.2f}" for m in METRICS)
            print(f"| {t} | {c} | {g(t, c, 'score_mean'):.1f} | {skills} | {g(t, c, 'group_shape_growth_mean'):.2f} |")


def table_jitter():
    fs = ["0.25", "0.5", "1.0", "2.0"]
    print("| target | density / unit | jitter RMS at f=1 | " + " | ".join(f"f={f}" for f in fs) + " |")
    print("|---|---|---|" + "---|" * len(fs))
    for t in TARGETS:
        for r, unit in (("1", "own"), ("4", "own"), ("4", "truth")):
            if not has(t, f"jitter_{r}x_1.0{unit}"):
                continue
            cells = []
            for f in fs:
                c = f"jitter_{r}x_{f}{unit}"
                cells.append(f"{g(t, c, 'dr_score_mean'):+.2f} ({g(t, c, 'dr_score_sd'):.2f}) "
                             f"{seeds(t, c, 'seeds_drscore_neg')}")
            jr = g(t, f"jitter_{r}x_1.0{unit}", "jitter_rms_mean")
            print(f"| {t} | {r}x / {unit} | {jr:.1f} | " + " | ".join(cells) + " |")


def table_jitter_channels():
    print("| target | condition | Δ score | Δ dice skill | Δ nbhd skill | Δ d2 skill | Δ scale skill | "
          "nbhd raw/ceiling |")
    print("|---|---|---|---|---|---|---|---|")
    for t in TARGETS:
        for c in ("jitter_1x_1.0own", "jitter_1x_2.0own", "jitter_4x_2.0truth"):
            if not has(t, c):
                continue
            d = [g(t, c, f"dr_skill_{m}_mean") for m in ("occupancy_dice", "neighborhood_mmd", "d2_shape",
                                                           "scale_log_ratio")]
            print(f"| {t} | {c} | {g(t, c, 'dr_score_mean'):+.2f} | " + " | ".join(f"{v:+.3f}" for v in d)
                  + f" | {g(t, c, 'nbhd_over_ceiling_mean'):.2f} |")


def table_anchors():
    print("| target | ref (floor) | truth_cells | " + " | ".join(f"{s} floor / ceiling" for s in SHORT) + " |")
    print("|---|---|---|" + "---|" * len(SHORT))
    for t in TARGETS:
        if not has(t, "random_1x"):
            continue
        cells = []
        for m in METRICS:
            cells.append(f"{g(t, 'random_1x', f'floor_{m}_mean'):.4g} / {g(t, 'random_1x', f'ceil_{m}_mean'):.4g}")
        print(f"| {t} | {REFS.get(t, '?')} | {rows[(t, 'random_1x')]['truth_cells']} | " + " | ".join(cells) + " |")


for name, fn in [("anchors", table_anchors), ("random", table_random), ("nbhd_raw", table_nbhd_raw),
                 ("refstage", table_refstage), ("shuffle50", table_shuffle50), ("max_cells_levels", table_max),
                 ("blocks", table_blocks), ("block_positions", table_block_positions),
                 ("block_metrics", table_block_metrics), ("jitter", table_jitter),
                 ("jitter_channels", table_jitter_channels)]:
    print(f"\n<!-- {name} -->")
    fn()
