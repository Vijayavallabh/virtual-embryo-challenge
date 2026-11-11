#!/usr/bin/env python
"""Price the rule fixed before the screened read of the Task 1 page statement (docs/PLAN.md, decision 2026-11-11 (c)).

Target-free: it reads only the committed per-seed payoffs of the T1 grid (practice task T1/E9.5->E10.5:chamber, 10
seeds, commit 90b25d3; `results/per_seed/t1-grid/`). For each of the 15 tissue mixes taken as the truth, it scores the
pair of uploads the rule would leave if the statement gave that mix's shares exactly, as the better of the two paired
by seed (as `payoff.py hedge` does), against today's pair (s1 bown_cown, s2 bown_c0.6):
  - cardiomyocyte share stated: the upload whose share is nearer moves to it (s1 if nearer E10.5's own 0.265 than
    0.60, ties to s2), the other stays;
  - blood share stated: both uploads move to it;
  - both stated: both take the blood share, then the nearer one the cardiomyocyte share.
`cm_stated` prices a cardiomyocyte share that converts exactly to the truth's share of heart cells (a share stated among
non-blood cells, or blood at E10.5's own). `cm_stated_of_all_cells` prices the rule as written for a share stated among
all the target's cells with blood unstated: converted with E10.5's own blood (p / (1 - 0.2409)), then snapped to the
nearest grid prediction, so it is approximate. Blood 0 and one-sided bounds are off the grid and unpriced.
Also priced and not taken: with the cardiomyocyte share stated, s1 at it and s2 at it with blood at 5% or 40%.

Usage (repo root): python -I -B experiments/s2_hedge/statement_rule.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = HERE / "results" / "per_seed" / "t1-grid"
OUT = HERE / "results" / "t1_grid"
OWN_CM, S2_CM = 0.265, 0.60  # E10.5's own cardiomyocyte share of heart cells, and today's s2
OWN_B, CUT = 0.24092736063377407, 0.4325  # E10.5's own blood share of chamber cells; s1 below CUT, s2 from it
GRID_CM = {"0.15": 0.15, "own": 0.26498350088251094, "0.45": 0.45, "0.6": 0.60, "0.75": 0.75}
GRID_B = {"own": OWN_B, "0.05": 0.05, "0.4": 0.40}


def load(runs: Path) -> tuple[dict, list[str], list[str], list[int]]:
    parts = [json.loads(q.read_text()) for d in sorted(runs.iterdir()) for q in sorted(d.glob("*.json"))]
    merged: dict[int, dict] = {}
    for part in parts:
        for r in part["rows"]:
            merged.setdefault(r["seed"], {}).update(r["scores"])
    tv, pv = parts[0]["truth_values"], parts[0]["pred_values"]
    seeds = sorted(merged)
    if any(set(merged[k]) != set(tv) for k in seeds):
        sys.exit("some seeds miss some truths")
    return merged, tv, pv, seeds


def parse(label: str) -> tuple[str, str]:
    b, c = label.split("_")
    return b[1:], c[1:]


def moved_cm(blood: str, c: str) -> tuple[str, str]:
    p = OWN_CM if c == "own" else float(c)
    s1_moves = abs(p - OWN_CM) < abs(p - S2_CM)
    return (f"b{blood}_c{c}", f"b{blood}_c0.6") if s1_moves else (f"b{blood}_cown", f"b{blood}_c{c}")


def cm_of_all_cells(t: str) -> tuple[str, str]:
    """The rule as written for a cardiomyocyte share stated among all the target's (chamber) cells, blood unstated."""
    b, c = parse(t)
    x = GRID_CM[c] * (1 - GRID_B[b]) / (1 - OWN_B)  # the stated share of all cells, converted with E10.5's own blood
    if x > 1:
        return ("bown_cown", "bown_c0.6")
    snap = min(GRID_CM, key=lambda g: abs(GRID_CM[g] - x))
    return (f"bown_c{snap}", "bown_c0.6") if x < CUT else ("bown_cown", f"bown_c{snap}")


RULES = {
    "today": lambda t: ("bown_cown", "bown_c0.6"),
    "cm_stated": lambda t: moved_cm("own", parse(t)[1]),
    "cm_stated_of_all_cells": cm_of_all_cells,
    "blood_stated": lambda t: (f"b{parse(t)[0]}_cown", f"b{parse(t)[0]}_c0.6"),
    "both_stated": lambda t: moved_cm(*parse(t)),
    "not_taken_cm_stated_s2_blood_0.05": lambda t: (f"bown_c{parse(t)[1]}", f"b0.05_c{parse(t)[1]}"),
    "not_taken_cm_stated_s2_blood_0.4": lambda t: (f"bown_c{parse(t)[1]}", f"b0.4_c{parse(t)[1]}"),
}


def main() -> None:
    merged, tv, pv, seeds = load(RUNS)

    def better(t: str, pair: tuple[str, str]) -> np.ndarray:
        for p in pair:
            if p not in pv:
                sys.exit(f"{p} is not a grid prediction")
        return np.array([max(merged[k][t][pair[0]], merged[k][t][pair[1]]) for k in seeds])

    # orientation: the README's blood-own slice (truth c0.15: s1 71.15; truth c0.6: bown_c0.6 71.22)
    assert abs(np.mean([merged[k]["bown_c0.15"]["bown_cown"] for k in seeds]) - 71.15) < 0.01
    assert abs(np.mean([merged[k]["bown_c0.6"]["bown_c0.6"] for k in seeds]) - 71.22) < 0.01
    priors = {"uniform": {t: 1 / len(tv) for t in tv},
              "own_third": {t: 1 / 3 if t == "bown_cown" else (2 / 3) / (len(tv) - 1) for t in tv},
              "blood_own_only": {t: 1 / 5 for t in tv if t.startswith("bown_")}}
    res: dict = {"source": "results/per_seed/t1-grid (commit 90b25d3)", "seeds": seeds, "priors": {}, "per_truth": {}}
    for name, w in priors.items():
        today = sum(w[t] * better(t, RULES["today"](t)) for t in w)
        row = {"today": round(float(today.mean()), 3)}
        for rule, pair in RULES.items():
            if rule == "today":
                continue
            d = sum(w[t] * better(t, pair(t)) for t in w) - today
            row[rule] = {"gain": round(float(d.mean()), 3), "sd": round(float(d.std(ddof=1)), 3),
                         "seeds_ahead": int((d > 1e-9).sum())}
        res["priors"][name] = row
    for t in tv:
        base = better(t, RULES["today"](t)).mean()
        res["per_truth"][t] = {"today": round(float(base), 2),
                               **{r: round(float(better(t, p(t)).mean() - base), 2) for r, p in RULES.items()
                                  if r != "today"}}
    lines = ["Gain of the statement rule over today's T1 pair (s1 bown_cown, s2 bown_c0.6), better of two paired by "
             f"seed, {len(seeds)} seeds, if the statement gives the truth's shares exactly",
             f"{'rule':>36s} " + " ".join(f"{n:>16s}" for n in priors)]
    for rule in [r for r in RULES if r != "today"]:
        lines.append(f"{rule:>36s} " + " ".join(f"{res['priors'][n][rule]['gain']:+16.2f}" for n in priors))
    lines.append(f"{'today (better of two)':>36s} " + " ".join(f"{res['priors'][n]['today']:16.2f}" for n in priors))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "statement_rule.json").write_text(json.dumps(res, indent=1) + "\n")
    (OUT / "statement_rule.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
