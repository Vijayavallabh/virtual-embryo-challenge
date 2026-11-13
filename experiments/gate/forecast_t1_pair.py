#!/usr/bin/env python
"""What tries T32 and T33's tables do to the forecast: the built T1 pair against the other pairs those tries priced,
inside experiments/gate/forecast.py's own model of the total. An aggregation of committed results, not a try; no test
score and no leaderboard figure enters it (forecast.py's comparison lines are applied to the finished draws only).

Pairs (first upload, second upload), each arm a column of try T30, T32 or T33's scenario files:
  built       lib3_base, s1_cm60    the pair of row R6 (forecast.py's files of 2026-11-13)
  lib3_pair   lib3_base, lib3_cm60  try T32
  mid_s2      lib3_base, mid_cm60   try T33 (only when its results are committed)
  mid_both    mid_base, mid_cm60    try T33

Model: forecast.py's files of 2026-11-13, unchanged for the heart, the embryo and T3. T1 as forecast.board_t1 builds
it: the oracle-base grid's truths, E12.5's depth world, and each upload moved, per seed, by k times the zero-gap
difference between its arm and the arm the grid stands for (s1_base for the first upload, s1_cm60 for the second) in
the scenario at that depth and the truth's mapped cardiomyocyte share. k and P(deep) are forecast.py's JUDGEMENT
values; the table gives each pair at their central, low and high values and under try T30's "latest" and "flat" depth
priors. One seed and one outcome order price every pair, so the differences are paired.

Run: python3 -I -B experiments/gate/forecast_t1_pair.py [--json PATH]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate as E  # noqa: E402
import forecast as F  # noqa: E402

T32_DIR = "experiments/t1_scope/results/t32_pair"
T33_DIR = "experiments/t1_scope/results/t33_mid"
E.READABLE = (*E.READABLE, f"{T32_DIR}/", f"{T33_DIR}/")
LINES = (203.00, 205.32, 209.0)
SHARES = ("own", "0.45", "0.60")
PAIRS = {"built": ("lib3_base", "s1_cm60"), "lib3_pair": ("lib3_base", "lib3_cm60"),
         "mid_s2": ("lib3_base", "mid_cm60"), "mid_both": ("mid_base", "mid_cm60")}
SOURCES = {"s1_base": ("t30", F.T30_DIR), "s1_cm60": ("t30", F.T30_DIR), "lib3_base": ("t30", F.T30_DIR),
           "lib3_cm60": ("t32", T32_DIR), "mid_base": ("t33", T33_DIR), "mid_cm60": ("t33", T33_DIR)}


def load_arms(inp: F.Inputs, arms: tuple[str, ...]) -> dict:
    """Per seed, each arm's score by (depth, truth share), from the scenario files of the try that ran it; the ten-seed
    means are checked against that try's pricing.json."""
    out = {a: {} for a in arms}
    for a in arms:
        stem, folder = SOURCES[a]
        summ = inp.door.json(f"{folder}/pricing.json")
        F.need(not summ.get("void", False), f"{stem}'s reading is not void", inp.checks)
        for dep in F.T30_DEPTHS:
            for sh in SHARES:
                raw = inp.door.json(f"{folder}/scenarios/{stem}_d{dep}_cm{sh}.json")
                rows = sorted(raw["rows"], key=lambda r: r["seed"])
                F.need([int(r["seed"]) for r in rows] == inp.grid["seeds"], f"{stem} d{dep} cm{sh}: the grid's seeds")
                v = np.array([float(r["scores"][a]) for r in rows])
                E.same([v.mean()], [summ["arm_by_scenario"][a][f"d{dep} cm{sh}"]],
                       f"{stem} d{dep} cm{sh}: {a}'s per-seed scores give pricing.json's mean", inp.checks)
                out[a][dep, sh] = v
    return out


def board_t1(inp: F.Inputs, prior_name: str, p_deep: float, k: float, arms: dict, pair: tuple[str, str],
             depth_w: dict | None = None) -> dict:
    """forecast.board_t1 for the new files, generalised to any pair of arms: the first upload moved by
    k x (its arm - s1_base), the second by k x (its arm - s1_cm60). depth_w: explicit depth weights in place of
    forecast.py's P(deep) split (try T30's "latest" and "flat" priors)."""
    sc, prior, step = inp.grid["scores"], inp.t1_hedge_priors[prior_name], inp.t1_steps["m1_minus_b0"]
    a1, a2 = pair

    def by_truth(d: float) -> list:
        parts = []
        for t, pt in prior.items():
            sh = F.T1_SHARE[t.split("_")[1]]
            s1 = sc[t][E.T1_S1] + k * (arms[a1][d, sh] - arms["s1_base"][d, sh])
            s2 = sc[t][E.T1_S2] + k * (arms[a2][d, sh] - arms["s1_cm60"][d, sh])
            parts.append((pt, E.grid_outcomes({t: s1}, {t: s2}, {t: 1.0}, shift=step), f"truth {t}"))
        return parts

    if depth_w is None:
        depth_w = {**{d: (1 - p_deep) * f for d, f in F.T30_SHALLOW_SPLIT.items()},
                   **{d: p_deep * f for d, f in F.T30_DEEP_SPLIT.items()}}
    F.need(abs(sum(depth_w.values()) - 1) < 1e-12, "depth weights sum to 1")
    return F.mix([(wd, F.mix(by_truth(d)), f"depth {d}") for d, wd in depth_w.items()])


def price(inp: F.Inputs, arms: dict, pairs: dict, p_deep: float, k: float, depth_w: dict | None = None) -> dict:
    par = {**F.central(), "t1_p_deep": p_deep, "t1_k": k}
    bd = F.boards(inp, "new", par)
    out = {}
    for name, pair in pairs.items():
        t1 = board_t1(inp, par["t1_prior"], p_deep, k, arms, pair, depth_w)
        inside = (t1["a"] >= 0) & (t1["a"] <= 100) & (t1["b"] >= 0) & (t1["b"] <= 100)
        F.need(bool(np.all(inside)), f"T1 {name}: scores in [0, 100]")
        if name == "built" and depth_w is None:
            F.need(F.same_dist(t1, bd["T1"]), "the built pair is forecast.py's T1 board for the new files")
        tot, _ = F.draw({**bd, "T1": t1})
        out[name] = {"t1_mean_better_of_two": E.reading(t1, "a"), "t1_mean_first_upload": E.reading(t1, "b"),
                     "total_mean": float(tot.mean()), **{f"p_{c:.2f}": float(np.mean(tot >= c)) for c in LINES}}
    out["minus_built"] = {n: {k2: out[n][k2] - out["built"][k2] for k2 in out["built"]} for n in pairs if n != "built"}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", metavar="PATH")
    a = ap.parse_args(argv)
    inp = F.Inputs()
    have_t33 = (F.REPO / T33_DIR / "pricing.json").exists()
    pairs = {n: p for n, p in PAIRS.items() if have_t33 or not n.startswith("mid")}
    arms = load_arms(inp, tuple(sorted({"s1_base", "s1_cm60", *(x for p in pairs.values() for x in p)})))
    F.need(all(np.allclose(arms["lib3_base"][key] - arms["s1_base"][key], inp.t30[key]) for key in inp.t30),
           "the first upload's move is forecast.py's (lib3_base - s1_base, per seed)", inp.checks)
    pd, kk = F.JUDGEMENT["t1_p_deep"], F.JUDGEMENT["t1_k"]
    latest = {0.35: 1 / 3, 0.68: 1 / 3, 1.0: 1 / 3}  # try T30's "latest" depth prior: E10.5's three libraries
    flat = {d: 0.25 for d in F.T30_DEPTHS}  # try T30's "flat" depth prior
    rows = [("central", pd[0], kk[0], None), ("P(deep) low", pd[1], kk[0], None), ("P(deep) high", pd[2], kk[0], None),
            ("k low", pd[0], kk[1], None), ("k high", pd[0], kk[2], None), ("P(deep) low, k high", pd[1], kk[2], None),
            ("T30 latest depths", 1 / 3, kk[0], latest), ("T30 latest, k high", 1 / 3, kk[2], latest),
            ("T30 flat depths", 0.5, kk[0], flat)]
    res = {"label": "the built T1 pair against the pairs tries T32 and T33 priced, inside forecast.py's model: "
                    "aggregation of committed results, not a try; the comparison lines are applied to the finished "
                    "draws only",
           "pairs": {n: list(p) for n, p in pairs.items()}, "draws": F.N_DRAWS, "rng_seed": F.SEED,
           "lines": {"203.00": F.LINES[203.00], "205.32": F.LINES[205.32], "209.00": F.LINES[209.0]}, "table": {}}
    others = [n for n in pairs if n != "built"]
    print("the built T1 pair, then each other pair as its difference from the built pair; forecast.py's model, files "
          "of 2026-11-13")
    print(f"  {'setting':<21} {'P(deep)':>7} {'k':>5} | {'pair':<9} {'T1 mean':>8} {'total':>8} {'P(>=203)':>9} "
          f"{'P(>=205.32)':>11} {'P(>=209)':>9}")
    for name, p, k, dw in rows:
        r = price(inp, arms, pairs, p, k, dw)
        res["table"][name] = {"p_deep": p, "k": k, "depth_weights": dw and {str(d): w for d, w in dw.items()}, **r}
        b = r["built"]
        print(f"  {name:<21} {p:7.3f} {k:5.2f} | {'built':<9} {b['t1_mean_better_of_two']:8.2f} {b['total_mean']:8.2f} "
              f"{100 * b['p_203.00']:8.2f}% {100 * b['p_205.32']:10.2f}% {100 * b['p_209.00']:8.2f}%")
        for n in others:
            d = r["minus_built"][n]
            print(f"  {'':<21} {'':>7} {'':>5} | {n:<9} {d['t1_mean_better_of_two']:+8.2f} {d['total_mean']:+8.2f} "
                  f"{100 * d['p_203.00']:+8.2f}  {100 * d['p_205.32']:+10.2f}  {100 * d['p_209.00']:+8.2f} ")
    res["checks"] = sorted(set(inp.checks))
    if a.json:
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.json).write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
