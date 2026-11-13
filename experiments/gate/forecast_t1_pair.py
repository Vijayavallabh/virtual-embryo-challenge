#!/usr/bin/env python
"""What try T32's table does to the forecast: the built T1 pair (first upload on library 3, second on the pooled base
with cardiomyocytes at 0.60) against the lib3 pair (both on library 3, the second with cardiomyocytes at 0.60), inside
experiments/gate/forecast.py's own model of the total. An aggregation of committed results, not a try; no test score
and no leaderboard figure enters it (forecast.py's comparison lines are applied to the finished draws only).

Model: forecast.py's files of 2026-11-13, unchanged for the heart, the embryo and T3. T1 as forecast.board_t1 builds
it (the oracle-base grid's truths; E12.5's depth world; the first upload moved by k x (lib3_base - s1_base) from try
T30's scenario at that depth and the truth's cardiomyocyte share), and for the lib3 pair the second upload moved the
same way by k x (lib3_cm60 - s1_cm60) from try T32's scenario. k and P(deep) are forecast.py's JUDGEMENT values; the
table gives each pair at their central, low and high values. The same draws price both pairs (one seed, the same
outcome order), so the differences are paired.

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
T32_SUMMARY = f"{T32_DIR}/pricing.json"
E.READABLE = (*E.READABLE, f"{T32_DIR}/")
LINES = (203.00, 205.32, 209.0)


def load_t32(inp: F.Inputs) -> dict:
    """Per seed, lib3_cm60 - s1_cm60 by (depth, truth share), from T32's scenario files; checked against its summary
    and against T30's stored s1_cm60 (T32's reproduction gate, restated here on the means)."""
    summ = inp.door.json(T32_SUMMARY)
    F.need(not summ["void"], "T32's reading is not void", inp.checks)
    out = {}
    for dep in F.T30_DEPTHS:
        for sh in ("own", "0.45", "0.60"):
            raw = inp.door.json(f"{T32_DIR}/scenarios/t32_d{dep}_cm{sh}.json")
            rows = sorted(raw["rows"], key=lambda r: r["seed"])
            F.need([int(r["seed"]) for r in rows] == inp.grid["seeds"], f"T32 d{dep} cm{sh}: seeds are the T1 grid's")
            new = np.array([float(r["scores"]["lib3_cm60"]) for r in rows])
            old = np.array([float(r["scores"]["s1_cm60"]) for r in rows])
            scen = f"d{dep} cm{sh}"
            E.same([new.mean(), old.mean()], [summ["arm_by_scenario"]["lib3_cm60"][scen],
                                              summ["arm_by_scenario"]["s1_cm60"][scen]],
                   f"T32 {scen}: per-seed files give pricing.json's arm means", inp.checks)
            out[dep, sh] = new - old
    return out


def board_t1(inp: F.Inputs, prior_name: str, p_deep: float, k: float, t32: dict | None,
             depth_w: dict | None = None) -> dict:
    """forecast.board_t1 for the new files, with the second upload also moved when t32 is given (the lib3 pair).
    depth_w: explicit depth weights in place of forecast.py's P(deep) split (try T30's "latest" and "flat" priors)."""
    sc, prior, step = inp.grid["scores"], inp.t1_hedge_priors[prior_name], inp.t1_steps["m1_minus_b0"]

    def by_truth(d: float) -> list:
        parts = []
        for t, pt in prior.items():
            sh = F.T1_SHARE[t.split("_")[1]]
            s1 = sc[t][E.T1_S1] + k * inp.t30[d, sh]
            s2 = sc[t][E.T1_S2] + (k * t32[d, sh] if t32 is not None else 0.0)
            parts.append((pt, E.grid_outcomes({t: s1}, {t: s2}, {t: 1.0}, shift=step), f"truth {t}"))
        return parts

    if depth_w is None:
        depth_w = {**{d: (1 - p_deep) * f for d, f in F.T30_SHALLOW_SPLIT.items()},
                   **{d: p_deep * f for d, f in F.T30_DEEP_SPLIT.items()}}
    F.need(abs(sum(depth_w.values()) - 1) < 1e-12, "depth weights sum to 1")
    return F.mix([(wd, F.mix(by_truth(d)), f"depth {d}") for d, wd in depth_w.items()])


def price(inp: F.Inputs, t32: dict, p_deep: float, k: float, depth_w: dict | None = None) -> dict:
    par = {**F.central(), "t1_p_deep": p_deep, "t1_k": k}
    bd = F.boards(inp, "new", par)
    out = {}
    for name, shift in (("built", None), ("lib3_pair", t32)):
        t1 = board_t1(inp, par["t1_prior"], p_deep, k, shift, depth_w)
        inside = (t1["a"] >= 0) & (t1["a"] <= 100) & (t1["b"] >= 0) & (t1["b"] <= 100)
        F.need(bool(np.all(inside)), f"T1 {name}: scores in [0, 100]")
        if shift is None and depth_w is None:
            F.need(F.same_dist(t1, bd["T1"]), "the built pair is forecast.py's T1 board for the new files")
        tot, _ = F.draw({**bd, "T1": t1})
        out[name] = {"t1_mean_better_of_two": E.reading(t1, "a"), "total_mean": float(tot.mean()),
                     **{f"p_{c:.2f}": float(np.mean(tot >= c)) for c in LINES}}
    out["lib3_pair_minus_built"] = {k2: out["lib3_pair"][k2] - out["built"][k2] for k2 in out["built"]}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", metavar="PATH")
    a = ap.parse_args(argv)
    inp = F.Inputs()
    t32 = load_t32(inp)
    pd, kk = F.JUDGEMENT["t1_p_deep"], F.JUDGEMENT["t1_k"]
    res = {"label": "the built T1 pair against the lib3 pair inside forecast.py's model: aggregation of committed "
                    "results, not a try; the comparison lines are applied to the finished draws only",
           "t32_commit": inp.door.json(T32_SUMMARY)["commit"], "draws": F.N_DRAWS, "rng_seed": F.SEED,
           "lines": {"203.00": F.LINES[203.00], "205.32": F.LINES[205.32], "209.00": F.LINES[209.0]}, "table": {}}
    latest = {0.35: 1 / 3, 0.68: 1 / 3, 1.0: 1 / 3}  # try T30's "latest" depth prior: E10.5's three libraries
    flat = {d: 0.25 for d in F.T30_DEPTHS}  # try T30's "flat" depth prior
    rows = [("central", pd[0], kk[0], None), ("P(deep) low", pd[1], kk[0], None), ("P(deep) high", pd[2], kk[0], None),
            ("k low", pd[0], kk[1], None), ("k high", pd[0], kk[2], None), ("P(deep) low, k high", pd[1], kk[2], None),
            ("T30 latest depths", 1 / 3, kk[0], latest), ("T30 latest, k high", 1 / 3, kk[2], latest),
            ("T30 flat depths", 0.5, kk[0], flat)]
    print("the built T1 pair / the lib3 pair (difference), forecast.py's model, files of 2026-11-13")
    print(f"  {'setting':<22} {'P(deep)':>7} {'k':>5} | {'T1 mean':>22} | {'total mean':>24} | "
          f"{'P(>=203)':>22} | {'P(>=205.32)':>22} | {'P(>=209)':>20}")
    for name, p, k, dw in rows:
        r = price(inp, t32, p, k, dw)
        res["table"][name] = {"p_deep": p, "k": k, "depth_weights": dw and {str(d): w for d, w in dw.items()}, **r}
        b, n, d = r["built"], r["lib3_pair"], r["lib3_pair_minus_built"]

        def cell(key: str, pct: bool = False, b=b, n=n, d=d) -> str:
            s = 100.0 if pct else 1.0
            return f"{s * b[key]:6.2f} / {s * n[key]:6.2f} ({s * d[key]:+5.2f})"
        print(f"  {name:<22} {p:7.3f} {k:5.2f} | {cell('t1_mean_better_of_two'):>22} | {cell('total_mean'):>24} | "
              f"{cell('p_203.00', True):>22} | {cell('p_205.32', True):>22} | {cell('p_209.00', True):>20}")
    res["checks"] = sorted(set(inp.checks))
    if a.json:
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.json).write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
