"""Upper-tail check of two options put to the captain (D3, R2) against today's pair, from committed results only.

Written for the captain's rulings of 2026-11-05 (docs/PLAN.md decision log 2026-11-05 (c)): he asked for the choices
that maximise our chances of a top-3 place, so each option is judged by P(Overall >= c) in reading A, not by its mean.

Aggregation with experiments/gate/estimate.py's own loaders and outcome rules, not a try: no target, no score computed,
no leaderboard figure. For each option, reading A (each board at the better of its two uploads) is drawn 200,000 times
with common random numbers (one outcome per board, the boards independent), and P(Overall >= c) is printed for every
c across our own upper range (the median to the 99.5th percentile of today's pair), so no outside threshold enters.
  D3: T1's s2 at cardiomyocytes 0.75 (the grid's tail column) instead of 0.60.
  R2: the embryo's s2 = fgw_aligned (E11, emb_aligned_v1) instead of fgw_shape0.
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate as E  # noqa: E402

door, checks = E.Door(REPO), []
t1, _, t1_steps = E.board_t1(door, checks)
heart, _, _ = E.board_heart(door, checks)
embryo = E.board_embryo(door, checks)
_, t3, _ = E.board_t3(door, checks)  # with K2a, as the record's reproducible reading

grid, hedge = E.load_runs(door, E.T1_GRID_RUNS), door.json(E.T1_GRID_HEDGE)
prior, sc = hedge["priors"][E.PRIOR], grid["scores"]
t1_d3 = E.grid_outcomes({t: sc[t][E.T1_S1] for t in prior}, {t: sc[t]["bown_c0.75"] for t in prior}, prior,
                        shift=t1_steps["m1_minus_b0"])

parts = []
for shape_rel, _fgw_rel in E.EMBRYO:
    s = E.load_task(door, shape_rel, ("floor", E.EMBRYO_S1, E.EMBRYO_S2))
    al = E.load_task(door, shape_rel.replace("emb_shape_v1", "emb_aligned_v1"), ("fgw_aligned",))
    if s["key"] != al["key"] or s["seeds"] != al["seeds"]:
        sys.exit(f"{shape_rel}: emb_aligned_v1 is not the same backtest and seeds")
    parts.append((1.0, E.grid_outcomes({"b": s[E.EMBRYO_S1]}, {"b": al["fgw_aligned"]}, {"b": 1.0})))
embryo_r2 = E.mixture(parts)

for a, b, what in ((t1, t1_d3, "T1"), (embryo, embryo_r2, "embryo")):
    if a["a"].shape != b["a"].shape or not np.allclose(a["w"], b["w"]):
        sys.exit(f"{what}: the option's outcomes are not paired with today's")

options = {"today": {"T1": t1, "heart": heart, "embryo": embryo, "T3": t3},
           "D3": {"T1": t1_d3, "heart": heart, "embryo": embryo, "T3": t3},
           "R2": {"T1": t1, "heart": heart, "embryo": embryo_r2, "T3": t3},
           "D3+R2": {"T1": t1_d3, "heart": heart, "embryo": embryo_r2, "T3": t3}}
rng = np.random.default_rng(E.RNG_SEED)
idx = {b: rng.choice(t1["w"].size if b == "T1" else options["today"][b]["w"].size, size=E.N_DRAWS,
                     p=options["today"][b]["w"]) for b in E.BOARDS}
tot = {n: sum(E.OVERALL_WEIGHT[b] * o[b]["a"][idx[b]] for b in E.BOARDS) for n, o in options.items()}

print("reading A means:", {n: round(E.readings(o)["A"]["Overall"], 2) for n, o in options.items()})
print(f"board A means:  T1 today {E.reading(t1, 'a'):.2f}, D3 {E.reading(t1_d3, 'a'):.2f}; "
      f"embryo today {E.reading(embryo, 'a'):.2f}, R2 {E.reading(embryo_r2, 'a'):.2f}")
qs = (50, 75, 90, 95, 99, 99.5)
for n, x in tot.items():
    print(f"{n:6s} percentiles", {q: round(float(np.percentile(x, q)), 1) for q in qs})
lo, hi = np.percentile(tot["today"], 50), np.percentile(tot["today"], 99.5)
cs = np.arange(np.floor(lo), np.ceil(hi) + 1, 2.0)
print("c     " + " ".join(f"{n:>8s}" for n in tot) + "   (P(Overall_A >= c), %)")
better = {n: 0 for n in tot if n != "today"}
worse = {n: 0 for n in tot if n != "today"}
for c in cs:
    p = {n: 100 * float(np.mean(x >= c)) for n, x in tot.items()}
    for n in better:
        better[n] += p[n] > p["today"] + 1e-9
        worse[n] += p[n] < p["today"] - 1e-9
    print(f"{c:5.0f} " + " ".join(f"{p[n]:8.2f}" for n in tot))
print("thresholds where the option beats / trails today:", {n: (better[n], worse[n]) for n in better},
      "of", len(cs))
print("self-checks run:", len(checks))
