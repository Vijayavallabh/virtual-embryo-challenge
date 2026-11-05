"""Price T3's second upload, today's refused s1 against the barred s1, under both readings of a floored severity.

Aggregation of committed results with stated judgement cells. Not a try: no target is read and no score is computed
(experiments/harness/TRIES.md row N10, registered with its decision rule before the run).

**Worlds.** U1 ~ Bern(0.85), anchor ~ Bern(p_anchor) and U2 ~ Bern(p_u2), independent:
- U1: the beta-catenin file was imaged apart from WT E8.75;
- anchor: the component shared by the two factors carries to the test;
- U2: the file shares Gata4's shallow depth.
(`experiments/research_2026_11/REPORT.md` section 3; task h's priors; the counter-review.)

**What each upload scores in each world.** "Committed" cells draw one Gata4-replicate seed row of the named stand-in
(paired by seed and replicate across uploads, as estimate.py pairs them). "Judgement" cells are task h's values: a
normal with SD 3 and a stated floor probability.

| Upload | World | Score |
|---|---|---|
| s1 | not U1 | judgement 46, floored with p 0.9 |
| s1 | U1, anchor, U2 | committed `xrep:F+depth` (sensitivity: judgement 66, floored with p 0) |
| s1 | U1, anchor, not U2 | judgement 50, floored with p 0.9 |
| s1 | U1, no anchor, U2 | committed `fac:depth` (the depth step without a working factor) |
| s1 | U1, no anchor, not U2 | judgement 46, floored with p 0.9 |
| Today's s2 | every world | committed `blind:G2_bcat_central_d50` |
| Barred s1 | not U1 | judgement 47, floored with p 0.6 |
| Barred s1 | U1, anchor | committed `transfer:mab21l2_all` |
| Barred s1 | U1, no anchor | committed `transfer:mab21l2_vs95` |

**The two readings.**
- "Half": a floored severity scores skill 0.5, as the committed values are.
- "Zero": every floored draw loses 0.25 x 0.5 x 100 = 12.5 T3 points.
- The mixture uses P(zero) = 0.35, the working answer of `docs/ORGANIZER_QUESTIONS.md`.

**Overall.** Overall = T1 + T3 + (heart + embryo)/2. The other boards come from estimate.py's reading A ("a": the
better of two; "b": s1, the consensus-like first upload), with common random numbers. The levels c are today's own
75th, 90th, 95th and 99th percentiles, so no outside threshold enters.

**Rank-aware pricing** (task e, rule K3): C3 = c + lambda (S1 - mean S1) + eta, with eta ~ N(0, 3^2) and S1 the
first-upload total (T3's s1 world score included). P_top3 = mean(Overall >= C3), for lambda in {0, 0.5, 1}.

**Decision rule (row N10).** At the band the checks select, the AI assistant recommends the swap only if both hold:
- the gain in E[T3] under the mixture is at least +0.5;
- the change in P_top3 is at least 0 at the 90th and 95th levels, for lambda 0.5 and 1, under the mixture.
Otherwise it recommends keeping today's s2. The captain decides either way.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate as E  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
P_U1 = 0.85
P_ZERO = 0.35
FLOOR_COST = 12.5  # 100 x weight 0.25 x skill 0.5
GRID_U2 = (0.5, 0.7, 0.85)
GRID_ANCHOR = (0.35, 0.5, 0.65)
LAMBDAS = (0.0, 0.5, 1.0)
LEVELS = (75, 90, 95, 99)
SD_J = 3.0
SEED = 20261106
COMMITTED = ("xrep:F+depth", "fac:depth", "blind:G2_bcat_central_d50", "transfer:mab21l2_all",
             "transfer:mab21l2_vs95")


def load_rows(door) -> dict[str, dict[str, np.ndarray]]:
    """Per variant: t3 and sev_floor over the 20 (replicate, seed) rows, in one shared row order."""
    rows = list(csv.DictReader(open(door.path(E.T3_RECORDS))))
    out = {}
    for v in COMMITTED:
        rs = sorted([r for r in rows if r["variant"] == v], key=lambda r: (r["target"], int(r["seed"])))
        if len(rs) != 20:
            raise SystemExit(f"{v}: {len(rs)} rows, expected 20")
        out[v] = {"key": [(r["target"], int(r["seed"])) for r in rs],
                  "t3": np.array([float(r["t3"]) for r in rs]),
                  "floor": np.array([r["sev_floor"] in ("True", "true", "1", "1.0") for r in rs])}
    keys = {tuple(map(tuple, out[v]["key"])) for v in COMMITTED}
    if len(keys) != 1:
        raise SystemExit("stand-ins are not paired on the same (replicate, seed) rows")
    return out


def simulate(rows, others, p_u2: float, p_anchor: float, best66: bool, rng) -> dict:
    n = E.N_DRAWS
    u1 = rng.random(n) < P_U1
    anc = rng.random(n) < p_anchor
    u2 = rng.random(n) < p_u2
    row = rng.integers(0, 20, n)  # one (replicate, seed) row per draw, shared by every committed cell

    def committed(v):
        return rows[v]["t3"][row], rows[v]["floor"][row]

    def judged(m, pf):
        return rng.normal(m, SD_J, n), rng.random(n) < pf

    s1, f1 = judged(46.0, 0.9)  # not U1, and U1/no anchor/not U2 (same values)
    for mask, (x, f) in (
            (u1 & anc & u2, judged(66.0, 0.0) if best66 else committed("xrep:F+depth")),
            (u1 & anc & ~u2, judged(50.0, 0.9)),
            (u1 & ~anc & u2, committed("fac:depth"))):
        s1, f1 = np.where(mask, x, s1), np.where(mask, f, f1)
    s2, f2 = committed("blind:G2_bcat_central_d50")
    bx, bf = judged(47.0, 0.6)
    for mask, (x, f) in ((u1 & anc, committed("transfer:mab21l2_all")),
                         (u1 & ~anc, committed("transfer:mab21l2_vs95"))):
        bx, bf = np.where(mask, x, bx), np.where(mask, f, bf)
    idx = {k: rng.choice(o["w"].size, size=n, p=o["w"]) for k, o in others.items()}
    rest_a = sum(E.OVERALL_WEIGHT[k] * o["a"][idx[k]] for k, o in others.items())
    rest_b = sum(E.OVERALL_WEIGHT[k] * o["b"][idx[k]] for k, o in others.items())
    eta = rng.normal(0.0, 3.0, n)
    res = {}
    for reading in ("half", "zero"):
        c = FLOOR_COST if reading == "zero" else 0.0
        a1, a2, ab = s1 - c * f1, s2 - c * f2, bx - c * bf
        t3_today, t3_swap = np.maximum(a1, a2), np.maximum(a1, ab)
        o_today, o_swap = rest_a + t3_today, rest_a + t3_swap
        s1tot = rest_b + a1
        lv = {q: float(np.percentile(o_today, q)) for q in LEVELS}
        p = {}
        for lam in LAMBDAS:
            thr_shift = lam * (s1tot - s1tot.mean()) + eta
            for q, cl in lv.items():
                pt, ps = np.mean(o_today >= cl + thr_shift), np.mean(o_swap >= cl + thr_shift)
                p[f"lambda{lam}_q{q}"] = {"today": float(pt), "swap": float(ps), "delta": float(ps - pt)}
        res[reading] = {"E_T3_today": float(t3_today.mean()), "E_T3_swap": float(t3_swap.mean()),
                        "dE_T3": float(t3_swap.mean() - t3_today.mean()),
                        "E_overall_today": float(o_today.mean()), "levels": lv, "P_top3": p,
                        "floor_rate_s1": float(f1.mean()), "floor_rate_s2": float(f2.mean()),
                        "floor_rate_barred": float(bf.mean())}
    mix = {"dE_T3": (1 - P_ZERO) * res["half"]["dE_T3"] + P_ZERO * res["zero"]["dE_T3"],
           "dP": {k: (1 - P_ZERO) * res["half"]["P_top3"][k]["delta"] + P_ZERO * res["zero"]["P_top3"][k]["delta"]
                  for k in res["half"]["P_top3"]}}
    rec = (mix["dE_T3"] >= 0.5 and all(mix["dP"][f"lambda{lam}_q{q}"] >= 0 for lam in (0.5, 1.0) for q in (90, 95)))
    return {"p_u2": p_u2, "p_anchor": p_anchor, "best66": best66, **res, "mixture": mix,
            "recommend_swap": bool(rec)}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="JSON output path")
    a = ap.parse_args(argv)
    door, checks = E.Door(REPO), []
    t1, _, _ = E.board_t1(door, checks)
    heart, _, _ = E.board_heart(door, checks)
    embryo = E.board_embryo(door, checks)
    rows = load_rows(door)
    stats = {v: {"mean_t3": float(rows[v]["t3"].mean()), "floored": int(rows[v]["floor"].sum())} for v in COMMITTED}
    others = {"T1": t1, "heart": heart, "embryo": embryo}
    cells = []
    for best66 in (False, True):
        for pu in GRID_U2:
            for pa in GRID_ANCHOR:
                cells.append(simulate(rows, others, pu, pa, best66, np.random.default_rng(SEED)))
    out = {"script": "experiments/gate/t3_backup_pricing.py", "self_checks": len(checks), "standins": stats,
           "p_u1": P_U1, "p_zero": P_ZERO, "cells": cells}
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"self-checks run: {len(checks)}; stand-ins: " +
          "; ".join(f"{v} {s['mean_t3']:.2f} ({s['floored']}/20 floored)" for v, s in stats.items()))
    print("best66 p_u2 p_anchor | dE_T3 half / zero / mix | dP(top3) mix at q90, q95 for lambda 0, 0.5, 1 | recommend")
    for c in cells:
        dp = c["mixture"]["dP"]
        ps = " ".join(f"{100 * dp[f'lambda{lam}_q{q}']:+.2f}" for lam in LAMBDAS for q in (90, 95))
        print(f"{c['best66']!s:5} {c['p_u2']:.2f} {c['p_anchor']:.2f} | {c['half']['dE_T3']:+.2f} / "
              f"{c['zero']['dE_T3']:+.2f} / {c['mixture']['dE_T3']:+.2f} | {ps} | {c['recommend_swap']}")


if __name__ == "__main__":
    main()
