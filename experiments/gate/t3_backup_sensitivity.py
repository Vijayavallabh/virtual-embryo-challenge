"""Sensitivity of T3's backup pricing (t3_backup_pricing.py, TRIES N10) to the one world that carries its gain.

TRIES N11; aggregation of committed results with stated judgement cells, not a try.

**Why.** A side agent's note to the captain (2026-11-06): nearly all of N10's gain comes from one world,
U1 & anchor & not-U2, with probability 0.85 x p_anchor x (1 - p_u2), about 1 in 6 at the selected band. In that world:
- s1's score is a judgement cell (50, floored with probability 0.9);
- the barred file's score is `transfer:mab21l2_all` as measured on the Gata4 replicates, whose files are shallow.
The beta-catenin file, in this world, is not.

**What varies, in that world only.** Every other world is exactly as in N10.
- s1 is a judgement cell (m, f) for (m, f) in (50, 0.9), (53, 0.6), (56, 0.3) and (59, 0.0).
- The barred file performs as measured with probability q, and as `transfer:mab21l2_vs95` otherwise (its gate fails;
  the counter-review's bimodal reading). q runs over 1.0, 0.75, 0.5 and 0.25.

**Output:** dE_T3 (half, zero, mixture), and dP_top3 under the mixture at the 90th and 95th levels for lambda 1. It is
reported at (p_u2, p_anchor) = (0.7, 0.65), N10's selected band, and at (0.7, 0.5), the band D6's raw placement points
to. Seeds and draws are as in N10.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate as E  # noqa: E402
import t3_backup_pricing as P  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
S1_WC = ((50.0, 0.9), (53.0, 0.6), (56.0, 0.3), (59.0, 0.0))
Q_BARRED = (1.0, 0.75, 0.5, 0.25)
BANDS = ((0.7, 0.65), (0.7, 0.5))


def simulate(rows, others, p_u2, p_anchor, s1_wc, q, rng) -> dict:
    n = E.N_DRAWS
    u1 = rng.random(n) < P.P_U1
    anc = rng.random(n) < p_anchor
    u2 = rng.random(n) < p_u2
    row = rng.integers(0, 20, n)

    def committed(v):
        return rows[v]["t3"][row], rows[v]["floor"][row]

    def judged(m, pf):
        return rng.normal(m, P.SD_J, n), rng.random(n) < pf

    wc = u1 & anc & ~u2
    s1, f1 = judged(46.0, 0.9)
    for mask, (x, f) in ((u1 & anc & u2, committed("xrep:F+depth")), (wc, judged(*s1_wc)),
                         (u1 & ~anc & u2, committed("fac:depth"))):
        s1, f1 = np.where(mask, x, s1), np.where(mask, f, f1)
    s2, f2 = committed("blind:G2_bcat_central_d50")
    bx, bf = judged(47.0, 0.6)
    ok = rng.random(n) < q
    ax, af = committed("transfer:mab21l2_all")
    vx, vf = committed("transfer:mab21l2_vs95")
    for mask, (x, f) in ((u1 & anc & ~wc, (ax, af)), (wc, (np.where(ok, ax, vx), np.where(ok, af, vf))),
                         (u1 & ~anc, (vx, vf))):
        bx, bf = np.where(mask, x, bx), np.where(mask, f, bf)
    idx = {k: rng.choice(o["w"].size, size=n, p=o["w"]) for k, o in others.items()}
    rest_a = sum(E.OVERALL_WEIGHT[k] * o["a"][idx[k]] for k, o in others.items())
    rest_b = sum(E.OVERALL_WEIGHT[k] * o["b"][idx[k]] for k, o in others.items())
    eta = rng.normal(0.0, 3.0, n)
    res = {}
    for reading in ("half", "zero"):
        c = P.FLOOR_COST if reading == "zero" else 0.0
        a1, a2, ab = s1 - c * f1, s2 - c * f2, bx - c * bf
        t_today, t_swap = np.maximum(a1, a2), np.maximum(a1, ab)
        o_today, o_swap = rest_a + t_today, rest_a + t_swap
        s1tot = rest_b + a1
        shift = 1.0 * (s1tot - s1tot.mean()) + eta
        dp = {}
        for q_ in (90, 95):
            cl = float(np.percentile(o_today, q_))
            dp[q_] = float(np.mean(o_swap >= cl + shift) - np.mean(o_today >= cl + shift))
        res[reading] = {"dE": float(t_swap.mean() - t_today.mean()), "dP_lambda1": dp,
                        "wc_gain": float((t_swap - t_today)[wc].mean()), "p_wc": float(wc.mean())}
    mix = {"dE": (1 - P.P_ZERO) * res["half"]["dE"] + P.P_ZERO * res["zero"]["dE"],
           "dP_lambda1": {k: (1 - P.P_ZERO) * res["half"]["dP_lambda1"][k] + P.P_ZERO * res["zero"]["dP_lambda1"][k]
                          for k in (90, 95)}}
    return {"p_u2": p_u2, "p_anchor": p_anchor, "s1_wc": s1_wc, "q_barred": q, **res, "mixture": mix}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    door, checks = E.Door(REPO), []
    t1, _, _ = E.board_t1(door, checks)
    heart, _, _ = E.board_heart(door, checks)
    embryo = E.board_embryo(door, checks)
    rows = P.load_rows(door)
    others = {"T1": t1, "heart": heart, "embryo": embryo}
    cells = [simulate(rows, others, pu, pa, s, q, np.random.default_rng(P.SEED))
             for pu, pa in BANDS for s in S1_WC for q in Q_BARRED]
    out = {"script": "experiments/gate/t3_backup_sensitivity.py", "self_checks": len(checks), "cells": cells}
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"self-checks run: {len(checks)}")
    print("p_u2 p_anchor | s1 in that world (mean, P floor) | q | dE half / zero / mixture | dP mix lambda1 q90, q95")
    for c in cells:
        m = c["mixture"]
        print(f"{c['p_u2']:.2f} {c['p_anchor']:.2f} | {c['s1_wc'][0]:.0f}, {c['s1_wc'][1]:.1f} | {c['q_barred']:.2f} | "
              f"{c['half']['dE']:+.2f} / {c['zero']['dE']:+.2f} / {m['dE']:+.2f} | "
              f"{100 * m['dP_lambda1'][90]:+.2f}, {100 * m['dP_lambda1'][95]:+.2f}")


if __name__ == "__main__":
    main()
