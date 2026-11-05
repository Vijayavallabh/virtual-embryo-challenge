"""T3's second upload: today's refused-branch s1 (no factor) against the barred-branch s1 (the Mab21l2 factor).

SUPERSEDED (2026-11-05; results/t3_swap_review_2026-11-05.md, docs/PLAN.md decision log 2026-11-05 (d)): the worlds
below treat s1's and the swap's outcomes as independent, though both factors are referenced to WT E8.75 and share
its E9.5-ward offset, and they floor s1's failure at K2a. Coupled and unfloored, the swap's expected T3 gain is about
-0.2 to +0.1, and it can lower the chance of middle totals. Kept as the record of what was priced.

Written for the captain's aim of 2026-11-05 (the choices that maximise the chance of a top-3 place). Aggregation of
committed results with experiments/gate/estimate.py's loaders, not a try: no target is read and no score is computed.
The two candidate files differ only by the Mab21l2 factor (docs/S2_RULES.md, T3: the same wild-type cells and the same
central composition at the Mab21l2 magnitude); s1, the allowed branch's s1, is unchanged.

The unseen knockout (beta-catenin) is priced as a stated mixture of three worlds, each built only from committed
per-seed scores, with the weights and the failure size varied rather than fitted:
  A   s1's own within-knockout pattern carries (weight a): on each Gata4 replicate and seed, s1 = xrep:F+depth, today's
      s2 = blind:G2_bcat_central_d50 (its stand-in in estimate.py) and the swap = transfer:mab21l2_all (stand-in B).
  KB  s1's factor fails across knockouts as K2a measured, and the knockout-series factor carries as D11b measured
      (weight b): s1 = K2a on Mab21l2, today's s2 = our WT copy there, the swap = that copy plus B's paired gain over
      today's s2 stand-in on Gata4 (every K2 seed with every Gata4 pair).
  KF  both factors fail (weight c): as KB, but the swap = the copy plus delta, delta = -0.66 (K2a minus the copy, a
      failed factor measured on Mab21l2) or -4.1 (the WT-series factor R against the null on Mab21l2,
      experiments/t3_decomp/rules.yaml note 3).
Each board keeps the better of its two uploads. Overall draws one outcome per board (the boards independent, common
random numbers), and P(Overall >= c) is printed across our own range, so no outside threshold enters.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate as E  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
WEIGHTS = [(0.5, 0.25, 0.25), (0.5, 0.125, 0.375), (0.5, 0.05, 0.45), (0.25, 0.375, 0.375), (0.25, 0.1875, 0.5625),
           (0.25, 0.075, 0.675)]
DELTAS = (-0.66, -4.1)


def t3_options(rec, k2, a, b, c, delta):
    """Two outcome sets with the same structure and weights: today's pair and the swap."""
    cur, swp, w = [], [], []
    for r in E.T3_REPLICATES:  # world A
        s1, wt, mb = rec[r][E.T3_STANDINS["A"]], rec[r][E.T3_S2], rec[r][E.T3_STANDINS["B"]]
        cur.append(np.maximum(s1, wt))
        swp.append(np.maximum(s1, mb))
        w.append(np.full(s1.size, a / (len(E.T3_REPLICATES) * s1.size)))
    gain = np.concatenate([rec[r][E.T3_STANDINS["B"]] - rec[r][E.T3_S2] for r in E.T3_REPLICATES])
    s1k, copy = np.repeat(k2[E.T3_K2A], gain.size), np.repeat(k2[E.T3_K2_S2], gain.size)
    cur.append(np.maximum(s1k, copy))  # world KB
    swp.append(np.maximum(s1k, copy + np.tile(gain, k2[E.T3_K2A].size)))
    w.append(np.full(s1k.size, b / s1k.size))
    s1f, copyf = k2[E.T3_K2A], k2[E.T3_K2_S2]
    cur.append(np.maximum(s1f, copyf))  # world KF
    swp.append(np.maximum(s1f, copyf + delta))
    w.append(np.full(s1f.size, c / s1f.size))
    w_ = np.concatenate(w)
    return ({"a": np.concatenate(cur), "w": w_ / w_.sum()}, {"a": np.concatenate(swp), "w": w_ / w_.sum()},
            float(gain.mean()), int((gain > 0).sum()), int(gain.size))


def main() -> None:
    door, checks = E.Door(REPO), []
    t1, _, _ = E.board_t1(door, checks)
    heart, _, _ = E.board_heart(door, checks)
    embryo = E.board_embryo(door, checks)
    rec = E.load_t3(door, (E.T3_STANDINS["A"], E.T3_STANDINS["B"], E.T3_S2), checks)
    k2 = E.load_task(door, E.T3_K2, ("floor", E.T3_K2A, E.T3_K2_S2))
    others = {"T1": t1, "heart": heart, "embryo": embryo}
    print(f"self-checks run: {len(checks)}; K2a {np.mean(k2[E.T3_K2A]):.2f}, our WT copy on its fold "
          f"{np.mean(k2[E.T3_K2_S2]):.2f}")
    first = True
    for a, b, c in WEIGHTS:
        for delta in DELTAS:
            cur, swp, g, ahead, n = t3_options(rec, k2, a, b, c, delta)
            if first:
                print(f"B minus today's s2 stand-in on Gata4, paired: mean {g:+.2f}, ahead in {ahead} of {n}")
                print("weights a/b/c, delta | T3 today -> swap (gain) | Overall today -> swap | P(Overall >= c), %, "
                      "today -> swap, at today's 75th/90th/95th/99th percentiles | levels where the swap trails")
                first = False
            rng = np.random.default_rng(E.RNG_SEED)
            idx = {k: rng.choice(o["w"].size, size=E.N_DRAWS, p=o["w"]) for k, o in others.items()}
            it3 = rng.choice(cur["w"].size, size=E.N_DRAWS, p=cur["w"])
            base = sum(E.OVERALL_WEIGHT[k] * o["a"][idx[k]] for k, o in others.items())
            tot_c, tot_s = base + cur["a"][it3], base + swp["a"][it3]
            qs = [float(np.percentile(tot_c, q)) for q in (75, 90, 95, 99)]
            ps = " ".join(f"{100 * np.mean(tot_c >= q):.1f}->{100 * np.mean(tot_s >= q):.1f}" for q in qs)
            levels = np.arange(np.floor(np.percentile(tot_c, 50)), np.ceil(np.percentile(tot_c, 99.5)) + 1, 2.0)
            trails = int(sum(np.mean(tot_s >= x) < np.mean(tot_c >= x) - 1e-9 for x in levels))
            t3c, t3s = E.reading(cur, "a"), E.reading(swp, "a")
            print(f"{a:.2f}/{b:.3f}/{c:.3f}, {delta:+.2f} | {t3c:.2f} -> {t3s:.2f} ({t3s - t3c:+.2f}) | "
                  f"{float(np.mean(tot_c)):.2f} -> {float(np.mean(tot_s)):.2f} | {ps} | {trails} of {len(levels)}")


if __name__ == "__main__":
    main()
