"""Which second uploads give the best chance of a high total: every committed candidate pair, judged jointly.

SUPERSEDED in part (2026-11-05; docs/PLAN.md decision log 2026-11-05 (d)): its T3 options use t3_swap.py's worlds,
which the review refuted, and its T1 top-end gain from blood at 0.4 rests on the grid's design prior (a third on
blood near 40%), which no allowed data supports. Kept as the record of what was priced.

Written for the captain's aim of 2026-11-05 (the best chance of a top-3 place). Aggregation of committed results with
experiments/gate/estimate.py's loaders, not a try: no target is read and no score is computed. Each board keeps the
better of its two uploads; Overall = T1 + (heart + embryo) / 2 + T3, drawn 200,000 times with common random numbers
(one outcome per board, the boards independent). Levels are our own percentiles of today's pair (75th, 90th, 95th,
99th), so no outside threshold enters; which level matters is the captain's judgement.

Candidates (s1 fixed on every board):
  T1      s2 = any column of the oracle-base grid (s1 bown_cown), under each of the grid's five stated priors. A column
          that moves the blood share is a composition change: it enters only by the captain's ruling on a prior he
          writes (experiments/harness/TRIES.md, T1 item 1).
  heart   today's s2 only: its other committed candidates change the scope, which the post-exposure block closes
          (docs/DATA_POLICY.md 2026-11-04 (b)-(d)).
  embryo  s2 = fgw_shape0 (today), fgw_aligned (R2, emb_aligned_v1) or fgw_shape1, on the two counted backtests.
  T3      s2 = today's refused-branch s1 or the barred-branch s1 (experiments/gate/t3_swap.py's worlds; a = 0.5,
          b = 0.25, c = 0.25 and delta = -0.66, then the pessimistic a = 0.25, b = 0.075, c = 0.675, delta = -4.1).
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate as E  # noqa: E402
import t3_swap  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
LEVELS = (75, 90, 95, 99)
T3_MODELS = {"base": (0.5, 0.25, 0.25, -0.66), "pessimistic": (0.25, 0.075, 0.675, -4.1)}
EMBRYO_S2 = ("fgw_shape0", "fgw_aligned", "fgw_shape1")


def embryo_options(door):
    out = {k: [] for k in EMBRYO_S2}
    for shape_rel, _ in E.EMBRYO:
        s = E.load_task(door, shape_rel, ("floor", "fgw", "fgw_shape0", "fgw_shape1"))
        al = E.load_task(door, shape_rel.replace("emb_shape_v1", "emb_aligned_v1"), ("fgw_aligned",))
        if s["key"] != al["key"] or s["seeds"] != al["seeds"]:
            sys.exit(f"{shape_rel}: emb_aligned_v1 is not the same backtest and seeds")
        cols = {"fgw_shape0": s["fgw_shape0"], "fgw_shape1": s["fgw_shape1"], "fgw_aligned": al["fgw_aligned"]}
        for k in EMBRYO_S2:
            out[k].append((1.0, E.grid_outcomes({"b": s["fgw"]}, {"b": cols[k]}, {"b": 1.0})))
    return {k: E.mixture(v) for k, v in out.items()}


def main() -> None:
    door, checks = E.Door(REPO), []
    _, _, steps = E.board_t1(door, checks)
    heart, _, _ = E.board_heart(door, checks)
    emb = embryo_options(door)
    if not all(np.allclose(emb[k]["w"], emb["fgw_shape0"]["w"]) for k in EMBRYO_S2):
        sys.exit("embryo candidates are not paired with today's s2")
    rec = E.load_t3(door, (E.T3_STANDINS["A"], E.T3_STANDINS["B"], E.T3_S2), checks)
    k2 = E.load_task(door, E.T3_K2, ("floor", E.T3_K2A, E.T3_K2_S2))
    grid, hedge = E.load_runs(door, E.T1_GRID_RUNS), door.json(E.T1_GRID_HEDGE)
    sc = grid["scores"]
    t1_s2 = [p for p in grid["preds"] if p not in ("floor", E.T1_S1)]
    print(f"self-checks run: {len(checks)}; candidates: T1 {len(t1_s2)}, embryo {len(EMBRYO_S2)}, T3 2; heart fixed")
    today = (E.T1_S2, "fgw_shape0", "refused")
    for (mname, (a, b, c, delta)), (pname, prior) in itertools.product(T3_MODELS.items(), hedge["priors"].items()):
        cur3, swp3, *_ = t3_swap.t3_options(rec, k2, a, b, c, delta)
        t3 = {"refused": cur3, "barred": swp3}
        t1 = {p: E.grid_outcomes({t: sc[t][E.T1_S1] for t in prior}, {t: sc[t][p] for t in prior}, prior,
                                 shift=steps["m1_minus_b0"]) for p in t1_s2}
        rng = np.random.default_rng(E.RNG_SEED)
        i1 = rng.choice(t1[E.T1_S2]["w"].size, size=E.N_DRAWS, p=t1[E.T1_S2]["w"])
        ih = rng.choice(heart["w"].size, size=E.N_DRAWS, p=heart["w"])
        ie = rng.choice(emb["fgw_shape0"]["w"].size, size=E.N_DRAWS, p=emb["fgw_shape0"]["w"])
        i3 = rng.choice(cur3["w"].size, size=E.N_DRAWS, p=cur3["w"])
        h = E.OVERALL_WEIGHT["heart"] * heart["a"][ih]
        tot = {}
        for p, e, s in itertools.product(t1_s2, EMBRYO_S2, t3):
            tot[p, e, s] = (h + t1[p]["a"][i1] + E.OVERALL_WEIGHT["embryo"] * emb[e]["a"][ie] + t3[s]["a"][i3])
        lv = [float(np.percentile(tot[today], q)) for q in LEVELS]
        print(f"\nT3 model {mname}, T1 prior {pname}: today's pair, mean {float(np.mean(tot[today])):.2f}, "
              f"P at its 75/90/95/99th = 25/10/5/1 %")
        for q, x in zip(LEVELS, lv, strict=True):
            best = max(tot, key=lambda k: (np.mean(tot[k] >= x), np.mean(tot[k])))
            pb, pt = 100 * np.mean(tot[best] >= x), 100 * np.mean(tot[today] >= x)
            print(f"  level {q}th ({x:.1f}): best {best[0]} + {best[1]} + T3 {best[2]}: {pt:.1f} -> {pb:.1f} %, "
                  f"mean {float(np.mean(tot[best])):.2f}")


if __name__ == "__main__":
    main()
