"""How far today's entries sit from our own upper percentiles, and which board's outcome carries the upper tail.

Written for the state-of-the-art research of 2026-11-06 (experiments/research_2026_11_sota/REPORT.md, section 2).
Aggregation with experiments/gate/estimate.py's own loaders and outcome rules, not a try: no target, no score computed,
no leaderboard figure used to choose anything. The captain's 209 (his upload gate, the leading Human total when he
set it; docs/PLAN.md section 8) is reported once, for information only: it is never a threshold or evidence for a
try, an adoption or a method choice (docs/DATA_POLICY.md 2026-11-05). Options are compared at our own percentiles.

Reading A (each board at the better of its two uploads, with K2a, the record's reproducible reading) is drawn
200,000 times, one outcome per board, the boards independent, with the estimator's seed. Printed:
  1. P(Overall >= c) for c at our own 50th-99th percentiles (and, for information, at 209).
  2. Per board: its mean overall, and its mean in the draws that reach our own 95th percentile (where the upper tail
     comes from).
  3. What-if: a sure gain of g Overall points (any mix of boards), and P(Overall + g >= c) at our own 90th and 95th
     percentiles of today (and, for information, at 209). An expected gain applied as a sure shift.
  4. Sensitivity: the weight of T3's stand-in A (the Gata4 replicate fold, a fold the adoption rule never counts) at
     0.5 (the record), 0.375, 0.25 and 0, the rest split evenly between stand-in B and K2a as in the record, and
     P(Overall >= our own 90th / 95th percentile of today).
"""
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate as E  # noqa: E402

TARGET = 209.0  # the captain's upload gate (results/upload_gate.json), for information only
OWN_Q = (50, 75, 90, 95, 99)
OVERALL_GAINS = (0.35, 0.8, 1, 2, 3, 5, 8, 10, 16)  # sure Overall gains, incl. the research's +0.35 and +0.8
A_WEIGHTS = (0.5, 0.375, 0.25, 0.0)  # T3 stand-in A's weight; the record's is 0.5 (estimate.T3_WEIGHTS_K2A)


def t3_with_weight(door, checks, w_a: float) -> dict:
    """T3's outcomes as estimate.board_t3 builds them, with stand-in A at weight w_a and B and K2a at (1 - w_a)/2."""
    rec = E.load_t3(door, (E.T3_STANDINS["A"], E.T3_STANDINS["B"], E.T3_S2), checks)
    k2 = E.load_task(door, E.T3_K2, ("floor", E.T3_K2A, E.T3_K2_S2))
    standins = {name: {r: rec[r][v] for r in E.T3_REPLICATES} for name, v in E.T3_STANDINS.items()}
    standins["K2a"] = (k2[E.T3_K2A], k2[E.T3_K2_S2])
    s2 = {r: rec[r][E.T3_S2] for r in E.T3_REPLICATES}
    rest = (1.0 - w_a) / 2.0
    weights = {k: v for k, v in (("A", w_a), ("B", rest), ("K2a", rest)) if v > 0}
    return E.t3_outcomes(standins, s2, weights)


def draw_total(boards: dict, seed: int) -> tuple[np.ndarray, dict]:
    rng = np.random.default_rng(seed)
    draws = {}
    for b in E.BOARDS:
        o = boards[b]
        draws[b] = o["a"][rng.choice(o["w"].size, size=E.N_DRAWS, p=o["w"])]
    return sum(E.OVERALL_WEIGHT[b] * draws[b] for b in E.BOARDS), draws


def main() -> int:
    door, checks = E.Door(REPO), []
    t1, _, _ = E.board_t1(door, checks)
    heart, _, _ = E.board_heart(door, checks)
    embryo = E.board_embryo(door, checks)
    _, t3, _ = E.board_t3(door, checks)  # with K2a, as the record's reproducible reading
    boards = {"T1": t1, "heart": heart, "embryo": embryo, "T3": t3}
    if abs(E.overall({b: E.reading(o, "a") for b, o in boards.items()}) - 191.51) > 0.01:
        raise SystemExit("path_to_target: reading A does not reproduce the record's 191.51")
    E.same(E.reading(t3_with_weight(door, [], 0.5), "a"), E.reading(t3, "a"),
           "path_to_target: T3 rebuilt at stand-in A weight 0.5 reproduces board_t3", checks)
    tot, draws = draw_total(boards, E.RNG_SEED)
    own = {q: float(np.percentile(tot, q)) for q in OWN_Q}
    out = {"label": "aggregation of committed results, not a try; 209 for information only",
           "mean_overall": float(tot.mean()),
           "own_percentiles": {str(q): v for q, v in own.items()},
           "p_at_least_own": {str(q): float(np.mean(tot >= v)) for q, v in own.items()},
           "for_information_p_at_least_209": float(np.mean(tot >= TARGET)),
           "board_means": {b: float(draws[b].mean()) for b in E.BOARDS}}
    hit = tot >= own[95]
    out["board_means_given_own_q95"] = {b: float(draws[b][hit].mean()) for b in E.BOARDS}
    out["sure_overall_gain"] = {
        str(g): {"p_own_q90": float(np.mean(tot + g >= own[90])), "p_own_q95": float(np.mean(tot + g >= own[95])),
                 "for_information_p_209": float(np.mean(tot + g >= TARGET))}
        for g in OVERALL_GAINS}
    sens = {}
    for w in A_WEIGHTS:
        t3w = t3_with_weight(door, [], w)
        tw, _ = draw_total({**boards, "T3": t3w}, E.RNG_SEED)
        sens[str(w)] = {"T3_mean": E.reading(t3w, "a"), "overall_mean": float(tw.mean()),
                        "p_own_q90": float(np.mean(tw >= own[90])), "p_own_q95": float(np.mean(tw >= own[95]))}
    out["t3_standin_A_weight_sensitivity"] = sens
    out["checks"] = checks
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
