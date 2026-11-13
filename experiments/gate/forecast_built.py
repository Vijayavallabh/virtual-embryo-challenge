#!/usr/bin/env python3
"""The forecast restated for the files as built at the third rehearsal, and what it would take to pass a line.

Aggregation of committed practice results only, not a try: no target, no data file, no scorer, no build, no new score.
Everything is read through experiments/gate/estimate.py's door with forecast.py's own loaders, world table and draws;
this script adds no world, weight or judgement of its own. No test score and no leaderboard figure enters any number;
the comparison lines are forecast.py's and are applied to the finished draws only.

Three parts, each written to the JSON and printed:

1. FILES. forecast.py's files of 2026-11-13, then the same model with the two rulings that followed it that day:
   T1's second upload on the mid base (row R7: forecast_t1_pair.py's pair `mid_s2`, first upload on library 3 and
   second on the mid base at 60% cardiomyocytes), and the heart's trend (row H11: both heart uploads moved by the
   committed paired gain of `b0_trend` over `b0`, the mean of the two counted backtests, since both files carry it).
2. WHAT IT TAKES. For the files as built: each board's better-of-two outcome (mean, percentiles, maximum), the sum of
   the boards' maxima and of their 95th percentiles in Overall units, and P(total >= line) after a constant shift of
   one board's better-of-two outcome, in every world or inside one labelled world. A shift is not a method: it shows
   how large a gain on a board would have to be, and where a gain counts most.
3. A SECOND FILE CHOSEN AFTER THE FIRST SCORE. docs/S2_RULES.md fixes every second upload before any test score. This
   prices the alternative on committed practice results: the second file picked from the built menu by a threshold on
   the first upload's board total, the thresholds fitted on half of the seeds and read on the other half, with an
   unknown offset between practice and test score levels (one draw per board, SD 0, 2, 4 and 6 board points). T1: the
   oracle-base grid's 15 truths and 15 files. Embryo: the nine practice tasks of try E28 (means; leave-one-task-out).

Usage, from the repository root:
  uv run --frozen --no-sync python -I -B experiments/gate/forecast_built.py \
      --json experiments/gate/results/forecast_built_2026-11-13.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import estimate as E  # noqa: E402
import forecast as F  # noqa: E402
import forecast_t1_pair as FP  # noqa: E402

TREND = "experiments/harness/results/trend_heart_v1/tasks/{}.json"
TREND_TASKS = ("T2h_E9.5_E10.5", "T2h_E8.75_E10.5")
TREND_ROW_GAIN = {"T2h_E9.5_E10.5": 0.36, "T2h_E8.75_E10.5": 0.48}  # TRIES H11; checked against the files
BUILT_PAIR = ("lib3_base", "mid_cm60")  # rows R6 and R7
SHIFTS = (2.0, 5.0, 8.0, 10.0)
WORLD_SHIFTS = (("T3", "stand-in A", 5.0), ("T3", "stand-in A", 8.0), ("T3", "stand-in B", 8.0),
                ("T3", "stand-in K2a", 8.0), ("T1", "depth 1.37", 5.0), ("T1", "depth 1.0", 5.0),
                ("embryo", "D:", 5.0), ("embryo", "O:", 5.0))
JOINT = ({"T1": 5.0, "T3": 5.0, "heart": 5.0, "embryo": 5.0}, {"T1": 8.0, "T3": 8.0, "heart": 6.0, "embryo": 6.0},
         {"T1": 10.0, "T3": 10.0, "heart": 8.0, "embryo": 8.0})
OFFSET_SDS = (0.0, 2.0, 4.0, 6.0)
OFFSET_REPS = 400
OFFSET_SEED = 20261113
NINE = "experiments/harness/results/e28_nine_tasks/nine_tasks.json"
EMB_S1, EMB_S2 = "fgw", "fgw_aligned_series"


def wq(o: dict, qs=(0.05, 0.25, 0.5, 0.75, 0.95)) -> list[float]:
    """Weighted quantiles of a board's better-of-two outcome."""
    i = np.argsort(o["a"])
    cw = np.cumsum(o["w"][i])
    return [float(o["a"][i][min(np.searchsorted(cw, q), i.size - 1)]) for q in qs]


def stats(bd: dict) -> dict:
    tot, _ = F.draw(bd)
    return {"mean": float(tot.mean()), "median": float(np.percentile(tot, 50)),
            "percentiles": {str(q): float(np.percentile(tot, q)) for q in F.QS},
            "boards_better_of_two": {b: E.reading(o, "a") for b, o in bd.items()},
            "p_at_least": {f"{c:.2f}": float(np.mean(tot >= c)) for c in F.LINES}}


def shifted(bd: dict, moves: dict[str, float], world: str | None = None) -> dict:
    out = {b: dict(o) for b, o in bd.items()}
    for b, delta in moves.items():
        a = out[b]["a"].copy()
        m = np.ones(a.size, bool) if world is None else np.array([world in str(x) for x in out[b]["lab"]])
        F.need(bool(m.any()), f"{b}: a world labelled {world!r} exists")
        a[m] = np.minimum(100.0, a[m] + delta)
        out[b]["a"] = a
    return out


# =========================================================================================== part 1: the file sets
def file_sets(inp: F.Inputs) -> tuple[dict, dict]:
    par = F.central()
    bd = F.boards(inp, "new", par)
    arms = FP.load_arms(inp, ("s1_base", "s1_cm60", *BUILT_PAIR))
    t1 = FP.board_t1(inp, par["t1_prior"], par["t1_p_deep"], par["t1_k"], arms, BUILT_PAIR)
    gains = {}
    for task in TREND_TASKS:
        s = E.load_task(inp.door, TREND.format(task), ("b0", "b0_trend"))
        g = s["b0_trend"] - s["b0"]
        F.need(abs(float(g.mean()) - TREND_ROW_GAIN[task]) < 0.005, f"heart {task}: the trend's gain is row H11's",
               inp.checks)
        gains[task] = {"mean": float(g.mean()), "ahead": int((g > 0).sum()), "seeds": int(g.size)}
    shift = float(np.mean([v["mean"] for v in gains.values()]))
    heart = dict(bd["heart"])
    heart["a"], heart["b"] = heart["a"] + shift, heart["b"] + shift
    built = {**bd, "T1": t1, "heart": heart}
    sets = {"files of 2026-11-13 (forecast.py, committed)": stats(bd),
            "with T1's second upload on the mid base (R7)": stats({**bd, "T1": t1}),
            "the files as built: R7 and the heart trend (H11)": stats(built)}
    return {"sets": sets, "heart_trend_gain": gains, "heart_shift": shift, "t1_pair": list(BUILT_PAIR)}, built


# ====================================================================================== part 2: what it would take
def what_it_takes(built: dict) -> dict:
    per = {b: {"mean": E.reading(o, "a"), "q05_q25_q50_q75_q95": wq(o), "max": float(o["a"].max())}
           for b, o in built.items()}
    out = {"boards": per,
           "sum_of_maxima": float(sum(E.OVERALL_WEIGHT[b] * built[b]["a"].max() for b in built)),
           "sum_of_95th": float(sum(E.OVERALL_WEIGHT[b] * wq(built[b], (0.95,))[0] for b in built)),
           "board_shift": {}, "world_shift": {}, "joint_shift": {}}
    for b in E.BOARDS:
        for d in SHIFTS:
            out["board_shift"][f"{b} +{d:g}"] = stats(shifted(built, {b: d}))
    for b, w, d in WORLD_SHIFTS:
        out["world_shift"][f"{b} [{w}] +{d:g}"] = stats(shifted(built, {b: d}, w))
    for moves in JOINT:
        out["joint_shift"][", ".join(f"{b} +{d:g}" for b, d in moves.items())] = stats(shifted(built, moves))
    return out


# ======================================================================== part 3: a second file chosen after a score
def _fit(x: np.ndarray, pay: np.ndarray, k: int, default: int) -> tuple[tuple[float, ...], list[int]]:
    """The best k-bin threshold rule: x the first upload's score per case, pay[case, file] the better of two."""
    cand = np.unique(np.quantile(x, np.linspace(0.1, 0.9, 17)))
    best = (-np.inf, (), [default] * k)
    for cuts in itertools.combinations(cand, k - 1):
        edges, tot, files = [-np.inf, *cuts, np.inf], 0.0, []
        for a, b in itertools.pairwise(edges):
            m = (x > a) & (x <= b)
            j = int(pay[m].sum(0).argmax()) if m.any() else default
            files.append(j)
            tot += float(pay[m, j].sum())
        if tot > best[0]:
            best = (tot, tuple(float(c) for c in cuts), files)
    return best[1], best[2]


def _read(x: np.ndarray, pay: np.ndarray, cuts, files, sd: float, rng) -> float:
    vals = []
    for _ in range(OFFSET_REPS if sd > 0 else 1):
        off = rng.normal(0.0, sd) if sd > 0 else 0.0  # one unknown level offset for the whole board
        choice = np.array(files)[np.searchsorted(np.array(cuts), x + off)]
        vals.append(float(pay[np.arange(x.size), choice].mean()))
    return float(np.mean(vals))


def second_file_t1(inp: F.Inputs) -> dict:
    sc, seeds = inp.grid["scores"], inp.grid["seeds"]
    truths = sorted(sc)
    files = sorted(p for p in sc[truths[0]] if p != "floor")
    S = np.array([[np.asarray(sc[t][p], float) for p in files] for t in truths])  # truth x file x seed
    i1, i2 = files.index(E.T1_S1), files.index(E.T1_S2)
    fixed = float(np.maximum(S[:, i1], S[:, i2]).mean())
    F.need(abs(fixed - 68.28) < 0.01, "T1 grid: the fixed pair's better of two is the estimator's 68.28", inp.checks)
    halves = ([i for i, s in enumerate(seeds) if s % 2 == 0], [i for i, s in enumerate(seeds) if s % 2 == 1])
    rng = np.random.default_rng(OFFSET_SEED)
    out = {"fixed_pair": fixed, "first_alone": float(S[:, i1].mean()),
           "oracle_file_per_truth": float(np.maximum(S[:, i1], S.max(1)).mean()), "rules": {}}
    for k in (2, 3):
        gain = {}
        for sd in OFFSET_SDS:
            g = []
            for tr, te in (halves, halves[::-1]):
                x_tr, x_te = S[:, i1][:, tr].ravel(), S[:, i1][:, te].ravel()
                pay = lambda idx: np.maximum(S[:, i1][:, None, idx], S[:, :, idx]).transpose(0, 2, 1).reshape(  # noqa: E731
                    -1, len(files))
                cuts, chosen = _fit(x_tr, pay(tr), k, i2)
                base = float(np.maximum(S[:, i1][:, te], S[:, i2][:, te]).mean())
                g.append(_read(x_te, pay(te), cuts, chosen, sd, rng) - base)
            gain[f"sd {sd:g}"] = float(np.mean(g))
        allx = S[:, i1].ravel()
        cuts, chosen = _fit(allx, np.maximum(S[:, i1][:, None, :], S).transpose(0, 2, 1).reshape(-1, len(files)), k, i2)
        out["rules"][f"{k} bins"] = {"cuts_all_seeds": list(cuts), "files_all_seeds": [files[j] for j in chosen],
                                     "held_out_gain_over_fixed_pair": gain}
    return out


def second_file_embryo(inp: F.Inputs) -> dict:
    single = inp.door.json(NINE)["single"]
    tasks = sorted(single[EMB_S1])
    files = sorted(f for f in single if f != EMB_S1)
    x = np.array([single[EMB_S1][t] for t in tasks])
    pay = np.array([[max(single[EMB_S1][t], single[f][t]) for f in files] for t in tasks])
    i2 = files.index(EMB_S2)
    rng = np.random.default_rng(OFFSET_SEED)
    out = {"fixed_pair_mean9": float(pay[:, i2].mean()), "oracle_file_per_task": float(pay.max(1).mean()), "rules": {}}
    gain = {}
    for sd in OFFSET_SDS:
        g = []
        for hold in range(len(tasks)):  # leave one task out: fit the cut on eight tasks, read it on the ninth
            m = np.arange(len(tasks)) != hold
            cuts, chosen = _fit(x[m], pay[m], 2, i2)
            g.append(_read(x[[hold]], pay[[hold]], cuts, chosen, sd, rng) - float(pay[hold, i2]))
        gain[f"sd {sd:g}"] = float(np.mean(g))
    cuts, chosen = _fit(x, pay, 2, i2)
    out["rules"]["2 bins"] = {"cuts_all_tasks": list(cuts), "files_all_tasks": [files[j] for j in chosen],
                              "leave_one_task_out_gain_over_fixed_pair": gain}
    return out


# =============================================================================================================== main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", metavar="PATH")
    a = ap.parse_args(argv)
    inp = F.Inputs()
    part1, built = file_sets(inp)
    part2 = what_it_takes(built)
    part3 = {"T1": second_file_t1(inp), "embryo": second_file_embryo(inp)}
    lines = [f"{c:.2f}" for c in F.LINES]

    def row(name: str, s: dict) -> str:
        return (f"  {name:<52} {s['mean']:7.2f} {s['median']:7.2f} "
                + " ".join(f"{100 * s['p_at_least'][c]:6.2f}%" for c in lines))

    print("label: aggregation of committed practice results, not a try; lines applied to the finished draws only")
    print(f"  {'':<52} {'mean':>7} {'median':>7} " + " ".join(f"{'>=' + c:>7}" for c in lines))
    print("1. file sets")
    for n, s in part1["sets"].items():
        print(row(n, s))
        print("      boards, better of two: " + ", ".join(f"{b} {v:.2f}" for b, v in s["boards_better_of_two"].items()))
    print("2. the files as built: each board's better of two (mean; 5/25/50/75/95th percentiles; maximum)")
    for b, v in part2["boards"].items():
        print(f"  {b:<7} {v['mean']:6.2f}; " + " / ".join(f"{x:.1f}" for x in v["q05_q25_q50_q75_q95"])
              + f"; {v['max']:.1f}")
    print(f"  sum of the boards' maxima in Overall units {part2['sum_of_maxima']:.1f}; of their 95th percentiles "
          f"{part2['sum_of_95th']:.1f}")
    for key in ("board_shift", "world_shift", "joint_shift"):
        print(f"  {key.replace('_', ' ')} (a constant added to the better of two; not a method)")
        for n, s in part2[key].items():
            print(row("  " + n, s))
    print("3. a second file chosen by a threshold on the first upload's board total (held-out; gain over the fixed "
          "pair by the SD of an unknown practice-to-test level offset)")
    t1 = part3["T1"]
    print(f"  T1 grid: fixed pair {t1['fixed_pair']:.2f}, first upload alone {t1['first_alone']:.2f}, the best file "
          f"per truth {t1['oracle_file_per_truth']:.2f}")
    for n, r in t1["rules"].items():
        print(f"    {n}: " + ", ".join(f"{k} {v:+.2f}" for k, v in r["held_out_gain_over_fixed_pair"].items())
              + f"   (all seeds: cuts {[round(c, 1) for c in r['cuts_all_seeds']]}, files {r['files_all_seeds']})")
    em = part3["embryo"]
    print(f"  embryo, nine practice tasks: fixed pair {em['fixed_pair_mean9']:.2f}, the best file per task "
          f"{em['oracle_file_per_task']:.2f}")
    for n, r in em["rules"].items():
        loo = r["leave_one_task_out_gain_over_fixed_pair"]
        print(f"    {n}: " + ", ".join(f"{k} {v:+.2f}" for k, v in loo.items()))
    res = {"label": "the forecast restated for the files as built at the third rehearsal (rows R6, R7, H11), what a "
                    "constant shift per board would do, and a second file chosen after the first score, priced on "
                    "practice results: aggregation of committed results, not a try; no test score or leaderboard "
                    "figure enters any number; the comparison lines are applied to the finished draws only",
           "draws": F.N_DRAWS, "rng_seed": F.SEED, "lines": {f"{c:.2f}": t for c, t in F.LINES.items()},
           "parameters": F.central(), "files": part1, "what_it_takes": part2, "second_file_after_a_score": part3,
           "checks": sorted(set(inp.checks)), "inputs_sha256": dict(sorted(inp.door.read.items()))}
    if a.json:
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.json).write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
