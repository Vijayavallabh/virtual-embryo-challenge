#!/usr/bin/env python
"""Reading N39 of experiments/harness/TRIES.md (registered before it was run; it can only withdraw): does the gain
that reading N37 credits to `fgw_trend_render` carry to the test's geometry?

Why. A cold adversarial check of request R-D held that the arm's gain over `fgw` on a practice task follows the
task's fit geometry and not its label. The trend rendering moves every cell to the profile a least-squares line
through the task's samples gives at the target's day. That profile is a weighted sum of the samples' own profiles,
with the line's weights; the plain mixture `fgw` holds the two brackets' profiles in the shares (1 - t) and t. Where
the truth is like one sample, the arm should gain if the line gives that sample more weight than the mixture gives
it, and lose if it gives it less. On the task that carries most of N37's value (the E7.25 backtest) the line gives
the truth-like sample more; at the test, in the world that task stands for, it gives it less.

What it reads. Stage days (vec.harness.STAGE_DAYS) and the nine practice tasks' definitions, which hold no data; and
per-seed board scores committed by counted tries: `fgw` from the stored runs of E24, E25 and E28 and
`fgw_trend_render` from the committed task files of try E35 (given by --runs). It opens no data file, runs no
candidate and no scorer: an aggregation of committed scores, not a try.

Definitions, fixed before the run.
  samples   a task's trend samples by row E35's rule: its two brackets plus every stage of E6.75, E7.25, E7.5, E8.0
            strictly outside the bracket interval.
  weights   the least-squares line's weights at the target's day x over the sample days x_i:
            w_i = 1/n + (x - mean)(x_i - mean)/sum((x_j - mean)^2). They sum to 1.
  t         the target's fraction between its brackets; the mixture's shares are 1 - t (earlier) and t (later).
  like      the sample the record says the truth is like (nine_tasks.py, diagnostic N23): the later bracket on half,
            late and new1; the earlier bracket on x1, x2, x3 and on the E7.25 backtest; none on new2 and on the E7.5
            backtest, whose truths have a profile of their own (they are reported and not fitted).
  g         weight(like) - share(like).
  y         the arm's own mean board score minus fgw's own, over the 10 seeds.
  fit       y = a + b g by least squares over the seven tasks that have a truth-like sample; and the same without
            the E7.25 backtest, the task that prompted the request.
  the test  brackets E7.5 and E8.0, target day 7.75, samples E6.75, E7.25, E7.5, E8.0, t = 1/2. Two worlds with a
            truth-like sample: like its earlier bracket (g_early = w(E7.5) - 1/2) and like its later bracket
            (g_late = w(E8.0) - 1/2).

Rule, fixed in the row before the run. The AI assistant's advice on request R-D changes from yes to no (the request
stays on record as reading N37's outcome) if the fit's prediction of y at g_early is at or below 0 in the fit over
all seven tasks OR in the fit without the E7.25 backtest. Otherwise the advice stands. The reading cannot strengthen
the request. Reported, deciding nothing: the correlation of y with g, the prediction at g_late, and N37's S(0.45)
recomputed with D set to 0 (the tasks arranged like the test, half, late and new1, give 0.00) with and without the
E7.25 backtest.

  python -I -B experiments/harness/n39_fit_geometry_check.py --runs DIR --out DIR
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nine_tasks as nt  # noqa: E402

from vec import harness  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ARM, BASE, SECOND = "fgw_trend_render", "fgw", nt.RA
SERIES = ("T2/E6.75.h5ad", "T2/E7.25.h5ad", "T2/E7.5.h5ad", "T2/E8.0.h5ad")
LIKE = {"half": "after", "late": "after", "new1": "after", "x1": "before", "x2": "before", "x3": "before",
        "e725": "before", "new2": None, "e75": None}
PROMPTED = "e725"  # the E7.25 backtest, where try E33's one-sided gain prompted the request
TEST = {"before": 7.5, "after": 8.0, "target": 7.75}
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
P_MID = 0.45


def line_weights(days, x: float) -> np.ndarray:
    """The least-squares line's weights at x over the sample days (they sum to 1)."""
    d = np.asarray(days, float)
    c = d - d.mean()
    sxx = float((c ** 2).sum())
    if sxx <= 0:
        raise ValueError("the samples share one day: no line")
    return 1.0 / d.size + (x - d.mean()) * c / sxx


def sample_days(before: float, after: float) -> list[float]:
    """Row E35's sample rule on days: the brackets plus the series' stages strictly outside the bracket interval."""
    series = [harness.STAGE_DAYS[f] for f in SERIES]
    return sorted({before, after, *[s for s in series if s < before or s > after]})


def geometry(before: float, after: float, target: float, like: str | None) -> dict:
    days = sample_days(before, after)
    w = line_weights(days, target)
    t = (target - before) / (after - before)
    out = {"days": days, "t": t, "weight_before": float(w[days.index(before)]),
           "weight_after": float(w[days.index(after)]), "like": like}
    if like is not None:
        weight = out[f"weight_{like}"]
        share = t if like == "after" else 1.0 - t
        out.update(weight_like=weight, share_like=share, g=weight - share)
    return out


def fit(g, y) -> dict:
    g, y = np.asarray(g, float), np.asarray(y, float)
    b, a = np.polyfit(g, y, 1)
    r = float(np.corrcoef(g, y)[0, 1]) if g.size > 2 and g.std() > 0 and y.std() > 0 else float("nan")
    return {"a": float(a), "b": float(b), "r": r, "n": int(g.size)}


def check(scores: dict, tasks: dict) -> dict:
    rows = {}
    for n, t in tasks.items():
        geo = geometry(harness.STAGE_DAYS[t.inputs["before"]], harness.STAGE_DAYS[t.inputs["after"]],
                       harness.STAGE_DAYS[t.target], LIKE[n])
        s = scores[t.key]
        own = float(s[ARM].mean() - s[BASE].mean())
        pair = float(np.maximum(s[ARM], s[SECOND]).mean() - np.maximum(s[BASE], s[SECOND]).mean())
        rows[n] = {**geo, "y": own, "ahead": int((s[ARM] > s[BASE]).sum()), "pair_value": pair}
    classed = [n for n in rows if rows[n]["like"] is not None]
    full = fit([rows[n]["g"] for n in classed], [rows[n]["y"] for n in classed])
    rest = [n for n in classed if n != PROMPTED]
    without = fit([rows[n]["g"] for n in rest], [rows[n]["y"] for n in rest])
    te = geometry(TEST["before"], TEST["after"], TEST["target"], "before")
    tl = geometry(TEST["before"], TEST["after"], TEST["target"], "after")
    pred = {k: {"g_early": te["g"], "g_late": tl["g"], "y_early": f["a"] + f["b"] * te["g"],
                "y_late": f["a"] + f["b"] * tl["g"]}
            for k, f in (("all_seven", full), (f"without_{PROMPTED}", without))}
    other = [rows[n]["pair_value"] for n in nt.OTHER]
    o_all = float(np.mean(other))
    o_rest = float(np.mean([rows[n]["pair_value"] for n in nt.OTHER if n != PROMPTED]))
    withdraw = bool(pred["all_seven"]["y_early"] <= 0 or pred[f"without_{PROMPTED}"]["y_early"] <= 0)
    return {"arm": ARM, "tasks": rows, "fit_all_seven": full, f"fit_without_{PROMPTED}": without,
            "test_geometry": {"days": te["days"], "weight_before": te["weight_before"],
                              "weight_after": te["weight_after"], "t": te["t"]},
            "prediction_at_the_test": pred,
            "n37_with_D_zero": {"S0.45": (1 - P_MID) * o_all, f"S0.45_without_{PROMPTED}": (1 - P_MID) * o_rest},
            "advice_changes_to_no": withdraw}


def reading(run_dirs) -> dict:
    scores, why, meta = nt._load([*nt.STORED, E28_RUNS, *run_dirs])
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (ARM, BASE, SECOND):
            if arm not in scores.get(t.key, {}):
                why.append(f"{n}: no scores for {arm}")
    if set(tasks) != set(LIKE):
        why.append("the nine tasks are not the ones this reading classes")
    if why:
        return {"void": True, "why": sorted(set(why))}
    return {"void": False, "why": [], "runs": meta, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    k = f"without_{PROMPTED}"
    lines = [f"Reading N39 (runs {', '.join(res['runs'])}): the fit geometry of {res['arm']} on the nine practice "
             "tasks and at the test; aggregation of committed scores and stage days, not a try",
             "advice on request R-D: " + ("changes to NO (the request stays on record as reading N37's outcome)"
                                          if res["advice_changes_to_no"] else "stands"),
             "  task   samples (days)                          t   like     weight  share      g    own-fgw  ahead"
             "  pair value"]
    for n, r in res["tasks"].items():
        days = " ".join(f"{d:g}" for d in r["days"])
        if r["like"] is None:
            mid = f"{'own':<7} {'':>6} {'':>6} {'':>6}"
        else:
            mid = f"{r['like']:<7} {r['weight_like']:6.3f} {r['share_like']:6.3f} {r['g']:+6.3f}"
        lines.append(f"  {n:<6} {days:<34} {r['t']:5.3f}  {mid}  {r['y']:+7.2f}  {r['ahead']:>2}/10  "
                     f"{r['pair_value']:+7.2f}")
    for name, f in (("all seven classed tasks", res["fit_all_seven"]), (f"without {PROMPTED}", res[f"fit_{k}"])):
        lines.append(f"  fit over {name}: own-fgw = {f['a']:+.2f} {f['b']:+.2f} g  (r {f['r']:+.2f}, n {f['n']})")
    tg, p = res["test_geometry"], res["prediction_at_the_test"]
    lines.append(f"  the test: samples {' '.join(f'{d:g}' for d in tg['days'])}, t {tg['t']:.2f}; the line's weight on "
                 f"E7.5 {tg['weight_before']:.3f} and on E8.0 {tg['weight_after']:.3f}")
    for name, q in p.items():
        lines.append(f"    {name:<16} truth like E7.5: g {q['g_early']:+.3f} -> own-fgw {q['y_early']:+.2f};   "
                     f"truth like E8.0: g {q['g_late']:+.3f} -> own-fgw {q['y_late']:+.2f}")
    z = res["n37_with_D_zero"]
    lines.append(f"  N37's S(0.45) with D set to 0: {z['S0.45']:+.2f}; and without {PROMPTED}: "
                 f"{z[f'S0.45_without_{PROMPTED}']:+.2f}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of try E35")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n39.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n39.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
