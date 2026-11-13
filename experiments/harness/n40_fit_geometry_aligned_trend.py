#!/usr/bin/env python
"""Reading N40 of experiments/harness/TRIES.md (registered before it was run; it can only withdraw): reading N39's
fit-geometry check, put to the other arm that passed reading N37 as a first upload, `fgw_aligned_trend`.

Why. Reading N37 found two arms that pass row N33's rule as the embryo's first upload: `fgw_trend_render` (request
R-D) and `fgw_aligned_trend`. Reading N39 showed that the first one's gain over `fgw` follows the trend line's fit
geometry and withdrew the AI assistant's yes. The decision log then said that nothing from the round is recommended,
without a reading of the second arm. An independent audit of the round named that gap. `fgw_aligned_trend` carries
the same trend rendering, at try E10's aligned positions, so the same question applies to it.

What it reads. Stage days and the nine practice tasks' definitions, which hold no data; and per-seed board scores
committed by counted tries: `fgw` and `fgw_aligned` from the stored runs of E24, E25 and E28, `fgw_aligned_trend` from
the committed task files of try E35 (given by --runs). It opens no data file, runs no candidate and no scorer: an
aggregation of committed scores, not a try.

Definitions: reading N39's (n39_fit_geometry_check.py), imported unchanged: the sample rule, the line's weights, the
mixture's shares, the truth-like sample of each task, g, the two fits (the seven tasks with a truth-like sample, and
the same without the E7.25 backtest) and the test's geometry. Two quantities are fitted on g:
  own    the arm's own mean board score minus fgw's own: reading N39's quantity, and the one this row's rule reads;
  trend  the arm's own mean minus fgw_aligned's own: the trend rendering's part at the aligned positions, the arm's
         one difference from an arm the ledger has already read (reported, deciding nothing).

Rule, fixed in the row before the run (reading N39's, with the arm changed). `fgw_aligned_trend` is recorded as not
put forward as a first upload if the fit's prediction of `own` at g_early (a test truth like its earlier bracket) is
at or below 0, in the fit over all seven tasks OR in the fit without the E7.25 backtest. Otherwise it stays on record
as an alternative that this lens does not withdraw, and goes to a cold adversarial check before any advice. The
reading cannot put the arm forward by itself.

  python -I -B experiments/harness/n40_fit_geometry_aligned_trend.py --runs DIR --out DIR
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n39_fit_geometry_check as n39  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec import harness  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ARM = "fgw_aligned_trend"
QUANTITIES = {"own": "fgw", "trend": "fgw_aligned"}  # the quantity's name -> the arm it is measured against
RULE_QUANTITY = "own"
PROMPTED = n39.PROMPTED


def lens(scores: dict, tasks: dict, arm: str, base: str) -> dict:
    """Reading N39's fits for one arm against one base: per task g and the paired mean gain, the two lines, and their
    values at the test's two geometries."""
    rows = {}
    for n, t in tasks.items():
        geo = n39.geometry(harness.STAGE_DAYS[t.inputs["before"]], harness.STAGE_DAYS[t.inputs["after"]],
                           harness.STAGE_DAYS[t.target], n39.LIKE[n])
        s = scores[t.key]
        rows[n] = {**geo, "y": float(s[arm].mean() - s[base].mean()), "ahead": int((s[arm] > s[base]).sum())}
    classed = [n for n in rows if rows[n]["like"] is not None]
    rest = [n for n in classed if n != PROMPTED]
    fits = {"all_seven": n39.fit([rows[n]["g"] for n in classed], [rows[n]["y"] for n in classed]),
            f"without_{PROMPTED}": n39.fit([rows[n]["g"] for n in rest], [rows[n]["y"] for n in rest])}
    te = n39.geometry(n39.TEST["before"], n39.TEST["after"], n39.TEST["target"], "before")
    tl = n39.geometry(n39.TEST["before"], n39.TEST["after"], n39.TEST["target"], "after")
    pred = {k: {"g_early": te["g"], "g_late": tl["g"], "y_early": f["a"] + f["b"] * te["g"],
                "y_late": f["a"] + f["b"] * tl["g"]} for k, f in fits.items()}
    return {"arm": arm, "base": base, "tasks": rows, "fits": fits, "prediction_at_the_test": pred,
            "test_geometry": {"days": te["days"], "weight_before": te["weight_before"],
                              "weight_after": te["weight_after"], "t": te["t"]}}


def check(scores: dict, tasks: dict) -> dict:
    out = {"arm": ARM, "quantities": {q: lens(scores, tasks, ARM, base) for q, base in QUANTITIES.items()}}
    p = out["quantities"][RULE_QUANTITY]["prediction_at_the_test"]
    out["not_put_forward"] = bool(p["all_seven"]["y_early"] <= 0 or p[f"without_{PROMPTED}"]["y_early"] <= 0)
    return out


def reading(run_dirs) -> dict:
    scores, why, meta = nt._load([*nt.STORED, n39.E28_RUNS, *run_dirs])
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (ARM, *QUANTITIES.values()):
            if arm not in scores.get(t.key, {}):
                why.append(f"{n}: no scores for {arm}")
    if set(tasks) != set(n39.LIKE):
        why.append("the nine tasks are not the ones this reading classes")
    if why:
        return {"void": True, "why": sorted(set(why))}
    return {"void": False, "why": [], "runs": meta, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    k = f"without_{PROMPTED}"
    lines = [f"Reading N40 (runs {', '.join(res['runs'])}): the fit geometry of {res['arm']} on the nine practice "
             "tasks and at the test; aggregation of committed scores and stage days, not a try",
             f"{res['arm']} as a first upload: " + ("NOT put forward (the fitted gain over fgw is at or below 0 at "
                                                   "the test's earlier-bracket geometry)"
                                                   if res["not_put_forward"]
                                                   else "stays on record as an alternative; a cold check comes next")]
    for q, v in res["quantities"].items():
        role = "the rule reads this one" if q == RULE_QUANTITY else "reported, deciding nothing"
        lines.append(f"  quantity '{q}': {v['arm']} minus {v['base']} ({role})")
        lines.append("    task   t      like     weight  share      g     gain   ahead")
        for n, r in v["tasks"].items():
            if r["like"] is None:
                mid = f"{'own':<7} {'':>6} {'':>6} {'':>6}"
            else:
                mid = f"{r['like']:<7} {r['weight_like']:6.3f} {r['share_like']:6.3f} {r['g']:+6.3f}"
            lines.append(f"    {n:<6} {r['t']:5.3f}  {mid}  {r['y']:+6.2f}  {r['ahead']:>2}/10")
        for name, key in (("all seven classed tasks", "all_seven"), (f"without {PROMPTED}", k)):
            f, p = v["fits"][key], v["prediction_at_the_test"][key]
            lines.append(f"    fit over {name}: gain = {f['a']:+.2f} {f['b']:+.2f} g  (r {f['r']:+.2f}, n {f['n']}); "
                         f"truth like E7.5: g {p['g_early']:+.3f} -> {p['y_early']:+.2f};  truth like E8.0: "
                         f"g {p['g_late']:+.3f} -> {p['y_late']:+.2f}")
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
    (out / "n40.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n40.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
