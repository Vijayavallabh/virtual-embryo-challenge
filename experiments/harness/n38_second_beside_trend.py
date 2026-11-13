#!/usr/bin/env python
"""Reading N38 of experiments/harness/TRIES.md (registered before it was run): the embryo's second upload beside the
first upload that reading N37 put to the captain.

Why. Reading N37 found that `fgw_trend_render` passes row N33's rule as a first upload beside the second upload as
built. docs/S2_RULES.md sends the second upload back to the captain whenever the first one changes, so the request
has to say which second upload goes with the new first one. Every second upload was chosen beside `fgw`; no table
prices them beside the trend rendering.

What it reads. Only per-seed board scores committed by counted tries: the stored runs of E24, E25 and E28
(`nine_tasks.STORED` and `results/e28_nine_tasks/tasks`) and the committed task files of tries E34 and E35 given by
`--runs`. It opens no data file, runs no candidate and no scorer: an aggregation of committed scores, not a try.

Value of a second upload X beside the first upload FIRST: on each task, the mean over seeds of max(FIRST, X) minus
the same for max(FIRST, fgw_aligned_series), the second upload as built. D, O and S(p) as in nine_tasks.py (DEEP,
OTHER, PS); mean9; worst9 = the lowest task mean of the pair's better of two, in board points; and S(0.45) with each
task left out in turn.

Rule, fixed in the row before the reading (rows E34 and E35's, with the first upload changed). Beside FIRST, another
second upload replaces fgw_aligned_series in the advice only if all hold:
  (a) S(0.45) > BAR (0.3);
  (b) S(0.3) >= 0 and S(0.6) >= 0;
  (c) its pair's worst9 is not lower than the pair (FIRST, fgw_aligned_series)'s by more than TOL;
  (d) S(0.45) stays above 0 whichever single task is left out.
If several pass, the one with the largest S(0.45). Otherwise the advice is to keep fgw_aligned_series as the second
upload. An arm without scores on every one of the nine tasks and seeds is reported as not readable and cannot pass.
Nothing changes an upload without the captain's ruling.

  python -I -B experiments/harness/n38_second_beside_trend.py --runs DIR [--runs DIR ...] --out DIR
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

from vec.provenance import repo_state  # noqa: E402

FIRST = "fgw_trend_render"  # the first upload reading N37 put to the captain
BASE_SECOND = nt.RA  # fgw_aligned_series, the second upload as built
BUILT_FIRST = "fgw"  # the first upload as built
SECONDS = ("fgw", "fgw_aligned", "fgw_shape0", "fgw_render_hi", "fgw_shape0_hi", "fgw_series", "fgw_shape0_series",
           "fgw_aligned_hi", "fgw_early_specialist", "fgw_aligned_trend")
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
LEAVE_OUT_P = 0.45


def _has(scores: dict, tasks: dict, arm: str) -> list[str]:
    """The tasks on which an arm has no complete, finite per-seed scores."""
    return [n for n in (*nt.DEEP, *nt.OTHER)
            if arm not in scores[tasks[n].key] or not np.isfinite(scores[tasks[n].key][arm]).all()]


def price(scores: dict, tasks: dict, first: str = FIRST, candidates=SECONDS) -> dict:
    """The rule on per-seed scores: scores[task key][arm] -> array over the seeds."""
    names = [*nt.DEEP, *nt.OTHER]

    def pair_mean(a: str, b: str, n: str) -> float:
        s = scores[tasks[n].key]
        return float(np.maximum(s[a], s[b]).mean())

    base_abs = {n: pair_mean(first, BASE_SECOND, n) for n in names}
    built_abs = {n: pair_mean(BUILT_FIRST, BASE_SECOND, n) for n in names}
    out = {"first": first, "base_second": BASE_SECOND, "base_pair_task_means": base_abs,
           "base_pair_worst9": min(base_abs.values()),
           "base_pair_against_the_pair_as_built": {n: base_abs[n] - built_abs[n] for n in names},
           "second_uploads": {}}
    passed = []
    for x in candidates:
        missing = _has(scores, tasks, x)
        if missing:
            out["second_uploads"][x] = {"readable": False, "missing_on": missing, "passes": False}
            continue
        absx = {n: pair_mean(first, x, n) for n in names}
        per = {n: absx[n] - base_abs[n] for n in names}
        d = float(np.mean([per[n] for n in nt.DEEP]))
        o = float(np.mean([per[n] for n in nt.OTHER]))
        s = {p: p * d + (1 - p) * o for p in nt.PS}
        loo = {}
        for drop in names:
            dd = [per[n] for n in nt.DEEP if n != drop]
            oo = [per[n] for n in nt.OTHER if n != drop]
            loo[drop] = LEAVE_OUT_P * float(np.mean(dd)) + (1 - LEAVE_OUT_P) * float(np.mean(oo))
        worst = min(absx.values())
        rule = {"a_S0.45_above_bar": bool(s[0.45] > nt.BAR),
                "b_S0.3_and_S0.6_not_below_0": bool(s[0.3] >= 0 and s[0.6] >= 0),
                "c_worst9_not_lower": bool(worst >= out["base_pair_worst9"] - nt.TOL),
                "d_leave_one_task_out_positive": bool(min(loo.values()) > 0)}
        ok = all(rule.values())
        if ok:
            passed.append(x)
        out["second_uploads"][x] = {"readable": True, "tasks": per, "D": d, "O": o,
                                    **{f"S{p}": v for p, v in s.items()},
                                    "mean9": float(np.mean(list(per.values()))), "worst9_abs": worst,
                                    "S0.45_leave_one_task_out_min": min(loo.values()),
                                    "S0.45_leave_one_task_out": loo, "rule": rule, "passes": ok}
    out["advised_second"] = (max(passed, key=lambda x: out["second_uploads"][x]["S0.45"]) if passed
                             else BASE_SECOND)
    return out


def reading(run_dirs) -> dict:
    scores, why, meta = nt._load([*nt.STORED, E28_RUNS, *run_dirs])
    tasks = nt.all_tasks()
    for arm in (FIRST, BASE_SECOND, BUILT_FIRST):
        for n in tasks:
            if tasks[n].key not in scores:
                why.append(f"{n}: no scores at all")
            elif arm not in scores[tasks[n].key]:
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": sorted(set(why))}
    return {"void": False, "why": [], "runs": meta, **price(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    lines = [f"Reading N38 (runs {', '.join(res['runs'])}): second uploads beside the first upload {res['first']}, "
             f"against the pair ({res['first']}, {res['base_second']}); aggregation of committed scores, not a try",
             f"advised second upload: {res['advised_second']}"
             + ("" if res["advised_second"] != res["base_second"] else " (the second upload as built stays)"),
             "  second upload          " + " ".join(f"{n:>6}" for n in names)
             + "  |      D      O  " + " ".join(f"S({p})" for p in nt.PS) + "   mean9  worst9  LOO-min  rule a b c d"]
    for x, v in res["second_uploads"].items():
        if not v["readable"]:
            lines.append(f"  {x:<22} not readable: no scores on {', '.join(v['missing_on'])}")
            continue
        r = v["rule"]
        lines.append(f"  {x:<22} " + " ".join(f"{v['tasks'][n]:+6.2f}" for n in names)
                     + f"  | {v['D']:+6.2f} {v['O']:+6.2f}  " + "  ".join(f"{v[f'S{p}']:+5.2f}" for p in nt.PS)
                     + f"   {v['mean9']:+5.2f}  {v['worst9_abs']:6.2f}  {v['S0.45_leave_one_task_out_min']:+7.2f}  "
                     + " ".join("y" if r[k] else "n" for k in sorted(r)))
    lines.append(f"  the pair ({res['first']}, {res['base_second']}): worst9 {res['base_pair_worst9']:.2f}; task means "
                 + " ".join(f"{res['base_pair_task_means'][n]:.2f}" for n in names))
    lines.append("  that pair against the pair as built:                    "
                 + " ".join(f"{res['base_pair_against_the_pair_as_built'][n]:+.2f}" for n in names))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", action="append", required=True, help="a committed tasks folder of try E34 or E35")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n38.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n38.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
