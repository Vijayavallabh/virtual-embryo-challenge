#!/usr/bin/env python
"""Reading N37 of experiments/harness/TRIES.md (registered before any score it uses was read): the embryo arms that
tries E34 and E35 score on the nine practice tasks, priced as a FIRST upload beside the second upload as built.

Why. Try E33 found the trend rendering one-sided on the two counted backtests: a large gain where the truth is like
its earlier, shallower bracket and a loss on the other. Rows E34 and E35 read such arms as second uploads, where they
would replace `fgw_aligned_series`, the file that pays where the hidden stage shares its deeper bracket's profile. The
other placement is the first upload: beside `fgw_aligned_series` a one-sided arm only has to beat `fgw` where the
second upload does not already lead. Row N33 priced the older arms this way; this script applies N33's rule, unchanged,
to the new ones.

What it reads. Only per-seed board scores committed by counted tries: the stored runs of E24, E25 and E28
(`nine_tasks.STORED` and `results/e28_nine_tasks/tasks`) and the committed task files of tries E34 and E35 given by
`--runs` (their runs score the new arms on the nine tasks; the task keys and seeds are the stored runs' own). It opens
no data file, runs no candidate and no scorer: an aggregation of committed scores, not a try.

Value of a first upload X: on each task, the mean over seeds of max(X, second) minus the same for max(fgw, second),
with second = fgw_aligned_series; D, O and S(p) as in nine_tasks.py (DEEP, OTHER, PS); mean9; worst9 = the lowest task
mean of the pair's better of two, in board points; and S(0.45) with each task left out in turn.

Rule, fixed in the row before the reading (row N33's, word for word). A first upload other than fgw is put to the
captain only if all hold:
  (a) S(0.45) > BAR (0.3);
  (b) S(0.3) >= 0 and S(0.6) >= 0;
  (c) its pair's worst9 is not lower than fgw's pair's worst9 by more than TOL;
  (d) S(0.45) stays above 0 whichever single task is left out.
If several pass, the one with the largest S(0.45). Otherwise nothing is requested and the first upload stays fgw. An
arm without scores on every one of the nine tasks and seeds is reported as not readable and can never pass. A first
upload that passes still fails the adoption rule on the two counted backtests (try E33's veto), so it could enter by
the captain's ruling only, and docs/S2_RULES.md then sends the second upload back to the captain.

  python -I -B experiments/harness/n37_first_upload_new_arms.py --runs DIR [--runs DIR ...] --out DIR
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

SECOND = nt.RA  # fgw_aligned_series, the second upload since row R4
BASE = "fgw"  # the first upload as built
NEW = ("fgw_trend_render", "fgw_aligned_trend", "fgw_aligned_hi", "fgw_early_specialist")
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
LEAVE_OUT_P = 0.45


def price(scores: dict, tasks: dict, candidates=NEW) -> dict:
    """N33's rule on per-seed scores: scores[task key][arm] -> array over the seeds. An arm missing on any task, or
    with a non-finite score, is not readable."""
    names = [*nt.DEEP, *nt.OTHER]

    def pair_mean(first: str, n: str) -> float:
        s = scores[tasks[n].key]
        return float(np.maximum(s[first], s[SECOND]).mean())

    base_abs = {n: pair_mean(BASE, n) for n in names}
    out = {"second": SECOND, "base_first": BASE, "base_pair_task_means": base_abs,
           "base_pair_worst9": min(base_abs.values()), "first_uploads": {}}
    passed = []
    for x in candidates:
        missing = [n for n in names if x not in scores[tasks[n].key]
                   or not np.isfinite(scores[tasks[n].key][x]).all()]
        if missing:
            out["first_uploads"][x] = {"readable": False, "missing_on": missing, "passes": False}
            continue
        absx = {n: pair_mean(x, n) for n in names}
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
                "b_S0.3_and_S0.6_not_lower": bool(s[0.3] >= 0 and s[0.6] >= 0),
                "c_worst9_not_lower": bool(worst >= out["base_pair_worst9"] - nt.TOL),
                "d_leave_one_task_out_positive": bool(min(loo.values()) > 0)}
        ok = all(rule.values())
        if ok:
            passed.append(x)
        own = {n: float(scores[tasks[n].key][x].mean() - scores[tasks[n].key][BASE].mean()) for n in names}
        out["first_uploads"][x] = {"readable": True, "tasks": per, "own_minus_fgw": own, "D": d, "O": o,
                                   **{f"S{p}": v for p, v in s.items()},
                                   "mean9": float(np.mean(list(per.values()))), "worst9_abs": worst,
                                   "S0.45_leave_one_task_out_min": min(loo.values()),
                                   "S0.45_leave_one_task_out": loo, "rule": rule, "passes": ok}
    out["requested_first_upload"] = max(passed, key=lambda x: out["first_uploads"][x]["S0.45"]) if passed else None
    return out


def reading(run_dirs) -> dict:
    scores, why, meta = nt._load([*nt.STORED, E28_RUNS, *run_dirs])
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (BASE, SECOND):
            if arm not in scores.get(t.key, {}):
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": why}
    return {"void": False, "why": [], "runs": meta, **price(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    lines = [f"Reading N37 (runs {', '.join(res['runs'])}): the new arms as first uploads beside the second upload "
             f"{res['second']}, against the pair as built ({res['base_first']}, {res['second']}); aggregation of "
             "committed scores, not a try",
             f"requested first upload: {res['requested_first_upload'] or 'none (the first upload stays fgw)'}",
             "  first upload           " + " ".join(f"{n:>6}" for n in names)
             + "  |      D      O  " + " ".join(f"S({p})" for p in nt.PS) + "   mean9  worst9  LOO-min  rule a b c d"]
    for x, v in res["first_uploads"].items():
        if not v["readable"]:
            lines.append(f"  {x:<22} not readable: no scores on {', '.join(v['missing_on'])}")
            continue
        r = v["rule"]
        lines.append(f"  {x:<22} " + " ".join(f"{v['tasks'][n]:+6.2f}" for n in names)
                     + f"  | {v['D']:+6.2f} {v['O']:+6.2f}  " + "  ".join(f"{v[f'S{p}']:+5.2f}" for p in nt.PS)
                     + f"   {v['mean9']:+5.2f}  {v['worst9_abs']:6.2f}  {v['S0.45_leave_one_task_out_min']:+7.2f}  "
                     + " ".join("y" if r[k] else "n" for k in sorted(r)))
        lines.append(f"  {'  its own mean - fgw':<22} " + " ".join(f"{v['own_minus_fgw'][n]:+6.2f}" for n in names))
    lines.append(f"  the pair as built: worst9 {res['base_pair_worst9']:.2f}; task means "
                 + " ".join(f"{res['base_pair_task_means'][n]:.2f}" for n in names))
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
    (out / "n37.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n37.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
