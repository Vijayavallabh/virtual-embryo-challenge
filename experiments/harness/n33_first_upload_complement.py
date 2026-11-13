#!/usr/bin/env python
"""Reading N33 of experiments/harness/TRIES.md (registered before this script was run): the embryo's first upload
priced as the complement of its second upload, on the nine practice tasks of try E28.

Why. The second upload changed on 2026-11-13 (row R4: `fgw_aligned_series`, which pays where the hidden stage shares
its deeper bracket's profile). The first upload, `fgw`, was chosen before that change, and tries E25 and E28 priced
second uploads beside `fgw` or `fgw_aligned` only. No table prices another first upload beside the new second one.

What it reads. Only per-seed board scores already committed by tries E24, E25 and E28 (every arm on every one of the
nine tasks was scored and counted there): `results/e24_pricing/tasks`, `results/e25_pairs/tasks` and
`results/e28_nine_tasks/tasks`, through `nine_tasks._load`, which checks seeds, clean trees and agreement between runs.
It opens no data file, runs no candidate and no scorer, and scores nothing new: an aggregation, not a try.

Value of a first upload X: on each task, the mean over seeds of max(X, second) minus the same for max(fgw, second),
with second = fgw_aligned_series; D, O and S(p) as in nine_tasks.py (DEEP, OTHER, PS); mean9; worst9 = the lowest
task mean of the pair's better of two, in board points; and S(0.45) with each task left out in turn.

Rule, fixed in the row before the reading. A first upload other than fgw is put to the captain only if all hold:
  (a) S(0.45) > BAR (0.3);
  (b) S(0.3) >= 0 and S(0.6) >= 0;
  (c) its pair's worst9 is not lower than fgw's pair's worst9 by more than TOL;
  (d) S(0.45) stays above 0 whichever single task is left out.
Otherwise nothing is requested and the first upload stays fgw. A first upload that passes still fails the adoption rule
on the two counted backtests unless it is fgw_aligned's kind of change, so it could enter by the captain's ruling only,
and docs/S2_RULES.md then sends the second upload back to the captain.

  python -I -B experiments/harness/n33_first_upload_complement.py --out experiments/harness/results/n33_first_upload
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
CANDIDATES = ("fgw_aligned", "fgw_shape0", "fgw_render_hi", "fgw_shape0_hi", "fgw_series", "fgw_shape0_series")
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"


def reading() -> dict:
    scores, why, meta = nt._load([*nt.STORED, E28_RUNS])
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (BASE, SECOND, *CANDIDATES):
            if arm not in scores.get(t.key, {}):
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": why}
    names = [*nt.DEEP, *nt.OTHER]

    def pair_mean(first: str, n: str) -> float:
        s = scores[tasks[n].key]
        return float(np.maximum(s[first], s[SECOND]).mean())

    base_abs = {n: pair_mean(BASE, n) for n in names}
    out = {"void": False, "why": [], "runs": meta, "second": SECOND, "base_first": BASE,
           "base_pair_task_means": base_abs, "base_pair_worst9": min(base_abs.values()), "first_uploads": {}}
    passed = []
    for x in CANDIDATES:
        absx = {n: pair_mean(x, n) for n in names}
        per = {n: absx[n] - base_abs[n] for n in names}
        d = float(np.mean([per[n] for n in nt.DEEP]))
        o = float(np.mean([per[n] for n in nt.OTHER]))
        s = {p: p * d + (1 - p) * o for p in nt.PS}
        loo = {}
        for drop in names:
            dd = [per[n] for n in nt.DEEP if n != drop]
            oo = [per[n] for n in nt.OTHER if n != drop]
            loo[drop] = 0.45 * float(np.mean(dd)) + 0.55 * float(np.mean(oo))
        worst = min(absx.values())
        rule = {"a_S0.45_above_bar": bool(s[0.45] > nt.BAR),
                "b_S0.3_and_S0.6_not_lower": bool(s[0.3] >= 0 and s[0.6] >= 0),
                "c_worst9_not_lower": bool(worst >= out["base_pair_worst9"] - nt.TOL),
                "d_leave_one_task_out_positive": bool(min(loo.values()) > 0)}
        ok = all(rule.values())
        if ok:
            passed.append(x)
        out["first_uploads"][x] = {"tasks": per, "D": d, "O": o, **{f"S{p}": v for p, v in s.items()},
                                   "mean9": float(np.mean(list(per.values()))), "worst9_abs": worst,
                                   "S0.45_leave_one_task_out_min": min(loo.values()),
                                   "S0.45_leave_one_task_out": loo, "rule": rule, "passes": ok}
    out["requested_first_upload"] = max(passed, key=lambda x: out["first_uploads"][x]["S0.45"]) if passed else None
    return out


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    lines = [f"Reading N33 (runs {', '.join(res['runs'])}): first uploads beside the second upload {res['second']}, "
             f"against the pair as built ({res['base_first']}, {res['second']}); aggregation of committed scores, not "
             "a try",
             f"requested first upload: {res['requested_first_upload'] or 'none (the first upload stays fgw)'}",
             "  first upload         " + " ".join(f"{n:>6}" for n in names)
             + "  |      D      O  " + " ".join(f"S({p})" for p in nt.PS) + "   mean9  worst9  LOO-min  rule a b c d"]
    for x, v in res["first_uploads"].items():
        r = v["rule"]
        lines.append(f"  {x:<20} " + " ".join(f"{v['tasks'][n]:+6.2f}" for n in names)
                     + f"  | {v['D']:+6.2f} {v['O']:+6.2f}  " + "  ".join(f"{v[f'S{p}']:+5.2f}" for p in nt.PS)
                     + f"   {v['mean9']:+5.2f}  {v['worst9_abs']:6.2f}  {v['S0.45_leave_one_task_out_min']:+7.2f}  "
                     + " ".join("y" if r[k] else "n" for k in sorted(r)))
    lines.append(f"  the pair as built: worst9 {res['base_pair_worst9']:.2f}; task means "
                 + " ".join(f"{res['base_pair_task_means'][n]:.2f}" for n in names))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading()
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n33.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n33.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
