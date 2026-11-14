#!/usr/bin/env python
"""Readings of tries E38a to E38d of experiments/harness/TRIES.md (registered before any run): the embryo arms with
a mirror-only correction of FGW's partners (the ledger's pruned line of that name: keep FGW's moves; put each partner
on its anchor's side of the symmetric direction) in two forms, scored on the nine practice tasks by
`n42_mirrorfix_arms.py` on the branch `campaign/emb-mirror`:
  `fgw_mirrorfix`, `fgw_mirrorfix_series`    every partner on the wrong side of the direction is mirrored;
  `fgw_mirrorgate`, `fgw_mirrorgate_series`  the same, only in a stage whose partners are mirrored along it
                                             (smallest eigenvalue below 0).

What it reads. Only per-seed board scores committed by counted tries: the stored runs of E24, E25 and E28
(`nine_tasks.STORED` and `results/e28_nine_tasks/tasks`) and the committed task files of tries E38 given by `--runs`.
It opens no data file, runs no candidate and no scorer.

Three readings, each an existing rule applied unchanged (the same three as tries E36), all fixed in the row before the
run:
  adoption       the adoption rule (vec.adoption.keep_change) for each first-upload form (`fgw_mirrorfix`,
                 `fgw_mirrorgate`) against the entry `fgw` on the two counted embryo backtests (the nine tasks' e75
                 and e725), seed for seed;
  first upload   row N33's rule (n37_first_upload_new_arms.price), each arm as a first upload beside the second
                 upload as built;
  second upload  rows E34 and E35's rule (n38_second_beside_trend.price with the first upload `fgw`), each arm as a
                 second upload beside `fgw`, against `fgw_aligned_series`.
A pass of the adoption rule leads to a ruling request. A pass of a pair rule leads to a cold adversarial check and
then, if it survives, to a priced request. Anything else closes the line. Reported, deciding nothing: each arm's own
mean minus its base's on every task (the arm it differs from by the correction alone: `fgw` and `fgw_series`).

Void on a missing task or seed, a dirty tree, or a floor that differs from the committed runs by more than
nine_tasks.TOL (nine_tasks._load checks every arm two runs share). An arm without scores on every task is reported as
not readable and cannot pass.

  python -I -B experiments/harness/n42_mirrorfix_readings.py --runs DIR --out DIR
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n37_first_upload_new_arms as n37  # noqa: E402
import n38_second_beside_trend as n38  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec.adoption import keep_change  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ENTRY = "fgw"  # the first upload as built
SECOND = nt.RA  # fgw_aligned_series, the second upload as built
BASE_OF = {"fgw_mirrorfix": ENTRY, "fgw_mirrorfix_series": "fgw_series",  # each arm without the correction
           "fgw_mirrorgate": ENTRY, "fgw_mirrorgate_series": "fgw_series"}
ARMS = tuple(BASE_OF)
ADOPTION_ARMS = ("fgw_mirrorfix", "fgw_mirrorgate")  # the first-upload forms: the entry with its partners corrected
BACKTESTS = ("e75", "e725")  # the two counted embryo backtests among the nine tasks
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"


def adoption(scores: dict, tasks: dict, arm: str, base: str = ENTRY) -> dict:
    """The adoption rule on the two counted backtests, seed for seed (rows are in seed order)."""
    paired = {}
    for n in BACKTESTS:
        s = scores[tasks[n].key]
        paired[n] = {i: (float(b), float(a)) for i, (b, a) in enumerate(zip(s[base], s[arm], strict=True))}
    return keep_change(paired)


def check(scores: dict, tasks: dict) -> dict:
    names = [*nt.DEEP, *nt.OTHER]
    readable = {a: all(a in scores[tasks[n].key] and np.isfinite(scores[tasks[n].key][a]).all() for n in names)
                for a in ARMS}
    own = {a: {n: float(scores[tasks[n].key][a].mean() - scores[tasks[n].key][BASE_OF[a]].mean()) for n in names}
           for a in ARMS if readable[a]}
    out = {"arms": list(ARMS), "readable": readable, "own_minus_base": own,
           "adoption": {a: (adoption(scores, tasks, a) if readable[a] else None) for a in ADOPTION_ARMS},
           "first_upload": n37.price(scores, tasks, candidates=ARMS),
           "second_upload": n38.price(scores, tasks, first=ENTRY, candidates=ARMS)}
    out["adopted"] = [a for a in ADOPTION_ARMS if out["adoption"][a] and out["adoption"][a]["keep"]]
    out["requested_first_upload"] = out["first_upload"]["requested_first_upload"]
    advised = out["second_upload"]["advised_second"]
    out["requested_second_upload"] = None if advised == SECOND else advised
    out["closed"] = not (out["adopted"] or out["requested_first_upload"] or out["requested_second_upload"])
    return out


def reading(run_dirs) -> dict:
    scores, why, meta = nt._load([*nt.STORED, E28_RUNS, *run_dirs])
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in sorted({ENTRY, SECOND, *BASE_OF.values()}):
            if arm not in scores.get(t.key, {}):
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": sorted(set(why))}
    return {"void": False, "why": [], "runs": meta, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    passed = [x for x in (("the adoption rule is passed by " + ", ".join(res["adopted"])) if res["adopted"] else "",
                          f"first upload put forward: {res['requested_first_upload']}"
                          if res["requested_first_upload"] else "",
                          f"second upload put forward: {res['requested_second_upload']}"
                          if res["requested_second_upload"] else "") if x]
    lines = [f"Readings of tries E38a to E38d (runs {', '.join(res['runs'])}): the mirror-only correction on the nine "
             "practice tasks; aggregation of committed scores",
             "outcome: " + ("closed: no rule is passed" if res["closed"] else "; ".join(passed))]
    for arm, a in res["adoption"].items():
        if a is None:
            lines.append(f"  adoption rule, {arm} against {ENTRY}: not readable")
            continue
        parts = [f"{n} {t['mean_gain']:+.2f} ({t['seeds_ahead']}/10)" for n, t in a["tasks"].items()]
        lines.append(f"  adoption rule, {arm} against {ENTRY}: " + ", ".join(parts)
                     + f" -> {'kept' if a['keep'] else 'not kept'} ({a['reason']})")
    for arm in res["arms"]:
        if not res["readable"][arm]:
            lines.append(f"  {arm}: not readable")
            continue
        lines.append(f"  {arm} minus {BASE_OF[arm]}, own means: "
                     + " ".join(f"{n} {res['own_minus_base'][arm][n]:+.2f}" for n in names))
    for title, block, key in (("as a first upload (row N33's rule)", res["first_upload"], "first_uploads"),
                              ("as a second upload beside fgw (rows E34 and E35's rule)", res["second_upload"],
                               "second_uploads")):
        lines.append(f"  {title}:" + " " * 6 + "D      O  " + " ".join(f"S({p})" for p in nt.PS)
                     + "  LOO-min  rule a b c d")
        for arm in res["arms"]:
            v = block[key][arm]
            if not v["readable"]:
                lines.append(f"    {arm:<32} not readable: no scores on {', '.join(v['missing_on'])}")
                continue
            r = v["rule"]
            lines.append(f"    {arm:<32} {v['D']:+6.2f} {v['O']:+6.2f}  "
                         + "  ".join(f"{v[f'S{p}']:+5.2f}" for p in nt.PS)
                         + f"  {v['S0.45_leave_one_task_out_min']:+7.2f}  "
                         + " ".join("y" if r[k] else "n" for k in sorted(r)))
            lines.append("      per task " + " ".join(f"{n} {v['tasks'][n]:+.2f}" for n in names))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of tries E38")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n42.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n42.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
