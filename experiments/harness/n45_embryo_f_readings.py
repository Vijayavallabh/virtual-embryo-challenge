#!/usr/bin/env python
"""Readings of tries E41a and E41b of experiments/harness/TRIES.md (registered before any run): two embryo arms the
record had left unbuilt, scored on the nine practice tasks by `n45_embryo_f_arms.py` on the branch `campaign/emb-f`:
  `fgw_potts`    the entry with a zero-temperature Potts sorting of which drawn cell sits at which position, stopped
                 at the same-group neighbour share interpolated between the two brackets;
  `fgw_partial`  the entry with a partial FGW coupling whose transported share is the overlap of the two brackets'
                 anchor label distributions.

What it reads. Only per-seed board scores committed by counted tries: the stored runs of E24, E25 and E28
(`nine_tasks.STORED` and `results/e28_nine_tasks/tasks`) and the committed task files of tries E41 given by `--runs`.
It opens no data file, runs no candidate and no scorer.

Validity, a check that can fail. The run of tries E41 rescored the entry `fgw` and the floor. Each must agree with
its committed scores seed for seed within nine_tasks.TOL on every task, or the reading is void. The largest
difference per arm is reported. The readings then use the committed scores of `fgw` and of the second upload, as
every earlier reading did, and the run's scores of the two new arms.

Three readings, each an existing rule applied unchanged (the same three as the rows of tries E36 and E38 onwards),
all fixed in the row before the run:
  adoption       the adoption rule (vec.adoption.keep_change) for each arm against the entry `fgw` on the two counted
                 embryo backtests (the nine tasks' e75 and e725), seed for seed;
  first upload   row N33's rule (n37_first_upload_new_arms.price), each arm as a first upload beside the second
                 upload as built;
  second upload  rows E34 and E35's rule (n38_second_beside_trend.price with the first upload `fgw`), each arm as a
                 second upload beside `fgw`, against `fgw_aligned_series`.
A pass of the adoption rule leads to a ruling request. A pass of a pair rule leads to a cold adversarial check and
then, if it survives, to a priced request. Anything else closes the line. Reported, deciding nothing: each arm's own
mean minus `fgw`'s on every task, and on every task the number of seeds in which the arm's score is `fgw`'s (within
nine_tasks.TOL).

Void on a missing task or seed, a dirty tree, or the rescored `fgw` or floor off its committed scores by more than
nine_tasks.TOL. An arm without scores on every task is reported as not readable and cannot pass.

  python -I -B experiments/harness/n45_embryo_f_readings.py --runs DIR --out DIR
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
BASE_OF = {"fgw_potts": ENTRY, "fgw_partial": ENTRY}  # each arm without its one idea
ARMS = tuple(BASE_OF)
ADOPTION_ARMS = ARMS  # both are the first-upload form: the entry with one step changed
RESCORED = ("floor", ENTRY)  # scored again by the run; must match the record
BACKTESTS = ("e75", "e725")  # the two counted embryo backtests among the nine tasks
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"


def adoption(scores: dict, tasks: dict, arm: str, base: str = ENTRY) -> dict:
    """The adoption rule on the two counted backtests, seed for seed (rows are in seed order)."""
    paired = {}
    for n in BACKTESTS:
        s = scores[tasks[n].key]
        paired[n] = {i: (float(b), float(a)) for i, (b, a) in enumerate(zip(s[base], s[arm], strict=True))}
    return keep_change(paired)


def rescored(committed: dict, run: dict, tasks: dict) -> tuple[dict, list[str]]:
    """Per rescored arm, the largest absolute difference between the run's per-seed scores and the committed ones over
    the nine tasks, and the reasons the run is void (an arm missing from the run, or off by more than nt.TOL)."""
    worst, why = {}, []
    for arm in RESCORED:
        big = 0.0
        for n, t in tasks.items():
            new, old = run.get(t.key, {}).get(arm), committed.get(t.key, {}).get(arm)
            if new is None:
                why.append(f"{n}: the run did not rescore {arm}")
                continue
            if old is None:
                why.append(f"{n}: no committed scores for {arm}")
                continue
            d = float(np.abs(np.asarray(new, float) - np.asarray(old, float)).max())
            big = max(big, d)
            if not d <= nt.TOL:
                why.append(f"{n}: {arm} differs from its committed scores by {d:.2e}")
        worst[arm] = big
    return worst, why


def check(scores: dict, tasks: dict) -> dict:
    names = [*nt.DEEP, *nt.OTHER]
    readable = {a: all(a in scores[tasks[n].key] and np.isfinite(scores[tasks[n].key][a]).all() for n in names)
                for a in ARMS}
    own = {a: {n: float(scores[tasks[n].key][a].mean() - scores[tasks[n].key][BASE_OF[a]].mean()) for n in names}
           for a in ARMS if readable[a]}
    same = {a: {n: int((np.abs(scores[tasks[n].key][a] - scores[tasks[n].key][BASE_OF[a]]) <= nt.TOL).sum())
                for n in names} for a in ARMS if readable[a]}
    out = {"arms": list(ARMS), "readable": readable, "own_minus_base": own, "seeds_equal_to_base": same,
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
    committed, why, meta = nt._load([*nt.STORED, E28_RUNS])
    run, why_run, meta_run = nt._load(list(run_dirs))
    why, meta = [*why, *why_run], sorted({*meta, *meta_run})
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in sorted({ENTRY, SECOND, *BASE_OF.values()}):
            if arm not in committed.get(t.key, {}):
                why.append(f"{n}: no committed scores for {arm}")
        if t.key not in run:
            why.append(f"{n}: no task file in the run")
    worst, why_rescored = rescored(committed, run, tasks)
    why = [*why, *why_rescored]
    if why:
        return {"void": True, "why": sorted(set(why)), "rescored_max_abs_difference": worst}
    scores = {t.key: {**committed[t.key], **{a: run[t.key][a] for a in ARMS if a in run[t.key]}}
              for t in tasks.values()}
    return {"void": False, "why": [], "runs": meta, "rescored_max_abs_difference": worst, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    passed = [x for x in (("the adoption rule is passed by " + ", ".join(res["adopted"])) if res["adopted"] else "",
                          f"first upload put forward: {res['requested_first_upload']}"
                          if res["requested_first_upload"] else "",
                          f"second upload put forward: {res['requested_second_upload']}"
                          if res["requested_second_upload"] else "") if x]
    lines = [f"Readings of tries E41a and E41b (runs {', '.join(res['runs'])}): the Potts sorting and the partial "
             "coupling on the nine practice tasks; aggregation of committed scores",
             "outcome: " + ("closed: no rule is passed" if res["closed"] else "; ".join(passed)),
             "  rescored arms against their committed scores, largest difference over tasks and seeds: "
             + ", ".join(f"{a} {res['rescored_max_abs_difference'][a]:.1e}" for a in RESCORED)
             + f" (void above {nt.TOL})"]
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
        lines.append(f"    seeds (of 10) in which {arm} scores as {BASE_OF[arm]}: "
                     + " ".join(f"{n} {res['seeds_equal_to_base'][arm][n]}" for n in names))
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
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of tries E41")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n45.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n45.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
