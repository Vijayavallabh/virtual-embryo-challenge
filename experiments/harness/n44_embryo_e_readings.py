#!/usr/bin/env python
"""Readings of tries E40a to E40d of experiments/harness/TRIES.md (registered before any run): four embryo arms
scored on the nine practice tasks by `n44_embryo_e_arms.py` on the branch `campaign/emb-e`:
  `fgw_aligned_mix`        one file at the second upload's cells and positions, a random half of its rows carrying
                           the second upload's series expression and the other half the aligned trend's (try H16's
                           within-file mixture, on this board);
  `fgw_aligned_series_mn`  the second upload as built with its render's counts redrawn per cell as one multinomial
                           (try H19's kernel, on this board);
  `fgw_quad_render`        the trend render of try E33 with a quadratic in stage day in place of its line;
  `fgw_aligned_quad`       the aligned trend of try E35b with the same quadratic.

What it reads. Only per-seed board scores committed by counted tries: the stored runs of E24, E25 and E28
(`nine_tasks.STORED` and `results/e28_nine_tasks/tasks`), the committed task files of try E35
(`results/n36_trend_second_upload/tasks`: `fgw_trend_render` and `fgw_aligned_trend`) and the committed task files of
tries E40 given by `--runs`. It opens no data file, runs no candidate and no scorer.

Validity, a check that can fail. The run of tries E40 rescored four arms the record already holds (`fgw`,
`fgw_aligned_series`, `fgw_trend_render`, `fgw_aligned_trend`) and the floor. Each must agree with its committed
scores seed for seed within nine_tasks.TOL on every task, or the reading is void. The largest difference per arm is
reported. The readings then use the committed scores of those arms, as every earlier reading did, and the run's
scores of the four new arms.

Three readings, each an existing rule applied unchanged (the same three as tries E36, E38 and E39), all fixed in the
row before the run:
  adoption       the adoption rule (vec.adoption.keep_change) for `fgw_quad_render` against the entry `fgw` on the
                 two counted embryo backtests (the nine tasks' e75 and e725), seed for seed;
  first upload   row N33's rule (n37_first_upload_new_arms.price), each arm as a first upload beside the second
                 upload as built;
  second upload  rows E34 and E35's rule (n38_second_beside_trend.price with the first upload `fgw`), each arm as a
                 second upload beside `fgw`, against `fgw_aligned_series`.
A pass of the adoption rule leads to a ruling request. A pass of a pair rule leads to a cold adversarial check and
then, if it survives, to a priced request. Anything else closes the line. Reported, deciding nothing: each arm's own
mean minus its base's on every task (the arm it differs from by its one idea: the second upload as built for the
mixture and the redraw, the line's arm for each quadratic), its own mean minus one more arm's (the mixture against
its other half `fgw_aligned_trend`, `fgw_quad_render` against `fgw`, `fgw_aligned_quad` against the second upload),
and for the mixture the number of tasks on which its own mean lies between its two halves' means.

Void on a missing task or seed, a dirty tree, or a rescored arm or floor off its committed scores by more than
nine_tasks.TOL. An arm without scores on every task is reported as not readable and cannot pass.

  python -I -B experiments/harness/n44_embryo_e_readings.py --runs DIR --out DIR
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
LINE_FIRST, LINE_SECOND = "fgw_trend_render", "fgw_aligned_trend"  # the line's arms (tries E33 and E35b)
MIX = "fgw_aligned_mix"
BASE_OF = {MIX: SECOND, "fgw_aligned_series_mn": SECOND,  # each arm without its one idea
           "fgw_quad_render": LINE_FIRST, "fgw_aligned_quad": LINE_SECOND}
ALSO_AGAINST = {MIX: LINE_SECOND, "fgw_quad_render": ENTRY, "fgw_aligned_quad": SECOND}  # reported only
ARMS = tuple(BASE_OF)
ADOPTION_ARMS = ("fgw_quad_render",)  # the first-upload form: the entry's cells with the quadratic's render
RESCORED = ("floor", ENTRY, SECOND, LINE_FIRST, LINE_SECOND)  # scored again by the run; must match the record
BACKTESTS = ("e75", "e725")  # the two counted embryo backtests among the nine tasks
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
N36_RUNS = HERE / "results" / "n36_trend_second_upload" / "tasks"


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

    def minus(a: str, b: str) -> dict:
        return {n: float(scores[tasks[n].key][a].mean() - scores[tasks[n].key][b].mean()) for n in names}

    own = {a: minus(a, BASE_OF[a]) for a in ARMS if readable[a]}
    also = {a: minus(a, b) for a, b in ALSO_AGAINST.items() if readable[a]}
    between = None
    if readable[MIX]:
        between = 0
        for n in names:
            s = scores[tasks[n].key]
            lo, hi = sorted((float(s[SECOND].mean()), float(s[LINE_SECOND].mean())))
            between += int(lo <= float(s[MIX].mean()) <= hi)
    out = {"arms": list(ARMS), "readable": readable, "own_minus_base": own, "own_minus_other": also,
           "mix_tasks_between_its_halves": between,
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
    committed, why, meta = nt._load([*nt.STORED, E28_RUNS, N36_RUNS])
    run, why_run, meta_run = nt._load(list(run_dirs))
    why, meta = [*why, *why_run], sorted({*meta, *meta_run})
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in sorted({ENTRY, SECOND, *BASE_OF.values(), *ALSO_AGAINST.values()}):
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
    lines = [f"Readings of tries E40a to E40d (runs {', '.join(res['runs'])}): the within-file mixture, the "
             "multinomial render and the quadratic trend on the nine practice tasks; aggregation of committed scores",
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
        if arm in res["own_minus_other"]:
            lines.append(f"    minus {ALSO_AGAINST[arm]}: "
                         + " ".join(f"{n} {res['own_minus_other'][arm][n]:+.2f}" for n in names))
    if res["mix_tasks_between_its_halves"] is not None:
        lines.append(f"  {MIX}: its own mean lies between its two halves' means ({SECOND}, {LINE_SECOND}) on "
                     f"{res['mix_tasks_between_its_halves']} of {len(names)} tasks")
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
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of tries E40")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n44.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n44.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
