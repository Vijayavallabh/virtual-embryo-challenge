#!/usr/bin/env python
"""Tries E35a and E35b of experiments/harness/TRIES.md (campaign branch campaign/emb-trend2; the row is registered
before any run): the embryo trend rendering read as a second upload on the nine practice tasks, each task given its
trend samples.

Why. fgw_trend_render (fgw's draws and positions; every drawn cell of both brackets moved by a per-gene factor to the
profile a least-squares trend through all the task's samples gives at the target's day, vec.t2_build.TREND_RENDER_RULE)
came out one-sided on the two counted embryo backtests (try E33), which have three input stages each. A one-sided arm
is a second-upload question, read on the nine practice tasks by n35's pair rule. Seven of the nine carry only their two
brackets, where the trend does not apply; this script gives every task its trend samples.

The sample rule. A task's trend samples are its two brackets plus every stage of the embryo board's own released series
(SERIES: E6.75, E7.25, E7.5, E8.0) that lies strictly outside the bracket interval and is not the task's target. A
stage strictly inside the bracket interval is never a sample (the hidden test has none). The two counted backtests get
exactly the inputs configs/pseudo_tasks.yaml gives them; the test pair (E7.5 + E8.0 -> E7.75) would get its brackets
plus E6.75 and E7.25. A stage added by the rule takes the role vec.harness.TREND_SAMPLE_PREFIX + its stem
("trend_E6.75"), so "before", "after", "earlier" and "later" keep their meaning, and the size rule
(vec.harness.target_rms) does not read it: fgw's positions and every bracket-only path are what they are on the
two-bracket task.

Arms, each a second upload beside the first upload fgw:
  fgw_trend_render   unchanged (E35a)
  fgw_aligned_trend  try E10's aligned moves with fgw_trend_render's expression (E35b)
Only the two arms are scored here, with the floor, which every harness run scores and which must agree with the
committed runs to 0.01 per seed. fgw and fgw_aligned_series are NOT rescored: their per-seed scores come from the
committed runs of E24, E25 and E28 through nine_tasks._load, as in n35.

Rule (n35's, NEW = the two arms), fixed in the row before the run. A second upload X is put forward only if the pair
(fgw, X) beats the pair (fgw, fgw_aligned_series) with all of:
  (a) S(0.45) > BAR (0.3);
  (b) S(0.3) >= 0 and S(0.6) >= 0;
  (c) the pair (fgw, X)'s worst9 is not lower than the pair (fgw, fgw_aligned_series)'s by more than TOL (0.01);
  (d) S(0.45) stays above 0 whichever single task is left out.
If both pass, the one with the larger S(0.45); otherwise nothing is requested. Nothing changes an upload without the
captain's ruling (docs/S2_RULES.md).

The reading is void on any missing task, arm or seed, a dirty tree, a run file whose inputs are not the sample rule's,
a new run without the floor, a floor that differs from the committed runs by more than 0.01, rescored fgw or
fgw_aligned_series, or more than one commit among the new run's task files. An arm whose build cannot be made on a
task (for example the trend's refusal "no extraembryonic label held by every sample") is recorded per task as not
buildable with the exception's message; the reading reports it as "not readable", never as passing, and still reads
the other arm.

Also written, deciding nothing: per task the two arms' own means and their paired difference from fgw, and
fgw_trend_render's per-seed scores on e75 and e725 (the main session compares them with try E33's committed run).

  run   --task NAME --out DIR   score the floor and the two arms on ONE task (half, late, new1, new2, e75, e725, x1,
                                x2, x3); start the nine together so that they record one commit
  read  --runs DIR --out DIR    apply the rule; write n36.json and n36.txt
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nine_tasks as nt  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

BASE_FIRST = "fgw"  # the first upload as built
BASE_SECOND = nt.RA  # fgw_aligned_series, the second upload since row R4
NEW = ("fgw_trend_render", "fgw_aligned_trend")
ARMS = (harness.FLOOR, *NEW)
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
LEAVE_OUT_P = 0.45
SERIES = ("T2/E6.75.h5ad", "T2/E7.25.h5ad", "T2/E7.5.h5ad", "T2/E8.0.h5ad")  # the embryo board's own released series
E33_CHECK = ("e75", "e725")  # fgw_trend_render's per-seed scores here are compared with try E33's by the main session
NOT_BUILDABLE = "not buildable on this task"
NOT_READABLE = "not readable"


# ------------------------------------------------------------------------------------------------ the sample rule
def stage_day(f: str, days: dict | None = None) -> float:
    """A stage's day from vec.harness.STAGE_DAYS (or `days`, for a stage that is not released); refuses a stage
    without one."""
    d = (days or {}).get(f, harness.STAGE_DAYS.get(f))
    if d is None:
        raise ValueError(f"{f}: no stage day (vec.harness.STAGE_DAYS)")
    return float(d)


def trend_sample_files(before: str, after: str, target: str, days: dict | None = None) -> list[str]:
    """The sample rule: the two brackets plus every SERIES stage strictly outside (day(before), day(after)) that is not
    the target, in day order. `days` gives a day to a stage STAGE_DAYS lacks (the test target, E7.75)."""
    lo, hi, tday = stage_day(before, days), stage_day(after, days), stage_day(target, days)
    if not lo < tday < hi:
        raise ValueError(f"{target} is not inside its bracket {before} - {after}")
    out = [before, after]
    for f in SERIES:
        d = stage_day(f, days)
        if f == target or f in out or lo < d < hi:
            continue
        out.append(f)
    return sorted(out, key=lambda f: (stage_day(f, days), f))


def trend_role(f: str) -> str:
    return harness.TREND_SAMPLE_PREFIX + harness.stem(f)


def extend(t: harness.PseudoTask) -> harness.PseudoTask:
    """t with its inputs extended by the sample rule; every other field is t's. A configured input besides the
    brackets must be one of the samples (it keeps its role); each added stage takes the role trend_<stem>."""
    files = trend_sample_files(t.inputs["before"], t.inputs["after"], t.target)
    stray = sorted(f for r, f in t.inputs.items() if f not in files)
    if stray:
        raise ValueError(f"{t.key}: configured inputs outside the sample rule: {stray}")
    clash = sorted(r for r in t.inputs if r.startswith(harness.TREND_SAMPLE_PREFIX))
    if clash:
        raise ValueError(f"{t.key}: roles {clash} already use the trend-sample prefix")
    extra = {trend_role(f): f for f in files if f not in t.inputs.values()}
    return harness.check_task(dataclasses.replace(t, inputs={**t.inputs, **extra}))


def all_tasks() -> dict:
    """The nine tasks of nine_tasks.all_tasks(), each with its inputs extended by the sample rule."""
    return {n: extend(t) for n, t in nt.all_tasks().items()}


# ------------------------------------------------------------------------------------------------ run
def _message(e: BaseException, key: str) -> str:
    """The exception's type and message with the task key named as such and every number masked: no data value."""
    msg = str(e).replace(key, "<task>")
    return f"{type(e).__name__}: " + re.sub(r"[-+]?\d+(\.\d+)?([eE][-+]?\d+)?", "#", msg)


def _guarded(c: harness.Candidate, failures: dict, key: str) -> harness.Candidate:
    """c with its build recording (arm, seed, message) in `failures` before the exception goes on: the same name,
    params and applicability."""
    inner = c.build

    def build(ctx):
        try:
            return inner(ctx)
        except Exception as e:
            failures.setdefault(c.name, {"status": NOT_BUILDABLE, "seed": int(ctx.seed - harness.BUILD_SEED_OFFSET),
                                         "message": _message(e, key)})
            raise

    return dataclasses.replace(c, build=build)


def _merge(parts: list[dict]) -> dict:
    """One task result from runs of the same task that each scored the floor and some arms."""
    res = parts[0]
    for p in parts[1:]:
        for row, other in zip(res["rows"], p["rows"], strict=True):
            if row["seed"] != other["seed"] or row["scores"][harness.FLOOR] != other["scores"][harness.FLOOR]:
                raise ValueError(f"{res['task']['key']}: the floor differs between the per-arm runs")
            for k in ("scores", "skills", "metrics", "n_cells", "groups", "reads"):
                row[k] = {**row[k], **other[k]}
        res["candidates"] = [*res["candidates"], *[c for c in p["candidates"] if c not in res["candidates"]]]
        res["params"] = {**res["params"], **p["params"]}
    return res


def score_task(t: harness.PseudoTask, data_dir: Path, runner=harness.run_task, seeds=SEEDS) -> dict:
    """The floor and both arms on t. A build that raises stops neither the other arm nor the floor: the task is then
    scored arm by arm, and an arm that fails is recorded under "not_buildable" (the seed and the exception's message)
    and left out of the rows. A failure outside an arm's build (data, floor, scoring) stops the task."""
    floor = harness.CANDIDATES[harness.FLOOR]
    failures: dict = {}
    try:
        res = runner(t, [floor, *[_guarded(harness.CANDIDATES[n], failures, t.key) for n in NEW]], seeds=seeds,
                     data_dir=data_dir)
        res["not_buildable"] = {}
        return res
    except Exception:
        if not failures:
            raise
    parts, not_buildable = [], {}
    for n in NEW:
        fails: dict = {}
        try:
            parts.append(runner(t, [floor, _guarded(harness.CANDIDATES[n], fails, t.key)], seeds=seeds,
                                data_dir=data_dir))
        except Exception:
            if n not in fails:
                raise
            not_buildable[n] = fails[n]
    if not parts:
        parts.append(runner(t, [floor], seeds=seeds, data_dir=data_dir))
    res = _merge(parts)
    res["not_buildable"] = not_buildable
    return res


def cmd_run(a) -> None:
    tasks = all_tasks()
    if a.task not in tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(tasks)}")
    t = tasks[a.task]
    missing = [n for n in NEW if not harness.CANDIDATES[n].applies(t)]
    if missing:
        sys.exit(f"{t.key}: arms that do not apply: {missing}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{sw.safe(t.key)}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = score_task(t, Path(a.data))
    res["task_source"] = "experiments/harness/n36_trend_second_upload.py"
    res["sample_rule"] = {"series": list(SERIES), "samples": trend_sample_files(t.inputs["before"], t.inputs["after"],
                                                                                 t.target)}
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def _arms_in(d: Path) -> set[str]:
    """Every arm a run directory's task files hold."""
    arms: set[str] = set()
    for p in sorted(Path(d).glob("T2e_*.json")):
        for row in json.loads(p.read_text())["rows"]:
            arms |= set(row["scores"])
    return arms


def evaluate(scores: dict, tasks: dict, new) -> dict:
    """n35's rule on per-seed scores (scores[task key][arm] -> array over SEEDS), with the second uploads `new`."""
    names = [*nt.DEEP, *nt.OTHER]

    def pair_mean(second: str, n: str) -> float:
        s = scores[tasks[n].key]
        return float(np.maximum(s[BASE_FIRST], s[second]).mean())

    base_abs = {n: pair_mean(BASE_SECOND, n) for n in names}
    out = {"first": BASE_FIRST, "base_second": BASE_SECOND, "base_pair_task_means": base_abs,
           "base_pair_worst9": min(base_abs.values()), "second_uploads": {}}
    passed = []
    for x in new:
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
    out["put_forward"] = max(passed, key=lambda x: out["second_uploads"][x]["S0.45"]) if passed else None
    return out


def describe(scores: dict, tasks: dict) -> dict:
    """Deciding nothing: per task each arm's own mean, fgw's, and the arm's paired difference from fgw; and
    fgw_trend_render's per-seed scores on E33_CHECK."""
    per = {}
    for n, t in tasks.items():
        s = scores.get(t.key, {})
        row = {}
        for x in NEW:
            if x not in s:
                row[x] = None
                continue
            d = s[x] - s[BASE_FIRST]
            row[x] = {"mean": float(s[x].mean()), "fgw_mean": float(s[BASE_FIRST].mean()),
                      "paired_diff_from_fgw": float(d.mean()), "ahead_of_fgw": int((d > 0).sum())}
        per[n] = row
    e33 = {n: ([float(v) for v in scores[tasks[n].key]["fgw_trend_render"]]
               if "fgw_trend_render" in scores.get(tasks[n].key, {}) else None) for n in E33_CHECK}
    return {"per_task": per, "fgw_trend_render_per_seed_for_the_E33_check": e33}


def reading(runs_dir: Path) -> dict:
    runs_dir = Path(runs_dir)
    why = []
    stored = [*nt.STORED, E28_RUNS]
    for d in stored:
        clash = _arms_in(d) & set(NEW)
        if clash:
            why.append(f"{d}: already holds {sorted(clash)}")
    rescored = _arms_in(runs_dir) & {BASE_FIRST, BASE_SECOND}
    if rescored:
        why.append(f"{runs_dir}: rescored {sorted(rescored)} (they must come from the committed runs only)")
    tasks = all_tasks()
    files = {}
    for p in sorted(runs_dir.glob("T2e_*.json")):
        r = json.loads(p.read_text())
        files[r["task"]["key"]] = r
    commits = sorted({r["provenance"]["commit"] for r in files.values()})
    if len(commits) > 1:
        why.append(f"{runs_dir}: the task files record more than one commit ({', '.join(c[:7] for c in commits)})")
    scores, why_load, meta = nt._load([*stored, runs_dir])
    why += why_load
    not_buildable: dict[str, dict] = {x: {} for x in NEW}
    for n, t in tasks.items():
        r = files.get(t.key)
        if r is None:
            why.append(f"{n}: no task file in {runs_dir}")
            continue
        if dict(r["task"]["inputs"]) != dict(t.inputs):
            why.append(f"{n}: the run's inputs are not the sample rule's")
        if any(harness.FLOOR not in row["scores"] for row in r["rows"]):
            why.append(f"{n}: the new run has no floor")
        for x, rec in (r.get("not_buildable") or {}).items():
            if x in not_buildable:
                not_buildable[x][n] = rec
        have = scores.get(t.key, {})
        for arm in (BASE_FIRST, BASE_SECOND):
            if arm not in have:
                why.append(f"{n}: no scores for {arm}")
        for x in NEW:
            if x not in have and n not in not_buildable[x]:
                why.append(f"{n}: no scores for {x}")
            if x in have and n in not_buildable[x]:
                why.append(f"{n}: {x} is both scored and recorded as not buildable")
    if why:
        return {"void": True, "why": why}
    readable = [x for x in NEW if not not_buildable[x]]
    res = {"void": False, "why": [], "runs": meta, **evaluate(scores, tasks, readable)}
    for x in NEW:
        if x not in readable:
            res["second_uploads"][x] = {"readable": False, "status": NOT_READABLE, "passes": False,
                                        "not_buildable_on": not_buildable[x]}
    res["second_uploads"] = {x: res["second_uploads"][x] for x in NEW}
    res["described"] = describe(scores, tasks)
    return res


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    lines = [f"n36, tries E35a and E35b (runs {', '.join(res['runs'])}): second uploads beside {res['first']}, against "
             f"the pair as built ({res['first']}, {res['base_second']})",
             f"put forward: {res['put_forward'] or 'none (the second upload stays ' + res['base_second'] + ')'}",
             "  second upload          " + " ".join(f"{n:>6}" for n in names)
             + "  |      D      O  " + " ".join(f"S({p})" for p in nt.PS) + "   mean9  worst9  LOO-min  rule a b c d"]
    for x, v in res["second_uploads"].items():
        if not v["readable"]:
            where = "; ".join(f"{n} (seed {r['seed']}: {r['message']})" for n, r in v["not_buildable_on"].items())
            lines.append(f"  {x:<22} {NOT_READABLE}: {NOT_BUILDABLE} on {where}")
            continue
        r = v["rule"]
        lines.append(f"  {x:<22} " + " ".join(f"{v['tasks'][n]:+6.2f}" for n in names)
                     + f"  | {v['D']:+6.2f} {v['O']:+6.2f}  " + "  ".join(f"{v[f'S{p}']:+5.2f}" for p in nt.PS)
                     + f"   {v['mean9']:+5.2f}  {v['worst9_abs']:6.2f}  {v['S0.45_leave_one_task_out_min']:+7.2f}  "
                     + " ".join("y" if r[k] else "n" for k in sorted(r)))
    lines.append(f"  the pair as built: worst9 {res['base_pair_worst9']:.2f}; task means "
                 + " ".join(f"{res['base_pair_task_means'][n]:.2f}" for n in names))
    lines.append("deciding nothing: each arm's own mean and its paired difference from fgw (seeds ahead)")
    for x in NEW:
        cells = []
        for n in names:
            v = res["described"]["per_task"][n][x]
            cells.append(f"{n} -" if v is None else
                         f"{n} {v['mean']:.2f} ({v['paired_diff_from_fgw']:+.2f}, {v['ahead_of_fgw']}/10)")
        lines.append(f"  {x}: " + "; ".join(cells))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} is not empty: refusing to overwrite")
    res = reading(Path(a.runs))
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n36.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n36.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
