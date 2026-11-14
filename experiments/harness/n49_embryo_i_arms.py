#!/usr/bin/env python
"""The embryo batch I scored on the nine practice tasks, each given its trend samples (campaign branch campaign/emb-i;
the experiments/harness/TRIES.md row is registered by the main session before any run).

Why. The research report's recipe R2 ("reference-residual factor for an interior reference") reads a sample's
deviation from the straight line through its two neighbouring stages as that sample's own run effect and removes part
of it. On the embryo the record read the trend render (fgw_trend_render: every drawn cell of both brackets moved to the
trend's profile at the TARGET's day), never R2's own form, which corrects each bracket at its own day and so keeps the
stage difference between the brackets. This arm is R2's form with no fitted number (vec.t2_build.RESID_RENDER_RULE).

Arms:
  fgw_resid         fgw's draws and positions; each bracket with an earlier and a later sample moved by the share
                    L = 1 / (1 + w1^2 + w2^2) of its residual from the line through those two samples
and fgw and fgw_trend_render, rescored under their own names, so that the run can be checked against the committed
scores (fgw on the extended task is what it is on the two-bracket task: the size rule does not read a trend sample).

Tasks: the nine of experiments/harness/nine_tasks.py, each with its inputs extended by n36's sample rule
(n36_trend_second_upload.all_tasks, imported unchanged), 10 seeds (vec.adoption.SEEDS). The task file is the
harness's own result (vec.harness.run_task), so nine_tasks._load reads it; it also holds the sample rule's samples, as
n36's does. fgw_resid's build record (per seed, per bracket: corrected, the neighbours' stage names, w1, w2, L, the
share of panel genes whose factor was clipped and the mean |log factor|) is written under "build_records". An arm
whose build cannot be made on a task is recorded under "not_buildable" with the seed and the exception's message,
every number masked, and left out of the rows; the floor and the other arms are still scored (n36's guard and merge,
imported unchanged). This script reads nothing.

  run   --task NAME --out DIR   score the floor, fgw, fgw_trend_render and fgw_resid on ONE task; start the nine
                                together so that they record one commit; refuses an existing task file
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n36_trend_second_upload as n36  # noqa: E402  (the sample rule, the per-arm guard, the merge)
import nine_tasks as nt  # noqa: E402, F401  (the task file loads through nt._load)
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

NEW = ("fgw_resid",)
RESCORED = ("fgw", "fgw_trend_render")  # under their own names: the run's check against the committed scores
ARMS = (harness.FLOOR, *RESCORED, *NEW)


def _recorded(c: harness.Candidate, records: dict) -> harness.Candidate:
    """c with its build record kept per seed in records[c.name][seed] (the build's `build_record` keyword): the same
    name, params and applicability."""
    inner = c.build

    def build(ctx):
        rec: dict = {}
        out = inner(ctx, build_record=rec)
        records.setdefault(c.name, {})[str(int(ctx.seed - harness.BUILD_SEED_OFFSET))] = rec
        return out

    return dataclasses.replace(c, build=build)


def _arm(name: str, failures: dict, records: dict, key: str) -> harness.Candidate:
    if name in NEW:
        return n36._guarded(_recorded(harness.CANDIDATES[name], records), failures, key)
    return n36._guarded(harness.CANDIDATES[name], failures, key)


def score_task(t: harness.PseudoTask, data_dir: Path, runner=None, seeds=SEEDS) -> dict:
    """The floor, fgw, fgw_trend_render and fgw_resid on t. A build of fgw_trend_render or fgw_resid that raises stops
    no other arm: the task is then scored as the floor with fgw, and the floor with each of the two, and an arm that
    fails is recorded under "not_buildable" (the seed and the exception's message, every number masked) and left out of
    the rows. A failure elsewhere (data, floor, fgw, scoring) stops the task."""
    runner = harness.run_task if runner is None else runner
    floor, entry = harness.CANDIDATES[harness.FLOOR], harness.CANDIDATES["fgw"]
    others = (*RESCORED[1:], *NEW)
    failures: dict = {}
    records: dict = {}
    try:
        res = runner(t, [floor, entry, *[_arm(n, failures, records, t.key) for n in others]], seeds=seeds,
                     data_dir=data_dir)
        res["not_buildable"] = {}
        res["build_records"] = records
        return res
    except Exception:
        if not failures:
            raise
    records = {}
    parts, not_buildable = [runner(t, [floor, entry], seeds=seeds, data_dir=data_dir)], {}
    for n in others:
        fails: dict = {}
        try:
            parts.append(runner(t, [floor, _arm(n, fails, records, t.key)], seeds=seeds, data_dir=data_dir))
        except Exception:
            if n not in fails:
                raise
            not_buildable[n] = fails[n]
            records.pop(n, None)
    res = n36._merge(parts)
    res["not_buildable"] = not_buildable
    res["build_records"] = records
    return res


def cmd_run(a) -> None:
    tasks = n36.all_tasks()
    if a.task not in tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(tasks)}")
    t = tasks[a.task]
    cands = [harness.CANDIDATES[n] for n in ARMS]
    missing = [c.name for c in cands if c.name != harness.FLOOR and not c.applies(t)]
    if missing:
        sys.exit(f"{t.key}: arms that do not apply: {missing}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{sw.safe(t.key)}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = score_task(t, Path(a.data))
    res["task_source"] = "experiments/harness/n49_embryo_i_arms.py"
    res["sample_rule"] = {"series": list(n36.SERIES),
                          "samples": n36.trend_sample_files(t.inputs["before"], t.inputs["after"], t.target)}
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    {"run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    main()
