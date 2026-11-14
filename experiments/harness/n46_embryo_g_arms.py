#!/usr/bin/env python
"""The embryo batch G scored on the nine practice tasks (campaign branch campaign/emb-g; the
experiments/harness/TRIES.md row is registered by the main session before any run).

Why. fgw_partial (tries E41b, branch campaign/emb-f) passed the adoption rule and the first-upload pair rule on the nine
practice tasks, but its row had warned that the arm is in large part a cloud that moves less: an anchor with no mass
does not move. The captain asked for "the same test with cells simply moved less". Three controls answer that, each the
entry fgw with ONE change to the displacement each anchor carries into the k-anchor mean, on the entry's own balanced
coupling, partners and rotation (vec.t2_build.interpolate_clouds, anchor_moves), s the share fgw_partial transports
(the overlap of the brackets' anchors' celltype distributions):

Arms:
  fgw_damp       every anchor's displacement multiplied by s
  fgw_part_rand  a uniformly random subset of round(s x its anchors) anchors of each bracket keeps its displacement,
                 the others get 0 (stream [seed, vec.harness.ANCHOR_STILL_STREAM])
  fgw_still      every anchor's displacement 0
and fgw and fgw_partial, rescored under their own names, so that the run can be checked against the committed scores.

Tasks: the nine of experiments/harness/nine_tasks.py (half, late, new1, new2, e75, e725, x1, x2, x3), each with exactly
the inputs nine_tasks.all_tasks() gives it, 10 seeds (vec.adoption.SEEDS). The task file is the harness's own result
(vec.harness.run_task, which records the commit and the tree state), so nine_tasks._load reads it. Each new arm's build
record (per seed: s where the arm uses it, the share of anchors with a non-zero displacement per bracket and the mean
norm of the anchors' displacements as a ratio to the entry's) is written under "build_records". This script reads
nothing: the reading is a separate script on main, registered in the ledger row by the main session.

An arm whose build cannot be made on a task (fgw_partial's solver failing, or the builder's own cloud guard refusing
the cloud) is recorded in the task file under "not_buildable" with the seed and the exception's message, every number
masked, and left out of the rows; fgw, the other arms and the floor are still scored (the per-arm guard, the merge and
the record are n36_trend_second_upload's own, imported unchanged), as in n45_embryo_f_arms.py.

  run   --task NAME --out DIR   score the floor, fgw, fgw_partial and the three arms on ONE task; start the nine
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
import n36_trend_second_upload as n36  # noqa: E402  (the per-arm guard, the merge and the not-buildable record)
import nine_tasks as nt  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

NEW = ("fgw_damp", "fgw_part_rand", "fgw_still")
ENTRY = "fgw"  # rescored under its own name: the run's check against the committed scores
PARTIAL = "fgw_partial"  # rescored under its own name; guarded, since its solver can refuse
GUARDED = (PARTIAL, *NEW)
ARMS = (harness.FLOOR, ENTRY, *GUARDED)


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
    c = harness.CANDIDATES[name]
    return n36._guarded(_recorded(c, records) if name in NEW else c, failures, key)


def score_task(t: harness.PseudoTask, data_dir: Path, runner=None, seeds=SEEDS) -> dict:
    """The floor, fgw, fgw_partial and the three arms on t. A build of fgw_partial or a new arm that raises stops
    neither the other arms, fgw nor the floor: the task is then scored as the floor with fgw, and the floor with each
    guarded arm, and an arm that fails is recorded under "not_buildable" (the seed and the exception's message, every
    number masked) and left out of the rows. A failure elsewhere (data, floor, fgw, scoring) stops the task."""
    runner = harness.run_task if runner is None else runner
    floor, entry = harness.CANDIDATES[harness.FLOOR], harness.CANDIDATES[ENTRY]
    failures: dict = {}
    records: dict = {}
    try:
        res = runner(t, [floor, entry, *[_arm(n, failures, records, t.key) for n in GUARDED]], seeds=seeds,
                     data_dir=data_dir)
        res["not_buildable"] = {}
        res["build_records"] = records
        return res
    except Exception:
        if not failures:
            raise
    records = {}
    parts, not_buildable = [runner(t, [floor, entry], seeds=seeds, data_dir=data_dir)], {}
    for n in GUARDED:
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
    tasks = nt.all_tasks()
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
    res["task_source"] = "experiments/harness/n46_embryo_g_arms.py"
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
