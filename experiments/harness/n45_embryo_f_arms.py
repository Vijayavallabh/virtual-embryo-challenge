#!/usr/bin/env python
"""The embryo batch F scored on the nine practice tasks (campaign branch campaign/emb-f; the
experiments/harness/TRIES.md row is registered by the main session before any run).

Why. Two lines of the ledger are closed on a label rather than a reading: row E31 ("Potts sorting ... is not built: no
update rule or constant is on record") and row E32 (a partial FGW coupling, "the share of mass it would leave unmatched
is its one constant and no record gives it"). Each is built here in a form with no free number, so that neither stays
closed on a label: a zero-temperature Potts sorting of which drawn cell sits where (vec.t2_build.potts_arrangement)
and a partial FGW coupling whose transported share is the overlap of the brackets' anchors' celltype distributions
(vec.t2_build.partial_fgw_coupling). The partial form rests on three choices no record fixes: the label level
(celltype), the overlap as the estimate of the share, and the zero displacement of an unmatched anchor.

Arms:
  fgw_potts    the entry with its drawn cells sorted among its own positions by the Potts rule
  fgw_partial  the entry with a partial FGW coupling (the entry's own coupling where the overlap is exactly 1)
and fgw, rescored under its own name, so that the run can be checked against the committed scores.

Tasks: the nine of experiments/harness/nine_tasks.py (half, late, new1, new2, e75, e725, x1, x2, x3), each with exactly
the inputs nine_tasks.all_tasks() gives it, 10 seeds (vec.adoption.SEEDS). The task file is the harness's own result
(vec.harness.run_task, which records the commit and the tree state), so nine_tasks._load reads it. Each new arm's build
record (per seed: the Potts sorting's shares, the number of accepted swaps and its stop; the partial coupling's s, the
transported mass and the share of anchors with no mass per side) is written under "build_records". This script reads
nothing: the reading is a separate script on main, registered in the ledger row by the main session.

An arm whose build cannot be made on a task (the partial solver failing, or the builder's own cloud guard refusing the
cloud) is recorded in the task file under "not_buildable" with the seed and the exception's message, every number
masked, and left out of the rows; fgw, the other arm and the floor are still scored (the per-arm guard, the merge and
the record are n36_trend_second_upload's own, imported unchanged), as in n43_embryo_d_arms.py.

  run   --task NAME --out DIR   score the floor, fgw and the two arms on ONE task; start the nine together so that they
                                record one commit; refuses an existing task file
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

NEW = ("fgw_potts", "fgw_partial")
ENTRY = "fgw"  # rescored under its own name: the run's check against the committed scores
ARMS = (harness.FLOOR, ENTRY, *NEW)


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
    return n36._guarded(_recorded(harness.CANDIDATES[name], records), failures, key)


def score_task(t: harness.PseudoTask, data_dir: Path, runner=None, seeds=SEEDS) -> dict:
    """The floor, fgw and the two arms on t. A build of a new arm that raises stops neither the other arm, fgw nor the
    floor: the task is then scored as the floor with fgw, and the floor with each arm, and an arm that fails is recorded
    under "not_buildable" (the seed and the exception's message, every number masked) and left out of the rows. A
    failure elsewhere (data, floor, fgw, scoring) stops the task."""
    runner = harness.run_task if runner is None else runner
    floor, entry = harness.CANDIDATES[harness.FLOOR], harness.CANDIDATES[ENTRY]
    failures: dict = {}
    records: dict = {}
    try:
        res = runner(t, [floor, entry, *[_arm(n, failures, records, t.key) for n in NEW]], seeds=seeds,
                     data_dir=data_dir)
        res["not_buildable"] = {}
        res["build_records"] = records
        return res
    except Exception:
        if not failures:
            raise
    records = {}
    parts, not_buildable = [runner(t, [floor, entry], seeds=seeds, data_dir=data_dir)], {}
    for n in NEW:
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
    res["task_source"] = "experiments/harness/n45_embryo_f_arms.py"
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
