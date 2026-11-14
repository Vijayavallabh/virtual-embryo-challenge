#!/usr/bin/env python
"""The embryo batch J scored on the nine practice tasks (campaign branch campaign/emb-j; the
experiments/harness/TRIES.md row is registered by the main session before any run).

Why. The movement-matched control fgw_damp_m (fix brief M2, campaign/emb-g) passed both keep rules and fgw_partial does
not lead it (tries E42), so by the ledger's rule fgw_damp_m is the arm put forward, on one condition: its amount of
movement comes from cell-type labels (the labels behind fgw_partial's share and so behind the partial coupling that
gives its two factors), and cell-type names are not harmonised across stages. Its twin at the lineage-group level is
built and read first.

Arm:
  fgw_damp_m_grp  fgw_damp_m with the labels behind its two factors the brackets' anchors' obs["group"] in place of
                  obs["celltype"], so the partial coupling that sets them is fgw_partial_grp's (campaign/emb-h); the
                  entry's coupling, directions, partners, rotation, draws, streams and expression are untouched
and fgw, fgw_partial_grp and fgw_damp_m, rescored under their own names, so that the run can be checked against the
committed scores.

Tasks: the nine of experiments/harness/nine_tasks.py (half, late, new1, new2, e75, e725, x1, x2, x3), each with exactly
the inputs nine_tasks.all_tasks() gives it, 10 seeds (vec.adoption.SEEDS). The task file is the harness's own result
(vec.harness.run_task, which records the commit and the tree state), so nine_tasks._load reads it. For fgw_partial_grp,
fgw_damp_m and fgw_damp_m_grp each seed's build record is written under "build_records" exactly as
n46_embryo_g_arms.py writes it (its own record-keeping, imported unchanged): the arm's own record fields, with
"movement_ratio": {"anchor": [b, a], "cell": [b, a]}, the arm's mean anchor shift and mean cell shift per bracket
divided by fgw's on the same seed and task (null where fgw's is 0), and for the two damped arms also "factors"
[before, after]. The raw shifts and fgw's own records are not written. This script reads nothing: the reading is a
separate script on main, registered in the ledger row by the main session.

An arm whose build cannot be made on a task (a bracket without obs["group"], the partial solver failing, or the
builder's own cloud guard refusing the cloud) is recorded in the task file under "not_buildable" with the seed and the
exception's message, every number masked, and left out of the rows; fgw, the other arms and the floor are still scored
(the per-arm guard, the merge and the record are n36_trend_second_upload's own, imported unchanged), as in
n46_embryo_g_arms.py.

  run   --task NAME --out DIR   score the floor, fgw, fgw_partial_grp, fgw_damp_m and the arm on ONE task; start the
                                nine together so that they record one commit; refuses an existing task file
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n36_trend_second_upload as n36  # noqa: E402  (the per-arm guard, the merge and the not-buildable record)
import n46_embryo_g_arms as n46  # noqa: E402  (the build records and their movement ratios)
import nine_tasks as nt  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

NEW = ("fgw_damp_m_grp",)
ENTRY = "fgw"  # rescored under its own name: the run's check against the committed scores
RESCORED = ("fgw_partial_grp", "fgw_damp_m")  # rescored under their own names; guarded, since each can refuse
GUARDED = (*RESCORED, *NEW)
ARMS = (harness.FLOOR, ENTRY, *GUARDED)


def score_task(t: harness.PseudoTask, data_dir: Path, runner=None, seeds=SEEDS) -> dict:
    """The floor, fgw, fgw_partial_grp, fgw_damp_m and the arm on t. A build of a guarded arm that raises stops neither
    the other arms, fgw nor the floor: the task is then scored as the floor with fgw, and the floor with each guarded
    arm, and an arm that fails is recorded under "not_buildable" (the seed and the exception's message, every number
    masked) and left out of the rows. A failure elsewhere (data, floor, fgw, scoring) stops the task."""
    runner = harness.run_task if runner is None else runner
    floor, entry = harness.CANDIDATES[harness.FLOOR], harness.CANDIDATES[ENTRY]
    failures: dict = {}
    records: dict = {}
    try:
        res = runner(t, [floor, n46._recorded(entry, records),
                         *[n46._arm(n, failures, records, t.key) for n in GUARDED]], seeds=seeds, data_dir=data_dir)
        res["not_buildable"] = {}
        res["build_records"] = n46._movement(records)
        return res
    except Exception:
        if not failures:
            raise
    records = {}
    parts, not_buildable = [runner(t, [floor, n46._recorded(entry, records)], seeds=seeds, data_dir=data_dir)], {}
    for n in GUARDED:
        fails: dict = {}
        try:
            parts.append(runner(t, [floor, n46._arm(n, fails, records, t.key)], seeds=seeds, data_dir=data_dir))
        except Exception:
            if n not in fails:
                raise
            not_buildable[n] = fails[n]
            records.pop(n, None)
    res = n36._merge(parts)
    res["not_buildable"] = not_buildable
    res["build_records"] = n46._movement(records)
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
    res["task_source"] = "experiments/harness/n50_embryo_j_arms.py"
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
