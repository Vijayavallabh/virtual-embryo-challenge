#!/usr/bin/env python
"""The embryo mirror-only partner correction scored on the nine practice tasks (campaign branch campaign/emb-mirror;
the experiments/harness/TRIES.md row is registered by the main session before any run).

Why. The ledger's "Pruned before running" holds a mirror-only correction of FGW's partners (keep FGW's moves; put
each partner on its anchor's side of the symmetric direction), pruned by analogy with try E11 and never read. Row E10
measured the defect: FGW's one-to-one partners sit at random or mirrored across the body's symmetric direction, so
the entry's moves pull both brackets towards its mid-plane and flatten the cloud. This line keeps FGW's own partners
and only moves a partner on the wrong side to its mirror position (vec.t2_build.mirror_fixed_partners). No constant
is new. Fix E2 adds the gated form: the same correction applied to a stage only where the smallest eigenvalue of its
pairs' symmetrised whitened cross-correlation is below 0 (row E10's sign of a mirrored direction), since the ungated
sign test also mirrors a share of partners that are right and only noisy.

Arms:
  fgw_mirrorfix         the entry with the mirror-only correction of FGW's partners
  fgw_mirrorfix_series  the second upload's expression (series) on FGW's corrected partners
  fgw_mirrorgate        the entry with the correction applied only to a stage whose smallest eigenvalue is below 0
  fgw_mirrorgate_series the second upload's expression (series) on FGW's partners with the gated correction

Tasks: the nine of experiments/harness/nine_tasks.py (half, late, new1, new2, e75, e725, x1, x2, x3), each with exactly
the inputs nine_tasks.all_tasks() gives it, 10 seeds (vec.adoption.SEEDS). Only the floor and the four arms are scored;
fgw and fgw_aligned_series are NOT rescored. The task file is the harness's own result (vec.harness.run_task, which
records the commit and the tree state), so nine_tasks._load reads it. This script reads nothing: the reading is a
separate script on main, registered in the ledger row by the main session.

An arm whose build cannot be made on a task (for example the builder's own cloud guard refusing the cloud) is
recorded in the task file under "not_buildable" with the seed and the exception's message, every number masked,
and left out of the rows; the other arms and the floor are still scored (the per-arm guard, the merge and the record
are n36_trend_second_upload's own, imported unchanged), as in n41_polar_arms.py.

  run   --task NAME --out DIR   score the floor and the four arms on ONE task; start the nine together so that they
                                record one commit; refuses an existing task file
"""

from __future__ import annotations

import argparse
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

NEW = ("fgw_mirrorfix", "fgw_mirrorfix_series", "fgw_mirrorgate", "fgw_mirrorgate_series")
ARMS = (harness.FLOOR, *NEW)
NOT_RESCORED = ("fgw", nt.RA)


def score_task(t: harness.PseudoTask, data_dir: Path, runner=None, seeds=SEEDS) -> dict:
    """The floor and the four arms on t (n36's score_task with this script's arms). A build that raises stops neither
    the other arms nor the floor: the task is then scored arm by arm, and an arm that fails is recorded under
    "not_buildable" (the seed and the exception's message, every number masked) and left out of the rows. A failure
    outside an arm's build (data, floor, scoring) stops the task."""
    runner = harness.run_task if runner is None else runner
    floor = harness.CANDIDATES[harness.FLOOR]
    failures: dict = {}
    try:
        res = runner(t, [floor, *[n36._guarded(harness.CANDIDATES[n], failures, t.key) for n in NEW]], seeds=seeds,
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
            parts.append(runner(t, [floor, n36._guarded(harness.CANDIDATES[n], fails, t.key)], seeds=seeds,
                                data_dir=data_dir))
        except Exception:
            if n not in fails:
                raise
            not_buildable[n] = fails[n]
    if not parts:
        parts.append(runner(t, [floor], seeds=seeds, data_dir=data_dir))
    res = n36._merge(parts)
    res["not_buildable"] = not_buildable
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
    res["task_source"] = "experiments/harness/n42_mirrorfix_arms.py"
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
