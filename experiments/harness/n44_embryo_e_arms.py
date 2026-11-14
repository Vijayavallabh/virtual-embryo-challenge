#!/usr/bin/env python
"""The embryo batch E scored on the nine practice tasks (campaign branch campaign/emb-e; the
experiments/harness/TRIES.md row is registered by the main session before any run).

Why. The ledger's "Pruned before running" holds "All boards: within-upload mixtures over shift sizes ... Each is about
0" and "negative-binomial resampling", read on the heart only (tries H16 and H19), and the round's record closed a
quadratic trend as a variant of the line. An independent audit named all three open. They are built here as campaign
arms; no constant is new besides the mixture's half, which is H16's.

Arms:
  fgw_aligned_mix        fgw_aligned_series' cells and positions, a uniformly random half of the rows (floor(n / 2))
                         with fgw_aligned_series' expression and the others with fgw_aligned_trend's
  fgw_aligned_series_mn  the second upload as built (fgw_aligned_series) with a whole-cell multinomial redraw of each
                         rendered cell (try H19's form)
  fgw_quad_render        fgw_trend_render with a least-squares quadratic in place of the trend's line
  fgw_aligned_quad       fgw_aligned_trend with the same quadratic
Rescored here under their own names, so that the run can be checked against the committed scores: fgw,
fgw_aligned_series, fgw_trend_render and fgw_aligned_trend.

Tasks: the nine of experiments/harness/nine_tasks.py (half, late, new1, new2, e75, e725, x1, x2, x3), each with its
trend samples by n36's sample rule (n36_trend_second_upload.all_tasks(): the brackets plus every released series stage
strictly outside the bracket interval that is not the target, as roles trend_<stem>, which the size rule does not read),
so that the trend arms apply; 10 seeds (vec.adoption.SEEDS). The task file is the harness's own result
(vec.harness.run_task, which records the commit and the tree state), so nine_tasks._load reads it. This script reads
nothing: the reading is a separate script, registered in the ledger row by the main session.

An arm whose build cannot be made on a task (the builder's own refusal, for example the quadratic's "at least three
distinct days" or the cloud guard) is recorded in the task file under "not_buildable" with the seed and the exception's
message, every number masked, and left out of the rows; the other arms and the floor are still scored (the per-arm
guard, the merge and the record are n36_trend_second_upload's own, imported unchanged), as in n43_embryo_d_arms.py.

  run   --task NAME --out DIR   score the floor, the four arms and the four rescored arms on ONE task; start the nine
                                together so that they record one commit; refuses an existing task file
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n36_trend_second_upload as n36  # noqa: E402  (the sample rule, the per-arm guard, the merge, the record)
import nine_tasks as nt  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

NEW = ("fgw_aligned_mix", "fgw_aligned_series_mn", "fgw_quad_render", "fgw_aligned_quad")
RESCORED = ("fgw", nt.RA, "fgw_trend_render", "fgw_aligned_trend")
SCORED = (*NEW, *RESCORED)
ARMS = (harness.FLOOR, *SCORED)


def all_tasks() -> dict:
    """The nine tasks with their trend samples (n36's sample rule)."""
    return n36.all_tasks()


def score_task(t: harness.PseudoTask, data_dir: Path, runner=None, seeds=SEEDS) -> dict:
    """The floor and the eight arms on t (n43's score_task with this script's arms). A build that raises stops neither
    the other arms nor the floor: the task is then scored arm by arm, and an arm that fails is recorded under
    "not_buildable" (the seed and the exception's message, every number masked) and left out of the rows. A failure
    outside an arm's build (data, floor, scoring) stops the task."""
    runner = harness.run_task if runner is None else runner
    floor = harness.CANDIDATES[harness.FLOOR]
    failures: dict = {}
    try:
        res = runner(t, [floor, *[n36._guarded(harness.CANDIDATES[n], failures, t.key) for n in SCORED]],
                     seeds=seeds, data_dir=data_dir)
        res["not_buildable"] = {}
        return res
    except Exception:
        if not failures:
            raise
    parts, not_buildable = [], {}
    for n in SCORED:
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
    tasks = all_tasks()
    if a.task not in tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(tasks)}")
    t = tasks[a.task]
    missing = [n for n in SCORED if not harness.CANDIDATES[n].applies(t)]
    if missing:
        sys.exit(f"{t.key}: arms that do not apply: {missing}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{sw.safe(t.key)}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = score_task(t, Path(a.data))
    res["task_source"] = "experiments/harness/n44_embryo_e_arms.py"
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
