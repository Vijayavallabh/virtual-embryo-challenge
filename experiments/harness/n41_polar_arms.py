#!/usr/bin/env python
"""The embryo polar-move arms scored on the nine practice tasks (campaign branch campaign/emb-polar; the
experiments/harness/TRIES.md row is registered by the main session before any run).

Why. experiments/research_2026_11_runs/REPORT.md (recipe R10) records that linear interpolation shrinks parts that
rotate; the remedy factors each local map A = R S and moves along R(t)((1 - t) I + t S). The entry moves every cell by
the mean displacement of its 8 nearest anchors; the polar move replaces that translation by the same anchors' local
linear map, interpolated along its polar factors (vec.t2_build.polar_moves). No constant is new.

Arms:
  fgw_polar                       the entry with the polar move on FGW's own partners
  fgw_aligned_polar_series        the second upload as built (fgw_aligned_series) with the polar move
  fgw_polar_frame                 the entry with the polar move in its frame form (fix C2: the patch's whole affine
                                  map of the common frame interpolated, translation included)
  fgw_aligned_polar_frame_series  the second upload as built with the frame form

Tasks: the nine of experiments/harness/nine_tasks.py (half, late, new1, new2, e75, e725, x1, x2, x3), each with exactly
the inputs nine_tasks.all_tasks() gives it, 10 seeds (vec.adoption.SEEDS). Only the floor and the four arms are scored;
fgw and fgw_aligned_series are NOT rescored. The task file is the harness's own result (vec.harness.run_task, which
records the commit and the tree state), so nine_tasks._load reads it. This script reads nothing: the reading is a
separate script on main, registered in the ledger row by the main session.

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
import nine_tasks as nt  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

NEW = ("fgw_polar", "fgw_aligned_polar_series", "fgw_polar_frame", "fgw_aligned_polar_frame_series")
ARMS = (harness.FLOOR, *NEW)
NOT_RESCORED = ("fgw", nt.RA)


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
    res = harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data))
    res["task_source"] = "experiments/harness/n41_polar_arms.py"
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
