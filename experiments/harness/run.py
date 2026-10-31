#!/usr/bin/env python
"""Pseudo-task harness runner (feat-008; vec.harness). DGX for `run` (the stages live under data/raw), anywhere for
`decide` (it reads only the run's aggregate JSONs).

  run     score candidates on a test board's pseudo-tasks (all of them, or --task KEY ...), seeds 0-9, writing one JSON
          per pseudo-task (per-seed scores, skills, anchors, cell counts, provenance) into --out
  decide  apply the adoption rule (vec.adoption.keep_change) to each --change against --base over the board's
          pseudo-tasks found in --runs; writes decisions.{json,txt} into --out and, with --register, one row per
          pseudo-task and candidate into results/registry.csv

Usage (DGX, repo root, BLAS threads capped as docs/COMPUTE.md says; long runs in tmux with a log):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 uv run python experiments/harness/run.py run \
      --board T2:heart:test_extrap --candidates floor,copy,pseudobulk_shift --out runs/harness/<run_id>
  uv run python experiments/harness/run.py decide --runs runs/harness/<run_id> --board T2:heart:test_extrap \
      --base floor --change copy --change pseudobulk_shift --out experiments/harness/results/<run_id> --register
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

from vec import harness, registry
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT


def safe(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", key).strip("_")


def parse_seeds(text: str) -> tuple[int, ...]:
    out = []
    for part in text.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return tuple(out)


def cmd_run(a) -> None:
    tasks = harness.load_tasks()
    if a.board not in tasks:
        sys.exit(f"unknown test board {a.board}; known: {sorted(tasks)}")
    chosen = [t for t in tasks[a.board] if not a.task or t.key in a.task]
    if a.task and len(chosen) != len(set(a.task)):
        sys.exit(f"unknown pseudo-task(s) among {a.task}; known: {[t.key for t in tasks[a.board]]}")
    names = [n.strip() for n in a.candidates.split(",") if n.strip()]
    unknown = [n for n in names if n not in harness.CANDIDATES]
    if unknown:
        sys.exit(f"unknown candidate(s) {unknown}; known: {sorted(harness.CANDIDATES)}")
    cands = [harness.CANDIDATES[n] for n in names]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for t in chosen:
        path = out / f"{safe(t.key)}.json"
        if path.exists() and not a.overwrite:
            sys.exit(f"{path} exists; pass --overwrite to replace it")
        t0 = time.time()
        res = harness.run_task(t, cands, seeds=parse_seeds(a.seeds), data_dir=Path(a.data), build=a.build)
        res["run_id"] = out.name
        res["seconds"] = round(time.time() - t0, 1)
        path.write_text(json.dumps(res, indent=1, default=float))
        print(f"wrote {path} ({res['seconds']} s)", flush=True)


def _summary(r: dict, name: str) -> tuple[float, float]:
    s = np.array([row["scores"][name] for row in r["rows"]], float)
    return float(s.mean()), float(s.std(ddof=1)) if len(s) > 1 else 0.0


def cmd_decide(a) -> None:
    runs = Path(a.runs)
    results = [json.loads(p.read_text()) for p in sorted(runs.glob("*.json")) if p.name != "run.json"]
    results = [r for r in results if r["task"]["test_board"] == a.board]
    if not results:
        sys.exit(f"no pseudo-task results for {a.board} in {runs}")
    decisions = [harness.decide(results, a.base, c) for c in a.change]
    lines = [f"{a.board}: {len(results)} pseudo-task(s) in {runs} (protocol {results[0]['protocol']}, build "
             f"{results[0]['build']}, anchors computed, seeds {results[0]['seeds']})"]
    for r in results:
        lines.append(f"  {r['task']['key']} on {r['board']['key']}: "
                     + "; ".join(f"{n} {m:.2f} +/- {sd:.2f}" for n in r["candidates"] for m, sd in [_summary(r, n)])
                     + f" (commit {r['provenance']['commit'][:7]}, clean {r['provenance']['tree_clean']})")
    for d in decisions:
        lines.append(f"{d['change']} against {d['base']}: {'KEEP' if d['keep'] else 'not kept'} ({d['reason']})")
        for k, v in d["tasks"].items():
            if v.get("complete"):
                lines.append(f"    {k}: paired mean gain {v['mean_gain']:+.2f}, ahead in {v['seeds_ahead']} of "
                             f"{len(SEEDS)} seeds, passes {v['passes']}, vetoes {v['vetoes']}")
            else:
                lines.append(f"    {k}: incomplete {v}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"decisions_{safe(a.board)}.json").write_text(json.dumps(
        {"board": a.board, "runs": str(runs), "decisions": decisions,
         "summaries": {r["task"]["key"]: {n: _summary(r, n) for n in r["candidates"]} for r in results}},
        indent=1, default=float))
    (out / f"decisions_{safe(a.board)}.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    if a.register:
        rows = []
        for r in results:
            for n in r["candidates"]:
                m, sd = _summary(r, n)
                rows.append({
                    "utc_time": r["provenance"]["utc"], "commit": r["provenance"]["commit"],
                    "tree_clean": r["provenance"]["tree_clean"], "source": "experiments/harness/run.py",
                    "run_id": r.get("run_id", runs.name), "test_board": a.board, "pseudo_task": r["task"]["key"],
                    "scoring_board": r["board"]["key"], "candidate": n, "params": json.dumps(r["params"].get(n, {})),
                    "seeds": ",".join(map(str, r["seeds"])), "score_mean": f"{m:.3f}", "score_sd": f"{sd:.3f}",
                    "protocol": r["protocol"], "build": r["build"], "anchors": r["anchors"],
                    "gate": "floor_normalised", "frac": r["frac"],
                    "scores_json": str((runs / f"{safe(r['task']['key'])}.json").relative_to(REPO_ROOT))
                    if (runs / f"{safe(r['task']['key'])}.json").resolve().is_relative_to(REPO_ROOT)
                    else str(runs / f"{safe(r['task']['key'])}.json")})
        print(f"registered {registry.append(rows)} rows in {registry.REGISTRY}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--board", required=True, help="test board key, e.g. T2:heart:test_extrap")
    r.add_argument("--task", action="append", help="pseudo-task key (repeatable); default: all of the board's")
    r.add_argument("--candidates", default="floor,copy", help="comma-separated names from vec.harness.CANDIDATES")
    r.add_argument("--seeds", default="0-9")
    r.add_argument("--build", default=harness.BUILD, choices=("main",), help="candidates are compared under main")
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--out", required=True)
    r.add_argument("--overwrite", action="store_true")
    d = sub.add_parser("decide")
    d.add_argument("--runs", required=True)
    d.add_argument("--board", required=True)
    d.add_argument("--base", default=harness.FLOOR)
    d.add_argument("--change", action="append", required=True)
    d.add_argument("--out", required=True)
    d.add_argument("--register", action="store_true")
    a = ap.parse_args(argv)
    (cmd_run if a.cmd == "run" else cmd_decide)(a)


if __name__ == "__main__":
    main()
