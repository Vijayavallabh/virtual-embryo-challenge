#!/usr/bin/env python
"""Pseudo-task harness runner (feat-008; vec.harness). DGX for `run` (the stages live under data/raw), anywhere for
`decide` (it reads only the run's aggregate JSONs).

  run     score candidates on a test board's pseudo-tasks (all of them, or --task KEY ...), seeds 0-9, writing one JSON
          per pseudo-task (per-seed scores, skills, anchors, cell counts, provenance) into --out
  decide  apply the adoption rule (vec.adoption.keep_change) to each --change against --base over the board's
          pseudo-tasks found in --runs; writes decisions.{json,txt} into --out and, with --register, one row per
          pseudo-task and candidate into results/registry.csv. On T3 each decision also carries the fit checks of
          tries K2a and K2b (vec.harness.t3_fit_checks: the margin by scored group, the top 10 genes' share of the
          de_direction margin, the seeds where the severity gate differs)

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

from vec import harness, latent_extrap, registry
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT


def _rel(p: Path) -> str:
    p = p.resolve()
    return str(p.relative_to(REPO_ROOT)) if p.is_relative_to(REPO_ROOT) else str(p)


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
    # --task takes a key or its shell-safe form (T1_E9.5_E10.5 for T1/E9.5->E10.5: '>' would redirect in a shell)
    chosen = [t for t in tasks[a.board] if not a.task or t.key in a.task or safe(t.key) in a.task]
    if a.task and len(chosen) != len(set(a.task)):
        sys.exit(f"unknown pseudo-task(s) among {a.task}; known: {[t.key for t in tasks[a.board]]}")
    names = [n.strip() for n in a.candidates.split(",") if n.strip()]
    unknown = [n for n in names if n not in harness.CANDIDATES]
    if unknown:
        sys.exit(f"unknown candidate(s) {unknown}; known: {sorted(harness.CANDIDATES)}")
    cands = [harness.CANDIDATES[n] for n in names]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    failed = []
    for t in chosen:
        path = out / f"{safe(t.key)}.json"
        if path.exists() and not a.overwrite:
            sys.exit(f"{path} exists; pass --overwrite to replace it")
        t0 = time.time()
        try:
            res = harness.run_task(t, cands, seeds=parse_seeds(a.seeds), data_dir=Path(a.data), build=a.build)
        except Exception as e:  # noqa: BLE001 (recorded, so decide sees the pseudo-task as missing, not absent)
            failed.append(t.key)
            (out / f"{safe(t.key)}.failed.txt").write_text(f"{type(e).__name__}: {e}\n")
            print(f"FAILED {t.key}: {type(e).__name__}: {e}", flush=True)
            continue
        res["run_id"] = out.name
        res["seconds"] = round(time.time() - t0, 1)
        path.write_text(json.dumps(res, indent=1, default=float))
        print(f"wrote {path} ({res['seconds']} s)", flush=True)
        latent_extrap.write_records(out / "latent_records" / f"{safe(t.key)}.json", t.key)  # latent arms only
    if failed:
        sys.exit(f"{len(failed)} pseudo-task(s) failed: {failed}")


def _summary(r: dict, name: str) -> tuple[float, float]:
    s = np.array([row["scores"][name] for row in r["rows"]], float)
    return float(s.mean()), float(s.std(ddof=1)) if len(s) > 1 else 0.0


def cmd_decide(a) -> None:
    runs = Path(a.runs).resolve()
    configured = harness.load_tasks().get(a.board)
    if not configured:
        sys.exit(f"unknown test board {a.board}")
    results = []
    for p in sorted(runs.glob("*.json")):
        r = json.loads(p.read_text())
        if isinstance(r, dict) and "task" in r and "rows" in r and r["task"]["test_board"] == a.board:
            results.append(r)
    if not results:
        sys.exit(f"no pseudo-task results for {a.board} in {runs}")
    sha = harness._sha256(harness.CONFIG)
    decisions = [harness.decide(results, configured, a.base, c, config_sha256=sha) for c in a.change]
    for d in decisions:  # T3 rows carry the parts of the response metrics (tries K2a and K2b's fit checks)
        d["fit_checks"] = {r["task"]["key"]: fc for r in results
                           if (fc := harness.t3_fit_checks(r, a.base, d["change"])) is not None}
    lines = [f"{a.board}: {len(results)} pseudo-task(s) in {runs} (protocol {results[0]['protocol']}, build "
             f"{results[0]['build']}, anchors computed, seeds {results[0]['seeds']})"]
    for r in results:
        lines.append(f"  {r['task']['key']} on {r['board']['key']}: "
                     + "; ".join(f"{n} {m:.2f} +/- {sd:.2f}" for n in r["candidates"] for m, sd in [_summary(r, n)])
                     + f" (commit {r['provenance']['commit'][:7]}, clean {r['provenance']['tree_clean']})")
    for d in decisions:
        lines.append(f"{d['change']} against {d['base']}: {'KEEP' if d['keep'] else 'not kept'} ({d['reason']}; "
                     f"decided on {d['decided_on']}; not counted: {d['excluded'] or 'none'})")
        for k, v in d["tasks"].items():
            if v.get("complete"):
                lines.append(f"    {k}: paired mean gain {v['mean_gain']:+.2f}, ahead in {v['seeds_ahead']} of "
                             f"{len(SEEDS)} seeds, passes {v['passes']}, vetoes {v['vetoes']}")
            else:
                lines.append(f"    {k}: incomplete {v}")
        for k, fc in d["fit_checks"].items():
            share = "n/a" if fc["top_share"] is None else f"{fc['top_share']:.2f}"
            lines.append(f"    {k} fit checks: margin by group "
                         + ", ".join(f"{g} {m:+.2f}" for g, m in fc["group_margins"].items())
                         + f"; de_direction margin {fc['de_direction_margin']:+.4f}, its top {fc['top']} genes' share "
                         f"{share}; severity gate differs in {len(fc['severity_gate_differs'])} of {len(fc['seeds'])} "
                         f"seeds" + ("" if fc["matches_panel"] else "; WARNING: the parts do not reproduce the panel"))
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
                    "scores_json": _rel(runs / f"{safe(r['task']['key'])}.json")})
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
