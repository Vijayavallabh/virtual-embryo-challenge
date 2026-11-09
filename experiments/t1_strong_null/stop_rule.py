#!/usr/bin/env python
"""Row T12 of experiments/harness/TRIES.md: the stop rule of the strong-null down-tail prior, applied to the harness's
run of its downside measurement. Pure arithmetic on the harness's own result files (no data file is read, no score is
computed here); nothing is written if the results are inconsistent (mixed commits, a dirty tree, other seeds;
vec.harness.check_results).

The downside measurement (the row's held-out evidence; `experiments/harness/run.py run --board T1:test --candidates
floor,b0,b0_scope,b0_strong_null,scope_strong_null`, 10 seeds, the harness's bake scoring):
  - b0_strong_null against b0 on the two E9.5 -> E10.5 tasks (the prior with its table calibrated on the other
    transition, E8.5 -> E9.5, the only legal one on them);
  - scope_strong_null against b0_scope on T1/E8.5->E9.5 (no table is legal there, so the same realisation with the genes
    chosen by the reference level alone: the prior's cost against a truth where highly expressed genes stay flat).

The stop rule, fixed in the row before any of this runs, in T1 board points (mean over the 10 seeds of the paired
difference, change minus base):
  S1: any task's mean difference below LOSS_BOUND, -0.66 (twice the approved price of -0.33 T1, the point where the
      bet's break-even moves from about 16% to about 27%, the top of the 0.2-0.3 range sheet C gave him); or
  S2: on the E8.5 -> E9.5 control (the one weak-null truth), the distribution and covariance groups together (mmd_u,
      weight 0.30; variogram, 0.20) losing more than -DIST_BOUND = 0.33 points on average, the whole approved price
      spent on terms the price did not cover.
If S1 or S2 holds, the verdict is STOP: the build stops and the numbers go to the captain. Otherwise the verdict is
"proceed": the main session may pin the table and flip the option on at a later build commit (the build still needs its
own target-free checks; apply_down_prior's `checks`). The candidates' group margins and metric changes are reported in
any case.

Outputs, in experiments/t1_strong_null/results/stop_rule/ (no test-board count, limit or anchor is read or written):
stop_rule.json and stop_rule.txt, each with a vec.provenance build stamp and the stamp fields of the row (commit,
clean-tree flag, UTC, the result files' sha256); it refuses to overwrite existing outputs (a rerun is a new row).

Usage (anywhere; it only reads the run's JSON files, from the repo root):
  uv run --frozen --no-sync python -I -B experiments/t1_strong_null/stop_rule.py --runs runs/harness/t12_v1
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, provenance
from vec.adoption import SEEDS
from vec.boards import GROUPS

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "stop_rule"
LOSS_BOUND = -0.66  # S1, T1 board points (twice the approved price of -0.33)
DIST_BOUND = -0.33  # S2, T1 board points of the distribution and covariance groups together, on the control task
CONTROL_TASK = "T1/E8.5->E9.5"
PAIRS = {"T1/E9.5->E10.5": ("b0", "b0_strong_null"), "T1/E9.5->E10.5:chamber": ("b0", "b0_strong_null"),
         CONTROL_TASK: ("b0_scope", "scope_strong_null")}
WEIGHTS = {name: weight for name, weight, _ in GROUPS["T1"]}
DISTRIBUTIONAL = ("distribution", "covariance")


def paired(result: Mapping, base: str, change: str) -> dict:
    """One task's paired differences (change minus base), per seed: the board score, and each scored group in T1 board
    points (100 x its weight x the difference of its group skill)."""
    rows = [r for r in result["rows"] if base in r["scores"] and change in r["scores"]]
    score = np.array([float(r["scores"][change]) - float(r["scores"][base]) for r in rows], dtype=np.float64)
    groups = {g: np.array([100.0 * w * (float(r["groups"][change][g]) - float(r["groups"][base][g])) for r in rows],
                          dtype=np.float64) for g, w in WEIGHTS.items()}
    metrics = {m: float(np.mean([float(r["metrics"][change][m]) - float(r["metrics"][base][m]) for r in rows]))
               for m in ("de_score", "de_direction", "mmd_u", "variogram")}
    return {"base": base, "change": change, "seeds": [int(r["seed"]) for r in rows], "n_seeds": len(rows),
            "mean": float(score.mean()) if len(rows) else float("nan"),
            "sd": float(score.std(ddof=1)) if len(rows) > 1 else 0.0,
            "ahead": int((score > 0).sum()), "per_seed": [float(x) for x in score],
            "groups_mean": {g: float(v.mean()) for g, v in groups.items()},
            "groups_per_seed": {g: [float(x) for x in v] for g, v in groups.items()},
            "metric_differences_mean": metrics}


def verdict(tasks: Mapping[str, Mapping]) -> dict:
    """T12's stop rule (module docstring) from `paired` per task; a task absent or without the full 10 seeds means the
    measurement is incomplete, and the verdict is STOP."""
    problems, s1, s2 = [], [], None
    for key in PAIRS:
        t = tasks.get(key)
        if t is None or t["n_seeds"] != len(SEEDS):
            problems.append(f"{key}: incomplete (needs all {len(SEEDS)} seeds)")
            continue
        if t["mean"] < LOSS_BOUND:
            s1.append(f"{key}: mean {t['mean']:+.3f} below {LOSS_BOUND}")
    control = tasks.get(CONTROL_TASK)
    if control is not None and control["n_seeds"] == len(SEEDS):
        s2 = float(sum(control["groups_mean"][g] for g in DISTRIBUTIONAL))
    stop = bool(s1) or bool(problems) or (s2 is not None and s2 < DIST_BOUND)
    return {"verdict": "STOP" if stop else "proceed", "S1_failures": s1, "problems": problems,
            "control_distributional_points": s2, "loss_bound": LOSS_BOUND, "distributional_bound": DIST_BOUND}


def report_text(res: dict) -> str:
    s = res["stamp"]
    lines = [f"T12 stop rule: {res['verdict']['verdict']}. S1 (any task's mean difference below {LOSS_BOUND}): "
             f"{res['verdict']['S1_failures'] or 'none'}. S2 (control's distribution + covariance points below "
             f"{DIST_BOUND}): {res['verdict']['control_distributional_points']}. "
             f"Problems: {res['verdict']['problems'] or 'none'}",
             f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}", ""]
    for key, t in res["tasks"].items():
        lines.append(f"{key}: {t['change']} minus {t['base']}: mean {t['mean']:+.3f} (sd {t['sd']:.3f}, ahead in "
                     f"{t['ahead']} of {t['n_seeds']}); groups " + ", ".join(
                         f"{g} {v:+.3f}" for g, v in t["groups_mean"].items()) + "; metrics " + ", ".join(
                         f"{m} {v:+.4f}" for m, v in t["metric_differences_mean"].items()))
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", required=True, help="the harness run's folder (one JSON per pseudo-task)")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("T12's stop rule has already written its outputs; a rerun needs a new TRIES row")
    runs = Path(a.runs)
    results = []
    files = {}
    for p in sorted(runs.glob("*.json")):
        r = json.loads(p.read_text())
        if isinstance(r, dict) and "task" in r and "rows" in r and r["task"]["test_board"] == "T1:test":
            results.append(r)
            files[p.name] = provenance._sha256_size(p)[0]
    configured = harness.load_tasks()["T1:test"]
    harness.check_results(results, configured, harness._sha256(harness.CONFIG))
    by_key = {r["task"]["key"]: r for r in results}
    tasks = {k: paired(by_key[k], *PAIRS[k]) for k in PAIRS if k in by_key}
    decision = verdict(tasks)
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    res = {"row": "T12 (experiments/harness/TRIES.md)", "verdict": decision, "tasks": tasks,
           "inputs_sha256": files, "runs": str(runs),
           "stamp": {"commit": commit, "tree_clean": unclean == [],
                     "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    (out / "stop_rule.json").write_text(json.dumps(res, indent=1) + "\n")
    provenance.write_stamp(out / "stop_rule.json")
    text = report_text(res)
    (out / "stop_rule.txt").write_text(text)
    provenance.write_stamp(out / "stop_rule.txt")
    print(text, end="")


if __name__ == "__main__":
    main()
