#!/usr/bin/env python
"""Rows T12a and T12b of experiments/harness/TRIES.md: the stop rule of the strong-null down-tail prior, applied to the
harness's run of its downside measurement. Pure arithmetic on the harness's own result files (no data file is read, no
score is computed here); nothing is written if the results are inconsistent (mixed commits, a dirty tree, other seeds;
vec.harness.check_results).

The downside measurement (the rows' held-out evidence; `experiments/harness/run.py run --board T1:test --candidates
floor,b0,b0_scope,b0_strong_null,scope_strong_null`, 10 seeds, the harness's bake scoring). Both candidates choose genes
as the test build does, by the reference's level alone (`vec.t1_scope.DownPriorTable.level_only`, nothing fitted):
  - T12a: b0_strong_null against b0 on the two E9.5 -> E10.5 tasks (a strong-null truth);
  - T12b: scope_strong_null against b0_scope on T1/E8.5->E9.5, the control (a weak-null truth: the prior's cost when
    highly expressed genes stay flat).

The stop rule, fixed in the rows before any of this runs (the captain's rulings of 2026-11-10 (a) and (d): the bet goes
into T1's files only if its measured cost comes in near the -0.33 he approved, otherwise the numbers go back to him
before any file changes). Means over the 10 seeds of the paired difference, change minus base:
  S1 (any task; T1 board points): the total difference, or the distribution and covariance groups together (mmd_u,
      weight 0.30; variogram, 0.20), below LOSS_BOUND, -0.66 (twice the approved price; break-even moves from about 16%
      to about 27%). The total alone is not enough: on a strong-null truth a de_score gain can hide a distribution loss.
  S2 (the control; raw metric units, because points depend on a skill map whose slope is the operating point's: one
      raw unit of mmd_u costs -52.4 points at b0_scope's operating point on the control, -68 to -94 near the floor for
      d_f 0.11 to 0.08, -312 near the floor at the control's own d_f of 0.024, and the test's d_f is hidden): the raw
      mmd_u increase above RAW_MMD_BOUND (0.0035, which is 0.33 points near the floor at d_f 0.08), the raw variogram
      increase above RAW_VARIOGRAM_BOUND (0.0004, which is 0.36 points near the floor at the E9.5 -> E10.5 tasks' d_f of
      0.0055), or the two groups together below DIST_BOUND, -0.33 points (the rule as first written, kept so nothing is
      loosened).
  Verdicts: "STOP" if S1 or S2 holds, or the run is incomplete or has a non-finite mean (the build stops, the numbers go
  to the captain). "captain" if no stop holds but any task's total or distribution-plus-covariance points lie below
  APPROVED_LOSS (-0.33): between -0.33 and -0.66 the build proceeds only on his explicit reading. "proceed" otherwise:
  the main session may set the approval switch (configs/t1_strong_null.yaml) at a later build commit, the build still
  needs its own target-free checks (apply_down_prior's `checks`) and the change still needs his ruling.
The implied points of each task's raw increases are reported for several values of d_f (near-floor slopes -7.5 / d_f for
mmd_u and -5 / d_f for the variogram), since the test's d_f is unknown. The candidates' group margins and metric
changes are reported in any case.

Outputs, in experiments/t1_strong_null/results/stop_rule/ (no test-board count, limit or anchor is read or written):
stop_rule.json and stop_rule.txt, each with a vec.provenance build stamp and the stamp fields of the rows (commit,
clean-tree flag, UTC, the result files' sha256); it refuses to overwrite existing outputs (a rerun is a new row).

Usage (anywhere; it only reads the run's JSON files, from the repo root):
  uv run --frozen --no-sync python -I -B experiments/t1_strong_null/stop_rule.py --runs runs/harness/t12_v1
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, provenance, t1_scope
from vec.adoption import SEEDS
from vec.boards import GROUPS

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "stop_rule"
LOSS_BOUND = -0.66  # S1, T1 board points (twice the approved price of -0.33): below it, STOP
APPROVED_LOSS = -0.33  # the approved price: a measured loss below this and not below LOSS_BOUND asks the captain
DIST_BOUND = -0.33  # S2, the control's distribution and covariance points (the rule as first written, kept)
RAW_MMD_BOUND = t1_scope.DP_MMD_BOUND  # S2, the control's raw mmd_u increase (0.33 points near the floor at d_f 0.08)
RAW_VARIOGRAM_BOUND = 0.0004  # S2, the control's raw variogram increase (0.36 points near the floor at d_f 0.0055)
CONTROL_TASK = "T1/E8.5->E9.5"
PAIRS = {"T1/E9.5->E10.5": ("b0", "b0_strong_null"), "T1/E9.5->E10.5:chamber": ("b0", "b0_strong_null"),
         CONTROL_TASK: ("b0_scope", "scope_strong_null")}
WEIGHTS = {name: weight for name, weight, _ in GROUPS["T1"]}
DISTRIBUTIONAL = ("distribution", "covariance")
DF_MMD = (0.024, 0.05, 0.08, 0.11)  # d_f values shown for mmd_u: the control's, two between, the E9.5 -> E10.5 tasks'
DF_VARIOGRAM = (0.0009, 0.003, 0.0055, 0.01)  # the same for the variogram: the control's, one between, E9.5 -> E10.5's


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


def distributional_points(task: Mapping) -> float:
    """A task's mean paired T1 points of the distribution and covariance groups together."""
    return float(sum(task["groups_mean"][g] for g in DISTRIBUTIONAL))


def implied_points(raw_mmd: float, raw_variogram: float) -> dict:
    """The T1 points a raw increase would cost near the floor, for several d_f (the skill map's slope there: -7.5 / d_f
    per unit of mmd_u, -5 / d_f per unit of the variogram; the test's d_f is unknown)."""
    return {"mmd_u": {str(d): -7.5 / d * raw_mmd for d in DF_MMD},
            "variogram": {str(d): -5.0 / d * raw_variogram for d in DF_VARIOGRAM}}


def _finite(task: Mapping) -> bool:
    numbers = [task["mean"], *task["groups_mean"].values(), *task["metric_differences_mean"].values()]
    return all(math.isfinite(float(x)) for x in numbers)


def verdict(tasks: Mapping[str, Mapping]) -> dict:
    """The rows' stop rule (module docstring) from `paired` per task. A task absent or without the full 10 seeds, or
    with a non-finite mean or group mean, means the measurement is incomplete, and the verdict is STOP."""
    problems, s1, s2, band = [], [], [], []
    for key in PAIRS:
        t = tasks.get(key)
        if t is None or t["n_seeds"] != len(SEEDS):
            problems.append(f"{key}: incomplete (needs all {len(SEEDS)} seeds)")
            continue
        if not _finite(t):
            problems.append(f"{key}: a non-finite mean or group mean")
            continue
        for what, value in (("total", t["mean"]), ("distribution+covariance", distributional_points(t))):
            if value < LOSS_BOUND:
                s1.append(f"{key}: {what} {value:+.3f} below {LOSS_BOUND}")
            elif value < APPROVED_LOSS:
                band.append(f"{key}: {what} {value:+.3f} below {APPROVED_LOSS}")
    control, raw, dist = tasks.get(CONTROL_TASK), None, None
    if control is not None and control["n_seeds"] == len(SEEDS) and _finite(control):
        dist = distributional_points(control)
        m = control["metric_differences_mean"]
        raw = {"mmd_u": m["mmd_u"], "variogram": m["variogram"]}
        if m["mmd_u"] > RAW_MMD_BOUND:
            s2.append(f"control: raw mmd_u increase {m['mmd_u']:+.5f} above {RAW_MMD_BOUND}")
        if m["variogram"] > RAW_VARIOGRAM_BOUND:
            s2.append(f"control: raw variogram increase {m['variogram']:+.5f} above {RAW_VARIOGRAM_BOUND}")
        if dist < DIST_BOUND:
            s2.append(f"control: distribution+covariance {dist:+.3f} below {DIST_BOUND}")
    stop = bool(s1) or bool(s2) or bool(problems)
    return {"verdict": "STOP" if stop else ("captain" if band else "proceed"), "S1_failures": s1, "S2_failures": s2,
            "captain_band": band, "problems": problems, "control_distributional_points": dist,
            "control_raw_increase": raw, "loss_bound": LOSS_BOUND, "approved_loss": APPROVED_LOSS,
            "distributional_bound": DIST_BOUND, "raw_mmd_bound": RAW_MMD_BOUND,
            "raw_variogram_bound": RAW_VARIOGRAM_BOUND}


def report_text(res: dict) -> str:
    s, v = res["stamp"], res["verdict"]
    lines = [f"T12a/T12b stop rule: {v['verdict']}. S1 (a task's total or distribution+covariance points below "
             f"{LOSS_BOUND}): {v['S1_failures'] or 'none'}. S2 (control: raw mmd_u above {RAW_MMD_BOUND}, raw "
             f"variogram above {RAW_VARIOGRAM_BOUND}, distribution+covariance points below {DIST_BOUND}): "
             f"{v['S2_failures'] or 'none'}. Between {APPROVED_LOSS} and {LOSS_BOUND} (only on the captain's explicit "
             f"reading): {v['captain_band'] or 'none'}. Problems: {v['problems'] or 'none'}",
             f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}", ""]
    for key, t in res["tasks"].items():
        lines.append(f"{key}: {t['change']} minus {t['base']}: mean {t['mean']:+.3f} (sd {t['sd']:.3f}, ahead in "
                     f"{t['ahead']} of {t['n_seeds']}); distribution+covariance {distributional_points(t):+.3f}; "
                     "groups "
                     + ", ".join(f"{g} {v2:+.3f}" for g, v2 in t["groups_mean"].items()) + "; metrics " + ", ".join(
                         f"{m} {v2:+.5f}" for m, v2 in t["metric_differences_mean"].items()))
        imp = implied_points(t["metric_differences_mean"]["mmd_u"], t["metric_differences_mean"]["variogram"])
        lines.append("  raw increases' implied points near the floor: mmd_u " + ", ".join(
            f"d_f {d}: {p:+.2f}" for d, p in imp["mmd_u"].items()) + "; variogram " + ", ".join(
            f"d_f {d}: {p:+.2f}" for d, p in imp["variogram"].items()))
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
    implied = {k: implied_points(t["metric_differences_mean"]["mmd_u"], t["metric_differences_mean"]["variogram"])
               for k, t in tasks.items()}
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    res = {"rows": "T12a and T12b (experiments/harness/TRIES.md)", "verdict": decision, "tasks": tasks,
           "implied_points_near_floor": implied, "inputs_sha256": files, "runs": str(runs),
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
