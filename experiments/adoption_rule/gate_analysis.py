#!/usr/bin/env python
"""How often does an adoption rule keep a change that does nothing? (2026-11-01, before any release-data adoption)

The plan keeps a model change if it "improves >= 2 of 3 leave-one-transition-out pseudo-tasks" (ASSUMPTIONS.md A17;
docs/PLAN.md section 4.1; feat-015), with paired comparisons over 10 seeds (feat-008). With no margin, a change with no
effect is ahead on each pseudo-task half the time, so it passes 2 of 3 half the time, however many seeds are averaged.
A red-team agent flagged this on 2026-10-31 and a side agent raised it with the captain on 2026-11-01.

This script (1) measures the seed-to-seed noise of a paired score difference (same seeds, same working copies) in the
team's committed results on released data and validation targets, and (2) simulates the chance that each candidate
rule adopts a change, for true per-pseudo-task gains from -0.5 to +2 board points, with the transitions allowed to
disagree (true gain varying across pseudo-tasks with sd tau). Rules:
  current   paired mean gain > 0 on at least 2 of 3 pseudo-tasks (all of them when a board has fewer than 3)
  rule7     the approved count rule's test applied per pseudo-task (paired mean gain > 0.3 points and ahead in at least
            8 of 10 seeds), on at least 2 of 3 pseudo-tasks (all when fewer than 3), and no pseudo-task with a paired
            mean loss greater than 0.3 points
Reads only committed aggregate tables; writes experiments/adoption_rule/results/gate_analysis.{json,txt}.
"""

from __future__ import annotations

import csv
import json
from itertools import combinations
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "experiments" / "adoption_rule" / "results"
SEEDS = 10
MARGIN, MIN_AHEAD = 0.3, 8


def paired_sds() -> dict[str, list[float]]:
    """sd over seeds of the paired score difference between two candidates scored on the same seeds."""
    out: dict[str, list[float]] = {}
    # density on the validation targets: two count levels of the same sampling pattern (the condition name carries
    # the level, so pairs are formed within target and pattern)
    rows = list(csv.DictReader(open(REPO / "experiments/density/results/release/results_long.csv")))
    by = {}
    for r in rows:
        by.setdefault((r["target"], r["pattern"], r["level"]), {})[int(r["seed"])] = float(r["score"])
    sds = []
    for (t, c, lvl), s in by.items():
        for (t2, c2, lvl2), s2 in by.items():
            if (t, c) == (t2, c2) and float(lvl2) > float(lvl):
                common = sorted(set(s) & set(s2))
                if len(common) >= 5:
                    d = np.array([s2[k] - s[k] for k in common])
                    sds.append(float(d.std(ddof=1)))
    out["density levels (validation targets)"] = sds
    # T1 scope pseudo-task on released stages: every pair of candidates
    for lib in ("pseudo_E8.5_E9.5_lib1", "pseudo_E8.5_E9.5_lib2"):
        rows = list(csv.DictReader(open(REPO / f"experiments/t1_scope/results/{lib}/scores.csv")))
        cand = {}
        for r in rows:
            cand.setdefault(r["candidate"], {})[int(r["seed"])] = float(r["score"])
        sds = []
        for a, b in combinations(sorted(cand), 2):
            common = sorted(set(cand[a]) & set(cand[b]))
            if len(common) >= 5:
                sds.append(float(np.std([cand[b][k] - cand[a][k] for k in common], ddof=1)))
        out[f"T1 pseudo-task candidate pairs ({lib})"] = sds
    # calibration: an organiser baseline against the floor (copy_last scores exactly 50 at every seed)
    sds = []
    for p in sorted((REPO / "experiments/calibration/results/boards").glob("*.json")):
        res = json.loads(p.read_text())
        for model in ("pseudobulk_shift", "shift_transfer"):
            sc = [r["score"] for r in res["rows"] if r["model"] == model and r["construction"] == "wc"
                  and r["anchors"] == "computed"]
            if len(sc) >= 5:
                sds.append(float(np.std(sc, ddof=1)))
    out["calibration baseline vs floor"] = sds
    return out


def adopt(diffs: np.ndarray, rule: str) -> np.ndarray:
    """diffs: (sims, tasks, seeds) paired gains. Returns a boolean adoption per simulation."""
    mean = diffs.mean(axis=2)
    k = diffs.shape[1]
    need = 2 if k >= 3 else k
    if rule == "current":
        return (mean > 0).sum(axis=1) >= need
    if rule == "rule7":
        ahead = (diffs > 0).sum(axis=2)
        passed = (mean > MARGIN) & (ahead >= MIN_AHEAD)
        return (passed.sum(axis=1) >= need) & ~(mean < -MARGIN).any(axis=1)
    raise ValueError(rule)


def simulate(gain: float, sigma: float, tau: float, tasks: int, rule: str, sims: int = 100_000, seed: int = 0):
    rng = np.random.default_rng(seed)
    true = gain + tau * rng.standard_normal((sims, tasks, 1))
    diffs = true + sigma * rng.standard_normal((sims, tasks, SEEDS))
    return float(adopt(diffs, rule).mean())


def main() -> None:
    noise = paired_sds()
    pooled = np.array([v for vs in noise.values() for v in vs])
    summary = {k: {"n_pairs": len(v), "median_sd": float(np.median(v)), "q25": float(np.percentile(v, 25)),
                   "q75": float(np.percentile(v, 75))} for k, v in noise.items() if v}
    summary["all"] = {"n_pairs": len(pooled), "median_sd": float(np.median(pooled)),
                      "q25": float(np.percentile(pooled, 25)), "q75": float(np.percentile(pooled, 75))}
    sigmas = sorted({0.2, round(float(np.median(pooled)), 2), 1.0})
    gains = (-0.5, 0.0, 0.3, 0.5, 1.0, 2.0)
    table = []
    for tasks in (3, 2):
        for tau in (0.0, 0.5):
            for sigma in sigmas:
                for rule in ("current", "rule7"):
                    table.append({"tasks": tasks, "tau": tau, "sigma": sigma, "rule": rule,
                                  **{f"gain {g:+.1f}": simulate(g, sigma, tau, tasks, rule) for g in gains}})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gate_analysis.json").write_text(json.dumps({"paired_sd": summary, "table": table}, indent=1))
    lines = ["Paired score-difference sd over seeds (board points), committed results:"]
    lines += [f"  {k}: median {v['median_sd']:.2f} (IQR {v['q25']:.2f}-{v['q75']:.2f}; {v['n_pairs']} pairs)"
              for k, v in summary.items()]
    lines.append("")
    lines.append("P(adopt) by true per-pseudo-task gain; tau = sd of the true gain across pseudo-tasks:")
    head = "tasks tau  sigma rule     " + "  ".join(f"{g:+.1f}" for g in gains)
    lines.append(head)
    for r in table:
        lines.append(f"{r['tasks']:>5} {r['tau']:.1f} {r['sigma']:>6.2f} {r['rule']:<8} "
                     + "  ".join(f"{r[f'gain {g:+.1f}']:.2f}" for g in gains))
    (OUT / "gate_analysis.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
