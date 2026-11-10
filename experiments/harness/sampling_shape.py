#!/usr/bin/env python
"""Tries E16 and E17 of experiments/harness/TRIES.md (the captain's go-ahead of 2026-11-11; docs/PLAN.md decision log
2026-11-11 (a)): should the embryo's second upload take its positions from the bracket sampled like the target?

The Data Protocol page, read through the screen (docs/DATA_NOTES.md, 2026-11-10), gives the test target E7.75 the
procedure of E8.0, its later bracket, and docs/PLAN.md 4.2 records E7.5's sampling as unlike E8.0's. Today's s2
(fgw_shape0) takes the earlier bracket's geometry, chosen on a backtest whose target was sampled like its earlier
bracket. Two diagnostic pseudo-tasks, defined here (configs/pseudo_tasks.yaml is unchanged) and checked by
vec.harness.check_task, copy the test's structure on released stages:

  main     T2e/E8.0+E8.5->E8.25_late: the target and the later bracket are heart-tube portions with the yolk sac; the
           earlier bracket, E8.0, runs from the yolk-sac boundary to where the cardiac crescent ends; t = 0.6
  control  T2e/E8.25_late+E8.75->E8.5: all three are heart-tube portions; t = 0.44

Arms: floor, fgw (s1), fgw_shape0 (today's s2, the earlier bracket's geometry) and fgw_shape1 (the later bracket's),
seeds 0-9, the harness's bake protocol and build main, on the embryo validation board's gene panel and cell limits.

Rule, fixed before the run: per seed on the main task, D = max(fgw, fgw_shape1) - max(fgw, fgw_shape0). The switch of
s2 to fgw_shape1 is recommended to the captain if mean D > 0.3 and D > 0 in at least 8 of 10 seeds; otherwise there is
no change. The control decides nothing and is reported. Neither task is counted by the adoption rule, and the frozen s2
rule changes only on the captain's say-so (docs/S2_RULES.md).

The reading is void, and recommends nothing, if a task, seed or arm is missing, the two results come from different
commits or a dirty tree, or the expression-only skills (de_score, de_direction, mmd_u) differ between fgw and either
shape variant in any seed: the three arms share their cells and differ only in positions.

  run   (DGX, repo root) score the arms on both tasks; one JSON per task under --out
  read  apply the rule to the JSONs in --runs; write reading.json and reading.txt under --out
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

from vec import harness
from vec.adoption import SEEDS

BOARD = "T2:embryo:test_interp"
PANEL = "T2:embryo:val_interp"
ARMS = ("floor", "fgw", "fgw_shape0", "fgw_shape1")
BAR = 0.3
MIN_AHEAD = 8
EXPRESSION_ONLY = ("de_score", "de_direction", "mmd_u")
TOL = 1e-9

MAIN = harness.check_task(harness.PseudoTask(
    key="T2e/E8.0+E8.5->E8.25_late", test_board=BOARD, task="T2", kind="interp", target="T2/E8.25_late.h5ad",
    reference="T2/E8.0.h5ad", inputs={"before": "T2/E8.0.h5ad", "after": "T2/E8.5.h5ad"}, twin_board=PANEL,
    setting="embryo", adoption=False,
    notes="Tries E16 and E17, main: the target is sampled like the later bracket and unlike the earlier, t = 0.6."))
CONTROL = harness.check_task(harness.PseudoTask(
    key="T2e/E8.25_late+E8.75->E8.5", test_board=BOARD, task="T2", kind="interp", target="T2/E8.5.h5ad",
    reference="T2/E8.25_late.h5ad", inputs={"before": "T2/E8.25_late.h5ad", "after": "T2/E8.75.h5ad"},
    twin_board=PANEL, setting="embryo", adoption=False,
    notes="Tries E16 and E17, control: all three stages sampled alike, t = 0.44."))
TASKS = (MAIN, CONTROL)
SOURCE = "experiments/harness/sampling_shape.py (in-script pseudo-tasks; configs/pseudo_tasks.yaml unchanged)"


def safe(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", key).strip("_")


def cmd_run(a) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cands = [harness.CANDIDATES[n] for n in ARMS]
    for t in TASKS:
        missing = [c.name for c in cands if c.name != harness.FLOOR and not c.applies(t)]
        if missing:
            sys.exit(f"{t.key}: arms that do not apply: {missing}")
    for t in TASKS:
        path = out / f"{safe(t.key)}.json"
        if path.exists():
            sys.exit(f"{path} exists")
        res = harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data))
        res["task_source"] = SOURCE
        res["run_id"] = out.name
        path.write_text(json.dumps(res, indent=1, default=float))
        print(f"wrote {path}", flush=True)


def _per_task(r: dict) -> dict:
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    s = {a: np.array([row["scores"][a] for row in rows], float) for a in ARMS}
    b0, b1 = np.maximum(s["fgw"], s["fgw_shape0"]), np.maximum(s["fgw"], s["fgw_shape1"])
    d = b1 - b0
    skills = {a: {m: float(np.mean([row["skills"][a][m] for row in rows])) for m in rows[0]["skills"][a]}
              for a in ARMS}
    return {
        "seeds": [row["seed"] for row in rows],
        "mean_score": {a: float(v.mean()) for a, v in s.items()},
        "shape0_minus_fgw": float((s["fgw_shape0"] - s["fgw"]).mean()),
        "shape0_ahead_of_fgw": int((s["fgw_shape0"] > s["fgw"]).sum()),
        "shape1_minus_fgw": float((s["fgw_shape1"] - s["fgw"]).mean()),
        "shape1_ahead_of_fgw": int((s["fgw_shape1"] > s["fgw"]).sum()),
        "pair_with_shape0": float(b0.mean()),
        "pair_with_shape1": float(b1.mean()),
        "D_mean": float(d.mean()),
        "D_ahead": int((d > 0).sum()),
        "D_per_seed": [float(x) for x in d],
        "mean_skills": skills,
    }


def void_reasons(results: dict[str, dict]) -> list[str]:
    """Why the reading may not be read; empty when it may."""
    why = []
    for t in TASKS:
        r = results.get(t.key)
        if r is None:
            why.append(f"missing task {t.key}")
            continue
        seeds = sorted(row["seed"] for row in r["rows"])
        if seeds != list(SEEDS):
            why.append(f"{t.key}: seeds {seeds}, not {list(SEEDS)}")
        for row in r["rows"]:
            gone = [a for a in ARMS if a not in row["scores"]]
            if gone:
                why.append(f"{t.key} seed {row['seed']}: arms missing {gone}")
                continue
            for m in EXPRESSION_ONLY:
                ref = row["skills"]["fgw"][m]
                for a in ("fgw_shape0", "fgw_shape1"):
                    if abs(row["skills"][a][m] - ref) > TOL:
                        why.append(f"{t.key} seed {row['seed']}: {a} {m} differs from fgw's")
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    return why


def reading(results: dict[str, dict]) -> dict:
    why = void_reasons(results)
    if why:
        return {"void": True, "why": why, "passes": False, "recommendation": "none: the reading is void"}
    main, control = _per_task(results[MAIN.key]), _per_task(results[CONTROL.key])
    passes = main["D_mean"] > BAR and main["D_ahead"] >= MIN_AHEAD
    rec = ("recommend to the captain: s2 = fgw_shape1 (the later bracket's geometry, E8.0's at test)" if passes
           else "no change: s2 stays fgw_shape0")
    return {"void": False, "why": [], "rule": {"bar": BAR, "min_ahead": MIN_AHEAD, "task": MAIN.key},
            "main": main, "control": control, "passes": bool(passes), "recommendation": rec,
            "commit": next(iter(results.values()))["provenance"]["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries E16 and E17 (commit {res['commit']}): {res['recommendation']}",
             f"rule: on {res['rule']['task']}, mean D > {BAR} and D > 0 in at least {MIN_AHEAD} of 10 seeds, "
             "D = max(fgw, fgw_shape1) - max(fgw, fgw_shape0)"]
    for name in ("main", "control"):
        m = res[name]
        sc = ", ".join(f"{a} {v:.2f}" for a, v in m["mean_score"].items())
        lines += [f"{name}: {sc}",
                  f"  fgw_shape0 - fgw {m['shape0_minus_fgw']:+.2f} (ahead {m['shape0_ahead_of_fgw']} of 10); "
                  f"fgw_shape1 - fgw {m['shape1_minus_fgw']:+.2f} (ahead {m['shape1_ahead_of_fgw']} of 10)",
                  f"  better of fgw and shape0 {m['pair_with_shape0']:.2f}; of fgw and shape1 "
                  f"{m['pair_with_shape1']:.2f}; D {m['D_mean']:+.2f}, ahead {m['D_ahead']} of 10"]
        shape = {a: {k: round(m["mean_skills"][a][k], 3) for k in ("d2_shape", "occupancy_dice", "neighborhood_mmd")}
                 for a in ("fgw", "fgw_shape0", "fgw_shape1")}
        lines.append(f"  shape skills (d2_shape, occupancy_dice, neighborhood_mmd): {shape}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    runs = Path(a.runs)
    results = {}
    for p in sorted(runs.glob("*.json")):
        r = json.loads(p.read_text())
        if isinstance(r, dict) and "rows" in r and r.get("task", {}).get("key") in {t.key for t in TASKS}:
            results[r["task"]["key"]] = r
    res = reading(results)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
