#!/usr/bin/env python
"""Try E27 of experiments/harness/TRIES.md (registered before the run): the confirmation of diagnostic N24b's pick.

N24 and N24b read the renderings' depth top-up on a target-free round trip (thin a deep released stage, render it
back, compare with held-out cells of the same stage). The top-up as built closes 0.92 to 0.96 of the gap; drawing each
cell's profile from 100 neighbours instead of 15 closes 0.96 to 0.98. That is a gain against thinned cells, so by
N24's reading it needs this confirmation on the practice tasks before it may enter a candidate.

Arms: floor, fgw, fgw_shape0, fgw_aligned_series (request R-A's candidate, the reference here) and
fgw_aligned_series_wide (the same file with vec.t2_build.RENDER_K_WIDE neighbours behind each rendered cell's profile).
Tasks: the six of experiments/harness/s2_pricing.py, one process each, started together.

Rule, fixed before the run. LIKE = (half, late, new1), the tasks whose target shares its later bracket's profile
(diagnostic N23). The wider neighbourhood replaces the one in request R-A's candidate only if
  (a) fgw_aligned_series_wide - fgw_aligned_series is above 0.3 on average over LIKE,
  (b) on at least two of LIKE it is above 0 with at least 8 of 10 seeds ahead, and
  (c) on none of the six tasks the better of fgw and the wide file is more than 0.3 below the better of fgw and
      the reference file.
Otherwise it is not kept, and the round trip's bound stands as the reason. One counted try (embryo 27 to 28).

  run   --task NAME --out DIR   one task (half, late, new1, new2, e75, e725)
  read  --runs DIR --out DIR    apply the rule; write reading.json and reading.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s2_pricing as sp  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

ARMS = ("floor", "fgw", "fgw_shape0", "fgw_aligned_series", "fgw_aligned_series_wide")
REF, NEW = "fgw_aligned_series", "fgw_aligned_series_wide"
LIKE = ("half", "late", "new1")
BAR = 0.3
MIN_AHEAD = 8


def cmd_run(a) -> None:
    all_tasks = sp.tasks()
    if a.task not in all_tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(all_tasks)}")
    t = all_tasks[a.task]
    cands = [harness.CANDIDATES[n] for n in ARMS]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{sw.safe(t.key)}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data))
    res["task_source"] = "experiments/harness/render_wide.py"
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def reading(results: dict) -> dict:
    all_tasks = sp.tasks()
    why = [f"missing task {n}" for n, t in all_tasks.items() if t.key not in results]
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    if why:
        return {"void": True, "why": why}
    per = {}
    for n, t in all_tasks.items():
        rows = sorted(results[t.key]["rows"], key=lambda x: x["seed"])
        if [x["seed"] for x in rows] != list(SEEDS) or any(a not in x["scores"] for x in rows for a in ARMS):
            return {"void": True, "why": [f"{n}: a seed or an arm is missing"]}
        s = {a: np.array([x["scores"][a] for x in rows], float) for a in ARMS}
        d = s[NEW] - s[REF]
        per[n] = {"mean_score": {a: float(v.mean()) for a, v in s.items()},
                  "wide_minus_ref": float(d.mean()), "ahead": int((d > 0).sum()),
                  "pair_ref": float(np.maximum(s["fgw"], s[REF]).mean()),
                  "pair_wide": float(np.maximum(s["fgw"], s[NEW]).mean()),
                  "skills": {a: {m: float(np.mean([x["skills"][a][m] for x in rows])) for m in rows[0]["skills"][a]}
                             for a in (REF, NEW)}}
    a_ok = bool(np.mean([per[n]["wide_minus_ref"] for n in LIKE]) > BAR)
    b_ok = bool(sum(per[n]["wide_minus_ref"] > 0 and per[n]["ahead"] >= MIN_AHEAD for n in LIKE) >= 2)
    c_ok = bool(all(per[n]["pair_wide"] >= per[n]["pair_ref"] - BAR for n in per))
    mean_like = float(np.mean([per[n]["wide_minus_ref"] for n in LIKE]))
    return {"void": False, "why": [], "tasks": per, "mean_over_like": mean_like, "a": a_ok, "b": b_ok, "c": c_ok,
            "kept": bool(a_ok and b_ok and c_ok), "commit": prov[0]["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Try E27 (commit {res['commit']}): the wider neighbourhood is kept: {res['kept']} (mean over half, "
             f"late, new1 {res['mean_over_like']:+.2f} against a bar of {BAR}; clauses a {res['a']}, b {res['b']}, "
             f"c {res['c']})"]
    for n, m in res["tasks"].items():
        lines.append(f"{n}: " + ", ".join(f"{a} {v:.2f}" for a, v in m["mean_score"].items())
                     + f"; wide - reference {m['wide_minus_ref']:+.2f} ({m['ahead']}/10); better of two with fgw "
                     f"{m['pair_wide']:.2f} against {m['pair_ref']:.2f}")
        lines.append("  skills reference / wide: " + ", ".join(
            f"{k} {m['skills'][REF][k]:.3f} / {m['skills'][NEW][k]:.3f}" for k in ("mmd_u", "variogram",
                                                                                 "neighborhood_mmd", "de_score",
                                                                                 "de_direction")))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    keys = {t.key for t in sp.tasks().values()}
    results = {}
    for p in sorted(Path(a.runs).glob("*.json")):
        r = json.loads(p.read_text())
        if isinstance(r, dict) and "rows" in r and r.get("task", {}).get("key") in keys:
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
    r.add_argument("--task", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
