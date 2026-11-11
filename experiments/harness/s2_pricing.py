#!/usr/bin/env python
"""Tries E24a, E24b and E24c of experiments/harness/TRIES.md (registered before the run): which second embryo upload
pays most beside fgw, priced over both worlds the released stages can stand in for.

Worlds. H: the target belongs to the later bracket's series of runs (the four pseudo-tasks of
experiments/harness/series_world.py: a heart-series truth, E7.5 as earlier bracket and reference). E: the target
belongs to the earlier bracket's series (the two counted embryo backtests of configs/pseudo_tasks.yaml).

Arms, each scored beside fgw (s1) in the same run and seeds:
  fgw_shape0          today's s2: the earlier bracket's geometry (the cover for fgw's flattened cloud on the test pair)
  fgw_render_hi       E19c: every cell at the deeper bracket's depth, fgw's positions
  fgw_shape0_hi       E24b: both of the above in one file
  fgw_series          E24a: the shallower bracket's cells moved by the extraembryonic per-gene factor, then to the
                      deeper bracket's depth (vec.t2_build.SERIES_RULE), fgw's positions
  fgw_shape0_series   E24c: fgw_series in the earlier bracket's geometry

Rules, fixed before the run.
  (1) fgw_series "adds to the top-up" only if it beats fgw_render_hi by more than 0.3 on average over the four H tasks
      and on at least three of them, ahead in at least 8 of 10 seeds on those. If it does not, the two series arms
      leave the menu.
  (2) For an arm X on a task, V = mean over seeds of max(fgw, X) - fgw: what X adds as a second upload beside fgw.
      m_H(X) and m_E(X) are its means over the H tasks and the E tasks, and G_X(p) = p m_H(X) + (1 - p) m_E(X) for the
      chance p that the test target is a heart-series sample (a judgement from the organisers' text, never measured).
  (3) Recommendation to the captain (no upload changes without the captain's ruling; the s2 rule is frozen): the menu
      arm with the largest G at p = 0.5, provided it beats fgw_shape0's G by more than 0.3 at p = 0.3 and at p = 0.7.
      A combined arm (geometry and depth in one file) is preferred to its fgw-positions twin when their G at p = 0.5
      differ by less than 0.3, because it keeps the cover for the flattened cloud measured on the test pair, which no
      task here can score. Otherwise there is no recommendation.
The reading is void if a task, seed or arm is missing, or results come from different commits or a dirty tree.

  run   --task NAME --out DIR   score the arms on ONE task (half, late, new1, new2, e75, e725); start the six
                                together so that they record one commit
  read  --runs DIR --out DIR    apply the rules; write pricing.json and pricing.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import series_world as sw  # noqa: E402  (the four H tasks, defined once)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

BOARD = "T2:embryo:test_interp"
ARMS = ("floor", "fgw", "fgw_shape0", "fgw_render_hi", "fgw_shape0_hi", "fgw_series", "fgw_shape0_series")
MENU = ARMS[2:]
SERIES = ("fgw_series", "fgw_shape0_series")
TWIN = {"fgw_shape0_hi": "fgw_render_hi", "fgw_shape0_series": "fgw_series"}  # combined arm -> fgw-positions twin
BAR = 0.3
MIN_AHEAD = 8
PS = (0.3, 0.5, 0.7, 0.9)
COUNTED = {"e75": "T2e/E7.25+E8.0->E7.5", "e725": "T2e/E6.75+E7.5->E7.25"}


def tasks() -> dict:
    counted = {t.key: t for t in harness.load_tasks()[BOARD]}
    return {**sw.TASKS, **{name: counted[key] for name, key in COUNTED.items()}}


H, E = tuple(sw.TASKS), tuple(COUNTED)


def cmd_run(a) -> None:
    all_tasks = tasks()
    if a.task not in all_tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(all_tasks)}")
    t = all_tasks[a.task]
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
    res["task_source"] = "experiments/harness/s2_pricing.py"
    res["run_id"] = out.name
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def _per_task(r: dict) -> dict:
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    s = {a: np.array([row["scores"][a] for row in rows], float) for a in ARMS}
    out = {"mean_score": {a: float(v.mean()) for a, v in s.items()}, "V": {}, "vs_fgw": {}}
    for x in MENU:
        d = s[x] - s["fgw"]
        out["vs_fgw"][x] = {"mean": float(d.mean()), "ahead": int((d > 0).sum())}
        out["V"][x] = float(np.maximum(s["fgw"], s[x]).mean() - s["fgw"].mean())
    d = s["fgw_series"] - s["fgw_render_hi"]
    out["series_vs_hi"] = {"mean": float(d.mean()), "ahead": int((d > 0).sum())}
    return out


def reading(results: dict[str, dict]) -> dict:
    all_tasks = tasks()
    why = []
    for name, t in all_tasks.items():
        r = results.get(t.key)
        if r is None:
            why.append(f"missing task {name}")
            continue
        if sorted(row["seed"] for row in r["rows"]) != list(SEEDS):
            why.append(f"{name}: seeds are not {list(SEEDS)}")
        if any(a not in row["scores"] for row in r["rows"] for a in ARMS):
            why.append(f"{name}: an arm is missing")
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    if why:
        return {"void": True, "why": why, "recommendation": "none: the reading is void"}
    per = {name: _per_task(results[t.key]) for name, t in all_tasks.items()}
    sv = [per[n]["series_vs_hi"] for n in H]
    series_adds = bool(np.mean([x["mean"] for x in sv]) > BAR
                       and sum(x["mean"] > BAR and x["ahead"] >= MIN_AHEAD for x in sv) >= 3)
    menu = [x for x in MENU if series_adds or x not in SERIES]
    m_h = {x: float(np.mean([per[n]["V"][x] for n in H])) for x in MENU}
    m_e = {x: float(np.mean([per[n]["V"][x] for n in E])) for x in MENU}
    g = {x: {str(p): p * m_h[x] + (1 - p) * m_e[x] for p in PS} for x in MENU}
    best = max(menu, key=lambda x: g[x]["0.5"])
    twin_of = {v: k for k, v in TWIN.items()}
    if best in twin_of and twin_of[best] in menu and g[best]["0.5"] - g[twin_of[best]]["0.5"] < BAR:
        best = twin_of[best]  # the combined file keeps the cover for the flattened cloud
    robust = all(g[best][str(p)] - g["fgw_shape0"][str(p)] > BAR for p in (0.3, 0.7))
    recommend = bool(best != "fgw_shape0" and robust)
    rec = (f"recommend to the captain: embryo s2 = {best}" if recommend else "no recommendation: s2 stays fgw_shape0")
    return {"void": False, "why": [], "tasks": per, "series_adds_to_the_top_up": series_adds, "menu": menu,
            "m_H": m_h, "m_E": m_e, "G": g, "best_at_half": best, "robust_at_0.3_and_0.7": bool(robust),
            "recommend": recommend, "recommendation": rec,
            "commit": next(iter(results.values()))["provenance"]["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries E24a-E24c (commit {res['commit']}): {res['recommendation']}",
             f"rule (1), fgw_series adds to the top-up: {res['series_adds_to_the_top_up']}",
             f"menu: {', '.join(res['menu'])}",
             "what each arm adds as a second upload beside fgw (mean of max(fgw, X) - fgw):",
             "  arm                 m_H    m_E   " + "  ".join(f"G(p={p})" for p in PS)]
    for x in MENU:
        lines.append(f"  {x:<18} {res['m_H'][x]:+6.2f} {res['m_E'][x]:+6.2f}   "
                     + "   ".join(f"{res['G'][x][str(p)]:+6.2f}" for p in PS))
    for name, m in res["tasks"].items():
        lines.append(f"{name}: " + ", ".join(f"{a} {v:.2f}" for a, v in m["mean_score"].items()))
        lines.append("  against fgw: " + ", ".join(f"{x} {m['vs_fgw'][x]['mean']:+.2f} ({m['vs_fgw'][x]['ahead']}/10)"
                                                    for x in MENU)
                     + f"; series - hi {m['series_vs_hi']['mean']:+.2f} ({m['series_vs_hi']['ahead']}/10)")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    keys = {t.key for t in tasks().values()}
    results = {}
    for p in sorted(Path(a.runs).glob("*.json")):
        r = json.loads(p.read_text())
        if isinstance(r, dict) and "rows" in r and r.get("task", {}).get("key") in keys:
            results[r["task"]["key"]] = r
    res = reading(results)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "pricing.json").write_text(json.dumps(res, indent=1))
    (out / "pricing.txt").write_text(text(res))
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
