#!/usr/bin/env python
"""Diagnostic E26 of experiments/harness/TRIES.md (registered before the run; existing candidates only, so no try is
counted): does the series rendering carry to a second, unrelated profile jump?

Diagnostic N23 found three neighbour steps at which a sample's per-gene detection profile changes. Tries E22 to E25 read
the rendering on one of them only (E7.5 against E8.0's group), the test pair's own. The step from E8.5 to E8.75 is a
second one: E8.0, E8.25_late and E8.5 share a profile and E8.75 has its own. Three practice tasks straddle it with a
target on the E8.0 side:

  x1  T2e/E8.0+E8.75->E8.25_late        t = 0.40
  x2  T2e/E8.0+E8.75->E8.5              t = 0.67
  x3  T2e/E8.25_late+E8.75->E8.5        t = 0.44   (try E16's control task; fgw and the shape arms were scored there)

Here the target shares the EARLIER bracket's profile, which is also the deeper one, so the rendering (the shallower
bracket's cells moved by the extraembryonic factor into the deeper bracket's profile and topped up to its depth) moves
E8.75's cells towards the target's profile: the same operation as on the test pair, on other samples. The published
baseline is a copy of the earlier bracket and shares the target's profile here, so every score sits near 50.

Arms: floor, fgw, fgw_shape0, fgw_render_hi, fgw_series, fgw_aligned_series.

Reading, fixed before the run. The mechanism "carries" if fgw_series beats fgw by more than 0.3 with at least 8 of
10 seeds ahead on at least two of the three tasks, and beats fgw_render_hi (the top-up without the factor) on at least
two. Otherwise the report says the rendering's gain is so far shown on the E7.5 / E8.0 pair only. Either way request
R-A is unchanged in form (the test pair is that pair); the captain is told the outcome with it. No rule, recommendation
or upload depends on this run.

  run   --task x1|x2|x3 --out DIR   one task per process; start the three together
  read  --runs DIR --out DIR        write reading.json and reading.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

ARMS = ("floor", "fgw", "fgw_shape0", "fgw_render_hi", "fgw_series", "fgw_aligned_series")
BAR = 0.3
MIN_AHEAD = 8


def _task(key: str, before: str, after: str, target: str):
    return harness.check_task(harness.PseudoTask(
        key=key, test_board=sw.BOARD, task="T2", kind="interp", target=f"T2/{target}.h5ad",
        reference=f"T2/{before}.h5ad", inputs={"before": f"T2/{before}.h5ad", "after": f"T2/{after}.h5ad"},
        twin_board=sw.PANEL, setting="embryo", adoption=False, notes="Diagnostic E26: the second profile jump."))


TASKS = {"x1": _task("T2e/E8.0+E8.75->E8.25_late", "E8.0", "E8.75", "E8.25_late"),
         "x2": _task("T2e/E8.0+E8.75->E8.5", "E8.0", "E8.75", "E8.5"),
         "x3": _task("T2e/E8.25_late+E8.75->E8.5", "E8.25_late", "E8.75", "E8.5")}


def cmd_run(a) -> None:
    t = TASKS[a.task]
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
    res["task_source"] = "experiments/harness/second_jump.py"
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def reading(results: dict) -> dict:
    why = [f"missing task {n}" for n, t in TASKS.items() if t.key not in results]
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    if why:
        return {"void": True, "why": why}
    per = {}
    for n, t in TASKS.items():
        rows = sorted(results[t.key]["rows"], key=lambda x: x["seed"])
        if [x["seed"] for x in rows] != list(SEEDS):
            return {"void": True, "why": [f"{n}: seeds are not {list(SEEDS)}"]}
        s = {a: np.array([x["scores"][a] for x in rows], float) for a in ARMS}

        def diff(x, y, s=s):
            d = s[x] - s[y]
            return {"mean": float(d.mean()), "ahead": int((d > 0).sum())}

        per[n] = {"mean_score": {a: float(v.mean()) for a, v in s.items()},
                  "series_vs_fgw": diff("fgw_series", "fgw"), "series_vs_hi": diff("fgw_series", "fgw_render_hi"),
                  "hi_vs_fgw": diff("fgw_render_hi", "fgw"), "aligned_series_vs_series":
                      diff("fgw_aligned_series", "fgw_series"), "shape0_vs_fgw": diff("fgw_shape0", "fgw"),
                  "mean_skills": {a: {m: float(np.mean([x["skills"][a][m] for x in rows]))
                                      for m in rows[0]["skills"][a]} for a in ARMS if a != "floor"}}
    over_fgw = sum(per[n]["series_vs_fgw"]["mean"] > BAR and per[n]["series_vs_fgw"]["ahead"] >= MIN_AHEAD for n in per)
    over_hi = sum(per[n]["series_vs_hi"]["mean"] > 0 for n in per)
    carries = bool(over_fgw >= 2 and over_hi >= 2)
    return {"void": False, "why": [], "tasks": per, "tasks_over_fgw": int(over_fgw), "tasks_over_hi": int(over_hi),
            "carries": carries, "commit": prov[0]["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Diagnostic E26 (commit {res['commit']}): the rendering carries to the second profile jump: "
             f"{res['carries']} (over fgw on {res['tasks_over_fgw']} of 3 tasks, over the plain top-up on "
             f"{res['tasks_over_hi']} of 3)"]
    for n, m in res["tasks"].items():
        lines.append(f"{n} {TASKS[n].key}: " + ", ".join(f"{a} {v:.2f}" for a, v in m["mean_score"].items()))
        lines.append("  " + "; ".join(f"{k.replace('_', ' ')} {v['mean']:+.2f} ({v['ahead']}/10)" for k, v in m.items()
                                        if isinstance(v, dict) and "ahead" in v))
        lines.append("  skills fgw / fgw_series: " + ", ".join(
            f"{k} {m['mean_skills']['fgw'][k]:.3f} / {m['mean_skills']['fgw_series'][k]:.3f}"
            for k in m["mean_skills"]["fgw"]))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    keys = {t.key for t in TASKS.values()}
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
    r.add_argument("--task", choices=sorted(TASKS), required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
