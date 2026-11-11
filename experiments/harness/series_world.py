#!/usr/bin/env python
"""Tries E23a and E23b of experiments/harness/TRIES.md (registered before the run): the embryo entry in the world where
the target belongs to the later bracket's series of runs.

E22 (experiments/harness/depth_series.py) scored two pseudo-tasks whose target, E8.0, is a heart-series stage between
the shallow E7.5 and a heart-series later bracket. fgw scored 65.99 and 80.53 there, far above its 61.12 and 58.78 on
the counted backtests, and topping the E7.5 cells up to the later bracket's depth added +3.86 and +1.55 (10 of 10 seeds
on each). E22's rule gave no recommendation (on one task 47% of the gain sat on the cell-state group, against the 50%
it asked for; the rest sat on neighbourhoods). Both tasks share one truth, and fgw_render_hi has now been seen on it.

This script asks the two questions E22 leaves open, on four pseudo-tasks defined here (configs/pseudo_tasks.yaml is
unchanged; checked by vec.harness.check_task), all with E7.5 as the earlier bracket and the reference, as at test:

  half   T2e/E7.5+E8.5->E8.0          t = 0.5    (seen in E22)
  late   T2e/E7.5+E8.25_late->E8.0    t = 0.625  (seen in E22)
  new1   T2e/E7.5+E8.5->E8.25_late    t = 0.8    a truth E22 did not read
  new2   T2e/E7.5+E8.75->E8.5         t = 0.8    a truth E22 did not read

Arms: floor, fgw (s1), fgw_shape0 (today's s2), fgw_render_hi and fgw_late (every cell from the later bracket, the
same-series bracket alone; vec.harness._fgw_late), seeds 0-9, bake protocol, build main, the embryo validation
board's panel and cell limits.

Rules, fixed before the run.
  E23b (confirmation of fgw_render_hi on truths it has not been scored against): it is confirmed only if on BOTH new1
       and new2 its paired mean gain over fgw is above 0.3 with at least 8 of 10 seeds ahead.
  E23a (fgw_late): reported on all four tasks against fgw and against fgw_render_hi. It is "the better heart-series
       entry" only if it beats fgw_render_hi by more than 0.3 on average over the four tasks and on at least three of
       them.
  Recommendation to the captain (no upload changes without the captain's ruling; the s2 rule is frozen): if E23b
       confirms, the candidate second upload for the heart-series world is fgw_render_hi, or fgw_late if E23a's
       condition holds; with X that candidate, it is recommended only if D = max(fgw, X) - max(fgw, fgw_shape0) has a
       mean above 0.3 on at least three of the four tasks. Otherwise there is no recommendation.
What each arm costs when the target is NOT in the later bracket's series is read from the counted backtests (E19 for
fgw_render_hi; the harness run of fgw_late registered in the same row). Neither kind of task is adoption evidence here.

The reading is void if a task, seed or arm is missing, the results come from different commits or a dirty tree, or
fgw_shape0's expression-only skills differ from fgw's in any seed.

  run   (repo root) score the arms on the four tasks; one JSON per task under --out
  read  apply the rules to the JSONs in --runs; write reading.json and reading.txt under --out
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
ARMS = ("floor", "fgw", "fgw_shape0", "fgw_render_hi", "fgw_late")
BAR = 0.3
MIN_AHEAD = 8
EXPRESSION_ONLY = ("de_score", "de_direction", "mmd_u")
TOL = 1e-9


def _task(key: str, target: str, after: str, note: str):
    return harness.check_task(harness.PseudoTask(
        key=key, test_board=BOARD, task="T2", kind="interp", target=target, reference="T2/E7.5.h5ad",
        inputs={"before": "T2/E7.5.h5ad", "after": after}, twin_board=PANEL, setting="embryo", adoption=False,
        notes=f"Tries E23a and E23b, {note}."))


TASKS = {
    "half": _task("T2e/E7.5+E8.5->E8.0", "T2/E8.0.h5ad", "T2/E8.5.h5ad", "half: t = 0.5 (seen in E22)"),
    "late": _task("T2e/E7.5+E8.25_late->E8.0", "T2/E8.0.h5ad", "T2/E8.25_late.h5ad", "late: t = 0.625 (seen in E22)"),
    "new1": _task("T2e/E7.5+E8.5->E8.25_late", "T2/E8.25_late.h5ad", "T2/E8.5.h5ad", "new1: t = 0.8"),
    "new2": _task("T2e/E7.5+E8.75->E8.5", "T2/E8.5.h5ad", "T2/E8.75.h5ad", "new2: t = 0.8"),
}
NEW = ("new1", "new2")
SOURCE = "experiments/harness/series_world.py (in-script pseudo-tasks; configs/pseudo_tasks.yaml unchanged)"


def safe(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", key).strip("_")


def cmd_run(a) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cands = [harness.CANDIDATES[n] for n in ARMS]
    for t in TASKS.values():
        missing = [c.name for c in cands if c.name != harness.FLOOR and not c.applies(t)]
        if missing:
            sys.exit(f"{t.key}: arms that do not apply: {missing}")
    for t in TASKS.values():
        path = out / f"{safe(t.key)}.json"
        if path.exists():
            sys.exit(f"{path} exists")
        res = harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data))
        res["task_source"] = SOURCE
        res["run_id"] = out.name
        path.write_text(json.dumps(res, indent=1, default=float))
        print(f"wrote {path}", flush=True)


def _pair(rows, change: str, base: str) -> dict:
    d = np.array([r["scores"][change] - r["scores"][base] for r in rows], float)
    return {"mean": float(d.mean()), "sd": float(d.std(ddof=1)), "ahead": int((d > 0).sum()), "seeds": int(len(d))}


def _per_task(r: dict) -> dict:
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    s = {a: np.array([row["scores"][a] for row in rows], float) for a in ARMS}
    out = {"mean_score": {a: float(v.mean()) for a, v in s.items()},
           "hi_vs_fgw": _pair(rows, "fgw_render_hi", "fgw"), "late_vs_fgw": _pair(rows, "fgw_late", "fgw"),
           "late_vs_hi": _pair(rows, "fgw_late", "fgw_render_hi"), "shape0_vs_fgw": _pair(rows, "fgw_shape0", "fgw"),
           "mean_skills": {a: {m: float(np.mean([row["skills"][a][m] for row in rows]))
                               for m in rows[0]["skills"][a]} for a in ARMS}}
    pair0 = np.maximum(s["fgw"], s["fgw_shape0"])
    for x in ("fgw_render_hi", "fgw_late"):
        d = np.maximum(s["fgw"], s[x]) - pair0
        out[f"D_{x}"] = {"mean": float(d.mean()), "ahead": int((d > 0).sum())}
    return out


def void_reasons(results: dict[str, dict]) -> list[str]:
    why = []
    for name, t in TASKS.items():
        r = results.get(t.key)
        if r is None:
            why.append(f"missing task {name}")
            continue
        if sorted(row["seed"] for row in r["rows"]) != list(SEEDS):
            why.append(f"{name}: seeds are not {list(SEEDS)}")
        for row in r["rows"]:
            gone = [a for a in ARMS if a not in row["scores"]]
            if gone:
                why.append(f"{name} seed {row['seed']}: arms missing {gone}")
                continue
            for m in EXPRESSION_ONLY:
                if abs(row["skills"]["fgw_shape0"][m] - row["skills"]["fgw"][m]) > TOL:
                    why.append(f"{name} seed {row['seed']}: fgw_shape0 {m} differs from fgw's")
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    return why


def reading(results: dict[str, dict]) -> dict:
    why = void_reasons(results)
    if why:
        return {"void": True, "why": why, "recommendation": "none: the reading is void"}
    per = {name: _per_task(results[t.key]) for name, t in TASKS.items()}
    confirmed = all(per[n]["hi_vs_fgw"]["mean"] > BAR and per[n]["hi_vs_fgw"]["ahead"] >= MIN_AHEAD for n in NEW)
    late_means = [per[n]["late_vs_hi"]["mean"] for n in TASKS]
    late_better = bool(np.mean(late_means) > BAR and sum(m > BAR for m in late_means) >= 3)
    cand = "fgw_late" if late_better else "fgw_render_hi"
    d_ok = sum(per[n][f"D_{cand}"]["mean"] > BAR for n in TASKS)
    recommend = bool(confirmed and d_ok >= 3)
    rec = (f"recommend to the captain: embryo s2 = {cand} for the heart-series world" if recommend
           else "no recommendation: s2 stays fgw_shape0")
    return {"void": False, "why": [], "tasks": per, "E23b_confirmed": bool(confirmed),
            "E23a_late_better_than_hi": late_better, "late_vs_hi_mean_over_tasks": float(np.mean(late_means)),
            "candidate": cand, "tasks_with_D_above_bar": int(d_ok), "recommend": recommend, "recommendation": rec,
            "commit": next(iter(results.values()))["provenance"]["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries E23a and E23b (commit {res['commit']}): {res['recommendation']}",
             f"E23b, fgw_render_hi confirmed on the two new truths: {res['E23b_confirmed']}",
             f"E23a, fgw_late better than fgw_render_hi: {res['E23a_late_better_than_hi']} "
             f"(mean over tasks {res['late_vs_hi_mean_over_tasks']:+.2f})",
             f"candidate {res['candidate']}, D above {BAR} on {res['tasks_with_D_above_bar']} of 4 tasks"]
    for name, m in res["tasks"].items():
        sc = ", ".join(f"{a} {v:.2f}" for a, v in m["mean_score"].items())
        lines += [f"{name} {TASKS[name].key}: {sc}"]
        for k in ("hi_vs_fgw", "late_vs_fgw", "late_vs_hi", "shape0_vs_fgw"):
            p = m[k]
            lines.append(f"  {k}: {p['mean']:+.2f} (sd {p['sd']:.2f}, ahead {p['ahead']} of {p['seeds']})")
        lines.append("  D (better of fgw and X, minus better of fgw and shape0): "
                     + ", ".join(f"{x} {m['D_' + x]['mean']:+.2f} (ahead {m['D_' + x]['ahead']})"
                                 for x in ("fgw_render_hi", "fgw_late")))
        lines.append("  skills: " + "; ".join(
            f"{a}: " + ", ".join(f"{k} {v:.3f}" for k, v in m["mean_skills"][a].items())
            for a in ("fgw", "fgw_render_hi", "fgw_late")))
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
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
