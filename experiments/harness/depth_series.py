#!/usr/bin/env python
"""Tries E22a and E22b of experiments/harness/TRIES.md (registered before the run): what the embryo's depth rendering
pays when the target is as deep as its later bracket and the earlier bracket is shallow, the test's likely situation.

Why. Inside every coarse lineage group, E8.0's cells hold about three times E7.5's counts (a target-free reading of the
test's own inputs, 2026-11-12), so the depth step between the two is a property of the runs, not of composition. The
organisers' Data Protocol text (docs/DATA_NOTES.md, released through the screen on 2026-11-10) places the test target
E7.75 in the MERFISH heart series with E8.0 ("sections from E7.75 to E9.5 were chosen to contain heart tissue"), whose
released stages are deep (212, 233, 166, 130 counts per cell), while E7.5 belongs to the shallow earlier series (30,
35, 62). If E7.75 is deep, half of the entry's mixture (its E7.5 cells) is about 3.4 times too shallow, and tries
E19a-E19c measured too shallow as the costly side (-7.07 on E7.5 for a prediction 1.8 times too shallow, against -0.95
for one 3.4 times too deep). No counted backtest has a deep target with a shallow earlier bracket. Two diagnostic
pseudo-tasks, defined here (configs/pseudo_tasks.yaml is unchanged) and checked by vec.harness.check_task, do:

  half  T2e/E7.5+E8.5->E8.0        earlier bracket and reference E7.5 (62), later bracket E8.5 (166), target E8.0
                                   (212); t = 0.5, the test's fraction
  late  T2e/E7.5+E8.25_late->E8.0  the same with E8.25_late (233) as the later bracket; t = 0.625

Arms: floor, fgw (s1), fgw_shape0 (today's s2), fgw_render_hi (E22a: every cell at the deeper bracket's depth, the
shallower bracket topped up) and fgw_render_mid (E22b: every cell at the log-linear depth), seeds 0-9, the harness's
bake protocol and build main, on the embryo validation board's gene panel and cell limits.

Rule, fixed before the run. Per seed and task, H = fgw_render_hi - fgw and D = max(fgw, fgw_render_hi) - max(fgw,
fgw_shape0). A switch of the embryo's s2 to fgw_render_hi is recommended to the captain only if, on BOTH tasks, mean
H > 1.0, mean D > 0.3 and D > 0 in at least 8 of 10 seeds, and at least half of H sits on the cell-state group (mmd_u
and variogram). Otherwise there is no recommendation. fgw_render_mid decides nothing: it is the dose check (if depth is
the cause, its cell-state skill lies between fgw's and fgw_render_hi's). Neither task is counted by the adoption rule,
and the frozen s2 rule changes only on the captain's say-so (docs/S2_RULES.md).

The reading is void if a task, seed or arm is missing, the results come from different commits or a dirty tree, or
fgw_shape0's expression-only skills differ from fgw's in any seed (they share their cells).

  run   (repo root) score the arms on both tasks; one JSON per task under --out
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
from vec.boards import GROUPS

BOARD = "T2:embryo:test_interp"
PANEL = "T2:embryo:val_interp"
ARMS = ("floor", "fgw", "fgw_shape0", "fgw_render_hi", "fgw_render_mid")
MIN_H = 1.0
BAR = 0.3
MIN_AHEAD = 8
MIN_STATE_SHARE = 0.5
EXPRESSION_ONLY = ("de_score", "de_direction", "mmd_u")
TOL = 1e-9

HALF = harness.check_task(harness.PseudoTask(
    key="T2e/E7.5+E8.5->E8.0", test_board=BOARD, task="T2", kind="interp", target="T2/E8.0.h5ad",
    reference="T2/E7.5.h5ad", inputs={"before": "T2/E7.5.h5ad", "after": "T2/E8.5.h5ad"}, twin_board=PANEL,
    setting="embryo", adoption=False,
    notes="Tries E22a and E22b, half: a deep target between a shallow earlier bracket and a deep later one, t = 0.5."))
LATE = harness.check_task(harness.PseudoTask(
    key="T2e/E7.5+E8.25_late->E8.0", test_board=BOARD, task="T2", kind="interp", target="T2/E8.0.h5ad",
    reference="T2/E7.5.h5ad", inputs={"before": "T2/E7.5.h5ad", "after": "T2/E8.25_late.h5ad"}, twin_board=PANEL,
    setting="embryo", adoption=False,
    notes="Tries E22a and E22b, late: the same with E8.25_late as the later bracket, t = 0.625."))
TASKS = (HALF, LATE)
SOURCE = "experiments/harness/depth_series.py (in-script pseudo-tasks; configs/pseudo_tasks.yaml unchanged)"


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


def _state_points(row: dict, arm: str) -> float:
    """The cell-state group's board points of one arm in one seed: 100 x group weight x its group skill."""
    weight = next(w for name, w, _ in GROUPS["T2"] if name == "cell_state")
    return 100.0 * weight * row["groups"][arm]["cell_state"]


def _per_task(r: dict) -> dict:
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    s = {a: np.array([row["scores"][a] for row in rows], float) for a in ARMS}
    h = s["fgw_render_hi"] - s["fgw"]
    d = np.maximum(s["fgw"], s["fgw_render_hi"]) - np.maximum(s["fgw"], s["fgw_shape0"])
    state = np.array([_state_points(row, "fgw_render_hi") - _state_points(row, "fgw") for row in rows])
    skills = {a: {m: float(np.mean([row["skills"][a][m] for row in rows])) for m in rows[0]["skills"][a]}
              for a in ARMS}
    cs = {a: 0.6 * skills[a]["mmd_u"] + 0.4 * skills[a]["variogram"]
          for a in ("fgw", "fgw_render_mid", "fgw_render_hi")}
    return {
        "seeds": [row["seed"] for row in rows],
        "mean_score": {a: float(v.mean()) for a, v in s.items()},
        "H_mean": float(h.mean()), "H_sd": float(h.std(ddof=1)), "H_ahead": int((h > 0).sum()),
        "H_on_cell_state": float(state.mean()),
        "H_state_share": float(state.mean() / h.mean()) if h.mean() > 0 else None,
        "mid_minus_fgw": float((s["fgw_render_mid"] - s["fgw"]).mean()),
        "shape0_minus_fgw": float((s["fgw_shape0"] - s["fgw"]).mean()),
        "pair_with_hi": float(np.maximum(s["fgw"], s["fgw_render_hi"]).mean()),
        "pair_with_shape0": float(np.maximum(s["fgw"], s["fgw_shape0"]).mean()),
        "D_mean": float(d.mean()), "D_ahead": int((d > 0).sum()), "D_per_seed": [float(x) for x in d],
        "cell_state_skill": cs,
        "dose_ordered": bool(min(cs["fgw"], cs["fgw_render_hi"]) <= cs["fgw_render_mid"]
                             <= max(cs["fgw"], cs["fgw_render_hi"])),
        "mean_skills": skills,
    }


def void_reasons(results: dict[str, dict]) -> list[str]:
    why = []
    for t in TASKS:
        r = results.get(t.key)
        if r is None:
            why.append(f"missing task {t.key}")
            continue
        if sorted(row["seed"] for row in r["rows"]) != list(SEEDS):
            why.append(f"{t.key}: seeds are not {list(SEEDS)}")
        for row in r["rows"]:
            gone = [a for a in ARMS if a not in row["scores"]]
            if gone:
                why.append(f"{t.key} seed {row['seed']}: arms missing {gone}")
                continue
            for m in EXPRESSION_ONLY:
                if abs(row["skills"]["fgw_shape0"][m] - row["skills"]["fgw"][m]) > TOL:
                    why.append(f"{t.key} seed {row['seed']}: fgw_shape0 {m} differs from fgw's")
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    return why


def _passes(m: dict) -> bool:
    share = m["H_state_share"]
    return bool(m["H_mean"] > MIN_H and m["D_mean"] > BAR and m["D_ahead"] >= MIN_AHEAD
                and share is not None and share >= MIN_STATE_SHARE)


def reading(results: dict[str, dict]) -> dict:
    why = void_reasons(results)
    if why:
        return {"void": True, "why": why, "passes": False, "recommendation": "none: the reading is void"}
    per = {t.key: _per_task(results[t.key]) for t in TASKS}
    ok = {k: _passes(m) for k, m in per.items()}
    passes = all(ok.values())
    rec = ("recommend to the captain: embryo s2 = fgw_render_hi (every cell at the deeper bracket's depth)" if passes
           else "no recommendation: s2 stays fgw_shape0")
    return {"void": False, "why": [], "rule": {"min_H": MIN_H, "bar_D": BAR, "min_ahead": MIN_AHEAD,
                                               "min_state_share": MIN_STATE_SHARE}, "tasks": per, "task_passes": ok,
            "passes": bool(passes), "recommendation": rec,
            "commit": next(iter(results.values()))["provenance"]["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries E22a and E22b (commit {res['commit']}): {res['recommendation']}",
             f"rule: on both tasks, mean H > {MIN_H}, mean D > {BAR}, D > 0 in at least {MIN_AHEAD} of 10 seeds and at "
             f"least {MIN_STATE_SHARE:.0%} of H on cell_state; H = fgw_render_hi - fgw, "
             "D = max(fgw, fgw_render_hi) - max(fgw, fgw_shape0)"]
    for key, m in res["tasks"].items():
        sc = ", ".join(f"{a} {v:.2f}" for a, v in m["mean_score"].items())
        share = "n/a" if m["H_state_share"] is None else f"{m['H_state_share']:.2f}"
        lines += [f"{key} ({'passes' if res['task_passes'][key] else 'does not pass'}): {sc}",
                  f"  H {m['H_mean']:+.2f} (sd {m['H_sd']:.2f}, ahead {m['H_ahead']} of 10), on cell_state "
                  f"{m['H_on_cell_state']:+.2f} (share {share}); mid - fgw {m['mid_minus_fgw']:+.2f}; "
                  f"shape0 - fgw {m['shape0_minus_fgw']:+.2f}",
                  f"  better of fgw and hi {m['pair_with_hi']:.2f}; of fgw and shape0 {m['pair_with_shape0']:.2f}; "
                  f"D {m['D_mean']:+.2f}, ahead {m['D_ahead']} of 10",
                  "  cell-state skill (0.6 mmd_u + 0.4 variogram): "
                  + ", ".join(f"{a} {v:.3f}" for a, v in m["cell_state_skill"].items())
                  + f"; dose ordered: {m['dose_ordered']}",
                  "  skills: " + "; ".join(
                      f"{a}: " + ", ".join(f"{k} {v:.3f}" for k, v in m["mean_skills"][a].items())
                      for a in ("fgw", "fgw_render_hi", "fgw_render_mid", "fgw_shape0"))]
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    results = {}
    for p in sorted(Path(a.runs).glob("*.json")):
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
