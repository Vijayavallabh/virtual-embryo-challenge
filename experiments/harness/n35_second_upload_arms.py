#!/usr/bin/env python
"""Two second-upload arms for the embryo, scored on the nine practice tasks of try E28 (campaign branch
campaign/emb-expr; the experiments/harness/TRIES.md row is registered before any run).

Arms, each a second upload beside the first upload fgw:
  fgw_aligned_hi        a top-up-only file at the aligned positions (row E26, its correction (v)): fgw_aligned's moves
                        with fgw_render_hi's expression (every cell at the deeper bracket's depth, no factor)
  fgw_early_specialist  the earlier-side specialist by gene-wise thinning (experiments/research_2026_11_runs/REPORT2.md,
                        recipe 5): fgw's cells and positions, the deeper bracket's cells thinned gene by gene to the
                        shallower bracket's profile and depth (vec.t2_build.EARLY_RULE)

Tasks: the nine of experiments/harness/nine_tasks.py (DEEP = half, late, new1, x1, x2, x3; OTHER = new2, e75, e725).
Only the two arms are scored here (with the floor, which every harness run scores). fgw and fgw_aligned_series are NOT
rescored: their per-seed scores are taken from the committed runs of E24, E25 and E28 (results/e24_pricing/tasks,
results/e25_pairs/tasks, results/e28_nine_tasks/tasks) through nine_tasks._load, which checks the seeds, clean trees and
agreement between runs to 0.01 per seed. A run directory that holds fgw or fgw_aligned_series makes the reading void.

Value of a second upload X on a task: the mean over seeds of max(fgw, X) minus the same for max(fgw,
fgw_aligned_series), the pair as built since row R4. D, O and S(p) = p D + (1 - p) O as in nine_tasks.py; worst9 = the
lowest task mean of a pair's better of two, in board points; S(0.45) with each task left out in turn.

Rule, fixed in the row before the run. A second upload X is put forward only if all hold:
  (a) S(0.45) > BAR (0.3);
  (b) S(0.3) >= 0 and S(0.6) >= 0;
  (c) the pair (fgw, X)'s worst9 is not lower than the pair (fgw, fgw_aligned_series)'s by more than TOL (0.01);
  (d) S(0.45) stays above 0 whichever single task is left out.
If both pass, the one with the larger S(0.45) is put forward; otherwise nothing is requested and the second upload stays
fgw_aligned_series. Nothing changes an upload without the captain's ruling (docs/S2_RULES.md).

  run   --task NAME --out DIR   score the two arms on ONE task (half, late, new1, new2, e75, e725, x1, x2, x3); start
                                the nine together so that they record one commit (a GPU is not needed)
  read  --runs DIR --out DIR    apply the rule; write n35.json and n35.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nine_tasks as nt  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

BASE_FIRST = "fgw"  # the first upload as built
BASE_SECOND = nt.RA  # fgw_aligned_series, the second upload since row R4
NEW = ("fgw_aligned_hi", "fgw_early_specialist")
ARMS = (harness.FLOOR, *NEW)
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
LEAVE_OUT_P = 0.45


def cmd_run(a) -> None:
    tasks = nt.all_tasks()
    if a.task not in tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(tasks)}")
    t = tasks[a.task]
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
    res["task_source"] = "experiments/harness/n35_second_upload_arms.py"
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def _arms_in(d: Path) -> set[str]:
    """Every arm a run directory's task files hold."""
    arms: set[str] = set()
    for p in sorted(Path(d).glob("T2e_*.json")):
        for row in json.loads(p.read_text())["rows"]:
            arms |= set(row["scores"])
    return arms


def evaluate(scores: dict, tasks: dict) -> dict:
    """The rule on per-seed scores: scores[task key][arm] -> array over SEEDS."""
    names = [*nt.DEEP, *nt.OTHER]

    def pair_mean(second: str, n: str) -> float:
        s = scores[tasks[n].key]
        return float(np.maximum(s[BASE_FIRST], s[second]).mean())

    base_abs = {n: pair_mean(BASE_SECOND, n) for n in names}
    out = {"first": BASE_FIRST, "base_second": BASE_SECOND, "base_pair_task_means": base_abs,
           "base_pair_worst9": min(base_abs.values()), "second_uploads": {}}
    passed = []
    for x in NEW:
        absx = {n: pair_mean(x, n) for n in names}
        per = {n: absx[n] - base_abs[n] for n in names}
        d = float(np.mean([per[n] for n in nt.DEEP]))
        o = float(np.mean([per[n] for n in nt.OTHER]))
        s = {p: p * d + (1 - p) * o for p in nt.PS}
        loo = {}
        for drop in names:
            dd = [per[n] for n in nt.DEEP if n != drop]
            oo = [per[n] for n in nt.OTHER if n != drop]
            loo[drop] = LEAVE_OUT_P * float(np.mean(dd)) + (1 - LEAVE_OUT_P) * float(np.mean(oo))
        worst = min(absx.values())
        rule = {"a_S0.45_above_bar": bool(s[0.45] > nt.BAR),
                "b_S0.3_and_S0.6_not_below_0": bool(s[0.3] >= 0 and s[0.6] >= 0),
                "c_worst9_not_lower": bool(worst >= out["base_pair_worst9"] - nt.TOL),
                "d_leave_one_task_out_positive": bool(min(loo.values()) > 0)}
        ok = all(rule.values())
        if ok:
            passed.append(x)
        out["second_uploads"][x] = {"tasks": per, "D": d, "O": o, **{f"S{p}": v for p, v in s.items()},
                                    "mean9": float(np.mean(list(per.values()))), "worst9_abs": worst,
                                    "S0.45_leave_one_task_out_min": min(loo.values()),
                                    "S0.45_leave_one_task_out": loo, "rule": rule, "passes": ok}
    out["put_forward"] = max(passed, key=lambda x: out["second_uploads"][x]["S0.45"]) if passed else None
    return out


def reading(runs_dir: Path) -> dict:
    why = []
    stored = [*nt.STORED, E28_RUNS]
    for d in stored:
        clash = _arms_in(d) & set(NEW)
        if clash:
            why.append(f"{d}: already holds {sorted(clash)}")
    rescored = _arms_in(runs_dir) & {BASE_FIRST, BASE_SECOND}
    if rescored:
        why.append(f"{runs_dir}: rescored {sorted(rescored)} (they must come from the committed runs only)")
    scores, why_load, meta = nt._load([*stored, runs_dir])
    why += why_load
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (BASE_FIRST, BASE_SECOND, *NEW):
            if arm not in scores.get(t.key, {}):
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": why}
    return {"void": False, "why": [], "runs": meta, **evaluate(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    lines = [f"n35 (runs {', '.join(res['runs'])}): second uploads beside {res['first']}, against the pair as built "
             f"({res['first']}, {res['base_second']})",
             f"put forward: {res['put_forward'] or 'none (the second upload stays ' + res['base_second'] + ')'}",
             "  second upload          " + " ".join(f"{n:>6}" for n in names)
             + "  |      D      O  " + " ".join(f"S({p})" for p in nt.PS) + "   mean9  worst9  LOO-min  rule a b c d"]
    for x, v in res["second_uploads"].items():
        r = v["rule"]
        lines.append(f"  {x:<22} " + " ".join(f"{v['tasks'][n]:+6.2f}" for n in names)
                     + f"  | {v['D']:+6.2f} {v['O']:+6.2f}  " + "  ".join(f"{v[f'S{p}']:+5.2f}" for p in nt.PS)
                     + f"   {v['mean9']:+5.2f}  {v['worst9_abs']:6.2f}  {v['S0.45_leave_one_task_out_min']:+7.2f}  "
                     + " ".join("y" if r[k] else "n" for k in sorted(r)))
    lines.append(f"  the pair as built: worst9 {res['base_pair_worst9']:.2f}; task means "
                 + " ".join(f"{res['base_pair_task_means'][n]:.2f}" for n in names))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} is not empty: refusing to overwrite")
    res = reading(Path(a.runs))
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n35.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n35.txt").write_text(text(res))
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
