#!/usr/bin/env python
"""Tries E28a to E28c of experiments/harness/TRIES.md (registered before the run): the embryo pair table over all nine
practice tasks, the second profile jump included.

Why. An independent check of 2026-11-13 found that the record prices request R-A (the second upload) on the six tasks
of E24 and E25 only, while diagnostic E26 scored the same kind of file on three more tasks (the second profile jump)
and found it mixed: against today's pair the R-A file gives +3.67, 0.00 and -3.09 there. E26 did not score the
top-up in E7.5's geometry (fgw_shape0_hi), E24's first pick (fgw_shape0_series) or the aligned first upload
(fgw_aligned) on those tasks, so no table covers every candidate on every task. This run scores them there.

Tasks. The six of experiments/harness/s2_pricing.py and the three of experiments/harness/second_jump.py:
  DEEP   = half, late, new1, x1, x2, x3   the target shares its DEEPER bracket's profile (the later bracket on the
                                          first three, the earlier one on the last three; diagnostic N23)
  OTHER  = new2 (deep target with a profile of its own), e75 (shallow target with a profile of its own),
           e725 (target like its shallower, earlier bracket)

Arms. First uploads: fgw (today), fgw_aligned. Second uploads: fgw_shape0 (today), fgw_render_hi, fgw_shape0_hi,
fgw_series, fgw_shape0_series, fgw_aligned_series. This run scores all of them on x1 to x3; on the six other tasks
their per-seed scores are taken from the committed runs of E24 and E25 (experiments/harness/results/e24_pricing,
e25_pairs). An arm scored in two of the runs must agree seed for seed to 0.01 of a board point, or the reading is
void (the rendered arms differ by up to 1e-3 between runs: their neighbour search follows the BLAS thread count).

Value of a pair (A, B) on a task: mean over seeds of max(A, B), minus the same for today's pair (fgw, fgw_shape0).
Summaries: D = mean over DEEP; O = mean over OTHER; S(p) = p D + (1 - p) O for the chance p that the hidden target
shares its deeper bracket's profile; also the mean and the worst value over the nine tasks.

Rules, fixed before the run (three arms' values on x1 to x3 were already known from E26: fgw_render_hi, fgw_series
and fgw_aligned_series; the other arms' were not).
  (1) Second upload beside fgw. The file request R-A names (fgw_aligned_series) is replaced in the advice by another
      second upload only if that one's S(0.45) is higher by more than 0.3 and its S(0.3) and S(0.6) are not lower.
      If no second upload has S(0.45) above 0.3, the advice is to keep today's.
  (2) First upload. With S2* the second upload of rule (1), the pair (fgw_aligned, S2*) is put to the captain as a
      request only if its S(0.45) exceeds that of (fgw, S2*) by more than 0.3 and its S(0.3) and S(0.6) are not
      lower. Otherwise both first uploads' values are reported and nothing is requested.
Counted: three tries (fgw_shape0_hi, fgw_shape0_series and fgw_aligned on truths they have not been scored against).

  run   --task x1|x2|x3 --out DIR   one task per process; start the three together
  read  --runs DIR --out DIR        build the nine-task table and apply the rules
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s2_pricing as sp  # noqa: E402
import second_jump as sj  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

FIRST = ("fgw", "fgw_aligned")
SECOND = ("fgw_shape0", "fgw_render_hi", "fgw_shape0_hi", "fgw_series", "fgw_shape0_series", "fgw_aligned_series")
ARMS = ("floor", *FIRST, *SECOND)
RA = "fgw_aligned_series"
DEEP = ("half", "late", "new1", "x1", "x2", "x3")
OTHER = ("new2", "e75", "e725")
PS = (0.3, 0.45, 0.6)
BAR = 0.3
TOL = 0.01  # board points per seed; the rendered arms move by up to 1e-3 between runs (BLAS threads)
STORED = (HERE / "results" / "e24_pricing" / "tasks", HERE / "results" / "e25_pairs" / "tasks")


def all_tasks() -> dict:
    return {**sp.tasks(), **sj.TASKS}


def cmd_run(a) -> None:
    t = sj.TASKS[a.task]
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
    res["task_source"] = "experiments/harness/nine_tasks.py"
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def _load(dirs) -> tuple[dict, list[str], list[str]]:
    """Per task key and arm, the per-seed scores, each arm checked across the runs that hold it."""
    scores: dict[str, dict[str, np.ndarray]] = {}
    why, meta = [], []
    for d in dirs:
        for p in sorted(Path(d).glob("T2e_*.json")):
            r = json.loads(p.read_text())
            key = r["task"]["key"]
            rows = sorted(r["rows"], key=lambda x: x["seed"])
            if [x["seed"] for x in rows] != list(SEEDS):
                why.append(f"{p}: seeds are not {list(SEEDS)}")
                continue
            if not r["provenance"]["tree_clean"]:
                why.append(f"{p}: a dirty tree")
            where = Path(d).parent.name if Path(d).name == "tasks" else Path(d).name
            meta.append(f"{where}:{r['provenance']['commit'][:7]}")
            for arm in rows[0]["scores"]:
                v = np.array([x["scores"][arm] for x in rows], float)
                old = scores.setdefault(key, {}).get(arm)
                if old is not None and float(np.abs(old - v).max()) > TOL:
                    why.append(f"{key}: {arm} differs between runs by {float(np.abs(old - v).max()):.2e}")
                scores[key][arm] = v
    return scores, why, sorted(set(meta))


def reading(runs_dir: Path) -> dict:
    scores, why, meta = _load([*STORED, runs_dir])
    tasks = all_tasks()
    for n, t in tasks.items():
        have = scores.get(t.key, {})
        for arm in (*FIRST, *SECOND):
            if arm not in have:
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": why}
    val = {}
    for a1 in FIRST:
        for a2 in SECOND:
            per = {}
            for n, t in tasks.items():
                s = scores[t.key]
                per[n] = float(np.maximum(s[a1], s[a2]).mean() - np.maximum(s["fgw"], s["fgw_shape0"]).mean())
            d, o = float(np.mean([per[n] for n in DEEP])), float(np.mean([per[n] for n in OTHER]))
            val[f"{a1}|{a2}"] = {"tasks": per, "D": d, "O": o, **{f"S{p}": p * d + (1 - p) * o for p in PS},
                                 "mean9": float(np.mean(list(per.values()))), "worst9": float(min(per.values()))}
    single = {a: {n: float(scores[t.key][a].mean()) for n, t in tasks.items()} for a in (*FIRST, *SECOND)}
    ra = val[f"fgw|{RA}"]
    better = [a2 for a2 in SECOND if a2 != RA and val[f"fgw|{a2}"]["S0.45"] - ra["S0.45"] > BAR
              and val[f"fgw|{a2}"]["S0.3"] >= ra["S0.3"] and val[f"fgw|{a2}"]["S0.6"] >= ra["S0.6"]]
    s2 = max(better, key=lambda a2: val[f"fgw|{a2}"]["S0.45"]) if better else RA
    keep_today = bool(max(val[f"fgw|{a2}"]["S0.45"] for a2 in SECOND) <= BAR)
    if keep_today:
        s2 = "fgw_shape0"
    f0, f1 = val[f"fgw|{s2}"], val[f"fgw_aligned|{s2}"]
    request_first = bool(f1["S0.45"] - f0["S0.45"] > BAR and f1["S0.3"] >= f0["S0.3"] and f1["S0.6"] >= f0["S0.6"])
    return {"void": False, "why": [], "runs": meta, "values": val, "single": single, "rule1_second_upload": s2,
            "rule1_replaces": better, "rule1_keep_today": keep_today, "rule2_request_first_upload": request_first,
            "rule2_gain_S0.45": f1["S0.45"] - f0["S0.45"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*DEEP, *OTHER]
    lines = [f"Tries E28a-E28c (runs {', '.join(res['runs'])}): rule (1) second upload beside fgw: "
             f"{res['rule1_second_upload']}"
             + (" (today's is kept)" if res["rule1_keep_today"] else "")
             + f"; rule (2) aligned first upload requested: {res['rule2_request_first_upload']} "
             f"(S(0.45) {res['rule2_gain_S0.45']:+.2f})",
             "better of two against today's pair (fgw, fgw_shape0), board points:",
             "  first, second                         " + " ".join(f"{n:>6}" for n in names)
             + "  |      D      O  " + " ".join(f"S({p})" for p in PS) + "   mean9  worst9"]
    for k, v in res["values"].items():
        a1, a2 = k.split("|")
        lines.append(f"  {a1 + ', ' + a2:<36} " + " ".join(f"{v['tasks'][n]:+6.2f}" for n in names)
                     + f"  | {v['D']:+6.2f} {v['O']:+6.2f}  " + "  ".join(f"{v[f'S{p}']:+5.2f}" for p in PS)
                     + f"   {v['mean9']:+5.2f}  {v['worst9']:+6.2f}")
    lines.append("single arms, mean board score:")
    lines.append("  arm                    " + " ".join(f"{n:>6}" for n in names))
    for a, v in res["single"].items():
        lines.append(f"  {a:<22} " + " ".join(f"{v[n]:6.2f}" for n in names))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    res = reading(Path(a.runs))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "nine_tasks.json").write_text(json.dumps(res, indent=1))
    (out / "nine_tasks.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", choices=sorted(sj.TASKS), required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
