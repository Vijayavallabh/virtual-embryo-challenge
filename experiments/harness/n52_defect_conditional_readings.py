#!/usr/bin/env python
"""Reading N52 (a reading of committed files, not a try): the committed per-seed scores of the moving-less arms on the
nine practice tasks, split by row N51's flag on each draw (the entry's coupling of the misshapen kind or not).

Why. On 2026-11-14 a family of arms that move cells less passed the team's keep rules on average over the ten draws.
Per draw the picture is split, and the real upload is one draw, of the flagged kind by the diagnostic of 2026-11-03.
Row N51 flags every draw of the nine tasks with that diagnostic's own measures and thresholds; this script reads the
committed scores by that flag. It decides nothing and puts nothing forward.

What it reads. `n51.json` (committed after its run) and the per-seed board scores committed by counted tries, through
`nine_tasks._load`: the stored runs (`nine_tasks.STORED`), results/e28_nine_tasks/tasks,
results/n45_embryo_f_arms/tasks, results/n46_embryo_g_arms/tasks, results/n47_embryo_h_arms/tasks and, if the folder
exists, results/n50_embryo_j_arms/tasks. It opens no data file, runs no candidate and no scorer.

Arms: fgw_partial, fgw_partial_grp, fgw_damp_m, fgw_damp_m_grp (if present), fgw_damp, fgw_part_rand, fgw_still; the
entry fgw; the second upload fgw_aligned_series. An arm without finite scores on all nine tasks is not read.

For each arm read:
  1. per task: the number of flagged draws; fgw's own mean board score on flagged and on unflagged draws; the arm's
     mean gain over fgw on flagged and on unflagged draws, with the number of draws ahead in each;
  2. the first-upload pair value by flag: per task and draw v = max(arm, second) - max(fgw, second), paired by seed;
     V_flag(task) the mean of v over the task's flagged draws, V_rest(task) over its unflagged draws. A task enters the
     flagged summary only with at least MIN_DRAWS flagged draws, the unflagged summary only with at least MIN_DRAWS
     unflagged draws; the tasks left out are named. D, O and S(p) = p D + (1 - p) O for p in nine_tasks.PS, over the
     DEEP and OTHER tasks that enter (none where no DEEP or no OTHER task enters). Beside them the all-draw values,
     which must equal row N33's (n37_first_upload_new_arms.price) within CHECK_TOL, or the reading is void;
  3. the same v on the single draw s = 1 of each task (the upload's seed), with its flag, deciding nothing.

  python -I -B experiments/harness/n52_defect_conditional_readings.py \
    --n51 experiments/harness/results/n51_coupling_defect_nine/n51.json \
    --out experiments/harness/results/n52_defect_conditional_readings
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import n37_first_upload_new_arms as n37  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ENTRY = "fgw"
SECOND = nt.RA  # fgw_aligned_series, the second upload as built
ARMS = ("fgw_partial", "fgw_partial_grp", "fgw_damp_m", "fgw_damp_m_grp", "fgw_damp", "fgw_part_rand", "fgw_still")
MIN_DRAWS = 3
UPLOAD_SEED = 1  # s = 1, the upload's seed
CHECK_TOL = 1e-9
N51 = HERE / "results" / "n51_coupling_defect_nine" / "n51.json"
RUNS = (*nt.STORED, HERE / "results" / "e28_nine_tasks" / "tasks", HERE / "results" / "n45_embryo_f_arms" / "tasks",
        HERE / "results" / "n46_embryo_g_arms" / "tasks", HERE / "results" / "n47_embryo_h_arms" / "tasks")
OPTIONAL_RUNS = (HERE / "results" / "n50_embryo_j_arms" / "tasks",)


def run_dirs() -> list[Path]:
    return [*RUNS, *(d for d in OPTIONAL_RUNS if d.is_dir())]


def flags(n51: dict, tasks: dict) -> tuple[dict, list[str]]:
    """Per task name, the flag per seed (a bool array in SEEDS order), and the reasons the flag table is unusable."""
    out, why = {}, []
    if not n51.get("clean", False):
        why.append("n51.json comes from a dirty tree")
    for n, t in tasks.items():
        r = n51.get("tasks", {}).get(n)
        if r is None:
            why.append(f"{n}: not in n51.json")
            continue
        if r["key"] != t.key:
            why.append(f"{n}: n51.json's task key is {r['key']}, not {t.key}")
            continue
        by_seed = {d["seed"]: bool(d["defect"]) for d in r["draws"]}
        if sorted(by_seed) != list(SEEDS):
            why.append(f"{n}: n51.json's seeds are not {list(SEEDS)}")
            continue
        out[n] = np.array([by_seed[s] for s in SEEDS], dtype=bool)
    return out, why


def _mean(x: np.ndarray) -> float | None:
    return float(x.mean()) if len(x) else None


def _summary(per: dict, counts: dict) -> dict:
    """D, O and S(p) over the tasks with at least MIN_DRAWS draws; the tasks left out named."""
    enter = [n for n in [*nt.DEEP, *nt.OTHER] if counts[n] >= MIN_DRAWS]
    deep, other = [n for n in nt.DEEP if n in enter], [n for n in nt.OTHER if n in enter]
    d = float(np.mean([per[n] for n in deep])) if deep else None
    o = float(np.mean([per[n] for n in other])) if other else None
    s = {f"S{p}": (p * d + (1 - p) * o if d is not None and o is not None else None) for p in nt.PS}
    return {"tasks": enter, "left_out": [n for n in [*nt.DEEP, *nt.OTHER] if n not in enter], "D": d, "O": o, **s}


def read_arm(scores: dict, tasks: dict, flag: dict, arm: str) -> dict:
    names = [*nt.DEEP, *nt.OTHER]
    per_task, v_flag, v_rest, v_all, n_flag, n_rest, upload = {}, {}, {}, {}, {}, {}, {}
    for n in names:
        s = scores[tasks[n].key]
        f = flag[n]
        gain = s[arm] - s[ENTRY]
        v = np.maximum(s[arm], s[SECOND]) - np.maximum(s[ENTRY], s[SECOND])
        per_task[n] = {"n_flagged": int(f.sum()),
                       "fgw_mean_flagged": _mean(s[ENTRY][f]), "fgw_mean_unflagged": _mean(s[ENTRY][~f]),
                       "gain_flagged": _mean(gain[f]), "ahead_flagged": int((gain[f] > 0).sum()),
                       "gain_unflagged": _mean(gain[~f]), "ahead_unflagged": int((gain[~f] > 0).sum()),
                       "V_flag": _mean(v[f]), "V_rest": _mean(v[~f]), "V_all": float(v.mean())}
        v_flag[n], v_rest[n], v_all[n] = per_task[n]["V_flag"], per_task[n]["V_rest"], per_task[n]["V_all"]
        n_flag[n], n_rest[n] = int(f.sum()), int((~f).sum())
        i = list(SEEDS).index(UPLOAD_SEED)
        upload[n] = {"defect": bool(f[i]), "v": float(v[i])}
    d = float(np.mean([v_all[n] for n in nt.DEEP]))
    o = float(np.mean([v_all[n] for n in nt.OTHER]))
    all_draws = {"D": d, "O": o, **{f"S{p}": p * d + (1 - p) * o for p in nt.PS}}
    return {"tasks": per_task, "flagged": _summary(v_flag, n_flag), "unflagged": _summary(v_rest, n_rest),
            "all_draws": all_draws, "upload_seed": upload}


def check_all_draws(arms: dict, scores: dict, tasks: dict) -> tuple[float, list[str]]:
    """The all-draw values against row N33's pricing (n37.price); the largest difference and the failures."""
    worst, why = 0.0, []
    if not arms:
        return worst, why
    price = n37.price(scores, tasks, candidates=tuple(arms))["first_uploads"]
    for arm, r in arms.items():
        ref = price[arm]
        if not ref.get("readable", False):
            why.append(f"{arm}: n37.price cannot read it")
            continue
        pairs = [(r["all_draws"][k], ref[k]) for k in ("D", "O", *(f"S{p}" for p in nt.PS))]
        pairs += [(r["tasks"][n]["V_all"], ref["tasks"][n]) for n in r["tasks"]]
        d = max(abs(a - b) for a, b in pairs)
        worst = max(worst, d)
        if not d <= CHECK_TOL:
            why.append(f"{arm}: the all-draw values differ from n37.price's by {d:.2e}")
    return worst, why


def compute(scores: dict, tasks: dict, n51: dict) -> dict:
    """The reading on per-seed scores (task key -> arm -> array over SEEDS) and the flag table."""
    names = [*nt.DEEP, *nt.OTHER]
    flag, why = flags(n51, tasks)
    for n in names:
        for arm in (ENTRY, SECOND):
            if arm not in scores.get(tasks[n].key, {}):
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": why}
    readable = {a: all(a in scores[tasks[n].key] and np.isfinite(scores[tasks[n].key][a]).all() for n in names)
                for a in ARMS}
    arms = {a: read_arm(scores, tasks, flag, a) for a in ARMS if readable[a]}
    worst, why = check_all_draws(arms, scores, tasks)
    if why:
        return {"void": True, "why": why, "check_max_abs_difference": worst}
    return {"void": False, "why": [], "min_draws": MIN_DRAWS, "upload_seed": UPLOAD_SEED,
            "n_flagged": {n: int(flag[n].sum()) for n in names},
            "not_read": [a for a in ARMS if not readable[a]], "check_max_abs_difference": worst, "arms": arms}


def reading(n51_path: Path, dirs=None) -> dict:
    scores, why, meta = nt._load(run_dirs() if dirs is None else dirs)
    if why:
        return {"void": True, "why": why}
    n51 = json.loads(Path(n51_path).read_text())
    res = compute(scores, nt.all_tasks(), n51)
    res["runs"] = meta
    res["n51"] = {"path": str(n51_path), "commit": n51.get("commit"), "clean": n51.get("clean")}
    return res


def _f(v, fmt="+6.2f") -> str:
    return "     -" if v is None else format(v, fmt)


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    lines = ["Reading N52: committed per-seed scores split by row N51's flag (the entry's coupling of the misshapen "
             "kind). It decides nothing and puts nothing forward.",
             f"runs {', '.join(res.get('runs', []))}; n51 at {(res.get('n51') or {}).get('commit')}",
             f"all-draw values against n37.price: largest difference {res['check_max_abs_difference']:.1e} "
             f"(void above {CHECK_TOL:.0e})",
             "flagged draws per task: " + ", ".join(f"{n} {res['n_flagged'][n]}" for n in names),
             f"not read (no scores on all nine tasks): {', '.join(res['not_read']) or 'none'}"]
    for arm, r in res["arms"].items():
        lines += ["", f"{arm} (pair value beside {SECOND}, against ({ENTRY}, {SECOND}); board points)",
                  "  task   flagged  fgw flag  fgw rest  gain flag (ahead)  gain rest (ahead)   V_flag  V_rest   V_all"]
        for n in names:
            t = r["tasks"][n]
            lines.append(f"  {n:<6} {t['n_flagged']:>4}/10  {_f(t['fgw_mean_flagged'], '8.2f')}  "
                         f"{_f(t['fgw_mean_unflagged'], '8.2f')}  {_f(t['gain_flagged'])} "
                         f"({t['ahead_flagged']:>2})         {_f(t['gain_unflagged'])} "
                         f"({t['ahead_unflagged']:>2})       {_f(t['V_flag'])}  {_f(t['V_rest'])}  {_f(t['V_all'])}")
        for k in ("flagged", "unflagged", "all_draws"):
            s = r[k]
            row = f"  {k:<10} D {_f(s['D'])}  O {_f(s['O'])}  " + "  ".join(
                f"S({p}) {_f(s[f'S{p}'])}" for p in nt.PS)
            if k != "all_draws":
                row += f"   left out (fewer than {MIN_DRAWS} draws): {', '.join(s['left_out']) or 'none'}"
            lines.append(row)
        lines.append(f"  s = {UPLOAD_SEED} (deciding nothing): " + ", ".join(
            f"{n} {r['upload_seed'][n]['v']:+.2f}{' flagged' if r['upload_seed'][n]['defect'] else ''}"
            for n in names))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n51", default=str(N51))
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if (out / "n52.json").exists() or (out / "n52.txt").exists():
        raise SystemExit(f"{out} already holds n52's outputs: refusing to overwrite")
    res = reading(Path(a.n51))
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n52.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n52.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
