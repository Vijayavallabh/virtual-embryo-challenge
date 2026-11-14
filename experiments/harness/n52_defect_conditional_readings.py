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
     DEEP and OTHER tasks that enter (none where no DEEP or no OTHER task enters). Beside them the all-draw values;
  3. the same v on the single draw s = 1 of each task (the upload's seed), with its flag, deciding nothing;
  4. each flagged draw of each task with its xcorr_min, thin_ratio and `by` from n51.json and the arm's v there, and per
     task the number of flagged draws with xcorr_min below 0 and at or above 0 (a sign split, not a kind).

Two checks, either of which voids the reading:
  committed    for every arm with a committed reading (COMMITTED_READINGS: rows N45, N46 and N47, their
               first_upload.first_uploads[arm]), the all-draw D, O and S(p) and the per-task V_all equal the committed
               values within CHECK_TOL: the scores read here are the ones those rows priced. An arm without a committed
               reading is named as unchecked;
  arithmetic   the all-draw values equal n37_first_upload_new_arms.price on the same scores within CHECK_TOL. This
               checks the arithmetic only: it reads the same scores, so it cannot tell whether they are the priced ones.

What the flag is not. It pools draws whose coupling partners sit at random (xcorr_min near 0, as on the committed
diagnostic's E7.5 backtest) with draws whose partners are mirrored (xcorr_min well below 0, as on its E7.25 backtest's
flagged draws). On every flagged draw of the committed diagnostic (both backtests and the test pair) both conditions
fire, so `by` reads "both" there and cannot tell the kinds apart. V_flag is a reading of draws past one of two
thresholds, not of draws like the upload's.

  uv run --frozen --no-sync python -I -B experiments/harness/n52_defect_conditional_readings.py \
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
COMMITTED_READINGS = tuple(HERE / "results" / f for f in ("n45_embryo_f_readings/n45.json",
                                                         "n46_embryo_g_readings/n46.json",
                                                         "n47_embryo_h_readings/n47.json"))
N51 = HERE / "results" / "n51_coupling_defect_nine" / "n51.json"
RUNS = (*nt.STORED, HERE / "results" / "e28_nine_tasks" / "tasks", HERE / "results" / "n45_embryo_f_arms" / "tasks",
        HERE / "results" / "n46_embryo_g_arms" / "tasks", HERE / "results" / "n47_embryo_h_arms" / "tasks")
OPTIONAL_RUNS = (HERE / "results" / "n50_embryo_j_arms" / "tasks",)


def run_dirs() -> list[Path]:
    return [*RUNS, *(d for d in OPTIONAL_RUNS if d.is_dir())]


def flags(n51: dict, tasks: dict) -> tuple[dict, dict, list[str]]:
    """Per task name, the flag per seed (a bool array in SEEDS order) and the draws' n51 records in SEEDS order, and
    the reasons the flag table is unusable."""
    out, draws, why = {}, {}, []
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
        rec = {d["seed"]: d for d in r["draws"]}
        draws[n] = [rec[s] for s in SEEDS]
    return out, draws, why


def flag_kinds(flag: dict, draws: dict) -> dict:
    """Per task, its flagged draws' measures from n51.json, the `by` tally and a sign split of xcorr_min."""
    out = {}
    for n, f in flag.items():
        fl = [{"seed": s, "xcorr_min": float(draws[n][i]["xcorr_min"]), "thin_ratio": float(draws[n][i]["thin_ratio"]),
               "by": draws[n][i]["by"]} for i, s in enumerate(SEEDS) if f[i]]
        by: dict = {}
        for d in fl:
            by[d["by"]] = by.get(d["by"], 0) + 1
        out[n] = {"flagged_draws": fl, "by": by,
                  "xcorr_below_0": sum(d["xcorr_min"] < 0 for d in fl),
                  "xcorr_at_or_above_0": sum(d["xcorr_min"] >= 0 for d in fl)}
    return out


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


def read_arm(scores: dict, tasks: dict, flag: dict, arm: str, draws: dict) -> dict:
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
                       "V_flag": _mean(v[f]), "V_rest": _mean(v[~f]), "V_all": float(v.mean()),
                       "flagged_draws": [{"seed": s, "xcorr_min": float(draws[n][i]["xcorr_min"]),
                                          "thin_ratio": float(draws[n][i]["thin_ratio"]), "v": float(v[i])}
                                         for i, s in enumerate(SEEDS) if f[i]]}
        v_flag[n], v_rest[n], v_all[n] = per_task[n]["V_flag"], per_task[n]["V_rest"], per_task[n]["V_all"]
        n_flag[n], n_rest[n] = int(f.sum()), int((~f).sum())
        i = list(SEEDS).index(UPLOAD_SEED)
        upload[n] = {"defect": bool(f[i]), "v": float(v[i])}
    d = float(np.mean([v_all[n] for n in nt.DEEP]))
    o = float(np.mean([v_all[n] for n in nt.OTHER]))
    all_draws = {"D": d, "O": o, **{f"S{p}": p * d + (1 - p) * o for p in nt.PS}}
    return {"tasks": per_task, "flagged": _summary(v_flag, n_flag), "unflagged": _summary(v_rest, n_rest),
            "all_draws": all_draws, "upload_seed": upload}


def _keys() -> tuple[str, ...]:
    return ("D", "O", *(f"S{p}" for p in nt.PS))


def load_committed(paths=None) -> tuple[dict, list[str]]:
    """Per arm, the committed first-upload records (source file name, record) of rows N45-N47 (COMMITTED_READINGS
    unless `paths` is given), and the reasons they cannot be used."""
    out, why = {}, []
    for path in map(Path, COMMITTED_READINGS if paths is None else paths):
        if not path.is_file():
            why.append(f"{path.name}: no committed reading")
            continue
        d = json.loads(path.read_text())
        fu = d.get("first_upload") or {}
        if d.get("void", True) or "first_uploads" not in fu:
            why.append(f"{path.name}: void or without first_upload")
            continue
        if (fu.get("base_first"), fu.get("second")) != (ENTRY, SECOND):
            why.append(f"{path.name}: its pair is ({fu.get('base_first')}, {fu.get('second')}), "
                       f"not ({ENTRY}, {SECOND})")
            continue
        for arm, rec in fu["first_uploads"].items():
            if rec.get("readable", False):
                out.setdefault(arm, []).append((path.name, rec))
    return out, why


def check_committed(arms: dict, committed: dict) -> tuple[float, dict, list[str], list[str]]:
    """The all-draw values and per-task V_all of every arm read against each committed reading of it: the largest
    difference, the sources checked per arm, the arms without a committed reading, and the failures."""
    worst, checked, unchecked, why = 0.0, {}, [], []
    for arm, r in arms.items():
        refs = committed.get(arm, [])
        if not refs:
            unchecked.append(arm)
            continue
        for src, ref in refs:
            pairs = [(r["all_draws"][k], ref[k]) for k in _keys()]
            pairs += [(r["tasks"][n]["V_all"], ref["tasks"][n]) for n in r["tasks"]]
            d = max(abs(a - b) for a, b in pairs)
            worst = max(worst, d)
            checked.setdefault(arm, []).append(src)
            if not d <= CHECK_TOL:
                why.append(f"{arm}: the all-draw values differ from the committed reading {src} by {d:.2e}")
    return worst, checked, unchecked, why


def check_arithmetic(arms: dict, scores: dict, tasks: dict) -> tuple[float, list[str]]:
    """The all-draw values against n37.price on the same scores, an arithmetic check only (it reads the same scores, so
    it cannot tell whether they are the priced ones); the largest difference and the failures."""
    worst, why = 0.0, []
    if not arms:
        return worst, why
    price = n37.price(scores, tasks, candidates=tuple(arms))["first_uploads"]
    for arm, r in arms.items():
        ref = price[arm]
        if not ref.get("readable", False):
            why.append(f"{arm}: n37.price cannot read it")
            continue
        pairs = [(r["all_draws"][k], ref[k]) for k in _keys()]
        pairs += [(r["tasks"][n]["V_all"], ref["tasks"][n]) for n in r["tasks"]]
        d = max(abs(a - b) for a, b in pairs)
        worst = max(worst, d)
        if not d <= CHECK_TOL:
            why.append(f"{arm}: the all-draw values differ from n37.price's on the same scores by {d:.2e}")
    return worst, why


def compute(scores: dict, tasks: dict, n51: dict, committed: dict) -> dict:
    """The reading on per-seed scores (task key -> arm -> array over SEEDS), the flag table and the committed
    first-upload records per arm (load_committed)."""
    names = [*nt.DEEP, *nt.OTHER]
    flag, draws, why = flags(n51, tasks)
    for n in names:
        for arm in (ENTRY, SECOND):
            if arm not in scores.get(tasks[n].key, {}):
                why.append(f"{n}: no scores for {arm}")
    if why:
        return {"void": True, "why": why}
    readable = {a: all(a in scores[tasks[n].key] and np.isfinite(scores[tasks[n].key][a]).all() for n in names)
                for a in ARMS}
    arms = {a: read_arm(scores, tasks, flag, a, draws) for a in ARMS if readable[a]}
    c_worst, checked, unchecked, c_why = check_committed(arms, committed)
    a_worst, a_why = check_arithmetic(arms, scores, tasks)
    checks = {"committed_check_max_abs_difference": c_worst, "committed_checked": checked,
              "no_committed_reading": unchecked, "arithmetic_check_max_abs_difference": a_worst}
    if c_why or a_why:
        return {"void": True, "why": c_why + a_why, **checks}
    return {"void": False, "why": [], "min_draws": MIN_DRAWS, "upload_seed": UPLOAD_SEED,
            "n_flagged": {n: int(flag[n].sum()) for n in names}, "flag_kinds": flag_kinds(flag, draws),
            "not_read": [a for a in ARMS if not readable[a]], **checks, "arms": arms}


def reading(n51_path: Path, dirs=None) -> dict:
    scores, why, meta = nt._load(run_dirs() if dirs is None else dirs)
    if why:
        return {"void": True, "why": why}
    committed, why = load_committed()
    if why:
        return {"void": True, "why": why}
    n51 = json.loads(Path(n51_path).read_text())
    res = compute(scores, nt.all_tasks(), n51, committed)
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
             "all-draw values against the committed readings of rows N45-N47: largest difference "
             f"{res['committed_check_max_abs_difference']:.1e} (void above {CHECK_TOL:.0e}); checked: "
             + (", ".join(f"{a} ({', '.join(s)})" for a, s in res["committed_checked"].items()) or "none")
             + f"; no committed reading, unchecked: {', '.join(res['no_committed_reading']) or 'none'}",
             "arithmetic check (n37.price on the same scores; it cannot tell whether they are the priced ones): "
             f"largest difference {res['arithmetic_check_max_abs_difference']:.1e}",
             "flagged draws per task: " + ", ".join(f"{n} {res['n_flagged'][n]}" for n in names),
             f"not read (no scores on all nine tasks): {', '.join(res['not_read']) or 'none'}",
             "",
             "What the flag is not: it pools draws whose partners sit at random (xcorr_min near 0) with mirrored ones "
             "(xcorr_min well below 0). On every flagged draw of the committed diagnostic (both backtests and the "
             "test pair) both conditions fire, so `by` reads 'both' and cannot tell the kinds apart; V_flag is a "
             "reading of draws past one of two thresholds, not of draws like the upload's.",
             "flagged draws per task (seed xcorr_min thin_ratio); the sign split of xcorr_min is not a kind:"]
    for n in names:
        k = res["flag_kinds"][n]
        lines.append(f"  {n:<6} below 0: {k['xcorr_below_0']}, at or above 0: {k['xcorr_at_or_above_0']}; by: "
                     + (", ".join(f"{b} {c}" for b, c in k["by"].items()) or "-") + "; "
                     + (", ".join(f"s{d['seed']} {d['xcorr_min']:+.2f} {d['thin_ratio']:.2f}"
                                  for d in k["flagged_draws"]) or "none"))
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
        lines.append("  v on each flagged draw: " + "; ".join(
            f"{n} " + " ".join(f"s{d['seed']} {d['v']:+.2f}" for d in r["tasks"][n]["flagged_draws"])
            for n in names if r["tasks"][n]["flagged_draws"]))
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
