#!/usr/bin/env python
"""Tries E25a, E25b and E25c of experiments/harness/TRIES.md (registered before the run): which PAIR of embryo uploads
pays most, now that E24 named a second upload for the heart-series reading, and which geometry goes with it.

Why a pair. On the test pair fgw's moved cloud is flattened in every draw (row N1: thinnest axis 0.21 to 0.25 of its
RMS radius; every released stage lies between 0.37 and 0.53), and try E10's aligned moves give 0.43 there. Today's s2
covers that with E7.5's geometry. If the s2 becomes the heart-series file (E24), the question of what geometry each
upload carries is open again, and the better of two uploads is what counts.

Worlds, as in s2_pricing.py. H: the four heart-series pseudo-tasks (series_world.py). E: the two counted backtests.

Arms (beside the floor), all in one run and seeds:
  fgw                 s1 today
  fgw_shape0          s2 today (E7.5's geometry at test)
  fgw_series          E24a
  fgw_shape0_series   E24c, E24's recommendation for s2
  fgw_aligned         try E11's candidate (not kept there: +2.13 on E7.5, -1.59 on E7.25), read here on four new truths
  fgw_aligned_series  E25a: the series rendering at the aligned moves' positions
  fgw_shape1_series   E25b: the series rendering in the later bracket's geometry
  fgw_late_shape1     E25c: the later bracket's cells alone at their own positions

Definitions. For uploads A and B on a task, P(A, B) = mean over seeds of max(score A, score B). m_H and m_E are the
means of P over the H tasks and the E tasks, and G(p) = p m_H + (1 - p) m_E for the chance p that the test target is a
heart-series sample (a judgement from the organisers' text, never measured).

Rules, fixed before the run.
  (1) The second upload beside fgw. X in (fgw_aligned_series, fgw_shape1_series, fgw_late_shape1) replaces
      fgw_shape0_series only if G(fgw, X) exceeds G(fgw, fgw_shape0_series) by more than 0.3 at p = 0.5, is not below
      it at p = 0.3 or at p = 0.7, and X scores above fgw_shape0_series on at least three of the four H tasks (mean
      paired difference above 0, ahead in at least 8 of 10 seeds). If several qualify, the largest G at p = 0.5.
      The outcome S* is that arm, or fgw_shape0_series.
  (2) The first upload. Of the pairs (fgw_aligned, S), S in (fgw_shape0, fgw_shape0_series, fgw_aligned_series,
      fgw_shape1_series, fgw_late_shape1), take the one with the largest G at p = 0.5. It is put to the captain as a
      change of both uploads only if its G exceeds G(fgw, S*) by more than 0.3 at p = 0.5 and is not below it at
      p = 0.3 or at p = 0.7. Otherwise the recommendation keeps s1 = fgw. Try E11's veto on E7.25 stands either way:
      a change of s1 can only be the captain's ruling, never an adoption.
  (3) Reported; it recommends nothing. The same tables over the defect seeds only: the seeds in which fgw's coupling
      shows row E10's defect on that task's bracket pair (smallest whitened cross-correlation below 0.3, or thinnest
      axis under 0.9 of the thinner bracket's; `flags`, target-free, diag_emb_coupling.py's own measurements at the
      harness's build seeds), as the test pair does in 10 of 10 draws. A task with no defect seed is left out of that
      table and named. If rule (2) does not put an aligned-first pair forward but on the defect seeds its G exceeds
      G(fgw, S*) by more than 0.3 at p = 0.3, 0.5 and 0.7, the first upload is reported to the captain as an open
      question, with both tables, and with no recommendation from this script.
The reading is void if a task, seed, arm or flag file is missing, or results come from different commits or a dirty
tree.

  flags --task NAME --out DIR   the defect flags of ONE task's bracket pair at the ten build seeds (reads the two
                                brackets only)
  run   --task NAME --out DIR   score the arms on ONE task (half, late, new1, new2, e75, e725); start the six
                                together so that they record one commit
  read  --runs DIR --out DIR    apply the rules; write pairs.json and pairs.txt
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s2_pricing as sp  # noqa: E402  (the six tasks and the two worlds, defined once)
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.prep import load_stage  # noqa: E402
from vec.t2_build import interpolate_clouds  # noqa: E402

ARMS = ("floor", "fgw", "fgw_shape0", "fgw_series", "fgw_shape0_series", "fgw_aligned", "fgw_aligned_series",
        "fgw_shape1_series", "fgw_late_shape1")
E24_S2 = "fgw_shape0_series"
NEW = ("fgw_aligned_series", "fgw_shape1_series", "fgw_late_shape1")
S2_MENU = ("fgw_shape0", E24_S2, *NEW)
S1S = ("fgw", "fgw_aligned")
PAIRS = tuple((a, b) for a in S1S for b in S2_MENU)
TODAY = ("fgw", "fgw_shape0")
BAR = sp.BAR
MIN_AHEAD = sp.MIN_AHEAD
PS = sp.PS
H, E = sp.H, sp.E


def _diag():
    spec = importlib.util.spec_from_file_location("diag_emb_coupling", HERE / "diag_emb_coupling.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cmd_flags(a) -> None:
    all_tasks = sp.tasks()
    if a.task not in all_tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(all_tasks)}")
    t = all_tasks[a.task]
    diag = _diag()
    view = harness.TaskView.of(t)
    lo, hi = view.role_days["before"], view.role_days["after"]
    frac = (view.target_day - lo) / (hi - lo)
    genes = list(load_boards()[t.panel_board].genes)
    data = Path(a.data)
    before, after = (load_stage(data / t.inputs[r], genes=genes, drop_duplicates=True) for r in ("before", "after"))
    Cb, Ca, Eb, Ea, lab_b, lab_a = diag._inputs(before, after, genes)
    rows = []
    for s in SEEDS:
        seed = harness.BUILD_SEED_OFFSET + s
        res = interpolate_clouds(Cb, Ca, Eb, Ea, frac, diag.N_DRAW, seed, target_rms=1.0, move_coupling="fgw",
                                 with_coupling=True)
        m = diag._measure(res, frac, lab_b, lab_a, diag._nn_rows(seed, len(Cb), len(Ca)))
        defect = bool(m["xcorr_min"] < diag.DEFECT_XCORR or m["thin_ratio"] < diag.DEFECT_THIN)
        rows.append({"seed": s, "build_seed": seed, "xcorr_min": m["xcorr_min"], "thin_ratio": m["thin_ratio"],
                     "thin_over_rms_at_t": m["shape_at_t"][-1], "thin_over_rms_before": m["shape_before"][-1],
                     "thin_over_rms_after": m["shape_after"][-1], "defect": defect})
        print(f"{t.key} seed {s}: xcorr {m['xcorr_min']:+.2f}, thin ratio {m['thin_ratio']:.2f}, defect {defect}",
              flush=True)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rec = {"task": t.key, "name": a.task, "t": frac, "before": t.inputs["before"], "after": t.inputs["after"],
           "gate": {"DEFECT_XCORR": diag.DEFECT_XCORR, "DEFECT_THIN": diag.DEFECT_THIN},
           "commit": harness._git("rev-parse", "HEAD"),
           "tree_clean": harness._git("status", "--porcelain", "--untracked-files=no") == "", "seeds": rows}
    (out / f"flags_{a.task}.json").write_text(json.dumps(rec, indent=1))


def cmd_run(a) -> None:
    all_tasks = sp.tasks()
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
    res["task_source"] = "experiments/harness/pair_pricing.py"
    res["run_id"] = out.name
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def _scores(r: dict) -> dict[str, np.ndarray]:
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    return {a: np.array([row["scores"][a] for row in rows], float) for a in ARMS}


def _g(per: dict, pair: tuple[str, str], key: str) -> dict | None:
    """G at each p for one pair from the per-task pair means under `key` ("P" or "P_defect"); None if a world has no
    task with a value."""
    k = "|".join(pair)
    h = [per[n][key][k] for n in H if per[n][key].get(k) is not None]
    e = [per[n][key][k] for n in E if per[n][key].get(k) is not None]
    if not h or not e:
        return None
    m_h, m_e = float(np.mean(h)), float(np.mean(e))
    return {"m_H": m_h, "m_E": m_e, **{str(p): p * m_h + (1 - p) * m_e for p in PS}}


def reading(results: dict[str, dict], flags: dict[str, dict]) -> dict:
    all_tasks = sp.tasks()
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
        f = flags.get(name)
        if f is None or f.get("task") != t.key or [x["seed"] for x in f["seeds"]] != list(SEEDS):
            why.append(f"{name}: its defect flags are missing or incomplete")
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    if why:
        return {"void": True, "why": why, "recommendation": "none: the reading is void"}
    per = {}
    for name, t in all_tasks.items():
        s = _scores(results[t.key])
        mask = np.array([x["defect"] for x in flags[name]["seeds"]], bool)
        m = {"mean_score": {a: float(v.mean()) for a, v in s.items()}, "n_defect": int(mask.sum()),
             "vs_fgw": {}, "vs_e24": {}, "P": {}, "P_defect": {}}
        for x in ARMS[2:]:
            d = s[x] - s["fgw"]
            m["vs_fgw"][x] = {"mean": float(d.mean()), "ahead": int((d > 0).sum()),
                              "mean_defect": float(d[mask].mean()) if mask.any() else None,
                              "mean_clear": float(d[~mask].mean()) if (~mask).any() else None}
        for x in NEW:
            d = s[x] - s[E24_S2]
            m["vs_e24"][x] = {"mean": float(d.mean()), "ahead": int((d > 0).sum())}
        for a_, b_ in PAIRS:
            best = np.maximum(s[a_], s[b_])
            m["P"][f"{a_}|{b_}"] = float(best.mean())
            m["P_defect"][f"{a_}|{b_}"] = float(best[mask].mean()) if mask.any() else None
        per[name] = m
    g = {"|".join(p): _g(per, p, "P") for p in PAIRS}
    g_def = {"|".join(p): _g(per, p, "P_defect") for p in PAIRS}
    ref = g[f"fgw|{E24_S2}"]
    # rule (1)
    qualifies = {}
    for x in NEW:
        gx = g[f"fgw|{x}"]
        tasks_ahead = sum(per[n]["vs_e24"][x]["mean"] > 0 and per[n]["vs_e24"][x]["ahead"] >= MIN_AHEAD for n in H)
        qualifies[x] = {"G_gain_at_0.5": gx["0.5"] - ref["0.5"], "not_below_at_0.3": bool(gx["0.3"] >= ref["0.3"]),
                        "not_below_at_0.7": bool(gx["0.7"] >= ref["0.7"]), "H_tasks_ahead": int(tasks_ahead),
                        "passes": bool(gx["0.5"] - ref["0.5"] > BAR and gx["0.3"] >= ref["0.3"]
                                       and gx["0.7"] >= ref["0.7"] and tasks_ahead >= 3)}
    passed = [x for x in NEW if qualifies[x]["passes"]]
    s_star = max(passed, key=lambda x: g[f"fgw|{x}"]["0.5"]) if passed else E24_S2
    base = g[f"fgw|{s_star}"]
    # rule (2)
    best_al = max(S2_MENU, key=lambda x: g[f"fgw_aligned|{x}"]["0.5"])
    ga = g[f"fgw_aligned|{best_al}"]
    request = bool(ga["0.5"] - base["0.5"] > BAR and ga["0.3"] >= base["0.3"] and ga["0.7"] >= base["0.7"])
    # rule (3): reported
    left_out = [n for n in (*H, *E) if per[n]["n_defect"] == 0]
    base_d = g_def[f"fgw|{s_star}"]
    al_d = {x: g_def[f"fgw_aligned|{x}"] for x in S2_MENU}
    open_q = False
    best_al_d = None
    if base_d is not None and all(v is not None for v in al_d.values()):
        best_al_d = max(S2_MENU, key=lambda x: al_d[x]["0.5"])
        open_q = bool(not request and all(al_d[best_al_d][str(p)] - base_d[str(p)] > BAR for p in (0.3, 0.5, 0.7)))
    s2_text = (f"embryo s2 = {s_star}" + (" (E24's recommendation stands)" if s_star == E24_S2
                                           else " (replaces E24's recommendation)"))
    if request:
        s1_text = f"put to the captain as a ruling request: s1 = fgw_aligned with s2 = {best_al}"
    elif open_q:
        s1_text = ("s1 stays fgw by the rule; on the defect seeds the aligned-first pair "
                   f"(fgw_aligned, {best_al_d}) is ahead: an open question for the captain, no recommendation here")
    else:
        s1_text = "s1 stays fgw"
    return {"void": False, "why": [], "tasks": per, "G": g, "G_defect_seeds": g_def, "defect_tasks_left_out": left_out,
            "rule1": qualifies, "s2": s_star, "rule2": {"best_aligned_pair_s2": best_al, "request": request,
                                                       "G_gain_at_0.5": ga["0.5"] - base["0.5"]},
            "rule3": {"best_aligned_pair_s2_defect": best_al_d, "open_question": open_q},
            "recommendation": f"{s2_text}; {s1_text}",
            "commit": next(iter(results.values()))["provenance"]["commit"]}


def _table(g: dict) -> list[str]:
    lines = ["  s1, s2                                  m_H     m_E    " + "  ".join(f"G(p={p})" for p in PS)]
    for k, v in g.items():
        a_, b_ = k.split("|")
        if v is None:
            lines.append(f"  {a_ + ', ' + b_:<38} not available")
            continue
        lines.append(f"  {a_ + ', ' + b_:<38} {v['m_H']:6.2f}  {v['m_E']:6.2f}   "
                     + "   ".join(f"{v[str(p)]:6.2f}" for p in PS))
    return lines


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries E25a-E25c (commit {res['commit']}): {res['recommendation']}",
             "pairs of uploads, the better of the two per seed (all seeds):", *_table(res["G"]),
             "rule (1), against fgw_shape0_series as s2 beside fgw: "
             + "; ".join(f"{x} G(0.5) {q['G_gain_at_0.5']:+.2f}, ahead on {q['H_tasks_ahead']} of 4 H tasks, "
                         f"passes {q['passes']}" for x, q in res["rule1"].items()),
             f"rule (2): best aligned-first pair (fgw_aligned, {res['rule2']['best_aligned_pair_s2']}), G(0.5) "
             f"{res['rule2']['G_gain_at_0.5']:+.2f} against (fgw, {res['s2']}); request {res['rule2']['request']}",
             "rule (3), reported: the same pairs over the defect seeds only"
             + (f" (left out, no defect seed: {', '.join(res['defect_tasks_left_out'])})"
                if res["defect_tasks_left_out"] else "") + ":", *_table(res["G_defect_seeds"]),
             f"  open question for the captain: {res['rule3']['open_question']}"]
    for name, m in res["tasks"].items():
        lines.append(f"{name} ({m['n_defect']} defect seeds of 10): "
                     + ", ".join(f"{a} {v:.2f}" for a, v in m["mean_score"].items()))
        lines.append("  against fgw (all; defect; clear seeds): " + ", ".join(
            f"{x} {v['mean']:+.2f} ({v['ahead']}/10; "
            + (f"{v['mean_defect']:+.2f}" if v["mean_defect"] is not None else "n/a") + "; "
            + (f"{v['mean_clear']:+.2f}" if v["mean_clear"] is not None else "n/a") + ")"
            for x, v in m["vs_fgw"].items()))
        lines.append("  against fgw_shape0_series: " + ", ".join(f"{x} {v['mean']:+.2f} ({v['ahead']}/10)"
                                                                 for x, v in m["vs_e24"].items()))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    keys = {t.key for t in sp.tasks().values()}
    results, flags = {}, {}
    for p in sorted(Path(a.runs).glob("*.json")):
        r = json.loads(p.read_text())
        if p.name.startswith("flags_"):
            flags[r["name"]] = r
        elif isinstance(r, dict) and "rows" in r and r.get("task", {}).get("key") in keys:
            results[r["task"]["key"]] = r
    res = reading(results, flags)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "pairs.json").write_text(json.dumps(res, indent=1))
    (out / "pairs.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("flags", "run"):
        r = sub.add_parser(name)
        r.add_argument("--task", required=True)
        r.add_argument("--out", required=True)
        r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"flags": cmd_flags, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
