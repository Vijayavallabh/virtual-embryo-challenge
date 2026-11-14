#!/usr/bin/env python
"""Reading of try E44 of experiments/harness/TRIES.md (registered before any run): `fgw_resid`, the entry `fgw` with
each bracket's drawn cells moved the share L of the way to the straight line through its two neighbouring samples
(the research report's recipe R2 in its form with no fitted number), scored on the nine practice tasks, each with
its trend samples, by `n49_embryo_i_arms.py` on the branch `campaign/emb-i`.

Why. The record read the trend render on the embryo (every drawn cell moved to a trend's value at the target's day)
and never R2's own form, which corrects a bracket at its own day. The seventh independent sweep of the record
(2026-11-14) named that a hard call.

What it reads. Only per-seed board scores committed by counted tries and the run's own build record: the stored runs
of E24, E25 and E28, the committed task files of tries E35 (`results/n36_trend_second_upload/tasks`:
`fgw_trend_render` on the tasks with their trend samples) and the committed task files of try E44 given by `--runs`
(their scores and, under "build_records", which bracket was corrected, its L, the share of panel genes whose factor
was clipped and the mean |log factor|). It opens no data file, runs no candidate and no scorer.

Validity, a check that can fail. The run rescored the floor, `fgw` and `fgw_trend_render`; each must agree with its
committed scores seed for seed within nine_tasks.TOL on every task, or the reading is void. For `fgw` that also
checks that a task extended by its trend samples builds the two-bracket entry.

Reading, fixed in the row before the run:
  first upload   the adoption rule (vec.adoption.keep_change) for `fgw_resid` against `fgw` on the two counted embryo
                 backtests, and row N33's first-upload rule beside the second upload as built. Put forward as a first
                 upload only if it passes both.
  second upload  rows E34 and E35's rule for `fgw_resid` as the second upload beside `fgw`, against the pair as
                 built. Put forward as a second upload only if it passes.
Whatever is put forward goes to an independent recomputation and a cold adversarial check before any request; the
reading requests nothing by itself. If neither passes, the line is closed by this reading.
Reported, deciding nothing: the arm's own mean minus `fgw`'s and minus `fgw_trend_render`'s on every task, and per
task the corrected bracket, its L, and the means over seeds of its clipped share and its mean |log factor|.

  python -I -B experiments/harness/n49_embryo_i_readings.py --runs DIR --out DIR
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
import n38_second_beside_trend as n38  # noqa: E402
import nine_tasks as nt  # noqa: E402

from vec.adoption import keep_change  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ENTRY = "fgw"  # the first upload as built
SECOND = nt.RA  # fgw_aligned_series, the second upload as built
TREND = "fgw_trend_render"  # tries E33 and E35a: every cell to the trend's value at the target's day
ARM = "fgw_resid"  # R2's form: each bracket the share L of the way to its neighbours' line, at its own day
RESCORED = ("floor", ENTRY, TREND)  # scored again by the run; must match the record
BACKTESTS = ("e75", "e725")  # the two counted embryo backtests among the nine tasks
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
N36_RUNS = HERE / "results" / "n36_trend_second_upload" / "tasks"


def adoption(scores: dict, tasks: dict, arm: str, base: str = ENTRY) -> dict:
    """The adoption rule on the two counted backtests, seed for seed (rows are in seed order)."""
    paired = {}
    for n in BACKTESTS:
        s = scores[tasks[n].key]
        paired[n] = {i: (float(b), float(a)) for i, (b, a) in enumerate(zip(s[base], s[arm], strict=True))}
    return keep_change(paired)


def rescored(committed: dict, run: dict, tasks: dict) -> tuple[dict, list[str]]:
    """Per rescored arm, the largest absolute difference between the run's per-seed scores and the committed ones over
    the nine tasks, and the reasons the run is void (an arm missing from the run, or off by more than nt.TOL)."""
    worst, why = {}, []
    for arm in RESCORED:
        big = 0.0
        for n, t in tasks.items():
            new, old = run.get(t.key, {}).get(arm), committed.get(t.key, {}).get(arm)
            if new is None:
                why.append(f"{n}: the run did not rescore {arm}")
                continue
            if old is None:
                why.append(f"{n}: no committed scores for {arm}")
                continue
            d = float(np.abs(np.asarray(new, float) - np.asarray(old, float)).max())
            big = max(big, d)
            if not d <= nt.TOL:
                why.append(f"{n}: {arm} differs from its committed scores by {d:.2e}")
        worst[arm] = big
    return worst, why


def corrections(run_dirs) -> dict:
    """Per task key, what the run recorded of the arm's correction: for each corrected bracket its neighbours and L
    (the first seed's; they follow from the task alone) and the means over the recorded seeds of its clipped share
    and its mean |log factor|. Empty where the task file holds no record."""
    out: dict = {}
    for d in run_dirs:
        for p in sorted(Path(d).glob("T2e_*.json")):
            r = json.loads(p.read_text())
            recs = r.get("build_records", {}).get(ARM, {})
            per_seed = [v.get("resid", {}).get("brackets", {}) for v in recs.values()]
            rec = {}
            for name in sorted({k for b in per_seed for k, v in b.items() if v.get("corrected")}):
                rows = [b[name] for b in per_seed if b.get(name, {}).get("corrected")]
                rec[name] = {"neighbours": rows[0]["neighbours"], "L": float(rows[0]["L"]), "seeds": len(rows),
                             "clipped_share": float(np.mean([x["clipped_share"] for x in rows])),
                             "mean_abs_log_factor": float(np.mean([x["mean_abs_log_factor"] for x in rows]))}
            out[r["task"]["key"]] = rec
    return out


def check(scores: dict, tasks: dict) -> dict:
    names = [*nt.DEEP, *nt.OTHER]
    readable = all(ARM in scores[tasks[n].key] and np.isfinite(scores[tasks[n].key][ARM]).all() for n in names)

    def minus(a: str, b: str) -> dict:
        return {n: float(scores[tasks[n].key][a].mean() - scores[tasks[n].key][b].mean()) for n in names}

    first = n37.price(scores, tasks, candidates=(ARM,))
    second = n38.price(scores, tasks, first=ENTRY, candidates=(ARM,))
    out = {"arm": ARM, "readable": readable,
           "own_minus_fgw": minus(ARM, ENTRY) if readable else None,
           "own_minus_trend": minus(ARM, TREND) if readable else None,
           "adoption": adoption(scores, tasks, ARM) if readable else None,
           "first_upload": first, "second_upload": second}
    out["first_passes"] = bool(readable and out["adoption"]["keep"] and first["first_uploads"][ARM]["passes"])
    out["second_passes"] = bool(readable and second["second_uploads"][ARM]["passes"])
    out["closed"] = not (out["first_passes"] or out["second_passes"])
    return out


def reading(run_dirs) -> dict:
    run_dirs = list(run_dirs)
    committed, why, meta = nt._load([*nt.STORED, E28_RUNS, N36_RUNS])
    run, why_run, meta_run = nt._load(run_dirs)
    why, meta = [*why, *why_run], sorted({*meta, *meta_run})
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (ENTRY, SECOND, TREND):
            if arm not in committed.get(t.key, {}):
                why.append(f"{n}: no committed scores for {arm}")
        if t.key not in run:
            why.append(f"{n}: no task file in the run")
    worst, why_rescored = rescored(committed, run, tasks)
    why = [*why, *why_rescored]
    if why:
        return {"void": True, "why": sorted(set(why)), "rescored_max_abs_difference": worst}
    scores = {t.key: {**committed[t.key], **({ARM: run[t.key][ARM]} if ARM in run[t.key] else {})}
              for t in tasks.values()}
    rec = corrections(run_dirs)
    return {"void": False, "why": [], "runs": meta, "rescored_max_abs_difference": worst,
            "corrections": {n: rec.get(t.key, {}) for n, t in tasks.items()}, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    put = [x for x in ("as a first upload" if res["first_passes"] else "",
                       "as a second upload" if res["second_passes"] else "") if x]
    lines = [f"Reading of try E44 (runs {', '.join(res['runs'])}): fgw with each bracket moved the share L of the way "
             "to its neighbours' line, on the nine practice tasks with their trend samples; aggregation of committed "
             "scores",
             (f"{ARM} is put forward {' and '.join(put)} (to a recomputation and a cold check before any request)"
              if put else f"{ARM} is not put forward: closed by this reading"),
             "  rescored arms against their committed scores, largest difference over tasks and seeds: "
             + ", ".join(f"{a} {res['rescored_max_abs_difference'][a]:.1e}" for a in RESCORED)
             + f" (void above {nt.TOL})"]
    if not res["readable"]:
        return "\n".join([*lines, f"  {ARM}: not readable"]) + "\n"
    a = res["adoption"]
    parts = [f"{n} {t['mean_gain']:+.2f} ({t['seeds_ahead']}/10)" for n, t in a["tasks"].items()]
    v, w = res["first_upload"]["first_uploads"][ARM], res["second_upload"]["second_uploads"][ARM]
    lines.append(f"  adoption against {ENTRY}: " + ", ".join(parts) + f" -> {'kept' if a['keep'] else 'not kept'}")
    for label, x in (("first upload beside " + SECOND + " (row N33's rule)", v),
                     ("second upload beside " + ENTRY + " (rows E34 and E35's rule)", w)):
        r = x["rule"]
        lines.append(f"  {label}: S(0.45) {x['S0.45']:+.2f}, clauses "
                     + " ".join("y" if r[k] else "n" for k in sorted(r))
                     + f" -> {'passes' if x['passes'] else 'does not pass'}")
    lines.append(f"  own mean minus {ENTRY}: " + " ".join(f"{n} {res['own_minus_fgw'][n]:+.2f}" for n in names))
    lines.append(f"  own mean minus {TREND}: " + " ".join(f"{n} {res['own_minus_trend'][n]:+.2f}" for n in names))
    lines.append("  corrected bracket, L, clipped share, mean |log factor|: "
                 + "; ".join(f"{n} " + (", ".join(f"{b} {c['L']:.3f} {c['clipped_share']:.3f} "
                                                  f"{c['mean_abs_log_factor']:.3f}"
                                                  for b, c in res["corrections"][n].items()) or "none")
                             for n in names))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of try E44")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n49.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n49.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
