#!/usr/bin/env python
"""Readings of tries E42a to E42d of experiments/harness/TRIES.md (registered before any run): four controls for
`fgw_partial` (try E41b), scored on the nine practice tasks by `n46_embryo_g_arms.py` on the branch
`campaign/emb-g`. Each is the entry `fgw` with one change to its anchors' displacements, on the entry's own coupling,
partners and rotation; none of them uses which anchors the partial coupling keeps:
  `fgw_still`      no anchor moves;
  `fgw_damp`       every anchor moves by the share s of its displacement, s the share `fgw_partial` transports;
  `fgw_part_rand`  the share s of each bracket's anchors, chosen at random, moves as in the entry; the rest stay;
  `fgw_damp_m`     every anchor's displacement is multiplied by one factor per bracket, set so that the bracket's
                   cells move, on average, exactly as far as under `fgw_partial` (the movement-matched control).

Why. `fgw_partial` passed the adoption rule and the first-upload pair rule, but 27 to 90% of its anchors do not
move, so its row could not say whether the gain comes from partial matching or only from moving the cloud less.
The captain asked for "the same test with cells simply moved less". The first three controls match the share of
anchors that moves and not the distance: their builder and their reviewer showed that they move the cloud about 1.3
to 1.8 times as far as `fgw_partial`, so `fgw_damp_m` was added before any run.

What it reads. Only per-seed board scores committed by counted tries and the run's own movement record: the stored
runs of E24, E25 and E28 (`nine_tasks.STORED` and `results/e28_nine_tasks/tasks`), the committed task files of tries
E41 (`results/n45_embryo_f_arms/tasks`: `fgw_partial`) and the committed task files of tries E42 given by `--runs`
(their scores and, under "build_records", each arm's mean displacement norm as a ratio to the entry's, per bracket,
at the anchors and at the cells). It opens no data file, runs no candidate and no scorer.

Validity, checks that can fail. (1) The run of tries E42 rescored `fgw`, `fgw_partial` and the floor: each must
agree with its committed scores seed for seed within nine_tasks.TOL on every task. (2) Every arm the run scored
holds its movement record for every seed. (3) The matched control is matched: on every task, seed and bracket its
cell-level ratio equals `fgw_partial`'s within MATCH_TOL. Otherwise the reading is void.

Readings, all fixed in the row before the run:
  rules     for each control and for `fgw_partial`: the adoption rule (vec.adoption.keep_change) against `fgw` on the
            two counted embryo backtests (the nine tasks' e75 and e725), and row N33's first-upload rule
            (n37_first_upload_new_arms.price) beside the second upload as built. An arm "passes" if it passes both.
  lead      `fgw_partial` against each control, paired by seed, by the adoption rule's own test on the two counted
            backtests (a paired mean above 0.3 with at least 8 of 10 seeds ahead on both, no loss beyond 0.3).
  verdict   "matching" if `fgw_partial` leads every one of the four controls by that test: no way of simply moving
            less reaches it. Otherwise "moving less": a control that uses no matching is not led, so the gain is
            not shown to need the partial coupling's choices. "Not led" is not "equal": the test has ten seeds.
  the arm   put forward: `fgw_partial` if the verdict is "matching". Otherwise the first arm that passes, in the
            order `fgw_still`, `fgw_damp`, `fgw_part_rand`, `fgw_damp_m`, `fgw_partial`: the arm with the fewest
            unrecorded choices first (`fgw_still` reads no label and uses no number; `fgw_damp_m` needs the partial
            solver for its two factors). None if no arm passes.
Any arm put forward goes to a cold adversarial check before a request; the reading requests nothing by itself.
Reported, deciding nothing: each arm's own mean minus `fgw`'s on every task, `fgw_partial`'s own mean minus each
control's on every task, the second-upload rule for the controls, and how far each arm moves its cells.

What the controls do not separate: `fgw_partial` also turns the later bracket by its own rotation and gives its
matched anchors other partners than the entry's. A verdict "matching" therefore means that the partial coupling's
choices as a whole are needed, not which of them.

  python -I -B experiments/harness/n46_embryo_g_readings.py --runs DIR --out DIR
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

from vec.adoption import SEEDS, keep_change  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ENTRY = "fgw"  # the first upload as built
SECOND = nt.RA  # fgw_aligned_series, the second upload as built
PARTIAL = "fgw_partial"  # try E41b, the arm under control
STILL, DAMP, RAND, MATCHED = "fgw_still", "fgw_damp", "fgw_part_rand", "fgw_damp_m"
CONTROLS = (STILL, DAMP, RAND, MATCHED)  # none uses which anchors the partial coupling keeps
ORDER = (*CONTROLS, PARTIAL)  # fewest unrecorded choices first
MOVED = (PARTIAL, *CONTROLS)  # the arms whose movement record the run holds
RESCORED = ("floor", ENTRY, PARTIAL)  # scored again by the run; must match the record
BACKTESTS = ("e75", "e725")  # the two counted embryo backtests among the nine tasks
MATCH_TOL = 1e-6  # the matched control's cell-level ratio against fgw_partial's (equal by construction)
MOVEMENT_KEY = "movement_ratio"  # build_records[arm][seed][MOVEMENT_KEY] = {"anchor": [b, a], "cell": [b, a]}
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
N45_RUNS = HERE / "results" / "n45_embryo_f_arms" / "tasks"


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


def movement(run_dirs) -> dict:
    """Per task key, arm and seed, the run's movement record: {"anchor": [before, after], "cell": [before, after]}, the
    arm's mean displacement norm as a ratio to the entry's on the same anchors and the same cells. None where the task
    file holds none."""
    out: dict = {}
    for d in run_dirs:
        for p in sorted(Path(d).glob("T2e_*.json")):
            r = json.loads(p.read_text())
            recs = r.get("build_records", {})
            out[r["task"]["key"]] = {a: {int(s): v.get(MOVEMENT_KEY) for s, v in recs.get(a, {}).items()}
                                     for a in MOVED}
    return out


def _pair(rec, level: str):
    """The record's two ratios at one level (before, after) as floats, or None if it does not hold two finite ones."""
    v = rec.get(level) if isinstance(rec, dict) else None
    if not isinstance(v, (list, tuple)) or len(v) != 2 or any(x is None for x in v):
        return None
    v = [float(x) for x in v]
    return v if all(np.isfinite(v)) else None


def moved(moves: dict, run: dict, tasks: dict) -> tuple[dict, float | None, list[str]]:
    """How far each arm moves, and the check that the matched control is matched.

    Returns (per arm and task, the mean over seeds of the cell-level and the anchor-level ratio per bracket; the
    largest absolute difference between the matched control's cell-level ratio and fgw_partial's over tasks, seeds
    and brackets; the reasons the run is void: an arm the run scored without its record for a seed, or a matched
    control that is off fgw_partial's cell-level ratio by more than MATCH_TOL)."""
    table: dict = {a: {} for a in MOVED}
    worst, why = None, []
    for n, t in tasks.items():
        have = moves.get(t.key, {})
        per: dict = {}
        for arm in MOVED:
            if arm not in run.get(t.key, {}):
                continue  # not scored on this task (not buildable): nothing to hold
            recs = {s: have.get(arm, {}).get(s) for s in SEEDS}
            pairs = {s: {lv: _pair(recs[s], lv) for lv in ("cell", "anchor")} for s in SEEDS}
            bad = [s for s in SEEDS if pairs[s]["cell"] is None or pairs[s]["anchor"] is None]
            if bad:
                why.append(f"{n}: no movement record for {arm} (seeds {bad})")
                continue
            per[arm] = pairs
            table[arm][n] = {lv: [float(np.mean([pairs[s][lv][b] for s in SEEDS])) for b in (0, 1)]
                             for lv in ("cell", "anchor")}
        if MATCHED in per and PARTIAL in per:
            d = max(abs(per[MATCHED][s]["cell"][b] - per[PARTIAL][s]["cell"][b]) for s in SEEDS for b in (0, 1))
            worst = d if worst is None else max(worst, d)
            if not d <= MATCH_TOL:
                why.append(f"{n}: {MATCHED} moves its cells {d:.2e} off {PARTIAL} (ratio to the entry's)")
    return table, worst, why


def check(scores: dict, tasks: dict) -> dict:
    names = [*nt.DEEP, *nt.OTHER]
    readable = {a: all(a in scores[tasks[n].key] and np.isfinite(scores[tasks[n].key][a]).all() for n in names)
                for a in ORDER}

    def minus(a: str, b: str) -> dict:
        return {n: float(scores[tasks[n].key][a].mean() - scores[tasks[n].key][b].mean()) for n in names}

    arms = [a for a in ORDER if readable[a]]
    first = n37.price(scores, tasks, candidates=ORDER)
    out = {"arms": list(ORDER), "readable": readable,
           "own_minus_fgw": {a: minus(a, ENTRY) for a in arms},
           "partial_minus_control": {c: minus(PARTIAL, c) for c in CONTROLS if readable[c] and readable[PARTIAL]},
           "adoption": {a: (adoption(scores, tasks, a) if readable[a] else None) for a in ORDER},
           "first_upload": first,
           "second_upload": n38.price(scores, tasks, first=ENTRY, candidates=CONTROLS),
           "lead": {c: (adoption(scores, tasks, PARTIAL, base=c) if readable[c] and readable[PARTIAL] else None)
                    for c in CONTROLS}}
    out["passes"] = {a: bool(readable[a] and out["adoption"][a]["keep"] and first["first_uploads"][a]["passes"])
                     for a in ORDER}
    leads = [out["lead"][c] is not None and out["lead"][c]["keep"] for c in CONTROLS]
    out["verdict"] = "matching" if all(leads) else "moving less"
    if out["verdict"] == "matching" and out["passes"][PARTIAL]:
        out["put_forward"] = PARTIAL
    else:
        out["put_forward"] = next((a for a in ORDER if out["passes"][a]), None)
    out["closed"] = out["put_forward"] is None
    return out


def reading(run_dirs) -> dict:
    run_dirs = list(run_dirs)
    committed, why, meta = nt._load([*nt.STORED, E28_RUNS, N45_RUNS])
    run, why_run, meta_run = nt._load(run_dirs)
    why, meta = [*why, *why_run], sorted({*meta, *meta_run})
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (ENTRY, SECOND, PARTIAL):
            if arm not in committed.get(t.key, {}):
                why.append(f"{n}: no committed scores for {arm}")
        if t.key not in run:
            why.append(f"{n}: no task file in the run")
    worst, why_rescored = rescored(committed, run, tasks)
    table, worst_match, why_moved = moved(movement(run_dirs), run, tasks)
    why = [*why, *why_rescored, *why_moved]
    head = {"rescored_max_abs_difference": worst, "movement": table, "matched_max_abs_difference": worst_match}
    if why:
        return {"void": True, "why": sorted(set(why)), **head}
    scores = {t.key: {**committed[t.key], **{a: run[t.key][a] for a in CONTROLS if a in run[t.key]}}
              for t in tasks.values()}
    return {"void": False, "why": [], "runs": meta, **head, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    match = res["matched_max_abs_difference"]
    lines = [f"Readings of tries E42a to E42d (runs {', '.join(res['runs'])}): four controls for fgw_partial on the "
             "nine practice tasks; aggregation of committed scores",
             f"verdict: {res['verdict']}; arm put forward: {res['put_forward'] or 'none'}"
             + (" (to a cold check before any request)" if res["put_forward"] else ""),
             "  rescored arms against their committed scores, largest difference over tasks and seeds: "
             + ", ".join(f"{a} {res['rescored_max_abs_difference'][a]:.1e}" for a in RESCORED)
             + f" (void above {nt.TOL})",
             "  how far the cells move, as a share of the entry's distance (mean over tasks and seeds; earlier "
             "bracket, later bracket): "
             + ", ".join(f"{a} {np.mean([v['cell'][0] for v in res['movement'][a].values()]):.2f} "
                         f"{np.mean([v['cell'][1] for v in res['movement'][a].values()]):.2f}"
                         for a in MOVED if res["movement"][a]),
             f"  {MATCHED} against {PARTIAL}, largest difference of that share over tasks, seeds and brackets: "
             + ("not readable" if match is None else f"{match:.1e}") + f" (void above {MATCH_TOL})"]
    for arm in res["arms"]:
        a = res["adoption"][arm]
        if a is None:
            lines.append(f"  {arm}: not readable")
            continue
        parts = [f"{n} {t['mean_gain']:+.2f} ({t['seeds_ahead']}/10)" for n, t in a["tasks"].items()]
        v = res["first_upload"]["first_uploads"][arm]
        r = v["rule"]
        lines.append(f"  {arm:<14} adoption against {ENTRY}: " + ", ".join(parts)
                     + f" -> {'kept' if a['keep'] else 'not kept'}; first upload S(0.45) {v['S0.45']:+.2f}, clauses "
                     + " ".join("y" if r[k] else "n" for k in sorted(r))
                     + f" -> {'passes both' if res['passes'][arm] else 'does not pass both'}")
        lines.append(f"    own mean minus {ENTRY}: "
                     + " ".join(f"{n} {res['own_minus_fgw'][arm][n]:+.2f}" for n in names))
    for c, a in res["lead"].items():
        if a is None:
            lines.append(f"  {PARTIAL} against {c}: not readable")
            continue
        parts = [f"{n} {t['mean_gain']:+.2f} ({t['seeds_ahead']}/10)" for n, t in a["tasks"].items()]
        lines.append(f"  {PARTIAL} against {c}: " + ", ".join(parts)
                     + f" -> {'leads' if a['keep'] else 'does not lead'}")
        lines.append("    on every task: "
                     + " ".join(f"{n} {res['partial_minus_control'][c][n]:+.2f}" for n in names))
    lines.append("  the controls as a second upload beside fgw (rows E34 and E35's rule), S(0.45): "
                 + ", ".join(f"{c} " + (f"{v['S0.45']:+.2f}" if v["readable"] else "not readable")
                             for c, v in res["second_upload"]["second_uploads"].items()))
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of tries E42")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n46.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n46.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
