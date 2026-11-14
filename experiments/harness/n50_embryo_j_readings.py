#!/usr/bin/env python
"""Reading of try E45 of experiments/harness/TRIES.md (registered before any run): `fgw_damp_m_grp`, the
lineage-level twin of the movement-matched control `fgw_damp_m`, scored on the nine practice tasks by
`n50_embryo_j_arms.py` on the branch `campaign/emb-j`.

Why. Row E42's reading put `fgw_damp_m` forward: the entry `fgw` with every anchor's displacement multiplied by one
factor per bracket, so that the bracket's cells move on average exactly as far as under `fgw_partial`. Its amount
of movement comes from labels (the cell-type labels behind the partial coupling that gives its factors), and
cell-type names are not harmonised across stages. Rows E42 and E43 fixed, before their runs, that the twin with
the share and the partial coupling at the lineage-group level is built and read first, under the same two rules.

What it reads. Only per-seed board scores committed by counted tries and the run's own movement record: the stored
runs of E24, E25 and E28, the committed task files of tries E42 (`results/n46_embryo_g_arms/tasks`: `fgw_damp_m`)
and E43 (`results/n47_embryo_h_arms/tasks`: `fgw_partial_grp`), and the committed task files of try E45 given by
`--runs` (their scores and, under "build_records", each moving arm's ratios to the entry's displacement and the two
damped arms' factors). It opens no data file, runs no candidate and no scorer.

Validity, checks that can fail. (1) The run rescored the floor, `fgw`, `fgw_partial_grp` and `fgw_damp_m`: each must
agree with its committed scores seed for seed within nine_tasks.TOL on every task. (2) The twin is scored on every
task: it uses the entry's coupling and `fgw_partial_grp`'s solver, which build there, so a refusal is a fault. (3)
`fgw_partial_grp` and the two damped arms hold their movement record for every seed. (4) Each damped arm's record
is what the arm is by construction: its anchors and its cells move by its recorded factors, within MATCH_TOL. (5)
The twin is matched: on every task, seed and bracket its cell-level ratio equals `fgw_partial_grp`'s within
MATCH_TOL. Otherwise the reading is void, decides nothing, and is rerun only after its cause is fixed.

Reading, fixed in the row before the run: the adoption rule (vec.adoption.keep_change) for the twin against `fgw` on
the two counted embryo backtests, and row N33's first-upload rule beside the second upload as built. If the twin
passes both, its file is the one a request would name, after an independent recomputation and a cold adversarial
check. If it does not, the line closes with no request and no other arm takes its place (rows E42 and E43).
Reported, deciding nothing: the twin's own mean minus `fgw`'s, minus `fgw_damp_m`'s and minus `fgw_partial_grp`'s on
every task; both damped arms' factors; and, on the E7.25 backtest, the twin's gain on the seeds whose coupling is
mirrored (3, 4, 5, 8 and 9; row E42's outcome) and on the other five.

  python -I -B experiments/harness/n50_embryo_j_readings.py --runs DIR --out DIR
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

from vec.adoption import SEEDS, keep_change  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ENTRY = "fgw"  # the first upload as built
SECOND = nt.RA  # fgw_aligned_series, the second upload as built
GRP = "fgw_partial_grp"  # try E43: the partial coupling with its share from the lineage-group labels
DAMP_M = "fgw_damp_m"  # try E42d: the control matched to fgw_partial's movement (cell-type labels)
ARM = "fgw_damp_m_grp"  # its twin, matched to fgw_partial_grp's movement
RESCORED = ("floor", ENTRY, GRP, DAMP_M)  # scored again by the run; must match the record
MOVED = (GRP, DAMP_M, ARM)  # the arms whose movement record the run holds
DAMPED = (DAMP_M, ARM)
BACKTESTS = ("e75", "e725")  # the two counted embryo backtests among the nine tasks
MIRRORED = ("e725", (3, 4, 5, 8, 9))  # the E7.25 backtest's seeds whose coupling sits mirrored (row E42's outcome)
MATCH_TOL = 1e-6
MOVEMENT_KEY = "movement_ratio"
E28_RUNS = HERE / "results" / "e28_nine_tasks" / "tasks"
N46_RUNS = HERE / "results" / "n46_embryo_g_arms" / "tasks"
N47_RUNS = HERE / "results" / "n47_embryo_h_arms" / "tasks"


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
    """Per task key, arm and seed, what the run recorded of the arm's movement: "ratio", {"anchor": [before, after],
    "cell": [before, after]}, the arm's mean displacement norm over the entry's on the same anchors and cells, and
    "factors", a damped arm's two factors. Each is None where the task file holds none."""
    out: dict = {}
    for d in run_dirs:
        for p in sorted(Path(d).glob("T2e_*.json")):
            r = json.loads(p.read_text())
            recs = r.get("build_records", {})
            out[r["task"]["key"]] = {a: {int(s): {"ratio": v.get(MOVEMENT_KEY), "factors": v.get("factors")}
                                         for s, v in recs.get(a, {}).items()} for a in MOVED}
    return out


def _two(v):
    """Two finite floats (before, after), or None."""
    if not isinstance(v, (list, tuple)) or len(v) != 2 or any(x is None for x in v):
        return None
    v = [float(x) for x in v]
    return v if all(np.isfinite(v)) else None


def moved(moves: dict, run: dict, tasks: dict) -> tuple[dict, dict, float | None, list[str]]:
    """The checks of the movement record. Returns (per arm and task, the mean over seeds of the cell-level ratio per
    bracket; per damped arm and task, its factors as [min, max] over seeds per bracket; the largest absolute
    difference between the twin's cell-level ratio and fgw_partial_grp's; the reasons the run is void)."""
    table: dict = {a: {} for a in MOVED}
    factors: dict = {a: {} for a in DAMPED}
    worst, why = None, []
    for n, t in tasks.items():
        have = moves.get(t.key, {})
        per: dict = {}
        for arm in MOVED:
            if arm not in run.get(t.key, {}):
                continue  # reading() voids a run that did not score it
            recs = {s: have.get(arm, {}).get(s) or {} for s in SEEDS}
            pairs = {s: {lv: _two((recs[s].get("ratio") or {}).get(lv)) for lv in ("cell", "anchor")} for s in SEEDS}
            bad = [s for s in SEEDS if pairs[s]["cell"] is None or pairs[s]["anchor"] is None]
            if bad:
                why.append(f"{n}: no movement record for {arm} (seeds {bad})")
                continue
            per[arm] = pairs
            table[arm][n] = [float(np.mean([pairs[s]["cell"][b] for s in SEEDS])) for b in (0, 1)]
            if arm in DAMPED:
                want = {s: _two(recs[s].get("factors")) for s in SEEDS}
                if any(w is None for w in want.values()):
                    why.append(f"{n}: {arm}'s movement record holds no factors")
                    continue
                off = max(abs(pairs[s][lv][b] - want[s][b]) for s in SEEDS for lv in ("cell", "anchor")
                          for b in (0, 1))
                if not off <= MATCH_TOL:
                    why.append(f"{n}: {arm}'s movement record is {off:.2e} off its factors")
                f = np.array([want[s] for s in SEEDS], dtype=float)
                factors[arm][n] = [[float(f[:, b].min()), float(f[:, b].max())] for b in (0, 1)]
        if ARM in per and GRP in per:
            d = max(abs(per[ARM][s]["cell"][b] - per[GRP][s]["cell"][b]) for s in SEEDS for b in (0, 1))
            worst = d if worst is None else max(worst, d)
            if not d <= MATCH_TOL:
                why.append(f"{n}: {ARM} moves its cells {d:.2e} off {GRP} (ratio to the entry's)")
    return table, factors, worst, why


def check(scores: dict, tasks: dict) -> dict:
    names = [*nt.DEEP, *nt.OTHER]

    def minus(a: str, b: str) -> dict:
        return {n: float(scores[tasks[n].key][a].mean() - scores[tasks[n].key][b].mean()) for n in names}

    first = n37.price(scores, tasks, candidates=(ARM, DAMP_M, GRP))
    task, seeds = MIRRORED
    gain = scores[tasks[task].key][ARM] - scores[tasks[task].key][ENTRY]
    others = [s for s in SEEDS if s not in seeds]
    out = {"arm": ARM,
           "own_minus_fgw": minus(ARM, ENTRY), "own_minus_damp_m": minus(ARM, DAMP_M),
           "own_minus_partial_grp": minus(ARM, GRP),
           "adoption": adoption(scores, tasks, ARM), "adoption_damp_m": adoption(scores, tasks, DAMP_M),
           "first_upload": first,
           "mirrored": {"task": task, "seeds": list(seeds), "gain_on_them": float(gain[list(seeds)].mean()),
                        "gain_on_the_others": float(gain[others].mean())}}
    out["passes"] = bool(out["adoption"]["keep"] and first["first_uploads"][ARM]["passes"])
    out["closed"] = not out["passes"]
    return out


def reading(run_dirs) -> dict:
    run_dirs = list(run_dirs)
    committed, why, meta = nt._load([*nt.STORED, E28_RUNS, N46_RUNS, N47_RUNS])
    run, why_run, meta_run = nt._load(run_dirs)
    why, meta = [*why, *why_run], sorted({*meta, *meta_run})
    tasks = nt.all_tasks()
    for n, t in tasks.items():
        for arm in (ENTRY, SECOND, GRP, DAMP_M):
            if arm not in committed.get(t.key, {}):
                why.append(f"{n}: no committed scores for {arm}")
        if t.key not in run:
            why.append(f"{n}: no task file in the run")
        elif ARM not in run[t.key] or not np.isfinite(run[t.key][ARM]).all():
            why.append(f"{n}: the run did not score {ARM} (it uses fgw's coupling and {GRP}'s solver, which build "
                       "there: a refusal is a fault)")
    worst, why_rescored = rescored(committed, run, tasks)
    table, factors, worst_match, why_moved = moved(movement(run_dirs), run, tasks)
    why = [*why, *why_rescored, *why_moved]
    if worst_match is None and not why:
        why.append(f"{ARM} was never compared with {GRP}")
    head = {"rescored_max_abs_difference": worst, "movement": table, "factors": factors,
            "matched_max_abs_difference": worst_match}
    if why:
        return {"void": True, "why": sorted(set(why)), **head}
    scores = {t.key: {**committed[t.key], ARM: run[t.key][ARM]} for t in tasks.values()}
    return {"void": False, "why": [], "runs": meta, **head, **check(scores, tasks)}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = [*nt.DEEP, *nt.OTHER]
    m = res["mirrored"]
    lines = [f"Reading of try E45 (runs {', '.join(res['runs'])}): the lineage-level twin of the movement-matched "
             "control, on the nine practice tasks; aggregation of committed scores",
             (f"{ARM} passes both rules: it goes to a recomputation and a cold check before any request"
              if res["passes"] else
              f"{ARM} does not pass both rules: the line closes with no request (rows E42 and E43)"),
             "  rescored arms against their committed scores, largest difference over tasks and seeds: "
             + ", ".join(f"{a} {res['rescored_max_abs_difference'][a]:.1e}" for a in RESCORED)
             + f" (void above {nt.TOL})",
             f"  {ARM} against {GRP}, largest difference of the cells' share of the entry's distance over tasks, "
             f"seeds and brackets: {res['matched_max_abs_difference']:.1e} (void above {MATCH_TOL})"]
    for arm, a in ((ARM, res["adoption"]), (DAMP_M, res["adoption_damp_m"])):
        parts = [f"{n} {t['mean_gain']:+.2f} ({t['seeds_ahead']}/10)" for n, t in a["tasks"].items()]
        v = res["first_upload"]["first_uploads"][arm]
        r = v["rule"]
        lines.append(f"  {arm:<15} adoption against {ENTRY}: " + ", ".join(parts)
                     + f" -> {'kept' if a['keep'] else 'not kept'}; first upload S(0.45) {v['S0.45']:+.2f}, clauses "
                     + " ".join("y" if r[k] else "n" for k in sorted(r)))
    for label, key in ((ENTRY, "own_minus_fgw"), (DAMP_M, "own_minus_damp_m"), (GRP, "own_minus_partial_grp")):
        lines.append(f"  {ARM} own mean minus {label}: " + " ".join(f"{n} {res[key][n]:+.2f}" for n in names))
    for arm in DAMPED:
        lines.append(f"  {arm}'s factors, earlier and later bracket, lowest to highest over seeds: "
                     + "; ".join(f"{n} {f[0][0]:.2f}-{f[0][1]:.2f} {f[1][0]:.2f}-{f[1][1]:.2f}"
                                 for n, f in res["factors"][arm].items()))
    lines.append(f"  on {m['task']}, the twin's gain over {ENTRY} on the mirrored seeds {m['seeds']}: "
                 f"{m['gain_on_them']:+.2f}; on the other five: {m['gain_on_the_others']:+.2f}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", action="append", required=True, help="the committed tasks folder of try E45")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty: refusing to overwrite")
    res = reading([Path(d) for d in a.runs])
    state = repo_state(HERE.parents[1])
    res["provenance"] = {"commit": state.commit, "tree_clean": state.unclean == (), "utc": state.utc}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n50.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / "n50.txt").write_text(text(res))
    print(text(res), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
