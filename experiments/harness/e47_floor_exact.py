#!/usr/bin/env python
"""Tries E47a and E47b of experiments/harness/TRIES.md (registered before the run): the embryo's exact floor-form arm
and its series twin on the nine practice tasks of try E28.

Why. The first upload (fgw) scored far below the floor on d2_shape and occupancy_dice (docs/PLAN.md decision log
2026-11-18 (j), (k)). The record's floor-form arm, fgw_shape0, keeps E7.5's geometry only in part: its earlier cells
stay where they are, but each later cell goes to the mean position of its partners, which can pull that half of the
cloud inward. fgw_shape0x (E47a) puts every cell on a real position of the earlier bracket: its own cells where they
are, the later bracket's cells on positions of earlier cells not drawn, by an exact assignment in expression
(vec.t2_build._floor_exact_positions). fgw_shape0x_series (E47b) is the same positions with the series rendering.
Each keeps its control's cells, order and expression (fgw_shape0 and fgw_shape0_series), so only positions differ.

Tasks. nine_tasks.all_tasks(): half, late, new1, x1, x2, x3, new2, e75, e725 (10 seeds, bake protocol, build main).
Tasks that share a truth are merged for every nine-task mean: the mean over truths of the mean over that truth's tasks
(E8.0: half, late; E8.25_late: new1, x1; E8.5: new2, x2, x3; E7.5: e75; E7.25: e725).

Arms. The floor; the reference arms fgw, fgw_shape0, fgw_shape0_series and fgw_aligned_series (committed on every task
by tries E24, E25 and E28, recomputed here only for the guard and the pairing, not counted again); the new arms.

Reading, fixed before the run.
  (i)   Per task and per metric, the new arm minus its control (fgw_shape0x - fgw_shape0, fgw_shape0x_series -
        fgw_shape0_series), paired by seed: mean, seeds ahead, and the seed SD of the paired difference.
  (ii)  Geometry. A new arm replaces its control in the question to the captain if (a) its nine-task mean board score
        is not lower than the control's by more than GEOMETRY_SCORE_BAR, (b) its nine-task mean neighborhood_mmd skill
        is not lower by more than GEOMETRY_NBHD_BAR and (c) on no task its d2_shape or occupancy_dice skill (the task's
        mean over seeds) is more than GEOMETRY_SHAPE_BAR below the control's. Otherwise the control stays.
  (iii) Nothing else is decided. The table prices every arm against fgw per task, both ways: the paired mean, the seeds
        ahead, the better-of-two value beside fgw (mean of max(fgw, arm) minus fgw's mean), the nine-task merged mean,
        the mean gain over the tasks where it gains and the mean loss over those where it loses, worst and best task.
Noise floor, reported: per task the seed SD of each arm's board score and of each paired difference; for the nine-task
merged mean of a paired difference its standard error over seeds (SD of the merged per-seed difference / sqrt(10)).

Void (nothing is read): a task, seed or arm missing; results from different commits or a dirty tree; a reference arm
off its committed scores on a task by more than GUARD_MEAN in the ten-seed mean or by more than GUARD_SEED in any seed;
the committed runs disagreeing among themselves (nine_tasks._load's own check); a new arm's de_score, de_direction,
mmd_u or variogram skill differing from its control's by more than EXPR_TOL in any seed (they hold the same cells).

  run   --task NAME --out DIR   one task per process (half, late, new1, x1, x2, x3, new2, e75, e725); start the nine
                                together so that they record one commit
  read  --runs DIR --out DIR    apply the reading; write e47.json and e47.txt

Refused: a run while TRIES.md has no row E47a, a dirty tree, an existing output.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))
import n12_hygiene_size as n12  # noqa: E402
import nine_tasks as nt  # noqa: E402
import series_world as sw  # noqa: E402

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW = "E47a"
TRIES = HERE / "TRIES.md"
REFERENCE = ("fgw", "fgw_shape0", "fgw_shape0_series", "fgw_aligned_series")
NEW = {"fgw_shape0x": "fgw_shape0", "fgw_shape0x_series": "fgw_shape0_series"}  # arm -> its control
ARMS = ("floor", *REFERENCE, *NEW)
PRICED = ("fgw_shape0", "fgw_shape0_series", "fgw_aligned_series", *NEW)
COMMITTED = (*nt.STORED, HERE / "results" / "e28_nine_tasks" / "tasks")
GUARD_MEAN = 0.1  # board points: a reference arm's ten-seed mean against its committed one
GUARD_SEED = 0.5  # board points: any one seed
EXPR_METRICS = ("de_score", "de_direction", "mmd_u", "variogram")  # they read expression only
EXPR_TOL = 1e-9
GEOMETRY_SCORE_BAR = 0.3  # board points, nine-task merged mean
GEOMETRY_NBHD_BAR = 0.02  # neighborhood_mmd skill, nine-task merged mean
GEOMETRY_SHAPE_BAR = 0.10  # d2_shape and occupancy_dice skill, per task
SHAPE_METRICS = ("d2_shape", "occupancy_dice")


def tasks() -> dict:
    return nt.all_tasks()


def truths(task_map: dict | None = None) -> dict[str, list[str]]:
    """Truth (the target file) -> the task names that score against it, in the tasks' order."""
    out: dict[str, list[str]] = {}
    for name, t in (task_map or tasks()).items():
        out.setdefault(t.target, []).append(name)
    return out


def merged(values: dict[str, float], groups: dict[str, list[str]]) -> float:
    """The nine-task mean with tasks that share a truth merged: the mean over truths of the mean over its tasks."""
    return float(np.mean([np.mean([values[n] for n in names]) for names in groups.values()]))


def merged_seeds(values: dict[str, np.ndarray], groups: dict[str, list[str]]) -> np.ndarray:
    """`merged` seed by seed (each value an array over the same seeds)."""
    return np.mean([np.mean([values[n] for n in names], axis=0) for names in groups.values()], axis=0)


# ------------------------------------------------------------------------------------------------ run
def _state() -> tuple[str, bool]:
    s = repo_state(REPO)
    return s.commit, s.unclean == ()


def cmd_run(a) -> None:
    all_tasks = tasks()
    if a.task not in all_tasks:
        sys.exit(f"unknown task {a.task}; known: {sorted(all_tasks)}")
    if ROW not in n12.tries_rows(TRIES):
        sys.exit(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    if not _state()[1]:
        sys.exit("the tree is not clean: tries E47a and E47b run on a committed tree")
    t = all_tasks[a.task]
    cands = [harness.CANDIDATES[n] for n in ARMS]
    missing = [c.name for c in cands if c.name != harness.FLOOR and not c.applies(t)]
    if missing:
        sys.exit(f"{t.key}: arms that do not apply: {missing}")
    out = Path(a.out)
    path = out / f"{sw.safe(t.key)}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    out.mkdir(parents=True, exist_ok=True)
    res = harness.run_task(t, cands, seeds=SEEDS, data_dir=Path(a.data))
    res["task_source"] = "experiments/harness/e47_floor_exact.py"
    res["run_id"] = out.name
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ read
def per_seed(r: dict) -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    """One task result's board scores and metric skills per arm, as arrays over the seeds in seed order."""
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    arms = list(rows[0]["scores"])
    scores = {a: np.array([x["scores"][a] for x in rows], float) for a in arms}
    skills = {a: {m: np.array([x["skills"][a][m] for x in rows], float) for m in rows[0]["skills"][a]}
              for a in arms if a in rows[0]["skills"]}
    return scores, skills


def committed_scores() -> tuple[dict[str, dict[str, np.ndarray]], list[str]]:
    """Task key -> arm -> committed per-seed board scores, from the runs of tries E24, E25 and E28 (each arm checked
    across the runs that hold it, nine_tasks._load), and the reasons that check fails."""
    scores, why, _ = nt._load(COMMITTED)
    return scores, why


def guard(run: dict[str, dict[str, np.ndarray]], committed: dict[str, dict[str, np.ndarray]],
          keys: list[str], arms=REFERENCE) -> tuple[dict, list[str]]:
    """The reference arms recomputed in this run against their committed scores, task by task: the largest absolute
    difference of the ten-seed mean and of any seed, and the reasons the run is void."""
    worst = {"mean": 0.0, "seed": 0.0}
    why = []
    for key in keys:
        for arm in arms:
            old, new = committed.get(key, {}).get(arm), run.get(key, {}).get(arm)
            if old is None or new is None:
                why.append(f"{key}: {arm} has no {'committed' if old is None else 'recomputed'} scores")
                continue
            if old.shape != new.shape:
                why.append(f"{key}: {arm} has {len(new)} seeds against {len(old)} committed")
                continue
            dm, ds = abs(float(new.mean() - old.mean())), float(np.abs(new - old).max())
            worst = {"mean": max(worst["mean"], dm), "seed": max(worst["seed"], ds)}
            if not dm <= GUARD_MEAN:
                why.append(f"{key}: {arm}'s mean is {dm:.3f} off its committed one (tolerance {GUARD_MEAN})")
            if not ds <= GUARD_SEED:
                why.append(f"{key}: {arm} is {ds:.3f} off its committed score in a seed (tolerance {GUARD_SEED})")
    return worst, why


def expression_identity(skills: dict[str, dict[str, dict[str, np.ndarray]]]) -> tuple[float, list[str]]:
    """Each new arm against its control on the metrics that read expression only: the largest difference and the
    reasons the run is void."""
    worst, why = 0.0, []
    for key, sk in skills.items():
        for arm, ctl in NEW.items():
            for m in EXPR_METRICS:
                d = float(np.abs(sk[arm][m] - sk[ctl][m]).max())
                worst = max(worst, d)
                if not d <= EXPR_TOL:
                    why.append(f"{key}: {arm}'s {m} differs from {ctl}'s by {d:.2e}")
    return worst, why


def _diff(x: np.ndarray, y: np.ndarray) -> dict:
    d = x - y
    return {"mean": float(d.mean()), "ahead": int((d > 0).sum()), "sd": float(d.std(ddof=1))}


def geometry_rule(per: dict, groups: dict[str, list[str]], arm: str) -> dict:
    """Rule (ii) for one new arm against its control, from the per-task paired differences of `paired`."""
    ctl = NEW[arm]
    board = merged({n: per[n]["paired"][arm]["board"]["mean"] for n in per}, groups)
    nbhd = merged({n: per[n]["paired"][arm]["skills"]["neighborhood_mmd"]["mean"] for n in per}, groups)
    shape = {n: {m: per[n]["paired"][arm]["skills"][m]["mean"] for m in SHAPE_METRICS} for n in per}
    worst_task, worst_metric = min(((n, m) for n in shape for m in SHAPE_METRICS), key=lambda nm: shape[nm[0]][nm[1]])
    a = bool(board >= -GEOMETRY_SCORE_BAR)
    b = bool(nbhd >= -GEOMETRY_NBHD_BAR)
    c = bool(all(v >= -GEOMETRY_SHAPE_BAR for s in shape.values() for v in s.values()))
    return {"control": ctl, "board_merged": board, "nbhd_merged": nbhd, "shape_per_task": shape,
            "worst_shape": {"task": worst_task, "metric": worst_metric, "value": shape[worst_task][worst_metric]},
            "a_board": a, "b_nbhd": b, "c_shape": c, "replaces": bool(a and b and c),
            "outcome": (f"{arm} replaces {ctl} in the question to the captain" if a and b and c
                        else f"{ctl} stays")}


def price(scores: dict[str, dict[str, np.ndarray]], groups: dict[str, list[str]], arm: str) -> dict:
    """Rule (iii)'s row for one arm against fgw, both ways."""
    per = {}
    for n, s in scores.items():
        d = _diff(s[arm], s["fgw"])
        per[n] = {**d, "pair_value": float(np.maximum(s[arm], s["fgw"]).mean() - s["fgw"].mean())}
    means = {n: v["mean"] for n, v in per.items()}
    gains, losses = [v for v in means.values() if v > 0], [v for v in means.values() if v < 0]
    return {"tasks": per, "merged": merged(means, groups),
            "pair_value_merged": merged({n: v["pair_value"] for n, v in per.items()}, groups),
            "gain_tasks": len(gains), "mean_gain": float(np.mean(gains)) if gains else 0.0,
            "loss_tasks": len(losses), "mean_loss": float(np.mean(losses)) if losses else 0.0,
            "worst": min(means.items(), key=lambda kv: kv[1]), "best": max(means.items(), key=lambda kv: kv[1])}


def reading(results: dict[str, dict], committed: tuple[dict, list[str]] | None = None) -> dict:
    """The reading of fixed rule (i) to (iii) from one run's task results (task key -> run_task's result)."""
    all_tasks = tasks()
    groups = truths(all_tasks)
    why = []
    for name, t in all_tasks.items():
        r = results.get(t.key)
        if r is None:
            why.append(f"missing task {name}")
            continue
        if sorted(row["seed"] for row in r["rows"]) != list(SEEDS):
            why.append(f"{name}: seeds are not {list(SEEDS)}")
        if any(a not in row["scores"] or a not in row["skills"] for row in r["rows"] for a in ARMS):
            why.append(f"{name}: an arm is missing")
    prov = [r["provenance"] for r in results.values()]
    if len({p["commit"] for p in prov}) > 1:
        why.append("the results come from different commits")
    if not all(p["tree_clean"] for p in prov):
        why.append("a result comes from a dirty tree")
    if why:
        return {"void": True, "why": why}
    data = {n: per_seed(results[t.key]) for n, t in all_tasks.items()}
    scores = {n: d[0] for n, d in data.items()}
    skills = {n: d[1] for n, d in data.items()}
    old, old_why = committed if committed is not None else committed_scores()
    worst_guard, guard_why = guard({t.key: scores[n] for n, t in all_tasks.items()}, old,
                                   [t.key for t in all_tasks.values()])
    worst_expr, expr_why = expression_identity(skills)
    why = [*old_why, *guard_why, *expr_why]
    if why:
        return {"void": True, "why": why, "guard": worst_guard, "expression_identity": worst_expr}
    per = {}
    for n in all_tasks:
        s, k = scores[n], skills[n]
        per[n] = {"mean_score": {a: float(v.mean()) for a, v in s.items()},
                  "seed_sd": {a: float(v.std(ddof=1)) for a, v in s.items()},
                  "paired": {arm: {"board": _diff(s[arm], s[ctl]),
                                   "skills": {m: _diff(k[arm][m], k[ctl][m]) for m in k[arm]}}
                             for arm, ctl in NEW.items()}}
    noise = {arm: {"merged_board_mean": merged({n: per[n]["paired"][arm]["board"]["mean"] for n in per}, groups),
                   "merged_board_se": float(merged_seeds({n: scores[n][arm] - scores[n][ctl] for n in per}, groups)
                                            .std(ddof=1) / np.sqrt(len(SEEDS)))}
             for arm, ctl in NEW.items()}
    rule2 = {arm: geometry_rule(per, groups, arm) for arm in NEW}
    pricing = {arm: price(scores, groups, arm) for arm in PRICED}
    return {"void": False, "why": [], "commit": prov[0]["commit"], "truths": groups, "tasks": per, "noise": noise,
            "rule2": rule2, "pricing": pricing, "guard": worst_guard, "expression_identity": worst_expr,
            "thresholds": {"GUARD_MEAN": GUARD_MEAN, "GUARD_SEED": GUARD_SEED, "EXPR_TOL": EXPR_TOL,
                           "GEOMETRY_SCORE_BAR": GEOMETRY_SCORE_BAR, "GEOMETRY_NBHD_BAR": GEOMETRY_NBHD_BAR,
                           "GEOMETRY_SHAPE_BAR": GEOMETRY_SHAPE_BAR}}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    names = list(res["tasks"])
    lines = [f"Tries E47a and E47b (commit {res['commit']}); guard: reference arms within {res['guard']['mean']:.3f} "
             f"(mean) and {res['guard']['seed']:.3f} (any seed) of their committed scores; expression metrics equal to "
             f"the controls' within {res['expression_identity']:.1e}.",
             "rule (ii), geometry:"]
    for arm, g in res["rule2"].items():
        w = g["worst_shape"]
        lines.append(f"  {arm} against {g['control']}: board {g['board_merged']:+.3f} (a {g['a_board']}), "
                     f"neighborhood_mmd {g['nbhd_merged']:+.4f} (b {g['b_nbhd']}), worst shape {w['value']:+.3f} "
                     f"({w['metric']} on {w['task']}; c {g['c_shape']}): {g['outcome']}")
        lines.append(f"    noise: merged paired mean {res['noise'][arm]['merged_board_mean']:+.3f}, its standard error "
                     f"over seeds {res['noise'][arm]['merged_board_se']:.3f}")
    lines.append("rule (i), the new arm minus its control per task (board mean, seeds ahead, seed SD; then skills):")
    for n in names:
        for arm, p in res["tasks"][n]["paired"].items():
            b = p["board"]
            lines.append(f"  {n:<5} {arm:<19} {b['mean']:+6.2f} ({b['ahead']:>2}/10, sd {b['sd']:.2f})  "
                         + " ".join(f"{m} {v['mean']:+.3f}" for m, v in p["skills"].items()))
    lines.append("rule (iii), each arm against fgw (paired mean per task; merged; gain / loss tasks and means; better "
                 "of two beside fgw, merged):")
    lines.append("  arm                  " + " ".join(f"{n:>6}" for n in names)
                 + "  merged  gains        losses       pair")
    for arm, p in res["pricing"].items():
        lines.append(f"  {arm:<20} " + " ".join(f"{p['tasks'][n]['mean']:+6.2f}" for n in names)
                     + f"  {p['merged']:+6.2f}  {p['gain_tasks']} {p['mean_gain']:+5.2f}   "
                     f"{p['loss_tasks']} {p['mean_loss']:+5.2f}   {p['pair_value_merged']:+5.2f}")
    lines.append("seed SD of each arm's board score per task:")
    for n in names:
        lines.append(f"  {n:<5} " + ", ".join(f"{a} {v:.2f}" for a, v in res["tasks"][n]["seed_sd"].items()))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    keys = {t.key for t in tasks().values()}
    results = {}
    for p in sorted(Path(a.runs).glob("T2e_*.json")):
        r = json.loads(p.read_text())
        if r.get("task", {}).get("key") in keys:
            results[r["task"]["key"]] = r
    out = Path(a.out)
    files = [out / "e47.json", out / "e47.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    res = reading(results)
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(res, indent=1, default=float))
    files[1].write_text(text(res))
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
