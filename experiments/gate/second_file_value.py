#!/usr/bin/env python3
"""What a second file is worth on each test board, from committed practice results only.

One official upload is left on each board, and a board keeps the better of a team's two. A second file is therefore
worth V = E[max(D, 0)], where D is its score minus the banked first file's. On a practice task D is read as a PAIRED
difference: the candidate minus the practice analogue of the file actually uploaded, same task, same seed. The banked
test score then drops out. No test score enters any number here: the first uploads only say which file is banked
(docs/PLAN.md decision log 2026-11-18 (j)), and the lines of entry (k) hold (no constant from a test value, no task
chosen or weighted by its closeness to one).

Not a try (experiments/harness/TRIES.md, rule 1): it aggregates committed per-seed practice scores. It loads no
competition file, builds nothing, runs no scorer and compares no predictor with a target. Every file is read through
`Reader`, which refuses boards/, data/, runs/, docs/ and any notes/ folder, and records each file's sha256.

Per board and candidate it gives:
  D per task          the seed-paired mean of candidate - first file over the task's seeds, its seed SD (n - 1) and
                      the number of seeds ahead;
  units               tasks that share one truth (the same target file) are merged into one unit; the unit's mean D is
                      the mean of its tasks' mean D;
  across units        the mean and the SD of the units' mean D (n - 1 with equal weights; the weighted population SD
                      under a stated prior);
  V, task level       the (weighted) mean over units of max(unit mean D, 0);
  V, seed level       the (weighted) mean over units of the mean over the unit's tasks and seeds of max(D, 0): the
                      expected gain of the better of two when one scoring draw of the unit is the test;
  share positive      the (weighted) share of units with mean D > 0.
Weights are equal unless a committed file states a prior; then every stated prior is shown and the file named. No
prior is made up here.

Boards:
  embryo  first file fgw; nine practice tasks of try E28 (per-seed scores from runs E24, E25 and E28, merged as
          experiments/harness/nine_tasks.py merges them); five truths. Also try E28's own weighting S(p) = p D + (1 -
          p) O over its unmerged groups, at its three stated p.
  T1      first file lib3_base; try T33's twelve zero-gap scenarios (T33's own arms, T32's stored s1_cm60 and
          lib3_cm60, as T33 reads them) under T30's nine stated priors, imported by reading t30_pricing.py's
          constants. The one-day reading that exists is the oracle-base grid (s2_hedge t1_grid, E9.5 -> E10.5): it
          has no library-3 or mid base, so it prices only the cardiomyocyte step on the pooled base (bown_c0.6 against
          bown_cown) under its hedge file's five stated priors. Try T31 (one day, bases only) is shown as context.
  heart   first file b0_trend, candidate copy, on the two counted backtests of trend_heart_v1 (one truth).
  T3      no paired analogue of the allowed s1 exists (the script checks where its stand-ins were scored); the
          candidates are given against the floor on K7's and K8's folds, and the first file's stand-ins against the
          floor on their own folds.

Usage, from the repository root (standard library only):
  python3 -I -B experiments/gate/second_file_value.py \
      --json experiments/gate/results/second_file_value_2026-11-18.json \
      --txt experiments/gate/results/second_file_value_2026-11-18.txt
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
import math
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SEEDS = list(range(10))
TOL_RUNS = 0.01  # nine_tasks.TOL: an arm scored in two runs agrees seed for seed to this (board points)
TOL_SUMMARY = 1e-9  # a committed summary recomputed from its per-seed files


# ====================================================================================================== the reader
class Reader:
    """Reads committed result files only, and records the sha256 of each."""

    ALLOWED = ("experiments/", "results/registry.csv")
    REFUSED = ("boards", "data", "runs", "notes", "research", "reference", ".venv", ".venv-gpu", ".git", ".claude")

    def __init__(self, root: Path = REPO):
        self.root = Path(root)
        self.read: dict[str, str] = {}

    def check(self, rel: str) -> str:
        rel = Path(rel).as_posix()
        parts = rel.split("/")
        if (rel.startswith("/") or ".." in parts or not rel.startswith(self.ALLOWED)
                or any(p in self.REFUSED for p in parts) or "index.json" in parts):
            raise PermissionError(f"second_file_value: {rel} is not a committed result this script may read")
        return rel

    def bytes(self, rel: str) -> bytes:
        rel = self.check(rel)
        b = (self.root / rel).read_bytes()
        self.read[rel] = hashlib.sha256(b).hexdigest()
        return b

    def json(self, rel: str):
        return json.loads(self.bytes(rel))

    def text(self, rel: str) -> str:
        return self.bytes(rel).decode()

    def csv(self, rel: str) -> list[dict]:
        return list(csv.DictReader(io.StringIO(self.text(rel))))

    def glob(self, rel_dir: str, pattern: str) -> list[str]:
        rel_dir = self.check(rel_dir.rstrip("/") + "/")
        return sorted(p.relative_to(self.root).as_posix() for p in (self.root / rel_dir).glob(pattern))


# ==================================================================================================== arithmetic
def need(cond: bool, what: str, checks: list[str] | None = None) -> None:
    """A self-check: stop if it fails, and record it."""
    if not cond:
        raise AssertionError(f"second_file_value: check failed: {what}")
    if checks is not None:
        checks.append(what)


def paired(first: list[float], cand: list[float]) -> list[float]:
    """D per seed: candidate minus first file, same seed."""
    need(len(first) == len(cand) and len(first) > 0, "two arms with the same seeds")
    return [c - f for f, c in zip(first, cand, strict=True)]


def seed_summary(d: list[float]) -> dict:
    return {"mean": statistics.fmean(d), "seed_sd": statistics.stdev(d) if len(d) > 1 else 0.0,
            "ahead": sum(x > 0 for x in d), "seeds": len(d)}


def truth_units(targets: dict[str, str]) -> dict[str, list[str]]:
    """Tasks that share one target file merged into one unit, keyed by the target, in first-seen order."""
    units: dict[str, list[str]] = {}
    for task, target in targets.items():
        units.setdefault(target, []).append(task)
    return units


def unit_values(d_by_task: dict[str, list[float]], units: dict[str, list[str]]) -> dict[str, dict]:
    """Per unit: the mean of its tasks' mean D and the mean over its tasks and seeds of max(D, 0)."""
    out = {}
    for u, tasks in units.items():
        means = [statistics.fmean(d_by_task[t]) for t in tasks]
        out[u] = {"tasks": list(tasks), "mean_d": statistics.fmean(means),
                  "seed_v": statistics.fmean([max(x, 0.0) for t in tasks for x in d_by_task[t]])}
    return out


def summarise(units: dict[str, dict], weights: dict[str, float] | None = None) -> dict:
    """Across units: mean and SD of mean D, V at task and seed level, and the share of units with mean D > 0."""
    names = list(units)
    if weights is None:
        w = {u: 1.0 / len(names) for u in names}
    else:
        tot = sum(weights[u] for u in names)
        need(tot > 0, "a prior with positive weight")
        w = {u: weights[u] / tot for u in names}
    m = [units[u]["mean_d"] for u in names]
    mean = sum(w[u] * units[u]["mean_d"] for u in names)
    if len(m) < 2:
        sd = None  # one unit: no spread across units
    elif weights is None:
        sd = statistics.stdev(m)
    else:
        sd = math.sqrt(sum(w[u] * (units[u]["mean_d"] - mean) ** 2 for u in names))
    return {"units": len(names), "mean_d": mean, "sd_d": sd, "sd_kind": "n-1" if weights is None else "weighted",
            "v_task": sum(w[u] * max(units[u]["mean_d"], 0.0) for u in names),
            "v_seed": sum(w[u] * units[u]["seed_v"] for u in names),
            "share_positive": sum(w[u] for u in names if units[u]["mean_d"] > 0)}


def candidate_table(first: dict[str, list[float]], cand: dict[str, list[float]], targets: dict[str, str]) -> dict:
    """One candidate against the first file on every task both hold: per task, per unit, and the summary."""
    tasks = [t for t in targets if t in first and t in cand]
    d = {t: paired(first[t], cand[t]) for t in tasks}
    units = unit_values(d, truth_units({t: targets[t] for t in tasks}))
    return {"per_task": {t: seed_summary(d[t]) for t in tasks}, "units": units, "summary": summarise(units), "d": d}


# ===================================================================================== constants of a committed file
def _const(node):
    """A literal, or arithmetic on literals (1 / 7, (0.25,) * 4, -x): the forms t30_pricing.py's priors use."""
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float, str)):
        return node.value
    if isinstance(node, (ast.Tuple, ast.List)):
        return tuple(_const(e) for e in node.elts)
    if isinstance(node, ast.Dict):
        return {_const(k): _const(v) for k, v in zip(node.keys, node.values, strict=True)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_const(node.operand)
    if isinstance(node, ast.BinOp):
        a, b = _const(node.left), _const(node.right)
        ops = {ast.Div: lambda x, y: x / y, ast.Mult: lambda x, y: x * y, ast.Add: lambda x, y: x + y,
               ast.Sub: lambda x, y: x - y}
        if type(node.op) in ops:
            return ops[type(node.op)](a, b)
    raise ValueError(f"not a constant expression: {ast.dump(node)}")


def module_constants(source: str, names: tuple[str, ...]) -> dict:
    """Top-level NAME = <constant expression> assignments of a module's source, without importing it."""
    out = {}
    for node in ast.parse(source).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in names:
                out[node.targets[0].id] = _const(node.value)
    missing = [n for n in names if n not in out]
    if missing:
        raise ValueError(f"constants not found: {missing}")
    return out


# =================================================================================================== task files
def task_scores(rd: Reader, rel: str, checks: list[str]) -> dict:
    """A committed harness or diagnostic task file: key, target, scorer build and each arm's score per seed."""
    r = rd.json(rel)
    rows = sorted(r["rows"], key=lambda x: x["seed"])
    need([int(x["seed"]) for x in rows] == SEEDS, f"{rel}: seeds 0 to 9", checks)
    prov = r.get("provenance", {})
    need(prov.get("tree_clean", True) is True, f"{rel}: a clean tree", checks)
    t = r["task"]
    return {"key": t["key"] if isinstance(t, dict) else t, "target": t.get("target") if isinstance(t, dict) else None,
            "reference": t.get("reference") if isinstance(t, dict) else None, "build": r.get("build"),
            "arms": {a: [float(x["scores"][a]) for x in rows] for a in rows[0]["scores"]}}


# ========================================================================================================= embryo
EMB_RUNS = ("experiments/harness/results/e24_pricing/tasks", "experiments/harness/results/e25_pairs/tasks",
            "experiments/harness/results/e28_nine_tasks/tasks")
EMB_SUMMARY = "experiments/harness/results/e28_nine_tasks/nine_tasks.json"
# task names: experiments/harness/series_world.py TASKS, s2_pricing.py COUNTED, second_jump.py TASKS
EMB_TASKS = {"half": "T2e/E7.5+E8.5->E8.0", "late": "T2e/E7.5+E8.25_late->E8.0", "new1": "T2e/E7.5+E8.5->E8.25_late",
             "x1": "T2e/E8.0+E8.75->E8.25_late", "x2": "T2e/E8.0+E8.75->E8.5", "x3": "T2e/E8.25_late+E8.75->E8.5",
             "new2": "T2e/E7.5+E8.75->E8.5", "e75": "T2e/E7.25+E8.0->E7.5", "e725": "T2e/E6.75+E7.5->E7.25"}
EMB_DEEP = ("half", "late", "new1", "x1", "x2", "x3")  # nine_tasks.py DEEP
EMB_OTHER = ("new2", "e75", "e725")  # nine_tasks.py OTHER
EMB_PS = (0.3, 0.45, 0.6)  # nine_tasks.py PS; TRIES E28; forecast.py JUDGEMENT emb_p_deep (0.45; 0.30 to 0.60)
EMB_FIRST = "fgw"
EMB_CANDS = ("fgw_aligned_series", "fgw_shape0", "fgw_shape0_series", "fgw_series", "fgw_render_hi", "fgw_shape0_hi",
             "fgw_aligned")
EMB_NOTES = {
    "fgw_aligned_series": "the built second file (builds/T2embryo_s2.h5ad, ruling R4); named by try E28's rule (1) on "
                          "these nine tasks, so all nine were used to select it; its values on x1 to x3 were known "
                          "from diagnostic E26 before E28's rules were written",
    "fgw_shape0": "the frozen s2 before ruling R4 (E7.5's geometry, the cover for fgw's flattened cloud)",
    "fgw_shape0_series": "try E24's first pick on E24's six tasks (half, late, new1, new2, e75, e725)",
    "fgw_series": "its values on x1 to x3 were known from diagnostic E26 before E28's rules were written",
    "fgw_render_hi": "its values on x1 to x3 were known from diagnostic E26 before E28's rules were written",
    "fgw_shape0_hi": "first scored on x1 to x3 by try E28 (counted)",
    "fgw_aligned": "try E11's candidate, not kept (E11's veto stands); read here as a second file; first scored on "
                   "x1 to x3 by try E28",
}
EMB_LIMITS = [
    "first file: fgw was adopted on the two counted backtests e75 and e725 (decision log 2026-11-04); the other seven "
    "tasks count for no adoption",
    "every score here was taken under the practice scorer build 'main'; decision log (k) names the live board's "
    "pre-fix build as the one recorded factor that acts on occupancy_dice alone, not established",
    "after the fact, decision log (k), changing no number here: at the test pair fgw's partners are mirrored in 10 of "
    "10 draws (rows E11, N1, N51), while the practice tasks carry that defect only in some seeds (row N52 reads them "
    "by the flag); occupancy_dice on the first upload is explained by no practice task ((k): a practice task with "
    "mirrored, flatter draws scored 50.1 on it). (k) names no embryo practice task as misleading on expression",
]


def board_embryo(rd: Reader, checks: list[str]) -> dict:
    scores: dict[str, dict[str, list[float]]] = {}
    meta: dict[str, dict] = {}
    worst = 0.0
    for d in EMB_RUNS:
        for rel in rd.glob(d, "T2e_*.json"):
            t = task_scores(rd, rel, checks)
            meta[t["key"]] = {"target": t["target"], "build": t["build"]}
            for arm, v in t["arms"].items():
                old = scores.setdefault(t["key"], {}).get(arm)
                if old is not None:
                    diff = max(abs(a - b) for a, b in zip(old, v, strict=True))
                    worst = max(worst, diff)
                    need(diff <= TOL_RUNS, f"embryo {t['key']} {arm}: runs agree to {TOL_RUNS}")
                scores[t["key"]][arm] = v  # the later run's scores, as nine_tasks._load keeps them
    summ = rd.json(EMB_SUMMARY)
    need(summ["void"] is False, "E28's reading is not void", checks)
    per: dict[str, dict[str, list[float]]] = {}
    for name, key in EMB_TASKS.items():
        need(key in scores, f"embryo: task {name} ({key}) is in the runs")
        per[name] = scores[key]
        for arm in (EMB_FIRST, *EMB_CANDS):
            need(abs(statistics.fmean(scores[key][arm]) - summ["single"][arm][name]) < TOL_SUMMARY,
                 f"embryo {name} {arm}: per-seed scores give nine_tasks.json's mean")
    checks.append("embryo: every arm on every task gives nine_tasks.json's mean (1e-9)")
    targets = {n: meta[k]["target"] for n, k in EMB_TASKS.items()}
    first = {n: per[n][EMB_FIRST] for n in EMB_TASKS}
    cands = {}
    for c in EMB_CANDS:
        tab = candidate_table(first, {n: per[n][c] for n in EMB_TASKS}, targets)
        s = {}
        for p in EMB_PS:  # try E28's weighting, over its unmerged groups
            deep = [max(tab["per_task"][n]["mean"], 0.0) for n in EMB_DEEP]
            other = [max(tab["per_task"][n]["mean"], 0.0) for n in EMB_OTHER]
            deep_s = [statistics.fmean([max(x, 0.0) for x in tab["d"][n]]) for n in EMB_DEEP]
            other_s = [statistics.fmean([max(x, 0.0) for x in tab["d"][n]]) for n in EMB_OTHER]
            s[f"{p:g}"] = {"v_task": p * statistics.fmean(deep) + (1 - p) * statistics.fmean(other),
                           "v_seed": p * statistics.fmean(deep_s) + (1 - p) * statistics.fmean(other_s)}
        tab.pop("d")
        cands[c] = {**tab, "e28_weighting": s, "note": EMB_NOTES[c]}
    return {"first": EMB_FIRST, "tasks": dict(EMB_TASKS), "targets": targets,
            "units": truth_units(targets), "builds": sorted({m["build"] for m in meta.values()}),
            "runs_largest_difference": worst, "candidates": cands, "limits": EMB_LIMITS,
            "merge_source": "the task files' target field (task.target); experiments/harness/series_world.py says "
                            "half and late share one truth; docs/PLAN.md decision log 2026-11-13 (gg), item 2, says "
                            "x2 and x3 do",
            "e28_weighting_source": "experiments/harness/nine_tasks.py (DEEP, OTHER, PS) and TRIES row E28: S(p) = "
                                    "p D + (1 - p) O over unmerged tasks; p 0.45 central, 0.30 and 0.60 its range "
                                    "(experiments/gate/forecast.py JUDGEMENT emb_p_deep)"}


# ============================================================================================================= T1
T30_SRC = "experiments/t1_scope/t30_pricing.py"
T33_DIR = "experiments/t1_scope/results/t33_mid"
T32_DIR = "experiments/t1_scope/results/t32_pair"
T1_FIRST = "lib3_base"  # row R6: the first upload on library 3's cells
T1_CANDS = ("mid_cm60", "s1_cm60", "lib3_cm60", "mid_base")
T1_STORED = ("s1_cm60", "lib3_cm60")  # t33_mid.py STORED: read from T32's scenario files
T1_PAIR_OF = {"mid_cm60": "mid_s2", "s1_cm60": "built", "lib3_cm60": "lib3_pair"}  # T33's pairs with lib3_base first
T1_NOTES = {
    "mid_cm60": "the built second file (builds/T1_s2.h5ad, row R7); try T33 named it by its rule on these twelve "
                "scenarios, so all twelve were used to select it",
    "s1_cm60": "the pooled base at cardiomyocytes 0.60, the second file of row R6 before R7",
    "lib3_cm60": "try T32's second file; its pair did not qualify under T33's rule",
    "mid_base": "the mid base at its own mix; T33 priced it only as a first upload (mid_both)",
}
T1_LIMITS = [
    "zero time gap: T30's own words, 'sizes are upper bounds for the test, where two days of development lie between "
    "base and truth'",
    "no arm carries the uploaded files' expression step (M1, Feng's direction): every arm draws real E10.5 cells, so "
    "D prices the base and the cell mix only",
    "the twelve scenarios are one truth (one half of E10.5's library 3) thinned to four depths and re-drawn to three "
    "cardiomyocyte shares: they are weighted by T30's priors, not merged",
    "after the fact, decision log (k), changing no number here: the counted practice cuts (not used here) could not "
    "see scope or composition and under-called the first file; (k) names none of T33's scenarios",
]
GRID_RUNS = "experiments/s2_hedge/results/per_seed/t1-grid"
GRID_HEDGE = "experiments/s2_hedge/results/t1_grid/hedge_bown_cown.json"
GRID_FIRST, GRID_CAND = "bown_cown", "bown_c0.6"  # experiments/gate/estimate.py T1_S1, T1_S2
T31 = "experiments/t1_scope/results/t31_backward/backward.json"


def _t1_arms(rd: Reader, checks: list[str]) -> tuple[dict, dict, dict]:
    t33 = rd.json(f"{T33_DIR}/pricing.json")
    need(t33["void"] is False, "T33's reading is not void", checks)
    c = module_constants(rd.text(T30_SRC), ("DEPTHS", "CMS", "DEPTH_PRIORS", "CM_PRIORS", "MAIN", "CHECKS"))
    arms: dict[str, dict[str, list[float]]] = {a: {} for a in (T1_FIRST, *T1_CANDS)}
    for d in c["DEPTHS"]:
        for cm in c["CMS"]:
            scen = f"d{d} cm{cm}"
            own = rd.json(f"{T33_DIR}/scenarios/t33_d{d}_cm{cm}.json")
            stored = rd.json(f"{T32_DIR}/scenarios/t32_d{d}_cm{cm}.json")
            for src in (own, stored):
                need([int(x["seed"]) for x in sorted(src["rows"], key=lambda x: x["seed"])] == SEEDS,
                     f"T1 {scen}: seeds 0 to 9")
                need(src["tree_clean"] is True, f"T1 {scen}: a clean tree")
            o = {x["seed"]: x["scores"] for x in own["rows"]}
            s = {x["seed"]: x["scores"] for x in stored["rows"]}
            for a in arms:
                arms[a][scen] = [float((s if a in T1_STORED else o)[k][a]) for k in SEEDS]
                need(abs(statistics.fmean(arms[a][scen]) - t33["arm_by_scenario"][a][scen]) < TOL_SUMMARY,
                     f"T1 {scen} {a}: per-seed scores give T33's pricing.json mean")
    checks.append("T1: every arm in every scenario gives T33's pricing.json mean (1e-9)")
    priors = {}
    for dp, dw in c["DEPTH_PRIORS"].items():
        for cp, cw in c["CM_PRIORS"].items():
            priors[f"{dp} x {cp}"] = {f"d{d} cm{cm}": dw[i] * cw[j] for i, d in enumerate(c["DEPTHS"])
                                      for j, cm in enumerate(c["CMS"])}
    return arms, priors, {"t33": t33, "constants": c}


def board_t1(rd: Reader, checks: list[str]) -> dict:
    arms, priors, ctx = _t1_arms(rd, checks)
    t33, c = ctx["t33"], ctx["constants"]
    scen = list(arms[T1_FIRST])
    targets = {s: s for s in scen}  # one unit per scenario; the priors weight them
    cands = {}
    for cand in T1_CANDS:
        tab = candidate_table(arms[T1_FIRST], arms[cand], targets)
        tab.pop("d")
        by_prior = {p: summarise(tab["units"], w) for p, w in priors.items()}
        pair = T1_PAIR_OF.get(cand)
        if pair:  # V at seed level + the first arm's expected score is T33's expected better of two for that pair
            for p, w in priors.items():
                first_exp = sum(w[s] * statistics.fmean(arms[T1_FIRST][s]) for s in scen)
                need(abs(by_prior[p]["v_seed"] + first_exp - t33["expected_pair"][pair][p]) < 1e-9,
                     f"T1 {cand} {p}: seed-level V plus lib3_base's mean is T33's expected {pair}")
            checks.append(f"T1 {cand}: seed-level V reproduces T33's expected better of two of {pair} (nine priors)")
        cands[cand] = {**tab, "by_prior": by_prior, "note": T1_NOTES[cand]}
    # one day: the oracle-base grid (pooled base only)
    merged: dict[int, dict] = {}
    head = None
    for rel in rd.glob(GRID_RUNS, "*/*.json"):
        g = rd.json(rel)
        this = (g["task"], tuple(g["truth_values"]), tuple(g["pred_values"]))
        need(head is None or this == head, f"{rel}: a part of one grid")
        head = this
        for row in g["rows"]:
            merged.setdefault(int(row["seed"]), {}).update(row["scores"])  # parts split the truths (estimate.load_runs)
    need(sorted(merged) == SEEDS, "T1 grid: seeds 0 to 9", checks)
    truths = list(head[1])
    first = {t: [float(merged[k][t][GRID_FIRST]) for k in SEEDS] for t in truths}
    cand = {t: [float(merged[k][t][GRID_CAND]) for k in SEEDS] for t in truths}
    gtab = candidate_table(first, cand, {t: t for t in truths})
    gtab.pop("d")
    hedge = rd.json(GRID_HEDGE)
    need(hedge["s1"] == GRID_FIRST and hedge["seeds"] == SEEDS, "T1 grid: the hedge file's s1 and seeds", checks)
    g_prior = {}
    for p, w in hedge["priors"].items():
        g_prior[p] = summarise({t: gtab["units"][t] for t in w}, w)
        first_exp = sum(w[t] * statistics.fmean(first[t]) for t in w) / sum(w.values())
        need(abs(g_prior[p]["v_seed"] + first_exp - hedge["candidates"][GRID_CAND][p]["mean"]) < 1e-9,
             f"T1 grid {p}: seed-level V plus bown_cown's mean is the hedge file's expected better of two")
    checks.append("T1 grid: seed-level V reproduces hedge_bown_cown.json's expected better of two (five priors)")
    t31 = rd.json(T31)
    rows = sorted(t31["rows"], key=lambda x: x["seed"])
    d31 = paired([x["scores"]["pooled_base"] for x in rows], [x["scores"]["lib3_base"] for x in rows])
    need(abs(statistics.fmean(d31) - t31["readings"]["lib3_base - pooled_base"]["mean"]) < TOL_SUMMARY,
         "T31: lib3_base - pooled_base is backward.json's reading", checks)
    return {"first": T1_FIRST, "scenarios": scen, "priors": priors, "main_prior": " x ".join(c["MAIN"]),
            "check_priors": [" x ".join(x) for x in c["CHECKS"]],
            "prior_source": f"{T30_SRC} DEPTH_PRIORS x CM_PRIORS (read from its source; T33 imports them unchanged)",
            "candidates": cands, "limits": T1_LIMITS,
            "one_day": {"source": f"{GRID_RUNS} and {GRID_HEDGE}", "task": head[0], "first": GRID_FIRST,
                        "candidate": GRID_CAND, "per_truth": gtab["per_task"], "by_prior": g_prior,
                        "summary_uniform": gtab["summary"],
                        "what": "the pooled base at the truth's own chamber mix (bown_cown, s1_base's analogue) "
                                "against the same at cardiomyocytes 0.60 (bown_c0.6, s1_cm60's analogue), over 15 "
                                "re-sampled mixes of one E10.5 truth; no library-3 or mid arm exists at one day, so "
                                "no one-day reading of the built pair exists"},
            "t31_context": {"source": T31, "what": "one-day backward gap, bases only, no cardiomyocyte arm: "
                                                   "lib3_base minus pooled_base", **seed_summary(d31),
                            "verdict": t31["verdict"]}}


# ========================================================================================================== heart
HEART_TASKS = "experiments/harness/results/trend_heart_v1/tasks"
HEART_FIRST, HEART_CAND = "b0_trend", "copy"
HEART_HEDGE = tuple(f"experiments/s2_hedge/results/per_seed/heart-scope{x}" for x in ("", "-e875", "-e825"))
REGISTRY = "results/registry.csv"
HEART_LIMITS = [
    "neither practice arm carried the uploaded files' size step or the half: the first upload had the pre-registered "
    "size (xy_grow_z_hold) and a uniform random half, the second has the raw-mode isotropic size and a half; the "
    "harness applies no size step (src/vec/harness.py _b0_t2), so D prices the expression difference only",
    "b0_trend was kept (row H11) on these same two backtests, so they are tuned-on for the first file",
    "the two backtests share one truth (E10.5): one unit",
    "after the fact, decision log (k), changing no number here: these two backtests did not predict the test for the "
    "expression step; its budget at the test was 0.5 day against 0.125 and 0.375 day on them, and it had been "
    "adopted on a margin of 0.06 on two backtests that share one truth",
]


def board_heart(rd: Reader, checks: list[str]) -> dict:
    tasks = {}
    for rel in rd.glob(HEART_TASKS, "*.json"):
        t = task_scores(rd, rel, checks)
        tasks[t["key"]] = t
    need(len(tasks) == 2, "heart: two backtests", checks)
    reg = [r for r in rd.csv(REGISTRY) if r["run_id"] == "trend_heart_v1"]
    for key, t in tasks.items():
        for arm in (HEART_FIRST, HEART_CAND):
            row = [r for r in reg if r["pseudo_task"] == key and r["candidate"] == arm]
            need(len(row) == 1 and abs(float(row[0]["score_mean"]) - statistics.fmean(t["arms"][arm])) <= 5e-4,
                 f"heart {key} {arm}: the registry's mean (written to three decimals)")
    checks.append("heart: both arms on both backtests give results/registry.csv's means (to its three decimals)")
    targets = {k: t["target"] for k, t in tasks.items()}
    tab = candidate_table({k: t["arms"][HEART_FIRST] for k, t in tasks.items()},
                          {k: t["arms"][HEART_CAND] for k, t in tasks.items()}, targets)
    tab.pop("d")
    hedge_arms = set()
    for d in HEART_HEDGE:
        for rel in rd.glob(d, "*/*.json"):
            hedge_arms |= set(rd.json(rel)["pred_values"])
    has = sorted(hedge_arms & {HEART_FIRST, HEART_CAND})
    return {"first": HEART_FIRST, "candidate": HEART_CAND, "targets": targets, "units": truth_units(targets),
            "builds": sorted({t["build"] for t in tasks.values()}), "table": tab, "limits": HEART_LIMITS,
            "s2_hedge_heart": {"arms": sorted(hedge_arms), "holds_first_or_candidate": has,
                               "what": "the heart-scope payoffs score scope variants scaled to an oracle size; they "
                                       "hold neither b0_trend nor copy, so they give no paired D"}}


# ============================================================================================================= T3
K7 = "experiments/t3_decomp/results/k7_wt_revert/tasks/{}.json"
K8 = "experiments/t3_decomp/results/k8_s2_form/tasks/{}.json"
K10 = "experiments/t3_decomp/results/k10_s1_form/tasks/{}.json"
K2A = "experiments/harness/results/t3_k2_v1/tasks/T3_Gata4_Mab21l2.json"
GATA4 = "experiments/t3_decomp/results/gata4/gata4_records.csv"
T3_FOLDS = ("g1", "g2", "m")
T3_KNOCKOUT = {"g1": "Gata4", "g2": "Gata4", "m": "Mab21l2"}  # g1 and g2: two replicates of one knockout, one slide
T3_CANDS = {  # name: ({fold: (source, arm)}, note)
    "refused form, WT-only factor per group (K8b; built file)": (
        {"g1": ("K8", "s2_form_revert_grp"), "g2": ("K8", "s2_form_revert_grp"), "m": ("K7", "wt_revert_grp")},
        "the built last file (t3_build/refused_s1.h5ad, ruling R5); on m the K7b form (no composition), as try K17 "
        "reads m; R5 was chosen on K7's and K8's readings of these three folds"),
    "refused form, pooled factor (K8a)": ({"g1": ("K8", "s2_form_revert"), "g2": ("K8", "s2_form_revert"),
                                           "m": ("K7", "wt_revert")}, "on m the K7a form"),
    "refused form, no factor (K8 base)": ({"g1": ("K8", "s2_form"), "g2": ("K8", "s2_form")},
                                          "the refused s1 as frozen before ruling R5; not scored on m"),
    "K7a wt_revert": ({f: ("K7", "wt_revert") for f in T3_FOLDS}, "wild-type cells, pooled factor, no composition"),
    "K7b wt_revert_grp": ({f: ("K7", "wt_revert_grp") for f in T3_FOLDS}, "wild-type cells, factor per group"),
    "K7 wt_draw": ({f: ("K7", "wt_draw") for f in T3_FOLDS}, "wild-type cells, no factor (K7's base)"),
    "K7c wt_revert_perm (control)": ({f: ("K7", "wt_revert_perm") for f in T3_FOLDS}, "K7a's factor permuted"),
}
T3_LIMITS = [
    "every candidate was read on these folds before ruling R5: K8's row says the factor's form is no longer blind to "
    "g1 and g2, and tries K16 and K17 read the exponent on all three",
    "m's reference is WT E9.5 (the K11 to K14 rows: reverts were flattered there); g1 and g2 are two replicates of one "
    "knockout on one slide",
    "after the fact, decision log (k), changing no number here: the same-knockout stand-in of the first file "
    "(stand-in A, xrep:F+depth) did not predict the test, while the cross-knockout readings did",
]


def _no_value(tab: dict) -> dict:
    """Against the floor, max(D, 0) is not the value of a second file: drop V from a T3 table."""
    units = {u: {k: v for k, v in x.items() if k != "seed_v"} for u, x in tab["units"].items()}
    summary = {k: v for k, v in tab["summary"].items() if k not in ("v_task", "v_seed")}
    return {**tab, "units": units, "summary": summary}


def board_t3(rd: Reader, checks: list[str]) -> dict:
    k7 = {f: task_scores(rd, K7.format(f), checks) for f in T3_FOLDS}
    k8 = {f: task_scores(rd, K8.format(f), checks) for f in ("g1", "g2")}
    for f in ("g1", "g2"):
        need(k7[f]["key"] == k8[f]["key"] and k7[f]["target"] == k8[f]["target"], f"T3 {f}: K7 and K8 one task", checks)
    for src in (*k7.values(), *k8.values()):
        need(all(x == 50.0 for x in src["arms"]["floor"]), f"T3 {src['key']}: the floor scores 50 in every seed")
    checks.append("T3: the floor is 50.0 in every seed on K7's and K8's folds")
    src = {"K7": k7, "K8": k8}
    targets = {f: k7[f]["target"] for f in T3_FOLDS}
    floor = {f: k7[f]["arms"]["floor"] for f in T3_FOLDS}
    cands = {}
    for name, (where, note) in T3_CANDS.items():
        sc = {f: src[s][f]["arms"][a] for f, (s, a) in where.items()}
        tab = candidate_table(floor, sc, targets)
        tab.pop("d")
        ko = unit_values({f: paired(floor[f], sc[f]) for f in sc},
                         truth_units({f: T3_KNOCKOUT[f] for f in sc}))
        cands[name] = {**_no_value(tab), "by_knockout": _no_value({"units": ko, "summary": summarise(ko)}),
                       "arms": {f: f"{s} {a}" for f, (s, a) in where.items()}, "note": note}
    # where the first file's stand-ins were scored, and why none pairs with a candidate
    rec = {}
    for r in rd.csv(GATA4):
        if r["t3"] != "":
            rec.setdefault((r["target"], r["variant"]), {})[int(r["seed"])] = float(r["t3"])
    stand = {}
    for rep, fold in (("rep1", "g1"), ("rep2", "g2")):
        base = rec[(rep, "_floor_wt_workcopy")]
        need(sorted(base) == SEEDS and all(v == 50.0 for v in base.values()), f"T3 {rep}: the records' floor is 50")
        for v in ("xrep:F+depth", "transfer:mab21l2_all"):
            stand[f"{v} on {rep} (gata4_records.csv)"] = seed_summary(paired([base[k] for k in SEEDS],
                                                                              [rec[(rep, v)][k] for k in SEEDS]))
        same = [rec[(rep, "blind:G2_bcat_central_d50")][k] - k8[fold]["arms"]["s2_form"][k] for k in SEEDS]
        stand[f"check {rep}: the frozen refused s1 in two runs, blind:G2_bcat_central_d50 (records, estimate.py "
              f"T3_S2) minus K8's s2_form, per seed"] = {
            "largest_abs": max(abs(x) for x in same), "mean": statistics.fmean(same)}
    for f in ("r12", "r21"):
        t = task_scores(rd, K10.format(f), checks)
        stand[f"s1x on {f} (K10; {t['key']})"] = seed_summary(paired(t["arms"]["floor"], t["arms"]["s1x"]))
    k2 = task_scores(rd, K2A, checks)
    need(k2["target"] == k7["m"]["target"] and k2["reference"] == k7["m"]["reference"],
         "T3: K2a's task has fold m's target and reference", checks)
    stand[f"t3_factor_xko on {k2['key']} (K2a)"] = seed_summary(paired(k2["arms"]["floor"],
                                                                      k2["arms"]["t3_factor_xko"]))
    xko = statistics.fmean(k2["arms"]["t3_factor_xko"])
    unpaired = {name: statistics.fmean(src[where["m"][0]]["m"]["arms"][where["m"][1]]) - xko
                for name, (where, _) in T3_CANDS.items() if "m" in where}
    return {"first": "allowed s1 (pooled Gata4 factor at exponent 1, the composition, Gata4's per-label depth ratios)",
            "paired_analogue": None,
            "why_none": "the allowed s1's stand-ins were scored on the replicate folds r12 and r21 (s1x, tries K10 "
                        "and K15), in the decomposition records on the Gata4 replicates (xrep:F+depth, another scoring "
                        "run: the same frozen file scored there and in K8 differs seed by seed, see the check lines) "
                        "and on T3/Gata4->Mab21l2 (K2a's t3_factor_xko, another run); no candidate was scored in any "
                        "of those runs, so no seed-paired D exists and no V is computed",
            "targets": targets, "candidates": cands, "first_standins_vs_floor": stand,
            "m_unpaired_cross_run": {"what": "on fold m's truth and reference only: the candidate's mean (K7) minus "
                                             "K2a's t3_factor_xko mean (t3_k2_v1), two runs, NOT seed-paired; no "
                                             "shared arm checks that the two runs score alike",
                                     "values": unpaired},
            "limits": T3_LIMITS}


# =========================================================================================================== text
def _sd(s: dict) -> str:
    return "  n/a" if s["sd_d"] is None else f"{s['sd_d']:5.2f}"


def _row(name: str, s: dict) -> str:
    return (f"  {name:<34} mean D {s['mean_d']:+6.2f}  SD {_sd(s)} ({s['sd_kind']}, {s['units']} units)  "
            f"V task {s['v_task']:5.2f}  V seed {s['v_seed']:5.2f}  share D>0 {s['share_positive']:.2f}")


def _tasks(tab: dict, order: list[str]) -> list[str]:
    return [f"    {t:<28} D {tab['per_task'][t]['mean']:+6.2f} (seed SD {tab['per_task'][t]['seed_sd']:4.2f}, ahead "
            f"{tab['per_task'][t]['ahead']} of {tab['per_task'][t]['seeds']})" for t in order if t in tab["per_task"]]


def text(res: dict) -> str:
    L = [res["label"], ""]
    e = res["embryo"]
    L += [f"EMBRYO. First file: {e['first']} (builds/T2embryo_s1.h5ad). Nine practice tasks, 10 seeds; scorer build "
          f"{', '.join(e['builds'])}; runs agree to {e['runs_largest_difference']:.4f} where an arm is in two.",
          "  units (one truth each): " + "; ".join(f"{u}: {', '.join(v)}" for u, v in e["units"].items()),
          f"  merges from: {e['merge_source']}"]
    for c, v in e["candidates"].items():
        L.append(_row(c, v["summary"]))
        L += _tasks(v, list(e["tasks"]))
        L.append("    E28's weighting S(p), unmerged: " + "; ".join(
            f"p {p}: V task {x['v_task']:.2f}, V seed {x['v_seed']:.2f}" for p, x in v["e28_weighting"].items()))
        L.append(f"    note: {v['note']}")
    L.append(f"  E28's weighting from: {e['e28_weighting_source']}")
    L += [f"  limit: {x}" for x in e["limits"]]
    t = res["T1"]
    L += ["", f"T1. First file: {t['first']} (builds/T1_s1.h5ad, row R6). Try T33's twelve zero-gap scenarios, 10 "
              f"seeds; priors: {t['prior_source']}; T33's main prior {t['main_prior']}, checks "
              f"{', '.join(t['check_priors'])}."]
    for c, v in t["candidates"].items():
        L.append(f"  {c} against {t['first']}: equal weights over the scenarios, then each stated prior")
        L.append(_row("    equal weights", v["summary"]))
        for p, s in v["by_prior"].items():
            L.append(_row(f"    {p}", s))
        L += _tasks(v, t["scenarios"])
        L.append(f"    note: {v['note']}")
    L.append("  stated priors, weight per scenario (" + ", ".join(t["scenarios"]) + "):")
    for p, w in t["priors"].items():
        L.append(f"    {p:<16} " + " ".join(f"{w[s]:.3f}" for s in t["scenarios"]))
    o = t["one_day"]
    L += [f"  one day ({o['task']}; {o['source']}): {o['what']}",
          _row(f"    {o['candidate']} vs {o['first']}, equal", o["summary_uniform"])]
    for p, s in o["by_prior"].items():
        L.append(_row(f"    {p}", s))
    L += _tasks({"per_task": o["per_truth"]}, list(o["per_truth"]))
    c31 = t["t31_context"]
    L.append(f"  context, try T31 ({c31['what']}): {c31['mean']:+.2f} (seed SD {c31['seed_sd']:.2f}, ahead "
             f"{c31['ahead']} of {c31['seeds']}); its verdict: {c31['verdict']}")
    L += [f"  limit: {x}" for x in t["limits"]]
    h = res["heart"]
    L += ["", f"HEART. First file: {h['first']} (the uploaded half of builds/T2heart_s1.h5ad). Candidate: "
              f"{h['candidate']} (the built half of the heart2 set). Two counted backtests, 10 seeds; scorer build "
              f"{', '.join(h['builds'])}.",
          "  units: " + "; ".join(f"{u}: {', '.join(v)}" for u, v in h["units"].items()),
          _row(h["candidate"], h["table"]["summary"])]
    L += _tasks(h["table"], list(h["targets"]))
    L.append(f"  s2_hedge heart-scope arms {h['s2_hedge_heart']['arms']}: {h['s2_hedge_heart']['what']}")
    L += [f"  limit: {x}" for x in h["limits"]]
    r = res["T3"]
    L += ["", f"T3. First file: {r['first']}. NO PAIRED ANALOGUE: {r['why_none']}.",
          "  The candidates against the floor (wt_identity, 50 in every seed), NOT against the banked file; "
          "three folds, one truth each (" + ", ".join(f"{f}: {x}" for f, x in r["targets"].items()) + "):"]
    for c, v in r["candidates"].items():
        s = v["summary"]
        L.append(f"  {c:<56} over floor {s['mean_d']:+6.2f}  SD {_sd(s)} ({s['units']} folds)  share > 0 "
                 f"{s['share_positive']:.2f}; by knockout {v['by_knockout']['summary']['mean_d']:+6.2f} "
                 f"({v['by_knockout']['summary']['units']} units)")
        L += _tasks(v, list(T3_FOLDS))
        L.append(f"    arms: {v['arms']}; note: {v['note']}")
    L.append("  the first file's stand-ins against their own floors, on their own folds (no candidate there):")
    for k, s in r["first_standins_vs_floor"].items():
        if "seeds" in s:
            L.append(f"    {k:<62} {s['mean']:+6.2f} (seed SD {s['seed_sd']:.2f}, ahead {s['ahead']} of "
                     f"{s['seeds']})")
        else:
            L.append(f"    {k}: largest {s['largest_abs']:.2f}, mean {s['mean']:+.2f}")
    u = r["m_unpaired_cross_run"]
    L.append(f"  {u['what']}:")
    L += [f"    {k:<56} {v:+6.2f}" for k, v in u["values"].items()]
    L += [f"  limit: {x}" for x in r["limits"]]
    L += ["", "checks passed: " + str(len(res["checks"])), "inputs (sha256):"]
    L += [f"  {k} {v}" for k, v in res["inputs_sha256"].items()]
    return "\n".join(L) + "\n"


# =========================================================================================================== main
LABEL = ("Second-file value table, 2026-11-18: V = E[max(D, 0)], D = candidate minus the practice analogue of the "
         "uploaded first file, seed-paired on committed practice tasks. Aggregation of committed results, not a try; "
         "no test score, leaderboard figure or made-up prior enters any number.")


def build(rd: Reader | None = None) -> dict:
    rd = rd or Reader()
    checks: list[str] = []
    res = {"label": LABEL, "embryo": board_embryo(rd, checks), "T1": board_t1(rd, checks),
           "heart": board_heart(rd, checks), "T3": board_t3(rd, checks)}
    res["checks"] = sorted(set(checks))
    res["inputs_sha256"] = dict(sorted(rd.read.items()))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", metavar="PATH")
    ap.add_argument("--txt", metavar="PATH")
    a = ap.parse_args(argv)
    res = build()
    out = text(res)
    print(out, end="")
    if a.json:
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.json).write_text(json.dumps(res, indent=1) + "\n")
    if a.txt:
        Path(a.txt).parent.mkdir(parents=True, exist_ok=True)
        Path(a.txt).write_text(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
