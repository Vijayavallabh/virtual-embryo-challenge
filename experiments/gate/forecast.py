#!/usr/bin/env python3
"""A forecast of our leaderboard total as a distribution, for the upload files of 2026-11-05 and of 2026-11-13.

Aggregation of committed practice results only, not a try: no target, no data file, no scorer, no build, no new score.
Every number is a mean, maximum, difference or percentile of board scores already committed under experiments/, read
through experiments/gate/estimate.py's door (which refuses boards/, panels/, docs/research/, data/, runs/ and any
index.json), with that script's loaders and outcome rules imported, not copied. No test score and no leaderboard figure
enters any world, weight, score distribution or choice. The comparison lines (LINES below) are printed against the
finished distribution only.

Total = T1 + T3 + (heart + embryo) / 2, each board at the better of its two uploads in the same world and seed
(estimate.py, reading A). A Monte Carlo draws one outcome per board, the boards independent (estimate.py's assumption;
no committed result ties two boards' hidden facts together), and reads percentiles and P(total >= c) off the draws.

File sets:
  2026-11-05  the files the committed upper reading prices (191.51; results/upload_gate.json estimate_reproducible):
              T1 s1 pooled base, s2 cardiomyocytes 0.60; heart s1 hold, s2 half_radial; embryo s1 fgw, s2 fgw_shape0;
              T3 first upload the refused branch's s1 without factor, second the allowed branch's s1.
  2026-11-13  rows R4, R5, R6 (experiments/harness/TRIES.md; docs/PLAN.md decision log 2026-11-13 (d)): embryo s2
              fgw_aligned_series; T3's first upload carries the wild-type-only factor per lineage group; T1 s1 on
              E10.5's library 3. Heart and the other files unchanged.

Usage, from the repository root:
  python3 -I -B experiments/gate/forecast.py [--files 2026-11-05|2026-11-13|both] [--json PATH]
  uv run --frozen --no-sync python -I -B experiments/gate/forecast.py --files both \
      --json experiments/gate/results/forecast_2026-11-13.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import estimate as E  # noqa: E402

# ==================================================================================================== THE WORLD TABLE
# Every world, weight and per-world score distribution. "OLD" = stated by estimate.py (the committed upper reading's
# procedure) and reused as is. "SOURCED" = a committed result file and the ledger row it comes from. "JUDGEMENT" = no
# committed result fixes it; low and high are the values the sensitivity table sets it to, one at a time.
#
# T1, all files. World: the oracle-base grid's truth, E10.5 at blood own/0.05/0.4 x cardiomyocytes
#   0.15/own/0.45/0.60/0.75 (15 truths). Weight: OLD, uniform 1/15 (estimate.PRIOR; the hedge file's four other stated
#   priors in the sensitivity table). Scores: OLD, the grid per seed (s2_hedge per_seed/t1-grid; s1 bown_cown, s2
#   bown_c0.6) plus M1's step over B0 (+1.29, t1_m1 diag_v1).
# T1, new files (R6: s1 on library 3). World: E12.5's depth relative to library 3, 0.35 or 0.68 (shallow), 1.0 or 1.37
#   (deep), independent of the grid truth. Weight: p_deep = P(1.0 or 1.37), JUDGEMENT 5/7 (low 1/3, high 6/7); 5/7 is
#   T30's census prior (its main), 1/3 its latest prior; within a side T30 census's split (1:1 shallow, 1:4 deep).
#   Scores: s1 moved by k x (lib3_base - s1_base), per seed, from T30's scenario at that depth and the truth's
#   cardiomyocyte share (t1_scope/results/t30_pricing/scenarios; TRIES T30); s2 unchanged. k = JUDGEMENT 0.2 (low 0.10,
#   high 0.33): T30 is at zero time gap, and R6 expects "about a tenth to a third of those sizes" at a real gap.
# heart, all files. World: transition E9.5 or E8.75 -> E10.5 x five scaled truths. Weight: OLD, transitions even,
#   uniform over truths (the four other stated priors in the sensitivity table). Scores: OLD, heart-scope payoffs moved
#   to the entry's d2_v1 level.
# embryo, old files. World: the two counted backtests (E7.5, E7.25). Weight: OLD, even. Scores: OLD, emb_shape_v1.
# embryo, new files (R4: s2 fgw_aligned_series). Worlds: D, E7.75 shares its deeper bracket's profile (measured like
#   E8.0; tasks half, late, new1, x1, x2, x3); O, it does not (new2, E7.5, E7.25). Weight: p = P(D), JUDGEMENT 0.45
#   (low 0.30, high 0.60): decision log 2026-11-13 (a) "0.4 to 0.5, a judgement"; TRIES E28 reads S at 0.3/0.45/0.6.
#   Scores: today's pair on the counted backtests (OLD level) plus the new pair's per-seed change against today's pair
#   on each task of the world (TRIES E24, E25, E28; harness/results/e24_pricing, e25_pairs, e28_nine_tasks), never
#   below the unchanged first upload in the same seed.
# T3, all files. World: the fold s1's bet meets: A (Gata4, other replicate), B (Mab21l2 factor on Gata4), K2a (Gata4
#   factor on Mab21l2). Weight: OLD, A 0.5, B 0.25, K2a 0.25 (estimate.T3_WEIGHTS_K2A); JUDGEMENT low 0.25, high 0.75
#   for A, the rest split evenly. Scores: OLD, t3_decomp/results/gata4/gata4_records.csv and t3_k2_v1.
# T3, new files (R5: the first upload carries the WT-only factor per lineage group). Worlds: carry, the factor's gain
#   carries to beta-catenin as on Gata4 (K8b); vanish, it acts like an unrelated factor (K7c, the permuted control).
#   Weight: q = P(carry), JUDGEMENT 0.6 (low 0.35, high 0.85): 0.6 is K8's own pre-run chance, 0.35 K7's. Scores: the
#   first upload's stand-in plus K8b's per-seed gain over s2_form (k8_s2_form) or K7c's over wt_draw (k7_wt_revert);
#   on K2a's fold, our WT copy plus the same gain averaged over the two replicates.
#
# Boards are independent (estimate.py). Within a board the new worlds are independent of the old ones (T1 depth of grid
# truth; T3 carry of the fold).
JUDGEMENT = {  # name: (central, low, high, what it is, file sets it applies to)
    "t3_standin_A": (0.5, 0.25, 0.75, "T3: weight of stand-in A (central = estimate.T3_WEIGHTS_K2A)", ("old", "new")),
    "emb_p_deep": (0.45, 0.30, 0.60, "embryo: P(E7.75 shares its deeper bracket's profile)", ("new",)),
    "t3_q_carry": (0.6, 0.35, 0.85, "T3: P(the WT-only factor's Gata4 gain carries to beta-catenin)", ("new",)),
    "t1_p_deep": (5 / 7, 1 / 3, 6 / 7, "T1: P(E12.5 at least as deep as library 3)", ("new",)),
    "t1_k": (0.2, 0.10, 0.33, "T1: real-gap size of T30's zero-gap base effect", ("new",)),
}
STATED_PRIORS = ("t1_prior", "heart_prior")  # OLD weights with stated alternatives (each hedge file's priors)

# T1: T30's cardiomyocyte shares for the grid's truths (nearest T30 level; 0.15 and 0.75 lie outside T30's range)
T1_SHARE = {"c0.15": "own", "cown": "own", "c0.45": "0.45", "c0.6": "0.60", "c0.75": "0.60"}
T30_DEPTHS = (0.35, 0.68, 1.0, 1.37)
T30_DEEP_SPLIT = {1.0: 0.2, 1.37: 0.8}  # T30 census: 1/7 and 4/7 of the deep side
T30_SHALLOW_SPLIT = {0.35: 0.5, 0.68: 0.5}  # T30 census: 1/7 and 1/7
T30_DIR = "experiments/t1_scope/results/t30_pricing"
T30_SUMMARY = f"{T30_DIR}/pricing.json"

# embryo: the nine practice tasks (TRIES E28; experiments/harness/nine_tasks.py DEEP and OTHER)
EMB_RUNS = ("experiments/harness/results/e24_pricing/tasks", "experiments/harness/results/e25_pairs/tasks",
            "experiments/harness/results/e28_nine_tasks/tasks")
EMB_SUMMARY = "experiments/harness/results/e28_nine_tasks/nine_tasks.json"
EMB_DEEP = ("half", "late", "new1", "x1", "x2", "x3")
EMB_OTHER = ("new2", "e75", "e725")
EMB_NEW_S2 = "fgw_aligned_series"

# T3: the first upload's factor (TRIES K7, K8, R5)
K8 = "experiments/t3_decomp/results/k8_s2_form/tasks/{}.json"
K7 = "experiments/t3_decomp/results/k7_wt_revert/tasks/{}.json"
T3_REP_TASK = {"rep1": "g1", "rep2": "g2"}
K8_ROW_GAIN = {"rep1": 4.72, "rep2": 7.12}  # TRIES K8, K8b against s2_form; checked against the files

# the committed upper reading, for the self-check of the old-files run
UPPER_JSON = "experiments/gate/results/estimate_2026-11-05.json"
PATH_TO_TARGET_MEAN = 191.49689560437218  # experiments/gate/results/path_to_target_2026-11-06.json mean_overall
REPRO_TOL = 0.5

# the forecast's draws
N_DRAWS = 1_000_000
SEED = E.RNG_SEED
QS = (5, 25, 50, 75, 90, 95)

# ================================================================= comparison lines: output only, never an input
# 203.00, 205.32, 215.41: the sixth, third and first Human-track overall totals of a screened leaderboard read,
# 2026-11-13 07:52 UTC (overall totals only, calibration only; the sixth is the last of six prize places).
# 209: the captain's fixed gate (results/upload_gate.json). 197 and 200: round lines.
LINES = {197.0: "round line", 200.0: "round line", 203.00: "6th Human total, screened read 2026-11-13 (last prize)",
         205.32: "3rd Human total, screened read 2026-11-13", 209.0: "the captain's gate",
         215.41: "1st Human total, screened read 2026-11-13"}
TAIL_LINES = (203.00, 209.0)

E.READABLE = (*E.READABLE, f"{T30_DIR}/", "experiments/t3_decomp/results/k7_wt_revert/",
              "experiments/t3_decomp/results/k8_s2_form/", UPPER_JSON)


def need(cond: bool, what: str, checks: list[str] | None = None) -> None:
    """A self-check: assert, and record it."""
    if not cond:
        raise AssertionError(f"forecast: check failed: {what}")
    if checks is not None:
        checks.append(what)


# ============================================================================================== labelled outcomes
def mix(parts: list[tuple[float, dict, str]]) -> dict:
    """estimate.mixture with a world label per outcome: parts are (weight, outcome set, label)."""
    parts = [p for p in parts if p[0] > 0]
    tot = sum(p for p, _, _ in parts)
    out = E.mixture([(p, o) for p, o, _ in parts])
    out["lab"] = np.concatenate([np.full(o["a"].size, lab, dtype=object) if "lab" not in o
                                 else np.array([f"{lab}|{x}" for x in o["lab"]], dtype=object) for _, o, lab in parts])
    need(abs(sum(p / tot for p, _, _ in parts) - 1) < 1e-12 and abs(out["w"].sum() - 1) < 1e-9, "weights sum to 1")
    return out


def same_dist(x: dict, y: dict) -> bool:
    """Two outcome sets hold the same weighted values (order free)."""
    def key(o, k):
        i = np.lexsort((o["w"], o[k]))
        return o[k][i], o["w"][i]
    return all(np.allclose(key(x, k)[0], key(y, k)[0]) and np.allclose(key(x, k)[1], key(y, k)[1]) for k in "ab")


# ======================================================================================================== boards
class Inputs:
    """Everything read, once, through estimate.py's door."""

    def __init__(self):
        self.door, self.checks = E.Door(REPO), []
        d, c = self.door, self.checks
        self.t1_old, self.t1_priors, self.t1_steps = E.board_t1(d, c)
        self.heart_old, self.heart_priors, _ = E.board_heart(d, c)
        self.embryo_old = E.board_embryo(d, c)
        _, self.t3_old_est, _ = E.board_t3(d, c)
        self.grid = E.load_runs(d, E.T1_GRID_RUNS)
        self.t1_hedge_priors = d.json(E.T1_GRID_HEDGE)["priors"]
        self._t30()
        self._embryo()
        self._t3()

    # ---- T1: T30's per-seed base effect, lib3_base - s1_base, by depth and share
    def _t30(self):
        summ = self.door.json(T30_SUMMARY)["arm_by_scenario"]
        self.t30 = {}
        for dep in T30_DEPTHS:
            for sh in ("own", "0.45", "0.60"):
                rel = f"{T30_DIR}/scenarios/t30_d{dep}_cm{sh}.json"
                raw = self.door.json(rel)
                rows = sorted(raw["rows"], key=lambda r: r["seed"])
                seeds = [int(r["seed"]) for r in rows]
                s1 = np.array([float(r["scores"]["s1_base"]) for r in rows])
                l3 = np.array([float(r["scores"]["lib3_base"]) for r in rows])
                need(seeds == self.grid["seeds"], f"T30 d{dep} cm{sh}: seeds are the T1 grid's (paired by seed)")
                scen = f"d{dep} cm{sh}"
                E.same([s1.mean(), l3.mean()], [summ["s1_base"][scen], summ["lib3_base"][scen]],
                       f"T30 {scen}: per-seed files give pricing.json's arm means", self.checks)
                self.t30[dep, sh] = l3 - s1

    # ---- embryo: per-seed change of the new pair against today's on each of the nine tasks
    def _embryo(self):
        scores: dict[str, dict[str, np.ndarray]] = {}
        for d in EMB_RUNS:
            for p in sorted((REPO / d).glob("T2e_*.json")):
                rel = p.relative_to(REPO).as_posix()
                t = self.door.json(rel)
                rows = sorted(t["rows"], key=lambda r: r["seed"])
                key = t["task"]["key"]
                for arm in ("fgw", "fgw_shape0", EMB_NEW_S2):
                    if arm not in rows[0]["scores"]:  # each run holds some of the arms (nine_tasks.py merges them)
                        continue
                    v =np.array([float(r["scores"][arm]) for r in rows])
                    old = scores.setdefault(key, {}).get(arm)
                    need(old is None or float(np.abs(old - v).max()) <= 0.01, f"embryo {key} {arm}: runs agree")
                    scores[key][arm] = v
        summ = self.door.json(EMB_SUMMARY)
        fgw_means, pair = summ["single"]["fgw"], summ["values"][f"fgw|{EMB_NEW_S2}"]["tasks"]
        self.emb_delta = {}
        for key, s in scores.items():
            names = [n for n, m in fgw_means.items() if abs(m - s["fgw"].mean()) < 1e-9]
            need(len(names) == 1, f"embryo {key}: one E28 task name by fgw's mean")
            delta = np.maximum(s["fgw"], s[EMB_NEW_S2]) - np.maximum(s["fgw"], s["fgw_shape0"])
            E.same(delta.mean(), pair[names[0]], f"embryo {names[0]}: per-seed change gives E28's value", self.checks)
            self.emb_delta[names[0]] = delta
        need(sorted(self.emb_delta) == sorted(EMB_DEEP + EMB_OTHER), "embryo: all nine E28 tasks found", self.checks)
        # today's pair on the counted backtests, per backtest and seed (estimate.py's own files)
        self.emb_level = []
        for shape_rel, _ in E.EMBRYO:
            s = E.load_task(self.door, shape_rel, (E.EMBRYO_S1, E.EMBRYO_S2))
            self.emb_level.append((np.maximum(s[E.EMBRYO_S1], s[E.EMBRYO_S2]), s[E.EMBRYO_S1]))
            need(len(s["seeds"]) == len(next(iter(self.emb_delta.values()))), "embryo: counted backtests' seeds")

    # ---- T3: stand-ins, and the first upload's factor gain per replicate
    def _t3(self):
        rec = E.load_t3(self.door, (E.T3_STANDINS["A"], E.T3_STANDINS["B"], E.T3_S2), self.checks)
        k2 = E.load_task(self.door, E.T3_K2, ("floor", E.T3_K2A, E.T3_K2_S2))
        self.t3_standins = {n: {r: rec[r][v] for r in E.T3_REPLICATES} for n, v in E.T3_STANDINS.items()}
        self.t3_k2a, self.t3_k2_copy = k2[E.T3_K2A], k2[E.T3_K2_S2]
        self.t3_s2 = {r: rec[r][E.T3_S2] for r in E.T3_REPLICATES}
        self.t3_gain = {"carry": {}, "vanish": {}}
        for r, g in T3_REP_TASK.items():
            k8 = E.load_task(self.door, K8.format(g), ("s2_form", "s2_form_revert_grp"))
            k7 = E.load_task(self.door, K7.format(g), ("wt_draw", "wt_revert_perm"))
            need(k8["seeds"] == k7["seeds"] and len(k8["seeds"]) == self.t3_s2[r].size, f"T3 {r}: K7/K8 seeds")
            self.t3_gain["carry"][r] = k8["s2_form_revert_grp"] - k8["s2_form"]
            self.t3_gain["vanish"][r] = k7["wt_revert_perm"] - k7["wt_draw"]
            need(abs(self.t3_gain["carry"][r].mean() - K8_ROW_GAIN[r]) < 0.01,
                 f"T3 {r}: K8b's gain over s2_form is the row's {K8_ROW_GAIN[r]}", self.checks)


def board_t1(inp: Inputs, files: str, prior_name: str, p_deep: float, k: float) -> dict:
    sc, prior, step = inp.grid["scores"], inp.t1_hedge_priors[prior_name], inp.t1_steps["m1_minus_b0"]

    def by_truth(shift_s1) -> list:
        parts = []
        for t, pt in prior.items():
            s1 = sc[t][E.T1_S1] + shift_s1(t)
            o = E.grid_outcomes({t: s1}, {t: sc[t][E.T1_S2]}, {t: 1.0}, shift=step)
            parts.append((pt, o, f"truth {t}"))
        return parts

    if files == "old":
        return mix(by_truth(lambda t: 0.0))
    depth_w = {**{d: (1 - p_deep) * f for d, f in T30_SHALLOW_SPLIT.items()},
               **{d: p_deep * f for d, f in T30_DEEP_SPLIT.items()}}
    return mix([(wd, mix(by_truth(lambda t, d=d: k * inp.t30[d, T1_SHARE[t.split("_")[1]]])), f"depth {d}")
                for d, wd in depth_w.items()])


def board_heart(inp: Inputs, prior_name: str) -> dict:
    o = dict(inp.heart_priors[prior_name])
    o["lab"] = np.full(o["a"].size, "heart (old procedure)", dtype=object)
    return o


def board_embryo(inp: Inputs, files: str, p: float) -> dict:
    if files == "old":
        o = dict(inp.embryo_old)
        o["lab"] = np.full(o["a"].size, "counted backtests", dtype=object)
        return o

    def world(tasks) -> dict:
        a, b = [], []
        for lev_a, lev_b in inp.emb_level:
            for n in tasks:
                # the first upload (fgw) is unchanged, so the pair never falls below it in the same seed
                a.append(np.maximum(lev_b, lev_a + inp.emb_delta[n]))
                b.append(lev_b)
        a, b = np.concatenate(a), np.concatenate(b)
        return {"a": a, "b": b, "w": np.full(a.size, 1.0 / a.size)}
    return mix([(p, world(EMB_DEEP), "D: like its deeper bracket"), (1 - p, world(EMB_OTHER), "O: other")])


def board_t3(inp: Inputs, files: str, w_a: float, q: float) -> dict:
    weights = {"A": w_a, "B": (1 - w_a) / 2, "K2a": (1 - w_a) / 2}

    def by_standin(gain: dict | None) -> dict:
        s2 = {r: inp.t3_s2[r] + (gain[r] if gain else 0.0) for r in inp.t3_s2}
        copy = inp.t3_k2_copy + (np.mean([gain[r] for r in gain], axis=0) if gain else 0.0)
        standins = {**inp.t3_standins, "K2a": (inp.t3_k2a, copy)}
        return mix([(wt, E.t3_outcomes(standins, s2, {n: 1.0}), f"stand-in {n}") for n, wt in weights.items()])

    if files == "old":
        return by_standin(None)
    return mix([(q, by_standin(inp.t3_gain["carry"]), "factor carries"),
                (1 - q, by_standin(inp.t3_gain["vanish"]), "factor vanishes")])


def boards(inp: Inputs, files: str, par: dict) -> dict:
    out = {"T1": board_t1(inp, files, par["t1_prior"], par["t1_p_deep"], par["t1_k"]),
           "heart": board_heart(inp, par["heart_prior"]),
           "embryo": board_embryo(inp, files, par["emb_p_deep"]),
           "T3": board_t3(inp, files, par["t3_standin_A"], par["t3_q_carry"])}
    for b, o in out.items():
        need(abs(o["w"].sum() - 1) < 1e-9, f"{b}: weights sum to 1")
        inside = (o["a"] >= 0) & (o["a"] <= 100) & (o["b"] >= 0) & (o["b"] <= 100)
        need(bool(np.all(inside)), f"{b}: scores in [0, 100]")
        need(bool(np.all(o["a"] >= o["b"] - 1e-12)), f"{b}: the better of two is at least the first upload")
    return out


# ===================================================================================================== the draws
def draw(bd: dict, n: int = N_DRAWS, seed: int = SEED) -> tuple[np.ndarray, dict]:
    """One outcome per board per draw, the boards independent (estimate.bands' order and generator)."""
    rng = np.random.default_rng(seed)
    idx = {b: rng.choice(bd[b]["w"].size, size=n, p=bd[b]["w"]) for b in E.BOARDS}
    tot = sum(E.OVERALL_WEIGHT[b] * bd[b]["a"][idx[b]] for b in E.BOARDS)
    return tot, idx


def summary(bd: dict, tot: np.ndarray) -> dict:
    return {"mean_exact": E.overall({b: E.reading(o, "a") for b, o in bd.items()}),
            "mean_draws": float(tot.mean()), "median": float(np.percentile(tot, 50)),
            "percentiles": {str(q): float(np.percentile(tot, q)) for q in QS},
            "board_mean_better_of_two": {b: E.reading(o, "a") for b, o in bd.items()},
            "board_mean_first_upload": {b: E.reading(o, "b") for b, o in bd.items()},
            "total_first_round_mean": E.overall({b: E.reading(o, "b") for b, o in bd.items()})}


def p_at_least(tot: np.ndarray) -> dict:
    return {f"{c:.2f}": float(np.mean(tot >= c)) for c in LINES}


def tail(bd: dict, tot: np.ndarray, idx: dict) -> dict:
    """Per board, each world's weight and its share of the draws at or above each tail line."""
    out = {}
    for c in TAIL_LINES:
        hit = tot >= c
        per = {}
        for b in E.BOARDS:
            labs = bd[b]["lab"]
            coarse = {}
            for lab in sorted(set(labs.tolist())):  # each world dimension on its own ("|" separates dimensions)
                m = labs == lab
                wt = float(bd[b]["w"][m].sum())
                sh = float(np.mean(m[idx[b][hit]])) if hit.any() else 0.0
                for key in lab.split("|"):
                    coarse.setdefault(key, [0.0, 0.0])
                    coarse[key][0] += wt
                    coarse[key][1] += sh
            per[b] = {"mean_in_tail": float(bd[b]["a"][idx[b][hit]].mean()) if hit.any() else None,
                      "worlds": {k: {"weight": v[0], "share_of_tail": v[1]} for k, v in coarse.items()}}
        out[f"{c:.2f}"] = {"draws": int(hit.sum()), "boards": per}
    return out


# ===================================================================================================== forecast
def central() -> dict:
    return {**{k: v[0] for k, v in JUDGEMENT.items()}, "t1_prior": E.PRIOR, "heart_prior": E.PRIOR}


def run_set(inp: Inputs, files: str) -> dict:
    par = central()
    bd = boards(inp, files, par)
    tot, idx = draw(bd)
    res = {"files": files, "parameters": par, **summary(bd, tot), "p_at_least": p_at_least(tot),
           "tail": tail(bd, tot, idx)}
    sens = {}
    for name, (c, lo, hi, what, applies) in JUDGEMENT.items():
        if ("old" if files == "old" else "new") not in applies:
            continue
        row = {"what": what, "central": c}
        for side, v in (("low", lo), ("high", hi)):
            t, _ = draw(boards(inp, files, {**par, name: v}))
            row[side] = {"value": v, "mean": float(t.mean()), "p_203": float(np.mean(t >= 203.0)),
                         "p_209": float(np.mean(t >= 209.0))}
        sens[name] = row
    for name, priors in (("t1_prior", inp.t1_hedge_priors), ("heart_prior", inp.heart_priors)):
        row = {"what": f"{name.split('_')[0]}: the hedge file's stated priors (OLD: {E.PRIOR})", "central": E.PRIOR}
        for alt in priors:
            if alt == E.PRIOR:
                continue
            t, _ = draw(boards(inp, files, {**par, name: alt}))
            row[alt] = {"mean": float(t.mean()), "p_203": float(np.mean(t >= 203.0)),
                        "p_209": float(np.mean(t >= 209.0))}
        sens[name] = row
    res["sensitivity"] = sens
    return res, bd, tot


def forecast(sets: tuple[str, ...]) -> dict:
    inp = Inputs()
    checks = inp.checks
    # the old procedure, reproduced: rebuilt old boards equal estimate.py's, and the committed readings hold
    old_bd = boards(inp, "old", central())
    est = {"T1": inp.t1_old, "heart": inp.heart_old, "embryo": inp.embryo_old, "T3": inp.t3_old_est}
    for b in E.BOARDS:
        need(same_dist(old_bd[b], est[b]), f"old files, {b}: rebuilt outcomes equal estimate.py's", checks)
    upper = inp.door.json(UPPER_JSON)["readings"]["A"]["Overall"]
    gate_upper = inp.door.json(E.GATE)["estimate_reproducible"]["better_of_two"]
    t200, _ = draw(est, E.N_DRAWS, E.RNG_SEED)
    need(abs(float(t200.mean()) - PATH_TO_TARGET_MEAN) < 1e-9,
         "old files: estimate.py's boards at 200,000 draws give path_to_target's committed mean exactly", checks)
    out = {"label": "forecast of our leaderboard total from committed practice results only: aggregation, not a try; "
                    "no test score and no leaderboard figure enters any world, weight or distribution",
           "total": "T1 + T3 + (heart + embryo) / 2, each board at the better of its two uploads",
           "correlation": "boards independent (estimate.py's assumption; no committed result links two boards' "
                          "hidden facts); within a board the new worlds independent of the old ones",
           "judgement_parameters": {k: {"central": v[0], "low": v[1], "high": v[2], "what": v[3]}
                                    for k, v in JUDGEMENT.items()},
           "comparison_lines": {f"{c:.2f}": v for c, v in LINES.items()},
           "draws": N_DRAWS, "rng_seed": SEED, "sets": {}}
    for s in sets:
        files = "old" if s == "2026-11-05" else "new"
        res, bd, tot = run_set(inp, files)
        if files == "old":
            for what, ref in (("estimate_2026-11-05.json readings.A.Overall", upper),
                              ("upload_gate.json estimate_reproducible.better_of_two", gate_upper)):
                need(abs(res["mean_exact"] - ref) <= REPRO_TOL and abs(res["mean_draws"] - ref) <= REPRO_TOL,
                     f"old files: the mean matches the committed upper reading ({what} {ref}) within {REPRO_TOL}",
                     checks)
        out["sets"][s] = res
    if "2026-11-13" in out["sets"]:
        new = out["sets"]["2026-11-13"]
        need(new["board_mean_better_of_two"]["heart"] == E.reading(inp.heart_old, "a"), "heart unchanged", checks)
        # what moved: the old files with one board's new file put in, one board at a time
        new_bd = boards(inp, "new", central())
        out["one_board_at_a_time"] = {}
        for b in ("T1", "embryo", "T3"):
            t, _ = draw({**old_bd, b: new_bd[b]})
            out["one_board_at_a_time"][b] = {
                "board_mean_change": E.reading(new_bd[b], "a") - E.reading(old_bd[b], "a"),
                "total_mean": float(t.mean()), "p_203": float(np.mean(t >= 203.0)), "p_209": float(np.mean(t >= 209.0))}
    out["checks"] = sorted(set(checks))
    hashes = dict(inp.door.read)
    for f in ("experiments/gate/estimate.py", "experiments/gate/path_to_target.py", "experiments/gate/forecast.py"):
        hashes[f] = hashlib.sha256((REPO / f).read_bytes()).hexdigest()
    out["inputs_sha256"] = dict(sorted(hashes.items()))
    return out


# ======================================================================================================= output
def text(res: dict) -> str:
    L = []
    for s, r in res["sets"].items():
        L.append(f"== files of {s} ==  mean {r['mean_exact']:.2f} (draws {r['mean_draws']:.2f}), median "
                 f"{r['median']:.2f}; percentiles " + ", ".join(f"{q}th {v:.1f}" for q, v in r["percentiles"].items()))
        L.append("   board means, better of two / first upload alone: " + ", ".join(
            f"{b} {r['board_mean_better_of_two'][b]:.2f}/{r['board_mean_first_upload'][b]:.2f}" for b in E.BOARDS)
            + f"; first round total {r['total_first_round_mean']:.2f}")
        L.append("   P(total >= c): " + ", ".join(f"{c} {100 * p:.2f}%" for c, p in r["p_at_least"].items()))
        for c, t in r["tail"].items():
            L.append(f"   tail >= {c} ({t['draws']} draws): " + "; ".join(
                f"{b} (mean {d['mean_in_tail']:.1f}) " + ", ".join(
                    f"{k} w{v['weight']:.2f}->{v['share_of_tail']:.2f}" for k, v in d["worlds"].items()
                    if v["weight"] < 1 and not k.startswith("truth"))
                for b, d in t["boards"].items()))
        L.append("   sensitivity (P>=203 / P>=209, %):")
        for n, row in r["sensitivity"].items():
            cells = [f"{k} {('' if 'value' not in v else format(v['value'], '.3g') + ': ')}"
                     f"{100 * v['p_203']:.2f}/{100 * v['p_209']:.2f}" for k, v in row.items() if isinstance(v, dict)]
            L.append(f"     {n:13s} " + "; ".join(cells))
    for b, v in res.get("one_board_at_a_time", {}).items():
        L.append(f"old files with only {b} new: board {v['board_mean_change']:+.2f}, total mean {v['total_mean']:.2f}, "
                 f"P>=203 {100 * v['p_203']:.2f}%, P>=209 {100 * v['p_209']:.2f}%")
    L.append(f"{len(res['checks'])} self-checks passed; {len(res['inputs_sha256'])} files hashed; {res['correlation']}")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--files", choices=("2026-11-05", "2026-11-13", "both"), default="both")
    ap.add_argument("--json", metavar="PATH")
    a = ap.parse_args(argv)
    sets = ("2026-11-05", "2026-11-13") if a.files == "both" else (a.files,)
    res = forecast(sets)
    again = forecast(sets)
    need(json.dumps(res, sort_keys=True) == json.dumps(again, sort_keys=True), "a fixed seed gives identical output")
    res["checks"].append("a fixed seed gives identical output (the whole forecast run twice)")
    print(text(res))
    if a.json:
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.json).write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
