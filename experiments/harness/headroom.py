#!/usr/bin/env python3
"""Baseline headroom per test board: what each entry gains over the published baselines, and where most is left.

Principle served: the captain's guiding principle (2026-11-04), in particular "Begin with the published copy-last-stage
and unchanged-wild-type baselines across the boards you intend to enter. Improve the weakest weighted component before
polishing what already works", and its question for every component: "What meaningful change beyond the baseline does
this predict, and what held-out evidence shows it?" This script answers the first half of that question with numbers
already committed; the second half is the label printed with every row.

For each test board, for the current s1 entry (and the s2 where the harness scored it), on each practice task the
adoption rule counts, it prints for every ranked metric:
- its weight in the board score (docs/CHALLENGE.md "Metric panel"; `GROUPS` below has the group ids and weights of
  src/vec/boards.py GROUPS). On T2 the three shape terms are multiplied by the structure gate min(1, 2 x
  skill(neighborhood_mmd)) of the same candidate and seed (src/vec/evaluate.py `aggregate`, gate "floor_normalised"),
  so their weight is (0.25 / 3) x the entry's mean gate;
- the entry's mean skill over the seeds; the published floor's (copy_last on T1 and T2, wt_identity on T3; 0.5 by the
  skill map); the copy row's (real reference cells at the count rule) where the run has one;
- the points gained, 100 x weight x (entry skill - baseline skill), over the floor and over the copy, as seed-paired
  means of each metric's contribution to the board score;
- the points left to the split-half ceiling, 100 x weight x (1 - entry skill).
Per scored group it sums its metrics' points left and gains. The weakest weighted component is the group with the most
points left (docs/PLAN.md section 3, item 9), ties (at 0.01 points) going to the one with the least gained over the
floor; the weakest metric, the same rule over single metrics, is printed for reference only. Across boards it ranks
the groups in Overall units (a T2 board point is half an Overall point).
Contributions are computed per seed, so they sum to the board score. On T2, neighborhood_mmd's own weight is 0.25; its
effect through the gate is booked on the shape terms, whose points left are 100 x (0.25 / 3) x (1 - gate x skill)
(equal to the formula above whenever the gate is 1), so the points left of all metrics sum to 100 - board score.

Where the harness cannot score an entry, a stand-in is used and labelled:
- T1 s1 (M1): `m1_feng` of the labelled diagnostic experiments/t1_m1/results/diag_v1/ (not an adoption result; the
  harness refuses Feng's E9.5 cells on every T1 pseudo-task, and the readout that chose Feng ran on this transition);
- T2 heart s1 (H1): the harness's `b0`, which is H1 without the pre-registered size step (src/vec/harness.py
  `_b0_t2`: coordinates stay the reference's; no backtest can check a size);
- T3 (both entries): the closest rows of the pre-registered decomposition on the Gata4 replicates
  (experiments/t3_decomp/results/gata4/gata4_records.csv, the T3 script's own scorer, 10 seeds), next to the
  harness's published-baseline rows on the two counted cross-knockout folds.

Inputs (committed practice-task results only, read through `_src`, which refuses everything else): the task and
decision files under experiments/harness/results/ (d2_v1, emb_fgw_v1, emb_shape_v1, baselines_v1),
experiments/t1_m1/results/diag_v1/, experiments/t3_decomp/results/gata4/gata4_records.csv and results/registry.csv (a
cross-check of the means). It opens no target, no data file, nothing under boards/ (index.json above all), no
test-board file and nothing under docs/research/. Each harness task file embeds the spec of the validation or pseudo
board it was scored on; only that board's key is used. From every file it uses only skills, board scores, seeds, task
and candidate names, decision outcomes and provenance stamps, and it writes no cell count, anchor or raw metric value.
It runs no candidate on any backtest, so it is not a try in the sense of experiments/harness/TRIES.md, and nothing in it
chooses an entry: a change it points to still needs a ledger row before its run and the adoption rule after it.

Usage (stdlib only, from the repository root):

    python3 -I experiments/harness/headroom.py [--out-dir experiments/harness/results/headroom]

Written by the AI assistant on 2026-11-04 at the captain's request to internalise his guiding principle.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT_DEFAULT = "experiments/harness/results/headroom"

PRINCIPLE = ("Begin with the published copy-last-stage and unchanged-wild-type baselines across the boards you intend "
             "to enter. Improve the weakest weighted component before polishing what already works.")
QUESTION = "What meaningful change beyond the baseline does this predict, and what held-out evidence shows it?"

# task -> [(group id, group weight, [(metric, within-group weight), ...])]: copied from src/vec/boards.py GROUPS, ids
# included (proposal Table 1; docs/CHALLENGE.md "Metric panel"). The self-check below recomputes every stored board
# score. On T2, shape_growth is the gated group.
GROUPS = {
    "T1": [("de_recovery", 0.25, [("de_score", 1.0)]),
           ("direction", 0.25, [("de_direction", 1.0)]),
           ("distribution", 0.30, [("mmd_u", 1.0)]),
           ("covariance", 0.20, [("variogram", 1.0)])],
    "T2": [("expression_change", 0.25, [("de_score", 0.5), ("de_direction", 0.5)]),
           ("cell_state", 0.25, [("mmd_u", 0.6), ("variogram", 0.4)]),
           ("shape_growth", 0.25, [("d2_shape", 1 / 3), ("occupancy_dice", 1 / 3), ("scale_log_ratio", 1 / 3)]),
           ("local_organisation", 0.25, [("neighborhood_mmd", 1.0)])],
    "T3": [("response_identity", 0.30, [("de_score", 1.0)]),
           ("response_direction", 0.25, [("de_direction", 1.0)]),
           ("response_magnitude", 0.25, [("severity_slope", 1.0)]),
           ("mutant_state", 0.20, [("mmd_u", 0.6), ("variogram", 0.4)])],
}
WEIGHT = {t: {m: gw * w for _, gw, ms in g for m, w in ms} for t, g in GROUPS.items()}
METRICS = {t: [m for _, _, ms in g for m, _ in ms] for t, g in GROUPS.items()}
SHAPE = ("d2_shape", "occupancy_dice", "scale_log_ratio")
# Overall standing = T1 + mean(T2 heart, T2 embryo) + T3 (docs/PLAN.md section 1): a T2 board point is half a point.
OVERALL = {"T1:test": 1.0, "T2:heart:test_extrap": 0.5, "T2:embryo:test_interp": 0.5, "T3:bcatenin": 1.0}
TASK_OF = {"T1:test": "T1", "T2:heart:test_extrap": "T2", "T2:embryo:test_interp": "T2", "T3:bcatenin": "T3"}

READABLE = ("experiments/harness/results/", "experiments/t1_m1/results/diag_v1/", "experiments/t3_decomp/results/",
            "results/registry.csv")
REFUSED = ("boards/", "docs/research/", "data/", "runs/")
_READ: dict[str, str] = {}  # path -> sha256 of every input read


def _src(rel: str) -> Path:
    """The one door to the file system: committed practice-task results only."""
    if rel.startswith(REFUSED) or "index.json" in rel or not rel.startswith(READABLE):
        raise SystemExit(f"headroom: refusing {rel}: not a committed practice-task result")
    path = (REPO / rel).resolve()
    if REPO not in path.parents:
        raise SystemExit(f"headroom: refusing {rel}: outside the repository")
    _READ[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def _num(x) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return math.nan
    return v


def _skill(x) -> float:
    """A skill as the server counts it: NaN or missing is 0."""
    v = _num(x)
    return v if math.isfinite(v) else 0.0


def _mean(xs) -> float:
    xs = list(xs)
    return sum(xs) / len(xs)


# ------------------------------------------------------------------------------------------------------------ loaders
def load_harness(rel: str) -> dict:
    """A committed harness task file: per seed, every candidate's ranked skills and board score."""
    d = json.loads(_src(rel).read_text())
    rows = sorted(d["rows"], key=lambda r: r["seed"])
    return {"rel": rel, "kind": "harness", "run_id": d["run_id"], "key": d["task"]["key"],
            "test_board": d["task"]["test_board"], "scoring": d["board"]["key"],
            "seeds": [r["seed"] for r in rows], "cands": list(d["candidates"]),
            "skills": {c: [r["skills"][c] for r in rows] for c in d["candidates"]},
            "scores": {c: [r["scores"][c] for r in rows] for c in d["candidates"]},
            "commit": d["provenance"]["commit"][:7], "clean": bool(d["provenance"]["tree_clean"]),
            "protocol": f"{d['protocol']}, build {d['build']}, anchors {d['anchors']}"}


def load_diag(rel: str) -> dict:
    """The labelled M1 diagnostic (experiments/t1_m1/diagnose.py): the same per-seed layout as a harness task."""
    d = json.loads(_src(rel).read_text())
    rows = sorted(d["rows"], key=lambda r: r["seed"])
    cands = list(rows[0]["skills"])
    return {"rel": rel, "kind": "diagnostic", "run_id": "t1_m1_diag_v1", "key": d["task"], "test_board": "T1:test",
            "scoring": None, "seeds": [r["seed"] for r in rows], "cands": cands,
            "skills": {c: [r["skills"][c] for r in rows] for c in cands},
            "scores": {c: [r["scores"][c] for r in rows] for c in cands},
            "commit": d["provenance"]["commit"][:7], "clean": bool(d["provenance"]["tree_clean"]),
            "protocol": "bake, build main, anchors computed (diag_v1 README)", "label": d["label"]}


def load_t3(rel: str, variants: tuple[str, ...]) -> dict:
    """Per-seed skills of the named variants of the pre-registered T3 decomposition, one 'run' per replicate target.
    Only the skill columns and the board score are read."""
    per: dict = {}
    with open(_src(rel), newline="") as fh:
        for r in csv.DictReader(fh):
            if r["variant"] not in variants:
                continue
            sk = {m: _num(r["skill_" + m]) for m in METRICS["T3"]}
            per.setdefault(r["target"], {}).setdefault(r["variant"], {})[int(r["seed"])] = (sk, _num(r["t3"]))
    runs = {}
    for target, by_var in per.items():
        seeds = sorted(next(iter(by_var.values())))
        for v, rows in by_var.items():
            if sorted(rows) != seeds:
                raise SystemExit(f"headroom: {rel} {target} {v}: seeds differ")
        runs[target] = {"rel": rel, "kind": "t3_decomposition", "run_id": "t3_decomp gata4 (2026-10-28)",
                        "key": f"Gata4 {target} (E8.75, against WT E8.75)", "test_board": "T3:bcatenin",
                        "scoring": None, "seeds": seeds, "cands": list(by_var),
                        "skills": {v: [rows[s][0] for s in seeds] for v, rows in by_var.items()},
                        "scores": {v: [rows[s][1] for s in seeds] for v, rows in by_var.items()},
                        "commit": "see gata4_meta.json", "clean": True,
                        "protocol": "experiments/t3_decomp/t3_decompose.py's own scorer, 10 seeds"}
    return runs


def load_decision(rel: str, base: str, change: str) -> dict:
    d = json.loads(_src(rel).read_text())
    for x in d["decisions"]:
        if x["base"] == base and x["change"] == change:
            return {"rel": rel, "base": base, "change": change, "keep": x["keep"], "counted": list(x["counted"]),
                    "excluded": dict(x["excluded"]), "decided_on": x["decided_on"]}
    raise SystemExit(f"headroom: no decision {base} -> {change} in {rel}")


# ------------------------------------------------------------------------------------------------------------ scoring
def per_seed(task: str, skills: dict) -> tuple[dict, float, dict, dict]:
    """Skills (NaN as 0), the structure gate, each metric's contribution to the board score and its points left."""
    s = {m: _skill(skills.get(m)) for m in METRICS[task]}
    gate = min(1.0, 2.0 * s["neighborhood_mmd"]) if task == "T2" else 1.0
    contrib = {m: 100.0 * WEIGHT[task][m] * (gate if m in SHAPE else 1.0) * s[m] for m in METRICS[task]}
    left = {m: 100.0 * WEIGHT[task][m] - contrib[m] for m in METRICS[task]}
    return s, gate, contrib, left


def board_score(task: str, skills: dict) -> float:
    return sum(per_seed(task, skills)[2].values())


def headroom(task: str, entry: list, floor: list, copy: list | None) -> dict:
    """Seed-paired headroom of one entry against the floor (and the copy) on one practice task."""
    E = [per_seed(task, r) for r in entry]
    F = [per_seed(task, r) for r in floor]
    C = [per_seed(task, r) for r in copy] if copy is not None else None
    gates = [e[1] for e in E]
    rows = {}
    for m in METRICS[task]:
        rows[m] = {
            "weight": WEIGHT[task][m] * (_mean(gates) if m in SHAPE else 1.0),
            "entry": _mean(e[0][m] for e in E),
            "floor": _mean(f[0][m] for f in F),
            "copy": _mean(c[0][m] for c in C) if C else None,
            "gain_floor": _mean(e[2][m] - f[2][m] for e, f in zip(E, F, strict=True)),
            "gain_copy": _mean(e[2][m] - c[2][m] for e, c in zip(E, C, strict=True)) if C else None,
            "left": _mean(e[3][m] for e in E),
        }
    groups = group_totals(task, rows)
    return {"metrics": rows, "gate_mean": _mean(gates), "gate_min": min(gates),
            "score": _mean(sum(e[2].values()) for e in E), "floor_score": _mean(sum(f[2].values()) for f in F),
            "copy_score": _mean(sum(c[2].values()) for c in C) if C else None,
            "groups": groups, "weakest_group": weakest(groups), "weakest_metric": weakest(rows)}


def group_totals(task: str, rows: dict) -> dict:
    out = {}
    for name, _, ms in GROUPS[task]:
        out[name] = {"left": sum(rows[m]["left"] for m, _ in ms),
                     "gain_floor": sum(rows[m]["gain_floor"] for m, _ in ms)}
    return out


def _ranked(rows: dict) -> list[str]:
    """Most points left first; ties at 0.01 points go to the one with the least gained over the floor (metrics or
    groups)."""
    return sorted(rows, key=lambda k: (-round(rows[k]["left"], 2), round(rows[k]["gain_floor"], 2)))


def weakest(rows: dict) -> str:
    return _ranked(rows)[0]


def aggregate(task: str, blocks: list[dict]) -> dict:
    """Means over the counted practice tasks of one entry."""
    rows = {}
    for m in METRICS[task]:
        rows[m] = {k: (_mean(b["metrics"][m][k] for b in blocks)
                       if all(b["metrics"][m][k] is not None for b in blocks) else None)
                   for k in ("weight", "entry", "floor", "copy", "gain_floor", "gain_copy", "left")}
    groups = group_totals(task, rows)
    return {"metrics": rows, "score": _mean(b["score"] for b in blocks), "n_tasks": len(blocks),
            "floor_score": _mean(b["floor_score"] for b in blocks), "groups": groups,
            "weakest_group": weakest(groups), "weakest_metric": weakest(rows)}


# ------------------------------------------------------------------------------------------------------------ boards
H = "experiments/harness/results/"
T3_REC = "experiments/t3_decomp/results/gata4/gata4_records.csv"
T3_STAND_INS = ("_floor_wt_workcopy", "fac:none", "xrep:F+depth", "transfer:mab21l2_all",
                "blind:G2_bcat_central_d50")

BOARDS = [
    {"board": "T1:test", "task": "T1", "floor_model": "copy_last",
     "what": ("E12.5 scRNA. s1 = M1: real E10.5 cells at E10.5's own chamber composition, Feng's E9.5 -> E14.5 "
              "direction at 2/5 for its four groups, half the E9.5 -> E10.5 per-type change for the rest "
              "(METHOD_SUMMARY.md). s2 = M1 at cardiomyocytes 0.60 of heart cells (docs/S2_RULES.md)."),
     "entries": [
         {"slot": "s1", "name": "M1", "cand": "m1_feng", "source": "diag",
          "decision": (H + "d2_v1/decisions_T1_test.json", "floor", "b0"),
          "label": ("STAND-IN, a labelled diagnostic, not an adoption result: m1_feng of experiments/t1_m1/results/"
                    "diag_v1 (B0 plus Feng's direction at 1/5, the step's share of Feng's bracket on this "
                    "transition). The harness refuses Feng's E9.5 cells on every T1 pseudo-task, and the "
                    "pre-registered readout that chose Feng ran on this same transition, so it is in-sample for "
                    "that choice. M1 has B0's scope, so B0's counted set against the floor is used.")}],
     "context": [("b0", "harness", "B0: M1 without Feng's direction (the harness's own candidate, d2_v1)")],
     "tasks": [
         {"key": "T1/E9.5->E10.5:chamber", "harness": H + "d2_v1/tasks/T1_E9.5_E10.5_chamber.json",
          "diag": "experiments/t1_m1/results/diag_v1/T1_E9.5_E10.5_chamber.json"},
         {"key": "T1/E9.5->E10.5", "harness": H + "d2_v1/tasks/T1_E9.5_E10.5.json",
          "diag": "experiments/t1_m1/results/diag_v1/T1_E9.5_E10.5.json"}],
     "not_scored": ("s2 (M1 at 0.60 cardiomyocytes) is not a harness candidate; its evidence is the oracle-base "
                    "sensitivity grid in experiments/s2_hedge/results/t1_grid/ (not an adoption result)."),
     "closed": ("Cell numbers, sizes, scope and T1's composition are closed to tries after the 2026-11-04 exposure, "
                "except births by the legal route (docs/DATA_POLICY.md 2026-11-04, (b), (c)); the cardiomyocyte and "
                "blood shares stay as the s1 and s2 rules set them. State shifts (direction) stay open under the "
                "block, but no backtest can judge a change to Feng's direction (Feng's E9.5 lies inside every T1 "
                "window).")},
    {"board": "T2:heart:test_extrap", "task": "T2", "floor_model": "copy_last",
     "what": ("E12.5 MERFISH, extrapolation. s1 = H1: B0 on real E10.5 tissue with the pre-registered size "
              "(METHOD_SUMMARY.md). s2 = H2 as the half trim (docs/S2_RULES.md)."),
     "entries": [
         {"slot": "s1", "name": "H1", "cand": "b0", "source": "harness",
          "decision": (H + "d2_v1/decisions_T2_heart_test_extrap.json", "floor", "b0"),
          "label": ("STAND-IN: the harness's b0 is H1 without the pre-registered size step; its coordinates and size "
                    "stay the reference stage's (src/vec/harness.py _b0_t2), so its shape and scale skills equal a "
                    "copy's and say nothing about the size the test build applies, which no backtest can check.")}],
     "context": [],
     "tasks": [
         {"key": "T2h/<=E9.5->E10.5", "harness": H + "d2_v1/tasks/T2h_E9.5_E10.5.json"},
         {"key": "T2h/<=E8.75->E10.5", "harness": H + "d2_v1/tasks/T2h_E8.75_E10.5.json"}],
     "not_scored": ("s2 (H2 as the half trim) is not a harness candidate; its evidence is the oracle-core-size "
                    "sensitivity analysis in experiments/s2_hedge/results/heart_scope*/ (not an adoption result)."),
     "closed": ("Cell numbers, the size and the scope are closed to tries after the 2026-11-04 exposure "
                "(docs/DATA_POLICY.md 2026-11-04 and (b)); the composition trend was dropped by the captain on "
                "2026-11-01 (docs/PLAN.md 4.3, pipeline step 2); the depth line is closed (TRIES.md, H3).")},
    {"board": "T2:embryo:test_interp", "task": "T2", "floor_model": "copy_last",
     "what": ("E7.75 MERFISH, interpolation. s1 = fgw: displacement interpolation between E7.5 and E8.0 with a "
              "mixture of real cells of both, at the embryo size rule. s2 = fgw_shape0: fgw's mixture in E7.5's "
              "geometry (docs/S2_RULES.md)."),
     "entries": [
         {"slot": "s1", "name": "fgw", "cand": "fgw", "source": "harness",
          "decision": (H + "emb_fgw_v1/fgw_vs_b0/decisions_T2_embryo_test_interp.json", "b0", "fgw"),
          "label": "the entry itself, scored by the harness (emb_fgw_v1); kept over B0 by the adoption rule."},
         {"slot": "s2", "name": "fgw_shape0", "cand": "fgw_shape0", "source": "harness_s2",
          "decision": (H + "emb_shape_v1/decisions_T2_embryo_test_interp.json", "fgw", "fgw_shape0"),
          "label": ("the entry itself, scored by the harness (emb_shape_v1); chosen as s2 by the frozen rule's "
                    "better-of-two, not kept as an s1 (vetoed on E7.25).")}],
     "context": [],
     "tasks": [
         {"key": "T2e/E7.25+E8.0->E7.5", "harness": H + "emb_fgw_v1/tasks/T2e_E7.25_E8.0_E7.5.json",
          "harness_s2": H + "emb_shape_v1/tasks/T2e_E7.25_E8.0_E7.5.json"},
         {"key": "T2e/E6.75+E7.5->E7.25", "harness": H + "emb_fgw_v1/tasks/T2e_E6.75_E7.5_E7.25.json",
          "harness_s2": H + "emb_shape_v1/tasks/T2e_E6.75_E7.5_E7.25.json"}],
     "not_scored": None,
     "closed": ("Cell numbers, the size and the scope are closed to tries after the 2026-11-04 exposure "
                "(docs/DATA_POLICY.md 2026-11-04 and (b)); the coupling line is closed for the s1 unless the captain "
                "rules otherwise (TRIES.md, after E11), and so is the depth line (TRIES.md, H3).")},
]

T3 = {
    "board": "T3:bcatenin", "task": "T3", "floor_model": "wt_identity",
    "what": ("beta-catenin knockout at E8.75. s1 = the allowed branch's s1: the pooled Gata4 per-gene factor "
             "(exponent 1), the depth ratio and the central composition at the Gata4 magnitude (log2 1.38). "
             "s2 = the refused branch's s1: wild-type cells with the central composition at the Mab21l2 magnitude "
             "(log2 0.718), no factor (docs/S2_RULES.md; experiments/t3_decomp/results/gata4/gata4_plan.json)."),
    "decision": (H + "baselines_v1/decisions_T3_bcatenin.json", "floor", "copy"),
    "harness_tasks": [
        {"key": "T3/Gata4->Mab21l2", "harness": H + "baselines_v1/tasks/T3_Gata4_Mab21l2.json"},
        {"key": "T3/Mab21l2->Gata4rep2", "harness": H + "baselines_v1/tasks/T3_Mab21l2_Gata4rep2.json"}],
    "stand_ins": [
        {"slot": "s1", "name": "allowed s1", "cand": "xrep:F+depth",
         "label": ("STAND-IN A, optimistic: the Gata4 factor and depth ratio borrowed from the other replicate, scored "
                   "on each replicate. Same knockout, so a replicate fold, which never counts under the adoption "
                   "rule; it lacks s1's composition term (D8 measured that term at +0.22 and -0.43 points on top of "
                   "the borrowed factor; gata4_decisions.json). It shows the factor carries between two runs of one "
                   "knockout, which a processing effect shared by one run could explain (rules.yaml, the 2026-10-19 "
                   "note that adds rule 6), not that it carries to another knockout.")},
        {"slot": "s1", "name": "allowed s1", "cand": "transfer:mab21l2_all",
         "label": ("STAND-IN B, the only cross-knockout reading of a per-gene factor: the Mab21l2 factor on WT E8.75 "
                   "cells, scored on Gata4 (D11b, whose pre-registered bar is a board score of at least 60 on every "
                   "target; it scored 59.40 on rep1, so it failed). A different factor from s1's, but the situation "
                   "of the test: one knockout's factor applied to another knockout.")},
        {"slot": "s2", "name": "refused s1", "cand": "blind:G2_bcat_central_d50",
         "label": ("STAND-IN: wild-type cells with s2's five depleted states at depletion 0.5 (log2 1.0, against "
                   "s2's 0.608, log2 0.718), no factor, scored on Gata4: the beta-catenin composition rule read "
                   "on another knockout, written before any Gata4 file was opened (rules.yaml, A11).")}],
    "not_scored": ("Neither T3 entry is a harness candidate: every fold is one they were selected or fitted on "
                   "(METHOD_SUMMARY.md, Task 3), so they stand on the pre-registered rules "
                   "(experiments/t3_decomp/rules.yaml; docs/PLAN.md decision log 2026-11-01 (c) and (d))."),
    "closed": ("T3's composition, cell numbers and sizes are closed to tries after the 2026-11-04 exposure "
               "(docs/DATA_POLICY.md 2026-11-04 and (b)); the entries are fixed by the frozen pre-registered rules; a "
               "new T3 factor change faces the adoption rule on the two cross-knockout folds (docs/PLAN.md decision "
               "log 2026-11-01 (d)) and then the captain, since only his ruling departs from the frozen rules' "
               "outcome; experiments/t3_decomp/rules.yaml is not edited."),
}


# ------------------------------------------------------------------------------------------------------------ checks
def score_check(run: dict, task: str) -> float:
    """Largest gap between a stored board score and the one recomputed here from the stored skills."""
    gap = 0.0
    for c in run["cands"]:
        for sk, sc in zip(run["skills"][c], run["scores"][c], strict=True):
            gap = max(gap, abs(board_score(task, sk) - sc))
    return gap


def registry_means() -> dict:
    out = {}
    with open(_src("results/registry.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            out[(r["run_id"], r["pseudo_task"], r["candidate"])] = float(r["score_mean"])
    return out


def same_draws(a: dict, b: dict, cand: str) -> float:
    if a["seeds"] != b["seeds"]:
        raise SystemExit(f"headroom: {a['rel']} and {b['rel']} have different seeds")
    return max(abs(x - y) for x, y in zip(a["scores"][cand], b["scores"][cand], strict=True))


# ------------------------------------------------------------------------------------------------------------ text
def fmt(v, spec=".3f", signed=False) -> str:
    if v is None:
        return "-"
    return format(v, ("+" if signed else "") + spec)


def table_lines(task: str, hr: dict, copy_name: str = "copy") -> list[str]:
    out = [f"    {'metric':<17}{'weight':>8}{'entry':>8}{'floor':>8}{copy_name:>8}{'+floor':>9}{'+' + copy_name:>9}"
           f"{'left':>8}"]
    for m in METRICS[task]:
        r = hr["metrics"][m]
        out.append(f"    {m:<17}{r['weight']:>8.4f}{fmt(r['entry']):>8}{fmt(r['floor']):>8}{fmt(r['copy']):>8}"
                   f"{fmt(r['gain_floor'], '.2f', True):>9}{fmt(r['gain_copy'], '.2f', True):>9}{r['left']:>8.2f}")
    tot = {k: sum(hr["metrics"][m][k] for m in METRICS[task]) if all(hr["metrics"][m][k] is not None
           for m in METRICS[task]) else None for k in ("weight", "gain_floor", "gain_copy", "left")}
    out.append(f"    {'board (sum)':<17}{tot['weight']:>8.4f}{hr['score']:>8.2f}{'':>16}"
               f"{fmt(tot['gain_floor'], '.2f', True):>9}{fmt(tot['gain_copy'], '.2f', True):>9}{tot['left']:>8.2f}")
    g = hr["groups"]
    out += wrap("groups (left / +floor): " + "; ".join(f"{k} {v['left']:.2f} / {v['gain_floor']:+.2f}"
                                                       for k, v in g.items()), "    ")
    wg, wm = hr["weakest_group"], hr["weakest_metric"]
    out += wrap(f"weakest weighted component (group): {wg} ({g[wg]['left']:.2f} points left, "
                f"{g[wg]['gain_floor']:+.2f} gained over the floor); weakest metric: {wm} "
                f"({hr['metrics'][wm]['left']:.2f} left)", "    ")
    return out


def wrap(text: str, indent: str = "  ", width: int = 118) -> list[str]:
    words, lines, cur = text.split(), [], indent
    for w in words:
        if len(cur) + len(w) + 1 > width and cur.strip():
            lines.append(cur.rstrip())
            cur = indent + "  "
        cur += w + " "
    if cur.strip():
        lines.append(cur.rstrip())
    return lines


def agg_text(agg: dict, board: str) -> str:
    wg, wm = agg["weakest_group"], agg["weakest_metric"]
    g = agg["groups"][wg]
    return (f"weakest weighted component (group) {wg}, {g['left']:.2f} points left "
            f"({g['left'] * OVERALL[board]:.2f} in Overall), {g['gain_floor']:+.2f} gained over the floor; weakest "
            f"metric {wm} ({agg['metrics'][wm]['left']:.2f} left)")


def ranked_lines(task: str, agg: dict) -> list[str]:
    out = wrap("groups, points left, most first: " + ", ".join(
        f"{k} {agg['groups'][k]['left']:.2f}" for k in _ranked(agg["groups"])), "    ")
    out += wrap("metrics, points left, most first: " + ", ".join(
        f"{m} {agg['metrics'][m]['left']:.2f}" for m in _ranked(agg["metrics"])), "    ")
    return out


# ------------------------------------------------------------------------------------------------------------ main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default=OUT_DEFAULT, help="repository-relative folder for headroom.txt and .json")
    args = ap.parse_args(argv)

    reg = registry_means()
    checks, txt, result = [], [], {"principle": PRINCIPLE, "question": QUESTION, "boards": {}}
    worst_score_gap, worst_reg_gap = 0.0, 0.0

    def check_run(run: dict, task: str) -> None:
        nonlocal worst_score_gap, worst_reg_gap
        worst_score_gap = max(worst_score_gap, score_check(run, task))
        for c in run["cands"]:
            k = (run["run_id"], run["key"], c)
            if k in reg:
                worst_reg_gap = max(worst_reg_gap, abs(_mean(run["scores"][c]) - reg[k]))
            elif run["kind"] == "diagnostic" and c == "floor":
                pass  # diagnose.py registers its b0 and m1_feng rows only; the floor is 50 by construction
            elif run["kind"] != "t3_decomposition":
                checks.append(f"not in results/registry.csv: {run['run_id']} {run['key']} {c}")
        for c in ("floor", "_floor_wt_workcopy"):
            if c in run["cands"] and any(abs(s - 50.0) > 1e-9 for s in run["scores"][c]):
                raise SystemExit(f"headroom: {run['rel']}: the floor is not 50 on every seed")

    summary = []
    for cfg in BOARDS:
        task, board = cfg["task"], cfg["board"]
        txt += ["", "=" * 118, f"{board}  (weight in Overall: {OVERALL[board]:g}; floor = {cfg['floor_model']})"]
        txt += wrap(cfg["what"])
        bres = {"what": cfg["what"], "floor_model": cfg["floor_model"], "overall_factor": OVERALL[board],
                "entries": [], "closed_levers": cfg["closed"], "not_scored": cfg["not_scored"]}
        runs = {}
        for t in cfg["tasks"]:
            r = {"harness": load_harness(t["harness"])}
            if t["harness"] and r["harness"]["key"] != t["key"]:
                raise SystemExit(f"headroom: {t['harness']} holds {r['harness']['key']}, not {t['key']}")
            for extra in ("diag", "harness_s2"):
                if extra in t:
                    r[extra] = load_diag(t[extra]) if extra == "diag" else load_harness(t[extra])
                    if r[extra]["key"] != t["key"]:
                        raise SystemExit(f"headroom: {t[extra]} holds {r[extra]['key']}, not {t['key']}")
            for run in r.values():
                check_run(run, task)
            if "diag" in r:
                gap = same_draws(r["harness"], r["diag"], "b0")
                checks.append(f"{t['key']}: diag_v1's b0 equals d2_v1's b0 seed for seed (max gap {gap:.1e} points), "
                              "so its m1_feng pairs with the harness's floor and copy")
                if gap > 1e-6:
                    raise SystemExit("headroom: the diagnostic's draws are not the harness's")
            if "harness_s2" in r:
                gap = same_draws(r["harness"], r["harness_s2"], "fgw")
                checks.append(f"{t['key']}: emb_shape_v1's fgw equals emb_fgw_v1's fgw seed for seed "
                              f"(max gap {gap:.1e} points)")
                if gap > 1e-6:
                    raise SystemExit("headroom: the s2 run's draws are not the s1 run's")
            runs[t["key"]] = r

        for ent in cfg["entries"]:
            dec = load_decision(*ent["decision"])
            counted = [k for k in runs if k in dec["counted"]]
            eres = {"slot": ent["slot"], "entry": ent["name"], "candidate": ent["cand"], "label": ent["label"],
                    "counted_by": f"{dec['rel']} ({dec['base']} -> {dec['change']}: {dec['decided_on']})",
                    "excluded": dec["excluded"], "tasks": {}}
            txt += ["", f"  {ent['slot']} = {ent['name']} (candidate {ent['cand']})"]
            txt += wrap(ent["label"], "    ")
            txt += wrap(f"Counted practice tasks from {dec['rel']} ({dec['base']} -> {dec['change']}): "
                        f"{', '.join(dec['counted'])}; not counted: "
                        + ("; ".join(f"{k} ({v})" for k, v in dec["excluded"].items()) or "none") + ".", "    ")
            blocks = []
            for key, r in runs.items():
                run = r[ent["source"]]
                hrun = r["harness"]
                copy = hrun["skills"].get("copy")
                hr = headroom(task, run["skills"][ent["cand"]], hrun["skills"]["floor"], copy)
                status = "counted" if key in counted else ("reported, not counted: " + dec["excluded"].get(key, "?"))
                where = f"{run['rel']} (commit {run['commit']}, clean {run['clean']}, seeds {run['seeds'][0]}-" \
                        f"{run['seeds'][-1]}, {run['protocol']})"
                if run is not hrun:
                    where += f"; floor and copy from {hrun['rel']} (commit {hrun['commit']})"
                txt += ["", f"   task {key} [{status}]"]
                txt += wrap(f"source: {where}", "    ")
                txt.append(f"    board: entry {hr['score']:.2f}, floor {hr['floor_score']:.2f}, copy "
                           f"{fmt(hr['copy_score'], '.2f')}" + (f"; structure gate mean {hr['gate_mean']:.3f}, "
                                                                f"min {hr['gate_min']:.3f}" if task == "T2" else ""))
                for cand, src_name, what in cfg["context"]:
                    cr = r[src_name]
                    chr_ = headroom(task, cr["skills"][cand], hrun["skills"]["floor"], copy)
                    txt += wrap(f"context: {what}: board {chr_['score']:.2f}; skills " + ", ".join(
                        f"{m} {chr_['metrics'][m]['entry']:.3f}" for m in METRICS[task]), "    ")
                    eres.setdefault("context", {})[key] = {"candidate": cand, "what": what, "score": chr_["score"],
                                                           "skills": {m: chr_["metrics"][m]["entry"]
                                                                      for m in METRICS[task]}}
                txt += table_lines(task, hr)
                eres["tasks"][key] = dict(hr, status=status, source=where)
                if key in counted:
                    blocks.append(hr)
            if blocks:
                agg = aggregate(task, blocks)
                eres["aggregate"] = agg
                txt += [""] + wrap(f"{ent['slot']} over the {agg['n_tasks']} counted task(s): board "
                                   f"{agg['score']:.2f}; " + agg_text(agg, board), "   ")
                txt += ranked_lines(task, agg)
                summary.append((board, ent["slot"], ent["name"], ent["label"].startswith("STAND-IN"), agg))
            bres["entries"].append(eres)
        txt += [""] + wrap("Not in the table: " + (cfg["not_scored"] or "nothing; both entries are harness "
                                                    "candidates."), "  ")
        txt += wrap("Closed levers: " + cfg["closed"], "  ")
        result["boards"][board] = bres

    # T3: the published baselines on the counted folds, then the labelled stand-ins on the Gata4 replicates.
    task, board = "T3", T3["board"]
    txt += ["", "=" * 118, f"{board}  (weight in Overall: {OVERALL[board]:g}; floor = {T3['floor_model']})"]
    txt += wrap(T3["what"])
    dec = load_decision(*T3["decision"])
    bres = {"what": T3["what"], "floor_model": T3["floor_model"], "overall_factor": OVERALL[board],
            "baselines_on_counted_folds": {}, "stand_ins": [], "closed_levers": T3["closed"],
            "not_scored": T3["not_scored"]}
    txt += [""] + wrap(f"Baseline rows on the counted cross-knockout folds ({dec['rel']}: "
                       f"{', '.join(dec['counted'])}; not counted: "
                       + "; ".join(f"{k} ({v})" for k, v in dec["excluded"].items()) + "). The published baseline, "
                       "wt_identity, is the floor (50 by construction). Rows for our copy (real WT cells at the count "
                       "rule, a fresh draw, not a published baseline) and, as context, the organisers' published "
                       "shift_transfer; neither is an entry, and neither entry is a harness candidate.", "  ")
    for t in T3["harness_tasks"]:
        run = load_harness(t["harness"])
        if run["key"] != t["key"]:
            raise SystemExit(f"headroom: {t['harness']} holds {run['key']}")
        check_run(run, task)
        status = "counted" if t["key"] in dec["counted"] else "not counted"
        txt += ["", f"   task {t['key']} [{status}]"]
        txt += wrap(f"source: {run['rel']} (commit {run['commit']}, clean {run['clean']}, seeds {run['seeds'][0]}-"
                    f"{run['seeds'][-1]}, {run['protocol']})", "    ")
        fold = {}
        for cand in ("copy", "shift_transfer"):
            hr = headroom(task, run["skills"][cand], run["skills"]["floor"], None)
            fold[cand] = hr
            txt += wrap(f"{cand}: board {hr['score']:.2f} (floor 50.00); skills " + ", ".join(
                f"{m} {hr['metrics'][m]['entry']:.3f}" for m in METRICS[task]) + "; gained over the floor " +
                ", ".join(f"{m} {hr['metrics'][m]['gain_floor']:+.2f}" for m in METRICS[task]), "    ")
        fl = headroom(task, run["skills"]["floor"], run["skills"]["floor"], None)
        txt += wrap("floor (wt_identity), points left: " + ", ".join(
            f"{m} {fl['metrics'][m]['left']:.2f}" for m in METRICS[task]) + f" (weakest group: {fl['weakest_group']})",
            "    ")
        bres["baselines_on_counted_folds"][t["key"]] = {"status": status, "copy": fold["copy"],
                                                        "shift_transfer": fold["shift_transfer"], "floor": fl}

    t3runs = load_t3(T3_REC, T3_STAND_INS)
    for run in t3runs.values():
        gap = 0.0
        for c in run["cands"]:
            for sk, sc in zip(run["skills"][c], run["scores"][c], strict=True):
                gap = max(gap, abs(board_score(task, sk) - sc))
        worst_score_gap = max(worst_score_gap, gap)
        if any(abs(s - 50.0) > 1e-9 for s in run["scores"]["_floor_wt_workcopy"]):
            raise SystemExit("headroom: the T3 decomposition's floor is not 50 on every seed")
    txt += [""] + wrap(f"Stand-ins for the entries ({T3_REC}: the pre-registered decomposition on the two Gata4 "
                       "replicates, scored against WT E8.75 by the T3 script's own scorer, seeds 0-9; not in "
                       "results/registry.csv and not harness folds). The copy column is fac:none, the "
                       "decomposition's WT cells resampled.", "  ")
    for si in T3["stand_ins"]:
        txt += ["", f"  {si['slot']} = {si['name']}, stand-in {si['cand']}"]
        txt += wrap(si["label"], "    ")
        blocks, sres = [], {"slot": si["slot"], "entry": si["name"], "candidate": si["cand"], "label": si["label"],
                            "targets": {}}
        for target, run in sorted(t3runs.items()):
            hr = headroom(task, run["skills"][si["cand"]], run["skills"]["_floor_wt_workcopy"],
                          run["skills"]["fac:none"])
            txt += ["", f"   target {run['key']} [stand-in, not counted]"]
            txt.append(f"    board: stand-in {hr['score']:.2f}, floor {hr['floor_score']:.2f}, "
                       f"fac:none {fmt(hr['copy_score'], '.2f')}")
            txt += table_lines(task, hr, copy_name="wt_rs")
            sres["targets"][target] = hr
            blocks.append(hr)
        agg = aggregate(task, blocks)
        sres["aggregate_over_replicates"] = agg
        txt += [""] + wrap(f"over both replicates: board {agg['score']:.2f}; " + agg_text(agg, board), "   ")
        txt += ranked_lines(task, agg)
        summary.append((board, si["slot"], f"{si['name']} ~ {si['cand']}", True, agg))
        bres["stand_ins"].append(sres)
    txt += [""] + wrap("Not in the table: " + T3["not_scored"], "  ")
    txt += wrap("Closed levers: " + T3["closed"], "  ")
    result["boards"][board] = bres

    # Summary across boards, in board points and in Overall points (a T2 board point is half an Overall point).
    head = ["Baseline headroom per test board (experiments/harness/headroom.py; written by the AI assistant, "
            "2026-11-04)", "",
            *wrap(f"The captain's guiding principle (2026-11-04): \"{PRINCIPLE}\" Asked of every component: "
                  f"\"{QUESTION}\"", ""),
            "",
            *wrap("Weights from docs/CHALLENGE.md \"Metric panel\" (src/vec/boards.py GROUPS). Skills: 0.5 = the "
                  "published floor (copy_last on T1 and T2, wt_identity on T3), 1 = the split-half ceiling. +floor and "
                  "+copy: points gained over the floor and the copy row, 100 x weight x (entry skill - baseline "
                  "skill), seed-paired. left: points left to the ceiling, 100 x weight x (1 - entry skill). T2 shape "
                  "weights carry the entry's mean structure gate min(1, 2 x skill(neighborhood_mmd)); a shape term's "
                  "points left include any gate shortfall, so the points left sum to 100 - board score. Weakest "
                  "weighted component: the scored group (src/vec/boards.py GROUPS id) with the most points left, "
                  "ties at 0.01 points to the least gained over the floor (docs/PLAN.md section 3, item 9); weakest "
                  "metric: the same rule over single metrics, for reference only. Inputs are committed practice-task "
                  "results only; no target, no test-board file, no cell count.", ""),
            "", *wrap("Summary, means over the counted practice tasks (n; for T3 the two Gata4 replicates). * marks "
                      "a stand-in. score: board score; +floor: its gain over the floor's 50. weakest group: the "
                      "weakest weighted component, its points left on the board and in Overall units (x 0.5 on a T2 "
                      "board), and what it gained over the floor.", ""),
            f"  {'board':<23}{'slot':<5}{'entry':<42}{'n':>2}{'score':>7}{'+floor':>8}  {'weakest group':<19}"
            f"{'left':>6}{'Overall':>8}{'gained':>7}"]
    result["summary"] = []
    for board, slot, name, stand, agg in summary:
        wg, wm = agg["weakest_group"], agg["weakest_metric"]
        g = agg["groups"][wg]
        head.append(f"  {board:<23}{slot:<5}{(name + (' *' if stand else '')):<42}{agg['n_tasks']:>2}"
                    f"{agg['score']:>7.2f}{agg['score'] - agg['floor_score']:>+8.2f}  {wg:<19}{g['left']:>6.2f}"
                    f"{g['left'] * OVERALL[board]:>8.2f}{g['gain_floor']:>+7.2f}")
        result["summary"].append({"board": board, "slot": slot, "entry": name, "stand_in": stand,
                                  "n_tasks": agg["n_tasks"], "board_score": agg["score"],
                                  "gain_over_floor": agg["score"] - agg["floor_score"],
                                  "left_overall_board": (100.0 - agg["score"]) * OVERALL[board],
                                  "weakest_group": wg, "weakest_group_left": g["left"],
                                  "weakest_group_left_overall": g["left"] * OVERALL[board],
                                  "weakest_group_gain_floor": g["gain_floor"],
                                  "weakest_metric": wm, "weakest_metric_left": agg["metrics"][wm]["left"],
                                  "left_by_group": {k: v["left"] for k, v in agg["groups"].items()},
                                  "left_by_metric": {m: agg["metrics"][m]["left"] for m in agg["metrics"]}})
    head += ["", "Weakest metric (single metrics, for reference), where it is not in the weakest weighted component:"]
    diff = [(b, s, n, a) for b, s, n, _, a in summary
            if not any(m == a["weakest_metric"] for m, _ in next(ms for g, _, ms in GROUPS[TASK_OF[b]]
                                                                if g == a["weakest_group"]))]
    head += [f"  {b} {s} {n}: {a['weakest_metric']} {a['metrics'][a['weakest_metric']]['left']:.2f} left, in "
             f"another group than {a['weakest_group']} ({a['groups'][a['weakest_group']]['left']:.2f})"
             for b, s, n, a in diff] or ["  none"]
    s1 = [x for x in summary if x[1] == "s1"]
    head += ["", "Points left in Overall units per board (100 - board score, x 0.5 on a T2 board), s1 or its "
             "stand-ins:"]
    head += [f"  {(100.0 - agg['score']) * OVERALL[b]:6.2f}  {b} ({n})" for b, _, n, _, agg in s1]
    result["overall_left_by_board_s1"] = [{"left_overall": (100.0 - agg["score"]) * OVERALL[b], "board": b,
                                           "entry": n} for b, _, n, _, agg in s1]
    gpairs = sorted(((agg["groups"][k]["left"] * OVERALL[b], b, n, k) for b, _, n, _, agg in s1
                     for k in agg["groups"]), key=lambda x: (-x[0], x[1], x[2], x[3]))
    head += [""] + wrap("Points left per scored group in Overall units, s1 (or its stand-ins; T3 has two), largest "
                        "first (groups compared across boards, docs/PLAN.md section 3, item 9):", "")
    head += [f"  {v:6.2f}  {b} {k} ({n})" for v, b, n, k in gpairs]
    result["overall_ranking_s1_by_group"] = [{"left_overall": v, "board": b, "entry": n, "group": k}
                                             for v, b, n, k in gpairs]
    pairs = sorted(((agg["metrics"][m]["left"] * OVERALL[b], b, n, m) for b, _, n, _, agg in s1
                    for m in agg["metrics"]), reverse=True)
    head += ["", "Points left per metric in Overall units, s1 (or its stand-ins), largest first (top 12; for "
             "reference):"]
    head += [f"  {v:6.2f}  {b} {m} ({n})" for v, b, n, m in pairs[:12]]
    result["overall_ranking_s1_by_metric"] = [{"left_overall": v, "board": b, "entry": n, "metric": m}
                                              for v, b, n, m in pairs]

    checks.insert(0, "every stored board score (harness, diagnostic and T3 decomposition rows) recomputed from its "
                     f"skills with these weights and the gate: max gap {worst_score_gap:.1e} points")
    checks.insert(1, "every per-candidate mean used agrees with results/registry.csv score_mean: max gap "
                     f"{worst_reg_gap:.4f} points (the registry rounds to 0.001)")
    if worst_score_gap > 1e-6:
        raise SystemExit(f"headroom: recomputed board scores differ from the stored ones by {worst_score_gap}")
    caveats = [
        "Points left are measured to the split-half ceiling, which no prediction reaches, and the skill map is "
        "1/(1 + x): halving a metric's remaining distance takes its skill from 0.5 to 0.667, so equal points left on "
        "two metrics do not mean equal reach.",
        "Where an entry sits near the floor (the M1 and H1 stand-ins, T3's s2 stand-in), the ranking of points left "
        "mostly restates the weights: at the floor every metric leaves 50 x weight points, so the heaviest metric "
        "comes first (T1 mmd_u 0.30, T2 neighborhood_mmd 0.25, T3 de_score 0.30). On T2 the four groups weigh 0.25 "
        "each and neighborhood_mmd alone carries a whole group's weight, so the weakest metric is read for reference "
        "only and the weakest weighted component is the group (docs/PLAN.md section 3, item 9).",
        "Where a board's counting backtests name different weakest groups, those groups count as tied, and so do "
        "gaps across boards that lie within the spread of a board's stand-ins (experiments/harness/TRIES.md rule 6); "
        "the per-task lines above show each backtest's own weakest group.",
        "These are backtests on released transitions, not the test: T1 on E9.5 -> E10.5 (test E10.5 -> E12.5); the "
        "heart one and 1.75 days ahead (test two); the embryo at t = 1/3 and 2/3 (test 1/2). Each backtest's geometry, "
        "sampling and depth differ from the test's (TRIES.md diagnostics).",
        "T1: the stand-in draws E9.5's own chamber mix (cardiomyocytes about 66% of heart cells) against E10.5's "
        "chamber-cut truth, while the test entry starts from E10.5's own mix, so mmd_u and the variogram here carry "
        "a composition mismatch the test entry does not carry in that form. B0 drawn to the truth's own chamber mix "
        "(an oracle base) scores 69.92 against that truth, against 49.68 for B0 at E9.5's mix "
        "(experiments/s2_hedge/README.md, t1-grid at 90b25d3, 10 seeds; a sensitivity analysis, not an adoption "
        "result). At test the open question is how E12.5's states "
        "and shares differ from E10.5's, which is where births by the legal route come in.",
        "T1: the whole-E10.5 task is shown but not counted (a different tissue mix from the floor's); E8.5 -> E9.5 "
        "cannot run the entry.",
        "T2 heart: the stand-in's shape and scale skills are a copy's (the reference stage's coordinates); the test "
        "entry's size comes from the pre-registered procedure, which no backtest scores, and sizes are closed to "
        "tries anyway.",
        "T3: no row is an entry. Stand-in A is a replicate fold (one knockout, two runs), stand-in B a different "
        "factor; the composition magnitudes differ from the entries'. The harness rows on the counted folds are "
        "baseline rows only: wt_identity (the floor), our copy of WT cells and the organisers' shift_transfer.",
        "The weights and the gate are the emulator's (src/vec/evaluate.py, gate floor_normalised). docs/PLAN.md "
        "section 1 says all three board formulas, the T2 gate included, were re-derived from over 1,000 public "
        "leaderboard rows to within 1e-4; the organisers have not published the gate's code.",
        "This table chooses nothing: a change aimed at a weak component is registered in TRIES.md before its run "
        "and kept only under the adoption rule (docs/PLAN.md section 3, item 8) and TRIES.md rule 2; the data rules, "
        "the frozen s2 rules and the 2026-11-04 block stand.",
    ]
    tail = ["", "=" * 118, "Checks:"] + [ln for c in checks for ln in wrap("- " + c, "  ")]
    tail += ["", "Caveats:"] + [ln for c in caveats for ln in wrap("- " + c, "  ")]
    tail += ["", "Inputs (sha256):"] + [f"  {k}  {v}" for k, v in sorted(_READ.items())]
    result.update({"checks": checks, "caveats": caveats, "inputs_sha256": dict(sorted(_READ.items())),
                   "weights": WEIGHT, "definitions": {
                       "skill": "0.5 = the published floor, 1 = the split-half ceiling; NaN counts as 0",
                       "gain": "seed-paired mean of 100 x weight x (entry skill - baseline skill)",
                       "left": "100 x weight x (1 - entry skill); T2 shape terms: 100 x (0.25/3) x (1 - gate x skill)",
                       "gate": "T2: min(1, 2 x skill(neighborhood_mmd)) per candidate and seed",
                       "weakest_group": ("the weakest weighted component: the scored group (src/vec/boards.py GROUPS "
                                         "id) with the most points left; ties at 0.01 points to the least gained over "
                                         "the floor (docs/PLAN.md section 3, item 9)"),
                       "weakest_metric": "the same rule over single metrics, for reference only",
                       "left_overall": "points left x the board's weight in Overall (1 on T1 and T3, 0.5 on a T2 "
                                       "board)"}})

    text = "\n".join(head + txt + tail) + "\n"
    out = (REPO / args.out_dir).resolve()
    if REPO not in out.parents:
        raise SystemExit("headroom: the output folder must be inside the repository")
    out.mkdir(parents=True, exist_ok=True)
    (out / "headroom.txt").write_text(text)
    (out / "headroom.json").write_text(json.dumps(result, indent=1, default=float) + "\n")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
