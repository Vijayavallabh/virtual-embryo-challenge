#!/usr/bin/env python
"""Try T32 (experiments/harness/TRIES.md, registered before the run): the one T1 pair try T30 left unpriced. A labelled
diagnostic for a ruling; no upload changes, the builder still refuses the file this would name, and no cell's depth
changes in any arm.

Why. T30 priced three pairs on its grid (today's, the base hedge, the switched pair). Since row R6 the built pair is
the switched one: the first upload on library 3's cells at library 3's own mix, the second on the pooled base with
cardiomyocytes at 0.60 of heart cells. That pair matches the truth fully only where the truth is as deep as library 3
with its own mix, or shallower and cardiomyocyte-rich. Where it is as deep as library 3 AND cardiomyocyte-rich, neither
file matches on both counts. No arm of T30 puts the cardiomyocyte hedge on library 3's cells.

Geometry, scenarios, truths, floor, priors and seeds: T30's, imported unchanged (experiments/t1_scope/t30_pricing.py,
which imports T27's experiments/t1_scope/kind_gap.py).

Arms (n cells by the count rule; real cells; none from the truth half):
  lib3_base   T27's: library 3's cells at library 3's own chamber mix (rerun; must reproduce T30's stored scores)
  s1_cm60     T30's: the pooled base's draw rule at cardiomyocytes 0.60 of heart cells (rerun; must reproduce)
  lib3_cm60   NEW: library 3's cells at library 3's own chamber mix with cardiomyocytes moved to 0.60 of heart cells
              by vec.t1_scope.with_share, the rule the s2 build uses (blood's share kept)

Pairs, the better of two per seed: built = (lib3_base, s1_cm60), T30's "switched", the pair of row R6;
lib3_pair = (lib3_base, lib3_cm60).

Rule, fixed before the run: T30's, with the built pair in the place of "today's". Main prior census x flat; checks
latest x flat and census x like.
  - The lib3 pair is put to the captain only if its expected better-of-two score exceeds the built pair's by more than
    1.0 under the main prior and is not below the built pair's under both checks.
  - The built pair is confirmed if it beats the lib3 pair by more than 1.0 under the main prior and is not below it
    under both checks.
  - Otherwise: no recommendation; the captain is shown the table.
Reported, deciding nothing: all nine prior products; the two pairs at library 3's depth alone (the depth-1.0 row, the
numbers an organiser's answer "prepared with library _3" or "the same protocol" would be read on); under each depth
prior with the flat cell-mix prior, the weight of scenarios whose mean better-of-two is 95 or more; how often the new
arm draws a cell more than once (its cardiomyocyte quota exceeds half of library 3's cardiomyocytes, the pool here).
Void if a scenario, seed or arm is missing, results come from different commits or a dirty tree, or a rerun arm
differs from T30's stored per-seed score by more than 0.02.

  run   --depth F --cm own|0.45|0.60 --out DIR   one scenario (start the twelve together: one commit)
  read  --runs DIR --t30 DIR --out DIR           apply the rule; write pricing.json and pricing.txt
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kind_gap as kg  # noqa: E402  (T27's geometry and draws, unchanged)
import t30_pricing as t30  # noqa: E402  (T30's scenarios, priors and rule constants, unchanged)

from vec import harness, t1_scope  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import RANKED_METRICS  # noqa: E402
from vec.core import quotas  # noqa: E402
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample  # noqa: E402

ARMS = ("lib3_base", "s1_cm60", "lib3_cm60")
RERUN = ("lib3_base", "s1_cm60")  # must reproduce T30's stored per-seed scores
NEW = "lib3_cm60"
LIB3_CM60_STREAM = 11  # the new arm's draw stream, apart from kind_gap's arm indices and T30's 10
PAIRS = {"built": ("lib3_base", "s1_cm60"), "lib3_pair": ("lib3_base", "lib3_cm60")}
REPRO_TOL = 0.02  # points per seed, the tolerance of the ledger's reproduction gates (T26)
HIGH = 95.0  # R6's reason (4): how often a pair reaches 95 or more


def prepare(data_dir: Path, depth: float, cm: str) -> dict:
    d = t30.prepare(data_dir, depth, cm)
    d["lib3_cm60_mix"] = t1_scope.with_share(d["lib3_mix"], "cm", t30.CM_HEDGE)
    return d


def arm_rows(arm: str, n: int, seed: int) -> np.ndarray:
    if arm != NEW:
        return t30.arm_rows(arm, n, seed)  # T30's and T27's draws, stream for stream
    d = kg._D
    rng = np.random.default_rng([seed, LIB3_CM60_STREAM])
    q = quotas(d["lib3_cm60_mix"], n)
    return np.sort(np.concatenate([kg._draw(d["pools"]["b"][c], int(q.get(c, 0)), rng) for c in d["chamber"]]))


def score_seed(seed: int) -> dict:
    d = kg._D
    t0 = time.time()
    X, board = d["X"], d["board"]
    target, ref = Stage(X=d["x_truth"]), Stage(X=X[d["ref_rows"]])
    f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, seed), target, ref, seed=seed, anchors="computed",
                    build=harness.BUILD, protocol=harness.PROTOCOL)
    if abs(f.score - 50.0) > 1e-9 or set(RANKED_METRICS["T1"]) - set(f.anchors):
        raise ValueError(f"seed {seed}: the floor scores {f.score}")
    n = harness.count_for(board, len(d["ref_rows"]))
    row = {"seed": seed, "scores": {}, "skills": {}, "median_library": {}, "repeat_share": {}}
    for arm in ARMS:
        rows = arm_rows(arm, n, harness.BUILD_SEED_OFFSET + seed)
        r = score_board(board, Stage(X=X[rows]), target, ref, seed=seed, anchors=f.anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        row["scores"][arm], row["skills"][arm] = r.score, r.skills
        row["median_library"][arm] = float(np.median(d["tot"][rows]))
        row["repeat_share"][arm] = float(1.0 - len(np.unique(rows)) / len(rows))
    row["seconds"] = round(time.time() - t0, 1)
    print(f"seed {seed}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["scores"].items()), flush=True)
    return row


def cmd_run(a) -> None:
    if a.depth not in t30.DEPTHS or a.cm not in t30.CMS:
        sys.exit(f"--depth must be one of {t30.DEPTHS} and --cm one of {t30.CMS}")
    commit = harness._git("rev-parse", "HEAD")
    clean = harness._git("status", "--porcelain", "--untracked-files=no") == ""
    kg._D.update(prepare(Path(a.data), a.depth, a.cm))
    d = kg._D
    if a.jobs > 1:
        with mp.get_context("fork").Pool(a.jobs) as pool:
            rows = pool.map(score_seed, list(SEEDS))
    else:
        rows = [score_seed(s) for s in SEEDS]
    rows = sorted(rows, key=lambda r: r["seed"])
    res = {"try": "T32", "commit": commit, "tree_clean": clean, "depth": a.depth, "cm": a.cm,
           "truth_cells": d["truth_cells"], "truth_cm_of_heart": d["truth_cm_of_heart"],
           "seeds": list(SEEDS), "arms": list(ARMS), "protocol": harness.PROTOCOL, "build": harness.BUILD,
           "frac": SUBSAMPLE_FRAC, "library3_chamber_mix": d["lib3_mix"], "cm60_mix": d["cm60_mix"],
           "lib3_cm60_mix": d["lib3_cm60_mix"],
           "library3_pool_cells": {c: int(len(d["pools"]["b"][c])) for c in d["chamber"]}, "rows": rows}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"t32_d{a.depth}_cm{a.cm}.json"
    path.write_text(json.dumps(res, indent=1))
    print(f"wrote {path}", flush=True)


def reproduction(runs: dict, stored: dict) -> dict:
    """The rerun arms against T30's stored per-seed scores: the largest absolute difference, and whether it is within
    REPRO_TOL in every scenario and seed."""
    worst, missing = 0.0, []
    for key, r in runs.items():
        s = stored.get(key)
        if s is None:
            missing.append(f"d{key[0]} cm{key[1]}")
            continue
        old = {x["seed"]: x["scores"] for x in s["rows"]}
        for x in r["rows"]:
            for arm in RERUN:
                worst = max(worst, abs(x["scores"][arm] - old[x["seed"]][arm]))
    return {"largest_difference": worst, "missing_in_t30": missing, "ok": bool(not missing and worst <= REPRO_TOL)}


def reading(runs: dict, stored: dict) -> dict:
    why = []
    for d in t30.DEPTHS:
        for c in t30.CMS:
            r = runs.get((d, c))
            if r is None:
                why.append(f"missing scenario depth {d}, cm {c}")
            elif sorted(x["seed"] for x in r["rows"]) != list(SEEDS) or any(
                    a not in x["scores"] for x in r["rows"] for a in ARMS):
                why.append(f"depth {d}, cm {c}: a seed or an arm is missing")
    if len({r["commit"] for r in runs.values()}) > 1:
        why.append("the results come from different commits")
    if not all(r["tree_clean"] for r in runs.values()):
        why.append("a result comes from a dirty tree")
    repro = reproduction(runs, stored) if not why else None
    if repro is not None and not repro["ok"]:
        why.append(f"a rerun arm does not reproduce T30 (largest difference {repro['largest_difference']:.4f}; "
                   f"missing {repro['missing_in_t30']})")
    if why:
        return {"void": True, "why": why, "reproduction": repro, "recommendation": "none: the reading is void"}
    arm, pair, ahead, repeats = {a: {} for a in ARMS}, {p: {} for p in PAIRS}, {}, {}
    for key, r in runs.items():
        rows = sorted(r["rows"], key=lambda x: x["seed"])
        s = {a: np.array([x["scores"][a] for x in rows], float) for a in ARMS}
        for a in ARMS:
            arm[a][key] = float(s[a].mean())
        for p, (x, y) in PAIRS.items():
            pair[p][key] = float(np.maximum(s[x], s[y]).mean())
        diff = np.maximum(s["lib3_base"], s[NEW]) - np.maximum(s["lib3_base"], s["s1_cm60"])
        ahead[key] = {"mean": float(diff.mean()), "ahead": int((diff > 0).sum()), "behind": int((diff < 0).sum())}
        repeats[key] = float(np.mean([x["repeat_share"][NEW] for x in rows]))
    priors = [(dp, cp) for dp in t30.DEPTH_PRIORS for cp in t30.CM_PRIORS]
    exp_pair = {p: {f"{dp} x {cp}": t30._expect(pair[p], dp, cp) for dp, cp in priors} for p in PAIRS}
    exp_arm = {a: {f"{dp} x {cp}": t30._expect(arm[a], dp, cp) for dp, cp in priors} for a in ARMS}
    main = " x ".join(t30.MAIN)
    checks = [" x ".join(c) for c in t30.CHECKS]
    lead = {k: exp_pair["lib3_pair"][k] - exp_pair["built"][k] for k in exp_pair["built"]}
    qualifies = bool(lead[main] > t30.BAR and all(lead[c] >= 0 for c in checks))
    confirmed = bool(-lead[main] > t30.BAR and all(lead[c] <= 0 for c in checks))
    if qualifies:
        rec = (f"put to the captain: the pair lib3_pair {PAIRS['lib3_pair']} in place of the built pair "
               f"{PAIRS['built']}")
    elif confirmed:
        rec = "the built pair is confirmed"
    else:
        rec = "no recommendation: the pairs are within the priors' reach; the captain is shown the table"
    high = {p: {dp: float(sum(t30.DEPTH_PRIORS[dp][i] * t30.CM_PRIORS["flat"][j] * (pair[p][(d, c)] >= HIGH)
                              for i, d in enumerate(t30.DEPTHS) for j, c in enumerate(t30.CMS)))
                for dp in t30.DEPTH_PRIORS} for p in PAIRS}
    at_lib3_depth = {c: {p: pair[p][(1.0, c)] for p in PAIRS} for c in t30.CMS}
    name = lambda k: f"d{k[0]} cm{k[1]}"  # noqa: E731
    return {"void": False, "why": [], "reproduction": repro,
            "arm_by_scenario": {a: {name(k): v for k, v in t.items()} for a, t in arm.items()},
            "pair_by_scenario": {p: {name(k): v for k, v in t.items()} for p, t in pair.items()},
            "lib3_pair_minus_built_by_scenario": {name(k): v for k, v in ahead.items()},
            "expected_pair": exp_pair, "expected_arm": exp_arm, "lib3_pair_minus_built": lead,
            "main_prior": main, "check_priors": checks, "qualifies": qualifies, "built_confirmed": confirmed,
            "recommendation": rec, "weight_at_or_above_95_flat_cm": high, "at_library3_depth": at_lib3_depth,
            "new_arm_repeat_share": {name(k): v for k, v in repeats.items()},
            "commit": next(iter(runs.values()))["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Try T32 (commit {res['commit']}): {res['recommendation']}",
             f"rerun arms against T30's stored scores: largest difference "
             f"{res['reproduction']['largest_difference']:.4f} points (tolerance {REPRO_TOL})",
             f"expected better-of-two score by prior (main: {res['main_prior']}; checks: "
             f"{', '.join(res['check_priors'])}):",
             "  prior (depth x cm)            built    lib3_pair   lib3_pair - built   | single arms: "
             + "  ".join(f"{a:>9}" for a in ARMS)]
    for k in res["expected_pair"]["built"]:
        built, lib3 = res["expected_pair"]["built"][k], res["expected_pair"]["lib3_pair"][k]
        lines.append(f"  {k:<22}  {built:11.2f}  {lib3:11.2f}"
                     f"  {res['lib3_pair_minus_built'][k]:+18.2f}   |              "
                     + "  ".join(f"{res['expected_arm'][a][k]:9.2f}" for a in ARMS))
    lines.append("mean score by scenario (truth depth relative to library 3, truth cardiomyocyte share); the pair "
                 "difference with the seeds ahead and behind; the new arm's share of repeated draws:")
    lines.append("  scenario        " + "  ".join(f"{a:>9}" for a in ARMS)
                 + "  |       built    lib3_pair   difference (ahead/behind)   repeats")
    for k in res["arm_by_scenario"][ARMS[0]]:
        dd = res["lib3_pair_minus_built_by_scenario"][k]
        lines.append(f"  {k:<14}  " + "  ".join(f"{res['arm_by_scenario'][a][k]:9.2f}" for a in ARMS) + "  | "
                     + f"{res['pair_by_scenario']['built'][k]:11.2f}  {res['pair_by_scenario']['lib3_pair'][k]:11.2f}"
                     f"   {dd['mean']:+8.2f} ({dd['ahead']}/{dd['behind']})            "
                     f"{res['new_arm_repeat_share'][k]:.2f}")
    lines.append("reported, deciding nothing:")
    lines.append("  at library 3's depth alone (depth 1.0), built / lib3_pair: "
                 + "; ".join(f"cm {c}: {v['built']:.2f} / {v['lib3_pair']:.2f}"
                             for c, v in res["at_library3_depth"].items()))
    lines.append("  weight of scenarios with a mean better-of-two of 95 or more (flat cell-mix prior), built / "
                 "lib3_pair: " + "; ".join(f"{dp}: {res['weight_at_or_above_95_flat_cm']['built'][dp]:.2f} / "
                                           f"{res['weight_at_or_above_95_flat_cm']['lib3_pair'][dp]:.2f}"
                                           for dp in t30.DEPTH_PRIORS))
    return "\n".join(lines) + "\n"


def _load(folder: Path, stem: str) -> dict:
    runs = {}
    for p in sorted(folder.glob(f"{stem}_d*_cm*.json")):
        r = json.loads(p.read_text())
        runs[(float(r["depth"]), str(r["cm"]))] = r
    return runs


def cmd_read(a) -> None:
    res = reading(_load(Path(a.runs), "t32"), _load(Path(a.t30), "t30"))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "pricing.json").write_text(json.dumps(res, indent=1))
    (out / "pricing.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--depth", type=float, required=True)
    r.add_argument("--cm", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--jobs", type=int, default=3)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--t30", default=str(Path(__file__).resolve().parent / "results" / "t30_pricing" / "scenarios"))
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
