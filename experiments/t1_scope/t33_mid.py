#!/usr/bin/env python
"""Try T33 (experiments/harness/TRIES.md, registered before the run): a base half-way between the pooled base and
library 3, priced on try T30's grid. A labelled diagnostic for a ruling; no upload changes, no builder option exists
for these files, and no cell's depth changes in any arm.

Why. Tries T30 and T32 priced second uploads on the pooled base (65% library 3's cells) and on library 3 alone. Which
of the two leads turns on the truth's depth, which no released file shows: T32's lib3 pair gains 4.5 to 9.2 where the
truth is as deep as library 3 or deeper and cardiomyocyte-rich, and loses 0.7 to 5.3 where it is shallower. The scored
distribution metrics are quadratic in a mixture's weights (the kernel distance and the variogram difference of a blend
of two populations against one of them fall with the square of the other's share), so a base half-way between the two
should give up about a quarter to a half of each gap instead of all of one. That is a hedge along the depth axis
inside one file, which two uploads spent on the cell mix cannot supply.

Geometry, scenarios, truths, floor, priors, seeds and the rule's constants: T30's, imported unchanged.

The mid base. In each chamber category c, with s_c library 3's share of the released file's cells of that category
(kind_gap's `lib3_share`), the cells are drawn at library 3's share (s_c + 1) / 2, the midpoint between the pooled
base's s_c and library 3's 1. A builder would get this by keeping every library-3 cell and a random fraction
s_c / (1 + s_c) of the other libraries' cells; the base's own chamber mix is then the pooled mix re-weighted by
2 s_c / (1 + s_c), which is the mix used here.

Arms (n cells by the count rule; real cells; none from the truth half):
  lib3_base   T27's (rerun; must reproduce T30's stored scores; the first upload of row R6)
  mid_base    NEW: the mid base at its own chamber mix
  mid_cm60    NEW: the mid base with cardiomyocytes moved to 0.60 of heart cells by vec.t1_scope.with_share
  s1_cm60, lib3_cm60   not rerun: their per-seed scores are read from T32's stored results (same machine, same seeds)

Pairs, the better of two per seed:
  built       (lib3_base, s1_cm60)     the pair of row R6
  lib3_pair   (lib3_base, lib3_cm60)   T32's
  mid_s2      (lib3_base, mid_cm60)    the first upload as built, the second on the mid base
  mid_both    (mid_base, mid_cm60)     both on the mid base

Rule, fixed before the run: T30's, with the built pair in the place of "today's". Main prior census x flat; checks
latest x flat and census x like.
  - A pair other than the built pair is put to the captain only if its expected better-of-two score exceeds the built
    pair's by more than 1.0 under the main prior and is not below the built pair's under both checks. If two qualify,
    the higher under the main prior; mid_both goes forward only if it also exceeds mid_s2 by more than 1.0 there, since
    it would change the first upload too (row R6).
  - The built pair is confirmed if it beats every other pair by more than 1.0 under the main prior and is not below
    any under both checks.
  - Otherwise: no recommendation; the captain is shown the table.
Reported, deciding nothing: all nine prior products; each pair's smallest lead over the nine; the largest loss and
gain by scenario; the weight of scenarios with a mean better-of-two of 95 or more; the new arms' shares of repeated
draws and median libraries (ten-seed means).
Void if a scenario, seed or arm is missing, results come from different commits or a dirty tree, or the rerun arm
differs from T30's stored per-seed score by more than 0.02.

  run   --depth F --cm own|0.45|0.60 --out DIR   one scenario (start the twelve together: one commit)
  read  --runs DIR --out DIR                     apply the rule; write pricing.json and pricing.txt
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

ARMS = ("lib3_base", "mid_base", "mid_cm60")  # the arms this try runs
RERUN = "lib3_base"
STORED = ("s1_cm60", "lib3_cm60")  # read from T32's stored per-seed scores
STREAMS = {"mid_base": 12, "mid_cm60": 13}  # apart from kind_gap's arm indices, T30's 10 and T32's 11
PAIRS = {"built": ("lib3_base", "s1_cm60"), "lib3_pair": ("lib3_base", "lib3_cm60"),
         "mid_s2": ("lib3_base", "mid_cm60"), "mid_both": ("mid_base", "mid_cm60")}
OTHERS = ("lib3_pair", "mid_s2", "mid_both")
REPRO_TOL = 0.02
HIGH = 95.0
HERE = Path(__file__).resolve().parent
T30_SCENARIOS = HERE / "results" / "t30_pricing" / "scenarios"
T32_SCENARIOS = HERE / "results" / "t32_pair" / "scenarios"


def mid_lambda(d: dict) -> dict[str, float]:
    """Library 3's share of each chamber category's draws on the mid base: the midpoint of its share in the file
    and 1."""
    return {c: (d["lib3_share"][c] + 1.0) / 2.0 for c in d["chamber"]}


def mid_mix(d: dict) -> dict[str, float]:
    """The mid base's own chamber mix: the pooled mix re-weighted by 2 s / (1 + s) per category (every library-3 cell
    kept and s / (1 + s) of the others) and renormalised."""
    w = {c: d["pooled_mix"][c] * 2.0 * d["lib3_share"][c] / (1.0 + d["lib3_share"][c]) for c in d["chamber"]}
    tot = sum(w.values())
    return {c: v / tot for c, v in w.items()}


def prepare(data_dir: Path, depth: float, cm: str) -> dict:
    d = t30.prepare(data_dir, depth, cm)
    d["mid_lambda"] = mid_lambda(d)
    d["mid_mix"] = mid_mix(d)
    d["mid_cm60_mix"] = t1_scope.with_share(d["mid_mix"], "cm", t30.CM_HEDGE)
    return d


def arm_rows(arm: str, n: int, seed: int) -> np.ndarray:
    if arm not in STREAMS:
        return t30.arm_rows(arm, n, seed)  # T30's and T27's draws, stream for stream
    d = kg._D
    rng = np.random.default_rng([seed, STREAMS[arm]])
    q = quotas(d["mid_mix"] if arm == "mid_base" else d["mid_cm60_mix"], n)
    out = []
    for c in d["chamber"]:
        k = int(q.get(c, 0))
        k3 = int(round(k * d["mid_lambda"][c]))
        out += [kg._draw(d["pools"]["b"][c], k3, rng), kg._draw(d["pools"]["o"][c], k - k3, rng)]
    return np.sort(np.concatenate(out))


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
    row = {"seed": seed, "scores": {}, "skills": {}, "median_library": {}, "repeat_share": {}, "lib3_draw_share": {}}
    lib3 = d["lib"] == kg.DISSECTED
    for arm in ARMS:
        rows = arm_rows(arm, n, harness.BUILD_SEED_OFFSET + seed)
        r = score_board(board, Stage(X=X[rows]), target, ref, seed=seed, anchors=f.anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        row["scores"][arm], row["skills"][arm] = r.score, r.skills
        row["median_library"][arm] = float(np.median(d["tot"][rows]))
        row["repeat_share"][arm] = float(1.0 - len(np.unique(rows)) / len(rows))
        row["lib3_draw_share"][arm] = float(lib3[rows].mean())
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
    res = {"try": "T33", "commit": commit, "tree_clean": clean, "depth": a.depth, "cm": a.cm,
           "truth_cells": d["truth_cells"], "truth_cm_of_heart": d["truth_cm_of_heart"],
           "seeds": list(SEEDS), "arms": list(ARMS), "protocol": harness.PROTOCOL, "build": harness.BUILD,
           "frac": SUBSAMPLE_FRAC, "library3_share_of_category": d["lib3_share"], "mid_lambda": d["mid_lambda"],
           "pooled_chamber_mix": d["pooled_mix"], "mid_mix": d["mid_mix"], "mid_cm60_mix": d["mid_cm60_mix"],
           "pool_cells": {k: {c: int(len(d["pools"][k][c])) for c in d["chamber"]} for k in ("b", "o")},
           "rows": rows}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"t33_d{a.depth}_cm{a.cm}.json"
    path.write_text(json.dumps(res, indent=1))
    print(f"wrote {path}", flush=True)


def _load(folder: Path, stem: str) -> dict:
    runs = {}
    for p in sorted(Path(folder).glob(f"{stem}_d*_cm*.json")):
        r = json.loads(p.read_text())
        runs[(float(r["depth"]), str(r["cm"]))] = r
    return runs


def reading(runs: dict, t30_runs: dict, t32_runs: dict) -> dict:
    why = []
    for d in t30.DEPTHS:
        for c in t30.CMS:
            r = runs.get((d, c))
            if r is None:
                why.append(f"missing scenario depth {d}, cm {c}")
            elif sorted(x["seed"] for x in r["rows"]) != list(SEEDS) or any(
                    a not in x["scores"] for x in r["rows"] for a in ARMS):
                why.append(f"depth {d}, cm {c}: a seed or an arm is missing")
            if (d, c) not in t30_runs or (d, c) not in t32_runs:
                why.append(f"depth {d}, cm {c}: T30's or T32's stored scenario is missing")
    if len({r["commit"] for r in runs.values()}) > 1:
        why.append("the results come from different commits")
    if not all(r["tree_clean"] for r in runs.values()):
        why.append("a result comes from a dirty tree")
    worst = 0.0
    if not why:
        for key, r in runs.items():
            old = {x["seed"]: x["scores"][RERUN] for x in t30_runs[key]["rows"]}
            worst = max(worst, max(abs(x["scores"][RERUN] - old[x["seed"]]) for x in r["rows"]))
        if worst > REPRO_TOL:
            why.append(f"the rerun arm does not reproduce T30 (largest difference {worst:.4f})")
    if why:
        return {"void": True, "why": why, "recommendation": "none: the reading is void"}
    names = (*ARMS, *STORED)
    arm, pair, extra = {a: {} for a in names}, {p: {} for p in PAIRS}, {}
    for key, r in runs.items():
        rows = sorted(r["rows"], key=lambda x: x["seed"])
        s = {a: np.array([x["scores"][a] for x in rows], float) for a in ARMS}
        stored = {x["seed"]: x["scores"] for x in t32_runs[key]["rows"]}
        for a in STORED:
            s[a] = np.array([stored[x["seed"]][a] for x in rows], float)
        for a in names:
            arm[a][key] = float(s[a].mean())
        for p, (x, y) in PAIRS.items():
            pair[p][key] = float(np.maximum(s[x], s[y]).mean())
        built = np.maximum(s["lib3_base"], s["s1_cm60"])
        extra[key] = {p: {"mean": float((np.maximum(s[x], s[y]) - built).mean()),
                          "ahead": int((np.maximum(s[x], s[y]) > built).sum()),
                          "behind": int((np.maximum(s[x], s[y]) < built).sum())}
                      for p, (x, y) in PAIRS.items() if p != "built"}
        extra[key]["new_arms"] = {a: {k: float(np.mean([x[k][a] for x in rows]))
                                      for k in ("repeat_share", "median_library", "lib3_draw_share")}
                                  for a in ("mid_base", "mid_cm60")}
    priors = [(dp, cp) for dp in t30.DEPTH_PRIORS for cp in t30.CM_PRIORS]
    exp_pair = {p: {f"{dp} x {cp}": t30._expect(pair[p], dp, cp) for dp, cp in priors} for p in PAIRS}
    exp_arm = {a: {f"{dp} x {cp}": t30._expect(arm[a], dp, cp) for dp, cp in priors} for a in names}
    main = " x ".join(t30.MAIN)
    checks = [" x ".join(c) for c in t30.CHECKS]
    lead = {p: {k: exp_pair[p][k] - exp_pair["built"][k] for k in exp_pair["built"]} for p in OTHERS}
    qualifies = {p: bool(lead[p][main] > t30.BAR and all(lead[p][c] >= 0 for c in checks)) for p in OTHERS}
    if qualifies["mid_both"] and not exp_pair["mid_both"][main] - exp_pair["mid_s2"][main] > t30.BAR:
        qualifies["mid_both"] = False  # it would change the first upload too: it must also beat mid_s2 by the bar
    confirmed = bool(all(-lead[p][main] > t30.BAR and all(lead[p][c] <= 0 for c in checks) for p in OTHERS))
    winners = [p for p in OTHERS if qualifies[p]]
    if winners:
        best = max(winners, key=lambda p: exp_pair[p][main])
        rec = f"put to the captain: the pair {best} {PAIRS[best]} in place of the built pair {PAIRS['built']}"
    elif confirmed:
        rec = "the built pair is confirmed"
    else:
        rec = "no recommendation: the pairs are within the priors' reach; the captain is shown the table"
    high = {p: {dp: float(sum(t30.DEPTH_PRIORS[dp][i] * t30.CM_PRIORS["flat"][j] * (pair[p][(d, c)] >= HIGH)
                              for i, d in enumerate(t30.DEPTHS) for j, c in enumerate(t30.CMS)))
                for dp in t30.DEPTH_PRIORS} for p in PAIRS}
    name = lambda k: f"d{k[0]} cm{k[1]}"  # noqa: E731
    by_scen = {p: {name(k): v[p] for k, v in extra.items()} for p in OTHERS}
    return {"void": False, "why": [], "reproduction_largest_difference": worst,
            "arm_by_scenario": {a: {name(k): v for k, v in t.items()} for a, t in arm.items()},
            "pair_by_scenario": {p: {name(k): v for k, v in t.items()} for p, t in pair.items()},
            "pair_minus_built_by_scenario": by_scen,
            "expected_pair": exp_pair, "expected_arm": exp_arm, "pair_minus_built": lead,
            "smallest_lead_over_nine_priors": {p: min(lead[p].values()) for p in OTHERS},
            "largest_loss_by_scenario": {p: min(v["mean"] for v in by_scen[p].values()) for p in OTHERS},
            "largest_gain_by_scenario": {p: max(v["mean"] for v in by_scen[p].values()) for p in OTHERS},
            "main_prior": main, "check_priors": checks, "qualifies": qualifies, "built_confirmed": confirmed,
            "recommendation": rec, "weight_at_or_above_95_flat_cm": high,
            "new_arms": {name(k): v["new_arms"] for k, v in extra.items()},
            "commit": next(iter(runs.values()))["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    P = list(PAIRS)
    lines = [f"Try T33 (commit {res['commit']}): {res['recommendation']}",
             f"rerun arm against T30's stored scores: largest difference {res['reproduction_largest_difference']:.4f} "
             f"points (tolerance {REPRO_TOL})",
             f"expected better-of-two score by prior (main: {res['main_prior']}; checks: "
             f"{', '.join(res['check_priors'])}); in brackets the lead over the built pair:",
             "  prior (depth x cm)      " + "  ".join(f"{p:>17}" for p in P)]
    for k in res["expected_pair"]["built"]:
        cells = [f"{res['expected_pair']['built'][k]:17.2f}"]
        cells += [f"{res['expected_pair'][p][k]:8.2f} ({res['pair_minus_built'][p][k]:+5.2f})" for p in OTHERS]
        lines.append(f"  {k:<22}  " + "  ".join(cells))
    lines.append("qualifies under the rule: " + ", ".join(f"{p} {res['qualifies'][p]}" for p in OTHERS)
                 + f"; built confirmed {res['built_confirmed']}")
    lines.append("smallest lead over the nine priors / largest scenario loss / largest scenario gain: " + "; ".join(
        f"{p} {res['smallest_lead_over_nine_priors'][p]:+.2f} / {res['largest_loss_by_scenario'][p]:+.2f} / "
        f"{res['largest_gain_by_scenario'][p]:+.2f}" for p in OTHERS))
    arms = list(res["arm_by_scenario"])
    lines.append("mean score by scenario (truth depth relative to library 3, truth cardiomyocyte share):")
    lines.append("  scenario        " + "  ".join(f"{a:>9}" for a in arms) + "  | " + "  ".join(f"{p:>9}" for p in P))
    for k in res["arm_by_scenario"][arms[0]]:
        lines.append(f"  {k:<14}  " + "  ".join(f"{res['arm_by_scenario'][a][k]:9.2f}" for a in arms) + "  | "
                     + "  ".join(f"{res['pair_by_scenario'][p][k]:9.2f}" for p in P))
    lines.append("pair minus built by scenario, mean (seeds ahead/behind):")
    for k in res["arm_by_scenario"][arms[0]]:
        lines.append(f"  {k:<14}  " + "  ".join(
            f"{p} {res['pair_minus_built_by_scenario'][p][k]['mean']:+6.2f} "
            f"({res['pair_minus_built_by_scenario'][p][k]['ahead']}/{res['pair_minus_built_by_scenario'][p][k]['behind']})"
            for p in OTHERS))
    lines.append("reported, deciding nothing:")
    lines.append("  weight of scenarios with a mean better-of-two of 95 or more (flat cell-mix prior): " + "; ".join(
        f"{dp}: " + " / ".join(f"{res['weight_at_or_above_95_flat_cm'][p][dp]:.2f}" for p in P)
        for dp in t30.DEPTH_PRIORS) + f"  (order: {' / '.join(P)})")
    first = next(iter(res["new_arms"].values()))
    lines.append("  new arms at the first scenario, ten-seed means: " + "; ".join(
        f"{a}: repeats {v['repeat_share']:.2f}, median library {v['median_library']:.0f}, library 3's share of draws "
        f"{v['lib3_draw_share']:.2f}" for a, v in first.items()))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    res = reading(_load(a.runs, "t33"), _load(a.t30, "t30"), _load(a.t32, "t32"))
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
    d.add_argument("--t30", default=str(T30_SCENARIOS))
    d.add_argument("--t32", default=str(T32_SCENARIOS))
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
