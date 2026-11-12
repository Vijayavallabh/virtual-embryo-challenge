#!/usr/bin/env python
"""Try T30 (experiments/harness/TRIES.md, registered before the run): which PAIR of T1 uploads pays most, priced on one
grid inside E10.5 at zero time gap. A labelled diagnostic for a ruling; no upload changes and no cell's depth changes
in any arm.

Why. A cold review of 2026-11-12 found the T1 half of the campaign's conclusions unpriced: no arm compares a second
upload that covers the base's library against today's second upload, which covers the cardiomyocyte share. And T29
probed only truths SHALLOWER than library 3, while four of the seven released libraries are deeper than it (median
counts per chamber cell: E8.5 14,755 and 13,990; E9.5 14,656 and 13,704; E10.5 4,030, 6,148 and 10,503).

Geometry: T27's (experiments/t1_scope/kind_gap.py, imported unchanged). The truth is one half of E10.5's library 3;
the reference is every other E10.5 cell; arms draw real cells outside the truth half.

Scenarios: 4 truth depths x 3 truth cardiomyocyte shares.
  depth 1.0      the truth half as released (T27's truth)
  depth 0.68, 0.35   its counts thinned binomially (T29's truths, the same draw)
  depth 1.37     a truth deeper than library 3, as the E8.5 and E9.5 libraries are: the truth stays as released and
                 every OTHER E10.5 cell (the reference, hence the floor, and every arm's cells) is thinned to 1/1.37,
                 so the truth is 1.37 times deeper relative to all of them. This thins a scenario, not a candidate:
                 all arms and the floor are thinned alike.
  cm "own"       the truth half's own mix (cardiomyocytes 0.25 of heart cells)
  cm 0.45, 0.60  the truth half re-drawn by quota so that cardiomyocytes are that share of heart cells, blood and the
                 non-chamber cells keeping their shares (vec.t1_scope.with_share, the rule the s2 build uses)

Arms (n cells by the count rule; real cells; none from the truth half):
  copy        a uniform draw of the reference pool
  s1_base     today's s1 base: chamber categories at the pooled E10.5 mix, libraries in the file's proportions (T27's)
  s1_cm60     today's s2 rule on that base: the same draw rule at cardiomyocytes 0.60 of heart cells
  lib3_base   T27's library-3 base: library 3's cells at library 3's own chamber mix

Pairs: P(A, B) on a scenario is the mean over seeds of the better of the two. today = (s1_base, s1_cm60);
base hedge = (s1_base, lib3_base); switched = (lib3_base, s1_cm60).

Priors, fixed before the run. Depth: "census" gives each released library one vote at its nearest grid depth
(0.35: 1/7, 0.68: 1/7, 1.0: 1/7, 1.37: 4/7); "latest" uses E10.5's three libraries only (1/3, 1/3, 1/3, 0); "flat" is
1/4 each. Cardiomyocyte share: "flat" (1/3 each); "like" (own 0.60, 0.45: 0.25, 0.60: 0.15: library 3 is the
like-for-like sample); "rich" (0, 1/2, 1/2: the s2 rule's own case). The joint prior is the product.

Rule, fixed before the run. Main prior: census x flat. Checks: latest x flat, and census x like.
  - A pair other than today's is put to the captain only if its expected score exceeds today's pair by more than 1.0
    under the main prior and is not below today's under both checks. If two qualify, the higher under the main prior.
  - Today's pair is confirmed if it beats both others by more than 1.0 under the main prior and is not below either
    under both checks.
  - Otherwise: no recommendation; the captain is shown the table.
All nine prior products are reported. Zero time gap: sizes are upper bounds for the test, where two days of
development lie between base and truth; what is read is which pair leads. The reading is void if a scenario, seed or
arm is missing, or results come from different commits or a dirty tree.

  run   --depth F --cm own|0.45|0.60 --out DIR   one scenario (start the twelve together: one commit)
  read  --runs DIR --out DIR                     apply the rule; write pricing.json and pricing.txt
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kind_gap as kg  # noqa: E402  (T27's geometry, draws and helpers, unchanged)

from vec import harness, t1_scope  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import RANKED_METRICS, load_boards  # noqa: E402
from vec.core import quotas  # noqa: E402
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample  # noqa: E402

ARMS = ("copy", "s1_base", "s1_cm60", "lib3_base")
CM_HEDGE = 0.60  # docs/S2_RULES.md: T1's s2
CM60_STREAM = 10  # s1_cm60's draw stream, apart from kind_gap's arm indices
DEPTHS = (0.35, 0.68, 1.0, 1.37)
CMS = ("own", "0.45", "0.60")
TRUTH_SEED = 20261114  # the truth half's quota re-draw
DEEP_SEED = 20261115  # the thinning of every non-truth cell in the deeper-truth scenario
PAIRS = {"today": ("s1_base", "s1_cm60"), "base_hedge": ("s1_base", "lib3_base"), "switched": ("lib3_base", "s1_cm60")}
DEPTH_PRIORS = {"census": (1 / 7, 1 / 7, 1 / 7, 4 / 7), "latest": (1 / 3, 1 / 3, 1 / 3, 0.0), "flat": (0.25,) * 4}
CM_PRIORS = {"flat": (1 / 3, 1 / 3, 1 / 3), "like": (0.60, 0.25, 0.15), "rich": (0.0, 0.5, 0.5)}
MAIN = ("census", "flat")
CHECKS = (("latest", "flat"), ("census", "like"))
BAR = 1.0
_D: dict = {}


def cm_truth_rows(d: dict, cm: float) -> np.ndarray:
    """The truth half re-drawn so that cardiomyocytes are `cm` of its heart cells: the largest quota draw its
    categories allow, the non-chamber cells cut by the same factor as the chamber cells."""
    cat, half = d["cat"], d["half_a"]
    chamber = d["chamber"]
    counts = {c: int((cat[half] == c).sum()) for c in chamber}
    total = sum(counts.values())
    target = t1_scope.with_share({c: counts[c] / total for c in chamber}, "cm", cm)
    n_ch = int(min(counts[c] / target[c] for c in chamber if target[c] > 0))
    rng = np.random.default_rng([TRUTH_SEED, int(round(cm * 100))])
    keep = []
    for c in chamber:
        pool = half[cat[half] == c]
        keep.append(pool[rng.choice(len(pool), min(len(pool), int(n_ch * target[c])), replace=False)])
    rest = half[~np.isin(cat[half], chamber)]
    keep.append(rest[rng.choice(len(rest), int(round(len(rest) * n_ch / total)), replace=False)])
    return np.sort(np.concatenate(keep))


def prepare(data_dir: Path, depth: float, cm: str) -> dict:
    d = kg.prepare(data_dir, 1.0)  # T27's split, pools, mixes and board, the truth as released
    rows = d["half_a"] if cm == "own" else cm_truth_rows(d, float(cm))
    x_truth = d["X"][rows]
    if depth < 1.0:  # T29's thinning of the truth (the same generator; with cm "own" it is T29's truth)
        x_truth = kg.lognorm(kg.thin(d["C"][rows], depth, np.random.default_rng(kg.THIN_SEED)))
    elif depth > 1.0:  # a deeper truth: every other cell thinned alike, the truth as released
        other = np.setdiff1d(np.arange(d["n_obs"]), d["half_a"])
        Ct = kg.thin(d["C"][other], 1.0 / depth, np.random.default_rng(DEEP_SEED))
        X = sparse.vstack([sparse.csr_matrix(d["X"][d["half_a"]]), kg.lognorm(Ct)], format="csr")
        order = np.argsort(np.r_[d["half_a"], other], kind="stable")  # back to the file's row order
        d["X"] = X[order]
        tot = d["tot"].astype(float).copy()
        tot[other] = np.asarray(Ct.sum(1)).ravel()
        d["tot"] = tot
    boards = load_boards()
    d["board"] = replace(boards[kg.TWIN], key=f"pseudo:T1/E10.5 pair pricing d{depth} cm{cm}", split="pseudo",
                         anchors={}, **harness.pseudo_counts(len(rows), len(d["ref_rows"])))
    d.update(x_truth=x_truth, truth_rows=rows, truth_depth=float(depth), truth_cm=cm,
             cm60_mix=t1_scope.with_share(d["pooled_mix"], "cm", CM_HEDGE))
    heart = list(t1_scope.load_spec().heart)
    tc = {c: int((d["cat"][rows] == c).sum()) for c in d["chamber"]}
    d["truth_cm_of_heart"] = tc["cm"] / max(sum(tc[c] for c in heart), 1)
    d["truth_cells"] = int(len(rows))
    return d


def arm_rows(arm: str, n: int, seed: int) -> np.ndarray:
    if arm != "s1_cm60":
        return kg.arm_rows(arm, n, seed)  # kind_gap's draws, stream for stream
    d = kg._D
    rng = np.random.default_rng([seed, CM60_STREAM])
    q = quotas(d["cm60_mix"], n)
    out = []
    for c in d["chamber"]:
        k = int(q.get(c, 0))
        k3 = int(round(k * d["lib3_share"][c]))
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
    row = {"seed": seed, "scores": {}, "skills": {}, "median_library": {}}
    for arm in ARMS:
        rows = arm_rows(arm, n, harness.BUILD_SEED_OFFSET + seed)
        r = score_board(board, Stage(X=X[rows]), target, ref, seed=seed, anchors=f.anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        row["scores"][arm], row["skills"][arm] = r.score, r.skills
        row["median_library"][arm] = float(np.median(d["tot"][rows]))
    row["seconds"] = round(time.time() - t0, 1)
    print(f"seed {seed}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["scores"].items()), flush=True)
    return row


def cmd_run(a) -> None:
    if a.depth not in DEPTHS or a.cm not in CMS:
        sys.exit(f"--depth must be one of {DEPTHS} and --cm one of {CMS}")
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
    res = {"try": "T30", "commit": commit, "tree_clean": clean, "depth": a.depth, "cm": a.cm,
           "truth_cells": d["truth_cells"], "truth_cm_of_heart": d["truth_cm_of_heart"],
           "seeds": list(SEEDS), "arms": list(ARMS), "protocol": harness.PROTOCOL, "build": harness.BUILD,
           "frac": SUBSAMPLE_FRAC, "pooled_chamber_mix": d["pooled_mix"], "cm60_mix": d["cm60_mix"],
           "library3_chamber_mix": d["lib3_mix"], "rows": rows}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"t30_d{a.depth}_cm{a.cm}.json"
    path.write_text(json.dumps(res, indent=1))
    print(f"wrote {path}", flush=True)


def _expect(table: dict, dp: str, cp: str) -> float:
    return float(sum(DEPTH_PRIORS[dp][i] * CM_PRIORS[cp][j] * table[(d, c)]
                     for i, d in enumerate(DEPTHS) for j, c in enumerate(CMS)))


def reading(runs: dict) -> dict:
    why = []
    for d in DEPTHS:
        for c in CMS:
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
    if why:
        return {"void": True, "why": why, "recommendation": "none: the reading is void"}
    arm, pair = {a: {} for a in ARMS}, {p: {} for p in PAIRS}
    for key, r in runs.items():
        rows = sorted(r["rows"], key=lambda x: x["seed"])
        s = {a: np.array([x["scores"][a] for x in rows], float) for a in ARMS}
        for a in ARMS:
            arm[a][key] = float(s[a].mean())
        for p, (x, y) in PAIRS.items():
            pair[p][key] = float(np.maximum(s[x], s[y]).mean())
    priors = [(dp, cp) for dp in DEPTH_PRIORS for cp in CM_PRIORS]
    exp_pair = {p: {f"{dp} x {cp}": _expect(pair[p], dp, cp) for dp, cp in priors} for p in PAIRS}
    exp_arm = {a: {f"{dp} x {cp}": _expect(arm[a], dp, cp) for dp, cp in priors} for a in ARMS}
    main = " x ".join(MAIN)
    checks = [" x ".join(c) for c in CHECKS]
    others = [p for p in PAIRS if p != "today"]
    qualifies = {p: bool(exp_pair[p][main] - exp_pair["today"][main] > BAR
                         and all(exp_pair[p][c] >= exp_pair["today"][c] for c in checks)) for p in others}
    confirmed = bool(all(exp_pair["today"][main] - exp_pair[p][main] > BAR
                         and all(exp_pair["today"][c] >= exp_pair[p][c] for c in checks) for p in others))
    winners = [p for p in others if qualifies[p]]
    if winners:
        best = max(winners, key=lambda p: exp_pair[p][main])
        rec = f"put to the captain: the pair {best} {PAIRS[best]} in place of today's {PAIRS['today']}"
    elif confirmed:
        rec = "today's pair is confirmed"
    else:
        rec = "no recommendation: the pairs are within the priors' reach; the captain is shown the table"
    return {"void": False, "why": [], "arm_by_scenario": {a: {f"d{k[0]} cm{k[1]}": v for k, v in t.items()}
                                                          for a, t in arm.items()},
            "pair_by_scenario": {p: {f"d{k[0]} cm{k[1]}": v for k, v in t.items()} for p, t in pair.items()},
            "expected_pair": exp_pair, "expected_arm": exp_arm, "main_prior": main, "check_priors": checks,
            "qualifies": qualifies, "today_confirmed": confirmed, "recommendation": rec,
            "truth": {f"d{k[0]} cm{k[1]}": {"cells": r["truth_cells"], "cm_of_heart": r["truth_cm_of_heart"]}
                      for k, r in runs.items()},
            "commit": next(iter(runs.values()))["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Try T30 (commit {res['commit']}): {res['recommendation']}",
             f"expected better-of-two score by prior (main: {res['main_prior']}; checks: "
             f"{', '.join(res['check_priors'])}):",
             "  prior (depth x cm)      " + "  ".join(f"{p:>11}" for p in PAIRS) + "   | single arms: "
             + "  ".join(f"{a:>9}" for a in ARMS[1:])]
    for k in res["expected_pair"]["today"]:
        lines.append(f"  {k:<22}  " + "  ".join(f"{res['expected_pair'][p][k]:11.2f}" for p in PAIRS)
                     + "   |              " + "  ".join(f"{res['expected_arm'][a][k]:9.2f}" for a in ARMS[1:]))
    lines.append("mean score by scenario (truth depth relative to library 3, truth cardiomyocyte share):")
    lines.append("  scenario        " + "  ".join(f"{a:>9}" for a in ARMS) + "  | "
                 + "  ".join(f"{p:>11}" for p in PAIRS))
    for k in res["arm_by_scenario"]["copy"]:
        lines.append(f"  {k:<14}  " + "  ".join(f"{res['arm_by_scenario'][a][k]:9.2f}" for a in ARMS) + "  | "
                     + "  ".join(f"{res['pair_by_scenario'][p][k]:11.2f}" for p in PAIRS))
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    runs = {}
    for p in sorted(Path(a.runs).glob("t30_d*_cm*.json")):
        r = json.loads(p.read_text())
        runs[(float(r["depth"]), str(r["cm"]))] = r
    res = reading(runs)
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
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
