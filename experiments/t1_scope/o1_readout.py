#!/usr/bin/env python
"""O1's reproduction gate and its readout (tries T8 to T10, experiments/harness/TRIES.md; the arms are built by
experiments/t1_scope/pseudo_task.py --o1-arms). A diagnostic: it decides only whether the within-type (G) line gets its
GPU, its days and the X002/X005 rulings; it never enters an upload.

Two steps, in order, each with the numbers fixed in the row:

  repro   Reproduce the committed E8.5 -> E9.5 runs before anything new is read: per truth library and seed, the arms
          copy_last, oracle and shift_1 of a fresh run (--only oracle,shift_1 --n-cells 5000 --strict-outputs, under a
          new tag) against the committed scores.csv. Pass: every |difference of score| <= 0.02 points and every
          |difference of a ranked metric's skill| <= 0.002. A failure stops the O1 run: it does not measure what the row
          says it measures, and that goes to the captain.
  decide  The bar, on the two truths' O1 runs (n_cells 2,500, seeds 0-9, both from one clean commit, repro passed):
          points of mmd_u + variogram = 100 x (0.30 x skill:mmd_u + 0.20 x skill:variogram), the arm minus shift_1.
          state_transport_1.0 clears the absolute bar if its mean over both truths and the 10 seeds is at least +1.5
          points and the per-seed difference (the mean over the two truths of each seed's difference) is above 0 in at
          least 8 of 10 seeds. Validity: state_transport_1.0's de_score and de_direction skills within 0.02 of shift_1's
          (mean over seeds, each truth), else the per-type pseudobulk was not held and the reading is void. Reported
          and deciding nothing: P (sibling_real minus shift_1, points), R (state_transport_1.0 minus shift_1, over P),
          the b0-like arm against shift_1, every arm's skills and points against the floor (skill 0.5), and cells
          reused per category.

Outputs go to experiments/t1_scope/results/<--out>/ (o1_repro.json and .txt, or o1_decision.json and .txt), each stamped
with the commit, the clean-tree flag (vec.provenance), the UTC time and the sha256 of every input file. The step refuses
a folder that already holds its outputs, and it reads only the run folders it is given (no data file, no board file).

A third step, gmom, reads try T11's two pseudo-task results (the harness's per-task JSON files of a run with b0 and
g_mom, experiments/harness/run.py) and writes what the adoption rule does not: per counted cut, the mean over seeds of
each scored group's margin of g_mom over b0 (100 x weight x difference of the group skill), which sum to the score
margin, and whether at least two thirds of the gain sits on the distribution and covariance groups (rule 2, the row
for T11). It decides nothing alone: T11 is kept only if the adoption rule keeps it on both cuts and this holds on both.

Usage (repo root, after the runs of the row; numbers only, runs on any machine with the results):
  python -I -B experiments/t1_scope/o1_readout.py repro --lib1 o1_repro_lib1 --lib2 o1_repro_lib2
  python -I -B experiments/t1_scope/o1_readout.py decide --lib1 o1_E8.5_E9.5_lib1 --lib2 o1_E8.5_E9.5_lib2
  python -I -B experiments/t1_scope/o1_readout.py gmom --runs runs/harness/t11_g_mom --out t11_g_mom_groups
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import pandas as pd

from vec import provenance
from vec.boards import GROUPS

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
SEEDS = list(range(10))
COMMITTED = {"1": "pseudo_E8.5_E9.5_lib1", "2": "pseudo_E8.5_E9.5_lib2"}
REPRO_ARMS = ("copy_last", "oracle", "shift_1")
RANKED_SKILLS = ("skill:de_score", "skill:de_direction", "skill:mmd_u", "skill:variogram")
REPRO_SCORE_TOL = 0.02  # points, per seed and arm
REPRO_SKILL_TOL = 0.002  # skill, per seed, arm and ranked metric
REPRO_N_CELLS = 5000  # the committed runs'
O1_N_CELLS = 2500
BAR_POINTS = 1.5
BAR_SEEDS = 8  # of 10
DE_SKILL_TOL = 0.02
T11_CUTS = {"T1/E9.5->E10.5": "T1_E9.5_E10.5.json", "T1/E9.5->E10.5:chamber": "T1_E9.5_E10.5_chamber.json"}
T11_SHAPE_SHARE = 2.0 / 3.0  # of the gain on the distribution and covariance groups, on each cut (rule 2, T11's row)
T11_BASE, T11_CHANGE = "b0", "g_mom"
BASE = "shift_1"
TRANSPORT, TRANSPORT_B0LIKE, SIBLING = "state_transport_1.0", "state_transport_b0like_1.0", "sibling_real"
ARMS = ("copy_last", "oracle", BASE, SIBLING, TRANSPORT, TRANSPORT_B0LIKE)
_WEIGHT = {m: w for _, w, members in GROUPS["T1"] for m, _ in members}
MMD_WEIGHT, VARIOGRAM_WEIGHT = _WEIGHT["mmd_u"], _WEIGHT["variogram"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stamp(inputs: list[Path]) -> dict:
    """The stamp every output carries: commit, the clean-tree flag (vec.provenance's rule), the UTC time, the inputs."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    return {"commit": commit, "tree_clean": unclean == [], "unclean_paths": unclean,
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "inputs": {str(p.relative_to(provenance.REPO_ROOT)) if p.is_relative_to(provenance.REPO_ROOT) else str(p):
                       sha256_file(p) for p in inputs}}


def refuse_existing(out: Path, names: tuple[str, ...]) -> None:
    """Refuses a folder that already holds one of this step's outputs."""
    if any((out / n).exists() for n in names):
        sys.exit(f"{out} already holds this step's outputs; a rerun needs a new --out and a new TRIES.md row")


def load_run(tag: str) -> tuple[dict, pd.DataFrame, list[Path]]:
    d = RESULTS / tag
    paths = [d / "run.json", d / "scores.csv"]
    if not all(p.exists() for p in paths):
        sys.exit(f"{d}: no run.json and scores.csv")
    return json.loads(paths[0].read_text()), pd.read_csv(paths[1]), paths


def check_seeds(df: pd.DataFrame, arms: tuple[str, ...], where: str) -> None:
    """Every arm present for each of seeds 0-9, once."""
    for arm in arms:
        got = sorted(df.loc[df["candidate"] == arm, "seed"].tolist())
        if got != SEEDS:
            sys.exit(f"{where}: arm {arm} has seeds {got}, not 0-9")


def compare_runs(new: pd.DataFrame, committed: pd.DataFrame, arms: tuple[str, ...] = REPRO_ARMS) -> dict:
    """Per arm and metric, the largest |difference| over the seeds between a fresh run and a committed run."""
    out = {}
    for arm in arms:
        a = new[new["candidate"] == arm].set_index("seed")
        b = committed[committed["candidate"] == arm].set_index("seed")
        out[arm] = {"score": float((a["score"] - b["score"]).abs().max()),
                    **{m: float((a[m] - b[m]).abs().max()) for m in RANKED_SKILLS}}
    return out


def repro_verdict(devs: dict) -> tuple[bool, list[str]]:
    """Pass if every deviation of every arm is within its tolerance."""
    problems = []
    for lib, per_arm in devs.items():
        for arm, d in per_arm.items():
            if d["score"] > REPRO_SCORE_TOL:
                problems.append(f"library {lib} {arm}: score differs by {d['score']:.4f} > {REPRO_SCORE_TOL}")
            for m in RANKED_SKILLS:
                if d[m] > REPRO_SKILL_TOL:
                    problems.append(f"library {lib} {arm}: {m} differs by {d[m]:.5f} > {REPRO_SKILL_TOL}")
    return not problems, problems


def cmd_repro(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("o1_repro.json", "o1_repro.txt"))
    res, inputs, runs = {}, [], {}
    for lib, tag in (("1", a.lib1), ("2", a.lib2)):
        run, new, paths = load_run(tag)
        committed_run, committed, cpaths = load_run(COMMITTED[lib])
        inputs += [*paths, *cpaths]
        if run["n_cells"] != REPRO_N_CELLS or committed_run["n_cells"] != REPRO_N_CELLS:
            sys.exit(f"{tag}: n_cells {run['n_cells']} (committed {committed_run['n_cells']}), not {REPRO_N_CELLS}")
        if run["truth"]["library"] != lib or list(run["seeds"]) != SEEDS or run.get("frac") != 0.1:
            sys.exit(f"{tag}: not truth library {lib}, seeds 0-9 and frac 0.1")
        check_seeds(new, REPRO_ARMS, tag)
        check_seeds(committed, REPRO_ARMS, COMMITTED[lib])
        res[lib] = compare_runs(new, committed)
        prov = run.get("provenance", {})
        runs[lib] = {"tag": tag, "committed": COMMITTED[lib], "commit_of_new_run": prov.get("commit"),
                     "tree_clean_of_new_run": prov.get("tree_clean")}
    ok, problems = repro_verdict(res)
    record = {"verdict": "pass" if ok else "fail", "tolerances": {"score_points": REPRO_SCORE_TOL,
                                                                  "ranked_skill": REPRO_SKILL_TOL},
              "max_abs_differences": res, "problems": problems, "runs": runs, "stamp": stamp(inputs)}
    out.mkdir(parents=True, exist_ok=True)
    (out / "o1_repro.json").write_text(json.dumps(record, indent=2) + "\n")
    lines = [f"O1 reproduction gate: {record['verdict']} (tolerances: score {REPRO_SCORE_TOL} points, skill "
             f"{REPRO_SKILL_TOL} per seed)"]
    for lib, per_arm in res.items():
        for arm, d in per_arm.items():
            lines.append(f"  library {lib} {arm}: max |d score| {d['score']:.4f}; "
                         + ", ".join(f"{m.split(':')[1]} {d[m]:.5f}" for m in RANKED_SKILLS))
    lines += [f"  PROBLEM: {p}" for p in problems]
    lines.append(f"Commit {record['stamp']['commit']}, tree clean {record['stamp']['tree_clean']}, UTC "
                 f"{record['stamp']['utc']}")
    text = "\n".join(lines) + "\n"
    (out / "o1_repro.txt").write_text(text)
    print(text, end="")
    if not ok:
        sys.exit("the reproduction gate failed: do not run or read the O1 arms; the row sends this to the captain")


def add_points(df: pd.DataFrame) -> pd.DataFrame:
    """Board points of mmd_u + variogram: 100 x (weight x skill) for each, as the board totals them."""
    return df.assign(points=100.0 * (MMD_WEIGHT * df["skill:mmd_u"] + VARIOGRAM_WEIGHT * df["skill:variogram"]))


def per_seed(df: pd.DataFrame, value: str) -> pd.DataFrame:
    return df.pivot(index="seed", columns="candidate", values=value).loc[SEEDS]


def read_truth(df: pd.DataFrame) -> dict:
    """One truth's numbers: mean over seeds of each arm's skills and points, and of the differences that matter."""
    pts = per_seed(df, "points")
    d_tr = pts[TRANSPORT] - pts[BASE]
    d_p = pts[SIBLING] - pts[BASE]
    out = {"seed_differences_transport_minus_shift_1": {int(s): float(v) for s, v in d_tr.items()},
           "mean_transport_minus_shift_1": float(d_tr.mean()), "seeds_ahead": int((d_tr > 0).sum()),
           "prize_P": float(d_p.mean()), "reach_R": float(d_tr.mean() / d_p.mean()) if d_p.mean() > 0 else float("nan")}
    arms = {}
    for arm in ARMS:
        row = {m.split(":")[1]: float(per_seed(df, m)[arm].mean()) for m in RANKED_SKILLS}
        row["points_mmd_u_variogram"] = float(pts[arm].mean())
        floor_points = 100.0 * 0.5 * (MMD_WEIGHT + VARIOGRAM_WEIGHT)  # both skills at the floor's 0.5
        row["points_above_floor_mmd_u_variogram"] = float(pts[arm].mean() - floor_points)
        arms[arm] = row
    out["arms"] = arms
    de_gap = {}
    for arm in (TRANSPORT, TRANSPORT_B0LIKE):
        de_gap[arm] = {m: float(per_seed(df, f"skill:{m}")[arm].mean() - per_seed(df, f"skill:{m}")[BASE].mean())
                       for m in ("de_score", "de_direction")}
    out["de_skill_gap_vs_shift_1"] = de_gap
    b0like = pts[TRANSPORT_B0LIKE]
    out["b0like"] = {"vs_shift_1_points": float((b0like - pts[BASE]).mean()),
                     "vs_transport_points": float((b0like - pts[TRANSPORT]).mean())}
    return out


def decide(lib1: dict, lib2: dict, seed_diffs: tuple[pd.Series, pd.Series]) -> dict:
    """The row's bar on two truths' numbers (read_truth) and their per-seed differences of the transport arm."""
    d1, d2 = seed_diffs
    both = (d1 + d2) / 2.0
    mean_gain = float(both.mean())
    ahead = int((both > 0).sum())
    de_ok = all(abs(t["de_skill_gap_vs_shift_1"][TRANSPORT][m]) <= DE_SKILL_TOL for t in (lib1, lib2)
                for m in ("de_score", "de_direction"))
    bar = mean_gain >= BAR_POINTS and ahead >= BAR_SEEDS
    p = (lib1["prize_P"] + lib2["prize_P"]) / 2.0
    verdict = "void" if not de_ok else "continue" if bar else "stop"
    return {"mean_transport_minus_shift_1": mean_gain, "seeds_ahead_of_10": ahead, "bar_points": BAR_POINTS,
            "bar_seeds": BAR_SEEDS, "bar_cleared": bar, "de_within_tolerance": de_ok, "de_skill_tol": DE_SKILL_TOL,
            "prize_P_mean_of_truths": p, "reach_R_ratio_of_means": mean_gain / p if p > 0 else float("nan"),
            "verdict": verdict}


def cmd_decide(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("o1_decision.json", "o1_decision.txt"))
    repro_path = RESULTS / a.repro / "o1_repro.json"
    if not repro_path.exists() or json.loads(repro_path.read_text()).get("verdict") != "pass":
        sys.exit(f"{repro_path}: no passed reproduction gate; the row runs it first")
    inputs: list[Path] = [repro_path]
    truths, runs, diffs = {}, {}, {}
    for lib, tag in (("1", a.lib1), ("2", a.lib2)):
        run, df, paths = load_run(tag)
        reuse_path = RESULTS / tag / "o1_reuse.json"
        inputs += [*paths, reuse_path]
        if run["n_cells"] != O1_N_CELLS or run["truth"]["library"] != lib or list(run["seeds"]) != SEEDS:
            sys.exit(f"{tag}: not truth library {lib} at n_cells {O1_N_CELLS}, seeds 0-9")
        prov = run.get("provenance") or {}
        if not run.get("o1") or not prov.get("tree_clean"):
            sys.exit(f"{tag}: not an O1 run from a clean tree (provenance: {prov})")
        check_seeds(df, ARMS, tag)
        df = add_points(df)
        truths[lib] = read_truth(df)
        pts = per_seed(df, "points")
        diffs[lib] = pts[TRANSPORT] - pts[BASE]
        runs[lib] = {"tag": tag, "commit": prov.get("commit"), "tree_clean": prov.get("tree_clean"),
                     "cells_reused_per_seed": json.loads(reuse_path.read_text())["cells_reused_per_seed"]}
    if runs["1"]["commit"] != runs["2"]["commit"]:
        sys.exit("the two truths' O1 runs come from different commits")
    dec = decide(truths["1"], truths["2"], (diffs["1"], diffs["2"]))
    record = {"decision": dec, "truths": truths, "runs": runs, "stamp": stamp(inputs),
              "row": "experiments/harness/TRIES.md, tries T8 to T10"}
    out.mkdir(parents=True, exist_ok=True)
    (out / "o1_decision.json").write_text(json.dumps(record, indent=2) + "\n")
    text = report_text(record)
    (out / "o1_decision.txt").write_text(text)
    print(text, end="")


def report_text(record: dict) -> str:
    dec = record["decision"]
    lines = [f"O1 readout: verdict {dec['verdict']} (bar {dec['bar_points']} points and {dec['bar_seeds']}/10 seeds; "
             f"mean gain {dec['mean_transport_minus_shift_1']:+.3f}, ahead in {dec['seeds_ahead_of_10']}/10; "
             f"de_score and de_direction within {dec['de_skill_tol']}: {dec['de_within_tolerance']})",
             f"Descriptive: P {dec['prize_P_mean_of_truths']:+.3f} points (mean of both truths), R "
             f"{dec['reach_R_ratio_of_means']:.3f}"]
    for lib, t in record["truths"].items():
        b0like = t["b0like"]
        lines.append(f"Truth library {lib}: transport - shift_1 {t['mean_transport_minus_shift_1']:+.3f} points "
                     f"(ahead {t['seeds_ahead']}/10), P {t['prize_P']:+.3f}, R {t['reach_R']:.3f}; b0-like arm vs "
                     f"shift_1 {b0like['vs_shift_1_points']:+.3f}, vs transport {b0like['vs_transport_points']:+.3f}")
        for arm, r in t["arms"].items():
            lines.append(f"    {arm:29s} mmd_u {r['mmd_u']:.3f}, variogram {r['variogram']:.3f}, de_score "
                         f"{r['de_score']:.3f}, de_direction {r['de_direction']:.3f}; points (mmd_u + variogram) "
                         f"{r['points_mmd_u_variogram']:.2f}, above floor "
                         f"{r['points_above_floor_mmd_u_variogram']:+.2f}")
        for arm, g in t["de_skill_gap_vs_shift_1"].items():
            lines.append(f"    {arm}: de_score {g['de_score']:+.4f}, de_direction {g['de_direction']:+.4f} skill vs "
                         "shift_1")
        for arm, per in record["runs"][lib]["cells_reused_per_seed"].items():
            worst = {}
            for seed_reuse in per.values():
                for cat, n in seed_reuse.items():
                    worst[cat] = max(worst.get(cat, 0), int(n))
            shown = ", ".join(f"{c} {n}" for c, n in sorted(worst.items()) if n > 0) or "none"
            lines.append(f"    cells reused, worst seed, {arm}: {shown}")
    s = record["stamp"]
    lines.append(f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}")
    return "\n".join(lines) + "\n"


def group_margins(result: dict, base: str = T11_BASE, change: str = T11_CHANGE) -> dict:
    """One pseudo-task's result (`vec.harness.run_task`'s JSON): the mean over seeds of each scored group's margin of
    `change` over `base` in board points (100 x weight x the difference of group skills, as `harness.t3_fit_checks`
    reads T3's), their sum, the score margin the harness scored, and the share of the sum on the distribution and
    covariance groups (None when the sum is not above 0, which `decide` would not keep either)."""
    groups = [(g, w) for g, w, _ in GROUPS["T1"]]
    per_seed = []
    for row in result["rows"]:
        pair = row["groups"][change], row["groups"][base]
        per_seed.append({g: 100.0 * w * (pair[0][g] - pair[1][g]) for g, w in groups})
    margins = {g: float(pd.Series([p[g] for p in per_seed]).mean()) for g, _ in groups}
    total = float(sum(margins.values()))
    score = float(pd.Series([row["scores"][change] - row["scores"][base] for row in result["rows"]]).mean())
    shape = margins["distribution"] + margins["covariance"]
    return {"margins": margins, "sum_of_margins": total, "score_margin": score, "n_seeds": len(result["rows"]),
            "share_on_distribution_and_covariance": shape / total if total > 0 else None}


def cmd_gmom(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("gmom_groups.json", "gmom_groups.txt"))
    runs = Path(a.runs)
    cuts, inputs = {}, []
    for key, fname in T11_CUTS.items():
        path = runs / fname
        if not path.exists():
            sys.exit(f"{path}: no result for {key}; run experiments/harness/run.py with --candidates floor,b0,g_mom")
        result = json.loads(path.read_text())
        inputs.append(path)
        missing = [c for c in (T11_BASE, T11_CHANGE) if c not in result["candidates"]]
        if missing or result["task"]["key"] != key:
            sys.exit(f"{path}: not {key} with {T11_BASE} and {T11_CHANGE}")
        cuts[key] = group_margins(result)
    shares = [c["share_on_distribution_and_covariance"] for c in cuts.values()]
    ok = all(s is not None and s >= T11_SHAPE_SHARE for s in shares)
    record = {"cuts": cuts, "shape_share_required": T11_SHAPE_SHARE, "rule2_gain_on_predicted_groups": ok,
              "row": "experiments/harness/TRIES.md, try T11", "stamp": stamp(inputs)}
    out.mkdir(parents=True, exist_ok=True)
    (out / "gmom_groups.json").write_text(json.dumps(record, indent=2, default=float) + "\n")
    lines = [f"g_mom against b0, group margins per cut (rule 2: at least {T11_SHAPE_SHARE:.0%} of the gain on "
             f"distribution and covariance, on both cuts): {ok}"]
    for key, c in cuts.items():
        share = c["share_on_distribution_and_covariance"]
        lines.append(f"  {key}: " + ", ".join(f"{g} {m:+.2f}" for g, m in c["margins"].items())
                     + f"; sum {c['sum_of_margins']:+.2f} (scored margin {c['score_margin']:+.2f}); share on "
                     f"distribution and covariance {'n/a' if share is None else f'{share:.2f}'}")
    st = record["stamp"]
    lines.append(f"Commit {st['commit']}, tree clean {st['tree_clean']}, UTC {st['utc']}")
    text = "\n".join(lines) + "\n"
    (out / "gmom_groups.txt").write_text(text)
    print(text, end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("repro")
    r.add_argument("--lib1", required=True, help="tag of the fresh run with truth library 1 (n_cells 5000)")
    r.add_argument("--lib2", required=True, help="tag of the fresh run with truth library 2 (n_cells 5000)")
    r.add_argument("--out", default="o1_repro")
    d = sub.add_parser("decide")
    d.add_argument("--lib1", required=True, help="tag of the O1 run with truth library 1 (n_cells 2500)")
    d.add_argument("--lib2", required=True, help="tag of the O1 run with truth library 2 (n_cells 2500)")
    d.add_argument("--repro", default="o1_repro", help="folder of the reproduction gate's outputs")
    d.add_argument("--out", default="o1_decision")
    g = sub.add_parser("gmom")
    g.add_argument("--runs", required=True, help="the harness run folder with the two T1 cuts' JSON files (try T11)")
    g.add_argument("--out", default="t11_g_mom_groups")
    a = ap.parse_args(argv)
    {"repro": cmd_repro, "decide": cmd_decide, "gmom": cmd_gmom}[a.cmd](a)


if __name__ == "__main__":
    main()
