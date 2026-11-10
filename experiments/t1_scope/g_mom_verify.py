#!/usr/bin/env python
"""Verification of T1 candidate g_mom (try T11, experiments/harness/TRIES.md) before it goes to the captain: its placebo
controls (tries T20 and T21), its no-development control and pseudobulk-held arm (tries T23 and T24), its composition
reading (try T22) and a target-free look at its factors at the test pair (diagnostic N21). Three steps, every number
fixed in their rows before their runs. None of them changes an upload.

  placebo      Tries T20, T21, T23 and T24. Reads a harness run (experiments/harness/run.py run with floor,b0,g_mom,
               g_mom_flip,g_mom_perm,g_mom_null,g_mom_held on the two counted E10.5 cuts, the same draws and base as
               g_mom) and applies the rows' rules: the momentum reading survives only if g_mom beats EACH control
               (the sign flip T20, the permutation T21 and the no-development arm T23) by a paired mean gain above 0.3
               points and in at least 8 of 10 seeds on both cuts (the adoption rule, `vec.harness.decide`, with the
               control as the base; a loss larger than 0.3 on either cut fails it too), and the pseudobulk-held arm
               T24 keeps a gain over b0 by the same rule on both cuts (`decide` with b0 as the base; if not, the gain
               is a per-type mean effect, consistent with O1, and not within-type momentum; its validity gate, its
               de_score and de_direction skills within 0.02 of b0's on each cut, voids the reading otherwise). It also
               reads the spread hypothesis from the harness's own metric panels (per arm and cut: variance_ratio, mmd_u,
               variogram, energy_distance, pb_rel_err), reports T20's sign-odd spread reading and T21's validity gate as
               a spread control (its variance_ratio within 0.003 of g_mom's on both cuts; otherwise T21 says only that
               the assignment matters). A pass of T20 says the effect is odd in the factor sign and a pass of T21 that
               the assignment matters; neither is a spread control or evidence of momentum. Before reading anything it
               requires the run's b0 to reproduce try T11's committed per-seed scores within 0.02 points and its g_mom
               within 0.05 points on both cuts (g_mom carries the fit's PCA and k-means, whose floats can differ between
               runs; b0 does not); otherwise the reading is void and goes to the captain.
  composition  Try T22. The same comparison of g_mom against b0, both built at the truth's own chamber mix, the
               oracle-base construction of experiments/s2_hedge/payoff.py t1-grid (the chamber cut's categories
               without its surrounding cells, renormalised), which removes the composition mismatch of the counted
               cuts (B0 draws E9.5's chamber mix). b0 at that mix (B0's population,
               `vec.t1_state_transport.b0_population`) and g_mom's factors of the fold on top (the fit on the
               reference's own chamber categories, as g_mom's) with R1's held labels unedited, built by the harness's
               builders; scored by the harness's scoring against the chamber cut (the rule, seeds 0-9) and, on the
               same scoring, the whole E10.5 cut (reported, deciding nothing). Rule: the gain must stay at or above
               +0.3 points in at least 8 of 10 seeds on the chamber cut; otherwise the gain is not shown to survive
               the known composition confound (TRIES rule 2's third condition) and g_mom is not proposed. A fail is
               inconclusive on power (its row says why). Before anything is decided the b0 arm must reproduce the
               committed t1-grid's oracle-base row (experiments/s2_hedge/results/t1_grid/payoff.json) within 0.02 points
               in the mean over seeds. Reported beside the rule, deciding nothing: the same oracle-mix population with
               the pseudobulk held (g_mom_held, the T24 reading at the oracle mix; O1's own construction was a held
               edit at an oracle mix), and a companion that does not shrink with the matched share: g_mom and b0 at the
               prediction's own chamber mix (the harness's arms) scored against the chamber truth re-sampled to that
               mix (payoff.stratified, the t1-cm axis's re-sampling, the truth's surrounding share kept). The drawn
               oracle-mix population's median library size against E10.5's is reported (docs/DATA_POLICY.md 2026-11-04
               (f)(3)).
  n21          Diagnostic N21, target-free, the test pair. g_mom's factors fitted exactly as the test would fit them
               (start E9.5, reference E10.5, chamber categories; R1's held labels excluded at application), and the
               backtest's fit (start E8.5, reference E9.5) on the categories in both. Reported per category, deciding
               nothing: the size of the factors, the share of cells and genes touched (E10.5's chamber cells, shares
               only), each edited cell's expected depth (kept after application, checked against the keep-depth
               constant `vec.core.apply_factors` returns), the concordance of the factor direction with the backtest's
               fit (the absolute cosine of the leading right singular vectors of the two per-metacell factor matrices,
               and the Spearman correlation of the per-gene sizes) with split-half ceilings (one random halving of each
               stage's cells per category, fixed seed; the same fit on each half: half-data estimates, so they
               under-estimate the full-data ceiling) and a no-development reference (the same statistics between two
               same-stage halves' fits, E10.5's and E9.5's), the amplification share (the regression slope and R^2
               of the finished factors on the metacell's own deviation) and the concordance after projecting that
               component out, and the share of the factor size on R2's union genes. It reads the released E8.5 and
               E9.5 files, the unlocked E10.5 file, the held-flat list, the chamber config and R2's union file; no
               board, no target, no score, no built file; it prints and stores no cell count.

Outputs go to experiments/t1_scope/results/<--out>/ (a JSON and a text file per step), stamped with the commit, the
clean-tree flag (vec.provenance's rule, from the start of the run to the write), the UTC time and the sha256 of every
input. A step refuses a folder that already holds its outputs; a rerun needs a new row.

Usage (DGX, repo root, BLAS at most 16 threads, after the row is committed and the harness run has finished for
`placebo`):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
      uv run --frozen --no-sync python -I -B experiments/t1_scope/g_mom_verify.py placebo \
      --runs runs/harness/t20_t24_g_mom_controls --out t20_t24_controls
  ... g_mom_verify.py composition --out t22_composition     (try T22; it reads data/raw, not the harness run)
  ... g_mom_verify.py n21 --out n21_g_mom_test_pair         (diagnostic N21)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

from vec import harness, provenance, t1_scope
from vec import t1_state_transport as ts
from vec.adoption import MAX_LOSS, MIN_AHEAD, MIN_GAIN, SEEDS
from vec.boards import GROUPS, RANKED_METRICS, REPO_ROOT, load_boards
from vec.core import apply_factors, counts_of, mean_cp10k
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
for _d in (HERE, REPO_ROOT / "experiments" / "s2_hedge", REPO_ROOT / "experiments" / "t1_hygiene"):
    sys.path.insert(0, str(_d))
import hygiene  # noqa: E402
import o1_readout as o1  # noqa: E402
import payoff  # noqa: E402

CUTS = {"T1/E9.5->E10.5": "T1_E9.5_E10.5.json", "T1/E9.5->E10.5:chamber": "T1_E9.5_E10.5_chamber.json"}
TASK_WHOLE, TASK_CHAMBER = tuple(CUTS)
PLACEBOS = {"g_mom_flip": "T20", "g_mom_perm": "T21", "g_mom_null": "T23"}  # candidate -> try; g_mom must beat each
HELD_ARM, HELD_TRY = "g_mom_held", "T24"  # the pseudobulk-held arm; it must keep a gain over b0
T11_TASKS = REPO_ROOT / "experiments" / "harness" / "results" / "t11_g_mom" / "tasks"
T1_GRID = REPO_ROOT / "experiments" / "s2_hedge" / "results" / "t1_grid" / "payoff.json"
ORACLE_ROW = "bown_cown"  # the grid's oracle-base prediction and truth: the chamber cut's own mix
REPRO_T11_TOL = {"b0": 0.02, "g_mom": 0.05}  # points per seed, T20 and T21's gate against try T11's committed scores
REPRO_SCORE_TOL = 0.02  # points, the mean over seeds, T22's gate against the committed t1-grid
T22_MIN_GAIN = MIN_GAIN  # 0.3 points, the rule's bar, per seed
T22_MIN_SEEDS = MIN_AHEAD  # 8 of the 10 seeds at or above it
SPREAD_METRICS = ("variance_ratio", "mmd_u", "variogram", "energy_distance", "pb_rel_err")
PANEL_METRICS = (*SPREAD_METRICS, "de_score", "de_direction", "library_size_ratio", "pseudobulk_pearson")
VR_TOL = 0.003  # T21's validity gate and T20's sign-odd reading: about 6 seed sds of the variance ratio (T11: 0.0005)
HELD_DE_TOL = 0.02  # T24's validity: its de_score and de_direction skills within this of b0's, as T9's gate
LIBRARY_SHIFT_TOL = 0.05  # docs/DATA_POLICY.md 2026-11-04 (f)(3): a larger shift goes to the captain
SPLIT_SEED = 20261110  # N21's halves: one generator per stage pair, default_rng([SPLIT_SEED, pair]) (test 0, back 1)
NO_DEV_STREAMS = {"test": 2, "backtest": 3}  # N21's no-development reference: the halves of E10.5 and of E9.5
N21_APPLY_STREAM = (SPLIT_SEED, 5)  # the one realised application N21 reports (depth only)
DEPTH_TOL = 1e-9  # the largest |expected depth after the edit / the cell's depth - 1| of an edited cell
TRUTH_RESAMPLE_SEED = 20261110  # T22's companion: the chamber truth re-sampled once to the prediction's own mix
LABEL_T22 = ("try T22: g_mom against b0 at the truth's own chamber mix (an oracle base, the construction of "
             "experiments/s2_hedge/payoff.py t1-grid); not an adoption result, it can only veto the proposal")
LABEL_N21 = ("diagnostic N21: g_mom's factors at the test pair (start E9.5, reference E10.5), against the backtest's "
             "(start E8.5, reference E9.5); target-free: no target, no score, no board, no built file, no cell count")


def _rel(path: Path) -> str:
    p = Path(path).resolve()
    return p.relative_to(REPO_ROOT).as_posix() if p.is_relative_to(REPO_ROOT) else str(path)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_state() -> dict:
    """HEAD and whether the tree is clean by vec.provenance's rule (None when git cannot tell), as N12 reads it."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    return {"commit": commit, "clean": None if unclean is None else not unclean}


def stamp(started: Mapping, inputs: Sequence[Path]) -> dict:
    """The stamp of a step's outputs: the commit, whether the tree was clean from the start of the run to the write and
    at one commit (vec.provenance's rule), the UTC time and the sha256 of every input file."""
    ended = tree_state()
    clean = (None if started["clean"] is None or ended["clean"] is None
             else started["clean"] and ended["clean"] and started["commit"] == ended["commit"])
    return {"commit": ended["commit"], "tree_clean": clean, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "inputs": {_rel(p): sha256_file(Path(p)) for p in inputs}}


def refuse_existing(out: Path, names: Sequence[str]) -> None:
    if any((out / n).exists() for n in names):
        sys.exit(f"{out} already holds this step's outputs; a rerun needs a new --out and a new TRIES.md row")


def parse_seeds(text: str) -> tuple[int, ...]:
    out: list[int] = []
    for part in text.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return tuple(out)


def write_outputs(out: Path, stem: str, record: dict, text: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{stem}.json").write_text(json.dumps(record, indent=1, default=float) + "\n")
    (out / f"{stem}.txt").write_text(text)
    print(text, end="")


# ---------------------------------------------------------------------------------------------------------------------
# Tries T20 and T21: the placebo controls


def load_run(runs: Path) -> dict[str, dict]:
    """The harness run's two per-cut results (key -> run_task JSON), with b0, g_mom, the controls and the held arm."""
    out = {}
    for key, fname in CUTS.items():
        path = runs / fname
        if not path.exists():
            sys.exit(f"{path}: no result for {key}; run experiments/harness/run.py with "
                     "floor,b0,g_mom,g_mom_flip,g_mom_perm,g_mom_null,g_mom_held")
        result = json.loads(path.read_text())
        missing = [c for c in ("b0", "g_mom", *PLACEBOS, HELD_ARM) if c not in result["candidates"]]
        if missing or result["task"]["key"] != key:
            sys.exit(f"{path}: not {key} with b0, g_mom, the controls and the held arm (missing {missing})")
        out[key] = result
    return out


def reproduce_t11(results: Mapping[str, dict]) -> dict:
    """Per cut, seed and arm (b0, g_mom), the absolute difference of score between this run and try T11's committed
    result; passes if every one is within its tolerance (REPRO_T11_TOL points: the arms are the same code, so a larger
    difference means the run is not the comparison T11 made; g_mom carries the fit's PCA and k-means, whose floats can
    differ between runs, b0 does not)."""
    worst = {arm: 0.0 for arm in REPRO_T11_TOL}
    per_cut = {}
    for key, fname in CUTS.items():
        t11 = json.loads((T11_TASKS / fname).read_text())
        by_seed = {row["seed"]: row["scores"] for row in t11["rows"]}
        diffs = {}
        for row in results[key]["rows"]:
            s = row["seed"]
            if s not in by_seed:
                sys.exit(f"{key}: seed {s} is not in try T11's committed result")
            for arm in REPRO_T11_TOL:
                d = abs(float(row["scores"][arm]) - float(by_seed[s][arm]))
                diffs[arm] = max(diffs.get(arm, 0.0), d)
                worst[arm] = max(worst[arm], d)
        per_cut[key] = diffs
    return {"tolerance_points": dict(REPRO_T11_TOL), "largest_difference": worst, "per_cut": per_cut,
            "reproduced": all(worst[arm] <= REPRO_T11_TOL[arm] for arm in worst)}


def arm_panels(results: Mapping[str, dict]) -> dict:
    """Per cut and arm the mean over seeds of each metric of PANEL_METRICS (the harness run's rows' `metrics`, the
    panel T11's own results carry); what the spread reading reads."""
    out = {}
    for key, result in results.items():
        out[key] = {arm: {m: float(np.mean([row["metrics"][arm][m] for row in result["rows"]])) for m in PANEL_METRICS}
                    for arm in result["candidates"]}
    return out


def spread_reading(panels: Mapping[str, Mapping[str, Mapping[str, float]]]) -> dict:
    """The rows' readings of the spread hypothesis on variance_ratio, per cut. T20 (sign flip): if the gain is odd in
    the factor sign, the flip's variance_ratio falls near b0's minus what g_mom added (first order; a second-order term
    of the factors' size lifts both arms alike). T21's validity gate as a spread control: the permuted arm's
    variance_ratio within VR_TOL of g_mom's on both cuts; otherwise it says only that the assignment matters. T23's arm
    is reported beside them (the no-development arm at g_mom's factor size). Nothing here decides anything."""
    per_cut = {}
    for key, arms in panels.items():
        vr = {a: arms[a]["variance_ratio"] for a in arms}
        added = vr["g_mom"] - vr["b0"]
        expected = vr["b0"] - added
        per_cut[key] = {"b0": vr["b0"], "g_mom": vr["g_mom"], "added_by_g_mom": added, "floor": vr.get("floor"),
                        "flip": vr["g_mom_flip"], "flip_if_sign_odd": expected,
                        "flip_minus_sign_odd_expectation": vr["g_mom_flip"] - expected,
                        "flip_near_sign_odd_expectation": abs(vr["g_mom_flip"] - expected) <= VR_TOL,
                        "perm": vr["g_mom_perm"], "perm_minus_g_mom": vr["g_mom_perm"] - vr["g_mom"],
                        "perm_within_tolerance": abs(vr["g_mom_perm"] - vr["g_mom"]) <= VR_TOL,
                        "null": vr["g_mom_null"], "null_minus_g_mom": vr["g_mom_null"] - vr["g_mom"],
                        "held": vr[HELD_ARM], "held_minus_g_mom": vr[HELD_ARM] - vr["g_mom"]}
    return {"tolerance": VR_TOL, "per_cut": per_cut,
            "t21_is_a_spread_control": all(v["perm_within_tolerance"] for v in per_cut.values()),
            "t20_sign_odd_reading_holds": all(v["flip_near_sign_odd_expectation"] for v in per_cut.values())}


def held_reading(results: Mapping[str, dict], configured, config_sha256: str | None) -> dict:
    """Try T24: the adoption rule for the pseudobulk-held arm against b0 (`harness.decide` with b0 as the base: a paired
    mean gain above 0.3 points, ahead in at least 8 of 10 seeds on both counted cuts, no cut losing more than 0.3), with
    its validity gate (the arm's de_score and de_direction skills, mean over seeds, within HELD_DE_TOL of b0's on each
    cut; otherwise the pseudobulk was not held and the reading is void). passes: the gain survives the hold, so it is
    not a per-type mean effect."""
    d = harness.decide(list(results.values()), configured, "b0", HELD_ARM, config_sha256=config_sha256)
    gaps = {}
    for key, result in results.items():
        gaps[key] = {m: float(np.mean([row["skills"][HELD_ARM][m] - row["skills"]["b0"][m] for row in result["rows"]]))
                     for m in ("de_score", "de_direction")}
    within = all(abs(x) <= HELD_DE_TOL for g in gaps.values() for x in g.values())
    return {"try": HELD_TRY, "decision": d,
            "group_margins": {k: o1.group_margins(results[k], base="b0", change=HELD_ARM) for k in CUTS},
            "validity": {"de_skill_gap_vs_b0": gaps, "tolerance": HELD_DE_TOL, "within_tolerance": within},
            "void": not within, "passes": bool(d["keep"] and within)}


def read_placebos(results: Mapping[str, dict], configured, config_sha256: str | None) -> dict:
    """For each control the adoption rule with the control as base and g_mom as change (`harness.decide`: a paired mean
    gain above 0.3 points, ahead in at least 8 of 10 seeds, on both counted cuts, no cut losing more than 0.3), with the
    margin by scored group per cut; the held arm's reading (`held_reading`), the metric panels and the spread reading
    (`arm_panels`, `spread_reading`). The momentum reading survives only if every control is beaten by that rule and the
    held arm keeps its gain."""
    out = {}
    for placebo, label in PLACEBOS.items():
        d = harness.decide(list(results.values()), configured, placebo, "g_mom", config_sha256=config_sha256)
        out[placebo] = {"try": label, "decision": d,
                        "group_margins": {k: o1.group_margins(results[k], base=placebo, change="g_mom") for k in CUTS}}
    held = held_reading(results, configured, config_sha256)
    panels = arm_panels(results)
    return {"placebos": out, "held": held, "panels": panels, "spread": spread_reading(panels),
            "survives": all(p["decision"]["keep"] for p in out.values()) and held["passes"]}


def placebo_text(record: Mapping) -> str:
    st = record["stamp"]
    lines = [f"Tries T20, T21, T23 and T24: g_mom against its controls on both E10.5 cuts (rule: g_mom beats each of "
             f"{', '.join(PLACEBOS)} by a paired mean gain above {MIN_GAIN} points and ahead in at least "
             f"{MIN_AHEAD} of 10 seeds on both cuts, no cut losing more than {MAX_LOSS}; and {HELD_ARM} keeps that "
             f"gain over b0, its validity gate holding): momentum reading "
             f"{'SURVIVES' if record['survives'] else 'DOES NOT SURVIVE'}"]
    rep = record["reproduction_of_t11"]
    lines.append(f"Reproduction of try T11's committed scores (per seed; b0 within {rep['tolerance_points']['b0']}, "
                 f"g_mom within {rep['tolerance_points']['g_mom']} points): "
                 f"{'yes' if rep['reproduced'] else 'NO, the reading is void'} (largest difference b0 "
                 f"{rep['largest_difference']['b0']:.4f}, g_mom {rep['largest_difference']['g_mom']:.4f})")
    for name, p in record["placebos"].items():
        d = p["decision"]
        lines.append(f"{p['try']} {name}: {'beaten' if d['keep'] else 'NOT beaten'} ({d['reason']})")
        for key, v in d["tasks"].items():
            if not v.get("complete"):
                lines.append(f"    {key}: incomplete {v}")
                continue
            m = p["group_margins"][key]
            lines.append(f"    {key}: g_mom minus {name}: paired mean gain {v['mean_gain']:+.2f}, ahead in "
                         f"{v['seeds_ahead']} of 10, passes {v['passes']}; by group "
                         + ", ".join(f"{g} {x:+.2f}" for g, x in m["margins"].items()))
    held, dh = record["held"], record["held"]["decision"]
    verdict = ("VOID (de_score or de_direction drifted past the gate: the pseudobulk was not held)" if held["void"]
               else "keeps its gain" if held["passes"]
               else "loses the gain: a per-type mean effect, consistent with O1")
    lines.append(f"{held['try']} {HELD_ARM} (the pseudobulk held per gene, as T9) against b0: {verdict} "
                 f"({dh['reason']})")
    for key, v in dh["tasks"].items():
        if v.get("complete"):
            m = held["group_margins"][key]
            lines.append(f"    {key}: {HELD_ARM} minus b0: paired mean gain {v['mean_gain']:+.2f}, ahead in "
                         f"{v['seeds_ahead']} of 10, passes {v['passes']}; by group "
                         + ", ".join(f"{g} {x:+.2f}" for g, x in m["margins"].items()))
    for key, g in held["validity"]["de_skill_gap_vs_b0"].items():
        lines.append(f"    validity on {key}: de_score skill gap {g['de_score']:+.3f}, de_direction skill gap "
                     f"{g['de_direction']:+.3f} (tolerance {held['validity']['tolerance']})")
    spread = record["spread"]
    lines.append(f"Spread reading (variance_ratio, deciding nothing): T21 "
                 f"{'may' if spread['t21_is_a_spread_control'] else 'may NOT'} be read as a spread control (its "
                 f"variance_ratio within {spread['tolerance']} of g_mom's on both cuts: "
                 f"{spread['t21_is_a_spread_control']}; otherwise it says only that the assignment matters); T20's "
                 f"sign-odd reading (the flip near b0 minus what g_mom added) holds on both cuts: "
                 f"{spread['t20_sign_odd_reading_holds']}. A pass of T20 says the effect is odd in the factor sign, a "
                 "pass of T21 that the assignment matters; neither is a spread control or evidence of momentum, which "
                 "is what T23 and T24 are for")
    for key, v in spread["per_cut"].items():
        lines.append(f"    {key}: variance_ratio floor {v['floor']:.4f}, b0 {v['b0']:.4f}, g_mom {v['g_mom']:.4f} "
                     f"(added {v['added_by_g_mom']:+.4f}), flip {v['flip']:.4f} (near {v['flip_if_sign_odd']:.4f} if "
                     f"sign-odd), perm {v['perm']:.4f} (g_mom {v['perm_minus_g_mom']:+.4f}), null {v['null']:.4f} "
                     f"(g_mom {v['null_minus_g_mom']:+.4f}), held {v['held']:.4f}")
    for key, arms in record["panels"].items():
        lines.append(f"    {key} panel (mean over seeds): "
                     + "; ".join(f"{a}: " + ", ".join(f"{m} {arms[a][m]:.4f}" for m in SPREAD_METRICS) for a in arms))
    lines.append("Limit (one transition): E8.5 -> E9.5 momentum on E9.5 -> E10.5 is one fit and one permutation draw; "
                 "the ten seeds are draw and scoring noise, so a control beaten in 8 of 10 is one comparison, not a "
                 "replication of the biology")
    lines.append(f"Decided on {record['decided_on']}; run {record['run']}; commit {st['commit']}, tree clean "
                 f"{st['tree_clean']}, UTC {st['utc']}")
    return "\n".join(lines) + "\n"


def cmd_placebo(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("placebo.json", "placebo.txt"))
    started = tree_state()
    runs = Path(a.runs)
    results = load_run(runs)
    repro = reproduce_t11(results)
    configured = harness.load_tasks()["T1:test"]
    read = read_placebos(results, configured, harness._sha256(harness.CONFIG))
    if not repro["reproduced"]:
        read["survives"] = False
    inputs = [runs / f for f in CUTS.values()] + [T11_TASKS / f for f in CUTS.values()]
    record = {"label": "tries T20, T21, T23 and T24: g_mom against its controls on both E10.5 cuts", "run": str(runs),
              "rule": {"min_mean_gain": MIN_GAIN, "min_seeds_ahead": MIN_AHEAD, "max_cut_loss": MAX_LOSS,
                       "every_control_on_both_cuts": True, "held_arm_over_b0": True, "held_de_skill_tol": HELD_DE_TOL,
                       "t21_variance_ratio_tol": VR_TOL},
              "reproduction_of_t11": repro, **read,
              "decided_on": next(iter(read["placebos"].values()))["decision"]["decided_on"],
              "stamp": stamp(started, inputs)}
    write_outputs(out, "placebo", record, placebo_text(record))
    if not repro["reproduced"]:
        sys.exit("the run does not reproduce try T11's committed b0 and g_mom: the reading is void; it goes to the "
                 "captain")


# ---------------------------------------------------------------------------------------------------------------------
# Try T22: the composition reading (the oracle base of the t1-grid)


def oracle_shares(cut_categories: Sequence[str]) -> dict[str, float]:
    """payoff.py t1-grid's oracle base: the chamber categories of the chamber cut (its categories without its
    surrounding cells), renormalised to 1. Nothing here reads the target's expression, only its categories."""
    mix = t1_scope.composition(list(cut_categories))
    chamber = set(t1_scope.load_spec().in_scope("chamber"))
    kept = {k: v for k, v in mix.items() if k in chamber and v > 0}
    total = sum(kept.values())
    return {k: v / total for k, v in kept.items()}


def oracle_arms(ctx: harness.BuildContext, shares: Mapping[str, float], held: Sequence[str]) -> tuple[dict, object]:
    """b0, g_mom and g_mom_held at the mix `shares`, on one build context of the reference (E9.5 with start E8.5, the
    task's own inputs): B0's population (`ts.b0_population`: payoff.py's b0 at that mix, cell for cell), and on the same
    population g_mom's edit: the fold's finished factors (`harness._g_mom_factors`: fitted once per loaded fold on the
    reference's own chamber categories, E9.5's, as g_mom's) applied through `ts.apply_state_factors` on stream
    [seed, 5], R1's held labels unedited, and the same edit with each category's pseudobulk held (rescale=True, as
    `harness._g_mom_held`). Returns ({arm: Stage}, as the harness would finish them, and the B0 population at the
    mix)."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    start, start_src = harness._shallow(ctx.inputs["start"]), ctx.task.role_stems["start"]
    own = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    pop = ts.b0_population(ref, start, shares, n, ctx.seed, list(held), ref_source=src, start_source=start_src,
                           damping=harness.B0_DAMPING)
    genes = ctx.board.genes
    b0 = harness.finish(ctx, ref, dense(pop[:, genes].X).astype(np.float32), None)
    factors = harness._g_mom_factors(ctx, ref, src, start, start_src, own)
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    no_edit = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    arms = {"b0": b0}
    for arm, rescale in (("g_mom", False), (HELD_ARM, True)):
        moved, _ = ts.apply_state_factors(pop, factors, list(held), [ctx.seed, 5], rescale=rescale, no_edit=no_edit)
        arms[arm] = harness.finish(ctx, ref, dense(moved[:, genes].X).astype(np.float32), None)
    return arms, pop


def companion_mix(own: Mapping[str, float], cut_mix: Mapping[str, float]) -> tuple[dict[str, float], list[str]]:
    """The mix T22's companion truth is re-sampled to: the prediction's own chamber mix `own` (its categories the cut
    holds, renormalised) scaled to the chamber cells' share of the cut, and the cut's surrounding cells keeping their
    share (the t1-cm axis's convention, experiments/s2_hedge/payoff.py). Returns (mix, the categories of `own` the cut
    lacks); nothing here reads the target's expression, only its categories' shares."""
    chamber = set(t1_scope.load_spec().in_scope("chamber"))
    kept = {k: v for k, v in own.items() if cut_mix.get(k, 0.0) > 0}
    dropped = sorted(k for k in own if k not in kept)
    total = sum(kept.values())
    chamber_share = sum(v for k, v in cut_mix.items() if k in chamber)
    mix = {k: v / total * chamber_share for k, v in kept.items()}
    mix.update({k: v for k, v in cut_mix.items() if k not in chamber and v > 0})
    return mix, dropped


def floor_anchors(t, board, truth: Stage, ref_stage: Stage, seed: int):
    """The floor's anchors for one truth and scoring seed, with the harness's checks (the floor scores exactly 50)."""
    f = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, seed), truth, ref_stage, seed=seed, anchors="computed",
                    build=harness.BUILD, protocol=harness.PROTOCOL)
    missing = sorted(set(RANKED_METRICS[t.task]) - set(f.anchors))
    if abs(f.score - 50.0) > 1e-9 or f.extras["inverted_anchors"] or missing:
        raise ValueError(f"{t.key} seed {seed}: the floor scores {f.score} (inverted {f.extras['inverted_anchors']}, "
                         f"missing anchors {missing})")
    return f


def score_arms(board, truth: Stage, ref_stage: Stage, seed: int, anchors, arms: Mapping[str, Stage]) -> dict:
    """The arms against one truth's shared anchors: score, group skills and the metric panel (PANEL_METRICS) of each."""
    out = {}
    for name, pred in arms.items():
        r = score_board(board, pred, truth, ref_stage, seed=seed, anchors=anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        out[name] = {"score": float(r.score), "groups": {g: float(v) for g, v in r.group_skill.items()},
                     "metrics": {m: float(r.metrics[m]) for m in PANEL_METRICS if m in r.metrics}}
    return out


def t22_reading(rows: Sequence[Mapping], key: str = "chamber", change: str = "g_mom") -> dict:
    """Per-seed gain (`change` minus b0, board points) on one cut (or the companion), its mean, the seeds at or above
    the rule's bar (+0.3) and above 0, the rule (at least 8 of 10 at or above the bar, on all ten seeds), the margin by
    scored group (mean over seeds of 100 x weight x the difference of group skills) and the mean over seeds of the
    metric panel of each arm in the row (the spread reading)."""
    gains = {int(r["seed"]): float(r[key][change]["score"] - r[key]["b0"]["score"]) for r in rows}
    groups = [(g, w) for g, w, _ in GROUPS["T1"]]
    margins = {g: float(np.mean([100.0 * w * (r[key][change]["groups"][g] - r[key]["b0"]["groups"][g]) for r in rows]))
               for g, w in groups}
    v = np.array(list(gains.values()))
    at_or_above = int((v >= T22_MIN_GAIN).sum())
    panel = {arm: {m: float(np.mean([r[key][arm]["metrics"][m] for r in rows])) for m in entry.get("metrics", {})}
             for arm, entry in rows[0][key].items() if isinstance(entry, dict)}
    return {"gains_per_seed": gains, "mean_gain": float(v.mean()), "seeds_at_or_above_bar": at_or_above,
            "seeds_above_zero": int((v > 0).sum()), "n_seeds": len(v), "group_margins": margins, "panel": panel,
            "passes": len(v) == len(SEEDS) and at_or_above >= T22_MIN_SEEDS}


def composition_text(record: Mapping) -> str:
    rule = record["rule"]
    ch, wh = record["chamber"], record["whole"]
    st = record["stamp"]
    rep = record["reproduction_of_oracle_b0"]
    held, companion, lib = record["chamber_held"], record["companion"], record["median_library_size"]
    verdict = ("incomplete (not all ten seeds)" if not record["complete_seeds"] else
               ("VOID: the oracle-base b0 does not reproduce the committed t1-grid; the reading goes to the captain"
                if rep["reproduced"] is False else
                ("PASSES" if ch["passes"] else
                 "FAILS: the gain is not shown to survive the composition match (TRIES rule 2's third condition); "
                 "g_mom is not proposed. A fail is inconclusive on power: the matched mix shrinks the edited share "
                 "(the row states why), and the companion below does not")))
    repro = {True: "yes", False: "NO, the reading is void", None: "not read (not all ten seeds)"}[rep["reproduced"]]
    lines = [f"Try T22: {LABEL_T22}",
             f"Rule (chamber cut): the gain at or above {rule['min_gain']:+.1f} points in at least {rule['min_seeds']} "
             f"of {len(SEEDS)} seeds: {verdict}",
             f"Reproduction of the committed t1-grid's oracle-base b0 ({ORACLE_ROW}) on the chamber cut: fresh mean "
             f"{rep['fresh_mean']:.3f}, committed {rep['committed_mean']:.3f} (tolerance {rep['tolerance_points']}): "
             f"{repro}",
             f"Chamber cut, g_mom minus b0 at the oracle mix: mean {ch['mean_gain']:+.2f}; seeds at or above the bar "
             f"{ch['seeds_at_or_above_bar']} of {ch['n_seeds']}, above 0 {ch['seeds_above_zero']}; per seed "
             + ", ".join(f"{g:+.2f}" for g in ch["gains_per_seed"].values())
             + "; by group " + ", ".join(f"{g} {x:+.2f}" for g, x in ch["group_margins"].items()),
             "    variance_ratio at the oracle mix (mean over seeds): "
             + ", ".join(f"{a} {ch['panel'][a]['variance_ratio']:.4f}" for a in ch["panel"]
                         if "variance_ratio" in ch["panel"][a]),
             f"Chamber cut, {HELD_ARM} minus b0 at the oracle mix (reported, deciding nothing; the pseudobulk held, "
             f"O1's construction): mean {held['mean_gain']:+.2f}; seeds at or above the bar "
             f"{held['seeds_at_or_above_bar']} of {held['n_seeds']}; by group "
             + ", ".join(f"{g} {x:+.2f}" for g, x in held["group_margins"].items()),
             f"Companion (reported, deciding nothing; g_mom minus b0 at the prediction's own chamber mix against the "
             f"chamber truth re-sampled to that mix, seed {record['companion_truth']['resample_seed']}; the edited "
             f"share does not shrink): mean {companion['mean_gain']:+.2f}; seeds at or above the bar "
             f"{companion['seeds_at_or_above_bar']} of {companion['n_seeds']}, above 0 "
             f"{companion['seeds_above_zero']}; by group "
             + ", ".join(f"{g} {x:+.2f}" for g, x in companion["group_margins"].items()),
             f"Whole cut (reported, deciding nothing: a whole-mix oracle base could be drawn from E9.5's non-chamber "
             f"cells, but it would not be the chamber-scoped candidate, so the whole cut is scored for the same two "
             f"predictions and keeps its non-chamber composition gap): mean {wh['mean_gain']:+.2f}; seeds at or above "
             f"the bar {wh['seeds_at_or_above_bar']} of {wh['n_seeds']}; by group "
             + ", ".join(f"{g} {x:+.2f}" for g, x in wh["group_margins"].items()),
             "Oracle mix (shares of the prediction): "
             + ", ".join(f"{k} {v:.3f}" for k, v in sorted(record["oracle_shares"].items())),
             f"Median library size (docs/DATA_POLICY.md 2026-11-04 (f)(3); first seed): the oracle-mix population "
             f"{lib['oracle_mix_population']:.0f}, E10.5's {lib['e10_5_released']:.0f}, ratio {lib['ratio']:.3f}"
             + (" (beyond 5%: this goes to the captain)" if lib["beyond_5_percent"] else " (within 5%)"),
             "Limit: one transition and one fit; the seeds are draw and scoring noise, not biological replication.",
             f"Commit {st['commit']}, tree clean {st['tree_clean']}, UTC {st['utc']}"]
    return "\n".join(lines) + "\n"


def median_library_size(population, released: Path) -> dict:
    """docs/DATA_POLICY.md 2026-11-04 (f)(3)'s report: the drawn population's median library size (total counts per
    cell of the population's own counts) against that of the released file (E10.5's), and whether it shifted by more
    than 5% (sends the rule to the captain)."""
    def med(adata) -> float:
        return float(np.median(np.asarray(counts_of(adata).sum(1)).ravel()))

    pop_median = med(population)
    e105_median = med(load_stage(released, lineage=False))
    ratio = pop_median / e105_median if e105_median > 0 else float("nan")
    return {"oracle_mix_population": pop_median, "e10_5_released": e105_median, "ratio": ratio,
            "beyond_5_percent": bool(abs(ratio - 1.0) > LIBRARY_SHIFT_TOL)}


def cmd_composition(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("composition.json", "composition.txt"))
    started = tree_state()
    seeds = parse_seeds(a.seeds)
    data = Path(a.data)
    tasks = {t.key: t for t in harness.load_tasks()["T1:test"]}
    t_ch, t_wh = tasks[TASK_CHAMBER], tasks[TASK_WHOLE]
    boards = load_boards()
    genes = boards[t_ch.panel_board].genes
    cut, cut_cat = payoff.chamber_cut(t_ch, data, genes)
    shares = oracle_shares(cut_cat)
    grid = json.loads(T1_GRID.read_text())
    want_mix = grid["pred_mix"][ORACLE_ROW]
    if set(want_mix) != set(shares) or max(abs(want_mix[k] - shares[k]) for k in shares) > 1e-9:
        sys.exit("the oracle mix is not the committed t1-grid's: the construction has changed; the reading is void")
    truth_ch = Stage.from_anndata(cut, genes)
    truth_wh = Stage.from_anndata(read_stage(data / t_wh.target, genes), genes)
    ref_stage = Stage.from_anndata(read_stage(data / t_ch.reference, genes), genes)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t_ch.inputs.items()}
    held = list(t1_scope.held_flat_genes())
    board_ch = harness.task_board(t_ch, boards, truth_ch.n, ref_stage.n)
    board_wh = harness.task_board(t_wh, boards, truth_wh.n, ref_stage.n)
    views = {"chamber": (t_ch, board_ch, truth_ch), "whole": (t_wh, board_wh, truth_wh)}
    ref_role = harness.TaskView.of(t_ch).reference_role
    own = harness.t1_chamber_shares(harness._shallow(inputs[ref_role]), harness.TaskView.of(t_ch).role_stems[ref_role])
    mix, dropped = companion_mix(own, t1_scope.composition(list(cut_cat)))
    resampled = payoff.stratified(cut_cat, mix, np.random.default_rng([TRUTH_RESAMPLE_SEED, 0]))
    truth_companion = Stage.from_anndata(cut[resampled], genes)
    board_companion = harness.task_board(t_ch, boards, truth_companion.n, ref_stage.n)
    rows, library = [], None
    for s in seeds:
        t0 = time.time()
        row: dict = {"seed": s}
        for key, (t, board, truth) in views.items():
            view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
            ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                       derived={"t1_held_flat": tuple(held)})
            arms, population = oracle_arms(ctx, shares, held)
            f = floor_anchors(t, board, truth, ref_stage, s)
            row[key] = {"floor": float(f.score), **score_arms(board, truth, ref_stage, s, f.anchors, arms)}
            if library is None:  # DATA_POLICY 2026-11-04 (f)(3); the first seed's draw
                library = median_library_size(population, data / t_wh.target)
            if key == "chamber":  # the companion: same seed and context, the prediction's own mix, re-sampled truth
                own_arms = {"b0": harness._b0_t1(harness.B0_DAMPING)(ctx), "g_mom": harness._g_mom(ctx)}
                fc = floor_anchors(t, board_companion, truth_companion, ref_stage, s)
                row["companion"] = {"floor": float(fc.score), **score_arms(board_companion, truth_companion, ref_stage,
                                                                          s, fc.anchors, own_arms)}
        row["seconds"] = round(time.time() - t0, 1)
        rows.append(row)
        print(f"T22 seed {s} ({row['seconds']} s): chamber gain "
              f"{row['chamber']['g_mom']['score'] - row['chamber']['b0']['score']:+.2f}, whole gain "
              f"{row['whole']['g_mom']['score'] - row['whole']['b0']['score']:+.2f}, companion gain "
              f"{row['companion']['g_mom']['score'] - row['companion']['b0']['score']:+.2f}", flush=True)
    chamber = t22_reading(rows, "chamber")
    whole = t22_reading(rows, "whole")
    chamber_held = t22_reading(rows, "chamber", change=HELD_ARM)
    companion = t22_reading(rows, "companion")
    fresh = float(np.mean([r["chamber"]["b0"]["score"] for r in rows]))
    committed = float(grid["mean"][ORACLE_ROW][ORACLE_ROW])
    full = tuple(sorted(seeds)) == tuple(SEEDS)
    reproduced = None if not full else abs(fresh - committed) <= REPRO_SCORE_TOL
    inputs_paths = [data / f for f in sorted({t_ch.target, t_ch.reference, *t_ch.inputs.values()})] + [
        t1_scope.HELD_FLAT_GENES, harness.CONFIG, T1_GRID]
    record = {"label": LABEL_T22, "rule": {"min_gain": T22_MIN_GAIN, "min_seeds": T22_MIN_SEEDS, "cut": TASK_CHAMBER},
              "oracle_shares": shares, "chamber": chamber, "whole": whole, "chamber_held": chamber_held,
              "companion": companion,
              "companion_truth": {"mix": mix, "categories_of_the_prediction_the_cut_lacks": dropped,
                                  "resample_seed": TRUTH_RESAMPLE_SEED},
              "median_library_size": library, "complete_seeds": full,
              "reproduction_of_oracle_b0": {"fresh_mean": fresh, "committed_mean": committed,
                                            "tolerance_points": REPRO_SCORE_TOL, "reproduced": reproduced},
              "seeds": list(seeds), "rows": rows, "stamp": stamp(started, inputs_paths)}
    write_outputs(out, "composition", record, composition_text(record))
    if reproduced is False:
        sys.exit("the oracle-base b0 does not reproduce the committed t1-grid: the reading is void; it goes to the "
                 "captain")


# ---------------------------------------------------------------------------------------------------------------------
# Diagnostic N21: the factors at the test pair


def free_genes(genes: Sequence[str], held: set[str]) -> np.ndarray:
    return np.array([g not in held for g in genes], dtype=bool)


def factor_size(cf: ts.CategoryFactors, free: np.ndarray) -> dict:
    """The size of a finished factor matrix on the genes not held flat: the mean of |log f| over metacells (cell-
    weighted) and genes, the mean over the entries not zeroed and its 90% quantile, the share of entries not zeroed,
    the share of genes touched in at least one metacell, and the share of the entries at the clip."""
    F = np.abs(cf.logf[:, free])
    w = cf.sizes / cf.sizes.sum()
    nz = F > 0
    return {"mean_abs_log_f": float((w * F.mean(1)).sum()),
            "mean_abs_log_f_nonzero": float(F[nz].mean()) if nz.any() else 0.0,
            "q90_abs_log_f_nonzero": float(np.quantile(F[nz], 0.9)) if nz.any() else 0.0,
            "share_entries_touched": float(nz.mean()), "share_genes_touched": float(nz.any(0).mean()),
            "share_entries_at_clip": float((F >= ts.CLIP_AT - 1e-12).mean()), "n_metacells": int(cf.n_metacells)}


def _columns(sf_genes: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """The positions of the genes `keep` (in that order) in the fit's gene order."""
    index = {g: i for i, g in enumerate(sf_genes)}
    return np.array([index[g] for g in keep], dtype=np.int64)


def leading_axis(cf: ts.CategoryFactors, genes: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """The per-category factor direction: the first right singular vector (gene loadings, unit norm, on the genes `keep`
    in that order) of the metacell factor matrix, each metacell weighted by the square root of its size; zeros when the
    matrix is all zero (no direction)."""
    M = np.sqrt(cf.sizes.astype(np.float64))[:, None] * cf.logf[:, _columns(genes, keep)]
    if not M.any():
        return np.zeros(len(keep))
    return np.linalg.svd(M, full_matrices=False)[2][0]


def axis_cosine(u: np.ndarray, v: np.ndarray) -> float:
    """|cosine| of two gene-space axes on the same genes (the sign of a singular vector is arbitrary); NaN when either
    is zero."""
    a, b = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    return abs(float(np.dot(u, v)) / (a * b)) if a > 0 and b > 0 else float("nan")


def gene_size(cf: ts.CategoryFactors, genes: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """The per-gene size of the factors: the cell-weighted standard deviation over metacells of each kept gene's log
    factor (the genes `keep`, in that order)."""
    F = cf.logf[:, _columns(genes, keep)]
    w = cf.sizes / cf.sizes.sum()
    mean = (w[:, None] * F).sum(0)
    return np.sqrt((w[:, None] * (F - mean) ** 2).sum(0))


def profile_spearman(x: np.ndarray, y: np.ndarray) -> float:
    """The Spearman correlation of two per-gene size profiles (NaN when either is constant)."""
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(spearmanr(x, y)[0])


def concordance(a: ts.StateFactors, b: ts.StateFactors, category: str, held: set[str]) -> dict:
    """The concordance of two fits on one category: on the genes both fits hold and the held-flat list leaves free (the
    same genes in the same order for both), the |cosine| of their leading axes (their absolute value: an axis has no
    sign) and the Spearman correlation of their per-gene sizes over the genes with a non-zero size in either fit (the
    many genes with no factor at all would otherwise tie and drive it)."""
    keep = sorted((set(a.genes) & set(b.genes)) - held)
    ca, cb = a.categories[category], b.categories[category]
    if not keep:
        return {"axis_abs_cosine": float("nan"), "size_spearman": float("nan")}
    sa, sb = gene_size(ca, a.genes, keep), gene_size(cb, b.genes, keep)
    moved = (sa > 0) | (sb > 0)
    spearman = profile_spearman(sa[moved], sb[moved]) if moved.sum() >= 3 else float("nan")
    return {"axis_abs_cosine": axis_cosine(leading_axis(ca, a.genes, keep), leading_axis(cb, b.genes, keep)),
            "size_spearman": spearman}


halves = ts.category_halves  # the halves of N21's split-half fits and no-development reference (the module's rule)


def classified(a, source: str, spec):
    c = harness._shallow(a)
    t1_scope.classify(c, source=source, spec=spec)
    return c


def split_half_fits(ref, ref_src: str, start, start_src: str, categories: Sequence[str], pair: int, spec):
    """The same fit on two random halves of each stage's cells of each category: halves drawn from
    default_rng([SPLIT_SEED, pair]), start stage first, then the reference, categories in sorted order. Returns (fit A,
    fit B), each `harness.g_mom_pair_factors` on the cells of its half, for the given categories."""
    rng = np.random.default_rng([SPLIT_SEED, pair])
    ms, mr = halves(classified(start, start_src, spec), rng), halves(classified(ref, ref_src, spec), rng)
    fit_a = harness.g_mom_pair_factors(ref[mr].copy(), ref_src, start[ms].copy(), start_src, categories=categories)[0]
    fit_b = harness.g_mom_pair_factors(ref[~mr].copy(), ref_src, start[~ms].copy(), start_src,
                                       categories=categories)[0]
    return fit_a, fit_b


def expected_depth_ratio(C, factors: np.ndarray, row: np.ndarray, held: np.ndarray, scale: np.ndarray) -> np.ndarray:
    """Per cell of counts C, the expected total count after `vec.core.apply_factors` (keep_depth) with these factors,
    over the cell's own total: c * f * s on the free genes, c on the held genes, with s = `scale`, the keep-depth
    constant `apply_factors` itself returned for this call (a change to that function's constant shows up here)."""
    C = sparse.csr_matrix(C, dtype=np.float64)
    n = C.shape[0]
    row_of = np.repeat(np.arange(n), np.diff(C.indptr))
    f = factors[row[row_of], C.indices]
    free = ~held[C.indices]
    expected = np.bincount(row_of, weights=C.data * np.where(free, f * scale[row_of], 1.0), minlength=n)
    depth = np.asarray(C.sum(1), dtype=np.float64).ravel()
    return np.divide(expected, depth, out=np.ones(n), where=depth > 0)


def n21_pool(ref, src: str, categories: Sequence[str], spec):
    """E10.5's chamber cells of the fitted categories as a population `apply_state_factors` can read (each cell its own
    source cell, its category and its R1 mask); no draw and no cell number leaves this function."""
    pop = classified(ref, src, spec)
    cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    pop = pop[np.isin(cats, list(categories))].copy()
    pop.obs["source_obs"] = np.asarray(pop.obs_names).astype(str)
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    held_labels = ts.r1_unbracketed_labels(np.unique(labels))
    return pop, labels, np.isin(labels, sorted(held_labels)), sorted(held_labels)


def n21_category(pop, labels, no_edit, factors: ts.StateFactors, held_genes: set[str],
                 category: str) -> tuple[dict, np.ndarray]:
    """One category of the test pair's fit on E10.5's chamber cells: the factor size, the shares of the category's cells
    touched and held by R1 (shares only) and the labels held; returns that record and the category's edited cells
    (boolean over the population, for the depth check)."""
    cf = factors.categories[category]
    cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    in_cat = cats == category
    src = pop.obs["source_obs"].astype(str).to_numpy()
    has_factor = np.fromiter((s in cf.index for s in src), dtype=bool, count=len(src))
    edited = in_cat & has_factor & ~no_edit
    out = {"factor_size": factor_size(cf, free_genes(factors.genes, held_genes)),
           "share_of_category_cells_touched": float(edited.sum() / max(in_cat.sum(), 1)),
           "share_of_category_cells_held_by_r1": float((in_cat & no_edit).sum() / max(in_cat.sum(), 1)),
           "held_labels_in_category": sorted(set(labels[in_cat & no_edit]))}
    return out, edited


def depth_checks(pop, no_edit, factors: ts.StateFactors, held_genes: set[str], edited_any: np.ndarray) -> dict:
    """Over the edited cells: the expected depth after the edit against the cell's own depth (`expected_depth_ratio`,
    with the keep-depth constant `vec.core.apply_factors` returns for these very factors: a move larger than DEPTH_TOL,
    which would mean the constant no longer keeps depth, stops the run), and the realised depth ratio after one
    application of `ts.apply_state_factors` on N21_APPLY_STREAM (median, 1% and 99% quantiles over the edited cells);
    whether every held cell's depth is identical; the share of the fitted categories' cells edited at all; and the
    expected and realised depth per category (per_category, over each category's edited cells)."""
    C = sparse.csr_matrix(pop.layers["counts"], dtype=np.int64)
    row, logf, _ = ts._factor_rows(pop, factors, t1_scope.CATEGORY_KEY, "source_obs")
    row = row.copy()
    row[no_edit] = logf.shape[0] - 1
    held = np.array([g in held_genes for g in factors.genes], dtype=bool)
    F = np.exp(logf)
    F[:, held] = 1.0
    _, scale, _ = apply_factors(C, F, row, held, keep_depth=True, seed=list(N21_APPLY_STREAM))
    ratio = expected_depth_ratio(C, F, row, held, scale)
    moved, _ = ts.apply_state_factors(pop, factors, sorted(held_genes), list(N21_APPLY_STREAM), rescale=False,
                                      no_edit=no_edit)
    C_new = sparse.csr_matrix(moved.layers["counts"], dtype=np.int64)
    depth = np.asarray(C.sum(1), dtype=np.float64).ravel()
    new_depth = np.asarray(C_new.sum(1), dtype=np.float64).ravel()
    realised = np.divide(new_depth, depth, out=np.ones_like(depth), where=depth > 0)
    idx = np.flatnonzero(edited_any)
    worst = float(np.abs(ratio[idx] - 1.0).max()) if len(idx) else 0.0
    if worst > DEPTH_TOL:
        sys.exit(f"an edited cell's expected depth moves by {worst:.2e} (allowed {DEPTH_TOL:g}): the edit does not "
                 "keep depth; do not use these factors")

    def spread(rows: np.ndarray) -> dict:
        return {"median": float(np.median(realised[rows])) if len(rows) else None,
                "q01": float(np.quantile(realised[rows], 0.01)) if len(rows) else None,
                "q99": float(np.quantile(realised[rows], 0.99)) if len(rows) else None}

    cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    per_category = {}
    for c in sorted(factors.categories):
        rows = np.flatnonzero(edited_any & (cats == c))
        largest = float(np.abs(ratio[rows] - 1.0).max()) if len(rows) else 0.0
        per_category[c] = {"max_abs_expected_ratio_minus_1": largest, "realised_depth_ratio": spread(rows)}
    return {"max_abs_expected_ratio_minus_1": worst,
            "share_of_fitted_cells_edited": float(edited_any.mean()) if pop.n_obs else 0.0,
            "realised_depth_ratio_edited": spread(idx),
            "realised_depth_identical_when_held": bool(np.array_equal(new_depth[no_edit], depth[no_edit])),
            "per_category": per_category}


def own_deviation(stage, cf: ts.CategoryFactors, genes: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """Per metacell (cf's, on the cells cf.index names in `stage`) and kept gene: the log of the metacell's own mean
    CP10k (the fit's pseudocount, the fit's gene universe for the cell totals) minus its cell-weighted mean over the
    category's metacells: the deviation that the factor's numerator (the metacell's own mean) contains, so that a factor
    correlated with it expands the state contrast whether or not anything developed."""
    names = list(cf.index)
    rows = stage.obs_names.get_indexer(names)
    if (rows < 0).any():
        raise KeyError("a metacell cell is missing from the stage")
    labels = np.fromiter((cf.index[n] for n in names), dtype=np.int64, count=len(names))
    C = counts_of(stage, list(genes))[rows]
    means = np.vstack([mean_cp10k(C[labels == m]) for m in range(cf.n_metacells)])
    logm = np.log(means + ts.PSEUDOCOUNT)[:, _columns(genes, keep)]
    w = cf.sizes / cf.sizes.sum()
    return logm - (w[:, None] * logm).sum(0)


def amplification(F: np.ndarray, dev: np.ndarray, sizes: np.ndarray) -> dict:
    """The regression of finished factors F on the metacell's own deviation `dev` over the metacells and the kept genes
    (through the origin, metacells weighted by cells, genes equally): the slope and the share of F's sum of squares it
    explains (R^2, the amplification share: how much of the factor is the metacell's own deviation expanded)."""
    w = (sizes / sizes.sum())[:, None]
    num, den = float((w * F * dev).sum()), float((w * dev * dev).sum())
    slope = num / den if den > 0 else float("nan")
    ss_tot = float((w * F * F).sum())
    r2 = 1.0 - float((w * (F - slope * dev) ** 2).sum()) / ss_tot if ss_tot > 0 and den > 0 else float("nan")
    return {"slope": slope, "r2": r2}


def projected_fit(sf: ts.StateFactors, stage, categories: Sequence[str],
                  held_genes: set[str]) -> tuple[ts.StateFactors, dict]:
    """A copy of `sf` with each named category's factors reduced on the free genes by that category's amplification
    component (slope x own deviation, `amplification`), and the per-category amplification (the slope and R^2).
    Held-flat genes keep their factors (they are not free)."""
    keep = [g for g in sf.genes if g not in held_genes]
    cats, amp = dict(sf.categories), {}
    for c in categories:
        cf = sf.categories[c]
        dev = own_deviation(stage, cf, sf.genes, keep)
        cols = _columns(sf.genes, keep)
        amp[c] = amplification(cf.logf[:, cols], dev, cf.sizes)
        logf = cf.logf.copy()
        if np.isfinite(amp[c]["slope"]):
            logf[:, cols] -= amp[c]["slope"] * dev
        cats[c] = replace(cf, logf=logf)
    return replace(sf, categories=cats), amp


def union_share(cf: ts.CategoryFactors, genes: Sequence[str], union: set[str], held_genes: set[str]) -> dict:
    """The share of a finished factor matrix's size carried by R2's union genes (`hygiene.load_union`, the technical
    modules R2 holds flat in Feng's direction), among the genes not held flat: the cell-weighted sum of squares (l2) and
    of absolute values (l1) over the metacells, and how many union genes the fit holds."""
    free = np.array([g not in held_genes for g in genes], dtype=bool)
    on = free & np.array([g in union for g in genes], dtype=bool)
    w = (cf.sizes / cf.sizes.sum())[:, None]
    l2, l2_on = float((w * cf.logf[:, free] ** 2).sum()), float((w * cf.logf[:, on] ** 2).sum())
    l1, l1_on = float((w * np.abs(cf.logf[:, free])).sum()), float((w * np.abs(cf.logf[:, on])).sum())
    return {"l2_share": l2_on / l2 if l2 > 0 else float("nan"), "l1_share": l1_on / l1 if l1 > 0 else float("nan"),
            "union_genes_in_fit": int(on.sum())}


def size_ratio(no_dev: ts.StateFactors, fit: ts.StateFactors, category: str, held_genes: set[str]) -> float | None:
    """How large the no-development fit's finished factors are against g_mom's, on one category (the ratio of their
    cell-weighted root mean squares over the metacells and the free genes): the size the no-development control (T23)
    has before it is scaled up to g_mom's, on these stages. None when either fit lacks the category or its size is 0."""
    if category not in no_dev.categories or category not in fit.categories:
        return None
    free_no = np.array([g not in held_genes for g in no_dev.genes], dtype=bool)
    free_fit = np.array([g not in held_genes for g in fit.genes], dtype=bool)
    a, b = no_dev.categories[category], fit.categories[category]
    denominator = ts._factor_rms(b.logf, b.sizes, free_fit)
    return ts._factor_rms(a.logf, a.sizes, free_no) / denominator if denominator > 0 else None


def no_development_fit(stage, src: str, categories: Sequence[str], pair: str, spec) -> ts.StateFactors:
    """g_mom's finished fit between two random halves of one stage (the halves from default_rng([SPLIT_SEED,
    NO_DEV_STREAMS[pair]]), metacells on one half, partners in the other): nothing develops between them, so the
    concordance of two such fits is what the concordance statistics read without development."""
    mask = halves(classified(stage, src, spec), np.random.default_rng([SPLIT_SEED, NO_DEV_STREAMS[pair]]))
    return harness.g_mom_pair_factors(stage[~mask].copy(), src, stage[mask].copy(), src, categories=categories)[0]


def n21_text(res: Mapping) -> str:
    def num(x, f="{:.3f}"):
        return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f.format(x)

    d = res["depth"]
    lines = [f"N21: {LABEL_N21}. Decides nothing; reported to the captain with the ruling.",
             f"R1's held labels (no edit): {', '.join(res['r1_held_labels']) or 'none'}. Cells edited: "
             f"{num(d['share_of_fitted_cells_edited'], '{:.2f}')} of the fitted categories' cells",
             f"Depth, all edited cells: largest |E[depth after] / depth - 1| {d['max_abs_expected_ratio_minus_1']:.2e} "
             f"(kept by construction); realised depth ratio median {num(d['realised_depth_ratio_edited']['median'])}, "
             f"1% to 99% {num(d['realised_depth_ratio_edited']['q01'])} to "
             f"{num(d['realised_depth_ratio_edited']['q99'])}; held cells' depth identical "
             f"{d['realised_depth_identical_when_held']}"]
    for c, r in res["test_pair"].items():
        fs = r["factor_size"]
        lines.append(f"{c}: {fs['n_metacells']} metacells; mean |log f| {num(fs['mean_abs_log_f'])} (not zeroed "
                     f"{num(fs['mean_abs_log_f_nonzero'])}, 90% quantile {num(fs['q90_abs_log_f_nonzero'])}); entries "
                     f"touched {num(fs['share_entries_touched'])}, genes touched {num(fs['share_genes_touched'])}, "
                     f"entries at the clip {num(fs['share_entries_at_clip'])}; cells touched "
                     f"{num(r['share_of_category_cells_touched'], '{:.2f}')} of the category (held by R1 "
                     f"{num(r['share_of_category_cells_held_by_r1'], '{:.2f}')}: "
                     f"{', '.join(r['held_labels_in_category']) or 'none'}); expected depth moves at most "
                     f"{d['per_category'][c]['max_abs_expected_ratio_minus_1']:.2e}, realised depth ratio median "
                     f"{num(d['per_category'][c]['realised_depth_ratio']['median'])}")
        amp = r["amplification"]
        lines.append(f"    amplification (the finished factors' regression on the metacell's own deviation, "
                     f"through the origin): test pair slope {num(amp['test']['slope'])}, R^2 {num(amp['test']['r2'])}"
                     + ("" if "backtest" not in amp else
                        f"; backtest pair slope {num(amp['backtest']['slope'])}, R^2 {num(amp['backtest']['r2'])}"))
        un = r["r2_union_share"]
        lines.append("    share of the factor size on R2's union genes (the technical modules Feng's direction holds "
                     "flat; g_mom's step does not hold them): "
                     + ("not available" if un is None else
                        f"l2 {num(un['test']['l2_share'])}, l1 {num(un['test']['l1_share'])} of "
                        f"{un['test']['union_genes_in_fit']} union genes (test pair)"))
        con = r.get("concordance_with_backtest")
        if con is None:
            lines.append("    not in the backtest's fit: no concordance")
            continue
        lines.append(f"    concordance with the backtest's fit: |cos| of the leading axes "
                     f"{num(con['axis_abs_cosine'])}, Spearman of the per-gene sizes {num(con['size_spearman'])}")
        proj = r["concordance_after_projecting_amplification"]
        lines.append(f"    the same after projecting the amplification component out of each fit: |cos| "
                     f"{num(proj['axis_abs_cosine'])}, Spearman {num(proj['size_spearman'])}")
        nod = r["no_development_reference"]
        lines.append("    no-development reference (the concordance of two same-stage halves' fits, E10.5's and "
                     "E9.5's: what these statistics read with nothing developing): "
                     + ("n/a" if nod is None else
                        f"|cos| {num(nod['axis_abs_cosine'])}, Spearman {num(nod['size_spearman'])}"))
        ratio = r["no_development_size_ratio"]
        lines.append("    size of the no-development fit against g_mom's fit (before T23 scales it up to g_mom's "
                     f"size): test pair {num(ratio['test'])}, backtest pair {num(ratio['backtest'])}")
        ceil = r["split_half_ceiling"]
        lines.append(f"    split-half ceilings (half-data estimates, each halving the cells of a fit and so "
                     f"under-estimating the full-data ceiling; test pair; backtest pair; their geometric mean): axis "
                     f"|cos| {num(ceil['axis_abs_cosine_test'])}; {num(ceil['axis_abs_cosine_backtest'])}; "
                     f"{num(ceil['axis_geometric_mean'])}; size Spearman {num(ceil['size_spearman_test'])}; "
                     f"{num(ceil['size_spearman_backtest'])}")
    st = res["stamp"]
    lines.append(f"Commit {st['commit']}, tree clean {st['tree_clean']}, UTC {st['utc']}")
    return "\n".join(lines) + "\n"


def run_n21(data: Path, started: Mapping) -> dict:
    spec = t1_scope.load_spec()
    files = {"E8.5": "T1/E8.5_RNA.h5ad", "E9.5": "T1/E9.5_RNA.h5ad", "E10.5": "T1/E10.5_RNA.h5ad"}
    stage = {k: load_stage(data / f, lineage=False) for k, f in files.items()}
    held_genes = set(t1_scope.held_flat_genes())
    test_cats = [c for c, s in harness.t1_chamber_shares(stage["E10.5"], "E10.5_RNA").items() if s > 0]
    back_cats = [c for c, s in harness.t1_chamber_shares(stage["E9.5"], "E9.5_RNA").items() if s > 0]
    test_f = harness.g_mom_pair_factors(stage["E10.5"], "E10.5_RNA", stage["E9.5"], "E9.5_RNA", categories=test_cats)[0]
    back_f = harness.g_mom_pair_factors(stage["E9.5"], "E9.5_RNA", stage["E8.5"], "E8.5_RNA", categories=back_cats)[0]
    both = sorted(set(test_f.categories) & set(back_f.categories))
    test_a, test_b = split_half_fits(stage["E10.5"], "E10.5_RNA", stage["E9.5"], "E9.5_RNA", test_cats, 0, spec)
    back_a, back_b = split_half_fits(stage["E9.5"], "E9.5_RNA", stage["E8.5"], "E8.5_RNA", back_cats, 1, spec)
    no_dev_test = no_development_fit(stage["E10.5"], "E10.5_RNA", test_cats, "test", spec)
    no_dev_back = no_development_fit(stage["E9.5"], "E9.5_RNA", back_cats, "backtest", spec)
    try:
        union = hygiene.load_union()
        union_genes: set[str] | None = set(union["genes"])
    except hygiene.Refused as e:
        print(f"R2's union file is not available ({e}): the union-gene shares are not reported", flush=True)
        union, union_genes = None, None
    pop, labels, no_edit, held_labels = n21_pool(stage["E10.5"], "E10.5_RNA", list(test_f.categories), spec)
    test_rows, edited_any = {}, np.zeros(pop.n_obs, dtype=bool)
    test_proj, amp_test = projected_fit(test_f, stage["E10.5"], sorted(test_f.categories), held_genes)
    back_proj, amp_back = projected_fit(back_f, stage["E9.5"], both, held_genes)
    for c in sorted(test_f.categories):
        row, edited = n21_category(pop, labels, no_edit, test_f, held_genes, c)
        edited_any |= edited
        row["amplification"] = {"test": amp_test[c], **({"backtest": amp_back[c]} if c in both else {})}
        row["r2_union_share"] = None if union_genes is None else {
            "test": union_share(test_f.categories[c], test_f.genes, union_genes, held_genes)}
        if c in both:
            row["concordance_with_backtest"] = concordance(test_f, back_f, c, held_genes)
            row["concordance_after_projecting_amplification"] = concordance(test_proj, back_proj, c, held_genes)
            ca = concordance(test_a, test_b, c, held_genes) if {c} <= set(test_a.categories) & set(test_b.categories) \
                else None
            cb = concordance(back_a, back_b, c, held_genes) if {c} <= set(back_a.categories) & set(back_b.categories) \
                else None
            at = None if ca is None else ca["axis_abs_cosine"]
            ab = None if cb is None else cb["axis_abs_cosine"]
            both_axes = at is not None and ab is not None and np.isfinite(at) and np.isfinite(ab)
            row["split_half_ceiling"] = {
                "axis_abs_cosine_test": at, "axis_abs_cosine_backtest": ab,
                "axis_geometric_mean": float(np.sqrt(at * ab)) if both_axes else None,
                "size_spearman_test": None if ca is None else ca["size_spearman"],
                "size_spearman_backtest": None if cb is None else cb["size_spearman"]}
            row["no_development_reference"] = (concordance(no_dev_test, no_dev_back, c, held_genes)
                                               if {c} <= set(no_dev_test.categories) & set(no_dev_back.categories)
                                               else None)
            row["no_development_size_ratio"] = {
                "test": size_ratio(no_dev_test, test_f, c, held_genes),
                "backtest": size_ratio(no_dev_back, back_f, c, held_genes)}
        test_rows[c] = row
    res = {"label": LABEL_N21, "pairs": {"test": ["E9.5", "E10.5"], "backtest": ["E8.5", "E9.5"]},
           "categories_test_pair": sorted(test_f.categories), "categories_backtest_pair": sorted(back_f.categories),
           "categories_in_both": both, "skipped_test_pair": sorted(test_f.skipped),
           "skipped_backtest_pair": sorted(back_f.skipped), "r1_held_labels": held_labels,
           "test_pair": test_rows, "depth": depth_checks(pop, no_edit, test_f, held_genes, edited_any),
           "settings": {"damping": harness.B0_DAMPING, "split_seed": SPLIT_SEED,
                        "no_development_streams": dict(NO_DEV_STREAMS), "apply_stream": list(N21_APPLY_STREAM),
                        "zero_below": ts.ZERO_BELOW, "clip_at": ts.CLIP_AT}}
    inputs = [data / f for f in files.values()] + [t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ts.R1_CONFIG]
    if union is not None:
        inputs.append(REPO_ROOT / union["file"])
    res["stamp"] = stamp(started, inputs)
    return res


def cmd_n21(a) -> None:
    out = RESULTS / a.out
    refuse_existing(out, ("n21.json", "n21.txt"))
    started = tree_state()
    res = run_n21(Path(a.data), started)
    write_outputs(out, "n21", res, n21_text(res))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("placebo", help="tries T20 and T21, from a harness run's results")
    p.add_argument("--runs", required=True, help="the harness run folder (experiments/harness/run.py run --out)")
    p.add_argument("--out", required=True, help="folder name under experiments/t1_scope/results/")
    p = sub.add_parser("composition", help="try T22: g_mom against b0 at the truth's own chamber mix")
    p.add_argument("--seeds", default="0-9")
    p.add_argument("--out", required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("n21", help="diagnostic N21: the factors at the test pair (target-free)")
    p.add_argument("--out", required=True)
    p.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    {"placebo": cmd_placebo, "composition": cmd_composition, "n21": cmd_n21}[a.cmd](a)


if __name__ == "__main__":
    main()
