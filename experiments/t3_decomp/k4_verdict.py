#!/usr/bin/env python
"""Tries K4 and K4c (experiments/harness/TRIES.md): the keep-or-stop rule of K4's row, applied by code to the two
`experiments/harness/run.py decide` outputs (K4 against the raw base and against the control). It reads aggregate
decision files only, no data and no score of a test board.

K4 is kept only if, in BOTH decide files:
  - the base is the file's named base (t3_factor_raw, then t3_factor_eb_ctrl) and the change is t3_factor_eb;
  - both counted folds were decided and the adoption rule keeps it (paired mean gain above 0.3 and ahead in at least
    8 of 10 seeds on both folds, no fold losing more than 0.3 on average);
  - on each fold, the fit checks' parts reproduce the panel (matches_panel);
  - on each fold, the margin sits on identity and direction: (the mean margin of response_magnitude + mutant_state) /
    (the mean score margin) is at most 0.40 (vec.harness.t3_fit_checks: group_margins, score_margin);
  - on each fold, the severity gate differs between the two arms in at most 2 of the 10 seeds.
Any other outcome is "not kept". The verdict records each condition; the row and the main session read it.

Usage (repo root, after both decide runs):
  uv run --frozen --no-sync python -I -B experiments/t3_decomp/k4_verdict.py \\
      --vs-raw experiments/harness/results/t3_k4_v1/vs_raw/decisions_T3_bcatenin.json \\
      --vs-ctrl experiments/harness/results/t3_k4_v1/vs_ctrl/decisions_T3_bcatenin.json \\
      --out experiments/harness/results/t3_k4_v1/verdict
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

FOLDS = ("T3/Gata4->Mab21l2", "T3/Mab21l2->Gata4rep2")
CHANGE = "t3_factor_eb"
BASES = {"vs_raw": "t3_factor_raw", "vs_ctrl": "t3_factor_eb_ctrl"}
MAX_OUTSIDE_SHARE = 0.40  # K4's row: at most this share of the margin on response_magnitude + mutant_state
MAX_GATE_DIFFERS = 2  # K4's row: the severity gate may differ in at most this many seeds


def comparison(dec: dict, base: str) -> dict:
    """One decide file's decision against K4's row: the problems (empty when it passes) and its fold fit readings."""
    problems = []
    if dec.get("base") != base or dec.get("change") != CHANGE:
        problems.append(f"it decides {dec.get('change')} against {dec.get('base')}, not {CHANGE} against {base}")
    if sorted(dec.get("counted", [])) != sorted(FOLDS) or dec.get("need") != 2:
        problems.append(f"counted folds {dec.get('counted')}, need {dec.get('need')}, not both folds, need 2")
    if not dec.get("keep"):
        problems.append(f"the adoption rule does not keep it: {dec.get('reason')}")
    folds = {}
    for fold in FOLDS:
        fc = (dec.get("fit_checks") or {}).get(fold)
        if fc is None:
            problems.append(f"{fold}: no fit checks")
            continue
        if not fc.get("matches_panel", False):
            problems.append(f"{fold}: the fit checks' parts do not reproduce the panel")
        m, score = fc["group_margins"], fc["score_margin"]
        outside = (m["response_magnitude"] + m["mutant_state"]) / score if score > 0 else None
        gate = len(fc["severity_gate_differs"])
        folds[fold] = {"score_margin": score, "group_margins": m, "share_outside_identity_direction": outside,
                       "severity_gate_differs_seeds": gate}
        if outside is None or outside > MAX_OUTSIDE_SHARE:
            problems.append(f"{fold}: the share of the margin outside identity and direction is {outside}, above "
                            f"{MAX_OUTSIDE_SHARE}")
        if gate > MAX_GATE_DIFFERS:
            problems.append(f"{fold}: the severity gate differs in {gate} seeds, above {MAX_GATE_DIFFERS}")
    return {"base": base, "keep_by_adoption_rule": bool(dec.get("keep")), "problems": problems, "folds": folds}


def verdict(raw_decision: dict, ctrl_decision: dict) -> dict:
    """K4's rule over the two decide files: kept only if both comparisons have no problem."""
    parts = {"vs_raw": comparison(raw_decision, BASES["vs_raw"]),
             "vs_ctrl": comparison(ctrl_decision, BASES["vs_ctrl"])}
    return {"kept": not parts["vs_raw"]["problems"] and not parts["vs_ctrl"]["problems"], "parts": parts}


def _load(path: str) -> tuple[dict, str]:
    raw = Path(path).read_bytes()
    obj = json.loads(raw)
    decisions = obj["decisions"]
    if len(decisions) != 1:
        sys.exit(f"{path}: expected one decision (one --change per decide run), found {len(decisions)}")
    return {**decisions[0], "runs": obj.get("runs")}, hashlib.sha256(raw).hexdigest()


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="K4's keep-or-stop rule over its two decide files")
    ap.add_argument("--vs-raw", required=True, help="decisions_T3_bcatenin.json of --base t3_factor_raw")
    ap.add_argument("--vs-ctrl", required=True, help="decisions_T3_bcatenin.json of --base t3_factor_eb_ctrl")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("k4_verdict has already written its outputs; a rerun needs a new TRIES row")

    from vec import provenance

    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    raw, raw_sha = _load(a.vs_raw)
    ctrl, ctrl_sha = _load(a.vs_ctrl)
    v = verdict(raw, ctrl)
    if raw.get("runs") != ctrl.get("runs"):
        v["kept"] = False
        v["runs_mismatch"] = {"vs_raw": raw.get("runs"), "vs_ctrl": ctrl.get("runs")}
    res = {"label": "K4's keep-or-stop rule over its two decide files (experiments/harness/TRIES.md)", **v,
           "inputs": {"vs_raw": {"path": a.vs_raw, "sha256": raw_sha},
                      "vs_ctrl": {"path": a.vs_ctrl, "sha256": ctrl_sha}},
           "stamp": {"commit": commit, "tree_clean": None if unclean is None else not unclean,
                     "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    (out / "k4_verdict.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    text = f"K4 kept: {v['kept']}\n" + "\n".join(
        f"{name}: adoption keep {p['keep_by_adoption_rule']}; problems: {p['problems'] or 'none'}"
        for name, p in v["parts"].items())
    (out / "k4_verdict.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
