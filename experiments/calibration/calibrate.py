#!/usr/bin/env python
"""Emulator calibration against the published reference rows of the validation boards (feat-005, assumption A2).

Pre-registered test (ASSUMPTIONS.md A2, 2026-10-18), run after the team unlocked the validation answers (2026-10-28):
on the validation boards, local scores of the organisers' reference baselines (copy_last / wt_identity,
pseudobulk_shift, shift_transfer, gene_ko) are within 2 points of the published reference rows (mean of 10 seeds),
and the local anchors bracket the published ones (boards/index.json; validation boards only).

Since runbook step 4.6 (2026-10-28) the validation bake is known exactly for the T2 boards (vec.evaluate
protocol="bake": numpy Generator draws at the seed for the target and the reference, veckit split_half(n, 0) with the
first half as the truth, metric seed = the seed). This script also asks the sharper question: at seed 0, do the
bake-protocol anchors equal the published anchors on all five validation boards (T1 and T3 included), and do the
scores equal the published rows?

How the organisers built each baseline prediction is published only as a description (src/vec/baselines.py). A
prediction must respect the board's cell limits, so it cannot be the whole last stage on most boards. The script
scores two constructions, fixed before any comparison (2026-10-31):
  wc          the baseline computed on the 10% working copies of all its input stages (numpy Generator draws at the
              seed, rows sorted, as vec.evaluate.subsample), so copy_last and wt_identity are exactly the floor
  full_delta  the deltas (per cell type on T1/T2, pseudobulk on T3) from the full input stages, applied to the working
              copy of the last stage (T3: of the matched wild type)
Verdict rules (fixed before the run):
  - anchors: at seed 0 under the bake protocol, every computed floor and ceiling of every ranked metric rounds to the
    published value at the published number of decimals ("exact"), else it is within 2% of it or 1e-4 ("close");
  - A2 passes if the anchors are exact or close on every board and, for every published reference row, the better of
    the two constructions is within 2 points of the published score by the mean over seeds 0-9 (bake protocol,
    computed anchors, scorer build "prefix", the live build per runbook step 4.6 and the screened research); it is
    "incomplete" when nothing failed but a board was not run or a published row of a rebuilt model (MODELS) could not
    be scored; rows of models this script does not rebuild (the organisers' ODE and flow models, the controls) are
    listed and do not count;
  - a row is "exact" where its seed-0 score with published anchors matches the published score to two decimals.
The comparison only reads published validation rows; nothing here changes a prediction or a decision rule.

Compliance (docs/DATA_POLICY.md): validation boards only; boards/index.json must be the repo's own (refreshed by
`vec boards update`) and must carry no anchors for any non-validation board; outputs are aggregate numbers (JSON, TXT).
Run on the DGX from the repo root, one board per process, with BLAS threads capped (docs/COMPUTE.md):

  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 uv run python experiments/calibration/calibrate.py \
      run --board T1:val --out runs/calibration
  uv run python experiments/calibration/calibrate.py summarise --rows experiments/calibration/reference_rows.json \
      --runs runs/calibration --out experiments/calibration/results
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from vec import baselines
from vec.boards import REPO_ROOT, VALIDATION_BOARDS, load_boards
from vec.evaluate import Stage, score_board, subsample
from vec.io import read_stage

SEEDS = tuple(range(10))
FRAC = 0.1
BUILD = "prefix"
PROTOCOL = "bake"
CONSTRUCTIONS = ("wc", "full_delta")
MODELS = ("copy_last", "wt_identity", "pseudobulk_shift", "shift_transfer", "gene_ko")  # the reference rows rebuilt
TOL_POINTS = 2.0

# Input files per board (data/raw names from docs/DATA_NOTES.md). The T3 truth is the Gata4 replicate whose 10%
# working copy halves to the board's truth_cells; the script checks the counts of every board before scoring.
BOARD_FILES = {
    "T1:val": {"target": "T1/E10.5_RNA.h5ad", "ref": "T1/E9.5_RNA.h5ad", "prev": "T1/E8.5_RNA.h5ad"},
    "T2:heart:val_extrap": {"target": "T2/E10.5.h5ad", "ref": "T2/E9.5.h5ad", "prev": "T2/E8.75.h5ad"},
    "T2:heart:val_interp": {"target": "T2/E8.5.h5ad", "ref": "T2/E8.25_late.h5ad", "prev": None},
    "T2:embryo:val_interp": {"target": "T2/E7.5.h5ad", "ref": "T2/E7.25.h5ad", "prev": "T2/E6.75.h5ad"},
    "T3:gata4": {"target": "T3/E8.75_gata4_ko_2.h5ad", "ref": "T2/E8.75.h5ad", "ko_train": "T3/E9.5_mab21l2_ko.h5ad",
                 "wt_train": "T2/E9.5.h5ad", "target_gene": "Gata4"},
}


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _decimals(x: float) -> int:
    s = repr(float(x))
    if "e" in s:
        mant, exp = s.split("e")
        return max(0, len(mant.split(".")[-1]) - int(exp))
    return len(s.split(".")[-1])


def anchor_match(computed: float | None, published: float) -> str:
    if computed is None or not math.isfinite(computed):
        return "missing"
    if round(computed, _decimals(published)) == round(published, _decimals(published)):
        return "exact"
    if abs(computed - published) <= max(1e-4, 0.02 * abs(published)):
        return "close"
    return "off"


def counts_ok(board, N_target: int, N_ref: int) -> dict:
    nt = max(2, int(round(N_target * FRAC)))
    nr = max(2, int(round(N_ref * FRAC)))
    return {"N_target": N_target, "N_ref": N_ref, "target_wc": nt, "truth_half": nt // 2, "ref_wc": nr,
            "truth_cells": board.truth_cells, "ref_cells": board.ref_cells,
            "ok": board.truth_cells in (nt // 2, nt - nt // 2) and board.ref_cells == nr}


def _type_shift(prev: Stage, last: Stage, apply_to: Stage) -> Stage:
    """pseudobulk_shift with the per-type deltas measured on (prev, last) and added to `apply_to` (clip at 0)."""
    X = apply_to.X.copy()
    prev_types, last_types = set(prev.celltype), set(last.celltype)
    for c in np.unique(apply_to.celltype):
        if c not in prev_types or c not in last_types:
            continue
        sel = apply_to.celltype == c
        X[sel] += last.X[last.celltype == c].mean(0) - prev.X[prev.celltype == c].mean(0)
    return Stage(X=np.clip(X, 0.0, None), coords=None if apply_to.coords is None else apply_to.coords.copy())


def predictions(task: str, stages: dict, files: dict, genes: list[str], seed: int) -> dict[str, dict[str, Stage]]:
    """model -> construction -> prediction, at one seed."""
    w = {k: subsample(s, FRAC, seed) for k, s in stages.items() if k != "target"}
    out: dict[str, dict[str, Stage]] = {}
    if task in ("T1", "T2"):
        floor = baselines.copy_last(w["ref"])
        out["copy_last"] = {"wc": floor, "full_delta": floor}
        if "prev" in stages:
            out["pseudobulk_shift"] = {
                "wc": baselines.pseudobulk_shift(w["prev"], w["ref"]),
                "full_delta": _type_shift(stages["prev"], stages["ref"], w["ref"]),
            }
    else:
        floor = baselines.wt_identity(w["ref"])
        out["wt_identity"] = {"wc": floor, "full_delta": floor}
        st_wc = baselines.shift_transfer(w["ko_train"], w["wt_train"], w["ref"])
        st_full = baselines.shift_transfer(stages["ko_train"], stages["wt_train"], w["ref"])
        out["shift_transfer"] = {"wc": st_wc, "full_delta": st_full}
        gk = {}
        for name, base in (("wc", st_wc), ("full_delta", st_full)):
            X = base.X.copy()
            if files["target_gene"] in genes:
                X[:, genes.index(files["target_gene"])] = 0.0
            gk[name] = Stage(X=X, coords=base.coords)
        out["gene_ko"] = gk
    return out


def check_index(boards: dict) -> None:
    leaking = [k for k, b in boards.items() if k not in VALIDATION_BOARDS and b.anchors]
    if leaking:  # test-board anchors encode hidden targets and are never read
        sys.exit(f"refused: {len(leaking)} non-validation board(s) in boards/index.json carry anchors")


def run_board(key: str, data: Path, out: Path) -> Path:
    boards = load_boards()
    check_index(boards)
    if key not in VALIDATION_BOARDS:
        sys.exit(f"refused: {key} is not a validation board")
    board, files = boards[key], BOARD_FILES[key]
    t0 = time.time()

    def load(rel: str) -> Stage:
        return Stage.from_anndata(read_stage(data / rel, board.genes), board.genes, need_coords=board.needs_coords)

    stages = {k: load(v) for k, v in files.items() if k != "target_gene" and v}
    counts = counts_ok(board, stages["target"].n, stages["ref"].n)
    if not counts["ok"]:
        sys.exit(f"refused: {key} working-copy counts do not match the board: {counts}")

    rows, anchors_seed0 = [], None
    for seed in SEEDS:
        preds = predictions(board.task, stages, files, board.genes, seed)
        for model, by_c in preds.items():
            for construction, pred in by_c.items():
                if construction == "full_delta" and model in ("copy_last", "wt_identity"):
                    continue  # identical to wc
                r = score_board(board, pred, stages["target"], stages["ref"], seed=seed, anchors="computed",
                                build=BUILD, protocol=PROTOCOL)
                row = {"model": model, "construction": construction, "seed": seed, "anchors": "computed",
                       "score": r.score, "skills": r.skills, "group_skill": r.group_skill, "pred_cells": pred.n}
                rows.append(row)
                if seed == 0:
                    if anchors_seed0 is None:
                        anchors_seed0 = r.anchors
                    p = score_board(board, pred, stages["target"], stages["ref"], seed=0, anchors="published",
                                    build=BUILD, protocol=PROTOCOL)
                    rows.append({**row, "anchors": "published", "score": p.score, "skills": p.skills,
                                 "group_skill": p.group_skill})
        print(f"{key}: seed {seed} done ({time.time() - t0:.0f} s)", flush=True)

    anchor_cmp = {m: {k: {"published": board.anchors[m][k], "computed": anchors_seed0.get(m, {}).get(k),
                          "match": anchor_match(anchors_seed0.get(m, {}).get(k), board.anchors[m][k])}
                      for k in ("floor", "ceiling")} for m in board.anchors}
    res = {"board": key, "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), "commit": _git("rev-parse", "HEAD"),
           "tree_clean": _git("status", "--porcelain", "--untracked-files=no") == "", "build": BUILD,
           "protocol": PROTOCOL, "seeds": SEEDS, "frac": FRAC, "counts": counts, "anchors_seed0": anchor_cmp,
           "rows": rows,
           "seconds": round(time.time() - t0, 1)}
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{key.replace(':', '__')}.json"
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}")
    return path


def summarise(rows_file: Path, runs: Path, out: Path) -> dict:
    all_rows = json.loads(rows_file.read_text())["rows"]
    published = [r for r in all_rows if r["model"] in MODELS]
    out_of_scope = sorted({(r["board"], r["model"]) for r in all_rows if r["model"] not in MODELS})
    per_board = {p.stem.replace("__", ":"): json.loads(p.read_text()) for p in sorted(runs.glob("*.json"))}
    missing = sorted(VALIDATION_BOARDS - set(per_board))
    comparisons, lines = [], []
    anchors_ok = True
    for key, res in per_board.items():
        bad = [(m, k, v["match"]) for m, a in res["anchors_seed0"].items() for k, v in a.items()
               if v["match"] not in ("exact", "close")]
        exact = all(v["match"] == "exact" for a in res["anchors_seed0"].values() for v in a.values())
        anchors_ok &= not bad
        lines.append(f"## {key}: counts ok {res['counts']['ok']}; anchors at seed 0 "
                     f"{'all exact' if exact else ('exact or close' if not bad else f'OFF {bad}')}; "
                     f"commit {res['commit'][:7]} clean {res['tree_clean']}")
        for m, a in res["anchors_seed0"].items():
            lines.append("  " + m + ": " + "; ".join(
                f"{k} published {v['published']} computed {v['computed']:.6g} [{v['match']}]"
                if v["computed"] is not None else f"{k} published {v['published']} computed - [missing]"
                for k, v in a.items()))
        for pub in [r for r in published if r["board"] == key]:
            model = pub["model"]
            mine = [r for r in res["rows"] if r["model"] == model]
            if not mine:
                comparisons.append({"board": key, "model": model, "published": pub["score"], "status": "not scored"})
                lines.append(f"  {model}: published {pub['score']} - not scored here")
                continue
            best = None
            for c in CONSTRUCTIONS:
                comp = [r["score"] for r in mine if r["construction"] == c and r["anchors"] == "computed"]
                row0 = [r for r in mine if r["construction"] == c and r["anchors"] == "published"]
                pub0 = [r["score"] for r in row0]
                if not comp:
                    continue
                d = {"board": key, "model": model, "construction": c, "published": pub["score"],
                     "mean_10": float(np.mean(comp)), "sd_10": float(np.std(comp, ddof=1)) if len(comp) > 1 else 0.0,
                     "seed0_published_anchors": pub0[0] if pub0 else None}
                d["diff_mean"] = d["mean_10"] - pub["score"]
                d["exact"] = pub0 != [] and round(pub0[0], 2) == round(float(pub["score"]), 2)
                # diagnostics only: seed-0 skills (0-100 scale, published anchors) against the published skills
                d["skill_diffs_seed0"] = {m: 100.0 * row0[0]["skills"][m] - v for m, v in pub.get("skills", {}).items()
                                          if row0 and row0[0]["skills"].get(m) is not None
                                          and math.isfinite(row0[0]["skills"][m])}
                comparisons.append(d)
                worst = max(d["skill_diffs_seed0"].items(), key=lambda kv: abs(kv[1]), default=None)
                lines.append(f"  {model} [{c}]: published {pub['score']}; mean of 10 seeds {d['mean_10']:.2f} "
                             f"+/- {d['sd_10']:.2f} (diff {d['diff_mean']:+.2f}); seed 0 with published anchors "
                             f"{d['seed0_published_anchors']:.2f}{' [exact]' if d['exact'] else ''}"
                             + (f"; largest seed-0 skill gap {worst[0]} {worst[1]:+.2f}" if worst else ""))
                if best is None or abs(d["diff_mean"]) < abs(best["diff_mean"]):
                    best = d
            if best is not None:
                best["best"] = True
    rows_ok = all(abs(c["diff_mean"]) <= TOL_POINTS for c in comparisons if c.get("best"))
    unscored = [c for c in comparisons if c.get("status") == "not scored"]
    if not (anchors_ok and rows_ok):
        verdict = "fail"
    else:  # nothing failed; a board not run or a published row of a reproduced model not scored leaves it open
        verdict = "incomplete" if missing or unscored else "pass"
    summary = {"utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), "verdict_A2": verdict,
               "anchors_ok": anchors_ok, "rows_within_2_points": rows_ok, "boards_missing": missing,
               "rows_not_scored": [[c["board"], c["model"]] for c in unscored],
               "rows_out_of_scope": [list(x) for x in out_of_scope], "comparisons": comparisons}
    out.mkdir(parents=True, exist_ok=True)
    (out / "calibration_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    head = [f"A2 verdict: {verdict} (anchors ok {anchors_ok}; best construction within {TOL_POINTS} points on every "
            f"published row {rows_ok}; boards missing {missing}; rows not scored {summary['rows_not_scored']}; "
            f"{len(out_of_scope)} published rows of models this script does not rebuild)"]
    (out / "calibration_summary.txt").write_text("\n".join(head + lines) + "\n")
    print("\n".join(head + lines))
    return summary


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="score the reference baselines on one validation board")
    r.add_argument("--board", choices=sorted(BOARD_FILES), required=True)
    r.add_argument("--data", type=Path, default=REPO_ROOT / "data" / "raw")
    r.add_argument("--out", type=Path, default=REPO_ROOT / "runs" / "calibration")
    s = sub.add_parser("summarise", help="compare the board runs with the published rows; write the verdict")
    s.add_argument("--rows", type=Path, required=True)
    s.add_argument("--runs", type=Path, required=True)
    s.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "summarise":
        summarise(a.rows, a.runs, a.out)
    else:
        run_board(a.board, a.data, a.out)


if __name__ == "__main__":
    main()
