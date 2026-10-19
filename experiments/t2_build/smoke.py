#!/usr/bin/env python
"""Smoke test of vec.t2_build on released data: copy-like T2 submissions through the count rule, the sampler, the
guard, the builder and the emulator under both scorer builds. DGX only; aggregates only.

Jobs. Board specs come from the repo's boards/index.json (refreshed only through vec.boards.update_boards); only the
gene panel and the cell limits are used, never an anchor.
- heart: E8.75 cells on the T2:heart:val_extrap spec, scored against the released stand-in target E9.5 (reference
  E8.75).
- embryo: E7.25 cells on the T2:embryo:val_interp spec, scored against the released stand-in target E8.0 (reference
  E7.25).
The stand-in targets are released training stages: no validation answer, test input or held-out stage is read, and
the board's own truth_cells and max_cells belong to its held-out stage (they set only the count, as on the real board).

Per job: load the source with drop_duplicates=True (a prediction input); n_cells (max_cells); sample_cells (rows repeat
only when the count exceeds the pool, and repeated copies are jittered apart); write the file with build_submission (it
must pass vec.io.check_submission); cloud_guard against the pool; read the file back as `vec eval` does and score it
with vec.evaluate.score_board_seeds under builds "main" and "prefix", anchors computed on the stand-in. Two labelled
controls: a contiguous mid-file block of pool rows (the guard must flag it) and the copy rescaled to the stand-in
target's RMS radius (checks rescale_to_rms through scale_log_ratio and the guard's expected ratio; an oracle on a
released stage, not a method).

Under "prefix" the parity coin also acts on the computed floor and ceiling, so a stand-in can show inverted Dice anchors
(a mirror-mode ceiling below the floor; the Dice skill of that seed is then not meaningful). They are recorded with
their seeds; the inverted-anchor check applies to "main" only.

Outputs (aggregates only: counts, records, guard numbers, scores and skills): <out>/smoke_summary.{json,txt}, which
stay on the DGX with the submission files, and <out>/smoke_summary_shareable.{json,txt}, an allow-list of the fields
from which the validation boards' published cell limits cannot be recovered without the released data (`_shareable`;
held-out-stage sizes stay out of committed files). Bring back only the shareable pair.

Usage (DGX, repo root, in tmux with a log under runs/logs/):
  OMP_NUM_THREADS=8 POT_BACKEND_DISABLE_JAX=1 XLA_PYTHON_CLIENT_PREALLOCATE=false \
    uv run python experiments/t2_build/smoke.py --out runs/t2_build/smoke
Exit status: 0 when every hard check passes, 1 when one fails, 2 when an input file is missing.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import typer
from rich.console import Console
from rich.table import Table

from vec.boards import RANKED_METRICS, REPO_ROOT, VALIDATION_BOARDS, get_board
from vec.evaluate import SCORER_BUILDS, Stage, score_board_seeds
from vec.io import COORD_KEY, read_stage
from vec.prep import load_stage
from vec.t2_build import build_submission, cloud_guard, n_cells, rescale_to_rms, rms_radius, sample_cells

app = typer.Typer(add_completion=False)
console = Console(emoji=False, highlight=False)
SHAPE_BUILD_FREE = ("d2_shape", "scale_log_ratio")  # never call _canonicalise: identical under both builds
ORACLE_MAX_ABS_SCALE = 0.05  # |scale_log_ratio| of the copy rescaled to the stand-in's RMS (truth half: a 5% draw)
SHAREABLE_NOTE = ("allow-list (smoke._shareable): only fields from which a validation board's cell limits cannot be "
                  "recovered without the released data, so no count, sample or block-control record, guard reference "
                  "statistic, file size or handedness count. Threat model set by the lead: nothing about a limit may "
                  "be recoverable from these files alone; data access is not assumed, so the rounded guard ratios and "
                  "scores of the submitted cells stay")


@dataclass(frozen=True)
class Job:
    board: str
    source: str  # prediction input, loaded with duplicates dropped
    target: str  # released stand-in target, loaded as released (the emulator path)
    ref: str  # DE reference and copy_last floor, loaded as released


JOBS = {
    "heart": Job("T2:heart:val_extrap", "T2/E8.75.h5ad", "T2/E9.5.h5ad", "T2/E8.75.h5ad"),
    "embryo": Job("T2:embryo:val_interp", "T2/E7.25.h5ad", "T2/E8.0.h5ad", "T2/E7.25.h5ad"),
}


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _sha256(path: Path) -> str:
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def _eval_stage(path: Path, genes: list[str]) -> Stage:
    """A file as the emulator reads it (`vec eval`): panel columns, float32 coordinates, labels if any."""
    return Stage.from_anndata(read_stage(path, genes), genes, need_coords=True)


def _summarise(res: dict) -> dict:
    """Aggregates of score_board_seeds output: score, mean skills, group skills and ranked metric values."""
    per = res["per_seed"]

    def mean(key: str, m: str) -> float:
        vals = [p[key].get(m) for p in per]
        vals = [float(v) for v in vals if v is not None and np.isfinite(v)]
        return float(np.mean(vals)) if vals else float("nan")

    return {
        "build": res["build"], "score_mean": res["score_mean"], "score_sd": res["score_sd"],
        "scores": [p["score"] for p in per], "seeds": [p["seed"] for p in per],
        "skills_mean": {m: mean("skills", m) for m in RANKED_METRICS["T2"]},
        "metrics_mean": {m: mean("metrics", m) for m in RANKED_METRICS["T2"]},
        "metrics_per_seed": {m: [p["metrics"].get(m) for p in per] for m in SHAPE_BUILD_FREE + ("occupancy_dice",)},
        "group_skill_mean": {g: float(np.mean([p["group_skill"][g] for p in per])) for g in per[0]["group_skill"]},
        "truth_cells_emulated": per[0]["extras"]["truth_cells"], "pred_cells": per[0]["extras"]["pred_cells"],
        "inverted_anchors": sorted({m for p in per for m in p["extras"].get("inverted_anchors", [])}),
        "inverted_anchor_seeds": [p["seed"] for p in per if p["extras"].get("inverted_anchors")],
        "handedness": res.get("handedness"),
    }


def _score(board, pred_path: Path, target: Stage, ref: Stage, seeds: tuple[int, ...]) -> dict:
    pred = _eval_stage(pred_path, board.genes)
    out = {}
    for build in SCORER_BUILDS:
        t0 = time.time()
        out[build] = _summarise(score_board_seeds(board, pred, target, ref, seeds=seeds, anchors="computed",
                                                  build=build))
        out[build]["seconds"] = round(time.time() - t0, 1)
    return out


def _run_job(name: str, job: Job, data: Path, out: Path, seeds: tuple[int, ...], sample_seed: int,
             score: bool) -> tuple[dict, dict]:
    board = get_board(job.board)
    if job.board not in VALIDATION_BOARDS:
        raise ValueError(f"{job.board} is not a validation board")
    paths = {k: data / getattr(job, k) for k in ("source", "target", "ref")}
    t0 = time.time()
    pool = load_stage(paths["source"], genes=board.genes, lineage=False, drop_duplicates=True)
    P = np.asarray(pool.obsm[COORD_KEY], dtype=np.float64)[:, :3]
    qc = pool.uns.get("qc", {})
    n, count_record = n_cells(board, pool.n_obs)
    s = sample_cells(pool, n, seed=sample_seed)
    genes = list(pool.var_names)
    files = {"copy": out / f"{name}_copy.h5ad"}
    _, rep = build_submission(board, pool.X[s.rows], s.coords, genes=genes, out=files["copy"])
    guard = {"copy": cloud_guard(s.coords, P, seed=sample_seed).as_dict()}

    # control 1: a contiguous mid-file block of pool rows (rows are stored section by section)
    m = min(n, pool.n_obs // 2)
    start = (pool.n_obs - m) // 2
    block = np.arange(start, start + m)
    guard["block_control"] = cloud_guard(P[block], P, seed=sample_seed).as_dict()
    files["block_control"] = out / f"{name}_block_control.h5ad"
    build_submission(board, pool.X[block], P[block], genes=genes, out=files["block_control"])

    # control 2: the copy rescaled to the stand-in target's RMS radius (an oracle on a released stage)
    target = _eval_stage(paths["target"], board.genes)
    target_rms = rms_radius(target.coords)
    scaled = rescale_to_rms(s.coords, target_rms)
    expected = rms_radius(scaled) / rms_radius(P)
    guard["oracle_size_control"] = cloud_guard(scaled, P, seed=sample_seed, expected_rms_ratio=expected).as_dict()
    files["oracle_size_control"] = out / f"{name}_oracle_size_control.h5ad"
    build_submission(board, pool.X[s.rows], scaled, genes=genes, out=files["oracle_size_control"])
    build_seconds = round(time.time() - t0, 1)

    scores = {}
    if score:
        ref = _eval_stage(paths["ref"], board.genes)
        scores = {cond: _score(board, f, target, ref, seeds) for cond, f in files.items()}

    result = {
        "board": job.board, "board_spec": {"min_cells": board.min_cells, "max_cells": board.max_cells,
                                           "truth_cells": board.truth_cells, "n_genes": board.n_genes},
        "inputs": {k: {"file": getattr(job, k), "sha256": _sha256(p)} for k, p in paths.items()},
        "pool": {"n_input": int(qc.get("n_input", pool.n_obs)), "n_dropped_duplicates": int(qc.get("n_dropped", 0)),
                 "n_pool": int(pool.n_obs), "rms_radius": rms_radius(P)},
        "count": count_record, "sample": s.record,
        "check_submission": {"ok": rep.ok, "errors": rep.errors, "warnings": rep.warnings, "n_cells": rep.n_cells},
        "block_control": {"n": int(m), "first_row": int(start)},
        "oracle_size_control": {"target_rms_radius": target_rms, "scale": target_rms / rms_radius(s.coords)},
        "guard": guard, "files_mb": {k: round(f.stat().st_size / 2 ** 20, 1) for k, f in files.items()},
        "scores": scores, "build_seconds": build_seconds,
    }
    return result, _checks(board, n, result, score)


def _checks(board, n: int, r: dict, score: bool) -> dict:
    """Hard checks of one job (True = pass)."""
    g = r["guard"]
    c = {
        "copy_count_is_max_cells": n == board.max_cells and r["check_submission"]["n_cells"] == n,
        "copy_passes_check_submission": r["check_submission"]["ok"],
        "copy_file_has_no_repeated_positions": not any("exact coordinates" in w
                                                       for w in r["check_submission"]["warnings"]),
        "copy_passes_guard": g["copy"]["ok"],
        "block_control_flagged_as_block": any(f.startswith("block") for f in g["block_control"]["flags"]),
        "oracle_size_control_passes_guard": g["oracle_size_control"]["ok"],
    }
    if score:
        s = r["scores"]
        c["scores_finite_both_builds"] = all(np.isfinite(s[cond][b]["score_mean"]) for cond in s for b in s[cond])
        c["build_free_metrics_identical_across_builds"] = all(
            s[cond]["main"]["metrics_per_seed"][m] == s[cond]["prefix"]["metrics_per_seed"][m]
            for cond in s for m in SHAPE_BUILD_FREE)
        c["oracle_scale_log_ratio_near_zero"] = all(
            abs(s["oracle_size_control"][b]["metrics_mean"]["scale_log_ratio"]) < ORACLE_MAX_ABS_SCALE
            for b in SCORER_BUILDS)
        # Under "prefix" the parity coin also acts on the computed floor and ceiling, so a mirror-mode ceiling can
        # fall below the floor (inverted anchors): recorded, not a failure. Under "main" it must not happen.
        c["no_inverted_anchors_main_build"] = not any(s[cond]["main"]["inverted_anchors"] for cond in s)
    return c


def _table(name: str, r: dict) -> Table:
    t = Table("condition", "build", "score", "dice", "d2", "scale", "nbhd", "guard", title=f"{name}: {r['board']}")
    for cond, by_build in r["scores"].items():
        flags = r["guard"][cond]["flags"]
        verdict = "ok" if not flags else ",".join(f.split(":")[0] for f in flags)
        for b, v in by_build.items():
            sk = v["skills_mean"]
            t.add_row(cond, b, f"{v['score_mean']:.2f} +/- {v['score_sd']:.2f}", f"{sk['occupancy_dice']:.3f}",
                      f"{sk['d2_shape']:.3f}", f"{sk['scale_log_ratio']:.3f}", f"{sk['neighborhood_mmd']:.3f}", verdict)
    return t


def _guard_view(g: dict) -> dict:
    v = {"ok": bool(g["ok"]), "cover_ratio": round(float(g["cover_ratio"]), 2),
         "rms_ratio": round(float(g["rms_ratio"]), 3)}
    if g.get("expected_rms_source") == "argument":  # the oracle's; a random reference's is a reference statistic
        v["expected_rms_ratio"] = round(float(g["expected_rms_ratio"]), 3)
    return v


def _score_view(v: dict) -> dict:
    return {"score_mean": round(float(v["score_mean"]), 2), "score_sd": round(float(v["score_sd"]), 2),
            "skills_mean": {m: round(float(x), 3) for m, x in v["skills_mean"].items()},
            "inverted_anchors": list(v["inverted_anchors"]),
            "inverted_anchor_seeds": [int(s) for s in v["inverted_anchor_seeds"]]}


def _shareable(summary: dict) -> dict:
    """The committed form of a summary: an allow-list of the fields from which a validation board's published cell
    limits cannot be recovered without the released data (committed docs and results keep held-out-stage sizes out:
    docs/DATA_POLICY.md). Kept: the run's time, commit, emulator seeds and checks; per job the board, the count rule's
    name (an override's argument can be a count), the pool size and dropped records, the input files and their sha256,
    per guard condition ok, the coverage ratio (2 decimals) and the RMS ratio (3 decimals), per scored condition and
    build the score mean and SD (2 decimals), the mean skills (3 decimals) and the inverted anchors, and
    check_submission's verdict and warnings. The expected RMS ratio (3 decimals) stays only where it was passed in (the
    oracle's): a random reference's is a reference statistic, and the copy and the block control draw theirs with one
    seed, so equal counts would show as equal values. Every other field is dropped, fields added to the summary later
    included."""
    jobs = {}
    for name, r in summary["jobs"].items():
        jobs[name] = {
            "board": r["board"],
            "count_rule": r["count"]["rule"].split()[0],
            "pool": {k: int(r["pool"][k]) for k in ("n_pool", "n_dropped_duplicates")},
            "inputs": {k: {"file": f["file"], "sha256": f["sha256"]} for k, f in r["inputs"].items()},
            "guard": {cond: _guard_view(g) for cond, g in r["guard"].items()},
            "scores": {cond: {b: _score_view(v) for b, v in bb.items()} for cond, bb in r["scores"].items()},
            "check_submission": {"ok": bool(r["check_submission"]["ok"]),
                                 "warnings": list(r["check_submission"]["warnings"])},
        }
    return {"generated_utc": summary["generated_utc"], "git_commit": summary["git_commit"],
            "seeds": [int(s) for s in summary["seeds"]], "shareable": SHAREABLE_NOTE, "jobs": jobs,
            "checks": {name: {k: bool(v) for k, v in c.items()} for name, c in summary["checks"].items()},
            "all_checks_pass": bool(summary["all_checks_pass"])}


def _lines(summary: dict, seeds: str) -> list[str]:
    """Text report of a shareable summary (`_shareable` output)."""
    lines = [f"t2_build smoke {summary['generated_utc']} commit {summary['git_commit']} seeds {seeds}"]
    for name, r in summary["jobs"].items():
        lines.append(f"\n[{name}] {r['board']}: count rule {r['count_rule']}, pool {r['pool']['n_pool']} after "
                     f"dropping {r['pool']['n_dropped_duplicates']} synthetic records")
        for cond, gd in r["guard"].items():
            expected = f"{gd['expected_rms_ratio']:.3f}" if "expected_rms_ratio" in gd else "from the random reference"
            lines.append(f"  guard {cond}: ok={gd['ok']} cover ratio {gd['cover_ratio']:.2f}, rms ratio "
                         f"{gd['rms_ratio']:.3f} (expected {expected})")
        for cond, by_build in r["scores"].items():
            for b, v in by_build.items():
                sk = v["skills_mean"]
                inv = (f"  inverted anchors {v['inverted_anchors']} in seeds {v['inverted_anchor_seeds']}"
                       if v["inverted_anchors"] else "")
                lines.append(f"  score {cond:20s} {b:6s} {v['score_mean']:6.2f} +/- {v['score_sd']:.2f}  dice "
                             f"{sk['occupancy_dice']:.3f} d2 {sk['d2_shape']:.3f} scale {sk['scale_log_ratio']:.3f} "
                             f"nbhd {sk['neighborhood_mmd']:.3f}{inv}")
        cs = r["check_submission"]
        lines.append(f"  check_submission: ok={cs['ok']}, warnings: {'; '.join(cs['warnings']) or 'none'}")
        checks = summary["checks"][name]
        lines.append("  checks: " + ", ".join(f"{k}={'pass' if v else 'FAIL'}" for k, v in checks.items()))
    lines.append(f"\nall checks pass: {summary['all_checks_pass']}")
    return lines


def _count_lines(summary: dict) -> list[str]:
    """The counts behind each job, for the full text report only (they depend on the boards' cell limits)."""
    lines = ["\ncounts (full report only)"]
    for name, r in summary["jobs"].items():
        c, s, b = r["count"], r["sample"], r["block_control"]
        lines.append(f"  [{name}] n={c['n']} (floor {c['floor']} from {c['floor_source']}), repeated rows "
                     f"{s['n_repeated_rows']}, jittered {s['n_jittered']} (max shift {s['max_shift']:.3g}, budget "
                     f"{s['jitter_radius_max']:.3g}), block control {b['n']} rows from row {b['first_row']}")
    return lines


def _write(out: Path, summary: dict, seeds: str) -> list[str]:
    """Write <out>/smoke_summary.{json,txt} (the full summary, stays on the DGX) and
    <out>/smoke_summary_shareable.{json,txt}; returns the shareable text lines."""
    shareable = _shareable(summary)
    lines = _lines(shareable, seeds)
    for suffix, data_, text in (("", summary, lines + _count_lines(summary)), ("_shareable", shareable, lines)):
        (out / f"smoke_summary{suffix}.json").write_text(json.dumps(data_, indent=2, default=float) + "\n")
        (out / f"smoke_summary{suffix}.txt").write_text("\n".join(text) + "\n")
    return lines


@app.command()
def main(
    data: Path = typer.Option(Path("data/raw"), help="released data root (holds T2/)"),
    out: Path = typer.Option(Path("runs/t2_build/smoke"), help="output folder (stays on the DGX)"),
    jobs: str = typer.Option("heart,embryo", help="comma-separated jobs: heart, embryo"),
    seeds: str = typer.Option("0,1,2", help="emulator subsampling seeds"),
    sample_seed: int = typer.Option(0, help="seed for sample_cells and the guard's reference draws"),
    score: bool = typer.Option(True, help="score with the emulator (both builds)"),
):
    """Build copy-like T2 submissions with vec.t2_build from released stages and score them on stand-in targets."""
    names = [j.strip() for j in jobs.split(",") if j.strip()]
    unknown = [j for j in names if j not in JOBS]
    if unknown:
        raise typer.BadParameter(f"unknown jobs {unknown}; known: {sorted(JOBS)}")
    missing = [str(data / getattr(JOBS[j], k)) for j in names for k in ("source", "target", "ref")
               if not (data / getattr(JOBS[j], k)).exists()]
    if missing:
        console.print(f"missing input files: {sorted(set(missing))}")
        raise typer.Exit(2)
    out.mkdir(parents=True, exist_ok=True)
    seed_tuple = tuple(int(x) for x in seeds.split(","))
    summary = {"generated_utc": datetime.now(UTC).isoformat(timespec="seconds"), "git_commit": _git_commit(),
               "host": platform.node(), "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
               "seeds": list(seed_tuple), "sample_seed": sample_seed, "anchors": "computed (stand-in targets)",
               "jobs": {}, "checks": {}}
    for name in names:
        console.print(f"[{name}] building and scoring ({JOBS[name].board})")
        t0 = time.time()
        res, checks = _run_job(name, JOBS[name], data, out, seed_tuple, sample_seed, score)
        res["seconds"] = round(time.time() - t0, 1)
        summary["jobs"][name], summary["checks"][name] = res, checks
        console.print(f"[{name}] n={res['count']['n']} pool={res['pool']['n_pool']} "
                      f"jittered={res['sample']['n_jittered']} guard copy ok={res['guard']['copy']['ok']} "
                      f"block flagged={checks['block_control_flagged_as_block']} ({res['seconds']} s)")
        if score:
            console.print(_table(name, res))
    passed = all(all(c.values()) for c in summary["checks"].values())
    summary["all_checks_pass"] = passed
    console.print("\n".join(_write(out, summary, seeds)))
    raise typer.Exit(0 if passed else 1)


if __name__ == "__main__":
    app()
