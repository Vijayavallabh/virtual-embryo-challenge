#!/usr/bin/env python
"""E12d (experiments/harness/TRIES.md; a diagnostic, not a try): the embryo s2 candidate fgw_shape0_aligned against
fgw_shape0 at the test's own mapped fraction (t = 1/2) with the test pair's mirrored partners, on the test pair's own
section halves, in a pseudo-truth world that keeps E7.5's geometry (the embryo design's C1,
experiments/research_2026_11_sota/notes/design/embryo.md, narrowed by E12's sharpened spec to that world at t' = 1/2).

Target-free: it reads only the test's two inputs, E7.5 (unlocked, A20) and E8.0, on the test board's gene panel (the
panel alone, through vec.boards, as N1 read it), no target and no score of one, and it prints and writes no cell count
(a stage's, a half's or a draw's). It is the only read at the test's mapped fraction: the backtests map a third of
their cells with partners random along the symmetric direction (E7.5) or in only five of ten seeds mirrored (E7.25),
the test pair maps half of them through mirrored partners in every draw (N1).

Method, fixed in row E12d before it runs:
  1. Each bracket's cells are split into halves A and B by alternating snapped z sections (vec.spatial.snap_z; even
     levels are A, odd levels are B).
  2. A pseudo-truth world from half A keeps E7.5's geometry at t' = 1/2: vec.t2_build.interpolate_clouds(t = 1/2,
     t_pos = 0) of A's E7.5 and E8.0 cells, E7.5's real cells where they are and A's E8.0 cells placed onto E7.5's
     geometry by a coupling. The world is built twice: the placement by FGW's coupling (the incumbent's, mirrored at
     test) and by the aligned one (the candidate's), which each favour their own candidate; both are read, at the size
     rule on A's two stages (PCHIP of log RMS, as the backtests' size rule).
  3. The predictions from half B are the harness's own builders: fgw_shape0 and fgw_shape0_aligned (t_pos = 0, same
     mixture of B's cells, the moves from FGW's or the aligned coupling), at the build seeds the harness uses.
  4. Each is scored by the calibrated emulator (vec.evaluate.score_board, protocol "bake", build "main", computed
     anchors, B's E7.5 as the floor and reference), over seeds 0-9 of the world, the same seed for its draws, for both
     worlds. The worlds' own couplings and draws take seeds 2,000,000 + s, apart from the builds'.
  5. Reported per world, for the candidate against fgw_shape0: the paired mean gain, the seeds ahead, the shape-skill
     changes; and per seed the partners' cross-correlation (row E10's measure) of both couplings and the thinnest axis
     of each candidate's pool (reported, not judged).

Reading, fixed in row E12d (it decides nothing by itself; its numbers go to the captain with E12's outcome):
  "supports" if the candidate's paired mean gain is above +0.3 with at least 8 seeds ahead in the aligned-placed world
  (the adoption rule's form) and not below -0.3 in the FGW-placed world (the adoption rule's loss bound);
  "against" if its paired mean gain is below -0.3 in the aligned-placed world (the one that favours it);
  otherwise "no clear read". A partial run gives no reading.
Degenerate anchors. In a world that keeps E7.5's geometry the floor (a copy of B's E7.5) can sit as close to the truth's
shape as the ceiling does, so on a shape term a seed's computed anchors can be inverted (the harness refuses such a
board). Such a seed gives no score for that world and is counted; a world with fewer than 8 usable seeds gives "no clear
read". The candidates' raw shape-metric changes (all seeds, degenerate or not) are reported beside the scores.

Outputs (--out, default experiments/harness/results/e12d_halves/, refused if it holds anything): e12d_halves.json and
.txt, each stamped (commit, clean-tree flag, UTC, the two input files' sha256 and the panel's) and with a vec.provenance
build stamp. Scores and skills per seed only, no per-cell data; a world's truth is a mixture of the test pair's own
halves, so nothing target-derived is written.

Usage (DGX, repo root, BLAS threads capped as docs/COMPUTE.md says; an estimated 15 to 30 minutes of one process at 16
threads, from E11's 30-odd seconds per seed row, here four couplings and two scoring rounds per seed):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
    uv run --frozen --no-sync python -I -B experiments/harness/diag_e12_halves.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np

from vec import adoption, harness, provenance, t2_build
from vec.boards import load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, panel_metrics, score_board, subsample, working_copies
from vec.harness import BUILD, BUILD_SEED_OFFSET, PROTOCOL
from vec.io import COORD_KEY, dense
from vec.prep import load_stage
from vec.spatial import snap_z

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import diag_emb_coupling as dec  # noqa: E402  (row E10's partner cross-correlation and thinnest-axis measures)

TEST_BOARD = "T2:embryo:test_interp"  # read for its gene panel alone (N1's precedent)
TWIN_BOARD = "T2:embryo:val_interp"  # the validation twin whose spec the pseudo-board reuses (never a test-board limit)
DAYS = {"E7.5": 7.5, "E8.0": 8.0}
TARGET_DAY = 7.75
T_PRIME = 0.5  # the world's mixing fraction, the test's
SEEDS = tuple(range(10))
TRUTH_SEED_OFFSET = 2_000_000  # the worlds' couplings and draws, apart from the builds' (BUILD_SEED_OFFSET + s)
FRAC = SUBSAMPLE_FRAC  # the scorer's 10% working copies (a module value so a small synthetic test can use 1.0)
WORLDS = {"fgw_placed": "fgw", "aligned_placed": "aligned"}  # the coupling that places a world's mapped E8.0 half
CANDIDATES = {"fgw_shape0": "fgw", "fgw_shape0_aligned": "aligned"}  # the two s2 candidates, as the harness builds them
ENTRY, CAND = "fgw_shape0", "fgw_shape0_aligned"
OUT = HERE / "results" / "e12d_halves"
NAMES = ("e12d_halves.json", "e12d_halves.txt")
SKILL_READ = ("d2_shape", "occupancy_dice", "scale_log_ratio", "neighborhood_mmd", "de_score", "de_direction", "mmd_u",
              "variogram")
RAW_READ = ("d2_shape", "occupancy_dice", "scale_log_ratio", "neighborhood_mmd")  # raw values, every seed
MIN_USABLE_SEEDS = 8  # of the 10: a world with fewer seeds with sound anchors gives no reading (the adoption rule's 8)


def section_halves(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Boolean masks of a stage's two halves by alternating snapped z sections (even levels A, odd levels B)."""
    index = snap_z(np.asarray(coords, dtype=float)[:, 2]).index
    even = index % 2 == 0
    return even, ~even


def split(a: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    a_mask, b_mask = section_halves(np.asarray(a.obsm[COORD_KEY], dtype=float))
    if not a_mask.any() or not b_mask.any():
        sys.exit("a stage has too few z levels for two halves")
    return a[a_mask].copy(), a[b_mask].copy()


def task_view() -> harness.TaskView:
    """The test geometry as the harness shows it to a builder (build_test.py's _view, for the halves)."""
    return harness.TaskView(key="E12d:halves", task="T2", kind="interp", setting="embryo", test_board=TEST_BOARD,
                            roles=("before", "after"), role_days={"before": DAYS["E7.5"], "after": DAYS["E8.0"]},
                            role_stems={"before": "E7.5", "after": "E8.0"}, reference_role="before",
                            target_day=TARGET_DAY)


def mixture(res: t2_build.Interpolation, Eb: np.ndarray, Ea: np.ndarray) -> np.ndarray:
    """A mixture of real cells, as harness._fgw builds it: the drawn rows' expression, before and after."""
    X = np.empty((len(res.rows), Eb.shape[1]), dtype=np.float32)
    X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    return X


def predict(ctx: harness.BuildContext, move_coupling: str) -> tuple[Stage, t2_build.Interpolation]:
    """One s2 candidate as the harness builds it (harness._fgw at t_pos = 0, its guards and checks included), keeping
    the coupling's internals for the mechanism read; pinned equal to the harness's builders by its test."""
    t, g = ctx.task, list(ctx.board.genes)
    before, after = ctx.inputs["before"], ctx.inputs["after"]
    frac = (t.target_day - t.role_days["before"]) / (t.role_days["after"] - t.role_days["before"])
    Eb, Ea = dense(before[:, g].X).astype(np.float32), dense(after[:, g].X).astype(np.float32)
    res = t2_build.interpolate_clouds(np.asarray(before.obsm[COORD_KEY], float)[:, :3],
                                      np.asarray(after.obsm[COORD_KEY], float)[:, :3], Eb, Ea, frac,
                                      harness.count_for(ctx.board, before.n_obs), ctx.seed,
                                      target_rms=harness.target_rms(ctx), t_pos=0.0, move_coupling=move_coupling,
                                      with_coupling=True)
    stage = harness.finish(ctx, before, mixture(res, Eb, Ea), res.coords, guard_source=res.pool_coords,
                           guard_weights=res.pool_weights)
    return stage, res


def world_truth(before: ad.AnnData, after: ad.AnnData, genes: list[str], n: int, seed: int,
                move_coupling: str) -> tuple[Stage, t2_build.Interpolation]:
    """A pseudo-truth world from half A at t' = 1/2 in E7.5's geometry (t_pos = 0): A's E7.5 cells where they are, A's
    E8.0 cells placed by `move_coupling`, at the size rule on A's own two stages; n cells with no row repeated (the
    caller sets n <= 2 x the smaller stage's cells)."""
    Cb, Ca = (np.asarray(a.obsm[COORD_KEY], float)[:, :3] for a in (before, after))
    Eb, Ea = dense(before[:, genes].X).astype(np.float32), dense(after[:, genes].X).astype(np.float32)
    rms = t2_build.pchip_rms([DAYS["E7.5"], DAYS["E8.0"]], [t2_build.rms_radius(Cb), t2_build.rms_radius(Ca)],
                             TARGET_DAY)
    res = t2_build.interpolate_clouds(Cb, Ca, Eb, Ea, T_PRIME, n, seed, target_rms=rms, t_pos=0.0,
                                      move_coupling=move_coupling, with_coupling=True)
    return Stage(X=mixture(res, Eb, Ea), coords=res.coords.astype(np.float32)), res


def coupling_read(res: t2_build.Interpolation) -> dict:
    """The mechanism read of one interpolation at t_pos = 0: the smallest whitened cross-correlation of anchors and
    their partners (row E10's measure, below 0 where partners are mirrored) and the pool's thinnest axis over its RMS
    radius (N1's measure of the s2's pool)."""
    c = res.coupling
    ib, ia, Tm = c["anchors_before"], c["anchors_after"], c["moves"]
    Zb, Za = c["before_frame"], c["after_frame"]
    partners = (Tm @ Za[ia]) / Tm.sum(1, keepdims=True)
    return {"xcorr_min": dec._xcorr(Zb[ib], partners)[0],
            "pool_thinnest": dec._shape(res.pool_coords, res.pool_weights)[-1]}


def pseudo_board(boards: dict, genes: list[str], n_truth: int, n_ref: int):
    """The board the worlds are scored on: the validation twin's spec with the test panel and a pseudo-task's counts
    (harness.pseudo_counts), no field taken from a test board."""
    return replace(boards[TWIN_BOARD], key="pseudo:E12d_halves", split="pseudo", anchors={}, genes=list(genes),
                   **harness.pseudo_counts(n_truth, n_ref))


def _raw(metrics: Mapping) -> dict[str, float | None]:
    return {m: (None if metrics.get(m) is None else float(metrics[m])) for m in RAW_READ}


def score_world(board, truth: Stage, ref: Stage, preds: dict[str, Stage], seed: int) -> dict:
    """The floor (B's E7.5, checked to score 50 with complete anchors, as vec.harness.score_seeds) and both candidates
    on one world's seed, against the same anchors; their raw shape metrics are kept in every case. A seed whose computed
    anchors are inverted or missing (the floor cannot then score 50) is returned as "degenerate", with no scores: in a
    world that keeps E7.5's geometry the floor can sit as close to the truth's shape as the ceiling does (the module
    docstring)."""
    floor_pred = subsample(ref, FRAC, seed)
    f = score_board(board, floor_pred, truth, ref, seed=seed, anchors="computed", frac=FRAC, build=BUILD,
                    protocol=PROTOCOL)
    missing = sorted(set(harness.RANKED_METRICS["T2"]) - set(f.anchors))
    inverted = list(f.extras["inverted_anchors"])
    out: dict = {"scores": {}, "skills": {}, "raw": {}}
    if inverted or missing:
        out["degenerate"] = {"inverted": inverted, "missing": missing}
        w_truth, _, w_ref = working_copies(truth, ref, seed, FRAC, PROTOCOL)
        for name, pred in preds.items():
            out["raw"][name] = _raw(panel_metrics("T2", pred, w_truth, w_ref, seed=seed, build=BUILD))
        return out
    if abs(f.score - 50.0) > 1e-9:
        raise ValueError(f"seed {seed}: the floor scores {f.score} with complete anchors")
    for name, pred in preds.items():
        r = score_board(board, pred, truth, ref, seed=seed, anchors=f.anchors, frac=FRAC, build=BUILD,
                        protocol=PROTOCOL)
        out["scores"][name], out["skills"][name] = float(r.score), {k: float(v) for k, v in r.skills.items()}
        out["raw"][name] = _raw(r.metrics)
    return out


def world_summary(rows: list[dict]) -> dict:
    """One world over the seeds: each candidate's mean score, the candidate's paired gain over fgw_shape0 and the seeds
    ahead, and the mean change of each ranked metric's skill, over the usable seeds; the mean raw change of the shape
    metrics over every seed (reported, not judged; below 0 is better for d2_shape, neighborhood_mmd and
    scale_log_ratio, above 0 for occupancy_dice)."""
    usable = [r for r in rows if "degenerate" not in r]
    raw = {m: float(np.mean([r["raw"][CAND][m] - r["raw"][ENTRY][m] for r in rows])) for m in RAW_READ}
    base = {"n_seeds": len(usable), "n_seeds_total": len(rows),
            "degenerate": [r["degenerate"] for r in rows if "degenerate" in r], "raw_change_mean": raw}
    if not usable:
        return {**base, "score_mean": None, "gain_mean": None, "gain_sd": None, "seeds_ahead": 0,
                "skill_change_mean": None}
    s = {n: np.array([r["scores"][n] for r in usable]) for n in (ENTRY, CAND)}
    gain = s[CAND] - s[ENTRY]
    changes = {m: float(np.mean([r["skills"][CAND][m] - r["skills"][ENTRY][m] for r in usable])) for m in SKILL_READ}
    return {**base, "score_mean": {n: float(v.mean()) for n, v in s.items()}, "gain_mean": float(gain.mean()),
            "gain_sd": float(gain.std(ddof=1)) if len(gain) > 1 else 0.0, "seeds_ahead": int((gain > 0).sum()),
            "skill_change_mean": changes}


def reading(summary: dict[str, dict], full: bool) -> str:
    """Row E12d's three-way reading, fixed before the run; it decides nothing by itself."""
    if not full:
        return "partial run: no reading"
    aligned, fgw = summary["aligned_placed"], summary["fgw_placed"]
    if min(aligned["n_seeds"], fgw["n_seeds"]) < MIN_USABLE_SEEDS:
        return "no clear read"
    if aligned["gain_mean"] < -adoption.MAX_LOSS:
        return "against"
    if (aligned["gain_mean"] > adoption.MIN_GAIN and aligned["seeds_ahead"] >= adoption.MIN_AHEAD
            and fgw["gain_mean"] >= -adoption.MAX_LOSS):
        return "supports"
    return "no clear read"


def _span(values) -> str:
    lo, hi = min(values), max(values)
    return f"{lo:.2f}" if f"{lo:.2f}" == f"{hi:.2f}" else f"{lo:.2f} to {hi:.2f}"


def _seed_line(world: str, scored: dict) -> str:
    if "degenerate" in scored:
        d = scored["degenerate"]
        return f"{world} degenerate anchors (inverted {d['inverted']}, missing {d['missing']})"
    return f"{world} gain {scored['gain']:+.3f}"


def report_text(rep: dict) -> str:
    lines = ["E12d (experiments/harness/TRIES.md; a diagnostic, not a try): fgw_shape0_aligned against fgw_shape0 at "
             "the test's mapped fraction (t = 1/2) with its mirrored partners,",
             "on the test pair's own section halves, in a pseudo-truth world keeping E7.5's geometry at t' = 1/2. "
             "Target-free: reads E7.5 and E8.0 only; no cell count is printed or written.", ""]
    for w, s in rep["summary"].items():
        coupling = WORLDS[w]
        lines.append(f"World {w} (the mapped E8.0 half of the pseudo-truth placed by the {coupling} coupling); usable "
                     f"seeds {s['n_seeds']} of {s['n_seeds_total']}" + (f" (degenerate anchors in the others: "
                     f"{s['degenerate']})" if s["degenerate"] else "") + ":")
        if s["score_mean"] is None:
            lines.append("  no usable seed: no scores")
        else:
            lines.append(f"  mean score fgw_shape0 {s['score_mean'][ENTRY]:.3f}, candidate "
                         f"{s['score_mean'][CAND]:.3f}; paired gain {s['gain_mean']:+.3f} (SD {s['gain_sd']:.3f}), "
                         f"ahead in {s['seeds_ahead']} of {s['n_seeds']} seeds")
            lines.append("  mean skill change (candidate minus fgw_shape0): "
                         + ", ".join(f"{m} {v:+.4f}" for m, v in s["skill_change_mean"].items()))
        lines.append("  mean raw change, every seed (candidate minus fgw_shape0; below 0 is better for d2_shape, "
                     "neighborhood_mmd and scale_log_ratio, above 0 for occupancy_dice): "
                     + ", ".join(f"{m} {v:+.5f}" for m, v in s["raw_change_mean"].items()))
    m = rep["mechanism"]
    lines += ["", "Mechanism (reported, not judged): smallest partner cross-correlation (below 0: mirrored) and the "
              "pool's thinnest axis over its RMS radius, per seed:",
              f"  half B, fgw_shape0's coupling (FGW): cross-correlation {_span(m['B']['fgw']['xcorr_min'])}; mirrored "
              f"in {sum(v < 0 for v in m['B']['fgw']['xcorr_min'])} of {rep['n_seeds']} seeds; pool thinnest axis "
              f"{_span(m['B']['fgw']['pool_thinnest'])}",
              f"  half B, the candidate's coupling (aligned): cross-correlation "
              f"{_span(m['B']['aligned']['xcorr_min'])}; pool thinnest axis "
              f"{_span(m['B']['aligned']['pool_thinnest'])}",
              f"  half A, the FGW-placed world's coupling: cross-correlation "
              f"{_span(m['A']['fgw_placed']['xcorr_min'])}; "
              f"mirrored in {sum(v < 0 for v in m['A']['fgw_placed']['xcorr_min'])} of {rep['n_seeds']} seeds",
              f"  half A, the aligned-placed world's coupling: cross-correlation "
              f"{_span(m['A']['aligned_placed']['xcorr_min'])}",
              "", f"Reading (row E12d; it decides nothing by itself): {rep['reading']}",
              ("Reading rule: supports if the candidate's paired mean gain is above +0.3 with at least 8 seeds ahead "
               "in the aligned-placed world and not below -0.3 in the FGW-placed world; against if it is below -0.3 "
               "in the aligned-placed world; otherwise no clear read (also when a world has fewer than "
               f"{MIN_USABLE_SEEDS} usable seeds)."),
              "", f"Inputs: {rep['inputs']}",
              f"Commit {rep['stamp']['commit']}, tree clean {rep['stamp']['tree_clean']}, UTC {rep['stamp']['utc']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--data", default=str(harness.DATA_DIR))
    ap.add_argument("--seeds", type=int, nargs="*", default=list(SEEDS))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} exists and is not empty: E12d has already written its outputs; a rerun needs a new TRIES row")
    paths = [out / n for n in NAMES]
    exempt = [*paths, *(provenance.stamp_path(p) for p in paths)]
    started = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    boards = load_boards()
    if TEST_BOARD not in boards:
        sys.exit(f"{TEST_BOARD} is not in boards/: run the board intake first")
    genes = list(boards[TEST_BOARD].genes)
    files = {s: Path(a.data) / "T2" / f"{s}.h5ad" for s in DAYS}
    stages = {s: load_stage(files[s], genes=genes, drop_duplicates=True) for s in DAYS}
    (a_b, b_b), (a_a, b_a) = split(stages["E7.5"]), split(stages["E8.0"])
    n_truth = 2 * min(a_b.n_obs, a_a.n_obs)  # the world's draw: no row repeated, half from each bracket
    board = pseudo_board(boards, genes, n_truth, b_b.n_obs)
    ref = Stage.from_anndata(b_b, genes, need_coords=True)
    view = task_view()
    per_seed, summary_rows = [], {w: [] for w in WORLDS}
    mech = {"A": {w: {"xcorr_min": [], "pool_thinnest": []} for w in WORLDS},
            "B": {c: {"xcorr_min": [], "pool_thinnest": []} for c in ("fgw", "aligned")}}
    for s in a.seeds:
        t0 = time.time()
        build_seed, truth_seed = BUILD_SEED_OFFSET + s, TRUTH_SEED_OFFSET + s
        ctx = harness.BuildContext(task=view, board=harness._candidate_board(board),
                                   inputs={"before": b_b, "after": b_a}, seed=build_seed)
        preds, row = {}, {"seed": s, "build_seed": build_seed, "truth_seed": truth_seed, "worlds": {}}
        for name, coupling in CANDIDATES.items():
            preds[name], res = predict(ctx, coupling)
            for k, v in coupling_read(res).items():
                mech["B"][coupling][k].append(v)
        for world, coupling in WORLDS.items():
            truth, res_t = world_truth(a_b, a_a, genes, n_truth, truth_seed, coupling)
            for k, v in coupling_read(res_t).items():
                mech["A"][world][k].append(v)
            scored = score_world(board, truth, ref, preds, s)
            if "degenerate" not in scored:
                scored["gain"] = scored["scores"][CAND] - scored["scores"][ENTRY]
            row["worlds"][world] = scored
            summary_rows[world].append(scored)
        row["seconds"] = round(time.time() - t0, 1)
        per_seed.append(row)
        print(f"seed {s}: " + "; ".join(_seed_line(w, row["worlds"][w]) for w in WORLDS)
              + f" ({row['seconds']} s)", flush=True)
    summary = {w: world_summary(rows) for w, rows in summary_rows.items()}
    full = sorted(a.seeds) == list(SEEDS)
    inputs = {str(p.relative_to(a.data) if p.is_relative_to(a.data) else p): harness._sha256(p) for p in files.values()}
    inputs["panel (test board's gene list)"] = hashlib.sha256("\n".join(genes).encode()).hexdigest()
    ended = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    clean = (None if started[1] is None or ended[1] is None
             else not started[1] and not ended[1] and started[0] == ended[0])
    stamp = {"commit": ended[0], "tree_clean": clean, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    rep = {"what": "E12d (TRIES.md; a diagnostic, not a try): fgw_shape0_aligned against fgw_shape0 at the test's "
                   "mapped fraction on the test pair's own section halves",
           "test_board": TEST_BOARD, "panel_board": TEST_BOARD, "scoring_board": "validation twin spec, pseudo counts",
           "worlds": WORLDS, "candidates": CANDIDATES, "t_prime": T_PRIME, "seeds": list(a.seeds),
           "n_seeds": len(a.seeds),
           "summary": summary, "reading": reading(summary, full), "mechanism": mech, "per_seed": per_seed,
           "inputs": inputs, "stamp": stamp}
    text = report_text(rep)
    out.mkdir(parents=True, exist_ok=True)
    (out / NAMES[0]).write_text(json.dumps(rep, indent=1, default=float) + "\n")
    (out / NAMES[1]).write_text(text)
    for p in paths:
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
