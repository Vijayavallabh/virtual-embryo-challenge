#!/usr/bin/env python
"""Try H28 of experiments/harness/TRIES.md: the heart's trend step (b0_trend, row H11) at strengths 1, 0.5 and 0.25 on
the full trim the captain named for the heart's last file, against that file's own copy expression, with a placebo of
the same factor sizes per strength. A labelled diagnostic on one truth (E10.5); it changes no file, rule or default.

The control is row H27's arm C at its forecast size f (`core_iso_half:f`): `vec.harness._copy` on a context whose
reference input is E10.5's pool under `radial_trim_mask(C, core, 0.0)` (earlier stages whole, as the heart builder loads
them), then `t2_build.apply_size_procedure` in raw mode with the core hypothesis on the in-memory size result row H27
forms from the input stages only, then row H25's half, then the checked builder (`h27.size_arm`, imported). Every arm
here is the control with ONE change, the expression step, made before the size step:
  trend_s{s}    b0_trend's expression (`vec.harness._b0_t2(ctx, trend=True)`, rebuilt step by step below) with every
                lineage's log factor multiplied by s: the logf_by_type mapping of `b0_trend_factors` times s, and the
                damping that `momentum_shift` gives a lineage the mapping leaves out (B0's own factor) times s. So s = 1
                is the harness's own b0_trend expression on the same cells (guard V5 checks it bit for bit in every
                seed), and s = 0 leaves every count as drawn (the reference trend_s0, not counted). At the test the
                step's time budget is B0_DAMPING x the last gap (E9.5 to E10.5) = 0.5 day; on the backtests it is
                0.375 day (E9.5 -> E10.5) and 0.125 day (E8.75 -> E10.5). Strength s scales all three alike.
  placebo_s{s}  the same factors at strength s with no developmental content: within each lineage, the genes' factors
                moved by one fixed derangement inside ten bins of the lineage's mean log level at the reference
                (PERM_SEED; a gene never keeps its own factor in a bin of two or more). The same magnitudes, the same
                count-space thinning and top-up, on other genes of like level.
The draw (`draw_copy`), coordinates, size step, half rows and cell number are the control's in every arm (guard V5).

Transitions: the two harness backtests on E10.5 (`T2h/<=E9.5->E10.5`, `T2h/<=E8.75->E10.5`; one truth, merged seed by
seed into fold E10.5). Row H27's third transition, E8.25_late -> E8.5, has one input stage: b0_trend is not defined
there (`trend_stages` raises below three), so it is not built and is printed as not buildable. Truths: the truth as
released (own) and row H26's CLEAN world (radial_0, thin_0.5, thin_0, equal weights), at payoff.TRUTH_SEED.

  run   one transition, a seed range, one JSON under --out. Refused before anything loads unless row H28 is in
        TRIES.md, the tree is clean by vec.provenance's rule and the output file does not exist. Prints no score.
  read  refused without the row or with an existing reading; applies the void conditions, the gate and the fixed
        statements before it prints anything; prints only the outcome words and what follows; writes the numbers to
        reading.json and reading.txt.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
import zlib
from dataclasses import replace
from pathlib import Path

import numpy as np

from vec import harness, prep, provenance, t1_scope, t2_build
from vec.adoption import SEEDS
from vec.boards import load_boards
from vec.core import counts_of, mean_cp10k
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parents[2]  # this checkout


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h27 = _load("h27_heart_forecast_size", HERE / "experiments" / "harness" / "h27_heart_forecast_size.py")
h26, h25, payoff = h27.h26, h27.h25, h27.payoff

ROW = "H28"
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h28-heart-trend-strength"
PART = "h28.json"
SOURCE = "experiments/harness/h28_heart_trend_strength.py (try H28)"
TASK_KEYS = h26.HARNESS_KEYS  # the two backtests on E10.5
NOT_BUILT = "T2h/E8.25_late->E8.5"
NOT_BUILT_LINE = ("fold E8.5 (T2h/E8.25_late->E8.5): not buildable: one input stage (b0_trend needs at least three); "
                  "no arm is built there")
OWN, CLEAN = h26.OWN, h26.CLEAN
TRUTHS = (OWN, *CLEAN)
WORLDS = (("own", (OWN,)), ("CLEAN", CLEAN))
FORM = "core_iso_half"  # C, the full trim
SIZING = "f"  # row H27's forecast size
H27_ARM = h27.arm_name(FORM, SIZING)
STRENGTHS = (1.0, 0.5, 0.25)  # counted, with a placebo each
ZERO_STRENGTH = 0.0
CONTROL = "control"
ZERO = "trend_s0"


def arm(s: float) -> str:
    return f"trend_s{s:g}"


def placebo(s: float) -> str:
    return f"placebo_s{s:g}"


COUNTED = tuple(n for s in STRENGTHS for n in (arm(s), placebo(s)))
BUILT = (CONTROL, ZERO, *COUNTED)
ALL = (harness.FLOOR, *BUILT)
PERM_SEED = 20261120  # the placebo's derangement: one per lineage, from this seed and the lineage's name
PERM_BINS = 10
TEST_BUDGET = harness.B0_DAMPING * (10.5 - 9.5)  # b0_trend's budget at the test: half the last gap E9.5 -> E10.5
PRACTICE_BUDGET = {"T2h/<=E9.5->E10.5": 0.375, "T2h/<=E8.75->E10.5": 0.125}  # V6 checks the records against these
GUARD_MEAN, GUARD_SEED = h25.GUARD_MEAN, h25.GUARD_SEED  # 0.1 and 0.5 board points (row H27's gate)
TOL = 1e-9
FACTOR_RTOL = 1e-12
DE_TERMS = ("de_score", "de_direction")
TERMS = ("de_score", "de_direction", "mmd_u", "variogram", "neighborhood_mmd")
SHAPE_TERMS = ("d2_shape", "occupancy_dice", "scale_log_ratio")  # must not move: the cells' positions are the control's
BAR = 0.3  # board points, paired mean, on fold E10.5
MIN_SEEDS = 8  # of 10 seeds above zero
TAG = "[practice diagnostic, released stages, one truth]"


# ------------------------------------------------------------------------------------------------ refusals
def tries_rows(path: Path = TRIES) -> set[str]:
    return h27.tries_rows(path)


def refuse(path: Path) -> None:
    """The refusals before anything loads: no row H28, an unclean tree (vec.provenance's rule), an existing output."""
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    _, unclean = provenance.git_state(HERE)
    if unclean is None or unclean:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    if Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


# ------------------------------------------------------------------------------------------------ the expression step
def _typed(a, tag: str):
    """`_b0_t2`'s typing: obs["merged"], with Untracked tagged by stage so it never pairs across stages."""
    m = np.asarray(a.obs["merged"]).astype(str)
    return harness._shallow(a, b0_type=np.where(m == harness.T2_UNTRACKED, f"{harness.T2_UNTRACKED}@{tag}", m))


def lineage_log_means(ctx: harness.BuildContext, role: str, genes, lineages) -> dict[str, np.ndarray]:
    """y = log(mean CP10k + eps) of each lineage's own cells at one input stage (b0_trend_factors' own y)."""
    eps = float(harness._momentum_default("eps"))
    a = ctx.inputs[role]
    C = counts_of(a, genes)
    merged = np.asarray(a.obs[harness.T2_LINEAGE_KEY]).astype(str)
    return {lin: np.log(mean_cp10k(C[np.flatnonzero(merged == lin)]) + eps) for lin in lineages}


def derangement(basis: np.ndarray, key: str, bins: int = PERM_BINS, seed: int = PERM_SEED) -> np.ndarray:
    """A permutation of the genes that moves every gene to another gene of the same level bin (ten equal bins of
    `basis` by rank); a bin of one gene keeps it. Fixed by `seed` and the lineage's name `key`."""
    basis = np.asarray(basis, dtype=np.float64)
    rng = np.random.default_rng([seed, zlib.crc32(key.encode())])
    perm = np.arange(len(basis))
    for b in np.array_split(np.argsort(basis, kind="stable"), bins):
        if len(b) < 2:
            continue
        o = rng.permutation(b)
        perm[o] = np.roll(o, -1)
    return perm


def trend_expression(ctx: harness.BuildContext, strength: float, placebo_arm: bool = False):
    """b0_trend's expression step at `strength` on the context's draw: `_b0_t2(ctx, trend=True)` line by line, with the
    trend mapping times strength and momentum_shift's damping times strength (B0's own factor for a lineage the mapping
    leaves out). placebo_arm: every shifted lineage's factor (the mapping's, or B0's own for a lineage it leaves out),
    at strength, with its genes moved by `derangement` on the lineage's reference level. Returns the stage through
    `harness.finish`, the shifted counts, the draw's rows and the target-free record."""
    harness._b0_trend_task_check(ctx.task)
    s = float(strength)
    pool, rows, _, coords = harness.draw_copy(ctx)
    start, end = ctx.inputs["start"], pool
    like = _typed(pool, "end")[rows]
    genes = list(like.var_names)
    types = sorted(set(np.asarray(like.obs["b0_type"]).astype(str)))
    logf, rec = harness.b0_trend_factors(ctx, genes, types)
    fallback = sorted(r["lineage"] for r in rec["lineages"] if r["own_factor"])
    mapping = {k: s * v for k, v in logf.items()}
    damping = harness.B0_DAMPING * s
    fixed_points = {}
    if placebo_arm:
        shifted_lins = sorted(r["lineage"] for r in rec["lineages"])
        y_ref = lineage_log_means(ctx, ctx.task.reference_role, genes, shifted_lins)
        y_start = lineage_log_means(ctx, "start", genes, fallback) if fallback else {}
        full = {**{k: s * (harness.B0_DAMPING * (y_ref[k] - y_start[k])) for k in fallback}, **mapping}
        out = {}
        for k in shifted_lins:
            p = derangement(y_ref[k], k)
            out[k] = np.asarray(full[k], dtype=np.float64)[p]
            fixed_points[k] = int((p == np.arange(len(p))).sum())
        mapping = out
    shifted, report = t1_scope.momentum_shift(_typed(start, "start"), _typed(end, "end"), like, damping, (),
                                              key="b0_type", seed=[ctx.seed, 3], logf_by_type=mapping)
    st = harness.finish(ctx, pool, dense(shifted[:, ctx.board.genes].X), coords)
    record = {"strength": s, "placebo": bool(placebo_arm), "budget_days": float(rec["budget_days"]),
              "effective_budget_days": s * float(rec["budget_days"]), "stages": rec["stages"], "days": rec["days"],
              "n_trend": int(rec["n_trend"]), "fallback": fallback,
              "lineages": [{k: r[k] for k in ("lineage", "stages", "n_cells", "own_factor", "cosine", "norm_ratio")}
                           for r in rec["lineages"]],
              "shifted": sorted(r["type"] for r in report.rows if r["status"] == "shifted"),
              "fixed_points": fixed_points, "ref_cells": int(ctx.inputs[ctx.task.reference_role].n_obs)}
    return st, counts_of(shifted, list(ctx.board.genes)), rows, record


def whole_pool_comparison(ctx: harness.BuildContext, whole_ref, genes, types) -> dict:
    """Target-free: per lineage, the trimmed-pool trend factor against the whole-pool one (cosine, norm ratio)."""
    trimmed, _ = harness.b0_trend_factors(ctx, genes, types)
    whole, _ = harness.b0_trend_factors(replace(ctx, inputs={**ctx.inputs, "ref": whole_ref}), genes, types)
    out = {}
    for k in sorted(set(trimmed) | set(whole)):
        if k in trimmed and k in whole:
            n = float(np.linalg.norm(whole[k]))
            out[k] = {"cosine": harness._cosine(trimmed[k], whole[k]),
                      "norm_ratio": float(np.linalg.norm(trimmed[k])) / n if n > 0 else float("nan")}
        else:
            out[k] = {"in_trimmed": k in trimmed, "in_whole": k in whole}
    return out


def build_arms(ctx: harness.BuildContext, size: dict, whole_ref=None) -> tuple[dict, dict]:
    """Every arm of one seed on the full trim's context: the control (`_copy`), trend_s0, the trend and placebo arms,
    each then sized by row H27's f step, halved and checked; and the identity checks of guard V5."""
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    pre = {CONTROL: (harness._copy(ctx), None, None, {"strength": None})}
    pre[ZERO] = trend_expression(ctx, ZERO_STRENGTH)
    for s in STRENGTHS:
        pre[arm(s)] = trend_expression(ctx, s)
        pre[placebo(s)] = trend_expression(ctx, s, placebo_arm=True)
    ref_b0 = harness._b0_t2(ctx, trend=True)  # the harness's own b0_trend on the same context (guard V5, not scored)
    pool, rows, _, _ = harness.draw_copy(ctx)
    drawn = counts_of(pool[rows], list(ctx.board.genes))
    ctrl = pre[CONTROL][0]
    checks = {"s1_equals_harness_b0_trend": bool(np.array_equal(pre[arm(1.0)][0].X, ref_b0.X)
                                                 and np.array_equal(pre[arm(1.0)][0].coords, ref_b0.coords)),
              "s0_counts_equal_drawn": bool((pre[ZERO][1] != drawn).nnz == 0),
              "s0_x_max_abs_diff_control": float(np.abs(np.asarray(pre[ZERO][0].X, np.float64)
                                                        - np.asarray(ctrl.X, np.float64)).max()),
              "same_draw": bool(all(np.array_equal(v[2], rows) for k, v in pre.items() if k != CONTROL)),
              "pre_size_coords_equal": bool(all(np.array_equal(v[0].coords, ctrl.coords) for v in pre.values()))}
    out, recs = {}, {}
    for name, (st, _, _, rec) in pre.items():
        X, C, srec = h27.size_arm(st.X, st.coords, FORM, SIZING, {}, size, {}, seed)
        sub, rep = t2_build.build_submission(ctx.board, X, C, genes=list(ctx.board.genes))
        if not rep.ok:
            raise ValueError(f"{ctx.task.key} {name}: the built file fails its check: {rep.errors}")
        out[name] = Stage.from_anndata(sub, ctx.board.genes, need_coords=True)
        recs[name] = {"size": srec, "trend": rec}
    checks["post_size_coords_equal"] = bool(all(np.array_equal(v.coords, out[CONTROL].coords) for v in out.values()))
    checks["n_cells"] = {k: int(v.n) for k, v in out.items()}
    if whole_ref is not None:
        like = _typed(pool, "end")[rows]
        checks["whole_pool"] = whole_pool_comparison(ctx, whole_ref, list(like.var_names),
                                                     sorted(set(np.asarray(like.obs["b0_type"]).astype(str))))
    return out, {"records": recs, "checks": checks}


# ------------------------------------------------------------------------------------------------ the run
def run_transition(key: str, seeds, out: Path, data_dir: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse(path)
    if key not in TASK_KEYS:
        sys.exit(f"{key}: not buildable (one input stage) or not a task of row {ROW}")
    t = payoff.heart_task(key)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    geoms, _ = h27.input_geometries(t, data_dir)  # the forecast first, from the input stages only (row H27's)
    size = h27.analogue_size_result(geoms, float(harness.STAGE_DAYS[t.target]))
    ftargets = h27.forecast_targets(size)
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    tm = h26.truth_masks(TC, tcore)
    truths = {q: Stage.from_anndata(ta[tm[q]], genes, need_coords=True) for q in TRUTHS}
    t_info = {q: h26.truth_info(TC, tcore, tm[q]) for q in TRUTHS}
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    whole = inputs["ref"]
    pcore = np.asarray(whole.obs["core"], dtype=bool)
    PC = np.asarray(whole.obsm[COORD_KEY], dtype=float)[:, :3]
    keep = h26.pool_masks(PC, pcore)[FORM]
    pool = whole[keep].copy()
    p_info = {FORM: h26.pool_info(PC, pcore, keep)}
    v7 = h27.v7_crosscheck(PC, size)
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        ctx = harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pool},
                                   seed=harness.BUILD_SEED_OFFSET + s)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            preds, extra = build_arms(ctx, size, whole_ref=whole)
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": extra["checks"]["n_cells"],
               "checks": extra["checks"], "records": extra["records"], "scores": {}, "skills": {}, "groups": {}}
        for q, truth in truths.items():
            board = replace(base_board, truth_cells=truth.n)
            fl = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, s), truth, ref_stage, seed=s,
                             anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
            row["scores"][q] = {harness.FLOOR: fl.score}
            row["skills"][q] = {harness.FLOOR: fl.skills}
            row["groups"][q] = {harness.FLOOR: dict(fl.group_skill)}
            for name, pred in preds.items():
                r = score_board(board, pred, truth, ref_stage, seed=s, anchors=fl.anchors, build=harness.BUILD,
                                protocol=harness.PROTOCOL)
                row["scores"][q][name], row["skills"][q][name] = r.score, r.skills
                row["groups"][q][name] = dict(r.group_skill)
        row["seconds"] = round(time.time() - t0, 1)
        rows_out.append(row)
        print(f"{key} seed {s} ({row['seconds']} s) done", flush=True)
    res = {"label": "try H28: the trend step at three strengths and its placebos on the full trim; a diagnostic on "
                    "one truth, not an adoption result", "axis": AXIS, "task": t.key, "arms": list(BUILT),
           "truth_values": list(TRUTHS), "forecast_targets": ftargets, "size_result": size, "v7": v7,
           "truth_info": t_info, "pred_info": p_info, "truth_seed": payoff.TRUTH_SEED,
           "half_seed_base": h25.HALF_SEED_BASE, "perm_seed": PERM_SEED, "seeds": list(seeds), "rows": rows_out,
           "source": SOURCE, "provenance": payoff.provenance()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


# ------------------------------------------------------------------------------------------------ loading
def load_parts(root: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for p in sorted(Path(root).rglob(PART)):
        r = json.loads(p.read_text())
        if r.get("axis") == AXIS:
            out.setdefault(r["task"], []).append(r)
    return out


H27_PARTS = HERE / "experiments" / "harness" / "results" / "h27_heart_forecast_size" / "parts"


def load_h27(root: Path = H27_PARTS) -> dict[str, list[dict]]:
    return h27.load_parts(root)


# ------------------------------------------------------------------------------------------------ the guards
def gate(key: str, rows: list[dict], h27_rows: list[dict], mine: str = CONTROL, theirs: str = H27_ARM,
         truths=TRUTHS) -> tuple[list[str], dict]:
    """V4, the gate: `mine` against row H27's committed `theirs`, per truth: the ten-seed mean within GUARD_MEAN and
    every seed within GUARD_SEED."""
    other = {r["seed"]: r for r in h27_rows}
    why, big = [], {"mean": 0.0, "seed": 0.0}
    for q in truths:
        try:
            d = np.array([r["scores"][q][mine] - other[r["seed"]]["scores"][q][theirs] for r in rows], float)
        except KeyError:
            why.append(f"V4: {key} {q}: no committed row H27 score to compare")
            continue
        big["mean"], big["seed"] = max(big["mean"], abs(float(d.mean()))), max(big["seed"], float(np.abs(d).max()))
        if abs(d.mean()) > GUARD_MEAN:
            why.append(f"V4: {key} {q}: the control's ten-seed mean is {d.mean():+.3f} from row H27's")
        if np.abs(d).max() > GUARD_SEED:
            why.append(f"V4: {key} {q}: a seed is {np.abs(d).max():.3f} from row H27's")
    return why, big


def identities(key: str, rows: list[dict]) -> list[str]:
    """V5: strength 1 is the harness's b0_trend on the same context, bit for bit; strength 0 leaves every drawn count
    as drawn; every arm has the control's draw, coordinates (before and after the size step and half) and cell
    number; and no shape term moves."""
    why = []
    for r in rows:
        c = r["checks"]
        tag = f"V5: {key} seed {r['seed']}"
        if c.get("s1_equals_harness_b0_trend") is not True:
            why.append(f"{tag}: strength 1 is not the harness's b0_trend expression on the same cells")
        if c.get("s0_counts_equal_drawn") is not True:
            why.append(f"{tag}: strength 0 changed a drawn count")
        for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
            if c.get(k) is not True:
                why.append(f"{tag}: {k.replace('_', ' ')} fails")
        if len(set(r["n_cells"].values())) != 1 or set(r["n_cells"]) != set(BUILT):
            why.append(f"{tag}: the arms' cell numbers differ")
        for q in TRUTHS:
            for name in BUILT:
                for m in SHAPE_TERMS:
                    if abs(r["skills"][q][name][m] - r["skills"][q][CONTROL][m]) > TOL:
                        why.append(f"{tag} {q}: {name}'s {m} is not the control's")
    return why


def trend_records(key: str, rows: list[dict], pool_cells: int) -> list[str]:
    """V6: each trend arm's record has the transition's own budget, the trimmed pool as its reference, the same
    shifted lineages at every strength, and every placebo a derangement (no gene keeps its factor)."""
    why = []
    want = PRACTICE_BUDGET[key]
    for r in rows:
        shifted = set()
        for name in BUILT:
            rec = r["records"][name]["trend"]
            if name == CONTROL:
                continue
            if abs(rec["budget_days"] - want) > FACTOR_RTOL * want:
                why.append(f"V6: {key} seed {r['seed']}: {name}'s budget is {rec['budget_days']}, not {want}")
            if rec["ref_cells"] != pool_cells:
                why.append(f"V6: {key} seed {r['seed']}: {name}'s reference is not the trimmed pool")
            if name != ZERO:
                shifted.add(tuple(rec["shifted"]))
            if rec["placebo"] and any(v for v in rec["fixed_points"].values()):
                why.append(f"V6: {key} seed {r['seed']}: {name} leaves a gene on its own factor")
        if len(shifted) != 1:
            why.append(f"V6: {key} seed {r['seed']}: the arms shift different lineages")
    return why


def void_reasons(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], commit: str | None = None
                 ) -> tuple[list[str], dict]:
    why: list[str] = []
    dist: dict = {"V4": {}}
    every = [p for ps in parts.values() for p in ps]
    complete = True
    for key in TASK_KEYS:  # V1
        ps = parts.get(key)
        if not ps:
            why.append(f"V1: missing transition {key}")
            complete = False
            continue
        rows, twice = h26.merge_rows(ps)
        if twice:
            why.append(f"V1: {key}: seeds {twice} occur in more than one part")
            complete = False
        if [r["seed"] for r in rows] != list(SEEDS):
            why.append(f"V1: {key}: seeds {[r['seed'] for r in rows]}, not {list(SEEDS)}")
            complete = False
        for r in rows:
            for field in ("scores", "skills", "groups"):
                gone = sorted({f"{q}/{a}" for q in TRUTHS for a in ALL if a not in r.get(field, {}).get(q, {})})
                if gone:
                    why.append(f"V1: {key} seed {r['seed']}: missing from {field}: {gone[:6]}")
                    complete = False
            if set(r.get("records", {})) != set(BUILT) or "checks" not in r:
                why.append(f"V1: {key} seed {r['seed']}: records or checks missing")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"V1: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not every:
        return why or ["V1: no part"], dist
    if len({p["provenance"]["commit"] for p in every}) > 1:  # V2
        why.append("V2: the parts come from different commits")
    if not all(p["provenance"]["tree_clean"] for p in every):
        why.append("V2: a part comes from a dirty tree")
    if not commit or any(not str(p["provenance"]["commit"]).startswith(commit) for p in every):
        why.append(f"V2: a part does not come from the row's commit {commit!r}")
    if not complete:
        return why, dist
    for key in TASK_KEYS:
        ps, theirs = parts[key], h27_parts.get(key)
        if not theirs:
            why.append(f"V3: {key}: no committed row H27 part")
            continue
        ref = theirs[0]
        for p in ps:  # V3
            for q in TRUTHS:
                if not h27._same(p["truth_info"].get(q), ref["truth_info"].get(q), "cells"):
                    why.append(f"V3: {key}: truth {q} differs from row H27's (cells or core share)")
            if not h27._same(p["pred_info"].get(FORM), ref["pred_info"].get(FORM), "pool_cells"):
                why.append(f"V3: {key}: the full trim's pool differs from row H27's (cells or core share)")
            f, g = p["forecast_targets"][FORM], ref["forecast_targets"][FORM]
            if not abs(f - g) <= FACTOR_RTOL * abs(g):
                why.append(f"V3: {key}: the forecast size differs from row H27's")
            if not p["v7"]["ok"]:
                why.append(f"V3: {key}: the reference cloud is not the cloud the size result measured")
        rows = h26.merge_rows(ps)[0]
        w4, dist["V4"][key] = gate(key, rows, h26.merge_rows(theirs)[0])
        why += w4 + identities(key, rows) + trend_records(key, rows, ps[0]["pred_info"][FORM]["pool_cells"])
    return why, dist


# ------------------------------------------------------------------------------------------------ the tables
def _stats(d) -> dict:
    d = np.asarray(d, float)
    return {"mean": float(d.mean()), "lo": float(d.min()), "hi": float(d.max()), "ahead": int((d > 0).sum()),
            "se": float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else float("nan")}


def series(rows: dict[str, list[dict]], keys, truths, name: str, what: str = "score") -> np.ndarray:
    """Per seed: `name`'s board score (or one skill) averaged over `truths`, then over the transitions `keys` seed by
    seed (the fold)."""
    def one(key, q):
        if what == "score":
            return np.array([r["scores"][q][name] for r in rows[key]], float)
        return np.array([r["skills"][q][name][what] for r in rows[key]], float)
    return np.mean([np.mean([one(key, q) for q in truths], axis=0) for key in keys], axis=0)


def compare(rows, keys, truths, a: str, b: str) -> dict:
    """a minus b, paired by seed: the board score and every term."""
    out = {"score": _stats(series(rows, keys, truths, a) - series(rows, keys, truths, b))}
    for m in (*TERMS, *SHAPE_TERMS):
        out[m] = float((series(rows, keys, truths, a, m) - series(rows, keys, truths, b, m)).mean())
    return out


PAIRS = (*((arm(s), CONTROL) for s in STRENGTHS), (ZERO, CONTROL), *((arm(s), placebo(s)) for s in STRENGTHS),
         *((placebo(s), CONTROL) for s in STRENGTHS), (arm(0.5), arm(1.0)), (arm(0.25), arm(1.0)),
         (arm(0.25), arm(0.5)))
UNITS = {**{k: (k,) for k in TASK_KEYS}, "fold E10.5": TASK_KEYS}


def pname(a: str, b: str) -> str:
    return f"{a} - {b}"


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: {pname(a, b): compare(rows, ks, ts, a, b) for a, b in PAIRS} for u, ks in UNITS.items()}
            for w, ts in WORLDS}


# ------------------------------------------------------------------------------------------------ the statements
def passes(s: dict) -> bool:
    """The bar: paired mean above +BAR board points with at least MIN_SEEDS of 10 seeds above zero."""
    return bool(s["mean"] > BAR and s["ahead"] >= MIN_SEEDS)


def statements(tab: dict) -> dict:
    """The fixed statements on fold E10.5, in both worlds (own and CLEAN), per strength."""
    def fold(w, a, b):
        return tab[w]["fold E10.5"][pname(a, b)]
    out = {}
    for s in STRENGTHS:
        ahead = all(passes(fold(w, arm(s), CONTROL)["score"]) for w, _ in WORLDS)
        over_placebo = all(passes(fold(w, arm(s), placebo(s))["score"]) for w, _ in WORLDS)
        on_de = all(fold(w, arm(s), CONTROL)[m] > 0 for w, _ in WORLDS for m in DE_TERMS)
        st = {"ahead": ahead, "over_placebo": over_placebo, "on_de": on_de,
              "qualified": bool(ahead and over_placebo and on_de),
              "least_world_mean": min(fold(w, arm(s), CONTROL)["score"]["mean"] for w, _ in WORLDS)}
        out[s] = st
    for s in STRENGTHS:
        out[s]["beats"] = {o: all(passes(fold(w, arm(s), arm(o))["score"]) for w, _ in WORLDS)
                           for o in STRENGTHS if o > s}
    return out


def output(st: dict) -> float | None:
    """The default is the largest qualified strength (1 first). A smaller qualified strength is named only if it
    beats that default itself, paired, by the bar in both worlds; among several, the one whose least-world mean over
    the copy is highest (ties: the larger strength)."""
    qual = [s for s in STRENGTHS if st[s]["qualified"]]
    if not qual:
        return None
    default = qual[0]
    smaller = [s for s in qual if s < default and st[s]["beats"].get(default)]
    if not smaller:
        return default
    return max(smaller, key=lambda s: (st[s]["least_world_mean"], s))


def outcome(st: dict | None, void: bool) -> str:
    """Exactly one of VOID, AHEAD, PLACEBO, OFF_DE, NONE."""
    if void or st is None:
        return "VOID"
    if output(st) is not None:
        return "AHEAD"
    ahead = [s for s in STRENGTHS if st[s]["ahead"]]
    if not ahead:
        return "NONE"
    if any(not st[s]["over_placebo"] for s in ahead):
        return "PLACEBO"
    return "OFF_DE"


# ------------------------------------------------------------------------------------------------ the fixed words
def budget_at_test(s: float) -> float:
    return s * TEST_BUDGET


def j6_line(s: float) -> str:
    ratios = ", ".join(f"{k}: practice {s * b:.4g} day, ratio {budget_at_test(s) / (s * b):.3g}"
                       for k, b in PRACTICE_BUDGET.items())
    return (f"strength {s:g}: test budget {budget_at_test(s):.4g} day; {ratios}; above 1 on both backtests, so J6's "
            f"first branch: the gain is shown as a range from zero to its practice value, never as a point")


FOLLOWS = {
    "VOID": ("VOID: {why}. No figure exists under row H28 and there is no rerun under it; the ruled file stands (the "
             "full trim with the copy expression)."),
    "AHEAD": ("AHEAD: the trend step at strength {s:g} (test budget {b:.4g} day) is ahead of the copy on the full "
              "trim, ahead of its placebo and on the DE terms, on fold E10.5 in both worlds (one truth; a diagnostic). "
              "What follows: ONE priced question goes to the captain: his ruled file (the full trim, expression copy) "
              "beside the full trim with the trend step at strength {s:g} (budget {b:.4g} day at the test), with the "
              "practice figures in both worlds as ranges over the two backtests and never a point; rule 9's J6 (the "
              "test's budget over each backtest's is 1.33 and 4, so the gain is a range from zero to its practice "
              "value); J1's caveat (the test's slope uses four stages and a last step across the change from "
              "sectioned embryo portions to dissected hearts, which no backtest has, compounded here with the trim's "
              "selection); and, as its own item, the reversal of the amendment of 2026-11-18 of docs/S2_RULES.md "
              "(the S2_RULES amendment is owed before any build). His ruled file stands unless he changes it; an "
              "off-by-default builder option and check builds follow only on his word."),
    "PLACEBO": ("MATCHED BY ITS PLACEBO: a strength ahead of the copy does not beat its own placebo by the bar in "
                "both worlds; the gain is count noise or depth, the heart's closed depth lever. What follows: no "
                "question; the captain is told in one line; the ruled file stands."),
    "OFF_DE": ("OFF THE DE TERMS: a strength ahead of the copy and of its placebo does not gain on both de_score and "
               "de_direction in both worlds; the gain sits where the step's reason did not put it (count noise or "
               "depth, a closed lever). What follows: no question; the captain is told in one line; the ruled file "
               "stands."),
    "NONE": ("NO STRENGTH AHEAD: no strength is ahead of the copy on the full trim by the bar in both worlds. What "
             "follows: the ruled file stands; the captain is told the result in one line."),
}
L0 = ("Try H28 (commit {commit}). A diagnostic on one truth (E10.5, two backtests merged seed by seed); every figure "
      "is a practice reading on released stages and none is to be read beside a test score. Nothing is built and no "
      "file, rule or default changes; the heart's scope and size are the captain's ruling and are not in question.")
L1 = ("Arms: the control is row H27's full trim at its forecast size (copy expression); trend_s<s> is b0_trend's "
      "expression with every lineage's log factor times s, before the size step; placebo_s<s> is the same factors "
      "with each lineage's genes deranged within ten level bins; trend_s0 is a reference (strength zero, counts as "
      "drawn, expression re-rendered from counts as the trend path renders it).")


def _st(s: dict) -> str:
    return f"{s['mean']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10 (se {s['se']:.3f})"


def _terms(c: dict) -> str:
    return ", ".join(f"{m} {c[m]:+.4f}" for m in (*TERMS, *SHAPE_TERMS))


def statement_lines(st: dict, tab: dict) -> list[str]:
    def fold(w, a, b):
        return tab[w]["fold E10.5"][pname(a, b)]["score"]
    lines = []
    def de(w, s, m):
        return tab[w]["fold E10.5"][pname(arm(s), CONTROL)][m]
    for s in STRENGTHS:
        x = st[s]
        lines.append(f"S1 strength {s:g} AHEAD of the copy on the full trim: {str(x['ahead']).upper()} (fold E10.5; "
                     + "; ".join(f"{w} {_st(fold(w, arm(s), CONTROL))}" for w, _ in WORLDS)
                     + f"; bar: mean above +{BAR} with at least {MIN_SEEDS} of 10 seeds above 0, in both worlds)")
        lines.append(f"S2 strength {s:g} AHEAD of its placebo: {str(x['over_placebo']).upper()} (fold E10.5; "
                     + "; ".join(f"{w} {_st(fold(w, arm(s), placebo(s)))}" for w, _ in WORLDS) + "; the same bar)")
        lines.append(f"S3 strength {s:g} gains on the DE terms: {str(x['on_de']).upper()} (fold E10.5, arm minus "
                     "copy; " + "; ".join(f"{w} de_score {de(w, s, 'de_score'):+.4f}, de_direction "
                                          f"{de(w, s, 'de_direction'):+.4f}" for w, _ in WORLDS)
                     + "; both above 0 in both worlds)")
        for o, ok in x["beats"].items():
            lines.append(f"S4 strength {s:g} AHEAD of strength {o:g}: {str(ok).upper()} (fold E10.5; "
                         + "; ".join(f"{w} {_st(fold(w, arm(s), arm(o)))}" for w, _ in WORLDS) + "; the same bar)")
    return lines


def follows(res: dict) -> str:
    o = res["outcome"]
    if o == "VOID":
        return FOLLOWS["VOID"].format(why="; ".join(res["why"]))
    if o == "AHEAD":
        return FOLLOWS["AHEAD"].format(s=res["output"], b=budget_at_test(res["output"]))
    return FOLLOWS[o]


def reading(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], commit: str | None) -> dict:
    try:
        why, dist = void_reasons(parts, h27_parts, commit)
    except (KeyError, TypeError, ValueError, IndexError) as e:  # a malformed part
        why, dist = [f"a part cannot be read ({type(e).__name__}: {e})"], {}
    if why:
        res = {"void": True, "why": why, "guards": dist, "commit": commit, "outcome": "VOID", "output": None}
        res["follows"] = follows(res)
        return res
    rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
    tab = tables(rows)
    st = statements(tab)
    res = {"void": False, "why": [], "row": ROW, "commit": parts[TASK_KEYS[0]][0]["provenance"]["commit"],
           "guards": dist, "tables": tab, "statements": {f"{s:g}": v for s, v in st.items()},
           "output": output(st), "outcome": outcome(st, False),
           "budgets": {f"{s:g}": {"test": budget_at_test(s), **{k: s * b for k, b in PRACTICE_BUDGET.items()}}
                       for s in STRENGTHS},
           "rendering_offset": {w: tab[w]["fold E10.5"][pname(ZERO, CONTROL)] for w, _ in WORLDS},
           "s0_x_max_abs_diff_control": max(r["checks"]["s0_x_max_abs_diff_control"] for k in TASK_KEYS
                                            for r in rows[k]),
           "trend_record": {k: {"seed0": rows[k][0]["records"][arm(1.0)]["trend"],
                                "whole_pool": rows[k][0]["checks"].get("whole_pool")} for k in TASK_KEYS},
           "not_built": NOT_BUILT_LINE}
    res["follows"] = follows(res)
    res["lines"] = statement_lines(st, tab)
    return res


def text(res: dict) -> str:
    if res["void"]:
        return res["follows"] + "\n"
    g = res["guards"]["V4"]
    lines = [L0.format(commit=res["commit"]), L1,
             f"GATE (the control against row H27's committed {H27_ARM}): largest distance "
             f"{max(v['mean'] for v in g.values()):.3f} (ten-seed mean) and {max(v['seed'] for v in g.values()):.3f} "
             f"(one seed), board points; void beyond {GUARD_MEAN} and {GUARD_SEED}.", res["not_built"],
             "J6 per strength: " + " | ".join(j6_line(s) for s in STRENGTHS),
             "RENDERING OFFSET (trend_s0 minus the control, fold E10.5; expression re-rendered from unchanged counts): "
             + "; ".join(f"{w} {_st(v['score'])}" for w, v in res["rendering_offset"].items())
             + f"; largest |X difference| {res['s0_x_max_abs_diff_control']:.3g}.",
             *res["lines"], f"OUTCOME: {res['outcome']}" + (f" (strength {res['output']:g})" if res["output"] else ""),
             res["follows"], "TABLES (paired by seed, mean [min, max] seeds ahead of 10 (se); per-term mean skill "
             "differences):"]
    for w, _ in WORLDS:
        for u in UNITS:
            for p, c in res["tables"][w][u].items():
                lines.append(f"{TAG} {w} {u} {p}: {_st(c['score'])}; {_terms(c)}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def cmd_run(a) -> None:
    refuse(Path(a.out) / PART)
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS} ({NOT_BUILT_LINE})")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data))


def cmd_read(a) -> None:
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its read)")
    out = Path(a.out)
    if (out / "reading.json").exists() or (out / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")
    res = reading(load_parts(Path(a.runs)), load_h27(Path(a.h27)), a.commit)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "reading.txt").write_text(text(res))
    print(f"OUTCOME: {res['outcome']}" + (f" (strength {res['output']:g})" if res.get("output") else ""))
    print(res["follows"])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", required=True, help=f"one of {TASK_KEYS}")
    r.add_argument("--seeds", required=True, help="a range such as 0-1")
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--h27", default=str(H27_PARTS))
    d.add_argument("--out", required=True)
    d.add_argument("--commit", required=True, help="row H28's own commit (a prefix of the parts' commit)")
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
