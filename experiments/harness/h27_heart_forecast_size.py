#!/usr/bin/env python
"""Try H27 of experiments/harness/TRIES.md: row H26's three heart forms, each sized as the builder sizes it at the
test (raw mode, by the size procedure's forecast formed with the target held out), on row H26's three released
transitions, five truths and ten seeds. A labelled practice analogue that prices; it decides nothing.

The forms are row H26's (h26_heart_scope_forms, imported, not copied): A full_iso_half (the whole reference stage),
B half_iso_half (radial_trim_mask at 0.5), C core_iso_half (radial_trim_mask at 0.0). Each is `vec.harness._copy` on a
context whose reference input is the form's pool, then a size step, then row H25's half (`h25.kept_rows`), then the
checked builder. Only the size step differs between the three sizings of one form, built in one process from one draw
and one set of half rows:
  o  ORACLE TWIN (reference, not counted): row H26's own arm, `h26.form` at row H26's core-ratio factor
     rms(truth core) / rms(reference core). The gate (it must reproduce row H26's committed scores) and the link to
     row H26's table. A convention of row H26, not the builder's size.
  p  PERFECT TWIN (counted): `t2_build.rescale_to_rms` to the size a perfect forecast of the as-released stage would
     give the builder's field: A and B the truth's whole-cloud RMS radius, C the truth's core RMS radius. An oracle
     control, never a centre.
  f  FORECAST ARM (counted): `t2_build.apply_size_procedure(coords, size, hyp, z_mode="raw", expect_board="heart")`,
     hyp H1 for A, H2 for B, "core" for C, on an in-memory size result formed per transition by the size procedure's
     own functions (experiments/coords/size_procedure.py, imported) on the transition's INPUT stages only:
     `measure_stage` at the size file's days (experiments/coords/run_p3.sh), then
       two or more inputs: `heart_size_from_geometry(z_mode="raw", min_start_t=-inf, tie_tol=inf)`, so all three
                           pre-registered rules are averaged, as at the test (where no pseudo-target is eligible);
       one input:          the copy_last analogue (`rule_copy_last` on the core, `scope_log_ratio` for H0, H1, H2).
Departures from the builder at the test (each also in reading.json): (a) the E9.0 rate cut-off is lifted; (b) rule
selection is forced to all three rules (tie_tol = inf); (c) copy_last on E8.25_late -> E8.5 (one input stage); (d) an
in-memory result through apply_size_procedure, size_or_hold's crosscheck redone by guard V7 on the untrimmed reference
stage; (e) the practice board's count rule sets the cell number; (f) the preregistered field holds the raw result.

  run   one transition, a seed range, one JSON under --out. Refused before anything loads unless row H27 is in
        TRIES.md, the tree is clean by vec.provenance's rule and the output file does not exist. Prints no score.
  read  applies the void conditions and the statements fixed in row H27, writes reading.json and reading.txt.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from vec import harness, prep, provenance, t2_build
from vec.adoption import SEEDS
from vec.boards import load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parents[2]  # this checkout


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h26 = _load("h26_heart_scope_forms", HERE / "experiments" / "harness" / "h26_heart_scope_forms.py")
h25 = h26.h25  # row H25's halving and guard tolerances, as row H26 imports them
payoff = h26.payoff
sp = h25.sp  # experiments/coords/size_procedure.py, as row H25 imports it

ROW = "H27"
ROW_ID = re.compile(r"[A-Z]{1,2}\d+[a-z]*")
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h27-heart-forecast-size"
PART = "h27.json"
SOURCE = "experiments/harness/h27_heart_forecast_size.py (try H27)"
TASK_KEYS, FOLDS, FOLD_NAMES = h26.TASK_KEYS, h26.FOLDS, h26.FOLD_NAMES
TRUTHS, OWN, CLEAN = h26.TRUTHS, h26.OWN, h26.CLEAN
FORMS = h26.ARMS  # full_iso_half (A), half_iso_half (B), core_iso_half (C)
LETTER = h26.LETTER
PAIRS = h26.PAIRS
SIZINGS = ("o", "p", "f")
COUNTED_SIZINGS = ("p", "f")
HYP = {"full_iso_half": "H1", "half_iso_half": "H2", "core_iso_half": t2_build.SIZE_CORE_SCOPE}
HEART_RULES = tuple(sp.PREREGISTERED["heart_rules"])  # latest_rate, half_rate, power_t
COPY_RULE = "copy_last"
MIN_START_T = -np.inf  # departure (a): the E9.0 cut-off lifted
TIE_TOL = np.inf  # departure (b): every pre-registered rule kept, as at the test
# the size file's days for the input stages (experiments/coords/run_p3.sh, step 4); a test parses run_p3.sh
SIZE_FILE_DAYS = {"T2/E8.25_late.h5ad": 8.25, "T2/E8.5.h5ad": 8.5, "T2/E8.75.h5ad": 8.75, "T2/E9.5.h5ad": 9.5}
ONE_STAGE_KEY = "T2h/E8.25_late->E8.5"  # payoff.HEART_EXTRA_TASKS: one input stage, the copy_last analogue
TEST_INTERVAL, TEST_GAP, TEST_H = ("E9.5", "E10.5"), 2.0, 2.0
ANALOGUE_RATES = "size procedure, E9.0 cut-off lifted, all three rules averaged"
ANALOGUE_COPY = "copy_last analogue, one input stage"
DEPARTURES = {
    "a": "the E9.0 rate cut-off is lifted (min_start_t = -inf): every practice rate interval starts before E9.0, so "
         "the pre-registered settings give valid=False and no figure (row H25)",
    "b": "rule selection is forced to all three pre-registered rules (tie_tol = inf), the test's outcome, reached at "
         "the test by no eligible pseudo-target",
    "c": "copy_last on E8.25_late -> E8.5, the one transition with one input stage",
    "d": "an in-memory size result applied through apply_size_procedure; size_or_hold's crosscheck (hardwired to "
         "E10.5) is redone by guard V7 on each transition's untrimmed reference stage",
    "e": "the practice board's count rule sets the cell number, not the test board's",
    "f": "the preregistered field holds the raw result; the builder passes --z-mode raw, so the recorded mode never "
         "enters",
}
H26_PARTS = HERE / "experiments" / "harness" / "results" / "h26_heart_scope_forms" / "parts"
GEOMETRY_CSV = HERE / "experiments" / "coords" / "results" / "stage_geometry.csv"
MODES_CSV = HERE / "experiments" / "coords" / "results" / "size_heart_modes.csv"
V6_STAGES = ("E8.25_late", "E8.75", "E9.5")
V6_FIELDS = ("rms_full", "rms_core", "log_scope")
V6_PSEUDO = ("E8.25_late", "E8.75", 9.5)  # inputs and target day of the committed raw no-cut-off row
GUARD_MEAN, GUARD_SEED = h25.GUARD_MEAN, h25.GUARD_SEED  # 0.1 and 0.5 board points
TOL = h26.TOL  # 1e-9
NBHD = "neighborhood_mmd"
NBHD_TOL = h25.NBHD_TWIN_TOL  # 1e-3
EXPRESSION_ONLY = h25.EXPRESSION_ONLY
FACTOR_RTOL = h26.FACTOR_RTOL  # 1e-12
V6_RTOL = 1e-9
SIZE_LOG_TOL = 1e-9
SHOWN_METRICS = ("scale_log_ratio", "d2_shape", "occupancy_dice")
TAG = "[practice analogue, released stages]"
LABELS = {"f": "practice analogue, forecast size", "p": "perfect size (oracle control)",
          "o": "row H26 core-ratio size"}


def arm_name(form: str, sizing: str) -> str:
    return f"{form}:{sizing}"


ARM_NAMES = tuple(arm_name(f, s) for f in FORMS for s in SIZINGS)
ALL = (harness.FLOOR, *ARM_NAMES)


# ------------------------------------------------------------------------------------------------ refusals
def tries_rows(path: Path = TRIES) -> set[str]:
    """The row IDs of TRIES.md's tables (each table row's first word, where it is a row ID)."""
    words = (ln[2:].split("|", 1)[0].split() for ln in Path(path).read_text(encoding="utf-8").splitlines()
             if ln[:2] == "| ")
    return {w[0] for w in words if w and ROW_ID.fullmatch(w[0])}


def refuse(path: Path) -> None:
    """The refusals before anything loads: no row H27, an unclean tree (vec.provenance's rule), an existing output."""
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    _, unclean = provenance.git_state(HERE)
    if unclean is None or unclean:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    if Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


# ------------------------------------------------------------------------------------------------ the forecast
def input_files(task) -> list[tuple[str, float]]:
    """The task's input stages (never its target) with the size file's days, earliest first."""
    files = sorted(set(task.inputs.values()), key=lambda f: SIZE_FILE_DAYS[f])
    if task.target in files:
        raise ValueError(f"{task.key}: the target is an input")
    return [(f, SIZE_FILE_DAYS[f]) for f in files]


def input_geometries(task, data_dir: Path) -> tuple[list, list[dict]]:
    """sp.measure_stage on each input stage (duplicated records dropped, as the size file measures). The target's path
    is never passed."""
    geoms, infos = [], []
    for f, day in input_files(task):
        g, info = sp.measure_stage(Path(data_dir) / f, day)
        geoms.append(g)
        infos.append(info)
    return geoms, infos


def _stage_record(g) -> dict:
    return asdict(g) | {"log_scope": g.log_scope, "log_scope_xy": g.log_scope_xy,
                        "core_thin_axis_angle_deg": g.core_thin_axis_angle_deg}


def copy_analogue(g, target_t: float) -> dict:
    """One input stage: the procedure's copy_last rule on the core and its scope hypotheses (H1 = H0 = s, H2 = s / 2),
    in the layout of heart_size_from_geometry's raw result."""
    core_log = sp.rule_copy_last(np.array([g.t]), np.array([np.log(g.rms_core)]), target_t)
    core = float(np.exp(core_log))
    scope = {}
    for h in sp.PREREGISTERED["heart_scope_hypotheses"]:
        s = sp.scope_log_ratio([g.log_scope], [g.t], target_t, h)
        full = core * float(np.exp(s))
        scope[h] = {"pred_log_scope": s, "pred_rms_full": full, "scale_vs_reference_full": full / g.rms_full}
    return {"target_t": float(target_t), "z_mode_requested": "raw", "z_mode": "raw", "rules": [COPY_RULE],
            "min_start_t": None, "tie_tol": None, "reference": g.name, "stages": [_stage_record(g)],
            "growth": {"core": {"valid": True, "applicable": [COPY_RULE], "chosen": [COPY_RULE],
                                "pred_log_target": core_log}},
            "scope": scope, "valid": True, "pred_rms_core": core,
            "pred_log_growth_core": float(np.log(core / g.rms_core))}


def analogue_size_result(geoms, target_t: float) -> dict:
    """The in-memory size result of one transition, in the layout apply_size_procedure reads. Never written as a size
    file."""
    geoms = sorted(geoms, key=lambda g: g.t)
    if len(geoms) >= 2:
        res = sp.heart_size_from_geometry(geoms, target_t, z_mode="raw", min_start_t=MIN_START_T, tie_tol=TIE_TOL)
        label = ANALOGUE_RATES
    elif len(geoms) == 1:
        res, label = copy_analogue(geoms[0], target_t), ANALOGUE_COPY
    else:
        raise ValueError("no input stage")
    return {"board": "heart", "valid": res["valid"], "target_t": float(target_t),
            "timestamp": "try H27 analogue (in memory, never a size file)", "preregistered": res,
            "z_modes": {"raw": res},
            "provenance": {"analogue": label, "git_commit": None, "size_procedure_sha256": None,
                           "departures": sorted(DEPARTURES)}}


def forecast_targets(size: dict) -> dict[str, float]:
    """The field each form's builder path reads in raw mode: z_modes.raw.scope.H1/H2.pred_rms_full, z_modes.raw
    .pred_rms_core."""
    r = size["z_modes"]["raw"]
    return {form: float(r["pred_rms_core"] if HYP[form] == t2_build.SIZE_CORE_SCOPE
                        else r["scope"][HYP[form]]["pred_rms_full"]) for form in FORMS}


def forecast_record(size: dict, target_t: float) -> dict:
    """What L2 prints of the analogue result: inputs, rules, rate interval, gap and the procedure's horizon h."""
    r = size["z_modes"]["raw"]
    st = r["stages"]
    names, ts = [g["name"] for g in st], [float(g["t"]) for g in st]
    two = len(st) >= 2
    return {"inputs": names, "days": ts, "chosen": list(r["growth"]["core"]["chosen"]),
            "analogue": size["provenance"]["analogue"], "min_start_t": r.get("min_start_t"),
            "tie_tol": r.get("tie_tol"), "interval": names[-2:] if two else None, "gap": float(target_t - ts[-1]),
            "h": float((target_t - ts[-1]) / (ts[-1] - ts[-2])) if two else None, "target_t": float(target_t)}


def perfect_targets(TC: np.ndarray, tcore: np.ndarray) -> dict[str, float]:
    """The p twins' sizes: a perfect forecast of the as-released stage, read through each form's field."""
    full, core = t2_build.rms_radius(TC), t2_build.rms_radius(TC[tcore])
    return {"full_iso_half": full, "half_iso_half": full, "core_iso_half": core}


def convention_offsets(ptargets: dict, even: float, pool_rms: dict, ref_core_rms: float) -> dict:
    """log(p target / o's nominal size), o's nominal size = even factor x the form's pool RMS radius. C's is
    target-free: log rms(ref core) - log rms(C pool)."""
    out = {form: float(np.log(ptargets[form] / (even * pool_rms[form]))) for form in FORMS}
    out["C_parts"] = {"log_rms_ref_core": float(np.log(ref_core_rms)),
                      "log_rms_c_pool": float(np.log(pool_rms["core_iso_half"]))}
    return out


def v7_crosscheck(ref_coords: np.ndarray, size: dict) -> dict:
    """size_or_hold's crosscheck on the transition's reference: the untrimmed, drop-duplicates reference cloud against
    the analogue result's measured rms_full of its reference stage, within SIZE_CROSSCHECK_LOG_TOL in log."""
    r = size["z_modes"]["raw"]
    g = {s["name"]: s for s in r["stages"]}[r["reference"]]
    mine = t2_build.rms_radius(ref_coords)
    d = float(np.log(mine / g["rms_full"]))
    return {"reference": r["reference"], "loaded_rms_full": mine, "size_rms_full": float(g["rms_full"]),
            "log_diff": d, "ok": bool(abs(d) <= t2_build.SIZE_CROSSCHECK_LOG_TOL)}


# ------------------------------------------------------------------------------------------------ the arms
def size_arm(X: np.ndarray, C: np.ndarray, form: str, sizing: str, factors: dict, size: dict, ptargets: dict,
             seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    """One sizing of one form's drawn cells, then row H25's half. o: row H26's form; p: rescale_to_rms to the perfect
    target; f: apply_size_procedure in raw mode with the form's hypothesis."""
    C = np.asarray(C, dtype=np.float64)
    before = t2_build.rms_radius(C)
    if sizing == "o":
        Xh, Ch = h26.form(X, C, form, factors, seed)
        return Xh, Ch, {"sizing": "o", "factor": float(factors["even"]), "rms_before": before,
                        "rms_after": before * float(factors["even"])}
    if sizing == "p":
        T = float(ptargets[form])
        Cs = t2_build.rescale_to_rms(C, T)
        after = t2_build.rms_radius(Cs)
        rec = {"sizing": "p", "target_rms_full": T, "rms_before": before, "rms_after": after,
               "log_rms_after_vs_target": float(np.log(after / T))}
    elif sizing == "f":
        Cs, r = t2_build.apply_size_procedure(C, size, HYP[form], z_mode="raw", expect_board="heart")
        rec = {"sizing": "f", **{k: r.get(k) for k in ("hypothesis", "z_mode_used", "target_rms_full", "rms_before",
                                                         "rms_after", "log_rms_after_vs_target", "target_source",
                                                         "warnings")}}
    else:
        raise ValueError(f"unknown sizing {sizing!r}")
    rows = h25.kept_rows(len(Cs), seed, 1)
    return np.asarray(X)[rows], Cs[rows], rec


def build_forms(form: str, ctx: harness.BuildContext, factors: dict, size: dict, ptargets: dict) -> dict:
    """One form on the context whose reference input is its pool: one harness draw, three sizings, the same half,
    each through the checked builder."""
    st = harness._copy(ctx)
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    out = {}
    for sizing in SIZINGS:
        X, C, rec = size_arm(st.X, st.coords, form, sizing, factors, size, ptargets, seed)
        sub, rep = t2_build.build_submission(ctx.board, X, C, genes=list(ctx.board.genes))
        if not rep.ok:
            raise ValueError(f"{ctx.task.key} {form} {sizing}: the built file fails its check: {rep.errors}")
        out[sizing] = (Stage.from_anndata(sub, ctx.board.genes, need_coords=True), rec)
    return out


# ------------------------------------------------------------------------------------------------ the run
def run_transition(key: str, seeds, out: Path, data_dir: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse(path)
    t = payoff.heart_task(key)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    # the forecast first, from the input stages only
    geoms, infos = input_geometries(t, data_dir)
    target_t = float(harness.STAGE_DAYS[t.target])
    size = analogue_size_result(geoms, target_t)
    ftargets = forecast_targets(size)
    frec = forecast_record(size, target_t)
    # the truths and pools, as row H26 loads them
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    tm = h26.truth_masks(TC, tcore)
    truths = {q: Stage.from_anndata(ta[m], genes, need_coords=True) for q, m in tm.items()}
    t_info = {q: h26.truth_info(TC, tcore, m) for q, m in tm.items()}
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    pool = inputs["ref"]
    pcore = np.asarray(pool.obs["core"], dtype=bool)
    PC = np.asarray(pool.obsm[COORD_KEY], dtype=float)[:, :3]
    pm = {form: m for form, m in h26.pool_masks(PC, pcore).items() if form in FORMS}
    pools = {form: (pool if h26.TRIM[form] is None else pool[m].copy()) for form, m in pm.items()}
    p_info = {form: h26.pool_info(PC, pcore, m) for form, m in pm.items()}
    pool_rms = {form: t2_build.rms_radius(PC[m]) for form, m in pm.items()}
    factors = h26.oracle_factors(TC, tcore, PC, pcore)
    ptargets = perfect_targets(TC, tcore)
    offsets = convention_offsets(ptargets, factors["even"], pool_rms, t2_build.rms_radius(PC[pcore]))
    v7 = v7_crosscheck(PC, size)
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        preds, recs, n_cells = {}, {}, {}
        for form in FORMS:
            ctx = harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pools[form]},
                                       seed=harness.BUILD_SEED_OFFSET + s)
            with harness.guarded_reads(harness.REPO_ROOT, data_dir):
                built = build_forms(form, ctx, factors, size, ptargets)
            for sizing, (st, rec) in built.items():
                preds[arm_name(form, sizing)] = st
                recs[arm_name(form, sizing)] = rec
                n_cells[arm_name(form, sizing)] = int(st.n)
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": n_cells, "size_records": recs,
               "scores": {}, "skills": {}, "groups": {}}
        for q, truth in truths.items():
            board = replace(base_board, truth_cells=truth.n)
            fl = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, s), truth, ref_stage, seed=s,
                             anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
            row["scores"][q] = {harness.FLOOR: fl.score}
            row["skills"][q] = {harness.FLOOR: fl.skills}
            row["groups"][q] = {harness.FLOOR: dict(fl.group_skill)}
            for arm, pred in preds.items():
                r = score_board(board, pred, truth, ref_stage, seed=s, anchors=fl.anchors, build=harness.BUILD,
                                protocol=harness.PROTOCOL)
                row["scores"][q][arm], row["skills"][q][arm] = r.score, r.skills
                row["groups"][q][arm] = dict(r.group_skill)
        row["seconds"] = round(time.time() - t0, 1)
        rows_out.append(row)
        print(f"{key} seed {s} ({row['seconds']} s) done", flush=True)
    res = {"label": "try H27: practice analogue (forecast size f, perfect size p, row H26 core-ratio size o); not an "
                    "adoption result", "axis": AXIS, "task": t.key, "arms": list(ARM_NAMES),
           "truth_values": list(TRUTHS),
           "size_result": size, "input_records": infos, "forecast": frec, "forecast_targets": ftargets,
           "perfect_targets": ptargets, "oracle_factors": factors, "pool_rms": pool_rms, "offsets": offsets,
           "truth_rms": {"full": ptargets["full_iso_half"], "core": ptargets["core_iso_half"]},
           "v7": v7, "truth_info": t_info, "pred_info": p_info, "truth_seed": payoff.TRUTH_SEED,
           "half_seed_base": h25.HALF_SEED_BASE, "departures": DEPARTURES, "seeds": list(seeds), "rows": rows_out,
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


def load_h26(root: Path = H26_PARTS) -> dict[str, list[dict]]:
    return h26.load_parts(root)


def load_committed() -> dict:
    """The committed pre-release size results V6 checks against: stage_geometry.csv's rows and size_heart_modes.csv's
    raw no-cut-off procedure row (predictions only; no released target is read from them)."""
    with open(GEOMETRY_CSV, newline="") as f:
        geom = {r["name"]: {k: float(r[k]) for k in V6_FIELDS} for r in csv.DictReader(f) if r["name"] in V6_STAGES}
    with open(MODES_CSV, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["z_mode"] == "raw" and r["settings"] == "no_rate_cutoff_diagnostic"
                and r["rule"].startswith("procedure")]
    if len(rows) != 1:
        raise ValueError("size_heart_modes.csv: expected one raw no_rate_cutoff_diagnostic procedure row")
    return {"geometry": geom, "modes_raw": {"pred_rms_core": float(rows[0]["pred_rms_core"]),
                                            "full_H1": float(rows[0]["full_H1"])}}


# ------------------------------------------------------------------------------------------------ the guards
def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else abs(a - b)


def _same(mine: dict | None, theirs: dict | None, cells: str) -> bool:
    """V3's equality: the same cell number and the same core share, exactly."""
    return bool(mine and theirs and mine[cells] == theirs[cells] and mine["core_share"] == theirs["core_share"])


def _geometry(d: dict):
    return sp.StageGeometry(**{k: d[k] for k in sp.StageGeometry.__dataclass_fields__})


def v6_reproduction(parts: dict[str, list[dict]], committed: dict) -> tuple[list[str], dict]:
    """V6: the parts' measured geometries of E8.25_late, E8.75 and E9.5 against stage_geometry.csv, and the procedure on
    [E8.25_late, E8.75] to t = 9.5 (raw, cut-off lifted) against size_heart_modes.csv's committed predictions, at the
    pre-registered tie_tol and at tie_tol = inf."""
    why, dist = [], {"geometry": 0.0, "procedure": None}
    seen: dict[str, dict] = {}
    for ps in parts.values():
        for p in ps:
            for g in p["size_result"]["z_modes"]["raw"]["stages"]:
                if g["name"] in V6_STAGES:
                    seen.setdefault(g["name"], g)
                    for k in V6_FIELDS:
                        e = _rel(float(g[k]), committed["geometry"][g["name"]][k])
                        dist["geometry"] = max(dist["geometry"], e)
                        if not e <= V6_RTOL:
                            why.append(f"V6: {g['name']} {k} differs from stage_geometry.csv (relative {e:.2e})")
    gone = [n for n in V6_STAGES if n not in seen]
    if gone:
        why.append(f"V6: no part measured {gone}")
        return why, dist
    a, b, t_star = V6_PSEUDO
    geoms = [_geometry(seen[a]), _geometry(seen[b])]
    want = committed["modes_raw"]
    errs = []
    for tie in (None, TIE_TOL):
        r = sp.heart_size_from_geometry(geoms, t_star, z_mode="raw", min_start_t=MIN_START_T, tie_tol=tie)
        errs += [_rel(r["pred_rms_core"], want["pred_rms_core"]), _rel(r["scope"]["H1"]["pred_rms_full"],
                                                                       want["full_H1"])]
    dist["procedure"] = max(errs)
    if not dist["procedure"] <= V6_RTOL:
        why.append(f"V6: the procedure on {a}, {b} does not reproduce size_heart_modes.csv (relative "
                   f"{dist['procedure']:.2e})")
    return why, dist


def v8_size_results(key: str, part: dict) -> list[str]:
    why = []
    size = part["size_result"]
    r = size["z_modes"]["raw"]
    if size.get("valid") is not True or r.get("valid") is not True:
        why.append(f"V8: {key}: the analogue result is not valid")
    chosen = list(r["growth"]["core"]["chosen"])
    want = [COPY_RULE] if key == ONE_STAGE_KEY else list(HEART_RULES)
    if chosen != want:
        why.append(f"V8: {key}: rules chosen {chosen}, not {want}")
    for name, vals in (("forecast", part["forecast_targets"]), ("perfect", part["perfect_targets"])):
        for form in FORMS:
            v = vals.get(form)
            if not (isinstance(v, int | float) and np.isfinite(v) and v > 0):
                why.append(f"V8: {key}: {name} target of {form} is not a positive finite number")
    for row in part["rows"]:
        for form in FORMS:
            fr, pr = row["size_records"][arm_name(form, "f")], row["size_records"][arm_name(form, "p")]
            if not abs(fr["log_rms_after_vs_target"]) <= SIZE_LOG_TOL:
                why.append(f"V8: {key} seed {row['seed']}: {form} f is not at its target")
            if not abs(float(fr["target_rms_full"]) - part["forecast_targets"][form]) <= SIZE_LOG_TOL * abs(
                    part["forecast_targets"][form]):
                why.append(f"V8: {key} seed {row['seed']}: {form} f's target is not the form's field")
            if not abs(pr["log_rms_after_vs_target"]) <= SIZE_LOG_TOL:
                why.append(f"V8: {key} seed {row['seed']}: {form} p is not at its target")
    return why


def v4_twins(key: str, rows: list[dict], h26_rows: list[dict]) -> tuple[list[str], dict]:
    """V4: each oracle twin against row H26's committed score, per truth: ten-seed mean within GUARD_MEAN, every seed
    within GUARD_SEED."""
    theirs = {r["seed"]: r for r in h26_rows}
    why, big = [], {"mean": 0.0, "seed": 0.0}
    for q in TRUTHS:
        for form in FORMS:
            try:
                d = np.array([r["scores"][q][arm_name(form, "o")] - theirs[r["seed"]]["scores"][q][form]
                              for r in rows], float)
            except KeyError:
                why.append(f"V4: {key} {q} {form}: no committed row H26 score to compare")
                continue
            big["mean"], big["seed"] = max(big["mean"], abs(d.mean())), max(big["seed"], float(np.abs(d).max()))
            if abs(d.mean()) > GUARD_MEAN:
                why.append(f"V4: {key} {q} {form}: the oracle twin's ten-seed mean is {d.mean():+.3f} from row H26's")
            if np.abs(d).max() > GUARD_SEED:
                why.append(f"V4: {key} {q} {form}: a seed is {np.abs(d).max():.3f} from row H26's")
    return why, big


def v5_identity(key: str, rows: list[dict]) -> tuple[list[str], dict]:
    """V5: each p and f arm against its o twin: the expression-only skills to TOL, neighborhood_mmd to NBHD_TOL, the
    same cell number."""
    why, big = [], {"expression": 0.0, "neighbourhood": 0.0}
    for r in rows:
        for form in FORMS:
            o = arm_name(form, "o")
            for sizing in COUNTED_SIZINGS:
                a = arm_name(form, sizing)
                if r["n_cells"][a] != r["n_cells"][o]:
                    why.append(f"V5: {key} seed {r['seed']}: {a} has another cell number than {o}")
                for q in TRUTHS:
                    mine, twin = r["skills"][q][a], r["skills"][q][o]
                    e = max(abs(mine[m] - twin[m]) for m in EXPRESSION_ONLY)
                    n = abs(mine[NBHD] - twin[NBHD])
                    big["expression"], big["neighbourhood"] = max(big["expression"], e), max(big["neighbourhood"], n)
                    if e > TOL:
                        why.append(f"V5: {key} seed {r['seed']} {q}: {a}'s expression skills differ from {o}'s")
                    if n > NBHD_TOL:
                        why.append(f"V5: {key} seed {r['seed']} {q}: {a}'s neighborhood_mmd differs from {o}'s")
    return why, big


def void_reasons(parts: dict[str, list[dict]], h26_parts: dict[str, list[dict]], committed: dict,
                 commit: str | None = None) -> tuple[list[str], dict]:
    """Why the read may not be read (empty when it may), and the guard distances."""
    why: list[str] = []
    dist: dict = {"V4": {}, "V5": {}, "V6": None, "V7": {}}
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
        seeds = [r["seed"] for r in rows]
        if seeds != list(SEEDS):
            why.append(f"V1: {key}: seeds {seeds}, not {list(SEEDS)}")
            complete = False
        for r in rows:
            for field in ("scores", "skills", "groups"):
                gone = sorted({f"{q}/{a}" for q in TRUTHS for a in ALL if a not in r.get(field, {}).get(q, {})})
                if gone:
                    why.append(f"V1: {key} seed {r['seed']}: missing from {field}: {gone[:6]}")
                    complete = False
            if set(r.get("size_records", {})) != set(ARM_NAMES) or set(r.get("n_cells", {})) != set(ARM_NAMES):
                why.append(f"V1: {key} seed {r['seed']}: size records or cell numbers missing")
                complete = False
    if not every:
        return why or ["V1: no part"], dist
    commits = {p["provenance"]["commit"] for p in every}  # V2
    if len(commits) > 1:
        why.append("V2: the parts come from different commits")
    if not all(p["provenance"]["tree_clean"] for p in every):
        why.append("V2: a part comes from a dirty tree")
    if not commit or any(not str(p["provenance"]["commit"]).startswith(commit) for p in every):
        why.append(f"V2: a part does not come from the row's commit {commit!r}")
    if not complete:
        return why, dist
    for key in TASK_KEYS:
        ps, theirs = parts[key], h26_parts.get(key)
        if not theirs:
            why.append(f"V3: {key}: no committed row H26 part")
            continue
        ref = theirs[0]
        for p in ps:  # V3
            for q in TRUTHS:
                if not _same(p["truth_info"].get(q), ref["truth_info"].get(q), "cells"):
                    why.append(f"V3: {key}: truth {q} differs from row H26's (cells or core share)")
            for form in FORMS:
                if not _same(p["pred_info"].get(form), ref["pred_info"].get(form), "pool_cells"):
                    why.append(f"V3: {key}: {form}'s pool differs from row H26's (cells or core share)")
            for k in ("even", "in_plane"):
                f, g = p["oracle_factors"][k], ref["oracle_factors"][k]
                if not abs(f - g) <= FACTOR_RTOL * abs(g):
                    why.append(f"V3: {key}: the {k} oracle factor differs from row H26's")
            if not p["v7"]["ok"]:  # V7
                why.append(f"V7: {key}: the reference cloud is not the cloud the analogue result measured")
            dist["V7"][key] = max(dist["V7"].get(key, 0.0), abs(p["v7"]["log_diff"]))
            why += v8_size_results(key, p)  # V8
            for field in ("forecast_targets", "perfect_targets"):  # one forecast and one perfect size per transition
                if p[field] != ps[0][field]:
                    why.append(f"V8: {key}: the parts disagree on the {field.replace('_', ' ')}")
        rows = h26.merge_rows(ps)[0]
        w4, dist["V4"][key] = v4_twins(key, rows, h26.merge_rows(theirs)[0])
        w5, dist["V5"][key] = v5_identity(key, rows)
        why += w4 + w5
    w6, dist["V6"] = v6_reproduction(parts, committed)
    return why + w6, dist


# ------------------------------------------------------------------------------------------------ the tables
UNITS = h26.UNITS  # each transition; fold E10.5; fold E8.5; mean of folds (last)
UNIT_LABEL = {**{k: k for k in TASK_KEYS}, "fold E10.5": "fold E10.5",
              "fold E8.5": "fold E8.5 (copy_last analogue, no growth step)",
              "mean of folds": "mean of folds (the two folds weighed equally)"}
WORLDS = (("own", (OWN,)), ("CLEAN", CLEAN), *((f"truth {q}", (q,)) for q in TRUTHS))
STATED = ("own", "CLEAN")


def _series(rows, folds, truths, form, sizing) -> np.ndarray:
    return h26.series(rows, folds, truths, arm_name(form, sizing))


def cell(rows: dict[str, list[dict]], folds, truths) -> dict:
    """One cell: per sizing the plain means of A, B, C and the paired differences; per form the size terms f-p, p-o,
    f-o; per pair the shifts of the difference (f minus o, f minus p)."""
    s = {(f, z): _series(rows, folds, truths, f, z) for f in FORMS for z in SIZINGS}
    out = {"mean": {z: {LETTER[f]: float(s[f, z].mean()) for f in FORMS} for z in SIZINGS},
           "diff": {z: {f"{LETTER[x]}-{LETTER[y]}": h26._stats(s[x, z] - s[y, z]) for x, y in PAIRS}
                    for z in SIZINGS},
           "size_term": {LETTER[f]: {"f-p": h26._stats(s[f, "f"] - s[f, "p"]), "p-o": h26._stats(s[f, "p"] - s[f, "o"]),
                                     "f-o": h26._stats(s[f, "f"] - s[f, "o"])} for f in FORMS},
           "shift": {f"{LETTER[x]}-{LETTER[y]}": {
               "f-o": h26._stats((s[x, "f"] - s[y, "f"]) - (s[x, "o"] - s[y, "o"])),
               "f-p": h26._stats((s[x, "f"] - s[y, "f"]) - (s[x, "p"] - s[y, "p"]))} for x, y in PAIRS}}
    return out


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: cell(rows, fs, ts) for u, fs in UNITS.items()} for w, ts in WORLDS}


def ranges(tab: dict) -> dict:
    """S0: over the three transitions, per stated world, the range of each difference at the forecast size and of each
    form's forecast error (f minus p)."""
    out = {}
    for w in STATED:
        cells = [tab[w][k] for k in TASK_KEYS]
        out[w] = {"diff_f": {d: [min(c["diff"]["f"][d]["mean"] for c in cells),
                                 max(c["diff"]["f"][d]["mean"] for c in cells)] for d in ("B-A", "C-A", "C-B")},
                  "f-p": {LETTER[f]: [min(c["size_term"][LETTER[f]]["f-p"]["mean"] for c in cells),
                                      max(c["size_term"][LETTER[f]]["f-p"]["mean"] for c in cells)] for f in FORMS}}
    return out


def groups_table(rows: dict[str, list[dict]]) -> dict:
    """Per transition and truth: each arm's mean group skills and its scale_log_ratio, d2_shape and occupancy_dice
    skills, and the size terms' mean differences on the same quantities."""
    out = {}
    for key in TASK_KEYS:
        for q in TRUTHS:
            def val(arm, what, q=q, key=key):
                src = "groups" if what in h26.GROUPS else "skills"
                return np.array([r[src][q][arm][what] for r in rows[key]], float)
            whats = (*h26.GROUPS, *SHOWN_METRICS)
            arms = {arm: {w: float(val(arm, w).mean()) for w in whats} for arm in ARM_NAMES}
            terms = {LETTER[f]: {t: {w: float((val(arm_name(f, t[0]), w) - val(arm_name(f, t[2]), w)).mean())
                                     for w in whats} for t in ("f-p", "p-o", "f-o")} for f in FORMS}
            out.setdefault(key, {})[q] = {"arms": arms, "size_terms": terms}
    return out


def released_figures(parts: dict[str, list[dict]]) -> dict:
    """L2, per transition, in the procedure's units: the forecast's log errors against this released target (core
    against the truth's core RMS, H1 and H2 full against the as-released truth's full RMS) and the convention
    offsets."""
    out = {}
    for key in TASK_KEYS:
        p = parts[key][0]
        r = p["size_result"]["z_modes"]["raw"]
        tr = p["truth_rms"]
        out[key] = {"forecast": p["forecast"], "forecast_targets": p["forecast_targets"],
                    "perfect_targets": p["perfect_targets"], "oracle_factors": p["oracle_factors"],
                    "log_error": {"core": float(np.log(r["pred_rms_core"] / tr["core"])),
                                  "full_H1": float(np.log(r["scope"]["H1"]["pred_rms_full"] / tr["full"])),
                                  "full_H2": float(np.log(r["scope"]["H2"]["pred_rms_full"] / tr["full"]))},
                    "offsets": p["offsets"], "size_result": p["size_result"], "departures": p["departures"]}
    return out


H26_LINE = ("Beside it, row H26 at its core-ratio oracle size, mean of folds: as released B-A -1.00, C-A -3.86, "
            "C-B -2.85; CLEAN B-A +5.78, C-A +11.79, C-B +6.01.")
H26_FIGURES = {"own": {"B-A": -1.00, "C-A": -3.86, "C-B": -2.85}, "CLEAN": {"B-A": 5.78, "C-A": 11.79, "C-B": 6.01}}
VOID_TAIL = ("No forecast-size figure exists under row H27 and there is no rerun under it; row H26's table stands "
             "alone, at its core-ratio oracle size.")


def reading(parts: dict[str, list[dict]], h26_parts: dict[str, list[dict]], committed: dict,
            commit: str | None) -> dict:
    try:
        why, dist = void_reasons(parts, h26_parts, committed, commit)
    except (KeyError, TypeError, ValueError, IndexError) as e:  # a malformed part
        why, dist = [f"a part cannot be read ({type(e).__name__}: {e})"], {}
    if why:
        return {"void": True, "why": why, "guards": dist, "commit": commit}
    rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
    tab = tables(rows)
    return {"void": False, "why": [], "commit": parts[TASK_KEYS[0]][0]["provenance"]["commit"], "row": ROW,
            "labels": LABELS, "decides": "nothing: the read prices; no pick, no outcome, no advice",
            "guards": dist, "transitions": released_figures(parts), "ranges": ranges(tab), "tables": tab,
            "groups": groups_table(rows), "h26_quoted": H26_FIGURES}


# ------------------------------------------------------------------------------------------------ the fixed words
L0 = ("This read prices; it decides nothing. Every figure is a practice analogue on released stages (rule 9, J1a), "
      "and none is to be read beside a test score. No file, rule, size or builder default changes; the file the "
      "captain named stands unless he names another. The AI assistant states no preference.")
L1 = ("Try H27 (commit {commit}). FORECAST SIZE (f): each arm sized as the builder sizes it at the test (raw mode; A "
      "to scope H1's full-cloud forecast, B to H2's, C to the core forecast), by the size procedure's own functions on "
      "the transition's input stages, the target held out, with all three pre-registered rules averaged as at the "
      "test. This is an analogue, not the test's step: the E9.0 rate cut-off is lifted (no released transition has a "
      "rate interval at or after E9.0 before its target), and on E8.25_late to E8.5 (one input stage) the forecast is "
      "the procedure's copy_last rule. PERFECT (p): the same arm at the size a perfect forecast of the as-released "
      "stage would give (an oracle control). ORACLE (o): row H26's arm at its core-ratio size (a convention of row "
      "H26, not the builder's size).")
S0_HEAD = "RANGES FIRST (rule 9 J1a and J6; a range, never a point for the test):"
S0_TAIL = ("Every practice rate interval lies before E9.0, where the size procedure says rates carry no information "
           "about later growth. The test's budget (2 days, h 2) lies outside all three practice budgets (1.0, 1.75, "
           "0.25 days); the nearest in days has its rate interval on the flat core, at h 7. A practice size term is "
           "shown between zero and its value, never as the test's.")
S1_HEAD = ("AT THE FORECAST SIZE, per world, plain means (A B C | paired difference mean [min, max] seeds ahead of "
           "10):")
S2_HEAD = ("SIZE TERM per form, paired by seed, mean [min, max] seeds above zero of 10: forecast error (f minus p), "
           "convention offset (p minus o), and their sum (f minus o, against row H26's core-ratio size, a convention "
           "and not the right size):")
S3_HEAD = ("SHIFT OF EACH COMPARISON (at the forecast size minus at row H26's core-ratio size; and at the forecast "
           "size minus at the perfect size):")
CLOSING = ("No pick and no advice. These tables go to the captain beside row H26's, in these words; the named file "
           "stands unless he names another.")


def _st(s: dict) -> str:
    return f"{s['mean']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10"


def _l2(key: str, t: dict) -> str:
    fc, e, o = t["forecast"], t["log_error"], t["offsets"]
    interval = (f"{fc['interval'][0]} to {fc['interval'][1]}" if fc["interval"] else "none (one input stage)")
    h = f"{fc['h']:.2f}" if fc["h"] is not None else "none"
    return (f"{key}: inputs {', '.join(fc['inputs'])}; rules {', '.join(fc['chosen'])}; rate interval {interval} "
            f"(test: E9.5 to E10.5); horizon {fc['gap']:.2f} days, h {h} (test: 2 days, h 2). Forecast log error on "
            f"this released target, in the procedure's units: core {e['core']:+.3f}, full H1 {e['full_H1']:+.3f}, "
            f"full H2 {e['full_H2']:+.3f}. Convention offset log(p/o): A {o['full_iso_half']:+.3f}, B "
            f"{o['half_iso_half']:+.3f}, C {o['core_iso_half']:+.3f} (C's from the reference alone: log rms(ref "
            f"core) {o['C_parts']['log_rms_ref_core']:.3f} minus log rms(C pool) "
            f"{o['C_parts']['log_rms_c_pool']:.3f}).")


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + ". " + VOID_TAIL + "\n"
    g = res["guards"]
    lines = [L0, L1.format(commit=res["commit"])]
    lines += [_l2(k, res["transitions"][k]) for k in TASK_KEYS]
    lines.append(f"ORACLE TWINS (row H26's arms rebuilt here): largest distance to row H26's committed scores "
                 f"{max(v['mean'] for v in g['V4'].values()):.3f} (ten-seed mean) and "
                 f"{max(v['seed'] for v in g['V4'].values()):.3f} (one seed), board points; void beyond 0.1 and 0.5.")
    lines.append(S0_HEAD)
    for w in STATED:
        r = res["ranges"][w]
        lines.append(f"{w}: over the three transitions, at the forecast size "
                     + ", ".join(f"{d} from {r['diff_f'][d][0]:+.2f} to {r['diff_f'][d][1]:+.2f}"
                                 for d in ("B-A", "C-A", "C-B"))
                     + "; forecast error (f-p) "
                     + ", ".join(f"{x} from {r['f-p'][x][0]:+.2f} to {r['f-p'][x][1]:+.2f}" for x in ("A", "B", "C"))
                     + ".")
    lines.append(S0_TAIL)
    tab = res["tables"]
    lines.append(S1_HEAD)
    for w, _ in WORLDS:
        for u in UNITS:
            c = tab[w][u]
            lines.append(f"{TAG} {w} {UNIT_LABEL[u]}: "
                         + " ".join(f"{x} {v:.2f}" for x, v in c["mean"]["f"].items()) + " | "
                         + "; ".join(f"{d} {_st(s)}" for d, s in c["diff"]["f"].items()))
    lines.append(S2_HEAD)
    for w, _ in WORLDS:
        for u in UNITS:
            c = tab[w][u]["size_term"]
            lines.append(f"{TAG} {w} {UNIT_LABEL[u]}: "
                         + "; ".join(f"{x} f-p {_st(t['f-p'])}, p-o {t['p-o']['mean']:+.2f}, "
                                     f"f-o {t['f-o']['mean']:+.2f}" for x, t in c.items()))
    lines.append(S3_HEAD)
    for w, _ in WORLDS:
        for u in UNITS:
            c = tab[w][u]["shift"]
            lines.append(f"{TAG} {w} {UNIT_LABEL[u]}: f-o "
                         + ", ".join(f"{d} {s['f-o']['mean']:+.2f}" for d, s in c.items()) + "; f-p "
                         + ", ".join(f"{d} {s['f-p']['mean']:+.2f}" for d, s in c.items()))
    lines += [H26_LINE, CLOSING]
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def cmd_run(a) -> None:
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS}")
    refuse(Path(a.out) / PART)
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data))


def cmd_read(a) -> None:
    out = Path(a.out)
    if (out / "reading.json").exists() or (out / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")
    res = reading(load_parts(Path(a.runs)), load_h26(Path(a.h26)), load_committed(), a.commit)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


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
    d.add_argument("--h26", default=str(H26_PARTS))
    d.add_argument("--out", required=True)
    d.add_argument("--commit", required=True, help="the commit row H27 names (a prefix of the parts' commit)")
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
