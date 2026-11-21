#!/usr/bin/env python
"""Try H29a of experiments/harness/TRIES.md: the practice read of two heart candidates, H-T1 (the whole E10.5 pool,
row H27's form A `full_iso_half`) and H-T1B (the half trim, row H27's form B `half_iso_half`), each with the heart's
trend step (b0_trend, row H11) at strength 1 and its placebo. A labelled diagnostic on one truth (rule 9 J5); it
changes no file, rule or default and builds no test file.

Row H28's script (h28_heart_trend_strength, imported, not copied) parameterised by form. For each form F of the two:
  F:control     row H27's arm F at its forecast size f (`F:f`): `vec.harness._copy` on a context whose reference input
                is F's pool (E10.5 whole for A, under `radial_trim_mask(C, core, 0.5)` for B; earlier stages whole),
                then `t2_build.apply_size_procedure` in raw mode with F's size hypothesis (h27.HYP: H1 for A, H2 for B)
                on row H27's in-memory size result from the input stages only, then row H25's half, then the checked
                builder (`h27.size_arm`, imported).
  F:trend_s1    the control with one change, the expression step, made before the size step: b0_trend's expression
                (`h28.trend_expression` at strength 1, which is `vec.harness._b0_t2(ctx, trend=True)` bit for bit; guard
                V5 checks it in every seed).
  F:placebo_s1  the same factors with each shifted lineage's genes moved by row H28's derangement within ten level bins
                (seed 20261120).
Row H28's trend_s0 and its strengths 0.5 and 0.25 are not built here. The draw, coordinates, size step, half rows and
cell number are the control's in every arm of a form (guard V5).

Transitions, truths and seeds are row H28's: the two harness backtests on E10.5 merged seed by seed into fold E10.5,
the truth as released (own) and row H26's CLEAN world (radial_0, thin_0.5, thin_0), payoff.TRUTH_SEED, seeds 0 to 9.

  run   one transition, a seed range, one JSON under --out. Refused before anything loads unless row H29a is in
        TRIES.md, the tree is clean by vec.provenance's rule, HEAD is --commit and the output file does not exist.
        Prints only "<task> seed <s> (<sec> s) done" and "wrote <path>".
  read  refused unless row H29a is in TRIES.md, the tree is clean and HEAD is --commit, and with an existing reading;
        applies the void conditions, the gates and the fixed statements; prints one outcome line and its fixed
        what-follows line per candidate and nothing else (no figure; void reasons stay in reading.txt); writes the
        numbers to reading.json and reading.txt.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from vec import harness, prep, provenance
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


h28 = _load("h28_heart_trend_strength", HERE / "experiments" / "harness" / "h28_heart_trend_strength.py")
h27, h26, h25, payoff = h28.h27, h28.h26, h28.h25, h28.payoff

ROW = "H29a"
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h29a-heart-trend-forms"
PART = "h29a.json"
SOURCE = "experiments/harness/h29a_heart_trend_forms.py (try H29a)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h29a_heart_trend_forms"
TASK_KEYS = h28.TASK_KEYS  # the two backtests on E10.5
NOT_BUILT, NOT_BUILT_LINE = h28.NOT_BUILT, h28.NOT_BUILT_LINE
OWN, CLEAN, TRUTHS, WORLDS = h28.OWN, h28.CLEAN, h28.TRUTHS, h28.WORLDS
FULL, HALF = "full_iso_half", "half_iso_half"  # row H27's A and B
FORMS = (FULL, HALF)
CANDIDATE = {FULL: "H-T1", HALF: "H-T1B"}
POOL_WORDS = {FULL: "the whole pool", HALF: "the half-trim pool"}
SIZING = h28.SIZING  # "f", row H27's forecast size
STRENGTH = 1.0
STRENGTHS = (STRENGTH,)
CONTROL = h28.CONTROL  # "control"
ARM = h28.arm(STRENGTH)  # "trend_s1"
PLACEBO = h28.placebo(STRENGTH)  # "placebo_s1"
LOCAL = (CONTROL, ARM, PLACEBO)
PERM_SEED = h28.PERM_SEED
PRACTICE_BUDGET = h28.PRACTICE_BUDGET
GUARD_MEAN, GUARD_SEED = h28.GUARD_MEAN, h28.GUARD_SEED  # 0.1 and 0.5 board points
TOL, FACTOR_RTOL = h28.TOL, h28.FACTOR_RTOL
DE_TERMS, TERMS, SHAPE_TERMS = h28.DE_TERMS, h28.TERMS, h28.SHAPE_TERMS
BAR, MIN_SEEDS = h28.BAR, h28.MIN_SEEDS  # +0.3 board points, 8 of 10 seeds
TAG = h28.TAG
QUALIFIED, NOT_QUALIFIED, VOID = "QUALIFIED", "NOT QUALIFIED", "VOID"


def name(form: str, local: str) -> str:
    return f"{form}:{local}"


def h27_arm(form: str) -> str:
    """The gate's committed arm: row H27's form at its forecast size."""
    return h27.arm_name(form, SIZING)


BUILT = tuple(name(f, a) for f in FORMS for a in LOCAL)
ALL = (harness.FLOOR, *BUILT)


# ------------------------------------------------------------------------------------------------ refusals
def tries_rows(path: Path = TRIES) -> set[str]:
    return h27.tries_rows(path)


def refuse_state(commit: str | None, what: str = "run") -> None:
    """No row H29a, an unclean tree (vec.provenance's rule), or a HEAD other than the run commit."""
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its {what})")
    head, unclean = provenance.git_state(HERE)
    if unclean is None or unclean:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    if not commit or head != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit {commit!r}")


def refuse(path: Path, commit: str | None) -> None:
    """The refusals before anything loads."""
    refuse_state(commit)
    if Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


# ------------------------------------------------------------------------------------------------ the run
def form_pools(whole, PC: np.ndarray, pcore: np.ndarray) -> tuple[dict, dict]:
    """Each form's pool as row H27 builds it: the whole reference stage itself for A, the trimmed copy for B."""
    pm = {f: m for f, m in h26.pool_masks(PC, pcore).items() if f in FORMS}
    pools = {f: (whole if h26.TRIM[f] is None else whole[m].copy()) for f, m in pm.items()}
    return pools, {f: h26.pool_info(PC, pcore, m) for f, m in pm.items()}


def build_form(form: str, ctx: harness.BuildContext, size: dict, whole_ref=None) -> tuple[dict, dict]:
    """One form's three arms on the context whose reference input is the form's pool (h28.build_arms with the form,
    strength 1 only and no trend_s0); keys are `name(form, local)`."""
    out, extra = h28.build_arms(ctx, size, whole_ref=whole_ref, form=form, strengths=STRENGTHS, zero=False)
    if set(out) != set(LOCAL):
        raise ValueError(f"{form}: arms {sorted(out)}, not {sorted(LOCAL)}")
    return ({name(form, k): v for k, v in out.items()},
            {"records": {name(form, k): v for k, v in extra["records"].items()}, "checks": extra["checks"]})


def run_transition(key: str, seeds, out: Path, data_dir: Path, commit: str | None) -> dict:
    out = Path(out)
    path = out / PART
    refuse(path, commit)
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
    pools, p_info = form_pools(whole, PC, pcore)
    v7 = h27.v7_crosscheck(PC, size)
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        preds, recs, checks, n_cells = {}, {}, {}, {}
        for form in FORMS:
            ctx = harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pools[form]},
                                       seed=harness.BUILD_SEED_OFFSET + s)
            with harness.guarded_reads(harness.REPO_ROOT, data_dir):
                p, extra = build_form(form, ctx, size, whole_ref=None if h26.TRIM[form] is None else whole)
            preds.update(p)
            recs.update(extra["records"])
            checks[form] = extra["checks"]
            n_cells.update({name(form, k): v for k, v in extra["checks"]["n_cells"].items()})
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": n_cells, "checks": checks,
               "records": recs, "scores": {}, "skills": {}, "groups": {}}
        for q, truth in truths.items():
            board = replace(base_board, truth_cells=truth.n)
            fl = score_board(board, subsample(ref_stage, SUBSAMPLE_FRAC, s), truth, ref_stage, seed=s,
                             anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
            row["scores"][q] = {harness.FLOOR: fl.score}
            row["skills"][q] = {harness.FLOOR: fl.skills}
            row["groups"][q] = {harness.FLOOR: dict(fl.group_skill)}
            for arm_name, pred in preds.items():
                r = score_board(board, pred, truth, ref_stage, seed=s, anchors=fl.anchors, build=harness.BUILD,
                                protocol=harness.PROTOCOL)
                row["scores"][q][arm_name], row["skills"][q][arm_name] = r.score, r.skills
                row["groups"][q][arm_name] = dict(r.group_skill)
        row["seconds"] = round(time.time() - t0, 1)
        rows_out.append(row)
        print(f"{key} seed {s} ({row['seconds']} s) done", flush=True)
    res = {"label": "try H29a: the trend step at strength 1 and its placebo on row H27's forms A (whole pool, H-T1) "
                    "and B (half trim, H-T1B); a diagnostic on one truth, not an adoption result", "axis": AXIS,
           "task": t.key, "arms": list(BUILT), "truth_values": list(TRUTHS), "forecast_targets": ftargets,
           "size_result": size, "v7": v7, "truth_info": t_info, "pred_info": p_info, "truth_seed": payoff.TRUTH_SEED,
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


H27_PARTS = h28.H27_PARTS


def load_h27(root: Path = H27_PARTS) -> dict[str, list[dict]]:
    return h27.load_parts(root)


# ------------------------------------------------------------------------------------------------ the guards
def gate(key: str, form: str, rows: list[dict], h27_rows: list[dict]) -> tuple[list[str], dict]:
    """V4: the form's control against row H27's committed `form:f` (h28.gate's code and tolerances)."""
    why, big = h28.gate(key, rows, h27_rows, mine=name(form, CONTROL), theirs=h27_arm(form))
    return [w.replace("V4: ", f"V4: {form} ", 1) for w in why], big


def identities(key: str, form: str, rows: list[dict]) -> list[str]:
    """V5: strength 1 is the harness's b0_trend on the same context, bit for bit; every arm has the control's draw,
    coordinates (before and after the size step and half) and cell number; no shape term moves (1e-9)."""
    why = []
    names = [name(form, a) for a in LOCAL]
    for r in rows:
        c = r["checks"][form]
        tag = f"V5: {key} {form} seed {r['seed']}"
        if c.get("s1_equals_harness_b0_trend") is not True:
            why.append(f"{tag}: strength 1 is not the harness's b0_trend expression on the same cells")
        for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
            if c.get(k) is not True:
                why.append(f"{tag}: {k.replace('_', ' ')} fails")
        n = c.get("n_cells") or {}
        mismatch = any(r["n_cells"].get(m) != n.get(a) for m, a in zip(names, LOCAL, strict=True))
        if set(n) != set(LOCAL) or len(set(n.values())) != 1 or mismatch:
            why.append(f"{tag}: the arms' cell numbers differ")
        for q in TRUTHS:
            for m in names:
                for t in SHAPE_TERMS:
                    if abs(r["skills"][q][m][t] - r["skills"][q][name(form, CONTROL)][t]) > TOL:
                        why.append(f"{tag} {q}: {m}'s {t} is not the control's")
    return why


def trend_records(key: str, form: str, rows: list[dict], pool_cells: int) -> list[str]:
    """V6: each trend arm's record has the transition's budget at strength 1, the form's pool as its reference, the
    arm and the placebo shift the same lineages, and the placebo leaves no gene on its own factor."""
    why = []
    want = PRACTICE_BUDGET[key]
    for r in rows:
        shifted = set()
        tag = f"V6: {key} {form} seed {r['seed']}"
        for local in (ARM, PLACEBO):
            rec = r["records"][name(form, local)]["trend"]
            if rec.get("strength") != STRENGTH or bool(rec.get("placebo")) != (local == PLACEBO):
                why.append(f"{tag}: {local}'s record is not strength 1 {'placebo' if local == PLACEBO else 'trend'}")
            if abs(rec["budget_days"] - want) > FACTOR_RTOL * want:
                why.append(f"{tag}: {local}'s budget is {rec['budget_days']}, not {want}")
            if rec["ref_cells"] != pool_cells:
                why.append(f"{tag}: {local}'s reference is not the form's pool")
            shifted.add(tuple(rec["shifted"]))
            if local == PLACEBO and any(v for v in rec["fixed_points"].values()):
                why.append(f"{tag}: {local} leaves a gene on its own factor")
        if len(shifted) != 1:
            why.append(f"{tag}: the arm and its placebo shift different lineages")
    return why


def void_reasons(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], commit: str | None
                 ) -> tuple[list[str], dict[str, list[str]], dict]:
    """(reasons that void the whole read, reasons that void one form, gate distances). V1 (completeness), V2 (commit
    and clean trees) and V3's truths and reference cloud are the whole read's; V3's pool and forecast, V4, V5 and V6
    are each form's."""
    why: list[str] = []
    per: dict[str, list[str]] = {f: [] for f in FORMS}
    dist: dict = {"V4": {f: {} for f in FORMS}}
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
            if set(r.get("records", {})) != set(BUILT) or set(r.get("checks", {})) != set(FORMS):
                why.append(f"V1: {key} seed {r['seed']}: records or checks missing")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"V1: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not every:
        return why or ["V1: no part"], per, dist
    if len({p["provenance"]["commit"] for p in every}) > 1:  # V2
        why.append("V2: the parts come from different commits")
    if not all(p["provenance"]["tree_clean"] for p in every):
        why.append("V2: a part comes from a dirty tree")
    if not commit or any(str(p["provenance"]["commit"]) != commit for p in every):
        why.append(f"V2: a part does not come from the row's commit {commit!r}")
    if not complete:
        return why, per, dist
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
            if not p["v7"]["ok"]:
                why.append(f"V3: {key}: the reference cloud is not the cloud the size result measured")
            for form in FORMS:
                if not h27._same(p["pred_info"].get(form), ref["pred_info"].get(form), "pool_cells"):
                    per[form].append(f"V3: {key}: {form}'s pool differs from row H27's (cells or core share)")
                f, g = p["forecast_targets"][form], ref["forecast_targets"][form]
                if not abs(f - g) <= FACTOR_RTOL * abs(g):
                    per[form].append(f"V3: {key}: {form}'s forecast size differs from row H27's")
        rows = h26.merge_rows(ps)[0]
        theirs_rows = h26.merge_rows(theirs)[0]
        for form in FORMS:
            w4, dist["V4"][form][key] = gate(key, form, rows, theirs_rows)
            per[form] += (w4 + identities(key, form, rows)
                          + trend_records(key, form, rows, ps[0]["pred_info"][form]["pool_cells"]))
    return why, per, dist


# ------------------------------------------------------------------------------------------------ the tables
_stats, series, passes, pname = h28._stats, h28.series, h28.passes, h28.pname


def compare(rows, keys, truths, a: str, b: str) -> dict:
    """a minus b, paired by seed: the board score and every term."""
    out = {"score": _stats(series(rows, keys, truths, a) - series(rows, keys, truths, b))}
    for m in (*TERMS, *SHAPE_TERMS):
        out[m] = float((series(rows, keys, truths, a, m) - series(rows, keys, truths, b, m)).mean())
    return out


def pairs(form: str) -> tuple[tuple[str, str], ...]:
    """Reported per form: arm and placebo minus the control; arm minus placebo."""
    return ((name(form, ARM), name(form, CONTROL)), (name(form, PLACEBO), name(form, CONTROL)),
            (name(form, ARM), name(form, PLACEBO)))


UNITS = {**{k: (k,) for k in TASK_KEYS}, "fold E10.5": TASK_KEYS}


def tables(rows: dict[str, list[dict]], form: str) -> dict:
    return {w: {u: {pname(a, b): compare(rows, ks, ts, a, b) for a, b in pairs(form)} for u, ks in UNITS.items()}
            for w, ts in WORLDS}


# ------------------------------------------------------------------------------------------------ the statements
def statements(tab: dict, form: str) -> dict:
    """Row H28's S1 to S3 for one form at strength 1, on fold E10.5, required in both worlds (own and CLEAN)."""
    def fold(w, a, b):
        return tab[w]["fold E10.5"][pname(name(form, a), name(form, b))]
    ahead = all(passes(fold(w, ARM, CONTROL)["score"]) for w, _ in WORLDS)
    over_placebo = all(passes(fold(w, ARM, PLACEBO)["score"]) for w, _ in WORLDS)
    on_de = all(fold(w, ARM, CONTROL)[m] > 0 for w, _ in WORLDS for m in DE_TERMS)
    return {"S1": bool(ahead), "S2": bool(over_placebo), "S3": bool(on_de), "qualified": bool(ahead and over_placebo
                                                                                             and on_de)}


def candidates(st: dict[str, dict | None], whole_void: bool = False) -> dict[str, dict]:
    """The fixed pass statements. st[form] is None when the form is void.
    H-T1 qualifies only if its own S1, S2 and S3 hold on the whole pool.
    H-T1B qualifies only if its own S1, S2 and S3 hold on the half-trim pool AND H-T1's S2 holds on the whole pool."""
    out = {}
    full, half = st.get(FULL), st.get(HALF)
    if whole_void:
        return {c: {"outcome": VOID, "why": "the whole read is void"} for c in ("H-T1", "H-T1B")}
    if full is None:
        out["H-T1"] = {"outcome": VOID, "why": "its form is void"}
    else:
        out["H-T1"] = {"outcome": QUALIFIED if full["qualified"] else NOT_QUALIFIED,
                       "why": "S1, S2 and S3 hold" if full["qualified"]
                       else "fails " + ", ".join(k for k in ("S1", "S2", "S3") if not full[k])}
    if half is None:
        out["H-T1B"] = {"outcome": VOID, "why": "its form is void"}
    else:
        own_ok = half["qualified"]
        t1_s2 = None if full is None else bool(full["S2"])
        fails = [k for k in ("S1", "S2", "S3") if not half[k]]
        if t1_s2 is None:
            fails.append("H-T1's S2 cannot be read (H-T1's form is void)")
        elif not t1_s2:
            fails.append("H-T1's S2 on the whole pool")
        ok = own_ok and t1_s2 is True
        out["H-T1B"] = {"outcome": QUALIFIED if ok else NOT_QUALIFIED,
                        "why": "its own S1, S2 and S3 and H-T1's S2 hold" if ok else "fails " + ", ".join(fails)}
    return out


# ------------------------------------------------------------------------------------------------ the fixed words
FOLLOWS = {
    ("H-T1", QUALIFIED): (
        "H-T1 QUALIFIED under row H29a: trend_s1 is ahead of the copy and of its placebo by the bar and gains on "
        "de_score and de_direction, on fold E10.5 in both worlds, on the whole pool against row H27's "
        "full_iso_half:f (a diagnostic on one truth). What follows: the file enters at its place in the order only if "
        "two DGX builds from one clean commit give identical sha256 and pass the format check, and the captain says "
        "yes to it by name; one priced question goes to him, with the practice figures as ranges over the two "
        "backtests from zero to their practice value (J6), never a point."),
    ("H-T1", NOT_QUALIFIED): (
        "H-T1 NOT QUALIFIED under row H29a ({why}): the file is skipped and the next file moves up. No priced "
        "question; no rerun under this row."),
    ("H-T1", VOID): (
        "H-T1 VOID under row H29a ({why}; reasons in reading.txt): no figure exists for it and there is no rerun "
        "under this row; it does not qualify and is skipped."),
    ("H-T1B", QUALIFIED): (
        "H-T1B QUALIFIED under row H29a: on the half-trim pool against row H27's half_iso_half:f, trend_s1 is ahead "
        "of the copy and of its placebo by the bar and gains on de_score and de_direction, on fold E10.5 in both "
        "worlds, and H-T1's S2 holds on the whole pool (a diagnostic on one truth). What follows: the file enters at "
        "its place in the order only if two DGX builds from one clean commit give identical sha256 and pass the "
        "format check, and the captain says yes to it by name, his yes read as also covering K1's scope for this one "
        "file and the question saying so; one priced question goes to him, with the practice figures as ranges over "
        "the two backtests from zero to their practice value (J6), never a point."),
    ("H-T1B", NOT_QUALIFIED): (
        "H-T1B NOT QUALIFIED under row H29a ({why}): the file is skipped and the next file moves up. No priced "
        "question; no rerun under this row."),
    ("H-T1B", VOID): (
        "H-T1B VOID under row H29a ({why}; reasons in reading.txt): no figure exists for it and there is no rerun "
        "under this row; it does not qualify and is skipped."),
}
L0 = ("Try H29a (commit {commit}). A diagnostic on one truth (E10.5, two backtests merged seed by seed; rule 9 J5); "
      "every figure is a practice reading on released stages and none is to be read beside a test score. Nothing is "
      "built and no file, rule or default changes.")
L1 = ("Arms per form: the control is row H27's form at its forecast size (copy expression; A full_iso_half, H1, the "
      "whole pool, for H-T1; B half_iso_half, H2, radial_trim_mask 0.5, for H-T1B); trend_s1 is b0_trend's expression "
      "before the size step; placebo_s1 is the same factors with each lineage's genes deranged within ten level bins "
      f"(seed {PERM_SEED}).")


def follows(cand: str, c: dict) -> str:
    return FOLLOWS[cand, c["outcome"]].format(why=c["why"])


def reading(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], commit: str | None) -> dict:
    try:
        why, per, dist = void_reasons(parts, h27_parts, commit)
    except (KeyError, TypeError, ValueError, IndexError) as e:  # a malformed part
        why, per, dist = [f"a part cannot be read ({type(e).__name__}: {e})"], {f: [] for f in FORMS}, {}
    res = {"row": ROW, "commit": commit, "why": why, "form_why": per, "guards": dist, "not_built": NOT_BUILT_LINE}
    if why:
        st = {f: None for f in FORMS}
        res.update(void=True, tables={}, statements={})
    else:
        rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
        res["void"] = False
        res["tables"] = {f: tables(rows, f) for f in FORMS if not per[f]}
        st = {f: (statements(res["tables"][f], f) if not per[f] else None) for f in FORMS}
        res["statements"] = {f: v for f, v in st.items() if v is not None}
        res["trend_record"] = {f: {k: rows[k][0]["records"][name(f, ARM)]["trend"] for k in TASK_KEYS}
                               for f in FORMS if not per[f]}
        res["whole_pool"] = {k: rows[k][0]["checks"][HALF].get("whole_pool") for k in TASK_KEYS}
    res["candidates"] = candidates(st, whole_void=bool(why))
    res["follows"] = {c: follows(c, v) for c, v in res["candidates"].items()}
    return res


def _st(s: dict) -> str:
    return h28._st(s)


def _terms(c: dict) -> str:
    return ", ".join(f"{m} {c[m]:+.4f}" for m in (*TERMS, *SHAPE_TERMS))


def statement_lines(res: dict, form: str) -> list[str]:
    tab, st, cand = res["tables"][form], res["statements"][form], CANDIDATE[form]
    pool = POOL_WORDS[form]

    def fold(w, a, b):
        return tab[w]["fold E10.5"][pname(name(form, a), name(form, b))]
    return [
        f"S1 {cand} strength 1 AHEAD of the copy on {pool}: {str(st['S1']).upper()} (fold E10.5; "
        + "; ".join(f"{w} {_st(fold(w, ARM, CONTROL)['score'])}" for w, _ in WORLDS)
        + f"; bar: mean above +{BAR} with at least {MIN_SEEDS} of 10 seeds above 0, in both worlds)",
        f"S2 {cand} strength 1 AHEAD of its placebo: {str(st['S2']).upper()} (fold E10.5; "
        + "; ".join(f"{w} {_st(fold(w, ARM, PLACEBO)['score'])}" for w, _ in WORLDS) + "; the same bar)",
        f"S3 {cand} strength 1 gains on the DE terms: {str(st['S3']).upper()} (fold E10.5, arm minus copy; "
        + "; ".join(f"{w} de_score {fold(w, ARM, CONTROL)['de_score']:+.4f}, de_direction "
                    f"{fold(w, ARM, CONTROL)['de_direction']:+.4f}" for w, _ in WORLDS)
        + "; both above 0 in both worlds)"]


def text(res: dict) -> str:
    lines = [L0.format(commit=res["commit"]), L1, res["not_built"], "J6: " + h28.j6_line(STRENGTH)]
    if res["why"]:
        lines.append("VOID (the whole read): " + "; ".join(res["why"]))
    for f in FORMS:
        if res["form_why"].get(f):
            lines.append(f"VOID ({f}, {CANDIDATE[f]}): " + "; ".join(res["form_why"][f]))
    for f in FORMS:
        g = (res.get("guards") or {}).get("V4", {}).get(f) or {}
        if g and not res["why"]:
            lines.append(f"GATE {f} (the control against row H27's committed {h27_arm(f)}): largest distance "
                         f"{max(v['mean'] for v in g.values()):.3f} (ten-seed mean) and "
                         f"{max(v['seed'] for v in g.values()):.3f} (one seed), board points; void beyond "
                         f"{GUARD_MEAN} and {GUARD_SEED}.")
    for f in FORMS:
        if f in res.get("statements", {}):
            lines += statement_lines(res, f)
    for c, v in res["candidates"].items():
        lines.append(f"OUTCOME {c}: {v['outcome']}")
        lines.append(res["follows"][c])
    if res.get("tables"):
        lines.append("TABLES (reported, deciding nothing; paired by seed, mean [min, max] seeds ahead of 10 (se); "
                     "per-term mean skill differences):")
        for tab in res["tables"].values():
            for w, _ in WORLDS:
                for u in UNITS:
                    for p, c in tab[w][u].items():
                        lines.append(f"{TAG} {w} {u} {p}: {_st(c['score'])}; {_terms(c)}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def cmd_run(a) -> None:
    refuse(Path(a.out) / PART, a.commit)
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS} ({NOT_BUILT_LINE})")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data), a.commit)


def cmd_read(a) -> None:
    refuse_state(a.commit, "read")
    out = Path(a.out)
    if (out / "reading.json").exists() or (out / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")
    res = reading(load_parts(Path(a.runs)), load_h27(Path(a.h27)), a.commit)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "reading.txt").write_text(text(res))
    for c in res["candidates"]:
        print(f"OUTCOME {c}: {res['candidates'][c]['outcome']}")
        print(res["follows"][c])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", required=True, help=f"one of {TASK_KEYS}")
    r.add_argument("--seeds", required=True, help="a range such as 0-1")
    r.add_argument("--out", required=True)
    r.add_argument("--commit", required=True, help="row H29a's run commit (full hash; HEAD must be it)")
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--h27", default=str(H27_PARTS))
    d.add_argument("--out", required=True)
    d.add_argument("--commit", required=True, help="row H29a's run commit (full hash; HEAD must be it)")
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
