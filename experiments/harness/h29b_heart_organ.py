#!/usr/bin/env python
"""Row H29b of experiments/harness/TRIES.md: the practice read of the heart's organ-cut file D (row R21; the
specification of the score-blind round, the cold challenger's final text governing), extending row H27's script
(h27_heart_forecast_size, imported, not copied). A labelled practice analogue that prices the downside and checks the
code; it cannot test D's bet (only the tautological organ truth can favour D over C), decides nothing and builds no
test file.

  gate  target-free and not a try: reads E10.5 only, the test's reference input, with no target and no score. D's pool
        share of E10.5, its core share, L, and the Jaccard overlap of D's rows with C's pool (radial_trim_mask 0) and
        with B's pool (radial_trim_mask 0.5). Stops fixed in the specification: an unmapped label refuses (the config
        is fixed first, in a new commit; exit 4); a Jaccard overlap with C or with B at or above 0.90 means D is not a
        new form: STOP (exit 3), nothing further runs and nothing counts. PASS exits 0. Prints one line.
  run   one transition, a seed range, one JSON under --out; refused before anything loads unless row H29b is in
        TRIES.md, the tree is clean by vec.provenance's rule, HEAD is --commit, the gate of --commit said PASS and the
        output does not exist. Arms, each through the checked builder after row H25's half:
          full_iso_half:f, half_iso_half:f, core_iso_half:f   controls A_f, B_f, C_f, built by row H27's code
                                (`harness._copy` on the form's pool, `h27.size_arm` at the forecast size f);
          organ:f   D_f (counted): `harness._copy` on D's pool (vec.t2_build.organ_scope on the transition's reference
                    stage), rescaled about the centroid to the size procedure's core result on the transition's inputs
                    (row H27's in-memory raw result, z_modes.raw.pred_rms_core) times rms(D pool) / rms(core) on the
                    reference stage (vec.t2_build.organ_core_rms_ratio): the arithmetic of
                    vec.t2_build.size_organ_scope, whose E10.5 checks are redone by row H27's V7;
          organ:p   D_p (counted; an oracle control, never a centre): the same draw at the target's core RMS radius
                    times the same ratio (a perfect forecast of the as-released stage through D's field);
          label:f   the pure-label mask (counted diagnostic): the chamber categories without the tissue step, sized as
                    D_f with its own ratio.
        Truths: row H27's five plus `organ`, the target cut by D's own operator on the target's labels (a tautology,
        as C shares radial_0's operator). A label of the transition's stages or target that the config does not list
        writes a part that names the labels and holds no row (the read is then void). Prints only
        "<task> seed <s> (<sec> s) done" and "wrote <path>".
  read  refused unless row H29b is in TRIES.md, the tree is clean and HEAD is --commit, and with an existing reading;
        applies the void conditions (V3 among them: a truth, control pool or forecast size different from row
        H27's committed parts, or a V7 crosscheck that fails) and the fixed statements; prints only "OUTCOME H29b: READ" or "OUTCOME H29b: VOID"
        and the fixed what-follows line; writes the numbers to reading.json and reading.txt.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from vec import harness, prep, provenance, t2_build
from vec.adoption import SEEDS
from vec.boards import load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parents[2]  # this checkout: the repository and the commit come from here


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h27 = _load("h27_heart_forecast_size", HERE / "experiments" / "harness" / "h27_heart_forecast_size.py")
h26, h25, payoff = h27.h26, h27.h25, h27.payoff

ROW = "H29b"
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h29b-heart-organ"
PART = "h29b.json"
GATE_JSON, GATE_TXT = "gate.json", "gate.txt"
SOURCE = "experiments/harness/h29b_heart_organ.py (row H29b)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h29b_heart_organ"
H27_RESULTS = HERE / "experiments" / "harness" / "results" / "h27_heart_forecast_size"
H27_PARTS = H27_RESULTS / "parts"
H27_READING = H27_RESULTS / "reading.json"
H26_PARTS = h27.H26_PARTS
REFERENCE_FILE = "T2/E10.5.h5ad"  # the gate's one stage: the test's reference input

TASK_KEYS = h27.TASK_KEYS  # row H27's three transitions
OWN, CLEAN = h27.OWN, h27.CLEAN
ORGAN = "organ"
H27_TRUTHS = tuple(h27.TRUTHS)  # own, radial_0.5, radial_0, thin_0.5, thin_0
TRUTHS = (*H27_TRUTHS, ORGAN)
FORMS = tuple(h27.FORMS)  # full_iso_half (A), half_iso_half (B), core_iso_half (C)
LETTER = h27.LETTER
SIZING = "f"
CONTROLS = tuple(h27.arm_name(f, SIZING) for f in FORMS)
D_F, D_P, LABEL_F = "organ:f", "organ:p", "label:f"
COUNTED = (D_F, D_P, LABEL_F)  # three tries per transition (rule 4), as row H27 counted its arms
BUILT = (*CONTROLS, *COUNTED)
ALL = (harness.FLOOR, *BUILT)
POOLS = ("organ", "label")  # D's pool and the pure-label pool
ARM_POOL = {D_F: "organ", D_P: "organ", LABEL_F: "label"}
CORE_FORM = "core_iso_half"  # C: its field (z_modes.raw.pred_rms_core) is D's base

GATE_STOP_JACCARD = 0.90  # fixed in the specification: at or above, D is not a new form
PASS, STOP, REFUSED = "PASS", "STOP", "REFUSED"
GATE_EXIT = {PASS: 0, STOP: 3, REFUSED: 4}
READ, VOID = "READ", "VOID"
SIZE_LOG_TOL = h27.SIZE_LOG_TOL  # 1e-9
TOL = h27.TOL  # 1e-9
FACTOR_RTOL = h27.FACTOR_RTOL  # 1e-12
EXPRESSION_ONLY = h27.EXPRESSION_ONLY
TAG = "[practice analogue, released stages]"

# per world (own, CLEAN, then each truth) and per unit (each transition, the folds, their mean): row H27's
WORLDS = (*h27.WORLDS, (f"truth {ORGAN}", (ORGAN,)))
H27_WORLDS = tuple(w for w, _ in h27.WORLDS)
UNITS = h27.UNITS
UNIT_LABEL = h27.UNIT_LABEL
HEADLINE = "own"  # J4: the least favourable committed world, the truth as released
STATED_PAIRS = ((D_F, "core_iso_half:f"), (D_F, "half_iso_half:f"), (D_F, "full_iso_half:f"))
DIAGNOSTIC_PAIRS = ((LABEL_F, D_F), (D_P, D_F))


def pname(a: str, b: str) -> str:
    return f"{a} - {b}"


# ------------------------------------------------------------------------------------------------ refusals
def tries_rows(path: Path = TRIES) -> set[str]:
    return h27.tries_rows(path)


def git_head_clean() -> tuple[str, bool]:
    head, unclean = provenance.git_state(HERE)
    return head, bool(unclean is not None and not unclean)


def refuse_state(commit: str | None, what: str = "run") -> None:
    """No row H29b, an unclean tree (vec.provenance's rule), or a HEAD other than the run commit."""
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its {what})")
    head, clean = git_head_clean()
    if not clean:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    if not commit or head != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit {commit!r}")


def gate_passed(gate_dir: Path, commit: str | None) -> bool:
    try:
        g = json.loads((Path(gate_dir) / GATE_JSON).read_text())
    except (OSError, ValueError):
        return False
    return bool(g.get("verdict") == PASS and commit and g.get("commit") == commit)


def refuse_run(path: Path, commit: str | None, gate_dir: Path) -> None:
    """The refusals before anything loads."""
    refuse_state(commit)
    if not gate_passed(gate_dir, commit):
        sys.exit(f"{ROW} is refused: the gate of {commit!r} has not said {PASS} ({Path(gate_dir) / GATE_JSON})")
    if Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


def run_provenance() -> dict:
    """The commit and cleanliness of this checkout (from the script's own location, never the working folder)."""
    def git(*a):
        return subprocess.run(["git", "-C", str(HERE), *a], capture_output=True, text=True).stdout.strip()
    return {"commit": git("rev-parse", "HEAD"),
            "tree_clean": git("status", "--porcelain", "--untracked-files=no") == "",
            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())}


# ------------------------------------------------------------------------------------------------ the masks
def label_mask(labels, spec: t2_build.T2ChamberSpec | None = None) -> np.ndarray:
    """The pure-label mask: the chamber categories of configs/t2_chamber_classes.yaml, without the tissue step. An
    unmapped label refuses (UnmappedT2LabelError)."""
    spec = spec or t2_build.load_t2_chamber_spec()
    return np.isin(t2_build.t2_chamber_categories(labels, spec), spec.kept)


def organ_pools(C: np.ndarray, labels, core: np.ndarray, spec=None) -> tuple[dict, dict]:
    """D's pool (vec.t2_build.organ_scope) and the pure-label pool on one stage, as row masks, and D's record."""
    spec = spec or t2_build.load_t2_chamber_spec()
    keep, rec = t2_build.organ_scope(C, labels, core, spec)
    return {"organ": keep, "label": label_mask(labels, spec) | np.asarray(core, bool)}, rec


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    u = int((a | b).sum())
    return float((a & b).sum() / u) if u else 1.0


def gate_verdict(j_c: float, j_b: float) -> str:
    return STOP if (j_c >= GATE_STOP_JACCARD or j_b >= GATE_STOP_JACCARD) else PASS


def gate_figures(C: np.ndarray, labels, core: np.ndarray, spec=None) -> dict:
    """The gate's figures on one stage (E10.5 at the gate): target-free. Raises UnmappedT2LabelError on a label the
    config does not list."""
    C = np.asarray(C, dtype=np.float64)[:, :3]
    core = np.asarray(core, dtype=bool)
    masks, rec = organ_pools(C, labels, core, spec)
    keep = masks["organ"]
    c_pool = t2_build.radial_trim_mask(C, core, 0.0)
    b_pool = t2_build.radial_trim_mask(C, core, 0.5)
    j_c, j_b = jaccard(keep, c_pool), jaccard(keep, b_pool)
    return {"pool_share": float(keep.mean()), "core_share": float(core[keep].mean()), "link_length": rec["link_length"],
            "link_percentile": rec["link_percentile"], "jaccard_C": j_c, "jaccard_B": j_b,
            "label_pool_share": float(masks["label"].mean()), "chamber_share": rec["chamber_share"],
            "C_pool_share": float(c_pool.mean()), "B_pool_share": float(b_pool.mean()),
            "components_kept": rec["components_kept"], "components_total": rec["components_total"],
            "label_config": rec["label_config"], "label_config_sha256": rec["label_config_sha256"],
            "stop_jaccard": GATE_STOP_JACCARD, "verdict": gate_verdict(j_c, j_b)}


GATE_LINES = {
    PASS: "GATE PASS: D is a new form (its overlap with C's pool and with B's pool is below 0.90); the run may follow "
          "at this commit.",
    STOP: "GATE STOP: D is not a new form (its overlap with C's or B's pool is 0.90 or more); nothing further runs and "
          "nothing counts.",
    REFUSED: "GATE REFUSED: a label of E10.5 has no category in configs/t2_chamber_classes.yaml; the config is fixed "
             "first, in a new commit, and nothing runs before.",
}


def gate_text(g: dict) -> str:
    if g["verdict"] == REFUSED:
        return GATE_LINES[REFUSED] + " Unmapped: " + ", ".join(f"{k} ({v})" for k, v in g["unmapped"].items()) + "\n"
    return (f"Row {ROW}'s target-free gate (commit {g['commit']}; E10.5 only, the test's reference input; no target, "
            f"no score; not a try). D's pool share of E10.5 {g['pool_share']:.4f}; its core share "
            f"{g['core_share']:.4f}; "
            f"L {g['link_length']:.4f} (the {g['link_percentile']:g}th percentile of the core's nearest-core "
            f"distances, stored units); Jaccard overlap of D's rows with C's pool {g['jaccard_C']:.4f} and with B's "
            f"pool {g['jaccard_B']:.4f} (stop at {GATE_STOP_JACCARD:.2f}). Beside them: the pure-label pool's share "
            f"{g['label_pool_share']:.4f}, C's {g['C_pool_share']:.4f}, B's {g['B_pool_share']:.4f}; components kept "
            f"{g['components_kept']} of {g['components_total']}; label config sha256 {g['label_config_sha256']}.\n"
            + GATE_LINES[g["verdict"]] + "\n")


def run_gate(out: Path, data_dir: Path, commit: str | None) -> int:
    """The gate step: refusals first, then E10.5 loaded as the builder loads it (duplicated records dropped)."""
    refuse_state(commit, "gate")
    out = Path(out)
    if (out / GATE_JSON).exists() or (out / GATE_TXT).exists():
        sys.exit(f"{ROW} is refused: {out} already holds a gate reading")
    ref = load_stage(Path(data_dir) / REFERENCE_FILE, drop_duplicates=True)
    C = np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3]
    core = np.asarray(ref.obs["core"], dtype=bool)
    labels = ref.obs[t2_build.T2_LABEL_KEY].to_numpy()
    base = {"row": ROW, "stage": REFERENCE_FILE, "commit": commit, "provenance": run_provenance(),
            "label": "target-free gate of row H29b: not a try; E10.5 only; no target and no score"}
    try:
        g = base | gate_figures(C, labels, core)
    except t2_build.UnmappedT2LabelError as e:
        g = base | {"verdict": REFUSED, "unmapped": e.labels}
    out.mkdir(parents=True, exist_ok=True)
    (out / GATE_JSON).write_text(json.dumps(g, indent=1, default=float))
    (out / GATE_TXT).write_text(gate_text(g))
    print(GATE_LINES[g["verdict"]].split(":")[0], flush=True)
    return GATE_EXIT[g["verdict"]]


# ------------------------------------------------------------------------------------------------ the arms
def size_to(C: np.ndarray, target: float, sizing: str, ratio: float, base: float, source: str
            ) -> tuple[np.ndarray, dict]:
    """An isotropic rescale about the centroid to `target` (= base x ratio), through the builder's rescale_to_rms."""
    C = np.asarray(C, dtype=np.float64)
    before = t2_build.rms_radius(C)
    Cs = t2_build.rescale_to_rms(C, float(target))
    after = t2_build.rms_radius(Cs)
    return Cs, {"sizing": sizing, "target_rms_full": float(target), "base_rms": float(base), "ratio": float(ratio),
                "target_source": source, "rms_before": before, "rms_after": after,
                "log_rms_after_vs_target": float(np.log(after / target))}


F_SOURCE = "z_modes.raw.pred_rms_core of row H27's in-memory size result x rms(pool) / rms(core) on the reference stage"
P_SOURCE = "rms(target core) x rms(pool) / rms(core) on the reference stage (perfect size, an oracle control)"


def d_targets(ftargets: dict, ptargets: dict, ratios: dict) -> dict:
    """The sizes of the three counted arms: D_f and the label arm at the core forecast times their own pool's ratio,
    D_p at the target's core RMS radius times D's ratio."""
    fc, pc = float(ftargets[CORE_FORM]), float(ptargets[CORE_FORM])
    return {D_F: {"target": fc * ratios["organ"], "base": fc, "ratio": ratios["organ"], "source": F_SOURCE},
            D_P: {"target": pc * ratios["organ"], "base": pc, "ratio": ratios["organ"], "source": P_SOURCE},
            LABEL_F: {"target": fc * ratios["label"], "base": fc, "ratio": ratios["label"], "source": F_SOURCE}}


def _checked(ctx: harness.BuildContext, X, C, what: str) -> Stage:
    sub, rep = t2_build.build_submission(ctx.board, X, C, genes=list(ctx.board.genes))
    if not rep.ok:
        raise ValueError(f"{ctx.task.key} {what}: the built file fails its check: {rep.errors}")
    return Stage.from_anndata(sub, ctx.board.genes, need_coords=True)


def build_control(form: str, ctx: harness.BuildContext, factors: dict, size: dict, ptargets: dict
                  ) -> tuple[Stage, dict]:
    """Row H27's arm `form:f` on the context whose reference input is the form's pool: one harness draw, row H27's
    size_arm at the forecast size (with its half), the checked builder; bit for bit h27.build_forms's f."""
    st = harness._copy(ctx)
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    X, C, rec = h27.size_arm(st.X, st.coords, form, SIZING, factors, size, ptargets, seed)
    return _checked(ctx, X, C, h27.arm_name(form, SIZING)), rec


def build_pool_arms(pool: str, ctx: harness.BuildContext, targets: dict) -> tuple[dict, dict, dict]:
    """The counted arms of one pool on the context whose reference input is that pool: one harness draw, the size step
    of each arm, row H25's half, the checked builder. Returns the stages, their size records and the draw record
    (the practice count drawn against the pool's rows: rows repeat with jitter when the count exceeds the pool)."""
    st = harness._copy(ctx)
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    stages, recs = {}, {}
    for arm in (a for a in COUNTED if ARM_POOL[a] == pool):
        t = targets[arm]
        Cs, rec = size_to(st.coords, t["target"], arm.split(":")[1], t["ratio"], t["base"], t["source"])
        rows = h25.kept_rows(len(Cs), seed, 1)
        stages[arm] = _checked(ctx, np.asarray(st.X)[rows], Cs[rows], arm)
        recs[arm] = rec
    n_pool = int(ctx.inputs["ref"].n_obs)
    draw = {"pool_cells": n_pool, "drawn": int(st.n), "short": bool(st.n > n_pool),
            "repeated_rows_at_least": int(max(0, st.n - n_pool))}
    return stages, recs, draw


# ------------------------------------------------------------------------------------------------ the run
def _unmapped_part(t, where: str, e: t2_build.UnmappedT2LabelError, seeds) -> dict:
    return {"label": f"row {ROW}: a label without a T2 chamber category; no arm built", "axis": AXIS, "task": t.key,
            "unmapped": {"where": where, "labels": e.labels}, "seeds": list(seeds), "rows": [], "source": SOURCE,
            "provenance": run_provenance()}


def run_transition(key: str, seeds, out: Path, data_dir: Path, commit: str | None, gate_dir: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse_run(path, commit, gate_dir)
    if key not in TASK_KEYS:
        sys.exit(f"{key}: not a transition of row {ROW}")
    t = payoff.heart_task(key)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    spec = t2_build.load_t2_chamber_spec()
    # the forecast first, from the input stages only (row H27's)
    geoms, _ = h27.input_geometries(t, data_dir)
    target_t = float(harness.STAGE_DAYS[t.target])
    size = h27.analogue_size_result(geoms, target_t)
    ftargets = h27.forecast_targets(size)
    # the truths, as row H27 loads them, and the organ truth by D's own operator on the target's labels
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    tm = h26.truth_masks(TC, tcore)
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    whole = inputs["ref"]
    pcore = np.asarray(whole.obs["core"], dtype=bool)
    PC = np.asarray(whole.obsm[COORD_KEY], dtype=float)[:, :3]
    try:
        tmask, trec = organ_pools(TC, ta.obs[t2_build.T2_LABEL_KEY].to_numpy(), tcore, spec)
    except t2_build.UnmappedT2LabelError as e:
        return _write(path, _unmapped_part(t, "target", e, seeds))
    try:
        pm, prec = organ_pools(PC, whole.obs[t2_build.T2_LABEL_KEY].to_numpy(), pcore, spec)
    except t2_build.UnmappedT2LabelError as e:
        return _write(path, _unmapped_part(t, "reference", e, seeds))
    tm[ORGAN] = tmask["organ"]
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    truths = {q: Stage.from_anndata(ta[tm[q]], genes, need_coords=True) for q in TRUTHS}
    t_info = {q: h26.truth_info(TC, tcore, tm[q]) for q in TRUTHS}
    # the controls' pools, factors and perfect targets, as row H27 builds them
    fm = {f: m for f, m in h26.pool_masks(PC, pcore).items() if f in FORMS}
    pools = {f: (whole if h26.TRIM[f] is None else whole[m].copy()) for f, m in fm.items()}
    pools.update({p: whole[pm[p]].copy() for p in POOLS})
    p_info = {f: h26.pool_info(PC, pcore, m) for f, m in fm.items()}
    p_info.update({p: h26.pool_info(PC, pcore, pm[p]) for p in POOLS})
    factors = h26.oracle_factors(TC, tcore, PC, pcore)
    ptargets = h27.perfect_targets(TC, tcore)
    ratios = {p: t2_build.organ_core_rms_ratio(PC, pm[p], pcore) for p in POOLS}
    targets = d_targets(ftargets, ptargets, ratios)
    v7 = h27.v7_crosscheck(PC, size)
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        preds, recs, draws = {}, {}, {}

        def ctx_of(pool, s=s):
            return harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pools[pool]},
                                        seed=harness.BUILD_SEED_OFFSET + s)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            for form in FORMS:
                st, rec = build_control(form, ctx_of(form), factors, size, ptargets)
                preds[h27.arm_name(form, SIZING)], recs[h27.arm_name(form, SIZING)] = st, rec
            for pool in POOLS:
                stages, rs, draws[pool] = build_pool_arms(pool, ctx_of(pool), targets)
                preds.update(stages)
                recs.update(rs)
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
               "n_cells": {a: int(p.n) for a, p in preds.items()},
               "size_records": recs, "draws": draws, "scores": {}, "skills": {}, "groups": {}}
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
    res = {"label": "row H29b: the organ cut D (forecast size f counted, perfect size p an oracle control, the pure "
                    "label mask a counted diagnostic) beside row H27's A, B and C at the forecast size; a practice "
                    "analogue, not an adoption result", "axis": AXIS, "task": t.key, "arms": list(BUILT),
           "truth_values": list(TRUTHS), "unmapped": None, "size_result": size, "forecast_targets": ftargets,
           "perfect_targets": ptargets, "oracle_factors": factors, "ratios": ratios, "d_targets": targets,
           "organ_record": {"reference": prec, "target": trec}, "v7": v7, "truth_info": t_info, "pred_info": p_info,
           "truth_seed": payoff.TRUTH_SEED, "half_seed_base": h25.HALF_SEED_BASE, "departures": h27.DEPARTURES,
           "seeds": list(seeds), "rows": rows_out, "source": SOURCE, "provenance": run_provenance()}
    return _write(path, res)


def _write(path: Path, res: dict) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def load_h27(root: Path = H27_PARTS) -> dict[str, list[dict]]:
    return h27.load_parts(root)


def load_h27_reading(path: Path = H27_READING) -> dict:
    return json.loads(Path(path).read_text())


def load_gate(gate_dir: Path) -> dict | None:
    try:
        return json.loads((Path(gate_dir) / GATE_JSON).read_text())
    except (OSError, ValueError):
        return None


# ------------------------------------------------------------------------------------------------ the guard (J7)
def guard_scores(key: str, rows: list[dict], theirs_rows: list[dict], pairs, truths=H27_TRUTHS) -> list[str]:
    """Bit for bit: each (mine, theirs) arm's score equals the committed one in every seed and truth, exactly."""
    other = {r["seed"]: r for r in theirs_rows}
    why = []
    for mine, theirs in pairs:
        for q in truths:
            for r in rows:
                o = other.get(r["seed"])
                try:
                    a, b = r["scores"][q][mine], o["scores"][q][theirs]
                except (KeyError, TypeError):
                    why.append(f"J7: {key} {q} seed {r['seed']}: no committed score of {theirs} to compare")
                    continue
                if a != b:
                    why.append(f"J7: {key} {q} seed {r['seed']}: {mine} is {a - b:+.3e} from the committed {theirs}")
    return why


def guard_reading(rows: dict[str, list[dict]], committed: dict) -> list[str]:
    """Bit for bit: the controls' plain means per world and unit, by row H27's own arithmetic (h27._series), equal row
    H27's committed reading.json (tables[world][unit].mean.f)."""
    why = []
    if committed.get("void") is not False:
        return ["J7: row H27's committed reading is void or unreadable"]
    for w, ts in h27.WORLDS:
        for u, fs in UNITS.items():
            want = committed["tables"][w][u]["mean"][SIZING]
            for f in FORMS:
                got = float(h27._series(rows, fs, ts, f, SIZING).mean())
                if got != want[LETTER[f]]:
                    why.append(f"J7: {w} {u} {LETTER[f]}_f mean {got!r} is not row H27's {want[LETTER[f]]!r}")
    return why


# ------------------------------------------------------------------------------------------------ the void
def void_reasons(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], h27_reading: dict,
                 gate: dict | None, commit: str | None) -> tuple[list[str], dict]:
    """Why the read is void (empty when it may be read), and the pools-short report (never a void)."""
    why: list[str] = []
    short: dict = {}
    if gate is None:  # G
        why.append("G: no gate reading")
    elif gate.get("verdict") != PASS:
        why.append(f"G: the gate said {gate.get('verdict')!r}, not {PASS}")
    elif gate.get("commit") != commit:
        why.append("G: the gate reading is not of the run commit")
    every = [p for ps in parts.values() for p in ps]
    complete = True
    for key in TASK_KEYS:  # V1
        ps = parts.get(key)
        if not ps:
            why.append(f"V1: missing transition {key}")
            complete = False
            continue
        for p in ps:
            if p.get("unmapped"):
                u = p["unmapped"]
                why.append(f"V1: {key}: labels without a T2 chamber category on the {u['where']} stage: "
                           f"{sorted(u['labels'])}")
                complete = False
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
            if set(r.get("size_records", {})) != set(BUILT) or set(r.get("n_cells", {})) != set(BUILT) \
                    or set(r.get("draws", {})) != set(POOLS):
                why.append(f"V1: {key} seed {r['seed']}: size records, cell numbers or draws missing")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"V1: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not every:
        return why or ["V1: no part"], short
    if len({p["provenance"]["commit"] for p in every}) > 1:  # V2
        why.append("V2: the parts come from different commits")
    if not all(p["provenance"]["tree_clean"] for p in every):
        why.append("V2: a part comes from a dirty tree")
    if not commit or any(str(p["provenance"]["commit"]) != commit for p in every):
        why.append(f"V2: a part does not come from the row's run commit {commit!r}")
    if not complete:
        return why, short
    merged = {}
    for key in TASK_KEYS:
        ps, theirs = parts[key], h27_parts.get(key)
        if not theirs:
            why.append(f"V3: {key}: no committed row H27 part")
            continue
        ref = theirs[0]
        for p in ps:  # V3
            for q in H27_TRUTHS:
                if not h27._same(p["truth_info"].get(q), ref["truth_info"].get(q), "cells"):
                    why.append(f"V3: {key}: truth {q} differs from row H27's (cells or core share)")
            for f in FORMS:
                if not h27._same(p["pred_info"].get(f), ref["pred_info"].get(f), "pool_cells"):
                    why.append(f"V3: {key}: {f}'s pool differs from row H27's (cells or core share)")
                a, b = p["forecast_targets"][f], ref["forecast_targets"][f]
                if not abs(a - b) <= FACTOR_RTOL * abs(b):
                    why.append(f"V3: {key}: {f}'s forecast size differs from row H27's")
            if not p["v7"]["ok"]:
                why.append(f"V3: {key}: the reference cloud is not the cloud the size result measured")
            for field in ("forecast_targets", "perfect_targets", "ratios", "d_targets"):
                if p[field] != ps[0][field]:
                    why.append(f"V3: {key}: the parts disagree on the {field.replace('_', ' ')}")
        rows = h26.merge_rows(ps)[0]
        merged[key] = rows
        why += guard_scores(key, rows, h26.merge_rows(theirs)[0], tuple(zip(CONTROLS, CONTROLS, strict=True)))
        why += size_checks(key, rows, ps[0]["d_targets"])
        short[key] = {pool: {"pool_cells": rows[0]["draws"][pool]["pool_cells"],
                             "drawn": sorted({r["draws"][pool]["drawn"] for r in rows}),
                             "short_in_seeds": sum(bool(r["draws"][pool]["short"]) for r in rows),
                             "repeated_rows_at_least": max(r["draws"][pool]["repeated_rows_at_least"] for r in rows)}
                      for pool in POOLS}
    if not why:
        why += guard_reading(merged, h27_reading)
    return why, short


def size_checks(key: str, rows: list[dict], targets: dict) -> list[str]:
    """Each counted arm at its target (log within 1e-9); D_p and D_f one draw: the same cell number and the same
    expression-only skills (1e-9) in every truth."""
    why = []
    for r in rows:
        for arm in COUNTED:
            rec = r["size_records"][arm]
            if not abs(rec["log_rms_after_vs_target"]) <= SIZE_LOG_TOL:
                why.append(f"V5: {key} seed {r['seed']}: {arm} is not at its target")
            if not abs(rec["target_rms_full"] - targets[arm]["target"]) <= SIZE_LOG_TOL * abs(targets[arm]["target"]):
                why.append(f"V5: {key} seed {r['seed']}: {arm}'s target is not its field")
        if r["n_cells"][D_P] != r["n_cells"][D_F]:
            why.append(f"V5: {key} seed {r['seed']}: {D_P} has another cell number than {D_F}")
        for q in TRUTHS:
            e = max(abs(r["skills"][q][D_P][m] - r["skills"][q][D_F][m]) for m in EXPRESSION_ONLY)
            if e > TOL:
                why.append(f"V5: {key} seed {r['seed']} {q}: {D_P}'s expression skills differ from {D_F}'s")
    return why


# ------------------------------------------------------------------------------------------------ the tables
def _stats(d) -> dict:
    d = np.asarray(d, float)
    return {"mean": float(d.mean()), "se": float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else float("nan"),
            "lo": float(d.min()), "hi": float(d.max()), "ahead": int((d > 0).sum())}


def series(rows: dict[str, list[dict]], folds, truths, arm: str) -> np.ndarray:
    return h26.series(rows, folds, truths, arm)


def cell(rows: dict[str, list[dict]], folds, truths) -> dict:
    s = {a: series(rows, folds, truths, a) for a in BUILT}
    return {"mean": {a: float(v.mean()) for a, v in s.items()},
            "stated": {pname(a, b): _stats(s[a] - s[b]) for a, b in STATED_PAIRS},
            "diagnostic": {pname(a, b): _stats(s[a] - s[b]) for a, b in DIAGNOSTIC_PAIRS}}


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: cell(rows, fs, ts) for u, fs in UNITS.items()} for w, ts in WORLDS}


def reading(parts: dict[str, list[dict]], h27_parts: dict[str, list[dict]], h27_reading: dict, gate: dict | None,
            commit: str | None) -> dict:
    try:
        why, short = void_reasons(parts, h27_parts, h27_reading, gate, commit)
    except (KeyError, TypeError, ValueError, IndexError) as e:  # a malformed part
        why, short = [f"a part cannot be read ({type(e).__name__}: {e})"], {}
    res = {"row": ROW, "commit": commit, "gate": gate, "why": why, "void": bool(why),
           "outcome": VOID if why else READ, "decides": DECIDES}
    if not why:
        rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
        res.update(tables=tables(rows), headline_world=HEADLINE, pools_short=short,
                   organ_records={k: parts[k][0]["organ_record"] for k in TASK_KEYS},
                   ratios={k: parts[k][0]["ratios"] for k in TASK_KEYS},
                   d_targets={k: parts[k][0]["d_targets"] for k in TASK_KEYS},
                   pool_info={k: parts[k][0]["pred_info"] for k in TASK_KEYS},
                   truth_info={k: parts[k][0]["truth_info"] for k in TASK_KEYS})
    res["follows"] = FOLLOWS[res["outcome"]]
    return res


# ------------------------------------------------------------------------------------------------ the fixed words
DECIDES = ("nothing: the read prices the downside and checks the code; it cannot test D's bet (only the tautological "
           "organ truth can favour D over C); no pick rule and no world weight; both readings go to the captain")
FOLLOWS = {
    READ: ("H29b READ: the figures go to fresh analysts, not the main session, who write the captain's question; the "
           "read builds nothing, and file D enters the heart's order only if its two DGX builds give one sha256 and "
           "pass the format check and the captain says yes to it by name."),
    VOID: ("H29b VOID (reasons in reading.txt): no figure exists under row H29b and there is no rerun under it; by the "
           "heart's order rule file D does not qualify and is skipped."),
}
L0 = ("Row H29b (commit {commit}). A practice analogue on released stages (rule 9, J1a); no figure is to be read "
      "beside a test score. It cannot test D's bet: only the organ truth can favour D over C, and that truth is a "
      "tautology (D is built by the operator that cuts it, as C shares radial_0's). It prices the downside and checks "
      "the code. Nothing is built and no file, rule or default changes. The main session states no preference.")
L1 = ("Arms: D_f (organ:f, counted) is the copy draw from D's pool on the transition's reference stage at the size "
      "procedure's core result times rms(D pool)/rms(core) there; D_p (organ:p, counted, an oracle control, never a "
      "centre) is the same draw at the target's core RMS radius times that ratio; the pure-label mask (label:f, a "
      "counted diagnostic) drops the tissue step. Controls A_f, B_f, C_f are row H27's arms at the forecast size, "
      "rebuilt by its code and checked against its committed reading bit for bit (J7).")
HEAD_LINE = ("HEADLINE (J4): the least favourable committed world, the truth as released ({w}); no pick rule and no "
             "world weight; both readings (own and CLEAN) go to the captain.")
J5_LINE = ("J5: a diagnostic: the two E10.5 transitions share one truth; fold E8.5 is the copy_last analogue with no "
           "growth step.")
CAVEAT = ("Rule 2 caveat: the designer saw E10.5's label composition, the truth of two practice transitions; the mask "
          "comes from the 2026-10-19 rule text and the lineage map, not from those shares, but E10.5's cells are not "
          "blind.")


def _st(s: dict) -> str:
    return f"{s['mean']:+.2f} (se {s['se']:.2f}) [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10"


def text(res: dict) -> str:
    if res["void"]:
        return f"OUTCOME {ROW}: {VOID}\n" + res["follows"] + "\nVOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"OUTCOME {ROW}: {READ}", res["follows"], L0.format(commit=res["commit"]), L1,
             HEAD_LINE.format(w=HEADLINE), J5_LINE, CAVEAT]
    if res.get("gate"):
        lines.append(gate_text(res["gate"]).strip())
    for k in TASK_KEYS:
        r, sh = res["ratios"][k], res["pools_short"][k]
        lines.append(f"{k}: ratio rms(pool)/rms(core) on the reference stage: D {r['organ']:.4f}, label "
                     f"{r['label']:.4f}; "
                     + "; ".join(f"{p} pool {v['pool_cells']} rows, drawn {v['drawn']}, short in {v['short_in_seeds']} "
                                 f"of 10 seeds (rows repeated with jitter, at least {v['repeated_rows_at_least']})"
                                 for p, v in sh.items()))
    tab = res["tables"]
    lines.append("STATED, per world and unit, paired by seed: D_f minus C_f, B_f and A_f, mean (se) [min, max] seeds "
                 "ahead of 10; then the diagnostics label_f minus D_f (what the tissue step does) and D_p minus D_f "
                 "(the forecast's size error):")
    for w, _ in WORLDS:
        for u in UNITS:
            c = tab[w][u]
            head = "HEADLINE " if w == HEADLINE else ""
            lines.append(f"{head}{TAG} {w} {UNIT_LABEL[u]}: "
                         + "; ".join(f"{d} {_st(s)}" for d, s in c["stated"].items()) + " | "
                         + "; ".join(f"{d} {_st(s)}" for d, s in c["diagnostic"].items()))
    lines.append("No pick and no advice; these figures go to fresh analysts, who write the captain's question.")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def cmd_gate(a) -> None:
    sys.exit(run_gate(Path(a.out), Path(a.data), a.commit))


def cmd_run(a) -> None:
    refuse_run(Path(a.out) / PART, a.commit, Path(a.gate))
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS}")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data), a.commit, Path(a.gate))


def cmd_read(a) -> None:
    refuse_state(a.commit, "read")
    out = Path(a.out)
    if (out / "reading.json").exists() or (out / "reading.txt").exists():
        sys.exit(f"{out} already holds a reading: there is one read under row {ROW}")
    res = reading(load_parts(Path(a.runs)), load_h27(Path(a.h27)), load_h27_reading(Path(a.h27_reading)),
                  load_gate(Path(a.gate)), a.commit)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "reading.txt").write_text(text(res))
    print(f"OUTCOME {ROW}: {res['outcome']}")
    print(res["follows"])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--commit", required=True, help="row H29b's run commit (full hash; HEAD must be it)")
    g.add_argument("--data", default=str(harness.DATA_DIR))
    r = sub.add_parser("run")
    r.add_argument("--task", required=True, help=f"one of {TASK_KEYS}")
    r.add_argument("--seeds", required=True, help="a range such as 0-1")
    r.add_argument("--out", required=True)
    r.add_argument("--gate", required=True, help="the gate step's folder (its gate.json must say PASS at --commit)")
    r.add_argument("--commit", required=True, help="row H29b's run commit (full hash; HEAD must be it)")
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--gate", required=True)
    d.add_argument("--h27", default=str(H27_PARTS))
    d.add_argument("--h27-reading", default=str(H27_READING))
    d.add_argument("--out", required=True)
    d.add_argument("--commit", required=True, help="row H29b's run commit (full hash; HEAD must be it)")
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
