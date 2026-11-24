#!/usr/bin/env python
"""Row H34 of experiments/harness/TRIES.md: the practice read that the heart's fourth search line requires before any
build (the line's final text: "Setting", "Practice read", "Target-free gate first", "Pass rule", "Void"): row H28's
trend arms on row H29b's D_f form. A screen, not evidence of development. It builds no test file and changes no file,
rule or default. No row on record had read the trend step on D's pool at any strength before this one.

Rows H28's and H29b's scripts are imported, not copied (and through them rows H27's and H26's): from
h28_heart_trend_strength the expression step and its placebo (trend_expression with its fixed derangement), the
series, the paired comparison, the S1 bar, the seeds and the worlds; from h29b_heart_organ D's pool (organ_pools),
D_f's size (d_targets, size_to), row H25's half and the checked builder.

The form is row H29b's D_f: `vec.harness._copy` on a context whose reference input is D's pool (vec.t2_build
.organ_scope on the transition's reference stage, E9.5 or E8.75, the analogue of E10.5 at the test; earlier stages
whole, as the heart builder loads them), rescaled about the centroid to the size procedure's core result on the
transition's input stages (row H27's in-memory raw result) times rms(D pool) / rms(core) on the reference stage, then
row H25's half, then the checked builder. Every arm is that control with ONE change, the expression step, made before
the size step, exactly as row H28 places it on form C:
  control        the copy D_f (a reference, not counted)
  trend_s1.0     row H28's step at strength 1: b0_trend's expression on D's context (`h28.trend_expression(ctx, 1)`;
                 the harness's own b0_trend bit for bit, checked in every seed)
  trend_s1.5     the same step with every lineage's log factor times 1.5: the trend mapping and the damping that
                 momentum_shift gives a lineage the mapping leaves out, both times 1.5
  placebo_s1.0, placebo_s1.5   the same factors at that strength with each shifted lineage's genes moved by row H28's
                 fixed derangement within ten level bins (seed 20261120)
The step's reference is D's pool: it is the context's reference input, so every lineage's last mean is taken on D's
cells. Counted under rule 4: eight tries (two strengths and two placebos on two transitions).

Transitions: the two heart backtests on E10.5 (`T2h/<=E9.5->E10.5`, `T2h/<=E8.75->E10.5`), merged seed by seed into
fold E10.5; seeds 0 to 9; worlds own (the truth as released) and CLEAN (radial_0, thin_0.5, thin_0, equal weights).

THE RULE (row H34's eighth cell), on fold E10.5, every difference paired by seed; a difference of exactly 0 is neither
above nor below 0. Clause (a) holds in a world if trend_s1.5 minus trend_s1.0 has a paired mean below 0 with at least
8 of the 10 seeds below 0 there. Row H28's S1 bar is met in a world if trend_s1.5 minus the control has a paired mean
above +0.3 with at least 8 of the 10 seeds above 0 there. Clause (b) holds unless the bar is met in own and also in
CLEAN. STOPPED if clause (a) holds in own or in CLEAN, or clause (b) holds; NOT STOPPED otherwise. VOID is decided
first, and under VOID the rule is not read. Arm minus placebo and the DE terms are printed as labelled diagnostics and
gate nothing.

  gate  first, counting nothing. No counted arm and no placebo is built or scored in it.
        G1: the control D_f, built and scored on both transitions and seeds 0 to 9 against the four truths, against
        row H29b's committed D_f. Row H29b committed no per-seed part: its committed D_f is the plain ten-seed means
        under `organ:f` in its reading.json (the file of sha256 H29B_READING_SHA256), tables[world][unit].mean for
        the worlds own, CLEAN and the four single truths and the units of the two transitions and fold E10.5:
        eighteen entries. Each, recomputed from the control's per-seed scores by row H29b's own arithmetic
        (h29b.series with the world's truths of h29b.WORLDS and the unit's entry of h29b.UNITS, then the mean over
        the seeds), must equal the committed number exactly, with no tolerance.
        G2: in every seed's context of each transition, the log-factor mapping and the damping that row H28's step
        gives vec.t1_scope.momentum_shift at strength 1.5 against those it gives at strength 1.0 (read through
        h28.trend_expression itself, stopped at that call, so that no count is touched): every lineage's vector
        equals 1.5 times the other, element by element, and the damping equals 1.5 times the other; exactly.
        Writes gate.json and gate.txt under --out; exit 0 PASS, 3 FAIL. A failed gate voids the read.
  run   one transition, a seed range, one JSON under --out: the four counted arms, beside a rebuilt control (the
        paired reference of the read). Refused unless the gate record of this commit says PASS. Prints only
        "<task> seed <s> (<sec> s) done" and "wrote <path>".
  read  the void conditions first (the eighth cell's v3 to v5, as far as code can hold them), then the rule; writes
        reading.json and reading.txt under --out and never overwrites a reading; prints the outcome line and what
        follows, and nothing else. It repeats G1 on the run's own control, the control its differences are paired
        with (void under v3 if it fails).
Every step is refused before anything loads unless row H34 is in TRIES.md, the tracked tree is clean, HEAD is the
commit given and vec is this checkout's own.

Beside the gate and the read, gating nothing: the transitions' target-free records (the truths' and D's pool's cell
numbers and core shares, D's record on the reference stage, both rms ratios, D_f's size, the draw, the size after
the size step) are held against row H29b's committed ones, and any difference is written as a note.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import importlib.util
import json
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2]  # this checkout
sys.path.insert(0, str(HERE / "src"))  # this checkout's vec first: a shared environment may name another checkout

from vec import harness, prep, t1_scope, t2_build  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.core import counts_of  # noqa: E402
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample  # noqa: E402
from vec.io import COORD_KEY, read_stage  # noqa: E402
from vec.prep import load_stage  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


h28 = _load("h28_heart_trend_strength", HERE / "experiments" / "harness" / "h28_heart_trend_strength.py")
h29b = _load("h29b_heart_organ", HERE / "experiments" / "harness" / "h29b_heart_organ.py")
h27, h26, h25, payoff = h29b.h27, h29b.h26, h29b.h25, h29b.payoff  # the form is row H29b's: its own imports

ROW = "H34"
TRIES = HERE / "experiments" / "harness" / "TRIES.md"
AXIS = "h34-heart-organ-trend"
GATE_AXIS = AXIS + ":gate"
PART = "h34.json"
GATE_JSON, GATE_TXT = "gate.json", "gate.txt"
SOURCE = "experiments/harness/h34_heart_organ_trend.py (row H34)"
RESULTS = HERE / "experiments" / "harness" / "results" / "h34_heart_organ_trend"
H29B_READING = h29b.RESULTS / "reading.json"  # row H29b's committed reading (its parts were not committed)
# the file row H29b's outcome cell records (and row H34's eighth cell names): G1 is held against this file only
H29B_READING_SHA256 = "e86754f7b4ca012886aee9d7fcadf969eda26d57f27eeaf928e469d8788a3784"

TASK_KEYS = h28.TASK_KEYS  # the two heart backtests on E10.5
UNITS = h28.UNITS  # each transition, and fold E10.5 (the two merged seed by seed)
FOLD = "fold E10.5"
SEEDS = h28.SEEDS  # 0 to 9
OWN, CLEAN, TRUTHS, WORLDS = h28.OWN, h28.CLEAN, h28.TRUTHS, h28.WORLDS  # worlds own and CLEAN
POOL = "organ"  # D's pool (h29b.organ_pools)
D_F = h29b.D_F  # "organ:f": row H29b's arm, this row's control
LO, HI = 1.0, 1.5  # row H28's strength 1 and the fourth line's 1.5
STRENGTHS = (LO, HI)
CONTROL = h28.CONTROL
PERM_SEED = h28.PERM_SEED  # the placebo's derangement (20261120)
PRACTICE_BUDGET = h28.PRACTICE_BUDGET
BAR, MIN_SEEDS = h28.BAR, h28.MIN_SEEDS  # row H28's S1 bar: +0.3 board points, 8 of 10 seeds
TOL, FACTOR_RTOL = h28.TOL, h28.FACTOR_RTOL  # 1e-9 on a shape term's skill; relative 1e-12 on a budget
SIZE_LOG_TOL = h29b.SIZE_LOG_TOL
DE_TERMS, TERMS, SHAPE_TERMS = h28.DE_TERMS, h28.TERMS, h28.SHAPE_TERMS
TAG = h28.TAG
FILE, RULING_ROW = "OTS15", "R32"  # the line's one file and the ruling row its builder option stands under
PASS, FAIL = "PASS", "FAIL"
GATE_EXIT = {PASS: 0, FAIL: 3}
VOID, STOPPED, NOT_STOPPED = "VOID", "STOPPED", "NOT STOPPED"


def arm(s: float) -> str:
    return f"trend_s{s:.1f}"


def placebo(s: float) -> str:
    return f"placebo_s{s:.1f}"


COUNTED = tuple(n for s in STRENGTHS for n in (arm(s), placebo(s)))  # per transition: eight tries on the two
BUILT = (CONTROL, *COUNTED)
ALL = (harness.FLOOR, *BUILT)
pname = h28.pname
RULE_PAIRS = ((arm(HI), arm(LO)), (arm(HI), CONTROL))  # clause (a), clause (b)
DIAGNOSTIC_PAIRS = (*((arm(s), placebo(s)) for s in STRENGTHS), (arm(LO), CONTROL),
                    *((placebo(s), CONTROL) for s in STRENGTHS))
PAIRS = (*RULE_PAIRS, *DIAGNOSTIC_PAIRS)
# G1: the worlds and units of row H29b's committed table that this row's truths and transitions can form
G1_WORLDS = tuple((w, ts) for w, ts in h29b.WORLDS if set(ts) <= set(TRUTHS))
G1_UNITS = (*TASK_KEYS, FOLD)
ORGAN_FIELDS = ("label_config_sha256", "kept_categories", "link_percentile", "link_length", "kept_share",
                "chamber_share", "core_share_of_kept", "components_kept", "components_total")
TARGET_FIELDS = ("target", "base", "ratio")


# ------------------------------------------------------------------------------------------------ refusals
def tries_rows(path: Path | None = None) -> set[str]:
    return h27.tries_rows(TRIES if path is None else path)


def git_state() -> dict:
    """This checkout's HEAD and whether its tracked files are as committed (from the script's own location)."""
    return h29b.run_provenance()


def vec_home() -> Path:
    """The checkout the imported vec package lives in."""
    return Path(harness.__file__).resolve().parents[2]


def refuse_state(commit: str | None, what: str = "run") -> None:
    """The refusals of every step, before anything loads: no row H34, vec from another checkout, tracked files that
    differ from the commit, or a HEAD other than the commit given."""
    if ROW not in tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its {what})")
    if vec_home() != HERE:
        sys.exit(f"{ROW} is refused: vec is imported from {vec_home()}, not from this checkout")
    state = git_state()
    if state.get("tree_clean") is not True:
        sys.exit(f"{ROW} is refused: tracked files differ from the commit")
    if not commit or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: HEAD is not the run commit {commit!r}")


def load_gate(gate_dir: Path) -> dict | None:
    try:
        return json.loads((Path(gate_dir) / GATE_JSON).read_text())
    except (OSError, ValueError):
        return None


def gate_record_reasons(gate: dict | None, commit: str | None) -> list[str]:
    """v3: the gate failed, or its record is missing or of another commit than the run commit (a record whose tree was
    not the commit's is not of the commit)."""
    if not isinstance(gate, dict) or gate.get("axis") != GATE_AXIS:
        return ["v3: no gate record"]
    if gate.get("verdict") != PASS:
        return [f"v3: the gate said {gate.get('verdict')!r}, not {PASS}"]
    prov = gate.get("provenance") or {}
    if not commit or gate.get("commit") != commit or prov.get("commit") != commit or prov.get("tree_clean") is not True:
        return ["v3: the gate record is not of the run commit"]
    return []


def gate_passed(gate_dir: Path, commit: str | None) -> bool:
    return not gate_record_reasons(load_gate(gate_dir), commit)


def refuse_run(path: Path, commit: str | None, gate_dir: Path) -> None:
    refuse_state(commit)
    if not gate_passed(gate_dir, commit):
        sys.exit(f"{ROW} is refused: the gate of {commit!r} has not said {PASS} ({Path(gate_dir) / GATE_JSON})")
    if Path(path).exists():
        sys.exit(f"{ROW} is refused: {path} exists")


# ------------------------------------------------------------------------------------------------ one transition
@dataclass
class Transition:
    """One backtest as rows H28 and H29b load it, with D's pool as the context's reference input."""

    key: str
    view: harness.TaskView
    cboard: object
    base_board: object
    inputs: dict  # the input stages by role; the reference role holds D's pool, earlier stages are whole
    ref_stage: Stage  # the reference stage as released: the scorer's reference and the floor's source
    truths: dict  # worlds own and CLEAN: own, radial_0, thin_0.5, thin_0
    target: dict  # D_f's size (h29b.d_targets): target, base, ratio, source
    facts: dict  # the target-free records a part and the gate record carry


def load_transition(key: str, data_dir: Path) -> Transition:
    """One backtest, loaded as row H29b loads it: the forecast first, from the input stages only; the target's cells
    enter the truths alone; D's pool by h29b.organ_pools on the reference stage, earlier stages whole. The truths are
    row H28's (own and row H26's CLEAN world)."""
    data_dir = Path(data_dir)
    t = payoff.heart_task(key)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    spec = t2_build.load_t2_chamber_spec()
    geoms, _ = h27.input_geometries(t, data_dir)
    size = h27.analogue_size_result(geoms, float(harness.STAGE_DAYS[t.target]))
    ftargets = h27.forecast_targets(size)
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    tm = h26.truth_masks(TC, tcore)
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    whole = inputs["ref"]
    pcore = np.asarray(whole.obs["core"], dtype=bool)
    PC = np.asarray(whole.obsm[COORD_KEY], dtype=float)[:, :3]
    pm, prec = h29b.organ_pools(PC, whole.obs[t2_build.T2_LABEL_KEY].to_numpy(), pcore, spec)
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    truths = {q: Stage.from_anndata(ta[tm[q]], genes, need_coords=True) for q in TRUTHS}
    ratios = {p: t2_build.organ_core_rms_ratio(PC, pm[p], pcore) for p in h29b.POOLS}
    # D_f's size by row H29b's own arithmetic. Its oracle control (the perfect size) is not formed here: no target
    # enters a size.
    target = h29b.d_targets(ftargets, {h29b.CORE_FORM: float("nan")}, ratios)[D_F]
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    facts = {"truth_info": {q: h26.truth_info(TC, tcore, tm[q]) for q in TRUTHS},
             "pred_info": {POOL: h26.pool_info(PC, pcore, pm[POOL])}, "organ_record": {"reference": prec},
             "ratios": ratios, "d_targets": {D_F: target}, "forecast_targets": ftargets, "size_result": size,
             "v7": h27.v7_crosscheck(PC, size)}
    return Transition(key=t.key, view=harness.TaskView.of(t), cboard=harness._candidate_board(base_board),
                      base_board=base_board, inputs={**inputs, "ref": whole[pm[POOL]].copy()}, ref_stage=ref_stage,
                      truths=truths, target=target, facts=facts)


def context(tr: Transition, s: int) -> harness.BuildContext:
    return harness.BuildContext(task=tr.view, board=tr.cboard, inputs=tr.inputs, seed=harness.BUILD_SEED_OFFSET + s)


# ------------------------------------------------------------------------------------------------ the arms
def size_and_half(ctx: harness.BuildContext, st: Stage, target: dict, what: str) -> tuple[Stage, dict]:
    """Row H29b's steps for D_f after the expression (h29b.build_pool_arms's own calls for its arm organ:f): the
    rescale about the centroid to D_f's size, row H25's half, the checked builder."""
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    Cs, rec = h29b.size_to(st.coords, target["target"], h29b.SIZING, target["ratio"], target["base"],
                           target["source"])
    rows = h25.kept_rows(len(Cs), seed, 1)
    return h29b._checked(ctx, np.asarray(st.X)[rows], Cs[rows], what), rec


def draw_record(ctx: harness.BuildContext, st: Stage) -> dict:
    """The draw as h29b.build_pool_arms records it: the practice count against the pool's rows."""
    n_pool = int(ctx.inputs["ref"].n_obs)
    return {"pool_cells": n_pool, "drawn": int(st.n), "short": bool(st.n > n_pool),
            "repeated_rows_at_least": int(max(0, st.n - n_pool))}


def build_control(ctx: harness.BuildContext, target: dict) -> tuple[Stage, dict, dict]:
    """The control alone: row H29b's D_f on D's context (`harness._copy`, then size_and_half)."""
    st = harness._copy(ctx)
    stage, rec = size_and_half(ctx, st, target, CONTROL)
    return stage, rec, draw_record(ctx, st)


def build_arms(ctx: harness.BuildContext, target: dict, strengths=STRENGTHS) -> tuple[dict, dict]:
    """Every arm of one seed on D's context: the control (`harness._copy`) and, per strength, row H28's step and its
    placebo (h28.trend_expression), each then through size_and_half; and the identity checks the read's condition v5
    holds the arms to. Strength 1 must be among the strengths: its arm is held against the harness's own b0_trend on
    the same context."""
    if LO not in strengths:
        raise ValueError("strength 1 must be among the strengths")
    pre = {CONTROL: (harness._copy(ctx), None, None, {"strength": None})}
    for s in strengths:
        pre[arm(s)] = h28.trend_expression(ctx, s)
        pre[placebo(s)] = h28.trend_expression(ctx, s, placebo_arm=True)
    ref_b0 = harness._b0_t2(ctx, trend=True)  # the harness's own b0_trend on the same context (a check, not scored)
    pool, rows, _, _ = harness.draw_copy(ctx)
    drawn = counts_of(pool[rows], list(ctx.board.genes))
    ctrl, one = pre[CONTROL][0], pre[arm(LO)][0]
    checks = {"s1_equals_harness_b0_trend": bool(np.array_equal(one.X, ref_b0.X)
                                                 and np.array_equal(one.coords, ref_b0.coords)),
              "same_draw": bool(all(np.array_equal(v[2], rows) for k, v in pre.items() if k != CONTROL)),
              "pre_size_coords_equal": bool(all(np.array_equal(v[0].coords, ctrl.coords) for v in pre.values())),
              # count entries the step changed among the drawn cells (0 at strength 0: the counts are the control's)
              "changed_counts": {k: int((v[1] != drawn).nnz) for k, v in pre.items() if k != CONTROL}}
    out, recs = {}, {}
    for name, (st, _, _, rec) in pre.items():
        out[name], srec = size_and_half(ctx, st, target, name)
        recs[name] = {"size": srec, "trend": rec}
    checks["post_size_coords_equal"] = bool(all(np.array_equal(v.coords, out[CONTROL].coords) for v in out.values()))
    checks["n_cells"] = {k: int(v.n) for k, v in out.items()}
    return out, {"records": recs, "checks": checks, "draw": draw_record(ctx, ctrl)}


def score_arms(tr: Transition, s: int, preds: dict) -> tuple[dict, dict, dict]:
    """Rows H28's and H29b's scoring of one seed: per truth, the floor and its anchors (shared by the arms), then
    every arm; the harness's bake protocol and scorer build main."""
    scores, skills, groups = {}, {}, {}
    for q, truth in tr.truths.items():
        board = replace(tr.base_board, truth_cells=truth.n)
        fl = score_board(board, subsample(tr.ref_stage, SUBSAMPLE_FRAC, s), truth, tr.ref_stage, seed=s,
                         anchors="computed", build=harness.BUILD, protocol=harness.PROTOCOL)
        scores[q] = {harness.FLOOR: fl.score}
        skills[q] = {harness.FLOOR: fl.skills}
        groups[q] = {harness.FLOOR: dict(fl.group_skill)}
        for name, pred in preds.items():
            r = score_board(board, pred, truth, tr.ref_stage, seed=s, anchors=fl.anchors, build=harness.BUILD,
                            protocol=harness.PROTOCOL)
            scores[q][name], skills[q][name] = r.score, r.skills
            groups[q][name] = dict(r.group_skill)
    return scores, skills, groups


# ------------------------------------------------------------------------------------------------ the gate's factors
class _StepRead(Exception):
    """Raised in place of momentum_shift once its arguments are read: the step stops before any count is touched."""


def step_arguments(ctx: harness.BuildContext, strength: float) -> tuple[float, dict[str, np.ndarray]]:
    """The damping and the log-factor mapping (lineage -> vector, `logf_by_type`) that row H28's step at `strength`
    gives vec.t1_scope.momentum_shift on this context: read through h28.trend_expression itself, which is stopped at
    that call. Both come from the input stages and the draw only. No count is touched and no arm is built;
    vec.t1_scope is as it was afterwards. The stand-in keeps momentum_shift's signature, whose defaults (eps,
    min_cells) b0_trend_factors reads before the call."""
    seen: list[tuple[float, dict[str, np.ndarray]]] = []
    own = t1_scope.momentum_shift

    @functools.wraps(own)
    def stop(start, ref, target_like, damping, *a, logf_by_type=None, **k):
        seen.append((float(damping), {t: np.array(v, dtype=np.float64) for t, v in (logf_by_type or {}).items()}))
        raise _StepRead

    t1_scope.momentum_shift = stop
    try:
        h28.trend_expression(ctx, strength)
    except _StepRead:
        pass
    finally:
        t1_scope.momentum_shift = own
    if len(seen) != 1:
        raise RuntimeError("the trend step did not reach momentum_shift exactly once")
    return seen[0]


def factor_identity(lo: tuple[float, dict], hi: tuple[float, dict], ratio: float = HI / LO) -> dict:
    """G2 in one seed's context, exactly and with no tolerance. lo and hi are step_arguments at the lower and the
    higher strength. For every lineage of the mapping, the log-factor vector given at the higher strength equals
    `ratio` times the vector given at the lower one, element by element; and the damping given at the higher strength
    equals `ratio` times the damping given at the lower one."""
    (d_lo, m_lo), (d_hi, m_hi) = lo, hi
    rec = {"ratio": float(ratio), "lineages": sorted(m_lo), "damping": [float(d_lo), float(d_hi)], "ok": False}
    if sorted(m_lo) != sorted(m_hi):
        return {**rec, "why": "the two strengths' mappings name different lineages"}
    off = [k for k in sorted(m_lo)
           if m_lo[k].shape != m_hi[k].shape or not np.array_equal(m_hi[k], ratio * m_lo[k])]
    if off:
        return {**rec, "why": f"the log-factor vectors of {off} are not {ratio:g} times the lower strength's"}
    if not d_hi == ratio * d_lo:
        return {**rec, "why": f"the damping {d_hi!r} is not {ratio:g} times the lower strength's {d_lo!r}"}
    return {**rec, "ok": True, "why": None}


# ------------------------------------------------------------------------------------------------ the committed D_f
def load_h29b_reading(path: Path | None = None) -> dict | None:
    """Row H29b's committed reading: the file of sha256 H29B_READING_SHA256 and no other (None otherwise)."""
    try:
        raw = Path(H29B_READING if path is None else path).read_bytes()
        return json.loads(raw) if hashlib.sha256(raw).hexdigest() == H29B_READING_SHA256 else None
    except (OSError, ValueError):
        return None


def committed_reasons(committed: dict | None, tag: str) -> list[str]:
    if (not isinstance(committed, dict) or committed.get("row") != h29b.ROW or committed.get("void") is not False
            or committed.get("outcome") != h29b.READ):
        return [f"{tag}: row H29b's committed reading is not at hand (missing, not the file its row records, or void)"]
    return []


def reproduces(rows: dict[str, list[dict]], committed: dict, tag: str) -> tuple[list[str], dict]:
    """G1, bit for bit: the control's plain mean per world and unit, by row H29b's own arithmetic (as h29b.cell forms
    a mean: h29b.series with the world's truths and the unit's entry of h29b.UNITS, then the mean over the seeds),
    equals row H29b's committed reading.json (tables[world][unit].mean under organ:f) exactly: eighteen entries, no
    tolerance."""
    why, n, worst = [], 0, 0.0
    for w, ts in G1_WORLDS:
        for u in G1_UNITS:
            want = committed["tables"][w][u]["mean"][D_F]
            got = float(h29b.series(rows, h29b.UNITS[u], ts, CONTROL).mean())
            n += 1
            if got != want:
                worst = max(worst, abs(got - want))
                why.append(f"{tag}: {w} {u}: the control's mean is {got - want:+.3e} from row H29b's committed D_f")
    return why, {"means": n, "unequal": len(why), "largest_abs_diff": worst}


def gate_reasons(rows: dict[str, list[dict]], committed: dict | None) -> tuple[list[str], dict]:
    """Why the gate fails (empty when it passes), and what it compared. G1: reproduces, on the control's rows of both
    transitions and every seed. G2: factor_identity in every seed's context of each transition."""
    dist: dict = {"G1": None, "G2": None}
    why = committed_reasons(committed, "G1")
    for key in TASK_KEYS:
        if key not in rows:
            why.append(f"G1: missing transition {key}")
        elif [r["seed"] for r in rows[key]] != list(SEEDS):
            why.append(f"G1: {key}: seeds {[r['seed'] for r in rows[key]]}, not {list(SEEDS)}")
    if why:
        return why, dist
    w1, dist["G1"] = reproduces(rows, committed, "G1")
    w2, lineages, n = [], set(), 0
    for key in TASK_KEYS:
        for r in rows[key]:
            f = r["factors"]
            n += 1
            lineages |= set(f.get("lineages") or [])
            if f.get("ok") is not True:
                w2.append(f"G2: {key} seed {r['seed']}: {f.get('why')}")
    dist["G2"] = {"contexts": n, "failed": len(w2), "lineages": sorted(lineages)}
    return w1 + w2, dist


# ------------------------------------------------------------------------------------------------ notes (gate nothing)
def _fields_differ(mine, theirs, fields) -> bool:
    return not (isinstance(mine, dict) and isinstance(theirs, dict)
                and all(k in mine and k in theirs and mine[k] == theirs[k] for k in fields))


def same_facts(key: str, f: dict, committed: dict, tag: str = "note") -> list[str]:
    """A transition's target-free records against row H29b's committed reading, exactly: the truths of the two
    worlds and D's pool (cell numbers and core shares, as row H27's V3 compares them), D's record on the reference
    stage, both rms ratios and D_f's size; and the size result's reference is the loaded cloud (row H27's V7)."""
    why = []
    for q in TRUTHS:
        if not h27._same(f["truth_info"].get(q), committed["truth_info"][key].get(q), "cells"):
            why.append(f"{tag}: {key}: truth {q} differs from row H29b's (cells or core share)")
    if not h27._same(f["pred_info"].get(POOL), committed["pool_info"][key].get(POOL), "pool_cells"):
        why.append(f"{tag}: {key}: D's pool differs from row H29b's (rows or core share)")
    if _fields_differ(f["organ_record"].get("reference"), committed["organ_records"][key].get("reference"),
                      ORGAN_FIELDS):
        why.append(f"{tag}: {key}: D's record on the reference stage differs from row H29b's")
    if _fields_differ(f["ratios"], committed["ratios"][key], h29b.POOLS):
        why.append(f"{tag}: {key}: rms(pool) / rms(core) differs from row H29b's")
    if _fields_differ(f["d_targets"].get(D_F), committed["d_targets"][key].get(D_F), TARGET_FIELDS):
        why.append(f"{tag}: {key}: D_f's size differs from row H29b's")
    if f["v7"].get("ok") is not True:
        why.append(f"{tag}: {key}: the reference cloud is not the cloud the size result measured")
    return why


def pool_draw(rows: list[dict]) -> dict:
    """The control's draw over the seeds, as row H29b's read reports a pool's (its pools_short)."""
    return {"pool_cells": rows[0]["draw"]["pool_cells"], "drawn": sorted({r["draw"]["drawn"] for r in rows}),
            "short_in_seeds": sum(bool(r["draw"]["short"]) for r in rows),
            "repeated_rows_at_least": max(r["draw"]["repeated_rows_at_least"] for r in rows)}


def same_draw(key: str, rows: list[dict], committed: dict, tag: str = "note") -> list[str]:
    if pool_draw(rows) != committed["pools_short"][key].get(POOL):
        return [f"{tag}: {key}: the control's draw from D's pool differs from row H29b's"]
    return []


def at_size(key: str, seed: int, name: str, rec: dict, target: dict, tag: str = "note") -> list[str]:
    """A build at D_f's size (h29b.size_checks's two conditions)."""
    why = []
    if not abs(rec["log_rms_after_vs_target"]) <= SIZE_LOG_TOL:
        why.append(f"{tag}: {key} seed {seed}: {name} is not at D_f's size")
    if not abs(rec["target_rms_full"] - target["target"]) <= SIZE_LOG_TOL * abs(target["target"]):
        why.append(f"{tag}: {key} seed {seed}: {name}'s size is not D_f's field")
    return why


def notes(facts: dict[str, dict], rows: dict[str, list[dict]], committed: dict | None, size_of) -> list[str]:
    """Gating nothing: each transition's target-free records against row H29b's committed ones (same_facts,
    same_draw) and the control at D_f's size in every seed (at_size; `size_of` gives a row's size record of the
    control). Every difference is a note; none is expected."""
    try:
        if committed_reasons(committed, "note"):
            return ["note: row H29b's committed reading is not at hand; its target-free records are not compared"]
        out = []
        for key in TASK_KEYS:
            out += same_facts(key, facts[key], committed) + same_draw(key, rows[key], committed)
            for r in rows[key]:
                out += at_size(key, r["seed"], "the control", size_of(r), facts[key]["d_targets"][D_F])
        return out
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:
        return [f"note: the target-free records cannot be compared ({type(e).__name__}: {e})"]


# ------------------------------------------------------------------------------------------------ the fixed words
SCREEN = "A screen, not evidence of development."
NO_PRIOR = "No row on record had read the trend step on D's pool at any strength before this one."
FOLLOWS = {  # row H34's eighth cell, "What follows, fixed now"
    NOT_STOPPED: (f"What follows: step (3), the builder option under ruling row {RULING_ROW}, is coded and reviewed, "
                  "then the identity builds (step 4), then the file (step 5), each as the line's Build passage has "
                  "it; NOT STOPPED uploads nothing, the file going in only under the line's Rule and on the "
                  "captain's word at the time."),
    STOPPED: "What follows: nothing is built.",
    VOID: "What follows: nothing is built.",
}
RULE = (f"RULE, fixed in the row before its run, on fold E10.5 (every difference paired by seed; a difference of "
        f"exactly 0 is neither above nor below 0). Clause (a) holds in a world if {arm(HI)} minus {arm(LO)} has a "
        f"paired mean below 0 with at least {MIN_SEEDS} of the 10 seeds below 0 in that world. The S1 bar (row "
        f"H28's) is met in a world if {arm(HI)} minus the control D_f has a paired mean above +{BAR} with at least "
        f"{MIN_SEEDS} of the 10 seeds above 0 in that world. Clause (b) holds unless the S1 bar is met in own and "
        "also in CLEAN; the bar failing in one world is enough for clause (b). STOPPED if clause (a) holds in own "
        "or in CLEAN, or clause (b) holds; NOT STOPPED otherwise. VOID is decided first, and under VOID the rule is "
        "not read. Arm minus placebo and the DE terms are labelled diagnostics and gate nothing.")
L0 = ("Row H34 (commit {commit}). The practice read the heart's fourth search line requires before any build: row "
      "H28's trend arms on row H29b's D_f form (D's pool at its practice forecast size, the half). One truth (E10.5; "
      "the two backtests merged seed by seed into fold E10.5), worlds own and CLEAN, seeds 0 to 9, released stages "
      "only. No figure here is to be read beside a test score, and no test score enters the row, its rule or any "
      "value. The read builds nothing.")
L1 = (f"Arms: the control is the copy D_f (row H29b's {D_F}: the copy draw from D's pool on the transition's "
      "reference stage, sized to the core forecast times rms(D pool) / rms(core) there, then the half); "
      f"{arm(LO)} and {arm(HI)} are row H28's expression step at strengths {LO:g} and {HI:g}, made before the size "
      "step, with D's pool as the step's reference and earlier stages whole; "
      f"{placebo(LO)} and {placebo(HI)} are the same factors with each shifted lineage's genes deranged within ten "
      f"level bins (seed {PERM_SEED}). Counted: eight tries (two strengths and two placebos on two transitions); the "
      "control is a reference.")
GATE_LINES = {
    PASS: (f"GATE PASS: G1 and G2 hold (the control D_f reproduces row H29b's committed D_f bit for bit; the "
           f"strength-{HI:g} factors equal {HI / LO:g} times the strength-{LO:g} factors); the run may follow at "
           "this commit."),
    FAIL: "GATE FAIL: G1 or G2 does not hold; the run step refuses, no try is counted and the outcome is VOID.",
}


def outcome_line(outcome: str, why=(), clauses=()) -> str:
    """The outcome line, wording fixed."""
    if outcome == VOID:
        return f"Try {ROW}: VOID: {'; '.join(why)}; nothing is built."
    if outcome == STOPPED:
        return f"Try {ROW}: STOPPED ({', '.join(clauses)}): {FILE} is not built."
    return f"Try {ROW}: NOT STOPPED."


def follows_line(outcome: str) -> str:
    return f"{SCREEN} {FOLLOWS[outcome]}"


# ------------------------------------------------------------------------------------------------ the gate step
def gate_text(g: dict) -> str:
    lines = [f"Row {ROW}'s gate (commit {g['commit']}): first, counting nothing. No counted arm and no placebo is "
             "built or scored in it: the one build it scores is the control D_f (the copy draw from D's pool at its "
             "practice forecast size, the half), a reference whose figures row H29b committed; the factors it checks "
             "come from the input stages and the draw only.",
             f"G1: the control, built and scored on both transitions and seeds 0 to 9 against the four truths, "
             f"against the plain ten-seed means under {D_F} in row H29b's reading.json (sha256 "
             f"{H29B_READING_SHA256}): eighteen entries, by row H29b's own arithmetic, exactly, with no tolerance.",
             f"G2: in every seed's context, every lineage's log-factor vector given to momentum_shift at strength "
             f"{HI:g} against {HI / LO:g} times the one given at strength {LO:g}, and the damping likewise; exactly."]
    d = g.get("distances") or {}
    if d.get("G1"):
        lines.append(f"G1 found: {d['G1']['means']} means compared, {d['G1']['unequal']} unequal (largest "
                     f"difference {d['G1']['largest_abs_diff']:.3e}).")
    if d.get("G2"):
        lines.append(f"G2 found: {d['G2']['contexts']} contexts, {d['G2']['failed']} failed; the mapping's lineages: "
                     f"{', '.join(d['G2']['lineages']) or 'none'}.")
    lines.append("Beside the gate, gating nothing (the target-free records against row H29b's committed ones): "
                 + ("no difference." if not g.get("notes") else "; ".join(g["notes"])))
    lines.append(GATE_LINES[g["verdict"]])
    if g["verdict"] == FAIL:
        lines += ["WHY: " + "; ".join(g["why"]), outcome_line(VOID, g["why"]), follows_line(VOID)]
    return "\n".join(lines) + "\n"


def gate_seed(tr: Transition, s: int, data_dir: Path) -> dict:
    """One seed of the gate: the control built and scored; the step's arguments at the two strengths read, no arm
    built."""
    t0 = time.time()
    ctx = context(tr, s)
    with harness.guarded_reads(harness.REPO_ROOT, data_dir):
        stage, rec, draw = build_control(ctx, tr.target)
        factors = factor_identity(step_arguments(ctx, LO), step_arguments(ctx, HI))
    scores, skills, groups = score_arms(tr, s, {CONTROL: stage})
    return {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": {CONTROL: int(stage.n)},
            "size_record": rec, "draw": draw, "factors": factors, "scores": scores, "skills": skills,
            "groups": groups, "seconds": round(time.time() - t0, 1)}


def run_gate(out: Path, data_dir: Path, commit: str | None, h29b_reading: Path | None = None) -> int:
    """The gate step: the refusals first, then both transitions and every seed in this one process."""
    refuse_state(commit, "gate")
    out = Path(out)
    if (out / GATE_JSON).exists() or (out / GATE_TXT).exists():
        sys.exit(f"{ROW} is refused: {out} already holds a gate record")
    committed = load_h29b_reading(h29b_reading)
    rows, facts, why = {}, {}, []
    for key in TASK_KEYS:
        try:
            tr = load_transition(key, data_dir)
        except t2_build.UnmappedT2LabelError as e:
            why.append(f"G1: {key}: the control cannot be built: labels without a T2 chamber category on the "
                       f"reference stage: {sorted(e.labels)}")
            continue
        facts[key], rows[key] = tr.facts, []
        for s in SEEDS:
            rows[key].append(gate_seed(tr, s, data_dir))
            print(f"{key} seed {s} ({rows[key][-1]['seconds']} s) done", flush=True)
        del tr
    state = git_state()
    if state.get("tree_clean") is not True or state.get("commit") != commit:
        sys.exit(f"{ROW} is refused: tracked files or HEAD changed while the gate ran; no gate record is written")
    dist: dict = {}
    if not why:
        try:
            why, dist = gate_reasons(rows, committed)
        except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:  # a malformed committed record
            why = [f"G1: row H29b's committed reading cannot be read ({type(e).__name__}: {e})"]
    g = {"row": ROW, "axis": GATE_AXIS, "commit": commit, "provenance": state,
         "label": "the gate of row H34: counting nothing; the control D_f alone is built and scored, and its figures "
                  "are row H29b's own", "verdict": FAIL if why else PASS, "why": why, "distances": dist,
         "notes": notes(facts, rows, committed, lambda r: r["size_record"]) if len(facts) == len(TASK_KEYS) else [],
         "h29b_reading_sha256": H29B_READING_SHA256, "strengths": list(STRENGTHS), "transitions": facts, "rows": rows,
         "source": SOURCE}
    out.mkdir(parents=True, exist_ok=True)
    (out / GATE_JSON).write_text(json.dumps(g, indent=1, default=float))
    (out / GATE_TXT).write_text(gate_text(g))
    print(f"GATE {g['verdict']}", flush=True)
    if why:
        print(outcome_line(VOID, why), flush=True)
        print(follows_line(VOID), flush=True)
    return GATE_EXIT[g["verdict"]]


# ------------------------------------------------------------------------------------------------ the run
def run_transition(key: str, seeds, out: Path, data_dir: Path, commit: str | None, gate_dir: Path) -> dict:
    out = Path(out)
    path = out / PART
    refuse_run(path, commit, gate_dir)
    if key not in TASK_KEYS:
        sys.exit(f"{key}: not a transition of row {ROW}")
    tr = load_transition(key, data_dir)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        ctx = context(tr, s)
        with harness.guarded_reads(harness.REPO_ROOT, data_dir):
            preds, extra = build_arms(ctx, tr.target)
        scores, skills, groups = score_arms(tr, s, preds)
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": extra["checks"]["n_cells"],
               "checks": extra["checks"], "records": extra["records"], "draw": extra["draw"], "scores": scores,
               "skills": skills, "groups": groups, "seconds": round(time.time() - t0, 1)}
        rows_out.append(row)
        print(f"{key} seed {s} ({row['seconds']} s) done", flush=True)
    res = {"label": "row H34: row H28's trend step at strengths 1 and 1.5, and a placebo of each, on row H29b's D_f "
                    "form; a screen on one truth, not an adoption result", "axis": AXIS, "task": tr.key,
           "arms": list(BUILT), "truth_values": list(TRUTHS), **tr.facts, "strengths": list(STRENGTHS),
           "truth_seed": payoff.TRUTH_SEED, "half_seed_base": h25.HALF_SEED_BASE, "perm_seed": PERM_SEED,
           "departures": h27.DEPARTURES, "seeds": list(seeds), "rows": rows_out, "source": SOURCE,
           "provenance": git_state()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


def load_parts(root: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for p in sorted(Path(root).rglob(PART)):
        r = json.loads(p.read_text())
        if r.get("axis") == AXIS:
            out.setdefault(r["task"], []).append(r)
    return out


# ------------------------------------------------------------------------------------------------ the void
def identities(key: str, rows: list[dict]) -> list[str]:
    """v5, the arm is the arm the third cell defines: trend_s1.0's expression and coordinates before the size step
    are `harness._b0_t2(ctx, trend=True)`'s on the same context, bit for bit; every arm's draw, coordinates (before
    and after the size step and the half) and cell number are the control's; and no arm's d2_shape, occupancy_dice or
    scale_log_ratio is more than 1e-9 from the control's in any seed and truth."""
    why = []
    for r in rows:
        c = r["checks"]
        tag = f"v5: {key} seed {r['seed']}"
        if c.get("s1_equals_harness_b0_trend") is not True:
            why.append(f"{tag}: {arm(LO)} before the size step is not the harness's b0_trend on the same context")
        for k in ("same_draw", "pre_size_coords_equal", "post_size_coords_equal"):
            if c.get(k) is not True:
                why.append(f"{tag}: {k.replace('_', ' ')} fails")
        if set(r["n_cells"]) != set(BUILT) or len(set(r["n_cells"].values())) != 1:
            why.append(f"{tag}: the arms' cell numbers are not the control's")
        for q in TRUTHS:
            for name in COUNTED:
                for m in SHAPE_TERMS:
                    if not abs(r["skills"][q][name][m] - r["skills"][q][CONTROL][m]) <= TOL:
                        why.append(f"{tag} {q}: {name}'s {m} is not the control's")
    return why


def trend_records(key: str, rows: list[dict], pool_cells: int) -> list[str]:
    """v5, the arms' records: each holds its transition's budget as it is before the strength is applied (relative
    1e-12) and D's pool as the step's reference; the arms shift the same lineages; no placebo leaves a gene on its
    own factor."""
    why = []
    want = PRACTICE_BUDGET[key]
    for r in rows:
        shifted = set()
        tag = f"v5: {key} seed {r['seed']}"
        for s in STRENGTHS:
            for name, is_placebo in ((arm(s), False), (placebo(s), True)):
                rec = r["records"][name]["trend"]
                if not abs(rec["budget_days"] - want) <= FACTOR_RTOL * want:
                    why.append(f"{tag}: {name}'s budget is {rec['budget_days']}, not {want}")
                if rec["ref_cells"] != pool_cells:
                    why.append(f"{tag}: {name}'s reference is not D's pool")
                shifted.add(tuple(rec["shifted"]))
                if is_placebo and any(v for v in rec["fixed_points"].values()):
                    why.append(f"{tag}: {name} leaves a gene on its own factor")
        if len(shifted) != 1:
            why.append(f"{tag}: the arms shift different lineages")
    return why


def void_reasons(parts: dict[str, list[dict]], committed: dict | None, gate: dict | None, commit: str | None
                 ) -> tuple[list[str], dict]:
    """Why the read is void (empty when it may be read), by the eighth cell's conditions that code can hold, and
    what G1 found on the run's own control.
      v3  the gate failed, or its record is missing or of another commit than the run commit; and G1 repeated on the
          run's own control, the control the rule's differences are paired with
      v4  the read is not the read registered: a transition, seed, truth or arm is missing or a seed occurs twice
          (or a part of another transition is among them); the parts come from more than one commit, from a commit
          other than the run commit or from an unclean tree
      v5  an arm is not the arm the third cell defines (identities, trend_records)
    (v1, v2 and v6 are not the read's to decide: the line's own void, the row before its steps, which every step's
    refusal holds, and no test score in the row.)"""
    dist: dict = {"G1": None}
    why = gate_record_reasons(gate, commit)
    every = [p for ps in parts.values() for p in ps]
    complete = True
    for key in TASK_KEYS:  # v4: completeness
        ps = parts.get(key)
        if not ps:
            why.append(f"v4: missing transition {key}")
            complete = False
            continue
        rows, twice = h26.merge_rows(ps)
        if twice:
            why.append(f"v4: {key}: seeds {twice} occur in more than one part")
            complete = False
        if [r["seed"] for r in rows] != list(SEEDS):
            why.append(f"v4: {key}: seeds {[r['seed'] for r in rows]}, not {list(SEEDS)}")
            complete = False
        for r in rows:
            for field in ("scores", "skills", "groups"):
                gone = sorted({f"{q}/{a}" for q in TRUTHS for a in ALL if a not in r.get(field, {}).get(q, {})})
                if gone:
                    why.append(f"v4: {key} seed {r['seed']}: missing from {field}: {gone[:6]}")
                    complete = False
            if set(r.get("records", {})) != set(BUILT) or "checks" not in r or "n_cells" not in r:
                why.append(f"v4: {key} seed {r['seed']}: records, checks or cell numbers missing")
                complete = False
    if set(parts) - set(TASK_KEYS):
        why.append(f"v4: parts of other transitions {sorted(set(parts) - set(TASK_KEYS))}")
    if not every:
        return why, dist
    if len({p["provenance"]["commit"] for p in every}) > 1:  # v4: the commit and the tree
        why.append("v4: the parts come from more than one commit")
    if not all(p["provenance"]["tree_clean"] is True for p in every):
        why.append("v4: a part comes from an unclean tree")
    if not commit or any(str(p["provenance"]["commit"]) != commit for p in every):
        why.append(f"v4: a part comes from a commit other than the run commit {commit!r}")
    if not complete:
        return why, dist
    merged = {}
    for key in TASK_KEYS:  # v5
        rows = merged[key] = h26.merge_rows(parts[key])[0]
        why += identities(key, rows) + trend_records(key, rows, parts[key][0]["pred_info"][POOL]["pool_cells"])
    g1 = committed_reasons(committed, "v3")  # v3: G1 on the run's own control
    if not g1:
        g1, dist["G1"] = reproduces(merged, committed, "v3: G1 on the run's own control")
    return why + g1, dist


# ------------------------------------------------------------------------------------------------ the tables
def compare(rows: dict[str, list[dict]], keys, truths, a: str, b: str) -> dict:
    """h28.compare (a minus b, paired by seed: the board score and every term), with the seeds below zero beside the
    seeds above it."""
    out = h28.compare(rows, keys, truths, a, b)
    d = h28.series(rows, keys, truths, a) - h28.series(rows, keys, truths, b)
    out["score"]["behind"] = int((d < 0).sum())
    return out


def tables(rows: dict[str, list[dict]]) -> dict:
    return {w: {u: {pname(a, b): compare(rows, ks, ts, a, b) for a, b in PAIRS} for u, ks in UNITS.items()}
            for w, ts in WORLDS}


# ------------------------------------------------------------------------------------------------ the rule
def below(s: dict) -> bool:
    """Clause (a) in one world: a paired mean below 0 with at least MIN_SEEDS of 10 seeds below 0."""
    return bool(s["mean"] < 0 and s["behind"] >= MIN_SEEDS)


def stop_rule(tab: dict) -> dict:
    """The rule on fold E10.5. Clause (a) holds in a world if trend_s1.5 minus trend_s1.0 is below there (`below`).
    Row H28's S1 bar is met in a world if trend_s1.5 minus the control passes there (h28.passes). Clause (b) holds
    unless the bar is met in own and also in CLEAN. STOPPED if clause (a) holds in own or in CLEAN, or clause (b)
    holds."""
    def fold(w, a, b):
        return tab[w][FOLD][pname(a, b)]["score"]
    a = {w: below(fold(w, arm(HI), arm(LO))) for w, _ in WORLDS}
    bar = {w: bool(h28.passes(fold(w, arm(HI), CONTROL))) for w, _ in WORLDS}
    holds = {"a": any(a.values()), "b": not all(bar.values())}
    clauses = [c for c in ("a", "b") if holds[c]]
    return {"a_in_world": a, "bar_met": bar, "holds": holds, "clauses": clauses, "stop": bool(clauses)}


def read_notes(parts: dict[str, list[dict]], committed: dict | None) -> list[str]:
    """Gating nothing: `notes` on the run's parts, and a transition whose arms shift no lineage."""
    try:
        rows = {k: h26.merge_rows(parts[k])[0] for k in TASK_KEYS}
        out = notes({k: parts[k][0] for k in TASK_KEYS}, rows, committed, lambda r: r["records"][CONTROL]["size"])
        return out + [f"note: {k} seed {r['seed']}: the arms shift no lineage" for k in TASK_KEYS for r in rows[k]
                      if not r["records"][arm(LO)]["trend"]["shifted"]]
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:
        return [f"note: the run's records cannot be compared ({type(e).__name__}: {e})"]


def reading(parts: dict[str, list[dict]], committed: dict | None, gate: dict | None, commit: str | None) -> dict:
    try:
        why, dist = void_reasons(parts, committed, gate, commit)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:  # a malformed part or record
        why, dist = [f"v4: a part or a committed record cannot be read ({type(e).__name__}: {e})"], {}
    g = gate if isinstance(gate, dict) else {}
    res = {"row": ROW, "commit": commit, "void": bool(why), "why": why, "guards": dist, "prior": NO_PRIOR,
           "screen": SCREEN, "rule": RULE, "h29b_reading_sha256": H29B_READING_SHA256,
           "gate": {k: g.get(k) for k in ("verdict", "commit", "why", "distances", "notes", "provenance")}}
    if why:
        res.update(outcome=VOID, clauses=[])
    else:
        rows = {key: h26.merge_rows(parts[key])[0] for key in TASK_KEYS}
        tab = tables(rows)
        st = stop_rule(tab)
        res.update(tables=tab, statements=st, clauses=st["clauses"], outcome=STOPPED if st["stop"] else NOT_STOPPED,
                   notes=read_notes(parts, committed),
                   budgets={f"{s:g}": {"test": h28.budget_at_test(s), **{k: s * b for k, b in PRACTICE_BUDGET.items()}}
                            for s in STRENGTHS},
                   trend_record={k: rows[k][0]["records"][arm(LO)]["trend"] for k in TASK_KEYS},
                   changed_counts={k: rows[k][0]["checks"].get("changed_counts") for k in TASK_KEYS},
                   d_targets={k: parts[k][0]["d_targets"] for k in TASK_KEYS},
                   pred_info={k: parts[k][0]["pred_info"] for k in TASK_KEYS},
                   truth_info={k: parts[k][0]["truth_info"] for k in TASK_KEYS})
    res["line"] = outcome_line(res["outcome"], why, res["clauses"])
    res["follows"] = follows_line(res["outcome"])
    return res


def _st(s: dict) -> str:
    return (f"{s['mean']:+.2f} [{s['lo']:+.2f}, {s['hi']:+.2f}] {s['ahead']}/10 above, {s['behind']}/10 below "
            f"(se {s['se']:.3f})")


def _terms(c: dict) -> str:
    return ", ".join(f"{m} {c[m]:+.4f}" for m in (*TERMS, *SHAPE_TERMS))


def rule_lines(res: dict) -> list[str]:
    tab, st = res["tables"], res["statements"]

    def fold(w, a, b):
        return tab[w][FOLD][pname(a, b)]
    lines = [RULE]
    for w, _ in WORLDS:
        lines.append(f"Clause (a) in {w}: {'HOLDS' if st['a_in_world'][w] else 'does not hold'} ({arm(HI)} minus "
                     f"{arm(LO)}, fold E10.5: {_st(fold(w, arm(HI), arm(LO))['score'])})")
    for w, _ in WORLDS:
        lines.append(f"The S1 bar in {w}: {'MET' if st['bar_met'][w] else 'not met'} ({arm(HI)} minus the control "
                     f"D_f, fold E10.5: {_st(fold(w, arm(HI), CONTROL)['score'])})")
    lines.append(f"Clause (a) holds in own or in CLEAN: {str(st['holds']['a']).upper()}. Clause (b) holds (the bar "
                 f"is not met in both worlds): {str(st['holds']['b']).upper()}.")
    for s in STRENGTHS:
        lines.append(f"DIAGNOSTIC, gates nothing: {arm(s)} minus {placebo(s)}, fold E10.5: "
                     + "; ".join(f"{w} {_st(fold(w, arm(s), placebo(s))['score'])}" for w, _ in WORLDS))
    for a, b in ((arm(LO), CONTROL), (arm(HI), CONTROL), (arm(HI), arm(LO))):
        lines.append(f"DIAGNOSTIC, gates nothing: the DE terms of {a} minus {b}, fold E10.5: "
                     + "; ".join(f"{w} " + ", ".join(f"{m} {fold(w, a, b)[m]:+.4f}" for m in DE_TERMS)
                                 for w, _ in WORLDS))
    return lines


def text(res: dict) -> str:
    lines = [res["line"], res["follows"], L0.format(commit=res["commit"]), NO_PRIOR, L1]
    g = res.get("gate") or {}
    lines.append("GATE (first, counting nothing): "
                 + (f"{g['verdict']} at commit {g.get('commit')}." if g.get("verdict") else "no record."))
    if res["void"]:
        lines += [RULE, "VOID, the conditions that hold, one by one:", *res["why"]]
        return "\n".join(lines) + "\n"
    g1 = res["guards"]["G1"]
    lines.append(f"G1 repeated on the run's own control, against row H29b's committed D_f: {g1['means']} plain "
                 f"means compared, {g1['unequal']} unequal.")
    lines.append("Beside the read, gating nothing (the target-free records against row H29b's committed ones): "
                 + ("no difference." if not res["notes"] else "; ".join(res["notes"])))
    lines.append("Budgets (B0's damping times the last gap, times the strength): "
                 + "; ".join(f"strength {s}: " + ", ".join(f"{k} {v:.4g} day" for k, v in b.items())
                             for s, b in res["budgets"].items()) + ".")
    lines += rule_lines(res)
    lines.append("TABLES, practice diagnostics that gate nothing (paired by seed: mean [min, max] seeds above and "
                 "below 0 of 10 (se); per-term mean skill differences):")
    for w, _ in WORLDS:
        for u in UNITS:
            for p, c in res["tables"][w][u].items():
                lines.append(f"{TAG} {w} {u} {p}: {_st(c['score'])}; {_terms(c)}")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------ command line
def _write_new(path: Path, content: str) -> None:
    """Write a file that must not exist: a reading is never overwritten."""
    with open(path, "x", encoding="utf-8") as f:
        f.write(content)


def cmd_gate(a) -> None:
    sys.exit(run_gate(Path(a.out), Path(a.data), a.commit, Path(a.h29b)))


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
    res = reading(load_parts(Path(a.runs)), load_h29b_reading(Path(a.h29b)), load_gate(Path(a.gate)), a.commit)
    out.mkdir(parents=True, exist_ok=True)
    _write_new(out / "reading.json", json.dumps(res, indent=1, default=float))
    _write_new(out / "reading.txt", text(res))
    print(res["line"])
    print(res["follows"])


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--commit", required=True, help="row H34's run commit (full hash; HEAD must be it)")
    g.add_argument("--data", default=str(harness.DATA_DIR))
    g.add_argument("--h29b", default=str(H29B_READING), help="row H29b's committed reading.json")
    r = sub.add_parser("run")
    r.add_argument("--task", required=True, help=f"one of {TASK_KEYS}")
    r.add_argument("--seeds", required=True, help="a range such as 0-1")
    r.add_argument("--out", required=True)
    r.add_argument("--gate", required=True, help="the gate step's folder (its gate.json must say PASS at --commit)")
    r.add_argument("--commit", required=True, help="row H34's run commit (full hash; HEAD must be it)")
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--gate", required=True)
    d.add_argument("--h29b", default=str(H29B_READING), help="row H29b's committed reading.json")
    d.add_argument("--out", default=str(RESULTS))
    d.add_argument("--commit", required=True, help="row H34's run commit (full hash; HEAD must be it)")
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
