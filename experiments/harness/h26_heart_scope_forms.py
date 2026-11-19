#!/usr/bin/env python
"""Try H26 of experiments/harness/TRIES.md: the heart's three candidate last-upload forms, which differ only in the
scope of the reference stage they draw from, scored on the s2 hedge study's three heart transitions against its five
truths each, with the scope differences the study ASSUMED to carry over now measured (docs/PLAN.md decision log
2026-11-19 (m), (n)).

The three forms (docs/S2_RULES.md "Amendment of 2026-11-18" and "T2 heart"):
  A  full_iso_half  the file in force: the copy's draw from the whole reference stage, expression unchanged, positions
                    scaled evenly on all three axes about the centroid, a uniform random half of the rows
  B  half_iso_half  A's form on the reference stage under vec.t2_build.radial_trim_mask(C, core, 0.5) (the half trim)
  C  core_iso_half  A's form under radial_trim_mask(C, core, 0.0) (the full trim)
Each is built as experiments/builds/build_test.py `build_heart` builds `--expression copy`: `vec.harness._copy` on a
context whose reference input is the pool (the harness's count rule for the practice board), then the rescale, then
the half as row H25's iso_half takes it (`h25_heart_forms.kept_rows`, seed HALF_SEED_BASE + scoring seed), then the
checked builder (`t2_build.build_submission`).
Reference, not counted:
  s1form_half       the first upload's form on the whole stage: `vec.harness._b0_trend` (row H25's s1form expression)
                    on the two harness backtests and the plain copy on E8.25_late -> E8.5 (no earlier stage there, as
                    in the study), an IN-PLANE rescale with z held, then the same half. It stands for the banked
                    first upload in the better-of-two.

The size is an ORACLE, shared by every arm and labelled so everywhere. One factor per transition, from the cores
(the five truths of a transition share one core: core cells are never trimmed): the even arms are scaled by
rms3d(truth core) / rms3d(reference core), and s1form_half in plane by the study's own factor
(xy RMS of the truth core over the reference core, z held; payoff.py's oracle_core_xy_factor). An arm's whole-cloud
size then follows from its own scope. A second size variant (each even arm at the oracle core size times its scope
hypothesis's full-over-core ratio) is LEFT OUT: the size procedure's code has no entry point that takes an outside
core size, so it would need a formula of this script's own (row H26 says so).

Truths and scoring are the study's (experiments/s2_hedge/payoff.py, imported as a module, `run_heart`): the five
truths per transition built by payoff's own `trim_mask` with its TRUTH_SEED; each truth's floor and anchors per seed,
shared by the arms; the harness's bake protocol and scorer build main; build seed BUILD_SEED_OFFSET + seed; 10 seeds.

The reading (`reading`), fixed in row H26 before the run. Void (nothing follows, no rerun under the row) if: a
transition, truth, seed or arm is missing or a seed repeats; the parts come from more than one commit, from a commit
other than the row's (read --commit) or from a dirty tree; a truth's cell number or core share differs from the
study's committed truth_info; a pool's cell number or core share differs from the study's committed pred_info (A:
hold, B: half_radial, C: core_radial); the in-plane oracle factor differs from the study's committed
oracle_core_xy_factor (relative 1e-12); or, on the two harness backtests with the truth as released, full_iso_half
differs from row H25's committed iso_half on the de_score, de_direction or mmd_u skill in any seed by more than 1e-9.
Then:
  R1  carry-over (descriptive, decides nothing): per transition and truth, half_iso_half minus full_iso_half and
      core_iso_half minus full_iso_half (mean, min, max, seeds ahead), beside the study's half_radial minus hold and
      core_radial minus hold; per trim, how many of the 15 cells agree in sign or both lie within 1.0 point of zero.
  R2  the better of two in practice space, no test score: V(X) = mean over seeds of max(s1form_half, X), per
      transition and truth; its expectation under the frozen section's five priors, per transition; the plain mean
      score; and s1form_half alone.
  R3  per metric group and metric: what the trims gain or lose against full_iso_half, per transition and truth.
  The pick rule, used only if the captain leaves the direction to the AI assistant: D3 = mean over the 15
  prior-by-transition cells of E V(half_iso_half) - E V(full_iso_half); DJ4 = the same with the two E10.5
  transitions merged into one fold first, then the two folds averaged. B if the run is valid and D3 >= 1.0 and
  DJ4 >= 1.0 board points; otherwise A, the file in force. C is reported and never picked.

  run   (DGX, repo root) one transition, a seed range; one JSON under --out
  read  apply the reading to every part under --runs; write reading.json and reading.txt under --out
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from vec import harness, prep, t2_build
from vec.adoption import SEEDS
from vec.boards import load_boards
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parents[2]  # this checkout (the experiments beside the script)
sys.path.insert(0, str(HERE / "experiments" / "harness"))
import h25_heart_forms as h25  # noqa: E402  (row H25's halving, imported, not copied)


def _load_payoff():
    spec = importlib.util.spec_from_file_location("s2_hedge_payoff", HERE / "experiments" / "s2_hedge" / "payoff.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


payoff = _load_payoff()

TASK_KEYS = ("T2h/<=E9.5->E10.5", "T2h/E8.25_late->E8.5", "T2h/<=E8.75->E10.5")  # payoff.HEART_TASK and its extras
HARNESS_KEYS = ("T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5")  # the two harness backtests (row H25's)
FOLDS = (("T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"), ("T2h/E8.25_late->E8.5",))  # one truth per fold (rule 9, J4)
TRUTHS = tuple(payoff.HEART_TRUTHS)  # own, radial_0.5, radial_0, thin_0.5, thin_0
OWN = "own"
REFERENCE = "s1form_half"  # not counted
ARMS = ("full_iso_half", "half_iso_half", "core_iso_half")  # counted: one try per arm and transition
BUILT = (REFERENCE, *ARMS)
ALL = (harness.FLOOR, *BUILT)
TRIM = {REFERENCE: None, "full_iso_half": None, "half_iso_half": 0.5, "core_iso_half": 0.0}  # radial_trim_mask frac
STUDY_TWIN = {"full_iso_half": "hold", "half_iso_half": "half_radial", "core_iso_half": "core_radial"}
STUDY_CONTROL = "hold"
CONTROL = "full_iso_half"  # A, the file in force
CANDIDATE = "half_iso_half"  # B, the only form the pick rule can return
TRIMMED = ("half_iso_half", "core_iso_half")
STUDY_RUNS = {"T2h/<=E9.5->E10.5": "experiments/s2_hedge/results/per_seed/heart-scope",
              "T2h/E8.25_late->E8.5": "experiments/s2_hedge/results/per_seed/heart-scope-e825",
              "T2h/<=E8.75->E10.5": "experiments/s2_hedge/results/per_seed/heart-scope-e875"}
H25_RUN = "experiments/harness/results/h25_heart_forms/tasks"
H25_TWIN = "iso_half"
EXPRESSION_ONLY = h25.EXPRESSION_ONLY  # de_score, de_direction, mmd_u
METRICS = h25.METRICS
GROUPS = h25.GROUPS
TOL = 1e-9
FACTOR_RTOL = 1e-12
NEAR_ZERO = 1.0  # R1: both within this many board points of zero count as agreeing
PICK_BAR = 1.0  # board points, for D3 and DJ4 alike
PRIORS = {  # docs/S2_RULES.md "T2 heart", "Why the amendment": the five priors over the truths
    "uniform": {q: 1 / len(TRUTHS) for q in TRUTHS},
    "own_third": {q: (1 / 3 if q == OWN else (2 / 3) / (len(TRUTHS) - 1)) for q in TRUTHS},
    "released": {OWN: 1.0},
    "radial_family": {q: 1 / 3 for q in ("own", "radial_0.5", "radial_0")},
    "thin_family": {q: 1 / 2 for q in ("thin_0.5", "thin_0")},
}
ORACLE = "oracle size"
AXIS = "h26-heart-scope-forms"
SOURCE = "experiments/harness/h26_heart_scope_forms.py (try H26)"


def safe(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9.]+", "_", key).strip("_")


def s1_expression(key: str) -> str:
    """The reference arm's expression step on a transition: b0_trend where the harness backtest has earlier stages,
    the plain copy on E8.25_late -> E8.5 (no earlier stage; the study built copies there)."""
    return "b0_trend" if key in HARNESS_KEYS else "copy"


# ------------------------------------------------------------------------------------------------ truths and pools
def truth_masks(TC: np.ndarray, tcore: np.ndarray) -> dict[str, np.ndarray]:
    """The study's five truths as row masks of the target: payoff.trim_mask at payoff.TRUTH_SEED, in the order and
    with the generators `payoff.run_heart` uses (truth i draws from default_rng([TRUTH_SEED, i]))."""
    out = {}
    for i, name in enumerate(payoff.HEART_TRUTHS):
        if name == OWN:
            out[name] = np.ones(len(TC), dtype=bool)
        else:
            mode, frac = name.split("_")
            out[name] = payoff.trim_mask(TC, tcore, float(frac), mode, np.random.default_rng([payoff.TRUTH_SEED, i]))
    return out


def pool_masks(PC: np.ndarray, pcore: np.ndarray) -> dict[str, np.ndarray]:
    """Each arm's pool as a row mask of the reference stage: whole, or vec.t2_build.radial_trim_mask (the builder's
    trim of the frozen rule's original second upload) to half the in-plane scope or to none."""
    return {arm: (np.ones(len(PC), dtype=bool) if TRIM[arm] is None
                  else t2_build.radial_trim_mask(PC, pcore, TRIM[arm])) for arm in BUILT}


def truth_info(C: np.ndarray, core: np.ndarray, keep: np.ndarray) -> dict:
    """As payoff.run_heart records a truth (cells, core share, in-plane scope)."""
    return {"cells": int(keep.sum()), "core_share": float(core[keep].mean()), "xy_scope": payoff.xy_scope(C[keep],
                                                                                                         core[keep])}


def pool_info(C: np.ndarray, core: np.ndarray, keep: np.ndarray) -> dict:
    """As payoff.run_heart records a prediction's pool (pool cells, core share, in-plane scope)."""
    return {"pool_cells": int(keep.sum()), "core_share": float(core[keep].mean()),
            "xy_scope": payoff.xy_scope(C[keep], core[keep])}


# ------------------------------------------------------------------------------------------------ the oracle size
def oracle_factors(TC: np.ndarray, tcore: np.ndarray, PC: np.ndarray, pcore: np.ndarray) -> dict:
    """The oracle size factors of one transition, both from the cores: `even` the 3-D RMS radius of the truth's core
    over the reference stage's core (the even arms); `in_plane` the study's own factor, the in-plane RMS of the
    truth's core over the reference's (s1form_half, z held; payoff.py's oracle_core_xy_factor)."""
    even = t2_build.rms_radius(TC[tcore]) / t2_build.rms_radius(PC[pcore])
    xy = payoff.xy_rms(TC[tcore]) / payoff.xy_rms(PC[pcore])
    return {"label": ORACLE, "even": float(even), "in_plane": float(xy)}


def even_rescale(C: np.ndarray, f: float) -> np.ndarray:
    """An even rescale about the centroid by f on all three axes, through the builder's own `rescale_to_rms`."""
    C = np.asarray(C, dtype=np.float64)
    return t2_build.rescale_to_rms(C, float(f) * t2_build.rms_radius(C))


def in_plane_rescale(C: np.ndarray, f: float) -> np.ndarray:
    """x and y scaled about their centroid by f, z untouched: the study's rescale (payoff.run_heart)."""
    C = np.asarray(C, dtype=np.float64).copy()
    c = C[:, :2].mean(0)
    C[:, :2] = c + (C[:, :2] - c) * float(f)
    return C


def form(X: np.ndarray, coords: np.ndarray, arm: str, factors: dict, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """An arm's rows from its harness build: the oracle rescale (even, or in plane for s1form_half), then the half of
    row H25 (`h25.kept_rows(n, seed, 1)`)."""
    C = (in_plane_rescale(coords, factors["in_plane"]) if arm == REFERENCE
         else even_rescale(coords, factors["even"]))
    rows = h25.kept_rows(len(C), seed, 1)
    return np.asarray(X)[rows], C[rows]


def build_arm(arm: str, ctx: harness.BuildContext, factors: dict, expression: str) -> Stage:
    """Arm `arm` on the context whose reference input is its pool: the harness build (`_copy`, or `_b0_trend` for
    s1form_half on a harness backtest), the form, then the checked builder."""
    base = harness._b0_trend if (arm == REFERENCE and expression == "b0_trend") else harness._copy
    st = base(ctx)
    seed = ctx.seed - harness.BUILD_SEED_OFFSET
    X, C = form(st.X, st.coords, arm, factors, seed)
    sub, rep = t2_build.build_submission(ctx.board, X, C, genes=list(ctx.board.genes))
    if not rep.ok:
        raise ValueError(f"{ctx.task.key} {arm}: the built file fails its check: {rep.errors}")
    return Stage.from_anndata(sub, ctx.board.genes, need_coords=True)


# ------------------------------------------------------------------------------------------------ the run
def run_transition(key: str, seeds, out: Path, data_dir: Path) -> dict:
    """One transition, as payoff.run_heart loads and scores it, with this row's four arms."""
    out = Path(out)
    path = out / "h26.json"
    if path.exists():
        sys.exit(f"{path} exists")
    t = payoff.heart_task(key)
    boards = load_boards()
    genes = boards[t.panel_board].genes
    ta = read_stage(data_dir / t.target, genes)
    prep.add_lineage(ta, source=harness.stem(t.target))
    tcore = np.asarray(ta.obs["core"], dtype=bool)
    TC = np.asarray(ta.obsm[COORD_KEY], dtype=float)[:, :3]
    ref_stage = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=True)
    tm = truth_masks(TC, tcore)
    truths = {q: Stage.from_anndata(ta[m], genes, need_coords=True) for q, m in tm.items()}
    t_info = {q: truth_info(TC, tcore, m) for q, m in tm.items()}
    inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    pool = inputs["ref"]
    pcore = np.asarray(pool.obs["core"], dtype=bool)
    PC = np.asarray(pool.obsm[COORD_KEY], dtype=float)[:, :3]
    pm = pool_masks(PC, pcore)
    pools = {arm: (pool if TRIM[arm] is None else pool[m].copy()) for arm, m in pm.items()}
    p_info = {arm: pool_info(PC, pcore, m) for arm, m in pm.items()}
    factors = oracle_factors(TC, tcore, PC, pcore)
    expression = s1_expression(key)
    base_board = harness.task_board(t, boards, truths[OWN].n, ref_stage.n)
    view, cboard = harness.TaskView.of(t), harness._candidate_board(base_board)
    rows_out = []
    for s in seeds:
        t0 = time.time()
        preds, n_cells = {}, {}
        for arm in BUILT:
            ctx = harness.BuildContext(task=view, board=cboard, inputs={**inputs, "ref": pools[arm]},
                                       seed=harness.BUILD_SEED_OFFSET + s)
            with harness.guarded_reads(harness.REPO_ROOT, data_dir):
                preds[arm] = build_arm(arm, ctx, factors, expression)
            n_cells[arm] = int(preds[arm].n)
        row = {"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "n_cells": n_cells, "scores": {},
               "skills": {}, "groups": {}}
        for q, truth in truths.items():
            board = replace(base_board, truth_cells=truth.n)  # as payoff.run_heart: the trimmed truth's own count
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
    res = {"label": f"try H26: {ORACLE} (an oracle shared by every arm), not an adoption result", "axis": AXIS,
           "task": t.key, "s1form_expression": expression, "pred_values": list(BUILT), "truth_values": list(TRUTHS),
           "truth_info": t_info, "pred_info": p_info, "oracle_factors": factors, "truth_seed": payoff.TRUTH_SEED,
           "half_seed_base": h25.HALF_SEED_BASE, "seeds": list(seeds), "rows": rows_out, "source": SOURCE,
           "provenance": payoff.provenance()}
    out.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)
    return res


# ------------------------------------------------------------------------------------------------ loading
def load_parts(root: Path) -> dict[str, list[dict]]:
    """Every part (one process's JSON) under root, by transition."""
    out: dict[str, list[dict]] = {}
    for p in sorted(Path(root).rglob("h26.json")):
        r = json.loads(p.read_text())
        if r.get("axis") == AXIS:
            out.setdefault(r["task"], []).append(r)
    return out


def merge_rows(parts: list[dict]) -> tuple[list[dict], list[int]]:
    """The parts' rows by seed, and the seeds that occur more than once."""
    by_seed: dict[int, dict] = {}
    twice = []
    for p in parts:
        for r in p["rows"]:
            if r["seed"] in by_seed:
                twice.append(r["seed"])
            by_seed[r["seed"]] = r
    return [by_seed[k] for k in sorted(by_seed)], sorted(set(twice))


def load_study(root: Path = HERE) -> dict[str, dict]:
    """The hedge study's committed per-seed files, merged per transition: rows, truth_info, pred_info and the
    in-plane oracle factor (each part's, all of which must agree)."""
    out = {}
    for key, d in STUDY_RUNS.items():
        parts = [json.loads(p.read_text()) for p in sorted((root / d).glob("*/heart-scope.json"))]
        rows, _ = merge_rows(parts)
        out[key] = {"rows": rows, "parts": parts, "truth_info": parts[0]["truth_info"],
                    "pred_info": parts[0]["pred_info"], "oracle_core_xy_factor": parts[0]["oracle_core_xy_factor"],
                    "task": parts[0]["task"]}
    return out


def load_h25(root: Path = HERE) -> dict[str, dict]:
    return h25.load_dir(root / H25_RUN)


# ------------------------------------------------------------------------------------------------ the guards
def void_reasons(parts: dict[str, list[dict]], study: dict[str, dict], committed_h25: dict[str, dict],
                 commit: str | None = None) -> list[str]:
    """Why the reading may not be read; empty when it may."""
    why = []
    every = [p for ps in parts.values() for p in ps]
    for key in TASK_KEYS:
        ps = parts.get(key)
        if not ps:
            why.append(f"missing transition {key}")
            continue
        rows, twice = merge_rows(ps)
        if twice:
            why.append(f"{key}: seeds {twice} occur in more than one part")
        seeds = [r["seed"] for r in rows]
        if seeds != list(SEEDS):
            why.append(f"{key}: seeds {seeds}, not {list(SEEDS)}")
        complete = True
        for r in rows:
            for field in ("scores", "skills", "groups"):
                gone_q = [q for q in TRUTHS if q not in r.get(field, {})]
                if gone_q:
                    why.append(f"{key} seed {r['seed']}: truths missing from {field}: {gone_q}")
                    complete = False
                    continue
                gone = sorted({a for q in TRUTHS for a in ALL if a not in r[field][q]})
                if gone:
                    why.append(f"{key} seed {r['seed']}: arms missing from {field}: {gone}")
                    complete = False
        st = study.get(key)
        if st is None:
            why.append(f"{key}: no committed study run to check the truths against")
            continue
        for p in ps:
            for q in TRUTHS:
                mine, theirs = p["truth_info"].get(q), st["truth_info"].get(q)
                if mine is None or theirs is None or mine["cells"] != theirs["cells"] \
                        or mine["core_share"] != theirs["core_share"]:
                    why.append(f"{key}: truth {q} differs from the study's (cells or core share)")
            for arm, twin in STUDY_TWIN.items():
                mine, theirs = p["pred_info"].get(arm), st["pred_info"].get(twin)
                if mine is None or theirs is None or mine["pool_cells"] != theirs["pool_cells"] \
                        or mine["core_share"] != theirs["core_share"]:
                    why.append(f"{key}: {arm}'s pool differs from the study's {twin} (cells or core share)")
            f, g = p["oracle_factors"]["in_plane"], st["oracle_core_xy_factor"]
            if not abs(f - g) <= FACTOR_RTOL * abs(g):
                why.append(f"{key}: the in-plane oracle factor differs from the study's")
        if key in HARNESS_KEYS and complete and seeds == list(SEEDS):
            ref = committed_h25.get(key)
            if ref is None:
                why.append(f"{key}: no committed row H25 run to check full_iso_half against")
            else:
                why += h25_identity(key, rows, ref)
    commits = {p["provenance"]["commit"] for p in every}
    if len(commits) > 1:
        why.append("the parts come from different commits")
    if not all(p["provenance"]["tree_clean"] for p in every):
        why.append("a part comes from a dirty tree")
    if commit is not None and (not commit or any(not str(p["provenance"]["commit"]).startswith(commit)
                                                 for p in every)):
        why.append(f"a part does not come from the row's commit {commit!r}")
    return why


def h25_identity(key: str, rows: list[dict], committed: dict) -> list[str]:
    """full_iso_half with the truth as released against row H25's committed iso_half (same draw, same rows, same
    expression; only the size factor differs): equal on the expression-only skills in every seed to TOL."""
    theirs = {r["seed"]: r["skills"].get(H25_TWIN) for r in committed["rows"]}
    why = []
    for r in rows:
        t = theirs.get(r["seed"])
        if t is None:
            why.append(f"{key} seed {r['seed']}: row H25's run has no {H25_TWIN}")
            continue
        mine = r["skills"][OWN][CONTROL]
        for m in EXPRESSION_ONLY:
            if abs(mine[m] - t[m]) > TOL:
                why.append(f"{key} seed {r['seed']}: full_iso_half {m} differs from row H25's {H25_TWIN}")
    return why


# ------------------------------------------------------------------------------------------------ the reading
def _stats(d) -> dict:
    d = np.asarray(d, float)
    return {"mean": float(d.mean()), "lo": float(d.min()), "hi": float(d.max()), "ahead": int((d > 0).sum()),
            "behind": int((d < 0).sum())}


def scores(rows: list[dict], q: str, arm: str) -> np.ndarray:
    return np.array([r["scores"][q][arm] for r in rows], float)


def paired(rows: list[dict], q: str, arm: str, base: str, with_parts: bool = False) -> dict:
    """arm minus base on truth q, paired by seed: the board score and, with_parts, every skill and group skill."""
    out = {"score": _stats(scores(rows, q, arm) - scores(rows, q, base))}
    if with_parts:
        for m in METRICS:
            out[m] = _stats([r["skills"][q][arm][m] - r["skills"][q][base][m] for r in rows])
        for g in GROUPS:
            out[g] = _stats([r["groups"][q][arm][g] - r["groups"][q][base][g] for r in rows])
    return out


def agrees(new: float, old: float, near: float = NEAR_ZERO) -> bool:
    """R1's agreement: the same sign, or both within `near` board points of zero."""
    return bool(np.sign(new) == np.sign(old) or (abs(new) <= near and abs(old) <= near))


def carry_over(rows: dict[str, list[dict]], study: dict[str, dict]) -> dict:
    """R1: per transition and truth, each trim minus full_iso_half beside the study's trim minus hold; per trim, the
    count of the 15 cells that agree."""
    cells, count = {}, {}
    for arm in TRIMMED:
        n = 0
        for key in TASK_KEYS:
            for q in TRUTHS:
                new = paired(rows[key], q, arm, CONTROL)["score"]
                old = paired(study[key]["rows"], q, STUDY_TWIN[arm], STUDY_CONTROL)["score"]
                ok = agrees(new["mean"], old["mean"])
                n += ok
                cells.setdefault(key, {}).setdefault(q, {})[arm] = {"new": new, "study": old, "agree": ok}
        count[arm] = {"agree": n, "of": len(TASK_KEYS) * len(TRUTHS)}
    return {"cells": cells, "agreement": count}


def better_of_two(rows: list[dict], q: str, arm: str, s1: str = REFERENCE) -> float:
    """V(arm) on truth q: the mean over seeds of max(score(s1), score(arm)), paired by seed."""
    return float(np.mean(np.maximum(scores(rows, q, s1), scores(rows, q, arm))))


def values(rows: dict[str, list[dict]]) -> dict:
    """R2: V and the plain mean score per transition, truth and arm (s1form_half's V is its own mean)."""
    out = {}
    for key in TASK_KEYS:
        for q in TRUTHS:
            out.setdefault(key, {})[q] = {
                "V": {arm: better_of_two(rows[key], q, arm) for arm in BUILT},
                "mean": {arm: float(scores(rows[key], q, arm).mean()) for arm in BUILT}}
    return out


def expected(vals: dict, what: str = "V") -> dict:
    """The expectation of V (or of the plain mean) under each prior, per transition and arm."""
    return {key: {name: {arm: float(sum(w * vals[key][q][what][arm] for q, w in pr.items())) for arm in BUILT}
                  for name, pr in PRIORS.items()} for key in TASK_KEYS}


def pick(ev: dict, valid: bool) -> dict:
    """The pick rule fixed in row H26: D3 over the 15 prior-by-transition cells (transitions weighted equally), DJ4
    with the two E10.5 transitions merged into one fold first; B if valid and both reach PICK_BAR, else A."""
    diff = {key: {name: ev[key][name][CANDIDATE] - ev[key][name][CONTROL] for name in PRIORS} for key in TASK_KEYS}
    d3 = float(np.mean([diff[key][name] for key in TASK_KEYS for name in PRIORS]))
    folds = [float(np.mean([np.mean([diff[key][name] for key in fold]) for name in PRIORS])) for fold in FOLDS]
    dj4 = float(np.mean(folds))
    b = bool(valid and d3 >= PICK_BAR and dj4 >= PICK_BAR)
    return {"D3": d3, "DJ4": dj4, "folds": dict(zip([" + ".join(f) for f in FOLDS], folds, strict=True)),
            "per_cell": diff, "bar": PICK_BAR,
            "pick": "B (half_iso_half)" if b else "A (full_iso_half, the file in force)",
            "rule": f"B if the run is valid and D3 >= {PICK_BAR} and DJ4 >= {PICK_BAR} board points; otherwise A. "
                    "C is reported and never picked. Used only if the captain leaves the direction to the AI "
                    "assistant; declining keeps A with no further step."}


def groups_table(rows: dict[str, list[dict]]) -> dict:
    """R3: per transition and truth, each trim minus full_iso_half on every group and metric skill."""
    return {key: {q: {arm: paired(rows[key], q, arm, CONTROL, with_parts=True) for arm in TRIMMED} for q in TRUTHS}
            for key in TASK_KEYS}


def reading(parts: dict[str, list[dict]], study: dict[str, dict], committed_h25: dict[str, dict],
            commit: str | None = None) -> dict:
    why = void_reasons(parts, study, committed_h25, commit)
    if why:
        return {"void": True, "why": why, "pick": "A (full_iso_half, the file in force): the reading is void",
                "follows": "nothing: the reading is void (no table is read, no rerun under this row)"}
    rows = {key: merge_rows(parts[key])[0] for key in TASK_KEYS}
    vals = values(rows)
    ev, em = expected(vals, "V"), expected(vals, "mean")
    first = {key: parts[key][0] for key in TASK_KEYS}
    return {
        "void": False, "why": [], "commit": first[TASK_KEYS[0]]["provenance"]["commit"],
        "labels": {ORACLE: "every score here: the size is an oracle from the truth's core, shared by every arm; a "
                           "sensitivity bound, never a file's score (rule 9, J3)"},
        "oracle_factors": {key: first[key]["oracle_factors"] for key in TASK_KEYS},
        "truth_info": {key: first[key]["truth_info"] for key in TASK_KEYS},
        "pred_info": {key: first[key]["pred_info"] for key in TASK_KEYS},
        "floor_vs_study": {key: _floor_distance(rows[key], study[key]["rows"]) for key in TASK_KEYS},
        "R1_carry_over": carry_over(rows, study),
        "R2_values": vals, "R2_expected_V": ev, "R2_expected_mean": em,
        "R3_groups": groups_table(rows),
        "pick": pick(ev, True),
        "follows": "nothing is built and no rule changes; the table and the rule's pick go to the captain",
    }


def _floor_distance(rows: list[dict], study_rows: list[dict]) -> float | None:
    """Reported only: the largest |floor score here - the study's| over seeds and truths (both 50 by construction)."""
    try:
        return float(max(abs(a["scores"][q][harness.FLOOR] - b["scores"][q]["floor"])
                         for a, b in zip(rows, study_rows, strict=True) for q in TRUTHS))
    except (KeyError, ValueError):
        return None


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\nPick: " + res["pick"] + "\n"
    lines = [f"Try H26 (commit {res['commit']}). Every score below is at an {ORACLE} (the truth's core, shared by "
             "every arm): a sensitivity bound, never a file's score.",
             "Oracle factors: " + "; ".join(f"{k} even {v['even']:.4f}, in plane {v['in_plane']:.4f}"
                                             for k, v in res["oracle_factors"].items())]
    p = res["pick"]
    lines.append(f"PICK RULE: D3 {p['D3']:+.3f}, DJ4 {p['DJ4']:+.3f} (bar {p['bar']}): {p['pick']}")
    lines.append("R1 carry-over (trim minus full_iso_half | the study's trim minus hold; mean [min, max] ahead/10):")
    for key in TASK_KEYS:
        for q in TRUTHS:
            c = res["R1_carry_over"]["cells"][key][q]
            lines.append(f"  {key} {q}: " + "; ".join(
                f"{arm} {c[arm]['new']['mean']:+.2f} [{c[arm]['new']['lo']:+.2f}, {c[arm]['new']['hi']:+.2f}] "
                f"{c[arm]['new']['ahead']} | study {c[arm]['study']['mean']:+.2f} {c[arm]['study']['ahead']}"
                f"{'' if c[arm]['agree'] else ' (disagree)'}" for arm in TRIMMED))
    lines.append("  agreement: " + ", ".join(f"{a} {v['agree']} of {v['of']}"
                                             for a, v in res["R1_carry_over"]["agreement"].items()))
    lines.append("R2 expected better of two, E V (plain mean score in brackets), per transition and prior:")
    for key in TASK_KEYS:
        for name in PRIORS:
            lines.append(f"  {key} {name}: " + ", ".join(
                f"{arm} {res['R2_expected_V'][key][name][arm]:.2f} ({res['R2_expected_mean'][key][name][arm]:.2f})"
                for arm in BUILT))
    lines.append("R3 groups (trim minus full_iso_half, mean group skill):")
    for key in TASK_KEYS:
        for q in TRUTHS:
            g = res["R3_groups"][key][q]
            lines.append(f"  {key} {q}: " + "; ".join(
                f"{arm} " + ", ".join(f"{gg} {g[arm][gg]['mean']:+.4f}" for gg in GROUPS) for arm in TRIMMED))
    return "\n".join(lines) + "\n"


def cmd_run(a) -> None:
    if a.task not in TASK_KEYS:
        sys.exit(f"--task must be one of {TASK_KEYS}")
    run_transition(a.task, payoff.parse_seeds(a.seeds), Path(a.out), Path(a.data))


def cmd_read(a) -> None:
    res = reading(load_parts(Path(a.runs)), load_study(), load_h25(), a.commit)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1))
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
    d.add_argument("--out", required=True)
    d.add_argument("--commit", required=True, help="the commit row H26 names (a prefix of the parts' commit)")
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
