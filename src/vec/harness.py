"""Pseudo-task harness (feat-008): paired backtests of candidates on released stages, decided by vec.adoption.

A pseudo-task (configs/pseudo_tasks.yaml) stands in for a test board's geometry using released stages only (the
validation answers included). For each scoring seed s in vec.adoption.SEEDS the harness
  1. draws the calibrated bake protocol's working copies (vec.evaluate, protocol="bake"): the target's 10% copy and
     its truth half, and the reference's 10% copy, which is both the floor and the DE reference. The floor therefore
     scores exactly 50, and the harness checks that it does and that its anchors cover every ranked metric;
  2. computes the floor and ceiling anchors once and scores every candidate against them, with the scorer build
     "main" (docs/PLAN.md section 3: candidates are compared without the pre-fix build's parity coin; "prefix" is for
     calibration only, and selecting on it would be the parity hedge). On the first seed one candidate is also scored
     with freshly computed anchors, and the two scores must agree;
  3. builds every candidate from the pseudo-task's inputs only, with the build seed BUILD_SEED_OFFSET + s. All
     candidates share it (common random numbers), and it is never a scoring stream: numpy's default_rng([s, 0]),
     which vec.t2_build.sample_cells draws rows with, is the same stream as default_rng(s), the bake protocol's
     reference draw, so a candidate built at seed s would contain the floor's own cells (checked 2026-11-01).
Adoption goes through vec.adoption.keep_change in `decide`, over a test board's configured pseudo-tasks (see `decide`
for which ones count). Results are aggregate numbers per seed; registry rows go to results/registry.csv
(vec.registry).

Leakage rules. `check_task`, applied when the config loads and again before a run: the target is never an input;
extrapolation inputs are no later than the reference, which is earlier than the target, whatever their modality or
genotype; interpolation inputs lie outside the "before"-"after" bracket and none on the target's day; knockout tasks
use no file of the target's genotype unless marked `replicate` (always a diagnostic); `forbidden` files are refused.
Candidates get a view of the task without file paths (TaskView) and of the board without truth_cells and ref_cells
(the test boards may not publish them), and they build inside `guarded_reads`, which refuses and logs any read of
data/raw (inputs arrive loaded), runs/ or experiments/*/results/ (where files derived from targets live, such as the
size and count decisions measured on E10.5 and E7.5). External data a candidate uses must be declared
(Candidate.externals, by stage day) and is refused inside the pseudo-task's version of the test's window. A committed
file computed from released stages (a derived input, DERIVED: the pre-registered T1 held-flat list) is read by the
harness, never by the candidate, and handed over only on pseudo-tasks whose own inputs include every stage it was
computed from (Candidate.derived; elsewhere the candidate does not apply).
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable, Mapping
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field, replace
from fnmatch import fnmatch
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import h5py
import numpy as np
import yaml

from vec import t2_build
from vec.adoption import SEEDS, keep_change
from vec.boards import GROUPS, RANKED_METRICS, REPO_ROOT, VALIDATION_BOARDS, Board, load_boards
from vec.core import draw_balanced, herd_rows
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample, t3_response_detail, working_copies
from vec.io import COORD_KEY, check_submission, dense, make_submission, read_stage, write_submission
from vec.prep import load_stage

CONFIG = REPO_ROOT / "configs" / "pseudo_tasks.yaml"
DATA_DIR = REPO_ROOT / "data" / "raw"
BUILD = "main"
PROTOCOL = "bake"
BUILD_SEED_OFFSET = 1_000_000
FLOOR = "floor"
TASKS = ("T1", "T2", "T3")
KINDS = ("extrap", "interp", "ko")
# read guard: never readable while a candidate builds (inputs arrive loaded; results and runs may hold target-derived
# numbers); boards/t1_composition.json publishes released-stage compositions, the targets' included
DENY_DIRS = ("data", "runs", "experiments/*/results", "docs/research")
DENY_FILES = ("boards/t1_composition.json",)

# E-day of every released file a pseudo-task may name (docs/DATA_NOTES.md; E8.25_late is about E8.3, Data Protocol).
STAGE_DAYS = {
    "T1/E8.5_RNA.h5ad": 8.5, "T1/E9.5_RNA.h5ad": 9.5, "T1/E10.5_RNA.h5ad": 10.5,
    "T2/E6.75.h5ad": 6.75, "T2/E7.25.h5ad": 7.25, "T2/E7.5.h5ad": 7.5, "T2/E8.0.h5ad": 8.0,
    "T2/E8.25_late.h5ad": 8.3, "T2/E8.5.h5ad": 8.5, "T2/E8.75.h5ad": 8.75, "T2/E9.5.h5ad": 9.5,
    "T2/E10.5.h5ad": 10.5,
    "T3/E8.75_gata4_ko_1.h5ad": 8.75, "T3/E8.75_gata4_ko_2.h5ad": 8.75, "T3/E9.5_mab21l2_ko.h5ad": 9.5,
}
GENOTYPE = {"T3/E8.75_gata4_ko_1.h5ad": "gata4", "T3/E8.75_gata4_ko_2.h5ad": "gata4",
            "T3/E9.5_mab21l2_ko.h5ad": "mab21l2"}
KNOCKOUTS = frozenset(GENOTYPE)


class LeakageError(ValueError):
    """A pseudo-task, input or read that could reveal a target."""


@dataclass(frozen=True)
class Derived:
    """A committed file computed from released stages that a candidate may declare (Candidate.derived) for pseudo-tasks
    of `task`. The harness reads it before the build (the read guard refuses it to candidates, as it lives under
    experiments/*/results/) and hands it over only where the pseudo-task's own inputs include every stage in
    `sources()`, so it carries nothing the inputs do not."""

    task: str
    path: Path
    sources: Callable[[], tuple[str, ...]]
    load: Callable[[Path], object]


@lru_cache
def _held_flat_sources() -> tuple[str, ...]:
    """The stages the held-flat list's pool rule was computed on, as its producing config names them."""
    cfg = yaml.safe_load((REPO_ROOT / "experiments" / "t1_direction" / "config.yaml").read_text())
    files = tuple(f"T1/{s}_RNA.h5ad" for s in cfg["held_flat"]["pool_rule"]["stages"])
    unknown = [f for f in files if f not in STAGE_DAYS]
    if unknown:
        raise ValueError(f"the held-flat pool rule names unknown stages {unknown}")
    return files


def _held_flat_load(path: Path) -> tuple[str, ...]:
    from vec import t1_scope

    return tuple(t1_scope.held_flat_genes(path))


DERIVED = {
    # the pre-registered list (experiments/t1_direction): genes named in its config (sex-linked, Y-linked, strain)
    # plus its pool rule on E8.5 and E9.5; committed 2026-10-19 (c992e7a), nine days before the E10.5 file reached us
    "t1_held_flat": Derived("T1", REPO_ROOT / "experiments" / "t1_direction" / "results" / "held_flat_genes.csv",
                            _held_flat_sources, _held_flat_load),
}


def derived_ok(name: str, t: PseudoTask) -> bool:
    """Whether derived input `name` may be handed over on pseudo-task t: it belongs to t's task and every stage it was
    computed from is one of t's own inputs."""
    d = DERIVED[name]
    return d.task == t.task and set(d.sources()) <= set(t.inputs.values())


@dataclass(frozen=True)
class PseudoTask:
    key: str
    test_board: str
    task: str
    kind: str
    target: str
    reference: str
    inputs: Mapping[str, str]
    validation_board: str | None = None
    twin_board: str | None = None
    setting: str | None = None
    forbidden: tuple[str, ...] = ()
    adoption: bool = True
    replicate: bool = False
    target_filter: Mapping | None = None
    notes: str = ""

    @property
    def panel_board(self) -> str:
        """The board whose gene panel and limits rule this pseudo-task uses."""
        return self.validation_board or self.twin_board

    @property
    def counted(self) -> bool:
        """Whether the adoption rule counts this pseudo-task (replicate folds and diagnostics never)."""
        return self.adoption and not self.replicate


def stem(path: str) -> str:
    return Path(path).name.removesuffix(".h5ad")


def check_task(t: PseudoTask) -> PseudoTask:
    """Refuse a pseudo-task that breaks a leakage rule or names an unknown file (see the module docstring)."""
    if t.task not in TASKS or t.kind not in KINDS:
        raise ValueError(f"{t.key}: task must be one of {TASKS} and kind one of {KINDS}")
    if (t.task == "T2") != (t.setting is not None):
        raise ValueError(f"{t.key}: a T2 pseudo-task names its setting, and only T2 does")
    if not t.panel_board or (t.validation_board and t.validation_board not in VALIDATION_BOARDS) or (
            t.twin_board and t.twin_board not in VALIDATION_BOARDS):
        raise ValueError(f"{t.key}: validation_board or twin_board must name a validation board")
    if t.target_filter is not None and (t.task != "T1" or t.validation_board):
        raise ValueError(f"{t.key}: a target filter is for T1 pseudo-tasks on a synthetic board")
    files = list(t.inputs.values())
    unknown = [f for f in [t.target, t.reference, *files, *t.forbidden] if f not in STAGE_DAYS]
    if unknown:
        raise ValueError(f"{t.key}: unknown files {unknown}")
    if t.target in files:
        raise LeakageError(f"{t.key}: the target is an input")
    if set(t.forbidden) & set(files):
        raise LeakageError(f"{t.key}: forbidden inputs {sorted(set(t.forbidden) & set(files))}")
    if t.reference not in files:
        raise ValueError(f"{t.key}: the reference must be one of the inputs")
    day, ref_day = STAGE_DAYS[t.target], STAGE_DAYS[t.reference]
    days = {r: STAGE_DAYS[f] for r, f in t.inputs.items()}
    if t.kind == "extrap":
        if not ref_day < day:
            raise LeakageError(f"{t.key}: the reference must be earlier than the target")
        late = sorted(r for r, d in days.items() if d > ref_day)
        if late:
            raise LeakageError(f"{t.key}: extrapolation inputs later than the reference: {late}")
    elif t.kind == "interp":
        if "before" not in t.inputs or "after" not in t.inputs or t.inputs["before"] != t.reference:
            raise ValueError(f"{t.key}: an interpolation task has 'before' (the reference) and 'after' inputs")
        lo, hi = days["before"], days["after"]
        if not lo < day < hi:
            raise LeakageError(f"{t.key}: the target is not inside its bracket")
        inside = sorted(r for r, d in days.items() if lo < d < hi or d == day)
        if inside:
            raise LeakageError(f"{t.key}: interpolation inputs inside the bracket: {inside}")
    else:
        if t.target not in KNOCKOUTS or t.reference in KNOCKOUTS:
            raise ValueError(f"{t.key}: a knockout pseudo-task targets a knockout file against a wild-type reference")
        sibling = sorted(r for r, f in t.inputs.items() if GENOTYPE.get(f) == GENOTYPE[t.target])
        if sibling and not t.replicate:
            raise LeakageError(f"{t.key}: inputs of the target's own genotype ({sibling}) need replicate: true")
    if t.replicate and t.kind != "ko":
        raise ValueError(f"{t.key}: only a knockout task can be a replicate fold")
    return t


def load_tasks(path: str | Path = CONFIG) -> dict[str, list[PseudoTask]]:
    """test board -> its pseudo-tasks, every one checked."""
    raw = yaml.safe_load(Path(path).read_text())
    out, seen = {}, set()
    for test_board, items in raw.items():
        tasks = []
        for it in items:
            it = dict(it)
            it["forbidden"] = tuple(it.get("forbidden") or ())
            it["notes"] = " ".join(str(it.get("notes", "")).split())
            t = check_task(PseudoTask(test_board=test_board, **it))
            if t.key in seen:
                raise ValueError(f"duplicate pseudo-task key {t.key}")
            seen.add(t.key)
            tasks.append(t)
        out[test_board] = tasks
    return out


def pseudo_counts(n_target: int, n_ref: int) -> dict[str, int]:
    """The organisers' cell counts for a board from its target and reference stages (rows as released), as they hold
    on all five validation boards: 10% working copies, the truth one half, min_cells = min(1000, truth_cells),
    max_cells = max(5000, 3 x the larger of truth_cells and ref_cells)."""
    n_w = max(2, int(round(n_target * SUBSAMPLE_FRAC)))
    truth = n_w // 2
    ref = max(2, int(round(n_ref * SUBSAMPLE_FRAC)))
    return {"truth_cells": truth, "ref_cells": ref, "min_cells": min(1000, truth),
            "max_cells": max(5000, 3 * max(truth, ref))}


def check_index(boards: Mapping[str, Board]) -> None:
    leaking = [k for k, b in boards.items() if k not in VALIDATION_BOARDS and b.anchors]
    if leaking:  # test-board anchors encode hidden targets and are never read
        raise LeakageError(f"{len(leaking)} non-validation board(s) in boards/index.json carry anchors")


def task_board(t: PseudoTask, boards: Mapping[str, Board], n_target: int, n_ref: int) -> Board:
    """The board a pseudo-task is scored on, always with computed anchors. A validation geometry keeps its published
    spec, and the files must reproduce its published counts; any other pseudo-task gets a synthetic board with the
    twin's gene panel, its own counts and its own knockout. Nothing is ever taken from a test board."""
    counts = pseudo_counts(n_target, n_ref)
    ko = GENOTYPE.get(t.target)
    if t.validation_board:
        b = boards[t.validation_board]
        off = sorted(k for k, v in counts.items() if getattr(b, k) != v)
        if off:
            raise ValueError(f"{t.key}: the files do not reproduce {b.key}'s published counts ({off})")
        if ko and (b.target_ko or "").lower() != ko:
            raise ValueError(f"{t.key}: {b.key} is the {b.target_ko} board, not {ko}")
        return replace(b, anchors={})
    twin = boards[t.twin_board]
    return replace(twin, key=f"pseudo:{t.key}", split="pseudo", anchors={}, target_ko=ko or twin.target_ko, **counts)


# ------------------------------------------------------------------------------------------------ candidates
@dataclass(frozen=True)
class TaskView:
    """What a candidate may know about its pseudo-task: no file paths, nothing about the target but its day."""

    key: str
    task: str
    kind: str
    setting: str | None
    test_board: str
    roles: tuple[str, ...]
    role_days: Mapping[str, float]
    role_stems: Mapping[str, str]
    reference_role: str
    target_day: float

    @classmethod
    def of(cls, t: PseudoTask) -> TaskView:
        return cls(key=t.key, task=t.task, kind=t.kind, setting=t.setting, test_board=t.test_board,
                   roles=tuple(t.inputs), role_days={r: STAGE_DAYS[f] for r, f in t.inputs.items()},
                   role_stems={r: stem(f) for r, f in t.inputs.items()}, reference_role=base_role(t),
                   target_day=STAGE_DAYS[t.target])


@dataclass
class BuildContext:
    """What a candidate sees: a path-free view of the task, the board without truth_cells and ref_cells, the loaded
    input stages by role (duplicates dropped on T2/T3; never the target), the build seed and the derived inputs it
    declared for this task (name -> content, read by the harness). Candidates must not modify the inputs."""

    task: TaskView
    board: Board
    inputs: Mapping[str, ad.AnnData]
    seed: int
    derived: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class Candidate:
    """A named builder. It applies to a pseudo-task of one of `tasks` when `requires` (default: every role in one of
    `needs`, a tuple of alternative role sets) says so. `tuned_on` lists pseudo-task keys the candidate was tuned or
    selected on, which `decide` does not count for it; `externals` declares external data as (EXTERNAL_DATA id,
    first day, last day), refused inside a pseudo-task's version of the test window; `derived` names DERIVED inputs
    it uses, and it does not apply where one of its task's cannot be handed over (`derived_ok`)."""

    name: str
    tasks: frozenset[str]
    needs: tuple[frozenset[str], ...]
    build: Callable[[BuildContext], Stage] | None
    params: Mapping = field(default_factory=dict)
    requires: Callable[[PseudoTask], bool] | None = None
    tuned_on: frozenset[str] = frozenset()
    externals: tuple[tuple[str, float, float], ...] = ()
    derived: tuple[str, ...] = ()

    def derived_for(self, t: PseudoTask) -> tuple[str, ...]:
        """The derived inputs this candidate uses on pseudo-task t (those of t's task)."""
        return tuple(n for n in self.derived if DERIVED[n].task == t.task)

    def applies(self, t: PseudoTask) -> bool:
        if t.task not in self.tasks:
            return False
        if not all(derived_ok(n, t) for n in self.derived_for(t)):
            return False
        if self.requires is not None:
            return self.requires(t)
        return any(n <= set(t.inputs) for n in self.needs)


def external_window(t: PseudoTask) -> tuple[float, float] | None:
    """The pseudo-task's version of the test's external-data window: extrapolation (reference day - 1, target day +
    1] (the test: after E9.5 up to and including E13.5 for E12.5 from E10.5); interpolation the open bracket."""
    if t.kind == "extrap":
        return STAGE_DAYS[t.reference] - 1.0, STAGE_DAYS[t.target] + 1.0
    if t.kind == "interp":
        return STAGE_DAYS[t.inputs["before"]], STAGE_DAYS[t.inputs["after"]]
    return None


def check_externals(c: Candidate, t: PseudoTask) -> None:
    if not c.externals:
        return
    w = external_window(t)
    if w is None:
        raise LeakageError(f"{c.name} declares external data; no rule covers it on knockout task {t.key}")
    lo, hi = w
    if t.kind == "extrap":  # (lo, hi]: an external stage range [first, last] overlaps it when last > lo, first <= hi
        bad = [e for e in c.externals if e[2] > lo and e[1] <= hi]
    else:  # (lo, hi), open on both ends
        bad = [e for e in c.externals if e[2] > lo and e[1] < hi]
    if bad:
        raise LeakageError(f"{c.name}: external data inside {t.key}'s window {w}: {[e[0] for e in bad]}")


def load_derived(c: Candidate, t: PseudoTask) -> dict[str, object]:
    """The derived inputs candidate c uses on pseudo-task t, read by the harness (outside the read guard). Refuses one
    computed from a stage that is not among t's own inputs (`applies` already leaves such a candidate out)."""
    out = {}
    for n in c.derived_for(t):
        if not derived_ok(n, t):
            raise LeakageError(f"{c.name}: derived input {n} was computed from {list(DERIVED[n].sources())}, not "
                               f"all of them inputs of {t.key}")
        out[n] = DERIVED[n].load(DERIVED[n].path)
    return out


def derived_record(c: Candidate, t: PseudoTask) -> dict[str, dict]:
    """What a result records about the derived inputs c used on t: the file and the stages it was computed from."""
    return {n: {"file": str(DERIVED[n].path.relative_to(REPO_ROOT)), "sources": list(DERIVED[n].sources())}
            for n in c.derived_for(t)}


def base_role(t: PseudoTask) -> str:
    """The input role holding the reference stage's own cells."""
    return next(role for role, f in t.inputs.items() if f == t.reference)


def count_for(board: Board, pool_size: int) -> int:
    """The approved count rule (max_cells, within the board's limits) on any board."""
    return t2_build.n_cells(board, pool_size)[0]


def draw_copy(ctx: BuildContext, herd: bool = False) -> tuple[ad.AnnData, np.ndarray, np.ndarray, np.ndarray | None]:
    """Real reference cells at the board's count: (pool, rows, X float32 in board gene order, coords or None). herd
    (coordinate boards; try H4 of experiments/harness/TRIES.md): the draw herded towards the pool's mean panel
    expression (`t2_build.sample_cells(..., herd_expr=...)`)."""
    pool = ctx.inputs[ctx.task.reference_role]
    n = count_for(ctx.board, pool.n_obs)
    if ctx.board.needs_coords:
        smp = t2_build.sample_cells(pool, n, seed=ctx.seed,
                                    herd_expr=pool[:, ctx.board.genes].X if herd else None)
        rows, coords = smp.rows, smp.coords
    else:
        rows, coords = draw_balanced(pool.n_obs, n, np.random.default_rng([ctx.seed, 0])), None
    X = dense(pool[:, ctx.board.genes].X[rows]).astype(np.float32)
    return pool, rows, X, coords


def finish(ctx: BuildContext, pool: ad.AnnData, X: np.ndarray, coords: np.ndarray | None,
           guard_source: np.ndarray | None = None, guard_weights: np.ndarray | None = None) -> Stage:
    """A candidate as the scorer reads its file, after the checks the server applies: written and checked with
    vec.io (T1), or through the T2 head's cloud guard and checked builder (coordinate boards). A failure raises.
    guard_source / guard_weights: the cloud the cells were drawn from, in the submission's frame, and the draw weights,
    when that is not `pool` as loaded (a size step or an interpolation moved it); the guard then checks the cells
    against it as it checks a copy against its pool."""
    X = np.clip(np.asarray(X, dtype=np.float32), 0.0, None)
    if not ctx.board.needs_coords:
        sub = make_submission(X, ctx.board.genes)
        with tempfile.TemporaryDirectory() as tmp:
            rep = check_submission(write_submission(sub, Path(tmp) / "submission.h5ad", stamp=False), ctx.board)
        if not rep.ok:
            raise ValueError(f"{ctx.task.key}: the built file fails its check: {rep.errors}")
        return Stage.from_anndata(sub, ctx.board.genes)
    source = np.asarray(pool.obsm[COORD_KEY], dtype=float)[:, :3] if guard_source is None else guard_source
    guard = t2_build.cloud_guard(coords, source, seed=ctx.seed, source_weights=guard_weights)
    if not guard.ok:
        raise ValueError(f"{ctx.task.key}: cloud guard flags {guard.flags}")
    sub, rep = t2_build.build_submission(ctx.board, X, coords, genes=list(ctx.board.genes))
    if not rep.ok:
        raise ValueError(f"{ctx.task.key}: the built file fails its check: {rep.errors}")
    return Stage.from_anndata(sub, ctx.board.genes, need_coords=True)


def _copy(ctx: BuildContext) -> Stage:
    pool, _, X, coords = draw_copy(ctx)
    return finish(ctx, pool, X, coords)


def _type_means(a: ad.AnnData, genes) -> dict[str, np.ndarray]:
    X, ct = a[:, genes].X, np.asarray(a.obs["celltype"]).astype(str)
    return {c: np.asarray(X[ct == c].mean(0)).ravel() for c in np.unique(ct)}


def start_role(t: PseudoTask | TaskView) -> str | None:
    """The stage before the reference, as the organisers' pseudobulk_shift uses it: `start` on extrapolation,
    `earlier` (before the earlier bracket) on interpolation."""
    roles = t.roles if isinstance(t, TaskView) else tuple(t.inputs)
    if t.kind == "extrap" and "start" in roles:
        return "start"
    if t.kind == "interp" and "earlier" in roles:
        return "earlier"
    return None


def _pseudobulk_shift(ctx: BuildContext) -> Stage:
    """The organisers' pseudobulk_shift on real reference cells at the board's count: each cell moves by its type's
    mean change from the stage before to the reference stage (raw celltype labels; types absent before stay put)."""
    pool, rows, X, coords = draw_copy(ctx)
    return finish(ctx, pool, _type_shift(ctx, pool, rows, X), coords)


def _type_shift(ctx: BuildContext, pool: ad.AnnData, rows: np.ndarray, X: np.ndarray) -> np.ndarray:
    """pseudobulk_shift's expression step: each drawn cell plus its raw type's mean change from the stage before."""
    ref_m = _type_means(pool, ctx.board.genes)
    start_m = _type_means(ctx.inputs[start_role(ctx.task)], ctx.board.genes)
    ct = np.asarray(pool.obs["celltype"]).astype(str)[rows]
    for c in np.unique(ct):
        if c in start_m:
            X[ct == c] += (ref_m[c] - start_m[c]).astype(np.float32)
    return X


def _shift_transfer(ctx: BuildContext) -> Stage:
    """The organisers' shift_transfer on real wild-type cells at the board's count: the training knockout's pseudobulk
    change against its matched wild type, added to every cell (both replicates pooled when a task gives two)."""
    pool, _, X, coords = draw_copy(ctx)
    g = ctx.board.genes
    kos = [ctx.inputs[r][:, g].X for r in ("ko_train", "ko_train_2") if r in ctx.inputs]
    ko_mean = np.asarray(np.vstack([dense(k) for k in kos]).mean(0)).ravel()
    wt_mean = np.asarray(ctx.inputs["wt_train"][:, g].X.mean(0)).ravel()
    return finish(ctx, pool, X + (ko_mean - wt_mean).astype(np.float32), coords)


# ------------------------------------------------------------------------------------------------ B0 (feat-009)
# Fixed on 2026-11-01, before any B0 result: the damping of B0's momentum shift, a small change of the right sign.
B0_DAMPING = 0.5


def _shallow(a: ad.AnnData, **obs_cols) -> ad.AnnData:
    """The same cells and matrices with a copied obs (plus columns), so a head may write obs without touching the
    inputs; the matrices are shared, not copied."""
    obs = a.obs.copy()
    for k, v in obs_cols.items():
        obs[k] = v
    return ad.AnnData(X=a.X, obs=obs, var=a.var, layers=dict(a.layers), obsm=dict(a.obsm))


def t1_chamber_shares(ref: ad.AnnData, source: str) -> dict[str, float]:
    """B0's T1 composition: the reference stage's own shares, restricted to the chamber categories of
    configs/t1_chamber_classes.yaml and renormalised (the E12.5 test is the heart chamber only; captain, 2026-11-01)."""
    from vec import t1_scope

    spec = t1_scope.load_spec()
    a = _shallow(ref)
    t1_scope.classify(a, source=source, spec=spec)
    comp = t1_scope.composition(a)
    keep = {c: v for c, v in comp.items() if c in spec.in_scope("chamber") and v > 0}
    total = sum(keep.values())
    return {c: v / total for c, v in keep.items()}


def _b0_t1(damping: float, herd: bool = False) -> Callable[[BuildContext], Stage]:
    def build(ctx: BuildContext) -> Stage:
        from vec import t1_scope

        role = ctx.task.reference_role
        ref, src = _shallow(ctx.inputs[role]), ctx.task.role_stems[role]
        shares = t1_chamber_shares(ref, src)
        start = _shallow(ctx.inputs["start"]) if damping else None
        held = list(ctx.derived["t1_held_flat"]) if damping else None  # handed over: the guard refuses the file
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, count_for(ctx.board, ref.n_obs), ctx.seed,
                                           start=start, damping=damping, held_flat=held, ref_source=src,
                                           start_source=ctx.task.role_stems["start"] if damping else None,
                                           herd=herd)
        return finish(ctx, ref, dense(sub.X), None)
    return build


def _copy_herd(ctx: BuildContext) -> Stage:
    """A copy of the reference with its draw herded towards the pool's mean expression (T1 try T6, a diagnostic of how
    much of a copy's de_score handicap is sampling noise; experiments/harness/TRIES.md)."""
    pool, rows, _, coords = draw_copy(ctx)
    g = ctx.board.genes
    rows, _ = herd_rows(pool[:, g].X, np.arange(pool.n_obs), rows, np.random.default_rng([ctx.seed, 9]))
    return finish(ctx, pool, dense(pool[:, g].X[rows]).astype(np.float32), coords)


def _b0_t2(ctx: BuildContext, depth: bool = False, herd: bool = False) -> Stage:
    """B0 on a T2 board: real reference cells at the count rule, moved in count space by the per-type change of the
    harmonised lineage level (obs["merged"]; groups without a stable cross-stage meaning held flat). Extrapolation:
    from `start` to the reference, damped by B0_DAMPING. Interpolation: from the earlier bracket towards the later one,
    by the target's fraction of the bracket (its day is part of the geometry). Coordinates and size stay the source's
    (the heart size cannot be backtested; the test build applies the pre-registered size). depth (extrapolation only;
    heart try H3 of experiments/harness/TRIES.md): the moved cells are then binomially thinned to the target depth that
    vec.t2_build.DEPTH_EXTRAP_RULE extrapolates at B0_DAMPING from the start's and the reference's median depths.
    herd (try H4): the copy's draw herded towards the pool's mean panel expression (`draw_copy`)."""
    from vec import t1_scope

    pool, rows, _, coords = draw_copy(ctx, herd=herd)
    t = ctx.task
    if t.kind == "extrap":
        start, end, damping = ctx.inputs["start"], pool, B0_DAMPING
    else:
        lo, hi = t.role_days["before"], t.role_days["after"]
        start, end, damping = pool, ctx.inputs["after"], (t.target_day - lo) / (hi - lo)

    def typed(a: ad.AnnData, tag: str) -> ad.AnnData:
        m = np.asarray(a.obs["merged"]).astype(str)
        return _shallow(a, b0_type=np.where(m == "Untracked", f"Untracked@{tag}", m))

    like = typed(pool, "start" if t.kind == "interp" else "end")[rows]
    shifted, _ = t1_scope.momentum_shift(typed(start, "start"), typed(end, "end"), like, damping, (), key="b0_type",
                                         seed=[ctx.seed, 3])
    g = ctx.board.genes
    if depth:
        if t.kind != "extrap":
            raise ValueError("the depth step of B0 is for extrapolation boards")
        d_start, d_ref = t2_build.median_depth(start, g), t2_build.median_depth(pool, g)
        target = t2_build.extrapolated_depth(d_start, d_ref, t.role_days["start"], t.role_days[t.reference_role],
                                             t.target_day, B0_DAMPING)
        X, _ = t2_build.thin_to_depth(shifted, g, d_ref, target, [ctx.seed, 8])
        return finish(ctx, pool, X, coords)
    return finish(ctx, pool, dense(shifted[:, g].X), coords)


# ------------------------------------------------------------------------------------------------ T2 embryo (feat-010)
def target_rms(ctx: BuildContext) -> float:
    """The pre-registered embryo size rule (raw mode) on the task's own input stages: PCHIP of log RMS radius against
    the stage day, read at the target's day (vec.t2_build.pchip_rms)."""
    days = [ctx.task.role_days[r] for r in ctx.inputs]
    rms = [t2_build.rms_radius(np.asarray(a.obsm[COORD_KEY], dtype=float)[:, :3]) for a in ctx.inputs.values()]
    return t2_build.pchip_rms(days, rms, ctx.task.target_day)


def _psb_size(ctx: BuildContext) -> Stage:
    """pseudobulk_shift (the embryo base since 2026-11-01) with the cloud scaled about the pool's centroid to the
    embryo size rule read on the task's own stages: the size step alone."""
    pool, rows, X, coords = draw_copy(ctx)
    X = _type_shift(ctx, pool, rows, X)
    P = np.asarray(pool.obsm[COORD_KEY], dtype=float)[:, :3]
    c0, f = P.mean(0), target_rms(ctx) / t2_build.rms_radius(P)
    return finish(ctx, pool, X, c0 + (coords - c0) * f, guard_source=c0 + (P - c0) * f)


def _fgw(ctx: BuildContext, t_pos: float | None = None, expression: str = "mixture", herd: bool = False,
         move_coupling: str = t2_build.MOVE_COUPLINGS[0]) -> Stage:
    """The T2 embryo pipeline (docs/PLAN.md 4.2): displacement interpolation between the brackets at the target's
    time fraction (vec.t2_build.interpolate_clouds: fused Gromov-Wasserstein coupling, proper rotation, McCann moves,
    a mixture of real cells of both brackets), scaled to the embryo size rule read on the task's own stages. t_pos:
    the s2 shape variants' position fraction (0: the earlier bracket's geometry, 1: the later's; default: the
    mixture's). expression: "mixture" (the cells as released), "depth" (vec.t2_build.depth_matched_expression) or
    "blend" (vec.t2_build.blend_expression), tries E7 and E8 of experiments/harness/TRIES.md, or "latent" (try E14:
    vec.t2_latent, the task's negative-binomial VAE trained once on both brackets' panel counts). move_coupling: "fgw"
    (the entry) or "aligned" (try E10: the moves from vec.core.aligned_coupling in the rotated common frame)."""
    t = ctx.task
    before, after = ctx.inputs["before"], ctx.inputs["after"]
    lo, hi = t.role_days["before"], t.role_days["after"]
    g = ctx.board.genes
    frac = (t.target_day - lo) / (hi - lo)
    Eb, Ea = dense(before[:, g].X).astype(np.float32), dense(after[:, g].X).astype(np.float32)
    res = t2_build.interpolate_clouds(np.asarray(before.obsm[COORD_KEY], dtype=float)[:, :3],
                                      np.asarray(after.obsm[COORD_KEY], dtype=float)[:, :3], Eb, Ea,
                                      frac, count_for(ctx.board, before.n_obs), ctx.seed,
                                      target_rms=target_rms(ctx), t_pos=t_pos,
                                      with_partners=expression in ("blend", "latent"),
                                      herd=herd, move_coupling=move_coupling, with_coupling=expression == "latent")
    if expression == "mixture":
        X = np.empty((len(res.rows), len(g)), dtype=np.float32)
        X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    elif expression == "depth":
        X, _ = t2_build.depth_matched_expression(before, after, res, g, frac, [ctx.seed, 7])
    elif expression == "blend":
        X, _ = t2_build.blend_expression(before, after, res, g, frac, [ctx.seed, 7])
    elif expression == "latent":  # try E14: t2_latent loads torch and scvi-tools only when it trains
        from vec import t2_latent

        X, _ = t2_latent.latent_midpoint_counts(t2_latent.task_model(before, after, g), before, after, res, g, frac,
                                                [ctx.seed, t2_latent.EXPRESSION_STREAM])
    else:
        raise ValueError(f"unknown expression step {expression!r}")
    return finish(ctx, before, X, res.coords, guard_source=res.pool_coords, guard_weights=res.pool_weights)


def _fgw_shape0(ctx: BuildContext) -> Stage:
    return _fgw(ctx, t_pos=0.0)


def _fgw_shape1(ctx: BuildContext) -> Stage:
    return _fgw(ctx, t_pos=1.0)


def _fgw_depth(ctx: BuildContext) -> Stage:
    return _fgw(ctx, expression="depth")


def _b0_depth(ctx: BuildContext) -> Stage:
    return _b0_t2(ctx, depth=True)


def _b0_herd(ctx: BuildContext) -> Stage:
    return _b0_t1(B0_DAMPING, herd=True)(ctx) if ctx.task.task == "T1" else _b0_t2(ctx, herd=True)


def _fgw_herd(ctx: BuildContext) -> Stage:
    return _fgw(ctx, herd=True)


def _t2_extrap_with_start(t: PseudoTask) -> bool:
    return t.task == "T2" and t.kind == "extrap" and "start" in t.inputs


def _fgw_blend(ctx: BuildContext) -> Stage:
    return _fgw(ctx, expression="blend")


def _fgw_aligned(ctx: BuildContext) -> Stage:
    return _fgw(ctx, move_coupling="aligned")


def _fgw_latent(ctx: BuildContext) -> Stage:
    return _fgw(ctx, expression="latent")


def _interp(t: PseudoTask) -> bool:
    return t.task == "T2" and t.kind == "interp" and {"before", "after"} <= set(t.inputs)


def _b0(ctx: BuildContext) -> Stage:
    return _b0_t1(B0_DAMPING)(ctx) if ctx.task.task == "T1" else _b0_t2(ctx)


def _b0_applies(t: PseudoTask) -> bool:
    if t.task == "T1" or (t.task == "T2" and t.kind == "extrap"):
        return "start" in t.inputs
    return t.task == "T2" and t.kind == "interp" and {"before", "after"} <= set(t.inputs)


# ------------------------------------------------------------------------------------------------ T3 tries K2a and K2b
# Fixed in their experiments/harness/TRIES.md rows (registered 2026-11-04 at 945de7d, before any run): s1's lead
# component, the pooled Gata4 factor of the pre-registered T3 run, on the Gata4 -> Mab21l2 fold's WT E9.5 (K2a,
# t3_factor_xko), and the same factor permuted across genes within deciles of WT E8.75 mean expression (K2b,
# t3_factor_xko_perm, the control). Built with the T3 head's own code, experiments/t3_decomp (t3_decompose.py is
# hash-pinned by build_submission.py; rules.yaml holds every setting the factor and the build read).
T3_DIR = REPO_ROOT / "experiments" / "t3_decomp"
K2_FACTOR_FILES = ("T3/E8.75_gata4_ko_1.h5ad", "T3/E8.75_gata4_ko_2.h5ad", "T2/E8.75.h5ad")  # ko_train(_2), wt_train
K2_FACTOR_PRESET = "gata4"  # rules.yaml preset of s1's factor: the replicates against WT E8.75, its label mode
K2_REFERENCE_PRESET = "mab21l2"  # rules.yaml preset whose reference is WT E9.5, the fold's WT: its label mode
K2_PERM_SEED = 20261104  # K2b: one permutation, this seed (its row)
K2_PERM_BINS = 10  # K2b: deciles
K2_STREAM = "t3_factor_xko"  # both arms build from this stream: the same cells drawn, then each arm's factor applied
_K2_PREP: dict = {}


def _t3_head():
    """experiments/t3_decomp/build_submission.py, and through it (`.t3d`) the pre-registered t3_decompose.py, imported
    as tests/test_hidden_counts.py imports them."""
    if str(T3_DIR) not in sys.path:
        sys.path.insert(0, str(T3_DIR))
    import build_submission

    return build_submission


def _t3_vocab(t3d, a: ad.AnnData, rules: Mapping) -> tuple[object, np.ndarray]:
    """(vocabulary, each cell's label index) of a pool labelled by t3_decompose.add_label_column, as
    build_submission.load_reference makes them (every label a class)."""
    labels = np.array(sorted(set(a.obs["t3_label"].astype(str))))
    V = t3d.make_vocab(a, "t3_label", SimpleNamespace(classes=labels), rules)
    pos = {lab: i for i, lab in enumerate(V.labels)}
    return V, np.array([pos[lab] for lab in a.obs["t3_label"].astype(str)])


def t3_gata4_factor(wt: ad.AnnData, reps: list[ad.AnnData], rules: Mapping, n_jobs: int = 1) -> np.ndarray:
    """Try K2a's factor, s1's lead component: the pooled non-Mesp1 Gata4 factor as the pre-registered T3 run computes
    it for s1 (t3_decompose.build_context, then submission_parameters, which writes the factor file
    build_submission.py loads for s1). WT E8.75's own labels in the gata4 preset's label mode; each replicate's cells
    (littermate wild types split off) called into that vocabulary by the frozen harmoniser, depth-matched as there;
    both replicates' cells pooled; t3_decompose.gene_factor over the strict non-Mesp1 labels with at least min_cells
    cells on both sides, weighted by the WT shares (rules.yaml components.factor: eps 1e-5, clip 0.25-4, 20 cells).
    wt and reps as the harness loads them (vec.prep.load_stage(drop_duplicates=True)), replicates in file order; none
    of them is modified."""
    t3d = _t3_head().t3d
    hp, fc = rules["harmonisation"], rules["components"]["factor"]
    w = _shallow(wt)
    t3d.add_label_column(w, rules["presets"][K2_FACTOR_PRESET]["label_mode"])
    V, wlab = _t3_vocab(t3d, w, rules)
    Cw = t3d.dense_counts(w)
    harm = t3d.Harmoniser(w, "t3_label", t3d.label_groups(w, "t3_label"), hp, n_jobs)
    Cs, labs = [], []
    for a in reps:
        a = _shallow(a)
        littermates, _, _ = t3d.littermate_split(a)
        if littermates.any():
            a = a[~littermates].copy()
        calls, _ = harm.calls(a, seed=hp["seed"])
        Cs.append(t3d.dense_counts(a))
        labs.append(t3d.calls_to_index(V, calls))
    C, lab = np.vstack(Cs), np.concatenate(labs)
    ok = ((np.bincount(lab[lab >= 0], minlength=V.L) >= fc["min_cells"])
          & (np.bincount(wlab, minlength=V.L) >= fc["min_cells"]))
    phi, _ = t3d.gene_factor(t3d.props_type_means(C, lab, V.L), t3d.props_type_means(Cw, wlab, V.L), V.p, ok & V.nm,
                             fc["eps"], tuple(fc["clip"]))
    return phi


def t3_mean_expression(a: ad.AnnData) -> np.ndarray:
    """K2b's decile basis: each gene's mean released expression (log1p CP10k, `.X`) over a pool's cells, in float64."""
    return dense(a.X).astype(np.float64).mean(0)


def t3_decile_permutation(phi, basis, seed: int = K2_PERM_SEED, bins: int = K2_PERM_BINS) -> np.ndarray:
    """Try K2b's factor: phi permuted across genes within bins of `basis`. Genes ranked by basis (a stable sort, so
    ties keep panel order) are cut into `bins` runs of consecutive ranks (sizes differing by at most one); one
    generator, default_rng(seed), permutes each run's values in turn, from the lowest basis up. Every gene gets a value
    of its own bin, and each bin keeps its values."""
    phi = np.asarray(phi, dtype=np.float64)
    rng = np.random.default_rng(seed)
    out = phi.copy()
    for genes in np.array_split(np.argsort(np.asarray(basis, dtype=np.float64), kind="stable"), bins):
        out[genes] = phi[genes[rng.permutation(len(genes))]]
    return out


def t3_reference(pool: ad.AnnData, preset: str, rules: Mapping):
    """A WT pool in t3_decompose's form, as build_submission.load_reference prepares a file without stage mixing: the
    preset's label mode, the vocabulary, int counts and the kNN-smoothed proportions that are the Poisson top-up
    profile (rules.yaml components.smoothing, note 5), from a pool the harness loaded."""
    bs = _t3_head()
    w = _shallow(pool)
    bs.t3d.add_label_column(w, rules["presets"][preset]["label_mode"])
    V, wlab = _t3_vocab(bs.t3d, w, rules)
    Cw = bs.t3d.dense_counts(w)
    sm = rules["components"]["smoothing"]
    return bs.Reference(preset=preset, rules=rules, genes=list(w.var_names), V=V,
                        pool=bs.t3d.Pool(Cw, wlab, bs.t3d.smoothed_props(Cw, sm["n_pcs"], sm["k"]), V.L),
                        coords=np.asarray(w.obsm[COORD_KEY], dtype=np.float64)[:, :3])


def _k2_prepared(ctx: BuildContext) -> dict:
    """K2a's factor, K2b's permuted factor and the fold's WT in t3_decompose's form. None depends on the seed, and the
    harness hands every seed the same input objects, so they are made once per loaded fold."""
    ref_role = ctx.task.reference_role
    objs = tuple(ctx.inputs[r] for r in ("wt_train", "ko_train", "ko_train_2", ref_role))
    hit = _K2_PREP.get("fold")
    if hit is None or any(a is not b for a, b in zip(hit["inputs"], objs, strict=True)):
        rules = _t3_head().load_rules()
        wt, pool = ctx.inputs["wt_train"], ctx.inputs[ref_role]
        if list(wt.var_names) != list(pool.var_names):
            raise ValueError(f"{ctx.task.key}: the factor's WT and the reference hold different genes")
        phi = t3_gata4_factor(wt, [ctx.inputs["ko_train"], ctx.inputs["ko_train_2"]], rules)
        _K2_PREP.clear()
        hit = _K2_PREP["fold"] = {"inputs": objs, "phi": phi,
                                  "phi_perm": t3_decile_permutation(phi, t3_mean_expression(wt)),
                                  "ref": t3_reference(pool, K2_REFERENCE_PRESET, rules)}
    return hit


def _t3_factor_xko(ctx: BuildContext, permuted: bool = False) -> Stage:
    """Tries K2a (t3_factor_xko) and K2b (permuted: t3_factor_xko_perm): the factor, at exponent 1, on the reference
    WT's real cells at its own composition, built as the T3 s1 build builds: build_submission.construct, i.e.
    t3_decompose.build (alloc and draw_rows per label, then apply_gene_factor with the kNN-smoothed top-up profile; no
    depth ratio, no composition term, no stage mixing), then build_submission.assemble (panel columns, random row
    order, each cell at its source cell's position). Never vec.core.apply_factors, whose top-up is from each cell's
    own counts."""
    bs = _t3_head()
    prep = _k2_prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    built = bs.construct(ref, count_for(ctx.board, pool.n_obs), bs.t3d.stable_rng(ctx.seed, K2_STREAM), f=ref.V.p,
                         phi=prep["phi_perm" if permuted else "phi"])
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{K2_STREAM}:assemble"))
    return finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def _t3_factor_xko_perm(ctx: BuildContext) -> Stage:
    return _t3_factor_xko(ctx, permuted=True)


def _k2_applies(t: PseudoTask) -> bool:
    """K2a and K2b run where the factor's own inputs are the task's: both Gata4 replicates and their matched WT E8.75
    (T3/Gata4->Mab21l2 in configs/pseudo_tasks.yaml)."""
    return t.task == "T3" and tuple(t.inputs.get(r) for r in ("ko_train", "ko_train_2", "wt_train")) == K2_FACTOR_FILES


K2_PARAMS = {"factor": "s1's lead component: Gata4 replicates pooled against WT E8.75, strict non-Mesp1 labels "
                       "(t3_decompose harmoniser calls and gene_factor; rules.yaml components.factor), exponent 1",
             "build": "build_submission.construct and assemble (t3_decompose.build: apply_gene_factor with the "
                      "kNN-smoothed top-up profile of rules.yaml note 5); no depth, no composition",
             "reference_labels": f"rules.yaml preset {K2_REFERENCE_PRESET}", "stream": K2_STREAM}

_ANY = (frozenset(),)
CANDIDATES = {
    FLOOR: Candidate(FLOOR, frozenset(TASKS), _ANY, None,
                     {"what": "the organisers' floor: the reference's 10% working copy (copy_last or wt_identity)",
                      "scope": "whole"}),
    "copy": Candidate("copy", frozenset(TASKS), _ANY, _copy,
                      {"what": "real reference cells at the count rule", "scope": "whole"}),
    "pseudobulk_shift": Candidate("pseudobulk_shift", frozenset({"T1", "T2"}), (), _pseudobulk_shift,
                                  {"what": "per-type mean change from the stage before, on copied reference cells",
                                   "scope": "whole"},
                                  requires=lambda t: start_role(t) is not None),
    "shift_transfer": Candidate("shift_transfer", frozenset({"T3"}), (frozenset({"ko_train", "wt_train"}),),
                                _shift_transfer, {"what": "training knockout's pseudobulk change on copied WT cells"}),
    "b0": Candidate("b0", frozenset({"T1", "T2"}), (), _b0,
                    {"what": "B0 (feat-009): T1 chamber-only copy plus damped momentum; T2 copy plus damped momentum "
                             "or bracket interpolation, both in count space", "damping": B0_DAMPING,
                     "scope": "chamber"}, requires=_b0_applies, derived=("t1_held_flat",)),
    "b0_scope": Candidate("b0_scope", frozenset({"T1"}), _ANY, _b0_t1(0.0),
                          {"what": "T1 chamber-only copy without a shift (B0's composition alone)",
                           "scope": "chamber"}),
    "psb_size": Candidate("psb_size", frozenset({"T2"}), (), _psb_size,
                          {"what": "pseudobulk_shift with the cloud scaled to the embryo size rule (PCHIP of log RMS "
                                   "radius through the task's own stages, read at the target's day)", "scope": "whole"},
                          requires=lambda t: _interp(t) and start_role(t) is not None),
    "fgw": Candidate("fgw", frozenset({"T2"}), (), _fgw,
                     {"what": "displacement interpolation between the brackets (fused Gromov-Wasserstein coupling, "
                              "proper rotation, McCann moves) with a mixture of real cells of both, scaled to the "
                              "embryo size rule", "alpha": t2_build.INTERP_ALPHA, "n_anchor": t2_build.INTERP_N_ANCHOR,
                      "n_pc": t2_build.INTERP_N_PC, "k": t2_build.INTERP_K, "scope": "whole"},
                     requires=_interp),
    # the embryo s2 shape variants (docs/S2_RULES.md): fgw's coupling and mixture, positions from one bracket's geometry
    "fgw_shape0": Candidate("fgw_shape0", frozenset({"T2"}), (), _fgw_shape0,
                            {"what": "fgw's mixture of real cells in the earlier bracket's geometry (t_pos = 0), "
                                     "scaled to the embryo size rule", "scope": "whole"}, requires=_interp),
    "fgw_shape1": Candidate("fgw_shape1", frozenset({"T2"}), (), _fgw_shape1,
                            {"what": "fgw's mixture of real cells in the later bracket's geometry (t_pos = 1), "
                                     "scaled to the embryo size rule", "scope": "whole"}, requires=_interp),
    # tries E7 and E8 (experiments/harness/TRIES.md, registered before their run): the target's detection depth
    "fgw_depth": Candidate("fgw_depth", frozenset({"T2"}), (), _fgw_depth,
                           {"what": "fgw with each bracket's cells binomially thinned to the interpolated stage's "
                                    "predicted depth (log-linear between the brackets' median counts per cell; a "
                                    "shallower bracket kept as released)", "rule": t2_build.DEPTH_RULE,
                            "scope": "whole"}, requires=_interp),
    "fgw_blend": Candidate("fgw_blend", frozenset({"T2"}), (), _fgw_blend,
                           {"what": "fgw's cells as binomial blends of each cell and a coupled partner of the other "
                                    "bracket, splitting the expected counts (1 - t) : t at the predicted depth",
                            "rule": "blend: " + t2_build.DEPTH_RULE, "scope": "whole"}, requires=_interp),
    # T1 tries T5 and T6 (experiments/harness/TRIES.md, registered before their run): herded draws
    "b0_herd": Candidate("b0_herd", frozenset({"T1", "T2"}), (), _b0_herd,
                         {"what": "B0 with its draw herded towards the pool's mean expression (vec.core.herd_rows): "
                                  "on T1 each category's draw, on T2 the copy's draw before B0's shift",
                          "damping": B0_DAMPING, "scope": "chamber"},
                         requires=_b0_applies, derived=("t1_held_flat",)),
    "fgw_herd": Candidate("fgw_herd", frozenset({"T2"}), (), _fgw_herd,
                          {"what": "fgw with each bracket's draw herded towards that bracket's mean panel expression "
                                   "(embryo try E9)", "scope": "whole"}, requires=_interp),
    # embryo try E10 (experiments/harness/TRIES.md, registered before its run; it runs only if its target-free gate
    # passes): FGW's rotation, the moves from an exact transport of the same anchors in the common frame
    "fgw_aligned": Candidate("fgw_aligned", frozenset({"T2"}), (), _fgw_aligned,
                             {"what": "fgw with the moves taken from an exact transport of the same anchors on a fused "
                                      "cost of expression and position in the rotated common frame "
                                      "(vec.core.aligned_coupling; FGW's rotation kept)",
                              "alpha": t2_build.INTERP_ALPHA, "scope": "whole"}, requires=_interp),
    # embryo try E14 (experiments/harness/TRIES.md, registered before its run; it runs only if its target-free fidelity
    # gate N18 passes): fgw's draws, partners and positions, expression from a count-model latent interpolation
    "fgw_latent": Candidate("fgw_latent", frozenset({"T2"}), (), _fgw_latent,
                            {"what": "fgw's draws, partners and positions with expression from a negative-binomial "
                                     "VAE (scvi-tools, n_latent 20, gene dispersion, observed library, no batch key, "
                                     "trained once per task on both brackets' panel counts): each coupled pair's "
                                     "latent means interpolated with a variance-preserving term, counts sampled "
                                     "(never a decoded mean) at the drawn real cell's panel library (vec.t2_latent)",
                             "scope": "whole"}, requires=_interp),
    "copy_herd": Candidate("copy_herd", frozenset({"T1"}), _ANY, _copy_herd,
                           {"what": "a copy of the reference with its draw herded towards the pool's mean expression "
                                    "(a diagnostic)", "scope": "whole"}),
    # heart try H3 (experiments/harness/TRIES.md, registered before its run)
    "b0_depth": Candidate("b0_depth", frozenset({"T2"}), (), _b0_depth,
                          {"what": "B0 on a T2 extrapolation board with its cells binomially thinned to the target's "
                                   "depth, extrapolated at B0's damping (half the latest log rate per day over the "
                                   "horizon)", "damping": B0_DAMPING, "rule": t2_build.DEPTH_EXTRAP_RULE,
                           "scope": "whole"}, requires=_t2_extrap_with_start),
    # T3 tries K2a and K2b (experiments/harness/TRIES.md, registered before their run): s1's Gata4 factor on the
    # Gata4 -> Mab21l2 fold, and its permuted control
    "t3_factor_xko": Candidate("t3_factor_xko", frozenset({"T3"}), (), _t3_factor_xko,
                               {"what": "try K2a: s1's pooled Gata4 factor on the reference WT's real cells at its "
                                        "own composition", **K2_PARAMS}, requires=_k2_applies),
    "t3_factor_xko_perm": Candidate("t3_factor_xko_perm", frozenset({"T3"}), (), _t3_factor_xko_perm,
                                    {"what": "try K2b: K2a with its factor permuted across genes within deciles of WT "
                                             "E8.75 mean expression (log1p CP10k over the deduplicated cells)",
                                     **K2_PARAMS, "permutation_seed": K2_PERM_SEED, "bins": K2_PERM_BINS},
                                    requires=_k2_applies),
}


# ------------------------------------------------------------------------------------------------ read guard
_GUARD: dict | None = None  # the active guard's rules and log, if a candidate is building
_HOOKED = False


def _under(rel: Path, pattern: str) -> bool:
    """Whether a repository-relative path lies at or under a pattern such as 'experiments/*/results'."""
    parts = Path(pattern).parts
    return len(rel.parts) >= len(parts) and all(fnmatch(a, b) for a, b in zip(rel.parts, parts, strict=False))


def _check_read(path) -> None:
    g = _GUARD
    if g is None or not isinstance(path, (str, bytes, os.PathLike)):
        return
    try:
        p = Path(os.fsdecode(path)).resolve()
    except (TypeError, ValueError, OSError):
        return
    try:
        rel = p.relative_to(g["root"])
    except ValueError:
        if g["data"] in p.parents or p == g["data"]:
            rel = None
        else:
            return
    if rel is not None and (rel.parts[:1] == (".venv",) or "__pycache__" in rel.parts or
                            p.suffix in (".py", ".pyc", ".so")):
        return
    shown = str(rel) if rel is not None else str(p)
    g["log"].append(shown)
    denied = rel is None or any(_under(rel, pat) for pat in DENY_DIRS)
    if denied or str(rel) in DENY_FILES:
        g["refused"].append(shown)
        raise PermissionError(f"read refused while a candidate builds: {shown}")


def _audit(event, args) -> None:
    if event == "open" and args:
        _check_read(args[0])


_H5_INIT = h5py.File.__init__


def _guarded_h5_init(self, name, *a, **kw):
    _check_read(name)
    return _H5_INIT(self, name, *a, **kw)


@contextmanager
def guarded_reads(root: Path = REPO_ROOT, data_dir: Path = DATA_DIR):
    """While active, refuse (PermissionError) and log every read of data/raw, runs/, experiments/*/results/,
    docs/research/ or boards/t1_composition.json, and log every other file opened under the repository. Covers
    Python file opens (audit hook) and HDF5 opens (h5py.File, which anndata uses)."""
    global _GUARD, _HOOKED
    if not _HOOKED:
        sys.addaudithook(_audit)
        _HOOKED = True
    _GUARD = {"root": Path(root).resolve(), "data": Path(data_dir).resolve(), "log": [], "refused": []}
    h5py.File.__init__ = _guarded_h5_init
    try:
        yield _GUARD
    finally:
        h5py.File.__init__ = _H5_INIT
        _GUARD = None


# ------------------------------------------------------------------------------------------------ running
def _candidate_board(board: Board) -> Board:
    """The board as a candidate sees it: limits only (the test boards may publish only their limits)."""
    return replace(board, truth_cells=None, ref_cells=None, anchors={})


def score_seeds(t: PseudoTask, board: Board, target: Stage, ref: Stage, inputs: Mapping[str, ad.AnnData],
                candidates: list[Candidate], seeds=SEEDS, build: str = BUILD, log=None,
                root: Path = REPO_ROOT, data_dir: Path = DATA_DIR) -> list[dict]:
    """Per seed: the floor (checked to score 50 with complete anchors) and every applicable candidate, all against
    one seed's anchors; candidates build inside the read guard. Each row keeps the full metric panel per candidate
    (`metrics`, de_score's audit trail included) and its group skills (`groups`), so a result can be explained without
    rebuilding it; on T3 also what the response metrics are made of (`t3`: vec.evaluate.t3_response_detail, the
    severity gate and de_direction's per-gene terms, for the fit checks of tries K2a and K2b)."""
    view, cboard = TaskView.of(t), _candidate_board(board)
    given: dict[str, dict[str, object]] = {}
    for c in candidates:
        if c.name != FLOOR and c.applies(t):
            check_externals(c, t)
            given[c.name] = load_derived(c, t)
    rows = []
    for i, s in enumerate(seeds):
        t0 = time.time()
        floor_pred = subsample(ref, SUBSAMPLE_FRAC, s)
        f = score_board(board, floor_pred, target, ref, seed=s, anchors="computed", build=build, protocol=PROTOCOL)
        missing = sorted(set(RANKED_METRICS[t.task]) - set(f.anchors))
        if abs(f.score - 50.0) > 1e-9 or f.extras["inverted_anchors"] or missing:
            raise ValueError(f"{t.key} seed {s}: the floor scores {f.score} (inverted {f.extras['inverted_anchors']},"
                             f" missing anchors {missing})")
        row = {"seed": s, "build_seed": BUILD_SEED_OFFSET + s, "scores": {FLOOR: f.score}, "skills": {FLOOR: f.skills},
               "metrics": {FLOOR: dict(f.metrics)}, "n_cells": {FLOOR: int(round(ref.n * SUBSAMPLE_FRAC))},
               "groups": {FLOOR: dict(f.group_skill)}, "anchors": f.anchors, "reads": {}}
        detail = None
        if t.task == "T3":
            from common.core_metrics import de_genes

            truth, _, ref_w = working_copies(target, ref, s, SUBSAMPLE_FRAC, PROTOCOL)
            de = de_genes(truth.X, ref_w.X)

            def detail(pred, metrics, truth=truth, ref_w=ref_w, de=de):
                return t3_response_detail(pred.X, truth.X, ref_w.X, de, metrics)

            row["t3"] = {FLOOR: detail(floor_pred, f.metrics)}
        cross_checked = False
        for c in candidates:
            if c.name == FLOOR or not c.applies(t):
                continue
            with guarded_reads(root, data_dir) as g:
                pred = c.build(BuildContext(task=view, board=cboard, inputs=inputs, seed=BUILD_SEED_OFFSET + s,
                                            derived=given[c.name]))
            row["reads"][c.name] = sorted(set(g["log"]))
            r = score_board(board, pred, target, ref, seed=s, anchors=f.anchors, build=build, protocol=PROTOCOL)
            if i == 0 and not cross_checked:
                full = score_board(board, pred, target, ref, seed=s, anchors="computed", build=build,
                                   protocol=PROTOCOL)
                if full.score != r.score:
                    raise ValueError(f"{t.key} seed {s}: shared anchors give {r.score}, fresh ones {full.score}")
                cross_checked = True
            row["scores"][c.name], row["skills"][c.name], row["n_cells"][c.name] = r.score, r.skills, pred.n
            row["metrics"][c.name], row["groups"][c.name] = dict(r.metrics), dict(r.group_skill)
            if detail is not None:
                row["t3"][c.name] = detail(pred, r.metrics)
        row["seconds"] = round(time.time() - t0, 1)
        rows.append(row)
        if log:
            log(f"{t.key} seed {s}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["scores"].items()))
    return rows


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def filter_target(t: PseudoTask, data_dir: Path, genes) -> tuple[Stage, dict]:
    """The released T1 target cut by the task's target filter (fixed in the config before any result): the cells of
    the kept categories of configs/t1_chamber_classes.yaml, plus a share of contamination cells drawn at the filter's
    seed from the named categories; returns the scoring Stage and a record (counts only)."""
    from vec import t1_scope

    spec = t1_scope.load_spec()
    a = read_stage(data_dir / t.target, genes)
    t1_scope.classify(a, source=stem(t.target))
    cat = np.asarray(a.obs[t1_scope.CATEGORY_KEY]).astype(str)
    f = t.target_filter
    keep = set(spec.in_scope(f["keep"]))
    kept = np.flatnonzero(np.isin(cat, sorted(keep)))
    pool = np.flatnonzero(np.isin(cat, list(f.get("contamination_from", ()))))
    c = float(f.get("contamination", 0.0))
    k = min(len(pool), int(round(c / (1.0 - c) * len(kept)))) if c > 0 else 0
    extra = np.random.default_rng(int(f["seed"])).choice(pool, k, replace=False) if k else np.array([], int)
    rows = np.sort(np.concatenate([kept, extra]))
    stage = Stage.from_anndata(a[rows], genes)
    return stage, {"kept": int(len(kept)), "contamination": int(k), "of": int(a.n_obs)}


def run_task(t: PseudoTask, candidates: list[Candidate], seeds=SEEDS, data_dir: Path = DATA_DIR,
             boards: Mapping[str, Board] | None = None, build: str = BUILD, log=print) -> dict:
    """Load one pseudo-task's stages and score the candidates on it (DGX: the files live under data/raw)."""
    commit, clean = _git("rev-parse", "HEAD"), _git("status", "--porcelain", "--untracked-files=no") == ""
    boards = boards if boards is not None else load_boards()
    check_index(boards)
    check_task(t)
    genes = boards[t.panel_board].genes
    need_coords = boards[t.panel_board].needs_coords
    if t.target_filter is not None:
        target, filt = filter_target(t, data_dir, genes)
    else:
        target, filt = Stage.from_anndata(read_stage(data_dir / t.target, genes), genes, need_coords=need_coords), None
    ref = Stage.from_anndata(read_stage(data_dir / t.reference, genes), genes, need_coords=need_coords)
    board = task_board(t, boards, target.n, ref.n)
    if t.task == "T1":
        inputs = {r: load_stage(data_dir / f, lineage=False) for r, f in t.inputs.items()}
    else:
        inputs = {r: load_stage(data_dir / f, genes=genes, drop_duplicates=True) for r, f in t.inputs.items()}
    for r, a in inputs.items():
        missing = [g for g in genes if g not in a.var_names]
        if missing:
            raise ValueError(f"{t.key}: input {r} lacks {len(missing)} panel genes")
    applicable = [FLOOR] + [c.name for c in candidates if c.name != FLOOR and c.applies(t)]
    rows = score_seeds(t, board, target, ref, inputs, candidates, seeds=seeds, build=build, log=log,
                       data_dir=data_dir)
    files = sorted({t.target, t.reference, *t.inputs.values()})
    used = [c for c in candidates if c.name in applicable[1:]]
    derived = sorted({n for c in used for n in c.derived_for(t)})
    return {
        "task": asdict(t) | {"inputs": dict(t.inputs), "forbidden": list(t.forbidden),
                             "target_filter": dict(t.target_filter) if t.target_filter else None},
        "board": {k: getattr(board, k) for k in ("key", "task", "setting", "split", "target_ko", "needs_coords",
                                                   "min_cells", "max_cells", "truth_cells", "ref_cells",
                                                   "floor_model")},
        "target_filter_record": filt,
        "candidates": applicable,
        "params": {c.name: dict(c.params) | {"tuned_on": sorted(c.tuned_on), "externals": list(c.externals)}
                   | ({"derived": derived_record(c, t)} if c.derived_for(t) else {})
                   for c in candidates if c.name in applicable},
        "protocol": PROTOCOL, "build": build, "anchors": "computed", "frac": SUBSAMPLE_FRAC,
        "seeds": list(seeds), "rows": rows,
        "provenance": {"commit": commit, "tree_clean": clean, "config_sha256": _sha256(CONFIG),
                       "files_sha256": {f: _sha256(data_dir / f) for f in files},
                       "derived_sha256": {n: _sha256(DERIVED[n].path) for n in derived},
                       "build_seed_offset": BUILD_SEED_OFFSET,
                       "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())},
    }


# ------------------------------------------------------------------------------------------------ deciding
class InconsistentResults(ValueError):
    """Results that may not be combined into one decision."""


def check_results(results: list[dict], configured: list[PseudoTask], config_sha256: str | None = None) -> None:
    """Refuse results that mix commits, configs, protocols, builds, anchors, fractions, seeds or build seeds, come
    from a dirty tree, carry an unknown or repeated pseudo-task, or disagree on a shared file's or derived input's
    checksum."""
    keys = [r["task"]["key"] for r in results]
    known = {t.key for t in configured}
    problems = []
    if len(keys) != len(set(keys)):
        problems.append(f"repeated pseudo-tasks {sorted({k for k in keys if keys.count(k) > 1})}")
    if set(keys) - known:
        problems.append(f"pseudo-tasks not in the config {sorted(set(keys) - known)}")
    for field_, want in (("protocol", PROTOCOL), ("build", BUILD), ("anchors", "computed"), ("frac", SUBSAMPLE_FRAC)):
        got = {r[field_] for r in results}
        if got != {want}:
            problems.append(f"{field_} {sorted(map(str, got))}, not {want}")
    if {tuple(r["seeds"]) for r in results} != {tuple(SEEDS)}:
        problems.append("seeds other than vec.adoption.SEEDS")
    commits = {r["provenance"]["commit"] for r in results}
    if len(commits) != 1:
        problems.append(f"{len(commits)} commits")
    if not all(r["provenance"]["tree_clean"] for r in results):
        problems.append("a run from a dirty tree")
    shas = {r["provenance"]["config_sha256"] for r in results}
    if len(shas) != 1 or (config_sha256 is not None and shas != {config_sha256}):
        problems.append("a pseudo-task config other than the current one")
    for r in results:
        if any(row["build_seed"] != BUILD_SEED_OFFSET + row["seed"] for row in r["rows"]):
            problems.append(f"{r['task']['key']}: build seeds off the offset")
    seen: dict[str, str] = {}
    for r in results:
        prov = r["provenance"]
        for f, h in [*prov["files_sha256"].items(), *prov.get("derived_sha256", {}).items()]:
            if seen.setdefault(f, h) != h:
                problems.append(f"{f}: two checksums")
    if problems:
        raise InconsistentResults("; ".join(problems))


def decide(results: list[dict], configured: list[PseudoTask], base: str, change: str,
           candidates: Mapping[str, Candidate] = CANDIDATES, config_sha256: str | None = None) -> dict:
    """The adoption rule for `change` against `base` on one test board, with the readings the captain approved on
    2026-11-01 ("yes for all"):
      - counted are the board's configured pseudo-tasks that are not diagnostics or replicate folds
        (PseudoTask.counted), on which `change` applies (a pseudo-task it cannot run on is not counted; the record
        says how many of the configured pseudo-tasks decided), on which `base` applies too (the captain's reading of
        2026-11-01, a fourth "yes for all": a pseudo-task the base cannot run on is not counted either, e.g. the
        embryo's E7.25 task for pseudobulk_shift, which has no stage before E6.75), and which `change` was not tuned
        on (if none remain, the change is not kept);
      - on T1 the tissue mix is a premise (the organisers' E12.5 dissection note), so candidates of the same declared
        scope (params["scope"]: "whole" for the reference's own mix, "chamber" for the heart chamber only) are judged
        on every counted pseudo-task, and candidates of different scopes (decision point D2: B0 against the floor)
        only on the pseudo-tasks whose truth is cut to the chamber; the whole-E10.5 result is still reported.
    A counted pseudo-task with no result, or without both candidates on every seed, is incomplete, so the change is
    not kept."""
    check_results(results, configured, config_sha256)
    by_key = {r["task"]["key"]: r for r in results}
    c = candidates.get(change)
    b = candidates.get(base)
    cross_scope = (c is not None and b is not None and c.params.get("scope") and b.params.get("scope")
                   and c.params["scope"] != b.params["scope"])
    counted, excluded = [], {}
    for t in configured:
        if not t.counted:
            excluded[t.key] = "replicate fold" if t.replicate else "diagnostic"
        elif c is not None and not c.applies(t):
            excluded[t.key] = "change does not apply"
        elif b is not None and not b.applies(t):
            excluded[t.key] = "base does not apply"
        elif c is not None and t.key in c.tuned_on:
            excluded[t.key] = "change was tuned on it"
        elif cross_scope and t.task == "T1" and (t.target_filter or {}).get("keep") != "chamber":
            excluded[t.key] = "different tissue mixes are judged on the chamber-cut truth"
        else:
            counted.append(t.key)
    paired = {}
    for k in counted:
        r = by_key.get(k)
        if r is None or base not in r["candidates"] or change not in r["candidates"]:
            paired[k] = {}
            continue
        paired[k] = {row["seed"]: (row["scores"][base], row["scores"][change]) for row in r["rows"]
                     if base in row["scores"] and change in row["scores"]}
    out = keep_change(paired)
    return {"base": base, "change": change, "counted": counted, "excluded": excluded,
            "configured": [t.key for t in configured],
            "decided_on": f"{len(counted)} of {len(configured)} configured pseudo-tasks", **out}


def t3_fit_checks(result: dict, base: str, change: str, top: int = 10) -> dict | None:
    """What the column "A fit, not an effect, if" of tries K2a and K2b (experiments/harness/TRIES.md) reads, for
    `change` against `base` on one T3 pseudo-task result, over the seeds whose rows hold both (None if none does):
      group_margins: per scored group (vec.boards.GROUPS), the mean over seeds of 100 x weight x (change's group skill
        - base's); they sum to the mean score margin (score_margin);
      de_direction_margin: the mean over seeds of the difference of de_direction's per-gene terms
        (vec.evaluate.t3_response_detail), summed over genes; top_share: the share of it carried by the `top` genes
        with the largest seed-mean differences (None unless the margin is positive);
      severity_gate_differs: the seeds where the severity gate (scored as the worst value: an r2 or beta gate, or no
        response) holds for one of the two only; severity_states: both states per seed;
      matches_panel: whether every detail used reproduced its candidate's panel values."""
    rows = [r for r in result["rows"] if {base, change} <= set(r.get("t3", {})) & set(r.get("groups", {}))]
    if not rows:
        return None
    margins = {name: float(np.mean([100.0 * w * (r["groups"][change][name] - r["groups"][base][name]) for r in rows]))
               for name, w, _ in GROUPS["T3"]}
    d = np.mean([np.asarray(r["t3"][change]["de_direction_terms"]) - np.asarray(r["t3"][base]["de_direction_terms"])
                 for r in rows], axis=0)
    margin = float(d.sum())
    sev = {r["seed"]: {n: r["t3"][n]["severity"] for n in (base, change)} for r in rows}
    return {"base": base, "change": change, "seeds": sorted(sev), "group_margins": margins,
            "score_margin": float(sum(margins.values())), "de_direction_margin": margin, "top": top,
            "top_share": float(np.sort(d)[::-1][:top].sum() / margin) if margin > 0 else None,
            "severity_gate_differs": [s for s, v in sev.items() if v[base]["gated"] != v[change]["gated"]],
            "severity_states": {s: {n: v[n]["state"] for n in v} for s, v in sev.items()},
            "matches_panel": all(all(r["t3"][n].get("matches_panel", {"none": False}).values())
                                 for r in rows for n in (base, change))}
