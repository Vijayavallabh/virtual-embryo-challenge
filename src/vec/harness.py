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
(Candidate.externals, by stage day) and is refused inside the pseudo-task's version of the test's window.
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
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import yaml

from vec import t2_build
from vec.adoption import SEEDS, keep_change
from vec.boards import RANKED_METRICS, REPO_ROOT, VALIDATION_BOARDS, Board, load_boards
from vec.core import draw_balanced
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
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
    input stages by role (duplicates dropped on T2/T3; never the target) and the build seed. Candidates must not
    modify the inputs."""

    task: TaskView
    board: Board
    inputs: Mapping[str, ad.AnnData]
    seed: int


@dataclass(frozen=True)
class Candidate:
    """A named builder. It applies to a pseudo-task of one of `tasks` when `requires` (default: every role in one of
    `needs`, a tuple of alternative role sets) says so. `tuned_on` lists pseudo-task keys the candidate was tuned or
    selected on, which `decide` does not count for it; `externals` declares external data as (EXTERNAL_DATA id,
    first day, last day), refused inside a pseudo-task's version of the test window."""

    name: str
    tasks: frozenset[str]
    needs: tuple[frozenset[str], ...]
    build: Callable[[BuildContext], Stage] | None
    params: Mapping = field(default_factory=dict)
    requires: Callable[[PseudoTask], bool] | None = None
    tuned_on: frozenset[str] = frozenset()
    externals: tuple[tuple[str, float, float], ...] = ()

    def applies(self, t: PseudoTask) -> bool:
        if t.task not in self.tasks:
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


def base_role(t: PseudoTask) -> str:
    """The input role holding the reference stage's own cells."""
    return next(role for role, f in t.inputs.items() if f == t.reference)


def count_for(board: Board, pool_size: int) -> int:
    """The approved count rule (max_cells, within the board's limits) on any board."""
    return t2_build.n_cells(board, pool_size)[0]


def draw_copy(ctx: BuildContext) -> tuple[ad.AnnData, np.ndarray, np.ndarray, np.ndarray | None]:
    """Real reference cells at the board's count: (pool, rows, X float32 in board gene order, coords or None)."""
    pool = ctx.inputs[ctx.task.reference_role]
    n = count_for(ctx.board, pool.n_obs)
    if ctx.board.needs_coords:
        smp = t2_build.sample_cells(pool, n, seed=ctx.seed)
        rows, coords = smp.rows, smp.coords
    else:
        rows, coords = draw_balanced(pool.n_obs, n, np.random.default_rng([ctx.seed, 0])), None
    X = dense(pool[:, ctx.board.genes].X[rows]).astype(np.float32)
    return pool, rows, X, coords


def finish(ctx: BuildContext, pool: ad.AnnData, X: np.ndarray, coords: np.ndarray | None) -> Stage:
    """A candidate as the scorer reads its file, after the checks the server applies: written and checked with
    vec.io (T1), or through the T2 head's cloud guard and checked builder (coordinate boards). A failure raises."""
    X = np.clip(np.asarray(X, dtype=np.float32), 0.0, None)
    if not ctx.board.needs_coords:
        sub = make_submission(X, ctx.board.genes)
        with tempfile.TemporaryDirectory() as tmp:
            rep = check_submission(write_submission(sub, Path(tmp) / "submission.h5ad", stamp=False), ctx.board)
        if not rep.ok:
            raise ValueError(f"{ctx.task.key}: the built file fails its check: {rep.errors}")
        return Stage.from_anndata(sub, ctx.board.genes)
    guard = t2_build.cloud_guard(coords, np.asarray(pool.obsm[COORD_KEY], dtype=float)[:, :3], seed=ctx.seed)
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
    ref_m = _type_means(pool, ctx.board.genes)
    start_m = _type_means(ctx.inputs[start_role(ctx.task)], ctx.board.genes)
    ct = np.asarray(pool.obs["celltype"]).astype(str)[rows]
    for c in np.unique(ct):
        if c in start_m:
            X[ct == c] += (ref_m[c] - start_m[c]).astype(np.float32)
    return finish(ctx, pool, X, coords)


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


def _b0_t1(damping: float) -> Callable[[BuildContext], Stage]:
    def build(ctx: BuildContext) -> Stage:
        from vec import t1_scope

        role = ctx.task.reference_role
        ref, src = _shallow(ctx.inputs[role]), ctx.task.role_stems[role]
        shares = t1_chamber_shares(ref, src)
        start = _shallow(ctx.inputs["start"]) if damping else None
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, count_for(ctx.board, ref.n_obs), ctx.seed,
                                           start=start, damping=damping, ref_source=src,
                                           start_source=ctx.task.role_stems["start"] if damping else None)
        return finish(ctx, ref, dense(sub.X), None)
    return build


def _b0_t2(ctx: BuildContext) -> Stage:
    """B0 on a T2 board: real reference cells at the count rule, moved in count space by the per-type change of the
    harmonised lineage level (obs["merged"]; groups without a stable cross-stage meaning held flat). Extrapolation:
    from `start` to the reference, damped by B0_DAMPING. Interpolation: from the earlier bracket towards the later one,
    by the target's fraction of the bracket (its day is part of the geometry). Coordinates and size stay the source's
    (the heart size cannot be backtested; the test build applies the pre-registered size)."""
    from vec import t1_scope

    pool, rows, _, coords = draw_copy(ctx)
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
    return finish(ctx, pool, dense(shifted[:, ctx.board.genes].X), coords)


def _b0(ctx: BuildContext) -> Stage:
    return _b0_t1(B0_DAMPING)(ctx) if ctx.task.task == "T1" else _b0_t2(ctx)


def _b0_applies(t: PseudoTask) -> bool:
    if t.task == "T1" or (t.task == "T2" and t.kind == "extrap"):
        return "start" in t.inputs
    return t.task == "T2" and t.kind == "interp" and {"before", "after"} <= set(t.inputs)


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
                     "scope": "chamber"}, requires=_b0_applies),
    "b0_scope": Candidate("b0_scope", frozenset({"T1"}), _ANY, _b0_t1(0.0),
                          {"what": "T1 chamber-only copy without a shift (B0's composition alone)",
                           "scope": "chamber"}),
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
    one seed's anchors; candidates build inside the read guard."""
    view, cboard = TaskView.of(t), _candidate_board(board)
    for c in candidates:
        if c.name != FLOOR and c.applies(t):
            check_externals(c, t)
    rows = []
    for i, s in enumerate(seeds):
        t0 = time.time()
        f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, s), target, ref, seed=s, anchors="computed",
                        build=build, protocol=PROTOCOL)
        missing = sorted(set(RANKED_METRICS[t.task]) - set(f.anchors))
        if abs(f.score - 50.0) > 1e-9 or f.extras["inverted_anchors"] or missing:
            raise ValueError(f"{t.key} seed {s}: the floor scores {f.score} (inverted {f.extras['inverted_anchors']},"
                             f" missing anchors {missing})")
        row = {"seed": s, "build_seed": BUILD_SEED_OFFSET + s, "scores": {FLOOR: f.score}, "skills": {FLOOR: f.skills},
               "n_cells": {FLOOR: int(round(ref.n * SUBSAMPLE_FRAC))}, "anchors": f.anchors, "reads": {}}
        cross_checked = False
        for c in candidates:
            if c.name == FLOOR or not c.applies(t):
                continue
            with guarded_reads(root, data_dir) as g:
                pred = c.build(BuildContext(task=view, board=cboard, inputs=inputs, seed=BUILD_SEED_OFFSET + s))
            row["reads"][c.name] = sorted(set(g["log"]))
            r = score_board(board, pred, target, ref, seed=s, anchors=f.anchors, build=build, protocol=PROTOCOL)
            if i == 0 and not cross_checked:
                full = score_board(board, pred, target, ref, seed=s, anchors="computed", build=build,
                                   protocol=PROTOCOL)
                if full.score != r.score:
                    raise ValueError(f"{t.key} seed {s}: shared anchors give {r.score}, fresh ones {full.score}")
                cross_checked = True
            row["scores"][c.name], row["skills"][c.name], row["n_cells"][c.name] = r.score, r.skills, pred.n
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
    return {
        "task": asdict(t) | {"inputs": dict(t.inputs), "forbidden": list(t.forbidden),
                             "target_filter": dict(t.target_filter) if t.target_filter else None},
        "board": {k: getattr(board, k) for k in ("key", "task", "setting", "split", "target_ko", "needs_coords",
                                                   "min_cells", "max_cells", "truth_cells", "ref_cells",
                                                   "floor_model")},
        "target_filter_record": filt,
        "candidates": applicable,
        "params": {c.name: dict(c.params) | {"tuned_on": sorted(c.tuned_on), "externals": list(c.externals)}
                   for c in candidates if c.name in applicable},
        "protocol": PROTOCOL, "build": build, "anchors": "computed", "frac": SUBSAMPLE_FRAC,
        "seeds": list(seeds), "rows": rows,
        "provenance": {"commit": commit, "tree_clean": clean, "config_sha256": _sha256(CONFIG),
                       "files_sha256": {f: _sha256(data_dir / f) for f in files},
                       "build_seed_offset": BUILD_SEED_OFFSET,
                       "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())},
    }


# ------------------------------------------------------------------------------------------------ deciding
class InconsistentResults(ValueError):
    """Results that may not be combined into one decision."""


def check_results(results: list[dict], configured: list[PseudoTask], config_sha256: str | None = None) -> None:
    """Refuse results that mix commits, configs, protocols, builds, anchors, fractions, seeds or build seeds, come
    from a dirty tree, carry an unknown or repeated pseudo-task, or disagree on a shared file's checksum."""
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
        for f, h in r["provenance"]["files_sha256"].items():
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
        says how many of the configured pseudo-tasks decided), and which `change` was not tuned on (if none remain,
        the change is not kept);
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
