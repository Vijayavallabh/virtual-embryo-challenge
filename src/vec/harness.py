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
computed from (Candidate.derived; elsewhere the candidate does not apply, unless it takes the input only where it can,
Candidate.derived_optional). External statistics a candidate declares in Candidate.externals are handed over the same
way, reduced to the declared stages (the T1 late-proxy tries, DERIVED "t1_late_proxy"). A T1 pseudo-task's target
filter may name a chamber category (`within`), so that a read is made within one cell type: the target, the reference
and the inputs are cut to it alike (the late-proxy tournament, configs/t1_late_proxy.yaml).
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable, Mapping, Sequence
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


# T1 late-proxy tries T13-T18 (experiments/harness/TRIES.md, registered before their run; the captain's approval of
# 2026-11-10, docs/PLAN.md decision log, and the X002 amendment in EXTERNAL_DATA.md): the stand-in direction's Feng
# bracket and the only stages of Feng's statistics a candidate is ever handed.
LP_BRACKET = ("E14.5", "E18")  # the file's own label for its last embryonic stage (N20's first run, row N20b)
LP_BRACKET_DAYS = (14.5, 18.0)
LP_EXTERNAL = ("X002", 14.5, 18.0)  # as declared in Candidate.externals; legal on every T1 pseudo-task's window
LP_STATS = REPO_ROOT / "data" / "external" / "feng2022_gse193346_geo" / "cd1_stats.npz"  # feng_dirs.clean_stats_path()


def _feng_dirs():
    """experiments/t1_m1/feng_dirs.py (M1's direction code; `.dl` the readout's dirlib, `.fio` the filter's
    manifest check) and through it the late-proxy tries' statistics, imported as tests and diagnose.py import them."""
    m1 = str(REPO_ROOT / "experiments" / "t1_m1")
    if m1 not in sys.path:
        sys.path.insert(0, m1)
    import feng_dirs

    return feng_dirs


def _t1_hygiene():
    """experiments/t1_hygiene/hygiene.py (R2's modules; it reads only the committed config and extracts)."""
    d = str(REPO_ROOT / "experiments" / "t1_hygiene")
    if d not in sys.path:
        sys.path.insert(0, d)
    import hygiene

    return hygiene


def _late_proxy_sources() -> tuple[str, ...]:
    return ()  # an external source: no released T1 stage went into it (Candidate.externals declares its window)


@lru_cache
def _late_proxy_load(path: Path) -> dict:
    """The late-proxy tries' hand-over, read by the harness outside the read guard (data/ is refused to a candidate):
    Feng's clean statistics (X002, checked against the filter's manifest) reduced to the bracket's two end stages, so
    no other Feng stage reaches a candidate; the T1 direction readout's config; the held-flat list's named genes (its
    sex-linked, Y-linked and strain names, the part that needs no T1 stage); and R2's four module symbol lists (D and
    N2 from the committed extracts, A1 and N1 from their patterns on Feng's gene symbols)."""
    fd = _feng_dirs()
    fd.fio.verify_manifest_output(path)
    st = fd.dl.load_stats(path, fields=fd.dl.STAT_FIELDS)
    cfg = fd.load_config()
    named = sorted({str(g) for key in ("sex_linked", "y_linked", "strain_named") for g in cfg["held_flat"][key]})
    return {"stats": fd.bracket_stats(st, LP_BRACKET), "cfg": cfg, "named": named,
            "modules": _t1_hygiene().build_modules([str(g) for g in st["genes"]])}


# The late-proxy hand-over on Feng's recomputed groups (experiments/external/feng_relabel.py, register row X002; off
# unless a candidate built by `_lp(..., feng_groups="recomputed")` declares it; no such candidate is registered here: a
# backtest on it is a scored try whose row comes first). The relabelling labels E9.5 and E14.5 only, so this hand-over
# refuses (bracket_stats) until a row adds E18 to feng_relabel.yaml's stages.
LP_STATS_RECOMPUTED = (REPO_ROOT / "data" / "derived" / "feng_relabel_c"
                       / "cd1_stats_relabel.npz")  # feng_dirs.recomputed_stats_path() (feng_relabel.yaml output.dir)


@lru_cache
def _late_proxy_load_recomputed(path: Path) -> dict:
    """As `_late_proxy_load`, on the recomputed statistics: the file checked by feng_dirs.check_recomputed (manifest,
    written by feng_relabel from the yaml as it stands, usable, its manifest committed at HEAD as the labeller's
    RECORD) and check_recomputed_keys (no author label), and the
    readout's config with feng_relabel.yaml's group map in place of the authors' labels; the yaml rides along for its
    fallback."""
    fd = _feng_dirs()
    fd.check_recomputed(path)
    st = fd.dl.load_stats(path, fields=fd.dl.STAT_FIELDS)
    fd.check_recomputed_keys(st)
    rcfg = fd.load_relabel_config()
    cfg = fd.load_config()
    named = sorted({str(g) for key in ("sex_linked", "y_linked", "strain_named") for g in cfg["held_flat"][key]})
    return {"stats": fd.bracket_stats(st, LP_BRACKET), "cfg": fd.groups_config(cfg, rcfg), "named": named,
            "modules": _t1_hygiene().build_modules([str(g) for g in st["genes"]]), "relabel": rcfg}


DERIVED = {
    # the pre-registered list (experiments/t1_direction): genes named in its config (sex-linked, Y-linked, strain)
    # plus its pool rule on E8.5 and E9.5; committed 2026-10-19 (c992e7a), nine days before the E10.5 file reached us
    "t1_held_flat": Derived("T1", REPO_ROOT / "experiments" / "t1_direction" / "results" / "held_flat_genes.csv",
                            _held_flat_sources, _held_flat_load),
    # the late-proxy tries' statistics (above): external, so it enters through this hand-over and Candidate.externals
    "t1_late_proxy": Derived("T1", LP_STATS, _late_proxy_sources, _late_proxy_load),
    # the same on Feng's recomputed groups (above): read only by a candidate that declares it; none does yet
    "t1_late_proxy_recomputed": Derived("T1", LP_STATS_RECOMPUTED, _late_proxy_sources, _late_proxy_load_recomputed),
}
LP_DERIVED = {"authors": "t1_late_proxy", "recomputed": "t1_late_proxy_recomputed"}  # --feng-groups' two values


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
    if t.target_filter is not None and t.target_filter.get("within") is not None:
        from vec import t1_scope

        if t.target_filter["within"] not in t1_scope.load_spec().categories:
            raise ValueError(f"{t.key}: `within` must name a category of configs/t1_chamber_classes.yaml")
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
    it uses, and it does not apply where one of its task's cannot be handed over (`derived_ok`); `derived_optional`
    names inputs it takes only where they can be handed over and does without elsewhere (the late-proxy tries: the
    held-flat list's pool rule, computed on E8.5 and E9.5, cannot be handed over on the E8.5 -> E9.5 task, whose target
    is E9.5)."""

    name: str
    tasks: frozenset[str]
    needs: tuple[frozenset[str], ...]
    build: Callable[[BuildContext], Stage] | None
    params: Mapping = field(default_factory=dict)
    requires: Callable[[PseudoTask], bool] | None = None
    tuned_on: frozenset[str] = frozenset()
    externals: tuple[tuple[str, float, float], ...] = ()
    derived: tuple[str, ...] = ()
    derived_optional: tuple[str, ...] = ()

    def derived_for(self, t: PseudoTask) -> tuple[str, ...]:
        """The derived inputs this candidate uses on pseudo-task t: those of t's task, and its optional inputs where
        they can be handed over."""
        return tuple(n for n in self.derived if DERIVED[n].task == t.task) + tuple(
            n for n in self.derived_optional if DERIVED[n].task == t.task and derived_ok(n, t))

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


def _b0_t1(damping: float, herd: bool = False,
           prior: Callable[[BuildContext], object] | None = None) -> Callable[[BuildContext], Stage]:
    """B0 on T1 (or, with damping 0, its chamber copy alone). prior (T1 tries T12a and T12b): a function of the context
    that gives the strong-null down-tail prior for `vec.t1_scope.build_submission`."""
    def build(ctx: BuildContext) -> Stage:
        from vec import t1_scope

        role = ctx.task.reference_role
        ref, src = _shallow(ctx.inputs[role]), ctx.task.role_stems[role]
        shares = t1_chamber_shares(ref, src)
        start = _shallow(ctx.inputs["start"]) if damping else None
        held = list(ctx.derived["t1_held_flat"]) if damping else None  # handed over: the guard refuses the file
        extra = {} if prior is None else {"down_prior": prior(ctx)}
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, count_for(ctx.board, ref.n_obs), ctx.seed,
                                           start=start, damping=damping, held_flat=held, ref_source=src,
                                           start_source=ctx.task.role_stems["start"] if damping else None,
                                           herd=herd, **extra)
        return finish(ctx, ref, dense(sub.X), None)
    return build


def _copy_herd(ctx: BuildContext) -> Stage:
    """A copy of the reference with its draw herded towards the pool's mean expression (T1 try T6, a diagnostic of how
    much of a copy's de_score handicap is sampling noise; experiments/harness/TRIES.md)."""
    pool, rows, _, coords = draw_copy(ctx)
    g = ctx.board.genes
    rows, _ = herd_rows(pool[:, g].X, np.arange(pool.n_obs), rows, np.random.default_rng([ctx.seed, 9]))
    return finish(ctx, pool, dense(pool[:, g].X[rows]).astype(np.float32), coords)


# T1 try T11 (experiments/harness/TRIES.md, registered before its run): a mean-free within-type momentum on top of B0
_G_MOM_PREP: dict = {}


def _g_mom_factors(ctx: BuildContext, ref: ad.AnnData, src: str, start: ad.AnnData, start_src: str,
                   shares: Mapping[str, float]):
    """Try T11's per-metacell factors for this fold, fitted once per pair of loaded inputs (the harness hands every
    seed the same objects): OT within each chamber category from the start stage to the reference in a PCA of the fold's
    inputs, metacells of the reference cells, the log ratio of each metacell's mean CP10k to its transported partners'
    (start-stage) mean, the category's mean removed, damped by B0_DAMPING, zeroed and clipped
    (`vec.t1_state_transport`). ref and start: classified copies of the inputs, as `_g_mom` makes them."""
    from vec import t1_scope
    from vec import t1_state_transport as ts

    objs = (ctx.inputs[ctx.task.reference_role], ctx.inputs["start"])
    hit = _G_MOM_PREP.get("fold")
    if hit is None or any(a is not b for a, b in zip(hit["inputs"], objs, strict=True)):
        spec = t1_scope.load_spec()
        fit = ts.fit_state_factors(start, ref, "later", [c for c, s in shares.items() if s > 0],
                                   earlier_fallback=spec.fallbacks(start_src), later_fallback=spec.fallbacks(src))
        _G_MOM_PREP.clear()
        hit = _G_MOM_PREP["fold"] = {"inputs": objs, "factors": ts.finish_momentum(fit, B0_DAMPING)}
    return hit["factors"]


def _g_mom(ctx: BuildContext) -> Stage:
    """Try T11: B0's cells and shift (the same draws and factors, `vec.t1_state_transport.b0_population`) with a
    within-type momentum per reference-stage metacell on top, applied through `vec.core.apply_factors`; cells of R1's
    held types (labels with no Feng direction) get no within-type factor. Backtests: start E8.5, ref E9.5; at test
    start would be E9.5 and ref E10.5 (promotion needs the captain's ruling)."""
    from vec import t1_scope
    from vec import t1_state_transport as ts

    role = ctx.task.reference_role
    ref, src = _shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    start, start_src = _shallow(ctx.inputs["start"]), ctx.task.role_stems["start"]
    shares = t1_chamber_shares(ref, src)
    held = list(ctx.derived["t1_held_flat"])
    pop = ts.b0_population(ref, start, shares, count_for(ctx.board, ref.n_obs), ctx.seed, held, ref_source=src,
                           start_source=start_src, damping=B0_DAMPING)
    factors = _g_mom_factors(ctx, ref, src, start, start_src, shares)
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    no_edit = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    moved, _ = ts.apply_state_factors(pop, factors, held, [ctx.seed, 5], rescale=False, no_edit=no_edit)
    return finish(ctx, ref, dense(moved[:, ctx.board.genes].X).astype(np.float32), None)


# T1 late-proxy tries T13-T18 (experiments/harness/TRIES.md): their constants, fixed in their rows before any run.
LP_GROUPS = ("CM", "Endocardium", "Endothelium", "Epicardium")  # M1's four direction groups (R2's guard: never Blood)
LP_PERM_SEED = 20261110  # the permuted control: one permutation, this seed
LP_BINS = 10  # the permuted control's and the shrinkage's expression deciles (K2b's bins)


def _late_proxy(hold: tuple[str, ...] = (), flip: bool = False, permute: bool = False,
                shrink: bool = False, feng_groups: str = "authors") -> Callable[[BuildContext], Stage]:
    """T1 tries T13-T18: B0's chamber composition (b0_scope's cells, at the same seed) with the cells of M1's four
    direction groups moved along Feng's E14.5 -> E18 per-group log ratio of mean CP10k (genes mapped by symbol, as
    feng_dirs.feng_directions does) by the step's share of that 3.5-day bracket (a one-day step takes 2/7), in count
    space (vec.t1_scope.external_shift, as M1). No momentum and no other change: the cells without a direction stay as
    drawn, as R1 keeps M1's. Held flat: the pre-registered list where the pseudo-task's inputs include its pool rule's
    stages (derived_optional; on E8.5 -> E9.5, whose target is E9.5, only its named genes), plus R2's `hold` modules.
    flip negates the direction (control Cflip); permute shuffles each group's values within deciles of the group's mean
    CP10k averaged over the two ends (control Cperm, K2b's procedure, one permutation at LP_PERM_SEED); shrink is
    feng_dirs.eb_shrink's reliability weighting. feng_groups (default "authors", unchanged): "recomputed" reads the
    hand-over on Feng's recomputed groups (LP_DERIVED) and applies feng_relabel.yaml's fallback; with a fallback fired
    it refuses permute and shrink, whose halves exist only for a group's own cells."""
    if feng_groups not in LP_DERIVED:
        raise ValueError(f"feng_groups must be one of {sorted(LP_DERIVED)}, got {feng_groups!r}")

    def build(ctx: BuildContext) -> Stage:
        from vec import t1_scope

        fd = _feng_dirs()
        bundle = ctx.derived[LP_DERIVED[feng_groups]]
        role = ctx.task.reference_role
        ref, src = _shallow(ctx.inputs[role]), ctx.task.role_stems[role]
        genes = list(ref.var_names)
        own = fd.directions_from_stats(bundle["stats"], genes, bundle["cfg"], LP_BRACKET)[0]
        if feng_groups != "authors":
            own, fired = fd.apply_fallback(own, bundle["relabel"])
            if (permute or shrink) and any(v["fired"] for v in fired.values()):
                raise ValueError(f"{ctx.task.key}: a recomputed group fell back, and its halves cannot be read")
        dirs = {g: d for g, d in own.items() if g in LP_GROUPS}
        if not dirs:
            raise ValueError(f"{ctx.task.key}: Feng's statistics give no direction for any of {list(LP_GROUPS)}")
        if permute or shrink:
            halves = fd.direction_halves(bundle["stats"], genes, bundle["cfg"], LP_BRACKET, groups=set(dirs))
        if shrink:
            dirs = {g: fd.eb_shrink(d, halves[g]["d_a"], halves[g]["d_b"], halves[g]["mean_cp"], LP_BINS)[0]
                    for g, d in dirs.items()}
        if permute:
            for g in dirs:
                d, ok = dirs[g].copy(), np.isfinite(dirs[g])
                d[ok] = t3_decile_permutation(d[ok], halves[g]["mean_cp"][ok], [LP_PERM_SEED, LP_GROUPS.index(g)],
                                              LP_BINS)
                dirs[g] = d
        if flip:
            dirs = {g: -d for g, d in dirs.items()}
        frac = fd.bracket_fraction(ctx.task.role_days[role], ctx.task.target_day, LP_BRACKET_DAYS)
        held = {*bundle["named"], *map(str, ctx.derived.get("t1_held_flat", ()))}  # the pool rule only where legal
        for module in hold:
            held |= set(bundle["modules"][module])
        groups = fd.label_groups(ref.obs["celltype"].astype(str).unique(), bundle["cfg"])
        shares = t1_chamber_shares(ref, src)
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, count_for(ctx.board, ref.n_obs), ctx.seed,
                                           held_flat=sorted(held), ref_source=src, directions=dirs,
                                           direction_fraction=frac, direction_groups=groups)
        return finish(ctx, ref, dense(sub.X), None)

    return build


def g_mom_pair_factors(ref: ad.AnnData, src: str, start: ad.AnnData, start_src: str,
                       categories: Sequence[str] | None = None) -> tuple[object, object]:
    """Try T11's factors for one pair of released stages outside a pseudo-task (the test pair of `experiments/builds/
    build_test.py t1 --g-mom`, diagnostic N21 of experiments/harness/TRIES.md): the steps of `_g_mom_factors` (OT within
    each chamber category from `start` to `ref` in a joint PCA, metacells of the `ref` cells, the log ratio of each
    metacell's mean CP10k to its transported partners', the per-category mean removed, damped by B0_DAMPING, zeroed and
    clipped), on classified shallow copies of the two stages, the start stage's genes aligned to the reference's order
    by name when the two files order them differently (as `vec.core.counts_of` reads a start's genes; a gene the start
    lacks raises KeyError). categories default to the chamber categories with a positive share at `ref` (B0's rule);
    src, start_src are the file stems. Returns (finished factors, the raw fit); `ref` and `start` are not modified."""
    from vec import t1_scope
    from vec import t1_state_transport as ts

    spec = t1_scope.load_spec()
    ref_c, start_c = _shallow(ref), _shallow(start)
    if list(start_c.var_names) != list(ref_c.var_names):
        start_c = start_c[:, list(ref_c.var_names)].copy()
    t1_scope.classify(ref_c, source=src, spec=spec)
    t1_scope.classify(start_c, source=start_src, spec=spec)
    cats = ([c for c, s in t1_chamber_shares(ref, src).items() if s > 0] if categories is None else list(categories))
    raw = ts.fit_state_factors(start_c, ref_c, "later", cats, earlier_fallback=spec.fallbacks(start_src),
                               later_fallback=spec.fallbacks(src))
    return ts.finish_momentum(raw, B0_DAMPING), raw


# T1 tries T20, T21, T23 and T24 (experiments/harness/TRIES.md, registered before their run): g_mom's control arms,
# built like g_mom on the same draws and base (the same stream), so a paired comparison reads the factors alone
G_MOM_PERM_SEED = 20261110  # T21: one permutation, this seed (its row)
G_MOM_NULL_SEED = 20261110  # T23: one split of the reference stage into two halves, this seed (its row)
_G_MOM_NULL_PREP: dict = {}


def _g_mom_arm(ctx: BuildContext, transform: Callable[[object], object], *, rescale: bool = False) -> Stage:
    """g_mom's own build (`_g_mom`, whose steps this repeats so that `_g_mom` stays as registered; tests/
    test_g_mom_verify.py pins the two together) with `transform` applied to the fold's finished factors before they are
    applied: the same reference, shares, B0 population, fit and random streams, so the arms differ from g_mom only in
    the factors. rescale: as `vec.t1_state_transport.apply_state_factors` (False is g_mom's; True, try T24, holds each
    category's pseudobulk per gene, on the same stream)."""
    from vec import t1_scope
    from vec import t1_state_transport as ts

    role = ctx.task.reference_role
    ref, src = _shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    start, start_src = _shallow(ctx.inputs["start"]), ctx.task.role_stems["start"]
    shares = t1_chamber_shares(ref, src)
    held = list(ctx.derived["t1_held_flat"])
    pop = ts.b0_population(ref, start, shares, count_for(ctx.board, ref.n_obs), ctx.seed, held, ref_source=src,
                           start_source=start_src, damping=B0_DAMPING)
    factors = transform(_g_mom_factors(ctx, ref, src, start, start_src, shares))
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    no_edit = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    moved, _ = ts.apply_state_factors(pop, factors, held, [ctx.seed, 5], rescale=rescale, no_edit=no_edit)
    return finish(ctx, ref, dense(moved[:, ctx.board.genes].X).astype(np.float32), None)


def g_mom_no_development_factors(ref: ad.AnnData, src: str, categories: Sequence[str], target, held_flat: Sequence[str],
                                 seed: int = G_MOM_NULL_SEED):
    """Try T23's factors, the no-development control of g_mom: `ref` (the reference stage, the file stem `src`) split at
    random within each category into two halves (`ts.category_halves`, default_rng(seed)), T11's fit
    (`g_mom_pair_factors`) from one half to the other and back with each half in the start role in turn, so every
    reference cell of the fitted categories has a factor from the direction in which its own half carries the
    metacells, and the raw factors finished per category at the size of `target` (g_mom's finished factors,
    `ts.finish_no_development`, over the genes not in held_flat). categories: g_mom's fitted categories, so that the
    control edits the same cells as g_mom. Nothing develops between the halves, so whatever such factors gain is
    state-aligned spread and noise."""
    from vec import t1_scope
    from vec import t1_state_transport as ts

    stage = _shallow(ref)
    t1_scope.classify(stage, source=src, spec=t1_scope.load_spec())
    in_a = ts.category_halves(stage, np.random.default_rng(seed))
    half_a, half_b = ref[in_a].copy(), ref[~in_a].copy()
    _, raw_on_b = g_mom_pair_factors(half_b, src, half_a, src, categories=categories)
    _, raw_on_a = g_mom_pair_factors(half_a, src, half_b, src, categories=categories)
    return ts.finish_no_development(raw_on_b, raw_on_a, target, held_flat, B0_DAMPING)


def _g_mom_null_factors(ctx: BuildContext, finished) -> object:
    """The no-development factors for one fold, fitted once per loaded reference and g_mom fit (the harness hands every
    seed the same objects), as `_g_mom_factors` is."""
    role = ctx.task.reference_role
    ref, src = ctx.inputs[role], ctx.task.role_stems[role]
    hit = _G_MOM_NULL_PREP.get("fold")
    if hit is None or hit["ref"] is not ref or hit["target"] is not finished:
        _G_MOM_NULL_PREP.clear()
        factors = g_mom_no_development_factors(ref, src, sorted(finished.categories), finished,
                                               list(ctx.derived["t1_held_flat"]))
        hit = _G_MOM_NULL_PREP["fold"] = {"ref": ref, "target": finished, "factors": factors}
    return hit["factors"]


def _g_mom_null(ctx: BuildContext) -> Stage:
    """Try T23: g_mom's draws, base and stream with the no-development factors (`g_mom_no_development_factors`: the same
    pipeline between two halves of the reference stage, at g_mom's size per category) in place of g_mom's."""
    return _g_mom_arm(ctx, lambda finished: _g_mom_null_factors(ctx, finished))


def _g_mom_held(ctx: BuildContext) -> Stage:
    """Try T24: g_mom with each category's pseudobulk held per gene (the rescale of T9, `apply_state_factors` with
    rescale=True): the same factors, draws and stream, so what is left of g_mom's gain is within-type shape."""
    return _g_mom_arm(ctx, lambda finished: finished, rescale=True)


def _g_mom_flip(ctx: BuildContext) -> Stage:
    """Try T20: g_mom with every within-type factor's sign reversed (`vec.t1_state_transport.flip_factors`)."""
    from vec import t1_state_transport as ts

    return _g_mom_arm(ctx, ts.flip_factors)


def _g_mom_perm(ctx: BuildContext) -> Stage:
    """Try T21: g_mom with each category's metacell factor rows re-assigned among its metacells by one fixed derangement
    (`vec.t1_state_transport.permute_factors`, seed G_MOM_PERM_SEED)."""
    from vec import t1_state_transport as ts

    return _g_mom_arm(ctx, lambda f: ts.permute_factors(f, G_MOM_PERM_SEED))


# Heart candidate b0_trend (row H11 of experiments/harness/TRIES.md; ported from campaign/trend 8bb23bf after the
# captain's ruling of 2026-11-13, docs/PLAN.md decision log (w)): B0 with each lineage's per-gene velocity taken as the
# least-squares trend over every input stage instead of the last step alone.
B0_TREND_RULE = ("B0 with each lineage's per-gene velocity taken as the least-squares slope of log(mean CP10k + "
                 "eps) on the stage day over every input stage where the lineage has at least min_cells cells, "
                 "times B0_DAMPING x the last gap in days; a lineage seen at only the last two stages keeps B0's "
                 "own factor; draw, stream, depth rule and coordinates are B0's")
T2_LINEAGE_KEY = "merged"
T2_UNTRACKED = "Untracked"


def trend_stages(task: TaskView) -> list[tuple[str, float]]:
    """b0_trend's input stages, (role, day) in day order, from `task.roles` and `task.role_days` only (never the target
    or its day). The reference must be the last and `start` the one before it, with at least three stages and
    increasing days; otherwise b0_trend is not defined and this raises."""
    roles = sorted(task.roles, key=lambda r: (float(task.role_days[r]), r))
    if len(roles) < 3:
        raise ValueError("b0_trend needs at least three input stages (two give one transition: B0's own step)")
    if roles[-1] != task.reference_role or roles[-2] != "start":
        raise ValueError(f"the reference must be the latest input stage and `start` the one before it, got {roles}")
    days = [float(task.role_days[r]) for r in roles]
    if not all(b > a for a, b in zip(days, days[1:], strict=False)):
        raise ValueError("the input stages' days must increase")
    return list(zip(roles, days, strict=True))


def _b0_trend_task_check(t: TaskView) -> None:
    """b0_trend is defined on heart extrapolation tasks with at least three input stages, `start` second to last and
    the reference last (`trend_stages` raises otherwise)."""
    if t.task != "T2" or t.setting != "heart" or t.kind != "extrap":
        raise ValueError(f"b0_trend is defined on heart extrapolation tasks only, not {t.task} {t.setting} {t.kind}")
    trend_stages(t)


def _momentum_default(name: str):
    """`vec.t1_scope.momentum_shift`'s own default for one keyword (eps, min_cells), the values B0 runs with."""
    import inspect

    from vec import t1_scope

    return inspect.signature(t1_scope.momentum_shift).parameters[name].default


def trend_transition_weights(days: Sequence[float]) -> np.ndarray:
    """The weights with which the least-squares slope over stages at `days` uses each consecutive transition's
    velocity (for the record only): with c_i = (d_i - dbar) / S_dd, slope = sum_k a_k (y_k - y_{k-1}) with
    a_k = -sum_{i<k} c_i, and w_k = a_k * gap_k, so slope = sum_k w_k (y_k - y_{k-1}) / gap_k and the w_k sum to 1."""
    d = np.asarray(days, dtype=np.float64)
    c = (d - d.mean()) / float(((d - d.mean()) ** 2).sum())
    return -np.cumsum(c)[:-1] * np.diff(d)


def trend_log_factor(log_means: Sequence[np.ndarray], days: Sequence[float], budget: float) -> np.ndarray:
    """Per gene, `budget` x the ordinary least-squares slope of y_k (one array per stage) on day_k:
    slope = sum_k (d_k - dbar) y_k / sum_k (d_k - dbar)^2. b0_trend's budget is B0_DAMPING x the last gap in days."""
    d = np.asarray(days, dtype=np.float64)
    dc = d - d.mean()
    Y = np.vstack([np.asarray(y, dtype=np.float64) for y in log_means])
    return budget * (dc @ Y) / float((dc ** 2).sum())


def _cosine(x: np.ndarray, y: np.ndarray) -> float:
    nx, ny = float(np.linalg.norm(x)), float(np.linalg.norm(y))
    return float(x @ y / (nx * ny)) if nx > 0 and ny > 0 else float("nan")


def b0_trend_factors(ctx: BuildContext, genes: Sequence[str], types: Sequence[str]) -> tuple[dict[str, np.ndarray],
                                                                                              dict]:
    """b0_trend's log factors (lineage -> damped log factor in `genes` order) and its target-free record. For each
    lineage among `types` (the drawn cells' B0 types; "Untracked" never) that B0 shifts on this task (at least
    min_cells own cells at `start` and at the reference), the usable stages are the input stages (`trend_stages`: day
    order, never the target) where it has at least min_cells own cells (obs["merged"]); with y_k = log(mean CP10k +
    eps) there, its factor is B0_DAMPING x (day_ref - day_start) x the least-squares slope of y_k on day_k. A lineage
    usable at `start` and the reference only is left out of the mapping, so momentum_shift gives it B0's own factor
    from B0's own expression. The record holds, per shifted lineage, the usable stages and their cell numbers, the
    transition weights, and the cosine, L2-norm ratio (trend / B0) and mean absolute log factor against B0's own factor
    on the same rows; and the medians of cosine and norm ratio over the shifted lineages. eps and min_cells are
    momentum_shift's defaults, as B0 uses them. (A caller that shows the record for a test board drops the cell
    numbers: experiments/builds/build_test.py.)"""
    from vec.core import counts_of, mean_cp10k

    t = ctx.task
    _b0_trend_task_check(t)
    stages = trend_stages(t)
    ref_role = t.reference_role
    eps, min_cells = float(_momentum_default("eps")), int(_momentum_default("min_cells"))
    budget = B0_DAMPING * (float(t.role_days[ref_role]) - float(t.role_days["start"]))
    C = {r: counts_of(ctx.inputs[r], genes) for r, _ in stages}
    merged = {r: np.asarray(ctx.inputs[r].obs[T2_LINEAGE_KEY]).astype(str) for r, _ in stages}
    out: dict[str, np.ndarray] = {}
    lineages = []
    for lin in types:
        if lin == T2_UNTRACKED or lin.startswith(T2_UNTRACKED + "@"):
            continue
        rows = {r: np.flatnonzero(merged[r] == lin) for r, _ in stages}
        if len(rows["start"]) < min_cells or len(rows[ref_role]) < min_cells:
            continue  # B0 does not shift it
        usable = [(r, d) for r, d in stages if len(rows[r]) >= min_cells]
        y = {r: np.log(mean_cp10k(C[r][rows[r]]) + eps) for r, _ in usable}
        own = B0_DAMPING * (y[ref_role] - y["start"])
        if [r for r, _ in usable] == ["start", ref_role]:
            logf, kept = own, True
        else:
            logf = trend_log_factor([y[r] for r, _ in usable], [d for _, d in usable], budget)
            out[lin], kept = logf, False
        n_own = float(np.linalg.norm(own))
        lineages.append({"lineage": lin, "stages": [r for r, _ in usable], "days": [d for _, d in usable],
                         "n_cells": [int(len(rows[r])) for r, _ in usable], "own_factor": kept,
                         "weights": [float(w) for w in trend_transition_weights([d for _, d in usable])],
                         "cosine": _cosine(logf, own),
                         "norm_ratio": float(np.linalg.norm(logf)) / n_own if n_own > 0 else float("nan"),
                         "mean_abs_logf_trend": float(np.abs(logf).mean()),
                         "mean_abs_logf_b0": float(np.abs(own).mean())})
    cos = [r["cosine"] for r in lineages]
    ratio = [r["norm_ratio"] for r in lineages]
    rec = {"rule": "vec.harness.B0_TREND_RULE", "stages": [r for r, _ in stages], "days": [d for _, d in stages],
           "budget_days": budget, "eps": eps, "min_cells": min_cells, "n_trend": len(out), "lineages": lineages,
           "median_cosine": float(np.nanmedian(cos)) if lineages else float("nan"),
           "median_norm_ratio": float(np.nanmedian(ratio)) if lineages else float("nan")}
    return out, rec


def _b0_t2(ctx: BuildContext, depth: bool = False, herd: bool = False, cfm: bool = False, trend: bool = False,
           trend_record: dict | None = None) -> Stage:
    """B0 on a T2 board: real reference cells at the count rule, moved in count space by the per-type change of the
    harmonised lineage level (obs["merged"]; groups without a stable cross-stage meaning held flat). Extrapolation:
    from `start` to the reference, damped by B0_DAMPING. Interpolation: from the earlier bracket towards the later one,
    by the target's fraction of the bracket (its day is part of the geometry). Coordinates and size stay the source's
    (the heart size cannot be backtested; the test build applies the pre-registered size). depth (extrapolation only;
    heart try H3 of experiments/harness/TRIES.md): the moved cells are then binomially thinned to the target depth that
    vec.t2_build.DEPTH_EXTRAP_RULE extrapolates at B0_DAMPING from the start's and the reference's median depths.
    herd (try H4): the copy's draw herded towards the pool's mean panel expression (`draw_copy`). cfm (extrapolation
    only; heart try H6): each lineage's momentum goes through vec.t2_build.covariance_filter before it is applied.
    trend (heart candidate b0_trend, row H11; heart extrapolation with at least three input stages only): each
    lineage's log factor is its least-squares trend over every input stage (`b0_trend_factors`, B0_TREND_RULE); the
    draw, typing, stream and everything else are B0's. No try read it with depth, herd or cfm, so those raise.
    trend_record: a dict that receives b0_trend_factors' record (with trend only; the build record's source)."""
    from vec import t1_scope

    if trend and (depth or herd or cfm):
        raise ValueError("the trend of B0 is not defined together with depth, herd or cfm (no try read them)")
    if trend_record is not None and not trend:
        raise ValueError("trend_record is for trend=True")
    if trend:
        _b0_trend_task_check(ctx.task)
    if cfm and ctx.task.kind != "extrap":
        raise ValueError("the covariance filter of B0 is for extrapolation boards")
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
    logf_by_type = None
    if trend:
        logf_by_type, rec = b0_trend_factors(ctx, list(like.var_names),
                                             sorted(set(np.asarray(like.obs["b0_type"]).astype(str))))
        if trend_record is not None:
            trend_record.update(rec)
    shifted, _ = t1_scope.momentum_shift(typed(start, "start"), typed(end, "end"), like, damping, (), key="b0_type",
                                         seed=[ctx.seed, 3], logf_filter=t2_build.covariance_filter if cfm else None,
                                         logf_by_type=logf_by_type)
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
         move_coupling: str = t2_build.MOVE_COUPLINGS[0], local_perm: bool = False,
         draw: str = t2_build.DRAWS[0], t_mix: float | None = None,
         positions: str = t2_build.POSITIONS[0], shape_form: str | None = None) -> Stage:
    """The T2 embryo pipeline (docs/PLAN.md 4.2): displacement interpolation between the brackets at the target's
    time fraction (vec.t2_build.interpolate_clouds: fused Gromov-Wasserstein coupling, proper rotation, McCann moves,
    a mixture of real cells of both brackets), scaled to the embryo size rule read on the task's own stages. t_pos:
    the s2 shape variants' position fraction (0: the earlier bracket's geometry, 1: the later's; default: the
    mixture's). expression: "mixture" (the cells as released), "depth" (vec.t2_build.depth_matched_expression) or
    "blend" (vec.t2_build.blend_expression), tries E7 and E8 of experiments/harness/TRIES.md, or "latent" (try E14:
    vec.t2_latent, the task's negative-binomial VAE trained once on both brackets' panel counts). move_coupling: "fgw"
    (the entry) or "aligned" (try E10: the moves from vec.core.aligned_coupling in the rotated common frame).
    local_perm (try E15): after the draw, the positions are re-arranged among the drawn cells within coarse type and a
    fixed radius (vec.t2_build.local_type_permutation); the expression rows stay where they are. draw: "uniform" (the
    entry) or "midpoint" (try E13, vec.t2_build._midpoint_draw; it needs the inputs' counts layer). positions: "moved"
    (the entry) or "floor_exact" (tries E47a and E47b, with t_pos = 0: every cell on a real position of the earlier
    bracket, vec.t2_build._floor_exact_positions). shape_form (row E49, as drafted; with floor_exact only): a form of
    vec.shape_interp on the exact-position cloud before the size rule; left at None, nothing changes."""
    t = ctx.task
    before, after = ctx.inputs["before"], ctx.inputs["after"]
    lo, hi = t.role_days["before"], t.role_days["after"]
    g = ctx.board.genes
    frac = (t.target_day - lo) / (hi - lo)
    Eb, Ea = dense(before[:, g].X).astype(np.float32), dense(after[:, g].X).astype(np.float32)
    log_depth = (None if draw == t2_build.DRAWS[0]
                 else (t2_build.log_panel_counts(before, g), t2_build.log_panel_counts(after, g)))
    res = t2_build.interpolate_clouds(np.asarray(before.obsm[COORD_KEY], dtype=float)[:, :3],
                                      np.asarray(after.obsm[COORD_KEY], dtype=float)[:, :3], Eb, Ea,
                                      frac, count_for(ctx.board, before.n_obs), ctx.seed,
                                      target_rms=target_rms(ctx), t_pos=t_pos,
                                      with_partners=expression in ("blend", "latent"),
                                      herd=herd, move_coupling=move_coupling, draw=draw, log_depth=log_depth,
                                      with_coupling=expression == "latent", t_mix=t_mix, positions=positions,
                                      **({} if shape_form is None else {"shape_form": shape_form}))
    if expression == "mixture":
        X = np.empty((len(res.rows), len(g)), dtype=np.float32)
        X[res.stage == 0], X[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    elif expression == "depth":
        X, _ = t2_build.depth_matched_expression(before, after, res, g, frac, [ctx.seed, 7])
    elif expression == "blend":
        X, _ = t2_build.blend_expression(before, after, res, g, frac, [ctx.seed, 7])
    elif expression == "series":  # try E24: the shallower bracket's cells in the deeper bracket's run profile
        X, _ = t2_build.series_rendered_expression(before, after, res, g, frac, [ctx.seed, 7])
    elif expression == "series_wide":  # try E27: the same with the wider neighbourhood of diagnostic N24b
        X, _ = t2_build.series_rendered_expression(before, after, res, g, frac, [ctx.seed, 7],
                                                   k=t2_build.RENDER_K_WIDE)
    elif expression in RENDER_EXPRESSIONS:  # tries E19a-E19c: every cell moved to one depth (t2_build.RENDER_RULE)
        X, _ = t2_build.depth_rendered_expression(before, after, res, g, frac, [ctx.seed, 7],
                                                  mode=RENDER_EXPRESSIONS[expression])
    elif expression == "latent":  # try E14: t2_latent loads torch and scvi-tools only when it trains
        from vec import t2_latent

        X, _ = t2_latent.latent_midpoint_counts(t2_latent.task_model(before, after, g), before, after, res, g, frac,
                                                [ctx.seed, t2_latent.EXPRESSION_STREAM])
    else:
        raise ValueError(f"unknown expression step {expression!r}")
    coords = res.coords
    if local_perm:
        coords = coords[nbhd_rearrangement(coords, t2_build.drawn_labels(res, before.obs["group"], after.obs["group"]),
                                           ctx.seed)]
    return finish(ctx, before, X, coords, guard_source=res.pool_coords, guard_weights=res.pool_weights)


RENDER_EXPRESSIONS = {f"render_{m}": m for m in t2_build.RENDER_MODES}


def _fgw_series(ctx: BuildContext) -> Stage:
    """Embryo try E24: fgw's draws and positions; the shallower bracket's cells moved by the extraembryonic per-gene
    factor and to the deeper bracket's depth (vec.t2_build.SERIES_RULE)."""
    return _fgw(ctx, expression="series")


def _series_applies(t: PseudoTask) -> bool:
    return _interp(t)


def _fgw_shape0_hi(ctx: BuildContext) -> Stage:
    """Embryo try E24b: today's s2 geometry (the earlier bracket's, t_pos = 0) with the depth top-up of E19c."""
    return _fgw(ctx, t_pos=0.0, expression="render_hi")


def _fgw_shape0_series(ctx: BuildContext) -> Stage:
    """Embryo try E24c: today's s2 geometry with the series rendering of E24a."""
    return _fgw(ctx, t_pos=0.0, expression="series")


def _fgw_render(mode: str) -> Callable[[BuildContext], Stage]:
    """Embryo tries E19a (mid), E19b (lo) and E19c (hi): fgw's draws and positions, every cell moved to one depth."""
    return lambda ctx: _fgw(ctx, expression=f"render_{mode}")


def _fgw_shape0x(ctx: BuildContext) -> Stage:
    """Embryo try E47a: fgw_shape0's rows, mixture and size rule, every cell on a real position of the earlier bracket
    (its own rows where they are, the later bracket's rows on unused positions by an exact expression assignment)."""
    return _fgw(ctx, t_pos=0.0, positions=t2_build.POSITIONS[1])


# Row E49 (as drafted; branch campaign/embryo-shape-interp): fgw_shape0x with one shape form (vec.shape_interp) on its
# exact-position cloud before the size rule. Three arms and four controls; the builder builds the arms only.
SHAPE_CANDIDATES = {"fgw_shape0x_sr": "stretch_ratio", "fgw_shape0x_sg": "stretch_registered",
                    "fgw_shape0x_ov": "overlay", "fgw_shape0x_sr_anti": "stretch_ratio_anti",
                    "fgw_shape0x_sg_anti": "stretch_registered_anti", "fgw_shape0x_sg_rot": "stretch_registered_rot",
                    "fgw_shape0x_ov_rot": "overlay_rot"}


def _fgw_shape0x_form(form: str) -> Callable[[BuildContext], Stage]:
    return lambda ctx: _fgw(ctx, t_pos=0.0, positions=t2_build.POSITIONS[1], shape_form=form)


def _fgw_shape0x_series(ctx: BuildContext) -> Stage:
    """Embryo try E47b: fgw_shape0x's positions with the series rendering of E24a (fgw_shape0_series' expression)."""
    return _fgw(ctx, t_pos=0.0, expression="series", positions=t2_build.POSITIONS[1])


def _fgw_late(ctx: BuildContext) -> Stage:
    """Embryo try E23: fgw's coupling, moves and size, every cell drawn from the later bracket (t_mix = 1)."""
    return _fgw(ctx, t_mix=1.0)


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


def _b0_cfm(ctx: BuildContext) -> Stage:
    return _b0_t2(ctx, cfm=True)


def _b0_trend(ctx: BuildContext) -> Stage:
    """Heart candidate b0_trend (row H11, B0_TREND_RULE): `_b0_t2` with trend=True."""
    return _b0_t2(ctx, trend=True)


def _fgw_herd(ctx: BuildContext) -> Stage:
    return _fgw(ctx, herd=True)


def _t2_extrap_with_start(t: PseudoTask) -> bool:
    return t.task == "T2" and t.kind == "extrap" and "start" in t.inputs


def _heart_extrap_with_start(t: PseudoTask) -> bool:
    """A heart extrapolation task with a `start` stage (b0_trend's tasks; the setting is named, so a future non-heart
    T2 extrapolation task does not match)."""
    return t.task == "T2" and t.kind == "extrap" and "start" in t.inputs and t.setting == "heart"


def _fgw_blend(ctx: BuildContext) -> Stage:
    return _fgw(ctx, expression="blend")


def _fgw_aligned(ctx: BuildContext) -> Stage:
    return _fgw(ctx, move_coupling="aligned")


def _fgw_shape0_aligned(ctx: BuildContext) -> Stage:
    return _fgw(ctx, t_pos=0.0, move_coupling="aligned")


def _fgw_aligned_series(ctx: BuildContext) -> Stage:
    """Embryo try E25a: try E10's aligned moves with the series rendering of E24a (fgw_series' cells and expression)."""
    return _fgw(ctx, move_coupling="aligned", expression="series")


def _fgw_aligned_series_wide(ctx: BuildContext) -> Stage:
    """Embryo try E27: fgw_aligned_series with each rendered cell's profile smoothed over RENDER_K_WIDE neighbours."""
    return _fgw(ctx, move_coupling="aligned", expression="series_wide")


def _fgw_shape1_series(ctx: BuildContext) -> Stage:
    """Embryo try E25b: the later bracket's geometry (t_pos = 1) with the series rendering of E24a."""
    return _fgw(ctx, t_pos=1.0, expression="series")


def _fgw_late_shape1(ctx: BuildContext) -> Stage:
    """Embryo try E25c: the later bracket's cells alone at their own positions (t_mix = 1, t_pos = 1), at the size
    rule: no cell and no move from the earlier bracket."""
    return _fgw(ctx, t_pos=1.0, t_mix=1.0)


# embryo try E15 (experiments/harness/TRIES.md; registered before its run, and run only if its target-free gate, the
# diagnostic N15, passes): every number below is fixed in the row, one value each, no sweep (a swept value is a try).
NBHD_K = 15  # cells per neighbourhood, self included, as the scorer's (T2 metrics `neighborhood_mmd`)
NBHD_RADIUS_FACTOR = 2.0  # how far a cell may end from where it started, in neighbourhood radii: try E15's fixed radius
NBHD_SWEEPS = 30  # rounds of the swap chain, a size of the computation (the mixing's plateau, tests/test_nbhd.py)
NBHD_STREAM = 20  # the swap chain's random stream, apart from the interpolation's and expression steps' 0-9


def nbhd_rearrangement(coords: np.ndarray, labels: np.ndarray, seed: int) -> np.ndarray:
    """E15's re-arrangement of one draw: `at` such that `coords[at]` holds the positions re-arranged among the drawn
    cells within coarse type and the row's fixed radius (vec.t2_build.local_type_permutation at NBHD_RADIUS_FACTOR x the
    draw's median 15-nearest-cell reach, NBHD_SWEEPS rounds, stream [seed, NBHD_STREAM]). The audit N15 (its reported
    after_over_built) and the candidate both call this, so the audit measures exactly the candidate's own step."""
    radius = NBHD_RADIUS_FACTOR * t2_build.neighbourhood_radius(coords, NBHD_K)
    return t2_build.local_type_permutation(coords, labels, radius, NBHD_SWEEPS,
                                           np.random.default_rng([seed, NBHD_STREAM]))


def _fgw_perm(ctx: BuildContext) -> Stage:
    return _fgw(ctx, local_perm=True)


def _fgw_midpoint(ctx: BuildContext) -> Stage:
    return _fgw(ctx, draw=t2_build.DRAWS[1])


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


# ------------------------------------------------------------------------------------------------ T1 tries T12a, T12b
# The strong-null down-tail prior's downside, measured before any test build (experiments/harness/TRIES.md, rows T12a
# and T12b, their constants in vec.t1_scope): the prior is applied by the same code as the test build
# (`vec.t1_scope.apply_down_prior`, the same gene budget, band and thinning rule) and with the same selection: genes
# ordered by the reference's level alone (`DownPriorTable.level_only`, nothing fitted on any truth). T12a:
# `b0_strong_null` on the B0 builds of the two E9.5 -> E10.5 tasks (a strong-null truth). T12b: `scope_strong_null` on
# B0's chamber copy of the E8.5 -> E9.5 task (a weak-null truth).
T12_GROUPS = ("CM", "Endocardium", "Endothelium", "Epicardium")  # the groups with a Feng direction at test (N12): the
# cells the test build moves along Feng's direction, and so the ones the prior may touch (R1 holds the others)


def _t12_labels(ctx: BuildContext) -> frozenset[str]:
    """The reference stage's labels of T12's cells: those whose group by experiments/t1_m1/feng_dirs.label_groups (the
    config's regexes, as build_t1 maps them) is one of T12_GROUPS. The groups are named here because a candidate may not
    read Feng's statistics (and a backtest's Feng cells lie inside its external window)."""
    d = REPO_ROOT / "experiments" / "t1_m1"
    if str(d) not in sys.path:
        sys.path.insert(0, str(d))
    import feng_dirs

    labels = ctx.inputs[ctx.task.reference_role].obs["celltype"].astype(str).unique()
    return frozenset(lab for lab, g in feng_dirs.label_groups(labels).items() if g in T12_GROUPS)


def _t12_prior(ctx: BuildContext):
    """The prior of tries T12a and T12b: the test build's own selection (the reference's level alone, no table) on the
    cells of T12_GROUPS."""
    from vec import t1_scope

    return t1_scope.DownPrior(table=t1_scope.DownPriorTable.level_only(), covered_labels=_t12_labels(ctx))


def _b0_strong_null(ctx: BuildContext) -> Stage:
    return _b0_t1(B0_DAMPING, prior=_t12_prior)(ctx)


def _scope_strong_null(ctx: BuildContext) -> Stage:
    return _b0_t1(0.0, prior=_t12_prior)(ctx)


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

# ------------------------------------------------------------------------------------------------ T3 tries K4 and K4c
# Fixed in their experiments/harness/TRIES.md rows (registered before any fold is read): K2a's factor estimator on each
# fold's own source (t3_factor_raw, the base arm), its gene-wise empirical-Bayes shrinkage (t3_factor_eb, try K4; the
# estimator is experiments/t3_decomp/eb_factor.py) and the control that permutes the shrinkage weights within the same
# deciles (t3_factor_eb_ctrl, try K4c), each built as K2a builds its factor. Decided twice by `run.py decide`, against
# the raw base and against the control; both must keep on both counted cross-knockout folds. The replicate fold is
# reported and never counted (`decide` excludes replicate folds).
K4_FOLDS = {  # pseudo-task key -> its source: knockout roles and files, matched WT file, rules.yaml presets
    "T3/Gata4->Mab21l2": {"ko_roles": ("ko_train", "ko_train_2"), "source_files": K2_FACTOR_FILES[:2],
                          "wt_file": K2_FACTOR_FILES[2], "source_preset": "gata4", "reference_preset": "mab21l2",
                          "counted": True},
    "T3/Mab21l2->Gata4rep2": {"ko_roles": ("ko_train",), "source_files": ("T3/E9.5_mab21l2_ko.h5ad",),
                              "wt_file": "T2/E9.5.h5ad", "source_preset": "mab21l2", "reference_preset": "gata4",
                              "counted": True},
    "T3/Gata4rep1->rep2": {"ko_roles": ("ko_train",), "source_files": K2_FACTOR_FILES[:1],
                           "wt_file": K2_FACTOR_FILES[2], "source_preset": "gata4", "reference_preset": "gata4",
                           "counted": False},
}
K4_CTRL_SEED = 20261109  # K4c: one permutation of the shrinkage weights, this seed (its row)
K4_BINS = 10  # deciles of the factor's denominator
K4_STREAM = K2_STREAM  # all three arms draw as K2a: t3_factor_raw on the Gata4 -> Mab21l2 fold is K2a's build
_K4_PREP: dict = {}


def _t3_eb():
    """experiments/t3_decomp/eb_factor.py, imported as `_t3_head` imports build_submission."""
    _t3_head()
    import eb_factor

    return eb_factor


def t3_source_frame(wt: ad.AnnData, kos: list[ad.AnnData], rules: Mapping, preset: str,
                    n_jobs: int = 1) -> SimpleNamespace:
    """What K2a's factor estimator computes for a fold's source, for any rules.yaml preset (t3_gata4_factor's steps, in
    its order, so `phi` is byte-identical to it on the Gata4 source), and the section blocks K4's jackknife needs. WT
    `wt`'s own labels in the preset's label mode; each knockout file's cells (littermate wild types split off) called
    into that vocabulary by the frozen harmoniser, depth-matched as there; all knockout cells pooled in file order;
    gene_factor over the strict non-Mesp1 labels with at least min_cells cells on both sides, weighted by the WT
    shares. Returns phi (K2a's factor), den, mu_wt, weights (V.p), mask, the pooled counts C and labels lab, block (each
    pooled cell's section: eb_factor.block_ids over the files' kept cells, replicate-by-level ids when files are
    pooled), n_levels, eps and clip. wt and kos as the harness loads them; none of them is modified."""
    t3d = _t3_head().t3d
    hp, fc = rules["harmonisation"], rules["components"]["factor"]
    w = _shallow(wt)
    t3d.add_label_column(w, rules["presets"][preset]["label_mode"])
    V, wlab = _t3_vocab(t3d, w, rules)
    Cw = t3d.dense_counts(w)
    harm = t3d.Harmoniser(w, "t3_label", t3d.label_groups(w, "t3_label"), hp, n_jobs)
    Cs, labs, zs = [], [], []
    for a in kos:
        a = _shallow(a)
        littermates, _, _ = t3d.littermate_split(a)
        if littermates.any():
            a = a[~littermates].copy()
        calls, _ = harm.calls(a, seed=hp["seed"])
        Cs.append(t3d.dense_counts(a))
        labs.append(t3d.calls_to_index(V, calls))
        zs.append(np.asarray(a.obsm[COORD_KEY], dtype=np.float64)[:, 2])
    C, lab = np.vstack(Cs), np.concatenate(labs)
    block, n_levels = _t3_eb().block_ids(zs)
    ok = ((np.bincount(lab[lab >= 0], minlength=V.L) >= fc["min_cells"])
          & (np.bincount(wlab, minlength=V.L) >= fc["min_cells"]))
    mask = ok & V.nm
    mu_wt = t3d.props_type_means(Cw, wlab, V.L)
    phi, den = t3d.gene_factor(t3d.props_type_means(C, lab, V.L), mu_wt, V.p, mask, fc["eps"], tuple(fc["clip"]))
    return SimpleNamespace(phi=phi, den=den, mu_wt=mu_wt, weights=V.p, mask=mask, C=C, lab=lab, block=block,
                           n_levels=n_levels, eps=fc["eps"], clip=tuple(fc["clip"]))


def t3_k4_fit(frame: SimpleNamespace, seed: int = K4_CTRL_SEED, bins: int = K4_BINS):
    """K4's and K4c's factors from a fold's source frame (`t3_source_frame`): eb_factor.jackknife, then
    eb_factor.eb_factors."""
    eb = _t3_eb()
    jack = eb.jackknife(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask, frame.eps)
    return eb.eb_factors(frame.phi, jack, clip=frame.clip, seed=seed, bins=bins)


def _k4_prepared(ctx: BuildContext) -> dict:
    """The fold's raw, EB and control factors and the fold's WT in t3_decompose's form, once per loaded fold (none
    depends on the seed; the harness hands every seed the same input objects)."""
    fold = K4_FOLDS[ctx.task.key]
    roles = ("wt_train", *fold["ko_roles"], ctx.task.reference_role)
    objs = tuple(ctx.inputs[r] for r in roles)
    hit = _K4_PREP.get("fold")
    if hit is None or hit["key"] != ctx.task.key or any(a is not b for a, b in zip(hit["inputs"], objs, strict=True)):
        rules = _t3_head().load_rules()
        wt, pool = ctx.inputs["wt_train"], ctx.inputs[ctx.task.reference_role]
        if list(wt.var_names) != list(pool.var_names):
            raise ValueError(f"{ctx.task.key}: the factor's WT and the reference hold different genes")
        fit = t3_k4_fit(t3_source_frame(wt, [ctx.inputs[r] for r in fold["ko_roles"]], rules, fold["source_preset"]))
        _K4_PREP.clear()
        hit = _K4_PREP["fold"] = {"key": ctx.task.key, "inputs": objs, "raw": fit.raw, "eb": fit.eb, "ctrl": fit.ctrl,
                                  "ref": t3_reference(pool, fold["reference_preset"], rules)}
    return hit


def _t3_factor_k4(arm: str) -> Callable[[BuildContext], Stage]:
    """Tries K4 and K4c and their raw base (arm "eb", "ctrl", "raw"): the fold's arm factor, at exponent 1, on the
    reference WT's real cells at its own composition, built as K2a builds its factor (`_t3_factor_xko`), from one
    stream for all three arms."""

    def build(ctx: BuildContext) -> Stage:
        bs = _t3_head()
        prep = _k4_prepared(ctx)
        ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
        built = bs.construct(ref, count_for(ctx.board, pool.n_obs), bs.t3d.stable_rng(ctx.seed, K4_STREAM), f=ref.V.p,
                             phi=prep[arm])
        X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes),
                                   bs.t3d.stable_rng(ctx.seed, f"{K4_STREAM}:assemble"))
        return finish(ctx, pool, X, coords if ctx.board.needs_coords else None)

    return build


def _k4_applies(t: PseudoTask) -> bool:
    """K4's arms run on the two counted cross-knockout folds and, reported only, on the replicate fold, each with its
    own source (K4_FOLDS); `decide` never counts the replicate fold."""
    fold = K4_FOLDS.get(t.key)
    return (fold is not None and t.task == "T3"
            and tuple(t.inputs.get(r) for r in (*fold["ko_roles"], "wt_train"))
            == (*fold["source_files"], fold["wt_file"]))


K4_PARAMS = {"factor": "K2a's estimator on the fold's own source against its matched WT (t3_source_frame): strict "
                       "non-Mesp1 labels, gene_factor (rules.yaml components.factor), exponent 1",
             "estimator": "experiments/t3_decomp/eb_factor.py", "bins": K4_BINS, "stream": K4_STREAM,
             "build": "build_submission.construct and assemble, as K2a; no depth, no composition"}

_ANY = (frozenset(),)
_LP_PARAMS = {"bracket": list(LP_BRACKET), "groups": list(LP_GROUPS), "fraction": "the step's share of the bracket "
              "(target day - reference day) / 4; no momentum, no other change",
              "held_flat": "the pre-registered list (named genes alone where its pool rule's stages are not inputs)",
              "scope": "chamber"}


def _lp(name: str, what: str, feng_groups: str = "authors", **build) -> Candidate:
    """A T1 late-proxy try (`_late_proxy`), with its external bracket and the hand-over declared. feng_groups "authors"
    (the default: every registered candidate, unchanged) or "recomputed" (Feng's permitted-only groups; no candidate
    uses it until a row registers one)."""
    if feng_groups == "authors":
        return Candidate(name, frozenset({"T1"}), _ANY, _late_proxy(**build), {"what": what, **_LP_PARAMS},
                         externals=(LP_EXTERNAL,), derived=("t1_late_proxy",), derived_optional=("t1_held_flat",))
    return Candidate(name, frozenset({"T1"}), _ANY, _late_proxy(feng_groups=feng_groups, **build),
                     {"what": what, **_LP_PARAMS, "feng_groups": feng_groups}, externals=(LP_EXTERNAL,),
                     derived=(LP_DERIVED[feng_groups],), derived_optional=("t1_held_flat",))


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
    # embryo try E24 (experiments/harness/TRIES.md, registered before its run): the later bracket's run profile
    "fgw_series": Candidate("fgw_series", frozenset({"T2"}), (), _fgw_series,
                            {"what": "fgw's draws and positions with the shallower bracket's cells moved by a per-gene "
                                     "factor measured on the extraembryonic cell types both brackets label alike, then "
                                     "to the deeper bracket's depth", "rule": t2_build.SERIES_RULE, "scope": "whole"},
                            requires=_series_applies),
    "fgw_shape0_hi": Candidate("fgw_shape0_hi", frozenset({"T2"}), (), _fgw_shape0_hi,
                               {"what": "fgw_shape0 (the earlier bracket's geometry) with every cell at the deeper "
                                        "bracket's depth: the shape cover of today's s2 and the depth top-up in one "
                                        "file", "rule": t2_build.RENDER_RULE, "scope": "whole"}, requires=_interp),
    "fgw_shape0_series": Candidate("fgw_shape0_series", frozenset({"T2"}), (), _fgw_shape0_series,
                                   {"what": "fgw_shape0 (the earlier bracket's geometry) with the series rendering",
                                    "rule": t2_build.SERIES_RULE, "scope": "whole"}, requires=_series_applies),
    # embryo tries E47a and E47b (experiments/harness/TRIES.md, registered before their run): fgw_shape0 and
    # fgw_shape0_series with the earlier bracket's geometry exact, a sample of its own positions
    "fgw_shape0x": Candidate("fgw_shape0x", frozenset({"T2"}), (), _fgw_shape0x,
                             {"what": "fgw_shape0's mixture of real cells, every cell on a real position of the "
                                      "earlier bracket: its own cells where they are, the later bracket's cells on "
                                      "positions of earlier cells not drawn, assigned one to one by the exact minimum "
                                      "squared distance in fgw's joint expression PCs; scaled to the embryo size rule",
                              "positions": t2_build.POSITIONS[1], "stream": t2_build.FLOOR_EXACT_STREAM,
                              "scope": "whole"}, requires=_interp),
    "fgw_shape0x_series": Candidate("fgw_shape0x_series", frozenset({"T2"}), (), _fgw_shape0x_series,
                                    {"what": "fgw_shape0x (the earlier bracket's own positions) with the series "
                                             "rendering", "rule": t2_build.SERIES_RULE,
                                     "positions": t2_build.POSITIONS[1], "stream": t2_build.FLOOR_EXACT_STREAM,
                                     "scope": "whole"}, requires=_series_applies),
    # row E49 (as drafted): fgw_shape0x with a shape form between the brackets, three arms and four controls
    **{name: Candidate(name, frozenset({"T2"}), (), _fgw_shape0x_form(form),
                       {"what": "fgw_shape0x (the earlier bracket's own positions) with the shape form "
                                f"vec.shape_interp {form} applied before the size rule",
                        "positions": t2_build.POSITIONS[1], "stream": t2_build.FLOOR_EXACT_STREAM, "shape_form": form,
                        "scope": "whole"}, requires=_interp)
       for name, form in SHAPE_CANDIDATES.items()},
    # embryo tries E25a-E25c (experiments/harness/TRIES.md, registered before their run): the geometry that goes with
    # the series rendering, and the later bracket alone
    "fgw_aligned_series": Candidate("fgw_aligned_series", frozenset({"T2"}), (), _fgw_aligned_series,
                                    {"what": "fgw_series (the series rendering) with the moves taken from the aligned "
                                             "coupling of try E10 (vec.core.aligned_coupling in the rotated common "
                                             "frame; FGW's rotation kept)", "rule": t2_build.SERIES_RULE,
                                     "alpha": t2_build.INTERP_ALPHA, "scope": "whole"}, requires=_series_applies),
    # embryo try E27 (experiments/harness/TRIES.md, registered before its run): the wider smoothing of N24b
    "fgw_aligned_series_wide": Candidate("fgw_aligned_series_wide", frozenset({"T2"}), (), _fgw_aligned_series_wide,
                                         {"what": "fgw_aligned_series with each rendered cell's profile taken from "
                                                  "its RENDER_K_WIDE nearest cells instead of RENDER_K",
                                          "rule": t2_build.SERIES_RULE, "k": t2_build.RENDER_K_WIDE,
                                          "alpha": t2_build.INTERP_ALPHA, "scope": "whole"},
                                         requires=_series_applies),
    "fgw_shape1_series": Candidate("fgw_shape1_series", frozenset({"T2"}), (), _fgw_shape1_series,
                                   {"what": "fgw_shape1 (the later bracket's geometry, t_pos = 1) with the series "
                                            "rendering", "rule": t2_build.SERIES_RULE, "scope": "whole"},
                                   requires=_series_applies),
    "fgw_late_shape1": Candidate("fgw_late_shape1", frozenset({"T2"}), (), _fgw_late_shape1,
                                 {"what": "the later bracket's cells alone at their own positions (t_mix = 1, "
                                          "t_pos = 1), scaled to the embryo size rule", "scope": "whole"},
                                 requires=_interp),
    # embryo try E23 (experiments/harness/TRIES.md, registered before its run): the later bracket's cells only
    "fgw_late": Candidate("fgw_late", frozenset({"T2"}), (), _fgw_late,
                          {"what": "fgw's coupling, moves and size rule with every cell drawn from the later bracket "
                                   "(t_mix = 1): the same-series bracket alone, at its moved positions",
                           "scope": "whole"}, requires=_interp),
    # embryo tries E19a-E19c (experiments/harness/TRIES.md, registered before their run): one depth for every cell
    **{f"fgw_render_{m}": Candidate(f"fgw_render_{m}", frozenset({"T2"}), (), _fgw_render(m),
                                    {"what": "fgw's draws and positions with every cell moved to one depth: the "
                                             "deeper bracket thinned, the shallower topped up from kNN-smoothed "
                                             "profiles", "mode": m, "rule": t2_build.RENDER_RULE, "scope": "whole"},
                                    requires=_interp) for m in t2_build.RENDER_MODES},
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
    # embryo try E12 (experiments/harness/TRIES.md, registered before its run): the s2 shape variant fgw_shape0 with
    # try E10's aligned moves, so its mapped later bracket lands on the side position says. Both embryo tasks were
    # selected on when fgw_shape0 was chosen, so neither counts for it under the adoption rule (decide excludes them)
    "fgw_shape0_aligned": Candidate("fgw_shape0_aligned", frozenset({"T2"}), (), _fgw_shape0_aligned,
                                    {"what": "fgw_shape0 (the earlier bracket's geometry, t_pos = 0) with the moves "
                                             "taken from the aligned coupling of try E10 (vec.core.aligned_coupling "
                                             "in the rotated common frame; FGW's rotation kept)",
                                     "alpha": t2_build.INTERP_ALPHA, "scope": "whole"}, requires=_interp,
                                    tuned_on=frozenset({"T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25"})),
    # embryo try E15 (experiments/harness/TRIES.md, registered before its run; it runs only if its target-free gate,
    # the diagnostic N15, passes): fgw's cells and expression rows, positions re-arranged among them
    "fgw_perm": Candidate("fgw_perm", frozenset({"T2"}), (), _fgw_perm,
                          {"what": "fgw with the positions re-arranged among the drawn cells within coarse type "
                                   "(the lineage map's group) and a fixed radius (vec.t2_build.local_type_permutation)",
                           "k": NBHD_K, "radius_factor": NBHD_RADIUS_FACTOR, "sweeps": NBHD_SWEEPS,
                           "scope": "whole"}, requires=_interp),
    # embryo try E13 (experiments/harness/TRIES.md, registered before its run; it runs only if its target-free gate N16
    # passes): fgw's cells and positions; which real cells are drawn is tilted towards the local displacement midpoint
    "fgw_midpoint": Candidate("fgw_midpoint", frozenset({"T2"}), (), _fgw_midpoint,
                              {"what": "fgw with each bracket's real cells drawn by a tilt towards the local "
                                       "displacement midpoint (w = exp(-s^2 / (2 var(s))), s the cell's PC offset "
                                       "from the midpoint projected on the local bracket-difference axis, within "
                                       "k-means strata of the joint PCs; every cell unedited at its moved position, "
                                       "exact bracket fractions)",
                               "k_strata": t2_build.MIDPOINT_K_STRATA, "k_anchors": t2_build.MIDPOINT_K_ANCHORS,
                               "kmeans_seed": t2_build.MIDPOINT_KMEANS_SEED,
                               "streams": list(t2_build.MIDPOINT_STREAMS), "scope": "whole"}, requires=_interp),
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
    # T1 try T11 (experiments/harness/TRIES.md, registered before its run): a mean-free within-type momentum on B0
    "g_mom": Candidate("g_mom", frozenset({"T1"}), (), _g_mom,
                       {"what": "B0 plus a mean-free within-type momentum per reference-stage metacell, from an OT "
                                "between the start stage and the reference (vec.t1_state_transport); R1's held types "
                                "keep B0's cells", "damping": B0_DAMPING, "scope": "chamber"},
                       requires=_b0_applies, derived=("t1_held_flat",)),
    # T1 tries T12a and T12b (experiments/harness/TRIES.md, registered before their run): the strong-null down-tail
    # prior
    "b0_strong_null": Candidate("b0_strong_null", frozenset({"T1"}), (), _b0_strong_null,
                                {"what": "T1 try T12a: B0 plus the strong-null down-tail prior, genes chosen by the "
                                         "reference level alone (the test build's selection, no table)",
                                 "damping": B0_DAMPING, "scope": "chamber"},
                                requires=_b0_applies, derived=("t1_held_flat",)),
    "scope_strong_null": Candidate("scope_strong_null", frozenset({"T1"}), (), _scope_strong_null,
                                   {"what": "T1 try T12b: B0's chamber copy plus the same prior, genes chosen by the "
                                            "reference level alone (a weak-null truth)",
                                    "scope": "chamber"}, requires=lambda t: "start" not in t.inputs),
    # T1 tries T20 and T21 (experiments/harness/TRIES.md, registered before their run): g_mom's placebo arms, on g_mom's
    # own draws, base and fit; judged by g_mom against each of them (a paired mean gain above 0.3 and 8 of 10 seeds on
    # both E10.5 cuts, experiments/t1_scope/g_mom_verify.py placebo)
    "g_mom_flip": Candidate("g_mom_flip", frozenset({"T1"}), (), _g_mom_flip,
                            {"what": "try T20, a placebo: g_mom with every within-type log factor's sign reversed "
                                     "(vec.t1_state_transport.flip_factors)", "damping": B0_DAMPING,
                             "scope": "chamber"}, requires=_b0_applies, derived=("t1_held_flat",)),
    "g_mom_perm": Candidate("g_mom_perm", frozenset({"T1"}), (), _g_mom_perm,
                            {"what": "try T21, a placebo: g_mom with each category's metacell factor rows re-assigned "
                                     "among its metacells by one fixed derangement (vec.t1_state_transport."
                                     "permute_factors)", "permutation_seed": G_MOM_PERM_SEED, "damping": B0_DAMPING,
                             "scope": "chamber"}, requires=_b0_applies, derived=("t1_held_flat",)),
    # T1 tries T23 and T24 (experiments/harness/TRIES.md, registered before their run): g_mom's no-development control
    # (the same pipeline between two halves of the reference stage, at g_mom's size) and its pseudobulk-held arm, on
    # g_mom's own draws, base and stream (experiments/t1_scope/g_mom_verify.py placebo)
    "g_mom_null": Candidate("g_mom_null", frozenset({"T1"}), (), _g_mom_null,
                            {"what": "try T23, a no-development control: g_mom's draws and stream with factors from "
                                     "the same pipeline fitted between two random halves of the reference stage "
                                     "(cross-fitted both ways), per category at g_mom's own root-mean-square size "
                                     "(vec.harness.g_mom_no_development_factors)", "split_seed": G_MOM_NULL_SEED,
                             "damping": B0_DAMPING, "scope": "chamber"}, requires=_b0_applies,
                            derived=("t1_held_flat",)),
    "g_mom_held": Candidate("g_mom_held", frozenset({"T1"}), (), _g_mom_held,
                            {"what": "try T24: g_mom with each category's pseudobulk held per gene (rescale=True of "
                                     "vec.t1_state_transport.apply_state_factors, T9's construction), the same "
                                     "factors, draws and stream", "damping": B0_DAMPING, "scope": "chamber"},
                            requires=_b0_applies, derived=("t1_held_flat",)),
    # heart try H3 (experiments/harness/TRIES.md, registered before its run)
    "b0_depth": Candidate("b0_depth", frozenset({"T2"}), (), _b0_depth,
                          {"what": "B0 on a T2 extrapolation board with its cells binomially thinned to the target's "
                                   "depth, extrapolated at B0's damping (half the latest log rate per day over the "
                                   "horizon)", "damping": B0_DAMPING, "rule": t2_build.DEPTH_EXTRAP_RULE,
                           "scope": "whole"}, requires=_t2_extrap_with_start),
    # heart try H6 (experiments/harness/TRIES.md, registered before its run and run only if its target-free gate N17
    # passes): each lineage's momentum filtered by the lineage's own within-stage covariance
    "b0_cfm": Candidate("b0_cfm", frozenset({"T2"}), (), _b0_cfm,
                        {"what": "B0 on a T2 extrapolation board with each lineage's momentum projected onto the "
                                 "directions of that lineage's own within-stage covariance above the Marchenko-Pastur "
                                 "edge and renormalised to its unfiltered L2 norm", "damping": B0_DAMPING,
                         "rule": t2_build.CFM_RULE, "scope": "whole"}, requires=_t2_extrap_with_start),
    # heart candidate b0_trend (row H11 of experiments/harness/TRIES.md; both heart files by the captain's ruling of
    # 2026-11-13, docs/PLAN.md decision log (w))
    "b0_trend": Candidate("b0_trend", frozenset({"T2"}), (), _b0_trend,
                          {"what": "B0 with each lineage's per-gene velocity the least-squares trend of log(mean CP10k "
                                   "+ eps) over every input stage where it has min_cells cells, read over B0's time "
                                   "budget (heart only; count space)", "damping": B0_DAMPING,
                           "rule": "vec.harness.B0_TREND_RULE", "scope": "chamber"},
                          requires=_heart_extrap_with_start),
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
    # T3 tries K4 and K4c and their raw base (experiments/harness/TRIES.md, registered before their run): on the two
    # counted cross-knockout folds (and, reported only, the replicate fold), K2a's factor estimator, its empirical-Bayes
    # shrinkage, and the shrinkage-matched control
    "t3_factor_raw": Candidate("t3_factor_raw", frozenset({"T3"}), (), _t3_factor_k4("raw"),
                               {"what": "K4's base arm: K2a's factor estimator on the fold's own source, on the "
                                        "reference WT's real cells at its own composition", **K4_PARAMS},
                               requires=_k4_applies),
    "t3_factor_eb": Candidate("t3_factor_eb", frozenset({"T3"}), (), _t3_factor_k4("eb"),
                              {"what": "try K4: the raw factor, each gene's log ratio shrunk towards f = 1 by its "
                                       "delete-one-section jackknife instability (DerSimonian-Laird within deciles "
                                       "of the factor's denominator, numerator side), then rescaled to the raw "
                                       "factor's robust SD by one constant", **K4_PARAMS}, requires=_k4_applies),
    "t3_factor_eb_ctrl": Candidate("t3_factor_eb_ctrl", frozenset({"T3"}), (), _t3_factor_k4("ctrl"),
                                   {"what": "try K4c: K4 with its shrinkage weights permuted across genes within the "
                                            "same deciles (K2b's logic, one permutation), the same constant",
                                    **K4_PARAMS, "permutation_seed": K4_CTRL_SEED}, requires=_k4_applies),
    # T1 tries T13-T18 (experiments/harness/TRIES.md, registered before their run): the late-proxy tournament
    "lp_proxy": _lp("lp_proxy", "try T13: Feng's E14.5 -> E18 direction at the step's share of its bracket, on "
                    "B0's chamber copy (the stand-in for M1's real E9.5 -> E14.5 direction)"),
    "lp_cflip": _lp("lp_cflip", "try T14, control Cflip: lp_proxy's direction negated", flip=True),
    "lp_cperm": _lp("lp_cperm", "try T15, control Cperm: lp_proxy's direction permuted across genes within deciles of "
                    "the group's mean CP10k (one permutation, seed 20261110)", permute=True),
    "lp_hyg": _lp("lp_hyg", "try T16: lp_proxy with R2's classes D and A1 held flat", hold=("D", "A1")),
    "lp_hyg_n": _lp("lp_hyg_n", "try T17: lp_proxy with R2's classes D, A1, N1 and N2 held flat (N1 and N2, the lines "
                    "the captain dropped)", hold=("D", "A1", "N1", "N2")),
    "lp_eb": _lp("lp_eb", "try T18: lp_proxy with empirical-Bayes reliability shrinkage per expression decile (from "
                 "Feng's two per-cell halves)", shrink=True),
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


def _within_rows(a: ad.AnnData, file: str, category: str) -> np.ndarray:
    """Rows of a released T1 stage whose chamber category (configs/t1_chamber_classes.yaml, by the file's own per-file
    rules) is `category`: a within-type pseudo-task cuts its target, reference and inputs the same way."""
    from vec import t1_scope

    return np.flatnonzero(t1_scope.classify(a, source=stem(file), key_added=None) == category)


def filter_target(t: PseudoTask, data_dir: Path, genes) -> tuple[Stage, dict]:
    """The released T1 target cut by the task's target filter (fixed in the config before any result): the cells of
    the kept categories of configs/t1_chamber_classes.yaml, plus a share of contamination cells drawn at the filter's
    seed from the named categories; or, for a within-type read (`within`: a category), only that category's cells.
    Returns the scoring Stage and a record (counts only)."""
    from vec import t1_scope

    spec = t1_scope.load_spec()
    a = read_stage(data_dir / t.target, genes)
    f = t.target_filter
    if f.get("within") is not None:
        rows = _within_rows(a, t.target, f["within"])
        return Stage.from_anndata(a[rows], genes), {"kept": int(len(rows)), "within": f["within"], "of": int(a.n_obs)}
    t1_scope.classify(a, source=stem(t.target))
    cat = np.asarray(a.obs[t1_scope.CATEGORY_KEY]).astype(str)
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
    within = (t.target_filter or {}).get("within")  # a within-type read cuts the reference and the inputs as the target
    if t.target_filter is not None:
        target, filt = filter_target(t, data_dir, genes)
    else:
        target, filt = Stage.from_anndata(read_stage(data_dir / t.target, genes), genes, need_coords=need_coords), None
    ref_cells = read_stage(data_dir / t.reference, genes)
    if within is not None:
        ref_cells = ref_cells[_within_rows(ref_cells, t.reference, within)]
    ref = Stage.from_anndata(ref_cells, genes, need_coords=need_coords)
    board = task_board(t, boards, target.n, ref.n)
    if t.task == "T1":
        inputs = {r: load_stage(data_dir / f, lineage=False) for r, f in t.inputs.items()}
        if within is not None:
            inputs = {r: a[_within_rows(a, t.inputs[r], within)].copy() for r, a in inputs.items()}
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
