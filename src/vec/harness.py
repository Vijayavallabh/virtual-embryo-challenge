"""Pseudo-task harness (feat-008): paired backtests of candidates on released stages, decided by vec.adoption.

A pseudo-task (configs/pseudo_tasks.yaml) stands in for a test board's geometry using released stages only (the
validation answers included). For each scoring seed s in vec.adoption.SEEDS the harness
  1. draws the calibrated bake protocol's working copies (vec.evaluate, protocol="bake"): the target's 10% copy and
     its truth half, and the reference's 10% copy, which is both the floor and the DE reference. The floor therefore
     scores exactly 50, and the harness checks that it does;
  2. computes the floor and ceiling anchors once and scores every candidate against them, with the scorer build
     "main" (docs/PLAN.md section 3: candidates are compared without the pre-fix build's parity coin; "prefix" is for
     calibration only, and selecting on it would be the parity hedge);
  3. builds every candidate from the pseudo-task's inputs only, with the build seed BUILD_SEED_OFFSET + s. All
     candidates share it (common random numbers), and it is never a scoring stream: numpy's default_rng([s, 0]),
     which vec.t2_build.sample_cells draws rows with, is the same stream as default_rng(s), the bake protocol's
     reference draw, so a candidate built at seed s would contain the floor's own cells (checked 2026-11-01).
Adoption goes through vec.adoption.keep_change over the pseudo-tasks of a test board on which both candidates apply
(`decide`). Results are aggregate numbers per seed; registry rows go to results/registry.csv (vec.registry).

Leakage rules (`check_task`, applied when the config loads): the target is never an input; extrapolation inputs are
all earlier than the target; interpolation inputs include a stage on each side of the target and none strictly
inside the bracket; a task's `forbidden` files are never inputs; and candidates receive only the loaded input stages,
never a path. A candidate must not read files derived from a pseudo-task's target either (the size JSONs measured
E10.5 and E7.5; the count decision was decided on them): the built-in candidates read nothing but their inputs.
"""

from __future__ import annotations

import hashlib
import subprocess
import time
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import anndata as ad
import numpy as np
import yaml

from vec import t2_build
from vec.adoption import SEEDS, keep_change
from vec.boards import REPO_ROOT, VALIDATION_BOARDS, Board, load_boards
from vec.core import draw_balanced
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample
from vec.io import COORD_KEY, dense, read_stage
from vec.prep import load_stage

CONFIG = REPO_ROOT / "configs" / "pseudo_tasks.yaml"
DATA_DIR = REPO_ROOT / "data" / "raw"
BUILD = "main"
PROTOCOL = "bake"
BUILD_SEED_OFFSET = 1_000_000
FLOOR = "floor"
TASKS = ("T1", "T2", "T3")
KINDS = ("extrap", "interp", "ko")

# E-day of every released file a pseudo-task may name (docs/DATA_NOTES.md; E8.25_late is about E8.3, Data Protocol).
STAGE_DAYS = {
    "T1/E8.5_RNA.h5ad": 8.5, "T1/E9.5_RNA.h5ad": 9.5, "T1/E10.5_RNA.h5ad": 10.5,
    "T2/E6.75.h5ad": 6.75, "T2/E7.25.h5ad": 7.25, "T2/E7.5.h5ad": 7.5, "T2/E8.0.h5ad": 8.0,
    "T2/E8.25_late.h5ad": 8.3, "T2/E8.5.h5ad": 8.5, "T2/E8.75.h5ad": 8.75, "T2/E9.5.h5ad": 9.5,
    "T2/E10.5.h5ad": 10.5,
    "T3/E8.75_gata4_ko_1.h5ad": 8.75, "T3/E8.75_gata4_ko_2.h5ad": 8.75, "T3/E9.5_mab21l2_ko.h5ad": 9.5,
}
KNOCKOUTS = frozenset(f for f in STAGE_DAYS if f.startswith("T3/"))


class LeakageError(ValueError):
    """A pseudo-task whose inputs could reveal its target."""


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
    notes: str = ""

    @property
    def panel_board(self) -> str:
        """The board whose gene panel and limits rule this pseudo-task uses."""
        return self.validation_board or self.twin_board


def check_task(t: PseudoTask) -> PseudoTask:
    """Refuse a pseudo-task that breaks a leakage rule or names an unknown file (see the module docstring)."""
    if t.task not in TASKS or t.kind not in KINDS:
        raise ValueError(f"{t.key}: task must be one of {TASKS} and kind one of {KINDS}")
    if (t.task == "T2") != (t.setting is not None):
        raise ValueError(f"{t.key}: a T2 pseudo-task names its setting, and only T2 does")
    if not t.panel_board or (t.validation_board and t.validation_board not in VALIDATION_BOARDS) or (
            t.twin_board and t.twin_board not in VALIDATION_BOARDS):
        raise ValueError(f"{t.key}: validation_board or twin_board must name a validation board")
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
    day = STAGE_DAYS[t.target]
    days = [STAGE_DAYS[f] for f in files]
    if t.kind == "extrap" and any(d >= day for d in days):
        raise LeakageError(f"{t.key}: an extrapolation input is not earlier than the target")
    if t.kind == "interp":
        if any(d == day for d in days) or not any(d < day for d in days) or not any(d > day for d in days):
            raise LeakageError(f"{t.key}: interpolation inputs need one stage on each side and none at the target")
        before, after = max(d for d in days if d < day), min(d for d in days if d > day)
        if STAGE_DAYS[t.reference] != before:
            raise ValueError(f"{t.key}: the reference must be the earlier bracket")
        if any(before < d < after for d in days):
            raise LeakageError(f"{t.key}: an input lies inside the bracket")
    if t.kind == "ko" and (t.target not in KNOCKOUTS or t.reference in KNOCKOUTS):
        raise ValueError(f"{t.key}: a knockout pseudo-task targets a knockout file against a wild-type reference")
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
    twin's gene panel and its own counts. Nothing is ever taken from a test board."""
    counts = pseudo_counts(n_target, n_ref)
    if t.validation_board:
        b = boards[t.validation_board]
        off = sorted(k for k, v in counts.items() if getattr(b, k) != v)
        if off:
            raise ValueError(f"{t.key}: the files do not reproduce {b.key}'s published counts ({off})")
        return replace(b, anchors={})
    twin = boards[t.twin_board]
    return replace(twin, key=f"pseudo:{t.key}", split="pseudo", anchors={}, **counts)


# ------------------------------------------------------------------------------------------------ candidates
@dataclass
class BuildContext:
    """What a candidate sees: the pseudo-task, its board, the loaded input stages (roles of the config; duplicates
    dropped on T2/T3; never the target) and the build seed. Candidates must not modify the inputs."""

    task: PseudoTask
    board: Board
    inputs: Mapping[str, ad.AnnData]
    seed: int


@dataclass(frozen=True)
class Candidate:
    """A named builder. `needs` lists alternative sets of input roles; the candidate applies to a pseudo-task of one
    of `tasks` that has every role of at least one set."""

    name: str
    tasks: frozenset[str]
    needs: tuple[frozenset[str], ...]
    build: Callable[[BuildContext], Stage] | None
    params: Mapping = field(default_factory=dict)

    def applies(self, t: PseudoTask) -> bool:
        return t.task in self.tasks and any(n <= set(t.inputs) for n in self.needs)


def base_role(t: PseudoTask) -> str:
    """The input role holding the reference stage's own cells."""
    return next(role for role, f in t.inputs.items() if f == t.reference)


def count_for(board: Board, pool_size: int) -> int:
    """The approved count rule (max_cells, within the board's limits) on any board."""
    return t2_build.n_cells(board, pool_size)[0]


def draw_copy(ctx: BuildContext) -> tuple[ad.AnnData, np.ndarray, np.ndarray, np.ndarray | None]:
    """Real reference cells at the board's count: (pool, rows, X float32 in board gene order, coords or None)."""
    pool = ctx.inputs[base_role(ctx.task)]
    n = count_for(ctx.board, pool.n_obs)
    if ctx.board.needs_coords:
        smp = t2_build.sample_cells(pool, n, seed=ctx.seed)
        rows, coords = smp.rows, smp.coords
    else:
        rows, coords = draw_balanced(pool.n_obs, n, np.random.default_rng([ctx.seed, 0])), None
    X = dense(pool[:, ctx.board.genes].X[rows]).astype(np.float32)
    return pool, rows, X, coords


def finish(ctx: BuildContext, pool: ad.AnnData, X: np.ndarray, coords: np.ndarray | None) -> Stage:
    """A candidate as the scorer reads its file: T1 as float32 X; coordinate boards through the T2 head's guard and
    checked builder (a failing guard or check raises)."""
    X = np.clip(np.asarray(X, dtype=np.float32), 0.0, None)
    if not ctx.board.needs_coords:
        return Stage(X=X)
    guard = t2_build.cloud_guard(coords, np.asarray(pool.obsm[COORD_KEY], dtype=float)[:, :3], seed=ctx.seed)
    if not guard.ok:
        raise ValueError(f"{ctx.task.key}: cloud guard flags {guard.flags}")
    sub, rep = t2_build.build_submission(ctx.board, X, coords, genes=list(ctx.board.genes))
    if not rep.ok:
        raise ValueError(f"{ctx.task.key}: the built file fails its check: {rep}")
    return Stage.from_anndata(sub, ctx.board.genes, need_coords=True)


def _copy(ctx: BuildContext) -> Stage:
    pool, _, X, coords = draw_copy(ctx)
    return finish(ctx, pool, X, coords)


def _type_means(a: ad.AnnData, genes) -> dict[str, np.ndarray]:
    X, ct = a[:, genes].X, np.asarray(a.obs["celltype"]).astype(str)
    return {c: np.asarray(X[ct == c].mean(0)).ravel() for c in np.unique(ct)}


def start_role(t: PseudoTask) -> str | None:
    """The stage before the reference: `start` on extrapolation, `earlier` (before the earlier bracket) on
    interpolation, as the organisers' pseudobulk_shift uses them."""
    return "start" if "start" in t.inputs else ("earlier" if t.kind == "interp" and "earlier" in t.inputs else None)


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


_ANY = (frozenset(),)
CANDIDATES = {
    FLOOR: Candidate(FLOOR, frozenset(TASKS), _ANY, None,
                     {"what": "the organisers' floor: the reference's 10% working copy (copy_last or wt_identity)"}),
    "copy": Candidate("copy", frozenset(TASKS), _ANY, _copy, {"what": "real reference cells at the count rule"}),
    "pseudobulk_shift": Candidate("pseudobulk_shift", frozenset({"T1", "T2"}),
                                  (frozenset({"start"}), frozenset({"earlier"})), _pseudobulk_shift,
                                  {"what": "per-type mean change from the stage before, on copied reference cells"}),
    "shift_transfer": Candidate("shift_transfer", frozenset({"T3"}), (frozenset({"ko_train", "wt_train"}),),
                                _shift_transfer, {"what": "training knockout's pseudobulk change on copied WT cells"}),
}


# ------------------------------------------------------------------------------------------------ running
def score_seeds(t: PseudoTask, board: Board, target: Stage, ref: Stage, inputs: Mapping[str, ad.AnnData],
                candidates: list[Candidate], seeds=SEEDS, build: str = BUILD, log=None) -> list[dict]:
    """Per seed: the floor (checked to score 50) and every applicable candidate, all against one seed's anchors."""
    rows = []
    for s in seeds:
        t0 = time.time()
        f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, s), target, ref, seed=s, anchors="computed",
                        build=build, protocol=PROTOCOL)
        if abs(f.score - 50.0) > 1e-9 or f.extras["inverted_anchors"]:
            raise ValueError(f"{t.key} seed {s}: the floor scores {f.score} (inverted {f.extras['inverted_anchors']})")
        row = {"seed": s, "build_seed": BUILD_SEED_OFFSET + s, "scores": {FLOOR: f.score}, "skills": {FLOOR: f.skills},
               "n_cells": {FLOOR: int(round(ref.n * SUBSAMPLE_FRAC))}, "anchors": f.anchors}
        for c in candidates:
            if c.name == FLOOR or not c.applies(t):
                continue
            pred = c.build(BuildContext(task=t, board=board, inputs=inputs, seed=BUILD_SEED_OFFSET + s))
            r = score_board(board, pred, target, ref, seed=s, anchors=f.anchors, build=build, protocol=PROTOCOL)
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


def run_task(t: PseudoTask, candidates: list[Candidate], seeds=SEEDS, data_dir: Path = DATA_DIR,
             boards: Mapping[str, Board] | None = None, build: str = BUILD, log=print) -> dict:
    """Load one pseudo-task's stages and score the candidates on it (DGX: the files live under data/raw)."""
    commit, clean = _git("rev-parse", "HEAD"), _git("status", "--porcelain", "--untracked-files=no") == ""
    boards = boards if boards is not None else load_boards()
    check_index(boards)
    check_task(t)
    genes = boards[t.panel_board].genes
    need_coords = boards[t.panel_board].needs_coords
    target = Stage.from_anndata(read_stage(data_dir / t.target, genes), genes, need_coords=need_coords)
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
    applicable = [c.name for c in candidates if c.applies(t)]
    rows = score_seeds(t, board, target, ref, inputs, candidates, seeds=seeds, build=build, log=log)
    files = sorted({t.target, t.reference, *t.inputs.values()})
    return {
        "task": asdict(t) | {"inputs": dict(t.inputs), "forbidden": list(t.forbidden)},
        "board": {k: getattr(board, k) for k in ("key", "task", "setting", "split", "needs_coords", "min_cells",
                                                   "max_cells", "truth_cells", "ref_cells", "floor_model")},
        "candidates": applicable, "params": {c.name: dict(c.params) for c in candidates if c.name in applicable},
        "protocol": PROTOCOL, "build": build, "anchors": "computed", "frac": SUBSAMPLE_FRAC,
        "seeds": list(seeds), "rows": rows,
        "provenance": {"commit": commit, "tree_clean": clean, "config_sha256": _sha256(CONFIG),
                       "files_sha256": {f: _sha256(data_dir / f) for f in files},
                       "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())},
    }


def decide(results: list[dict], base: str, change: str) -> dict:
    """The adoption rule for `change` against `base` over the pseudo-tasks (one test board's results) where both
    were scored."""
    paired = {r["task"]["key"]: {row["seed"]: (row["scores"][base], row["scores"][change]) for row in r["rows"]}
              for r in results if base in r["candidates"] and change in r["candidates"]}
    out = keep_change(paired)
    return {"base": base, "change": change, "pseudo_tasks": sorted(paired), **out}
