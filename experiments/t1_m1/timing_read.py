#!/usr/bin/env python
"""Try T67 (experiments/harness/TRIES.md): per-gene timing of the recomputed Feng direction from Feng's post-window
slope, the moonshot's kill test on row T62's five practice tasks.

The recomputed-groups file moves each recomputed Feng group by a fraction of the E9.5 -> E14.5 log chord D for every
gene, as if each gene changed at a steady rate. The arm replaces D, per group, by the cubic Hermite's shape through the
bracket (m0 fixed at D / 5; m1 each gene's per-day slope over Feng's E14.5 to E17.5 cells of the same recomputed group;
limiter, split-half and late-cell fallbacks), norm-matched to D (`vec.t1_timing`, the builder's own function: J1). The
practice arm applies the test's gene profile at the practice fraction (1/5), not the practice span's Hermite. This
script imports recomputed_read.py (T62) and, through it, bets_read.py (T61) and gmom_combo.py (T26), and copies none of
their arithmetic; it imports double_step_read.py (T63) for its source guard only.

Inputs, all target-free: the recomputed statistics of row N54's run (sha256 rc.STATS_SHA256; D, the E14.5 point and the
E9.5 level), the late statistics of feng_relabel_late.py (data/derived/feng_relabel_late/, the frozen labeller on
E15.5-E17.5, its manifest committed at HEAD as its RECORD, its sha256 named in row T67), the frozen settings file
(sha256 rc.YAML_SHA256) and the T1 gene names.

Tasks, draws and seeds: row T62's five (the two E9.5 -> E10.5 cuts, chamber primary, start E8.5, build seed 1,000,000 +
seed; T60's three library tasks, split 20261112, copy_last scored first), seeds 0 to 9. One truth, E10.5: a labelled
in-sample diagnostic.

Arms per seed, from the same build seed and streams:
  m1rc_gmom              gate, not counted: T62's candidate, rebuilt; within GATE_TOL of T62's committed value on
                         every task and seed
  m1rc_timing_gmom       counted: m1rc_gmom with each group's recomputed direction replaced by the norm-matched timing
                         direction (`vec.t1_timing.group_timing`'s "direction")
  m1rc_timing_perm_gmom  counted: the same with the deviations m1 - D / 5 permuted within 20 bins of |D| crossed with
                         Feng E9.5 level quintiles (seed 20261122), then the same limiter and norm match
On the first seed of each task the arm is also built with m1 := D / 5 for every gene and must equal m1rc_gmom bit for
bit (its prediction is not scored).

  pregate  Target-free, not counted: no T1 cell, no truth. On the T1 gene names: each group's ||direction - D|| / ||D||,
           the fallback share, the late cell counts, the steady-rate identity. Idle when every group's ratio is below
           0.05: exit IDLE_EXIT, the row is pruned before any scoring and no try is counted. Writes tasks/pregate.json
           and tasks/pregate.txt; prints one step line.
  run-cut  One E10.5 cut. Refused before anything loads unless row T67 is in TRIES.md, the output file does not exist,
           the tree is clean, the source guards pass, the relabelling's and the late labelling's checks pass and the
           pre-gate of this commit says not idle. Prints no score.
  run-lib  One library task; the same refusals. Prints no score.
  read     T67's statement, applied before anything is printed; writes t67.json and t67.txt; prints the outcome words.

Usage: on the DGX campaign checkout ~/v/vec-campaign-d, through experiments/t1_m1/timing.sh (one step per call).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
from vec import t1_timing as tt
from vec.adoption import SEEDS
from vec.boards import get_board, load_boards
from vec.evaluate import Stage
from vec.io import dense, read_stage
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
TRIES = ROOT / "experiments" / "harness" / "TRIES.md"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments" / "external"))
import double_step_read as ds  # noqa: E402  (T63: its source guard only)
import feng_relabel_late as frl  # noqa: E402  (the late labelling's names and record)
import recomputed_read as rc  # noqa: E402  (T62: its arms' pieces, guards and checks)

br, gc = rc.br, rc.gc
bt, gv, pt, rr = rc.bt, rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs
dl = feng_dirs.dl

ROW = "T67"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"
CANDIDATE = rc.CANDIDATE  # m1rc_gmom: the gate
TIMING = "m1rc_timing_gmom"
PERM = "m1rc_timing_perm_gmom"
STEADY = "m1rc_timing_steady"  # the arm with m1 := D / 5: built on the first seed, never scored
ARMS = (CANDIDATE, TIMING, PERM)
COUNTED = (TIMING, PERM)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first
LIB_ARMS = (LIB_REFERENCE, *ARMS)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 10
GATE_TOL = 0.02
PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2
FRACTION_TOL = 1e-12
PREGATE_FILES = ("pregate.json", "pregate.txt")
OUT_FILES = ("t67.json", "t67.txt")
IDLE_EXIT = 4
CHAMBER_BAR = 0.3
WHOLE_BAR = 0.0
LIB_BAR = 0.3
LIB_FLOOR = -0.3
LIB3 = "3"
PERM_CHAMBER_BAR = 0.2
PERM_LIB_BAR = 0.3
MIN_SEEDS = 8
SHARE_BAR = 0.6
DE_GROUPS = ("de_recovery", "direction")
GROUP_WEIGHTS = rc.GROUP_WEIGHTS
LATE_PATH = ROOT / frl.LATE_DIR / frl.STATS
LATE_MANIFEST = ROOT / frl.LATE_DIR / "MANIFEST.json"
LATE_RECORD = frl.RECORD

LABEL = ("try T67, a labelled in-sample diagnostic, never an adoption result: the per-gene Hermite timing of the "
         "recomputed Feng direction from Feng's post-window slope (m1rc_timing_gmom) against m1rc_gmom (T62's "
         "candidate) and against its timing-permuted control, on row T62's five tasks; every task scores one truth "
         "(E10.5) over one day; the practice arm carries the test's gene profile at the practice size")
AHEAD_WORDS = "AHEAD"
NOT_AHEAD_WORDS = ("NOT AHEAD: this closes 'per-gene Hermite timing from Feng's post-window slope' on T1")


# ---------------------------------------------------------------------------------------------------------------------
# Guards


def source_guard(path: Path | None = None) -> list[str]:
    """T63's guard on this file (recomputed_read's guard, no import of a script forming the atlas direction on another
    route, no dynamic import)."""
    return ds.source_guard(Path(__file__ if path is None else path))


def row_line(row: str = ROW, path: Path | None = None) -> str | None:
    return frl.row_line(row, TRIES if path is None else path)


def late_checks(manifest: Path | None = None, record: Path | None = None, stats_path: Path | None = None,
                line: str | None = None, need_stats: bool = True, record_at_head: bytes | None = None) -> dict:
    """The late labelling's guards: its manifest (written by feng_relabel_late from the frozen settings, the late stages
    only, authors_labels_read false, usable, not stopped, a clean gate of an N row), its copy committed at HEAD as the
    RECORD byte for byte, the statistics file's sha256 (hashed when need_stats) equal to the manifest's, and that sha256
    named in row T67's line. Returns the values and the failed checks."""
    manifest = LATE_MANIFEST if manifest is None else Path(manifest)
    record = LATE_RECORD if record is None else Path(record)
    reasons = []
    raw = manifest.read_bytes() if manifest.exists() else None
    man = json.loads(raw) if raw else {}
    if not man:
        reasons.append(f"no late labelling manifest at {manifest}")
    head = frl.fr.committed_at_head(record) if record_at_head is None else record_at_head
    if raw is None or head != raw:
        reasons.append("the late labelling's manifest is not committed at HEAD as its RECORD, byte for byte")
    if man.get("source") != frl.SOURCE:
        reasons.append(f"the late manifest's source is {man.get('source')!r}, not {frl.SOURCE!r}")
    if (man.get("config") or {}).get("sha256") != rc.YAML_SHA256:
        reasons.append("the late labelling ran under another settings file")
    if man.get("authors_labels_read") is not False:
        reasons.append(f"the late manifest's authors_labels_read is {man.get('authors_labels_read')!r}, not false")
    if man.get("stages") != list(frl.LATE_STAGES):
        reasons.append(f"the late manifest's stages are {man.get('stages')}, not {list(frl.LATE_STAGES)}")
    if man.get("stopped") or man.get("usable") is not True:
        reasons.append("the late labelling stopped or failed the frozen yaml's unassigned cap (usable not true)")
    gate = man.get("gate") or {}
    if gate.get("clean") is not True or not frl.N_ROW.fullmatch(str(gate.get("row"))):
        reasons.append("the late labelling has no clean gate of an N row")
    sha = ((man.get("outputs") or {}).get(frl.STATS) or {}).get("sha256")
    hashed = None
    if need_stats:
        path = LATE_PATH if stats_path is None else Path(stats_path)
        hashed = rc.sha256_file(path) if path.exists() else None
        if sha is None or hashed != sha:
            reasons.append(f"the late statistics file's sha256 is {hashed}, not the manifest's {sha}")
    line = row_line() if line is None else line
    if not sha or line is None or sha not in line:
        reasons.append(f"row {ROW}'s line does not name the late statistics' sha256 {sha}")
    return {"stats_sha256": sha, "stats_sha256_hashed": hashed, "yaml_sha256": (man.get("config") or {}).get("sha256"),
            "authors_labels_read": man.get("authors_labels_read"), "late_row": gate.get("row"),
            "usable": man.get("usable"), "passes": not reasons, "reasons": reasons}


def refuse(out: Path, names: Sequence[str], need_pregate: bool = True) -> tuple[dict, dict, dict]:
    """The refusals before anything loads: the row, an existing output file, an unclean tree, the source guards, the
    relabelling's and the late labelling's checks, and (for a run) the pre-gate of this commit saying not idle."""
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, tuple(names))
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    checks = rc.relabel_checks()
    if not checks["passes"]:
        sys.exit(f"{ROW} is refused: {checks['reasons']}")
    late = late_checks()
    if not late["passes"]:
        sys.exit(f"{ROW} is refused: {late['reasons']}")
    if not frl.x002_names_this_use():
        sys.exit(f"{ROW} is refused: EXTERNAL_DATA.md's X002 row does not name {frl.X002_TOKEN!r}")
    if need_pregate:
        why = pregate_refusal(out / PREGATE_FILES[0], started["commit"])
        if why:
            sys.exit(f"{ROW} is refused: {why}")
    return started, checks, late


def pregate_refusal(path: Path, commit: str) -> str | None:
    """Why a run may not start on the pre-gate at `path` (None when it may): it exists, says not idle, holds the
    steady-rate identity and was written on a clean tree at `commit`."""
    if not Path(path).exists():
        return "the pre-gate has not run"
    p = json.loads(Path(path).read_text())
    if p.get("idle") is not False:
        return "the pre-gate found the arm idle (pruned before any scoring)"
    if p.get("steady_identity") is not True:
        return "the pre-gate did not hold the steady-rate identity"
    st = p.get("stamp") or {}
    if st.get("commit") != commit or st.get("tree_clean") is not True:
        return "the pre-gate was not written on a clean tree at this commit"
    return None


# ---------------------------------------------------------------------------------------------------------------------
# The weights (target-free)


def _levels(st: Mapping, stage: str, labels, idx: np.ndarray, eps: float, halves) -> tuple[np.ndarray, int]:
    sel = dl.feng_sel(st, stage, labels, halves=halves)
    if sel is None:
        return np.full(len(idx), np.nan), 0
    return dl.on_t1(tt.log_level(sel["cp"], eps), idx), int(sel["n"])


def timing_kit(genes: Sequence[str], directions: Mapping[str, np.ndarray], fallbacks: Mapping, st_e: Mapping,
               st_l: Mapping, cfg: Mapping) -> dict:
    """Per group with a recomputed direction: the late slope over SLOPE_STAGES (E14.5 from the recomputed statistics
    `st_e`, the later stages from the late statistics `st_l`) from both halves pooled and from each half alone, on the
    labels of the group the direction came from (the yaml's fallback, `fallbacks[g]['takes']`), the group's Feng E9.5
    mean CP10k, the late cell counts, and `vec.t1_timing.group_timing`'s result."""
    eps = float(cfg["metrics"]["pseudocount_cp10k"])
    idx = {"e": dl.sym_index(st_e["genes"], list(genes)), "l": dl.sym_index(st_l["genes"], list(genes))}
    out = {}
    for g in sorted(directions):
        if g not in tt.GROUP_ORDER:
            raise ValueError(f"{g} has a recomputed direction but is not one of {tt.GROUP_ORDER} (R2's guard)")
        src = (fallbacks.get(g) or {}).get("takes") or g
        labels = (cfg["external_labels"].get(src) or [None])[0]
        if not labels:
            raise ValueError(f"{g}: its direction's group {src} has no recomputed labels")
        slopes, cells = {}, {}
        for key, halves in (("AB", ("A", "B")), ("A", ("A",)), ("B", ("B",))):
            rows = []
            for stage in tt.SLOPE_STAGES:
                which = "e" if stage == tt.SLOPE_STAGES[0] else "l"
                lev, n = _levels(st_e if which == "e" else st_l, stage, labels, idx[which], eps, halves)
                rows.append(lev)
                if key == "AB":
                    cells[stage] = n
            slopes[key] = tt.ols_slope(np.vstack(rows))
        e95 = dl.feng_sel(st_e, feng_dirs.BRACKET[0], labels)
        level = dl.on_t1(e95["cp"], idx["e"]) if e95 is not None else np.full(len(genes), np.nan)
        res = tt.group_timing(g, directions[g], slopes["AB"], slopes["A"], slopes["B"], level,
                              {s: cells[s] for s in tt.LATE_STAGES})
        res.update({"source_group": src, "labels": list(labels), "cells": cells, "m1_raw": slopes["AB"]})
        out[g] = res
    return out


def arm_directions(directions: Mapping[str, np.ndarray], kit: Mapping) -> dict[str, dict[str, np.ndarray]]:
    """arm -> group -> the direction the build applies: m1rc_gmom's D, the timing arm's, the control's, and the
    steady-rate arm's (m1 := D / 5)."""
    return {CANDIDATE: dict(directions), TIMING: {g: kit[g]["direction"] for g in directions},
            PERM: {g: kit[g]["perm_direction"] for g in directions},
            STEADY: {g: tt.timing_direction(d, tt.steady_m1(d)) for g, d in directions.items()}}


def direction_sha(dirs: Mapping[str, np.ndarray]) -> dict[str, str]:
    return {g: hashlib.sha256(np.ascontiguousarray(d, dtype=np.float64).tobytes()).hexdigest()
            for g, d in sorted(dirs.items())}


def steady_holds(directions: Mapping[str, np.ndarray], steady: Mapping[str, np.ndarray]) -> bool:
    return set(directions) == set(steady) and direction_sha(directions) == direction_sha(steady)


def kit_record(kit: Mapping, dirs: Mapping[str, Mapping]) -> dict:
    """What a task file and the pre-gate keep: per group the ratios, the fallbacks, the cells, and each arm's
    direction sha256."""
    groups = {g: {"idle_ratio": r["idle_ratio"], "perm_idle_ratio": r["perm_idle_ratio"],
                  "fallback_share": r["fallback_share"], "group_ok": r["group_ok"], "genes": r["genes"],
                  "cells": r["cells"], "source_group": r["source_group"]} for g, r in kit.items()}
    return {"groups": groups, "idle": tt.is_idle({g: r["idle_ratio"] for g, r in kit.items()}),
            "steady_identity": steady_holds(dirs[CANDIDATE], dirs[STEADY]),
            "direction_sha256": {a: direction_sha(d) for a, d in dirs.items()}}


def load_kit(genes: Sequence[str]) -> tuple[dict, dict, dict, dict]:
    """From the files: the recomputed directions on `genes` (rc.recomputed_inputs, feng_dirs' checks), the recomputed
    and late statistics (each through its checks), and the timing kit. Returns (arm directions, kit, its record,
    the recomputed record)."""
    directions, ext = rc.recomputed_inputs(list(genes))
    path = feng_dirs.recomputed_stats_path()
    feng_dirs.check_recomputed(path)
    st_e = dl.load_stats(path)
    feng_dirs.check_recomputed_keys(st_e)
    frl.fio.verify_manifest_output(LATE_PATH)
    st_l = dl.load_stats(LATE_PATH)
    feng_dirs.check_recomputed_keys(st_l)
    stray = sorted({k.split("|")[0] for k in st_l["keys"]} - set(frl.LATE_STAGES))
    if stray:
        raise ValueError(f"the late statistics hold stages {stray} outside {frl.LATE_STAGES}")
    cfg = feng_dirs.groups_config(feng_dirs.load_config(), feng_dirs.load_relabel_config())
    kit = timing_kit(genes, directions, ext["fallbacks"], st_e, st_l, cfg)
    dirs = arm_directions(directions, kit)
    return dirs, kit, kit_record(kit, dirs), ext


def _late_fields(late: Mapping) -> dict:
    return {k: late[k] for k in ("stats_sha256", "stats_sha256_hashed", "yaml_sha256", "authors_labels_read",
                                 "late_row", "usable")}


def pregate_text(res: Mapping) -> str:
    lines = [f"{ROW}'s target-free pre-gate (no T1 cell, no truth; it counts no try): idle "
             f"{'yes: the row is pruned before any scoring' if res['idle'] else 'no'}; steady-rate identity "
             f"{res['steady_identity']}"]
    for g, r in res["groups"].items():
        lines.append(f"{g} (labels of {r['source_group']}): ||w_scaled - 0.4 D|| / ||0.4 D|| {r['idle_ratio']:.4f} "
                     f"(control {r['perm_idle_ratio']:.4f}); fallback share {r['fallback_share']:.3f}; late-cell "
                     f"minimum met {r['group_ok']}; cells per stage {r['cells']}")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def cmd_pregate(a) -> int:
    out = RESULTS / a.out / "tasks"
    started, checks, late = refuse(out, PREGATE_FILES, need_pregate=False)
    genes = list(get_board(rr.BOARD).genes)  # the gene names alone
    _, _, record, ext = load_kit(genes)
    res = {"label": f"{ROW}'s pre-gate", "row": ROW, "bar": tt.IDLE_BAR, **record, "relabel":
           rc._relabel_fields(ext, checks), "late": _late_fields(late),
           "stamp": gv.stamp(started, [feng_dirs.recomputed_stats_path(), LATE_PATH, LATE_MANIFEST,
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, LATE_RECORD])}
    out.mkdir(parents=True, exist_ok=True)
    (out / PREGATE_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / PREGATE_FILES[1]).write_text(pregate_text(res))
    if res["idle"]:
        print("pre-gate: the arm is IDLE; the row is pruned before any scoring and no try is counted")
        return IDLE_EXIT
    print("pre-gate: not idle; the runs may start")
    return 0


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def prediction_sha(stage) -> str:
    return ds.prediction_sha(stage)


def cut_arms(ctx: harness.BuildContext, factors, dirs: Mapping[str, Mapping], fraction: float,
             held_m1: Sequence[str], steady: bool = False) -> dict[str, Stage]:
    """Every arm of one seed on a harness cut, each T62's candidate build (`t1_scope.build_submission` with the
    candidate's g_mom hook and mask) on its own directions at `fraction`, finished as the harness finishes a candidate;
    with `steady`, also the steady-rate arm (STEADY, not scored)."""
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, dirs[CANDIDATE])
    hook = rc.rc_hook(factors, ctx.seed, held, with_direction, covered)

    def build(d) -> Stage:
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src, post_cells=hook,
                                           held_flat=list(held_m1), directions=d, direction_fraction=fraction,
                                           direction_groups=groups)
        return harness.finish(ctx, ref, dense(sub.X), None)

    names = (*ARMS, STEADY) if steady else ARMS
    return {a: build(dirs[a]) for a in names}


def library_arms(e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], dirs: Mapping[str, Mapping], fraction: float,
                 n_cells: int = rr.N_CELLS, steady: bool = False) -> dict:
    """Every arm of one seed on a library task, from T60's oracle draw (recomputed_read.library_pair's route:
    `gmom_combo.moved_by_feng` on stream [seed, 2], then the candidate's g_mom), as AnnData."""
    base, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    labels = e95.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, dirs[CANDIDATE])
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    names = (*ARMS, STEADY) if steady else ARMS
    return {a: hook(gc.moved_by_feng(base, dirs[a], fraction, groups, list(held_m1), [seed, 2])) for a in names}


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def _check_fraction(fraction: float, where: str) -> None:
    if abs(fraction - PRACTICE_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")


def _task_kit(genes: Sequence[str]) -> tuple[dict, dict, dict]:
    dirs, _, record, ext = load_kit(genes)
    if not record["steady_identity"]:
        sys.exit(f"{ROW}: with m1 := D / 5 the timing direction is not D bit for bit; nothing is run")
    return dirs, record, ext


def cmd_run_cut(a) -> None:
    out = RESULTS / a.out / "tasks"
    started, checks, late = refuse(out, (CUT_FILES[a.cut],))
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    t = {x.key: x for x in harness.load_tasks()["T1:test"]}[CUTS[a.cut]]
    boards = load_boards()
    genes = boards[t.panel_board].genes
    target, _ = (harness.filter_target(t, data, genes) if t.target_filter is not None else
                 (Stage.from_anndata(read_stage(data / t.target, genes), genes), None))
    ref_stage = Stage.from_anndata(read_stage(data / t.reference, genes), genes)
    board = harness.task_board(t, boards, target.n, ref_stage.n)
    inputs = {r: load_stage(data / f, lineage=False) for r, f in t.inputs.items()}
    preregistered = t1_scope.held_flat_genes()
    view, cboard = harness.TaskView.of(t), harness._candidate_board(board)
    ref_role = harness.base_role(t)
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    _check_fraction(fraction, t.key)
    held_m1, union = rc._held_m1(preregistered)
    dirs, record, ext = _task_kit(list(inputs[ref_role].var_names))
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for i, s in enumerate(seeds):
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms = cut_arms(ctx, factors, dirs, fraction, held_m1, steady=i == 0)
        shas = {k: prediction_sha(v) for k, v in arms.items()}
        scored = {k: v for k, v in arms.items() if k != STEADY}
        row = {"seed": s, "floor": float(f.score), "arms": gv.score_arms(board, target, ref_stage, s, f.anchors,
                                                                         scored),
               "sha256": {k: v for k, v in shas.items() if k != STEADY}, "seconds": round(time.time() - t0, 1)}
        if i == 0:
            row["steady_equal"] = shas[STEADY] == shas[CANDIDATE]
        rows.append(row)
        del arms, scored
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "relabel": rc._relabel_fields(ext, checks), "late": _late_fields(late), "timing": record,
           "constants": {"edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1)},
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       LATE_PATH, LATE_MANIFEST, feng_dirs.RELABEL_CONFIG, rc.RECORD, LATE_RECORD,
                                       ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / CUT_FILES[a.cut]).write_text(json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_lib(a) -> None:
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started, checks, late = refuse(out, (name,))
    seeds = gv.parse_seeds(a.seeds)
    data = Path(a.data)
    board = get_board(rr.BOARD)
    spec = t1_scope.load_spec()
    chamber = spec.in_scope("chamber")
    stage = {k: load_stage(data / f, lineage=False) for k, f in rr.FILES.items()}
    for k, s in stage.items():
        if list(s.var_names) != board.genes:
            sys.exit(f"{k}: genes differ from {rr.BOARD}'s gene order")
        t1_scope.classify(s, source=f"{k}_RNA", spec=spec)
    e85, e95, e105 = stage["E8.5"], stage["E9.5"], stage["E10.5"]
    within = {"surrounding": list(rr.CONTAMINATION_LABELS)}
    fb = {k: spec.fallbacks(f"{k}_RNA") for k in rr.FILES}
    held = t1_scope.held_flat_genes()
    ref_stage = pt.to_stage(e95)
    half_a, _ = rr.split_library(e105.obs_names, a.library)
    truth = e105[rr.chamber_rows(e105, half_a, chamber)].copy()
    shares = t1_scope.load_shares(t1_scope.composition(truth), spec)
    truth_stage = pt.to_stage(truth)
    factors = br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    fraction = feng_dirs.bracket_fraction(9.5, 10.5)
    _check_fraction(fraction, f"library {a.library}")
    held_m1, union = rc._held_m1(held)
    dirs, record, ext = _task_kit(list(e95.var_names))
    rows = []
    for i, s in enumerate(seeds):
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms = library_arms(e95, shares, factors, held, held_m1, s, within, fb, dirs, fraction, steady=i == 0)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        shas = {k: prediction_sha(v) for k, v in stages.items()}
        scored = {k: v for k, v in stages.items() if k != STEADY}
        row = {"seed": s, **rr.score_arms(board, scored, truth_stage, ref_stage, s),
               "sha256": {k: v for k, v in shas.items() if k != STEADY}, "seconds": round(time.time() - t0, 1)}
        if i == 0:
            row["steady_equal"] = shas[STEADY] == shas[CANDIDATE]
        rows.append(row)
        del arms, stages, scored
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows, "fraction": fraction,
           "relabel": rc._relabel_fields(ext, checks), "late": _late_fields(late), "timing": record,
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD,
                         "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES), "held_m1_genes": len(held_m1)},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], LATE_PATH, LATE_MANIFEST,
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, LATE_RECORD, ROOT / union["file"]])}
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def task_files() -> dict[str, str]:
    return rc.task_files()


def committed_t62(name: str, tasks: Path | None = None) -> dict[int, float]:
    """seed -> m1rc_gmom's score in row T62's committed task file."""
    rows = json.loads(((T62_TASKS if tasks is None else tasks) / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][CANDIDATE]["score"]) for r in rows}


def guards(loaded: Mapping[str, Mapping | None], committed: Mapping[str, Mapping[int, float]],
           pregate: Mapping | None, checks: Mapping, late: Mapping, source_bad: Sequence[str],
           seeds: Sequence[int] = SEEDS) -> list[str]:
    """Every void condition of row T67; an empty list means none fired."""
    reasons = []
    commits = set()
    if pregate is None:
        reasons.append("no pre-gate file")
    else:
        commits.add((pregate.get("stamp") or {}).get("commit"))
        if (pregate.get("stamp") or {}).get("tree_clean") is not True:
            reasons.append("the pre-gate's tree was not clean")
        if pregate.get("idle") is not False:
            reasons.append("the pre-gate did not say not idle")
    for k, r in loaded.items():
        if r is None:
            reasons.append(f"{k}: no task file")
            continue
        rows = r.get("rows") or []
        if sorted(int(x["seed"]) for x in rows) != list(seeds):
            reasons.append(f"{k}: seeds {[x['seed'] for x in rows]}, not {list(seeds)}")
        want = LIB_ARMS if k in LIBRARIES else ARMS
        absent = sorted({a for x in rows for a in want if a not in x["arms"] or a not in (x.get("sha256") or {})})
        if absent or not rows:
            reasons.append(f"{k}: arms {absent} missing")
            continue
        stamp = r.get("stamp") or {}
        commits.add(stamp.get("commit"))
        if stamp.get("tree_clean") is not True:
            reasons.append(f"{k}: its tree was not clean at one commit")
        ref = committed.get(k, {})
        worst = max((abs(br.score_of(x, CANDIDATE) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                     for x in rows), default=float("inf"))
        if not worst <= GATE_TOL:
            reasons.append(f"{k}: the gate fails: m1rc_gmom differs from row T62's committed m1rc_gmom by "
                           f"{worst:.3g} (tolerance {GATE_TOL:g})")
        if not abs(float(r.get("fraction", float("nan"))) - PRACTICE_FRACTION) <= FRACTION_TOL:
            reasons.append(f"{k}: the practice fraction is {r.get('fraction')}, not {PRACTICE_FRACTION}")
        rel = r.get("relabel") or {}
        if rel.get("sha256") != rc.STATS_SHA256 or rel.get("stats_sha256_hashed") != rc.STATS_SHA256:
            reasons.append(f"{k}: ran on another recomputed statistics file")
        if rel.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: ran under another settings file")
        if rel.get("route") != rc.GROUPS_OPTION or rel.get("record_route") != rc.GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        if rel.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the relabelling record did not say authors_labels_read false")
        lt = r.get("late") or {}
        want_sha = late.get("stats_sha256")
        if not want_sha or lt.get("stats_sha256") != want_sha or lt.get("stats_sha256_hashed") != want_sha:
            reasons.append(f"{k}: ran on late statistics other than the row's")
        if lt.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: the late labelling ran under another settings file")
        if lt.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the late record did not say authors_labels_read false")
        tm = r.get("timing") or {}
        if tm.get("steady_identity") is not True or rows[0].get("steady_equal") is not True:
            reasons.append(f"{k}: with m1 := D / 5 the arm is not m1rc_gmom bit for bit")
    if len(commits) > 1:
        reasons.append(f"the task files and the pre-gate ran at more than one commit: {sorted(map(str, commits))}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"at the read: {x}" for x in late.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def library_mean(tasks: Mapping[str, Mapping]) -> float:
    return float(np.mean([tasks[lib]["mean"] for lib in LIBRARIES]))


def clause_i(vs: Mapping[str, Mapping]) -> dict:
    """(i) The arm minus m1rc_gmom: a chamber mean at least +0.3 with at least 8 of 10 seeds ahead; a whole-cut mean
    above 0; the libraries' mean at least +0.3 with library 3 ahead on at least 8 of 10; no library mean below -0.3."""
    c, w = vs["chamber"], vs["whole"]
    parts = {"chamber": bool(c["mean"] >= CHAMBER_BAR and c["ahead"] >= MIN_SEEDS),
             "whole": bool(w["mean"] > WHOLE_BAR),
             "library_mean": bool(library_mean(vs) >= LIB_BAR),
             "library_3_seeds": bool(vs[LIB3]["ahead"] >= MIN_SEEDS),
             "no_library_below_floor": not any(vs[lib]["mean"] < LIB_FLOOR for lib in LIBRARIES)}
    return {"holds": all(parts.values()), "parts": parts, "library_mean": library_mean(vs)}


def clause_ii(vp: Mapping[str, Mapping]) -> dict:
    """(ii) The arm minus the permuted control: a chamber mean at least +0.2 with at least 8 of 10 seeds ahead, and the
    libraries' mean at least +0.3."""
    c = vp["chamber"]
    parts = {"chamber": bool(c["mean"] >= PERM_CHAMBER_BAR and c["ahead"] >= MIN_SEEDS),
             "library_mean": bool(library_mean(vp) >= PERM_LIB_BAR)}
    return {"holds": all(parts.values()), "parts": parts, "library_mean": library_mean(vp)}


def clause_iii(vs: Mapping[str, Mapping]) -> dict:
    """(iii) Over the five tasks, at least 60 percent of (i)'s gain (the arm minus m1rc_gmom) lies in de_recovery plus
    direction: the five tasks' summed points in those two groups over their summed paired means; it cannot hold on a
    gain that is not positive."""
    total = float(sum(p["mean"] for p in vs.values()))
    de = float(sum(p["groups"][g] for p in vs.values() for g in DE_GROUPS))
    share = de / total if total > 0 else float("nan")
    return {"holds": bool(total > 0 and share >= SHARE_BAR), "share": share, "gain": total, "de_points": de}


def t67_reading(loaded: Mapping[str, Mapping | None], committed: Mapping, pregate: Mapping | None, checks: Mapping,
                late: Mapping, source_bad: Sequence[str] = ()) -> dict:
    """Row T67's statement, applied to the task files."""
    reasons = guards(loaded, committed, pregate, checks, late, source_bad)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons}
    rows = {k: r["rows"] for k, r in loaded.items()}
    vs = {k: br.paired(v, TIMING, CANDIDATE) for k, v in rows.items()}
    vp = {k: br.paired(v, TIMING, PERM) for k, v in rows.items()}
    i, ii, iii = clause_i(vs), clause_ii(vp), clause_iii(vs)
    ahead = bool(i["holds"] and ii["holds"] and iii["holds"])
    reported = {"perm_minus_m1rc": {k: br.paired(v, PERM, CANDIDATE) for k, v in rows.items()},
                "mean_scores": {k: {a: float(np.mean([br.score_of(x, a) for x in v]))
                                    for a in (LIB_ARMS if k in LIBRARIES else ARMS)} for k, v in rows.items()},
                "inverted_anchor_libraries": [lib for lib in LIBRARIES if not rc.library_valid(loaded[lib])],
                "pregate": {g: {k: x[k] for k in ("idle_ratio", "perm_idle_ratio", "fallback_share", "group_ok")}
                            for g, x in (pregate.get("groups") or {}).items()}}
    return {"outcome": AHEAD_WORDS if ahead else NOT_AHEAD_WORDS, "ahead": ahead, "i": i, "ii": ii, "iii": iii,
            "vs_m1rc": vs, "vs_perm": vp, "reported": reported, "j6_range": sorted([0.0, vs["chamber"]["mean"]])}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + "); no rerun under row T67 and its tries count"]
    lines = [f"Outcome: {r['outcome']}",
             f"(i) against m1rc_gmom: {'holds' if r['i']['holds'] else 'fails'}; (ii) against the permuted control: "
             f"{'holds' if r['ii']['holds'] else 'fails'}; (iii) the gain's share in de_recovery plus direction: "
             f"{'holds' if r['iii']['holds'] else 'fails'}"]
    lines.append("What follows: one priced question to the captain (this file beside R13's)" if r["ahead"] else
                 "What follows: nothing is built; the line is closed on T1")
    return lines


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} of {p['of']}, per seed "
            f"{p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}")


def t67_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        for title, per in (("m1rc_timing_gmom - m1rc_gmom", r["vs_m1rc"]),
                           ("m1rc_timing_gmom - m1rc_timing_perm_gmom", r["vs_perm"]),
                           ("reported: m1rc_timing_perm_gmom - m1rc_gmom", r["reported"]["perm_minus_m1rc"])):
            lines.append(title)
            lines += [_line(k, p) for k, p in per.items()]
        lines.append(f"(i) parts {r['i']['parts']}, library mean {r['i']['library_mean']:+.2f}; (ii) parts "
                     f"{r['ii']['parts']}, library mean {r['ii']['library_mean']:+.2f}; (iii) share "
                     f"{r['iii']['share']:.3f} of a summed gain {r['iii']['gain']:+.2f}")
        lines.append("Mean scores: " + "; ".join(f"{k}: " + ", ".join(f"{a} {v:.2f}" for a, v in m.items())
                                                 for k, m in r["reported"]["mean_scores"].items()))
        lines.append(f"Libraries with inverted anchors (reported only): "
                     f"{', '.join(r['reported']['inverted_anchor_libraries']) or 'none'}")
        lines.append(f"Pre-gate (target-free): {r['reported']['pregate']}")
        lines.append(f"J6 range (0 to the chamber cut's gain; the test's gap is twice every task's, and the practice "
                     f"arm carries the test's gene profile, so mid-window curvature is unpriced): "
                     f"{r['j6_range'][0]:+.2f} to {r['j6_range'][1]:+.2f}")
    lines.append(f"Tries counted: {res['tries']}")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path) -> dict:
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md")
    gv.refuse_existing(out, OUT_FILES)
    started = gv.tree_state()
    paths = {k: out / "tasks" / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    pg = out / "tasks" / PREGATE_FILES[0]
    pregate = json.loads(pg.read_text()) if pg.exists() else None
    committed = {k: committed_t62(f) for k, f in task_files().items()}
    checks = rc.relabel_checks(need_stats=False)  # the stats files' hashes are in each task file
    late = late_checks(need_stats=False)
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    reading = t67_reading(loaded, committed, pregate, checks, late, bad)
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT}
    res["stamp"] = gv.stamp(started, [p for p in [*paths.values(), pg] if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    (out / OUT_FILES[0]).write_text(json.dumps(res, indent=1, default=float) + "\n")
    (out / OUT_FILES[1]).write_text(t67_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pregate", help="try T67's target-free idle pre-gate (counts no try)")
    p.add_argument("--out", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T67")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T67")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T67's reading")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out)
        return 0
    if a.cmd == "pregate":
        return cmd_pregate(a)
    {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
