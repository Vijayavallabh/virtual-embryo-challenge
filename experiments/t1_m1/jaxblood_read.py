#!/usr/bin/env python
"""Try T66 (experiments/harness/TRIES.md): the kill test of the moonshot "held blood moved along a same-modality nuclei
atlas", read on row T62's five practice tasks. Ruling R1 keeps the drawn Blood cells at their counts as drawn; this
moves them by the per-gene change of the erythroid group of JAX/Qiu 2024 (register row X003) recomputed from raw counts
by experiments/external/jax_blood.py (row N54's marker-score machinery, UMI-downsampled nuclei, equal weight per embryo;
early end embryos 1-12, late end embryos 59-63).

The practice analogue of R13's file is row T62's `m1rc_gmom`, taken unchanged from `recomputed_read.cut_pair` and
`recomputed_read.library_pair` (their CANDIDATE): the gate and the base, not counted. The two counted arms are that
base with one change on the same draw, before g_mom (`recomputed_read.rc_hook`, the same factors, stream [seed, 5], the
same mask, which holds Blood, so g_mom draws nothing there):
  m1rc_gmom_jaxblood  the drawn Blood cells (on practice R1's held set is Blood alone) moved by
                      log f = FRACTION x direction, FRACTION = 1/7.5, |log f| clipped at log 3, f = 1 on R2's held genes
                      (the direction step's list, `recomputed_read._held_m1`) and on genes with no JAX value,
                      keep_depth,
                      through `t1_scope.external_shift` on stream [seed, JAX_STREAM] (`move_blood`)
  m1rc_gmom_jaxperm   the same, with the direction's gene-to-value map permuted within PERM_BINS equal-count bins of the
                      early end's erythroid mean level (held and unmapped genes left out of the permutation), one
                      permutation fixed by default_rng(PERM_SEED) (`permuted`): the same values, the same per-cell size
Every cell not labelled Blood keeps the base's bytes in both arms; each run checks that per seed before scoring
(`heldmom_read.held_record` for R1's held set, `blood_record` for Blood) and stops otherwise. The comparator, not
counted and not rerun: row T64's committed `m1rc_gmom_heldmom` per seed (experiments/t1_m1/results/t64_heldmom/tasks).
It imports recomputed_read.py (T62) and heldmom_read.py (T64) and copies none of their arithmetic.

  pregate  Target-free: no T1 cell is loaded or scored and no truth is read. From the reducer's output (its manifest,
           checked) and the T1 gene names alone: the nuclei assigned to the recomputed erythroid group at each end and
           per embryo (the specification's count, before the UMI floor); the share of the direction's squared norm on
           R2's held genes. Passes when each end has at least MIN_END_NUCLEI assigned nuclei and at least MIN_EMBRYOS
           embryos with at least MIN_EMBRYO_NUCLEI assigned nuclei each, the share is below MAX_HELD_SHARE, and the
           reducer's record has none of manifest_problems (a direction among them). Reports, deciding nothing: the
           nuclei entering the direction (at least the UMI target), the per-gene log f summary and each end's
           experimental_batch and day labels. Writes tasks/pregate.json and pregate.txt; prints the outcome words only;
           exits 0 when it passes and PREGATE_FAIL when it fails (then nothing is scored and no try is counted).
  run-cut  One E10.5 cut. Refused before anything loads unless row T66 is in TRIES.md, HEAD is --commit, the tree is
           clean by vec.provenance's rule, the task file does not exist, the source guards pass, the relabelling's
           committed checks pass, and the pre-gate passed at this commit on the same reducer manifest. Prints no score.
  run-lib  One library task; the same refusals. Prints no score.
  read     Row T66's statement, applied before anything is printed; writes t66.json and t66.txt; prints the outcome
           words and what follows only. The same refusals on the row, the commit and the tree.

Usage: through experiments/t1_m1/jaxblood.sh on the DGX campaign checkout ~/v/vec-campaign-d (it alone holds
data/derived/feng_relabel_c/), at the row's commit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from vec import harness, t1_scope
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
import heldmom_read as hmr  # noqa: E402  (T64: the held-cell record and byte check, the source guard)
import jax_blood as jb  # noqa: E402  (the reducer: its constants, whitelist and embryo lists)

rc = hmr.rc  # recomputed_read (T62): its candidate's arms, guards and checks
br, gc = rc.br, rc.gc
bt, gv, pt, rr = rc.bt, rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs
fr, fio = jb.fr, jb.fio

ROW = "T66"
OUT = "t66_jaxblood"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"
T64_TASKS = RESULTS / "t64_heldmom" / "tasks"
T64_ARM = hmr.ARM  # m1rc_gmom_heldmom, the comparator
OUT_FILES = ("t66.json", "t66.txt")
PREGATE_FILES = ("pregate.json", "pregate.txt")

BASE = rc.CANDIDATE  # m1rc_gmom: the gate and the base, not counted
ARM = "m1rc_gmom_jaxblood"
PERM = "m1rc_gmom_jaxperm"
ARMS = (BASE, ARM, PERM)
COUNTED = (ARM, PERM)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first (it sets the anchors)
LIB_ARMS = (LIB_REFERENCE, *ARMS)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 10

BLOOD = "Blood"  # the label the arms move
FRACTION = 1.0 / 7.5  # one day over (mean late day 16.0 - 8.5)
CLIP = math.log(3.0)
JAX_STREAM = 1  # the move's stream [seed, 1] (T64's held-cell stream; the base draws nothing on it)
PERM_SEED = 20261122
PERM_BINS = 10
FENG_FRACTION = rc.PRACTICE_FRACTION  # 0.2: the base's recomputed direction step, T62's
REPRO_TOL = br.REPRO_TOL  # 0.02 board points: m1rc_gmom against row T62's committed m1rc_gmom, per task and seed
FRACTION_TOL = 1e-12
# The pre-gate (fixed in the row)
MIN_END_NUCLEI = 300
MIN_EMBRYO_NUCLEI = 30
MIN_EMBRYOS = 3
MAX_HELD_SHARE = 0.5
PREGATE_FAIL = 3
# The statement (fixed in the row)
LIB_BAR = 0.3  # (i) and (ii): library-task mean
MIN_AHEAD = 8  # seeds of 10
AHEAD_LIBS_ARM = ("1", "3")  # (i)
AHEAD_LIBS_PERM = ("3",)  # (ii)
MIN_T64_LIBS = 2  # (iii)
CUT_FLOOR = -0.05  # (iv)
GROUP_WEIGHTS = rc.GROUP_WEIGHTS
UP = math.log(1.1)

AHEAD_WORDS = {True: "AHEAD", False: "NOT AHEAD"}
FOLLOWS_AHEAD = ("What follows: one priced question to the captain: amend R1 for Blood by this direction, or keep R1; "
                 "its default is keep R1, and only his ruling naming R1 changes a file. Nothing is built or changed by "
                 "this read.")
FOLLOWS_END = ("What follows: NOT AHEAD closes 'JAX erythroid direction on the held Blood cells' on T1. No question "
               "goes to the captain, nothing is built, and nothing is rerun under another row.")
FOLLOWS_VOID = ("What follows: the tries still count; nothing is rerun under row T66, no question goes to the captain "
                "and nothing is built.")
FOLLOWS_PRUNED = ("What follows: nothing is scored and no try is counted under row T66; the arm is pruned at the "
                  "pre-gate.")
BESIDE = ("J2: the October readout already compared JAX's blood direction with the E10.5 truth, so this read carries "
          "no forecast weight.",
          "The arm moves the expression of the drawn Blood cells only; R1 is the captain's ruling, and only his ruling "
          "naming R1 puts the arm in a file. Macrophage and Early VIC, held at the test, are unpriced here.")
LABEL = ("try T66, a labelled in-sample diagnostic, never an adoption result: the drawn Blood cells moved along the "
         "recomputed erythroid direction of JAX/Qiu 2024 nuclei (m1rc_gmom_jaxblood) and its permuted control "
         "(m1rc_gmom_jaxperm) against m1rc_gmom (row T62's) on T62's five tasks; every task scores one truth (E10.5) "
         "over one day")


# ---------------------------------------------------------------------------------------------------------------------
# Guards on this source


def source_guard(path: Path | None = None) -> list[str]:
    """heldmom_read's guard on this file (recomputed_read's guard, the banned imports, no dynamic import)."""
    return hmr.source_guard(Path(__file__ if path is None else path))


def head_is(commit: str) -> dict:
    """The start state, refused unless HEAD is `commit` and the tree is clean by vec.provenance's rule."""
    started = gv.tree_state()
    if started["commit"] != commit:
        sys.exit(f"{ROW} is refused: HEAD is {started['commit']}, not the run commit {commit}")
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    return started


# ---------------------------------------------------------------------------------------------------------------------
# The reducer's output


def reducer_paths() -> tuple[Path, Path]:
    cfg = jb.load_config()
    dst = fio.REPO / cfg["output"]["dir"]
    return dst / "MANIFEST.json", dst / cfg["output"]["direction"]


def manifest_problems(man: Mapping) -> list[str]:
    """Every way the reducer's manifest is not the row's: the yaml's sha256, the gate, the columns read (the whitelist
    check), an embryo outside the row's lists, a direction missing."""
    bad = []
    if (man.get("config") or {}).get("sha256") != jb.YAML_SHA256:
        bad.append("the reducer ran under another yaml than the row's")
    gate = man.get("gate") or {}
    if gate.get("row") != ROW or gate.get("clean") is not True:
        bad.append("the reducer's gate is not row T66's on a clean tree")
    cols = list(man.get("columns_read") or [])
    if not cols or any(c not in jb.CELL_WHITELIST or c in jb.REFUSED_COLUMNS or fr.ANNOTATION_PATTERN.search(c)
                       for c in cols):
        bad.append(f"the reducer read columns outside the whitelist: {cols}")
    if man.get("authors_labels_read") is not False or man.get("stats_npz_read") is not False:
        bad.append("the reducer's record does not say that no author label and no stats file was read")
    want = sorted(n for e in jb.ENDS for n in jb.ROW_EMBRYOS[e])
    if sorted(int(x) for x in man.get("embryos_read") or []) != want:
        bad.append("an embryo outside the row's lists entered, or one is missing")
    lists = man.get("embryo_lists") or {}
    if {e: sorted(lists.get(e) or []) for e in jb.ENDS} != {e: list(jb.ROW_EMBRYOS[e]) for e in jb.ENDS}:
        bad.append("the reducer's embryo lists are not the row's")
    used = man.get("embryos_used") or {}
    if any(int(n) not in jb.ROW_EMBRYOS[e] for e in jb.ENDS for n in used.get(e) or []):
        bad.append("an embryo outside the row's lists entered the direction")
    if man.get("has_direction") is not True:
        bad.append("the reducer wrote no direction")
    return bad


def load_reducer(manifest: Path | None = None, direction: Path | None = None) -> tuple[dict, dict, str]:
    """The reducer's manifest, its direction arrays (checked against the manifest) and the manifest's sha256."""
    if manifest is None or direction is None:
        manifest, direction = reducer_paths()
    man = json.loads(Path(manifest).read_text())
    fio.verify_manifest_output(direction)
    with np.load(direction, allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    return man, arrays, fio.sha256_file(manifest)


def on_t1(arrays: Mapping, genes: Sequence[str], key: str = "direction") -> np.ndarray:
    """A per-JAX-gene array on T1's gene order through JAX's unique symbols (`feng_relabel.unique_index`, as
    dirlib.jax_index drops a symbol the gene table holds twice); NaN where a T1 gene has no JAX value."""
    uidx = fr.unique_index([str(s) for s in arrays["symbol"]])
    vals = np.asarray(arrays[key], dtype=np.float64)
    return np.array([vals[uidx[g]] if g in uidx else np.nan for g in genes], dtype=np.float64)


def permuted(d: np.ndarray, level: np.ndarray, held: np.ndarray, seed: int = PERM_SEED,
             bins: int = PERM_BINS) -> np.ndarray:
    """The direction with its gene-to-value map permuted within `bins` equal-count bins of `level` (the early end's
    erythroid mean CP10k): the genes with a finite value and a finite level, off the held list, ranked by level (ties by
    gene position), cut by np.array_split; one generator, default_rng(seed), permutes each bin in order. Held genes and
    genes with no value are left as they are."""
    d = np.asarray(d, dtype=np.float64)
    level = np.asarray(level, dtype=np.float64)
    ok = np.isfinite(d) & np.isfinite(level) & ~np.asarray(held, dtype=bool)
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(level[idx], kind="stable")]
    rng = np.random.default_rng(seed)
    out = d.copy()
    for b in np.array_split(order, bins):
        out[b] = d[b[rng.permutation(b.size)]]
    return out


def log_factors(d: np.ndarray, held: np.ndarray) -> np.ndarray:
    """The factor external_shift applies (for the summaries): FRACTION x d, NaN as 0, clipped at log 3, 0 on held."""
    logf = np.clip(FRACTION * np.nan_to_num(np.asarray(d, dtype=np.float64), nan=0.0), -CLIP, CLIP)
    logf[np.asarray(held, dtype=bool)] = 0.0
    return logf


def write_new(path: Path, text: str) -> None:
    """Exclusive create (mode 'x'): a second copy of a step that finishes later stops with FileExistsError and never
    replaces the file the first wrote."""
    with open(path, "x") as f:
        f.write(text)


def vector_sha(v: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(v, dtype=np.float64).tobytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate (target-free)


def pregate_reading(man: Mapping, d: np.ndarray, held: np.ndarray, genes: Sequence[str]) -> dict:
    """The pre-gate's fixed checks and its reported tables (deciding nothing), from the manifest and the direction on
    T1's genes."""
    per = man.get("per_embryo") or {}
    ends = {}
    for e in jb.ENDS:
        rows = {int(n): r for n, r in per.items() if int(n) in jb.ROW_EMBRYOS[e]}
        assigned = {n: int(r.get("erythroid_assigned", 0)) for n, r in sorted(rows.items())}
        entering = {n: int(r.get("erythroid_entering", 0)) for n, r in sorted(rows.items())}
        # the specification's count: nuclei assigned to the recomputed erythroid group (before the UMI floor)
        ends[e] = {"assigned": assigned, "entering": entering,
                   "total": int(sum(assigned.values())),
                   "embryos_at_floor": int(sum(v >= MIN_EMBRYO_NUCLEI for v in assigned.values())),
                   "entering_total": int(sum(entering.values())),
                   "experimental_batch": sorted({b for r in rows.values() for b in r.get("experimental_batch", [])}),
                   "days": sorted({x for r in rows.values() for x in r.get("days", [])})}
    fin = np.isfinite(d)
    sq = np.where(fin, d, 0.0) ** 2
    total = float(sq.sum())
    share = float(sq[np.asarray(held, dtype=bool)].sum() / total) if total > 0 else float("nan")
    failed = [f"problem with the reducer's record: {x}" for x in manifest_problems(man)]
    for e in jb.ENDS:
        if ends[e]["total"] < MIN_END_NUCLEI:
            failed.append(f"fewer than {MIN_END_NUCLEI} assigned erythroid nuclei at the {e} end")
        if ends[e]["embryos_at_floor"] < MIN_EMBRYOS:
            failed.append(f"fewer than {MIN_EMBRYOS} embryos with at least {MIN_EMBRYO_NUCLEI} assigned erythroid "
                          f"nuclei at the {e} end")
    if not share < MAX_HELD_SHARE:
        failed.append(f"the share of the direction's squared norm on R2's held genes is not below {MAX_HELD_SHARE:g} "
                      "(the arm is near-idle)")
    logf = log_factors(d, held)
    mapped = fin & ~np.asarray(held, dtype=bool)
    order = np.argsort(logf, kind="stable")
    g = np.asarray(genes, dtype=str)
    summary = {"mapped_genes": int(fin.sum()), "moved_genes": int(mapped.sum()), "held_genes": int(np.sum(held)),
               "quantiles": {str(q): float(np.quantile(logf[mapped], q)) if mapped.any() else None
                             for q in (0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0)},
               "n_up": int((logf > UP).sum()), "n_down": int((logf < -UP).sum()),
               "n_clipped": int((np.abs(FRACTION * np.where(fin, d, 0.0)) > CLIP).sum()),
               "top_up": [[str(g[i]), float(logf[i])] for i in order[::-1][:15]],
               "top_down": [[str(g[i]), float(logf[i])] for i in order[:15]]}
    return {"passes": not failed, "failed": failed, "ends": ends, "held_share": share, "logf_summary": summary}


def pregate_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], pregate_outcome(r)]
    for e, v in r["ends"].items():
        lines.append(f"{e} end: {v['total']} assigned erythroid nuclei, {v['embryos_at_floor']} embryos at the floor; "
                     f"per embryo {v['assigned']}; entering the direction (at least the UMI target) "
                     f"{v['entering_total']}, per embryo {v['entering']}; "
                     f"experimental_batch {v['experimental_batch']}; "
                     f"days {v['days']}")
    s = r["logf_summary"]
    lines.append(f"R2-held share of the squared direction {r['held_share']:.3f}; log f over {s['moved_genes']} moved "
                 f"genes (mapped {s['mapped_genes']}, held {s['held_genes']}): quantiles {s['quantiles']}; up "
                 f"{s['n_up']}, down {s['n_down']}, clipped {s['n_clipped']}")
    lines.append(f"top up {s['top_up']}")
    lines.append(f"top down {s['top_down']}")
    lines.append("The run confound between the early end (run_4) and the late runs is declared, not removed.")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def pregate_outcome(r: Mapping) -> str:
    if r["passes"]:
        return "Pre-gate: passes"
    return "Pre-gate: FAILS (" + "; ".join(r["failed"]) + ")"


def cmd_pregate(a) -> int:
    out = RESULTS / a.out / "tasks"
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, PREGATE_FILES)
    started = head_is(a.commit)
    bad = source_guard()
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    genes = list(get_board(rr.BOARD).genes)  # the gene names alone
    man, arrays, man_sha = load_reducer()
    held_m1, _ = rc._held_m1(t1_scope.held_flat_genes())
    held = np.isin(np.asarray(genes, dtype=str), held_m1)
    reading = pregate_reading(man, on_t1(arrays, genes), held, genes)
    mpath, dpath = reducer_paths()
    res = {"label": LABEL, "row": ROW, "reading": reading, "manifest_sha256": man_sha,
           "yaml_sha256": (man.get("config") or {}).get("sha256"), "fraction": FRACTION, "clip": CLIP,
           "stamp": gv.stamp(started, [mpath, dpath, t1_scope.HELD_FLAT_GENES])}
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / PREGATE_FILES[0], json.dumps(res, indent=1, default=float) + "\n")
    write_new(out / PREGATE_FILES[1], pregate_text(res))
    print("Pre-gate: passes" if reading["passes"] else f"Pre-gate: FAILS\n{FOLLOWS_PRUNED}")  # words only
    return 0 if reading["passes"] else PREGATE_FAIL


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def move_blood(pop, d: np.ndarray, held_m1: Sequence[str], seed) -> tuple:
    """The population with its Blood cells moved by log f = FRACTION x d (clipped at log 3, f = 1 on `held_m1` and where
    d has no value), keep_depth, through `t1_scope.external_shift` on the Blood rows alone (stream `seed`), reassembled
    by `t1_scope._stack_rows`; every other row is the population's own. Returns the population and the shift's report
    (None when no Blood cell was drawn)."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    blood = labels == BLOOD
    if not blood.any():
        return pop, None
    idx, rest = np.flatnonzero(blood), np.flatnonzero(~blood)
    moved, rep = t1_scope.external_shift(pop[idx].copy(), {BLOOD: np.asarray(d, dtype=np.float64)}, FRACTION,
                                         {BLOOD: BLOOD}, list(held_m1), label_key=t1_scope.LABEL_KEY,
                                         max_abs_log_factor=CLIP, keep_depth=True, seed=seed)
    parts = [(idx, moved)] + ([(rest, pop[rest].copy())] if rest.size else [])
    return t1_scope._stack_rows(parts, pop.n_obs), rep


def blood_record(pop, base_X, arm_X, name: str) -> dict:
    """After the byte check: every row whose label (pop's obs, in row order) is not Blood has the base's bytes in the
    arm, or the run stops."""
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    blood = labels == BLOOD
    a = np.ascontiguousarray(dense(base_X), dtype=np.float32)
    b = np.ascontiguousarray(dense(arm_X), dtype=np.float32)
    if a.shape != b.shape or a.shape[0] != len(labels):
        raise SystemExit(f"{ROW}: {BASE} and {name} are not the same cells ({a.shape}, {b.shape}, {len(labels)})")
    equal = a[~blood].tobytes() == b[~blood].tobytes()
    if not equal:
        raise SystemExit(f"{ROW}: a cell not labelled Blood differs between {BASE} and {name}; the run stops")
    return {"blood_cells": int(blood.sum()), "blood_cells_changed": int((a[blood] != b[blood]).any(axis=1).sum()),
            "non_blood_equal": bool(equal)}


def _records(pop, groups, directions, base_X, arms: Mapping) -> dict:
    return {name: {**hmr.held_record(pop, groups, directions, base_X, X, None), **blood_record(pop, base_X, X, name)}
            for name, X in arms.items()}


def cut_arms(ctx: harness.BuildContext, factors, directions: Mapping, fraction: float, held_m1: Sequence[str],
             d_jax: np.ndarray, d_perm: np.ndarray) -> tuple[dict[str, Stage], dict]:
    """m1rc_gmom (recomputed_read.cut_pair's CANDIDATE, as called there) and the two counted arms, each one
    build_submission call with cut_pair's keywords for its candidate and the Blood move before g_mom; each finished as
    the harness finishes a candidate, and the per-arm records."""
    base = rc.cut_pair(ctx, factors, directions, fraction, held_m1)[BASE]
    role = ctx.task.reference_role
    ref, src = harness._shallow(ctx.inputs[role]), ctx.task.role_stems[role]
    spec = t1_scope.load_spec()
    t1_scope.classify(ref, source=src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    n = harness.count_for(ctx.board, ref.n_obs)
    held = list(ctx.derived["t1_held_flat"])
    labels = ref.obs[t1_scope.LABEL_KEY].astype(str).unique()
    groups, covered, with_direction = rc.direction_kit(labels, directions)
    hook = rc.rc_hook(factors, ctx.seed, held, with_direction, covered)
    out, seen, reports = {BASE: base}, {}, {}
    for name, d in ((ARM, d_jax), (PERM, d_perm)):
        def post(pop, name=name, d=d):
            seen[name] = pop  # the drawn and direction-moved population before the Blood move and g_mom
            moved, rep = move_blood(pop, d, held_m1, [ctx.seed, JAX_STREAM])
            reports[name] = None if rep is None else rep.as_dict()
            return hook(moved)

        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src,
                                           held_flat=list(held_m1), directions=directions,
                                           direction_fraction=fraction, direction_groups=groups, post_cells=post)
        out[name] = harness.finish(ctx, ref, dense(sub.X), None)
    if list(seen[ARM].obs_names) != list(seen[PERM].obs_names):
        raise SystemExit(f"{ROW}: the two arms were not built on one draw")
    rec = _records(seen[ARM], groups, directions, base.X, {ARM: out[ARM].X, PERM: out[PERM].X})
    for name in COUNTED:
        rec[name]["shift"] = reports[name]
    return out, rec


def library_population(e95, half_shares: Mapping[str, float], held_m1: Sequence[str], seed: int, within: Mapping,
                       fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                       n_cells: int = rr.N_CELLS) -> tuple:
    """The base's population before g_mom on a library task, as recomputed_read.library_pair forms it: T60's oracle
    draw and the recomputed direction on the covered rows (`gmom_combo.moved_by_feng`, stream [seed, 2]). Returns the
    population, the label -> group map, the covered labels and the groups with a direction."""
    base0, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    groups, covered, with_direction = rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), directions)
    moved = gc.moved_by_feng(base0, directions, fraction, groups, list(held_m1), [seed, 2])
    return moved, groups, covered, with_direction


def library_arms(e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                 d_jax: np.ndarray, d_perm: np.ndarray, n_cells: int = rr.N_CELLS) -> tuple[dict, dict]:
    """m1rc_gmom (recomputed_read.library_pair's CANDIDATE, as called there) and the two counted arms on the same draw,
    as AnnData, and the per-arm records."""
    base = rc.library_pair(e95, half_shares, factors, held, held_m1, seed, within, fb, directions, fraction,
                           n_cells=n_cells)[BASE]
    pop, groups, covered, with_direction = library_population(e95, half_shares, held_m1, seed, within, fb,
                                                              directions, fraction, n_cells=n_cells)
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    out, reports = {BASE: base}, {}
    for name, d in ((ARM, d_jax), (PERM, d_perm)):
        moved, rep = move_blood(pop, d, held_m1, [seed, JAX_STREAM])
        out[name] = hook(moved)
        reports[name] = None if rep is None else rep.as_dict()
    rec = _records(pop, groups, directions, base.X, {ARM: out[ARM].X, PERM: out[PERM].X})
    for name in COUNTED:
        rec[name]["shift"] = reports[name]
    return out, rec


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def refuse(out: Path, name: str, commit: str) -> tuple[dict, dict, dict]:
    """The refusals before anything loads: the row, an existing output file, HEAD and the tree, the source guards, the
    relabelling's committed checks and the pre-gate (passed, at this commit, on the reducer's present manifest).
    Returns the start state, the relabelling's checks and the pre-gate record."""
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, (name,))
    started = head_is(commit)
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    checks = rc.relabel_checks()
    if not checks["passes"]:
        sys.exit(f"{ROW} is refused: {checks['reasons']}")
    pg_path = out / PREGATE_FILES[0]
    pg = json.loads(pg_path.read_text()) if pg_path.exists() else None
    mpath, _ = reducer_paths()
    if pg is None or (pg.get("reading") or {}).get("passes") is not True:
        sys.exit(f"{ROW} is refused: the pre-gate has not passed (nothing is scored after a failed pre-gate)")
    if (pg.get("stamp") or {}).get("commit") != commit or pg.get("manifest_sha256") != fio.sha256_file(mpath):
        sys.exit(f"{ROW} is refused: the pre-gate passed at another commit or on another reducer manifest")
    return started, checks, pg


def jax_inputs(genes: Sequence[str], held_m1: Sequence[str]) -> tuple[np.ndarray, np.ndarray, dict]:
    """The JAX direction and its permuted control on `genes`, and the fields every task file records."""
    man, arrays, man_sha = load_reducer()
    bad = manifest_problems(man)
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    held = np.isin(np.asarray(genes, dtype=str), list(held_m1))
    d = on_t1(arrays, genes)
    d_perm = permuted(d, on_t1(arrays, genes, "m_early"), held)
    _, dpath = reducer_paths()
    fields = {"manifest_sha256": man_sha, "yaml_sha256": (man.get("config") or {}).get("sha256"),
              "direction_sha256": fio.sha256_file(dpath), "embryos_read": man.get("embryos_read"),
              "embryo_lists": man.get("embryo_lists"), "columns_read": man.get("columns_read"),
              "authors_labels_read": man.get("authors_labels_read"), "stats_npz_read": man.get("stats_npz_read"),
              "d_sha256": vector_sha(d), "d_perm_sha256": vector_sha(d_perm),
              "perm_same_values": bool(np.array_equal(np.sort(d[np.isfinite(d)]),
                                                      np.sort(d_perm[np.isfinite(d_perm)])))}
    return d, d_perm, fields


def _constants(held_m1: Sequence[str]) -> dict:
    return {"fraction": FRACTION, "clip": CLIP, "jax_stream": JAX_STREAM, "perm_seed": PERM_SEED,
            "perm_bins": PERM_BINS, "feng_fraction": FENG_FRACTION, "held_m1_genes": len(held_m1),
            "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES)}


def cmd_run_cut(a) -> None:
    out = RESULTS / a.out / "tasks"
    started, checks, pg = refuse(out, CUT_FILES[a.cut], a.commit)
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
    directions, ext = rc.recomputed_inputs(list(inputs[ref_role].var_names))
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    hmr.check_fraction(fraction, t.key)
    held_m1, union = rc._held_m1(preregistered)
    d_jax, d_perm, jfields = jax_inputs(list(inputs[ref_role].var_names), held_m1)
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms, rec = cut_arms(ctx, factors, directions, fraction, held_m1, d_jax, d_perm)
        shas = {k: hmr.prediction_sha(v.X) for k, v in arms.items()}
        rows.append({"seed": s, "floor": float(f.score), "arms": gv.score_arms(board, target, ref_stage, s, f.anchors,
                                                                               arms),
                     "sha256": shas, "records": rec, "seconds": round(time.time() - t0, 1)})
        del arms
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    mpath, dpath = reducer_paths()
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "feng_fraction": fraction,
           "relabel": rc._relabel_fields(ext, checks), "jax": jfields, "pregate_sha256": pg["manifest_sha256"],
           "constants": _constants(held_m1),
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, ROOT / union["file"], mpath, dpath,
                                       jb.CONFIG])}
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / CUT_FILES[a.cut], json.dumps(res, indent=1, default=float) + "\n")


def cmd_run_lib(a) -> None:
    out = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started, checks, pg = refuse(out, name, a.commit)
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
    directions, ext = rc.recomputed_inputs(list(e95.var_names))
    fraction = feng_dirs.bracket_fraction(9.5, 10.5)
    hmr.check_fraction(fraction, f"library {a.library}")
    held_m1, union = rc._held_m1(held)
    d_jax, d_perm, jfields = jax_inputs(list(e95.var_names), held_m1)
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms, rec = library_arms(e95, shares, factors, held, held_m1, s, within, fb, directions, fraction, d_jax,
                                 d_perm)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        shas = {k: hmr.prediction_sha(v.X) for k, v in stages.items()}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s), "sha256": shas,
                     "records": rec, "seconds": round(time.time() - t0, 1)})
        del arms, stages
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    mpath, dpath = reducer_paths()
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows,
           "feng_fraction": fraction, "relabel": rc._relabel_fields(ext, checks), "jax": jfields,
           "pregate_sha256": pg["manifest_sha256"],
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD, **_constants(held_m1)},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], feng_dirs.RELABEL_CONFIG,
                                       rc.RECORD, ROOT / union["file"], mpath, dpath, jb.CONFIG])}
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / name, json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def task_files() -> dict[str, str]:
    return rc.task_files()


def committed(name: str, tasks: Path, arm: str) -> dict[int, float]:
    """seed -> score of `arm` in a committed task file."""
    rows = json.loads((tasks / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][arm]["score"]) for r in rows}


def library_valid(r: Mapping) -> bool:
    return not any(x.get("inverted_anchors") for x in r["rows"])


def guards(loaded: Mapping[str, Mapping | None], t62: Mapping[str, Mapping[int, float]],
           t64: Mapping[str, Mapping[int, float]], checks: Mapping, pregate: Mapping | None, manifest_sha: str | None,
           source_bad: Sequence[str], seeds: Sequence[int] = SEEDS, read_commit: str | None = None) -> list[str]:
    """Every void condition of row T66; an empty list means none fired. `read_commit`, when given (the read passes its
    own HEAD), must be the task files' commit."""
    reasons = []
    if pregate is None or (pregate.get("reading") or {}).get("passes") is not True:
        reasons.append("the pre-gate record is missing or did not pass")
    if pregate is not None and (pregate.get("stamp") or {}).get("tree_clean") is not True:
        reasons.append("the pre-gate's tree was not clean")
    commits = set()
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
            reasons.append(f"{k}: arms or sha256 {absent} missing")
            continue
        stamp = r.get("stamp") or {}
        commits.add(stamp.get("commit"))
        if stamp.get("tree_clean") is not True:
            reasons.append(f"{k}: its tree was not clean at one commit")
        ref = t62.get(k, {})
        worst = max((abs(br.score_of(x, BASE) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                     for x in rows), default=float("inf"))
        if not worst <= REPRO_TOL:
            reasons.append(f"{k}: the gate failed: {BASE} differs from row T62's committed {BASE} by up to "
                           f"{worst:.3g} (tolerance {REPRO_TOL:g})")
        if k in LIBRARIES and sorted(t64.get(k, {})) != list(seeds):
            reasons.append(f"{k}: row T64's committed {T64_ARM} is missing a seed")
        if not abs(float(r.get("feng_fraction", float("nan"))) - FENG_FRACTION) <= FRACTION_TOL:
            reasons.append(f"{k}: the base's fraction is {r.get('feng_fraction')}, not {FENG_FRACTION}")
        c = r.get("constants") or {}
        if (c.get("fraction"), c.get("clip"), c.get("jax_stream"), c.get("perm_seed"), c.get("perm_bins")) != (
                FRACTION, CLIP, JAX_STREAM, PERM_SEED, PERM_BINS):
            reasons.append(f"{k}: the move's constants are not the row's")
        rel = r.get("relabel") or {}
        if rel.get("sha256") != rc.STATS_SHA256 or rel.get("stats_sha256_hashed") != rc.STATS_SHA256:
            reasons.append(f"{k}: ran on another recomputed statistics file")
        if rel.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: ran under another relabelling settings file")
        if rel.get("route") != rc.GROUPS_OPTION or rel.get("record_route") != rc.GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        j = r.get("jax") or {}
        if j.get("yaml_sha256") != jb.YAML_SHA256:
            reasons.append(f"{k}: the reducer's yaml sha256 is not the row's")
        if not manifest_sha or j.get("manifest_sha256") != manifest_sha or r.get("pregate_sha256") != manifest_sha:
            reasons.append(f"{k}: the reducer's manifest sha256 is not the pre-gate's and the present one")
        cols = list(j.get("columns_read") or [])
        if not cols or any(x not in jb.CELL_WHITELIST or x in jb.REFUSED_COLUMNS for x in cols) \
                or j.get("authors_labels_read") is not False or j.get("stats_npz_read") is not False:
            reasons.append(f"{k}: a label column or a label-keyed file was read")
        want_e = sorted(n for e in jb.ENDS for n in jb.ROW_EMBRYOS[e])
        if sorted(int(x) for x in j.get("embryos_read") or []) != want_e:
            reasons.append(f"{k}: an embryo outside the row's lists entered, or one is missing")
        if j.get("perm_same_values") is not True:
            reasons.append(f"{k}: the permuted control does not hold the direction's values")
        for x in rows:
            for name in COUNTED:
                h = (x.get("records") or {}).get(name) or {}
                if h.get("non_blood_equal") is not True or h.get("non_held_equal") is not True:
                    reasons.append(f"{k}: seed {x['seed']}: a non-Blood cell is not byte-equal in {name}")
                if sorted(h.get("held_set") or []) != sorted(h.get("uncovered_set") or []) or "held_set" not in h:
                    reasons.append(f"{k}: seed {x['seed']}: the held set is not build_submission's uncovered set")
                elif list(h.get("held_set") or []) != [BLOOD]:
                    reasons.append(f"{k}: seed {x['seed']}: the held set drawn is {h.get('held_set')}, not Blood alone "
                                   "(the arm moves Blood only, so it would not be the R1-held set)")
    if pregate is not None and commits and (pregate.get("stamp") or {}).get("commit") not in commits:
        reasons.append("the pre-gate ran at another commit than the task files")
    if pregate is not None and pregate.get("manifest_sha256") != manifest_sha:
        reasons.append("the pre-gate ran on another reducer manifest")
    if len(commits) > 1:
        reasons.append(f"the task files ran at more than one commit: {sorted(map(str, commits))}")
    if read_commit is not None and commits and commits != {read_commit}:
        reasons.append(f"the task files ran at another commit than the read's HEAD {read_commit}")
    reasons += [f"at the read: {x}" for x in checks.get("reasons", [])]
    reasons += [f"source: {x}" for x in source_bad]
    return reasons


def paired(rows: Sequence[Mapping], change: str, base: str) -> dict:
    """change - base per seed (`bets_read.paired`) with the seeds below 0."""
    p = br.paired(rows, change, base)
    p["behind"] = int(sum(br.score_of(r, change) - br.score_of(r, base) < 0 for r in rows))
    return p


def statement(arm: Mapping[str, Mapping], ap: Mapping[str, Mapping], arm_means: Mapping[str, float],
              t64_means: Mapping[str, float]) -> tuple[bool, list[str], dict]:
    """Row T66's statement: (i) to (iv), with d = arm - m1rc_gmom (`arm`) and arm - jaxperm (`ap`) paired per seed.
    Returns AHEAD, the clauses that failed (words only) and the values the clauses read."""
    lib_d = float(np.mean([arm[lib]["mean"] for lib in LIBRARIES]))
    lib_ap = float(np.mean([ap[lib]["mean"] for lib in LIBRARIES]))
    at_t64 = [lib for lib in LIBRARIES if arm_means[lib] >= t64_means[lib]]
    failed = []
    if not lib_d >= LIB_BAR:
        failed.append(f"(i) the library-task mean of the arm minus {BASE} is below +{LIB_BAR:g}")
    for lib in AHEAD_LIBS_ARM:
        if not arm[lib]["ahead"] >= MIN_AHEAD:
            failed.append(f"(i) library {lib} is ahead of {BASE} in fewer than {MIN_AHEAD} of 10 seeds")
    if not lib_ap >= LIB_BAR:
        failed.append(f"(ii) the library-task mean of the arm minus {PERM} is below +{LIB_BAR:g}")
    for lib in AHEAD_LIBS_PERM:
        if not ap[lib]["ahead"] >= MIN_AHEAD:
            failed.append(f"(ii) library {lib} is ahead of {PERM} in fewer than {MIN_AHEAD} of 10 seeds")
    if len(at_t64) < MIN_T64_LIBS:
        failed.append(f"(iii) the arm is at or above row T64's {T64_ARM} on fewer than {MIN_T64_LIBS} of 3 libraries")
    for cut in CUTS:
        if not arm[cut]["mean"] >= CUT_FLOOR:
            failed.append(f"(iv) the {cut} cut's mean of the arm minus {BASE} is below {CUT_FLOOR:g}")
    return not failed, failed, {"library_mean_d": lib_d, "library_mean_arm_minus_perm": lib_ap,
                                "libraries_at_or_above_t64": at_t64}


def t66_reading(loaded: Mapping[str, Mapping | None], t62: Mapping[str, Mapping[int, float]],
                t64: Mapping[str, Mapping[int, float]], checks: Mapping, pregate: Mapping | None,
                manifest_sha: str | None, source_bad: Sequence[str] = (), read_commit: str | None = None) -> dict:
    """Row T66's reading: the void conditions (the gate among them) first, then the statement."""
    reasons = guards(loaded, t62, t64, checks, pregate, manifest_sha, source_bad, read_commit=read_commit)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons, "follows": FOLLOWS_VOID}
    arm = {k: paired(r["rows"], ARM, BASE) for k, r in loaded.items()}
    perm = {k: paired(r["rows"], PERM, BASE) for k, r in loaded.items()}
    ap = {k: paired(r["rows"], ARM, PERM) for k, r in loaded.items()}
    means = {k: {a: float(np.mean([br.score_of(x, a) for x in r["rows"]])) for a in ARMS}
             for k, r in loaded.items()}
    t64_means = {lib: float(np.mean(list(t64[lib].values()))) for lib in LIBRARIES}
    ahead, failed, values = statement(arm, ap, {lib: means[lib][ARM] for lib in LIBRARIES}, t64_means)
    idle = [k for k, r in loaded.items() if all(x["sha256"][ARM] == x["sha256"][BASE] for x in r["rows"])]
    blood = {k: {name: [int((x["records"][name]).get("blood_cells_changed", 0)) for x in r["rows"]]
                 for name in COUNTED} for k, r in loaded.items()}
    drawn = {k: [int(x["records"][ARM].get("blood_cells", 0)) for x in r["rows"]] for k, r in loaded.items()}
    return {"outcome": AHEAD_WORDS[ahead], "ahead": ahead, "failed": failed, "values": values,
            "tasks": {"arm_minus_base": arm, "perm_minus_base": perm, "arm_minus_perm": ap}, "mean_scores": means,
            "t64_means": t64_means, "valid_libraries": [lib for lib in LIBRARIES if library_valid(loaded[lib])],
            "idle_tasks": idle, "blood_drawn": drawn, "blood_changed": blood,
            "beside": list(BESIDE), "follows": FOLLOWS_AHEAD if ahead else FOLLOWS_END}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]]
    first = f"Outcome: {r['outcome']}" + ("" if r["ahead"] else " (" + "; ".join(r["failed"]) + ")")
    lines = [first, *r["beside"]]
    lines += [f"the arm moved no cell on {k}: nothing is priced there" for k in r["idle_tasks"]]
    return lines + [r["follows"]]


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} behind {p['behind']} of "
            f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}")


def t66_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        v = r["values"]
        lines.append(f"Clause values: library mean of the arm minus {BASE} {v['library_mean_d']:+.2f}; library mean "
                     f"of the arm minus {PERM} {v['library_mean_arm_minus_perm']:+.2f}; libraries at or above row "
                     f"T64's {T64_ARM}: {', '.join(v['libraries_at_or_above_t64']) or 'none'}")
        for title, key in ((f"{ARM} - {BASE}", "arm_minus_base"), (f"{PERM} - {BASE}", "perm_minus_base"),
                           (f"{ARM} - {PERM}", "arm_minus_perm")):
            lines.append(f"{title}, paired per seed on the same draw:")
            lines += [_line(k, p) for k, p in r["tasks"][key].items()]
        lines.append("Mean scores (reported, deciding nothing beyond clause (iii)):")
        for k, m in r["mean_scores"].items():
            t64 = f"; row T64's {T64_ARM} {r['t64_means'][k]:.2f}" if k in r["t64_means"] else ""
            lines.append(f"  {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in m.items()) + t64)
        for k in r["blood_drawn"]:
            lines.append(f"  {k}: Blood cells drawn per seed {r['blood_drawn'][k]}; changed per seed "
                         f"{r['blood_changed'][k]}")
        lines.append(f"Libraries without inverted anchors (reported only): {', '.join(r['valid_libraries']) or 'none'}")
    lines.append(f"Tries counted: {res['tries']} (two counted arms on five tasks; the row numbers them)")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, commit: str) -> dict:
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md")
    gv.refuse_existing(out, OUT_FILES)
    started = head_is(commit)
    paths = {k: out / "tasks" / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    t62 = {k: committed(f, T62_TASKS, BASE) for k, f in task_files().items()}
    t64 = {lib: committed(f"lib{lib}.json", T64_TASKS, T64_ARM) for lib in LIBRARIES}
    pg_path = out / "tasks" / PREGATE_FILES[0]
    pregate = json.loads(pg_path.read_text()) if pg_path.exists() else None
    mpath, _ = reducer_paths()
    manifest_sha = fio.sha256_file(mpath) if mpath.exists() else None
    checks = rc.relabel_checks(need_stats=False)  # the stats file's hash is in each task file, checked by `guards`
    bad = source_guard() + [f"recomputed_read.py: {x}" for x in rc.source_guard()]
    reading = t66_reading(loaded, t62, t64, checks, pregate, manifest_sha, bad, read_commit=started["commit"])
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT}
    res["stamp"] = gv.stamp(started, [p for p in [*paths.values(), pg_path] if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / OUT_FILES[0], json.dumps(res, indent=1, default=float) + "\n")
    write_new(out / OUT_FILES[1], t66_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pregate", help="try T66's target-free pre-gate")
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True)
    p = sub.add_parser("run-cut", help="one E10.5 cut of try T66")
    p.add_argument("--cut", required=True, choices=list(CUTS))
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("run-lib", help="one library task of try T66")
    p.add_argument("--library", required=True, choices=list(LIBRARIES))
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True)
    p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T66's reading")
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out, a.commit)
        return 0
    if a.cmd == "pregate":
        return cmd_pregate(a)
    {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
