#!/usr/bin/env python
"""Try T71 (experiments/harness/TRIES.md): the kill test of candidate C1A, "the CM direction blended with a
skeletal-guarded single-nucleus atlas direction", read on row T62's five practice tasks.

The recomputed Feng step moves the cardiomyocyte group along Feng's E9.5 -> E14.5 direction (`feng_dirs.
recomputed_directions`). This read replaces that group's entry, and only it, by the per-gene mean of Feng's direction and
the nuclei-atlas CM direction of experiments/external/jax_cm.py (row T71's reducer, practice ends: embryos 1-12 and
59-63), rescaled to Feng's RMS over the non-held finite genes, so the step's size is unchanged and only which genes move,
and by how much, changes (`blend`). Endocardium, Endothelium and Epicardium keep Feng's direction byte for byte.

The bases (gates, not counted): `m1rc_gmom` (row T62's candidate at 1/5 of the bracket) and `m1rc_gmom_x2` (row T63's
doubled step, 2/5), each rebuilt with the calls T62 and T63 make: the cuts through `t1_scope.build_submission`, the
libraries through `gmom_combo.moved_by_feng`, stream [seed, 2], then g_mom (`recomputed_read.rc_hook`, stream [seed, 5],
the same mask). Counted, all on one draw, at each fraction f in (1/5, 2/5):
  cmjax_x<k>      the base with the CM entry of the direction dict replaced by the blend
  cmjax_plc_x<k>  the same with the residual-permuting placebo in place of JAX's direction (`placebo`): D_plc = a x
                  D_feng + r_perm on G, G the non-held genes finite in both, a = <D_jax, D_feng>_G / <D_feng, D_feng>_G,
                  r = D_jax - a x D_feng permuted within 10 equal-count bins of the early end's JAX CM mean CP10k
                  (np.array_split on a stable rank; one permutation from default_rng(PLC_SEED)); NaN where D_jax is
                  NaN; then the same blend and rescale. It keeps what JAX shares with Feng and the residual's values,
                  and destroys only the gene identity of JAX's part beyond Feng.
How a counted arm is spliced, so that every non-CM cell keeps the base's bytes (the specification's byte check): the
direction step runs with the arm's dict on the same draw and stream; its CM rows replace the base's CM rows before g_mom
(`recomputed_read.hold_group` with group CM); g_mom runs as the base's on that population; its CM rows replace the
base's final CM rows. vec.core.apply_factors draws one stream over every moved cell, so a plain rebuild would redraw
non-CM cells whenever a CM factor crosses 1; the splice keeps them. With D_jax equal to D_feng every arm is its base,
byte for byte (tested). Each run checks the non-CM rows per seed and arm and stops otherwise (void as missing).

  pregate  Target-free: no T1 cell is loaded or scored and no truth is read. From the reducer's manifest and direction,
           the recomputed Feng statistics (CM means at each end and the CM direction) and the T1 gene names alone. It
           passes only if all hold: (a) at least MIN_END_NUCLEI CM nuclei enter at each end; (b) at the late end the
           share of r_cm nuclei whose skm score is at least the yaml's purity_score is at most MAX_SKM_SHARE; (c)
           Spearman(D_jax, D_feng_CM) over genes finite in both with Feng mean CP10k at least FENG_MIN_CP at either end
           lies in [SPEARMAN_LOW, SPEARMAN_HIGH]; (d) the share of D_jax's squared norm on T1 genes on R2's held genes is
           below MAX_HELD_SHARE; and the reducer's record is the row's. Reported only: each end's mean CP10k of the five
           skm genes in JAX CM nuclei beside Feng's E14.5 CM mean of them, nuclei per embryo, day labels. Exits 0 on a
           pass and PREGATE_FAIL on a failure (nothing is scored and no try is counted).
  run-cut  One E10.5 cut. Refused before anything loads unless row T71 is in TRIES.md, HEAD is --commit, the tree is
           clean, the task file does not exist, the source guard passes, the relabelling's committed checks pass and
           the pre-gate passed at this commit on the present reducer manifest. Prints no score.
  run-lib  One library task; the same refusals. Prints no score.
  read     Row T71's statement, applied before anything is printed; writes t71.json and t71.txt; prints the outcome
           words and what follows only.

Usage: through experiments/t1_m1/cmjax.sh, from the root of a checkout that holds data/derived/feng_relabel_c/, at the
row's commit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

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
import jax_cm as jc  # noqa: E402  (row T71's reducer: its constants, lists and yaml)
import recomputed_read as rc  # noqa: E402  (T62: its arms' pieces, guards and checks)

br, gc = rc.br, rc.gc
bt, gv, pt, rr = rc.bt, rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs
jb, fr, fio = jc.jb, jc.fr, jc.fio

ROW = "T71"
OUT = "t71_cmjax"
MODE = "practice"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T62_TASKS = RESULTS / "t62_recomputed" / "tasks"
T63_TASKS = RESULTS / "t63_double_step" / "tasks"
OUT_FILES = ("t71.json", "t71.txt")
PREGATE_FILES = ("pregate.json", "pregate.txt")
CM = "CM"  # the direction group whose entry the arms replace

BASE1 = rc.CANDIDATE  # m1rc_gmom (row T62's)
BASE2 = "m1rc_gmom_x2"  # row T63's doubled step
ARM1, PLC1, ARM2, PLC2 = "cmjax_x1", "cmjax_plc_x1", "cmjax_x2", "cmjax_plc_x2"
FRACTIONS = {"x1": (BASE1, ARM1, PLC1, 1.0), "x2": (BASE2, ARM2, PLC2, 2.0)}  # base, arm, placebo, multiplier
ARMS = (BASE1, ARM1, PLC1, BASE2, ARM2, PLC2)
GATED = (BASE1, BASE2)
COUNTED = (ARM1, PLC1, ARM2, PLC2)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first
LIB_ARMS = (LIB_REFERENCE, *ARMS)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 20

PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2
TEST_FRACTION = 0.4
FRACTION_TOL = 1e-12
REPRO_TOL = br.REPRO_TOL  # 0.02 board points
PLC_SEED = 20261123
PLC_BINS = 10
BLEND_WEIGHT = 0.5
# The pre-gate (fixed in the row)
MIN_END_NUCLEI = 200
MAX_SKM_SHARE = 0.10
SPEARMAN_LOW, SPEARMAN_HIGH = 0.2, 0.95
FENG_MIN_CP = 1.0
MAX_HELD_SHARE = 0.5
PREGATE_FAIL = 3
# The statement (fixed in the row), per fraction
CHAMBER_BAR = 0.3
LIB_BAR = 0.3
LIB_FLOOR = -0.3
PLC_CHAMBER_BAR = 0.2
PLC_LIB_BAR = 0.2
MIN_AHEAD = 8
AHEAD_LIB = "3"
BANNED_IMPORTS = ("perm_read", "eb_read", "step_read", "deep_read", "over_read", "lib_step_read", "ctl_range_read",
                  "n48_m1_direction_shrink", "n12_hygiene_size", "double_step_read", "heldmom_read", "jaxblood_read")
BANNED_PREFIX = "n5" + "4"

FOLLOWS_AHEAD = ("What follows: one priced question to the captain for each fraction read AHEAD ({fractions}) and for "
                 "that fraction only; its default is no change. Nothing is built or changed by this read.")
FOLLOWS_END = ("What follows: NOT AHEAD at both fractions closes 'JAX nuclei CM co-direction' on T1. No question goes to "
               "the captain, nothing is built, and nothing is rerun under another row.")
FOLLOWS_VOID = ("What follows: the tries still count; nothing is rerun under row T71, no question goes to the captain "
                "and nothing is built.")
FOLLOWS_PRUNED = "What follows: nothing is scored and no try is counted under row T71; the line closes at the pre-gate."
BESIDE = ("J2: Feng's E9.5 is inside every practice window.",
          "J5: a labelled in-sample diagnostic, never an adoption result.",
          "J6: the test spans two days, so at the test the gain counts zero at the centre.",
          "The practice early end is run_4 (E8.5) and the test's is run_15 (E8.75-E9.5); the run confound is declared, "
          "not removed.")
LABEL = ("try T71, a labelled in-sample diagnostic, never an adoption result: the CM group's recomputed Feng direction "
         "blended with the JAX/Qiu 2024 nuclei CM direction (RMS-matched) and its residual-permuting placebo, at 1/5 "
         "and 2/5, against m1rc_gmom (row T62's) and m1rc_gmom_x2 (row T63's) on T62's five tasks; every task scores "
         "one truth (E10.5) over one day")


# ---------------------------------------------------------------------------------------------------------------------
# Guards


def source_guard(path: Path | None = None) -> list[str]:
    """recomputed_read's guard on this file (every feng_directions call names the recomputed groups; no token of the
    other groups, their statistics or row N54's listing), no import of a script that forms the atlas direction on
    another route or of another row's read (BANNED_IMPORTS), of any row-N54 script, or any dynamic import."""
    path = Path(__file__ if path is None else path)
    bad = rc.source_guard(path)
    text = path.read_text(encoding="utf-8")
    for m in re.finditer(r"^\s*(?:from\s+(\S+)\s+import\s+([^\n#]+)|import\s+([^\n#]+))", text, re.M):
        names = [x for x in re.split(r"[\s,.()]+", " ".join(x for x in m.groups() if x)) if x]
        hit = sorted(set(names) & set(BANNED_IMPORTS)) + sorted(x for x in names if x.startswith(BANNED_PREFIX))
        if hit:
            bad.append(f"the source imports {hit}")
    for token in ("import" + "lib", "__imp" + "ort__", "spec_from" + "_file_location", "T1_s2" + ".h5ad"):
        if token in text:
            bad.append(f"the source names {token}")
    return bad


def head_is(commit: str) -> dict:
    started = gv.tree_state()
    if started["commit"] != commit:
        sys.exit(f"{ROW} is refused: HEAD is {started['commit']}, not the run commit {commit}")
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean by vec.provenance's rule")
    return started


def write_new(path: Path, text: str) -> None:
    """Exclusive create: a second copy of a step never replaces the file the first wrote."""
    with open(path, "x") as f:
        f.write(text)


def vector_sha(v: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(v, dtype=np.float64).tobytes()).hexdigest()


def prediction_sha(X) -> str:
    """The sha256 of an arm's scored matrix (cells x genes, float32, C order)."""
    return hashlib.sha256(np.ascontiguousarray(dense(X), dtype=np.float32).tobytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# The reducer's output


def reducer_paths() -> tuple[Path, Path]:
    cfg = jc.load_config()
    dst = jc.output_dir(cfg, MODE)
    return dst / "MANIFEST.json", dst / cfg["output"]["direction"]


def manifest_problems(man: Mapping) -> list[str]:
    """Every way the reducer's manifest is not the row's practice reduction."""
    bad = []
    if (man.get("config") or {}).get("sha256") != jc.YAML_SHA256:
        bad.append("the reducer ran under another yaml than the row's")
    gate = man.get("gate") or {}
    if gate.get("row") != ROW or gate.get("clean") is not True or gate.get("mode") != MODE:
        bad.append(f"the reducer's gate is not row {ROW}'s {MODE} run on a clean tree")
    if man.get("mode") != MODE:
        bad.append(f"the reducer's output is not its {MODE} mode's")
    if man.get("gene_universe") != jc.GENE_UNIVERSE:
        bad.append("the reducer's gene universe is not the gene table's genes")
    cols = list(man.get("columns_read") or [])
    if not cols or any(c not in jc.CELL_WHITELIST or fr.ANNOTATION_PATTERN.search(c) for c in cols):
        bad.append(f"the reducer read columns outside the whitelist: {cols}")
    if man.get("authors_labels_read") is not False or man.get("stats_npz_read") is not False:
        bad.append("the reducer's record does not say that no author label and no stats file was read")
    lists = jc.ROW_EMBRYOS[MODE]
    if sorted(int(x) for x in man.get("embryos_read") or []) != sorted(n for e in jc.ENDS for n in lists[e]):
        bad.append("an embryo outside the row's lists entered, or one is missing")
    got = man.get("embryo_lists") or {}
    if {e: sorted(got.get(e) or []) for e in jc.ENDS} != {e: list(lists[e]) for e in jc.ENDS}:
        bad.append("the reducer's embryo lists are not the row's")
    used = man.get("embryos_used") or {}
    if any(int(n) not in lists[e] for e in jc.ENDS for n in used.get(e) or []):
        bad.append("an embryo outside the row's lists entered the direction")
    if man.get("has_direction") is not True:
        bad.append("the reducer wrote no direction (an end had too few embryos at the floor)")
    return bad


def load_reducer(manifest: Path | None = None, direction: Path | None = None) -> tuple[dict, dict, str]:
    if manifest is None or direction is None:
        manifest, direction = reducer_paths()
    man = json.loads(Path(manifest).read_text())
    fio.verify_manifest_output(direction)
    with np.load(direction, allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    return man, arrays, fio.sha256_file(manifest)


def on_t1(arrays: Mapping, genes: Sequence[str], key: str = "direction") -> np.ndarray:
    """A per-JAX-gene array on `genes` through JAX's unique symbols (a symbol the gene table holds twice gets none)."""
    uidx = fr.unique_index([str(s) for s in arrays["symbol"]])
    vals = np.asarray(arrays[key], dtype=np.float64)
    return np.array([vals[uidx[g]] if g in uidx else np.nan for g in genes], dtype=np.float64)


# ---------------------------------------------------------------------------------------------------------------------
# The blend and the placebo (target-free)


def _rms(v: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(v)))) if v.size else float("nan")


def blend(d_feng: np.ndarray, d_jax: np.ndarray, free: np.ndarray) -> tuple[np.ndarray, float]:
    """D_b = 0.5 D_feng + 0.5 D_jax where both are finite, else D_feng; then D_b x RMS(D_feng) / RMS(D_b) over the
    non-held (`free`) genes where D_b is finite. Returns D_b and the factor. Refused if RMS(D_b) is 0 there."""
    d_feng, d_jax = np.asarray(d_feng, dtype=np.float64), np.asarray(d_jax, dtype=np.float64)
    both = np.isfinite(d_feng) & np.isfinite(d_jax)
    out = d_feng.copy()
    out[both] = BLEND_WEIGHT * d_feng[both] + (1.0 - BLEND_WEIGHT) * d_jax[both]
    m = np.asarray(free, dtype=bool) & np.isfinite(out)
    den = _rms(out[m])
    if not den > 0:
        raise ValueError("the blend has no size on the non-held genes")
    k = _rms(d_feng[m]) / den
    return out * k, k


def placebo(d_feng: np.ndarray, d_jax: np.ndarray, level: np.ndarray, free: np.ndarray, seed: int = PLC_SEED,
            bins: int = PLC_BINS) -> tuple[np.ndarray, float, np.ndarray]:
    """The residual-permuting placebo before the blend: G = the non-held genes with D_jax, D_feng and the level finite;
    a = <D_jax, D_feng>_G / <D_feng, D_feng>_G; r = D_jax - a D_feng on G, permuted within `bins` equal-count bins of
    `level` (the early end's JAX CM mean CP10k; ranked stably, ties by gene position; np.array_split); one generator,
    default_rng(seed), permutes each bin in order. D_plc = a D_feng + r_perm on G and D_jax elsewhere (NaN where D_jax
    is NaN). Returns D_plc, a and the permuted residual (on G; r elsewhere)."""
    d_feng, d_jax = np.asarray(d_feng, dtype=np.float64), np.asarray(d_jax, dtype=np.float64)
    level = np.asarray(level, dtype=np.float64)
    g = np.asarray(free, dtype=bool) & np.isfinite(d_jax) & np.isfinite(d_feng) & np.isfinite(level)
    idx = np.flatnonzero(g)
    den = float(np.sum(d_feng[idx] ** 2))
    if not den > 0:
        raise ValueError("D_feng has no size on G")
    a = float(np.sum(d_jax[idx] * d_feng[idx])) / den
    r = d_jax - a * d_feng
    out = d_jax.copy()
    order = idx[np.argsort(level[idx], kind="stable")]
    rng = np.random.default_rng(seed)
    perm = r.copy()
    for b in np.array_split(order, bins):
        perm[b] = r[b[rng.permutation(b.size)]]
    out[idx] = a * d_feng[idx] + perm[idx]
    return out, a, perm


def arm_directions(directions: Mapping, d_cm: np.ndarray) -> dict:
    """The direction dict with the CM entry replaced (every other entry the same object, so byte for byte)."""
    if CM not in directions:
        raise ValueError("the recomputed directions have no CM entry")
    return {g: (np.asarray(d_cm, dtype=np.float64) if g == CM else v) for g, v in directions.items()}


def jax_cm_inputs(directions: Mapping, arrays: Mapping, genes: Sequence[str], held_m1: Sequence[str]) -> dict:
    """D_jax and the early level on `genes`; the arm's and the placebo's CM entries; the fields the task file records."""
    free = ~np.isin(np.asarray(genes, dtype=str), list(held_m1))
    d_feng = np.asarray(directions[CM], dtype=np.float64)
    d_jax = on_t1(arrays, genes)
    level = on_t1(arrays, genes, "m_early")
    d_arm, k_arm = blend(d_feng, d_jax, free)
    d_plc0, a, r_perm = placebo(d_feng, d_jax, level, free)
    d_plc, k_plc = blend(d_feng, d_plc0, free)
    g = free & np.isfinite(d_jax) & np.isfinite(d_feng) & np.isfinite(level)
    r = d_jax - a * d_feng
    fields = {"d_jax_sha256": vector_sha(d_jax), "d_arm_sha256": vector_sha(d_arm), "d_plc_sha256": vector_sha(d_plc),
              "rescale_arm": k_arm, "rescale_plc": k_plc, "a": a, "g_genes": int(g.sum()),
              "residual_same_values": bool(np.array_equal(np.sort(r[g]), np.sort(r_perm[g]))),
              "plc_nan_where_jax_nan": bool(np.array_equal(np.isnan(d_plc0), np.isnan(d_jax)))}
    return {"d_jax": d_jax, "arm": d_arm, "plc": d_plc, "fields": fields}


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate (target-free)


def feng_cm_means(genes: Sequence[str]) -> dict[str, np.ndarray]:
    """Feng's recomputed CM group mean CP10k at each bracket end on `genes` (any symbol list), through feng_dirs' own
    checks of the recomputed file and its group map."""
    rcfg = feng_dirs.load_relabel_config()
    path = feng_dirs.recomputed_stats_path(rcfg)
    feng_dirs.check_recomputed(path)
    st = feng_dirs.dl.load_stats(path)
    feng_dirs.check_recomputed_keys(st)
    cfg = feng_dirs.groups_config(feng_dirs.load_config(), rcfg)
    labels = cfg["external_labels"][CM][0]
    idx = feng_dirs.dl.sym_index(st["genes"], list(genes))
    out = {}
    for stage in feng_dirs.BRACKET:
        sel = feng_dirs.dl.feng_sel(st, stage, labels)
        if sel is None:
            raise ValueError(f"Feng's recomputed CM group has no cell at {stage}")
        out[stage] = feng_dirs.dl.on_t1(sel["cp"], idx)
    return out


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(x, y).statistic)


def pregate_reading(man: Mapping, d_jax: np.ndarray, d_feng: np.ndarray, feng_cp: Mapping[str, np.ndarray],
                    held: np.ndarray, skm: Mapping | None = None) -> dict:
    """The pre-gate's clauses (a) to (d) and the record check, from the manifest, D_jax and D_feng on T1's genes, Feng's
    CM means at each end on T1's genes and R2's held mask; `skm` is the reported skm table (deciding nothing)."""
    per = man.get("per_embryo") or {}
    used = man.get("embryos_used") or {}
    ends = {}
    for e in jc.ENDS:
        rows = {int(n): r for n, r in per.items() if int(n) in jc.ROW_EMBRYOS[MODE][e]}
        contrib = {int(n) for n in used.get(e) or []}
        ends[e] = {"assigned": {n: int(r.get("cm_assigned", 0)) for n, r in sorted(rows.items())},
                   "entering": {n: int(r.get("cm_entering", 0)) for n, r in sorted(rows.items())},
                   "skm_high": {n: int(r.get("cm_skm_high", 0)) for n, r in sorted(rows.items())},
                   "days": {n: r.get("day") for n, r in sorted(rows.items())},
                   "contributing": sorted(contrib),
                   "entering_total": int(sum(int(r.get("cm_entering", 0)) for n, r in rows.items() if n in contrib))}
        a = sum(ends[e]["assigned"].values())
        ends[e]["skm_share"] = (sum(ends[e]["skm_high"].values()) / a) if a > 0 else float("nan")
    failed = [f"problem with the reducer's record: {x}" for x in manifest_problems(man)]
    for e in jc.ENDS:
        if not ends[e]["entering_total"] >= MIN_END_NUCLEI:
            failed.append(f"(a) fewer than {MIN_END_NUCLEI} CM nuclei enter at the {e} end")
    share = ends["late"]["skm_share"]
    if not share <= MAX_SKM_SHARE:
        failed.append(f"(b) the late end's share of r_cm nuclei with an skm score at the purity bar is not at most "
                      f"{MAX_SKM_SHARE:g}")
    expressed = np.zeros(len(d_jax), dtype=bool)
    for v in feng_cp.values():
        expressed |= np.nan_to_num(np.asarray(v, dtype=np.float64), nan=-1.0) >= FENG_MIN_CP
    both = np.isfinite(d_jax) & np.isfinite(d_feng) & expressed
    rho = _spearman(d_jax[both], d_feng[both])
    if not SPEARMAN_LOW <= rho <= SPEARMAN_HIGH:
        failed.append(f"(c) Spearman(D_jax, D_feng_CM) is not within [{SPEARMAN_LOW:g}, {SPEARMAN_HIGH:g}]")
    sq = np.where(np.isfinite(d_jax), d_jax, 0.0) ** 2
    total = float(sq.sum())
    held_share = float(sq[np.asarray(held, dtype=bool)].sum() / total) if total > 0 else float("nan")
    if not held_share < MAX_HELD_SHARE:
        failed.append(f"(d) the share of D_jax's squared norm on R2's held genes is not below {MAX_HELD_SHARE:g}")
    return {"passes": not failed, "failed": failed, "ends": ends, "late_skm_share": share, "spearman": rho,
            "spearman_genes": int(both.sum()), "held_share": held_share, "delta_days": man.get("delta_days"),
            "end_left_out_umi_share": man.get("end_left_out_umi_share"), "skm_means": skm}


def skm_table(arrays: Mapping, feng_late: np.ndarray) -> dict:
    """Reported only: each end's mean CP10k of the five skm genes in JAX CM nuclei, beside Feng's E14.5 CM mean."""
    early, late = on_t1(arrays, jc.SKM_GENES, "m_early"), on_t1(arrays, jc.SKM_GENES, "m_late")
    return {g: {"jax_early": float(early[i]), "jax_late": float(late[i]), "feng_e14_5": float(feng_late[i])}
            for i, g in enumerate(jc.SKM_GENES)}


def pregate_outcome(r: Mapping) -> str:
    return "Pre-gate: passes" if r["passes"] else "Pre-gate: FAILS (" + "; ".join(r["failed"]) + ")"


def pregate_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], pregate_outcome(r)]
    for e, v in r["ends"].items():
        lines.append(f"{e} end: CM nuclei entering (contributing embryos {v['contributing']}) {v['entering_total']}; "
                     f"assigned per embryo {v['assigned']}; entering per embryo {v['entering']}; with an skm score at "
                     f"the bar per embryo {v['skm_high']} (share {v['skm_share']:.3f}); day labels {v['days']}")
    lines.append(f"Spearman(D_jax, D_feng_CM) {r['spearman']:.3f} over {r['spearman_genes']} genes; R2-held share of "
                 f"D_jax's squared norm {r['held_share']:.3f}; Delta {r['delta_days']} days; each end's UMI share on "
                 f"genes outside the gene table {r['end_left_out_umi_share']} (reported only)")
    lines.append(f"skm genes, mean CP10k in JAX CM nuclei and Feng's E14.5 CM (reported only): {r['skm_means']}")
    lines.append(BESIDE[3])
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


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
    directions, _ = rc.recomputed_inputs(genes)
    held_m1, _ = rc._held_m1(t1_scope.held_flat_genes())
    held = np.isin(np.asarray(genes, dtype=str), held_m1)
    feng_cp = feng_cm_means(genes)
    skm = skm_table(arrays, feng_cm_means(list(jc.SKM_GENES))[feng_dirs.BRACKET[1]])
    reading = pregate_reading(man, on_t1(arrays, genes), np.asarray(directions[CM], dtype=np.float64), feng_cp, held,
                              skm)
    mpath, dpath = reducer_paths()
    res = {"label": LABEL, "row": ROW, "reading": reading, "manifest_sha256": man_sha,
           "yaml_sha256": (man.get("config") or {}).get("sha256"),
           "stamp": gv.stamp(started, [mpath, dpath, t1_scope.HELD_FLAT_GENES, feng_dirs.recomputed_stats_path()])}
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / PREGATE_FILES[0], json.dumps(res, indent=1, default=float) + "\n")
    write_new(out / PREGATE_FILES[1], pregate_text(res))
    print("Pre-gate: passes" if reading["passes"] else f"Pre-gate: FAILS\n{FOLLOWS_PRUNED}")  # words only
    return 0 if reading["passes"] else PREGATE_FAIL


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def cm_record(labels: np.ndarray, groups: Mapping[str, str], base_X, arm_X, name: str) -> dict:
    """After the build: every row whose label is not in the CM group has the base's bytes in the arm, or the run
    stops."""
    cm = np.array([groups.get(str(lab)) == CM for lab in labels], dtype=bool)
    a = np.ascontiguousarray(dense(base_X), dtype=np.float32)
    b = np.ascontiguousarray(dense(arm_X), dtype=np.float32)
    if a.shape != b.shape or a.shape[0] != len(labels):
        raise SystemExit(f"{ROW}: the base and {name} are not the same cells ({a.shape}, {b.shape}, {len(labels)})")
    if a[~cm].tobytes() != b[~cm].tobytes():
        raise SystemExit(f"{ROW}: a cell outside the CM group differs between its base and {name}; the run stops")
    return {"cm_cells": int(cm.sum()), "cm_cells_changed": int((a[cm] != b[cm]).any(axis=1).sum()),
            "non_cm_equal": True}


def splice(base_moved, base_final, arm_moved, hook, groups: Mapping[str, str]):
    """The arm's population: the arm's direction-moved CM rows in the base's moved population, then g_mom (`hook`),
    then those CM rows in the base's final population (`recomputed_read.hold_group` with group CM)."""
    if list(arm_moved.obs_names) != list(base_moved.obs_names):
        raise SystemExit(f"{ROW}: an arm and its base were not built on one draw")
    moved = rc.hold_group(base_moved, arm_moved, groups, CM)
    return rc.hold_group(base_final, hook(moved), groups, CM)


def cut_arms(ctx: harness.BuildContext, factors, directions: Mapping, fraction: float, held_m1: Sequence[str],
             d_arm: np.ndarray, d_plc: np.ndarray) -> tuple[dict[str, Stage], dict]:
    """The six arms of one seed on a harness cut, from one draw, each finished as the harness finishes a candidate,
    and the per-arm records."""
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

    def build(dirs, f, post) -> Stage:
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src,
                                           held_flat=list(held_m1), directions=dirs, direction_fraction=f,
                                           direction_groups=groups, post_cells=post)
        return harness.finish(ctx, ref, dense(sub.X), None)

    out, rec = {}, {}
    for base_name, arm, plc, mult in FRACTIONS.values():
        f = mult * float(fraction)
        cap = {}

        def base_post(pop, cap=cap):
            cap["moved"] = pop
            cap["final"] = hook(pop)
            return cap["final"]

        out[base_name] = build(directions, f, base_post)
        drawn_labels = cap["moved"].obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
        for name, d in ((arm, d_arm), (plc, d_plc)):
            def post(pop, cap=cap):
                return splice(cap["moved"], cap["final"], pop, hook, groups)

            out[name] = build(arm_directions(directions, d), f, post)
            rec[name] = cm_record(drawn_labels, groups, out[base_name].X, out[name].X, name)
    return out, rec


def library_arms(e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                 d_arm: np.ndarray, d_plc: np.ndarray, n_cells: int = rr.N_CELLS) -> tuple[dict, dict]:
    """The six arms of one seed on a library task from T60's oracle draw (recomputed_read.library_pair's route), as
    AnnData, and the per-arm records."""
    base0, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    groups, covered, with_direction = rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), directions)
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    labels = base0.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    out, rec = {}, {}
    for base_name, arm, plc, mult in FRACTIONS.values():
        f = mult * float(fraction)
        moved = gc.moved_by_feng(base0, directions, f, groups, list(held_m1), [seed, 2])
        out[base_name] = hook(moved)
        for name, d in ((arm, d_arm), (plc, d_plc)):
            arm_moved = gc.moved_by_feng(base0, arm_directions(directions, d), f, groups, list(held_m1), [seed, 2])
            out[name] = splice(moved, out[base_name], arm_moved, hook, groups)
            rec[name] = cm_record(labels, groups, out[base_name].X, out[name].X, name)
    return out, rec


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def refuse(out: Path, name: str, commit: str) -> tuple[dict, dict, dict]:
    """The refusals before anything loads; returns the start state, the relabelling's checks and the pre-gate."""
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    gv.refuse_existing(out, (name,))
    started = head_is(commit)
    bad = source_guard()
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


def task_inputs(directions: Mapping, genes: Sequence[str], held_m1: Sequence[str]) -> dict:
    man, arrays, man_sha = load_reducer()
    bad = manifest_problems(man)
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    got = jax_cm_inputs(directions, arrays, genes, held_m1)
    _, dpath = reducer_paths()
    got["fields"].update({"manifest_sha256": man_sha, "yaml_sha256": (man.get("config") or {}).get("sha256"),
                          "direction_sha256": fio.sha256_file(dpath), "mode": man.get("mode"),
                          "embryos_read": man.get("embryos_read"), "columns_read": man.get("columns_read"),
                          "authors_labels_read": man.get("authors_labels_read"),
                          "stats_npz_read": man.get("stats_npz_read")})
    return got


def constants(held_m1: Sequence[str]) -> dict:
    return {"multipliers": {k: v[3] for k, v in FRACTIONS.items()}, "plc_seed": PLC_SEED, "plc_bins": PLC_BINS,
            "blend_weight": BLEND_WEIGHT, "held_m1_genes": len(held_m1),
            "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES)}


def check_fractions(fraction: float, where: str) -> None:
    if abs(fraction - PRACTICE_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")
    if abs(feng_dirs.bracket_fraction(10.5, 12.5) - TEST_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the test's fraction is not {TEST_FRACTION}")


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
    genes_ref = list(inputs[ref_role].var_names)
    directions, ext = rc.recomputed_inputs(genes_ref)
    fraction = feng_dirs.bracket_fraction(harness.STAGE_DAYS[t.reference], harness.STAGE_DAYS[t.target])
    check_fractions(fraction, t.key)
    held_m1, union = rc._held_m1(preregistered)
    ji = task_inputs(directions, genes_ref, held_m1)
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for s in seeds:
        t0 = time.time()
        f = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms, rec = cut_arms(ctx, factors, directions, fraction, held_m1, ji["arm"], ji["plc"])
        shas = {k: prediction_sha(v.X) for k, v in arms.items()}
        rows.append({"seed": s, "floor": float(f.score), "arms": gv.score_arms(board, target, ref_stage, s, f.anchors,
                                                                               arms),
                     "sha256": shas, "records": rec, "seconds": round(time.time() - t0, 1)})
        del arms
        print(f"{t.key} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    mpath, dpath = reducer_paths()
    res = {"label": LABEL, "row": ROW, "task": t.key, "seeds": list(seeds), "rows": rows, "feng_fraction": fraction,
           "relabel": rc._relabel_fields(ext, checks), "jax": ji["fields"], "pregate_sha256": pg["manifest_sha256"],
           "constants": constants(held_m1),
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, ROOT / union["file"], mpath, dpath,
                                       jc.CONFIG])}
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
    check_fractions(fraction, f"library {a.library}")
    held_m1, union = rc._held_m1(held)
    ji = task_inputs(directions, list(e95.var_names), held_m1)
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms, rec = library_arms(e95, shares, factors, held, held_m1, s, within, fb, directions, fraction, ji["arm"],
                                 ji["plc"])
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        shas = {k: prediction_sha(v.X) for k, v in stages.items()}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s), "sha256": shas,
                     "records": rec, "seconds": round(time.time() - t0, 1)})
        del arms, stages
        print(f"library {a.library} seed {s}: {time.time() - t0:.0f} s", flush=True)  # no score
    mpath, dpath = reducer_paths()
    res = {"label": LABEL, "row": ROW, "library": a.library, "seeds": list(seeds), "rows": rows,
           "feng_fraction": fraction, "relabel": rc._relabel_fields(ext, checks), "jax": ji["fields"],
           "pregate_sha256": pg["manifest_sha256"],
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD, **constants(held_m1)},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], feng_dirs.RELABEL_CONFIG,
                                       rc.RECORD, ROOT / union["file"], mpath, dpath, jc.CONFIG])}
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / name, json.dumps(res, indent=1, default=float) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def task_files() -> dict[str, str]:
    return rc.task_files()


def committed(name: str, tasks: Path, arm: str) -> dict[int, float]:
    rows = json.loads((tasks / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][arm]["score"]) for r in rows}


def guards(loaded: Mapping[str, Mapping | None], gates: Mapping[str, Mapping[str, Mapping[int, float]]],
           checks: Mapping, pregate: Mapping | None, manifest_sha: str | None, source_bad: Sequence[str],
           seeds: Sequence[int] = SEEDS, read_commit: str | None = None) -> list[str]:
    """Every void condition of row T71; `gates` maps each gated base to its committed scores per task and seed (T62's
    m1rc_gmom, T63's m1rc_gmom_x2)."""
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
        for base in GATED:
            ref = gates.get(base, {}).get(k, {})
            worst = max((abs(br.score_of(x, base) - ref[int(x["seed"])]) if int(x["seed"]) in ref else float("inf")
                         for x in rows), default=float("inf"))
            if not worst <= REPRO_TOL:
                reasons.append(f"{k}: the gate failed: {base} differs from its committed scores by up to {worst:.3g} "
                               f"(tolerance {REPRO_TOL:g})")
        if not abs(float(r.get("feng_fraction", float("nan"))) - PRACTICE_FRACTION) <= FRACTION_TOL:
            reasons.append(f"{k}: the base's fraction is {r.get('feng_fraction')}, not {PRACTICE_FRACTION}")
        c = r.get("constants") or {}
        if (c.get("multipliers"), c.get("plc_seed"), c.get("plc_bins"), c.get("blend_weight")) != (
                {k2: v[3] for k2, v in FRACTIONS.items()}, PLC_SEED, PLC_BINS, BLEND_WEIGHT):
            reasons.append(f"{k}: a constant is not the row's")
        rel = r.get("relabel") or {}
        if rel.get("sha256") != rc.STATS_SHA256 or rel.get("stats_sha256_hashed") != rc.STATS_SHA256:
            reasons.append(f"{k}: ran on another recomputed statistics file")
        if rel.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: ran under another relabelling settings file")
        if rel.get("route") != rc.GROUPS_OPTION or rel.get("record_route") != rc.GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        j = r.get("jax") or {}
        if j.get("yaml_sha256") != jc.YAML_SHA256 or j.get("mode") != MODE:
            reasons.append(f"{k}: the reducer's yaml sha256 or mode is not the row's")
        if not manifest_sha or j.get("manifest_sha256") != manifest_sha or r.get("pregate_sha256") != manifest_sha:
            reasons.append(f"{k}: the reducer's manifest sha256 is not the pre-gate's and the present one")
        cols = list(j.get("columns_read") or [])
        if not cols or any(x not in jc.CELL_WHITELIST for x in cols) or j.get("authors_labels_read") is not False \
                or j.get("stats_npz_read") is not False:
            reasons.append(f"{k}: a label column or a label-keyed file was read")
        want_e = sorted(n for e in jc.ENDS for n in jc.ROW_EMBRYOS[MODE][e])
        if sorted(int(x) for x in j.get("embryos_read") or []) != want_e:
            reasons.append(f"{k}: an embryo outside the row's lists entered, or one is missing")
        if j.get("residual_same_values") is not True or j.get("plc_nan_where_jax_nan") is not True:
            reasons.append(f"{k}: the placebo does not keep the residual's values or D_jax's missing genes")
        for x in rows:
            for name in COUNTED:
                h = (x.get("records") or {}).get(name) or {}
                if h.get("non_cm_equal") is not True:
                    reasons.append(f"{k}: seed {x['seed']}: a cell outside the CM group is not byte-equal in {name}")
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
    p = br.paired(rows, change, base)
    p["behind"] = int(sum(br.score_of(r, change) - br.score_of(r, base) < 0 for r in rows))
    return p


def fraction_statement(ab: Mapping[str, Mapping], ap: Mapping[str, Mapping]) -> tuple[bool, list[str], dict]:
    """One fraction: (i) arm minus base, (ii) arm minus placebo (paired per seed). AHEAD when both hold."""
    lib_ab = float(np.mean([ab[lib]["mean"] for lib in LIBRARIES]))
    lib_ap = float(np.mean([ap[lib]["mean"] for lib in LIBRARIES]))
    failed = []
    if not (ab["chamber"]["mean"] >= CHAMBER_BAR and ab["chamber"]["ahead"] >= MIN_AHEAD):
        failed.append(f"(i) the chamber cut's mean of the arm minus its base is below +{CHAMBER_BAR:g} or fewer than "
                      f"{MIN_AHEAD} of 10 seeds are ahead")
    if not ab["whole"]["mean"] > 0:
        failed.append("(i) the whole cut's mean of the arm minus its base is not above 0")
    if not lib_ab >= LIB_BAR:
        failed.append(f"(i) the library mean of the arm minus its base is below +{LIB_BAR:g}")
    if not ab[AHEAD_LIB]["ahead"] >= MIN_AHEAD:
        failed.append(f"(i) library {AHEAD_LIB} is ahead of its base in fewer than {MIN_AHEAD} of 10 seeds")
    low = [lib for lib in LIBRARIES if not ab[lib]["mean"] >= LIB_FLOOR]
    if low:
        failed.append(f"(i) libraries {low} have a mean of the arm minus its base below {LIB_FLOOR:g}")
    if not (ap["chamber"]["mean"] >= PLC_CHAMBER_BAR and ap["chamber"]["ahead"] >= MIN_AHEAD):
        failed.append(f"(ii) the chamber cut's mean of the arm minus its placebo is below +{PLC_CHAMBER_BAR:g} or "
                      f"fewer than {MIN_AHEAD} of 10 seeds are ahead")
    if not lib_ap >= PLC_LIB_BAR:
        failed.append(f"(ii) the library mean of the arm minus its placebo is below +{PLC_LIB_BAR:g}")
    return not failed, failed, {"library_mean_arm_minus_base": lib_ab, "library_mean_arm_minus_placebo": lib_ap}


def t71_reading(loaded: Mapping[str, Mapping | None], gates: Mapping, checks: Mapping, pregate: Mapping | None,
                manifest_sha: str | None, source_bad: Sequence[str] = (), read_commit: str | None = None) -> dict:
    reasons = guards(loaded, gates, checks, pregate, manifest_sha, source_bad, read_commit=read_commit)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons, "follows": FOLLOWS_VOID}
    per, tasks = {}, {}
    for key, (base, arm, plc, _) in FRACTIONS.items():
        ab = {k: paired(r["rows"], arm, base) for k, r in loaded.items()}
        ap = {k: paired(r["rows"], arm, plc) for k, r in loaded.items()}
        pb = {k: paired(r["rows"], plc, base) for k, r in loaded.items()}
        ahead, failed, values = fraction_statement(ab, ap)
        idle = [k for k, r in loaded.items() if all(x["sha256"][arm] == x["sha256"][base] for x in r["rows"])]
        per[key] = {"ahead": ahead, "failed": failed, "values": values, "idle_tasks": idle}
        tasks[key] = {"arm_minus_base": ab, "arm_minus_placebo": ap, "placebo_minus_base": pb}
    passing = [k for k, v in per.items() if v["ahead"]]
    means = {k: {a: float(np.mean([br.score_of(x, a) for x in r["rows"]])) for a in ARMS} for k, r in loaded.items()}
    outcome = "; ".join(f"{k}: {'AHEAD' if v['ahead'] else 'NOT AHEAD'}" for k, v in per.items())
    return {"outcome": outcome, "fractions": per, "passing": passing, "tasks": tasks, "mean_scores": means,
            "cm_changed": {k: {n: [int(x["records"][n].get("cm_cells_changed", 0)) for x in r["rows"]]
                               for n in COUNTED} for k, r in loaded.items()},
            "beside": list(BESIDE),
            "follows": FOLLOWS_AHEAD.format(fractions=", ".join(passing)) if passing else FOLLOWS_END}


def outcome_lines(r: Mapping) -> list[str]:
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]]
    lines = []
    for k, v in r["fractions"].items():
        lines.append(f"Outcome at {k}: " + ("AHEAD" if v["ahead"] else "NOT AHEAD (" + "; ".join(v["failed"]) + ")"))
        lines += [f"at {k} the arm moved no cell on {t}: nothing is priced there" for t in v["idle_tasks"]]
    return [*lines, *r["beside"], r["follows"]]


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} behind {p['behind']} of "
            f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}")


def t71_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], *outcome_lines(r)]
    if r["outcome"] != "VOID":
        for key, t in r["tasks"].items():
            v = r["fractions"][key]["values"]
            lines.append(f"{key}: library mean of the arm minus its base {v['library_mean_arm_minus_base']:+.2f}; of "
                         f"the arm minus its placebo {v['library_mean_arm_minus_placebo']:+.2f}")
            for title, part in t.items():
                lines.append(f"{key} {title.replace('_', ' ')}, paired per seed on the same draw:")
                lines += [_line(k, p) for k, p in part.items()]
        lines.append("Mean scores (reported):")
        for k, m in r["mean_scores"].items():
            lines.append(f"  {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in m.items()))
        for k, c in r["cm_changed"].items():
            lines.append(f"  {k}: CM cells changed per seed {c}")
    lines.append(f"Tries counted: {res['tries']} (four counted arms on five tasks; the row numbers them)")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, commit: str) -> dict:
    if ROW not in br.tries_rows(TRIES):
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md")
    gv.refuse_existing(out, OUT_FILES)
    started = head_is(commit)
    paths = {k: out / "tasks" / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    gates = {BASE1: {k: committed(f, T62_TASKS, BASE1) for k, f in task_files().items()},
             BASE2: {k: committed(f, T63_TASKS, BASE2) for k, f in task_files().items()}}
    pg_path = out / "tasks" / PREGATE_FILES[0]
    pregate = json.loads(pg_path.read_text()) if pg_path.exists() else None
    mpath, _ = reducer_paths()
    manifest_sha = fio.sha256_file(mpath) if mpath.exists() else None
    checks = rc.relabel_checks(need_stats=False)
    reading = t71_reading(loaded, gates, checks, pregate, manifest_sha, source_guard(), read_commit=started["commit"])
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT}
    res["stamp"] = gv.stamp(started, [p for p in [*paths.values(), pg_path] if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / OUT_FILES[0], json.dumps(res, indent=1, default=float) + "\n")
    write_new(out / OUT_FILES[1], t71_text(res))
    print("\n".join(outcome_lines(reading)) + f"\n(numbers in {out / OUT_FILES[1]})")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pregate", help="try T71's target-free pre-gate")
    p.add_argument("--out", required=True)
    p.add_argument("--commit", required=True)
    for cmd, flag, choices in (("run-cut", "--cut", list(CUTS)), ("run-lib", "--library", list(LIBRARIES))):
        p = sub.add_parser(cmd)
        p.add_argument(flag, required=True, choices=choices)
        p.add_argument("--out", required=True)
        p.add_argument("--commit", required=True)
        p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
        p.add_argument("--data", default=str(harness.DATA_DIR))
    p = sub.add_parser("read", help="try T71's reading")
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
