#!/usr/bin/env python
"""XPC, the heart's run-profile table (row H36 of experiments/harness/TRIES.md; ruling row R56, the captain's named
exception for file OXP125 and his A15 exception, bound to row H36's read): E10.5's own MERFISH run profile measured
against the organisers' Multiome RNA at the stages both instruments hold (E8.5, E9.5, E10.5), the biology cancelled by
stage differences (the step-change round's HEART MODEL design, proposal P1, with every fix of its cold challenge).

Only row H36's steps call this module (experiments/harness/h36_heart_run_profile.py, through
experiments/harness/run_h36.sh): its entry (E), its practice profile (Q: E8.5 and E9.5 only, no E10.5 file) and,
only after the practice read's outcome is committed, its test profile (T), which writes PROFILE, the file the heart
builder's `--run-profile <sha256>` takes. It has no command line of its own: reading E10.5 here is a counted reading
(the challenge's fix 3), so it runs only where row H36 puts it.

Definitions, fixed before any read:
  categories k   cm, endocardium, progenitor (the challenge's fix 2: no stand-ins). MERFISH cells by
                 configs/t2_chamber_classes.yaml (vec.t2_build.t2_chamber_categories); RNA cells by
                 configs/t1_chamber_classes.yaml with the file's own per-file overrides (vec.t1_scope.classify; no
                 file_fallbacks, which only re-weighting reads)
  cells          MERFISH: the released file as the heart builder loads it (vec.prep.load_stage, the board's panel,
                 duplicate records dropped); RNA: the released file (load_stage, lineage=False), each cell's counts
                 restricted to the panel genes. E10.5's RNA is library 3 alone (cell names end in _3: the dissected
                 hearts); E8.5's and E9.5's both libraries. A cell with no panel count is left out
  M_skg, R_skg   log(mean over the group's cells of the gene's panel share + 1e-5), the share a gene's count over the
                 cell's panel count; MERFISH and RNA
  D_skg          M_skg - R_skg: each gene's constant platform term cancels
  delta_k        test: D_E10.5,k - (D_E8.5,k + D_E9.5,k) / 2; practice: D_E9.5,k - D_E8.5,k (no E10.5 file)
  delta          the median over the categories; 0 where the gene's mean CP10k at the reference stage (E10.5 for the
                 test, E9.5 for the practice profile) over the cells of the three categories is under 0.5 on either
                 instrument (the eligible genes are the others); clipped to |delta| <= log 3
  entry          every category at least 200 cells on both instruments at E8.5, E9.5 and E10.5 (library 3)
  halves         every stage x instrument x category group split in two at seed 20261125, one stream per group
                 (default_rng([seed, crc32(group name)])); delta^A and delta^B by the same rule on each half, with the
                 whole data's eligible genes
  checks         over the eligible genes (the challenge's fix 4 in place of the design's G0(d)):
                   noise     SD(delta) >= 2 sigma_n, sigma_n = SD(delta^A - delta^B) / 2
                   halves    Spearman(delta^A, delta^B) >= 0.50
                   lineages  the median pairwise Spearman of the delta_k >= 0.40
                   library   test only: Spearman of delta against the same rule with E10.5's RNA from libraries 1,
                             2 and 3 >= 0.70
A profile passes when it holds three categories with at least 200 cells in every group and every check holds.
"""

from __future__ import annotations

import hashlib
import json
import math
import zlib
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import rankdata

from vec import t1_scope, t2_build
from vec.boards import REPO_ROOT, load_boards
from vec.core import counts_of, mean_cp10k
from vec.prep import TARGET_SUM, load_stage

AXIS = "heart-xpc-profile"
HEART_BOARD = "T2:heart:test_extrap"  # the panel's genes and their order (the board's own, as the builder uses them)
CATEGORIES = ("cm", "endocardium", "progenitor")
STAGES = ("E8.5", "E9.5", "E10.5")
MERFISH, RNA = "merfish", "rna"
INSTRUMENTS = (MERFISH, RNA)
FILES = {MERFISH: {"E8.5": "T2/E8.5.h5ad", "E9.5": "T2/E9.5.h5ad", "E10.5": "T2/E10.5.h5ad"},
         RNA: {"E8.5": "T1/E8.5_RNA.h5ad", "E9.5": "T1/E9.5_RNA.h5ad", "E10.5": "T1/E10.5_RNA.h5ad"}}
RNA_STEMS = {"E8.5": "E8.5_RNA", "E9.5": "E9.5_RNA", "E10.5": "E10.5_RNA"}  # t1_scope's per-file override keys
DISSECTED = "3"  # E10.5's dissected-heart library (cell names end in _3)
ALL_E105_LIBRARIES = ("1", "2", "3")  # the library check's alternative
MIN_CELLS = 200
SHARE_EPS = 1e-5
FLOOR_CP10K = 0.5
CLIP = math.log(3.0)
HALF_SEED = 20261125
NOISE_RATIO = 2.0
HALF_RHO = 0.50
LINEAGE_RHO = 0.40
LIBRARY_RHO = 0.70
PRACTICE, TEST = "practice", "test"
KINDS = (PRACTICE, TEST)
REFERENCE = {PRACTICE: "E9.5", TEST: "E10.5"}
WEIGHTS = {PRACTICE: {"E9.5": 1.0, "E8.5": -1.0}, TEST: {"E10.5": 1.0, "E8.5": -0.5, "E9.5": -0.5}}
KIND_STAGES = {PRACTICE: ("E8.5", "E9.5"), TEST: STAGES}
PROFILE = REPO_ROOT / "experiments" / "heart_profile" / "results" / "xpc_profile.json"  # what the builder reads
RULE = ("XPC profile: categories cm, endocardium, progenitor (T2 and T1 chamber configs, no stand-ins); M and R = "
        "log(mean panel share + 1e-5) over MERFISH and RNA cells (E10.5 RNA library 3; E8.5 and E9.5 both libraries); "
        "D = M - R; delta_k = D_E10.5 - (D_E8.5 + D_E9.5) / 2 (practice: D_E9.5 - D_E8.5); delta = median over k, 0 "
        "below 0.5 CP10k on either instrument at the reference stage, clipped to |delta| <= log 3; entry: 200 cells "
        "per group; checks: SD(delta) >= 2 sigma_n (halves at seed 20261125), Spearman of the halves >= 0.50, median "
        "pairwise Spearman of the delta_k >= 0.40, test only Spearman against E10.5 RNA libraries 1 to 3 >= 0.70")

Group = tuple[str, str, str]  # (stage, instrument, category)


# ------------------------------------------------------------------------------------------------ reading the files
def panel_genes() -> list[str]:
    """The heart panel in the heart board's own order (the builder's genes)."""
    return list(load_boards()[HEART_BOARD].genes)


def libraries_of(names) -> np.ndarray:
    """Each RNA cell's library, the last `_`-separated part of its name; refused unless every one is a number."""
    lib = np.array([str(x).rsplit("_", 1)[-1] if "_" in str(x) else "" for x in names])
    if not all(k.isdigit() for k in set(lib.tolist())):
        raise ValueError("RNA cell names must end in _<library number>")
    return lib


def read_files(data_dir: Path, stages: Sequence[str], genes: Sequence[str],
               load: Callable = load_stage) -> dict:
    """The released files of `stages`, as raw groups: per instrument, stage and category the cells' panel counts (CSR
    int64, cells x genes in `genes`' order), with each RNA cell's library. Every cell of a category is kept here (the
    cells with no panel count leave in `groups`)."""
    genes = list(genes)
    spec2 = t2_build.load_t2_chamber_spec()
    out: dict = {MERFISH: {}, RNA: {}}
    for s in stages:
        m = load(Path(data_dir) / FILES[MERFISH][s], genes=genes, drop_duplicates=True)
        cat = t2_build.t2_chamber_categories(m.obs[t2_build.T2_LABEL_KEY].to_numpy(), spec2)
        C = counts_of(m, genes)
        out[MERFISH][s] = {k: (C[np.flatnonzero(cat == k)], None) for k in CATEGORIES}
        r = load(Path(data_dir) / FILES[RNA][s], lineage=False)
        cat = np.asarray(t1_scope.classify(r, source=RNA_STEMS[s], key_added=None)).astype(str)
        lib = libraries_of(r.obs_names)
        C = counts_of(r, genes)
        out[RNA][s] = {k: (C[np.flatnonzero(cat == k)], lib[cat == k]) for k in CATEGORIES}
        del m, r
    return out


def groups(raw: dict, stages: Sequence[str],
           e105_libraries: Sequence[str] = (DISSECTED,)) -> dict[Group, sparse.csr_matrix]:
    """The groups a profile reads: per stage, instrument and category the panel counts of the cells that enter (E10.5's
    RNA from `e105_libraries` only; cells with no panel count left out)."""
    out = {}
    for s in stages:
        for inst in INSTRUMENTS:
            for k in CATEGORIES:
                C, lib = raw[inst][s][k]
                rows = np.ones(C.shape[0], dtype=bool)
                if inst == RNA and s == "E10.5":
                    rows &= np.isin(lib, list(e105_libraries))
                rows &= np.asarray(C.sum(1)).ravel() > 0
                out[(s, inst, k)] = sparse.csr_matrix(C[np.flatnonzero(rows)])
    return out


def cells(gs: Mapping[Group, sparse.csr_matrix]) -> dict[str, int]:
    """Cells per group, keyed "stage/instrument/category" (released files' cells; no test board's count)."""
    return {"/".join(g): int(C.shape[0]) for g, C in sorted(gs.items())}


def entry(gs: Mapping[Group, sparse.csr_matrix], stages: Sequence[str] = STAGES) -> dict:
    """The entry: every category at least MIN_CELLS cells on both instruments at every stage."""
    short = sorted("/".join(g) for g in ((s, i, k) for s in stages for i in INSTRUMENTS for k in CATEGORIES)
                   if g not in gs or gs[g].shape[0] < MIN_CELLS)
    return {"ok": not short, "short": short, "min_cells": MIN_CELLS, "cells": cells(gs)}


# ------------------------------------------------------------------------------------------------ the profile
def log_level(C) -> np.ndarray:
    """log(mean over the cells of each gene's panel share + SHARE_EPS)."""
    return np.log(mean_cp10k(C) / TARGET_SUM + SHARE_EPS)


def reference_cp10k(gs: Mapping[Group, sparse.csr_matrix], kind: str) -> dict[str, np.ndarray]:
    """Each instrument's mean CP10k per gene at the kind's reference stage, over the cells of the three categories
    together."""
    s = REFERENCE[kind]
    return {inst: mean_cp10k(sparse.vstack([gs[(s, inst, k)] for k in CATEGORIES]).tocsr()) for inst in INSTRUMENTS}


def eligible_genes(gs: Mapping[Group, sparse.csr_matrix], kind: str) -> np.ndarray:
    """Genes at or above FLOOR_CP10K mean CP10k on both instruments at the kind's reference stage, over the cells of
    the three categories together."""
    lv = reference_cp10k(gs, kind)
    return (lv[MERFISH] >= FLOOR_CP10K) & (lv[RNA] >= FLOOR_CP10K)


def reference_level(gs: Mapping[Group, sparse.csr_matrix], kind: str) -> np.ndarray:
    """log(mean panel share + 1e-5) on MERFISH at the reference stage over the three categories: the level the
    placebos of row H36 are binned on (as row H28 bins its derangement on the reference level)."""
    return np.log(reference_cp10k(gs, kind)[MERFISH] / TARGET_SUM + SHARE_EPS)


def category_profiles(gs: Mapping[Group, sparse.csr_matrix], kind: str) -> dict[str, np.ndarray]:
    """delta_k: the weighted stage differences of D = M - R, per category."""
    out = {}
    for k in CATEGORIES:
        d = sum(w * (log_level(gs[(s, MERFISH, k)]) - log_level(gs[(s, RNA, k)])) for s, w in WEIGHTS[kind].items())
        out[k] = np.asarray(d, dtype=np.float64)
    return out


def combine(by_category: Mapping[str, np.ndarray], eligible: np.ndarray) -> np.ndarray:
    """delta: the median over the categories, 0 off the eligible genes, clipped to |delta| <= CLIP."""
    d = np.median(np.vstack([by_category[k] for k in CATEGORIES]), axis=0)
    d = np.where(eligible, d, 0.0)
    return np.clip(d, -CLIP, CLIP)


def profile(gs: Mapping[Group, sparse.csr_matrix], kind: str, eligible: np.ndarray | None = None) -> dict:
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    e = eligible_genes(gs, kind) if eligible is None else np.asarray(eligible, dtype=bool)
    by_k = category_profiles(gs, kind)
    return {"delta": combine(by_k, e), "by_category": by_k, "eligible": e}


def halves(gs: Mapping[Group, sparse.csr_matrix], seed: int = HALF_SEED) -> tuple[dict, dict]:
    """Every group split in two at `seed`, one stream per group (its name), so that a group's split does not depend
    on any other group."""
    a, b = {}, {}
    for g in sorted(gs):
        n = gs[g].shape[0]
        p = np.random.default_rng([seed, zlib.crc32("/".join(g).encode())]).permutation(n)
        a[g] = gs[g][np.sort(p[: n // 2])]
        b[g] = gs[g][np.sort(p[n // 2:])]
    return a, b


def spearman(x, y) -> float:
    """Spearman's rank correlation (average ranks); nan when either side is constant or has under two values."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if x.size < 2 or x.size != y.size:
        return float("nan")
    rx, ry = rankdata(x), rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _sd(x) -> float:
    x = np.asarray(x, dtype=np.float64)
    return float(x.std(ddof=1)) if x.size > 1 else float("nan")


def measure(gs: Mapping[Group, sparse.csr_matrix], kind: str, alternative=None) -> dict:
    """The profile and its checks. `alternative`: for the test profile, the groups with E10.5's RNA from libraries 1,
    2 and 3 (the library check). A check that cannot be formed (nan) does not hold."""
    ent = entry(gs, KIND_STAGES[kind])
    if not ent["ok"]:
        return {"kind": kind, "entry": ent, "passed": False, "failed": [f"entry: {', '.join(ent['short'])}"]}
    full = profile(gs, kind)
    e = full["eligible"]
    ga, gb = halves(gs)
    da, db = profile(ga, kind, e)["delta"], profile(gb, kind, e)["delta"]
    d = full["delta"]
    sd, sigma = _sd(d[e]), _sd((da - db)[e]) / 2.0
    rho_half = spearman(da[e], db[e])
    pairs = [(a, b) for i, a in enumerate(CATEGORIES) for b in CATEGORIES[i + 1:]]
    rho_k = {f"{a}~{b}": spearman(full["by_category"][a][e], full["by_category"][b][e]) for a, b in pairs}
    rho_lineage = float(np.median(list(rho_k.values()))) if not any(math.isnan(v) for v in rho_k.values()) \
        else float("nan")
    checks = {"noise": {"sd": sd, "sigma_n": sigma, "ratio": NOISE_RATIO, "holds": bool(sd >= NOISE_RATIO * sigma)},
              "halves": {"spearman": rho_half, "bar": HALF_RHO, "holds": bool(rho_half >= HALF_RHO)},
              "lineages": {"pairs": rho_k, "median": rho_lineage, "bar": LINEAGE_RHO,
                           "holds": bool(rho_lineage >= LINEAGE_RHO)}}
    if kind == TEST:
        if alternative is None:
            raise ValueError("the test profile's library check needs the groups with E10.5's RNA from libraries 1-3")
        alt = profile(alternative, kind, e)["delta"]
        rho_lib = spearman(d[e], alt[e])
        checks["library"] = {"spearman": rho_lib, "bar": LIBRARY_RHO, "holds": bool(rho_lib >= LIBRARY_RHO),
                             "alternative_cells": cells(alternative)}
    failed = [name for name, c in checks.items() if not c["holds"]]
    return {"kind": kind, "entry": ent, "delta": d, "by_category": full["by_category"], "eligible": e,
            "level": reference_level(gs, kind), "halves": {"seed": HALF_SEED, "a": da, "b": db}, "checks": checks,
            "failed": failed, "passed": not failed}


# ------------------------------------------------------------------------------------------------ the record
def _floats(v) -> list[float]:
    return [float(x) for x in np.asarray(v, dtype=np.float64)]


def record(res: dict, genes: Sequence[str], files: Mapping | None = None) -> dict:
    """The JSON record of a profile: what the builder reads (kind, passed, genes, delta) and how it was made; no cell
    number of any test board (the cells are the released files')."""
    out = {"axis": AXIS, "kind": res["kind"], "passed": bool(res["passed"]), "failed": list(res["failed"]),
           "rule": RULE, "genes": list(genes), "categories": list(CATEGORIES), "reference": REFERENCE[res["kind"]],
           "weights": WEIGHTS[res["kind"]], "files": dict(files or {}), "entry": res["entry"],
           "constants": {"min_cells": MIN_CELLS, "share_eps": SHARE_EPS, "floor_cp10k": FLOOR_CP10K, "clip": CLIP,
                         "half_seed": HALF_SEED, "noise_ratio": NOISE_RATIO, "half_rho": HALF_RHO,
                         "lineage_rho": LINEAGE_RHO, "library_rho": LIBRARY_RHO, "dissected_library": DISSECTED}}
    if "delta" in res:
        out.update(delta=_floats(res["delta"]), eligible=[bool(x) for x in res["eligible"]],
                   level=_floats(res["level"]), by_category={k: _floats(v) for k, v in res["by_category"].items()},
                   halves={"seed": res["halves"]["seed"], "a": _floats(res["halves"]["a"]),
                           "b": _floats(res["halves"]["b"])},
                   checks=res["checks"])
    return out


def dumps(rec: dict) -> str:
    return json.dumps(rec, indent=1, sort_keys=True, default=float) + "\n"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_profile(path: Path, sha256: str, genes: Sequence[str], kind: str = TEST) -> np.ndarray:
    """The profile a step applies: the file at `path` must have sha256 `sha256`, be of `kind`, have passed, name
    `genes` in their order and hold a finite delta within |delta| <= CLIP. Returns delta (raises ValueError
    otherwise)."""
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != sha256:
        raise ValueError(f"{path} is not the file of sha256 {sha256}")
    rec = json.loads(raw)
    if rec.get("axis") != AXIS or rec.get("kind") != kind or rec.get("passed") is not True:
        raise ValueError(f"{path} is not a passed {kind} profile of {AXIS}")
    if rec.get("genes") != list(genes):
        raise ValueError(f"{path} does not name the panel's genes in the board's order")
    d = np.asarray(rec.get("delta"), dtype=np.float64)
    if d.shape != (len(genes),) or not np.all(np.isfinite(d)) or np.abs(d).max(initial=0.0) > CLIP * (1 + 1e-12):
        raise ValueError(f"{path}: delta is not finite and within log 3 on every panel gene")
    return d
