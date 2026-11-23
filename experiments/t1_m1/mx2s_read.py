#!/usr/bin/env python
"""Try T72 (experiments/harness/TRIES.md): T1 candidate MX2s, a within-study sign-agreement consensus of three heart
atlases (register rows X002, X019, X020) for the CM, Endocardium and Endothelium directions of the T1 direction step,
at the ruled step's size; read on row T62's five practice tasks at the test's fraction 2/5. Practice read only (ruling
row R27). Task 1 only: nothing here enters vec.core, T2 or T3 code (E8.5 lies inside the Task 2 heart window).

The row governs; this script implements its cells 3 and 8 and chooses no constant and no clause. Where the row leaves a
point open, the reading that refuses more is taken and named in the module's OPEN list.

Members (a member is an atlas in one group's consensus): CM has X002, X019 and X020; Endocardium and Endothelium have
X002 and X019; Epicardium keeps X002's direction byte for byte. X002 is row N54's recomputed direction
(`recomputed_read.recomputed_inputs`, i.e. `feng_dirs.feng_directions(genes, groups="recomputed")`); X019 and X020 are
read only from the data side's tables (engineer_notes section B):
  data/derived/t1_ext_t72/tables/<ACC>/<group>.npz   genes (T1 order), m_early, m_late, D, n_cells_early,
                                                     n_cells_late (and per-sample means, not read here)
  data/derived/t1_ext_t72/tables/MANIFEST.json       see `manifest_problems` for every field this script reads

  pregate   Target-free: no T1 cell is loaded or scored and no truth is read. Member clauses (m1) to (m4) drop a member
            from that group's consensus only; then, per changed group (a group left with at least two members), the
            sign-rule consensus (`sign_consensus`), the rescale to X002's L2 size over S_g (`rescale`), the two
            controls (`mask_control`, `perm_member`), and the line clauses (p1) to (p6) and the record. Writes the
            consensus and control tables to data/derived/t1_ext_t72/consensus/ and pregate.json (with its own sha256
            beside it); prints only the fixed words; exit 3 on FAIL. It refuses before it reads anything if pregate.json
            or the consensus folder exists, and creates both exclusively; after its HEAD and source guards it refuses
            (non-zero, never 3, no pre-gate words) on an unfilled or wrong atlases.yaml constant before it creates
            anything.
  run-cut   One E10.5 cut (chamber or whole); run-lib one library task (1, 2 or 3). The base m1rc_gmom_x2 (row T63's
            doubled step: the practice fraction 0.2 times 2.0) and each counted arm by the splice on the same draw;
            every cell outside the changed groups is checked byte for byte per seed and arm, or the run stops and
            writes no file. After the pre-load guards pass and before anything loads, an exclusive start record
            <task>.started. Prints no score.
  check-pregate
            The launcher's `runs` entry's check (no output beside an exit code and a refusal's reason): pregate.json
            passed at --commit and the arm tables' sha256 equal pregate.json's.
  read      Row T72's statement after the gate against row T63's committed per-seed m1rc_gmom_x2 within 1e-9 and the
            void list; prints only the fixed words; writes t72.txt (with a line beginning OUTCOME) and t72.json. It
            refuses, writing nothing, when pregate.json is missing or did not pass, or when no start record exists
            (no try was counted, so no read may say the tries count).

Usage: through experiments/t1_m1/mx2s.sh, from the root of a checkout that holds data/derived/feng_relabel_c/, the
released T1 stages and the row's data folders, at the row's commit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
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
import cmjax_read as cj  # noqa: E402  (T71: its splice, byte record and Spearman)
import recomputed_read as rc  # noqa: E402  (T62: the recomputed direction, its hook, guards and checks)

br, gc = rc.br, rc.gc
bt, gv, pt, rr = rc.bt, rc.gv, rc.pt, rc.rr
feng_dirs = rc.feng_dirs

ROW = "T72"
OUT = "t72_mx2s"
CUTS, CUT_FILES, LIBRARIES = rc.CUTS, rc.CUT_FILES, rc.LIBRARIES
T63_TASKS = RESULTS / "t63_double_step" / "tasks"
OUT_FILES = ("t72.json", "t72.txt")
PREGATE_FILE = "pregate.json"
PREGATE_SHA_FILE = "pregate.json.sha256"
STARTED_SUFFIX = ".started"

BASE = "m1rc_gmom_x2"  # row T63's doubled step: the base and the gate (not counted)
CONS, MASK, PERM = "mx_cons", "mx_mask", "mx_perm"
COUNTED = (CONS, MASK, PERM)
ARMS = (BASE, *COUNTED)
LIB_REFERENCE = rc.LIB_REFERENCE  # copy_last, scored first
LIB_ARMS = (LIB_REFERENCE, *ARMS)
TRIES_COUNT = len(COUNTED) * (len(CUTS) + len(LIBRARIES))  # 15

MULTIPLIER = 2.0
PRACTICE_FRACTION = rc.PRACTICE_FRACTION  # 0.2
TEST_FRACTION = 0.4
FRACTION_TOL = 1e-12
REPRO_TOL = 1e-9  # the gate against T63's committed m1rc_gmom_x2 (not br.REPRO_TOL)
MASK_SEED = 20261123  # plus the group index
PERM_SEED = 20261124  # plus 10 x the member index plus the group index
BINS = 10

MEMBERS = ("X002", "X019", "X020")  # the member index is the position here
NEW_MEMBERS = ("X019", "X020")
ACCESSION = {"X019": "GSE244783", "X020": "GSE230531"}
GROUP_MEMBERS = {"CM": ("X002", "X019", "X020"), "Endocardium": ("X002", "X019"), "Endothelium": ("X002", "X019")}
CONSENSUS_GROUPS = tuple(GROUP_MEMBERS)  # Epicardium keeps X002's direction byte for byte
X002_STAGES = feng_dirs.BRACKET  # ("E9.5", "E14.5")

# The pre-gate's bars (fixed in the row)
MIN_CELLS = {"CM": 200, "Endocardium": 50, "Endothelium": 50}
SKM_SCORE = 0.3
MAX_SKM_SHARE = 0.10
SPEARMAN_LOW, SPEARMAN_HIGH = 0.2, 0.95
X002_MIN_CP = 1.0
ZERO_LOW, ZERO_HIGH = 0.05, 0.60
MAX_CONS_SPEARMAN = 0.98
MAX_HELD_SHARE = 0.5
PREGATE_FAIL = 3
# The data side's constants, checked in its manifest (fixed in the row)
DIRECTION_SCALE = 5.0 / 8.0
PSEUDOCOUNT = 0.5
MIN_TABLE_CELLS = 10
CELL_CALLS = {"GSE230531": {"min_umi": 1000, "min_genes": 500, "max_mito": 0.20},
              "GSE244783": {"min_umi": 500, "min_genes": 300, "max_mito": 0.05}}
# The statement's bars (fixed in the row)
CHAMBER_BAR = 0.3
LIB_BAR = 0.3
LIB_FLOOR = -0.3
MIN_AHEAD = 8
LIB3 = "3"
LIB3_MIN_AHEAD = 7
CTRL_CHAMBER_BAR = 0.2
CTRL_LIB_BAR = 0.2

# Samples (cell 3): the pins and the never-fetched list, by name
PINNED = {"GSE244783": {"early": ("GSM7829699",), "late": ("GSM7829703",)},
          "GSE230531": {"early": ("GSM7226268", "GSM7226269"), "late": ("GSM7226274", "GSM7226276")}}
NEVER_FETCHED = ("GSM7829701", "GSM7829700", "GSM7829702", "GSM7829704", "GSE244783_RAW.tar",
                 "GSM7226279", "GSM7226280", "GSM7226270", "GSM7226271", "GSM7226272", "GSM7226273",
                 "GSE230531_RAW.tar")
# FILLED AT PLACEMENT (cell 3): the frozen atlases.yaml's sha256. Until it is filled (or while the file differs from
# it) the pre-gate, the runs and the runs entry's check refuse before they create anything; at the read it is a void
# reason (cell 8), reached only after a passed pre-gate and a start record.
ATLASES_YAML_SHA256: str | None = None
ATLASES_YAML = ROOT / "experiments" / "t1_ext" / "atlases.yaml"

DERIVED = ROOT / "data" / "derived" / "t1_ext_t72"
TABLES = DERIVED / "tables"
TABLES_MANIFEST = TABLES / "MANIFEST.json"
CONSENSUS = DERIVED / "consensus"

BANNED_IMPORTS = ("perm_read", "eb_read", "step_read", "deep_read", "over_read", "lib_step_read", "ctl_range_read",
                  "double_step_read", "heldmom_read", "jaxblood_read", "n12_hygiene_size")
BANNED_PREFIXES = ("n4" + "8", "n5" + "4")

LINE = "'multi-atlas sign consensus'"
FOLLOWS_PREGATE = ("What follows: nothing is scored and no try is counted under row T72; the line "
                   f"{LINE} closes on T1 at the pre-gate.")
AHEAD_WORDS = "Outcome: MX2s is AHEAD (clauses (i), (ii) and (iii) hold at the test's fraction 2/5)"
FOLLOWS_AHEAD = ("What follows: this read names no file and builds nothing; one priced question goes to the captain, "
                 "default no change; a file would need a builder option, the register amendments and the captain's "
                 "ruling, none of which this row opens.")
FOLLOWS_END = (f"What follows: NOT AHEAD closes {LINE} on T1. No question goes to the captain, nothing is built, "
               "and nothing is rerun under another row.")
FOLLOWS_VOID = ("What follows: the tries count; nothing is built; no question goes to the captain; the line "
                f"{LINE} closes on T1, and nothing is rerun under this row or another.")
BESIDE = ("J2: X002's E9.5, and the new atlases' E8.5 under a wide plug convention, sit in the practice window; no "
          "forecast weight.",
          "J5: a labelled in-sample diagnostic on one truth, E10.5.",
          "J6: the test spans two days and every practice task one (horizon ratio 2), so any gain counts 0 at the "
          "centre and is shown as 0 to the practice value.",
          "The base is R13's file's practice analogue; any step the file in force adds is unpriced (J1).")
IDLE_WORDS = "the arm moved no cell on {task}: nothing is priced there"
LABEL = ("try T72, a labelled in-sample diagnostic, never an adoption result: T1 candidate MX2s, the sign-agreement "
         "consensus of X002, X019 and X020 for CM, Endocardium and Endothelium at X002's L2 size, with its mask and "
         "member-permutation controls, against m1rc_gmom_x2 (row T63's) at the test's fraction 2/5 on T62's five "
         "tasks; every task scores one truth (E10.5) over one day")

# Points the row leaves open, and the reading taken (the one that refuses more)
OPEN = (
    "a member with no table (fewer than 10 cells at an end) is dropped under (m1); a table the manifest names that is "
    "missing or off its sha256 fails the record",
    "a control whose rescale factor does not exist fails (p6) for that arm",
    "a task file without its start record is VOID, as a start record without its task file",
    "every launcher step runs once: its exit record blocks a second call (the runs entry is the one restart)",
    "the idle-task words are printed for every counted arm, each led by the arm's name",
    "the read refuses, writing nothing, when pregate.json is missing or did not pass or no start record exists (no "
    "read may say the tries count for a row with no counted try)",
    "the pre-gate, the runs and the runs entry's check refuse (non-zero, never 3) on an unfilled or wrong atlases.yaml "
    "constant before they create anything",
)


# ---------------------------------------------------------------------------------------------------------------------
# Guards


def source_guard(path: Path | None = None) -> list[str]:
    """recomputed_read's guard on this file (every feng_directions call names the recomputed groups; no token of the
    authors' groups, their statistics or row N54's listing), no import of another row's read (BANNED_IMPORTS) or of
    any row-N48 or row-N54 script, and no dynamic import."""
    path = Path(__file__ if path is None else path)
    bad = rc.source_guard(path)
    text = path.read_text(encoding="utf-8")
    for m in re.finditer(r"^\s*(?:from\s+(\S+)\s+import\s+([^\n#]+)|import\s+([^\n#]+))", text, re.M):
        names = [x for x in re.split(r"[\s,.()]+", " ".join(x for x in m.groups() if x)) if x]
        hit = sorted(set(names) & set(BANNED_IMPORTS)) + sorted(x for x in names if x.startswith(BANNED_PREFIXES))
        if hit:
            bad.append(f"the source imports {hit}")
    for token in ("import" + "lib", "__imp" + "ort__", "spec_from" + "_file_location", "T1_s2" + ".h5ad",
                  "n54" + "_relabel_listing"):
        if token in text:
            bad.append(f"the source names {token}")
    return bad


def all_source_bad() -> list[str]:
    """This script's guard and recomputed_read.py's own."""
    return [*source_guard(), *(f"recomputed_read.py: {x}" for x in rc.source_guard())]


def has_row(path: Path | None = None) -> bool:
    """Exactly one ledger row T72."""
    text = Path(TRIES if path is None else path).read_text(encoding="utf-8")
    return len(re.findall(rf"^\| {ROW} ", text, re.M)) == 1


def need_row() -> None:
    if not has_row():
        sys.exit(f"{ROW} is refused: experiments/harness/TRIES.md does not hold exactly one row {ROW} (a row comes "
                 "before its run)")


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


def sha256_file(path: Path) -> str:
    return rc.sha256_file(Path(path))


def vector_sha(v: np.ndarray) -> str:
    return cj.vector_sha(v)


def yaml_problems(yaml_path: Path | None = None) -> list[str]:
    """The frozen atlases.yaml is the row's: its constant is filled and the file hashes to it."""
    yaml_path = ATLASES_YAML if yaml_path is None else Path(yaml_path)
    if not ATLASES_YAML_SHA256:
        return ["the frozen atlases.yaml's sha256 is not filled in this script (a placement fill)"]
    got = sha256_file(yaml_path) if yaml_path.exists() else None
    return [] if got == ATLASES_YAML_SHA256 else [f"atlases.yaml's sha256 is {got}, not the row's"]


# ---------------------------------------------------------------------------------------------------------------------
# The consensus, the rescale and the controls (target-free, on T1's genes)


def sign_consensus(d_x002: np.ndarray, others: Sequence[np.ndarray]) -> np.ndarray:
    """Per gene i: V = the members (X002 first, then `others`) with a finite, non-zero value; n = |V|. NaN where
    D_X002 is NaN. Else, if n >= 2 and one sign is held by at least 2 members and by more than n/2, c_i = (1/n) x the
    sum of that sign's values; otherwise 0."""
    d_x002 = np.asarray(d_x002, dtype=np.float64)
    stack = np.vstack([d_x002, *[np.asarray(o, dtype=np.float64) for o in others]])
    with np.errstate(invalid="ignore"):
        valid = np.isfinite(stack) & (stack != 0.0)
        pos, neg = valid & (stack > 0), valid & (stack < 0)
    n = valid.sum(axis=0)
    n_pos, n_neg = pos.sum(axis=0), neg.sum(axis=0)
    s_pos = np.where(pos, stack, 0.0).sum(axis=0)
    s_neg = np.where(neg, stack, 0.0).sum(axis=0)
    c = np.zeros(stack.shape[1], dtype=np.float64)
    keep_pos = (n >= 2) & (n_pos >= 2) & (n_pos > n / 2.0)
    keep_neg = (n >= 2) & (n_neg >= 2) & (n_neg > n / 2.0)
    c[keep_pos] = s_pos[keep_pos] / n[keep_pos]
    c[keep_neg] = s_neg[keep_neg] / n[keep_neg]
    c[np.isnan(d_x002)] = np.nan
    return c


def s_genes(d_x002: np.ndarray, free: np.ndarray) -> np.ndarray:
    """S_g: T1's genes not on R2's held list at which D_X002 is finite."""
    return np.asarray(free, dtype=bool) & np.isfinite(np.asarray(d_x002, dtype=np.float64))


def rescale(c: np.ndarray, d_x002: np.ndarray, s: np.ndarray) -> tuple[np.ndarray, float | None]:
    """c x k, k = the L2 norm of D_X002 over S_g divided by the L2 norm of c over S_g (the whole vector multiplied).
    A zero or non-finite norm of c gives no factor: (c unchanged, None)."""
    c = np.asarray(c, dtype=np.float64)
    s = np.asarray(s, dtype=bool)
    den = float(np.sqrt(np.sum(np.square(c[s]))))
    num = float(np.sqrt(np.sum(np.square(np.asarray(d_x002, dtype=np.float64)[s]))))
    if not (np.isfinite(den) and den > 0 and np.isfinite(num)):
        return c, None
    k = num / den
    return c * k, k


def zero_share(c: np.ndarray, s: np.ndarray) -> float:
    """The share of S_g at which the consensus is 0 (NaN when S_g is empty)."""
    s = np.asarray(s, dtype=bool)
    return float(np.mean(np.asarray(c)[s] == 0.0)) if s.any() else float("nan")


def held_share(d: np.ndarray, held: np.ndarray) -> float:
    """The share of the squared L2 norm over T1's genes with finite values that falls on R2's held genes."""
    d = np.asarray(d, dtype=np.float64)
    sq = np.where(np.isfinite(d), d, 0.0) ** 2
    total = float(sq.sum())
    return float(sq[np.asarray(held, dtype=bool)].sum() / total) if total > 0 else float("nan")


def bins_of(idx: np.ndarray, level: np.ndarray, bins: int = BINS) -> list[np.ndarray]:
    """`idx` in `bins` equal-count bins of `level` (np.array_split over a stable rank, ties by gene position)."""
    order = idx[np.argsort(np.asarray(level, dtype=np.float64)[idx], kind="stable")]
    return np.array_split(order, bins)


def mask_seed(group_index: int) -> int:
    """The mask's generator seed: 20261123 plus the group index."""
    return MASK_SEED + int(group_index)


def perm_seed(member_index: int, group_index: int) -> int:
    """A permutation's generator seed: 20261124 plus 10 x the member index plus the group index."""
    return PERM_SEED + 10 * int(member_index) + int(group_index)


def mask_control(d_x002: np.ndarray, c: np.ndarray, level_e95: np.ndarray, free: np.ndarray, group_index: int,
                 bins: int = BINS) -> tuple[np.ndarray, list[dict]]:
    """X002's direction with, in each of `bins` equal-count bins of X002's E9.5 group mean CP10k over S'_g (S_g with
    D_X002 non-zero and that level finite), as many genes set to 0 as the consensus set to 0 in that bin, drawn
    uniformly without replacement, bins in order, one generator default_rng(MASK_SEED + group index). Returns the
    masked vector (before the rescale) and per-bin records."""
    d = np.asarray(d_x002, dtype=np.float64)
    c = np.asarray(c, dtype=np.float64)
    level = np.asarray(level_e95, dtype=np.float64)
    s = s_genes(d, free)
    with np.errstate(invalid="ignore"):
        sp = s & (d != 0.0) & np.isfinite(level)
    rng = np.random.default_rng(mask_seed(group_index))
    out = d.copy()
    rec = []
    for b in bins_of(np.flatnonzero(sp), level, bins):
        k_b = int(np.sum(c[b] == 0.0))
        if k_b:
            out[b[rng.choice(b.size, size=k_b, replace=False)]] = 0.0
        rec.append({"genes": int(b.size), "consensus_zeros": k_b, "mask_zeros": int(np.sum(out[b] == 0.0))})
    return out, rec


def perm_member(d: np.ndarray, m_early: np.ndarray, free: np.ndarray, member_index: int, group_index: int,
                bins: int = BINS) -> tuple[np.ndarray, dict]:
    """A new member's D permuted within `bins` equal-count bins of its own early-end mean CP10k over its non-held genes
    with D and that level finite (NaN stays at its gene; held genes untouched), bins in order, one generator
    default_rng(PERM_SEED + 10 x member index + group index). Returns the permuted vector and its record."""
    d = np.asarray(d, dtype=np.float64)
    level = np.asarray(m_early, dtype=np.float64)
    g = np.asarray(free, dtype=bool) & np.isfinite(d) & np.isfinite(level)
    seed = perm_seed(member_index, group_index)
    rng = np.random.default_rng(seed)
    out = d.copy()
    same = True
    for b in bins_of(np.flatnonzero(g), level, bins):
        out[b] = d[b[rng.permutation(b.size)]]
        same &= bool(np.array_equal(np.sort(out[b]), np.sort(d[b])))
    nan_kept = bool(np.array_equal(np.isnan(out), np.isnan(d)))
    untouched = bool(np.array_equal(out[~g], d[~g], equal_nan=True))
    return out, {"seed": seed, "genes": int(g.sum()), "bins_keep_values": same, "nan_kept": nan_kept,
                 "outside_untouched": untouched}


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return cj._spearman(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate (target-free)


def x002_group_means(genes: Sequence[str], group: str) -> dict[str, np.ndarray]:
    """X002's recomputed group mean CP10k at each bracket end on `genes`, through feng_dirs' own checks of the
    recomputed file and its group map."""
    rcfg = feng_dirs.load_relabel_config()
    path = feng_dirs.recomputed_stats_path(rcfg)
    feng_dirs.check_recomputed(path)
    st = feng_dirs.dl.load_stats(path)
    feng_dirs.check_recomputed_keys(st)
    cfg = feng_dirs.groups_config(feng_dirs.load_config(), rcfg)
    labels = cfg["external_labels"][group][0]
    if not labels:
        raise ValueError(f"X002's recomputed groups have no {group} group")
    idx = feng_dirs.dl.sym_index(st["genes"], list(genes))
    out = {}
    for stage in X002_STAGES:
        sel = feng_dirs.dl.feng_sel(st, stage, labels)
        if sel is None:
            raise ValueError(f"X002's recomputed {group} group has no cell at {stage}")
        out[stage] = feng_dirs.dl.on_t1(sel["cp"], idx)
    return out


def x002_counts(record: Path | None = None) -> dict:
    """X002's labelled cells per group and end, from row N54's committed record."""
    rec = json.loads(Path(rc.RECORD if record is None else record).read_text())
    return {g: {s: int(v) for s, v in d.items()} for g, d in (rec.get("group_cells") or {}).items()}


def data_constants_ok(c: Mapping | None) -> bool:
    """The data side's constants as its manifest records them are the row's: the scale 5/8, the pseudocount 0.5,
    at least 10 cells at each end for a table, the cell calls, and the skm score bar of the purity share."""
    c = c or {}
    return (c.get("direction_scale"), c.get("pseudocount"), c.get("min_cells"), c.get("skm_score"),
            json.loads(json.dumps(c.get("cell_calls")))) == (DIRECTION_SCALE, PSEUDOCOUNT, MIN_TABLE_CELLS, SKM_SCORE,
                                                             json.loads(json.dumps(CELL_CALLS)))


def manifest_problems(man: Mapping, commit: str, genes_ok: bool = True) -> tuple[list[str], list[str]]:
    """The tables manifest against the row: returns (p1 failures, record failures). Fields read:
    tasks, row, commit, clean, yaml_sha256, author_labels_read (no author label, cluster or embedding),
    constants {direction_scale, pseudocount, min_cells, skm_score, cell_calls {ACC: {min_umi, min_genes,
    max_mito}}},
    atlases {ACC: {register, status ('use' or 'dropped'), samples {GSM: {end, pinned, loaded}},
    tables {group: {file, sha256, n_cells_early, n_cells_late}}, skm_late_share (CM)}}."""
    p1, record = [], []
    if list(man.get("tasks") or []) != ["T1"]:
        record.append("the tables manifest does not say tasks: [T1]")
    if man.get("row") != ROW:
        record.append(f"the tables were not written under row {ROW}")
    if man.get("commit") != commit or man.get("clean") is not True:
        record.append("the tables were not written on a clean tree at this row's commit")
    if not ATLASES_YAML_SHA256 or man.get("yaml_sha256") != ATLASES_YAML_SHA256:
        record.append("the tables were written under another atlases.yaml than the row's")
    if man.get("author_labels_read") is not False:
        record.append("the record does not say that no author label, cluster or embedding was read")
    if not data_constants_ok(man.get("constants")):
        record.append("a constant of the tables is not the row's")
    atlases = man.get("atlases") or {}
    if set(atlases) != set(PINNED):
        record.append(f"the manifest's atlases {sorted(atlases)} are not the row's {sorted(PINNED)}")
    for acc, a in atlases.items():
        if acc not in PINNED:
            continue
        if a.get("register") != {v: k for k, v in ACCESSION.items()}[acc]:
            record.append(f"{acc}: the register row is not the row's")
        status = a.get("status")
        if status not in ("use", "dropped"):
            record.append(f"{acc}: status {status!r} is neither use nor dropped")
        samples = a.get("samples") or {}
        pins = {g: e for e, gs in PINNED[acc].items() for g in gs}
        for gsm, s in samples.items():
            if gsm in NEVER_FETCHED:
                record.append(f"{acc}: a never-fetched accession entered: {gsm}")
            elif gsm not in pins or s.get("end") != pins[gsm] or s.get("pinned") is not True:
                record.append(f"{acc}: an unpinned sample entered: {gsm}")
        if status == "use":
            missing = sorted(g for g in pins if (samples.get(g) or {}).get("loaded") is not True)
            if missing:
                p1.append(f"(p1) {acc}: pinned samples not loaded under their gates: {missing}")
        for g, t in (a.get("tables") or {}).items():
            if not t.get("sha256") or not t.get("file"):
                record.append(f"{acc} {g}: a table is named without its file or sha256")
    if not genes_ok:
        record.append("a table's genes are not T1's genes in T1's order")
    return p1, record


def read_tables(man: Mapping, genes: Sequence[str], tables: Path | None = None) -> tuple[dict, list[str], bool]:
    """The member tables the manifest names for the row's members and groups, each checked against its sha256.
    Returns {(register, group): arrays}, record failures and whether every table's genes are T1's."""
    tables = TABLES if tables is None else Path(tables)
    out, bad, genes_ok = {}, [], True
    for reg in NEW_MEMBERS:
        acc = ACCESSION[reg]
        a = (man.get("atlases") or {}).get(acc) or {}
        if a.get("status") != "use":
            continue
        for g, t in (a.get("tables") or {}).items():
            if g not in GROUP_MEMBERS or reg not in GROUP_MEMBERS[g]:
                continue
            path = tables / str(t.get("file"))
            if not path.is_file() or sha256_file(path) != t.get("sha256"):
                bad.append(f"{acc} {g}: the table is missing or not the one the manifest names")
                continue
            with np.load(path, allow_pickle=False) as z:
                arr = {k: z[k] for k in ("genes", "m_early", "m_late", "D", "n_cells_early", "n_cells_late")}
            if [str(x) for x in arr["genes"]] != list(genes):
                genes_ok = False
                continue
            out[(reg, g)] = {"D": np.asarray(arr["D"], dtype=np.float64),
                             "m_early": np.asarray(arr["m_early"], dtype=np.float64),
                             "n_early": int(arr["n_cells_early"]), "n_late": int(arr["n_cells_late"]),
                             "file": str(t.get("file")), "sha256": t.get("sha256")}
    return out, bad, genes_ok


def member_clauses(reg: str, g: str, table: Mapping | None, skm_late_share, d_x002: np.ndarray,
                   x002_cp: Mapping[str, np.ndarray], held: np.ndarray) -> dict:
    """(m1) to (m4) for one new member in one group; every figure is recorded, every failed clause is named."""
    failed, fig = [], {}
    if table is None:
        failed.append("(m1)")
        fig["n_early"] = fig["n_late"] = None
        return {"failed": failed, "figures": fig}
    fig["n_early"], fig["n_late"] = table["n_early"], table["n_late"]
    if not (table["n_early"] >= MIN_CELLS[g] and table["n_late"] >= MIN_CELLS[g]):
        failed.append("(m1)")
    if g == "CM":
        share = float(skm_late_share) if isinstance(skm_late_share, (int, float)) else float("nan")
        fig["skm_late_share"] = share
        if not share <= MAX_SKM_SHARE:
            failed.append("(m2)")
    d = table["D"]
    expressed = np.zeros(len(d), dtype=bool)
    for v in x002_cp.values():
        expressed |= np.nan_to_num(np.asarray(v, dtype=np.float64), nan=-1.0) >= X002_MIN_CP
    both = np.isfinite(d) & np.isfinite(d_x002) & expressed
    rho = spearman(d[both], d_x002[both])
    fig["spearman"], fig["spearman_genes"] = rho, int(both.sum())
    if not SPEARMAN_LOW <= rho <= SPEARMAN_HIGH:
        failed.append("(m3)")
    hs = held_share(d, held)
    fig["held_share"] = hs
    if not hs < MAX_HELD_SHARE:
        failed.append("(m4)")
    fig["genes_finite"] = int(np.isfinite(d).sum())
    return {"failed": failed, "figures": fig}


def group_arms(g: str, d_x002: np.ndarray, members: Mapping[str, Mapping], x002_e95: np.ndarray, free: np.ndarray,
               held: np.ndarray) -> dict:
    """One changed group: the consensus, the mask control and the permutation control, each rescaled, with their
    figures and validity records. `members` maps each surviving new member to its table."""
    gi = harness.LP_GROUPS.index(g)
    s = s_genes(d_x002, free)
    names = [m for m in NEW_MEMBERS if m in members]
    c0 = sign_consensus(d_x002, [members[m]["D"] for m in names])
    cons, k_cons = rescale(c0, d_x002, s)
    masked0, mask_bins = mask_control(d_x002, c0, x002_e95, free, gi)
    mask, k_mask = rescale(masked0, d_x002, s)
    perms, perm_rec = [], {}
    for m in names:
        p, r = perm_member(members[m]["D"], members[m]["m_early"], free, MEMBERS.index(m), gi)
        perms.append(p)
        perm_rec[m] = r
    p0 = sign_consensus(d_x002, perms)
    perm, k_perm = rescale(p0, d_x002, s)
    out = {"group_index": gi, "members": ["X002", *names], "s_genes": int(s.sum()), "arms": {}, "mask_bins": mask_bins,
           "mask_seed": mask_seed(gi), "perm": perm_rec, "vectors": {CONS: cons, MASK: mask, PERM: perm}}
    for arm, v, k in ((CONS, cons, k_cons), (MASK, mask, k_mask), (PERM, perm, k_perm)):
        out["arms"][arm] = {"k_g": k, "zero_share": zero_share(v, s), "spearman_x002": spearman(v[s], d_x002[s]),
                            "held_share": held_share(v, held), "vector_sha256": vector_sha(v)}
    out["member_gene_counts"] = {"X002": int(np.isfinite(d_x002).sum()),
                                 **{m: int(np.isfinite(members[m]["D"]).sum()) for m in names}}
    return out


def pregate_reading(genes: Sequence[str], d_x002: Mapping[str, np.ndarray], x002_cp: Mapping[str, Mapping],
                    x002_n: Mapping[str, Mapping[str, int]], held: np.ndarray, man: Mapping, tables: Mapping,
                    table_bad: Sequence[str], genes_ok: bool, commit: str, yaml_bad: Sequence[str] = ()) -> dict:
    """Every clause of the pre-gate in the row's order: the record and (p1); the member clauses and the drops; the
    consensus, the rescale and the controls of each changed group; (p2) to (p6)."""
    held = np.asarray(held, dtype=bool)
    free = ~held
    p1, record = manifest_problems(man, commit, genes_ok)
    record = [*record, *table_bad, *yaml_bad]
    failed = [*p1, *(f"record: {x}" for x in record)]
    # (m1) on X002 (cell 8): its CM end counts, from N54's record, are recorded in pregate.json; X002 is never
    # dropped, so a count below 200 at either end fails the pre-gate
    x002_cm = {e: (x002_n.get("CM") or {}).get(e) for e in X002_STAGES}
    x002_m1_failed = [e for e in X002_STAGES
                      if not (isinstance(x002_cm[e], int) and x002_cm[e] >= MIN_CELLS["CM"])]
    failed += [f"(m1) X002 CM at {e}: below {MIN_CELLS['CM']} labelled cells in N54's record" for e in x002_m1_failed]
    atlases = man.get("atlases") or {}
    members, drops = {}, []
    for g, regs in GROUP_MEMBERS.items():
        members[g] = {}
        for reg in regs[1:]:
            a = atlases.get(ACCESSION[reg]) or {}
            if a.get("status") != "use":
                drops.append({"member": reg, "group": g, "clauses": ["screens"], "figures": {}})
                continue
            r = member_clauses(reg, g, tables.get((reg, g)), a.get("skm_late_share"), d_x002[g], x002_cp[g], held)
            if r["failed"]:
                drops.append({"member": reg, "group": g, "clauses": r["failed"], "figures": r["figures"]})
            else:
                members[g][reg] = {**tables[(reg, g)], "figures": r["figures"]}
    changed = [g for g in CONSENSUS_GROUPS if members[g]]
    groups = {}
    if not changed:
        failed.append("(p2) no changed group")
    for g in changed:
        ga = group_arms(g, d_x002[g], members[g], x002_cp[g][X002_STAGES[0]], free, held)
        groups[g] = ga
        c = ga["arms"][CONS]
        if not ZERO_LOW <= c["zero_share"] <= ZERO_HIGH:
            failed.append(f"(p3) {g}")
        if not c["spearman_x002"] < MAX_CONS_SPEARMAN:
            failed.append(f"(p4) {g}")
        if not c["held_share"] < MAX_HELD_SHARE:
            failed.append(f"(p5) {g}")
        for arm in COUNTED:  # OPEN[2]: a control without a factor fails (p6) for that arm
            if ga["arms"][arm]["k_g"] is None:
                failed.append(f"(p6) {g}" + ("" if arm == CONS else f" {arm}"))
    return {"passes": not failed, "failed": failed, "changed_groups": changed, "drops": drops, "groups": groups,
            "members": {g: {m: v["figures"] for m, v in d.items()} for g, d in members.items()},
            "x002_m1": {"source": "row N54's record (group_cells)", "cm_cells": x002_cm, "bar": MIN_CELLS["CM"],
                        "failed_ends": x002_m1_failed}}


def pregate_words(r: Mapping) -> list[str]:
    if r["passes"]:
        drops = "; ".join(f"{d['member']}, {d['group']}, {' '.join(d['clauses'])}" for d in r["drops"]) or "none"
        return [f"Pre-gate: PASSES (changed groups: {', '.join(r['changed_groups'])}; members dropped: {drops})"]
    return [f"Pre-gate: FAILS ({'; '.join(r['failed'])})", FOLLOWS_PREGATE]


def write_arm_tables(dest: Path, genes: Sequence[str], groups: Mapping) -> dict:
    """consensus/<group>_{cons,mask,perm}.npz, each created exclusively; returns their paths and sha256."""
    out = {}
    for g, ga in groups.items():
        out[g] = {}
        for arm, v in ga["vectors"].items():
            path = dest / f"{g}_{arm.split('_')[1]}.npz"
            with open(path, "xb") as f:
                np.savez(f, genes=np.asarray(list(genes), dtype=str), direction=np.asarray(v, dtype=np.float64))
            out[g][arm] = {"file": path.name, "sha256": sha256_file(path), "vector_sha256": vector_sha(v)}
    return out


def pregate_constants() -> dict:
    return {"min_cells": MIN_CELLS, "skm_score": SKM_SCORE, "max_skm_share": MAX_SKM_SHARE,
            "spearman": [SPEARMAN_LOW, SPEARMAN_HIGH], "x002_min_cp": X002_MIN_CP, "zero_share": [ZERO_LOW, ZERO_HIGH],
            "max_cons_spearman": MAX_CONS_SPEARMAN, "max_held_share": MAX_HELD_SHARE, "mask_seed": MASK_SEED,
            "perm_seed": PERM_SEED, "bins": BINS, "direction_scale": DIRECTION_SCALE, "pseudocount": PSEUDOCOUNT,
            "min_table_cells": MIN_TABLE_CELLS, "cell_calls": CELL_CALLS}


def _jsonable(r: Mapping) -> dict:
    """The reading without the vectors (they are in the consensus tables)."""
    out = dict(r)
    out["groups"] = {g: {k: v for k, v in ga.items() if k != "vectors"} for g, ga in r["groups"].items()}
    return out


def cmd_pregate(a) -> int:
    tasks = RESULTS / a.out / "tasks"
    if (tasks / PREGATE_FILE).exists() or CONSENSUS.exists():
        sys.exit(f"{ROW} is refused: pregate.json or the consensus folder exists (the pre-gate runs once)")
    need_row()
    started = head_is(a.commit)
    bad = all_source_bad()
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    bad = yaml_problems()  # a placement slip stops the step (never exit 3): nothing is created, no pre-gate words
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    tasks.mkdir(parents=True, exist_ok=True)
    CONSENSUS.parent.mkdir(parents=True, exist_ok=True)
    os.mkdir(CONSENSUS)  # exclusive: a second pre-gate is refused from here on
    genes = list(get_board(rr.BOARD).genes)  # the gene names alone
    man = json.loads(TABLES_MANIFEST.read_text())
    tables, table_bad, genes_ok = read_tables(man, genes)
    directions, _ = rc.recomputed_inputs(genes)
    held_m1, _ = rc._held_m1(t1_scope.held_flat_genes())
    held = np.isin(np.asarray(genes, dtype=str), held_m1)
    x002_cp = {g: x002_group_means(genes, g) for g in CONSENSUS_GROUPS}
    reading = pregate_reading(genes, {g: np.asarray(directions[g], dtype=np.float64) for g in CONSENSUS_GROUPS},
                              x002_cp, x002_counts(), held, man, tables, table_bad, genes_ok, a.commit, yaml_problems())
    arm_tables = write_arm_tables(CONSENSUS, genes, reading["groups"])
    res = {"label": LABEL, "row": ROW, "reading": _jsonable(reading), "arm_tables": arm_tables,
           "member_tables": {f"{m} {g}": {"file": t["file"], "sha256": t["sha256"]} for (m, g), t in tables.items()},
           "tables_manifest_sha256": sha256_file(TABLES_MANIFEST), "yaml_sha256": man.get("yaml_sha256"),
           "samples": {acc: sorted(x.get("samples") or {}) for acc, x in (man.get("atlases") or {}).items()},
           "author_labels_read": man.get("author_labels_read"), "data_constants": man.get("constants"),
           "constants": pregate_constants(), "held_m1_genes": len(held_m1), "open": list(OPEN),
           "stamp": gv.stamp(started, [TABLES_MANIFEST, t1_scope.HELD_FLAT_GENES, feng_dirs.recomputed_stats_path(),
                                       rc.RECORD, *(TABLES / t["file"] for t in tables.values())])}
    text = json.dumps(res, indent=1, default=float) + "\n"
    write_new(tasks / PREGATE_FILE, text)
    write_new(tasks / PREGATE_SHA_FILE, hashlib.sha256(text.encode()).hexdigest() + "  " + PREGATE_FILE + "\n")
    print("\n".join(pregate_words(reading)))  # the fixed words only; every figure is in pregate.json
    return 0 if reading["passes"] else PREGATE_FAIL


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def arm_directions(directions: Mapping, entries: Mapping[str, np.ndarray]) -> dict:
    """The direction dict with each changed group's entry replaced (every other entry the same object)."""
    missing = sorted(set(entries) - set(directions))
    if missing:
        raise ValueError(f"the recomputed directions have no entry for {missing}")
    return {g: (np.asarray(entries[g], dtype=np.float64) if g in entries else v) for g, v in directions.items()}


def changed_map(groups: Mapping[str, str], changed: Sequence[str]) -> dict[str, str]:
    """The label -> group map with every changed group's labels mapped to cmjax_read's CM, so that its splice and
    byte record act on the changed groups together; every other label keeps a name no changed group has."""
    return {lab: (cj.CM if g in changed else f"unchanged:{g}") for lab, g in groups.items()}


def splice(base_moved, base_final, arm_moved, hook, groups: Mapping[str, str], changed: Sequence[str]):
    """cmjax_read's splice generalised to the changed groups: the arm's direction-moved rows of those groups in the
    base's moved population, then g_mom (`hook`), then those rows in the base's final population."""
    try:
        return cj.splice(base_moved, base_final, arm_moved, hook, changed_map(groups, changed))
    except SystemExit as e:
        raise SystemExit(f"{ROW}: {e}") from None


def byte_record(labels: np.ndarray, groups: Mapping[str, str], changed: Sequence[str], base_X, arm_X,
                name: str) -> dict:
    """Every row whose label is outside the changed groups has the base's bytes in the arm, or the run stops."""
    try:
        r = cj.cm_record(labels, changed_map(groups, changed), base_X, arm_X, name)
    except SystemExit as e:
        raise SystemExit(f"{ROW}: {str(e).replace('the CM group', 'the changed groups')}") from None
    return {"changed_cells": r["cm_cells"], "cells_changed": r["cm_cells_changed"], "outside_equal": r["non_cm_equal"]}


def cut_arms(ctx: harness.BuildContext, factors, directions: Mapping, fraction: float, held_m1: Sequence[str],
             entries: Mapping[str, Mapping[str, np.ndarray]]) -> tuple[dict[str, Stage], dict]:
    """The base and the three counted arms of one seed on a harness cut, from one draw, each finished as the harness
    finishes a candidate; `entries` maps each counted arm to its changed groups' directions."""
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
    f = MULTIPLIER * float(fraction)

    def build(dirs, post) -> Stage:
        sub, _ = t1_scope.build_submission(ref, ctx.board, shares, n, ctx.seed, ref_source=src,
                                           held_flat=list(held_m1), directions=dirs, direction_fraction=f,
                                           direction_groups=groups, post_cells=post)
        return harness.finish(ctx, ref, dense(sub.X), None)

    cap = {}

    def base_post(pop):
        cap["moved"] = pop
        cap["final"] = hook(pop)
        return cap["final"]

    out = {BASE: build(directions, base_post)}
    drawn_labels = cap["moved"].obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    rec = {}
    for name in COUNTED:
        changed = sorted(entries[name])

        def post(pop, changed=changed):
            return splice(cap["moved"], cap["final"], pop, hook, groups, changed)

        out[name] = build(arm_directions(directions, entries[name]), post)
        rec[name] = byte_record(drawn_labels, groups, changed, out[BASE].X, out[name].X, name)
    return out, rec


def library_arms(e95, half_shares: Mapping[str, float], factors, held: Sequence[str], held_m1: Sequence[str],
                 seed: int, within: Mapping, fb: Mapping[str, Mapping], directions: Mapping, fraction: float,
                 entries: Mapping[str, Mapping[str, np.ndarray]], n_cells: int = rr.N_CELLS) -> tuple[dict, dict]:
    """The base and the three counted arms of one seed on a library task from T60's oracle draw
    (recomputed_read.library_pair's route), as AnnData, and the per-arm records."""
    base0, _ = t1_scope.reweight(e95, half_shares, n_cells, seed, within=within, fallback=fb["E9.5"])
    groups, covered, with_direction = rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), directions)
    hook = rc.rc_hook(factors, seed, held, with_direction, covered)
    labels = base0.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    f = MULTIPLIER * float(fraction)
    moved = gc.moved_by_feng(base0, directions, f, groups, list(held_m1), [seed, 2])
    out = {BASE: hook(moved)}
    rec = {}
    for name in COUNTED:
        changed = sorted(entries[name])
        arm_moved = gc.moved_by_feng(base0, arm_directions(directions, entries[name]), f, groups, list(held_m1),
                                     [seed, 2])
        out[name] = splice(moved, out[BASE], arm_moved, hook, groups, changed)
        rec[name] = byte_record(labels, groups, changed, out[BASE].X, out[name].X, name)
    return out, rec


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def started_name(task_file: str) -> str:
    return Path(task_file).stem + STARTED_SUFFIX


def load_pregate(tasks: Path) -> tuple[dict | None, str | None, list[str]]:
    """pregate.json, its sha256 and what is wrong with it (missing, its own sha256 off)."""
    p, s = tasks / PREGATE_FILE, tasks / PREGATE_SHA_FILE
    if not p.exists():
        return None, None, ["pregate.json is missing"]
    sha = sha256_file(p)
    bad = []
    if not s.exists() or s.read_text().split()[:1] != [sha]:
        bad.append("pregate.json's sha256 is not the one kept beside it")
    return json.loads(p.read_text()), sha, bad


def arm_table_problems(pg: Mapping, consensus: Path | None = None) -> list[str]:
    """Every arm table pregate.json names is present with its sha256."""
    consensus = CONSENSUS if consensus is None else Path(consensus)
    bad = []
    at = pg.get("arm_tables") or {}
    changed = list((pg.get("reading") or {}).get("changed_groups") or [])
    if not changed or sorted(at) != sorted(changed):
        bad.append("pregate.json names no arm tables for its changed groups")
    for g, arms in at.items():
        if sorted(arms) != sorted(COUNTED):
            bad.append(f"{g}: pregate.json does not name the three arms' tables")
        for arm, t in arms.items():
            path = consensus / str(t.get("file"))
            if not path.is_file() or sha256_file(path) != t.get("sha256"):
                bad.append(f"{g} {arm}: the arm table is not the one pregate.json names")
    return bad


def pregate_problems(pg: Mapping | None, commit: str) -> list[str]:
    if pg is None:
        return ["the pre-gate has not run"]
    bad = []
    if (pg.get("reading") or {}).get("passes") is not True:
        bad.append("the pre-gate has not passed (nothing is scored after a failed pre-gate)")
    st = pg.get("stamp") or {}
    if st.get("commit") != commit or st.get("tree_clean") is not True:
        bad.append("the pre-gate ran at another commit or on an unclean tree")
    return bad


def refuse(tasks: Path, task_file: str, commit: str) -> tuple[dict, dict, dict, str]:
    """The refusals before anything loads; returns the start state, the relabelling's checks, pregate.json and its
    sha256."""
    need_row()
    if (tasks / task_file).exists():
        sys.exit(f"{ROW} is refused: {tasks / task_file} exists (nothing is ever overwritten)")
    if (tasks / started_name(task_file)).exists():
        sys.exit(f"{ROW} is refused: this task's start record exists (a task that started is not started again)")
    started = head_is(commit)
    bad = all_source_bad() + yaml_problems()
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    checks = rc.relabel_checks()
    if not checks["passes"]:
        sys.exit(f"{ROW} is refused: {checks['reasons']}")
    pg, pg_sha, bad = load_pregate(tasks)
    bad += pregate_problems(pg, commit) if pg is not None else []
    bad += arm_table_problems(pg) if pg is not None else []
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    check_fractions(feng_dirs.bracket_fraction(9.5, 10.5), "the practice bracket")
    return started, checks, pg, pg_sha


def write_started(tasks: Path, task_file: str, commit: str) -> Path:
    """The exclusive start record, after the pre-load guards and before anything loads."""
    path = tasks / started_name(task_file)
    write_new(path, json.dumps({"row": ROW, "task_file": task_file, "commit": commit,
                                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}) + "\n")
    return path


def arm_entries(pg: Mapping, genes: Sequence[str], consensus: Path | None = None) -> dict[str, dict[str, np.ndarray]]:
    """Each counted arm's changed-group directions on `genes`, from the tables pregate.json names (by gene name)."""
    consensus = CONSENSUS if consensus is None else Path(consensus)
    out = {arm: {} for arm in COUNTED}
    for g, arms in (pg.get("arm_tables") or {}).items():
        for arm, t in arms.items():
            with np.load(consensus / t["file"], allow_pickle=False) as z:
                tg, v = [str(x) for x in z["genes"]], np.asarray(z["direction"], dtype=np.float64)
            if vector_sha(v) != t.get("vector_sha256"):
                raise SystemExit(f"{ROW}: {g} {arm}: the table's direction is not the one pregate.json names")
            pos = {x: i for i, x in enumerate(tg)}
            if len(pos) != len(tg) or set(genes) != set(tg):
                raise SystemExit(f"{ROW}: {g} {arm}: the table's genes are not the task's genes")
            out[arm][g] = v[[pos[x] for x in genes]]
    return out


def constants(held_m1: Sequence[str], pg: Mapping) -> dict:
    return {"multiplier": MULTIPLIER, "practice_fraction": PRACTICE_FRACTION, "test_fraction": TEST_FRACTION,
            **pregate_constants(), "pregate_constants": pg.get("constants"), "data_constants": pg.get("data_constants"),
            "held_m1_genes": len(held_m1), "edit_categories": list(bt.G_MOM_EDIT_CATEGORIES)}


def check_fractions(fraction: float, where: str) -> None:
    if abs(fraction - PRACTICE_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the practice fraction is {fraction}, not {PRACTICE_FRACTION}")
    if abs(feng_dirs.bracket_fraction(10.5, 12.5) - TEST_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the test's fraction is not {TEST_FRACTION}")
    if abs(MULTIPLIER * fraction - TEST_FRACTION) > FRACTION_TOL:
        sys.exit(f"{where}: the applied fraction is not {TEST_FRACTION}")


def task_record(pg: Mapping, pg_sha: str, ext: Mapping, checks: Mapping, held_m1, fraction: float) -> dict:
    return {"label": LABEL, "row": ROW, "feng_fraction": fraction, "multiplier": MULTIPLIER,
            "applied_fraction": MULTIPLIER * fraction, "relabel": rc._relabel_fields(ext, checks),
            "pregate_sha256": pg_sha, "arm_tables": pg.get("arm_tables"), "yaml_sha256": pg.get("yaml_sha256"),
            "changed_groups": (pg.get("reading") or {}).get("changed_groups"), "constants": constants(held_m1, pg)}


def cmd_run_cut(a) -> None:
    tasks = RESULTS / a.out / "tasks"
    name = CUT_FILES[a.cut]
    started, checks, pg, pg_sha = refuse(tasks, name, a.commit)
    write_started(tasks, name, a.commit)
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
    entries = arm_entries(pg, genes_ref)
    factors = br.gmom_fits(inputs[ref_role], view.role_stems[ref_role], inputs["start"],
                           view.role_stems["start"])["file"]
    rows = []
    for s in seeds:
        t0 = time.time()
        fl = gv.floor_anchors(t, board, target, ref_stage, s)
        ctx = harness.BuildContext(task=view, board=cboard, inputs=inputs, seed=harness.BUILD_SEED_OFFSET + s,
                                   derived={"t1_held_flat": preregistered})
        arms, rec = cut_arms(ctx, factors, directions, fraction, held_m1, entries)
        shas = {k: cj.prediction_sha(v.X) for k, v in arms.items()}
        rows.append({"seed": s, "floor": float(fl.score),
                     "arms": gv.score_arms(board, target, ref_stage, s, fl.anchors, arms), "sha256": shas,
                     "records": rec, "seconds": round(time.time() - t0, 1)})
        del arms
    res = {**task_record(pg, pg_sha, ext, checks, held_m1, fraction), "task": t.key, "seeds": list(seeds),
           "rows": rows,
           "stamp": gv.stamp(started, [data / t.target, data / t.reference, *[data / f for f in t.inputs.values()],
                                       t1_scope.HELD_FLAT_GENES, t1_scope.CHAMBER_CLASSES, ROOT / ext["file"],
                                       feng_dirs.RELABEL_CONFIG, rc.RECORD, ROOT / union["file"],
                                       tasks / PREGATE_FILE])}
    write_new(tasks / name, json.dumps(res, indent=1, default=float) + "\n")
    print(tasks / name)  # the path only: no score


def cmd_run_lib(a) -> None:
    tasks = RESULTS / a.out / "tasks"
    name = f"lib{a.library}.json"
    started, checks, pg, pg_sha = refuse(tasks, name, a.commit)
    write_started(tasks, name, a.commit)
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
    entries = arm_entries(pg, list(e95.var_names))
    rows = []
    for s in seeds:
        t0 = time.time()
        pick = np.sort(np.random.default_rng([s, 99]).choice(e95.n_obs, rr.N_CELLS, replace=False))
        arms, rec = library_arms(e95, shares, factors, held, held_m1, s, within, fb, directions, fraction, entries)
        stages = {LIB_REFERENCE: pt.to_stage(e95[pick]), **{k: pt.to_stage(v) for k, v in arms.items()}}
        shas = {k: cj.prediction_sha(v.X) for k, v in stages.items()}
        rows.append({"seed": s, **rr.score_arms(board, stages, truth_stage, ref_stage, s), "sha256": shas,
                     "records": rec, "seconds": round(time.time() - t0, 1)})
        del arms, stages
    res = {**task_record(pg, pg_sha, ext, checks, held_m1, fraction), "library": a.library, "seeds": list(seeds),
           "rows": rows,
           "constants": {"split_seed": rr.SPLIT_SEED, "n_cells": rr.N_CELLS, "board": rr.BOARD,
                         **constants(held_m1, pg)},
           "stamp": gv.stamp(started, [*(data / f for f in rr.FILES.values()), t1_scope.HELD_FLAT_GENES,
                                       t1_scope.CHAMBER_CLASSES, ROOT / ext["file"], feng_dirs.RELABEL_CONFIG,
                                       rc.RECORD, ROOT / union["file"], tasks / PREGATE_FILE])}
    write_new(tasks / name, json.dumps(res, indent=1, default=float) + "\n")
    print(tasks / name)  # the path only: no score


def cmd_check_pregate(a) -> int:
    """The `runs` entry's check: pregate.json passed at --commit with its own sha256, and the arm tables are the ones
    it names; no start record and no task file exists."""
    tasks = RESULTS / a.out / "tasks"
    need_row()
    pg, _, bad = load_pregate(tasks)
    bad += pregate_problems(pg, a.commit) if pg is not None else []
    bad += arm_table_problems(pg) if pg is not None else []
    bad += yaml_problems()
    for f in task_files().values():
        if (tasks / f).exists() or (tasks / started_name(f)).exists():
            bad.append(f"a task file or start record exists for {f}")
    if bad:
        sys.exit(f"{ROW} is refused: {bad}")
    return 0


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def task_files() -> dict[str, str]:
    return rc.task_files()


def committed(name: str, tasks: Path = T63_TASKS, arm: str = BASE) -> dict[int, float]:
    rows = json.loads((tasks / name).read_text())["rows"]
    return {int(r["seed"]): float(r["arms"][arm]["score"]) for r in rows}


def _finite(x) -> bool:
    try:
        return bool(np.isfinite(float(x)))
    except (TypeError, ValueError):
        return False


def row_task_constants(c: Mapping) -> bool:
    """A task file's constants are the row's."""
    return (c.get("multiplier") == MULTIPLIER and c.get("practice_fraction") == PRACTICE_FRACTION
            and all(c.get(k) == v for k, v in pregate_constants().items())
            and json.loads(json.dumps(c.get("pregate_constants"))) == json.loads(json.dumps(pregate_constants())))


def seed_problems(g: str, ga: Mapping) -> list[str]:
    """The mask's and each permutation's seed in pregate.json, recomputed here from the row's formulas, written out
    (the mask 20261123 plus the group index; a permutation 20261124 plus 10 x the member index plus the group index;
    the group index from harness.LP_GROUPS, the member index from X002, X019, X020)."""
    if g not in GROUP_MEMBERS or g not in harness.LP_GROUPS:
        return [f"{g}: not a group of the row's consensus"]
    gi = harness.LP_GROUPS.index(g)
    bad = []
    if ga.get("group_index") != gi:
        bad.append(f"{g}: the recorded group index is not harness.LP_GROUPS's")
    if ga.get("mask_seed") != 20261123 + gi:
        bad.append(f"{g}: the mask's seed is not 20261123 plus the group index")
    perm = ga.get("perm") or {}
    if not perm or sorted(perm) != sorted((ga.get("members") or [])[1:]):
        bad.append(f"{g}: the permutation records are not one per new member")
    for m, p in perm.items():
        if m not in NEW_MEMBERS or (p or {}).get("seed") != 20261124 + 10 * MEMBERS.index(m) + gi:
            bad.append(f"{g}: the permutation seed of {m} is not 20261124 plus 10 x the member index plus the group "
                       "index")
    return bad


def guards(loaded: Mapping[str, Mapping | None], gate: Mapping[str, Mapping[int, float]], checks: Mapping,
           pregate: Mapping | None, pregate_sha: str | None, pregate_bad: Sequence[str], started: Mapping[str, bool],
           source_bad: Sequence[str], yaml_bad: Sequence[str] = (), seeds: Sequence[int] = SEEDS,
           read_commit: str | None = None) -> list[str]:
    """Every void condition of row T72; `gate` maps each task to T63's committed m1rc_gmom_x2 per seed; `started`
    maps each task to whether its start record exists."""
    reasons = [*pregate_bad]
    pr = (pregate or {}).get("reading") or {}
    if pregate is None or pr.get("passes") is not True:
        reasons.append("pregate.json is missing or did not pass")
    if pregate is not None and (pregate.get("stamp") or {}).get("tree_clean") is not True:
        reasons.append("the pre-gate ran on an unclean tree")
    if pregate is not None:
        if (pregate.get("constants") is None or json.loads(json.dumps(pregate.get("constants")))
                != json.loads(json.dumps(pregate_constants()))):
            reasons.append("a constant of the pre-gate is not the row's")
        if not data_constants_ok(pregate.get("data_constants")):
            reasons.append("a constant of the tables (the scale, the pseudocount or a cell call) is not the row's")
        if not ATLASES_YAML_SHA256 or pregate.get("yaml_sha256") != ATLASES_YAML_SHA256:
            reasons.append("atlases.yaml's sha256 is not the row's")
        for acc, gsms in (pregate.get("samples") or {}).items():
            pins = {g for e in (PINNED.get(acc) or {}).values() for g in e}
            if set(gsms) & set(NEVER_FETCHED):
                reasons.append(f"{acc}: a never-fetched accession entered")
            if set(gsms) - pins - set(NEVER_FETCHED):
                reasons.append(f"{acc}: an unpinned sample entered")
        if pregate.get("author_labels_read") is not False:
            reasons.append("an author label, cluster or embedding entered")
        for g, ga in (pr.get("groups") or {}).items():
            reasons += seed_problems(g, ga)
            if any(b.get("mask_zeros") != b.get("consensus_zeros") for b in ga.get("mask_bins") or []):
                reasons.append(f"{g}: the mask's zero count in a bin differs from the consensus's")
            for m, p in (ga.get("perm") or {}).items():
                if p.get("bins_keep_values") is not True or p.get("nan_kept") is not True:
                    reasons.append(f"{g}: the permutation of {m} does not keep its bins' values or its NaN positions")
            arms = ga.get("arms") or {}
            vs = {a: (arms.get(a) or {}).get("vector_sha256") for a in COUNTED}
            for a in (MASK, PERM):
                if vs[a] is None or vs[a] == vs[CONS]:
                    reasons.append(f"{g}: {a}'s direction equals mx_cons's (or is missing)")
    reasons += [f"at the read: {x}" for x in yaml_bad]
    commits = set()
    for k, r in loaded.items():
        if r is None:
            reasons.append(f"{k}: no task file" + (" (its start record exists)" if started.get(k) else ""))
            continue
        if not started.get(k):
            reasons.append(f"{k}: a task file without its start record")
        rows = r.get("rows") or []
        if sorted(int(x["seed"]) for x in rows) != list(seeds):
            reasons.append(f"{k}: seeds {[x['seed'] for x in rows]}, not {list(seeds)}")
        want = LIB_ARMS if k in LIBRARIES else ARMS
        absent = sorted({a for x in rows for a in want if a not in x["arms"] or a not in (x.get("sha256") or {})})
        if absent or not rows:
            reasons.append(f"{k}: arms or sha256 {absent} missing")
            continue
        if not all(_finite((x["arms"][a] or {}).get("score")) for x in rows for a in want):
            reasons.append(f"{k}: a score is not finite")
            continue
        stamp = r.get("stamp") or {}
        commits.add(stamp.get("commit"))
        if stamp.get("tree_clean") is not True:
            reasons.append(f"{k}: its tree was not clean at one commit")
        ref = gate.get(k) or {}
        fresh = {int(x["seed"]): float(x["arms"][BASE]["score"]) for x in rows}
        miss = [s for s in seeds if s not in fresh or s not in ref]
        worst = max((abs(fresh[s] - ref[s]) for s in seeds if s not in miss), default=float("inf"))
        if miss or not worst <= REPRO_TOL:
            reasons.append(f"{k}: the gate failed: {BASE} differs from row T63's committed values by up to "
                           f"{worst:.3g}" + (f", seeds missing {miss}" if miss else "") + f" (tolerance {REPRO_TOL:g})")
        if r.get("multiplier") != MULTIPLIER:
            reasons.append(f"{k}: the multiplier is {r.get('multiplier')}, not {MULTIPLIER:g}")
        if not (abs(float(r.get("feng_fraction", float("nan"))) - PRACTICE_FRACTION) <= FRACTION_TOL
                and abs(float(r.get("applied_fraction", float("nan"))) - TEST_FRACTION) <= FRACTION_TOL):
            reasons.append(f"{k}: a fraction is off")
        if not row_task_constants(r.get("constants") or {}):
            reasons.append(f"{k}: a constant is not the row's")
        rel = r.get("relabel") or {}
        if rel.get("sha256") != rc.STATS_SHA256 or rel.get("stats_sha256_hashed") != rc.STATS_SHA256:
            reasons.append(f"{k}: the relabelling's statistics sha256 is not N54's")
        if rel.get("yaml_sha256") != rc.YAML_SHA256:
            reasons.append(f"{k}: the relabelling's settings sha256 is not N54's")
        if rel.get("route") != rc.GROUPS_OPTION or rel.get("record_route") != rc.GROUPS_OPTION:
            reasons.append(f"{k}: the route was not recomputed")
        if rel.get("authors_labels_read") is not False:
            reasons.append(f"{k}: the relabelling's record does not say authors_labels_read false")
        if r.get("pregate_sha256") != pregate_sha or pregate_sha is None:
            reasons.append(f"{k}: ran on another pregate.json than the present one")
        if pregate is not None and json.loads(json.dumps(r.get("arm_tables"))) != json.loads(
                json.dumps(pregate.get("arm_tables"))):
            reasons.append(f"{k}: names other arm tables than pregate.json's")
        if not ATLASES_YAML_SHA256 or r.get("yaml_sha256") != ATLASES_YAML_SHA256:
            reasons.append(f"{k}: atlases.yaml's sha256 is not the row's")
        for x in rows:
            for name in COUNTED:
                h = (x.get("records") or {}).get(name) or {}
                if h.get("outside_equal") is not True:
                    reasons.append(f"{k}: seed {x['seed']}: the byte-check record of {name} is not clean")
    if pregate is not None and commits and (pregate.get("stamp") or {}).get("commit") not in commits:
        reasons.append("pregate.json ran at another commit than the task files")
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


def lib_mean(p: Mapping[str, Mapping]) -> float:
    """The library mean: the mean of the three libraries' ten-seed means."""
    return float(np.mean([p[lib]["mean"] for lib in LIBRARIES]))


def clause_i(ab: Mapping[str, Mapping]) -> list[str]:
    """Clause (i) on a comparison against the base; returns each failed part, named."""
    failed = []
    if not (ab["chamber"]["mean"] >= CHAMBER_BAR and ab["chamber"]["ahead"] >= MIN_AHEAD):
        failed.append(f"(i) the chamber cut's mean is below +{CHAMBER_BAR:g} or fewer than {MIN_AHEAD} of 10 seeds "
                      "are ahead")
    if not ab["whole"]["mean"] > 0:
        failed.append("(i) the whole cut's mean is not above 0")
    if not lib_mean(ab) >= LIB_BAR:
        failed.append(f"(i) the mean over the three libraries is below +{LIB_BAR:g}")
    if not (ab[LIB3]["mean"] > 0 and ab[LIB3]["ahead"] >= LIB3_MIN_AHEAD):
        failed.append(f"(i) library {LIB3}'s mean is not above 0 or fewer than {LIB3_MIN_AHEAD} of 10 seeds are ahead")
    low = [lib for lib in LIBRARIES if not ab[lib]["mean"] >= LIB_FLOOR]
    if low:
        failed.append(f"(i) libraries {low} have a mean below {LIB_FLOOR:g}")
    return failed


def statement(ab: Mapping, am: Mapping, ap: Mapping) -> tuple[bool, list[str]]:
    """(i) mx_cons minus the base, (ii) minus mx_mask, (iii) minus mx_perm; AHEAD when all hold. Every comparison is
    on unrounded values."""
    failed = clause_i(ab)
    if not (am["chamber"]["mean"] >= CTRL_CHAMBER_BAR and am["chamber"]["ahead"] >= MIN_AHEAD):
        failed.append(f"(ii) the chamber cut's mean of mx_cons minus mx_mask is below +{CTRL_CHAMBER_BAR:g} or fewer "
                      f"than {MIN_AHEAD} of 10 seeds are ahead")
    if not lib_mean(am) >= CTRL_LIB_BAR:
        failed.append(f"(ii) the library mean of mx_cons minus mx_mask is below +{CTRL_LIB_BAR:g}")
    if not (ap["chamber"]["mean"] >= CTRL_CHAMBER_BAR and ap["chamber"]["ahead"] >= MIN_AHEAD):
        failed.append(f"(iii) the chamber cut's mean of mx_cons minus mx_perm is below +{CTRL_CHAMBER_BAR:g} or fewer "
                      f"than {MIN_AHEAD} of 10 seeds are ahead")
    return not failed, failed


def t72_reading(loaded: Mapping[str, Mapping | None], gate: Mapping, checks: Mapping, pregate: Mapping | None,
                pregate_sha: str | None, pregate_bad: Sequence[str] = (), started: Mapping[str, bool] | None = None,
                source_bad: Sequence[str] = (), yaml_bad: Sequence[str] = (), read_commit: str | None = None) -> dict:
    started = {k: True for k in loaded} if started is None else started
    reasons = guards(loaded, gate, checks, pregate, pregate_sha, pregate_bad, started, source_bad, yaml_bad,
                     read_commit=read_commit)
    if reasons:
        return {"outcome": "VOID", "reasons": reasons, "follows": FOLLOWS_VOID}
    comp = {f"{a}_minus_{BASE}": (a, BASE) for a in COUNTED}
    comp.update({f"{CONS}_minus_{MASK}": (CONS, MASK), f"{CONS}_minus_{PERM}": (CONS, PERM)})
    tasks = {name: {k: paired(r["rows"], a, b) for k, r in loaded.items()} for name, (a, b) in comp.items()}
    ab, am, ap = (tasks[f"{CONS}_minus_{BASE}"], tasks[f"{CONS}_minus_{MASK}"], tasks[f"{CONS}_minus_{PERM}"])
    ahead, failed = statement(ab, am, ap)
    controls = {a: not clause_i(tasks[f"{a}_minus_{BASE}"]) for a in (MASK, PERM)}  # deciding nothing
    idle = {a: [k for k, r in loaded.items() if all(x["sha256"][a] == x["sha256"][BASE] for x in r["rows"])]
            for a in COUNTED}
    means = {k: {a: float(np.mean([br.score_of(x, a) for x in r["rows"]])) for a in ARMS} for k, r in loaded.items()}
    lib = {name: lib_mean(t) for name, t in tasks.items()}
    j6 = sorted([0.0, ab["chamber"]["mean"]])
    return {"outcome": "AHEAD" if ahead else "NOT AHEAD", "failed": failed, "tasks": tasks, "library_means": lib,
            "controls_meet_clause_i": controls, "idle_tasks": idle, "mean_scores": means, "j6_range": j6,
            "cells_changed": {k: {n: [int(x["records"][n].get("cells_changed", 0)) for x in r["rows"]]
                                  for n in COUNTED} for k, r in loaded.items()},
            "beside": list(BESIDE), "follows": FOLLOWS_AHEAD if ahead else FOLLOWS_END}


def outcome_lines(r: Mapping) -> list[str]:
    """The printed words: the outcome, the idle tasks of every counted arm, what is printed beside, the controls' standing
    (deciding nothing) and what follows. No figure."""
    if r["outcome"] == "VOID":
        return ["Outcome: VOID (" + "; ".join(r["reasons"]) + ")", r["follows"]]
    lines = [AHEAD_WORDS if r["outcome"] == "AHEAD" else "Outcome: MX2s is NOT AHEAD (" + "; ".join(r["failed"]) + ")"]
    lines += [f"{a}: " + IDLE_WORDS.format(task=t) for a in COUNTED for t in r["idle_tasks"][a]]
    lines += list(r["beside"])
    lines.append("Deciding nothing: " + "; ".join(
        f"{a} would {'meet' if ok else 'not meet'} clause (i) against the base"
        for a, ok in r["controls_meet_clause_i"].items()))
    return [*lines, r["follows"]]


def _line(name: str, p: Mapping) -> str:
    g = ", ".join(f"{x} {v:+.2f}" for x, v in p["groups"].items())
    return (f"  {name}: {p['mean']:+.2f} points, sd {p['sd']:.2f}, ahead {p['ahead']} behind {p['behind']} of "
            f"{p['of']}, per seed {p['range'][0]:+.2f} to {p['range'][1]:+.2f}; by group: {g}")


def _fmt(x) -> str:
    return f"{x:.4g}" if isinstance(x, (int, float)) else str(x)


def t72_text(res: Mapping) -> str:
    r = res["reading"]
    lines = [res["label"], f"OUTCOME {r['outcome']}", *outcome_lines(r)]
    if r["outcome"] != "VOID":
        for name, t in r["tasks"].items():
            lines.append(f"{name.replace('_', ' ')}, paired per seed on the same draw (library mean "
                         f"{r['library_means'][name]:+.2f}):")
            lines += [_line(k, p) for k, p in t.items()]
        lines.append("Mean scores (reported):")
        for k, m in r["mean_scores"].items():
            lines.append(f"  {k}: " + ", ".join(f"{a} {x:.2f}" for a, x in m.items()))
        for k, c in r["cells_changed"].items():
            lines.append(f"  {k}: cells changed per seed {c}")
        lo, hi = r["j6_range"]
        lines.append(f"J6's range at the test: [{lo:+.2f}, {hi:+.2f}]")
        pg = res.get("pregate_reading") or {}
        lines.append(f"Changed groups {pg.get('changed_groups')}; drops "
                     f"{[(d['member'], d['group'], d['clauses']) for d in pg.get('drops') or []]}")
        for g, ga in (pg.get("groups") or {}).items():
            lines.append(f"  {g}: members {ga.get('members')}, gene counts {ga.get('member_gene_counts')}")
            for a, v in (ga.get("arms") or {}).items():
                lines.append(f"    {a}: k_g {_fmt(v.get('k_g'))}, zero share {_fmt(v.get('zero_share'))}, Spearman "
                             f"with X002 over S_g {_fmt(v.get('spearman_x002'))}")
        for g, d in (pg.get("members") or {}).items():
            for m, f in d.items():
                lines.append(f"  member {m} {g}: Spearman with X002 {f.get('spearman')}, cells {f.get('n_early')} and "
                             f"{f.get('n_late')}")
    lines.append(f"Tries counted: {res['tries']} (three counted arms on five tasks; the row numbers them)")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def read(out: Path, commit: str) -> dict:
    need_row()
    gv.refuse_existing(out, OUT_FILES)
    started = head_is(commit)
    tasks = out / "tasks"
    paths = {k: tasks / f for k, f in task_files().items()}
    loaded = {k: json.loads(p.read_text()) if p.exists() else None for k, p in paths.items()}
    begun = {k: (tasks / started_name(f)).exists() for k, f in task_files().items()}
    pregate, pg_sha, pg_bad = load_pregate(tasks)
    if pregate is None or (pregate.get("reading") or {}).get("passes") is not True:
        sys.exit(f"{ROW} is refused: pregate.json is missing or did not pass, so no try was counted and no read "
                 "exists (nothing is written)")
    if not any(begun.values()):
        sys.exit(f"{ROW} is refused: no start record exists, so no task run loaded and no try was counted (nothing "
                 "is written)")
    gate = {k: committed(f) for k, f in task_files().items()}
    if pregate is not None:
        pg_bad = [*pg_bad, *arm_table_problems(pregate)]
    checks = rc.relabel_checks(need_stats=False)
    reading = t72_reading(loaded, gate, checks, pregate, pg_sha, pg_bad, begun, all_source_bad(), yaml_problems(),
                          read_commit=started["commit"])
    res = {"label": LABEL, "row": ROW, "reading": reading, "tries": TRIES_COUNT, "pregate_sha256": pg_sha,
           "pregate_reading": (pregate or {}).get("reading"), "open": list(OPEN)}
    res["stamp"] = gv.stamp(started, [p for p in [*paths.values(), tasks / PREGATE_FILE] if p.exists()])
    out.mkdir(parents=True, exist_ok=True)
    write_new(out / OUT_FILES[0], json.dumps(res, indent=1, default=float) + "\n")
    write_new(out / OUT_FILES[1], t72_text(res))
    print("\n".join(outcome_lines(reading)))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for cmd in ("pregate", "check-pregate", "read"):
        p = sub.add_parser(cmd)
        p.add_argument("--out", required=True)
        p.add_argument("--commit", required=True)
    for cmd, flag, choices in (("run-cut", "--cut", list(CUTS)), ("run-lib", "--library", list(LIBRARIES))):
        p = sub.add_parser(cmd)
        p.add_argument(flag, required=True, choices=choices)
        p.add_argument("--out", required=True)
        p.add_argument("--commit", required=True)
        p.add_argument("--seeds", default=",".join(map(str, SEEDS)))
        p.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    if a.cmd == "read":
        read(RESULTS / a.out, a.commit)
        return 0
    if a.cmd == "pregate":
        return cmd_pregate(a)
    if a.cmd == "check-pregate":
        return cmd_check_pregate(a)
    {"run-cut": cmd_run_cut, "run-lib": cmd_run_lib}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
