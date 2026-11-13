#!/usr/bin/env python
"""Row N28 of experiments/harness/TRIES.md (T1; a diagnostic, not a try): does the library-3 base already carry the
strong-null bet of tries T12a and T12b?

Since rows R6 and R7 (2026-11-13) T1's first upload is drawn from E10.5's library 3 alone, the dissected-heart library,
while the scorer's reference stage is expected to be E10.5 as released, all three libraries (assumption A3). If the
kind change between embryo portions and a dissected heart is what made E9.5 -> E10.5 strong-null (mostly highly
expressed genes going down), the library-3 base already carries it in its own dp = pb(base) - pb(reference), and the
prior of T12a and T12b (vec.t1_scope.apply_down_prior) has little left to add. N28 measures that overlap.

What it reads: the released E10.5 file only (T1/E10.5_RNA.h5ad under the data directory, loaded as
experiments/builds/build_test.py t1 loads it; E10.5 is an input of the test task, assumption A20), the T1 panel's gene
names (vec.boards, names only), and the two committed gene lists the builder holds flat (the pre-registered list and
R2's union, through the builder's own `_t1_held`). It reads no target of any fold, backtest or board, compares nothing
with a truth, computes no score and builds no file. The pools are the builder's own: library 3 as build_test.py t1
identifies it (cell names ending in _<T1_DISSECTED_LIBRARY>), the chamber categories of harness.t1_chamber_shares, and
for base B the builder's `_t1_mid_rows` at the builder's default build seed.

Definition (every constant below is the row's):
  pb_ref: the mean over all released E10.5 cells of the file's log1p CP10k values, per panel gene.
  Base A (the first upload's, `--base lib3`): the library-3 cells of the chamber categories at that pool's own mix
    (harness.t1_chamber_shares of the pool), before any draw, direction step or R1 to R3b edit; pb_A the
    category-weighted mean of its category means (equal to the plain mean over those cells).
  Base B (the second upload's, `--base mid --cm-share 0.60`; reported only): the mid pool of `_t1_mid_rows` at the
    builder's default seed, its chamber shares with the cardiomyocyte share of heart cells set to T1_MID_CM_SHARE
    (t1_scope.with_share, as the builder does); pb_B the category-weighted mean of its category means.
  dp = pb_base - pb_ref; held: the builder's held-flat genes; level, change = down_prior_bins(pb_ref, dp,
    DownPriorTable.level_only()).
  eligible: level bin among the first DP_PROMOTED_LEVELS and not the table's last, not held, pb_ref >= DP_MIN_PB.
  top: the DP_BUDGET eligible genes with the highest reference level.
  redundancy: the share of top whose rank by dp (most negative first) is below DP_BOTTOM.
  GATE (fixed in N28's row before it runs): "kill" if redundancy on base A >= REDUNDANCY_KILL, otherwise "open".

Reported, deciding nothing (per base): how many genes down_prior_order returns (capped at DP_BUDGET) and their median
dp; the share of dp's bottom DP_BOTTOM genes whose reference level rank is below LEVEL_RANK_CUT; the Spearman
correlation of dp with -pb_ref over genes with pb_ref >= DP_MIN_PB; and base B's redundancy. No gene name, no per-gene
value, no share of any pool and no cell number of any stage or sample is written.

Outputs, in experiments/t1_strong_null/results/n28_lib3_overlap/: n28.json and n28.txt, each with the commit, the
clean-tree flag (vec.provenance), the UTC time and the inputs' sha256, and a vec.provenance stamp each. N28 runs once:
it refuses when its output folder already holds anything (a rerun needs a new TRIES row).

Usage (DGX, repo root; after N28's row is committed):
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/t1_strong_null/n28_lib3_overlap.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

from vec import harness, t1_scope
from vec.boards import REPO_ROOT

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "n28_lib3_overlap"
sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
import build_test  # noqa: E402  (the builder's pools, library mark, mid draw and held-flat genes, reused as they are)

E105_FILE = "T1/E10.5_RNA.h5ad"  # the released reference stage, the one file read
E105_SOURCE = "E10.5_RNA"  # its stem, for the chamber config's per-file rules (as the builder passes it)
PANEL_BOARD = "T1:test"  # the panel the prior ranks (its gene names only; no limit of the board is read or written)
REDUNDANCY_KILL = 0.70  # the gate: base A already carries the bet if at least this share of top lies in dp's bottom
LEVEL_RANK_CUT = t1_scope.DP_LEVEL_EDGES[t1_scope.DP_PROMOTED_LEVELS - 1]  # 800: the promoted levels' rank edge
MID_SEED = harness.BUILD_SEED_OFFSET + 1  # build_test.py's default build seed, for the mid base's draw
VERDICTS = ("kill", "open")


def verdict(redundancy_a: float) -> str:
    """N28's gate: "kill" if base A's redundancy is at least REDUNDANCY_KILL, otherwise "open"."""
    return VERDICTS[0] if float(redundancy_a) >= REDUNDANCY_KILL else VERDICTS[1]


def _ranks_desc(values) -> np.ndarray:
    """Each gene's rank when sorted from the largest down, ties by position (down_prior_bins' level rank)."""
    order = np.argsort(-np.asarray(values, dtype=np.float64), kind="stable")
    out = np.empty(order.size, dtype=np.int64)
    out[order] = np.arange(order.size)
    return out


def _ranks_asc(values) -> np.ndarray:
    """Each gene's rank when sorted from the most negative up, ties by position (down_prior_bins' change rank)."""
    order = np.argsort(np.asarray(values, dtype=np.float64), kind="stable")
    out = np.empty(order.size, dtype=np.int64)
    out[order] = np.arange(order.size)
    return out


def top_genes(pb_ref, dp, held) -> np.ndarray:
    """Positions of `top`: the DP_BUDGET eligible genes with the highest reference level (level bin among the first
    DP_PROMOTED_LEVELS and not the table's last, not held, pb_ref >= DP_MIN_PB), highest level first."""
    table = t1_scope.DownPriorTable.level_only()
    pb_ref = np.asarray(pb_ref, dtype=np.float64)
    level, _ = t1_scope.down_prior_bins(pb_ref, dp, table)
    ok = ((level < min(len(table.level_edges), t1_scope.DP_PROMOTED_LEVELS)) & (level != len(table.level_edges))
          & ~np.asarray(held, dtype=bool) & (pb_ref >= t1_scope.DP_MIN_PB))
    idx = np.flatnonzero(ok)
    return idx[np.argsort(_ranks_desc(pb_ref)[idx], kind="stable")][:t1_scope.DP_BUDGET]


def _spearman(a, b) -> float | None:
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.size < 3 or a.std() == 0 or b.std() == 0:
        return None
    return float(spearmanr(a, b).correlation)


def overlap(pb_ref, pb_base, held, genes: Sequence[str]) -> dict:
    """N28's numbers on one base (the module docstring): the redundancy and the reported quantities. Summaries only."""
    pb_ref, pb_base = np.asarray(pb_ref, dtype=np.float64), np.asarray(pb_base, dtype=np.float64)
    held = np.asarray(held, dtype=bool)
    dp = pb_base - pb_ref
    table = t1_scope.DownPriorTable.level_only()
    level, change = t1_scope.down_prior_bins(pb_ref, dp, table)
    top = top_genes(pb_ref, dp, held)
    dp_rank = _ranks_asc(dp)
    redundancy = float((dp_rank[top] < t1_scope.DP_BOTTOM).mean()) if top.size else None
    order = t1_scope.down_prior_order(table, level, change, pb_ref, dp, genes, held)[:t1_scope.DP_BUDGET]
    bottom = np.flatnonzero(dp_rank < t1_scope.DP_BOTTOM)
    expressed = pb_ref >= t1_scope.DP_MIN_PB
    return {"redundancy": redundancy, "top_genes": int(top.size),
            "prior_order_genes": int(order.size),
            "prior_order_median_dp": float(np.median(dp[order])) if order.size else None,
            "bottom_share_high_level": float((_ranks_desc(pb_ref)[bottom] < LEVEL_RANK_CUT).mean()) if bottom.size
            else None,
            "spearman_dp_vs_minus_pb_ref": _spearman(dp[expressed], -pb_ref[expressed])}


def _panel_matrix(a, genes: Sequence[str]):
    cols = np.asarray(a.var_names.get_indexer(list(genes)))
    if (cols < 0).any():
        sys.exit(f"{int((cols < 0).sum())} panel genes are missing from the E10.5 file (names not shown)")
    return sparse.csr_matrix(a.X[:, cols], dtype=np.float64)


def weighted_pseudobulk(pool, genes: Sequence[str], shares: Mapping[str, float]) -> np.ndarray:
    """The category-weighted mean of the pool's per-category means of log1p CP10k, at `shares` (categories with a
    positive share; each must hold cells in the pool). The pool is classified as the builder's pools are."""
    b = harness._shallow(pool)
    cats = np.asarray(t1_scope.classify(b, source=E105_SOURCE, spec=t1_scope.load_spec(), key_added=None)).astype(str)
    X = _panel_matrix(pool, genes)
    keep = {c: float(s) for c, s in shares.items() if s > 0}
    total = sum(keep.values())
    out = np.zeros(len(genes), dtype=np.float64)
    for c, s in sorted(keep.items()):
        rows = np.flatnonzero(cats == c)
        if rows.size == 0:
            raise ValueError(f"the pool holds no cell of the category {c} it is weighted by")
        out += (s / total) * np.asarray(X[rows].mean(axis=0)).ravel()
    return out


def library3_mask(ref) -> np.ndarray:
    """Library 3 exactly as build_test.py t1 marks it for --base lib3 and mid (its lines and its refusal, with its
    constant T1_DISSECTED_LIBRARY): cell names ending in _<library number>."""
    names = [str(x) for x in ref.obs_names]
    library = np.array([x.rsplit("_", 1)[-1] if "_" in x else "" for x in names])
    kinds = set(library.tolist())
    if not all(k.isdigit() for k in kinds) or build_test.T1_DISSECTED_LIBRARY not in kinds or len(kinds) < 2:
        sys.exit("the reference's cell names do not separate library 3 from the others (every name must end in "
                 "_<library number>, and library 3 must be one of at least two)")
    return library == build_test.T1_DISSECTED_LIBRARY


def run(ref, genes: Sequence[str], held_genes, seed: int = MID_SEED) -> dict:
    """N28 on the loaded E10.5 file: both bases, the gate on base A. No cell number enters the result."""
    genes = [str(g) for g in genes]
    held_set = {str(g) for g in held_genes}
    held = np.array([g in held_set for g in genes], dtype=bool)
    pb_ref = np.asarray(_panel_matrix(ref, genes).mean(axis=0)).ravel()
    in3 = library3_mask(ref)
    pool_a = ref[np.flatnonzero(in3)].copy()
    shares_a = harness.t1_chamber_shares(harness._shallow(pool_a), E105_SOURCE)  # the pool's own chamber mix
    base_a = overlap(pb_ref, weighted_pseudobulk(pool_a, genes, shares_a), held, genes)
    rows, _ = build_test._t1_mid_rows(ref, in3, seed)  # the builder's draw; its record (shares) is not kept
    pool_b = ref[rows].copy()
    shares_b = t1_scope.with_share(harness.t1_chamber_shares(harness._shallow(pool_b), E105_SOURCE), "cm",
                                   build_test.T1_MID_CM_SHARE)
    base_b = overlap(pb_ref, weighted_pseudobulk(pool_b, genes, shares_b), held, genes)
    if base_a["redundancy"] is None:
        raise ValueError("no eligible gene on base A: the gate is undefined")
    return {"verdict": verdict(base_a["redundancy"]),
            "gate": {"statistic": "base A's redundancy", "kill_if_at_least": REDUNDANCY_KILL},
            "base_a": {"pool": "build_test.py t1 --base lib3: library-3 cells of the chamber categories at the pool's "
                               "own mix, before any draw or edit",
                       "categories": sorted(c for c, s in shares_a.items() if s > 0), **base_a},
            "base_b_reported_only": {"pool": "build_test.py t1 --base mid --cm-share 0.60: _t1_mid_rows at the "
                                             "builder's default build seed, cardiomyocytes at the ruled share of "
                                             "heart cells", "seed": int(seed), "cm_share_of_heart":
                                     build_test.T1_MID_CM_SHARE,
                                     "categories": sorted(c for c, s in shares_b.items() if s > 0), **base_b},
            "held_genes_in_panel": int(held.sum())}


def report_text(res: dict) -> str:
    def f(x) -> str:
        return "n/a" if x is None else f"{x:.3f}"

    a, b = res["base_a"], res["base_b_reported_only"]
    lines = ["N28 (experiments/harness/TRIES.md): does the library-3 base already carry the strong-null bet?",
             f"Verdict: {res['verdict']} (kill if base A's redundancy >= {REDUNDANCY_KILL:.2f}, otherwise open)", ""]
    for name, r in (("Base A (first upload, --base lib3)", a),
                    ("Base B (second upload, --base mid; reported only)", b)):
        lines += [f"{name}: redundancy {f(r['redundancy'])} (share of the top {r['top_genes']} eligible genes by "
                  f"reference level whose dp rank is below {t1_scope.DP_BOTTOM})",
                  f"  down_prior_order genes (capped at {t1_scope.DP_BUDGET}): {r['prior_order_genes']}, their median "
                  f"dp {f(r['prior_order_median_dp'])}",
                  f"  share of dp's bottom {t1_scope.DP_BOTTOM} genes with reference level rank below "
                  f"{LEVEL_RANK_CUT}: {f(r['bottom_share_high_level'])}",
                  f"  Spearman of dp with -pb_ref (pb_ref >= {t1_scope.DP_MIN_PB}): "
                  f"{f(r['spearman_dp_vs_minus_pb_ref'])}", ""]
    s = res["stamp"]
    lines += [f"Inputs (sha256): {res['inputs']}",
              f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="N28: does the library-3 base already carry the strong-null bet?")
    ap.add_argument("--data", default=None, help="the released data (default: data/raw)")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N28 has already written its outputs; a rerun needs a new TRIES row")

    from vec import provenance
    from vec.boards import load_boards
    from vec.prep import load_stage

    started = provenance.git_state(provenance.REPO_ROOT)
    data = Path(a.data) if a.data else harness.DATA_DIR
    genes = list(load_boards()[PANEL_BOARD].genes)  # gene names only
    union = build_test.hygiene.load_union()  # R2's union, as build_test.py t1 gets it (refused unless pinned)
    held_genes, _ = build_test._t1_held(union, (), genes)  # the builder's held-flat genes (no Feng group is read)
    ref = load_stage(data / E105_FILE, lineage=False)
    res = run(ref, genes, held_genes, MID_SEED)
    inputs = {E105_FILE: harness._sha256(data / E105_FILE),
              str(Path(t1_scope.HELD_FLAT_GENES).relative_to(REPO_ROOT)): harness._sha256(Path(
                  t1_scope.HELD_FLAT_GENES)),
              union["file"]: union["sha256"]}
    ended = provenance.git_state(provenance.REPO_ROOT)
    clean = None if started[1] is None or ended[1] is None else (not started[1] and not ended[1]
                                                                 and started[0] == ended[0])
    res = {"label": "N28 (experiments/harness/TRIES.md): the library-3 base's overlap with the strong-null bet; the "
                    "released E10.5 file only, no target, no score, no file built",
           **res, "inputs": inputs,
           "settings": {"redundancy_kill": REDUNDANCY_KILL, "dp_budget": t1_scope.DP_BUDGET,
                        "dp_bottom": t1_scope.DP_BOTTOM, "dp_min_pb": t1_scope.DP_MIN_PB,
                        "dp_promoted_levels": t1_scope.DP_PROMOTED_LEVELS, "level_rank_cut": LEVEL_RANK_CUT,
                        "level_edges": list(t1_scope.DP_LEVEL_EDGES), "change_edges": list(t1_scope.DP_CHANGE_EDGES),
                        "panel_board": PANEL_BOARD, "mid_seed": MID_SEED},
           "stamp": {"commit": ended[0], "tree_clean": clean,
                     "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    text = report_text(res)
    paths = {out / "n28.json": json.dumps(res, indent=1) + "\n", out / "n28.txt": text}
    for p, body in paths.items():
        p.write_text(body)
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
