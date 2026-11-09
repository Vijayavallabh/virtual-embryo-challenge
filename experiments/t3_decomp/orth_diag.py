#!/usr/bin/env python
"""Tries U3 and U3d of experiments/harness/TRIES.md (T3; one row, one code path), built on the pre-registered T3 script
(t3_decompose.py and build_submission.py, hash-pinned: imported, never edited).

U3, the stage-orthogonal reading. The barred file's factor (`transfer:mab21l2_all`, the Mab21l2 knockout against WT
E8.75) scores 60.77 on the Gata4 replicates; the same knockout against its own WT E9.5 (`transfer:mab21l2_vs95`) 49.14.
N8 put corr(f_M, d) at 0.42 for the WT developmental axis d (E9.5 against E8.25_late, in the WT E8.75 vocabulary, built
without WT E8.75). Nobody has separated that stage part from the rest by score. U3 removes it from the barred factor
and scores what is left:
  phi_perp = exp(log phi_B - b d), b the OLS slope of log phi_B on d over the panel's genes, then gene_factor's rescale
  (the same label mask as the barred factor's: labels with >= 100 cells on both sides) and gene_factor's clip.
It starts from the hashed factor file (sha256 in gata4_decisions.json), so B-perp is exactly the barred factor minus
its stage part. B-perp's reading: q_hat = (mean B-perp - mean vs95) / (mean B - mean vs95) over the 20 (replicate,
seed) rows, the means being the committed ones (60.7745 and 49.1371). It is a share of score on Gata4, not N11's
probability q: B keeps its stage term.

U3d, the depth-matched Gata4 reading: both cells N11 had to judge (s1 and the barred file in the world where the target
is not Gata4-shallow), measured. In each (replicate, seed) the WT E8.75 reference working copy and the WT pool are
thinned to the replicate's per-label depth with the T3 script's rescale_depth; then `transfer:mab21l2_all` (B, no depth
step) and `xrep:F+depth` (s1's stand-in, its depth step applied on top of the already matched pool, so it overshoots
as s1 would against a not-shallow target) are scored with anchors recomputed on the matched reference. `n11` then
reruns N11 (experiments/gate/t3_backup_sensitivity.py) with those measured rows in the place of the judgement cell and
the shallow-file barred value.

Scoring is the T3 script's: score_seed's truth half, ceiling half, reference working copy, anchors, record fields and
streams. B and B-perp share the stream `stable_rng(seed, "transfer:mab21l2_all")`, so they draw the same cells; a
rerun of B is compared with the committed rows, and if it cannot be matched B counts as one more try.

Reads released and unlocked files only: WT E8.25_late, E8.5, E8.75 and E9.5; the Mab21l2 knockout; the two Gata4
replicates as truths, as D11b read them. Never a board file, a board anchor, a test file or a leaderboard. No number
derived from Gata4 enters B-perp's construction: the gate below builds its context without any Gata4 file. Nothing
here changes an upload: B-perp could only ever be the slot-2 swap file, by the captain's ruling.

The gate comes before any scoring: `axis` reads no Gata4 file and no score. It stops (exit status 1, outputs written)
if the axis fails its self-check against experiments/research_2026_11/results/d5d6.json (corr_fM_d within 0.01 of
0.4225) or the barred factor is not reproduced from the run's own label mask. If the removal is a no-op (SD of log phi
changes by less than 0.02), U3 is inconclusive and B-perp is not scored (no try), while U3d, which does not depend on
the removal, still runs. `score` refuses to run unless that gate passed with this script. No second axis is scored
(the secondary axis is reported only) and no follow-up variant is run.

Usage (DGX, repository root; the row is committed first and the captain signs off on its plain summary before the run;
at most 16 BLAS threads, NUMPY_MADVISE_HUGEPAGE=0):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
      uv run --frozen --no-sync python -I -B experiments/t3_decomp/orth_diag.py axis
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
      uv run --frozen --no-sync python -I -B experiments/t3_decomp/orth_diag.py score --workers 8 --threads 2
  uv run --frozen --no-sync python -I -B experiments/t3_decomp/orth_diag.py n11

Outputs in experiments/t3_decomp/results/orth_diag/ (scores on the Gata4 validation files; no test-board count, limit
or anchor is read or written): axis_gate.json and .txt, orth_records.csv, orth_decision.json and .txt, n11_measured.json
and .txt, each stamped (commit, clean-tree flag by vec.provenance's rule, UTC, input hashes) and with a
vec.provenance build stamp. B-perp's per-gene factor goes to runs/t3_decomp/orth_diag/factors/ (the per-gene factors
stay with the runs), its sha256 in the gate. Every command refuses to overwrite its outputs unless --rerun-row names
a row of experiments/harness/TRIES.md, other than U3 and U3d, that registers the rerun (the files of a run that
stopped early are removed by hand).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import re
import sys
import tempfile
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in (REPO / "src", HERE, REPO / "experiments" / "gate"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import build_submission as bs  # noqa: E402  (first: it records this build's start through vec.provenance)
import estimate as E  # noqa: E402
import t3_backup_pricing as P  # noqa: E402
import t3_backup_sensitivity as S  # noqa: E402
import t3_decompose as t3d  # noqa: E402

from vec import provenance  # noqa: E402

DATA = REPO / "data" / "raw"
RUN_DIR = REPO / "runs" / "t3_decomp" / "gata4_20261028T122425Z"
DECISION = HERE / "results" / "gata4" / "gata4_decisions.json"
RECORDS = HERE / "results" / "gata4" / "gata4_records.csv"
META = HERE / "results" / "gata4" / "gata4_meta.json"
D5D6 = REPO / "experiments" / "research_2026_11" / "results" / "d5d6.json"
TRIES = REPO / "experiments" / "harness" / "TRIES.md"
OUT = HERE / "results" / "orth_diag"
FACTOR_DIR = REPO / "runs" / "t3_decomp" / "orth_diag" / "factors"
OWN_ROWS = ("U3", "U3d")  # this script's rows: a rerun is registered under another row
ROW_ID = re.compile(r"[A-Z][0-9]+[a-z]?")  # a TRIES row's ID (R2, N12, K2b, ...)
WT_FILES = {"E8.25_late": "T2/E8.25_late.h5ad", "E8.5": "T2/E8.5.h5ad", "E8.75": "T2/E8.75.h5ad",
            "E9.5": "T2/E9.5.h5ad"}
KO_FILES = {"rep1": "T3/E8.75_gata4_ko_1.h5ad", "rep2": "T3/E8.75_gata4_ko_2.h5ad",
            "mab21l2": "T3/E9.5_mab21l2_ko.h5ad"}

AXIS_JSON, AXIS_TXT = "axis_gate.json", "axis_gate.txt"
PERP_CSV = "transfer_mab21l2_all_perp.csv"
RECORDS_CSV, DECISION_JSON, DECISION_TXT = "orth_records.csv", "orth_decision.json", "orth_decision.txt"
N11_JSON, N11_TXT = "n11_measured.json", "n11_measured.txt"

BARRED = "transfer_mab21l2_all"  # the factor file; its scored variant is "transfer:mab21l2_all"
B_STREAM = "transfer:mab21l2_all"
S1_STREAM = "xrep:F+depth"
V_B = "transfer:mab21l2_all"                # the rerun of B (the pairing check)
V_PERP = "transfer:mab21l2_all_perp"        # U3
V_B_DM = "dm:transfer:mab21l2_all"          # U3d: B on the depth-matched pool, against the matched reference
V_S1_DM = "dm:xrep:F+depth"                 # U3d: s1's stand-in, its depth step on top of the matched pool
BARRED_MIN_CELLS = 100  # t3_decompose.transfer_factors: a label enters the barred factor with >= 100 cells on each side
AXIS_MIN_CELLS = 100    # t3_checks.py d5d6: a label enters the axis with >= 100 cells in both stages
SELF_CHECK_TOL = 0.01   # the axis's corr(f_M, d) against the committed D5 value (TRIES U3: 0.4225 +- 0.01)
REPRO_TOL_LOG = 0.01    # the barred factor recomputed from this run's label mask against the hashed file, in |log phi|.
#                         It guards the mask of the rescale; a wrong mask moves some genes by far more than this, while
#                         a few relabelled Mab21l2 cells stay well inside it
NOOP_SD = 0.02          # U3's no-op: SD of log phi changes by less than this
REPLAY_TOL = 0.02       # a rerun row of B within this of the committed row (T3 points) and with the same floor flag
COMMITTED_MEANS = {"transfer:mab21l2_all": 60.7745, "transfer:mab21l2_vs95": 49.1371}  # means over the 20 rows
MEANS_TOL = 0.005       # the committed means, as the row states them, against the records file


# ---------------------------------------------------------------------------------------------------------------------
# Stamps


def _rel(path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO).as_posix()
    except ValueError:
        return str(path)


def tree_state(exempt: Sequence[Path]) -> dict:
    """HEAD and whether the tree is clean by vec.provenance's rule, the outputs exempt (None when git cannot tell)."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    return {"commit": commit, "clean": None if unclean is None else not unclean}


def run_stamp(started: dict, exempt: Sequence[Path], rerun_row: str | None) -> dict:
    ended = tree_state(exempt)
    clean = (None if started["clean"] is None or ended["clean"] is None
             else started["clean"] and ended["clean"] and started["commit"] == ended["commit"])
    return {"commit": ended["commit"], "tree_clean": clean, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "rerun_row": rerun_row}


def tries_rows(path: Path = TRIES) -> set[str]:
    """The row IDs in experiments/harness/TRIES.md's tables: each table row's first word, where it is a row ID."""
    words = (ln[2:].split("|", 1)[0].split() for ln in path.read_text(encoding="utf-8").splitlines() if ln[:2] == "| ")
    return {w[0] for w in words if w and ROW_ID.fullmatch(w[0])}


def refuse_existing(paths: Sequence[Path], rerun_row: str | None) -> None:
    """A command runs once: its outputs are not rewritten, unless the row that registers the rerun is named. A rerun's
    row must be a row of experiments/harness/TRIES.md other than this script's own (registered before the rerun);
    the files of a run that stopped early are removed by hand."""
    if rerun_row is not None and (rerun_row in OWN_ROWS or not ROW_ID.fullmatch(rerun_row)
                                  or rerun_row not in tries_rows()):
        sys.exit(f"--rerun-row {rerun_row}: not a row of experiments/harness/TRIES.md other than "
                 f"{' and '.join(OWN_ROWS)} (a rerun is registered there before it runs)")
    existing = [str(p) for p in paths if Path(p).exists()]
    if existing and not rerun_row:
        sys.exit(f"orth_diag: refusing to overwrite {existing} (a rerun needs a new experiments/harness/TRIES.md row, "
                 "named with --rerun-row)")


def input_hashes(paths: Sequence[Path]) -> dict:
    return {_rel(p): t3d.sha256(p) for p in paths}


def write_outputs(items: Sequence[tuple[Path, str]]) -> None:
    """Write each (path, text), then its vec.provenance build stamp."""
    for path, text in items:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        provenance.write_stamp(path)


# ---------------------------------------------------------------------------------------------------------------------
# The axis and the orthogonalisation (U3)


def context_paths(data: Path, with_targets: bool) -> dict:
    """The T3 script's `paths` for build_context: the gata4 run's (with the two replicates) or, for the target-free
    gate, without any Gata4 file (as t3_checks.py d5d6)."""
    paths: dict = {"wt": str(data / WT_FILES["E8.75"]), "targets": [],
                   "stage:E9.5": str(data / WT_FILES["E9.5"]), "stage:E8.25_late": str(data / WT_FILES["E8.25_late"])}
    if with_targets:
        paths["targets"] = [("rep1", str(data / KO_FILES["rep1"])), ("rep2", str(data / KO_FILES["rep2"]))]
        paths.update({"series:E8.25_late": paths["stage:E8.25_late"], "series:E9.5": paths["stage:E9.5"],
                      "series:E8.5": str(data / WT_FILES["E8.5"]), "transfer_ko": str(data / KO_FILES["mab21l2"])})
    return paths


def centred_log(phi) -> np.ndarray:
    v = np.log(np.asarray(phi, dtype=np.float64))
    return v - v.mean()


def stage_axis(mu: dict, cnt: dict, V, fc: dict, late: str, early: str) -> np.ndarray:
    """d5d6's axis: the centred log of t3_decompose.gene_factor(later stage against earlier, label weights V.p, the
    strict non-Mesp1 labels with >= 100 cells in both stages), in the WT E8.75 vocabulary."""
    mask = np.asarray(V.nm, dtype=bool) & (np.asarray(cnt[late]) >= AXIS_MIN_CELLS) & (
        np.asarray(cnt[early]) >= AXIS_MIN_CELLS)
    phi, _ = t3d.gene_factor(mu[late], mu[early], V.p, mask, fc["eps"], tuple(fc["clip"]))
    return centred_log(phi)


def ladder_means(ctx, e85_path) -> tuple[dict, dict]:
    """Mean proportions and label cell counts per WT stage in the WT E8.75 vocabulary, exactly as t3_checks.py d5d6
    reads them: the T3 context's stage pools (E9.5, E8.25_late) and, for the secondary axis, WT E8.5 called into the
    vocabulary with no confidence filter."""
    from vec.prep import load_stage

    V, hp = ctx.V, ctx.rules["harmonisation"]
    mu = {sp.name: t3d.props_type_means(sp.pool.C, sp.pool.lab, V.L) for sp in ctx.stages}
    cnt = {sp.name: np.asarray(sp.count) for sp in ctx.stages}
    if "E8.5" not in mu:
        a = load_stage(e85_path, genes=ctx.genes, drop_duplicates=True)
        calls, _ = ctx.harm.calls(a, seed=hp["seed"])
        lab = t3d.calls_to_index(V, calls)
        mu["E8.5"] = t3d.props_type_means(t3d.dense_counts(a), lab, V.L)
        cnt["E8.5"] = np.bincount(lab[lab >= 0], minlength=V.L)
    return mu, cnt


def barred_denominator(ctx, mab_path) -> tuple[np.ndarray, np.ndarray, np.ndarray, object]:
    """t3_decompose.transfer_factors' computation of the barred factor in this context: the label mask (labels with
    >= 100 Mab21l2 cells and >= 100 WT cells), the factor it gives (to compare with the hashed file) and the denominator
    the factor's rescale uses (gene_factor's den: the mask's WT-share-weighted mean proportions)."""
    from vec.prep import load_stage

    hp, V, wt, fc = ctx.rules["harmonisation"], ctx.V, ctx.wt, ctx.rules["components"]["factor"]
    a = load_stage(mab_path, genes=ctx.genes, drop_duplicates=True)
    calls, key = ctx.harm.calls(a, seed=hp["seed"])
    lab = t3d.calls_to_index(V, calls)
    C = t3d.dense_counts(a)
    mu = t3d.props_type_means(C, lab, V.L)
    cnt = np.bincount(lab[lab >= 0], minlength=V.L)
    both = (cnt >= BARRED_MIN_CELLS) & (wt.count >= BARRED_MIN_CELLS)
    phi, den = t3d.gene_factor(mu, wt.mu, V.p, both, fc["eps"], tuple(fc["clip"]))
    return both, phi, den, key


def remove_axis(phi_b, d, den, clip) -> dict:
    """B-perp: exp(log phi_B - b d), b the OLS slope of log phi_B (centred) on the centred axis d over all genes, then
    gene_factor's rescale (den @ phi = den.sum()) and clip. d must be centred, as stage_axis returns it, so the removal
    keeps the mean of log phi and the rescale only moves it back to unit total proportion."""
    phi_b = np.asarray(phi_b, dtype=np.float64)
    d = np.asarray(d, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    v = np.log(phi_b)
    vc = v - v.mean()
    b = float(vc @ d / (d @ d))
    raw = np.exp(v - b * d)
    scaled = raw / (den @ raw / max(den.sum(), 1e-12))
    phi = np.clip(scaled, *clip)
    return {"phi_perp": phi, "b": b, "corr_before_clip_with_axis": float(np.corrcoef(centred_log(scaled), d)[0, 1]),
            "share_of_genes_clipped": float(np.mean((scaled < clip[0]) | (scaled > clip[1])))}


def noop(sd_b: float, sd_perp: float) -> bool:
    return abs(sd_perp - sd_b) < NOOP_SD


def gate_outcome(corr_fm_d: float, corr_committed: float, repro_max_log: float, labels: int, labels_committed: int,
                 sd_b: float, sd_perp: float) -> dict:
    """The gate of U3's row: the axis's self-check and the barred factor's reproduction decide whether `score` may run
    (`proceed`); a no-op removal leaves U3 inconclusive with B-perp unscored while U3d still runs (`noop`)."""
    self_check = bool(np.isfinite(corr_fm_d) and abs(corr_fm_d - corr_committed) <= SELF_CHECK_TOL)
    repro = bool(np.isfinite(repro_max_log) and repro_max_log <= REPRO_TOL_LOG and labels == labels_committed)
    return {"self_check_passed": self_check, "barred_reproduced": repro,
            "noop": bool(self_check and repro and noop(sd_b, sd_perp)), "proceed": bool(self_check and repro)}


def axis_gate(ctx, dec: dict, run_dir: Path, data: Path, d5d6_path: Path, meta_path: Path) -> tuple[dict, np.ndarray]:
    """The gate's numbers: from WT files, the Mab21l2 file and the hashed Mab21l2 factors; no Gata4 file, no score."""
    fc = ctx.rules["components"]["factor"]
    clip = tuple(fc["clip"])
    genes = list(ctx.genes)
    phi_b, rec_b = bs.load_factor(run_dir, BARRED, genes, dec)
    phi_95, _ = bs.load_factor(run_dir, "transfer_mab21l2_vs95", genes, dec)
    mu, cnt = ladder_means(ctx, data / WT_FILES["E8.5"])
    d1 = stage_axis(mu, cnt, ctx.V, fc, "E9.5", "E8.25_late")
    d2 = stage_axis(mu, cnt, ctx.V, fc, "E9.5", "E8.5")
    f_b, f_95 = centred_log(phi_b), centred_log(phi_95)
    committed = json.loads(Path(d5d6_path).read_text())["D5"]
    corr_fm_d = float(np.corrcoef(f_b, d1)[0, 1])
    corr_committed = float(committed["primary"]["corr_fM_d"])
    both, phi_re, den, key = barred_denominator(ctx, data / KO_FILES["mab21l2"])
    repro_max = float(np.abs(np.log(phi_re) - np.log(phi_b)).max())
    labels_committed = int(json.loads(Path(meta_path).read_text())["transfer_ko"]["labels_used_all"])
    res = remove_axis(phi_b, d1, den, clip)
    phi_perp = res.pop("phi_perp")
    sd_b, sd_perp = float(np.log(phi_b).std()), float(np.log(phi_perp).std())
    second = {"axis": "E9.5 vs E8.5", "corr_fM_d": float(np.corrcoef(f_b, d2)[0, 1]),
              "corr_fM_d_committed": float(committed["secondary"]["corr_fM_d"]), "scored": False}
    gate = {"corr_fM_d": corr_fm_d, "corr_fM_d_committed": corr_committed, "tolerance": SELF_CHECK_TOL,
            "reported_only": {
                "corr_all_minus_vs95_d": float(np.corrcoef(f_b - f_95, d1)[0, 1]),
                "corr_all_minus_vs95_d_committed": float(committed["context"]["corr_d_primary_with_log_all_over_vs95"]),
                "secondary_axis": second},
            "barred": {"factor": rec_b["name"], "sha256": rec_b["sha256"], "labels_in_mask": int(both.sum()),
                       "labels_in_mask_committed": labels_committed, "max_abs_diff_log_phi": repro_max,
                       "harmoniser_depth_target": key},
            "removal": {**res, "sd_log_phi_B": sd_b, "sd_log_phi_perp": sd_perp, "sd_change": sd_perp - sd_b,
                        "corr_log_phi_perp_with_axis": float(np.corrcoef(centred_log(phi_perp), d1)[0, 1]),
                        "corr_log_phi_perp_with_barred": float(np.corrcoef(centred_log(phi_perp), f_b)[0, 1])},
            **gate_outcome(corr_fm_d, corr_committed, repro_max, int(both.sum()), labels_committed, sd_b, sd_perp)}
    return gate, phi_perp


def axis_text(res: dict) -> str:
    g, r = res["gate"], res["gate"]["removal"]
    ro, sec = g["reported_only"], g["reported_only"]["secondary_axis"]
    verdict = ("STOP and fix the code before any scoring" if not g["proceed"] else
               ("the removal is a no-op: U3 is inconclusive and B-perp is not scored; U3d still runs" if g["noop"]
                else "proceed to score U3 and U3d"))
    lines = ["U3's gate (experiments/harness/TRIES.md; target-free: no Gata4 file, no score)",
             f"self-check: corr(f_M, d) = {g['corr_fM_d']:.4f} against the committed D5 value "
             f"{g['corr_fM_d_committed']:.4f} (tolerance {g['tolerance']}): "
             f"{'passed' if g['self_check_passed'] else 'FAILED'}",
             f"barred factor recomputed from the run's label mask: max |diff log phi| "
             f"{g['barred']['max_abs_diff_log_phi']:.2e}, labels in the mask {g['barred']['labels_in_mask']} "
             f"(committed {g['barred']['labels_in_mask_committed']}): "
             f"{'reproduced' if g['barred_reproduced'] else 'NOT REPRODUCED'}",
             f"removal: b = {r['b']:+.4f}, SD of log phi {r['sd_log_phi_B']:.4f} -> {r['sd_log_phi_perp']:.4f} "
             f"(change {r['sd_change']:+.4f}; no-op below {NOOP_SD}), corr with the axis after the removal "
             f"{r['corr_log_phi_perp_with_axis']:+.4f}, corr with the barred factor "
             f"{r['corr_log_phi_perp_with_barred']:+.4f}, genes at the clip {r['share_of_genes_clipped']:.3f}",
             f"reported only: secondary axis ({sec['axis']}) corr(f_M, d) {sec['corr_fM_d']:.4f} "
             f"(committed {sec['corr_fM_d_committed']:.4f}), not scored; corr(f_M - f_M95, d) "
             f"{ro['corr_all_minus_vs95_d']:.4f} (committed {ro['corr_all_minus_vs95_d_committed']:.4f})",
             f"decision: {verdict}",
             f"B-perp factor: {res['factor_perp']['file']} sha256 {res['factor_perp']['sha256']}",
             f"commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, utc {res['stamp']['utc']}"]
    return "\n".join(lines) + "\n"


def cmd_axis(a) -> None:
    out, fdir = Path(a.out), Path(a.factor_dir)
    json_p, txt_p, csv_p = out / AXIS_JSON, out / AXIS_TXT, fdir / PERP_CSV
    outputs = [json_p, txt_p, csv_p]
    exempt = [*outputs, *(provenance.stamp_path(p) for p in outputs)]
    refuse_existing(outputs, a.rerun_row)
    started = tree_state(exempt)
    rules = bs.load_rules()
    decision = Path(a.decision)
    dec = bs.load_decision(decision, rules)
    hashes = bs.check_hashes(dec, allow_mismatch=False)
    data = Path(a.data)
    run_dir = Path(a.run_dir)
    paths = context_paths(data, with_targets=False)
    with bs.thread_limits(a.threads), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context(dec["preset"], paths, rules, n_jobs=a.threads, threads=a.threads, out=Path(tmp))
        gate, phi_perp = axis_gate(ctx, dec, run_dir, data, Path(a.d5d6), Path(a.meta))
    fdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": list(ctx.genes), "factor": phi_perp}).to_csv(csv_p, index=False)
    provenance.write_stamp(csv_p)
    factor_rec = {"file": _rel(csv_p), "sha256": t3d.sha256(csv_p)}
    factors = [run_dir / "factors" / f"{n}.csv" for n in (BARRED, "transfer_mab21l2_vs95")]
    inputs = input_hashes([data / WT_FILES[k] for k in ("E8.25_late", "E8.5", "E8.75", "E9.5")]
                          + [data / KO_FILES["mab21l2"], decision, Path(a.d5d6), Path(a.meta), *factors])
    res = {"what": "U3's gate", "script_sha256": t3d.sha256(Path(__file__)), "t3_hashes": hashes, "gate": gate,
           "factor_perp": factor_rec, "inputs_sha256": inputs,
           "stamp": run_stamp(started, exempt, a.rerun_row)}
    write_outputs([(json_p, json.dumps(t3d.json_safe(res), indent=1) + "\n"), (txt_p, axis_text(res))])
    print(axis_text(res), end="")
    if not gate["proceed"]:
        sys.exit("orth_diag axis: the gate failed; stop and fix the code before any scoring")


def load_gate(out: Path, script_sha: str) -> dict:
    p = out / AXIS_JSON
    if not p.is_file():
        sys.exit(f"orth_diag: {p} not found: run `orth_diag.py axis` first (the gate comes before any scoring)")
    gate = json.loads(p.read_text())
    if not gate["gate"]["proceed"]:
        sys.exit(f"orth_diag: the gate in {p} did not pass; no scoring (see {AXIS_TXT})")
    if gate["script_sha256"] != script_sha:
        sys.exit(f"orth_diag: {p} was written by another version of this script; rerun the gate first")
    return gate


def load_phi_perp(gate: dict, genes: list[str]) -> np.ndarray | None:
    """B-perp's factor, checked against the gate's hash; None when the removal was a no-op (B-perp is not scored)."""
    if gate["gate"]["noop"]:
        return None
    p = REPO / gate["factor_perp"]["file"]  # the gate's path: relative to the repository, or absolute outside it
    h = t3d.sha256(p)
    if h != gate["factor_perp"]["sha256"]:
        sys.exit(f"orth_diag: {p}: sha256 {h[:16]} differs from the gate's ({gate['factor_perp']['sha256'][:16]})")
    df = pd.read_csv(p, float_precision="round_trip")
    if list(df["gene"].astype(str)) != list(genes):
        sys.exit(f"orth_diag: {p}: genes differ from the WT reference's genes (order included)")
    return df["factor"].to_numpy(np.float64)


# ---------------------------------------------------------------------------------------------------------------------
# U3d: the depth-matched world


def matched_counts(C, lab, ratio, props, rng) -> np.ndarray:
    """Counts of cells with their label's depth ratio (1.0 for a cell outside the vocabulary, lab -1) applied with
    t3_decompose.rescale_depth (binomial thinning below 1, Poisson top-up from `props` above 1)."""
    r = np.ones(len(lab), dtype=np.float64)
    known = lab >= 0
    r[known] = np.asarray(ratio, dtype=np.float64)[lab[known]]
    return t3d.rescale_depth(C, r, props, rng)


def reference_counts(ctx, wt_path, label_mode: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """The reference as released (the rows of ctx.ref_X, in file order): dense counts (interpolated records rounded, as
    t3_decompose.dense_counts does), each row's label index (-1 outside the vocabulary) and kNN-smoothed proportions for
    the top-up of a label whose knockout is deeper than WT. Checked against ctx.ref_X before use: over the first 5,000
    rows holding integer counts, the largest gap between the log1p CP10k of the counts and the released expression must
    be about 0 (the file's own convention; only the rounded interpolated records can differ)."""
    from vec.prep import load_stage

    ref_ad = load_stage(wt_path, genes=ctx.genes, drop_duplicates=False)
    if ref_ad.n_obs != len(ctx.ref_X):
        sys.exit("orth_diag: the reference counts do not line up with the context's reference rows")
    t3d.add_label_column(ref_ad, label_mode)
    pos = {t: i for i, t in enumerate(ctx.V.labels)}
    lab = np.array([pos.get(t, -1) for t in ref_ad.obs["t3_label"].astype(str).to_numpy()])
    C = t3d.dense_counts(ref_ad)
    M = sparse.csr_matrix(ref_ad.layers["counts"], dtype=np.float64)
    fractional = np.abs(M.data - np.rint(M.data)) > 1e-6
    has_fraction = np.zeros(M.shape[0], dtype=bool)
    has_fraction[np.repeat(np.arange(M.shape[0]), np.diff(M.indptr))[fractional]] = True
    rows = np.flatnonzero(~has_fraction)[:5000]
    gap = float(np.abs(t3d.lognorm(C[rows]) - np.asarray(ctx.ref_X)[rows]).max()) if rows.size else float("nan")
    if not np.isfinite(gap) or gap > 1e-3:
        sys.exit(f"orth_diag: the released reference's counts do not reproduce its expression (largest gap {gap:.3g}); "
                 "stop and look")
    sm = ctx.rules["components"]["smoothing"]
    props = t3d.smoothed_props(C, sm["n_pcs"], sm["k"])
    return C, lab, props, gap


# ---------------------------------------------------------------------------------------------------------------------
# Scoring one (target, seed), as t3_decompose.score_seed does, for the four variants of this try


_STATE: dict = {}


def _task(task) -> list[dict]:
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=_STATE["threads"]):
        return score_row(_STATE, *task)


def score_row(st: dict, ti: int, seed: int) -> list[dict]:
    """The T3 script's truth half, ceiling half, reference working copy and n_pred for (target ti, seed), its per-seed
    estimate (for the replicate's per-label depth ratios), then B and B-perp on the committed stream against the
    reference as released, and B and s1's stand-in on the depth-matched pool against the matched reference. Anchors
    are computed as score_seed computes them, once per reference. Records hold score_seed's fields, with the world."""
    from common.core_metrics import SEVERITY_LOG_WORST, de_genes

    from vec.boards import RANKED_METRICS
    from vec.evaluate import Stage, aggregate, panel_metrics, skills_from_anchors, split_half

    ctx, phi_b, phi_perp = st["ctx"], st["phi_b"], st["phi_perp"]
    T, pr, V = ctx.targets[ti], ctx.rules["protocol"], ctx.V
    t0 = time.time()
    rows_t = t3d.working_copy_rows(len(T.X_rel), pr["working_copy_fraction"], seed)
    rows_r = t3d.working_copy_rows(len(ctx.ref_X), pr["working_copy_fraction"], seed + 1)
    truth, other = split_half(Stage(X=T.X_rel[rows_t]), seed=pr["split_seed"])
    ref = Stage(X=ctx.ref_X[rows_r])
    n_pred = pr["prediction_cells_multiplier"] * max(truth.n, ref.n)
    pool_rows = np.flatnonzero(~np.isin(T.hash_pool, T.hash_rel[rows_t]))
    est = t3d.estimate(ctx, T, pool_rows, t3d.stable_rng(seed, f"estimate:{T.name}"))
    floor_cut = SEVERITY_LOG_WORST + pr["severity_floor_tolerance"]

    def anchors_of(reference) -> dict:
        floor_v = panel_metrics("T3", reference, truth, reference, seed=seed)
        ceil_v = panel_metrics("T3", other, truth, reference, seed=seed)
        a = {m: {"floor": floor_v.get(m), "ceiling": ceil_v.get(m)} for m in RANKED_METRICS["T3"]}
        return {m: x for m, x in a.items() if x["floor"] is not None and x["ceiling"] is not None}

    def record(world: str, name: str, X, reference, anchors: dict, up, dn, extra: dict | None = None) -> dict:
        vals = panel_metrics("T3", Stage(X=X), truth, reference, seed=seed)
        skills = skills_from_anchors("T3", vals, anchors)
        score, groups = aggregate("T3", skills)
        sev = vals.get("severity_slope")
        rec = {"target": T.name, "seed": seed, "world": world, "variant": name, "n_pred": int(X.shape[0]), "t3": score,
               **{f"skill_{m}": s for m, s in skills.items()}, **{f"group_{g}": s for g, s in groups.items()},
               **{k: vals.get(k) for k in ("de_score", "de_direction", "severity_slope", "mmd_u", "variogram",
                                           "energy_distance", "pb_rel_err", "pseudobulk_pearson", "_de_chance",
                                           "_slope_r2")},
               "sev_floor": bool(sev is None or sev <= floor_cut), **(extra or {})}
        rec.update(t3d._diagnostics(X, truth.X, reference.X, up, dn))
        return rec

    def stream(name: str):
        return t3d.stable_rng(seed, name)

    # the reference as released (the committed rows' world)
    anchors = anchors_of(ref)
    up, dn, _ = de_genes(truth.X, ref.X)
    X_b = t3d.build(ctx, n_pred, stream(B_STREAM), f=V.p, phi=phi_b)
    recs = [record("as_released", V_B, X_b, ref, anchors, up, dn)]
    if phi_perp is not None:  # U3: B-perp, built on B's stream (the same cells); not scored after a no-op removal
        X_perp = t3d.build(ctx, n_pred, stream(B_STREAM), f=V.p, phi=phi_perp)
        recs.append(record("as_released", V_PERP, X_perp, ref, anchors, up, dn))
    # the depth-matched world: the WT reference working copy and the WT pool at the replicate's per-label depth
    wt = ctx.wt.pool
    C_ref = matched_counts(st["ref_counts"][rows_r], st["ref_lab"][rows_r], est.ratio, st["ref_props"][rows_r],
                           stream(f"dm_ref:{T.name}"))
    C_pool = matched_counts(wt.C, wt.lab, est.ratio, wt.P, stream(f"dm_pool:{T.name}"))
    ref_m = Stage(X=t3d.lognorm(C_ref))
    src = t3d.Pool(C_pool, wt.lab, wt.P, V.L)
    anchors_m = anchors_of(ref_m)
    up_m, dn_m, _ = de_genes(truth.X, ref_m.X)
    extra = {"dm_ref_depth_kept": float(C_ref.sum() / max(st["ref_counts"][rows_r].sum(), 1)),
             "dm_pool_depth_kept": float(C_pool.sum() / max(wt.C.sum(), 1))}
    S1 = ctx.targets[1 - ti].full
    X_bm = t3d.build(ctx, n_pred, stream(B_STREAM), f=V.p, src=src, phi=phi_b)
    X_s1 = t3d.build(ctx, n_pred, stream(S1_STREAM), f=V.p, src=src, phi=S1.phi_nm, ratio=S1.ratio)
    recs.append(record("depth_matched", V_B_DM, X_bm, ref_m, anchors_m, up_m, dn_m, extra))
    recs.append(record("depth_matched", V_S1_DM, X_s1, ref_m, anchors_m, up_m, dn_m, extra))
    print(f"[{T.name} seed {seed}] U3 and U3d: {len(recs)} variants, {time.time() - t0:.0f}s", flush=True)
    return recs


# ---------------------------------------------------------------------------------------------------------------------
# The decision of U3's row


def committed_standins(path, check: bool = True) -> dict:
    """The committed (replicate, seed) rows of the barred factor and of `vs95` and their means, checked (20 rows, the
    means the row states) unless `check` is False (a code check on synthetic files)."""
    rec = pd.read_csv(path)
    out = {}
    for v, stated in COMMITTED_MEANS.items():
        rows = rec[rec["variant"] == v][["target", "seed", "t3", "sev_floor"]].sort_values(["target", "seed"])
        if check and len(rows) != 20:
            sys.exit(f"orth_diag: {v}: {len(rows)} committed rows, expected 20")
        mean = float(rows["t3"].mean())
        if check and abs(mean - stated) > MEANS_TOL:
            sys.exit(f"orth_diag: {v}: committed mean {mean:.4f} differs from the {stated} the row states")
        out[v] = {"rows": rows.reset_index(drop=True), "mean": mean}
    return out


def bool_col(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "1.0"))


def replay_check(recs: pd.DataFrame, committed: dict) -> dict:
    """The rerun of B against the committed rows: the streams match if every row's t3 is within REPLAY_TOL of the
    committed one and every severity-floor flag is the same."""
    rerun = recs[recs["variant"] == V_B][["target", "seed", "t3", "sev_floor"]].rename(
        columns={"t3": "t3_rerun", "sev_floor": "floor_rerun"})
    merged = committed[V_B]["rows"].merge(rerun, on=["target", "seed"], how="inner")
    if merged.empty:
        return {"rows": 0, "matched": False, "note": "no rerun row pairs with a committed row"}
    d = (merged["t3_rerun"] - merged["t3"]).abs()
    floors_equal = bool((bool_col(merged["sev_floor"]) == bool_col(merged["floor_rerun"])).all())
    return {"rows": int(len(merged)), "max_abs_diff_t3": float(d.max()), "floor_flags_equal": floors_equal,
            "matched": bool(len(merged) == len(rerun) and d.max() <= REPLAY_TOL and floors_equal),
            "tolerance": REPLAY_TOL}


def u3_outcome(q_hat: float, floored: int) -> str:
    """U3's rule: (ii) first (a low share or a floored severity reports q <= 0.5), then (i), else (iii)."""
    if not np.isfinite(q_hat):
        return "iii"
    if q_hat <= 0.5 or floored >= 5:
        return "ii"
    if q_hat >= 0.75 and floored <= 2:
        return "i"
    return "iii"


U3_MEANING = {
    "i": "q-hat >= 0.75 and at most 2 of 20 rows at the severity floor: q >= 0.75 reported for N11's grid; the swap is "
         "recommended for the no-answer and zero-answer branches only, subject to the N11 rerun at q-hat "
         "(no_answer_rule); under an explicit half answer N11's number at q-hat is reported and ruling (b) stays. "
         "B-perp is a candidate slot-2 swap file only, by the captain's ruling; it cannot fill the full-bar switch.",
    "ii": "q-hat <= 0.5 or 5 or more of 20 rows at the severity floor: q <= 0.5 reported; ruling (b) stays as it is.",
    "iii": "neither: no change, reported only.",
    "no-op": "the removal is a no-op (SD of log phi changes by less than 0.02): U3 is inconclusive and B-perp was not "
             "scored (no try for U3); U3d's rows were scored; N11 is rerun at the grid's q values only."}


def summarise_u3(recs: pd.DataFrame, committed: dict, expect_rows: int | None = 20) -> dict:
    perp = recs[recs["variant"] == V_PERP].sort_values(["target", "seed"])
    if expect_rows is not None and len(perp) != expect_rows:
        sys.exit(f"orth_diag: {len(perp)} rows of {V_PERP}, expected {expect_rows}")
    mean_b, mean_95 = committed[V_B]["mean"], committed["transfer:mab21l2_vs95"]["mean"]
    if perp.empty:
        return {"variant": V_PERP, "mean_t3": None, "per_replicate": {}, "n_rows": 0, "rows_at_severity_floor": None,
                "committed_B_mean": mean_b, "committed_vs95_mean": mean_95, "q_hat": None,
                "paired_gain_over_B_rerun": None, "outcome": "no-op", "reading": U3_MEANING["no-op"]}
    mean_perp = float(perp["t3"].mean())
    floored = int(bool_col(perp["sev_floor"]).sum())
    q_hat = (mean_perp - mean_95) / (mean_b - mean_95)
    b_rerun = recs[recs["variant"] == V_B].sort_values(["target", "seed"])
    paired = (perp.set_index(["target", "seed"])["t3"] - b_rerun.set_index(["target", "seed"])["t3"]).dropna()
    per_rep = perp.groupby("target")["t3"].mean().to_dict()
    outcome = u3_outcome(q_hat, floored)
    return {"variant": V_PERP, "mean_t3": mean_perp, "per_replicate": {k: float(v) for k, v in per_rep.items()},
            "n_rows": int(len(perp)), "rows_at_severity_floor": floored, "committed_B_mean": mean_b,
            "committed_vs95_mean": mean_95, "q_hat": float(q_hat),
            "paired_gain_over_B_rerun": float(paired.mean()) if len(paired) else None,
            "outcome": outcome, "reading": U3_MEANING[outcome]}


def summarise_u3d(recs: pd.DataFrame) -> dict:
    out = {}
    for v in (V_B_DM, V_S1_DM):
        g = recs[recs["variant"] == v]
        out[v] = {"mean_t3": float(g["t3"].mean()),
                  "per_replicate": {k: float(x) for k, x in g.groupby("target")["t3"].mean().items()},
                  "rows_at_severity_floor": int(bool_col(g["sev_floor"]).sum()), "n_rows": int(len(g)),
                  "mean_dm_ref_depth_kept": float(g["dm_ref_depth_kept"].mean()),
                  "mean_dm_pool_depth_kept": float(g["dm_pool_depth_kept"].mean())}
    return out


def decision_text(res: dict) -> str:
    u3, ud, rep = res["u3"], res["u3d"], res["replay_of_B"]
    head = (f"U3: B-perp {u3['mean_t3']:.3f} over {u3['n_rows']} rows ({u3['rows_at_severity_floor']} at the severity "
            f"floor), committed B {u3['committed_B_mean']:.4f} and vs95 {u3['committed_vs95_mean']:.4f}: q-hat "
            f"{u3['q_hat']:.3f} (a share of score on Gata4, not N11's probability q: B keeps its stage term, so if "
            "the stage does not carry B scores below what B-perp scored; fed to N11 as a stated stand-in)."
            if u3["n_rows"] else "U3: the removal was a no-op; B-perp was not scored and q-hat was not measured.")
    lines = ["U3 and U3d (experiments/harness/TRIES.md): reading of the barred factor without its stage part, and the "
             "two N11 cells measured on a depth-matched Gata4 target",
             head,
             f"  rule outcome ({u3['outcome']}): {u3['reading']}",
             "U3d, the depth-matched world (reference working copy and pool thinned to each replicate's per-label "
             "depth; anchors recomputed):"]
    for v, s in ud.items():
        per = ", ".join(f"{k} {x:.2f}" for k, x in s["per_replicate"].items())
        lines.append(f"  {v}: mean {s['mean_t3']:.3f} ({per}); at the severity floor in "
                     f"{s['rows_at_severity_floor']} of {s['n_rows']} rows; depth kept reference "
                     f"{s['mean_dm_ref_depth_kept']:.3f}, pool {s['mean_dm_pool_depth_kept']:.3f}")
    tries = res["tries"]
    lines += [f"Rerun of B against the committed rows: max |diff t3| {rep.get('max_abs_diff_t3')}, floor flags equal "
              f"{rep.get('floor_flags_equal')}, streams matched: {rep['matched']}. If not matched, B counts as one "
              "more try on the new streams (the row).",
              f"Tries used by this run (TRIES rule 4): U3 {tries['U3']}, U3d {tries['U3d']}, the rerun of B "
              f"{tries['rerun_of_B']}.",
              "N11 with the measured cells: run `orth_diag.py n11`.",
              f"commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, utc {res['stamp']['utc']}"]
    return "\n".join(lines) + "\n"


def cmd_score(a) -> None:
    out = Path(a.out)
    rec_p, dec_p, txt_p = out / RECORDS_CSV, out / DECISION_JSON, out / DECISION_TXT
    outputs = [rec_p, dec_p, txt_p]
    exempt = [*outputs, *(provenance.stamp_path(p) for p in outputs)]
    refuse_existing(outputs, a.rerun_row)
    script_sha = t3d.sha256(Path(__file__))
    gate = load_gate(out, script_sha)  # before anything else: the gate comes before any scoring
    started = tree_state(exempt)
    rules = bs.load_rules()
    decision = Path(a.decision)
    dec = bs.load_decision(decision, rules)
    hashes = bs.check_hashes(dec, allow_mismatch=False)
    data, run_dir, records = Path(a.data), Path(a.run_dir), Path(a.records)
    seeds = list(rules["protocol"]["seeds"]) if a.seeds is None else t3d.parse_seeds(a.seeds)
    registered = seeds == list(rules["protocol"]["seeds"]) and not a.code_check
    committed = committed_standins(records, check=not a.code_check)
    paths = context_paths(data, with_targets=True)
    n_par = a.workers * a.threads
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=Path(tmp))
        t3d.selfcheck(ctx, seeds[0])
        phi_b, _ = bs.load_factor(run_dir, BARRED, list(ctx.genes), dec)
        phi_perp = load_phi_perp(gate, list(ctx.genes))
        ref_C, ref_lab, ref_P, gap = reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        _STATE.update(ctx=ctx, phi_b=phi_b, phi_perp=phi_perp, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P,
                      threads=a.threads)
        tasks = [(ti, s) for ti in range(len(ctx.targets)) for s in seeds]
        rows: list[dict] = []
        if a.workers > 1:
            if not bs.openmp_limited():
                raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
            with mp.get_context("fork").Pool(min(a.workers, len(tasks))) as pool:
                for r in pool.imap_unordered(_task, tasks):
                    rows.extend(r)
        else:
            for t in tasks:
                rows.extend(_task(t))
    recs = pd.DataFrame(rows).sort_values(["world", "variant", "target", "seed"]).reset_index(drop=True)
    replay = replay_check(recs, committed)
    with_perp = phi_perp is not None
    u3 = summarise_u3(recs, committed, expect_rows=None if a.code_check else (20 if with_perp else 0))
    u3d = summarise_u3d(recs)
    inputs = input_hashes([data / WT_FILES[k] for k in ("E8.25_late", "E8.5", "E8.75", "E9.5")]
                          + [data / KO_FILES["rep1"], data / KO_FILES["rep2"], data / KO_FILES["mab21l2"], decision,
                             records, run_dir / "factors" / f"{BARRED}.csv", out / AXIS_JSON,
                             *([REPO / gate["factor_perp"]["file"]] if with_perp else [])])
    tries = {"U3": 1 if with_perp else 0, "U3d": 2, "rerun_of_B": 0 if replay["matched"] else 1}
    res = {"what": "U3 and U3d", "seeds": seeds, "protocol_as_registered": registered,
           "script_sha256": script_sha, "t3_hashes": hashes, "reference_counts_largest_gap": gap,
           "replay_of_B": replay, "u3": u3, "u3d": u3d, "tries": tries, "inputs_sha256": inputs}
    res["stamp"] = run_stamp(started, exempt, a.rerun_row)
    out.mkdir(parents=True, exist_ok=True)
    recs.to_csv(rec_p, index=False)
    provenance.write_stamp(rec_p)
    write_outputs([(dec_p, json.dumps(t3d.json_safe(res), indent=1) + "\n"), (txt_p, decision_text(res))])
    print(decision_text(res), end="")


# ---------------------------------------------------------------------------------------------------------------------
# N11 rerun with U3d's measured cells


def measured_rows(rec: pd.DataFrame, key_order: list[tuple]) -> dict:
    """t3 and severity-floor flag of the two measured variants, in N11's row order ((target, seed) sorted)."""
    out = {}
    for v in (V_B_DM, V_S1_DM):
        g = rec[rec["variant"] == v]
        t3 = {(t, int(s)): float(x) for t, s, x in zip(g["target"], g["seed"], g["t3"], strict=True)}
        fl = {(t, int(s)): bool(x) for t, s, x in zip(g["target"], g["seed"], bool_col(g["sev_floor"]), strict=True)}
        if len(g) != len(key_order) or set(t3) != set(key_order):
            sys.exit(f"orth_diag: {v}: rows do not pair with N11's (replicate, seed) rows")
        out[v] = {"t3": np.array([t3[k] for k in key_order], dtype=np.float64),
                  "floor": np.array([fl[k] for k in key_order], dtype=bool)}
    return out


def simulate_measured(rows, measured, others, p_u2, p_anchor, q, rng, lambdas=(0.5, 1.0), levels=(90, 95)) -> dict:
    """N11's simulation (t3_backup_sensitivity.simulate) with the world "not Gata4-shallow, shared signal carries"
    (U1, anchor, not U2) changed in two cells: s1 there is U3d's measured `dm:xrep:F+depth` row, in the place of the
    judgement cell (drawn and ignored, so that every other world keeps N11's random numbers), and the barred file is
    U3d's measured `dm:transfer:mab21l2_all` row with probability q, otherwise the committed vs95 row as in N11. Every
    other world is exactly N11's. dP is reported for each lambda and level, as N10's rule reads it."""
    n = E.N_DRAWS
    u1 = rng.random(n) < P.P_U1
    anc = rng.random(n) < p_anchor
    u2 = rng.random(n) < p_u2
    row = rng.integers(0, 20, n)

    def committed(v):
        return rows[v]["t3"][row], rows[v]["floor"][row]

    def got(v):
        return measured[v]["t3"][row], measured[v]["floor"][row]

    def judged(m, pf):
        return rng.normal(m, P.SD_J, n), rng.random(n) < pf

    wc = u1 & anc & ~u2
    s1, f1 = judged(46.0, 0.9)
    judged(50.0, 0.9)  # N11's draw for the cell measured here: consumed, so the other worlds keep N11's numbers
    for mask, (x, f) in ((u1 & anc & u2, committed("xrep:F+depth")), (wc, got(V_S1_DM)),
                         (u1 & ~anc & u2, committed("fac:depth"))):
        s1, f1 = np.where(mask, x, s1), np.where(mask, f, f1)
    s2, f2 = committed("blind:G2_bcat_central_d50")
    bx, bf = judged(47.0, 0.6)
    ok = rng.random(n) < q
    ax, af = committed("transfer:mab21l2_all")
    vx, vf = committed("transfer:mab21l2_vs95")
    mx, mf = got(V_B_DM)
    for mask, (x, f) in ((u1 & anc & ~wc, (ax, af)), (wc, (np.where(ok, mx, vx), np.where(ok, mf, vf))),
                         (u1 & ~anc, (vx, vf))):
        bx, bf = np.where(mask, x, bx), np.where(mask, f, bf)
    idx = {k: rng.choice(o["w"].size, size=n, p=o["w"]) for k, o in others.items()}
    rest_a = sum(E.OVERALL_WEIGHT[k] * o["a"][idx[k]] for k, o in others.items())
    rest_b = sum(E.OVERALL_WEIGHT[k] * o["b"][idx[k]] for k, o in others.items())
    eta = rng.normal(0.0, 3.0, n)
    res = {}
    for reading in ("half", "zero"):
        c = P.FLOOR_COST if reading == "zero" else 0.0
        a1, a2, ab = s1 - c * f1, s2 - c * f2, bx - c * bf
        t_today, t_swap = np.maximum(a1, a2), np.maximum(a1, ab)
        o_today, o_swap = rest_a + t_today, rest_a + t_swap
        s1tot = rest_b + a1
        dp = {}
        for lam in lambdas:
            shift = lam * (s1tot - s1tot.mean()) + eta
            for lv in levels:
                cl = float(np.percentile(o_today, lv))
                dp[f"lambda{lam}_q{lv}"] = float(np.mean(o_swap >= cl + shift) - np.mean(o_today >= cl + shift))
        res[reading] = {"dE": float(t_swap.mean() - t_today.mean()), "dP": dp,
                        "wc_gain": float((t_swap - t_today)[wc].mean()), "p_wc": float(wc.mean())}
    mix = {"dE": (1 - P.P_ZERO) * res["half"]["dE"] + P.P_ZERO * res["zero"]["dE"],
           "dP": {k: (1 - P.P_ZERO) * v + P.P_ZERO * res["zero"]["dP"][k] for k, v in res["half"]["dP"].items()}}
    return {"p_u2": p_u2, "p_anchor": p_anchor, "q_barred": q, "half": res["half"], "zero": res["zero"], "mixture": mix}


def n10_bar_met(cell: dict) -> bool:
    """N10's rule on one cell: mixture dE_T3 >= +0.5 and mixture dP >= 0 at the 90th and 95th levels for lambda 0.5
    and 1."""
    m = cell["mixture"]
    return bool(m["dE"] >= 0.5 and all(m["dP"][f"lambda{lam}_q{q}"] >= 0 for lam in (0.5, 1.0) for q in (90, 95)))


def no_answer_rule(u3_rule: str, cells_at_q_hat: Sequence[dict]) -> dict:
    """The row's rule for ruling (b)'s no-answer default, fixed before the run: the swap is recommended for the
    no-answer branch only if U3's rule (i) fired and N10's bar is met at q = q-hat with the measured cells at both
    N11 bands; otherwise the default stands (the numbers go to the captain, who decides either way). Under an explicit
    half answer ruling (b) stays and N11's number at q-hat is reported; on a zero answer the swap stands as ruled."""
    bar = bool(cells_at_q_hat) and all(n10_bar_met(c) for c in cells_at_q_hat)
    change = bool(u3_rule == "i" and bar)
    return {"u3_rule": u3_rule, "n10_bar_met_at_both_bands": bar, "change_no_answer_default_to_swap": change}


def n11_text(res: dict) -> str:
    qh = ("no q-hat (U3's removal was a no-op): the grid's q values only" if res["q_hat"] is None else
          f"q-hat {res['q_hat']:.3f} (clipped to {res['q_hat_clipped']:.3f}) as a stated stand-in for N11's "
          "probability q")
    lines = [f"N11 with U3d's measured cells (s1 in the not-shallow world = measured {V_S1_DM}; barred there = "
             f"measured {V_B_DM} with probability q, else committed vs95), {len(res['cells'])} cells, {qh}; "
             f"{res['self_checks']} estimator self-checks.",
             "p_u2 p_anchor | q | dE half / zero / mixture | wc gain half / zero | dP mixture at q90, q95, lambda 0.5 "
             "and 1"]
    for c in res["cells"]:
        m = c["mixture"]
        dp = " ".join(f"{100 * m['dP'][f'lambda{lam}_q{lv}']:+.2f}" for lam in (0.5, 1.0) for lv in (90, 95))
        lines.append(f"{c['p_u2']:.2f} {c['p_anchor']:.2f} | {c['q_barred']:.2f} | {c['half']['dE']:+.2f} / "
                     f"{c['zero']['dE']:+.2f} / {m['dE']:+.2f} | {c['half']['wc_gain']:+.2f} / "
                     f"{c['zero']['wc_gain']:+.2f} | {dp}")
    r = res["no_answer_rule"]
    lines += [f"The row's no-answer rule: U3's outcome {r['u3_rule']}; N10's bar at q-hat with the measured cells at "
              f"both bands: {r['n10_bar_met_at_both_bands']}. Change the no-answer default to the swap: "
              f"{r['change_no_answer_default_to_swap']} (the captain decides; an explicit half answer leaves ruling "
              "(b) as it is; B-perp can only ever be the slot-2 swap file, by his ruling).",
              f"commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, utc {res['stamp']['utc']}"]
    return "\n".join(lines) + "\n"


def cmd_n11(a) -> None:
    out = Path(a.out)
    json_p, txt_p = out / N11_JSON, out / N11_TXT
    outputs = [json_p, txt_p]
    exempt = [*outputs, *(provenance.stamp_path(p) for p in outputs)]
    refuse_existing(outputs, a.rerun_row)
    started = tree_state(exempt)
    dec_p, rec_p = out / DECISION_JSON, out / RECORDS_CSV
    if not dec_p.is_file() or not rec_p.is_file():
        sys.exit(f"orth_diag: {dec_p} and {rec_p} not found: run `orth_diag.py score` first")
    dec = json.loads(dec_p.read_text())
    rec = pd.read_csv(rec_p)
    q_hat = None if dec["u3"]["q_hat"] is None else float(dec["u3"]["q_hat"])  # None after a no-op removal
    q_clip = None if q_hat is None else float(np.clip(q_hat, 0.0, 1.0))
    door, checks = E.Door(REPO), []
    t1, _, _ = E.board_t1(door, checks)
    heart, _, _ = E.board_heart(door, checks)
    embryo = E.board_embryo(door, checks)
    rows = P.load_rows(door)
    measured = measured_rows(rec, rows["transfer:mab21l2_all"]["key"])
    others = {"T1": t1, "heart": heart, "embryo": embryo}
    grid = sorted(set(S.Q_BARRED) | ({q_clip} if q_clip is not None else set()), reverse=True)
    cells = [simulate_measured(rows, measured, others, pu, pa, q, np.random.default_rng(P.SEED))
             for pu, pa in S.BANDS for q in grid]
    at_qhat = [] if q_clip is None else [c for c in cells if c["q_barred"] == q_clip]
    res = {"what": "N11 rerun with U3d's measured cells", "q_hat": q_hat, "q_hat_clipped": q_clip,
           "self_checks": len(checks), "bands": list(S.BANDS), "grid_q": grid, "cells": cells,
           "no_answer_rule": no_answer_rule(dec["u3"]["outcome"], at_qhat),
           "inputs_sha256": {**input_hashes([dec_p, rec_p]), **door.read},
           "stamp": run_stamp(started, exempt, a.rerun_row)}
    write_outputs([(json_p, json.dumps(t3d.json_safe(res), indent=1) + "\n"), (txt_p, n11_text(res))])
    print(n11_text(res), end="")


# ---------------------------------------------------------------------------------------------------------------------
# Command line


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("--out", default=str(OUT), help="output folder (default: the registered one)")
        p.add_argument("--rerun-row", default=None, help="the experiments/harness/TRIES.md row that registers a "
                                                         "rerun, recorded in every output; without it a command "
                                                         "refuses to overwrite its outputs")

    def inputs(p):
        p.add_argument("--data", default=str(DATA), help="the released data (default: data/raw)")
        p.add_argument("--run-dir", default=str(RUN_DIR), help="the gata4 run with factors/")
        p.add_argument("--decision", default=str(DECISION), help="the run's decision file (factor sha256s)")

    pa = sub.add_parser("axis", help="the gate: the axis's self-check, the barred factor's reproduction and B-perp, "
                                     "target-free (no Gata4 file, no score)")
    common(pa)
    inputs(pa)
    pa.add_argument("--factor-dir", default=str(FACTOR_DIR), help="where B-perp's per-gene factor goes (runs/, not "
                                                                  "results/)")
    pa.add_argument("--d5d6", default=str(D5D6), help="the committed D5/D6 results")
    pa.add_argument("--meta", default=str(META), help="the gata4 run's committed meta (labels in the barred factor's "
                                                      "mask)")
    pa.add_argument("--threads", type=int, default=16, help="BLAS threads (the T3 run's 16 in the context)")
    ps = sub.add_parser("score", help="the counted run: U3 and U3d on the Gata4 replicates (the gate first)")
    common(ps)
    inputs(ps)
    ps.add_argument("--records", default=str(RECORDS), help="the committed gata4 records (rows of B and vs95)")
    ps.add_argument("--workers", type=int, default=8)
    ps.add_argument("--threads", type=int, default=2, help="BLAS threads per worker (the T3 run used 8 x 2)")
    ps.add_argument("--seeds", default=None, help="override the protocol's seeds (a code check only)")
    ps.add_argument("--code-check", action="store_true", help="a code check on synthetic files: the committed rows "
                                                              "are not held to the row's numbers (never for the run)")
    pn = sub.add_parser("n11", help="N11 rerun with U3d's measured cells, committed results only (after score)")
    common(pn)
    a = ap.parse_args(argv)
    {"axis": cmd_axis, "score": cmd_score, "n11": cmd_n11}[a.command](a)


if __name__ == "__main__":
    main()
