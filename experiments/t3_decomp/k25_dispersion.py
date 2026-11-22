#!/usr/bin/env python
"""Try K25 (experiments/harness/TRIES.md; the row is committed before any run, and every subcommand refuses without
it): M-B, carried within-label dispersion. The first T3 arm that changes cell-to-cell spread and leaves the mean. A
labelled reading for a ruling, not the adoption rule. The specification is the score-blind judge's text for M-B; its
fixed statements (the gate, the arms, the control, the pass line, the void conditions, the try count) are implemented
here as written; where it was unclear the narrowest reading is taken and listed in the hand-back.

Estimator (every constant fixed in the row; `estimate_ratio`, a pure function of counts and labels):
  for each label with at least MIN_LABEL_CELLS (50) cells in both source files (the knockout's cells called into the
  wild type's vocabulary by the frozen harmoniser, littermates split off, as vec.harness.t3_source_frame does it; the
  wild type's own labels in the source preset's label mode), y = log1p(1e4 c / N) with N each cell's total count, each
  gene residualised on ln N within the label (OLS with an intercept, in each file separately); v_KO and v_WT the
  residual variances (sum of squares over n - 2); ln q_Lg = ln(v_KO / v_WT) where both are positive.
  Per gene, ln r_g = the knockout-cell-count-weighted median of ln q_Lg over the labels, times w_g = n_g / (n_g +
  SHRINK) (SHRINK 200), n_g the knockout's cells in those labels where the gene is detected (count > 0) in at least
  DETECT (10%) of the label's knockout cells; r_g clipped to [0.5, 2]. A gene with no label holds ln r_g = 0 exactly.
  Sources: g1, g2 (both worlds): the Mab21l2 knockout against WT E9.5 (preset mab21l2). x2: the two Gata4 replicates
  pooled against WT E8.75 (preset gata4). A fold's truth never enters its own ratio.

Apply (`stretch`): after the base's factor step, on the stepped counts the base's own build made (captured, never
recomputed: `captured_counts`), for each label L and gene g, c' = m + sqrt(r_g) (c - m), m the label-gene mean of the
stepped counts; negatives clipped to 0 and the clipped mass removed proportionally from the positive cells of the same
label and gene, so the mean is exact (asserted to 1e-9 per label before rounding); rounded stochastically with
default_rng(ROUND_SEED = 20261122). A gene with r_g = 1 is not touched, so r = 1 returns the base byte for byte.

Cells: plain g1, plain g2, dm g1, dm g2 exactly as K18 and K21 to K24 define them (k18's worlds, tasks, stream and
drawn rows; base barred_form, K18's arm value for value, reproduced against K18's committed rows within K17's tolerance
as a guard), and x2 (K9's fold "T3/Gata4->Mab21l2", plain, base K9's `raw`, reproduced against K9's committed x2 rows
within the same tolerance as a guard). K9's module (experiments/t3_decomp/k9_residual.py) is the file of its campaign
branch byte for byte (sha256 pinned below). Arms: `disp` (the try) and `disp_perm` (the counted control: r permuted
across genes within deciles of the reference's mean expression, vec.harness.t3_decile_permutation, one permutation,
seed PERM_SEED = 20261113, K21's). Inside one process k18's ARMS, build_arm and fold_factors, and k9's arm_built, are
replaced by this try's and restored on exit, as K21 to K24 do; no existing file changes.

  pregate --commit C --out DIR           the source-only pre-gate: DIR/pregate.json (no fold target read, no try);
                                         refused first (exit 6, no pregate.json) when a source's ratios miss a gene
                                         of a reference they are mapped into
  run --commit C --world plain --task g1|g2 [--seeds 0-4] --out DIR    one plain part, seed files
  run --commit C --world dm [--workers 8 --threads 2] --out DIR        the depth-matched world, seed files
  run --commit C --world x2 [--seeds 0-4] --out DIR                    K9's x2 fold, seed files
  read --commit C --runs DIR --out DIR   the fixed reading (never over an existing one: exit 4)

Exit codes: 0 done; 1 a crash (masked traceback; seed files written so far stay); 2 refused (no row, a dirty tree,
HEAD not C, a pinned file off its sha256, pregate.json missing, failing or of another commit or machine, a pregate.json
already there); 3 read: pregate.json missing (not a reading, nothing counted); 4 read: the reading exists already;
5 pregate: the pre-gate fails (nothing runs, nothing counts); 6 pregate: a source's ratios do not cover every gene of
the references they are mapped into (the Mab21l2 source: WT E8.75's genes as the dm world loads them and the panel
genes of g1's and g2's boards; the Gata4 source: the panel genes of x2's board): no pregate.json, nothing runs, nothing
counts.
The script prints step lines only ("K25 <part> seed <s>: scored", "wrote <path>", "K25 pregate: passes|fails",
"K25 read: wrote <path>"), never a figure or a count.

Tries (the row): 2 (disp and disp_perm), each counted once across the five cells; a VOID reading counts 2; none when
the pre-gate fails.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import tempfile
import traceback
from contextlib import contextmanager
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (_paired, the series panel)
import k8_s2_form as k8  # noqa: E402  (the upload's central composition, for the record)
import k9_residual as k9  # noqa: E402  (K9's folds, base arm and build; its campaign branch's file, unchanged)
import k18_built_files as k18  # noqa: E402  (the worlds, tasks, stream, seeds, build path and scoring; unchanged)
import k21_factor_forms as k21  # noqa: E402  (K18's originals, the seed parser, gated count, K18's committed rows)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import GROUPS  # noqa: E402

SCRIPT = "experiments/t3_decomp/k25_dispersion.py"
ROW_LABEL = "K25"
ROW_SCRIPT = f"`{SCRIPT}`"
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"

# The estimator's constants (the row)
MIN_LABEL_CELLS = 50
TARGET_SUM = 1e4
SHRINK = 200.0
DETECT = 0.10
R_CLIP = (0.5, 2.0)
# The pre-gate (the row)
PREGATE_MIN_GENES = 50
PREGATE_LN = math.log(1.2)
# The apply step and the control
ROUND_SEED = 20261122
MEAN_TOL = 1e-9
PERM_SEED = k7.PERM_SEED  # 20261113, K21's
PERM_BINS = harness.K2_PERM_BINS  # deciles
# The reading (the row; the adoption rule's constants)
MIN_GAIN, MIN_AHEAD = k7.MIN_GAIN, k7.MIN_AHEAD  # 0.3 and 8 of 10
MAX_EXTRA_GATED = 1
GUARD_MEAN, GUARD_SEED = k18.GUARD_MEAN, k18.GUARD_SEED  # K17's 0.1 and 0.5
TRIES = 2

ARM, CONTROL = "disp", "disp_perm"
K18_BASE = k18.BARRED  # "barred_form"
X2_BASE = k9.BASE  # "raw"
K18_ARMS = (K18_BASE, ARM, CONTROL)
X2_ARMS = (X2_BASE, ARM, CONTROL)
STREAM = k18.STREAM
TASKS = k18.TASKS  # g1, g2
WORLDS = k18.WORLDS  # plain, dm
X2_TASK = k9.TASKS["x2"]
CELLS = ("plain g1", "plain g2", "dm g1", "dm g2", "x2")
# The only harness tasks this try runs (practice folds; no test board): K18's g1 and g2 and K9's x2
PRACTICE_KEYS = frozenset({"T3w/WT->Gata4rep1", "T3w/WT->Gata4rep2", "T3k9/x2"})


def practice_keys() -> frozenset:
    return frozenset({*(t.key for t in TASKS.values()), X2_TASK.key})

# Sources: (fold of vec.harness.K4_FOLDS whose source this is, the cells it serves)
SOURCES = {"mab21l2": "T3/Mab21l2->Gata4rep2", "gata4": "T3/Gata4->Mab21l2"}
SOURCE_OF_CELL = {"plain g1": "mab21l2", "plain g2": "mab21l2", "dm g1": "mab21l2", "dm g2": "mab21l2",
                  "x2": "gata4"}

# Pinned files (committed; a test pins each)
K9_MODULE_SHA256 = "1c58a0908b20cbb196bd3d1f64a522a7858e1d74b73fb1f73bf743c893af1f1c"  # campaign/k9 30ede1f
K9_X2 = HERE / "results" / "k9_residual" / "tasks" / "x2.json"
K9_X2_SHA256 = "8be12f798aeab27e68c44fd54fe3106ee949d8d2b4e4694bd9b5db0bfb5458ba"
K18_TASKS = k21.K18_TASKS
K18_SHA256 = dict(k21.K18_SHA256)

EXIT_REFUSED, EXIT_INCOMPLETE, EXIT_EXISTS, EXIT_PREGATE_FAILS, EXIT_UNCOVERED = 2, 3, 4, 5, 6
MUTANT_STATE = next(ms for g, _, ms in GROUPS["T3"] if g == "mutant_state")  # [("mmd_u", 0.6), ("variogram", 0.4)]
_DISP: dict = {}  # pregate.json's genes and ratios, loaded before any arm builds


def refuse(msg: str) -> None:
    print(f"refused: {msg}", file=sys.stderr, flush=True)
    sys.exit(EXIT_REFUSED)


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# The estimator (pure functions)


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """The lower weighted median: values sorted (stable), the first whose cumulative weight reaches half the total."""
    v = np.asarray(values, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    if v.size == 0 or w.sum() <= 0:
        raise ValueError("a weighted median needs at least one value with a positive weight")
    o = np.argsort(v, kind="stable")
    c = np.cumsum(w[o])
    return float(v[o][int(np.searchsorted(c, 0.5 * c[-1], side="left"))])


def residual_variance(C: np.ndarray) -> np.ndarray:
    """Per gene, the residual variance of y = log1p(1e4 c / N) on ln N (OLS with an intercept) over the rows of C
    (cells with N > 0), as the sum of squared residuals over n - 2; NaN when fewer than 3 cells."""
    C = np.asarray(C, dtype=np.float64)
    N = C.sum(1)
    C, N = C[N > 0], N[N > 0]
    n = C.shape[0]
    if n < 3:
        return np.full(C.shape[1], np.nan)
    y = np.log1p(TARGET_SUM * C / N[:, None])
    x = np.log(N)
    x = x - x.mean()
    yc = y - y.mean(0)
    sxx = float(x @ x)
    beta = (x @ yc) / sxx if sxx > 0 else np.zeros(C.shape[1])
    res = yc - np.outer(x, beta)
    return (res ** 2).sum(0) / (n - 2)


def estimate_ratio(C_wt: np.ndarray, lab_wt: np.ndarray, C_ko: np.ndarray, lab_ko: np.ndarray, L: int) -> dict:
    """The row's estimator from one source: wild-type and knockout counts (cells x genes, one gene order) and label
    indices into one vocabulary of L labels (-1: no label, left out). Returns ln_r (shrunk, before the clip), r
    (clipped), n_g, the labels used and per-gene flags."""
    C_wt, C_ko = np.asarray(C_wt), np.asarray(C_ko)
    lab_wt, lab_ko = np.asarray(lab_wt), np.asarray(lab_ko)
    G = C_wt.shape[1]
    if C_ko.shape[1] != G:
        raise ValueError("the two source files hold different genes")
    used = [t for t in range(L) if (lab_wt == t).sum() >= MIN_LABEL_CELLS and (lab_ko == t).sum() >= MIN_LABEL_CELLS]
    lnq = np.full((len(used), G), np.nan)
    n_ko = np.zeros(len(used))
    detected = np.zeros((len(used), G), dtype=bool)
    for i, t in enumerate(used):
        ko = C_ko[lab_ko == t]
        v_ko, v_wt = residual_variance(ko), residual_variance(C_wt[lab_wt == t])
        ok = (v_ko > 0) & (v_wt > 0) & np.isfinite(v_ko) & np.isfinite(v_wt)
        lnq[i, ok] = np.log(v_ko[ok] / v_wt[ok])
        n_ko[i] = ko.shape[0]
        detected[i] = (ko > 0).mean(0) >= DETECT
    med = np.zeros(G)
    has = np.zeros(G, dtype=bool)
    for g in range(G):
        ok = np.isfinite(lnq[:, g])
        if ok.any():
            med[g] = weighted_median(lnq[ok, g], n_ko[ok])
            has[g] = True
    n_g = (detected * n_ko[:, None]).sum(0)
    w = n_g / (n_g + SHRINK)
    ln_r = np.where(has, med * w, 0.0)
    lo, hi = math.log(R_CLIP[0]), math.log(R_CLIP[1])
    clipped = (ln_r < lo) | (ln_r > hi)
    r = np.where(ln_r == 0.0, 1.0, np.exp(np.clip(ln_r, lo, hi)))
    return {"ln_r": ln_r, "r": r, "n_g": n_g, "labels_used": used, "label_cells_ko": n_ko.tolist(),
            "has_label": has, "clipped": clipped}


def pregate_figures(est: dict) -> dict:
    """The pre-gate's figures for one source: genes with |ln r_g| > ln 1.2 (after the clip), the n_g distribution
    (quantiles and a binned histogram) and the share clipped."""
    r = np.asarray(est["r"], dtype=np.float64)
    big = int((np.abs(np.log(r)) > PREGATE_LN).sum())
    n_g = np.asarray(est["n_g"], dtype=np.float64)
    edges = [0, 1, 50, 100, 200, 500, 1000, 2000, 5000, 10000, math.inf]
    hist = {f"[{edges[i]}, {edges[i + 1]})": int(((n_g >= edges[i]) & (n_g < edges[i + 1])).sum())
            for i in range(len(edges) - 1)}
    return {"genes": int(r.size), "genes_big": big, "passes": big >= PREGATE_MIN_GENES,
            "share_clipped": float(np.mean(est["clipped"])) if r.size else 0.0,
            "genes_without_label": int((~np.asarray(est["has_label"])).sum()),
            "genes_r_exactly_1": int((r == 1.0).sum()), "labels_used": len(est["labels_used"]),
            "n_g_quantiles": {q: float(np.quantile(n_g, float(q))) for q in ("0", "0.1", "0.25", "0.5", "0.75",
                                                                            "0.9", "1")} if n_g.size else {},
            "n_g_histogram": hist, "sd_ln_r": float(np.std(np.log(r))) if r.size else 0.0}


# ---------------------------------------------------------------------------------------------------------------------
# The apply step (pure functions)


def stretch(C: np.ndarray, lab: np.ndarray, r: np.ndarray, seed: int = ROUND_SEED) -> np.ndarray:
    """c' = m + sqrt(r_g) (c - m) per label and gene (m the label-gene mean of C), negatives clipped to 0 and the
    clipped mass taken proportionally from the positive cells of the same label and gene, the label-gene means asserted
    equal to MEAN_TOL before rounding, then stochastic rounding with default_rng(seed). Genes with r_g == 1 are not
    touched; with every r_g == 1 the input array itself comes back (byte for byte)."""
    C = np.asarray(C)
    lab = np.asarray(lab)
    r = np.asarray(r, dtype=np.float64)
    if r.shape != (C.shape[1],) or lab.shape != (C.shape[0],):
        raise ValueError("the ratio must hold one value per gene and the labels one per cell")
    if not np.all(np.isfinite(r) & (r > 0)):
        raise ValueError("the ratio must be finite and positive")
    cols = np.flatnonzero(r != 1.0)
    if cols.size == 0:
        return C
    s = np.sqrt(r[cols])
    F = C[:, cols].astype(np.float64)
    for t in np.unique(lab):
        rows = np.flatnonzero(lab == t)
        B = F[rows]
        m = B.mean(0)
        new = m + s * (B - m)
        neg = new < 0
        deficit = np.where(neg, -new, 0.0).sum(0)
        new[neg] = 0.0
        tot = new.sum(0)
        take = (deficit > 0) & (tot > 0)
        new[:, take] *= 1.0 - deficit[take] / tot[take]
        if np.max(np.abs(new.mean(0) - m), initial=0.0) > MEAN_TOL:
            raise AssertionError(f"label {t}: a gene mean moved by more than {MEAN_TOL} before rounding")
        F[rows] = new
    rng = np.random.default_rng(seed)
    fl = np.floor(F)
    out = C.copy()
    out[:, cols] = (fl + (rng.random(F.shape) < (F - fl))).astype(C.dtype)
    return out


@contextmanager
def captured_counts(t3d):
    """While it is open, every call of t3d.lognorm records its argument and returns what the original returns, so the
    stepped counts a build makes are read, never recomputed. t3d.lognorm is restored on exit."""
    original = t3d.lognorm
    calls: list = []

    def spy(C):
        calls.append(np.asarray(C))
        return original(C)

    t3d.lognorm = spy
    try:
        yield calls
    finally:
        t3d.lognorm = original


def stretched(bs, base_builder, r: np.ndarray):
    """The base's Built (base_builder() under captured_counts), then its stepped counts stretched by r and log-normed
    by the frozen script's lognorm; draws, labels and source rows unchanged."""
    with captured_counts(bs.t3d) as calls:
        built = base_builder()
    if len(calls) != 1:
        raise RuntimeError(f"the base's build called lognorm {len(calls)} times, not once")
    C = calls[0]
    if not np.array_equal(bs.t3d.lognorm(C), built.X):
        raise RuntimeError("the captured counts are not the base's stepped counts")
    C2 = stretch(C, built.lab, r)
    if C2 is C:
        return built
    return bs.Built(X=bs.t3d.lognorm(C2), lab=built.lab, source=built.source, row=built.row)


def permuted_ratio(r: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """disp_perm's ratio: r permuted across genes within deciles of `basis`, one permutation, PERM_SEED."""
    return harness.t3_decile_permutation(np.asarray(r, dtype=np.float64), basis, PERM_SEED, PERM_BINS)


def ratio_in_order(source: str, genes: list[str]) -> np.ndarray:
    if not _DISP.get(source):
        raise RuntimeError("the ratios are not loaded from pregate.json")
    d = _DISP[source]
    return k18.in_gene_order(np.asarray(d["r"], dtype=np.float64), d["genes"], list(genes))


# ---------------------------------------------------------------------------------------------------------------------
# The arms: K18's worlds (g1, g2) and K9's x2


_K18_ORIGINAL = dict(k21._K18_ORIGINAL)  # k18's own ARMS, build_arm and fold_factors
_K9_ORIGINAL = {"arm_built": k9.arm_built}


def fold_factors(bs, ref, per: dict, phi: np.ndarray, group: np.ndarray, basis: np.ndarray, phi_barred: np.ndarray,
                 spec: dict) -> dict:
    """K18's fold_factors (the original, called unchanged) plus the Mab21l2 source's ratio in the reference's gene
    order and its permutation within deciles of `basis` (the world's basis as k18 passes it)."""
    F = _K18_ORIGINAL["fold_factors"](bs, ref, per, phi, group, basis, phi_barred, spec)
    F["r"] = ratio_in_order("mab21l2", list(ref.genes))
    F["r_perm"] = permuted_ratio(F["r"], basis)
    return F


def build_arm(bs, ref, arm: str, n: int, seed: int, F: dict):
    """barred_form: K18's arm, the original build_arm, unchanged. disp and disp_perm: the same build, then the stepped
    counts stretched by r (or its permutation)."""
    if arm not in K18_ARMS:
        raise RuntimeError(f"unknown arm {arm!r}")
    base = lambda: _K18_ORIGINAL["build_arm"](bs, ref, K18_BASE, n, seed, F)  # noqa: E731
    if arm == K18_BASE:
        return base()
    return stretched(bs, base, F["r"] if arm == ARM else F["r_perm"])


@contextmanager
def k18_arms():
    """For this process only: k18's ARMS, build_arm and fold_factors replaced by this try's, restored on exit (to what
    they were on entry: K18's own)."""
    saved = (k18.ARMS, k18.build_arm, k18.fold_factors)
    k18.ARMS, k18.build_arm, k18.fold_factors = K18_ARMS, build_arm, fold_factors
    k18._FOLD.clear()
    try:
        yield
    finally:
        k18.ARMS, k18.build_arm, k18.fold_factors = saved
        k18._FOLD.clear()


def k18_candidates() -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: k18._build_plain(ctx, arm),
                                 params={"form": arm, "stream": STREAM, "composition": k8.spec()["composition"],
                                         "factor": "the barred s1's", "dispersion": None if arm == K18_BASE else arm,
                                         "round_seed": None if arm == K18_BASE else ROUND_SEED,
                                         "perm_seed": PERM_SEED if arm == CONTROL else None},
                                 requires=lambda t: t.key in k18.KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in K18_ARMS)]


_X2: dict = {}


def _x2_ratios(ctx) -> dict:
    """The Gata4 source's ratio in x2's reference gene order and its permutation within deciles of the reference's
    mean expression (K9's basis), once per loaded fold."""
    pool = ctx.inputs[ctx.task.reference_role]
    if _X2.get("pool") is not pool:
        prep = k9._prepared(ctx)
        r = ratio_in_order("gata4", list(prep["ref"].genes))
        _X2.clear()
        _X2.update(pool=pool, r=r, r_perm=permuted_ratio(r, harness.t3_mean_expression(pool)))
    return _X2


@contextmanager
def k9_stretch(r: np.ndarray):
    """For one build: k9's arm_built replaced by K9's own `raw` build (the original, unchanged) followed by the
    stretch; restored on exit."""
    def arm_built(bs, ref, n, seed, f, arm, phi_raw, rho, group, clip):
        if arm != X2_BASE:
            raise RuntimeError("the stretch follows K9's raw build only")
        return stretched(bs, lambda: _K9_ORIGINAL["arm_built"](bs, ref, n, seed, f, X2_BASE, phi_raw, rho, group,
                                                               clip), r)

    saved = k9.arm_built
    k9.arm_built = arm_built
    try:
        yield
    finally:
        k9.arm_built = saved


def _build_x2(ctx, arm: str):
    if arm == X2_BASE:
        return k9._build(ctx, X2_BASE)
    rr = _x2_ratios(ctx)
    with k9_stretch(rr["r"] if arm == ARM else rr["r_perm"]):
        return k9._build(ctx, X2_BASE)


def x2_candidates() -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build_x2(ctx, arm),
                                 params=harness.K4_PARAMS | {"dispersion": None if arm == X2_BASE else arm,
                                                             "round_seed": None if arm == X2_BASE else ROUND_SEED,
                                                             "perm_seed": PERM_SEED if arm == CONTROL else None},
                                 requires=lambda t: t.key == X2_TASK.key)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in X2_ARMS)]


# ---------------------------------------------------------------------------------------------------------------------
# Preconditions


# What the row must state, word for word (the rules this code fixes and that decide the reading). The script refuses
# under a row that leaves any of them out, so the committed row and the code cannot part.
ROW_MUST_STATE = (
    # the estimator
    "at least 50 cells in both source files",
    "y = log1p(1e4 c / N), N each cell's total count over every gene of the file as loaded",
    "OLS on ln N with an intercept, in each file separately",
    "the sum of squared residuals over n - 2",
    "a label is left out for a gene when either variance is 0",
    "the lower weighted median, weighted by the knockout's cells in the label",
    "n_g the knockout's cells in the qualifying labels where the gene is detected (count > 0) in at least 10% of the "
    "label's knockout cells",
    "w_g = n_g / (n_g + 200)",
    "clipped to [0.5, 2]",
    "a gene with no label left has ln r_g = 0 exactly (r_g = 1)",
    # the apply step
    "c' = m + sqrt(r_g) (c - m)",
    "default_rng(20261122), created fresh for every build",
    "label-gene means equal to 1e-9 before rounding",
    # the pre-gate
    "at least 50 genes with abs(ln r_g) > ln 1.2 in each source",
    "counted over every gene of the source's wild-type file, after the clip",
    "exit 6, launcher 26",
    # the arms, the control, the guards
    "seed 20261113, 10 bins",
    "within K17's tolerance (ten-seed mean within 0.1, every seed within 0.5)",
    # the pass line
    "AHEAD: paired mean > 0.3 with at least 8 of 10 seeds positive; BEHIND: mean < -0.3",
    "gated: r2_gate, beta_gate, no_response",
    "mutant-state skill up: the ten-seed mean of the paired difference in 0.6 mmd_u + 0.4 variogram skill is above 0",
    # the void list and the count
    "made on another pregate.json",
    "two tries",
)


def row_present(path: Path | None = None) -> str:
    """The K25 row naming this script must be in TRIES.md, exactly one row, and it must state every rule of
    ROW_MUST_STATE word for word."""
    path = TRIES_MD if path is None else Path(path)
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.startswith("| ") and ROW_SCRIPT in ln]
    ids = sorted({ln[2:].split(" ", 1)[0] for ln in rows})
    if ids != [ROW_LABEL]:
        refuse(f"row {ROW_LABEL} naming {ROW_SCRIPT} is not (alone) in {Path(path).name} ({ids}); the row is committed "
               "before any run")
    unstated = [p for p in ROW_MUST_STATE if not any(p in ln for ln in rows)]
    if unstated:
        refuse(f"row {ROW_LABEL} does not state the rules the code fixes: {unstated}")
    return ids[0]


def pinned_ok() -> list[str]:
    why = []
    if sha256_file(HERE / "k9_residual.py") != K9_MODULE_SHA256:
        why.append("k9_residual.py is not its campaign branch's file (sha256)")
    if not K9_X2.is_file() or sha256_file(K9_X2) != K9_X2_SHA256:
        why.append("K9's committed x2 task file is missing or off its sha256")
    return why


def preconditions(commit: str) -> dict:
    row_present()
    st = k18.stamp()
    if not st["tree_clean"]:
        refuse("the tree is not clean (tracked files); not started, nothing written")
    if st["commit"] != commit:
        refuse(f"HEAD is {st['commit']}, not the run commit {commit}")
    why = pinned_ok()
    if why:
        refuse("; ".join(why))
    return st


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate (source only; no fold target)


def source_ratio(source: str, data: Path) -> tuple[dict, dict]:
    """One source's ratio from the released files of vec.harness.K4_FOLDS[SOURCES[source]], loaded as the harness
    loads inputs (vec.prep.load_stage, duplicates dropped; the knockouts in the wild type's genes): the knockout cells
    and labels through vec.harness.t3_source_frame (littermates split off, the frozen harmoniser), the wild type's
    labels in the preset's label mode."""
    from vec.prep import load_stage

    fold = harness.K4_FOLDS[SOURCES[source]]
    bs = harness._t3_head()
    t3d = bs.t3d
    rules = bs.load_rules()
    wt = load_stage(Path(data) / fold["wt_file"], drop_duplicates=True)
    genes = list(wt.var_names)
    kos = [load_stage(Path(data) / f, genes=genes, drop_duplicates=True) for f in fold["source_files"]]
    frame = harness.t3_source_frame(wt, kos, rules, fold["source_preset"])
    w = harness._shallow(wt)
    t3d.add_label_column(w, rules["presets"][fold["source_preset"]]["label_mode"])
    V, wlab = harness._t3_vocab(t3d, w, rules)
    Cw = t3d.dense_counts(w)
    if not np.array_equal(t3d.props_type_means(Cw, wlab, V.L), frame.mu_wt):
        raise RuntimeError("the wild type's labels are not the source frame's")
    est = estimate_ratio(Cw, wlab, frame.C, frame.lab, V.L)
    files = [fold["wt_file"], *fold["source_files"]]
    rec = {"fold": SOURCES[source], "preset": fold["source_preset"], "files_sha256": {
        f: t3d.sha256(Path(data) / f) for f in files}, "labels_used": [str(V.labels[t]) for t in est["labels_used"]]}
    return est, {"genes": genes, "r": [float(x) for x in est["r"]], "record": rec}


def needed_genes(ratios: dict, panel: dict[str, list[str]]) -> dict[str, list[str]]:
    """The genes each source's ratio is mapped into (k18.in_gene_order raises on a missing one): the Mab21l2 source
    serves plain g1 and g2 (their boards' panel genes) and the dm world (WT E8.75's genes as build_context loads them,
    load_stage(drop_duplicates=True), which are the Gata4 source's wild-type genes here); the Gata4 source serves x2
    (its board's panel genes)."""
    mab = [*panel["g1"], *panel["g2"], *ratios["gata4"]["genes"]]
    return {"mab21l2": list(dict.fromkeys(mab)), "gata4": list(dict.fromkeys(panel["x2"]))}


def panel_genes() -> dict[str, list[str]]:
    """The panel genes of g1's, g2's and x2's boards, as vec.harness.run_task loads them (validation boards only)."""
    boards = harness.load_boards()
    return {**{n: list(boards[t.panel_board].genes) for n, t in TASKS.items()},
            "x2": list(boards[X2_TASK.panel_board].genes)}


def uncovered(ratios: dict, need: dict[str, list[str]]) -> list[str]:
    """The sources whose ratio misses a gene it is mapped into."""
    return [s for s in SOURCES if not set(need[s]) <= set(ratios[s]["genes"])]


def cmd_pregate(a) -> None:
    st = preconditions(a.commit)
    out = Path(a.out)
    path = out / "pregate.json"
    if path.exists():
        refuse(f"{path} exists: the pre-gate is made once")
    out.mkdir(parents=True, exist_ok=True)
    try:
        figs, ratios = {}, {}
        for s in SOURCES:
            est, ratios[s] = source_ratio(s, Path(a.data))
            figs[s] = pregate_figures(est) | {"record": ratios[s]["record"]}
        short = uncovered(ratios, needed_genes(ratios, panel_genes()))
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if short:
        print("K25 pregate: refused: a source's ratios miss a gene of a reference they are mapped into; no "
              "pregate.json; nothing runs and nothing counts", flush=True)
        sys.exit(EXIT_UNCOVERED)
    try:
        passes = all(f["passes"] for f in figs.values())
        res = {"row": ROW_LABEL, "passes": passes, "figures": figs, "ratios": ratios,
               "constants": {"min_label_cells": MIN_LABEL_CELLS, "shrink": SHRINK, "detect": DETECT,
                             "r_clip": list(R_CLIP), "pregate_min_genes": PREGATE_MIN_GENES,
                             "pregate_ln": PREGATE_LN, "round_seed": ROUND_SEED, "perm_seed": PERM_SEED,
                             "perm_bins": PERM_BINS}, "tries_if_run": TRIES if passes else 0, "provenance": st}
        k18.write_task_file(path, res)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    print(f"K25 pregate: {'passes' if passes else 'fails; nothing runs and nothing counts'}", flush=True)
    if not passes:
        sys.exit(EXIT_PREGATE_FAILS)


def load_pregate(out: Path, st: dict) -> str:
    """pregate.json of this commit and machine, passing; fills _DISP. Returns its sha256."""
    p = Path(out) / "pregate.json"
    if not p.is_file():
        refuse(f"{p} is missing: the pre-gate runs first")
    raw = p.read_bytes()
    body = json.loads(raw)
    prov = body.get("provenance") or {}
    if prov.get("commit") != st["commit"] or prov.get("machine") != st["machine"]:
        refuse("pregate.json is from another commit or machine")
    if body.get("passes") is not True:
        refuse("the pre-gate failed: nothing runs and nothing counts")
    _DISP.clear()
    _DISP.update({s: {"genes": list(v["genes"]), "r": list(v["r"])} for s, v in body["ratios"].items()})
    return hashlib.sha256(raw).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# Running


def progress_log(part: str):
    def log(line: str) -> None:
        m = re.match(r".* seed (\d+):", str(line))
        print(f"K25 {part} seed {m.group(1) if m else '?'}: scored", flush=True)

    return log


def _check_stamp(res: dict, st: dict) -> None:
    got = res["provenance"]
    if got.get("tree_clean") is not True or got.get("commit") != st["commit"]:
        raise k21._Unclean(f"the harness's stamp says commit {got.get('commit')}, clean {got.get('tree_clean')}: no "
                           "task file written")
    k21._stamp_ok(st)


def _write_seed(path: Path, res: dict, row: dict, st: dict, extra: dict) -> None:
    body = {k: v for k, v in res.items() if k != "rows"} | {"rows": [row], "seed": row["seed"]} | extra
    body["provenance"] = dict(res["provenance"]) | {"machine": st["machine"]}
    prov = k18.write_task_file(path, body)
    print(f"wrote {path}: commit {prov['commit']}, clean {prov['tree_clean']}", flush=True)


def run_plain(a, st: dict, out: Path, pre_sha: str) -> None:
    folder = out / f"plain_{a.task}"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    rec = k18.load_barred_for_plain(Path(a.data), Path(a.run_dir), Path(a.decision))
    with k18_arms():
        for s in k21.parse_seeds(a.seeds):
            path = folder / f"seed_{s}.json"
            if path.exists():
                print(f"K25 plain {a.task} seed {s}: present, skipped", flush=True)
                continue
            res = harness.run_task(TASKS[a.task], k18_candidates(), seeds=[s], data_dir=Path(a.data),
                                   log=progress_log(f"plain {a.task}"))
            _check_stamp(res, st)
            _write_seed(path, res, res["rows"][0], st, {
                "world": "plain", "task_source": SCRIPT, "pregate_sha256": pre_sha, "barred_factor": rec,
                "factor_record": k7._PREP["fold"]["record"]})


def run_x2(a, st: dict, out: Path, pre_sha: str) -> None:
    folder = out / "x2"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    seeds = [s for s in k21.parse_seeds(a.seeds) if not (folder / f"seed_{s}.json").exists()]
    if not seeds:
        print("K25 x2: every seed present, skipped", flush=True)
        return
    srcs = (_DISP_RECORD.get("gata4") or {}).get("files_sha256") or {}
    roles = (*harness.K4_FOLDS[SOURCES["gata4"]]["ko_roles"], "wt_train")
    if sorted(srcs) != sorted(X2_TASK.inputs[r] for r in roles):
        raise RuntimeError("x2's source inputs are not the Gata4 source's files of the pre-gate")
    bs = harness._t3_head()
    for f, h in srcs.items():
        if bs.t3d.sha256(Path(a.data) / f) != h:
            raise RuntimeError(f"{f} is not the file the pre-gate read")
    res = harness.run_task(X2_TASK, x2_candidates(), seeds=seeds, data_dir=Path(a.data), log=progress_log("x2"))
    _check_stamp(res, st)
    for row in res["rows"]:
        _write_seed(folder / f"seed_{row['seed']}.json", res, row, st, {
            "world": "x2", "task_source": SCRIPT, "pregate_sha256": pre_sha,
            "factor_record": k9._PREP["fold"]["record"]})


_DISP_RECORD: dict = {}


def run_dm(a, st: dict, out: Path, pre_sha: str) -> None:
    """The depth-matched world set up as k21.run_dm sets it up, each (replicate, seed) row through k18.dm_row."""
    import multiprocessing as mp

    import orth_diag as od

    folder = out / "dm"
    folder.mkdir(parents=True, exist_ok=True)
    k18.check_folder(folder, st)
    seeds = k21.parse_seeds(a.seeds)
    bs = harness._t3_head()
    t3d = bs.t3d
    spec, dec, rules = k18.barred_spec(Path(a.decision))
    if list(rules["protocol"]["seeds"]) != list(SEEDS):
        raise RuntimeError("the frozen script's seeds are not vec.adoption.SEEDS")
    paths = od.context_paths(Path(a.data), with_targets=True)
    n_par = a.workers * a.threads
    by_name = {v: k for k, v in k18.DM_TARGETS.items()}
    with bs.thread_limits(n_par), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=Path(tmp))
        if sorted(T.name for T in ctx.targets) != sorted(k18.DM_TARGETS.values()):
            raise RuntimeError(f"targets {[T.name for T in ctx.targets]}, expected {sorted(k18.DM_TARGETS.values())}")
        k18.quiet_selfcheck(t3d, ctx)
        genes = list(ctx.genes)
        phi_b, frec = k18.barred_factor(Path(a.run_dir), dec, spec, genes)
        ref_C, ref_lab, ref_P, gap = od.reference_counts(ctx, paths["wt"], ctx.preset["label_mode"])
        wt_path = Path(paths["wt"])
        series = bs.wt_revert_series(wt_path, {n: Path(a.data) / "T2" / f"{n}.h5ad" for n in k7.FULL_PANEL
                                               if n != wt_path.stem})
        phi, per, group, rrec = bs.wt_revert_factors(wt_path, series, genes, ctx.wt.pool.C)
        refw = bs.Reference(preset=ctx.preset_name, rules=rules, genes=genes, V=ctx.V, pool=ctx.wt.pool,
                            coords=np.zeros((ctx.wt.pool.C.shape[0], 3)))
        with k18_arms():
            F = k18.fold_factors(bs, refw, per, phi, group, np.asarray(ctx.wt.X, dtype=np.float64).mean(0), phi_b,
                                 spec)
            k18._STATE.clear()
            k18._STATE.update(ctx=ctx, phi_b=phi_b, ref_counts=ref_C, ref_lab=ref_lab, ref_props=ref_P, F=F,
                              threads=a.threads)
            todo = [(ti, s) for ti, T in enumerate(ctx.targets) for s in seeds
                    if not (folder / f"{by_name[T.name]}_seed_{s}.json").exists()]

            def write(row: dict) -> None:
                k21._stamp_ok(st)
                path = folder / f"{by_name[row['target']]}_seed_{row['seed']}.json"
                k18.write_task_file(path, row | {"world": "dm", "task": by_name[row["target"]],
                                                 "arms": list(K18_ARMS), "pregate_sha256": pre_sha,
                                                 "barred_factor": frec, "reference_counts_largest_gap": gap,
                                                 "provenance": st, "task_source": SCRIPT})
                print(f"wrote {path}", flush=True)

            if a.workers > 1 and len(todo) > 1:
                if not bs.openmp_limited():
                    raise RuntimeError("the OpenMP runtime is not limited to one thread: forked workers would deadlock")
                with mp.get_context("fork").Pool(min(a.workers, len(todo))) as p:
                    for row in p.imap_unordered(k18._dm_task, todo):
                        write(row)
            else:
                for t in todo:
                    write(k18._dm_task(t))


def cmd_run(a) -> None:
    if a.world == "plain" and a.task not in TASKS:
        refuse("--world plain needs --task g1 or g2")
    if practice_keys() != PRACTICE_KEYS:
        refuse("a task of this try is not one of its five practice cells")
    st = preconditions(a.commit)
    out = Path(a.out)
    pre_sha = load_pregate(out, st)
    body = json.loads((out / "pregate.json").read_text())
    _DISP_RECORD.clear()
    _DISP_RECORD.update({s: v["record"] for s, v in body["ratios"].items()})
    try:
        {"plain": run_plain, "dm": run_dm, "x2": run_x2}[a.world](a, st, out, pre_sha)
    except Exception:  # noqa: BLE001  (any crash: masked traceback, exit 1; seed files written so far stay)
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def expected_files() -> list[str]:
    return ([f"plain_{n}/seed_{s}.json" for n in TASKS for s in SEEDS]
            + [f"dm/{n}_seed_{s}.json" for n in TASKS for s in SEEDS] + [f"x2/seed_{s}.json" for s in SEEDS])


def cell_of(f: str) -> tuple[str, int]:
    if f.startswith("plain_"):
        return f"plain {f[6:8]}", int(f.rsplit("_", 1)[1][:-5])
    if f.startswith("dm/"):
        return f"dm {f[3:5]}", int(f.rsplit("_", 1)[1][:-5])
    return "x2", int(f.rsplit("_", 1)[1][:-5])


def base_of(cell: str) -> str:
    return X2_BASE if cell == "x2" else K18_BASE


def group_skill(skills: dict) -> float:
    return float(sum(w * skills[m] for m, w in MUTANT_STATE))


def summarise(rows: list[dict], cell: str) -> dict:
    """Per arm: mean and sd, mean skills, severity states, seeds gated, the mutant-state group skill. Paired
    (k7._paired): disp and disp_perm against the base, disp against disp_perm, the base against the floor where it is
    there. Per term: each arm's mean per-seed skill difference against the base. The mutant-state group: each arm's
    mean per-seed difference against the base."""
    base = base_of(cell)
    arms = (base, ARM, CONTROL)
    with_floor = all(harness.FLOOR in q["scores"] for q in rows)
    out: dict = {"base": base, "arms": {}, "paired": {}, "per_term": {}, "mutant_state_diff": {}}
    for x in ((harness.FLOOR,) if with_floor else ()) + arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "mutant_state": float(np.mean([group_skill(q["skills"][x]) for q in rows])),
            "severity_states": {v: sum(z["state"] == v for z in sev) for v in sorted({z["state"] for z in sev})},
            "severity_gated": k21._gated(sev)}
    pairs = [(ARM, base), (CONTROL, base), (ARM, CONTROL)] + ([(base, harness.FLOOR)] if with_floor else [])
    out["paired"] = {f"{x} - {y}": k7._paired(rows, x, y) for x, y in pairs}
    for x in (ARM, CONTROL):
        out["per_term"][x] = {m: float(np.mean([q["skills"][x][m] - q["skills"][base][m] for q in rows]))
                              for m in rows[0]["skills"][base]}
        d = [group_skill(q["skills"][x]) - group_skill(q["skills"][base]) for q in rows]
        out["mutant_state_diff"][x] = {"mean": float(np.mean(d)), "up_seeds": int(sum(v > 0 for v in d))}
    return out


def _ahead(v: dict) -> bool:
    return bool(v["mean"] > MIN_GAIN and v["ahead"] >= MIN_AHEAD)


def decide(per: dict) -> dict:
    """The pass line fixed in the row, on summarise()'s output per cell (all five present)."""
    d = {c: per[c]["paired"][f"{ARM} - {base_of(c)}"] for c in CELLS}
    ahead = {c: _ahead(v) for c, v in d.items()}
    behind = {c: bool(v["mean"] < -MIN_GAIN) for c, v in d.items()}
    gated = {c: (per[c]["arms"][ARM]["severity_gated"], per[c]["arms"][base_of(c)]["severity_gated"]) for c in CELLS}
    up = [c for c in CELLS if ahead[c]]
    a = bool(ahead["x2"] and any(ahead[c] for c in ("dm g1", "dm g2"))
             and any(ahead[c] for c in ("plain g1", "plain g2")))
    b = not any(behind.values())
    cc = all(g <= g0 + MAX_EXTRA_GATED for g, g0 in gated.values())
    dd = all(_ahead(per[c]["paired"][f"{ARM} - {CONTROL}"]) for c in up)
    ee = all(per[c]["mutant_state_diff"][ARM]["mean"] > 0 for c in up)
    passes = bool(a and b and cc and dd and ee)
    return {"passes": passes, "a": a, "b": b, "c": cc, "d": dd, "e": ee,
            "ahead_cells": [c for c in CELLS if ahead[c]], "behind_cells": [c for c in CELLS if behind[c]],
            "gated_arm_base": {c: list(v) for c, v in gated.items()},
            "outcome": (PASS_TEXT if passes else FAIL_TEXT)}


PASS_TEXT = ("PASS: the dispersion arm passes K25's kill test; a labelled reading for a ruling, not the adoption rule. "
             "The build path of the specification follows only by the captain's ruling")
FAIL_TEXT = "FAIL: the dispersion arm does not pass K25's kill test; nothing follows"
VOID_TEXT = "VOID: no file; nothing follows; the two tries count"
NOT_A_READING = "not a reading: pregate.json is missing, so nothing ran; nothing counted"


def guards(rows: dict, k18_refs: dict, k9_rows: list[dict]) -> dict:
    """The base reproduced within K17's tolerance: barred_form against K18's committed rows in the four K18 cells,
    raw against K9's committed x2 rows."""
    out = {}
    for c in CELLS:
        ref = k9_rows if c == "x2" else k18_refs.get(tuple(c.split(" ")), [])
        dist = k18.distance(rows.get(c, []), base_of(c), ref, base_of(c))
        out[c] = {"distance": dist, "passes": k18.guard_passes(dist)}
    return out


def k9_committed_rows(path: Path = K9_X2) -> tuple[list[dict], list[str]]:
    if not Path(path).is_file() or sha256_file(path) != K9_X2_SHA256:
        return [], ["K9's committed x2 task file is missing or off its sha256"]
    return json.loads(Path(path).read_text())["rows"], []


def reading(runs: Path, k18_dir: Path = K18_TASKS, k9_x2: Path = K9_X2) -> dict:
    runs = Path(runs)
    if not (runs / "pregate.json").is_file():
        return {"incomplete": True, "void": False, "tries_counted": 0, "outcome": NOT_A_READING}
    pre_raw = (runs / "pregate.json").read_bytes()
    pre = json.loads(pre_raw)
    pre_sha = hashlib.sha256(pre_raw).hexdigest()
    prov = {"pregate.json": pre.get("provenance") or {}}
    common = {"incomplete": False, "pregate": {"passes": pre.get("passes"), "figures": pre.get("figures")},
              "pregate_sha256": pre_sha}
    if pre.get("passes") is not True:
        return common | {"void": False, "tries_counted": 0, "outcome": "the pre-gate failed: nothing ran, nothing "
                         "counts", "worlds": {}}
    why: list[str] = []
    rows: dict = {c: [] for c in CELLS}
    sha: dict = {}
    for f in expected_files():
        p = runs / f
        if not p.is_file():
            why.append(f"missing seed file {f}")
            continue
        r, sha[f] = k21._read_json(p)
        prov[f] = r.get("provenance") or {}
        cell, seed = cell_of(f)
        if r.get("pregate_sha256") != pre_sha:
            why.append(f"{f}: made on another pregate.json")
        if cell.startswith("dm"):
            if r.get("task") != cell[3:] or r.get("seed") != seed:
                why.append(f"{f}: another task or seed")
                continue
            row = r
        else:
            key = X2_TASK.key if cell == "x2" else TASKS[cell[6:]].key
            if r["task"]["key"] != key or [q["seed"] for q in r["rows"]] != [seed]:
                why.append(f"{f}: another task or seed")
                continue
            row = r["rows"][0]
        want = (base_of(cell), ARM, CONTROL)
        miss = [x for x in want if x not in row.get("scores", {}) or x not in row.get("skills", {})
                or x not in row.get("t3", {})]
        if miss:
            why.append(f"{f}: missing arms {miss}")
            continue
        rows[cell].append(row)
    rows = {c: sorted(v, key=lambda q: q["seed"]) for c, v in rows.items()}
    for f, p in prov.items():
        if p.get("tree_clean") is not True:
            why.append(f"{f}: a dirty tree")
    commits = {p.get("commit") for p in prov.values()}
    machines = {p.get("machine") for p in prov.values()}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(map(str, commits))}")
    if len(machines) > 1:
        why.append(f"results from different machines {sorted(map(str, machines))}")
    saved = dict(k21.K18_SHA256)
    k21.K18_SHA256.clear()
    k21.K18_SHA256.update(K18_SHA256)
    try:
        k18_refs, kwhy = k21.k18_rows(k18_dir)
    finally:
        k21.K18_SHA256.clear()
        k21.K18_SHA256.update(saved)
    k9_rows, k9why = k9_committed_rows(k9_x2)
    why += kwhy + k9why
    per: dict = {}
    for c in CELLS:
        if [q["seed"] for q in rows[c]] != list(SEEDS):
            why.append(f"{c}: not the ten seeds")
            continue
        per[c] = summarise(rows[c], c)
    g = guards(rows, k18_refs, k9_rows)
    for c, v in g.items():
        if not v["passes"]:
            why.append(f"guard: {c}: the base is off its committed rows (limits {GUARD_MEAN} and {GUARD_SEED})")
    dm = {n: rows[f"dm {n}"] for n in TASKS}
    check = k18.dm_check(dm) if all(len(v) == len(SEEDS) for v in dm.values()) else None
    base = common | {"files_sha256": sha, "provenance": prov, "guards": g, "dm_check_u3d_report_only": check,
                     "tries_counted": TRIES, "worlds": per}
    if why:
        return base | {"void": True, "why": why, "outcome": VOID_TEXT}
    return base | {"void": False, "why": [], "commit": next(iter(commits)), "machine": next(iter(machines))} | \
        decide(per)


def _dist(d: dict | None) -> str:
    return k21._dist(d)


def text(res: dict) -> str:
    lines = [f"Try K25: {res['outcome']}."]
    if res.get("incomplete"):
        return "\n".join(lines) + "\n"
    lines.append(f"tries counted: {res['tries_counted']}")
    pg = res.get("pregate") or {}
    for s, f in (pg.get("figures") or {}).items():
        lines.append(f"pre-gate {s}: passes {f['passes']}; genes with |ln r| > ln 1.2: {f['genes_big']} of "
                     f"{f['genes']}; share clipped {f['share_clipped']:.4f}; sd ln r {f['sd_ln_r']:.4f}; labels used "
                     f"{f['labels_used']}; genes without a label {f['genes_without_label']}; n_g quantiles "
                     f"{f['n_g_quantiles']}; n_g histogram {f['n_g_histogram']}")
    if not (res.get("worlds") or res.get("void")):
        return "\n".join(lines) + "\n"
    for c, v in (res.get("guards") or {}).items():
        lines.append(f"guard {c} ({base_of(c)}): passes {v['passes']}; {_dist(v['distance'])}")
    if res["void"]:
        lines.append("VOID: " + "; ".join(res["why"]))
    else:
        lines.append(f"commit {res['commit']}, machine {res['machine']}. (a) {res['a']}, (b) {res['b']}, (c) "
                     f"{res['c']}, (d) {res['d']}, (e) {res['e']}; ahead in {res['ahead_cells'] or 'none'}; behind in "
                     f"{res['behind_cells'] or 'none'}; gated (arm, base) {res['gated_arm_base']}")
    chk = res.get("dm_check_u3d_report_only")
    if chk:
        lines.append("U3d's check arm in the dm world, reported only: " + "; ".join(
            f"{n} {_dist(c['distance'])}" for n, c in chk.items()))
    lines.append("Reported tables, deciding nothing (g1 and g2 are two replicates of one knockout run):")
    for c, v in (res.get("worlds") or {}).items():
        lines.append(f"{c} (base {v['base']}):")
        for a, s in v["arms"].items():
            lines.append(f"    {a:<12} {s['mean']:.2f} (sd {s['sd']:.2f}); gated {s['severity_gated']} "
                         f"{s['severity_states']}; mutant-state group skill {s['mutant_state']:.4f}; skills "
                         + ", ".join(f"{m} {x:.3f}" for m, x in s["skills"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
        for a, terms in v["per_term"].items():
            ms = v["mutant_state_diff"][a]
            lines.append(f"    per-term skill, {a} - base: " + ", ".join(f"{m} {x:+.4f}" for m, x in terms.items())
                         + f"; mutant-state group {ms['mean']:+.4f} (up in {ms['up_seeds']} seeds)")
    if res.get("files_sha256"):
        lines.append("task files: " + ", ".join(f"{f} {h[:12]}" for f, h in res["files_sha256"].items()))
    if res.get("reading_code"):
        rc = res["reading_code"]
        lines.append(f"read by {rc['script']} {rc['script_sha256'][:12]} at commit {rc['commit']}, clean "
                     f"{rc['tree_clean']}, machine {rc['machine']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    if any((out / f).exists() for f in ("reading.json", "reading.txt")):
        print(f"refused: {out} holds a reading already; it is never overwritten", file=sys.stderr, flush=True)
        sys.exit(EXIT_EXISTS)
    st = preconditions(a.commit)
    try:
        res = reading(Path(a.runs))
        if not res.get("incomplete"):
            res["reading_code"] = {"script": SCRIPT, "script_sha256": sha256_file(Path(__file__)),
                                   "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}
        t = text(res)
        js = json.dumps(res, indent=1, default=float)
    except Exception:  # noqa: BLE001
        print(k18.masked(traceback.format_exc()), file=sys.stderr, end="", flush=True)
        sys.exit(1)
    if res.get("incomplete"):
        print("K25 read: not a reading (pregate.json missing)", flush=True)
        sys.exit(EXIT_INCOMPLETE)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(js)
    (out / "reading.txt").write_text(t)
    print(f"K25 read: wrote {out / 'reading.json'}", flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pregate")
    p.add_argument("--out", required=True)
    r = sub.add_parser("run")
    r.add_argument("--world", choices=(*WORLDS, "x2"), required=True)
    r.add_argument("--task", choices=sorted(TASKS))
    r.add_argument("--seeds", default=None)
    r.add_argument("--out", required=True)
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--threads", type=int, default=8)
    r.add_argument("--run-dir", default=str(k18.RUN_DIR))
    r.add_argument("--decision", default=str(k18.DECISION))
    for q in (p, r):
        q.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    for q in (p, r, d):
        q.add_argument("--commit", required=True)
    a = ap.parse_args(argv)
    {"pregate": cmd_pregate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
