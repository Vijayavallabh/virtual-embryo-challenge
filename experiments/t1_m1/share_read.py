#!/usr/bin/env python
"""Try T47 (experiments/harness/TRIES.md): how far the released days move along Feng's direction, by kind of sample.
A labelled diagnostic; it scores nothing and calls no scorer.

M1 applies Feng's E9.5 -> E14.5 change (per group, the log ratio of mean CP10k, feng_dirs.py) in proportion to days:
two fifths at the test, E10.5 to E12.5. Try T44 found that on the released day E9.5 -> E10.5 twice the proportional
step still scores better. Three things can make the scoring step differ from a fifth: the kind of sample (E10.5's
library 3 is a dissected heart, libraries 1 and 2 and the earlier stages embryo portions), the pace (the first day may
hold more than a fifth of the five days' change) and the scale (the same contrast may be larger or smaller in the
challenge's data than in Feng's). This read measures the three with the entry's own formula on the three released T1
files and Feng's statistics at the bracket's two ends only (feng_dirs.bracket_stats).

Definitions (the entry's own wherever one exists):
  groups   M1's direction groups (feng_dirs.directions_from_stats on the bracket); a group's cells at a T1 stage are the
           cells whose label maps to it (feng_dirs.label_groups), else, at a stage without the group, those of its
           config fallback (dirlib.antecedent_labels, as the readout uses it). A cell set under the config's
           min_cells_t1 gives no value.
  genes    per group: Feng's direction finite and the gene not in M1's held list (try T44's: gmom_combo.held_feng).
  change   dirlib.lfc_ratio of the two sets' mean CP10k with the config's pseudocount; d = Feng's direction.
  r_d      per group, the reliability of Feng's direction, computed here from the bracket's two halves
           (feng_dirs.direction_halves on the same bracket statistics: d_a from half A's cells at both ends, d_b from
           half B's) over the share's own gene mask: r_d = max(0, mean(d^2) - mean(((d_a - d_b) / 2)^2)) / mean(d^2)
           over the genes in the mask where d, d_a and d_b are finite (contrast_reliability, N48's js formula; NaN,
           not N48's 1.0, where mean(d^2) = 0 or no gene has the three values). N48's committed 'js_without_held' and
           'js' scales (whose gene set leaves out only the pre-registered held list, not R2's union) are reported
           beside it and enter nothing.
  share    s = sum(change x d) / sum(d x d); beside it the Pearson correlation, the slope with an intercept, the share
           over the half of the genes with the larger |d|, share_corrected = s / r_d (where r_d > 0), the effective
           gene number n_eff = (sum d^2)^2 / sum d^4, the share without the top-|d| genes that carry half of sum d^2,
           and the Pearson correlation over the genes with d != 0 (the share itself is the same over those genes: a
           gene with d = 0 adds nothing).
  combined the groups' shares averaged, weighted by the groups' cells in E10.5's library 3, and unweighted.
  (a) E8.5 -> E9.5; (b) E9.5 -> E10.5 library 1, 2, 3 and pooled; (c) E10.5 libraries 1+2 -> 3 and 1 -> 2;
  (d) per pair of groups, the slope through the origin of the challenge's E9.5 contrast (g against h) on Feng's E9.5
      contrast, over the genes both groups keep; L = the median over the pairs whose Feng label sets do not overlap
      (the median over every pair beside it); the pair's Feng reliability r_c = max(0, mean(c^2) - mean(((c_a -
      c_b) / 2)^2)) / mean(c^2) from Feng's two halves (N48's formula on the contrast), L_corr = the median of slope
      / r_c over the same pairs where r_c > 0; beside it E10.5 library 3 on Feng E14.5. The pairs each of L, L_corr
      and L over every pair rests on (a finite value at the point) are written by name and number.
  calibration  M1's own one-day step (t1_scope.external_shift, fraction 0.2, stream (47, 1)) applied to E10.5
           library 3's cells, read before -> after with the same share code: the share M1's 0.2 produces.
Every value has a bootstrap over the challenge's cells (200 draws, cells resampled within label, stage and library;
Feng's statistics, r_d, r_c and the calibration are fixed) with its 5th and 95th percentiles and its number of
undefined draws.

The reading (fixed in row T47 before the run; written with the numbers): A = combined share of (b) libraries 1 and 2,
averaged; K = combined share of (b) library 3 - A; A_corr and K_corr the same from share_corrected; all four over one
set of groups, those with a share and a share_corrected in all three libraries (the groups each condition drops are
written); P_own_cells = combined share of (a) over the groups whose cells are their own at E8.5 and E9.5 (the reading's
P), and beside it P = the same over every group with an (a) value; L = the scale of (d). A is attenuated by the noise
in Feng's direction and L by the noise in Feng's E9.5 contrast; the pace and the entry-fraction statements are read on
the corrected pair (A_corr, L_corr: the courses in units of Feng's true direction) and on the raw pair, and are "not
separated" where the two disagree (which of the two causes applies is written). Statements on the kind (K against 0.1),
the pace (A / L against 0.2 and 0.3) and the entry's fraction 0.4 against f_max = L - A, f_even = (L - A) / 2 and
f_decay = L (x - x^3) / (1 - x^5) with (1 - x) / (1 - x^5) = A / L (f_decay = 0 where A / L >= 1, no decaying course
where A / L <= 0). "above the bound" is read from f_max alone; a statement whose quantity is undefined in any draw reads
"not computable (n of N draws undefined)", and the corrected pace and fraction statements read "not computable" where
L_corr rests on fewer pairs than L. The courses times the weighted r_d (f_*_corr_noisy_d) are reported only. The read
puts nothing forward by itself; its result can only go to the team captain. The output holds shares, slopes,
correlations, gene numbers and percentiles only: no cell count.

Ledger sentences (printed in t47.txt as notes):
- A is attenuated by the noise in Feng's direction and L by the noise in Feng's E9.5 contrast; the pace and the
  entry-fraction statements are read on the corrected pair and are not separated where raw and corrected disagree.
- A share is a d^2-weighted average dominated by the few genes with the largest Feng change; n_eff is reported beside
  it.
- The corrected fraction statement compares 0.4 with the courses in units of Feng's true direction; the courses times
  r_d (the best multiple of the noisy d, which includes N48's shrinkage) are reported beside it and decide nothing.
- r_d is computed in this read from the bracket's two halves over the share's own genes; N48's scales are reported
  beside it and enter nothing.
- The intervals hold only the challenge's cell sampling; Feng's own noise, of the same size as the cell sampling on the
  reviewer's synthetic check, is not in them, so a separated statement is weaker than its interval reads.
- K counts any difference between the heart library and the embryo libraries that lies along Feng's direction, ambient
  heart RNA included (2% soup gives about 0.13 on synthetic data); a depth difference alone gives none. A K above 0.1
  does not show that dissected-heart cells are further along in development.
- L is the size of between-type contrasts at E9.5 in the challenge's data relative to Feng's; using it for the
  within-type change over five days assumes one scale for both, which the read cannot test; the E10.5 library 3
  against E14.5 slope and each pair's correlation are reported as the only check.
- A share is what the genes with the largest Feng change say; where share and share_without_top_d2_half differ by more
  than their intervals, the statements describe a handful of genes, not the direction.
- Not separated may mean that correcting for Feng's noise changes the statement; the corrected statement is the
  estimate and the raw statement shows how much the correction moves it.
- The courses assume the same E9.5 and E14.5 in both datasets, one scale for all genes and groups and a course that
  never turns back; they are two courses among many, and a statement about 0.4 says nothing about scores (T44 found
  larger steps score better).
- P over own cells is the reading; P with fallback groups compares one group's cells at E8.5 with another group's at
  E9.5 and is shown only beside it.
- No statement of T47 is a reason to change M1's fraction; any change goes as a counted try on the E10.5 cuts with its
  own row.

  uv run --frozen --no-sync python -I -B experiments/t1_m1/share_read.py run --out t47_share_read

Writes experiments/t1_m1/results/<out>/t47.{json,txt}; refused without row T47 in TRIES.md, on an existing output and
on a dirty tree.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Mapping, Sequence
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in ("t1_hygiene", "t1_scope", "t1_direction", "external"):
    sys.path.insert(0, str(REPO / "experiments" / _p))
sys.path.insert(0, str(HERE))
import dirlib as dl  # noqa: E402
import feng_dirs  # noqa: E402
import filterio as fio  # noqa: E402
import g_mom_verify as gv  # noqa: E402
import gmom_combo as gc  # noqa: E402
import hygiene  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402
import t1_stats  # noqa: E402

from vec import t1_scope  # noqa: E402

ROW = "T47"
STEM = "t47"
RESULTS = HERE / "results"
T1_STAGES = ("E8.5", "E9.5", "E10.5")
T1_FILES = {s: REPO / "data" / "raw" / "T1" / f"{s}_RNA.h5ad" for s in T1_STAGES}
GENE_STAGE = "E9.5"  # the T1 gene order (the build's: its start file)
N48_JSON = HERE / "results" / "n48_m1_direction_shrink" / "n48.json"
BRACKET = feng_dirs.BRACKET  # ("E9.5", "E14.5")
N_BOOT = 200
BOOT_SEED = (47, 0)  # row T47's stream: no build seed applies to this read; no existing step draws on it
CAL_STREAM = (47, 1)  # the calibration arm's draw (external_shift's depth restoration)
CAL_FRACTION = 0.2  # M1's one-day step in the E9.5 -> E14.5 bracket (feng_dirs.bracket_fraction(9.5, 10.5))
PERCENTILES = (5.0, 95.0)
LIB_EMBRYO = ("1", "2")  # E10.5's embryo-portion libraries (docs/DATA_NOTES.md)
LIB_HEART = "3"  # E10.5's dissected-heart library
ENTRY_FRACTION = 0.4  # M1's fraction at the test, E10.5 -> E12.5 in the E9.5 -> E14.5 bracket
KIND_AT = 0.1
PACE_FIFTH, PACE_UPPER = 0.2, 0.3

# the comparisons: name -> ((stage, libraries or None for all), (stage, libraries or None))
COMPARISONS = {
    "a_E8.5_E9.5": (("E8.5", None), ("E9.5", None)),
    "b_E9.5_E10.5_lib1": (("E9.5", None), ("E10.5", ("1",))),
    "b_E9.5_E10.5_lib2": (("E9.5", None), ("E10.5", ("2",))),
    "b_E9.5_E10.5_lib3": (("E9.5", None), ("E10.5", (LIB_HEART,))),
    "b_E9.5_E10.5_pooled": (("E9.5", None), ("E10.5", None)),
    "c_E10.5_lib12_lib3": (("E10.5", LIB_EMBRYO), ("E10.5", (LIB_HEART,))),
    "c_E10.5_lib1_lib2": (("E10.5", ("1",)), ("E10.5", ("2",))),
}
AK_COMPARISONS = ("b_E9.5_E10.5_lib1", "b_E9.5_E10.5_lib2", "b_E9.5_E10.5_lib3")  # A and K read these together
SCALE_SETS = {"E9.5": (("E9.5", None), BRACKET[0]),  # (challenge cells, Feng stage)
              "E10.5_lib3_on_feng_E14.5": (("E10.5", (LIB_HEART,)), BRACKET[1])}
COMBINED_FIELDS = ("share", "share_corrected")
SHARE_FIELDS = ("share", "pearson", "slope_intercept", "share_top_half", "share_corrected", "n_eff",
                "share_without_top_d2_half", "pearson_d_nonzero", "genes_d_nonzero")
COURSE_KEYS = ("f_max", "f_even", "f_decay")
NOISY_D = "_noisy_d"  # suffix of the corrected courses times the weighted r_d: reported only, read by no statement

KIND_ALONG = "the change of kind lies along Feng's direction"
KIND_NOT = "it does not"
NOT_SEPARATED = "not separated"
PACE_MORE = "the day holds more than a fifth"
PACE_FIFTH_OR_LESS = "about a fifth or less"
ABOVE_BOUND = "above the bound"
INSIDE = "inside"
ABOVE_BOTH = "above both courses"
BELOW_BOTH = "below both courses"
NOT_COMPUTABLE = "not computable"
LABEL_T47 = ("try T47, a labelled diagnostic, counted, never an adoption result: how far the released days move "
             "along Feng's E9.5 -> E14.5 direction, by kind of sample, with M1's own formula; it compares a predictor "
             "with the cuts' target stage, scores nothing and calls no scorer; the read puts nothing forward by itself "
             "and its result can only go to the team captain")
CORRECTION_CHANGES = "the correction changes the statement"
INTERVAL_STRADDLES = "the interval does not separate"


def _ledger_notes(doc: str | None) -> tuple[str, ...]:
    """The ledger sentences of the module docstring (its 'Ledger sentences' list), one string each."""
    if not doc or "Ledger sentences" not in doc:
        raise RuntimeError("the module docstring holds the ledger sentences; run without -OO")
    lines = doc.split("Ledger sentences", 1)[1].split("\n")[1:]
    notes: list[str] = []
    for line in lines:
        if line.startswith("- "):
            notes.append(line[2:].strip())
        elif line.startswith("  ") and notes:
            notes[-1] += " " + line.strip()
        else:
            break
    return tuple(notes)


LEDGER_NOTES = _ledger_notes(__doc__)


# ---------------------------------------------------------------------------------------------------------------------
# Cell tables and their (resampled) block sums


def stage_table(adata: ad.AnnData, genes: Sequence[str], keep_labels, label_key: str = "celltype") -> dict:
    """One T1 stage as blocks of (label, library): the CP10k matrix (expm1 of .X, as t1_stats) of the cells whose label
    is in keep_labels, on `genes`' order, and each cell's block. Library = the obs-name suffix
    (t1_stats.library_tags)."""
    pos = adata.var_names.get_indexer(list(genes))
    if (pos < 0).any():
        raise ValueError(f"a T1 stage lacks {int((pos < 0).sum())} of the reference stage's genes")
    labels = adata.obs[label_key].astype(str).to_numpy()
    libs, _ = t1_stats.library_tags(adata.obs_names)
    keep = np.isin(labels, list(keep_labels))
    X = adata.X[np.flatnonzero(keep)]
    X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    X = X[:, pos].astype(np.float32).tocsr()
    X.data = np.expm1(X.data)
    X.eliminate_zeros()
    keys = [f"{lab}|{li}" for lab, li in zip(labels[keep], np.asarray(libs)[keep], strict=True)]
    uk, inv = np.unique(np.asarray(keys, dtype=str), return_inverse=True)
    order = np.argsort(inv, kind="stable")
    bounds = np.searchsorted(inv[order], np.arange(len(uk) + 1))
    members = [order[bounds[i]:bounds[i + 1]] for i in range(len(uk))]
    return {"cp": X, "inv": inv, "keys": [tuple(k.split("|")) for k in uk], "members": members,
            "n_block": np.bincount(inv, minlength=len(uk)).astype(np.float64)}


def resample_weights(tab: dict, rng: np.random.Generator) -> np.ndarray:
    """Per cell, how many times a bootstrap draw takes it: within each (label, library) block, as many draws with
    replacement as the block has cells (blocks in sorted key order, one call per block)."""
    n = len(tab["inv"])
    picks = [m[rng.integers(0, len(m), len(m))] for m in tab["members"]]
    return np.bincount(np.concatenate(picks) if picks else np.zeros(0, int), minlength=n).astype(np.float64)


def block_sums(tab: dict, weights: np.ndarray | None = None) -> np.ndarray:
    """Per block, the (weighted) sum of the cells' CP10k on every gene."""
    n = len(tab["inv"])
    w = np.ones(n) if weights is None else np.asarray(weights, dtype=np.float64)
    Ind = sparse.csr_matrix((w, (tab["inv"], np.arange(n))), shape=(len(tab["keys"]), n))
    return np.asarray((Ind @ tab["cp"]).todense(), dtype=np.float64)


def cell_set(tab: dict, sums: np.ndarray, labels, libs, min_cells: float) -> dict | None:
    """Mean CP10k over the blocks whose label is in `labels` and library in `libs` (None: every library); None under
    min_cells cells. A within-block resample keeps every block's size, so the set's cell number is fixed."""
    ls = set(labels)
    m = np.array([lab in ls and (libs is None or li in libs) for lab, li in tab["keys"]], dtype=bool)
    n = float(tab["n_block"][m].sum()) if m.any() else 0.0
    if n == 0 or n < min_cells:
        return None
    return {"n": n, "cp": sums[m].sum(0) / n}


# ---------------------------------------------------------------------------------------------------------------------
# Shares, scales and the courses


def share_stats(change: np.ndarray, d: np.ndarray, mask: np.ndarray, scale: float | None) -> dict:
    """The multiple of d that fits `change` over the genes in mask (least squares through the origin), and beside it
    the Pearson correlation, the slope with an intercept, the share over the half of the genes with the larger |d|
    (ceil(n / 2) genes, ties in gene order), the share over the N48 scale, the effective gene number (sum d^2)^2 /
    sum d^4, the share without the top-|d| genes that together carry half of sum d^2 (the fewest such genes, ties in
    gene order), and the Pearson correlation and gene number over the genes with d != 0."""
    x, y = np.asarray(d, np.float64)[mask], np.asarray(change, np.float64)[mask]
    n = len(x)
    nan = float("nan")
    out = dict.fromkeys(SHARE_FIELDS, nan)
    if n == 0 or float(x @ x) == 0.0:
        return out
    ss = float(x @ x)
    out["share"] = float(y @ x / ss)
    if n >= 3 and np.std(x) > 0 and np.std(y) > 0:
        out["pearson"] = float(np.corrcoef(x, y)[0, 1])
    if n >= 2 and np.var(x) > 0:
        out["slope_intercept"] = float(np.mean((x - x.mean()) * (y - y.mean())) / np.var(x))
    order = np.argsort(-np.abs(x), kind="stable")
    top = order[: n - n // 2]
    if float(x[top] @ x[top]) > 0:
        out["share_top_half"] = float(y[top] @ x[top] / (x[top] @ x[top]))
    if scale is not None and np.isfinite(scale) and scale > 0:
        out["share_corrected"] = out["share"] / float(scale)
    out["n_eff"] = ss ** 2 / float(np.sum(x ** 4))
    cum = np.cumsum(x[order] ** 2)
    k = int(np.searchsorted(cum, 0.5 * ss)) + 1  # the fewest top genes whose d^2 reaches half the sum
    rest = order[k:]
    if rest.size and float(x[rest] @ x[rest]) > 0:
        out["share_without_top_d2_half"] = float(y[rest] @ x[rest] / (x[rest] @ x[rest]))
    nz = x != 0
    out["genes_d_nonzero"] = float(nz.sum())
    if nz.sum() >= 3 and np.std(x[nz]) > 0 and np.std(y[nz]) > 0:
        out["pearson_d_nonzero"] = float(np.corrcoef(x[nz], y[nz])[0, 1])
    return out


def origin_slope(y: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    """Least-squares slope through the origin of y on x, and their Pearson correlation (NaN where undefined)."""
    nan = float("nan")
    if len(x) == 0 or float(x @ x) == 0.0:
        return nan, nan
    r = float(np.corrcoef(x, y)[0, 1]) if len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0 else nan
    return float(y @ x / (x @ x)), r


def contrast_reliability(c: np.ndarray, c_a: np.ndarray, c_b: np.ndarray, mask: np.ndarray) -> float:
    """N48's scalar James-Stein formula on a contrast: over the genes in mask where c and both halves' c_a and c_b are
    finite, max(0, mean(c^2) - mean(((c_a - c_b) / 2)^2)) / mean(c^2); NaN where no gene or mean(c^2) = 0."""
    c, c_a, c_b = (np.asarray(v, np.float64) for v in (c, c_a, c_b))
    ok = np.asarray(mask, bool) & np.isfinite(c) & np.isfinite(c_a) & np.isfinite(c_b)
    if not ok.any():
        return float("nan")
    msq = float(np.mean(c[ok] ** 2))
    if msq == 0:
        return float("nan")
    return max(msq - float(np.mean(((c_a[ok] - c_b[ok]) / 2.0) ** 2)), 0.0) / msq


def pair_median(slopes: Sequence[float]) -> float:
    """The scale: the median over the pairs' finite slopes (NaN where none)."""
    fin = [s for s in slopes if np.isfinite(s)]
    return float(np.median(fin)) if fin else float("nan")


def combine(values: Mapping[str, float], weights: Mapping[str, float]) -> tuple[float, float]:
    """(weighted mean, unweighted mean) over the groups with a finite value; the weighted one over those of them with a
    positive weight, renormalised; NaN where none."""
    ok = [g for g, v in values.items() if np.isfinite(v)]
    unw = float(np.mean([values[g] for g in ok])) if ok else float("nan")
    wok = [g for g in ok if weights.get(g, 0.0) > 0]
    tot = float(sum(weights[g] for g in wok))
    wtd = float(sum(weights[g] * values[g] for g in wok) / tot) if tot > 0 else float("nan")
    return wtd, unw


def combined_groups(values: Mapping[str, float], weights: Mapping[str, float]) -> tuple[list[str], list[str]]:
    """The groups `combine` uses: (weighted, unweighted)."""
    ok = sorted(g for g, v in values.items() if np.isfinite(v))
    return [g for g in ok if weights.get(g, 0.0) > 0], ok


def decay_x(ratio: float) -> float:
    """x > 0 with (1 - x) / (1 - x^5) = ratio, i.e. 1 + x + x^2 + x^3 + x^4 = 1 / ratio (x = 1 at 0.2); x = 0 where
    ratio >= 1 (the day holds the whole change or more: nothing is left for a decaying course); NaN where ratio <= 0
    (no decaying course) or not finite."""
    if not np.isfinite(ratio) or ratio <= 0:
        return float("nan")
    if ratio >= 1:
        return 0.0
    roots = np.roots([1.0, 1.0, 1.0, 1.0, 1.0 - 1.0 / ratio])
    real = [float(r.real) for r in roots if abs(r.imag) < 1e-9 and r.real > 0]
    return min(real, key=lambda r: abs(1 + r + r ** 2 + r ** 3 + r ** 4 - 1 / ratio))


def courses(A: float, L: float) -> dict:
    """The test's applied fraction under the two simple courses and the bound, from the day's share A and the scale L:
    f_max = L - A, f_even = (L - A) / 2 and f_decay = L (x - x^3) / (1 - x^5), written L x (1 + x) / (1 + x + x^2 +
    x^3 + x^4) (the same, without the 0/0 at x = 1); f_decay = 0 where A / L >= 1, NaN where A / L <= 0."""
    nan = float("nan")
    ratio = A / L if np.isfinite(A) and np.isfinite(L) and L != 0 else nan
    x = decay_x(ratio)
    f_decay = L * x * (1 + x) / (1 + x + x ** 2 + x ** 3 + x ** 4) if np.isfinite(x) else nan
    return {"A_over_L": ratio, "x": x, "f_max": L - A, "f_even": (L - A) / 2.0, "f_decay": f_decay}


# ---------------------------------------------------------------------------------------------------------------------
# One measurement (the point or one bootstrap draw)


def common_support(out: Mapping[str, float], groups: Sequence[str], names: Sequence[str], field: str) -> list[str]:
    """The groups with a finite `field` in every comparison of `names`."""
    return [g for g in groups if all(np.isfinite(out[f"{n}/{g}/{field}"]) for n in names)]


def reading_values(out: Mapping[str, float], setup: Mapping, info: dict | None = None) -> dict:
    """A, K, P, L and the courses, raw and corrected, from one measurement's per-group values. A, K, A_corr and K_corr
    over one set of groups: those with a share and a share_corrected in libraries 1, 2 and 3 (the intersection of the
    two common supports). P over the groups with an (a) value and, beside it, over those whose cells are their own at
    E8.5 and E9.5. The corrected courses are in units of Feng's true direction; times the weighted r_d (NOISY_D) they
    are reported only. info, when given, receives the group lists."""
    groups, w = list(setup["dirs"]), setup["weights"]

    def wmean(name, field, support):
        return combine({g: out[f"{name}/{g}/{field}"] for g in support}, w)[0]

    rv: dict[str, float] = {}
    sup_raw = common_support(out, groups, AK_COMPARISONS, "share")
    sup_corr = common_support(out, groups, AK_COMPARISONS, "share_corrected")
    sup = [g for g in sup_raw if g in sup_corr]
    for field, tag in (("share", ""), ("share_corrected", "_corr")):
        A = (wmean(AK_COMPARISONS[0], field, sup) + wmean(AK_COMPARISONS[1], field, sup)) / 2.0
        rv[f"A{tag}"] = A
        rv[f"K{tag}"] = wmean(AK_COMPARISONS[2], field, sup) - A
    p_sup = common_support(out, groups, ("a_E8.5_E9.5",), "share")
    p_own = [g for g in p_sup if all(setup["source"].get(s, {}).get(g) == "own" for s in ("E8.5", "E9.5"))]
    rv["P"] = wmean("a_E8.5_E9.5", "share", p_sup)
    rv["P_own_cells"] = wmean("a_E8.5_E9.5", "share", p_own)
    sc = setup["scales"]
    rv["r_d_weighted"] = combine({g: (sc.get(g) if sc.get(g) is not None else float("nan")) for g in sup}, w)[0]
    L = out["scale_E9.5/median_slope"]
    L_corr = out["scale_E9.5/median_slope_corrected"]
    rv["L"], rv["L_all_pairs"], rv["L_corr"] = L, out["scale_E9.5/median_slope_all_pairs"], L_corr
    for k, v in courses(rv["A"], L).items():
        rv[k] = v
    cc = courses(rv["A_corr"], L_corr)
    rv["A_over_L_corr"], rv["x_corr"] = cc["A_over_L"], cc["x"]
    for k in COURSE_KEYS:
        rv[f"{k}_corr"] = cc[k]  # in units of Feng's true direction: what the corrected statement reads
        rv[f"{k}_corr{NOISY_D}"] = cc[k] * rv["r_d_weighted"]  # the best multiple of the noisy d: reported only
    if info is not None:
        info.update({"A_K_support": sup, "A_K_excluded": [g for g in groups if g not in sup],
                     "A_K_excluded_no_share": [g for g in groups if g not in sup_raw],
                     "A_K_excluded_no_share_corrected": [g for g in groups if g in sup_raw and g not in sup_corr],
                     "P_support": p_sup, "P_own_cells_support": p_own,
                     "P_fallback_groups": [g for g in p_sup if g not in p_own]})
    return rv


def measure(tabs: Mapping[str, dict], sums: Mapping[str, np.ndarray], setup: Mapping,
            info: dict | None = None) -> dict:
    """Every scalar of the read from one set of block sums: a flat mapping name -> value."""
    eps, mc = setup["eps"], setup["min_cells"]
    dirs, masks, labels, scales = setup["dirs"], setup["masks"], setup["labels"], setup["scales"]
    out: dict[str, float] = {}

    def mean_of(spec, g):
        stage, libs = spec
        labs = labels.get(stage, {}).get(g)
        if not labs or stage not in tabs:
            return None
        return cell_set(tabs[stage], sums[stage], labs, libs, mc)

    for name, (sa, sb) in COMPARISONS.items():
        per = {}
        for g in dirs:
            a, b = mean_of(sa, g), mean_of(sb, g)
            if a is None or b is None:
                st = dict.fromkeys(SHARE_FIELDS, float("nan"))
            else:
                st = share_stats(dl.lfc_ratio(a, b, eps), dirs[g], masks[g], scales.get(g))
            for k, v in st.items():
                out[f"{name}/{g}/{k}"] = v
            per[g] = st
        for f in COMBINED_FIELDS:
            wtd, unw = combine({g: per[g][f] for g in per}, setup["weights"])
            out[f"{name}/combined_weighted/{f}"] = wtd
            out[f"{name}/combined_unweighted/{f}"] = unw

    for sname, (spec, _) in SCALE_SETS.items():
        kept, every, corrected = [], [], []
        for g, h in setup["pairs"]:
            a, b = mean_of(spec, h), mean_of(spec, g)
            fx = setup["feng_contrast"][sname].get((g, h))
            m = masks[g] & masks[h]
            if a is None or b is None or fx is None:
                s, r = float("nan"), float("nan")
            else:
                s, r = origin_slope(dl.lfc_ratio(a, b, eps)[m], fx[m])
            rc = setup["feng_reliability"][sname].get((g, h), float("nan"))
            sc = s / rc if np.isfinite(rc) and rc > 0 else float("nan")
            out[f"scale_{sname}/{g}~{h}/slope"] = s
            out[f"scale_{sname}/{g}~{h}/pearson"] = r
            out[f"scale_{sname}/{g}~{h}/slope_corrected"] = sc
            every.append(s)
            if (g, h) not in setup["overlap_pairs"]:
                kept.append(s)
                corrected.append(sc)
        out[f"scale_{sname}/median_slope"] = pair_median(kept)
        out[f"scale_{sname}/median_slope_all_pairs"] = pair_median(every)
        out[f"scale_{sname}/median_slope_corrected"] = pair_median(corrected)
        if info is not None:  # the pairs each median rests on (a finite value)
            names = [f"{g}~{h}" for g, h in setup["pairs"]]
            kept_names = [p for p, (g, h) in zip(names, setup["pairs"], strict=True)
                          if (g, h) not in setup["overlap_pairs"]]
            info.setdefault("pairs", {})[sname] = {
                "L": [p for p in kept_names if np.isfinite(out[f"scale_{sname}/{p}/slope"])],
                "L_corr": [p for p in kept_names if np.isfinite(out[f"scale_{sname}/{p}/slope_corrected"])],
                "L_all_pairs": [p for p in names if np.isfinite(out[f"scale_{sname}/{p}/slope"])]}

    for k, v in reading_values(out, setup, info).items():
        out[f"reading/{k}"] = v
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The statements


def _interval(rec: Mapping) -> tuple[float, float]:
    return rec["p5"], rec["p95"]


def kind_statement(lo: float, hi: float) -> str:
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return NOT_COMPUTABLE
    if lo > KIND_AT:
        return KIND_ALONG
    if hi < KIND_AT:
        return KIND_NOT
    return NOT_SEPARATED


def pace_statement(lo: float, hi: float) -> str:
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return NOT_COMPUTABLE
    if lo > PACE_FIFTH:
        return PACE_MORE
    if hi < PACE_UPPER and lo <= PACE_FIFTH:  # below 0.3 and including 0.2 or lying below it
        return PACE_FIFTH_OR_LESS
    return NOT_SEPARATED


def fraction_statement(f_max: Mapping, f_even: Mapping, f_decay: Mapping, fraction: float = ENTRY_FRACTION) -> str:
    """The statement on the entry's fraction: "above the bound" if it is above f_max's upper percentile (read from f_max
    alone); else against the range from the lower of the two courses' lower percentiles to the higher of their upper
    percentiles."""
    if not np.isfinite(f_max["p95"]):
        return NOT_COMPUTABLE
    if fraction > f_max["p95"]:
        return ABOVE_BOUND
    vals = [f_even["p5"], f_even["p95"], f_decay["p5"], f_decay["p95"]]
    if not all(np.isfinite(v) for v in vals):
        return NOT_COMPUTABLE
    lo, hi = min(f_even["p5"], f_decay["p5"]), max(f_even["p95"], f_decay["p95"])
    if lo <= fraction <= hi:
        return INSIDE
    return ABOVE_BOTH if fraction > hi else BELOW_BOTH


def undefined(recs: Sequence[Mapping], n_boot: int) -> str | None:
    """'not computable (n of N draws undefined)' when a quantity a statement reads is undefined in any draw or at the
    point; None when every one is defined."""
    n = max(int(r["nonfinite_draws"]) for r in recs)
    if n > 0 or not all(np.isfinite(r["value"]) for r in recs):
        return f"{NOT_COMPUTABLE} ({n} of {n_boot} draws undefined)"
    return None


def joined(raw: str, corrected: str) -> str:
    """The statement on the raw and the corrected pair: the shared one, a 'not computable' one as it stands, else
    'not separated'."""
    for s in (raw, corrected):
        if s.startswith(NOT_COMPUTABLE):
            return s
    return raw if raw == corrected else NOT_SEPARATED


def not_separated_cause(raw: str, corrected: str) -> str | None:
    """Why a joined statement reads 'not separated': the correction changes a separated statement, or an interval
    straddles a threshold (raw, corrected or both); None where the joined statement is not 'not separated'."""
    if joined(raw, corrected) != NOT_SEPARATED:
        return None
    if raw != NOT_SEPARATED and corrected != NOT_SEPARATED:
        return CORRECTION_CHANGES
    which = [lab for lab, s in (("raw", raw), ("corrected", corrected)) if s == NOT_SEPARATED]
    return f"{INTERVAL_STRADDLES} ({' and '.join(which)})"


def entry_fraction_statement(recs: Sequence[Mapping], n_boot: int) -> str:
    """The fraction statement with undefined draws: 'above the bound' needs only f_max defined in every draw; the
    comparison with the two courses needs f_even and f_decay defined in every draw (f_decay has no value where A / L
    <= 0), else it reads 'not computable (n of N draws undefined)'."""
    f_max, f_even, f_decay = recs
    u = undefined([f_max], n_boot)
    if u is not None:
        return u
    if fraction_statement(f_max, f_even, f_decay) == ABOVE_BOUND:
        return ABOVE_BOUND
    return undefined([f_even, f_decay], n_boot) or fraction_statement(f_max, f_even, f_decay)


def statements(r: Mapping[str, Mapping], n_boot: int, pairs: Mapping[str, Sequence[str]] | None = None) -> dict:
    """The three statements from the reading's records (value, p5, p95, nonfinite_draws). pairs: the pairs L and
    L_corr rest on; where L_corr rests on fewer, the corrected pace and fraction statements are 'not computable'."""
    out = {"kind": undefined([r["K"]], n_boot) or kind_statement(*_interval(r["K"]))}
    short = None
    if pairs is not None and len(pairs["L_corr"]) < len(pairs["L"]):
        short = f"{NOT_COMPUTABLE} (L_corr rests on {len(pairs['L_corr'])} of the {len(pairs['L'])} pairs of L)"
    for tag in ("", "_corr"):
        lab = "raw" if tag == "" else "corrected"
        q = r[f"A_over_L{tag}"]
        recs = [r[f"{k}{tag}"] for k in COURSE_KEYS]
        if tag and short:
            out[f"pace_{lab}"] = out[f"entry_fraction_{lab}"] = short
            continue
        out[f"pace_{lab}"] = undefined([q], n_boot) or pace_statement(*_interval(q))
        out[f"entry_fraction_{lab}"] = entry_fraction_statement(recs, n_boot)
    for s in ("pace", "entry_fraction"):
        out[s] = joined(out[f"{s}_raw"], out[f"{s}_corrected"])
        out[f"{s}_why_not_separated"] = not_separated_cause(out[f"{s}_raw"], out[f"{s}_corrected"])
    return out


def summarise(point: Mapping[str, float], draws: Mapping[str, np.ndarray], n_boot: int) -> dict:
    """Per quantity: the point value, the percentiles over its finite draws and the number of undefined draws."""
    values = {}
    for k, v in point.items():
        d = np.asarray(draws[k], np.float64)
        fin = d[np.isfinite(d)]
        lo, hi = (np.percentile(fin, PERCENTILES) if fin.size else (np.nan, np.nan))
        values[k] = {"value": float(v), "p5": float(lo), "p95": float(hi), "nonfinite_draws": int(n_boot - fin.size)}
    return values


# ---------------------------------------------------------------------------------------------------------------------
# The calibration arm: M1's own 0.2 step read with the same share code


def calibration(adata: ad.AnnData, genes: Sequence[str], dirs: Mapping[str, np.ndarray], masks: Mapping,
                labels: Mapping[str, Sequence[str]], held: Sequence[str], scales: Mapping, weights: Mapping,
                eps: float, min_cells: float, label_key: str = "celltype") -> dict:
    """t1_scope.external_shift with fraction CAL_FRACTION (stream CAL_STREAM, its defaults otherwise, as M1's build
    calls it) on library 3's cells of the direction groups at E10.5, on the file's own gene set; then, per group, the
    share of before -> after (mean CP10k on `genes`), and the combined shares."""
    libs, _ = t1_stats.library_tags(adata.obs_names)
    lab = adata.obs[label_key].astype(str).to_numpy()
    group_of = {lb: g for g in sorted(dirs) for lb in labels.get(g) or []}
    cells = np.flatnonzero((np.asarray(libs) == LIB_HEART) & np.isin(lab, list(group_of)))
    pos = adata.var_names.get_indexer(list(genes))
    full = {}
    for g, v in dirs.items():
        w = np.full(adata.n_vars, np.nan)
        w[pos] = v
        full[g] = w
    nan_rec = dict.fromkeys(SHARE_FIELDS, float("nan"))
    per = {g: dict(nan_rec) for g in dirs}
    if cells.size:
        sub = adata[cells].copy()
        moved, _ = t1_scope.external_shift(sub, full, CAL_FRACTION, group_of, held, label_key=label_key,
                                           seed=list(CAL_STREAM))
        sl = lab[cells]

        def cp(a, m):
            X = a.X[np.flatnonzero(m)]
            X = X.toarray() if sparse.issparse(X) else np.asarray(X)
            return {"cp": np.expm1(np.asarray(X, np.float64)[:, pos]).mean(0)}

        for g in dirs:
            m = np.isin(sl, list(labels.get(g) or []))
            if m.sum() == 0 or m.sum() < min_cells:
                continue
            per[g] = share_stats(dl.lfc_ratio(cp(sub, m), cp(moved, m), eps), dirs[g], masks[g], scales.get(g))
    shares = {g: per[g]["share"] for g in per}
    wtd, unw = combine(shares, weights)
    return {"fraction": CAL_FRACTION, "stream": list(CAL_STREAM), "groups": per, "combined_weighted": wtd,
            "combined_unweighted": unw, "combined_groups": combined_groups(shares, weights)[0]}


# ---------------------------------------------------------------------------------------------------------------------
# The read


def group_labels_by_stage(stages: Mapping[str, ad.AnnData], groups: Sequence[str], cfg: Mapping,
                          label_key: str) -> tuple[dict, dict]:
    """stage -> group -> labels (own, else the config fallback's), and stage -> group -> 'own', 'fallback:<group>' or
    'none'."""
    labels, source = {}, {}
    for s, a in stages.items():
        mapping = feng_dirs.label_groups(a.obs[label_key].astype(str).unique(), cfg)
        labels[s], source[s] = {}, {}
        for g in groups:
            own = dl.group_labels(mapping, g)
            labs = dl.antecedent_labels(mapping, g, (), cfg.get("fallback"))
            labels[s][g] = labs
            source[s][g] = "own" if own else (f"fallback:{cfg['fallback'][g]}" if labs else "none")
    return labels, source


def feng_pair_contrasts(feng: dict, pairs, cfg: Mapping, idx: np.ndarray, masks: Mapping, eps: float):
    """Per scale set and pair (g, h): Feng's contrast (h's cells against g's, both halves pooled) on the T1 genes, and
    its reliability r_c from the two halves' contrasts (contrast_reliability over the genes both groups keep)."""
    contrast, reliability = {}, {}
    for sname, (_, fstage) in SCALE_SETS.items():
        contrast[sname], reliability[sname] = {}, {}
        for g, h in pairs:
            lg, lh = cfg["external_labels"][g][0], cfg["external_labels"][h][0]
            fg, fh = dl.feng_sel(feng, fstage, lg), dl.feng_sel(feng, fstage, lh)
            if fg is None or fh is None:
                contrast[sname][(g, h)] = None
                reliability[sname][(g, h)] = float("nan")
                continue
            c = dl.on_t1(dl.lfc_ratio(fh, fg, eps), idx)
            contrast[sname][(g, h)] = c
            halves = []
            for half in ("A", "B"):
                hg, hh = dl.feng_sel(feng, fstage, lg, halves=(half,)), dl.feng_sel(feng, fstage, lh, halves=(half,))
                halves.append(None if hg is None or hh is None else dl.on_t1(dl.lfc_ratio(hh, hg, eps), idx))
            reliability[sname][(g, h)] = (float("nan") if any(x is None for x in halves)
                                          else contrast_reliability(c, halves[0], halves[1], masks[g] & masks[h]))
    return contrast, reliability


def direction_reliability(feng: dict, genes: Sequence[str], cfg: Mapping, dirs: Mapping[str, np.ndarray],
                          masks: Mapping[str, np.ndarray]) -> dict[str, float]:
    """Per group, r_d: contrast_reliability of Feng's direction d against its two halves' directions d_a and d_b
    (feng_dirs.direction_halves on the bracket, the same statistics) over the group's share mask."""
    halves = feng_dirs.direction_halves(feng, list(genes), dict(cfg), BRACKET, groups=list(dirs))
    return {g: (contrast_reliability(dirs[g], halves[g]["d_a"], halves[g]["d_b"], masks[g]) if g in halves
                else float("nan")) for g in dirs}


def reading(stages: Mapping[str, ad.AnnData], feng: dict, held: Sequence[str], n48: Mapping[str, float] | None,
            cfg: Mapping, n_boot: int = N_BOOT, seed=BOOT_SEED, calibrate: bool = False) -> dict:
    """The whole read on loaded inputs: T1 stages (stage -> AnnData), Feng's statistics (only the bracket's two ends),
    M1's held list and N48's scales (reported beside r_d only; r_d is computed here from the bracket's halves).
    calibrate: run the calibration arm (it needs counts recoverable from E10.5's .X, as the released files have)."""
    label_key = cfg.get("label_key", "celltype")
    genes = list(stages[GENE_STAGE].var_names)
    dirs, _, idx = feng_dirs.directions_from_stats(feng, genes, cfg, BRACKET)
    groups = sorted(dirs)
    dirs = {g: dirs[g] for g in groups}
    held_set = set(map(str, held))
    not_held = np.array([g not in held_set for g in genes], dtype=bool)
    masks = {g: np.isfinite(dirs[g]) & not_held for g in groups}
    r_d = direction_reliability(feng, genes, cfg, dirs, masks)
    n48 = n48 or {}
    labels, source = group_labels_by_stage(stages, groups, cfg, label_key)
    keep = {s: {lab for g in groups for lab in labels[s][g]} for s in stages}
    tabs = {s: stage_table(stages[s], genes, keep[s], label_key) for s in stages}
    eps = cfg["metrics"]["pseudocount_cp10k"]
    mc = cfg["metrics"]["min_cells_t1"]
    pairs = list(combinations(groups, 2))
    feng_contrast, feng_rel = feng_pair_contrasts(feng, pairs, cfg, idx, masks, eps)
    overlap = [(g, h) for g, h in pairs
               if set(cfg["external_labels"][g][0]) & set(cfg["external_labels"][h][0])]
    heart = tabs.get("E10.5")
    weights = {}
    for g in groups:
        labs = set(labels.get("E10.5", {}).get(g) or [])
        weights[g] = (float(sum(n for (lab, li), n in zip(heart["keys"], heart["n_block"], strict=True)
                                if lab in labs and li == LIB_HEART)) if heart is not None else 0.0)
    setup = {"eps": eps, "min_cells": mc, "dirs": dirs, "masks": masks, "labels": labels, "source": source,
             "scales": r_d, "weights": weights, "pairs": pairs,
             "feng_contrast": feng_contrast, "feng_reliability": feng_rel, "overlap_pairs": set(overlap)}
    info: dict = {}
    point = measure(tabs, {s: block_sums(t) for s, t in tabs.items()}, setup, info)
    rng = np.random.default_rng(list(seed))
    draws = {k: np.empty(n_boot) for k in point}
    t0 = time.time()
    for i in range(n_boot):
        sums = {s: block_sums(t, resample_weights(t, rng)) for s, t in tabs.items()}
        m = measure(tabs, sums, setup)
        for k, v in m.items():
            draws[k][i] = v
    seconds = time.time() - t0
    values = summarise(point, draws, n_boot)
    r = {k.split("/", 1)[1]: values[k] for k in values if k.startswith("reading/")}
    groups_of = {}
    for name in COMPARISONS:
        for f in COMBINED_FIELDS:
            wg, ug = combined_groups({g: point[f"{name}/{g}/{f}"] for g in groups}, weights)
            groups_of[f"{name}/combined_weighted/{f}"] = wg
            groups_of[f"{name}/combined_unweighted/{f}"] = ug
    pair_genes = {f"{g}~{h}": int((masks[g] & masks[h]).sum()) for g, h in pairs}
    scale_pairs = info.pop("pairs")
    res = {"groups": groups, "genes_per_group": {g: int(masks[g].sum()) for g in groups},
           "genes_per_pair": pair_genes, "label_source": source, "r_d": r_d,
           "n48_scale": {g: n48.get(g) for g in groups}, "values": values, "reading": r,
           "statements": statements(r, n_boot, scale_pairs["E9.5"]), "support": info, "combined_groups": groups_of,
           "scale_pairs": scale_pairs,
           "excluded_pairs_shared_feng_labels": [f"{g}~{h}" for g, h in overlap],
           "feng_reliability": {sname: {f"{g}~{h}": v for (g, h), v in d.items()} for sname, d in feng_rel.items()},
           "n_boot": n_boot, "seed": list(seed), "percentiles": list(PERCENTILES),
           "bootstrap_seconds": round(seconds, 1)}
    if calibrate:
        cal = calibration(stages["E10.5"], genes, dirs, masks, labels["E10.5"], held, setup["scales"], weights,
                          eps, mc, label_key)
        s_cal = cal["combined_weighted"]
        ok = np.isfinite(s_cal) and s_cal != 0
        cal["A_in_M1_fraction"] = r["A"]["value"] * CAL_FRACTION / s_cal if ok else float("nan")
        cal["K_in_M1_fraction"] = r["K"]["value"] * CAL_FRACTION / s_cal if ok else float("nan")
        res["calibration"] = cal
    return res


# ---------------------------------------------------------------------------------------------------------------------
# Inputs and the run


def feng_inputs(path: Path | None = None) -> tuple[dict, dict]:
    """Feng's clean statistics (manifest-checked) cut to the bracket's two ends (feng_dirs.bracket_stats): no other
    Feng stage reaches the read."""
    path = feng_dirs.clean_stats_path() if path is None else Path(path)
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    return feng_dirs.bracket_stats(st, BRACKET), {"file": gv._rel(path), "sha256": entry["sha256"]}


def n48_scales(path: Path = N48_JSON, field: str = "js_without_held") -> dict[str, float]:
    """Each group's scalar James-Stein scale from N48's record, reported beside the read's own r_d and entering
    nothing: by default 'js_without_held' (N48 leaves out only the pre-registered held list, not R2's union, so its
    gene set is not the shares'); field='js' gives N48's rule's own record (held genes included)."""
    rep = json.loads(Path(path).read_text())
    return {g: r[field]["scale"] for g, r in rep["groups"].items() if r[field]["scale"] is not None}


def _fmt(rec: Mapping, n_boot: int | None = None) -> str:
    s = f"{rec['value']:+.3f} [{rec['p5']:+.3f}, {rec['p95']:+.3f}]"
    if n_boot is not None:
        s += f" ({rec['nonfinite_draws']} of {n_boot} draws undefined)"
    return s


def text(res: Mapping) -> str:
    r = res["read"]
    rd, st, nb = r["reading"], r["statements"], r["n_boot"]
    sup = r["support"]

    def f(k):
        return _fmt(rd[k], nb)

    def why(s):
        cause = st.get(f"{s}_why_not_separated")
        return f" [{cause}]" if cause else ""

    sp = r["scale_pairs"]["E9.5"]
    lines = [res["label"],
             f"Kind: {st['kind']} (K = {f('K')}, against {KIND_AT})",
             f"Pace: {st['pace']}{why('pace')} (raw: {st['pace_raw']}, A / L = {f('A_over_L')}; corrected: "
             f"{st['pace_corrected']}, A_corr / L_corr = {f('A_over_L_corr')}; against {PACE_FIFTH} and {PACE_UPPER})",
             f"Entry's fraction {ENTRY_FRACTION}: {st['entry_fraction']}{why('entry_fraction')} (two courses among "
             f"many, not a bracket)",
             f"    raw: {st['entry_fraction_raw']}; f_max {f('f_max')}; f_even {f('f_even')}; f_decay {f('f_decay')}",
             f"    corrected (units of Feng's true direction): {st['entry_fraction_corrected']}; "
             f"f_max {f('f_max_corr')}; f_even {f('f_even_corr')}; f_decay {f('f_decay_corr')}",
             f"    reported only, the corrected courses times the weighted r_d {f('r_d_weighted')} (units of the noisy "
             f"d): f_max {f('f_max_corr' + NOISY_D)}; f_even {f('f_even_corr' + NOISY_D)}; f_decay "
             f"{f('f_decay_corr' + NOISY_D)}",
             f"A {f('A')}; K {f('K')}; A_corr {f('A_corr')}; K_corr {f('K_corr')}",
             f"    A, K, A_corr and K_corr over {sup['A_K_support']} (excluded: {sup['A_K_excluded']}; no share in "
             f"every library: {sup['A_K_excluded_no_share']}; no corrected share in every library: "
             f"{sup['A_K_excluded_no_share_corrected']})",
             f"P over own cells {f('P_own_cells')} over {sup['P_own_cells_support']}; beside it P {f('P')} over "
             f"{sup['P_support']} (fallback groups: {sup['P_fallback_groups']})",
             f"L {f('L')} (pairs with shared Feng labels excluded: {r['excluded_pairs_shared_feng_labels']}); "
             f"L over every pair {f('L_all_pairs')}; L_corr {f('L_corr')}",
             *(f"    {k} rests on {len(sp[k])} pairs: {sp[k]}" for k in ("L", "L_corr", "L_all_pairs")),
             f"(value [p{PERCENTILES[0]:g}, p{PERCENTILES[1]:g}], {nb} draws, stream {r['seed']})",
             "Groups (genes kept; r_d; N48 js_without_held, beside it only; label source E8.5 / E9.5 / E10.5): "
             + "; ".join(f"{g} ({r['genes_per_group'][g]}; {r['r_d'][g]:.3f}; {r['n48_scale'][g]}; "
                         + " / ".join(r["label_source"].get(s, {}).get(g, "-") for s in T1_STAGES) + ")"
                         for g in r["groups"])]
    if res.get("n48_js") is not None:
        lines.append(f"N48 'js' (held genes included), beside it only: {res['n48_js']}")
    cal = r.get("calibration")
    if cal is not None:
        lines.append(f"Calibration: M1's step of {cal['fraction']} on E10.5 library 3 reads a combined weighted share "
                     f"of {cal['combined_weighted']:+.3f} (unweighted {cal['combined_unweighted']:+.3f}; groups "
                     f"{cal['combined_groups']}); A in M1's fraction {cal['A_in_M1_fraction']:+.3f}, K "
                     f"{cal['K_in_M1_fraction']:+.3f}")
        for g, s in cal["groups"].items():
            lines.append(f"    {g}: share {s['share']:+.3f}")
    v, cg = r["values"], r["combined_groups"]
    for name in COMPARISONS:
        kw, ku = f"{name}/combined_weighted/share", f"{name}/combined_unweighted/share"
        lines.append(f"{name}: combined weighted {_fmt(v[kw], nb)} over {cg[kw]}, unweighted {_fmt(v[ku], nb)} over "
                     f"{cg[ku]}")
        for g in r["groups"]:
            def one(field, g=g, name=name):
                return v[f"{name}/{g}/{field}"]["value"]
            lines.append(f"    {g}: share {_fmt(v[f'{name}/{g}/share'], nb)}, r {one('pearson'):+.3f}, "
                         f"slope+int {one('slope_intercept'):+.3f}, top half {one('share_top_half'):+.3f}, "
                         f"corrected {one('share_corrected'):+.3f}, n_eff {one('n_eff'):.1f}, without the top-d^2 "
                         f"half {one('share_without_top_d2_half'):+.3f}, r over d != 0 {one('pearson_d_nonzero'):+.3f} "
                         f"({one('genes_d_nonzero'):.0f} genes)")
    for sname in SCALE_SETS:
        lines.append(f"scale {sname}: median slope {_fmt(v[f'scale_{sname}/median_slope'], nb)}; every pair "
                     f"{_fmt(v[f'scale_{sname}/median_slope_all_pairs'], nb)}; corrected "
                     f"{_fmt(v[f'scale_{sname}/median_slope_corrected'], nb)}")
        for pair, n in r["genes_per_pair"].items():
            lines.append(f"    {pair} ({n} genes): slope {_fmt(v[f'scale_{sname}/{pair}/slope'], nb)}, r "
                         f"{v[f'scale_{sname}/{pair}/pearson']['value']:+.3f}, r_c "
                         f"{r['feng_reliability'][sname][pair]:.3f}, corrected "
                         f"{v[f'scale_{sname}/{pair}/slope_corrected']['value']:+.3f}")
    lines.extend(f"Note: {s}" for s in LEDGER_NOTES)
    lines.append("Limit: a labelled diagnostic; never adoption evidence; puts nothing forward by itself; it can only "
                 "go to the team captain")
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(lines) + "\n"


def run(out_name: str, files: Mapping[str, Path] = T1_FILES) -> dict:
    if ROW not in n12.tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    out = RESULTS / out_name
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean (a counted read runs on a committed tree)")
    cfg = feng_dirs.load_config()
    union = hygiene.load_union()
    preregistered = t1_scope.held_flat_genes()
    held = gc.held_feng(union, preregistered)
    feng, feng_rec = feng_inputs()
    n48 = n48_scales()  # reported beside r_d only
    js = n48_scales(field="js")
    stages = {s: ad.read_h5ad(files[s]) for s in T1_STAGES}
    read = reading(stages, feng, held, n48, cfg, calibrate=True)
    res = {"label": LABEL_T47, "read": read, "feng": feng_rec, "held_genes": len(held), "n48_js": js}
    res["stamp"] = gv.stamp(started, [*files.values(), feng_dirs.clean_stats_path(), t1_scope.HELD_FLAT_GENES,
                                      REPO / union["file"], N48_JSON])
    gv.write_outputs(out, STEM, res, text(res))
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="try T47's read")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    run(a.out)


if __name__ == "__main__":
    main()
