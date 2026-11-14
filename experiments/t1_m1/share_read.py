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
  share    s = sum(change x d) / sum(d x d); beside it the Pearson correlation, the slope with an intercept, the share
           over the half of the genes with the larger |d|, and s / the group's N48 James-Stein scale.
  combined the groups' shares averaged, weighted by the groups' cells in E10.5's library 3, and unweighted.
  (a) E8.5 -> E9.5; (b) E9.5 -> E10.5 library 1, 2, 3 and pooled; (c) E10.5 libraries 1+2 -> 3 and 1 -> 2;
  (d) per pair of groups, the slope through the origin of the challenge's E9.5 contrast (g against h) on Feng's E9.5
      contrast, over the genes both groups keep; L = the median over pairs; beside it E10.5 library 3 on Feng E14.5.
Every value has a bootstrap over the challenge's cells (200 draws, cells resampled within label, stage and library;
Feng's statistics are fixed) with its 5th and 95th percentiles.

The reading (fixed in row T47 before the run; written with the numbers): A = combined share of (b) libraries 1 and 2,
averaged; K = combined share of (b) library 3 - A; P = combined share of (a); L = the scale of (d). Statements on the
kind (K against 0.1), the pace (A / L against 0.2 and 0.3) and the entry's fraction 0.4 against f_max = L - A,
f_even = (L - A) / 2 and f_decay = L (x - x^3) / (1 - x^5) with (1 - x) / (1 - x^5) = A / L. The read puts nothing
forward by itself; its result can only go to the team captain. The output holds shares, slopes, correlations, gene
numbers and percentiles only: no cell count.

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
BOOT_SEED = (47, 0)  # row T47's stream: no existing step draws on it
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
SCALE_SETS = {"E9.5": (("E9.5", None), BRACKET[0]),  # (challenge cells, Feng stage)
              "E10.5_lib3_on_feng_E14.5": (("E10.5", (LIB_HEART,)), BRACKET[1])}
COMBINED_FIELDS = ("share", "share_corrected")
SHARE_FIELDS = ("share", "pearson", "slope_intercept", "share_top_half", "share_corrected")

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
    (ceil(n / 2) genes, ties in gene order) and the share over the N48 scale."""
    x, y = np.asarray(d, np.float64)[mask], np.asarray(change, np.float64)[mask]
    n = len(x)
    nan = float("nan")
    out = dict.fromkeys(SHARE_FIELDS, nan)
    if n == 0 or float(x @ x) == 0.0:
        return out
    out["share"] = float(y @ x / (x @ x))
    if n >= 3 and np.std(x) > 0 and np.std(y) > 0:
        out["pearson"] = float(np.corrcoef(x, y)[0, 1])
    if n >= 2 and np.var(x) > 0:
        out["slope_intercept"] = float(np.mean((x - x.mean()) * (y - y.mean())) / np.var(x))
    top = np.argsort(-np.abs(x), kind="stable")[: n - n // 2]
    if float(x[top] @ x[top]) > 0:
        out["share_top_half"] = float(y[top] @ x[top] / (x[top] @ x[top]))
    if scale is not None and np.isfinite(scale) and scale > 0:
        out["share_corrected"] = out["share"] / float(scale)
    return out


def origin_slope(y: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    """Least-squares slope through the origin of y on x, and their Pearson correlation (NaN where undefined)."""
    nan = float("nan")
    if len(x) == 0 or float(x @ x) == 0.0:
        return nan, nan
    r = float(np.corrcoef(x, y)[0, 1]) if len(x) >= 3 and np.std(x) > 0 and np.std(y) > 0 else nan
    return float(y @ x / (x @ x)), r


def combine(values: Mapping[str, float], weights: Mapping[str, float]) -> tuple[float, float]:
    """(weighted mean, unweighted mean) over the groups with a finite value; the weighted one over those of them with a
    positive weight, renormalised; NaN where none."""
    ok = [g for g, v in values.items() if np.isfinite(v)]
    unw = float(np.mean([values[g] for g in ok])) if ok else float("nan")
    wok = [g for g in ok if weights.get(g, 0.0) > 0]
    tot = float(sum(weights[g] for g in wok))
    wtd = float(sum(weights[g] * values[g] for g in wok) / tot) if tot > 0 else float("nan")
    return wtd, unw


def decay_x(ratio: float) -> float:
    """x > 0 with (1 - x) / (1 - x^5) = ratio, i.e. 1 + x + x^2 + x^3 + x^4 = 1 / ratio (x = 1 at 0.2, x = 0 at 1);
    NaN outside 0 < ratio <= 1."""
    if not np.isfinite(ratio) or ratio <= 0 or ratio > 1:
        return float("nan")
    if ratio == 1:
        return 0.0
    roots = np.roots([1.0, 1.0, 1.0, 1.0, 1.0 - 1.0 / ratio])
    real = [float(r.real) for r in roots if abs(r.imag) < 1e-9 and r.real > 0]
    return min(real, key=lambda r: abs(1 + r + r ** 2 + r ** 3 + r ** 4 - 1 / ratio))


def courses(A: float, L: float) -> dict:
    """The test's applied fraction under the two simple courses and the bound, from the day's share A and the scale L:
    f_max = L - A, f_even = (L - A) / 2 and f_decay = L (x - x^3) / (1 - x^5), written L x (1 + x) / (1 + x + x^2 +
    x^3 + x^4) (the same, without the 0/0 at x = 1)."""
    nan = float("nan")
    ratio = A / L if np.isfinite(A) and np.isfinite(L) and L != 0 else nan
    x = decay_x(ratio)
    f_decay = L * x * (1 + x) / (1 + x + x ** 2 + x ** 3 + x ** 4) if np.isfinite(x) else nan
    return {"A_over_L": ratio, "x": x, "f_max": L - A, "f_even": (L - A) / 2.0, "f_decay": f_decay}


# ---------------------------------------------------------------------------------------------------------------------
# One measurement (the point or one bootstrap draw)


def measure(tabs: Mapping[str, dict], sums: Mapping[str, np.ndarray], setup: Mapping) -> dict:
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
        slopes = []
        for g, h in setup["pairs"]:
            a, b = mean_of(spec, h), mean_of(spec, g)
            fx = setup["feng_contrast"][sname].get((g, h))
            m = masks[g] & masks[h]
            if a is None or b is None or fx is None:
                s, r = float("nan"), float("nan")
            else:
                s, r = origin_slope(dl.lfc_ratio(a, b, eps)[m], fx[m])
            out[f"scale_{sname}/{g}~{h}/slope"] = s
            out[f"scale_{sname}/{g}~{h}/pearson"] = r
            slopes.append(s)
        fin = [s for s in slopes if np.isfinite(s)]
        out[f"scale_{sname}/median_slope"] = float(np.median(fin)) if fin else float("nan")

    cw = "combined_weighted/share"
    A = (out[f"b_E9.5_E10.5_lib1/{cw}"] + out[f"b_E9.5_E10.5_lib2/{cw}"]) / 2.0
    L = out["scale_E9.5/median_slope"]
    out["reading/A"] = A
    out["reading/K"] = out[f"b_E9.5_E10.5_lib3/{cw}"] - A
    out["reading/P"] = out[f"a_E8.5_E9.5/{cw}"]
    out["reading/L"] = L
    for k, v in courses(A, L).items():
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
    """The statement on the entry's fraction: "above the bound" if it is above f_max's upper percentile; else against
    the range from the lower percentile of the smaller course (by point value) to the upper percentile of the
    larger."""
    vals = [f_max["p95"], f_even["value"], f_decay["value"], f_even["p5"], f_even["p95"], f_decay["p5"],
            f_decay["p95"]]
    if not all(np.isfinite(v) for v in vals):
        return NOT_COMPUTABLE
    if fraction > f_max["p95"]:
        return ABOVE_BOUND
    small, large = (f_even, f_decay) if f_even["value"] <= f_decay["value"] else (f_decay, f_even)
    lo, hi = small["p5"], large["p95"]
    if lo <= fraction <= hi:
        return INSIDE
    return ABOVE_BOTH if fraction > hi else BELOW_BOTH


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


def reading(stages: Mapping[str, ad.AnnData], feng: dict, held: Sequence[str], scales: Mapping[str, float],
            cfg: Mapping, n_boot: int = N_BOOT, seed=BOOT_SEED) -> dict:
    """The whole read on loaded inputs: T1 stages (stage -> AnnData), Feng's statistics (only the bracket's two ends),
    M1's held list and each group's N48 scale."""
    label_key = cfg.get("label_key", "celltype")
    genes = list(stages[GENE_STAGE].var_names)
    dirs, _, idx = feng_dirs.directions_from_stats(feng, genes, cfg, BRACKET)
    groups = sorted(dirs)
    held_set = set(map(str, held))
    not_held = np.array([g not in held_set for g in genes], dtype=bool)
    masks = {g: np.isfinite(dirs[g]) & not_held for g in groups}
    labels, source = group_labels_by_stage(stages, groups, cfg, label_key)
    keep = {s: {lab for g in groups for lab in labels[s][g]} for s in stages}
    tabs = {s: stage_table(stages[s], genes, keep[s], label_key) for s in stages}
    eps = cfg["metrics"]["pseudocount_cp10k"]
    mc = cfg["metrics"]["min_cells_t1"]
    pairs = list(combinations(groups, 2))
    feng_contrast = {}
    for sname, (_, fstage) in SCALE_SETS.items():
        feng_contrast[sname] = {}
        for g, h in pairs:
            fg = dl.feng_sel(feng, fstage, cfg["external_labels"][g][0])
            fh = dl.feng_sel(feng, fstage, cfg["external_labels"][h][0])
            feng_contrast[sname][(g, h)] = None if fg is None or fh is None else dl.on_t1(dl.lfc_ratio(fh, fg, eps),
                                                                                          idx)
    heart = tabs.get("E10.5")
    weights = {}
    for g in groups:
        labs = set(labels.get("E10.5", {}).get(g) or [])
        weights[g] = (float(sum(n for (lab, li), n in zip(heart["keys"], heart["n_block"], strict=True)
                                if lab in labs and li == LIB_HEART)) if heart is not None else 0.0)
    setup = {"eps": eps, "min_cells": mc, "dirs": dirs, "masks": masks, "labels": labels,
             "scales": {g: scales.get(g) for g in groups}, "weights": weights, "pairs": pairs,
             "feng_contrast": feng_contrast}
    point = measure(tabs, {s: block_sums(t) for s, t in tabs.items()}, setup)
    rng = np.random.default_rng(list(seed))
    draws = {k: np.empty(n_boot) for k in point}
    t0 = time.time()
    for i in range(n_boot):
        sums = {s: block_sums(t, resample_weights(t, rng)) for s, t in tabs.items()}
        m = measure(tabs, sums, setup)
        for k, v in m.items():
            draws[k][i] = v
    seconds = time.time() - t0
    values = {}
    for k, v in point.items():
        fin = draws[k][np.isfinite(draws[k])]
        lo, hi = (np.percentile(fin, PERCENTILES) if fin.size else (np.nan, np.nan))
        values[k] = {"value": float(v), "p5": float(lo), "p95": float(hi),
                     "nonfinite_draws": int(n_boot - fin.size)}
    r = {k.split("/", 1)[1]: values[k] for k in values if k.startswith("reading/")}
    statements = {"kind": kind_statement(*_interval(r["K"])), "pace": pace_statement(*_interval(r["A_over_L"])),
                  "entry_fraction": fraction_statement(r["f_max"], r["f_even"], r["f_decay"])}
    pair_genes = {f"{g}~{h}": int((masks[g] & masks[h]).sum()) for g, h in pairs}
    return {"groups": groups, "genes_per_group": {g: int(masks[g].sum()) for g in groups},
            "genes_per_pair": pair_genes, "label_source": source,
            "n48_scale": {g: scales.get(g) for g in groups}, "values": values, "reading": r,
            "statements": statements, "n_boot": n_boot, "seed": list(seed), "percentiles": list(PERCENTILES),
            "bootstrap_seconds": round(seconds, 1)}


# ---------------------------------------------------------------------------------------------------------------------
# Inputs and the run


def feng_inputs(path: Path | None = None) -> tuple[dict, dict]:
    """Feng's clean statistics (manifest-checked) cut to the bracket's two ends (feng_dirs.bracket_stats): no other
    Feng stage reaches the read."""
    path = feng_dirs.clean_stats_path() if path is None else Path(path)
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path, fields=dl.STAT_FIELDS)
    return feng_dirs.bracket_stats(st, BRACKET), {"file": gv._rel(path), "sha256": entry["sha256"]}


def n48_scales(path: Path = N48_JSON) -> dict[str, float]:
    """Each group's scalar James-Stein scale as N48's rule read it (the 'js' record, the held-flat genes included)."""
    rep = json.loads(Path(path).read_text())
    return {g: r["js"]["scale"] for g, r in rep["groups"].items() if r["js"]["scale"] is not None}


def _fmt(rec: Mapping) -> str:
    return f"{rec['value']:+.3f} [{rec['p5']:+.3f}, {rec['p95']:+.3f}]"


def text(res: Mapping) -> str:
    r = res["read"]
    rd = r["reading"]
    lines = [res["label"],
             f"Kind: {r['statements']['kind']} (K = {_fmt(rd['K'])}, against {KIND_AT})",
             f"Pace: {r['statements']['pace']} (A / L = {_fmt(rd['A_over_L'])}, against {PACE_FIFTH} and {PACE_UPPER})",
             f"Entry's fraction {ENTRY_FRACTION}: {r['statements']['entry_fraction']} (f_max {_fmt(rd['f_max'])}; "
             f"f_even {_fmt(rd['f_even'])}; f_decay {_fmt(rd['f_decay'])}; two courses among many, not a bracket)",
             f"A {_fmt(rd['A'])}; K {_fmt(rd['K'])}; P {_fmt(rd['P'])}; L {_fmt(rd['L'])} "
             f"(value [p{PERCENTILES[0]:g}, p{PERCENTILES[1]:g}], {r['n_boot']} draws)",
             "Groups (genes kept; N48 scale; label source E8.5 / E9.5 / E10.5): " + "; ".join(
                 f"{g} ({r['genes_per_group'][g]}; {r['n48_scale'][g]}; "
                 + " / ".join(r["label_source"].get(s, {}).get(g, "-") for s in T1_STAGES) + ")"
                 for g in r["groups"])]
    v = r["values"]
    for name in COMPARISONS:
        lines.append(f"{name}: combined weighted {_fmt(v[f'{name}/combined_weighted/share'])}, unweighted "
                     f"{_fmt(v[f'{name}/combined_unweighted/share'])}")
        for g in r["groups"]:
            lines.append(f"    {g}: share {_fmt(v[f'{name}/{g}/share'])}, r {v[f'{name}/{g}/pearson']['value']:+.3f}, "
                         f"slope+int {v[f'{name}/{g}/slope_intercept']['value']:+.3f}, top half "
                         f"{v[f'{name}/{g}/share_top_half']['value']:+.3f}, corrected "
                         f"{v[f'{name}/{g}/share_corrected']['value']:+.3f}")
    for sname in SCALE_SETS:
        lines.append(f"scale {sname}: median slope {_fmt(v[f'scale_{sname}/median_slope'])}")
        for pair, n in r["genes_per_pair"].items():
            lines.append(f"    {pair} ({n} genes): slope {_fmt(v[f'scale_{sname}/{pair}/slope'])}, r "
                         f"{v[f'scale_{sname}/{pair}/pearson']['value']:+.3f}")
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
    scales = n48_scales()
    stages = {s: ad.read_h5ad(files[s]) for s in T1_STAGES}
    read = reading(stages, feng, held, scales, cfg)
    res = {"label": LABEL_T47, "read": read, "feng": feng_rec, "held_genes": len(held)}
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
