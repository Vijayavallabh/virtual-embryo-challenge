#!/usr/bin/env python
"""Task 3 response decomposition under the official metrics (pre-registered for the Gata4 release, 2026-10-20).

Question. Which parts of a knockout-vs-matched-WT MERFISH response carry signal under the official Task 3 metrics,
and which of them carry over from one knockout embryo to another? Five components, each built from real wild-type
cells in count space, are scored alone and in every combination (a full factorial, summarised by Shapley values):

  cov     section coverage: the knockout's strict non-Mesp1 label shares (endoderm, surface ectoderm, neural tube,
          neural crest; branchial arch counts as mixed); every other label keeps its WT relative share
  mesp1   lineage composition: the knockout's relative shares of every other label
  stage   stage axis: within each label a fraction of the cells comes from other WT stages (oracle NNLS mixing)
  depth   per-label sequencing depth of the knockout
  factor  one per-gene factor shared across cell types, estimated from the strict non-Mesp1 compartments only

Compositions come from vec.harmonise (knockout cells relabelled into the WT vocabulary), as a top-call/EM bracket,
bias-matched against the WT's own out-of-fold calls. Prediction pools come from vec.prep.load_stage with
drop_duplicates=True; the emulator (vec.evaluate: official veckit metrics, 10% working copies, split-half ceiling,
computed anchors) scores against the target and the reference as released. Oracle components are re-estimated for
every seed from knockout cells outside that seed's working copy. On Gata4, components borrowed from the other
replicate (the construction a beta-catenin prediction can use), the Mab21l2 factor and the WT-series factor R are
scored as well, and the rules in rules.yaml turn the scores into the decision file: which components enter
submissions 1 and 2, with which weights.

Every analysis choice and every rule is in rules.yaml next to this file; only file paths come from the command line.
No board anchor, boards/index.json or leaderboard score is read.

Commands (DGX, repository root; about 10 min for mab21l2 and 20-30 min for gata4 with 8 workers x 2 threads):

  uv run python experiments/t3_decomp/t3_decompose.py mab21l2 --out runs/t3_decomp/mab21l2
  uv run python experiments/t3_decomp/t3_decompose.py gata4 REP1.h5ad REP2.h5ad data/raw/T2/E8.75.h5ad \\
      --out runs/t3_decomp/gata4

gata4 also reads WT E9.5, WT E8.25_late and the Mab21l2 KO (defaults under data/raw/; see --help) and adds WT E8.5 to
the R series when data/raw/T2/E8.5.h5ad (or --e85) exists. It exits with status 0 and a message, without loading
anything, when REP1, REP2 or WT E8.75 is missing (before the release). run_gata4.sh wraps it with logging.
"""

from __future__ import annotations

import os

os.environ.setdefault("POT_BACKEND_DISABLE_JAX", "1")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

import argparse  # noqa: E402
import ast  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import re  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import zlib  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from itertools import product  # noqa: E402
from math import factorial  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy.optimize import nnls  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RULES_PATH = HERE / "rules.yaml"
SCENARIOS_PATH = REPO / "configs" / "t3_scenarios.yaml"
DATA = REPO / "data" / "raw"
TARGET_SUM = 1e4
EST = ("top", "em")
EST_COL = {"top": "top_call", "em": "em"}
COMPONENTS = ("cov", "mesp1", "stage", "depth", "factor")
WT_LIKE = re.compile(r"(?:^|[^a-z0-9])(?:wt|wild[ _-]?type|control|ctrl)(?:[^a-z0-9]|$)", re.IGNORECASE)
VEC_MODULES = ("evaluate", "prep", "harmonise", "spatial", "boards", "io")
# Names of the prototype (screened report, 2026-10-18) for the same constructions, for side-by-side tables.
PROTOTYPE_NAMES = {
    "fac:none": "wt_resample", "fac[top]:cov": "wt_cov", "fac[top]:mesp1": "wt_mesp1",
    "fac[top]:cov+mesp1": "wt_comp", "fac:stage": "wt_stage", "fac[top]:cov+mesp1+stage": "wt_comp_stage",
    "fac:depth": "wt_dm", "fac[top]:cov+depth": "wt_cov_dm", "fac[top]:cov+mesp1+depth": "wt_comp_dm",
    "fac[top]:cov+mesp1+stage+depth": "wt_comp_stage_dm", "fac:factor": "wt_gf_nm",
    "fac[top]:cov+factor": "wt_cov_gf_nm", "fac[top]:cov+depth+factor": "wt_cov_gf_nm_dm",
    "fac[top]:cov+mesp1+factor": "wt_comp_gf_nm", "fac[top]:cov+mesp1+stage+factor": "wt_comp_stage_gf_nm",
    "fac[top]:cov+mesp1+stage+depth+factor": "wt_comp_stage_gf_nm_dm", "oracle:comp+factor_all": "wt_comp_gf_all",
    "oracle:ko_draw": "ko_draw", "oracle:ko_within": "ko_within", "oracle:ko_within_wtdepth": "ko_within_dm",
    "plan_N_as_written": "ko_sys", "cov_broad[top]": "wt_cov_broad", "R": "wt_R",
}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# Count-space operations (dense int32 counts, cells x genes)


def lognorm(C) -> np.ndarray:
    """log1p CP10k of dense counts as float32 (the released convention; same values as vec.prep.log1p_cp10k)."""
    C = np.asarray(C, dtype=np.float64)
    n = C.sum(1, keepdims=True)
    return np.log1p(np.divide(TARGET_SUM * C, n, out=np.zeros_like(C), where=n > 0)).astype(np.float32)


def thin(C, p, rng) -> np.ndarray:
    """Binomial thinning with per-cell keep probability p; a cell that would become empty keeps its counts."""
    p = np.clip(np.broadcast_to(np.asarray(p, dtype=np.float64).reshape(-1, 1), (C.shape[0], 1)), 0.0, 1.0)
    out = rng.binomial(C.astype(np.int64), p).astype(np.int32)
    empty = out.sum(1) == 0
    out[empty] = C[empty]
    return out


def rescale_depth(C, ratio, P, rng) -> np.ndarray:
    """Per-cell depth change by `ratio`: thinning below 1, Poisson top-up from smoothed proportions P above 1."""
    ratio = np.asarray(ratio, dtype=np.float64)
    out = thin(C, np.minimum(ratio, 1.0), rng)
    up = ratio > 1.0
    if up.any():
        n = C[up].sum(1).astype(np.float64)
        out[up] = out[up] + rng.poisson(P[up].astype(np.float64) * (n * (ratio[up] - 1.0))[:, None]).astype(np.int32)
    return out


def apply_gene_factor(C, phi, P, rng) -> np.ndarray:
    """Per-gene factor in count space: binomial thinning where phi < 1, Poisson top-up from P where phi > 1."""
    out = C.copy()
    down = phi < 1.0
    if down.any():
        out[:, down] = rng.binomial(C[:, down].astype(np.int64), phi[down]).astype(np.int32)
    up = phi > 1.0
    if up.any():
        n = C.sum(1, keepdims=True).astype(np.float64)
        out[:, up] += rng.poisson(P[:, up].astype(np.float64) * n * (phi[up] - 1.0)).astype(np.int32)
    empty = out.sum(1) == 0
    out[empty] = C[empty]
    return out


def smoothed_props(C, n_pcs=30, k=15, seed=0, chunk=4000) -> np.ndarray:
    """kNN-smoothed within-cell proportions (mean over the k nearest cells in a PCA of log1p CP10k)."""
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    X = lognorm(C)
    Z = PCA(n_components=min(n_pcs, X.shape[1] - 1, X.shape[0] - 1), random_state=seed).fit_transform(X)
    _, ind = NearestNeighbors(n_neighbors=min(k, X.shape[0])).fit(Z).kneighbors(Z)
    Pc = C / np.maximum(C.sum(1, keepdims=True), 1)
    P = np.empty(C.shape, dtype=np.float32)
    for i in range(0, C.shape[0], chunk):
        P[i:i + chunk] = Pc[ind[i:i + chunk]].mean(1)
    return P


def alloc(n, f) -> np.ndarray:
    """Largest-remainder allocation of n cells to shares f."""
    raw = n * np.asarray(f, dtype=np.float64)
    k = np.floor(raw).astype(int)
    rem = int(n - k.sum())
    if rem > 0:
        k[np.argsort(-(raw - k))[:rem]] += 1
    return k


def draw_rows(pool_idx, k, rng) -> np.ndarray:
    m = len(pool_idx)
    if k <= 0 or m == 0:
        return np.array([], dtype=int)
    if k <= m:
        return rng.choice(pool_idx, k, replace=False)
    return np.concatenate([pool_idx, rng.choice(pool_idx, k - m, replace=True)])


def indicator(lab, L) -> sparse.csr_matrix:
    ok = lab >= 0
    return sparse.csr_matrix((np.ones(int(ok.sum())), (lab[ok], np.flatnonzero(ok))), shape=(L, len(lab)))


def type_means(X, lab, L) -> np.ndarray:
    """Per-label mean of the rows of X (labels -1 are left out; empty labels give zeros)."""
    cnt = np.bincount(lab[lab >= 0], minlength=L).astype(np.float64)
    return np.asarray(indicator(lab, L) @ np.asarray(X, dtype=np.float64)) / np.maximum(cnt, 1)[:, None]


def props_type_means(C, lab, L) -> np.ndarray:
    """Per-label mean within-cell proportions c / n (depth-free in expectation)."""
    return type_means(C / np.maximum(C.sum(1, keepdims=True), 1), lab, L)


def label_medians(lab, values, L) -> np.ndarray:
    s = pd.Series(np.asarray(values, dtype=np.float64)[lab >= 0]).groupby(lab[lab >= 0]).median()
    out = np.full(L, np.nan)
    out[s.index.to_numpy()] = s.to_numpy()
    return out


def gene_factor(mu_ko, mu_wt, w, mask, eps=1e-5, clip=(0.25, 4.0)) -> tuple[np.ndarray, np.ndarray]:
    """Shared per-gene factor: ratio of w-weighted mean proportions over the labels in `mask`, rescaled so the
    weighted total proportion is unchanged, clipped. Returns (phi, denominator)."""
    ww = (np.asarray(w, dtype=np.float64) * np.asarray(mask, dtype=np.float64))[:, None]
    num, den = (ww * mu_ko).sum(0), (ww * mu_wt).sum(0)
    if den.sum() <= 0:
        return np.ones(mu_wt.shape[1]), den
    phi = (num + eps) / (den + eps)
    phi = phi / (den @ phi / max(den.sum(), 1e-12))
    return np.clip(phi, *clip), den


def stable_rng(seed: int, name: str) -> np.random.Generator:
    """Generator for one (seed, variant): reproducible whatever the order or the number of workers."""
    return np.random.default_rng([int(seed), zlib.crc32(name.encode())])


def working_copy_rows(n: int, frac: float, seed: int) -> np.ndarray:
    """The rows vec.evaluate.subsample(stage, frac, seed) keeps (checked against it at start-up)."""
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, max(2, int(round(n * frac))), replace=False))


def row_hashes_layer(a) -> np.ndarray:
    """64-bit content hash of every row of layers['counts'] (vec.spatial.row_hashes on the exact float64 values, so
    integer and half-integer records hash as stored)."""
    from vec.spatial import row_hashes

    M = sparse.csr_matrix(a.layers["counts"], dtype=np.float64)
    M.sum_duplicates()
    M.eliminate_zeros()
    return row_hashes(M)


# ---------------------------------------------------------------------------------------------------------------------
# Small helpers


def sha256(path, chunk=1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def parse_seeds(s: str) -> list[int]:
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",") if x.strip()]


def json_safe(o):
    if isinstance(o, dict):
        return {str(k): json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return json_safe(o.tolist())
    return o


_ALLOWED_NODES = (ast.Expression, ast.BoolOp, ast.And, ast.Or, ast.UnaryOp, ast.Not, ast.Name, ast.Load, ast.Constant)


def eval_condition(expr, names: dict) -> bool:
    """A rule condition such as "D2 and not D10": and/or/not over rule names (true/false allowed). An outcome that
    could not be evaluated (None) counts as False; an unknown name is an error."""
    if isinstance(expr, bool):
        return expr
    tree = ast.parse(str(expr).strip(), mode="eval")
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES) or (isinstance(node, ast.Constant)
                                                   and not isinstance(node.value, bool)):
            raise ValueError(f"condition {expr!r}: only and/or/not over rule names are allowed")

    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.BoolOp):
            vals = [ev(v) for v in n.values]
            return all(vals) if isinstance(n.op, ast.And) else any(vals)
        if isinstance(n, ast.UnaryOp):
            return not ev(n.operand)
        if isinstance(n, ast.Constant):
            return bool(n.value)
        if n.id.lower() in ("true", "false"):
            return n.id.lower() == "true"
        if n.id not in names:
            raise KeyError(f"unknown name {n.id!r} in condition {expr!r}")
        return bool(names[n.id])

    return ev(tree)


def first_match(table: list[dict], names: dict) -> dict:
    for row in table:
        if eval_condition(row["when"], names):
            return row
    raise ValueError("no row matched (every table must end with when: true)")


# ---------------------------------------------------------------------------------------------------------------------
# Data


@dataclass
class Vocab:
    """The reference (matched WT) label vocabulary with lineage flags and the WT pool composition."""

    labels: list[str]
    group: np.ndarray
    mesp1: np.ndarray
    nm: np.ndarray            # strict non-Mesp1 (bool)
    broad: np.ndarray         # not Mesp1 "yes" (bool)
    model_idx: np.ndarray     # positions of the classifier's classes, in model.classes order
    p: np.ndarray             # WT pool annotation composition
    count: np.ndarray

    @property
    def L(self) -> int:
        return len(self.labels)

    def index(self, names) -> np.ndarray:
        pos = {t: i for i, t in enumerate(self.labels)}
        return np.array([pos[n] for n in names if n in pos], dtype=int)


@dataclass
class Pool:
    """Cells to resample from: counts, label index in the vocabulary (-1: none), smoothed proportions."""

    C: np.ndarray
    lab: np.ndarray
    P: np.ndarray | None
    L: int
    idx: list = field(default_factory=list)

    def __post_init__(self):
        order = np.argsort(self.lab, kind="stable")
        bounds = np.searchsorted(self.lab[order], np.arange(self.L + 1))
        self.idx = [order[bounds[t]:bounds[t + 1]] for t in range(self.L)]

    @property
    def count(self) -> np.ndarray:
        return np.array([len(i) for i in self.idx])


def add_label_column(a, mode: str) -> None:
    ct = a.obs["celltype"].astype(str).to_numpy() if "celltype" in a.obs else np.array(["Unknown"] * a.n_obs)
    lab = ct
    if mode == "combined" and "cm_celltype" in a.obs:
        cm = a.obs["cm_celltype"].astype(str).to_numpy()
        lab = np.where(cm != "Unknown", cm, ct)
    a.obs["t3_label"] = pd.Categorical(lab)


def littermate_split(a) -> tuple[np.ndarray, str | None, dict]:
    """Cells marked wild type by obs['condition'] or obs['genotype'] in a file that also holds other levels."""
    for col in ("condition", "genotype"):
        if col in a.obs:
            lev = a.obs[col].astype(str)
            wt = lev.str.contains(WT_LIKE).to_numpy()
            if wt.any() and (~wt).any():
                return wt, col, lev.value_counts().to_dict()
            return np.zeros(a.n_obs, dtype=bool), col, lev.value_counts().to_dict()
    return np.zeros(a.n_obs, dtype=bool), None, {}


def dense_counts(a) -> np.ndarray:
    """Dense int32 counts of a prediction pool. Deduplicated pools hold no interpolated (non-integer) records; any
    non-integer value left is rounded and counted in uns['t3_rounded_entries']."""
    C = a.layers["counts"]
    C = C.toarray() if sparse.issparse(C) else np.asarray(C)
    if np.issubdtype(C.dtype, np.floating):
        frac = np.abs(C - np.rint(C)) > 1e-6
        a.uns["t3_rounded_entries"] = int(frac.sum())
        C = np.rint(C)
    return C.astype(np.int32)


def file_record(path, a) -> dict:
    qc = a.uns.get("qc", {})
    return {"path": str(path), "sha256": sha256(path), "n_cells": int(a.n_obs),
            "n_in_file": int(qc.get("n_input", a.n_obs)), "duplicates_dropped": int(qc.get("n_dropped", 0)),
            **{k: int(qc[k]) for k in ("n_copies", "n_interpolated") if k in qc},
            "unmapped_labels": list(map(str, a.uns.get("unmapped_labels", [])))}


# ---------------------------------------------------------------------------------------------------------------------
# Harmonisation (vec.harmonise, settings frozen in rules.yaml)


class Harmoniser:
    """Label transfer into the reference vocabulary. One model per reference depth: a reference deeper than the
    query is refitted thinned to the query's median depth; a deeper query is thinned by `predict` (5 draws)."""

    def __init__(self, ref, label_key, group_map, hp, n_jobs):
        from vec.harmonise import median_depth

        self.ref, self.key, self.gm, self.hp, self.n_jobs = ref, label_key, group_map, hp, n_jobs
        self.ref_depth = median_depth(ref)
        self.models, self.oofs, self.oof_pos = {}, {}, {}

    def depth_key(self, query):
        from vec.harmonise import median_depth

        d = median_depth(query)
        if d is None or self.ref_depth is None or d >= self.ref_depth * (1.0 - self.hp["depth_tol"]):
            return None
        return float(round(d, 1))

    def model(self, key):
        from vec.harmonise import fit_transfer

        if key not in self.models:
            t0 = time.time()
            hp = self.hp
            self.models[key] = fit_transfer(
                self.ref, self.key, classifier=hp["classifier"], n_pcs=hp["n_pcs"], C=hp["C"], clip=hp["clip"],
                threshold=hp["threshold"], cv_folds=hp["cv_folds"], exclude_labels=tuple(hp["exclude_labels"]),
                min_class_cells=hp["min_class_cells"], group_map=self.gm, depth=key, depth_tol=hp["depth_tol"],
                seed=hp["seed"], n_jobs=self.n_jobs)
            log(f"  harmoniser fitted (reference depth target {key}): {len(self.models[key].classes)} classes, "
                f"{time.time() - t0:.0f}s")
        return self.models[key]

    def oof(self, key):
        """WT out-of-fold calls (vec.harmonise.cross_validate) with the model's settings, and their pool rows."""
        from vec.harmonise import cross_validate, transfer_metrics

        if key not in self.oofs:
            t0 = time.time()
            hp = self.hp
            df = cross_validate(
                self.ref, self.key, method=hp["classifier"], n_folds=hp["cv_folds"], n_pcs=hp["n_pcs"], C=hp["C"],
                clip=hp["clip"], threshold=hp["threshold"], exclude_labels=tuple(hp["exclude_labels"]),
                min_class_cells=hp["min_class_cells"], group_map=self.gm, depth=key, depth_tol=hp["depth_tol"],
                novelty=False, seed=hp["seed"], n_jobs=self.n_jobs)
            m = self.model(key)
            cols = [f"p:{c}" for c in m.classes]
            missing = [c for c in cols if c not in df.columns]
            if missing:
                raise RuntimeError(f"out-of-fold classes differ from the model's: {missing[:5]}")
            self.oofs[key] = df
            self.oof_pos[key] = self.ref.obs_names.get_indexer(df.index)
            tm = transfer_metrics(df)
            log(f"  WT out-of-fold calls (depth target {key}): accuracy {tm['accuracy']:.3f}, macro-F1 "
                f"{tm['macro_f1']:.3f}, {time.time() - t0:.0f}s")
        return self.oofs[key]

    def calls(self, query, seed):
        from vec.harmonise import predict

        key = self.depth_key(query)
        df = predict(self.model(key), query, threshold=self.hp["threshold"], depth_match=True,
                     depth_draws=self.hp["depth_draws"], depth_tol=self.hp["depth_tol"], seed=seed, proba=True)
        return df, key


def bracket(calls, model, rows=None) -> dict:
    """Top-call and EM label shares (vec.harmonise.composition_bracket), in model.classes order."""
    from vec.harmonise import composition_bracket

    df = calls if rows is None else calls.iloc[rows]
    br = composition_bracket(df, model, level="label", estimators=("top_call", "em"))
    br = br.reindex([str(c) for c in model.classes])
    return {e: br[EST_COL[e]].to_numpy(dtype=np.float64) for e in EST}


def calls_to_index(V: Vocab, calls, conf=None) -> np.ndarray:
    pos = {t: i for i, t in enumerate(V.labels)}
    lab = np.array([pos.get(str(x), -1) for x in calls["label_raw"].to_numpy()])
    if conf is not None:
        lab = np.where(calls["label_conf"].to_numpy() >= conf, lab, -1)
    return lab


# ---------------------------------------------------------------------------------------------------------------------
# Compositions over the full vocabulary


def bias_matched(V: Vocab, q, p_oof, p=None) -> np.ndarray:
    """WT annotation shares plus the knockout-minus-WT difference of call shares (clipped at 0, renormalised over the
    classifier's labels); labels outside the classifier keep their WT shares."""
    p = V.p if p is None else p
    pm = p[V.model_idx]
    s = pm.sum()
    if s <= 0:
        return p.copy()
    fm = np.clip(pm / s + (np.asarray(q) - np.asarray(p_oof)), 0.0, None)
    fm = fm / fm.sum()
    f = p.copy()
    f[V.model_idx] = fm * s
    return f


def coverage_comp(V: Vocab, f, mask=None) -> np.ndarray:
    """Shares of the masked (non-Mesp1) labels from f; every other label keeps its WT relative share."""
    nm = V.nm if mask is None else mask
    out = V.p.copy()
    out[nm] = f[nm]
    rest_p, rest_target = V.p[~nm].sum(), 1.0 - f[nm].sum()
    out[~nm] = V.p[~nm] * (rest_target / rest_p) if rest_p > 0 else 0.0
    return out / out.sum()


def mesp1_comp(V: Vocab, f) -> np.ndarray:
    """Relative shares of every label outside the strict non-Mesp1 set from f; non-Mesp1 shares stay at WT."""
    out = V.p.copy()
    rest_f, rest_target = f[~V.nm].sum(), 1.0 - V.p[V.nm].sum()
    out[~V.nm] = f[~V.nm] * (rest_target / rest_f) if rest_f > 0 else V.p[~V.nm]
    return out / out.sum()


def rule_comp(V: Vocab, rule: dict, lists: dict, log2_magnitude: float | None = None) -> tuple[np.ndarray, list]:
    """A blind rule or scenario on the WT composition: depleted states x (1 - fraction) or x 2^-m, enriched x factor."""
    f = V.p.copy()
    dep = dict(rule.get("deplete", {}))
    for name in lists.get(rule.get("deplete_list"), []):
        dep[name] = rule.get("fraction", 0.5)
    missing = []
    for t, d in dep.items():
        i = V.index([t])
        if not i.size:
            missing.append(t)
            continue
        f[i] *= (2.0 ** -log2_magnitude) if log2_magnitude is not None else (1.0 - d)
    for t, e in rule.get("enrich", {}).items():
        i = V.index([t])
        if not i.size:
            missing.append(t)
            continue
        f[i] *= e
    return f / f.sum(), missing


def within(f, states) -> np.ndarray:
    x = np.asarray(f, dtype=np.float64)[states]
    return x / max(x.sum(), 1e-300)


# ---------------------------------------------------------------------------------------------------------------------
# Component estimates for one knockout pool (rows of its deduplicated cells)


@dataclass
class Estimate:
    q: dict
    f: dict                    # f[est] = {"cov": ..., "mesp1": ..., "comp": ...} over the vocabulary
    q_hard: np.ndarray
    ratio: np.ndarray          # per-label depth ratio knockout / WT
    phi_nm: np.ndarray
    phi_all: np.ndarray
    alphas: np.ndarray
    mu: np.ndarray
    count: np.ndarray
    pool: Pool
    nm_delta: np.ndarray


def estimate(ctx, T, rows, rng) -> Estimate:
    V, wt, cfg = ctx.V, ctx.wt, ctx.rules["components"]
    C, lab = T.C[rows], T.lab[rows]
    model = ctx.harm.model(T.key)
    q = bracket(T.calls, model, rows)
    f = {}
    for e in EST:
        comp = bias_matched(V, q[e], ctx.p_oof[T.key][e])
        f[e] = {"comp": comp, "cov": coverage_comp(V, comp), "mesp1": mesp1_comp(V, comp)}
    count = np.bincount(lab[lab >= 0], minlength=V.L)
    q_hard = count / max(count.sum(), 1)
    dmin = cfg["depth"]["min_cells"]
    med = label_medians(lab, C.sum(1), V.L)
    ok_d = (count >= dmin) & (wt.count >= dmin) & np.isfinite(med) & (wt.med > 0)
    ratio = np.where(ok_d, med / np.where(wt.med > 0, wt.med, 1.0), 1.0)
    fc = cfg["factor"]
    mu = props_type_means(C, lab, V.L)
    ok = (count >= fc["min_cells"]) & (wt.count >= fc["min_cells"])
    phi_nm, _ = gene_factor(mu, wt.mu, V.p, ok & V.nm, fc["eps"], tuple(fc["clip"]))
    phi_all, _ = gene_factor(mu, wt.mu, V.p, ok, fc["eps"], tuple(fc["clip"]))
    sm = ctx.rules["components"]["smoothing"]
    P = smoothed_props(C, sm["n_pcs"], sm["k"])
    pool = Pool(C, lab, P, V.L)
    alphas = np.zeros(len(ctx.stages))
    if ctx.stages:
        C_dm = rescale_depth(C, np.where(ratio > 0, 1.0 / ratio, 1.0)[np.maximum(lab, 0)], P, rng)
        M_dm = type_means(lognorm(C_dm), lab, V.L)
        resid = q_hard @ (M_dm - wt.M)
        S = stage_axes(ctx, q_hard)
        alphas, _ = nnls(S.T, resid)
        if alphas.sum() > 1:
            alphas = alphas / alphas.sum()
    nmk, nmw = V.nm[np.maximum(lab, 0)] & (lab >= 0), V.nm[wt.pool.lab]
    nm_delta = (lognorm(C[nmk]).mean(0).astype(np.float64) - wt.X[nmw].mean(0)) if nmk.any() and nmw.any() \
        else np.zeros(C.shape[1])
    return Estimate(q, f, q_hard, ratio, phi_nm, phi_all, alphas, mu, count, pool, nm_delta)


def stage_axes(ctx, f) -> np.ndarray:
    """Within-type stage axis per stage pool: sum_t f_t (M_stage,t - M_WT,t) over labels with enough stage cells."""
    mn = ctx.rules["components"]["stage"]["min_cells"]
    return np.array([((f * (sp.count >= mn))[:, None] * (sp.M - ctx.wt.M)).sum(0) for sp in ctx.stages])


# ---------------------------------------------------------------------------------------------------------------------
# Prediction builder


def build(ctx, n, rng, f, src=None, src_by_label=None, alphas=None, phi=None, w=1.0, ratio=None) -> np.ndarray:
    """Resample cells to composition f (WT pool unless `src` / `src_by_label`), optionally replace a fraction alpha_j
    of each label's cells with depth-matched cells of stage pool j, apply phi ** w, then per-label depth ratios.
    Returns log1p CP10k (float32)."""
    src = ctx.wt.pool if src is None else src
    pick = (lambda t: src) if src_by_label is None else src_by_label
    f = np.asarray(f, dtype=np.float64).copy()
    avail = np.array([len(pick(t).idx[t]) > 0 for t in range(ctx.V.L)])
    f = np.where(avail, f, 0.0)
    f /= f.sum()
    k = alloc(n, f)
    mn = ctx.rules["components"]["stage"]["min_cells"]
    Cs, Ps, labs = [], [], []
    for t in np.flatnonzero(k):
        kt, used = int(k[t]), 0
        if alphas is not None and np.sum(alphas) > 0:
            for sp, a in zip(ctx.stages, alphas, strict=True):
                if a <= 0 or len(sp.pool.idx[t]) < mn:
                    continue
                kk = min(int(round(kt * a)), kt - used)
                r = draw_rows(sp.pool.idx[t], kk, rng)
                if len(r):
                    Cs.append(sp.pool.C[r])
                    Ps.append(sp.pool.P[r])
                    labs.append(np.full(len(r), t))
                    used += len(r)
        s = pick(t)
        r = draw_rows(s.idx[t], kt - used, rng)
        if len(r):
            Cs.append(s.C[r])
            Ps.append(s.P[r])
            labs.append(np.full(len(r), t))
    C, P, lab = np.vstack(Cs), np.vstack(Ps), np.concatenate(labs)
    if phi is not None:
        C = apply_gene_factor(C, np.asarray(phi) ** w if w != 1.0 else np.asarray(phi), P, rng)
    if ratio is not None:
        C = rescale_depth(C, np.asarray(ratio)[lab], P, rng)
    return lognorm(C)


def fac_name(est: str, on) -> str:
    on = [c for c in COMPONENTS if c in set(on)]
    body = "+".join(on) or "none"
    return f"fac[{est}]:{body}" if {"cov", "mesp1"} & set(on) else f"fac:{body}"


def fac_parse(name: str) -> tuple[str | None, tuple]:
    """Inverse of fac_name: (estimator or None, components)."""
    head, body = name.split(":", 1)
    est = head[4:-1] if head.startswith("fac[") else None
    return est, () if body == "none" else tuple(body.split("+"))


def factorial_cells(comps) -> list[tuple[str | None, tuple]]:
    """(estimator or None, components on) for every factorial cell; cells without a composition component are
    shared by both estimators."""
    out = []
    for bits in product((0, 1), repeat=len(comps)):
        on = tuple(c for c, b in zip(comps, bits, strict=True) if b)
        if {"cov", "mesp1"} & set(on):
            out += [(e, on) for e in EST]
        else:
            out.append((None, on))
    return out


def variant_specs(ctx, ti: int, est: Estimate):
    """(name, builder(n, rng)) for every variant scored against target ti."""
    V, wt, preset = ctx.V, ctx.wt, ctx.preset
    specs = []

    def add(name, **kw):
        specs.append((name, lambda n, rng, kw=kw: build(ctx, n, rng, **kw)))

    for e, on in factorial_cells(ctx.comps):
        fe = est.f[e or "top"]
        if "cov" in on and "mesp1" in on:
            f = fe["comp"]
        elif "cov" in on:
            f = fe["cov"]
        elif "mesp1" in on:
            f = fe["mesp1"]
        else:
            f = V.p
        add(fac_name(e or "top", on), f=f, alphas=est.alphas if "stage" in on else None,
            phi=est.phi_nm if "factor" in on else None, ratio=est.ratio if "depth" in on else None)
    # oracle rows and the plan's component N as first written
    ko = est.pool
    add("oracle:ko_draw", f=est.q_hard, src=ko)
    pm = np.zeros(V.L)
    pm[V.model_idx] = V.p[V.model_idx]
    add("oracle:ko_within", f=pm / pm.sum(), src=ko)
    add("oracle:ko_within_wtdepth", f=pm / pm.sum(), src=ko, ratio=np.where(est.ratio > 0, 1.0 / est.ratio, 1.0))
    add("oracle:factor_all", f=V.p, phi=est.phi_all)
    add("oracle:comp+factor_all", f=est.f["top"]["comp"], phi=est.phi_all)
    add("plan_N_as_written", f=est.f["top"]["cov"], src_by_label=lambda t: ko if V.nm[t] else wt.pool)
    add("cov_broad[top]", f=coverage_comp(V, est.f["top"]["comp"], V.broad))
    for name in preset["blind"]:
        add(f"blind:{name}", f=ctx.blind[name])
        add(f"factor+blind:{name}", f=ctx.blind[name], phi=est.phi_nm)
    if ctx.phi_R is not None:
        add("R", f=V.p, phi=ctx.phi_R)
        add(f"R+blind:{preset['primary_blind']}", f=ctx.blind[preset["primary_blind"]], phi=ctx.phi_R)
        if ctx.scenario_f is not None:
            add("R+scenario@m", f=ctx.scenario_f, phi=ctx.phi_R)
    T = ctx.targets[ti]
    if T.littermate is not None:
        lm = T.littermate
        add("littermate:draw", f=lm.q_hard, src=lm.pool)
        add("littermate:factor", f=V.p, phi=lm.phi)
        add("littermate:factor+comp", f=est.f["top"]["comp"], phi=lm.phi)
    # components borrowed from the other replicate (the construction a beta-catenin prediction can borrow)
    if len(ctx.targets) == 2:
        S = ctx.targets[1 - ti].full
        add("xrep:F", f=V.p, phi=S.phi_nm)
        add("xrep:F+depth", f=V.p, phi=S.phi_nm, ratio=S.ratio)
        if ctx.stages:
            add("xrep:F+stage", f=V.p, phi=S.phi_nm, alphas=S.alphas)
        add("xrep:F+cov", f=S.f["top"]["cov"], phi=S.phi_nm)
        add("xrep:N", f=S.f["top"]["cov"], phi=S.phi_nm, ratio=S.ratio)
        add("xrep:comp+F", f=S.f["top"]["comp"], phi=S.phi_nm)
        add("xrep:draw", f=S.q_hard, src=S.pool)
        add("xrep:sys", f=S.f["top"]["cov"], src_by_label=lambda t, S=S: S.pool if V.nm[t] else wt.pool)
        if ctx.scenario_f is not None:
            add("xrep:F+scenario@m", f=ctx.scenario_f, phi=S.phi_nm)
            add(f"xrep:F+{preset['primary_blind']}@m", f=ctx.blind_m, phi=S.phi_nm)
        W = ctx.rules["weights"]["factor_exponent"]
        for w in W["grid"]:
            add(W["variant_depth_off"].format(w=w), f=V.p, phi=S.phi_nm, w=float(w))
            add(W["variant_depth_on"].format(w=w), f=V.p, phi=S.phi_nm, w=float(w), ratio=S.ratio)
    for name, phi in ctx.transfer.items():
        add(f"transfer:{name}", f=V.p, phi=phi)
    return specs


# ---------------------------------------------------------------------------------------------------------------------
# Scoring one (target, seed): runs in a worker process


_CTX = None


def _diagnostics(Xp, truth, ref, up, dn) -> dict:
    dp = Xp.mean(0, dtype=np.float64) - ref.mean(0, dtype=np.float64)
    dt = truth.mean(0, dtype=np.float64) - ref.mean(0, dtype=np.float64)
    out = {"dp_sd": float(dp.std())}
    idx = np.concatenate([up, dn]).astype(int)
    if len(idx):
        a, b = dp[idx], dt[idx]
        beta = float(b @ a) / max(float(b @ b), 1e-12)
        resid = a - beta * b
        out.update(beta=beta, r2=1.0 - float(resid @ resid) / max(float(((a - a.mean()) ** 2).sum()), 1e-12),
                   sign_agree=float(np.mean(np.sign(a) == np.sign(b))),
                   sign_agree_up=float(np.mean(dp[up] > 0)) if len(up) else np.nan,
                   sign_agree_dn=float(np.mean(dp[dn] < 0)) if len(dn) else np.nan,
                   offset_share=float(a.mean() ** 2 / max(np.mean(a ** 2), 1e-12)))
    return out


def run_task(task):
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=_CTX.threads):
        return score_seed(_CTX, *task)


def score_seed(ctx, ti: int, seed: int) -> list[dict]:
    from common.core_metrics import SEVERITY_LOG_WORST, de_genes

    from vec.boards import RANKED_METRICS
    from vec.evaluate import Stage, aggregate, panel_metrics, skills_from_anchors, split_half

    t0 = time.time()
    T, pr = ctx.targets[ti], ctx.rules["protocol"]
    rows_t = working_copy_rows(len(T.X_rel), pr["working_copy_fraction"], seed)
    rows_r = working_copy_rows(len(ctx.ref_X), pr["working_copy_fraction"], seed + 1)
    truth, other = split_half(Stage(X=T.X_rel[rows_t]), seed=pr["split_seed"])
    ref = Stage(X=ctx.ref_X[rows_r])
    floor_v = panel_metrics("T3", ref, truth, ref, seed=seed)
    ceil_v = panel_metrics("T3", other, truth, ref, seed=seed)
    anchors = {m: {"floor": floor_v.get(m), "ceiling": ceil_v.get(m)} for m in RANKED_METRICS["T3"]}
    anchors = {m: a for m, a in anchors.items() if a["floor"] is not None and a["ceiling"] is not None}
    up, dn, _ = de_genes(truth.X, ref.X)
    n_pred = pr["prediction_cells_multiplier"] * max(truth.n, ref.n)
    pool_rows = np.flatnonzero(~np.isin(T.hash_pool, T.hash_rel[rows_t]))
    est = estimate(ctx, T, pool_rows, stable_rng(seed, f"estimate:{T.name}"))
    floor_cut = SEVERITY_LOG_WORST + pr["severity_floor_tolerance"]
    recs = []

    def record(name, vals, n_cells, X=None):
        skills = skills_from_anchors("T3", vals, anchors)
        score, groups = aggregate("T3", skills)
        sev = vals.get("severity_slope")
        rec = {"target": T.name, "seed": seed, "variant": name, "n_pred": n_cells, "t3": score,
               **{f"skill_{m}": s for m, s in skills.items()}, **{f"group_{g}": s for g, s in groups.items()},
               **{k: vals.get(k) for k in ("de_score", "de_direction", "severity_slope", "mmd_u", "variogram",
                                          "energy_distance", "pb_rel_err", "pseudobulk_pearson", "_de_chance",
                                          "_slope_r2")},
               "sev_floor": bool(sev is None or sev <= floor_cut)}
        if X is not None:
            rec.update(_diagnostics(X, truth.X, ref.X, up, dn))
        recs.append(rec)

    for name, builder in variant_specs(ctx, ti, est):
        X = builder(n_pred, stable_rng(seed, name))
        record(name, panel_metrics("T3", Stage(X=X), truth, ref, seed=seed), X.shape[0], X)
    record("_floor_wt_workcopy", floor_v, ref.n)
    record("_ceiling_split_half", ceil_v, other.n, other.X)
    meta = {"target": T.name, "seed": seed, "variant": "_meta", "n_truth": truth.n, "n_ref": ref.n, "n_up": len(up),
            "n_dn": len(dn), "n_pool": int(len(pool_rows)), "alphas": json.dumps(np.round(est.alphas, 4).tolist()),
            "seconds": round(time.time() - t0, 1)}
    recs.append(meta)
    print(f"[{T.name} seed {seed}] {len(recs) - 3} variants, pool {len(pool_rows)}, truth {truth.n}, DE up {len(up)} "
          f"down {len(dn)}, alphas {np.round(est.alphas, 3).tolist()}, {time.time() - t0:.0f}s", flush=True)
    return recs


# ---------------------------------------------------------------------------------------------------------------------
# Context: everything a worker needs, built once in the parent


@dataclass
class Target:
    name: str
    path: str
    X_rel: np.ndarray
    hash_rel: np.ndarray
    ad: object
    C: np.ndarray
    lab: np.ndarray
    calls: pd.DataFrame
    key: object
    hash_pool: np.ndarray
    littermate: object = None
    full: Estimate | None = None
    record: dict = field(default_factory=dict)


@dataclass
class StagePool:
    name: str
    pool: Pool
    M: np.ndarray
    count: np.ndarray
    record: dict


@dataclass
class WT:
    pool: Pool
    X: np.ndarray
    M: np.ndarray
    mu: np.ndarray
    med: np.ndarray
    count: np.ndarray


@dataclass
class Context:
    preset_name: str
    preset: dict
    rules: dict
    V: Vocab
    wt: WT
    ref_X: np.ndarray
    harm: Harmoniser
    p_oof: dict
    targets: list
    stages: list
    comps: tuple
    blind: dict
    phi_R: np.ndarray | None
    transfer: dict
    scenario_f: np.ndarray | None = None
    blind_m: np.ndarray | None = None
    threads: int = 2
    genes: list = field(default_factory=list)


def label_groups(a, key) -> dict[str, str]:
    df = pd.DataFrame({"label": a.obs[key].astype(str).to_numpy(), "group": a.obs["group"].astype(str).to_numpy()})
    return df.groupby("label")["group"].agg(lambda v: v.value_counts().index[0]).to_dict()


def make_vocab(a, key, model, rules) -> Vocab:
    from vec.prep import lineage_map

    spec, lin = lineage_map()["groups"], rules["lineage"]
    lab = a.obs[key].astype(str).to_numpy()
    vc = pd.Series(lab).value_counts()
    labels = sorted(vc.index, key=lambda t: (-vc[t], t))
    gm = label_groups(a, key)
    group = np.array([gm[t] for t in labels])
    mesp1 = np.array([str(spec.get(g, {}).get("mesp1", "partial")) for g in group])
    nm = np.isin(group, lin["nonmesp1_groups"]) & ~np.isin(labels, lin["mixed_labels"])
    pos = {t: i for i, t in enumerate(labels)}
    count = np.array([vc[t] for t in labels])
    return Vocab(labels=list(labels), group=group, mesp1=mesp1, nm=nm, broad=mesp1 != "yes",
                 model_idx=np.array([pos[str(c)] for c in model.classes]), p=count / count.sum(), count=count)


def git_state() -> dict:
    def run(*cmd):
        try:
            return subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception as e:  # noqa: BLE001
            return f"unavailable: {e!r}"
    return {"commit": run("git", "rev-parse", "HEAD"),
            "status": run("git", "status", "--porcelain", "--", "src/vec", "experiments/t3_decomp", "configs")}


def provenance(paths: dict) -> dict:
    from importlib.metadata import version

    import sklearn

    files = {"script": Path(__file__), "rules": RULES_PATH, "t3_scenarios": SCENARIOS_PATH,
             "lineage_map": REPO / "configs" / "lineage_map.yaml"}
    files.update({f"vec.{m}": REPO / "src" / "vec" / f"{m}.py" for m in VEC_MODULES})
    veckit = None
    try:
        from importlib.metadata import distribution

        veckit = json.loads(distribution("veckit").read_text("direct_url.json") or "{}")
    except Exception as e:  # noqa: BLE001
        veckit = f"unavailable: {e!r}"
    return {"sha256": {k: sha256(v) for k, v in files.items() if Path(v).exists()},
            "inputs": {k: str(v) for k, v in paths.items() if v is not None},
            "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__,
                         "scipy": __import__("scipy").__version__, "sklearn": sklearn.__version__,
                         "anndata": version("anndata"), "veckit": veckit},
            "git": git_state()}


def check_scenario(rules) -> dict:
    """The frozen scenario states must match configs/t3_scenarios.yaml `central`."""
    sc = yaml.safe_load(SCENARIOS_PATH.read_text())["scenarios"]["central"]
    mine = rules["presets"]["gata4"]["scenario"]
    ok = (list(sc["depleted_states"]) == list(mine["depleted_states"])
          and list(sc["spared_states"]) == list(mine["spared_states"]))
    if not ok:
        raise SystemExit("rules.yaml scenario states differ from configs/t3_scenarios.yaml `central`: the "
                         "pre-registration needs a dated note before it can run")
    return {"depleted_states": mine["depleted_states"], "spared_states": mine["spared_states"]}


def build_context(preset_name, paths, rules, n_jobs, threads, out) -> tuple[Context, dict]:
    from vec.harmonise import median_depth
    from vec.prep import load_stage

    preset, hp = rules["presets"][preset_name], rules["harmonisation"]
    meta = {"files": {}, "targets": {}, "stage_pools": {}, "series": {}}
    log(f"reference {paths['wt']}")
    wt_ad = load_stage(paths["wt"], drop_duplicates=True)
    add_label_column(wt_ad, preset["label_mode"])
    genes = list(wt_ad.var_names)
    meta["files"]["wt_pool"] = file_record(paths["wt"], wt_ad)
    ref_ad = load_stage(paths["wt"], genes=genes, drop_duplicates=False)
    ref_X = ref_ad.X.toarray().astype(np.float32)
    meta["files"]["wt_reference_as_released"] = {"n_cells": int(ref_ad.n_obs)}
    del ref_ad
    gm = label_groups(wt_ad, "t3_label")
    harm = Harmoniser(wt_ad, "t3_label", gm, hp, n_jobs)
    model0 = harm.model(None)
    V = make_vocab(wt_ad, "t3_label", model0, rules)
    meta["vocabulary"] = pd.DataFrame({"label": V.labels, "group": V.group, "mesp1": V.mesp1, "nonmesp1": V.nm,
                                       "wt_cells": V.count}).to_dict("records")
    log(f"vocabulary: {V.L} labels, strict non-Mesp1 {[t for t, b in zip(V.labels, V.nm, strict=True) if b]}")
    pos = {t: i for i, t in enumerate(V.labels)}
    wlab = np.array([pos[t] for t in wt_ad.obs["t3_label"].astype(str)])
    Cw = dense_counts(wt_ad)
    sm = rules["components"]["smoothing"]
    wpool = Pool(Cw, wlab, smoothed_props(Cw, sm["n_pcs"], sm["k"]), V.L)
    Xw = lognorm(Cw).astype(np.float64)
    wt = WT(pool=wpool, X=Xw, M=type_means(Xw, wlab, V.L), mu=props_type_means(Cw, wlab, V.L),
            med=label_medians(wlab, Cw.sum(1), V.L), count=wpool.count)
    cache: dict = {}

    def other_sample(name, path, conf=True):
        if path in cache:
            return cache[path]
        a = load_stage(path, genes=genes, drop_duplicates=True)
        calls, key = harm.calls(a, seed=hp["seed"])
        lab = calls_to_index(V, calls, hp["confident_call"] if conf else None)
        C = dense_counts(a)
        rec = {**file_record(path, a), "median_depth": median_depth(a), "harmoniser_depth_target": key,
               "confident_frac": float((lab >= 0).mean())}
        cache[path] = (a, C, lab, calls, rec)
        log(f"  {name}: {a.n_obs} cells ({rec['duplicates_dropped']} duplicates dropped), "
            f"{rec['confident_frac']:.2f} confident")
        return cache[path]

    stages = []
    for name in preset["stage_pools"]:
        a, C, lab, calls, rec = other_sample(name, paths[f"stage:{name}"])
        mn = rules["components"]["depth"]["min_cells"]
        cnt = np.bincount(lab[lab >= 0], minlength=V.L)
        med = label_medians(lab, C.sum(1), V.L)
        ok = (cnt >= mn) & (wt.count >= mn) & np.isfinite(med) & (med > 0)
        r = np.where(ok, wt.med / np.where(med > 0, med, 1.0), 1.0)
        keep = lab >= 0
        P = smoothed_props(C, sm["n_pcs"], sm["k"])
        Cdm = rescale_depth(C[keep], r[lab[keep]], P[keep], stable_rng(0, f"stage:{name}"))
        sp_pool = Pool(Cdm, lab[keep], P[keep], V.L)
        stages.append(StagePool(name, sp_pool, type_means(lognorm(Cdm), lab[keep], V.L), sp_pool.count, rec))
        meta["stage_pools"][name] = rec
    comps = tuple(c for c in COMPONENTS if c != "stage" or stages)

    # blind rules on the WT composition
    lists = {k: v for k, v in preset.items() if isinstance(v, list)}
    blind, missing = {}, {}
    for name, rule in preset["blind"].items():
        blind[name], miss = rule_comp(V, rule, lists)
        if miss:
            missing[name] = miss
    meta["blind_states_absent_from_vocabulary"] = missing

    # component R from the WT series (no knockout data)
    phi_R, meta["R"] = None, None
    series = {n: t for n, t in preset["series"].items() if paths.get(f"series:{n}")}
    if len(series) >= 2:
        sc = rules["components"]["series_factor"]
        mus, cnts, times = [], [], []
        for n, t in series.items():
            a, C, lab, calls, rec = other_sample(n, paths[f"series:{n}"])
            mus.append(props_type_means(C, lab, V.L))
            cnts.append(np.bincount(lab[lab >= 0], minlength=V.L))
            times.append(float(t))
            meta["series"][n] = rec
        ok = np.all([c >= sc["min_cells"] for c in cnts], axis=0) & (wt.count >= sc["min_cells"])
        tt = np.array(times)
        Y = np.log(np.stack(mus) + 1e-6)                          # samples x labels x genes
        A = np.c_[np.ones_like(tt), tt - tt.mean()]
        coef = np.linalg.lstsq(A, Y.reshape(len(tt), -1), rcond=None)[0]
        x0 = preset["reference_time"] - tt.mean()
        mu_trend = np.exp(coef[0] + coef[1] * x0).reshape(Y.shape[1:]) - 1e-6
        fc = rules["components"]["factor"]
        phi_R, _ = gene_factor(np.clip(mu_trend, 0, None), wt.mu, V.p, ok, fc["eps"], tuple(fc["clip"]))
        meta["R"] = {"series": series, "labels_used": int(ok.sum()), "sd_log_phi": float(np.log(phi_R).std())}
        log(f"component R from {list(series)} at {preset['reference_time']}: {int(ok.sum())} labels")

    # targets: as released (scoring) and deduplicated (pools, calls)
    targets = []
    for name, path in paths["targets"]:
        rel = load_stage(path, genes=genes, drop_duplicates=False)
        pool_ad = load_stage(path, genes=genes, drop_duplicates=True)
        lm_rel, col, levels = littermate_split(rel)
        lm_pool, _, _ = littermate_split(pool_ad)
        rec = {**file_record(path, pool_ad), "as_released_cells": int(rel.n_obs), "genotype_column": col,
               "genotype_levels": levels, "littermate_wt_cells": int(lm_pool.sum())}
        hash_rel = row_hashes_layer(rel)
        X_rel = rel.X.toarray().astype(np.float32)
        lm_obj = None
        if lm_pool.any():
            lm_ad = pool_ad[lm_pool].copy()
            pool_ad = pool_ad[~lm_pool].copy()
            X_rel, hash_rel = X_rel[~lm_rel], hash_rel[~lm_rel]
            lm_obj = littermate(V, wt, harm, rules, lm_ad)
            log(f"  {name}: {int(lm_pool.sum())} littermate wild-type cells split off ({col}: {levels})")
        add_label_column(pool_ad, preset["label_mode"])
        calls, key = harm.calls(pool_ad, seed=hp["seed"])
        C = dense_counts(pool_ad)
        lab = calls_to_index(V, calls)
        rec.update(median_depth=median_depth(pool_ad), harmoniser_depth_target=key,
                   rounded_count_entries=int(pool_ad.uns.get("t3_rounded_entries", 0)))
        own = pool_ad.obs["t3_label"].astype(str).to_numpy()
        xw = pd.crosstab(pd.Series(own, name="own_label"), pd.Series(calls["label_raw"].astype(str).to_numpy(),
                                                                       name="wt_label"), normalize="index")
        xw.round(4).to_csv(out / f"{preset_name}_{name}_crosswalk.csv")
        targets.append(Target(name=name, path=str(path), X_rel=X_rel, hash_rel=hash_rel, ad=pool_ad, C=C, lab=lab,
                              calls=calls, key=key, hash_pool=row_hashes_layer(pool_ad), littermate=lm_obj,
                              record=rec))
        meta["targets"][name] = rec
        log(f"target {name}: {rel.n_obs} cells as released, {pool_ad.n_obs} in the pool "
            f"({rec['duplicates_dropped']} duplicates dropped), median depth {rec['median_depth']:.0f} "
            f"(WT {harm.ref_depth:.0f})")
        del rel

    # WT out-of-fold call shares per harmoniser depth target
    p_oof = {}
    for T in targets:
        if T.key not in p_oof:
            p_oof[T.key] = bracket(harm.oof(T.key), harm.model(T.key))

    ctx = Context(preset_name=preset_name, preset=preset, rules=rules, V=V, wt=wt, ref_X=ref_X, harm=harm,
                  p_oof=p_oof, targets=targets, stages=stages, comps=comps, blind=blind, phi_R=phi_R, transfer={},
                  threads=threads, genes=genes)
    # full-pool estimates: the analytic decomposition and the source of every borrowed component
    for T in targets:
        T.full = estimate(ctx, T, np.arange(len(T.lab)), stable_rng(0, f"full:{T.name}"))
    return ctx, meta


@dataclass
class Littermate:
    pool: Pool
    q_hard: np.ndarray
    phi: np.ndarray


def littermate(V, wt, harm, rules, a) -> Littermate:
    """Littermate wild-type cells inside a knockout file: their composition and their per-gene factor against the
    reference (a batch-matched reference, if the files carry one)."""
    hp, fc, sm = rules["harmonisation"], rules["components"]["factor"], rules["components"]["smoothing"]
    calls, _ = harm.calls(a, seed=hp["seed"])
    lab = calls_to_index(V, calls)
    C = dense_counts(a)
    count = np.bincount(lab[lab >= 0], minlength=V.L)
    ok = (count >= 100) & (wt.count >= 100)
    phi, _ = gene_factor(props_type_means(C, lab, V.L), wt.mu, V.p, ok, fc["eps"], tuple(fc["clip"]))
    return Littermate(Pool(C, lab, smoothed_props(C, sm["n_pcs"], sm["k"]), V.L), count / max(count.sum(), 1), phi)


# ---------------------------------------------------------------------------------------------------------------------
# Composition calibration (D7 / D7m), constructive decomposition, transfer factors


def calibrate(V: Vocab, harm: Harmoniser, wlab, targets, cal: dict) -> tuple[dict, pd.DataFrame]:
    """Within-lineage log2 shifts of the knockout(s) against the reference, per estimator and z-trim; returns the
    magnitude summary and the per-state table (see rules.yaml, calibration). targets: (name, AnnData, calls, depth
    key) per knockout embryo; wlab: label index of every reference pool cell."""
    from vec.harmonise import default_trims, trim_mask, z_levels

    states = np.array([i for i in V.model_idx if V.mesp1[i] == "yes" and V.count[i] >= cal["min_wt_cells"]])
    if states.size == 0:
        return {"m": 0.0, "m_raw": 0.0, "n_counted": 0, "states": []}, pd.DataFrame()
    trims = default_trims(k=cal["trims"]["k"], both=tuple(cal["trims"]["both"]))
    ps = cal["pseudo"]
    wt_zl = z_levels(harm.ref)
    names = [t[0] for t in targets]
    ell = {}          # (target, est, trim) -> log2 shift over states
    frac = {}         # (target, est) -> untrimmed within-lineage knockout fractions
    skipped = set()
    for name, ad_k, calls_k, key in targets:
        model = harm.model(key)
        oof = harm.oof(key)
        opos = harm.oof_pos[key]
        ko_zl = z_levels(ad_k)
        for tname, (al, ah, bl, bh) in trims.items():
            wmask = trim_mask(wt_zl, al, ah)
            kmask = trim_mask(ko_zl, bl, bh)
            if (np.unique(np.asarray(wt_zl.index)[wmask]).size < 2
                    or np.unique(np.asarray(ko_zl.index)[kmask]).size < 2):
                skipped.add(tname)      # as vec.harmonise.trim_sensitivity: fewer than two levels left
                continue
            cnt = np.bincount(wlab[wmask], minlength=V.L).astype(np.float64)
            p_tr = cnt / cnt.sum()
            br_w = bracket(oof, model, np.flatnonzero(wmask[opos]))
            br_k = bracket(calls_k, model, np.flatnonzero(kmask))
            for e in EST:
                f = bias_matched(V, br_k[e], br_w[e], p=p_tr)
                x, y = within(f, states), within(p_tr, states)
                ell[(name, e, tname)] = np.log2((x + ps) / (y + ps))
                if tname == "none":
                    frac[(name, e)] = x
    keys = [k for k in ell if k[2] not in skipped]
    L = np.stack([ell[k] for k in keys])                       # combos x states
    sign_ok = np.all(np.sign(L) == np.sign(L[:1]), axis=0) & (np.abs(L).min(0) > 0)
    un = np.stack([ell[(n, e, "none")] for n in names for e in EST])
    noise = np.zeros(states.size)
    if len(targets) == 2:
        a, b = names
        noise = np.max([np.abs(np.log2((frac[(a, e)] + ps) / (frac[(b, e)] + ps))) for e in EST], axis=0)
        counted = sign_ok & (np.abs(un).min(0) > noise)
    else:
        counted = sign_ok
    mean_abs = np.abs(un).mean(0)
    table = pd.DataFrame({"label": [V.labels[i] for i in states], "group": V.group[states],
                          "wt_cells": V.count[states], "wt_within_share": within(V.p, states),
                          **{f"l2_{n}_{e}": ell[(n, e, "none")] for n in names for e in EST},
                          "l2_min_over_trims": L.min(0), "l2_max_over_trims": L.max(0),
                          "sign_consistent": sign_ok, "replicate_noise_bar": noise,
                          "mean_abs_l2": mean_abs, "counted": counted})
    def med(v):
        return float(np.median(v)) if len(v) else 0.0

    summary = {"m": med(mean_abs[counted]), "m_lo": med(np.abs(un).min(0)[counted]),
               "m_hi": med(np.abs(un).max(0)[counted]), "m_raw": med(mean_abs[sign_ok]),
               "n_states": int(states.size), "n_counted": int(counted.sum()), "n_sign_consistent": int(sign_ok.sum()),
               "counted_states": [V.labels[i] for i in states[counted]],
               "directions": {V.labels[i]: ("up" if un[:, j].mean() > 0 else "down")
                              for j, i in enumerate(states) if counted[j]},
               "n_trims": len(trims) - len(skipped), "skipped_trims": sorted(skipped), "replicates": names}
    return summary, table


def constructive(ctx, T) -> tuple[pd.DataFrame, dict, np.ndarray]:
    """Full-file decomposition of the pseudobulk delta into WT-frame components (coverage, Mesp1 composition, stage,
    depth, factor) and a residual; shares, DE genes carried, a joint non-negative fit and within-type sharing."""
    from common.core_metrics import de_genes

    V, wt, E = ctx.V, ctx.wt, T.full
    rng = stable_rng(0, f"constructive:{T.name}")
    Xk = lognorm(T.C).astype(np.float64)
    dt = Xk.mean(0) - wt.X.mean(0)
    up, dn, _ = de_genes(lognorm(T.C), lognorm(wt.pool.C))
    q = E.q_hard
    comps = {"coverage": (E.f["top"]["cov"] - V.p) @ wt.M,
             "mesp1_composition": (E.f["top"]["comp"] - E.f["top"]["cov"]) @ wt.M,
             "stage": (E.alphas @ stage_axes(ctx, q)) if ctx.stages else np.zeros_like(dt)}
    Cdk = rescale_depth(wt.pool.C, E.ratio[wt.pool.lab], wt.pool.P, rng)
    comps["depth"] = q @ (type_means(lognorm(Cdk), wt.pool.lab, V.L) - wt.M)
    Cgf = apply_gene_factor(wt.pool.C, E.phi_nm, wt.pool.P, rng)
    comps["factor"] = q @ (type_means(lognorm(Cgf), wt.pool.lab, V.L) - wt.M)
    comps["residual"] = dt - sum(comps.values())
    idx = np.concatenate([up, dn]).astype(int)
    rows = []
    for name, c in comps.items():
        r = {"component": name, "share_of_dt": float(c @ dt / max(dt @ dt, 1e-12)),
             "spearman_dt": float(spearmanr(c, dt).correlation) if np.std(c) > 0 else np.nan}
        if len(idx):
            a, b = c[idx], dt[idx]
            r.update(sign_agree_DE=float(np.mean(np.sign(a) == np.sign(b))),
                     n_DE_carried=int(np.sum((np.sign(a) == np.sign(b)) & (np.abs(a) >= 0.5 * np.abs(b)))))
        rows.append(r)
    names = [k for k in comps if k != "residual"]
    A = np.vstack([comps[k] for k in names]).T
    tot = float(((dt - dt.mean()) ** 2).sum())
    w, _ = nnls(A, dt)
    r2 = 1 - float(((dt - A @ w) ** 2).sum()) / tot
    drop = {}
    for i, k in enumerate(names):
        keep = [j for j in range(len(names)) if j != i]
        wi, _ = nnls(A[:, keep], dt)
        drop[k] = r2 - (1 - float(((dt - A[:, keep] @ wi) ** 2).sum()) / tot)
    # is the within-type change shared across labels?
    big = np.flatnonzero((E.count >= 300) & (wt.count >= 300))
    Mk = type_means(Xk, T.lab, V.L)
    W = (Mk - wt.M)[big]
    R = np.corrcoef(W) if len(big) > 1 else np.ones((1, 1))
    fc = ctx.rules["components"]["factor"]
    ok = (E.count >= fc["min_cells"]) & (wt.count >= fc["min_cells"])
    phi_m1, den_m1 = gene_factor(E.mu, wt.mu, V.p, ok & (V.mesp1 == "yes"), fc["eps"], tuple(fc["clip"]))
    _, den_nm = gene_factor(E.mu, wt.mu, V.p, ok & V.nm, fc["eps"], tuple(fc["clip"]))
    expr = (den_nm > 1e-4 * max(den_nm.sum(), 1e-12)) & (den_m1 > 1e-4 * max(den_m1.sum(), 1e-12))
    info = {"n_DE_up": int(len(up)), "n_DE_down": int(len(dn)),
            "joint_nnls_weights": dict(zip(names, map(float, w), strict=True)), "joint_r2": r2,
            "drop_one_r2_loss": drop,
            "residual_share": float(comps["residual"] @ dt / max(dt @ dt, 1e-12)),
            "within_type_mean_pairwise_r": float(R[~np.eye(len(big), dtype=bool)].mean()) if len(big) > 1 else None,
            "log_factor_nonmesp1_vs_mesp1_spearman": float(spearmanr(np.log(E.phi_nm[expr]),
                                                                     np.log(phi_m1[expr])).correlation)
            if expr.sum() > 2 else None,
            "dt_spearman_with_wt_pseudobulk": float(spearmanr(dt, wt.X.mean(0)).correlation),
            "full_pool_alphas": E.alphas.tolist(), "stage_pools": [s.name for s in ctx.stages]}
    return pd.DataFrame(rows), info, comps["residual"]


def diag_projection(ctx, T, name, path) -> dict:
    """A WT sample used only as a diagnostic axis: projection of the knockout's depth-free within-type change on that
    sample's within-type offset from the reference (never enters a prediction)."""
    from vec.prep import load_stage

    hp, V, wt, E = ctx.rules["harmonisation"], ctx.V, ctx.wt, T.full
    a = load_stage(path, genes=ctx.genes, drop_duplicates=True)
    calls, _ = ctx.harm.calls(a, seed=hp["seed"])
    lab = calls_to_index(V, calls, hp["confident_call"])
    C = dense_counts(a)
    keep = lab >= 0
    mn = ctx.rules["components"]["stage"]["min_cells"]
    cnt = np.bincount(lab[keep], minlength=V.L)
    med = label_medians(lab, C.sum(1), V.L)
    ok = (cnt >= 20) & (wt.count >= 20) & np.isfinite(med) & (med > 0)
    r = np.where(ok, wt.med / np.where(med > 0, med, 1.0), 1.0)
    sm = ctx.rules["components"]["smoothing"]
    P = smoothed_props(C[keep], sm["n_pcs"], sm["k"])
    Cdm = rescale_depth(C[keep], r[lab[keep]], P, stable_rng(0, f"diag:{name}"))
    Mo = type_means(lognorm(Cdm), lab[keep], V.L)
    S = (((E.q_hard * (cnt >= mn))[:, None]) * (Mo - wt.M)).sum(0)
    Ck = rescale_depth(E.pool.C, np.where(E.ratio > 0, 1.0 / E.ratio, 1.0)[np.maximum(E.pool.lab, 0)], E.pool.P,
                       stable_rng(1, f"diag:{name}"))
    wtc = E.q_hard @ (type_means(lognorm(Ck), E.pool.lab, V.L) - wt.M)
    return {"file": str(path), "ols_coefficient": float(S @ wtc / max(S @ S, 1e-12)),
            "cosine": float(S @ wtc / max(np.linalg.norm(S) * np.linalg.norm(wtc), 1e-12)),
            "labels_used": int((cnt >= mn).sum())}


def transfer_factors(ctx, path, meta) -> dict:
    """The released Mab21l2 KO's per-gene factor against the reference (all shared labels; non-Mesp1 labels) and
    against the E9.5 stage pool (its knockout-specific efficiency), for transfer variants on reference cells."""
    from vec.prep import load_stage

    hp, V, wt, fc = ctx.rules["harmonisation"], ctx.V, ctx.wt, ctx.rules["components"]["factor"]
    a = load_stage(path, genes=ctx.genes, drop_duplicates=True)
    calls, key = ctx.harm.calls(a, seed=hp["seed"])
    lab = calls_to_index(V, calls)
    C = dense_counts(a)
    mu = props_type_means(C, lab, V.L)
    cnt = np.bincount(lab[lab >= 0], minlength=V.L)
    both = (cnt >= 100) & (wt.count >= 100)
    out = {"mab21l2_all": gene_factor(mu, wt.mu, V.p, both, fc["eps"], tuple(fc["clip"]))[0],
           "mab21l2_nm": gene_factor(mu, wt.mu, V.p, both & V.nm, fc["eps"], tuple(fc["clip"]))[0]}
    s95 = [s for s in ctx.stages if s.name == "E9.5"]
    if s95:
        sp = s95[0]
        mu95 = props_type_means(sp.pool.C, sp.pool.lab, V.L)
        both95 = (cnt >= 100) & (sp.count >= 100)
        out["mab21l2_vs95"] = gene_factor(mu, mu95, V.p, both95, fc["eps"], tuple(fc["clip"]))[0]
    meta["transfer_ko"] = {**file_record(path, a), "labels_used_all": int(both.sum()), "harmoniser_depth_target": key}
    log(f"transfer knockout {path}: {a.n_obs} cells, factor from {int(both.sum())} shared labels")
    return out


def mab21l2_calibration(rules, paths, n_jobs) -> tuple[dict, pd.DataFrame]:
    """D7m for the gata4 run: Mab21l2 vs WT E9.5 in the E9.5 vocabulary (the mab21l2 preset's own calibration)."""
    from vec.prep import load_stage

    preset, hp = rules["presets"]["mab21l2"], rules["harmonisation"]
    wt_ad = load_stage(paths["stage:E9.5"], drop_duplicates=True)
    add_label_column(wt_ad, preset["label_mode"])
    harm = Harmoniser(wt_ad, "t3_label", label_groups(wt_ad, "t3_label"), hp, n_jobs)
    V = make_vocab(wt_ad, "t3_label", harm.model(None), rules)
    pos = {t: i for i, t in enumerate(V.labels)}
    wlab = np.array([pos[t] for t in wt_ad.obs["t3_label"].astype(str)])
    ko = load_stage(paths["transfer_ko"], genes=list(wt_ad.var_names), drop_duplicates=True)
    calls, key = harm.calls(ko, seed=hp["seed"])
    return calibrate(V, harm, wlab, [("mab21l2", ko, calls, key)], rules["calibration"])


# ---------------------------------------------------------------------------------------------------------------------
# Summaries and rules


def summarise(rec: pd.DataFrame) -> pd.DataFrame:
    df = rec[rec["variant"] != "_meta"].copy()
    num = [c for c in df.columns if c not in ("target", "seed", "variant") and pd.api.types.is_numeric_dtype(df[c])]
    g = df.groupby(["target", "variant"])
    out = g[num].mean().join(g[num].std(ddof=1), rsuffix="_sd")
    out["sev_floor_frac"] = g["sev_floor"].mean()
    out["n_seeds"] = g["seed"].nunique()
    out = out.reset_index()
    out["prototype_name"] = out["variant"].map(PROTOTYPE_NAMES)
    return out


def shapley_table(rec: pd.DataFrame, comps, metrics=("t3", "de_direction", "de_score", "skill_severity_slope",
                                                      "skill_mmd_u", "skill_variogram")) -> pd.DataFrame:
    """Per target, estimator and component: mean (sd over seeds) of its effect added alone, its Shapley value over
    the full factorial and its effect added last."""
    df = rec[rec["variant"].str.startswith("fac")]
    n = len(comps)
    rows = []
    for (tgt, seed), g in df.groupby(["target", "seed"]):
        val = g.set_index("variant")
        for e in EST:
            def v(on, m, e=e, val=val):
                return float(val.at[fac_name(e, on), m])
            for k in comps:
                others = [c for c in comps if c != k]
                for m in metrics:
                    sh = 0.0
                    for r in range(n):
                        for S in _subsets(others, r):
                            wgt = factorial(len(S)) * factorial(n - len(S) - 1) / factorial(n)
                            sh += wgt * (v(S + (k,), m) - v(S, m))
                    rows.append({"target": tgt, "seed": seed, "estimator": e, "component": k, "metric": m,
                                 "alone": v((k,), m) - v((), m), "shapley": sh,
                                 "last": v(tuple(comps), m) - v(tuple(others), m)})
    t = pd.DataFrame(rows)
    g = t.groupby(["target", "estimator", "component", "metric"])[["alone", "shapley", "last"]]
    return g.mean().join(g.std(ddof=1), rsuffix="_sd").reset_index()


def _subsets(items, r):
    from itertools import combinations

    return [tuple(c) for c in combinations(items, r)]


def evaluate_rules(rules, rec, summ, shap, cal, cal_m, analytic, ctx, extra) -> dict:
    """Apply rules.yaml to the scores. Every outcome is True, False or None (not applicable to this run)."""
    R = rules["rules"]
    tags = [T.name for T in ctx.targets]
    S = summ.set_index(["target", "variant"])
    have = set(summ["variant"])
    prim = ctx.preset["primary_blind"]
    out: dict = {"targets": tags}

    def m(t, v, c):
        try:
            x = S.at[(t, v), c]
            return float(x) if pd.notna(x) else None
        except KeyError:
            return None

    def gain(t, a, b):
        piv = rec[(rec["target"] == t) & rec["variant"].isin([a, b])].pivot_table(index="seed", columns="variant",
                                                                                 values="t3")
        if a not in piv or b not in piv:
            return None, None
        d = (piv[a] - piv[b]).dropna()
        return float(d.mean()), (float(d.std(ddof=1)) if len(d) > 1 else 0.0)

    def leq(x, thr):
        return x is not None and x <= thr

    # D1
    r1 = R["D1"]

    def sign_ok(t, v):
        d, b, r2, fl = (m(t, v, c) for c in ("de_direction", "beta", "r2", "sev_floor_frac"))
        if None in (d, b, r2, fl):
            return None
        return bool(d > r1["de_direction_gt"] and b > r1["beta_gt"] and r2 > r1["r2_gt"]
                    and fl <= r1["sev_floor_frac_le"])
    v1 = r1["variant"].format(primary=prim)
    per = {t: sign_ok(t, v1) for t in tags}
    out["D1"] = {"variant": v1, "per_target": per,
                 "outcome": None if None in per.values() else all(per.values()),
                 "D1b_all_blind_rules": {f"blind:{b}": {t: sign_ok(t, f"blind:{b}") for t in tags}
                                         for b in ctx.preset["blind"]}}
    # D2 and D2u
    two = len(tags) == 2
    if "xrep:F" in have:
        r = R["D2"]
        per = {t: {"t3": m(t, r["variant"], "t3"), "sev_floor_frac": m(t, r["variant"], "sev_floor_frac")}
               for t in tags}
        rho = extra.get("nonmesp1_delta_spearman")
        ok = all(p["t3"] is not None and p["t3"] >= r["t3_ge"] and leq(p["sev_floor_frac"], r["sev_floor_frac_le"])
                 for p in per.values())
        out["D2"] = {"per_target": per, "nonmesp1_delta_spearman": rho,
                     "also_reported_xrep_N": {t: m(t, "xrep:N", "t3") for t in tags},
                     "outcome": bool(ok and rho is not None and rho >= r["nonmesp1_delta_spearman_ge"])}
        r = R["D2u"]
        per = {}
        for t in tags:
            g_, sd = gain(t, r["variant"], r["baseline"])
            per[t] = {"t3_gain": g_, "t3_gain_sd": sd, "sev_floor_frac": m(t, r["variant"], "sev_floor_frac")}
        out["D2u"] = {"per_target": per, "outcome": all(p["t3_gain"] is not None and p["t3_gain"] >= r["t3_gain_ge"]
                                                         and leq(p["sev_floor_frac"], r["sev_floor_frac_le"])
                                                         for p in per.values())}
    else:
        out["D2"] = out["D2u"] = {"outcome": None, "note": "needs two replicates"}
    # D3
    r = R["D3"]
    d3 = {}
    sh = shap.set_index(["target", "estimator", "component", "metric"])
    full = fac_name("top", ctx.comps)
    for k in ctx.comps:
        per = {}
        carries = True
        for t in tags:
            for e in EST:
                dd, tt = sh.at[(t, e, k, "de_direction"), "shapley"], sh.at[(t, e, k, "t3"), "shapley"]
                fl_full = m(t, fac_name(e, ctx.comps), "sev_floor_frac")
                fl_alone = m(t, fac_name(e, (k,)), "sev_floor_frac")
                per[f"{t}|{e}"] = {"shapley_de_direction": float(dd), "shapley_t3": float(tt),
                                   "alone_t3": float(sh.at[(t, e, k, "t3"), "alone"]),
                                   "last_t3": float(sh.at[(t, e, k, "t3"), "last"]),
                                   "alone_sev_floor_frac": fl_alone, "full_set_sev_floor_frac": fl_full}
                carries &= bool(dd >= r["shapley_de_direction_ge"] and tt >= r["shapley_t3_ge"]
                                and leq(fl_alone, r["sev_floor_frac_le"]) and leq(fl_full, r["sev_floor_frac_le"]))
        d3[k] = {"per_target_estimator": per, "carries": carries}
    out["D3"] = {"components": d3, "full_set_variant": full, "outcome": {k: v["carries"] for k, v in d3.items()}}

    # D4-D6 (transfer of each extra component on top of the borrowed factor)
    def extra_rule(rid, oracle_check=None):
        r = R[rid]
        if r["variant"] not in have:
            return {"outcome": None, "note": "needs two replicates"}
        per, ok = {}, True
        for t in tags:
            g_, sd = gain(t, r["variant"], r["baseline"])
            fl = m(t, r["variant"], "sev_floor_frac")
            per[t] = {"t3_gain": g_, "t3_gain_sd": sd, "sev_floor_frac": fl}
            ok &= bool(g_ is not None and g_ >= r["t3_gain_ge"] and leq(fl, r["sev_floor_frac_le"]))
        res = {"per_target": per}
        if oracle_check is not None:
            chk, cok = oracle_check(r)
            res["oracle_check"] = chk
            ok &= cok
        res["outcome"] = bool(ok)
        return res

    def depth_oracle(r):
        chk = {f"{t}|{e}": m(t, fac_name(e, ctx.comps), "sev_floor_frac") for t in tags for e in EST}
        return chk, all(leq(x, r["oracle_full_set_sev_floor_frac_le"]) for x in chk.values())

    def stage_oracle(r):
        chk = {t: gain(t, "fac:stage+factor", "fac:factor")[0] for t in tags}
        return chk, all(x is not None and x >= r["oracle_gain_after_factor_ge"] for x in chk.values())

    out["D4"] = extra_rule("D4", depth_oracle)
    out["D5"] = extra_rule("D5", stage_oracle) if ctx.stages else {"outcome": False, "note": "no stage pools"}
    out["D6"] = extra_rule("D6")
    # depth floor check on the oracle factorial alone (prototype D5), reported for every run
    flags = []
    for t in tags:
        for v in sorted(have):
            if not v.startswith("fac") or "depth" not in fac_parse(v)[1]:
                continue
            e, on = fac_parse(v)
            base = fac_name(e or "top", [c for c in on if c != "depth"])
            fa, fb = m(t, base, "sev_floor_frac"), m(t, v, "sev_floor_frac")
            if fa is not None and fb is not None and fb >= 0.5 > fa:
                flags.append(f"{t}:{v}")
    out["depth_floor_flags"] = sorted(flags)
    # D7, D7m
    out["D7"] = ({**cal, "outcome": cal["n_counted"] >= R["D7"]["min_counted_states"]} if cal is not None
                 else {"outcome": None, "note": "needs two replicates"})
    out["D7m"] = ({**cal_m, "outcome": cal_m["n_counted"] >= R["D7m"]["min_counted_states"]} if cal_m is not None
                  else {"outcome": None})
    out["m_raw_positive"] = {"outcome": bool(cal is not None and cal.get("m_raw", 0) > 0)}
    # D8
    r = R["D8"]
    if r["variant"] in have:
        per, fires = {}, False
        for t in tags:
            g_, sd = gain(t, r["variant"], r["baseline"])
            fl = m(t, r["variant"], "sev_floor_frac")
            per[t] = {"t3_change": g_, "t3_change_sd": sd, "sev_floor_frac": fl}
            fires |= bool((g_ is not None and g_ <= r["fires_if_t3_change_le"])
                          or (fl is not None and fl > r["fires_if_sev_floor_frac_gt"]))
        out["D8"] = {"per_target": per, "scenario_log2_magnitude": extra.get("scenario_log2"),
                     "phenotype_rule_at_m": {t: gain(t, f"xrep:F+{prim}@m", "xrep:F")[0] for t in tags},
                     "outcome": bool(fires)}
    else:
        out["D8"] = {"outcome": None, "note": "needs two replicates"}
    # D9
    r = R["D9"]
    shares = {t: analytic[t]["residual_share"] for t in tags}
    rho = extra.get("residual_spearman")
    ok = all(s is not None and s >= r["residual_share_ge"] for s in shares.values())
    out["D9"] = {"residual_share": shares, "replicate_residual_spearman": rho,
                 "outcome": bool(ok and (rho >= r["replicate_residual_spearman_ge"] if two else True))}
    # D10
    r = R["D10"]
    if r["variant"] in have:
        per = {t: {"t3": m(t, r["variant"], "t3"), "sev_floor_frac": m(t, r["variant"], "sev_floor_frac")}
               for t in tags}
        out["D10"] = {"per_target": per, "series": extra.get("R_series"),
                      "outcome": all(p["t3"] is not None and p["t3"] >= r["t3_ge"]
                                     and leq(p["sev_floor_frac"], r["sev_floor_frac_le"]) for p in per.values())}
    else:
        out["D10"] = {"outcome": None, "note": "no WT series"}
    # D11
    r = R["D11a"]
    rho = extra.get("log_factor_spearman")
    out["D11a"] = {"replicate_log_factor_spearman": rho,
                   "outcome": None if rho is None else bool(rho >= r["replicate_log_factor_spearman_ge"])}
    r = R["D11b"]
    if r["variant"] in have:
        per = {t: {"t3": m(t, r["variant"], "t3"), "sev_floor_frac": m(t, r["variant"], "sev_floor_frac")}
               for t in tags}
        out["D11b"] = {"per_target": per, "outcome": all(p["t3"] is not None and p["t3"] >= r["t3_ge"]
                                                          and leq(p["sev_floor_frac"], r["sev_floor_frac_le"])
                                                          for p in per.values())}
    else:
        out["D11b"] = {"outcome": None}
    # W1
    W = rules["weights"]["factor_exponent"]
    tmpl = W["variant_depth_on"] if out["D4"].get("outcome") else W["variant_depth_off"]
    grid = []
    for w in W["grid"]:
        v = tmpl.format(w=w)
        if v not in have:
            continue
        t3s = [m(t, v, "t3") for t in tags]
        fls = [m(t, v, "sev_floor_frac") for t in tags]
        grid.append({"w": float(w), "variant": v, "mean_t3": float(np.mean(t3s)) if None not in t3s else None,
                     "feasible": all(leq(x, W["constraint_sev_floor_frac_le"]) for x in fls)})
    feas = [g for g in grid if g["feasible"] and g["mean_t3"] is not None]
    if feas:
        best = max(g["mean_t3"] for g in feas)
        near = [g for g in feas if g["mean_t3"] >= best - W["tie_window_t3"]]
        chosen = min(near, key=lambda g: (abs(g["w"] - 1.0), g["w"]))["w"]
        rest = sorted([g for g in feas if g["w"] != chosen], key=lambda g: -g["mean_t3"])
        runner = rest[0]["w"] if rest else None
    else:
        chosen, runner = float(W["default"]), None
    out["W1"] = {"grid": grid, "chosen": chosen if grid else None, "runner_up": runner,
                 "applicable": bool(grid), "variant_template": tmpl}
    return out


def decide_submissions(rules, dec, ctx, extra) -> dict:
    """Both submissions for each organiser answer, from the tables in rules.yaml."""
    sub = rules["submissions"]
    names = {k: v.get("outcome") for k, v in dec.items() if isinstance(v, dict) and "outcome" in v
             and not isinstance(v.get("outcome"), dict)}
    w1 = dec.get("W1", {})
    cal = dec.get("D7", {})
    cal_m = dec.get("D7m", {})
    scen = rules["presets"]["gata4"]["scenario"]
    wts = rules["weights"]
    out = {}
    for branch in ("allowed", "refused"):
        b = sub[branch]
        row = first_match(b["factor_table"], names)
        mag = cal_m if b.get("composition_magnitude") == "D7m" else cal

        def comp_spec(kind, mag=mag, b=b):
            if kind == "none":
                return {"kind": "none"}
            base = mag.get("m_raw" if kind == "mild_raw" else "m", 0.0) or 0.0
            l2 = float(base * wts["composition"][kind])
            return {"kind": kind, "depleted_states": scen["depleted_states"], "spared_states": scen["spared_states"],
                    "log2_magnitude": l2, "depletion_factor": float(2.0 ** -l2),
                    "magnitude_source": "Mab21l2 (D7m)" if mag is cal_m else "Gata4 replicates (D7)"}

        def spec(source, comp_kind, exponent=None, b=b):
            exps = b["factor_exponent"]
            chosen = w1.get("chosen") if w1.get("chosen") is not None else 1.0
            e = exponent if exponent is not None else (
                chosen if exps.get(source) == "W1" else float(exps.get(source, 1.0)))
            s = {"factor": {"source": source, "exponent": e}, "composition": comp_spec(comp_kind)}
            for comp, rid in b.get("gata4_extras", {}).items():
                s[comp] = bool(source == "gata4" and dec.get(rid, {}).get("outcome"))
            return s

        c1 = first_match(b["composition_s1"], names)["value"]
        tab2 = b.get("composition_s2_other_factor") if row["s2"] != row["s1"] else None
        c2 = first_match(tab2 or b["composition_s2_same_factor"], names)["value"]
        s1, s2 = spec(row["s1"], c1), spec(row["s2"], c2)
        note = None
        if s1 == s2:
            alt = w1.get("runner_up") if branch == "allowed" and row["s1"] == "gata4" else None
            alt = alt if alt is not None else 0.5 * (s1["factor"]["exponent"] or 0.0)
            s2 = spec(row["s2"], c2, exponent=alt)
            note = f"s2 equalled s1: s2 takes factor exponent {alt}"
        warn = None
        if s1["factor"]["source"] == "none" and s1["composition"]["kind"] == "none":
            warn = "no validated component: s1 is plain WT resampling; review before any upload"
        out[branch] = {"factor_row": row, "s1": s1, "s2": s2, "note": note, "warning": warn}
    out["parameters"] = extra.get("submission_parameters", {})
    out["rule_outcomes"] = names
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Runs


def run(preset_name, paths, args, rules):
    from threadpoolctl import threadpool_limits

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    seeds = parse_seeds(args.seeds) if args.seeds else list(rules["protocol"]["seeds"])
    meta = {"preset": preset_name, "argv": sys.argv, "started": time.strftime("%Y-%m-%d %H:%M:%S"),
            "seeds": seeds, "protocol_as_registered": seeds == list(rules["protocol"]["seeds"]),
            "workers": args.workers, "threads_per_worker": args.threads, **provenance(paths)}
    if preset_name == "gata4":
        meta["scenario"] = check_scenario(rules)
    log(f"T3 decomposition, preset {preset_name}; rules sha256 {meta['sha256']['rules'][:16]}, script "
        f"{meta['sha256']['script'][:16]}; seeds {seeds}")
    n_par = args.workers * args.threads
    extra: dict = {}
    # The parent never starts an OpenMP thread pool before forking the workers (libgomp is not fork-safe; OpenBLAS
    # is), and uses at most workers x threads BLAS threads on the shared machine.
    with threadpool_limits(limits=n_par, user_api="blas"), threadpool_limits(limits=1, user_api="openmp"):
        ctx, cmeta = build_context(preset_name, paths, rules, n_jobs=n_par, threads=args.threads, out=out)
        meta.update(cmeta)
        selfcheck(ctx, seeds[0])
        meta["selfcheck"] = "passed"
        if preset_name == "gata4":
            ctx.transfer = transfer_factors(ctx, paths["transfer_ko"], meta)
        # composition calibration: D7 with two replicates; the single-embryo form (D7m) otherwise
        kos = [(T.name, T.ad, T.calls, T.key) for T in ctx.targets]
        cal = None
        if len(ctx.targets) == 2:
            cal, cal_tab = calibrate(ctx.V, ctx.harm, ctx.wt.pool.lab, kos, rules["calibration"])
            cal_tab.to_csv(out / f"{preset_name}_calibration.csv", index=False)
            log(f"D7 composition magnitude: m {cal['m']:.3f} ({cal['n_counted']} of {cal['n_states']} states "
                f"counted), m_raw {cal['m_raw']:.3f}")
            cal_m, cal_m_tab = mab21l2_calibration(rules, paths, n_par)
            cal_m_tab.to_csv(out / f"{preset_name}_calibration_mab21l2.csv", index=False)
        else:
            cal_m, cal_m_tab = calibrate(ctx.V, ctx.harm, ctx.wt.pool.lab, kos, rules["calibration"])
            cal_m_tab.to_csv(out / f"{preset_name}_calibration.csv", index=False)
        log(f"D7m (one embryo, sign over every trim and estimator): m {cal_m['m']:.3f} over {cal_m['n_counted']} "
            f"of {cal_m['n_states']} states")
        if preset_name == "gata4":
            sc = rules["presets"]["gata4"]["scenario"]
            l2 = cal["m"] if cal and cal["m"] > 0 else rules["calibration"]["scenario_fallback_log2"]
            ctx.scenario_f, _ = rule_comp(ctx.V, {"deplete": {s_: 0.0 for s_ in sc["depleted_states"]}}, {}, l2)
            lists = {k: v for k, v in ctx.preset.items() if isinstance(v, list)}
            ctx.blind_m, _ = rule_comp(ctx.V, ctx.preset["blind"][ctx.preset["primary_blind"]], lists, l2)
            extra["scenario_log2"] = l2
        extra["R_series"] = meta.get("R")
        # full-file analytic decomposition per target
        analytic, residuals, rows = {}, {}, []
        for T in ctx.targets:
            dec_tab, info, res = constructive(ctx, T)
            dec_tab.insert(0, "target", T.name)
            rows.append(dec_tab)
            analytic[T.name], residuals[T.name] = info, res
        pd.concat(rows).to_csv(out / f"{preset_name}_decomposition.csv", index=False)
        meta["analytic"] = analytic
        for spec in args.diag or []:
            name, path = spec.split("=", 1)
            meta.setdefault("diagnostic_axes", {})[name] = {T.name: diag_projection(ctx, T, name, path)
                                                            for T in ctx.targets}
        if len(ctx.targets) == 2:
            a, b = ctx.targets
            extra["nonmesp1_delta_spearman"] = float(spearmanr(a.full.nm_delta, b.full.nm_delta).correlation)
            extra["residual_spearman"] = float(spearmanr(residuals[a.name], residuals[b.name]).correlation)
            fc = rules["components"]["factor"]
            dens = []
            for T in (a, b):
                ok = (T.full.count >= fc["min_cells"]) & (ctx.wt.count >= fc["min_cells"]) & ctx.V.nm
                dens.append(gene_factor(T.full.mu, ctx.wt.mu, ctx.V.p, ok, fc["eps"], tuple(fc["clip"]))[1])
            expr = (dens[0] > 1e-4 * dens[0].sum()) & (dens[1] > 1e-4 * dens[1].sum())
            extra["log_factor_spearman"] = float(spearmanr(np.log(a.full.phi_nm[expr]),
                                                           np.log(b.full.phi_nm[expr])).correlation)
        composition_table(ctx).to_csv(out / f"{preset_name}_composition.csv", index=False)
        if preset_name == "gata4":
            extra["submission_parameters"] = submission_parameters(ctx, out)
        # scored variants: one task per (target, seed)
        global _CTX
        _CTX = ctx
        tasks = [(ti, s_) for ti in range(len(ctx.targets)) for s_ in seeds]
        log(f"scoring {len(tasks)} (target, seed) tasks on {args.workers} workers x {args.threads} threads")
        recs = []
        if args.workers > 1:
            with mp.get_context("fork").Pool(min(args.workers, len(tasks))) as pool:
                for r in pool.imap_unordered(run_task, tasks):
                    recs.extend(r)
        else:
            for t in tasks:
                recs.extend(run_task(t))
    rec = pd.DataFrame(recs).sort_values(["target", "variant", "seed"])
    rec.to_csv(out / f"{preset_name}_records.csv", index=False)
    summ = summarise(rec)
    summ.to_csv(out / f"{preset_name}_summary.csv", index=False)
    shap = shapley_table(rec, ctx.comps)
    shap.to_csv(out / f"{preset_name}_factorial.csv", index=False)
    dec = evaluate_rules(rules, rec, summ, shap, cal, cal_m, analytic, ctx, extra)
    decisions = {"preset": preset_name, "rules_sha256": meta["sha256"]["rules"],
                 "script_sha256": meta["sha256"]["script"], "protocol_as_registered": meta["protocol_as_registered"],
                 "seeds": seeds, "rules": dec}
    if preset_name == "gata4":
        decisions["submissions"] = decide_submissions(rules, dec, ctx, extra)
    else:
        decisions["submissions"] = "not applicable (smoke test on the training knockout; the gata4 run decides)"
    (out / f"{preset_name}_decisions.json").write_text(json.dumps(json_safe(decisions), indent=1))
    text = decision_text(decisions, summ, shap, ctx)
    (out / f"{preset_name}_decisions.txt").write_text(text)
    meta["runtime_min"] = round((time.time() - t_start) / 60, 1)
    (out / f"{preset_name}_meta.json").write_text(json.dumps(json_safe(meta), indent=1))
    print(text, flush=True)
    log(f"done in {meta['runtime_min']} min; outputs in {out}")
    remote = str(out) if out.is_absolute() else f"v/virtual-embryo-challenge/{out}"
    log("bring back the aggregates (CSV/JSON/TXT; factors/ stays): rsync -az --include='*.csv' --include='*.json' "
        f"--include='*.txt' --exclude='*' PrakashDGX_H2:{remote}/ experiments/t3_decomp/results/{preset_name}/")


def selfcheck(ctx, seed) -> None:
    from vec.boards import Board
    from vec.evaluate import Stage, score_board, subsample

    T = ctx.targets[0]
    frac = ctx.rules["protocol"]["working_copy_fraction"]
    rows = working_copy_rows(len(T.X_rel), frac, seed)
    assert np.array_equal(subsample(Stage(X=T.X_rel), frac, seed).X, T.X_rel[rows]), "working-copy rows differ"
    est = T.full
    X = build(ctx, 3000, stable_rng(seed, "selfcheck"), f=ctx.V.p, phi=est.phi_nm)
    board = Board(key="T3:emulated", task="T3", genes=ctx.genes, needs_coords=False)
    ref = score_board(board, Stage(X=X), Stage(X=T.X_rel), Stage(X=ctx.ref_X), seed=seed)
    # the same scoring as score_seed
    from vec.boards import RANKED_METRICS
    from vec.evaluate import aggregate, panel_metrics, skills_from_anchors, split_half

    truth, other = split_half(Stage(X=T.X_rel[rows]), seed=ctx.rules["protocol"]["split_seed"])
    rw = Stage(X=ctx.ref_X[working_copy_rows(len(ctx.ref_X), frac, seed + 1)])
    fv, cv = panel_metrics("T3", rw, truth, rw, seed=seed), panel_metrics("T3", other, truth, rw, seed=seed)
    anchors = {m: {"floor": fv.get(m), "ceiling": cv.get(m)} for m in RANKED_METRICS["T3"]
               if fv.get(m) is not None and cv.get(m) is not None}
    score, _ = aggregate("T3", skills_from_anchors("T3", panel_metrics("T3", Stage(X=X), truth, rw, seed=seed),
                                                   anchors))
    if abs(score - ref.score) > 1e-9:
        raise RuntimeError(f"self-check failed: factored scoring {score} vs vec.evaluate.score_board {ref.score}")
    log(f"self-check passed: factored scoring equals vec.evaluate.score_board ({score:.4f})")


def composition_table(ctx) -> pd.DataFrame:
    V = ctx.V
    df = pd.DataFrame({"label": V.labels, "group": V.group, "mesp1": V.mesp1, "nonmesp1": V.nm, "wt_cells": V.count,
                       "wt_share": V.p, "wt_median_depth": ctx.wt.med})
    for key, br in ctx.p_oof.items():
        for e in EST:
            col = np.full(V.L, np.nan)
            col[V.model_idx] = br[e]
            df[f"wt_oof_{e}" + ("" if key is None else f"@{key}")] = col
    for T in ctx.targets:
        E = T.full
        for e in EST:
            col = np.full(V.L, np.nan)
            col[V.model_idx] = E.q[e]
            df[f"{T.name}_calls_{e}"] = col
            df[f"{T.name}_f_{e}"] = E.f[e]["comp"]
        df[f"{T.name}_hard_share"] = E.q_hard
        df[f"{T.name}_depth_ratio"] = E.ratio
    for s in ctx.stages:
        df[f"stage_{s.name}_confident_cells"] = s.count
    return df


def submission_parameters(ctx, out) -> dict:
    """What a beta-catenin builder needs from this run: the pooled Gata4 factor and depth ratios, the stage alphas and
    coverage shares averaged over replicates, R and the Mab21l2 factor. Per-gene factors go to factors/*.csv (the
    output of committed code on released data; keep them with the run)."""
    V, wt, fc = ctx.V, ctx.wt, ctx.rules["components"]["factor"]
    a, b = ctx.targets
    C = np.vstack([a.C, b.C])
    lab = np.concatenate([a.lab, b.lab])
    cnt = np.bincount(lab[lab >= 0], minlength=V.L)
    ok = (cnt >= fc["min_cells"]) & (wt.count >= fc["min_cells"])
    phi, _ = gene_factor(props_type_means(C, lab, V.L), wt.mu, V.p, ok & V.nm, fc["eps"], tuple(fc["clip"]))
    med = label_medians(lab, C.sum(1), V.L)
    dmin = ctx.rules["components"]["depth"]["min_cells"]
    okd = (cnt >= dmin) & (wt.count >= dmin) & np.isfinite(med) & (wt.med > 0)
    ratio = np.where(okd, med / np.where(wt.med > 0, wt.med, 1.0), 1.0)
    fdir = out / "factors"
    fdir.mkdir(exist_ok=True)
    factors = {"gata4_pooled": phi, "reference_R": ctx.phi_R, **{f"transfer_{k}": v for k, v in ctx.transfer.items()}}
    files = {}
    for k, v in factors.items():
        if v is None:
            continue
        p = fdir / f"{k}.csv"
        pd.DataFrame({"gene": ctx.genes, "factor": v}).to_csv(p, index=False)
        files[k] = {"file": str(p), "sha256": sha256(p)}
    cov = 0.5 * (a.full.f["top"]["cov"] + b.full.f["top"]["cov"])
    return {"factor_files": files,
            "depth_ratio_pooled": {V.labels[i]: float(ratio[i]) for i in range(V.L) if ratio[i] != 1.0},
            "stage_alphas_mean": dict(zip([s.name for s in ctx.stages],
                                          map(float, 0.5 * (a.full.alphas + b.full.alphas)), strict=True)),
            "coverage_shares_mean": {V.labels[i]: float(cov[i]) for i in np.flatnonzero(V.nm)}}


def decision_text(decisions, summ, shap, ctx) -> str:
    """A readable summary of the decision file."""
    lines = [f"T3 decomposition, preset {decisions['preset']}",
             f"rules sha256 {decisions['rules_sha256']}", f"script sha256 {decisions['script_sha256']}",
             f"protocol as registered: {decisions['protocol_as_registered']}", ""]
    show = ["t3", "t3_sd", "de_score", "de_direction", "skill_severity_slope", "sev_floor_frac", "beta", "r2",
            "sign_agree", "offset_share", "mmd_u", "variogram"]
    pd.set_option("display.width", 250)
    for t in summ["target"].unique():
        s = summ[summ["target"] == t].set_index("variant")
        lines.append(f"== {t}: means over {int(s['n_seeds'].max())} seeds")
        lines.append(s[[c for c in show if c in s]].round(3).to_string())
        lines.append("")
    if len(shap):
        sh = shap[shap["metric"].isin(["t3", "de_direction"])].pivot_table(
            index=["target", "estimator", "component"], columns="metric", values=["alone", "shapley", "last"])
        lines.append("== factorial: component effects (alone / Shapley / added last)")
        lines.append(sh.round(3).to_string())
        lines.append("")
    lines.append("== rule outcomes")
    for k, v in decisions["rules"].items():
        if isinstance(v, dict) and "outcome" in v:
            lines.append(f"{k}: {v['outcome']}")
    if isinstance(decisions.get("submissions"), dict):
        for br in ("allowed", "refused"):
            s = decisions["submissions"][br]
            lines.append(f"\n== submissions if the organisers {br} Gata4 files inside the beta-catenin prediction")
            lines.append(f"s1: {json.dumps(json_safe(s['s1']))}")
            lines.append(f"s2: {json.dumps(json_safe(s['s2']))}")
            if s.get("note"):
                lines.append(f"note: {s['note']}")
            if s.get("warning"):
                lines.append(f"WARNING: {s['warning']}")
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="preset", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--out", required=True, help="output directory (full outputs; aggregates are CSV/JSON/TXT)")
    common.add_argument("--seeds", default=None, help="override the registered seeds (e.g. 0-1 for a code check; "
                        "the decision file then says protocol_as_registered: false)")
    common.add_argument("--workers", type=int, default=8)
    common.add_argument("--threads", type=int, default=2, help="BLAS/OpenMP threads per worker")
    common.add_argument("--diag", action="append", default=None,
                        help="NAME=PATH: a WT sample used only as a diagnostic stage axis (e.g. E10.5=...)")
    m = sub.add_parser("mab21l2", parents=[common], help="Mab21l2 KO vs WT E9.5 (smoke test)")
    m.add_argument("--wt", default=str(DATA / "T2" / "E9.5.h5ad"))
    m.add_argument("--ko", default=str(DATA / "T3" / "E9.5_mab21l2_ko.h5ad"))
    m.add_argument("--e875", default=str(DATA / "T2" / "E8.75.h5ad"))
    m.add_argument("--e825", default=str(DATA / "T2" / "E8.25_late.h5ad"))
    g = sub.add_parser("gata4", parents=[common], help="Gata4 rep1 and rep2 vs WT E8.75 (release day)")
    g.add_argument("rep1")
    g.add_argument("rep2")
    g.add_argument("wt_e875")
    g.add_argument("--e95", default=str(DATA / "T2" / "E9.5.h5ad"))
    g.add_argument("--e825", default=str(DATA / "T2" / "E8.25_late.h5ad"))
    g.add_argument("--e85", default=None, help="WT E8.5 for the R series (default data/raw/T2/E8.5.h5ad if present)")
    g.add_argument("--mab21l2", default=str(DATA / "T3" / "E9.5_mab21l2_ko.h5ad"))
    args = ap.parse_args(argv)
    rules = yaml.safe_load(RULES_PATH.read_text())
    if args.preset == "gata4":
        absent = [p for p in (args.rep1, args.rep2, args.wt_e875) if not Path(p).is_file()]
        if absent:
            print("t3_decompose gata4: input files not found, nothing was run: " + ", ".join(absent) + "\n"
                  "(expected before the 2026-10-20 release; pass the two Gata4 replicate files and WT E8.75)",
                  file=sys.stderr)
            return 0
        e85 = args.e85 or (str(DATA / "T2" / "E8.5.h5ad") if (DATA / "T2" / "E8.5.h5ad").is_file() else None)
        paths = {"wt": args.wt_e875, "targets": [("rep1", args.rep1), ("rep2", args.rep2)],
                 "stage:E9.5": args.e95, "stage:E8.25_late": args.e825, "series:E8.25_late": args.e825,
                 "series:E9.5": args.e95, "series:E8.5": e85, "transfer_ko": args.mab21l2}
    else:
        paths = {"wt": args.wt, "targets": [("mab21l2", args.ko)], "stage:E8.75": args.e875,
                 "series:E8.25_late": args.e825, "series:E8.75": args.e875}
    need = [p for k, p in paths.items() if k != "targets" and p is not None and not Path(p).is_file()]
    need += [p for _, p in paths["targets"] if not Path(p).is_file()]
    if need:
        print("t3_decompose: required released files not found: " + ", ".join(need), file=sys.stderr)
        return 2
    for spec in args.diag or []:
        if "=" not in spec or not Path(spec.split("=", 1)[1]).is_file():
            print(f"t3_decompose: --diag {spec}: expected NAME=PATH to an existing file", file=sys.stderr)
            return 2
    run(args.preset, paths, args, rules)
    return 0


if __name__ == "__main__":
    sys.exit(main())
