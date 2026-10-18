"""Library for the pre-registered T1 external-direction readout (experiments/t1_direction).

Sufficient statistics (per group: cells, per-gene sums of log1p-CP10k, CP10k and detections), gene maps, change
vectors, vector versions of veckit's de_score and de_direction, held-flat genes, label harmonisation, population
predictions and the decision rule. Ported from the screened research bundle (common.py, o1_oct20.py; 2026-10-19).

Named dirlib, not common: veckit installs a top-level package called `common`, which this module imports.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from common.core_metrics import _signed_overlap  # pinned veckit (git main cc3859a): the scorer's own tie handling
from scipy.stats import rankdata, spearmanr

# ------------------------------------------------------------------ sufficient statistics


STAT_FIELDS = ("genes", "keys", "n", "S_log", "S_cp", "S_nz", "depth_sum")


def load_stats(path: str | Path, fields: Sequence[str] | None = None) -> dict:
    """Load a statistics npz (T1 from t1_stats.py; Feng and JAX from data/external, written by the filters).
    `fields` limits what is decompressed (the JAX files also hold raw counts, which are not needed)."""
    z = np.load(path, allow_pickle=True)  # the JAX statistics from the 2026-10-18 run hold object arrays
    d = {k: z[k] for k in (z.files if fields is None else [f for f in fields if f in z.files])}
    d["keys"] = np.array([str(k) for k in d["keys"]], dtype=object)
    d["genes"] = np.array([str(g) for g in d["genes"]], dtype=object)
    return d


def agg(st: dict, pred: Callable[[list[str]], bool]) -> dict | None:
    """Mean log1p-CP10k, mean CP10k, detection rate and mean depth over the rows whose key parts satisfy `pred`."""
    m = np.array([bool(pred(k.split("|"))) for k in st["keys"]])
    n = float(st["n"][m].sum()) if m.any() else 0.0
    if n == 0:
        return None
    return {"n": n, "log": st["S_log"][m].sum(0) / n, "cp": st["S_cp"][m].sum(0) / n,
            "det": st["S_nz"][m].sum(0) / n, "depth": float(st["depth_sum"][m].sum() / n)}


def t1_sel(st: dict, labels: Iterable[str] | None, halves: Sequence[str] = ("A", "B")) -> dict | None:
    """T1 keys are 'celltype|library|half'; labels=None selects every cell type."""
    ls = None if labels is None else set(labels)
    return agg(st, lambda p: (ls is None or p[0] in ls) and p[2] in halves)


def feng_sel(st: dict, stage: str, labels: Iterable[str], halves: Sequence[str] = ("A", "B")) -> dict | None:
    """Feng keys are 'stage|zone|cell_type_spec|half'."""
    ls = set(labels)
    return agg(st, lambda p: p[0] == stage and p[2] in ls and p[3] in halves)


def load_jax(stats_dir: str | Path, embryos: Iterable[int]) -> dict:
    """Stack per-embryo statistics into one table with keys 'embryo|label|half' (gene order checked)."""
    tabs, genes = [], None
    for e in sorted(set(embryos)):
        st = load_stats(Path(stats_dir) / f"embryo_{e}.npz", fields=STAT_FIELDS)
        if genes is None:
            genes = st["genes"]
        elif not np.array_equal(genes, st["genes"]):
            idx = pd.Index(st["genes"]).get_indexer(genes)
            if (idx < 0).any():
                raise ValueError(f"embryo {e}: gene set differs from the first embryo")
            for k in ("S_log", "S_cp", "S_nz"):
                st[k] = st[k][:, idx]
        st["keys"] = np.array([f"{e}|{k}" for k in st["keys"]], dtype=object)
        tabs.append(st)
    out = {"genes": genes, "keys": np.concatenate([t["keys"] for t in tabs]),
           "n": np.concatenate([t["n"] for t in tabs]), "depth_sum": np.concatenate([t["depth_sum"] for t in tabs])}
    for k in ("S_log", "S_cp", "S_nz"):
        out[k] = np.vstack([t[k] for t in tabs])
    return out


def jax_sel(st: dict, embryos: Iterable[int], labels: Iterable[str] | None,
            halves: Sequence[str] = ("A", "B")) -> dict | None:
    es = {str(e) for e in embryos}
    ls = None if labels is None else set(labels)
    return agg(st, lambda p: p[0] in es and (ls is None or p[1] in ls) and p[2] in halves)


# ------------------------------------------------------------------ gene maps (unmapped -> -1)


def sym_index(src_genes: Sequence[str], t1_genes: Sequence[str]) -> np.ndarray:
    """Column of each T1 gene in a symbol-keyed source; symbols used more than once in the source are dropped."""
    s = pd.Series(list(src_genes))
    dup = set(s[s.duplicated(keep=False)])
    pos = {g: i for i, g in enumerate(src_genes) if g not in dup}
    return np.array([pos.get(g, -1) for g in t1_genes])


def jax_index(jax_genes: Sequence[str], t1_genes: Sequence[str], gene_table: str | Path) -> np.ndarray:
    """JAX columns are Ensembl IDs: map T1 symbols through the JAX gene table, dropping ambiguous symbols."""
    dg = pd.read_csv(gene_table)
    dg = dg[~dg.gene_short_name.duplicated(keep=False)]
    sym2id = dict(zip(dg.gene_short_name, dg.gene_id, strict=True))
    pos = {g: i for i, g in enumerate(jax_genes)}
    return np.array([pos.get(sym2id.get(g, ""), -1) for g in t1_genes])


def on_t1(vec: np.ndarray, idx: np.ndarray) -> np.ndarray:
    out = np.full(len(idx), np.nan)
    ok = idx >= 0
    out[ok] = vec[idx[ok]]
    return out


# ------------------------------------------------------------------ change vectors


def lfc_log(a: dict, b: dict) -> np.ndarray:
    """Change in mean log1p-CP10k (the challenge's pseudobulk definition)."""
    return b["log"] - a["log"]


def lfc_ratio(a: dict, b: dict, eps: float = 0.5) -> np.ndarray:
    """Log ratio of mean CP10k with a pseudocount (robust to depth changes between external groups)."""
    return np.log((b["cp"] + eps) / (a["cp"] + eps))


def zscale(v: np.ndarray) -> np.ndarray:
    """v / sd over finite entries (0 where not finite); used only for the informational z-sum combinations."""
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    out = np.zeros_like(v)
    if ok.sum() > 10 and np.std(v[ok]) > 0:
        out[ok] = v[ok] / np.std(v[ok])
    return out


# ------------------------------------------------------------------ metrics (vector forms of veckit core_metrics)


def de_direction_vec(dp: np.ndarray, dt: np.ndarray, ref: np.ndarray) -> float:
    """veckit de_direction on pseudobulk change vectors: rank-partial correlation of dp and dt given ref."""
    if np.std(dt) < 1e-12 or np.std(dp) < 1e-2 * np.std(dt):
        return 0.0
    C = np.corrcoef(np.vstack([rankdata(dp), rankdata(dt), rankdata(ref)]))
    if not np.all(np.isfinite(C)) or 1.0 - abs(C[0, 2]) < 1e-6:
        return 0.0
    num = C[0, 1] - C[0, 2] * C[1, 2]
    den = np.sqrt(max((1 - C[0, 2] ** 2) * (1 - C[1, 2] ** 2), 1e-12))
    r = num / den
    return float(np.clip(r, -1.0, 1.0)) if np.isfinite(r) else 0.0


def de_score_vec(dp: np.ndarray, dt: np.ndarray, ref: np.ndarray, thr: float = 0.25) -> tuple[float, float]:
    """veckit de_score on pseudobulk change vectors; DE genes are |dt| >= thr (veckit's docstring: its BH step is inert
    on T1 once the effect floor applies). Returns (score, chance)."""
    up, dn = np.where(dt >= thr)[0], np.where(dt <= -thr)[0]
    if len(up) + len(dn) == 0:
        return float("nan"), float("nan")
    chance = max(_signed_overlap(ref, up, dn)[0], _signed_overlap(-ref, up, dn)[0])
    if np.std(dt) < 1e-12:
        return float("nan"), float(chance)
    if np.std(dp) < 1e-2 * np.std(dt):
        return 0.0, float(chance)
    raw = _signed_overlap(dp, up, dn)[0]
    return (float((raw - chance) / (1 - chance)) if chance < 1 else float("nan")), float(chance)


def concord(d_pred: np.ndarray, d_true: np.ndarray, ref: np.ndarray, thr: float = 0.25, min_genes: int = 5) -> dict:
    """Agreement of a predicted change vector with a target change vector, both on T1 gene order. NaN in d_pred marks
    a gene the predictor says nothing about (no external counterpart, or held flat): it is left out of the sign and
    rank statistics on the target's DE genes S, and counts as a zero change in the full-gene metrics."""
    ok = np.isfinite(d_pred) & np.isfinite(d_true)
    S = ok & (np.abs(d_true) >= thr)
    r = {"n_mapped": int(ok.sum()), "n_S": int(S.sum()),
         "n_S_all": int((np.isfinite(d_true) & (np.abs(d_true) >= thr)).sum())}
    if S.sum() >= min_genes:
        r["sign_agree_S"] = float(np.mean(np.sign(d_pred[S]) == np.sign(d_true[S])))
        r["spearman_S"] = float(spearmanr(d_pred[S], d_true[S]).correlation)
    P = ok & (np.abs(d_pred) >= thr)
    r["n_P"] = int(P.sum())
    if P.sum() >= min_genes:
        r["precision_P"] = float(np.mean(np.sign(d_true[P]) == np.sign(d_pred[P])))
    r["de_direction_mapped"] = de_direction_vec(d_pred[ok], d_true[ok], ref[ok]) if ok.sum() > 2 else float("nan")
    okt = np.isfinite(d_true)
    dp_full = np.where(np.isfinite(d_pred), d_pred, 0.0)
    r["de_direction_full"] = de_direction_vec(dp_full[okt], d_true[okt], ref[okt])
    r["de_score_full"], r["de_score_chance"] = de_score_vec(dp_full[okt], d_true[okt], ref[okt], thr)
    return r


def mean_rows(rows: Sequence[Mapping]) -> dict:
    """Fold mean of the numeric fields of concord() results (NaN-aware; a field missing in a fold counts as NaN)."""
    keys = sorted({k for r in rows for k, v in r.items() if isinstance(v, int | float)})
    out = {}
    for k in keys:
        vals = np.array([float(r.get(k, np.nan)) for r in rows])
        out[k] = float(np.nanmean(vals)) if np.isfinite(vals).any() else float("nan")
    return out


# ------------------------------------------------------------------ held-flat genes


def per_label_detection(st: dict) -> dict[str, tuple[float, np.ndarray]]:
    labels = np.array([k.split("|")[0] for k in st["keys"]])
    out = {}
    for lab in np.unique(labels):
        m = labels == lab
        n = float(st["n"][m].sum())
        out[str(lab)] = (n, st["S_nz"][m].sum(0) / n)
    return out


def pool_effect_genes(st_a: dict, st_b: dict, min_cells: int, min_abs_delta: float,
                      min_type_fraction: float) -> tuple[np.ndarray, dict]:
    """Genes whose detection rate moves by at least `min_abs_delta` with one sign in at least `min_type_fraction` of
    the cell types present at both stages (exact labels, >= `min_cells` cells at each). Developmental changes are
    type-specific; a shift shared by almost every type is a property of the pooled sample (sex ratio, background)."""
    pa, pb = per_label_detection(st_a), per_label_detection(st_b)
    shared = sorted(t for t in set(pa) & set(pb) if min(pa[t][0], pb[t][0]) >= min_cells)
    if len(shared) < 3:
        raise ValueError(f"pool rule needs at least 3 shared cell types with >= {min_cells} cells; got {shared}")
    D = np.vstack([pb[t][1] - pa[t][1] for t in shared])
    k = math.ceil(min_type_fraction * len(shared) - 1e-9)
    up = (D >= min_abs_delta).sum(0) >= k
    dn = (D <= -min_abs_delta).sum(0) >= k
    return up | dn, {"shared_types": shared, "min_types": k, "median_delta": np.median(D, axis=0)}


def held_flat_genes(cfg: Mapping, genes: Sequence[str], stats_by_stage: Mapping[str, dict]) -> pd.DataFrame:
    """The held-flat list: named sex-linked, Y-linked and strain-variable genes present in T1, plus the pool rule on
    the two released stages. One row per gene with every reason that applies."""
    g = np.asarray(genes, dtype=object)
    present = set(g)
    reasons: dict[str, list[str]] = {}
    for key, tag in (("sex_linked", "sex-linked"), ("y_linked", "Y-linked"), ("strain_named", "strain (named)")):
        for name in cfg[key]:
            if name in present:
                reasons.setdefault(name, []).append(tag)
    rule = cfg["pool_rule"]
    a, b = (stats_by_stage[s] for s in rule["stages"])
    mask, info = pool_effect_genes(a, b, rule["min_cells"], rule["min_abs_delta_detection"], rule["min_type_fraction"])
    med = dict(zip(g, info["median_delta"], strict=True))
    for name in g[mask]:
        reasons.setdefault(str(name), []).append("pool rule")
    rows = [{"gene": n, "reasons": "; ".join(rs), "median_delta_detection": round(float(med[n]), 4)}
            for n, rs in reasons.items()]
    df = pd.DataFrame(rows, columns=["gene", "reasons", "median_delta_detection"])
    df.attrs["pool_rule"] = {"stages": list(rule["stages"]), "shared_types": info["shared_types"],
                             "min_types": info["min_types"]}
    return df.sort_values("gene").reset_index(drop=True)


def flat_mask(genes: Sequence[str], held: Iterable[str]) -> np.ndarray:
    hs = set(held)
    return np.array([x in hs for x in genes])


# ------------------------------------------------------------------ label harmonisation


def harmonise(labels: Iterable[str], rules: Sequence[Sequence[str]],
              override: Mapping[str, str] | None = None) -> dict[str, str]:
    """label -> group: an override wins, else the first matching regex, else 'other:<label>'."""
    override = override or {}
    out = {}
    for lab in sorted(set(labels)):
        if lab in override:
            out[lab] = override[lab]
            continue
        out[lab] = next((grp for grp, rx in rules if re.search(rx, lab)), f"other:{lab}")
    return out


def group_labels(mapping: Mapping[str, str], group: str, cm_includes: Iterable[str] = ()) -> list[str]:
    ls = [lab for lab, g in mapping.items() if g == group]
    if group == "CM":
        inc = set(cm_includes)
        ls += [lab for lab, g in mapping.items() if g in inc]
    return sorted(ls)


def antecedent_labels(mapping: Mapping[str, str], group: str, cm_includes: Iterable[str] = (),
                      fallback: Mapping[str, str] | None = None) -> list[str]:
    """The group's own labels at a stage, else those of its fallback antecedent group (e.g. Endothelium at E8.5)."""
    ls = group_labels(mapping, group, cm_includes)
    if not ls and fallback and group in fallback:
        ls = group_labels(mapping, fallback[group], cm_includes)
    return ls


# ------------------------------------------------------------------ population predictions and the decision


def population_prediction(vectors: Mapping[str, np.ndarray], weights: Mapping[str, float],
                          groups: Iterable[str], flat: np.ndarray) -> np.ndarray:
    """sum over groups of weight x per-type change (genes without a value contribute 0; held-flat genes stay 0)."""
    groups = list(groups)
    if not groups:
        raise ValueError("no groups in the population prediction")
    pred = np.zeros(len(flat))
    for g in groups:
        pred += weights[g] * np.nan_to_num(vectors[g], nan=0.0)
    pred[flat] = 0.0
    return pred


def decide(rule: Mapping, readings: Mapping) -> dict:
    """Apply the pre-registered rule. `readings` holds the fold-mean population metrics of the primary predictor and
    the baseline under the primary design, and the CM sign agreement of the primary predictor. Missing or non-finite
    readings fail their check (fail closed). A reading passes when it is >= its threshold, within 1e-9 for float
    rounding (0.25 - 0.20 is 0.04999999999999999 in binary)."""

    def check(value, threshold):
        ok = value is not None and np.isfinite(value) and value >= threshold - 1e-9
        return {"value": None if value is None or not np.isfinite(value) else float(value),
                "threshold": float(threshold), "pass": bool(ok)}

    p, b = readings.get("primary", {}), readings.get("baseline", {})
    d_score = p.get("de_score", np.nan) - b.get("de_score", np.nan)
    d_dir = p.get("de_direction", np.nan) - b.get("de_direction", np.nan)
    checks = {"delta_de_score": check(d_score, rule["min_delta_de_score"]),
              "delta_de_direction": check(d_dir, rule["min_delta_de_direction"]),
              "cm_sign_agreement": check(readings.get("cm_sign_agreement", np.nan), rule["min_cm_sign_agreement"])}
    passed = all(c["pass"] for c in checks.values())
    return {"verdict": rule["pass_verdict"] if passed else rule["fail_verdict"], "passed": passed, "checks": checks,
            "failed": [k for k, c in checks.items() if not c["pass"]]}
