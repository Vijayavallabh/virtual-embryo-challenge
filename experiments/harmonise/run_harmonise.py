"""Label harmonisation across MERFISH vocabularies (run on the DGX from the repo root).

    OMP_NUM_THREADS=16 uv run python experiments/harmonise/run_harmonise.py --out runs/harmonise_v2

Stages (--stages, default all, in this order):
  ref       WT E9.5 and WT E8.75 against their own labels: classifier candidates (stratified 5-fold CV), stratified 20%
            hold-out through fit_transfer / predict, section-grouped and slab-grouped 5-fold CV, direct group
            classifier vs summed fine-label probabilities, confidence-threshold trade-off.
  ko        Mab21l2 KO vs WT E9.5, file against file (one embryo each, sectioned differently). The KO is mapped into
            the WT vocabulary; the WT side is its own out-of-fold calls. Top-call/EM bracket (plus mean probability)
            at label, group and merged level and within CM, with 95% section-resampling intervals (within-embryo
            sampling only, not genotype effects); the 8 z-trims of the review; z-coverage profiles of both stacks;
            cross-tab against the KO's own labels (lineage map with the KO's file overrides); reverse mapping and
            mutual best label matches. Primary summary: the within-CM shift with its trim range (ko_summary.json).
  shift     composition estimators under simulated prior shift on WT E9.5 hold-out cells (needs ref).
  cross     cross-stage transfer E8.75 <-> E9.5 at group and merged level, with depth-handling variants.
  contrast  in-distribution estimator check: lower vs upper z-half of an embryo scored with a model trained on cells
            from both halves (hold-out cells; needs ref and cross).
  shifted   shifted-query estimator check, the analogue of WT-trained calls on a knockout: train on one z-half (its
            out-of-fold calls are the reference side) and predict the other half; label, group and merged level,
            E9.5 and E8.75, both directions.
  release   2026-10-20 release, input paths as arguments: Gata4 KO rep1 vs rep2 and each vs WT E8.75 (same-genotype
            embryo-to-embryo variation, the noise bar a composition shift must clear), and cross-stage transfer along
            a stage series with E8.5 and E10.5 added (adjacent pairs, both directions, group and merged level). Parts
            whose files are missing, or carry labels configs/lineage_map.yaml does not map, are skipped with exit
            status 0 and listed in release_status.json.
  figs      PNG figures from the saved tables.

Every file is loaded with vec.prep.load_stage(drop_duplicates=True): redundant copies of records are removed (counts
per file in the manifest). Tables (CSV/JSON) and figures go to --out; per-cell predictions and fitted models stay in
--out/cells (DGX only). Only released files under data/raw are read; nothing from boards/index.json is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import platform
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import classification_report, f1_score
from sklearn.model_selection import train_test_split

from vec.harmonise import (
    BRACKET,
    ESTIMATORS,
    UNASSIGNED,
    UNTRACKED,
    bootstrap_em,
    bracket_difference,
    composition,
    composition_bracket,
    cross_validate,
    default_trims,
    em_priors,
    fit_transfer,
    group_order,
    median_depth,
    predict,
    section_ids,
    shared_genes,
    subset_levels,
    transfer_metrics,
    trim_sensitivity,
    trim_summary,
    z_levels,
    z_location,
    z_profile,
)
from vec.io import COORD_KEY
from vec.prep import lineage_map, load_stage

FILES = {"E9.5": "data/raw/T2/E9.5.h5ad", "E8.75": "data/raw/T2/E8.75.h5ad", "KO": "data/raw/T3/E9.5_mab21l2_ko.h5ad"}
SERIES = ("E8.25_late=data/raw/T2/E8.25_late.h5ad,E8.5=data/raw/T2/E8.5.h5ad,E8.75=data/raw/T2/E8.75.h5ad,"
          "E9.5=data/raw/T2/E9.5.h5ad,E10.5=data/raw/T2/E10.5.h5ad")
GATA4 = {"rep1": "data/raw/T3/E8.75_gata4_ko_rep1.h5ad", "rep2": "data/raw/T3/E8.75_gata4_ko_rep2.h5ad",
         "wt": "data/raw/T2/E8.75.h5ad"}
CANDIDATES = {
    "logreg_C1_pc50": dict(method="logreg", C=1.0, n_pcs=50),
    "knn_k15_pc50": dict(method="knn", k=15, n_pcs=50),
    "logreg_C0.1_pc50": dict(method="logreg", C=0.1, n_pcs=50),
    "logreg_C1_pc30": dict(method="logreg", C=1.0, n_pcs=30),
    "logreg_C1_pc100": dict(method="logreg", C=1.0, n_pcs=100),
    "knn_k30_pc50": dict(method="knn", k=30, n_pcs=50),
}
DEFAULTS = ("logreg_C1_pc50", "knn_k15_pc50")  # the pair fit_transfer(classifier="auto") compares
THRESHOLDS = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
METRIC_KEYS = ["accuracy", "macro_f1", "group_accuracy", "group_macro_f1"]
SPEC_FILES = {"label": "labels", "group": "groups", "merged": "merged", "within_cm": "within_cm"}
KO_DEPTH = {"match": dict(depth_match=True, depth_draws=1), "match5": dict(depth_match=True, depth_draws=5),
            "none": dict(depth_match=False)}
INTERVAL_NOTE = ("95% percentile intervals from resampling whole sections (snapped z levels) within each embryo, "
                 "independently per population. They cover within-embryo section sampling only: with one embryo per "
                 "genotype they are not intervals for genotype effects, and they carry no between-embryo or "
                 "sectioning-coverage variance (see the z-trim table and the z profiles).")
T0 = time.time()
MAX_CELLS: int | None = None  # --max-cells: random subsample per file, for smoke tests only
_CACHE: dict = {}
LOADED: dict = {}


def log(msg: str) -> None:
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def jsonable(o):
    """Plain Python types for json.dumps (numpy scalars and arrays, nan -> None)."""
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, list | tuple | np.ndarray | pd.Index):
        return [jsonable(v) for v in list(o)]
    if isinstance(o, bool | np.bool_):
        return bool(o)
    if isinstance(o, int | np.integer):
        return int(o)
    if isinstance(o, float | np.floating):
        return float(o) if np.isfinite(o) else None
    return o


def load_path(path: str):
    """vec.prep.load_stage(path, drop_duplicates=True), cached per path (optionally subsampled for smoke tests)."""
    path = str(path)
    if path not in _CACHE:
        a = load_stage(path, drop_duplicates=True)
        qc = a.uns.get("qc", {})
        LOADED[path] = {"n_input": int(qc.get("n_input", a.n_obs)), "n_dropped_duplicates": int(qc.get("n_dropped", 0)),
                        "cells": int(a.n_obs), "unmapped_labels": [str(x) for x in a.uns.get("unmapped_labels", [])],
                        "sha256_16": sha256(path)}
        if MAX_CELLS is not None and a.n_obs > MAX_CELLS:
            keep = np.sort(np.random.default_rng(0).choice(a.n_obs, MAX_CELLS, replace=False))
            a = a[keep].copy()
            LOADED[path]["subsampled_to"] = int(MAX_CELLS)
        log(f"loaded {path}: {LOADED[path]['n_input']} records, {LOADED[path]['n_dropped_duplicates']} duplicate "
            f"copies dropped, {a.n_obs} cells used")
        _CACHE[path] = a
    return _CACHE[path]


def load(name: str):
    return load_path(FILES[name])


def fold_sd(oof: pd.DataFrame) -> dict:
    per = pd.DataFrame([transfer_metrics(g) for _, g in oof.groupby("fold")])
    return {f"{k}_sd": float(per[k].std(ddof=1)) for k in METRIC_KEYS}


def slab_ids(a, n_slabs: int = 5) -> np.ndarray:
    """Contiguous blocks of section levels along z (leave-one-slab-out is a harsh regional split)."""
    zl = z_levels(a)
    n = zl.levels.size
    slab_of_level = np.minimum((np.arange(n) * n_slabs) // n, n_slabs - 1)
    return np.array([f"slab{s}" for s in slab_of_level[zl.index]])


def group_confusion(true, pred, order=None) -> pd.DataFrame:
    cats = list(order) if order is not None else list(group_order(set(true) | set(pred)))
    extra = [c for c in pd.unique(np.concatenate([np.asarray(true), np.asarray(pred)])) if c not in cats]
    cats = cats + extra
    return pd.crosstab(pd.Categorical(true, cats), pd.Categorical(pred, cats), dropna=False,
                       rownames=["true"], colnames=["predicted"])


def per_group_scores(true, pred, reference_groups, order=None) -> pd.DataFrame:
    """Recall, precision and F1 per true category over the given cells."""
    true, pred = np.asarray(true).astype(str), np.asarray(pred).astype(str)
    rows = []
    for g in (order if order is not None else group_order(set(true))):
        t, p = true == g, pred == g
        if not t.any():
            continue
        rec = float((p & t).sum() / t.sum())
        prec = float((p & t).sum() / p.sum()) if p.sum() else np.nan
        f1 = 2 * rec * prec / (rec + prec) if np.isfinite(prec) and (rec + prec) > 0 else 0.0
        top = pd.Series(pred[t]).value_counts(normalize=True)
        rows.append({"group": g, "n": int(t.sum()), "in_reference": g in set(reference_groups), "recall": rec,
                     "precision": prec, "f1": f1, "top_predicted": top.index[0], "top_share": float(top.iloc[0])})
    return pd.DataFrame(rows)


def _sec(df: pd.DataFrame) -> pd.DataFrame:
    """Name bootstrap intervals as section-resampling intervals (`_sec_lo`, `_sec_hi`, `sec_excludes_0`)."""
    ren = {}
    for c in df.columns:
        n = c
        if n.endswith("_lo"):
            n = n[:-3] + "_sec_lo"
        elif n.endswith("_hi"):
            n = n[:-3] + "_sec_hi"
        if n.endswith("excludes_0"):
            n = n[: -len("excludes_0")] + "sec_excludes_0"
        ren[c] = n
    return df.rename(columns=ren)


def merged_order_of(cats) -> list[str]:
    from vec.harmonise import merged_order

    return list(merged_order(cats))


# ---------------------------------------------------------------------------------------------------------------------
# Stage ref


def stage_ref(out: Path, nj: int, seed: int) -> None:
    summary, direct, cand_all = [], [], []
    for name in ("E9.5", "E8.75"):
        a = load(name)
        sec, slab = section_ids(a), slab_ids(a)
        labels = a.obs["celltype"].astype(str).to_numpy()
        log(f"{name}: {a.n_obs} cells, {len(np.unique(sec))} sections, {len(np.unique(labels))} labels, "
            f"median depth {median_depth(a):.0f}")

        # 1. candidates, stratified 5-fold CV (scaling, PCA and classifier refitted per fold)
        for cname, kw in CANDIDATES.items():
            t = time.time()
            oof = cross_validate(a, n_folds=5, novelty=False, seed=seed, n_jobs=nj, **kw)
            row = {"stage": name, "candidate": cname, **kw, **transfer_metrics(oof), **fold_sd(oof),
                   "seconds": time.time() - t, "dropped": json.dumps(oof.attrs["dropped"])}
            cand_all.append(row)
            log(f"  {cname}: acc {row['accuracy']:.4f} macroF1 {row['macro_f1']:.4f} "
                f"group acc {row['group_accuracy']:.4f} ({row['seconds']:.0f}s)")
        cand = pd.DataFrame([r for r in cand_all if r["stage"] == name]).set_index("candidate")
        best = max(DEFAULTS, key=lambda c: cand.loc[c, "macro_f1"])
        method = CANDIDATES[best]["method"]
        log(f"  auto pair picks {best}")

        # 2. stratified 5-fold OOF of the picked default, with novelty (cached for the ko stage)
        oof = cross_validate(a, method=method, n_folds=5, seed=seed, n_jobs=nj)
        oof["section"] = pd.Series(sec, index=a.obs_names).loc[oof.index].to_numpy()
        oof.attrs["genes"] = list(a.var_names)
        oof.to_pickle(out / "cells" / f"oof_{name}.pkl")
        strat = {"stage": name, "split": "stratified 5-fold CV", "method": method, **transfer_metrics(oof),
                 **fold_sd(oof), "novelty_gt1": float((oof["novelty"] > 1).mean()),
                 "duplicates_dropped": LOADED[FILES[name]]["n_dropped_duplicates"]}
        summary.append(strat)

        # 3. grouped splits: random whole sections, and contiguous slabs of sections
        for split, grp in (("section-grouped 5-fold CV", sec), ("slab-grouped 5-fold CV (leave one of 5 z-blocks out)",
                                                                  slab)):
            g = cross_validate(a, method=method, n_folds=5, groups=grp, seed=seed, n_jobs=nj, novelty=False)
            summary.append({"stage": name, "split": split, "method": method, **transfer_metrics(g), **fold_sd(g)})
            log(f"  {split}: acc {summary[-1]['accuracy']:.4f} macroF1 {summary[-1]['macro_f1']:.4f}")

        # 4. stratified 20% hold-out through the public API
        vc = pd.Series(labels).value_counts()
        idx = np.flatnonzero((labels != "Unknown") & np.isin(labels, vc.index[vc >= 10]))
        tr, te = train_test_split(idx, test_size=0.2, stratify=labels[idx], random_state=seed)
        model = fit_transfer(a[tr], "celltype", classifier="auto", cv_folds=5, seed=seed, n_jobs=nj)
        model.cv.assign(stage=name).to_csv(out / f"ref_holdout_auto_cv_{name}.csv", index=False)
        pred = predict(model, a[te])
        pred.insert(0, "true_merged", a.obs["merged"].astype(str).to_numpy()[te])
        pred.insert(0, "true_group", a.obs["group"].astype(str).to_numpy()[te])
        pred.insert(0, "true", labels[te])
        pred["section"] = sec[te]
        met = transfer_metrics(pred)
        summary.append({"stage": name, "split": "stratified 20% hold-out", "method": model.method, **met,
                        "novelty_gt1": float((pred["novelty"] > 1).mean()), "n_train": len(tr)})
        log(f"  hold-out ({model.method}): acc {met['accuracy']:.4f} macroF1 {met['macro_f1']:.4f} "
            f"group acc {met['group_accuracy']:.4f} group macroF1 {met['group_macro_f1']:.4f} "
            f"merged acc {met['merged_accuracy']:.4f}")
        pred.to_pickle(out / "cells" / f"holdout_{name}.pkl")
        with open(out / "cells" / f"model_holdout_{name}.pkl", "wb") as fh:
            pickle.dump(model, fh)

        rep = classification_report(pred["true"], pred["label_raw"], labels=sorted(set(pred["true"])),
                                    output_dict=True, zero_division=0)
        per = pd.DataFrame({k: v for k, v in rep.items() if k in set(pred["true"])}).T
        per.index.name = "label"
        per["group"] = [model.group_of.get(lab, "Unknown") for lab in per.index]
        per.sort_values("f1-score").to_csv(out / f"ref_perclass_{name}.csv")
        group_confusion(pred["true_group"], pred["group_raw"]).to_csv(out / f"ref_group_confusion_{name}.csv")
        pg = per_group_scores(pred["true_group"], pred["group_raw"], model.groups)
        pg.to_csv(out / f"ref_pergroup_{name}.csv", index=False)

        thr_rows = []
        hit, ghit = pred["label_raw"] == pred["true"], pred["group_raw"] == pred["true_group"]
        for t in THRESHOLDS:
            ok, gok = pred["label_conf"] >= t, pred["group_conf"] >= t
            thr_rows.append({"stage": name, "threshold": t, "coverage": float(ok.mean()),
                             "assigned_accuracy": float(hit[ok].mean()), "group_coverage": float(gok.mean()),
                             "group_assigned_accuracy": float(ghit[gok].mean())})
        pd.DataFrame(thr_rows).to_csv(out / f"ref_threshold_{name}.csv", index=False)

        # 5. coarse groups: direct group classifier vs summed fine-label probabilities (same method, 5-fold)
        gmap = {g: g for g in lineage_map()["groups"]}
        dg = cross_validate(a, label_key="group", method=method, n_folds=5, seed=seed, n_jobs=nj, novelty=False,
                            group_map=gmap)
        direct += [
            {"stage": name, "approach": "summed fine-label probabilities", "group_accuracy": strat["group_accuracy"],
             "group_macro_f1": strat["group_macro_f1"]},
            {"stage": name, "approach": "direct group classifier",
             "group_accuracy": float((dg["label_raw"] == dg["true"]).mean()),
             "group_macro_f1": float(f1_score(dg["true"], dg["label_raw"], labels=np.unique(dg["true"]),
                                              average="macro", zero_division=0))},
        ]
    pd.DataFrame(cand_all).to_csv(out / "ref_candidates.csv", index=False)
    pd.DataFrame(summary).to_csv(out / "ref_summary.csv", index=False)
    pd.DataFrame(direct).to_csv(out / "ref_direct_vs_summed.csv", index=False)


# ---------------------------------------------------------------------------------------------------------------------
# File-vs-file comparison machinery (ko and release stages)


def own_composition(labels, cats) -> pd.Series:
    return pd.Series(np.asarray(labels).astype(str)).value_counts(normalize=True).reindex(cats, fill_value=0.0)


def compare_populations(out: Path, prefix: str, model, pops: dict, pairs: list, n_boot: int, seed: int,
                        cm_labels, extra_cols: dict | None = None, nj: int = 1) -> tuple[dict, dict, dict]:
    """Composition tables for populations mapped into one model's vocabulary.

    pops: {name: (predictions, ZLevels or None)}; pairs: [(a, b)] differences b - a. Writes
    {prefix}_composition_{labels,groups,merged,within_cm}.csv with, per population, top call, EM and mean probability
    (with section-resampling intervals and the top-call/EM bracket) and, per pair, the differences.
    """
    with_ci = all(zl is not None for _, zl in pops.values())
    em, seeds = {}, {}
    for i, (name, (pred, zl)) in enumerate(pops.items()):
        seeds[name] = seed + 10 * i
        t0 = time.time()
        P = pred[[f"p:{c}" for c in model.classes]].to_numpy()
        em[name] = (bootstrap_em(P, model.priors, zl.index, n_boot, seeds[name], n_jobs=nj) if with_ci
                    else (em_priors(P, model.priors)[0], None))
        log(f"{prefix}: EM for {name} ({len(pred)} cells, {n_boot if with_ci else 0} section-resampling draws) "
            f"in {time.time() - t0:.1f}s")
    specs = {"label": ("label", None), "group": ("group", None), "merged": ("merged", None)}
    if cm_labels:
        specs["within_cm"] = ("label", list(cm_labels))
    tables = {}
    for spec, (level, subset) in specs.items():
        parts = []
        if extra_cols and spec in extra_cols:
            parts.append(extra_cols[spec])
        for name, (pred, zl) in pops.items():
            b = composition_bracket(pred, model, level, zl.index if with_ci else None, subset=subset,
                                    estimators=ESTIMATORS, n_boot=n_boot, seed=seeds[name], em=em[name])
            parts.append(_sec(b).add_prefix(f"{name}_"))
        for a, b_ in pairs:
            (pa, za), (pb, zb) = pops[a], pops[b_]
            d = bracket_difference(pa, pb, model, level, za.index if with_ci else None,
                                   zb.index if with_ci else None, subset=subset, estimators=ESTIMATORS,
                                   n_boot=n_boot, seeds=(seeds[a], seeds[b_]), em_a=em[a], em_b=em[b_])
            d = d[[c for c in d.columns if not c.startswith(("a_", "b_"))]]
            parts.append(_sec(d).add_prefix(f"{b_}_vs_{a}_"))
        idx = parts[-1].index
        t = pd.concat([p.reindex(idx) for p in parts], axis=1)
        t.index.name = spec
        t.to_csv(out / f"{prefix}_composition_{SPEC_FILES[spec]}.csv")
        tables[spec] = t
    return tables, em, seeds


def trims_for(out: Path, prefix: str, model, pops: dict, pairs: list, n_boot: int, seed: int,
              cm_labels, nj: int = 1) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The 8 z-trims of the review for each pair (group, merged and within-CM; top call and EM)."""
    specs = {"group": ("group", None), "merged": ("merged", None)}
    if cm_labels:
        specs["within_cm"] = ("label", list(cm_labels))
    tabs, sums = [], []
    for a, b in pairs:
        (pa, za), (pb, zb) = pops[a], pops[b]
        if za is None or zb is None:
            continue
        t0 = time.time()
        t = trim_sensitivity(pa, za, pb, zb, model, trims=default_trims(a, b), specs=specs, estimators=BRACKET,
                             n_boot=n_boot, seed=seed, n_jobs=nj)
        s = trim_summary(t)
        s["skipped_trims"] = ";".join(t.attrs.get("skipped_trims", []))
        t.insert(0, "pair", f"{b}_vs_{a}")
        s.insert(0, "pair", f"{b}_vs_{a}")
        tabs.append(_sec(t))
        sums.append(s.rename(columns={"lo_min": "sec_lo_min", "hi_max": "sec_hi_max",
                                      "excludes_0_all_trims": "sec_excludes_0_all_trims",
                                      "n_trims_excluding_0": "n_trims_sec_excluding_0"}))
        log(f"{prefix}: z-trims {b} vs {a} ({len(t)} rows) in {time.time() - t0:.1f}s")
    T = pd.concat(tabs, ignore_index=True) if tabs else pd.DataFrame()
    S = pd.concat(sums, ignore_index=True) if sums else pd.DataFrame()
    T.to_csv(out / f"{prefix}_trim_sensitivity.csv", index=False)
    S.to_csv(out / f"{prefix}_trim_summary.csv", index=False)
    return T, S


def profiles_for(out: Path, prefix: str, model, pops: dict, pairs: list, n_bins: int = 10) -> pd.DataFrame:
    """z-coverage profiles: group composition (top call) along each stack, where each group sits, and the per-bin
    composition distance between stacks (relative z position, 0 = lowest level)."""
    prof, loc = [], []
    for name, (pred, zl) in pops.items():
        if zl is None:
            continue
        p = z_profile(pred, zl, "group", n_bins, categories=model.groups)
        p.insert(0, "population", name)
        prof.append(p)
        lc = z_location(pred, zl, "group").reindex(list(model.groups))
        lc.insert(0, "population", name)
        loc.append(lc)
    if not prof:
        return pd.DataFrame()
    P = pd.concat(prof, ignore_index=True)
    L = pd.concat(loc)
    L.index.name = "group"
    groups = list(model.groups)
    rows = []
    for a, b in pairs:
        if a not in set(P["population"]) or b not in set(P["population"]):
            continue
        pa = P[P["population"] == a].reset_index(drop=True)
        pb = P[P["population"] == b].reset_index(drop=True)
        tv = 0.5 * np.abs(pa[groups].to_numpy() - pb[groups].to_numpy()).sum(1)
        for i in range(n_bins):
            rows.append({"pair": f"{b}_vs_{a}", "bin_lo": pa.at[i, "bin_lo"], "bin_hi": pa.at[i, "bin_hi"],
                         f"cells_{a}": int(pa.at[i, "n_cells"]), f"cells_{b}": int(pb.at[i, "n_cells"]),
                         "total_variation": float(tv[i])})
        la, lb = L[L["population"] == a], L[L["population"] == b]
        shift = (lb["mean_pos"] - la["mean_pos"]).rename(f"{b}_minus_{a}_mean_pos")
        L = L.join(shift, how="left") if shift.name not in L else L
    D = pd.DataFrame(rows)
    P.to_csv(out / f"{prefix}_zprofile.csv", index=False)
    L.to_csv(out / f"{prefix}_zlocation.csv")
    D.to_csv(out / f"{prefix}_zprofile_distance.csv", index=False)
    return D


def concordance(own_label, own_group, own_merged, pred, model) -> pd.DataFrame:
    """Each own label of a query file against its calls in the model vocabulary (label, group and merged level)."""
    rows = []
    for lab in pd.unique(own_label):
        m = own_label == lab
        fine = pd.Series(pred["label_raw"].to_numpy()[m]).value_counts(normalize=True)
        grp = pd.Series(pred["group_raw"].to_numpy()[m]).value_counts(normalize=True)
        mrg = pd.Series(pred["merged_raw"].to_numpy()[m]).value_counts(normalize=True)
        g_own = pd.Series(own_group[m]).mode().iloc[0]
        m_own = pd.Series(own_merged[m]).mode().iloc[0]
        rows.append({
            "own_label": lab, "n": int(m.sum()), "lineage_group": g_own, "lineage_merged": m_own,
            "group_in_reference_vocab": g_own in set(model.groups),
            "merged_tracked_in_reference": m_own != UNTRACKED and m_own in set(model.merged),
            "top_label": fine.index[0], "top_label_share": float(fine.iloc[0]),
            "second_label": fine.index[1] if len(fine) > 1 else "",
            "second_share": float(fine.iloc[1]) if len(fine) > 1 else 0.0,
            "top_group": grp.index[0], "top_group_share": float(grp.iloc[0]),
            "agrees_with_lineage_map": float((pred["group_raw"].to_numpy()[m] == g_own).mean()),
            "top_merged": mrg.index[0], "top_merged_share": float(mrg.iloc[0]),
            "agrees_with_merged": float((pred["merged_raw"].to_numpy()[m] == m_own).mean()),
            "unassigned_at_0.5": float((pred["label"].to_numpy()[m] == UNASSIGNED).mean()),
            "novelty_gt1": float((pred["novelty"].to_numpy()[m] > 1).mean()),
            "median_conf": float(np.median(pred["label_conf"].to_numpy()[m])),
        })
    return pd.DataFrame(rows).sort_values("n", ascending=False)


def agreement(own_group, own_merged, pred, model) -> dict:
    gcomp = np.isin(own_group, list(model.groups))
    mcomp = (own_merged != UNTRACKED) & np.isin(own_merged, list(model.merged))
    return {"group_comparable_fraction": float(gcomp.mean()),
            "group_agreement": float((pred["group_raw"].to_numpy()[gcomp] == own_group[gcomp]).mean()),
            "merged_comparable_fraction": float(mcomp.mean()),
            "merged_agreement": float((pred["merged_raw"].to_numpy()[mcomp] == own_merged[mcomp]).mean())}


def within_cm_block(tab: pd.DataFrame, tsum: pd.DataFrame, a: str, b: str) -> dict:
    """The within-CM shift of b vs a with its section intervals and its range over the z-trims."""
    pair = f"{b}_vs_{a}"
    out = {}
    for lab in tab.index:
        r = tab.loc[lab]
        ts = tsum[(tsum["pair"] == pair) & (tsum["spec"] == "within_cm") & (tsum["category"] == lab)]
        ts = ts.set_index("estimator")
        rec = {f"{a}_top_call": r[f"{a}_top_call"], f"{b}_top_call": r[f"{b}_top_call"],
               f"{a}_em": r[f"{a}_em"], f"{b}_em": r[f"{b}_em"]}
        for e in BRACKET:
            rec[f"{e}_diff"] = r[f"{pair}_{e}_diff"]
            if f"{pair}_{e}_diff_sec_lo" in r:
                rec[f"{e}_diff_sec_ci"] = [r[f"{pair}_{e}_diff_sec_lo"], r[f"{pair}_{e}_diff_sec_hi"]]
            if e in ts.index:
                rec[f"{e}_trim_range"] = [ts.at[e, "diff_min"], ts.at[e, "diff_max"]]
                rec[f"{e}_sec_excludes_0_all_trims"] = bool(ts.at[e, "sec_excludes_0_all_trims"])
                rec[f"{e}_sign_flips_over_trims"] = bool(ts.at[e, "sign_flips"])
        rec["bracket"] = [r[f"{pair}_bracket_low"], r[f"{pair}_bracket_high"]]
        out[lab] = rec
    return out


def group_block(tab: pd.DataFrame, tsum: pd.DataFrame, a: str, b: str, spec: str = "group") -> dict:
    pair = f"{b}_vs_{a}"
    out = {}
    for g in tab.index:
        r = tab.loc[g]
        ts = tsum[(tsum["pair"] == pair) & (tsum["spec"] == spec) & (tsum["category"] == g)].set_index("estimator")
        rec = {"bracket": [r[f"{pair}_bracket_low"], r[f"{pair}_bracket_high"]],
               "top_call_diff": r[f"{pair}_top_call_diff"], "em_diff": r[f"{pair}_em_diff"],
               "bracket_sec_excludes_0": bool(r.get(f"{pair}_bracket_sec_excludes_0", False))}
        for e in BRACKET:
            if e in ts.index:
                rec[f"{e}_trim_range"] = [ts.at[e, "diff_min"], ts.at[e, "diff_max"]]
                rec[f"{e}_sec_excludes_0_all_trims"] = bool(ts.at[e, "sec_excludes_0_all_trims"])
        out[g] = rec
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Stage ko


def stage_ko(out: Path, nj: int, seed: int, n_boot: int, n_boot_trim: int, ko_depth: str = "match5") -> None:
    wt, ko = load("E9.5"), load("KO")
    genes = shared_genes(wt, ko)
    d_wt, d_ko = median_depth(wt, genes), median_depth(ko, genes)
    log(f"WT E9.5 {wt.n_obs} cells (depth {d_wt:.0f}); KO {ko.n_obs} cells (depth {d_ko:.0f}); "
        f"{len(genes)} shared genes, same order: {genes == list(wt.var_names) == list(ko.var_names)}")
    model = fit_transfer(wt, "celltype", genes=genes, classifier="auto", cv_folds=5, seed=seed, n_jobs=nj)
    model.cv.to_csv(out / "ko_model_cv.csv", index=False)
    with open(out / "cells" / "model_E9.5_full.pkl", "wb") as fh:
        pickle.dump(model, fh)
    log(f"WT E9.5 model: {model.method}, {len(model.classes)} labels, groups {list(model.groups)}, "
        f"merged {list(model.merged)}, dropped {model.dropped}")

    variants = {k: predict(model, ko, seed=seed, **kw) for k, kw in KO_DEPTH.items()}
    pred = variants[ko_depth]
    log(f"KO primary variant '{ko_depth}': kept {pred.attrs['depth_frac']:.3f} of its molecules "
        f"(target depth {model.depth:.0f}, {pred.attrs['depth_draws']} draw(s))")
    pred.assign(ko_label=ko.obs["celltype"].astype(str).to_numpy()).to_pickle(out / "cells" / "ko_pred.pkl")

    path = out / "cells" / "oof_E9.5.pkl"
    oof = pd.read_pickle(path) if path.exists() else None
    if oof is None or oof.attrs.get("method") != model.method or oof.attrs.get("genes") != genes:
        oof = cross_validate(wt, method=model.method, genes=genes, n_folds=5, seed=seed, n_jobs=nj)
    zl_wt_all, zl_ko = z_levels(wt), z_levels(ko)
    pos = wt.obs_names.get_indexer(oof.index)
    zl_wt = subset_levels(zl_wt_all, pos)
    log(f"sections: WT {zl_wt_all.levels.size} (step {zl_wt_all.step:.3f}), KO {zl_ko.levels.size} "
        f"(step {zl_ko.step:.3f}, gaps {zl_ko.gaps()}); WT OOF cells {len(oof)}")

    cm_labels = [c for c in model.classes if model.group_of[c] == "CM"]
    wt_label = wt.obs["celltype"].astype(str).to_numpy()[pos]
    wt_group = wt.obs["group"].astype(str).to_numpy()[pos]
    wt_merged = wt.obs["merged"].astype(str).to_numpy()[pos]
    ko_label = ko.obs["celltype"].astype(str).to_numpy()
    ko_group = ko.obs["group"].astype(str).to_numpy()
    ko_merged = ko.obs["merged"].astype(str).to_numpy()
    in_cm = np.isin(wt_label, cm_labels)
    extra = {
        "label": pd.DataFrame({"WT_annotated": own_composition(wt_label, model.classes)}),
        "group": pd.DataFrame({"WT_annotated": own_composition(wt_group, model.groups),
                               "KO_own_labels_lineage_map": own_composition(ko_group, model.groups)}),
        "merged": pd.DataFrame({"WT_annotated": own_composition(wt_merged, model.merged),
                                "KO_own_labels_lineage_map": own_composition(ko_merged, model.merged)}),
        "within_cm": pd.DataFrame({"WT_annotated": own_composition(wt_label[in_cm], cm_labels)}),
    }
    pops = {"WT": (oof, zl_wt), "KO": (pred, zl_ko)}
    pairs = [("WT", "KO")]
    t0 = time.time()
    tables, em, seeds = compare_populations(out, "ko", model, pops, pairs, n_boot, seed, cm_labels, extra, nj)
    log(f"composition tables in {time.time() - t0:.0f}s")

    # sensitivity columns: KO without depth matching (point estimates)
    for spec, (level, subset) in {"label": ("label", None), "group": ("group", None), "merged": ("merged", None),
                                  "within_cm": ("label", cm_labels)}.items():
        nd = composition_bracket(variants["none"], model, level, subset=subset)
        t = tables[spec]
        t["KO_no_depth_match_top_call"] = nd["top_call"].reindex(t.index)
        t["KO_no_depth_match_em"] = nd["em"].reindex(t.index)
        t.to_csv(out / f"ko_composition_{SPEC_FILES[spec]}.csv")
    for spec in ("label", "group", "within_cm"):
        t = tables[spec]
        log(f"{spec}: section-resampling 95% CI of KO-WT excludes 0 for "
            f"{int(t['KO_vs_WT_top_call_sec_excludes_0'].sum())} (top call) / "
            f"{int(t['KO_vs_WT_em_sec_excludes_0'].sum())} (EM) / "
            f"{int(t['KO_vs_WT_bracket_sec_excludes_0'].sum())} (both) of {len(t)} categories")

    T, S = trims_for(out, "ko", model, pops, pairs, n_boot_trim, seed, cm_labels, nj)
    D = profiles_for(out, "ko", model, pops, pairs)

    # cross-tab against the KO's own labels (lineage map with the KO's file overrides, from load_stage)
    labels = list(model.classes)
    ct = pd.crosstab(pd.Series(ko_label, name="ko_label"), pd.Series(pred["label"].to_numpy(), name="predicted"))
    ct = ct[[c for c in labels + [UNASSIGNED] if c in ct.columns]]
    ct.to_csv(out / "ko_crosstab_counts.csv")
    ct.div(ct.sum(1), axis=0).to_csv(out / "ko_crosstab_rowfrac.csv")
    ctt = pd.crosstab(pd.Series(ko_label, name="ko_label"), pd.Series(pred["label_raw"].to_numpy(), name="top_call"))
    ctt = ctt[[c for c in labels if c in ctt.columns]]
    ctt.to_csv(out / "ko_crosstab_top_counts.csv")
    ctt.div(ctt.sum(1), axis=0).to_csv(out / "ko_crosstab_top_rowfrac.csv")
    group_confusion(ko_group, pred["group_raw"].to_numpy()).to_csv(out / "ko_group_confusion.csv")
    group_confusion(ko_merged, pred["merged_raw"].to_numpy(),
                    order=merged_order_of(set(ko_merged) | set(pred["merged_raw"]))).to_csv(
        out / "ko_merged_confusion.csv")
    conc = concordance(ko_label, ko_group, ko_merged, pred, model)
    conc.rename(columns={"own_label": "ko_label"}).to_csv(out / "ko_label_concordance.csv", index=False)

    summary = {
        "files": {k: LOADED[FILES[k]] for k in ("E9.5", "KO")},
        "wt_depth": d_wt, "ko_depth": d_ko, "ko_depth_frac": variants["match"].attrs["depth_frac"],
        "ko_primary_variant": ko_depth, "method": model.method, "n_boot": n_boot, "n_boot_em": n_boot,
        "n_boot_trim": n_boot_trim, "wt_sections": int(zl_wt_all.levels.size), "ko_sections": int(zl_ko.levels.size),
        "wt_z_step": zl_wt_all.step, "ko_z_step": zl_ko.step, "ko_z_gaps": zl_ko.gaps(),
        "interval_note": INTERVAL_NOTE,
        "reading": ("File against file: one WT embryo and one KO embryo whose z-stacks cover different anatomy. "
                    "Group-level differences are not genotype effects or T3 magnitudes. The primary summary is the "
                    "within-CM shift with its z-trim range. Compositions are quoted as a top-call/EM bracket."),
        "primary_within_cm": within_cm_block(tables["within_cm"], S, "WT", "KO"),
        "file_vs_file_groups": group_block(tables["group"], S, "WT", "KO", "group"),
        "file_vs_file_merged": group_block(tables["merged"], S, "WT", "KO", "merged"),
        "zprofile_mean_total_variation": float(D["total_variation"].mean()) if len(D) else None,
    }
    for tag, p in variants.items():
        ag = agreement(ko_group, ko_merged, p, model)
        for k, v in ag.items():
            summary[f"{k}_{tag}"] = v
        summary[f"unassigned_{tag}"] = float((p["label"] == UNASSIGNED).mean())
        summary[f"novelty_gt1_{tag}"] = float((p["novelty"] > 1).mean())
        summary[f"median_label_conf_{tag}"] = float(p["label_conf"].median())
        summary[f"label_change_vs_primary_{tag}"] = float((p["label_raw"] != pred["label_raw"]).mean())
    log("KO own-group agreement on comparable cells: " + ", ".join(
        f"{k} {summary[f'group_agreement_{k}']:.3f} (merged {summary[f'merged_agreement_{k}']:.3f})"
        for k in variants))

    # reverse mapping: WT cells into the KO vocabulary (KO thinned to WT depth), and mutual best matches
    rev = fit_transfer(ko, "celltype", genes=genes, classifier=model.method, depth=d_wt, seed=seed, n_jobs=nj)
    rpred = predict(rev, wt)
    wt_label_all = wt.obs["celltype"].astype(str).to_numpy()
    rct = pd.crosstab(pd.Series(wt_label_all, name="wt_label"), pd.Series(rpred["label_raw"].to_numpy(),
                                                                           name="predicted_ko_label"))
    rct.to_csv(out / "ko_reverse_crosstab_counts.csv")
    rfrac = rct.div(rct.sum(1), axis=0)
    rfrac.to_csv(out / "ko_reverse_crosstab_rowfrac.csv")
    ffrac = ctt.div(ctt.sum(1), axis=0)
    matches = []
    for lab in ffrac.index:
        w = ffrac.loc[lab].idxmax()
        back = rfrac.loc[w].idxmax() if w in rfrac.index else ""
        matches.append({"ko_label": lab, "best_wt_label": w, "forward_share": float(ffrac.loc[lab, w]),
                        "best_ko_label_for_that_wt_label": back,
                        "reverse_share": float(rfrac.loc[w, back]) if back else np.nan, "mutual": back == lab})
    pd.DataFrame(matches).to_csv(out / "ko_vocabulary_matches.csv", index=False)
    summary["reverse_model_depth"] = rev.depth
    summary["reverse_method"] = rev.method
    summary["reverse_model_groups"] = {k: rev.group_of[k] for k in ("LPM", "DPM") if k in rev.group_of}
    (out / "ko_summary.json").write_text(json.dumps(jsonable(summary), indent=2))
    for lab, r in summary["primary_within_cm"].items():
        tr, te = r["top_call_trim_range"], r["em_trim_range"]
        log(f"within-CM {lab}: WT {r['WT_top_call']:.3f} -> KO {r['KO_top_call']:.3f}; top-call diff "
            f"{r['top_call_diff']:+.3f} (trims {tr[0]:+.3f}..{tr[1]:+.3f}), "
            f"EM diff {r['em_diff']:+.3f} (trims {te[0]:+.3f}..{te[1]:+.3f})")


# ---------------------------------------------------------------------------------------------------------------------
# Stage shift


def stage_shift(out: Path, seed: int, n_rep: int = 200, n_draw: int = 5000, sigma: float = 0.5) -> None:
    pred = pd.read_pickle(out / "cells" / "holdout_E9.5.pkl")
    with open(out / "cells" / "model_holdout_E9.5.pkl", "rb") as fh:
        model = pickle.load(fh)
    classes, M = list(model.classes), model.membership
    K = len(classes)
    P = pred[[f"p:{c}" for c in classes]].to_numpy()
    yi = np.searchsorted(model.classes, pred["true"].to_numpy())
    ri = np.searchsorted(model.classes, pred["label_raw"].to_numpy())
    gi = np.argmax(M[yi], axis=1)
    rng = np.random.default_rng(seed)

    def est(sel):
        return {"true": np.bincount(yi[sel], minlength=K) / len(sel),
                "top call": np.bincount(ri[sel], minlength=K) / len(sel),
                "mean probability": P[sel].mean(0),
                "EM": em_priors(P[sel], model.priors, tol=1e-7)[0]}

    base = est(np.arange(len(yi)))
    recs = []
    for regime in ("group-level shift", "label-level shift"):
        for r in range(n_rep):
            w = (rng.lognormal(0.0, sigma, M.shape[1])[gi] if regime.startswith("group")
                 else rng.lognormal(0.0, sigma, K)[yi])
            sel = rng.choice(len(yi), n_draw, replace=True, p=w / w.sum())
            e = est(sel)
            for level in ("label", "group"):
                f = (lambda v: v @ M) if level == "group" else (lambda v: v)
                T, Tb = f(e["true"]), f(base["true"])
                for m in ("top call", "mean probability", "EM"):
                    E, Eb = f(e[m]), f(base[m])
                    recs.append({"regime": regime, "level": level, "estimator": m, "rep": r,
                                 "mae": float(np.abs(E - T).mean()), "max_err": float(np.abs(E - T).max()),
                                 "num": float(((E - Eb) * (T - Tb)).sum()), "den": float(((T - Tb) ** 2).sum())})
    df = pd.DataFrame(recs)
    summ = df.groupby(["regime", "level", "estimator"]).agg(
        mae=("mae", "mean"), max_err=("max_err", "mean"), num=("num", "sum"), den=("den", "sum")).reset_index()
    summ["delta_slope"] = summ["num"] / summ["den"]
    summ.drop(columns=["num", "den"]).to_csv(out / "shift_estimators.csv", index=False)
    log("shift: " + "; ".join(f"{r.regime}/{r.level}/{r.estimator}: mae {r.mae:.4f} slope {r.delta_slope:.3f}"
                              for r in summ.itertuples()))


# ---------------------------------------------------------------------------------------------------------------------
# Cross-stage transfer (cross and release stages)


def transfer_scores(tag: str, m, q, p) -> tuple[dict, pd.DataFrame, dict]:
    """Group- and merged-level agreement of a transfer with the query's own labels through the lineage map."""
    true = q.obs["group"].astype(str).to_numpy()
    true_m = q.obs["merged"].astype(str).to_numpy()
    known = true != "Unknown"
    transferable = known & np.isin(true, m.groups)
    m_transferable = known & (true_m != UNTRACKED) & np.isin(true_m, m.merged)
    g_hat, m_hat = p["group_raw"].to_numpy(), p["merged_raw"].to_numpy()

    def f1(t, h, sel):
        if not sel.any():  # e.g. a release file whose labels are all Unknown or untracked
            return float("nan")
        return float(f1_score(t[sel], h[sel], labels=np.unique(t[sel]), average="macro", zero_division=0))

    row = {
        "direction": tag, "method": m.method, "depth_frac": p.attrs["depth_frac"],
        "depth_draws": p.attrs["depth_draws"], "reference_depth": m.depth, "n_genes": len(m.genes),
        "n_query": int(q.n_obs), "known_group_cells": int(known.sum()),
        "transferable_fraction": float(transferable.sum() / max(1, known.sum())),
        "group_accuracy_transferable": float((g_hat[transferable] == true[transferable]).mean()),
        "group_macro_f1_transferable": f1(true, g_hat, transferable),
        "group_accuracy_all_known": float((g_hat[known] == true[known]).mean()),
        "merged_transferable_fraction": float(m_transferable.sum() / max(1, known.sum())),
        "merged_accuracy_transferable": float((m_hat[m_transferable] == true_m[m_transferable]).mean()),
        "merged_macro_f1_transferable": f1(true_m, m_hat, m_transferable),
        "merged_accuracy_all_known": float((m_hat[known] == true_m[known]).mean()),
        "unassigned_label_0.5": float((p["label"] == UNASSIGNED).mean()),
        "unassigned_group_0.5": float((p["group"] == UNASSIGNED).mean()),
        "novelty_gt1": float((p["novelty"] > 1).mean()),
    }
    parts = []
    for level, t, h, cats in (("group", true, g_hat, m.groups), ("merged", true_m, m_hat, m.merged)):
        order = (group_order(set(t[known])) if level == "group" else merged_order_of(set(t[known])))
        pg = per_group_scores(t[known], h[known], cats, order=order)
        sub = pd.DataFrame({"g": t[known], "nov": p["novelty"].to_numpy()[known],
                            "un": p["group"].to_numpy()[known] == UNASSIGNED})
        pg["novelty_gt1"] = pg["group"].map(sub.groupby("g")["nov"].apply(lambda v: float((v > 1).mean())))
        pg["unassigned_group_0.5"] = pg["group"].map(sub.groupby("g")["un"].mean())
        pg.insert(0, "level", level)
        pg.insert(0, "direction", tag)
        parts.append(pg)
    conf = {"group": group_confusion(true, g_hat),
            "merged": group_confusion(true_m, m_hat, order=merged_order_of(set(true_m) | set(m_hat)))}
    return row, pd.concat(parts, ignore_index=True), conf


def verdict_table(pg: pd.DataFrame, fwd: str, bwd: str, cats_fwd, cats_bwd, level: str) -> pd.DataFrame:
    """Per category: F1 each way and a verdict (clean >= 0.8, partial >= 0.5, else poor; one stage only)."""
    a = pg[(pg["direction"] == fwd) & (pg["level"] == level)].set_index("group")
    b = pg[(pg["direction"] == bwd) & (pg["level"] == level)].set_index("group")
    present = set(a.index) | set(b.index)
    allg = group_order(present) if level == "group" else merged_order_of(present)
    rows = []
    for g in allg:
        in_ref_a, in_ref_b = g in set(cats_fwd), g in set(cats_bwd)
        fa = a["f1"].get(g, np.nan) if in_ref_a else np.nan
        fb = b["f1"].get(g, np.nan) if in_ref_b else np.nan
        if g == UNTRACKED:
            v = "untracked (no cross-stage meaning)"
        elif not (in_ref_a and in_ref_b) or g not in a.index or g not in b.index:
            v = "one stage only"
        elif min(fa, fb) >= 0.8:
            v = "clean"
        elif min(fa, fb) >= 0.5:
            v = "partial"
        else:
            v = "poor"
        rows.append({"level": level, "category": g, "forward": fwd, "backward": bwd,
                     "f1_forward": fa, "recall_forward": a["recall"].get(g, np.nan),
                     "precision_forward": a["precision"].get(g, np.nan),
                     "f1_backward": fb, "recall_backward": b["recall"].get(g, np.nan),
                     "precision_backward": b["precision"].get(g, np.nan), "verdict": v})
    return pd.DataFrame(rows)


def safe_tag(tag: str) -> str:
    return tag.replace("->", "_to_").replace("(", "_").replace(")", "").replace(" ", "_")


def stage_cross(out: Path, nj: int, seed: int) -> None:
    e875, e95 = load("E8.75"), load("E9.5")
    genes = shared_genes(e875, e95)
    d875, d95 = median_depth(e875, genes), median_depth(e95, genes)
    log(f"E8.75 depth {d875:.0f}, E9.5 depth {d95:.0f}, {len(genes)} shared genes")
    mA = fit_transfer(e875, "celltype", genes=genes, classifier="auto", cv_folds=5, seed=seed, n_jobs=nj)
    mA_thin = fit_transfer(e875, "celltype", genes=genes, classifier=mA.method, depth=d95, seed=seed, n_jobs=nj)
    mB = fit_transfer(e95, "celltype", genes=genes, classifier="auto", cv_folds=5, seed=seed, n_jobs=nj)
    log(f"E8.75 model {mA.method} (dropped {mA.dropped}); E9.5 model {mB.method}")
    with open(out / "cells" / "models_cross.pkl", "wb") as fh:
        pickle.dump({"E8.75": mA, "E8.75_thinned": mA_thin, "E9.5": mB}, fh)
    variants = {
        "E8.75->E9.5": (mA, e95, {}),
        "E8.75(thinned)->E9.5": (mA_thin, e95, {}),
        "E9.5->E8.75(1 thinning draw)": (mB, e875, {"depth_draws": 1}),
        "E9.5->E8.75(5 thinning draws)": (mB, e875, {"depth_draws": 5}),
        "E9.5->E8.75(no depth match)": (mB, e875, {"depth_match": False}),
    }
    rows, pergroup = [], []
    for tag, (m, q, kw) in variants.items():
        p = predict(m, q, seed=seed, **kw)
        row, pg, conf = transfer_scores(tag, m, q, p)
        rows.append(row)
        pergroup.append(pg)
        safe = safe_tag(tag)
        keep_cols = ["label_raw", "group_raw", "merged_raw", "label_conf", "group_conf", "novelty"] + [
            c for c in p.columns if c.startswith(("p:", "pg:", "pm:"))]
        p[keep_cols].assign(true_label=q.obs["celltype"].astype(str).to_numpy(),
                            true_group=q.obs["group"].astype(str).to_numpy(),
                            true_merged=q.obs["merged"].astype(str).to_numpy(),
                            section=section_ids(q)).to_pickle(out / "cells" / f"cross_{safe}.pkl")
        conf["group"].to_csv(out / f"cross_group_confusion_{safe}.csv")
        conf["merged"].to_csv(out / f"cross_merged_confusion_{safe}.csv")
        ct = pd.crosstab(pd.Series(q.obs["celltype"].astype(str).to_numpy(), name="query_label"),
                         pd.Series(p["label_raw"].to_numpy(), name="reference_label"))
        ct.div(ct.sum(1), axis=0).to_csv(out / f"cross_crosstab_rowfrac_{safe}.csv")
        log(f"{tag}: transferable {row['transferable_fraction']:.3f}, group acc "
            f"{row['group_accuracy_transferable']:.4f}, macroF1 {row['group_macro_f1_transferable']:.4f}; merged acc "
            f"{row['merged_accuracy_transferable']:.4f}, macroF1 {row['merged_macro_f1_transferable']:.4f}")
    pd.DataFrame(rows).to_csv(out / "cross_summary.csv", index=False)
    pg = pd.concat(pergroup, ignore_index=True)
    pg[pg["level"] == "group"].drop(columns="level").to_csv(out / "cross_pergroup.csv", index=False)
    pg[pg["level"] == "merged"].drop(columns="level").to_csv(out / "cross_pergroup_merged.csv", index=False)
    fwd, bwd = "E8.75(thinned)->E9.5", "E9.5->E8.75(5 thinning draws)"
    vg = verdict_table(pg, fwd, bwd, mA.groups, mB.groups, "group")
    vg.to_csv(out / "cross_verdict.csv", index=False)
    vm = verdict_table(pg, fwd, bwd, mA.merged, mB.merged, "merged")
    vm.to_csv(out / "cross_verdict_merged.csv", index=False)
    log("merged verdicts: " + ", ".join(f"{r.category} {r.verdict}" for r in vm.itertuples()))


# ---------------------------------------------------------------------------------------------------------------------
# Estimator checks: in-distribution contrast and shifted query


def _contrast_rows(setting: str, pred_a, truth_a: dict, pred_b, truth_b: dict, model, levels, extra=None) -> list:
    """Composition change from population a to b: estimated (top call, mean probability, EM) vs true.

    truth_*: {level: array of true categories}. Top call at a level is the argmax of the summed probabilities
    (`{level}_raw`), as in the knockout comparison.
    """
    rows = []
    for level in levels:
        ta, tb = np.asarray(truth_a[level]).astype(str), np.asarray(truth_b[level]).astype(str)
        cats = sorted(set(ta) | set(tb) | {str(c) for c in model.categories(level)})
        Ta, Tb = own_composition(ta, cats).to_numpy(), own_composition(tb, cats).to_numpy()
        dT = Tb - Ta
        for name, method in (("top call", "raw"), ("mean probability", "soft"), ("EM", "em")):
            Ea = composition(pred_a, level, method, model).reindex(cats, fill_value=0.0).to_numpy()
            Eb = composition(pred_b, level, method, model).reindex(cats, fill_value=0.0).to_numpy()
            dE = Eb - Ea
            rows.append({"setting": setting, "level": level, "estimator": name, "n_a": len(ta), "n_b": len(tb),
                         "true_contrast_tv": float(np.abs(dT).sum() / 2),
                         "delta_slope": float(dE @ dT / (dT @ dT)), "delta_mae": float(np.abs(dE - dT).mean()),
                         "delta_max_err": float(np.abs(dE - dT).max()),
                         "composition_mae_a": float(np.abs(Ea - Ta).mean()),
                         "composition_mae_b": float(np.abs(Eb - Tb).mean()), **(extra or {})})
    return rows


def _halves(sections) -> np.ndarray:
    """True for the upper half of the section levels present (z above the median level)."""
    z = np.array([float(s[1:]) for s in sections])
    return z > np.median(np.unique(z))


def stage_contrast(out: Path) -> None:
    rows = []
    for st in ("E9.5", "E8.75"):
        pred = pd.read_pickle(out / "cells" / f"holdout_{st}.pkl")
        with open(out / "cells" / f"model_holdout_{st}.pkl", "rb") as fh:
            model = pickle.load(fh)
        pred = pred[np.isin(pred["true"], model.classes)]
        up = _halves(pred["section"].to_numpy())
        truth = {"label": pred["true"].to_numpy(), "group": pred["true_group"].to_numpy(),
                 "merged": pred["true_merged"].to_numpy()}
        rows += _contrast_rows(f"in-distribution within {st}: 20% hold-out cells, model trained on both halves",
                               pred[~up], {k: v[~up] for k, v in truth.items()}, pred[up],
                               {k: v[up] for k, v in truth.items()}, model, ("label", "group", "merged"))
    with open(out / "cells" / "models_cross.pkl", "rb") as fh:
        models = pickle.load(fh)
    for tag, key in (("E8.75_thinned_to_E9.5", "E8.75_thinned"), ("E9.5_to_E8.75_5_thinning_draws", "E9.5")):
        p = pd.read_pickle(out / "cells" / f"cross_{tag}.pkl")
        m = models[key]
        p = p[np.isin(p["true_group"], m.groups)]  # oracle restriction to groups the reference has
        up = _halves(p["section"].to_numpy())
        tg = p["true_group"].to_numpy()
        rows += _contrast_rows(f"cross-stage {tag}: transferable cells", p[~up], {"group": tg[~up]}, p[up],
                               {"group": tg[up]}, m, ("group",))
    df = pd.DataFrame(rows)
    df.to_csv(out / "contrast_estimators.csv", index=False)
    log("contrast: " + "; ".join(f"{r.setting.split(':')[0]}/{r.level}/{r.estimator}: slope {r.delta_slope:.3f} "
                                 f"mae {r.delta_mae:.4f}" for r in df.itertuples()))


def stage_shifted(out: Path, nj: int, seed: int) -> None:
    """Train on one z-half (out-of-fold calls on it), predict the other half: label, group, merged level."""
    rows, info = [], []
    for st in ("E9.5", "E8.75"):
        a = load(st)
        zl = z_levels(a)
        upper = zl.levels[zl.index] > np.median(zl.levels)
        lab = a.obs["celltype"].astype(str).to_numpy()
        ok = lab != "Unknown"
        for train_half, query_half in (("lower", "upper"), ("upper", "lower")):
            trm = (upper if train_half == "upper" else ~upper) & ok
            qm = (~upper if train_half == "upper" else upper) & ok
            t0 = time.time()
            ref = a[trm]
            oof = cross_validate(ref, method="logreg", n_folds=5, seed=seed, n_jobs=nj, novelty=False)
            model = fit_transfer(ref, "celltype", classifier="logreg", seed=seed, n_jobs=nj)
            q = a[qm]
            pq = predict(model, q, seed=seed)
            truth_q = {"label": lab[qm], "group": q.obs["group"].astype(str).to_numpy(),
                       "merged": q.obs["merged"].astype(str).to_numpy()}
            truth_tr = {"label": oof["true"].to_numpy(), "group": oof["true_group"].to_numpy(),
                        "merged": oof["true_merged"].to_numpy()}
            untrained = float((~np.isin(lab[qm], model.classes)).mean())
            inside = np.isin(lab[qm], model.classes)
            q_acc = float((pq["label_raw"].to_numpy()[inside] == lab[qm][inside]).mean())
            extra = {"stage": st, "train_half": train_half, "query_half": query_half,
                     "query_cells_of_untrained_labels": untrained, "query_accuracy": q_acc,
                     "depth_frac": pq.attrs["depth_frac"]}
            setting = f"shifted query {st}: train {train_half} half (out-of-fold), predict {query_half} half"
            rows += _contrast_rows(setting, oof, truth_tr, pq, truth_q, model, ("label", "group", "merged"), extra)
            info.append({**extra, "n_train": int(len(oof)), "n_query": int(qm.sum()),
                         "dropped_labels": json.dumps(model.dropped), "seconds": time.time() - t0})
            log(f"shifted {st} train {train_half}: query accuracy {q_acc:.3f}, untrained-label share {untrained:.3f}, "
                f"depth frac {pq.attrs['depth_frac']:.3f} ({time.time() - t0:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(out / "contrast_shifted_query.csv", index=False)
    pd.DataFrame(info).to_csv(out / "contrast_shifted_query_settings.csv", index=False)
    summ = df.groupby(["level", "estimator"]).agg(
        slope_mean=("delta_slope", "mean"), slope_min=("delta_slope", "min"), slope_max=("delta_slope", "max"),
        mae_mean=("delta_mae", "mean"), n_settings=("delta_slope", "size")).reset_index()
    best = df.loc[df.groupby(["setting", "level"])["delta_mae"].idxmin(), ["level", "estimator"]].value_counts()
    summ["settings_with_lowest_mae"] = [int(best.get((r.level, r.estimator), 0)) for r in summ.itertuples()]
    summ.to_csv(out / "contrast_shifted_query_summary.csv", index=False)
    log("shifted query: " + "; ".join(f"{r.level}/{r.estimator}: slope {r.slope_mean:.3f} "
                                      f"[{r.slope_min:.3f}, {r.slope_max:.3f}] mae {r.mae_mean:.4f}"
                                      for r in summ.itertuples()))


# ---------------------------------------------------------------------------------------------------------------------
# Stage release (2026-10-20)


def discover(words) -> list[str]:
    """h5ad files under data/raw whose names contain any of `words` (case-insensitive), for the status file."""
    root = Path("data/raw")
    if not root.exists():
        return []
    return sorted(str(p) for p in root.rglob("*.h5ad") if any(w.lower() in p.name.lower() for w in words))


def parse_series(spec: str) -> list[tuple[str, str]]:
    out = []
    for item in spec.split(","):
        name, _, path = item.partition("=")
        if not path:
            raise SystemExit(f"--series entries are NAME=PATH, got {item!r}")
        out.append((name.strip(), path.strip()))
    return out


def run_gata4(out: Path, adatas: dict, nj: int, seed: int, n_boot: int, n_boot_trim: int) -> dict:
    wt, r1, r2 = adatas["WT"], adatas["rep1"], adatas["rep2"]
    genes = shared_genes(wt, r1, r2)
    depths = {k: median_depth(a, genes) for k, a in adatas.items()}
    target = min(d for d in depths.values() if d is not None)
    depth = target if depths["WT"] is not None and depths["WT"] > 1.05 * target else None
    log(f"gata4: {len(genes)} shared genes (WT panel {wt.n_vars}); depths {depths}; WT reference thinned to "
        f"{depth if depth is not None else 'no thinning'}")
    model = fit_transfer(wt, "celltype", genes=genes, classifier="auto", cv_folds=5, depth=depth, seed=seed,
                         n_jobs=nj)
    oof = cross_validate(wt, method=model.method, genes=genes, depth=depth, n_folds=5, seed=seed, n_jobs=nj)
    preds = {k: predict(model, adatas[k], seed=seed) for k in ("rep1", "rep2")}
    for k, p in preds.items():
        p.to_pickle(out / "cells" / f"gata4_{k}_pred.pkl")
    coords = all(COORD_KEY in a.obsm for a in adatas.values())
    pos = wt.obs_names.get_indexer(oof.index)
    zl = {"WT": subset_levels(z_levels(wt), pos) if coords else None,
          "rep1": z_levels(r1) if coords else None, "rep2": z_levels(r2) if coords else None}
    pops = {"WT": (oof, zl["WT"]), "rep1": (preds["rep1"], zl["rep1"]), "rep2": (preds["rep2"], zl["rep2"])}
    pairs = [("WT", "rep1"), ("WT", "rep2"), ("rep1", "rep2")]
    cm_labels = [c for c in model.classes if model.group_of[c] == "CM"]
    wt_label = wt.obs["celltype"].astype(str).to_numpy()[pos]
    in_cm = np.isin(wt_label, cm_labels)
    extra = {"label": pd.DataFrame({"WT_annotated": own_composition(wt_label, model.classes)}),
             "group": pd.DataFrame({"WT_annotated": own_composition(wt.obs["group"].astype(str).to_numpy()[pos],
                                                                    model.groups)}),
             "merged": pd.DataFrame({"WT_annotated": own_composition(wt.obs["merged"].astype(str).to_numpy()[pos],
                                                                     model.merged)}),
             "within_cm": pd.DataFrame({"WT_annotated": own_composition(wt_label[in_cm], cm_labels)})}
    tables, _, _ = compare_populations(out, "gata4", model, pops, pairs, n_boot, seed, cm_labels, extra, nj)
    T, S = trims_for(out, "gata4", model, pops, pairs, n_boot_trim, seed, cm_labels, nj) if coords else (None, None)
    D = profiles_for(out, "gata4", model, pops, pairs) if coords else pd.DataFrame()

    # noise bar: does each rep-vs-WT shift exceed the rep1-vs-rep2 difference, with the same sign in both reps?
    bars = []
    for spec, t in tables.items():
        for cat in t.index:
            d = {f"{b}_vs_{a}_{e}": t.at[cat, f"{b}_vs_{a}_{e}_diff"] for a, b in pairs for e in BRACKET}
            noise = max(abs(d["rep2_vs_rep1_top_call"]), abs(d["rep2_vs_rep1_em"]))
            shifts = [d[f"{r}_vs_WT_{e}"] for r in ("rep1", "rep2") for e in BRACKET]
            same = bool(np.all(np.sign(shifts) == np.sign(shifts[0])) and shifts[0] != 0)
            smallest = float(min(abs(s) for s in shifts)) if same else 0.0
            sec = [f"{r}_vs_WT_bracket_sec_excludes_0" for r in ("rep1", "rep2")]
            bars.append({"spec": spec, "category": cat, **d, "rep_to_rep_noise": noise,
                         "smallest_shift_same_sign": smallest, "same_sign_both_reps_both_estimators": same,
                         "clears_rep_noise": bool(same and smallest > noise),
                         "ratio_shift_to_noise": smallest / noise if noise > 0 else np.nan,
                         "both_reps_bracket_sec_excludes_0": bool(all(t.at[cat, c] for c in sec))
                         if all(c in t for c in sec) else None})
    B = pd.DataFrame(bars)
    B.to_csv(out / "gata4_noise_bar.csv", index=False)

    status = {"ran": True, "method": model.method, "genes": len(genes), "depths": depths,
              "reference_depth": model.depth, "interval_note": INTERVAL_NOTE,
              "coords": coords, "n_boot": n_boot, "n_boot_trim": n_boot_trim}
    for k, p in preds.items():
        a = adatas[k]
        status[f"{k}_depth_frac"] = p.attrs["depth_frac"]
        if "celltype" in a.obs:
            c = concordance(a.obs["celltype"].astype(str).to_numpy(), a.obs["group"].astype(str).to_numpy(),
                            a.obs["merged"].astype(str).to_numpy(), p, model)
            c.to_csv(out / f"gata4_{k}_label_concordance.csv", index=False)
            status[f"{k}_agreement"] = agreement(a.obs["group"].astype(str).to_numpy(),
                                                 a.obs["merged"].astype(str).to_numpy(), p, model)
    if S is not None and len(S):
        status["within_cm"] = {f"{b}_vs_{a}": within_cm_block(tables["within_cm"], S, a, b) for a, b in pairs}
    for spec in ("group", "merged"):
        sub = B[B["spec"] == spec]
        status[f"{spec}_clearing_rep_noise"] = sub.loc[sub["clears_rep_noise"], "category"].tolist()
        status[f"{spec}_median_rep_to_rep_noise"] = float(sub["rep_to_rep_noise"].median())
    status["zprofile_mean_total_variation"] = ({p: float(g["total_variation"].mean()) for p, g in D.groupby("pair")}
                                               if len(D) else None)
    (out / "gata4_summary.json").write_text(json.dumps(jsonable(status), indent=2))
    log(f"gata4: categories clearing the rep-to-rep noise: group {status['group_clearing_rep_noise']}, "
        f"merged {status['merged_clearing_rep_noise']}")
    return status


def run_series(out: Path, stages: list[tuple[str, str]], nj: int, seed: int, pairs_mode: str) -> dict:
    loaded = {n: load_path(p) for n, p in stages}
    names = [n for n, _ in stages]
    pairs = (list(zip(names[:-1], names[1:], strict=True)) if pairs_mode == "adjacent"
             else list(combinations(names, 2)))
    models: dict = {}
    rows, pergroup, verdicts = [], [], []
    for x, y in pairs:
        tags = []
        for rn, qn in ((x, y), (y, x)):
            ref, q = loaded[rn], loaded[qn]
            genes = shared_genes(ref, q)
            d_ref, d_q = median_depth(ref, genes), median_depth(q, genes)
            depth = d_q if (d_ref is not None and d_q is not None and d_ref > 1.05 * d_q) else None
            key = (rn, hashlib.sha1("|".join(genes).encode()).hexdigest(), None if depth is None else round(depth, 3))
            if key not in models:
                t0 = time.time()
                models[key] = fit_transfer(ref, "celltype", genes=genes, classifier="auto", cv_max_cells=20000,
                                           depth=depth, seed=seed, n_jobs=nj)
                log(f"series: {rn} model ({models[key].method}, depth {models[key].depth:.0f}, {len(genes)} genes) "
                    f"in {time.time() - t0:.0f}s")
            m = models[key]
            p = predict(m, q, seed=seed)
            tag = f"{rn}->{qn}"
            row, pg, conf = transfer_scores(tag, m, q, p)
            row.update(pair=f"{x}|{y}", reference=rn, query=qn, reference_depth_raw=d_ref, query_depth=d_q)
            rows.append(row)
            pergroup.append(pg)
            conf["group"].to_csv(out / f"release_cross_group_confusion_{safe_tag(tag)}.csv")
            conf["merged"].to_csv(out / f"release_cross_merged_confusion_{safe_tag(tag)}.csv")
            tags.append((tag, m))
            log(f"series {tag}: group acc {row['group_accuracy_transferable']:.4f} (transferable "
                f"{row['transferable_fraction']:.3f}), merged acc {row['merged_accuracy_transferable']:.4f}")
        pg_all = pd.concat(pergroup, ignore_index=True)
        (fwd, mf), (bwd, mb) = tags
        for level, cf, cb in (("group", mf.groups, mb.groups), ("merged", mf.merged, mb.merged)):
            v = verdict_table(pg_all, fwd, bwd, cf, cb, level)
            v.insert(0, "pair", f"{x}|{y}")
            verdicts.append(v)
    pd.DataFrame(rows).to_csv(out / "release_cross_summary.csv", index=False)
    pd.concat(pergroup, ignore_index=True).to_csv(out / "release_cross_pergroup.csv", index=False)
    V = pd.concat(verdicts, ignore_index=True)
    V.to_csv(out / "release_cross_verdict.csv", index=False)
    return {"ran": True, "stages": names, "pairs": [f"{x}|{y}" for x, y in pairs],
            "files": {n: LOADED.get(str(p), {}) for n, p in stages}}


def stage_release(out: Path, nj: int, seed: int, n_boot: int, n_boot_trim: int, args) -> None:
    status: dict = {"checked": time.strftime("%Y-%m-%d %H:%M:%S"), "inputs": {}}
    errors: list[str] = []
    g_paths = {"WT": args.gata4_wt, "rep1": args.gata4_rep1, "rep2": args.gata4_rep2}
    series = parse_series(args.series)
    for k, p in g_paths.items():
        status["inputs"][f"gata4_{k}"] = {"path": p, "exists": Path(p).exists()}
    for n, p in series:
        status["inputs"][f"series_{n}"] = {"path": p, "exists": Path(p).exists()}

    def unmapped_of(paths: dict) -> dict:
        return {k: LOADED[str(p)]["unmapped_labels"] for k, p in paths.items() if LOADED[str(p)]["unmapped_labels"]}

    missing = [k for k, p in g_paths.items() if not Path(p).exists()]
    if missing:
        status["gata4"] = {"ran": False, "reason": f"input files not present: {missing}",
                           "candidates_under_data_raw": discover(["gata4"])}
        log(f"release/gata4 skipped: {status['gata4']['reason']} (candidates: "
            f"{status['gata4']['candidates_under_data_raw']})")
    else:
        adatas = {k: load_path(p) for k, p in g_paths.items()}
        unmapped = unmapped_of(g_paths)
        if unmapped and not args.allow_unmapped:
            status["gata4"] = {"ran": False, "unmapped_labels": unmapped,
                               "reason": "labels not in configs/lineage_map.yaml; the pre-registered rule says map "
                                         "and commit them before running anything on the file (or pass "
                                         "--allow-unmapped)"}
            log(f"release/gata4 skipped: unmapped labels {unmapped}")
        else:
            try:
                status["gata4"] = run_gata4(out, adatas, nj, seed, n_boot, n_boot_trim)
                status["gata4"]["unmapped_labels"] = unmapped
                status["gata4"]["files"] = {k: {"path": str(p), **LOADED[str(p)]} for k, p in g_paths.items()}
            except Exception as e:  # record and keep going, so release_status.json is always written
                status["gata4"] = {"ran": False, "error": repr(e)}
                errors.append(f"gata4: {e!r}")

    present = [(n, p) for n, p in series if Path(p).exists()]
    absent = [n for n, p in series if not Path(p).exists()]
    if len(present) < 2:
        status["series"] = {"ran": False, "reason": f"fewer than two stage files present (missing {absent})"}
    else:
        for _, p in present:
            load_path(p)
        unmapped = unmapped_of(dict(present))
        if unmapped and not args.allow_unmapped:
            keep = [(n, p) for n, p in present if n not in unmapped]
            status["series_unmapped_skipped"] = unmapped
            log(f"release/series: leaving out stages with unmapped labels {unmapped}")
            present = keep
        if len(present) < 2:
            status["series"] = {"ran": False, "reason": "fewer than two usable stage files"}
        else:
            try:
                status["series"] = run_series(out, present, nj, seed, args.release_pairs)
            except Exception as e:
                status["series"] = {"ran": False, "error": repr(e)}
                errors.append(f"series: {e!r}")
    status["series_missing"] = absent
    status["release_files_present"] = {"gata4": not missing,
                                       "series": [n for n, _ in present if n in ("E8.5", "E10.5")]}
    (out / "release_status.json").write_text(json.dumps(jsonable(status), indent=2))
    log(f"release: gata4 ran {status['gata4']['ran']}; series ran {status['series']['ran']} "
        f"({status['series'].get('pairs', [])}); missing series files {absent}")
    if errors:
        raise SystemExit(f"release stage errors (recorded in release_status.json): {errors}")


# ---------------------------------------------------------------------------------------------------------------------
# Figures

INK, INK2, MUTED, GRID, BASE, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def _style(plt):
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": BASE, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": INK2,
        "text.color": INK, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
        "font.family": ["DejaVu Sans"],
    })


def _heatmap(ax, frac: pd.DataFrame, title: str, xlabel: str, ylabel: str, annotate_min: float = 0.05):
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ)
    im = ax.imshow(frac.to_numpy(dtype=float), cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(frac.shape[1]), frac.columns, rotation=60, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(frac.shape[0]), frac.index)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    vals = frac.to_numpy(dtype=float)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            if np.isfinite(v) and v >= annotate_min:
                ax.text(j, i, f"{v:.2f}".lstrip("0") if v < 1 else "1", ha="center", va="center", fontsize=6.5,
                        color="white" if v > 0.55 else INK)
    return im


def _grid(ax, axis="x"):
    ax.grid(axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _bracket_dots(ax, t, pops, colors, labels, dy=0.16):
    """Top call (filled dot) and EM (open diamond) per population with section-resampling intervals."""
    y = np.arange(len(t))
    offs = np.linspace(dy, -dy, len(pops)) if len(pops) > 1 else [0.0]
    for p, c, lab, o in zip(pops, colors, labels, offs, strict=True):
        ax.hlines(y + o, t[f"{p}_top_call_sec_lo"], t[f"{p}_top_call_sec_hi"], color=c, linewidth=2)
        ax.hlines(y + o, t[[f"{p}_top_call", f"{p}_em"]].min(1), t[[f"{p}_top_call", f"{p}_em"]].max(1), color=c,
                  linewidth=1, alpha=0.5)
        ax.scatter(t[f"{p}_top_call"], y + o, s=40, color=c, edgecolor=SURFACE, linewidth=1.5, zorder=3,
                   label=f"{lab}: top call")
        ax.scatter(t[f"{p}_em"], y + o, s=40, marker="D", facecolor=SURFACE, edgecolor=c, linewidth=1.5, zorder=4,
                   label=f"{lab}: EM")
    ax.set_yticks(y, t.index)


def stage_figs(out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _style(plt)

    def rowfrac(df):
        return df.div(df.sum(1).replace(0, np.nan), axis=0)

    # 1. hold-out group confusion
    paths = [out / f"ref_group_confusion_{s}.csv" for s in ("E9.5", "E8.75")]
    if all(p.exists() for p in paths):
        fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), constrained_layout=True)
        for ax, p, s in zip(axes, paths, ("E9.5", "E8.75"), strict=True):
            cm = pd.read_csv(p, index_col=0)
            cm = cm.loc[cm.sum(1) > 0, :]
            cm = cm.loc[:, (cm.sum(0) > 0) | cm.columns.isin(cm.index)]
            im = _heatmap(ax, rowfrac(cm), f"WT {s}: 20% hold-out, lineage groups (row = true)",
                          "predicted group", "true group")
        fig.colorbar(im, ax=axes, shrink=0.5, label="fraction of the true group's cells", aspect=30)
        fig.savefig(out / "fig_ref_group_confusion.png", dpi=160)
        plt.close(fig)

    # 2. KO vs WT group composition: top-call/EM bracket with section-resampling intervals
    p = out / "ko_composition_groups.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        t = t.loc[t[["WT_top_call", "KO_top_call", "WT_annotated"]].max(1) > 0.002].sort_values("WT_top_call")
        fig, ax = plt.subplots(figsize=(8, 0.46 * len(t) + 1.8), constrained_layout=True)
        _grid(ax)
        _bracket_dots(ax, t, ["WT", "KO"], [BLUE, ORANGE], ["WT E9.5 (out-of-fold)", "Mab21l2 KO"])
        ax.scatter(t["WT_annotated"], np.arange(len(t)) + 0.16, s=30, facecolor="none", edgecolor=INK2,
                   linewidth=1.0, zorder=5, label="WT E9.5: its own labels")
        ax.set_xlabel("fraction of cells (bars: 95% section-resampling interval of the top call, within embryo)")
        ax.set_title("WT E9.5 vs Mab21l2 KO, WT vocabulary\n(one embryo per genotype: file vs file)", loc="left")
        ax.legend(frameon=False, loc="lower right", fontsize=7.5)
        fig.savefig(out / "fig_ko_vs_wt_groups.png", dpi=160)
        plt.close(fig)

    # 3. KO own labels vs top-call WT labels
    p = out / "ko_crosstab_top_rowfrac.csv"
    if p.exists():
        ct = pd.read_csv(p, index_col=0)
        counts = pd.read_csv(out / "ko_crosstab_top_counts.csv", index_col=0)
        ct = ct.loc[counts.sum(1).sort_values(ascending=False).index]
        fig, ax = plt.subplots(figsize=(11, 6.2), constrained_layout=True)
        im = _heatmap(ax, ct, "Mab21l2 KO: own labels (rows) mapped to WT E9.5 labels (top call, no threshold)",
                      "WT E9.5 label (top call)", "KO file label")
        fig.colorbar(im, ax=ax, shrink=0.5, label="fraction of the KO label's cells", aspect=30)
        fig.savefig(out / "fig_ko_crosstab.png", dpi=160)
        plt.close(fig)

    # 4. cross-stage group confusion
    tags = ["E8.75_thinned_to_E9.5", "E9.5_to_E8.75_5_thinning_draws"]
    paths = [out / f"cross_group_confusion_{t}.csv" for t in tags]
    if all(p.exists() for p in paths):
        fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), constrained_layout=True)
        titles = ["E8.75 model (thinned to E9.5 depth) on E9.5",
                  "E9.5 model on E8.75 (thinned to E9.5 depth, 5 draws)"]
        for ax, p, title in zip(axes, paths, titles, strict=True):
            cm = pd.read_csv(p, index_col=0)
            cm = cm.drop(index=[i for i in cm.index if i == "Unknown"], errors="ignore")
            cm = cm.loc[cm.sum(1) > 0, cm.sum(0) > 0]
            im = _heatmap(ax, rowfrac(cm), title, "predicted group (reference vocabulary)",
                          "true group (query labels)")
        fig.colorbar(im, ax=axes, shrink=0.5, label="fraction of the true group's cells", aspect=30)
        fig.savefig(out / "fig_cross_stage_groups.png", dpi=160)
        plt.close(fig)

    # 5. z-coverage profiles of the WT and KO stacks
    p = out / "ko_zprofile.csv"
    if p.exists():
        prof = pd.read_csv(p)
        groups = [c for c in prof.columns if c not in ("population", "bin_lo", "bin_hi", "n_cells", "n_levels")]
        ncol = 4
        nrow = int(np.ceil(len(groups) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(12, 2.3 * nrow + 0.8), sharex=True, constrained_layout=True)
        for ax, g in zip(np.ravel(axes), groups, strict=False):
            _grid(ax, "y")
            for pop, c, lab in (("WT", BLUE, "WT E9.5"), ("KO", ORANGE, "Mab21l2 KO")):
                s = prof[prof["population"] == pop]
                ax.plot((s["bin_lo"] + s["bin_hi"]) / 2, s[g], color=c, linewidth=2, marker="o", markersize=4,
                        label=lab)
            ax.set_title(g, loc="left", fontsize=9)
            ax.set_ylim(bottom=0)
        for ax in np.ravel(axes)[len(groups):]:
            ax.set_visible(False)
        np.ravel(axes)[0].legend(frameon=False, fontsize=7.5)
        fig.supxlabel("relative z position in the stack (0 = lowest section level, 1 = highest)", color=INK2,
                      fontsize=9)
        fig.supylabel("fraction of cells in the z bin (top call)", color=INK2, fontsize=9)
        fig.suptitle("Where each group sits along the WT and KO stacks (coverage check before comparing files)",
                     x=0.01, ha="left", fontweight="bold", fontsize=10)
        fig.savefig(out / "fig_ko_zprofile.png", dpi=160)
        plt.close(fig)

    # 6. within-CM KO-WT differences under the 8 z-trims
    p = out / "ko_trim_sensitivity.csv"
    if p.exists():
        t = pd.read_csv(p)
        t = t[t["spec"] == "within_cm"]
        labs = list(pd.unique(t["category"]))
        trims = list(pd.unique(t["trim"]))
        fig, axes = plt.subplots(1, len(labs), figsize=(3.2 * len(labs) + 1.5, 3.6), sharey=True,
                                 constrained_layout=True)
        for ax, lab in zip(np.atleast_1d(axes), labs, strict=True):
            _grid(ax)
            ax.axvline(0, color=BASE, linewidth=1)
            for e, c, dy, mk in (("top_call", BLUE, 0.13, "o"), ("em", ORANGE, -0.13, "D")):
                s = t[(t["category"] == lab) & (t["estimator"] == e)].set_index("trim").reindex(trims)
                y = np.arange(len(trims)) + dy
                ax.hlines(y, s["diff_sec_lo"], s["diff_sec_hi"], color=c, linewidth=2)
                ax.scatter(s["diff"], y, s=36, marker=mk, color=c, edgecolor=SURFACE, linewidth=1.2, zorder=3,
                           label={"top_call": "top call", "em": "EM"}[e])
            ax.set_title(f"{lab} (within CM)", loc="left", fontsize=9)
            ax.set_yticks(np.arange(len(trims)), trims)
            ax.invert_yaxis()
        np.atleast_1d(axes)[0].legend(frameon=False, fontsize=7.5, loc="lower left")
        fig.supxlabel("KO - WT fraction within cardiomyocytes (bars: 95% section-resampling interval)", color=INK2,
                      fontsize=9)
        fig.savefig(out / "fig_ko_within_cm_trims.png", dpi=160)
        plt.close(fig)

    # 7. estimator slopes: in-distribution vs shifted query
    pa, pb = out / "contrast_estimators.csv", out / "contrast_shifted_query.csv"
    if pa.exists() and pb.exists():
        a = pd.read_csv(pa)
        a = a[a["setting"].str.startswith("in-distribution")].assign(kind="in distribution")
        b = pd.read_csv(pb).assign(kind="shifted query")
        d = pd.concat([a, b], ignore_index=True)
        levels = [lv for lv in ("label", "group", "merged") if lv in set(d["level"])]
        fig, axes = plt.subplots(1, len(levels), figsize=(4.2 * len(levels), 3.2), sharey=True,
                                 constrained_layout=True)
        for ax, lv in zip(np.atleast_1d(axes), levels, strict=True):
            _grid(ax)
            ax.axvline(1, color=BASE, linewidth=1)
            for j, (kind, c) in enumerate((("in distribution", BLUE), ("shifted query", ORANGE))):
                for i, e in enumerate(("top call", "mean probability", "EM")):
                    s = d[(d["level"] == lv) & (d["kind"] == kind) & (d["estimator"] == e)]["delta_slope"]
                    yy = np.full(len(s), i + (0.15 if j == 0 else -0.15))
                    ax.scatter(s, yy, s=30, color=c, edgecolor=SURFACE, linewidth=1, zorder=3,
                               label=kind if i == 0 else None)
            ax.set_yticks(range(3), ["top call", "mean probability", "EM"])
            ax.set_title(f"{lv} level", loc="left", fontsize=9)
            ax.invert_yaxis()
        np.atleast_1d(axes)[0].legend(frameon=False, fontsize=7.5, loc="lower left")
        fig.supxlabel("estimated / true composition change (slope; 1 = unbiased; one dot per z-half setting)",
                      color=INK2, fontsize=9)
        fig.savefig(out / "fig_estimator_slopes.png", dpi=160)
        plt.close(fig)

    # 8. Gata4 replicates (release day)
    p = out / "gata4_composition_groups.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        t = t.loc[t[["WT_top_call", "rep1_top_call", "rep2_top_call"]].max(1) > 0.002].sort_values("WT_top_call")
        fig, ax = plt.subplots(figsize=(8, 0.5 * len(t) + 1.8), constrained_layout=True)
        _grid(ax)
        ax.axvline(0, color=BASE, linewidth=1)
        y = np.arange(len(t))
        for (pair, c, lab), o in zip((("rep1_vs_WT", BLUE, "rep1 - WT"), ("rep2_vs_WT", ORANGE, "rep2 - WT"),
                                      ("rep2_vs_rep1", AQUA, "rep2 - rep1 (embryo-to-embryo)")),
                                     (0.2, 0.0, -0.2), strict=True):
            if f"{pair}_top_call_diff_sec_lo" in t:
                ax.hlines(y + o, t[f"{pair}_envelope_sec_lo"], t[f"{pair}_envelope_sec_hi"], color=c, linewidth=1.5)
            ax.hlines(y + o, t[f"{pair}_bracket_low"], t[f"{pair}_bracket_high"], color=c, linewidth=4)
            ax.scatter(t[f"{pair}_top_call_diff"], y + o, s=30, color=c, edgecolor=SURFACE, linewidth=1, zorder=3,
                       label=lab)
        ax.set_yticks(y, t.index)
        ax.set_xlabel("difference in group fraction (dot: top call; thick: top-call/EM bracket; thin: interval)")
        ax.set_title("Gata4 KO replicates vs WT E8.75 and vs each other", loc="left")
        ax.legend(frameon=False, loc="lower right", fontsize=7.5)
        fig.savefig(out / "fig_gata4_groups.png", dpi=160)
        plt.close(fig)
    log("figures written")


# ---------------------------------------------------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stages", default="ref,ko,shift,cross,contrast,shifted,release,figs")
    ap.add_argument("--out", default="runs/harmonise_v2")
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--n-boot", type=int, default=2000, help="section-resampling draws (top call, EM, mean prob)")
    ap.add_argument("--n-boot-trim", type=int, default=None, help="draws per z-trim (default: --n-boot)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-cells", type=int, default=None, help="subsample every file (smoke tests only)")
    ap.add_argument("--ko-depth", choices=sorted(KO_DEPTH), default="match5",
                    help="depth handling for the KO predictions used in the composition tables")
    ap.add_argument("--gata4-rep1", default=GATA4["rep1"])
    ap.add_argument("--gata4-rep2", default=GATA4["rep2"])
    ap.add_argument("--gata4-wt", default=GATA4["wt"], help="matched wild type for the Gata4 replicates")
    ap.add_argument("--series", default=SERIES, help="NAME=PATH,... stage files in developmental order")
    ap.add_argument("--release-pairs", choices=["adjacent", "all"], default="adjacent")
    ap.add_argument("--allow-unmapped", action="store_true",
                    help="run release parts on files with labels the lineage map does not cover")
    args = ap.parse_args()
    global MAX_CELLS
    MAX_CELLS = args.max_cells
    n_boot_trim = args.n_boot if args.n_boot_trim is None else args.n_boot_trim
    out = Path(args.out)
    (out / "cells").mkdir(parents=True, exist_ok=True)
    manifest = {
        "started": time.strftime("%Y-%m-%d %H:%M:%S"), "stages": args.stages, "seed": args.seed,
        "n_boot": args.n_boot, "n_boot_em": args.n_boot, "n_boot_trim": n_boot_trim, "max_cells": args.max_cells,
        "ko_depth": args.ko_depth, "load": "vec.prep.load_stage(drop_duplicates=True)",
        "python": platform.python_version(), "sklearn": sklearn.__version__,
        "harmonise_py_sha256_16": sha256("src/vec/harmonise.py"), "spatial_py_sha256_16": sha256("src/vec/spatial.py"),
        "prep_py_sha256_16": sha256("src/vec/prep.py"), "lineage_map_sha256_16": sha256("configs/lineage_map.yaml"),
        "script_sha256_16": sha256(__file__),
        "release_args": {"gata4_rep1": args.gata4_rep1, "gata4_rep2": args.gata4_rep2, "gata4_wt": args.gata4_wt,
                         "series": args.series, "release_pairs": args.release_pairs,
                         "allow_unmapped": args.allow_unmapped},
    }
    for stage in args.stages.split(","):
        log(f"=== stage {stage}")
        t0 = time.time()
        if stage == "ref":
            stage_ref(out, args.n_jobs, args.seed)
        elif stage == "ko":
            stage_ko(out, args.n_jobs, args.seed, args.n_boot, n_boot_trim, args.ko_depth)
        elif stage == "shift":
            stage_shift(out, args.seed)
        elif stage == "cross":
            stage_cross(out, args.n_jobs, args.seed)
        elif stage == "contrast":
            stage_contrast(out)
        elif stage == "shifted":
            stage_shifted(out, args.n_jobs, args.seed)
        elif stage == "release":
            stage_release(out, args.n_jobs, args.seed, args.n_boot, n_boot_trim, args)
        elif stage == "figs":
            stage_figs(out)
        else:
            raise SystemExit(f"unknown stage {stage}")
        manifest.setdefault("stage_seconds", {})[stage] = round(time.time() - t0, 1)
    manifest["files"] = LOADED
    manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    manifest["seconds"] = round(time.time() - T0, 1)
    (out / f"manifest_{args.stages.replace(',', '-')}.json").write_text(json.dumps(jsonable(manifest), indent=2))
    log("done")


if __name__ == "__main__":
    main()
