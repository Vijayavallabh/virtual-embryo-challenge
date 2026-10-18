"""Label harmonisation across MERFISH vocabularies (run on the DGX from the repo root).

    uv run python experiments/harmonise/run_harmonise.py --stages ref,ko,shift,cross,figs --out runs/harmonise

Stages:
  ref    WT E9.5 and WT E8.75 against their own labels: classifier candidates (stratified 5-fold CV), stratified 20%
         hold-out through the public API (fit_transfer / predict), section-grouped and slab-grouped 5-fold CV,
         direct group classifier vs summed fine-label probabilities, confidence-threshold trade-off.
  ko     Mab21l2 KO mapped into the WT E9.5 vocabulary: harmonised composition KO vs WT (fine labels, groups,
         within-CM) with section-bootstrap intervals, cross-tab against the KO's own labels, depth-matching
         ablation, reverse mapping (WT cells into the KO vocabulary) and mutual best label matches.
  shift  composition estimators (top call, mean probability, EM) under simulated prior shift on WT E9.5 hold-out
         cells (needs ref).
  cross  cross-stage transfer at the lineage-group level: E8.75 model on E9.5 and E9.5 model on E8.75, with
         depth-handling variants.
  contrast  composition estimators on a real composition contrast with known truth: the two z-halves of one embryo
         (within-stage hold-out cells, and cross-stage transferable cells; needs ref and cross).
  figs   PNG figures from the saved tables.

Tables (CSV/JSON) and figures go to --out; per-cell predictions and fitted models stay in --out/cells (DGX only).
Only released training files are read; nothing from boards/index.json is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import classification_report, f1_score
from sklearn.model_selection import train_test_split

from vec.harmonise import (
    UNASSIGNED,
    bootstrap_composition,
    bootstrap_difference,
    composition,
    cross_validate,
    em_priors,
    fit_transfer,
    group_order,
    one_hot,
    predict,
    section_ids,
    transfer_metrics,
)
from vec.prep import lineage_map, load_stage

FILES = {"E9.5": "data/raw/T2/E9.5.h5ad", "E8.75": "data/raw/T2/E8.75.h5ad", "KO": "data/raw/T3/E9.5_mab21l2_ko.h5ad"}
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
T0 = time.time()
MAX_CELLS: int | None = None  # --max-cells: random subsample per file, for smoke tests only


def load(name: str):
    """vec.prep.load_stage on a released file (optionally subsampled at random for smoke tests)."""
    a = load_stage(FILES[name])
    if MAX_CELLS is not None and a.n_obs > MAX_CELLS:
        keep = np.sort(np.random.default_rng(0).choice(a.n_obs, MAX_CELLS, replace=False))
        a = a[keep].copy()
    return a


def log(msg: str) -> None:
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def fold_sd(oof: pd.DataFrame) -> dict:
    per = pd.DataFrame([transfer_metrics(g) for _, g in oof.groupby("fold")])
    return {f"{k}_sd": float(per[k].std(ddof=1)) for k in METRIC_KEYS}


def slab_ids(a, n_slabs: int = 5) -> np.ndarray:
    """Contiguous blocks of sections along z (leave-one-slab-out is a harsh regional split)."""
    z = np.round(np.asarray(a.obsm["spatial_3D"])[:, 2].astype(float), 0)
    levels = np.unique(z)
    slab_of_level = np.minimum((np.arange(len(levels)) * n_slabs) // len(levels), n_slabs - 1)
    return np.array([f"slab{s}" for s in slab_of_level[np.searchsorted(levels, z)]])


def median_depth(a) -> float:
    return float(np.median(np.asarray(a.layers["counts"].sum(1)).ravel()))


def group_confusion(true, pred, order=None) -> pd.DataFrame:
    cats = list(order) if order is not None else list(group_order(set(true) | set(pred)))
    extra = [c for c in pd.unique(np.concatenate([np.asarray(true), np.asarray(pred)])) if c not in cats]
    cats = cats + extra
    return pd.crosstab(pd.Categorical(true, cats), pd.Categorical(pred, cats), dropna=False,
                       rownames=["true"], colnames=["predicted"])


def per_group_scores(true, pred, reference_groups) -> pd.DataFrame:
    """Recall, precision and F1 per true group over all cells with a known group."""
    true, pred = np.asarray(true).astype(str), np.asarray(pred).astype(str)
    rows = []
    for g in group_order(set(true)):
        t, p = true == g, pred == g
        rec = float((p & t).sum() / t.sum()) if t.sum() else np.nan
        prec = float((p & t).sum() / p.sum()) if p.sum() else np.nan
        f1 = 2 * rec * prec / (rec + prec) if (rec + prec) and np.isfinite(rec) and np.isfinite(prec) else 0.0
        top = pd.Series(pred[t]).value_counts(normalize=True)
        rows.append({"group": g, "n": int(t.sum()), "in_reference": g in set(reference_groups), "recall": rec,
                     "precision": prec, "f1": f1, "top_predicted": top.index[0], "top_share": float(top.iloc[0])})
    return pd.DataFrame(rows)


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

        # 2. stratified 5-fold OOF of the picked default, with novelty (cached for the ko/shift stages)
        oof = cross_validate(a, method=method, n_folds=5, seed=seed, n_jobs=nj)
        oof["section"] = pd.Series(sec, index=a.obs_names).loc[oof.index].to_numpy()
        oof.to_pickle(out / "cells" / f"oof_{name}.pkl")
        strat = {"stage": name, "split": "stratified 5-fold CV", "method": method, **transfer_metrics(oof),
                 **fold_sd(oof), "novelty_gt1": float((oof["novelty"] > 1).mean())}
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
        pred.insert(0, "true_group", a.obs["group"].astype(str).to_numpy()[te])
        pred.insert(0, "true", labels[te])
        pred["section"] = sec[te]
        met = transfer_metrics(pred)
        summary.append({"stage": name, "split": "stratified 20% hold-out", "method": model.method, **met,
                        "novelty_gt1": float((pred["novelty"] > 1).mean()), "n_train": len(tr)})
        log(f"  hold-out ({model.method}): acc {met['accuracy']:.4f} macroF1 {met['macro_f1']:.4f} "
            f"group acc {met['group_accuracy']:.4f} group macroF1 {met['group_macro_f1']:.4f}")
        pred.to_pickle(out / "cells" / f"holdout_{name}.pkl")
        with open(out / "cells" / f"model_holdout_{name}.pkl", "wb") as fh:
            pickle.dump(model, fh)

        rep = classification_report(pred["true"], pred["label_raw"], labels=sorted(set(pred["true"])),
                                    output_dict=True, zero_division=0)
        per = pd.DataFrame({k: v for k, v in rep.items() if k in set(pred["true"])}).T
        per.index.name = "label"
        per["group"] = [lineage_map()["labels"].get(lab, "Unknown") for lab in per.index]
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
# Stage ko


def _weights(pred: pd.DataFrame, kind: str, level: str, cats: list[str]) -> pd.DataFrame:
    """Per-cell weights for a composition: top calls ("raw"), thresholded calls with Unassigned ("hard") or
    probabilities ("soft")."""
    if kind == "raw":
        return one_hot(pred[f"{level}_raw"].astype(str).to_numpy(), cats + [UNASSIGNED])
    if kind == "hard":
        return one_hot(pred[level].astype(str).to_numpy(), cats + [UNASSIGNED])
    prefix = "pg:" if level == "group" else "p:"
    W = pred[[f"{prefix}{c}" for c in cats]].copy()
    W.columns = cats
    return W.reset_index(drop=True)


def _comp_block(named: dict, n_boot: int, seed: int) -> pd.DataFrame:
    cols = {}
    for name, (W, cl) in named.items():
        b = bootstrap_composition(W, cl, n_boot=n_boot, seed=seed)
        cols[name], cols[f"{name}_lo"], cols[f"{name}_hi"] = b["estimate"], b["lo"], b["hi"]
    return pd.DataFrame(cols)


KO_DEPTH = {"match": dict(depth_match=True, depth_draws=1), "match5": dict(depth_match=True, depth_draws=5),
            "none": dict(depth_match=False)}


def em_draws(P, priors, clusters, n_boot: int, seed: int):
    """EM priors on all cells and on cluster-bootstrap resamples (weighted, warm-started EM)."""
    q0, _ = em_priors(P, priors)
    ids, inv = np.unique(np.asarray(clusters), return_inverse=True)
    mult = np.random.default_rng(seed).multinomial(len(ids), np.full(len(ids), 1.0 / len(ids)), size=n_boot)
    D = np.empty((n_boot, P.shape[1]))
    for b in range(n_boot):
        w = mult[b][inv].astype(float)
        keep = w > 0
        D[b] = em_priors(P[keep], priors, weights=w[keep], init=q0, tol=1e-7)[0]
    return q0, D


def stage_ko(out: Path, nj: int, seed: int, n_boot: int, ko_depth: str = "match5") -> None:
    wt, ko = load("E9.5"), load("KO")
    d_wt, d_ko = median_depth(wt), median_depth(ko)
    log(f"WT E9.5 {wt.n_obs} cells (depth {d_wt:.0f}); KO {ko.n_obs} cells (depth {d_ko:.0f}); "
        f"same genes in same order: {list(wt.var_names) == list(ko.var_names)}")
    model = fit_transfer(wt, "celltype", classifier="auto", cv_folds=5, seed=seed, n_jobs=nj)
    model.cv.to_csv(out / "ko_model_cv.csv", index=False)
    with open(out / "cells" / "model_E9.5_full.pkl", "wb") as fh:
        pickle.dump(model, fh)
    log(f"WT E9.5 model: {model.method}, {len(model.classes)} labels, dropped {model.dropped}")

    variants = {k: predict(model, ko, seed=seed, **kw) for k, kw in KO_DEPTH.items()}
    pred = variants[ko_depth]
    pred_nodm = variants["none"]
    log(f"KO primary variant '{ko_depth}': kept {pred.attrs['depth_frac']:.3f} of its molecules "
        f"(target depth {model.depth:.0f}, {pred.attrs['depth_draws']} draw(s))")
    pred.assign(ko_label=ko.obs["celltype"].astype(str).to_numpy()).to_pickle(out / "cells" / "ko_pred.pkl")

    oof = pd.read_pickle(out / "cells" / "oof_E9.5.pkl") if (out / "cells" / "oof_E9.5.pkl").exists() else None
    if oof is None or oof.attrs.get("method") != model.method:
        oof = cross_validate(wt, method=model.method, n_folds=5, seed=seed, n_jobs=nj)
        oof["section"] = pd.Series(section_ids(wt), index=wt.obs_names).loc[oof.index].to_numpy()
    sec_wt_all = section_ids(wt)
    sec_wt = oof["section"].to_numpy()
    sec_ko = section_ids(ko)
    log(f"sections: WT {len(np.unique(sec_wt_all))}, KO {len(np.unique(sec_ko))}; WT OOF cells {len(oof)}")

    labels = list(model.classes)
    groups = list(model.groups)
    wt_label = wt.obs["celltype"].astype(str).to_numpy()
    wt_group = wt.obs["group"].astype(str).to_numpy()
    summary = {"wt_depth": d_wt, "ko_depth": d_ko, "ko_depth_frac": variants["match"].attrs["depth_frac"],
               "ko_primary_variant": ko_depth,
               "method": model.method, "n_boot": n_boot, "wt_sections": int(len(np.unique(sec_wt_all))),
               "ko_sections": int(len(np.unique(sec_ko)))}

    n_boot_em = min(n_boot, 1000)
    t_em = time.time()
    em_wt, D_wt = em_draws(oof[[f"p:{c}" for c in labels]].to_numpy(), model.priors, sec_wt, n_boot_em, seed)
    em_ko, D_ko = em_draws(pred[[f"p:{c}" for c in labels]].to_numpy(), model.priors, sec_ko, n_boot_em, seed + 1)
    summary["n_boot_em"] = n_boot_em
    log(f"EM cluster bootstrap ({n_boot_em} draws per embryo) in {time.time() - t_em:.0f}s")

    for level, cats, wt_true in (("label", labels, wt_label), ("group", groups, wt_group)):
        named = {
            "wt_annotated": (one_hot(wt_true, cats + [UNASSIGNED]), sec_wt_all),
            "wt_pred_raw": (_weights(oof, "raw", level, cats), sec_wt),
            "ko_pred_raw": (_weights(pred, "raw", level, cats), sec_ko),
            "ko_pred_raw_no_depth_match": (_weights(pred_nodm, "raw", level, cats), sec_ko),
            "wt_pred_hard": (_weights(oof, "hard", level, cats), sec_wt),
            "wt_pred_soft": (_weights(oof, "soft", level, cats).reindex(columns=cats + [UNASSIGNED], fill_value=0.0),
                             sec_wt),
            "ko_pred_hard": (_weights(pred, "hard", level, cats), sec_ko),
            "ko_pred_soft": (_weights(pred, "soft", level, cats).reindex(columns=cats + [UNASSIGNED],
                                                                          fill_value=0.0), sec_ko),
            "ko_pred_soft_no_depth_match": (_weights(pred_nodm, "soft", level, cats).reindex(
                columns=cats + [UNASSIGNED], fill_value=0.0), sec_ko),
        }
        table = _comp_block(named, n_boot, seed)
        f = (lambda v: v @ model.membership) if level == "group" else (lambda v: v)
        ew, ek, Dw, Dk = f(em_wt), f(em_ko), f(D_wt), f(D_ko)
        lr = np.log2((Dk + 1e-3) / (Dw + 1e-3))
        em = pd.DataFrame({
            "wt_pred_em": ew, "wt_pred_em_lo": np.quantile(Dw, 0.025, 0), "wt_pred_em_hi": np.quantile(Dw, 0.975, 0),
            "ko_pred_em": ek, "ko_pred_em_lo": np.quantile(Dk, 0.025, 0), "ko_pred_em_hi": np.quantile(Dk, 0.975, 0),
            "em_diff": ek - ew, "em_diff_lo": np.quantile(Dk - Dw, 0.025, 0),
            "em_diff_hi": np.quantile(Dk - Dw, 0.975, 0), "em_log2_ratio": np.log2((ek + 1e-3) / (ew + 1e-3)),
            "em_log2_lo": np.quantile(lr, 0.025, 0), "em_log2_hi": np.quantile(lr, 0.975, 0),
        }, index=cats)
        np.testing.assert_allclose(ek, composition(pred, level, "em", model).reindex(cats).to_numpy(), atol=1e-6)
        table = table.join(em).fillna({c: 0.0 for c in em.columns})
        table["em_ci_excludes_0"] = (table["em_diff_lo"] > 0) | (table["em_diff_hi"] < 0)
        for kind in ("raw", "hard", "soft"):
            d = bootstrap_difference(named[f"wt_pred_{kind}"][0], sec_wt, named[f"ko_pred_{kind}"][0], sec_ko,
                                     n_boot=n_boot, seed=seed)
            for c in ("diff", "diff_lo", "diff_hi", "log2_ratio", "log2_lo", "log2_hi", "ci_excludes_0"):
                table[f"{kind}_{c}"] = d[c]
        d = bootstrap_difference(named["wt_annotated"][0], sec_wt_all, named["ko_pred_soft"][0], sec_ko,
                                 n_boot=n_boot, seed=seed)
        table["soft_vs_annotated_diff"], table["soft_vs_annotated_lo"], table["soft_vs_annotated_hi"] = (
            d["diff"], d["diff_lo"], d["diff_hi"])
        table.index.name = level
        table.to_csv(out / f"ko_composition_{level}s.csv")
        log(f"composition by {level}: 95% CI of the KO-WT difference excludes 0 for "
            f"{int(table['raw_ci_excludes_0'].sum())} (top call) / {int(table['soft_ci_excludes_0'].sum())} (soft) "
            f"/ {int(table['em_ci_excludes_0'].sum())} (EM) of {len(cats)} categories")

    # within-cardiomyocyte composition (CM-called cells only), less sensitive to dissection scope
    cm_labels = [c for c in labels if model.group_of[c] == "CM"]
    # top call: cells whose top fine label is a CM label, split by that label
    wt_cm = np.isin(oof["label_raw"].to_numpy(), cm_labels)
    ko_cm = np.isin(pred["label_raw"].to_numpy(), cm_labels)
    d_raw = bootstrap_difference(one_hot(oof["label_raw"].to_numpy()[wt_cm], cm_labels), sec_wt[wt_cm],
                                 one_hot(pred["label_raw"].to_numpy()[ko_cm], cm_labels), sec_ko[ko_cm],
                                 n_boot=n_boot, seed=seed)
    # soft: CM-group cells, CM-label probabilities renormalised within CM
    wt_cm = oof["group_raw"].to_numpy() == "CM"
    ko_cm = pred["group_raw"].to_numpy() == "CM"
    W_wt = oof.loc[wt_cm, [f"p:{c}" for c in cm_labels]].to_numpy()
    W_ko = pred.loc[ko_cm, [f"p:{c}" for c in cm_labels]].to_numpy()
    W_wt = pd.DataFrame(W_wt / W_wt.sum(1, keepdims=True), columns=cm_labels)
    W_ko = pd.DataFrame(W_ko / W_ko.sum(1, keepdims=True), columns=cm_labels)
    d_soft = bootstrap_difference(W_wt, sec_wt[wt_cm], W_ko, sec_ko[ko_cm], n_boot=n_boot, seed=seed)
    ann = wt_group == "CM"
    within = pd.concat({"raw": d_raw.rename(columns={"a": "wt", "b": "ko"}),
                        "soft": d_soft.rename(columns={"a": "wt", "b": "ko"})}, axis=1)
    within.columns = [f"{k}_{c}" for k, c in within.columns]
    within.insert(0, "wt_annotated", bootstrap_composition(one_hot(wt_label[ann], cm_labels), sec_wt_all[ann],
                                                           n_boot=n_boot, seed=seed)["estimate"])
    within.index.name = "cm_label"
    within.to_csv(out / "ko_composition_within_cm.csv")

    # cross-tab against the KO's own labels
    ko_label = ko.obs["celltype"].astype(str).to_numpy()
    ko_group = ko.obs["group"].astype(str).to_numpy()
    ct = pd.crosstab(pd.Series(ko_label, name="ko_label"), pd.Series(pred["label"].to_numpy(), name="predicted"))
    ct = ct[[c for c in labels + [UNASSIGNED] if c in ct.columns]]
    ct.to_csv(out / "ko_crosstab_counts.csv")
    ct.div(ct.sum(1), axis=0).to_csv(out / "ko_crosstab_rowfrac.csv")
    gct = group_confusion(ko_group, pred["group_raw"].to_numpy())
    gct.to_csv(out / "ko_group_confusion.csv")

    rows = []
    shared_groups = set(groups)
    for lab in pd.unique(ko_label):
        m = ko_label == lab
        fine = pd.Series(pred["label_raw"].to_numpy()[m]).value_counts(normalize=True)
        grp = pd.Series(pred["group_raw"].to_numpy()[m]).value_counts(normalize=True)
        own = lineage_map()["labels"].get(lab, "Unknown")
        rows.append({
            "ko_label": lab, "n": int(m.sum()), "lineage_group": own, "group_in_wt_vocab": own in shared_groups,
            "top_wt_label": fine.index[0], "top_wt_label_share": float(fine.iloc[0]),
            "second_wt_label": fine.index[1] if len(fine) > 1 else "", "second_share": float(fine.iloc[1])
            if len(fine) > 1 else 0.0,
            "top_group": grp.index[0], "top_group_share": float(grp.iloc[0]),
            "agrees_with_lineage_map": float((pred["group_raw"].to_numpy()[m] == own).mean()),
            "unassigned_at_0.5": float((pred["label"].to_numpy()[m] == UNASSIGNED).mean()),
            "novelty_gt1": float((pred["novelty"].to_numpy()[m] > 1).mean()),
            "median_conf": float(np.median(pred["label_conf"].to_numpy()[m])),
        })
    conc = pd.DataFrame(rows).sort_values("n", ascending=False)
    conc.to_csv(out / "ko_label_concordance.csv", index=False)
    comparable = np.isin(ko_group, list(shared_groups))
    for tag, p in variants.items():
        g_hat = p["group_raw"].to_numpy()
        summary[f"group_agreement_{tag}"] = float((g_hat[comparable] == ko_group[comparable]).mean())
        summary[f"unassigned_{tag}"] = float((p["label"] == UNASSIGNED).mean())
        summary[f"novelty_gt1_{tag}"] = float((p["novelty"] > 1).mean())
        summary[f"median_label_conf_{tag}"] = float(p["label_conf"].median())
        summary[f"label_change_vs_primary_{tag}"] = float((p["label_raw"] != pred["label_raw"]).mean())
    summary["comparable_fraction"] = float(comparable.mean())
    log("KO own-group agreement on comparable cells: " + ", ".join(
        f"{k} {summary[f'group_agreement_{k}']:.3f}" for k in variants))

    # reverse mapping: WT cells into the KO vocabulary (KO thinned to WT depth), and mutual best matches
    rev = fit_transfer(ko, "celltype", classifier=model.method, depth=d_wt, seed=seed, n_jobs=nj)
    rpred = predict(rev, wt)
    rct = pd.crosstab(pd.Series(wt_label, name="wt_label"), pd.Series(rpred["label_raw"].to_numpy(),
                                                                       name="predicted_ko_label"))
    rct.to_csv(out / "ko_reverse_crosstab_counts.csv")
    rfrac = rct.div(rct.sum(1), axis=0)
    rfrac.to_csv(out / "ko_reverse_crosstab_rowfrac.csv")
    fct = pd.crosstab(pd.Series(ko_label, name="ko_label"), pd.Series(pred["label_raw"].to_numpy(), name="wt"))
    ffrac = fct.div(fct.sum(1), axis=0)
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
    (out / "ko_summary.json").write_text(json.dumps(summary, indent=2))


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
# Stage cross


def stage_cross(out: Path, nj: int, seed: int) -> None:
    e875, e95 = load("E8.75"), load("E9.5")
    d875, d95 = median_depth(e875), median_depth(e95)
    log(f"E8.75 depth {d875:.0f}, E9.5 depth {d95:.0f}")
    mA = fit_transfer(e875, "celltype", classifier="auto", cv_folds=5, seed=seed, n_jobs=nj)
    mA_thin = fit_transfer(e875, "celltype", classifier=mA.method, depth=d95, seed=seed, n_jobs=nj)
    mB = fit_transfer(e95, "celltype", classifier="auto", cv_folds=5, seed=seed, n_jobs=nj)
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
        true = q.obs["group"].astype(str).to_numpy()
        true_lab = q.obs["celltype"].astype(str).to_numpy()
        known = true != "Unknown"
        transferable = known & np.isin(true, m.groups)
        g_hat = p["group_raw"].to_numpy()
        rows.append({
            "direction": tag, "method": m.method, "depth_frac": p.attrs["depth_frac"],
            "depth_draws": p.attrs["depth_draws"], "reference_depth": m.depth,
            "n_query": int(q.n_obs), "known_group_cells": int(known.sum()),
            "transferable_fraction": float(transferable.sum() / known.sum()),
            "group_accuracy_transferable": float((g_hat[transferable] == true[transferable]).mean()),
            "group_macro_f1_transferable": float(f1_score(true[transferable], g_hat[transferable],
                                                          labels=np.unique(true[transferable]), average="macro",
                                                          zero_division=0)),
            "group_accuracy_all_known": float((g_hat[known] == true[known]).mean()),
            "unassigned_label_0.5": float((p["label"] == UNASSIGNED).mean()),
            "unassigned_group_0.5": float((p["group"] == UNASSIGNED).mean()),
            "novelty_gt1": float((p["novelty"] > 1).mean()),
        })
        pg = per_group_scores(true[known], g_hat[known], m.groups)
        sub = pd.DataFrame({"g": true[known], "nov": p["novelty"].to_numpy()[known],
                            "un": p["group"].to_numpy()[known] == UNASSIGNED})
        pg["novelty_gt1"] = pg["group"].map(sub.groupby("g")["nov"].apply(lambda v: float((v > 1).mean())))
        pg["unassigned_group_0.5"] = pg["group"].map(sub.groupby("g")["un"].mean())
        pg.insert(0, "direction", tag)
        pergroup.append(pg)
        safe = tag.replace("->", "_to_").replace("(", "_").replace(")", "").replace(" ", "_")
        keep_cols = ["label_raw", "group_raw", "label_conf", "group_conf", "novelty"] + [
            c for c in p.columns if c.startswith("p:")]
        p[keep_cols].assign(true_label=true_lab, true_group=true, section=section_ids(q)).to_pickle(
            out / "cells" / f"cross_{safe}.pkl")
        group_confusion(true, g_hat).to_csv(out / f"cross_group_confusion_{safe}.csv")
        ct = pd.crosstab(pd.Series(true_lab, name="query_label"), pd.Series(p["label_raw"].to_numpy(),
                                                                            name="reference_label"))
        ct.div(ct.sum(1), axis=0).to_csv(out / f"cross_crosstab_rowfrac_{safe}.csv")
        log(f"{tag}: transferable {rows[-1]['transferable_fraction']:.3f}, group acc "
            f"{rows[-1]['group_accuracy_transferable']:.4f}, macroF1 {rows[-1]['group_macro_f1_transferable']:.4f}")
    pd.DataFrame(rows).to_csv(out / "cross_summary.csv", index=False)
    pg = pd.concat(pergroup, ignore_index=True)
    pg.to_csv(out / "cross_pergroup.csv", index=False)

    a = pg[pg["direction"] == "E8.75(thinned)->E9.5"].set_index("group")
    b = pg[pg["direction"] == "E9.5->E8.75(5 thinning draws)"].set_index("group")
    allg = group_order(set(a.index) | set(b.index))
    verdict = []
    for g in allg:
        in875 = g in set(mA.groups)
        in95 = g in set(mB.groups)
        fa = a["f1"].get(g, np.nan) if in875 else np.nan
        fb = b["f1"].get(g, np.nan) if in95 else np.nan
        if not (in875 and in95):
            v = "one stage only"
        elif min(fa, fb) >= 0.8:
            v = "clean"
        elif min(fa, fb) >= 0.5:
            v = "partial"
        else:
            v = "poor"
        verdict.append({"group": g, "labelled_in_E8.75": in875, "labelled_in_E9.5": in95,
                        "f1_E8.75_to_E9.5": fa, "recall_E8.75_to_E9.5": a["recall"].get(g, np.nan),
                        "precision_E8.75_to_E9.5": a["precision"].get(g, np.nan),
                        "f1_E9.5_to_E8.75": fb, "recall_E9.5_to_E8.75": b["recall"].get(g, np.nan),
                        "precision_E9.5_to_E8.75": b["precision"].get(g, np.nan), "verdict": v})
    pd.DataFrame(verdict).to_csv(out / "cross_verdict.csv", index=False)


# ---------------------------------------------------------------------------------------------------------------------
# Stage contrast


def _contrast_rows(setting, pred, model, true_label, true_group, sections, levels):
    """Composition change between the lower and upper halves of an embryo's sections: estimated vs true."""
    z = np.array([float(sec[1:]) for sec in sections])
    upper = z > np.median(np.unique(z))
    classes, groups, M = list(model.classes), list(model.groups), model.membership
    P = pred[[f"p:{c}" for c in classes]].to_numpy()
    raw = np.searchsorted(model.classes, pred["label_raw"].to_numpy())
    rows = []
    for level in levels:
        if level == "label":
            ti, K, f = np.searchsorted(model.classes, true_label), len(classes), (lambda v: v)
        else:
            gi = {g: i for i, g in enumerate(groups)}
            ti, K, f = np.array([gi[g] for g in true_group]), len(groups), (lambda v: v @ M)

        def comps(sel, ti=ti, K=K, f=f):
            n = sel.sum()
            return {"true": np.bincount(ti[sel], minlength=K) / n,
                    "top call": f(np.bincount(raw[sel], minlength=len(classes)) / n),
                    "mean probability": f(P[sel].mean(0)),
                    "EM": f(em_priors(P[sel], model.priors)[0])}

        A, B = comps(~upper), comps(upper)
        dT = B["true"] - A["true"]
        for m in ("top call", "mean probability", "EM"):
            dE = B[m] - A[m]
            rows.append({"setting": setting, "level": level, "estimator": m, "n_lower": int((~upper).sum()),
                         "n_upper": int(upper.sum()), "true_contrast_tv": float(np.abs(dT).sum() / 2),
                         "delta_slope": float(dE @ dT / (dT @ dT)), "delta_mae": float(np.abs(dE - dT).mean()),
                         "delta_max_err": float(np.abs(dE - dT).max()),
                         "composition_mae": float((np.abs(A[m] - A["true"]).mean()
                                                   + np.abs(B[m] - B["true"]).mean()) / 2)})
    return rows


def stage_contrast(out: Path) -> None:
    rows = []
    for st in ("E9.5", "E8.75"):
        pred = pd.read_pickle(out / "cells" / f"holdout_{st}.pkl")
        with open(out / "cells" / f"model_holdout_{st}.pkl", "rb") as fh:
            model = pickle.load(fh)
        ok = np.isin(pred["true"], model.classes)
        pred = pred[ok]
        rows += _contrast_rows(f"within {st}: 20% hold-out cells", pred, model, pred["true"].to_numpy(),
                               pred["true_group"].to_numpy(), pred["section"].to_numpy(), ("label", "group"))
    with open(out / "cells" / "models_cross.pkl", "rb") as fh:
        models = pickle.load(fh)
    for tag, key in (("E8.75_thinned_to_E9.5", "E8.75_thinned"), ("E9.5_to_E8.75_5_thinning_draws", "E9.5")):
        p = pd.read_pickle(out / "cells" / f"cross_{tag}.pkl")
        m = models[key]
        p = p[np.isin(p["true_group"], m.groups)]  # oracle restriction to groups the reference has
        rows += _contrast_rows(f"cross-stage {tag}: transferable cells", p, m, None, p["true_group"].to_numpy(),
                               p["section"].to_numpy(), ("group",))
    df = pd.DataFrame(rows)
    df.to_csv(out / "contrast_estimators.csv", index=False)
    log("contrast: " + "; ".join(f"{r.setting.split(':')[0]}/{r.level}/{r.estimator}: slope {r.delta_slope:.3f} "
                                 f"mae {r.delta_mae:.4f}" for r in df.itertuples()))


# ---------------------------------------------------------------------------------------------------------------------
# Figures

INK, INK2, MUTED, GRID, BASE, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"
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


def stage_figs(out: Path, estimator: str = "raw") -> None:
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

    # 2. KO vs WT group composition with section-bootstrap intervals
    p = out / "ko_composition_groups.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        e = estimator
        t = t.loc[(t[[f"wt_pred_{e}", f"ko_pred_{e}", "wt_annotated"]].max(1) > 0.002)]
        t = t.drop(index=[i for i in t.index if i == UNASSIGNED])
        t = t.sort_values(f"wt_pred_{e}")
        y = np.arange(len(t))
        fig, ax = plt.subplots(figsize=(7.5, 0.42 * len(t) + 1.4), constrained_layout=True)
        ax.grid(axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for col, color, dy, lab in (
            (f"wt_pred_{e}", BLUE, 0.14, "WT E9.5 (classifier, out-of-fold)"),
            (f"ko_pred_{e}", ORANGE, -0.14, "Mab21l2 KO (classifier)"),
        ):
            lo, hi = f"{col}_lo", f"{col}_hi"
            ax.hlines(y + dy, t[lo], t[hi], color=color, linewidth=2)
            ax.scatter(t[col], y + dy, s=36, color=color, edgecolor=SURFACE, linewidth=1.5, zorder=3, label=lab)
        ax.scatter(t["wt_annotated"], y + 0.14, s=34, facecolor="none", edgecolor=INK2, linewidth=1.0, zorder=4,
                   label="WT E9.5 (its own labels)")
        ax.set_yticks(y, t.index)
        how = {"raw": "top call per cell", "soft": "mean class probability", "em": "EM prior re-estimation"}[e]
        ax.set_xlabel(f"fraction of cells ({how}; bars: 95% section bootstrap)")
        ax.set_title("Lineage-group composition in the WT E9.5 vocabulary", loc="left")
        ax.legend(frameon=False, loc="lower right", fontsize=8)
        fig.savefig(out / "fig_ko_vs_wt_groups.png", dpi=160)
        plt.close(fig)

    # 3. KO own labels vs predicted WT labels
    p = out / "ko_crosstab_rowfrac.csv"
    if p.exists():
        ct = pd.read_csv(p, index_col=0)
        counts = pd.read_csv(out / "ko_crosstab_counts.csv", index_col=0)
        ct = ct.loc[counts.sum(1).sort_values(ascending=False).index]
        fig, ax = plt.subplots(figsize=(11, 6.2), constrained_layout=True)
        im = _heatmap(ax, ct, "Mab21l2 KO: own labels (rows) mapped to WT E9.5 labels (threshold 0.5)",
                      "predicted WT E9.5 label", "KO file label")
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
    log("figures written")


# ---------------------------------------------------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stages", default="ref,ko,shift,cross,contrast,figs")
    ap.add_argument("--out", default="runs/harmonise")
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-cells", type=int, default=None, help="subsample every file (smoke tests only)")
    ap.add_argument("--ko-depth", choices=sorted(KO_DEPTH), default="match5",
                    help="depth handling for the KO predictions used in the composition tables")
    ap.add_argument("--fig-estimator", choices=["raw", "soft", "em"], default="raw",
                    help="composition estimator plotted in fig_ko_vs_wt_groups")
    args = ap.parse_args()
    global MAX_CELLS
    MAX_CELLS = args.max_cells
    out = Path(args.out)
    (out / "cells").mkdir(parents=True, exist_ok=True)
    manifest = {
        "started": time.strftime("%Y-%m-%d %H:%M:%S"), "stages": args.stages, "seed": args.seed,
        "n_boot": args.n_boot, "max_cells": args.max_cells, "ko_depth": args.ko_depth,
        "fig_estimator": args.fig_estimator,
        "python": platform.python_version(),
        "sklearn": sklearn.__version__,
        "harmonise_py_sha256_16": sha256("src/vec/harmonise.py"), "script_sha256_16": sha256(__file__),
        "data_sha256_16": {k: sha256(v) for k, v in FILES.items()},
    }
    for stage in args.stages.split(","):
        log(f"=== stage {stage}")
        if stage == "ref":
            stage_ref(out, args.n_jobs, args.seed)
        elif stage == "ko":
            stage_ko(out, args.n_jobs, args.seed, args.n_boot, args.ko_depth)
        elif stage == "shift":
            stage_shift(out, args.seed)
        elif stage == "cross":
            stage_cross(out, args.n_jobs, args.seed)
        elif stage == "contrast":
            stage_contrast(out)
        elif stage == "figs":
            stage_figs(out, args.fig_estimator)
        else:
            raise SystemExit(f"unknown stage {stage}")
    manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    manifest["seconds"] = round(time.time() - T0, 1)
    (out / f"manifest_{args.stages.replace(',', '-')}.json").write_text(json.dumps(manifest, indent=2))
    log("done")


if __name__ == "__main__":
    main()
