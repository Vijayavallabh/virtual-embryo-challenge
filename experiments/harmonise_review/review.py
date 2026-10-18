"""Adversarial review of experiments/harmonise. Run from the repo root on the DGX:

    OMP_NUM_THREADS=16 uv run python runs/harmonise_review/review.py --parts A,B,C,D,E,F,G,H

A  independent re-implementation of the 20% hold-out (sklearn directly on .X), 3 split seeds
B  independent KO-vs-WT composition (own OOF, own thinning, own bootstrap), LR and kNN, with/without depth match
C  where along z each group sits in the WT and KO stacks (coverage confound)
D  section-bootstrap calibration (odd vs even sections), trimming/jackknife sensitivity of KO-WT differences
E  composition estimators when the query is shifted: train on one z-half, OOF on it, predict the other half
F  seed robustness of the cross-stage depth claims (vec.harmonise API)
G  within-CM composition with a CM-only classifier
H  gene-level technical shift between the KO and WT files in matched non-Mesp1 populations, and a corrected re-call
Writes small CSVs to runs/harmonise_review/. Reads only the three released files plus runs/harmonise/cells (DGX).
"""
from __future__ import annotations

import argparse
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler

from vec import harmonise as H
from vec.prep import lineage_map, load_stage

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 60)
pd.set_option("display.max_rows", 200)
OUT = Path("runs/harmonise_review")
HR = Path("runs/harmonise")
FILES = {"E9.5": "data/raw/T2/E9.5.h5ad", "E8.75": "data/raw/T2/E8.75.h5ad", "KO": "data/raw/T3/E9.5_mab21l2_ko.h5ad"}
GM = lineage_map()["labels"]
GORDER = list(lineage_map()["groups"])
CM_LABELS = ["A-CM", "OFT-CM", "SV-CM", "V-CM"]
T0 = time.time()
_cache: dict = {}
STATE: dict = {}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def load(name):
    if name not in _cache:
        _cache[name] = load_stage(FILES[name])
    return _cache[name]


def clip10(X):
    return np.clip(X, -10, 10)


def pipe(method="logreg", n_pcs=50, C=1.0, k=15, seed=0):
    steps = [("scale", StandardScaler()), ("clip", FunctionTransformer(clip10)),
             ("pca", PCA(n_pcs, random_state=seed))]
    if method == "logreg":
        steps += [("pcs", StandardScaler()), ("clf", LogisticRegression(C=C, max_iter=2000))]
    else:
        steps += [("clf", KNeighborsClassifier(k, n_jobs=16))]
    return Pipeline(steps)


def lognorm(counts):
    C = sparse.csr_matrix(counts, dtype=np.float64)
    lib = np.asarray(C.sum(1)).ravel()
    C = sparse.diags(1e4 / np.maximum(lib, 1)) @ C
    C.data = np.log1p(C.data)
    return C.toarray().astype(np.float32)


def thin(counts, frac, rng):
    C = sparse.csr_matrix(counts, dtype=np.int64, copy=True)
    C.data = rng.binomial(C.data, frac)
    C.eliminate_zeros()
    return C


def group_matrix(classes):
    g = [GM[c] for c in classes]
    order = [x for x in GORDER if x in set(g)]
    return order, np.array([[gi == o for o in order] for gi in g], float)


def oof_proba(X, y, method="logreg", n_folds=5, seed=7):
    cls = np.unique(y)
    P = np.zeros((len(y), len(cls)))
    for tr, te in StratifiedKFold(n_folds, shuffle=True, random_state=seed).split(X, y):
        p = pipe(method).fit(X[tr], y[tr])
        P[np.ix_(te, np.searchsorted(cls, p.classes_))] = p.predict_proba(X[te])
    return P, cls


def my_boot(W, cl, n_boot=1000, seed=0):
    """Cluster bootstrap of a composition (own implementation; resample cluster indices with replacement)."""
    W = np.asarray(W, float)
    ids, inv = np.unique(np.asarray(cl), return_inverse=True)
    S = np.zeros((len(ids), W.shape[1]))
    np.add.at(S, inv, W)
    n = np.bincount(inv).astype(float)
    draws = np.random.default_rng(seed).integers(0, len(ids), size=(n_boot, len(ids)))
    return np.stack([S[d].sum(0) / n[d].sum() for d in draws])


def onehot(idx, K):
    W = np.zeros((len(idx), K))
    W[np.arange(len(idx)), idx] = 1.0
    return W


def diff_table(Wa, cla, Wb, clb, names, n_boot=1000, seed=0):
    Da, Db = my_boot(Wa, cla, n_boot, seed), my_boot(Wb, clb, n_boot, seed + 1)
    d = Db - Da
    t = pd.DataFrame({"a": np.asarray(Wa).mean(0), "b": np.asarray(Wb).mean(0)}, index=names)
    t["diff"] = t["b"] - t["a"]
    t["lo"], t["hi"] = np.quantile(d, 0.025, 0), np.quantile(d, 0.975, 0)
    t["excl0"] = (t["lo"] > 0) | (t["hi"] < 0)
    return t


# ---------------------------------------------------------------------------------------------------------------------


def part_a():
    rows = []
    for name in ("E9.5", "E8.75"):
        a = load(name)
        X = a.X.toarray()
        lab = a.obs["celltype"].astype(str).to_numpy()
        vc = pd.Series(lab).value_counts()
        idx = np.flatnonzero((lab != "Unknown") & np.isin(lab, vc.index[vc >= 10]))
        for s in (0, 1, 2):
            tr, te = train_test_split(idx, test_size=0.2, stratify=lab[idx], random_state=s)
            p = pipe().fit(X[tr], lab[tr])
            P = p.predict_proba(X[te])
            cls = p.classes_
            order, M = group_matrix(cls)
            y = lab[te]
            yhat = cls[P.argmax(1)]
            ghat = np.array(order)[(P @ M).argmax(1)]
            gtrue = np.array([GM[x] for x in y])
            rows.append(dict(stage=name, split_seed=s, n_train=len(tr), n_test=len(te),
                             accuracy=accuracy_score(y, yhat),
                             macro_f1=f1_score(y, yhat, labels=np.unique(y), average="macro", zero_division=0),
                             group_accuracy=accuracy_score(gtrue, ghat),
                             group_macro_f1=f1_score(gtrue, ghat, labels=np.unique(gtrue), average="macro",
                                                     zero_division=0),
                             coverage_0p5=float((P.max(1) >= 0.5).mean())))
            log(f"A {name} seed {s}: " + ", ".join(f"{k} {v:.4f}" for k, v in rows[-1].items()
                                                  if isinstance(v, float)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "A_holdout_independent.csv", index=False)
    print(df.round(4).to_string())


def part_b():
    wt, ko = load("E9.5"), load("KO")
    Xw = wt.X.toarray()
    lw = wt.obs["celltype"].astype(str).to_numpy()
    sec_w, sec_k = H.section_ids(wt), H.section_ids(ko)
    dw = float(np.median(np.asarray(wt.layers["counts"].sum(1)).ravel()))
    Ck = sparse.csr_matrix(ko.layers["counts"])
    dk = float(np.median(np.asarray(Ck.sum(1)).ravel()))
    frac = dw / dk
    rng = np.random.default_rng(12345)
    Xk_draws = [lognorm(thin(Ck, frac, rng)) for _ in range(5)]
    Xk_raw = ko.X.toarray()
    log(f"B: WT depth {dw}, KO depth {dk}, frac {frac:.4f}; max |lognorm(KO counts) - KO.X| = "
        f"{np.abs(lognorm(Ck) - Xk_raw).max():.2e}")
    STATE.update(Xw=Xw, lw=lw, sec_w=sec_w, sec_k=sec_k, Xk_draws=Xk_draws, Xk_raw=Xk_raw, dw=dw, dk=dk)
    group_rows, label_rows, cm_rows = [], [], []
    for method in ("logreg", "knn"):
        P_oof, cls = oof_proba(Xw, lw, method)
        full = pipe(method).fit(Xw, lw)
        assert np.array_equal(full.classes_, cls)
        Pk = np.mean([full.predict_proba(X) for X in Xk_draws], axis=0)
        Pk_nodm = full.predict_proba(Xk_raw)
        order, M = group_matrix(cls)
        STATE[method] = dict(P_oof=P_oof, Pk=Pk, Pk_nodm=Pk_nodm, cls=cls, order=order, M=M, model=full)
        for tag, PK in (("match5", Pk), ("none", Pk_nodm)):
            gw, gk = (P_oof @ M).argmax(1), (PK @ M).argmax(1)
            t = diff_table(onehot(gw, len(order)), sec_w, onehot(gk, len(order)), sec_k, order)
            t.insert(0, "variant", f"{method}/{tag}")
            group_rows.append(t)
            lwi, lki = P_oof.argmax(1), PK.argmax(1)
            tl = diff_table(onehot(lwi, len(cls)), sec_w, onehot(lki, len(cls)), sec_k, list(cls))
            tl.insert(0, "variant", f"{method}/{tag}")
            label_rows.append(tl)
            cmi = [list(cls).index(c) for c in CM_LABELS]
            mw, mk = np.isin(lwi, cmi), np.isin(lki, cmi)
            Ww = np.stack([lwi[mw] == i for i in cmi], 1)
            Wk = np.stack([lki[mk] == i for i in cmi], 1)
            tc = diff_table(Ww.astype(float), sec_w[mw], Wk.astype(float), sec_k[mk], CM_LABELS)
            tc.insert(0, "variant", f"{method}/{tag}")
            cm_rows.append(tc)
        log(f"B {method} done")
    G = pd.concat(group_rows)
    L = pd.concat(label_rows)
    C = pd.concat(cm_rows)
    G.to_csv(OUT / "B_ko_groups_independent.csv")
    L.to_csv(OUT / "B_ko_labels_independent.csv")
    C.to_csv(OUT / "B_ko_within_cm_independent.csv")
    print(G.round(4).to_string())
    print(C.round(4).to_string())
    print("labels with CI excluding 0:", L.groupby("variant")["excl0"].sum().to_dict())
    # KO's own labels through the lineage map (no classifier) vs WT annotation through the lineage map
    kol = ko.obs["celltype"].astype(str).to_numpy()
    own = pd.Series([GM.get(x, "Unknown") for x in kol]).value_counts(normalize=True)
    wta = pd.Series([GM.get(x, "Unknown") for x in lw]).value_counts(normalize=True)
    print("group fractions by each file's OWN labels via the lineage map:")
    print(pd.DataFrame({"wt_own": wta, "ko_own": own}).fillna(0).round(4).to_string())


def part_c():
    """Group fractions by section-rank quintile, WT (OOF top calls, LR) and KO (top calls, LR, 5 draws)."""
    sec_w, sec_k = STATE["sec_w"], STATE["sec_k"]
    s = STATE["logreg"]
    order, M = s["order"], s["M"]
    out = []
    for emb, P, sec in (("WT", s["P_oof"], sec_w), ("KO", s["Pk"], sec_k)):
        z = np.array([float(x[1:]) for x in sec])
        levels = np.unique(z)
        rank = np.searchsorted(levels, z) / (len(levels) - 1)
        q = np.minimum((rank * 5).astype(int), 4)
        g = np.array(order)[(P @ M).argmax(1)]
        ct = pd.crosstab(q, g, normalize="index").reindex(columns=order, fill_value=0.0)
        ct.insert(0, "n_cells", np.bincount(q))
        ct.insert(0, "embryo", emb)
        out.append(ct)
        # share of each group's cells in the outer 10% of sections at each end
        lo, hi = rank <= 0.1, rank >= 0.9
        ends = pd.DataFrame({grp: [float((lo & (g == grp)).sum() / max(1, (g == grp).sum())),
                                   float((hi & (g == grp)).sum() / max(1, (g == grp).sum()))] for grp in order},
                            index=[f"{emb}_low10%", f"{emb}_high10%"])
        out.append(ends.assign(embryo=emb))
        if emb == "KO":
            gaps = np.diff(levels)
            j = int(np.argmax(gaps))
            big = pd.Series(sec).value_counts()
            log(f"C KO: largest z gap {gaps[j]:.0f} between section ranks {j} and {j + 1} of {len(levels)}; "
                f"largest section {big.index[0]} ({big.iloc[0]} cells) at rank "
                f"{int(np.searchsorted(levels, float(big.index[0][1:])))}")
    df = pd.concat(out)
    df.to_csv(OUT / "C_group_by_z_quintile.csv")
    print(df.round(3).to_string())


def part_d():
    sec_w, sec_k = STATE["sec_w"], STATE["sec_k"]
    s = STATE["logreg"]
    order, M, cls = s["order"], s["M"], s["cls"]
    gw, gk = (s["P_oof"] @ M).argmax(1), (s["Pk"] @ M).argmax(1)
    Gw, Gk = onehot(gw, len(order)), onehot(gk, len(order))
    Lw, Lk = onehot(s["P_oof"].argmax(1), len(cls)), onehot(s["Pk"].argmax(1), len(cls))
    rows = []
    # 1. null: odd vs even sections inside one embryo (same anatomy, interleaved)
    for emb, Gx, Lx, sec in (("WT", Gw, Lw, sec_w), ("KO", Gk, Lk, sec_k)):
        z = np.array([float(x[1:]) for x in sec])
        rank = np.searchsorted(np.unique(z), z)
        odd = rank % 2 == 1
        for lvl, W, names in (("group", Gx, order), ("label", Lx, list(cls))):
            t = diff_table(W[~odd], sec[~odd], W[odd], sec[odd], names, n_boot=1000, seed=3)
            rows.append(dict(test="odd vs even sections", embryo=emb, level=lvl, n_cat=len(names),
                             n_excl0=int(t["excl0"].sum()), max_abs_diff=float(t["diff"].abs().max())))
            # also blocks of 4 consecutive sections alternating (closer to a different-coverage null)
            blk = (rank // 4) % 2 == 1
            t2 = diff_table(W[~blk], sec[~blk], W[blk], sec[blk], names, n_boot=1000, seed=4)
            rows.append(dict(test="alternating 4-section blocks", embryo=emb, level=lvl, n_cat=len(names),
                             n_excl0=int(t2["excl0"].sum()), max_abs_diff=float(t2["diff"].abs().max())))
    null = pd.DataFrame(rows)
    null.to_csv(OUT / "D_bootstrap_null.csv", index=False)
    print(null.to_string())
    # 2. trimming sensitivity of KO - WT group differences
    def trimmed(sec, lo_k, hi_k):
        z = np.array([float(x[1:]) for x in sec])
        lv = np.unique(z)
        r = np.searchsorted(lv, z)
        return (r >= lo_k) & (r <= len(lv) - 1 - hi_k)

    trows = []
    for tag, (wl, wh, kl, kh) in {"none": (0, 0, 0, 0), "both ends -2": (2, 2, 2, 2), "both ends -4": (4, 4, 4, 4),
                                  "both ends -6": (6, 6, 6, 6), "WT low -6": (6, 0, 0, 0), "WT high -6": (0, 6, 0, 0),
                                  "KO low -6": (0, 0, 6, 0), "KO high -6": (0, 0, 0, 6)}.items():
        mw, mk = trimmed(sec_w, wl, wh), trimmed(sec_k, kl, kh)
        t = diff_table(Gw[mw], sec_w[mw], Gk[mk], sec_k[mk], order, n_boot=1000, seed=5)
        trows.append(t["diff"].rename(tag))
        trows.append(t["excl0"].rename(tag + " excl0"))
    T = pd.concat(trows, axis=1)
    T.to_csv(OUT / "D_trim_sensitivity.csv")
    print(T.round(4).to_string())
    # 3. single-section jackknife: largest change of each group's KO-WT difference when one section is dropped
    base = Gk.mean(0) - Gw.mean(0)
    worst = np.zeros(len(order))
    for emb, Gx, sec in (("WT", Gw, sec_w), ("KO", Gk, sec_k)):
        for s_id in np.unique(sec):
            m = sec != s_id
            d = (Gk[m].mean(0) - Gw.mean(0)) if emb == "KO" else (Gk.mean(0) - Gw[m].mean(0))
            worst = np.maximum(worst, np.abs(d - base))
    J = pd.DataFrame({"diff": base, "max_jackknife_shift": worst}, index=order)
    J.to_csv(OUT / "D_jackknife.csv")
    print(J.round(4).to_string())


def _est(P, priors):
    return {"top call": np.bincount(P.argmax(1), minlength=P.shape[1]) / len(P), "mean probability": P.mean(0),
            "EM": H.em_priors(P, priors, tol=1e-8)[0]}


def part_e():
    rows = []
    for name in ("E9.5", "E8.75"):
        a = load(name)
        X = a.X.toarray()
        lab = a.obs["celltype"].astype(str).to_numpy()
        ok = lab != "Unknown"
        z = np.round(np.asarray(a.obsm["spatial_3D"])[:, 2].astype(float))
        upper = z > np.median(np.unique(z))
        U = np.unique(lab[ok])
        Ug, Mu = group_matrix(U)
        for tr_h, q_h in (("lower", "upper"), ("upper", "lower")):
            trm = (~upper if tr_h == "lower" else upper) & ok
            qm = (upper if tr_h == "lower" else ~upper) & ok
            vc = pd.Series(lab[trm]).value_counts()
            kept = vc.index[vc >= 10]
            trm &= np.isin(lab, kept)
            Xtr, ytr = X[trm], lab[trm]
            P_oof, cls = oof_proba(Xtr, ytr, "logreg")
            m = pipe().fit(Xtr, ytr)
            Pq = m.predict_proba(X[qm])
            priors = np.array([(ytr == c).mean() for c in cls])
            pos = np.searchsorted(U, cls)
            Tq = np.array([(lab[qm] == u).mean() for u in U])
            Ttr = np.array([(ytr == u).mean() for u in U])
            Etr, Eq = _est(P_oof, priors), _est(Pq, priors)
            absent = 1 - np.isin(lab[qm], cls).mean()
            for est in ("top call", "mean probability", "EM"):
                e_tr, e_q = np.zeros(len(U)), np.zeros(len(U))
                e_tr[pos], e_q[pos] = Etr[est], Eq[est]
                for lvl, f in (("label", lambda v: v), ("group", lambda v: v @ Mu)):
                    dT, dE = f(Tq) - f(Ttr), f(e_q) - f(e_tr)
                    rows.append(dict(stage=name, train_half=tr_h, level=lvl, estimator=est,
                                     query_cells_of_untrained_labels=float(absent),
                                     true_tv=float(np.abs(dT).sum() / 2), delta_slope=float(dE @ dT / (dT @ dT)),
                                     delta_mae=float(np.abs(dE - dT).mean()), delta_max=float(np.abs(dE - dT).max()),
                                     query_comp_mae=float(np.abs(f(e_q) - f(Tq)).mean()),
                                     train_comp_mae=float(np.abs(f(e_tr) - f(Ttr)).mean())))
            acc_q = float((cls[Pq.argmax(1)] == lab[qm]).mean())
            log(f"E {name} train {tr_h}: query accuracy {acc_q:.3f}, untrained-label share {absent:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "E_estimators_shifted_query.csv", index=False)
    print(df.round(4).to_string())
    print(df.groupby(["level", "estimator"])[["delta_slope", "delta_mae"]].mean().round(4).to_string())


def part_f():
    e875, e95 = load("E8.75"), load("E9.5")
    d95 = H.median_depth(e95)
    rows = []

    def acc(m, q, p, tag, s):
        true = q.obs["group"].astype(str).to_numpy()
        tf = (true != "Unknown") & np.isin(true, m.groups)
        g = p["group_raw"].to_numpy()
        rows.append(dict(direction=tag, seed=s, group_accuracy=float((g[tf] == true[tf]).mean()),
                         group_macro_f1=float(f1_score(true[tf], g[tf], labels=np.unique(true[tf]), average="macro",
                                                       zero_division=0)),
                         novelty_gt1=float((p["novelty"] > 1).mean())))
        log(f"F {tag} seed {s}: {rows[-1]}")

    for s in (0, 1, 2):
        mt = H.fit_transfer(e875, "celltype", classifier="logreg", depth=d95, seed=s, n_jobs=16)
        acc(mt, e95, H.predict(mt, e95, seed=s), "E8.75(thinned)->E9.5", s)
    for s in (0, 1):
        m = H.fit_transfer(e875, "celltype", classifier="logreg", seed=s, n_jobs=16)
        acc(m, e95, H.predict(m, e95, seed=s), "E8.75->E9.5", s)
    mb = H.fit_transfer(e95, "celltype", classifier="logreg", seed=0, n_jobs=16)
    for s in (0, 1, 2):
        acc(mb, e875, H.predict(mb, e875, depth_draws=1, seed=100 * s), "E9.5->E8.75 1 draw", s)
        acc(mb, e875, H.predict(mb, e875, depth_draws=5, seed=100 * s), "E9.5->E8.75 5 draws", s)
    acc(mb, e875, H.predict(mb, e875, depth_match=False), "E9.5->E8.75 no depth match", 0)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "F_cross_seeds.csv", index=False)
    print(df.round(4).to_string())


def part_g():
    Xw, lw, sec_w, sec_k = STATE["Xw"], STATE["lw"], STATE["sec_w"], STATE["sec_k"]
    s = STATE["logreg"]
    ko = load("KO")
    kol = ko.obs["celltype"].astype(str).to_numpy()
    cmw = np.isin(lw, CM_LABELS)
    P_oof, cls = oof_proba(Xw[cmw], lw[cmw], "logreg")
    m = pipe().fit(Xw[cmw], lw[cmw])
    kcm = (s["Pk"] @ s["M"]).argmax(1) == s["order"].index("CM")
    Pk = np.mean([m.predict_proba(X[kcm]) for X in STATE["Xk_draws"]], axis=0)
    t = diff_table(onehot(P_oof.argmax(1), 4), sec_w[cmw], onehot(Pk.argmax(1), 4), sec_k[kcm], list(cls))
    t.to_csv(OUT / "G_within_cm_cm_only_classifier.csv")
    print("CM-only classifier: WT annotated-CM cells (OOF) vs KO cells called CM by the full LR")
    print(t.round(4).to_string())
    print("WT annotation within CM:", pd.Series(lw[cmw]).value_counts(normalize=True).round(4).to_dict())
    own = pd.crosstab(kol[kcm], np.asarray(cls)[Pk.argmax(1)], normalize="index")
    print("KO own label (rows) -> CM-only call (cols), KO cells called CM:")
    print(own.round(3).to_string())
    ownfrac = pd.Series(kol[np.isin(kol, ["vCM", "SV-CM", "OFT/RV-CM"])]).value_counts(normalize=True)
    print("KO own CM labels:", ownfrac.round(4).to_dict())


def part_h():
    wt, ko = load("E9.5"), load("KO")
    lw = wt.obs["celltype"].astype(str).to_numpy()
    kol = ko.obs["celltype"].astype(str).to_numpy()
    genes = np.asarray(wt.var_names)
    Cw = sparse.csr_matrix(wt.layers["counts"]).astype(np.float64)
    Ck = thin(ko.layers["counts"], STATE["dw"] / STATE["dk"], np.random.default_rng(99)).astype(np.float64)

    def cp10k_mean(C, mask):
        Cm = C[mask]
        lib = np.asarray(Cm.sum(1)).ravel()
        return np.asarray((sparse.diags(1e4 / lib) @ Cm).mean(0)).ravel()

    pairs = [("Hepatocyte", "Hepatocytes"), ("Surface Ectoderm", "SE"), ("Foregut", "Foregut"),
             ("Branchial Arch", "NCC"), ("vCM", "V-CM"), ("Pericardium", "Peri")]
    R = {}
    for kl, wl in pairs:
        mk, mw = kol == kl, lw == wl
        a, b = cp10k_mean(Ck, mk), cp10k_mean(Cw, mw)
        ok = (a >= 5) & (b >= 5)
        r = np.full(len(genes), np.nan)
        r[ok] = np.log2(a[ok] / b[ok])
        R[f"{kl}|{wl}"] = r
        log(f"H {kl} ({mk.sum()}) vs {wl} ({mw.sum()}): {ok.sum()} genes >=5 CP10k in both; median |log2 KO/WT| "
            f"{np.nanmedian(np.abs(r)):.3f}; share |log2|>1 {np.nanmean(np.abs(r) > 1):.3f}")
    Rdf = pd.DataFrame(R, index=genes)
    corr = Rdf.corr(method="spearman", min_periods=20)
    n = Rdf.notna().astype(int).T @ Rdf.notna().astype(int)
    print("Spearman correlation of per-gene log2(KO/WT) between matched populations (pairwise complete genes):")
    print(corr.round(3).to_string())
    print("genes in common per pair:")
    print(n.to_string())
    Rdf.to_csv(OUT / "H_gene_log2_ratios.csv")
    # corrected re-call: divide KO CP10k by the per-gene ratio estimated from the non-Mesp1 pairs (median over pairs)
    nonm = Rdf[["Hepatocyte|Hepatocytes", "Surface Ectoderm|SE", "Foregut|Foregut", "Branchial Arch|NCC"]]
    est = nonm.median(1, skipna=True)
    cnt = nonm.notna().sum(1)
    est[(cnt < 2)] = 0.0
    est = est.fillna(0.0).clip(-2, 2)
    log(f"H correction: {int((est != 0).sum())} genes with an estimate; median |log2| {np.median(np.abs(est[est != 0])):.3f}")
    s = STATE["logreg"]
    fac = 2.0 ** est.to_numpy()
    Xc = [np.log1p(np.expm1(X) / fac[None, :]).astype(np.float32) for X in STATE["Xk_draws"]]
    Pc = np.mean([s["model"].predict_proba(X) for X in Xc], axis=0)
    order, M, cls = s["order"], s["M"], s["cls"]
    gk0, gkc = (s["Pk"] @ M).argmax(1), (Pc @ M).argmax(1)
    gw = (s["P_oof"] @ M).argmax(1)
    comp = pd.DataFrame({"wt_oof": np.bincount(gw, minlength=len(order)) / len(gw),
                         "ko": np.bincount(gk0, minlength=len(order)) / len(gk0),
                         "ko_gene_corrected": np.bincount(gkc, minlength=len(order)) / len(gkc)}, index=order)
    cmi = [list(cls).index(c) for c in CM_LABELS]
    li0, lic, liw = s["Pk"].argmax(1), Pc.argmax(1), s["P_oof"].argmax(1)
    within = pd.DataFrame({"wt_oof": [np.mean(liw[np.isin(liw, cmi)] == i) for i in cmi],
                           "ko": [np.mean(li0[np.isin(li0, cmi)] == i) for i in cmi],
                           "ko_gene_corrected": [np.mean(lic[np.isin(lic, cmi)] == i) for i in cmi]}, index=CM_LABELS)
    print("group composition, KO recalled after per-gene correction:")
    print(comp.round(4).to_string())
    print("within-CM composition after correction:")
    print(within.round(4).to_string())
    print(f"label changes after correction: {np.mean(li0 != lic):.3f} of KO cells")
    comp.to_csv(OUT / "H_corrected_groups.csv")
    within.to_csv(OUT / "H_corrected_within_cm.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="A,B,C,D,E,F,G,H")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for p in args.parts.split(","):
        log(f"===== part {p}")
        try:
            {"A": part_a, "B": part_b, "C": part_c, "D": part_d, "E": part_e, "F": part_f, "G": part_g,
             "H": part_h}[p]()
        except Exception:
            traceback.print_exc()
    log("done")


if __name__ == "__main__":
    main()
