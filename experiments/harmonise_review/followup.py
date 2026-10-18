"""Follow-up checks for the harmonise review (run from the repo root on the DGX).

1. within-CM KO-WT differences under z trimming, and within-CM composition per z quintile (harmonise's own per-cell
   outputs: runs/harmonise/cells/oof_E9.5.pkl and ko_pred.pkl)
2. direct group classifier vs summed fine-label probabilities over more CV seeds
"""
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

from vec import harmonise as H
from vec.prep import lineage_map, load_stage

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 60)
CM = ["A-CM", "OFT-CM", "SV-CM", "V-CM"]
OUT = "runs/harmonise_review"

wt = load_stage("data/raw/T2/E9.5.h5ad")
ko = load_stage("data/raw/T3/E9.5_mab21l2_ko.h5ad")
oof = pd.read_pickle("runs/harmonise/cells/oof_E9.5.pkl")
kp = pd.read_pickle("runs/harmonise/cells/ko_pred.pkl")
assert list(kp.index) == list(ko.obs_names) and list(oof.index) == list(wt.obs_names)
sec_w, sec_k = oof["section"].to_numpy(), H.section_ids(ko)


def rank_of(sec):
    z = np.array([float(s[1:]) for s in sec])
    lv = np.unique(z)
    return np.searchsorted(lv, z), len(lv)


rw, nw = rank_of(sec_w)
rk, nk = rank_of(sec_k)
lw, lk = oof["label_raw"].to_numpy(), kp["label_raw"].to_numpy()
cw, ck = np.isin(lw, CM), np.isin(lk, CM)

rows = []
for tag, (wl, wh, kl, kh) in {"none": (0, 0, 0, 0), "both ends -4": (4, 4, 4, 4), "both ends -6": (6, 6, 6, 6),
                              "WT low -6": (6, 0, 0, 0), "WT high -6": (0, 6, 0, 0), "KO low -6": (0, 0, 6, 0),
                              "KO high -6": (0, 0, 0, 6), "KO low -8": (0, 0, 8, 0)}.items():
    mw = cw & (rw >= wl) & (rw <= nw - 1 - wh)
    mk = ck & (rk >= kl) & (rk <= nk - 1 - kh)
    d = H.bootstrap_difference(H.one_hot(lw[mw], CM), sec_w[mw], H.one_hot(lk[mk], CM), sec_k[mk], n_boot=1000)
    for c in CM:
        rows.append(dict(trim=tag, label=c, wt=d.loc[c, "a"], ko=d.loc[c, "b"], diff=d.loc[c, "diff"],
                         lo=d.loc[c, "diff_lo"], hi=d.loc[c, "diff_hi"], excl0=bool(d.loc[c, "ci_excludes_0"])))
T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/F1_within_cm_trim.csv", index=False)
print("within-CM KO-WT differences under z trimming (top call, harmonise per-cell outputs):")
print(T.pivot(index="trim", columns="label", values="diff").round(3).to_string())
print(T.pivot(index="trim", columns="label", values="excl0").to_string())

q_rows = []
for emb, lab, r, n, m in (("WT", lw, rw, nw, cw), ("KO", lk, rk, nk, ck)):
    q = np.minimum((r / (n - 1) * 5).astype(int), 4)
    ct = pd.crosstab(q[m], lab[m], normalize="index").reindex(columns=CM, fill_value=0.0)
    ct.insert(0, "cm_cells", np.bincount(q[m], minlength=5))
    ct.insert(0, "embryo", emb)
    q_rows.append(ct)
Q = pd.concat(q_rows)
Q.to_csv(f"{OUT}/F1_within_cm_by_quintile.csv")
print("within-CM composition by z quintile (0 = lowest z, the liver/septum-transversum end in both stacks):")
print(Q.round(3).to_string())

gmap = {g: g for g in lineage_map()["groups"]}
rows = []
for name, f in (("E9.5", "data/raw/T2/E9.5.h5ad"), ("E8.75", "data/raw/T2/E8.75.h5ad")):
    a = wt if name == "E9.5" else load_stage(f)
    for s in (0, 1, 2):
        o = H.cross_validate(a, method="logreg", n_folds=5, seed=s, n_jobs=16, novelty=False)
        m = H.transfer_metrics(o)
        dg = H.cross_validate(a, label_key="group", method="logreg", n_folds=5, seed=s, n_jobs=16, novelty=False,
                              group_map=gmap)
        rows.append(dict(stage=name, seed=s, summed_group_acc=m["group_accuracy"], summed_group_f1=m["group_macro_f1"],
                         direct_group_acc=float((dg["label_raw"] == dg["true"]).mean()),
                         direct_group_f1=float(f1_score(dg["true"], dg["label_raw"], labels=np.unique(dg["true"]),
                                                        average="macro", zero_division=0)),
                         macro_f1=m["macro_f1"], accuracy=m["accuracy"]))
        print(rows[-1], flush=True)
D = pd.DataFrame(rows)
D["f1_gain_summed"] = D["summed_group_f1"] - D["direct_group_f1"]
D.to_csv(f"{OUT}/F2_direct_vs_summed_seeds.csv", index=False)
print(D.round(4).to_string())
