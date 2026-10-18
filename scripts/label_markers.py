"""Mean expression of lineage marker genes per cell-type label, per file (run on the DGX).

Used to place each stage-specific label into a coarse lineage group (configs/lineage_map.yaml) from data rather
than from label names alone. Output: runs/label_markers.tsv (label x marker means, log1p CP10k) and a short print.
"""

import glob
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

MARKERS = {
    "CM": ["Myl7", "Myl4", "Ttn", "Myh7", "Nppa", "Actc1", "Tnni1", "Myl3"],
    "SHF": ["Isl1", "Tbx1", "Fgf10", "Fgf8", "Hand2"],
    "FHF_LV": ["Hand1", "Tbx5", "Hcn4"],
    "Endothelium": ["Cdh5", "Pecam1", "Emcn", "Kdr", "Sox17", "Npr3", "Plvap"],
    "Epicardium": ["Wt1", "Tbx18", "Upk3b", "Tcf21", "Aldh1a2"],
    "Endoderm": ["Foxa2", "Sox17", "Cdh1", "Epcam", "Nkx2-1"],
    "Liver_VE": ["Afp", "Apoa1", "Apoa2", "Ttr", "Alb"],
    "Neural": ["Sox2", "Sox1", "Pax6", "Otx2", "Sox3"],
    "NCC": ["Sox10", "Foxd3", "Tfap2a", "Tfap2b", "Twist1"],
    "Ectoderm": ["Krt8", "Krt18", "Tfap2a", "Wnt6", "Grhl2"],
    "Paraxial": ["Meox1", "Tcf15", "Tbx6", "Pax3", "Foxc1"],
    "Blood": ["Gata1", "Hba-x", "Hbb-y", "Klf1", "Runx1"],
    "ExEM": ["Hand1", "Pitx1", "Postn", "Ahnak", "Tbx20"],
    "Epiblast_PS": ["Pou5f1", "T", "Mixl1", "Nanog", "Fgf5"],
    "Prolif": ["Top2a", "Cdk1", "Mki67", "Ccna2"],
}


def main(root="data/raw", out="runs/label_markers.tsv"):
    rows = []
    for f in sorted(glob.glob(f"{root}/*/*.h5ad")):
        a = ad.read_h5ad(f)
        genes = list(a.var_names)
        X = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(np.asarray(a.X))
        labels = a.obs["celltype"].astype(str).values
        present = {g for gs in MARKERS.values() for g in gs if g in genes}
        idx = {g: genes.index(g) for g in present}
        for lab in pd.unique(labels):
            m = labels == lab
            sub = X[m]
            row = {"file": f.split("raw/")[1], "label": lab, "n": int(m.sum())}
            for g, j in idx.items():
                row[g] = float(sub[:, j].mean())
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    score_cols = []
    for grp, gs in MARKERS.items():
        cols = [g for g in gs if g in df.columns]
        if cols:
            df[f"S_{grp}"] = df[cols].mean(1)
            score_cols.append(f"S_{grp}")
    best = df[score_cols].idxmax(1).str[2:]
    second = df[score_cols].apply(lambda r: r.nlargest(2).index[1][2:], axis=1)
    for (_, r), b, s in zip(df.iterrows(), best, second, strict=True):
        print(f"{r['file']:28s} {r['label']:24s} n={r['n']:6d}  top={b:12s} 2nd={s:12s} " +
              " ".join(f"{c[2:]}={r[c]:.2f}" for c in score_cols if r[c] > 0.3))
    print(f"wrote {out}")


if __name__ == "__main__":
    main(*sys.argv[1:])
