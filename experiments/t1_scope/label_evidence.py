"""Evidence for the chamber-category placements in configs/t1_chamber_classes.yaml (run on the DGX).

For every label of two T1 stages: mean log1p-CP10k and detection rate of lineage marker genes, label-transfer calls
across the two stages (multinomial logistic regression with vec.harmonise on highly variable genes of the training
stage, both directions), and the Spearman correlation of label pseudobulks across stages. Aggregates only (one row per
label); nothing per cell leaves the DGX. Stage names come from the file names (E8.5_RNA.h5ad -> E8.5).

Usage (DGX, repo root):
  uv run python experiments/t1_scope/label_evidence.py data/raw/T1/E8.5_RNA.h5ad data/raw/T1/E9.5_RNA.h5ad
After the release, for the new labels (results go to results/label_evidence_E9.5_E10.5/):
  uv run python experiments/t1_scope/label_evidence.py data/raw/T1/E9.5_RNA.h5ad data/raw/T1/E10.5_RNA.h5ad \
      --out experiments/t1_scope/results/label_evidence_E9.5_E10.5
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import typer
from rich.console import Console
from scipy import sparse
from scipy.stats import spearmanr

from vec.harmonise import fit_transfer, predict
from vec.prep import load_stage

app = typer.Typer(add_completion=False)
console = Console(emoji=False, highlight=False)
OUT = Path(__file__).resolve().parent / "results" / "label_evidence"

MARKERS = {
    "pan_endothelial": ["Cdh5", "Pecam1", "Kdr", "Emcn", "Tek", "Sox18"],
    "endocardial": ["Npr3", "Nfatc1", "Hapln1", "Tmem100", "Cdh11"],
    "vascular": ["Cldn5", "Gja5", "Dll4", "Sox17", "Gja4", "Efnb2", "Plvap", "Lyve1", "Etv2", "Runx1"],
    "cardiomyocyte": ["Myl7", "Tnnt2", "Actc1", "Nppa", "Myh6", "Myl2", "Myl4"],
    "shf": ["Isl1", "Tbx1", "Fgf10", "Fgf8"],
    "epicardial": ["Wt1", "Tbx18", "Tcf21", "Upk3b", "Aldh1a2"],
    "mesothelial": ["Msln", "Krt19", "Upk1b"],
    "septum_hepatic": ["Hlx", "Lhx2", "Afp", "Alb", "Hnf4a", "Ttr"],
    "foregut": ["Foxa2", "Sox2", "Nkx2-1", "Shh", "Cdh1", "Epcam"],
    "neural_crest": ["Sox10", "Foxd3", "Tfap2a", "Twist1", "Prrx1", "Dlx2", "Dlx5"],
    "mesenchyme": ["Postn", "Sox9", "Col1a1", "Pdgfra"],
    "blood": ["Hba-x", "Hbb-y", "Gata1"],
}


def marker_table(a, stage: str) -> pd.DataFrame:
    genes = [g for gs in MARKERS.values() for g in gs if g in a.var_names]
    idx = a.var_names.get_indexer(genes)
    X = sparse.csr_matrix(a.X)[:, idx]
    labels = a.obs["celltype"].astype(str).to_numpy()
    rows = []
    for lab in sorted(set(labels)):
        m = labels == lab
        sub = X[m]
        mean = np.asarray(sub.mean(0)).ravel()
        det = np.asarray((sub > 0).mean(0)).ravel()
        row = {"stage": stage, "label": lab, "n": int(m.sum())}
        row.update({f"mean:{g}": round(float(v), 3) for g, v in zip(genes, mean, strict=True)})
        row.update({f"det:{g}": round(float(v), 3) for g, v in zip(genes, det, strict=True)})
        for name, gs in MARKERS.items():
            cols = [f"mean:{g}" for g in gs if g in genes]
            row[f"score:{name}"] = round(float(np.mean([row[c] for c in cols])), 3) if cols else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def hvg(a, n_top: int) -> list[str]:
    """Highest dispersion (variance / mean of log1p-CP10k) among genes detected in >= 1% of cells."""
    X = sparse.csr_matrix(a.X)
    mean = np.asarray(X.mean(0)).ravel()
    sq = np.asarray(X.multiply(X).mean(0)).ravel()
    var = sq - mean**2
    det = np.asarray((X > 0).mean(0)).ravel()
    disp = np.where((det >= 0.01) & (mean > 0), var / np.maximum(mean, 1e-9), -np.inf)
    order = np.argsort(-disp)[:n_top]
    return [str(g) for g in a.var_names[np.sort(order)]]


def transfer(train, query, genes: list[str], seed: int) -> tuple[pd.DataFrame, dict]:
    model = fit_transfer(train[:, genes].copy(), "celltype", genes=genes, classifier="logreg", exclude_labels=(),
                         seed=seed, n_jobs=8)
    calls = predict(model, query[:, genes].copy(), proba=False, seed=seed)
    tab = pd.crosstab(query.obs["celltype"].astype(str).to_numpy(), calls["label_raw"].to_numpy(),
                      normalize="index").round(3)
    tab.index.name = "query_label"
    return tab, {"method": model.method, "n_genes": len(genes), "n_train": int(model.params["n_train"]),
                 "depth_frac": float(calls.attrs.get("depth_frac", 1.0))}


def pseudobulk_corr(a, b, genes: list[str]) -> pd.DataFrame:
    def pb(x):
        X = sparse.csr_matrix(x[:, genes].X)
        labels = x.obs["celltype"].astype(str).to_numpy()
        return {lab: np.asarray(X[labels == lab].mean(0)).ravel() for lab in sorted(set(labels))}

    pa, pbb = pb(a), pb(b)
    out = pd.DataFrame(index=list(pa), columns=list(pbb), dtype=float)
    for la, va in pa.items():
        for lb, vb in pbb.items():
            out.loc[la, lb] = round(float(spearmanr(va, vb).correlation), 3)
    return out


def stage_name(path: Path) -> str:
    return path.name.split("_")[0].removesuffix(".h5ad")


@app.command()
def main(first: Path = typer.Argument(..., help="earlier T1 stage h5ad (e.g. E8.5_RNA.h5ad)"),
         second: Path = typer.Argument(..., help="later T1 stage h5ad (e.g. E9.5_RNA.h5ad)"),
         n_hvg: int = typer.Option(2000, help="highly variable genes per training stage"),
         seed: int = typer.Option(0), out: Path = typer.Option(OUT, help="output directory (aggregates only)")):
    """Write marker means, cross-stage label-transfer calls and pseudobulk correlations per T1 label."""
    out.mkdir(parents=True, exist_ok=True)
    sa, sb = stage_name(first), stage_name(second)
    a = load_stage(first, lineage=False)
    b = load_stage(second, lineage=False)
    markers = pd.concat([marker_table(a, sa), marker_table(b, sb)], ignore_index=True)
    markers.to_csv(out / "markers.csv", index=False)
    info = {"stages": [sa, sb]}
    genes_b = hvg(b, n_hvg)
    tab_ab, info[f"{sa}_on_{sb}_model"] = transfer(b, a, genes_b, seed)
    tab_ab.to_csv(out / f"calls_{sa}_cells_on_{sb}_labels.csv")
    genes_a = hvg(a, n_hvg)
    tab_ba, info[f"{sb}_on_{sa}_model"] = transfer(a, b, genes_a, seed)
    tab_ba.to_csv(out / f"calls_{sb}_cells_on_{sa}_labels.csv")
    union = sorted(set(genes_a) | set(genes_b))
    pseudobulk_corr(a, b, union).to_csv(out / f"pseudobulk_spearman_{sa}_vs_{sb}.csv")
    info["hvg_union"] = len(union)
    (out / "info.json").write_text(json.dumps(info, indent=2) + "\n")
    for stage, tab in ((sa, tab_ab), (sb, tab_ba)):
        for lab in tab.index:
            top = tab.loc[lab].sort_values(ascending=False).head(3)
            console.print(f"{stage} {lab:18s} -> " + ", ".join(f"{k} {v:.2f}" for k, v in top.items()))
    console.print(f"wrote {out}")


if __name__ == "__main__":
    app()
