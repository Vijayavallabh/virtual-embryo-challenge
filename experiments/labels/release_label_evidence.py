"""Classifier evidence for the cell-type labels of the released validation answers (runbook step 3, run on the DGX).

The pre-registered rule in `configs/lineage_map.yaml` (2026-10-19) maps every label first seen in a new file before
anything else runs on that file; an ambiguous name is checked against its marker genes (`scripts/label_markers.py`)
and its classifier calls with `vec.harmonise` against the nearest released stage, and defaults to shell. This script
gives the classifier half: for each new file it fits a transfer model on the nearest released stage (its own labels,
duplicates dropped as for prediction inputs) and reports, for every label of the new file, where its cells land in the
reference vocabulary (fine labels and lineage groups). It changes no prediction and writes only label-level shares.

    uv run python experiments/labels/release_label_evidence.py [--out experiments/labels/results]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from vec import harmonise as hz
from vec.prep import lineage_map, load_stage

RAW = Path("data/raw")
# new file -> nearest released stage of the same assay and setting
PAIRS = {
    "T1/E10.5_RNA.h5ad": "T1/E9.5_RNA.h5ad",
    "T2/E10.5.h5ad": "T2/E9.5.h5ad",
    "T3/E8.75_gata4_ko_1.h5ad": "T2/E8.75.h5ad",
    "T3/E8.75_gata4_ko_2.h5ad": "T2/E8.75.h5ad",
}
T1_GENES = 3000  # T1: the reference's most variable genes (log1p CP10k) among the shared ones


def gene_set(ref, query) -> list[str]:
    genes = hz.shared_genes(ref, query)
    if len(genes) <= 1000:  # a MERFISH panel: all of it
        return genes
    X = sparse.csr_matrix(ref[:, genes].X)  # log1p CP10k as released
    mean = np.asarray(X.mean(0)).ravel()
    var = np.asarray(X.multiply(X).mean(0)).ravel() - mean ** 2
    return [genes[i] for i in np.argsort(-var)[:T1_GENES]]


def label_table(query_labels: np.ndarray, pred: pd.DataFrame, mapped: set[str]) -> pd.DataFrame:
    rows = []
    for lab in pd.unique(query_labels):
        m = query_labels == lab
        fine = pd.Series(pred["label_raw"].to_numpy()[m]).value_counts(normalize=True)
        grp = pd.Series(pred["group_raw"].to_numpy()[m]).value_counts(normalize=True)
        rows.append({
            "label": lab, "n": int(m.sum()), "in_map_before": lab in mapped,
            **{f"label_{i + 1}": f"{fine.index[i]} {fine.iloc[i]:.2f}" if len(fine) > i else "" for i in range(3)},
            **{f"group_{i + 1}": f"{grp.index[i]} {grp.iloc[i]:.2f}" if len(grp) > i else "" for i in range(2)},
            "novelty_gt1": round(float((pred["novelty"].to_numpy()[m] > 1).mean()), 3),
            "median_conf": round(float(np.median(pred["label_conf"].to_numpy()[m])), 3),
        })
    return pd.DataFrame(rows).sort_values("n", ascending=False)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", type=Path, default=Path("experiments/labels/results"))
    p.add_argument("--n-jobs", type=int, default=8)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    mapped = set(lineage_map()["labels"])
    summary = {}
    for qf, rf in PAIRS.items():
        ref, query = load_stage(RAW / rf, drop_duplicates=True), load_stage(RAW / qf, drop_duplicates=True)
        genes = gene_set(ref, query)
        qd, rd = hz.median_depth(query, genes), hz.median_depth(ref, genes)
        model = hz.fit_transfer(ref, "celltype", genes=genes, depth=qd if rd and qd and rd > qd else None,
                                n_jobs=args.n_jobs)
        pred = hz.predict(model, query)
        t = label_table(query.obs["celltype"].astype(str).to_numpy(), pred, mapped)
        name = f"{Path(qf).stem}__vs__{Path(rf).stem}"
        t.to_csv(args.out / f"{name}.csv", index=False)
        summary[name] = {"query": qf, "reference": rf, "genes": len(genes), "method": model.method,
                         "query_depth": qd, "reference_depth": rd, "query_depth_frac": pred.attrs.get("depth_frac"),
                         "new_labels": int((~t["in_map_before"]).sum())}
        print(f"== {name}: {summary[name]}")
        print(t.to_string(index=False))
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
