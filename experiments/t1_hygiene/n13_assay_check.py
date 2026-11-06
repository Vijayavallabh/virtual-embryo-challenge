"""N13 (experiments/harness/TRIES.md): does the released T1 E9.5 file look like single nuclei or whole cells?

R2's lines N1 (mitochondrially encoded genes) and N2 (nucleus-retained ncRNAs) rest on T1 being nuclei. The organisers'
Data-page text says so (docs/DATA_NOTES.md lines 152-153: the 10x Chromium Nuclei Isolation Kit, and doublets filtered
with Scrublet and AMULET, which needs ATAC fragments from nuclei). The 2026-10-18 loader note (line 93) guessed whole
cells from library size alone.

N13 measures, per cell of the released E9.5 file, the share of its UMIs in the mt- genes, Malat1, R2's N2 module and the
ribosomal-protein genes. It compares them with the same pooled shares in Feng's E9.5 whole cells (X002), over all cells
and over M1's CM group, and applies the reading its row states. It decides nothing by itself: the reading goes to the
captain with R2's sign-off. It reads no E10.5 file, target, score, board or built file, and it prints shares only,
never a cell count.

Usage (DGX, repo root):
    OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
        uv run --frozen --no-sync python -I -B experiments/t1_hygiene/n13_assay_check.py
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in ("t1_m1", "t1_direction", "external"):
    sys.path.insert(0, str(REPO / "experiments" / _p))
sys.path.insert(0, str(HERE))
import dirlib as dl  # noqa: E402
import feng_dirs  # noqa: E402
import filterio as fio  # noqa: E402
import hygiene  # noqa: E402

from vec import provenance  # noqa: E402
from vec.prep import load_stage  # noqa: E402

STAGE = "E9.5"
T1_FILE = REPO / "data" / "raw" / "T1" / f"{STAGE}_RNA.h5ad"
OUT = HERE / "results" / "n13_assay_check"
MT = re.compile(r"^mt-")
RP = re.compile(r"^Rp[ls]\d")
SETS = ("mt", "Malat1", "N2", "Rp")


def gene_sets(genes, n2) -> dict[str, np.ndarray]:
    """Column indices of each set in a gene list: mt- genes, Malat1, R2's N2 module, ribosomal-protein genes."""
    g = [str(x) for x in genes]
    n2s = set(n2)
    return {"mt": np.flatnonzero([bool(MT.match(x)) for x in g]),
            "Malat1": np.flatnonzero([x == "Malat1" for x in g]),
            "N2": np.flatnonzero([x in n2s for x in g]),
            "Rp": np.flatnonzero([bool(RP.match(x)) for x in g])}


def cell_shares(X, sets: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Per cell, the share of its UMIs in each set; cells without a UMI are left out."""
    X = sparse.csr_matrix(X)
    total = np.asarray(X.sum(axis=1)).ravel()
    keep = total > 0
    out = {}
    for k, idx in sets.items():
        s = np.asarray(X[:, idx].sum(axis=1)).ravel() if len(idx) else np.zeros(X.shape[0])
        out[k] = s[keep] / total[keep]
    return out


def summarise(v: np.ndarray) -> dict:
    if v.size == 0:
        return {"mean": None, "median": None, "p10": None, "p90": None}
    return {"mean": float(v.mean()), "median": float(np.median(v)), "p10": float(np.quantile(v, 0.1)),
            "p90": float(np.quantile(v, 0.9))}


def pooled_shares(cp: np.ndarray, sets: dict[str, np.ndarray]) -> dict[str, float]:
    """Feng: a set's summed mean CP10k over 10,000, which is the mean per-cell share."""
    return {k: (float(cp[idx].sum()) / 1e4 if len(idx) else 0.0) for k, idx in sets.items()}


def reading(t1_mt: float, feng_mt: float, t1_malat: float, feng_malat: float, n_mt_genes: int) -> str:
    """N13's reading as its row states it before the run, on the means over all cells."""
    if n_mt_genes == 0:
        if t1_malat >= 2 * feng_malat:
            return "nuclei"
        if t1_malat <= feng_malat:
            return "whole cells"
        return "unclear"
    if t1_mt < feng_mt / 5 and t1_malat >= 2 * feng_malat:
        return "nuclei"
    if t1_mt >= feng_mt / 2:
        return "whole cells"
    return "unclear"


def _pct(x: float | None) -> str:
    return "n/a" if x is None else f"{100 * x:.3f}%"


def report_text(res: dict) -> str:
    lines = [f"N13: does the released T1 {STAGE} file look like nuclei or whole cells? Reading: {res['reading']}",
             f"T1 genes: {res['t1_genes']['mt']} mt- genes, Malat1 present {res['t1_genes']['Malat1'] > 0}, "
             f"{res['t1_genes']['N2']} N2 genes, {res['t1_genes']['Rp']} ribosomal-protein genes",
             f"Feng genes: {res['feng_genes']['mt']} mt- genes, Malat1 present {res['feng_genes']['Malat1'] > 0}",
             "", "Share of UMIs (T1: per-cell mean [median, p10-p90]; Feng: pooled mean per-cell share)"]
    for scope in ("all", "CM"):
        lines.append(f"  {scope} cells")
        for k in SETS:
            t = res["t1"][scope][k]
            f = res["feng"][scope][k] if res["feng"][scope] else None
            lines.append(f"    {k:7s} T1 {_pct(t['mean'])} [{_pct(t['median'])}, {_pct(t['p10'])}-{_pct(t['p90'])}]"
                         f"   Feng {_pct(f)}")
    s = res["stamp"]
    lines += ["", f"Inputs: {res['inputs']}", f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="N13: nuclei or whole cells in the released T1 E9.5 file")
    ap.add_argument("--t1-file", default=str(T1_FILE))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N13 has already written its outputs; a rerun needs a new TRIES row")
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    cfg = feng_dirs.load_config()
    n2 = hygiene.read_extract("N2")

    ref = load_stage(Path(a.t1_file), lineage=False)
    t1_sets = gene_sets(ref.var_names, n2)
    labels = ref.obs["celltype"].astype(str).to_numpy()
    groups = feng_dirs.label_groups(np.unique(labels), cfg)
    cm = np.array([groups.get(lab) == "CM" for lab in labels])
    t1_all, t1_cm = cell_shares(ref.X, t1_sets), cell_shares(ref.X[cm], t1_sets)

    path = feng_dirs.clean_stats_path()
    entry = fio.verify_manifest_output(path)
    st = dl.load_stats(path)
    f_sets = gene_sets(st["genes"], n2)
    f_all = dl.agg(st, lambda p: p[0] == STAGE)
    f_cm = dl.feng_sel(st, STAGE, cfg["external_labels"]["CM"][0])

    res = {"t1": {"all": {k: summarise(v) for k, v in t1_all.items()},
                  "CM": {k: summarise(v) for k, v in t1_cm.items()}},
           "feng": {"all": pooled_shares(f_all["cp"], f_sets),
                    "CM": pooled_shares(f_cm["cp"], f_sets) if f_cm else None},
           "t1_genes": {k: int(len(v)) for k, v in t1_sets.items()},
           "feng_genes": {k: int(len(v)) for k, v in f_sets.items()}}
    res["reading"] = reading(res["t1"]["all"]["mt"]["mean"] or 0.0, res["feng"]["all"]["mt"],
                             res["t1"]["all"]["Malat1"]["mean"] or 0.0, res["feng"]["all"]["Malat1"],
                             res["t1_genes"]["mt"])
    res["inputs"] = {"t1_file": str(Path(a.t1_file).relative_to(REPO)) if Path(a.t1_file).is_relative_to(REPO)
                     else str(a.t1_file), "t1_file_sha256": hygiene.sha256_file(a.t1_file),
                     "feng_stats": str(path.relative_to(REPO)), "feng_stats_sha256": entry["sha256"],
                     "n2_extract_sha256": hygiene.list_sha256(n2)}
    res["stamp"] = {"commit": commit, "tree_clean": unclean == [],
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n13_assay_check.json").write_text(json.dumps(res, indent=1) + "\n")
    text = report_text(res)
    (out / "n13_assay_check.txt").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
