"""Validate vec.prep on every released file (run on the DGX): count recovery, counts-layer consistency, lineage
coverage and coarse composition. Prints one block per file."""

import glob

import anndata as ad
import numpy as np
from scipy import sparse

from vec.prep import load_stage, log1p_cp10k, recover_counts

for f in sorted(glob.glob("data/raw/*/*.h5ad")):
    raw = ad.read_h5ad(f)
    X = sparse.csr_matrix(raw.X, dtype=np.float32)
    msg = []
    if "counts" in raw.layers:
        Xc = log1p_cp10k(raw.layers["counts"])
        diff = abs(Xc - X)
        msg.append(f"counts-layer->X max|diff|={diff.max():.2e}")
    else:
        _, lib, resid = recover_counts(X)
        msg.append(f"recovered counts: max residual {resid:.2e}, median library {np.median(lib):.0f}, "
                   f"p10-p90 {np.percentile(lib, 10):.0f}-{np.percentile(lib, 90):.0f}")
    a = load_stage(f)
    lib = np.asarray(a.layers["counts"].sum(1)).ravel()
    comp = a.obs["group"].value_counts(normalize=True).round(3)
    print(f"== {f.split('raw/')[1]}  cells={a.n_obs} genes={a.n_vars}  counts/cell median={np.median(lib):.0f}")
    print("   " + "; ".join(msg))
    print(f"   unmapped labels: {a.uns['unmapped_labels']}  core fraction={a.obs['core'].mean():.3f}  "
          f"mesp1 yes/partial/no={[round(float((a.obs['mesp1'] == k).mean()), 3) for k in ('yes', 'partial', 'no')]}")
    print("   groups: " + ", ".join(f"{k}={v}" for k, v in comp.items()))
