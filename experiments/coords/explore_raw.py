"""Quick structural look at the coordinates of every released MERFISH file (run on the DGX).

Prints obsm keys, raw obs columns and the first facts about obsm["spatial_3D"]; no expression values.
Usage: uv run python experiments/coords/explore_raw.py
"""

import anndata as ad
import numpy as np

FILES = [
    "data/raw/T2/E6.75.h5ad",
    "data/raw/T2/E7.25.h5ad",
    "data/raw/T2/E8.0.h5ad",
    "data/raw/T2/E8.25_late.h5ad",
    "data/raw/T2/E8.75.h5ad",
    "data/raw/T2/E9.5.h5ad",
    "data/raw/T3/E9.5_mab21l2_ko.h5ad",
]

for f in FILES:
    a = ad.read_h5ad(f, backed="r")
    print(f"== {f}  n={a.n_obs}")
    print("   obs columns:", list(a.obs.columns))
    print("   obsm keys:", {k: (np.asarray(v).shape, np.asarray(v).dtype) for k, v in a.obsm.items()})
    C = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)
    for j, name in enumerate("xyz"):
        v = C[:, j]
        u = np.unique(v)
        frac_int = float(np.mean(np.isclose(v, np.round(v))))
        print(f"   {name}: min {v.min():.3f} max {v.max():.3f} range {np.ptp(v):.3f} unique {u.size} "
              f"frac_integer {frac_int:.3f}")
    z = np.unique(C[:, 2])
    print("   first 12 z levels:", np.round(z[:12], 4).tolist())
    d = np.diff(z)
    vals, cnt = np.unique(np.round(d, 3), return_counts=True)
    order = np.argsort(-cnt)[:10]
    print("   z-step counts (top 10):", [(float(vals[i]), int(cnt[i])) for i in order])
    for col in a.obs.columns:
        s = a.obs[col]
        if s.dtype.kind in "iuf":
            print(f"   obs[{col}] numeric: min {s.min()} max {s.max()} unique {s.nunique()}")
        else:
            print(f"   obs[{col}] categorical/str: unique {s.nunique()} e.g. {list(map(str, s.unique()[:6]))}")
    other = [k for k in a.obsm.keys() if k != "spatial_3D"]
    for k in other:
        M = np.asarray(a.obsm[k], dtype=np.float64)
        if M.ndim == 2 and M.shape[1] <= 3:
            print(f"   obsm[{k}] ranges:", [(round(float(M[:, j].min()), 2), round(float(M[:, j].max()), 2))
                                          for j in range(M.shape[1])])
    a.file.close()
