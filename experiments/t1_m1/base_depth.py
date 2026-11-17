"""Target-free listing (decision log 2026-11-17 (e) and 2026-11-18 (a)): how deep the cells of Task 1's second upload's
base (the mid base) are against the cells of its first upload's base (library 3 of the released E10.5 file), by the
builder's own row selection. Reads the released E10.5 file only; builds nothing, reads no board, no target and no
score, and compares no prediction with anything. Prints ratios only.

Usage, repository root:  uv run --frozen --no-sync python -I -B experiments/t1_m1/base_depth.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("build_test_mod", REPO / "experiments/builds/build_test.py")
bt = importlib.util.module_from_spec(spec)
sys.modules["build_test_mod"] = bt
spec.loader.exec_module(bt)
from vec import harness, t1_scope  # noqa: E402

ref = bt.load_stage(bt.DATA / "T1/E10.5_RNA.h5ad", lineage=False)
names = [str(x) for x in ref.obs_names]
library = np.array([n.rsplit("_", 1)[-1] for n in names])
in3 = library == bt.T1_DISSECTED_LIBRARY
C = ref.layers["counts"] if "counts" in ref.layers else None
assert C is not None, list(ref.layers)
tot = np.asarray(C.sum(1)).ravel() if sparse.issparse(C) else np.asarray(C).sum(1)
seed = harness.BUILD_SEED_OFFSET + 1
rows, record = bt._t1_mid_rows(ref, in3, seed)
b = harness._shallow(ref)
t1_scope.classify(b, source="E10.5_RNA", spec=t1_scope.load_spec())
cat = np.asarray(b.obs[t1_scope.CATEGORY_KEY]).astype(str)
cats = sorted(record)
mid = np.zeros(len(tot), bool)
mid[rows] = True
print("mean counts per cell, mid base over library 3, per chamber category (ratios only):")
r_all = []
for c in cats:
    m3 = in3 & (cat == c)
    mm = mid & (cat == c)
    r = float(tot[mm].mean() / tot[m3].mean())
    r_all.append(r)
    print(f"  {c}: library 3's share of the pool {record[c]['library3_share_in_pool']:.3f}; depth ratio {r:.3f}")
ch = np.isin(cat, cats)
pooled = float(tot[mid & ch].mean() / tot[in3 & ch].mean())
print(f"all chamber categories, cell-weighted over each pool: {pooled:.3f}; "
      f"unweighted mean of the categories' ratios {float(np.mean(r_all)):.3f}")
print(f"the other libraries' chamber cells over library 3's: {float(tot[~in3 & ch].mean() / tot[in3 & ch].mean()):.3f}")
