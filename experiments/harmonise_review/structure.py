"""Structure checks for the harmonise review: z levels -> sections, sample ids, coverage. Structure only."""
import numpy as np
import pandas as pd

from vec.harmonise import section_ids
from vec.prep import load_stage

pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 200)
FILES = {"E9.5": "data/raw/T2/E9.5.h5ad", "E8.75": "data/raw/T2/E8.75.h5ad", "KO": "data/raw/T3/E9.5_mab21l2_ko.h5ad"}

for name, f in FILES.items():
    a = load_stage(f)
    xyz = np.asarray(a.obsm["spatial_3D"], dtype=float)
    z = xyz[:, 2]
    uz, cnt = np.unique(z, return_counts=True)
    sec = section_ids(a)
    zr = np.round(z, 0)
    ur = np.unique(zr)
    print(f"===== {name}: {a.n_obs} cells; obs columns {a.obs.columns.tolist()}")
    print(f"raw z levels {len(uz)}; rounded levels {len(ur)}; section_ids {len(np.unique(sec))}")
    d = np.diff(ur)
    print("rounded-level steps:", dict(zip(*np.unique(d, return_counts=True), strict=True)))
    frac = np.abs(uz - np.round(uz))
    print("raw levels within 0.1 of an x.5 boundary:", int((np.abs(frac - 0.5) < 0.1).sum()))
    df = pd.DataFrame({"z": z, "zr": zr})
    g = df.groupby("zr")["z"].agg(["min", "max", "count", "nunique"])
    g["span"] = g["max"] - g["min"]
    print("per rounded level: raw-z span max", round(float(g["span"].max()), 4),
          "; levels with >1 raw value:", int((g["nunique"] > 1).sum()))
    print("cells per section: min", int(g["count"].min()), "median", int(g["count"].median()), "max",
          int(g["count"].max()))
    if name == "E9.5":
        print("raw z values (value:count), first 98:")
        print(", ".join(f"{v:.3f}:{c}" for v, c in zip(uz, cnt, strict=True)))
    if "sample_id" in a.obs:
        sid = a.obs["sample_id"].astype(str).to_numpy()
        ct = pd.crosstab(sid, sec)
        print(f"sample_id: {len(np.unique(sid))} ids; sections per sample_id: "
              f"{dict(zip(*np.unique((ct > 0).sum(1), return_counts=True), strict=True))}; "
              f"sample_ids per section: {dict(zip(*np.unique((ct > 0).sum(0), return_counts=True), strict=True))}")
        print("example sample_id values:", sorted(np.unique(sid))[:5])
    for col in ("genotype", "condition"):
        if col in a.obs:
            print(col, a.obs[col].astype(str).value_counts().to_dict())
    # xy layout per section: centroid drift and spread (multiple embryos would show as large spread / bimodality)
    xy = pd.DataFrame({"x": xyz[:, 0], "y": xyz[:, 1], "sec": zr})
    s = xy.groupby("sec").agg(xc=("x", "mean"), yc=("y", "mean"), xsd=("x", "std"), ysd=("y", "std"), n=("x", "size"))
    print("xy: overall sd", np.round(xyz[:, :2].std(0), 1), "; per-section sd median", np.round(
        s[["xsd", "ysd"]].median().to_numpy(), 1), "; centroid sd across sections",
        np.round(s[["xc", "yc"]].std().to_numpy(), 1))
