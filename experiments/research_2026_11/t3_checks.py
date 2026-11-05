"""T3 checks D1, D2, D5 and D6 behind the T3 second-upload ruling (experiments/research_2026_11/REPORT.md, 3 and 7).

Diagnostics, not tries (experiments/harness/TRIES.md rows N6-N9, registered with their decision rules before any run).
None computes a score. None reads a test file or a test-board limit.

  d1    Run fingerprint: per released T2/T3 file, its stream and metadata. Covers the matrix type, layers and obs
        columns; the slide, section, embryo, sample and batch identifiers; z levels and step; synthetic records; and
        per-object counts and genes. Readout: whether the Gata4 and Mab21l2 knockout files share a slide or batch.
  d2    Areal capture: per file, within each z level holding at least 200 objects, the median 2D nearest-neighbour
        distance between objects (an area-free proxy for objects per area) and the median counts per object. Readout:
        the object-density ratio R_obj = (NN_WT / NN_KO)^2, the counts-per-object ratio R_cpo, and counts per area
        R_cpa = R_obj x R_cpo. Pairs: each Gata4 replicate against WT E8.75; Mab21l2 against WT E9.5.
  d5d6  D5 asks whether the component shared by s1's factor (gata4_pooled, Gata4 against WT E8.75) and the barred
        branch's factor (transfer_mab21l2_all, Mab21l2 against WT E8.75) lies along a WT developmental axis built
        without WT E8.75, or across it.
        - The axis d is E9.5 against E8.25_late (primary) or E9.5 against E8.5 (secondary), from the T3 script's own
          estimator (gene_factor: non-Mesp1 labels holding 100 or more cells in both, label weights V.p) in the WT
          E8.75 vocabulary.
        - On centred log-factors: f = a d_hat + r, and <f_G, f_M> = a_G a_M + <r_G, r_M>.
        D6 places WT E8.75 on the ladder.
        - It uses the centred log of the label-weighted, pooled non-Mesp1 mean proportions of each WT file, which is
          depth-free and composition-free.
        - Positions are projections on the E8.25_late -> E9.5 axis.
        - Known-lag check: E8.5 lands strictly between the ends (0 < s85 < 1), and E8.25_late lands before E8.5 on the
          E8.5 -> E9.5 axis.
        - The nominal position of E8.75, calibrated on E8.5, is s85 + 0.25 (1 - s85).

Run on the DGX from the repo root at a clean, committed tree:
    uv run --frozen --no-sync python -I -B experiments/research_2026_11/t3_checks.py d1 --out DIR
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments" / "t3_decomp"))
sys.path.insert(0, str(REPO / "experiments" / "research_2026_11"))

from diagnostics import provenance  # noqa: E402

from vec import harness  # noqa: E402
from vec.prep import load_stage  # noqa: E402

DATA = harness.DATA_DIR
WT = {"E8.25_late": "T2/E8.25_late.h5ad", "E8.5": "T2/E8.5.h5ad", "E8.75": "T2/E8.75.h5ad", "E9.5": "T2/E9.5.h5ad",
      "E10.5": "T2/E10.5.h5ad"}
KO = {"Mab21l2 KO (E9.5)": "T3/E9.5_mab21l2_ko.h5ad", "Gata4 KO rep1 (E8.75)": "T3/E8.75_gata4_ko_1.h5ad",
      "Gata4 KO rep2 (E8.75)": "T3/E8.75_gata4_ko_2.h5ad"}
ID_COLS = ("slide", "section", "embryo", "sample_id", "perturb_batch", "batch", "condition", "genotype")
COORD = "spatial_3D"
RUN_DIR = REPO / "runs/t3_decomp/gata4_20261028T122425Z"
DECISION = REPO / "experiments/t3_decomp/results/gata4/gata4_decisions.json"


def counts_of(a) -> np.ndarray:
    from scipy import sparse
    c = a.layers["counts"] if "counts" in a.layers else a.X
    c = c.tocsr() if sparse.issparse(c) else np.asarray(c)
    tot = np.asarray(c.sum(axis=1)).ravel()
    ng = np.asarray((c > 0).sum(axis=1)).ravel()
    return tot, ng


def levels(coords: np.ndarray) -> np.ndarray:
    return np.round(coords[:, 2].astype(np.float64), 4)


# ---------------------------------------------------------------------------------------------------------------------
def d1() -> tuple[dict, list[Path]]:
    import anndata as ad
    from scipy import sparse
    out, inputs = {}, []
    for name, rel in {**{f"WT {k}": v for k, v in WT.items()}, **KO}.items():
        p = DATA / rel
        inputs.append(p)
        a = ad.read_h5ad(p)
        ded = load_stage(p, drop_duplicates=True)
        tot, ng = counts_of(ded)
        rec = {"file": rel, "matrix": "sparse" if sparse.issparse(a.X) else "dense", "dtype": str(a.X.dtype),
               "layers": sorted(map(str, a.layers.keys())), "obs_columns": sorted(map(str, a.obs.columns)),
               "obsm": sorted(map(str, a.obsm.keys())), "objects_released": int(a.n_obs),
               "synthetic_records_dropped": int(a.n_obs - ded.n_obs)}
        for col in ID_COLS:
            if col in a.obs.columns:
                vc = a.obs[col].astype(str).value_counts()
                rec[f"id:{col}"] = {"n_unique": int(vc.size), "top": {str(k): int(v) for k, v in vc.head(12).items()}}
        if COORD in a.obsm:
            z = np.unique(levels(np.asarray(a.obsm[COORD])))
            rec["z_levels"] = int(z.size)
            rec["z_step_median"] = float(np.median(np.diff(z))) if z.size > 1 else None
        rec["counts_per_object"] = {"min": float(tot.min()), "median": float(np.median(tot))}
        rec["genes_per_object"] = {"min": int(ng.min()), "median": float(np.median(ng))}
        out[name] = rec

    def ids(name, col):
        return set(out[name].get(f"id:{col}", {}).get("top", {}))
    shared = {}
    for col in ID_COLS:
        g = ids("Gata4 KO rep1 (E8.75)", col) | ids("Gata4 KO rep2 (E8.75)", col)
        m = ids("Mab21l2 KO (E9.5)", col)
        if g and m:
            shared[col] = sorted(g & m)
    return {"files": out, "gata4_mab21l2_shared_ids": shared}, inputs


# ---------------------------------------------------------------------------------------------------------------------
def d2() -> tuple[dict, list[Path]]:
    from sklearn.neighbors import NearestNeighbors
    files = {"WT E8.25_late": WT["E8.25_late"], "WT E8.5": WT["E8.5"], "WT E8.75": WT["E8.75"], "WT E9.5": WT["E9.5"],
             **KO}
    out, inputs = {}, []
    for name, rel in files.items():
        p = DATA / rel
        inputs.append(p)
        a = load_stage(p, drop_duplicates=True)
        xyz = np.asarray(a.obsm[COORD], dtype=np.float64)[:, :3]
        tot, _ = counts_of(a)
        lv = levels(xyz)
        nn_med, per_level = [], []
        for z in np.unique(lv):
            i = np.flatnonzero(lv == z)
            if i.size < 200:
                continue
            d, _ = NearestNeighbors(n_neighbors=2).fit(xyz[i, :2]).kneighbors(xyz[i, :2])
            nn_med.append(float(np.median(d[:, 1])))
            per_level.append(int(i.size))
        out[name] = {"file": rel, "levels_used": len(nn_med), "nn2d_median": float(np.median(nn_med)),
                     "objects_per_level_median": float(np.median(per_level)),
                     "counts_per_object_median": float(np.median(tot))}

    def ratio(ko, wt):
        r_obj = (out[wt]["nn2d_median"] / out[ko]["nn2d_median"]) ** 2
        r_cpo = out[ko]["counts_per_object_median"] / out[wt]["counts_per_object_median"]
        r_cpa = r_obj * r_cpo
        if r_obj >= 1.25 and 0.75 <= r_cpa <= 1.25:
            cls = "granularity"
        elif 0.8 <= r_obj < 1.25 and r_cpa < 0.75:
            cls = "capture"
        else:
            cls = "mixed"
        units_ok = 1 / 3 <= out[wt]["nn2d_median"] / out[ko]["nn2d_median"] <= 3
        return {"R_obj": r_obj, "R_cpo": r_cpo, "R_cpa": r_cpa, "class": cls if units_ok else "units differ",
                "nn_ratio_wt_over_ko": out[wt]["nn2d_median"] / out[ko]["nn2d_median"]}
    pairs = {"Gata4 rep1 vs WT E8.75": ratio("Gata4 KO rep1 (E8.75)", "WT E8.75"),
             "Gata4 rep2 vs WT E8.75": ratio("Gata4 KO rep2 (E8.75)", "WT E8.75"),
             "Mab21l2 vs WT E9.5": ratio("Mab21l2 KO (E9.5)", "WT E9.5")}
    absolute = {k: out[k]["counts_per_object_median"] / out["WT E8.75"]["counts_per_object_median"] for k in KO}
    return {"files": out, "pairs": pairs, "ko_depth_over_wt_e875": absolute}, inputs


# ---------------------------------------------------------------------------------------------------------------------
def d5d6(threads: int) -> tuple[dict, list[Path]]:
    import build_submission as bs
    import t3_decompose as t3d

    rules = bs.load_rules()
    dec = bs.load_decision(DECISION, rules)
    preset = dec["preset"]
    pools = rules["presets"][preset]["stage_pools"]
    import tempfile
    paths = {"wt": str(DATA / WT["E8.75"]), "targets": [],
             **{f"stage:{n}": str(DATA / "T2" / f"{n}.h5ad") for n in pools}}
    with bs.thread_limits(threads), tempfile.TemporaryDirectory() as tmp:
        ctx, _ = t3d.build_context(preset, paths, rules, n_jobs=threads, threads=threads, out=Path(tmp))
    V = ctx.V
    fc = rules["components"]["factor"]
    hp = rules["harmonisation"]
    mu = {"E8.75": ctx.wt.mu}
    cnt = {"E8.75": np.asarray(ctx.wt.count)}
    for sp in ctx.stages:
        mu[sp.name] = t3d.props_type_means(sp.pool.C, sp.pool.lab, V.L)
        cnt[sp.name] = np.asarray(sp.count)
    if "E8.5" not in mu:
        a = load_stage(DATA / WT["E8.5"], genes=ctx.genes, drop_duplicates=True)
        calls, _ = ctx.harm.calls(a, seed=hp["seed"])
        lab = t3d.calls_to_index(V, calls)
        C = t3d.dense_counts(a)
        mu["E8.5"] = t3d.props_type_means(C, lab, V.L)
        cnt["E8.5"] = np.bincount(lab[lab >= 0], minlength=V.L)
    need = ("E8.25_late", "E8.5", "E8.75", "E9.5")
    missing = [s for s in need if s not in mu]
    if missing:
        raise SystemExit(f"stages missing from the T3 context: {missing}")
    nm = np.asarray(V.nm, dtype=bool)

    def axis(late: str, early: str) -> np.ndarray:
        mask = nm & (cnt[late] >= 100) & (cnt[early] >= 100)
        phi, _ = t3d.gene_factor(mu[late], mu[early], V.p, mask, fc["eps"], tuple(fc["clip"]))
        v = np.log(phi)
        return v - v.mean()

    def logf(name: str) -> np.ndarray:
        phi, _ = bs.load_factor(RUN_DIR, name, list(ctx.genes), dec)
        v = np.log(phi)
        return v - v.mean()

    fG, fM = logf("gata4_pooled"), logf("transfer_mab21l2_all")
    fM95, fR = logf("transfer_mab21l2_vs95"), logf("reference_R")

    def decompose(d: np.ndarray) -> dict:
        dh = d / np.linalg.norm(d)
        aG, aM = float(fG @ dh), float(fM @ dh)
        rG, rM = fG - aG * dh, fM - aM * dh
        tot = float(fG @ fM)
        return {"inner_fG_fM": tot, "along_axis": aG * aM, "orthogonal": float(rG @ rM),
                "orthogonal_share": float(rG @ rM) / tot if tot > 0 else None,
                "corr_fG_fM": float(np.corrcoef(fG, fM)[0, 1]), "corr_fG_d": float(np.corrcoef(fG, d)[0, 1]),
                "corr_fM_d": float(np.corrcoef(fM, d)[0, 1]), "corr_rG_rM": float(np.corrcoef(rG, rM)[0, 1])}

    d_primary, d_secondary = axis("E9.5", "E8.25_late"), axis("E9.5", "E8.5")
    D5 = {"primary_axis": "E9.5 vs E8.25_late", "primary": decompose(d_primary),
          "secondary_axis": "E9.5 vs E8.5", "secondary": decompose(d_secondary),
          "context": {"corr_fM_fM95": float(np.corrcoef(fM, fM95)[0, 1]),
                      "corr_d_primary_with_log_all_over_vs95": float(np.corrcoef(d_primary, fM - fM95)[0, 1]),
                      "corr_fG_R": float(np.corrcoef(fG, fR)[0, 1])}}
    s = D5["primary"]["orthogonal_share"]
    D5["band"] = "low" if (s is None or s < 0.4) else ("middle" if s < 0.6 else "high")

    mask_all = nm.copy()
    for k in need:
        mask_all &= cnt[k] >= 100
    w = (np.asarray(V.p, dtype=np.float64) * mask_all)[:, None]

    def L(k: str) -> np.ndarray:
        x = np.log((w * mu[k]).sum(0) + fc["eps"])
        return x - x.mean()
    Ls = {k: L(k) for k in need}

    def place(k: str, a: str, b: str) -> float:
        v = Ls[b] - Ls[a]
        return float((Ls[k] - Ls[a]) @ v / (v @ v))
    s85 = place("E8.5", "E8.25_late", "E9.5")
    s875 = place("E8.75", "E8.25_late", "E9.5")
    s825_on_85_95 = place("E8.25_late", "E8.5", "E9.5")
    s875_on_85_95 = place("E8.75", "E8.5", "E9.5")
    check = bool(0 < s85 < 1 and s825_on_85_95 < 0)
    nominal = s85 + 0.25 * (1 - s85)
    if not check:
        pos = "check failed: uninformative"
    elif s875 < nominal - 0.1:
        pos = "young"
    elif s875 > nominal + 0.1:
        pos = "old"
    else:
        pos = "nominal"
    D6 = {"s_E8.5_on_E8.25late_E9.5": s85, "s_E8.75_on_E8.25late_E9.5": s875, "nominal_E8.75": nominal,
          "s_E8.25late_on_E8.5_E9.5": s825_on_85_95, "s_E8.75_on_E8.5_E9.5": s875_on_85_95,
          "known_lag_check_passed": check, "placement": pos, "labels_used": int(mask_all.sum())}
    shift = {"young": 1, "old": -1}.get(pos, 0)
    bands = ["low", "middle", "high"]
    final = bands[min(2, max(0, bands.index(D5["band"]) + shift))]
    p_anchor = {"low": 0.35, "middle": 0.5, "high": 0.65}[final]
    inputs = [DATA / WT[k] for k in need] + [DECISION] + [RUN_DIR / "factors" / f"{n}.csv" for n in
                                                          ("gata4_pooled", "transfer_mab21l2_all",
                                                           "transfer_mab21l2_vs95", "reference_R")]
    return {"D5": D5, "D6": D6, "p_anchor_band": final, "p_anchor": p_anchor}, inputs


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("which", choices=("d1", "d2", "d5d6"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res, inputs = {"d1": d1, "d2": d2, "d5d6": lambda: d5d6(a.threads)}[a.which]()
    res = {"check": a.which, "provenance": provenance(inputs), **res}
    (out / f"{a.which}.json").write_text(json.dumps(res, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps(res, indent=1, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
