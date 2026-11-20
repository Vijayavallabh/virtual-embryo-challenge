"""Factor tables in the frozen factor's form, and the reference statistics they are finished against.

A table is kept in LOG space while it is built (one value per lineage group and panel gene, NaN where the gene is
absent from the outside data or the group did not pass). It is finished, for every table, in one fixed order (the
challenger's change 2): exp, rescale so the reference's weighted panel total is unchanged (phi / (den . phi / sum den),
as t3_decompose.gene_factor does with its denominator), clip to [0.25, 4], and only then set absent genes to exactly 1.

Two forms are written per table:
  pooled     <out>/factors/<name>.csv with columns gene, factor in the WT E8.75 reference's gene order (the layout of
             transfer_mab21l2_all.csv) and <out>/factors/factors.json holding its sha256 under 'files'. The pooled log
             value is the reference-group-share-weighted mean over the groups that passed; den is the same
             share-weighted sum of the reference's group mean proportions over those groups.
  per group  <out>/<name>_per_group.csv, long form (group, gene, factor), one block per group that has its own
             factor; den is that group's reference mean proportions. A group without a row gets factor 1 (no change)
             in any builder, never the pooled factor.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import t3x_common as cm

Refused = cm.Refused


def finish(L, den) -> np.ndarray:
    """phi from log values L (NaN = absent) and the denominator profile den: exp, rescale, clip, absent = 1 last."""
    L = np.asarray(L, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    absent = ~np.isfinite(L)
    phi = np.exp(np.where(absent, 0.0, L))
    if den.sum() > 0:
        phi = phi / (den @ phi / den.sum())
    phi = np.clip(phi, *cm.CLIP)
    phi[absent] = 1.0
    return phi


def pool_groups(Lg: dict[str, np.ndarray], shares: dict[str, float]) -> tuple[np.ndarray, list[str]]:
    """Share-weighted mean over groups, per gene over the groups that hold a value; NaN where none does."""
    groups = [g for g in Lg if np.isfinite(Lg[g]).any()]
    if not groups:
        return np.full(len(next(iter(Lg.values()))) if Lg else 0, np.nan), []
    V = np.vstack([Lg[g] for g in groups])
    w = np.array([shares[g] for g in groups], dtype=np.float64)[:, None] * np.isfinite(V)
    num = np.nansum(np.where(np.isfinite(V), V, 0.0) * w, 0)
    tot = w.sum(0)
    return np.where(tot > 0, num / np.where(tot > 0, tot, 1.0), np.nan), groups


def finish_pooled(Lg: dict[str, np.ndarray], ref: dict) -> tuple[np.ndarray, list[str]]:
    L, groups = pool_groups(Lg, ref["shares"])
    if not groups:
        return np.ones(len(ref["genes"])), []
    den = sum(ref["shares"][g] * ref["mu"][g] for g in groups)
    return finish(L, den), groups


def finish_per_group(Lg: dict[str, np.ndarray], ref: dict) -> dict[str, np.ndarray]:
    return {g: finish(L, ref["mu"][g]) for g, L in Lg.items() if np.isfinite(L).any()}


def write_table(out_dir, name: str, genes: list[str], phi: np.ndarray, per_group: dict[str, np.ndarray],
                manifest: dict) -> dict:
    """Write both forms and the manifest outside the repository; returns the file records with sha256."""
    out_dir = cm.outside_repo(out_dir)
    fdir = out_dir / "factors"
    fdir.mkdir(parents=True, exist_ok=True)
    phi = np.asarray(phi, dtype=np.float64)
    if len(phi) != len(genes) or not np.all(np.isfinite(phi) & (phi > 0)):
        raise Refused(f"{name}: factor values must be finite, positive and one per reference gene")
    p = fdir / f"{name}.csv"
    pd.DataFrame({"gene": list(genes), "factor": phi}).to_csv(p, index=False)
    rows = [(g, gene, float(v)) for g in sorted(per_group) for gene, v in zip(genes, per_group[g], strict=True)]
    pg = out_dir / f"{name}_per_group.csv"
    pd.DataFrame(rows, columns=["group", "gene", "factor"]).to_csv(pg, index=False)
    files = {name: {"file": str(p), "sha256": cm.sha256_file(p)},
             f"{name}_per_group": {"file": str(pg), "sha256": cm.sha256_file(pg)}}
    fj = fdir / "factors.json"
    old = json.loads(fj.read_text()) if fj.is_file() else {"files": {}, "parameters": {}}
    old["files"].update({name: files[name]})
    old["note"] = ("outside-knockout factor tables (experiments/t3_ext); not a run of t3_decompose; the frozen "
                   "builder's decision file does not name them, so only a thin builder that checks these sha256 can "
                   "load them")
    cm.dump_json(old, fj)
    cm.dump_json({**manifest, "name": name, "files": files,
                  "pooled_summary": {"n_up": int((phi > 1).sum()), "n_down": int((phi < 1).sum()),
                                     "n_at_clip": int(((phi <= cm.CLIP[0]) | (phi >= cm.CLIP[1])).sum())},
                  "groups_with_own_factor": sorted(per_group)}, out_dir / f"{name}_manifest.json")
    return files


# ---------------------------------------------------------------------------------------------------------------------
# Reference statistics (WT E8.75): genes, lineage groups, their shares and mean within-cell proportions


def ref_stats_from_counts(C, groups_of_cells, genes: list[str]) -> dict:
    """Group shares and per-group mean within-cell proportions (t3_decompose.props_type_means) over the panel;
    'Unknown' is left out."""
    C = np.asarray(C.toarray() if hasattr(C, "toarray") else C, dtype=np.float64)
    g = np.asarray(groups_of_cells).astype(str)
    keep = g != "Unknown"
    names = sorted(set(g[keep]))
    P = C / np.maximum(C.sum(1, keepdims=True), 1)
    mu = {n: P[g == n].mean(0) for n in names}
    cnt = {n: int((g == n).sum()) for n in names}
    tot = sum(cnt.values())
    return {"genes": list(genes), "groups": names, "shares": {n: cnt[n] / tot for n in names}, "mu": mu,
            "count": cnt}


def save_ref_stats(st: dict, path, wt_sha256: str) -> None:
    path = cm.outside_repo(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, genes=np.array(st["genes"]), groups=np.array(st["groups"]),
             shares=np.array([st["shares"][g] for g in st["groups"]]),
             mu=np.vstack([st["mu"][g] for g in st["groups"]]),
             count=np.array([st["count"][g] for g in st["groups"]]), wt_sha256=np.array(wt_sha256))


def load_ref_stats(path) -> dict:
    z = np.load(path, allow_pickle=False)
    groups = [str(g) for g in z["groups"]]
    return {"genes": [str(g) for g in z["genes"]], "groups": groups,
            "shares": dict(zip(groups, map(float, z["shares"]), strict=True)),
            "mu": {g: z["mu"][i] for i, g in enumerate(groups)},
            "count": dict(zip(groups, map(int, z["count"]), strict=True)), "wt_sha256": str(z["wt_sha256"]),
            "file_sha256": cm.sha256_file(path)}
