#!/usr/bin/env python
"""One fetched outside-knockout dataset to its factor tables and one QC file (experiments/t3_ext/README.md).

  reference --wt WT_E8.75.h5ad --out D/reference/refstats.npz
      the WT E8.75 reference's panel genes (its gene order), lineage groups (obs['group'] of vec.prep.load_stage,
      'Unknown' left out), their shares and their mean within-cell proportions. Read by run and by combine.
  run ACC --dest D --wt WT_E8.75.h5ad --refstats D/reference/refstats.npz --out D/processed
      reads D/meta/ACC/plan.json and D/raw/ACC/MANIFEST.json (every file's sha256 checked). It touches ONLY the
      samples in datasets.yaml's pinned block: it refuses a dataset with no pins, a manifest fetched under other pins,
      a pinned accession absent from its own series' plan table, a pinned sample whose stored metadata no longer
      hashes to its pin, and screens each pinned sample's record again with t3x_screen BEFORE any count is read (a
      sample that does not pass refuses the run). Arms, genes and stages come from the pins only; no file of an
      unpinned sample is opened. It reads each pinned sample's counts restricted to the panel genes, labels every
      cell with a lineage group by transfer
      from the WT reference, and computes one log-ratio table per knockout gene with the fixed guards, minimum-cells
      rule and shrinkage below. It writes per gene <out>/<ACC>/<ACC>__<gene>.unit.json (log space, for combine) and
      its finished tables (pooled and per group), and <out>/<ACC>/QC.json. It prints no gene-level value.

The rules (constants in t3x_common.py, fixed 2026-11-20, none from a score):
  * Cells: stage capped at E9.5 (the T1 window applied to every cell); fewer than 200 detected genes is not a cell;
    a per-sample matrix of more than 60,000 barcodes is refused as unfiltered.
  * Groups: a logistic-regression transfer (vec.harmonise.fit_transfer, the frozen settings: 50 PCs, C 1.0, clip 10,
    5 depth draws, seed 0) trained on the reference's lineage GROUPS directly, so the outside call and the reference
    cell's own obs['group'] mean the same thing (the challenger's change 2). Confident: label_conf >= 0.6.
  * A stage drops if fewer than half of either arm's cells are confident. A (group, stage) drops if: fewer than 100
    confident cells in any arm; the group's share of the arm's cells differs across arms by more than 2x (a knockout
    that moves cells between calls is not read as expression change); the control arm's group profile does not have
    the reference group as its best Spearman match, ahead of the next by 0.05; or, where author labels and a committed
    dictionary exist, fewer than half the control cells' dictionary groups agree.
  * Each sample's confident and unconfident shares go into QC.json; a knockout gene's table is refused (status
    refused_low_confidence, nothing written) if any of its samples has a confident share under datasets.yaml's
    qc_floors.min_confident_share_per_sample.
  * Log ratio per gene on summed UMIs with a pseudocount of 10: log((x_ko+10)/N_ko) - log((x_ctl+10)/N_ctl); NaN (not
    measured) for a gene with fewer than 20 UMIs in any arm; clipped to log [0.25, 4]. Chimeras only with a
    wild-type-chimera second difference: [mut vs host] in knockout chimeras minus [mut vs host] in wild-type
    chimeras.
  * Shrinkage L = k * log ratio, k = min(n_h/(n_h+200), m/(m+1)) * m/(m+1) * d (t3x_common.shrink_k).
  * Per group, its passing stages are averaged with equal weight. No group passing = the dataset is left out
    (status 'left_out'), never given a pooled fallback.
"""

from __future__ import annotations

import gzip
import io
import json
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import fetch as fe  # noqa: E402  (experiments/t3_ext/fetch.py: the plan loader and the pin check)
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402

Refused = cm.Refused
CONTROL_ARMS = ("control", "ko_host", "wt_host")
PLAIN_ARMS = ("ko", "control")
CHIMERA_ARMS = ("ko_mut", "ko_host", "wt_mut", "wt_host")


# ---------------------------------------------------------------------------------------------------------------------
# Reading one sample's counts (no extraction into any folder)


def _read_text(raw: bytes, name: str) -> list[str]:
    if name.endswith(".gz"):
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", errors="replace").splitlines()


def _from_parts(mtx: bytes, mtx_name: str, feats: list[str], bcs: list[str]):
    data = gzip.decompress(mtx) if mtx_name.endswith(".gz") else mtx
    M = sparse.csr_matrix(spio.mmread(io.BytesIO(data), spmatrix=False).T)  # cells x features
    rows = [r.split("\t") for r in feats if r.strip()]
    symbols = [r[1] if len(r) > 1 else r[0] for r in rows]
    types = [r[2] if len(r) > 2 else "Gene Expression" for r in rows]
    barcodes = [b.split("\t")[0] for b in bcs if b.strip()]
    if M.shape != (len(barcodes), len(symbols)):
        raise Refused(f"matrix shape {M.shape} does not match {len(barcodes)} barcodes x {len(symbols)} features")
    return M, symbols, types, barcodes


def _pick(names: list[str], key: str) -> str:
    hits = [n for n in names if Path(n).name.lower().endswith((key, key + ".gz"))]
    if len(hits) != 1:
        raise Refused(f"expected exactly one '{key}' file, found {len(hits)}")
    return hits[0]


def read_sample(paths: list[Path]):
    """(cells x features CSR, symbols, feature types, barcodes) from a sample's mtx triple, its 10x h5, or a tar.gz
    holding a mtx triple (read member by member in memory)."""
    paths = [Path(p) for p in paths]
    h5 = [p for p in paths if p.name.endswith(".h5")]
    tars = [p for p in paths if p.name.endswith((".tar.gz", ".tgz", ".tar"))]
    if h5:
        if len(h5) != 1:
            raise Refused(f"{len(h5)} h5 files for one sample")
        import h5py

        with h5py.File(h5[0], "r") as f:
            g = f["matrix"]
            shape = tuple(g["shape"][()])
            M = sparse.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), shape=shape).T.tocsr()

            def dec(a):
                return [x.decode() if isinstance(x, bytes) else str(x) for x in a]

            symbols = dec(g["features"]["name"][()])
            types = dec(g["features"]["feature_type"][()])
            barcodes = dec(g["barcodes"][()])
        return M, symbols, types, barcodes
    if tars:
        if len(tars) != 1:
            raise Refused(f"{len(tars)} archives for one sample")
        with tarfile.open(tars[0], "r:*") as t:
            members = {m.name: m for m in t.getmembers() if m.isfile()}
            names = list(members)
            mtx = _pick(names, "matrix.mtx")
            feat = [n for n in names if Path(n).name.lower().startswith(("features.tsv", "genes.tsv"))
                    or Path(n).name.lower().endswith(("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"))]
            if len(feat) != 1:
                raise Refused(f"{tars[0].name}: expected one features/genes file, found {len(feat)}")
            bc = _pick(names, "barcodes.tsv")

            def rd(n):
                return t.extractfile(members[n]).read()

            return _from_parts(rd(mtx), mtx, _read_text(rd(feat[0]), feat[0]), _read_text(rd(bc), bc))
    names = [str(p) for p in paths]
    mtx = _pick(names, "matrix.mtx")
    feat = [n for n in names if "features.tsv" in Path(n).name.lower() or "genes.tsv" in Path(n).name.lower()]
    if len(feat) != 1:
        raise Refused(f"expected one features/genes file, found {len(feat)}")
    bc = _pick(names, "barcodes.tsv")
    return _from_parts(Path(mtx).read_bytes(), mtx, _read_text(Path(feat[0]).read_bytes(), feat[0]),
                       _read_text(Path(bc).read_bytes(), bc))


def to_panel(M, symbols: list[str], types: list[str], genes: list[str]) -> tuple[sparse.csr_matrix, np.ndarray, dict]:
    """Counts on the panel genes in the reference's order (absent columns zero), the presence mask and a record.
    Only gene-expression features; exact symbol first, then the committed alias table; duplicates summed."""
    gex = np.array([("gene" in t.lower() and "expression" in t.lower()) or t == "" for t in types])
    pos: dict[str, list[int]] = {}
    for i, s in enumerate(symbols):
        if gex[i]:
            pos.setdefault(s, []).append(i)
    rows, cols, how = [], [], {"exact": 0, "alias": 0, "duplicate_symbols": 0}
    for j, g in enumerate(genes):
        idx = pos.get(g)
        if idx:
            how["exact"] += 1
        elif g in cm.ALIASES and pos.get(cm.ALIASES[g]):
            idx = pos[cm.ALIASES[g]]
            how["alias"] += 1
        if idx:
            how["duplicate_symbols"] += int(len(idx) > 1)
            rows += idx
            cols += [j] * len(idx)
    Map = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(symbols), len(genes)))
    present = np.zeros(len(genes), bool)
    present[cols] = True
    C = sparse.csr_matrix(M @ Map).astype(np.int64)
    return C, present, {**how, "matched": int(present.sum()), "of": len(genes)}


def cell_filter(M, types) -> np.ndarray:
    gex = np.array([("gene" in t.lower() and "expression" in t.lower()) or t == "" for t in types])
    return np.asarray((M[:, np.flatnonzero(gex)] > 0).sum(1)).ravel() >= cm.MIN_GENES_PER_CELL


# ---------------------------------------------------------------------------------------------------------------------
# The log ratios (pure numpy; synthetic tests run these)


def _nanmean_rows(A: np.ndarray) -> np.ndarray:
    held = np.isfinite(A)
    n = held.sum(0)
    return np.where(n > 0, np.where(held, A, 0.0).sum(0) / np.maximum(n, 1), np.nan)


def _sum(C, rows) -> np.ndarray:
    return np.asarray(C[rows].sum(0)).ravel().astype(np.float64)


def log_ratio(x_a, x_b, present) -> np.ndarray:
    """log((x_a+a)/N_a) - log((x_b+a)/N_b) on present genes, clipped; NaN for an absent gene and for a gene with
    fewer than MIN_UMI in either arm (not measured: it holds no value in any mean and ends as factor 1 if no group
    or study holds it)."""
    x_a, x_b = np.asarray(x_a, float), np.asarray(x_b, float)
    Na, Nb = x_a[present].sum(), x_b[present].sum()
    out = np.full(len(x_a), np.nan)
    if Na <= 0 or Nb <= 0:
        return out
    L = np.log((x_a + cm.PSEUDOCOUNT) / Na) - np.log((x_b + cm.PSEUDOCOUNT) / Nb)
    L = np.where(np.minimum(x_a, x_b) < cm.MIN_UMI, np.nan, np.clip(L, *cm.LOG_CLIP))
    out[present] = L[present]
    return out


def chimera_log_ratio(x: dict[str, np.ndarray], present) -> np.ndarray:
    """[ko_mut vs ko_host] - [wt_mut vs wt_host], each with the pseudocount; NaN where any arm has < MIN_UMI."""
    a = log_ratio(x["ko_mut"], x["ko_host"], present)
    b = log_ratio(x["wt_mut"], x["wt_host"], present)
    low = np.min(np.vstack([x[k] for k in CHIMERA_ARMS]), 0) < cm.MIN_UMI
    L = np.where(low, np.nan, np.clip(a - b, *cm.LOG_CLIP))
    out = np.full(len(L), np.nan)
    out[present] = L[present]
    return out


def spearman_guard(profile, ref: dict, g: str, present) -> tuple[bool, dict]:
    """The control group's mean-proportion profile must match reference group g best, ahead of the next by the
    fixed margin (computed on control cells and the reference only: not a try)."""
    r = rankdata(np.asarray(profile)[present])
    rho = {}
    for h in ref["groups"]:
        q = rankdata(np.asarray(ref["mu"][h])[present])
        rho[h] = float(np.corrcoef(r, q)[0, 1]) if r.std() > 0 and q.std() > 0 else float("nan")
    order = sorted((v, h) for h, v in rho.items() if np.isfinite(v))[::-1]
    if not order or order[0][1] != g:
        return False, {"best": order[0][1] if order else None, "rho_own": rho.get(g)}
    margin = order[0][0] - (order[1][0] if len(order) > 1 else -1.0)
    return margin >= cm.SPEARMAN_MARGIN, {"best": g, "rho_own": rho[g], "margin": margin}


def unit_log_ratios(C, present, obs: pd.DataFrame, ref: dict, d: float, chimera: bool = False) -> tuple[dict, dict]:
    """One knockout gene's table: {group: L (P,) averaged over passing stages, NaN absent}, and its QC record.

    obs columns: arm (ko/control, or ko_mut/ko_host/wt_mut/wt_host for chimeras), stage (float), sample, call
    (reference group or ''), conf; optional author_group (the committed dictionary's group, or 'none')."""
    arms = CHIMERA_ARMS if chimera else PLAIN_ARMS
    if set(obs["arm"]) - set(arms):
        raise Refused(f"arms {sorted(set(obs['arm']) - set(arms))} are not of this design ({arms})")
    if chimera and set(arms) - set(obs["arm"]):
        raise Refused("a chimera needs knockout and wild-type chimeras (mut and host each): no second difference, "
                      "refused")
    C = sparse.csr_matrix(C)
    conf_ok = obs["conf"].to_numpy(float) >= cm.CONFIDENT
    per_stage: dict[float, dict] = {}
    qc = {"stages": {}, "groups": {}}
    for s in sorted(set(obs["stage"])):
        cm.stage_ok(f"E{s:g}")
        at = obs["stage"].to_numpy() == s
        if not all((at & (obs["arm"] == a).to_numpy()).any() for a in arms):
            qc["stages"][f"E{s:g}"] = {"skipped": "not every arm has cells"}
            continue
        tot = {a: int((at & (obs["arm"] == a).to_numpy()).sum()) for a in arms}
        cshare = {a: float((at & (obs["arm"] == a).to_numpy() & conf_ok).sum() / tot[a]) for a in arms}
        rec = {"cells": tot, "confident_share": cshare, "groups": {}}
        qc["stages"][f"E{s:g}"] = rec
        if min(cshare.values()) < cm.MIN_CONF_SHARE:
            rec["skipped"] = "fewer than half of an arm's cells called confidently"
            continue
        for g in ref["groups"]:
            sel = {a: np.flatnonzero(at & (obs["arm"] == a).to_numpy() & conf_ok & (obs["call"] == g).to_numpy())
                   for a in arms}
            n = {a: len(v) for a, v in sel.items()}
            share = {a: n[a] / tot[a] for a in arms}
            gr = {"cells": n, "call_share": share}
            rec["groups"][g] = gr
            if min(n.values()) < cm.MIN_CELLS:
                gr["skipped"] = f"fewer than {cm.MIN_CELLS} confident cells in an arm"
                continue
            if max(share.values()) > cm.SHARE_RATIO_MAX * min(share.values()):
                gr["skipped"] = "call share differs across arms by more than 2x (relabelling)"
                continue
            ctl = np.concatenate([sel[a] for a in arms if a in CONTROL_ARMS])
            Pc = C[ctl].toarray().astype(float)
            prof = (Pc / np.maximum(Pc.sum(1, keepdims=True), 1)).mean(0)
            ok, info = spearman_guard(prof, ref, g, present)
            gr["profile_guard"] = info
            if not ok:
                gr["skipped"] = "control profile does not match its reference group by the margin"
                continue
            if "author_group" in obs.columns:
                agree = float((obs["author_group"].to_numpy()[ctl] == g).mean())
                gr["author_agreement"] = agree
                if agree < cm.AUTHOR_AGREE:
                    gr["skipped"] = "author-label agreement under one half"
                    continue
            x = {a: _sum(C, sel[a]) for a in arms}
            L = chimera_log_ratio(x, present) if chimera else log_ratio(x["ko"], x["control"], present)
            m = min(len(set(obs["sample"].to_numpy()[sel[a]].astype(str))) for a in arms)
            k = cm.shrink_k(list(n.values()), m, d)
            gr.update(k=k, m=m, passed=True)
            per_stage.setdefault(s, {})[g] = k * L
    out = {}
    for g in ref["groups"]:
        vals = [per_stage[s][g] for s in per_stage if g in per_stage[s]]
        if vals:
            out[g] = _nanmean_rows(np.vstack(vals))
            qc["groups"][g] = {"stages": [f"E{s:g}" for s in per_stage if g in per_stage[s]]}
    qc["status"] = "ok" if out else "left_out"
    if not out:
        qc["reason"] = "no group passes the guards (fail closed: left out, no pooled fallback)"
    return out, qc


# ---------------------------------------------------------------------------------------------------------------------
# Classification (DGX only; synthetic tests pass their own calls)


def classify(wt_path: Path, C, present, genes: list[str], names: list[str], threads: int = 8):
    """Group calls and confidences for outside cells by transfer from the WT E8.75 reference's lineage groups."""
    import anndata as ad

    from vec import harmonise as hm
    from vec.prep import load_stage, log1p_cp10k

    pres = [g for g, p in zip(genes, present, strict=True) if p]
    ref = load_stage(wt_path, genes=genes, drop_duplicates=True)
    ref = ref[:, pres].copy()
    Cq = sparse.csr_matrix(C)[:, np.flatnonzero(present)]
    q = ad.AnnData(X=log1p_cp10k(Cq),
                   obs=pd.DataFrame(index=names), var=pd.DataFrame(index=pres))
    q.layers["counts"] = Cq.astype(np.int32)
    qd = float(np.median(np.asarray(Cq.sum(1)).ravel()))
    rd = hm.median_depth(ref, pres)
    depth = qd if (rd is not None and rd > qd * 1.05) else None
    model = hm.fit_transfer(ref, "group", genes=pres, classifier="logreg", n_pcs=50, C=1.0, clip=10.0,
                            threshold=cm.CONFIDENT, min_class_cells=10, exclude_labels=("Unknown",), depth=depth,
                            depth_tol=0.05, seed=0, n_jobs=threads)
    df = hm.predict(model, q, depth_draws=5, depth_tol=0.05, seed=0, proba=False)
    return df["label_raw"].astype(str).to_numpy(), df["label_conf"].to_numpy(float), {
        "reference_depth": rd, "query_depth": qd, "refit_depth": depth, "genes_used": len(pres)}


# ---------------------------------------------------------------------------------------------------------------------
# Orchestration


def write_unit(out_dir: Path, acc: str, gene: str, Lg: dict, qc: dict, ref: dict, meta: dict) -> dict:
    name = f"{acc}__{gene}".replace("+", "_and_").replace("/", "_")
    out_dir = cm.outside_repo(out_dir)
    unit = {"accession": acc, "gene": gene, "genes": ref["genes"], "refstats_sha256": ref["file_sha256"],
            "L": {g: [None if not np.isfinite(v) else float(v) for v in L] for g, L in Lg.items()}, **meta}
    up = out_dir / f"{name}.unit.json"
    cm.dump_json(unit, up)
    rec = {"unit_file": str(up), "unit_sha256": cm.sha256_file(up), "status": qc["status"]}
    if qc["status"] == "ok":
        phi, used = tb.finish_pooled(Lg, ref)
        rec["tables"] = tb.write_table(out_dir / "tables" / name, name, ref["genes"], phi,
                                       tb.finish_per_group(Lg, ref), {**meta, "groups_pooled": used})
    return rec


def run(acc: str, dest: Path, wt: Path, refstats: Path, out: Path, datasets: Path | None = None,
        threads: int = 8, classifier=None) -> dict:
    ds = cm.dataset(acc, datasets)
    floor = cm.qc_floors(datasets)["min_confident_share_per_sample"]
    if ds.get("cell_labels"):
        raise Refused(f"{acc}: per-cell genotype designs (cell-level labels inside one sample) are excluded in this "
                      "stage")
    if not ds.get("pinned"):
        raise Refused(f"{acc}: the pinned block is empty: nothing is processed")
    dest, out = cm.outside_repo(dest), cm.outside_repo(out)
    pl = fe._load_plan(acc, dest, datasets)
    mp = dest / "raw" / acc / "MANIFEST.json"
    if not mp.is_file():
        raise Refused(f"{mp}: no manifest; fetch first")
    man = json.loads(mp.read_text())
    if man.get("plan_sha256") != cm.sha256_file(dest / "meta" / acc / "plan.json"):
        raise Refused(f"{acc}: the manifest was not fetched from this plan")
    if man.get("pins_sha256") != cm.pins_sha256(ds["pinned"]) or man.get("pinned") != ds["pinned"]:
        raise Refused(f"{acc}: the samples, arms, genes or stages are not exactly as pinned when the files were "
                      "fetched: fetch again")
    # 1. the pins against the stored metadata, and the screen again on every pinned sample, before any count is read
    pinned = fe.verify_pins(ds, pl)
    units = cm.pinned_units(ds["pinned"])
    for s in pinned:
        cm.stage_ok(f"E{s['stage']:g}")
    screen = [s["accession"] for s in pinned]
    ref = tb.load_ref_stats(refstats)
    if ref["wt_sha256"] != cm.sha256_file(wt):
        raise Refused("the WT file is not the one the reference statistics were computed from")
    genes = ref["genes"]
    # 2. verify the pinned files, then read counts on the panel (no other sample's file is opened)
    blocks, obs_parts, presence, panel_rec, dropped = [], [], [], {}, {}
    for s in pinned:
        gsm = s["accession"]
        files = []
        for f in s["files"]:
            key = f"{gsm}/{f['name']}"
            p = dest / "raw" / acc / gsm / f["name"]
            if key not in man["files"] or not p.is_file() or cm.sha256_file(p) != man["files"][key]["sha256"]:
                raise Refused(f"{key}: missing or its sha256 differs from the manifest")
            files.append(p)
        M, symbols, types, bcs = read_sample(files)
        if M.shape[0] > cm.MAX_BARCODES:
            raise Refused(f"{gsm}: {M.shape[0]} barcodes: an unfiltered matrix, refused")
        keep = cell_filter(M, types)
        dropped[gsm] = int((~keep).sum())
        Cp, pres, prec = to_panel(M[np.flatnonzero(keep)], symbols, types, genes)
        panel_rec[gsm] = prec
        presence.append(pres)
        blocks.append(Cp)
        obs_parts.append(pd.DataFrame({"arm": s["arm"], "gene": ";".join(s["genes"]), "stage": s["stage"],
                                       "sample": gsm},
                                      index=[f"{gsm}:{b}" for b, k in zip(bcs, keep, strict=True) if k]))
    present = np.logical_and.reduce(presence)
    C = sparse.vstack(blocks).tocsr()
    obs = pd.concat(obs_parts)
    qc = {"accession": acc, "tier": ds["tier"], "study": ds["study"], "d": ds["d"], "screened_samples": screen,
          "pinned": ds["pinned"], "pins_sha256": cm.pins_sha256(ds["pinned"]),
          "screen_record": pl["screen"], "cells_dropped_low_genes": dropped, "panel_match": panel_rec,
          "panel_genes_present": int(present.sum()), "panel_genes": len(genes),
          "created": time.strftime("%Y-%m-%d %H:%M:%S %Z"), "git": _git_head()}
    adir = out / acc
    if present.mean() < cm.MIN_PANEL_MATCH:
        qc.update(status="failed", reason=f"only {present.sum()} of {len(genes)} panel genes matched")
        cm.dump_json(qc, adir / "QC.json")
        return qc
    # 3. group calls
    calls, conf, crec = (classifier or classify)(wt, C, present, genes, list(obs.index), threads)
    obs["call"], obs["conf"] = calls, conf
    qc["classifier"] = crec
    # 3b. the transfer has never run on real cells: each sample's confident share goes into QC.json, and a table
    # with a sample under the floor fixed in datasets.yaml is refused
    by_sample = {}
    for gsm, sub in obs.groupby("sample"):
        share = float((sub["conf"].to_numpy(float) >= cm.CONFIDENT).mean()) if len(sub) else 0.0
        by_sample[str(gsm)] = {"cells": int(len(sub)), "confident_share": share, "unconfident_share": 1.0 - share}
    qc["confidence_by_sample"] = by_sample
    qc["min_confident_share_per_sample"] = floor
    # 4. one table per knockout gene
    qc["units"] = {}
    for u in units:
        low = sorted(g for g in u["ko"] + u["control"] if by_sample.get(g, {}).get("confident_share", 0.0) < floor)
        if low:
            qc["units"][u["gene"]] = {"status": "refused_low_confidence", "samples_under_floor": low,
                                      "reason": f"confident share under {floor} in {low}: no table written"}
            continue
        rows = obs["sample"].isin(u["ko"] + u["control"]).to_numpy()
        sub = obs[rows].copy()
        sub["arm"] = np.where(sub["sample"].isin(u["ko"]), "ko", "control")
        Lg, uqc = unit_log_ratios(C[np.flatnonzero(rows)], present, sub, ref, ds["d"])
        cells = sub.groupby(["stage", "arm"]).size()
        uqc["cells_per_stage_and_arm"] = {f"E{s:g}:{a}": int(v) for (s, a), v in cells.items()}
        meta = {"tier": ds["tier"], "study": ds["study"], "d": ds["d"], "ko_samples": u["ko"],
                "control_samples": u["control"], "present": present.tolist(), "git": qc["git"]}
        uqc.update(write_unit(adir, acc, u["gene"], Lg, uqc, ref, meta))
        qc["units"][u["gene"]] = uqc
    qc["status"] = "ok" if any(v["status"] == "ok" for v in qc["units"].values()) else "left_out"
    cm.dump_json(qc, adir / "QC.json")
    return qc


def _git_head() -> dict:
    import subprocess

    def g(*a):
        try:
            return subprocess.run(["git", "-C", str(cm.REPO), *a], capture_output=True, text=True, timeout=30,
                                  env={"GIT_OPTIONAL_LOCKS": "0", "PATH": "/usr/bin:/bin"}).stdout.strip()
        except Exception:  # noqa: BLE001
            return "unknown"

    return {"commit": g("rev-parse", "HEAD"), "dirty": bool(g("status", "--porcelain", "--untracked-files=no"))}


def reference(wt: Path, out: Path) -> dict:
    from vec.prep import load_stage

    out = cm.outside_repo(out)
    a = load_stage(wt, drop_duplicates=True)
    st = tb.ref_stats_from_counts(a.layers["counts"], a.obs["group"].astype(str).to_numpy(), list(a.var_names))
    tb.save_ref_stats(st, out, cm.sha256_file(wt))
    return st


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="step", required=True)
    r = sub.add_parser("reference")
    r.add_argument("--wt", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    u = sub.add_parser("run")
    u.add_argument("accession")
    u.add_argument("--dest", type=Path, required=True)
    u.add_argument("--wt", type=Path, required=True)
    u.add_argument("--refstats", type=Path, required=True)
    u.add_argument("--out", type=Path, required=True)
    u.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    try:
        if a.step == "reference":
            st = reference(a.wt, a.out)
            print(f"reference: {len(st['genes'])} genes, {len(st['groups'])} groups -> {a.out}")
            return 0
        qc = run(a.accession, a.dest, a.wt, a.refstats, a.out, threads=a.threads)
    except Refused as e:
        print(f"refused: {e}")
        return 2
    print(f"{a.accession}: status {qc['status']}; panel genes present {qc['panel_genes_present']}/{qc['panel_genes']}")
    for gsm, v in (qc.get("confidence_by_sample") or {}).items():
        print(f"  {gsm}: {v['cells']} cells, confident share {v['confident_share']:.3f}")
    for gname, v in (qc.get("units") or {}).items():
        passed = sorted(v.get("groups", {}))
        print(f"  {gname}: {v['status']}; groups with a factor: {passed}")
    return 0 if qc["status"] != "failed" else 1


if __name__ == "__main__":
    sys.exit(main())
