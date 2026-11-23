#!/usr/bin/env python
"""Recomputed cell groups of the outside atlases, by the team's fixed marker scores, and the cross-atlas checks.

The labeller is experiments/external/feng_relabel.relabel, IMPORTED, never copied or edited, on the frozen base yaml
of row N54 (atlases.yaml relabel.base, its sha256 pinned there). Each atlas's config is derived in code
(derived_config): the base unchanged except its `stages` (the atlas's two bracket stages; relabel() refuses any other
stage) and, for an atlas whose early end is a whole embryo, the extra competing classes of atlases.yaml
extra_competing, whose genes the released E8.5 data decide (derive_extra, run once by `groups.py derive` before any
atlas file is opened; its output is committed and pinned). Adding a class adds its genes to every set's control
exclusions, so the controls are not those of the base: each derived config's sha256 goes to every manifest. A
competing class only removes cells from the target groups. The derived config keeps the base's terminal section, so
relabel() runs on it unchanged (its route output is ignored here: no builder route follows from an outside atlas).

The checks (atlases.yaml checks; all fixed before any atlas file is opened; each drops an atlas-group pair, never
relaxes a threshold):
  (a) profile_check: a group's mean log1p CP10k profile at the atlas's early end correlates (Spearman, genes with mean
      CP10k >= 1 in either profile) best with a released same-stage label of its own category;
  (b) the same at the late end against the atlas in use's recomputed E14.5 group profiles;
  (c) the share check is deleted with the composition block (the strict form, ruling row R27): no share is computed;
  purity_gate: for every competing class with released labels at the early stage, the group's mean class score lies
      below the midpoint between the released CM labels' and the class's own released labels' mean scores;
  embryo_gene_gate: in a whole-embryo early end, a gene enters the atlas's CM direction only if its mean CP10k in the
      CM group is above its mean in the same atlas's early non-heart cells (ratio > 1);
  min_cells: a group with fewer than min_cells_ext cells at either end has no direction.
Nothing here prints a gene-level value.

  uv run --frozen --no-sync python -I -B experiments/t1_ext/groups.py derive --released-e85 <E8.5_RNA.h5ad> --out F
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
for _p in (HERE, HERE.parents[1] / "experiments" / "external"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
import feng_relabel as fr  # noqa: E402
import t1x_common as cm  # noqa: E402

Refused = cm.Refused
CHAMBER = cm.REPO / "configs" / "t1_chamber_classes.yaml"
EXTRA_SOURCE = "released_e85_extra"
TARGET_LABELS = ("r_cm", "r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous", "r_epicardial")
CP10K = 1e4


# ------------------------------------------------------------------------------------------------ configs


def base_config(cfg: Mapping) -> dict:
    path = cm.REPO / cfg["relabel"]["base"]
    if cm.sha256_file(path) != cfg["relabel"]["base_sha256"]:
        raise Refused(f"{cfg['relabel']['base']} is not the frozen labeller of row N54 (sha256 differs)")
    return fr.load_config(path)


def load_extra(cfg: Mapping) -> dict:
    """The derived extra classes, from the committed file pinned in atlases.yaml."""
    ex = cfg["extra_competing"]
    f = cm.REPO / ex["derived_file"]
    if ex.get("derived_sha256") in (None, cm.PENDING) or not f.exists() or cm.sha256_file(f) != ex["derived_sha256"]:
        raise Refused("the extra competing classes are not derived, committed and pinned")
    return json.loads(f.read_text())


def derived_config(a: Mapping, cfg: Mapping, extra: Mapping | None, early: str, late: str) -> tuple[dict, str]:
    """The labeller's config for one atlas: the base with stages [early, late] and, for a whole-embryo early end, the
    extra classes. Checked by feng_relabel.check_config. Returns it and its sha256."""
    base = base_config(cfg)
    out = copy.deepcopy(base)
    out["stages"] = [early, late]
    if a.get("early_tissue") == "whole_embryo":
        if extra is None:
            raise Refused(f"{a.get('accession')}: a whole-embryo early end needs the derived extra classes")
        out["sources"][EXTRA_SOURCE] = (
            "Extra competing classes for whole-embryo early ends (experiments/t1_ext/atlases.yaml extra_competing): "
            "candidates proposed by the AI assistant from general knowledge, members decided by the released T1 E8.5 "
            "label means by the rule stated there.")
        for name, spec in extra["classes"].items():
            if name in out["sets"]:
                raise Refused(f"extra class {name} collides with a base set")
            if len(spec["genes"]) < int(out["score"]["min_genes_found"]):
                continue
            out["sets"][name] = {"role": "competing", "label": spec["label"],
                                 "released_labels": list(spec["released_labels"]),
                                 "genes": {g: EXTRA_SOURCE for g in spec["genes"]}}
    fr.check_config(out)
    return out, cm.sha256_json(out)


def derive_extra(label_means: Mapping[str, Mapping[str, float]], cfg: Mapping) -> dict:
    """atlases.yaml's extra-class rule on released per-label means (label -> gene -> mean log1p CP10k, one stage):
    a candidate joins its class if its highest mean over the class's released labels is >= min_mean and >= min_ratio
    times its highest mean over every other label; a gene the base yaml names never joins; a class with fewer than
    min_genes_found members is dropped."""
    ex = cfg["extra_competing"]
    base = base_config(cfg)
    named = fr.named_genes(base)
    need = int(base["score"]["min_genes_found"])
    out, record = {}, {}
    for name, spec in ex["classes"].items():
        own = [lab for lab in spec["released_labels"] if lab in label_means]
        if not own:
            record[name] = {"dropped": "none of its released labels is in the released file"}
            continue
        others = [lab for lab in label_means if lab not in spec["released_labels"]]
        genes, why = [], {}
        for g in spec["candidates"]:
            if g in named:
                why[g] = "named by the base yaml"
                continue
            hi = max(label_means[lab].get(g, 0.0) for lab in own)
            oth = max((label_means[lab].get(g, 0.0) for lab in others), default=0.0)
            if hi >= float(ex["min_mean"]) and hi >= float(ex["min_ratio"]) * oth:
                genes.append(g)
            else:
                why[g] = "fails the rule"
        record[name] = {"members": genes, "left_out": why}
        if len(genes) >= need:
            out[name] = {"label": spec["label"], "released_labels": own, "genes": genes}
        else:
            record[name]["dropped"] = f"fewer than {need} members"
    return {"stage": ex["stage"], "classes": out, "record": record, "rule": {k: ex[k] for k in ("min_mean",
                                                                                                    "min_ratio")}}


# ------------------------------------------------------------------------------------------------ released profiles


def chamber_categories() -> dict[str, str]:
    return dict(yaml.safe_load(CHAMBER.read_text())["labels"])


def label_profiles(C: sparse.spmatrix, genes: Sequence[str], labels: np.ndarray) -> dict:
    """Per label: cells, mean CP10k and mean log1p CP10k by gene symbol (duplicated symbols dropped)."""
    C = sparse.csr_matrix(C, dtype=np.float64)
    lib = np.asarray(C.sum(1)).ravel()
    if (lib <= 0).any():
        raise Refused("a released cell has no counts")
    cp = sparse.csr_matrix(sparse.diags(CP10K / lib) @ C)
    lg = cp.copy()
    lg.data = np.log1p(lg.data)
    uidx = fr.unique_index(list(genes))
    cols = np.array(sorted(uidx.values()))
    names = np.asarray(genes, dtype=str)[cols]
    out = {}
    for lab in sorted(set(labels)):
        rows = np.flatnonzero(labels == lab)
        out[lab] = {"n": int(rows.size), "genes": names,
                    "cp": np.asarray(cp[rows][:, cols].mean(0)).ravel(),
                    "log": np.asarray(lg[rows][:, cols].mean(0)).ravel()}
    return out


def released_stage(path: Path) -> tuple[sparse.csr_matrix, list[str], np.ndarray]:
    """A released T1 training stage: counts, gene names and the organisers' celltype column (DGX only)."""
    from vec.core import counts_of
    from vec.prep import load_stage

    ref = load_stage(path, lineage=False)
    return counts_of(ref), [str(g) for g in ref.var_names], ref.obs["celltype"].astype(str).to_numpy()


# ------------------------------------------------------------------------------------------------ checks


def _on(genes_from: Sequence[str], values: np.ndarray, genes_to: Sequence[str]) -> np.ndarray:
    pos = {g: i for i, g in enumerate(genes_from)}
    return np.array([values[pos[g]] if g in pos else np.nan for g in genes_to])


def profile_check(group: Mapping, refs: Mapping[str, Mapping], own: Sequence[str], min_cp: float) -> dict:
    """Spearman of the group's mean log1p CP10k with each reference profile, over shared genes with mean CP10k >= min_cp
    in either profile; passes when the best reference is one of `own`."""
    rho = {}
    for lab, r in refs.items():
        g_cp = _on(group["genes"], group["cp"], r["genes"])
        g_lg = _on(group["genes"], group["log"], r["genes"])
        ok = np.isfinite(g_cp) & ((g_cp >= min_cp) | (np.asarray(r["cp"]) >= min_cp))
        rho[lab] = float(spearmanr(g_lg[ok], np.asarray(r["log"])[ok])[0]) if ok.sum() >= 10 else float("nan")
    finite = {k: v for k, v in rho.items() if np.isfinite(v)}
    best = max(finite, key=finite.get) if finite else None
    return {"rho": rho, "best": best, "own": list(own), "passed": best is not None and best in own}


def purity_gate(group_mean: Mapping[str, float], released_cm: Mapping[str, float],
                released_own: Mapping[str, float | None]) -> dict:
    """Per competing class with a released mean of its own: passed when the group's mean class score is below the
    midpoint of the released CM labels' and the class's own labels' mean scores."""
    per, ok = {}, True
    for c, own in released_own.items():
        if own is None or c not in released_cm or c not in group_mean:
            per[c] = {"checked": False}
            continue
        mid = (released_cm[c] + own) / 2.0
        p = bool(group_mean[c] < mid)
        per[c] = {"checked": True, "passed": p}
        ok &= p
    return {"passed": ok, "classes": per}


def embryo_gene_gate(cm_cp: np.ndarray, nonheart_cp: np.ndarray) -> np.ndarray:
    """True where a gene may enter a whole-embryo atlas's CM direction: CM mean CP10k above the non-heart mean."""
    cm_cp, nh = np.asarray(cm_cp, dtype=np.float64), np.asarray(nonheart_cp, dtype=np.float64)
    return np.isfinite(cm_cp) & np.isfinite(nh) & (cm_cp > nh)


# ------------------------------------------------------------------------------------------------ the labelling


def label_atlas(red: Mapping, cfg_a: Mapping) -> dict:
    """relabel() on one reduced atlas, and the per-sample diagnostics of row 0 (label counts, the unassigned split into
    'no class reached min_score' and 'margin failed', the target classes' score quantiles). Deciding nothing."""
    out = fr.relabel(red["C"], red["genes"], red["stage"], cfg_a)
    labels, scores = out["labels"], out["scores"]
    classes = fr.classes(cfg_a)
    low = np.vstack([scores[n] for n in classes]).max(0) < float(cfg_a["assign"]["min_score"])
    un = cfg_a["assign"]["unassigned"]
    per = {}
    for g in sorted(set(red["gsm"])):
        m = red["gsm"] == g
        per[g] = {"stage": str(red["stage"][m][0]), "cells": int(m.sum()),
                  "labels": {lab: int(np.sum(m & (labels == lab))) for lab in fr.all_labels(cfg_a)},
                  "unassigned_no_class": int(np.sum(m & low)),
                  "unassigned_margin": int(np.sum(m & (labels == un) & ~low)),
                  "score_quantiles": {n: [float(q) for q in np.quantile(scores[n][m], [0.1, 0.5, 0.9])]
                                      for n, s in cfg_a["sets"].items() if s["role"] == "target"}}
    return {"labels": labels, "scores": scores, "per_sample": per, "group_cells": out["group_cells"],
            "label_counts": out["label_counts"], "unassigned_share": out["unassigned_share"],
            "sets": out["sets"]}


def class_means(scores: Mapping[str, np.ndarray], mask: np.ndarray, cfg_a: Mapping) -> dict[str, float]:
    return {n: float(np.mean(scores[n][mask])) for n, s in cfg_a["sets"].items()
            if s["role"] == "competing" and mask.any()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("derive", help="the extra classes from the released E8.5 (target-free; before any atlas file)")
    d.add_argument("--released-e85", required=True, type=Path)
    d.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    cfg = cm.load_atlases()
    C, genes, labels = released_stage(a.released_e85)
    prof = label_profiles(C, genes, labels)
    means = {lab: dict(zip(p["genes"], p["log"].tolist(), strict=True)) for lab, p in prof.items()}
    res = derive_extra(means, cfg)
    res["released_sha256"] = cm.sha256_file(a.released_e85)
    cm.dump_json(res, a.out)
    print(json.dumps({k: {"members": len(v.get("members", [])), "dropped": v.get("dropped")}
                      for k, v in res["record"].items()}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
