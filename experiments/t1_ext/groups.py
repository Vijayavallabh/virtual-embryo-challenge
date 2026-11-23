#!/usr/bin/env python
"""Recomputed cell groups of the outside atlases of row T72 (candidate MX2s; Task 1 only), by the team's fixed
marker scores.

The labeller is experiments/external/feng_relabel.relabel, IMPORTED, never copied or edited, on the frozen base yaml
of row N54 (atlases.yaml relabel.base, its sha256 pinned there): its sets and thresholds unchanged (log1p CP10k;
control sets of 25 bins and 10 genes per set gene, seed 20261119; min_genes_found 2; min_score 0.3, min_margin 0.25,
split_margin 0.1) and its label-to-group map. Each atlas's configs are derived in code (label_configs), the base
unchanged except:
  - `stages`: the stages of the cells the config labels;
  - the competing class skm (Tnnt3, Tnni2, Mylpf, Myl1, Tnnc2; label r_skm; source X010_T71, the skeletal-muscle
    competing class of row T71), on both new atlases. Its genes are set genes, so relabel() never draws one as a
    control gene (feng_relabel.relabel excludes every set gene);
  - on X020's E8.5 samples ONLY (an atlas whose early end is a whole embryo; atlases.yaml extra_competing.stage), the
    derived neural, paraxial and EXEM classes of `groups.py derive` (computed on the released T1 E8.5 file, a training
    stage, before any atlas file is opened; committed at experiments/t1_ext/results/extra_competing.json and pinned by
    its sha256 in atlases.yaml).
So X019 (heart at both ends) is labelled in one run over both stages (as row N54 labelled its two stages), and X020
in two runs: its E8.5 cells with skm and the derived classes, its E16.5 cells with skm only (the derived classes are
the E8.5 samples' only). Each config's sha256 goes to the tables' manifest. A competing class only removes cells from
the target groups. No author label, cluster or embedding is read; nothing here prints a gene-level value.

  uv run --frozen --no-sync python -I -B experiments/t1_ext/groups.py derive --released-e85 data/raw/T1/E8.5_RNA.h5ad \\
      --out experiments/t1_ext/results/extra_competing.json
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
for _p in (HERE, HERE.parents[1] / "experiments" / "external"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
import feng_relabel as fr  # noqa: E402
import t1x_common as cm  # noqa: E402

Refused = cm.Refused
EXTRA_SOURCE = "released_e85_extra"
SKM_SET = "skm"
SKM_LABEL = "r_skm"
SKM_GENES = ("Tnnt3", "Tnni2", "Mylpf", "Myl1", "Tnnc2")
SKM_SOURCE = "X010_T71"
SKM_PURITY_SCORE = 0.3  # (m2): a late-end r_cm cell whose skm score is at least this counts against purity
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


def derived_config(base: Mapping, stages: Sequence[str], extra: Mapping | None) -> tuple[dict, str]:
    """The base with `stages`, the skm class and, when `extra` is given, the derived classes. Checked by
    feng_relabel.check_config. Returns it and its sha256."""
    out = copy.deepcopy(dict(base))
    out["stages"] = list(stages)
    if SKM_SET in out["sets"]:
        raise Refused("the base already holds an skm set")
    out["sources"][SKM_SOURCE] = ("The skeletal-muscle competing class of row T71 (register row X010's amendment): "
                                  "gene symbols only.")
    out["sets"][SKM_SET] = {"role": "competing", "label": SKM_LABEL, "genes": {g: SKM_SOURCE for g in SKM_GENES}}
    if extra is not None:
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


def label_configs(a: Mapping, cfg: Mapping, extra: Mapping | None, early: str, late: str) -> list[dict]:
    """[{stages, config, sha256}]: one run over both ends, or, for an atlas whose early end is a whole embryo, the early
    end (with the derived classes, which must be of that stage) and the late end (without them)."""
    base = base_config(cfg)
    if a.get("early_tissue") != "whole_embryo":
        c, h = derived_config(base, [early, late], None)
        return [{"stages": [early, late], "config": c, "sha256": h}]
    if extra is None:
        raise Refused(f"{a.get('accession')}: a whole-embryo early end needs the derived extra classes")
    if early != cfg["extra_competing"]["stage"] or extra.get("stage") != early:
        raise Refused(f"{a.get('accession')}: the derived classes are of {extra.get('stage')}, not {early}")
    ce, he = derived_config(base, [early], extra)
    cl, hl = derived_config(base, [late], None)
    return [{"stages": [early], "config": ce, "sha256": he}, {"stages": [late], "config": cl, "sha256": hl}]


def derive_extra(label_means: Mapping[str, Mapping[str, float]], cfg: Mapping) -> dict:
    """atlases.yaml's extra-class rule on released per-label means (label -> gene -> mean log1p CP10k, one stage):
    a candidate joins its class if its highest mean over the class's released labels is >= min_mean and >= min_ratio
    times its highest mean over every other label; a gene the base yaml names never joins; a class with fewer than
    min_genes_found members is dropped."""
    ex = cfg["extra_competing"]
    base = base_config(cfg)
    named = fr.named_genes(base) | set(SKM_GENES)
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
                why[g] = "named by the base yaml or the skm class"
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


# ------------------------------------------------------------------------------------------------ the labelling


def label_atlas(red: Mapping, runs: Sequence[Mapping]) -> dict:
    """relabel() on one reduced atlas, one run per config (each on the cells of its stages). Returns per cell the
    label and the skm score, and per run its config's sha256, the sets found and the control genes' count."""
    stage = np.asarray(red["stage"], dtype=str)
    n = stage.size
    labels = np.full(n, "", dtype=object)
    skm = np.full(n, np.nan)
    rec = []
    done = np.zeros(n, dtype=bool)
    for r in runs:
        m = np.isin(stage, r["stages"])
        if not m.any():
            raise Refused(f"no cell at {r['stages']}")
        rows = np.flatnonzero(m)
        out = fr.relabel(sparse.csr_matrix(red["C"])[rows], red["genes"], stage[rows], r["config"])
        if set(out["sets"][SKM_SET]["controls"]) & set(SKM_GENES) or any(
                set(s["controls"]) & set(SKM_GENES) for s in out["sets"].values()):
            raise Refused("an skm gene was drawn as a control gene")
        labels[rows] = out["labels"]
        skm[rows] = out["scores"][SKM_SET]
        done[rows] = True
        rec.append({"stages": list(r["stages"]), "config_sha256": r["sha256"], "cells": int(rows.size),
                    "sets_found": {k: len(v["found"]) for k, v in out["sets"].items()},
                    "unassigned_share": out["unassigned_share"]})
    if not done.all():
        raise Refused("a cell was not labelled")
    return {"labels": labels.astype(str), "skm": skm, "runs": rec}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("derive", help="the extra classes from the released E8.5 (target-free; before any atlas file)")
    d.add_argument("--released-e85", required=True, type=Path)
    d.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    cfg = cm.load_atlases()
    want = (cm.REPO / cfg["extra_competing"]["derived_file"]).resolve()
    if a.out.resolve() != want:
        raise SystemExit(f"--out must be {cfg['extra_competing']['derived_file']} (the file atlases.yaml pins)")
    if Path(a.released_e85).name != f"{cfg['extra_competing']['stage']}_RNA.h5ad":
        raise SystemExit("--released-e85 must be the released E8.5 file (a training stage)")
    C, genes, labels = released_stage(a.released_e85)
    prof = label_profiles(C, genes, labels)
    means = {lab: dict(zip(p["genes"], p["log"].tolist(), strict=True)) for lab, p in prof.items()}
    res = derive_extra(means, cfg)
    res["released_sha256"] = cm.sha256_file(a.released_e85)
    want.parent.mkdir(parents=True, exist_ok=True)
    want.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: {"members": len(v.get("members", [])), "dropped": v.get("dropped")}
                      for k, v in res["record"].items()}, indent=1))
    print(f"wrote {cfg['extra_competing']['derived_file']} sha256 {cm.sha256_file(want)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
