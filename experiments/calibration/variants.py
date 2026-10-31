#!/usr/bin/env python
"""Post-hoc diagnostic (2026-10-31): which pseudobulk_shift construction did the organisers score on a T2 board?

Added after calibrate.py's fixed constructions reproduced the embryo board's pseudobulk_shift row exactly (55.22) but
not the heart-extrapolation row (50.32 against 47.76; every anchor exact). On that board E8.75 and E9.5 use different
cell-type annotations (5 of 33 / 22 `celltype` labels shared; `cm_celltype` shares only "Unknown"), so the per-type
deltas depend on which labels define a type. This script scores, at seed 0 under the bake protocol with the published
anchors and the prefix build, the working-copy construction with each label choice below, and reports the score and
skills next to the published row. It is a search over the organisers' baseline construction, not over the scorer:
a variant counts as found only if it matches the published headline to its printed precision. Nothing here changes
calibrate.py's verdict rules; a found variant is added there as a construction with this provenance.

Usage (DGX, repo root): uv run python experiments/calibration/variants.py --board T2:heart:val_extrap \
    --out runs/calibration_variants
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np

from vec import baselines
from vec.boards import REPO_ROOT, load_boards
from vec.evaluate import Stage, score_board, subsample
from vec.io import read_stage

ALT_PREV = {"T2:heart:val_extrap": "T2/E8.25_late.h5ad"}  # the earlier heart stage, if E8.75 was not in its series
FILES = {"T2:heart:val_extrap": ("T2/E10.5.h5ad", "T2/E9.5.h5ad", "T2/E8.75.h5ad"),
         "T2:embryo:val_interp": ("T2/E7.5.h5ad", "T2/E7.25.h5ad", "T2/E6.75.h5ad"),
         "T1:val": ("T1/E10.5_RNA.h5ad", "T1/E9.5_RNA.h5ad", "T1/E8.5_RNA.h5ad")}


def labels(path: Path, how: str) -> np.ndarray:
    obs = ad.read_h5ad(path, backed="r").obs
    ct = obs["celltype"].astype(str).to_numpy()
    if how == "celltype":
        return ct
    if how == "cm_celltype":
        return obs["cm_celltype"].astype(str).to_numpy() if "cm_celltype" in obs else ct
    if how == "cm_or_celltype":  # cm_celltype where it names a cardiomyocyte state, else celltype
        if "cm_celltype" not in obs:
            return ct
        cm = obs["cm_celltype"].astype(str).to_numpy()
        return np.where(cm == "Unknown", ct, cm)
    if how == "single_type":  # one global delta for every cell
        return np.array(["all"] * len(ct))
    raise ValueError(how)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--board", choices=sorted(FILES), required=True)
    ap.add_argument("--data", type=Path, default=REPO_ROOT / "data" / "raw")
    ap.add_argument("--rows", type=Path, default=REPO_ROOT / "experiments" / "calibration" / "reference_rows.json")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "runs" / "calibration_variants")
    ap.add_argument("--round2", action="store_true", help="also vary the base cells and the earlier stage")
    a = ap.parse_args(argv)
    board = load_boards()[a.board]
    t, r, p = (a.data / f for f in FILES[a.board])

    def load(path: Path, how: str) -> Stage:
        s = Stage.from_anndata(read_stage(path, board.genes), board.genes, need_coords=board.needs_coords)
        return Stage(X=s.X, coords=s.coords, celltype=labels(path, how))

    pub = next(x for x in json.loads(a.rows.read_text())["rows"]
               if x["board"] == a.board and x["model"] == "pseudobulk_shift")
    target = load(t, "celltype")
    out = {"board": a.board, "published": pub, "variants": {}}
    # (label choice, prev file or None for the default, base the shift is applied to); round 2 added after round 1
    # (label choices on the working copies) found no match on the heart-extrapolation board
    variants = [(how, None, "wc") for how in ("celltype", "cm_celltype", "cm_or_celltype", "single_type")]
    if a.round2:
        alt_prev = ALT_PREV.get(a.board)
        variants += [("celltype", None, "full"), ("celltype", None, "wc_seed1"), ("celltype", None, "max_cells"),
                     ("celltype", None, "full_deltas_full_base")]
        if alt_prev:
            variants += [("celltype", alt_prev, "wc"), ("celltype", alt_prev, "full")]
    for how, prev_file, base in variants:
        ref = load(r, how)
        prev = load(a.data / prev_file if prev_file else p, how)
        if base == "wc":
            pred = baselines.pseudobulk_shift(subsample(prev, 0.1, 0), subsample(ref, 0.1, 0))
        elif base == "wc_seed1":
            pred = baselines.pseudobulk_shift(subsample(prev, 0.1, 1), subsample(ref, 0.1, 1))
        elif base in ("full", "full_deltas_full_base"):
            pred = baselines.pseudobulk_shift(prev, ref)
        elif base == "max_cells":
            full = baselines.pseudobulk_shift(prev, ref)
            pred = full.take(np.sort(np.random.default_rng(0).choice(full.n, min(full.n, board.max_cells),
                                                                       replace=False)))
        else:
            raise ValueError(base)
        res = score_board(board, pred, target, ref, seed=0, anchors="published", build="prefix", protocol="bake")
        skills = {m: round(100 * v, 2) for m, v in res.skills.items()}
        name = f"{how}|prev={prev_file or 'default'}|base={base}"
        out["variants"][name] = {"score": res.score, "skills": skills, "pred_cells": pred.n,
                                 "matches_headline": round(res.score, 2) == round(float(pub["score"]), 2)}
        print(f"{a.board} {name}: {res.score:.4f} (published {pub['score']}; {pred.n} cells) {skills}", flush=True)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / f"{a.board.replace(':', '__')}.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
