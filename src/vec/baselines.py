"""The organisers' published reference baselines, re-implemented from their descriptions on
https://virtualembryo.ai/challenge/baselines. All inputs are `Stage`s already aligned to the same gene panel.
"""

from __future__ import annotations

import numpy as np

from vec.evaluate import Stage


def copy_last(last: Stage) -> Stage:
    """T1/T2 floor: the last observed stage, verbatim (coordinates included)."""
    return Stage(X=last.X.copy(), coords=None if last.coords is None else last.coords.copy())


def wt_identity(wt: Stage) -> Stage:
    """T3 floor: the matched wild type, unchanged."""
    return copy_last(wt)


def pseudobulk_shift(prev: Stage, last: Stage) -> Stage:
    """Per-cell-type constant velocity: add mean(last|c) - mean(prev|c) to each real last-stage cell of type c,
    clipped at zero. Types absent from `prev` are held fixed; coordinates are carried over unchanged."""
    if prev.celltype is None or last.celltype is None:
        raise ValueError("pseudobulk_shift needs cell-type labels on both stages")
    X = last.X.copy()
    prev_types = set(prev.celltype)
    for c in np.unique(last.celltype):
        if c not in prev_types:
            continue
        sel = last.celltype == c
        delta = last.X[sel].mean(0) - prev.X[prev.celltype == c].mean(0)
        X[sel] += delta
    return Stage(X=np.clip(X, 0.0, None), coords=None if last.coords is None else last.coords.copy())


def shift_transfer(ko_train: Stage, wt_train: Stage, wt_target: Stage) -> Stage:
    """T3: assume every knockout shares one response; add the observed KO's pseudobulk shift to the target WT."""
    delta = ko_train.X.mean(0) - wt_train.X.mean(0)
    X = np.clip(wt_target.X + delta, 0.0, None)
    return Stage(X=X, coords=None if wt_target.coords is None else wt_target.coords.copy())


def gene_ko(ko_train: Stage, wt_train: Stage, wt_target: Stage, genes: list[str], target_gene: str) -> Stage:
    """shift_transfer, then zero the knocked-out gene if it is on the panel (Ctnnb1 is not on the MERFISH panel)."""
    pred = shift_transfer(ko_train, wt_train, wt_target)
    if target_gene in genes:
        pred.X[:, genes.index(target_gene)] = 0.0
    return pred
