#!/usr/bin/env python
"""Validation anchor check for the two scorer builds. Pre-registered 2026-10-19, before the validation answers exist.

From 2026-10-20, with the released T2 validation answers, this reproduces the published shape anchors of the three T2
validation boards (d2_shape, occupancy_dice and scale_log_ratio, floor and ceiling) under both scorer builds of
vec.evaluate ("main" = pinned veckit git cc3859a, "prefix" = veckit 0.1.1) and classifies each published
occupancy_dice anchor as a proper-mode or a mirror-mode value.

Background. The published anchors predate veckit PR #9 (right-handed frame). Under the pre-fix build a Dice value is
the proper-mode value ("main") when the prediction's and the truth's SVD bases have the same parity, and the
mirror-mode value (the prediction against the truth's mirror image) otherwise. d2_shape and scale_log_ratio never
depend on the build, so they fingerprint the exact 10% draw, the half split and the floor prediction.

Compliance (docs/DATA_POLICY.md):
- Validation boards only. --board keys must be T2 validation boards (or synthetic "TEST:" keys made by `self-test`).
- The index must be the repo's boards/index.json, refreshed only through vec.boards.update_boards. The run stops,
  printing counts only, if the index sits under a "panels" directory or if any non-validation board in it still
  carries anchors (test-board anchors encode hidden targets and are never read). Nothing is fetched from the network.
- Published validation anchors are compared with recomputed values; they set no number in any prediction.
- When no validation answer file is present the script exits with status 0 before it opens the index.
- Outputs are aggregate numbers only (JSON and TXT). Run it on the DGX; bring back only those files.

Per board:
 1. Counts: truth_cells must be one half of the 10% working copy of the released target (first half n//2 or second
    half n - n//2), and ref_cells the 10% working copy of the reference stage.
 2. Exact-protocol search: plausible implementations of the 10% draws (numpy Generator / RandomState / permutation,
    separate or shared generators, a list of seeds, the reference drawn with the same seed or seed + 1, rows sorted
    or as drawn), of the split (veckit split_half or the emulator's sorted halves; seed 0 or the draw seed; truth as
    the first or second half), the floor prediction (the reference's working copy or the whole reference stage) and
    the metric seed (0 or the draw seed). A hit reproduces d2_shape (5 dp) and scale_log_ratio (4 dp) for both floor
    and ceiling exactly. For each hit the Dice anchors are recomputed under both builds and in mirror mode.
 3. Seed distribution under the emulator's protocol (vec.evaluate.subsample / split_half, metric seed 0) for the
    floor (working copy and whole stage) and the ceiling: Dice main / mirror / prefix, d2_shape, scale_log_ratio,
    and the z-score and percentile of each published value.
 4. Classification of each Dice anchor (rules fixed 2026-10-19):
    - with exact hits: "proper" if the published value equals the hit's main (proper-mode) Dice and not its
      mirror-mode Dice, "mirror" for the converse, "tie" if both are equal, "neither" if neither is; "ambiguous" if
      the hits disagree;
    - without a hit: from the seed distribution, "proper" if |z_proper| <= 2.5 and |z_mirror| >= 4, "mirror" for the
      converse, otherwise "undetermined"; "undetermined (fingerprints off)" if the published d2_shape or
      scale_log_ratio of that anchor sits more than 4 sd from its distribution. The floor uses whichever floor
      prediction (working copy or whole stage) fits the fingerprints better.
 5. Theory checks on every computed row: prefix Dice = the proper value when the parities agree and the mirror value
    otherwise, and a z-mirrored prediction scores exactly like the original under prefix.

What the result can say (screened research, 2026-10-18): a Dice anchor at a mirror-mode value shows the bake ran a
pre-fix build; whether the server runs it today is a separate question (leaderboard consistency, organiser answer).
With --index-before (an older repo copy of boards/index.json, e.g. from git), a change in the Dice anchors without a
change in d2_shape points to a re-bake with another canonicalisation.

Usage (DGX, repo root, from 2026-10-20; long runs in tmux with a log under runs/logs/):
  uv run python experiments/scorer_build/anchor_check.py check --out runs/scorer_build/anchor_check
  uv run python experiments/scorer_build/anchor_check.py check \
      --board T2:heart:val_extrap PATH_TO_E10.5_ANSWERS data/raw/T2/E9.5.h5ad
  uv run python experiments/scorer_build/anchor_check.py self-test --target data/raw/T2/E9.5.h5ad \
      --ref data/raw/T2/E8.75.h5ad --out runs/scorer_build/anchor_selftest
Answer files are found under --data (default data/raw) as T2/.../E8.5*.h5ad, E10.5*.h5ad and E7.5*.h5ad; pass
--board KEY TARGET REF when the release names them differently.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from vec.boards import REPO_ROOT, VALIDATION_BOARDS, VALIDATION_SPLITS
from vec.evaluate import (
    SUBSAMPLE_FRAC,
    Stage,
    handedness,
    mirror_z,
    scorer_build,
    split_half,
    subsample,
)

SHAPE_METRICS = ("d2_shape", "occupancy_dice", "scale_log_ratio")
DECIMALS = {"d2_shape": 5, "occupancy_dice": 4, "scale_log_ratio": 4}  # as score_task2_v2 rounds them
# board -> (released stage of the answers, reference stage of copy_last); ref_cells in boards/index.json confirm the
# reference (5872 = 10% of E8.25_late, 5374 = 10% of E9.5, 1330 = 10% of E7.25)
T2_BOARDS = {
    "T2:heart:val_interp": ("E8.5", "E8.25_late"),
    "T2:heart:val_extrap": ("E10.5", "E9.5"),
    "T2:embryo:val_interp": ("E7.5", "E7.25"),
}
T2_VALIDATION_BOARDS = frozenset(k for k in VALIDATION_BOARDS if k.startswith("T2:"))
assert T2_VALIDATION_BOARDS == frozenset(T2_BOARDS), "keep T2_BOARDS in sync with vec.boards.VALIDATION_BOARDS"
SYNTHETIC_PREFIX = "TEST:"
DRAW_KINDS = ("gen", "gen_noshuffle", "legacy", "perm")
SHARED_KINDS = ("gen", "legacy")
DEFAULT_SEARCH_SEEDS = (0, 1, 2, 3, 7, 42, 123, 2024, 2025, 2026, 12345)
MAX_HITS = 20
# distribution quantity -> the published metric it is compared with
METRIC_OF = {"dice_main": "occupancy_dice", "dice_mirror": "occupancy_dice", "dice_prefix": "occupancy_dice",
             "d2_shape": "d2_shape", "scale_log_ratio": "scale_log_ratio"}
# pre-registered classification thresholds (2026-10-19)
Z_INSIDE, Z_OUTSIDE, Z_FINGERPRINT = 2.5, 4.0, 4.0


# --------------------------------------------------------------------------------------------- guards and inputs
def check_board_keys(keys) -> None:
    """Only T2 validation boards or synthetic TEST: boards may be analysed."""
    for k in keys:
        if not (k in T2_VALIDATION_BOARDS or str(k).startswith(SYNTHETIC_PREFIX)):
            raise SystemExit(f"--board {k!r} is not a T2 validation board (or a synthetic TEST: board): refused.")


def screen_index(index: dict, keys, index_path) -> None:
    """Stop before any anchor is used unless the index is a stripped repo index and every requested board is a T2
    validation board (or synthetic). Prints counts only, never values."""
    if "panels" in Path(index_path).resolve().parts:
        raise SystemExit("refusing an index under a 'panels' directory: the public board file may carry test-board "
                         "anchors. Use the repo's boards/index.json (refreshed only via vec.boards.update_boards).")
    unstripped = sum(1 for k, e in index.items()
                     if isinstance(e, dict) and e.get("anchors") and not str(k).startswith(SYNTHETIC_PREFIX)
                     and not (k in VALIDATION_BOARDS and e.get("split") in VALIDATION_SPLITS))
    if unstripped:
        raise SystemExit(f"the index carries anchors for {unstripped} non-validation board(s), so it is not the "
                         "stripped repo index. Refresh boards only via vec.boards.update_boards. Nothing was printed.")
    check_board_keys(keys)
    for k in keys:
        e = index.get(k)
        if not isinstance(e, dict):
            raise SystemExit(f"board {k!r} is not in the index")
        if not str(k).startswith(SYNTHETIC_PREFIX) and e.get("split") not in VALIDATION_SPLITS:
            raise SystemExit(f"board {k!r} is not marked as a validation split in the index: refused.")
        missing = [m for m in SHAPE_METRICS if m not in (e.get("anchors") or {})]
        if missing:
            raise SystemExit(f"board {k!r} has no published anchors for {missing}")


def load_index_entries(index_path, keys) -> dict:
    """The requested boards' entries from a screened index; only their shape anchors and counts are kept."""
    path = Path(index_path)
    if "panels" in path.resolve().parts:  # checked before the file is opened
        screen_index({}, keys, path)
    index = json.loads(path.read_text())
    screen_index(index, keys, path)
    return {k: {"truth_cells": index[k].get("truth_cells"), "ref_cells": index[k].get("ref_cells"),
                "anchors": {m: dict(index[k]["anchors"][m]) for m in SHAPE_METRICS}} for k in keys}


def find_answer(data_dir: Path, stage: str) -> Path | None:
    """The released answer file for a T2 validation stage (name starting with the stage, then no digit)."""
    pat = re.compile(rf"^{re.escape(stage)}(?![0-9])")
    hits = sorted(p for p in Path(data_dir).rglob("*.h5ad") if "T2" in p.parts and pat.match(p.name))
    if len(hits) > 1:
        raise SystemExit(f"several candidate answer files for {stage}: {[str(p) for p in hits]}; pass --board")
    return hits[0] if hits else None


def load_coords(path) -> np.ndarray:
    """obsm['spatial_3D'][:, :3] as float32 (how the scorer reads submitted coordinates)."""
    import h5py

    with h5py.File(path, "r") as f:
        o = f["obsm"]["spatial_3D"]
        if isinstance(o, h5py.Dataset):
            C = o[()]
        else:
            C = np.stack([o[c][()] for c in list(o.attrs.get("column-order"))], 1)
    return np.asarray(C, dtype=np.float32)[:, :3]


# --------------------------------------------------------------------------------------------- metric helpers
def _sm():
    import common.shape_metrics as sm

    return sm


def d2(P, A, metric_seed: int) -> float:
    return _sm().d2_distance(P, A, seed=metric_seed)


def scale(P, A) -> float:
    return _sm().scale_log_ratio(P, A)


def matches(value: float, published, decimals: int) -> bool:
    """Equal to the published value at the panel's rounding (or exactly, if the published value is unrounded)."""
    if published is None or not np.isfinite(value):
        return False
    return abs(round(value, decimals) - published) < 1e-9 or abs(value - published) < 1e-9


def dice_row(P, A, metric_seed: int = 0) -> dict:
    """Dice of prediction P against truth A under both builds and in mirror mode, plus the theory checks."""
    h = handedness(P, A, seed=metric_seed)
    with scorer_build("prefix"):
        pre = _sm().occupancy_dice(P, A, seed=metric_seed)[0]
        pre_twin = _sm().occupancy_dice(mirror_z(P), A, seed=metric_seed)[0]
    return {"dice_main": h["dice_proper"], "dice_mirror": h["dice_mirror"], "dice_prefix": round(pre, 4),
            "parity_pred": h["parity_pred"], "parity_truth": h["parity_truth"], "prefix_mode": h["prefix_mode"],
            "rule_ok": round(pre, 4) == h["dice_prefix"], "twin_ok": abs(pre_twin - pre) < 1e-12}


def shape_row(P, A, metric_seed: int = 0, both_builds: bool = False) -> dict:
    """All three shape metrics of P against A. both_builds: also compute d2_shape and scale_log_ratio under the
    prefix build and record whether they agree with main (they never call _canonicalise)."""
    row = dice_row(P, A, metric_seed)
    vals = {}
    for build in ("main", "prefix") if both_builds else ("main",):
        with scorer_build(build):
            vals[build] = (round(d2(P, A, metric_seed), 5), round(scale(P, A), 4))
    row["d2_shape"], row["scale_log_ratio"] = vals["main"]
    if both_builds:
        row["d2_scale_build_independent"] = vals["main"] == vals["prefix"]
    return row


# --------------------------------------------------------------------------------------------- candidate protocols
def _draw(kind: str, seed: int, N: int, n: int) -> np.ndarray:
    if kind == "gen":
        return np.random.default_rng(seed).choice(N, n, replace=False)
    if kind == "gen_noshuffle":
        return np.random.default_rng(seed).choice(N, n, replace=False, shuffle=False)
    if kind == "legacy":
        return np.random.RandomState(seed).choice(N, n, replace=False)
    if kind == "perm":
        return np.random.default_rng(seed).permutation(N)[:n]
    raise ValueError(kind)


def _second_draw(kind: str, seed: int, first: tuple[int, int], second: tuple[int, int]) -> np.ndarray:
    """The second of two draws from one generator (a shared-generator implementation)."""
    g = np.random.default_rng(seed) if kind == "gen" else np.random.RandomState(seed)
    g.choice(first[0], first[1], replace=False)
    return g.choice(second[0], second[1], replace=False)


def target_draws(Nt: int, nts, Nr: int, nrs, seeds, truth_cells) -> dict[str, tuple[int, np.ndarray]]:
    """label -> (draw seed, row indices) for the target's 10% working copy."""
    out = {}
    for nt in nts:
        for s in seeds:
            for kind in DRAW_KINDS:
                out[f"{kind}(seed={s}, n={nt})"] = (s, _draw(kind, s, Nt, nt))
            for kind in SHARED_KINDS:
                for nr in nrs:
                    out[f"{kind}-shared(seed={s}, n={nt}, reference drawn first n={nr})"] = (
                        s, _second_draw(kind, s, (Nr, nr), (Nt, nt)))
    if truth_cells and Nt in (2 * truth_cells, 2 * truth_cells + 1):
        out["released file is the working copy"] = (0, np.arange(Nt))
    return out


def ref_draws(Nt: int, nts, Nr: int, nrs, seeds) -> dict[str, np.ndarray | None]:
    """label -> row indices for the reference's 10% working copy (None = the whole reference stage)."""
    out: dict[str, np.ndarray | None] = {"whole reference stage": None}
    for nr in nrs:
        for s in sorted(set(seeds) | {s + 1 for s in seeds}):
            for kind in DRAW_KINDS:
                out[f"{kind}(seed={s}, n={nr})"] = _draw(kind, s, Nr, nr)
        for s in seeds:
            for kind in SHARED_KINDS:
                for nt in nts:
                    out[f"{kind}-shared(seed={s}, n={nr}, target drawn first n={nt})"] = _second_draw(
                        kind, s, (Nt, nt), (Nr, nr))
    return out


def _halves(n: int, method: str, seed: int) -> tuple[np.ndarray, np.ndarray]:
    from common.core_metrics import split_half as veckit_split_half

    a, b = veckit_split_half(n, seed=seed)
    return (a, b) if method == "veckit" else (np.sort(a), np.sort(b))


@dataclass(frozen=True)
class Protocol:
    """One complete implementation of the server's anchor bake."""

    target_draw: str  # a target_draws label
    target_sorted: bool
    split: str  # "veckit" (unsorted halves) or "emulator" (sorted halves)
    split_seed: int
    truth_first: bool
    floor: str  # a ref_draws label, plus ", sorted" / ", as drawn" for working copies
    metric_seed: int


def realise(protocol: Protocol, T, R, tdraws: dict, rdraws: dict):
    """(truth A, ceiling prediction B, floor prediction P) for a protocol."""
    _, ti = tdraws[protocol.target_draw]
    T10 = T[np.sort(ti) if protocol.target_sorted else ti]
    ia, ib = _halves(len(T10), protocol.split, protocol.split_seed)
    A, B = (T10[ia], T10[ib]) if protocol.truth_first else (T10[ib], T10[ia])
    if protocol.floor == "whole reference stage":
        P = R
    else:
        label, order = protocol.floor.rsplit(", ", 1)
        ri = rdraws[label]
        P = R[np.sort(ri) if order == "sorted" else ri]
    return A, B, P


def _counts(N: int, frac: float) -> list[int]:
    """Possible 10% working-copy sizes of an N-cell stage (rounding conventions differ)."""
    return sorted({math.floor(frac * N), math.ceil(frac * N), int(round(frac * N))})


def exact_search(T, R, published: dict, truth_cells, ref_cells, seeds, frac: float = SUBSAMPLE_FRAC) -> dict:
    """Protocols that reproduce the published d2_shape and scale_log_ratio anchors (floor and ceiling) exactly."""
    Nt, Nr = len(T), len(R)
    nts = [n for n in _counts(Nt, frac) if truth_cells is None or truth_cells in (n // 2, n - n // 2)]
    nrs = [n for n in _counts(Nr, frac) if n == ref_cells] or _counts(Nr, frac)
    tdraws = target_draws(Nt, nts, Nr, nrs, seeds, truth_cells)
    rdraws = ref_draws(Nt, nts, Nr, nrs, seeds)
    pub_d2, pub_sc = published["d2_shape"], published["scale_log_ratio"]
    # stage 1: truth half and ceiling (scale first: cheap and independent of row order)
    halves, n_tried, n_scale = [], 0, 0
    for t_label, (s, ti) in tdraws.items():
        for t_sorted in (False, True):
            T10 = T[np.sort(ti) if t_sorted else ti]
            for split in ("veckit", "emulator"):
                for split_seed in sorted({0, s}):
                    ia, ib = _halves(len(T10), split, split_seed)
                    for truth_first in (True, False):
                        A, B = (T10[ia], T10[ib]) if truth_first else (T10[ib], T10[ia])
                        if truth_cells is not None and len(A) != truth_cells:
                            continue
                        n_tried += 1
                        if not matches(scale(B, A), pub_sc["ceiling"], DECIMALS["scale_log_ratio"]):
                            continue
                        n_scale += 1
                        for ms in sorted({0, s}):
                            if matches(d2(B, A, ms), pub_d2["ceiling"], DECIMALS["d2_shape"]):
                                halves.append((t_label, t_sorted, split, split_seed, truth_first, ms, A))
    # stage 2: floor prediction for each surviving truth half
    hits = []
    for t_label, t_sorted, split, split_seed, truth_first, ms, A in halves:
        for r_label, ri in rdraws.items():
            if ri is None:
                floors = [("whole reference stage", R)]
            else:
                floors = [(f"{r_label}, sorted", R[np.sort(ri)]), (f"{r_label}, as drawn", R[ri])]
            if not matches(scale(floors[0][1], A), pub_sc["floor"], DECIMALS["scale_log_ratio"]):
                continue
            for f_label, P in floors:
                if matches(d2(P, A, ms), pub_d2["floor"], DECIMALS["d2_shape"]):
                    hits.append(Protocol(t_label, t_sorted, split, split_seed, truth_first, f_label, ms))
    return {"hits": hits, "n_truth_halves_tried": n_tried, "n_ceiling_scale_matches": n_scale,
            "n_ceiling_matches": len(halves), "n_target_draws": len(tdraws), "n_ref_draws": len(rdraws),
            "working_copy_sizes": {"target": nts, "reference": nrs}, "_tdraws": tdraws, "_rdraws": rdraws}


# --------------------------------------------------------------------------------------------- distribution
def _stage(C) -> Stage:
    return Stage(X=np.zeros((len(C), 1), np.float32), coords=np.asarray(C, np.float32))


def distribution(T, R, n_seeds: int, frac: float = SUBSAMPLE_FRAC, metric_seed: int = 0,
                 ref_frac: float | None = None) -> list[dict]:
    """Shape rows over seeds of the emulator's protocol (vec.evaluate.score_board's working copies and split).
    frac applies to the target (1.0 when the released file is already the working copy), ref_frac to the reference."""
    tgt, ref = _stage(T), _stage(R)
    rows = []
    for s in range(n_seeds):
        truth, other = split_half(subsample(tgt, frac, s), seed=0)
        ref_w = subsample(ref, frac if ref_frac is None else ref_frac, s + 1)
        A = truth.coords
        rows.append({"seed": s, "floor": shape_row(ref_w.coords, A, metric_seed),
                     "floor_whole_stage": shape_row(R, A, metric_seed),
                     "ceiling": shape_row(other.coords, A, metric_seed)})
    return rows


def where(value, xs) -> dict:
    xs = np.asarray(xs, float)
    m, sd = float(xs.mean()), (float(xs.std(ddof=1)) if len(xs) > 1 else 0.0)
    out = {"mean": round(m, 5), "sd": round(sd, 5), "min": round(float(xs.min()), 5), "max": round(float(xs.max()), 5)}
    if value is not None and np.isfinite(value):
        out["z_of_published"] = round((value - m) / sd, 2) if sd > 0 else (0.0 if abs(value - m) < 1e-12 else None)
        out["pct_of_published"] = round(float((xs <= value).mean()) * 100, 1)
    return out


def _absz(d: dict) -> float:
    z = d.get("z_of_published")
    return float("inf") if z is None else abs(z)


# --------------------------------------------------------------------------------------------- classification
def classify_exact(pub: float, rows: list[dict]) -> str:
    labels = set()
    for r in rows:
        eq_p = matches(r["dice_main"], pub, DECIMALS["occupancy_dice"])
        eq_m = matches(r["dice_mirror"], pub, DECIMALS["occupancy_dice"])
        labels.add("tie" if eq_p and eq_m else "proper" if eq_p else "mirror" if eq_m else "neither")
    return labels.pop() if len(labels) == 1 else "ambiguous"


def classify_distribution(pub_dice: float, stats: dict) -> str:
    fp = max(_absz(stats["d2_shape"]), _absz(stats["scale_log_ratio"]))
    if fp > Z_FINGERPRINT:
        return "undetermined (fingerprints off)"
    zp, zm = _absz(stats["dice_main"]), _absz(stats["dice_mirror"])
    if zp <= Z_INSIDE and zm >= Z_OUTSIDE:
        return "proper"
    if zm <= Z_INSIDE and zp >= Z_OUTSIDE:
        return "mirror"
    return "undetermined"


def check_board(key: str, T, R, entry: dict, seeds=DEFAULT_SEARCH_SEEDS, n_seeds: int = 50,
                frac: float = SUBSAMPLE_FRAC) -> dict:
    """Steps 1-5 for one board. entry: {"truth_cells", "ref_cells", "anchors": {metric: {floor, ceiling}}}."""
    pub = entry["anchors"]
    tc, rc = entry.get("truth_cells"), entry.get("ref_cells")
    Nt, Nr = len(T), len(R)
    nts = _counts(Nt, frac)
    out = {"board": key, "N_target": Nt, "N_ref": Nr, "index_truth_cells": tc, "index_ref_cells": rc,
           "counts": {"target_working_copy": nts,
                      "truth_is_first_half": [n for n in nts if tc == n // 2],
                      "truth_is_second_half": [n for n in nts if tc == n - n // 2],
                      "ref_cells_match": rc in _counts(Nr, frac) if rc is not None else None,
                      "released_file_is_working_copy": bool(tc) and Nt in (2 * tc, 2 * tc + 1),
                      "released_file_is_truth_half": bool(tc) and Nt == tc},
           "published": {m: pub[m] for m in SHAPE_METRICS}}
    t0 = time.time()
    found = exact_search(T, R, pub, tc, rc, seeds, frac)
    tdraws, rdraws = found.pop("_tdraws"), found.pop("_rdraws")
    hits = found["hits"]
    # identical arrays reached through different labels (e.g. choice with and without shuffle, then sorted) are one hit
    unique: dict[bytes, list[Protocol]] = {}
    for p in hits:
        A, B, P = realise(p, T, R, tdraws, rdraws)
        digest = hashlib.sha256(b"".join(np.ascontiguousarray(x).tobytes() for x in (A, B, P))
                                + str(p.metric_seed).encode()).digest()
        unique.setdefault(digest, []).append(p)
    hit_rows = []
    for group in list(unique.values())[:MAX_HITS]:
        p = group[0]
        A, B, P = realise(p, T, R, tdraws, rdraws)
        hit_rows.append({"protocol": asdict(p), "equivalent": [asdict(q) for q in group],
                         "floor": shape_row(P, A, p.metric_seed, both_builds=True),
                         "ceiling": shape_row(B, A, p.metric_seed, both_builds=True)})
    out["exact_search"] = {k: v for k, v in found.items() if k != "hits"}
    out["exact_search"]["n_hits"] = len(unique)
    out["exact_search"]["n_hit_labels"] = len(hits)
    out["exact_search"]["sec"] = round(time.time() - t0, 1)
    out["exact_hits"] = hit_rows

    target_frac = frac
    if out["counts"]["released_file_is_working_copy"] and not (out["counts"]["truth_is_first_half"]
                                                               or out["counts"]["truth_is_second_half"]):
        target_frac = 1.0  # the released file is the server's working copy itself: split it, do not subsample it
    rows = distribution(T, R, n_seeds, target_frac, ref_frac=frac)
    dist = {}
    for kind in ("floor", "floor_whole_stage", "ceiling"):
        akey = "ceiling" if kind == "ceiling" else "floor"
        dist[kind] = {q: where(pub[m][akey], [r[kind][q] for r in rows]) for q, m in METRIC_OF.items()}
        gaps = np.array([r[kind]["dice_main"] - r[kind]["dice_mirror"] for r in rows])
        dist[kind]["proper_minus_mirror"] = {"mean": round(float(gaps.mean()), 4),
                                            "sd": round(float(gaps.std(ddof=1)), 4) if len(gaps) > 1 else 0.0}
        dist[kind]["prefix_mirror_mode_frac"] = round(float(np.mean([r[kind]["prefix_mode"] == "mirror"
                                                                     for r in rows])), 3)
    out["distribution"] = {"n_seeds": n_seeds, "target_frac": target_frac,
                           "protocol": "vec.evaluate subsample(target, s) / subsample(ref, s+1) / split_half(seed 0), "
                                       "metric seed 0", **dist}

    classes = {}
    for akey in ("floor", "ceiling"):
        pub_d = pub["occupancy_dice"][akey]
        if hit_rows:
            label, basis = classify_exact(pub_d, [h[akey] for h in hit_rows]), "exact"
            prefix_match = all(matches(h[akey]["dice_prefix"], pub_d, 4) for h in hit_rows)
            kind = akey
        else:
            if akey == "floor":
                kind = min(("floor", "floor_whole_stage"),
                           key=lambda k: _absz(dist[k]["d2_shape"]) + _absz(dist[k]["scale_log_ratio"]))
            else:
                kind = "ceiling"
            label, basis, prefix_match = classify_distribution(pub_d, dist[kind]), "distribution", None
        classes[akey] = {"label": label, "basis": basis, "published": pub_d, "compared_with": kind,
                         "published_equals_prefix_value_here": prefix_match}
    out["classification"] = classes

    # reproduction table: published vs recomputed under each build (exact hit if unique, else seed mean +/- sd)
    rep = {}
    for m in SHAPE_METRICS:
        rep[m] = {}
        for akey in ("floor", "ceiling"):
            if len(hit_rows) == 1:
                h = hit_rows[0][akey]
                main = h["dice_main"] if m == "occupancy_dice" else h[m]
                prefix = h["dice_prefix"] if m == "occupancy_dice" else h[m]
                rep[m][akey] = {"published": pub[m][akey], "main": main, "prefix": prefix, "basis": "exact hit"}
            else:
                kind = classes[akey]["compared_with"]
                qm, qp = ("dice_main", "dice_prefix") if m == "occupancy_dice" else (m, m)
                rep[m][akey] = {"published": pub[m][akey],
                                "main": [dist[kind][qm]["mean"], dist[kind][qm]["sd"]],
                                "prefix": [dist[kind][qp]["mean"], dist[kind][qp]["sd"]],
                                "basis": f"seed mean, sd ({kind})"}
    out["reproduction"] = rep

    computed = [r[k] for r in rows for k in ("floor", "floor_whole_stage", "ceiling")]
    computed += [h[k] for h in hit_rows for k in ("floor", "ceiling")]
    both = [h[k] for h in hit_rows for k in ("floor", "ceiling")]
    out["theory"] = {"rows": len(computed),
                     "prefix_equals_parity_rule": sum(r["rule_ok"] for r in computed),
                     "mirror_twin_unchanged_under_prefix": sum(r["twin_ok"] for r in computed),
                     "d2_scale_rows_under_both_builds": len(both),
                     "d2_scale_build_independent": sum(r["d2_scale_build_independent"] for r in both)}
    return out


def compare_indices(before: dict, after: dict, keys) -> dict:
    """Which T2 validation shape anchors changed between two screened repo copies of the index."""
    out = {}
    for k in keys:
        if k not in before or k not in after:
            continue
        ch = {}
        for m in SHAPE_METRICS:
            for akey in ("floor", "ceiling"):
                b = (before[k].get("anchors") or {}).get(m, {}).get(akey)
                a = (after[k].get("anchors") or {}).get(m, {}).get(akey)
                if b != a:
                    ch[f"{m}.{akey}"] = [b, a]
        out[k] = ch
    return out


def summary_lines(report: dict) -> list[str]:
    lines = [f"anchor check {report['generated']} ({report['mode']})"]
    for key, b in report["boards"].items():
        c = b["counts"]
        lines.append(f"## {key}: N_target {b['N_target']}, N_ref {b['N_ref']}, truth_cells {b['index_truth_cells']} "
                     f"(first half {c['truth_is_first_half']}, second half {c['truth_is_second_half']}), "
                     f"ref_cells ok {c['ref_cells_match']}")
        es = b["exact_search"]
        lines.append(f"  exact protocol: {es['n_hits']} hit(s) from {es['n_truth_halves_tried']} truth halves "
                     f"({es['n_ceiling_matches']} ceiling matches), {es['sec']} s")
        for h in b["exact_hits"][:3]:
            lines.append(f"    {h['protocol']}")
            for akey in ("floor", "ceiling"):
                r = h[akey]
                lines.append(f"      {akey}: dice main {r['dice_main']} mirror {r['dice_mirror']} prefix "
                             f"{r['dice_prefix']} (parities {r['parity_pred']:+d}/{r['parity_truth']:+d})")
        for m, per in b["reproduction"].items():
            lines.append(f"  {m}: " + "; ".join(f"{akey} published {v['published']} main {v['main']} prefix "
                                                  f"{v['prefix']} [{v['basis']}]" for akey, v in per.items()))
        for kind in ("floor", "floor_whole_stage", "ceiling"):
            d = b["distribution"][kind]
            lines.append(f"  {kind}: proper-minus-mirror {d['proper_minus_mirror']['mean']} "
                         f"(sd {d['proper_minus_mirror']['sd']}), z(published Dice) proper "
                         f"{d['dice_main'].get('z_of_published')} mirror {d['dice_mirror'].get('z_of_published')}")
        for akey, cl in b["classification"].items():
            lines.append(f"  Dice {akey} anchor {cl['published']}: {cl['label'].upper()} ({cl['basis']}, "
                         f"{cl['compared_with']}; equals the prefix value here: "
                         f"{cl['published_equals_prefix_value_here']})")
        t = b["theory"]
        lines.append(f"  theory: parity rule {t['prefix_equals_parity_rule']}/{t['rows']}, mirror twin unchanged "
                     f"{t['mirror_twin_unchanged_under_prefix']}/{t['rows']}, d2/scale build-independent "
                     f"{t['d2_scale_build_independent']}/{t['d2_scale_rows_under_both_builds']}")
    labels = [cl["label"] for b in report["boards"].values() for cl in b["classification"].values()]
    n_other = len(labels) - sum(labels.count(x) for x in ("proper", "mirror", "tie"))
    lines.append(f"Dice anchors: {len(labels)}; proper {labels.count('proper')}, mirror {labels.count('mirror')}, "
                 f"tie {labels.count('tie')}, other {n_other}")
    if "mirror" in labels:
        lines.append("Reading: at least one Dice anchor is a mirror-mode value, so the bake ran a pre-fix build.")
    elif labels and all(x == "proper" for x in labels):
        lines.append("Reading: every Dice anchor is a proper-mode value (right-handed bake, or pre-fix with lucky "
                     "parities); ask the organisers.")
    if report.get("index_changes"):
        lines.append(f"index changes against --index-before: {report['index_changes']}")
    return lines


def write_report(report: dict, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "anchor_check.json").write_text(json.dumps(report, indent=1, default=float) + "\n")
    lines = summary_lines(report)
    (out / "anchor_check.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


# --------------------------------------------------------------------------------------------- self-test
def bake_entry(T, R, protocol: Protocol, build: str, seeds=DEFAULT_SEARCH_SEEDS, frac: float = SUBSAMPLE_FRAC):
    """A synthetic index entry baked with a known protocol and build, and the true mode of each Dice anchor."""
    Nt, Nr = len(T), len(R)
    nts, nrs = _counts(Nt, frac), _counts(Nr, frac)
    tdraws = target_draws(Nt, nts, Nr, nrs, seeds, None)
    rdraws = ref_draws(Nt, nts, Nr, nrs, seeds)
    A, B, P = realise(protocol, T, R, tdraws, rdraws)
    anchors, modes = {m: {} for m in SHAPE_METRICS}, {}
    for akey, pred in (("floor", P), ("ceiling", B)):
        r = shape_row(pred, A, protocol.metric_seed)
        anchors["d2_shape"][akey] = r["d2_shape"]
        anchors["scale_log_ratio"][akey] = r["scale_log_ratio"]
        anchors["occupancy_dice"][akey] = r["dice_main"] if build == "main" else r["dice_prefix"]
        if r["dice_main"] == r["dice_mirror"]:
            modes[akey] = "tie"
        else:
            modes[akey] = "proper" if build == "main" or r["prefix_mode"] == "proper" else "mirror"
    ref_cells = len(P) if protocol.floor != "whole reference stage" else int(round(frac * Nr))
    return {"truth_cells": len(A), "ref_cells": ref_cells, "split": "val_interp", "anchors": anchors}, modes


SELF_TEST_CASES = (
    # (build, protocol): a non-default implementation baked with the pre-fix build, and the emulator's own protocol
    ("prefix", Protocol("legacy(seed=7, n={nt})", False, "veckit", 7, False, "legacy(seed=8, n={nr}), as drawn", 0)),
    ("main", Protocol("gen(seed=3, n={nt})", True, "emulator", 0, True, "gen(seed=4, n={nr}), sorted", 3)),
)


def _fill(protocol: Protocol, Nt: int, Nr: int, frac: float) -> Protocol:
    nt, nr = int(round(frac * Nt)), int(round(frac * Nr))
    return Protocol(protocol.target_draw.format(nt=nt), protocol.target_sorted, protocol.split, protocol.split_seed,
                    protocol.truth_first, protocol.floor.format(nr=nr), protocol.metric_seed)


def mirror_case(T, R, seeds=DEFAULT_SEARCH_SEEDS, frac: float = SUBSAMPLE_FRAC) -> Protocol | None:
    """An emulator-style protocol whose ceiling the prefix build scores in mirror mode (so a pre-fix bake of it has a
    mirror-mode Dice anchor), or None if no search seed gives one."""
    Nt, Nr = len(T), len(R)
    nts, nrs = _counts(Nt, frac), _counts(Nr, frac)
    tdraws, rdraws = target_draws(Nt, nts, Nr, nrs, seeds, None), ref_draws(Nt, nts, Nr, nrs, seeds)
    for s in seeds:
        proto = _fill(Protocol(f"gen(seed={s}, n={{nt}})", True, "emulator", 0, True,
                               f"gen(seed={s + 1}, n={{nr}}), sorted", 0), Nt, Nr, frac)
        A, B, _ = realise(proto, T, R, tdraws, rdraws)
        h = handedness(B, A, seed=0)
        if h["prefix_mode"] == "mirror" and h["dice_proper"] != h["dice_mirror"]:
            return proto
    return None


def self_test(T, R, seeds=DEFAULT_SEARCH_SEEDS, n_seeds: int = 20, frac: float = SUBSAMPLE_FRAC,
              cases=None) -> dict:
    """Bake synthetic anchors with known protocols and builds, run check_board on them and record whether each
    protocol is recovered and each Dice anchor classified as baked. cases: (build, Protocol) pairs; by default the
    SELF_TEST_CASES plus a pre-fix bake with a mirror-mode ceiling (mirror_case)."""
    if cases is None:
        cases = list(SELF_TEST_CASES)
        auto = mirror_case(T, R, seeds, frac)
        if auto is not None:
            cases.append(("prefix", auto))
    results = {}
    for i, (build, proto) in enumerate(cases):
        proto = _fill(proto, len(T), len(R), frac)
        entry, modes = bake_entry(T, R, proto, build, seeds, frac)
        key = f"{SYNTHETIC_PREFIX}{i}:{build}"
        res = check_board(key, T, R, entry, seeds, n_seeds, frac)
        recovered = [q for h in res["exact_hits"] for q in h["equivalent"]]
        labels = {akey: res["classification"][akey]["label"] for akey in ("floor", "ceiling")}
        ok_protocol = asdict(proto) in recovered
        ok_modes = all(labels[a] == modes[a] for a in modes)
        res["self_test"] = {"build": build, "baked_protocol": asdict(proto), "true_modes": modes,
                            "classified": labels, "protocol_recovered": ok_protocol,
                            "n_hits": res["exact_search"]["n_hits"],
                            "modes_recovered": ok_modes, "passed": ok_protocol and ok_modes}
        results[key] = res
    return results


# --------------------------------------------------------------------------------------------- command line
def resolve_boards(board_args, data_dir: Path) -> list[tuple[str, Path | None, Path]]:
    if board_args:
        return [(k, Path(t) if Path(t).exists() else None, Path(r)) for k, t, r in board_args]
    out = []
    for key, (stage, ref_stage) in T2_BOARDS.items():
        out.append((key, find_answer(data_dir, stage), data_dir / "T2" / f"{ref_stage}.h5ad"))
    return out


def _seeds(text: str) -> tuple[int, ...]:
    return tuple(int(x) for x in text.split(",") if x.strip())


RELEASE_UTC = "2026-10-20"  # validation answers become training data; no anchor is read before this date


def cmd_check(a) -> int:
    data = Path(a.data)
    boards = resolve_boards(a.board, data)
    synthetic = [k for k, _, _ in boards if str(k).startswith(SYNTHETIC_PREFIX)]
    if synthetic:
        raise SystemExit(f"synthetic boards {synthetic} are for self-test only: refused in check mode.")
    check_board_keys([k for k, _, _ in boards])
    present = [(k, t, r) for k, t, r in boards if t is not None]
    if not present:
        print(f"No T2 validation answers found (looked under {data} for "
              f"{', '.join(s for s, _ in T2_BOARDS.values())}; released 2026-10-20). Nothing was read. Exiting.")
        return 0
    if time.strftime("%Y-%m-%d", time.gmtime()) < RELEASE_UTC:
        print(f"Before the {RELEASE_UTC} release the validation anchors stay unread, even with stand-in files "
              "(use self-test). Exiting.")
        return 0
    missing = [k for k, t, _ in boards if t is None]
    entries = load_index_entries(a.index, [k for k, _, _ in present])
    report = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "mode": "check", "index": Path(a.index).name,
              "boards": {}, "missing_answers": missing,
              "thresholds": {"z_inside": Z_INSIDE, "z_outside": Z_OUTSIDE, "z_fingerprint": Z_FINGERPRINT}}
    if a.index_before:
        before = load_index_entries(a.index_before, [k for k, _, _ in present])
        report["index_changes"] = compare_indices(before, entries, [k for k, _, _ in present])
    for key, tpath, rpath in present:
        T, R = load_coords(tpath), load_coords(rpath)
        report["boards"][key] = check_board(key, T, R, entries[key], _seeds(a.search_seeds), a.seeds)
        report["boards"][key]["files"] = {"target": Path(tpath).name, "reference": Path(rpath).name}
        print(f"{key} done", flush=True)
        write_report(report, Path(a.out))
    return 0


def cmd_self_test(a) -> int:
    T, R = load_coords(a.target), load_coords(a.ref)
    results = self_test(T, R, _seeds(a.search_seeds), a.seeds)
    report = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"), "mode": "self-test",
              "files": {"target": Path(a.target).name, "reference": Path(a.ref).name}, "boards": results}
    write_report(report, Path(a.out))
    verdicts = {k: v["self_test"] for k, v in results.items()}
    print(json.dumps(verdicts, indent=1))
    return 0 if all(v["passed"] for v in verdicts.values()) else 1


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in ("check", "self-test", "-h", "--help"):
        argv = ["check", *argv]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="reproduce and classify the published T2 validation shape anchors")
    c.add_argument("--index", default=str(REPO_ROOT / "boards" / "index.json"),
                   help="the repo's boards/index.json, refreshed only via vec.boards.update_boards")
    c.add_argument("--index-before", help="an older screened repo copy of boards/index.json (e.g. from git)")
    c.add_argument("--data", default=str(REPO_ROOT / "data" / "raw"))
    c.add_argument("--board", nargs=3, action="append", metavar=("KEY", "TARGET_H5AD", "REF_H5AD"))
    c.add_argument("--seeds", type=int, default=50, help="seeds of the emulator protocol for the distributions")
    c.add_argument("--search-seeds", default=",".join(map(str, DEFAULT_SEARCH_SEEDS)))
    c.add_argument("--out", default=str(REPO_ROOT / "runs" / "scorer_build" / "anchor_check"))
    s = sub.add_parser("self-test", help="bake synthetic anchors from released stages with known protocols and "
                                         "check that they are recovered and classified")
    s.add_argument("--target", required=True)
    s.add_argument("--ref", required=True)
    s.add_argument("--seeds", type=int, default=20)
    s.add_argument("--search-seeds", default=",".join(map(str, DEFAULT_SEARCH_SEEDS)))
    s.add_argument("--out", default=str(REPO_ROOT / "runs" / "scorer_build" / "anchor_selftest"))
    a = ap.parse_args(argv)
    return cmd_check(a) if a.cmd == "check" else cmd_self_test(a)


if __name__ == "__main__":
    raise SystemExit(main())
