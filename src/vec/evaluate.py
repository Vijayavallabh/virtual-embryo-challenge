"""Local board emulator: score a prediction the way the challenge server does, as far as it is published.

What is reproduced exactly: the metric code (the pinned veckit functions are called directly), the skill map, and the
group weights of proposal Table 1. What is emulated: the server subsamples every stage to 10% and splits the target in
two halves at a fixed seed (one half is the truth, the other half is the ceiling prediction); the floor is the
copy_last / wt_identity prediction scored against the same truth. The server's subsampling seeds are not published,
so scores are reported per seed and averaged. Runbook step 4.6 (2026-10-28) recovered the validation bake exactly from
the T2 shape anchors; `score_board(..., protocol="bake")` reproduces it (at seed 0, the server's own draw).

Known gap: Task 2's shape group is "gated" by neighbourhood-structure skill on the live board, but the exact gate is not
published. `gate` selects an assumption; calibrate it against the public reference rows once the validation answers
are out.

Scorer builds (`build`). veckit's `occupancy_dice` (ranked on T2) and `sliced_wasserstein` (unranked) canonicalise
each cloud with `common.shape_metrics._canonicalise` (centre, PCA basis from `np.linalg.svd`, RMS scale) and then
maximise over the four proper axis-sign flips. The two builds differ only in that function:

- "main" (default): the pinned git cc3859a (PR #9, issue #7). The PCA basis is made right-handed, so a proper rotation
  of a prediction scores the same and a mirror image is compared as a mirror image.
- "prefix": veckit 0.1.0 / 0.1.1 (git 46d41e6; PyPI 0.1.2 also lacks the flip, per the screened research). The basis
  is used as LAPACK returns it (det +1 or -1). When the prediction's and the truth's bases have the same parity, the
  comparison is the proper one ("proper mode", the "main" value); otherwise the prediction is compared with the
  truth's mirror image ("mirror mode"). The parity follows the SVD sign pattern, which follows the rows and their
  order, so it acts as a coin per file and truth half, and a z-mirrored copy of a file scores exactly like the file.

The screened research of 2026-10-18 found the live validation board behaving like "prefix"; the default stays "main"
until the lead decides. `d2_shape` and `scale_log_ratio` never call `_canonicalise`, so they are identical under both
builds. Whatever the build, `score_board` reports for a T2 candidate its proper-mode Dice, its mirror-mode Dice, which
of the two the prefix build returns on the emulated truth, and the file's local parity under the prefix build
(`extras["handedness"]`, see `handedness`).
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np

from vec.boards import ABS_METRICS, GROUPS, LOWER_IS_BETTER, RANKED_METRICS, Board
from vec.io import COORD_KEY, dense

SUBSAMPLE_FRAC = 0.1
SCORER_BUILDS = ("main", "prefix")
DEFAULT_BUILD = "main"
PROTOCOLS = ("emulator", "bake")  # how score_board draws the working copies and the truth half


def _veckit_modules(task: str):
    """(metrics, metrics_v2) for a task, loaded through veckit's own loader so the code path is the official one."""
    import score_h5ad

    return score_h5ad._load_task_metrics(task)


def _shape_metrics():
    import common.shape_metrics as sm

    return sm


def _skill(value, floor, ceiling, lower_is_better):
    from common.core_metrics import skill

    return skill(value, floor, ceiling, lower_is_better=lower_is_better)


# ---------------------------------------------------------------------------------------------------------------------
# Pre-fix canonicalisation: veckit 0.1.1 (git 46d41e6, https://github.com/aristoteleo/veckit),
# common/shape_metrics.py::_canonicalise, copied and reformatted to this repo's lint rules with every array operation
# unchanged (so LAPACK makes the same sign choices). The pinned git cc3859a differs from it only by the two-line
# right-handed flip. Licence of the copied code:
#
#   MIT License
#
#   Copyright (c) 2026 Virtual Embryo Challenge organizers
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
#   copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
# ---------------------------------------------------------------------------------------------------------------------
def _canonicalise_prefix(C):
    """Centre, rotate onto principal axes, scale by RMS radius. Returns (Z, rms_radius). (veckit 0.1.1, pre-fix.)

    PCA leaves each axis' sign undetermined; callers minimise over the sign combinations with determinant +1. This
    build does not make the basis right-handed, so its handedness is whatever LAPACK returns.
    """
    C = np.asarray(C, float)
    mu = C.mean(0)
    X = C - mu
    rms = float(np.sqrt((X**2).sum(1).mean()))
    if rms < 1e-12:
        return X, rms
    _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    return (X @ Vt.T) / rms, rms


_PINNED_CANONICALISE = None


def _pinned_canonicalise():
    """veckit's own `_canonicalise` as installed (the pinned git cc3859a), captured before any swap."""
    global _PINNED_CANONICALISE
    if _PINNED_CANONICALISE is None:
        current = _shape_metrics()._canonicalise
        if current is _canonicalise_prefix:
            raise RuntimeError("common.shape_metrics._canonicalise was swapped before the pinned one was recorded")
        _PINNED_CANONICALISE = current
    return _PINNED_CANONICALISE


@contextmanager
def scorer_build(build: str = DEFAULT_BUILD):
    """Run veckit's metrics under one scorer build: "main" runs the pinned code untouched, "prefix" swaps in the
    pre-fix `_canonicalise` for the duration. The previous function is always restored (nesting works). The swap is
    a module-global monkeypatch, so do not score from several threads at once."""
    if build not in SCORER_BUILDS:
        raise ValueError(f"build must be one of {SCORER_BUILDS}, got {build!r}")
    sm = _shape_metrics()
    pinned = _pinned_canonicalise()
    previous = sm._canonicalise
    sm._canonicalise = pinned if build == "main" else _canonicalise_prefix
    try:
        yield
    finally:
        sm._canonicalise = previous


def mirror_z(coords) -> np.ndarray:
    """The z-mirrored copy of a cloud (an improper transform; same rows, same order)."""
    C = np.array(coords, copy=True)
    C[:, 2] = -C[:, 2]
    return C


def svd_parity(coords) -> int:
    """Parity (+1 or -1) of the PCA basis the pre-fix `_canonicalise` builds for these rows, computed with the same
    array operations so LAPACK makes the same sign choices. 0 for a degenerate cloud (fewer than 3 points or no
    spread), where the pre-fix function does not rotate."""
    C = np.asarray(coords, float)
    X = C - C.mean(0)
    if X.shape[0] < 3 or float(np.sqrt((X**2).sum(1).mean())) < 1e-12:
        return 0
    _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    return int(np.sign(np.linalg.det(Vt)))


def match_rows(pred_coords, true_coords, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """The rows veckit's `_match_n` compares: both clouds subsampled (or permuted) to a common count by
    default_rng(seed), the prediction drawn first."""
    A, B, _ = _shape_metrics()._match_n(np.asarray(pred_coords, float), np.asarray(true_coords, float), seed)
    return A, B


def prefix_parity(pred_coords, n_true: int, seed: int = 0) -> int:
    """A file's local parity under the prefix build: the parity of the rows `_match_n` takes from the prediction when
    the truth has `n_true` cells, at metric seed `seed`. It depends on the prediction's rows, their order and count,
    n_true and the seed, never on the truth's coordinates. Pass coordinates as the scorer reads them (float32).

    For the live board: seed 0 and n_true = the board's truth_cells, assuming the server's metric seed is fixed at 0
    and its LAPACK picks SVD signs like numpy/OpenBLAS (screened research of 2026-10-18; neither is shown).
    """
    A, _ = match_rows(pred_coords, np.zeros((int(n_true), 3)), seed)
    return svd_parity(A)


def handedness(pred_coords, true_coords, seed: int = 0, live_n_true: int | None = None, live_seed: int = 0) -> dict:
    """Where a prediction stands under the two builds against one truth (T2 occupancy_dice, rounded as the panel does).

    dice_proper: the "main" Dice, i.e. the comparison with the truth as it is. dice_mirror: the prediction against
    the truth's mirror image (the "main" Dice of the z-mirrored prediction). The prefix build returns dice_proper when
    the two clouds' SVD bases (the rows `_match_n` hands over at metric seed `seed`) have the same parity and
    dice_mirror otherwise: prefix_mode and dice_prefix. local_parity is the file's own parity as the live board would
    compute it (`prefix_parity` at `live_seed` with `live_n_true` truth cells, default this truth's count).
    """
    sm = _shape_metrics()
    with scorer_build("main"):
        proper = sm.occupancy_dice(pred_coords, true_coords, seed=seed)[0]
        mirror = sm.occupancy_dice(mirror_z(pred_coords), true_coords, seed=seed)[0]
    A, B = match_rows(pred_coords, true_coords, seed)
    parity_pred, parity_truth = svd_parity(A), svd_parity(B)
    mode = "proper" if parity_pred == parity_truth else "mirror"
    n_live = int(live_n_true) if live_n_true is not None else int(np.asarray(true_coords).shape[0])
    return {
        "dice_proper": round(proper, 4),
        "dice_mirror": round(mirror, 4),
        "dice_gap": round(proper - mirror, 4),
        "parity_pred": parity_pred,
        "parity_truth": parity_truth,
        "prefix_mode": mode,
        "dice_prefix": round(proper if mode == "proper" else mirror, 4),
        "metric_seed": int(seed),
        "local_parity": prefix_parity(pred_coords, n_live, live_seed),
        "local_parity_n_true": n_live,
        "local_parity_seed": int(live_seed),
    }


class _NoProbe:
    """Stand-in for the frozen cell-type probe when composition_JSD is not wanted (it is not a ranked metric)."""

    def __init__(self, classes):
        self.classes_ = np.asarray(classes)

    def predict_proba(self, X):
        return np.full((X.shape[0], len(self.classes_)), 1.0 / len(self.classes_))


@dataclass
class Stage:
    """Expression (cells x board genes, log-normalised), optional coordinates and cell-type labels."""

    X: np.ndarray
    coords: np.ndarray | None = None
    celltype: np.ndarray | None = None

    @classmethod
    def from_anndata(cls, a, genes, need_coords=False):
        a = a[:, genes]
        coords = np.asarray(a.obsm[COORD_KEY], dtype=np.float32)[:, :3] if need_coords else None
        ct = np.asarray(a.obs["celltype"]).astype(str) if "celltype" in a.obs else None
        return cls(X=dense(a.X).astype(np.float32), coords=coords, celltype=ct)

    @property
    def n(self) -> int:
        return self.X.shape[0]

    def take(self, idx) -> Stage:
        return Stage(
            X=self.X[idx],
            coords=None if self.coords is None else self.coords[idx],
            celltype=None if self.celltype is None else self.celltype[idx],
        )


def subsample(stage: Stage, frac: float, seed: int) -> Stage:
    rng = np.random.default_rng(seed)
    n = max(2, int(round(stage.n * frac)))
    return stage.take(np.sort(rng.choice(stage.n, n, replace=False)))


def split_half(stage: Stage, seed: int = 0) -> tuple[Stage, Stage]:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(stage.n)
    return stage.take(np.sort(idx[: stage.n // 2])), stage.take(np.sort(idx[stage.n // 2:]))


def veckit_halves(stage: Stage, seed: int = 0) -> tuple[Stage, Stage]:
    """veckit's own split (`common.core_metrics.split_half`), rows kept in veckit's order (not sorted)."""
    from common.core_metrics import split_half as veckit_split_half

    a, b = veckit_split_half(stage.n, seed=seed)
    return stage.take(a), stage.take(b)


def panel_metrics(task: str, pred: Stage, truth: Stage, ref: Stage, seed: int = 0, probe="skip",
                  build: str = DEFAULT_BUILD) -> dict:
    """The full official metric panel for one prediction against one truth, under scorer build `build`."""
    metrics, metrics_v2 = _veckit_modules(task)
    true_ct = truth.celltype if truth.celltype is not None else np.array(["NA"] * truth.n)
    if probe == "skip":
        probe_obj = _NoProbe(np.unique(true_ct))
    elif probe == "train":
        probe_obj = metrics.train_frozen_probe(truth.X, true_ct)
    else:
        probe_obj = probe
    pred_ct = np.array(["NA"] * pred.n)
    with scorer_build(build):
        if task == "T1":
            out = metrics_v2.score_task1_v2(pred.X, pred_ct, truth.X, true_ct, ref.X, probe=probe_obj, seed=seed)
        elif task == "T2":
            out = metrics_v2.score_task2_v2(pred.X, pred.coords, pred_ct, truth.X, truth.coords, true_ct, ref.X,
                                            probe=probe_obj, seed=seed)
        elif task == "T3":
            out = metrics_v2.score_task3_v2(pred.X, pred.coords, pred_ct, truth.X, truth.coords, true_ct, ref.X,
                                            probe=probe_obj, seed=seed)
        else:
            raise ValueError(f"unknown task {task!r}")
    if probe == "skip":
        out["composition_JSD"] = None
    return out


def _metric_value(m: str, values: dict):
    v = values.get(m)
    if v is None or not np.isfinite(v):
        return float("nan")
    return abs(v) if m in ABS_METRICS else float(v)


def skills_from_anchors(task: str, values: dict, anchors: dict) -> dict:
    out = {}
    for m in RANKED_METRICS[task]:
        a = anchors.get(m)
        if a is None:
            out[m] = float("nan")
            continue
        floor, ceiling = a["floor"], a["ceiling"]
        if m in ABS_METRICS:
            floor, ceiling = abs(floor), abs(ceiling)
        out[m] = _skill(_metric_value(m, values), floor, ceiling, m in LOWER_IS_BETTER)
    return out


def inverted_anchors(anchors: dict) -> list[str]:
    """Metrics whose ceiling is no better than their floor. The skill map then gives the floor 1.0, so scores on
    such a board are not meaningful (it happens on tiny samples; the published boards have none)."""
    bad = []
    for m, a in anchors.items():
        f, c = (abs(a["floor"]), abs(a["ceiling"])) if m in ABS_METRICS else (a["floor"], a["ceiling"])
        if (m in LOWER_IS_BETTER and c >= f) or (m not in LOWER_IS_BETTER and c <= f):
            bad.append(m)
    return bad


def aggregate(task: str, skills: dict, gate: str = "floor_normalised") -> tuple[float, dict]:
    """Task score on the 0-100 board scale and per-group skills. NaN skills count as 0, as on the server.

    gate (Task 2 shape group only): "none" = ungated; "multiplicative" = shape * neighbourhood skill (the floor would
    then sit below 50); "floor_normalised" = shape * min(1, neighbourhood skill / 0.5), which keeps the floor at 50.
    """
    group_skill = {}
    for name, _, members in GROUPS[task]:
        s = sum(w * (0.0 if not np.isfinite(skills.get(m, np.nan)) else skills[m]) for m, w in members)
        group_skill[name] = s
    if task == "T2" and gate != "none":
        nb = skills.get("neighborhood_mmd", np.nan)
        nb = 0.0 if not np.isfinite(nb) else nb
        factor = nb if gate == "multiplicative" else min(1.0, nb / 0.5)
        group_skill["shape_growth"] *= factor
    score = 100.0 * sum(w * group_skill[name] for name, w, _ in GROUPS[task])
    return score, group_skill


@dataclass
class BoardScore:
    score: float
    group_skill: dict
    skills: dict
    metrics: dict
    anchors: dict
    seed: int
    extras: dict = field(default_factory=dict)


def score_board(board: Board, pred: Stage, target: Stage, ref: Stage, floor_pred: Stage | None = None,
                seed: int = 0, anchors: str = "computed", gate: str = "floor_normalised",
                probe="skip", frac: float = SUBSAMPLE_FRAC, build: str = DEFAULT_BUILD,
                protocol: str = "emulator") -> BoardScore:
    """Score `pred` on a board, emulating the server's 10% working copies and half-split truth.

    target: the full held-out stage (e.g. released validation answers). ref: the full DE reference stage (preceding
    stage for T1/T2, matched wild type for T3). floor_pred: the full floor prediction; defaults to `ref`, which is
    what copy_last and wt_identity return. anchors: "computed" measures floor and ceiling on these working copies;
    "published" uses the board's anchors from index.json (validation boards only). frac: working-copy fraction;
    0.1 is what the server does, 1.0 is useful for small files and low-noise diagnostics. build: the scorer build
    ("main" or "prefix", see the module docstring) for the candidate and for computed anchors. protocol: how the
    working copies and the truth half are drawn. "emulator" draws the target with `seed` and the reference with
    `seed + 1` and takes sorted halves of a seeded permutation. "bake" is the validation bake that runbook step 4.6
    recovered exactly on all three T2 validation boards (2026-10-28, experiments/scorer_build/results/release/): the
    target and the reference drawn with the same `seed`, the truth the first half of veckit's split_half(n, 0) in
    veckit's row order; at seed 0 these are the server's own working copies and truth half.

    On T2 the candidate's handedness diagnostics go to extras["handedness"] whatever the build (dice_proper,
    dice_mirror, prefix_mode on this truth half, and local_parity at metric seed 0 with the board's truth_cells);
    with computed anchors, the floor's and ceiling's go to extras["anchor_handedness"].
    """
    if build not in SCORER_BUILDS:
        raise ValueError(f"build must be one of {SCORER_BUILDS}, got {build!r}")
    if protocol not in PROTOCOLS:
        raise ValueError(f"protocol must be one of {PROTOCOLS}, got {protocol!r}")
    task = board.task
    ref_seed = seed + 1 if protocol == "emulator" else seed
    target_w = subsample(target, frac, seed)
    ref_w = subsample(ref, frac, ref_seed)
    truth, other_half = split_half(target_w, seed=0) if protocol == "emulator" else veckit_halves(target_w)
    values = panel_metrics(task, pred, truth, ref_w, seed=seed, probe=probe, build=build)

    extras = {"build": build, "protocol": protocol}
    floor_w = None
    if anchors == "computed":
        floor_w = ref_w if floor_pred is None else subsample(floor_pred, frac, ref_seed)
        floor_vals = panel_metrics(task, floor_w, truth, ref_w, seed=seed, build=build)
        ceil_vals = panel_metrics(task, other_half, truth, ref_w, seed=seed, build=build)
        anchor_map = {m: {"floor": floor_vals.get(m), "ceiling": ceil_vals.get(m)} for m in RANKED_METRICS[task]}
        anchor_map = {m: a for m, a in anchor_map.items()
                      if a["floor"] is not None and a["ceiling"] is not None}
        extras["floor_metrics"], extras["ceiling_metrics"] = floor_vals, ceil_vals
    elif anchors == "published":
        if not board.anchors:
            raise ValueError(f"{board.key}: no published anchors (test boards are stripped at intake)")
        anchor_map = board.anchors
    else:
        raise ValueError("anchors must be 'computed' or 'published'")

    if task == "T2" and pred.coords is not None and truth.coords is not None:
        live_n = board.truth_cells if board.truth_cells else truth.n
        hand = handedness(pred.coords, truth.coords, seed=seed, live_n_true=live_n)
        # the panel's Dice must be the value this build returns: proper mode under main, the parity rule under prefix
        hand["panel_dice_matches_build"] = values.get("occupancy_dice") == hand[
            "dice_prefix" if build == "prefix" else "dice_proper"]
        extras["handedness"] = hand
        if floor_w is not None:
            extras["anchor_handedness"] = {
                kind: handedness(P.coords, truth.coords, seed=seed, live_n_true=live_n)
                for kind, P in (("floor", floor_w), ("ceiling", other_half)) if P.coords is not None}

    skills = skills_from_anchors(task, values, anchor_map)
    score, group_skill = aggregate(task, skills, gate=gate)
    extras.update(truth_cells=truth.n, ref_cells=ref_w.n, pred_cells=pred.n,
                  inverted_anchors=inverted_anchors(anchor_map))
    return BoardScore(score=score, group_skill=group_skill, skills=skills, metrics=values,
                      anchors=anchor_map, seed=seed, extras=extras)


def _handedness_summary(hands: list[dict]) -> dict:
    gap = np.array([h["dice_gap"] for h in hands], float)
    return {
        "dice_proper_mean": float(np.mean([h["dice_proper"] for h in hands])),
        "dice_mirror_mean": float(np.mean([h["dice_mirror"] for h in hands])),
        "dice_gap_mean": float(gap.mean()),
        "dice_gap_sd": float(gap.std(ddof=1)) if len(gap) > 1 else 0.0,
        "prefix_mirror_mode_frac": float(np.mean([h["prefix_mode"] == "mirror" for h in hands])),
        "local_parity": sorted({h["local_parity"] for h in hands}),
        "local_parity_n_true": hands[0]["local_parity_n_true"],
    }


def score_board_seeds(board: Board, pred: Stage, target: Stage, ref: Stage, seeds=(0, 1, 2), **kw) -> dict:
    """Score over several subsampling seeds; returns per-seed results plus mean and sd of the board score (and, on T2,
    a summary of the candidate's handedness diagnostics over the seeds)."""
    runs = [score_board(board, pred, target, ref, seed=s, **kw) for s in seeds]
    scores = np.array([r.score for r in runs])
    out = {
        "board": board.key,
        "score_mean": float(scores.mean()),
        "score_sd": float(scores.std(ddof=1)) if len(scores) > 1 else 0.0,
        "per_seed": [
            {"seed": r.seed, "score": r.score, "group_skill": r.group_skill, "skills": r.skills,
             "metrics": r.metrics, "anchors": r.anchors, "extras": {k: v for k, v in r.extras.items()
                                                                     if not k.endswith("_metrics")}}
            for r in runs
        ],
        "build": runs[0].extras.get("build", DEFAULT_BUILD),
    }
    hands = [r.extras.get("handedness") for r in runs]
    if hands and all(h is not None for h in hands):
        out["handedness"] = _handedness_summary(hands)
    return out
