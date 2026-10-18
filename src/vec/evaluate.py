"""Local board emulator: score a prediction the way the challenge server does, as far as it is published.

What is reproduced exactly: the metric code (the pinned veckit functions are called directly), the skill map, and the
group weights of proposal Table 1. What is emulated: the server subsamples every stage to 10% and splits the target in
two halves at a fixed seed (one half is the truth, the other half is the ceiling prediction); the floor is the
copy_last / wt_identity prediction scored against the same truth. The server's subsampling seeds are not published,
so scores are reported per seed and averaged.

Known gap: Task 2's shape group is "gated" by neighbourhood-structure skill on the live board, but the exact gate is not
published. `gate` selects an assumption; calibrate it against the public reference rows once the validation answers
are out.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from vec.boards import ABS_METRICS, GROUPS, LOWER_IS_BETTER, RANKED_METRICS, Board
from vec.io import COORD_KEY, dense

SUBSAMPLE_FRAC = 0.1


def _veckit_modules(task: str):
    """(metrics, metrics_v2) for a task, loaded through veckit's own loader so the code path is the official one."""
    import score_h5ad

    return score_h5ad._load_task_metrics(task)


def _skill(value, floor, ceiling, lower_is_better):
    from common.core_metrics import skill

    return skill(value, floor, ceiling, lower_is_better=lower_is_better)


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


def panel_metrics(task: str, pred: Stage, truth: Stage, ref: Stage, seed: int = 0, probe="skip") -> dict:
    """The full official metric panel for one prediction against one truth."""
    metrics, metrics_v2 = _veckit_modules(task)
    true_ct = truth.celltype if truth.celltype is not None else np.array(["NA"] * truth.n)
    if probe == "skip":
        probe_obj = _NoProbe(np.unique(true_ct))
    elif probe == "train":
        probe_obj = metrics.train_frozen_probe(truth.X, true_ct)
    else:
        probe_obj = probe
    pred_ct = np.array(["NA"] * pred.n)
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
                probe="skip", frac: float = SUBSAMPLE_FRAC) -> BoardScore:
    """Score `pred` on a board, emulating the server's 10% working copies and half-split truth.

    target: the full held-out stage (e.g. released validation answers). ref: the full DE reference stage (preceding
    stage for T1/T2, matched wild type for T3). floor_pred: the full floor prediction; defaults to `ref`, which is
    what copy_last and wt_identity return. anchors: "computed" measures floor and ceiling on these working copies;
    "published" uses the board's anchors from index.json (validation boards only). frac: working-copy fraction;
    0.1 is what the server does, 1.0 is useful for small files and low-noise diagnostics.
    """
    task = board.task
    target_w = subsample(target, frac, seed)
    ref_w = subsample(ref, frac, seed + 1)
    truth, other_half = split_half(target_w, seed=0)
    values = panel_metrics(task, pred, truth, ref_w, seed=seed, probe=probe)

    extras = {}
    if anchors == "computed":
        floor_w = ref_w if floor_pred is None else subsample(floor_pred, frac, seed + 1)
        floor_vals = panel_metrics(task, floor_w, truth, ref_w, seed=seed)
        ceil_vals = panel_metrics(task, other_half, truth, ref_w, seed=seed)
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

    skills = skills_from_anchors(task, values, anchor_map)
    score, group_skill = aggregate(task, skills, gate=gate)
    extras.update(truth_cells=truth.n, ref_cells=ref_w.n, pred_cells=pred.n,
                  inverted_anchors=inverted_anchors(anchor_map))
    return BoardScore(score=score, group_skill=group_skill, skills=skills, metrics=values,
                      anchors=anchor_map, seed=seed, extras=extras)


def score_board_seeds(board: Board, pred: Stage, target: Stage, ref: Stage, seeds=(0, 1, 2), **kw) -> dict:
    """Score over several subsampling seeds; returns per-seed results plus mean and sd of the board score."""
    runs = [score_board(board, pred, target, ref, seed=s, **kw) for s in seeds]
    scores = np.array([r.score for r in runs])
    return {
        "board": board.key,
        "score_mean": float(scores.mean()),
        "score_sd": float(scores.std(ddof=1)) if len(scores) > 1 else 0.0,
        "per_seed": [
            {"seed": r.seed, "score": r.score, "group_skill": r.group_skill, "skills": r.skills,
             "metrics": r.metrics, "anchors": r.anchors, "extras": {k: v for k, v in r.extras.items()
                                                                     if not k.endswith("_metrics")}}
            for r in runs
        ],
    }
