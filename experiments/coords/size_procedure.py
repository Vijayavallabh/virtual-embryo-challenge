"""Pre-registered size procedure for the two T2 test boards: heart E12.5 (extrapolation) and embryo E7.75
(interpolation).

Rules this file obeys:
- Every number is computed from the released training files passed in. There are no size numbers in this file; the
  only constants are the forms of the rules (e.g. "half" the latest rate) and the selection settings in PREREGISTERED.
- It never reads boards/index.json anchors, leaderboard scores, the PyPI veckit 0.1.2 docstrings or veckit PR #6, and
  it takes no external data. Nominal stage times (E-days) are inputs, taken from the released file names.
- RMS radius is the scorer's definition (veckit common.shape_metrics.scale_log_ratio): sqrt(mean squared distance to
  the centroid) over all three raw coordinates, in each file's own convention.

Why z gets special treatment (experiments/coords/REPORT.md, parts 1-2): in every released file the cells fill the
section stack uniformly along z (SD of z = 0.91-1.02 x range/sqrt(12)) and reach the end sections, and every stack has
30-40 sections whatever the stage. So the z extent is set by how many sections were imaged and the stored step, not
by tissue growth, and in-plane (xy) growth is the anatomical signal when the sectioning plane is comparable.

Heart E12.5 (`heart_size_procedure`):
  1. For each heart training stage: geometry of the core-cardiac cells (obs["core"] from vec.prep.load_stage:
     cardiomyocytes and endocardium) and of the full cloud: raw 3D RMS, in-plane RMS (x, y), SD of z, and how close the
     core's thinnest principal axis is to z.
  2. z mode (PREREGISTERED["heart_z_mode"] = "auto"): "xy_grow_z_hold" when the core's thinnest principal axis is
     within 30 degrees of z in both of the last two stages (comparable sectioning plane), otherwise "raw".
       raw:            grow the core's 3D RMS; scope ratio on 3D RMS.
       xy_grow_z_hold: grow the core's in-plane RMS; scope ratio on in-plane RMS; hold the core and full-cloud z SD at
                       the reference stage's values; recombine in quadrature.
       common_step:    as raw after rescaling every stage's z to the reference stage's step (one-pitch hypothesis).
       split:          grow the in-plane RMS and the z SD with separately selected rules.
  3. Candidate growth rules (RULES) on the grown series against t. A rule may only use rate intervals that start at
     or after heart_min_rate_start_t (E9.0): before that the core is being remodelled (crescent -> tube -> looping), its
     RMS is flat and then jumps, so rates measured there say nothing about later growth (REPORT.md, part 4). Rules
     that would need an earlier interval are not applicable (decaying_rate needs its last two intervals in range; the
     others need the last one). If no rule is applicable the result is flagged valid=False. Applicable rules are
     scored on leave-future-out pseudo-targets from the same stages (predict stage j from the stages before it, and
     from the stages before j-1) whose last training interval also starts at or after E9.0; rules within tie_tol of
     the lowest mean |log error| are averaged in log space; with no eligible pseudo-target (the case with the stages
     released by 2026-10-20) every applicable rule is averaged.
  4. Scope s = log(full / core) at the target under three pre-registered hypotheses:
       H1 "trend": s decays exponentially at the rate fitted on the last two stages (held if not decreasing);
       H0 "hold":  the reference (latest) stage's s;
       H2 "half":  half the reference stage's s (heart-focused scope, shell influence halved).
  5. Output: predicted core and full-cloud RMS per scope hypothesis. A submission rescales its coordinates about the
     centroid to the full-cloud RMS of the scope hypothesis it implements (`rescale_to_rms`; proper scaling only).
Embryo E7.75 (`embryo_size_procedure`):
  PCHIP (monotone cubic) interpolation of log full-cloud RMS against t through every embryo training stage (z mode
  "raw"), with leave-one-out errors on the interior stages as the check.

Refit on 2026-10-20: call the same functions with the newly released stages added (heart: E8.5, E10.5; embryo:
E7.5). The forms, candidate rules, selection settings, z-mode switch and scope hypotheses stay exactly as here.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

PREREGISTERED = {
    "heart_rules": ["latest_rate", "half_rate", "power_t", "decaying_rate"],
    "heart_min_rate_start_t": 9.0,
    "tie_tol": 0.15,
    "heart_scope_hypotheses": ["H1", "H0", "H2"],
    "heart_z_mode": "auto",
    "auto_thin_axis_max_angle_deg": 30.0,
    "embryo_interpolant": "pchip",
    "embryo_z_mode": "raw",
    "core_key": "core",
}


# ----------------------------------------------------------------------------------------------- geometry
def rms_radius(C: np.ndarray) -> float:
    X = np.asarray(C, float) - np.asarray(C, float).mean(0)
    return float(np.sqrt((X ** 2).sum(1).mean()))


def z_levels(z: np.ndarray, tol: float = 1e-3) -> np.ndarray:
    """Distinct section levels (z values closer than tol are one level: floating-point duplicates exist)."""
    u = np.unique(z)
    return u[np.concatenate([[True], np.diff(u) > tol])]


def thin_axis_cos_z(C: np.ndarray) -> float:
    """|cos| between the cloud's smallest principal axis and the z axis."""
    X = C - C.mean(0)
    w, V = np.linalg.eigh(X.T @ X / X.shape[0])
    return float(abs(V[2, int(np.argmin(w))]))


@dataclass
class StageGeometry:
    name: str
    t: float
    n_cells: int
    n_core: int
    n_levels: int
    z_step: float
    z_range: float
    rms_full: float
    rms_core: float
    rms_xy_full: float
    sd_z_full: float
    rms_xy_core: float
    sd_z_core: float
    core_thin_axis_cos_z: float

    @property
    def log_scope(self) -> float:
        return float(np.log(self.rms_full / self.rms_core))

    @property
    def log_scope_xy(self) -> float:
        return float(np.log(self.rms_xy_full / self.rms_xy_core))


def geometry_from_arrays(name: str, t: float, C: np.ndarray, core: np.ndarray) -> StageGeometry:
    C = np.asarray(C, float)[:, :3]
    core = np.asarray(core, bool)
    lv = z_levels(C[:, 2])
    ok = core.sum() >= 4
    nan = float("nan")
    return StageGeometry(
        name=name, t=float(t), n_cells=int(C.shape[0]), n_core=int(core.sum()), n_levels=int(lv.size),
        z_step=float(np.median(np.diff(lv))) if lv.size > 1 else nan, z_range=float(np.ptp(C[:, 2])),
        rms_full=rms_radius(C), rms_core=rms_radius(C[core]) if ok else nan,
        rms_xy_full=rms_radius(C[:, :2]), sd_z_full=float(C[:, 2].std()),
        rms_xy_core=rms_radius(C[core, :2]) if ok else nan, sd_z_core=float(C[core, 2].std()) if ok else nan,
        core_thin_axis_cos_z=thin_axis_cos_z(C[core]) if ok else nan)


def stage_geometry(path: str | Path, t: float, name: str | None = None, core_key: str = "core",
                   allow_unmapped: bool = False) -> StageGeometry:
    """Load a released file with vec.prep.load_stage and measure it.

    Refuses files whose cell-type labels are missing from configs/lineage_map.yaml: an unmapped label silently
    becomes non-core, which would change the core-cardiac RMS. Map every new label (and commit the map) first."""
    from vec.prep import load_stage

    a = load_stage(path)
    unmapped = list(a.uns.get("unmapped_labels", []))
    if unmapped and not allow_unmapped:
        raise ValueError(f"{path}: labels not in configs/lineage_map.yaml: {unmapped}")
    C = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)[:, :3]
    return geometry_from_arrays(name or Path(path).stem, t, C, np.asarray(a.obs[core_key], dtype=bool))


def with_common_step(g: StageGeometry, step_ref: float) -> StageGeometry:
    """Rescale a stage's z to the reference stage's step (one-pitch hypothesis); xy is unchanged."""
    f = step_ref / g.z_step
    out = StageGeometry(**asdict(g))
    out.sd_z_full, out.sd_z_core, out.z_range = g.sd_z_full * f, g.sd_z_core * f, g.z_range * f
    out.rms_full = float(np.hypot(g.rms_xy_full, out.sd_z_full))
    out.rms_core = float(np.hypot(g.rms_xy_core, out.sd_z_core))
    out.z_step = step_ref
    return out


# ----------------------------------------------------------------------------------------------- growth rules
# Each rule maps training times ts, values ys (log size) and a target time to a predicted log size at the target.
def _rate(ts, ys, i, j):
    return (ys[j] - ys[i]) / (ts[j] - ts[i])


def rule_copy_last(ts, ys, t_star):
    return float(ys[-1])


def rule_latest_rate(ts, ys, t_star):
    """Log-linear growth at the rate of the last interval."""
    return float(ys[-1] + _rate(ts, ys, -2, -1) * (t_star - ts[-1]))


def rule_half_rate(ts, ys, t_star):
    """Half the latest log rate (strong deceleration)."""
    return float(ys[-1] + 0.5 * _rate(ts, ys, -2, -1) * (t_star - ts[-1]))


def rule_power_t(ts, ys, t_star):
    """log size linear in log t through the last two stages (relative growth rate falls as 1/t)."""
    b = (ys[-1] - ys[-2]) / (np.log(ts[-1]) - np.log(ts[-2]))
    return float(ys[-1] + b * (np.log(t_star) - np.log(ts[-1])))


def rule_decaying_rate(ts, ys, t_star):
    """Exponentially decaying log rate (Gompertz), decay fitted on the last two interval rates.

    Falls back to the latest rate with fewer than three stages or when the rates do not show positive, decelerating
    growth (r_prev > r_last > 0 fails)."""
    if len(ts) < 3:
        return rule_latest_rate(ts, ys, t_star)
    ra, rb = _rate(ts, ys, -3, -2), _rate(ts, ys, -2, -1)
    if not (ra > rb > 0):
        return rule_latest_rate(ts, ys, t_star)
    ma, mb = 0.5 * (ts[-3] + ts[-2]), 0.5 * (ts[-2] + ts[-1])
    lam = np.log(ra / rb) / (mb - ma)
    # integral of rb * exp(-lam (t - mb)) from ts[-1] to t_star, written with expm1 so lam -> 0 is exact
    growth = rb * np.exp(-lam * (ts[-1] - mb)) * (-np.expm1(-lam * (t_star - ts[-1]))) / lam
    return float(ys[-1] + growth)


RULES: dict[str, Callable] = {
    "copy_last": rule_copy_last,
    "latest_rate": rule_latest_rate,
    "half_rate": rule_half_rate,
    "power_t": rule_power_t,
    "decaying_rate": rule_decaying_rate,
}


def applicable_rules(ts: Sequence[float], rules: Sequence[str], min_start_t: float) -> list[str]:
    """Rules whose rate intervals all start at or after min_start_t (copy_last is always applicable)."""
    ts = np.asarray(ts, float)
    need = {"copy_last": 0, "latest_rate": 2, "half_rate": 2, "power_t": 2, "decaying_rate": 3}
    return [r for r in rules if need[r] == 0 or (len(ts) >= need[r] and ts[-need[r]] >= min_start_t)]


def pseudo_targets(ts: Sequence[float], ys: Sequence[float], rules: Sequence[str], min_start_t: float = -np.inf,
                   names: Sequence[str] | None = None) -> list[dict]:
    """Leave-future-out errors: predict stage j from stages < j ("next") and from stages < j-1 ("skip").
    A (pseudo-target, rule) pair is eligible for selection when the rule is applicable on that training set, i.e.
    every rate interval it uses starts at or after min_start_t."""
    ts, ys = np.asarray(ts, float), np.asarray(ys, float)
    names = list(names) if names is not None else [f"t{t:g}" for t in ts]
    out = []
    for j in range(2, len(ts)):
        for kind, end in (("next", j), ("skip", j - 1)):
            if end < 2:
                continue
            for r in rules:
                pred = RULES[r](ts[:end], ys[:end], ts[j])
                out.append({"target": names[j], "t_target": float(ts[j]), "kind": kind,
                            "train": ",".join(names[:end]), "rule": r, "pred": pred, "obs": float(ys[j]),
                            "error": pred - float(ys[j]),
                            "eligible": bool(applicable_rules(ts[:end], [r], min_start_t))})
    return out


def select_rules(pt: list[dict], rules: Sequence[str], tie_tol: float) -> tuple[list[str], dict]:
    """Among rules with at least one eligible pseudo-target, those within tie_tol of the lowest mean |error|.
    If no rule could be scored, every rule in `rules` is returned (averaged by the caller)."""
    mae = {}
    for r in rules:
        errs = [abs(p["error"]) for p in pt if p["rule"] == r and p["eligible"]]
        mae[r] = float(np.mean(errs)) if errs else float("nan")
    scored = [r for r in rules if np.isfinite(mae[r])]
    if not scored:
        return list(rules), mae
    best = min(mae[r] for r in scored)
    return [r for r in scored if mae[r] <= best + tie_tol], mae


def predict_log(ts, ys, t_star, chosen: Sequence[str]) -> float:
    return float(np.mean([RULES[r](np.asarray(ts, float), np.asarray(ys, float), t_star) for r in chosen]))


# ----------------------------------------------------------------------------------------------- scope
def scope_log_ratio(s: Sequence[float], ts: Sequence[float], t_star: float, hypothesis: str) -> float:
    """Scope log-ratio log(full/core) at t_star from the per-stage values s (sorted by time)."""
    s, ts = np.asarray(s, float), np.asarray(ts, float)
    if hypothesis == "H0":
        return float(s[-1])
    if hypothesis == "H2":
        return float(0.5 * s[-1])
    if hypothesis == "H1":
        if len(s) >= 2 and s[-2] > s[-1] > 0:
            kappa = np.log(s[-2] / s[-1]) / (ts[-1] - ts[-2])
            return float(s[-1] * np.exp(-kappa * (t_star - ts[-1])))
        return float(s[-1])
    raise ValueError(f"unknown scope hypothesis {hypothesis!r}")


# ----------------------------------------------------------------------------------------------- procedures
def resolve_heart_z_mode(geoms: Sequence[StageGeometry], z_mode: str) -> str:
    if z_mode != "auto":
        return z_mode
    lim = np.cos(np.radians(PREREGISTERED["auto_thin_axis_max_angle_deg"]))
    last_two = sorted(geoms, key=lambda g: g.t)[-2:]
    return "xy_grow_z_hold" if all(g.core_thin_axis_cos_z >= lim for g in last_two) else "raw"


def _grow(ts, ys, names, target_t, rules, min_start_t, tie_tol, tag, out):
    pt = pseudo_targets(ts, ys, list(RULES), min_start_t, names)
    cand = applicable_rules(ts, rules, min_start_t)
    out["valid"] = bool(cand)
    out[f"applicable_{tag}"] = cand
    chosen, mae = select_rules(pt, cand or list(rules), tie_tol)
    out[f"pseudo_targets_{tag}"] = pt
    out[f"mae_{tag}"] = mae
    out[f"chosen_{tag}"] = chosen
    out[f"all_rule_predictions_{tag}"] = {r: RULES[r](ts, ys, target_t) for r in RULES}
    return predict_log(ts, ys, target_t, chosen)


def heart_size_from_geometry(geoms: Sequence[StageGeometry], target_t: float, z_mode: str | None = None,
                             rules: Sequence[str] | None = None, min_start_t: float | None = None,
                             tie_tol: float | None = None) -> dict:
    """The heart procedure on measured stage geometries (any order; the latest stage is the reference).
    Check result["valid"]: False means no candidate rule had a rate interval inside the growth regime, and the
    numbers (an average of all rules, kept for diagnosis) must not be used."""
    rules = list(rules or PREREGISTERED["heart_rules"])
    min_start_t = PREREGISTERED["heart_min_rate_start_t"] if min_start_t is None else min_start_t
    tie_tol = PREREGISTERED["tie_tol"] if tie_tol is None else tie_tol
    geoms = sorted(geoms, key=lambda g: g.t)
    requested = z_mode or PREREGISTERED["heart_z_mode"]
    mode = resolve_heart_z_mode(geoms, requested)
    if mode == "common_step":
        geoms = [with_common_step(g, geoms[-1].z_step) for g in geoms]
    ts = np.array([g.t for g in geoms])
    names = [g.name for g in geoms]
    ref = geoms[-1]
    out = {"target_t": target_t, "z_mode_requested": requested, "z_mode": mode, "rules": rules,
           "min_start_t": min_start_t, "tie_tol": tie_tol, "reference": ref.name,
           "stages": [asdict(g) | {"log_scope": g.log_scope, "log_scope_xy": g.log_scope_xy} for g in geoms]}
    hyps = PREREGISTERED["heart_scope_hypotheses"]
    out["scope"] = {}
    if mode in ("raw", "common_step"):
        y = _grow(ts, np.log([g.rms_core for g in geoms]), names, target_t, rules, min_start_t, tie_tol, "core", out)
        core = float(np.exp(y))
        for h in hyps:
            s = scope_log_ratio([g.log_scope for g in geoms], ts, target_t, h)
            out["scope"][h] = {"log_scope": s, "pred_rms_full": core * float(np.exp(s))}
    elif mode == "xy_grow_z_hold":
        y = _grow(ts, np.log([g.rms_xy_core for g in geoms]), names, target_t, rules, min_start_t, tie_tol,
                  "core_xy", out)
        core_xy = float(np.exp(y))
        core = float(np.hypot(core_xy, ref.sd_z_core))
        out["pred_rms_xy_core"] = core_xy
        for h in hyps:
            s = scope_log_ratio([g.log_scope_xy for g in geoms], ts, target_t, h)
            full_xy = core_xy * float(np.exp(s))
            out["scope"][h] = {"log_scope_xy": s, "pred_rms_xy_full": full_xy,
                               "pred_rms_full": float(np.hypot(full_xy, ref.sd_z_full))}
    elif mode == "split":
        yx = _grow(ts, np.log([g.rms_xy_core for g in geoms]), names, target_t, rules, min_start_t, tie_tol,
                   "core_xy", out)
        yz = _grow(ts, np.log([g.sd_z_core for g in geoms]), names, target_t, rules, min_start_t, tie_tol,
                   "core_z", out)
        core = float(np.hypot(np.exp(yx), np.exp(yz)))
        for h in hyps:
            s = scope_log_ratio([g.log_scope for g in geoms], ts, target_t, h)
            out["scope"][h] = {"log_scope": s, "pred_rms_full": core * float(np.exp(s))}
    else:
        raise ValueError(f"unknown z_mode {mode!r}")
    out["pred_rms_core"] = core
    out["pred_log_growth_core"] = float(np.log(core / ref.rms_core))
    for v in out["scope"].values():
        v["scale_vs_reference_full"] = v["pred_rms_full"] / ref.rms_full
    return out


def heart_size_procedure(stages: Sequence[tuple[str | Path, float]], target_t: float, z_mode: str | None = None,
                         **kw) -> dict:
    """Heart E12.5 size from heart training files. stages: [(path, nominal E-day), ...], any order."""
    geoms = [stage_geometry(p, t, core_key=PREREGISTERED["core_key"]) for p, t in stages]
    return heart_size_from_geometry(geoms, target_t, z_mode=z_mode, **kw)


def embryo_size_from_geometry(geoms: Sequence[StageGeometry], target_t: float, z_mode: str | None = None) -> dict:
    z_mode = z_mode or PREREGISTERED["embryo_z_mode"]
    geoms = sorted(geoms, key=lambda g: g.t)
    if z_mode == "common_step":
        geoms = [with_common_step(g, geoms[-1].z_step) for g in geoms]
    elif z_mode != "raw":
        raise ValueError("embryo procedure supports z_mode 'raw' or 'common_step'")
    ts = np.array([g.t for g in geoms])
    ys = np.log([g.rms_full for g in geoms])
    if not ts.min() < target_t < ts.max():
        raise ValueError("the embryo target must lie strictly inside the training stages")
    f = PchipInterpolator(ts, ys, extrapolate=False)
    lo = int(np.searchsorted(ts, target_t)) - 1
    out = {"target_t": target_t, "z_mode": z_mode, "interpolant": PREREGISTERED["embryo_interpolant"],
           "stages": [asdict(g) for g in geoms], "pred_rms_full": float(np.exp(f(target_t))),
           "bracket": [geoms[lo].name, geoms[lo + 1].name]}
    out["scale_vs_lower_bracket"] = out["pred_rms_full"] / geoms[lo].rms_full
    out["scale_vs_upper_bracket"] = out["pred_rms_full"] / geoms[lo + 1].rms_full
    loo = []
    for j in range(1, len(ts) - 1):
        keep = np.arange(len(ts)) != j
        fj = PchipInterpolator(ts[keep], ys[keep], extrapolate=False)
        loo.append({"left_out": geoms[j].name, "t": float(ts[j]), "pred": float(fj(ts[j])), "obs": float(ys[j]),
                    "error": float(fj(ts[j]) - ys[j])})
    out["loo_interior"] = loo
    return out


def embryo_size_procedure(stages: Sequence[tuple[str | Path, float]], target_t: float, z_mode: str | None = None
                          ) -> dict:
    """Embryo E7.75 size from embryo training files. stages: [(path, nominal E-day), ...], any order."""
    geoms = [stage_geometry(p, t, core_key=PREREGISTERED["core_key"]) for p, t in stages]
    return embryo_size_from_geometry(geoms, target_t, z_mode=z_mode)


def rescale_to_rms(C: np.ndarray, target_rms: float) -> np.ndarray:
    """Scale coordinates about their centroid so their raw 3D RMS radius equals target_rms (no rotation, no mirror)."""
    C = np.asarray(C, float)
    mu = C.mean(0)
    return mu + (C - mu) * (target_rms / rms_radius(C))


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Run the pre-registered size procedure on released training files.")
    ap.add_argument("--board", choices=["heart", "embryo"], required=True)
    ap.add_argument("--stage", nargs=2, action="append", metavar=("PATH", "E_DAY"), required=True)
    ap.add_argument("--target", type=float, required=True)
    ap.add_argument("--z-mode", default=None)
    ap.add_argument("--out", default=None, help="write the full result as JSON here")
    a = ap.parse_args()
    stages = [(p, float(t)) for p, t in a.stage]
    fn = heart_size_procedure if a.board == "heart" else embryo_size_procedure
    res = fn(stages, a.target, z_mode=a.z_mode)
    txt = json.dumps(res, indent=2, default=float)
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(txt)
    keys = ("valid", "z_mode", "pred_rms_core", "pred_log_growth_core", "scope", "pred_rms_full", "bracket",
            "loo_interior")
    print(json.dumps({k: res[k] for k in keys if k in res}, indent=2, default=float))
