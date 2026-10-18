"""Pre-registered size procedure for the two T2 test boards: heart E12.5 (extrapolation) and embryo E7.75
(interpolation).

Rules this file obeys:
- Every number is computed from the released training files passed in. There are no size numbers and no assumed growth
  rates in this file: the numeric constants are E-days, an angle, tolerances and the forms of the rules (for example
  "half" the latest rate). test_size_procedure.py checks every numeric literal against an explicit allow-list.
- It never reads boards/index.json anchors, leaderboard scores, the PyPI veckit docstrings or the veckit pull request
  that lists radii, and it takes no external data. Nominal stage times (E-days) are inputs, from the file names.
- RMS radius is the scorer's definition (veckit common.shape_metrics.scale_log_ratio): sqrt(mean squared distance to
  the centroid) over all three raw coordinates, in each file's own convention.
- Duplicated records (identical count rows at the same xy in the same or an adjacent z level) are dropped before any
  geometry is measured (vec.prep.load_stage(..., drop_duplicates=True)); the numbers dropped are logged per stage.
- Core-cardiac cells are obs["core"] from vec.prep.load_stage at the committed configs/lineage_map.yaml
  (cardiomyocytes and endocardium). A file with labels missing from the map is refused.
- The Mab21l2 knockout (T3) is not an input: cross-task use is an open question.

Why z gets special treatment (experiments/coords/REPORT.md, parts 1-2): every released stack holds a similar number of
sections whatever the stage, cells fill it along z, and tissue reaches the end sections. So the z extent is set by how
many sections were imaged and by the stored step, not by tissue growth, and in-plane (xy) growth is the anatomical
signal when the sectioning plane is comparable.

Heart E12.5 (`heart_size_procedure`):
  1. For each heart training stage: geometry of the core-cardiac cells and of the full cloud: raw 3D RMS, in-plane RMS
     (x, y), SD of z, and how close the core's thinnest principal axis is to z.
  2. z mode (PREREGISTERED["heart_z_mode"] = "auto"): "xy_grow_z_hold" when the core's thinnest principal axis is
     within 30 degrees of z in both of the last two stages (comparable sectioning plane), otherwise "raw". All four
     modes are computed and logged next to the pre-registered one:
       raw:            grow the core's 3D RMS; scope ratio on 3D RMS.
       xy_grow_z_hold: grow the core's in-plane RMS; scope ratio on in-plane RMS; hold the core and full-cloud z SD at
                       the reference stage's values; recombine in quadrature.
       common_step:    as raw after rescaling every stage's z to the reference stage's step (one-pitch hypothesis).
       split:          grow the in-plane RMS and the z SD separately; recombine in quadrature.
  3. Growth: the fixed rules latest_rate, half_rate and power_t (PREREGISTERED["heart_rules"]), each using only rate
     intervals that start at or after heart_min_rate_start_t (E9.0). Before that the core is remodelled (crescent,
     tube, looping): its RMS is flat and then jumps, so rates measured there say nothing about later growth (REPORT.md,
     part 4). If no rule is applicable the result is flagged valid=False. Applicable rules are compared on
     leave-future-out pseudo-targets (predict stage j from the stages before it, and from the stages before j-1),
     using only pseudo-targets on which every applicable rule is eligible, i.e. whose last training interval also
     starts at or after E9.0. Rules within tie_tol of the lowest mean |log error| are averaged in log space; with no
     such pseudo-target every applicable rule is averaged. With the stages released by 2026-10-20 the only interval
     in range is E9.5 to E10.5, no pseudo-target is eligible, and the prediction is the average of the three rules.

     Why these three, averaged (rule_table(); results/size_rule_table.csv). Each candidate rule predicts the log growth
     from the reference stage to the target as m times the observed log growth G of the last interval. Let rho be the
     later mean log rate relative to the observed one; the true growth is h rho G, with h = (target - reference) /
     (reference - previous), so the bias is (m - h rho) G. Two properties follow without assuming any rate value: the
     worst-case bias coefficient, max over rho in [0.5, 1] of |m - h rho| (later growth no faster than observed, and at
     least half as fast), and the noise amplification sqrt((1 + m)^2 + m^2 + 1), the SD of the prediction error in
     units of the per-embryo SD when the previous, reference and target embryos vary independently. For the
     2026-10-20 layout (E9.5 and E10.5 to E12.5, h = 2):

         rule          m       worst-case bias   noise amplification
         copy_last     0       2                 1.414
         half_rate     1       1                 2.449
         power_t       1.742   0.742             3.399
         latest_rate   2       1                 3.742
         average       1.581   0.581             3.187

     The average of half_rate, power_t and latest_rate has the smallest worst-case bias of these rules, at a higher
     noise amplification than half_rate or copy_last. Which term dominates depends on the later rate and the
     per-embryo SD, and neither is assumed.
  4. Scope s = log(full / core) at the target under three pre-registered hypotheses:
       H1 "trend": s decays exponentially at the rate fitted on the last two stages (held if not decreasing);
       H0 "hold":  the reference (latest) stage's s;
       H2 "half":  half the reference stage's s (heart-focused scope, shell influence halved).
  Output: predicted core and full-cloud RMS per scope hypothesis. A submission rescales its coordinates about the
  centroid to the full-cloud RMS of the scope hypothesis it implements (`rescale_to_rms`; proper scaling only).
Embryo E7.75 (`embryo_size_procedure`):
  PCHIP (monotone cubic) interpolation of log full-cloud RMS against t through every embryo training stage (z mode
  "raw"), with leave-one-out errors on the interior stages as the check. The other three z modes are logged:
  common_step as for the heart, xy_grow_z_hold (in-plane RMS interpolated, z SD held at the lower bracketing stage) and
  split (in-plane RMS and z SD interpolated separately). Extrapolation is refused.

Command line (`run`): writes one JSON with the pre-registered result, every z mode (and, for the heart, every scope
hypothesis), the valid flag, the stage geometry and the duplicates dropped per stage, the git commit, the sha256 of
every input file, of this file, of the lineage map and of the loader code, and the timestamp given with --timestamp.
`--redact-forecasts` replaces every target forecast with a marker (for dry runs on incomplete stage sets, which are code
checks, not predictions).

Refit on 2026-10-20 (experiments/coords/run_p3.sh): the same functions with the newly released stages added (heart:
E8.5, E10.5; embryo: E7.5). The forms, rules, selection settings, z-mode switch and scope hypotheses stay as here.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

PREREGISTERED = {
    "heart_rules": ["latest_rate", "half_rate", "power_t"],
    "heart_min_rate_start_t": 9.0,
    "tie_tol": 0.15,
    "heart_scope_hypotheses": ["H1", "H0", "H2"],
    "heart_z_mode": "auto",
    "auto_thin_axis_max_angle_deg": 30.0,
    "embryo_interpolant": "pchip",
    "embryo_z_mode": "raw",
    "core_key": "core",
    "drop_duplicates": True,
}
Z_MODES = ("raw", "common_step", "xy_grow_z_hold", "split")
RHO_RANGE = (0.5, 1.0)  # later mean log rate / observed log rate: no faster than observed, at least half as fast
TABLE_RULES = ("copy_last", "half_rate", "power_t", "latest_rate")
MIN_CORE_CELLS = 4
REPO = Path(__file__).resolve().parents[2]
PROVENANCE_FILES = ("configs/lineage_map.yaml", "src/vec/prep.py", "src/vec/spatial.py",
                    "experiments/coords/handedness.py")
REDACTED = "[redacted: dry-run forecast]"


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

    @property
    def core_thin_axis_angle_deg(self) -> float:
        return float(np.degrees(np.arccos(min(self.core_thin_axis_cos_z, 1.0))))


def geometry_from_arrays(name: str, t: float, C: np.ndarray, core: np.ndarray) -> StageGeometry:
    C = np.asarray(C, float)[:, :3]
    core = np.asarray(core, bool)
    lv = z_levels(C[:, 2])
    ok = core.sum() >= MIN_CORE_CELLS
    nan = float("nan")
    return StageGeometry(
        name=name, t=float(t), n_cells=int(C.shape[0]), n_core=int(core.sum()), n_levels=int(lv.size),
        z_step=float(np.median(np.diff(lv))) if lv.size > 1 else nan, z_range=float(np.ptp(C[:, 2])),
        rms_full=rms_radius(C), rms_core=rms_radius(C[core]) if ok else nan,
        rms_xy_full=rms_radius(C[:, :2]), sd_z_full=float(C[:, 2].std()),
        rms_xy_core=rms_radius(C[core, :2]) if ok else nan, sd_z_core=float(C[core, 2].std()) if ok else nan,
        core_thin_axis_cos_z=thin_axis_cos_z(C[core]) if ok else nan)


def measure_stage(path: str | Path, t: float, name: str | None = None, core_key: str | None = None,
                  allow_unmapped: bool = False) -> tuple[StageGeometry, dict]:
    """Load a released file with vec.prep.load_stage (duplicated records dropped) and measure it.

    Refuses files whose cell-type labels are missing from configs/lineage_map.yaml: an unmapped label silently
    becomes non-core, which would change the core-cardiac RMS. Map every new label (and commit the map) first.
    Returns (geometry, record of the load: duplicates dropped, core labels and their cell counts)."""
    from vec.prep import load_stage

    core_key = core_key or PREREGISTERED["core_key"]
    a = load_stage(path, drop_duplicates=PREREGISTERED["drop_duplicates"])
    unmapped = list(a.uns.get("unmapped_labels", []))
    if unmapped and not allow_unmapped:
        raise ValueError(f"{path}: labels not in configs/lineage_map.yaml: {unmapped}")
    C = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)[:, :3]
    core = np.asarray(a.obs[core_key], dtype=bool)
    g = geometry_from_arrays(name or Path(path).stem, t, C, core)
    qc = a.uns.get("qc") or {}
    labels = a.obs["celltype"].astype(str) if "celltype" in a.obs else None
    core_labels = {} if labels is None else {str(k): int(v) for k, v in labels[core].value_counts().items()}
    info = {
        "name": g.name, "path": str(path), "t": float(t),
        "n_input": int(qc.get("n_input", a.n_obs)), "n_dropped_duplicates": int(qc.get("n_dropped", 0)),
        "dropped_by_level_gap": {str(k): int(v) for k, v in dict(qc.get("by_gap", {})).items()},
        "dropped_max_xy_offset": float(qc.get("max_xy_offset", 0.0)),
        "duplicate_matches_below_min_library": int(qc.get("below_min_library", 0)),
        "n_cells_measured": g.n_cells, "n_core": g.n_core, "core_labels": core_labels, "unmapped_labels": unmapped,
    }
    return g, info


def stage_geometry(path: str | Path, t: float, name: str | None = None, core_key: str = "core",
                   allow_unmapped: bool = False) -> StageGeometry:
    """Geometry only (see measure_stage)."""
    return measure_stage(path, t, name=name, core_key=core_key, allow_unmapped=allow_unmapped)[0]


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
    """Diagnostic only (not a pre-registered rule): exponentially decaying log rate, decay fitted on the last two
    interval rates. Falls back to the latest rate with fewer than three stages or without positive, decelerating
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
RULE_NEEDS = {"copy_last": 0, "latest_rate": 2, "half_rate": 2, "power_t": 2, "decaying_rate": 3}


def applicable_rules(ts: Sequence[float], rules: Sequence[str], min_start_t: float) -> list[str]:
    """Rules whose rate intervals all start at or after min_start_t (copy_last is always applicable)."""
    ts = np.asarray(ts, float)
    return [r for r in rules if RULE_NEEDS[r] == 0 or (len(ts) >= RULE_NEEDS[r] and ts[-RULE_NEEDS[r]] >= min_start_t)]


def rule_multiplier(rule: str, t_prev: float, t_ref: float, t_target: float) -> float:
    """m such that the rule's predicted log growth from t_ref to t_target is m times the observed log growth from
    t_prev to t_ref. Exact for the candidate rules, which are linear in the observed log sizes."""
    ts = np.array([t_prev, t_ref])
    return float(RULES[rule](ts, np.array([0.0, 1.0]), t_target) - RULES[rule](ts, np.array([0.0, 0.0]), t_target)
                 - 1.0)


def rule_table(t_prev: float, t_ref: float, t_target: float, rules: Sequence[str] | None = None,
               rho_range: tuple[float, float] = RHO_RANGE) -> list[dict]:
    """Rate-free comparison of the candidate rules (see the module docstring): per rule and for the log-space average
    of `rules` (default PREREGISTERED["heart_rules"]): multiplier m on the observed log growth, horizon ratio h,
    worst-case bias coefficient max over rho in rho_range of |m - h rho|, and noise amplification
    sqrt((1 + m)^2 + m^2 + 1) in units of the per-embryo SD."""
    rules = list(rules or PREREGISTERED["heart_rules"])
    h = (t_target - t_ref) / (t_ref - t_prev)
    ms = {r: rule_multiplier(r, t_prev, t_ref, t_target) for r in TABLE_RULES}
    ms["average"] = float(np.mean([rule_multiplier(r, t_prev, t_ref, t_target) for r in rules]))
    out = []
    for name, m in ms.items():
        out.append({"rule": name, "averaged_rules": "+".join(rules) if name == "average" else name,
                    "t_prev": float(t_prev), "t_ref": float(t_ref), "t_target": float(t_target), "horizon_ratio": h,
                    "multiplier_m": m, "rho_min": float(rho_range[0]), "rho_max": float(rho_range[1]),
                    "worst_case_bias_coef": float(max(abs(m - h * rho_range[0]), abs(m - h * rho_range[1]))),
                    "noise_amplification": float(np.sqrt((1 + m) ** 2 + m ** 2 + 1))})
    return out


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


def select_rules(pt: list[dict], rules: Sequence[str], tie_tol: float) -> tuple[list[str], dict, list[str]]:
    """Compare `rules` on the pseudo-targets where every one of them is eligible (one common set, so each mean |error|
    is over the same targets) and keep those within tie_tol of the lowest. With no common eligible pseudo-target every
    rule is kept (averaged by the caller). Returns (kept rules, mean |error| per rule, the common pseudo-targets)."""
    rules = list(rules)
    elig = {(p["target"], p["kind"], p["rule"]): p for p in pt if p["eligible"]}
    keys = sorted({(p["target"], p["kind"]) for p in pt})
    common = [k for k in keys if all((k[0], k[1], r) in elig for r in rules)]
    mae = {r: float(np.mean([abs(elig[(k[0], k[1], r)]["error"]) for k in common])) if common else float("nan")
           for r in rules}
    labels = [f"{k[0]}:{k[1]}" for k in common]
    if not common:
        return rules, mae, labels
    best = min(mae.values())
    return [r for r in rules if mae[r] <= best + tie_tol], mae, labels


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


def _grow(ts, ys, names, target_t, rules, min_start_t, tie_tol) -> dict:
    """Growth of one series: its own validity, applicable and chosen rules, and the predicted log value."""
    pt = pseudo_targets(ts, ys, list(RULES), min_start_t, names)
    cand = applicable_rules(ts, rules, min_start_t)
    chosen, mae, common = select_rules(pt, cand or list(rules), tie_tol)
    return {"valid": bool(cand), "applicable": cand, "chosen": chosen, "mae_common": mae,
            "common_pseudo_targets": common, "pseudo_targets": pt,
            "all_rule_predictions": {r: RULES[r](ts, ys, target_t) for r in RULES},
            "pred_log_target": predict_log(ts, ys, target_t, chosen)}


def heart_size_from_geometry(geoms: Sequence[StageGeometry], target_t: float, z_mode: str | None = None,
                             rules: Sequence[str] | None = None, min_start_t: float | None = None,
                             tie_tol: float | None = None) -> dict:
    """The heart procedure on measured stage geometries (any order; the latest stage is the reference).
    Check result["valid"]: False means some grown series had no candidate rule with a rate interval inside the
    growth regime, and the numbers (an average of all rules, kept for diagnosis) must not be used."""
    rules = list(rules or PREREGISTERED["heart_rules"])
    min_start_t = PREREGISTERED["heart_min_rate_start_t"] if min_start_t is None else min_start_t
    tie_tol = PREREGISTERED["tie_tol"] if tie_tol is None else tie_tol
    geoms = sorted(geoms, key=lambda g: g.t)
    requested = z_mode or PREREGISTERED["heart_z_mode"]
    mode = resolve_heart_z_mode(geoms, requested)
    if mode not in Z_MODES:
        raise ValueError(f"unknown z_mode {mode!r}")
    if mode == "common_step":
        geoms = [with_common_step(g, geoms[-1].z_step) for g in geoms]
    ts = np.array([g.t for g in geoms])
    names = [g.name for g in geoms]
    ref = geoms[-1]
    out = {"target_t": target_t, "z_mode_requested": requested, "z_mode": mode, "rules": rules,
           "min_start_t": min_start_t, "tie_tol": tie_tol, "reference": ref.name,
           "stages": [asdict(g) | {"log_scope": g.log_scope, "log_scope_xy": g.log_scope_xy,
                                   "core_thin_axis_angle_deg": g.core_thin_axis_angle_deg} for g in geoms],
           "growth": {}, "scope": {}}

    def grow(tag: str, values) -> float:
        out["growth"][tag] = _grow(ts, np.log(values), names, target_t, rules, min_start_t, tie_tol)
        return float(np.exp(out["growth"][tag]["pred_log_target"]))

    hyps = PREREGISTERED["heart_scope_hypotheses"]
    if mode in ("raw", "common_step"):
        core = grow("core", [g.rms_core for g in geoms])
        for h in hyps:
            s = scope_log_ratio([g.log_scope for g in geoms], ts, target_t, h)
            out["scope"][h] = {"pred_log_scope": s, "pred_rms_full": core * float(np.exp(s))}
    elif mode == "xy_grow_z_hold":
        core_xy = grow("core_xy", [g.rms_xy_core for g in geoms])
        core = float(np.hypot(core_xy, ref.sd_z_core))
        out["pred_rms_xy_core"] = core_xy
        for h in hyps:
            s = scope_log_ratio([g.log_scope_xy for g in geoms], ts, target_t, h)
            full_xy = core_xy * float(np.exp(s))
            out["scope"][h] = {"pred_log_scope_xy": s, "pred_rms_xy_full": full_xy,
                               "pred_rms_full": float(np.hypot(full_xy, ref.sd_z_full))}
    else:  # split
        core_xy = grow("core_xy", [g.rms_xy_core for g in geoms])
        core_z = grow("core_z", [g.sd_z_core for g in geoms])
        core = float(np.hypot(core_xy, core_z))
        out["pred_rms_xy_core"], out["pred_sd_z_core"] = core_xy, core_z
        for h in hyps:
            s = scope_log_ratio([g.log_scope for g in geoms], ts, target_t, h)
            out["scope"][h] = {"pred_log_scope": s, "pred_rms_full": core * float(np.exp(s))}
    out["valid"] = all(g["valid"] for g in out["growth"].values())
    out["pred_rms_core"] = core
    out["pred_log_growth_core"] = float(np.log(core / ref.rms_core))
    for v in out["scope"].values():
        v["scale_vs_reference_full"] = v["pred_rms_full"] / ref.rms_full
    return out


def _mode_summary(res: dict, preregistered_mode: str) -> dict:
    return {"z_mode": res["z_mode"], "preregistered": res["z_mode"] == preregistered_mode, "valid": res["valid"],
            "chosen": {k: g["chosen"] for k, g in res.get("growth", {}).items()},
            "pred_rms_core": res.get("pred_rms_core"), "pred_log_growth_core": res.get("pred_log_growth_core"),
            "pred_rms_full": ({h: v["pred_rms_full"] for h, v in res["scope"].items()} if "scope" in res
                              else res.get("pred_rms_full"))}


def heart_report(geoms: Sequence[StageGeometry], target_t: float) -> dict:
    """The pre-registered heart result plus every z mode (each with all scope hypotheses) and the rule table."""
    pre = heart_size_from_geometry(geoms, target_t)
    modes = {m: heart_size_from_geometry(geoms, target_t, z_mode=m) for m in Z_MODES}
    ts = sorted(g.t for g in geoms)
    table = rule_table(ts[-2], ts[-1], target_t) if len(ts) >= 2 else []
    return {"valid": pre["valid"], "preregistered": pre,
            "summary": [_mode_summary(modes[m], pre["z_mode"]) for m in Z_MODES], "z_modes": modes,
            "rule_table": table}


def heart_size_procedure(stages: Sequence[tuple[str | Path, float]], target_t: float, z_mode: str | None = None,
                         **kw) -> dict:
    """Heart E12.5 size from heart training files. stages: [(path, nominal E-day), ...], any order."""
    geoms = [stage_geometry(p, t, core_key=PREREGISTERED["core_key"]) for p, t in stages]
    return heart_size_from_geometry(geoms, target_t, z_mode=z_mode, **kw)


def _embryo_interp(geoms: Sequence[StageGeometry], t: float, z_mode: str) -> float:
    ts = np.array([g.t for g in geoms])

    def pchip(values) -> float:
        return float(np.exp(PchipInterpolator(ts, np.log(values), extrapolate=False)(t)))

    if z_mode in ("raw", "common_step"):
        return pchip([g.rms_full for g in geoms])
    xy = pchip([g.rms_xy_full for g in geoms])
    if z_mode == "xy_grow_z_hold":
        return float(np.hypot(xy, geoms[int(np.searchsorted(ts, t)) - 1].sd_z_full))
    return float(np.hypot(xy, pchip([g.sd_z_full for g in geoms])))


def embryo_size_from_geometry(geoms: Sequence[StageGeometry], target_t: float, z_mode: str | None = None) -> dict:
    z_mode = z_mode or PREREGISTERED["embryo_z_mode"]
    if z_mode not in Z_MODES:
        raise ValueError(f"unknown z_mode {z_mode!r}")
    geoms = sorted(geoms, key=lambda g: g.t)
    if z_mode == "common_step":
        geoms = [with_common_step(g, geoms[-1].z_step) for g in geoms]
    ts = np.array([g.t for g in geoms])
    if not ts.min() < target_t < ts.max():
        raise ValueError("the embryo target must lie strictly inside the training stages")
    lo = int(np.searchsorted(ts, target_t)) - 1
    out = {"target_t": target_t, "z_mode": z_mode, "interpolant": PREREGISTERED["embryo_interpolant"], "valid": True,
           "stages": [asdict(g) for g in geoms], "pred_rms_full": _embryo_interp(geoms, target_t, z_mode),
           "bracket": [geoms[lo].name, geoms[lo + 1].name]}
    out["scale_vs_lower_bracket"] = out["pred_rms_full"] / geoms[lo].rms_full
    out["scale_vs_upper_bracket"] = out["pred_rms_full"] / geoms[lo + 1].rms_full
    loo = []
    for j in range(1, len(ts) - 1):
        rest = [g for i, g in enumerate(geoms) if i != j]
        pred = _embryo_interp(rest, ts[j], z_mode)
        loo.append({"left_out": geoms[j].name, "t": float(ts[j]), "pred": float(np.log(pred)),
                    "obs": float(np.log(geoms[j].rms_full)), "error": float(np.log(pred / geoms[j].rms_full))})
    out["loo_interior"] = loo
    return out


def embryo_report(geoms: Sequence[StageGeometry], target_t: float) -> dict:
    """The pre-registered embryo result plus every z mode."""
    pre = embryo_size_from_geometry(geoms, target_t)
    modes = {m: embryo_size_from_geometry(geoms, target_t, z_mode=m) for m in Z_MODES}
    return {"valid": pre["valid"], "preregistered": pre,
            "summary": [_mode_summary(modes[m], pre["z_mode"]) | {"loo_interior": modes[m]["loo_interior"]}
                        for m in Z_MODES], "z_modes": modes}


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


# ----------------------------------------------------------------------------------------------- run and log
def sha256_file(path: str | Path) -> str:
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def _git(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=REPO, text=True, stderr=subprocess.DEVNULL).rstrip("\n")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def provenance(paths: Sequence[str | Path]) -> dict:
    """git commit, uncommitted changes to the files that define the procedure, and sha256 of inputs and code."""
    here = Path(__file__).resolve()
    tracked = [str(here.relative_to(REPO)), *PROVENANCE_FILES]
    return {"git_commit": _git("rev-parse", "HEAD"),
            "git_uncommitted": [ln for ln in _git("status", "--porcelain", "--", *tracked).splitlines() if ln],
            "size_procedure_sha256": sha256_file(here),
            "code_sha256": {p: sha256_file(REPO / p) for p in PROVENANCE_FILES if (REPO / p).exists()},
            "inputs_sha256": {str(p): sha256_file(p) for p in paths}}


def _stage_names(stages: Sequence[tuple[str | Path, float]]) -> list[str]:
    stems = [Path(p).stem for p, _ in stages]
    return [s if stems.count(s) == 1 else f"{s}@{t:g}" for s, (_, t) in zip(stems, stages, strict=True)]


def run(board: str, stages: Sequence[tuple[str | Path, float]], target_t: float, timestamp: str) -> dict:
    """Measure every stage, run the board's procedure with all z modes, and attach provenance."""
    measured = [measure_stage(p, t, name=n) for (p, t), n in zip(stages, _stage_names(stages), strict=True)]
    geoms = [m[0] for m in measured]
    report = heart_report(geoms, target_t) if board == "heart" else embryo_report(geoms, target_t)
    return {"board": board, "target_t": float(target_t), "timestamp": timestamp, "valid": report["valid"],
            "preregistered_settings": dict(PREREGISTERED), "provenance": provenance([p for p, _ in stages]),
            "inputs": [m[1] for m in measured], **report}


def redact_forecasts(obj):
    """Replace every target forecast (keys starting with pred_ or scale_vs_, and all_rule_predictions) with a marker.
    Pseudo-target and leave-one-out records of released stages ("pred", "obs", "error") are kept."""
    if isinstance(obj, dict):
        return {k: REDACTED if (k.startswith(("pred_", "scale_vs_")) or k == "all_rule_predictions")
                else redact_forecasts(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [redact_forecasts(v) for v in obj]
    return obj


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, list | tuple):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_jsonable(v) for v in o.tolist()]
    if isinstance(o, bool | np.bool_):
        return bool(o)
    if isinstance(o, int | np.integer):
        return int(o)
    if isinstance(o, float | np.floating):
        return float(o)
    return o


def dumps(res: dict) -> str:
    return json.dumps(_jsonable(res), indent=2)


def main(argv: Sequence[str] | None = None) -> dict:
    import argparse

    ap = argparse.ArgumentParser(description="Run the pre-registered size procedure on released training files.")
    ap.add_argument("--board", choices=["heart", "embryo"], required=True)
    ap.add_argument("--stage", nargs=2, action="append", metavar=("PATH", "E_DAY"), required=True)
    ap.add_argument("--target", type=float, required=True)
    ap.add_argument("--timestamp", required=True, help="run time to record in the JSON (an ISO time string)")
    ap.add_argument("--out", default=None, help="write the full result as JSON here")
    ap.add_argument("--redact-forecasts", action="store_true",
                    help="replace target forecasts with a marker (dry runs on incomplete stage sets)")
    a = ap.parse_args(argv)
    res = run(a.board, [(p, float(t)) for p, t in a.stage], a.target, a.timestamp)
    if a.redact_forecasts:
        res = redact_forecasts(res)
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(dumps(res))
    short = {"valid": res["valid"], "board": res["board"], "target_t": res["target_t"], "timestamp": res["timestamp"],
             "git_commit": res["provenance"]["git_commit"], "git_uncommitted": res["provenance"]["git_uncommitted"],
             "duplicates_dropped": {i["name"]: i["n_dropped_duplicates"] for i in res["inputs"]},
             "preregistered_z_mode": res["preregistered"]["z_mode"], "summary": res["summary"]}
    print(dumps(short))
    return res


if __name__ == "__main__":
    main()
