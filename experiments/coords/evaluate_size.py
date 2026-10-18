"""Part 4: evaluate candidate size rules on pseudo-targets built from released training stages only (DGX).

Writes to runs/coords/:
  stage_geometry.csv        per-file geometry used by the size procedure (raw, common-step and isotropy-corrected z)
  size_heart_pseudo.csv     heart core-RMS rules on every pseudo-target available now, with sensitivity variants
  size_heart_scope.csv      scope log-ratio hypotheses on the E9.5 pseudo-target
  size_heart_full.csv       full-cloud pseudo-target errors (core rule x scope hypothesis) and scale skill
  size_heart_modes.csv      the procedure end to end on the E9.5 pseudo-target for every z mode
  size_embryo.csv           embryo interpolation checks (leave-one-out) and scale skill
  size_noise.csv            embryo-to-embryo spread (E9.5 WT vs Mab21l2 KO) and rule noise multipliers
  size_count_density.csv    diagnostic only: RMS / N^(1/3) consistency (needs an organiser ruling before any use)
  size_ceiling.csv          split-half |scale_log_ratio| (the metric's ceiling) per released stage, 10 seeds
  size_dryrun_heart.json, size_dryrun_embryo.json   end-to-end dry runs of size_procedure on current stages
Usage: uv run python experiments/coords/evaluate_size.py
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import size_procedure as sp  # noqa: E402
from coordlib import FILES  # noqa: E402

from vec.evaluate import Stage, _skill, split_half, subsample  # noqa: E402
from vec.prep import load_stage  # noqa: E402

T = {"E6.75": 6.75, "E7.25": 7.25, "E8.0": 8.0, "E8.25_late": 8.25, "E8.75": 8.75, "E9.5": 9.5,
     "E9.5_mab21l2_ko": 9.5}
SEEDS = range(10)


def scale_log_ratio(pred_C, true_C):
    from common.shape_metrics import scale_log_ratio as slr

    return slr(pred_C, true_C)


def iso_geometry(g: sp.StageGeometry, s: float) -> sp.StageGeometry:
    """z divided by the isotropy-collapse stretch s (diagnostic variant only)."""
    out = sp.StageGeometry(**asdict(g))
    out.sd_z_full, out.sd_z_core = g.sd_z_full / s, g.sd_z_core / s
    out.rms_full = float(np.sqrt(g.rms_xy_full ** 2 + out.sd_z_full ** 2))
    out.rms_core = float(np.sqrt(g.rms_xy_core ** 2 + out.sd_z_core ** 2))
    return out


def scale_skill(pred_rms: float, target: Stage, floor_stage: Stage, seeds=SEEDS) -> dict:
    """Server-faithful scale_log_ratio skill of a prediction with raw RMS pred_rms: target subsampled to 10% and
    split in half (truth / ceiling), floor = copy of the floor stage's 10% working copy."""
    sk, vals, fl, ce = [], [], [], []
    for s in seeds:
        tw = subsample(target, 0.1, s)
        truth, other = split_half(tw, seed=0)
        fw = subsample(floor_stage, 0.1, s + 1)
        v = float(np.log(pred_rms / sp.rms_radius(truth.coords)))
        f = scale_log_ratio(fw.coords, truth.coords)
        c = scale_log_ratio(other.coords, truth.coords)
        sk.append(_skill(abs(v), abs(f), abs(c), lower_is_better=True))
        vals.append(v)
        fl.append(f)
        ce.append(c)
    return {"skill_mean": float(np.mean(sk)), "skill_sd": float(np.std(sk, ddof=1)), "value_mean": float(np.mean(vals)),
            "floor_mean": float(np.mean(fl)), "ceiling_abs_mean": float(np.mean(np.abs(ce)))}


def main(out: Path) -> None:
    stages, geoms = {}, {}
    for name, path in FILES.items():
        a = load_stage(path)
        C = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)[:, :3]
        core = np.asarray(a.obs["core"], dtype=bool)
        stages[name] = Stage(X=np.zeros((C.shape[0], 1), np.float32), coords=C)
        geoms[name] = sp.geometry_from_arrays(name, T[name], C, core)
    zfit = pd.read_csv(out / "zpitch_fit.csv")
    s_iso = zfit[zfit.measure == "mism_group"].set_index("stage")["s_best"].to_dict()

    geo_rows = []
    for name, g in geoms.items():
        r = asdict(g) | {"log_scope": g.log_scope if np.isfinite(g.rms_core) else np.nan, "s_iso": s_iso.get(name)}
        gi = iso_geometry(g, s_iso[name])
        r |= {"rms_full_iso": gi.rms_full, "rms_core_iso": gi.rms_core}
        geo_rows.append(r)
    pd.DataFrame(geo_rows).to_csv(out / "stage_geometry.csv", index=False)

    # ---------------- heart: core-RMS rules on E9.5 from (E8.25_late, E8.75) ----------------
    rules = list(sp.RULES)
    rows = []
    heart = ["E8.25_late", "E8.75", "E9.5"]

    def add_rows(variant, gl, target="E9.5", t_override=None):
        ts = np.array([t_override.get(g.name, g.t) if t_override else g.t for g in gl])
        for key, ys in (("core", np.log([g.rms_core for g in gl])), ("full", np.log([g.rms_full for g in gl])),
                        ("xy_core", np.log([g.rms_xy_core for g in gl]))):
            for r in rules:
                pred = sp.RULES[r](ts[:-1], ys[:-1], ts[-1])
                rows.append({"variant": variant, "series": key, "target": target, "rule": r,
                             "train": ",".join(g.name for g in gl[:-1]), "pred_log": pred, "obs_log": ys[-1],
                             "error": pred - ys[-1], "pred_rms": float(np.exp(pred)), "obs_rms": float(np.exp(ys[-1]))})

    add_rows("raw", [geoms[k] for k in heart])
    add_rows("raw_t8.25late=8.375", [geoms[k] for k in heart], t_override={"E8.25_late": 8.375})
    add_rows("raw_t8.25late=8.5", [geoms[k] for k in heart], t_override={"E8.25_late": 8.5})
    gl = [geoms[k] for k in heart]
    add_rows("common_step", [sp.with_common_step(g, gl[-1].z_step) for g in gl])
    add_rows("iso_z", [iso_geometry(geoms[k], s_iso[k]) for k in heart])
    add_rows("raw_with_E8.0", [geoms[k] for k in ["E8.0"] + heart])
    # E8.75 from (E8.0, E8.25_late): the remodelling window, for context
    add_rows("raw_context", [geoms[k] for k in ["E8.0", "E8.25_late", "E8.75"]], target="E8.75")
    hp = pd.DataFrame(rows)
    hp.to_csv(out / "size_heart_pseudo.csv", index=False)

    # every z mode through the procedure itself: E9.5 from (E8.25_late, E8.75), all rules and scope hypotheses
    mode_rows = []
    for mode in ("raw", "common_step", "xy_grow_z_hold", "split"):
        res = sp.heart_size_from_geometry([geoms[k] for k in heart[:2]], 9.5, z_mode=mode)
        ref = geoms["E8.75"]
        tag = "core" if mode in ("raw", "common_step") else "core_xy"
        for r, y in res[f"all_rule_predictions_{tag}"].items():
            if mode == "xy_grow_z_hold":
                core_r = float(np.hypot(np.exp(y), ref.sd_z_core))
            elif mode == "split":
                yz = res["all_rule_predictions_core_z"][r]
                core_r = float(np.hypot(np.exp(y), np.exp(yz)))
            else:
                core_r = float(np.exp(y))
            mode_rows.append({"z_mode": mode, "rule": r, "pred_rms_core": core_r,
                              "core_error": float(np.log(core_r / geoms["E9.5"].rms_core))})
        mode_rows.append({"z_mode": mode, "rule": "procedure(avg of " + "+".join(res[f"chosen_{tag}"]) + ")",
                          "pred_rms_core": res["pred_rms_core"],
                          "core_error": float(np.log(res["pred_rms_core"] / geoms["E9.5"].rms_core))}
                         | {f"full_{h}": v["pred_rms_full"] for h, v in res["scope"].items()}
                         | {f"full_error_{h}": float(np.log(v["pred_rms_full"] / geoms["E9.5"].rms_full))
                            for h, v in res["scope"].items()}
                         | {f"full_skill_{h}": scale_skill(v["pred_rms_full"], stages["E9.5"], stages["E8.75"])[
                             "skill_mean"] for h, v in res["scope"].items()})
    pd.DataFrame(mode_rows).to_csv(out / "size_heart_modes.csv", index=False)

    # ---------------- heart: scope hypotheses on E9.5 ----------------
    scope_rows = []
    for h in ("H1", "H0", "H2"):
        s = sp.scope_log_ratio([geoms["E8.25_late"].log_scope, geoms["E8.75"].log_scope], [8.25, 8.75], 9.5, h)
        scope_rows.append({"hypothesis": h, "pred_log_scope": s, "obs_log_scope": geoms["E9.5"].log_scope,
                           "error": s - geoms["E9.5"].log_scope})
    pd.DataFrame(scope_rows).to_csv(out / "size_heart_scope.csv", index=False)

    # ---------------- heart: full-cloud pseudo-target with scale skill (floor = copy of E8.75) ----------------
    full_rows = []
    core_pred = hp[(hp.variant == "raw") & (hp.series == "core")].set_index("rule")["pred_rms"].to_dict()
    for r in rules:
        for sr in scope_rows:
            pr = core_pred[r] * float(np.exp(sr["pred_log_scope"]))
            sk = scale_skill(pr, stages["E9.5"], stages["E8.75"])
            full_rows.append({"core_rule": r, "scope": sr["hypothesis"], "pred_rms_full": pr,
                              "obs_rms_full": geoms["E9.5"].rms_full,
                              "error": float(np.log(pr / geoms["E9.5"].rms_full))} | sk)
    for label, pr in (("direct_full_latest_rate", float(np.exp(hp[(hp.variant == "raw") & (hp.series == "full") &
                                                                   (hp.rule == "latest_rate")].pred_log.iloc[0]))),
                      ):
        sk = scale_skill(pr, stages["E9.5"], stages["E8.75"])
        full_rows.append({"core_rule": label, "scope": "-", "pred_rms_full": pr, "obs_rms_full": geoms["E9.5"].rms_full,
                          "error": float(np.log(pr / geoms["E9.5"].rms_full))} | sk)
    pd.DataFrame(full_rows).to_csv(out / "size_heart_full.csv", index=False)

    # ---------------- embryo: interpolation checks ----------------
    emb = ["E6.75", "E7.25", "E8.0"]
    emb_rows = []
    for variant in ("raw", "common_step"):
        gl = [geoms[k] for k in emb]
        if variant == "common_step":
            gl = [sp.with_common_step(g, gl[-1].z_step) for g in gl]
        res = sp.embryo_size_from_geometry(gl, 7.75, z_mode="raw")
        for lo in res["loo_interior"]:
            pr = float(np.exp(lo["pred"]))
            # the board's floor for an interpolation target is a copy of the previous stage
            sk = scale_skill(pr, stages[lo["left_out"]], stages["E6.75"]) if variant == "raw" else {}
            emb_rows.append({"variant": variant, "check": "loo_interior", "target": lo["left_out"],
                             "pred_rms": pr, "obs_rms": float(np.exp(lo["obs"])), "error": lo["error"]} | sk)
        # extrapolative checks for context (latest rate)
        ts = np.array([g.t for g in gl])
        ys = np.log([g.rms_full for g in gl])
        p = sp.rule_latest_rate(ts[:2], ys[:2], ts[2])
        emb_rows.append({"variant": variant, "check": "extrap_latest_rate", "target": gl[2].name,
                         "pred_rms": float(np.exp(p)), "obs_rms": float(np.exp(ys[2])), "error": p - ys[2]})
        # also xy-only and z-only views of the same interpolation (raw only)
        if variant == "raw":
            for key, ser in (("xy_full", [g.rms_xy_full for g in gl]), ("sdz_full", [g.sd_z_full for g in gl])):
                y = np.log(ser)
                p = y[0] + (y[2] - y[0]) * (ts[1] - ts[0]) / (ts[2] - ts[0])
                emb_rows.append({"variant": key, "check": "loo_interior", "target": "E7.25",
                                 "pred_rms": float(np.exp(p)), "obs_rms": float(np.exp(y[1])), "error": p - y[1]})
    pd.DataFrame(emb_rows).to_csv(out / "size_embryo.csv", index=False)

    # ---------------- embryo-to-embryo spread and rule noise multipliers ----------------
    wt, ko = geoms["E9.5"], geoms["E9.5_mab21l2_ko"]
    noise = []
    for key in ("rms_core", "rms_full", "rms_xy_core", "rms_xy_full", "sd_z_core", "sd_z_full"):
        d = float(np.log(getattr(ko, key) / getattr(wt, key)))
        noise.append({"quantity": key, "log_ko_over_wt": d, "sigma_e_if_independent": abs(d) / np.sqrt(2)})
    # linear-rule noise multipliers for the test layout (reference E10.5, previous E9.5, target E12.5): the SD of the
    # prediction error in units of the per-embryo SD sigma_e, counting the reference, previous and target embryos
    t1, t2, t3 = 9.5, 10.5, 12.5
    k_latest = (t3 - t2) / (t2 - t1)
    k_power = np.log(t3 / t2) / np.log(t2 / t1)
    for name, k in (("copy_last", 0.0), ("half_rate", 0.5 * k_latest), ("power_t", k_power),
                    ("latest_rate", k_latest)):
        mult = float(np.sqrt((1 + k) ** 2 + k ** 2 + 1))
        noise.append({"quantity": f"multiplier_{name}", "log_ko_over_wt": np.nan, "sigma_e_if_independent": mult})
    # decaying rate: Monte Carlo around a decelerating baseline built from rates in the released data scale
    pd.DataFrame(noise).to_csv(out / "size_noise.csv", index=False)

    # ---------------- count-density diagnostic (not part of the procedure) ----------------
    cd = []
    hs = ["E8.25_late", "E8.75", "E9.5", "E9.5_mab21l2_ko"]
    c = {k: geoms[k].rms_full / geoms[k].n_cells ** (1 / 3) for k in FILES}
    for k in hs:
        others = [c[o] for o in hs if o != k]
        pred = float(np.median(others)) * geoms[k].n_cells ** (1 / 3)
        cd.append({"stage": k, "c_rms_over_cuberoot_n": c[k], "loo_pred_rms": pred, "obs_rms": geoms[k].rms_full,
                   "error": float(np.log(pred / geoms[k].rms_full))})
    for k in ("E6.75", "E7.25", "E8.0"):
        cd.append({"stage": k, "c_rms_over_cuberoot_n": c[k], "loo_pred_rms": np.nan, "obs_rms": geoms[k].rms_full,
                   "error": np.nan})
    pd.DataFrame(cd).to_csv(out / "size_count_density.csv", index=False)

    # ---------------- ceiling noise of the metric ----------------
    ce = []
    for k, st in stages.items():
        v = []
        for s in SEEDS:
            truth, other = split_half(subsample(st, 0.1, s), seed=0)
            v.append(scale_log_ratio(other.coords, truth.coords))
        ce.append({"stage": k, "truth_cells": int(round(st.n * 0.1)) // 2,
                   "ceiling_abs_mean": float(np.mean(np.abs(v))), "ceiling_abs_max": float(np.max(np.abs(v)))})
    pd.DataFrame(ce).to_csv(out / "size_ceiling.csv", index=False)

    # ---------------- dry runs of the procedure on the stages released so far ----------------
    dh = sp.heart_size_from_geometry([geoms[k] for k in heart], 12.5)
    de = sp.embryo_size_from_geometry([geoms[k] for k in emb], 7.75)
    (out / "size_dryrun_heart.json").write_text(json.dumps(dh, indent=2, default=float))
    (out / "size_dryrun_embryo.json").write_text(json.dumps(de, indent=2, default=float))

    pd.set_option("display.width", 220)
    print(hp[(hp.series == "core")].pivot(index="rule", columns="variant", values="error").round(3))
    print(pd.DataFrame(mode_rows).round(3).to_string(index=False))
    print(pd.DataFrame(scope_rows).round(3))
    print(pd.DataFrame(full_rows)[["core_rule", "scope", "pred_rms_full", "error", "skill_mean", "floor_mean",
                                   "ceiling_abs_mean"]].round(3).to_string(index=False))
    print(pd.DataFrame(emb_rows).round(3).to_string(index=False))
    print(pd.DataFrame(noise).round(3).to_string(index=False))
    print(pd.DataFrame(cd).round(3).to_string(index=False))
    print(pd.DataFrame(ce).round(4).to_string(index=False))
    print("dry run heart:", {k: v for k, v in dh.items() if k in ("z_mode", "pred_rms_core", "pred_log_growth_core")
                             or k.startswith("chosen_")},
          {h: round(v["pred_rms_full"], 1) for h, v in dh["scope"].items()})
    print("dry run embryo:", de["pred_rms_full"], de["bracket"], de["loo_interior"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    main(Path(ap.parse_args().out))
