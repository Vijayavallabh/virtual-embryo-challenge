"""Part 4: evaluate candidate size rules on pseudo-targets built from released training stages only (DGX).

Geometry comes from size_procedure.measure_stage (duplicated records dropped, unmapped labels refused), exactly as the
procedure measures it. Scorer emulation (scale skill, ceilings) uses the files as released, because the scorer draws
its truth halves from them. No rate value is assumed anywhere: the rule comparison is size_procedure.rule_table.

Writes to --out:
  stage_geometry.csv         per-file geometry used by the size procedure, with the duplicates dropped per file
  size_heart_pseudo.csv      heart rules on the E9.5 pseudo-target (and the E8.75 context target), with variants
  size_heart_modes.csv       the procedure on the E9.5 pseudo-target for every z mode; every row carries the run's valid
                             flag: pre-registered settings (valid=False, no rate interval starts at or after E9.0, so
                             the numbers are diagnostic only) and the same rules with the cut-off removed
  size_heart_scope.csv       scope log-ratio hypotheses on the E9.5 pseudo-target
  size_heart_full.csv        full-cloud pseudo-target errors (core rule x scope hypothesis) and scale skill
  size_rule_table.csv        the rate-free rule table for the 2026-10-20 layout (E9.5, E10.5 -> E12.5)
  size_embryo.csv            embryo leave-one-out checks for every z mode, scale skill, and an extrapolation check
  size_noise.csv             informational: E9.5 WT vs Mab21l2 KO spread (the KO is not used by the size procedure)
  size_count_density.csv     informational: RMS / N^(1/3) (needs an organiser ruling before any use)
  size_dedup_sensitivity.csv informational: RMS as released, de-duplicated (what the procedure uses), and with copies
                             whose xy moved and interpolated (non-integer) records also removed
  size_ceiling.csv           split-half |scale_log_ratio| (the metric's ceiling) per released stage, 10 seeds
  size_dryrun_heart.json, size_dryrun_embryo.json   dry runs on the current stages with every target forecast
                             redacted (code checks, not predictions); the unredacted JSONs go to --out/dryrun_full/ and
                             stay on the DGX
Usage: uv run python experiments/coords/evaluate_size.py [--out runs/coords] [--timestamp TEXT]
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import size_procedure as sp  # noqa: E402
from coordlib import FILES, nominal_t, raw_coords  # noqa: E402

from vec.evaluate import Stage, _skill, split_half, subsample  # noqa: E402

T = {k: nominal_t(k) for k in FILES}
SEEDS = range(10)
HEART = ["E8.25_late", "E8.75", "E9.5"]
EMBRYO = ["E6.75", "E7.25", "E8.0"]
INFO = "informational: not used by size_procedure"


def scale_log_ratio(pred_C, true_C):
    from common.shape_metrics import scale_log_ratio as slr

    return slr(pred_C, true_C)


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


def dedup_sensitivity(name: str, path: str) -> dict:
    """RMS of the full cloud and the core under four record sets (informational)."""
    import anndata as ad
    from scipy import sparse

    from vec.prep import add_lineage
    from vec.spatial import find_duplicate_records

    a = ad.read_h5ad(path)
    add_lineage(a, source=Path(path).name.removesuffix(".h5ad"))
    C = np.asarray(a.obsm["spatial_3D"], float)[:, :3]
    core = np.asarray(a.obs["core"], bool)
    m_dup, _ = find_duplicate_records(a)
    m_any, _ = find_duplicate_records(a, xy_tol=np.inf)
    M = sparse.csr_matrix(a.layers["counts"], dtype=np.float64)
    frac = np.abs(M.data - np.rint(M.data)) > 1e-6
    nonint = np.zeros(a.n_obs, bool)
    nonint[np.unique(np.repeat(np.arange(a.n_obs), np.diff(M.indptr))[frac])] = True
    row = {"stage": name, "n_as_released": int(a.n_obs), "n_duplicates": int(m_dup.sum()),
           "n_duplicates_any_xy": int(m_any.sum()), "n_non_integer_rows": int(nonint.sum()), "note": INFO}
    for tag, drop in (("as_released", np.zeros(a.n_obs, bool)), ("dedup", m_dup), ("dedup_any_xy", m_any),
                      ("dedup_any_xy_no_nonint", m_any | nonint)):
        keep = ~drop
        row[f"rms_full_{tag}"] = sp.rms_radius(C[keep])
        row[f"rms_core_{tag}"] = sp.rms_radius(C[keep & core]) if (keep & core).sum() >= 4 else np.nan
    for tag in ("as_released", "dedup_any_xy", "dedup_any_xy_no_nonint"):
        for part in ("full", "core"):
            row[f"log_{part}_{tag}_vs_dedup"] = float(np.log(row[f"rms_{part}_{tag}"] / row[f"rms_{part}_dedup"]))
    return row


def main(out: Path, timestamp: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    geoms, infos, stages = {}, {}, {}
    for name, path in FILES.items():
        geoms[name], infos[name] = sp.measure_stage(path, T[name], name=name)
        C = raw_coords(path)
        stages[name] = Stage(X=np.zeros((C.shape[0], 1), np.float32), coords=C)
        print(f"{name}: {infos[name]['n_input']} records, dropped {infos[name]['n_dropped_duplicates']} duplicates; "
              f"core {geoms[name].n_core}", flush=True)

    geo_rows = []
    for name, g in geoms.items():
        geo_rows.append(asdict(g) | {"log_scope": g.log_scope if np.isfinite(g.rms_core) else np.nan,
                                     "log_scope_xy": g.log_scope_xy if np.isfinite(g.rms_core) else np.nan,
                                     "core_thin_axis_angle_deg": g.core_thin_axis_angle_deg,
                                     "n_input": infos[name]["n_input"],
                                     "n_dropped_duplicates": infos[name]["n_dropped_duplicates"],
                                     "rms_full_as_released": sp.rms_radius(stages[name].coords),
                                     "setting_note": INFO if name == "E9.5_mab21l2_ko" else ""})
    pd.DataFrame(geo_rows).to_csv(out / "stage_geometry.csv", index=False)

    # ---------------- heart: rules on E9.5 from (E8.25_late, E8.75) ----------------
    rules = list(sp.RULES)
    rows = []

    def add_rows(variant, gl, target="E9.5", t_override=None):
        ts = np.array([t_override.get(g.name, g.t) if t_override else g.t for g in gl])
        for key, ys in (("core", np.log([g.rms_core for g in gl])), ("full", np.log([g.rms_full for g in gl])),
                        ("xy_core", np.log([g.rms_xy_core for g in gl]))):
            for r in rules:
                pred = sp.RULES[r](ts[:-1], ys[:-1], ts[-1])
                rows.append({"variant": variant, "series": key, "target": target, "rule": r,
                             "preregistered_rule": r in sp.PREREGISTERED["heart_rules"],
                             "train": ",".join(g.name for g in gl[:-1]), "pred_log": pred, "obs_log": ys[-1],
                             "error": pred - ys[-1], "pred_rms": float(np.exp(pred)), "obs_rms": float(np.exp(ys[-1]))})

    gl = [geoms[k] for k in HEART]
    add_rows("raw", gl)
    add_rows("raw_t8.25late=8.375", gl, t_override={"E8.25_late": 8.375})
    add_rows("raw_t8.25late=8.5", gl, t_override={"E8.25_late": 8.5})
    add_rows("common_step_ref", [sp.with_common_step(g, gl[-2].z_step) for g in gl])  # the procedure's choice
    add_rows("common_step_target", [sp.with_common_step(g, gl[-1].z_step) for g in gl])  # oracle: the target's step
    add_rows("raw_context", [geoms[k] for k in ["E8.0", "E8.25_late", "E8.75"]], target="E8.75")
    hp = pd.DataFrame(rows)
    hp.to_csv(out / "size_heart_pseudo.csv", index=False)

    # ---------------- every z mode through the procedure itself, with the valid flag on every row ----------------
    mode_rows = []
    obs_core, obs_full = geoms["E9.5"].rms_core, geoms["E9.5"].rms_full
    for mode in sp.Z_MODES:
        for settings, mst in (("preregistered", None), ("no_rate_cutoff_diagnostic", -np.inf)):
            res = sp.heart_size_from_geometry([geoms[k] for k in HEART[:2]], 9.5, z_mode=mode, min_start_t=mst)
            ref = sorted([geoms[k] for k in HEART[:2]], key=lambda g: g.t)[-1]
            tag = "core" if mode in ("raw", "common_step") else "core_xy"
            base = {"z_mode": mode, "settings": settings, "min_start_t": res["min_start_t"], "valid": res["valid"]}
            for r, y in res["growth"][tag]["all_rule_predictions"].items():
                if mode == "xy_grow_z_hold":
                    core_r = float(np.hypot(np.exp(y), ref.sd_z_core))
                elif mode == "split":
                    core_r = float(np.hypot(np.exp(y), np.exp(res["growth"]["core_z"]["all_rule_predictions"][r])))
                else:
                    core_r = float(np.exp(y))
                mode_rows.append(base | {"rule": r, "pred_rms_core": core_r,
                                         "core_error": float(np.log(core_r / obs_core))})
            chosen = "+".join(res["growth"][tag]["chosen"])
            row = base | {"rule": f"procedure: average of {chosen}", "pred_rms_core": res["pred_rms_core"],
                          "core_error": float(np.log(res["pred_rms_core"] / obs_core))}
            for h, v in res["scope"].items():
                row[f"full_{h}"] = v["pred_rms_full"]
                row[f"full_error_{h}"] = float(np.log(v["pred_rms_full"] / obs_full))
                row[f"full_skill_{h}"] = scale_skill(v["pred_rms_full"], stages["E9.5"], stages["E8.75"])["skill_mean"]
            mode_rows.append(row)
    pd.DataFrame(mode_rows).to_csv(out / "size_heart_modes.csv", index=False)

    # ---------------- heart: scope hypotheses on E9.5 ----------------
    scope_rows = []
    for h in sp.PREREGISTERED["heart_scope_hypotheses"]:
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
            full_rows.append({"core_rule": r, "scope": sr["hypothesis"], "pred_rms_full": pr, "obs_rms_full": obs_full,
                              "error": float(np.log(pr / obs_full))} | sk)
    pr = float(np.exp(hp[(hp.variant == "raw") & (hp.series == "full") & (hp.rule == "latest_rate")].pred_log.iloc[0]))
    full_rows.append({"core_rule": "direct_full_latest_rate", "scope": "-", "pred_rms_full": pr,
                      "obs_rms_full": obs_full, "error": float(np.log(pr / obs_full))}
                     | scale_skill(pr, stages["E9.5"], stages["E8.75"]))
    pd.DataFrame(full_rows).to_csv(out / "size_heart_full.csv", index=False)

    # ---------------- the rate-free rule table for the release layout ----------------
    table = pd.DataFrame(sp.rule_table(9.5, 10.5, 12.5))
    table.to_csv(out / "size_rule_table.csv", index=False)

    # ---------------- embryo: leave-one-out for every z mode ----------------
    emb_rows = []
    rep = sp.embryo_report([geoms[k] for k in EMBRYO], 7.75)
    for mode, res in rep["z_modes"].items():
        for lo in res["loo_interior"]:
            pr = float(np.exp(lo["pred"]))
            sk = scale_skill(pr, stages[lo["left_out"]], stages["E6.75"]) if mode == "raw" else {}
            emb_rows.append({"z_mode": mode, "preregistered": mode == rep["preregistered"]["z_mode"],
                             "check": "loo_interior", "target": lo["left_out"], "pred_rms": pr,
                             "obs_rms": float(np.exp(lo["obs"])), "error": lo["error"]} | sk)
    ts = np.array([geoms[k].t for k in EMBRYO])
    ys = np.log([geoms[k].rms_full for k in EMBRYO])
    p = sp.rule_latest_rate(ts[:2], ys[:2], ts[2])
    emb_rows.append({"z_mode": "raw", "preregistered": False, "check": "extrap_latest_rate", "target": EMBRYO[2],
                     "pred_rms": float(np.exp(p)), "obs_rms": float(np.exp(ys[2])), "error": p - ys[2]})
    pd.DataFrame(emb_rows).to_csv(out / "size_embryo.csv", index=False)

    # ---------------- informational: embryo-to-embryo spread (the KO is not used by the procedure) ----------------
    wt, ko = geoms["E9.5"], geoms["E9.5_mab21l2_ko"]
    noise = []
    for key in ("rms_core", "rms_full", "rms_xy_core", "rms_xy_full", "sd_z_core", "sd_z_full"):
        d = float(np.log(getattr(ko, key) / getattr(wt, key)))
        noise.append({"quantity": key, "log_ko_over_wt": d, "sigma_e_if_independent": abs(d) / np.sqrt(2),
                      "note": INFO})
    pd.DataFrame(noise).to_csv(out / "size_noise.csv", index=False)

    # ---------------- informational: count-density diagnostic ----------------
    cd = []
    hs = ["E8.25_late", "E8.75", "E9.5", "E9.5_mab21l2_ko"]
    c = {k: geoms[k].rms_full / geoms[k].n_cells ** (1 / 3) for k in FILES}
    for k in hs:
        pred = float(np.median([c[o] for o in hs if o != k])) * geoms[k].n_cells ** (1 / 3)
        cd.append({"stage": k, "c_rms_over_cuberoot_n": c[k], "loo_pred_rms": pred, "obs_rms": geoms[k].rms_full,
                   "error": float(np.log(pred / geoms[k].rms_full)), "note": INFO})
    for k in EMBRYO:
        cd.append({"stage": k, "c_rms_over_cuberoot_n": c[k], "loo_pred_rms": np.nan, "obs_rms": geoms[k].rms_full,
                   "error": np.nan, "note": INFO})
    pd.DataFrame(cd).to_csv(out / "size_count_density.csv", index=False)

    # ---------------- informational: how much the record set moves the RMS ----------------
    ds = pd.DataFrame([dedup_sensitivity(k, p) for k, p in FILES.items()])
    ds.to_csv(out / "size_dedup_sensitivity.csv", index=False)

    # ---------------- ceiling noise of the metric (files as released) ----------------
    ce = []
    for k, st in stages.items():
        v = []
        for s in SEEDS:
            truth, other = split_half(subsample(st, 0.1, s), seed=0)
            v.append(scale_log_ratio(other.coords, truth.coords))
        ce.append({"stage": k, "truth_cells": int(round(st.n * 0.1)) // 2,
                   "ceiling_abs_mean": float(np.mean(np.abs(v))), "ceiling_abs_max": float(np.max(np.abs(v)))})
    pd.DataFrame(ce).to_csv(out / "size_ceiling.csv", index=False)

    # ---------------- dry runs on the stages released so far (code checks; forecasts redacted) ----------------
    (out / "dryrun_full").mkdir(exist_ok=True)
    dry = {}
    for board, names, target in (("heart", HEART, 12.5), ("embryo", EMBRYO, 7.75)):
        res = sp.run(board, [(FILES[k], T[k]) for k in names], target, timestamp)
        (out / "dryrun_full" / f"size_dryrun_{board}.json").write_text(sp.dumps(res))
        dry[board] = sp.redact_forecasts(res)
        (out / f"size_dryrun_{board}.json").write_text(sp.dumps(dry[board]))

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(hp[(hp.series == "core")].pivot(index="rule", columns="variant", values="error").round(3))
    print(pd.DataFrame(mode_rows).round(3).to_string(index=False))
    print(pd.DataFrame(scope_rows).round(3))
    print(pd.DataFrame(full_rows)[["core_rule", "scope", "pred_rms_full", "error", "skill_mean", "floor_mean",
                                   "ceiling_abs_mean"]].round(3).to_string(index=False))
    print(table.round(3).to_string(index=False))
    print(pd.DataFrame(emb_rows).round(3).to_string(index=False))
    print(pd.DataFrame(noise).round(3).to_string(index=False))
    print(pd.DataFrame(cd).round(3).to_string(index=False))
    print(ds.round(4).to_string(index=False))
    print(pd.DataFrame(ce).round(4).to_string(index=False))
    for board, d in dry.items():
        print(f"dry run {board} (forecasts redacted): valid={d['valid']} preregistered z mode "
              f"{d['preregistered']['z_mode']}; per mode:",
              [(s["z_mode"], s["valid"], s.get("chosen")) for s in d["summary"]])
        if board == "embryo":
            print("  leave-one-out:", [(s["z_mode"], [(x["left_out"], round(x["error"], 3)) for x in s["loo_interior"]])
                                       for s in d["summary"]])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    ap.add_argument("--timestamp", default="unset", help="recorded in the dry-run JSONs")
    args = ap.parse_args()
    main(Path(args.out), args.timestamp)
