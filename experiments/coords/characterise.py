"""Part 1: characterise obsm["spatial_3D"] in every MERFISH file (run on the DGX).

Duplicated records (vec.spatial.find_duplicate_records: identical count rows at the same xy in the same or an adjacent
level) are dropped before anything is measured (vec.prep.load_stage(..., drop_duplicates=True)); the numbers dropped
are recorded. Per-section counts are normalised by the largest physical section: a level holding two sections (the
Mab21l2 KO stores sample_ids 30 and 35 at one z) counts as two when the file records obs["sample_id"]. Vertical
neighbours are found in the level one step above by z difference, not by index.

Writes to --out (geometry and counts only):
  file_summary.csv     one row per file: duplicates dropped, extents, z levels and steps, spacing, density, RMS,
                       principal SDs, core subset, end-section fractions
  z_levels.csv         one row per (file, z level): cells, core cells, dropped duplicates, sections, spacing, density
  z_raw_values.csv     every raw unique z value: its snapped level, deviation, cell count and xy footprint
Usage: uv run python experiments/coords/characterise.py [--out runs/coords] [--stage NAME=PATH ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coordlib import (  # noqa: E402
    FILES,
    SETTING,
    coords_of,
    dropped_per_level,
    inplane_nn,
    level_partners,
    level_step,
    load,
    local_density_2d,
    multi_section_levels,
    occupied_area,
    principal_sds,
    rms_radius,
    rms_xy,
    section_ids,
    snap_z,
    stage_specs,
    summarize,
    vertical_offset,
)


def shape_block(C: np.ndarray, prefix: str) -> dict:
    out = {}
    if C.shape[0] < 4:
        return out
    sds, V = principal_sds(C)
    out[f"{prefix}rms_3d"] = rms_radius(C)
    out[f"{prefix}rms_xy"] = rms_xy(C)
    out[f"{prefix}sd_z"] = float(C[:, 2].std())
    out[f"{prefix}frac_var_z"] = float(C[:, 2].var() / max(rms_radius(C) ** 2, 1e-12))
    for i in range(3):
        out[f"{prefix}pc{i + 1}_sd"] = float(sds[i])
        out[f"{prefix}pc{i + 1}_cos_z"] = float(abs(V[2, i]))
    out[f"{prefix}sd_ratio_12"] = float(sds[0] / sds[1])
    out[f"{prefix}sd_ratio_23"] = float(sds[1] / sds[2])
    out[f"{prefix}sd_ratio_13"] = float(sds[0] / sds[2])
    out[f"{prefix}thin_axis_angle_z_deg"] = float(np.degrees(np.arccos(min(abs(V[2, 2]), 1.0))))
    return out


def end_fractions(per_section: pd.Series, per_level: np.ndarray) -> tuple[float, float, float]:
    """(peak cells in one physical section, first-level count / peak, last-level count / peak)."""
    peak = float(per_section.max()) if len(per_section) else 0.0
    if peak <= 0:
        return peak, np.nan, np.nan
    return peak, float(per_level[0] / peak), float(per_level[-1] / peak)


def main(out_dir: Path, files: dict[str, str]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows, level_rows, raw_rows = [], [], []
    for stage, path in files.items():
        a, dd = load(path, drop_duplicates=True)
        C = coords_of(a)
        core = np.asarray(a.obs["core"], dtype=bool)
        z = C[:, 2]
        lev, levels, spread = snap_z(z)
        steps = np.diff(levels)
        med_step = level_step(levels)
        vals, cnt = np.unique(np.round(steps, 3), return_counts=True)
        step_str = "; ".join(f"{v:g}x{c}" for v, c in zip(vals, cnt, strict=True))
        p1 = level_partners(levels, 1, med_step)
        multi = multi_section_levels(a, lev)
        sec = section_ids(a, lev)

        nn_ip = inplane_nn(C, lev, k=1)
        dens = local_density_2d(C, lev, k=10)
        d3, i3 = cKDTree(C).query(C, k=2)
        same_level = lev[i3[:, 1]] == lev
        vo = vertical_offset(C, lev, levels, dk=1, step=med_step)
        dup = int(C.shape[0] - np.unique(C, axis=0).shape[0])
        per_level_n = np.bincount(lev, minlength=levels.size)
        per_level_core = np.bincount(lev[core], minlength=levels.size)
        per_section = pd.Series(1, index=sec).groupby(level=0).sum()
        per_section_core = pd.Series(1, index=sec[core]).groupby(level=0).sum()
        drop_lev = dropped_per_level(a, levels)

        # per-level occupied area and density on two pixel sizes
        area20 = np.array([occupied_area(C[lev == g, :2], 20.0) for g in range(levels.size)])
        area40 = np.array([occupied_area(C[lev == g, :2], 40.0) for g in range(levels.size)])
        ok = per_level_n >= 50
        rho20 = per_level_n[ok] / area20[ok]
        rho40 = per_level_n[ok] / area40[ok]

        rho_loc = summarize(dens)["median"]
        peak, first_f, last_f = end_fractions(per_section, per_level_n)
        peak_c, first_c, last_c = end_fractions(per_section_core, per_level_core)
        row = {
            "stage": stage, "setting": SETTING.get(stage, "unknown"), "n_input": dd["n_input"],
            "n_dropped_duplicates": dd["n_dropped"],
            "dropped_by_gap": ";".join(f"gap{k}:{v}" for k, v in dd["dropped_by_gap"].items()),
            "n_cells": int(C.shape[0]), "n_core": int(core.sum()),
            "core_frac": float(core.mean()), "shell_frac": float(1 - core.mean()),
            "n_duplicate_coords": dup,
        }
        for j, ax in enumerate("xyz"):
            row[f"{ax}_min"], row[f"{ax}_max"] = float(C[:, j].min()), float(C[:, j].max())
            row[f"{ax}_range"] = float(np.ptp(C[:, j]))
            row[f"{ax}_mean"] = float(C[:, j].mean())
        row.update({
            "z_unique_raw": int(np.unique(z).size), "z_levels": int(levels.size),
            "z_within_level_spread": spread, "z_step_median": med_step, "z_step_min": float(steps.min()),
            "z_step_max": float(steps.max()), "z_steps": step_str,
            "z_gaps": int((steps > 1.5 * med_step).sum()),
            "z_levels_without_partner_above": int((p1[:-1] < 0).sum()),
            "multi_section_levels": ";".join(f"{levels[g]:g}" for g in multi),
            "n_sections": int(per_section.size),
            "cells_per_level_median": float(np.median(per_level_n)), "cells_per_level_max": int(per_level_n.max()),
            "peak_section_cells": peak, "first_section_frac": first_f, "last_section_frac": last_f,
            "peak_section_core": peak_c, "first_section_core_frac": first_c, "last_section_core_frac": last_c,
            "inplane_nn_median": summarize(nn_ip)["median"], "inplane_nn_p10": summarize(nn_ip)["p10"],
            "inplane_nn_p90": summarize(nn_ip)["p90"],
            "nn3d_median": float(np.median(d3[:, 1])), "frac_nn3d_same_level": float(same_level.mean()),
            "vertical_offset_median": summarize(vo)["median"],
            "frac_vertical_offset_lt2": float(np.nanmean(vo < 2.0)),
            "frac_vertical_offset_lt4": float(np.nanmean(vo < 4.0)),
            "rho2d_local_per_1000": rho_loc * 1000.0,
            "rho2d_grid20_per_1000": float(np.median(rho20)) * 1000.0,
            "rho2d_grid40_per_1000": float(np.median(rho40)) * 1000.0,
            "step_over_inplane_nn": med_step / summarize(nn_ip)["median"],
            "sd_z_over_uniform": float(z.std() / (np.ptp(z) / np.sqrt(12.0))),
        })
        row.update(shape_block(C, ""))
        if core.sum() >= 4:
            row.update(shape_block(C[core], "core_"))
            row["rms_full_over_core"] = row["rms_3d"] / row["core_rms_3d"]
            row["log_scope"] = float(np.log(row["rms_full_over_core"]))
            row["core_centroid_offset"] = float(np.linalg.norm(C[core].mean(0) - C.mean(0)))
            row["core_levels_ge5"] = int((per_level_core >= 5).sum())
            row["core_inplane_nn_median"] = summarize(inplane_nn(C[core], lev[core], k=1))["median"]
        rows.append(row)

        n_sec_level = pd.Series(sec).groupby(lev).nunique().reindex(range(levels.size), fill_value=0).to_numpy()
        for g in range(levels.size):
            m = lev == g
            level_rows.append({
                "stage": stage, "level": g, "z": float(levels[g]), "n_cells": int(m.sum()),
                "n_core": int((m & core).sum()), "n_dropped_duplicates": int(drop_lev[g]),
                "n_sections": int(n_sec_level[g]), "partner_above": int(p1[g]),
                "inplane_nn_median": summarize(nn_ip[m])["median"],
                "rho2d_local_median_per_1000": summarize(dens[m])["median"] * 1000.0,
                "area_px20": float(area20[g]), "area_px40": float(area40[g]),
                "x_range": float(np.ptp(C[m, 0])) if m.any() else np.nan,
                "y_range": float(np.ptp(C[m, 1])) if m.any() else np.nan,
                "n_raw_z_values": int(np.unique(z[m]).size),
            })

        for v in np.unique(z):
            m = z == v
            g = int(lev[np.flatnonzero(m)[0]])
            raw_rows.append({
                "stage": stage, "z_raw": repr(float(v)), "level": g, "level_z": float(levels[g]),
                "deviation": float(v - levels[g]), "n_cells": int(m.sum()),
                "x_centroid": float(C[m, 0].mean()), "y_centroid": float(C[m, 1].mean()),
                "x_range": float(np.ptp(C[m, 0])), "y_range": float(np.ptp(C[m, 1])),
            })
        print(f"{stage}: n={C.shape[0]} (dropped {dd['n_dropped']} of {dd['n_input']}) levels={levels.size} "
              f"(raw unique {np.unique(z).size}) step {step_str} spread {spread:.2e} sections {per_section.size} "
              f"first/last {first_f:.3f}/{last_f:.3f} core {first_c:.3f}/{last_c:.3f} "
              f"nn_ip {row['inplane_nn_median']:.2f} rms3d {row['rms_3d']:.1f}", flush=True)

    pd.DataFrame(rows).to_csv(out_dir / "file_summary.csv", index=False)
    pd.DataFrame(level_rows).to_csv(out_dir / "z_levels.csv", index=False)
    pd.DataFrame(raw_rows).to_csv(out_dir / "z_raw_values.csv", index=False)
    print("wrote", out_dir)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    ap.add_argument("--stage", action="append", help="NAME=PATH (repeatable; default: every released file)")
    args = ap.parse_args()
    main(Path(args.out), stage_specs(args.stage, FILES))
