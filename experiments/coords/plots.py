"""Figures for REPORT.md, drawn from the CSVs in runs/coords (DGX). Writes three PNGs next to the CSVs.

Revision 2026-10-19: counts are after dropping duplicated records, and each panel of the z profiles is normalised by
the largest physical section (the Mab21l2 KO's z = -320 level holds two sections and no longer sets the scale); the
pitch figure uses the main k = 1-3 fit.

Palette: reference slots blue / orange / aqua (validated all-pairs, light mode) plus muted gray for context; aqua is
below 3:1 on the surface, so every figure carries a legend and distinct marker shapes.
Usage: uv run python experiments/coords/plots.py [--dir runs/coords]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BLUE, ORANGE, AQUA, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
SURFACE, INK, INK2, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
ORDER = ["E6.75", "E7.25", "E8.0", "E8.25_late", "E8.75", "E9.5", "E9.5_mab21l2_ko"]
LABEL = {"E9.5_mab21l2_ko": "E9.5 Mab21l2 KO"}


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.8, ls="-")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.title.set_color(INK)


def plain_log_ticks(ax, ticks):
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator

    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_yticklabels([str(t) for t in ticks])


def fig_z_profiles(d: Path):
    L = pd.read_csv(d / "z_levels.csv")
    F = pd.read_csv(d / "file_summary.csv").set_index("stage")
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.2), facecolor=SURFACE)
    for ax, st in zip(axes.flat, ORDER, strict=False):
        g = L[L.stage == st]
        f = F.loc[st]
        ax.plot(g.level, g.n_cells / f.peak_section_cells, color=BLUE, lw=1.8, solid_capstyle="round",
                label="all cells")
        if g.n_core.max() > 0:
            ax.plot(g.level, g.n_core / f.peak_section_core, color=ORANGE, lw=1.8, solid_capstyle="round",
                    label="core-cardiac cells")
        ax.set_title(f"{LABEL.get(st, st)}: {int(f.z_levels)} levels, step {f.z_step_median:g}", fontsize=9)
        ax.set_ylim(0, max(1.05, float(g.n_cells.max() / f.peak_section_cells) + 0.05))
        ax.set_xlabel("section index (stack order)", fontsize=8)
        style(ax)
    axes[0, 0].set_ylabel("cells per level / largest section", fontsize=8)
    axes[1, 0].set_ylabel("cells per level / largest section", fontsize=8)
    leg = axes.flat[-1]
    leg.axis("off")
    h = [plt.Line2D([], [], color=BLUE, lw=2), plt.Line2D([], [], color=ORANGE, lw=2)]
    leg.legend(h, ["all cells", "core-cardiac cells (CM + endocardium)"], loc="center left", frameon=False,
               fontsize=9, labelcolor=INK)
    leg.text(0.0, 0.12, "Duplicated records dropped. The KO level at z = -320\nholds two sections, so it rises "
             "above 1.", fontsize=8.5, color=INK2, transform=leg.transAxes)
    fig.suptitle("Cells per z level, every released MERFISH file", color=INK, fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(d / "fig_z_profiles.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_pitch(d: Path):
    xy = pd.read_csv(d / "variogram_xy.csv")
    vz = pd.read_csv(d / "variogram_z.csv")
    fit = pd.read_csv(d / "zpitch_fit.csv")
    fit = fit[(fit.measure == "mism_group") & (fit.k_set == "1-3")].set_index("stage")
    stages = ["E8.25_late", "E8.75", "E9.5", "E9.5_mab21l2_ko"]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), facecolor=SURFACE, sharey=True)
    for ax, st in zip(axes, stages, strict=True):
        g = xy[(xy.stage == st) & (xy.n_pairs >= 200) & (xy.d_mean >= 4) & (xy.d_mean <= 120)]
        ax.plot(g.d_mean, g.mism_group, color=MUTED, lw=1.8, label="same section, xy distance")
        z = vz[vz.stage == st]
        step, p = float(fit.loc[st, "step"]), float(fit.loc[st, "p_best"])
        r_nom = np.sqrt(z.mean_dxy2 + (z.k * step) ** 2)
        r_fit = np.sqrt(z.mean_dxy2 + (z.k * p) ** 2)
        ax.plot(r_nom, z.mism_group, "o", ms=6, color=BLUE, mec=SURFACE, mew=1.2,
                label="k sections apart, at k x stored step")
        ax.plot(r_fit, z.mism_group, "s", ms=6, color=ORANGE, mec=SURFACE, mew=1.2,
                label="k sections apart, at k x fitted pitch")
        ax.set_title(f"{LABEL.get(st, st)}: step {step:g}, fitted pitch {p:g}", fontsize=9)
        ax.set_xlabel("distance (xy units)", fontsize=8)
        ax.set_xlim(0, 160)
        style(ax)
    axes[0].set_ylabel("P(different lineage group)", fontsize=8)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, fontsize=8.5, labelcolor=INK)
    fig.suptitle("Isotropy check: pairs k levels apart (paired by z difference) against the in-plane curve",
                 color=INK, fontsize=11, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(d / "fig_pitch_collapse.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_sizes(d: Path):
    G = pd.read_csv(d / "stage_geometry.csv").set_index("name")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.8), facecolor=SURFACE)
    heart = ["E8.25_late", "E8.75", "E9.5"]
    t = G.loc[heart, "t"]
    a1.plot(t, G.loc[heart, "rms_full"], "-o", color=BLUE, lw=1.8, ms=6, mec=SURFACE, mew=1.2, label="full cloud")
    a1.plot(t, G.loc[heart, "rms_core"], "-s", color=ORANGE, lw=1.8, ms=6, mec=SURFACE, mew=1.2,
            label="core-cardiac cells")
    ko = G.loc["E9.5_mab21l2_ko"]
    a1.plot([9.5], [ko.rms_full], "o", ms=7, mfc="none", mec=BLUE, mew=1.6,
            label="Mab21l2 KO (second E9.5 embryo; informational)")
    a1.plot([9.5], [ko.rms_core], "s", ms=7, mfc="none", mec=ORANGE, mew=1.6)
    hp = pd.read_csv(d / "size_heart_pseudo.csv")
    pr = hp[(hp.variant == "raw") & (hp.series == "core") & (hp.rule != "copy_last")]
    a1.plot([9.5] * len(pr), pr.pred_rms, "^", ms=6, color=MUTED, label="E9.5 core predicted from E8.25_late, E8.75")
    a1.set_yscale("log")
    a1.set_xlabel("nominal stage (E-day)", fontsize=8)
    a1.set_ylabel("raw 3D RMS radius (log scale)", fontsize=8)
    a1.set_title("Heart setting: core flat over E8.25_late-E8.75, then a jump", fontsize=9)
    a1.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    plain_log_ticks(a1, [150, 200, 250, 300, 350])
    style(a1)
    emb = ["E6.75", "E7.25", "E8.0"]
    a2.plot(G.loc[emb, "t"], G.loc[emb, "rms_full"], "-o", color=BLUE, lw=1.8, ms=6, mec=SURFACE, mew=1.2,
            label="full cloud (raw 3D)")
    a2.plot(G.loc[emb, "t"], G.loc[emb, "rms_xy_full"], "-s", color=ORANGE, lw=1.8, ms=6, mec=SURFACE, mew=1.2,
            label="in-plane (xy) RMS")
    a2.plot(G.loc[emb, "t"], G.loc[emb, "sd_z_full"], "-D", color=AQUA, lw=1.8, ms=5.5, mec=SURFACE, mew=1.2,
            label="SD of z")
    se = pd.read_csv(d / "size_embryo.csv")
    lo = se[(se.z_mode == "raw") & (se.check == "loo_interior")].iloc[0]
    a2.plot([7.25], [lo.pred_rms], "^", ms=7, color=MUTED, label="E7.25 interpolated from E6.75, E8.0")
    a2.set_yscale("log")
    a2.set_xlabel("nominal stage (E-day)", fontsize=8)
    gxy = G.loc[emb[-1], "rms_xy_full"] / G.loc[emb[0], "rms_xy_full"]
    gz = G.loc[emb[-1], "sd_z_full"] / G.loc[emb[0], "sd_z_full"]
    a2.set_title(f"Embryo setting: xy RMS x{gxy:.1f}, z SD x{gz:.2f} from {emb[0]} to {emb[-1]}", fontsize=9)
    a2.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    plain_log_ticks(a2, [70, 100, 150, 200, 300])
    style(a2)
    fig.tight_layout()
    fig.savefig(d / "fig_size_trajectories.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="runs/coords")
    d = Path(ap.parse_args().dir)
    fig_z_profiles(d)
    fig_pitch(d)
    fig_sizes(d)
    print("wrote figures to", d)
