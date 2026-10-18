"""How was each 3D stack built? Per section, fit spatial_3D[:, :2] from the file's 2D coordinates (DGX).

For every physical section we fit (a) a pure translation, (b) a rigid transform (rotation + translation) and (c) a
similarity transform (adds one isotropic scale), by least squares with known cell correspondence. A rigid fit with
~zero residual means the section was placed by a rigid transform (xy keeps the 2D units); a large rigid residual means
the section was deformed (non-rigid alignment). The fitted scale says whether 3D xy and 2D share units.

A physical section is a z level, split by obs["sample_id"] where the file records it (the Mab21l2 KO stores two
sections at z = -320). Duplicated records are dropped first (vec.prep.load_stage(..., drop_duplicates=True)).
`rigid` flags sections whose rigid residual is at most RIGID_TOL.

2D source: obsm["spatial"] (slide coordinates) when present (Mab21l2 KO), else obsm["spatial_2D"].
Writes runs/coords/section_transforms.csv.
Usage: uv run python experiments/coords/section_transforms.py [--out runs/coords] [--stage NAME=PATH ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coordlib import FILES, load, snap_z, stage_specs  # noqa: E402

RIGID_TOL = 0.01


def fit(src: np.ndarray, dst: np.ndarray, mode: str) -> tuple[float, float, float]:
    """Least-squares fit dst ~ s R src + t. Returns (rms residual, scale, rotation angle in degrees)."""
    ms, md = src.mean(0), dst.mean(0)
    A, B = src - ms, dst - md
    if mode == "translation":
        return float(np.sqrt(((A - B) ** 2).sum(1).mean())), 1.0, 0.0
    U, S, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, d])
    R = (U @ D @ Vt).T  # proper rotation (a reflection is never fitted; see resid_rigid_mirrored)
    scale = float((S * np.diag(D)).sum() / (A ** 2).sum()) if mode == "similarity" else 1.0
    res = B - scale * A @ R.T
    ang = float(np.degrees(np.arctan2(R[1, 0], R[0, 0])))
    return float(np.sqrt((res ** 2).sum(1).mean())), scale, ang


def reflection_needed(src: np.ndarray, dst: np.ndarray) -> bool:
    A, B = src - src.mean(0), dst - dst.mean(0)
    U, _, Vt = np.linalg.svd(A.T @ B)
    return bool(np.linalg.det(U @ Vt) < 0)


def main(out: Path, files: dict[str, str]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for stage, path in files.items():
        a, dd = load(path, drop_duplicates=True)
        C3 = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)
        key = "spatial" if "spatial" in a.obsm else "spatial_2D"
        C2 = np.asarray(a.obsm[key], dtype=np.float64)[:, :2]
        lev, levels, _ = snap_z(C3[:, 2])
        sid = np.asarray(a.obs["sample_id"]).astype(str) if "sample_id" in a.obs else np.full(len(lev), "-")
        for (g, s), idx in pd.Series(np.arange(len(lev))).groupby([lev, sid]).groups.items():
            idx = np.asarray(idx)
            if idx.size < 10:
                continue
            src, dst = C2[idx], C3[idx, :2]
            r_t, _, _ = fit(src, dst, "translation")
            r_r, _, ang = fit(src, dst, "rigid")
            r_s, sc, _ = fit(src, dst, "similarity")
            # mirrored rigid fit (flip src x) to detect a reflection between 2D and 3D
            r_rm, _, _ = fit(src * np.array([-1.0, 1.0]), dst, "rigid")
            rows.append({"stage": stage, "source_2d": key, "level": int(g), "z": float(levels[g]), "sample_id": s,
                         "n": int(idx.size), "rms_dst": float(np.sqrt(((dst - dst.mean(0)) ** 2).sum(1).mean())),
                         "resid_translation": r_t, "resid_rigid": r_r, "resid_rigid_mirrored": r_rm,
                         "resid_similarity": r_s, "scale_similarity": sc, "rotation_deg": ang,
                         "rigid": bool(r_r <= RIGID_TOL), "svd_wants_reflection": reflection_needed(src, dst)})
        df = pd.DataFrame([r for r in rows if r["stage"] == stage])
        nonrigid = [(r.z, r.sample_id, round(r.resid_rigid, 2), round(r.scale_similarity, 3))
                    for r in df[~df.rigid].itertuples()]
        print(f"{stage} ({key}): sections {len(df)} (dropped {dd['n_dropped']} duplicates); rigid {int(df.rigid.sum())}"
              f"/{len(df)}; median resid translation {df.resid_translation.median():.3g}, rigid "
              f"{df.resid_rigid.median():.3g}, similarity {df.resid_similarity.median():.3g}; scale median "
              f"{df.scale_similarity.median():.4f} [{df.scale_similarity.min():.4f}, {df.scale_similarity.max():.4f}]; "
              f"rotation range {df.rotation_deg.min():.1f}..{df.rotation_deg.max():.1f}; mirrored-better sections "
              f"{int((df.resid_rigid_mirrored < df.resid_rigid).sum())}; non-rigid (z, sample, resid, scale): "
              f"{nonrigid}", flush=True)
    pd.DataFrame(rows).to_csv(out / "section_transforms.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    ap.add_argument("--stage", action="append", help="NAME=PATH (repeatable; default: every released file)")
    args = ap.parse_args()
    main(Path(args.out), stage_specs(args.stage, FILES))
