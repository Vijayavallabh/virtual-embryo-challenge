"""How was each 3D stack built? Per section, fit spatial_3D[:, :2] from the file's 2D coordinates (DGX).

For every z level we fit (a) a pure translation, (b) a rigid transform (rotation + translation) and (c) a similarity
transform (adds one isotropic scale), by least squares with known cell correspondence. A rigid fit with ~zero residual
means the section was placed by a rigid transform (xy keeps the 2D units); a large rigid residual means the section
was deformed (non-rigid alignment). The fitted scale says whether 3D xy and 2D share units.

2D source: obsm["spatial"] (slide coordinates) when present (Mab21l2 KO), else obsm["spatial_2D"].
Writes runs/coords/section_transforms.csv. Usage: uv run python experiments/coords/section_transforms.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coordlib import FILES, snap_z  # noqa: E402


def fit(src: np.ndarray, dst: np.ndarray, mode: str) -> tuple[float, float, float]:
    """Least-squares fit dst ~ s R src + t. Returns (rms residual, scale, rotation angle in degrees)."""
    ms, md = src.mean(0), dst.mean(0)
    A, B = src - ms, dst - md
    if mode == "translation":
        return float(np.sqrt(((A - B) ** 2).sum(1).mean())), 1.0, 0.0
    U, S, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, d])
    R = (U @ D @ Vt).T  # proper rotation (reflection allowed only if d < 0 is forced: report it via angle)
    scale = float((S * np.diag(D)).sum() / (A ** 2).sum()) if mode == "similarity" else 1.0
    res = B - scale * A @ R.T
    ang = float(np.degrees(np.arctan2(R[1, 0], R[0, 0])))
    return float(np.sqrt((res ** 2).sum(1).mean())), scale, ang


def reflection_needed(src: np.ndarray, dst: np.ndarray) -> bool:
    A, B = src - src.mean(0), dst - dst.mean(0)
    U, _, Vt = np.linalg.svd(A.T @ B)
    return bool(np.linalg.det(U @ Vt) < 0)


def main(out: Path) -> None:
    rows = []
    for stage, path in FILES.items():
        a = ad.read_h5ad(path, backed="r")
        C3 = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)
        key = "spatial" if "spatial" in a.obsm else "spatial_2D"
        C2 = np.asarray(a.obsm[key], dtype=np.float64)[:, :2]
        a.file.close()
        lev, levels, _ = snap_z(C3[:, 2])
        for g in range(levels.size):
            m = lev == g
            if m.sum() < 10:
                continue
            src, dst = C2[m], C3[m, :2]
            r_t, _, _ = fit(src, dst, "translation")
            r_r, _, ang = fit(src, dst, "rigid")
            r_s, sc, _ = fit(src, dst, "similarity")
            # mirrored rigid fit (flip src x) to detect a reflection between 2D and 3D
            src_m = src * np.array([-1.0, 1.0])
            r_rm, _, _ = fit(src_m, dst, "rigid")
            rows.append({"stage": stage, "source_2d": key, "level": g, "z": float(levels[g]), "n": int(m.sum()),
                         "rms_dst": float(np.sqrt(((dst - dst.mean(0)) ** 2).sum(1).mean())),
                         "resid_translation": r_t, "resid_rigid": r_r, "resid_rigid_mirrored": r_rm,
                         "resid_similarity": r_s, "scale_similarity": sc, "rotation_deg": ang,
                         "svd_wants_reflection": reflection_needed(src, dst)})
        df = pd.DataFrame([r for r in rows if r["stage"] == stage])
        print(f"{stage} ({key}): levels {len(df)}; median resid translation {df.resid_translation.median():.3g}, "
              f"rigid {df.resid_rigid.median():.3g}, rigid-mirrored {df.resid_rigid_mirrored.median():.3g}, "
              f"similarity {df.resid_similarity.median():.3g}; scale median {df.scale_similarity.median():.4f} "
              f"[{df.scale_similarity.min():.4f}, {df.scale_similarity.max():.4f}]; rotation range "
              f"{df.rotation_deg.min():.1f}..{df.rotation_deg.max():.1f}; mirrored-better levels "
              f"{int((df.resid_rigid_mirrored < df.resid_rigid).sum())}", flush=True)
    pd.DataFrame(rows).to_csv(out / "section_transforms.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/coords")
    main(Path(ap.parse_args().out))
