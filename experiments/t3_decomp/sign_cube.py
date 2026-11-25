"""The T3 sign cube's algebra (ledger row R72, the captain's order; T3_MOON5_design.md section 2.2 as the fixes F1
to F10 of its two cold checks word it): pure functions of arrays, shared by cube_tables.py (the tables and G0),
k31_cube_read.py (row K31's counted read) and build_sign_cube.py (the builder). Nothing here loads a file.

Tables (each a t3_decompose.gene_factor on released files, so already rescaled and clipped to [0.25, 4]):
  l_M = ln transfer_mab21l2_all, l_G = ln gata4_pooled (C2's two pinned tables);
  S_M, Y_M = ln (WT E9.5, WT E8.0 over WT E8.75) on l_M's den; S_G, Y_G = the same on l_G's den.
Axes (design 2.2; c is C2's, SD(l_M) / SD(l_M + l_G)):
  X_S = c (S_M - Y_M) / 2           the stage bit: older pole E9.5, younger pole E8.0
  X_M = c (l_M - S_M)               the Mab21l2 knockout against WT E9.5, on l_M's labels
  X_G = c (l_G - (S_G + Y_G) / 2)   the Gata4 knockout against the stage-neutral WT midpoint
  A   = c [(S_M + Y_M) / 2 + (S_G + Y_G) / 2],  A + X_S + X_M + X_G = c (l_M + l_G) = ln phi_C2
A vertex v = (s_S, s_M, s_G), letters p (C2's sign) or n (reflected), F1 and F7:
  ln phi_v = ln phi_C2 + sum_j (k_v s_j - 1) Q(X_j),  Q(X) = X - B (B'B)^-1 B'X,
  B's columns m_L * phi_C2 for every label L with a cell in the draw (m_L: L's count-weighted mean within-cell
  proportions there, the counts as drawn, before the factor step; row R72's words, "every label's expected library
  is held within 1% of C2's", bind: F1 named the labels holding >= 1% of the drawn cells, which left the smaller
  labels free to move once there are two of them or the den is not the drawn cells' mean); k_v > 0 holds SD(ln phi_v) at
  SD(ln phi_C2) (the root in [0.5, 2] nearest 1; refused when there is none; |k_v - 1| <= 1e-9 is taken as
  k_v = 1);
  the clip gene by gene to [min(0.25, phi_C2), max(4, phi_C2)] (C2's builder never re-clips its combined factor);
  one constant so that den . phi_v = den . phi_C2; then rounds of a Newton correction within B's span, taken jointly
  with the den constant, each followed by the clip and the exact den constant, until every label's expected library
  m_L . phi is C2's to 1e-6 and the den total C2's to 1e-9, at most 20 rounds; else the vertex is not built (a label
  too small to be held is refused, never passed). Every label's library is then checked within 1% (LIB_BAND).
The identity vertex ppp is C2's own array by construction (every coefficient k s_j - 1 is exactly 0, so the factor is
phi_C2 times exp(0); the clip, the den constant and the hold then leave it unchanged), with no special case.
G1's analogue (row K31) uses the same function with the base phi_M, the axes X_S' = (S_M - Y_M) / 2 and
X_M' = l_M - S_M, and B's columns m_L * phi_M for every label of K18's drawn cells.
"""

from __future__ import annotations

import hashlib

import numpy as np

AXES = ("S", "M", "G")
VERTICES = tuple(a + b + c for a in "pn" for b in "pn" for c in "pn")  # ppp, ppn, pnp, pnn, npp, npn, nnp, nnn
IDENTITY = "ppp"
NAMED = tuple(v for v in VERTICES if v != IDENTITY)  # the seven files of row R72
CLIP = (0.25, 4.0)
K_RANGE = (0.5, 2.0)  # k_v refused outside it (design 2.2)
K_SNAP = 1e-9  # F3: |k - 1| <= 1e-9 is taken as k = 1
LIB_TOL = 1e-6  # F7: each held label's expected library at C2's to 1e-6
DEN_TOL = 1e-9  # F1: den . phi_v = den . phi_C2 within 1e-9
MAX_ROUNDS = 20  # F7
LABEL_SHARE = 0.0  # every label with a drawn cell is held (row R72's "every label"; F1's 1% left smaller ones free)
LIB_BAND = 0.01  # F1's G0d bar and row R72: every label's expected library within 1% of C2's
RANK_TOL = 1e-12  # B's rank: singular values above this share of the largest


class CubeRefused(Exception):
    """A vertex (or an arm) whose factor cannot be made by the rule: not built."""


def signs(v: str) -> tuple[int, ...]:
    """The signs of a vertex string of p and n letters (p: +1, C2's sign; n: -1, reflected)."""
    if not isinstance(v, str) or not v or any(ch not in "pn" for ch in v):
        raise CubeRefused(f"vertex {v!r} is not a string of p and n")
    return tuple(1 if ch == "p" else -1 for ch in v)


def sd(v) -> float:
    """Population SD over the genes (ddof 0), as build_two_knockout's c."""
    return float(np.std(np.asarray(v, dtype=np.float64), ddof=0))


def factor_sha256(phi) -> str:
    """The sha256 of a factor's float64 bytes (build_two_knockout.factor_sha256's rule)."""
    return hashlib.sha256(np.ascontiguousarray(phi, dtype=np.float64).tobytes()).hexdigest()


def cube_axes(l_m, l_g, s_m, y_m, s_g, y_g, c: float) -> dict[str, np.ndarray]:
    """The anchor and the three axes of design 2.2 from the six log tables and C2's c."""
    l_m, l_g, s_m, y_m, s_g, y_g = (np.asarray(x, dtype=np.float64) for x in (l_m, l_g, s_m, y_m, s_g, y_g))
    return {"A": c * ((s_m + y_m) / 2.0 + (s_g + y_g) / 2.0), "S": c * (s_m - y_m) / 2.0, "M": c * (l_m - s_m),
            "G": c * (l_g - (s_g + y_g) / 2.0)}


def analogue_axes(l_m, s_m, y_m) -> dict[str, np.ndarray]:
    """G1's analogue cube without a Gata4 table (design 2.6): A' + X_S' + X_M' = l_M."""
    l_m, s_m, y_m = (np.asarray(x, dtype=np.float64) for x in (l_m, s_m, y_m))
    return {"A": (s_m + y_m) / 2.0, "S": (s_m - y_m) / 2.0, "M": l_m - s_m}


def label_rows(counts, lab, share: float = LABEL_SHARE) -> tuple[np.ndarray, np.ndarray]:
    """The labels holding >= `share` of the drawn cells (by default every label with a drawn cell) and their rows m_L:
    the label's counts summed over its cells over its total count (the count-weighted mean within-cell proportions).
    Returns (labels, M), M labels x genes."""
    C = np.asarray(counts, dtype=np.float64)
    lab = np.asarray(lab)
    if C.ndim != 2 or lab.shape != (C.shape[0],) or C.shape[0] == 0:
        raise CubeRefused("the drawn cells' counts and labels do not match")
    labels, n = np.unique(lab, return_counts=True)
    held = labels[n >= share * C.shape[0]]
    rows = []
    for t in held:
        tot = C[lab == t].sum(0)
        s = float(tot.sum())
        if not s > 0:
            raise CubeRefused(f"label {t!r} holds no counts")
        rows.append(tot / s)
    if not rows:
        raise CubeRefused("no label holds the share of the drawn cells")
    return held, np.vstack(rows)


def basis(M, phi) -> np.ndarray:
    """B: genes x labels, columns m_L * phi (F1)."""
    return (np.asarray(M, dtype=np.float64) * np.asarray(phi, dtype=np.float64)[None, :]).T


def projector(B):
    """Q(X) = X - B (B'B)^-1 B'X, the orthogonal projection off B's column span (through an orthonormal basis of the
    span from the SVD, rank by RANK_TOL)."""
    B = np.asarray(B, dtype=np.float64)
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    r = int((s > s[0] * RANK_TOL).sum()) if s.size and s[0] > 0 else 0
    U = U[:, :r]

    def Q(X):
        X = np.asarray(X, dtype=np.float64)
        return X - U @ (U.T @ X)

    return Q


def solve_k(b, qs, sg, k_range=K_RANGE) -> tuple[float, dict]:
    """k > 0 with SD(b + sum_j (k s_j - 1) q_j) = SD(b): ln phi(k) = u + k w, u = b - sum q_j, w = sum s_j q_j, a
    quadratic in k. The positive root inside k_range nearest 1 (|k - 1| <= K_SNAP is taken as 1); CubeRefused when
    no positive root lies inside k_range."""
    b = np.asarray(b, dtype=np.float64)
    qs = [np.asarray(q, dtype=np.float64) for q in qs]
    u = b - np.sum(qs, axis=0)
    w = np.sum([s * q for s, q in zip(sg, qs, strict=True)], axis=0)
    uc, wc = u - u.mean(), w - w.mean()
    vw, cuw, vu, vb = float(wc @ wc) / w.size, float(uc @ wc) / w.size, float(uc @ uc) / u.size, sd(b) ** 2
    rec = {"var_w": vw, "cov_uw": cuw, "var_u": vu, "var_b": vb}
    if not vw > 0:
        raise CubeRefused("the flipped axes are zero: k is undefined")
    disc = cuw * cuw - vw * (vu - vb)
    if 0 > disc >= -1e-12 * max(cuw * cuw, vw * abs(vu - vb), 1e-300):
        disc = 0.0  # a double root (rounding only)
    if disc < 0:
        raise CubeRefused("no real k holds SD(ln phi) at the base's")
    roots = sorted({(-cuw + q * np.sqrt(disc)) / vw for q in (1.0, -1.0)})
    pos = [1.0 if abs(r - 1.0) <= K_SNAP else float(r) for r in roots if r > 0]
    rec["roots"] = [float(r) for r in roots]
    if not pos:
        raise CubeRefused("no positive k holds SD(ln phi) at the base's")
    inside = [r for r in pos if k_range[0] <= r <= k_range[1]]
    if not inside:
        k = min(pos, key=lambda r: (abs(r - 1.0), r))
        raise CubeRefused(f"k = {k:.6g} lies outside [{k_range[0]}, {k_range[1]}]")
    k = min(inside, key=lambda r: (abs(r - 1.0), r))  # the admissible root nearest 1
    rec["k"] = k
    return k, rec


def libraries(M, phi) -> np.ndarray:
    """Each held label's expected library m_L . phi (relative to its drawn total)."""
    return np.asarray(M, dtype=np.float64) @ np.asarray(phi, dtype=np.float64)


def hold(phi, phi_ref, M, den, B, lo, hi, rounds: int = MAX_ROUNDS) -> tuple[np.ndarray, dict]:
    """F7's hold: until every held label's m_L . phi is phi_ref's to LIB_TOL and den . phi is phi_ref's to DEN_TOL,
    rounds of one Newton step in log space within B's span with the den constant as one more unknown (a uniform
    scale, solved jointly so that the den step does not undo the labels; the Jacobian over the genes strictly inside
    their clip bounds, since the clip holds the others), then the clip to [lo, hi], then the exact den constant; at
    most `rounds` rounds. A factor already held is returned unchanged (round 0)."""
    phi = np.asarray(phi, dtype=np.float64).copy()
    M = np.asarray(M, dtype=np.float64)
    den = np.asarray(den, dtype=np.float64)
    e0, d0 = libraries(M, phi_ref), float(den @ np.asarray(phi_ref, dtype=np.float64))
    L = M.shape[0]
    for r in range(rounds + 1):
        lib = libraries(M, phi) / e0 - 1.0
        dd = float(den @ phi) / d0 - 1.0
        if float(np.abs(lib).max()) <= LIB_TOL and abs(dd) <= DEN_TOL:
            return phi, {"rounds": r, "converged": True, "max_lib_dev": float(np.abs(lib).max()), "den_dev": dd}
        if r == rounds:
            break
        # the genes the clip holds at a bound do not move under the step (an active set): the Jacobian counts only
        # the genes strictly inside their bounds, so a step is not undone by the clip that follows it
        free = phi * ((phi > lo * (1.0 + 1e-12)) & (phi < hi * (1.0 - 1e-12)))
        J = np.empty((L + 1, L + 1))
        J[:L, :L] = (M * free[None, :]) @ B / e0[:, None]
        J[:L, L] = libraries(M, free) / e0
        J[L, :L] = ((den * free) @ B) / d0
        J[L, L] = float(den @ free) / d0
        sol = np.linalg.lstsq(J, -np.concatenate([lib, [dd]]), rcond=None)[0]
        phi = phi * np.exp(B @ sol[:L] + sol[L])
        phi = np.clip(phi, lo, hi)
        phi = phi * (d0 / float(den @ phi))
    lib = libraries(M, phi) / e0 - 1.0
    return phi, {"rounds": rounds, "converged": False, "max_lib_dev": float(np.abs(lib).max()),
                 "den_dev": float(den @ phi) / d0 - 1.0}


def flip_factor(phi_base, qs, sg, M, den, B=None, k_range=K_RANGE, rounds: int = MAX_ROUNDS) -> tuple[np.ndarray,
                                                                                                         dict]:
    """The factor of one sign pattern on a base (C2's factor for the cube; phi_M for G1's analogue), F1 and F7:
    phi_base * exp(sum_j (k s_j - 1) q_j) with q_j = Q(X_j) already projected and k from solve_k; the clip to
    [min(0.25, phi_base), max(4, phi_base)]; the den constant; then hold(). CubeRefused when k has no admissible root
    or the hold does not converge. With every sign +1 the result is phi_base itself (k = 1, every coefficient 0)."""
    base = np.asarray(phi_base, dtype=np.float64)
    if not np.all(np.isfinite(base) & (base > 0)):
        raise CubeRefused("the base factor is not finite and positive")
    den = np.asarray(den, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)
    B = basis(M, base) if B is None else np.asarray(B, dtype=np.float64)
    lb = np.log(base)
    k, krec = solve_k(lb, qs, sg, k_range)
    coef = [k * s - 1.0 for s in sg]
    delta = np.zeros_like(lb)
    for a, q in zip(coef, qs, strict=True):
        delta = delta + a * np.asarray(q, dtype=np.float64)
    phi = base * np.exp(delta)
    lo, hi = np.minimum(CLIP[0], base), np.maximum(CLIP[1], base)
    pre_sd = sd(np.log(phi))
    n_lo, n_hi = int((phi < lo).sum()), int((phi > hi).sum())
    phi = np.clip(phi, lo, hi)
    phi = phi * (float(den @ base) / float(den @ phi))
    phi, hrec = hold(phi, base, M, den, B, lo, hi, rounds)
    if not np.all(np.isfinite(phi) & (phi > 0)):
        raise CubeRefused("the factor is not finite and positive")
    lib = libraries(M, phi) / libraries(M, base) - 1.0
    rec = {"k": k, "k_solve": krec, "coefficients": coef, "sd_ln_base": sd(lb), "sd_ln_before_clip": pre_sd,
           "sd_ln_factor": sd(np.log(phi)), "clipped_low": n_lo, "clipped_high": n_hi, "hold": hrec,
           "held_labels": int(M.shape[0]), "max_abs_library_dev": float(np.abs(lib).max()),
           "den_dev": float(den @ phi) / float(den @ base) - 1.0, "min": float(phi.min()), "max": float(phi.max()),
           "factor_sha256": factor_sha256(phi)}
    if not hrec["converged"]:
        raise CubeRefused(f"the hold did not converge in {rounds} rounds (largest label deviation "
                          f"{hrec['max_lib_dev']:.3g}, den {hrec['den_dev']:.3g})")
    return phi, rec


def vertex_factor(phi_c2, q_axes: dict, v: str, M, den, B=None) -> tuple[np.ndarray, dict]:
    """The cube's vertex v on C2's factor: flip_factor with the three projected axes in AXES order."""
    sg = signs(v)
    if len(sg) != len(AXES) or v not in VERTICES:
        raise CubeRefused(f"vertex {v!r} is not one of {list(VERTICES)}")
    phi, rec = flip_factor(phi_c2, [q_axes[a] for a in AXES], sg, M, den, B)
    return phi, rec | {"vertex": v, "signs": dict(zip(AXES, sg, strict=True))}


def within_band(rec: dict, band: float = LIB_BAND) -> bool:
    """G0d's bar (F1) and row R72: every held label's expected library (every label with a drawn cell, label_rows'
    default) within `band` of the base's."""
    return bool(rec["max_abs_library_dev"] <= band)


def stage_positions(tables: dict, young: str = "E8.0", old: str = "E9.5") -> dict[str, float]:
    """F4's G0g: s(x) = <ln x - ln young, ln old - ln young> / |ln old - ln young|^2 for each table x (log tables over
    WT E8.75 on l_M's den; WT E8.75's own table is 0)."""
    y, o = np.asarray(tables[young], dtype=np.float64), np.asarray(tables[old], dtype=np.float64)
    d = o - y
    dd = float(d @ d)
    if not dd > 0:
        raise CubeRefused("the stage axis is zero")
    return {x: float((np.asarray(t, dtype=np.float64) - y) @ d / dd) for x, t in tables.items()}


def g0g_passes(pos: dict) -> bool:
    """F4's bar: 0 < s(E8.25_late) < s(E8.75) < 1, |s(E8.75) - 0.5| <= 0.25, and E8.5 between E8.25_late and E8.75
    where its table is there."""
    a, c = pos.get("E8.25_late"), pos.get("E8.75")
    if a is None or c is None or not (np.isfinite(a) and np.isfinite(c)):
        return False
    ok = 0.0 < a < c < 1.0 and abs(c - 0.5) <= 0.25
    b = pos.get("E8.5")
    if b is not None:
        ok = ok and bool(np.isfinite(b)) and a < b < c
    return bool(ok)


def kept_vertices(kept_axes, allowed=NAMED) -> list[str]:
    """The named vertices whose every reflected axis is in `kept_axes` (an axis held at '+' is never reflected)."""
    keep = set(kept_axes)
    return [v for v in allowed if all(ch == "p" or ax in keep for ch, ax in zip(v, AXES, strict=True))]
