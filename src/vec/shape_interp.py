"""Shape interpolation of the embryo's exact-position cloud between its two released brackets (row E49 of
experiments/harness/TRIES.md, as drafted; branch campaign/embryo-shape-interp).

Every function here is pure and float64. Inputs are the two released bracket clouds and the exact-position cloud, each
in the earlier bracket's centred, RMS-1 frame as `vec.t2_build.interpolate_clouds` holds them (Zb, Za, Z), the
calendar fraction t, and, for the registered forms, the build's own coupling rotation R. Outputs are new coordinates
for the same rows, in the same order, and the source cloud the cloud guard checks them against. The forms move
positions only: rows, order, count and expression are the caller's and are never touched. The size rule's RMS is set
by the caller after the form, so the scale term cannot move.

Forms (the builder exposes these three only):
  stretch_ratio       F1: a stretch along the cloud's own principal axes to the calendar-fraction log interpolation of
                      the two brackets' trace-normalised principal SDs (paired by descending order, the scorer's own
                      pairing). Adjacent axes are stretched apart only if both brackets identify them (`identifiable`);
                      otherwise they form a block that keeps the cloud's own shape inside it (one factor).
  stretch_registered  F2: the symmetric linear optimal-transport map that gives the cloud the trace-normalised
                      log-Euclidean interpolation of the earlier bracket's covariance and the later bracket's covariance
                      in the earlier frame (Za @ R). Refused if its smallest or largest normalised SD leaves the range
                      of the brackets' own (`ShapeRefused`).
  overlay             F3: every later-bracket row at its own real position in the earlier frame ((Za @ R)[row]), the
                      earlier-bracket rows where they are, at the file's own share.
Controls (harness candidates only; the builder never builds one; no bracket-range guard):
  stretch_ratio_anti       F1 with the step mirrored about the earlier bracket: (1 + t) log s_b - t log s_a.
  stretch_registered_anti  F2 with the covariance step mirrored: expm((1 + t) logm S_b - t logm S_a), the same R.
  stretch_registered_rot   F2's interpolated covariance turned by Q (Q S* Q^T): the same spectrum, a random
                           orientation; Q is a uniform proper rotation from default_rng([seed, ROT_STREAM]).
  overlay_rot              F3 with Q in place of R.
Row E53 (as specified; branch campaign/moon-embryo-geodesic) adds two harness-only forms, GEODESIC_FORMS, outside
ALL_FORMS (the builder does not offer them): geodesic_overlay (F3's rows on one LDDMM geodesic between the brackets,
vec.lddmm) and its placebo geodesic_overlay_rot (the fit started from Q in place of R). They need every cell's joint PCs
(apply_form's `features`).
Row E55 (as specified; branch campaign/emb-osr) adds COMPOSED_FORMS (overlay_stretch, F3 then F1's stretch on the
overlay's union) and its harness-only placebos COMPOSED_CONTROLS (overlay_stretch_anti, overlay_stretch_axes), also
outside ALL_FORMS; the builder offers overlay_stretch only.
Row R71 (the captain's order of 2026-11-25; branch campaign/embryo-bf) adds BAYES_FORMS (overlay_bayes, BF: the
overlay's union carried by a smoothed optimal-transport map onto a consensus canonical density of released specimens,
EMBRYO_MOON2_design.md's P1 with its cold challenge's fixes) and its harness-only controls BAYES_CONTROLS
(overlay_bayes_brackets: the prior is the task's two brackets only; overlay_bayes_noflip: no flip copies, each specimen
turned onto the union's canonical frame by the best of the four proper flips), outside ALL_FORMS; the builder offers
overlay_bayes only, and only under row R71. No step of BF branches on the cell count.

Fixed constants (protocol, not tuned): IDENT_DRAWS draws of IDENT_FRAC of a bracket's cells (the emulator's 10 percent
subsample halved, the truth half's share), from the seed-free streams default_rng([IDENT_STREAM, h]), so a pair of
released brackets has one block structure at every build seed and at the test. ROT_STREAM = 50, which no other build
step uses (48 is the row order's COPY_PARITY_STREAM, 47 the floor-exact draw).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

FORMS = ("stretch_ratio", "stretch_registered", "overlay")
CONTROLS = ("stretch_ratio_anti", "stretch_registered_anti", "stretch_registered_rot", "overlay_rot")
ALL_FORMS = (*FORMS, *CONTROLS)
CONTROL_OF = {"stretch_ratio": ("stretch_ratio_anti",),
              "stretch_registered": ("stretch_registered_anti", "stretch_registered_rot"),
              "overlay": ("overlay_rot",)}
IDENT_DRAWS = 20
IDENT_FRAC = 0.05  # vec.evaluate.SUBSAMPLE_FRAC / 2: the truth half's share of a stage (tests pin it)
IDENT_STREAM = 49
ROT_STREAM = 50
MIN_CLOUD = 4  # vec.t2_build.MIN_CLOUD: fewest points for a PCA frame (tests pin it)
ROTATION_TOL = 1e-6
GUARD_SLACK = 1e-12  # relative: the bracket-range guard's rounding allowance


class ShapeRefused(ValueError):
    """A form refuses to build (F2's bracket-range guard, a non-finite or improper rotation)."""


@dataclass
class ShapeResult:
    """`apply_form` result: Z the new positions (same rows, same order), pool the source cloud for the cloud guard (in
    the same frame), pool_weights its draw weights (None: uniform over pool's rows), record flat values only."""

    Z: np.ndarray
    pool: np.ndarray
    pool_weights: np.ndarray | None
    record: dict


# ------------------------------------------------------------------------------------------------ second moments
def principal(Z) -> tuple[np.ndarray, np.ndarray]:
    """The scorer's principal frame of a cloud (shape_metrics._canonicalise without the sign fix): the SVD of the
    centred cloud. Returns (sigma, Vt): sigma the principal SDs (s / sqrt(N)) renormalised so that sum sigma^2 = 1,
    descending; Vt the axes as rows."""
    Z = np.asarray(Z, dtype=np.float64)
    X = Z - Z.mean(0)
    _, s, Vt = np.linalg.svd(X, full_matrices=False)
    sig = s / np.sqrt(len(Z))
    return sig / np.sqrt(float((sig ** 2).sum())), Vt


def normalised_sds(Z) -> np.ndarray:
    """The principal SDs of a cloud with sum of squares 1, descending."""
    return principal(Z)[0]


def _check_cloud(Z, name: str) -> np.ndarray:
    Z = np.asarray(Z, dtype=np.float64)
    if Z.ndim != 2 or Z.shape[1] != 3 or len(Z) < MIN_CLOUD or not np.isfinite(Z).all():
        raise ValueError(f"{name}: a finite cloud of at least {MIN_CLOUD} points in 3-D")
    return Z


def identifiable(Z, draws: int = IDENT_DRAWS, frac: float = IDENT_FRAC,
                 stream: int = IDENT_STREAM) -> tuple[bool, bool]:
    """Whether the axis pairs (1, 2) and (2, 3) of a bracket are identifiable at the scorer's count: in each of `draws`
    draws of max(MIN_CLOUD, round(frac N)) cells without replacement (stream default_rng([stream, h]), no build seed),
    each draw axis j is given its nearest full-stage axis, argmax_k |<v_j(draw), v_k(full)>|; a pair is stretched apart
    only if, in every draw, the draw's axes on each side of the cut are nearest full axes on the same side. It certifies
    only that no draw moves an axis across the cut, not that an axis is anatomical.

    Why sides of the cut and not "both axes matched": with a round (2, 3) plane, draw axis 2 is nearest full axis 3
    about half the time, so "both axes of (1, 2) matched" would fail too and merge all three axes (the identity),
    although axis 1 is well separated from the plane. The sides of the cut give the blocks {1}, {2, 3} the row's
    fallback describes."""
    Z = _check_cloud(Z, "bracket")
    _, V = principal(Z)
    m = min(len(Z), max(MIN_CLOUD, int(round(frac * len(Z)))))
    ok = [True, True]
    for h in range(draws):
        rows = np.random.default_rng([stream, h]).choice(len(Z), m, replace=False)
        _, v = principal(Z[rows])
        near = np.argmax(np.abs(v @ V.T), axis=1)
        for j in (0, 1):
            if not (set(near[:j + 1].tolist()) <= set(range(j + 1))
                    and set(near[j + 1:].tolist()) <= set(range(j + 1, 3))):
                ok[j] = False
    return bool(ok[0]), bool(ok[1])


def blocks_from(pairable: tuple[bool, bool]) -> tuple[tuple[int, ...], ...]:
    """The blocks of axes (0-based) a stretch acts on as one: adjacent axes whose pair is not pairable merge."""
    out = [[0]]
    for j, ok in enumerate(pairable, start=1):
        if ok:
            out.append([j])
        else:
            out[-1].append(j)
    return tuple(tuple(b) for b in out)


def bracket_blocks(Zb, Za) -> tuple[tuple[int, ...], ...]:
    """F1's blocks for a pair of brackets: a pair of adjacent axes is pairable if it is identifiable in both."""
    ib, ia = identifiable(Zb), identifiable(Za)
    return blocks_from((ib[0] and ia[0], ib[1] and ia[1]))


def blocks_text(blocks) -> str:
    """Blocks as flat text, 1-based: '1|2,3' (axes 2 and 3 move together)."""
    return "|".join(",".join(str(j + 1) for j in b) for b in blocks)


def log_interp_sds(s_b, s_a, t: float, anti: bool = False) -> np.ndarray:
    """log s* = (1 - t) log s_b + t log s_a (or, anti, (1 + t) log s_b - t log s_a), renormalised to sum of squares 1.
    Index by index: both inputs descending, as the scorer pairs them."""
    lb, la = np.log(np.asarray(s_b, np.float64)), np.log(np.asarray(s_a, np.float64))
    ls = (1.0 + t) * lb - t * la if anti else (1.0 - t) * lb + t * la
    s = np.exp(ls)
    return s / np.sqrt(float((s ** 2).sum()))


def stretch_factors(s_own, s_target, blocks) -> np.ndarray:
    """Per own axis: s_target / s_own outside blocks; inside a block one factor sqrt(sum s_target^2 / sum s_own^2),
    which keeps the cloud's own shape in that block. One block of all three axes gives factor 1 (the identity)."""
    s_own, s_target = np.asarray(s_own, np.float64), np.asarray(s_target, np.float64)
    f = np.empty(3)
    for b in blocks:
        idx = list(b)
        if len(idx) == 3:
            f[idx] = 1.0
        else:
            f[idx] = np.sqrt(float((s_target[idx] ** 2).sum()) / float((s_own[idx] ** 2).sum()))
    return f


def stretch_matrix(Vt, f) -> np.ndarray:
    """The symmetric linear map that scales the axes Vt (rows) by f: Vt^T diag(f) Vt (row vectors: x @ M)."""
    Vt = np.asarray(Vt, np.float64)
    return Vt.T @ np.diag(np.asarray(f, np.float64)) @ Vt


def cov_trace(Z) -> np.ndarray:
    """The cloud's covariance (population form) divided by its trace."""
    Z = np.asarray(Z, dtype=np.float64)
    X = Z - Z.mean(0)
    S = X.T @ X / len(Z)
    return S / float(np.trace(S))


def _sym_fn(S, fn) -> np.ndarray:
    S = np.asarray(S, np.float64)
    w, V = np.linalg.eigh((S + S.T) / 2.0)
    if not (np.isfinite(w).all() and (w > 0).all()):
        raise ShapeRefused("a covariance is not positive definite")
    return (V * fn(w)) @ V.T


def spd_log(S) -> np.ndarray:
    return _sym_fn(S, np.log)


def spd_exp(L) -> np.ndarray:
    L = np.asarray(L, np.float64)
    w, V = np.linalg.eigh((L + L.T) / 2.0)
    return (V * np.exp(w)) @ V.T


def spd_power(S, p: float) -> np.ndarray:
    return _sym_fn(S, lambda w: w ** p)


def log_euclid(S_b, S_a, t: float, anti: bool = False) -> np.ndarray:
    """Trace-normalised expm((1 - t) logm S_b + t logm S_a) (anti: (1 + t) logm S_b - t logm S_a)."""
    lb, la = spd_log(S_b), spd_log(S_a)
    L = (1.0 + t) * lb - t * la if anti else (1.0 - t) * lb + t * la
    S = spd_exp(L)
    return S / float(np.trace(S))


def brenier(S_o, S_star) -> np.ndarray:
    """The unique symmetric positive-definite A with A S_o A = S_star: S_o^-1/2 (S_o^1/2 S_star S_o^1/2)^1/2
    S_o^-1/2 (the linear optimal-transport map between centred Gaussians)."""
    r = spd_power(S_o, 0.5)
    ri = spd_power(S_o, -0.5)
    A = ri @ spd_power(r @ S_star @ r, 0.5) @ ri
    return (A + A.T) / 2.0


def check_rotation(R, name: str = "R") -> np.ndarray:
    """A finite proper rotation, or ShapeRefused."""
    R = np.asarray(R, dtype=np.float64)
    if (R.shape != (3, 3) or not np.isfinite(R).all() or np.abs(R.T @ R - np.eye(3)).max() > ROTATION_TOL
            or np.linalg.det(R) <= 0):
        raise ShapeRefused(f"{name} is not a finite proper rotation")
    return R


def random_rotation(seed: int) -> np.ndarray:
    """The controls' uniform proper rotation Q: scipy Rotation.random from default_rng([seed, ROT_STREAM])."""
    from scipy.spatial.transform import Rotation

    return check_rotation(Rotation.random(random_state=np.random.default_rng([int(seed), ROT_STREAM])).as_matrix(),
                          "Q")


def range_guard(sd_new, sd_b, sd_a) -> bool:
    """F2's arm guard: the smallest and the largest normalised principal SD of the new cloud each lie within [min,
    max] of the two brackets' corresponding SD (a relative GUARD_SLACK for rounding)."""
    ok = True
    for i in (-1, 0):
        lo, hi = min(sd_b[i], sd_a[i]), max(sd_b[i], sd_a[i])
        ok &= bool(lo * (1 - GUARD_SLACK) <= sd_new[i] <= hi * (1 + GUARD_SLACK))
    return ok


# ------------------------------------------------------------------------------------------------ the forms
def _affine(Z, m, M) -> np.ndarray:
    return m + (np.asarray(Z, np.float64) - m) @ M


def stretch_ratio(Z, Zb, Za, t: float, anti: bool = False) -> ShapeResult:
    """F1 (anti: its control). Applied to Z about Z's mean; the guard pool is Zb under the same map."""
    Z, Zb, Za = _check_cloud(Z, "Z"), _check_cloud(Zb, "Zb"), _check_cloud(Za, "Za")
    s_b, s_a = normalised_sds(Zb), normalised_sds(Za)
    blocks = bracket_blocks(Zb, Za)
    target = log_interp_sds(s_b, s_a, t, anti=anti)
    s_o, V = principal(Z)
    f = stretch_factors(s_o, target, blocks)
    M = stretch_matrix(V, f)
    m = Z.mean(0)
    Zn = _affine(Z, m, M)
    rec = {"blocks": blocks_text(blocks), "identity": bool(len(blocks) == 1),
           **{f"f{j + 1}": float(f[j]) for j in range(3)},
           **{f"target_sd{j + 1}": float(target[j]) for j in range(3)},
           "target_descending": bool(np.all(np.diff(target) <= 0)), "thin_sd": float(normalised_sds(Zn)[-1])}
    return ShapeResult(Z=Zn, pool=_affine(Zb, m, M), pool_weights=None, record=rec)


def stretch_registered(Z, Zb, Za, t: float, R, variant: str = "arm", Q=None) -> ShapeResult:
    """F2 (variant "arm"), or its controls ("anti": the mirrored covariance step, same R; "rot": the interpolated
    covariance turned by Q). Only the arm carries the bracket-range guard."""
    Z, Zb, Za = _check_cloud(Z, "Z"), _check_cloud(Zb, "Zb"), _check_cloud(Za, "Za")
    if variant not in ("arm", "anti", "rot"):
        raise ValueError(f"unknown variant {variant!r}")
    R = check_rotation(R)
    S_b, S_a = cov_trace(Zb), cov_trace(Za @ R)
    S_star = log_euclid(S_b, S_a, t, anti=variant == "anti")
    if variant == "rot":
        Q = check_rotation(Q, "Q")
        S_star = Q @ S_star @ Q.T
    A = brenier(cov_trace(Z), S_star)
    if not np.isfinite(A).all():
        raise ShapeRefused("the linear map is not finite")
    m = Z.mean(0)
    Zn = _affine(Z, m, A)
    sd = normalised_sds(Zn)
    eig = np.sort(np.linalg.eigvalsh(S_star))[::-1]
    rec = {**{f"eig{j + 1}": float(eig[j]) for j in range(3)}, "thin_sd": float(sd[-1]), "wide_sd": float(sd[0])}
    if variant == "arm":
        ok = range_guard(sd, normalised_sds(Zb), normalised_sds(Za))
        rec["guard"] = "pass" if ok else "refused"
        if not ok:
            raise ShapeRefused("stretch_registered: the cloud's smallest or largest normalised SD leaves the brackets' "
                               "range")
    return ShapeResult(Z=Zn, pool=_affine(Zb, m, A), pool_weights=None, record=rec)


def overlay(Zb, Za, rows_b, rows_a, perm, R, t_mix: float) -> ShapeResult:
    """F3 (and, with Q for R, its control): Z = vstack(Zb[rows_b], (Za @ R)[rows_a])[perm]; pool vstack(Zb, Za @ R)
    with weights (1 - t_mix) / N_b and t_mix / N_a."""
    Zb, Za = _check_cloud(Zb, "Zb"), _check_cloud(Za, "Za")
    R = check_rotation(R)
    Zr = Za @ R
    rows_b, rows_a = np.asarray(rows_b, np.int64), np.asarray(rows_a, np.int64)
    Zn = np.vstack([Zb[rows_b], Zr[rows_a]])[np.asarray(perm, np.int64)]
    w = np.r_[np.full(len(Zb), (1.0 - t_mix) / len(Zb)), np.full(len(Za), t_mix / len(Za))]
    return ShapeResult(Z=Zn, pool=np.vstack([Zb, Zr]), pool_weights=w,
                       record={"thin_sd": float(normalised_sds(Zn)[-1])})


def geodesic_overlay(Zb, Za, rows_b, rows_a, perm, R, t_mix: float, *, t: float, seed: int, Fb, Fa) -> ShapeResult:
    """Row E53 (as specified; branch campaign/moon-embryo-geodesic): overlay()'s rows, order and pool weights, with
    positions on one LDDMM geodesic (vec.lddmm) from the earlier cloud Zb to the later cloud in the earlier frame,
    Za @ R, fitted at the build seed on the brackets' positions and expression features (Fb, Fa: every cell's joint
    PCs). Each earlier row moves forward from 0 to t, each later row back from 1 to t. With Q for R (the placebo) the
    fit starts from Q's turn of the later cloud. The pool for the cloud guard is every cell of both brackets carried
    the same way. Record (flat): the fit's momentum and data facts (lddmm_*), the smallest Jacobian determinant over
    the drawn rows (jac_min), and the thinnest normalised SD of the moved union and of the two brackets."""
    from vec import lddmm

    Zb, Za = _check_cloud(Zb, "Zb"), _check_cloud(Za, "Za")
    R = check_rotation(R)
    if not 0.0 <= float(t) <= 1.0:
        raise ValueError(f"t must lie in [0, 1], got {t}")
    Fb, Fa = np.asarray(Fb, np.float64), np.asarray(Fa, np.float64)
    if Fb.ndim != 2 or Fa.ndim != 2 or len(Fb) != len(Zb) or len(Fa) != len(Za) or Fb.shape[1] != Fa.shape[1]:
        raise ValueError("features: one row per cell of each bracket, the same columns")
    Zr = Za @ R
    try:
        fit = lddmm.fit_cached(Zb, Zr, Fb, Fa, seed)
    except ValueError as e:
        raise ShapeRefused(f"geodesic_overlay: {e}") from e
    t = float(t)
    mb, ma = fit.carry(Zb, 0.0, t), fit.carry(Zr, 1.0, t)
    if not (np.isfinite(mb).all() and np.isfinite(ma).all()):
        raise ShapeRefused("geodesic_overlay: a carried position is not finite")
    rows_b, rows_a = np.asarray(rows_b, np.int64), np.asarray(rows_a, np.int64)
    Zn = np.vstack([mb[rows_b], ma[rows_a]])[np.asarray(perm, np.int64)]
    jac = []
    if len(rows_b):
        jac.append(fit.jacobian_det(Zb[np.unique(rows_b)], 0.0, t).min())
    if len(rows_a):
        jac.append(fit.jacobian_det(Zr[np.unique(rows_a)], 1.0, t).min())
    w = np.r_[np.full(len(Zb), (1.0 - t_mix) / len(Zb)), np.full(len(Za), t_mix / len(Za))]
    thin, thin_b, thin_a = (float(normalised_sds(c)[-1]) for c in (Zn, Zb, Za))
    rec = {**{f"lddmm_{k}": v for k, v in fit.record.items()}, "t": t, "jac_min": float(min(jac)),
           "thin_sd": thin, "thin_sd_before": thin_b, "thin_sd_after": thin_a,
           "thin_between": bool(min(thin_b, thin_a) <= thin <= max(thin_b, thin_a))}
    return ShapeResult(Z=Zn, pool=np.vstack([mb, ma]), pool_weights=w, record=rec)


GEODESIC_FORMS = ("geodesic_overlay", "geodesic_overlay_rot")  # row E53: harness candidates only, not in ALL_FORMS

# Row E55 (as specified by the cold challenger, round_1124_t2; branch campaign/emb-osr): F3 then F1's stretch, on the
# overlay's union. In their own tuples, outside FORMS, CONTROLS and ALL_FORMS: the builder offers overlay_stretch only,
# never a control. NEW_STREAM = 56 (alias AXES_STREAM), a stream no build step uses (47 the floor-exact draw, 48 the
# row order, 49 the identifiability draws, 50 ROT_STREAM, 51 to 52 vec.lddmm, 53 to 55 vec.t2_mmflow).
COMPOSED_FORMS = ("overlay_stretch",)
COMPOSED_CONTROLS = ("overlay_stretch_anti", "overlay_stretch_axes")
COMPOSED_CONTROL_OF = {"overlay_stretch": COMPOSED_CONTROLS}
NEW_STREAM = 56
AXES_STREAM = NEW_STREAM  # the same constant, by the name the code first used


def axes_rotation(seed: int) -> np.ndarray:
    """overlay_stretch_axes' uniform proper rotation Q: scipy Rotation.random from default_rng([seed, NEW_STREAM])."""
    from scipy.spatial.transform import Rotation

    return check_rotation(Rotation.random(random_state=np.random.default_rng([int(seed), NEW_STREAM])).as_matrix(),
                          "Q")


def overlay_stretch(Zb, Za, rows_b, rows_a, perm, R, t_mix: float, *, t: float, variant: str = "arm",
                    Q=None) -> ShapeResult:
    """Row E55 (OSR): U = overlay(...).Z, the union in the earlier bracket's normalised frame; blocks =
    bracket_blocks(Zb, Za); target = log_interp_sds(normalised_sds(Zb), normalised_sds(Za), t) (variant "anti":
    anti=True); s_o, V = principal(U); f = stretch_factors(s_o, target, blocks); M = stretch_matrix(V, f) (variant
    "axes": the arm's own f along the turned frame V @ Q, M_ax = (V Q)^T diag(f) (V Q)); Z' = mean(U) + (U - mean(U))
    @ M. The pool is the same affine map of overlay's pool, with overlay's pool weights. The caller applies the size
    rule, the jitter and the cloud guard after it, unchanged. Record (flat): the blocks, f1 to f3, g = log(f1 / f2),
    the target SDs and the thin SD of Z' (and of U, thin_sd_union)."""
    if variant not in ("arm", "anti", "axes"):
        raise ValueError(f"unknown variant {variant!r}")
    Zb, Za = _check_cloud(Zb, "Zb"), _check_cloud(Za, "Za")
    ov = overlay(Zb, Za, rows_b, rows_a, perm, R, t_mix)
    U = _check_cloud(ov.Z, "U")
    blocks = bracket_blocks(Zb, Za)
    target = log_interp_sds(normalised_sds(Zb), normalised_sds(Za), t, anti=variant == "anti")
    s_o, V = principal(U)
    f = stretch_factors(s_o, target, blocks)
    if variant == "axes":
        Q = check_rotation(Q, "Q")
        M = stretch_matrix(V @ Q, f)
    else:
        M = stretch_matrix(V, f)
    m = U.mean(0)
    Zn = _affine(U, m, M)
    rec = {"blocks": blocks_text(blocks), **{f"f{j + 1}": float(f[j]) for j in range(3)},
           "g": float(np.log(f[0] / f[1])), **{f"target_sd{j + 1}": float(target[j]) for j in range(3)},
           "thin_sd": float(normalised_sds(Zn)[-1]), "thin_sd_union": float(s_o[-1])}
    return ShapeResult(Z=Zn, pool=_affine(ov.pool, m, M), pool_weights=ov.pool_weights, record=rec)


# Row R71 (the captain's order of 2026-11-25: the embryo's consensus-form files BF on O and on SO only) and its gate row
# (experiments/coords/bf_read.py), as EMBRYO_MOON2_design.md's P1 specifies BF with the cold challenge's fixes 3, 4 and
# 7: the union's rows are carried by a smoothed optimal-transport map onto q, the equal-weight mean of the released
# specimens' canonical densities. Constants fixed before any read, none tuned. Streams 61 and 62 are used by no other
# build step (47 to 56 above and in vec.t2_build, vec.lddmm and vec.t2_mmflow).
BAYES_FORMS = ("overlay_bayes",)  # BF, buildable only under row R71
BAYES_CONTROLS = ("overlay_bayes_brackets", "overlay_bayes_noflip")  # BF2 and BFns, the gate's rule-7 controls
BAYES_CONTROL_OF = {"overlay_bayes": BAYES_CONTROLS}
BAYES_VARIANT = {"overlay_bayes": "arm", "overlay_bayes_brackets": "brackets", "overlay_bayes_noflip": "noflip"}
BAYES_PRIOR_STAGES = ("E6.75", "E7.25", "E7.5", "E8.0", "E8.25_late", "E8.75", "E9.5")  # the test's prior, in day order
BAYES_GRID_N = 48  # the density's cube over [-BAYES_EXTENT, BAYES_EXTENT] (canonical units)
BAYES_EXTENT = 3.0
BAYES_SIGMA = 0.1875  # the Gaussian smoothing's SD in canonical units (1.5 voxels of the 48-cube)
BAYES_ANCHORS = 4000  # rows of the union and draws from q in the exact assignment, at every cell count
BAYES_SMOOTH_K = 30  # each displacement is averaged over the cell's 30 nearest cells
BAYES_Y_STREAM = 61  # the draws from q
BAYES_X_STREAM = 62  # the union's drawn rows
NB_SHARE_K = 15  # G0's neighbour share: each cell's 15 nearest cells
OCC_GRID_N = 16  # the scorer's occupancy_dice grid (vec.t2_build.GRID_N)
OCC_EXTENT = 3.0
OCC_EDGE_EPS = 1e-9  # the scorer's clip below the upper edge
AXIS_SIGNS = tuple((a, b, c) for a in (1, -1) for b in (1, -1) for c in (1, -1))  # eight sign copies, identity first
PROPER_FLIPS = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))  # the scorer's _PROPER_FLIPS (det +1), in its order


@dataclass(frozen=True)
class BayesPrior:
    """The released specimens a BF form pools, each as its canonical cloud (`canonical_frame`), in BAYES_PRIOR_STAGES
    order, equal weight. Built by `bayes_prior` from coordinates the caller loaded (the harness's read guard keeps
    every stage file from a candidate; the caller leaves the task's truth out, J2)."""

    names: tuple[str, ...]
    clouds: tuple[np.ndarray, ...]


def canonical_frame(C) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """veckit main's `_canonicalise` (centre, PCA basis from np.linalg.svd of the centred cloud made right-handed,
    divide by the RMS radius), with the frame kept for the way back: (Z, mu, Vt, rms) with Z = (C - mu) @ Vt.T / rms,
    so C = mu + rms * Z @ Vt (`from_canonical`)."""
    C = _check_cloud(C, "cloud")
    mu = C.mean(0)
    X = C - mu
    rms = float(np.sqrt((X ** 2).sum(1).mean()))
    if not rms > 1e-12:
        raise ShapeRefused("a cloud with no spread has no canonical frame")
    _, _, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    if np.linalg.det(Vt) < 0:
        Vt[-1] *= -1
    return (X @ Vt.T) / rms, mu, Vt, rms


def from_canonical(Z, mu, Vt, rms: float) -> np.ndarray:
    """The inverse of `canonical_frame`'s map."""
    return np.asarray(mu, np.float64) + float(rms) * (np.asarray(Z, np.float64) @ np.asarray(Vt, np.float64))


def bayes_prior(coords: dict) -> BayesPrior:
    """A BayesPrior from {stage: coordinates (rows, 3)}: every key one of BAYES_PRIOR_STAGES, at least one, each
    canonicalised. The order is BAYES_PRIOR_STAGES', whatever the dict's."""
    names = [str(k) for k in coords]
    unknown = sorted(set(names) - set(BAYES_PRIOR_STAGES))
    if unknown or not names or len(set(names)) != len(names):
        raise ValueError(f"a BF prior takes stages of {BAYES_PRIOR_STAGES} once each, got {names}")
    order = tuple(s for s in BAYES_PRIOR_STAGES if s in coords)
    return BayesPrior(names=order, clouds=tuple(canonical_frame(np.asarray(coords[s], np.float64)[:, :3])[0]
                                                for s in order))


def _cells(Z, n_grid: int, extent: float) -> np.ndarray:
    """Integer cell indices of points on an n_grid cube over [-extent, extent], clipped into it as the scorer's
    occupancy_dice bins them."""
    return np.floor((np.clip(np.asarray(Z, np.float64), -extent, extent - OCC_EDGE_EPS) + extent) / (2 * extent)
                    * n_grid).astype(np.int64)


def specimen_density(Z, signs=AXIS_SIGNS) -> np.ndarray:
    """One canonical specimen's density on the 48-cube: the cloud stacked in each axis-sign copy of `signs` (all eight:
    flip-symmetric; a single sign triple: no flip copies), binned as the scorer bins (clipped into the cube),
    Gaussian-smoothed at SD BAYES_SIGMA (scipy.ndimage.gaussian_filter, zero outside the cube) and normalised to mass
    1."""
    from scipy.ndimage import gaussian_filter

    Z = _check_cloud(Z, "specimen")
    G = BAYES_GRID_N
    H = np.zeros(G ** 3)
    for s in signs:
        v = _cells(Z * np.asarray(s, np.float64), G, BAYES_EXTENT)
        H += np.bincount((v[:, 0] * G + v[:, 1]) * G + v[:, 2], minlength=G ** 3)
    H = gaussian_filter(H.reshape(G, G, G), sigma=BAYES_SIGMA / (2 * BAYES_EXTENT / G), mode="constant", cval=0.0)
    H = np.clip(H, 0.0, None)
    return H / float(H.sum())


def occupancy(Z) -> np.ndarray:
    """Occupied voxels of canonical coordinates on the scorer's 16-cube over [-3, 3], as sorted unique ids."""
    v = _cells(Z, OCC_GRID_N, OCC_EXTENT)
    return np.unique((v[:, 0] * OCC_GRID_N + v[:, 1]) * OCC_GRID_N + v[:, 2])


def dice(a, b) -> float:
    """2 |a & b| / (|a| + |b|) of two sorted unique id arrays."""
    a, b = np.asarray(a), np.asarray(b)
    return 2.0 * len(np.intersect1d(a, b, assume_unique=True)) / max(len(a) + len(b), 1)


def best_proper_flip(S, X) -> tuple[tuple[int, int, int], float]:
    """The proper flip f of the scorer's four (first best in their order) that maximises the 16-cube occupancy Dice of
    S * f against X (two canonical clouds, all their rows), and that Dice."""
    ox = occupancy(X)
    best, flip = -1.0, PROPER_FLIPS[0]
    for f in PROPER_FLIPS:
        d = dice(occupancy(np.asarray(S, np.float64) * np.asarray(f, np.float64)), ox)
        if d > best:
            best, flip = d, f
    return flip, float(best)


def consensus_density(prior: BayesPrior, variant: str = "arm", X=None) -> tuple[np.ndarray, tuple]:
    """q, the equal-weight mean of the prior's specimen densities: with their eight flip copies ("arm", "brackets"), or
    ("noflip") each specimen without flip copies, turned onto X's canonical frame by its best proper flip. Returns
    (q, the flips used: () unless "noflip")."""
    if variant == "noflip":
        if X is None:
            raise ValueError("the no-flip density needs the union's canonical cloud X")
        flips = tuple(best_proper_flip(S, X)[0] for S in prior.clouds)
        dens = [specimen_density(np.asarray(S) * np.asarray(f, np.float64), signs=((1, 1, 1),))
                for S, f in zip(prior.clouds, flips, strict=True)]
        return np.mean(dens, axis=0), flips
    return np.mean([specimen_density(S) for S in prior.clouds], axis=0), ()


def density_draws(q, k: int, rng: np.random.Generator) -> np.ndarray:
    """k points from q: a cell with probability its mass (a multinomial draw over the cells), then uniform within it."""
    q = np.asarray(q, np.float64)
    G = q.shape[0]
    edge = 2 * BAYES_EXTENT / G
    cells = rng.choice(q.size, size=int(k), p=q.ravel())
    ijk = np.stack(np.unravel_index(cells, q.shape), axis=1).astype(np.float64)
    return -BAYES_EXTENT + (ijk + rng.random((int(k), 3))) * edge


def transport_displacements(X, q, seed: int) -> np.ndarray:
    """Step 4 as the challenge's fix 3 sets it, with no branch on the cell count: BAYES_ANCHORS rows of X drawn with
    replacement (stream [seed, BAYES_X_STREAM]) and BAYES_ANCHORS points drawn from q (stream [seed, BAYES_Y_STREAM]);
    the exact assignment on squared distance (scipy linear_sum_assignment on float64 costs summed coordinate by
    coordinate, no BLAS); a drawn row's displacement is the mean over its draws of partner minus row; every row of X
    takes the displacement of its nearest drawn row (cKDTree). Returns the raw displacements, one per row of X."""
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial import cKDTree

    X = _check_cloud(X, "X")
    K = BAYES_ANCHORS
    Y = density_draws(q, K, np.random.default_rng([int(seed), BAYES_Y_STREAM]))
    rows = np.random.default_rng([int(seed), BAYES_X_STREAM]).integers(0, len(X), K)
    A = X[rows]
    cost = (A[:, 0, None] - Y[None, :, 0]) ** 2
    cost += (A[:, 1, None] - Y[None, :, 1]) ** 2
    cost += (A[:, 2, None] - Y[None, :, 2]) ** 2
    r, c = linear_sum_assignment(cost)
    d = Y[c] - A[r]
    drawn, inv = np.unique(rows[r], return_inverse=True)
    D = np.zeros((len(drawn), 3))
    np.add.at(D, inv, d)
    D /= np.bincount(inv, minlength=len(drawn))[:, None]
    near = cKDTree(X[drawn]).query(X, k=1)[1]
    return D[near]


def smoothed(X, d, at=None) -> np.ndarray:
    """The Gaussian average of X's displacements d over the BAYES_SMOOTH_K nearest rows of X: at X's own rows (`at`
    None; each row's own first) or at the points `at` (the pool), with weights exp(-r^2 / (2 h^2)), h the median over
    X's rows of the distance to their BAYES_SMOOTH_K-th nearest row."""
    from scipy.spatial import cKDTree

    X = _check_cloud(X, "X")
    if len(X) < BAYES_SMOOTH_K:
        raise ValueError(f"the smoothing needs at least {BAYES_SMOOTH_K} cells")
    tree = cKDTree(X)
    r, idx = tree.query(X, k=BAYES_SMOOTH_K)
    h = float(np.median(r[:, -1]))
    if not h > 0:
        raise ShapeRefused("the 30-nearest-cell reach is zero: no smoothing scale")
    if at is not None:
        r, idx = tree.query(np.asarray(at, np.float64), k=BAYES_SMOOTH_K)
    w = np.exp(-0.5 * (r / h) ** 2)
    return (w[..., None] * np.asarray(d, np.float64)[idx]).sum(1) / w.sum(1)[:, None]


def knn_sets(P, k: int = NB_SHARE_K) -> np.ndarray:
    """Each row's k nearest rows of P, itself left out (a coincident row may stand first; then the k + 1-th drops)."""
    from scipy.spatial import cKDTree

    P = np.asarray(P, np.float64)
    idx = cKDTree(P).query(P, k=k + 1)[1]
    own = idx == np.arange(len(P))[:, None]
    return np.take_along_axis(idx, np.argsort(own, axis=1, kind="stable"), axis=1)[:, :k]


def neighbour_share(P, P2, k: int = NB_SHARE_K) -> float:
    """G0's statistic: the median over rows of the share of each row's k nearest rows in P that stay among its k
    nearest in P2 (the same rows, moved)."""
    A, B = knn_sets(P, k), knn_sets(P2, k)
    return float(np.median((A[:, :, None] == B[:, None, :]).any(2).sum(1) / k))


def overlay_bayes(Zb, Za, rows_b, rows_a, perm, R, t_mix: float, *, seed: int, prior: BayesPrior | None = None,
                  variant: str = "arm") -> ShapeResult:
    """Row R71's BF (variant "arm") and the gate's controls ("brackets": the prior is the task's two brackets, Zb and
    Za, with flip copies; "noflip": the prior's specimens without flip copies, each turned onto the union's frame by its
    best proper flip). U = overlay(...).Z; X = U canonicalised (`canonical_frame`); q = `consensus_density`; raw
    displacements by `transport_displacements`; X' = X + `smoothed`(X, d); Z' = X' taken back to U's frame. The pool is
    overlay's pool moved by the same smoothed field (each pool cell takes the Gaussian 30-nearest-cell average
    displacement of X's rows), with overlay's pool weights (the challenge's fix 7). The caller applies the size rule,
    the jitter and the cloud guard after it, unchanged. Rows, order, count and expression are overlay's. Record (flat,
    no count): the variant, the prior's stages, the flips (noflip), the constants and the thin SDs of Z' and U."""
    if variant not in ("arm", "brackets", "noflip"):
        raise ValueError(f"unknown variant {variant!r}")
    if variant == "brackets":
        if prior is not None:
            raise ValueError("overlay_bayes_brackets takes no prior: its prior is the task's two brackets")
        prior = BayesPrior(names=("before", "after"),
                           clouds=(canonical_frame(Zb)[0], canonical_frame(Za)[0]))
    elif not isinstance(prior, BayesPrior) or not prior.clouds or len(prior.clouds) != len(prior.names):
        raise ValueError(f"overlay_bayes ({variant}) needs a BayesPrior")
    ov = overlay(Zb, Za, rows_b, rows_a, perm, R, t_mix)
    U = _check_cloud(ov.Z, "U")
    X, mu, Vt, rms = canonical_frame(U)
    q, flips = consensus_density(prior, "noflip" if variant == "noflip" else "arm", X)
    d = transport_displacements(X, q, seed)
    Xn = X + smoothed(X, d)
    P = (np.asarray(ov.pool, np.float64) - mu) @ Vt.T / rms
    Pn = P + smoothed(X, d, at=P)
    Zn = from_canonical(Xn, mu, Vt, rms)
    if not (np.isfinite(Zn).all() and np.isfinite(Pn).all()):
        raise ShapeRefused("overlay_bayes: a moved position is not finite")
    rec = {"variant": variant, "prior": "|".join(prior.names), "flip_copies": 1 if variant == "noflip" else 8,
           "anchors": BAYES_ANCHORS, "y_stream": BAYES_Y_STREAM, "x_stream": BAYES_X_STREAM, "grid": BAYES_GRID_N,
           "sigma": BAYES_SIGMA, "smooth_k": BAYES_SMOOTH_K, "thin_sd": float(normalised_sds(Zn)[-1]),
           "thin_sd_union": float(normalised_sds(U)[-1])}
    if variant == "noflip":
        rec["flips"] = "|".join(",".join(str(int(v)) for v in f) for f in flips)
    return ShapeResult(Z=Zn, pool=from_canonical(Pn, mu, Vt, rms), pool_weights=ov.pool_weights, record=rec)


def apply_form(form: str, *, Z, Zb, Za, t: float, R, seed: int, rows_b, rows_a, perm, t_mix: float,
               features=None, prior: BayesPrior | None = None) -> ShapeResult:
    """One named form or control on the exact-position cloud Z (Zb, Za: the brackets in the earlier frame before R;
    R: the build's coupling rotation; seed: the build seed, used by the controls' Q only). features (row E53's
    geodesic forms only): (Fb, Fa), every cell's joint PCs of each bracket. prior (row R71's overlay_bayes and its
    no-flip control only): the BayesPrior of released specimens the caller loaded."""
    if form in BAYES_FORMS or form in BAYES_CONTROLS:  # row R71: dispatched before every other branch
        if features is not None:
            raise ValueError(f"{form} takes no features")
        out = overlay_bayes(Zb, Za, rows_b, rows_a, perm, R, t_mix, seed=seed, prior=prior,
                            variant=BAYES_VARIANT[form])
        out.record = {"form": form, **out.record}
        return out
    if prior is not None:
        raise ValueError(f"{form} takes no prior")
    if form in GEODESIC_FORMS:
        if features is None:
            raise ValueError(f"{form} needs features=(Fb, Fa)")
        Q = R if form == "geodesic_overlay" else random_rotation(seed)
        out = geodesic_overlay(Zb, Za, rows_b, rows_a, perm, Q, t_mix, t=t, seed=seed, Fb=features[0],
                               Fa=features[1])
        out.record = {"form": form, **out.record}
        if form.endswith("_rot"):
            out.record["rot_stream"] = ROT_STREAM
        return out
    if features is not None:
        raise ValueError(f"{form} takes no features")
    if form in COMPOSED_FORMS or form in COMPOSED_CONTROLS:  # row E55: dispatched before the FORMS branch
        variant = {"overlay_stretch": "arm", "overlay_stretch_anti": "anti", "overlay_stretch_axes": "axes"}[form]
        out = overlay_stretch(Zb, Za, rows_b, rows_a, perm, R, t_mix, t=t, variant=variant,
                              Q=axes_rotation(seed) if variant == "axes" else None)
        out.record = {"form": form, **out.record, "ident_stream": IDENT_STREAM}
        if variant == "axes":
            out.record["axes_stream"] = NEW_STREAM
        return out
    if form == "stretch_ratio":
        out = stretch_ratio(Z, Zb, Za, t)
    elif form == "stretch_ratio_anti":
        out = stretch_ratio(Z, Zb, Za, t, anti=True)
    elif form == "stretch_registered":
        out = stretch_registered(Z, Zb, Za, t, R)
    elif form == "stretch_registered_anti":
        out = stretch_registered(Z, Zb, Za, t, R, variant="anti")
    elif form == "stretch_registered_rot":
        out = stretch_registered(Z, Zb, Za, t, R, variant="rot", Q=random_rotation(seed))
    elif form == "overlay":
        out = overlay(Zb, Za, rows_b, rows_a, perm, R, t_mix)
    elif form == "overlay_rot":
        out = overlay(Zb, Za, rows_b, rows_a, perm, random_rotation(seed), t_mix)
    else:
        raise ValueError(f"unknown shape form {form!r}; known: {ALL_FORMS}")
    out.record = {"form": form, **out.record}
    if form.endswith("_rot"):
        out.record["rot_stream"] = ROT_STREAM
    if form.startswith("stretch_ratio"):
        out.record["ident_stream"] = IDENT_STREAM
    return out
