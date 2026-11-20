"""The shape forms of row E49 (vec.shape_interp) and their hook in vec.t2_build.interpolate_clouds, on synthetic clouds
only. Nothing is read from data/ and nothing is scored.

Pinned: F1 reaches a hand-computed profile (an ellipsoid stretched to known ratios) and keeps the cloud's own ratio
inside a block when two variances are too close to identify; F2's map is the symmetric positive-definite map that gives
the interpolated covariance, and its guard refuses a cloud outside the brackets' range; F3 puts each later row at its
own rotated position; the controls are what the row says; every form keeps rows, order, count and the RMS radius the
size rule sets (to machine precision) and introduces no mirror image; the default call is bit-identical."""

from __future__ import annotations

import numpy as np
import pytest
from test_occupancy_frame import _stage

from vec import evaluate as ev
from vec import harness
from vec import shape_interp as si
from vec import t2_build as tb


def _rot(seed: int) -> np.ndarray:
    from scipy.spatial.transform import Rotation

    return Rotation.random(random_state=np.random.default_rng(seed)).as_matrix()


def _ellipsoid(n: int, sds, seed: int, R=None) -> np.ndarray:
    """A centred Gaussian cloud whose principal SDs are exactly `sds` along the rows of R (whitened, then scaled)."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 3))
    X -= X.mean(0)
    w, V = np.linalg.eigh(X.T @ X / n)
    X = X @ V @ np.diag(w ** -0.5) @ V.T  # exactly white
    X = X * np.asarray(sds, float)
    return X @ (np.eye(3) if R is None else R)


# ------------------------------------------------------------------------------------------------ constants
def test_the_constants_are_the_scorers_and_the_builders():
    assert si.IDENT_FRAC == ev.SUBSAMPLE_FRAC / 2 and si.MIN_CLOUD == tb.MIN_CLOUD
    assert si.ROT_STREAM == 50 and si.ROT_STREAM not in (tb.COPY_PARITY_STREAM, tb.FLOOR_EXACT_STREAM)
    assert si.IDENT_STREAM == 49 and si.IDENT_DRAWS == 20
    assert si.FORMS == ("stretch_ratio", "stretch_registered", "overlay")
    assert set(si.CONTROLS) == {c for cs in si.CONTROL_OF.values() for c in cs}
    assert set(harness.SHAPE_CANDIDATES.values()) == set(si.ALL_FORMS)


# ------------------------------------------------------------------------------------------------ F1
def test_log_interpolation_of_the_profile_by_hand():
    s_b = np.array([3.0, 2.0, 1.0]) / np.sqrt(14.0)
    s_a = np.array([4.0, 2.0, 0.5]) / np.sqrt(20.25)
    want = np.sqrt(np.array([12.0, 4.0, 0.5]) / 16.5)  # geometric means of the raw SDs, renormalised
    np.testing.assert_allclose(si.log_interp_sds(s_b, s_a, 0.5), want, rtol=0, atol=1e-15)
    # the mirrored step about the earlier bracket: s_b^1.5 / s_a^0.5, renormalised
    raw = np.array([3.0, 2.0, 1.0]) ** 1.5 / np.array([4.0, 2.0, 0.5]) ** 0.5
    np.testing.assert_allclose(si.log_interp_sds(s_b, s_a, 0.5, anti=True), raw / np.linalg.norm(raw), atol=1e-15)


def test_f1_stretches_an_ellipsoid_to_the_hand_computed_profile():
    Zb = _ellipsoid(20000, [3.0, 2.0, 1.0], 1, _rot(1))
    Za = _ellipsoid(20000, [4.0, 2.0, 0.5], 2, _rot(2))
    Z = _ellipsoid(3000, [3.0, 2.0, 1.0], 3, _rot(1))  # a sample shaped as the earlier bracket
    assert si.bracket_blocks(Zb, Za) == ((0,), (1,), (2,))
    out = si.stretch_ratio(Z, Zb, Za, 0.5)
    want = np.sqrt(np.array([12.0, 4.0, 0.5]) / 16.5)
    np.testing.assert_allclose(si.normalised_sds(out.Z), want, atol=1e-12)
    assert out.record["blocks"] == "1|2|3" and not out.record["identity"]
    # the map is symmetric with positive factors: no mirror image, the mean is kept
    M = np.linalg.lstsq(Z - Z.mean(0), out.Z - out.Z.mean(0), rcond=None)[0]
    assert np.linalg.det(M) > 0 and np.allclose(M, M.T, atol=1e-10)
    np.testing.assert_allclose(out.Z.mean(0), Z.mean(0), atol=1e-12)
    # the guard pool is the earlier bracket under the same map
    np.testing.assert_allclose(out.pool - Z.mean(0), (Zb - Z.mean(0)) @ M, atol=1e-9)


def test_f1_falls_back_inside_a_plane_whose_two_variances_are_too_close():
    Zb = _ellipsoid(20000, [3.0, 1.0, 1.0], 4, _rot(4))  # a round 2-3 plane: no axis in it is identifiable
    Za = _ellipsoid(20000, [4.0, 2.0, 0.5], 5, _rot(5))
    assert si.identifiable(Zb) == (True, False) and si.identifiable(Za) == (True, True)
    assert si.bracket_blocks(Zb, Za) == ((0,), (1, 2))
    Z = Zb[:3000]
    s0 = si.normalised_sds(Z)
    out = si.stretch_ratio(Z, Zb, Za, 0.5)
    s1 = si.normalised_sds(out.Z)
    assert s1[1] / s1[2] == pytest.approx(s0[1] / s0[2], rel=1e-12)  # the plane keeps its own shape
    assert out.record["f2"] == out.record["f3"] and out.record["blocks"] == "1|2,3"
    target = si.log_interp_sds(si.normalised_sds(Zb), si.normalised_sds(Za), 0.5)
    assert s1[0] == pytest.approx(target[0], rel=1e-12)  # axis 1 reaches its target; the plane takes the rest


def test_f1_is_the_identity_when_no_axis_is_identifiable():
    Zb = _ellipsoid(20000, [1.0, 1.0, 1.0], 6)
    Za = _ellipsoid(20000, [4.0, 2.0, 0.5], 7, _rot(7))
    assert si.bracket_blocks(Zb, Za) == ((0, 1, 2),)
    Z = Zb[:2000]
    out = si.stretch_ratio(Z, Zb, Za, 0.5)
    assert out.record["identity"] and out.record["f1"] == out.record["f2"] == out.record["f3"] == 1.0
    np.testing.assert_allclose(out.Z, Z, atol=1e-13)


def test_a_round_plane_still_lets_its_well_separated_axis_stretch():
    """The sides-of-the-cut criterion: a round (2, 3) plane merges axes 2 and 3 only; "both axes matched" would merge
    all."""
    Zb = _ellipsoid(20000, [3.0, 1.0, 1.0], 9, _rot(9))
    ident = si.identifiable(Zb)
    assert ident == (True, False)


def _identifiable_with_draw_axes(monkeypatch, Vt_draw) -> tuple[bool, bool]:
    """identifiable() on a stub cloud whose full frame is the identity and whose every draw has the axes Vt_draw."""
    calls = {"n": 0}

    def fake_principal(Z):
        calls["n"] += 1
        return np.full(3, 3 ** -0.5), (np.eye(3) if calls["n"] == 1 else np.asarray(Vt_draw, float))

    monkeypatch.setattr(si, "principal", fake_principal)
    Z = np.random.default_rng(0).normal(size=(max(si.MIN_CLOUD, 100), 3))
    return si.identifiable(Z)


def test_a_draw_exchanging_the_first_and_third_axes_fails_both_pairs(monkeypatch):
    """near = [2, 1, 0]: no adjacent index pair is swapped, yet axes cross both cuts, so both pairs fail."""
    assert _identifiable_with_draw_axes(monkeypatch, np.eye(3)[[2, 1, 0]]) == (False, False)


@pytest.mark.parametrize("perm, want", [([0, 1, 2], (True, True)), ([1, 0, 2], (False, True)),
                                        ([0, 2, 1], (True, False)), ([1, 2, 0], (False, False)),
                                        ([2, 0, 1], (False, False)), ([0, 2, 2], (True, False)),
                                        ([1, 1, 2], (False, True))])
def test_the_sides_of_the_cut_by_hand(monkeypatch, perm, want):
    Vt = np.eye(3)[perm] + 1e-3 * np.eye(3)  # a repeated nearest axis: the first of the pair stays nearest
    assert _identifiable_with_draw_axes(monkeypatch, Vt) == want


def test_the_hand_brackets_still_give_their_blocks():
    assert si.blocks_text(si.blocks_from(si.identifiable(_ellipsoid(20000, [3.0, 1.0, 1.0], 4, _rot(4))))) == "1|2,3"
    assert si.blocks_text(si.blocks_from(si.identifiable(_ellipsoid(20000, [1.0, 1.0, 1.0], 6)))) == "1,2,3"
    assert si.blocks_text(si.blocks_from(si.identifiable(_ellipsoid(20000, [3.0, 2.0, 1.0], 1, _rot(1))))) == "1|2|3"


def test_the_identifiability_draws_are_seed_free():
    Zb = _ellipsoid(5000, [3.0, 1.05, 1.0], 8, _rot(8))
    assert si.identifiable(Zb) == si.identifiable(Zb.copy())
    assert si.blocks_text(((0,), (1, 2))) == "1|2,3" and si.blocks_from((False, False)) == ((0, 1, 2),)


# ------------------------------------------------------------------------------------------------ F2
def test_f2_gives_the_interpolated_covariance_by_a_symmetric_positive_map():
    Zb = _ellipsoid(20000, [3.0, 2.0, 1.0], 11, _rot(11))
    Za = _ellipsoid(20000, [3.5, 2.0, 0.8], 12, _rot(12))
    R = _rot(12).T @ _rot(11)  # registers the later bracket onto the earlier one's axes
    Z = Zb[:4000]
    out = si.stretch_registered(Z, Zb, Za, 0.5, R)
    S_star = si.log_euclid(si.cov_trace(Zb), si.cov_trace(Za @ R), 0.5)
    np.testing.assert_allclose(si.cov_trace(out.Z), S_star, atol=1e-12)
    A = si.brenier(si.cov_trace(Z), S_star)
    np.testing.assert_allclose(A, A.T, atol=1e-14)
    assert (np.linalg.eigvalsh(A) > 0).all()
    np.testing.assert_allclose(A @ si.cov_trace(Z) @ A, S_star, atol=1e-12)
    assert out.record["guard"] == "pass"
    # commuting covariances: the log-Euclidean mean is the element-wise geometric mean
    Sb, Sa = np.diag([9.0, 4.0, 1.0]), np.diag([16.0, 4.0, 0.25])
    np.testing.assert_allclose(si.log_euclid(Sb, Sa, 0.5), np.diag([12.0, 4.0, 0.5]) / 16.5, atol=1e-14)


def test_f2s_guard_refuses_a_cloud_outside_the_brackets_range():
    # two brackets with their thin axes crossed: the mean is rounder than both, so its smallest SD leaves their range
    Zb = _ellipsoid(20000, [3.0, 2.0, 1.0], 14)
    Za = _ellipsoid(20000, [1.0, 2.0, 3.0], 15)
    with pytest.raises(si.ShapeRefused, match="leaves the brackets' range"):
        si.stretch_registered(Zb[:3000], Zb, Za, 0.5, np.eye(3))
    # the controls carry no bracket-range guard
    assert si.stretch_registered(Zb[:3000], Zb, Za, 0.5, np.eye(3), variant="anti").Z.shape == (3000, 3)
    # flatter than both brackets is refused by the guard itself
    sd_b, sd_a = np.sqrt([0.6, 0.3, 0.1]), np.sqrt([0.62, 0.3, 0.08])
    assert si.range_guard(np.sqrt([0.61, 0.3, 0.09]), sd_b, sd_a)
    assert not si.range_guard(np.sqrt([0.61, 0.33, 0.06]), sd_b, sd_a)
    assert not si.range_guard(np.sqrt([0.7, 0.21, 0.09]), sd_b, sd_a)


def test_f2s_controls_keep_its_size_and_change_its_direction():
    Zb = _ellipsoid(20000, [3.0, 2.0, 1.0], 16, _rot(16))
    Za = _ellipsoid(20000, [3.5, 2.0, 0.8], 17, _rot(17))
    R, Z = _rot(17).T @ _rot(16), Zb[:3000]
    arm = si.stretch_registered(Z, Zb, Za, 0.5, R)
    Q = si.random_rotation(1_000_001)
    rot = si.stretch_registered(Z, Zb, Za, 0.5, R, variant="rot", Q=Q)
    np.testing.assert_allclose(si.normalised_sds(rot.Z), si.normalised_sds(arm.Z), atol=1e-12)  # same spectrum
    assert not np.allclose(rot.Z, arm.Z)
    anti = si.stretch_registered(Z, Zb, Za, 0.5, R, variant="anti")
    want = si.log_euclid(si.cov_trace(Zb), si.cov_trace(Za @ R), 0.5, anti=True)
    np.testing.assert_allclose(si.cov_trace(anti.Z), want, atol=1e-12)
    np.testing.assert_array_equal(Q, si.random_rotation(1_000_001))  # stream [seed, 50]: reproducible
    assert np.linalg.det(Q) == pytest.approx(1.0)
    with pytest.raises(si.ShapeRefused):
        si.stretch_registered(Z, Zb, Za, 0.5, np.diag([1.0, 1.0, -1.0]))  # a mirror image is refused
    with pytest.raises(si.ShapeRefused):
        si.stretch_registered(Z, Zb, Za, 0.5, np.full((3, 3), np.nan))


# ------------------------------------------------------------------------------------------------ F3
def test_f3_puts_each_later_row_at_its_own_rotated_position():
    Zb = _ellipsoid(500, [3.0, 2.0, 1.0], 21)
    Za = _ellipsoid(300, [3.5, 2.0, 0.8], 22)
    R = _rot(23)
    rows_b, rows_a = np.arange(0, 100), np.arange(50, 150)
    perm = np.random.default_rng(0).permutation(200)
    out = si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.5)
    stage = np.r_[np.zeros(100, int), np.ones(100, int)][perm]
    rows = np.r_[rows_b, rows_a][perm]
    np.testing.assert_array_equal(out.Z[stage == 0], Zb[rows[stage == 0]])
    np.testing.assert_allclose(out.Z[stage == 1], (Za @ R)[rows[stage == 1]], atol=0)
    assert out.pool.shape == (800, 3) and out.pool_weights.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(out.pool_weights[:500], 0.5 / 500)


# ------------------------------------------------------------------------------------------------ the hook
GENES_N, N_FILE, SEED, TARGET = 12, 400, 1_000_001, 6.0


@pytest.fixture(scope="module")
def stages():
    return _stage(1500, 1, 0.0), _stage(900, 2, 0.3)


def _interp(stages, **kw):
    b, a = stages
    return tb.interpolate_clouds(np.asarray(b.obsm["spatial_3D"], float), np.asarray(a.obsm["spatial_3D"], float),
                                 b.X.toarray(), a.X.toarray(), 0.5, N_FILE, SEED, target_rms=TARGET, t_pos=0.0,
                                 positions="floor_exact", n_anchor=150, **kw)


@pytest.fixture(scope="module")
def built(stages):
    out = {None: _interp(stages, with_coupling=True)}
    for f in si.ALL_FORMS:
        out[f] = _interp(stages, shape_form=f, with_coupling=True)
    return out


def test_the_default_call_is_bit_identical(stages, built):
    plain = _interp(stages, with_coupling=True)
    none = _interp(stages, shape_form=None, with_coupling=True)
    for r in (plain, none):
        np.testing.assert_array_equal(r.coords, built[None].coords)
        np.testing.assert_array_equal(r.pool_coords, built[None].pool_coords)
        assert r.record == built[None].record
    assert not any(k.startswith("shape_interp") for k in built[None].record)


@pytest.mark.parametrize("form", si.ALL_FORMS)
def test_every_form_keeps_rows_order_count_and_the_size_rule(built, form):
    base, r = built[None], built[form]
    np.testing.assert_array_equal(r.rows, base.rows)
    np.testing.assert_array_equal(r.stage, base.stage)
    assert r.record["n_jittered"] == 0 == base.record["n_jittered"]
    assert tb.rms_radius(r.coords) == pytest.approx(TARGET, rel=1e-13)  # the size rule's RMS, set after the form
    assert tb.rms_radius(r.coords) == pytest.approx(tb.rms_radius(base.coords), rel=1e-13)
    np.testing.assert_allclose(r.coords.mean(0), base.coords.mean(0), atol=1e-12)  # placed at the same centroid
    assert not np.allclose(r.coords, base.coords)
    assert r.record["shape_interp"] == form
    assert r.record["floor_exact"] == base.record["floor_exact"]  # stream 47's draw ran as in the default
    assert r.pool_weights.sum() == pytest.approx(1.0)


@pytest.mark.parametrize("form", ["stretch_ratio", "stretch_registered", "stretch_ratio_anti",
                                  "stretch_registered_anti", "stretch_registered_rot"])
def test_the_linear_forms_add_no_mirror_image(built, form):
    base, r = built[None], built[form]
    X0, X1 = base.coords - base.coords.mean(0), r.coords - r.coords.mean(0)
    M = np.linalg.lstsq(X0, X1, rcond=None)[0]
    np.testing.assert_allclose(X0 @ M, X1, atol=1e-9)  # one linear map
    assert np.linalg.det(M) > 0


@pytest.mark.parametrize("form", ["overlay", "overlay_rot"])
def test_the_overlay_sits_on_the_rotated_later_brackets_own_positions(stages, built, form):
    r = built[form]
    Zb, Za_r = r.coupling["before_frame"], r.coupling["after_frame"]  # Za_r = Za @ R
    R = r.coupling["rotation"]
    if form == "overlay_rot":
        Za_r = Za_r @ R.T @ si.random_rotation(SEED)
    Z = np.empty((len(r.rows), 3))
    Z[r.stage == 0], Z[r.stage == 1] = Zb[r.rows[r.stage == 0]], Za_r[r.rows[r.stage == 1]]
    scale = r.record["scale_from_canonical"]
    np.testing.assert_allclose(r.coords - r.coords.mean(0), (Z - Z.mean(0)) * scale, atol=1e-10)


def test_the_hook_refuses_outside_the_exact_positions(stages):
    b, a = stages
    with pytest.raises(ValueError, match="floor_exact"):
        tb.interpolate_clouds(np.asarray(b.obsm["spatial_3D"], float), np.asarray(a.obsm["spatial_3D"], float),
                              b.X.toarray(), a.X.toarray(), 0.5, N_FILE, SEED, target_rms=TARGET, n_anchor=150,
                              shape_form="stretch_ratio")
    with pytest.raises(ValueError, match="shape_form must be one of"):
        _interp(stages, shape_form="squash")


def test_the_record_holds_flat_values_and_no_count(built):
    for f in si.ALL_FORMS:
        rec = {k: v for k, v in built[f].record.items() if k.startswith("shape_interp")}
        assert rec and all(not isinstance(v, (dict, list)) for v in rec.values())
        assert not any(v == N_FILE for v in rec.values() if isinstance(v, int) and not isinstance(v, bool))
