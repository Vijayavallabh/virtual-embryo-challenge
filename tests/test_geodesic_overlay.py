"""The geodesic overlay forms of vec.shape_interp (row E53, as specified) through apply_form, interpolate_clouds and the
harness candidates, on synthetic arrays only (CPU). Nothing is read from data/ and nothing is scored.

Pinned: E49's forms, ALL_FORMS, FORMS and SHAPE_CANDIDATES are untouched (the builder offers neither new form); with
zero momenta the geodesic overlay is overlay() exactly; earlier rows go forward from 0 to t and later rows back from 1
to t on one fit of Zb against Za @ R (the arm) or Za @ Q, Q from stream [seed, 50] (the placebo); the pool is both
brackets carried alike, with overlay's weights; the record is flat and holds the gate's quantities; through
interpolate_clouds the rows, order, count and expression are fgw_shape0x's and only positions move; a real (small) fit
runs end to end through the harness candidates."""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from vec import harness, lddmm  # noqa: E402
from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402


def _pair(seed=0, nb=120, na=150):
    rng = np.random.default_rng(seed)
    Zb = rng.normal(size=(nb, 3)) * [1.5, 1.0, 0.5]
    Za = rng.normal(size=(na, 3)) * [1.2, 1.1, 0.7]
    Zb, Za = ((Z - Z.mean(0)) / np.sqrt(((Z - Z.mean(0)) ** 2).sum(1).mean()) for Z in (Zb, Za))
    return Zb, Za, rng.normal(size=(nb, 4)), rng.normal(size=(na, 4))


def _draw(nb, na, n=60, seed=1):
    rng = np.random.default_rng(seed)
    rows_b, rows_a = rng.choice(nb, n // 2), rng.choice(na, n - n // 2)
    return rows_b, rows_a, rng.permutation(n)


class _Recorder:
    """A stand-in fit: records the clouds it is fitted on and each carry's times; moves by a fixed shift per unit time
    (the earlier cloud towards the later), so forward and back are distinguishable."""

    def __init__(self, shift=(0.0, 0.0, 0.0)):
        self.fits, self.carries, self.shift = [], [], np.asarray(shift, float)

    def __call__(self, Zb, Zt, Fb, Ft, seed):
        self.fits.append((np.array(Zb), np.array(Zt), seed))
        rec = self

        class Fit:
            record = {"data_ratio": 0.3, "p0_norm_mean": 0.1, "lbfgs_n_iter": 7, "device": "cpu"}

            def carry(self, X, a, b):
                rec.carries.append((a, b))
                return np.asarray(X) + (b - a) * rec.shift

            def jacobian_det(self, X, a, b):
                return np.ones(len(X))

        return Fit()


def test_e49s_forms_and_candidates_are_untouched():
    assert si.FORMS == ("stretch_ratio", "stretch_registered", "overlay")
    assert si.ALL_FORMS == (*si.FORMS, "stretch_ratio_anti", "stretch_registered_anti", "stretch_registered_rot",
                            "overlay_rot")
    assert si.GEODESIC_FORMS == ("geodesic_overlay", "geodesic_overlay_rot")
    assert not set(si.GEODESIC_FORMS) & set(si.ALL_FORMS)
    assert len(harness.SHAPE_CANDIDATES) == 7 and "fgw_shape0x_geo" not in harness.SHAPE_CANDIDATES
    assert harness.GEODESIC_CANDIDATES == {"fgw_shape0x_geo": "geodesic_overlay",
                                           "fgw_shape0x_geo_rot": "geodesic_overlay_rot"}
    for name, form in harness.GEODESIC_CANDIDATES.items():
        c = harness.CANDIDATES[name]
        assert c.params["shape_form"] == form and c.params["positions"] == tb.POSITIONS[1]
        assert c.requires is harness.CANDIDATES["fgw_shape0x_ov"].requires


def test_zero_momenta_give_overlay_exactly(monkeypatch):
    Zb, Za, Fb, Fa = _pair()
    rows_b, rows_a, perm = _draw(len(Zb), len(Za))
    R = si.random_rotation(7)
    zero = lddmm.GeodesicFit(q0=Zb[:5].copy(), p0=np.zeros((5, 3)), record={"data_ratio": 1.0})
    monkeypatch.setattr(lddmm, "fit_cached", lambda *a, **k: zero)
    g = si.geodesic_overlay(Zb, Za, rows_b, rows_a, perm, R, 0.4, t=0.5, seed=3, Fb=Fb, Fa=Fa)
    o = si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.4)
    np.testing.assert_allclose(g.Z, o.Z, atol=1e-12)
    np.testing.assert_allclose(g.pool, o.pool, atol=1e-12)
    np.testing.assert_array_equal(g.pool_weights, o.pool_weights)
    assert g.record["jac_min"] == pytest.approx(1.0)
    assert all(not isinstance(v, (list, dict, tuple, np.ndarray)) for v in g.record.values())


@pytest.mark.parametrize("form", si.GEODESIC_FORMS)
def test_earlier_rows_go_forward_later_rows_back_on_one_fit(monkeypatch, form):
    Zb, Za, Fb, Fa = _pair()
    rows_b, rows_a, perm = _draw(len(Zb), len(Za))
    R = si.random_rotation(11)
    rec = _Recorder(shift=(1.0, 0.0, 0.0))
    monkeypatch.setattr(lddmm, "fit_cached", rec)
    out = si.apply_form(form, Z=None, Zb=Zb, Za=Za, t=0.3, R=R, seed=1_000_005, rows_b=rows_b, rows_a=rows_a,
                        perm=perm, t_mix=0.3, features=(Fb, Fa))
    assert len(rec.fits) == 1
    turn = R if form == "geodesic_overlay" else si.random_rotation(1_000_005)
    np.testing.assert_allclose(rec.fits[0][0], Zb)
    np.testing.assert_allclose(rec.fits[0][1], Za @ turn)
    assert rec.fits[0][2] == 1_000_005
    assert (0.0, 0.3) in rec.carries and (1.0, 0.3) in rec.carries
    want = np.vstack([Zb[rows_b] + [0.3, 0, 0], (Za @ turn)[rows_a] - [0.7, 0, 0]])[perm]
    np.testing.assert_allclose(out.Z, want)
    np.testing.assert_allclose(out.pool, np.vstack([Zb + [0.3, 0, 0], Za @ turn - [0.7, 0, 0]]))
    assert out.record["form"] == form and out.record["lddmm_data_ratio"] == 0.3
    assert (out.record.get("rot_stream") == si.ROT_STREAM) == form.endswith("_rot")
    for k in ("jac_min", "thin_sd", "thin_sd_before", "thin_sd_after", "thin_between"):
        assert k in out.record


def test_the_forms_need_features_and_the_old_forms_refuse_them():
    Zb, Za, Fb, Fa = _pair()
    rows_b, rows_a, perm = _draw(len(Zb), len(Za))
    kw = {"Z": Zb[:60], "Zb": Zb, "Za": Za, "t": 0.5, "R": np.eye(3), "seed": 1, "rows_b": rows_b, "rows_a": rows_a,
          "perm": perm, "t_mix": 0.5}
    with pytest.raises(ValueError, match="needs features"):
        si.apply_form("geodesic_overlay", **kw)
    with pytest.raises(ValueError, match="takes no features"):
        si.apply_form("overlay", **kw, features=(Fb, Fa))
    with pytest.raises(ValueError, match="one row per cell"):
        si.geodesic_overlay(Zb, Za, rows_b, rows_a, perm, np.eye(3), 0.5, t=0.5, seed=1, Fb=Fb[:-1], Fa=Fa)
    with pytest.raises(si.ShapeRefused):
        si.geodesic_overlay(Zb, Za, rows_b, rows_a, perm, -np.eye(3), 0.5, t=0.5, seed=1, Fb=Fb, Fa=Fa)


def test_a_failed_fit_is_a_refusal(monkeypatch):
    Zb, Za, Fb, Fa = _pair()
    rows_b, rows_a, perm = _draw(len(Zb), len(Za))

    def bad(*a, **k):
        raise ValueError("the data term at the rigid start is not positive")

    monkeypatch.setattr(lddmm, "fit_cached", bad)
    with pytest.raises(si.ShapeRefused, match="rigid start"):
        si.geodesic_overlay(Zb, Za, rows_b, rows_a, perm, np.eye(3), 0.5, t=0.5, seed=1, Fb=Fb, Fa=Fa)


def _small(monkeypatch):
    monkeypatch.setattr(lddmm, "SUBSAMPLE", 120)
    monkeypatch.setattr(lddmm, "N_CONTROL", 30)
    monkeypatch.setattr(lddmm, "LBFGS_ITERS", 4)
    monkeypatch.setattr(lddmm, "_CACHE", {})
    monkeypatch.setenv(lddmm.DEVICE_ENV, "cpu")


def test_through_interpolate_clouds_only_positions_move(monkeypatch):
    """interpolate_clouds with a geodesic form: fgw_shape0x's rows, order, stage labels and count; the fit sees every
    cell's joint PCs of the coupling's anchors."""
    _small(monkeypatch)
    rng = np.random.default_rng(5)
    Cb = rng.normal(size=(200, 3)) * [3.0, 2.0, 1.0] + 10
    Ca = rng.normal(size=(240, 3)) * [3.5, 2.0, 1.4] - 4
    Eb, Ea = rng.gamma(2.0, 1.0, size=(200, 12)), rng.gamma(2.0, 1.2, size=(240, 12))
    seen = []
    real = lddmm.fit_cached

    def spy(Zb, Zt, Fb, Ft, seed):
        seen.append((Fb.shape, Ft.shape))
        return real(Zb, Zt, Fb, Ft, seed)

    monkeypatch.setattr(lddmm, "fit_cached", spy)
    kw = {"target_rms": 2.0, "t_pos": 0.0, "positions": tb.POSITIONS[1], "n_anchor": 80, "n_pc": 5}
    base = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 90, 1_000_002, **kw)
    ov = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 90, 1_000_002, shape_form="overlay", **kw)
    geo = tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 90, 1_000_002, shape_form="geodesic_overlay", **kw)
    assert seen == [((200, 5), (240, 5))]
    for r in (ov, geo):
        np.testing.assert_array_equal(r.rows, base.rows)
        np.testing.assert_array_equal(r.stage, base.stage)
    assert geo.coords.shape == base.coords.shape and not np.allclose(geo.coords, ov.coords)
    assert geo.record["shape_interp"] == "geodesic_overlay"
    assert 0 < geo.record["shape_interp_lddmm_data_ratio"] and geo.record["shape_interp_jac_min"] > 0
    assert tb.rms_radius(geo.coords) == pytest.approx(2.0, rel=1e-6)  # the size rule runs after the form
    with pytest.raises(ValueError, match="shape_form must be one of"):
        tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, 90, 1, shape_form="geodesic", **kw)


def test_the_harness_candidates_build_end_to_end_on_a_synthetic_pair(monkeypatch):
    """Both candidates on test_e47_floor_exact's synthetic embryo pair (a small fit): fgw_shape0x's expression rows
    and count, other positions, and a cloud guard that either passes or refuses with its own message."""
    from test_e47_floor_exact import _ctx as e47_ctx

    _small(monkeypatch)
    ctx = e47_ctx()
    c0 = harness.CANDIDATES["fgw_shape0x"].build(ctx)
    for name in harness.GEODESIC_CANDIDATES:
        try:
            out = harness.CANDIDATES[name].build(ctx)
        except ValueError as e:
            assert "cloud guard flags" in str(e)
            continue
        np.testing.assert_array_equal(out.X, c0.X)
        assert out.n == c0.n and not np.allclose(out.coords, c0.coords)
