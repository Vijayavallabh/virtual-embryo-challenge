"""Synthetic checks of size_procedure.py (no competition data). Run: uv run pytest experiments/coords -q"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import size_procedure as sp  # noqa: E402


def cloud(sd, n=4000, seed=0, core_frac=0.5, step=15.0):
    """Gaussian cloud with per-axis SDs `sd`, z snapped to sections of the given step; first half is core."""
    rng = np.random.default_rng(seed)
    C = rng.normal(size=(n, 3)) * np.asarray(sd, float)
    C[:, 2] = np.round(C[:, 2] / step) * step
    core = np.zeros(n, bool)
    core[: int(n * core_frac)] = True
    return C, core


def geom(name, t, scale, sd=(100.0, 60.0, 30.0), seed=0):
    C, core = cloud(np.asarray(sd) * scale, seed=seed)
    return sp.geometry_from_arrays(name, t, C, core)


def test_rules_on_exact_log_linear_growth():
    ts = np.array([9.0, 9.5, 10.5])
    ys = 5.0 + 0.4 * (ts - 9.0)
    assert sp.rule_latest_rate(ts, ys, 12.5) == pytest.approx(5.0 + 0.4 * 3.5)
    assert sp.rule_half_rate(ts, ys, 12.5) == pytest.approx(ys[-1] + 0.2 * 2.0)
    # constant rates: no deceleration evidence, decaying_rate falls back to the latest rate
    assert sp.rule_decaying_rate(ts, ys, 12.5) == pytest.approx(sp.rule_latest_rate(ts, ys, 12.5))
    assert sp.rule_copy_last(ts, ys, 12.5) == ys[-1]


def test_decaying_rate_decelerates_and_power_law_is_between():
    ts = np.array([9.0, 9.5, 10.5])
    ys = np.array([5.0, 5.4, 5.8])  # rates 0.8/day then 0.4/day
    d = sp.rule_decaying_rate(ts, ys, 12.5)
    assert ys[-1] < d < sp.rule_latest_rate(ts, ys, 12.5)
    p = sp.rule_power_t(ts, ys, 12.5)
    assert sp.rule_half_rate(ts, ys, 12.5) < p < sp.rule_latest_rate(ts, ys, 12.5)


def test_applicability_and_validity():
    early = [geom("a", 8.25, 1.0), geom("b", 8.75, 1.0, seed=1), geom("c", 9.5, 1.8, seed=2)]
    res = sp.heart_size_from_geometry(early, 12.5, z_mode="raw")
    assert res["valid"] is False  # the only rate interval starts at 8.75 < 9.0
    late = [geom("b", 8.75, 1.0, seed=1), geom("c", 9.5, 1.8, seed=2), geom("d", 10.5, 2.4, seed=3)]
    res = sp.heart_size_from_geometry(late, 12.5, z_mode="raw")
    assert res["valid"] is True
    assert set(res["applicable_core"]) == {"latest_rate", "half_rate", "power_t"}  # decaying needs 8.75 >= 9.0
    assert set(res["chosen_core"]) == set(res["applicable_core"])  # no eligible pseudo-target: average all
    expected = np.mean([res["all_rule_predictions_core"][r] for r in res["chosen_core"]])
    assert np.log(res["pred_rms_core"]) == pytest.approx(expected)


def test_scope_hypotheses():
    s, ts = [0.8, 0.4], [9.5, 10.5]
    assert sp.scope_log_ratio(s, ts, 12.5, "H0") == pytest.approx(0.4)
    assert sp.scope_log_ratio(s, ts, 12.5, "H2") == pytest.approx(0.2)
    assert sp.scope_log_ratio(s, ts, 12.5, "H1") == pytest.approx(0.4 * 0.5 ** 2)
    assert sp.scope_log_ratio([0.3, 0.4], ts, 12.5, "H1") == pytest.approx(0.4)  # increasing: held


def test_auto_z_mode_uses_in_plane_growth_when_the_thin_axis_is_z():
    flat = [geom("c", 9.5, 1.0, sd=(100, 80, 30), seed=2), geom("d", 10.5, 1.3, sd=(100, 80, 30), seed=3)]
    assert sp.resolve_heart_z_mode(flat, "auto") == "xy_grow_z_hold"
    tall = [geom("c", 9.5, 1.0, sd=(30, 80, 100), seed=2), geom("d", 10.5, 1.3, sd=(30, 80, 100), seed=3)]
    assert sp.resolve_heart_z_mode(tall, "auto") == "raw"
    res = sp.heart_size_from_geometry([geom("b", 9.0, 0.8, sd=(100, 80, 30), seed=1)] + flat, 12.5)
    assert res["z_mode"] == "xy_grow_z_hold"
    # z is held at the reference: predicted core RMS^2 = predicted xy^2 + reference z SD^2
    ref = flat[-1]
    assert res["pred_rms_core"] == pytest.approx(np.hypot(res["pred_rms_xy_core"], ref.sd_z_core))


def test_embryo_pchip_and_loo():
    gs = [geom("e1", 6.75, 1.0), geom("e2", 7.25, 1.2, seed=1), geom("e3", 7.5, 1.5, seed=2),
          geom("e4", 8.0, 2.4, seed=3)]
    res = sp.embryo_size_from_geometry(gs, 7.75)
    assert gs[2].rms_full < res["pred_rms_full"] < gs[3].rms_full  # monotone data: stays inside the bracket
    assert res["bracket"] == ["e3", "e4"]
    assert [r["left_out"] for r in res["loo_interior"]] == ["e2", "e3"]
    with pytest.raises(ValueError):
        sp.embryo_size_from_geometry(gs, 8.5)  # extrapolation is refused


def test_rescale_to_rms_is_a_proper_scaling():
    C, _ = cloud((50, 40, 30), seed=5)
    D = sp.rescale_to_rms(C, 123.0)
    assert sp.rms_radius(D) == pytest.approx(123.0)
    np.testing.assert_allclose(D.mean(0), C.mean(0), atol=1e-9)
    f = sp.rms_radius(D) / sp.rms_radius(C)
    np.testing.assert_allclose(D - D.mean(0), f * (C - C.mean(0)), atol=1e-8)


def test_no_size_numbers_in_preregistered_settings():
    num = [v for v in sp.PREREGISTERED.values() if isinstance(v, (int, float))]
    # only time cut-offs (E-days), an angle and a log tolerance; nothing on the scale of an RMS radius
    assert all(v < 100 for v in num)
