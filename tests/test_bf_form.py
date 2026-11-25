"""Row R71's form BF (vec.shape_interp.overlay_bayes) and its gate's controls, on synthetic clouds only: the canonical
frame is the scorer's (veckit main), the density is flip-symmetric with mass 1, no step branches on the cell count,
rows, order, count and expression are the base's, the RMS radius lands on the size target, two builds are equal, the
pool moves with the field, and the truth never enters a practice task's prior. Nothing is read from data/ and nothing
is scored."""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import shape_interp as si
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def fewer_anchors(monkeypatch):
    """The exact assignment on 600 draws in place of 4,000 (the same code path; the constant is pinned below by
    the module's text), so the tests run in seconds."""
    monkeypatch.setattr(si, "BAYES_ANCHORS", 600)


def _cloud(n, scale, seed, shift=(0.0, 0.0, 0.0)):
    return np.random.default_rng(seed).normal(size=(n, 3)) * np.asarray(scale) + np.asarray(shift)


def _prior(names=("E6.75", "E8.0", "E9.5"), n=900):
    return si.bayes_prior({s: _cloud(n, (1.6 - 0.1 * k, 0.9, 0.35 + 0.05 * k), 40 + k) for k, s in enumerate(names)})


def _overlay_args(nb=700, na=800, n=650, seed=3):
    rng = np.random.default_rng(seed)
    Zb, Za = _cloud(nb, (1.5, 1.0, 0.4), 1), _cloud(na, (1.3, 1.1, 0.5), 2, (0.2, 0.0, 0.0))
    rows_b, rows_a = rng.integers(0, nb, n // 2), rng.integers(0, na, n - n // 2)
    return Zb, Za, rows_b, rows_a, rng.permutation(n)


def test_the_canonical_frame_is_the_scorers_and_inverts():
    import common.shape_metrics as sm

    Q = np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))[0]
    C = _cloud(500, (3.0, 2.0, 0.5), 7, (10.0, -4.0, 2.0)) @ Q
    Z, mu, Vt, rms = si.canonical_frame(C)
    want, want_rms = sm._canonicalise(C)
    np.testing.assert_allclose(Z, want, atol=1e-12)
    assert rms == pytest.approx(want_rms) and np.linalg.det(Vt) > 0
    np.testing.assert_allclose(si.from_canonical(Z, mu, Vt, rms), C, atol=1e-9)
    with pytest.raises(si.ShapeRefused):
        si.canonical_frame(np.ones((10, 3)))


def test_occupancy_and_dice_are_the_scorers():
    import common.shape_metrics as sm

    A, B = _cloud(400, (1.4, 1.0, 0.3), 5), _cloud(400, (1.2, 0.8, 0.5), 6)
    Za, Zb = sm._canonicalise(A)[0], sm._canonicalise(B)[0]
    assert si.dice(si.occupancy(Zb), si.occupancy(Zb)) == 1.0
    np.testing.assert_array_equal(si.occupancy(Za), tb._voxel_ids(Za))
    best = max(si.dice(si.occupancy(Za * np.asarray(f, float)), si.occupancy(Zb)) for f in si.PROPER_FLIPS)
    assert best == pytest.approx(sm.occupancy_dice(A, B)[0])  # equal counts: _match_n keeps every row
    assert si.PROPER_FLIPS == tuple(tuple(f) for f in sm._PROPER_FLIPS) and len(si.AXIS_SIGNS) == 8


def test_the_density_is_flip_symmetric_with_mass_one_and_the_no_flip_one_is_not():
    S = si.canonical_frame(_cloud(3000, (1.5, 0.8, 0.3), 9, (0.4, 0.2, 0.0)))[0] + [0.3, 0.1, 0.05]
    q = si.specimen_density(S)
    assert q.shape == (48, 48, 48) and q.sum() == pytest.approx(1.0) and (q >= 0).all()
    for ax in range(3):
        np.testing.assert_allclose(q, np.flip(q, axis=ax), atol=1e-15)
    q1 = si.specimen_density(S, signs=((1, 1, 1),))
    assert q1.sum() == pytest.approx(1.0) and not np.allclose(q1, np.flip(q1, axis=0))
    prior = _prior()
    qc, flips = si.consensus_density(prior)
    assert flips == () and qc.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(qc, np.mean([si.specimen_density(c) for c in prior.clouds], axis=0))


def test_the_best_proper_flip_recovers_a_turned_specimen():
    X = si.canonical_frame(_cloud(2000, (1.5, 1.0, 0.4), 11))[0]
    X = X + np.where(X[:, :1] > 0, [0.0, 0.6, 0.2], [0.0, -0.2, 0.0])  # an asymmetric cloud
    for f in si.PROPER_FLIPS:
        got, d = si.best_proper_flip(X * np.asarray(f, float), X)
        assert tuple(np.asarray(got) * np.asarray(f)) == (1, 1, 1) and d == 1.0


def test_no_step_branches_on_the_cell_count():
    """The challenge's fix 3: 4,000 rows drawn with replacement and 4,000 draws from q at every count, one code path;
    the design's 12,000-cell switch is gone."""
    src = inspect.getsource(si.transport_displacements) + inspect.getsource(si.overlay_bayes)
    assert "12000" not in src and "12_000" not in src and "min(n" not in src and "len(X) >" not in src
    assert "entropic" not in src.lower() and "sinkhorn" not in src.lower()
    q = si.consensus_density(_prior())[0]
    for n in (40, 300, 5000):  # below and above the anchor count: the same path
        X = si.canonical_frame(_cloud(n, (1.4, 1.0, 0.4), n))[0]
        d = si.transport_displacements(X, q, 7)
        assert d.shape == (n, 3) and np.isfinite(d).all()
        np.testing.assert_array_equal(d, si.transport_displacements(X, q, 7))
    text = (REPO / "src/vec/shape_interp.py").read_text()
    assert "\nBAYES_ANCHORS = 4000  #" in text and (si.BAYES_Y_STREAM, si.BAYES_X_STREAM) == (61, 62)


def test_the_form_keeps_rows_order_count_and_pool_weights_and_is_deterministic():
    Zb, Za, rows_b, rows_a, perm = _overlay_args()
    R = np.eye(3)
    ov = si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.5)
    prior = _prior()
    out = si.overlay_bayes(Zb, Za, rows_b, rows_a, perm, R, 0.5, seed=1000001, prior=prior)
    again = si.overlay_bayes(Zb, Za, rows_b, rows_a, perm, R, 0.5, seed=1000001, prior=prior)
    assert out.Z.shape == ov.Z.shape and out.pool.shape == ov.pool.shape
    np.testing.assert_array_equal(out.Z, again.Z)
    np.testing.assert_array_equal(out.pool, again.pool)
    np.testing.assert_array_equal(out.pool_weights, ov.pool_weights)
    assert not np.allclose(out.Z, ov.Z)
    other = si.overlay_bayes(Zb, Za, rows_b, rows_a, perm, R, 0.5, seed=1000002, prior=prior)
    assert not np.array_equal(other.Z, out.Z)  # the build seed sets streams 61 and 62
    rec = out.record
    assert rec["variant"] == "arm" and rec["prior"] == "E6.75|E8.0|E9.5" and rec["flip_copies"] == 8
    assert not any(v == len(ov.Z) for v in rec.values() if isinstance(v, int) and not isinstance(v, bool))


def test_the_pool_moves_with_the_field_and_a_union_cell_moves_like_its_pool_twin():
    Zb, Za, rows_b, rows_a, perm = _overlay_args(nb=300, na=300, n=300)
    rows_b, rows_a = np.arange(150), np.arange(150)  # every drawn row is a distinct pool cell
    ov = si.overlay(Zb, Za, rows_b, rows_a, perm, np.eye(3), 0.5)
    out = si.overlay_bayes(Zb, Za, rows_b, rows_a, perm, np.eye(3), 0.5, seed=5, prior=_prior())
    union_rows = np.r_[rows_b, len(Zb) + rows_a][perm]  # pool row of each union row
    np.testing.assert_allclose(out.pool[union_rows], out.Z, atol=1e-9)
    assert not np.allclose(out.pool, ov.pool)


def test_the_controls_take_their_own_priors():
    Zb, Za, rows_b, rows_a, perm = _overlay_args()
    kw = dict(Z=None, Zb=Zb, Za=Za, t=0.5, R=np.eye(3), seed=3, rows_b=rows_b, rows_a=rows_a, perm=perm, t_mix=0.5)
    b2 = si.apply_form("overlay_bayes_brackets", **kw)
    assert b2.record["prior"] == "before|after" and b2.record["form"] == "overlay_bayes_brackets"
    with pytest.raises(ValueError, match="takes no prior"):
        si.apply_form("overlay_bayes_brackets", **kw, prior=_prior())
    ns = si.apply_form("overlay_bayes_noflip", **kw, prior=_prior())
    assert ns.record["flip_copies"] == 1 and len(ns.record["flips"].split("|")) == 3
    for form in ("overlay_bayes", "overlay_bayes_noflip"):
        with pytest.raises(ValueError, match="needs a BayesPrior"):
            si.apply_form(form, **kw)
    for form in ("overlay", "overlay_stretch", "stretch_ratio"):
        with pytest.raises(ValueError, match="takes no prior"):
            si.apply_form(form, **kw, prior=_prior())
    with pytest.raises(ValueError, match="takes no features"):
        si.apply_form("overlay_bayes", **kw, prior=_prior(), features=(np.zeros((1, 1)), np.zeros((1, 1))))
    assert set(si.BAYES_CONTROLS).isdisjoint(si.ALL_FORMS) and set(si.BAYES_FORMS).isdisjoint(si.ALL_FORMS)


def test_the_prior_takes_released_stages_once_each_in_day_order():
    p = si.bayes_prior({"E9.5": _cloud(50, (1, 1, 1), 1), "E6.75": _cloud(50, (1, 1, 1), 2)})
    assert p.names == ("E6.75", "E9.5")
    for bad in ({}, {"E7.75": _cloud(50, (1, 1, 1), 1)}, {"E8.5": _cloud(50, (1, 1, 1), 1)},
                {"E10.5": _cloud(50, (1, 1, 1), 1)}):
        with pytest.raises(ValueError):
            si.bayes_prior(bad)
    assert si.BAYES_PRIOR_STAGES == ("E6.75", "E7.25", "E7.5", "E8.0", "E8.25_late", "E8.75", "E9.5")


def test_the_truth_never_enters_a_practice_tasks_prior():
    sys.path.insert(0, str(REPO / "experiments" / "coords"))
    import bf_read

    assert bf_read.prior_stages(None) == si.BAYES_PRIOR_STAGES
    for truth, stage in (("T2/E8.0.h5ad", "E8.0"), ("T2/E7.5.h5ad", "E7.5"), ("T2/E7.25.h5ad", "E7.25"),
                         ("T2/E8.25_late.h5ad", "E8.25_late")):
        got = bf_read.prior_stages(truth)
        assert stage not in got and len(got) == 6
    assert bf_read.prior_stages("T2/E8.5.h5ad") == si.BAYES_PRIOR_STAGES  # E8.5 is never in the prior


def _interp(shape_form, prior=None, n=500, seed=1000001):
    rng = np.random.default_rng(0)
    Cb, Ca = _cloud(700, (6.0, 4.0, 1.5), 1, (50, 20, 5)), _cloud(800, (5.0, 4.5, 2.0), 2, (-10, 3, 1))
    Eb, Ea = rng.poisson(2.0, (700, 6)).astype(np.float32), rng.poisson(2.5, (800, 6)).astype(np.float32)
    kw = {} if prior is None else {"shape_prior": prior}
    return tb.interpolate_clouds(Cb, Ca, Eb, Ea, 0.5, n, seed, target_rms=7.0, t_pos=0.0, positions="floor_exact",
                                 shape_form=shape_form, n_anchor=150, **kw)


def test_through_interpolate_clouds_only_positions_change_and_the_size_rule_holds():
    o = _interp("overlay")
    bf = _interp("overlay_bayes", _prior())
    again = _interp("overlay_bayes", _prior())
    np.testing.assert_array_equal(bf.rows, o.rows)
    np.testing.assert_array_equal(bf.stage, o.stage)
    np.testing.assert_array_equal(bf.coords, again.coords)
    assert not np.allclose(bf.coords, o.coords)
    if bf.record["n_jittered"] == 0:
        assert tb.rms_radius(bf.coords) == pytest.approx(7.0, rel=1e-6)
    assert bf.record["target_rms"] == o.record["target_rms"] == 7.0
    assert bf.record["shape_interp"] == "overlay_bayes" and bf.record["shape_interp_variant"] == "arm"
    b2 = _interp("overlay_bayes_brackets")
    assert b2.record["shape_interp_prior"] == "before|after"


def test_shape_prior_goes_with_the_bayes_forms_only():
    with pytest.raises(ValueError, match="shape_prior goes with"):
        _interp("overlay", _prior())
    with pytest.raises(ValueError, match="shape_prior goes with"):
        _interp(None, _prior())
    with pytest.raises(ValueError, match="takes no prior"):
        _interp("overlay_bayes_brackets", _prior())


def test_the_neighbour_share_of_g0():
    P = _cloud(400, (1.0, 1.0, 1.0), 3)
    assert si.neighbour_share(P, P) == 1.0
    assert si.neighbour_share(P, P * 2.0 + 1.0) == 1.0  # a similarity keeps every neighbourhood
    shuffled = P[np.random.default_rng(1).permutation(len(P))]
    assert si.neighbour_share(P, shuffled) < 0.2
    sets = si.knn_sets(P, 15)
    assert sets.shape == (400, 15) and not (sets == np.arange(400)[:, None]).any()
