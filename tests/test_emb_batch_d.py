"""The embryo batch D (campaign branch campaign/emb-d): coupling on tissue profiles (vec.t2_build.tissue_profiles,
interpolate_clouds(..., features="tissue"), the arms fgw_tissue and fgw_tissue_series) and the series render with an
exact per-cell total (series_rendered_expression(..., kernel="exact_total"), the arm fgw_aligned_series_xt), and the
runner experiments/harness/n43_embryo_d_arms.py. Synthetic clouds and stages only; no data file, committed result or
scorer is read.

Tolerances stated here: the expected counts against a Monte Carlo mean, 5 standard errors of that mean per entry (and
1e-12 absolute where the draw cannot vary); the mean rendered library of the exact-total render against the
independent render's, 4 standard errors of the independent mean plus 1 count (the two roundings, panel and outside)."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from vec import core, harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
BASE = "217520d"  # campaign/emb-mirror, the branch's base
ARMS = ("fgw_tissue", "fgw_tissue_series", "fgw_aligned_series_xt")
K = harness.NBHD_K


def _helpers(name: str):
    spec = importlib.util.spec_from_file_location(f"{name}_helpers_d", REPO / "tests" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def polar():
    """tests/test_emb_polar.py's synthetic pair, base loader and comparison (no data file)."""
    return _helpers("test_emb_polar")


@pytest.fixture(scope="module")
def emb_expr():
    """The synthetic embryo stages and contexts of tests/test_emb_expr.py (no data file)."""
    return _helpers("test_emb_expr")


def _hand_profile(C: np.ndarray, E: np.ndarray, i: int, k: int) -> np.ndarray:
    """The cell itself and its k - 1 nearest other cells by a full sort of the distances (no tree)."""
    d = np.linalg.norm(C - C[i], axis=1)
    d[i] = -1.0
    near = np.argsort(d, kind="stable")[:k]
    assert near[0] == i
    return E[near].mean(0)


# ================================================================================ coupling on tissue profiles
# ------------------------------------------------------------------------------------------ T1. unset is the entry
@pytest.mark.parametrize("move_coupling", tb.MOVE_COUPLINGS)
def test_features_unset_is_the_entry_bit_for_bit(tmp_path, polar, move_coupling):
    A, B, E, EB = polar._pair()
    base = polar._base_module(tmp_path, BASE)
    old = base.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True, **polar.KW)
    new = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True, **polar.KW)
    own = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True,
                                features="own", **polar.KW)
    polar._same(old, new)
    polar._same(new, own)
    assert np.array_equal(old.partners, new.partners) and np.array_equal(new.partners, own.partners)
    assert "features" not in new.record
    for kw in ({"partners": "mirror_fixed"}, {"partners": "mirror_gated"}, {"moves": "polar"}):  # the base's arms
        polar._same(base.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, **kw, **polar.KW),
                    tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, **kw, **polar.KW))
    tissue = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=move_coupling, with_partners=True,
                                   features="tissue", tissue_k=K, **polar.KW)
    assert not np.array_equal(tissue.coords, new.coords)  # the idea is not idle here


def test_features_and_tissue_k_are_checked(polar):
    A, B, E, EB = polar._pair()
    with pytest.raises(ValueError, match="features"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, features="niche", **polar.KW)
    with pytest.raises(ValueError, match="tissue_k"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, features="tissue", **polar.KW)
    with pytest.raises(ValueError, match="tissue_k"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, tissue_k=K, **polar.KW)


# ------------------------------------------------------------------------------------------ T2. the profile by hand
def test_a_cells_profile_is_the_mean_over_itself_and_its_14_nearest_cells():
    rng = np.random.default_rng(7)
    C = rng.normal(size=(40, 3)) * [3.0, 2.0, 1.0]
    E = rng.gamma(2.0, 1.0, size=(40, 6))
    rows = np.array([0, 5, 17, 39, 5])
    got = tb.tissue_profiles(C, E, rows, K)
    assert K == 15 and got.shape == (len(rows), 6)
    for j, i in enumerate(rows):
        np.testing.assert_allclose(got[j], _hand_profile(C, E, i, K), rtol=1e-12, atol=0.0)
    # coincident positions: the cell itself is always one of the 15
    D = C.copy()
    D[:20] = D[0]
    got = tb.tissue_profiles(D, E, np.arange(20), K)
    for i in range(20):
        dist = np.linalg.norm(D - D[i], axis=1)
        others = np.flatnonzero(dist == 0.0)
        assert len(others) == 20  # 20 cells tie at distance 0; the profile holds the cell and 14 of the others
        near = np.argsort(np.where(np.arange(40) == i, -1.0, dist), kind="stable")[:K]
        assert i in near
    tree_rows = np.asarray(tb.cKDTree(D).query(D[:20], k=K)[1])
    lost = ~(tree_rows == np.arange(20)[:, None]).any(1)
    expect = tree_rows.copy()
    expect[lost, -1] = np.arange(20)[lost]
    np.testing.assert_array_equal(got, E[expect].mean(1))
    assert all(i in expect[i] for i in range(20))


def test_a_bracket_smaller_than_15_cells_gives_the_mean_over_all_of_them(polar):
    rng = np.random.default_rng(8)
    C, E = rng.normal(size=(9, 3)), rng.gamma(2.0, 1.0, size=(9, 5))
    got = tb.tissue_profiles(C, E, np.arange(9), K)
    np.testing.assert_allclose(got, np.broadcast_to(E.mean(0), got.shape), rtol=1e-12)
    one = tb.tissue_profiles(C[:1], E[:1], np.array([0]), K)
    assert np.array_equal(one, E[:1])
    # through interpolate_clouds: an after bracket of 12 cells
    A, B, EA, EB = polar._pair()
    res = tb.interpolate_clouds(A, B[:12], EA, EB[:12], 0.5, 200, features="tissue", tissue_k=K, **polar.KW)
    assert np.isfinite(res.coords).all() and res.record["features"] == "tissue"


# ------------------------------------------------------------------------------------------ T3. the spy
@pytest.mark.parametrize("move_coupling", tb.MOVE_COUPLINGS)
def test_tissue_features_go_into_the_joint_pcs_and_nothing_else_changes(polar, monkeypatch, move_coupling):
    """With "tissue", `_joint_pcs` receives the anchors' tissue profiles, the FGW coupling (and the aligned coupling)
    receives those joint PCs with the entry's structure, alpha and iteration arguments; with the couplings held at the
    entry's, every later step gives the entry's output bit for bit."""
    A, B, E, EB = polar._pair()
    kw = dict(move_coupling=move_coupling, with_partners=True, with_coupling=True, **polar.KW)
    real_pcs, real_fgw, real_al = tb._joint_pcs, tb.fgw_coupling, tb.aligned_coupling
    calls: dict[str, list] = {"pcs": [], "fgw": [], "al": []}

    def pcs(Eb, Ea, n_pc):
        out = real_pcs(Eb, Ea, n_pc)
        calls["pcs"].append((np.array(Eb), np.array(Ea), n_pc, out))
        return out

    def fgw(*a):
        out = real_fgw(*a)
        calls["fgw"].append(([np.array(x) if isinstance(x, np.ndarray) else x for x in a], out))
        return out

    def al(*a):
        out = real_al(*a)
        calls["al"].append(([np.array(x) if isinstance(x, np.ndarray) else x for x in a], out))
        return out

    monkeypatch.setattr(tb, "_joint_pcs", pcs)
    monkeypatch.setattr(tb, "fgw_coupling", fgw)
    monkeypatch.setattr(tb, "aligned_coupling", al)
    entry = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, **kw)
    tissue = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, features="tissue", tissue_k=K, **kw)
    (_, _, n0, _), (Pb, Pa, n1, (Fb, Fa)) = calls["pcs"]
    ib, ia = tissue.coupling["anchors_before"], tissue.coupling["anchors_after"]
    assert np.array_equal(ib, entry.coupling["anchors_before"]) and np.array_equal(ia, entry.coupling["anchors_after"])
    assert n0 == n1 == tb.INTERP_N_PC
    np.testing.assert_array_equal(Pb, tb.tissue_profiles(A, E, ib, K))
    np.testing.assert_array_equal(Pa, tb.tissue_profiles(B, EB, ia, K))
    (a0, T0), (a1, T1) = calls["fgw"]
    assert np.array_equal(a1[0], Fb) and np.array_equal(a1[1], Fa)  # the profiles' joint PCs
    for x, y in zip(a0[2:], a1[2:], strict=True):  # structure, alpha, max_iter, emd_max_iter: the entry's
        assert np.array_equal(x, y)
    assert not np.array_equal(T0, T1)
    if move_coupling == "aligned":
        (b0, M0), (b1, M1) = calls["al"]
        assert np.array_equal(b1[0], Fb) and np.array_equal(b1[1], Fa)
        assert np.array_equal(b0[2], b1[2]) and all(np.array_equal(x, y) for x, y in zip(b0[4:], b1[4:], strict=True))
    else:
        assert calls["al"] == []
    # the couplings held at the entry's: the tissue run is the entry, bit for bit, the record gains "features" only
    monkeypatch.setattr(tb, "fgw_coupling", lambda *a: T0)
    if move_coupling == "aligned":
        monkeypatch.setattr(tb, "aligned_coupling", lambda *a: calls["al"][0][1])
    held = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, features="tissue", tissue_k=K, **kw)
    polar._same(dataclasses.replace(held, record={k: v for k, v in held.record.items() if k != "features"}), entry)
    assert np.array_equal(held.partners, entry.partners)
    assert set(held.record) - set(entry.record) == {"features"} and held.record["features"] == "tissue"


# ------------------------------------------------------------------------------------------ T4. rows and expression
def test_the_drawn_rows_are_the_entrys_and_only_positions_differ(polar):
    A, B, E, EB = polar._pair()
    for mc in tb.MOVE_COUPLINGS:
        entry = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=mc, **polar.KW)
        tissue = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, move_coupling=mc, features="tissue", tissue_k=K,
                                       **polar.KW)
        assert np.array_equal(entry.rows, tissue.rows) and np.array_equal(entry.stage, tissue.stage)
        assert np.array_equal(entry.pool_weights, tissue.pool_weights)
        assert not np.array_equal(entry.coords, tissue.coords)


@pytest.mark.parametrize("arm, base", [("fgw_tissue", "fgw"), ("fgw_tissue_series", "fgw_series")])
def test_the_tissue_arms_build_on_a_synthetic_embryo_and_move_only_the_positions(emb_expr, arm, base):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    a, again, b = harness.CANDIDATES[arm].build(c), harness.CANDIDATES[arm].build(c), harness.CANDIDATES[base].build(c)
    assert emb_expr.digest(a) == emb_expr.digest(again)  # same seed, same output
    assert np.array_equal(a.X, b.X)  # the same draws and each cell's expression
    assert not np.array_equal(a.coords, b.coords) and np.isfinite(a.coords).all()


# ------------------------------------------------------------------------------------------ T5. reproducible
def test_tissue_is_reproducible(polar):
    A, B, E, EB = polar._pair()
    for mc in tb.MOVE_COUPLINGS:
        a = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, move_coupling=mc, features="tissue", tissue_k=K, **polar.KW)
        b = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, move_coupling=mc, features="tissue", tissue_k=K, **polar.KW)
        polar._same(a, b)
    other = tb.interpolate_clouds(A, B, E, EB, 0.5, 300, features="tissue", tissue_k=K, **{**polar.KW, "seed": 2})
    assert not np.array_equal(a.coords, other.coords)


# ================================================================================ the exact-total series render
def _render(emb_expr, monkeypatch, kernel=None, seed=(9, 7), repair=(9, harness.SERIES_XT_STREAM)):
    """The series render on the synthetic before/after pair at the aligned positions, with the rendered cells' counts
    and library sizes captured as `_log_normalise` receives them and the expected counts as computed."""
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    res, _, _, frac = emb_expr._interp(c, "aligned")
    seen: dict = {}
    real_norm, real_exp, real_fix = tb._log_normalise, tb.series_expected_counts, tb._exact_series_totals

    def norm(C, lib):
        seen["counts"], seen["lib"] = sparse.csr_matrix(C).toarray(), np.asarray(lib, dtype=np.float64)
        return real_norm(C, lib)

    def exp(*a):
        out = real_exp(*a)
        seen["expected"] = out
        return out

    def fix(C1, oth, *a):
        seen["drawn"] = (np.array(C1), np.array(oth))
        return real_fix(C1, oth, *a)

    monkeypatch.setattr(tb, "_log_normalise", norm)
    monkeypatch.setattr(tb, "series_expected_counts", exp)
    monkeypatch.setattr(tb, "_exact_series_totals", fix)
    kw = {} if kernel is None else {"kernel": kernel, "repair_seed": list(repair)}
    X, rec = tb.series_rendered_expression(s["before"], s["after"], res, list(emb_expr.PANEL), frac, list(seed), **kw)
    monkeypatch.undo()
    return X, rec, seen, res, s


# ------------------------------------------------------------------------------------------ X1. unset
def test_kernel_unset_is_fgw_aligned_series_bit_for_bit(tmp_path, polar, emb_expr, monkeypatch):
    base = polar._base_module(tmp_path, BASE)
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    res, _, _, frac = emb_expr._interp(c, "aligned")
    g = list(emb_expr.PANEL)
    X0, r0 = base.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7])
    X1, r1 = tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7])
    X2, r2 = tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], kernel="independent")
    assert np.array_equal(X0, X1) and np.array_equal(X1, X2)
    assert json.dumps(r0, sort_keys=True) == json.dumps(r1, sort_keys=True) == json.dumps(r2, sort_keys=True)
    # the built candidate itself is pinned by sha256 in tests/test_emb_expr.py (PINS_AT_E709655)
    built, xt = (harness.CANDIDATES[n].build(c) for n in ("fgw_aligned_series", "fgw_aligned_series_xt"))
    assert np.array_equal(built.coords, xt.coords)  # the coupling, moves, draws and cells of the second upload
    assert not np.array_equal(built.X, xt.X)
    with pytest.raises(ValueError, match="kernel"):
        tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], kernel="multinomial")
    with pytest.raises(ValueError, match="repair_seed"):
        tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], kernel="exact_total")
    with pytest.raises(ValueError, match="repair_seed"):
        tb.series_rendered_expression(s["before"], s["after"], res, g, frac, [9, 7], repair_seed=[9, 21])


# ------------------------------------------------------------------------------------------ X2. exact totals
def test_every_rendered_cell_is_at_its_rounded_expected_totals(emb_expr, monkeypatch):
    X, rec, seen, res, _ = _render(emb_expr, monkeypatch, kernel="exact_total")
    expected, expected_other = seen["expected"]
    C, lib = seen["counts"], seen["lib"]
    assert len(C) == int((res.stage == 0).sum()) > 0  # the shallower bracket (before) is rendered
    np.testing.assert_array_equal(C.sum(1), np.rint(expected.sum(1)))
    np.testing.assert_array_equal(lib - C.sum(1), np.rint(expected_other))
    assert ((C == 0) | (expected > 0)).all()
    drawn, oth = seen["drawn"]
    assert 0.0 < rec["xt_share_repaired"] <= 1.0
    size = np.abs(drawn.sum(1) - np.rint(expected.sum(1))) + np.abs(oth - np.rint(expected_other))
    assert rec["xt_share_repaired"] == pytest.approx(float((size > 0).mean()))
    assert rec["xt_mean_abs_repair"] == pytest.approx(float(size.mean()))
    # the render's own draws are untouched: what the repair receives is the independent render's output
    _, rec0, seen0, _, _ = _render(emb_expr, monkeypatch)
    np.testing.assert_array_equal(drawn, seen0["counts"])
    np.testing.assert_array_equal(oth, seen0["lib"] - seen0["counts"].sum(1))
    assert set(rec) - set(rec0) == {"xt_share_repaired", "xt_mean_abs_repair"}
    same = set(rec0) - {"rendered_median"}  # the rendered cells' median library follows the repair
    assert {k: rec[k] for k in same} == {k: rec0[k] for k in same}
    assert rec["rendered_median"] == float(np.median(lib))


# ------------------------------------------------------------------------------------------ X3. Monte Carlo
def _mc_render(C0, phi, P, top, other, ratio, reps, seed):
    rng = np.random.default_rng(seed)
    panel = np.zeros(C0.shape)
    sq = np.zeros(C0.shape)
    out_o = np.zeros(len(C0))
    sq_o = np.zeros(len(C0))
    for _ in range(reps):
        C1 = core.apply_gene_factor(C0.astype(np.int32), phi, P, rng).astype(np.int64)
        oth = other
        if ratio > 1.0:
            C1 = C1 + rng.poisson(top)
            oth = other + rng.poisson((ratio - 1.0) * other)
        panel += C1
        sq += C1.astype(np.float64) ** 2
        out_o += oth
        sq_o += np.asarray(oth, dtype=np.float64) ** 2
    m, mo = panel / reps, out_o / reps
    return m, np.sqrt(np.maximum(sq / reps - m * m, 0.0) / reps), mo, np.sqrt(np.maximum(sq_o / reps - mo * mo, 0.0)
                                                                              / reps)


def test_the_expected_counts_match_a_monte_carlo_mean():
    """Small cells chosen so that the restore rule matters (a cell whose only count sits on a thinned gene) next to
    ordinary cells; 5 standard errors of the Monte Carlo mean per entry (1e-12 absolute where nothing varies)."""
    phi = np.array([0.25, 0.6, 1.0, 1.7, 3.0])
    C0 = np.array([[1, 0, 0, 0, 0],  # empty unless the thinning keeps the one count: Pr 0.75 of the restore
                   [2, 1, 0, 0, 0],
                   [0, 0, 3, 0, 0],  # phi = 1 only: nothing moves
                   [4, 3, 2, 1, 5],
                   [0, 0, 0, 2, 0]])
    rng = np.random.default_rng(1)
    S = rng.dirichlet(np.ones(5), size=5) * 0.4
    S[0, :] = [0.4, 0.0, 0.0, 0.0, 0.0]  # no top-up from P on the raised genes: the restore rule alone keeps it
    S[2, :] = [0.0, 0.0, 0.4, 0.0, 0.0]
    P = S / S.sum(1, keepdims=True)
    lib = np.array([4.0, 6.0, 5.0, 30.0, 4.0])
    other = lib - C0.sum(1)
    for ratio in (1.0, 2.5):
        if ratio > 1.0:
            Sp = S * phi
            Sp *= (S.sum(1) / np.maximum(Sp.sum(1), 1e-12))[:, None]
            top = Sp * ((ratio - 1.0) * lib)[:, None]
        else:
            top = np.zeros(C0.shape)
        exp_p, exp_o = tb.series_expected_counts(C0, phi, P, top, other, ratio)
        m, se, mo, se_o = _mc_render(C0, phi, P, top, other, ratio, 40_000, 3)
        assert (np.abs(m - exp_p) <= 5.0 * se + 1e-12).all(), (ratio, m - exp_p, se)
        assert (np.abs(mo - exp_o) <= 5.0 * se_o + 1e-12).all()
    # by hand, ratio 1: the one-count cell keeps its count in expectation (0.25 kept + 0.75 restored)
    exp_p, exp_o = tb.series_expected_counts(C0, phi, P, np.zeros(C0.shape), other, 1.0)
    assert exp_p[0].tolist() == [1.0, 0.0, 0.0, 0.0, 0.0]
    assert exp_p[1, 0] > 2 * 0.25 and exp_p[1, 1] > 0.6  # the restore rule raises what thinning alone gives
    assert exp_p[2].tolist() == [0.0, 0.0, 3.0, 0.0, 0.0]
    np.testing.assert_array_equal(exp_o, other)


# ------------------------------------------------------------------------------------------ X4. the deeper bracket
def test_the_deeper_brackets_cells_are_bit_for_bit_as_before(emb_expr, monkeypatch):
    X0, _, _, res, _ = _render(emb_expr, monkeypatch)
    X1, _, _, _, _ = _render(emb_expr, monkeypatch, kernel="exact_total")
    hi = res.stage == 1
    assert hi.any() and np.array_equal(X0[hi], X1[hi])
    assert not np.array_equal(X0[~hi], X1[~hi])


# ------------------------------------------------------------------------------------------ X5. cells at their totals
def test_a_cell_already_at_its_totals_is_left_as_drawn_and_zeros_stay_zero():
    expected = np.array([[2.0, 1.0, 0.0, 1.0],  # total 4, a zero the render cannot raise (gene 2)
                         [1.5, 1.5, 0.0, 3.0],  # total 6
                         [0.5, 0.25, 0.0, 0.25]])  # total 1
    expected_other = np.array([10.0, 20.2, 3.0])
    C1 = np.array([[3, 0, 0, 1],  # at its total: left as drawn
                   [0, 0, 0, 2],  # a shortfall of 4
                   [0, 2, 0, 3]])  # a surplus of 4
    oth = np.array([10.0, 25.0, 3.0])
    rng = np.random.default_rng(4)
    state = rng.bit_generator.state
    out, o, rec = tb._exact_series_totals(C1[:1], oth[:1], expected[:1], expected_other[:1], rng)
    assert np.array_equal(out, C1[:1]) and np.array_equal(o, oth[:1])
    assert rng.bit_generator.state == state  # nothing drawn
    assert rec == {"xt_share_repaired": 0.0, "xt_mean_abs_repair": 0.0}
    out, o, rec = tb._exact_series_totals(C1, oth, expected, expected_other, np.random.default_rng(4))
    assert np.array_equal(out[0], C1[0])
    np.testing.assert_array_equal(out.sum(1), np.rint(expected.sum(1)))
    np.testing.assert_array_equal(o, [10.0, 20.0, 3.0])
    assert (out[:, 2] == 0).all()  # the zero the render cannot raise stays zero
    assert (out[2] <= C1[2]).all()  # a surplus only removes drawn units
    assert (out[1] >= C1[1]).all()
    assert rec["xt_share_repaired"] == pytest.approx(2 / 3)
    assert rec["xt_mean_abs_repair"] == pytest.approx((0 + 4 + 5 + 4) / 3)
    with pytest.raises(ValueError, match="expected count 0"):
        tb._exact_series_totals(np.array([[1, 0, 1, 0]]), oth[:1], expected[:1], expected_other[:1],
                                np.random.default_rng(0))


def test_zeros_the_render_cannot_raise_stay_zero_through_the_render(emb_expr, monkeypatch):
    X, _, seen, _, _ = _render(emb_expr, monkeypatch, kernel="exact_total")
    expected = seen["expected"][0]
    assert (seen["counts"][expected == 0.0] == 0).all()
    X0, _, seen0, _, _ = _render(emb_expr, monkeypatch)
    assert (seen0["counts"][expected == 0.0] == 0).all()


# ------------------------------------------------------------------------------------------ X6. the mean library
def test_the_mean_rendered_library_is_unchanged_within_its_tolerance(emb_expr, monkeypatch):
    _, _, seen0, _, _ = _render(emb_expr, monkeypatch)
    _, _, seen1, _, _ = _render(emb_expr, monkeypatch, kernel="exact_total")
    l0, l1 = seen0["lib"], seen1["lib"]
    tol = 4.0 * l0.std() / np.sqrt(len(l0)) + 1.0
    assert abs(l1.mean() - l0.mean()) <= tol, (l1.mean(), l0.mean(), tol)
    expected, expected_other = seen1["expected"]
    assert abs(l1.mean() - (expected.sum(1) + expected_other).mean()) <= 1.0  # the two roundings


# ------------------------------------------------------------------------------------------ X7. reproducible
def test_the_exact_total_render_is_reproducible(emb_expr, monkeypatch):
    a, ra, _, _, _ = _render(emb_expr, monkeypatch, kernel="exact_total")
    b, rb, _, _, _ = _render(emb_expr, monkeypatch, kernel="exact_total")
    assert np.array_equal(a, b) and ra == rb
    c, _, _, _, _ = _render(emb_expr, monkeypatch, kernel="exact_total", repair=(9, 22))
    assert not np.array_equal(a, c)  # the repair has its own stream
    s = emb_expr.stages()
    ctx = emb_expr.ctx(inputs=s, roles=("before", "after"))
    x, y = (harness.CANDIDATES["fgw_aligned_series_xt"].build(ctx) for _ in range(2))
    assert emb_expr.digest(x) == emb_expr.digest(y)


def test_the_repair_stream_is_new():
    assert harness.SERIES_XT_STREAM not in {*range(10), harness.NBHD_STREAM, *tb.HERD_STREAMS, *tb.MIDPOINT_STREAMS}
    from vec import t2_latent

    assert harness.SERIES_XT_STREAM != t2_latent.EXPRESSION_STREAM


# ================================================================================ the harness arms
def test_the_three_arms_are_registered_for_t2_and_call_fgw_with_their_keywords(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    for name, kw in (("fgw_tissue", {"features": "tissue"}),
                     ("fgw_tissue_series", {"features": "tissue", "expression": "series"}),
                     ("fgw_aligned_series_xt", {"move_coupling": "aligned", "expression": "series_xt"})):
        c = harness.CANDIDATES[name]
        assert c.tasks == frozenset({"T2"}) and c.name == name
        assert c.build(name) == kw
    assert harness.CANDIDATES["fgw_tissue"].requires is harness.CANDIDATES["fgw"].requires
    assert harness.CANDIDATES["fgw_tissue_series"].requires is harness.CANDIDATES["fgw_series"].requires
    assert harness.CANDIDATES["fgw_aligned_series_xt"].requires is harness.CANDIDATES["fgw_aligned_series"].requires
    assert harness.CANDIDATES["fgw_tissue"].params["tissue_k"] == harness.NBHD_K == 15


def test_fgw_passes_nbhd_k_as_the_tissue_size(emb_expr, monkeypatch):
    seen = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        seen.append((kw["features"], kw["tissue_k"]))
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    c = emb_expr.ctx(inputs=emb_expr.stages(), roles=("before", "after"))
    harness.CANDIDATES["fgw"].build(c)
    harness.CANDIDATES["fgw_tissue"].build(c)
    assert seen == [("own", None), ("tissue", harness.NBHD_K)]


# ================================================================================ the runner
@pytest.fixture(scope="module")
def n43():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n43_embryo_d_arms", HARNESS / "n43_embryo_d_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_and_the_three_arms_on_the_nine_tasks(n43):
    assert n43.NEW == ARMS and n43.ARMS == (harness.FLOOR, *ARMS)
    assert not set(n43.ARMS) & {"fgw", "fgw_series", "fgw_aligned_series"}  # never rescored
    assert set(n43.NOT_RESCORED) == {"fgw", "fgw_series", "fgw_aligned_series"}
    tasks = n43.nt.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in n43.NEW:
        assert all(harness.CANDIDATES[name].applies(t) for t in tasks.values())


def test_the_runners_task_file_loads_through_nine_tasks(tmp_path, monkeypatch, n43):
    calls = []

    def fake_run_task(t, cands, seeds, data_dir):  # synthetic scores: no data, no scorer
        calls.append((t, [c.name for c in cands], list(seeds)))
        rows = [{"seed": s, "scores": {c.name: 50.0 + 0.1 * s + i for i, c in enumerate(cands)}} for s in seeds]
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "provenance": {"commit": "c" * 40, "tree_clean": True}}

    monkeypatch.setattr(n43.harness, "run_task", fake_run_task)
    out = tmp_path / "runs"
    n43.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    t = n43.nt.all_tasks()["e75"]
    assert calls == [(t, list(n43.ARMS), list(n43.SEEDS))]  # one call with the floor and the three arms
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n43_embryo_d_arms.py" and data["not_buildable"] == {}
    scores, why, meta = n43.nt._load([out])
    assert why == [] and meta == ["runs:ccccccc"]
    assert set(scores[t.key]) == set(n43.ARMS) and all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n43.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n43.main(["run", "--task", "nope", "--out", str(out)])


def _fake_runner(score: dict[str, float]):
    """run_task's shape without data or scorer: each candidate is built (its build may raise) and given a constant."""
    calls = []

    def runner(t, cands, seeds, data_dir):
        calls.append([c.name for c in cands])
        rows = []
        for s in seeds:
            row = {"seed": s, "scores": {"floor": 50.0}, "skills": {"floor": {}}, "metrics": {"floor": {}},
                   "n_cells": {"floor": 1}, "groups": {"floor": {}}, "reads": {}}
            for c in cands:
                if c.name == "floor":
                    continue
                c.build(dataclasses.make_dataclass("Ctx", ["seed"])(seed=harness.BUILD_SEED_OFFSET + s))
                row["scores"][c.name] = score[c.name]
                for k in ("skills", "metrics", "groups"):
                    row[k][c.name] = {}
                row["n_cells"][c.name], row["reads"][c.name] = 1, []
            rows.append(row)
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "candidates": [c.name for c in cands], "params": {c.name: {} for c in cands},
                "provenance": {"commit": "a" * 40, "tree_clean": True}}

    return runner, calls


def _patch_builds(monkeypatch, n43, fail_on: set[str], key: str, seed_fails: int = 0):
    """Every new arm's build replaced: the arms in fail_on raise the builder's cloud-guard refusal from seed_fails."""
    cands = harness.CANDIDATES

    def ok(ctx):
        return "built"

    def bad(ctx):
        if ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
            raise ValueError(f"{key}: cloud guard flags ['rms: RMS radius is 0.000 x the source, expected 0.436']")
        return "built"

    for x in n43.NEW:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=bad if x in fail_on else ok))


def test_an_arm_the_builder_refuses_is_recorded_and_the_other_arms_are_still_scored(tmp_path, monkeypatch, n43):
    t = n43.nt.all_tasks()["e725"]
    score = {"fgw_tissue": 53.0, "fgw_tissue_series": 0.0, "fgw_aligned_series_xt": 52.0}
    _patch_builds(monkeypatch, n43, {"fgw_tissue_series"}, t.key, seed_fails=2)
    runner, calls = _fake_runner(score)
    res = n43.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n43.NEW], *[["floor", n] for n in n43.NEW]]  # together first, then arm by arm
    assert set(res["not_buildable"]) == {"fgw_tissue_series"}
    rec = res["not_buildable"]["fgw_tissue_series"]
    assert rec["status"] == n43.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == ("ValueError: <task>: cloud guard flags ['rms: RMS radius is # x the source, "
                              "expected #']")  # the task key named as such, every number masked
    scored = {"floor", "fgw_tissue", "fgw_aligned_series_xt"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert [r["seed"] for r in res["rows"]] == list(n43.SEEDS)
    monkeypatch.setattr(n43.harness, "run_task", runner)
    out = tmp_path / "runs"
    n43.main(["run", "--task", "e725", "--out", str(out), "--data", str(tmp_path / "no_data")])
    scores, why, _ = n43.nt._load([out])
    assert why == [] and set(scores[t.key]) == scored
    (path,) = out.glob("T2e_*.json")
    assert set(json.loads(path.read_text())["not_buildable"]) == {"fgw_tissue_series"}
    # every arm failing still leaves the floor scored
    _patch_builds(monkeypatch, n43, set(n43.NEW), t.key)
    runner, calls = _fake_runner(score)
    res = n43.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n43.NEW) and calls[-1] == ["floor"]
    assert all(set(r["scores"]) == {"floor"} for r in res["rows"])
    # nothing fails: one run, nothing recorded
    _patch_builds(monkeypatch, n43, set(), t.key)
    runner, calls = _fake_runner(score)
    res = n43.score_task(t, tmp_path, runner=runner)
    assert calls == [["floor", *n43.NEW]] and res["not_buildable"] == {}


def test_a_failure_outside_a_build_stops_the_task(tmp_path, n43):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n43.score_task(n43.nt.all_tasks()["x2"], tmp_path, runner=runner)
