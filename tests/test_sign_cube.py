"""The sign cube's algebra (experiments/t3_decomp/sign_cube.py; row R72, T3_MOON5_design.md 2.2 with fixes F1 to F10) on
synthetic arrays only: the identity A + sum X = c (l_M + l_G), the projection off B, k's SD hold, the clip, the den
constant and the Newton hold of each held label's expected library, ppp as C2's own array, G0g's positions and the
vertices an axis held at '+' leaves. No competition file is read."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t3_decomp"))
import sign_cube as sc  # noqa: E402

LN4 = np.log(4.0)


def _world(seed: int = 0, G: int = 400, L: int = 10, N: int = 2500):
    rng = np.random.default_rng(seed)
    lab = rng.integers(0, L, N)
    lab[:4] = L  # a label under 1% of the cells: never held
    prof = rng.dirichlet(np.ones(G) * 0.3, size=L + 1)
    C = np.vstack([rng.poisson(prof[t] * 150) for t in lab])
    # tables shaped like the record's: l_M mostly a stage offset (N8, K18 [9]), the WT poles around it
    l_m, l_g = rng.normal(0, 0.5, G), rng.normal(0, 0.4, G)
    raw = {"l_M": l_m, "l_G": l_g, "S_M": 0.5 * l_m + rng.normal(0, 0.05, G),
           "Y_M": -0.5 * l_m + rng.normal(0, 0.05, G), "S_G": 0.3 * l_g + rng.normal(0, 0.05, G),
           "Y_G": -0.3 * l_g + rng.normal(0, 0.05, G)}
    tab = {k: np.clip(v, -LN4, LN4) for k, v in raw.items()}
    c = sc.sd(tab["l_M"]) / sc.sd(tab["l_M"] + tab["l_G"])
    phi_m, phi_g = np.exp(tab["l_M"]), np.exp(tab["l_G"])
    phi_c2 = (phi_m ** c) * (phi_g ** c)  # build_two_knockout.combine('c2', 1.0)'s operations
    held, M = sc.label_rows(C, lab)
    den = C.sum(0) / C.sum()
    B = sc.basis(M, phi_c2)
    Q = sc.projector(B)
    ax = sc.cube_axes(tab["l_M"], tab["l_G"], tab["S_M"], tab["Y_M"], tab["S_G"], tab["Y_G"], c)
    q = {j: Q(ax[j]) for j in sc.AXES}
    return dict(tab=tab, c=c, phi_c2=phi_c2, held=held, M=M, den=den, B=B, Q=Q, ax=ax, q=q, C=C, lab=lab)


def test_the_anchor_and_axes_sum_to_c2s_log_factor():
    w = _world()
    t, c, ax = w["tab"], w["c"], w["ax"]
    tot = ax["A"] + ax["S"] + ax["M"] + ax["G"]
    assert np.abs(tot - c * (t["l_M"] + t["l_G"])).max() <= 1e-12
    assert np.abs(tot - np.log(w["phi_c2"])).max() <= 1e-12
    an = sc.analogue_axes(t["l_M"], t["S_M"], t["Y_M"])
    assert np.abs(an["A"] + an["S"] + an["M"] - t["l_M"]).max() <= 1e-12  # G1's base is l_M


def test_a_planted_den_deviation_cancels_in_every_axis_but_the_anchor():
    """The reference's own deviation (a shift d of ln den in every table) cancels in X_S, X_M and X_G and stays in A
    (design 2.2: the den's log cancels in every axis, off the clip)."""
    w = _world(1)
    t, c = w["tab"], w["c"]
    d = np.zeros_like(t["l_M"])
    d[7] = -0.3  # every table is ln(num) - ln(den); a den deviation enters each with a minus sign
    shifted = {k: v + d for k, v in t.items()}
    a0 = sc.cube_axes(*(t[k] for k in ("l_M", "l_G", "S_M", "Y_M", "S_G", "Y_G")), c)
    a1 = sc.cube_axes(*(shifted[k] for k in ("l_M", "l_G", "S_M", "Y_M", "S_G", "Y_G")), c)
    for j in sc.AXES:
        assert np.abs(a1[j] - a0[j]).max() <= 1e-12
    assert a1["A"][7] - a0["A"][7] == pytest.approx(2 * c * -0.3, abs=1e-12)


@pytest.mark.parametrize("which", ["E9.5", "E8.0"])
def test_f3_planted_gene_moves_the_axes_by_the_fixed_amounts(which):
    """F3, in the axes' algebra: a gene doubled in WT E9.5 only moves S_M and S_G by ln 2, so X_S +c ln2/2, X_M -c ln2,
    X_G -c ln2/2, A +c ln2; halved in WT E8.0 only moves Y_M and Y_G by -ln 2, so X_S +c ln2/2, X_G +c ln2/2,
    A -c ln2 (X_M unchanged). The tables' own response to the planted gene is tested in test_cube_tables.py."""
    w = _world(2)
    t, c, g = dict(w["tab"]), w["c"], 11
    a0 = sc.cube_axes(*(t[k] for k in ("l_M", "l_G", "S_M", "Y_M", "S_G", "Y_G")), c)
    if which == "E9.5":
        for k in ("S_M", "S_G"):
            t[k] = t[k].copy()
            t[k][g] += np.log(2)
        want = {"S": c * np.log(2) / 2, "M": -c * np.log(2), "G": -c * np.log(2) / 2, "A": c * np.log(2)}
    else:
        for k in ("Y_M", "Y_G"):
            t[k] = t[k].copy()
            t[k][g] -= np.log(2)
        want = {"S": c * np.log(2) / 2, "M": 0.0, "G": c * np.log(2) / 2, "A": -c * np.log(2)}
    a1 = sc.cube_axes(*(t[k] for k in ("l_M", "l_G", "S_M", "Y_M", "S_G", "Y_G")), c)
    for k, v in want.items():
        assert a1[k][g] - a0[k][g] == pytest.approx(v, abs=1e-12)
        others = np.delete(a1[k] - a0[k], g)
        assert np.abs(others).max() <= 1e-12


def test_the_projection_is_off_b_and_idempotent():
    w = _world(3)
    for j in sc.AXES:
        q = w["q"][j]
        assert np.abs(w["B"].T @ q).max() <= 1e-12 * max(1.0, np.abs(w["ax"][j]).max())
        assert np.abs(w["Q"](q) - q).max() <= 1e-12


def test_label_rows_hold_labels_with_one_percent_of_the_cells():
    C = np.array([[1, 3], [2, 2], [0, 4]] + [[5, 5]] * 197)
    lab = np.array([0, 0, 1] + [2] * 197)
    held, M = sc.label_rows(C, lab)
    assert list(held) == [0, 2]  # label 1 holds 1 of 200 cells, under 1%
    assert np.allclose(M[0], [3 / 8, 5 / 8]) and np.allclose(M.sum(1), 1)


def test_ppp_is_c2s_own_array_with_k_one():
    w = _world(4)
    phi, rec = sc.vertex_factor(w["phi_c2"], w["q"], "ppp", w["M"], w["den"], w["B"])
    assert rec["k"] == 1.0 and rec["coefficients"] == [0.0, 0.0, 0.0] and rec["hold"]["rounds"] == 0
    assert np.array_equal(phi, w["phi_c2"]) and phi.tobytes() == w["phi_c2"].tobytes()


@pytest.mark.parametrize("v", sc.NAMED)
def test_each_vertex_holds_size_den_and_every_held_label_library(v):
    w = _world(5)
    phi, rec = sc.vertex_factor(w["phi_c2"], w["q"], v, w["M"], w["den"], w["B"])
    lb = np.log(w["phi_c2"])
    assert rec["sd_ln_before_clip"] == pytest.approx(sc.sd(lb), rel=1e-9)  # k holds SD(ln phi) before the clip
    assert sc.K_RANGE[0] <= rec["k"] <= sc.K_RANGE[1]
    lib = (w["M"] @ phi) / (w["M"] @ w["phi_c2"]) - 1
    assert np.abs(lib).max() <= sc.LIB_TOL and sc.within_band(rec)
    assert abs(w["den"] @ phi / (w["den"] @ w["phi_c2"]) - 1) <= sc.DEN_TOL
    lo, hi = np.minimum(0.25, w["phi_c2"]), np.maximum(4.0, w["phi_c2"])
    # the clip, then the exact den constant: a clipped gene may sit off its bound by that one constant only
    excess = max(float((lo / phi).max()), float((phi / hi).max()))
    assert excess <= 1.0 + 1e-3
    assert not np.array_equal(phi, w["phi_c2"])
    flipped = [j for j, ch in zip(sc.AXES, v, strict=True) if ch == "n"]
    # the move lies along the flipped axes (first order): its correlation with -sum of the flipped Q(X_j) is high
    move = np.log(phi) - lb
    want = -np.sum([w["q"][j] for j in flipped], axis=0)
    assert np.corrcoef(move, want)[0, 1] > 0.8


def test_k_is_the_positive_root_nearest_one_and_refused_outside_the_range():
    rng = np.random.default_rng(6)
    b = rng.normal(0, 0.5, 300)
    q = [0.4 * b + rng.normal(0, 0.05, 300), 0.3 * b + rng.normal(0, 0.05, 300)]
    k, rec = sc.solve_k(b, q, (-1, 1), k_range=(1e-9, 1e9))
    lnphi = b + sum((k * s - 1) * x for s, x in zip((-1, 1), q, strict=True))
    assert sc.sd(lnphi) == pytest.approx(sc.sd(b), rel=1e-10) and k > 0
    assert k == min([r for r in rec["roots"] if r > 0], key=lambda r: abs(r - 1))
    k1, _ = sc.solve_k(b, q, (1, 1))
    assert k1 == 1.0  # the identity's root, snapped
    with pytest.raises(sc.CubeRefused, match="outside"):
        sc.solve_k(b, q, (-1, 1), k_range=(k + 0.01, k + 1.0))
    with pytest.raises(sc.CubeRefused):
        sc.solve_k(b, [np.zeros(300), np.zeros(300)], (-1, 1))
    with pytest.raises(sc.CubeRefused, match="no real k"):  # axes orthogonal to the base and too large to hold
        sc.solve_k(b, [rng.normal(0, 0.5, 300), rng.normal(0, 0.5, 300)], (-1, 1))


def test_a_hold_that_cannot_converge_refuses_the_vertex():
    w = _world(7)
    with pytest.raises(sc.CubeRefused, match="did not converge"):
        sc.flip_factor(w["phi_c2"], [w["q"][j] for j in sc.AXES], (-1, -1, -1), w["M"], w["den"], w["B"], rounds=0)


def test_signs_and_vertex_names():
    assert sc.VERTICES[0] == "ppp" and len(sc.VERTICES) == 8 and len(sc.NAMED) == 7 and "ppp" not in sc.NAMED
    assert sc.signs("pnp") == (1, -1, 1)
    for bad in ("", "pxp", 3, "PPN"):
        with pytest.raises(sc.CubeRefused):
            sc.signs(bad)
    w = _world(8)
    with pytest.raises(sc.CubeRefused):
        sc.vertex_factor(w["phi_c2"], w["q"], "pn", w["M"], w["den"], w["B"])


def test_kept_vertices_never_reflect_a_held_axis():
    assert sc.kept_vertices(["S", "M", "G"]) == list(sc.NAMED)
    assert sc.kept_vertices(["M", "G"]) == ["ppn", "pnp", "pnn"]
    assert sc.kept_vertices(["S"]) == ["npp"]
    assert sc.kept_vertices([]) == []


def test_g0g_positions_and_bar():
    young, old = np.zeros(5), np.ones(5)
    tabs = {"E8.0": young, "E9.5": old, "E8.25_late": np.full(5, 0.2), "E8.5": np.full(5, 0.35),
            "E8.75": np.full(5, 0.5)}
    pos = sc.stage_positions(tabs)
    assert pos["E8.25_late"] == pytest.approx(0.2) and pos["E8.75"] == pytest.approx(0.5)
    assert sc.g0g_passes(pos)
    assert not sc.g0g_passes({**pos, "E8.5": 0.6})  # E8.5 not between E8.25_late and E8.75
    assert not sc.g0g_passes({**pos, "E8.75": 0.8})  # |s(E8.75) - 0.5| > 0.25
    assert not sc.g0g_passes({**pos, "E8.25_late": -0.1})
    assert sc.g0g_passes({k: v for k, v in pos.items() if k != "E8.5"})
    with pytest.raises(sc.CubeRefused):
        sc.stage_positions({"E8.0": young, "E9.5": young})
