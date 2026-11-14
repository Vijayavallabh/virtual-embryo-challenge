"""The embryo batch F (campaign branch campaign/emb-f): a zero-temperature Potts sorting of which drawn cell sits where
(vec.t2_build.potts_arrangement, the arm fgw_potts) and a partial FGW coupling whose transported share is the overlap of
the brackets' anchors' label distributions (vec.t2_build.partial_fgw_coupling, interpolate_clouds(...,
coupling="partial"), the arm fgw_partial), and the runner experiments/harness/n45_embryo_f_arms.py. Synthetic clouds and
stages only; no data file, committed result or scorer is read.

Tolerances stated here: the partial coupling's total mass against s and its marginals against the uniform masses,
1e-9 absolute (POT's exact network simplex; its own convergence threshold is tol=1e-8 on the loss, not on the mass)."""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
HARNESS = REPO / "experiments" / "harness"
BASE = "f29a3eb"  # campaign/emb-d, the branch's base
ARMS = ("fgw_potts", "fgw_partial")
K = harness.NBHD_K


def _helpers(name: str):
    spec = importlib.util.spec_from_file_location(f"{name}_helpers_f", REPO / "tests" / f"{name}.py")
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


def _rows_multiset(X: np.ndarray) -> np.ndarray:
    X = np.ascontiguousarray(X)
    return X[np.lexsort(X.T[::-1])]


# ================================================================================ the Potts sorting
def _cloud(n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).uniform(0.0, 10.0, size=(n, 3))


def _mixed(n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).choice(np.array(["A", "B"]), n)


def _sorted(C: np.ndarray) -> np.ndarray:
    return np.where(C[:, 0] < np.median(C[:, 0]), "A", "B")


# ------------------------------------------------------------------------------------------ P1. the statistic by hand
def test_the_neighbour_lists_and_the_share_match_a_brute_force_count():
    C = _cloud(60, 1)
    lab = _mixed(60, 2)
    nb = tb.knn_others(C, K)
    assert nb.shape == (60, K - 1)
    D = np.linalg.norm(C[:, None] - C[None], axis=2)
    np.fill_diagonal(D, np.inf)
    by_hand = np.argsort(D, axis=1, kind="stable")[:, :K - 1]
    assert all(set(nb[i]) == set(by_hand[i]) for i in range(60))
    share = float(np.mean([(lab[by_hand[i]] == lab[i]).mean() for i in range(60)]))
    assert tb.same_group_share(C, lab, K) == pytest.approx(share, abs=1e-15)
    # fewer than K cells: every other cell
    assert tb.knn_others(C[:9], K).shape == (9, 8)
    with pytest.raises(ValueError, match="two cells"):
        tb.same_group_share(C[:1], lab[:1], K)


# ------------------------------------------------------------------------------------------ P2. the incremental H
def test_the_incremental_h_equals_a_full_recomputation_after_every_accepted_swap():
    C, lab = _cloud(80, 3), _mixed(80, 4)
    trace: list = []
    at, rec = tb.potts_sort(C, lab, 0.9, K, np.random.default_rng(5), trace=trace)
    assert len(trace) == rec["accepted_swaps"] > 5
    for at_k, h in trace:
        assert tb.same_group_share(C, lab[at_k], K) == h  # exact: both are the same integer count over n (K - 1)
    assert np.array_equal(trace[-1][0], at) and rec["H_final"] == trace[-1][1]
    assert sorted(at.tolist()) == list(range(80))


# ------------------------------------------------------------------------------------------ P3. rises and stops
def test_a_mixed_cloud_rises_towards_the_brackets_target_and_stops_at_or_just_past_it():
    n = 200
    C, lab = _cloud(n, 6), _mixed(n, 7)
    Cb, Ca = _cloud(300, 8), _cloud(150, 9)
    brackets = ((Cb, _mixed(300, 22)), (Ca, _sorted(Ca)))  # the later bracket sorted in two halves
    at, rec = tb.potts_arrangement(C, lab, brackets, 0.1, K, np.random.default_rng(10), np.random.default_rng(11))
    assert rec["H_fgw"] < rec["H_target"]
    assert rec["stop"] == "reached" and rec["H_final"] >= rec["H_target"]
    trace: list = []
    at2, rec2 = tb.potts_sort(C, lab, rec["H_target"], K, np.random.default_rng(11), trace=trace)
    assert np.array_equal(at, at2) and rec2["H_final"] == rec["H_final"]
    hs = [rec["H_fgw"], *[h for _, h in trace]]
    assert all(b > a for a, b in zip(hs, hs[1:], strict=False))  # every accepted swap raises H
    assert hs[-2] < rec["H_target"] <= hs[-1]  # the last swap reaches or crosses the target
    assert abs(hs[-1] - rec["H_target"]) < abs(hs[-2] - rec["H_target"])  # and ends closer to it than before
    assert rec["H_target"] == 0.9 * rec["H_before"] + 0.1 * rec["H_after"]
    assert set(rec) == {"H_fgw", "H_before", "H_after", "H_target", "H_final", "accepted_swaps", "stop"}


def test_a_target_far_above_a_mixed_cloud_freezes_short_of_it():
    """The zero-temperature rule is a greedy descent: from a fully mixed cloud with both brackets sorted (target near
    0.9) it freezes in a local optimum well below the target and ends on the stall rule, H still risen at every swap.
    Recorded here as a property of the specified rule, not a fault of the build."""
    n = 200
    C, lab = _cloud(n, 6), _mixed(n, 7)
    Cb, Ca = _cloud(300, 8), _cloud(150, 9)
    brackets = ((Cb, _sorted(Cb)), (Ca, _sorted(Ca)))
    at, rec = tb.potts_arrangement(C, lab, brackets, 0.5, K, np.random.default_rng(10), np.random.default_rng(11))
    assert rec["stop"] == "stalled" and rec["H_fgw"] < rec["H_final"] < rec["H_target"]
    assert rec["H_target"] - rec["H_final"] > 0.1


def test_next_to_the_target_the_strictly_closer_rule_can_end_on_the_stall_within_one_swap():
    """H moves in steps of whole same-group pairs over n (K - 1): when every swap left would land farther past the
    target than H sits below it, none is accepted and the run ends on the stall rule just short of the target (the
    same file and brackets as the rising test with t = 0.25, and a swap stream that meets this case)."""
    n = 200
    C, lab = _cloud(n, 6), _mixed(n, 7)
    Cb, Ca = _cloud(300, 8), _cloud(150, 9)
    brackets = ((Cb, _mixed(300, 22)), (Ca, _sorted(Ca)))
    _, rec = tb.potts_arrangement(C, lab, brackets, 0.25, K, np.random.default_rng(10), np.random.default_rng(11))
    assert rec["stop"] == "stalled" and rec["H_final"] < rec["H_target"]
    assert rec["H_target"] - rec["H_final"] < 2.0 / (n * (K - 1))  # within less than one swap's step of two pairs


def test_a_cloud_sorted_beyond_the_target_falls_towards_it():
    n = 200
    C = _cloud(n, 12)
    lab = _sorted(C)
    Cb, Ca = _cloud(300, 13), _cloud(300, 14)
    brackets = ((Cb, _mixed(300, 15)), (Ca, _mixed(300, 16)))
    at, rec = tb.potts_arrangement(C, lab, brackets, 0.25, K, np.random.default_rng(1), np.random.default_rng(2))
    assert rec["H_fgw"] > rec["H_target"]
    assert rec["stop"] == "reached" and rec["H_final"] <= rec["H_target"]
    trace: list = []
    tb.potts_sort(C, lab, rec["H_target"], K, np.random.default_rng(2), trace=trace)
    hs = [rec["H_fgw"], *[h for _, h in trace]]
    assert all(b < a for a, b in zip(hs, hs[1:], strict=False))


# ------------------------------------------------------------------------------------------ P4. idle cases
def test_at_the_target_nothing_is_drawn_and_nothing_moves():
    C, lab = _cloud(70, 17), _mixed(70, 18)
    h = tb.same_group_share(C, lab, K)
    rng = np.random.default_rng(3)
    state = rng.bit_generator.state
    at, rec = tb.potts_sort(C, lab, h, K, rng)
    assert np.array_equal(at, np.arange(70)) and rec["stop"] == "at_target" and rec["accepted_swaps"] == 0
    assert rng.bit_generator.state == state


def test_one_group_cannot_move_and_stalls_after_n_proposals():
    C = _cloud(50, 19)
    lab = np.full(50, "A")
    rng = np.random.default_rng(4)
    at, rec = tb.potts_sort(C, lab, 0.3, K, rng)
    assert np.array_equal(at, np.arange(50)) and rec["stop"] == "stalled" and rec["H_final"] == 1.0
    ref = np.random.default_rng(4)  # exactly n proposals: one block of n positions and n neighbour ranks
    ref.integers(0, 50, size=50)
    ref.integers(0, K - 1, size=50)
    assert rng.bit_generator.state == ref.bit_generator.state


def test_a_file_whose_cells_all_share_one_group_is_fgw_bit_for_bit(emb_expr):
    s = emb_expr.stages()
    for a in s.values():
        a.obs["group"] = "Mesoderm"
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    rec: dict = {}
    arm = harness._fgw_potts(c, build_record=rec)
    assert emb_expr.digest(arm) == emb_expr.digest(harness.CANDIDATES["fgw"].build(c))
    assert rec["potts"]["accepted_swaps"] == 0 and rec["potts"]["H_fgw"] == 1.0
    # with the brackets mixed the target is below 1, and the one-group file stalls without a swap
    s2 = emb_expr.stages()
    c2 = emb_expr.ctx(inputs=s2, roles=("before", "after"))
    real = tb.drawn_labels
    rec2: dict = {}
    try:
        tb.drawn_labels = lambda res, b, a: np.full(len(res.rows), "one")  # the file's cells all in one group
        arm2 = harness._fgw_potts(c2, build_record=rec2)
    finally:
        tb.drawn_labels = real
    assert rec2["potts"]["stop"] == "stalled" and rec2["potts"]["H_target"] < 1.0
    assert emb_expr.digest(arm2) == emb_expr.digest(harness.CANDIDATES["fgw"].build(c2))


def test_with_h_fgw_equal_to_the_target_the_arm_is_fgw_bit_for_bit(emb_expr, monkeypatch):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    fgw = harness.CANDIDATES["fgw"].build(c)
    first: dict = {}
    harness._fgw_potts(c, build_record=first)
    h = first["potts"]["H_fgw"]
    t = (c.task.target_day - c.task.role_days["before"]) / (c.task.role_days["after"] - c.task.role_days["before"])
    assert (1.0 - t) * h + t * h == h  # the precondition: both brackets at H(fgw) give the target H(fgw) exactly
    monkeypatch.setattr(tb, "bracket_share", lambda *a, **kw: h)
    rec: dict = {}
    arm = harness._fgw_potts(c, build_record=rec)
    assert rec["potts"]["stop"] == "at_target" and rec["potts"]["H_target"] == h
    assert emb_expr.digest(arm) == emb_expr.digest(fgw)


# ------------------------------------------------------------------------------------------ P5. the arm
def test_the_arm_keeps_fgws_cells_rows_and_positions_and_moves_only_the_assignment(emb_expr):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    fgw = harness.CANDIDATES["fgw"].build(c)
    rec: dict = {}
    arm = harness._fgw_potts(c, build_record=rec)
    again = harness.CANDIDATES["fgw_potts"].build(c)
    assert np.array_equal(arm.coords, fgw.coords)  # the positions array, bit for bit
    assert arm.n == fgw.n and arm.X.shape == fgw.X.shape
    assert np.array_equal(_rows_multiset(arm.X), _rows_multiset(fgw.X))  # a permutation of fgw's expression rows
    assert not np.array_equal(arm.X, fgw.X)  # the idea is not idle here
    assert emb_expr.digest(arm) == emb_expr.digest(again)  # same seed, same output
    # the permutation is the sorting's own: fgw's rows taken in the order `at`
    res, Eb, Ea, frac = emb_expr._interp(c)
    at, rec2 = harness.potts_order(res, s["before"], s["after"], frac, c.seed)
    assert np.array_equal(arm.X, fgw.X[at]) and rec2 == rec["potts"]
    assert rec["potts"]["accepted_swaps"] > 0 and rec["potts"]["stop"] in tb.POTTS_STOPS
    other = harness._fgw_potts(dataclasses.replace(c, seed=c.seed + 1))
    assert not np.array_equal(other.X, arm.X)


def test_the_arm_refuses_a_file_without_group_labels(emb_expr):
    s = emb_expr.stages()
    del s["after"].obs["group"]
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    with pytest.raises(ValueError, match="group"):
        harness.CANDIDATES["fgw_potts"].build(c)
    with pytest.raises(ValueError, match="local_perm"):
        harness._fgw(c, potts=True, local_perm=True)


def test_a_bracket_with_fewer_cells_than_the_file_uses_all_of_them():
    Cb, lab = _cloud(40, 20), _mixed(40, 21)
    rng = np.random.default_rng(6)
    state = rng.bit_generator.state
    assert tb.bracket_share(Cb, lab, 40, K, rng) == tb.same_group_share(Cb, lab, K)
    assert tb.bracket_share(Cb, lab, 100, K, rng) == tb.same_group_share(Cb, lab, K)
    assert rng.bit_generator.state == state  # no draw
    rows = np.random.default_rng(6).choice(40, 25, replace=False)
    assert tb.bracket_share(Cb, lab, 25, K, rng) == tb.same_group_share(Cb[rows], lab[rows], K)


# ------------------------------------------------------------------------------------------ P6. the streams
def test_the_two_potts_streams_are_new_and_differ_from_each_other():
    used = {*range(13), harness.NBHD_STREAM, harness.SERIES_XT_STREAM, *range(20, 26), *tb.HERD_STREAMS,
            *tb.MIDPOINT_STREAMS}
    from vec import t2_latent

    used.add(t2_latent.EXPRESSION_STREAM)
    a, b = harness.POTTS_STREAMS
    assert a != b and not {a, b} & used
    seed = harness.BUILD_SEED_OFFSET + 4
    draws = {s: np.random.default_rng([seed, s]).random(64) for s in sorted(used | {a, b})}
    for s in sorted(used):
        assert not np.array_equal(draws[a], draws[s]) and not np.array_equal(draws[b], draws[s])
    assert not np.array_equal(draws[a], draws[b])


def test_the_build_draws_on_the_two_streams_and_nothing_else_changes(emb_expr, monkeypatch):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    seen = []
    real = tb.potts_arrangement

    def spy(coords, labels, brackets, t, k, bracket_rng, swap_rng):
        seen.append((bracket_rng.random(4), swap_rng.random(4), k))
        return real(coords, labels, brackets, t, k, bracket_rng, swap_rng)

    monkeypatch.setattr(tb, "potts_arrangement", spy)
    harness.CANDIDATES["fgw_potts"].build(c)
    (br, sw, k), = seen
    assert k == harness.NBHD_K
    assert np.array_equal(br, np.random.default_rng([c.seed, harness.POTTS_STREAMS[0]]).random(4))
    assert np.array_equal(sw, np.random.default_rng([c.seed, harness.POTTS_STREAMS[1]]).random(4))


# ================================================================================ the partial coupling
def _labels(A: np.ndarray, B: np.ndarray, flip: int = 120) -> tuple[np.ndarray, np.ndarray]:
    la = np.where(A[:, 0] > np.median(A[:, 0]), "p", "q")
    lb = np.where(B[:, 0] > np.median(B[:, 0]), "p", "q")
    lb[:flip] = "r"
    return la, lb


# ------------------------------------------------------------------------------------------ Q1. s by hand
def test_the_overlap_by_hand():
    s = tb.label_overlap(["a", "a", "b", "c"], ["a", "b", "b", "b", "d"])
    # p_b: a 1/2, b 1/4, c 1/4; p_a: a 1/5, b 3/5, d 1/5; min: a 1/5, b 1/4
    assert s == Fraction(1, 5) + Fraction(1, 4) == Fraction(9, 20)
    assert tb.label_overlap(["x", "y", "y"], ["y", "x", "y"]) == 1
    assert tb.label_overlap(["x"], ["y"]) == 0


# ------------------------------------------------------------------------------------------ Q2. unset and s = 1
def test_coupling_unset_is_the_entry_bit_for_bit(tmp_path, polar):
    A, B, E, EB = polar._pair()
    base = polar._base_module(tmp_path, BASE)
    old = base.interpolate_clouds(A, B, E, EB, 0.4, 300, **polar.KW)
    new = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, **polar.KW)
    fgw = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, coupling="fgw", **polar.KW)
    polar._same(old, new)
    polar._same(new, fgw)


def test_an_overlap_of_exactly_one_is_the_entry_bit_for_bit(polar, monkeypatch):
    A, B, E, EB = polar._pair()
    entry = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, **polar.KW)

    def never(*a, **kw):
        raise AssertionError("the partial solver is not called at s = 1")

    monkeypatch.setattr(tb, "partial_fgw_coupling", never)
    one = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, coupling="partial",
                                coupling_labels=(np.full(len(A), "x"), np.full(len(B), "x")), **polar.KW)
    polar._same(dataclasses.replace(one, record={k: v for k, v in one.record.items()
                                                 if k not in ("coupling", "partial")}), entry)
    assert one.record["partial"] == {"s": 1.0, "transported_mass": None, "no_mass_share": [0.0, 0.0]}


def test_the_arm_with_equal_label_shares_is_fgw_bit_for_bit(emb_expr):
    c = emb_expr.ctx(inputs=emb_expr.stages(), roles=("before", "after"))
    rec: dict = {}
    arm = harness._fgw_partial(c, build_record=rec)
    assert rec["partial"]["s"] == 1.0  # the synthetic brackets hold each label in the same share, every cell an anchor
    assert emb_expr.digest(arm) == emb_expr.digest(harness.CANDIDATES["fgw"].build(c))


# ------------------------------------------------------------------------------------------ Q3. the solver's call
def test_the_solver_gets_the_entrys_costs_alpha_and_cap_and_mass_s(polar, monkeypatch):
    import ot

    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    calls: dict = {"fgw": [], "pfgw": []}
    real_fgw, real_p = ot.gromov.fused_gromov_wasserstein, ot.gromov.partial_fused_gromov_wasserstein

    def fgw(*a, **kw):
        calls["fgw"].append((a, kw))
        return real_fgw(*a, **kw)

    def pfgw(*a, **kw):
        calls["pfgw"].append((a, kw))
        return real_p(*a, **kw)

    monkeypatch.setattr(ot.gromov, "fused_gromov_wasserstein", fgw)
    monkeypatch.setattr(ot.gromov, "partial_fused_gromov_wasserstein", pfgw)
    tb.interpolate_clouds(A, B, E, EB, 0.4, 300, **polar.KW)
    res = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, coupling="partial", coupling_labels=(la, lb),
                                with_coupling=True, **polar.KW)
    (ea, ekw), = calls["fgw"]
    (pa, pkw), = calls["pfgw"]
    for x, y in zip(ea[:5], pa[:5], strict=True):  # M, C1, C2, p, q: the entry's
        assert np.array_equal(x, y)
    ib, ia = res.coupling["anchors_before"], res.coupling["anchors_after"]
    s = tb.label_overlap(la[ib], lb[ia])
    assert 0 < s < 1
    assert pkw == {"m": float(s), "alpha": tb.INTERP_ALPHA, "numItermax": tb.INTERP_EMD_ITER}
    assert ekw["alpha"] == pkw["alpha"] and ekw["loss_fun"] == "square_loss"
    assert res.record["coupling"] == "partial" and res.record["partial"]["s"] == float(s)


# ------------------------------------------------------------------------------------------ Q4. mass and marginals
def test_below_one_the_coupling_carries_mass_s_within_the_uniform_masses(polar, monkeypatch):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    rot = []
    real_rot = tb.coupling_rotation
    monkeypatch.setattr(tb, "coupling_rotation", lambda a, b, T: rot.append(np.array(T)) or real_rot(a, b, T))
    res = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, coupling="partial", coupling_labels=(la, lb),
                                with_coupling=True, **polar.KW)
    T = res.coupling["fgw"]
    s = res.record["partial"]["s"]
    n_b, n_a = T.shape
    assert abs(T.sum() - s) <= 1e-9 and res.record["partial"]["transported_mass"] == float(T.sum())
    assert (T >= 0).all() and (T.sum(1) <= 1.0 / n_b + 1e-9).all() and (T.sum(0) <= 1.0 / n_a + 1e-9).all()
    assert res.record["partial"]["no_mass_share"] == [float((T.sum(1) <= 0).mean()), float((T.sum(0) <= 0).mean())]
    (R_in,) = rot  # the rotation is coupling_rotation on this coupling
    assert np.array_equal(R_in, T)
    entry = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, **polar.KW)
    assert np.array_equal(entry.rows, res.rows) and np.array_equal(entry.stage, res.stage)  # the same draws
    assert not np.array_equal(entry.coords, res.coords)  # the idea is not idle here


# ------------------------------------------------------------------------------------------ Q5. the displacements
def test_an_anchor_without_mass_stays_and_a_partly_matched_one_uses_its_own_mass():
    T = np.array([[0.0, 0.0, 0.0],
                  [0.1, 0.3, 0.0],
                  [0.0, 0.0, 0.05]])
    own = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    other = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    out = tb.partial_partners(T, own, other)
    assert np.array_equal(out[0], own[0])  # displacement exactly 0
    np.testing.assert_allclose(out[1], (0.1 * other[0] + 0.3 * other[1]) / 0.4, rtol=1e-15)  # not / (1 / 3)
    np.testing.assert_allclose(out[2], other[2], rtol=1e-15)
    cols = tb.partial_partners(T.T, other, own)  # the after side: column 2 of T.T, i.e. row 2 of T, partly matched
    np.testing.assert_allclose(cols[0], own[1], rtol=1e-15)  # column 0 of T: only row 1 holds mass


def test_in_the_build_unmatched_anchors_do_not_move(polar, monkeypatch):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B, flip=200)
    seen = []
    real = tb.partial_partners

    def spy(T, own, other):
        out = real(T, own, other)
        seen.append((np.array(T), np.array(own), out))
        return out

    monkeypatch.setattr(tb, "partial_partners", spy)
    res = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, coupling="partial", coupling_labels=(la, lb), **polar.KW)
    assert len(seen) == 2 and max(res.record["partial"]["no_mass_share"]) > 0.0
    for T, own, out in seen:
        none = T.sum(1) <= 0.0
        assert np.array_equal(out[none] - own[none], np.zeros_like(own[none]))


# ------------------------------------------------------------------------------------------ Q6. refusals
def test_a_solver_failure_or_a_non_finite_coupling_refuses(polar, monkeypatch, emb_expr):
    import ot

    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)

    def boom(*a, **kw):
        raise RuntimeError("no convergence")

    monkeypatch.setattr(ot.gromov, "partial_fused_gromov_wasserstein", boom)
    with pytest.raises(ValueError, match="partial FGW solver failed"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 300, coupling="partial", coupling_labels=(la, lb), **polar.KW)
    monkeypatch.setattr(ot.gromov, "partial_fused_gromov_wasserstein",
                        lambda M, *a, **kw: np.full(M.shape, np.nan))
    with pytest.raises(ValueError, match="not finite"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 300, coupling="partial", coupling_labels=(la, lb), **polar.KW)
    # the arm refuses too, and does not fall back to the entry
    s = emb_expr.stages()
    ct = s["after"].obs["celltype"].astype(str).to_numpy().copy()
    ct[:50] = "PS"
    s["after"].obs["celltype"] = ct
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    monkeypatch.setattr(ot.gromov, "partial_fused_gromov_wasserstein", boom)
    with pytest.raises(ValueError, match="partial FGW solver failed"):
        harness.CANDIDATES["fgw_partial"].build(c)


def test_the_coupling_arguments_are_checked(polar):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    with pytest.raises(ValueError, match="coupling must be"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, coupling="unbalanced", **polar.KW)
    with pytest.raises(ValueError, match="coupling_labels"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, coupling="partial", **polar.KW)
    with pytest.raises(ValueError, match="coupling_labels"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, coupling_labels=(la, lb), **polar.KW)
    with pytest.raises(ValueError, match="one label per cell"):
        tb.interpolate_clouds(A, B, E, EB, 0.4, 30, coupling="partial", coupling_labels=(la[:5], lb), **polar.KW)
    for kw in ({"move_coupling": "aligned"}, {"with_partners": True}):
        with pytest.raises(ValueError, match="entry's coupling only"):
            tb.interpolate_clouds(A, B, E, EB, 0.4, 30, coupling="partial", coupling_labels=(la, lb), **kw,
                                  **polar.KW)


# ------------------------------------------------------------------------------------------ Q7. reproducible
def test_the_partial_build_is_reproducible(polar, emb_expr):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    kw = dict(coupling="partial", coupling_labels=(la, lb))
    a = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, **kw, **polar.KW)
    b = tb.interpolate_clouds(A, B, E, EB, 0.4, 300, **kw, **polar.KW)
    polar._same(a, b)
    s = emb_expr.stages()
    ct = s["after"].obs["celltype"].astype(str).to_numpy().copy()
    ct[:50] = "PS"
    s["after"].obs["celltype"] = ct
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    rec: dict = {}
    x = harness._fgw_partial(c, build_record=rec)
    y = harness.CANDIDATES["fgw_partial"].build(c)
    f = harness.CANDIDATES["fgw"].build(c)
    assert emb_expr.digest(x) == emb_expr.digest(y)  # same seed, same output
    assert 0.0 < rec["partial"]["s"] < 1.0 and abs(rec["partial"]["transported_mass"] - rec["partial"]["s"]) <= 1e-9
    assert np.array_equal(x.X, f.X) and not np.array_equal(x.coords, f.coords)  # only the positions move


# ================================================================================ the harness arms
def test_the_two_arms_are_registered_for_t2_and_call_fgw_with_their_keywords(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    for name, kw in (("fgw_potts", {"potts": True, "build_record": None}),
                     ("fgw_partial", {"coupling": "partial", "build_record": None})):
        c = harness.CANDIDATES[name]
        assert c.tasks == frozenset({"T2"}) and c.name == name
        assert c.build(name) == kw
        assert c.requires is harness.CANDIDATES["fgw"].requires
    assert harness.CANDIDATES["fgw_potts"].params["k"] == harness.NBHD_K == 15
    assert harness.CANDIDATES["fgw_potts"].params["label"] == "group"
    assert harness.CANDIDATES["fgw_partial"].params["label"] == "celltype"


def test_fgw_passes_celltype_to_the_partial_coupling_only(emb_expr, monkeypatch):
    seen = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        seen.append((kw.get("coupling"), kw.get("coupling_labels")))
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    harness.CANDIDATES["fgw"].build(c)
    harness.CANDIDATES["fgw_potts"].build(c)
    harness.CANDIDATES["fgw_partial"].build(c)
    assert seen[0] == (None, None) and seen[1] == (None, None)  # fgw and fgw_potts: the entry's call
    assert seen[2][0] == "partial"
    assert np.array_equal(np.asarray(seen[2][1][0]), np.asarray(s["before"].obs["celltype"]))
    assert np.array_equal(np.asarray(seen[2][1][1]), np.asarray(s["after"].obs["celltype"]))


# ================================================================================ the runner
@pytest.fixture(scope="module")
def n45():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n45_embryo_f_arms", HARNESS / "n45_embryo_f_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_fgw_and_the_two_arms_on_the_nine_tasks(n45):
    assert n45.NEW == ARMS and n45.ARMS == (harness.FLOOR, "fgw", *ARMS)
    tasks = n45.nt.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in ("fgw", *n45.NEW):
        assert all(harness.CANDIDATES[name].applies(t) for t in tasks.values())


def _patch_builds(monkeypatch, n45, fail_on: set[str], key: str, seed_fails: int = 0):
    """Every new arm's build replaced (it writes a record): the arms in fail_on raise the builder's refusal from
    seed_fails."""
    cands = harness.CANDIDATES

    def ok(ctx, build_record=None):
        build_record["seen"] = 1
        return "built"

    def bad(ctx, build_record=None):
        if ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
            raise ValueError(f"{key}: the partial FGW solver failed: ValueError: Error in the EMD resolution 3")
        build_record["seen"] = 1
        return "built"

    for x in n45.NEW:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=bad if x in fail_on else ok))
    monkeypatch.setitem(cands, "fgw", dataclasses.replace(cands["fgw"], build=lambda ctx: "built"))


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


def test_the_runners_task_file_loads_through_nine_tasks_with_fgw_under_its_own_name(tmp_path, monkeypatch, n45):
    score = {"fgw": 51.0, "fgw_potts": 52.0, "fgw_partial": 53.0}
    t = n45.nt.all_tasks()["e75"]
    _patch_builds(monkeypatch, n45, set(), t.key)
    runner, calls = _fake_runner(score)
    monkeypatch.setattr(n45.harness, "run_task", runner)
    out = tmp_path / "runs"
    n45.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert calls == [list(n45.ARMS)]  # one call with the floor, fgw and the two arms
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n45_embryo_f_arms.py" and data["not_buildable"] == {}
    assert data["build_records"] == {a: {str(s): {"seen": 1} for s in n45.SEEDS} for a in ARMS}
    scores, why, meta = n45.nt._load([out])
    assert why == [] and meta == ["runs:aaaaaaa"]
    assert set(scores[t.key]) == set(n45.ARMS) and "fgw" in scores[t.key]
    assert all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n45.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n45.main(["run", "--task", "nope", "--out", str(out)])


def test_an_arm_the_builder_refuses_is_recorded_and_the_others_are_still_scored(tmp_path, monkeypatch, n45):
    t = n45.nt.all_tasks()["e725"]
    score = {"fgw": 51.0, "fgw_potts": 53.0, "fgw_partial": 0.0}
    _patch_builds(monkeypatch, n45, {"fgw_partial"}, t.key, seed_fails=2)
    runner, calls = _fake_runner(score)
    res = n45.score_task(t, tmp_path, runner=runner)
    assert calls == [list(n45.ARMS), ["floor", "fgw"], ["floor", "fgw_potts"], ["floor", "fgw_partial"]]
    assert set(res["not_buildable"]) == {"fgw_partial"}
    rec = res["not_buildable"]["fgw_partial"]
    assert rec["status"] == n45.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == ("ValueError: <task>: the partial FGW solver failed: ValueError: Error in the EMD "
                              "resolution #")
    scored = {"floor", "fgw", "fgw_potts"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert set(res["build_records"]) == {"fgw_potts"}
    # every arm failing still leaves the floor and fgw scored
    _patch_builds(monkeypatch, n45, set(n45.NEW), t.key)
    runner, calls = _fake_runner(score)
    res = n45.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n45.NEW)
    assert all(set(r["scores"]) == {"floor", "fgw"} for r in res["rows"])


def test_a_failure_outside_a_new_arms_build_stops_the_task(tmp_path, monkeypatch, n45):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n45.score_task(n45.nt.all_tasks()["x2"], tmp_path, runner=runner)
    t = n45.nt.all_tasks()["x1"]
    _patch_builds(monkeypatch, n45, set(), t.key)

    def bad_fgw(ctx):
        raise ValueError("fgw fails")

    monkeypatch.setitem(harness.CANDIDATES, "fgw", dataclasses.replace(harness.CANDIDATES["fgw"], build=bad_fgw))
    runner, _ = _fake_runner({"fgw": 1.0, "fgw_potts": 1.0, "fgw_partial": 1.0})
    with pytest.raises(ValueError, match="fgw fails"):
        n45.score_task(t, tmp_path, runner=runner)
