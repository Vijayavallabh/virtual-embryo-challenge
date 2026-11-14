"""The embryo batch G (campaign branch campaign/emb-g): three controls for the partial coupling of campaign/emb-f, each
the entry fgw with one change to the displacement each anchor carries into the k-anchor mean
(vec.t2_build.interpolate_clouds(..., anchor_moves=...)): "damped" (fgw_damp: every anchor's displacement times s),
"random_still" (fgw_part_rand: a uniformly random subset of round(s x n) anchors per bracket keeps it, the others 0)
and "still" (fgw_still: every displacement 0), s fgw_partial's share; and the runner
experiments/harness/n46_embryo_g_arms.py. Synthetic clouds and stages only; no data file, committed result or scorer
is read.

Tolerances stated here: the damped arm's displacements against s times the entry's, 0 (bit for bit: the arm computes
exactly that product); its record's ratio of mean norms against s, 1e-12 relative (a ratio of two float sums); the
random arm's ratio against s, 0.1 absolute on 200 anchors (a random subset's share of the summed norms; the exact value
is checked separately)."""

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
BASE = "fcc1f6b"  # campaign/emb-f, the branch's base
ARMS = ("fgw_damp", "fgw_part_rand", "fgw_still")
MODES = {"fgw_damp": "damped", "fgw_part_rand": "random_still", "fgw_still": "still"}
RUN_ARMS = (*ARMS, "fgw_damp_m")  # fix brief M2: the runner's four controls
T = 0.4
N = 300


def _helpers(name: str):
    spec = importlib.util.spec_from_file_location(f"{name}_helpers_g", REPO / "tests" / f"{name}.py")
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


def _labels(A: np.ndarray, B: np.ndarray, flip: int = 120) -> tuple[np.ndarray, np.ndarray]:
    la = np.where(A[:, 0] > np.median(A[:, 0]), "p", "q")
    lb = np.where(B[:, 0] > np.median(B[:, 0]), "p", "q")
    lb[:flip] = "r"
    return la, lb


def _rng(seed: int = 1) -> np.random.Generator:
    return np.random.default_rng([seed, harness.ANCHOR_STILL_STREAM])


def _kw(mode: str, la, lb, rng_seed: int = 1) -> dict:
    kw: dict = {"anchor_moves": mode}
    if mode in ("damped", "random_still"):
        kw["anchor_labels"] = (la, lb)
    if mode == "random_still":
        kw["anchor_rng"] = _rng(rng_seed)
    return kw


def _entry_displacements(res) -> tuple[np.ndarray, np.ndarray]:
    """The entry's anchor displacements, recomputed from its coupling internals (its balanced barycentric partners)."""
    c = res.coupling
    ib, ia, Tm, Zb, Za_b = c["anchors_before"], c["anchors_after"], c["moves"], c["before_frame"], c["after_frame"]
    return tb.barycentric(Tm, Za_b[ia]) - Zb[ib], tb.barycentric(Tm.T, Zb[ib]) - Za_b[ia]


def _spy(monkeypatch) -> list:
    seen = []
    real = tb.controlled_displacements

    def spy(db, da, mode, s, rng):
        out = real(db, da, mode, s, rng)
        seen.append((np.array(db), np.array(da), s, out))
        return out

    monkeypatch.setattr(tb, "controlled_displacements", spy)
    return seen


def _strip(res):
    return dataclasses.replace(res, record={k: v for k, v in res.record.items()
                                            if k not in ("anchor_moves", "anchor_move")})


# ================================================================================ unset and the existing arms
def test_anchor_moves_unset_is_the_base_bit_for_bit_and_fgw_partial_is_unchanged(tmp_path, polar):
    A, B, E, EB = polar._pair()
    base = polar._base_module(tmp_path, BASE)
    old = base.interpolate_clouds(A, B, E, EB, T, N, **polar.KW)
    new = tb.interpolate_clouds(A, B, E, EB, T, N, **polar.KW)
    coupled = tb.interpolate_clouds(A, B, E, EB, T, N, anchor_moves="coupled", **polar.KW)
    polar._same(old, new)
    polar._same(new, coupled)
    la, lb = _labels(A, B)
    kw = dict(coupling="partial", coupling_labels=(la, lb))
    polar._same(base.interpolate_clouds(A, B, E, EB, T, N, **kw, **polar.KW),
                tb.interpolate_clouds(A, B, E, EB, T, N, **kw, **polar.KW))


def test_fgw_and_fgw_partial_call_the_builder_without_the_new_keywords(emb_expr, monkeypatch):
    seen = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        seen.append(kw)
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    c = emb_expr.ctx(inputs=emb_expr.stages(), roles=("before", "after"))
    harness.CANDIDATES["fgw"].build(c)
    harness.CANDIDATES["fgw_partial"].build(c)
    for kw in seen:
        assert not {"anchor_moves", "anchor_labels", "anchor_rng"} & set(kw)


# ================================================================================ the damped arm
def test_damped_every_displacement_is_s_times_the_entrys(polar, monkeypatch):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, with_coupling=True, **polar.KW)
    db0, da0 = _entry_displacements(entry)
    seen = _spy(monkeypatch)
    res = tb.interpolate_clouds(A, B, E, EB, T, N, with_coupling=True, **_kw("damped", la, lb), **polar.KW)
    (db_in, da_in, s, (db, da)), = seen
    assert np.array_equal(db_in, db0) and np.array_equal(da_in, da0)  # the entry's displacements go in
    assert 0 < s < 1
    assert np.array_equal(db, db0 * float(s)) and np.array_equal(da, da0 * float(s))
    # the entry's coupling, rotation and anchors; the cells move by the k-anchor mean of these displacements
    for key in ("anchors_before", "anchors_after", "fgw", "moves", "rotation", "before_frame", "after_frame"):
        assert np.array_equal(res.coupling[key], entry.coupling[key]), key
    Zb, Za_b = res.coupling["before_frame"], res.coupling["after_frame"]
    near_b = tb.nearest_anchors(Zb, Zb[res.coupling["anchors_before"]], tb.INTERP_K)[1]
    near_a = tb.nearest_anchors(Za_b, Za_b[res.coupling["anchors_after"]], tb.INTERP_K)[1]
    assert np.array_equal(res.coupling["shift_before"], db[near_b].mean(1))
    assert np.array_equal(res.coupling["shift_after"], da[near_a].mean(1))
    # s is fgw_partial's s on the same anchors
    part = tb.interpolate_clouds(A, B, E, EB, T, N, coupling="partial", coupling_labels=(la, lb), **polar.KW)
    ib, ia = res.coupling["anchors_before"], res.coupling["anchors_after"]
    assert s == tb.label_overlap(la[ib], lb[ia]) and part.record["partial"]["s"] == float(s)
    rec = res.record["anchor_move"]
    assert res.record["anchor_moves"] == "damped" and rec["s"] == float(s)
    for r in rec["shift_ratio"]:
        assert r == pytest.approx(float(s), rel=1e-12)
    assert rec["moving_share"] == [float((db0 != 0).any(1).mean()), float((da0 != 0).any(1).mean())]


def test_an_overlap_of_exactly_one_is_the_entry_bit_for_bit(polar):
    A, B, E, EB = polar._pair()
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, **polar.KW)
    one = (np.full(len(A), "x"), np.full(len(B), "x"))
    for mode in ("damped", "random_still"):
        res = tb.interpolate_clouds(A, B, E, EB, T, N, **_kw(mode, *one), **polar.KW)
        polar._same(_strip(res), entry)
        assert res.record["anchor_move"]["s"] == 1.0
        assert res.record["anchor_move"]["shift_ratio"] == [1.0, 1.0]


# ================================================================================ the random arm
def test_random_exactly_round_s_n_anchors_keep_the_entrys_displacement_and_the_rest_are_zero(polar, monkeypatch):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, with_coupling=True, **polar.KW)
    db0, da0 = _entry_displacements(entry)
    seen = _spy(monkeypatch)
    res = tb.interpolate_clouds(A, B, E, EB, T, N, **_kw("random_still", la, lb), **polar.KW)
    (db_in, da_in, s, (db, da)), = seen
    assert np.array_equal(db_in, db0) and np.array_equal(da_in, da0)
    assert 0 < s < 1
    ref = _rng()  # the subsets: the before bracket's draw first, then the after bracket's, on the arm's own stream
    for d, d0 in ((db, db0), (da, da0)):
        n = len(d0)
        m = round(s * n)
        assert m == int(np.rint(float(s) * n)) and 0 < m < n
        idx = ref.choice(n, m, replace=False)
        keep = np.zeros(n, dtype=bool)
        keep[idx] = True
        assert np.array_equal(d[keep], d0[keep])  # bit for bit
        assert np.array_equal(d[~keep], np.zeros_like(d0[~keep]))
        assert (d0[keep] != 0).any(1).all()  # the entry's kept displacements are not 0 here
    rec = res.record["anchor_move"]
    assert rec["s"] == float(s)
    assert rec["moving_share"] == [round(s * len(db0)) / len(db0), round(s * len(da0)) / len(da0)]
    for r, d, d0 in zip(rec["shift_ratio"], (db, da), (db0, da0), strict=True):
        exact = float(np.linalg.norm(d, axis=1).mean()) / float(np.linalg.norm(d0, axis=1).mean())
        assert r == exact and abs(r - float(s)) < 0.1


def test_the_random_stream_is_new_and_another_seed_gives_another_subset(polar, monkeypatch):
    from vec import t2_latent

    used = {*range(16), *range(20, 28), harness.NBHD_STREAM, harness.SERIES_XT_STREAM, *harness.POTTS_STREAMS,
            *tb.HERD_STREAMS, *tb.MIDPOINT_STREAMS, t2_latent.EXPRESSION_STREAM}
    a = harness.ANCHOR_STILL_STREAM
    assert a not in used
    seed = harness.BUILD_SEED_OFFSET + 4
    draws = {s: np.random.default_rng([seed, s]).random(64) for s in sorted(used | {a})}
    for s in sorted(used):
        assert not np.array_equal(draws[a], draws[s])
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    seen = _spy(monkeypatch)
    x = tb.interpolate_clouds(A, B, E, EB, T, N, **_kw("random_still", la, lb, rng_seed=1), **polar.KW)
    y = tb.interpolate_clouds(A, B, E, EB, T, N, **_kw("random_still", la, lb, rng_seed=2), **polar.KW)
    (_, _, _, (db1, da1)), (_, _, _, (db2, da2)) = seen
    assert not np.array_equal((db1 != 0).any(1), (db2 != 0).any(1))
    assert not np.array_equal((da1 != 0).any(1), (da2 != 0).any(1))
    assert np.array_equal(x.rows, y.rows) and not np.array_equal(x.coords, y.coords)


def test_the_arm_draws_its_subset_on_its_own_stream(emb_expr, monkeypatch):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    got = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        if "anchor_rng" in kw:
            got.append(kw["anchor_rng"].bit_generator.state)
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    harness.CANDIDATES["fgw_part_rand"].build(c)
    assert got == [np.random.default_rng([c.seed, harness.ANCHOR_STILL_STREAM]).bit_generator.state]


# ================================================================================ the still arm
def test_still_every_displacement_is_zero_and_no_label_is_read(polar, monkeypatch):
    A, B, E, EB = polar._pair()
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, with_coupling=True, **polar.KW)
    seen = _spy(monkeypatch)
    res = tb.interpolate_clouds(A, B, E, EB, T, N, anchor_moves="still", with_coupling=True, **polar.KW)
    (_, _, s, (db, da)), = seen
    assert s is None and not db.any() and not da.any()
    assert res.record["anchor_move"] == {"moving_share": [0.0, 0.0], "shift_ratio": [0.0, 0.0]}
    # no anchor moves: each bracket's cells at their own normalised positions, the after bracket rotated by the
    # entry's rotation into the before frame
    assert np.array_equal(res.coupling["rotation"], entry.coupling["rotation"])
    Zb, Za_b = res.coupling["before_frame"], res.coupling["after_frame"]
    Za = (B - B.mean(0)) / tb.rms_radius(B)
    assert np.array_equal(Za_b, Za @ entry.coupling["rotation"])
    pool = np.vstack([Zb, Za_b])
    Z = np.empty((N, 3))
    Z[res.stage == 0], Z[res.stage == 1] = Zb[res.rows[res.stage == 0]], Za_b[res.rows[res.stage == 1]]
    np.testing.assert_allclose(res.pool_coords, A.mean(0) + (pool - Z.mean(0)) * res.record["scale_from_canonical"],
                               rtol=0, atol=1e-12)
    with pytest.raises(ValueError, match="anchor_labels"):
        tb.interpolate_clouds(A, B, E, EB, T, N, anchor_moves="still", anchor_labels=_labels(A, B), **polar.KW)


def test_the_still_arm_does_not_depend_on_the_labels(emb_expr):
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    x = harness.CANDIDATES["fgw_still"].build(c)
    s2 = emb_expr.stages()
    for a in s2.values():
        a.obs["celltype"] = np.random.default_rng(0).permutation(a.obs["celltype"].astype(str).to_numpy())
    y = harness.CANDIDATES["fgw_still"].build(emb_expr.ctx(inputs=s2, roles=("before", "after")))
    assert emb_expr.digest(x) == emb_expr.digest(y)


# ================================================================================ all three
def _partly_shared(emb_expr):
    s = emb_expr.stages()
    ct = s["after"].obs["celltype"].astype(str).to_numpy().copy()
    ct[:50] = "PS"
    s["after"].obs["celltype"] = ct
    return emb_expr.ctx(inputs=s, roles=("before", "after"))


@pytest.mark.parametrize("arm", ARMS)
def test_each_arm_keeps_fgws_cells_and_rows_and_moves_only_positions(emb_expr, arm):
    c = _partly_shared(emb_expr)
    fgw = harness.CANDIDATES["fgw"].build(c)
    rec: dict = {}
    x = getattr(harness, f"_{arm}")(c, build_record=rec)
    y = harness.CANDIDATES[arm].build(c)
    assert emb_expr.digest(x) == emb_expr.digest(y)  # same seed, same output
    assert np.array_equal(x.X, fgw.X) and x.n == fgw.n  # the drawn cells' expression rows, bit for bit
    assert not np.array_equal(x.coords, fgw.coords)  # the idea is not idle here
    r = rec["anchor_move"]
    assert r["anchor_moves"] == MODES[arm]
    if arm == "fgw_still":
        assert "s" not in r and r["shift_ratio"] == [0.0, 0.0] and r["moving_share"] == [0.0, 0.0]
    else:
        assert 0.0 < r["s"] < 1.0
        part: dict = {}
        harness._fgw_partial(c, build_record=part)
        assert r["s"] == part["partial"]["s"]
        tol = 1e-12 if arm == "fgw_damp" else 0.1
        assert all(abs(q - r["s"]) <= tol * (r["s"] if arm == "fgw_damp" else 1.0) for q in r["shift_ratio"])
    assert set(r) <= {"anchor_moves", "s", "moving_share", "shift_ratio"}


@pytest.mark.parametrize("mode", ["damped", "random_still", "still"])
def test_each_mode_keeps_the_entrys_draws_and_is_reproducible(polar, mode):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, **polar.KW)
    a = tb.interpolate_clouds(A, B, E, EB, T, N, **_kw(mode, la, lb), **polar.KW)
    b = tb.interpolate_clouds(A, B, E, EB, T, N, **_kw(mode, la, lb), **polar.KW)
    polar._same(a, b)
    assert np.array_equal(a.rows, entry.rows) and np.array_equal(a.stage, entry.stage)
    assert np.array_equal(a.pool_weights, entry.pool_weights)
    assert not np.array_equal(a.coords, entry.coords)
    assert a.record["n_jittered"] >= 0 and a.record["rotation_deg"] == entry.record["rotation_deg"]


def test_the_arm_with_equal_label_shares_is_fgw_bit_for_bit(emb_expr):
    c = emb_expr.ctx(inputs=emb_expr.stages(), roles=("before", "after"))
    fgw = emb_expr.digest(harness.CANDIDATES["fgw"].build(c))
    for arm in ("fgw_damp", "fgw_part_rand"):
        rec: dict = {}
        x = getattr(harness, f"_{arm}")(c, build_record=rec)
        assert rec["anchor_move"]["s"] == 1.0  # the synthetic brackets hold each label in the same share
        assert emb_expr.digest(x) == fgw


def test_the_anchor_moves_arguments_are_checked(polar):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    with pytest.raises(ValueError, match="anchor_moves must be"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_moves="half", **polar.KW)
    with pytest.raises(ValueError, match="anchor_labels"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_moves="damped", **polar.KW)
    with pytest.raises(ValueError, match="anchor_labels"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_labels=(la, lb), **polar.KW)
    with pytest.raises(ValueError, match="anchor_rng"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_moves="random_still", anchor_labels=(la, lb), **polar.KW)
    with pytest.raises(ValueError, match="anchor_rng"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_moves="damped", anchor_labels=(la, lb), anchor_rng=_rng(),
                              **polar.KW)
    with pytest.raises(ValueError, match="one label per cell"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_moves="damped", anchor_labels=(la[:5], lb), **polar.KW)
    for kw in ({"coupling": "partial", "coupling_labels": (la, lb)}, {"move_coupling": "aligned"},
               {"partners": "mirror_fixed"}, {"moves": "polar"}):
        with pytest.raises(ValueError, match="anchor displacements only"):
            tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_moves="still", **kw, **polar.KW)


def test_the_still_anchor_subset_rounds_exactly():
    rng = np.random.default_rng(0)
    assert tb.still_anchor_subset(10, Fraction(1, 4), rng).sum() == 2  # 2.5 to the even 2, as rint
    assert tb.still_anchor_subset(10, Fraction(7, 20), rng).sum() == 4  # 3.5 to 4
    assert tb.still_anchor_subset(10, Fraction(1), rng).all()
    assert not tb.still_anchor_subset(10, Fraction(0), rng).any()


# ================================================================================ the harness arms
def test_the_three_arms_are_registered_for_t2_and_call_fgw_with_their_keywords(monkeypatch):
    seen = {}
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: seen.setdefault(ctx, kw))
    for name in ARMS:
        c = harness.CANDIDATES[name]
        assert c.tasks == frozenset({"T2"}) and c.name == name
        assert c.build(name) == {"anchor_moves": MODES[name], "build_record": None}
        assert c.requires is harness.CANDIDATES["fgw"].requires
    assert harness.CANDIDATES["fgw_damp"].params["label"] == "celltype"
    assert harness.CANDIDATES["fgw_part_rand"].params["stream"] == harness.ANCHOR_STILL_STREAM
    assert "label" not in harness.CANDIDATES["fgw_still"].params


def test_fgw_passes_celltype_to_the_damped_and_random_arms_only(emb_expr, monkeypatch):
    seen = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        seen.append(kw)
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    s = emb_expr.stages()
    c = emb_expr.ctx(inputs=s, roles=("before", "after"))
    for name in ARMS:
        harness.CANDIDATES[name].build(c)
    for kw, name in zip(seen, ARMS, strict=True):
        assert kw["anchor_moves"] == MODES[name]
        if name == "fgw_still":
            assert "anchor_labels" not in kw and "anchor_rng" not in kw
        else:
            assert np.array_equal(np.asarray(kw["anchor_labels"][0]), np.asarray(s["before"].obs["celltype"]))
            assert np.array_equal(np.asarray(kw["anchor_labels"][1]), np.asarray(s["after"].obs["celltype"]))
        assert ("anchor_rng" in kw) == (name == "fgw_part_rand")
        assert "coupling" not in kw


# ================================================================================ the runner
@pytest.fixture(scope="module")
def n46():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n46_embryo_g_arms", HARNESS / "n46_embryo_g_arms.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_scores_the_floor_fgw_fgw_partial_and_the_three_arms_on_the_nine_tasks(n46):
    # fix brief M2: the fourth control fgw_damp_m joins the three
    assert n46.NEW == RUN_ARMS and n46.ARMS == (harness.FLOOR, "fgw", "fgw_partial", *RUN_ARMS)
    tasks = n46.nt.all_tasks()
    assert set(tasks) == {"half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"}
    for name in ("fgw", "fgw_partial", *n46.NEW):
        assert all(harness.CANDIDATES[name].applies(t) for t in tasks.values())


FAKE_SHIFT = {"fgw": {"anchor": [2.0, 4.0], "cell": [1.0, 0.0]}, "arm": {"anchor": [1.0, 1.0], "cell": [0.5, 0.0]}}


def _patch_builds(monkeypatch, n46, fail_on: set[str], key: str, seed_fails: int = 0):
    """Every guarded arm's build replaced (each new arm writes a record; fix brief M2: fgw, fgw_partial and every
    new arm also write their shifts, fgw_damp_m its factors): the arms in fail_on raise the builder's refusal from
    seed_fails."""
    cands = harness.CANDIDATES

    def make(name):
        def build(ctx, build_record=None):
            if name in fail_on and ctx.seed - harness.BUILD_SEED_OFFSET >= seed_fails:
                raise ValueError(f"{key}: the partial FGW solver failed: ValueError: Error in the EMD resolution 3")
            if name in n46.NEW:
                build_record["seen"] = 1
            if name == "fgw_damp_m":
                build_record["anchor_move"] = {"factors": [0.5, 0.25]}
            build_record["shift"] = FAKE_SHIFT["arm"]
            return "built"

        return build

    def fgw(ctx, build_record=None):
        build_record["shift"] = FAKE_SHIFT["fgw"]
        return "built"

    for x in n46.GUARDED:
        monkeypatch.setitem(cands, x, dataclasses.replace(cands[x], build=make(x)))
    monkeypatch.setitem(cands, "fgw", dataclasses.replace(cands["fgw"], build=fgw))


def _fake_runner(score: dict[str, float], ctx_of=None):
    """run_task's shape without data or scorer: each candidate is built (its build may raise) and given a constant.
    ctx_of(seed) gives the build's context (default: a context holding the seed only)."""
    calls = []
    if ctx_of is None:
        def ctx_of(seed):
            return dataclasses.make_dataclass("Ctx", ["seed"])(seed=seed)

    def runner(t, cands, seeds, data_dir):
        calls.append([c.name for c in cands])
        rows = []
        for s in seeds:
            row = {"seed": s, "scores": {"floor": 50.0}, "skills": {"floor": {}}, "metrics": {"floor": {}},
                   "n_cells": {"floor": 1}, "groups": {"floor": {}}, "reads": {}}
            for c in cands:
                if c.name == "floor":
                    continue
                c.build(ctx_of(harness.BUILD_SEED_OFFSET + s))
                row["scores"][c.name] = score[c.name]
                for k in ("skills", "metrics", "groups"):
                    row[k][c.name] = {}
                row["n_cells"][c.name], row["reads"][c.name] = 1, []
            rows.append(row)
        return {"task": {"key": t.key, "inputs": dict(t.inputs)}, "rows": rows,
                "candidates": [c.name for c in cands], "params": {c.name: {} for c in cands},
                "provenance": {"commit": "a" * 40, "tree_clean": True}}

    return runner, calls


SCORE = {"fgw": 51.0, "fgw_partial": 52.0, "fgw_damp": 53.0, "fgw_part_rand": 54.0, "fgw_still": 55.0,
         "fgw_damp_m": 56.0}
FAKE_RATIO = {"anchor": [0.5, 0.25], "cell": [0.5, None]}


def test_the_runners_task_file_loads_through_nine_tasks_with_fgw_and_fgw_partial_under_their_names(tmp_path,
                                                                                                  monkeypatch, n46):
    t = n46.nt.all_tasks()["e75"]
    _patch_builds(monkeypatch, n46, set(), t.key)
    runner, calls = _fake_runner(SCORE)
    monkeypatch.setattr(n46.harness, "run_task", runner)
    out = tmp_path / "runs"
    n46.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert calls == [list(n46.ARMS)]  # one call with the floor, fgw, fgw_partial and the three arms
    (path,) = out.glob("T2e_*.json")
    data = json.loads(path.read_text())
    assert data["task_source"] == "experiments/harness/n46_embryo_g_arms.py" and data["not_buildable"] == {}
    # fix brief M2: fgw_partial and the four controls carry movement_ratio (fgw's record and raw shifts not written)
    want = {a: {str(s): {"seen": 1, "movement_ratio": FAKE_RATIO} for s in n46.SEEDS} for a in RUN_ARMS}
    want["fgw_partial"] = {str(s): {"movement_ratio": FAKE_RATIO} for s in n46.SEEDS}
    for s in n46.SEEDS:
        want["fgw_damp_m"][str(s)].update(anchor_move={"factors": [0.5, 0.25]}, factors=[0.5, 0.25])
    assert data["build_records"] == want
    scores, why, meta = n46.nt._load([out])
    assert why == [] and meta == ["runs:aaaaaaa"]
    assert set(scores[t.key]) == set(n46.ARMS) and {"fgw", "fgw_partial"} <= set(scores[t.key])
    assert all(len(v) == 10 for v in scores[t.key].values())
    with pytest.raises(SystemExit, match="exists"):
        n46.main(["run", "--task", "e75", "--out", str(out), "--data", str(tmp_path / "no_data")])
    assert len(calls) == 1
    with pytest.raises(SystemExit, match="unknown task"):
        n46.main(["run", "--task", "nope", "--out", str(out)])


def test_an_arm_the_builder_refuses_is_recorded_and_the_others_are_still_scored(tmp_path, monkeypatch, n46):
    t = n46.nt.all_tasks()["e725"]
    _patch_builds(monkeypatch, n46, {"fgw_partial", "fgw_damp"}, t.key, seed_fails=2)
    runner, calls = _fake_runner(SCORE)
    res = n46.score_task(t, tmp_path, runner=runner)
    assert calls == [list(n46.ARMS), ["floor", "fgw"], ["floor", "fgw_partial"], ["floor", "fgw_damp"],
                     ["floor", "fgw_part_rand"], ["floor", "fgw_still"], ["floor", "fgw_damp_m"]]
    assert set(res["not_buildable"]) == {"fgw_partial", "fgw_damp"}
    rec = res["not_buildable"]["fgw_damp"]
    assert rec["status"] == n46.n36.NOT_BUILDABLE and rec["seed"] == 2
    assert rec["message"] == ("ValueError: <task>: the partial FGW solver failed: ValueError: Error in the EMD "
                              "resolution #")
    scored = {"floor", "fgw", "fgw_part_rand", "fgw_still", "fgw_damp_m"}
    assert all(set(r["scores"]) == scored for r in res["rows"]) and set(res["candidates"]) == scored
    assert set(res["build_records"]) == {"fgw_part_rand", "fgw_still", "fgw_damp_m"}
    # every guarded arm failing still leaves the floor and fgw scored
    _patch_builds(monkeypatch, n46, set(n46.GUARDED), t.key)
    runner, calls = _fake_runner(SCORE)
    res = n46.score_task(t, tmp_path, runner=runner)
    assert set(res["not_buildable"]) == set(n46.GUARDED)
    assert all(set(r["scores"]) == {"floor", "fgw"} for r in res["rows"])


def test_a_failure_outside_a_guarded_arms_build_stops_the_task(tmp_path, monkeypatch, n46):
    def runner(t, cands, seeds, data_dir):
        raise OSError("the data are not there")

    with pytest.raises(OSError, match="not there"):
        n46.score_task(n46.nt.all_tasks()["x2"], tmp_path, runner=runner)
    t = n46.nt.all_tasks()["x1"]
    _patch_builds(monkeypatch, n46, set(), t.key)

    def bad_fgw(ctx, build_record=None):
        raise ValueError("fgw fails")

    monkeypatch.setitem(harness.CANDIDATES, "fgw", dataclasses.replace(harness.CANDIDATES["fgw"], build=bad_fgw))
    runner, _ = _fake_runner(SCORE)
    with pytest.raises(ValueError, match="fgw fails"):
        n46.score_task(t, tmp_path, runner=runner)


# ================================================================================ fix brief M2: fgw_damp_m
# Tolerances: fgw_damp_m's mean cell shift against fgw_partial's, 1e-9 relative (the same partial displacements; the
# after bracket's nearest-anchor indices are the entry's against fgw_partial's own, the same up to the rounding of a
# rotation); a cell-level ratio against the factor it was scaled by, 1e-12 relative.
DIGESTS_4E3C6EC = {  # the five arms on _partly_shared, computed at 4e3c6ec before fix M2 (tests/test_emb_expr digest)
    "fgw": "4a53a67d0be4e5195b3ddf8e554b5ac2ac0b506bfb8808f747ab48cea1d8164f",
    "fgw_partial": "3c6e99b0548240ff9e578ee1c6380417ebcbe0e0b01385616e1a952989afb77a",
    "fgw_damp": "8e29579079a375747cdaced71865d7b83f038ab24048c355b02cd77d30dda7a1",
    "fgw_part_rand": "7743a92233c1b3121a80f8726ebd579db62fa748827bb14558fd4dfb62dd5eb0",
    "fgw_still": "9754173a91bc162a9d59e10eb87b282f75b8877731ca82c63185b139babdcac4",
}


def _spy_m(monkeypatch) -> list:
    seen = []
    real = tb.controlled_displacements

    def spy(db, da, mode, s, rng, factors=None):
        out = real(db, da, mode, s, rng, factors=factors)
        seen.append((np.array(db), np.array(da), factors, out))
        return out

    monkeypatch.setattr(tb, "controlled_displacements", spy)
    return seen


def _m_kw(la, lb) -> dict:
    return {"anchor_moves": "damped_matched", "anchor_labels": (la, lb)}


def test_the_five_arms_keep_their_4e3c6ec_files(emb_expr):
    c = _partly_shared(emb_expr)
    for name, d in DIGESTS_4E3C6EC.items():
        assert emb_expr.digest(harness.CANDIDATES[name].build(c)) == d, name


def test_damp_m_scales_each_brackets_displacements_by_its_factor_and_matches_fgw_partials_cell_shift(polar,
                                                                                                     monkeypatch):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, with_coupling=True, **polar.KW)
    part = tb.interpolate_clouds(A, B, E, EB, T, N, coupling="partial", coupling_labels=(la, lb), **polar.KW)
    db0, da0 = _entry_displacements(entry)
    seen = _spy_m(monkeypatch)
    res = tb.interpolate_clouds(A, B, E, EB, T, N, with_coupling=True, **_m_kw(la, lb), **polar.KW)
    (db_in, da_in, factors, (db, da)), = seen
    assert np.array_equal(db_in, db0) and np.array_equal(da_in, da0)
    fb, fa = factors
    assert 0 < fb < 1 and 0 < fa < 1 and fb != fa
    assert np.array_equal(db, db0 * fb) and np.array_equal(da, da0 * fa)
    rec = res.record["anchor_move"]
    assert res.record["anchor_moves"] == "damped_matched" and rec["factors"] == [fb, fa]
    assert rec["s"] == part.record["partial"]["s"] and 0 < rec["s"] < 1
    assert rec["moving_share"] == [float((db0 != 0).any(1).mean()), float((da0 != 0).any(1).mean())]
    # each bracket's cells move on average as far as under fgw_partial; its anchors in general do not
    np.testing.assert_allclose(res.movement["mean_cell_shift"], part.movement["mean_cell_shift"], rtol=1e-9, atol=0)
    assert not np.allclose(res.movement["mean_anchor_shift"], part.movement["mean_anchor_shift"], rtol=1e-3)
    # the factors are the cell-level ratios of fgw_partial's shifts to the entry's
    for i, f in enumerate((fb, fa)):
        assert f == pytest.approx(part.movement["mean_cell_shift"][i] / entry.movement["mean_cell_shift"][i],
                                  rel=1e-9)
    # the entry's coupling, rotation, anchors and partners (directions); only the lengths change
    for key in ("anchors_before", "anchors_after", "fgw", "moves", "rotation", "before_frame", "after_frame"):
        assert np.array_equal(res.coupling[key], entry.coupling[key]), key
    assert res.record["rotation_deg"] == entry.record["rotation_deg"]
    np.testing.assert_allclose(res.coupling["shift_before"], fb * entry.coupling["shift_before"], rtol=1e-12, atol=0)
    np.testing.assert_allclose(res.coupling["shift_after"], fa * entry.coupling["shift_after"], rtol=1e-12, atol=0)


def test_damp_m_with_an_overlap_of_one_is_fgw_bit_for_bit_without_the_partial_solver(polar, monkeypatch):
    A, B, E, EB = polar._pair()
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, **polar.KW)

    def boom(*a, **kw):
        raise AssertionError("the partial solver is called")

    monkeypatch.setattr(tb, "partial_fgw_coupling", boom)
    one = (np.full(len(A), "x"), np.full(len(B), "x"))
    res = tb.interpolate_clouds(A, B, E, EB, T, N, **_m_kw(*one), **polar.KW)
    polar._same(_strip(res), entry)
    assert res.movement == entry.movement
    assert res.record["anchor_move"]["factors"] == [1.0, 1.0] and res.record["anchor_move"]["s"] == 1.0


def test_damp_m_refuses_on_a_solver_failure_and_a_zero_denominator(polar, emb_expr, monkeypatch):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)

    def fail(*a, **kw):
        raise ValueError("the partial FGW solver failed: ValueError: stub")

    with monkeypatch.context() as m:
        m.setattr(tb, "partial_fgw_coupling", fail)
        with pytest.raises(ValueError, match="partial FGW solver failed"):
            tb.interpolate_clouds(A, B, E, EB, T, N, **_m_kw(la, lb), **polar.KW)
        with pytest.raises(ValueError, match="partial FGW solver failed"):
            harness.CANDIDATES["fgw_damp_m"].build(_partly_shared(emb_expr))
    with pytest.raises(ValueError, match="mean cell shift of a bracket is 0"):
        tb.movement_factors([0.0, 1.0], [0.5, 0.5])
    assert tb.movement_factors([2.0, 4.0], [1.0, 1.0]) == (0.5, 0.25)
    with pytest.raises(ValueError, match="anchor_labels"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, anchor_moves="damped_matched", **polar.KW)
    with pytest.raises(ValueError, match="anchor displacements only"):
        tb.interpolate_clouds(A, B, E, EB, T, 30, coupling="partial", coupling_labels=(la, lb), **_m_kw(la, lb),
                              **polar.KW)


def test_damp_m_keeps_fgws_draws_and_rows_and_is_reproducible(polar, emb_expr):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    entry = tb.interpolate_clouds(A, B, E, EB, T, N, **polar.KW)
    a = tb.interpolate_clouds(A, B, E, EB, T, N, **_m_kw(la, lb), **polar.KW)
    b = tb.interpolate_clouds(A, B, E, EB, T, N, **_m_kw(la, lb), **polar.KW)
    polar._same(a, b)
    assert np.array_equal(a.rows, entry.rows) and np.array_equal(a.stage, entry.stage)
    assert np.array_equal(a.pool_weights, entry.pool_weights) and not np.array_equal(a.coords, entry.coords)
    c = _partly_shared(emb_expr)
    fgw = harness.CANDIDATES["fgw"].build(c)
    rec: dict = {}
    x = harness._fgw_damp_m(c, build_record=rec)
    y = harness.CANDIDATES["fgw_damp_m"].build(c)
    assert emb_expr.digest(x) == emb_expr.digest(y)
    assert np.array_equal(x.X, fgw.X) and x.n == fgw.n and not np.array_equal(x.coords, fgw.coords)
    assert rec["anchor_move"]["anchor_moves"] == "damped_matched" and len(rec["anchor_move"]["factors"]) == 2
    c1 = emb_expr.ctx(inputs=emb_expr.stages(), roles=("before", "after"))  # equal label shares: s = 1
    assert emb_expr.digest(harness.CANDIDATES["fgw_damp_m"].build(c1)) == emb_expr.digest(
        harness.CANDIDATES["fgw"].build(c1))


def test_damp_m_is_registered_like_fgw_partial_and_gets_the_celltype_labels(emb_expr, monkeypatch):
    c = harness.CANDIDATES["fgw_damp_m"]
    assert c.tasks == frozenset({"T2"}) and c.requires is harness.CANDIDATES["fgw_partial"].requires
    assert c.params["label"] == "celltype" and c.params["anchor_moves"] == "damped_matched"
    seen = []
    real = tb.interpolate_clouds

    def spy(*a, **kw):
        seen.append(kw)
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    s = emb_expr.stages()
    c.build(emb_expr.ctx(inputs=s, roles=("before", "after")))
    (kw,) = seen
    assert kw["anchor_moves"] == "damped_matched" and "anchor_rng" not in kw and "coupling" not in kw
    assert np.array_equal(np.asarray(kw["anchor_labels"][0]), np.asarray(s["before"].obs["celltype"]))
    assert np.array_equal(np.asarray(kw["anchor_labels"][1]), np.asarray(s["after"].obs["celltype"]))


# ================================================================================ fix brief M2: the movement record
@pytest.mark.parametrize("kw", ["entry", "partial", "damped", "random_still", "still", "damped_matched"])
def test_mean_cell_shift_is_present_for_every_arm_and_equals_a_by_hand_computation(polar, kw):
    A, B, E, EB = polar._pair()
    la, lb = _labels(A, B)
    fixed = {"entry": {}, "partial": {"coupling": "partial", "coupling_labels": (la, lb)},
             "damped_matched": _m_kw(la, lb)}
    extra = fixed[kw] if kw in fixed else _kw(kw, la, lb)
    res = tb.interpolate_clouds(A, B, E, EB, T, N, with_coupling=True, **extra, **polar.KW)
    by_hand = [float(np.mean([np.sqrt((v ** 2).sum()) for v in res.coupling[key]]))
               for key in ("shift_before", "shift_after")]
    np.testing.assert_allclose(res.movement["mean_cell_shift"], by_hand, rtol=1e-12, atol=0)
    assert res.movement["mean_anchor_shift"] == res.record["mean_anchor_shift"]
    assert "mean_cell_shift" not in res.record  # every record stays as it was


def test_mean_cell_shift_on_a_small_case_by_hand():
    d = np.array([[1.0, 0, 0], [0, 2.0, 0], [0, 0, 0], [3.0, 4.0, 0]])
    near = np.array([[0, 1], [2, 3], [0, 2]])
    # the cells' mean displacements (0.5, 1, 0), (1.5, 2, 0), (0.5, 0, 0): norms sqrt(1.25), 2.5, 0.5
    assert tb.mean_cell_shift(d, near) == pytest.approx((np.sqrt(1.25) + 2.5 + 0.5) / 3, rel=1e-15)


def test_the_runner_writes_movement_ratios_from_real_builds(tmp_path, emb_expr, n46):
    c = _partly_shared(emb_expr)

    def ctx_of(seed):
        return dataclasses.replace(c, seed=seed)

    runner, _ = _fake_runner(SCORE, ctx_of=ctx_of)
    res = n46.score_task(n46.nt.all_tasks()["e75"], tmp_path, runner=runner, seeds=(0, 1))
    br = res["build_records"]
    assert set(br) == {"fgw_partial", *RUN_ARMS} and res["not_buildable"] == {}
    for arm in br:
        for seed in ("0", "1"):
            r = br[arm][seed]
            assert "shift" not in r and set(r["movement_ratio"]) == {"anchor", "cell"}
            assert all(len(v) == 2 for v in r["movement_ratio"].values())
    for seed in ("0", "1"):
        s = br["fgw_damp"][seed]["anchor_move"]["s"]
        assert 0 < s < 1
        np.testing.assert_allclose(br["fgw_damp"][seed]["movement_ratio"]["cell"], [s, s], rtol=1e-12)
        np.testing.assert_allclose(br["fgw_damp"][seed]["movement_ratio"]["anchor"], [s, s], rtol=1e-12)
        assert br["fgw_still"][seed]["movement_ratio"] == {"anchor": [0.0, 0.0], "cell": [0.0, 0.0]}
        m = br["fgw_damp_m"][seed]
        np.testing.assert_allclose(m["movement_ratio"]["cell"], m["factors"], rtol=1e-12)
        np.testing.assert_allclose(m["movement_ratio"]["cell"], br["fgw_partial"][seed]["movement_ratio"]["cell"],
                                   rtol=1e-9)
        assert m["factors"] == m["anchor_move"]["factors"]
        assert "factors" not in br["fgw_damp"][seed]
    json.dumps(br)  # floats and None only
