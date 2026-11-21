"""Row E52 as specified (the captain's order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (v)): the copy-rows
file at the team's pre-registered size, vec.t2_build.copy_rows_sized and the harness arm fgw_copyrows_sized, on
synthetic arrays and stages only. No released stage, board file, size file or test count is read."""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from test_t2_build import _labelled_pair

from vec import evaluate as ev
from vec import harness
from vec import t2_build as tb
from vec.boards import Board
from vec.harness import CANDIDATES

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import e52_copy_rows_sized_read as e52  # noqa: E402

SEED = 1_000_003


def _bits(C):
    return np.ascontiguousarray(np.asarray(C, dtype=np.float32)).view(np.uint32)


def _cloud(n, seed, scales=(5.0, 2.0, 0.7), offset=(100.0, -50.0, 20.0), dup_frac=0.02):
    """An anisotropic, offset float32 cloud with about dup_frac of its rows exact copies of other rows, and the
    (copy, original) pairs."""
    rng = np.random.default_rng(seed)
    C = (rng.normal(size=(n, 3)) * np.asarray(scales) + np.asarray(offset)).astype(np.float32)
    k = max(1, int(round(dup_frac * n)))
    copies = rng.choice(np.arange(1, n), k, replace=False)
    sources = np.setdiff1d(np.arange(n), copies)  # an original is never itself a copy
    originals = rng.choice(sources, k)
    C[copies] = C[originals]
    return C, list(zip(copies.tolist(), originals.tolist(), strict=True))


# ------------------------------------------------------------------------------------------------ 1. arithmetic
@pytest.mark.parametrize("n", [2000, 9000, 40000])
@pytest.mark.parametrize("t", [0.37, 3.25, 41.0])
def test_the_arithmetic_is_the_specifications_bit_for_bit(n, t):
    floor32, pairs = _cloud(n, n)
    out = tb.copy_rows_sized(floor32, t)
    assert out.dtype == np.float64 and out.shape == floor32.shape
    S32 = out.astype(np.float32)
    C64 = floor32.astype(np.float64)
    mu = C64.mean(0)
    want = (mu + (C64 - mu) * (t / tb.rms_radius(C64))).astype(np.float32)
    np.testing.assert_array_equal(_bits(S32), _bits(want))  # shape, row order and float32 bits
    assert abs(tb.rms_radius(S32.astype(np.float64)) / t - 1) <= 1e-6
    for c, o in pairs:  # exact duplicates stay exact duplicates
        assert np.array_equal(S32[c].view(np.uint32), S32[o].view(np.uint32))
    assert len(np.unique(S32, axis=0)) == len(np.unique(floor32, axis=0))
    np.testing.assert_array_equal(_bits(tb.copy_rows_sized(floor32.astype(np.float64), t).astype(np.float32)),
                                  _bits(S32))  # float32 or float64 input: the same float64 arithmetic


@pytest.mark.parametrize("target", [0.0, -1.0, float("inf"), float("-inf"), float("nan")])
def test_a_target_that_is_not_finite_and_positive_raises(target):
    with pytest.raises(ValueError, match="finite and positive"):
        tb.copy_rows_sized(_cloud(2000, 1)[0], target)


def test_a_degenerate_cloud_raises():
    with pytest.raises(ValueError, match="no spread"):
        tb.copy_rows_sized(np.full((2000, 3), 4.5, np.float32), 3.0)


def test_the_rule_text_and_the_builder_rule_name_the_step():
    assert "rescale_to_rms" in tb.COPY_ROWS_SIZE_RULE and "centroid" in tb.COPY_ROWS_SIZE_RULE
    assert "no jitter" in tb.COPY_ROWS_SIZE_RULE


# ------------------------------------------------------------------------------------------------ 2. the scorer
@pytest.fixture(scope="module")
def scored():
    floor32, _ = _cloud(2000, 7)
    truth = (np.random.default_rng(8).normal(size=(3000, 3)) * np.array([4.0, 2.5, 1.0]) + 7.0).astype(np.float32)
    return floor32, truth


@pytest.mark.parametrize("build", ["main", "prefix"])
@pytest.mark.parametrize("s", [0, 1])
@pytest.mark.parametrize("t", [0.8, 3.25])
def test_the_scorer_terms_are_the_floors_except_the_size(scored, build, s, t):
    floor32, truth = scored
    S32 = tb.copy_rows_sized(floor32, t).astype(np.float32)
    sm = ev._shape_metrics()
    with ev.scorer_build(build):
        d2a = sm.d2_distance(S32, truth, seed=s)
        d2f = sm.d2_distance(floor32, truth, seed=s)
        fv = e52.frame_and_voxels(S32, floor32, truth, s)
        scale = sm.scale_log_ratio(S32, truth)
    assert abs(d2a - d2f) <= 1e-6
    assert fv["frame_dev"] <= 1e-5
    assert fv["dice_arm"] == fv["dice_floor"] or abs(fv["dice_arm"] - fv["dice_floor"]) <= e52.dice_bound(
        fv["v"], fv["occupied"])
    assert scale == pytest.approx(np.log(t / tb.rms_radius(truth.astype(np.float64))), abs=1e-6)


@pytest.mark.parametrize("s", [0, 1])
def test_the_prefix_parity_is_the_floors_at_every_count(scored, s):
    floor32, _ = scored
    S32 = tb.copy_rows_sized(floor32, 3.25).astype(np.float32)
    for n in range(4, len(floor32) + 1):
        z = np.zeros((n, 3))
        assert ev.svd_parity(ev.match_rows(S32, z, s)[0]) == ev.svd_parity(ev.match_rows(floor32, z, s)[0]), n
    assert e52.parity_flips(S32, floor32, s) == []
    M = S32.copy()
    M[:, 2] = -M[:, 2]  # a mirror flips the parity at every count
    assert len(e52.parity_flips(M, floor32, s, n_max=40)) == 37


def test_the_voxel_bound_sees_changed_voxels():
    """A cloud moved by more than rounding changes voxels; the bound counts them and stays a bound."""
    floor32, _ = _cloud(2000, 3)
    truth = _cloud(2500, 4)[0]
    moved = floor32.copy()
    moved[:50] += np.float32(0.4)
    with ev.scorer_build("main"):
        fv = e52.frame_and_voxels(moved, floor32, truth, 0)
    assert fv["v"] > 0 and abs(fv["dice_arm"] - fv["dice_floor"]) <= e52.dice_bound(fv["v"], fv["occupied"])


# ------------------------------------------------------------------------------------------------ 3. the harness arm
VIEW = harness.TaskView(key="syn:E52", task="T2", kind="interp", setting="embryo", test_board="T2:embryo:test_interp",
                        roles=("before", "after"), role_days={"before": 7.5, "after": 8.0},
                        role_stems={"before": "E7.5", "after": "E8.0"}, reference_role="before", target_day=7.75)


@pytest.fixture(scope="module")
def stages():
    """A released earlier bracket with exact copies inserted (every 7th cell repeated after itself), the same bracket
    after the duplicate filter (with uns["qc"]["dropped_rows"]), a later bracket and the panel (as test_copy_rows)."""
    before, after, panel = _labelled_pair(np.random.default_rng(61))
    order, dropped = [], []
    for i in range(before.n_obs):
        order.append(i)
        if i % 7 == 3:
            dropped.append(len(order))
            order.append(i)
    released = before[order].copy()
    released.obs_names_make_unique()
    deduped = released[np.setdiff1d(np.arange(len(order)), dropped)].copy()
    deduped.uns["qc"] = {"drop_duplicates": True, "dropped_rows": np.asarray(dropped, np.int64)}
    return released, deduped, after, list(panel)


def _board(panel):
    return Board(key="pseudo:synthetic_embryo", task="T2", genes=panel, needs_coords=True, setting="embryo",
                 min_cells=10, max_cells=5000)


def _ctx(stages, seed=SEED):
    """The context run_task hands a needs_reference_cells candidate: the roles plus the reference cells."""
    rel, deduped, after, panel = stages
    inputs = {"before": deduped, "after": after, harness.REFERENCE_RELEASED: rel[:, panel],
              harness.REFERENCE_RELEASED_ROWS: np.arange(rel.n_obs)}
    return harness.BuildContext(task=VIEW, board=_board(panel), inputs=inputs, seed=seed)


def _roles_target(ctx) -> float:
    """fgw_shape0x's target: target_rms on a context holding only the task's roles."""
    return harness.target_rms(replace(ctx, inputs={r: ctx.inputs[r] for r in ctx.task.roles}))


def test_the_unamended_target_fails_on_the_reference_augmented_inputs(stages):
    """Why the specification was amended: target_rms on the inputs run_task hands this arm raises (the reference cells
    have no role day)."""
    with pytest.raises(KeyError):
        harness.target_rms(_ctx(stages))


@pytest.mark.parametrize("seed", [SEED, SEED + 4])
def test_the_sized_arm_is_the_copy_rows_arm_rescaled_to_the_roles_target(stages, seed):
    ctx = _ctx(stages, seed)
    unsized = harness._copy_rows(ctx, sized=False)
    sized = harness._copy_rows(ctx, sized=True)
    assert sized.X.dtype == unsized.X.dtype and sized.X.tobytes() == unsized.X.tobytes()
    T = _roles_target(ctx)
    want = tb.copy_rows_sized(unsized.coords, T).astype(np.float32)
    np.testing.assert_array_equal(_bits(sized.coords), _bits(want))
    assert abs(tb.rms_radius(np.asarray(sized.coords, np.float64)) / T - 1) <= 1e-6
    assert abs(T / tb.rms_radius(np.asarray(unsized.coords, np.float64)) - 1) > 1e-3  # the size really changes
    floor = ev.subsample(ev.Stage.from_anndata(stages[0], stages[3], need_coords=True), 0.1,
                         seed - harness.BUILD_SEED_OFFSET)  # the floor copy the read recomputes
    np.testing.assert_array_equal(_bits(sized.coords), _bits(e52.expected_coords(floor.coords, T)))
    built = CANDIDATES["fgw_copyrows_sized"].build(ctx)
    np.testing.assert_array_equal(_bits(built.coords), _bits(sized.coords))
    assert built.X.tobytes() == sized.X.tobytes()


def test_the_sized_arm_refuses_as_the_copy_rows_arm(stages):
    rel, deduped, after, panel = stages
    small = harness.BuildContext(task=VIEW, board=replace(_board(panel), min_cells=4000),
                                 inputs=_ctx(stages).inputs, seed=SEED)
    msgs = []
    for sized in (False, True):
        with pytest.raises(ValueError) as e:
            harness._copy_rows(small, sized=sized)
        msgs.append(str(e.value))
    assert msgs[0] == msgs[1] and e52.refusal_category(msgs[0]) == "count_outside_limits"
    bare = harness.BuildContext(task=VIEW, board=_board(panel), inputs={"before": deduped, "after": after}, seed=SEED)
    with pytest.raises(ValueError, match="reference cells"):
        CANDIDATES["fgw_copyrows_sized"].build(bare)


def test_the_arm_is_registered_as_specified():
    c = CANDIDATES["fgw_copyrows_sized"]
    assert c.needs_reference_cells and c.requires is CANDIDATES["fgw_copyrows"].requires
    assert c.tasks == frozenset({"T2"}) and c.params["stream"] == 51
    assert c.params["what"] == "fgw_copyrows rescaled about the copy centroid to the embryo size rule (row E52)"
    tasks = {t.key: t for ts in harness.load_tasks().values() for t in ts}
    for t in tasks.values():
        assert c.applies(t) == CANDIDATES["fgw_copyrows"].applies(t)
