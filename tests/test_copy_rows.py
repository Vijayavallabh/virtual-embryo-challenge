"""Row E50 as specified (the captain's order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (h)): the copy-rows
file, vec.t2_build.copy_rows_cloud and the harness arms fgw_copyrows and fgw_copyrows_series, on synthetic arrays and
stages only. With the option off every output is what it was."""

from __future__ import annotations

import itertools
import re
from dataclasses import replace

import numpy as np
import pytest
from conftest import GENES, simulate_stage
from scipy.spatial.distance import cdist
from test_t2_build import _labelled_pair

from vec import evaluate as ev
from vec import harness
from vec import t2_build as tb
from vec.boards import RANKED_METRICS, Board
from vec.core import draw_balanced
from vec.evaluate import Stage
from vec.harness import CANDIDATES
from vec.io import COORD_KEY

SEED = 1_000_003


def _bits(C):
    return np.ascontiguousarray(np.asarray(C, dtype=np.float32)).view(np.uint32)


def _world(rng, n=400, n_after=300, n_genes=8, n_dup=6):
    """A released earlier bracket of n rows with n_dup exact copies (a copy repeats the row before it: position and
    expression), the rows the duplicate filter drops, the bracket after the filter, and a later bracket."""
    C = rng.normal(size=(n, 3)).astype(np.float32)
    X = rng.gamma(2.0, 1.0, size=(n, n_genes)).astype(np.float32)
    dup = np.sort(rng.choice(np.arange(1, n), n_dup, replace=False)) if n_dup else np.array([], np.int64)
    for j in dup:
        C[j], X[j] = C[j - 1], X[j - 1]
    keep = np.setdiff1d(np.arange(n), dup)
    Ea = (rng.gamma(2.0, 1.0, size=(n_after, n_genes)) + 0.3).astype(np.float32)
    return C, X, dup, X[keep], Ea


def _copy_index(n, copy_seed, frac=0.1):
    st = ev.subsample(Stage(X=np.zeros((n, 0), np.float32), coords=np.zeros((n, 3), np.float32),
                            celltype=np.arange(n)), frac, copy_seed)
    return np.asarray(st.celltype)


# ------------------------------------------------------------------------------------------------ 1. coordinates
@pytest.mark.parametrize("n", [45, 125, 400, 1234])
@pytest.mark.parametrize("copy_seed", [0, 3])
def test_the_coordinates_are_the_floor_copy_bit_for_bit(n, copy_seed):
    rng = np.random.default_rng(n)
    C, X, dup, Eb, Ea = _world(rng, n=n, n_dup=2 if n < 100 else 12)
    assert len(np.unique(C, axis=0)) < n  # duplicate positions in the reference
    res = tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, copy_seed, SEED, n_anchor=60)
    copy = ev.subsample(Stage(X=X, coords=C), 0.1, copy_seed)
    assert res.coords.dtype == np.float64 and res.coords.shape == copy.coords.shape
    np.testing.assert_array_equal(_bits(res.coords), _bits(copy.coords))  # rows, order, count, float32 bits
    assert len(res.rows) == max(2, int(round(n * 0.1)))
    assert not any(isinstance(v, int) and not isinstance(v, bool) and v == len(res.rows) for v in res.record.values())


def test_one_count_ends_in_a_half():
    assert any(abs((0.1 * n) % 1 - 0.5) < 1e-9 for n in (45, 125, 400, 1234))


# ------------------------------------------------------------------------------------------------ 2. the scorer
@pytest.fixture(scope="module")
def scored_world():
    rng = np.random.default_rng(5)
    means = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    later = means.copy()
    later[:, :8] *= 3.0
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float)
    ref = simulate_stage(3000, means, [0.4, 0.35, 0.25], rng, coords_scale=1.5, coords_shift=centres)
    ref.coords[1::97] = ref.coords[0:-1:97][: len(ref.coords[1::97])]  # some exact duplicate positions
    targets = {n: simulate_stage(n, later, [0.3, 0.4, 0.3], rng, coords_scale=1.8, coords_shift=centres * 1.3)
               for n in (8000, 4000)}  # truth halves of 400 and 200 against a copy of 300
    board = Board(key="pseudo:copy_rows", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  split="pseudo", min_cells=10, max_cells=20000)
    return ref, targets, board


@pytest.mark.parametrize("build", ["main", "prefix"])
@pytest.mark.parametrize("s", [0, 1])
def test_the_position_terms_are_the_floors_under_both_builds(scored_world, build, s):
    ref, targets, board = scored_world
    rng = np.random.default_rng(s + 10)
    Ea = rng.gamma(2.0, 1.0, size=(500, len(GENES))).astype(np.float32)
    res = tb.copy_rows_cloud(ref.coords, [], ref.X, Ea, ref.X, 0.5, s, SEED + s, n_anchor=200)
    X = rng.gamma(2.0, 1.0, size=(len(res.rows), len(GENES))).astype(np.float32)  # expression unlike the floor's
    pred = Stage(X=X, coords=res.coords.astype(np.float32))
    for n_target, target in targets.items():
        r = ev.score_board(board, pred, target, ref, seed=s, anchors="computed", build=build, protocol="bake")
        floor = r.extras["floor_metrics"]
        terms = ["d2_shape", "occupancy_dice", "scale_log_ratio"] + [m for m in ("sliced_wasserstein",)
                                                                     if m in floor and m in r.metrics]
        for m in terms:
            assert r.metrics[m] == floor[m], (n_target, m)
        assert r.metrics["mmd_u"] != floor["mmd_u"]  # the expression is not the floor's
        assert r.extras["handedness"]["prefix_mode"] == r.extras["anchor_handedness"]["floor"]["prefix_mode"]
    assert {targets[n].n // 20 for n in targets} == {400, 200}


# ------------------------------------------------------------------------------------------------ 3. expression
def test_own_rows_and_later_rows():
    rng = np.random.default_rng(3)
    C, X, dup, Eb, Ea = _world(rng, n=600, n_dup=20)
    res = tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 2, SEED, n_anchor=80)
    idx = _copy_index(600, 2)
    keep = np.setdiff1d(np.arange(600), dup)
    b, a = res.stage == 0, res.stage == 1
    np.testing.assert_array_equal(keep[res.rows[b]], idx[b])  # each earlier slot keeps its own cell
    np.testing.assert_array_equal(Eb[res.rows[b]], X[idx[b]])
    n = len(idx)
    n_a = n - int(round(0.5 * n))
    assert int(a.sum()) == n_a
    want = draw_balanced(len(Ea), n_a, np.random.default_rng([SEED, 1]))
    assert sorted(res.rows[a].tolist()) == sorted(want.tolist())


@pytest.mark.parametrize("t, n_rel", [(0.5, 470), (1 / 3, 530), (0.5, 455)])
def test_the_split_is_interpolate_clouds(t, n_rel):
    rng = np.random.default_rng(4)
    C, X, dup, Eb, Ea = _world(rng, n=n_rel, n_dup=0)
    res = tb.copy_rows_cloud(C, dup, Eb, Ea, X, t, 0, SEED, n_anchor=60)
    n = len(res.rows)
    ic = tb.interpolate_clouds(C, rng.normal(size=(len(Ea), 3)), Eb, Ea, t, n, SEED, target_rms=1.0, n_anchor=60,
                               k=4)
    assert (int((res.stage == 0).sum()), int((res.stage == 1).sum())) == (ic.record["n_before"],
                                                                          ic.record["n_after"])


def test_the_assignment_is_the_brute_force_minimum_on_six():
    rng = np.random.default_rng(6)
    C, X, dup, Eb, Ea = _world(rng, n=120, n_dup=0)
    res = tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 1, SEED, n_anchor=50, n_pc=4)
    idx = _copy_index(120, 1)
    slots = np.flatnonzero(res.stage == 1)
    assert len(slots) == 6
    rows_a = draw_balanced(len(Ea), 6, np.random.default_rng([SEED, 1]))
    r4 = np.random.default_rng([SEED, 4])
    ib = np.sort(r4.choice(len(Eb), min(50, len(Eb)), replace=False))
    ia = np.sort(r4.choice(len(Ea), min(50, len(Ea)), replace=False))
    centre, L = tb._joint_pc_loadings(Eb[ib], Ea[ia], 4)
    cost = cdist((Ea[rows_a] - centre) @ L, (X[idx[slots]] - centre) @ L, "sqeuclidean")
    best = min(itertools.permutations(range(6)), key=lambda p: sum(cost[i, p[i]] for i in range(6)))
    want = np.empty(6, np.int64)
    for i, j in enumerate(best):
        want[j] = rows_a[i]
    np.testing.assert_array_equal(res.rows[slots], want)


def test_stream_51_is_stable_at_a_seed_and_moves_with_it():
    rng = np.random.default_rng(7)
    C, X, dup, Eb, Ea = _world(rng, n=400, n_dup=0)
    one = tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 0, SEED, n_anchor=60)
    two = tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 0, SEED, n_anchor=60)
    other = tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 0, SEED + 1, n_anchor=60)
    s1, s2, s3 = (np.flatnonzero(r.stage == 1) for r in (one, two, other))
    np.testing.assert_array_equal(s1, s2)
    assert not np.array_equal(s1, s3)
    n = len(one.rows)
    want = np.sort(np.random.default_rng([SEED, tb.COPY_ROWS_STREAM]).choice(np.arange(n), n - round(0.5 * n),
                                                                              replace=False))
    np.testing.assert_array_equal(s1, want)  # no dropped slot: the later slots are stream 51's draw alone
    assert tb.COPY_ROWS_STREAM == 51 and 51 not in (tb.FLOOR_EXACT_STREAM, tb.COPY_PARITY_STREAM, 50)


# ------------------------------------------------------------------------------------------------ 4. dropped slots
def test_every_dropped_slot_takes_a_later_cell():
    rng = np.random.default_rng(8)
    C, X, dup, Eb, Ea = _world(rng, n=800, n_dup=120)
    res = tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 0, SEED, n_anchor=80)
    idx = _copy_index(800, 0)
    dropped_slots = np.flatnonzero(np.isin(idx, dup))
    assert len(dropped_slots) > 0
    assert (res.stage[dropped_slots] == 1).all()


def test_refusals_show_no_value():
    rng = np.random.default_rng(9)
    C, X, dup, Eb, Ea = _world(rng, n=400, n_dup=330)  # most rows are copies: more dropped slots than half
    with pytest.raises(ValueError) as e:
        tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 0, SEED, n_anchor=60)
    assert "share" in str(e.value) and not re.search(r"\d", str(e.value))
    C, X, dup, Eb, Ea = _world(rng, n=20, n_dup=0)
    with pytest.raises(ValueError) as e:  # a copy of two rows
        tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 0, SEED)
    assert not re.search(r"\d", str(e.value))
    C, X, dup, Eb, Ea = _world(rng, n=200, n_dup=0)
    C[5, 0] = np.nan
    with pytest.raises(ValueError) as e:
        tb.copy_rows_cloud(C, dup, Eb, Ea, X, 0.5, 0, SEED)
    assert "non-finite" in str(e.value) and not re.search(r"\d", str(e.value))


# ------------------------------------------------------------------------------------------------ the harness
VIEW = harness.TaskView(key="syn:E50", task="T2", kind="interp", setting="embryo", test_board="T2:embryo:test_interp",
                        roles=("before", "after"), role_days={"before": 7.5, "after": 8.0},
                        role_stems={"before": "E7.5", "after": "E8.0"}, reference_role="before", target_day=7.75)


@pytest.fixture(scope="module")
def stages():
    """A released earlier bracket with exact copies inserted (every 7th cell repeated after itself), the same bracket
    after the duplicate filter (with uns["qc"]["dropped_rows"]), a later bracket and the panel."""
    before, after, panel = _labelled_pair(np.random.default_rng(61))
    order, dropped = [], []
    for i in range(before.n_obs):
        order.append(i)
        if i % 7 == 3:
            dropped.append(len(order))
            order.append(i)
    released = before[order].copy()
    released.obs_names_make_unique()
    keep = np.setdiff1d(np.arange(len(order)), dropped)
    deduped = released[keep].copy()
    deduped.uns["qc"] = {"drop_duplicates": True, "dropped_rows": np.asarray(dropped, np.int64)}
    return released, deduped, after, list(panel)


def _board(panel, **kw):
    return Board(key="pseudo:synthetic_embryo", task="T2", genes=panel, needs_coords=True, setting="embryo",
                 min_cells=10, max_cells=5000, **kw)


def _ctx(stages, seed=SEED, rows=None, released=None):
    rel, deduped, after, panel = stages
    rel = rel if released is None else released
    rows = np.arange(rel.n_obs) if rows is None else rows
    inputs = {"before": deduped, "after": after, harness.REFERENCE_RELEASED: rel[:, panel],
              harness.REFERENCE_RELEASED_ROWS: rows}
    return harness.BuildContext(task=VIEW, board=_board(panel), inputs=inputs, seed=seed)


def test_the_arms_are_registered_off_by_default():
    for name in ("fgw_copyrows", "fgw_copyrows_series"):
        c = CANDIDATES[name]
        assert c.needs_reference_cells and c.params["stream"] == 51 and c.tasks == frozenset({"T2"})
    assert [n for n, c in CANDIDATES.items() if c.needs_reference_cells] == ["fgw_copyrows", "fgw_copyrows_series"]
    tasks = {t.key: t for ts in harness.load_tasks().values() for t in ts}
    for key in ("T2e/E6.75+E7.5->E7.25", "T2e/E7.25+E8.0->E7.5"):
        assert CANDIDATES["fgw_copyrows"].applies(tasks[key]) == CANDIDATES["fgw_shape0x"].applies(tasks[key])
        assert (CANDIDATES["fgw_copyrows_series"].applies(tasks[key])
                == CANDIDATES["fgw_shape0x_series"].applies(tasks[key]))


def test_the_arm_is_the_floor_copy_with_own_and_later_expression(stages):
    rel, deduped, after, panel = stages
    ctx = _ctx(stages)
    st = CANDIDATES["fgw_copyrows"].build(ctx)
    floor = ev.subsample(Stage.from_anndata(rel, panel, need_coords=True), 0.1, SEED - harness.BUILD_SEED_OFFSET)
    np.testing.assert_array_equal(_bits(st.coords), _bits(floor.coords))
    Xr = np.asarray(rel[:, panel].X.toarray(), np.float32)
    Ea = np.asarray(after[:, panel].X.toarray(), np.float32)
    idx = _copy_index(rel.n_obs, SEED - harness.BUILD_SEED_OFFSET)
    own = np.array([np.array_equal(st.X[i], Xr[idx[i]]) for i in range(len(idx))])
    later = np.array([any(np.array_equal(st.X[i], e) for e in Ea) for i in range(len(idx))])
    assert own.sum() == int(round(0.5 * len(idx))) and (own | later).all()
    dropped = np.isin(idx, deduped.uns["qc"]["dropped_rows"])
    assert dropped.any() and not own[dropped].any()  # a slot the filter drops never keeps its own row


def test_series_passes_the_stage_and_rows_through_unchanged(stages, monkeypatch):
    rel, deduped, after, panel = stages
    seen = {}

    def fake(before, after_, interp, genes, t, seed, **kw):
        seen.update(before=before, after=after_, stage=interp.stage.copy(), rows=interp.rows.copy(), t=t, seed=seed)
        return np.full((len(interp.rows), len(genes)), 0.5, np.float32), {}

    monkeypatch.setattr(tb, "series_rendered_expression", fake)
    CANDIDATES["fgw_copyrows_series"].build(_ctx(stages))
    Eb = np.asarray(deduped[:, panel].X.toarray(), np.float32)
    Ea = np.asarray(after[:, panel].X.toarray(), np.float32)
    Xr = np.asarray(rel[:, panel].X.toarray(), np.float32)
    res = tb.copy_rows_cloud(rel.obsm[COORD_KEY], deduped.uns["qc"]["dropped_rows"], Eb, Ea, Xr, 0.5,
                             SEED - harness.BUILD_SEED_OFFSET, SEED)
    np.testing.assert_array_equal(seen["stage"], res.stage)
    np.testing.assert_array_equal(seen["rows"], res.rows)
    assert seen["before"] is deduped and seen["after"] is after and seen["t"] == 0.5 and seen["seed"] == [SEED, 7]


@pytest.mark.parametrize("target_day", [7.75, 7.625, 7.875])
def test_the_arm_splits_at_the_tasks_own_t(stages, target_day):
    """Row E50, the analyst's decision (b): round((1 - t) n) slots keep their own cell's expression, t the task's own
    time fraction as fgw_shape0x uses it (half at t = 0.5); the rule texts say the same."""
    rel, deduped, after, panel = stages
    ctx = replace(_ctx(stages), task=replace(VIEW, target_day=target_day))
    st = CANDIDATES["fgw_copyrows"].build(ctx)
    t = (target_day - 7.5) / (8.0 - 7.5)
    Xr = np.asarray(rel[:, panel].X.toarray(), np.float32)
    idx = _copy_index(rel.n_obs, SEED - harness.BUILD_SEED_OFFSET)
    own = sum(np.array_equal(st.X[i], Xr[idx[i]]) for i in range(len(idx)))
    assert own == int(round((1 - t) * len(idx)))
    assert "round((1 - t) n)" in tb.COPY_ROWS_RULE
    assert "round((1 - t) n) slots (half at t = 0.5)" in CANDIDATES["fgw_copyrows"].params["what"]


@pytest.mark.parametrize("name", ["fgw_copyrows", "fgw_copyrows_series"])
def test_the_arms_call_no_cloud_guard(stages, monkeypatch, name):
    """Row E50, the analyst's decision: no cloud guard on the copy-rows path (reading (i) checks the coordinates bit
    for bit); a patched guard that raises is never reached."""
    def never(*a, **kw):
        raise AssertionError("the cloud guard ran on the copy-rows path")

    monkeypatch.setattr(tb, "cloud_guard", never)
    if name.endswith("series"):
        monkeypatch.setattr(tb, "series_rendered_expression",
                            lambda before, after, interp, genes, t, seed, **kw: (
                                np.full((len(interp.rows), len(genes)), 0.5, np.float32), {}))
    st = CANDIDATES[name].build(_ctx(stages))
    assert len(st.coords) > 0


def test_without_the_handover_the_arm_refuses(stages):
    rel, deduped, after, panel = stages
    ctx = harness.BuildContext(task=VIEW, board=_board(panel), inputs={"before": deduped, "after": after}, seed=SEED)
    with pytest.raises(ValueError, match="reference cells"):
        CANDIDATES["fgw_copyrows"].build(ctx)


class _Score:
    def __init__(self):
        self.score, self.skills, self.metrics, self.group_skill = 50.0, {}, {}, {}
        self.anchors = {m: {"floor": 1.0, "ceiling": 0.0} for m in RANKED_METRICS["T2"]}
        self.extras = {"inverted_anchors": []}


def _task():
    return harness.PseudoTask(key="syn", test_board="T2:embryo:test_interp", task="T2", setting="embryo",
                              kind="interp", target="T2/E8.0.h5ad", reference="T2/E7.5.h5ad",
                              inputs={"before": "T2/E7.5.h5ad", "after": "T2/E8.25_late.h5ad"},
                              twin_board="T2:embryo:val_interp")


def test_with_a_within_cut_the_arm_is_score_seeds_floor_copy_and_others_see_nothing_new(stages, monkeypatch,
                                                                                         tmp_path):
    rel, deduped, after, panel = stages
    cut = np.sort(np.random.default_rng(12).choice(rel.n_obs, rel.n_obs - 150, replace=False))
    ref_cells = rel[cut][:, panel]
    ref = Stage.from_anndata(ref_cells, panel, need_coords=True)
    floors, preds, roles = {}, {}, []
    orig = harness.subsample

    def spy_subsample(stage, frac, seed):
        out = orig(stage, frac, seed)
        floors[seed] = out
        return out

    def stub_score(board, pred, target, ref_, seed=0, **kw):
        preds.setdefault(seed, []).append(pred)
        return _Score()

    def spy(ctx):
        roles.append(set(ctx.inputs))
        return Stage(X=np.zeros((5, len(panel)), np.float32), coords=np.zeros((5, 3), np.float32))

    monkeypatch.setattr(harness, "subsample", spy_subsample)
    monkeypatch.setattr(harness, "score_board", stub_score)
    other = harness.Candidate("other", frozenset({"T2"}), harness._ANY, spy)
    arm = CANDIDATES["fgw_copyrows"]
    board = _board(panel, truth_cells=100, ref_cells=ref.n, split="pseudo")
    harness.score_seeds(_task(), board, ref, ref, {"before": deduped, "after": after},
                        [CANDIDATES[harness.FLOOR], other, arm], seeds=(0, 1), root=tmp_path, data_dir=tmp_path,
                        reference_cells=ref_cells, reference_rows=cut)
    for s in (0, 1):
        built = preds[s][-1]  # the arm, scored last
        np.testing.assert_array_equal(_bits(built.coords), _bits(floors[s].coords))
    assert roles == [{"before", "after"}] * 2  # a candidate without the flag gets no extra input
    with pytest.raises(ValueError, match="needs the reference cells"):
        harness.score_seeds(_task(), board, ref, ref, {"before": deduped, "after": after},
                            [CANDIDATES[harness.FLOOR], arm], seeds=(0,), root=tmp_path, data_dir=tmp_path)


def test_run_task_hands_the_reference_cells_over_only_for_the_arms(stages, monkeypatch, tmp_path):
    rel, deduped, after, panel = stages
    t = _task()
    calls = []
    ref_a = rel[:, panel].copy()

    def read(path, genes):
        return ref_a

    monkeypatch.setattr(harness, "read_stage", read)
    monkeypatch.setattr(harness, "load_stage", lambda path, **kw: deduped if str(path).endswith("E7.5.h5ad")
                        else after)
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    board = _board(panel, split="pseudo")
    monkeypatch.setattr(harness, "task_board", lambda *a: board)

    def stub(t, board, target, ref, inputs, candidates, **kw):
        calls.append(kw)
        return []

    monkeypatch.setattr(harness, "score_seeds", stub)
    boards = {"T2:embryo:val_interp": replace(board, key="T2:embryo:val_interp", split="val")}
    harness.run_task(t, [CANDIDATES["copy"]], seeds=(0,), data_dir=tmp_path, boards=boards, log=None)
    harness.run_task(t, [CANDIDATES["copy"], CANDIDATES["fgw_copyrows"]], seeds=(0,), data_dir=tmp_path,
                     boards=boards, log=None)
    assert "reference_cells" not in calls[0] and "reference_rows" not in calls[0]
    assert calls[1]["reference_cells"] is ref_a
    np.testing.assert_array_equal(calls[1]["reference_rows"], np.arange(ref_a.n_obs))

