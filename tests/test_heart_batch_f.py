"""The heart batch F (campaign branch campaign/heart-f): the Potts sorting of which cell sits where, carried from the
embryo's fgw_potts (row E31) to the heart's entry b0_trend: b0_trend's file to the last detail, then the drawn cells
(rows and obs["merged"] lineage labels) sorted among its fixed positions by vec.t2_build.potts_sort towards H*, with H*
the input stages' least-squares line read at the target's day (b0_trend_potts) or the file's own share plus the line's
slope over B0's time budget (b0_trend_potts_d). Synthetic heart stages only (tests/test_b0_trend.py's); no data file,
committed result or scorer is read."""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec import t2_build as tb
from vec.harness import B0_DAMPING, CANDIDATES, load_tasks

REPO = Path(__file__).resolve().parents[1]
ARMS = {"b0_trend_potts": "line", "b0_trend_potts_d": "damped"}
K = harness.NBHD_K


def _helpers(name: str):
    spec = importlib.util.spec_from_file_location(f"{name}_helpers_heart_f", REPO / "tests" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bt = _helpers("test_b0_trend")


@pytest.fixture(scope="module")
def emb_expr():
    return _helpers("test_emb_expr")


def _three() -> dict:
    return {"ref": bt._stage(1600, seed=2, shift=0.0, tiny=60), "start": bt._stage(1100, seed=3, shift=0.8, tiny=10),
            "earlier": bt._stage(900, seed=5, shift=1.6, tiny=10)}


@pytest.fixture(scope="module", params=["three", "four"])
def case(request):
    """(inputs, view, days): the synthetic heart task with three and with four input stages."""
    three = _three()
    if request.param == "three":
        return three, bt._view(), [8.5, 8.75, 9.5]
    four = {**three, "earliest": bt._stage(1000, seed=7, shift=2.2, tiny=10)}
    return four, bt._view(inputs=bt.FOUR_INPUTS), [8.3, 8.5, 8.75, 9.5]


def _ctx(case, seed: int = harness.BUILD_SEED_OFFSET + 1):
    inputs, view, _ = case
    return bt._ctx(inputs, view=view, seed=seed)


def _rows_multiset(X: np.ndarray) -> np.ndarray:
    X = np.ascontiguousarray(X)
    return X[np.lexsort(X.T[::-1])]


def _sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float32).tobytes()).hexdigest()


def _spy_build(ctx, form: str, monkeypatch):
    """The arm built with spies on the draw (its rows), potts_sort (labels, target, `at`) and the record."""
    seen: dict = {}
    real_draw, real_sort = harness.draw_copy, tb.potts_sort

    def draw(c, herd=False):
        out = real_draw(c, herd=herd)
        seen["rows"] = out[1].copy()
        return out

    def sort(coords, labels, target, k, rng, trace=None):
        at, rec = real_sort(coords, labels, target, k, rng, trace=trace)
        seen.update(coords=np.asarray(coords).copy(), labels=np.asarray(labels).copy(), target=target, k=k, at=at,
                    sort_rec=rec)
        return at, rec

    monkeypatch.setattr(harness, "draw_copy", draw)
    monkeypatch.setattr(tb, "potts_sort", sort)
    rec: dict = {}
    st = harness._b0_t2(ctx, trend=True, potts=form, potts_record=rec)
    monkeypatch.setattr(harness, "draw_copy", real_draw)
    monkeypatch.setattr(tb, "potts_sort", real_sort)
    return st, rec, seen


def _shares_by_hand(ctx, n: int) -> list[float]:
    """Each input stage's H by hand: a uniform draw of n of its cells without replacement on [seed, 29] (stages in
    day order, one generator), the same-label share of the 14 nearest other drawn cells."""
    rng = np.random.default_rng([ctx.seed, 29])
    out = []
    for r, _ in harness.trend_stages(ctx.task):
        a = ctx.inputs[r]
        C = np.asarray(a.obsm["spatial_3D"], dtype=float)[:, :3]
        lab = np.asarray(a.obs["merged"]).astype(str)
        rows = np.arange(len(C)) if len(C) <= n else rng.choice(len(C), n, replace=False)
        D = np.linalg.norm(C[rows][:, None] - C[rows][None], axis=2)
        np.fill_diagonal(D, np.inf)
        nb = np.argsort(D, axis=1, kind="stable")[:, :K - 1]
        out.append(float((lab[rows][nb] == lab[rows][:, None]).mean()))
    return out


# ================================================================================ the bases are unchanged
# fgw_potts on tests/test_emb_expr.py's synthetic pair (digest of X and coords), computed at the base fcc1f6b
FGW_POTTS_DIGEST = "e16d16d8f3ff6d681647ff9d239d49044e2a91b593810be177fad56ca5ec1e49"


def test_b0_trend_and_fgw_potts_are_unchanged_bit_for_bit(emb_expr):
    three = _three()
    st = CANDIDATES["b0_trend"].build(bt._ctx(three))
    assert _sha(st.X) == bt.B0_TREND_SHA256
    np.testing.assert_array_equal(harness._b0_t2(bt._ctx(three), trend=True).X, st.X)
    c = emb_expr.ctx(inputs=emb_expr.stages(), roles=("before", "after"))
    assert emb_expr.digest(CANDIDATES["fgw_potts"].build(c)) == FGW_POTTS_DIGEST


# ================================================================================ the arms
@pytest.mark.parametrize("name", list(ARMS))
def test_the_arm_is_b0_trends_file_with_its_cells_sorted_among_the_positions(case, name, monkeypatch):
    ctx = _ctx(case)
    base = CANDIDATES["b0_trend"].build(ctx)
    st, rec, seen = _spy_build(ctx, ARMS[name], monkeypatch)
    # positions, number and size: b0_trend's bit for bit
    np.testing.assert_array_equal(st.coords, base.coords)
    assert st.n == base.n and st.X.shape == base.X.shape
    # the expression rows: a permutation of b0_trend's, the permutation potts_sort's `at`
    at = seen["at"]
    assert sorted(at.tolist()) == list(range(base.n))
    np.testing.assert_array_equal(st.X, base.X[at])
    np.testing.assert_array_equal(_rows_multiset(st.X), _rows_multiset(base.X))
    assert rec["accepted_swaps"] > 0 and not np.array_equal(st.X, base.X)  # the idea is not idle here
    # the labels are the drawn cells' obs["merged"] (Untracked kept as a label) and move with the rows
    merged = np.asarray(case[0]["ref"].obs["merged"]).astype(str)
    np.testing.assert_array_equal(seen["labels"], merged[seen["rows"]])
    assert "Untracked" in set(seen["labels"])
    np.testing.assert_array_equal(seen["coords"], base.coords)
    assert seen["k"] == K
    assert tb.same_group_share(base.coords, seen["labels"][at], K) == rec["H_final"]
    assert tb.same_group_share(base.coords, seen["labels"], K) == rec["H_file"] == seen["sort_rec"]["H_start"]
    # the final H is not farther from the target than the start
    assert abs(rec["H_final"] - rec["H_target"]) <= abs(rec["H_file"] - rec["H_target"])
    assert rec["stop"] in tb.POTTS_STOPS


@pytest.mark.parametrize("name", list(ARMS))
def test_the_line_and_the_target_by_hand(case, name, monkeypatch):
    ctx = _ctx(case)
    days = case[2]
    _, rec, seen = _spy_build(ctx, ARMS[name], monkeypatch)
    n = len(seen["labels"])
    h = _shares_by_hand(ctx, n)
    assert rec["days"] == days and rec["stages"] == [r for r, _ in harness.trend_stages(ctx.task)]
    np.testing.assert_allclose(rec["H_stages"], h, atol=1e-15, rtol=0)
    d = np.asarray(days)
    A = np.vstack([d, np.ones_like(d)]).T
    slope, intercept = np.linalg.lstsq(A, np.asarray(h), rcond=None)[0]
    assert rec["slope"] == pytest.approx(slope, abs=1e-12) and rec["intercept"] == pytest.approx(intercept, abs=1e-12)
    if ARMS[name] == "line":
        want = intercept + slope * ctx.task.target_day
    else:
        want = rec["H_file"] + slope * B0_DAMPING * (9.5 - 8.75)
    assert rec["H_target"] == pytest.approx(want, abs=1e-12) and seen["target"] == rec["H_target"]


def test_potts_line_by_hand_and_its_refusals():
    s, i = harness.potts_line([1.0, 3.0], [0.2, 0.6])
    assert s == pytest.approx(0.2) and i == pytest.approx(0.0, abs=1e-15)  # two stages: the line through them
    s, i = harness.potts_line([0.0, 1.0, 2.0], [0.0, 1.0, 1.0])
    assert s == pytest.approx(0.5) and i == pytest.approx(1.0 / 6.0)
    with pytest.raises(ValueError, match="two stages"):
        harness.potts_line([1.0], [0.3])
    with pytest.raises(ValueError, match="two different days"):
        harness.potts_line([1.0, 1.0], [0.3, 0.4])


# ================================================================================ idle cases
@pytest.mark.parametrize("name", list(ARMS))
def test_with_the_files_h_at_the_target_the_arm_is_b0_trend_bit_for_bit(case, name, monkeypatch):
    """Every stage's H set to H of the file as drawn: the slope is 0, both targets equal H of the file exactly,
    nothing is drawn on the swap stream and the output is b0_trend's bit for bit."""
    ctx = _ctx(case)
    base = CANDIDATES["b0_trend"].build(ctx)
    first: dict = {}
    harness._b0_t2(ctx, trend=True, potts=ARMS[name], potts_record=first)
    h = first["H_file"]
    monkeypatch.setattr(tb, "bracket_share", lambda *a, **kw: h)
    swap_draws = []
    real_rng = np.random.default_rng

    def rng(seed=None):
        g = real_rng(seed)
        if isinstance(seed, list) and seed[-1] == harness.HEART_POTTS_STREAMS[1]:
            swap_draws.append(g)
        return g

    monkeypatch.setattr(harness.np.random, "default_rng", rng)
    rec: dict = {}
    st = harness._b0_t2(ctx, trend=True, potts=ARMS[name], potts_record=rec)
    monkeypatch.setattr(harness.np.random, "default_rng", real_rng)
    assert rec["slope"] == 0.0 and rec["H_target"] == h and rec["stop"] == "at_target" and rec["accepted_swaps"] == 0
    (g,) = swap_draws
    assert g.bit_generator.state == real_rng([ctx.seed, harness.HEART_POTTS_STREAMS[1]]).bit_generator.state
    np.testing.assert_array_equal(st.X, base.X)
    np.testing.assert_array_equal(st.coords, base.coords)


def test_with_a_zero_slope_the_damped_arm_is_b0_trend_bit_for_bit_whatever_the_intercept(case, monkeypatch):
    ctx = _ctx(case)
    base = CANDIDATES["b0_trend"].build(ctx)
    monkeypatch.setattr(harness, "potts_line", lambda days, shares: (0.0, 0.123))
    rec: dict = {}
    st = harness._b0_t2(ctx, trend=True, potts="damped", potts_record=rec)
    assert rec["H_target"] == rec["H_file"] and rec["stop"] == "at_target"
    np.testing.assert_array_equal(st.X, base.X)
    line: dict = {}
    harness._b0_t2(ctx, trend=True, potts="line", potts_record=line)  # the same zero slope: line's target is 0.123
    assert line["H_target"] == 0.123 and line["accepted_swaps"] > 0


def test_a_target_outside_0_1_runs_until_it_stalls_and_a_non_finite_one_refuses(case, monkeypatch):
    ctx = _ctx(case)
    monkeypatch.setattr(harness, "potts_line", lambda days, shares: (0.0, 5.0))
    rec: dict = {}
    harness._b0_t2(ctx, trend=True, potts="line", potts_record=rec)
    assert rec["H_target"] == 5.0 and rec["stop"] == "stalled" and rec["H_final"] > rec["H_file"]
    monkeypatch.setattr(harness, "potts_line", lambda days, shares: (float("nan"), 0.5))
    for form in ("line", "damped"):
        with pytest.raises(ValueError, match="finite"):
            harness._b0_t2(ctx, trend=True, potts=form)


def test_the_arms_need_the_trend_and_a_known_form():
    ctx = bt._ctx(_three())
    with pytest.raises(ValueError, match="needs trend"):
        harness._b0_t2(ctx, potts="line")
    with pytest.raises(ValueError, match="needs trend"):
        harness._b0_t2(ctx, trend=True, potts="other")
    with pytest.raises(ValueError, match="potts_record"):
        harness._b0_t2(ctx, trend=True, potts_record={})
    two = {k: v for k, v in _three().items() if k != "earlier"}
    with pytest.raises(ValueError, match="three"):
        CANDIDATES["b0_trend_potts"].build(bt._ctx(two, view=bt._view(inputs={"ref": "T2/E9.5.h5ad",
                                                                              "start": "T2/E8.75.h5ad"})))


# ================================================================================ seeds and streams
@pytest.mark.parametrize("name", list(ARMS))
def test_same_seed_same_output_another_seed_another_draw_of_the_stages_h(case, name):
    ctx = _ctx(case)
    a, b = CANDIDATES[name].build(ctx), CANDIDATES[name].build(ctx)
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.coords, b.coords)
    r1: dict = {}
    harness._b0_t2(ctx, trend=True, potts=ARMS[name], potts_record=r1)
    r2: dict = {}
    harness._b0_t2(dataclasses.replace(ctx, seed=ctx.seed + 1), trend=True, potts=ARMS[name], potts_record=r2)
    assert r1["H_stages"] != r2["H_stages"]


def test_the_two_streams_are_new_and_differ_from_every_stream_the_build_uses():
    used = {*range(16), *range(20, 29), harness.NBHD_STREAM, harness.SERIES_XT_STREAM, *harness.POTTS_STREAMS,
            *tb.HERD_STREAMS, *tb.MIDPOINT_STREAMS}
    a, b = harness.HEART_POTTS_STREAMS
    assert (a, b) == (29, 30) and a != b and not {a, b} & used
    seed = harness.BUILD_SEED_OFFSET + 1
    draws = {s: np.random.default_rng([seed, s]).random(64) for s in sorted(used | {a, b})}
    for s in sorted(used):
        assert not np.array_equal(draws[a], draws[s]) and not np.array_equal(draws[b], draws[s])
    assert not np.array_equal(draws[a], draws[b])
    assert not np.array_equal(draws[a], np.random.default_rng(seed).random(64))


def test_the_build_draws_the_stages_and_the_swaps_on_the_two_streams(case, monkeypatch):
    ctx = _ctx(case)
    seen = []
    real_share, real_sort = tb.bracket_share, tb.potts_sort

    def share(coords, labels, n, k, rng):
        seen.append(("stage", rng.bit_generator.state["state"]["state"]))
        return real_share(coords, labels, n, k, rng)

    def sort(coords, labels, target, k, rng, trace=None):
        seen.append(("swap", rng.bit_generator.state["state"]["state"]))
        return real_sort(coords, labels, target, k, rng, trace=trace)

    monkeypatch.setattr(tb, "bracket_share", share)
    monkeypatch.setattr(tb, "potts_sort", sort)
    harness._b0_t2(ctx, trend=True, potts="line")
    stages = [s for kind, s in seen if kind == "stage"]
    (swap,) = [s for kind, s in seen if kind == "swap"]
    assert len(stages) == len(case[2])
    assert stages[0] == np.random.default_rng([ctx.seed, 29]).bit_generator.state["state"]["state"]
    assert len(set(stages)) == len(stages)  # one generator, advanced from stage to stage
    assert swap == np.random.default_rng([ctx.seed, 30]).bit_generator.state["state"]["state"]


# ================================================================================ record, check line, registration
@pytest.mark.parametrize("name", list(ARMS))
def test_the_record_is_kept_with_the_trend_arms_and_printed_as_one_check_line(name, capsys):
    ctx = bt._ctx(_three())
    capsys.readouterr()
    CANDIDATES[name].build(ctx)
    out = [ln for ln in capsys.readouterr().out.splitlines() if ln.startswith(f"{name} check:")]
    assert len(out) == 1
    rec = harness.TREND_ARM_RECORDS[name]
    assert rec["task"] == ctx.task.key and rec["seed"] == ctx.seed and rec["form"] == ARMS[name]
    assert set(rec) == {"task", "seed", "form", "rule", "stages", "days", "H_stages", "slope", "intercept", "H_target",
                        "H_file", "H_final", "accepted_swaps", "stop"}
    line = out[0]
    for key in ("H_stages=", "slope=", "intercept=", "H_target=", "H_file=", "H_final=", "accepted_swaps=", "stop="):
        assert key in line
    assert f"accepted_swaps={rec['accepted_swaps']} stop={rec['stop']}" in line


def test_registration_is_b0_trends():
    tasks = {t.key: t for ts_ in load_tasks().values() for t in ts_}
    base = CANDIDATES["b0_trend"]
    for name, form in ARMS.items():
        c = CANDIDATES[name]
        assert c.tasks == base.tasks == frozenset({"T2"}) and c.needs == base.needs and c.derived == ()
        assert c.requires is base.requires is harness._heart_extrap_with_start
        assert {k for k, t in tasks.items() if c.applies(t)} == {"T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
        assert c.params["scope"] == base.params["scope"] == "chamber"
        assert c.params["form"] == form and c.params["streams"] == [29, 30] and c.params["k"] == K
        assert c.params["label"] == "merged"


# ================================================================================ review of b62db9a
# (1) Which scored quantities a permutation of cells over fixed positions can change. The pinned scorer subsamples the
# prediction's rows by index in variogram and energy_distance (1500 rows) and mmd_u (2000 rows), so above those sizes
# a permutation re-indexes the rows they read: "every metric that reads expression without positions equals
# b0_trend's exactly" is false there. The check below is the corrected one.
EXACT = ("d2_shape", "occupancy_dice", "scale_log_ratio", "sliced_wasserstein", "count_log_ratio")
ROUNDED = {"de_score": 4, "de_direction": 4, "pb_rel_err": 4, "pseudobulk_pearson": 4, "library_size_ratio": 3,
           "variance_ratio": 3, "composition_JSD": 4}  # decimals the pinned T2 scorer rounds them to
SUBSAMPLED = {"variogram": 1500, "energy_distance": 1500, "mmd_u": 2000}  # rows the scorer draws from the prediction
BY_DESIGN = ("neighborhood_mmd", "morans_I_agreement")


def _validity_check(base, arm, mb: dict, ma: dict) -> set[str]:
    """The corrected validity check of an arm against its base (same positions, cells permuted): (a) coords bytes and
    n identical, (b) the lexsorted row multiset identical, (c) the shape metrics equal exactly, (d) the order-invariant
    expression metrics equal within one unit of their rounding, (e) the subsampled ones equal where the file holds no
    more rows than the scorer draws. Returns the subsampled metrics that differ (allowed only above their size)."""
    assert np.asarray(base.coords).tobytes() == np.asarray(arm.coords).tobytes() and base.n == arm.n
    assert _sha(_rows_multiset(base.X)) == _sha(_rows_multiset(arm.X))
    for m in EXACT:
        assert mb[m] == ma[m], m
    for m, d in ROUNDED.items():
        if mb[m] is None or ma[m] is None:
            assert mb[m] is None and ma[m] is None, m
        else:
            assert abs(mb[m] - ma[m]) <= 10.0 ** -d * (1 + 1e-9), m
    changed = {m for m in SUBSAMPLED if mb[m] != ma[m]}
    for m in changed:
        assert base.n > SUBSAMPLED[m], m
    for m in BY_DESIGN:
        assert m in mb and m in ma
    return changed


def _metrics(st, truth, ref):
    from vec.evaluate import panel_metrics

    return panel_metrics("T2", st, truth, ref, seed=0)


def _cloud(n: int, seed: int):
    from vec.evaluate import Stage

    r = np.random.default_rng(seed)
    C = r.normal(size=(n, 3)).astype(np.float32)
    g = len(bt.GENES)
    X = np.log1p(r.poisson(np.exp(C[:, :1] * 0.8 + r.normal(size=(1, g))), size=(n, g))).astype(np.float32)
    return Stage(X=X, coords=C, celltype=np.array(["NA"] * n))


def test_the_arms_pass_the_corrected_validity_check_against_b0_trend():
    """At the synthetic file's size (at most 1500 rows) every expression-only metric is b0_trend's."""
    ctx = bt._ctx(_three())
    truth, ref = _cloud(900, 1), _cloud(900, 2)
    base = CANDIDATES["b0_trend"].build(ctx)
    assert base.n <= min(SUBSAMPLED.values())
    mb = _metrics(base, truth, ref)
    for name in ARMS:
        st = CANDIDATES[name].build(ctx)
        ma = _metrics(st, truth, ref)
        assert _validity_check(base, st, mb, ma) == set()
        assert mb["neighborhood_mmd"] != ma["neighborhood_mmd"]  # the arm's idea moves the local organisation


def test_above_the_scorers_subsample_sizes_a_permutation_changes_the_subsampled_metrics():
    """The brief's premise (expression-only metrics equal exactly) fails above 1500 rows: variogram and energy_distance
    change by subsample re-indexing; mmd_u, which draws 2000 rows, still reads every row at 1800."""
    from vec.evaluate import Stage

    truth, ref = _cloud(1500, 1), _cloud(1500, 2)
    p = _cloud(1800, 3)
    q = Stage(X=p.X[np.random.default_rng(9).permutation(p.n)], coords=p.coords, celltype=p.celltype)
    mp, mq = _metrics(p, truth, ref), _metrics(q, truth, ref)
    assert _validity_check(p, q, mp, mq) == {"variogram", "energy_distance"}


# (2) The two arms share one swap path (stream [seed, 30]) and differ only in where they stop.
def _traced(ctx, form: str, line: tuple[float, float], monkeypatch):
    """The arm built with a fixed line (slope, intercept) and potts_sort's trace kept: (stage, record, trace)."""
    trace: list = []
    real_sort, real_line = tb.potts_sort, harness.potts_line

    def sort(coords, labels, target, k, rng, trace_=None):
        return real_sort(coords, labels, target, k, rng, trace=trace)

    monkeypatch.setattr(tb, "potts_sort", sort)
    monkeypatch.setattr(harness, "potts_line", lambda days, shares: line)
    rec: dict = {}
    st = harness._b0_t2(ctx, trend=True, potts=form, potts_record=rec)
    monkeypatch.setattr(tb, "potts_sort", real_sort)
    monkeypatch.setattr(harness, "potts_line", real_line)
    return st, rec, trace


def _small_ref() -> dict:
    """A reference stage with fewer cells than the file: the file repeats rows (drawn with replacement, jittered)."""
    return {**_three(), "ref": bt._stage(500, seed=2, shift=0.0, tiny=60)}


@pytest.fixture(scope="module", params=["three", "four", "small_ref"])
def nest_case(request):
    if request.param == "small_ref":
        return _small_ref(), bt._view(), [8.5, 8.75, 9.5]
    three = _three()
    if request.param == "three":
        return three, bt._view(), [8.5, 8.75, 9.5]
    four = {**three, "earliest": bt._stage(1000, seed=7, shift=2.2, tiny=10)}
    return four, bt._view(inputs=bt.FOUR_INPUTS), [8.3, 8.5, 8.75, 9.5]


def _gap(ctx) -> float:
    t = ctx.task
    return float(t.role_days[t.reference_role]) - float(t.role_days["start"])


def test_the_two_arms_are_nested_on_one_swap_path(nest_case, monkeypatch):
    """Both targets above H of the file (damped H_file + 0.03, line H_file + 0.08): the nearer target's accepted swaps
    are the farther one's until the first swap the nearer one rejects as an overshoot (the proposals are the same), and
    the nearer one stops there or earlier. One target reached and the other farther: the files differ."""
    ctx = _ctx(nest_case)
    first: dict = {}
    harness._b0_t2(ctx, trend=True, potts="damped", potts_record=first)
    h0 = first["H_file"]
    s = 0.03 / (B0_DAMPING * _gap(ctx))
    line = (s, h0 + 0.08 - s * float(ctx.task.target_day))
    st_d, rec_d, tr_d = _traced(ctx, "damped", line, monkeypatch)
    st_l, rec_l, tr_l = _traced(ctx, "line", line, monkeypatch)
    near, far = rec_d["H_target"], rec_l["H_target"]
    assert rec_d["H_file"] == rec_l["H_file"] == h0 and h0 < near < far
    m = 0
    while m < min(len(tr_d), len(tr_l)) and np.array_equal(tr_d[m][0], tr_l[m][0]):
        m += 1
    if m < len(tr_d):  # a divergence: the farther run accepted a swap the nearer one rejected as an overshoot
        h_prev = tr_d[m - 1][1] if m else h0
        assert m < len(tr_l) and not abs(tr_l[m][1] - near) < abs(h_prev - near)
        assert tr_l[m][1] >= near
    assert rec_d["stop"] == "reached" and rec_l["accepted_swaps"] > rec_d["accepted_swaps"]
    assert _sha(st_d.X) != _sha(st_l.X)
    np.testing.assert_array_equal(st_d.coords, st_l.coords)


def test_with_both_targets_beyond_the_stall_the_two_arms_are_the_same_file(nest_case, monkeypatch):
    ctx = _ctx(nest_case)
    line = (10.0, 0.0)  # damped: H_file + 10 x B0_DAMPING x the gap; line: 10 x the target day; both far above 1
    st_d, rec_d, tr_d = _traced(ctx, "damped", line, monkeypatch)
    st_l, rec_l, tr_l = _traced(ctx, "line", line, monkeypatch)
    assert rec_d["H_target"] > 1.0 and rec_l["H_target"] > 1.0 and rec_d["H_target"] != rec_l["H_target"]
    assert rec_d["stop"] == rec_l["stop"] == "stalled"
    assert _sha(st_d.X) == _sha(st_l.X) and rec_d["H_final"] == rec_l["H_final"]
    assert rec_d["accepted_swaps"] == rec_l["accepted_swaps"] == len(tr_d) == len(tr_l)
    assert all(np.array_equal(a[0], b[0]) for a, b in zip(tr_d, tr_l, strict=True))


def test_a_reference_smaller_than_the_file(monkeypatch):
    """The reference holds fewer cells than the file: b0_trend repeats rows; the arm is still a permutation of
    b0_trend's rows over its positions, and the reference stage's H is read on all its cells (no draw)."""
    inputs = _small_ref()
    ctx = bt._ctx(inputs)
    base = CANDIDATES["b0_trend"].build(ctx)
    for name, form in ARMS.items():
        st, rec, seen = _spy_build(ctx, form, monkeypatch)
        assert base.n > inputs["ref"].n_obs and len(np.unique(seen["rows"])) < base.n
        np.testing.assert_array_equal(st.coords, base.coords)
        np.testing.assert_array_equal(st.X, base.X[seen["at"]])
        np.testing.assert_array_equal(_rows_multiset(st.X), _rows_multiset(base.X))
        i = rec["stages"].index("ref")
        a = inputs["ref"]
        want = tb.same_group_share(np.asarray(a.obsm["spatial_3D"], dtype=float)[:, :3],
                                   np.asarray(a.obs["merged"]).astype(str), K)
        assert rec["H_stages"][i] == want, name
        np.testing.assert_allclose(rec["H_stages"], _shares_by_hand(ctx, base.n), atol=1e-15, rtol=0)
        assert abs(rec["H_final"] - rec["H_target"]) <= abs(rec["H_file"] - rec["H_target"])
