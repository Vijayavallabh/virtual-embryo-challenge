"""Heart batch b (campaign/heart-b): four arms on `b0_trend`, each with one idea read from the record and nothing else
changed: `b0_trend_smooth` (the top-up from kNN-smoothed expected counts; TRIES.md "Pruned before running",
b0_smooth), `b0_trend_tail` (the heart tail decoder; TRIES.md N3, N5), `b0_trend_pava` (R9's persistence gate and
combiner) and `b0_trend_r8` (R8's run basis; experiments/research_2026_11_runs/REPORT.md, R8 to R10). Synthetic data
built here only; nothing scored, no prediction compared with any target."""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace

import numpy as np
import pytest
from scipy import sparse
from test_b0_trend import (
    BOARD,
    FOUR_INPUTS,
    TYPES,
    _ctx,
    _h,
    _only_start_and_ref,
    _stage,
    _view,
)

from vec import core, harness, t1_scope
from vec.harness import B0_DAMPING, CANDIDATES, load_tasks

ARMS = {"b0_trend_smooth": "smooth", "b0_trend_tail": "tail", "b0_trend_pava": "pava", "b0_trend_r8": "r8"}
RULES = {"smooth": "B0_TREND_SMOOTH_RULE", "tail": "B0_TREND_TAIL_RULE", "pava": "B0_TREND_PAVA_RULE",
         "r8": "B0_TREND_R8_RULE"}


@pytest.fixture(scope="module")
def three():
    return {"ref": _stage(1600, seed=2, shift=0.0, tiny=60), "start": _stage(1100, seed=3, shift=0.8, tiny=10),
            "earlier": _stage(900, seed=5, shift=1.6, tiny=10)}


@pytest.fixture(scope="module")
def four(three):
    return {**three, "earliest": _stage(1000, seed=7, shift=2.2, tiny=10)}


def _case(three, four, case):
    return (three, _view()) if case == "three" else (four, _view(inputs=FOUR_INPUTS))


def _rec_sha(rec: dict) -> str:
    return hashlib.sha256(json.dumps(rec, sort_keys=True).encode()).hexdigest()


# ---------------------------------------------------------------------------------------------------- pins
# Taken on this branch BEFORE any edit (a0f3b5b, the synthetic stages of tests/test_b0_trend.py, BUILD_SEED_OFFSET + 1):
# sha256 of each existing heart candidate's output X (float32); every one shares b0's coordinates; and of the H13 arms'
# records (json.dumps(rec, sort_keys=True)).
PIN_X = {
    "three": {"b0": "7f9b4113111857bdc351c027f81b80bbde24153a02e07ff99c9cbd31f764e137",
              "b0_trend": "8320700817f003ae9b5daff16ac5375fcb8442313f3454e69a048e7bfe87855c",
              "b0_trend_eb": "6bde13d2e2362877b1dbf56f2fd6d340cacc8e5822490b4036dea1f8601e2f86",
              "b0_trend_ts": "a9db9a36b36e5704cdb4abdce049babfd225a4453b4b175438e90b167df5f1b7",
              "b0_trend_h": "0e652f3362444d931a72ca6a386cc67c7a313a5eb0f26a2053100ccddf138e3e",
              "b0_depth": "7f9b4113111857bdc351c027f81b80bbde24153a02e07ff99c9cbd31f764e137",
              "b0_cfm": "ddb738ea85cdcecd367054f7ecb5b794f68dbf4a4a2b6f76f645b5082d638c2f"},
    "four": {"b0": "7f9b4113111857bdc351c027f81b80bbde24153a02e07ff99c9cbd31f764e137",
             "b0_trend": "2f176a19017393e88078138a572d5baf5bf7dad4589613aafd45fec6609916dd",
             "b0_trend_eb": "ac5bb900bbfdeaa6b8fb71b293430268a0d4f0010ae62d35b54919cef27c53e9",
             "b0_trend_ts": "dc92334d5fe0ab3651eae2958b8d1a280001d537a2966a69a1a3c8fb2d290ea8",
             "b0_trend_h": "daa3daa37539e369b51eaefe3596da89452f5d4849cde6167978b9e422d16d85",
             "b0_depth": "7f9b4113111857bdc351c027f81b80bbde24153a02e07ff99c9cbd31f764e137",
             "b0_cfm": "ddb738ea85cdcecd367054f7ecb5b794f68dbf4a4a2b6f76f645b5082d638c2f"}}
PIN_COORDS = "1704b9c9251aad9a5e6df8d0fe2a5aab8a42c4647c58af60536c86915664f9ee"
PIN_REC = {"three": {"b0_trend_eb": "4f8a4ba0201672423457f2db0b9ea523ac08e6ac8572db662afd0298002de883",
                     "b0_trend_ts": "1c2d5230ee5d16c5d954cc21c1e561e98f0cb251eee03b67bbe4fc9d29c406ab",
                     "b0_trend_h": "35146133a50117e2b1abc2191356d9a6abd3061714cad97028daec43709f35d1"},
           "four": {"b0_trend_eb": "3c6103330e60da8fdca9388684919e76ab9b131f26913f2d07316b115855b988",
                    "b0_trend_ts": "7ac201379600bd595d875ee67548d7f4d2bba1d63e28d9c7044d0becd2497715",
                    "b0_trend_h": "12032f2ab4942c679159472651bfffbf8c7f500ef85412e6191eb84e71ec352f"}}


@pytest.mark.parametrize("case", ["three", "four"])
def test_existing_heart_candidates_are_unchanged(three, four, case):
    """Every existing heart candidate builds what it built before this batch (output and coordinates; the H13 arms'
    records too); b0_trend's own record is pinned in tests/test_b0_trend.py."""
    inputs, view = _case(three, four, case)
    for name, want in PIN_X[case].items():
        st = CANDIDATES[name].build(_ctx(inputs, view=view))
        assert _h(st.X) == want, name
        assert _h(st.coords) == PIN_COORDS, name
        if name in PIN_REC[case]:
            assert _rec_sha(harness.TREND_ARM_RECORDS[name]) == PIN_REC[case][name], name


def test_momentum_shift_with_apply_factors_passed_is_the_default(three):
    """factor_apply=None is momentum_shift as before; passing vec.core.apply_factors itself gives the same counts."""
    from test_b0_trend import PRE_CHANGE_SHA256, _adata, _sha

    rng = np.random.default_rng(21)
    genes = [f"g{i}" for i in range(50)]
    mu = rng.gamma(2.0, 2.0, 50)
    start = _adata(rng.poisson(mu, (300, 50)), ["A"] * 200 + ["B"] * 100, genes, "s")
    ref = _adata(rng.poisson(mu * np.linspace(0.5, 2.0, 50), (400, 50)), ["A"] * 250 + ["B"] * 120 + ["C"] * 30,
                 genes, "r")
    like = ref[np.arange(0, 400, 3)].copy()
    a, _ = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3])
    b, _ = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], factor_apply=core.apply_factors)
    assert _sha(a) == _sha(b) == PRE_CHANGE_SHA256


# ---------------------------------------------------------------------------------------------------- every arm
@pytest.mark.parametrize("name", sorted(ARMS))
def test_each_arm_is_registered_like_b0_trend(name):
    tasks = [t for ts_ in load_tasks().values() for t in ts_]
    c, base = CANDIDATES[name], CANDIDATES["b0_trend"]
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if base.applies(t)]
    assert {t.key for t in tasks if c.applies(t)} == {"T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    assert c.tasks == base.tasks and c.derived == () and c.requires is harness._heart_extrap_with_start
    assert c.params["scope"] == "chamber" and c.params["damping"] == B0_DAMPING
    rule = RULES[ARMS[name]]
    assert c.params["rule"] == f"vec.harness.{rule}" and isinstance(getattr(harness, rule), str)


@pytest.mark.parametrize("name", sorted(ARMS))
def test_each_arm_refuses_what_b0_trend_refuses(three, name):
    idea = ARMS[name]
    for kw in ({"depth": True}, {"herd": True}, {"cfm": True}):
        with pytest.raises(ValueError, match="depth, herd or cfm"):
            harness._b0_t2(_ctx(three), trend=True, trend_idea=idea, **kw)
    with pytest.raises(ValueError, match="trend_idea"):
        harness._b0_t2(_ctx(three), trend_idea=idea)
    with pytest.raises(ValueError, match="trend_idea"):
        harness._b0_t2(_ctx(three), trend=True, trend_form="eb", trend_idea=idea)
    c = CANDIDATES[name]
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), kind="interp")))
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), setting="embryo")))
    two = {k: three[k] for k in ("ref", "start")}
    with pytest.raises(ValueError, match="three"):
        c.build(_ctx(two, view=_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})))


def test_an_unknown_idea_is_refused(three):
    with pytest.raises(ValueError, match="unknown trend idea"):
        harness._b0_t2(_ctx(three), trend=True, trend_idea="blend")


@pytest.mark.parametrize("name", sorted(ARMS))
@pytest.mark.parametrize("case", ["three", "four"])
def test_each_arm_keeps_b0_trends_rows_coordinates_and_cell_number(three, four, name, case, monkeypatch):
    """The same seed draws the same rows as b0_trend, with its coordinates and cell number; cells B0 does not shift
    (Untracked, Tiny) are b0_trend's; a trend lineage moves differently; the record is b0_trend's plus an "idea"
    block naming the arm's rule."""
    inputs, view = _case(three, four, case)
    seen = []
    real = harness.draw_copy

    def spy(ctx, herd=False):
        out = real(ctx, herd=herd)
        seen.append(out[1].copy())
        return out

    monkeypatch.setattr(harness, "draw_copy", spy)
    harness.TREND_ARM_RECORDS.pop(name, None)
    a = CANDIDATES[name].build(_ctx(inputs, view=view))
    base_rec: dict = {}
    b = harness._b0_t2(_ctx(inputs, view=view), trend=True, trend_record=base_rec)
    assert len(seen) == 2 and np.array_equal(seen[0], seen[1])
    assert a.n == b.n == BOARD.max_cells
    np.testing.assert_array_equal(a.coords, b.coords)
    merged = np.asarray(inputs["ref"].obs["merged"]).astype(str)[seen[0]]
    flat = (merged == "Untracked") | (merged == "Tiny")
    assert flat.any() and np.array_equal(a.X[flat], b.X[flat])
    assert not np.array_equal(a.X[merged == "CM"], b.X[merged == "CM"])
    rec = harness.TREND_ARM_RECORDS[name]
    assert {k: v for k, v in rec.items() if k != "idea"} == base_rec
    assert rec["idea"]["rule"] == f"vec.harness.{RULES[ARMS[name]]}"


@pytest.mark.parametrize("name", sorted(ARMS))
def test_no_arm_reads_the_target(three, name):
    """Built with the target's day NaN under the read guard: nothing refused, and the output is the same."""
    blind = replace(_view(), target_day=float("nan"))
    with harness.guarded_reads() as g:
        a = CANDIDATES[name].build(_ctx(three, view=blind))
    assert not g["refused"]
    np.testing.assert_array_equal(a.X, CANDIDATES[name].build(_ctx(three)).X)


@pytest.mark.parametrize("name", sorted(ARMS))
def test_with_nothing_shifted_each_arm_is_b0_trend_bit_for_bit(three, name):
    """No lineage has min_cells cells at `start` (B0 shifts nothing): every arm has nothing to do and equals b0_trend
    bit for bit."""
    start = three["start"].copy()
    start.obs["merged"] = np.where(np.arange(start.n_obs) % 2 == 0, "Elsewhere", "Untracked")
    inputs = {**three, "start": start}
    np.testing.assert_array_equal(CANDIDATES[name].build(_ctx(inputs)).X, CANDIDATES["b0_trend"].build(_ctx(inputs)).X)


@pytest.mark.parametrize("name", ["b0_trend_pava", "b0_trend_r8"])
def test_series_arms_with_no_trend_lineage_are_b0_trend(three, name):
    """No lineage usable before `start` (every shifted lineage keeps B0's own factor): the gate and the run basis have
    no series to act on, and the output is b0_trend's bit for bit."""
    inputs = _only_start_and_ref(three)
    np.testing.assert_array_equal(CANDIDATES[name].build(_ctx(inputs)).X, CANDIDATES["b0_trend"].build(_ctx(inputs)).X)


def test_series_out_changes_nothing_and_holds_the_trend_series(three):
    ctx = _ctx(three)
    series: dict = {}
    a, rec_a = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES)
    b, rec_b = harness.b0_trend_factors(ctx, list(BOARD.genes), TYPES, series_out=series)
    assert a.keys() == b.keys() == series.keys() and all(np.array_equal(a[k], b[k]) for k in a) and rec_a == rec_b
    for lin, (ys, ds) in series.items():
        assert ds == [8.5, 8.75, 9.5]
        np.testing.assert_array_equal(harness.trend_log_factor(ys, ds, rec_a["budget_days"]), a[lin])


# ---------------------------------------------------------------------------------------------------- smooth
def _counts(seed: int = 31):
    rng = np.random.default_rng(seed)
    C = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 2.0, 40), (60, 40)))
    return C, rng.integers(0, 3, 60), rng


def test_smooth_up_is_apply_factors_when_no_factor_is_above_one():
    """Factors at or below 1 (keep_depth off, so no constant lifts them; held genes too): bit for bit apply_factors,
    whatever the expected counts."""
    C, ti, rng = _counts()
    F = np.exp(-np.abs(rng.normal(0, 0.6, (3, 40))))
    F[0, :5] = 1.0
    held = np.zeros(40, dtype=bool)
    held[[2, 5]] = True
    E = rng.gamma(1.0, 3.0, C.shape)
    a, sa, da = core.apply_factors(C, F, ti, held, keep_depth=False, seed=[5, 3])
    b, sb, db = harness.apply_factors_smooth_up(C, F, ti, held, keep_depth=False, seed=[5, 3], expected=E)
    np.testing.assert_array_equal(a.toarray(), b.toarray())
    np.testing.assert_array_equal(sa, sb)
    np.testing.assert_array_equal(da, db)


def test_smooth_up_thins_as_apply_factors_and_tops_up_from_the_expectation():
    """With keep_depth: the same keep-depth constant and, entry for entry, the same thinned counts where f < 1 (the same
    first draws); held genes untouched; where f > 1 the mean gain over many draws is (f - 1) x expected, zeros
    included (detection can rise), where the own-count top-up never lifts a zero."""
    C, ti, rng = _counts()
    F = np.exp(rng.normal(0, 0.6, (3, 40)))
    held = np.zeros(40, dtype=bool)
    held[[2, 5]] = True
    E = rng.gamma(1.0, 3.0, C.shape)
    a, sa, _ = core.apply_factors(C, F, ti, held, keep_depth=True, seed=[5, 3])
    b, sb, _ = harness.apply_factors_smooth_up(C, F, ti, held, keep_depth=True, seed=[5, 3], expected=E)
    np.testing.assert_array_equal(sa, sb)
    Fd = F[ti] * sa[:, None]
    Cd = C.toarray()
    down = (Fd < 1.0) & ~held[None, :] & (Cd > 0)
    np.testing.assert_array_equal(a.toarray()[down], b.toarray()[down])
    np.testing.assert_array_equal(b.toarray()[:, held], Cd[:, held])
    up = (Fd > 1.0) & ~held[None, :]
    gains = np.mean([harness.apply_factors_smooth_up(C, F, ti, held, keep_depth=True, seed=s, expected=E)[0].toarray()
                     - Cd for s in range(400)], axis=0)
    want = (Fd - 1.0) * E
    assert abs(gains[up].sum() / want[up].sum() - 1.0) < 0.02
    zeros_up = up & (Cd == 0) & (E > 1.0)
    assert zeros_up.any() and gains[zeros_up].mean() > 0.5
    assert (a.toarray()[Cd == 0] == 0).all()


def test_smooth_up_refuses_bad_expectations():
    C, ti, _ = _counts()
    F = np.ones((3, 40))
    with pytest.raises(ValueError, match="expected"):
        harness.apply_factors_smooth_up(C, F, ti, expected=np.ones((2, 2)))
    with pytest.raises(ValueError, match="expected"):
        harness.apply_factors_smooth_up(C, F, ti, expected=-np.ones(C.shape))


def test_smooth_expected_counts_is_the_render_routine(three):
    """The expectation is t2_build._smoothed_props over the reference pool (RENDER_K neighbours, RENDER_N_PC
    components), times each drawn cell's library size: per cell it sums to the library's expected panel count."""
    from vec import t2_build

    pool = three["ref"]
    rows = np.arange(0, pool.n_obs, 7)
    genes = list(pool.var_names)
    E, C, lib, _ = t2_build._stage_counts(pool, genes)
    want = t2_build._smoothed_props(E, C, lib, rows, t2_build.RENDER_K) * lib[rows][:, None]
    got = harness.smooth_expected_counts(pool, genes, rows)
    np.testing.assert_array_equal(got, want)
    assert t2_build.RENDER_K == CANDIDATES["b0_trend_smooth"].params["k"] == 15
    assert t2_build.RENDER_N_PC == CANDIDATES["b0_trend_smooth"].params["n_pc"] == 30


def test_smooth_arm_calls_the_smoothed_applier_once_with_b0_trends_factors(three, monkeypatch):
    """The arm hands momentum_shift exactly b0_trend's factors (same mapping, same stream); only the applier differs:
    with the applier replaced by apply_factors the arm is b0_trend bit for bit."""
    calls = []

    def as_base(C, F, tidx, held=None, keep_depth=True, seed=0, *, expected):
        calls.append(expected.shape)
        return core.apply_factors(C, F, tidx, held, keep_depth=keep_depth, seed=seed)

    monkeypatch.setattr(harness, "apply_factors_smooth_up", as_base)
    a = CANDIDATES["b0_trend_smooth"].build(_ctx(three))
    assert calls == [(BOARD.max_cells, len(BOARD.genes))]
    np.testing.assert_array_equal(a.X, CANDIDATES["b0_trend"].build(_ctx(three)).X)


# ---------------------------------------------------------------------------------------------------- tail
def test_tail_with_its_threshold_at_zero_is_b0_trend(three, monkeypatch):
    """Strength zero: with the feasibility threshold at 0 no gene is infeasible, and the arm is b0_trend bit for bit."""
    monkeypatch.setattr(harness, "TAIL_PB_MIN", 0.0)
    np.testing.assert_array_equal(CANDIDATES["b0_trend_tail"].build(_ctx(three)).X,
                                  CANDIDATES["b0_trend"].build(_ctx(three)).X)
    assert harness.TREND_ARM_RECORDS["b0_trend_tail"]["idea"]["n_infeasible"] == 0


def test_tail_mask_is_the_scorers_pseudobulk_below_a_quarter(three):
    pool = three["ref"]
    genes = list(pool.var_names)
    pb = np.asarray(pool.X.mean(axis=0)).ravel()
    np.testing.assert_array_equal(harness.tail_infeasible(pool, genes), pb < 0.25)
    assert harness.TAIL_PB_MIN == CANDIDATES["b0_trend_tail"].params["pb_min"] == 0.25
    assert harness.tail_infeasible(pool, genes).any() and not harness.tail_infeasible(pool, genes).all()


def test_tail_filter_zeroes_only_negative_factors_on_infeasible_genes():
    infeasible = np.array([True, True, False, False, True])
    logf = np.array([-0.3, 0.2, -0.4, 0.1, 0.0])
    out, rec = harness.tail_filter(infeasible)(logf, None, np.zeros(5, dtype=bool))
    np.testing.assert_array_equal(out, [0.0, 0.2, -0.4, 0.1, 0.0])
    assert rec["n_compressed"] == 1 and logf[0] == -0.3


def test_tail_arm_filters_every_shifted_lineage_and_only_those_genes(three, monkeypatch):
    """The filter sees every shifted lineage (trend and own factor alike) and nothing else changes: with the filter
    made the identity, the arm is b0_trend bit for bit."""
    seen = []
    real = harness.tail_filter

    def spy(mask):
        f = real(np.zeros_like(mask))

        def g(logf, counts, held):
            seen.append(len(logf))
            return f(logf, counts, held)

        return g

    monkeypatch.setattr(harness, "tail_filter", spy)
    a = CANDIDATES["b0_trend_tail"].build(_ctx(three))
    assert seen == [len(BOARD.genes)] * 2  # CM and Endothelial
    np.testing.assert_array_equal(a.X, CANDIDATES["b0_trend"].build(_ctx(three)).X)


# ---------------------------------------------------------------------------------------------------- pava
@pytest.mark.parametrize("y, unpooled", [
    ([1.0, 2.0, 3.0], True), ([3.0, 1.0, 2.5], True), ([3.0, 1.0, 2.0], False), ([1.0, 3.0, 2.0], False),
    ([5.0, 1.0, 2.0], False), ([1.0, 2.0, 2.0], False), ([2.0, 1.0, 1.5, 3.0], True), ([1.0, 4.0, 2.0, 3.0], False),
    ([0.0, 4.0, 1.0, 2.0, 5.0], True)])
def test_pava_last_unpooled(y, unpooled):
    """By hand: (5, 1, 2) pools 5 and 1 to 3, then 3 >= 2 pools all; (3, 1, 2): 3 and 1 pool to 2, which ties 2 (a
    tie pools); (1, 4, 2, 3): 4 and 2 pool to 3, which ties 3."""
    assert harness.pava_last_unpooled(y) is unpooled


def test_persistence_gate_reads_the_direction_of_the_last_step():
    genes = np.array([[1.0, 3.0, 2.0], [2.0, 1.0, 3.0], [3.0, 2.0, 1.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]])
    open_, s = harness.persistence_gate(list(genes.T))  # one array per stage
    np.testing.assert_array_equal(s, [-1.0, 1.0, -1.0, 0.0, 0.0])
    # gene 0 (1, 3, 2) falls last: the falling fit pools everything; gene 1 (2, 1, 3) rises last: 2 and 1 pool, 3
    # stays alone; gene 2 (3, 2, 1) falls throughout: unpooled; genes 3, 4: no last step
    np.testing.assert_array_equal(open_, [False, True, True, False, False])


def test_trend_pava_with_one_monotone_lineage_is_the_identity():
    """One lineage (no analogous series) whose series move monotonically: every gate is open and the mapping is
    b0_trend's bit for bit."""
    rng = np.random.default_rng(3)
    days = [8.3, 8.5, 8.75, 9.5]
    a, b = rng.uniform(2.0, 4.0, 60), rng.choice([-1.0, 1.0], 60) * rng.uniform(0.2, 1.5, 60)
    ys = [a + b * (d - 9.5) for d in days]
    logf = {"CM": harness.trend_log_factor(ys, days, 0.375)}
    out, rec = harness.trend_pava(logf, {"CM": (ys, days)}, 0.375)
    assert out["CM"] is not logf["CM"] and np.array_equal(out["CM"], logf["CM"])
    assert rec["lineages"][0]["n_closed"] == 0 and rec["lineages"][0]["n_combined"] == 0


def test_trend_pava_gates_and_combines_by_hand():
    """Two lineages: a gene whose last step reverses its run is closed (0); an open gene agreeing with the other
    lineage's slope gets budget x the mean of the two slopes; an open gene disagreeing keeps b0_trend's factor."""
    days = [8.5, 8.75, 9.5]
    # rows are stages, columns genes
    cm = np.array([[1.0, 1.0, 1.0], [2.0, 1.5, 3.0], [3.0, 2.5, 2.9]])  # g0, g1 rise; g2 rises then dips (pooled)
    en = np.array([[1.0, 3.0, 2.0], [2.0, 2.0, 1.0], [4.0, 1.0, 0.5]])  # g0 rises; g1, g2 fall
    budget = B0_DAMPING * 0.75
    series = {"CM": (list(cm), days), "Endothelial": (list(en), days)}
    logf = {k: harness.trend_log_factor(v[0], days, budget) for k, v in series.items()}
    out, rec = harness.trend_pava(logf, series, budget)
    b_cm, b_en = (harness.trend_log_factor(v[0], days, 1.0) for v in series.values())
    assert b_en[0] > 0 and b_en[1] < 0 and b_cm[2] > 0
    # CM: g0 open and agrees (combined), g1 open and disagrees (own), g2 closed
    assert out["CM"][0] == budget * 0.5 * (b_cm[0] + b_en[0])
    assert out["CM"][1] == logf["CM"][1]
    assert out["CM"][2] == 0.0
    # Endothelial: g0 open and agrees (combined), g1 open and disagrees (own), g2 open (falls throughout) and its
    # analogue (CM's g2 slope) rises: own
    assert out["Endothelial"][0] == budget * 0.5 * (b_en[0] + b_cm[0])
    assert out["Endothelial"][1] == logf["Endothelial"][1] and out["Endothelial"][2] == logf["Endothelial"][2]
    rows = {r["lineage"]: r for r in rec["lineages"]}
    assert (rows["CM"]["n_closed"], rows["CM"]["n_combined"], rows["CM"]["n_own"]) == (1, 1, 1)
    en_row = rows["Endothelial"]
    assert (en_row["n_closed"], en_row["n_combined"], en_row["n_own"]) == (0, 1, 2)
    assert rows["CM"]["n_analogous"] == 1


def test_trend_pava_keeps_an_open_gene_that_disagrees():
    days = [8.5, 8.75, 9.5]
    cm = [np.array([1.0]), np.array([1.5]), np.array([3.0])]  # rises throughout
    en = [np.array([3.0]), np.array([2.0]), np.array([1.0])]  # falls throughout
    budget = 0.375
    series = {"CM": (cm, days), "Endothelial": (en, days)}
    logf = {k: harness.trend_log_factor(v[0], days, budget) for k, v in series.items()}
    out, _ = harness.trend_pava(logf, series, budget)
    np.testing.assert_array_equal(out["CM"], logf["CM"])
    np.testing.assert_array_equal(out["Endothelial"], logf["Endothelial"])


def test_pava_arm_hands_momentum_shift_the_gated_mapping(three, monkeypatch):
    """Through the candidate: with the gate and combiner replaced by the identity, the arm is b0_trend bit for bit
    (the idea is the only difference)."""
    monkeypatch.setattr(harness, "trend_pava", lambda logf, series, budget: (dict(logf), {"rule": "id"}))
    np.testing.assert_array_equal(CANDIDATES["b0_trend_pava"].build(_ctx(three)).X,
                                  CANDIDATES["b0_trend"].build(_ctx(three)).X)


# ---------------------------------------------------------------------------------------------------- r8
def test_run_contrasts_are_r2s_interior_residuals():
    """For days 8.5, 8.75, 9.5 the residual at 8.75 is y - 0.75 y(8.5) - 0.25 y(9.5) (REPORT.md R2's r_875); for 8.75,
    9.5, 10.5 at 9.5 it is y - 0.571 y(8.75) - 0.429 y(10.5)."""
    y = [np.array([1.0, 2.0]), np.array([3.0, 5.0]), np.array([2.0, 1.0])]
    M = harness.run_contrasts({"A": (y, [8.5, 8.75, 9.5])})
    np.testing.assert_allclose(M, [y[1] - 0.75 * y[0] - 0.25 * y[2]], atol=1e-12)
    M2 = harness.run_contrasts({"A": (y, [8.75, 9.5, 10.5])})
    np.testing.assert_allclose(M2, [y[1] - (1.0 / 1.75) * y[0] - (0.75 / 1.75) * y[2]], atol=1e-12)
    four = harness.run_contrasts({"A": (y + [np.array([0.0, 0.0])], [8.3, 8.5, 8.75, 9.5]),
                                  "B": (y, [8.5, 8.75, 9.5])})
    assert four.shape == (3, 2)


def test_r8_on_exactly_linear_series_leaves_the_factor():
    """Exactly linear series (integer arithmetic): every contrast is 0, there is no basis, the mapping is b0_trend's."""
    days = [8.0, 9.0, 10.0, 11.0]
    ys = [np.array([1.0, 2.0, -1.0]) * d for d in days]
    logf = {"CM": harness.trend_log_factor(ys, days, 0.5)}
    out, rec = harness.trend_r8(logf, {"CM": (ys, days)}, np.array([0.1, 2.0, 0.7]))
    assert not rec["basis"] and rec["n_contrasts"] == 2
    np.testing.assert_array_equal(out["CM"], logf["CM"])


def _tied_levels_and_profile(rng, G: int = 200, n_groups: int = 4):
    """Reference levels in n_groups tied groups (genes shuffled), and a unit run profile p with zero mean inside each
    group: p has no part monotone in the level (its isotonic fit, a function of the level, is 0)."""
    level = rng.permutation(np.repeat(np.arange(1.0, n_groups + 1.0) * 0.4, G // n_groups))
    p = rng.normal(size=G)
    for v in np.unique(level):
        p[level == v] -= p[level == v].mean()
    return level, p / np.linalg.norm(p)


def test_r8_removes_a_planted_run_profile():
    """Two lineages, linear trends plus one shared run profile p with a different weight at each stage, p with no
    part monotone in the reference level: the rank-1 basis is p (to sign), and each factor loses exactly its
    component along p, as before the level step."""
    rng = np.random.default_rng(41)
    G = 200
    level, p = _tied_levels_and_profile(rng, G)
    days = [8.3, 8.5, 8.75, 9.5]
    load = [0.0, 2.0, -1.5, 1.0]
    series = {}
    for lin in ("CM", "Endothelial"):
        a, b = rng.normal(2.0, 1.0, G), rng.normal(0.0, 1.0, G)
        series[lin] = ([a + b * d + w * p for d, w in zip(days, load, strict=True)], days)
    logf = {k: harness.trend_log_factor(v[0], days, 0.375) for k, v in series.items()}
    out, rec = harness.trend_r8(logf, series, level)
    assert rec["basis"] and rec["n_contrasts"] == 4 and rec["basis_energy_share"] > 0.999
    assert all(r["removed_energy_share"] < 1e-20 for r in rec["level_monotone"]["contrasts"])
    for lin in series:
        np.testing.assert_allclose(out[lin], logf[lin] - (logf[lin] @ p) * p, atol=1e-9)
        assert abs(out[lin] @ p) < 1e-9
    assert harness.R8_RANK == CANDIDATES["b0_trend_r8"].params["rank"] == 1


def test_r8_contrasts_come_from_the_tasks_own_stages(three, monkeypatch):
    """The arm's contrasts are those of b0_trend_factors' trend series on the task's input stages (one interior stage
    per lineage here: two contrasts)."""
    seen = []
    real = harness.run_contrasts

    def spy(series):
        seen.append({k: v[1] for k, v in series.items()})
        return real(series)

    monkeypatch.setattr(harness, "run_contrasts", spy)
    CANDIDATES["b0_trend_r8"].build(_ctx(three))
    assert seen == [{"CM": [8.5, 8.75, 9.5], "Endothelial": [8.5, 8.75, 9.5]}]
    assert harness.TREND_ARM_RECORDS["b0_trend_r8"]["idea"]["n_contrasts"] == 2


# ---------------------------------------------------------------------------------------------------- r8 level step
# Fix of campaign/heart-b's b0_trend_r8 (REPORT.md section 6, common rules: "Before taking residuals (R2, R3, R8),
# remove the part monotone in reference level"). Taken at e4b66d8 BEFORE the fix, on the synthetic stages above:
# sha256 of the other three batch-b arms' output X and records (json.dumps(rec, sort_keys=True)), and of b0_trend_r8's
# own output X as e4b66d8 built it.
PIN_B_X = {"three": {"b0_trend_smooth": "1072a14e87248b1012dac2daec85b786dbe0cb2bf88661e2771ac178f5eab51b",
                     "b0_trend_tail": "6e88977fd2d3fb7b54c65d6fffcb9ab97cce7bc83648e9c7e094744e7a091e8f",
                     "b0_trend_pava": "c473f6dc29a93c719796e4128b7149885a7d64f1bcdd35d6a4c370913c8ad28a"},
           "four": {"b0_trend_smooth": "62d3e8d50e4380054ec4734661c2abf5ff0bc08e350804fbf9a250ba7e650b36",
                    "b0_trend_tail": "1f6af82912d039db1811531532e52197c58da18fe388c0e0f1cdf44a15eb693b",
                    "b0_trend_pava": "8ca827e4728cd259a37bf8a8d02d8ba83243197035695168d5f0d4cd6754bf18"}}
PIN_B_REC = {"three": {"b0_trend_smooth": "ff741cd0024b24298253dedb6cba24c91bb7fd4c8fcf3cb150f4e620a372db4d",
                       "b0_trend_tail": "4b9a036dcefb88d4989bad745e68a3e8c088a11102a4753253305e9066e05073",
                       "b0_trend_pava": "2dfd87ee8e3a82233bf37b016eed21101176c6e6dbc25bb7af3f43176e649093"},
             "four": {"b0_trend_smooth": "49d7346fe2980095c10cd0254bc9452865aef0c195fd52cdef019bfd0652a7dc",
                      "b0_trend_tail": "052c32aee1ef129b9b88a21daed30d300166a1ee7b79adebbcbdb585501d8cdf",
                      "b0_trend_pava": "4d6fdeaa3ba5993c1429f30471aac45e2e94ffae0a87227e1f7b742d6fe083b8"}}
PIN_R8_E4B66D8_X = {"three": "a7b8215f6cc8f045b34e8e015243e12a3171d36186c3bb60a6a71360f35bbc4b",
                    "four": "98e919eccdfaf03072c8d9eef4b0a38866162d14ff6d7e8c53b5accbd4442044"}


@pytest.mark.parametrize("case", ["three", "four"])
def test_the_other_batch_b_arms_are_unchanged_by_the_r8_fix(three, four, case):
    """b0_trend_smooth, b0_trend_tail and b0_trend_pava build what they built at e4b66d8 (output, coordinates and
    record); the older candidates are pinned in test_existing_heart_candidates_are_unchanged."""
    inputs, view = _case(three, four, case)
    for name, want in PIN_B_X[case].items():
        st = CANDIDATES[name].build(_ctx(inputs, view=view))
        assert _h(st.X) == want, name
        assert _h(st.coords) == PIN_COORDS, name
        assert _rec_sha(harness.TREND_ARM_RECORDS[name]) == PIN_B_REC[case][name], name


@pytest.mark.parametrize("case", ["three", "four"])
def test_r8_with_the_level_fit_made_zero_is_its_e4b66d8_output(three, four, case, monkeypatch):
    """The level step is the only change: with the level-monotone fit replaced by zeros, b0_trend_r8 builds its
    e4b66d8 output bit for bit."""
    inputs, view = _case(three, four, case)
    monkeypatch.setattr(harness, "level_monotone_fit", lambda v, level: (np.zeros_like(v), "up"))
    assert _h(CANDIDATES["b0_trend_r8"].build(_ctx(inputs, view=view)).X) == PIN_R8_E4B66D8_X[case]


def _r8_e4b66d8(logf_by_type, series, rank=1):
    """b0_trend_r8's mapping as e4b66d8 built it (no level step): the uncentred rank-1 SVD basis of the raw
    contrasts, projected off each trend lineage's factor."""
    M = harness.run_contrasts(series)
    out = dict(logf_by_type)
    if M.shape[0] == 0 or not np.any(M):
        return out
    _, sv, Vt = np.linalg.svd(M, full_matrices=False)
    V = Vt[:rank][sv[:rank] > 0]
    for lin in series:
        logf = np.asarray(logf_by_type[lin], dtype=np.float64)
        out[lin] = logf - V.T @ (V @ logf)
    return out


def test_r8_with_contrasts_that_have_no_monotone_part_is_its_e4b66d8_mapping():
    """Contrasts with no part monotone in the reference level (zero mean within each tied level group): the
    isotonic fits are 0, and the mapping equals e4b66d8's to 1e-9 (two profiles, so the rank-1 basis is a mix)."""
    rng = np.random.default_rng(43)
    G = 240
    level, p = _tied_levels_and_profile(rng, G)
    q = rng.normal(size=G)
    for v in np.unique(level):
        q[level == v] -= q[level == v].mean()
    q /= np.linalg.norm(q)
    days = [8.3, 8.5, 8.75, 9.5]
    series = {}
    for lin, (wp, wq) in {"CM": ([0.0, 2.0, -1.5, 1.0], [0.0, 0.5, 0.3, 0.0]),
                          "Endothelial": ([0.0, -1.0, 0.8, 0.0], [1.0, 0.0, 0.7, -0.4])}.items():
        a, b = rng.normal(2.0, 1.0, G), rng.normal(0.0, 1.0, G)
        series[lin] = ([a + b * d + x * p + y * q for d, x, y in zip(days, wp, wq, strict=True)], days)
    for row in harness.run_contrasts(series):
        fit, _ = harness.level_monotone_fit(row, level)
        assert np.abs(fit).max() < 1e-12
    logf = {k: harness.trend_log_factor(v[0], days, 0.375) for k, v in series.items()}
    out, rec = harness.trend_r8(logf, series, level)
    want = _r8_e4b66d8(logf, series)
    assert rec["basis"]
    for lin in series:
        assert not np.allclose(want[lin], logf[lin])
        np.testing.assert_allclose(out[lin], want[lin], atol=1e-9, rtol=0)


@pytest.mark.parametrize("tied", [False, True])
def test_r8_basis_ignores_a_planted_level_monotone_component(tied):
    """Contrasts = a non-monotone run profile p plus a large component m that rises with the reference level, with
    the same sign in every lineage's contrasts. Without the level step (e4b66d8) the basis leans on m; with it the
    basis is orthogonal to m within tolerance (exactly, when p has zero mean within tied level groups) and p is still
    what each factor loses."""
    rng = np.random.default_rng(47)
    G = 300
    if tied:
        level, p = _tied_levels_and_profile(rng, G, n_groups=6)
    else:
        level = rng.gamma(1.0, 1.0, G)
        p = rng.normal(size=G)
        p /= np.linalg.norm(p)
    m = np.log1p(level) - np.log1p(level).mean()
    m /= np.linalg.norm(m)
    days = [8.3, 8.5, 8.75, 9.5]
    load = [0.0, 2.0, -1.5, 1.0]
    series = {}
    for amp, lin in ((6.0, "CM"), (7.0, "Endothelial")):
        a, b = rng.normal(2.0, 1.0, G), rng.normal(0.0, 1.0, G)
        series[lin] = ([a + b * d + w * p + (amp * m if k in (1, 2) else 0.0)
                        for k, (d, w) in enumerate(zip(days, load, strict=True))], days)
    M = harness.run_contrasts(series)
    assert ((M @ m) > 0).all()  # the planted level component has one sign in every contrast
    assert abs(np.linalg.svd(M, full_matrices=False)[2][0] @ m) > 0.5  # e4b66d8's basis leans on it
    logf = {k: harness.trend_log_factor(v[0], days, 0.375) for k, v in series.items()}
    out, rec = harness.trend_r8(logf, series, level)
    assert rec["basis"]
    assert {r["direction"] for r in rec["level_monotone"]["contrasts"]} == {"up"}
    for lin in series:
        d = logf[lin] - out[lin]
        basis = d / np.linalg.norm(d)  # rank 1: what a factor loses lies along the basis
        assert abs(basis @ m) < (1e-9 if tied else 0.05)
        assert abs(basis @ p) > (1 - 1e-9 if tied else 0.95)
        if tied:
            np.testing.assert_allclose(out[lin], logf[lin] - (logf[lin] @ p) * p, atol=1e-9)


@pytest.mark.parametrize("ties", [False, True])
def test_level_monotone_residual_of_a_monotone_vector_is_zero(ties):
    """A vector monotone in the level (rising or falling, genes in shuffled order, constant on tied levels) is its
    own fit: the isotonic residual is 0; the direction names the fit used."""
    rng = np.random.default_rng(53)
    level = rng.permutation(np.repeat(np.linspace(0.0, 3.0, 25), 4)) if ties else rng.gamma(1.0, 1.0, 100)
    for v, direction in ((np.sqrt(level) - 0.7, "up"), (-2.0 * level ** 2 + 1.0, "down")):
        fit, got = harness.level_monotone_fit(v, level)
        np.testing.assert_allclose(v - fit, 0.0, atol=1e-12)
        assert got == direction
    constant = np.full(100, 0.4)
    fit, got = harness.level_monotone_fit(constant, level)
    np.testing.assert_allclose(fit, constant, atol=1e-15)
    assert got == "up"


def test_level_monotone_fit_by_hand():
    """Levels (0.1, 0.5, 0.2, 0.9) order the genes 0, 2, 1, 3, so v (1, 2, 3, 4) reads (1, 3, 2, 4) in level order:
    the rising fit pools 3 and 2 to 2.5 (residual sum of squares 0.5), the falling one is the mean 2.5 throughout
    (5.0), and the rising one is used. Tied levels share one value."""
    fit, d = harness.level_monotone_fit(np.array([1.0, 2.0, 3.0, 4.0]), np.array([0.1, 0.5, 0.2, 0.9]))
    np.testing.assert_allclose(fit, [1.0, 2.5, 2.5, 4.0])
    assert d == "up"
    fit, d = harness.level_monotone_fit(np.array([4.0, 0.0, 1.0]), np.array([0.3, 0.3, 0.1]))
    np.testing.assert_allclose(fit, [2.0, 2.0, 1.0])
    assert d == "up"
    fit, d = harness.level_monotone_fit(np.array([3.0, 1.0, 2.0]), np.array([0.1, 0.2, 0.3]))
    np.testing.assert_allclose(fit, [3.0, 1.5, 1.5])
    assert d == "down"
    with pytest.raises(ValueError, match="1-d"):
        harness.level_monotone_fit(np.zeros(3), np.zeros(4))


def test_pava_blocks_with_unit_counts_are_the_gates_routine():
    """pava_last_unpooled reads pava_blocks with unit counts; counts weight a position as that many equal values."""
    assert [b[:2] for b in harness.pava_blocks([3.0, 1.0, 2.5])] == [[4.0, 2.0], [2.5, 1.0]]
    assert [b[:2] for b in harness.pava_blocks([6.0, 1.0], counts=[3, 1])] == [[7.0, 4.0]]  # means 2 > 1: pooled
    assert [b[:2] for b in harness.pava_blocks([2.0, 1.0], counts=[3, 1])] == [[2.0, 3.0], [1.0, 1.0]]  # 2/3 < 1
    assert [b[2] for b in harness.pava_blocks([1.0, 4.0, 2.0, 3.0])] == [0.0, 1.0]


def test_r8_arm_orders_genes_by_the_tail_arms_reference_pb(three, monkeypatch):
    """The arm's level is the reference pool's pb (mean log1p CP10k over its cells), the pb b0_trend_tail holds against
    0.25, in the factors' gene order; the record names the level step for each contrast."""
    seen = []
    real = harness.trend_r8

    def spy(logf, series, ref_level, rank=None):
        seen.append(np.array(ref_level))
        return real(logf, series, ref_level, rank)

    monkeypatch.setattr(harness, "trend_r8", spy)
    CANDIDATES["b0_trend_r8"].build(_ctx(three))
    pool = three["ref"]
    genes = list(BOARD.genes)
    pb = np.asarray(pool[:, genes].X.mean(axis=0)).ravel().astype(np.float64)
    assert len(seen) == 1
    np.testing.assert_array_equal(seen[0], pb)
    np.testing.assert_array_equal(harness.reference_pb(pool, genes), pb)
    np.testing.assert_array_equal(harness.tail_infeasible(pool, genes), harness.reference_pb(pool, genes) < 0.25)
    rows = harness.TREND_ARM_RECORDS["b0_trend_r8"]["idea"]["level_monotone"]["contrasts"]
    assert len(rows) == 2
    assert all(r["direction"] in ("up", "down") and 0.0 <= r["removed_energy_share"] <= 1.0 for r in rows)
    assert "monotone in reference level" in harness.B0_TREND_R8_RULE
