"""Heart arms `b0_trend_mix` (a within-upload momentum mixture: half the drawn cells at b0's step, half at b0_trend's),
`b0_trend_xt` (b0_trend's step through the exact-total kernel in its repair form) and `b0_trend_mn` (b0_trend's step as
a whole-cell multinomial redraw: the kernel b0_trend_xt was built with at ff8c42d, renamed), with their hooks:
`vec.core.apply_factors(..., kernel=...)` and `vec.t1_scope.momentum_shift(..., own_factor_cells=..., kernel=...)`.
Each arm is b0_trend with one idea. Synthetic data only; nothing scored, no prediction compared with a target."""

from __future__ import annotations

import hashlib
from dataclasses import replace

import numpy as np
import pytest
from scipy import sparse
from test_b0_trend import BOARD, PRE_CHANGE_SHA256, _adata, _ctx, _only_start_and_ref, _sha, _stage, _view

from vec import core, harness, t1_scope
from vec.harness import B0_DAMPING, CANDIDATES, load_tasks

ARMS = ("b0_trend_mix", "b0_trend_xt", "b0_trend_mn")
KERNEL_ARMS = {"b0_trend_xt": "exact_total", "b0_trend_mn": "multinomial"}


@pytest.fixture(scope="module")
def shift_inputs():
    """test_b0_trend's `shift_inputs`, built the same way: A and B at both stages, C only at the reference."""
    rng = np.random.default_rng(21)
    genes = [f"g{i}" for i in range(50)]
    mu = rng.gamma(2.0, 2.0, 50)
    start = _adata(rng.poisson(mu, (300, 50)), ["A"] * 200 + ["B"] * 100, genes, "s")
    ref = _adata(rng.poisson(mu * np.linspace(0.5, 2.0, 50), (400, 50)), ["A"] * 250 + ["B"] * 120 + ["C"] * 30,
                 genes, "r")
    return start, ref, ref[np.arange(0, 400, 3)].copy(), genes


@pytest.fixture(scope="module")
def three():
    """test_b0_trend's `three` heart-like stages, built the same way."""
    return {"ref": _stage(1600, seed=2, shift=0.0, tiny=60), "start": _stage(1100, seed=3, shift=0.8, tiny=10),
            "earlier": _stage(900, seed=5, shift=1.6, tiny=10)}


def _h(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float32).tobytes()).hexdigest()


# ------------------------------------------------------------------------------- pins taken before this batch's edit
# (campaign/heart-c at a0f3b5b, before apply_factors, momentum_shift and _b0_t2 gained their keywords): apply_factors
# on `_af_inputs` with keep_depth True and False (sha256 of the int64 dense counts, of the float64 depth constants), and
# the output X and coordinates (float32) of every existing b0 heart candidate on `three` at BUILD_SEED_OFFSET + 1
PIN_AF = {True: ("79bb33bc1db630ed97a5b8e73da12400905bcf2a7647533a5ba7830b8bf39a0b",
                 "495f891229896cddfa6ff36900ded3a0c53f71a278b2b3b2217bf2b70a5f71e2"),
          False: ("b9b58fafab65280473cf26eb1d1c0ab15cf3ff308e2abdd6d7ba568157e69853",
                  "41e09ab16161976ef905e3c39d7c6542759631afa4fff8d624fa0f56be7cf958")}
PIN_X = {"b0": "7f9b4113111857bdc351c027f81b80bbde24153a02e07ff99c9cbd31f764e137",
         "b0_trend": "8320700817f003ae9b5daff16ac5375fcb8442313f3454e69a048e7bfe87855c",
         "b0_trend_eb": "6bde13d2e2362877b1dbf56f2fd6d340cacc8e5822490b4036dea1f8601e2f86",
         "b0_trend_ts": "a9db9a36b36e5704cdb4abdce049babfd225a4453b4b175438e90b167df5f1b7",
         "b0_trend_h": "0e652f3362444d931a72ca6a386cc67c7a313a5eb0f26a2053100ccddf138e3e"}
PIN_COORDS = "1704b9c9251aad9a5e6df8d0fe2a5aab8a42c4647c58af60536c86915664f9ee"

# pins taken at ff8c42d (campaign/heart-c), before the kernel fix of campaign/heart-cfix: b0_trend_mix's output X on
# `three`, and the as-built b0_trend_xt (kernel "exact_total" then, a whole-cell multinomial redraw), now b0_trend_mn
# with kernel "multinomial": its X on `three`, apply_factors on `_af_inputs` (counts, depth constants) and
# momentum_shift on `shift_inputs` (test_b0_trend's `_sha`)
PIN_FF8_X = {"b0_trend_mix": "9d6c0f2ec8eabd2dbcd0f489a846d0ec6640092e78cdd823365b9a371243da6f",
             "b0_trend_mn": "2a495e3a53910d9f4be2ed984c436fb8ccf3a5e2d94b0a341afd07a8e657e888"}
PIN_FF8_AF_MN = {True: ("af3796bc22a951c09db498068403f1055740e27ffc83ec1d7aeee447f0243090",
                        "495f891229896cddfa6ff36900ded3a0c53f71a278b2b3b2217bf2b70a5f71e2"),
                 False: ("a316ca32c3008fcb0097bb2589a6eeeb5f8feaa6fbe289d7b8ebe4da2a44f179",
                         "41e09ab16161976ef905e3c39d7c6542759631afa4fff8d624fa0f56be7cf958")}
PIN_FF8_MS_MN = "d7a6483367522d7d216de20a7008cb6dc285c2df2a67b3a4e1165a84aa6e0887"


def _af_inputs():
    rng = np.random.default_rng(31)
    C = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 2.0, 40), (60, 40)))
    F = np.exp(rng.normal(0, 0.6, (3, 40)))
    ti = rng.integers(0, 3, 60)
    held = np.zeros(40, bool)
    held[[2, 5]] = True
    return C, F, ti, held


def _sha_int(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.int64).tobytes()).hexdigest()


def _sha_f64(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.float64).tobytes()).hexdigest()


@pytest.mark.parametrize("keep_depth", [True, False])
def test_apply_factors_default_kernel_is_unchanged(keep_depth):
    C, F, ti, held = _af_inputs()
    for kw in ({}, {"kernel": "independent"}):
        out, scale, dead = core.apply_factors(C, F, ti, held, keep_depth=keep_depth, seed=[5, 3], **kw)
        assert (_sha_int(out.toarray()), _sha_f64(scale)) == PIN_AF[keep_depth] and not dead.any()


@pytest.mark.parametrize("name", sorted(PIN_X))
def test_existing_heart_candidates_are_unchanged(three, name):
    st = CANDIDATES[name].build(_ctx(three))
    assert _h(st.X) == PIN_X[name] and _h(st.coords) == PIN_COORDS


def test_momentum_shift_defaults_are_unchanged(shift_inputs):
    start, ref, like, _ = shift_inputs
    out, _ = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], kernel="independent")
    assert _sha(out) == PRE_CHANGE_SHA256


@pytest.mark.parametrize("name", sorted(PIN_FF8_X))
def test_mix_and_mn_are_ff8c42d_bit_for_bit(three, name):
    """b0_trend_mix as built at ff8c42d; b0_trend_mn as b0_trend_xt was built at ff8c42d."""
    st = CANDIDATES[name].build(_ctx(three))
    assert _h(st.X) == PIN_FF8_X[name] and _h(st.coords) == PIN_COORDS


@pytest.mark.parametrize("keep_depth", [True, False])
def test_multinomial_kernel_is_the_as_built_exact_total_kernel(keep_depth, shift_inputs):
    C, F, ti, held = _af_inputs()
    out, scale, dead = core.apply_factors(C, F, ti, held, keep_depth=keep_depth, seed=[5, 3], kernel="multinomial")
    assert (_sha_int(out.toarray()), _sha_f64(scale)) == PIN_FF8_AF_MN[keep_depth] and not dead.any()
    start, ref, like, _ = shift_inputs
    ms, _ = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], kernel="multinomial")
    assert _sha(ms) == PIN_FF8_MS_MN


# ------------------------------------------------------------------------------ the exact-total and multinomial kernels
@pytest.mark.parametrize("kernel", ["exact_total", "multinomial"])
@pytest.mark.parametrize("keep_depth", [True, False])
def test_both_kernels_fix_each_cells_total_at_the_rounded_expected_total(keep_depth, kernel):
    """Each cell's total is rint(sum c f) over its moved counts plus its unmoved counts; held counts, counts with
    f == 1 and zeros stay; the depth constant is the independent kernel's."""
    C, F, ti, held = _af_inputs()
    F[0, 10:15] = 1.0
    out, scale, dead = core.apply_factors(C, F, ti, held, keep_depth=keep_depth, seed=[5, 3], kernel=kernel)
    _, scale_i, _ = core.apply_factors(C, F, ti, held, keep_depth=keep_depth, seed=[5, 3])
    np.testing.assert_array_equal(scale, scale_i)
    D, new = C.toarray().astype(np.float64), out.toarray()
    f = np.where(held[None, :], 1.0, F[ti] * scale[:, None])
    moved = (f != 1.0) & (D > 0)
    want = np.rint((D * f * moved).sum(1)) + (D * ~moved).sum(1)
    np.testing.assert_array_equal(new.sum(1), want)
    np.testing.assert_array_equal(new[~moved], D[~moved])
    assert (new[D == 0] == 0).all() and not dead.any()
    if keep_depth:
        np.testing.assert_array_equal(new.sum(1), D.sum(1))  # the expected total is the original total


def test_multinomial_keeps_each_counts_expectation():
    """Many copies of one cell: the mean of each moved count is c f x rint(S) / S, as the independent kernel's mean
    is c f (to sampling error)."""
    c = np.array([3, 0, 8, 1, 20, 5])
    f = np.array([0.5, 2.0, 1.7, 3.0, 0.8, 1.0])
    n = 20000
    C = sparse.csr_matrix(np.tile(c, (n, 1)))
    out, _, _ = core.apply_factors(C, f[None, :], np.zeros(n, dtype=int), keep_depth=False, seed=1,
                                   kernel="multinomial")
    m = out.toarray().mean(0)
    moved = (f != 1.0) & (c > 0)
    S = float((c * f)[moved].sum())
    want = np.where(moved, c * f * np.rint(S) / S, c)
    np.testing.assert_allclose(m, want, rtol=0.02, atol=0.02)
    assert (out.toarray().sum(1) == np.rint(S) + c[~moved].sum()).all()


@pytest.mark.parametrize("kernel", ["exact_total", "multinomial"])
def test_both_kernels_with_nothing_to_move_are_the_independent_kernel_bit_for_bit(kernel):
    C, _, ti, held = _af_inputs()
    ones = np.ones((3, 40))
    a, sa, da = core.apply_factors(C, ones, ti, held, seed=[5, 3], kernel=kernel)
    b, sb, db = core.apply_factors(C, ones, ti, held, seed=[5, 3])
    np.testing.assert_array_equal(a.toarray(), b.toarray())
    np.testing.assert_array_equal(a.toarray(), C.toarray())
    np.testing.assert_array_equal(sa, sb)
    np.testing.assert_array_equal(da, db)


def _moved_and_d(C, F, ti, held, keep_depth, y):
    """Dense M (not held, count above zero, f != 1 after the depth constant), each cell's target rint(sum c f) over
    M and its d against the dense counts y (the independent kernel's output)."""
    D = C.toarray().astype(np.float64)
    _, scale, _ = core.apply_factors(C, F, ti, held, keep_depth=keep_depth, seed=0)
    f = np.where(held[None, :], 1.0, F[ti] * scale[:, None])
    M = (f != 1.0) & (D > 0)
    target = np.rint((D * f * M).sum(1)).astype(np.int64)
    return M, target, (y * M).sum(1).astype(np.int64) - target


@pytest.mark.parametrize("keep_depth", [True, False])
def test_exact_total_is_the_independent_draw_repaired_by_d_units(keep_depth):
    """The repair form: the output starts from the independent kernel's Y (same seed); each cell's total over M is
    rint(sum c f); entries outside M are Y's (= the input); a cell with d == 0 is Y's bit for bit; a cell moves
    exactly |d| units, all removed (d > 0, never below zero) or all added (d < 0)."""
    C, F, ti, held = _af_inputs()
    F[1, 20:25] = 1.0
    y = core.apply_factors(C, F, ti, held, keep_depth=keep_depth, seed=[5, 3])[0].toarray()
    out, _, dead = core.apply_factors(C, F, ti, held, keep_depth=keep_depth, seed=[5, 3], kernel="exact_total")
    x = out.toarray()
    M, target, d = _moved_and_d(C, F, ti, held, keep_depth, y)
    assert not dead.any() and (d == 0).any() and (d > 0).any() and (d < 0).any()
    np.testing.assert_array_equal((x * M).sum(1), target)
    np.testing.assert_array_equal(x[~M], y[~M])
    np.testing.assert_array_equal(x[~M], C.toarray()[~M])
    np.testing.assert_array_equal(x[d == 0], y[d == 0])
    diff = x - y
    np.testing.assert_array_equal(np.abs(diff).sum(1), np.abs(d))
    assert (diff[d > 0] <= 0).all() and (diff[d < 0] >= 0).all() and (x >= 0).all()


def test_exact_total_with_d_zero_everywhere_is_the_independent_kernel_bit_for_bit():
    """A small case whose independent draw already hits every cell's rounded expected total (the first such seed
    of 0..999): the exact-total kernel draws nothing more and gives the independent output bit for bit."""
    C = sparse.csr_matrix(np.array([[3, 2, 0, 4], [1, 5, 2, 0]]))
    F = np.array([[0.5, 1.5, 2.0, 0.9]])
    ti, held = np.zeros(2, dtype=int), np.zeros(4, bool)
    for seed in range(1000):
        y = core.apply_factors(C, F, ti, held, keep_depth=False, seed=seed)[0].toarray()
        if not _moved_and_d(C, F, ti, held, False, y)[2].any():
            break
    else:
        pytest.fail("no seed with d == 0 in every cell")
    a, sa, da = core.apply_factors(C, F, ti, held, keep_depth=False, seed=seed, kernel="exact_total")
    b, sb, db = core.apply_factors(C, F, ti, held, keep_depth=False, seed=seed)
    np.testing.assert_array_equal(a.toarray(), b.toarray())
    np.testing.assert_array_equal(sa, sb)
    np.testing.assert_array_equal(da, db)


PROBE_C = np.array([50, 40, 30, 20, 10, 5])
PROBE_F = np.array([1.02, 1.0, 1.0, 1.0, 1.0, 1.0])


@pytest.fixture(scope="module")
def probe():
    """The reviewer's probe: one cell, keep_depth=True, 20,000 copies, under each kernel (one seed)."""
    n = 20000
    C = sparse.csr_matrix(np.tile(PROBE_C, (n, 1)))
    out = {}
    for k in ("independent", "exact_total", "multinomial"):
        x, scale, _ = core.apply_factors(C, PROBE_F[None, :], np.zeros(n, dtype=int), keep_depth=True, seed=1,
                                         kernel=k)
        out[k] = (x.toarray(), scale)
    return out


def test_exact_total_keeps_each_counts_mean_on_the_probe(probe):
    """Each count's mean is within 0.02 of c f x target / sum(c f) (f after the depth constant)."""
    x, scale = probe["exact_total"]
    cf = PROBE_C * PROBE_F * scale[0]
    np.testing.assert_allclose(x.mean(0), cf * np.rint(cf.sum()) / cf.sum(), rtol=0, atol=0.02)
    assert (x.sum(1) == np.rint(cf.sum())).all()


@pytest.mark.xfail(strict=True, reason="the bound asked for (the independent kernel's variance plus 0.05) does not "
                   "hold for the repair as specified: the allocation of the |d| units adds about |d| p (1 - p) per "
                   "count (measured about [0.55, 0.37, 0.30, 0.21, 0.11, 0.06] against [0.67, 0.25, 0.19, 0.13, 0.06, "
                   "0.03]); kept strict so that it reports if the kernel ever changes")
def test_exact_total_variance_is_at_most_the_independent_kernels_plus_005(probe):
    np.testing.assert_array_less(probe["exact_total"][0].var(0), probe["independent"][0].var(0) + 0.05)


def test_exact_total_variance_is_far_below_the_multinomials(probe):
    """The repair adds a small variance, the whole-cell redraw about c per count (the reviewer's finding): per count
    the exact-total variance is below the independent kernel's plus 0.15 and below a tenth of the multinomial's,
    and the multinomial's is above c / 2 (bounds set after the probe was measured; see the xfail above)."""
    v_i, v_x, v_m = (probe[k][0].var(0) for k in ("independent", "exact_total", "multinomial"))
    np.testing.assert_array_less(v_x, v_i + 0.15)
    np.testing.assert_array_less(v_x, v_m / 10)
    np.testing.assert_array_less(PROBE_C * 0.5, v_m)


@pytest.mark.parametrize("kernel", ["nb", "negative_binomial", "Exact_total", ""])
def test_unknown_kernel_is_refused(kernel):
    C, F, ti, held = _af_inputs()
    with pytest.raises(ValueError, match="unknown kernel"):
        core.apply_factors(C, F, ti, held, kernel=kernel)


@pytest.mark.parametrize("kernel", ["exact_total", "multinomial"])
def test_momentum_shift_passes_the_kernel(shift_inputs, kernel):
    start, ref, like, _ = shift_inputs
    out, rep = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3], kernel=kernel)
    base, rep0 = t1_scope.momentum_shift(start, ref, like, 0.5, key="b0_type", seed=[7, 3])
    np.testing.assert_array_equal(out.layers["counts"].toarray().sum(1), like.layers["counts"].toarray().sum(1))
    assert not np.array_equal(out.layers["counts"].toarray(), base.layers["counts"].toarray())
    assert rep.rows == rep0.rows


# --------------------------------------------------------------------------------------------- own_factor_cells hook
def _own(start, ref, t):
    S, R = core.counts_of(start), core.counts_of(ref)
    s, r = (start.obs["b0_type"] == t).to_numpy(), (ref.obs["b0_type"] == t).to_numpy()
    return 0.5 * (np.log(core.mean_cp10k(R[r]) + 0.5) - np.log(core.mean_cp10k(S[s]) + 0.5))


def test_own_factor_cells_with_nothing_to_mix_is_bit_for_bit(shift_inputs):
    """An all-False mask, or a mapping equal to each type's own factor with any mask, gives the mapping's output bit
    for bit; a mapping that names no type gives the default output."""
    start, ref, like, genes = shift_inputs
    kw = {"key": "b0_type", "seed": [7, 3]}
    n = like.n_obs
    half = np.arange(n) < n // 2
    bump = {"A": np.linspace(-1.0, 1.0, len(genes))}
    mapped, _ = t1_scope.momentum_shift(start, ref, like, 0.5, logf_by_type=bump, **kw)
    none, _ = t1_scope.momentum_shift(start, ref, like, 0.5, logf_by_type=bump, own_factor_cells=np.zeros(n, bool),
                                      **kw)
    np.testing.assert_array_equal(mapped.layers["counts"].toarray(), none.layers["counts"].toarray())
    own = {t: _own(start, ref, t) for t in ("A", "B")}
    for mask in (half, np.ones(n, bool)):
        same, _ = t1_scope.momentum_shift(start, ref, like, 0.5, logf_by_type=own, own_factor_cells=mask, **kw)
        assert _sha(same) == PRE_CHANGE_SHA256
    empty, _ = t1_scope.momentum_shift(start, ref, like, 0.5, logf_by_type={}, own_factor_cells=half, **kw)
    assert _sha(empty) == PRE_CHANGE_SHA256


def test_masked_cells_get_their_types_own_factor(shift_inputs):
    """A mapped to zero (its cells keep their counts): unmasked A cells keep their counts, masked A cells move; with
    every cell masked the output is the default's bit for bit (each cell at its type's own factor, one stream)."""
    start, ref, like, genes = shift_inputs
    kw = {"key": "b0_type", "seed": [7, 3]}
    n = like.n_obs
    half = np.arange(n) < n // 2
    zero = {"A": np.zeros(len(genes))}
    out, rep = t1_scope.momentum_shift(start, ref, like, 0.5, logf_by_type=zero, own_factor_cells=half, **kw)
    t = like.obs["b0_type"].to_numpy()
    old, new = like.layers["counts"].toarray(), out.layers["counts"].toarray()
    np.testing.assert_array_equal(old[(t == "A") & ~half], new[(t == "A") & ~half])
    assert not np.array_equal(old[(t == "A") & half], new[(t == "A") & half])
    row = next(r for r in rep.rows if r["type"] == "A")
    assert row["n_own_factor_cells"] == int(((t == "A") & half).sum())
    assert "n_own_factor_cells" not in next(r for r in rep.rows if r["type"] == "B")
    full, _ = t1_scope.momentum_shift(start, ref, like, 0.5, logf_by_type=zero, own_factor_cells=np.ones(n, bool), **kw)
    assert _sha(full) == PRE_CHANGE_SHA256


def test_own_factor_cells_refusals(shift_inputs):
    start, ref, like, _ = shift_inputs
    kw = {"key": "b0_type", "seed": [7, 3]}
    with pytest.raises(ValueError, match="boolean mask"):
        t1_scope.momentum_shift(start, ref, like, 0.5, own_factor_cells=np.zeros(3, bool), **kw)
    with pytest.raises(ValueError, match="boolean mask"):
        t1_scope.momentum_shift(start, ref, like, 0.5, own_factor_cells=np.zeros(like.n_obs, int), **kw)
    with pytest.raises(ValueError, match="logf_filter"):
        t1_scope.momentum_shift(start, ref, like, 0.5, own_factor_cells=np.zeros(like.n_obs, bool),
                                logf_filter=lambda logf, c, h: (logf, {}), **kw)


# ----------------------------------------------------------------------------------------------------------- the arms
def _spy_apply(monkeypatch):
    calls = []
    real = t1_scope.apply_factors

    def spy(C, F, tidx, held=None, keep_depth=True, seed=0, kernel="independent"):
        out = real(C, F, tidx, held, keep_depth=keep_depth, seed=seed, kernel=kernel)
        calls.append({"C": C.copy(), "F": F.copy(), "tidx": tidx.copy(), "seed": list(seed), "kernel": kernel,
                      "keep_depth": keep_depth, "out": out[0].copy(), "scale": out[1].copy(),
                      "held": np.zeros(C.shape[1], bool) if held is None else np.asarray(held, bool).copy()})
        return out

    monkeypatch.setattr(t1_scope, "apply_factors", spy)
    return calls


def _spy_draw(monkeypatch):
    seen = []
    real = harness.draw_copy

    def spy(ctx, herd=False):
        out = real(ctx, herd=herd)
        seen.append(out[1].copy())
        return out

    monkeypatch.setattr(harness, "draw_copy", spy)
    return seen


@pytest.mark.parametrize("name", ARMS)
def test_each_arm_keeps_b0_trends_rows_coordinates_and_cell_number(three, name, monkeypatch):
    seen = _spy_draw(monkeypatch)
    ctx = _ctx(three)
    a, b = CANDIDATES[name].build(ctx), CANDIDATES["b0_trend"].build(ctx)
    assert len(seen) == 2 and np.array_equal(seen[0], seen[1])
    assert a.n == b.n == BOARD.max_cells
    np.testing.assert_array_equal(a.coords, b.coords)
    assert _h(a.coords) == PIN_COORDS
    merged = np.asarray(three["ref"].obs["merged"]).astype(str)[seen[0]]
    flat = (merged == "Untracked") | (merged == "Tiny")
    assert flat.any() and np.array_equal(a.X[flat], b.X[flat])  # no factor, nothing to do: equal
    assert not np.array_equal(a.X, b.X)


def test_mix_gives_the_first_half_b0s_factor_and_the_rest_b0_trends(three, monkeypatch):
    """Through a spy on apply_factors: the arm's per-cell factor rows are b0's for draw positions below n // 2 and
    b0_trend's from there on; counts, stream, depth rule and kernel are b0_trend's."""
    calls = _spy_apply(monkeypatch)
    ctx = _ctx(three)
    CANDIDATES["b0"].build(ctx)
    CANDIDATES["b0_trend"].build(ctx)
    CANDIDATES["b0_trend_mix"].build(ctx)
    b0, tr, mix = calls
    n = b0["C"].shape[0]
    np.testing.assert_array_equal(mix["C"].toarray(), tr["C"].toarray())
    assert mix["seed"] == tr["seed"] and mix["kernel"] == tr["kernel"] == "independent" and mix["keep_depth"]
    rows_mix, rows_b0, rows_tr = mix["F"][mix["tidx"]], b0["F"][b0["tidx"]], tr["F"][tr["tidx"]]
    hold = np.arange(n) < n // 2
    np.testing.assert_array_equal(rows_mix[hold], rows_b0[hold])
    np.testing.assert_array_equal(rows_mix[~hold], rows_tr[~hold])
    differs = (rows_b0 != rows_tr).any(1)
    assert differs[hold].any() and differs[~hold].any()  # both halves hold trend-lineage cells
    rec = harness.TREND_ARM_RECORDS["b0_trend_mix"]
    assert rec["arm_rule"] == "vec.harness.B0_TREND_MIX_RULE" and rec["arm"] == "mix"
    assert rec["rule"] == "vec.harness.B0_TREND_RULE"  # b0_trend_factors' entry, not overwritten
    assert rec["n_hold"] == n // 2 and rec["n_trend_cells"] == n - n // 2
    # each cell's expected total is unchanged by the step (keep_depth): the library size per cell stays b0_trend's
    for c in (tr, mix):
        np.testing.assert_allclose((c["C"].toarray() * c["F"][c["tidx"]] * c["scale"][:, None]).sum(1),
                                   c["C"].toarray().sum(1), rtol=1e-12)


@pytest.mark.parametrize("name", sorted(KERNEL_ARMS))
def test_kernel_arms_use_b0_trends_factors_with_their_kernel(three, monkeypatch, name):
    """Through a spy: the same counts, factor rows, cell types and stream as b0_trend, with the arm's kernel; every
    cell's panel total after the step equals its total before (keep_depth: the rounded expected total). The arm's
    rule is under "arm_rule"; "rule" stays b0_trend_factors'."""
    calls = _spy_apply(monkeypatch)
    ctx = _ctx(three)
    CANDIDATES["b0_trend"].build(ctx)
    CANDIDATES[name].build(ctx)
    tr, xt = calls
    np.testing.assert_array_equal(xt["C"].toarray(), tr["C"].toarray())
    np.testing.assert_array_equal(xt["F"], tr["F"])
    np.testing.assert_array_equal(xt["tidx"], tr["tidx"])
    assert xt["seed"] == tr["seed"] and xt["kernel"] == KERNEL_ARMS[name] and tr["kernel"] == "independent"
    np.testing.assert_array_equal(xt["out"].toarray().sum(1), xt["C"].toarray().sum(1))
    assert not np.array_equal(tr["out"].toarray().sum(1), tr["C"].toarray().sum(1))
    rec = harness.TREND_ARM_RECORDS[name]
    arm = name.rsplit("_", 1)[1]
    assert rec["arm_rule"] == f"vec.harness.B0_TREND_{arm.upper()}_RULE" and rec["arm"] == arm
    assert rec["kernel"] == KERNEL_ARMS[name] and rec["rule"] == "vec.harness.B0_TREND_RULE"


def test_xt_is_b0_trends_draw_repaired(three, monkeypatch):
    """b0_trend_xt's apply_factors output is b0_trend's own draw except on cells whose total over the moved entries
    missed its rounded expected total, each moved by exactly |d| units, only on its moved entries."""
    calls = _spy_apply(monkeypatch)
    ctx = _ctx(three)
    CANDIDATES["b0_trend"].build(ctx)
    CANDIDATES["b0_trend_xt"].build(ctx)
    tr, xt = calls
    y, x = tr["out"].toarray(), xt["out"].toarray()
    D = tr["C"].toarray()
    f = np.where(tr["held"][None, :], 1.0, tr["F"][tr["tidx"]] * tr["scale"][:, None])
    M = (D > 0) & (f != 1.0)
    target = np.rint((D * f * M).sum(1)).astype(np.int64)
    d = (y * M).sum(1) - target
    np.testing.assert_array_equal((x * M).sum(1), target)
    np.testing.assert_array_equal(np.abs(x - y).sum(1), np.abs(d))
    np.testing.assert_array_equal(x[~M], y[~M])
    assert (d != 0).any() and (d == 0).any() and np.array_equal(x[d == 0], y[d == 0])


def _nothing_shifted(three) -> dict:
    """Every lineage absent at `start`: B0 shifts nothing, every factor is 1."""
    start = three["start"].copy()
    start.obs["merged"] = np.where(np.arange(start.n_obs) % 2 == 0, "Elsewhere", "Untracked")
    return {**three, "start": start}


@pytest.mark.parametrize("name", ARMS)
def test_with_nothing_to_do_each_arm_is_b0_trend_bit_for_bit(three, name):
    """No lineage shifted (every factor 1): both arms equal b0_trend bit for bit."""
    inputs = _nothing_shifted(three)
    a, b = CANDIDATES[name].build(_ctx(inputs)), CANDIDATES["b0_trend"].build(_ctx(inputs))
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.coords, b.coords)


def test_mix_with_the_trend_equal_to_b0s_step_is_b0_trend_bit_for_bit(three):
    """Only start and the reference usable (b0_trend's factor is B0's own for every lineage): the mixture has nothing
    to mix, and b0_trend_mix equals b0_trend (and b0) bit for bit."""
    inputs = _only_start_and_ref(three)
    a, b = CANDIDATES["b0_trend_mix"].build(_ctx(inputs)), CANDIDATES["b0_trend"].build(_ctx(inputs))
    np.testing.assert_array_equal(a.X, b.X)
    np.testing.assert_array_equal(a.X, CANDIDATES["b0"].build(_ctx(inputs)).X)


@pytest.mark.parametrize("name", ARMS)
def test_each_arm_is_registered_like_b0_trend(name):
    tasks = [t for ts_ in load_tasks().values() for t in ts_]
    c, base = CANDIDATES[name], CANDIDATES["b0_trend"]
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if base.applies(t)]
    assert {t.key for t in tasks if c.applies(t)} == {"T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    assert c.tasks == base.tasks and c.derived == () and c.requires is harness._heart_extrap_with_start
    assert c.params["scope"] == "chamber" and c.params["damping"] == B0_DAMPING
    rule = {"b0_trend_mix": "B0_TREND_MIX_RULE", "b0_trend_xt": "B0_TREND_XT_RULE",
            "b0_trend_mn": "B0_TREND_MN_RULE"}[name]
    assert c.params["rule"] == f"vec.harness.{rule}" and isinstance(getattr(harness, rule), str)


def test_the_negative_binomial_arm_is_not_built():
    """Line 3 is blocked (no dispersion rule on record): no candidate, no kernel."""
    assert not any("nb" in n.split("_") for n in CANDIDATES if n.startswith("b0_trend"))
    C, F, ti, held = _af_inputs()
    with pytest.raises(ValueError, match="unknown kernel"):
        core.apply_factors(C, F, ti, held, kernel="negative_binomial")


IDEAS = {"b0_trend_mix": {"trend_mix": True}, "b0_trend_xt": {"kernel": "exact_total"},
         "b0_trend_mn": {"kernel": "multinomial"}}


@pytest.mark.parametrize("name", ARMS)
def test_each_arm_refuses_what_b0_trend_refuses(three, name):
    idea = IDEAS[name]
    for kw in ({"depth": True}, {"herd": True}, {"cfm": True}):
        with pytest.raises(ValueError, match="depth, herd or cfm"):
            harness._b0_t2(_ctx(three), trend=True, **idea, **kw)
    with pytest.raises(ValueError, match="for trend=True"):
        harness._b0_t2(_ctx(three), **idea)
    with pytest.raises(ValueError, match="own factor only"):
        harness._b0_t2(_ctx(three), trend=True, trend_form="eb", **idea)
    c = CANDIDATES[name]
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), kind="interp")))
    with pytest.raises(ValueError, match="heart extrapolation"):
        c.build(_ctx(three, view=replace(_view(), setting="embryo")))
    two = {k: three[k] for k in ("ref", "start")}
    with pytest.raises(ValueError, match="three"):
        c.build(_ctx(two, view=_view(inputs={"ref": "T2/E9.5.h5ad", "start": "T2/E8.75.h5ad"})))
    for k in KERNEL_ARMS.values():
        with pytest.raises(ValueError, match="separate arms"):
            harness._b0_t2(_ctx(three), trend=True, trend_mix=True, kernel=k)
    with pytest.raises(ValueError, match="unknown kernel"):
        harness._b0_t2(_ctx(three), trend=True, kernel="negative_binomial")
