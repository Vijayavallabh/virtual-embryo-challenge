"""T1 arms `b0_mn`, `b0_half` and `b0_half_ctrl` (campaign/t1-f): b0 on T1 with one change each to its momentum step,
B0's only count edit on T1 (a whole-cell multinomial redraw; the step on a random half of the drawn cells; the
mixture's mean-matched control), with their check lines and records. Synthetic data only; nothing scored, no
prediction compared with a target."""

from __future__ import annotations

import hashlib
import re

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import core, harness, t1_scope
from vec.boards import Board
from vec.harness import B0_DAMPING, CANDIDATES, BuildContext, PseudoTask, TaskView, load_tasks
from vec.prep import log1p_cp10k

GENES = [f"g{i}" for i in range(50)]
LABELS = ["V-CM", "Endocardium", "Blood", "Macrophage"]
HELD = ("g49",)
TASK = PseudoTask(key="T1/E9.5->E10.5", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                  reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                  validation_board="T1:val")
ARMS = ("b0_mn", "b0_half", "b0_half_ctrl")
REL = 1e-9  # the brief's tolerance on expectations computed from the factors (not a constant of any arm)


def _labelled(n, seed):
    rng = np.random.default_rng(seed)
    mu = rng.gamma(2.0, 8.0, len(GENES))
    C = rng.poisson(mu, (n, len(GENES))).astype(np.int32)
    C[np.arange(n), rng.integers(0, len(GENES), n)] += 1
    C[: n // 3, 5] = 0  # a gene with zeros in many cells
    a = ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), layers={"counts": sparse.csr_matrix(C)},
                   obs=pd.DataFrame({"celltype": np.resize(LABELS, n)}, index=[f"c{i}" for i in range(n)]))
    a.var_names = GENES
    return a


@pytest.fixture(scope="module")
def stages():
    return {"ref": _labelled(1500, 1), "start": _labelled(1200, 2)}


def _ctx(inputs, seed=harness.BUILD_SEED_OFFSET + 1, max_cells=700):
    board = Board(key="pseudo", task="T1", genes=list(GENES), needs_coords=False, min_cells=50, max_cells=max_cells)
    return BuildContext(task=TaskView.of(TASK), board=harness._candidate_board(board), inputs=inputs, seed=seed,
                        derived={"t1_held_flat": HELD})


def _flat(stages):
    """start = the reference's own cells: every type's momentum is exactly 0, every factor exactly 1."""
    return {"ref": stages["ref"], "start": stages["ref"].copy()}


class _RecordingGenerator(np.random.Generator):
    """numpy's Generator on the same bit generator (so the same stream), recording the parameters of the
    independent kernel's draws: binomial keep probabilities and Poisson means."""

    def __init__(self, bitgen, log):
        super().__init__(bitgen)
        self.log = log

    def binomial(self, n, p, size=None):
        self.log.setdefault("binomial", []).append((np.array(n, copy=True), np.array(p, dtype=np.float64, copy=True)))
        return super().binomial(n, p, size)

    def poisson(self, lam=1.0, size=None):
        self.log.setdefault("poisson", []).append(np.array(lam, dtype=np.float64, copy=True))
        return super().poisson(lam, size)


def _spy(monkeypatch):
    """Each apply_factors call: its inputs and outputs, and the parameters of every draw it made (the generator
    `default_rng` gives inside the call is recorded; the stream is the same)."""
    calls = []
    real = t1_scope.apply_factors
    real_rng = np.random.default_rng

    def spy(C, F, tidx, held=None, keep_depth=True, seed=0, **kw):
        log: dict = {}
        monkeypatch.setattr(core.np.random, "default_rng",
                            lambda s=None: _RecordingGenerator(np.random.PCG64(s), log))
        try:
            out = real(C, F, tidx, held, keep_depth=keep_depth, seed=seed, **kw)
        finally:
            monkeypatch.setattr(core.np.random, "default_rng", real_rng)
        calls.append({"C": sparse.csr_matrix(C, dtype=np.int64), "F": F.copy(), "tidx": tidx.copy(),
                      "seed": list(np.atleast_1d(seed)), "kw": dict(kw), "keep_depth": keep_depth,
                      "out": out[0].copy(), "scale": out[1].copy(), "dead": out[2].copy(), "draws": log,
                      "held": np.zeros(C.shape[1], bool) if held is None else np.asarray(held, bool).copy()})
        return out

    monkeypatch.setattr(t1_scope, "apply_factors", spy)
    return calls


def _b0_factor(call):
    """B0's per-entry factor F (dense, cells x genes): the type's factor times the cell's keep-depth constant, 1 on
    held genes (apply_factors' own arithmetic)."""
    return np.where(call["held"][None, :], 1.0, call["F"][call["tidx"]] * call["scale"][:, None])


def _digest(stage) -> str:
    X = stage.X.toarray() if sparse.issparse(stage.X) else np.asarray(stage.X)
    return hashlib.sha256(np.ascontiguousarray(X).tobytes()).hexdigest()


# ------------------------------------------------------------------------------------- 1. b0 and b0_xt are unchanged
# Taken at 80e553a (campaign/t1-d, before this batch's edit) on this file's synthetic stages and contexts.
PINS = {"b0": "74d3aab2745e1bc2b7d9bcb4c843dd8a12f9cb319c5e8083440ea3f58688d671",
        "b0_xt": "5a56a916d8249047ac2c316d1547640820e88ca781da0085dd57620bdae10ee7",
        "b0_flat": "027a9475fbdf982fa38d2f0d2931b27d1d2e7d2cc18b323603a6e8d99effcbfc"}


def test_b0_and_b0_xt_are_unchanged_bit_for_bit(stages):
    assert _digest(CANDIDATES["b0"].build(_ctx(stages))) == PINS["b0"]
    assert _digest(CANDIDATES["b0_xt"].build(_ctx(stages))) == PINS["b0_xt"]
    assert _digest(CANDIDATES["b0"].build(_ctx(_flat(stages)))) == PINS["b0_flat"]


def test_new_keywords_unset_change_nothing_in_the_core():
    rng = np.random.default_rng(31)
    C = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 2.0, 40), (60, 40)))
    F = np.exp(rng.normal(0, 0.6, (3, 40)))
    ti = rng.integers(0, 3, 60)
    for kernel in ("independent", "exact_total", "multinomial"):
        a = core.apply_factors(C, F, ti, seed=[5, 3], kernel=kernel)
        b = core.apply_factors(C, F, ti, seed=[5, 3], kernel=kernel, step_cells=None, step_weight=None)
        for x, y in zip(a, b, strict=True):
            np.testing.assert_array_equal(*(z.toarray() if sparse.issparse(z) else z for z in (x, y)))
    empty: dict = {}
    core.apply_factors(C, F, ti, seed=[5, 3], repair_record=empty)
    assert empty == {}  # the independent kernel without a step weight records nothing


# -------------------------------------------------------------------------------------------------------- 2. b0_mn
def test_b0_mn_is_a_by_hand_multinomial_redraw_of_b0s_step(stages, monkeypatch):
    calls = _spy(monkeypatch)
    ctx = _ctx(stages)
    CANDIDATES["b0"].build(ctx)
    CANDIDATES["b0_mn"].build(ctx)
    b0, mn = calls
    for k in ("F", "tidx", "held", "scale"):
        np.testing.assert_array_equal(mn[k], b0[k])
    np.testing.assert_array_equal(mn["C"].toarray(), b0["C"].toarray())
    assert mn["seed"] == b0["seed"] and mn["keep_depth"] and mn["kw"]["kernel"] == "multinomial"
    assert b0["kw"].get("kernel", "independent") == "independent"
    D = b0["C"].toarray().astype(np.float64)
    f = _b0_factor(b0)
    M = (f != 1.0) & (D > 0)
    assert M.any()
    want = D.copy()
    rng = np.random.default_rng(b0["seed"])
    for i in np.flatnonzero(M.any(1)):  # cells in row order, entries in gene order (CSR data order)
        w = (D[i] * f[i])[M[i]]
        want[i, M[i]] = rng.multinomial(int(np.rint(w.sum())), w / w.sum())
    assert not mn["dead"].any()
    np.testing.assert_array_equal(mn["out"].toarray(), want)
    # the expected totals are b0's: keep_depth gives each cell its drawn total, the redraw its rounded value
    expected = (D * f).sum(1)
    np.testing.assert_allclose(expected, D.sum(1), rtol=REL)
    off = mn["out"].toarray().sum(1) - expected
    assert np.abs(off).max() <= 0.5 + 1e-9  # rint(S) against S: at most half a count in general
    # with keep_depth S is the integer sum of the moved counts (up to float rounding), so the realised total is the
    # drawn total exactly: the arm moves no cell's library
    np.testing.assert_array_equal(mn["out"].toarray().sum(1), D.sum(1))
    rec = harness.T1_ARM_RECORDS["b0_mn"]
    moved = M.any(1)
    assert rec["moved_share"] == float(moved.mean())
    assert rec["mean_abs_total_off"] == pytest.approx(float(np.abs(off[moved]).mean()), rel=1e-12, abs=1e-12)
    assert rec["mean_abs_total_off"] < 1e-6


# ------------------------------------------------------------------------------------------------------ 3. b0_half
def _rows(A):
    return {tuple(r) for r in A}


def test_b0_half_moves_exactly_the_first_half_of_the_draw_by_b0s_factors(stages, monkeypatch):
    calls = _spy(monkeypatch)
    ctx = _ctx(stages)
    CANDIDATES["b0"].build(ctx)
    stage = CANDIDATES["b0_half"].build(ctx)
    b0, half = calls
    for k in ("F", "tidx", "held", "scale"):  # b0's factors and keep-depth constant, computed on all drawn cells
        np.testing.assert_array_equal(half[k], b0[k])
    np.testing.assert_array_equal(half["C"].toarray(), b0["C"].toarray())
    assert half["seed"] == b0["seed"] and half["kw"].get("kernel") == b0["kw"].get("kernel") == "independent"
    D, Y = half["C"].toarray(), half["out"].toarray()
    n = D.shape[0]
    k = n // 2
    step = half["kw"]["step_cells"]
    np.testing.assert_array_equal(step, np.arange(n) < k)
    differ = (Y != D).any(1)
    assert differ.sum() == k and np.array_equal(differ, step)
    np.testing.assert_array_equal(Y[~step], D[~step])  # the other cells: the drawn counts, bit for bit
    np.testing.assert_array_equal(np.asarray(stage.X)[~step], log1p_cp10k(sparse.csr_matrix(D[~step])).toarray())
    # the moved cells' draws are at b0's factors: binomial keep probabilities and Poisson means, entry for entry
    f = _b0_factor(b0)
    on = (D > 0) & step[:, None]
    (bn, bp), = half["draws"]["binomial"]
    lam, = half["draws"]["poisson"]
    Cs = half["C"]
    rows = np.repeat(np.arange(n), np.diff(Cs.indptr))
    fd = f[rows, Cs.indices]
    st = step[rows]
    np.testing.assert_array_equal(bp, fd[(fd < 1.0) & st])
    np.testing.assert_array_equal(lam, (fd[(fd > 1.0) & st] - 1.0) * Cs.data[(fd > 1.0) & st])
    assert on.any() and (f[on] != 1.0).any()


def test_b0_half_uses_the_draw_order_and_another_seed_moves_another_half(stages, monkeypatch):
    """reweight's draw ends in a uniform permutation of the drawn cells, so its first n // 2 are a uniformly random
    half; no new stream is drawn (the arm's apply_factors call is b0's, on b0's stream)."""
    calls = _spy(monkeypatch)
    CANDIDATES["b0_half"].build(_ctx(stages))
    CANDIDATES["b0_half"].build(_ctx(stages, seed=harness.BUILD_SEED_OFFSET + 2))
    a, b = calls
    k = a["C"].shape[0] // 2
    A, B = a["C"].toarray(), b["C"].toarray()
    assert _rows(A[:k]) != _rows(B[:k])
    # b0's draw (reweight, as build_submission calls it) in its row order: not grouped by category
    src = TaskView.of(TASK).role_stems["ref"]
    spec = t1_scope.load_spec()
    ref = harness._shallow(stages["ref"])
    shares = harness.t1_chamber_shares(ref, src)
    t1_scope.classify(ref, source=src, spec=spec)
    pop, _ = t1_scope.reweight(ref, shares, A.shape[0], harness.BUILD_SEED_OFFSET + 1, fallback=spec.fallbacks(src),
                               spec=spec)
    np.testing.assert_array_equal(t1_scope.counts_of(pop).toarray(), A)  # the same cells in the same order
    cats = pop.obs["category"].astype(str).to_numpy()
    assert len(set(cats)) > 1 and set(cats[:k]) == set(cats) == set(cats[k:])
    order = {c: i for i, c in enumerate(spec.categories)}
    assert not (np.diff([order[c] for c in cats]) >= 0).all()


# ------------------------------------------------------------------------------------------------- 4. b0_half_ctrl
def test_b0_half_ctrl_factor_is_one_plus_f_over_two_and_matches_the_mixture_in_expectation(stages, monkeypatch):
    calls = _spy(monkeypatch)
    ctx = _ctx(stages)
    CANDIDATES["b0"].build(ctx)
    CANDIDATES["b0_half_ctrl"].build(ctx)
    b0, ctrl = calls
    for k in ("F", "tidx", "held", "scale"):  # nothing recomputed
        np.testing.assert_array_equal(ctrl[k], b0[k])
    np.testing.assert_array_equal(ctrl["C"].toarray(), b0["C"].toarray())
    n = b0["C"].shape[0]
    assert n % 2 == 0 and ctrl["kw"]["step_weight"] == 0.5
    Cs = b0["C"]
    rows = np.repeat(np.arange(n), np.diff(Cs.indptr))
    F = _b0_factor(b0)[rows, Cs.indices]  # b0's per-entry factor, CSR order
    G = (1.0 + F) / 2.0
    (_, bp), = ctrl["draws"]["binomial"]
    lam, = ctrl["draws"]["poisson"]
    np.testing.assert_allclose(bp, G[F < 1.0], rtol=1e-15, atol=0)  # the spied per-entry factor is (1 + F) / 2
    np.testing.assert_allclose(lam, (G[F > 1.0] - 1.0) * Cs.data[F > 1.0], rtol=1e-15, atol=0)
    g_used = np.ones_like(F)
    g_used[F < 1.0] = bp
    g_used[F > 1.0] = lam / Cs.data[F > 1.0] + 1.0
    c = Cs.data.astype(np.float64)
    k = n // 2
    # expectations from the factors: per gene over the drawn cells, the mixture's (each cell moved with chance k / n)
    mix_gene = np.bincount(Cs.indices, weights=c * (1.0 + (k / n) * (F - 1.0)), minlength=Cs.shape[1]) / n
    ctrl_gene = np.bincount(Cs.indices, weights=c * g_used, minlength=Cs.shape[1]) / n
    np.testing.assert_allclose(ctrl_gene, mix_gene, rtol=REL)
    # per cell: the control's expected total is b0's
    b0_cell = np.bincount(rows, weights=c * F, minlength=n)
    ctrl_cell = np.bincount(rows, weights=c * g_used, minlength=n)
    np.testing.assert_allclose(ctrl_cell, b0_cell, rtol=REL)
    rec = harness.T1_ARM_RECORDS["b0_half_ctrl"]
    on = (F != 1.0) & (c > 0)
    assert rec["step_ratio"] == pytest.approx(float((np.abs(np.log(G[on])) / np.abs(np.log(F[on]))).mean()),
                                              rel=1e-12)
    assert 0.0 < rec["step_ratio"] < 1.0


def test_odd_count_control_weight_is_the_mixtures_share_and_the_means_match():
    """With an odd count the mixture moves k = n // 2 of n cells, so its per-gene expected mean weights the step by
    k / n; the control's weight is that share (1/2 only for an even count), so the means still match exactly."""
    rng = np.random.default_rng(7)
    n = 61
    C = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 3.0, 30), (n, 30)))
    F = np.exp(rng.normal(0, 0.5, (2, 30)))
    ti = rng.integers(0, 2, n)
    _, scale, _ = core.apply_factors(C, F, ti, seed=[3, 1])
    D = C.toarray().astype(np.float64)
    f = F[ti] * scale[:, None]
    w = (n // 2) / n
    log: dict = {}
    real = np.random.default_rng
    try:
        core.np.random.default_rng = lambda s=None: _RecordingGenerator(np.random.PCG64(s), log)
        core.apply_factors(C, F, ti, seed=[3, 1], step_weight=w)
    finally:
        core.np.random.default_rng = real
    G = 1.0 + w * (f - 1.0)
    (_, bp), = log["binomial"]
    Cs = sparse.csr_matrix(C, dtype=np.int64)
    rows = np.repeat(np.arange(n), np.diff(Cs.indptr))
    Gd = G[rows, Cs.indices]
    np.testing.assert_allclose(bp, Gd[Gd < 1.0], rtol=1e-15, atol=0)
    mix = (D * (1.0 + (n // 2) / n * (f - 1.0))).mean(0)
    np.testing.assert_allclose((D * G).mean(0), mix, rtol=REL)


# ------------------------------------------------------------------------------------------------- 5. an idle step
@pytest.mark.parametrize("name", ARMS)
def test_idle_step_is_b0_bit_for_bit(stages, name):
    inputs = _flat(stages)
    a, b = CANDIDATES[name].build(_ctx(inputs)), CANDIDATES["b0"].build(_ctx(inputs))
    np.testing.assert_array_equal(np.asarray(a.X), np.asarray(b.X))


def test_idle_records(stages):
    inputs = _flat(stages)
    for name in ARMS:
        CANDIDATES[name].build(_ctx(inputs))
    assert harness.T1_ARM_RECORDS["b0_mn"]["moved_share"] == 0.0
    assert harness.T1_ARM_RECORDS["b0_mn"]["mean_abs_total_off"] == 0.0
    assert harness.T1_ARM_RECORDS["b0_half_ctrl"]["step_ratio"] == 0.0


# --------------------------------------------------------------------------------------------------------- 6. seeds
@pytest.mark.parametrize("name", ARMS)
def test_same_seed_same_output_another_seed_another_draw(stages, name):
    a = CANDIDATES[name].build(_ctx(stages))
    b = CANDIDATES[name].build(_ctx(stages))
    c = CANDIDATES[name].build(_ctx(stages, seed=harness.BUILD_SEED_OFFSET + 2))
    np.testing.assert_array_equal(np.asarray(a.X), np.asarray(b.X))
    assert np.asarray(a.X).shape == np.asarray(c.X).shape and not np.array_equal(np.asarray(a.X), np.asarray(c.X))


# --------------------------------------------------------------------------------------------------- 7. check lines
CHECKS = {"b0_mn": (r"moved_share=(\d+\.\d{4}) mean_abs_total_off=(\d+\.\d{4})",
                    {"arm", "arm_rule", "kernel", "task", "seed", "moved_share", "mean_abs_total_off"}),
          "b0_half": (r"moved_share=(\d+\.\d{2})", {"arm", "arm_rule", "task", "seed", "moved_share"}),
          "b0_half_ctrl": (r"step_ratio=(\d+\.\d{4})", {"arm", "arm_rule", "task", "seed", "step_ratio"})}


@pytest.mark.parametrize("name", ARMS)
def test_check_line_holds_only_the_allowed_fields(stages, capsys, name):
    ctx = _ctx(stages)
    CANDIDATES[name].build(ctx)
    lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.startswith(f"{name} check:")]
    assert len(lines) == 1
    pat, keys = CHECKS[name]
    m = re.fullmatch(rf"{name} check: task=(\S+) seed=(\d+) " + pat, lines[0])
    assert m and m.group(1) == TASK.key and int(m.group(2)) == ctx.seed
    rec = harness.T1_ARM_RECORDS[name]
    assert set(rec) == keys and rec["arm_rule"] == CANDIDATES[name].params["rule"]
    assert not any(g in lines[0] for g in GENES) and "700" not in lines[0]
    if name == "b0_half":
        assert rec["moved_share"] == 0.5


# ---------------------------------------------------------------------------------------------------- registration
@pytest.mark.parametrize("name", ARMS)
def test_registered_for_t1_with_b0_xts_requirements(name):
    tasks = [t for ts_ in load_tasks().values() for t in ts_]
    c, base = CANDIDATES[name], CANDIDATES["b0_xt"]
    assert c.tasks == frozenset({"T1"}) and c.requires is base.requires and c.derived == base.derived
    assert c.needs == base.needs
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if base.applies(t)]
    assert {"T1/E9.5->E10.5", "T1/E9.5->E10.5:chamber"} <= {t.key for t in tasks if c.applies(t)}
    assert c.params["scope"] == "chamber" and c.params["damping"] == B0_DAMPING
    rule = c.params["rule"]
    assert rule.startswith("vec.harness.") and isinstance(getattr(harness, rule.rsplit(".", 1)[1]), str)


def test_the_arms_refuse_other_combinations():
    with pytest.raises(ValueError, match="unknown T1 arm"):
        harness._b0_t1(B0_DAMPING, arm="mix")
    with pytest.raises(ValueError, match="B0 alone"):
        harness._b0_t1(0.0, arm="mn")
    with pytest.raises(ValueError, match="B0 alone"):
        harness._b0_t1(B0_DAMPING, herd=True, arm="half")
    with pytest.raises(ValueError, match="B0 alone"):
        harness._b0_t1(B0_DAMPING, prior=lambda ctx: None, arm="half_ctrl")
    with pytest.raises(ValueError, match="B0 alone"):
        harness._b0_t1(B0_DAMPING, kernel="exact_total", arm="mn")
    C = sparse.csr_matrix(np.ones((4, 3), dtype=np.int64))
    with pytest.raises(ValueError, match="step_cells"):
        core.apply_factors(C, np.ones((1, 3)), np.zeros(4, int), step_cells=np.ones(3, bool))
    with pytest.raises(ValueError, match="step_weight"):
        core.apply_factors(C, np.ones((1, 3)), np.zeros(4, int), step_weight=0.0)
