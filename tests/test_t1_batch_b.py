"""T1 batch t1-b (experiments/harness/TRIES.md, registered before any run), on synthetic data only: the tail decoder
(`b0_tail`, `vec.t1_scope.apply_tail_decoder`) and the posterior-predictive top-up (`b0_topup`,
`vec.t1_scope.topup_profile`, `vec.core.apply_factors_profile`). Each arm is B0 with one idea; with the idea idle it is
B0 bit for bit, and the code paths every existing candidate takes are pinned by checksums taken before the edit."""

from __future__ import annotations

import hashlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_t1_strong_null import GENES, LABELS, _ctx, _labelled

from vec import core, harness
from vec import t1_scope as ts
from vec.harness import PseudoTask
from vec.prep import log1p_cp10k

TASK = PseudoTask(key="T1/E9.5->E10.5", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                  reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                  validation_board="T1:val")
DERIVED = {"t1_held_flat": ("g49",)}
B0_TOPUP_SHA = "5f4403e64f88e5e0cf4eeb0a7f925d54d060998093e4d92d17bce0b94a1ef604"  # at 90a50d5, before t1-bfix


def _sha(a) -> str:
    a = np.ascontiguousarray(np.asarray(a.toarray() if sparse.issparse(a) else a))
    return hashlib.sha256(a.tobytes()).hexdigest()


def _factor_case():
    rng = np.random.default_rng(5)
    C = sparse.csr_matrix(rng.poisson(rng.gamma(0.5, 2.0, 40), (60, 40)))
    F = np.exp(rng.normal(0, 0.5, (3, 40)))
    ti = rng.integers(0, 3, 60)
    held = np.zeros(40, bool)
    held[3] = True
    return C, F, ti, held


def _build(name, ref, start, derived=DERIVED):
    return harness.CANDIDATES[name].build(_ctx(TASK, {"ref": ref, "start": start}, derived))


# ---------------------------------------------------------------------------------------------------------------------
# Today's paths, pinned before the edit (checksums taken at f63f964, before any line of this batch)

def test_apply_factors_without_a_profile_is_todays_code():
    C, F, ti, held = _factor_case()
    out, s, _ = core.apply_factors(C, F, ti, held, keep_depth=True, seed=[7, 1])
    assert _sha(out.astype(np.int64)) == "9bf920c4f8d40b3224bba8584145538a91efbf6b7a52d81ff7139a1f40477468"
    assert _sha(s) == "a726299fface640107ef74bda15dcf2b8738cf181c424e2d7d0cac9c5b2d9607"


def test_b0_is_todays_build():
    ref, start = _labelled(1500, LABELS, 1, GENES), _labelled(1200, LABELS, 2, GENES)
    st = _build("b0", ref, start)
    assert _sha(np.asarray(st.X, dtype=np.float32)) == (
        "23fab0fea7d206ed9d5f35050274d70e9c23ef5941888e869e8faf2a17a82508")


# ---------------------------------------------------------------------------------------------------------------------
# Registration and refusals

def test_both_arms_are_registered_on_b0_strong_nulls_tasks():
    tasks = harness.load_tasks()["T1:test"]
    ref = harness.CANDIDATES["b0_strong_null"]
    for name in ("b0_tail", "b0_topup"):
        c = harness.CANDIDATES[name]
        assert c.tasks == ref.tasks and c.derived == ref.derived == ("t1_held_flat",)
        assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if ref.applies(t)]
        assert c.params["damping"] == harness.B0_DAMPING and c.params["scope"] == "chamber"
        assert c.externals == () and c.tuned_on == frozenset()
    keys = {t.key for t in tasks if harness.CANDIDATES["b0_tail"].applies(t)}
    assert {"T1/E9.5->E10.5", "T1/E9.5->E10.5:chamber"} <= keys and "T1/E8.5->E9.5" not in keys


def test_the_constants_are_the_records():
    assert ts.TAIL_MIN_PB == ts.DP_MIN_PB == 0.25 and ts.TAIL_BOTTOM == ts.DP_BOTTOM == 400
    assert ts.TAIL_TARGET_DP == 0.0 and ts.TAIL_F_MAX == 1.0 / ts.DP_F_MIN == 2.0  # the bound again (t1-bfix2)
    assert (ts.TOPUP_N_PCS, ts.TOPUP_K) == (30, 15)
    p_tail, p_top = harness.CANDIDATES["b0_tail"].params, harness.CANDIDATES["b0_topup"].params
    assert (p_tail["bottom"], p_tail["min_pb"], p_tail["target_dp"], p_tail["f_max"]) == (400, 0.25, 0.0, 2.0)
    assert "more than twice their counts" in p_tail["what"] and "skipped" in p_tail["what"]
    assert (p_top["n_pcs"], p_top["k"]) == (30, 15)


def test_profile_refusals():
    C, F, ti, held = _factor_case()
    with pytest.raises(ValueError, match="shape"):
        core.apply_factors(C, F, ti, held, profile=sparse.csr_matrix((5, 40)))
    bad = sparse.csr_matrix(-np.ones((60, 40)))
    with pytest.raises(ValueError, match="non-negative"):
        core.apply_factors(C, F, ti, held, profile=bad)
    ref = _labelled(200, LABELS, 3, GENES)
    with pytest.raises(KeyError, match="source_obs"):
        ts.topup_profile(ref, ref)
    cells = ref[:5].copy()
    cells.obs["source_obs"] = ["nope"] * 5
    with pytest.raises(KeyError, match="not rows"):
        ts.topup_profile(ref, cells)


def test_the_tail_decoder_refuses_a_missing_panel_gene():
    pop = _labelled(100, LABELS, 4, GENES)
    with pytest.raises(KeyError, match="missing"):
        ts.apply_tail_decoder(pop, np.zeros(2), ["g0", "nope"], (), seed=1)


# ---------------------------------------------------------------------------------------------------------------------
# The exact top-up (b0_topup)

def test_the_depth_constant_is_exact_and_reduces_to_b0s_with_the_own_counts():
    rng = np.random.default_rng(0)
    for _ in range(50):
        c = rng.poisson(2.0, 30).astype(float)
        f = np.exp(rng.normal(0, 0.7, 30))
        m = rng.gamma(1.0, 2.0, 30)
        s = core._profile_depth_constant(c, f, m)
        fe = f * s
        total = np.where(fe <= 1, c * fe, c + (fe - 1) * m).sum()
        assert total == pytest.approx(c.sum(), rel=1e-12)
        s_own = core._profile_depth_constant(c, f, c)
        assert s_own == pytest.approx(c.sum() / (c * f).sum(), rel=1e-12)


def test_with_the_own_counts_as_profile_the_law_is_b0s():
    """m = c: each entry's law is apply_factors' (the same constant, binomial below 1, c + Poisson((f - 1) c) above), so
    the means over many draws agree, and no zero is ever raised."""
    C, F, ti, held = _factor_case()
    a = np.zeros(C.shape[1])
    b = np.zeros(C.shape[1])
    for k in range(200):
        a += np.asarray(core.apply_factors(C, F, ti, held, seed=[k])[0].sum(0)).ravel()
        out = core.apply_factors(C, F, ti, held, seed=[k], profile=C.astype(np.float64))[0].toarray()
        assert not ((C.toarray() == 0) & (out > 0)).any()
        b += out.sum(0)
    big = a > 2000
    assert big.sum() > 10 and np.abs(b[big] / a[big] - 1).max() < 0.05
    _, s0, _ = core.apply_factors(C, F, ti, held)
    _, s1, _ = core.apply_factors(C, F, ti, held, profile=C.astype(np.float64))
    assert np.allclose(s0, s1, rtol=1e-12)


def test_the_profile_raises_detection_and_realises_a_sparse_rise():
    rng = np.random.default_rng(1)
    n, g = 4000, 20
    C = sparse.csr_matrix(rng.poisson(np.r_[0.05, np.full(g - 1, 5.0)], (n, g)))
    F = np.ones((1, g))
    F[0, 0] = 3.0  # the sparse gene is meant to triple
    ti = np.zeros(n, dtype=np.int64)
    m = sparse.csr_matrix(np.tile(np.r_[0.05, np.full(g - 1, 5.0)], (n, 1)))  # its posterior-predictive mean
    own, _, _ = core.apply_factors(C, F, ti, keep_depth=True, seed=2)
    pp, s_pp, _ = core.apply_factors(C, F, ti, keep_depth=True, seed=2, profile=m)
    detected = (C[:, 0] > 0).sum()
    assert (own[:, 0] > 0).sum() <= detected  # the own-count kernel never raises detection
    assert (pp[:, 0] > 0).sum() > 1.8 * detected  # the posterior-predictive top-up does
    expected = C[:, 0].toarray().ravel() + (3.0 * s_pp - 1.0) * 0.05
    assert pp[:, 0].toarray().mean() == pytest.approx(expected.mean(), rel=0.15)
    assert pp[:, 0].toarray().mean() > 2.4 * C[:, 0].toarray().mean()  # most of the intended tripling is realised
    # the expected library is kept exactly, and the realised one within 1% on average
    lib = np.asarray(C.sum(1)).ravel()
    assert np.abs(np.asarray(pp.sum(1)).ravel().mean() / lib.mean() - 1) < 0.01


def test_held_genes_and_flat_factors_are_untouched_by_the_profile():
    C, F, ti, held = _factor_case()
    rng = np.random.default_rng(3)
    m = sparse.csr_matrix(rng.gamma(1.0, 1.0, C.shape))
    out, _, _ = core.apply_factors(C, F, ti, held, profile=m, seed=4)
    assert (out[:, 3] != C[:, 3]).nnz == 0
    ones = np.ones_like(F)
    out1, s1, _ = core.apply_factors(C, ones, ti, held, profile=m, seed=4)
    assert (out1 != C).nnz == 0 and np.all(s1 == 1.0)


def _pair(seed_ref=1, seed_start=2, n_ref=1500, n_start=1200):
    return _labelled(n_ref, LABELS, seed_ref, GENES), _labelled(n_start, LABELS, seed_start, GENES)


def test_topup_arm_is_b0_bit_for_bit_when_the_momentum_is_flat():
    ref, _ = _pair()
    start = ref.copy()  # no change between the stages: every factor is exactly 1, the idea has nothing to do
    a, b = _build("b0", ref, start), _build("b0_topup", ref, start)
    assert _sha(np.asarray(a.X, np.float32)) == _sha(np.asarray(b.X, np.float32))


def test_topup_arm_keeps_b0s_cells_and_depth_and_changes_only_the_rises(monkeypatch):
    ref, start = _pair()
    seen = []
    real = ts.momentum_shift

    def spy(*args, **kw):
        out, rep = real(*args, **kw)
        seen.append((kw.get("topup_source") is not None, out.obs.copy(),
                     np.asarray(out.layers["counts"].sum(1)).ravel(), np.asarray(ts.counts_of(args[2]).sum(1)).ravel()))
        return out, rep

    monkeypatch.setattr(ts, "momentum_shift", spy)
    a, b = _build("b0", ref, start), _build("b0_topup", ref, start)
    assert a.X.shape == b.X.shape and (np.asarray(a.X) != np.asarray(b.X)).any()
    (fa, obs_a, lib_a, in_a), (fb, obs_b, lib_b, in_b) = seen
    assert not fa and fb
    pd.testing.assert_frame_equal(obs_a, obs_b)  # the same drawn cells, categories and order
    assert np.array_equal(in_a, in_b)
    assert abs(lib_b.mean() / in_b.mean() - 1) < 0.01 and abs(lib_a.mean() / in_a.mean() - 1) < 0.01
    assert np.mean(np.abs(lib_b / in_b - 1)) < 0.05


def test_topup_profile_is_the_kernels_mean_of_neighbour_proportions():
    ref = _labelled(300, LABELS, 6, GENES)
    rows = np.array([0, 5, 5, 17])
    cells = ts._take(ref, rows)  # as reweight draws them: repeats allowed, obs names cell_0, cell_1, ...
    cells.obs["source_obs"] = np.asarray(ref.obs_names)[rows]
    P = ts.topup_profile(ref, cells).toarray()
    C = ts.counts_of(ref).toarray().astype(float)
    lib = C.sum(1)
    assert np.allclose(P.sum(1), lib[rows])  # expected counts add up to the cell's own library
    assert np.allclose(P[1], P[2])  # the same source row, the same profile
    assert (P[np.asarray(ts.counts_of(cells).toarray()) > 0] > 0).all()  # the cell itself is among its neighbours
    assert np.allclose(ts.topup_profile(ref, cells).toarray(), P)  # deterministic


# ---------------------------------------------------------------------------------------------------------------------
# The tail decoder (b0_tail)

def _tail_population(n=3000, seed=0):
    """Panel g0..g39; g0-g4 sparse (pb well below 0.25) and predicted down against pb_ref; the rest feasible."""
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(40)]
    mu = np.r_[np.full(5, 0.02), rng.gamma(4.0, 3.0, 35)]
    C = rng.poisson(mu, (n, 40))
    C[np.arange(n), 10] += 1
    pop = ad.AnnData(X=log1p_cp10k(sparse.csr_matrix(C)), layers={"counts": sparse.csr_matrix(C).astype(np.int32)},
                     obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    pop.var_names = genes
    pb = np.asarray(pop.X.mean(0)).ravel().astype(np.float64)
    pb_ref = pb.copy()
    pb_ref[:5] += 0.005  # the sparse genes fell (a scope-driven fall), still infeasible: pb_ref < 0.25
    pb_ref[5:10] += 0.5  # five feasible genes fell further
    return pop, genes, pb_ref


def test_the_decoder_lifts_infeasible_genes_to_flat_and_keeps_the_depth():
    pop, genes, pb_ref = _tail_population()
    out, rep = ts.apply_tail_decoder(pop, pb_ref, genes, ("g20",), seed=[3, 10])
    assert set(rep["lifted"]) <= {f"g{i}" for i in range(5)} and len(rep["lifted"]) >= 3
    assert rep["n_infeasible_in_bottom"] >= 5
    dp = np.asarray(out.X.mean(0)).ravel() - pb_ref
    lifted = [genes.index(g) for g in rep["lifted"]]
    assert np.all(np.abs(dp[lifted]) < 0.01)  # flat, within sampling
    assert np.all(dp[5:10] < -0.3)  # feasible genes keep their fall
    C0, C1 = ts.counts_of(pop), ts.counts_of(out)
    assert (C1[:, 20] != C0[:, 20]).nnz == 0  # held untouched
    assert abs(rep["depth_ratio_mean"] - 1) < 0.01 and rep["expected_depth_max_error"] < 1e-9
    assert out.n_obs == pop.n_obs and list(out.obs_names) == list(pop.obs_names)
    out2, _ = ts.apply_tail_decoder(pop, pb_ref, genes, ("g20",), seed=[3, 10])
    assert (ts.counts_of(out2) != C1).nnz == 0  # deterministic


def test_the_decoder_does_nothing_without_an_infeasible_gene_in_the_tail():
    pop, genes, pb_ref = _tail_population()
    pb_ref = pb_ref.copy()
    pb_ref[:5] = np.asarray(pop.X.mean(0)).ravel()[:5] - 0.01  # the sparse genes rise: none is in the down tail
    out, rep = ts.apply_tail_decoder(pop, pb_ref, genes, (), seed=1)
    assert out is pop and rep["n_lifted"] == 0
    pb_ref2 = np.asarray(pop.X.mean(0)).ravel() + 0.3  # every gene falls and every gene is feasible
    out, rep = ts.apply_tail_decoder(pop, pb_ref2, genes, (), seed=1)
    assert out is pop and rep["n_lifted"] == 0 and rep["n_infeasible_in_bottom"] == 0


def test_a_lift_out_of_reach_is_skipped_and_counted():
    pop, genes, pb_ref = _tail_population()
    pb_ref = pb_ref.copy()
    pb_ref[0] += 0.05  # g0 would need far more than f_max
    _, rep = ts.apply_tail_decoder(pop, pb_ref, genes, (), seed=1)
    assert "g0" not in rep["lifted"] and rep["n_lifted"] >= 3
    assert rep["n_skipped"] == rep["n_infeasible_in_bottom"] - rep["n_lifted"] >= 1
    assert 1.0 <= rep["median_factor"] <= rep["max_factor"] <= ts.TAIL_F_MAX
    C = ts.counts_of(pop).tolil()
    C[:, 1] = 0  # g1 has no counts in the file: slope zero, no answer
    pop0 = pop.copy()
    pop0.layers["counts"] = sparse.csr_matrix(C).astype(np.int32)
    pop0.X = log1p_cp10k(sparse.csr_matrix(C))
    pb0 = pb_ref.copy()
    pb0[1] = 0.01  # still fell and infeasible
    _, rep0 = ts.apply_tail_decoder(pop0, pb0, genes, (), seed=1)
    assert not {"g0", "g1"} & set(rep0["lifted"]) and rep0["n_skipped"] >= 2
    assert rep0["n_lifted"] + rep0["n_skipped"] == rep0["n_infeasible_in_bottom"]


def _sparse_case():
    """g0 very sparse (needs far more than twice its counts), g1-g4 modestly sparse and modestly fallen, the rest flat
    and well expressed."""
    rng = np.random.default_rng(11)
    n, g = 3000, 40
    genes = [f"g{i}" for i in range(g)]
    mu = np.r_[0.003, np.full(4, 0.05), rng.gamma(4.0, 3.0, g - 5)]
    C = rng.poisson(mu, (n, g))
    C[np.arange(n), 10] += 1
    Cs = sparse.csr_matrix(C)
    pop = ad.AnnData(X=log1p_cp10k(Cs), layers={"counts": Cs.astype(np.int32)},
                     obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    pop.var_names = genes
    pb_ref = np.asarray(pop.X.mean(0)).ravel().astype(np.float64)  # every other gene flat: dp = 0, not lifted
    pb_ref[0] += 0.01  # g0 fell: still infeasible (pb_ref < 0.25), out of reach
    pb_ref[1:5] += 0.003  # g1-g4 fell a little: within reach
    return pop, genes, pb_ref


def test_a_very_sparse_gene_is_skipped_and_the_modest_ones_are_lifted():
    pop, genes, pb_ref = _sparse_case()
    out, rep = ts.apply_tail_decoder(pop, pb_ref, genes, (), seed=[5, 10])
    assert sorted(rep["lifted"]) == ["g1", "g2", "g3", "g4"]
    assert rep["n_infeasible_in_bottom"] == 5 and rep["n_lifted"] == 4 and rep["n_skipped"] == 1
    assert 1.0 < rep["median_factor"] <= rep["max_factor"] <= ts.TAIL_F_MAX
    C0, new = ts.counts_of(pop), ts.counts_of(out)
    vals = new.toarray()
    assert np.issubdtype(new.dtype, np.integer) and np.all(vals >= 0)
    assert np.all(np.isfinite(np.asarray(out.X.toarray())))
    lib0, lib1 = np.asarray(C0.sum(1)).ravel(), np.asarray(new.sum(1)).ravel()
    assert np.mean(np.abs(lib1 / lib0 - 1)) < 0.01
    dp = np.asarray(out.X.mean(0)).ravel() - pb_ref
    assert np.all(np.abs(dp[1:5]) < 0.01)  # the modest ones land flat
    line = harness.tail_check_line(TASK.key, 1, rep)
    assert "infeasible_in_bottom=5 lifted=4 skipped=1 " in line and not any(g in line.split() for g in genes)


# Checksums (counts, X) of apply_tail_decoder's output on the synthetic cases, taken with 90a50d5's code before the
# t1-bfix2 edit: the bounded arithmetic is 90a50d5's bit for bit.
DECODER_PINS = {
    "tail": ("6398ef7ad80f439bdc904ca2d70046b1d59207c1be71ab4c1d2bc2a09bdb0a08",
             "3f5a28ac305555a3318582830fffe74fde4c0ccfe6d82334149167f00f7683aa"),
    "out_of_reach": ("a3abe7b81b922cbd23b215a738cdfc04eeea5eab4e59aed71521818e77f0daba",
                     "9c1d2bc87cd03a8897d53c3c998a3d0df12af44084a596987cb3987cb9ea7e0c"),
    "sparse": ("4b327c4188b882a7cc8c308f1b01a432578205440018e4f85a47df433a4001a7",
               "bf5217ed7c2ab518129cd6077facbdb66de6539d93df237b9568f5f6b8e0aa8c"),
}


def _decoder_sha(out):
    return _sha(ts.counts_of(out).astype(np.int64)), _sha(np.asarray(out.X.toarray(), np.float32))


def test_the_decoder_is_90a50d5s_bit_for_bit():
    pop, genes, pb_ref = _tail_population()
    assert _decoder_sha(ts.apply_tail_decoder(pop, pb_ref, genes, ("g20",), seed=[3, 10])[0]) == DECODER_PINS["tail"]
    p2 = pb_ref.copy()
    p2[0] += 0.05
    assert _decoder_sha(ts.apply_tail_decoder(pop, p2, genes, (), seed=1)[0]) == DECODER_PINS["out_of_reach"]
    sp, sg, sr = _sparse_case()
    assert _decoder_sha(ts.apply_tail_decoder(sp, sr, sg, (), seed=[5, 10])[0]) == DECODER_PINS["sparse"]


def test_the_check_line_holds_numbers_only():
    rep = {"n_infeasible_in_bottom": 7, "n_lifted": 5, "n_skipped": 2, "max_factor": 1.98765, "median_factor": 1.5,
           "lifted": ["g0"]}
    line = harness.tail_check_line("T1/E9.5->E10.5", 1001, rep)
    assert line == ("b0_tail check: task=T1/E9.5->E10.5 seed=1001 infeasible_in_bottom=7 lifted=5 skipped=2 "
                    "max_factor=1.99 median_factor=1.5")
    assert "g0" not in line
    idle = harness.tail_check_line("k", 1, {"n_infeasible_in_bottom": 0, "n_lifted": 0, "n_skipped": 0,
                                            "max_factor": float("nan"), "median_factor": float("nan")})
    assert idle.endswith("lifted=0 skipped=0 max_factor=nan median_factor=nan")


def test_tail_arm_is_b0_bit_for_bit_when_no_gene_is_infeasible(monkeypatch, capsys):
    ref, start = _pair()
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 0.0)  # strength zero: no gene can be infeasible
    a = _build("b0", ref, start)
    capsys.readouterr()
    b = _build("b0_tail", ref, start)
    assert _sha(np.asarray(a.X, np.float32)) == _sha(np.asarray(b.X, np.float32))
    lines = [x for x in capsys.readouterr().out.splitlines() if x.startswith("b0_tail check:")]
    assert lines == [f"b0_tail check: task={TASK.key} seed={harness.BUILD_SEED_OFFSET + 1} infeasible_in_bottom=0 "
                     "lifted=0 skipped=0 max_factor=nan median_factor=nan"]


def test_b0_topup_and_b0_tail_are_pinned_at_90a50d5(monkeypatch, capsys):
    """Checksums taken at 90a50d5, before the t1-bfix and t1-bfix2 edits: b0_topup is untouched, b0_tail with nothing
    to lift (the synthetic genes are all well expressed) is B0, and b0_tail with the tail made to count (TAIL_MIN_PB
    raised) is 90a50d5's build."""
    ref, start = _pair()
    assert _sha(np.asarray(_build("b0_topup", ref, start).X, np.float32)) == B0_TOPUP_SHA
    assert _sha(np.asarray(_build("b0_tail", ref, start).X, np.float32)) == (
        "23fab0fea7d206ed9d5f35050274d70e9c23ef5941888e869e8faf2a17a82508")
    assert "b0_tail check: " in capsys.readouterr().out
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 10.0)
    assert _sha(np.asarray(_build("b0_tail", ref, start).X, np.float32)) == (
        "f642c8d48e1ac0157b35469714d35650418f75954b7f5d73c0e259bb4324b629")


def test_tail_arm_keeps_b0s_cells_and_moves_only_the_decoders_genes(monkeypatch, capsys):
    ref, start = _pair()
    a = _build("b0", ref, start)
    capsys.readouterr()
    seen = []
    real = ts.apply_tail_decoder

    def spy(pop, pb_ref, genes, held, seed):
        out, rep = real(pop, pb_ref, genes, held, seed)
        seen.append((pop, out, rep, held, seed))
        return out, rep

    monkeypatch.setattr(ts, "apply_tail_decoder", spy)
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 10.0)  # synthetic genes are all well expressed: make the bottom tail count
    b = _build("b0_tail", ref, start)
    (pop, out, rep, held, seed), = seen
    assert list(held) == ["g49"] and seed == [harness.BUILD_SEED_OFFSET + 1, harness.TAIL_STREAM]
    assert rep["n_lifted"] > 0
    assert _sha(np.asarray(pop.X[:, :].toarray(), np.float32)) == _sha(np.asarray(a.X, np.float32))  # b0 before it
    assert a.X.shape == b.X.shape and list(out.obs_names) == list(pop.obs_names)
    lib0, lib1 = (np.asarray(ts.counts_of(x).sum(1)).ravel() for x in (pop, out))
    assert abs(lib1.mean() / lib0.mean() - 1) < 0.01
    assert (ts.counts_of(out)[:, GENES.index("g49")] != ts.counts_of(pop)[:, GENES.index("g49")]).nnz == 0
    lines = [x for x in capsys.readouterr().out.splitlines() if x.startswith("b0_tail check:")]
    assert lines == [harness.tail_check_line(TASK.key, seed[0], rep)]
    assert f"lifted={rep['n_lifted']} " in lines[0] and not any(g in lines[0] for g in rep["lifted"])
