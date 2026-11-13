"""T1 batch t1-c, on synthetic data only: candidate b0_tail_profile, the tail decoder of b0_tail with one change, the
kernel (`vec.t1_scope.apply_tail_decoder_profile`: c + Poisson((f - 1) m), m the neighbour profile of
`vec.t1_scope.topup_profile`, through `vec.core.apply_factors_profile`). With nothing to lift it is B0 bit for bit."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_t1_batch_b import TASK, _build, _pair, _sha, _sparse_case, _tail_population

from vec import harness
from vec import t1_scope as ts
from vec.prep import log1p_cp10k


def _own_source(pop: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """pop as its own source stage: each cell is the source row of the same name (as reweight records it)."""
    pop = pop.copy()
    pop.obs["source_obs"] = np.asarray(pop.obs_names)
    return pop, pop.copy()


def _spy_selection(monkeypatch) -> list:
    seen = []
    real = ts._tail_selection

    def spy(*args):
        res = real(*args)
        seen.append(np.array(res[-1]))
        return res

    monkeypatch.setattr(ts, "_tail_selection", spy)
    return seen


def _expected_selection(pop, pb_ref, genes, held) -> np.ndarray:
    """The selection of apply_tail_decoder's docstring, computed here from scratch."""
    cols = [list(pop.var_names).index(g) for g in genes]
    dp = np.asarray(sparse.csr_matrix(pop.X[:, cols], dtype=np.float64).sum(0)).ravel() / pop.n_obs - pb_ref
    bottom = np.argsort(dp, kind="stable")[:ts.TAIL_BOTTOM]
    held_panel = np.array([g in set(held) for g in genes])
    return bottom[(pb_ref[bottom] < ts.TAIL_MIN_PB) & (dp[bottom] < ts.TAIL_TARGET_DP) & ~held_panel[bottom]]


# ---------------------------------------------------------------------------------------------------------------------
# 1. The selection is b0_tail's

def test_the_selected_genes_are_the_own_count_decoders(monkeypatch):
    seen = _spy_selection(monkeypatch)
    for pop, genes, pb_ref, held in ((*_tail_population(), ("g20",)), (*_sparse_case(), ())):
        pop, src = _own_source(pop)
        _, rep_own = ts.apply_tail_decoder(pop, pb_ref, genes, held, seed=1)
        _, rep_pro = ts.apply_tail_decoder_profile(pop, pb_ref, genes, held, 1, src)
        own_sel, pro_sel = seen[-2], seen[-1]
        assert np.array_equal(own_sel, pro_sel) and np.array_equal(pro_sel, _expected_selection(pop, pb_ref, genes,
                                                                                                held))
        assert rep_own["n_infeasible_in_bottom"] == rep_pro["n_infeasible_in_bottom"] == len(pro_sel) > 0
        assert set(rep_pro["lifted"]) <= {genes[i] for i in pro_sel}
        assert rep_pro["n_lifted"] + rep_pro["n_skipped"] == rep_pro["n_infeasible_in_bottom"]


# ---------------------------------------------------------------------------------------------------------------------
# 2. A sparse gene the own-count decoder skips is lifted by the profile kernel, to flat

def _profile_case(n=6000, seed=21, mu0=0.05):
    """g0 sparse (Poisson mu0 per cell, so zero own counts in about 95% of cells at 0.05) in a source stage of n cells
    in five clusters over 200 genes (library about 1900), so that the 30 PCs follow the clusters and not g0; the cells
    are the source's own rows, and the kNN profile spreads g0's mass to many cells without counts of it."""
    rng = np.random.default_rng(seed)
    g, k = 200, 5
    genes = [f"g{i}" for i in range(g)]
    mu = rng.gamma(2.0, 5.0, (k, g))[rng.integers(0, k, n)]
    mu[:, 0] = mu0
    C = sparse.csr_matrix(rng.poisson(mu))
    pop = ad.AnnData(X=log1p_cp10k(C), layers={"counts": C.astype(np.int32)},
                     obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    pop.var_names = genes
    pop, src = _own_source(pop)
    pb = np.asarray(pop.X.mean(0)).ravel().astype(np.float64)
    return pop, src, genes, pb


def _slopes(pop, src, gene):
    """The own-count and the profile kernel's first-order slopes of the gene's mean log1p CP10k at f = 1."""
    C = ts.counts_of(pop)
    j = list(pop.var_names).index(gene)
    c = C[:, j].toarray().ravel().astype(float)
    a = 1e4 / np.asarray(C.sum(1)).ravel()
    m = ts.topup_profile(src, pop)[:, j].toarray().ravel()
    step = np.log1p(a * (c + 1)) - np.log1p(a * c)
    return float((c * step).mean()), float((m * step).mean())


def test_a_sparse_gene_the_own_counts_skip_is_lifted_by_the_profile_to_flat():
    pop, src, genes, pb = _profile_case()
    s_own, s_pro = _slopes(pop, src, "g0")
    c0 = ts.counts_of(pop)[:, 0].toarray().ravel()
    m0 = ts.topup_profile(src, pop)[:, 0].toarray().ravel()
    assert (c0 == 0).mean() > 0.9 and (m0[c0 == 0] > 0).mean() > 0.3  # zero own counts in most, profile mass in many
    fall = 0.5 * np.log(ts.TAIL_F_MAX) * s_pro  # half the profile's reach
    assert fall > np.log(ts.TAIL_F_MAX) * s_own  # beyond the own-count kernel's reach at the same bound
    pb_ref = pb.copy()
    pb_ref[0] += fall
    assert pb_ref[0] < ts.TAIL_MIN_PB  # still infeasible
    _, rep_own = ts.apply_tail_decoder(pop, pb_ref, genes, (), seed=[3, 15])
    assert rep_own["n_infeasible_in_bottom"] == 1 and rep_own["n_lifted"] == 0 and rep_own["n_skipped"] == 1
    out, rep = ts.apply_tail_decoder_profile(pop, pb_ref, genes, (), [3, 15], src)
    assert rep["lifted"] == ["g0"] and rep["n_lifted"] == 1 and rep["n_skipped"] == 0
    assert 1.0 < rep["factors"]["g0"] <= ts.TAIL_F_MAX and rep["kernel"] == "profile"
    new0 = ts.counts_of(out)[:, 0].toarray().ravel()
    assert ((c0 == 0) & (new0 > 0)).sum() > 50  # zeros raised: what the own-count kernel cannot do
    dp0 = float(np.asarray(out.X[:, 0].mean(0)).ravel()[0] - pb_ref[0])
    # Tolerance: the sampling noise of the realised mean (4 standard errors of the added log1p mass, which is about
    # Poisson in each cell) plus the second-order term the first-order solve leaves (the rise is at most fall, and
    # E log1p(a (c + P(lam))) departs from its tangent by a relative O(lam), lam = (f - 1) m < 0.1 here).
    lam = (rep["factors"]["g0"] - 1.0) * m0
    se = np.sqrt((lam * np.log1p(1e4 / np.asarray(ts.counts_of(pop).sum(1)).ravel()) ** 2).sum()) / pop.n_obs
    tol = 4 * se + fall * lam.max()
    assert abs(dp0) < tol < fall  # flat, within a tolerance that still tells a lift from no lift (dp = -fall)


# ---------------------------------------------------------------------------------------------------------------------
# 3. Skips: no profile mass, or a need beyond the bound

def test_no_profile_mass_and_out_of_reach_are_skipped_and_counted():
    pop, src, genes, pb = _profile_case(n=2000, mu0=0.01)  # sparser, so three times the reach stays infeasible
    C = ts.counts_of(pop).tolil()
    C[:, 1] = 0  # g1 has no counts in the source, so no profile mass anywhere: slope zero
    C = sparse.csr_matrix(C)
    pop.layers["counts"] = C.astype(np.int32)
    pop.X = log1p_cp10k(C)
    src = pop.copy()
    pb = np.asarray(pop.X.mean(0)).ravel().astype(np.float64)
    s_own, s_pro = _slopes(pop, src, "g0")
    pb_ref = pb.copy()
    pb_ref[1] = 0.01  # g1 fell (from 0.01 to 0) and is infeasible
    pb_ref[0] += 3.0 * np.log(ts.TAIL_F_MAX) * s_pro  # g0 needs three times the bound's reach
    assert pb_ref[0] < ts.TAIL_MIN_PB
    out, rep = ts.apply_tail_decoder_profile(pop, pb_ref, genes, (), 1, src)
    assert rep["n_infeasible_in_bottom"] == 2 and rep["n_lifted"] == 0 and rep["n_skipped"] == 2
    assert out is pop and np.isnan(rep["max_factor"]) and rep["factors"] == {}
    pb_ref[0] = pb[0] + 0.5 * np.log(ts.TAIL_F_MAX) * s_pro  # g0 back within reach: g1 alone is skipped
    _, rep = ts.apply_tail_decoder_profile(pop, pb_ref, genes, (), 1, src)
    assert rep["lifted"] == ["g0"] and rep["n_skipped"] == 1 and rep["n_infeasible_in_bottom"] == 2
    assert 1.0 <= rep["median_factor"] <= rep["max_factor"] <= ts.TAIL_F_MAX


# ---------------------------------------------------------------------------------------------------------------------
# 4 and 5. Depth kept; held genes untouched; every other free gene moves only by the depth constant

def test_depth_is_kept_held_genes_are_untouched_and_others_only_thinned(monkeypatch):
    pop, src, genes, pb = _profile_case()
    s_own, s_pro = _slopes(pop, src, "g0")
    pb_ref = pb.copy()
    pb_ref[0] += 0.5 * np.log(ts.TAIL_F_MAX) * s_pro
    held = ("g7",)
    scales = []
    real = ts.apply_factors_profile

    def spy(*args, **kw):
        res = real(*args, **kw)
        scales.append(res[1])
        return res

    monkeypatch.setattr(ts, "apply_factors_profile", spy)
    out, rep = ts.apply_tail_decoder_profile(pop, pb_ref, genes, held, [2, 15], src)
    assert rep["lifted"] == ["g0"]
    C0, C1 = ts.counts_of(pop), ts.counts_of(out)
    lib0, lib1 = np.asarray(C0.sum(1)).ravel(), np.asarray(C1.sum(1)).ravel()
    assert abs(rep["depth_ratio_mean"] - 1) < 0.01 and abs(np.mean(lib1 / lib0) - 1) < 0.01
    assert rep["expected_depth_max_error"] < 1e-9 and rep["no_cell_was_restored"]
    assert (C1[:, 7] != C0[:, 7]).nnz == 0  # held: untouched
    (scale,) = scales
    assert np.all(scale <= 1.0) and np.all(scale > 0.99)
    others = [j for j in range(len(genes)) if j not in (0, 7)]
    a, b = C0[:, others].toarray(), C1[:, others].toarray()
    assert np.all(b <= a)  # Binomial(c, s): never raised, never a zero raised
    expected = (scale[:, None] * a).sum()
    assert abs(b.sum() / expected - 1) < 0.001  # their mass is the depth constant's thinning, nothing more
    assert list(out.obs_names) == list(pop.obs_names) and out.n_obs == pop.n_obs


# ---------------------------------------------------------------------------------------------------------------------
# 6. Nothing to lift: the input, bit for bit

def test_nothing_to_lift_returns_the_population():
    pop, src, genes, pb = _profile_case(n=1000)
    out, rep = ts.apply_tail_decoder_profile(pop, pb - 0.01, genes, (), 1, src)  # every gene rose
    assert out is pop and rep["n_infeasible_in_bottom"] == 0 and rep["n_lifted"] == 0
    out, rep = ts.apply_tail_decoder_profile(pop, pb + 0.3, genes, (), 1, src)  # all fell, all feasible
    assert out is pop and rep["n_infeasible_in_bottom"] == 0 and np.isnan(rep["max_factor"])


def test_profile_arm_is_b0_bit_for_bit_when_nothing_is_infeasible(monkeypatch, capsys):
    ref, start = _pair()
    a = _build("b0", ref, start)
    capsys.readouterr()
    b = _build("b0_tail_profile", ref, start)  # synthetic genes all well expressed: none is infeasible
    assert _sha(np.asarray(a.X, np.float32)) == _sha(np.asarray(b.X, np.float32))
    assert _sha(np.asarray(b.X, np.float32)) == "23fab0fea7d206ed9d5f35050274d70e9c23ef5941888e869e8faf2a17a82508"
    lines = [x for x in capsys.readouterr().out.splitlines() if "check:" in x]
    assert lines == [f"b0_tail_profile check: task={TASK.key} seed={harness.BUILD_SEED_OFFSET + 1} "
                     "infeasible_in_bottom=0 lifted=0 skipped=0 max_factor=nan median_factor=nan"]
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 0.0)  # strength zero
    c = _build("b0_tail_profile", ref, start)
    assert _sha(np.asarray(c.X, np.float32)) == _sha(np.asarray(a.X, np.float32))


# ---------------------------------------------------------------------------------------------------------------------
# 7. Seeds

def test_same_seed_same_output_other_seed_other_draw():
    pop, src, genes, pb = _profile_case(n=2000)
    _, s_pro = _slopes(pop, src, "g0")
    pb_ref = pb.copy()
    pb_ref[0] += 0.5 * np.log(ts.TAIL_F_MAX) * s_pro
    o1, r1 = ts.apply_tail_decoder_profile(pop, pb_ref, genes, (), [4, 15], src)
    o2, r2 = ts.apply_tail_decoder_profile(pop, pb_ref, genes, (), [4, 15], src)
    o3, r3 = ts.apply_tail_decoder_profile(pop, pb_ref, genes, (), [5, 15], src)
    assert r1["n_lifted"] == 1 and r1["factors"] == r2["factors"] == r3["factors"]
    assert (ts.counts_of(o1) != ts.counts_of(o2)).nnz == 0
    assert (ts.counts_of(o1) != ts.counts_of(o3)).nnz > 0


# ---------------------------------------------------------------------------------------------------------------------
# 8. The check line; registration; the harness wiring

def test_the_check_line_holds_only_the_allowed_fields():
    rep = {"n_infeasible_in_bottom": 31, "n_lifted": 4, "n_skipped": 27, "max_factor": 1.98765,
           "median_factor": 1.5, "lifted": ["g0"], "factors": {"g0": 1.9}, "depth_ratio_mean": 0.9999}
    line = harness.tail_check_line("T1/E9.5->E10.5", 1001, rep, name="b0_tail_profile")
    assert line == ("b0_tail_profile check: task=T1/E9.5->E10.5 seed=1001 infeasible_in_bottom=31 lifted=4 "
                    "skipped=27 max_factor=1.99 median_factor=1.5")
    keys = [t.split("=")[0] for t in line.split()[2:]]
    assert keys == ["task", "seed", "infeasible_in_bottom", "lifted", "skipped", "max_factor", "median_factor"]
    assert line.split(" check: ")[1] == harness.tail_check_line("T1/E9.5->E10.5", 1001, rep).split(" check: ")[1]
    assert "g0" not in line


def test_registered_like_b0_tail():
    tasks = harness.load_tasks()["T1:test"]
    c, ref = harness.CANDIDATES["b0_tail_profile"], harness.CANDIDATES["b0_tail"]
    assert c.tasks == ref.tasks and c.derived == ref.derived and c.externals == () and c.tuned_on == frozenset()
    assert [t.key for t in tasks if c.applies(t)] == [t.key for t in tasks if ref.applies(t)]
    for k in ("damping", "scope", "bottom", "min_pb", "target_dp", "f_max"):
        assert c.params[k] == ref.params[k]
    assert (c.params["n_pcs"], c.params["k"]) == (ts.TOPUP_N_PCS, ts.TOPUP_K)
    assert harness.TAIL_PROFILE_STREAM not in (harness.TAIL_STREAM, 0, 1, 2, 3, 4, 5)


def test_profile_arm_takes_b0s_cells_and_the_reference_as_source(monkeypatch, capsys):
    ref, start = _pair()
    a = _build("b0", ref, start)
    capsys.readouterr()
    seen = []
    real = ts.apply_tail_decoder_profile

    def spy(pop, pb_ref, genes, held, seed, source):
        out, rep = real(pop, pb_ref, genes, held, seed, source)
        seen.append((pop, out, rep, held, seed, source))
        return out, rep

    monkeypatch.setattr(ts, "apply_tail_decoder_profile", spy)
    monkeypatch.setattr(ts, "TAIL_MIN_PB", 10.0)  # make the bottom tail count on well-expressed synthetic genes
    b = _build("b0_tail_profile", ref, start)
    (pop, out, rep, held, seed, source), = seen
    assert list(held) == ["g49"] and seed == [harness.BUILD_SEED_OFFSET + 1, harness.TAIL_PROFILE_STREAM]
    assert list(source.obs_names) == list(ref.obs_names) and source.n_vars == ref.n_vars
    assert _sha(np.asarray(pop.X.toarray(), np.float32)) == _sha(np.asarray(a.X, np.float32))  # b0 before it
    assert rep["n_infeasible_in_bottom"] > 0 and a.X.shape == b.X.shape
    if rep["n_lifted"]:
        lib0, lib1 = (np.asarray(ts.counts_of(x).sum(1)).ravel() for x in (pop, out))
        assert abs(lib1.mean() / lib0.mean() - 1) < 0.01
        assert (ts.counts_of(out)[:, 49] != ts.counts_of(pop)[:, 49]).nnz == 0
    lines = [x for x in capsys.readouterr().out.splitlines() if "check:" in x]
    assert lines == [harness.tail_check_line(TASK.key, seed[0], rep, name="b0_tail_profile")]


@pytest.mark.parametrize("name", ["b0_tail", "b0_topup"])
def test_the_t1b_arms_are_unchanged(name, capsys):
    from test_t1_batch_b import B0_TOPUP_SHA

    ref, start = _pair()
    want = {"b0_tail": "23fab0fea7d206ed9d5f35050274d70e9c23ef5941888e869e8faf2a17a82508", "b0_topup": B0_TOPUP_SHA}
    assert _sha(np.asarray(_build(name, ref, start).X, np.float32)) == want[name]
