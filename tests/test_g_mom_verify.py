"""The g_mom verification (tries T20 to T24 and diagnostic N21 of experiments/harness/TRIES.md) on synthetic data: the
placebo and no-development factors, the pseudobulk-held arm, their harness arms, the pair fit shared by N21 and the test
build, the controls' readout, the composition arms at the truth's own chamber mix, N21, the --g-mom option of the T1
test build, and a pair of synthetic worlds (with and without within-type momentum) that the controls must separate."""

from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from test_hidden_counts import T1_GENES, T1_HELD, UNION, _reference_file, bt, t1_env  # noqa: F401
from test_t1_scope import GENES, HELD, _profiles
from test_t1_state_transport import _gmom_ctx, _t1_task, _toy_population, states_stage

from vec import harness, t1_scope
from vec import t1_state_transport as ts
from vec.boards import REPO_ROOT, Board
from vec.evaluate import Stage
from vec.harness import CANDIDATES
from vec.prep import log1p_cp10k

SCRIPT = REPO_ROOT / "experiments" / "t1_scope" / "g_mom_verify.py"


def _small_inputs(seed=7):
    """E9.5- and E8.5-like stages with within-type structure, as `_t1_inputs` makes them but small (one Sinkhorn fit
    per category costs little), and an E10.5-like target with a real change, every label classified by the config."""
    ref = states_stage(seed, {"V-CM": 420, "Endocardium": 170, "Blood": 60, "Foregut": 100}, later=True)
    start = states_stage(seed + 1, {"V-CM": 300, "Endothelium": 130, "Blood": 50, "Foregut": 90})
    target = states_stage(seed + 2, {"V-CM": 500, "Endocardium": 220, "Blood": 80, "Foregut": 120}, later=True,
                          target=True)
    return ref, start, target


@pytest.fixture(scope="module")
def world():
    """One small set of stages; their fold factors are fitted once (harness._G_MOM_PREP) and shared by the arms."""
    harness._G_MOM_PREP.clear()
    return _small_inputs()


@pytest.fixture
def vm(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("g_mom_verify_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.provenance, "git_state", lambda *a, **k: ("abc123", []))
    return mod


def _factors(k=6, genes=12, seed=0, sizes=None):
    rng = np.random.default_rng(seed)
    sizes = np.array(sizes if sizes is not None else [20, 30, 40, 10, 50, 25][:k], dtype=np.int64)
    logf = rng.normal(0, 0.4, size=(k, genes))
    logf[np.abs(logf) < 0.1] = 0.0
    index = {f"c{i}": int(i % k) for i in range(int(sizes.sum()))}
    cf = ts.CategoryFactors("cm", index, sizes, logf, record={"finish": "x"})
    other = ts.CategoryFactors("endocardium", {"e0": 0, "e1": 1}, np.array([3, 4]), rng.normal(0, 0.3, size=(2, genes)))
    return ts.StateFactors(genes=[f"g{i}" for i in range(genes)], metacells_on="later",
                           categories={"cm": cf, "endocardium": other}, skipped={}, record={"finish": "momentum"})


# ------------------------------------------------------------------------------------------------ the placebo factors
def test_flip_factors_reverses_every_sign_and_leaves_the_input_alone():
    sf = _factors()
    flipped = ts.flip_factors(sf)
    for c, cf in sf.categories.items():
        assert np.array_equal(flipped.categories[c].logf, -cf.logf)
        assert np.array_equal(sf.categories[c].logf != 0, flipped.categories[c].logf != 0)  # zeros stay zero
    assert ts.flip_factors(flipped).categories["cm"].logf.tolist() == sf.categories["cm"].logf.tolist()
    assert "sign reversed" in flipped.categories["cm"].record["finish"] and "T20" in flipped.record["finish"]
    assert not np.array_equal(flipped.categories["cm"].logf, sf.categories["cm"].logf)  # the input is not modified


def test_permute_factors_is_a_fixed_derangement_of_each_categorys_rows():
    sf = _factors()
    a = ts.permute_factors(sf, 20261110)
    b = ts.permute_factors(sf, 20261110)
    c = ts.permute_factors(sf, 7)
    for cat, cf in sf.categories.items():
        pa = np.array(a.record["permutation"]["per_category"][cat])
        assert np.array_equal(a.categories[cat].logf, cf.logf[pa])
        assert np.array_equal(b.categories[cat].logf, a.categories[cat].logf)  # the seed fixes the permutation
        assert not np.any(pa == np.arange(cf.n_metacells))  # a derangement: no metacell keeps its own vector
        # every row (the gene programme of a metacell) survives intact, only the metacells it sits on change
        before = {tuple(r) for r in np.round(cf.logf, 12)}
        after = {tuple(r) for r in np.round(a.categories[cat].logf, 12)}
        assert before == after and sorted(np.abs(cf.logf).ravel()) == pytest.approx(
            sorted(np.abs(a.categories[cat].logf).ravel()))
    assert not np.array_equal(c.categories["cm"].logf, a.categories["cm"].logf)
    assert a.record["permutation"]["seed"] == 20261110 and "T21" in a.record["finish"]
    # the permutation depends on the categories' sizes and the seed, not on the dictionary order of the categories
    reordered = ts.StateFactors(genes=sf.genes, metacells_on=sf.metacells_on,
                                categories=dict(reversed(list(sf.categories.items()))), skipped={}, record={})
    assert np.array_equal(ts.permute_factors(reordered, 20261110).categories["cm"].logf, a.categories["cm"].logf)


def test_permute_factors_leaves_a_single_metacell_as_it_is():
    cf = ts.CategoryFactors("cm", {"x": 0}, np.array([5]), np.array([[0.2, 0.0]]))
    sf = ts.StateFactors(genes=["a", "b"], metacells_on="later", categories={"cm": cf}, skipped={}, record={})
    assert ts.permute_factors(sf, 1).categories["cm"].logf.tolist() == [[0.2, 0.0]]


def test_permute_factors_changes_nothing_for_identical_rows_only_the_assignment():
    """On a matrix whose rows differ only in which metacell holds them, the permuted applied factors move the same cells
    by another metacell's vector: the rows are re-assigned, so the sum over metacells of each gene's factor is the
    same whenever the metacells weigh the same."""
    logf = np.array([[0.5, -0.5, 0.0], [-0.5, 0.5, 0.0], [0.2, 0.0, -0.2], [0.0, -0.2, 0.2], [0.1, 0.1, -0.2]])
    cf = ts.CategoryFactors("cm", {f"c{i}": i for i in range(5)}, np.array([10] * 5), logf)
    sf = ts.StateFactors(genes=["a", "b", "c"], metacells_on="later", categories={"cm": cf}, skipped={}, record={})
    out = ts.permute_factors(sf, 3).categories["cm"]
    assert np.allclose(out.logf.sum(0), logf.sum(0))  # equal weights: the permuted per-gene mean stays exactly 0
    assert not np.array_equal(out.logf, logf)


# ------------------------------------------------------------------------------------------ the no-development factors
def _raw_fit(prefix, seed, n_meta=6, genes=40):
    """A raw fit on its own cells (names prefix_i), the sizes unequal so the cell weighting matters."""
    rng = np.random.default_rng(seed)
    sizes = rng.integers(5, 40, size=n_meta)
    names = [f"{prefix}{i}" for i in range(int(sizes.sum()))]
    index = dict(zip(names, np.repeat(np.arange(n_meta), sizes).tolist(), strict=True))
    cf = ts.CategoryFactors("cm", index, sizes, rng.normal(0.0, 0.3, size=(n_meta, genes)))
    return ts.StateFactors(genes=[f"g{i}" for i in range(genes)], metacells_on="later", categories={"cm": cf},
                           skipped={}, record={})


def test_no_development_finish_is_cross_fitted_mean_free_per_direction_and_matched_to_g_moms_size():
    held = ["g0", "g1"]
    first, second = _raw_fit("a", 1), _raw_fit("b", 2)
    target = ts.finish_momentum(_raw_fit("x", 3), 0.5)
    out = ts.finish_no_development(first, second, target, held, 0.5)
    cf = out.categories["cm"]
    # the two directions' metacells are stacked: every cell of both halves has a metacell, the second's shifted
    ka = first.categories["cm"].n_metacells
    assert cf.n_metacells == ka + second.categories["cm"].n_metacells and len(cf.index) == (
        len(first.categories["cm"].index) + len(second.categories["cm"].index))
    assert all(m >= ka for n, m in cf.index.items() if n.startswith("b"))
    free = np.array([g not in held for g in out.genes])
    # the size is g_mom's per category: the cell-weighted root mean square over the metacells and the free genes
    want = ts._factor_rms(target.categories["cm"].logf, target.categories["cm"].sizes, free)
    assert ts._factor_rms(cf.logf, cf.sizes, free) == pytest.approx(want, rel=0.02)
    rec = cf.record["no_development"]
    assert rec["rms_target"] == pytest.approx(want) and rec["relative_error"] < 0.02 and rec["damping"] > 0
    # each direction's mean was removed before the stacking, and the finish is the zero-and-clip of g_mom
    assert np.abs(cf.logf).max() <= ts.CLIP_AT + 1e-12
    assert np.all((np.abs(cf.logf) >= ts.ZERO_BELOW - 1e-12) | (cf.logf == 0))
    # what is matched is the size, not the direction: another target size scales the factors, the zero set aside
    bigger = ts.finish_momentum(_raw_fit("x", 3), 1.5)
    scaled = ts.finish_no_development(first, second, bigger, held, 0.5).categories["cm"]
    assert ts._factor_rms(scaled.logf, scaled.sizes, free) > ts._factor_rms(cf.logf, cf.sizes, free) * 1.3


def test_no_development_finish_refuses_shared_cells_and_skips_a_category_it_cannot_match():
    first = _raw_fit("a", 1)
    with pytest.raises(ValueError, match="share cells"):
        ts.finish_no_development(first, first, ts.finish_momentum(_raw_fit("x", 3), 0.5), [], 0.5)
    second = _raw_fit("b", 2)
    # a category g_mom has no factors for gets none; a category fitted in one direction only is skipped and recorded
    target = ts.finish_momentum(_raw_fit("x", 3), 0.5)
    only_one = ts.StateFactors(genes=second.genes, metacells_on="later", categories={}, skipped={"cm": "too few"},
                               record={})
    out = ts.finish_no_development(first, only_one, target, [], 0.5)
    assert not out.categories and "cm" in out.skipped


# ------------------------------------------------------------------------------------------------ the harness arms
def test_g_mom_arm_with_the_identity_is_byte_for_byte_g_mom(world):
    """g_mom itself is unchanged, and the arm repeats its steps: with an identity transform the arm's file is g_mom."""
    ref, start, _ = world
    ctx = _gmom_ctx(ref, start)
    gm = harness._g_mom(ctx)
    arm = harness._g_mom_arm(ctx, lambda f: f)
    assert np.array_equal(arm.X, gm.X) and arm.X.dtype == gm.X.dtype


def test_the_placebo_arms_share_g_moms_base_and_draws_and_differ_only_in_the_factors(world):
    ref, start, _ = world
    ctx = _gmom_ctx(ref, start)
    b0, gm = harness._b0_t1(harness.B0_DAMPING)(ctx), harness._g_mom(ctx)
    flip = CANDIDATES["g_mom_flip"].build(ctx)
    perm = CANDIDATES["g_mom_perm"].build(ctx)
    assert not np.array_equal(flip.X, gm.X) and not np.array_equal(perm.X, gm.X) and not np.array_equal(flip.X, perm.X)
    # the same population: cells of R1's held labels (Blood) keep B0's counts in every arm
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E9.5_RNA")
    pop = ts.b0_population(harness._shallow(ref), harness._shallow(start), shares, 600, ctx.seed, list(HELD),
                           ref_source="E9.5_RNA", start_source="E8.5_RNA", damping=harness.B0_DAMPING)
    blood = pop.obs["celltype"].astype(str).to_numpy() == "Blood"
    assert blood.any()
    for arm in (gm, flip, perm):
        assert np.array_equal(arm.X[blood], b0.X[blood])
    # double sign reversal is g_mom again: the flip is exactly the negation of the factors
    twice = harness._g_mom_arm(ctx, lambda f: ts.flip_factors(ts.flip_factors(f)))
    assert np.array_equal(twice.X, gm.X)


def test_the_flip_uses_the_same_stream_and_the_permutation_one_fixed_seed(world):
    ref, start, _ = world
    ctx = _gmom_ctx(ref, start)
    one = CANDIDATES["g_mom_perm"].build(ctx)
    two = CANDIDATES["g_mom_perm"].build(ctx)  # the permutation does not depend on the build seed's stream
    assert np.array_equal(one.X, two.X) and harness.G_MOM_PERM_SEED == 20261110
    other = harness._g_mom_arm(_gmom_ctx(ref, start, seed=1000001),
                               lambda f: ts.permute_factors(f, harness.G_MOM_PERM_SEED))
    assert not np.array_equal(one.X, other.X)  # another build seed: another draw of cells and counts


def test_the_no_development_arm_uses_g_moms_draws_and_stream_and_factors_for_every_reference_cell(world):
    """T23: the same base and draws as g_mom (R1's held labels keep B0's counts), cross-fitted so that every cell of
    g_mom's fitted categories has a factor from the direction in which its own half carries the metacells, at g_mom's
    factor size per category, and fitted once per fold (the object stays the same for another seed)."""
    ref, start, _ = world
    ctx = _gmom_ctx(ref, start)
    b0, gm = harness._b0_t1(harness.B0_DAMPING)(ctx), harness._g_mom(ctx)
    null = CANDIDATES["g_mom_null"].build(ctx)
    assert not np.array_equal(null.X, gm.X)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E9.5_RNA")
    pop = ts.b0_population(harness._shallow(ref), harness._shallow(start), shares, 600, ctx.seed, list(HELD),
                           ref_source="E9.5_RNA", start_source="E8.5_RNA", damping=harness.B0_DAMPING)
    blood = pop.obs["celltype"].astype(str).to_numpy() == "Blood"
    assert blood.any() and np.array_equal(null.X[blood], b0.X[blood])
    finished = harness._G_MOM_PREP["fold"]["factors"]
    factors = harness._G_MOM_NULL_PREP["fold"]["factors"]
    assert harness._G_MOM_NULL_PREP["fold"]["ref"] is ctx.inputs["ref"]
    assert set(factors.categories) <= set(finished.categories)  # g_mom's categories, and the control edits only them
    t1_scope_cats = harness._shallow(ref)
    t1_scope.classify(t1_scope_cats, source="E9.5_RNA")
    cats = t1_scope_cats.obs["category"].astype(str).to_numpy()
    for c, cf in factors.categories.items():
        # the two halves' metacells together cover the category's own cells (the cells g_mom's fit covers)
        assert set(cf.index) == set(finished.categories[c].index) or cf.index.keys() >= set(
            np.asarray(ref.obs_names)[cats == c].tolist())
    # the size is g_mom's, the record says so, and a second build (another seed) reuses the fold's fit
    free = np.array([g not in set(HELD) for g in factors.genes])
    for c, cf in factors.categories.items():
        goal = finished.categories[c]
        assert ts._factor_rms(cf.logf, cf.sizes, free) == pytest.approx(ts._factor_rms(goal.logf, goal.sizes, free),
                                                                        rel=0.1)
        assert cf.record["no_development"]["rms_target"] > 0
    CANDIDATES["g_mom_null"].build(_gmom_ctx(ref, start, seed=1000002))
    assert harness._G_MOM_NULL_PREP["fold"]["factors"] is factors
    assert CANDIDATES["g_mom_null"].params["split_seed"] == harness.G_MOM_NULL_SEED == 20261110


def test_the_halves_are_disjoint_and_cover_the_reference(world):
    ref, _, _ = world
    stage = harness._shallow(ref)
    t1_scope.classify(stage, source="E9.5_RNA")
    a = ts.category_halves(stage, np.random.default_rng(harness.G_MOM_NULL_SEED))
    b = ts.category_halves(stage, np.random.default_rng(harness.G_MOM_NULL_SEED))
    assert np.array_equal(a, b)  # the seed fixes the halves
    cats = stage.obs["category"].astype(str).to_numpy()
    for c in set(cats):
        assert a[cats == c].sum() == (cats == c).sum() // 2


def test_the_held_arm_is_g_mom_through_the_rescale_on_the_same_factors_and_stream(world, monkeypatch):
    """T24: the same factors, draws and stream as g_mom through rescale=True (T9's construction): the category's mean
    CP10k is held per gene in expectation, and R1's held labels keep B0's counts."""
    ref, start, _ = world
    ctx = _gmom_ctx(ref, start)
    b0, gm = harness._b0_t1(harness.B0_DAMPING)(ctx), harness._g_mom(ctx)
    calls = []
    real = ts.apply_state_factors

    def spy(pop, sf, held_flat, seed, **kw):
        calls.append((sf, list(seed), kw["rescale"]))
        return real(pop, sf, held_flat, seed, **kw)

    monkeypatch.setattr(ts, "apply_state_factors", spy)
    held_arm = CANDIDATES["g_mom_held"].build(ctx)
    monkeypatch.setattr(ts, "apply_state_factors", real)
    assert len(calls) == 1 and calls[0][1] == [ctx.seed, 5] and calls[0][2] is True
    assert calls[0][0] is harness._G_MOM_PREP["fold"]["factors"]  # g_mom's own finished factors, untransformed
    assert not np.array_equal(held_arm.X, gm.X) and not np.array_equal(held_arm.X, b0.X)
    twin = harness._g_mom_arm(ctx, lambda f: f, rescale=True)
    assert np.array_equal(twin.X, held_arm.X)  # the stream is fixed: the build is repeatable
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E9.5_RNA")
    pop = ts.b0_population(harness._shallow(ref), harness._shallow(start), shares, 600, ctx.seed, list(HELD),
                           ref_source="E9.5_RNA", start_source="E8.5_RNA", damping=harness.B0_DAMPING)
    blood = pop.obs["celltype"].astype(str).to_numpy() == "Blood"
    assert blood.any() and np.array_equal(held_arm.X[blood], b0.X[blood])


def test_the_placebo_arms_build_inside_the_read_guard(monkeypatch, tmp_path):
    """The harness's own run of g_mom and its controls on a T1 pseudo-task (one seed): nothing a candidate must not read
    is read, and every arm is scored against the shared anchors."""
    ref, start, tgt = _small_inputs(seed=47)
    harness._G_MOM_PREP.clear()
    task = _t1_task()
    board = Board(key="pseudo:t1gm", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=400,
                  truth_cells=60, ref_cells=40)
    stages = {"target": Stage(X=tgt.X.toarray()), "ref": Stage(X=ref.X.toarray())}

    def read(path, genes):
        s = stages["target"] if str(path).endswith(task.target) else stages["ref"]
        return ad.AnnData(X=np.asarray(s.X), var=pd.DataFrame(index=GENES))

    by_file = {"T1/E9.5_RNA.h5ad": ref, "T1/E8.5_RNA.h5ad": start}
    monkeypatch.setattr(harness, "read_stage", read)
    monkeypatch.setattr(harness, "load_stage", lambda path, **kw: by_file[str(path).split(str(tmp_path) + "/")[1]])
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    monkeypatch.setattr(harness, "derived_ok", lambda name, t: True)
    names = ["b0", "g_mom", "g_mom_flip", "g_mom_perm", "g_mom_null", "g_mom_held"]
    res = harness.run_task(task, [CANDIDATES[n] for n in names], seeds=(0,), data_dir=tmp_path,
                           boards={"T1:val": replace(board, key="T1:val")}, log=None)
    row = res["rows"][0]
    assert res["candidates"] == ["floor", *names] and row["scores"]["floor"] == 50.0
    for name in names:
        assert np.isfinite(row["scores"][name]) and row["n_cells"][name] == 400
        assert "variance_ratio" in row["metrics"][name] and "de_score" in row["skills"][name]  # the readout's inputs
    for name in ("g_mom", "g_mom_flip", "g_mom_perm", "g_mom_null", "g_mom_held"):
        assert not any("results" in p or p.startswith(("data/", "runs/")) for p in row["reads"][name])


def test_the_control_candidates_apply_where_g_mom_does():
    tasks = {t.key: t for ts_ in harness.load_tasks().values() for t in ts_}
    for name in ("g_mom_flip", "g_mom_perm", "g_mom_null", "g_mom_held"):
        c = CANDIDATES[name]
        assert c.applies(tasks["T1/E9.5->E10.5"]) and c.applies(tasks["T1/E9.5->E10.5:chamber"])
        assert not c.applies(tasks["T1/E8.5->E9.5"]) and not c.applies(tasks["T2h/<=E9.5->E10.5"])
        assert c.params["scope"] == CANDIDATES["g_mom"].params["scope"] and c.derived == ("t1_held_flat",)
    assert CANDIDATES["g_mom_perm"].params["permutation_seed"] == harness.G_MOM_PERM_SEED
    assert CANDIDATES["g_mom_held"].params["scope"] == "chamber"
    for name in ("g_mom_null", "g_mom_held"):
        assert CANDIDATES[name].derived == ("t1_held_flat",)


# ------------------------------------------------------------------------------------------------ two synthetic worlds
def _world_stage(seed, sizes, step, momentum, target=False):
    """A T1-like stage whose V-CM cells have two states (half of them with g30-g39 up 8x at every stage). With momentum,
    the second state also moves along its own programme g40-g44 by a factor 6 per stage (step 0, 1, 2 = start,
    reference, target), so the reference's own state contrast barely shows on those genes at the start and the true
    within-type change keeps going past the reference; without momentum nothing within a type changes between stages,
    only the whole-population shift the target carries (target=True, a broader version of `states_stage`'s)."""
    rng = np.random.default_rng(seed)
    prof = _profiles()
    blocks, labels = [], []
    for lab, n in sizes.items():
        mu = np.tile(prof[lab], (n, 1)) * rng.gamma(10.0, 0.1, size=(n, 1))
        if lab == "V-CM":
            second = np.arange(n) < n // 2
            mu[second, 30:40] *= 8.0
            if momentum:
                mu[second, 40:45] *= 6.0 ** step
        if target:  # the whole population's change, over enough genes that the scorer's DE sets are stable
            mu[:, :15] *= 4.0
            mu[:, 15:30] *= 0.25
        C = rng.poisson(mu)
        C[np.arange(n), rng.integers(0, C.shape[1], n)] += 1
        blocks.append(C)
        labels += [lab] * n
    C = sparse.csr_matrix(np.vstack(blocks).astype(np.int32))
    a = ad.AnnData(X=log1p_cp10k(C))
    a.var_names = GENES
    a.obs_names = [f"w{seed}_{i}" for i in range(a.n_obs)]
    a.obs["celltype"] = labels
    a.layers["counts"] = C
    return a


WORLD_ARMS = ["b0", "g_mom", "g_mom_flip", "g_mom_perm", "g_mom_null", "g_mom_held"]


def _run_world(monkeypatch, tmp_path, momentum, seed=71, seeds=(0, 1, 2)):
    """One world through the harness's own run on a T1 pseudo-task and the real scorer: per arm the mean over seeds of
    the score, the scores per seed, and the factors fitted (g_mom's and the no-development arm's)."""
    ref = _world_stage(seed, {"V-CM": 1800, "Endocardium": 600, "Blood": 200, "Foregut": 300}, 1, momentum)
    start = _world_stage(seed + 1, {"V-CM": 1400, "Endothelium": 500, "Blood": 180, "Foregut": 240}, 0, momentum)
    tgt = _world_stage(seed + 2, {"V-CM": 2200, "Endocardium": 800, "Blood": 240, "Foregut": 400}, 2, momentum,
                       target=True)
    harness._G_MOM_PREP.clear()
    harness._G_MOM_NULL_PREP.clear()
    task = _t1_task()
    board = Board(key="pseudo:t1gm", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=1200,
                  truth_cells=tgt.n_obs, ref_cells=ref.n_obs)

    def read(path, genes):
        stage = tgt if str(path).endswith(task.target) else ref
        return ad.AnnData(X=np.asarray(stage.X.todense()), var=pd.DataFrame(index=GENES))

    by_file = {"T1/E9.5_RNA.h5ad": ref, "T1/E8.5_RNA.h5ad": start}
    monkeypatch.setattr(harness, "read_stage", read)
    monkeypatch.setattr(harness, "load_stage", lambda path, **kw: by_file[str(path).split(str(tmp_path) + "/")[1]])
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    monkeypatch.setattr(harness, "derived_ok", lambda name, t: True)
    res = harness.run_task(task, [CANDIDATES[n] for n in WORLD_ARMS], seeds=seeds, data_dir=tmp_path,
                           boards={"T1:val": replace(board, key="T1:val")}, log=None)
    scores = {arm: [row["scores"][arm] for row in res["rows"]] for arm in WORLD_ARMS}
    return scores, res["rows"], (harness._G_MOM_PREP["fold"]["factors"], harness._G_MOM_NULL_PREP["fold"]["factors"])


def _gains(scores, arm, base):
    return float(np.mean(np.array(scores[arm]) - np.array(scores[base])))


def test_the_controls_separate_a_world_with_within_type_momentum_from_one_without(monkeypatch, tmp_path):
    """Through the real pipeline and scorer (toy size; no project data): without within-type momentum g_mom gains little
    over its no-development control, with real momentum (the second state keeps moving along a programme the reference's
    state contrast barely shows at the start) it gains clearly over it, the sign flip and the permutation; and the
    no-development arm, cross-fitted between halves of the reference, carries no g40-g44 momentum."""
    still, _, _ = _run_world(monkeypatch, tmp_path, momentum=False)
    moving, _, (moving_gm, moving_null) = _run_world(monkeypatch, tmp_path, momentum=True)
    # no within-type development: g_mom has nothing to gain over its no-development control, or over b0
    assert abs(_gains(still, "g_mom", "g_mom_null")) < 0.3 and _gains(still, "g_mom", "b0") < 0.3
    # real momentum: g_mom gains over b0 and over every control, and the gain survives the pseudobulk hold
    assert _gains(moving, "g_mom", "b0") > 0.5
    for control in ("g_mom_null", "g_mom_flip", "g_mom_perm"):
        assert _gains(moving, "g_mom", control) > 0.5, control
    assert _gains(moving, "g_mom_held", "b0") > 0.3
    # the no-development arm, which sees no development, cannot be the source of g_mom's gain in the moving world: on
    # the programme genes its factors are the smaller, though both arms have the same factor size per category
    cm_gm, cm_null = moving_gm.categories["cm"], moving_null.categories["cm"]
    idx = [GENES.index(f"g{i}") for i in range(40, 45)]
    w_gm, w_null = cm_gm.sizes / cm_gm.sizes.sum(), cm_null.sizes / cm_null.sizes.sum()
    size_gm = float((w_gm * np.abs(cm_gm.logf[:, idx]).mean(1)).sum())
    size_null = float((w_null * np.abs(cm_null.logf[:, idx]).mean(1)).sum())
    assert size_null < 0.8 * size_gm


# ------------------------------------------------------------------------------------------------ the pair fit
def test_pair_factors_are_the_folds_factors_and_align_the_start_stages_genes(world):
    ref, start, _ = world
    r0, s0 = harness._shallow(ref), harness._shallow(start)
    shares = harness.t1_chamber_shares(r0, "E9.5_RNA")
    ts.b0_population(r0, s0, shares, 400, 1000000, list(HELD), ref_source="E9.5_RNA", start_source="E8.5_RNA",
                     damping=harness.B0_DAMPING)
    ctx = _gmom_ctx(ref, start)
    fold = harness._g_mom_factors(ctx, r0, "E9.5_RNA", s0, "E8.5_RNA", shares)
    pair, raw = harness.g_mom_pair_factors(ref, "E9.5_RNA", start, "E8.5_RNA")
    assert set(pair.categories) == set(fold.categories) and raw.metacells_on == "later"
    for c, cf in fold.categories.items():
        assert np.array_equal(pair.categories[c].logf, cf.logf) and pair.categories[c].index == cf.index
    # the start stage's genes in another order are read by name: the same factors
    order = list(np.random.default_rng(1).permutation(start.var_names))
    shuffled = start[:, order].copy()
    pair_shuffled, _ = harness.g_mom_pair_factors(ref, "E9.5_RNA", shuffled, "E8.5_RNA")
    for c, cf in pair.categories.items():
        assert np.array_equal(pair_shuffled.categories[c].logf, cf.logf)
    # a gene the start lacks is refused (as counts_of reads a start without a reference gene)
    with pytest.raises(KeyError):
        harness.g_mom_pair_factors(ref, "E9.5_RNA", start[:, list(start.var_names[1:])].copy(), "E8.5_RNA")
    # the inputs are not modified
    assert "category" not in start.obs and "category" not in ref.obs


# ------------------------------------------------------------------------------------------------ the placebo readout
ARMS = ("floor", "b0", "g_mom", "g_mom_flip", "g_mom_perm", "g_mom_null", "g_mom_held")
PANEL = {"variance_ratio": 0.70, "mmd_u": 0.12, "variogram": 0.008, "energy_distance": 3.5, "pb_rel_err": 0.28,
         "de_score": -0.3, "de_direction": 0.17, "library_size_ratio": 1.0, "pseudobulk_pearson": 0.9}


def _run_rows(b0, g_mom, flip, perm, null=None, held=None, seeds=range(10), t11=None, vr=None, de_skill=None):
    """Per-seed rows of a harness run: scores, group skills, the metric panel and the skills per arm. null, held:
    default below g_mom by 1.0 and above b0 by 0.8 (the rule passes); t11: scores for b0 and g_mom to use instead; vr:
    arm -> variance_ratio; de_skill: arm -> (de_score, de_direction) skill."""
    rows = []
    vr = {"b0": 0.70, "g_mom": 0.712, "g_mom_flip": 0.688, "g_mom_perm": 0.711, "g_mom_null": 0.712,
          "g_mom_held": 0.711, "floor": 0.77, **(vr or {})}
    de_skill = {**{a: (0.5, 0.5) for a in ARMS}, **(de_skill or {})}
    for i, s in enumerate(seeds):
        scores = {"floor": 50.0, "b0": b0[i], "g_mom": g_mom[i], "g_mom_flip": flip[i], "g_mom_perm": perm[i],
                  "g_mom_null": g_mom[i] - 1.0 if null is None else null[i],
                  "g_mom_held": b0[i] + 0.8 if held is None else held[i]}
        if t11 is not None:
            scores["b0"], scores["g_mom"] = t11
        groups = {arm: {"de_recovery": 0.6, "direction": 0.6, "distribution": 0.7 + (v - 50) / 100,
                        "covariance": 0.7 + (v - 50) / 200} for arm, v in scores.items()}
        metrics = {arm: {**PANEL, "variance_ratio": vr[arm]} for arm in ARMS}
        skills = {arm: {"de_score": de_skill[arm][0], "de_direction": de_skill[arm][1]} for arm in ARMS}
        rows.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "scores": scores, "groups": groups,
                     "metrics": metrics, "skills": skills})
    return rows


def _write_run(tmp_path, folder, rows_by_cut, sha):
    d = tmp_path / folder
    d.mkdir(parents=True)
    for key, name in (("T1/E9.5->E10.5", "T1_E9.5_E10.5.json"),
                      ("T1/E9.5->E10.5:chamber", "T1_E9.5_E10.5_chamber.json")):
        (d / name).write_text(json.dumps(
            {"task": {"key": key, "test_board": "T1:test"},
             "candidates": list(ARMS), "rows": rows_by_cut[key],
             "protocol": "bake", "build": "main", "anchors": "computed", "frac": 0.1,
             "seeds": list(range(10)), "provenance": {"commit": "run1", "tree_clean": True, "config_sha256": sha,
                                                      "files_sha256": {"T1/E9.5_RNA.h5ad": "f"}}}))
    return d


@pytest.fixture
def t11_dir(vm, tmp_path, monkeypatch):
    d = tmp_path / "t11_tasks"
    d.mkdir()
    rows = [{"seed": s, "scores": {"floor": 50.0, "b0": 49.0 + 0.1 * s, "g_mom": 50.0 + 0.1 * s}} for s in range(10)]
    for name in vm.CUTS.values():
        (d / name).write_text(json.dumps({"rows": rows}))
    monkeypatch.setattr(vm, "T11_TASKS", d)
    return d


def _matching_rows(gain_flip, gain_perm, gain_null=None, held_gain=None, **kw):
    """Rows whose b0 and g_mom match try T11's committed scores; each control's gain per seed is g_mom minus it."""
    b0 = [49.0 + 0.1 * s for s in range(10)]
    g_mom = [50.0 + 0.1 * s for s in range(10)]
    flip = [g - gain_flip[i] for i, g in enumerate(g_mom)]
    perm = [g - gain_perm[i] for i, g in enumerate(g_mom)]
    null = None if gain_null is None else [g - gain_null[i] for i, g in enumerate(g_mom)]
    held = None if held_gain is None else [b + held_gain[i] for i, b in enumerate(b0)]
    return _run_rows(b0, g_mom, flip, perm, null=null, held=held, **kw)


def test_placebo_survives_only_when_each_control_is_beaten_and_the_held_arm_keeps_its_gain(vm, t11_dir, tmp_path):
    sha = vm.harness._sha256(vm.harness.CONFIG)
    good = {k: _matching_rows([1.0] * 10, [0.8] * 10) for k in vm.CUTS}
    run = _write_run(tmp_path, "good", good, sha)
    vm.main(["placebo", "--runs", str(run), "--out", "p_good"])
    rec = json.loads((vm.RESULTS / "p_good" / "placebo.json").read_text())
    assert rec["survives"] is True and rec["reproduction_of_t11"]["reproduced"] is True
    assert set(rec["placebos"]) == {"g_mom_flip", "g_mom_perm", "g_mom_null"}
    assert rec["held"]["try"] == "T24" and rec["held"]["passes"] is True and rec["held"]["void"] is False
    assert rec["held"]["decision"]["base"] == "b0" and rec["held"]["decision"]["change"] == "g_mom_held"
    assert rec["held"]["decision"]["tasks"]["T1/E9.5->E10.5"]["mean_gain"] == pytest.approx(0.8)
    dec = rec["placebos"]["g_mom_perm"]["decision"]
    assert dec["keep"] and dec["base"] == "g_mom_perm" and dec["change"] == "g_mom"
    assert dec["tasks"]["T1/E9.5->E10.5"]["mean_gain"] == pytest.approx(0.8)
    assert sorted(dec["counted"]) == sorted(vm.CUTS) and "T1/E8.5->E9.5" in dec["excluded"]
    assert rec["placebos"]["g_mom_flip"]["group_margins"]["T1/E9.5->E10.5"]["n_seeds"] == 10
    assert rec["stamp"]["commit"] == "abc123" and rec["stamp"]["tree_clean"] is True
    assert len(rec["stamp"]["inputs"]) == 4
    # the metric panels are in the record per arm and cut, and the spread reading reads them
    panel = rec["panels"]["T1/E9.5->E10.5"]
    assert set(panel) == set(ARMS) and panel["g_mom"]["variance_ratio"] == pytest.approx(0.712)
    assert panel["b0"]["mmd_u"] == pytest.approx(0.12) and "energy_distance" in panel["g_mom_null"]
    assert "pb_rel_err" in panel["g_mom_perm"]
    spread = rec["spread"]
    assert spread["t21_is_a_spread_control"] is True and spread["tolerance"] == 0.003
    cut = spread["per_cut"]["T1/E9.5->E10.5"]
    # the flip's variance ratio against b0 minus what g_mom added (the sign-odd reading of the row)
    assert cut["flip_if_sign_odd"] == pytest.approx(0.70 - 0.012) and cut["flip_near_sign_odd_expectation"]
    assert cut["perm_within_tolerance"] and spread["t20_sign_odd_reading_holds"] is True
    text = (vm.RESULTS / "p_good" / "placebo.txt").read_text()
    assert "SURVIVES" in text and "T20 g_mom_flip: beaten" in text and "T21 g_mom_perm: beaten" in text
    assert "T23 g_mom_null: beaten" in text and "T24 g_mom_held" in text and "keeps its gain" in text
    assert "variance_ratio" in text and "T21 may be read as a spread control" in text
    assert "neither is a spread control or evidence of momentum" in text
    assert "one fit and one permutation draw" in text
    with pytest.raises(SystemExit, match="already holds"):
        vm.main(["placebo", "--runs", str(run), "--out", "p_good"])


@pytest.mark.parametrize("flip_gain, perm_gain, null_gain, held_gain, survives", [
    ([1.0] * 10, [0.1] * 10, None, None, False),  # the permutation is not beaten by more than 0.3
    ([0.25] * 10, [1.0] * 10, None, None, False),  # the flip is within 0.3 of g_mom
    ([1.0] * 7 + [-0.5] * 3, [1.0] * 10, None, None, False),  # mean above 0.3 but only 7 of 10 seeds ahead
    ([1.0] * 10, [0.9] * 10, [0.1] * 10, None, False),  # the no-development arm is not beaten (T23)
    ([1.0] * 10, [0.9] * 10, None, [0.1] * 10, False),  # the held arm keeps no gain (T24): a per-type mean effect
    ([1.0] * 10, [0.9] * 10, None, None, True),
])
def test_placebo_rule_needs_mean_and_seeds_on_every_cut(vm, t11_dir, tmp_path, flip_gain, perm_gain, null_gain,
                                                        held_gain, survives):
    sha = vm.harness._sha256(vm.harness.CONFIG)
    rows = {k: _matching_rows(flip_gain, perm_gain, null_gain, held_gain) for k in vm.CUTS}
    run = _write_run(tmp_path, "r", rows, sha)
    vm.main(["placebo", "--runs", str(run), "--out", "p"])
    rec = json.loads((vm.RESULTS / "p" / "placebo.json").read_text())
    assert rec["survives"] is survives
    # one cut failing is enough: g_mom beaten only on the whole cut
    if survives:
        worse = {**rows, "T1/E9.5->E10.5:chamber": _matching_rows([0.1] * 10, [0.1] * 10, null_gain, held_gain)}
        run2 = _write_run(tmp_path, "r2", worse, sha)
        vm.main(["placebo", "--runs", str(run2), "--out", "p2"])
        assert json.loads((vm.RESULTS / "p2" / "placebo.json").read_text())["survives"] is False


def test_the_held_arms_reading_is_void_when_its_de_skills_drift(vm, t11_dir, tmp_path):
    """T24's validity, as T9's: the held arm's de_score and de_direction skills within 0.02 of b0's on each cut (mean
    over seeds), otherwise the pseudobulk was not held and the reading is void, whatever the score."""
    sha = vm.harness._sha256(vm.harness.CONFIG)
    rows = {k: _matching_rows([1.0] * 10, [1.0] * 10, de_skill={"g_mom_held": (0.53, 0.5)}) for k in vm.CUTS}
    run = _write_run(tmp_path, "void_held", rows, sha)
    vm.main(["placebo", "--runs", str(run), "--out", "p_void_held"])
    rec = json.loads((vm.RESULTS / "p_void_held" / "placebo.json").read_text())
    assert rec["held"]["void"] is True and rec["held"]["passes"] is False and rec["survives"] is False
    assert rec["held"]["validity"]["de_skill_gap_vs_b0"]["T1/E9.5->E10.5"]["de_score"] == pytest.approx(0.03)
    assert "VOID" in (vm.RESULTS / "p_void_held" / "placebo.txt").read_text()


def test_t21_may_be_read_as_a_spread_control_only_within_its_variance_gate(vm, t11_dir, tmp_path):
    """The permuted arm's variance_ratio within 0.003 of g_mom's on both cuts, or else T21 says only that the assignment
    matters; a pass of the rule is unaffected."""
    sha = vm.harness._sha256(vm.harness.CONFIG)
    near = {k: _matching_rows([1.0] * 10, [1.0] * 10, vr={"g_mom_perm": 0.7115}) for k in vm.CUTS}
    run = _write_run(tmp_path, "near", near, sha)
    vm.main(["placebo", "--runs", str(run), "--out", "p_near"])
    assert json.loads((vm.RESULTS / "p_near" / "placebo.json").read_text())["spread"]["t21_is_a_spread_control"] is True
    far = {k: _matching_rows([1.0] * 10, [1.0] * 10, vr={"g_mom_perm": 0.702}) for k in vm.CUTS}
    run2 = _write_run(tmp_path, "far", far, sha)
    vm.main(["placebo", "--runs", str(run2), "--out", "p_far"])
    rec = json.loads((vm.RESULTS / "p_far" / "placebo.json").read_text())
    assert rec["spread"]["t21_is_a_spread_control"] is False and rec["survives"] is True
    assert rec["spread"]["per_cut"]["T1/E9.5->E10.5"]["perm_within_tolerance"] is False
    assert "T21 may NOT be read as a spread control" in (vm.RESULTS / "p_far" / "placebo.txt").read_text()
    # and the flip's reading is only reported: a flip far from b0 minus what g_mom added is said so
    odd = {k: _matching_rows([1.0] * 10, [1.0] * 10, vr={"g_mom_flip": 0.70}) for k in vm.CUTS}
    run3 = _write_run(tmp_path, "odd", odd, sha)
    vm.main(["placebo", "--runs", str(run3), "--out", "p_odd"])
    spread_odd = json.loads((vm.RESULTS / "p_odd" / "placebo.json").read_text())["spread"]
    assert spread_odd["t20_sign_odd_reading_holds"] is False


def test_placebo_is_void_when_t11_does_not_reproduce(vm, t11_dir, tmp_path):
    sha = vm.harness._sha256(vm.harness.CONFIG)
    off = _matching_rows([1.0] * 10, [1.0] * 10)
    off[3]["scores"]["g_mom"] += 0.5  # one seed's g_mom differs from the committed try T11 by more than 0.02
    rows = {k: (off if k == "T1/E9.5->E10.5" else _matching_rows([1.0] * 10, [1.0] * 10)) for k in vm.CUTS}
    run = _write_run(tmp_path, "void", rows, sha)
    with pytest.raises(SystemExit, match="does not reproduce"):
        vm.main(["placebo", "--runs", str(run), "--out", "p_void"])
    rec = json.loads((vm.RESULTS / "p_void" / "placebo.json").read_text())
    assert rec["survives"] is False and rec["reproduction_of_t11"]["reproduced"] is False
    assert rec["reproduction_of_t11"]["largest_difference"]["g_mom"] == pytest.approx(0.5)
    assert "NO, the reading is void" in (vm.RESULTS / "p_void" / "placebo.txt").read_text()


def test_placebo_reproduction_allows_only_a_small_float_difference_in_g_mom(vm, t11_dir, tmp_path):
    """b0 must reproduce exactly (within 0.02); g_mom, whose fit carries PCA and k-means floats, within 0.05."""
    sha = vm.harness._sha256(vm.harness.CONFIG)
    ok = _matching_rows([1.0] * 10, [1.0] * 10)
    ok[0]["scores"]["g_mom"] += 0.04  # g_mom's own float wobble: within 0.05
    ok[0]["scores"]["g_mom_flip"] += 0.04  # the placebo moves with the run's own fit
    ok[0]["scores"]["g_mom_perm"] += 0.04
    run = _write_run(tmp_path, "ok", {k: ok for k in vm.CUTS}, sha)
    vm.main(["placebo", "--runs", str(run), "--out", "p_ok"])
    assert json.loads((vm.RESULTS / "p_ok" / "placebo.json").read_text())["reproduction_of_t11"]["reproduced"] is True
    off = _matching_rows([1.0] * 10, [1.0] * 10)
    off[0]["scores"]["b0"] += 0.03  # b0 has no float-sensitive fit: 0.02 is its bar
    run2 = _write_run(tmp_path, "off", {k: off for k in vm.CUTS}, sha)
    with pytest.raises(SystemExit, match="does not reproduce"):
        vm.main(["placebo", "--runs", str(run2), "--out", "p_off"])


def test_placebo_refuses_a_missing_cut_or_a_run_without_the_placebos(vm, t11_dir, tmp_path):
    sha = vm.harness._sha256(vm.harness.CONFIG)
    run = _write_run(tmp_path, "r", {k: _matching_rows([1.0] * 10, [1.0] * 10) for k in vm.CUTS}, sha)
    (run / "T1_E9.5_E10.5.json").unlink()
    with pytest.raises(SystemExit, match="no result for"):
        vm.main(["placebo", "--runs", str(run), "--out", "p3"])
    run2 = _write_run(tmp_path, "r2", {k: _matching_rows([1.0] * 10, [1.0] * 10) for k in vm.CUTS}, sha)
    for f in run2.iterdir():
        r = json.loads(f.read_text())
        r["candidates"] = ["floor", "b0", "g_mom", "g_mom_flip", "g_mom_perm"]  # no no-development or held arm
        f.write_text(json.dumps(r))
    with pytest.raises(SystemExit, match="the controls and the held arm"):
        vm.main(["placebo", "--runs", str(run2), "--out", "p4"])


# ------------------------------------------------------------------------------------------------ the composition arms
MIXFREE_ARMS = ("floor", "b0", "g_mom", "g_mom_null")


def _mixfree_result(key, b0, g_mom, null=None):
    """One within-type read's harness result (per-seed scores, group skills and the metric panel)."""
    rows = []
    for i, s in enumerate(range(10)):
        scores = {"floor": 50.0, "b0": b0[i], "g_mom": g_mom[i],
                  "g_mom_null": g_mom[i] - 1.0 if null is None else null[i]}
        groups = {arm: {"de_recovery": 0.6, "direction": 0.6, "distribution": 0.7 + (v - 50) / 100,
                        "covariance": 0.7 + (v - 50) / 200} for arm, v in scores.items()}
        metrics = {arm: {**PANEL, "variance_ratio": 0.7} for arm in MIXFREE_ARMS}
        rows.append({"seed": s, "scores": scores, "groups": groups, "metrics": metrics})
    return {"task": {"key": key}, "candidates": list(MIXFREE_ARMS), "rows": rows}


def _mixfree_results(vm, cm_gain, endo_gain):
    b0 = [49.0 + 0.1 * s for s in range(10)]
    return {vm.MIXFREE_CM: _mixfree_result(vm.MIXFREE_CM, b0, [b + g for b, g in zip(b0, cm_gain, strict=True)]),
            vm.MIXFREE_ENDO: _mixfree_result(vm.MIXFREE_ENDO, b0, [b + g for b, g in zip(b0, endo_gain, strict=True)])}


@pytest.mark.parametrize("cm_gain, endo_gain, passes", [
    ([0.5] * 10, [0.0] * 10, True),  # a gain on cardiomyocytes, endocardium flat: passes
    ([0.5] * 10, [0.5] * 10, True),
    ([0.5] * 10, [-0.2] * 10, True),  # an endocardium loss within the veto
    ([0.5] * 10, [-0.4] * 10, False),  # the endocardium veto
    ([0.2] * 10, [0.5] * 10, False),  # the mean on cardiomyocytes is below the bar
    ([1.0] * 7 + [-0.1] * 3, [0.5] * 10, False),  # ahead in 7 of 10 seeds only
])
def test_t25_needs_the_cardiomyocyte_gain_and_no_endocardium_veto(vm, cm_gain, endo_gain, passes):
    read = vm.t25_reading(_mixfree_results(vm, cm_gain, endo_gain))
    assert read["passes"] is passes
    assert read["vs_b0"][vm.MIXFREE_CM]["mean_gain"] == pytest.approx(np.mean(cm_gain))
    assert read["vs_null_reported"][vm.MIXFREE_CM]["mean_gain"] == pytest.approx(1.0)
    assert set(read["panels"][vm.MIXFREE_ENDO]) == set(MIXFREE_ARMS)


def test_t25_treats_a_read_with_missing_seeds_as_not_passing(vm):
    results = _mixfree_results(vm, [0.5] * 10, [0.0] * 10)
    results[vm.MIXFREE_ENDO]["rows"] = results[vm.MIXFREE_ENDO]["rows"][:5]
    read = vm.t25_reading(results)
    assert read["passes"] is False and read["vs_b0"][vm.MIXFREE_ENDO]["complete"] is False
    assert read["vs_b0"][vm.MIXFREE_ENDO]["group_margins"] is None


def _mixfree_io(vm, monkeypatch, tmp_path, results, b0_shift=0.0, fail=None):
    """The step's inputs without any file of the project: the two reads, run_task's results and the committed b0."""
    gate = tmp_path / "gate"
    gate.mkdir()
    for key, fname in vm.MIXFREE_READS.items():
        rows = [{"seed": r["seed"], "scores": {"b0": r["scores"]["b0"] + b0_shift}} for r in results[key]["rows"]]
        (gate / fname).write_text(json.dumps({"rows": rows}))
    monkeypatch.setattr(vm, "MIXFREE_B0", gate)
    original = vm.reproduce_late_proxy_b0
    monkeypatch.setattr(vm, "reproduce_late_proxy_b0", lambda res, folder=gate: original(res, folder))
    config = tmp_path / "reads.yaml"
    config.write_text("reads\n")
    monkeypatch.setattr(vm, "MIXFREE_CONFIG", config)

    class Read:
        def __init__(self, key):
            self.key = key

    monkeypatch.setattr(vm.harness, "load_tasks", lambda path=None: {vm.MIXFREE_BOARD: [Read(k) for k in results]})
    seen = {}

    def run(task, candidates, seeds, data_dir, build):
        if task.key == fail:
            raise ValueError("cannot build this read")
        seen[task.key] = ([c.name for c in candidates], tuple(seeds), build)
        return results[task.key]

    monkeypatch.setattr(vm.harness, "run_task", run)
    return seen


def test_mixfree_runs_the_arms_on_both_reads_and_records_the_rule(vm, monkeypatch, tmp_path):
    results = _mixfree_results(vm, [0.5] * 10, [0.1] * 10)
    seen = _mixfree_io(vm, monkeypatch, tmp_path, results)
    vm.main(["mixfree", "--out", "t25"])
    assert set(seen) == set(vm.MIXFREE_READS)
    assert all(v == (list(MIXFREE_ARMS), tuple(range(10)), vm.harness.BUILD) for v in seen.values())
    rec = json.loads((vm.RESULTS / "t25" / "mixfree.json").read_text())
    assert rec["passes"] is True and rec["reproduction_of_late_proxy_b0"]["reproduced"] is True
    assert rec["rule"]["primary"] == vm.MIXFREE_CM and rec["stamp"]["tree_clean"] is True
    assert len(rec["stamp"]["inputs"]) == 5
    assert sorted(p.name for p in (vm.RESULTS / "t25" / "tasks").iterdir()) == sorted(vm.MIXFREE_READS.values())
    text = (vm.RESULTS / "t25" / "mixfree.txt").read_text()
    assert "PASSES" in text and "g_mom minus g_mom_null" in text and "one comparison, not a replication" in text
    with pytest.raises(SystemExit, match="already holds"):
        vm.main(["mixfree", "--out", "t25"])


def test_mixfree_is_void_when_the_committed_b0_does_not_reproduce(vm, monkeypatch, tmp_path):
    _mixfree_io(vm, monkeypatch, tmp_path, _mixfree_results(vm, [0.5] * 10, [0.1] * 10), b0_shift=0.05)
    with pytest.raises(SystemExit, match="does not reproduce"):
        vm.main(["mixfree", "--out", "t25_void"])
    rec = json.loads((vm.RESULTS / "t25_void" / "mixfree.json").read_text())
    assert rec["passes"] is False and rec["reproduction_of_late_proxy_b0"]["reproduced"] is False
    assert "VOID" in (vm.RESULTS / "t25_void" / "mixfree.txt").read_text()


def test_mixfree_stops_void_and_uncounted_when_a_read_cannot_be_built(vm, monkeypatch, tmp_path):
    results = _mixfree_results(vm, [0.5] * 10, [0.1] * 10)
    _mixfree_io(vm, monkeypatch, tmp_path, results, fail=vm.MIXFREE_ENDO)
    with pytest.raises(SystemExit, match="void and uncounted"):
        vm.main(["mixfree", "--out", "t25_failed"])
    assert not (vm.RESULTS / "t25_failed" / "mixfree.json").exists()
    assert (vm.RESULTS / "t25_failed" / "tasks" / f"{vm.MIXFREE_READS[vm.MIXFREE_ENDO]}.failed.txt").exists()


def test_the_mixfree_reads_exist_in_the_committed_config_with_a_start_stage(vm):
    tasks = {t.key: t for t in harness.load_tasks(vm.MIXFREE_CONFIG)[vm.MIXFREE_BOARD]}
    for key in vm.MIXFREE_READS:
        assert CANDIDATES["g_mom"].applies(tasks[key]) and CANDIDATES["b0"].applies(tasks[key])
        assert CANDIDATES["g_mom_null"].applies(tasks[key])


def test_oracle_shares_are_the_cuts_chamber_categories_renormalised(vm):
    cats = ["cm"] * 2 + ["endocardium"] * 3 + ["blood"] * 3 + ["surrounding"] * 2
    assert vm.oracle_shares(cats) == pytest.approx({"cm": 0.25, "endocardium": 0.375, "blood": 0.375})


def _gain_rows(gains, key="chamber"):
    base = {"de_recovery": 0.6, "direction": 0.6, "distribution": 0.7, "covariance": 0.7}
    # b0 at 0.0, so that a gain of exactly +0.3 is exactly 0.3 in floating point (the bar is inclusive)
    return [{"seed": s, key: {"b0": {"score": 0.0, "groups": dict(base)},
                              "g_mom": {"score": g, "groups": {**base, "distribution": 0.7 + g / 30}}}}
            for s, g in zip(range(len(gains)), gains, strict=True)]


@pytest.mark.parametrize("gains, passes, at_bar", [
    ([0.3] * 8 + [0.0] * 2, True, 8),  # the bar is inclusive: exactly +0.3 counts, eight seeds suffice
    ([0.29] * 8 + [0.9] * 2, False, 2),  # just under the bar in eight seeds
    ([1.0] * 7 + [0.2] * 3, False, 7),  # only seven seeds at or above the bar
])
def test_t22_reading_counts_seeds_at_or_above_the_bar(vm, gains, passes, at_bar):
    r = vm.t22_reading(_gain_rows(gains))
    assert r["passes"] is passes and r["seeds_at_or_above_bar"] == at_bar and r["n_seeds"] == 10
    assert r["group_margins"]["distribution"] == pytest.approx(30 * (np.mean(gains) / 30))
    assert r["mean_gain"] == pytest.approx(np.mean(gains))


def test_t22_reading_does_not_pass_on_fewer_seeds(vm):
    assert vm.t22_reading(_gain_rows([1.0] * 5))["passes"] is False


T22_BOARD = Board(key="pseudo:t22", task="T1", genes=list(GENES), needs_coords=False, min_cells=50, max_cells=400,
                  truth_cells=60, ref_cells=40)


def _t22_io(vm, monkeypatch, tmp_path, full=False, committed_mean=None, seeds_patch=None):
    """A run of `composition` on synthetic stages: the tasks are the real ones (a board-free config read); every data,
    board and committed-grid read is replaced. Returns the fixture objects."""
    ref, start, tgt = _small_inputs(seed=31)
    t1_scope.classify(tgt, source="E10.5_RNA")
    cat = tgt.obs["category"].astype(str).to_numpy()
    keep = np.flatnonzero(np.isin(cat, ["cm", "endocardium", "blood"]))
    extra = np.flatnonzero(cat == "surrounding")[:len(keep) // 19]
    rows = np.sort(np.concatenate([keep, extra]))
    cut, cut_cat = tgt[rows].copy(), cat[rows]
    shares = vm.oracle_shares(cut_cat)
    grid = {"pred_mix": {vm.ORACLE_ROW: shares}, "mean": {vm.ORACLE_ROW: {vm.ORACLE_ROW: committed_mean or 0.0}}}
    grid_path = tmp_path / "payoff.json"
    grid_path.write_text(json.dumps(grid))
    monkeypatch.setattr(vm, "T1_GRID", grid_path)
    monkeypatch.setattr(vm.payoff, "chamber_cut", lambda t, data, genes: (cut, cut_cat))
    monkeypatch.setattr(vm, "read_stage", lambda path, genes: (tgt if "E10.5" in str(path) else ref)[:, genes].copy())

    def load(path, lineage=False):
        return (tgt if "E10.5" in str(path) else ref if "E9.5" in str(path) else start).copy()

    monkeypatch.setattr(vm, "load_stage", load)
    monkeypatch.setattr(vm, "load_boards", lambda: {"T1:val": T22_BOARD})
    monkeypatch.setattr(vm.harness, "task_board", lambda t, boards, n_t, n_r: T22_BOARD)
    monkeypatch.setattr(vm.t1_scope, "held_flat_genes", lambda *a, **k: list(HELD))
    monkeypatch.setattr(vm, "sha256_file", lambda path: "0" * 64)
    if seeds_patch is not None:
        monkeypatch.setattr(vm, "SEEDS", seeds_patch)
    harness._G_MOM_PREP.clear()
    return ref, start, tgt, shares


def test_composition_builds_both_arms_at_the_oracle_mix_and_reads_the_rule(vm, monkeypatch, tmp_path):
    ref, start, tgt, shares = _t22_io(vm, monkeypatch, tmp_path)
    vm.main(["composition", "--seeds", "0", "--out", "t22", "--data", str(tmp_path)])
    rec = json.loads((vm.RESULTS / "t22" / "composition.json").read_text())
    assert rec["oracle_shares"] == pytest.approx(shares) and rec["complete_seeds"] is False
    row = rec["rows"][0]
    assert set(row) >= {"seed", "chamber", "whole", "companion"} and row["chamber"]["floor"] == pytest.approx(50.0)
    for key in ("chamber", "whole"):
        assert row[key]["b0"]["score"] != row[key]["g_mom"]["score"]
        assert set(row[key]["g_mom"]["groups"]) == {"de_recovery", "direction", "distribution", "covariance"}
        assert set(row[key]) >= {"b0", "g_mom", "g_mom_held", "floor"}
        # the metric panel of every arm is kept (the spread reading reads variance_ratio)
        panel = set(row[key]["g_mom"]["metrics"])
        assert {"variance_ratio", "mmd_u", "variogram", "energy_distance", "pb_rel_err"} <= panel
    assert rec["chamber"]["passes"] is False and rec["reproduction_of_oracle_b0"]["reproduced"] is None
    assert "variance_ratio" in rec["chamber"]["panel"]["g_mom_held"]
    # the companion: the arms at the prediction's own mix against the chamber truth re-sampled to that mix
    assert set(row["companion"]) >= {"floor", "b0", "g_mom"} and row["companion"]["floor"] == pytest.approx(50.0)
    assert rec["companion"]["n_seeds"] == 1 and "group_margins" in rec["companion"]
    assert rec["companion_truth"]["resample_seed"] == 20261110
    want = harness.t1_chamber_shares(harness._shallow(ref), "E9.5_RNA")
    mix = rec["companion_truth"]["mix"]
    chamber_part = {k: v for k, v in mix.items() if k in want}
    assert sum(mix.values()) == pytest.approx(1.0) and set(want) <= set(mix)
    assert {k: v / sum(chamber_part.values()) for k, v in chamber_part.items()} == pytest.approx(want)
    assert rec["chamber_held"]["n_seeds"] == 1
    # the median library size report, DATA_POLICY 2026-11-04 (f)(3): shares and the two medians, no cell number
    lib = rec["median_library_size"]
    assert lib["ratio"] == pytest.approx(lib["oracle_mix_population"] / lib["e10_5_released"])
    assert lib["beyond_5_percent"] == (abs(lib["ratio"] - 1.0) > 0.05)
    text = (vm.RESULTS / "t22" / "composition.txt").read_text()
    assert "not all ten seeds" in text and "Oracle mix" in text and "Whole cut" in text
    assert "Companion (reported, deciding nothing" in text and "g_mom_held minus b0" in text
    assert "Median library size" in text and "one transition and one fit" in text
    assert "a whole-mix oracle base could be drawn" in text and "chamber-scoped candidate" in text
    assert rec["stamp"]["commit"] == "abc123" and any("E10.5" in k for k in rec["stamp"]["inputs"])
    with pytest.raises(SystemExit, match="already holds"):
        vm.main(["composition", "--seeds", "0", "--out", "t22", "--data", str(tmp_path)])


def test_composition_b0_arm_is_the_b0_population_and_the_g_mom_arm_edits_it(vm, monkeypatch, tmp_path):
    ref, start, tgt, shares = _t22_io(vm, monkeypatch, tmp_path)
    ctx = harness.BuildContext(task=harness.TaskView.of(next(t for t in harness.load_tasks()["T1:test"]
                                                             if t.key == vm.TASK_CHAMBER)),
                               board=harness._candidate_board(T22_BOARD), inputs={"ref": ref, "start": start},
                               seed=1000003, derived={"t1_held_flat": tuple(HELD)})
    arms, population = vm.oracle_arms(ctx, shares, HELD)
    b0, gm, held_arm = arms["b0"], arms["g_mom"], arms["g_mom_held"]
    # the b0 arm is B0's own population at the mix, as payoff.py builds it through build_submission
    n = harness.count_for(harness._candidate_board(T22_BOARD), ref.n_obs)
    want, _ = t1_scope.build_submission(harness._shallow(ref), T22_BOARD, shares, n, 1000003,
                                        start=harness._shallow(start), damping=harness.B0_DAMPING, held_flat=list(HELD),
                                        ref_source="E9.5_RNA", start_source="E8.5_RNA")
    assert np.array_equal(b0.X, np.asarray(want.X)) and b0.X.shape[0] == n
    assert gm.X.shape == b0.X.shape and not np.array_equal(gm.X, b0.X)
    assert held_arm.X.shape == b0.X.shape and not np.array_equal(held_arm.X, gm.X)
    # the g_mom arm's edit is the population at the oracle mix through g_mom's own fold factors (R1's held labels stay)
    r1, s1 = harness._shallow(ref), harness._shallow(start)
    pop = ts.b0_population(r1, s1, shares, n, 1000003, list(HELD), ref_source="E9.5_RNA", start_source="E8.5_RNA",
                           damping=harness.B0_DAMPING)
    blood = pop.obs["celltype"].astype(str).to_numpy() == "Blood"
    assert blood.any() and np.array_equal(gm.X[blood], b0.X[blood]) and np.array_equal(held_arm.X[blood], b0.X[blood])
    assert population.n_obs == n and "counts" in population.layers
    # the held arm is the same factors through the rescale: it holds each category's pseudobulk per gene in expectation
    # the fold's factors are fitted on the reference's own chamber categories, as g_mom's
    assert harness._G_MOM_PREP["fold"]["factors"].metacells_on == "later"


def test_companion_mix_is_the_predictions_chamber_mix_with_the_cuts_surrounding_share(vm):
    own = {"cm": 0.6, "endocardium": 0.2, "blood": 0.2}
    cut_mix = {"cm": 0.2, "endocardium": 0.3, "blood": 0.25, "epicardium": 0.2, "surrounding": 0.05}
    mix, dropped = vm.companion_mix(own, cut_mix)
    assert dropped == [] and mix["surrounding"] == pytest.approx(0.05) and sum(mix.values()) == pytest.approx(1.0)
    assert mix["cm"] == pytest.approx(0.6 * 0.95) and "epicardium" not in mix
    # a category of the prediction that the cut lacks is dropped and renormalised, and named
    mix2, dropped2 = vm.companion_mix({"cm": 0.5, "vascular_endothelium": 0.5}, {"cm": 0.95, "surrounding": 0.05})
    assert dropped2 == ["vascular_endothelium"] and mix2["cm"] == pytest.approx(0.95)


def test_composition_is_void_when_the_oracle_base_does_not_reproduce(vm, monkeypatch, tmp_path):
    ref, start, tgt, shares = _t22_io(vm, monkeypatch, tmp_path, committed_mean=12.0, seeds_patch=(0,))
    with pytest.raises(SystemExit, match="does not reproduce the committed t1-grid"):
        vm.main(["composition", "--seeds", "0", "--out", "t22v", "--data", str(tmp_path)])
    rec = json.loads((vm.RESULTS / "t22v" / "composition.json").read_text())
    assert rec["reproduction_of_oracle_b0"]["reproduced"] is False and rec["complete_seeds"] is True
    text = (vm.RESULTS / "t22v" / "composition.txt").read_text()
    assert "VOID" in text and "NO, the reading is void" in text and "PASSES" not in text.split("\n")[1]


def test_composition_refuses_another_oracle_mix(vm, monkeypatch, tmp_path):
    ref, start, tgt, shares = _t22_io(vm, monkeypatch, tmp_path)
    bad = {"pred_mix": {vm.ORACLE_ROW: {**shares, "cm": shares["cm"] + 0.1}},
           "mean": {vm.ORACLE_ROW: {vm.ORACLE_ROW: 0.0}}}
    (tmp_path / "payoff.json").write_text(json.dumps(bad))
    with pytest.raises(SystemExit, match="oracle mix is not the committed"):
        vm.main(["composition", "--seeds", "0", "--out", "t22b", "--data", str(tmp_path)])


# ------------------------------------------------------------------------------------------------ diagnostic N21
def test_leading_axis_and_cosine_recover_a_planted_direction(vm):
    rng = np.random.default_rng(3)
    loading = rng.normal(size=30)
    loading /= np.linalg.norm(loading)
    scores = np.array([2.0, -1.0, 1.5, -2.5, 0.5, -0.5])
    logf = scores[:, None] * loading[None, :]
    genes = [f"g{i}" for i in range(30)]
    cf = ts.CategoryFactors("cm", {}, np.ones(6, dtype=np.int64), logf)
    other = ts.CategoryFactors("cm", {}, np.array([3, 1, 2, 2, 1, 3]),
                               2.0 * logf + rng.normal(0, 1e-9, size=logf.shape))
    u = vm.leading_axis(cf, genes, genes)
    assert abs(u @ loading) == pytest.approx(1.0, abs=1e-6)
    assert vm.axis_cosine(u, -u) == pytest.approx(1.0)  # an axis has no sign
    v = vm.leading_axis(other, genes, genes)
    assert vm.axis_cosine(u, v) == pytest.approx(1.0, abs=1e-6)
    orth = rng.normal(size=30)
    orth -= (orth @ loading) * loading
    assert vm.axis_cosine(loading, orth) == pytest.approx(0.0, abs=1e-9)
    zero = ts.CategoryFactors("cm", {}, np.ones(3, dtype=np.int64), np.zeros((3, 30)))
    assert np.isnan(vm.axis_cosine(vm.leading_axis(zero, genes, genes), u))


def _two_fits(vm, same=True, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(40)]
    axis_a, axis_b = rng.normal(size=40), rng.normal(size=40)
    scores = rng.normal(size=8)
    axis2 = axis_a if same else axis_b
    unit = np.ones(8, dtype=np.int64)
    a = ts.CategoryFactors("cm", {}, unit, scores[:, None] * axis_a[None, :] / np.linalg.norm(axis_a))
    b = ts.CategoryFactors("cm", {}, unit,
                           scores[:, None] * axis2[None, :] / np.linalg.norm(axis2) + rng.normal(0, 0.01, size=(8, 40)))
    def fit(cf):
        return ts.StateFactors(genes=genes, metacells_on="later", categories={"cm": cf}, skipped={})

    return fit(a), fit(b)


def test_concordance_is_high_on_the_same_direction_and_low_on_another(vm):
    a, b = _two_fits(vm, same=True)
    near = vm.concordance(a, b, "cm", set())
    a2, b2 = _two_fits(vm, same=False)
    far = vm.concordance(a2, b2, "cm", set())
    assert near["axis_abs_cosine"] > 0.99 and far["axis_abs_cosine"] < 0.6
    assert near["size_spearman"] > 0.9 and far["size_spearman"] < 0.9
    # the held-flat genes leave the comparison, and the genes are read by name
    held = {"g0", "g3"}
    swapped = ts.StateFactors(genes=list(reversed(b.genes)), metacells_on="later", skipped={},
                              categories={"cm": ts.CategoryFactors("cm", {}, np.ones(8, dtype=np.int64),
                                                                   b.categories["cm"].logf[:, ::-1])})
    assert vm.concordance(a, swapped, "cm", held)["axis_abs_cosine"] == pytest.approx(
        vm.concordance(a, b, "cm", held)["axis_abs_cosine"], abs=1e-12)


def test_gene_size_and_profile_spearman(vm):
    cf = ts.CategoryFactors("cm", {}, np.array([1, 1, 2]), np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 0.0]]))
    size = vm.gene_size(cf, ["a", "b"], ["a", "b"])
    assert size[0] == pytest.approx(np.sqrt(0.5)) and size[1] == 0.0
    assert np.isnan(vm.profile_spearman(np.zeros(3), np.array([1.0, 2.0, 3.0])))
    assert vm.profile_spearman(np.array([1.0, 2.0, 3.0]), np.array([1.0, 3.0, 2.0])) == pytest.approx(0.5)


def test_expected_depth_ratio_checks_the_constant_apply_factors_returns(vm):
    """The check reads the keep-depth constant the real `apply_factors` returned for these factors: it is 1 for that
    constant and not 1 for another, so a change to the function's constant would be caught."""
    from vec.core import apply_factors

    rng = np.random.default_rng(1)
    C = sparse.csr_matrix(rng.poisson(1.0, size=(100, 30)).astype(np.int64))
    C = C + sparse.csr_matrix(np.eye(100, 30, dtype=np.int64))  # every cell has a count
    F = np.exp(rng.normal(0, 0.5, size=(4, 30)))
    row = rng.integers(0, 4, size=100)
    held = np.zeros(30, dtype=bool)
    held[:5] = True
    F[:, held] = 1.0
    _, scale, _ = apply_factors(C, F, row, held, keep_depth=True, seed=0)
    assert np.allclose(vm.expected_depth_ratio(C, F, row, held, scale), 1.0, atol=1e-12)
    assert not np.allclose(vm.expected_depth_ratio(C, F, row, held, scale * 1.01), 1.0, atol=1e-6)
    # apply_factors really keeps it in expectation: the realised totals agree on average
    depth = np.asarray(C.sum(1)).ravel().astype(np.float64)
    many = np.mean([np.asarray(apply_factors(C, F, row, held, keep_depth=True, seed=s)[0].sum(1)).ravel() / depth
                    for s in range(30)], axis=0)
    assert abs(float(many.mean()) - 1.0) < 0.02


def _amplified_factors(vm, seed=3, K=6, n_per=30, genes=30, slope=0.4, noise=0.05):
    """A stage whose metacells differ (so the own deviation is not zero) and factors that are a known multiple of that
    deviation plus noise, to read the amplification regression on."""
    rng = np.random.default_rng(seed)
    gene_names = [f"g{i}" for i in range(genes)]
    C = sparse.csr_matrix(rng.poisson(5.0 * np.exp(rng.normal(0, 1.0, size=(K, genes)))[np.repeat(np.arange(K), n_per)]
                                      + 1e-9, size=(K * n_per, genes)).astype(np.int64))
    names = [f"c{i}" for i in range(K * n_per)]
    stage = ad.AnnData(X=log1p_cp10k(C), var=pd.DataFrame(index=gene_names), obs=pd.DataFrame(index=names))
    stage.layers["counts"] = C
    index = {n: int(i // n_per) for i, n in enumerate(names)}
    sizes = np.full(K, n_per)
    cf = ts.CategoryFactors("cm", index, sizes, np.zeros((K, genes)))
    dev = vm.own_deviation(stage, cf, gene_names, gene_names)
    logf = slope * dev + rng.normal(0, noise, size=dev.shape)
    return stage, ts.StateFactors(genes=gene_names, metacells_on="later", skipped={},
                                  categories={"cm": replace(cf, logf=logf)}), dev


def test_own_deviation_and_amplification_read_a_planted_self_amplifier(vm):
    stage, sf, dev = _amplified_factors(vm)
    amp = vm.amplification(sf.categories["cm"].logf, dev, sf.categories["cm"].sizes)
    assert amp["slope"] == pytest.approx(0.4, abs=0.05) and amp["r2"] > 0.9
    # cell-weighted, the deviation of the metacells has mean zero per gene
    w = sf.categories["cm"].sizes / sf.categories["cm"].sizes.sum()
    assert np.allclose((w[:, None] * dev).sum(0), 0.0, atol=1e-12)
    # projecting the component out removes what the regression explained; the held-flat genes are left alone
    proj, record = vm.projected_fit(sf, stage, ["cm"], {"g0"})
    assert record["cm"]["slope"] == pytest.approx(0.4, abs=0.05)
    assert np.array_equal(proj.categories["cm"].logf[:, 0], sf.categories["cm"].logf[:, 0])
    keep = sf.genes[1:]
    dev_keep = vm.own_deviation(stage, sf.categories["cm"], sf.genes, keep)
    again = vm.amplification(proj.categories["cm"].logf[:, 1:], dev_keep, sf.categories["cm"].sizes)
    assert abs(again["slope"]) < 0.02 and again["r2"] < 0.02


def test_amplification_of_zero_factors_is_not_a_number(vm):
    dev = np.ones((4, 3))
    amp = vm.amplification(np.zeros((4, 3)), dev, np.ones(4, dtype=np.int64))
    assert amp["slope"] == 0.0 and np.isnan(amp["r2"])
    amp2 = vm.amplification(np.ones((4, 3)), np.zeros((4, 3)), np.ones(4, dtype=np.int64))
    assert np.isnan(amp2["slope"]) and np.isnan(amp2["r2"])


def test_union_share_reads_the_factor_size_on_the_union_genes(vm):
    logf = np.array([[2.0, 0.0, 1.0, 1.0], [0.0, 2.0, 1.0, -1.0]])
    cf = ts.CategoryFactors("cm", {}, np.array([1, 1]), logf)
    genes = ["a", "b", "c", "held"]
    out = vm.union_share(cf, genes, {"c", "held", "gone"}, {"held"})  # held-flat genes are out of the comparison
    assert out["union_genes_in_fit"] == 1
    assert out["l2_share"] == pytest.approx(1.0 / 5.0, abs=1e-12)
    assert out["l1_share"] == pytest.approx(1.0 / 3.0, abs=1e-12)
    zero = ts.CategoryFactors("cm", {}, np.array([1, 1]), np.zeros((2, 4)))
    assert np.isnan(vm.union_share(zero, genes, {"c"}, set())["l2_share"])


def _n21_stages(seed=11):
    """E8.5-, E9.5- and E10.5-like stages with within-type structure on the chamber labels (the config classifies
    them), small so that the full and the split-half fits stay quick; every half keeps 20 cells of each category."""
    e85 = states_stage(seed, {"V-CM": 320, "Endothelium": 140, "Blood": 60, "Foregut": 80})
    e95 = states_stage(seed + 1, {"V-CM": 420, "Endocardium": 170, "Blood": 60, "Foregut": 100}, later=True)
    e105 = states_stage(seed + 2, {"V-CM": 460, "Endocardium": 190, "Blood": 70, "Foregut": 100}, later=True,
                        target=True)
    return {"E8.5": e85, "E9.5": e95, "E10.5": e105}


def test_n21_splits_each_category_in_two_by_its_own_seed(vm):
    stages = _n21_stages()
    a = vm.classified(stages["E10.5"], "E10.5_RNA", t1_scope.load_spec())
    rng = np.random.default_rng([vm.SPLIT_SEED, 0])
    m = vm.halves(a, rng)
    cats = a.obs["category"].astype(str).to_numpy()
    for c in set(cats):
        assert m[cats == c].sum() == (cats == c).sum() // 2
    assert np.array_equal(vm.halves(a, np.random.default_rng([vm.SPLIT_SEED, 0])), m)
    assert not np.array_equal(vm.halves(a, np.random.default_rng([vm.SPLIT_SEED, 1])), m)


def test_n21_category_reports_shares_only_and_holds_r1(vm):
    pop, sf = _toy_population(seed=2)
    pop.obs["celltype"] = ["V-CM"] * 200 + ["Blood"] * 100
    labels = pop.obs["celltype"].astype(str).to_numpy()
    no_edit = labels == "Blood"
    held = set()
    row, edited = vm.n21_category(pop, labels, no_edit, sf, held, "cm")
    assert edited.sum() == 200 and row["share_of_category_cells_touched"] == pytest.approx(1.0)
    assert row["share_of_category_cells_held_by_r1"] == 0.0 and row["held_labels_in_category"] == []
    row_b, edited_b = vm.n21_category(pop, labels, no_edit, sf, held, "blood")
    assert edited_b.sum() == 0 and row_b["share_of_category_cells_held_by_r1"] == pytest.approx(1.0)
    assert row_b["held_labels_in_category"] == ["Blood"] and row_b["factor_size"]["n_metacells"] == 2


def test_depth_checks_keep_the_expected_depth_and_leave_held_cells_alone(vm):
    pop, sf = _toy_population(seed=3)
    pop.obs["celltype"] = ["V-CM"] * 200 + ["Blood"] * 100
    no_edit = pop.obs["celltype"].astype(str).to_numpy() == "Blood"
    edited = ~no_edit
    held = {"g0", "g7"}
    out = vm.depth_checks(pop, no_edit, sf, held, edited)
    assert out["max_abs_expected_ratio_minus_1"] < 1e-9 and out["realised_depth_identical_when_held"] is True
    assert out["share_of_fitted_cells_edited"] == pytest.approx(2 / 3)
    assert 0.5 < out["realised_depth_ratio_edited"]["median"] < 1.5  # the realised depth moves little (no cell number)
    # per category too: the edited category's depth is reported, and blood (held here) has no edited cell
    assert out["per_category"]["cm"]["max_abs_expected_ratio_minus_1"] < 1e-9
    assert 0.5 < out["per_category"]["cm"]["realised_depth_ratio"]["median"] < 1.5
    assert out["per_category"]["blood"]["realised_depth_ratio"]["median"] is None


def test_run_n21_on_synthetic_stages_reports_no_cell_number(vm, monkeypatch, tmp_path):
    stages = _n21_stages()

    def load(path, lineage=False):
        for k, a in stages.items():
            if f"/{k}_RNA" in f"/{Path(path).name}":
                return a.copy()
        raise AssertionError(f"unexpected read {path}")

    monkeypatch.setattr(vm, "load_stage", load)
    monkeypatch.setattr(vm.t1_scope, "held_flat_genes", lambda *a, **k: list(HELD))
    monkeypatch.setattr(vm.hygiene, "load_union", lambda *a, **k: UNION)
    monkeypatch.setattr(vm, "sha256_file", lambda path: "0" * 64)
    vm.main(["n21", "--out", "n21", "--data", str(tmp_path)])
    rec = json.loads((vm.RESULTS / "n21" / "n21.json").read_text())
    assert set(rec["test_pair"]) >= {"cm", "endocardium", "blood"} and "cm" in rec["categories_in_both"]
    cm = rec["test_pair"]["cm"]
    assert cm["factor_size"]["n_metacells"] >= 5 and 0.0 <= cm["factor_size"]["share_genes_touched"] <= 1.0
    assert cm["share_of_category_cells_touched"] > 0 and "concordance_with_backtest" in cm
    assert 0.0 <= cm["concordance_with_backtest"]["axis_abs_cosine"] <= 1.0 and "split_half_ceiling" in cm
    assert rec["depth"]["max_abs_expected_ratio_minus_1"] < 1e-9 and rec["depth"]["realised_depth_identical_when_held"]
    assert rec["depth"]["per_category"]["cm"]["realised_depth_ratio"]["median"] is not None
    # blood's cells are R1's held labels: none is edited, yet its factors are fitted and reported
    assert rec["test_pair"]["blood"]["share_of_category_cells_touched"] == 0.0 and "Blood" in rec["r1_held_labels"]
    # the additions of the review: the amplification share, the projected concordance, the no-development reference
    # and the share on R2's union genes, per category (both fits' amplification where the category is in both)
    assert set(cm["amplification"]) == {"test", "backtest"} and np.isfinite(cm["amplification"]["test"]["r2"])
    assert 0.0 <= cm["concordance_after_projecting_amplification"]["axis_abs_cosine"] <= 1.0
    assert cm["no_development_reference"] is not None and "size_spearman" in cm["no_development_reference"]
    assert cm["no_development_size_ratio"]["test"] > 0 and cm["no_development_size_ratio"]["backtest"] > 0
    assert cm["r2_union_share"]["test"]["union_genes_in_fit"] == 2  # g5 and g7 are free genes here (held_a is not held)
    assert 0.0 <= cm["r2_union_share"]["test"]["l2_share"] <= 1.0
    assert rec["settings"]["no_development_streams"] == {"test": 2, "backtest": 3}
    text = (vm.RESULTS / "n21" / "n21.txt").read_text()
    assert "Decides nothing" in text and "split-half" in text
    assert "half-data estimates" in text and "no-development reference" in text
    assert "amplification (the finished factors' regression" in text and "R2's union genes" in text
    assert "n_cells" not in json.dumps(rec) and "truth_cells" not in json.dumps(rec) and "n_obs" not in text
    assert rec["stamp"]["commit"] == "abc123" and rec["stamp"]["tree_clean"] is True
    assert any(k.endswith("hygiene_genes.csv") for k in rec["stamp"]["inputs"])  # the union file is an input
    with pytest.raises(SystemExit, match="already holds"):
        vm.main(["n21", "--out", "n21", "--data", str(tmp_path)])


def test_n21_without_the_union_file_reports_no_union_share(vm, monkeypatch, tmp_path):
    stages = _n21_stages(seed=29)

    def load(path, lineage=False):
        for k, a in stages.items():
            if f"/{k}_RNA" in f"/{Path(path).name}":
                return a.copy()
        raise AssertionError(f"unexpected read {path}")

    def refuse(*a, **k):
        raise vm.hygiene.Refused("no union file")

    monkeypatch.setattr(vm, "load_stage", load)
    monkeypatch.setattr(vm.t1_scope, "held_flat_genes", lambda *a, **k: list(HELD))
    monkeypatch.setattr(vm.hygiene, "load_union", refuse)
    monkeypatch.setattr(vm, "sha256_file", lambda path: "0" * 64)
    vm.main(["n21", "--out", "n21b", "--data", str(tmp_path)])
    rec = json.loads((vm.RESULTS / "n21b" / "n21.json").read_text())
    assert rec["test_pair"]["cm"]["r2_union_share"] is None
    assert "not available" in (vm.RESULTS / "n21b" / "n21.txt").read_text()


def test_n21_fit_is_the_candidates_fit(vm):
    """N21 fits through harness.g_mom_pair_factors, the function the test build's --g-mom calls: pinned to the
    candidate's own fold factors above; here it needs only the shape of its result."""
    stages = _n21_stages(seed=21)
    f, raw = vm.harness.g_mom_pair_factors(stages["E10.5"], "E10.5_RNA", stages["E9.5"], "E9.5_RNA")
    assert f.metacells_on == "later" and "cm" in f.categories
    cm = f.categories["cm"]
    assert np.abs(np.average(cm.logf, axis=0, weights=cm.sizes)).mean() < 0.3 * np.abs(cm.logf).mean()


# ------------------------------------------------------------------------------------------------ the T1 test build
def test_g_mom_option_is_off_by_default_and_off_reproduces_todays_file(bt, t1_env, tmp_path, capsys):  # noqa: F811
    """--g-mom is off by default: the file is byte-for-byte the one build_t1 wrote before this option (R1 on, R2 off
    here, so that no R2 input is needed), its record has no g_mom field, E9.5 is not read, and the cells the record
    lists are unchanged."""
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "runs/x/T1_s1.h5ad"]).g_mom is False
    assert ap.parse_args(["t1", "--g-mom", "--out", "runs/x/T1_s1.h5ad"]).g_mom is True
    today = _reference_file(bt, t1_env, tmp_path / "today.h5ad")
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--out", str(out)])
    assert out.read_bytes() == today.read_bytes()
    assert t1_env.reads == ["E10.5", ("Feng", list(T1_GENES))]  # no E9.5 read when the option is off
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert "g_mom" not in rec


def _g_mom_loads(bt, t1_env, monkeypatch):  # noqa: F811
    """The E9.5 file too, for runs with --g-mom (the fixture's load allows only E10.5)."""
    calls = []

    def load(path, **kw):
        calls.append(Path(path).name)
        return (t1_env.ref if Path(path).name.startswith("E10.5") else t1_env.start).copy()

    monkeypatch.setattr(bt, "load_stage", load)
    return calls


def test_g_mom_on_adds_the_factors_after_fengs_direction_and_keeps_r1(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    today = _reference_file(bt, t1_env, tmp_path / "today.h5ad")
    calls = _g_mom_loads(bt, t1_env, monkeypatch)
    # record what the hook receives: the population after Feng's direction (the same as with the option off)
    received = []
    real = bt.t1_scope.build_submission

    def spy(*a, **k):
        hook = k.get("post_cells")

        def wrapped(pop):
            received.append((pop.copy(), np.array(pop.layers["counts"].todense())))
            return hook(pop)

        return real(*a, **{**k, "post_cells": wrapped}) if hook else real(*a, **k)

    monkeypatch.setattr(bt.t1_scope, "build_submission", spy)
    out = tmp_path / "T1_s1_gmom.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--out", str(out)])
    assert calls == ["E10.5_RNA.h5ad", "E9.5_RNA.h5ad"]  # the E9.5 file is read for the fit, once
    on, off = ad.read_h5ad(out), ad.read_h5ad(today)
    assert on.shape == off.shape and list(on.var_names) == list(off.var_names)
    assert not np.array_equal(np.asarray(on.X), np.asarray(off.X))
    pop_before, counts_before = received[0]
    labels = pop_before.obs["celltype"].astype(str).to_numpy()
    held_labels = sorted(ts.r1_unbracketed_labels(np.unique(labels)))
    held = np.isin(labels, held_labels)
    assert held.any() and (~held).any() and {"Blood", "Macrophage", "Early VIC"} & set(held_labels)
    x_on, x_off = np.asarray(on.X), np.asarray(off.X)
    assert np.array_equal(x_on[held], x_off[held])  # R1's cells keep exactly what the option-off build gives them
    assert any(not np.array_equal(x_on[i], x_off[i]) for i in np.flatnonzero(~held))  # the others are edited
    # the record names the rule and no cell number; shares, count and scope are the off build's
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["g_mom"]["categories_with_factors"] and rec["g_mom"]["stream"] == [bt.harness.BUILD_SEED_OFFSET + 1, 5]
    text = json.dumps(rec["g_mom"])
    assert "n_cells" not in text and "cells_a" not in text and rec["n_cells"] == "max_cells"
    assert t1_env.held_calls == [T1_HELD]  # the pre-registered list goes to the build, as before (R2 off here)
    own = bt.harness.t1_chamber_shares(bt.harness._shallow(t1_env.ref), "E10.5_RNA")
    assert rec["shares"] == pytest.approx(own)  # the mix is E10.5's own chamber mix: the option touches no share
    assert counts_before.shape[0] == x_on.shape[0]  # the cell number is the option-off build's: no cell added


def test_g_mom_on_keeps_each_cells_expected_depth_and_the_shares(bt, t1_env, monkeypatch):  # noqa: F811
    _g_mom_loads(bt, t1_env, monkeypatch)
    ref = t1_env.ref.copy()
    cats = [c for c, s in bt.harness.t1_chamber_shares(bt.harness._shallow(ref), "E10.5_RNA").items() if s > 0]
    covered = {lab for lab, g in bt.feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique()).items()
               if g in ("CM", "Endocardium", "Endothelium", "Epicardium")}
    hook, rec = bt._t1_g_mom(ref, cats, 1000001, covered)
    assert rec["categories_with_factors"] and set(rec["categories_with_factors"]) <= set(cats)
    # a population as build_submission would hand it over: real cells with obs and counts
    t1_scope.classify(ref, source="E10.5_RNA")
    pop, _ = t1_scope.reweight(ref, {c: 1 / len(cats) for c in cats}, 200, 5,
                               fallback=t1_scope.load_spec().fallbacks("E10.5_RNA"))
    out = hook(pop)
    assert out.n_obs == pop.n_obs and list(out.var_names) == list(pop.var_names)
    before = np.asarray(pop.layers["counts"].sum(1)).ravel().astype(np.float64)
    after = np.asarray(sparse.csr_matrix(out.layers["counts"]).sum(1)).ravel().astype(np.float64)
    labels = pop.obs["celltype"].astype(str).to_numpy()
    held = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    edited = ~held & np.isin(pop.obs["category"].astype(str).to_numpy(), rec["categories_with_factors"])
    assert edited.any() and np.array_equal(before[held], after[held])
    assert np.abs(after[edited] / before[edited] - 1.0).mean() < 0.15  # expected depth kept: a noise-sized change
    assert rec["median_depth_scale"] > 0 and "n_edited" not in rec


def test_the_hook_refuses_a_label_the_builds_feng_direction_does_not_cover(bt, t1_env, monkeypatch):  # noqa: F811
    """The hook's own R1 mask and the build's Feng-covered labels are derived separately; a label edited outside the
    covered set stops the build rather than editing a label Feng's direction never moved."""
    _g_mom_loads(bt, t1_env, monkeypatch)
    ref = t1_env.ref.copy()
    cats = [c for c, s in bt.harness.t1_chamber_shares(bt.harness._shallow(ref), "E10.5_RNA").items() if s > 0]
    hook, _ = bt._t1_g_mom(ref, cats, 1000001, {"V-CM"})  # only the ventricular CM label is covered here
    t1_scope.classify(ref, source="E10.5_RNA")
    pop, _ = t1_scope.reweight(ref, {c: 1 / len(cats) for c in cats}, 200, 5,
                               fallback=t1_scope.load_spec().fallbacks("E10.5_RNA"))
    assert set(pop.obs["celltype"].astype(str)) - {"V-CM"}  # the draw holds other edited labels
    with pytest.raises(ValueError, match="does not cover"):
        hook(pop)


def test_g_mom_with_r2_keeps_the_candidates_held_flat_list(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    """With R2 on (the default), Feng's direction holds the union and g_mom's step only the pre-registered list, as the
    candidate's row: the option changes neither R1 nor R2, and adds its own list to neither."""
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: UNION)
    _g_mom_loads(bt, t1_env, monkeypatch)
    seen = []
    real = ts.apply_state_factors

    def spy(pop, sf, held_flat, seed, **kw):
        seen.append((list(held_flat), seed, kw["rescale"]))
        return real(pop, sf, held_flat, seed, **kw)

    monkeypatch.setattr(ts, "apply_state_factors", spy)
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", "--g-mom", "--out", str(out)])
    assert seen == [(list(T1_HELD), [bt.harness.BUILD_SEED_OFFSET + 1, 5], False)]  # one edit, the candidate's stream
    union_held = sorted(set(T1_HELD) | set(UNION["genes"]))
    assert t1_env.held_calls == [union_held]  # Feng's direction keeps R2's union, as without the option
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["hold_unbracketed"] is True and rec["feng_hygiene"]["genes"] == len(UNION["genes"])
    assert "pre-registered" in rec["g_mom"]["held_flat"]


def test_g_mom_is_refused_without_r1(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    monkeypatch.setattr(bt, "load_stage", lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    with pytest.raises(SystemExit, match="--g-mom is refused with --no-hold-unbracketed"):
        bt.main(["t1", "--g-mom", "--no-hold-unbracketed", "--no-feng-hygiene", "--out", str(tmp_path / "T1_s1.h5ad")])
    assert not (tmp_path / "T1_s1.h5ad").exists()


def test_the_hook_runs_last_and_sees_the_direction_applied_on_build_submission():
    """t1_scope.build_submission's new post_cells: the identity hook changes nothing, a hook is applied after the
    directions and before the panel is taken."""
    ref, start, _ = _small_inputs(seed=41)
    board = Board(key="pseudo:x", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600)
    shares = harness.t1_chamber_shares(harness._shallow(ref), "E9.5_RNA")
    groups = {"V-CM": "CM", "Endocardium": "Endocardium"}
    rng = np.random.default_rng(0)
    direction = np.zeros(len(GENES))
    direction[GENES.index("g0")] = 4.0  # a large direction on one gene, so its effect cannot be missed
    dirs = {"CM": direction, "Endocardium": rng.normal(0.0, 0.1, len(GENES))}
    kw = dict(start=None, held_flat=[], ref_source="E9.5_RNA", directions=dirs, direction_fraction=0.5,
              direction_groups=groups)
    plain, _ = t1_scope.build_submission(harness._shallow(ref), board, shares, 400, 1000011, **kw)
    same, _ = t1_scope.build_submission(harness._shallow(ref), board, shares, 400, 1000011,
                                        post_cells=lambda p: p, **kw)
    assert np.array_equal(plain.X, same.X)
    seen = {}

    def hook(pop):
        seen["pop"] = pop
        out = pop.copy()
        C = sparse.csr_matrix(pop.layers["counts"], dtype=np.int32).toarray()
        C[:, GENES.index("g1")] = 0
        out.layers["counts"] = sparse.csr_matrix(C)
        out.X = log1p_cp10k(sparse.csr_matrix(C))
        return out

    edited, _ = t1_scope.build_submission(harness._shallow(ref), board, shares, 400, 1000011, post_cells=hook, **kw)
    assert (edited.X[:, GENES.index("g1")] == 0).all() and not np.array_equal(edited.X, plain.X)
    pop = seen["pop"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    # the hook saw Feng's direction already applied: the V-CM cells' g0 went up from the drawn cells' g0
    classified_ref = harness._shallow(ref)
    t1_scope.classify(classified_ref, source="E9.5_RNA")
    drawn, _ = t1_scope.reweight(classified_ref, shares, 400, 1000011, label_key="celltype",
                                 fallback=t1_scope.load_spec().fallbacks("E9.5_RNA"))
    i = GENES.index("g0")
    cm = labels == "V-CM"
    before = sparse.csr_matrix(drawn.layers["counts"]).toarray()[cm][:, i].sum()
    after = sparse.csr_matrix(pop.layers["counts"]).toarray()[cm][:, i].sum()
    assert after > 2 * max(before, 1)
