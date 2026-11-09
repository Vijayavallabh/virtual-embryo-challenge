"""The g_mom verification (tries T20 to T22 and diagnostic N21 of experiments/harness/TRIES.md) on synthetic data: the
placebo factors, their harness arms, the pair fit shared by N21 and the test build, the placebo readout, the
composition arms at the truth's own chamber mix, N21 and the --g-mom option of the T1 test build."""

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
from test_t1_scope import GENES, HELD
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


def test_the_placebo_arms_build_inside_the_read_guard(monkeypatch, tmp_path):
    """The harness's own run of g_mom and its placebos on a T1 pseudo-task (one seed): nothing a candidate must not read
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
    names = ["b0", "g_mom", "g_mom_flip", "g_mom_perm"]
    res = harness.run_task(task, [CANDIDATES[n] for n in names], seeds=(0,), data_dir=tmp_path,
                           boards={"T1:val": replace(board, key="T1:val")}, log=None)
    row = res["rows"][0]
    assert res["candidates"] == ["floor", *names] and row["scores"]["floor"] == 50.0
    for name in names:
        assert np.isfinite(row["scores"][name]) and row["n_cells"][name] == 400
    for name in ("g_mom", "g_mom_flip", "g_mom_perm"):
        assert not any("results" in p or p.startswith(("data/", "runs/")) for p in row["reads"][name])


def test_the_placebo_candidates_apply_where_g_mom_does():
    tasks = {t.key: t for ts_ in harness.load_tasks().values() for t in ts_}
    for name in ("g_mom_flip", "g_mom_perm"):
        c = CANDIDATES[name]
        assert c.applies(tasks["T1/E9.5->E10.5"]) and c.applies(tasks["T1/E9.5->E10.5:chamber"])
        assert not c.applies(tasks["T1/E8.5->E9.5"]) and not c.applies(tasks["T2h/<=E9.5->E10.5"])
        assert c.params["scope"] == CANDIDATES["g_mom"].params["scope"] and c.derived == ("t1_held_flat",)
    assert CANDIDATES["g_mom_perm"].params["permutation_seed"] == harness.G_MOM_PERM_SEED


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
def _run_rows(b0, g_mom, flip, perm, seeds=range(10), t11=None):
    """Per-seed rows of a harness run: scores and group skills per arm; t11: scores for b0 and g_mom to use instead."""
    rows = []
    for i, s in enumerate(seeds):
        scores = {"floor": 50.0, "b0": b0[i], "g_mom": g_mom[i], "g_mom_flip": flip[i], "g_mom_perm": perm[i]}
        if t11 is not None:
            scores["b0"], scores["g_mom"] = t11
        groups = {arm: {"de_recovery": 0.6, "direction": 0.6, "distribution": 0.7 + (v - 50) / 100,
                        "covariance": 0.7 + (v - 50) / 200} for arm, v in scores.items()}
        rows.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s, "scores": scores, "groups": groups})
    return rows


def _write_run(tmp_path, folder, rows_by_cut, sha):
    d = tmp_path / folder
    d.mkdir(parents=True)
    for key, name in (("T1/E9.5->E10.5", "T1_E9.5_E10.5.json"),
                      ("T1/E9.5->E10.5:chamber", "T1_E9.5_E10.5_chamber.json")):
        (d / name).write_text(json.dumps(
            {"task": {"key": key, "test_board": "T1:test"},
             "candidates": ["floor", "b0", "g_mom", "g_mom_flip", "g_mom_perm"], "rows": rows_by_cut[key],
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


def _matching_rows(gain_flip, gain_perm):
    b0 = [49.0 + 0.1 * s for s in range(10)]
    g_mom = [50.0 + 0.1 * s for s in range(10)]
    flip = [g - gain_flip[i] for i, g in enumerate(g_mom)]
    perm = [g - gain_perm[i] for i, g in enumerate(g_mom)]
    return _run_rows(b0, g_mom, flip, perm)


def test_placebo_survives_only_when_each_placebo_is_beaten_on_both_cuts(vm, t11_dir, tmp_path):
    sha = vm.harness._sha256(vm.harness.CONFIG)
    good = {k: _matching_rows([1.0] * 10, [0.8] * 10) for k in vm.CUTS}
    run = _write_run(tmp_path, "good", good, sha)
    vm.main(["placebo", "--runs", str(run), "--out", "p_good"])
    rec = json.loads((vm.RESULTS / "p_good" / "placebo.json").read_text())
    assert rec["survives"] is True and rec["reproduction_of_t11"]["reproduced"] is True
    assert set(rec["placebos"]) == {"g_mom_flip", "g_mom_perm"}
    dec = rec["placebos"]["g_mom_perm"]["decision"]
    assert dec["keep"] and dec["base"] == "g_mom_perm" and dec["change"] == "g_mom"
    assert dec["tasks"]["T1/E9.5->E10.5"]["mean_gain"] == pytest.approx(0.8)
    assert sorted(dec["counted"]) == sorted(vm.CUTS) and "T1/E8.5->E9.5" in dec["excluded"]
    assert rec["placebos"]["g_mom_flip"]["group_margins"]["T1/E9.5->E10.5"]["n_seeds"] == 10
    assert rec["stamp"]["commit"] == "abc123" and rec["stamp"]["tree_clean"] is True
    assert len(rec["stamp"]["inputs"]) == 4
    text = (vm.RESULTS / "p_good" / "placebo.txt").read_text()
    assert "SURVIVES" in text and "T20 g_mom_flip: beaten" in text and "T21 g_mom_perm: beaten" in text
    with pytest.raises(SystemExit, match="already holds"):
        vm.main(["placebo", "--runs", str(run), "--out", "p_good"])


@pytest.mark.parametrize("flip_gain, perm_gain, survives", [
    ([1.0] * 10, [0.1] * 10, False),  # the permutation is not beaten by more than 0.3
    ([0.25] * 10, [1.0] * 10, False),  # the flip is within 0.3 of g_mom
    ([1.0] * 7 + [-0.5] * 3, [1.0] * 10, False),  # mean above 0.3 but only 7 of 10 seeds ahead
    ([1.0] * 10, [0.9] * 10, True),
])
def test_placebo_rule_needs_mean_and_seeds_on_every_cut(vm, t11_dir, tmp_path, flip_gain, perm_gain, survives):
    sha = vm.harness._sha256(vm.harness.CONFIG)
    rows = {k: _matching_rows(flip_gain, perm_gain) for k in vm.CUTS}
    run = _write_run(tmp_path, "r", rows, sha)
    vm.main(["placebo", "--runs", str(run), "--out", "p"])
    rec = json.loads((vm.RESULTS / "p" / "placebo.json").read_text())
    assert rec["survives"] is survives
    # one cut failing is enough: g_mom beaten only on the whole cut
    if survives:
        worse = {**rows, "T1/E9.5->E10.5:chamber": _matching_rows([0.1] * 10, [0.1] * 10)}
        run2 = _write_run(tmp_path, "r2", worse, sha)
        vm.main(["placebo", "--runs", str(run2), "--out", "p2"])
        assert json.loads((vm.RESULTS / "p2" / "placebo.json").read_text())["survives"] is False


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
        r["candidates"] = ["floor", "b0", "g_mom"]
        f.write_text(json.dumps(r))
    with pytest.raises(SystemExit, match="both placebos"):
        vm.main(["placebo", "--runs", str(run2), "--out", "p4"])


# ------------------------------------------------------------------------------------------------ the composition arms
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
    monkeypatch.setattr(vm, "load_stage",
                        lambda path, lineage=False: ref.copy() if "E9.5" in str(path) else start.copy())
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
    assert set(row) >= {"seed", "chamber", "whole"} and row["chamber"]["floor"] == pytest.approx(50.0)
    for key in ("chamber", "whole"):
        assert row[key]["b0"]["score"] != row[key]["g_mom"]["score"]
        assert set(row[key]["g_mom"]["groups"]) == {"de_recovery", "direction", "distribution", "covariance"}
    assert rec["chamber"]["passes"] is False and rec["reproduction_of_oracle_b0"]["reproduced"] is None
    text = (vm.RESULTS / "t22" / "composition.txt").read_text()
    assert "not all ten seeds" in text and "Oracle mix" in text and "Whole cut" in text
    assert rec["stamp"]["commit"] == "abc123" and any("E10.5" in k for k in rec["stamp"]["inputs"])
    with pytest.raises(SystemExit, match="already holds"):
        vm.main(["composition", "--seeds", "0", "--out", "t22", "--data", str(tmp_path)])


def test_composition_b0_arm_is_the_b0_population_and_the_g_mom_arm_edits_it(vm, monkeypatch, tmp_path):
    ref, start, tgt, shares = _t22_io(vm, monkeypatch, tmp_path)
    ctx = harness.BuildContext(task=harness.TaskView.of(next(t for t in harness.load_tasks()["T1:test"]
                                                             if t.key == vm.TASK_CHAMBER)),
                               board=harness._candidate_board(T22_BOARD), inputs={"ref": ref, "start": start},
                               seed=1000003, derived={"t1_held_flat": tuple(HELD)})
    b0, gm = vm.oracle_arms(ctx, shares, HELD)
    # the b0 arm is B0's own population at the mix, as payoff.py builds it through build_submission
    n = harness.count_for(harness._candidate_board(T22_BOARD), ref.n_obs)
    want, _ = t1_scope.build_submission(harness._shallow(ref), T22_BOARD, shares, n, 1000003,
                                        start=harness._shallow(start), damping=harness.B0_DAMPING, held_flat=list(HELD),
                                        ref_source="E9.5_RNA", start_source="E8.5_RNA")
    assert np.array_equal(b0.X, np.asarray(want.X)) and b0.X.shape[0] == n
    assert gm.X.shape == b0.X.shape and not np.array_equal(gm.X, b0.X)
    # the g_mom arm's edit is the population at the oracle mix through g_mom's own fold factors (R1's held labels stay)
    r1, s1 = harness._shallow(ref), harness._shallow(start)
    pop = ts.b0_population(r1, s1, shares, n, 1000003, list(HELD), ref_source="E9.5_RNA", start_source="E8.5_RNA",
                           damping=harness.B0_DAMPING)
    blood = pop.obs["celltype"].astype(str).to_numpy() == "Blood"
    assert blood.any() and np.array_equal(gm.X[blood], b0.X[blood])
    # the fold's factors are fitted on the reference's own chamber categories, as g_mom's
    assert harness._G_MOM_PREP["fold"]["factors"].metacells_on == "later"


def test_composition_is_void_when_the_oracle_base_does_not_reproduce(vm, monkeypatch, tmp_path):
    ref, start, tgt, shares = _t22_io(vm, monkeypatch, tmp_path, committed_mean=12.0, seeds_patch=(0,))
    with pytest.raises(SystemExit, match="does not reproduce the committed t1-grid"):
        vm.main(["composition", "--seeds", "0", "--out", "t22v", "--data", str(tmp_path)])
    rec = json.loads((vm.RESULTS / "t22v" / "composition.json").read_text())
    assert rec["reproduction_of_oracle_b0"]["reproduced"] is False and rec["complete_seeds"] is True


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


def test_expected_depth_ratio_is_one_on_cells_with_free_genes_and_real_counts(vm):
    rng = np.random.default_rng(1)
    C = sparse.csr_matrix(rng.poisson(1.0, size=(100, 30)).astype(np.int64))
    C = C + sparse.csr_matrix(np.eye(100, 30, dtype=np.int64))  # every cell has a count
    F = np.exp(rng.normal(0, 0.5, size=(4, 30)))
    row = rng.integers(0, 4, size=100)
    held = np.zeros(30, dtype=bool)
    held[:5] = True
    F[:, held] = 1.0
    ratio = vm.expected_depth_ratio(C, F, row, held)
    assert np.allclose(ratio, 1.0, atol=1e-12)
    # apply_factors really keeps it in expectation: the realised totals agree on average
    from vec.core import apply_factors
    depth = np.asarray(C.sum(1)).ravel().astype(np.float64)
    many = np.mean([np.asarray(apply_factors(C, F, row, held, keep_depth=True, seed=s)[0].sum(1)).ravel() / depth
                    for s in range(30)], axis=0)
    assert abs(float(many.mean()) - 1.0) < 0.02


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
    text = (vm.RESULTS / "n21" / "n21.txt").read_text()
    assert "Decides nothing" in text and "split-half" in text
    assert "n_cells" not in json.dumps(rec) and "truth_cells" not in json.dumps(rec) and "n_obs" not in text
    assert rec["stamp"]["commit"] == "abc123" and rec["stamp"]["tree_clean"] is True
    with pytest.raises(SystemExit, match="already holds"):
        vm.main(["n21", "--out", "n21", "--data", str(tmp_path)])


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
    hook, rec = bt._t1_g_mom(ref, cats, 1000001)
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
