"""experiments/heart_profile/xpc_table.py (row H36; ruling row R56): the XPC run-profile table, on synthetic arrays
only.

Synthetic groups: three categories at E8.5, E9.5 and E10.5 on two instruments, Poisson counts from shares that are a
category's base times a stage-and-category biology both instruments see, times each instrument's own per-gene
efficiency (a constant platform term), with or without a run profile planted in E10.5's MERFISH run. A planted shared
profile is recovered and passes every check; a lineage-specific shift is not (the lineage check stops it); constant
platform terms cancel and pure Poisson draws read STOP; the practice profile reads no E10.5 group; the entry, the
halves, the floor, the clip, the loader's categories and libraries, and load_profile's refusals. The run-profile step
itself (vec.t2_build.apply_run_profile): a zero profile changes nothing and draws nothing, depth is kept. Nothing is
read from data/."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
G = 200
GENES = [f"g{i}" for i in range(G)]


@pytest.fixture(scope="module")
def xpc():
    spec = importlib.util.spec_from_file_location("xpc_table_for_tests",
                                                  REPO / "experiments" / "heart_profile" / "xpc_table.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def synth(xpc, seed=0, n=600, depth=400.0, profile=None, lineage=None, platform=1.0, biology=1.0, stages=None,
          e95=None):
    """Groups (stage, instrument, category) -> CSR counts. The biology of a stage and category is the same on both
    instruments; each instrument has its own constant per-gene efficiency; `profile` is planted in E10.5's MERFISH run
    for every category, `e95` in E9.5's, `lineage` in E10.5's MERFISH cm cells only."""
    rng = np.random.default_rng(seed)
    base = {k: rng.dirichlet(np.full(G, 3.0)) for k in xpc.CATEGORIES}
    eff = {xpc.MERFISH: np.exp(platform * rng.normal(0, 1, G)), xpc.RNA: np.exp(platform * rng.normal(0, 1, G))}
    gs = {}
    for s in stages or xpc.STAGES:
        for k in xpc.CATEGORIES:
            b = base[k] * np.exp(biology * rng.normal(0, 0.5, G))
            for inst in xpc.INSTRUMENTS:
                w = b * eff[inst]
                if inst == xpc.MERFISH and s == "E10.5":
                    if profile is not None:
                        w = w * np.exp(profile)
                    if lineage is not None and k == "cm":
                        w = w * np.exp(lineage)
                if inst == xpc.MERFISH and s == "E9.5" and e95 is not None:
                    w = w * np.exp(e95)
                w = w / w.sum()
                C = rng.poisson(depth * w, size=(n, G))
                gs[(s, inst, k)] = sparse.csr_matrix(C[C.sum(1) > 0])
    return gs, eff


# ------------------------------------------------------------------------------------------------ the profile
def test_a_planted_shared_profile_is_recovered_and_passes_every_check(xpc):
    p = np.random.default_rng(5).normal(0, 0.5, G)
    gs, _ = synth(xpc, profile=p)
    res = xpc.measure(gs, xpc.TEST, alternative=gs)
    e = res["eligible"]
    assert e.sum() > G // 2
    assert xpc.spearman(res["delta"][e], p[e]) > 0.9
    # the planted profile less one constant (each run's normalisation): its spread is the planted one's
    resid = res["delta"][e] - p[e]
    assert np.std(resid) < 0.15 and np.std(res["delta"][e]) > 0.35
    assert res["passed"] and res["failed"] == []
    c = res["checks"]
    assert c["noise"]["sd"] >= 2 * c["noise"]["sigma_n"] and c["halves"]["spearman"] >= 0.5
    assert c["lineages"]["median"] >= 0.4 and c["library"]["spearman"] == pytest.approx(1.0)
    assert np.all(res["delta"][~e] == 0.0) and np.abs(res["delta"]).max() <= math.log(3) + 1e-15


def test_a_lineage_specific_shift_is_not_a_run_profile(xpc):
    shift = np.random.default_rng(6).normal(0, 0.6, G)
    gs, _ = synth(xpc, lineage=shift)
    res = xpc.measure(gs, xpc.TEST, alternative=gs)
    e = res["eligible"]
    assert xpc.spearman(res["by_category"]["cm"][e], shift[e]) > 0.9  # it is there, in cm alone
    assert res["checks"]["lineages"]["median"] < 0.4
    assert "lineages" in res["failed"] and not res["passed"]
    # the median over the categories does not carry it: what is left is the size of the other categories' noise
    assert np.std(res["delta"][e]) < 0.3 * np.std(shift[e])


def test_constant_platform_terms_cancel_and_read_stop(xpc):
    gs, eff = synth(xpc, platform=1.5, biology=1.0, depth=4000.0)
    res = xpc.measure(gs, xpc.TEST, alternative=gs)
    e = res["eligible"]
    platform = np.log(eff[xpc.MERFISH]) - np.log(eff[xpc.RNA])
    assert np.std(platform[e]) > 1.5
    assert np.std(res["delta"][e]) < 0.1  # what is left is a constant and the draws' noise
    assert not res["passed"] and "noise" in res["failed"]


def test_pure_poisson_draws_read_stop(xpc):
    for seed in (1, 2, 3):
        gs, _ = synth(xpc, seed=seed, platform=0.0, biology=0.0)
        for kind in xpc.KINDS:
            use = {g: v for g, v in gs.items() if g[0] in xpc.KIND_STAGES[kind]}
            res = xpc.measure(use, kind, alternative=use if kind == xpc.TEST else None)
            assert not res["passed"] and {"noise", "halves"} <= set(res["failed"]), (seed, kind, res["failed"])


def test_the_formulas(xpc):
    gs, _ = synth(xpc, seed=9, profile=np.linspace(-2.0, 2.0, G))
    lv = {g: xpc.log_level(C) for g, C in gs.items()}
    D = {(s, k): lv[(s, xpc.MERFISH, k)] - lv[(s, xpc.RNA, k)] for s in xpc.STAGES for k in xpc.CATEGORIES}
    t = xpc.category_profiles(gs, xpc.TEST)
    p = xpc.category_profiles(gs, xpc.PRACTICE)
    for k in xpc.CATEGORIES:
        np.testing.assert_allclose(t[k], D[("E10.5", k)] - (D[("E8.5", k)] + D[("E9.5", k)]) / 2, atol=1e-12)
        np.testing.assert_allclose(p[k], D[("E9.5", k)] - D[("E8.5", k)], atol=1e-12)
    C = gs[("E8.5", xpc.MERFISH, "cm")]
    share = np.asarray((sparse.diags(1.0 / np.asarray(C.sum(1)).ravel()) @ C).mean(0)).ravel()
    np.testing.assert_allclose(xpc.log_level(C), np.log(share + 1e-5), atol=1e-12)
    res = xpc.profile(gs, xpc.TEST)
    med = np.median(np.vstack([t[k] for k in xpc.CATEGORIES]), axis=0)
    e = res["eligible"]
    np.testing.assert_array_equal(res["delta"], np.clip(np.where(e, med, 0.0), -math.log(3), math.log(3)))
    assert (np.abs(med[e]) > math.log(3)).any()  # the clip acts here
    lvl = xpc.reference_cp10k(gs, xpc.TEST)
    np.testing.assert_array_equal(e, (lvl[xpc.MERFISH] >= 0.5) & (lvl[xpc.RNA] >= 0.5))


def test_the_floor_zeroes_genes_under_half_a_cp10k(xpc):
    gs, _ = synth(xpc, seed=4, profile=np.random.default_rng(1).normal(0, 0.5, G))
    gs = {g: (C.tolil() if g[0] == "E10.5" and g[1] == xpc.RNA else C) for g, C in gs.items()}
    for g in list(gs):
        if g[0] == "E10.5" and g[1] == xpc.RNA:
            gs[g][:, :5] = 0  # five genes never seen by the RNA at E10.5
            gs[g] = sparse.csr_matrix(gs[g])
    res = xpc.profile(gs, xpc.TEST)
    assert not res["eligible"][:5].any() and np.all(res["delta"][:5] == 0.0)
    assert xpc.profile(gs, xpc.PRACTICE)["eligible"][:5].all()  # the practice floor is E9.5's: no E10.5 group


def test_the_practice_profile_reads_no_e105_group(xpc):
    gs, _ = synth(xpc, stages=("E8.5", "E9.5"), profile=None)
    res = xpc.measure(gs, xpc.PRACTICE)
    assert res["entry"]["ok"] and set(res["checks"]) == {"noise", "halves", "lineages"}
    test = xpc.measure(gs, xpc.TEST, alternative=gs)  # the test profile needs E10.5: it stops at its entry
    assert not test["passed"] and test["failed"][0].startswith("entry: E10.5/merfish/cm")


def test_the_entry_wants_200_cells_in_every_group(xpc):
    gs, _ = synth(xpc, n=300)
    gs[("E10.5", xpc.RNA, "progenitor")] = gs[("E10.5", xpc.RNA, "progenitor")][:199]
    ent = xpc.entry(gs)
    assert not ent["ok"] and ent["short"] == ["E10.5/rna/progenitor"]
    res = xpc.measure(gs, xpc.TEST, alternative=gs)
    assert not res["passed"] and res["failed"] == ["entry: E10.5/rna/progenitor"] and "delta" not in res
    del gs[("E8.5", xpc.MERFISH, "cm")]
    assert xpc.entry(gs)["short"] == ["E10.5/rna/progenitor", "E8.5/merfish/cm"]


def test_the_halves_are_disjoint_fixed_and_one_stream_per_group(xpc):
    gs, _ = synth(xpc, n=301)
    a, b = xpc.halves(gs)
    a2, _ = xpc.halves(gs)
    one = {g: gs[g] for g in list(gs)[:3]}
    a3, _ = xpc.halves(one)
    for g in gs:
        n = gs[g].shape[0]
        assert a[g].shape[0] == n // 2 and b[g].shape[0] == n - n // 2
        assert (a[g] != a2[g]).nnz == 0
        both = sparse.vstack([a[g], b[g]])
        assert both.sum() == gs[g].sum()
    for g in one:
        assert (a[g] != a3[g]).nnz == 0  # a group's split does not depend on the other groups


def test_spearman(xpc):
    assert xpc.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert xpc.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert math.isnan(xpc.spearman([1, 1, 1], [1, 2, 3])) and math.isnan(xpc.spearman([1], [2]))


# ------------------------------------------------------------------------------------------------ the files
def _merfish(labels, rng):
    n = len(labels)
    C = rng.poisson(3.0, size=(n, G)).astype(np.float32)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / C.sum(1, keepdims=True))).astype(np.float32),
                   obs={"celltype": np.asarray(labels)})
    a.var_names = GENES
    a.obs_names = [f"m{i}" for i in range(n)]
    a.layers["counts"] = sparse.csr_matrix(C)
    return a


def _rna(labels, libs, rng, zero=()):
    n = len(labels)
    genes = [*GENES, *[f"x{i}" for i in range(20)]]
    C = rng.poisson(5.0, size=(n, len(genes))).astype(np.float32)
    for i in zero:  # a cell with no panel count
        C[i, :G] = 0
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / C.sum(1, keepdims=True))).astype(np.float32),
                   obs={"celltype": np.asarray(labels)})
    a.var_names = genes
    a.obs_names = [f"c{i}_{lib}" for i, lib in enumerate(libs)]
    a.layers["counts"] = sparse.csr_matrix(C)
    return a


def test_the_files_are_read_by_both_chamber_configs_with_library_3_at_e105(xpc, tmp_path):
    rng = np.random.default_rng(0)
    labels = np.array(["V-CM", "Endocardium", "aSHF", "Blood"] * 6)
    files, calls = {}, []
    for s in xpc.STAGES:
        files[xpc.FILES[xpc.MERFISH][s]] = _merfish(labels, rng)
        libs = (["1", "2", "3"] * 8) if s == "E10.5" else (["1", "2"] * 12)
        files[xpc.FILES[xpc.RNA][s]] = _rna(labels, libs, rng, zero=(0,) if s == "E9.5" else ())

    def load(path, genes=None, drop_duplicates=False, lineage=True):
        rel = str(Path(path).relative_to(tmp_path))
        calls.append((rel, None if genes is None else len(genes), drop_duplicates, lineage))
        return files[rel].copy()

    raw = xpc.read_files(tmp_path, xpc.STAGES, GENES, load=load)
    assert sorted(calls) == sorted([(xpc.FILES[xpc.MERFISH][s], G, True, True) for s in xpc.STAGES]
                                   + [(xpc.FILES[xpc.RNA][s], None, False, False) for s in xpc.STAGES])
    g3 = xpc.groups(raw, xpc.STAGES)
    g123 = xpc.groups(raw, xpc.STAGES, xpc.ALL_E105_LIBRARIES)
    for k, lab in zip(xpc.CATEGORIES, ("V-CM", "Endocardium", "aSHF"), strict=True):
        for s in xpc.STAGES:
            assert g3[(s, xpc.MERFISH, k)].shape == (6, G)
        rna = files[xpc.FILES[xpc.RNA]["E10.5"]]
        rows = np.flatnonzero((np.asarray(rna.obs["celltype"]) == lab)
                              & np.array([n.endswith("_3") for n in rna.obs_names]))
        np.testing.assert_array_equal(g3[("E10.5", xpc.RNA, k)].toarray(),
                                      np.asarray(rna.layers["counts"][rows][:, :G].todense(), dtype=np.int64))
        assert g123[("E10.5", xpc.RNA, k)].shape[0] == 6 and g3[("E10.5", xpc.RNA, k)].shape[0] == 2
        assert g3[("E8.5", xpc.RNA, k)].shape[0] == 6
    assert g3[("E9.5", xpc.RNA, "cm")].shape[0] == 5  # its first cell has no panel count: left out
    assert xpc.cells(g3)["E10.5/rna/cm"] == 2
    bad = _rna(labels, ["1", "x"] * 12, rng)
    with pytest.raises(ValueError):
        xpc.libraries_of(bad.obs_names)


# ------------------------------------------------------------------------------------------------ the record
def _profile_file(xpc, tmp_path, **over) -> tuple[Path, str]:
    gs, _ = synth(xpc, profile=np.random.default_rng(5).normal(0, 0.5, G))
    res = xpc.measure(gs, xpc.TEST, alternative=gs)
    rec = {**xpc.record(res, GENES, {"x": "y"}), **over}
    text = xpc.dumps(rec)
    path = tmp_path / "xpc_profile.json"
    path.write_text(text)
    return path, hashlib.sha256(text.encode()).hexdigest()


def test_the_record_and_load_profile(xpc, tmp_path):
    path, sha = _profile_file(xpc, tmp_path)
    rec = json.loads(path.read_text())
    assert rec["axis"] == xpc.AXIS and rec["kind"] == "test" and rec["passed"] is True and rec["genes"] == GENES
    assert set(rec) >= {"delta", "eligible", "level", "by_category", "halves", "checks", "rule", "constants"}
    d = xpc.load_profile(path, sha, GENES)
    np.testing.assert_array_equal(d, np.asarray(rec["delta"]))
    assert xpc.sha256_text(path.read_text()) == sha
    with pytest.raises(ValueError, match="sha256"):
        xpc.load_profile(path, "0" * 64, GENES)
    with pytest.raises(ValueError, match="order"):
        xpc.load_profile(path, sha, GENES[::-1])
    with pytest.raises(ValueError, match="passed practice"):
        xpc.load_profile(path, sha, GENES, kind=xpc.PRACTICE)
    for over, msg in (({"passed": False}, "passed"), ({"kind": "practice"}, "passed test"),
                      ({"delta": [2.0] * G}, "within log 3"), ({"delta": [float("nan")] * G}, "finite"),
                      ({"delta": [0.0] * (G - 1)}, "finite"), ({"axis": "other"}, "passed test")):
        p, s = _profile_file(xpc, tmp_path, **over)
        with pytest.raises(ValueError, match=msg):
            xpc.load_profile(p, s, GENES)


# ------------------------------------------------------------------------------------------------ the step
def test_apply_run_profile_changes_nothing_at_zero_and_keeps_depth():
    from vec import t2_build
    from vec.prep import log1p_cp10k

    rng = np.random.default_rng(3)
    C = sparse.csr_matrix(rng.poisson(2.0, size=(400, G)))
    new, X, rec = t2_build.apply_run_profile(C, np.zeros(G), [1000001, 37])
    assert (new != C).nnz == 0 and X.dtype == np.float32
    np.testing.assert_array_equal(X, log1p_cp10k(C).toarray())
    assert rec["stream"] == 37 and set(rec) == {"rule", "stream", "median_depth_scale"}  # no cell number
    d = rng.normal(0, 0.5, G)
    new, X, rec = t2_build.apply_run_profile(C, d, [1000001, 37])
    again, X2, _ = t2_build.apply_run_profile(C, d, [1000001, 37])
    assert (new != again).nnz == 0 and np.array_equal(X, X2)
    lib0, lib1 = np.asarray(C.sum(1)).ravel(), np.asarray(new.sum(1)).ravel()
    assert abs(lib1.mean() / lib0.mean() - 1) < 0.02  # expected depth kept
    up, down = d > 0.4, d < -0.4
    ratio = (np.asarray(new.sum(0)).ravel() + 1) / (np.asarray(C.sum(0)).ravel() + 1)
    assert ratio[down].mean() > 1.0 > ratio[up].mean()  # exp(-delta): genes the run over-detected go down
    with pytest.raises(ValueError):
        t2_build.apply_run_profile(C, np.zeros(G - 1), [1, 37])
    with pytest.raises(ValueError):
        t2_build.apply_run_profile(C, np.full(G, np.inf), [1, 37])
