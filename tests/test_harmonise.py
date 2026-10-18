import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.harmonise import (
    UNASSIGNED,
    UNTRACKED,
    bootstrap_composition,
    bootstrap_difference,
    bootstrap_draws,
    bootstrap_em,
    bracket_difference,
    composition,
    composition_bracket,
    cross_validate,
    default_trims,
    em_priors,
    em_priors_many,
    fit_transfer,
    match_depth,
    median_depth,
    one_hot,
    predict,
    section_ids,
    shared_genes,
    subset_levels,
    thin_counts,
    transfer_metrics,
    trim_mask,
    trim_sensitivity,
    trim_summary,
    z_levels,
    z_location,
    z_position,
    z_profile,
)
from vec.prep import log1p_cp10k
from vec.spatial import ZLevels

N_GENES = 40
# Real lineage-map labels: V-CM and A-CM share the CM group, so group calls can be right when fine calls are not.
LABELS = ["V-CM", "A-CM", "Peri", "NCC", "Neural Tube"]


def _profiles():
    """Per-type gene fractions. V-CM and A-CM share a CM block and differ by a weaker block; type 4 (Neural Tube)
    expresses a block that no other type uses and none of theirs."""
    rng = np.random.default_rng(42)
    p = rng.gamma(0.5, 1.0, size=(5, N_GENES)) + 0.05
    p[:2, 0:6] *= 8.0
    p[0, 6:9] *= 3.0
    p[1, 9:12] *= 3.0
    p[2, 12:18] *= 8.0
    p[3, 18:24] *= 8.0
    p[4, 24:30] *= 12.0
    return p / p.sum(1, keepdims=True)


def make_adata(n, props, seed, depth=80.0, n_sections=6):
    """Poisson counts from the type profiles, with log-normal library sizes, a counts layer and section z levels."""
    rng = np.random.default_rng(seed)
    props = np.asarray(props, dtype=float)
    props = np.concatenate([props, np.zeros(5 - len(props))])
    types = rng.choice(5, size=n, p=props / props.sum())
    lib = rng.lognormal(np.log(depth), 0.3, size=n)
    counts = rng.poisson(_profiles()[types] * lib[:, None])
    counts[np.arange(n), rng.integers(0, N_GENES, n)] += 1  # no empty cells
    C = sparse.csr_matrix(counts.astype(np.int32))
    a = ad.AnnData(X=log1p_cp10k(C))
    a.layers["counts"] = C
    a.var_names = [f"g{i}" for i in range(N_GENES)]
    a.obs_names = [f"s{seed}_c{i}" for i in range(n)]
    a.obs["celltype"] = np.asarray(LABELS, dtype=object)[types]
    z = rng.integers(0, n_sections, n) * 10.0 + rng.normal(0.0, 1e-9, n)
    a.obsm["spatial_3D"] = np.column_stack([rng.normal(size=n), rng.normal(size=n), z])
    return a


@pytest.fixture(scope="module")
def ref():
    return make_adata(800, [0.35, 0.25, 0.2, 0.2], seed=1)


@pytest.fixture(scope="module")
def model(ref):
    return fit_transfer(ref, "celltype", cv_folds=3, n_pcs=10)


def test_auto_classifier_is_picked_by_cv(model):
    assert set(model.cv["method"]) == {"logreg", "knn"}
    assert model.method == model.cv.loc[model.cv["macro_f1"].idxmax(), "method"]
    assert list(model.classes) == sorted(LABELS[:4])
    assert model.group_of["V-CM"] == "CM" and model.group_of["Peri"] == "Pericardium"
    assert list(model.groups) == ["CM", "Pericardium", "NCC"]  # lineage-map order
    assert model.depth == pytest.approx(80.0, rel=0.15)
    np.testing.assert_allclose(model.priors.sum(), 1.0)


def test_predict_recovers_fine_labels_and_groups(model):
    q = make_adata(400, [0.25, 0.25, 0.25, 0.25], seed=2)
    df = predict(model, q)
    assert list(df.index) == list(q.obs_names)
    true = q.obs["celltype"].to_numpy()
    assert (df["label_raw"].to_numpy() == true).mean() > 0.8
    true_group = np.array([model.group_of[t] for t in true])
    assert (df["group_raw"].to_numpy() == true_group).mean() > 0.97
    P = df[[f"p:{c}" for c in model.classes]].to_numpy()
    PG = df[[f"pg:{g}" for g in model.groups]].to_numpy()
    np.testing.assert_allclose(P.sum(1), 1.0, atol=1e-6)
    np.testing.assert_allclose(PG.sum(1), 1.0, atol=1e-6)
    np.testing.assert_allclose(df["pg:CM"], df["p:V-CM"] + df["p:A-CM"], atol=1e-12)
    assert (df["group_conf"] >= df["label_conf"] - 1e-12).all()


def test_threshold_leaves_low_confidence_cells_unassigned(model):
    q = make_adata(200, [0.25, 0.25, 0.25, 0.25], seed=3)
    assert (predict(model, q, threshold=1.01)["label"] == UNASSIGNED).all()
    df = predict(model, q, threshold=0.0)
    assert (df["label"] == df["label_raw"]).all() and (df["group"] == df["group_raw"]).all()
    df = predict(model, q, threshold=0.9, group_threshold=0.0)
    low = df["label_conf"].to_numpy() < 0.9
    assert (df["label"].to_numpy()[low] == UNASSIGNED).all()
    assert (df["label"].to_numpy()[~low] != UNASSIGNED).all()
    assert (df["group"] != UNASSIGNED).all()


def test_query_gene_order_and_missing_genes(model, ref):
    q = make_adata(150, [0.25, 0.25, 0.25, 0.25], seed=4)
    reordered = q[:, list(q.var_names[::-1])].copy()
    pd.testing.assert_frame_equal(predict(model, q), predict(model, reordered))
    partial = q[:, q.var_names[1:]].copy()
    with pytest.raises(KeyError):
        predict(model, partial)
    genes = shared_genes(ref, partial)
    assert genes == list(q.var_names[1:])
    m = fit_transfer(ref, genes=genes, classifier="logreg", n_pcs=10, cv_folds=3)
    assert len(predict(m, partial)) == 150


def test_thin_counts_and_match_depth():
    rng = np.random.default_rng(0)
    C = sparse.csr_matrix(rng.poisson(5.0, size=(300, 50)))
    T = thin_counts(C, 0.5, seed=1)
    assert T.dtype == np.int32 and (T.data > 0).all()
    assert (T.toarray() <= C.toarray()).all()
    assert abs(T.sum() / C.sum() - 0.5) < 0.02
    assert (thin_counts(C, 0.5, seed=1) != T).nnz == 0
    with pytest.raises(ValueError):
        thin_counts(C, 0.0)
    med = float(np.median(np.asarray(C.sum(1)).ravel()))
    thinned, frac = match_depth(C, target=med / 2)
    assert frac == pytest.approx(0.5)
    assert np.median(np.asarray(thinned.sum(1)).ravel()) == pytest.approx(med / 2, rel=0.05)
    same, frac = match_depth(C, target=med * 1.02)  # within tolerance: untouched
    assert frac == 1.0 and (same != C).nnz == 0


def test_depth_matching_thins_only_deeper_queries(model):
    deep = make_adata(300, [0.25, 0.25, 0.25, 0.25], seed=5, depth=240.0)
    df = predict(model, deep)
    q_med = median_depth(deep)
    assert q_med == float(np.median(np.asarray(deep.layers["counts"].sum(1)).ravel()))
    assert df.attrs["depth_frac"] == pytest.approx(model.depth / q_med)
    assert df.attrs["depth_draws"] == 5  # default: average over five thinning draws
    true_group = np.array([model.group_of[t] for t in deep.obs["celltype"]])
    assert (df["group_raw"].to_numpy() == true_group).mean() > 0.97
    assert predict(model, deep, depth_match=False).attrs["depth_frac"] == 1.0
    single = predict(model, deep, depth_draws=1)
    assert single.attrs["depth_draws"] == 1 and single.attrs["depth_frac"] == df.attrs["depth_frac"]
    np.testing.assert_allclose(single[[f"p:{c}" for c in model.classes]].sum(1), 1.0, atol=1e-6)
    assert (single["group_raw"].to_numpy() == true_group).mean() > 0.97
    assert not np.allclose(single["label_conf"], df["label_conf"])  # different draws, averaged
    shallow = make_adata(300, [0.25, 0.25, 0.25, 0.25], seed=6, depth=40.0)
    assert predict(model, shallow).attrs["depth_frac"] == 1.0
    assert predict(model, shallow, depth_draws=4).attrs["depth_draws"] == 1  # nothing to thin


def test_reference_can_be_thinned_to_a_shallower_query(ref):
    m = fit_transfer(ref, classifier="logreg", n_pcs=10, cv_folds=3, depth=40.0)
    assert m.depth == pytest.approx(40.0, rel=0.1)


def test_novelty_flags_a_cell_type_absent_from_the_reference(model):
    q = make_adata(500, [0.2, 0.2, 0.2, 0.2, 0.2], seed=7)
    df = predict(model, q, novelty_threshold=1.0)
    novel = q.obs["celltype"].to_numpy() == "Neural Tube"
    nov = df["novelty"].to_numpy()
    assert np.median(nov[novel]) > 1.0 > np.median(nov[~novel])
    assert np.median(nov[novel]) > 1.4 * np.median(nov[~novel])
    assert (nov[novel] > 1.0).mean() > 0.9 and (nov[~novel] > 1.0).mean() < 0.05
    assert (df["label"].to_numpy()[novel] == UNASSIGNED).mean() > 0.9
    np.testing.assert_allclose(nov, np.fmax(df["novelty_knn"], df["novelty_residual"]))


def test_excluded_and_rare_labels_are_dropped(ref):
    r = ref.copy()
    lab = r.obs["celltype"].astype(str).to_numpy().copy()
    lab[:3] = "Rare"
    lab[3:40] = "Unknown"
    r.obs["celltype"] = lab
    m = fit_transfer(r, classifier="knn", n_pcs=10, cv_folds=3)
    assert m.dropped["Unknown"] == "excluded" and "Rare" in m.dropped
    assert "Rare" not in m.classes and "Unknown" not in m.classes
    assert m.params["n_train"] == r.n_obs - 40


def test_cross_validate_out_of_fold_and_grouped(ref):
    oof = cross_validate(ref, method="logreg", n_folds=3, n_pcs=10)
    assert len(oof) == ref.n_obs and set(oof["fold"]) == {0, 1, 2}
    assert (oof["true"].to_numpy() == ref.obs["celltype"].to_numpy()).all()
    m = transfer_metrics(oof)
    assert m["accuracy"] > 0.8 and m["group_accuracy"] > 0.97 and 0.0 < m["coverage"] <= 1.0
    assert np.isfinite(m["log_loss"]) and np.isfinite(oof["novelty"]).all()
    sections = section_ids(ref)
    grouped = cross_validate(ref, method="knn", n_folds=3, n_pcs=10, groups=sections, novelty=False)
    per_section = pd.Series(grouped["fold"].to_numpy()).groupby(sections).nunique()
    assert (per_section == 1).all()
    assert grouped.attrs["grouped"] and grouped["novelty"].isna().all()




def test_em_priors_correct_classify_and_count_under_prior_shift():
    rng = np.random.default_rng(0)
    n = 20000
    y = rng.random(n) < 0.8  # query prior 0.8 for class 1; training priors 0.5 / 0.5
    x = np.where(y, 1.0, -1.0) + rng.normal(size=n)
    p1 = 1.0 / (1.0 + np.exp(-2.0 * x))  # exact posterior for N(+-1, 1) under equal priors
    q, W = em_priors(np.column_stack([1 - p1, p1]), [0.5, 0.5])
    assert abs(q[1] - 0.8) < 0.02
    assert abs((p1 > 0.5).mean() - 0.8) > 0.05  # classify-and-count is pulled toward the training prior
    np.testing.assert_allclose(W.sum(1), 1.0)
    # integer weights behave like duplicated cells; a warm start converges to the same priors
    P = np.column_stack([1 - p1, p1])[:500]
    mult = rng.integers(0, 4, size=500)
    dup = np.repeat(P, mult, axis=0)
    q_dup, _ = em_priors(dup, [0.5, 0.5], tol=1e-12)
    q_w, _ = em_priors(P, [0.5, 0.5], weights=mult, tol=1e-12)
    q_warm, _ = em_priors(P, [0.5, 0.5], weights=mult, init=q_w + [0.05, -0.05], tol=1e-12)
    np.testing.assert_allclose(q_w, q_dup, atol=1e-9)
    np.testing.assert_allclose(q_warm, q_w, atol=1e-7)


def test_composition_estimators_sum_to_one(model):
    df = predict(model, make_adata(300, [0.4, 0.1, 0.3, 0.2], seed=9))
    for method in ("hard", "raw", "soft", "em"):
        for level in ("group", "label"):
            assert composition(df, level=level, method=method, model=model).sum() == pytest.approx(1.0)
    with pytest.raises(ValueError):
        composition(df, method="em")


def test_cluster_bootstrap_composition_and_difference():
    clusters = np.repeat(np.arange(10), 10)
    W = one_hot(np.array(["A", "B"] * 50), ["A", "B"])  # every cluster is exactly half A
    bc = bootstrap_composition(W, clusters, n_boot=200)
    np.testing.assert_allclose(bc["estimate"], [0.5, 0.5])
    np.testing.assert_allclose(bc["lo"], bc["hi"])
    rng = np.random.default_rng(1)
    lab = np.where(rng.random(100) < np.repeat(rng.random(10), 10), "A", "B")  # cluster-specific composition
    W2 = one_hot(lab, ["A", "B", "C"])
    b2 = bootstrap_composition(W2, clusters, n_boot=500)
    assert (b2["lo"] <= b2["estimate"]).all() and (b2["estimate"] <= b2["hi"]).all()
    assert (b2.loc[["A", "B"], "hi"] > b2.loc[["A", "B"], "lo"]).all()
    d = bootstrap_difference(W, clusters, W2, clusters, n_boot=300)
    assert list(d.index) == ["A", "B", "C"]
    np.testing.assert_allclose(d["diff"], d["b"] - d["a"])
    assert d.loc["C", "diff"] == 0.0 and not d.loc["C", "ci_excludes_0"]


def _stack(z):
    a = ad.AnnData(X=np.zeros((len(z), 1), np.float32))
    a.obsm["spatial_3D"] = np.column_stack([np.zeros(len(z)), np.zeros(len(z)), np.asarray(z, dtype=float)])
    return a


def test_section_ids_snap_float_noise_without_integer_rounding():
    a = _stack([-15.000000001, -14.99999999, -1e-12, 30.2])
    assert list(section_ids(a)) == ["z-15", "z-15", "z0", "z30.2"]
    # noise straddling an x.5 value stays one section (integer rounding would split 0.5 - e and 0.5 + e)
    b = _stack([0.5 - 1e-12, 0.5 + 1e-12, 10.5, 10.5 + 1e-13, 20.5])
    assert list(section_ids(b)) == ["z0.5", "z0.5", "z10.5", "z10.5", "z20.5"]
    # levels off the integer grid (E8.75: -240.473939 + 14k) keep their values
    rng = np.random.default_rng(0)
    z = -240.473939 + 14.0 * np.repeat(np.arange(5), 3) + rng.normal(0, 1e-11, 15)
    ids = section_ids(_stack(z))
    assert len(set(ids)) == 5 and ids[0] == "z-240.473939" and ids[-1] == "z-184.473939"


def test_z_levels_positions_and_trims_respect_gaps():
    a = _stack([0, 0, 10, 20, 20, 50, 60])  # step 10, levels 30 and 40 missing
    zl = z_levels(a)
    assert list(zl.levels) == [0, 10, 20, 50, 60] and zl.step == pytest.approx(10.0)
    assert zl.gaps()[0]["missing"] == 2
    np.testing.assert_allclose(z_position(zl), [0, 0, 1 / 6, 1 / 3, 1 / 3, 5 / 6, 1])
    assert list(trim_mask(zl, 1, 0)) == [False, False, True, True, True, True, True]
    assert list(trim_mask(zl, 0, 2)) == [True, True, True, True, True, False, False]
    trims = default_trims("WT", "KO")
    assert len(trims) == 8 and trims["none"] == (0, 0, 0, 0) and trims["KO high -6"] == (0, 0, 0, 6)
    assert trims["both ends -4"] == (4, 4, 4, 4) and trims["WT low -6"] == (6, 0, 0, 0)


def test_merged_level_sums_group_probabilities(model):
    assert list(model.merged) == ["CM", "NCC", UNTRACKED]  # Pericardium has no merged category
    assert model.merged_of["Pericardium"] == UNTRACKED and model.merged_of["CM"] == "CM"
    df = predict(model, make_adata(300, [0.25, 0.25, 0.25, 0.25], seed=11))
    np.testing.assert_allclose(df["pm:CM"], df["pg:CM"], atol=1e-12)
    np.testing.assert_allclose(df[f"pm:{UNTRACKED}"], df["pg:Pericardium"], atol=1e-12)
    PM = df[[f"pm:{m}" for m in model.merged]].to_numpy()
    assert (df["merged_raw"].to_numpy() == np.asarray(model.merged)[PM.argmax(1)]).all()
    assert (df["merged_conf"] >= df["group_conf"] - 1e-12).all()
    for method in ("hard", "raw", "soft", "em"):
        assert composition(df, level="merged", method=method, model=model).sum() == pytest.approx(1.0)
    np.testing.assert_allclose(model.projection("merged").sum(1), 1.0)


def test_cross_validate_reports_merged_truth_and_metrics(ref):
    oof = cross_validate(ref, method="logreg", n_folds=3, n_pcs=10, novelty=False)
    assert set(oof["true_merged"]) == {"CM", "NCC", UNTRACKED}
    m = transfer_metrics(oof)
    assert m["merged_accuracy"] >= m["group_accuracy"] - 1e-12 and 0 < m["merged_macro_f1"] <= 1


def test_novelty_excludes_self_matches_on_request(model, ref):
    with pytest.warns(UserWarning, match="exclude_self"):
        biased = predict(model, ref)
    fair = predict(model, ref, exclude_self=True)
    fresh = predict(model, make_adata(800, [0.35, 0.25, 0.2, 0.2], seed=12))
    assert fair.attrs["exclude_self"] and not biased.attrs["exclude_self"]
    assert (fair["novelty_knn"].to_numpy() >= biased["novelty_knn"].to_numpy() - 1e-12).all()
    assert np.median(fair["novelty_knn"]) > 1.03 * np.median(biased["novelty_knn"])
    ratio = np.median(fair["novelty_knn"]) / np.median(fresh["novelty_knn"])
    assert 0.85 < ratio < 1.15  # leave-self-out novelty of reference cells matches that of new cells
    np.testing.assert_allclose(fair["novelty_residual"], biased["novelty_residual"])
    np.testing.assert_allclose(fair["label_conf"], biased["label_conf"])


def test_groups_follow_the_reference_obs_annotation(ref):
    r = ref.copy()
    lab = r.obs["celltype"].astype(str).to_numpy().copy()
    lab[lab == "Peri"] = "LPM"  # a file whose "LPM" is pericardium (as a lineage-map file override would say)
    r.obs["celltype"] = lab
    m = fit_transfer(r, classifier="logreg", n_pcs=10, cv_folds=3)
    assert m.group_of["LPM"] == "LPM"  # lineage map without overrides
    r.obs["group"] = np.where(lab == "LPM", "Pericardium", [m.group_of.get(x, "Unknown") for x in lab])
    m2 = fit_transfer(r, classifier="logreg", n_pcs=10, cv_folds=3)
    assert m2.group_of["LPM"] == "Pericardium" and "Pericardium" in m2.groups
    oof = cross_validate(r, method="logreg", n_folds=3, n_pcs=10, novelty=False)
    assert set(oof.loc[oof["true"] == "LPM", "true_group"]) == {"Pericardium"}


def test_batched_em_matches_the_single_population_loop():
    rng = np.random.default_rng(3)
    P = rng.dirichlet(np.full(4, 0.7), size=400)
    w0 = np.array([0.4, 0.3, 0.2, 0.1])
    W = rng.integers(0, 4, size=(400, 7)).astype(float)
    many = em_priors_many(P, w0, W, tol=1e-12, block=3)
    loop = np.array([em_priors(P, w0, weights=W[:, b], tol=1e-12)[0] for b in range(W.shape[1])])
    np.testing.assert_allclose(many, loop, atol=1e-8)
    clusters = rng.integers(0, 12, 400)
    q0, D = bootstrap_em(P, w0, clusters, n_boot=40, seed=5, tol=1e-12, block=16)
    np.testing.assert_allclose(q0, em_priors(P, w0)[0], atol=1e-7)
    ids, inv = np.unique(clusters, return_inverse=True)
    mult = np.random.default_rng(5).multinomial(len(ids), np.full(len(ids), 1 / len(ids)), size=40)
    for b in (0, 17, 39):
        np.testing.assert_allclose(D[b], em_priors(P, w0, weights=mult[b][inv], tol=1e-12)[0], atol=1e-7)
    q0t, Dt = bootstrap_em(P, w0, clusters, n_boot=40, seed=5, tol=1e-12, block=16, n_jobs=3)  # threads: same draws
    np.testing.assert_array_equal(Dt, D)
    np.testing.assert_array_equal(q0t, q0)


def test_bootstrap_em_pairs_with_top_call_draws():
    # with certain (one-hot) posteriors EM returns the top-call fractions, draw by draw on the same resamples
    rng = np.random.default_rng(4)
    y = rng.integers(0, 3, 300)
    P = np.eye(3)[y]
    clusters = rng.integers(0, 10, 300)
    q0, D = bootstrap_em(P, [0.2, 0.3, 0.5], clusters, n_boot=50, seed=9)
    np.testing.assert_allclose(q0, np.bincount(y, minlength=3) / 300, atol=1e-6)
    np.testing.assert_allclose(D, bootstrap_draws(P, clusters, n_boot=50, seed=9), atol=1e-6)


def test_composition_bracket_and_within_subset(model):
    q = make_adata(600, [0.4, 0.1, 0.3, 0.2], seed=13)
    df = predict(model, q)
    sec = section_ids(q)
    b = composition_bracket(df, model, "group", sec, n_boot=200, estimators=("top_call", "em", "mean_prob"))
    assert list(b.index) == list(model.groups)
    for e in ("top_call", "em", "mean_prob"):
        assert b[e].sum() == pytest.approx(1.0)
        assert (b[f"{e}_lo"] <= b[f"{e}_hi"]).all()
    assert (b["bracket_low"] <= b[["top_call", "em"]].min(1) + 1e-12).all()
    assert (b["bracket_high"] >= b[["top_call", "em"]].max(1) - 1e-12).all()
    assert (b["envelope_lo"] <= b["bracket_low"] + 1e-9).all() or (b["envelope_lo"] <= b["top_call"]).all()
    pd.testing.assert_series_equal(b["em"], composition(df, "group", "em", model).rename("em"), check_names=False)
    cm = ["A-CM", "V-CM"]
    w = composition_bracket(df, model, "label", sec, subset=cm, n_boot=100)
    assert list(w.index) == cm and w["top_call"].sum() == pytest.approx(1.0) and w["em"].sum() == pytest.approx(1.0)
    inside = df["label_raw"].isin(cm)
    assert w.loc["V-CM", "top_call"] == pytest.approx((df.loc[inside, "label_raw"] == "V-CM").mean())
    q_em = composition(df, "label", "em", model)
    assert w.loc["V-CM", "em"] == pytest.approx(q_em["V-CM"] / q_em[cm].sum())
    with pytest.raises(ValueError):
        composition_bracket(df, model, "group", subset=cm)
    point = composition_bracket(df, model, "merged")
    assert "top_call_lo" not in point and point["top_call"].sum() == pytest.approx(1.0)


def test_bracket_difference_detects_a_composition_change(model):
    qa = make_adata(700, [0.15, 0.25, 0.35, 0.25], seed=14)
    qb = make_adata(700, [0.45, 0.25, 0.05, 0.25], seed=15)  # far more V-CM, far less Peri
    pa, pb = predict(model, qa), predict(model, qb)
    d = bracket_difference(pa, pb, model, "group", section_ids(qa), section_ids(qb), n_boot=300)
    assert d.loc["CM", "top_call_diff"] > 0.15 and d.loc["CM", "em_diff"] > 0.15
    assert d.loc["Pericardium", "bracket_high"] < 0 and d.loc["Pericardium", "bracket_excludes_0"]
    np.testing.assert_allclose(d["bracket_low"], d[["top_call_diff", "em_diff"]].min(1))
    np.testing.assert_allclose(d["top_call_diff"], d["b_top_call"] - d["a_top_call"])
    assert d["top_call_diff"].sum() == pytest.approx(0.0, abs=1e-9)
    nob = bracket_difference(pa, pb, model, "label")
    assert "top_call_diff_lo" not in nob and "bracket_excludes_0" not in nob


def test_trim_sensitivity_and_profiles(model):
    qa = make_adata(600, [0.3, 0.2, 0.3, 0.2], seed=16, n_sections=16)
    qb = make_adata(600, [0.4, 0.2, 0.2, 0.2], seed=17, n_sections=18)
    pa, pb = predict(model, qa), predict(model, qb)
    za, zb = z_levels(qa), z_levels(qb)
    specs = {"group": ("group", None), "within_CM": ("label", ["A-CM", "V-CM"])}
    t = trim_sensitivity(pa, za, pb, zb, model, specs=specs, n_boot=100)
    assert t["trim"].nunique() == 8
    assert len(t) == 8 * (len(model.groups) + 2) * 2
    none = t[(t["trim"] == "none") & (t["spec"] == "group") & (t["estimator"] == "top_call")].set_index("category")
    direct = bracket_difference(pa, pb, model, "group", za.index, zb.index, n_boot=100, seeds=(0, 1))
    np.testing.assert_allclose(none["diff"], direct.loc[none.index, "top_call_diff"])
    assert (t.loc[t["trim"] == "both ends -2", "levels_a"] == 12).all() and not t.attrs["skipped_trims"]
    s = trim_summary(t)
    assert len(s) == (len(model.groups) + 2) * 2
    assert (s["diff_min"] <= s["diff_untrimmed"] + 1e-12).all() and (s["diff_max"] >= s["diff_untrimmed"] - 1e-12).all()
    prof = z_profile(pa, za, "group", n_bins=5, categories=model.groups)
    assert prof["n_cells"].sum() == qa.n_obs and len(prof) == 5
    np.testing.assert_allclose(prof[list(model.groups)].sum(1), 1.0)
    loc = z_location(pa, za, "group")
    assert loc["n"].sum() == qa.n_obs and loc["mean_pos"].between(0, 1).all()
    sub = subset_levels(za, np.arange(10))
    assert isinstance(sub, ZLevels) and sub.index.size == 10 and sub.levels.size == za.levels.size == 16


def test_trim_sensitivity_skips_trims_that_empty_a_stack(model):
    qa = make_adata(300, [0.3, 0.2, 0.3, 0.2], seed=18, n_sections=10)
    qb = make_adata(300, [0.3, 0.2, 0.3, 0.2], seed=19, n_sections=16)
    pa, pb = predict(model, qa), predict(model, qb)
    with pytest.warns(UserWarning, match="skipped"):
        t = trim_sensitivity(pa, z_levels(qa), pb, z_levels(qb), model, n_boot=50, estimators=("top_call",))
    assert t.attrs["skipped_trims"] == ["both ends -6"] and t["trim"].nunique() == 7
