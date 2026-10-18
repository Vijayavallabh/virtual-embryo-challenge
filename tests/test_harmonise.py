import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.harmonise import (
    UNASSIGNED,
    bootstrap_composition,
    bootstrap_difference,
    composition,
    cross_validate,
    em_priors,
    fit_transfer,
    match_depth,
    median_depth,
    one_hot,
    predict,
    section_ids,
    shared_genes,
    thin_counts,
    transfer_metrics,
)
from vec.prep import log1p_cp10k

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


def test_section_ids_round_away_float_noise():
    a = ad.AnnData(X=np.zeros((4, 1), np.float32))
    a.obsm["spatial_3D"] = np.array([[0, 0, -15.000000001], [0, 0, -14.99999999], [0, 0, -1e-12], [0, 0, 30.2]])
    assert list(section_ids(a)) == ["z-15", "z-15", "z0", "z30"]


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
