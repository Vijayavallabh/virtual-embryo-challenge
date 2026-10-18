import anndata as ad
import numpy as np
from scipy import sparse

from vec.prep import add_lineage, lineage_map, log1p_cp10k, recover_counts


def _counts(rng, n=200, g=300):
    C = rng.poisson(rng.gamma(0.3, 2.0, size=(n, g))).astype(np.int32)
    C[np.arange(n), rng.integers(0, g, n)] = 1  # guarantee a single count in every cell
    return sparse.csr_matrix(C)


def test_counts_round_trip_exactly():
    rng = np.random.default_rng(0)
    C = _counts(rng)
    X = log1p_cp10k(C)
    rec, lib, resid = recover_counts(X)
    assert resid < 1e-3
    assert (rec != C).nnz == 0
    np.testing.assert_allclose(lib, np.asarray(C.sum(1)).ravel(), rtol=1e-4)


def test_log1p_cp10k_convention():
    rng = np.random.default_rng(1)
    X = log1p_cp10k(_counts(rng)).toarray().astype(np.float64)
    np.testing.assert_allclose(np.expm1(X).sum(1), 1e4, rtol=1e-4)


def test_lineage_map_is_consistent():
    spec = lineage_map()
    assert set(spec["labels"].values()) <= set(spec["groups"])
    for g in spec["groups"].values():
        assert g["mesp1"] in {"yes", "partial", "no"}
        assert isinstance(g["core"], bool)


def test_add_lineage_flags_unmapped_labels():
    a = ad.AnnData(X=np.zeros((3, 2), np.float32))
    a.obs["celltype"] = ["V-CM", "Neural Tube", "Mystery"]
    unmapped = add_lineage(a)
    assert unmapped == ["Mystery"]
    assert list(a.obs["group"]) == ["CM", "Neural", "Unknown"]
    assert list(a.obs["core"]) == [True, False, False]
