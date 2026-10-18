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


def test_core_excludes_epicardium_and_merged_level_is_defined():
    spec = lineage_map()
    assert spec["groups"]["EpicardiumST"]["core"] is False
    assert {g for g, v in spec["groups"].items() if v["core"]} == {"CM", "Endocardium"}
    assert spec["groups"]["CardiacProgenitor"]["merged"] is None


def test_file_overrides_apply_only_to_their_file():
    def annotate(source):
        a = ad.AnnData(X=np.zeros((3, 2), np.float32))
        a.obs["celltype"] = ["LPM", "DPM", "V-CM"]
        add_lineage(a, source=source)
        return list(a.obs["group"]), list(a.obs["merged"])

    assert annotate("E8.25_late") == (["LPM", "CardiacProgenitor", "CM"], ["Untracked", "Untracked", "CM"])
    assert annotate("E9.5_mab21l2_ko") == (["Pericardium", "PharyngealMesoderm", "CM"],
                                           ["Untracked", "CranialMesoderm", "CM"])


def _merfish_file(path, rng, n=300, g=30):
    """A small MERFISH-like file with two redundant records appended: one copied a level down, one within its level."""
    counts = rng.poisson(4.0, size=(n, g)).astype(np.float32)
    coords = np.c_[rng.uniform(-100, 100, size=(n, 2)), np.repeat([30.0, 15.0, 0.0], n // 3)]
    counts = np.vstack([counts, counts[[5, 7]]])
    coords = np.vstack([coords, coords[5] - [0, 0, 15.0], coords[7] + [0.001, 0.0, 0.0]])
    a = ad.AnnData(X=log1p_cp10k(counts))
    a.layers["counts"] = sparse.csr_matrix(counts)
    a.obsm["spatial_3D"] = coords
    a.obs["celltype"] = "V-CM"
    a.var_names = [f"g{i}" for i in range(g)]
    a.write_h5ad(path)
    return n


def test_load_stage_keeps_duplicates_by_default_and_drops_them_on_request(tmp_path):
    from vec.prep import load_stage

    path = tmp_path / "E0.0.h5ad"
    n = _merfish_file(path, np.random.default_rng(5))
    full = load_stage(path, lineage=False)
    assert full.n_obs == n + 2 and "qc" not in full.uns
    dedup = load_stage(path, genes=["g3", "g1"], lineage=False, drop_duplicates=True)
    assert dedup.n_obs == n and list(dedup.var_names) == ["g3", "g1"]
    qc = dedup.uns["qc"]
    assert qc["drop_duplicates"] and qc["n_input"] == n + 2 and qc["n_dropped"] == 2
    assert qc["dropped_rows"].tolist() == [n, n + 1]
    assert qc["by_gap"] == {"0": 1, "1": 1} and int(qc["by_level_pair"]["n"].sum()) == 2
    dedup.write_h5ad(tmp_path / "dedup.h5ad")  # the qc record must stay writable
    assert int(ad.read_h5ad(tmp_path / "dedup.h5ad").uns["qc"]["n_dropped"]) == 2


def test_load_stage_drop_duplicates_without_coordinates_is_a_noted_no_op(tmp_path):
    from vec.prep import load_stage

    a = ad.AnnData(X=log1p_cp10k(_counts(np.random.default_rng(6), n=50, g=40)))
    a.var_names = [f"g{i}" for i in range(40)]
    a.write_h5ad(tmp_path / "rna.h5ad")
    out = load_stage(tmp_path / "rna.h5ad", lineage=False, drop_duplicates=True)
    assert out.n_obs == 50 and out.uns["qc"]["n_dropped"] == 0 and "note" in out.uns["qc"]


def test_load_stage_drops_moved_copies_and_interpolated_records_and_never_truncates(tmp_path):
    from vec.prep import load_stage

    rng = np.random.default_rng(9)
    n, g = 300, 30
    counts = rng.poisson(4.0, size=(n, g)).astype(np.float32)
    coords = np.c_[rng.uniform(-100, 100, size=(n, 2)), np.repeat([30.0, 15.0, 0.0], n // 3)]
    moved = counts[3]                                     # copy one level down, shifted 6 units in xy
    interp = (counts[110] + counts[220]) / 2              # average of two cells from neighbouring levels
    assert np.any(interp % 1)
    counts = np.vstack([counts, moved, interp])
    coords = np.vstack([coords, coords[3] + [6.0, 0.0, -15.0], [0.0, 0.0, 15.0]])
    a = ad.AnnData(X=log1p_cp10k(counts))
    a.layers["counts"] = sparse.csr_matrix(counts)
    a.obsm["spatial_3D"] = coords
    a.obs["celltype"] = "V-CM"
    a.var_names = [f"g{i}" for i in range(g)]
    a.write_h5ad(tmp_path / "E0.0.h5ad")

    full = load_stage(tmp_path / "E0.0.h5ad", lineage=False)
    assert full.n_obs == n + 2 and full.layers["counts"].dtype == np.float32
    np.testing.assert_allclose(full.layers["counts"][-1].toarray().ravel(), interp)   # not truncated
    dedup = load_stage(tmp_path / "E0.0.h5ad", lineage=False, drop_duplicates=True)
    qc = dedup.uns["qc"]
    assert dedup.n_obs == n and qc["n_copies"] == 1 and qc["n_interpolated"] == 1
    assert qc["dropped_rows"].tolist() == [n, n + 1] and dedup.layers["counts"].dtype == np.int32
