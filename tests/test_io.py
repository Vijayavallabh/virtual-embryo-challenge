import anndata as ad
import numpy as np
from conftest import GENES

from vec.io import check_submission, make_submission, write_submission


def test_valid_submission_passes(tmp_path, boards, own_repo):
    X = np.random.default_rng(0).gamma(1.0, 1.0, size=(50, len(GENES)))
    coords = np.random.default_rng(1).normal(size=(50, 3))
    path = write_submission(make_submission(X, GENES, coords), tmp_path / "pred.h5ad", repo=own_repo)
    rep = check_submission(path, boards["T2"])
    assert rep.ok, rep.errors


def test_wrong_order_and_missing_coords_rejected(tmp_path, boards):
    X = np.ones((50, len(GENES)), np.float32)
    a = ad.AnnData(X=X)
    a.var_names = GENES[::-1]
    a.write_h5ad(tmp_path / "bad.h5ad")
    rep = check_submission(tmp_path / "bad.h5ad", boards["T2"])
    assert any("different order" in e for e in rep.errors)
    assert any("spatial_3D" in e for e in rep.errors)


def test_negative_values_and_cell_count_rejected(tmp_path, boards):
    X = -np.ones((5, len(GENES)), np.float32)
    a = ad.AnnData(X=X)
    a.var_names = GENES
    a.write_h5ad(tmp_path / "neg.h5ad")
    rep = check_submission(tmp_path / "neg.h5ad", boards["T1"])
    assert any("negative" in e for e in rep.errors)
    assert any("outside" in e for e in rep.errors)


def test_duplicate_coordinates_warned(tmp_path, boards, own_repo):
    X = np.random.default_rng(2).gamma(1.0, 1.0, size=(50, len(GENES)))
    coords = np.random.default_rng(3).normal(size=(50, 3))
    ok = write_submission(make_submission(X, GENES, coords), tmp_path / "ok.h5ad", repo=own_repo)
    rep = check_submission(ok, boards["T2"])
    assert not any("exact coordinates" in w for w in rep.warnings)
    coords[[10, 11]] = coords[3]
    dup = write_submission(make_submission(X, GENES, coords), tmp_path / "dup.h5ad", repo=own_repo)
    rep = check_submission(dup, boards["T2"])
    assert rep.ok
    assert any(w.startswith("2 cells repeat the exact coordinates") for w in rep.warnings)
