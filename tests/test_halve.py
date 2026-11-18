"""`experiments/builds/halve.py` (the halving that made the heart file uploaded on 2026-11-18, docs/PLAN.md decision
log 2026-11-18 (j)) on a synthetic AnnData: the index rule, a repeatable file, nothing else changed, the seed
refusals, and a printout with no cell count. Nothing is read from data/ or runs/."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.harness import BUILD_SEED_OFFSET

REPO = Path(__file__).resolve().parents[1]
N = 1237  # odd, so n // 2 is checked


@pytest.fixture(scope="module")
def hv():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    return halve


def _synthetic(path: Path) -> ad.AnnData:
    rng = np.random.default_rng(3)
    X = sparse.random(N, 40, density=0.2, format="csr", random_state=4, dtype=np.float32)
    obs = pd.DataFrame({"celltype": rng.choice(["a", "b", "c"], N)}, index=[f"cell{i}" for i in range(N)])
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(40)]))
    a.obsm["spatial_3D"] = rng.normal(size=(N, 3)).astype(np.float32)
    a.write_h5ad(path, compression="gzip")
    return a


def test_halve_keeps_the_sorted_half_and_changes_nothing_else(hv, tmp_path):
    src = _synthetic(tmp_path / "in.h5ad")
    seed = BUILD_SEED_OFFSET + 1
    assert hv.DEFAULT_SEED == seed
    digest = hv.halve(str(tmp_path / "in.h5ad"), str(tmp_path / "out.h5ad"), seed)
    idx = np.sort(np.random.default_rng(seed).choice(N, size=N // 2, replace=False))
    np.testing.assert_array_equal(hv.half_rows(N, seed), idx)
    got = ad.read_h5ad(tmp_path / "out.h5ad")
    assert got.n_obs == N // 2 and list(got.obs_names) == [f"cell{i}" for i in idx]
    assert list(got.var_names) == list(src.var_names) and list(got.obs.columns) == ["celltype"]
    assert sparse.issparse(got.X) and got.X.dtype == np.float32
    np.testing.assert_array_equal(got.X.toarray(), src.X[idx].toarray())
    np.testing.assert_array_equal(got.obsm["spatial_3D"], src.obsm["spatial_3D"][idx])
    assert list(got.obs["celltype"]) == list(src.obs["celltype"].iloc[idx])
    assert digest == hashlib.sha256((tmp_path / "out.h5ad").read_bytes()).hexdigest()


def test_halve_is_repeatable_byte_for_byte_and_follows_the_seed(hv, tmp_path):
    _synthetic(tmp_path / "in.h5ad")
    a = hv.halve(str(tmp_path / "in.h5ad"), str(tmp_path / "a.h5ad"))
    b = hv.halve(str(tmp_path / "in.h5ad"), str(tmp_path / "b.h5ad"))
    c = hv.halve(str(tmp_path / "in.h5ad"), str(tmp_path / "c.h5ad"), BUILD_SEED_OFFSET + 2)
    assert a == b and (tmp_path / "a.h5ad").read_bytes() == (tmp_path / "b.h5ad").read_bytes()
    assert c != a


def test_halve_gives_the_bytes_of_the_hand_written_operations(hv, tmp_path):
    # the operations that made the 18 Nov file, written out by hand: the tool must add nothing to them
    _synthetic(tmp_path / "in.h5ad")
    seed = BUILD_SEED_OFFSET + 1
    hv.halve(str(tmp_path / "in.h5ad"), str(tmp_path / "tool.h5ad"), seed)
    a = ad.read_h5ad(tmp_path / "in.h5ad")
    a[np.sort(np.random.default_rng(seed).choice(a.n_obs, a.n_obs // 2, replace=False))].copy().write_h5ad(
        tmp_path / "hand.h5ad", compression="gzip")
    assert (tmp_path / "tool.h5ad").read_bytes() == (tmp_path / "hand.h5ad").read_bytes()


def test_halve_refuses_out_equal_to_in(hv, tmp_path):
    _synthetic(tmp_path / "in.h5ad")
    before = (tmp_path / "in.h5ad").read_bytes()
    (tmp_path / "sub").mkdir()
    for out in (tmp_path / "in.h5ad", tmp_path / "sub" / ".." / "in.h5ad"):
        with pytest.raises(SystemExit):
            hv.main([str(tmp_path / "in.h5ad"), str(out)])
    link = tmp_path / "link.h5ad"
    link.symlink_to(tmp_path / "in.h5ad")
    with pytest.raises(SystemExit):
        hv.halve(str(link), str(tmp_path / "in.h5ad"))
    assert (tmp_path / "in.h5ad").read_bytes() == before


@pytest.mark.parametrize("seed", [0, 1, BUILD_SEED_OFFSET - 1, -5])
def test_halve_refuses_seed_zero_and_seeds_below_the_offset(hv, tmp_path, seed):
    _synthetic(tmp_path / "in.h5ad")
    with pytest.raises(SystemExit):
        hv.main([str(tmp_path / "in.h5ad"), str(tmp_path / "o.h5ad"), "--seed", str(seed)])
    with pytest.raises(SystemExit):
        hv.halve(str(tmp_path / "in.h5ad"), str(tmp_path / "o.h5ad"), seed)
    assert not (tmp_path / "o.h5ad").exists()


def test_halve_prints_ok_and_the_sha256_only(hv, tmp_path, capsys):
    _synthetic(tmp_path / "in.h5ad")
    hv.main([str(tmp_path / "in.h5ad"), str(tmp_path / "o.h5ad")])
    lines = capsys.readouterr().out.splitlines()
    # exactly two lines, the word and a 64-character hex digest: no room for a cell count
    assert lines == ["ok", hashlib.sha256((tmp_path / "o.h5ad").read_bytes()).hexdigest()]
    assert len(lines[1]) == 64 and set(lines[1]) <= set("0123456789abcdef")
