"""N13 (experiments/harness/TRIES.md), on synthetic data only: its gene sets, per-cell shares and stated reading."""
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "t1_hygiene"))
import n13_assay_check as n13  # noqa: E402


def test_gene_sets_find_each_class_by_pattern():
    genes = ["mt-Co1", "mt-Nd1", "Malat1", "Meg3", "Rpl13", "Rps6", "Rplp0", "Actb"]
    s = n13.gene_sets(genes, ["Meg3"])
    assert list(s["mt"]) == [0, 1]
    assert list(s["Malat1"]) == [2]
    assert list(s["N2"]) == [3]
    assert list(s["Rp"]) == [4, 5]  # Rplp0 has no digit after Rp[ls], so it is not counted


def test_cell_shares_are_per_cell_fractions_and_skip_empty_cells():
    X = sparse.csr_matrix(np.array([[1, 3, 0], [0, 0, 0], [2, 2, 4]], dtype=np.float32))
    out = n13.cell_shares(X, {"a": np.array([0]), "none": np.array([], dtype=int)})
    assert np.allclose(out["a"], [0.25, 0.25])
    assert np.allclose(out["none"], [0.0, 0.0])


def test_pooled_shares_divide_summed_mean_cp10k_by_ten_thousand():
    cp = np.array([500.0, 1500.0, 8000.0])
    assert n13.pooled_shares(cp, {"x": np.array([0, 1])})["x"] == 0.2


def test_reading_follows_the_row_at_its_boundaries():
    assert n13.reading(0.009, 0.05, 0.02, 0.01, 13) == "nuclei"          # below a fifth, Malat1 twice Feng's
    assert n13.reading(0.010, 0.05, 0.02, 0.01, 13) == "unclear"         # exactly a fifth is not below it
    assert n13.reading(0.009, 0.05, 0.019, 0.01, 13) == "unclear"        # Malat1 under twice Feng's
    assert n13.reading(0.025, 0.05, 0.05, 0.01, 13) == "whole cells"     # at least half of Feng's mt- share
    assert n13.reading(0.0, 0.05, 0.02, 0.01, 0) == "nuclei"             # no mt- gene: Malat1 alone
    assert n13.reading(0.0, 0.05, 0.01, 0.01, 0) == "whole cells"
    assert n13.reading(0.0, 0.05, 0.015, 0.01, 0) == "unclear"
