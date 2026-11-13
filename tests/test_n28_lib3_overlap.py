"""experiments/t1_strong_null/n28_lib3_overlap.py (row N28 of experiments/harness/TRIES.md, a target-free diagnostic of
T1), on synthetic pseudobulks and a synthetic E10.5-like file only: its redundancy, its eligible genes, its gate, its
refusal of a written output folder, and that its result holds no cell number."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import t1_scope

REPO = Path(__file__).resolve().parents[1]
N_GENES = 4000


@pytest.fixture(scope="module")
def n28():
    spec = importlib.util.spec_from_file_location("t1_strong_null_n28",
                                                  REPO / "experiments" / "t1_strong_null" / "n28_lib3_overlap.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _genes(n: int = N_GENES) -> list[str]:
    return [f"g{i}" for i in range(n)]


def _pb_ref(n: int = N_GENES) -> np.ndarray:
    return np.linspace(5.0, 0.3, n)  # decreasing: gene i has level rank i, every gene above DP_MIN_PB


def test_redundancy_is_one_when_the_base_puts_exactly_the_top_level_genes_at_the_bottom(n28):
    pb_ref = _pb_ref()
    dp = np.zeros(N_GENES)
    dp[: t1_scope.DP_BUDGET] = -1.0 - np.linspace(0, 0.1, t1_scope.DP_BUDGET)
    r = n28.overlap(pb_ref, pb_ref + dp, np.zeros(N_GENES, bool), _genes())
    assert r["top_genes"] == t1_scope.DP_BUDGET
    assert r["redundancy"] == pytest.approx(1.0)
    assert n28.verdict(r["redundancy"]) == "kill"


def test_redundancy_is_zero_when_the_base_puts_none_of_them_at_the_bottom(n28):
    pb_ref = _pb_ref()
    dp = np.zeros(N_GENES)
    dp[: t1_scope.DP_BUDGET] = 0.5
    dp[1000: 1000 + t1_scope.DP_BOTTOM] = -1.0  # dp's bottom tail: genes of low reference level
    r = n28.overlap(pb_ref, pb_ref + dp, np.zeros(N_GENES, bool), _genes())
    assert r["redundancy"] == pytest.approx(0.0)
    assert r["bottom_share_high_level"] == pytest.approx(0.0)
    assert n28.verdict(r["redundancy"]) == "open"


def test_held_genes_and_genes_under_min_pb_are_never_in_top(n28):
    pb_ref = _pb_ref()
    pb_ref[150:] = np.linspace(0.24, 0.0, N_GENES - 150)  # ranks 150 to 799 are high-level bins but under DP_MIN_PB
    held = np.zeros(N_GENES, bool)
    held[:50] = True
    dp = -np.linspace(1.0, 0.0, N_GENES)  # every top-level gene at the bottom: they would all count if eligible
    top = n28.top_genes(pb_ref, dp, held)
    assert top.size == 100
    assert not held[top].any()
    assert (pb_ref[top] >= t1_scope.DP_MIN_PB).all()
    r = n28.overlap(pb_ref, pb_ref + dp, held, _genes())
    assert r["top_genes"] == 100 and r["redundancy"] == pytest.approx(1.0)


def test_top_is_capped_at_the_budget_and_ordered_by_reference_level(n28):
    pb_ref = _pb_ref()
    held = np.zeros(N_GENES, bool)
    held[[3, 7]] = True
    top = n28.top_genes(pb_ref, np.zeros(N_GENES), held)
    assert top.size == t1_scope.DP_BUDGET
    assert list(top) == [i for i in range(t1_scope.DP_BUDGET + 2) if i not in (3, 7)]


@pytest.mark.parametrize("k, expected", [(140, "kill"), (139, "open")])
def test_verdict_follows_the_gate_on_both_sides_through_overlap(n28, k, expected):
    pb_ref = _pb_ref()
    dp = np.zeros(N_GENES)
    dp[:k] = -1.0  # k of the top genes in the bottom tail
    dp[k: t1_scope.DP_BUDGET] = 0.5
    dp[1000: 1000 + t1_scope.DP_BOTTOM] = -0.5  # the rest of the tail filled with low-level genes
    r = n28.overlap(pb_ref, pb_ref + dp, np.zeros(N_GENES, bool), _genes())
    assert r["redundancy"] == pytest.approx(k / t1_scope.DP_BUDGET)
    assert n28.verdict(r["redundancy"]) == expected


@pytest.mark.parametrize("x, expected", [(0.70, "kill"), (0.6999, "open"), (1.0, "kill"), (0.0, "open")])
def test_verdict_gate_constant(n28, x, expected):
    assert n28.REDUNDANCY_KILL == 0.70
    assert n28.verdict(x) == expected


def test_refuses_a_non_empty_output_folder(n28, tmp_path):
    (tmp_path / "n28.json").write_text("{}")
    with pytest.raises(SystemExit) as e:
        n28.main(["--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])
    assert "already written" in str(e.value)


# a synthetic E10.5-like file: three libraries, chamber and surrounding labels, log1p CP10k in X

LIB_CELLS = {"1": 173, "2": 191, "3": 137}  # stand for cell numbers that must never be written


def _labels() -> dict[str, str]:
    """One released label per chamber category and one surrounding label, from the chamber config."""
    spec = t1_scope.load_spec()
    mapping = spec.mapping("E10.5_RNA")
    out = {}
    for cat in [*sorted(spec.in_scope("chamber")), "surrounding"]:
        out[cat] = sorted(lab for lab, c in mapping.items() if c == cat)[0]
    return out


def _ref(n_genes: int = 1200, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    labels = _labels()
    cats = list(labels)
    names, labs = [], []
    for lib, n in LIB_CELLS.items():
        for i in range(n):
            names.append(f"c{lib}x{i}_{lib}")
            labs.append(labels[cats[i % len(cats)]])
    level = np.exp(rng.normal(0, 1.5, n_genes))
    counts = rng.poisson(level[None, :] * rng.uniform(0.5, 1.5, (len(names), 1)))
    lib_size = np.maximum(counts.sum(1, keepdims=True), 1)
    X = np.log1p(1e4 * counts / lib_size).astype(np.float32)
    a = ad.AnnData(X=sparse.csr_matrix(X), obs=pd.DataFrame({"celltype": labs}, index=names),
                   var=pd.DataFrame(index=_genes(n_genes)))
    return a


def test_run_on_a_synthetic_file_writes_summaries_without_cell_numbers(n28):
    ref = _ref()
    genes = list(ref.var_names)
    res = n28.run(ref, genes, held_genes=genes[:5], seed=n28.MID_SEED)
    assert res["verdict"] in n28.VERDICTS
    assert res["verdict"] == n28.verdict(res["base_a"]["redundancy"])
    assert "surrounding" not in res["base_a"]["categories"]
    for key in ("redundancy", "prior_order_genes", "prior_order_median_dp", "bottom_share_high_level",
                "spearman_dp_vs_minus_pb_ref"):
        assert key in res["base_a"] and key in res["base_b_reported_only"]
    text = json.dumps(res)
    forbidden = {*LIB_CELLS.values(), sum(LIB_CELLS.values())}
    ints = set()

    def walk(x):
        if isinstance(x, dict):
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, int) and not isinstance(x, bool):
            ints.add(x)

    walk(res)
    assert not (ints & forbidden), text
    assert not any(f'"{g}"' in text for g in genes)  # no gene name
    report = n28.report_text({**res, "inputs": {}, "stamp": {"commit": "x", "tree_clean": True, "utc": "t"}})
    assert f"Verdict: {res['verdict']}" in report
    assert not any(str(n) in report for n in forbidden)


def test_base_a_is_the_mean_over_library_3_chamber_cells(n28):
    ref = _ref()
    genes = list(ref.var_names)
    in3 = n28.library3_mask(ref)
    pool = ref[np.flatnonzero(in3)].copy()
    shares = n28.harness.t1_chamber_shares(n28.harness._shallow(pool), n28.E105_SOURCE)
    pb = n28.weighted_pseudobulk(pool, genes, shares)
    cats = t1_scope.classify(n28.harness._shallow(pool), source=n28.E105_SOURCE, key_added=None)
    chamber = np.isin(np.asarray(cats).astype(str), list(shares))
    expected = np.asarray(pool.X[np.flatnonzero(chamber)].mean(axis=0)).ravel()
    np.testing.assert_allclose(pb, expected, rtol=1e-6, atol=1e-8)


def test_library_mark_refuses_names_without_a_library(n28):
    ref = _ref()
    ref.obs_names = [f"cell{i}" for i in range(ref.n_obs)]
    with pytest.raises(SystemExit):
        n28.library3_mask(ref)
