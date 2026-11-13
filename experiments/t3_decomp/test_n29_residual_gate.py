"""Tests for diagnostic N29 (experiments/t3_decomp/n29_residual_gate.py), on synthetic factors and synthetic stages."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "n29_residual_gate.py"
spec = importlib.util.spec_from_file_location("n29_residual_gate", SOURCE)
n29 = importlib.util.module_from_spec(spec)
sys.modules["n29_residual_gate"] = n29
spec.loader.exec_module(n29)

RNG = np.random.default_rng(0)


def test_f_is_zero_when_every_group_is_the_pool():
    pool = RNG.normal(size=50)
    r = {g: pool.copy() for g in ("A", "B", "C")}
    w = n29.weights({"A": 300, "B": 500, "C": 200}, ["A", "B", "C"])
    assert n29.fraction(r, pool, w) == pytest.approx(0.0, abs=1e-12)
    assert n29.verdict(n29.fraction(r, pool, w)) == "prune"


def test_f_is_one_when_the_pool_is_zero():
    r = {g: RNG.normal(size=50) for g in ("A", "B", "C")}
    w = n29.weights({"A": 300, "B": 500, "C": 200}, ["A", "B", "C"])
    assert n29.fraction(r, np.zeros(50), w) == pytest.approx(1.0)


def test_weights_follow_the_cell_shares():
    w = n29.weights({"A": 300, "B": 500, "C": 200, "D": 999}, ["A", "B", "C"])
    assert w == pytest.approx({"A": 0.3, "B": 0.5, "C": 0.2})
    # one group off the pool: F^2 is that group's share of the weighted |r|^2
    pool = np.ones(10)
    r = {"A": pool.copy(), "B": pool.copy(), "C": pool + np.r_[3.0, np.zeros(9)]}
    num = 0.2 * 9.0
    den = 0.3 * 10 + 0.5 * 10 + 0.2 * (16 + 9)
    assert n29.fraction(r, pool, w) == pytest.approx(np.sqrt(num / den))
    # moving the same cells' share changes F the way the weights say
    w2 = n29.weights({"A": 300, "B": 200, "C": 500}, ["A", "B", "C"])
    assert n29.fraction(r, pool, w2) == pytest.approx(np.sqrt(0.5 * 9 / (0.3 * 10 + 0.2 * 10 + 0.5 * 25)))


def test_verdict_follows_the_gate_on_both_sides():
    assert n29.F_PRUNE == 0.20
    assert n29.verdict(0.1999) == "prune"
    assert n29.verdict(0.20) == "open"
    assert n29.verdict(0.35) == "open"
    assert n29.verdict(None) == "open"
    unit = np.r_[1.0, np.zeros(9)]
    r = {"A": unit, "B": unit}
    w = {"A": 0.5, "B": 0.5}
    f_low, f_high = n29.fraction(r, 0.85 * unit, w), n29.fraction(r, 0.75 * unit, w)
    assert f_low == pytest.approx(0.15) and n29.verdict(f_low) == "prune"
    assert f_high == pytest.approx(0.25) and n29.verdict(f_high) == "open"


def test_lineage_sets_split_by_the_flag():
    flags = {"CM": "yes", "Neural": "no", "Blood": "partial", "EpiblastStreak": "no"}
    sets = n29.lineage_sets(["Blood", "CM", "EpiblastStreak", "Neural"], flags, ["Neural", "Endoderm"])
    assert sets == {"mesp1": ["CM"], "nonmesp1": ["Neural"], "other": ["Blood", "EpiblastStreak"]}


def test_source_names_no_knockout_file():
    src = SOURCE.read_text()
    for name in ("gata4_ko", "mab21l2_ko", "E8.75_gata4", "E9.5_mab21l2", "_ko_1", "_ko_2", "T3/"):
        assert name not in src, name
    assert not re.search(r"_ko\b", src)
    assert not re.findall(r'"[^"\s{}]+\.h5ad"', src)  # no literal file path: every one is built from a T2 stage name
    assert 'f"T2/{stage}.h5ad"' in src
    assert n29.stage_file("E8.75") == "T2/E8.75.h5ad"


def test_refuses_a_non_empty_output_folder(tmp_path):
    (tmp_path / "n29.json").write_text("{}")
    with pytest.raises(SystemExit, match="already written"):
        n29.main(["--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])


def _stage(shift: dict[str, np.ndarray], n: int = 220, genes: int = 30) -> ad.AnnData:
    rows, grp = [], []
    for g, s in shift.items():
        lam = np.exp(s)[None, :] * np.ones((n, genes))
        rows.append(RNG.poisson(5 * lam))
        grp += [g] * n
    C = sparse.csr_matrix(np.vstack(rows).astype(np.float32))
    a = ad.AnnData(X=C.copy(), layers={"counts": C}, obs=pd.DataFrame({"group": grp}))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.var_names = [f"g{i}" for i in range(genes)]
    return a


def test_factors_and_summary_on_synthetic_stages():
    genes = 30
    base = {g: RNG.normal(scale=0.3, size=genes) for g in ("CM", "Neural", "Blood")}
    days = [7.5, 8.0, 8.3, 8.5, 9.5, 10.5]
    others = [_stage(base, genes=genes) for _ in days]
    dev = RNG.normal(scale=0.5, size=genes)  # the reference's own deviation, shared by every group
    ref = _stage({g: v + dev for g, v in base.items()}, genes=genes)
    fx = n29.factors(ref, others, 8.75, days)
    assert fx["groups"] == ["Blood", "CM", "Neural"]
    flags = {"CM": "yes", "Neural": "no", "Blood": "partial"}
    res = n29.summarise(fx, flags, ["Neural"])
    assert res["F"] is not None and 0 <= res["F"] < 0.5
    assert res["verdict"] == ("prune" if res["F"] < n29.F_PRUNE else "open")
    rep = res["reported"]
    assert rep["n_groups"] == 3 and rep["n_genes"] == genes
    assert set(rep["per_group"]) == {"Blood", "CM", "Neural"}
    assert all(v["cosine_with_pool"] > 0.5 for v in rep["per_group"].values())
    assert rep["F_nonmesp1_pool"] is not None
    out = {"verdict": res["verdict"], **{k: v for k, v in res.items() if k != "verdict"}, "inputs": {},
           "stamp": {"commit": "x", "tree_clean": True, "utc": "t"}}
    assert json.loads(json.dumps(out))["verdict"] == res["verdict"]  # plain JSON: no array is written
    text = n29.report_text(out)
    assert f"Verdict: {res['verdict']}" in text and "g0" not in text  # no gene name
    # a group-specific deviation raises F
    ref2 = _stage({"CM": base["CM"] + dev + RNG.normal(scale=0.8, size=genes), "Neural": base["Neural"] + dev,
                   "Blood": base["Blood"] + dev}, genes=genes)
    res2 = n29.summarise(n29.factors(ref2, others, 8.75, days), flags, ["Neural"])
    assert res2["F"] > res["F"]
    assert res2["reported"]["per_group"]["CM"]["residual_ratio"] > res["reported"]["per_group"]["CM"]["residual_ratio"]
