"""Tests for experiments/t3_decomp/k25_dispersion.py (try K25, M-B: a mean-preserving stretch of each gene's counts
around its label mean by a per-gene variance ratio from a source knockout) and its launcher run_k25.sh. Synthetic
arrays and synthetic files only: no competition file is opened, nothing is scored against a fold target, and no
committed score is read (the reading's guards are tested on hand-made reference rows)."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec.adoption import SEEDS
from vec.boards import Board
from vec.evaluate import SEVERITY_GATED
from vec.prep import load_stage, log1p_cp10k

REPO = Path(__file__).resolve().parents[1]
T3 = REPO / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import k9_residual as k9  # noqa: E402
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402
import k25_dispersion as k25  # noqa: E402

LAUNCHER = T3 / "run_k25.sh"
ROW_FILE = T3 / "k25_row.txt"  # the row text handed to the main session, word for word


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# Constants, pins and scope


def test_the_rows_constants():
    assert (k25.MIN_LABEL_CELLS, k25.TARGET_SUM, k25.SHRINK, k25.DETECT, k25.R_CLIP) == (50, 1e4, 200.0, 0.10,
                                                                                         (0.5, 2.0))
    assert k25.PREGATE_MIN_GENES == 50 and k25.PREGATE_LN == math.log(1.2)
    assert (k25.ROUND_SEED, k25.PERM_SEED, k25.PERM_BINS) == (20261122, 20261113, 10)
    assert (k25.MIN_GAIN, k25.MIN_AHEAD, k25.MAX_EXTRA_GATED) == (0.3, 8, 1)
    assert (k25.GUARD_MEAN, k25.GUARD_SEED) == (0.1, 0.5) and k25.TRIES == 2
    assert k25.MUTANT_STATE == [("mmd_u", 0.6), ("variogram", 0.4)]
    assert k25.ROW_LABEL == "K25" and k25.ROW_SCRIPT == "`experiments/t3_decomp/k25_dispersion.py`"


def test_the_cells_bases_and_sources():
    assert k25.CELLS == ("plain g1", "plain g2", "dm g1", "dm g2", "x2")
    assert k25.K18_BASE == "barred_form" and k25.X2_BASE == "raw" and k25.STREAM == k18.STREAM
    assert k25.SOURCE_OF_CELL == {"plain g1": "mab21l2", "plain g2": "mab21l2", "dm g1": "mab21l2",
                                  "dm g2": "mab21l2", "x2": "gata4"}
    mab = harness.K4_FOLDS[k25.SOURCES["mab21l2"]]
    assert mab["source_files"] == ("T3/E9.5_mab21l2_ko.h5ad",) and mab["wt_file"] == "T2/E9.5.h5ad"
    gata = harness.K4_FOLDS[k25.SOURCES["gata4"]]
    assert gata["source_files"] == ("T3/E8.75_gata4_ko_1.h5ad", "T3/E8.75_gata4_ko_2.h5ad")
    assert gata["wt_file"] == "T2/E8.75.h5ad"
    # a fold's truth never enters its own ratio: x2's truth is the Mab21l2 knockout, g1's and g2's the Gata4 replicates
    assert k25.X2_TASK.target == "T3/E9.5_mab21l2_ko.h5ad" and k9.FOLDS["x2"] == k25.SOURCES["gata4"]
    for n in k25.TASKS:
        assert k25.TASKS[n].target not in (*mab["source_files"], mab["wt_file"])


def test_only_practice_tasks_run():
    assert k25.practice_keys() == k25.PRACTICE_KEYS
    assert k25.PRACTICE_KEYS == frozenset({"T3w/WT->Gata4rep1", "T3w/WT->Gata4rep2", "T3k9/x2"})
    src = (T3 / "k25_dispersion.py").read_text()
    assert "build_test" not in src and "bcatenin" not in src and "index.json" not in src


def test_the_pinned_files():
    assert _sha(T3 / "k9_residual.py") == k25.K9_MODULE_SHA256
    assert _sha(k25.K9_X2) == k25.K9_X2_SHA256 and k25.pinned_ok() == []
    assert k25.K18_SHA256 == k21.K18_SHA256


def test_existing_files_are_not_changed_by_this_try():
    # K18's, K21's and K9's arms are restored after each swap, and the modules keep their own names
    assert k18.ARMS == k21._K18_ORIGINAL["ARMS"] and k18.build_arm is k21._K18_ORIGINAL["build_arm"]
    assert k9.arm_built is k25._K9_ORIGINAL["arm_built"]


# ---------------------------------------------------------------------------------------------------------------------
# The estimator


def test_weighted_median_is_the_lower_weighted_median():
    assert k25.weighted_median(np.array([3.0, 1.0, 2.0]), np.array([1.0, 1.0, 1.0])) == 2.0
    assert k25.weighted_median(np.array([1.0, 2.0]), np.array([1.0, 1.0])) == 1.0
    assert k25.weighted_median(np.array([1.0, 2.0, 10.0]), np.array([1.0, 1.0, 5.0])) == 10.0
    with pytest.raises(ValueError):
        k25.weighted_median(np.array([]), np.array([]))


def test_residual_variance_removes_depth():
    rng = np.random.default_rng(0)
    N = rng.integers(200, 2000, 400)
    p = np.full(5, 0.02)
    C = rng.poisson(np.outer(N, p))
    C = np.c_[C, N - C.sum(1)]
    v = k25.residual_variance(C)
    y = np.log1p(1e4 * C / C.sum(1, keepdims=True))
    X = np.c_[np.ones(len(N)), np.log(C.sum(1))]
    res = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    assert v.shape == (6,) and np.allclose(v, (res ** 2).sum(0) / (len(N) - 2), rtol=1e-10)
    # a gene that follows depth exactly keeps no residual variance
    D = np.c_[N, 2 * N]
    assert np.allclose(k25.residual_variance(D), 0.0, atol=1e-20)
    assert np.all(np.isnan(k25.residual_variance(C[:2])))
    Z = np.vstack([np.zeros((3, 6), dtype=int), C[:5]])  # empty cells are left out
    assert np.allclose(k25.residual_variance(Z), k25.residual_variance(C[:5]))


def _source(rng, n_wt, n_ko, L, G, sd_wt, sd_ko, depth=400):
    """Cells of L labels with gene log-proportions spread by sd (log-normal multiplicative noise per cell)."""
    def make(n, sd):
        lab = np.arange(n) % L
        base = rng.lognormal(0, 0.5, (L, G))
        mu = base[lab] * np.exp(rng.normal(0, sd, (n, G)))
        C = rng.poisson(mu / mu.sum(1, keepdims=True) * depth * rng.lognormal(0, 0.2, (n, 1)))
        return C, lab
    rng_state = rng.bit_generator.state
    Cw, lw = make(n_wt, sd_wt)
    rng.bit_generator.state = rng_state  # the same label profiles on both sides
    Ck, lk = make(n_ko, sd_ko)
    return Cw, lw, Ck, lk


def test_the_ratio_follows_the_knockouts_spread_and_is_shrunk():
    rng = np.random.default_rng(1)
    Cw, lw, Ck, lk = _source(rng, 600, 600, 3, 30, 0.3, 0.8, depth=3000)
    est = k25.estimate_ratio(Cw, lw, Ck, lk, 3)
    assert est["labels_used"] == [0, 1, 2] and np.all(est["n_g"] == 600)
    assert np.median(est["ln_r"]) > math.log(1.2) and np.all(est["r"] <= 2.0)
    # shrinkage n / (n + 200) by hand, from the weighted median of the per-label ln q
    g = 0
    lnq = [math.log(k25.residual_variance(Ck[lk == t])[g] / k25.residual_variance(Cw[lw == t])[g]) for t in range(3)]
    assert est["ln_r"][g] == pytest.approx(k25.weighted_median(np.array(lnq), np.full(3, 200.0)) * 600 / 800)
    same = k25.estimate_ratio(Cw, lw, Cw.copy(), lw.copy(), 3)
    assert np.allclose(same["ln_r"], 0.0)


def test_labels_under_fifty_cells_on_either_side_are_left_out():
    rng = np.random.default_rng(2)
    Cw, lw, Ck, lk = _source(rng, 300, 300, 3, 10, 0.3, 0.6)
    lk = lk.copy()
    lk[np.flatnonzero(lk == 2)[49:]] = -1  # 49 knockout cells left in label 2
    est = k25.estimate_ratio(Cw, lw, Ck, lk, 3)
    assert est["labels_used"] == [0, 1]
    lw2 = lw.copy()
    lw2[np.flatnonzero(lw2 == 0)[:60]] = -1
    assert k25.estimate_ratio(Cw, lw2, Ck, lk, 3)["labels_used"] == [1]
    none = k25.estimate_ratio(Cw[:40], lw[:40], Ck, lk, 3)
    assert none["labels_used"] == [] and np.all(none["r"] == 1.0) and not none["has_label"].any()


def test_undetected_genes_get_no_weight_and_exactly_one():
    rng = np.random.default_rng(3)
    Cw, lw, Ck, lk = _source(rng, 300, 300, 2, 8, 0.3, 0.9)
    Ck = Ck.copy()
    Ck[:, 0] = 0
    Ck[rng.random(Ck.shape[0]) < 0.05, 0] = 1  # detected in about 5% of the knockout's cells
    est = k25.estimate_ratio(Cw, lw, Ck, lk, 2)
    assert est["n_g"][0] == 0 and est["ln_r"][0] == 0.0 and est["r"][0] == 1.0


def test_the_ratio_is_clipped_to_half_and_two():
    rng = np.random.default_rng(4)
    Cw, lw, Ck, lk = _source(rng, 400, 4000, 2, 6, 0.05, 1.5, depth=5000)
    est = k25.estimate_ratio(Cw, lw, Ck, lk, 2)
    assert est["clipped"].any() and np.max(est["r"]) == 2.0
    low = k25.estimate_ratio(Ck, lk, Cw, lw, 2)
    assert np.min(low["r"]) >= 0.5


def test_the_pregate_needs_fifty_genes_beyond_ln_1_2():
    r = np.ones(200)
    r[:50] = 1.21
    est = {"r": r, "n_g": np.full(200, 300.0), "clipped": np.zeros(200, bool), "has_label": np.ones(200, bool),
           "labels_used": [0, 1]}
    f = k25.pregate_figures(est)
    assert f["genes_big"] == 50 and f["passes"] and f["genes"] == 200 and f["share_clipped"] == 0.0
    r[0] = 1.2
    f = k25.pregate_figures(est | {"r": r})
    assert f["genes_big"] == 49 and not f["passes"]
    r[0] = 1 / 1.21
    assert k25.pregate_figures(est | {"r": r})["passes"]
    assert sum(f["n_g_histogram"].values()) == 200


# ---------------------------------------------------------------------------------------------------------------------
# The stretch


def _counts(rng, n=500, G=12, L=4):
    lab = rng.integers(0, L, n)
    C = rng.poisson(rng.gamma(1.5, 3.0, (L, G))[lab] * rng.lognormal(0, 0.6, (n, G))).astype(np.int32)
    return C, lab


def test_r_one_returns_the_input_itself():
    rng = np.random.default_rng(5)
    C, lab = _counts(rng)
    out = k25.stretch(C, lab, np.ones(C.shape[1]))
    assert out is C and out.tobytes() == C.tobytes()


def test_genes_at_one_are_untouched_and_the_output_is_integer_and_non_negative():
    rng = np.random.default_rng(6)
    C, lab = _counts(rng)
    r = np.ones(C.shape[1])
    r[[1, 4]] = [2.0, 0.5]
    out = k25.stretch(C, lab, r)
    assert out.dtype == C.dtype and out.shape == C.shape and out.min() >= 0
    keep = [j for j in range(C.shape[1]) if j not in (1, 4)]
    assert np.array_equal(out[:, keep], C[:, keep])
    assert not np.array_equal(out[:, 1], C[:, 1]) and not np.array_equal(out[:, 4], C[:, 4])
    assert np.array_equal(out, k25.stretch(C, lab, r))  # one fixed seed: reproducible


def test_the_stretch_moves_the_spread_and_keeps_the_label_means():
    rng = np.random.default_rng(7)
    lab = np.repeat([0, 1], 4000)
    C = np.c_[rng.poisson(50, 8000), rng.poisson(50, 8000)].astype(np.int32)
    out = k25.stretch(C, lab, np.array([1.69, 0.64]))
    for t in (0, 1):
        a, b = C[lab == t].astype(float), out[lab == t].astype(float)
        assert abs(b.mean(0) - a.mean(0)).max() < 0.05  # stochastic rounding only
        assert b[:, 0].var() / a[:, 0].var() == pytest.approx(1.69, rel=0.05)
        assert b[:, 1].var() / a[:, 1].var() == pytest.approx(0.64, rel=0.08)


def test_negatives_are_clipped_and_their_mass_taken_from_the_positive_cells():
    lab = np.zeros(4, dtype=int)
    C = np.array([[0], [0], [0], [12]], dtype=np.int32)
    seen = {}
    real = np.random.default_rng

    class NoNoise:
        def __init__(self, seed):
            seen["seed"] = seed

        def random(self, shape):
            return np.ones(shape)  # never rounds up: floor only

    k25.np.random.default_rng = NoNoise
    try:
        out = k25.stretch(C, lab, np.array([4.0]))
    finally:
        k25.np.random.default_rng = real
    # m = 3; m + 2 (c - m) = -3, -3, -3, 21 -> 0, 0, 0, 21 - 9 = 12: the mean stays 3
    assert seen["seed"] == k25.ROUND_SEED and out.ravel().tolist() == [0, 0, 0, 12]


def test_the_mean_check_raises_before_rounding(monkeypatch):
    rng = np.random.default_rng(8)
    C, lab = _counts(rng)
    monkeypatch.setattr(k25, "MEAN_TOL", -1.0)
    with pytest.raises(AssertionError):
        k25.stretch(C, lab, np.full(C.shape[1], 1.5))


def test_the_stretch_refuses_bad_ratios():
    C, lab = _counts(np.random.default_rng(9))
    for r in (np.ones(3), np.full(C.shape[1], np.nan), np.zeros(C.shape[1])):
        with pytest.raises(ValueError):
            k25.stretch(C, lab, r)


def test_the_control_is_the_decile_permutation_at_k21s_seed():
    rng = np.random.default_rng(10)
    r = np.exp(rng.normal(0, 0.3, 100))
    basis = rng.gamma(2.0, size=100)
    p = k25.permuted_ratio(r, basis)
    assert np.array_equal(p, harness.t3_decile_permutation(r, basis, 20261113, 10))
    assert sorted(p) == sorted(r) and not np.array_equal(p, r)
    for genes in np.array_split(np.argsort(basis, kind="stable"), 10):
        assert sorted(p[genes]) == sorted(r[genes])


# ---------------------------------------------------------------------------------------------------------------------
# The capture and the swaps


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


def test_captured_counts_records_and_restores(bs):
    real = bs.t3d.lognorm
    C = np.array([[1, 2], [3, 4]], dtype=np.int32)
    with k25.captured_counts(bs.t3d) as calls:
        X = bs.t3d.lognorm(C)
    assert bs.t3d.lognorm is real and len(calls) == 1 and np.array_equal(calls[0], C)
    assert np.array_equal(X, real(C))


def _built(bs, C, lab):
    return bs.Built(X=bs.t3d.lognorm(C), lab=lab, source=np.full(len(lab), "wt", dtype=object),
                    row=np.arange(len(lab)))


def test_stretched_is_the_base_with_its_counts_stretched(bs):
    rng = np.random.default_rng(11)
    C, lab = _counts(rng)
    base = lambda: _built(bs, C, lab)  # noqa: E731
    r = np.ones(C.shape[1])
    assert k25.stretched(bs, base, r).X.tobytes() == base().X.tobytes()
    r[2] = 1.8
    got = k25.stretched(bs, base, r)
    assert np.array_equal(got.X, bs.t3d.lognorm(k25.stretch(C, lab, r)))
    assert np.array_equal(got.lab, lab) and np.array_equal(got.row, np.arange(len(lab)))
    with pytest.raises(RuntimeError):
        k25.stretched(bs, lambda: (bs.t3d.lognorm(C), _built(bs, C, lab))[1], r)  # two lognorm calls


def test_k18s_swap_builds_the_base_unchanged_and_the_arms_from_its_counts(bs, monkeypatch):
    rng = np.random.default_rng(12)
    C, lab = _counts(rng)
    calls = []

    def fake_original(bs_, ref, arm, n, seed, F):
        calls.append(arm)
        return _built(bs_, C, lab)

    monkeypatch.setitem(k25._K18_ORIGINAL, "build_arm", fake_original)
    F = {"r": np.full(C.shape[1], 1.6), "r_perm": np.full(C.shape[1], 0.7)}
    with k25.k18_arms():
        assert k18.ARMS == k25.K18_ARMS and k18.build_arm is k25.build_arm
        base = k18.build_arm(bs, None, "barred_form", 0, 3, F)
        disp = k18.build_arm(bs, None, "disp", 0, 3, F)
        ctl = k18.build_arm(bs, None, "disp_perm", 0, 3, F)
        with pytest.raises(RuntimeError):
            k18.build_arm(bs, None, "refused_form", 0, 3, F)
    assert calls == ["barred_form"] * 3
    assert np.array_equal(base.X, bs.t3d.lognorm(C))
    assert np.array_equal(disp.X, bs.t3d.lognorm(k25.stretch(C, lab, F["r"])))
    assert np.array_equal(ctl.X, bs.t3d.lognorm(k25.stretch(C, lab, F["r_perm"])))
    assert k18.ARMS == k21._K18_ORIGINAL["ARMS"] and k18.build_arm is k21._K18_ORIGINAL["build_arm"]
    with pytest.raises(ValueError), k25.k18_arms():
        raise ValueError("restored on an exception too")
    assert k18.fold_factors is k21._K18_ORIGINAL["fold_factors"]


def test_fold_factors_adds_the_mab21l2_ratio_in_the_references_order(monkeypatch):
    monkeypatch.setitem(k25._K18_ORIGINAL, "fold_factors", lambda *a: {"k18": True})
    monkeypatch.setitem(k25._DISP, "mab21l2", {"genes": ["a", "b", "c"], "r": [1.0, 2.0, 0.5]})
    ref = SimpleNamespace(genes=["c", "a", "b"])
    basis = np.array([3.0, 1.0, 2.0])
    F = k25.fold_factors(None, ref, {}, None, None, basis, None, {})
    assert F["k18"] and F["r"].tolist() == [0.5, 1.0, 2.0]
    assert np.array_equal(F["r_perm"], harness.t3_decile_permutation(F["r"], basis, 20261113, 10))


# ---------------------------------------------------------------------------------------------------------------------
# End to end on synthetic files: the source estimator through the frozen harmoniser, and x2's arms through K9's build


E875 = [("V-CM", "OFT/RV-CM", 0.15), ("V-CM", "vCM1", 0.15), ("aSHF", "Unknown", 0.10), ("Endo", "Unknown", 0.08),
        ("NCC", "Unknown", 0.14), ("SE", "Unknown", 0.13), ("Foregut", "Unknown", 0.12),
        ("Neural Tube", "Unknown", 0.13)]
E95 = [("V-CM", None, 0.30), ("A-CM", None, 0.10), ("SE", None, 0.15), ("Foregut", None, 0.15), ("NCC", None, 0.20),
       ("Endo", None, 0.10)]
GENES = [f"G{j}" for j in range(40)]


def _label(cell) -> str:
    return cell[1] if cell[1] not in (None, "Unknown") else cell[0]


def _cells(cells, prof, n, rng, depth, factor=None, name="c", spread=0.0) -> ad.AnnData:
    w = np.array([c[2] for c in cells])
    lab = rng.choice(len(cells), n, p=w / w.sum())
    f = np.ones(len(GENES)) if factor is None else factor
    C = np.asarray([rng.poisson(prof[_label(cells[i])] * f * np.exp(rng.normal(0, spread, len(GENES)))
                                / (prof[_label(cells[i])] * f).sum() * rng.lognormal(np.log(depth), 0.3))
                    for i in lab], dtype=np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [cells[i][0] for i in lab]}, index=[f"{name}_{i}" for i in range(n)])
    if any(c[1] is not None for c in cells):
        obs["cm_celltype"] = [cells[i][1] or "Unknown" for i in lab]
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 500, n), rng.uniform(0, 500, n), rng.integers(0, 20, n) * 10.0]
    return a


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    """Synthetic files at the release's paths (x2's roles and the two sources), written under a temporary data dir."""
    d = tmp_path_factory.mktemp("k25data")
    (d / "T2").mkdir()
    (d / "T3").mkdir()
    rng = np.random.default_rng(0)
    prof = {k: rng.lognormal(0, 1.0, len(GENES)) for k in sorted({_label(c) for c in E875 + E95})}
    files = {"T2/E8.75.h5ad": _cells(E875, prof, 900, rng, 150, name="e875", spread=0.2),
             "T2/E9.5.h5ad": _cells(E95, prof, 900, rng, 120, name="e95", spread=0.2),
             "T3/E8.75_gata4_ko_1.h5ad": _cells(E875, prof, 600, rng, 140, rng.lognormal(0, 0.4, len(GENES)),
                                                name="g1", spread=0.6),
             "T3/E8.75_gata4_ko_2.h5ad": _cells(E875, prof, 600, rng, 140, rng.lognormal(0, 0.4, len(GENES)),
                                                name="g2", spread=0.6),
             "T3/E9.5_mab21l2_ko.h5ad": _cells(E95, prof, 700, rng, 110, rng.lognormal(0, 0.4, len(GENES)),
                                               name="mab", spread=0.6)}
    for s in ("E7.5", "E8.0", "E8.25_late", "E8.5", "E10.5"):
        files[f"T2/{s}.h5ad"] = _cells(E875, prof, 300, rng, 120, rng.lognormal(0, 0.1, len(GENES)), name=s)
    for f, a in files.items():
        a.write_h5ad(d / f)
    return d


def test_the_source_estimator_runs_through_the_harmoniser_and_reads_no_target(data, monkeypatch):
    opened = []
    real = load_stage

    def spy(path, *a, **kw):
        opened.append(str(Path(path).relative_to(data)))
        return real(path, *a, **kw)

    monkeypatch.setattr("vec.prep.load_stage", spy)
    for s, fold in k25.SOURCES.items():
        opened.clear()
        est, rec = k25.source_ratio(s, data)
        f = harness.K4_FOLDS[fold]
        assert sorted(opened) == sorted([f["wt_file"], *f["source_files"]])
        assert len(rec["genes"]) == len(GENES) and len(rec["r"]) == len(GENES)
        assert rec["record"]["files_sha256"] == {p: _sha(data / p) for p in opened}
        assert est["labels_used"] and np.all((np.asarray(rec["r"]) >= 0.5) & (np.asarray(rec["r"]) <= 2.0))
        assert np.median(est["ln_r"]) > 0  # the synthetic knockouts carry more spread
        assert json.loads(json.dumps(rec)) == rec


def test_the_dm_worlds_genes_are_the_gata4_sources_wild_type_genes():
    # needed_genes takes the dm world's reference genes from the Gata4 source's wild type: the same file, loaded the
    # same way (load_stage, duplicates dropped) as t3_decompose.build_context loads orth_diag's "wt"
    import orth_diag as od

    assert harness.K4_FOLDS[k25.SOURCES["gata4"]]["wt_file"] == od.WT_FILES["E8.75"]
    src = (T3 / "t3_decompose.py").read_text()
    assert 'wt_ad = load_stage(paths["wt"], drop_duplicates=True)' in src and "genes = list(wt_ad.var_names)" in src
    mine = (T3 / "k25_dispersion.py").read_text()
    assert 'wt = load_stage(Path(data) / fold["wt_file"], drop_duplicates=True)' in mine


def test_needed_genes_and_uncovered():
    ratios = {"mab21l2": {"genes": ["a", "b", "c", "d"]}, "gata4": {"genes": ["a", "b", "c"]}}
    panel = {"g1": ["a", "b"], "g2": ["b", "c"], "x2": ["a", "c"]}
    need = k25.needed_genes(ratios, panel)
    assert need == {"mab21l2": ["a", "b", "c"], "gata4": ["a", "c"]}
    assert k25.uncovered(ratios, need) == []
    assert k25.uncovered(ratios, k25.needed_genes(ratios, panel | {"x2": ["a", "d"]})) == ["gata4"]
    short = {"mab21l2": {"genes": ["a", "b"]}, "gata4": {"genes": ["a", "b", "c"]}}  # misses the dm world's "c"
    assert k25.uncovered(short, k25.needed_genes(short, {"g1": ["a"], "g2": ["b"], "x2": ["a"]})) == ["mab21l2"]
    p2 = {"mab21l2": {"genes": ["a", "b", "c"]}, "gata4": {"genes": ["a", "b"]}}
    assert k25.uncovered(p2, k25.needed_genes(p2, {"g1": ["a", "z"], "g2": ["b"], "x2": ["a"]})) == ["mab21l2"]


def _pregate_world(monkeypatch, genes: dict, panel: dict, big: int = 60) -> dict:
    st = {"commit": "a" * 40, "tree_clean": True, "machine": "m"}
    monkeypatch.setattr(k25, "preconditions", lambda commit: st)

    def fake_source(s, data):
        G = len(genes[s])
        r = np.ones(G)
        r[: min(big, G)] = 1.5
        est = {"r": r, "ln_r": np.log(r), "n_g": np.full(G, 300.0), "labels_used": [0],
               "has_label": np.ones(G, dtype=bool), "clipped": np.zeros(G, dtype=bool)}
        return est, {"genes": list(genes[s]), "r": [float(x) for x in r], "record": {"files_sha256": {}}}

    monkeypatch.setattr(k25, "source_ratio", fake_source)
    monkeypatch.setattr(k25, "panel_genes", lambda: panel)
    return st


def test_the_pregate_refuses_at_no_try_when_a_reference_gene_is_not_covered(tmp_path, monkeypatch, capsys):
    full = [f"G{j}" for j in range(80)]
    out = tmp_path / "o"
    # the dm world needs all 80 of WT E8.75's genes (the Gata4 source's); the Mab21l2 source holds 70
    _pregate_world(monkeypatch, {"mab21l2": full[:70], "gata4": full}, {"g1": full[:60], "g2": full[:60],
                                                                        "x2": full[:60]})
    with pytest.raises(SystemExit) as e:
        k25.main(["pregate", "--commit", "a" * 40, "--out", str(out)])
    assert e.value.code == k25.EXIT_UNCOVERED == 6
    assert not (out / "pregate.json").exists() and not (out / "pregate.json.prov.json").exists()
    printed = capsys.readouterr().out
    assert "refused" in printed and not re.search(r"\d", printed.replace("K25", ""))
    # x2's panel not covered by the Gata4 source
    _pregate_world(monkeypatch, {"mab21l2": full, "gata4": full[:50]}, {"g1": full[:40], "g2": full[:40],
                                                                        "x2": full[40:60]})
    with pytest.raises(SystemExit) as e:
        k25.main(["pregate", "--commit", "a" * 40, "--out", str(out)])
    assert e.value.code == k25.EXIT_UNCOVERED and not (out / "pregate.json").exists()
    # covered: written, and its pass or fail is the 50-gene rule in each source
    _pregate_world(monkeypatch, {"mab21l2": full, "gata4": full}, {"g1": full, "g2": full, "x2": full})
    k25.main(["pregate", "--commit", "a" * 40, "--out", str(out)])
    body = json.loads((out / "pregate.json").read_text())
    assert body["passes"] is True and body["tries_if_run"] == 2
    _pregate_world(monkeypatch, {"mab21l2": full, "gata4": full}, {"g1": full, "g2": full, "x2": full}, big=49)
    with pytest.raises(SystemExit) as e:
        k25.main(["pregate", "--commit", "a" * 40, "--out", str(tmp_path / "o49")])
    assert e.value.code == k25.EXIT_PREGATE_FAILS
    assert json.loads((tmp_path / "o49" / "pregate.json").read_text())["tries_if_run"] == 0


def test_the_pregate_needs_fifty_genes_in_each_source(tmp_path, monkeypatch):
    full = [f"G{j}" for j in range(80)]
    st = _pregate_world(monkeypatch, {"mab21l2": full, "gata4": full}, {"g1": full, "g2": full, "x2": full})
    fake = k25.source_ratio

    def one_short(s, data):
        est, rec = fake(s, data)
        if s == "gata4":
            est["r"][:] = 1.0
        return est, rec

    monkeypatch.setattr(k25, "source_ratio", one_short)
    with pytest.raises(SystemExit) as e:
        k25.main(["pregate", "--commit", st["commit"], "--out", str(tmp_path / "o")])
    assert e.value.code == k25.EXIT_PREGATE_FAILS
    body = json.loads((tmp_path / "o" / "pregate.json").read_text())
    assert body["figures"]["mab21l2"]["passes"] and not body["figures"]["gata4"]["passes"] and not body["passes"]


def test_panel_genes_reads_the_tasks_boards(monkeypatch):
    class B:
        def __init__(self, k):
            self.genes = [f"{k}:x"]

    class Boards(dict):
        def __missing__(self, k):
            return B(k)

    monkeypatch.setattr(harness, "load_boards", lambda: Boards())
    got = k25.panel_genes()
    assert set(got) == {"g1", "g2", "x2"} and got["x2"] == [f"{k25.X2_TASK.panel_board}:x"]
    assert got["g1"] == [f"{k25.TASKS['g1'].panel_board}:x"] and got["g2"] == [f"{k25.TASKS['g2'].panel_board}:x"]


def _ctx(task, inputs, seed=harness.BUILD_SEED_OFFSET + 3):
    board = Board(key=f"pseudo:{task.key}", task="T3", genes=list(GENES), needs_coords=True, split="pseudo",
                  target_ko="gata4", min_cells=100, max_cells=900, truth_cells=120, ref_cells=150)
    return harness.BuildContext(task=harness.TaskView.of(task), board=harness._candidate_board(board), inputs=inputs,
                                seed=seed)


def test_x2s_arms_are_k9s_raw_build_then_the_stretch(data, monkeypatch):
    import k7_wt_revert as k7

    monkeypatch.setattr(k7, "MIN_CELLS", 20)
    t = k25.X2_TASK
    inputs = {r: load_stage(data / f, genes=GENES, drop_duplicates=True) for r, f in t.inputs.items()}
    k9._PREP.clear()
    k25._X2.clear()
    rng = np.random.default_rng(3)
    monkeypatch.setitem(k25._DISP, "gata4", {"genes": list(GENES), "r": list(np.exp(rng.normal(0, 0.3, 40)))})
    ctx = _ctx(t, inputs)
    with harness.guarded_reads() as g:
        out = {a: k25._build_x2(ctx, a) for a in k25.X2_ARMS}
        raw = k9._build(ctx, k9.BASE)
    assert not g["refused"]
    assert np.array_equal(out["raw"].X, raw.X) and np.array_equal(out["raw"].coords, raw.coords)
    for a in ("disp", "disp_perm"):
        assert np.array_equal(out[a].coords, raw.coords) and not np.array_equal(out[a].X, raw.X)
    assert not np.array_equal(out["disp"].X, out["disp_perm"].X)
    assert k9.arm_built is k25._K9_ORIGINAL["arm_built"]
    pool = inputs[ctx.task.reference_role]
    assert np.array_equal(k25._X2["r_perm"], k25.permuted_ratio(k25._X2["r"], harness.t3_mean_expression(pool)))
    # r = 1: the base byte for byte
    k25._X2.clear()
    monkeypatch.setitem(k25._DISP, "gata4", {"genes": list(GENES), "r": [1.0] * 40})
    with harness.guarded_reads():
        one = k25._build_x2(ctx, "disp")
    assert one.X.tobytes() == raw.X.tobytes() and one.coords.tobytes() == raw.coords.tobytes()
    k9._PREP.clear()
    k25._X2.clear()


# ---------------------------------------------------------------------------------------------------------------------
# The reading


SKILLS = {"de_score": 0.5, "de_direction": 0.5, "severity_slope": 0.5, "mmd_u": 0.5, "variogram": 0.5}
PROV = {"commit": "c" * 40, "tree_clean": True, "machine": "dgx"}


def _row(seed, base, d_arm, d_ctl, *, floor=True, gated=(0, 0, 0), ms=0.01, cell=""):
    b = 52.0 + 0.1 * seed
    sc = {base: b, "disp": b + d_arm[seed], "disp_perm": b + d_ctl[seed]}
    if floor:
        sc[harness.FLOOR] = 50.0
    sk = {a: dict(SKILLS) for a in sc}
    sk["disp"] = SKILLS | {"mmd_u": 0.5 + ms, "variogram": 0.5 + ms}
    st = {}
    for a, k in zip((base, "disp", "disp_perm"), gated, strict=True):
        st[a] = {"severity": {"state": "r2_gate" if seed < k else "scored", "beta": 1.0, "r2": 0.5}}
    if floor:
        st[harness.FLOOR] = {"severity": {"state": "scored"}}
    return {"seed": seed, "scores": sc, "skills": sk, "t3": st}


def _const(v):
    return [v] * len(SEEDS)


def _runs(tmp: Path, spec: dict | None = None, prov: dict | None = None, drop: str | None = None,
          pre: dict | None = None) -> tuple[Path, Path, Path]:
    """A synthetic run folder, K18's and K9's reference rows (matching the bases) and their pins."""
    spec = spec or {}
    runs = tmp / "runs"
    runs.mkdir(parents=True)
    pre_body = {"passes": True, "figures": {}, "ratios": {}, "provenance": PROV} | (pre or {})
    raw = json.dumps(pre_body).encode()
    (runs / "pregate.json").write_bytes(raw)
    psha = hashlib.sha256(raw).hexdigest()
    for c in k25.CELLS:
        d_arm, d_ctl, gated, ms = spec.get(c, (_const(1.0), _const(0.0), (0, 0, 0), 0.01))
        for s in SEEDS:
            row = _row(s, k25.base_of(c), d_arm, d_ctl, floor=not c.startswith("dm"), gated=gated, ms=ms)
            p = (prov or {}).get((c, s), PROV)
            if c.startswith("dm"):
                f = runs / "dm" / f"{c[3:]}_seed_{s}.json"
                body = row | {"task": c[3:], "provenance": p, "pregate_sha256": psha}
            else:
                key = k25.X2_TASK.key if c == "x2" else k25.TASKS[c[6:]].key
                f = runs / ("x2" if c == "x2" else f"plain_{c[6:]}") / f"seed_{s}.json"
                body = {"task": {"key": key}, "rows": [row], "provenance": p, "pregate_sha256": psha}
            if drop == f"{c} {s}":
                continue
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text(json.dumps(body))
    refs = tmp / "k18"
    refs.mkdir()
    rows = lambda: [{"seed": s, "scores": {"barred_form": 52.0 + 0.1 * s}} for s in SEEDS]  # noqa: E731
    (refs / "plain_g1.json").write_text(json.dumps({"rows": rows()}))
    (refs / "plain_g2.json").write_text(json.dumps({"rows": rows()}))
    (refs / "dm.json").write_text(json.dumps({"rows": {"g1": rows(), "g2": rows()}}))
    x2 = tmp / "x2.json"
    x2.write_text(json.dumps({"rows": [{"seed": s, "scores": {"raw": 52.0 + 0.1 * s}} for s in SEEDS]}))
    return runs, refs, x2


@pytest.fixture
def read(tmp_path, monkeypatch):
    monkeypatch.setattr(k18, "dm_check", lambda rows, *a: {})

    def go(name="r", **kw):
        runs, refs, x2 = _runs(tmp_path / name, **kw)
        monkeypatch.setattr(k25, "K18_SHA256", {f: _sha(refs / f) for f in ("plain_g1.json", "plain_g2.json",
                                                                             "dm.json")})
        monkeypatch.setattr(k25, "K9_X2_SHA256", _sha(x2))
        return k25.reading(runs, refs, x2)

    return go


def test_a_pass_needs_every_clause(read):
    res = read()
    assert not res["void"] and res["passes"] and res["tries_counted"] == 2
    assert res["ahead_cells"] == list(k25.CELLS) and res["outcome"] == k25.PASS_TEXT
    assert res["worlds"]["x2"]["base"] == "raw" and res["worlds"]["dm g1"]["base"] == "barred_form"


def test_k21s_pins_are_restored_after_a_read(read):
    before = dict(k21.K18_SHA256)
    read()
    assert k21.K18_SHA256 == before


def _flat(cells, d_arm, d_ctl=None, gated=(0, 0, 0), ms=0.01):
    return {c: (d_arm, d_ctl if d_ctl is not None else _const(0.0), gated, ms) for c in cells}


def test_clause_a_needs_x2_and_a_dm_cell_and_a_plain_cell(read):
    assert not read("x2", spec=_flat(["x2"], _const(0.0)))["a"]
    assert not read("dm", spec=_flat(["dm g1", "dm g2"], _const(0.0)))["a"]
    assert not read("pl", spec=_flat(["plain g1", "plain g2"], _const(0.0)))["a"]
    one_each = read("one", spec=_flat(["dm g2", "plain g1"], _const(0.0)))
    assert one_each["a"] and one_each["passes"] and one_each["ahead_cells"] == ["plain g2", "dm g1", "x2"]


def test_ahead_needs_a_mean_above_0_3_and_eight_of_ten(read):
    seven = [1.0] * 7 + [-0.1] * 3  # mean 0.67, seven seeds ahead
    res = read("seven", spec=_flat(["x2"], seven))
    assert "x2" not in res["ahead_cells"] and not res["a"]
    eight = [0.5] * 8 + [-0.5] * 2  # mean 0.30: not above
    assert "x2" not in read("eight", spec=_flat(["x2"], eight))["ahead_cells"]


def test_clause_b_no_cell_behind(read):
    res = read("behind", spec=_flat(["plain g2"], _const(-0.31)))
    assert not res["b"] and res["behind_cells"] == ["plain g2"] and not res["passes"]
    assert read("level", spec=_flat(["plain g2"], _const(-0.3)))["b"]


def test_clause_c_gated_seeds_at_most_the_bases_plus_one(read):
    assert read("plus1", spec=_flat(["dm g1"], _const(1.0), gated=(1, 2, 0)))["c"]
    res = read("plus2", spec=_flat(["dm g1"], _const(1.0), gated=(1, 3, 0)))
    assert not res["c"] and not res["passes"]
    # a cell not ahead counts too
    assert not read("plus2b", spec=_flat(["dm g2"], _const(0.0), gated=(0, 2, 0)))["c"]


def test_clause_d_beats_the_control_in_every_ahead_cell(read):
    res = read("ctl", spec=_flat(["plain g1"], _const(1.0), _const(0.8)))
    assert not res["d"] and not res["passes"]
    # a cell that is not ahead need not beat the control
    assert read("ctl2", spec=_flat(["plain g1"], _const(0.0), _const(0.8)))["d"]


def test_clause_e_mutant_state_up_in_every_ahead_cell(read):
    res = read("ms", spec=_flat(["x2"], _const(1.0), ms=-0.01))
    assert not res["e"] and not res["passes"]
    assert not read("ms0", spec=_flat(["x2"], _const(1.0), ms=0.0))["e"]
    assert read("ms2", spec=_flat(["plain g1"], _const(0.0), ms=-0.01))["e"]


def test_guards_void_the_reading_when_a_base_is_off(tmp_path, monkeypatch):
    monkeypatch.setattr(k18, "dm_check", lambda rows, *a: {})
    runs, refs, x2 = _runs(tmp_path / "a")
    body = json.loads(x2.read_text())
    body["rows"][3]["scores"]["raw"] += 0.6  # one seed off by more than 0.5
    x2.write_text(json.dumps(body))
    monkeypatch.setattr(k25, "K18_SHA256", {f: _sha(refs / f) for f in ("plain_g1.json", "plain_g2.json", "dm.json")})
    monkeypatch.setattr(k25, "K9_X2_SHA256", _sha(x2))
    res = k25.reading(runs, refs, x2)
    assert res["void"] and res["tries_counted"] == 2 and not res["guards"]["x2"]["passes"]
    assert res["outcome"] == k25.VOID_TEXT and "passes" not in res
    p = refs / "dm.json"
    b = json.loads(p.read_text())
    for q in b["rows"]["g2"]:
        q["scores"]["barred_form"] += 0.15  # the mean off by more than 0.1
    p.write_text(json.dumps(b))
    monkeypatch.setattr(k25, "K18_SHA256", {f: _sha(refs / f) for f in ("plain_g1.json", "plain_g2.json", "dm.json")})
    res = k25.reading(runs, refs, x2)
    assert res["void"] and not res["guards"]["dm g2"]["passes"]
    monkeypatch.setattr(k25, "K9_X2_SHA256", "0" * 64)  # a pinned file off its sha256
    assert any("K9's committed x2" in w for w in k25.reading(runs, refs, x2)["why"])


def test_void_on_missing_seed_dirty_tree_other_commit_or_machine(read):
    res = read("miss", drop="dm g1 4")
    assert res["void"] and res["tries_counted"] == 2 and any("missing seed file" in w for w in res["why"])
    dirty = read("dirty", prov={("plain g1", 2): PROV | {"tree_clean": False}})
    assert dirty["void"] and any("dirty" in w for w in dirty["why"])
    other = read("commit", prov={("x2", 0): PROV | {"commit": "d" * 40}})
    assert other["void"] and any("different commits" in w for w in other["why"])
    mach = read("mach", prov={("dm g2", 9): PROV | {"machine": "local"}})
    assert mach["void"] and any("different machines" in w for w in mach["why"])


def test_a_seed_file_of_another_pregate_is_void(tmp_path, monkeypatch):
    monkeypatch.setattr(k18, "dm_check", lambda rows, *a: {})
    runs, refs, x2 = _runs(tmp_path / "p")
    (runs / "pregate.json").write_text(json.dumps({"passes": True, "figures": {}, "ratios": {}, "provenance": PROV,
                                                   "other": 1}))
    monkeypatch.setattr(k25, "K18_SHA256", {f: _sha(refs / f) for f in ("plain_g1.json", "plain_g2.json", "dm.json")})
    monkeypatch.setattr(k25, "K9_X2_SHA256", _sha(x2))
    res = k25.reading(runs, refs, x2)
    assert res["void"] and any("another pregate.json" in w for w in res["why"])


def test_no_pregate_is_not_a_reading_and_a_failed_pregate_counts_nothing(tmp_path, read):
    (tmp_path / "empty").mkdir()
    res = k25.reading(tmp_path / "empty")
    assert res["incomplete"] and res["tries_counted"] == 0
    failed = read("fail", pre={"passes": False})
    assert not failed["void"] and failed["tries_counted"] == 0 and "pre-gate failed" in failed["outcome"]


def test_the_reading_text_holds_no_error_and_the_tables(read):
    for res in (read("t1"), read("t2", drop="x2 0"), read("t3", pre={"passes": False})):
        t = k25.text(res)
        assert t.startswith("Try K25: ")


def test_read_never_overwrites_a_reading(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "reading.txt").write_text("x")
    with pytest.raises(SystemExit) as e:
        k25.main(["read", "--commit", "0" * 40, "--runs", str(tmp_path), "--out", str(out)])
    assert e.value.code == k25.EXIT_EXISTS


# ---------------------------------------------------------------------------------------------------------------------
# Preconditions


def _row_text() -> str:
    return ROW_FILE.read_text().strip()


def test_the_handed_row_has_eleven_cells_no_inner_pipe_and_only_two_placeholders():
    row = _row_text()
    assert "\n" not in row and row.startswith("| K25 (") and row.endswith(" |")
    cells = row[2:-2].split(" | ")
    assert len(cells) == 11 and row.count("|") == 12
    assert cells[1] == "T3" and cells[9] == "<this row's commit>" and cells[10] == "<outcome>"
    assert re.findall(r"<[^>]*>", row) == ["<this row's commit>", "<outcome>"]
    assert k25.ROW_SCRIPT in cells[2]


def test_the_row_states_every_rule_the_code_fixes(tmp_path):
    # the reviewer's blocker: the row states AHEAD, BEHIND, the gated states, "skill up", "in each source", the
    # estimator's constants and the rounding seed; the script refuses under a row that leaves any of them out
    row = _row_text()
    for phrase in k25.ROW_MUST_STATE:
        assert phrase in row, phrase
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n\n{row}\n")
    assert k25.row_present(p) == "K25"
    for phrase in k25.ROW_MUST_STATE:
        p.write_text(row.replace(phrase, "") + "\n")
        with pytest.raises(SystemExit) as e:
            k25.row_present(p)
        assert e.value.code == k25.EXIT_REFUSED, phrase
    # the first engineer's row text (no AHEAD definition, no gated states, no estimator constants) is refused
    old = ("| K25 (one row) | T3 | `experiments/t3_decomp/k25_dispersion.py` | Pre-gate, source-only, before any score: "
           "at least 50 genes with abs(ln r) > ln 1.2. PASS: AHEAD on x2 and in at least one dm and one plain cell |\n")
    p.write_text(old)
    with pytest.raises(SystemExit) as e:
        k25.row_present(p)
    assert e.value.code == k25.EXIT_REFUSED


def test_the_rows_numbers_are_the_codes():
    row = _row_text()
    assert f"default_rng({k25.ROUND_SEED})" in row and f"seed {k25.PERM_SEED}, {k25.PERM_BINS} bins" in row
    assert f"at least {k25.MIN_LABEL_CELLS} cells" in row and f"n_g / (n_g + {k25.SHRINK:.0f})" in row
    assert f"at least {k25.DETECT:.0%} of" in row and f"[{k25.R_CLIP[0]}, {k25.R_CLIP[1]:.0f}]" in row
    assert f"at least {k25.PREGATE_MIN_GENES} genes" in row and k25.MEAN_TOL == 1e-9 and "equal to 1e-9 before" in row
    assert f"mean > {k25.MIN_GAIN} with at least {k25.MIN_AHEAD} of 10" in row and f"mean < -{k25.MIN_GAIN}" in row
    assert f"plus {k25.MAX_EXTRA_GATED}" in row
    assert f"mean within {k25.GUARD_MEAN}, every seed within {k25.GUARD_SEED}" in row
    m = re.search(r"gated: ([a-z0-9_, ]+) \(severity states\)", row)
    assert m and sorted(m.group(1).split(", ")) == sorted(SEVERITY_GATED)
    ms = " + ".join(f"{w} {m}" for m, w in k25.MUTANT_STATE)
    assert f"{ms} skill is above 0" in row and "two tries" in row and k25.TRIES == 2
    for sha in (k25.K9_MODULE_SHA256, k25.K9_X2_SHA256, *k25.K18_SHA256.values()):
        assert sha[:12] in row
    assert f"exit {k25.EXIT_PREGATE_FAILS}, launcher 25" in row and f"exit {k25.EXIT_UNCOVERED}, launcher 26" in row


def test_the_row_must_be_k25_and_alone(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text(_row_text() + "\n")
    assert k25.row_present(p) == "K25"
    for text_ in ("", "| K26 (x) | T3 | `experiments/t3_decomp/k25_dispersion.py` |\n",
                  "| K25 (x) | T3 | `experiments/t3_decomp/k25_dispersion.py` |\n"
                  "| K26 (y) | T3 | `experiments/t3_decomp/k25_dispersion.py` |\n"):
        p.write_text(text_)
        with pytest.raises(SystemExit) as e:
            k25.row_present(p)
        assert e.value.code == k25.EXIT_REFUSED


def test_every_subcommand_refuses_another_commit_or_no_row(tmp_path, monkeypatch):
    monkeypatch.setattr(k25, "TRIES_MD", tmp_path / "TRIES.md")
    (tmp_path / "TRIES.md").write_text(_row_text() + "\n")
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "a" * 40, "tree_clean": True, "machine": "m"})
    for argv in (["pregate", "--out", str(tmp_path / "o")], ["run", "--world", "x2", "--out", str(tmp_path / "o")],
                 ["read", "--runs", str(tmp_path), "--out", str(tmp_path / "o2")]):
        with pytest.raises(SystemExit) as e:
            k25.main([*argv, "--commit", "b" * 40])
        assert e.value.code == k25.EXIT_REFUSED
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "b" * 40, "tree_clean": False, "machine": "m"})
    with pytest.raises(SystemExit) as e:
        k25.main(["pregate", "--out", str(tmp_path / "o"), "--commit", "b" * 40])
    assert e.value.code == k25.EXIT_REFUSED and not (tmp_path / "o").exists()
    (tmp_path / "TRIES.md").write_text("")
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "b" * 40, "tree_clean": True, "machine": "m"})
    with pytest.raises(SystemExit) as e:
        k25.main(["pregate", "--out", str(tmp_path / "o"), "--commit", "b" * 40])
    assert e.value.code == k25.EXIT_REFUSED


def test_a_run_refuses_a_missing_failing_or_foreign_pregate(tmp_path):
    st = {"commit": "a" * 40, "tree_clean": True, "machine": "m"}
    with pytest.raises(SystemExit):
        k25.load_pregate(tmp_path, st)
    for body in ({"passes": False, "provenance": st, "ratios": {}},
                 {"passes": True, "provenance": st | {"commit": "b" * 40}, "ratios": {}},
                 {"passes": True, "provenance": st | {"machine": "x"}, "ratios": {}}):
        (tmp_path / "pregate.json").write_text(json.dumps(body))
        with pytest.raises(SystemExit) as e:
            k25.load_pregate(tmp_path, st)
        assert e.value.code == k25.EXIT_REFUSED
    (tmp_path / "pregate.json").write_text(json.dumps({"passes": True, "provenance": st, "ratios": {
        "mab21l2": {"genes": ["a"], "r": [1.5], "record": {}}}}))
    assert len(k25.load_pregate(tmp_path, st)) == 64 and k25._DISP["mab21l2"]["r"] == [1.5]
    k25._DISP.clear()


def test_the_script_prints_step_lines_only():
    src = (T3 / "k25_dispersion.py").read_text()
    allowed = {"msg", "path", "prov['commit']", "prov['tree_clean']", "a.task", "s", "part", "out",
               "m.group(1) if m else '?'", "out / 'reading.json'",
               "'passes' if passes else 'fails; nothing runs and nothing counts'"}
    calls = re.findall(r'print\((f"[^"]*")', src)
    exprs = {e for c in calls for e in re.findall(r"\{([^{}]+)\}", c)}
    assert calls and exprs <= allowed, exprs - allowed
    assert not re.search(r"print\((?!f\"|\"|k18\.masked)", src)  # nothing else is printed


# ---------------------------------------------------------------------------------------------------------------------
# The launcher


def _launch(home: Path, *argv: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(home)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO, env=env, capture_output=True, text=True)


def test_the_launcher_deletes_nothing_and_quotes_nothing_in_a_usage_message():
    sh = LAUNCHER.read_text()
    code = [line.split("#")[0] for line in sh.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py|experiments/t3_decomp/test_[a-z0-9_]+\.py", sh)
    assert "tests/test_k25_dispersion.py" in listed and all((REPO / f).exists() for f in listed)
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and "OMP_NUM_THREADS=$1" in sh
    assert "experiments/t3_decomp/results/k25_dispersion" in sh


def test_the_launcher_usage_and_unknown_step(tmp_path):
    out = _launch(tmp_path)
    assert out.returncode == 2 and (tmp_path / "v" / "scratch" / "k25" / "usage.exit").read_text().strip() == "2"
    out = _launch(tmp_path, "bogus", "0" * 40)
    assert out.returncode == 97 and (tmp_path / "v" / "scratch" / "k25" / "bogus.exit").read_text().strip() == "97"


def test_each_step_refuses_unless_the_step_before_exited_0_at_the_same_commit(tmp_path):
    S = tmp_path / "v" / "scratch" / "k25"
    for step in ("pregate", "run", "read"):
        out = _launch(tmp_path, step, "0" * 40)
        assert out.returncode == 96 and (S / f"{step}.exit").read_text().strip() == "96", (step, out.stderr)
    (S / "tests.exit").write_text("0\n")
    (S / "tests.head").write_text("1" * 40 + "\n")  # exit 0, but at another commit
    out = _launch(tmp_path, "pregate", "0" * 40)
    assert out.returncode == 96 and (S / "pregate.exit").read_text().strip() == "96"
    (S / "tests.head").write_text("0" * 40 + "\n")
    (S / "tests.exit").write_text("10\n")  # the same commit, a failed step
    assert _launch(tmp_path, "pregate", "0" * 40).returncode == 96


def test_a_wrong_commit_is_refused_with_an_exit_file_and_nothing_run(tmp_path):
    S = tmp_path / "v" / "scratch" / "k25"
    out = _launch(tmp_path, "tests", "0" * 40)
    assert out.returncode == 99 and "HEAD is not" in out.stderr
    assert (S / "tests.exit").read_text().strip() == "99" and not (S / "tests.log").exists()
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    out = _launch(tmp_path, "tests", head[:12])  # a prefix of the right commit is never enough
    assert out.returncode == 99


def test_the_launcher_at_head_refuses_without_the_row(tmp_path):
    tries = (REPO / "experiments" / "harness" / "TRIES.md").read_text()
    if re.search(r"^\| K25 .*`experiments/t3_decomp/k25_dispersion\.py`", tries, flags=re.M):
        pytest.skip("the row is placed: this path would start the tests step itself")
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    out = _launch(tmp_path, "tests", head)
    S = tmp_path / "v" / "scratch" / "k25"
    assert out.returncode in (95, 98), out.stderr  # 98 only while the tree holds uncommitted changes
    assert (S / "tests.exit").read_text().strip() == str(out.returncode) and not (S / "tests.log").exists()


def test_the_launchers_pregate_passes_check(tmp_path):
    m = re.search(r"^PASSES_CHECK='([^']+)'$", LAUNCHER.read_text(), flags=re.M)
    assert m
    c = "c" * 40
    prov = {"provenance": {"commit": c}}
    for body, want in (({"passes": True} | prov, "pass"), ({"passes": False} | prov, "fail"), (prov, "fail"),
                       ({"passes": "true"} | prov, "fail"), ({"passes": True}, "foreign"),
                       ({"passes": True, "provenance": {"commit": "d" * 40}}, "foreign"),
                       ({"passes": False, "provenance": {"commit": "d" * 40}}, "foreign")):
        p = tmp_path / "p.json"
        p.write_text(json.dumps(body))
        got = subprocess.run([sys.executable, "-I", "-B", "-c", m.group(1), str(p), c], capture_output=True,
                             text=True)
        assert got.stdout.strip() == want, body


def test_the_launcher_maps_each_script_exit():
    sh = LAUNCHER.read_text()
    pre = sh[sh.index("  pregate)\n    if"):sh.index("  run)\n    mkdir")]
    assert '"$PASSES_CHECK" "$R/pregate.json" "$COMMIT"' in pre
    assert 'if [ "$got" = foreign ]; then finish 21; fi' in pre
    for rc, code in ((0, 0), (5, 25), (6, 26), (2, 22)):
        assert f"if [ $rc = {rc} ]; then finish {code}; fi" in pre
    assert pre.rstrip().endswith("finish 20 ;;")
    rd = sh[sh.index("  read)\n    if"):]
    for rc, code in ((0, 0), (3, 43), (4, 44), (2, 42)):
        assert f"if [ $rc = {rc} ]; then finish {code}; fi" in rd
    assert rd.rstrip().endswith("finish 41 ;;\nesac")
    assert (k25.EXIT_REFUSED, k25.EXIT_INCOMPLETE, k25.EXIT_EXISTS, k25.EXIT_PREGATE_FAILS, k25.EXIT_UNCOVERED) == (
        2, 3, 4, 5, 6)
