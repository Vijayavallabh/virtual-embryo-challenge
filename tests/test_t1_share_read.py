"""experiments/t1_m1/share_read.py (try T47 of experiments/harness/TRIES.md) on synthetic statistics and synthetic
AnnData stages built here. No released file, statistics file or score is read."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "share_read.py"
EPS = 0.5
N_CELLS = 21  # per label and library: above the config's min_cells_t1
FENG_LABELS = ("ventricular_cm", "endocardial_ec", "vascular_ec", "epicardial")
T1_LABEL = {"CM": "V-CM", "Endocardium": "Endocardium", "Endothelium": "BEC", "Epicardium": "Epicardium"}


@pytest.fixture
def sr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("share_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


@pytest.fixture
def cfg(sr):
    return sr.feng_dirs.load_config()


# ---------------------------------------------------------------------------------------------------------------------
# Synthetic inputs


def feng_stats(genes, means: dict, stages=("E9.5", "E14.5"), n=30) -> dict:
    """Feng-format statistics ('stage|zone|cell_type_spec|half'): means[(stage, label)] is the mean CP10k vector, or
    means[(stage, label, half)] a half's own (the two halves have equal cells, so the pooled mean is their average)."""
    keys, ns, S = [], [], []
    for st in stages:
        for lab in FENG_LABELS:
            for h in "AB":
                keys.append(f"{st}|z|{lab}|{h}")
                ns.append(n)
                S.append(means.get((st, lab, h), means[(st, lab)]) * n)
    S = np.array(S, dtype=np.float64)
    return {"genes": np.array(genes, dtype=object), "keys": np.array(keys, dtype=object), "n": np.array(ns, float),
            "S_log": np.log1p(S / n) * n, "S_cp": S, "S_nz": np.full_like(S, n),
            "depth_sum": np.full(len(keys), 1e4 * n)}


def stage(blocks: dict, genes, sizes: dict | None = None) -> ad.AnnData:
    """blocks[(label, library)] = mean CP10k vector; N_CELLS (or sizes[(label, library)]) identical cells each, named
    cell<i>_<library>."""
    rows, labs, names = [], [], []
    i = 0
    for (lab, lib), cp in blocks.items():
        for _ in range((sizes or {}).get((lab, lib), N_CELLS)):
            rows.append(np.log1p(cp))
            labs.append(lab)
            names.append(f"c{i}_{lib}")
            i += 1
    a = ad.AnnData(X=sparse.csr_matrix(np.array(rows, dtype=np.float32)),
                   obs=pd.DataFrame({"celltype": labs}, index=names), var=pd.DataFrame(index=list(genes)))
    return a


def world(seed=0, G=60, extra_stage=False, half_noise=0.0):
    """Feng statistics at E9.5 and E14.5 (and E12.5 with extra_stage) on genes g1..; gene 'g0' is absent from Feng (a
    non-finite direction). half_noise: the E9.5 halves differ (A = a (1 + e), B = a (1 - e), e ~ N(0, half_noise)),
    their pooled mean staying a."""
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(G)]
    feng_genes = genes[1:]
    means = {}
    for lab in FENG_LABELS:
        a = rng.uniform(2.0, 6.0, G - 1)
        means[("E9.5", lab)] = a
        if half_noise:
            e = np.clip(rng.normal(0, half_noise, G - 1), -0.9, 0.9)
            means[("E9.5", lab, "A")], means[("E9.5", lab, "B")] = a * (1 + e), a * (1 - e)
        means[("E14.5", lab)] = (a + EPS) * np.exp(rng.normal(0, 0.3, G - 1)) - EPS
        if extra_stage:
            means[("E12.5", lab)] = np.full(G - 1, 1e3)
    stages_f = ("E9.5", "E12.5", "E14.5") if extra_stage else ("E9.5", "E14.5")
    fst = feng_stats(feng_genes, means, stages_f)
    return genes, fst


def build_t1(sr, cfg, genes, fst, c_lib=(0.2, 0.3, 0.9), c_prev=0.25, k=1.5, held_gene=None, c_group=None,
             n_lib3=None):
    """Three T1 stages where, for every group, E9.5 -> E10.5 library i moves by c_lib[i] (times c_group[g]) times the
    group's Feng direction, E8.5 -> E9.5 by c_prev times it (CM; at E8.5 Endocardium takes Endothelium's cells and
    Epicardium has none), and the E9.5 contrasts between groups are k times Feng's E9.5 contrasts. n_lib3[g]: the
    group's cells in E10.5's library 3 (default N_CELLS)."""
    c_group = c_group or {}
    dirs, _, idx = sr.feng_dirs.directions_from_stats(fst, genes, cfg, sr.BRACKET)
    e95 = {}
    for g in T1_LABEL:
        fl = cfg["external_labels"][g][0]
        fm = sr.dl.feng_sel(fst, "E9.5", fl)["cp"]
        e95[g] = np.exp(k * np.log(sr.dl.on_t1(fm, idx) + EPS)) - EPS
        e95[g][~np.isfinite(e95[g])] = 5.0
    d0 = {g: np.nan_to_num(dirs[g]) for g in dirs}

    def moved(g, c):
        v = (e95[g] + EPS) * np.exp(c * d0[g]) - EPS
        if held_gene is not None:
            v[held_gene] = e95[g][held_gene] * 50.0  # a large change on a held gene, which must not enter
        return v

    E95 = {(T1_LABEL[g], li): e95[g] for g in T1_LABEL for li in ("1", "2")}
    E95[("Foregut", "1")] = np.full(len(genes), 3.0)
    E105 = {(T1_LABEL[g], li): moved(g, c * c_group.get(g, 1.0)) for g in T1_LABEL
            for li, c in zip(("1", "2", "3"), c_lib, strict=True)}
    sizes = {(T1_LABEL[g], "3"): n for g, n in (n_lib3 or {}).items()}
    E85 = {("V-CM", "1"): (e95["CM"] + EPS) * np.exp(-c_prev * d0["CM"]) - EPS,
           ("BEC", "1"): e95["Endothelium"], ("Foregut", "1"): np.full(len(genes), 3.0)}
    return {"E8.5": stage(E85, genes), "E9.5": stage(E95, genes), "E10.5": stage(E105, genes, sizes)}, dirs


# ---------------------------------------------------------------------------------------------------------------------
# The share


def test_a_change_made_as_c_times_the_direction_has_share_c_and_c_over_scale_corrected(sr):
    d = np.random.default_rng(1).normal(size=200)
    m = np.ones(200, bool)
    st = sr.share_stats(0.37 * d, d, m, 0.8)
    assert st["share"] == pytest.approx(0.37) and st["share_corrected"] == pytest.approx(0.37 / 0.8)
    assert st["pearson"] == pytest.approx(1.0) and st["slope_intercept"] == pytest.approx(0.37)
    assert st["share_top_half"] == pytest.approx(0.37)


def test_a_change_orthogonal_to_the_direction_has_share_zero(sr):
    d = np.random.default_rng(2).normal(size=100)
    v = np.random.default_rng(3).normal(size=100)
    orth = v - (v @ d) / (d @ d) * d
    assert sr.share_stats(orth, d, np.ones(100, bool), 1.0)["share"] == pytest.approx(0.0, abs=1e-12)


def test_the_top_half_takes_the_genes_with_the_larger_absolute_direction(sr):
    d = np.array([3.0, -2.0, 0.1, 0.2])
    change = np.array([3.0, -2.0, 5.0, -5.0])  # share 1 on the two large |d| genes; the small ones disagree
    assert sr.share_stats(change, d, np.ones(4, bool), None)["share_top_half"] == pytest.approx(1.0)
    assert np.isnan(sr.share_stats(change, d, np.ones(4, bool), None)["share_corrected"])


def test_held_genes_and_genes_without_a_finite_direction_do_not_enter(sr, cfg):
    genes, fst = world()
    stages, dirs = build_t1(sr, cfg, genes, fst, held_gene=5)
    clean = sr.reading(stages, fst, [genes[5]], {}, cfg, n_boot=2)
    for g in clean["groups"]:
        assert clean["genes_per_group"][g] == len(genes) - 2  # g0 has no Feng value, g5 is held
        assert clean["values"][f"b_E9.5_E10.5_lib1/{g}/share"]["value"] == pytest.approx(0.2, abs=1e-4)
    dirty = sr.reading(stages, fst, [], {}, cfg, n_boot=2)
    assert dirty["values"]["b_E9.5_E10.5_lib1/CM/share"]["value"] != pytest.approx(0.2, abs=1e-3)


def test_the_library_selection_uses_the_cell_name_suffix_and_the_reading_recovers_A_K_P(sr, cfg):
    genes, fst = world()
    stages, _ = build_t1(sr, cfg, genes, fst, c_lib=(0.2, 0.3, 0.9), c_prev=0.25)
    r = sr.reading(stages, fst, [], {"CM": 0.5}, cfg, n_boot=3)
    v = r["values"]
    for g in r["groups"]:
        for name, c in (("lib1", 0.2), ("lib2", 0.3), ("lib3", 0.9)):
            assert v[f"b_E9.5_E10.5_{name}/{g}/share"]["value"] == pytest.approx(c, abs=1e-4)
        assert v[f"c_E10.5_lib1_lib2/{g}/share"]["value"] == pytest.approx(0.1, abs=1e-4)
    # identical Feng halves: r_d = 1, so the corrected share is the share; N48's 0.5 is reported only
    assert r["r_d"]["CM"] == pytest.approx(1.0) and r["n48_scale"]["CM"] == 0.5
    assert v["b_E9.5_E10.5_lib3/CM/share_corrected"]["value"] == pytest.approx(0.9, abs=1e-4)
    assert v["b_E9.5_E10.5_lib3/Endocardium/share_corrected"]["value"] == pytest.approx(0.9, abs=1e-4)
    rd = r["reading"]
    assert rd["A"]["value"] == pytest.approx(0.25, abs=1e-4)
    assert rd["K"]["value"] == pytest.approx(0.65, abs=1e-4)
    assert v["a_E8.5_E9.5/CM/share"]["value"] == pytest.approx(0.25, abs=1e-4)
    # E8.5 has no Endocardium label (fallback to Endothelium's cells) and no Epicardium (no value)
    assert r["label_source"]["E8.5"]["Endocardium"] == "fallback:Endothelium"
    assert r["label_source"]["E8.5"]["Epicardium"] == "none"
    assert np.isnan(v["a_E8.5_E9.5/Epicardium/share"]["value"])
    assert r["label_source"]["E9.5"]["CM"] == "own"


def test_the_scale_is_recovered_when_the_challenge_contrasts_are_a_multiple_of_feng(sr, cfg):
    genes, fst = world()
    stages, _ = build_t1(sr, cfg, genes, fst, k=1.5)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=2)
    assert r["reading"]["L"]["value"] == pytest.approx(1.5, abs=1e-4)
    for pair in r["genes_per_pair"]:
        assert r["values"][f"scale_E9.5/{pair}/slope"]["value"] == pytest.approx(1.5, abs=1e-4)
        assert r["values"][f"scale_E9.5/{pair}/pearson"]["value"] == pytest.approx(1.0, abs=1e-5)


def test_the_combined_shares_weight_by_library_3_cells_and_skip_missing_groups(sr):
    wtd, unw = sr.combine({"a": 1.0, "b": 3.0, "c": float("nan")}, {"a": 3.0, "b": 1.0, "c": 5.0})
    assert wtd == pytest.approx(1.5) and unw == pytest.approx(2.0)
    wtd, unw = sr.combine({"a": 1.0, "b": 3.0}, {"a": 0.0, "b": 0.0})
    assert np.isnan(wtd) and unw == pytest.approx(2.0)


# ---------------------------------------------------------------------------------------------------------------------
# The bootstrap


def _noisy_stage():
    rng = np.random.default_rng(4)
    labs = ["A"] * 7 + ["A"] * 5 + ["B"] * 6
    libs = ["1"] * 7 + ["2"] * 5 + ["1"] * 6
    X = sparse.csr_matrix(np.log1p(rng.uniform(0, 5, (18, 4))).astype(np.float32))
    return ad.AnnData(X=X, obs=pd.DataFrame({"celltype": labs}, index=[f"c{i}_{li}" for i, li in enumerate(libs)]),
                      var=pd.DataFrame(index=["w", "x", "y", "z"]))


def test_the_bootstrap_is_seeded_and_resamples_within_label_stage_and_library(sr):
    a = _noisy_stage()
    tab = sr.stage_table(a, list(a.var_names), {"A", "B"})
    assert sorted(tab["keys"]) == [("A", "1"), ("A", "2"), ("B", "1")]
    w1 = sr.resample_weights(tab, np.random.default_rng([47, 0]))
    w2 = sr.resample_weights(tab, np.random.default_rng([47, 0]))
    w3 = sr.resample_weights(tab, np.random.default_rng([47, 1]))
    assert np.array_equal(w1, w2) and not np.array_equal(w1, w3)
    for b, members in enumerate(tab["members"]):  # every block draws its own size from its own cells
        assert w1[members].sum() == len(members)
        assert set(np.flatnonzero(tab["inv"] == b)) == set(members)
    assert w1.sum() == a.n_obs
    sums = sr.block_sums(tab, w1)
    one = sr.cell_set(tab, sums, ["A"], ("1",), 0)
    m = tab["members"][tab["keys"].index(("A", "1"))]
    cp = np.expm1(a.X.toarray().astype(np.float32)).astype(np.float64)
    assert one["n"] == 7 and np.allclose(one["cp"], (w1[m, None] * cp[m]).sum(0) / 7, atol=1e-5)
    assert sr.cell_set(tab, sums, ["A"], ("1",), 8) is None


def test_the_reading_bootstrap_is_seeded(sr, cfg):
    genes, fst = world()
    stages, _ = build_t1(sr, cfg, genes, fst)
    rng = np.random.default_rng(9)
    for a in stages.values():  # cell-level noise so that draws differ
        a.X = sparse.csr_matrix(np.log1p(np.expm1(a.X.toarray()) * rng.uniform(0.5, 1.5, a.shape)).astype(np.float32))
    r1 = sr.reading(stages, fst, [], {}, cfg, n_boot=4, seed=(47, 0))
    r2 = sr.reading(stages, fst, [], {}, cfg, n_boot=4, seed=(47, 0))
    r3 = sr.reading(stages, fst, [], {}, cfg, n_boot=4, seed=(47, 1))
    assert json.dumps(r1["values"], sort_keys=True) == json.dumps(r2["values"], sort_keys=True)
    assert r1["values"]["reading/A"]["p5"] != r3["values"]["reading/A"]["p5"]


# ---------------------------------------------------------------------------------------------------------------------
# The courses and the statements


@pytest.mark.parametrize(("ratio", "f_max", "f_even", "f_decay"),
                         [(0.2, 0.8, 0.4, 0.4), (0.3, 0.7, 0.35, 0.4285), (0.5, 0.5, 0.25, 0.3939)])
@pytest.mark.parametrize("L", [1.0, 2.0])
def test_the_courses_by_hand(sr, ratio, f_max, f_even, f_decay, L):
    c = sr.courses(ratio * L, L)
    assert c["f_max"] == pytest.approx(f_max * L)
    assert c["f_even"] == pytest.approx(f_even * L)
    assert c["f_decay"] == pytest.approx(f_decay * L, abs=1e-3 * L)
    x = c["x"]
    if abs(x - 1) > 1e-6:  # the brief's form, where it is defined
        assert (1 - x) / (1 - x ** 5) == pytest.approx(ratio)
        assert c["f_decay"] == pytest.approx(L * (x - x ** 3) / (1 - x ** 5))
    else:
        assert ratio == pytest.approx(0.2)


def test_the_decay_course_at_the_ends_of_its_range(sr):
    # no decaying course where the day holds nothing or goes back; nothing left where it holds the whole change or more
    assert np.isnan(sr.courses(-0.1, 1.0)["f_decay"]) and np.isnan(sr.courses(0.0, 1.0)["f_decay"])
    assert sr.courses(1.0, 1.0)["f_decay"] == 0.0 and sr.courses(1.2, 1.0)["f_decay"] == 0.0
    assert sr.courses(2.4, 2.0)["x"] == 0.0
    assert sr.courses(1e-4, 1.0)["f_decay"] < 0.02  # f_decay -> 0 as A / L -> 0+
    assert np.isnan(sr.courses(0.3, 0.0)["A_over_L"])


def test_each_kind_statement_is_reached(sr):
    assert sr.kind_statement(0.15, 0.4) == sr.KIND_ALONG
    assert sr.kind_statement(-0.2, 0.05) == sr.KIND_NOT
    assert sr.kind_statement(0.05, 0.2) == sr.NOT_SEPARATED


def test_each_pace_statement_is_reached(sr):
    assert sr.pace_statement(0.21, 0.5) == sr.PACE_MORE
    assert sr.pace_statement(0.15, 0.25) == sr.PACE_FIFTH_OR_LESS  # includes 0.2
    assert sr.pace_statement(0.05, 0.15) == sr.PACE_FIFTH_OR_LESS  # below 0.2
    assert sr.pace_statement(0.15, 0.35) == sr.NOT_SEPARATED


def _rec(v, lo, hi):
    return {"value": v, "p5": lo, "p95": hi}


def test_each_fraction_statement_is_reached(sr):
    big = _rec(1.0, 0.8, 1.2)
    assert sr.fraction_statement(_rec(0.3, 0.2, 0.35), _rec(0.15, 0.1, 0.2), _rec(0.2, 0.1, 0.3)) == sr.ABOVE_BOUND
    assert sr.fraction_statement(big, _rec(0.35, 0.3, 0.38), _rec(0.42, 0.39, 0.45)) == sr.INSIDE
    assert sr.fraction_statement(big, _rec(0.25, 0.2, 0.3), _rec(0.3, 0.25, 0.35)) == sr.ABOVE_BOTH
    assert sr.fraction_statement(big, _rec(0.5, 0.45, 0.55), _rec(0.6, 0.5, 0.7)) == sr.BELOW_BOTH
    # the range runs from the lower of the two lower percentiles to the higher of the two upper ones, whatever the
    # point values: here f_even's point is the smaller but its lower percentile lies above f_decay's upper one
    assert sr.fraction_statement(big, _rec(0.3, 0.38, 0.5), _rec(0.35, 0.2, 0.37)) == sr.INSIDE
    assert sr.fraction_statement(big, _rec(0.6, 0.2, 0.7), _rec(0.3, 0.41, 0.45)) == sr.INSIDE
    assert sr.fraction_statement(big, _rec(0.6, 0.45, 0.7), _rec(0.3, 0.41, 0.45)) == sr.BELOW_BOTH


def test_the_statements_come_out_of_the_reading(sr, cfg):
    genes, fst = world()
    stages, _ = build_t1(sr, cfg, genes, fst, c_lib=(0.18, 0.18, 0.6), k=1.0)
    r = sr.reading(stages, fst, [], dict.fromkeys(T1_LABEL, 1.0), cfg, n_boot=2)  # r_d = r_c = 1: corrected = raw
    assert r["statements"]["pace_raw"] == r["statements"]["pace_corrected"]
    assert r["statements"]["kind"] == sr.KIND_ALONG  # K = 0.42
    assert r["statements"]["pace"] == sr.PACE_FIFTH_OR_LESS  # A / L = 0.18
    assert r["statements"]["entry_fraction"] == sr.INSIDE  # f_decay 0.389, f_even 0.41
    txt = sr.text({"label": sr.LABEL_T47, "read": r, "stamp": {"commit": "abc", "tree_clean": True, "utc": "t"}})
    assert sr.KIND_ALONG in txt and sr.PACE_FIFTH_OR_LESS in txt and "Entry's fraction 0.4: inside" in txt
    assert "(0 of 2 draws undefined)" in txt  # beside each reading value
    assert "over ['CM', 'Endocardium', 'Endothelium', 'Epicardium']" in txt  # each combined value's groups
    assert "pairs with shared Feng labels excluded: ['Endocardium~Endothelium']" in txt
    assert all(note in txt for note in sr.LEDGER_NOTES)
    json.dumps(r)  # the record serialises (no tuple keys)
    r["calibration"] = {"fraction": 0.2, "combined_weighted": 0.19, "combined_unweighted": 0.18,
                        "combined_groups": ["CM"], "A_in_M1_fraction": 0.19, "K_in_M1_fraction": 0.44,
                        "groups": {"CM": {"share": 0.19}}}
    txt = sr.text({"label": sr.LABEL_T47, "read": r, "n48_js": {"CM": 0.7},
                   "stamp": {"commit": "abc", "tree_clean": True, "utc": "t"}})
    assert "Calibration: M1's step of 0.2" in txt and "N48 'js' (held genes included)" in txt


# ---------------------------------------------------------------------------------------------------------------------
# Refusals and the Feng stages


def test_refused_without_the_row(sr, monkeypatch):
    monkeypatch.setattr(sr.n12, "tries_rows", lambda *a: {"T44"})
    with pytest.raises(SystemExit, match="no row"):
        sr.run("x")


def test_refused_on_an_existing_output(sr, monkeypatch):
    monkeypatch.setattr(sr.n12, "tries_rows", lambda *a: {"T47"})
    out = sr.RESULTS / "x"
    out.mkdir(parents=True)
    (out / "t47.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        sr.run("x")


def test_refused_on_a_dirty_tree(sr, monkeypatch):
    monkeypatch.setattr(sr.n12, "tries_rows", lambda *a: {"T47"})
    monkeypatch.setattr(sr.gv, "tree_state", lambda: {"commit": "abc", "clean": False})
    monkeypatch.setattr(sr.feng_dirs, "load_config", lambda *a: pytest.fail("read past the refusal"))
    with pytest.raises(SystemExit, match="not clean"):
        sr.run("x")


def test_no_feng_stage_other_than_the_bracket_ends_is_read(sr, cfg, monkeypatch, tmp_path):
    genes, fst = world(extra_stage=True)
    monkeypatch.setattr(sr.fio, "verify_manifest_output", lambda p: {"sha256": "0" * 64})
    monkeypatch.setattr(sr.dl, "load_stats", lambda p, fields=None: fst)
    handed, rec = sr.feng_inputs(tmp_path / "cd1_stats.npz")
    assert {k.split("|")[0] for k in handed["keys"]} == set(sr.BRACKET)
    assert {k.split("|")[0] for k in fst["keys"]} == {"E9.5", "E12.5", "E14.5"}
    seen = []
    real = sr.dl.feng_sel

    def spy(st, stage, labels, halves=("A", "B")):
        seen.append((stage, tuple(halves)))
        return real(st, stage, labels, halves)

    monkeypatch.setattr(sr.dl, "feng_sel", spy)
    stages, _ = build_t1(sr, cfg, genes, handed)
    sr.reading(stages, handed, [], {}, cfg, n_boot=1)
    assert seen and {s for s, _ in seen} <= set(sr.BRACKET)
    # r_d's halves are read at both bracket ends, and at no other stage
    assert {s for s, h in seen if h == ("A",)} == set(sr.BRACKET)


# ---------------------------------------------------------------------------------------------------------------------
# Review fixes: corrected units, undefined draws, common support, weights, pinned choices, calibration


def test_the_contrast_reliability_is_n48s_formula_on_the_contrast(sr):
    c = np.array([1.0, -2.0, 0.5, np.nan, 3.0])
    c_a = np.array([1.4, -2.2, 0.1, 1.0, np.nan])
    c_b = np.array([0.6, -1.8, 0.9, 1.0, 1.0])
    m = np.ones(5, bool)
    ok = np.array([0, 1, 2])
    msq = np.mean(c[ok] ** 2)
    noise = np.mean(((c_a[ok] - c_b[ok]) / 2) ** 2)
    assert sr.contrast_reliability(c, c_a, c_b, m) == pytest.approx((msq - noise) / msq)
    assert sr.contrast_reliability(c, c + 10, c - 10, m) == 0.0  # noise above the signal: floored at 0
    assert sr.contrast_reliability(c, c, c, m) == pytest.approx(1.0)
    assert np.isnan(sr.contrast_reliability(c, c_a, c_b, np.zeros(5, bool)))


def test_L_is_corrected_by_each_pairs_feng_reliability_and_drops_pairs_with_shared_feng_labels(sr, cfg):
    genes, fst = world(half_noise=0.08)
    stages, _ = build_t1(sr, cfg, genes, fst, k=1.5)
    r = sr.reading(stages, fst, [], dict.fromkeys(T1_LABEL, 1.0), cfg, n_boot=2)
    assert r["excluded_pairs_shared_feng_labels"] == ["Endocardium~Endothelium"]
    rel = r["feng_reliability"]["E9.5"]
    # by hand for one pair, from Feng's E9.5 halves
    dirs, _, idx = sr.feng_dirs.directions_from_stats(fst, genes, cfg, sr.BRACKET)
    lg, lh = cfg["external_labels"]["CM"][0], cfg["external_labels"]["Epicardium"][0]

    def con(halves):
        fg, fh = sr.dl.feng_sel(fst, "E9.5", lg, halves), sr.dl.feng_sel(fst, "E9.5", lh, halves)
        return sr.dl.on_t1(np.log((fh["cp"] + EPS) / (fg["cp"] + EPS)), idx)

    c, c_a, c_b = con(("A", "B")), con(("A",)), con(("B",))
    ok = np.isfinite(dirs["CM"]) & np.isfinite(dirs["Epicardium"])
    msq = np.mean(c[ok] ** 2)
    assert rel["CM~Epicardium"] == pytest.approx((msq - np.mean(((c_a[ok] - c_b[ok]) / 2) ** 2)) / msq)
    assert 0 < rel["CM~Epicardium"] < 1
    kept = [p for p in r["genes_per_pair"] if p != "Endocardium~Endothelium"]
    v = r["values"]
    assert r["reading"]["L"]["value"] == pytest.approx(1.5, abs=1e-4)
    assert r["reading"]["L_corr"]["value"] == pytest.approx(np.median([1.5 / rel[p] for p in kept]), rel=1e-4)
    for p in kept:
        assert v[f"scale_E9.5/{p}/slope_corrected"]["value"] == pytest.approx(1.5 / rel[p], rel=1e-4)


def test_the_pair_median_is_a_median_over_finite_slopes(sr):
    assert sr.pair_median([1.0, 1.0, 4.0, float("nan")]) == 1.0
    assert np.isnan(sr.pair_median([float("nan")]))


def test_the_corrected_pair_and_the_fractions_in_true_d_units(sr, cfg):
    genes, fst = world(half_noise=0.1)
    stages, _ = build_t1(sr, cfg, genes, fst, c_lib=(0.2, 0.3, 0.9), k=1.0)
    r = sr.reading(stages, fst, [], {"CM": 0.5}, cfg, n_boot=2)
    rd = r["reading"]
    s = np.array([r["r_d"][g] for g in r["groups"]])  # equal library-3 cells: the weighted means are plain means
    assert np.all((s > 0) & (s < 1))
    A_corr = np.mean(0.25 / s)
    assert rd["A"]["value"] == pytest.approx(0.25, abs=1e-4)
    assert rd["A_corr"]["value"] == pytest.approx(A_corr, abs=1e-4)
    assert rd["K_corr"]["value"] == pytest.approx(np.mean(0.9 / s) - A_corr, abs=1e-4)
    assert rd["r_d_weighted"]["value"] == pytest.approx(s.mean())
    L_corr = rd["L_corr"]["value"]
    c = sr.courses(A_corr, L_corr)
    assert rd["A_over_L_corr"]["value"] == pytest.approx(A_corr / L_corr, abs=1e-4)
    for k in sr.COURSE_KEYS:
        assert rd[f"{k}_corr"]["value"] == pytest.approx(c[k], abs=1e-4)  # read by the corrected statement
        assert rd[f"{k}_corr_noisy_d"]["value"] == pytest.approx(c[k] * s.mean(), abs=1e-4)  # reported only


def test_raw_and_corrected_statements_are_joined(sr):
    assert sr.joined(sr.INSIDE, sr.INSIDE) == sr.INSIDE
    assert sr.joined(sr.PACE_MORE, sr.PACE_FIFTH_OR_LESS) == sr.NOT_SEPARATED
    nc = "not computable (3 of 200 draws undefined)"
    assert sr.joined(sr.INSIDE, nc) == nc and sr.joined(nc, sr.INSIDE) == nc


def _draw_records(sr, A_draws, L_draws, A0=0.5, L0=1.0):
    """The reading's records from draws of A and L (K fixed at 0.3; the corrected pair as the raw one)."""
    def one(a, el):
        return {**sr.courses(a, el), "K": 0.3}

    point = one(A0, L0)
    rows = [one(a, el) for a, el in zip(A_draws, L_draws, strict=True)]
    vals = sr.summarise(point, {k: np.array([row[k] for row in rows]) for k in point}, len(rows))
    return {**vals, **{f"{k}_corr": v for k, v in vals.items()}}


def test_a_third_of_the_draws_with_the_day_above_the_whole_change(sr):
    A = [0.5] * 20 + [1.5] * 10  # A / L > 1 in a third of the draws: f_decay = 0 there, nothing dropped
    r = _draw_records(sr, A, [1.0] * 30)
    assert r["f_decay"]["nonfinite_draws"] == 0 and r["f_decay"]["p5"] == 0.0
    st = sr.statements(r, 30)
    assert st["entry_fraction"] == sr.ABOVE_BOTH  # 0.4 above f_decay's 0.394 and f_even's 0.25
    assert st["pace"] == sr.PACE_MORE and st["kind"] == sr.KIND_ALONG


def test_a_statement_with_undefined_draws_is_not_computable(sr):
    r = _draw_records(sr, [0.5] * 20 + [-0.2] * 10, [1.0] * 30)  # no decaying course in a third of the draws
    st = sr.statements(r, 30)
    assert st["entry_fraction"] == "not computable (10 of 30 draws undefined)"
    assert st["entry_fraction_raw"] == st["entry_fraction_corrected"] == st["entry_fraction"]
    assert st["pace"] == sr.NOT_SEPARATED  # A / L is defined in every draw
    r = _draw_records(sr, [0.5] * 30, [1.0] * 25 + [0.0] * 5)  # L = 0: A / L undefined
    assert sr.statements(r, 30)["pace"] == "not computable (5 of 30 draws undefined)"
    r = _draw_records(sr, [0.5] * 30, [1.0] * 30)
    r["K"] = {**r["K"], "nonfinite_draws": 2}
    assert sr.statements(r, 30)["kind"] == "not computable (2 of 30 draws undefined)"


def test_the_weights_are_library_3_cells_and_A_K_P_share_their_groups(sr, cfg):
    genes, fst = world()
    cg = {"CM": 1.0, "Endocardium": 2.0, "Endothelium": 1.0, "Epicardium": 3.0}
    n3 = {"CM": 40, "Endocardium": 21, "Endothelium": 30, "Epicardium": 10}  # Epicardium under min_cells_t1 there
    stages, _ = build_t1(sr, cfg, genes, fst, c_lib=(0.1, 0.2, 0.5), c_prev=0.25, c_group=cg, n_lib3=n3)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=2)
    v, rd = r["values"], r["reading"]
    assert np.isnan(v["b_E9.5_E10.5_lib3/Epicardium/share"]["value"])  # 10 cells: no value
    # the comparison-level combined values: weights from library 3, every group with a value
    assert v["b_E9.5_E10.5_lib1/combined_unweighted/share"]["value"] == pytest.approx(0.1 * 7 / 4, abs=1e-4)
    assert v["b_E9.5_E10.5_lib1/combined_weighted/share"]["value"] == pytest.approx(
        0.1 * (40 + 42 + 30 + 30) / 101, abs=1e-4)
    assert r["combined_groups"]["b_E9.5_E10.5_lib3/combined_weighted/share"] == ["CM", "Endocardium", "Endothelium"]
    # A and K over the common support of libraries 1, 2 and 3
    assert r["support"]["A_K_support"] == ["CM", "Endocardium", "Endothelium"]
    assert r["support"]["A_K_excluded"] == ["Epicardium"]
    wc = (40 * 1 + 21 * 2 + 30 * 1) / 91
    assert rd["A"]["value"] == pytest.approx(0.15 * wc, abs=1e-4)
    assert rd["K"]["value"] == pytest.approx(0.35 * wc, abs=1e-4)
    # P with and without the groups that stand on another group's cells at E8.5
    p = {g: v[f"a_E8.5_E9.5/{g}/share"]["value"] for g in ("CM", "Endocardium", "Endothelium")}
    assert p["CM"] == pytest.approx(0.25, abs=1e-4) and p["Endothelium"] == pytest.approx(0.0, abs=1e-6)
    assert r["support"]["P_fallback_groups"] == ["Endocardium"]
    assert rd["P"]["value"] == pytest.approx((40 * p["CM"] + 21 * p["Endocardium"] + 30 * p["Endothelium"]) / 91,
                                             abs=1e-4)
    assert rd["P_own_cells"]["value"] == pytest.approx((40 * p["CM"] + 30 * p["Endothelium"]) / 70, abs=1e-4)


def test_the_brief_fixed_choices(sr):
    assert sr.BOOT_SEED == (47, 0) and sr.N_BOOT == 200 and sr.PERCENTILES == (5.0, 95.0)
    assert sr.CAL_STREAM == (47, 1) and sr.CAL_FRACTION == pytest.approx(sr.feng_dirs.bracket_fraction(9.5, 10.5))
    assert sr.KIND_AT == 0.1 and (sr.PACE_FIFTH, sr.PACE_UPPER) == (0.2, 0.3) and sr.ENTRY_FRACTION == 0.4
    assert sr.pace_statement(0.2, 0.25) == sr.PACE_FIFTH_OR_LESS  # "above 0.2" is strict


def test_the_n48_scale_is_read_without_the_held_genes(sr, tmp_path):
    p = tmp_path / "n48.json"
    p.write_text(json.dumps({"groups": {
        "CM": {"js": {"scale": 0.7}, "js_without_held": {"scale": 0.6}},
        "Blood": {"js": {"scale": None}, "js_without_held": {"scale": None}}}}))
    assert sr.n48_scales(p) == {"CM": 0.6}
    assert sr.n48_scales(p, field="js") == {"CM": 0.7}


def test_the_share_reports_concentration_beside_it(sr):
    d = np.array([4.0, 1.0, 1.0, 1.0, 0.0])
    y = np.array([8.0, 1.0, 1.0, 1.0, 5.0])
    st = sr.share_stats(y, d, np.ones(5, bool), None)
    assert st["share"] == pytest.approx(35 / 19)
    assert st["n_eff"] == pytest.approx(19 ** 2 / 259)
    assert st["share_without_top_d2_half"] == pytest.approx(1.0)  # gene 0 alone carries 16 of 19
    assert st["genes_d_nonzero"] == 4 and st["pearson_d_nonzero"] == pytest.approx(1.0)
    assert st["pearson"] < 0.9  # the d = 0 gene pulls the all-gene correlation


def test_the_calibration_arm_reads_m1s_own_step(sr):
    rng = np.random.default_rng(5)
    G, n = 300, 400
    genes = [f"g{i}" for i in range(G)]
    C = rng.poisson(rng.uniform(5, 60, G), size=(n, G)).astype(np.int64)
    C[:, 0] = 1  # a count of one in every cell: counts are recoverable from .X
    libs = ["3"] * 300 + ["1"] * 100
    a = ad.AnnData(X=sr.t1_scope.log1p_cp10k(sparse.csr_matrix(C)),
                   obs=pd.DataFrame({"celltype": ["V-CM"] * n}, index=[f"c{i}_{li}" for i, li in enumerate(libs)]),
                   var=pd.DataFrame(index=genes))
    d = rng.normal(0, 0.5, G)
    d[0] = np.nan
    dirs, masks = {"CM": d}, {"CM": np.isfinite(d)}
    cal = sr.calibration(a, genes, dirs, masks, {"CM": ["V-CM"]}, [], {}, {"CM": 300.0}, EPS, 20)
    # a count-space step of 0.2 reads as a share within 0.02 of 0.2 (the stated tolerance)
    assert cal["groups"]["CM"]["share"] == pytest.approx(0.2, abs=0.02)
    assert cal["combined_weighted"] == cal["groups"]["CM"]["share"]
    again = sr.calibration(a, genes, dirs, masks, {"CM": ["V-CM"]}, [], {}, {"CM": 300.0}, EPS, 20)
    assert again["combined_weighted"] == cal["combined_weighted"]  # seeded


def test_the_run_reads_n48_without_held_and_runs_the_calibration(sr, monkeypatch):
    monkeypatch.setattr(sr.n12, "tries_rows", lambda *a: {"T47"})
    monkeypatch.setattr(sr.gv, "tree_state", lambda: {"commit": "abc", "clean": True})
    monkeypatch.setattr(sr.hygiene, "load_union", lambda *a: {"file": "x", "genes": []})
    monkeypatch.setattr(sr.t1_scope, "held_flat_genes", lambda *a: [])
    monkeypatch.setattr(sr, "feng_inputs", lambda *a: ({}, {}))
    monkeypatch.setattr(sr, "n48_scales", lambda path=None, field="js_without_held": {"CM": field})
    monkeypatch.setattr(sr.ad, "read_h5ad", lambda p: None)
    seen = {}

    def spy(stages, feng, held, scales, cfg, **kw):
        seen.update(scales=scales, **kw)
        raise RuntimeError("stop")

    monkeypatch.setattr(sr, "reading", spy)
    with pytest.raises(RuntimeError, match="stop"):
        sr.run("x")
    assert seen == {"scales": {"CM": "js_without_held"}, "calibrate": True, "self_check": True}


# ---------------------------------------------------------------------------------------------------------------------
# Second review: the corrected fraction in units of the true d, r_d from the share's own genes, the fraction statement
# with undefined draws, one support for raw and corrected, boundaries, the pairs of L and L_corr, unchanged quantities


def _fake_measure(sr, shares: dict, corrected: dict, L=1.0, L_corr=1.0):
    """One measurement's flat values for reading_values: per group the same share in libraries 1, 2 and 3."""
    out = {"scale_E9.5/median_slope": L, "scale_E9.5/median_slope_corrected": L_corr,
           "scale_E9.5/median_slope_all_pairs": L}
    for g in shares:
        for name in sr.AK_COMPARISONS:
            out[f"{name}/{g}/share"] = shares[g]
            out[f"{name}/{g}/share_corrected"] = corrected[g]
        out[f"a_E8.5_E9.5/{g}/share"] = 0.1
    return out


def _constant_records(sr, rv: dict, n=10):
    return sr.summarise(rv, {k: np.full(n, v) for k, v in rv.items()}, n)


def test_the_corrected_fraction_reads_the_courses_in_units_of_the_true_direction(sr):
    # true A / L = 0.2 and L = 1, Feng's d noisy with r_d = 0.6: the raw share is 0.12, the corrected one 0.2
    setup = {"dirs": {"CM": None}, "weights": {"CM": 1.0}, "scales": {"CM": 0.6},
             "source": {"E8.5": {"CM": "own"}, "E9.5": {"CM": "own"}}}
    rv = sr.reading_values(_fake_measure(sr, {"CM": 0.12}, {"CM": 0.2}), setup)
    assert rv["A_corr"] == pytest.approx(0.2) and rv["r_d_weighted"] == pytest.approx(0.6)
    assert rv["f_even_corr"] == pytest.approx(0.4) and rv["f_decay_corr"] == pytest.approx(0.4)
    assert rv["f_max_corr_noisy_d"] == pytest.approx(0.48) and rv["f_even_corr_noisy_d"] == pytest.approx(0.24)
    st = sr.statements(_constant_records(sr, rv), 10)
    assert st["entry_fraction_corrected"] == sr.INSIDE  # the r_d-multiplied courses (0.24) would read "above both"
    assert st["pace_corrected"] == sr.PACE_FIFTH_OR_LESS


def test_r_d_is_computed_from_the_bracket_halves_over_the_shares_own_genes(sr, cfg):
    genes, fst = world(half_noise=0.1)
    stages, _ = build_t1(sr, cfg, genes, fst)
    held = [genes[3], genes[9], genes[20]]
    r = sr.reading(stages, fst, held, {"CM": 0.5, "Epicardium": 0.9}, cfg, n_boot=1)
    dirs, _, idx = sr.feng_dirs.directions_from_stats(fst, genes, cfg, sr.BRACKET)
    not_held = np.array([g not in held for g in genes])
    for g in r["groups"]:
        fl = cfg["external_labels"][g][0]

        def half(h, fl=fl):
            a, b = sr.dl.feng_sel(fst, "E9.5", fl, (h,)), sr.dl.feng_sel(fst, "E14.5", fl, (h,))
            return sr.dl.on_t1(np.log((b["cp"] + EPS) / (a["cp"] + EPS)), idx)

        d, d_a, d_b = dirs[g], half("A"), half("B")
        ok = np.isfinite(d) & not_held & np.isfinite(d_a) & np.isfinite(d_b)
        msq = np.mean(d[ok] ** 2)
        hand = max(msq - np.mean(((d_a[ok] - d_b[ok]) / 2) ** 2), 0.0) / msq
        assert r["r_d"][g] == pytest.approx(hand, rel=1e-12)
        every = np.isfinite(d) & np.isfinite(d_a) & np.isfinite(d_b)  # the held genes in: a different value
        msq_all = np.mean(d[every] ** 2)
        assert r["r_d"][g] != max(msq_all - np.mean(((d_a[every] - d_b[every]) / 2) ** 2), 0) / msq_all
        v = r["values"][f"b_E9.5_E10.5_lib1/{g}/share"]["value"]
        assert r["values"][f"b_E9.5_E10.5_lib1/{g}/share_corrected"]["value"] == pytest.approx(v / r["r_d"][g])
    assert r["n48_scale"] == {"CM": 0.5, "Endocardium": None, "Endothelium": None, "Epicardium": 0.9}  # beside only


def test_above_the_bound_is_read_from_f_max_even_where_f_decay_is_undefined(sr):
    # f_max's upper percentile 0.27 < 0.4 in every draw; A < 0 in 5 of 30 draws leaves f_decay undefined there
    r = _draw_records(sr, [0.05] * 25 + [-0.02] * 5, [0.25] * 30, A0=0.05, L0=0.25)
    assert r["f_decay"]["nonfinite_draws"] == 5
    st = sr.statements(r, 30)
    assert st["entry_fraction_raw"] == st["entry_fraction"] == sr.ABOVE_BOUND
    nan = float("nan")
    assert sr.fraction_statement(_rec(0.2, 0.18, 0.27), _rec(0.1, 0.09, 0.13), _rec(nan, nan, nan)) == sr.ABOVE_BOUND
    # f_max itself undefined in a draw: not computable
    r = _draw_records(sr, [0.05] * 30, [0.25] * 28 + [nan] * 2, A0=0.05, L0=0.25)
    assert sr.statements(r, 30)["entry_fraction"] == "not computable (2 of 30 draws undefined)"


def test_raw_and_corrected_A_and_K_share_one_set_of_groups(sr):
    nan = float("nan")
    setup = {"dirs": {"CM": None, "Endothelium": None, "Epicardium": None},
             "weights": {"CM": 1.0, "Endothelium": 3.0, "Epicardium": 1.0},
             "scales": {"CM": 0.5, "Endothelium": nan, "Epicardium": 0.8}, "source": {"E8.5": {}, "E9.5": {}}}
    out = _fake_measure(sr, {"CM": 0.1, "Endothelium": 0.4, "Epicardium": nan},
                        {"CM": 0.2, "Endothelium": nan, "Epicardium": nan})
    info = {}
    rv = sr.reading_values(out, setup, info)
    assert rv["A"] == pytest.approx(0.1) and rv["A_corr"] == pytest.approx(0.2)  # both over CM alone
    assert rv["r_d_weighted"] == pytest.approx(0.5)
    assert info["A_K_support"] == ["CM"] and info["A_K_excluded"] == ["Endothelium", "Epicardium"]
    assert info["A_K_excluded_no_share"] == ["Epicardium"]
    assert info["A_K_excluded_no_share_corrected"] == ["Endothelium"]


def test_the_threshold_boundaries(sr):
    assert sr.pace_statement(0.15, 0.3) == sr.NOT_SEPARATED  # "below 0.3" is strict
    assert sr.kind_statement(0.1, 0.4) == sr.NOT_SEPARATED  # "above 0.1" is strict
    assert sr.kind_statement(-0.1, 0.1) == sr.NOT_SEPARATED  # "below 0.1" is strict


def test_r_d_weighted_is_the_library_3_weighted_mean(sr, cfg):
    genes, fst = world(half_noise=0.15, seed=3)
    n3 = {"CM": 60, "Endocardium": 21, "Endothelium": 35, "Epicardium": 25}
    stages, _ = build_t1(sr, cfg, genes, fst, n_lib3=n3)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=1)
    s = r["r_d"]
    assert len({round(v, 6) for v in s.values()}) == 4  # unequal scales
    assert r["support"]["A_K_support"] == ["CM", "Endocardium", "Endothelium", "Epicardium"]
    wtd = sum(n3[g] * s[g] for g in n3) / sum(n3.values())
    assert r["reading"]["r_d_weighted"]["value"] == pytest.approx(wtd)
    assert r["reading"]["r_d_weighted"]["value"] != pytest.approx(np.mean(list(s.values())))


def test_the_pairs_of_L_and_L_corr_are_written_and_fewer_make_the_corrected_statements_not_computable(
        sr, cfg, monkeypatch):
    genes, fst = world(half_noise=0.05)
    stages, _ = build_t1(sr, cfg, genes, fst, c_lib=(0.18, 0.18, 0.6), k=1.0)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=2)
    sp = r["scale_pairs"]["E9.5"]
    kept = [p for p in r["genes_per_pair"] if p != "Endocardium~Endothelium"]
    assert sp["L"] == sp["L_corr"] == kept and sp["L_all_pairs"] == list(r["genes_per_pair"])
    assert not r["statements"]["pace_corrected"].startswith(sr.NOT_COMPUTABLE)
    real = sr.feng_pair_contrasts

    def floored(*a, **k):  # one pair's Feng contrast all noise: r_c = 0
        contrast, rel = real(*a, **k)
        rel["E9.5"][("CM", "Epicardium")] = 0.0
        return contrast, rel

    monkeypatch.setattr(sr, "feng_pair_contrasts", floored)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=2)
    sp = r["scale_pairs"]["E9.5"]
    assert sp["L_corr"] == [p for p in kept if p != "CM~Epicardium"] and sp["L"] == kept
    msg = f"not computable (L_corr rests on {len(kept) - 1} of the {len(kept)} pairs of L)"
    st = r["statements"]
    assert st["pace_corrected"] == st["entry_fraction_corrected"] == st["pace"] == st["entry_fraction"] == msg
    assert not st["pace_raw"].startswith(sr.NOT_COMPUTABLE)
    txt = sr.text({"label": sr.LABEL_T47, "read": r, "stamp": {"commit": "abc", "tree_clean": True, "utc": "t"}})
    assert f"L_corr rests on {len(kept) - 1} pairs: " in txt and f"L rests on {len(kept)} pairs: " in txt


def test_why_a_joined_statement_is_not_separated(sr):
    assert sr.not_separated_cause(sr.INSIDE, sr.INSIDE) is None
    assert sr.not_separated_cause(sr.PACE_MORE, sr.PACE_FIFTH_OR_LESS) == sr.CORRECTION_CHANGES
    assert sr.not_separated_cause(sr.PACE_MORE, sr.NOT_SEPARATED) == "the interval does not separate (corrected)"
    assert sr.not_separated_cause(sr.NOT_SEPARATED, sr.NOT_SEPARATED) == \
        "the interval does not separate (raw and corrected)"
    assert sr.not_separated_cause("not computable (1 of 2 draws undefined)", sr.INSIDE) is None


def test_the_ledger_sentences_are_in_the_docstring_and_printed(sr, cfg):
    flat = " ".join(sr.__doc__.split())
    assert len(sr.LEDGER_NOTES) == 21 and all(n in flat for n in sr.LEDGER_NOTES)
    for start in ("The corrected fraction statement compares 0.4 with the courses in units of Feng's true direction",
                  "K counts any difference between the heart library", "L is the size of between-type contrasts",
                  "A share is what the genes with the largest Feng change say", "Not separated may mean",
                  "The courses assume the same E9.5 and E14.5", "P over own cells is the reading",
                  "No statement of T47 is a reason to change M1's fraction", "The intervals hold only"):
        assert any(n.startswith(start) for n in sr.LEDGER_NOTES)
    genes, fst = world(half_noise=0.1)
    stages, _ = build_t1(sr, cfg, genes, fst)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=2)
    txt = sr.text({"label": sr.LABEL_T47, "read": r, "stamp": {"commit": "abc", "tree_clean": True, "utc": "t"}})
    assert all(f"Note: {n}" in txt for n in sr.LEDGER_NOTES)
    assert "corrected (units of Feng's true direction)" in txt and "reported only, the corrected courses" in txt


# Values at the base commit 82bdece24f on this synthetic world, every quantity no finding changes (sha256 of their
# json): the per-comparison shares and fields except share_corrected, the scales and pair slopes, A, K, P, L and the
# raw courses; and the raw statements.
GOLDEN_SHA = "e0acf3898ab13e5efc99c122df6b391ef6c9f5de0674296185f73801d0436e10"
GOLDEN_N = 292
GOLDEN_TABLE = REPO_ROOT / "tests" / "golden" / "t1_share_read_unchanged.json"
CHANGED = ("share_corrected", "reading/A_corr", "reading/K_corr", "reading/r_d_weighted", "reading/A_over_L_corr",
           "reading/x_corr", "reading/f_max_corr", "reading/f_even_corr", "reading/f_decay_corr")


def test_every_unchanged_quantity_keeps_its_output_bit_for_bit(sr, cfg):
    genes, fst = world(seed=11, half_noise=0.08)
    cg = {"CM": 1.0, "Endocardium": 1.5, "Endothelium": 0.7, "Epicardium": 1.2}
    n3 = {"CM": 40, "Endocardium": 25, "Endothelium": 30, "Epicardium": 22}
    stages, _ = build_t1(sr, cfg, genes, fst, c_lib=(0.15, 0.25, 0.6), c_prev=0.2, k=1.3, c_group=cg, n_lib3=n3)
    rng = np.random.default_rng(12)
    for a in stages.values():
        a.X = sparse.csr_matrix(np.log1p(np.expm1(a.X.toarray()) * rng.uniform(0.5, 1.5, a.shape)).astype(np.float32))
    r = sr.reading(stages, fst, [genes[7]], {"CM": 0.5, "Endothelium": 0.8}, cfg, n_boot=3)
    kept = {k: v for k, v in r["values"].items() if not any(c in k for c in CHANGED)}
    assert len(kept) == GOLDEN_N
    # The exact hash (GOLDEN_SHA) holds on the machine the table was made on. The stages' matrices are float32, and
    # another machine's expm1 and log1p differ from this one's at about 1e-8 relative (found on 2026-11-17, when the
    # full suite ran on the DGX). So the values are compared with the stored table to five digits, and the hash is
    # kept only as the table's own check.
    golden = json.loads(GOLDEN_TABLE.read_text())
    assert hashlib.sha256(json.dumps(golden, sort_keys=True).encode()).hexdigest() == GOLDEN_SHA
    assert sorted(kept) == sorted(golden)
    for key, rec in kept.items():
        assert sorted(rec) == sorted(golden[key]), key
        for field, value in rec.items():
            want = golden[key][field]
            if isinstance(want, float) and math.isnan(want):
                assert isinstance(value, float) and math.isnan(value), (key, field)
            else:
                assert value == pytest.approx(want, rel=1e-5, abs=1e-6), (key, field)
    assert {s: r["statements"][s] for s in ("kind", "pace_raw", "entry_fraction_raw")} == {
        "kind": sr.KIND_ALONG, "pace_raw": sr.PACE_FIFTH_OR_LESS, "entry_fraction_raw": sr.BELOW_BOTH}


# ---------------------------------------------------------------------------------------------------------------------
# Third review: the self-check arm and the widened statements, the calibration's r_d, undefined at the point, K's
# support, own cells at E9.5, two labels in the calibration, an empty Feng half, the new ledger sentences


def count_world(sr, cfg, n=80, depth=60.0, seed=21, c_lib=(0.2, 0.3, 0.9), half_noise=0.0):
    """Three T1 stages of Poisson counts (X = log1p CP10k of them, counts recoverable), n cells per label and library:
    E9.5 at Feng's E9.5 means (libraries 1 and 2), E10.5 libraries 1, 2 and 3 moved by c_lib times each group's
    direction, E8.5 with V-CM and BEC (library 1)."""
    genes, fst = world(seed=seed, half_noise=half_noise)
    dirs, _, idx = sr.feng_dirs.directions_from_stats(fst, genes, cfg, sr.BRACKET)
    rng = np.random.default_rng(seed + 1)
    e95 = {}
    for g in T1_LABEL:
        v = sr.dl.on_t1(sr.dl.feng_sel(fst, "E9.5", cfg["external_labels"][g][0])["cp"], idx)
        v[~np.isfinite(v)] = 5.0
        e95[g] = v
    d0 = {g: np.nan_to_num(dirs[g]) for g in dirs}

    def build(blocks):
        rows, labs, names = [], [], []
        for (lab, lib), mean in blocks.items():
            C = rng.poisson(np.maximum(mean, 0.0) * depth, size=(n, len(genes)))
            rows.append(C)
            labs += [lab] * n
            names += [f"c{len(names) + i}_{lib}" for i in range(n)]
        C = np.vstack(rows).astype(np.int64)
        C[:, 0] = 1  # g0 (no Feng value): a count of one in every cell, so counts are recoverable from .X
        C = sparse.csr_matrix(C)
        return ad.AnnData(X=sr.t1_scope.log1p_cp10k(C), obs=pd.DataFrame({"celltype": labs}, index=names),
                          var=pd.DataFrame(index=list(genes)))

    E95 = {(T1_LABEL[g], li): e95[g] for g in T1_LABEL for li in ("1", "2")}
    E105 = {(T1_LABEL[g], li): (e95[g] + EPS) * np.exp(c * d0[g]) - EPS for g in T1_LABEL
            for li, c in zip(("1", "2", "3"), c_lib, strict=True)}
    E85 = {("V-CM", "1"): e95["CM"], ("BEC", "1"): e95["Endothelium"]}
    return genes, fst, {"E8.5": build(E85), "E9.5": build(E95), "E10.5": build(E105)}


def test_the_self_check_split_is_seeded_and_splits_every_block_into_quarters_quarters_and_half(sr):
    a = _noisy_stage()
    parts = sr.split_cells(a)
    again = sr.split_cells(a)
    assert all(np.array_equal(parts[p], again[p]) for p in parts)
    other = sr.split_cells(a, seed=(47, 4))
    assert not all(np.array_equal(parts[p], other[p]) for p in parts)
    every = np.concatenate([parts["H1a"], parts["H1b"], parts["H2"]])
    assert sorted(every) == list(range(a.n_obs))  # a partition of the cells
    labs = a.obs["celltype"].to_numpy()
    libs = np.array([n.split("_")[-1] for n in a.obs_names])
    for lab, lib in (("A", "1"), ("A", "2"), ("B", "1")):  # blocks of 7, 5 and 6 cells
        m = (labs == lab) & (libs == lib)
        n = int(m.sum())
        cnt = {p: int(m[parts[p]].sum()) for p in parts}
        assert cnt == {"H1a": n // 4, "H1b": n // 2 - n // 4, "H2": n - n // 2}
    assert sr.SPLIT_STREAM == (47, 3) and sr.SELF_STREAM == (47, 2)


def test_the_self_check_world_and_its_bias(sr, cfg, monkeypatch):
    genes, fst, stages = count_world(sr, cfg)
    seen = {}
    real = sr.t1_scope.external_shift

    def spy(adata, dirs, fraction, group_of, held, **kw):
        seen.setdefault("calls", []).append({"n": adata.n_obs, "fraction": fraction, "seed": kw.get("seed"),
                                             "held": list(held), "group_of": dict(group_of)})
        return real(adata, dirs, fraction, group_of, held, **kw)

    monkeypatch.setattr(sr.t1_scope, "external_shift", spy)
    held = [genes[4]]
    r = sr.reading(stages, fst, held, {}, cfg, n_boot=2, self_check=True)
    sc = r["self_check"]
    parts = sr.split_cells(stages["E9.5"])
    (call,) = seen["calls"]  # one shift: H1's cells, M1's 0.2 step on stream (47, 2), the held genes flat
    assert call["n"] == len(parts["H1a"]) + len(parts["H1b"]) and call["fraction"] == 0.2
    assert call["seed"] == [47, 2] and call["held"] == held
    assert call["group_of"] == {"V-CM": "CM", "Endocardium": "Endocardium", "BEC": "Endothelium",
                                "Epicardium": "Epicardium"}
    known = {"A": 0.2, "L": 1.0, "A_over_L": 0.2, "f_max": 0.8, "f_even": 0.4, "f_decay": 0.4}
    assert sc["known"] == {**known, **{f"{k}_corr": v for k, v in known.items()}}
    for k, v in sc["bias"].items():
        assert v == sc["measured"][k] - sc["known"][k]
    # high counts: the world recovers its known answer closely, and r_d = 1 makes A_corr the raw A
    assert sc["measured"]["A"] == pytest.approx(0.2, abs=0.02)
    assert sc["measured"]["A_corr"] == sc["measured"]["A"]
    assert sc["measured"]["L_corr"] == pytest.approx(1.0, abs=0.05)
    assert sc["support"] == r["support"]["A_K_support"]
    # measured by hand: the courses from the measured (A, L)
    c = sr.courses(sc["measured"]["A"], sc["measured"]["L"])
    assert sc["measured"]["f_even"] == c["f_even"] and sc["measured"]["f_decay"] == c["f_decay"]
    for p, rec in sc["pairs"].items():
        if p != "Endocardium~Endothelium":
            assert 0 < rec["r_c"] <= 1 and rec["slope_corrected"] == pytest.approx(rec["slope"] / rec["r_c"])
    for s in ("pace_raw", "pace_corrected", "entry_fraction_raw", "entry_fraction_corrected"):
        assert r["statements"][f"{s}_decided_by"] in (sr.BOOT_WIDTH, sr.BIAS_WIDTH)
    txt = sr.text({"label": sr.LABEL_T47, "read": r, "stamp": {"commit": "abc", "tree_clean": True, "utc": "t"}})
    assert "Self-check (a known world from E9.5's cells" in txt and "[decided by the " in txt


def test_the_self_check_scale_by_hand_on_its_own_world(sr, cfg):
    """One pair's self-check slope by hand: the challenge's E9.5 contrast from H2's cells on the stand-in Feng contrast
    from H1's cells, both g against h, over the genes both groups keep."""
    genes, fst, stages = count_world(sr, cfg, seed=31)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=1, self_check=True)
    sc = r["self_check"]
    a95 = stages["E9.5"]
    parts = sr.split_cells(a95)
    h1 = np.sort(np.concatenate([parts["H1a"], parts["H1b"]]))
    lab = a95.obs["celltype"].to_numpy()

    def mean(idx, g):
        m = idx[lab[idx] == T1_LABEL[g]]
        return np.expm1(a95.X[m].toarray().astype(np.float32)).astype(np.float64).mean(0)

    def contrast(idx):
        return np.log((mean(idx, "CM") + EPS) / (mean(idx, "Epicardium") + EPS))

    dirs = sr.feng_dirs.directions_from_stats(fst, genes, cfg, sr.BRACKET)[0]
    m = np.isfinite(dirs["CM"]) & np.isfinite(dirs["Epicardium"])
    x, y = contrast(h1)[m], contrast(parts["H2"])[m]
    assert sc["pairs"]["CM~Epicardium"]["slope"] == pytest.approx(float(y @ x / (x @ x)), rel=1e-5)
    xa, xb = contrast(parts["H1a"])[m], contrast(parts["H1b"])[m]
    rc = (np.mean(x ** 2) - np.mean(((xa - xb) / 2) ** 2)) / np.mean(x ** 2)
    assert sc["pairs"]["CM~Epicardium"]["r_c"] == pytest.approx(rc, rel=1e-4)


def _bias(value=0.0, **kw):
    keys = [f"{k}{t}" for t in ("", "_corr") for k in ("A", "L", "A_over_L", "f_max", "f_even", "f_decay")]
    return {**dict.fromkeys(keys, value), **kw}


def test_the_fraction_statement_is_not_separated_where_0_4_lies_inside_the_bias_widening(sr):
    # the reviewer's case: true 0.4, both courses read about 0.38 with intervals of +-0.0015: "above both courses" on
    # the bootstrap interval; a self-check bias of -0.03 on the courses widens them past 0.4
    rv = {"K": 0.3, **sr.courses(0.19, 0.95)}  # f_even 0.38, f_decay 0.38
    recs = {k: {"value": v, "p5": v - 0.0015, "p95": v + 0.0015, "nonfinite_draws": 0} for k, v in rv.items()}
    r = {**recs, **{f"{k}_corr": v for k, v in recs.items()}}
    plain = sr.statements(r, 10)
    assert plain["entry_fraction"] == sr.ABOVE_BOTH and "entry_fraction_raw_decided_by" not in plain
    st = sr.statements(r, 10, bias=_bias(0.0, f_even=-0.03, f_decay=-0.03, f_even_corr=-0.03, f_decay_corr=-0.03))
    assert st["entry_fraction_raw"] == st["entry_fraction_corrected"] == st["entry_fraction"] == sr.NOT_SEPARATED
    assert st["entry_fraction_raw_decided_by"] == st["entry_fraction_corrected_decided_by"] == sr.BIAS_WIDTH
    assert st["entry_fraction_why_not_separated"] == ("the interval does not separate (raw widened by the self-check's "
                                                      "bias and corrected widened by the self-check's bias)")
    # a bias smaller than the margin leaves the statement, decided by the bootstrap interval
    st = sr.statements(r, 10, bias=_bias(0.005))
    assert st["entry_fraction"] == sr.ABOVE_BOTH and st["entry_fraction_raw_decided_by"] == sr.BOOT_WIDTH
    # the sign of the bias does not matter: each side widens by its absolute value
    assert sr.statements(r, 10, bias=_bias(0.0, f_even=0.03, f_decay=0.03, f_even_corr=0.03, f_decay_corr=0.03))[
        "entry_fraction"] == sr.NOT_SEPARATED
    # "inside" on the bootstrap interval stays "inside"; "above the bound" widened by f_max's own bias
    assert sr.statements(r, 10, bias=_bias(0.0, f_max=0.05))["entry_fraction"] == sr.ABOVE_BOTH
    rv2 = {"K": 0.3, **sr.courses(0.2, 1.0)}
    r2 = {k: {"value": v, "p5": v - 0.01, "p95": v + 0.01, "nonfinite_draws": 0} for k, v in rv2.items()}
    r2 = {**r2, **{f"{k}_corr": v for k, v in r2.items()}}
    assert sr.statements(r2, 10, bias=_bias(0.2))["entry_fraction"] == sr.INSIDE
    rv3 = {"K": 0.3, **sr.courses(0.05, 0.42)}  # f_max 0.37 < 0.4: above the bound
    r3 = {k: {"value": v, "p5": v - 0.005, "p95": v + 0.005, "nonfinite_draws": 0} for k, v in rv3.items()}
    r3 = {**r3, **{f"{k}_corr": v for k, v in r3.items()}}
    assert sr.statements(r3, 10, bias=_bias(0.0))["entry_fraction"] == sr.ABOVE_BOUND
    st = sr.statements(r3, 10, bias=_bias(0.0, f_max=-0.03, f_max_corr=0.03))
    assert st["entry_fraction"] == sr.NOT_SEPARATED and st["entry_fraction_raw_decided_by"] == sr.BIAS_WIDTH
    # an undefined bias makes the statement not computable
    st = sr.statements(r3, 10, bias=_bias(0.0, f_decay_corr=float("nan")))
    assert st["entry_fraction_corrected"] == "not computable (the self-check's bias of f_decay_corr is undefined)"
    assert st["entry_fraction_corrected_decided_by"] is None and st["entry_fraction"] == st["entry_fraction_corrected"]


def test_the_pace_statement_widens_by_the_bias_of_a_over_l(sr):
    r = _draw_records(sr, [0.21] * 10, [1.0] * 10, A0=0.21)  # A / L = 0.21 in every draw: "more than a fifth"
    assert sr.statements(r, 10)["pace"] == sr.PACE_MORE
    st = sr.statements(r, 10, bias=_bias(0.0, A_over_L=0.02))  # 0.2 inside the widened raw interval
    assert st["pace_raw"] == sr.NOT_SEPARATED and st["pace_raw_decided_by"] == sr.BIAS_WIDTH
    assert st["pace_corrected"] == sr.PACE_MORE and st["pace_corrected_decided_by"] == sr.BOOT_WIDTH
    assert st["pace"] == sr.NOT_SEPARATED
    assert st["pace_why_not_separated"] == "the interval does not separate (raw widened by the self-check's bias)"
    r = _draw_records(sr, [0.18] * 10, [1.0] * 10, A0=0.18)  # "about a fifth or less" until 0.3 enters the widening
    assert sr.statements(r, 10, bias=_bias(0.1))["pace"] == sr.PACE_FIFTH_OR_LESS
    assert sr.statements(r, 10, bias=_bias(0.13))["pace"] == sr.NOT_SEPARATED
    # the bias of the quantity the statement reads, and no other
    assert sr.statements(r, 10, bias=_bias(0.0, A=0.5, L=0.5, A_corr=0.5, L_corr=0.5))["pace"] == \
        sr.PACE_FIFTH_OR_LESS


def test_the_self_check_leaves_every_existing_value_bit_for_bit(sr, cfg):
    genes, fst, stages = count_world(sr, cfg, seed=41, half_noise=0.05)
    off = sr.reading(stages, fst, [genes[7]], {}, cfg, n_boot=3, calibrate=True)
    on = sr.reading(stages, fst, [genes[7]], {}, cfg, n_boot=3, calibrate=True, self_check=True)
    assert json.dumps(off["values"], sort_keys=True) == json.dumps(on["values"], sort_keys=True)
    for k in ("reading", "r_d", "support", "scale_pairs", "feng_reliability", "combined_groups"):
        assert json.dumps(off[k], sort_keys=True, default=str) == json.dumps(on[k], sort_keys=True, default=str)
    assert json.dumps(off["calibration"], sort_keys=True) == json.dumps(on["calibration"], sort_keys=True)
    assert on["statements"]["kind"] == off["statements"]["kind"]
    assert "self_check" in on and "self_check" not in off


def test_no_feng_stage_other_than_the_bracket_ends_is_read_with_the_self_check(sr, cfg, monkeypatch):
    genes, fst, stages = count_world(sr, cfg)
    seen, state = [], {"in_arm": False}
    real_sel, real_arm = sr.dl.feng_sel, sr.self_check_arm

    def spy(st, stage, labels, halves=("A", "B")):
        seen.append((stage, state["in_arm"]))
        return real_sel(st, stage, labels, halves)

    def arm(*a, **k):
        state["in_arm"] = True
        try:
            return real_arm(*a, **k)
        finally:
            state["in_arm"] = False

    monkeypatch.setattr(sr.dl, "feng_sel", spy)
    monkeypatch.setattr(sr, "self_check_arm", arm)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=1, self_check=True)
    assert "self_check" in r and seen
    assert {s for s, _ in seen} <= set(sr.BRACKET)
    assert not any(inside for _, inside in seen)  # the arm reads no Feng statistics: only d, from the setup


def test_the_calibration_divides_by_the_reads_r_d(sr, cfg):
    genes, fst, stages = count_world(sr, cfg, half_noise=0.1, seed=51)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=1, calibrate=True)
    cal = r["calibration"]
    for g, rec in cal["groups"].items():
        assert 0 < r["r_d"][g] < 1
        assert rec["share_corrected"] == pytest.approx(rec["share"] / r["r_d"][g], rel=1e-12)
    s_cal = cal["combined_weighted"]
    for k in ("f_max", "f_even", "f_decay", "f_max_corr", "f_even_corr", "f_decay_corr"):
        assert cal[f"{k}_in_M1_fraction"] == pytest.approx(r["reading"][k]["value"] * 0.2 / s_cal)


def test_the_calibration_share_corrected_uses_the_scales_handed_over(sr):
    rng = np.random.default_rng(5)
    G, n = 200, 120
    genes = [f"g{i}" for i in range(G)]
    C = rng.poisson(rng.uniform(5, 60, G), size=(n, G)).astype(np.int64)
    C[:, 0] = 1  # counts recoverable from .X
    a = ad.AnnData(X=sr.t1_scope.log1p_cp10k(sparse.csr_matrix(C)),
                   obs=pd.DataFrame({"celltype": ["V-CM"] * n}, index=[f"c{i}_3" for i in range(n)]),
                   var=pd.DataFrame(index=genes))
    d = rng.normal(0, 0.5, G)
    cal = sr.calibration(a, genes, {"CM": d}, {"CM": np.ones(G, bool)}, {"CM": ["V-CM"]}, [], {"CM": 0.5},
                         {"CM": 1.0}, EPS, 20)
    assert cal["groups"]["CM"]["share_corrected"] == cal["groups"]["CM"]["share"] / 0.5


def test_the_calibration_takes_every_label_of_a_group(sr):
    rng = np.random.default_rng(6)
    G, n = 200, 120
    genes = [f"g{i}" for i in range(G)]
    C = np.vstack([rng.poisson(rng.uniform(5, 60, G), size=(n, G)), rng.poisson(rng.uniform(1, 9, G), size=(n, G))])
    C[:, 0] = 1  # counts recoverable from .X
    two = ["V-CM"] * n + ["A-CM"] * n
    names = [f"c{i}_3" for i in range(2 * n)]

    def adata(labels):
        return ad.AnnData(X=sr.t1_scope.log1p_cp10k(sparse.csr_matrix(C.astype(np.int64))),
                          obs=pd.DataFrame({"celltype": labels}, index=names), var=pd.DataFrame(index=genes))

    d = rng.normal(0, 0.5, G)
    args = ({"CM": d}, {"CM": np.ones(G, bool)})
    both = sr.calibration(adata(two), genes, *args, {"CM": ["V-CM", "A-CM"]}, [], {}, {"CM": 1.0}, EPS, 20)
    one_label = sr.calibration(adata(["V-CM"] * (2 * n)), genes, *args, {"CM": ["V-CM"]}, [], {}, {"CM": 1.0}, EPS,
                               20)
    assert both["groups"]["CM"]["share"] == one_label["groups"]["CM"]["share"]  # both labels' cells enter
    first = sr.calibration(adata(two), genes, *args, {"CM": ["V-CM"]}, [], {}, {"CM": 1.0}, EPS, 20)
    assert first["groups"]["CM"]["share"] != both["groups"]["CM"]["share"]


def test_a_statement_undefined_at_the_point_says_so(sr):
    nan = float("nan")
    k = {"value": nan, "p5": 0.2, "p95": 0.3, "nonfinite_draws": 0}
    assert sr.undefined([k], 200) == "not computable (undefined at the point; 0 of 200 draws undefined)"
    r = _draw_records(sr, [0.5] * 30, [1.0] * 30)
    r["K"] = k
    assert sr.statements(r, 30)["kind"] == "not computable (undefined at the point; 0 of 30 draws undefined)"
    assert sr.undefined([{**k, "value": 0.1, "nonfinite_draws": 3}], 200) == "not computable (3 of 200 draws undefined)"
    assert sr.undefined([{**k, "nonfinite_draws": 3}], 200) == \
        "not computable (undefined at the point; 3 of 200 draws undefined)"


def test_K_is_read_over_the_joint_support(sr):
    nan = float("nan")
    setup = {"dirs": {"CM": None, "Endothelium": None, "Epicardium": None},
             "weights": {"CM": 1.0, "Endothelium": 3.0, "Epicardium": 1.0},
             "scales": {"CM": 0.5, "Endothelium": nan, "Epicardium": 0.8}, "source": {"E8.5": {}, "E9.5": {}}}
    out = _fake_measure(sr, {"CM": 0.1, "Endothelium": 0.4, "Epicardium": nan},
                        {"CM": 0.2, "Endothelium": nan, "Epicardium": nan})
    out["b_E9.5_E10.5_lib3/CM/share"], out["b_E9.5_E10.5_lib3/CM/share_corrected"] = 0.7, 1.4
    out["b_E9.5_E10.5_lib3/Endothelium/share"] = 0.9  # outside the joint support: must not enter K
    rv = sr.reading_values(out, setup)
    assert rv["K"] == pytest.approx(0.7 - 0.1) and rv["K_corr"] == pytest.approx(1.4 - 0.2)


def test_a_group_that_falls_back_at_e9_5_leaves_p_over_own_cells(sr, cfg):
    genes, fst = world()
    stages, _ = build_t1(sr, cfg, genes, fst)
    e95 = stages["E9.5"]
    endo = e95[e95.obs["celltype"] == "Endocardium"].copy()
    endo.obs_names = [f"e{i}_1" for i in range(endo.n_obs)]
    stages["E8.5"] = ad.concat([stages["E8.5"], endo])  # Endocardium's own cells at E8.5
    stages["E9.5"] = e95[e95.obs["celltype"] != "Endocardium"].copy()  # and none at E9.5
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=1)
    assert r["label_source"]["E8.5"]["Endocardium"] == "own"
    assert r["label_source"]["E9.5"]["Endocardium"] == "fallback:Endothelium"
    assert "Endocardium" in r["support"]["P_support"]
    assert "Endocardium" not in r["support"]["P_own_cells_support"]
    assert "Endocardium" in r["support"]["P_fallback_groups"]


def test_an_empty_feng_half_gives_that_group_no_r_d_and_the_read_goes_on(sr, cfg):
    genes, fst = world(half_noise=0.1)
    keys = list(fst["keys"])
    i = keys.index("E14.5|z|epicardial|B")
    j = keys.index("E14.5|z|epicardial|A")
    for f in ("S_log", "S_cp", "S_nz"):  # half B's cells moved into half A: the pooled mean is unchanged
        fst[f][j] = fst[f][j] + fst[f][i]
        fst[f][i] = 0.0
    fst["n"][j], fst["n"][i] = fst["n"][j] + fst["n"][i], 0.0
    with pytest.raises(ValueError, match="a half of its Feng cells is empty"):
        sr.feng_dirs.direction_halves(fst, genes, cfg, sr.BRACKET)
    stages, _ = build_t1(sr, cfg, genes, fst)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=1)
    assert np.isnan(r["r_d"]["Epicardium"])
    assert all(0 < r["r_d"][g] <= 1 for g in ("CM", "Endocardium", "Endothelium"))
    assert "Epicardium" in r["support"]["A_K_excluded_no_share_corrected"]
    assert "Epicardium" not in r["support"]["A_K_support"]
    # the other groups' r_d as a joint call gives them
    dirs, _, _ = sr.feng_dirs.directions_from_stats(fst, genes, cfg, sr.BRACKET)
    halves = sr.feng_dirs.direction_halves(fst, genes, cfg, sr.BRACKET, groups=["CM"])
    m = np.isfinite(dirs["CM"])
    assert r["r_d"]["CM"] == sr.contrast_reliability(dirs["CM"], halves["CM"]["d_a"], halves["CM"]["d_b"], m)


REVIEW3_SENTENCES = (
    "On synthetic data with the truth at 0.4, the cell-sampling intervals are about +-0.002 wide and L_corr reads 3-6% "
    "low (12-15% at scale 1.5) from the pseudocount and the two datasets' different depths; a fraction statement whose "
    "margin is smaller than that bias says nothing about 0.4.",
    "At a pace near 0.2 the corrected pace statement turns on differences of 1-3%, inside the read's own systematic "
    "error on synthetic data; 'the day holds more than a fifth' with p5 within 0.01 of 0.2 is not a finding.",
    "The intervals are a few thousandths wide at the released cell numbers; every statement is decided by the read's "
    "systematic errors (pseudocount, label purity, one scale for between-type and within-type change), which the "
    "intervals do not hold.",
    "The raw pace is biased by the ratio of the direction's to the contrast's reliability; a 'not separated' that "
    "comes from raw against corrected is that bias, not uncertainty.",
    "A is read on the embryo libraries alone; any difference between E9.5's and E10.5's embryo portions (dissection, "
    "composition within a label) that lies along Feng's direction is read as pace and subtracted from K.",
    "L is lowered by any impurity of the challenge's labels or any pooling of Feng's labels (10% mislabelled cells "
    "lowered it by a fifth on synthetic data), and A and L rest on different groups and weights; an 'above both "
    "courses' can come from L alone.",
    "0.4 is compared with the courses in share units; M1's own 0.2 step reads as s_cal in those units, so the fraction "
    "M1 would need is the course times 0.2 / s_cal, which is reported and decides nothing; neither unit says what "
    "scores.",
    "No statement of T47 measures the test's interval, sample kind or scoring; 'above both courses' means only that "
    "the released embryo days, at the E9.5 between-type scale, move less along Feng's direction than five even days "
    "would, and it is consistent with T44's finding that larger steps score better.",
)


def test_the_third_reviews_ledger_sentences_are_in_the_docstring_in_its_words(sr):
    for s in REVIEW3_SENTENCES:
        assert s in sr.LEDGER_NOTES
    assert sr.LEDGER_NOTES[-1].startswith("The self-check stands Feng in with the challenge's own E9.5 cells")
