"""experiments/t1_m1/share_read.py (try T47 of experiments/harness/TRIES.md) on synthetic statistics and synthetic
AnnData stages built here. No released file, statistics file or score is read."""

from __future__ import annotations

import importlib.util
import json
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
    """Feng-format statistics ('stage|zone|cell_type_spec|half'): means[(stage, label)] is the mean CP10k vector."""
    keys, ns, S = [], [], []
    for st in stages:
        for lab in FENG_LABELS:
            for h in "AB":
                keys.append(f"{st}|z|{lab}|{h}")
                ns.append(n)
                S.append(means[(st, lab)] * n)
    S = np.array(S, dtype=np.float64)
    return {"genes": np.array(genes, dtype=object), "keys": np.array(keys, dtype=object), "n": np.array(ns, float),
            "S_log": np.log1p(S / n) * n, "S_cp": S, "S_nz": np.full_like(S, n),
            "depth_sum": np.full(len(keys), 1e4 * n)}


def stage(blocks: dict, genes) -> ad.AnnData:
    """blocks[(label, library)] = mean CP10k vector; N_CELLS identical cells each, named cell<i>_<library>."""
    rows, labs, names = [], [], []
    i = 0
    for (lab, lib), cp in blocks.items():
        for _ in range(N_CELLS):
            rows.append(np.log1p(cp))
            labs.append(lab)
            names.append(f"c{i}_{lib}")
            i += 1
    a = ad.AnnData(X=sparse.csr_matrix(np.array(rows, dtype=np.float32)),
                   obs=pd.DataFrame({"celltype": labs}, index=names), var=pd.DataFrame(index=list(genes)))
    return a


def world(seed=0, G=60, extra_stage=False):
    """Feng statistics at E9.5 and E14.5 (and E12.5 with extra_stage) on genes g1..; gene 'g0' is absent from Feng (a
    non-finite direction)."""
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(G)]
    feng_genes = genes[1:]
    means = {}
    for lab in FENG_LABELS:
        a = rng.uniform(2.0, 6.0, G - 1)
        means[("E9.5", lab)] = a
        means[("E14.5", lab)] = (a + EPS) * np.exp(rng.normal(0, 0.3, G - 1)) - EPS
        if extra_stage:
            means[("E12.5", lab)] = np.full(G - 1, 1e3)
    stages_f = ("E9.5", "E12.5", "E14.5") if extra_stage else ("E9.5", "E14.5")
    fst = feng_stats(feng_genes, means, stages_f)
    return genes, fst


def build_t1(sr, cfg, genes, fst, c_lib=(0.2, 0.3, 0.9), c_prev=0.25, k=1.5, held_gene=None):
    """Three T1 stages where, for every group, E9.5 -> E10.5 library i moves by c_lib[i] times the group's Feng
    direction, E8.5 -> E9.5 by c_prev times it (CM; at E8.5 Endocardium takes Endothelium's cells and Epicardium has
    none), and the E9.5 contrasts between groups are k times Feng's E9.5 contrasts."""
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
    E105 = {(T1_LABEL[g], li): moved(g, c) for g in T1_LABEL for li, c in zip(("1", "2", "3"), c_lib, strict=True)}
    E85 = {("V-CM", "1"): (e95["CM"] + EPS) * np.exp(-c_prev * d0["CM"]) - EPS,
           ("BEC", "1"): e95["Endothelium"], ("Foregut", "1"): np.full(len(genes), 3.0)}
    return {"E8.5": stage(E85, genes), "E9.5": stage(E95, genes), "E10.5": stage(E105, genes)}, dirs


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
    assert v["b_E9.5_E10.5_lib3/CM/share_corrected"]["value"] == pytest.approx(0.9 / 0.5, abs=1e-3)
    assert np.isnan(v["b_E9.5_E10.5_lib3/Endocardium/share_corrected"]["value"])
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


def test_the_decay_course_is_not_defined_outside_its_range(sr):
    assert np.isnan(sr.courses(-0.1, 1.0)["f_decay"]) and np.isnan(sr.courses(1.2, 1.0)["f_decay"])
    assert sr.courses(1.0, 1.0)["f_decay"] == 0.0


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
    # the range runs from the smaller course's lower percentile to the larger's upper one
    # (the smaller and the larger by point value: here f_decay is the smaller)
    assert sr.fraction_statement(big, _rec(0.6, 0.2, 0.7), _rec(0.3, 0.35, 0.38)) == sr.INSIDE
    assert sr.fraction_statement(big, _rec(0.6, 0.2, 0.7), _rec(0.3, 0.41, 0.45)) == sr.BELOW_BOTH


def test_the_statements_come_out_of_the_reading(sr, cfg):
    genes, fst = world()
    stages, _ = build_t1(sr, cfg, genes, fst, c_lib=(0.18, 0.18, 0.6), k=1.0)
    r = sr.reading(stages, fst, [], {}, cfg, n_boot=2)
    assert r["statements"]["kind"] == sr.KIND_ALONG  # K = 0.42
    assert r["statements"]["pace"] == sr.PACE_FIFTH_OR_LESS  # A / L = 0.18
    assert r["statements"]["entry_fraction"] == sr.INSIDE  # f_decay 0.389, f_even 0.41
    txt = sr.text({"label": sr.LABEL_T47, "read": r, "stamp": {"commit": "abc", "tree_clean": True, "utc": "t"}})
    assert sr.KIND_ALONG in txt and sr.PACE_FIFTH_OR_LESS in txt and "Entry's fraction 0.4: inside" in txt
    json.dumps(r)  # the record serialises (no tuple keys)


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
        seen.append(stage)
        return real(st, stage, labels, halves)

    monkeypatch.setattr(sr.dl, "feng_sel", spy)
    stages, _ = build_t1(sr, cfg, genes, handed)
    sr.reading(stages, handed, [], {}, cfg, n_boot=1)
    assert seen and set(seen) <= set(sr.BRACKET)
