"""experiments/t1_ext groups.py and tables.py (row T72, candidate MX2s) on synthetic data only: the group score against
a hand case; the labeller's configs (row N54's frozen yaml unchanged, plus the skm class on both new atlases, plus the
derived classes on X020's E8.5 samples only), relabel() running on them, no skm gene ever a control; the extra-class
rule; the member tables against hand computations (per-sample CP10k on the cell's Gene Expression total, equal-weight
samples, pseudocount 0.5, the factor 5/8, unique-symbol mapping, ^Hb[ab]- and ^mt- NaN, at least 10 cells at each
end); the skm purity figure; and that tables.py forms no consensus, variant, rescale, control, composition or share.
"""

from __future__ import annotations

import copy
import importlib.util
import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t1_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t1x_common as cm  # noqa: E402


def _load(name):
    key = f"t1_ext_mod_{name}"
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[key] = m
    spec.loader.exec_module(m)
    return m


gr = _load("groups")
tb = _load("tables")
fr = gr.fr  # the labeller of row N54, as groups.py imports it
EXTRA = {"stage": "E8.5",
         "classes": {"neural": {"label": "r_neural", "released_labels": ["Neural Tube"], "genes": ["Sox1", "Pax6"]},
                     "paraxial": {"label": "r_paraxial", "released_labels": ["Paraxial Mesoderm"],
                                  "genes": ["Meox1"]}}}


# ------------------------------------------------------------------------------------------------ the labeller


def test_set_score_hand_case():
    C = sparse.csr_matrix(np.array([[9.0, 1.0, 0.0, 0.0], [0.0, 0.0, 5.0, 5.0]]))
    L = fr.log1p_cp10k(C)
    s = fr.set_score(L, [0], [2, 3])
    want = [math.log1p(9e3) - 0.0, 0.0 - math.log1p(5e3)]
    assert np.allclose(s, want)


def test_label_configs_keep_n54_and_add_skm_and_the_derived_classes_at_e85_only():
    cfg = cm.load_atlases()
    base = gr.base_config(cfg)
    heart = gr.label_configs({"accession": "X019", "early_tissue": "heart"}, cfg, EXTRA, "E8.5", "E16.5")
    emb = gr.label_configs({"accession": "X020", "early_tissue": "whole_embryo"}, cfg, EXTRA, "E8.5", "E16.5")
    assert [r["stages"] for r in heart] == [["E8.5", "E16.5"]]  # one run over both ends, as row N54's
    assert [r["stages"] for r in emb] == [["E8.5"], ["E16.5"]]
    for r in (*heart, *emb):
        c = r["config"]
        assert c["sets"]["skm"] == {"role": "competing", "label": "r_skm",
                                    "genes": {g: "X010_T71" for g in ("Tnnt3", "Tnni2", "Mylpf", "Myl1", "Tnnc2")}}
        for k in ("score", "assign", "groups", "terminal", "unassigned_cap", "min_cells_ext", "fallback"):
            assert c[k] == base[k], k  # N54's thresholds, map and rules unchanged
        assert {n: s for n, s in c["sets"].items() if n in base["sets"]} == base["sets"]
    assert "neural" in emb[0]["config"]["sets"] and "paraxial" not in emb[0]["config"]["sets"]  # under min_genes_found
    assert not {"neural", "paraxial"} & set(emb[1]["config"]["sets"])  # the late end: skm only
    assert not {"neural", "paraxial"} & set(heart[0]["config"]["sets"])
    assert len({r["sha256"] for r in (*heart, *emb)}) == 3
    with pytest.raises(cm.Refused):
        gr.label_configs({"accession": "x", "early_tissue": "whole_embryo"}, cfg, None, "E8.5", "E16.5")
    with pytest.raises(cm.Refused):  # derived classes of another stage
        gr.label_configs({"accession": "x", "early_tissue": "whole_embryo"}, cfg, {**EXTRA, "stage": "E9.5"},
                         "E8.5", "E16.5")
    assert gr.SKM_PURITY_SCORE == 0.3


def synthetic_reduction(seed=3, n_per=60, genes_extra=()):
    cfg = cm.load_atlases()
    runs = gr.label_configs({"accession": "X020", "early_tissue": "whole_embryo"}, cfg, EXTRA, "E8.5", "E16.5")
    named = sorted({g for r in runs for s in r["config"]["sets"].values() for g in s["genes"]})
    genes = named + [f"Fill{i}" for i in range(400)] + list(genes_extra)
    rng = np.random.default_rng(seed)
    n = 4 * n_per
    X = rng.poisson(0.3, size=(n, len(genes))).astype(np.float64)
    cmi = [genes.index(g) for g in runs[0]["config"]["sets"]["cm"]["genes"]]
    skm = [genes.index(g) for g in gr.SKM_GENES]
    for k in range(4):
        X[k * n_per:k * n_per + 25][:, cmi] += 40  # 25 cardiomyocytes per sample
    X[3 * n_per + 25:3 * n_per + 35][:, skm] += 60  # ten skeletal-muscle-like cells in a late sample
    stage = np.array(["E8.5"] * (2 * n_per) + ["E16.5"] * (2 * n_per))
    gsm = np.repeat(["GSM1", "GSM2", "GSM3", "GSM4"], n_per)
    red = {"C": sparse.csr_matrix(X), "genes": genes, "stage": stage, "gsm": gsm,
           "total_ge": X.sum(1) + 5.0}  # the cell's total on all its Gene Expression rows (more than these genes)
    return red, runs


def test_relabel_runs_on_each_config_and_no_skm_gene_is_a_control():
    red, runs = synthetic_reduction()
    lab = gr.label_atlas(red, runs)
    for k in range(4):
        assert np.mean(lab["labels"][k * 60:k * 60 + 25] == "r_cm") > 0.9
    assert np.all(np.isfinite(lab["skm"]))
    assert np.mean(lab["labels"][3 * 60 + 25:3 * 60 + 35] == "r_skm") > 0.9
    assert [r["stages"] for r in lab["runs"]] == [["E8.5"], ["E16.5"]]
    out = fr.relabel(red["C"], red["genes"], red["stage"], gr.label_configs(
        {"accession": "X019", "early_tissue": "heart"}, cm.load_atlases(), None, "E8.5", "E16.5")[0]["config"])
    assert not any(set(s["controls"]) & set(gr.SKM_GENES) for s in out["sets"].values())


def test_extra_class_rule():
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["extra_competing"]["classes"] = {"neural": {"label": "r_neural", "released_labels": ["Neural Tube"],
                                                    "candidates": ["Sox1", "Pax6", "Hes5", "Myl7"]}}
    means = {"Neural Tube": {"Sox1": 1.0, "Pax6": 0.5, "Hes5": 0.2, "Myl7": 2.0},
             "Foregut": {"Sox1": 0.1, "Pax6": 0.3, "Hes5": 0.0, "Myl7": 0.0}}
    r = gr.derive_extra(means, cfg)
    assert "neural" not in r["classes"]  # one member: under min_genes_found
    rec = r["record"]["neural"]
    assert rec["members"] == ["Sox1"]  # Pax6 fails the ratio, Hes5 the mean, Myl7 is named by the base
    assert rec["left_out"]["Myl7"].startswith("named by the base yaml") and "dropped" in rec


# ------------------------------------------------------------------------------------------------ the tables


def test_tables_hand_case():
    # two early samples, one late sample; genes A, B, mt-X, Hbb-y; T1 genes in another order with a duplicate and one
    # gene the atlas lacks
    genes = ["A", "B", "mt-X", "Hbb-y", "Dup"]
    rows = {"GSM1": [[10, 0, 0, 0, 1], [0, 10, 0, 0, 1]],
            "GSM2": [[5, 5, 0, 0, 1]] * 10,
            "GSM3": [[0, 20, 2, 2, 1]] * 12}
    X, gsm, stage, tot = [], [], [], []
    for g, rr in rows.items():
        for r in rr:
            X.append(r)
            gsm.append(g)
            stage.append("E16.5" if g == "GSM3" else "E8.5")
            tot.append(100.0)  # every cell's Gene Expression total
    X = np.array(X, dtype=np.float64)
    red = {"C": sparse.csr_matrix(X), "genes": genes, "stage": np.array(stage), "gsm": np.array(gsm),
           "total_ge": np.array(tot)}
    lab = {"labels": np.array(["r_cm"] * len(gsm)), "skm": np.zeros(len(gsm))}
    t1 = ["B", "Zz", "A", "Dup", "Dup", "mt-X", "Hbb-y"]
    out = tb.member_tables(red, lab, {"CM": ["r_cm"]}, ["CM"], "E8.5", "E16.5", t1)
    t = out["CM"]["table"]
    assert t is not None and list(t) == list(tb.TABLE_KEYS)
    # GSM1: 2 cells, CP10k A = (1000 + 0) / 2 = 500, B = 500; GSM2: A = B = 500; equal weight: m_early A = B = 500
    # GSM3: A = 0, B = 2000
    k = 5.0 / 8.0
    want_A = math.log((0 + 0.5) / (500 + 0.5)) * k
    want_B = math.log((2000 + 0.5) / (500 + 0.5)) * k
    assert t["D"][2] == pytest.approx(want_A) and t["D"][0] == pytest.approx(want_B)
    assert np.isnan(t["D"][1]) and np.isnan(t["D"][3]) and np.isnan(t["D"][4])  # lacked, T1 duplicate
    assert np.isnan(t["D"][5]) and np.isnan(t["D"][6])  # ^mt- and ^Hb[ab]- held
    assert t["m_late"][5] == pytest.approx(200.0)  # the level itself is kept
    assert (int(t["n_cells_early"]), int(t["n_cells_late"])) == (12, 12)
    assert t["per_sample_early"].shape == (2, len(t1)) and list(t["samples_early"]) == ["GSM1", "GSM2"]
    assert t["per_sample_early"][0][2] == pytest.approx(500.0) and t["per_sample_early"][1][2] == pytest.approx(500.0)
    # equal weight per sample, not per cell: GSM1's two cells weigh as much as GSM2's ten
    lab2 = {"labels": np.array(["r_cm"] * len(gsm)), "skm": np.zeros(len(gsm))}
    X2 = X.copy()
    X2[0] = [30, 0, 0, 0, 1]  # GSM1's first cell: A 3000 -> GSM1 A = 1500; m_early A = (1500 + 500) / 2 = 1000
    red2 = {**red, "C": sparse.csr_matrix(X2)}
    t2 = tb.member_tables(red2, lab2, {"CM": ["r_cm"]}, ["CM"], "E8.5", "E16.5", t1)["CM"]["table"]
    assert t2["m_early"][2] == pytest.approx(1000.0)


def test_a_group_needs_ten_cells_at_each_end_and_a_cell_in_every_sample():
    genes = ["A", "B"]
    X = np.ones((30, 2))
    stage = np.array(["E8.5"] * 20 + ["E16.5"] * 10)
    gsm = np.array(["GSM1"] * 10 + ["GSM2"] * 10 + ["GSM3"] * 10)
    red = {"C": sparse.csr_matrix(X), "genes": genes, "stage": stage, "gsm": gsm, "total_ge": np.full(30, 2.0)}
    labs = np.array(["r_cm"] * 30, dtype=object)
    labs[20:21] = "r_epicardial"  # the late end keeps 9 CM cells
    out = tb.member_tables(red, {"labels": labs.astype(str), "skm": np.zeros(30)}, {"CM": ["r_cm"]}, ["CM"],
                           "E8.5", "E16.5", genes)
    assert out["CM"]["table"] is None and "fewer than 10" in out["CM"]["record"]["no_table"][0]
    labs = np.array(["r_cm"] * 30, dtype=object)
    labs[10:20] = "r_epicardial"  # GSM2 has no CM cell
    out = tb.member_tables(red, {"labels": labs.astype(str), "skm": np.zeros(30)}, {"CM": ["r_cm"]}, ["CM"],
                           "E8.5", "E16.5", genes)
    assert out["CM"]["table"] is None and any("a sample without" in w for w in out["CM"]["record"]["no_table"])
    labs = np.array(["r_cm"] * 30, dtype=object)
    out = tb.member_tables(red, {"labels": labs.astype(str), "skm": np.zeros(30)}, {"CM": ["r_cm"]}, ["CM"],
                           "E8.5", "E16.5", genes)
    assert out["CM"]["table"] is not None and np.allclose(out["CM"]["table"]["D"], 0.0)
    with pytest.raises(cm.Refused):  # a cell without a Gene Expression total
        tb.cp10k_mean(sparse.csr_matrix(np.ones((1, 2))), np.array([0.0]))


def test_the_skm_purity_figure():
    lab = {"labels": np.array(["r_cm", "r_cm", "r_cm", "r_skm", "r_cm"]), "skm": np.array([0.31, 0.3, 0.29, 2.0, 5.0])}
    stage = np.array(["E16.5", "E16.5", "E16.5", "E16.5", "E8.5"])
    p = tb.skm_purity(lab, stage, "E16.5", 0.3)
    assert p == {"late_r_cm_cells": 3, "late_r_cm_skm_at_least": 2, "threshold": 0.3, "share": pytest.approx(2 / 3)}
    assert tb.skm_purity(lab, stage, "E17.5", 0.3)["share"] is None


def test_the_constants_are_the_rows():
    assert (tb.PSEUDOCOUNT, tb.FENG_BRACKET_DAYS, tb.MIN_CELLS, tb.CP10K) == (0.5, 5.0, 10, 1e4)
    assert tb.HELD == (r"^Hb[ab]-", r"^mt-")
    assert tb.direction(np.array([0.0]), np.array([0.0]), 8.5, 16.5)[0] == 0.0
    with pytest.raises(cm.Refused):  # an end inside the T1 window
        tb.direction(np.array([0.0]), np.array([0.0]), 8.5, 12.5)


FORBIDDEN = re.compile(r"\bdef (?:consensus|norm_match|group_tables|rescale|mask|perm)\b|VARIANTS|noX002|noE85start|"
                       r"composition|eb_shrink|deciding")


def test_tables_forms_no_consensus_variant_control_or_share():
    src = (EXT / "tables.py").read_text()
    code = "\n".join(ln for ln in src.splitlines() if not ln.lstrip().startswith("#"))
    body = code[code.index('"""', code.index('"""') + 3) + 3:]  # after the module docstring
    assert not FORBIDDEN.search(body)
    assert not hasattr(tb, "consensus") and not hasattr(tb, "VARIANTS")
    gsrc = (EXT / "groups.py").read_text()
    assert "share_check" not in gsrc and "profile_check" not in gsrc
