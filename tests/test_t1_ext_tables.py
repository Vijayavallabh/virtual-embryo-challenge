"""experiments/t1_ext groups.py and tables.py on synthetic data: the group score against a hand case, relabel()
running on every derived config, the extra-class rule, the checks, the consensus sign rule and weights (zeros count as
no value; k = 1 groups keep the ruled vector byte for byte), the norm match, and the composition rules on hand-checked
shares, including an atlas that lacks a group."""

from __future__ import annotations

import copy
import importlib.util
import math
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
NAN = np.nan
EXTRA = {"classes": {"neural": {"label": "r_neural", "released_labels": ["Neural Tube"], "genes": ["Sox1", "Pax6"]},
                     "paraxial": {"label": "r_paraxial", "released_labels": ["Paraxial Mesoderm"],
                                  "genes": ["Meox1"]}}}


# ------------------------------------------------------------------------------------------------ groups


def test_set_score_hand_case():
    C = sparse.csr_matrix(np.array([[9.0, 1.0, 0.0, 0.0], [0.0, 0.0, 5.0, 5.0]]))
    L = fr.log1p_cp10k(C)
    s = fr.set_score(L, [0], [2, 3])
    want = [math.log1p(9e3) - 0.0, 0.0 - math.log1p(5e3)]
    assert np.allclose(s, want)


def derived(early_tissue):
    cfg = cm.load_atlases()
    a = {"accession": "GSE999003", "early_tissue": early_tissue}
    return gr.derived_config(a, cfg, EXTRA, "E8.5", "E14.5")


def test_derived_configs_keep_the_base_and_add_extra_classes_only_for_embryos():
    heart, h_sha = derived("heart")
    emb, e_sha = derived("whole_embryo")
    base = gr.base_config(cm.load_atlases())
    assert heart["stages"] == ["E8.5", "E14.5"] and heart["sets"] == base["sets"]
    assert heart["terminal"] == base["terminal"] and heart["assign"] == base["assign"]
    assert "neural" in emb["sets"] and "paraxial" not in emb["sets"]  # one gene: under min_genes_found
    assert h_sha != e_sha
    with pytest.raises(cm.Refused):
        gr.derived_config({"accession": "x", "early_tissue": "whole_embryo"}, cm.load_atlases(), None, "E8.5",
                          "E14.5")


@pytest.mark.parametrize("tissue", ["heart", "whole_embryo"])
def test_relabel_runs_on_each_derived_config(tissue):
    cfg_a, _ = derived(tissue)
    named = sorted({g for s in cfg_a["sets"].values() for g in s["genes"]})
    fillers = [f"Fill{i}" for i in range(400)]
    genes = named + fillers
    rng = np.random.default_rng(3)
    n = 120
    X = rng.poisson(0.3, size=(n, len(genes))).astype(np.float64)
    cm_idx = [genes.index(g) for g in cfg_a["sets"]["cm"]["genes"]]
    X[:40][:, cm_idx] += 40  # forty cardiomyocytes
    stage = np.array(["E8.5"] * 60 + ["E14.5"] * 60)
    red = {"C": sparse.csr_matrix(X), "genes": genes, "stage": stage, "gsm": np.where(stage == "E8.5", "G1", "G2")}
    lab = gr.label_atlas(red, cfg_a)
    assert np.mean(lab["labels"][:40] == "r_cm") > 0.9
    assert set(lab["per_sample"]) == {"G1", "G2"}
    ps = lab["per_sample"]["G1"]
    assert ps["unassigned_no_class"] + ps["unassigned_margin"] == ps["labels"]["r_unassigned"]


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
    assert rec["left_out"]["Myl7"] == "named by the base yaml" and "dropped" in rec


def test_checks_hand_cases():
    genes = [f"G{i}" for i in range(30)]
    base = np.linspace(0.5, 30, 30)
    grp = {"genes": genes, "cp": base, "log": np.log1p(base)}
    refs = {"V-CM": {"genes": genes, "cp": base * 1.1, "log": np.log1p(base * 1.1)},
            "Neural Tube": {"genes": genes, "cp": base[::-1], "log": np.log1p(base[::-1])}}
    r = gr.profile_check(grp, refs, ["V-CM"], 1.0)
    assert r["passed"] and r["best"] == "V-CM"
    assert not gr.profile_check(grp, refs, ["Neural Tube"], 1.0)["passed"]
    p = gr.purity_gate({"hepatic": 0.1, "neural": 0.6}, {"hepatic": 0.0, "neural": 0.0},
                       {"hepatic": 1.0, "neural": 1.0, "mural": None})
    assert not p["passed"] and p["classes"]["hepatic"]["passed"] and not p["classes"]["neural"]["passed"]
    assert not p["classes"]["mural"]["checked"]
    assert gr.share_check(0.5, 0.5, 1.0)["passed"] and not gr.share_check(0.9, 0.5, 1.0)["passed"]
    assert gr.embryo_gene_gate(np.array([2.0, 1.0, NAN]), np.array([1.0, 1.0, 0.0])).tolist() == [True, False, False]


# ------------------------------------------------------------------------------------------------ consensus


def test_sign_rule_and_equal_weights():
    m = {"A": np.array([1.0, 1.0, 1.0, 0.0, NAN, 1.0, -1.0]),
         "B": np.array([2.0, -1.0, 1.0, 0.0, NAN, NAN, -3.0]),
         "C": np.array([3.0, -1.0, NAN, 0.0, NAN, NAN, 2.0]),
         "D": np.array([NAN, -1.0, NAN, 2.0, NAN, NAN, NAN])}
    r = tb.consensus(m)
    # gene 0: k 3, all agree: mean of 1, 2, 3 over k = 2.0
    # gene 1: k 4, three negative agree (> k/2): (-3) / 4; the dissent counts as zero
    # gene 2: k 2, both agree: 1.0
    # gene 3: zeros are no value: k 1, one cannot agree alone: 0.0
    # gene 4: no finite value: NaN; gene 5: one member: 0.0
    # gene 6: k 3, two negative agree (> 1.5): (-1 - 3) / 3
    want = [2.0, -0.75, 1.0, 0.0, NAN, 0.0, -4.0 / 3.0]
    assert np.allclose(r["c"], want, equal_nan=True)
    assert r["k"].tolist() == [3, 4, 2, 1, 0, 1, 3]
    two = tb.consensus({"A": np.array([1.0, 1.0]), "B": np.array([1.0, -1.0])})
    assert two["c"].tolist() == [1.0, 0.0]  # k = 2: both must agree


def test_norm_match_and_passthrough_byte_for_byte():
    ref = np.array([3.0, 4.0, NAN])
    # the second review's rule: the ruled norm over genes where c is finite and non-zero (and the ruled finite)
    c, rec = tb.norm_match(np.array([1.0, 0.0, 5.0]), ref)
    assert np.allclose(c, [3.0, 0.0, 15.0]) and rec["factor"] == 3.0 and rec["genes_in_norm"] == 1
    assert math.isclose(rec["zero_share"], 1 / 3)
    # the reviewer's hand case: one surviving gene (0.5 and 0.4 in two members) no longer takes the whole norm
    ruled0 = np.array([3.0, 2.0, 0.45])
    c2, rec2 = tb.norm_match(np.array([0.0, 0.0, 0.45]), ruled0, {"A": np.array([0.1, -0.2, 0.5]),
                                                                    "B": np.array([NAN, 0.1, 0.4])})
    assert np.allclose(c2, [0.0, 0.0, 0.45]) and rec2["factor"] == 1.0
    assert math.isclose(rec2["max_ratio_to_largest_member"], 0.9)
    ruled = {"CM": np.array([0.1, -0.2, NAN]), "Epicardium": np.array([0.3, 0.0, 0.7])}
    members = {"GSE193346": {"CM": np.array([0.1, -0.1, NAN]), "Epicardium": np.array([0.2, 0.0, 0.5])},
               "X1": {"CM": np.array([0.2, -0.1, 0.4])}, "X2": {"CM": np.array([0.3, 0.1, 0.4])}}
    early = {"GSE193346": 9.5, "X1": 8.5, "X2": 9.5}
    out = tb.group_tables(members, ruled, early)
    for v in tb.VARIANTS:
        assert out["tables"][v]["Epicardium"].tobytes() == ruled["Epicardium"].tobytes()
        assert out["record"][v]["Epicardium"]["passthrough"]
    raw = out["tables"]["raw"]["CM"]
    assert np.allclose(raw, [0.2, -0.1 * 2 / 3, 0.4], equal_nan=True)
    dec = out["tables"]["deciding"]["CM"]
    m = np.isfinite(dec) & (dec != 0) & np.isfinite(ruled["CM"])
    assert math.isclose(np.linalg.norm(dec[m]), np.linalg.norm(ruled["CM"][m]))
    r = out["record"]["deciding"]["CM"]
    assert r["norm_match"]["max_ratio_to_largest_member"] is not None and "genes_zero_one_member" in r
    # without the E8.5-start member, two members remain; without the atlas in use, two remain
    assert out["record"]["noE85start"]["CM"]["members"] == ["GSE193346", "X2"]
    assert out["record"]["noX002"]["CM"]["members"] == ["X1", "X2"]


def test_held_mask():
    m = tb.held_mask(["Hba-a1", "Hbb-bs", "mt-Co1", "Myl7", "Afp"], ["^Hb[ab]-", "^mt-"], {"Afp"})
    assert m.tolist() == [True, True, True, False, True]


# ------------------------------------------------------------------------------------------------ composition

RULES = {"heart_labels": ["r_cm", "r_endocardial_ec", "r_vascular_ec", "r_ec_ambiguous", "r_lymphatic",
                          "r_epicardial", "r_mural"], "mesenchyme_label": "r_mesenchyme", "mesenchyme_from_day": 14.0,
         "level_days": [14.0, 15.0], "level_halfrange_max": 1.0, "min_level_atlases": 2}


def lg(p):
    return math.log(p / (1 - p))


def test_cm_share_change_form_uses_one_label_set_and_level_form_adds_mesenchyme_late():
    counts = {"r_cm": 30, "r_endocardial_ec": 10, "r_mesenchyme": 20, "r_blood": 99, "r_unassigned": 50}
    assert tb.cm_share(counts, 9.5, RULES) == 0.75
    assert tb.cm_share(counts, 14.5, RULES) == 0.75           # change form: never r_mesenchyme
    assert tb.cm_share(counts, 14.5, RULES, "level") == 0.5   # level form: r_mesenchyme from E14
    assert tb.cm_share(counts, 9.5, RULES, "level") == 0.75
    assert tb.cm_share({"r_blood": 5}, 14.5, RULES) is None
    with pytest.raises(cm.Refused):
        tb.cm_share(counts, 14.5, RULES, "other")


def test_levels_use_the_level_share_when_given():
    e = [{**entries()[0], "late_share_level": 0.35}, {**entries()[1], "late_share_level": 0.25}]
    t = tb.composition_table(e, RULES)
    assert math.isclose(t["L"], (lg(0.35) + lg(0.25)) / 2)
    assert math.isclose(t["slopes_per_day"]["F"], (lg(0.4) - lg(0.5)) / 5)  # the change form keeps late_share


def entries():
    return [{"atlas": "F", "early_day": 9.5, "late_day": 14.5, "early_share": 0.5, "late_share": 0.4,
             "level_eligible": True, "change_eligible": True},
            {"atlas": "A", "early_day": 8.5, "late_day": 14.5, "early_share": None, "late_share": 0.3,
             "level_eligible": True, "change_eligible": False},
            # an atlas that lacks the CM group at its late end: left out of the levels
            {"atlas": "B", "early_day": 9.5, "late_day": 14.5, "early_share": None, "late_share": None,
             "level_eligible": True, "change_eligible": False},
            # outside the level days
            {"atlas": "N", "early_day": 8.5, "late_day": 16.5, "early_share": 0.6, "late_share": 0.2,
             "level_eligible": True, "change_eligible": True}]


def test_composition_rules_on_hand_checked_shares():
    t = tb.composition_table(entries(), RULES)
    assert set(t["levels_logit"]) == {"F", "A"} and "B level" in t["notes"]
    assert math.isclose(t["L"], (lg(0.4) + lg(0.3)) / 2)
    assert math.isclose(t["half_range"], (lg(0.4) - lg(0.3)) / 2) and t["level_ok"]
    sF, sN = (lg(0.4) - lg(0.5)) / 5, (lg(0.2) - lg(0.6)) / 8
    assert math.isclose(t["slope"], (sF + sN) / 2)
    f = tb.forecast(t, 0.6, "test")
    assert math.isclose(f["x_change"], 1 / (1 + math.exp(-(lg(0.6) + 2 * t["slope"]))))
    assert math.isclose(f["x_level"], 1 / (1 + math.exp(-(lg(0.6) + 0.5 * (t["L"] - lg(0.6))))))
    assert f["offered"]  # both steps negative
    p = tb.forecast(t, 0.6, "practice")
    assert math.isclose(p["step_change"], t["slope"]) and math.isclose(p["step_level"], 0.2 * (t["L"] - lg(0.6)))
    up = tb.forecast(t, 0.2, "test")  # base below L: the level form rises, the change form falls
    assert not up["offered"] and "disagree in sign" in " ".join(up["reasons"])


def test_composition_gate_refuses_one_level_or_a_wide_range():
    e = entries()[:1]
    t = tb.composition_table(e, RULES)
    assert not t["level_ok"] and t["change_label"] == "one atlas"
    assert not tb.forecast(t, 0.6, "test")["offered"]
    wide = [{**entries()[0], "late_share": 0.8}, entries()[1]]
    t2 = tb.composition_table(wide, RULES)
    assert t2["half_range"] > 1.0 and not t2["level_ok"]
    with pytest.raises(cm.Refused):
        tb.logit(1.0)


def test_atlas_member_on_a_synthetic_whole_embryo_atlas():
    cfg = copy.deepcopy(cm.load_atlases())
    acc = "GSE999004"
    cfg["atlases"][acc] = {"status": "use", "series": [acc], "form": "tenx_mtx", "early_tissue": "whole_embryo",
                           "early": {"labels": ["E8.5"], "tissue": "embryo"},
                           "late": {"labels": ["E14.5"], "tissue": "heart"}, "groups_formed": ["CM"],
                           "qc": {"min_umi": 1, "min_genes": 1, "max_mito": 1.0}, "level_eligible": True}
    cfg_a, _ = gr.derived_config({**cfg["atlases"][acc], "accession": acc}, cfg, EXTRA, "E8.5", "E14.5")
    named = sorted({g for s in cfg_a["sets"].values() for g in s["genes"]})
    genes = named + [f"Fill{i}" for i in range(400)]
    rng = np.random.default_rng(5)
    stage = np.array(["E8.5"] * 80 + ["E14.5"] * 80)
    X = rng.poisson(0.5, size=(160, len(genes))).astype(np.float64)
    cmi = [genes.index(g) for g in cfg_a["sets"]["cm"]["genes"]]
    is_cm = np.zeros(160, dtype=bool)
    is_cm[:40] = True
    is_cm[80:120] = True
    X[np.ix_(np.flatnonzero(is_cm), cmi)] += 30
    X[80:120, genes.index("Fill7")] += 20  # a gene that rises in CM between the ends
    red = {"C": sparse.csr_matrix(X), "genes": genes, "stage": stage,
           "gsm": np.array(["G1", "G2"] * 40 + ["G3", "G4"] * 40), "cell": np.array([f"c{i}" for i in range(160)])}
    lab = gr.label_atlas(red, cfg_a)
    rel = np.where(is_cm, "V-CM", "Neural Tube")
    early_p = gr.label_profiles(red["C"][:80], genes, rel[:80])
    late_p = gr.label_profiles(red["C"][80:], genes, np.where(is_cm[80:], "CM", "Endothelium"))
    comp = [n for n, s in cfg_a["sets"].items() if s["role"] == "competing"]
    checks_in = {"categories": {"V-CM": "cm", "Neural Tube": "surrounding"}, "released_profiles": {"E8.5": early_p},
                 "released_share": {}, "feng_late_profiles": late_p,
                 "purity": {"E8.5": {"cm": {c: 0.0 for c in comp}, "own": {c: 100.0 for c in comp}}}}
    t1_genes = ["Fill7", "Myl7", "Afp", "Hbb-y", "Sox1", "Fill1", "NotInAtlas"]
    m = tb.atlas_member(acc, cfg, red, lab, cfg_a, t1_genes, checks_in)
    assert m["record"]["groups_kept"] == ["CM"], m["record"]["checks"]
    d = m["directions"]["CM"]
    assert np.isnan(d[[2, 3, 4, 6]]).all()  # competing-class genes held in a whole-embryo atlas; unmapped NaN
    assert m["record"]["fractions"]["rescale_to_5_days"] == 5.0 / 6.0
    assert math.isclose(m["record"]["fractions"]["test"], 2.0 / 6.0)
    assert m["composition"]["early_share"] is None and m["composition"]["late_share"] > 0.5
    assert m["composition"]["late_share_level"] is not None
    assert set(m["record"]["purity_unchecked_classes"]) == {"CM"}


def rec_zero(out):
    import json

    return json.loads((out / "directions.json").read_text())["consensus"]["raw"]["CM"]["genes_zero_one_member"]


def test_run_wiring_on_stubs(tmp_path, monkeypatch):
    """run() with stub inputs (no N54 file, no released data): a group only the atlas in use forms passes through
    byte for byte in every variant, the consensus group is written, the fallback and norm records are present, and
    the manifest says tasks ['T1']. Every output goes outside the repository."""
    import json

    ld = _load("load")
    acc = "GSE999005"
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["atlases"][acc] = {"status": "use", "series": [acc], "form": "tenx_mtx", "early_tissue": "whole_embryo",
                           "early": {"labels": ["E8.5"], "tissue": "embryo"},
                           "late": {"labels": ["E14.5"], "tissue": "heart"}, "groups_formed": ["CM"],
                           "qc": {"min_umi": 1, "min_genes": 1, "max_mito": 1.0}}
    genes = ["G0", "G1", "G2", "G3"]
    ruled = {"CM": np.array([0.1, -0.2, 0.3, NAN]), "Epicardium": np.array([0.3, 0.0, 0.7, -0.1])}
    comp_f = {"atlas": tb.IN_USE, "early_day": 9.5, "late_day": 14.5, "early_share": 0.5, "late_share": 0.4,
              "late_share_level": 0.35, "level_eligible": True, "change_eligible": True, "share_check": None}
    feng = {"atlas": tb.IN_USE, "early": "E9.5", "late": "E14.5", "early_day": 9.5, "late_day": 14.5,
            "ruled": ruled, "late_profiles": {},
            "directions": {"CM": np.array([0.1, -0.1, 0.2, NAN]), "Epicardium": np.array([0.2, 0.0, 0.5, -0.1])},
            "record": {"ruled_record": {"fallbacks": {"Epicardium": {"takes": "Endothelium"}}}},
            "composition": comp_f}
    member = {"atlas": acc, "early": "E8.5", "late": "E14.5", "early_day": 8.5, "late_day": 14.5,
              "directions": {"CM": np.array([0.2, -0.1, 0.4, 0.1])}, "record": {"groups_kept": ["CM"]},
              "composition": {**comp_f, "atlas": acc, "early_day": 8.5, "early_share": None,
                              "change_eligible": False}}
    monkeypatch.setattr(cm, "require_frozen", lambda c: None)
    monkeypatch.setattr(gr, "load_extra", lambda c: {})
    monkeypatch.setattr(gr, "derived_config", lambda a, c, extra, e, late: ({"stages": [e, late]}, "b" * 64))
    monkeypatch.setattr(gr, "label_atlas", lambda red, c: {})
    monkeypatch.setattr(ld, "read_reduced", lambda a, d: {
        "manifest": {"chosen": {"early": "E8.5", "late": "E14.5"}, "inputs": {"fetch_manifest_sha256": "f" * 64}},
        "manifest_sha256": "a" * 64, "gsm": np.array(["GSM1", "GSM2"])})
    monkeypatch.setattr(tb, "t1_genes_from", lambda p: genes)
    monkeypatch.setattr(tb, "released_inputs", lambda c, d, cfgs: {"released_share": {"E9.5": 0.5}, "purity": {},
                                                                    "categories": {}, "released_profiles": {}})
    monkeypatch.setattr(tb, "feng_member", lambda c, g, s: feng)
    monkeypatch.setattr(tb, "atlas_member", lambda *a, **k: member)
    res = tb.run(cfg, tmp_path, [acc], tmp_path / "released", "stub")
    out = tmp_path / "tables" / "stub"
    z = np.load(out / "directions.npz")
    for v in tb.VARIANTS:
        assert z[f"{v}|Epicardium"].tobytes() == ruled["Epicardium"].tobytes()
    assert np.allclose(z["raw|CM"], [0.15, -0.1, 0.3, 0.0])  # G3: one member has a value, so 0
    assert rec_zero(out) == 1
    rec = json.loads((out / "directions.json").read_text())
    assert rec["consensus"]["deciding"]["CM"]["in_use_member_fallback"] is None
    assert rec["consensus"]["deciding"]["Epicardium"]["in_use_member_fallback"] == {"takes": "Endothelium"}
    assert rec["consensus"]["deciding"]["CM"]["norm_match"]["factor"] is not None
    man = json.loads((out / "MANIFEST.json").read_text())
    assert man["tasks"] == ["T1"] and man["atlases"] == [tb.IN_USE, acc]
    assert res["kept"][acc] == ["CM"]
    with pytest.raises(cm.Refused):
        tb.run(cfg, ROOT / "data" / "derived" / "t1_ext_test", [acc], tmp_path / "released", "stub")
