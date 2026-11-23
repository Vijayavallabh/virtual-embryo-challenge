"""experiments/t1_ext groups.py and tables.py on synthetic data: the group score against a hand case, relabel()
running on every derived config, the extra-class rule, the checks, the consensus sign rule and weights (zeros count as
no value; k = 1 groups keep the ruled vector byte for byte) and the norm match. The strict form (ruling row R27)
deleted the composition block: its tests went with it, and no output of the member or the run holds a composition or
share key (test_no_composition_key_in_any_output)."""

from __future__ import annotations

import copy
import importlib.util
import json
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
    assert not hasattr(gr, "share_check")  # deleted with the composition block
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


# keys that would carry a composition (a share of cells, a level or change form, label counts): none may be written
COMPOSITION_KEY = re.compile(r"composition|(?<!zero_)share|level|change_eligible|label_counts|forecast|logit")


def keys_of(obj, out=None):
    out = set() if out is None else out
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.add(str(k))
            keys_of(v, out)
    elif isinstance(obj, list | tuple):
        for v in obj:
            keys_of(v, out)
    return out


def composition_keys(obj):
    return sorted(k for k in keys_of(obj) if COMPOSITION_KEY.search(k))


def synthetic_member():
    cfg = copy.deepcopy(cm.load_atlases())
    acc = "GSE999004"
    cfg["atlases"][acc] = {"status": "use", "series": [acc], "form": "tenx_mtx", "early_tissue": "whole_embryo",
                           "early": {"labels": ["E8.5"], "tissue": "embryo"},
                           "late": {"labels": ["E14.5"], "tissue": "heart"}, "groups_formed": ["CM"],
                           "qc": {"min_umi": 1, "min_genes": 1, "max_mito": 1.0}}
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
                 "feng_late_profiles": late_p,
                 "purity": {"E8.5": {"cm": {c: 0.0 for c in comp}, "own": {c: 100.0 for c in comp}}}}
    t1_genes = ["Fill7", "Myl7", "Afp", "Hbb-y", "Sox1", "Fill1", "NotInAtlas"]
    return tb.atlas_member(acc, cfg, red, lab, cfg_a, t1_genes, checks_in)


def test_atlas_member_on_a_synthetic_whole_embryo_atlas():
    m = synthetic_member()
    assert m["record"]["groups_kept"] == ["CM"], m["record"]["checks"]
    d = m["directions"]["CM"]
    assert np.isnan(d[[2, 3, 4, 6]]).all()  # competing-class genes held in a whole-embryo atlas; unmapped NaN
    assert m["record"]["fractions"]["rescale_to_5_days"] == 5.0 / 6.0
    assert math.isclose(m["record"]["fractions"]["test"], 2.0 / 6.0)
    assert set(m["record"]["purity_unchecked_classes"]) == {"CM"}


def rec_zero(out):
    return json.loads((out / "directions.json").read_text())["consensus"]["raw"]["CM"]["genes_zero_one_member"]


def test_run_wiring_on_stubs(tmp_path, monkeypatch):
    """run() with stub inputs (no N54 file, no released data): a group only the atlas in use forms passes through
    byte for byte in every variant, the consensus group is written, the fallback and norm records are present, and
    the manifest says tasks ['T1']. Every output goes outside the repository."""

    ld = _load("load")
    acc = "GSE999005"
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["atlases"][acc] = {"status": "use", "series": [acc], "form": "tenx_mtx", "early_tissue": "whole_embryo",
                           "early": {"labels": ["E8.5"], "tissue": "embryo"},
                           "late": {"labels": ["E14.5"], "tissue": "heart"}, "groups_formed": ["CM"],
                           "qc": {"min_umi": 1, "min_genes": 1, "max_mito": 1.0}}
    genes = ["G0", "G1", "G2", "G3"]
    ruled = {"CM": np.array([0.1, -0.2, 0.3, NAN]), "Epicardium": np.array([0.3, 0.0, 0.7, -0.1])}
    feng = {"atlas": tb.IN_USE, "early": "E9.5", "late": "E14.5", "early_day": 9.5, "late_day": 14.5,
            "ruled": ruled, "late_profiles": {},
            "directions": {"CM": np.array([0.1, -0.1, 0.2, NAN]), "Epicardium": np.array([0.2, 0.0, 0.5, -0.1])},
            "record": {"ruled_record": {"fallbacks": {"Epicardium": {"takes": "Endothelium"}}}}}
    member = {"atlas": acc, "early": "E8.5", "late": "E14.5", "early_day": 8.5, "late_day": 14.5,
              "directions": {"CM": np.array([0.2, -0.1, 0.4, 0.1])}, "record": {"groups_kept": ["CM"]}}
    monkeypatch.setattr(cm, "require_frozen", lambda c: None)
    monkeypatch.setattr(gr, "load_extra", lambda c: {})
    monkeypatch.setattr(gr, "derived_config", lambda a, c, extra, e, late: ({"stages": [e, late]}, "b" * 64))
    monkeypatch.setattr(gr, "label_atlas", lambda red, c: {})
    monkeypatch.setattr(ld, "read_reduced", lambda a, d, c=None: {
        "manifest": {"chosen": {"early": "E8.5", "late": "E14.5"}, "inputs": {"fetch_manifest_sha256": "f" * 64}},
        "manifest_sha256": "a" * 64, "gsm": np.array(["GSM1", "GSM2"])})
    monkeypatch.setattr(tb, "t1_genes_from", lambda p: genes)
    monkeypatch.setattr(tb, "released_inputs", lambda c, d, cfgs: {"purity": {}, "categories": {},
                                                                    "released_profiles": {}})
    monkeypatch.setattr(tb, "feng_member", lambda c, g: feng)
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
    # no composition: no composition file, and no composition or share key in any output or in the return value
    assert sorted(p.name for p in out.iterdir()) == ["MANIFEST.json", "directions.json", "directions.npz"]
    for f in ("MANIFEST.json", "directions.json"):
        assert composition_keys(json.loads((out / f).read_text())) == [], f
    assert composition_keys(res) == [] and not any(COMPOSITION_KEY.search(k) for k in z.files)
    with pytest.raises(cm.Refused):
        tb.run(cfg, ROOT / "data" / "derived" / "t1_ext_test", [acc], tmp_path / "released", "stub")


def test_no_composition_key_in_any_output():
    """The strict form deleted the composition block (ruling row R27): the member record of an outside atlas and of the
    atlas in use carry no share, level, change form or label counts, the tables module has no composition function,
    and atlases.yaml has no composition block (the run's outputs are checked in test_run_wiring_on_stubs)."""
    m = synthetic_member()
    assert "composition" not in m and composition_keys(m) == []
    for name in ("cm_share", "composition_table", "forecast", "logit", "MODES"):
        assert not hasattr(tb, name), name
    cfg = cm.load_atlases()
    assert "composition" not in cfg and "share_logit_max" not in cfg["checks"]
    assert composition_keys(cfg["atlases"]) == []
    import inspect

    assert list(inspect.signature(tb.feng_member).parameters) == ["cfg", "t1_genes"]
    src = (EXT / "tables.py").read_text()
    assert "composition.json" not in src and "released_share" not in src
