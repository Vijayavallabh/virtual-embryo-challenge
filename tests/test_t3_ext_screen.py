"""The held-out-gene screen and the stage cap of experiments/t3_ext, on synthetic labels only."""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import pandas as pd
import pytest

EXT = Path(__file__).resolve().parents[1] / "experiments" / "t3_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t3x_common as cm  # noqa: E402
import t3x_screen as sc  # noqa: E402

HELD = [
    # knockout, conditional, knockdown, guide, degron, other alleles and spellings
    "Ctnnb1 KO", "Mesp1-Cre; Ctnnb1 F/F", "Ctnnb1 fl/fl", "ctnnb1-/-", "Ctnnb1 cKO", "CTNNB1 null",
    "Ctnnb exon3 deletion", "Catnb lox", "beta-catenin knockdown", "β-catenin shRNA", "β-cat cKO", "ß-catenin",
    "ß-cat", "bcat KO", "b-cat fl", "B_catenin", "beta catenin", "sgCtnnb1", "guide: sgRNA-Ctnnb1-2",
    "Ctnnb1-dTAG degron", "Ctnnb1 siRNA", "Gata4 KO", "GATA-4 null", "gata 4 siRNA", "Gata4 flox/flox", "sgGata4",
    "Gata4-AID degron", "shGata4", "Gata4_cKO",
    # drugs naming the protein or acting on it (the phenocopy clause)
    "beta-catenin inhibitor", "Gata4 inhibitor", "CHIR99021", "CHIR-99021 3uM", "XAV939", "XAV-939", "IWR-1",
    "IWR1-endo", "ICG-001", "PNU-74654", "iCRT14", "IWP-2", "IWP2", "LGK974", "LGK-974", "Wnt-C59", "WntC59",
]
REVIEW = ["Gata6 F/+", "GATA-6 het", "Apc min", "Axin2 KO", "Axin1/2", "Gsk3b inhibitor", "Gsk-3", "Wnt agonist",
          "CHIR treated", "alpha-catenin KO", "Porcn cKO", "inhibitor of Wnt secretion"]
CLEAN = ["Tbx1Cre/flox; Rosa26-GFPflox/+", "wild type", "E8.5_CRISPR_T_WT", "Hand2 null", "Mef2c KO",
         "subcategory", "Dnmt3a+3b", "Eed+Cdkn2a", "Gata1 chimera", "whole embryo"]


@pytest.mark.parametrize("text", HELD)
def test_every_held_out_spelling_is_caught(text):
    assert sc.match(text)[0] == "held_out"
    assert sc.screen_record("GSM1", {"characteristics_ch1": ["genotype: " + text]})["verdict"] == "held_out"


@pytest.mark.parametrize("text", REVIEW)
def test_review_terms_refuse_for_review(text):
    assert sc.match(text)[0] == "review"


@pytest.mark.parametrize("text", CLEAN)
def test_clean_labels_pass(text):
    assert sc.match(text) is None


def test_record_screens_every_field_and_refuses_an_empty_record():
    rec = {"title": ["E8.5 embryo 3"], "source_name_ch1": ["whole embryo"],
           "supplementary_file_1": ["ftp://x/GSM1_Gata4KO_matrix.mtx.gz"]}
    assert sc.screen_record("GSM1", rec)["verdict"] == "held_out"
    with pytest.raises(cm.Refused):
        sc.screen_record("GSM2", {"title": ["  "]})
    with pytest.raises(cm.Refused):
        sc.screen_record("GSM3", {})


def test_review_clearing_is_per_sample_and_never_clears_a_held_out_term():
    rec = {"characteristics_ch1": ["genotype: Gata6 F/+"]}
    ok = [{"sample": "GSM1", "term": "gata6", "reason": "a reason"}]
    assert sc.screen_record("GSM1", rec, review_cleared=ok)["verdict"] == "pass"
    assert sc.screen_record("GSM2", rec, review_cleared=ok)["verdict"] == "review"
    for term in ("gata4", "ctnnb1", "beta-catenin", "chir99021"):
        with pytest.raises(cm.Refused):
            sc.check_review_cleared([{"sample": "GSM1", "term": term, "reason": "x"}])


def test_extra_refuse_terms_only_add():
    rec = {"title": ["Wls fl/fl"]}
    assert sc.screen_record("GSM1", rec)["verdict"] == "pass"
    assert sc.screen_record("GSM1", rec, extra_refuse=["wls"])["verdict"] == "held_out"


def test_the_screen_has_no_switch(tmp_path):
    for f in (sc.screen_record, sc.screen_cells, sc.match):
        names = " ".join(inspect.signature(f).parameters).lower()
        assert not any(w in names for w in ("enable", "skip", "off", "disable", "bypass", "force"))
    y = tmp_path / "d.yaml"
    y.write_text("version: 1\ndatasets:\n  - {accession: GSE1, status: cleared, tier: wnt, study: S, design: "
                 "germline, control_relation: littermate, stages_used: [E8.5], screen: off}\n")
    with pytest.raises(cm.Refused, match="no switch"):
        cm.load_datasets(y)
    y.write_text("version: 1\ndatasets:\n  - {accession: GSE1, status: cleared, tier: wnt, study: S, design: "
                 "germline, control_relation: littermate, stages_used: [E8.5], review_cleared: "
                 "[{sample: GSM1, term: gata4, reason: x}]}\n")
    with pytest.raises(cm.Refused):
        cm.load_datasets(y)


def test_cells_are_removed_by_label_before_counts_and_unscreenable_tables_refused():
    obs = pd.DataFrame({
        "genotype": ["Wls fl/+", "Wls fl/fl", "Ctnnb1 fl/fl", "", None, "Wls fl/fl", "Wls fl/+"],
        "hashing": ["Bar1", "Bar2", "Bar3", "Bar1", "Bar2", "Negative", "Doublet"],
        "condition": ["ctrl", "ko", "ko", "ctrl", "ko", "ko", "ctrl"],
    }, index=[f"c{i}" for i in range(7)])
    keep, rec = sc.screen_cells(obs, ["genotype"], ["hashing"])
    assert keep.tolist() == [True, True, False, False, False, False, False]
    assert rec["removed"]["held_out:genotype"] == 1
    with pytest.raises(cm.Refused):
        sc.screen_cells(obs, [])
    with pytest.raises(cm.Refused):
        sc.screen_cells(obs.drop(columns="genotype"), ["genotype"])
    drug = obs.assign(treatment=["DMSO"] * 6 + ["XAV939"])
    keep2, _ = sc.screen_cells(drug, ["genotype"])
    assert not keep2.iloc[6]
    with pytest.raises(cm.Refused, match="review"):
        sc.screen_cells(obs.assign(note=["Gata6 F/+"] * 7), ["genotype"])


def test_guides_need_a_map_and_a_held_out_guide_is_removed():
    obs = pd.DataFrame({"genotype": ["pool"] * 3, "sgRNA": ["g1", "g2", "g3"]}, index=["a", "b", "c"])
    with pytest.raises(cm.Refused, match="guide"):
        sc.screen_cells(obs, ["genotype"])
    with pytest.raises(cm.Refused, match="not in the map"):
        sc.screen_cells(obs, ["genotype"], guide_map={"g1": "Hand2", "g2": "Tbx1"})
    keep, _ = sc.screen_cells(obs, ["genotype"], guide_map={"g1": "Hand2", "g2": "Ctnnb1", "g3": "non-targeting"})
    assert keep.tolist() == [True, False, True]


@pytest.mark.parametrize("label,ok", [("E8.5", 8.5), ("E9.25", 9.25), ("E9.5", 9.5), ("E7.75", 7.75)])
def test_stage_cap_keeps_stages_at_or_before_e95(label, ok):
    assert cm.stage_ok(label) == ok


@pytest.mark.parametrize("label", ["E9.75", "E10.5", "E13.5", "E14.5", "P0", "E9.5-E10.5", "E9.5 to E10.5",
                                   "embryo", "", "E8.5/E9.5x"])
def test_stage_cap_fails_closed(label):
    with pytest.raises(cm.Refused):
        cm.stage_ok(label)


def test_the_committed_dataset_list_loads_and_holds_back_what_the_plan_holds_back():
    ds = cm.load_datasets()
    cleared = {a for a, d in ds.items() if d["status"] == "cleared"}
    assert cleared == {"GSE205117", "GSE168768", "GSE186934", "GSE137337", "GSE204906"}
    assert {a for a, d in ds.items() if d["status"] == "held_phenocopy"} == {"GSE165300", "GSE126128", "GSE280543"}
    assert {a for a, d in ds.items() if d["status"] == "excluded_design"} == {"GSE259342", "E-MTAB-7325"}
    assert ds["GSE205117"]["d"] == 0.5 and ds["GSE137337"]["d"] == 0.5 and ds["GSE168768"]["d"] == 0.5
    for d in ds.values():
        for s in d["stages_used"]:
            assert cm.stage_ok(s) <= 9.5
    with pytest.raises(cm.Refused):
        cm.dataset("GSE165300")
    with pytest.raises(cm.Refused):
        cm.dataset("GSE000001")
