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
    # knockout, conditional, knockdown, guide, degron, other alleles and spellings (stage 1)
    "Ctnnb1 KO", "Mesp1-Cre; Ctnnb1 F/F", "Ctnnb1 fl/fl", "ctnnb1-/-", "Ctnnb1 cKO", "CTNNB1 null",
    "Ctnnb exon3 deletion", "Catnb lox", "beta-catenin knockdown", "β-catenin shRNA", "β-cat cKO", "ß-catenin",
    "ß-cat", "bcat KO", "b-cat fl", "B_catenin", "beta catenin", "sgCtnnb1", "guide: sgRNA-Ctnnb1-2",
    "Ctnnb1-dTAG degron", "Ctnnb1 siRNA", "Gata4 KO", "GATA-4 null", "gata 4 siRNA", "Gata4 flox/flox", "sgGata4",
    "Gata4-AID degron", "shGata4", "Gata4_cKO",
    "beta-catenin inhibitor", "Gata4 inhibitor", "CHIR99021", "CHIR-99021 3uM", "XAV939", "XAV-939", "IWR-1",
    "IWR1-endo", "ICG-001", "PNU-74654", "iCRT14", "IWP-2", "IWP2", "LGK974", "LGK-974", "Wnt-C59", "WntC59",
    # the reviewer's blocker strings (stage 2)
    "GATA binding protein 4", "GATA-binding factor 4", "sgbcat", "sgβcat", "shBcat",
    # symbols in any case, with separators, orthologues, allele and line names
    "CTNNB1", "ctnnb1", "Ctnn-b1", "CTNN B1", "ctnnb2 morphant", "GATA4", "gata4", "Gata-4", "GATA 4", "gata_4",
    "GATA4/GATA-4", "Ctnnb1tm2Kem", "Ctnnb1tm1Mmt", "Catnb^lox(ex3)", "Ctnnb1 lox(ex3)/+", "Ctnnb1 del(ex3)",
    "Ctnnb1 delta ex3", "Ctnnb1(ex3)fl", "Bfc/+ (batface)",
    # full names
    "catenin beta 1", "Catenin beta-1", "catenin (cadherin-associated protein), beta 1", "catenin, beta",
    "cadherin associated protein beta", "cadherin-associated protein, beta 1", "GATA binding protein 4 null",
    "transcription factor GATA-4", "GATA-binding protein 4", "gata binding factor-4", "GATA IV",
    # beta-catenin short names and encodings
    "beta-catenin", "β-catenin", "ß-catenin", "b-catenin", "bcat", "β-cat", "Β-CATENIN", "betacatenin", "βcatenin",
    "beta-cat", "bCat", "B-Cat", "ss-catenin", "ꞵ-catenin", "ϐ-catenin", "beta—catenin", "beta‑catenin",
    "β​catenin", "be­ta-catenin", "Î²-catenin", "%CE%B2-catenin", "&beta;-catenin", "&#946;-catenin",
    "beta-catenins", "anti-β-catenin antibody", "nuclear β-cat",
    # stabilised or dominant-active forms
    "deltaN beta-catenin", "ΔN89 β-catenin", "β-catenin S33Y", "S33Y", "S37A mutant", "T41A", "S45F", "S45del",
    "stabilised β-catenin", "dominant-active bcat", "Ctnnb1 ca",
    # guides, hairpins, siRNA, morpholinos, degrons, CRISPRi / a, glued prefixes and suffixes
    "sgbcat_1", "gbcat", "shbcat", "sibcat", "crbcat", "crisprbcat", "crispribcat", "crisprabcat", "MObcat",
    "bcat MO", "β-cat MO", "bcatMO", "bcatKO", "bcatfl/fl", "bcatcKO", "β-cat-AID", "bcat-dTAG", "βcatKD",
    "CRISPRi-Ctnnb1", "CRISPRa Gata4", "dCas9-KRAB sgGata4", "siGata4", "siCtnnb1", "MO-ctnnb1", "Gata4 MO",
    "crGata4", "gRNA_Gata4_1", "Gata4-dTAG", "Gata4 degron", "shRNA_Ctnnb1",
    # compound genotypes
    "Gata4;Gata6", "Gata4/6 dKO", "Gata4+/-; Gata6+/-", "GATA5/4", "Gata6, 4 double", "Gata 5/6/4",
    "Ctnnb1;Apc", "Apc; Ctnnb1 fl/fl", "Hand2;Gata4",
    # numerals, fullwidth and superscript forms
    "Ｇａｔａ４", "Gata⁴", "GATA４", "Ｃｔｎｎｂ１",
    # Gata4 alleles and the reporter knock-in (held out: fails closed)
    "Gata4 flox", "Gata4 fl/+", "G4 fl/fl", "G4 KO", "G4 -/-", "G4 knockout", "Gata4-H2B-GFP knock-in reporter",
    "Gata4tm1.1Sad", "Gata4 G295S",
    # identifiers
    "ENSMUSG00000006932", "ENSG00000168036", "ENSMUSG00000021944", "ENSG00000136574.15",
    "Entrez gene ID: 12387", "GeneID:1499", "ncbi gene 14463", "gene_id=2626",
    # drugs on beta-catenin or the TCF complex
    "CHIR98014", "PRI-724", "FH535", "PKF115-584", "CWP232291", "tegavivint", "BC2059", "KYA1797K", "MSAB",
    "WNT974",
    # the fly orthologue
    "armadillo", "arm (armadillo) RNAi",
    # stage 3: every string the two reviewers of 1172f53 found passing (or reaching review only) that is held out
    "bcatGOF", "bcatLOF", "bcatCA", "bcatf/f", "bcatiKO", "bcatiCKO", "bcatKI", "bcatcre", "β-catGOF", "βcatCA",
    "βcatLOF", "Mesp1Cre;bcatGOF", "CAbcat", "CAβcat", "deltaNbcat", "ΔNbcat", "iBcat", "shRNAbcat", "sgRNAbcat",
    "GATA binding 4", "GATA-binding 4", "Gata1-4", "Gata 2-6 TKO", "Gata four KO", "CT99021", "CT-99021",
    "b-catenine", "beta-catenine", "beta catenine", "bet-catenin", "G4 KD", "G4 siRNA", "G4 MO", "G4 deficient",
    "G4 mutant", "plate B cat", "CTNNB1;GATA4;SOX17 library", "guides: Sox2,Ctnnb1,T", "targets=Gata4",
    "lib_Gata4_g1", "βcatf/f", "bcatF/F", "bcatfx/fx", "bcatfloxed", "β-catfloxed", "bcatknockout", "bcatmutant",
    "bcatdeficient", "bcatlof", "bcatexon3", "bcatcKOs", "caBcat", "oeBcat", "lentiBcat", "batface", "Batface allele",
    "G4f/f", "G4fx", "GATA binding protein four", "Gata-four", "βcateninfloxed", "sgG4", "shG4", "Gata 1 to 4",
    "G4 morphant", "ibeta-catenin", "TGFbeta-catenin",
    # homoglyphs: a Cyrillic 'а' and 'С'
    "Gаta4", "Сtnnb1 KO",
]
REVIEW = ["Gata6 F/+", "GATA-6 het", "Gata5 KO", "GATA VI", "G4 ESC chimera", "Apc min", "ApcMin/+", "Apc fl",
          "adenomatous polyposis coli", "Axin2 KO", "Axin1/2", "conductin", "Gsk3b inhibitor", "Gsk-3",
          "glycogen synthase kinase 3", "LiCl 10 mM", "lithium chloride", "kenpaullone", "SB216763", "SB-415286",
          "TWS119", "AZD1080", "LY2090314", "BIO 2uM", "6-bromoindirubin-3'-oxime", "BIO-acetoxime", "Wnt agonist",
          "CHIR treated", "IWP", "IWR", "XAV", "LGK", "C59", "C-59", "tankyrase inhibitor", "Tnks1/2", "TNKS2",
          "G007-LK", "G244-LM", "NVP-TNKS656", "JW55", "Porcn cKO", "porcupine inhibitor", "ETC-159", "GNF-6231",
          "Wnt signalling inhibition", "inhibitor of Wnt secretion", "Wnt pathway blocked", "Wnt loss",
          "loss of Wnt signalling", "Wnt3a", "Wnt3 null", "WNT5A KO", "R-spondin", "Rspo3", "Dkk1", "Dickkopf",
          "Sfrp1", "Lrp5", "Lrp6 KO", "Lrp5/6", "Tcf7 KO", "Tcf7l1", "Tcf7l2 null", "Tcf1/Lef1 dKO", "Lef1",
          "TCF/LEF", "Tcf3", "Wls fl/fl", "Gpr177", "Wntless", "phenocopy", "phenocopies the mutant",
          "lox(ex3)", "delta ex3", "Δex3", "ex3 flox", "exon 3 del", "ΔN", "deltaN", "pannier",
          "alpha-catenin KO",
          # stage 3: the reviewers' review-level strings, and the over-refusals that must stay clearable
          "Ctnb1 KO", "Cttnb1 KO", "Ctnbb1 fl/fl", "Cnntb1", "Ctnnnb1", "Gat4 KO", "Gtaa4", "Chiron 3 µM", "2i/LIF",
          "2i medium", "3i", "Chi pulse", "3 uM Chi", "dnTCF", "dnTcf4", "dnLef1", "dnTCF7L2", "ΔNTcf3", "caTcf",
          "dominant negative Tcf", "pyrvinium", "NSC668036", "niclosamide", "Wnt-off", "Wnt low", "pangolin RNAi",
          "arm RNAi", "arm mutant", "pangolin", "dTCF", "sFRP", "CP21R7", "IM-12", "Chip B, Cat# 1000073",
          "Buffer B (Cat. No. 123)", "Kit B cat no 1000", "LB catalogue", "subcategory", "bobcat",
          "TGFbeta (Cat# 1)", "APC-conjugated anti-CD45", "APC-Cy7", "Wnt1-Cre", "ICG staining", "catenin alpha",
          "delta-catenin", "Gata6", "GATA binding protein 6", "Atg4 KO", "G4 ES cells"]
CLEAN = ["Tbx1Cre/flox; Rosa26-GFPflox/+", "wild type", "E8.5_CRISPR_T_WT", "E8.5_CRISPR_T_KO", "Hand2 null",
         "Mef2c KO", "Dnmt3a+3b", "Eed+Cdkn2a", "Gata1 chimera", "whole embryo", "C57BL/6",
         "bio replicate 1", "biological replicate 2", "conducting", "ΔNp63", "abcam", "Chromium Single Cell 3' v2",
         "Mesp1Cre; Tbx1 flox/flox", "Dnmt1 KO E8.5 rep2", "armadillo repeat protein Armc5", "category",
         "embryonic day 8.5", "10x Genomics", "L3mbtl2 KO", "Kmt2b null", "Rnf2", "G9a / Ehmt2",
         # stage 3: the reviewers' clean strings that stay clean
         "Gata1 KO", "Gata2", "Gata3", "Gata2/3 DKO", "Gata1, Gata2", "Gata 1-3", "Ctnna1", "Ctnnd1", "Mesp1-Cre",
         "Isl1 KO", "Tbx1 KO", "Hand1", "Nkx2-5 KO", "chirality", "E8.5 wild type", "centrifuged at 500 g at 4 C"]
# glued affixes (no separator) for the second property test (the second reviewer's list)
AFFIXES = ["f/f", "fx", "floxed", "mutant", "knockout", "lof", "deficient", "i", "ca", "oe", "lenti", "deltaN", "ΔN",
           "shRNA", "sgRNA"]
# the separators an insertion uses in the property test (glued to letters, a short name cannot be told from a word)
JOINS = [" ", "_", "; ", ", ", "-", "/", " | ", "\t", "(", ")"]
CLEAN_RECORD = {"title": ["E8.5_embryo_3"], "source_name_ch1": ["whole embryo"],
                "characteristics_ch1": ["genotype: wild type", "developmental stage: E8.5"],
                "description": ["10x Chromium 3' v3 library"], "extract_protocol_ch1": ["dissociated with TrypLE"],
                "supplementary_file_1": ["ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1_matrix.mtx.gz"]}


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


def test_the_clean_record_is_clean():
    assert sc.screen_record("GSM1", CLEAN_RECORD)["verdict"] == "pass"


def _insertions(value, spelling, join):
    mid = len(value) // 2
    return [spelling + join + value, value + join + spelling, value[:mid] + join + spelling + join + value[mid:]]


@pytest.mark.parametrize("i", range(len(HELD)))
def test_property_a_held_out_spelling_in_any_field_is_caught(i):
    """Any held-out spelling inserted into any field of a clean record (start, end, middle; rotating separators)
    makes the record held_out, and in a cell table removes that cell only."""
    spelling = HELD[i]
    for j, field in enumerate(CLEAN_RECORD):
        for join in (JOINS[(i + j) % len(JOINS)], JOINS[(i + j + 5) % len(JOINS)]):
            for text in _insertions(CLEAN_RECORD[field][0], spelling, join):
                r = sc.screen_record("GSM1", {**CLEAN_RECORD, field: [text]})
                assert r["verdict"] == "held_out", (field, text)
    for col in ("genotype", "note"):
        for text in _insertions("wild type", spelling, JOINS[i % len(JOINS)]):
            obs = pd.DataFrame({"genotype": ["wild type"] * 2, "note": ["wild type"] * 2}, index=["a", "b"])
            obs.loc["a", col] = text
            keep, _ = sc.screen_cells(obs, ["genotype"])
            assert keep.tolist() == [False, True], (col, text)


# G4 forms glued to letters are REVIEW, not held out: a letter before G4 is as often a gene that ends in g4 (Atg4,
# Gng4) and letters after a short allele word are ordinary text ('G4 mouse'); a held-out verdict there could never be
# cleared. They must still refuse.
GLUE_REVIEW_ONLY = {s for s in HELD if s.lower().startswith(("g4", "sgg4", "shg4"))}


@pytest.mark.parametrize("i", range(len(HELD)))
def test_property_a_held_out_spelling_glued_to_any_affix_is_held_out(i):
    """The second reviewer's extension: each held spelling glued with no separator to a rotating affix (f/f, fx,
    floxed, mutant, knockout, lof, deficient, i, ca, oe, lenti, deltaN, ΔN, shRNA, sgRNA), before and after, is held
    out (the G4 forms: refused, at least for review)."""
    spelling, n = HELD[i], len(AFFIXES)
    for a in (AFFIXES[i % n], AFFIXES[(i + 4) % n]):
        for text in (a + spelling, spelling + a):
            m = sc.match(text)
            if spelling in GLUE_REVIEW_ONLY:
                assert m is not None, text
            else:
                assert m is not None and m[0] == "held_out", (text, m)
            r = sc.screen_record("GSM1", {**CLEAN_RECORD, "title": [text]})
            assert r["verdict"] != "pass", text


def test_catalogue_wording_is_review_and_clearable_but_a_glued_short_name_is_not():
    for text, term in (("Chip B, Cat# 1000073", "b, cat#"), ("Buffer B (Cat. No. 123)", "b (cat."),
                       ("Kit B cat no 1000", "b cat no"), ("subcategory", "subcat"), ("bobcat", "bobcat")):
        rec = {"extract_protocol_ch1": [text]}
        assert sc.screen_record("GSM1", rec)["verdict"] == "review", text
        ok = [{"sample": "GSM1", "term": term, "reason": "catalogue wording in a kit name"}]
        assert sc.screen_record("GSM1", rec, review_cleared=ok)["verdict"] == "pass", text
    with pytest.raises(cm.Refused):
        sc.check_review_cleared([{"sample": "GSM1", "term": "bcat", "reason": "x"}])


def test_an_extended_form_voids_every_review_clearance_of_its_sample():
    ok = [{"sample": "GSM1", "term": "gata6", "reason": "a neighbour only"}]
    assert sc.screen_record("GSM1", {"title": ["Gata6 F/+"]}, review_cleared=ok)["verdict"] == "pass"
    for other in ("Ctnb1 KO", "Gtaa4 null", "deltaNbeta-catenin control", "b-catenine"):
        r = sc.screen_record("GSM1", {"title": ["Gata6 F/+"], "description": [other]}, review_cleared=ok)
        assert r["verdict"] != "pass" and r["clearance_void"], other
    for term in ("ctnb1", "gtaa4", "catenine", "ctnnnb1"):
        with pytest.raises(cm.Refused):
            sc.check_review_cleared([{"sample": "GSM1", "term": term, "reason": "a typo of something else"}])
    with pytest.raises(cm.Refused):
        sc.check_review_cleared([{"sample": "", "term": "gata6", "reason": "x"}])
    with pytest.raises(cm.Refused):
        sc.check_review_cleared([{"sample": "GSM1", "term": "gata6"}])


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
    for term in ("gata4", "ctnnb1", "beta-catenin", "chir99021", "gata binding protein 4", "bcat", "s33y",
                 "ensmusg00000006932", "not a term at all"):
        with pytest.raises(cm.Refused):
            sc.check_review_cleared([{"sample": "GSM1", "term": term, "reason": "x"}])
    with pytest.raises(cm.Refused):
        sc.check_review_cleared([{"sample": "GSM1", "term": "gata6", "reason": " "}])
    # clearing one review term never passes a second one in the same text, nor a held-out term beside it
    two = {"characteristics_ch1": ["genotype: Gata6 F/+; Axin2 KO"]}
    assert sc.screen_record("GSM1", two, review_cleared=ok)["verdict"] == "review"
    both = {"characteristics_ch1": ["genotype: Gata6 F/+; Gata4 F/+"]}
    assert sc.screen_record("GSM1", both, review_cleared=ok)["verdict"] == "held_out"


def test_a_column_name_with_a_screen_term_refuses_the_cell_table():
    obs = pd.DataFrame({"genotype": ["wt", "wt"], "Ctnnb1_fl": [0, 1]}, index=["a", "b"])
    with pytest.raises(cm.Refused, match="column name"):
        sc.screen_cells(obs, ["genotype"])
    obs = pd.DataFrame({"genotype": ["wt", "wt"], "note": ["Gata6 F/+", "x"]}, index=["a", "b"])
    cleared = [{"sample": "GSM1", "term": "gata6", "reason": "a reason"}]
    with pytest.raises(cm.Refused, match="review"):
        sc.screen_cells(obs, ["genotype"], review_cleared=cleared, sample_id="GSM2")
    keep, _ = sc.screen_cells(obs, ["genotype"], review_cleared=cleared, sample_id="GSM1")
    assert keep.tolist() == [True, True]


def test_extra_refuse_terms_only_add():
    rec = {"title": ["Hoxz fl/fl"]}
    assert sc.screen_record("GSM1", rec)["verdict"] == "pass"
    assert sc.screen_record("GSM1", rec, extra_refuse=["hoxz"])["verdict"] == "held_out"
    assert sc.screen_record("GSM1", {"title": ["Wls fl/fl"]})["verdict"] == "review"   # a stage-2 review term


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
        "genotype": ["Hoxz fl/+", "Hoxz fl/fl", "Ctnnb1 fl/fl", "", None, "Hoxz fl/fl", "Hoxz fl/+"],
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
    for acc in ("GSE168768", "GSE186934"):
        assert "is not a phenocopy of Gata4 loss" in ds[acc]["status_reason"]
        assert "enter no model" not in ds[acc]["status_reason"]
    assert ds["GSE205117"]["exclude_match_on"] == "title"
    assert 0 < cm.qc_floors()["min_confident_share_per_sample"] <= 1
    with pytest.raises(cm.Refused):
        cm.dataset("GSE165300")
    with pytest.raises(cm.Refused):
        cm.dataset("GSE000001")
    # the held and excluded accessions are exactly NEVER_FETCH; every cleared dataset has a pinned block (empty until
    # the two readers fill it)
    assert {a for a, d in ds.items() if d["status"] != "cleared"} == cm.NEVER_FETCH
    for a in cleared:
        assert ds[a]["pinned"] == [] or all(p["arm"] in cm.PIN_ARMS for p in ds[a]["pinned"])
    assert "Nomaru et al. (per the scout)" in ds["GSE186934"]["citation"]


def test_the_committed_dataset_list_hashes_to_its_pin_in_t3x_common():
    """Every edit of datasets.yaml (a pin, a clearance, a selector) must update t3x_common.DATASETS_SHA256 in the SAME
    commit: fetch, process and combine refuse the committed file otherwise. `fetch.py datasets-sha` prints both."""
    import hashlib

    got = hashlib.sha256(cm.DATASETS.read_bytes()).hexdigest()
    assert got == cm.DATASETS_SHA256, (
        f"experiments/t3_ext/datasets.yaml hashes to {got} but t3x_common.DATASETS_SHA256 is {cm.DATASETS_SHA256}: "
        "the pin-and-commit step must put the new value into t3x_common.DATASETS_SHA256 in the same commit")


def test_a_wrong_hash_pin_refuses_the_committed_file(monkeypatch):
    monkeypatch.setattr(cm, "DATASETS_SHA256", "0" * 64)
    for f in (cm.load_datasets, cm.qc_floors, cm.datasets_sha256):
        with pytest.raises(cm.Refused, match="same commit"):
            f()
    with pytest.raises(cm.Refused, match="same commit"):
        cm.load_datasets(cm.DATASETS)   # the committed file by its path is checked too
