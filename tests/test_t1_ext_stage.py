"""experiments/t1_ext: the stage rule (second version, after the two reviews of 2026-11-20) and the plan and fetch on
synthetic GEO records (no network). The reviewers' own cases are replayed: a window stage under a key that does not
name a stage, in the description or in pooled wording (blockers 1 and 6); stage tokens in forms the patterns do not
parse (blocker 2); the plug = E1 convention and unrecognised plug sentences (blocker 3); a processing sentence no
negative pattern catches, refused until its text is cleared by hash (blocker 4). Also: edges, integer labels, the
widening when no convention is stated, Theiler stages, somite counts, postnatal ages, a missing stage, conflicting
fields, the canonical processing text, and that only USED samples' files are fetched, after a live re-plan.

Third version (2026-11-21): the pins are the first line. Every case of the verifier's "These end USED" list now ends
WINDOW or UNKNOWN by the stage rule alone; and, independently, a sample the rule reads as permitted is not fetched or
loaded without a pin, is refused when its live record differs from its pin by one character, and is refused when
pinned for a stage its record does not state. Clean records of permitted stages are USED once pinned."""

from __future__ import annotations

import copy
import gzip
import importlib.util
import json
import sys
import textwrap
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t1_ext"


def _load(name):
    key = f"t1_ext_mod_{name}"
    if key in sys.modules:
        return sys.modules[key]
    if str(EXT) not in sys.path:
        sys.path.insert(0, str(EXT))
    spec = importlib.util.spec_from_file_location(key, EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[key] = m
    spec.loader.exec_module(m)
    return m


if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t1x_common as cm  # noqa: E402

fetch = _load("fetch")
E05 = ["Timed matings: the day of the vaginal plug was designated E0.5."]


def mark(fields, files=(), protocol=E05):
    return cm.stage_reading(fields, files, protocol)["mark"]


# ------------------------------------------------------------------------------------------------ the stage rule


@pytest.mark.parametrize("label,expected", [
    ("E9.5", "PERMITTED"),          # exactly E9.5: allowed (open lower edge), with plug = E0.5 stated
    ("E9.51", "WINDOW"),
    ("E13.5", "WINDOW"),            # the closed upper edge
    ("E13.6", "PERMITTED"),
    ("E14.5", "PERMITTED"),
    ("E8.5", "PERMITTED"),
    ("E12.5", "WINDOW"),
    ("embryonic day 10.5", "WINDOW"),
    ("14.5 dpc", "PERMITTED"),
    ("dpc 11.5", "WINDOW"),
    ("E14.5.", "PERMITTED"),        # a sentence's full stop after the number
    ("E12,5", "WINDOW"),            # a decimal comma
    ("e14p5", "PERMITTED"),
])
def test_half_day_labels_with_e05_stated(label, expected):
    assert mark({"characteristics_ch1": [f"developmental stage: {label}"]}) == expected


def test_lower_edge_is_pinned():
    assert not cm.touches_window(9.5, 9.5)
    assert cm.touches_window(9.5, 9.51) and cm.touches_window(13.5, 13.5) and not cm.touches_window(13.51, 14.0)


@pytest.mark.parametrize("label,expected", [("E9", "PERMITTED"), ("E10", "WINDOW"), ("E13", "WINDOW"),
                                            ("E14", "PERMITTED")])
def test_integer_label_is_a_half_day_interval(label, expected):
    r = cm.stage_reading({"title": [label]}, (), E05)
    assert r["mark"] == expected
    if expected == "PERMITTED":
        assert r["nominal_high"] - r["nominal_low"] == 0.5


@pytest.mark.parametrize("label,expected", [("E8.5", "PERMITTED"), ("E9.5", "WINDOW"), ("E14.5", "WINDOW"),
                                            ("E15.5", "PERMITTED"), ("E16.5", "PERMITTED"), ("P3", "PERMITTED")])
def test_unstated_convention_widens_by_minus_one_to_plus_half(label, expected):
    r = cm.stage_reading({"title": [f"{label} sample"]}, (), [])
    assert r["convention"]["name"] == "unstated" and r["convention"]["shift"] == cm.UNSTATED == (-1.0, 0.5)
    assert r["mark"] == expected
    if label == "E14.5":
        assert (r["low"], r["high"]) == (13.5, 15.0) and r["nominal_low"] == 14.5


def test_plug_e0_convention():
    proto = ["Timed matings: the morning of the vaginal plug was designated E0."]
    assert cm.convention(proto)["shift"] == (0.0, 0.5)
    assert mark({"title": ["E9.5 embryo"]}, protocol=proto) == "WINDOW"   # [9.5, 10.0]
    r = cm.stage_reading({"title": ["E14.5 heart"]}, (), proto)
    assert r["mark"] == "PERMITTED" and (r["low"], r["high"]) == (14.5, 15.0) and r["nominal_low"] == 14.5
    r = cm.stage_reading({"title": ["E13"]}, (), ["the morning of the plug was counted as day 0"])
    assert r["mark"] == "WINDOW"


@pytest.mark.parametrize("sentence", ["The day of plug was designated as E1.",
                                      "Noon on the day a vaginal plug was found was counted as day 1.",
                                      "The plug day was taken as the first day of gestation."])
def test_plug_e1_convention_is_not_permitted_at_e14_5(sentence):
    """Blocker 3: under plug = E1, E14.5 is [13.5, 14.0] on the E0.5 scale and meets the window."""
    c = cm.convention([sentence])
    assert c["name"] == "E1" and c["shift"] == (-1.0, -0.5)
    r = cm.stage_reading({"title": ["E14.5 heart"]}, (), [sentence])
    assert r["mark"] != "PERMITTED" and r["mark"] == "WINDOW"


@pytest.mark.parametrize("sentence", [
    "Embryos were collected at E14.5 after detection of a vaginal plug.",
    "Females with a plug were dissected 14 days later.",
    "The plug day was E0 in the first cohort and E0.5 in the second.",
])
def test_other_plug_sentences_and_two_conventions_refuse(sentence):
    assert cm.convention([sentence])["shift"] is None
    assert mark({"title": ["E16.5 heart"]}, protocol=[sentence]) == "UNKNOWN"


def test_plug_sentence_without_a_day_is_ignored_and_two_texts_conflict():
    assert cm.convention(["Females were checked for vaginal plugs."])["name"] == "unstated"
    assert cm.convention([*E05, "The plug day was designated E0."])["name"] == "conflict"


@pytest.mark.parametrize("text,expected", [
    ("TS13", "PERMITTED"), ("TS14", "WINDOW"), ("Theiler stage 22", "WINDOW"), ("TS24", "PERMITTED"),
    ("15 somites", "PERMITTED"), ("25 somites", "WINDOW"), ("somites: 30", "WINDOW"), ("18-24 somites", "WINDOW"),
])
def test_theiler_and_somites(text, expected):
    assert mark({"characteristics_ch1": [f"stage: {text}"]}) == expected


def test_somite_table_last_row_reaches_15():
    r = cm.stage_reading({"characteristics_ch1": ["stage: 45 somites"]}, (), E05)
    assert r["nominal_high"] == 15.0 and r["mark"] == "WINDOW"


def test_somite_count_overrides_a_permitted_label():
    assert mark({"title": ["E9.5 embryo"], "characteristics_ch1": ["stage: 26 somites"]}) == "WINDOW"
    assert mark({"title": ["E9.5 embryo"], "characteristics_ch1": ["stage: 16 somites"]}) == "PERMITTED"


@pytest.mark.parametrize("text,expected", [
    ("P3", "PERMITTED"), ("postnatal day 10", "PERMITTED"), ("neonatal", "PERMITTED"), ("8 weeks old", "PERMITTED"),
    ("6-8 weeks", "UNKNOWN"), ("2 weeks", "UNKNOWN"), ("24 months", "PERMITTED"), ("adult heart", "PERMITTED"),
    ("10 days", "WINDOW"), ("E9.5-E10.5", "WINDOW"), ("E14.5-E15.5", "PERMITTED"), ("E13.5 to E14.5", "WINDOW"),
])
def test_postnatal_ages_ranges_weeks_and_bare_days(text, expected):
    assert mark({"characteristics_ch1": [f"age: {text}"]}) == expected


@pytest.mark.parametrize("fields", [{"title": ["FetalHeart_1"]}, {"title": ["embryo"]}, {}])
def test_missing_stage_is_unknown(fields):
    assert mark(fields) == "UNKNOWN"


# blocker 1: a window stage under another key or in the description, beside a title saying E14.5
@pytest.mark.parametrize("extra", [
    {"characteristics_ch1": ["tissue: E12.5 heart"]},
    {"characteristics_ch1": ["gestation: E12.5"]},
    {"description": ["E12.5 hearts, library 2"]},
])
def test_window_stage_anywhere_in_the_record_refuses(extra):
    r = cm.stage_reading({"title": ["E14.5 heart rep2"], **extra}, (), E05)
    assert r["mark"] in ("UNKNOWN", "WINDOW") and r["window_reading"]


# blocker 6: pooled-stage wording
def test_pooled_stage_wording_reads_as_the_interval_it_covers():
    r = cm.stage_reading({"characteristics_ch1": ["tissue: hearts pooled from E12.5 and E14.5 embryos"]}, (), E05)
    assert r["mark"] == "WINDOW" and (r["nominal_low"], r["nominal_high"]) == (12.5, 14.5)
    r = cm.stage_reading({"title": ["FetalHeart_rep1"],
                          "characteristics_ch1": ["tissue: hearts pooled from E12.5 and E14.5 embryos", "age: E14.5"],
                          "description": ["Hearts were collected at E12.5"]}, (), E05)
    assert r["mark"] == "UNKNOWN" and r["window_reading"]
    assert mark({"title": ["Hearts, E14.5 and E16.5 pooled"]}) == "PERMITTED"
    assert cm.stage_reading({"title": ["Hearts, E14.5 and E16.5 pooled"]}, (), E05)["nominal_high"] == 16.5


# blocker 2: forms the patterns do not parse fail closed (beside an E14.5 title)
@pytest.mark.parametrize("fields,files", [
    ({"characteristics_ch1": ["age: 12.5"]}, ()),
    ({"characteristics_ch1": ["age: E-12.5"]}, ()),
    ({"characteristics_ch1": ["age: pc12.5"]}, ()),
    ({"characteristics_ch1": ["tissue: heart E-12.5"]}, ()),
    ({}, ("GSM1_HeartE12_5_matrix.mtx.gz",)),
    ({}, ("GSM1_E125_matrix.mtx.gz",)),
    ({}, ("GSM1_TS21_matrix.mtx.gz",)),
    ({"description": ["collected on embryonic day twelve"]}, ()),
    ({"description": ["gestational week 2"]}, ()),
    ({"characteristics_ch1": ["age (weeks): 2"]}, ()),
])
def test_unparsed_stage_tokens_fail_closed(fields, files):
    r = cm.stage_reading({"title": ["E14.5 heart"], **fields}, files, E05)
    assert r["mark"] in ("UNKNOWN", "WINDOW"), r


def test_bare_number_in_a_stage_keyed_field_is_an_eday():
    r = cm.stage_reading({"characteristics_ch1": ["age: 14.5"]}, (), E05)
    assert r["mark"] == "PERMITTED" and r["nominal_low"] == 14.5
    assert mark({"characteristics_ch1": ["time point: 12"]}) == "WINDOW"
    # a bare number under a key that names no stage is not read and is not a leftover
    assert mark({"title": ["E14.5 heart"], "characteristics_ch1": ["replicate: 2", "batch: 12"]}) == "PERMITTED"


def test_ordinary_record_text_does_not_trip_the_leftover_rule():
    fields = {"title": ["E14.5 heart rep2, PE150, GSE123456"],
              "characteristics_ch1": ["strain: CD1", "tissue: Heart", "developmental stage: E14.5",
                                      "genotype: wild type", "cell type: Epcam+ and Tie2+ cells"],
              "description": ["1E6 cells were loaded; Cre reporter R26mTmG; 10x v3"]}
    files = ["GSM7_barcodes.tsv.gz", "GSM7_features.tsv.gz", "GSM7_matrix.mtx.gz"]
    assert mark(fields, files) == "PERMITTED"


def test_conflicting_fields_and_file_tokens_refuse():
    assert mark({"title": ["E14.5 heart"], "characteristics_ch1": ["age: E16.5"]}) == "UNKNOWN"
    assert mark({"title": ["E8.5 embryo"]}, files=["GSM1_E8_5_matrix.mtx.gz"]) == "PERMITTED"
    assert mark({"title": ["E8.5 embryo"]}, files=["GSM1_E16_5_matrix.mtx.gz"]) == "UNKNOWN"
    assert mark({"title": ["E8.5 embryo"]}, files=["GSM1_E10_5_matrix.mtx.gz"]) == "UNKNOWN"  # never hidden


def test_bracket_ends_checked_with_windows_module():
    cm.bracket_ok(8.5, 14.5)
    cm.bracket_ok(9.5, 14.5)
    with pytest.raises(cm.Refused):
        cm.bracket_ok(9.5, 12.5)
    with pytest.raises(cm.Refused):
        cm.bracket_ok(14.5, 8.5)


# ------------------------------------------------------------------------------------------------ the verifier's list
#
# Every case of the second verifier's "These end USED" list, beside an E14.5 title with plug = E0.5 stated: the stage
# rule alone now gives WINDOW or UNKNOWN (blockers 1, 2 and 3; the "not realistic" three included).

VERIFIER_USED = [
    {"characteristics_ch1": ["tissue: hearts 13.5-14.5 dpc"]},
    {"characteristics_ch1": ["tissue: hearts 13.5 to 14.5 dpc"]},
    {"characteristics_ch1": ["tissue: hearts (12.5 and 14.5 dpc)"]},
    {"description": ["12.5–14.5 dpc hearts"]},
    {"description": ["hearts at 12.5 and 14.5 days post-coitum"]},
    {"characteristics_ch1": ["tissue: hearts from E14.5 and 13.5 embryos"]},
    {"description": ["E14.5, 13.5 hearts"]},
    {"description": ["heart E14.5 (and 12.5)"]},
    {"description": ["ED12.5 heart"]},
    {"description": ["ED 12.5 heart"]},
    {"characteristics_ch1": ["tissue: GD12.5 heart"]},
    {"characteristics_ch1": ["sample group: 12.5d heart"]},
    {"description": ["EmbryonicDay12.5"]},
    {"description": ["Embryonic12.5"]},
    {"description": ["heart_stage12.5"]},
    {"description": ["e.d. 12.5"]},
    {"description": ["midgestation heart"]},
    {"description": ["Ｅ12.5 heart"]},  # a full-width E
    {"extract_protocol_ch1": ["Hearts were dissected from E12.5 embryos."]},
    {"growth_protocol_ch1": ["Embryos were harvested at 12.5 dpc."]},
    {"data_processing": ["E12.5 and E14.5 libraries were processed with Cell Ranger count"]},
    {"label_ch1": ["E12.5"]},
    {"relation": ["BioSample: E12.5 heart"]},
    {"description": ["em12.5"]},
]


@pytest.mark.parametrize("extra", VERIFIER_USED)
def test_every_case_the_verifier_found_used_is_refused_by_the_rule(extra):
    r = cm.record_stage({"title": ["E14.5 heart"], **extra}, E05)
    assert r["mark"] in ("WINDOW", "UNKNOWN"), (extra, r["reason"])


@pytest.mark.parametrize("extra", [
    {"description": ["hearts 13.5-14.5 dpc"]},
    {"description": ["hearts at 12.5 and 14.5 days post-coitum"]},
    {"description": ["ED12.5 heart"]},
    {"description": ["EmbryonicDay12.5"]},
    {"description": ["Ｅ12.5 heart"]},
    {"extract_protocol_ch1": ["Hearts were dissected from E12.5 embryos."]},
    {"data_processing": ["E12.5 and E14.5 libraries were processed with Cell Ranger count"]},
])
def test_a_window_reading_alone_reads_window(extra):
    """Without the E14.5 title the same wordings read WINDOW: they are readings, not only leftovers."""
    r = cm.record_stage({"title": ["heart"], **extra}, E05)
    assert r["mark"] == "WINDOW" and r["window_reading"], r["reason"]


def test_a_protocol_reading_counts_like_any_other():
    """Blocker 3: the sample's own protocol text is read for the stage, sentence by sentence; the convention's own
    tokens are blanked first, so a stated plug convention is not itself a stage."""
    proto = ["Timed matings were set up; the day of the vaginal plug was designated E0.5.",
             "E14.5 hearts were digested in 0.5 mg/ml collagenase for 15 min at 37 C."]
    r = cm.record_stage({"title": ["E14.5 heart"], "extract_protocol_ch1": proto}, [])
    assert r["mark"] == "PERMITTED" and r["convention"]["name"] == "E0.5", r["reason"]
    assert {x["field"] for x in r["readings"]} == {"title", "extract_protocol_ch1"}
    # a protocol shared by every sample that lists every stage of the series refuses each sample (the default until
    # the main session sets another rule)
    shared = ["Hearts were collected at E8.5, E10.5, E14.5 and P3."]
    assert cm.record_stage({"title": ["E14.5 heart"], "extract_protocol_ch1": shared}, E05)["mark"] in (
        "WINDOW", "UNKNOWN")
    # a stage in the sample's own protocol that differs from its label refuses even outside the window
    assert cm.record_stage({"title": ["E14.5 heart"], "extract_protocol_ch1": ["Hearts were collected at E16.5."]},
                           E05)["mark"] == "UNKNOWN"


def test_ordinary_record_still_reads_permitted():
    """Ordinary text in every field the rule now reads does not trip it: accessions, read lengths, versions,
    quantities with units, relations, the instrument, the release date."""
    rec = {"title": ["E14.5 heart rep2, PE150, GSE123456"],
           "geo_accession": ["GSM7000001"], "status": ["Public on Jun 15 2023"],
           "submission_date": ["Apr 12 2023"], "last_update_date": ["Jun 15 2023"], "type": ["SRA"],
           "channel_count": ["1"], "source_name_ch1": ["heart"], "organism_ch1": ["Mus musculus"],
           "taxid_ch1": ["10090"],
           "characteristics_ch1": ["strain: CD1", "tissue: Heart", "developmental stage: E14.5",
                                   "genotype: wild type", "cell type: Epcam+ and Tie2+ cells"],
           "molecule_ch1": ["polyA RNA"],
           "extract_protocol_ch1": ["E14.5 hearts were digested in 0.5 mg/ml collagenase for 15 min at 37 C; cells "
                                    "were resuspended in 0.04% BSA and 1.5 ml tubes.",
                                    "10x Genomics Chromium Single Cell 3' v3 libraries."],
           "description": ["1E6 cells were loaded; Cre reporter R26mTmG; 10x v3", "Cell Ranger v3.0.2"],
           "data_processing": ["Cell Ranger 3.0.2 count per sample against mm10-3.0.0; Seurat resolution 0.5",
                               "Genome_build: mm10"],
           "platform_id": ["GPL24247"], "contact_name": ["A,,Person"], "contact_zip/postal_code": ["92093"],
           "instrument_model": ["Illumina NovaSeq 6000"], "library_selection": ["cDNA"],
           "library_source": ["transcriptomic single cell"], "library_strategy": ["RNA-Seq"],
           "relation": ["BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN12345678",
                        "SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX1234567"],
           "supplementary_file_1": [url("GSM7000001", "barcodes.tsv.gz")], "series_id": ["GSE123456"],
           "data_row_count": ["0"]}
    r = cm.record_stage(rec, E05)
    assert r["mark"] == "PERMITTED" and r["nominal_low"] == r["nominal_high"] == 14.5, r["reason"]


def test_wider_prefixes_and_ranges_read_as_stages():
    for text, span in (("ED14.5", (14.5, 14.5)), ("GD 16.5", (16.5, 16.5)), ("e.d. 8.5", (8.5, 8.5)),
                       ("EmbryonicDay16.5", (16.5, 16.5)), ("14.5-16.5 dpc", (14.5, 16.5)),
                       ("14.5 and 16.5 days post coitum", (14.5, 16.5)), ("Ｅ１４.5", (14.5, 14.5))):
        r = cm.stage_reading({"title": [text]}, (), E05)
        assert r["mark"] == "PERMITTED" and (r["nominal_low"], r["nominal_high"]) == span, (text, r["reason"])


@pytest.mark.parametrize("text", ["E14.5 heart, 13 hearts pooled", "E14 and 12 hearts", "E14.5 hearts (n=12)"])
def test_a_companion_integer_beside_a_reading_refuses(text):
    assert mark({"title": [text]}) == "UNKNOWN"


@pytest.mark.parametrize("sentence", ["The morning after mating was designated E1.",
                                      "The day of conception was counted as day 1.",
                                      "Fertilisation was taken as E0."])
def test_a_second_convention_in_a_mating_or_conception_sentence_refuses(sentence):
    """The verifier's lesser point 2: sentences about mating, coitus, copulation, conception or fertilisation are read
    for the convention too; E0.5 with such a sentence giving E1 (or E0) is UNKNOWN."""
    assert cm.convention([*E05, sentence])["name"] == "conflict"
    assert mark({"title": ["E16.5 heart"]}, protocol=[*E05, sentence]) == "UNKNOWN"
    # "post coitum" / "post coitus" is a dpc unit, not a convention sentence
    assert cm.convention(["Embryos were collected at 14.5 days post-coitus."])["name"] == "unstated"


# ------------------------------------------------------------------------------------------------ canonical texts


def test_record_text_canonical_form():
    """Pins the canonical record form stated in t1x_common and the README: a pin's record_sha256 is of these bytes."""
    rec = {"title": ["E14.5 heart"], "characteristics_ch1": ["tissue: heart", "age: E14.5\tx"],
           "extract_protocol_ch1": ["Hearts – dissected"], "relation": ["SRA: x"], "label_ch1": [],
           "supplementary_file_1": ["ftp://h/GSM1_matrix.mtx.gz"], "_other_lines": ["#a = b"]}
    text = cm.record_text("SAMPLE", "GSM1", rec)
    assert text == (
        "t1x-record-v1 SAMPLE GSM1\n"
        "title[1]: E14.5 heart\n"
        "characteristics_ch1[1]: tissue: heart\n"
        "characteristics_ch1[2]: age: E14.5\tx\n"
        "characteristics_ch1[2] escaped: age: E14.5\\tx\n"
        "extract_protocol_ch1[1]: Hearts – dissected\n"
        "extract_protocol_ch1[1] escaped: Hearts \\u2013 dissected\n"
        "relation[1]: SRA: x\n"
        "label_ch1: (none)\n"
        "supplementary_file_1[1]: ftp://h/GSM1_matrix.mtx.gz\n"
        "_other_lines[1]: #a = b\n")
    import hashlib
    assert cm.record_sha256("SAMPLE", "GSM1", rec) == hashlib.sha256(text.encode("utf-8")).hexdigest()
    # stored as pairs (the JSON writer sorts keys): the same text
    assert cm.record_text("SAMPLE", "GSM1", json.loads(json.dumps(cm.stored_record(rec), sort_keys=True))) == text
    # the order of the fields is part of the text
    assert cm.record_text("SAMPLE", "GSM1", dict(reversed(list(rec.items())))) != text
    assert cm.record_text("SERIES", "GSE2", {"title": ["x"]}) == "t1x-record-v1 SERIES GSE2\ntitle[1]: x\n"
    with pytest.raises(cm.Refused):
        cm.record_text("SAMPLE", "GSE2", rec)


def test_processing_text_canonical_form():
    """Pins the canonical form stated in t1x_common and the README: the sha256 a person writes is of these bytes.
    Version 2 adds the series' supplementary file names (the verifier's lesser point 4)."""
    series = [("GSE1", {"summary": ["Hearts."], "overall_design": [],
                        "supplementary_file": ["ftp://h/GSE1_aggr_filtered_feature_bc_matrix.h5",
                                               "ftp://h/GSE1_RAW.tar"]})]
    samples = {"GSM20": {"series": "GSE1", "text": {"data_processing": ["Cell Ranger count\tper sample", "v3"],
                                                     "extract_protocol_ch1": ["Hearts were dissected – 10x"]}},
               "GSM10": {"series": "GSE1", "text": {}}}
    text = cm.processing_text("GSE9", series, samples)
    assert text == (
        "t1x-processing-v2 GSE9\n"
        "SERIES GSE1\n"
        "summary[1]: Hearts.\n"
        "overall_design: (none)\n"
        "series_supplementary_file[1]: GSE1_aggr_filtered_feature_bc_matrix.h5\n"
        "series_supplementary_file[2]: GSE1_RAW.tar\n"
        "SAMPLE GSM10 (GSE1)\n"
        "data_processing: (none)\nextract_protocol_ch1: (none)\ndescription: (none)\n"
        "growth_protocol_ch1: (none)\ntreatment_protocol_ch1: (none)\n"
        "SAMPLE GSM20 (GSE1)\n"
        "data_processing[1]: Cell Ranger count\tper sample\n"
        "data_processing[1] escaped: Cell Ranger count\\tper sample\n"
        "data_processing[2]: v3\n"
        "extract_protocol_ch1[1]: Hearts were dissected – 10x\n"
        "extract_protocol_ch1[1] escaped: Hearts were dissected \\u2013 10x\n"
        "description: (none)\ngrowth_protocol_ch1: (none)\ntreatment_protocol_ch1: (none)\n")
    import hashlib
    assert cm.text_sha256(text) == hashlib.sha256(text.encode("utf-8")).hexdigest()
    assert cm.clearance({"accession": "GSE9"}, cm.text_sha256(text))[0] is False
    assert cm.clearance({"processing_cleared": {"sha256": cm.text_sha256(text)}}, cm.text_sha256(text))[0]
    assert not cm.clearance({"processing_cleared": {"sha256": "0" * 64}}, cm.text_sha256(text))[0]


# ------------------------------------------------------------------------------------------------ plan, pins and fetch

ACC = "GSE999001"
PROC = "Cell Ranger count per sample"


def url(gsm, name):
    return f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{gsm}_{name}"


def sample(gsm, title, stage, tissue, files=("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"),
           tag=None, extra=(), strategy="RNA-Seq", processing=PROC):
    lines = [f"^SAMPLE = {gsm}", f"!Sample_title = {title}", f"!Sample_source_name_ch1 = {tissue}"]
    if stage is not None:
        lines.append(f"!Sample_characteristics_ch1 = developmental stage: {stage}")
    lines += [f"!Sample_characteristics_ch1 = tissue: {tissue}", f"!Sample_library_strategy = {strategy}",
              "!Sample_molecule_ch1 = polyA RNA", "!Sample_extract_protocol_ch1 = 10x Genomics 3' v2 libraries",
              f"!Sample_data_processing = {processing}", *extra]
    for i, f in enumerate(files, 1):
        lines.append(f"!Sample_supplementary_file_{i} = {url(gsm, (tag + '_' if tag else '') + f)}")
    return lines


def soft(samples, design="One library per sample. The day of the vaginal plug was designated E0.5."):
    head = [f"^SERIES = {ACC}", "!Series_title = synthetic", "!Series_summary = Hearts and embryos.",
            f"!Series_overall_design = {design}",
            f"!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE999nnn/{ACC}/suppl/"
            f"{ACC}_RAW.tar"]
    return "\n".join(head + [ln for s in samples for ln in s]) + "\n"


def config(extra=None, early=("E8.5",), late=("E14.5", "E16.5")):
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["atlases"][ACC] = {"status": "use", "wave": 1, "series": [ACC], "form": "tenx_mtx",
                           "early": {"labels": list(early), "tissue": "embryo"},
                           "late": {"labels": list(late), "tissue": "heart"},
                           "assay_rule": {"library_strategy": "RNA-Seq"}, "early_tissue": "whole_embryo",
                           "groups_formed": ["CM"], "qc": {"min_umi": 1, "min_genes": 1, "max_mito": 1.0},
                           "protocol_terms_cleared": {}, "protocol_gate_cleared": {}, **(extra or {})}
    return cfg


BASE = [sample("GSM9990001", "E8.5 embryo rep1", "E8.5", "whole embryo", tag="E8_5"),
        sample("GSM9990002", "E14.5 heart rep1", "E14.5", "heart", tag="E14_5"),
        sample("GSM9990003", "E12.5 heart rep1", "E12.5", "heart", tag="E12_5"),
        # inside a permitted-looking accession range, but its own record states an in-window stage
        sample("GSM9990004", "embryo rep2", "E10.5", "whole embryo"),
        sample("GSM9990005", "P3 heart", "P3", "heart", tag="P3")]
PAIR = ("GSM9990001", "GSM9990002")


def getter_for(text, log):
    def get(u, dst):
        log.append(u)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if u.endswith("_family.soft.gz"):
            dst.write_bytes(gzip.compress(text.encode()))
        else:
            dst.write_bytes(b"x" * 10)
    return get


def from_template(text, clear=True):
    """What two readers paste after reading: the template's YAML with the readers and the date filled in."""
    body = "\n".join(ln for ln in text.splitlines() if not ln.startswith("#"))
    body = body.replace("READER_1", "Reader A").replace("READER_2", "Reader B").replace("YYYY-MM-DD", "2026-11-21")
    d = yaml.safe_load(textwrap.dedent(body))
    if not clear:
        d.pop("processing_cleared")
    return d


def pinned(tmp_path, samples, gsms=PAIR, cfg=None, clear=True, **kw):
    """Plan without pins, print the pin template for the named samples, paste it: the pinned (and, with clear, cleared)
    config."""
    cfg = cfg or config()
    fetch.plan(ACC, tmp_path, getter_for(soft(samples, **kw), []), cfg)
    t = fetch.pin_template(ACC, tmp_path, list(gsms), cfg)
    cfg2 = copy.deepcopy(cfg)
    cfg2["atlases"][ACC].update(from_template(t, clear))
    cm.check_pins(ACC, cm.atlas(ACC, cfg2))
    return cfg2


def plan_pinned(tmp_path, samples, gsms=PAIR, cfg=None, clear=True, **kw):
    cfg2 = pinned(tmp_path, samples, gsms, cfg, clear, **kw)
    return fetch.plan(ACC, tmp_path, getter_for(soft(samples, **kw), []), cfg2), cfg2


def test_plan_without_pins_uses_nothing_and_prints_every_permitted_record(tmp_path):
    log = []
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(BASE), log), config())
    rows = {r["gsm"]: r for r in pl["samples"]}
    assert pl["status"] == "refused" and not pl["chosen"]
    assert [rows[g]["mark"] for g in sorted(rows)] == ["NOT_USED", "NOT_USED", "WINDOW", "WINDOW", "NOT_USED"]
    assert "not pinned" in rows["GSM9990001"]["reason"] and rows["GSM9990002"]["candidate"] == {"end": "late",
                                                                                                "label": "E14.5"}
    assert any("do not cover both ends" in r for r in pl["refusals"])
    assert any(f"series {ACC} is not pinned" in r for r in pl["refusals"])
    assert all(u.endswith("_family.soft.gz") for u in log)  # the plan downloads metadata only
    mdir = tmp_path / "meta" / ACC
    recs = (mdir / "records.txt").read_text()
    for g in ("GSM9990001", "GSM9990002", "GSM9990005"):  # every sample the rule lets through: its full record
        assert rows[g]["record_sha256"] == cm.sha256_file(mdir / "records" / f"{g}.txt")
        assert (mdir / "records" / f"{g}.txt").read_text() in recs
        assert f"supplementary_file_1[1]: {url(g, '')}" in recs
    assert "!Sample" not in recs and "extract_protocol_ch1[1]: 10x Genomics 3' v2 libraries" in recs
    assert rows["GSM9990001"]["record_sha256"] in recs
    assert cm.sha256_file(mdir / "records" / f"{ACC}.txt") == pl["series_record_sha256"][ACC]
    for g in ("GSM9990003", "GSM9990004"):  # window rows: no record, only the fields that gave a reading
        assert "record" not in rows[g] and not (mdir / "records" / f"{g}.txt").exists()
        assert rows[g]["tissue_fields"] == fetch.WITHHELD and "convention_text" not in rows[g]
        assert set(rows[g]["stage_fields"]) == ({"characteristics_ch1", "title"} if g == "GSM9990003"
                                                else {"characteristics_ch1"})
        assert g not in recs
    table = (mdir / "plan_table.tsv").read_text()
    assert recs in table and table.endswith((mdir / "processing_text.txt").read_text())


def test_pinned_plan_passes_and_marks_used(tmp_path):
    pl, cfg = plan_pinned(tmp_path, BASE)
    rows = {r["gsm"]: r for r in pl["samples"]}
    assert pl["status"] == "passed", pl["refusals"]
    assert [rows[g]["mark"] for g in sorted(rows)] == ["USED", "USED", "WINDOW", "WINDOW", "NOT_USED"]
    assert pl["chosen"] == {"early": "E8.5", "late": "E14.5"}
    pins = cfg["atlases"][ACC]["pinned"]
    assert [p["accession"] for p in pins["samples"]] == list(PAIR) and pins["samples"][1]["stage"] == "E14.5"
    assert pins["samples"][0]["read_by"] == ["Reader A", "Reader B"] and pins["series"][0]["accession"] == ACC
    assert cfg["atlases"][ACC]["processing_cleared"]["sha256"] == pl["processing_sha256"]
    assert "series_supplementary_file[1]: GSE999001_RAW.tar" in pl["processing_text"]


def test_pin_form_is_checked():
    cfg = config()
    good = {"accession": "GSM1", "end": "late", "stage": "E14.5", "record_sha256": "a" * 64,
            "read_by": ["Reader A", "Reader B"], "date": "2026-11-21"}

    def check(samples, series=(), **kw):
        a = {**cm.atlas(ACC, cfg), **kw}
        a["pinned"] = {"samples": samples, "series": list(series)}
        return cm.check_pins(ACC, a)

    assert check([good])["samples"]["GSM1"]["stage"] == "E14.5"
    for bad in ({"read_by": ["READER_1", "READER_2"]}, {"read_by": ["Reader A"]}, {"read_by": ["A", "a"]},
                {"date": "YYYY-MM-DD"}, {"stage": "E12.5"}, {"end": "middle"}, {"record_sha256": "A" * 64},
                {"accession": "GSE1"}):
        with pytest.raises(cm.Refused):
            check([{**good, **bad}])
    with pytest.raises(cm.Refused):  # two stages at one end
        check([good, {**good, "accession": "GSM2", "stage": "E16.5"}])
    with pytest.raises(cm.Refused):  # an extra key
        check([{**good, "note": "x"}])
    with pytest.raises(cm.Refused):  # a series that is not the atlas's
        check([good], [{"accession": "GSE1", "record_sha256": "a" * 64, "read_by": ["A", "B"], "date": "2026-11-21"}])
    with pytest.raises(cm.Refused):  # a held atlas holds no pins
        check([good], status="held")


def test_a_permitted_sample_is_not_fetched_without_a_pin(tmp_path):
    """Independently of the rule: a sample it reads as permitted at a bracket stage is not USED, not fetched and not
    loaded until it is pinned; pinning only one end fetches nothing."""
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(BASE), []), config())
    assert all(r["mark"] != "USED" for r in pl["samples"])
    with pytest.raises(cm.Refused):
        fetch.fetch(ACC, tmp_path, cm.sha256_file(tmp_path / "meta" / ACC / "plan.json"),
                    getter_for(soft(BASE), []), config())
    # only the late end pinned (pasted by hand from a two-end template): refused at plan and at fetch
    cfg = pinned(tmp_path / "b", BASE)
    cfg["atlases"][ACC]["pinned"]["samples"] = cfg["atlases"][ACC]["pinned"]["samples"][1:]
    pl = fetch.plan(ACC, tmp_path / "b", getter_for(soft(BASE), []), cfg)
    assert pl["status"] == "refused" and any("do not cover both ends" in r for r in pl["refusals"])
    log = []
    with pytest.raises(cm.Refused):
        fetch.fetch(ACC, tmp_path / "b", cm.sha256_file(tmp_path / "b" / "meta" / ACC / "plan.json"),
                    getter_for(soft(BASE), log), cfg)
    assert not [u for u in log if not u.endswith("_family.soft.gz")]
    with pytest.raises(cm.Refused, match="not both"):  # the template itself refuses one end
        fetch.pin_template(ACC, tmp_path / "b", ["GSM9990002"], config())


def test_pinned_for_a_stage_its_record_does_not_state_is_refused(tmp_path):
    cfg = pinned(tmp_path, BASE)
    pins = cfg["atlases"][ACC]["pinned"]["samples"]
    pins[1]["stage"] = "E16.5"   # the record states E14.5
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(BASE), []), cfg)
    assert pl["status"] == "refused" and any("pinned at E16.5, but its record states E14.5" in r
                                             for r in pl["refusals"])
    assert {r["gsm"]: r["mark"] for r in pl["samples"]}["GSM9990002"] == "NOT_USED"
    pins[1]["stage"], pins[1]["end"] = "E8.5", "early"   # pinned at the other end, at a stage it does not state
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(BASE), []), cfg)
    assert pl["status"] == "refused" and any("pinned for the early end" in r for r in pl["refusals"])


def test_a_pinned_sample_the_rule_refuses_refuses_the_atlas(tmp_path):
    """The rule is a second line that can only refuse: a pinned sample whose live record gains a window reading is
    refused even though it is pinned (its hash changes too, but the refusal names the rule's mark)."""
    cfg = pinned(tmp_path, BASE)
    moved = copy.deepcopy(BASE)
    moved[1] = sample("GSM9990002", "E14.5 heart rep1", "E14.5", "heart", tag="E14_5",
                      extra=["!Sample_extract_protocol_ch1 = Hearts were dissected from E12.5 embryos."])
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(moved), []), cfg)
    assert pl["status"] == "refused" and any("pinned sample GSM9990002 is not usable" in r for r in pl["refusals"])


def test_unstated_convention_drops_the_e14_5_end(tmp_path):
    s = [*BASE, sample("GSM9990007", "E16.5 heart", "E16.5", "heart", tag="E16_5")]
    design = "One library per sample."
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(s, design=design), []), config())
    rows = {r["gsm"]: r for r in pl["samples"]}
    assert rows["GSM9990002"]["mark"] == "WINDOW"
    with pytest.raises(cm.Refused):  # a window sample cannot be pinned
        fetch.pin_template(ACC, tmp_path, ["GSM9990001", "GSM9990002"], config())
    pl, _ = plan_pinned(tmp_path / "b", s, ("GSM9990001", "GSM9990007"), design=design)
    assert pl["status"] == "passed" and pl["chosen"] == {"early": "E8.5", "late": "E16.5"}


@pytest.mark.parametrize("extra", [
    ["!Sample_characteristics_ch1 = gestation: E12.5"],
    ["!Sample_description = E12.5 hearts, library 2"],
    ["!Sample_characteristics_ch1 = tissue: hearts pooled from E12.5 and E14.5 embryos"],
    ["!Sample_characteristics_ch1 = tissue: hearts 13.5-14.5 dpc"],
    ["!Sample_description = hearts from E14.5 and 13.5 embryos"],
    ["!Sample_description = ED12.5 heart"],
    ["!Sample_extract_protocol_ch1 = Hearts were dissected from E12.5 embryos."],
    ["!Sample_growth_protocol_ch1 = Embryos were harvested at 12.5 dpc."],
    ["!Sample_data_processing = E12.5 and E14.5 libraries were processed with Cell Ranger count"],
])
def test_reviewers_records_never_become_used(tmp_path, extra):
    s = [BASE[0], sample("GSM9990002", "E14.5 heart rep2", "E14.5", "heart", tag="E14_5", extra=extra)]
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(s), []), config())
    row = next(r for r in pl["samples"] if r["gsm"] == "GSM9990002")
    assert row["mark"] in ("UNKNOWN", "WINDOW") and pl["status"] == "refused"
    if row["mark"] == "UNKNOWN":
        assert any("UNKNOWN" in r for r in pl["refusals"])
    assert "record" not in row
    with pytest.raises(cm.Refused):  # and it cannot be pinned
        fetch.pin_template(ACC, tmp_path, list(PAIR), config())


def test_plug_e1_in_the_extract_protocol(tmp_path):
    s = [BASE[0], sample("GSM9990002", "E14.5 heart", "E14.5", "heart", tag="E14_5",
                         extra=["!Sample_growth_protocol_ch1 = The day of plug was designated as E1."])]
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(s, design="One library per sample."), []), config())
    row = next(r for r in pl["samples"] if r["gsm"] == "GSM9990002")
    assert row["mark"] == "WINDOW" and pl["status"] == "refused"


def test_unknown_samples_refuse_only_when_the_atlas_would_use_them(tmp_path):
    other_assay = sample("GSM9990008", "heart ATAC", None, "heart", strategy="ATAC-seq")
    other_tissue = sample("GSM9990009", "liver", None, "liver")
    pl, _ = plan_pinned(tmp_path, [*BASE, other_assay, other_tissue])
    assert pl["status"] == "passed", pl["refusals"]
    rows = {r["gsm"]: r for r in pl["samples"]}
    assert rows["GSM9990008"]["mark"] == "UNKNOWN" and not rows["GSM9990008"]["would_use"]
    assert rows["GSM9990009"]["mark"] == "UNKNOWN" and "tissue" in rows["GSM9990009"]["would_use_reason"]
    s = [*BASE, sample("GSM9990006", "FetalHeart_1", None, "heart")]
    pl2, _ = plan_pinned(tmp_path / "b", s)
    assert pl2["status"] == "refused" and any("UNKNOWN" in r for r in pl2["refusals"])


@pytest.mark.parametrize("early,late,samples", [
    ("E9.5", "E14.5", [sample("GSM9990011", "E9.5 embryo", "E9.5", "whole embryo", tag="E9_5"), BASE[1]]),
    ("E8.5", "E16.5", [BASE[0], sample("GSM9990012", "E16.5 heart", "E16.5", "heart", tag="E16_5")]),
    ("E9.5", "P3", [sample("GSM9990011", "E9.5 embryo", "E9.5", "whole embryo", tag="E9_5"), BASE[4]]),
])
def test_clean_permitted_records_are_used_once_pinned(tmp_path, early, late, samples):
    """E9.5 exactly with E0.5 stated, E14.5, E16.5 and P3: USED once pinned, fetched after a live re-plan."""
    cfg = config(early=(early,), late=(late,))
    gsms = [samples[0][0].split(" = ")[1], samples[1][0].split(" = ")[1]]
    pl, cfg2 = plan_pinned(tmp_path, samples, gsms, cfg)
    assert pl["status"] == "passed", pl["refusals"]
    assert pl["chosen"] == {"early": early, "late": late}
    assert sorted(r["gsm"] for r in pl["samples"] if r["mark"] == "USED") == sorted(gsms)
    man = fetch.fetch(ACC, tmp_path, cm.sha256_file(tmp_path / "meta" / ACC / "plan.json"),
                      getter_for(soft(samples), []), cfg2)
    assert sorted(man["samples"]) == sorted(gsms)


def test_fallback_late_end_when_the_first_is_absent(tmp_path):
    s = [BASE[0], sample("GSM9990007", "E16.5 heart", "E16.5", "heart", tag="E16_5"), BASE[2]]
    pl, _ = plan_pinned(tmp_path, s, ("GSM9990001", "GSM9990007"))
    assert pl["status"] == "passed" and pl["chosen"]["late"] == "E16.5"


def test_processing_text_is_cleared_positively(tmp_path):
    """Blocker 4: a joint-processing sentence no term catches is printed, refused until cleared by its hash, and a
    changed text after clearance refuses again."""
    joint = ("Libraries from all stages were combined and normalized together in Seurat; cells passing joint QC "
             "were kept")
    s = copy.deepcopy(BASE)
    s[1] = sample("GSM9990002", "E14.5 heart rep1", "E14.5", "heart", tag="E14_5", processing=joint)
    pl, cfg = plan_pinned(tmp_path, s, clear=False)
    assert pl["status"] == "refused" and pl["refusals"] == [pl["refusals"][0]] and "not cleared" in pl["refusals"][0]
    assert joint in (tmp_path / "meta" / ACC / "processing_text.txt").read_text()
    assert joint in (tmp_path / "meta" / ACC / "processing_text_pinned.txt").read_text()
    pl2, cfg2 = plan_pinned(tmp_path / "b", s)
    assert pl2["status"] == "passed"
    s[1] = sample("GSM9990002", "E14.5 heart rep1", "E14.5", "heart", tag="E14_5", processing=joint + " (v2)")
    pl3 = fetch.plan(ACC, tmp_path / "c", getter_for(soft(s), []), cfg2)
    assert pl3["status"] == "refused" and any("now hashes" in r for r in pl3["refusals"])


def test_protocol_terms_and_gate_still_refuse_until_cleared(tmp_path):
    s = copy.deepcopy(BASE)
    s[1] = sample("GSM9990002", "E14.5 heart rep1", "E14.5", "heart", tag="E14_5",
                  processing="Libraries were combined with cellranger aggr")
    pl, cfg2 = plan_pinned(tmp_path, s)
    assert pl["status"] == "refused" and any("aggr" in r for r in pl["refusals"])
    cfg3 = copy.deepcopy(cfg2)
    cfg3["atlases"][ACC]["protocol_terms_cleared"] = {"aggr": "read: per-sample matrices"}
    pl = fetch.plan(ACC, tmp_path / "b", getter_for(soft(s), []), cfg3)
    assert pl["status"] == "passed", pl["refusals"]
    s2 = copy.deepcopy(BASE)
    s2[0] = [ln.replace("10x Genomics 3' v2 libraries", "MACS-enriched CD31 cells, 10x 3' v2") for ln in s2[0]]
    pl, _ = plan_pinned(tmp_path / "c", s2)
    assert pl["status"] == "refused" and any("enrichment" in r for r in pl["refusals"])


def test_fetch_needs_the_plan_sha_the_clearance_the_pins_and_an_unchanged_live_record(tmp_path):
    pl, cfg = plan_pinned(tmp_path, BASE)
    p = tmp_path / "meta" / ACC / "plan.json"
    log = []
    get = getter_for(soft(BASE), log)
    with pytest.raises(cm.Refused):
        fetch.fetch(ACC, tmp_path, "0" * 64, get, cfg)
    no_clear = copy.deepcopy(cfg)
    no_clear["atlases"][ACC].pop("processing_cleared")
    with pytest.raises(cm.Refused):  # the clearance removed after the plan
        fetch.fetch(ACC, tmp_path, cm.sha256_file(p), get, no_clear)
    no_pins = copy.deepcopy(cfg)
    no_pins["atlases"][ACC].pop("pinned")
    with pytest.raises(cm.Refused):  # the pins removed after the plan
        fetch.fetch(ACC, tmp_path, cm.sha256_file(p), get, no_pins)
    for change in (
            ("tissue: heart", "tissue: Heart"),                       # one character, no stage in it
            ("10x Genomics 3' v2 libraries", "10x Genomics 3' v2 libraries."),
            (PROC, PROC + "."),                                       # also changes the processing text
    ):
        changed = copy.deepcopy(BASE)
        changed[1] = [ln.replace(*change) for ln in changed[1]]
        log2 = []
        with pytest.raises(cm.Refused, match="no longer plans|differs"):
            fetch.fetch(ACC, tmp_path, cm.sha256_file(p), getter_for(soft(changed), log2), cfg)
        assert all(u.endswith("_family.soft.gz") for u in log2)  # no data file was requested
    series_changed = soft(BASE).replace("Hearts and embryos.", "Hearts and embryos!")
    with pytest.raises(cm.Refused):  # the series record changed
        fetch.fetch(ACC, tmp_path, cm.sha256_file(p), getter_for(series_changed, []), cfg)
    moved = copy.deepcopy(BASE)
    moved[1] = sample("GSM9990002", "E14.5 heart rep1", "E14.5", "heart", tag="E14_5",
                      extra=["!Sample_characteristics_ch1 = note: E12.5 littermates"])
    with pytest.raises(cm.Refused):  # the live record now holds a window reading
        fetch.fetch(ACC, tmp_path, cm.sha256_file(p), getter_for(soft(moved), []), cfg)
    log.clear()
    man = fetch.fetch(ACC, tmp_path, cm.sha256_file(p), get, cfg)
    assert sorted(man["samples"]) == list(PAIR)
    data = [u for u in log if not u.endswith("_family.soft.gz")]
    assert len(data) == 6 and all(("GSM9990001" in u or "GSM9990002" in u) for u in data)
    assert not any("RAW.tar" in u for u in log)
    m = json.loads((tmp_path / "raw" / ACC / "MANIFEST.json").read_text())
    s2 = m["samples"]["GSM9990002"]
    assert s2["record_sha256"] == cm.record_sha256("SAMPLE", "GSM9990002", s2["record"])
    assert s2["record_sha256"] == m["pins"]["samples"]["GSM9990002"]["record_sha256"]
    assert m["processing_sha256"] == pl["processing_sha256"] and s2["series_protocol"]


def test_pin_template_names_samples_and_never_bulk(tmp_path):
    fetch.plan(ACC, tmp_path, getter_for(soft(BASE), []), config())
    with pytest.raises(cm.Refused):
        fetch.pin_template(ACC, tmp_path, [], config())
    with pytest.raises(cm.Refused):  # a permitted sample that is not a candidate at a bracket end
        fetch.pin_template(ACC, tmp_path, ["GSM9990001", "GSM9990005"], config())
    with pytest.raises(cm.Refused):  # a window sample
        fetch.pin_template(ACC, tmp_path, ["GSM9990001", "GSM9990003"], config())
    t = fetch.pin_template(ACC, tmp_path, list(PAIR), config())
    assert "READER_1" in t and t.count("record_sha256") == 3 and "processing_cleared" in t
    with pytest.raises(cm.Refused):  # pasted without filling in the readers: the atlas check refuses it
        d = yaml.safe_load(textwrap.dedent("\n".join(ln for ln in t.splitlines() if not ln.startswith("#"))))
        cm.check_pins(ACC, {**cm.atlas(ACC, config()), **d})


def test_held_in_use_and_unlisted_accessions_are_refused(tmp_path):
    cfg = config({"status": "held"})
    get = getter_for(soft(BASE), [])
    pl = fetch.plan(ACC, tmp_path, get, cfg)
    assert pl["status"] == "refused"
    with pytest.raises(cm.Refused):
        fetch.fetch(ACC, tmp_path, cm.sha256_file(tmp_path / "meta" / ACC / "plan.json"), get, cfg)
    with pytest.raises(cm.Refused):
        fetch.pin_template(ACC, tmp_path, list(PAIR), cfg)
    with pytest.raises(cm.Refused):
        fetch.plan("GSE193346", tmp_path, get)
    with pytest.raises(cm.Refused):
        fetch.plan("GSE000000", tmp_path, get)
    assert json.loads((tmp_path / "meta" / ACC / "plan.json").read_text())["status"] == "refused"


def test_clearance_is_outside_the_entry_hash_and_pins_are_inside():
    cfg = config()
    c2 = copy.deepcopy(cfg)
    c2["atlases"][ACC]["processing_cleared"] = {"sha256": "1" * 64}
    assert cm.entry_sha256(ACC, cfg) == cm.entry_sha256(ACC, c2)
    c2["atlases"][ACC]["pinned"] = {"samples": [], "series": []}
    c3 = copy.deepcopy(c2)
    c3["atlases"][ACC]["pinned"]["samples"] = [{"accession": "GSM1", "end": "late", "stage": "E14.5",
                                                "record_sha256": "a" * 64, "read_by": ["A", "B"],
                                                "date": "2026-11-21"}]
    assert cm.entry_sha256(ACC, c2) != cm.entry_sha256(ACC, c3)


def test_committed_atlases_yaml_holds_no_pins_yet():
    cfg = cm.load_atlases()
    for acc, a in cfg["atlases"].items():
        assert not (a.get("pinned") or {}).get("samples"), acc


def test_plan_and_fetch_never_write_inside_the_repository():
    with pytest.raises(cm.Refused):
        fetch.plan(ACC, ROOT / "data" / "derived" / "t1_ext_test", getter_for(soft(BASE), []), config())
    with pytest.raises(cm.Refused):
        cm.dump_json({}, ROOT / "experiments" / "t1_ext" / "x.json")
