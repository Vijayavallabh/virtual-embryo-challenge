"""experiments/t1_ext in the strict form of candidate MX2s (ruling row R27), on synthetic records and synthetic files
only (no network, no released or outside data):

- the stage-parser exception, exactly as R27's third cell words it: a protocol or processing field byte-identical across
  every sample record of a series is read as series text (no stage reading), pooling words in it still refuse, and the
  sample's stage must then come from its own title and characteristics, which must agree; the parser stays the second
  line (it never makes an unpinned sample usable);
- the named samples and the never-fetch lists (the fetch takes per-sample folders only and refuses a series archive and
  every never-fetch accession by name);
- the SCREENS (a) to (d) on barcode and feature files before any matrix is fetched, the screens record, and the fetch
  of matrices refusing a member that is not recorded kept;
- the parked branch's known blocker, closed: a reduction is read only under the pins, clearance and screen in force;
- no code outside experiments/t1_ext (the shared core, the T2 heads, the builders such as build_test.py, T3) names it.
"""

from __future__ import annotations

import copy
import gzip
import importlib.util
import json
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
import yaml

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


fetch = _load("fetch")
ld = _load("load")
scr = _load("screens")
E05 = ["Timed matings: the day of the vaginal plug was designated E0.5."]
SHARED = "E8.5 and E10.5 embryos were dissected; hearts were collected at E12.5, E14.5, E16.5 and P3."


# ------------------------------------------------------------------------------------------------ the exception


def recs(n=3, shared=SHARED, own=None):
    """n sample records of one series, each with the same extract protocol; own: GSM -> fields to override."""
    out = {}
    stages = ["E8.5", "E16.5", "E12.5", "E16.5"]
    for i in range(n):
        g = f"GSM88{i:02d}"
        st = stages[i % len(stages)]
        out[g] = {"title": [f"{st} sample {i + 1}"], "characteristics_ch1": [f"developmental stage: {st}"],
                  "extract_protocol_ch1": [shared], "data_processing": [f"Cell Ranger count of library {i + 1}"]}
        out[g].update((own or {}).get(g, {}))
    return out


def test_the_exception_lets_through_a_sample_whose_only_obstacle_was_the_series_wide_sentence():
    rs = recs()
    st = cm.series_wide_fields(rs)
    assert st == {"extract_protocol_ch1": [SHARED]}  # data_processing differs per sample: not series text
    late = rs["GSM8801"]
    assert cm.record_stage(late, E05)["mark"] in ("WINDOW", "UNKNOWN")  # without the exception: refused
    r = cm.record_stage(late, E05, st)
    assert r["mark"] == "PERMITTED" and r["nominal_low"] == r["nominal_high"] == 16.5, r["reason"]
    assert r["series_text_fields"] == ["extract_protocol_ch1"]
    assert {x["field"] for x in r["readings"]} == {"title", "characteristics_ch1"}
    # the sample's own (non-series) processing text is still read: a stage there still refuses
    own = {**late, "data_processing": ["E12.5 libraries were processed with Cell Ranger count"]}
    assert cm.record_stage(own, E05, st)["mark"] in ("WINDOW", "UNKNOWN")
    # the series text still feeds the plug-day convention: a plug sentence in it with an unread day refuses
    rs2 = recs(shared="Females were checked for a plug on day 2 after mating; " + SHARED)
    st2 = cm.series_wide_fields(rs2)
    assert st2 and cm.record_stage(rs2["GSM8801"], [], st2)["mark"] == "UNKNOWN"


def test_a_field_identical_across_all_samples_but_one_is_not_series_text():
    rs = recs(4, own={"GSM8803": {"extract_protocol_ch1": [SHARED + " "]}})  # one character more in one sample
    # (parse_soft strips values, so a trailing space cannot come from a SOFT file; any other character counts too)
    rs["GSM8803"]["extract_protocol_ch1"] = [SHARED.replace("dissected", "dissected.")]
    assert cm.series_wide_fields(rs) == {}
    assert cm.record_stage(rs["GSM8801"], E05, cm.series_wide_fields(rs))["mark"] in ("WINDOW", "UNKNOWN")
    missing = recs(3)
    del missing["GSM8802"]["extract_protocol_ch1"]  # a field absent from one sample is not series-wide either
    assert cm.series_wide_fields(missing) == {}
    order = recs(3, own={"GSM8800": {"extract_protocol_ch1": ["Step two.", SHARED]},
                         "GSM8801": {"extract_protocol_ch1": [SHARED, "Step two."]},
                         "GSM8802": {"extract_protocol_ch1": [SHARED, "Step two."]}})
    assert cm.series_wide_fields(order) == {}  # the same values in another order are another field
    # a series text that is not this sample's own field gives it nothing
    other = {"extract_protocol_ch1": ["Another sentence."]}
    assert cm.record_stage(rs["GSM8801"], E05, other)["mark"] in ("WINDOW", "UNKNOWN")


ZW = "\u200b"  # a zero-width space (Unicode category Cf): NFKC keeps it


@pytest.mark.parametrize("sentence", ["Hearts were pooled from E8.5, E10.5 and E16.5 embryos.",
                                      "Libraries were POOLED and sequenced.",
                                      "Cells of all stages were merged before loading.",
                                      # the staged set: a pooling word in a sentence that names stages (the review of
                                      # 290da3f, its one blocker)
                                      "Embryos from E8.5 and E10.5 were combined before loading.",
                                      "Embryos from E8.5 and E10.5 were mixed before loading.",
                                      "Embryos from E8.5 and E10.5 were multiplexed in one lane.",
                                      "Embryos from E8.5 and E10.5 were superloaded together.",
                                      "Embryos from E8.5 and E10.5 were aggregated before loading.",
                                      "Libraries of E8.5 and E10.5 embryos were multiplexed and sequenced.",
                                      "Samples from all stages were mixed.",  # no stage label: the word stages
                                      "Cells from every time point were combined.",
                                      # a zero-width character inside the word
                                      f"Hearts were poo{ZW}led from E8.5, E10.5 and E16.5 embryos.",
                                      f"Libraries were poo{ZW}led and sequenced.",
                                      f"Embryos from E8.5 and E10.5 were com{ZW}bined before loading.",
                                      f"Hearts of all ages were me{ZW}rged.",
                                      # the review of 75356ec, blocker 1: noun forms and glued prefixes
                                      "A mixture of E8.5 and E10.5 embryos was loaded.",
                                      "A combination of E8.5 and E10.5 embryos was loaded.",
                                      "E8.5 and E10.5 embryos were loaded in combination.",
                                      "E8.5 and E10.5 embryos were admixed.",
                                      "E8.5 and E10.5 embryos were intermixed.",
                                      "E8.5 and E10.5 embryos were premixed.",
                                      "E8.5 and E10.5 embryos were recombined.",
                                      "E8.5 and E10.5 embryos were co‐mixed.",  # a Unicode hyphen
                                      # blocker 2: the hyphenated time point and the time course
                                      "Embryos at multiple time-points were mixed.",
                                      "Embryos across the timecourse were combined.",
                                      "Embryos across the time course were combined.",
                                      f"Embryos at multiple time{ZW}points were mixed.",
                                      # blocker 3: stages in one sentence, the pooling word in the next
                                      "Embryos were collected at E8.5 and E10.5; all were then mixed.",
                                      "Embryos were collected at E8.5 and E10.5. They were then combined and loaded.",
                                      "Embryos were collected at several stages. Cells were mixed and loaded."])
def test_pooling_words_in_series_text_still_refuse(sentence):
    rs = recs(shared=sentence)
    st = cm.series_wide_fields(rs)
    assert st == {"extract_protocol_ch1": [sentence]}
    r = cm.record_stage(rs["GSM8801"], E05, st)
    assert r["mark"] == "UNKNOWN" and "pooling word" in r["reason"], r["reason"]


@pytest.mark.parametrize("sentence", ["cDNA was combined with RT Master Mix.",
                                      "Cells were mixed with Trypan blue and counted.",
                                      "Libraries were prepared with the Single Cell Multiome kit.",
                                      f"cDNA was com{ZW}bined with RT Master Mix.",
                                      # two sentences, neither naming a stage: the whole value names none
                                      "Hearts were dissociated. cDNA was combined with RT Master Mix."])
def test_a_staged_pooling_word_in_a_value_naming_no_stage_does_not_refuse(sentence):
    rs = recs(shared=sentence)
    st = cm.series_wide_fields(rs)
    assert st == {"extract_protocol_ch1": [sentence]}
    r = cm.record_stage(rs["GSM8801"], E05, st)
    assert r["mark"] == "PERMITTED" and r["nominal_low"] == r["nominal_high"] == 16.5, r["reason"]


@pytest.mark.parametrize("sentence", [
    # benign prose that the whole-value rule now refuses (fails closed; listed in the README): a stage in one sentence
    # and reagent prose in another, a day count beside a mixing word, the word "stage" for a step of the protocol
    "Hearts were dissected at E16.5. cDNA was combined with RT Master Mix.",
    "cDNA was combined with RT Master Mix and incubated for 2 days.",
    "Samples were mixed by pipetting at each stage of the protocol."])
def test_benign_values_that_the_whole_value_rule_refuses(sentence):
    rs = recs(shared=sentence)
    st = cm.series_wide_fields(rs)
    assert st == {"extract_protocol_ch1": [sentence]}
    r = cm.record_stage(rs["GSM8801"], E05, st)
    assert r["mark"] == "UNKNOWN" and "pooling word" in r["reason"], r["reason"]


def test_a_sample_with_no_series_wide_field_reads_as_before():
    rs = recs(own={"GSM8800": {"extract_protocol_ch1": ["Hearts were dissected."]},
                   "GSM8801": {"extract_protocol_ch1": ["Hearts were dissected and dissociated."]}})
    st = cm.series_wide_fields(rs)
    assert st == {}
    r = cm.record_stage(rs["GSM8801"], E05, st)
    assert r["mark"] == "PERMITTED" and r["series_text_fields"] == [], r["reason"]


def test_a_cf_character_never_hides_a_pooling_word_that_matches_as_written():
    # as written, the zero-width space is a boundary and "pooled" matches; with it removed "xpooled" would not: a match
    # in either reading counts, so removing Cf characters never refuses less
    assert cm.pooling_in_series_text({"f": [f"x{ZW}pooled"]})
    assert cm.pooling_in_series_text({"f": [f"Embryos of E8.5 and E10.5 were x{ZW}mixed."]})
    assert not cm.pooling_in_series_text({"f": ["Hearts were dissected."]})


def test_the_pooling_patterns_are_in_the_entry_hash():
    src = (EXT / "t1x_common.py").read_text()
    body = src[src.index("def entry_sha256"):src.index("def require_frozen")]
    assert "POOLING_WORDS.pattern" in body and "POOLING_WORDS_STAGED.pattern" in body
    assert "POOLING_STAGE_WORDS.pattern" in body


@pytest.mark.parametrize("title,chars", [
    ("E16.5 heart", "developmental stage: E16"),   # [16.5] against [16, 16.5]: overlap, but not agreement
    ("heart rep1", "developmental stage: E16.5"),  # the title gives no stage
    ("E16.5 heart", "tissue: heart"),              # the characteristics give no stage
])
def test_title_and_characteristics_that_do_not_agree_refuse(title, chars):
    rs = recs(own={"GSM8801": {"title": [title], "characteristics_ch1": [chars]}})
    st = cm.series_wide_fields(rs)
    r = cm.record_stage(rs["GSM8801"], E05, st)
    assert r["mark"] == "UNKNOWN" and "title and characteristics" in r["reason"], r["reason"]
    # the requirement applies only when a field was read as series text: the same title and characteristics with no
    # series-wide field are read by the rule as before (here PERMITTED or UNKNOWN, never the exception's reason)
    plain = {"title": [title], "characteristics_ch1": [chars]}
    assert "title and characteristics" not in cm.record_stage(plain, E05)["reason"]


@pytest.mark.parametrize("title,chars,want", [
    ("E12.5 heart", "developmental stage: E12.5", "WINDOW"),
    ("E16.5 heart", "developmental stage: E10.5", ("WINDOW", "UNKNOWN")),
    ("E16.5 heart, 13 hearts", "developmental stage: E16.5", "UNKNOWN"),   # a companion integer
    ("E16.5 heart", "developmental stage: E16.5 (day 3 of culture)", ("WINDOW", "UNKNOWN")),  # a day count
    ("Heart TS25", "developmental stage: E16.5", "UNKNOWN"),               # Theiler 25 is not E16.5 exactly
])
def test_own_title_or_characteristics_in_the_window_or_unknown_win_whatever_the_protocol(title, chars, want):
    clean = "Hearts were dissected and dissociated with collagenase."
    for shared in (SHARED, clean):
        rs = recs(shared=shared, own={"GSM8801": {"title": [title], "characteristics_ch1": [chars]}})
        st = cm.series_wide_fields(rs)
        assert st == {"extract_protocol_ch1": [shared]}
        r = cm.record_stage(rs["GSM8801"], E05, st)
        assert r["mark"] in ((want,) if isinstance(want, str) else want), (title, chars, shared, r["reason"])


def test_a_series_with_one_sample_has_no_series_wide_field():
    one = recs(1)
    assert cm.series_wide_fields(one) == {}
    assert cm.series_wide_fields({}) == {}
    r = cm.record_stage(one["GSM8800"], E05, cm.series_wide_fields(one))
    assert r["series_text_fields"] == [] and r["mark"] in ("WINDOW", "UNKNOWN")


def test_only_protocol_and_processing_fields_can_be_series_text():
    rs = recs(own={g: {"description": ["E10.5 littermates"], "source_name_ch1": ["heart"]}
                   for g in ("GSM8800", "GSM8801", "GSM8802")})
    st = cm.series_wide_fields(rs)
    assert set(st) == {"extract_protocol_ch1"}  # description and source name are read as the sample's own
    assert cm.record_stage(rs["GSM8801"], E05, st)["mark"] in ("WINDOW", "UNKNOWN")


# ------------------------------------------------------------------------------------------------ plan, pins, fetch

ACC = "GSE999011"
EARLY, LATE, WIN, ATAC = "GSM9991101", "GSM9991102", "GSM9991103", "GSM9991104"


def url(gsm, name):
    return f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{gsm}_{name}"


def sample(gsm, title, stage, tissue, strategy="RNA-Seq", proto=SHARED, processing="Cell Ranger count per sample"):
    lines = [f"^SAMPLE = {gsm}", f"!Sample_title = {title}", f"!Sample_source_name_ch1 = {tissue}",
             f"!Sample_characteristics_ch1 = developmental stage: {stage}",
             f"!Sample_characteristics_ch1 = tissue: {tissue}", f"!Sample_library_strategy = {strategy}",
             "!Sample_molecule_ch1 = nuclear RNA", f"!Sample_extract_protocol_ch1 = {proto}",
             f"!Sample_data_processing = {processing}"]
    for i, f in enumerate(("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"), 1):
        lines.append(f"!Sample_supplementary_file_{i} = {url(gsm, f)}")
    return lines


def soft(samples):
    head = [f"^SERIES = {ACC}", "!Series_title = synthetic", "!Series_summary = Hearts.",
            "!Series_overall_design = One library per sample. The day of the vaginal plug was designated E0.5.",
            f"!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE999nnn/{ACC}/suppl/{ACC}_RAW.tar"]
    return "\n".join(head + [ln for s in samples for ln in s]) + "\n"


SAMPLES = [sample(EARLY, "E8.5-scRNA-Library", "E8.5", "heart"), sample(LATE, "E16.5-scRNA-Library", "E16.5", "heart"),
           sample(WIN, "E10.5-scRNA-Library", "E10.5", "heart"),
           sample(ATAC, "E16.5-ATAC-Library", "E16.5", "heart", strategy="ATAC-seq")]


def config(**extra):
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["atlases"][ACC] = {"status": "use", "series": [ACC], "form": "tenx_mtx", "multiome": True,
                           "early": {"labels": ["E8.5"], "tissue": "heart", "gsms": [EARLY]},
                           "late": {"labels": ["E16.5"], "tissue": "heart", "gsms": [LATE]},
                           "never_fetch": {"samples": [WIN, ATAC], "files": [f"{ACC}_RAW.tar"]},
                           "screens": ["suffix", "own_peaks", "no_hashing"],
                           "assay_rule": {"library_strategy": "RNA-Seq"}, "early_tissue": "heart",
                           "groups_formed": ["CM"], "qc": {"min_umi": 1, "min_genes": 1, "max_mito": 1.0},
                           "protocol_terms_cleared": {}, "protocol_gate_cleared": {}, **extra}
    return cfg


BARCODES = "".join(f"AAACGT{i:04d}-1\n" for i in range(12))


def features(peaks):
    rows = [f"ENSMUSG{i:05d}\tGene{i}\tGene Expression\tchr1\t{i}\t{i + 1}" for i in range(20)]
    rows += [f"chr1:{a}-{a + 500}\tchr1:{a}-{a + 500}\tPeaks\tchr1\t{a}\t{a + 500}" for a in peaks]
    return "\n".join(rows) + "\n"


def getter_for(text, log, files=None):
    """Serves the family SOFT and, per file name, synthetic gzipped content (files: name suffix -> text)."""
    files = files or {}

    def get(u, dst):
        log.append(u)
        dst.parent.mkdir(parents=True, exist_ok=True)
        if u.endswith("_family.soft.gz"):
            dst.write_bytes(gzip.compress(text.encode(), mtime=0))
            return
        name = u.rsplit("/", 1)[-1]
        body = next((v for k, v in files.items() if name.endswith(k)), "x\n")
        dst.write_bytes(gzip.compress(body.encode(), mtime=0))  # fixed bytes: a refetch is the same file
    return get


def from_template(text):
    body = "\n".join(ln for ln in text.splitlines() if not ln.startswith("#"))
    body = body.replace("READER_1", "Reader A").replace("READER_2", "Reader B").replace("YYYY-MM-DD", "2026-11-24")
    return yaml.safe_load(textwrap.dedent(body))


def pinned(tmp_path, cfg=None, samples=SAMPLES):
    cfg = cfg or config()
    fetch.plan(ACC, tmp_path, getter_for(soft(samples), []), cfg)
    t = fetch.pin_template(ACC, tmp_path, [EARLY, LATE], cfg)
    cfg2 = copy.deepcopy(cfg)
    cfg2["atlases"][ACC].update(from_template(t))
    cm.check_pins(ACC, cm.atlas(ACC, cfg2))
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(samples), []), cfg2)
    return pl, cfg2


def test_plan_reads_the_series_wide_sentence_as_series_text_and_still_needs_pins(tmp_path):
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(SAMPLES), []), config())
    rows = {r["gsm"]: r for r in pl["samples"]}
    assert pl["series_text"][ACC] == {"extract_protocol_ch1": [SHARED],
                                      "data_processing": ["Cell Ranger count per sample"]}
    # the parser lets the two named samples through but they are NOT used: the pins are the first line
    assert rows[EARLY]["mark"] == rows[LATE]["mark"] == "NOT_USED" and "not pinned" in rows[LATE]["reason"]
    assert rows[WIN]["mark"] == "WINDOW" and pl["status"] == "refused"
    recs_txt = (tmp_path / "meta" / ACC / "records.txt").read_text()
    assert "SERIES-WIDE FIELDS of GSE999011" in recs_txt and "data_processing, extract_protocol_ch1" in recs_txt
    pl2, cfg2 = pinned(tmp_path / "b")
    assert pl2["status"] == "passed", pl2["refusals"]
    assert sorted(r["gsm"] for r in pl2["samples"] if r["mark"] == "USED") == [EARLY, LATE]
    text = pl2["processing_text"]
    assert f"series_text extract_protocol_ch1[1]: {SHARED}" in text  # the readers' clearance covers it
    # without the series-wide reading (one sample's sentence differs), the same pins are refused by the rule
    changed = copy.deepcopy(SAMPLES)
    changed[2] = sample(WIN, "E10.5-scRNA-Library", "E10.5", "heart", proto=SHARED.replace("P3", "P3 hearts"))
    pl3 = fetch.plan(ACC, tmp_path / "b", getter_for(soft(changed), []), cfg2)
    assert pl3["status"] == "refused" and any("not usable" in r for r in pl3["refusals"])


def test_the_plan_prints_the_sha256_of_every_protocol_field_of_every_sample_and_nothing_else_of_a_window_sample(
        tmp_path):
    own = copy.deepcopy(SAMPLES)
    own[2] = sample(WIN, "E10.5-scRNA-Library", "E10.5", "heart", processing="E10.5 library processed alone")
    pl = fetch.plan(ACC, tmp_path, getter_for(soft(own), []), config())
    fs = pl["protocol_field_sha256"][ACC]
    assert set(fs) == {EARLY, LATE, WIN, ATAC}  # every sample record of the series, the window sample included
    assert all(set(v) == {"extract_protocol_ch1", "data_processing"} for v in fs.values())
    assert fs[EARLY]["extract_protocol_ch1"] == fs[WIN]["extract_protocol_ch1"] == \
        cm.field_sha256("extract_protocol_ch1", [SHARED])
    assert fs[WIN]["data_processing"] != fs[EARLY]["data_processing"]  # the readers see the difference
    lines = (tmp_path / "meta" / ACC / "protocol_field_sha256.tsv").read_text().splitlines()
    assert lines == fetch.field_sha_lines(pl) and len(lines) == 8
    assert all(ln.split("\t")[0] == "field_sha256" and len(ln.split("\t")) == 5 for ln in lines)
    win_lines = [ln for ln in lines if f"\t{WIN}\t" in ln]
    assert all("E10.5" not in ln and "alone" not in ln for ln in win_lines)  # nothing else of a window sample
    # the field hash is the hash of the field's lines in the canonical record text
    rec = [["data_processing", ["a", "bé"]]]
    body = cm.record_text("SAMPLE", "GSM1", rec).split("\n", 1)[1]
    assert cm.field_sha256("data_processing", ["a", "bé"]) == cm.text_sha256(body)


def test_a_stage_only_in_the_shared_sentence_refuses():
    rs = recs(own={"GSM8801": {"title": ["heart sample 2"], "characteristics_ch1": ["tissue: heart"]}},
              shared="Hearts were collected at E16.5.")
    st = cm.series_wide_fields(rs)
    assert st == {"extract_protocol_ch1": ["Hearts were collected at E16.5."]}
    r = cm.record_stage(rs["GSM8801"], E05, st)
    assert r["mark"] == "UNKNOWN"  # no stage of its own: refused (the shared sentence gives none)


def test_named_samples_and_never_fetch_lists(tmp_path):
    cfg = config()
    a = cm.atlas(ACC, cfg)
    good = {"accession": LATE, "end": "late", "stage": "E16.5", "record_sha256": "a" * 64,
            "read_by": ["Reader A", "Reader B"], "date": "2026-11-24"}
    for bad in ({"accession": WIN}, {"accession": "GSM9991199"}, {"accession": EARLY}):
        with pytest.raises(cm.Refused):  # never fetched, not named, named at the other end
            cm.check_pins(ACC, {**a, "pinned": {"samples": [{**good, **bad}], "series": []}})
    assert cm.pins_complete(a, {"samples": {LATE: good}}) == [f"named samples not pinned: {[EARLY]} (the strict "
                                                                f"form uses every named sample or none)"]
    fetch.plan(ACC, tmp_path, getter_for(soft(SAMPLES), []), cfg)
    with pytest.raises(cm.Refused, match="exactly the samples"):
        fetch.pin_template(ACC, tmp_path, [LATE, EARLY, WIN], cfg)
    for u, why in ((f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE999nnn/{ACC}/suppl/{ACC}_RAW.tar", "never"),
                   (url(WIN, "matrix.mtx.gz"), WIN), (url(ATAC, "atac_fragments.tsv.gz"), ATAC)):
        assert why in (cm.never_fetch_refusal(a, url=u) or "")
    assert cm.never_fetch_refusal(a, url=url(LATE, "matrix.mtx.gz")) is None
    with pytest.raises(cm.Refused, match="never-fetch"):
        fetch._get_data(a, WIN, {"name": f"{WIN}_matrix.mtx.gz", "url": url(WIN, "matrix.mtx.gz")}, tmp_path,
                        getter_for("", []))
    with pytest.raises(cm.Refused):  # not the sample's own folder
        fetch._get_data(a, LATE, {"name": f"{ACC}_RAW.tar", "url": f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/"
                                                                   f"GSE999nnn/{ACC}/suppl/{ACC}_RAW.tar"},
                        tmp_path, getter_for("", []))
    clash = copy.deepcopy(a)
    clash["never_fetch"]["samples"].append(LATE)
    with pytest.raises(cm.Refused, match="both named and never fetched"):
        cm.check_never_fetch(ACC, clash)
    # a never-fetch sample the rule would let through is never a candidate
    pl = fetch.plan(ACC, tmp_path / "c", getter_for(soft(SAMPLES), []),
                    config(never_fetch={"samples": [EARLY], "files": []},
                           early={"labels": ["E8.5"], "tissue": "heart", "gsms": ["GSM9991109"]}))
    row = next(r for r in pl["samples"] if r["gsm"] == EARLY)
    assert row["mark"] == "NOT_USED" and "never-fetch" in row["reason"]


def test_committed_atlases_yaml_is_the_strict_form(tmp_path):
    cfg = cm.load_atlases()
    at = cfg["atlases"]
    assert cfg["tasks"] == ["T1"] and "composition" not in cfg
    # holds before and after the placement commit (pins, screens and the freeze change no assertion here)
    st = {k: a["status"] for k, a in at.items()}
    assert st["GSE193346"] == "in_use" and st["GSE176063"] == st["GSE216970"] == "held"
    assert st["GSE244783"] in ("use", "dropped") and st["GSE230531"] in ("use", "dropped") and len(st) == 5
    assert {k: at[k].get("register") for k in ("GSE193346", "GSE244783", "GSE230531")} == cm.REGISTER_ROWS
    g4, g3 = at["GSE244783"], at["GSE230531"]
    assert (g4["early"]["labels"], g4["early"]["gsms"]) == (["E8.5"], ["GSM7829699"])
    assert (g4["late"]["labels"], g4["late"]["gsms"]) == (["E16.5"], ["GSM7829703"])
    assert g4["groups_formed"] == ["CM", "Endocardium", "Endothelium"] and g4["multiome"]
    assert set(g4["never_fetch"]["samples"]) == {"GSM7829701", "GSM7829700", "GSM7829702", "GSM7829704"}
    assert g4["never_fetch"]["files"] == ["GSE244783_RAW.tar"]
    assert g4["qc"] == {"min_umi": 500, "min_genes": 300, "max_mito": 0.05}
    assert set(g4["screens"]) == {"suffix", "own_peaks", "no_hashing"}
    assert (g3["early"]["gsms"], g3["late"]["gsms"]) == (["GSM7226268", "GSM7226269"], ["GSM7226274", "GSM7226276"])
    assert (g3["early"]["labels"], g3["late"]["labels"]) == (["E8.5"], ["E16.5"]) and g3["groups_formed"] == ["CM"]
    assert set(g3["never_fetch"]["samples"]) == {"GSM7226279", "GSM7226280", "GSM7226270", "GSM7226271",
                                                 "GSM7226272", "GSM7226273"}
    assert g3["never_fetch"]["files"] == ["GSE230531_RAW.tar"]
    assert g3["qc"] == {"min_umi": 1000, "min_genes": 500, "max_mito": 0.20}
    assert set(g3["screens"]) == {"suffix", "no_hashing", "whitelist_sized"}
    for e in ("early", "late"):  # every bracket end outside the T1 window under any plug convention
        for a in (g3, g4):
            d = cm.label_day(a[e]["labels"][0])
            assert not cm.touches_window(d + cm.UNSTATED[0], d + cm.UNSTATED[1])
    # a placeholder pin is refused by the atlas check
    with pytest.raises(cm.Refused):
        cm.check_pins("GSE244783", {**cm.atlas("GSE244783", cfg), "pinned": {"samples": [
            {"accession": "GSM7829699", "end": "early", "stage": "E8.5", "record_sha256": "PENDING",
             "read_by": ["READER_1", "READER_2"], "date": "YYYY-MM-DD"}], "series": []}})
    # a fetchable 10x atlas without screens, named samples or a never-fetch list is refused by the loader of the yaml
    raw = yaml.safe_load((EXT / "atlases.yaml").read_text())
    pending = {"verdict": "PENDING"}

    def without(key):
        def f(a):
            a.pop(key)
        return f

    for i, change in enumerate((without("screens"), without("never_fetch"), without("screen"),
                                lambda a: a.update(register="X018"),
                                lambda a: a["never_fetch"]["samples"].remove("GSM7829700"),
                                lambda a: a["late"].update(gsms=["GSM7829703", "GSM7829705"]),
                                lambda a: a.update(screen={"verdict": "READERS_VERDICT"}))):
        bad = copy.deepcopy(raw)
        a4 = bad["atlases"]["GSE244783"]
        a4.update(status="use", screen=pending, pinned={"samples": [], "series": []})  # the draft form, then one change
        change(a4)
        p = tmp_path / f"atlases_{i}.yaml"
        p.write_text(yaml.safe_dump(bad))
        with pytest.raises(cm.Refused):
            cm.load_atlases(p)
    ok = copy.deepcopy(raw)
    ok["atlases"]["GSE244783"].update(status="use", screen=pending, pinned={"samples": [], "series": []})
    (tmp_path / "ok.yaml").write_text(yaml.safe_dump(ok))
    cm.load_atlases(tmp_path / "ok.yaml")


# ------------------------------------------------------------------------------------------------ the screens


FILES = {f"{EARLY}_barcodes.tsv.gz": BARCODES, f"{LATE}_barcodes.tsv.gz": BARCODES,
         f"{EARLY}_features.tsv.gz": features(range(0, 5000, 1000)),
         f"{LATE}_features.tsv.gz": features(range(300, 6000, 700))}
WANT = {ACC: {EARLY: "early", LATE: "late"}}


def meta_ready(tmp_path, files=None, cfg=None):
    """The plan (no pins yet), then the screens' fetch of the named samples' barcodes and features into tmp_path (the
    screens' folder of the test)."""
    cfg = cfg or config()
    fetch.plan(ACC, tmp_path, getter_for(soft(SAMPLES), []), cfg)
    log = []
    fm = fetch.fetch_meta(ACC, tmp_path, cm.sha256_file(tmp_path / "meta" / ACC / "plan.json"),
                          getter_for(soft(SAMPLES), log, files or FILES), cfg, screen_root=tmp_path, samples=WANT)
    return fm, cfg, log


def test_fetch_meta_takes_barcodes_and_features_only_and_the_screens_pass(tmp_path):
    fm, cfg, log = meta_ready(tmp_path)
    data = [u for u in log if not u.endswith("_family.soft.gz")]
    assert len(data) == 4 and not any("matrix" in u for u in data) and all(EARLY in u or LATE in u for u in data)
    assert all("/geo/samples/" in u and "/suppl/" in u for u in data)  # each from the sample's own folder
    assert set(fm["samples"]) == {EARLY, LATE} and fm["tasks"] == ["T1"]
    assert all(set(s["files"]) == {"barcodes", "features"} for s in fm["samples"].values())
    assert sorted(p.name for p in (tmp_path / "screen" / ACC).rglob("*.gz")) == sorted(FILES)
    res = scr.screen(ACC, tmp_path, cfg, screen_root=tmp_path)
    ch = res["checks"]
    assert ch["passed"] and {k: ch[k]["passed"] for k in ("suffix", "own_peaks", "no_hashing")} == {
        "suffix": True, "own_peaks": True, "no_hashing": True}
    assert ch["per_file"][LATE]["barcodes"] == {"wc_l": 12, "not_suffix": 0}
    assert ch["per_file"][LATE]["features"]["types"] == {"Gene Expression": 20, "Peaks": 9}
    pair = ch["own_peaks"]["pairs"][f"{EARLY}|{LATE}"]
    assert pair["different_sha256"] and pair["peaks_partly_overlap"] and pair["gene_expression_rows_same_in_order"]
    block = scr.paste_block(res)
    assert "verdict: READERS_VERDICT" in block and fm["samples"][LATE]["files"]["features"]["sha256"] in block
    filled = yaml.safe_load(textwrap.dedent(block.replace("READERS_VERDICT", "kept").replace(
        "READER_1", "Reader A").replace("READER_2", "Reader B").replace("YYYY-MM-DD", "2026-11-24")))
    cm.check_screen_block(ACC, {**cm.atlas(ACC, cfg), **filled})  # the block fills the yaml's form
    with pytest.raises(cm.Refused):  # the placeholder verdict is refused
        cm.check_screen_block(ACC, {**cm.atlas(ACC, cfg), **yaml.safe_load(textwrap.dedent(block))})


def test_fetch_meta_refuses_before_any_request(tmp_path):
    cfg = config()
    fetch.plan(ACC, tmp_path, getter_for(soft(SAMPLES), []), cfg)
    sha = cm.sha256_file(tmp_path / "meta" / ACC / "plan.json")
    log = []

    def run(**kw):
        args = {"dest": tmp_path, "plan_sha256": sha, "getter": getter_for(soft(SAMPLES), log, FILES), "cfg": cfg,
                "screen_root": tmp_path, "samples": WANT, **kw}
        return fetch.fetch_meta(ACC, **args)

    for kw, why in (({"plan_sha256": "0" * 64}, "plan sha256"),
                    ({"dest": tmp_path / "elsewhere"}, "write only under"),
                    ({"samples": {ACC: {EARLY: "early"}}}, "names"),
                    ({"cfg": config(qc={"min_umi": 2, "min_genes": 1, "max_mito": 1.0})}, "another atlas entry"),
                    ({"cfg": config(screens=[])}, "two 10x members")):
        with pytest.raises(cm.Refused, match=why):
            run(**kw)
    assert log == []  # nothing was requested, not even the SOFT
    # a named sample the plan reads in the window is not a candidate
    wcfg = config(late={"labels": ["E16.5"], "tissue": "heart", "gsms": [WIN]},
                  never_fetch={"samples": [ATAC], "files": []})
    fetch.plan(ACC, tmp_path / "w", getter_for(soft(SAMPLES), []), wcfg)
    with pytest.raises(cm.Refused, match="not a candidate"):
        fetch.fetch_meta(ACC, tmp_path / "w", cm.sha256_file(tmp_path / "w" / "meta" / ACC / "plan.json"),
                         getter_for(soft(SAMPLES), log, FILES), wcfg, screen_root=tmp_path / "w",
                         samples={ACC: {EARLY: "early", WIN: "late"}})
    assert log == []
    a = cm.atlas(ACC, cfg)
    for role, f, why in (("matrix", {"name": f"{LATE}_matrix.mtx.gz", "url": url(LATE, "matrix.mtx.gz")}, "role"),
                         ("features", {"name": f"{LATE}_matrix.mtx.gz", "url": url(LATE, "matrix.mtx.gz")}, "matrix"),
                         ("features", {"name": f"{ACC}_RAW.tar", "url": f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/"
                                                                      f"GSE999nnn/{ACC}/suppl/{ACC}_RAW.tar"}, None),
                         ("barcodes", {"name": f"{WIN}_barcodes.tsv.gz", "url": url(WIN, "barcodes.tsv.gz")}, None),
                         ("barcodes", {"name": "GSM7829701_barcodes.tsv.gz",
                                       "url": url("GSM7829701", "barcodes.tsv.gz")}, "never-fetched")):
        with pytest.raises(cm.Refused, match=why):
            fetch._requested(a, f["name"].split("_")[0], role, f, r"matrix\.mtx(\.gz)?$")
    # a changed live record refuses after the SOFT and before any data file
    changed = copy.deepcopy(SAMPLES)
    changed[1] = sample(LATE, "E16.5-scRNA-Library", "E16.5", "heart", processing="Cell Ranger count, library 2")
    log.clear()
    with pytest.raises(cm.Refused, match="live record"):
        run(getter=getter_for(soft(changed), log, FILES))
    assert all(u.endswith("_family.soft.gz") for u in log)


def test_the_screens_folder_and_the_literal_lists():
    assert cm.SCREEN_ROOT == Path.home() / "v" / "scratch" / "t72_screen"
    assert cm.SCREEN_SAMPLES == {"GSE244783": {"GSM7829699": "early", "GSM7829703": "late"},
                                 "GSE230531": {"GSM7226268": "early", "GSM7226269": "early", "GSM7226274": "late",
                                               "GSM7226276": "late"}}
    assert cm.NEVER_FETCHED_ALL == {"GSM7829701", "GSM7829700", "GSM7829702", "GSM7829704", "GSE244783_RAW.tar",
                                    "GSM7226279", "GSM7226280", "GSM7226270", "GSM7226271", "GSM7226272",
                                    "GSM7226273", "GSE230531_RAW.tar"}
    for n in cm.NEVER_FETCHED_ALL:
        assert cm.never_fetched_refusal(n) and cm.never_fetched_refusal(f"https://x/{n}_matrix.mtx.gz")
    assert cm.never_fetched_refusal("https://x/GSE1_raw.TAR")
    assert cm.never_fetched_refusal("GSM78297010") is None and cm.never_fetched_refusal("GSM7829699") is None


@pytest.mark.parametrize("change,check", [
    ({f"{LATE}_barcodes.tsv.gz": BARCODES.replace("0003-1", "0003-2")}, "suffix"),
    ({f"{LATE}_features.tsv.gz": features(range(0, 5000, 1000))}, "own_peaks"),          # one file: same sha256
    ({f"{LATE}_features.tsv.gz": features(range(0, 3000, 1000))}, "own_peaks"),          # a subset: no partial overlap
    ({f"{LATE}_features.tsv.gz": features(range(20000, 26000, 1000))}, "own_peaks"),     # disjoint peak sets
    ({f"{LATE}_features.tsv.gz": features([])}, "own_peaks"),                             # no Peaks rows
    ({f"{LATE}_features.tsv.gz": features(range(300, 6000, 700)).replace("Gene3\t", "GeneX\t")}, "own_peaks"),
    ({f"{LATE}_barcodes.tsv.gz": BARCODES.replace("-1\n", "-1\tE16.5\n")}, "no_hashing"),  # a per-cell column
    ({f"{LATE}_features.tsv.gz": features(range(300, 6000, 700)) +
      "HTO_A\tHTO_A\tAntibody Capture\tNA\t0\t0\n"}, "no_hashing"),
])
def test_a_failing_screen_says_no(tmp_path, change, check):
    _, cfg, _ = meta_ready(tmp_path, {**FILES, **change})
    ch = scr.screen(ACC, tmp_path, cfg, screen_root=tmp_path)["checks"]
    assert not ch["passed"] and not ch[check]["passed"]
    assert all(ch[k]["passed"] for k in ("suffix", "own_peaks", "no_hashing") if k != check), ch
    if check == "suffix":
        assert ch["suffix"]["barcodes_not_ending_in_-1"] == 1


def test_whitelist_sized_is_the_v2_whitelist(tmp_path):
    assert scr.WHITELIST_V2 == 737280
    p = tmp_path / "b.tsv.gz"
    p.write_bytes(gzip.compress("".join(f"AAAC{i:07d}-1\n" for i in range(737280)).encode()))
    assert scr.check_d({"g": scr.barcode_facts(p)})["passed"]
    p.write_bytes(gzip.compress("".join(f"AAAC{i:07d}-1\n" for i in range(737281)).encode()))
    assert not scr.check_d({"g": scr.barcode_facts(p)})["passed"]
    q = tmp_path / "c.tsv"
    q.write_text("A-1\nB-1")  # no newline at the end: wc -l says 1
    assert scr.newline_count(q) == 1 and scr.barcode_facts(q)["barcodes"] == 2


def test_partial_overlap_hand_cases():
    a = [("chr1", 100, 200)]
    assert not scr.partial_overlap(a, [("chr1", 100, 200)])     # identical
    assert scr.partial_overlap(a, [("chr1", 150, 250)])          # partial
    assert scr.partial_overlap(a, [("chr1", 120, 180)])          # contained, not identical
    assert not scr.partial_overlap(a, [("chr1", 200, 300)])      # touching, half-open
    assert not scr.partial_overlap(a, [("chr2", 150, 250)])      # another chromosome
    assert scr.partial_overlap(a + [("chr1", 5000, 5100)], [("chr1", 0, 99), ("chr1", 4000, 5050)])


def test_fetch_ties_the_fetched_files_to_the_screened_copies(tmp_path):
    fm, cfg, _ = meta_ready(tmp_path / "screen")
    files = {g: {r: {"name": f["name"], "sha256": f["sha256"]} for r, f in s["files"].items()}
             for g, s in fm["samples"].items()}
    kept = {"verdict": "kept", "files": files, "read_by": ["Reader A", "Reader B"], "date": "2026-11-24"}
    for scr_block, ok in (({"verdict": "PENDING"}, False), (kept, True)):
        pl, cfg2 = pinned(tmp_path / ("row" if ok else "pending"), config(screen=scr_block))
        sha = cm.sha256_file(tmp_path / ("row" if ok else "pending") / "meta" / ACC / "plan.json")
        log = []
        if not ok:
            with pytest.raises(cm.Refused, match="PENDING"):
                fetch.fetch(ACC, tmp_path / "pending", sha, getter_for(soft(SAMPLES), log, FILES), cfg2)
            assert not any("matrix" in u for u in log)
            continue
        man = fetch.fetch(ACC, tmp_path / "row", sha, getter_for(soft(SAMPLES), log, FILES), cfg2)
        assert man["screen"] == files
        assert man["samples"][LATE]["files"]["features"]["sha256"] == files[LATE]["features"]["sha256"]
        assert man["series_text"][ACC] and man["samples"][LATE]["series_text"] == man["series_text"][ACC]
    # a features file that is not the screened copy refuses the fetch
    bad = {**FILES, f"{LATE}_features.tsv.gz": features(range(300, 6000, 700)) + "\n"}
    with pytest.raises(cm.Refused, match="not the file the readers screened"):
        fetch.fetch(ACC, tmp_path / "row", sha, getter_for(soft(SAMPLES), [], bad), cfg2)
    # the yaml's screen block form
    a = cm.atlas(ACC, cfg2)
    for blk in ({**kept, "files": {LATE: files[LATE]}}, {**kept, "read_by": ["Reader A"]},
                {**kept, "files": {**files, LATE: {"barcodes": files[LATE]["barcodes"]}}},
                {"verdict": "dropped", "failed": ["b"]}, {"verdict": "maybe"},
                {"verdict": "PENDING", "files": files}):
        with pytest.raises(cm.Refused):
            cm.check_screen_block(ACC, {**a, "screen": blk})
    cm.check_screen_block(ACC, {**a, "status": "dropped", "screen": {"verdict": "dropped", "failed": ["b", "4"]}})
    with pytest.raises(cm.Refused):
        cm.check_screen_block(ACC, {**a, "status": "dropped", "screen": {"verdict": "dropped", "failed": ["e"]}})
    assert cm.screen_kept(ACC, {**a, "screen": kept}) == files
    assert cm.screen_kept("X", {"screens": []}) is None


# ------------------------------------------------------------------------------------------------ the known blocker

LACC = "GSE999012"


def lconfig():
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["atlases"][LACC] = {"status": "use", "series": [LACC], "form": "tenx_mtx",
                            "early": {"labels": ["E8.5"], "tissue": "embryo"},
                            "late": {"labels": ["E16.5"], "tissue": "heart"},
                            "assay_rule": {"library_strategy": "RNA-Seq"}, "early_tissue": "whole_embryo",
                            "groups_formed": ["CM"], "qc": {"min_umi": 1, "min_genes": 1, "max_mito": 1.0},
                            "protocol_terms_cleared": {}, "protocol_gate_cleared": {}}
    return cfg


SERIES = {LACC: [["title", ["synthetic"]], ["summary", ["Hearts."]], ["overall_design", E05]]}


def reduction(tmp_path):
    """A loaded and written reduction of a synthetic atlas, made under pins and a clearance (load.py end to end, with a
    synthetic reader instead of files)."""
    from scipy import sparse

    samples = {}
    for g, title, label, end in (("GSM31", "E8.5 embryo", "E8.5", "early"), ("GSM32", "E16.5 heart", "E16.5", "late")):
        rec = [["title", [title]], ["data_processing", ["Cell Ranger count per sample"]]]
        samples[g] = {"series": LACC, "end": end, "label": label, "record": rec,
                      "record_sha256": cm.record_sha256("SAMPLE", g, rec), "series_protocol": E05, "files": {}}
    cfg = lconfig()
    who = {"read_by": ["Reader A", "Reader B"], "date": "2026-11-24"}
    cfg["atlases"][LACC]["pinned"] = {
        "samples": [{"accession": g, "end": s["end"], "stage": s["label"], "record_sha256": s["record_sha256"], **who}
                    for g, s in samples.items()],
        "series": [{"accession": LACC, "record_sha256": cm.record_sha256("SERIES", LACC, SERIES[LACC]), **who}]}
    text = cm.processing_text(LACC, [(LACC, SERIES[LACC])],
                              {g: {"series": LACC, "text": cm.sample_text(s["record"])} for g, s in samples.items()})
    cfg["atlases"][LACC]["processing_cleared"] = {"sha256": cm.text_sha256(text), "read_by": "test"}
    rdir = tmp_path / "raw" / LACC
    rdir.mkdir(parents=True)
    (rdir / "MANIFEST.json").write_text(json.dumps({
        "accession": LACC, "atlas_entry_sha256": cm.entry_sha256(LACC, cfg), "samples": samples, "tasks": ["T1"],
        "series_records": SERIES, "chosen": {"early": "E8.5", "late": "E16.5"}}))
    rng = np.random.default_rng(0)

    def reader(acc, form, a, s, folder, rd):
        X = sparse.csr_matrix(rng.poisson(3.0, size=(6, 8)).astype(np.int64))
        return {"C": X, "barcodes": np.array([f"B{i}-1" for i in range(6)]),
                "genes": np.array([f"G{i}" for i in range(8)]), "gates": {}}

    red = ld.load_atlas(LACC, tmp_path, cfg, check_frozen=False, reader=reader)
    ld.write_reduced(LACC, tmp_path, red, cfg)
    return cfg


def test_the_tables_step_reads_no_reduction_made_under_other_pins(tmp_path):
    """The parked branch's known blocker: tables.py read an earlier reduction without a pin check. read_reduced now
    refuses a reduction unless the pins, the atlas entry, the clearance and the screen in force are the ones it was
    made under."""
    cfg = reduction(tmp_path)
    assert ld.read_reduced(LACC, tmp_path, cfg)["C"].shape[0] == 12
    with pytest.raises(cm.Refused):  # the committed yaml does not hold this atlas: no reduction is read blind
        ld.read_reduced(LACC, tmp_path)
    repinned = copy.deepcopy(cfg)
    repinned["atlases"][LACC]["pinned"]["samples"][1]["record_sha256"] = "b" * 64  # a reader re-pinned another record
    with pytest.raises(cm.Refused, match="another atlas entry"):
        ld.read_reduced(LACC, tmp_path, repinned)
    unpinned = copy.deepcopy(cfg)
    unpinned["atlases"][LACC]["pinned"] = {"samples": [], "series": []}
    with pytest.raises(cm.Refused):
        ld.read_reduced(LACC, tmp_path, unpinned)
    held = copy.deepcopy(cfg)
    held["atlases"][LACC]["status"] = "held"
    with pytest.raises(cm.Refused, match="not read"):
        ld.read_reduced(LACC, tmp_path, held)
    lapsed = copy.deepcopy(cfg)
    lapsed["atlases"][LACC]["processing_cleared"] = {"sha256": "c" * 64}
    with pytest.raises(cm.Refused, match="now hashes"):
        ld.read_reduced(LACC, tmp_path, lapsed)
    # the same entry, but the manifest's records no longer match the pins (hand-edited after the load)
    out = tmp_path / "reduced" / LACC
    man = json.loads((out / "MANIFEST.json").read_text())
    man["stage_records"]["GSM32"]["record"][0][1] = ["E16.5 hearts"]
    (out / "MANIFEST.json").write_text(json.dumps(man))
    with pytest.raises(cm.Refused, match="not its pin"):
        ld.read_reduced(LACC, tmp_path, cfg)


def test_a_screened_member_dropped_after_its_load_is_not_read(tmp_path, monkeypatch):
    cfg = reduction(tmp_path)
    scfg = copy.deepcopy(cfg)
    scfg["atlases"][LACC]["screens"] = ["suffix"]
    files = {g: {r: {"name": f"{g}_{r}.tsv.gz", "sha256": "a" * 64} for r in cm.SCREEN_ROLES}
             for g in ("GSM31", "GSM32")}
    man = {"screen": copy.deepcopy(files), "samples": {g: {"files": copy.deepcopy(f)} for g, f in files.items()}}
    with pytest.raises(cm.Refused, match="None"):  # no screen record for it
        ld.check_screen(LACC, cm.atlas(LACC, scfg), man)
    scfg["atlases"][LACC]["screen"] = {"verdict": "kept", "files": files}
    assert ld.check_screen(LACC, cm.atlas(LACC, scfg), man) == files
    other = copy.deepcopy(man)
    other["samples"]["GSM32"]["files"]["features"]["sha256"] = "b" * 64  # fetched, but not the screened copy
    with pytest.raises(cm.Refused, match="not the file the readers screened"):
        ld.check_screen(LACC, cm.atlas(LACC, scfg), other)
    scfg["atlases"][LACC]["screen"] = {"verdict": "dropped", "failed": ["a"]}
    with pytest.raises(cm.Refused, match="dropped"):
        ld.check_screen(LACC, cm.atlas(LACC, scfg), man)


def test_tables_reads_each_reduction_with_the_atlas_config(monkeypatch, tmp_path):
    tb = _load("tables")
    seen = []

    def stop(acc, dest, cfg=None):
        seen.append(cfg)
        raise cm.Refused("stop")

    monkeypatch.setattr(ld, "read_reduced", stop)
    monkeypatch.setattr(cm, "require_frozen", lambda c: None)
    gr = _load("groups")
    monkeypatch.setattr(gr, "load_extra", lambda c: {})
    monkeypatch.setattr(tb, "t1_genes_from", lambda p: ["G0"])
    cfg = cm.load_atlases()
    with pytest.raises(cm.Refused, match="stop"):
        tb.run(cfg, tmp_path, ["GSE244783"], tmp_path, "x")
    assert seen == [cfg]


# ------------------------------------------------------------------------------------------------ Task 1 only


NAMES_T1_EXT = re.compile(r"experiments[/.]t1_ext\b|(?<![\w])t1_ext[/.]|[\"']t1_ext[\"']|t1x_common")


def test_no_code_outside_t1_ext_names_it():
    """E8.5 lies inside the T2 heart window: no T2 path, the shared core, the builders (build_test.py among them) or
    T3 code may import the outside-atlas tables."""
    for ref in ('EXT = ROOT / "experiments" / "t1_ext"', "experiments/t1_ext/tables.py", "import t1x_common",
                "from experiments.t1_ext import tables", "sys.path.insert(0, str(HERE / 't1_ext'))"):
        assert NAMES_T1_EXT.search(ref), ref
    assert not NAMES_T1_EXT.search("tests/test_t1_ext_mx2s.py tests/test_t1_ext_load.py")
    roots = [ROOT / "src", ROOT / "experiments"]
    hits = []
    for r in roots:
        for p in [*r.rglob("*.py"), *r.rglob("*.sh")]:
            if EXT in p.parents or "results" in p.relative_to(ROOT).parts:
                continue  # results folders hold outputs, some never opened (the brief's never-open list)
            t = p.read_text(errors="replace")
            if NAMES_T1_EXT.search(t):  # a path into the folder or its module; a test file's name is not a use
                hits.append(str(p.relative_to(ROOT)))
    assert hits == []
    assert "t1_ext" not in (ROOT / "experiments" / "builds" / "build_test.py").read_text()
