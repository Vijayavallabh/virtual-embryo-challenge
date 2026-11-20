"""experiments/t3_ext: the per-dataset `stage_fields` key (where a sample's stage is read from), on synthetic records
only. Nothing is downloaded."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t3x_common as cm  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(f"t3_ext_sf_{name}", EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


fe = _load("fetch")

ACC = "GSE900001"
PIN = "    pinned: []\n"
SF = "    stage_fields: [characteristics_ch1, source_name_ch1, description]\n"
DS = """version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
  - accession: GSE900001
    status: cleared
    tier: unrelated
    study: GSE900001
    design: germline
    control_relation: same_series
    repository: geo
    screen_series: [GSE900001]
    stages_used: [E8.5]
    arms:
      control: {on: text, all: ['(?<![a-z])wt(?![a-z])']}
      ko:
        - {gene: Hoxz, on: text, all: ['hoxz null']}
    sample_file_patterns: ['(matrix\\.mtx|barcodes\\.tsv|features\\.tsv)(\\.gz)?$']
""" + PIN
FILES = ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")


def _soft(samples: dict[str, dict]) -> str:
    lines = [f"^SERIES = {ACC}", "!Series_title = A series", "!Series_overall_design = embryos",
             "!Series_summary = a summary", "!Series_contact_name = A,,Person"]
    for gsm, f in samples.items():
        lines.append(f"^SAMPLE = {gsm}")
        for k, vals in {**f, "series_id": ACC}.items():
            for v in (vals if isinstance(vals, list) else [vals]):
                lines.append(f"!Sample_{k} = {v}")
    return "\n".join(lines) + "\n"


def _rec(gsm, source="whole embryo", desc=None, chars=("genotype: WT",), title=None):
    sup = [f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/{gsm}/suppl/{gsm}_{f}" for f in FILES]
    r = {"title": title or gsm, "source_name_ch1": source, "characteristics_ch1": list(chars)}
    if desc is not None:
        r["description"] = desc
    return {**r, **{f"supplementary_file_{i + 1}": s for i, s in enumerate(sup)}}


# GSE204906's form: the stage only in the source name and the description
SAMPLES = {
    "GSM1": _rec("GSM1", "E8.5 embryo", "Hoxz null E8.5 embryos", ("genotype: Hoxz null",)),
    "GSM2": _rec("GSM2", "E8.5 embryo", "WT E8.5 embryos"),
    "GSM3": _rec("GSM3", "E8.5 embryo"),                                    # source name only
    "GSM4": _rec("GSM4", "whole embryo", "WT E8.5 embryos"),               # description only
    "GSM5": _rec("GSM5", "E8.5 embryo", "WT E9.5 embryos"),                # the two disagree
    "GSM6": _rec("GSM6", "whole embryo", "WT embryos", title="WT E8.5"),   # none listed states one (title unlisted)
    "GSM7": _rec("GSM7", "E10.5 embryo", "WT E10.5 embryos"),              # after the cap
    "GSM8": _rec("GSM8", "E8.5 embryo", "WT E8.5 embryos", ("genotype: WT", "developmental stage: E8.5")),
    "GSM9": _rec("GSM9", "E8.5 embryo", "WT E8.5 embryos", ("genotype: WT", "developmental stage: E9.5")),
    "GSM10": _rec("GSM10", "E8.5-E9.5 embryos", "WT E8.5 embryos"),        # a range
}


def _ds_file(tmp_path, stage_fields=SF, name="datasets.yaml"):
    p = tmp_path / name
    p.write_text(DS.replace(PIN, stage_fields + PIN))
    return p


def _plan(dsy, dest, samples=SAMPLES):
    return fe.plan(ACC, dest, dsy, metas={ACC: fe.parse_soft(_soft(samples))})


def _rows(pl):
    return {r["accession"]: r for r in pl["tables"][ACC]}


def _table(dest):
    return (dest / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()


def test_the_stage_is_read_from_the_source_name_and_the_description(tmp_path):
    pl = _plan(_ds_file(tmp_path), tmp_path / "d")
    assert pl["status"] == "ok", pl.get("reason")
    rows = _rows(pl)
    assert (rows["GSM1"]["stage"], rows["GSM1"]["stage_from"]) == (8.5, "read from source_name_ch1, description")
    assert (rows["GSM3"]["stage"], rows["GSM3"]["stage_from"]) == (8.5, "read from source_name_ch1")
    assert (rows["GSM4"]["stage"], rows["GSM4"]["stage_from"]) == (8.5, "read from description")
    assert (rows["GSM8"]["stage"], rows["GSM8"]["stage_from"]) == (
        8.5, "read from characteristics_ch1, source_name_ch1, description")
    txt = _table(tmp_path / "d")
    assert "parsed stage E8.5 (read from source_name_ch1, description)" in txt


@pytest.mark.parametrize("gsm,why", [
    ("GSM5", "listed stage fields disagree: source_name_ch1 gives E8.5; description gives E9.5"),
    ("GSM9", "listed stage fields disagree: characteristics_ch1 gives E9.5; source_name_ch1 gives E8.5; "
             "description gives E8.5"),
    ("GSM6", "no listed stage field states a stage (characteristics_ch1, source_name_ch1, description)"),
    ("GSM10", "range in source_name_ch1 'E8.5-E9.5 embryos'"),
])
def test_disagreement_absence_or_a_range_give_none_and_a_refusal_to_pin(tmp_path, gsm, why):
    dsy = _ds_file(tmp_path)
    pl = _plan(dsy, tmp_path / "d")
    r = _rows(pl)[gsm]
    assert r["stage"] is None and r["stage_from"] == why
    assert f"parsed stage none ({why})" in _table(tmp_path / "d")
    with pytest.raises(cm.Refused, match="no single stage parses"):
        fe.pin_template(ACC, tmp_path / "d", "control", "Hoxz", [gsm], datasets=dsy)


def test_a_stage_after_e9_5_read_from_a_listed_field_is_still_capped(tmp_path):
    dsy = _ds_file(tmp_path)
    r = _rows(_plan(dsy, tmp_path / "d"))["GSM7"]
    assert r["stage"] is None and "after E9.5" in r["stage_from"] and "read from source_name_ch1, description" in \
        r["stage_from"]
    with pytest.raises(cm.Refused, match="no single stage parses"):
        fe.pin_template(ACC, tmp_path / "d", "control", "Hoxz", ["GSM7"], datasets=dsy)
    assert fe.sample_stage({"source_name_ch1": ["E9.5 embryo"]}, ("source_name_ch1",)) == (9.5, "read from "
                                                                                                 "source_name_ch1")
    assert fe.sample_stage({"source_name_ch1": ["E9.75 embryo"]}, ("source_name_ch1",))[0] is None


def test_a_listed_field_absent_from_the_record_states_no_stage():
    rec = {"source_name_ch1": ["E8.5 embryo"]}
    assert fe.sample_stage(rec, ("characteristics_ch1", "source_name_ch1", "description", "no_such_field")) == \
        (8.5, "read from source_name_ch1")


def test_the_default_dataset_is_unchanged(tmp_path):
    """No stage_fields (or exactly [characteristics_ch1]) is the rule from before: the characteristics, else the
    title, with the same 'how' text; the plan rows and table text are identical either way."""
    rec_ch = {"title": ["x"], "characteristics_ch1": ["genotype: wt", "developmental stage: E8.5"]}
    rec_title = {"title": ["WT E8.5 rep1"], "characteristics_ch1": ["genotype: wt"],
                 "source_name_ch1": ["E9.5 embryo"]}
    assert fe.sample_stage(rec_ch) == (8.5, "characteristics")
    assert fe.sample_stage(rec_title) == (8.5, "characteristics")         # the old 'how' quirk is kept as it was
    assert fe.sample_stage({"title": ["WT E8.5"]}) == (8.5, "title")
    assert fe.sample_stage({"title": ["WT"], "source_name_ch1": ["E8.5 embryo"]}) == (None, "0 stage values")
    assert fe.sample_stage(rec_title, ("characteristics_ch1",)) == fe.sample_stage(rec_title)
    default = _ds_file(tmp_path, stage_fields="", name="a.yaml")
    explicit = _ds_file(tmp_path, stage_fields="    stage_fields: [characteristics_ch1]\n", name="b.yaml")
    assert cm.stage_fields(cm.dataset(ACC, default)) == cm.DEFAULT_STAGE_FIELDS
    pa, pb = _plan(default, tmp_path / "a"), _plan(explicit, tmp_path / "b")
    assert json.dumps(pa["tables"], sort_keys=True).encode() == json.dumps(pb["tables"], sort_keys=True).encode()
    ta, tb = _table(tmp_path / "a").splitlines()[1:], _table(tmp_path / "b").splitlines()[1:]   # line 1: the time
    assert ta == tb
    rows = _rows(pa)
    # GSE204906's form gives none under the default rule, as a second reader found on the real records
    assert rows["GSM1"]["stage"] is None and rows["GSM1"]["stage_from"] == "0 stage values"
    assert (rows["GSM8"]["stage"], rows["GSM8"]["stage_from"]) == (8.5, "characteristics")
    assert (rows["GSM6"]["stage"], rows["GSM6"]["stage_from"]) == (8.5, "characteristics")   # the title fallback


def test_a_pinned_row_is_byte_identical_under_the_default(tmp_path):
    """A plan row of a synthetic record, serialised, is the same bytes with stage_fields absent and with the default
    written out; and its stage fields are exactly those of the rule from before."""
    s = {"GSM8": SAMPLES["GSM8"]}
    a = _plan(_ds_file(tmp_path, stage_fields="", name="a.yaml"), tmp_path / "a", s)
    b = _plan(_ds_file(tmp_path, stage_fields="    stage_fields: [characteristics_ch1]\n", name="b.yaml"),
              tmp_path / "b", s)
    ra, rb = a["tables"][ACC][0], b["tables"][ACC][0]
    assert json.dumps(ra, sort_keys=True) == json.dumps(rb, sort_keys=True)
    assert (ra["stage"], ra["stage_from"]) == (8.5, "characteristics")


@pytest.mark.parametrize("bad", ["    stage_fields: []\n", "    stage_fields: source_name_ch1\n",
                                 "    stage_fields: [source_name_ch1, source_name_ch1]\n",
                                 "    stage_fields: ['source name']\n", "    stage_fields: ['Sample_x/y']\n",
                                 "    stage_fields: [1]\n", "    stage_fields: ['contact_zip/postal_code']\n"])
def test_only_plain_field_names_are_accepted(tmp_path, bad):
    with pytest.raises(cm.Refused, match="stage_fields"):
        cm.load_datasets(_ds_file(tmp_path, stage_fields=bad))


def test_pins_verify_under_the_same_rule_and_a_changed_rule_refuses(tmp_path):
    dsy = _ds_file(tmp_path)
    dest = tmp_path / "d"
    _plan(dsy, dest)
    block = "".join(fe.pin_template(ACC, dest, arm, "Hoxz", gsms, datasets=dsy)
                    for arm, gsms in (("ko", ["GSM1"]), ("control", ["GSM2", "GSM3"])))
    assert "        stage: E8.5" in block
    series = fe.series_pin_template(ACC, dest, datasets=dsy)
    dsy.write_text(dsy.read_text().replace(PIN, "    pinned_series:\n" + series + "    pinned:\n" + block))
    pl = _plan(dsy, dest)                                   # the live check of the pins at plan time
    assert pl["status"] == "ok", pl.get("reason")
    got = fe.verify_pins(cm.dataset(ACC, dsy), json.loads((dest / "meta" / ACC / "plan.json").read_text()))
    assert {s["accession"]: s["stage"] for s in got} == {"GSM1": 8.5, "GSM2": 8.5, "GSM3": 8.5}
    # the same pins under the default rule: the stage no longer parses, so verify_pins (fetch, process) refuses
    stored = json.loads((dest / "meta" / ACC / "plan.json").read_text())
    other = tmp_path / "other.yaml"
    other.write_text(dsy.read_text().replace(SF, ""))
    with pytest.raises(cm.Refused, match="pinned stage E8.5 but the record gives 0 stage values"):
        fe.verify_pins(cm.dataset(ACC, other), stored)
    # and pin-template refuses a stored plan made under another stage rule
    with pytest.raises(cm.Refused, match="plan again"):
        fe.pin_template(ACC, dest, "control", "Hoxz", ["GSM2"], datasets=other)


@pytest.mark.parametrize("with_sf", [True, False])
@pytest.mark.parametrize("where", ["source_name_ch1", "description"])
def test_a_listed_field_that_names_a_held_out_term_is_screened_as_before(tmp_path, with_sf, where):
    """The screen reads every field whatever stage_fields says: a held-out term in a listed stage field refuses the
    dataset and withholds its text, with or without stage_fields."""
    bad = {**SAMPLES["GSM2"], where: "Ctnnb1 null E8.5 embryos"}
    dsy = _ds_file(tmp_path, stage_fields=SF if with_sf else "")
    pl = _plan(dsy, tmp_path / "d", {**SAMPLES, "GSM2": bad})
    assert pl["status"] == "refused" and pl["withheld"] and "held-out" in pl["reason"]
    assert _rows(pl)["GSM2"]["verdict"] == "held_out"
    assert "Ctnnb1" not in _table(tmp_path / "d") and "Ctnnb1" not in \
        (tmp_path / "d" / "meta" / ACC / "plan.json").read_text()
    with pytest.raises(cm.Refused, match="withheld|held"):
        fe.pin_template(ACC, tmp_path / "d", "control", "Hoxz", ["GSM2"], datasets=dsy)


def test_the_committed_list_gives_stage_fields_to_gse204906_only():
    allds = cm.load_datasets()
    assert allds["GSE204906"]["stage_fields"] == ["characteristics_ch1", "source_name_ch1", "description"]
    assert [a for a, ds in allds.items() if "stage_fields" in ds] == ["GSE204906"]
    assert all(cm.stage_fields(ds) == cm.DEFAULT_STAGE_FIELDS for a, ds in allds.items() if a != "GSE204906")
