"""fetch.py's plan tables, pins and fetch (no network: synthetic SOFT text and a stub download), process.py's readers
and a synthetic end-to-end run on pinned samples only, combine.py's weights and refusals, the held datasets that no flag
can fetch, and the rule that no step writes inside the repository or feeds a T1/T2 builder."""

from __future__ import annotations

import builtins
import gzip
import hashlib
import importlib.util
import io
import json
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy import io as spio
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
T3 = ROOT / "experiments" / "t3_decomp"
for p in (EXT, T3):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(f"t3_ext_{name}", EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


fe, pr, co = _load("fetch"), _load("process"), _load("combine")

PIN1 = "    pinned: []  # GSE900001\n"
PIN2 = "    pinned: []  # GSE900002\n"
DS_YAML = """version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
  - accession: GSE900001
    status: cleared
    tier: cardiac
    study: GSE900000
    design: conditional
    control_relation: littermate
    repository: geo
    screen_series: [GSE900001]
    stages_used: [E8.5, E9.5]
    arms:
      control: {on: text, all: ['genotype: wt']}
      ko:
        - {gene: Hoxz, on: text, all: ['genotype: hoxz ko']}
    exclude_match: ['atac']
    expect: {E8.5: {ko: 1, control: 1}}
    sample_file_patterns: ['(matrix\\.mtx|barcodes\\.tsv|features\\.tsv)(\\.gz)?$']
""" + PIN1 + """  - accession: GSE900002
    status: cleared
    tier: unrelated
    study: GSE900002
    design: f0_crispr
    control_relation: same_series
    repository: geo
    screen_series: [GSE900002]
    stages_used: [E8.5]
    chemistry: true
    arms:
      control: {on: title, all: ['(?<![a-z])wt(?![a-z])']}
      ko:
        - {gene: Aaa, on: title, all: ['aaa']}
        - {gene: Bbb, on: title, all: ['bbb']}
    sample_file_patterns: ['(matrix\\.mtx|barcodes\\.tsv|features\\.tsv)(\\.gz)?$']
""" + PIN2 + """  - accession: GSE900003
    status: held_phenocopy
    tier: wnt
    study: GSE900003
    design: conditional
    control_relation: same_capture
    repository: geo
    screen_series: [GSE900003]
    stages_used: [E8.5]
"""


def _soft(samples: dict[str, dict], title="A series", design="embryos", summary="a summary", acc="GSE900001") -> str:
    lines = [f"^SERIES = {acc}", f"!Series_title = {title}", f"!Series_overall_design = {design}",
             f"!Series_summary = {summary}"]
    for gsm, f in samples.items():
        lines.append(f"^SAMPLE = {gsm}")
        for k, vals in f.items():
            for v in (vals if isinstance(vals, list) else [vals]):
                lines.append(f"!Sample_{k} = {v}")
    return "\n".join(lines) + "\n"


def _sample(title, geno, stage="E8.5", files=("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"), extra=None):
    gsm = title
    sup = [f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9nnn/{gsm}/suppl/{gsm}_{f}" for f in files]
    return {"title": title, "source_name_ch1": "whole embryo",
            "characteristics_ch1": [f"genotype: {geno}", f"developmental stage: {stage}"] + (extra or []),
            **{f"supplementary_file_{i + 1}": s for i, s in enumerate(sup)}}


@pytest.fixture()
def dsy(tmp_path):
    p = tmp_path / "datasets.yaml"
    p.write_text(DS_YAML)
    return p


def _plan(dsy, dest, samples, acc="GSE900001", extra_series=None, **kw):
    metas = {acc: fe.parse_soft(_soft(samples, acc=acc, **kw))}
    metas.update(extra_series or {})
    return fe.plan(acc, dest, dsy, metas=metas)


def _pin(dsy, dest, acc, groups):
    """Paste pin-template output (one call per (arm, gene, gsms)) into the YAML, as the readers would."""
    marker = PIN1 if acc == "GSE900001" else PIN2
    block = "".join(fe.pin_template(acc, dest, arm, gene, gsms, datasets=dsy) for arm, gene, gsms in groups)
    dsy.write_text(dsy.read_text().replace(marker, "    pinned:\n" + block))


TWO = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT")}


# ---------------------------------------------------------------------------------------------------------------------
# the plan tables


def test_plan_writes_a_table_of_every_sample_with_verbatim_fields_hash_and_verdict(dsy, tmp_path):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT"),
         "GSM3": _sample("GSM3", "Hoxz KO", stage="E9.5"), "GSM4": _sample("GSM4", "WT", stage="E9.5"),
         "GSM5": _sample("GSM5", "Hoxz KO", stage="E7.5"),
         "GSM6": _sample("GSM6", "WT", extra=["assay: ATAC"], files=("atac_fragments.tsv.gz",))}
    pl = _plan(dsy, tmp_path / "d", s)
    assert pl["status"] == "ok", pl.get("reason")
    rows = {r["accession"]: r for r in pl["tables"]["GSE900001"]}
    assert sorted(rows) == ["GSM1", "GSM2", "GSM3", "GSM4", "GSM5", "GSM6"]   # every sample, chosen or not
    parsed = fe.parse_soft(_soft(s))["samples"]
    for g, r in rows.items():
        assert r["fields"] == parsed[g] and r["verdict"] == "pass"
        assert r["metadata_sha256"] == cm.metadata_sha256(g, parsed[g])
    assert rows["GSM1"]["suggested"] == {"arm": "ko", "gene": "Hoxz"} and rows["GSM6"]["suggested"]["arm"] is None
    assert [f["name"] for f in rows["GSM1"]["files"]] == ["GSM1_matrix.mtx.gz", "GSM1_barcodes.tsv.gz",
                                                          "GSM1_features.tsv.gz"]
    assert rows["GSM6"]["files"] == [] and rows["GSM6"]["other_files"] == ["GSM6_atac_fragments.tsv.gz"]
    assert all(f["url"].startswith("https://ftp.ncbi.nlm.nih.gov/") for f in rows["GSM1"]["files"])
    txt = (tmp_path / "d" / "meta" / "GSE900001" / "plan_table_GSE900001.txt").read_text()
    assert "genotype: Hoxz KO" in txt and rows["GSM5"]["metadata_sha256"] in txt and "verdict pass" in txt
    assert "pinned" not in pl and (tmp_path / "d" / "meta" / "GSE900001" / "plan.json").is_file()


def test_the_canonical_serialisation_is_exactly_as_documented():
    fields = {"title": ["Wt  e8.5 "], "characteristics_ch1": ["genotype: wt", "é stage: E8.5"]}
    want = ('{"accession":"GSM7","fields":{"characteristics_ch1":["genotype: wt","é stage: E8.5"],'
            '"title":["Wt  e8.5 "]}}').encode()
    assert cm.canonical_metadata("GSM7", fields) == want
    assert cm.metadata_sha256("GSM7", fields) == hashlib.sha256(want).hexdigest()
    # the order of values and every character count; the order of keys does not
    assert cm.metadata_sha256("GSM7", {**fields, "characteristics_ch1": fields["characteristics_ch1"][::-1]}) != \
        cm.metadata_sha256("GSM7", fields)
    assert cm.metadata_sha256("GSM7", {**fields, "title": ["Wt e8.5 "]}) != cm.metadata_sha256("GSM7", fields)
    assert cm.metadata_sha256("GSM7", dict(reversed(fields.items()))) == cm.metadata_sha256("GSM7", fields)
    assert cm.metadata_sha256("GSM8", fields) != cm.metadata_sha256("GSM7", fields)
    with pytest.raises(cm.Refused):
        cm.canonical_metadata("GSM7", {"title": "not a list"})


@pytest.mark.parametrize("bad", ["Ctnnb1 fl/fl", "Gata4 KO", "beta-catenin shRNA", "XAV939 treated", "bcatGOF"])
def test_a_held_out_hit_refuses_the_dataset_and_withholds_all_its_text(dsy, tmp_path, bad):
    s = {**TWO, "GSM9": _sample("GSM9", bad)}
    pl = _plan(dsy, tmp_path / "d", s)
    assert pl["status"] == "refused" and pl["withheld"] and "held-out" in pl["reason"]
    written = (tmp_path / "d" / "meta" / "GSE900001" / "plan.json").read_text() + \
        (tmp_path / "d" / "meta" / "GSE900001" / "plan_table_GSE900001.txt").read_text()
    assert bad not in written and "Hoxz" not in written and "whole embryo" not in written
    assert {r["accession"]: r["verdict"] for r in pl["tables"]["GSE900001"]}["GSM9"] == "held_out"
    with pytest.raises(cm.Refused, match="held-out"):
        fe.pin_template("GSE900001", tmp_path / "d", "ko", "Hoxz", ["GSM1"], datasets=dsy)


def test_plan_refuses_review_terms_and_held_out_series_titles(dsy, tmp_path):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT; Gata6 F/+")}
    pl = _plan(dsy, tmp_path / "d", s)
    assert pl["status"] == "refused" and not pl["withheld"]          # a review refusal keeps the table readable
    assert _plan(dsy, tmp_path / "d", TWO, title="Gata4 and Hoxz in the heart")["withheld"]


@pytest.mark.parametrize("key", ["title", "design", "summary"])
@pytest.mark.parametrize("text,kind", [("Gata6 and Hoxz", "review"), ("Wnt signalling inhibition", "review"),
                                       ("GATA binding protein 4", "held-out"), ("sgbcat embryos", "held-out")])
def test_series_title_design_and_summary_refuse_alike(dsy, tmp_path, key, text, kind):
    pl = _plan(dsy, tmp_path / "d", TWO, **{key: text})
    assert pl["status"] == "refused" and kind in pl["reason"]


def test_a_series_review_term_is_cleared_only_by_a_named_series_entry(dsy, tmp_path):
    ok = DS_YAML.replace("    exclude_match: ['atac']\n", "    exclude_match: ['atac']\n    review_cleared:\n"
                         "      - {sample: GSE900001, term: gata6, reason: names a neighbour only}\n", 1)
    dsy.write_text(ok)
    assert _plan(dsy, tmp_path / "d", TWO, summary="Gata6 is mentioned")["status"] == "ok"
    assert _plan(dsy, tmp_path / "d", TWO, summary="Gata6 and Gata4 are mentioned")["status"] == "refused"
    assert _plan(dsy, tmp_path / "d", TWO, summary="Gata6 and Ctnb1 are mentioned")["status"] == "refused"  # void
    dsy.write_text(ok.replace("sample: GSE900001", "sample: GSM1"))
    assert _plan(dsy, tmp_path / "d", TWO, summary="Gata6 is mentioned")["status"] == "refused"


def test_exclude_match_on_title_ignores_the_other_fields(dsy):
    ds = {"exclude_match": ["(?<![a-z])atac(?![a-z])"], "exclude_match_on": "title",
          "arms": {"control": {"on": "text", "all": ["genotype: wt"]}, "ko": []}}
    rna = {**_sample("E8.5_WT_rep1", "WT", extra=["library: multiome RNA (paired ATAC)"]), "title": ["E8.5_WT_rep1"]}
    atac = {**_sample("E8.5_WT_rep1_ATAC", "WT"), "title": ["E8.5_WT_rep1_ATAC"]}
    assert fe.resolve_arm(ds, rna)[0] == "control"
    assert fe.resolve_arm(ds, atac)[0] is None
    assert fe.resolve_arm({**ds, "exclude_match_on": "text"}, rna)[0] is None


# ---------------------------------------------------------------------------------------------------------------------
# pinning


def test_pin_template_prints_a_block_that_loads_and_never_writes_the_yaml(dsy, tmp_path):
    _plan(dsy, tmp_path / "d", TWO)
    before = dsy.read_bytes()
    out = fe.pin_template("GSE900001", tmp_path / "d", "ko", "Hoxz", ["GSM1"], datasets=dsy)
    assert dsy.read_bytes() == before
    pins = yaml.safe_load("pinned:\n" + out)["pinned"]
    rows = {r["accession"]: r for r in json.loads((tmp_path / "d" / "meta" / "GSE900001" / "plan.json")
                                                  .read_text())["tables"]["GSE900001"]}
    assert pins == [{"accession": "GSM1", "arm": "ko", "gene": "Hoxz", "stage": "E8.5",
                     "files": ["GSM1_matrix.mtx.gz", "GSM1_barcodes.tsv.gz", "GSM1_features.tsv.gz"],
                     "metadata_sha256": rows["GSM1"]["metadata_sha256"]}]
    assert fe.pin_template("GSE900001", tmp_path / "d", "control", "Hoxz", suggested=True, datasets=dsy).count(
        "accession: GSM2") == 1
    assert "GSM1" in fe.pin_template("GSE900001", tmp_path / "d", "ko", "Hoxz", title="^gsm1$", datasets=dsy)
    for kw in ({"arm": "KO"}, {"gsms": None}, {"gsms": ["GSM404"]}):
        with pytest.raises(cm.Refused):
            fe.pin_template("GSE900001", tmp_path / "d", kw.get("arm", "ko"), "Hoxz", kw.get("gsms", ["GSM1"]),
                            datasets=dsy)
    _pin(dsy, tmp_path / "d", "GSE900001", [("ko", "Hoxz", ["GSM1"]), ("control", "Hoxz", ["GSM2"])])
    ds = cm.dataset("GSE900001", dsy)
    assert [p["accession"] for p in ds["pinned"]] == ["GSM1", "GSM2"]
    assert cm.pinned_units(ds["pinned"]) == [{"gene": "Hoxz", "ko": ["GSM1"], "control": ["GSM2"], "stages": [8.5]}]


def test_pin_template_refuses_a_sample_without_one_stage(dsy, tmp_path):
    s = {**TWO, "GSM3": _sample("GSM3", "Hoxz KO", stage="E9.5-E10.5")}
    _plan(dsy, tmp_path / "d", s)
    with pytest.raises(cm.Refused, match="stage"):
        fe.pin_template("GSE900001", tmp_path / "d", "ko", "Hoxz", ["GSM3"], datasets=dsy)


@pytest.mark.parametrize("edit,match", [
    (lambda p: {**p, "arm": "KO"}, "exactly"), (lambda p: {**p, "arm": "wt"}, "exactly"),
    (lambda p: {**p, "stage": "E10.5"}, "after E9.5"), (lambda p: {**p, "stage": "E7.5"}, "stages_used"),
    (lambda p: {**p, "metadata_sha256": "abc"}, "hex"), (lambda p: {**p, "files": []}, "files"),
    (lambda p: {**p, "accession": "SRR1"}, "GSM"), (lambda p: {k: v for k, v in p.items() if k != "gene"}, "keys"),
    (lambda p: {**p, "note": "x"}, "keys")])
def test_the_pinned_block_form_is_checked(dsy, edit, match):
    good = {"accession": "GSM1", "arm": "ko", "gene": "Hoxz", "stage": "E8.5", "files": ["a.mtx.gz"],
            "metadata_sha256": "a" * 64}
    ds = cm.load_datasets(dsy)["GSE900001"]
    assert cm.check_pins("GSE900001", {**ds, "pinned": [good]}) == [good]
    with pytest.raises(cm.Refused, match=match):
        cm.check_pins("GSE900001", {**ds, "pinned": [edit(good)]})
    with pytest.raises(cm.Refused, match="two different"):
        cm.check_pins("GSE900001", {**ds, "pinned": [good, {**good, "gene": "Other", "stage": "E9.5"}]})
    with pytest.raises(cm.Refused, match="only a cleared"):
        cm.check_pins("GSE900003", {**ds, "status": "held_phenocopy", "pinned": [good]})
    with pytest.raises(cm.Refused, match="per-cell"):
        cm.check_pins("GSE900001", {**ds, "cell_labels": {"genotype_columns": ["g"]}, "pinned": [good]})
    with pytest.raises(cm.Refused, match="without a knockout"):
        cm.pinned_units([{**good, "arm": "control"}])
    with pytest.raises(cm.Refused, match="every stage needs both arms"):
        cm.pinned_units([good, {**good, "accession": "GSM2", "arm": "control", "stage": "E9.5"}])


def test_a_sample_pinned_under_two_datasets_is_refused(dsy, tmp_path):
    p = ("    pinned:\n      - {accession: GSM1, arm: ko, gene: Aaa, stage: E8.5, files: [x.gz], metadata_sha256: "
         f"'{'a' * 64}'}}\n      - {{accession: GSM2, arm: control, gene: Aaa, stage: E8.5, files: [y.gz], "
         f"metadata_sha256: '{'b' * 64}'}}\n")
    dsy.write_text(DS_YAML.replace(PIN2, p).replace(PIN1, p.replace("Aaa", "Hoxz")))
    with pytest.raises(cm.Refused, match="pinned under"):
        cm.load_datasets(dsy)


# ---------------------------------------------------------------------------------------------------------------------
# fetch: pinned samples only


def _stage_files(staging: Path, gsms, rng, mu, folds):
    for gsm in gsms:
        _write_sample(staging / gsm, gsm, {"CM": 160, "Neural": 160}, mu, folds.get(gsm, np.ones(len(PANEL))), rng)


def _fake_get(staging: Path, calls: list):
    def get(url, dst, max_bytes):
        name = Path(url).name
        gsm = name.split("_")[0]
        calls.append(url)
        dst = cm.outside_repo(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(staging / gsm / name, dst)
        return dst.stat().st_size, cm.sha256_file(dst)
    return get


def test_fetch_refuses_unpinned_unlisted_held_or_unplanned_datasets_and_paths_inside_the_repo(dsy, tmp_path):
    with pytest.raises(cm.Refused):
        fe.plan("GSE999999", tmp_path / "d", dsy, metas={})
    with pytest.raises(cm.Refused, match="held_phenocopy"):
        fe.plan("GSE900003", tmp_path / "d", dsy, metas={})
    with pytest.raises(cm.Refused, match="inside the repository"):
        fe.plan("GSE900001", ROOT / "scratch_dl", dsy, metas={"GSE900001": fe.parse_soft(_soft({}))})
    assert _plan(dsy, tmp_path / "d", TWO)["status"] == "ok"
    with pytest.raises(cm.Refused, match="pinned block is empty"):
        fe.fetch("GSE900001", tmp_path / "d", dsy)
    _pin(dsy, tmp_path / "d", "GSE900001", [("ko", "Hoxz", ["GSM1"]), ("control", "Hoxz", ["GSM2"])])
    with pytest.raises(cm.Refused, match="changed since the plan"):
        fe.fetch("GSE900001", tmp_path / "d", dsy)                  # the YAML changed: plan again
    with pytest.raises(cm.Refused, match="no plan"):
        fe.fetch("GSE900001", tmp_path / "empty", dsy)
    pl = _plan(dsy, tmp_path / "d", {**TWO, "GSM9": _sample("GSM9", "Ctnnb1 KO")})
    assert pl["status"] == "refused"
    with pytest.raises(cm.Refused, match="did not pass"):
        fe.fetch("GSE900001", tmp_path / "d", dsy)
    for url in ("https://example.org/x.mtx.gz", "http://ftp.ncbi.nlm.nih.gov/x", "file:///etc/passwd"):
        with pytest.raises(cm.Refused):
            fe.https_url(url)
    assert fe.DENY.search("GSM1_atac_fragments.tsv.gz") and fe.DENY.search("x_R1.fastq.gz")
    assert fe.DENY.search("GSE1_RAW_raw.tar") and not fe.DENY.search("GSM1_matrix.mtx.gz")


def test_fetch_requests_only_the_pinned_files_of_pinned_samples(dsy, tmp_path, monkeypatch):
    s = {**TWO, "GSM3": _sample("GSM3", "Hoxz KO")}
    _plan(dsy, tmp_path / "d", s)
    _pin(dsy, tmp_path / "d", "GSE900001", [("ko", "Hoxz", ["GSM1"]), ("control", "Hoxz", ["GSM2"])])
    assert _plan(dsy, tmp_path / "d", s)["status"] == "ok"
    rng = np.random.default_rng(3)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    _stage_files(tmp_path / "staging", ["GSM1", "GSM2", "GSM3"], rng, mu, {})
    calls: list = []
    monkeypatch.setattr(fe, "_get", _fake_get(tmp_path / "staging", calls))
    man = fe.fetch("GSE900001", tmp_path / "d", dsy)
    assert sorted({Path(u).name.split("_")[0] for u in calls}) == ["GSM1", "GSM2"] and len(calls) == 6
    assert man["pins_sha256"] == cm.pins_sha256(cm.dataset("GSE900001", dsy)["pinned"])
    assert not (tmp_path / "d" / "raw" / "GSE900001" / "GSM3").exists()


def _pinned_plan(dsy, tmp_path, s=None, extra_series=None):
    s = s or TWO
    _plan(dsy, tmp_path / "d", s, extra_series=extra_series)
    _pin(dsy, tmp_path / "d", "GSE900001", [("ko", "Hoxz", ["GSM1"]), ("control", "Hoxz", ["GSM2"])])
    return _plan(dsy, tmp_path / "d", s, extra_series=extra_series)


def test_a_pinned_sample_whose_live_or_stored_metadata_changed_is_refused(dsy, tmp_path):
    assert _pinned_plan(dsy, tmp_path)["status"] == "ok"
    changed = {"GSM1": {**_sample("GSM1", "Hoxz KO"), "description": "a corrected title"}, "GSM2": TWO["GSM2"]}
    live = _plan(dsy, tmp_path / "d2", changed)
    assert live["status"] == "refused" and "pins:" in live["reason"] and "live metadata" in live["reason"]
    pp = tmp_path / "d" / "meta" / "GSE900001" / "plan.json"
    pl = json.loads(pp.read_text())
    pl["tables"]["GSE900001"][0]["fields"]["title"] = ["GSM1 edited"]
    pp.write_text(json.dumps(pl))
    with pytest.raises(cm.Refused, match="stored metadata"):
        fe.fetch("GSE900001", tmp_path / "d", dsy)


def test_a_pinned_accession_must_be_in_its_own_series_table(dsy, tmp_path):
    other = {"GSE900009": fe.parse_soft(_soft({"GSM8": _sample("GSM8", "WT")}, acc="GSE900009"))}
    assert _pinned_plan(dsy, tmp_path, extra_series=other)["status"] == "ok"
    ds = cm.dataset("GSE900001", dsy)
    pl = json.loads((tmp_path / "d" / "meta" / "GSE900001" / "plan.json").read_text())
    p8 = {**ds["pinned"][1], "accession": "GSM8"}
    with pytest.raises(cm.Refused, match="belongs to another series"):
        fe.verify_pins({**ds, "pinned": [ds["pinned"][0], p8]}, pl)
    with pytest.raises(cm.Refused, match="not in the plan table"):
        fe.verify_pins({**ds, "pinned": [ds["pinned"][0], {**p8, "accession": "GSM77"}]}, pl)
    with pytest.raises(cm.Refused, match="pinned stage"):
        fe.verify_pins({**ds, "pinned": [{**p, "stage": "E9.5"} for p in ds["pinned"]]}, pl)
    with pytest.raises(cm.Refused, match="not per-sample files"):
        fe.verify_pins({**ds, "pinned": [{**ds["pinned"][0], "files": ["GSM2_matrix.mtx.gz"]}, ds["pinned"][1]]}, pl)


def test_pinned_chemistry_must_match_its_control(dsy, tmp_path):
    def smp(t, chem):
        r = _sample(t, "x")
        r["extract_protocol_ch1"] = f"Chromium Single Cell 3' {chem} chemistry"
        return r
    s = {"GSM11": smp("GSM11", "v3"), "GSM12": smp("GSM12", "v2"), "GSM13": smp("GSM13", "v3")}
    _plan(dsy, tmp_path / "d", s, acc="GSE900002")
    marker_yaml = dsy.read_text()
    _pin(dsy, tmp_path / "d", "GSE900002", [("ko", "Aaa", ["GSM11"]), ("control", "Aaa", ["GSM12"])])
    assert "chemistry" in _plan(dsy, tmp_path / "d", s, acc="GSE900002")["reason"]
    dsy.write_text(marker_yaml)
    _pin(dsy, tmp_path / "d", "GSE900002", [("ko", "Aaa", ["GSM11"]), ("control", "Aaa", ["GSM13"])])
    assert _plan(dsy, tmp_path / "d", s, acc="GSE900002")["status"] == "ok"


# ---------------------------------------------------------------------------------------------------------------------
# held datasets: no flag fetches them


HELD_CLEARED = """version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
  - {accession: GSE165300, status: cleared, tier: wnt, study: GSE165300, design: conditional,
     control_relation: same_capture, repository: geo, screen_series: [GSE165300], stages_used: [E8.5]}
"""


def test_a_never_fetch_accession_is_refused_whatever_a_yaml_says(tmp_path):
    y = tmp_path / "alt.yaml"
    y.write_text(HELD_CLEARED)
    with pytest.raises(cm.Refused, match="NEVER_FETCH"):
        cm.load_datasets(y)
    for acc in sorted(cm.NEVER_FETCH):
        with pytest.raises(cm.Refused, match="NEVER_FETCH"):
            cm.dataset(acc, y)
        with pytest.raises(cm.Refused, match="NEVER_FETCH"):
            fe.plan(acc, tmp_path / "d", y, metas={})
        with pytest.raises(cm.Refused, match="NEVER_FETCH"):
            fe.fetch(acc, tmp_path / "d", y)
        with pytest.raises(cm.Refused, match="NEVER_FETCH"):
            pr.run(acc, tmp_path / "d", tmp_path / "wt", tmp_path / "r", tmp_path / "o", y)
    assert cm.NEVER_FETCH == {"GSE165300", "GSE126128", "GSE280543", "GSE259342", "E-MTAB-7325"}
    for status in ("held_phenocopy", "excluded_design"):
        y.write_text(HELD_CLEARED.replace("status: cleared", f"status: {status}"))
        assert cm.load_datasets(y)["GSE165300"]["status"] == status


def test_no_cli_takes_another_dataset_list(tmp_path, capsys):
    y = tmp_path / "alt.yaml"
    y.write_text(HELD_CLEARED)
    for main, argv in ((fe.main, ["plan", "GSE165300", "--dest", str(tmp_path / "d"), "--datasets", str(y)]),
                       (fe.main, ["fetch", "GSE165300", "--dest", str(tmp_path / "d"), "--datasets", str(y)]),
                       (pr.main, ["run", "GSE165300", "--dest", "d", "--wt", "w", "--refstats", "r", "--out", "o",
                                  "--datasets", str(y)]),
                       (co.main, ["--processed", "p", "--refstats", "r", "--out", "o", "--datasets", str(y)])):
        with pytest.raises(SystemExit) as e:
            main(argv)
        assert e.value.code == 2
    assert fe.main(["plan", "GSE165300", "--dest", str(tmp_path / "d")]) == 2
    assert fe.main(["fetch", "GSE126128", "--dest", str(tmp_path / "d")]) == 2
    assert pr.main(["run", "GSE280543", "--dest", str(tmp_path / "d"), "--wt", "w", "--refstats", "r",
                    "--out", str(tmp_path / "o")]) == 2
    assert "NEVER_FETCH" in capsys.readouterr().out
    y.write_text(HELD_CLEARED.replace("status: cleared", "status: held_phenocopy"))
    q = tmp_path / "processed" / "GSE165300"
    q.mkdir(parents=True)
    (q / "QC.json").write_text(json.dumps({"accession": "GSE165300", "status": "ok", "units": {}}))
    with pytest.raises(cm.Refused, match="NEVER_FETCH"):
        co.load_units(tmp_path / "processed", y, {})


def test_the_pinned_gate_and_the_datasets_sha_command(dsy, tmp_path, capsys):
    assert fe.check_gate(["GSE900001"], dsy) == ["GSE900001: the pinned block is empty"]
    assert "NEVER_FETCH" in fe.check_gate(["GSE165300"], dsy)[0]
    _plan(dsy, tmp_path / "d", TWO)
    _pin(dsy, tmp_path / "d", "GSE900001", [("ko", "Hoxz", ["GSM1"]), ("control", "Hoxz", ["GSM2"])])
    assert fe.check_gate(["GSE900001"], dsy) == []
    # the committed file: its five cleared datasets have no pins yet, so the gate refuses every one
    assert fe.main(["check-pins", "GSE205117", "GSE168768"]) == 2
    assert fe.main(["datasets-sha"]) == 0
    assert cm.DATASETS_SHA256 in capsys.readouterr().out


# ---------------------------------------------------------------------------------------------------------------------
# readers


def _mtx_bytes(M) -> bytes:
    b = io.BytesIO()
    spio.mmwrite(b, sparse.coo_matrix(M))
    return gzip.compress(b.getvalue())


def test_read_sample_mtx_triple_and_tar_and_panel_mapping(tmp_path):
    feats = ["ENS1\tTbxt\tGene Expression", "ENS2\tHand1\tGene Expression", "ENS3\tHand1\tGene Expression",
             "chr1:1-100\tchr1:1-100\tPeaks", "ENS4\tOther\tGene Expression"]
    M = np.array([[1, 2, 3, 50, 4], [0, 1, 0, 9, 1]]).T  # features x cells
    (tmp_path / "s").mkdir()
    (tmp_path / "s" / "GSM1_matrix.mtx.gz").write_bytes(_mtx_bytes(M))
    (tmp_path / "s" / "GSM1_features.tsv.gz").write_bytes(gzip.compress("\n".join(feats).encode()))
    (tmp_path / "s" / "GSM1_barcodes.tsv.gz").write_bytes(gzip.compress(b"AAA-1\nCCC-1\n"))
    X, sym, types, bcs = pr.read_sample(sorted((tmp_path / "s").iterdir()))
    assert X.shape == (2, 5) and bcs == ["AAA-1", "CCC-1"]
    C, present, rec = pr.to_panel(X, sym, types, ["T", "Hand1", "chr1:1-100", "Absent"])
    assert present.tolist() == [True, True, False, False]
    assert C.toarray().tolist() == [[1, 5, 0, 0], [0, 1, 0, 0]]   # alias T <- Tbxt; duplicate Hand1 summed; no peaks
    assert rec["alias"] == 1 and rec["duplicate_symbols"] == 1
    tar = tmp_path / "GSM2_sample.tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        for f in sorted((tmp_path / "s").iterdir()):
            t.add(f, arcname=f"GSM2/{f.name.replace('GSM1_', '')}")
    X2, sym2, _, bcs2 = pr.read_sample([tar])
    assert (X2 != X).nnz == 0 and sym2 == sym and bcs2 == bcs
    assert not any(p.name.startswith("GSM2") and p.is_dir() for p in tmp_path.iterdir())  # nothing unpacked


# ---------------------------------------------------------------------------------------------------------------------
# a synthetic end-to-end run on pinned samples, and combine


PANEL = [f"g{i}" for i in range(40)]
GROUPS = ["CM", "Endoderm", "Neural"]


def _write_sample(dirp: Path, gsm: str, groups_n: dict, mu: dict, fold, rng) -> list[str]:
    extra = 260
    rows, names = [], []
    for g, n in groups_n.items():
        p = mu[g] * fold
        panel = rng.multinomial(2000, p / p.sum(), size=n)
        rows.append(np.hstack([panel, rng.poisson(2.0, size=(n, extra))]))
        names += [f"{g}_{i}" for i in range(n)]
    M = np.vstack(rows).T
    dirp.mkdir(parents=True, exist_ok=True)
    feats = [f"E{i}\t{s}\tGene Expression" for i, s in enumerate(PANEL + [f"x{j}" for j in range(extra)])]
    (dirp / f"{gsm}_matrix.mtx.gz").write_bytes(_mtx_bytes(M))
    (dirp / f"{gsm}_features.tsv.gz").write_bytes(gzip.compress("\n".join(feats).encode()))
    (dirp / f"{gsm}_barcodes.tsv.gz").write_bytes(gzip.compress("\n".join(names).encode()))
    return names


def _synthetic_run(tmp_path, dsy, monkeypatch, unpinned=False):
    """Plan, pin GSM1 (ko) and GSM2 (control), plan again, fetch through a stub download. With unpinned, the series
    also holds GSM3 (a Hoxz knockout nobody pinned) whose files lie in the raw folder and the manifest."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    wt = tmp_path / "wt.h5ad"
    wt.write_bytes(b"stand-in for the WT reference")
    st = {"genes": PANEL, "groups": GROUPS, "shares": {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5}, "mu": mu,
          "count": {"CM": 20, "Endoderm": 30, "Neural": 50}}
    refp = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats(st, refp, cm.sha256_file(wt))
    dest = tmp_path / "dest"
    s = {**TWO, **({"GSM3": _sample("GSM3", "Hoxz KO")} if unpinned else {})}
    _plan(dsy, dest, s)
    _pin(dsy, dest, "GSE900001", [("ko", "Hoxz", ["GSM1"]), ("control", "Hoxz", ["GSM2"])])
    assert _plan(dsy, dest, s)["status"] == "ok"
    fold = np.ones(len(PANEL))
    fold[:4] = 3.0
    staging = tmp_path / "staging"
    _stage_files(staging, list(s), rng, mu, {"GSM1": fold, "GSM3": fold})
    monkeypatch.setattr(fe, "_get", _fake_get(staging, []))
    fe.fetch("GSE900001", dest, dsy)
    if unpinned:   # the unpinned sample's files are present beside the pinned ones, and listed in the manifest
        mp = dest / "raw" / "GSE900001" / "MANIFEST.json"
        man = json.loads(mp.read_text())
        for p in (staging / "GSM3").iterdir():
            q = dest / "raw" / "GSE900001" / "GSM3" / p.name
            q.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p, q)
            man["files"][f"GSM3/{p.name}"] = {"sha256": cm.sha256_file(q), "size": q.stat().st_size}
        mp.write_text(json.dumps(man))

    def classifier(wt_path, C, present, genes, names, threads):
        calls = np.array([n.split(":")[1].split("_")[0] for n in names])
        return calls, np.full(len(names), 0.9), {"stub": True}

    return dest, wt, refp, classifier


def test_end_to_end_synthetic_run_writes_only_outside_and_combine_takes_it(tmp_path, dsy, monkeypatch):
    before = sorted(p.name for p in EXT.iterdir())
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy, monkeypatch)
    out = tmp_path / "processed"
    qc = pr.run("GSE900001", dest, wt, refp, out, dsy, classifier=clf)
    assert qc["status"] == "ok" and qc["panel_genes_present"] == 40
    assert qc["screened_samples"] == ["GSM1", "GSM2"] and qc["pins_sha256"]
    u = qc["units"]["Hoxz"]
    assert u["status"] == "ok" and set(u["groups"]) == {"CM", "Neural"}
    assert u["cells_per_stage_and_arm"] == {"E8.5:control": 320, "E8.5:ko": 320}
    import build_submission as bs

    phi, _ = bs.load_factor(Path(u["tables"]["GSE900001__Hoxz"]["file"]).parents[1], "GSE900001__Hoxz", PANEL, {})
    assert np.all((phi >= 0.25) & (phi <= 4)) and phi[:4].min() > 1
    r = co.combine(out, refp, tmp_path / "tables", dsy)
    assert set(r["tables"]) == {"study_GSE900000", "cardiac", "all"}
    man = json.loads((tmp_path / "tables" / "cardiac" / "cardiac_manifest.json").read_text())
    assert man["sign_rule"] == "single study" and man["studies"] == ["GSE900000"]
    assert sorted(p.name for p in EXT.iterdir()) == before
    with pytest.raises(cm.Refused, match="inside the repository"):
        pr.run("GSE900001", dest, wt, refp, EXT / "out", dsy, classifier=clf)
    with pytest.raises(cm.Refused, match="inside the repository"):
        co.combine(out, refp, ROOT / "tables_here", dsy)


def test_process_never_opens_a_count_file_of_an_unpinned_sample(tmp_path, dsy, monkeypatch):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy, monkeypatch, unpinned=True)
    assert (dest / "raw" / "GSE900001" / "GSM3").is_dir()
    real_open, real_io_open = builtins.open, io.open

    def guarded(opener):
        def f(file, *a, **k):
            if "GSM3" in str(file):
                raise AssertionError(f"opened a file of the unpinned sample: {file}")
            return opener(file, *a, **k)
        return f

    monkeypatch.setattr(builtins, "open", guarded(real_open))
    monkeypatch.setattr(io, "open", guarded(real_io_open))
    qc = pr.run("GSE900001", dest, wt, refp, tmp_path / "processed", dsy, classifier=clf)
    monkeypatch.setattr(builtins, "open", real_open)
    monkeypatch.setattr(io, "open", real_io_open)
    assert qc["status"] == "ok" and qc["screened_samples"] == ["GSM1", "GSM2"]
    assert set(qc["confidence_by_sample"]) == {"GSM1", "GSM2"}
    # the guard is live: opening the unpinned file under it raises
    monkeypatch.setattr(builtins, "open", guarded(real_open))
    with pytest.raises(AssertionError, match="unpinned"):
        open(next((dest / "raw" / "GSE900001" / "GSM3").iterdir()), "rb")
    monkeypatch.setattr(builtins, "open", real_open)


def test_process_refuses_arms_or_samples_that_are_not_exactly_as_pinned(tmp_path, dsy, monkeypatch):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy, monkeypatch)
    mp = dest / "raw" / "GSE900001" / "MANIFEST.json"
    man = json.loads(mp.read_text())
    good = json.dumps(man)
    man["pinned"][0]["arm"] = "control"
    mp.write_text(json.dumps(man))
    with pytest.raises(cm.Refused, match="not exactly as pinned"):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o", dsy, classifier=clf)
    mp.write_text(good)
    pin_yaml = dsy.read_text()
    dsy.write_text(pin_yaml.replace("arm: control", "arm: ko"))   # the YAML changed after the fetch
    with pytest.raises(cm.Refused):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o", dsy, classifier=clf)
    dsy.write_text(pin_yaml.replace("    pinned:\n", "    pinned: []\n    unused:\n"))
    with pytest.raises(cm.Refused):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o", dsy, classifier=clf)


def test_confident_shares_are_recorded_and_a_table_under_the_floor_is_refused(tmp_path, dsy, monkeypatch):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy, monkeypatch)

    def low(wt_path, C, present, genes, names, threads):
        calls, conf, rec = clf(wt_path, C, present, genes, names, threads)
        conf = np.where([n.startswith("GSM1:") and i % 5 < 3 for i, n in enumerate(names)], 0.3, conf)
        return calls, conf, rec

    out = tmp_path / "processed"
    qc = pr.run("GSE900001", dest, wt, refp, out, dsy, classifier=low)
    cs = qc["confidence_by_sample"]
    assert cs["GSM1"]["confident_share"] == pytest.approx(0.4) and cs["GSM2"]["confident_share"] == 1.0
    assert cs["GSM1"]["unconfident_share"] == pytest.approx(0.6)
    assert qc["units"]["Hoxz"]["status"] == "refused_low_confidence" and qc["status"] == "left_out"
    assert not list((out / "GSE900001").glob("*.unit.json"))
    saved = json.loads((out / "GSE900001" / "QC.json").read_text())
    assert saved["confidence_by_sample"]["GSM1"]["cells"] == 320
    with pytest.raises(cm.Refused, match="no unit passed"):
        co.combine(out, refp, tmp_path / "t", dsy)
    dsy.write_text(dsy.read_text().replace("qc_floors:\n  min_confident_share_per_sample: 0.5\n", ""))
    with pytest.raises(cm.Refused, match="qc_floors"):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o2", dsy, classifier=clf)


def test_process_refuses_changed_raw_files_and_a_pinned_sample_that_fails_the_screen(tmp_path, dsy, monkeypatch):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy, monkeypatch)
    p = next((dest / "raw" / "GSE900001" / "GSM1").glob("*barcodes*"))
    p.write_bytes(gzip.compress(b"changed"))
    with pytest.raises(cm.Refused, match="sha256"):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o", dsy, classifier=clf)
    dsy2 = tmp_path / "b" / "datasets.yaml"
    dsy2.parent.mkdir(parents=True)
    dsy2.write_text(DS_YAML)
    dest, wt, refp, clf = _synthetic_run(tmp_path / "b", dsy2, monkeypatch)
    # a screen term added later refuses the pinned sample before any count is read
    dsy2.write_text(dsy2.read_text().replace("    exclude_match: ['atac']\n",
                                             "    exclude_match: ['atac']\n    screen_extra_refuse: ['hoxz']\n", 1))
    pl_ds = cm.dataset("GSE900001", dsy2)
    pl = json.loads((dest / "meta" / "GSE900001" / "plan.json").read_text())
    with pytest.raises(cm.Refused, match="screen verdict"):
        fe.verify_pins(pl_ds, pl)
    with pytest.raises(cm.Refused):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o2", dsy2, classifier=lambda *a: pytest.fail("read"))


def test_sign_rule_and_equal_study_weights():
    V = np.array([[0.4], [0.2], [-0.1]])
    assert co.sign_rule_mean(V, 3)[0][0] == pytest.approx(0.5 / 3)
    assert co.sign_rule_mean(np.array([[0.4], [-0.3], [0.1]]), 3)[0][0] == 0.0     # opposite study too large
    assert co.sign_rule_mean(np.array([[0.4], [np.nan]]), 2)[0][0] == 0.0         # held by one study of two
    assert co.sign_rule_mean(np.array([[0.4]]), 1)[0][0] == pytest.approx(0.4)    # a one-study table
    assert np.isnan(co.sign_rule_mean(np.array([[np.nan], [np.nan]]), 2)[0][0])
    P = 2
    units = [{"accession": "A1", "gene": "g1", "study": "SA", "tier": "unrelated", "sha256": "1",
              "L": {"CM": np.array([0.2, 0.1])}},
             {"accession": "A1", "gene": "g2", "study": "SA", "tier": "unrelated", "sha256": "2",
              "L": {"CM": np.array([0.6, np.nan])}},
             {"accession": "B1", "gene": "h1", "study": "SB", "tier": "unrelated", "sha256": "3",
              "L": {"CM": np.array([0.1, 0.3])}}]
    st = co.study_tables(units, ["CM", "Neural"], P)
    assert st["SA"]["L"]["CM"].tolist() == pytest.approx([0.4, 0.1])   # equal weight per gene inside the study
    ref = {"genes": ["a", "b"], "groups": ["CM", "Neural"], "shares": {"CM": 0.4, "Neural": 0.6},
           "mu": {"CM": np.array([0.5, 0.5]), "Neural": np.array([0.9, 0.1])}}
    per_group, pooled, agree = co.table_from_studies(st, ["SA", "SB"], ref)
    assert per_group["CM"].tolist() == pytest.approx([0.25, 0.2])       # equal weight per study, not per gene
    assert pooled.tolist() == pytest.approx([0.25, 0.2]) and agree["same_sign_studies"] == [2, 2]


def test_combine_refuses_a_failed_qc_a_changed_unit_or_other_pins(tmp_path, dsy, monkeypatch):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy, monkeypatch)
    out = tmp_path / "processed"
    pr.run("GSE900001", dest, wt, refp, out, dsy, classifier=clf)
    qcp = out / "GSE900001" / "QC.json"
    good = qcp.read_text()
    qc = json.loads(good)
    qc["pins_sha256"] = "0" * 64
    qcp.write_text(json.dumps(qc))
    with pytest.raises(cm.Refused, match="other pins"):
        co.combine(out, refp, tmp_path / "t", dsy)
    qcp.write_text(good)
    unit = next((out / "GSE900001").glob("*.unit.json"))
    unit.write_text(unit.read_text().replace('"gene"', '"gene" ', 1))
    with pytest.raises(cm.Refused, match="changed"):
        co.combine(out, refp, tmp_path / "t", dsy)
    qc = json.loads(good)
    qc["status"] = "failed"
    qcp.write_text(json.dumps(qc))
    with pytest.raises(cm.Refused, match="QC failed"):
        co.combine(out, refp, tmp_path / "t", dsy)


# ---------------------------------------------------------------------------------------------------------------------
# the launcher and the guards


def test_the_launcher_plans_by_default_and_gates_on_pins_before_any_fetch():
    sh = EXT / "run_dgx.sh"
    assert subprocess.run(["bash", "-n", str(sh)], capture_output=True).returncode == 0
    text = sh.read_text()
    assert 'MODE=${2:-plan}' in text
    assert """trap 'echo 99 > "$D/logs/ALL.exit"; exit 99' TERM INT HUP""" in text
    plan_part, all_part = text.split('if [ "$MODE" = plan ]; then', 1)[1].split("\nfi\n", 1)
    assert "fetch.py plan" in plan_part and "fetch.py fetch" not in plan_part and "process.py" not in plan_part
    gate = all_part.index("fetch.py check-pins")
    assert gate < all_part.index("fetch.py fetch") and gate < all_part.index("process.py")
    assert 'echo 6 > "$D/logs/ALL.exit"' in all_part[:all_part.index("fetch.py fetch")]
    assert all_part.index('step "replan_$A"') < all_part.index('step "fetch_$A"')
    assert "bundle" not in text and "worktree" not in text


def test_the_launcher_refuses_the_reference_checkout(tmp_path):
    home = tmp_path / "home"
    ref = home / "v" / "virtual-embryo-challenge"
    ref.mkdir(parents=True)
    env = {"HOME": str(home), "PATH": "/usr/bin:/bin"}
    r = subprocess.run(["bash", str(EXT / "run_dgx.sh"), str(ref), "plan"], capture_output=True, text=True, env=env,
                       timeout=60)
    assert r.returncode == 3
    logs = home / "v" / "vec-ext" / "t3" / "logs"
    assert (logs / "ALL.exit").read_text().strip() == "3"
    assert "reference checkout" in (logs / "start_plan.txt").read_text()
    assert sorted(p.name for p in ref.iterdir()) == []


def test_outside_repo_guard(tmp_path, monkeypatch):
    with pytest.raises(cm.Refused):
        cm.outside_repo(ROOT / "x" / "y")
    with pytest.raises(cm.Refused):
        cm.outside_repo(ROOT)
    with pytest.raises(cm.Refused):
        tb.write_table(ROOT / "t", "n", ["a"], np.ones(1), {}, {})
    assert cm.outside_repo(tmp_path / "x") == (tmp_path / "x").resolve()
    monkeypatch.setattr(cm, "REFERENCE_CHECKOUT", tmp_path / "refco")
    with pytest.raises(cm.Refused, match="reference checkout"):
        cm.outside_repo(tmp_path / "refco" / ".git" / "objects")


def test_no_t1_or_t2_builder_reads_the_outside_knockout_code():
    files = [*(ROOT / "src" / "vec").glob("*.py"), *(ROOT / "experiments" / "builds").glob("*.py"),
             *(ROOT / "experiments" / "builds").glob("*.sh"), *(ROOT / "scripts").glob("*.py"),
             *(ROOT / "scripts").glob("*.sh")]
    for d in sorted((ROOT / "experiments").glob("t[12]_*")):
        files += [*d.glob("*.py"), *d.glob("*.sh")]
    assert files
    for f in files:
        assert "t3_ext" not in f.read_text(errors="replace"), f
