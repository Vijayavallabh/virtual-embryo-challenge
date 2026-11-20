"""Stage 6 of experiments/t3_ext (campaign/t3-ext-archives): series-level files (one archive or a loose 10x triple of a
series' OWN supplementary list, pinned to the samples whose cells it holds, with a split rule) and a control from a
sister series of the same study (screened, shown, hashed and pinned like the dataset's own series; d fixed at 0.5).
Synthetic SOFT text, synthetic counts and a stub download only: nothing here touches the network or real data."""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import re
import shutil
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
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(f"t3_ext_arch_{name}", EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


fe, pr, co = _load("fetch"), _load("process"), _load("combine")

ACC, SIS = "GSE910001", "GSE910002"
HOST = "ftp://ftp.ncbi.nlm.nih.gov/geo"
PIN = "    pinned: []  # GSE910001\n"
PLACEHOLDER_CS = """    control_series:
      accession: null
      screen_series: []
      relation: {same_lab: "TO FILL", same_protocol: "TO FILL", same_publication: "TO FILL"}
      reader_note: "TO FILL"
"""
NAMED_CS = """    control_series:
      accession: GSE910002
      screen_series: [GSE910002]
      relation: {same_lab: "same contact and institute on both series records",
                 same_protocol: "same extract and library protocol text, 10x 3' v3",
                 same_publication: "both series cite the same publication record"}
      reader_note: "both readers read both series' blocks and every sample of both tables"
"""
DS = """version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
  - accession: GSE910001
    status: cleared
    tier: unrelated
    study: GSE910001
    design: germline
    control_relation: littermate
    repository: geo
    screen_series: [GSE910001]
    stages_used: [E8.5]
    sample_file_patterns: ['(matrix\\.mtx|barcodes\\.tsv|features\\.tsv)(\\.gz)?$']
    series_files:
      - series: GSE910001
        patterns: ['^GSE910001_[A-Za-z0-9]+_E8\\.5_[A-Za-z0-9]+\\.tar\\.gz$',
                   '^GSE910001_[A-Za-z0-9]+_(matrix\\.mtx|barcodes\\.tsv|features\\.tsv)\\.gz$']
CONTROL_SERIES
""" + PIN

ARCHIVES = ["GSE910001_AAA_E8.5_1ab.tar.gz", "GSE910001_CCC_E8.5_1ab.tar.gz", "GSE910001_DDD_E8.5_1.tar.gz"]
LOOSE = ["GSE910001_BBB_barcodes.tsv.gz", "GSE910001_BBB_features.tsv.gz", "GSE910001_BBB_matrix.mtx.gz"]
OTHER_SERIES_FILES = ["GSE910001_RAW.tar", "GSE910001_notes.txt"]


def _yaml(tmp_path, named=False) -> Path:
    p = tmp_path / "datasets.yaml"
    p.write_text(DS.replace("CONTROL_SERIES\n", NAMED_CS if named else PLACEHOLDER_CS))
    return p


def _sample(gsm, title, geno, series, own_files=False):
    stub = "GSM" + gsm[3:-3] + "nnn"
    rec = {"title": title, "source_name_ch1": "whole embryo",
           "characteristics_ch1": [f"genotype: {geno}", "developmental stage: E8.5"],
           "extract_protocol_ch1": "Chromium Single Cell 3' v3", "series_id": series}
    if own_files:
        for i, f in enumerate(("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")):
            rec[f"supplementary_file_{i + 1}"] = f"{HOST}/samples/{stub}/{gsm}/suppl/{gsm}_{f}"
    else:
        rec["supplementary_file_1"] = "NONE"
    return rec


def _soft(acc, samples, files) -> str:
    stub = acc[:-3] + "nnn"
    lines = [f"^SERIES = {acc}", f"!Series_title = A synthetic series {acc}", "!Series_summary = embryos",
             "!Series_overall_design = one embryo pool per sample", "!Series_contact_name = A,,Person",
             # stage 7: a sister series is accepted only when both records share the institute and a contact
             "!Series_contact_institute = A Synthetic Institute"]
    lines += [f"!Series_supplementary_file = {HOST}/series/{stub}/{acc}/suppl/{f}" for f in files]
    for gsm, rec in samples.items():
        lines.append(f"^SAMPLE = {gsm}")
        for k, vals in rec.items():
            for v in (vals if isinstance(vals, list) else [vals]):
                lines.append(f"!Sample_{k} = {v}")
    return "\n".join(lines) + "\n"


OWN = {"GSM9101": _sample("GSM9101", "Aaa KO E8.5 1a", "Aaa -/-", ACC),
       "GSM9102": _sample("GSM9102", "Aaa KO E8.5 1b", "Aaa -/-", ACC),
       "GSM9103": _sample("GSM9103", "Bbb KO E8.5 1", "Bbb -/-", ACC),
       "GSM9104": _sample("GSM9104", "Bbb KO E8.5 2", "Bbb -/-", ACC),
       "GSM9106": _sample("GSM9106", "Ccc KO E8.5 1a", "Ccc -/-", ACC),
       "GSM9107": _sample("GSM9107", "Ccc KO E8.5 1b", "Ccc -/-", ACC),
       "GSM9105": _sample("GSM9105", "WT E8.5", "wild type", ACC, own_files=True)}
SISTER = {"GSM9201": _sample("GSM9201", "WT E8.5 1", "wild type", SIS, own_files=True),
          "GSM9202": _sample("GSM9202", "WT E8.5 2", "wild type", SIS, own_files=True)}


def _metas(named=False, own=None, sister=None, files=None, sister_title=None):
    m = {ACC: fe.parse_soft(_soft(ACC, own or OWN, files or ARCHIVES + LOOSE + OTHER_SERIES_FILES))}
    if named:
        t = _soft(SIS, sister or SISTER, [])
        if sister_title:
            t = t.replace(f"A synthetic series {SIS}", sister_title)
        m[SIS] = fe.parse_soft(t)
    return m


def _plan(y, dest, named=False, **kw):
    return fe.plan(ACC, dest, y, metas=_metas(named, **kw))


def _paste(y, block):
    y.write_text(y.read_text().replace(PIN, block))


# the pins the readers would make for the end-to-end case (own controls; the placeholder control series)
def _pin_all(y, dest, sister=False):
    ctl = ["GSM9201"] if sister else ["GSM9105"]
    block = "".join([
        fe.pin_template_from_series_file(ACC, dest, "ko", "Aaa", ["GSM9101", "GSM9102"], datasets=y),
        fe.pin_template_from_series_file(ACC, dest, "ko", "Bbb", ["GSM9103", "GSM9104"], datasets=y),
        fe.pin_template_from_series_file(ACC, dest, "ko", "Ccc", ["GSM9106", "GSM9107"], datasets=y),
        *[fe.pin_template(ACC, dest, "control", g, ctl, datasets=y) for g in ("Aaa", "Bbb", "Ccc")]])
    units = "".join([
        fe.series_file_template(ACC, dest, ACC, [ARCHIVES[0]], ["GSM9101", "GSM9102"], "folder",
                                ["AAA_1a=GSM9101", "AAA_1b=GSM9102"], datasets=y),
        fe.series_file_template(ACC, dest, ACC, LOOSE, ["GSM9103", "GSM9104"], "barcode_suffix",
                                ["-1=GSM9103", "-2=GSM9104"], datasets=y),
        fe.series_file_template(ACC, dest, ACC, [ARCHIVES[1]], ["GSM9106", "GSM9107"], "none", datasets=y)])
    series = fe.series_pin_template(ACC, dest, datasets=y)
    _paste(y, "    pinned_series:\n" + series + "    pinned_series_files:\n" + units + "    pinned:\n" + block)


# ---------------------------------------------------------------------------------------------------------------------
# synthetic counts

PANEL = [f"g{i}" for i in range(40)]
GROUPS = ["CM", "Endoderm", "Neural"]


def _triple(rng, mu, fold, n=160, suffix="", tag="") -> tuple[bytes, bytes, bytes]:
    extra = 260
    rows, names = [], []
    for g in ("CM", "Neural"):
        p = mu[g] * fold
        rows.append(np.hstack([rng.multinomial(2000, p / p.sum(), size=n), rng.poisson(2.0, size=(n, extra))]))
        names += [f"{g}_{tag}{i}{suffix}" for i in range(n)]
    M = np.vstack(rows).T
    b = io.BytesIO()
    spio.mmwrite(b, sparse.coo_matrix(M))
    feats = [f"E{i}\t{s}\tGene Expression" for i, s in enumerate(PANEL + [f"x{j}" for j in range(extra)])]
    return (gzip.compress(b.getvalue()), gzip.compress("\n".join(feats).encode()),
            gzip.compress("\n".join(names).encode()))


def _join(parts):
    """Two triples' bytes as one triple (barcodes concatenated): a loose aggregated matrix of two libraries."""
    mats, feats, bcs = [], None, []
    for m, f, b in parts:
        mats.append(sparse.csr_matrix(spio.mmread(io.BytesIO(gzip.decompress(m)), spmatrix=False)))
        feats = f
        bcs += gzip.decompress(b).decode().split("\n")
    b = io.BytesIO()
    spio.mmwrite(b, sparse.coo_matrix(sparse.hstack(mats)))
    return gzip.compress(b.getvalue()), feats, gzip.compress("\n".join(bcs).encode())


def _tar(path: Path, members: dict[str, bytes]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as t:
        for name, data in members.items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            t.addfile(ti, io.BytesIO(data))
    return path


def _stage(staging: Path, rng, mu):
    fold = np.ones(len(PANEL))
    fold[:4] = 3.0
    a1, a2 = _triple(rng, mu, fold, tag="a"), _triple(rng, mu, fold, tag="b")
    _tar(staging / ARCHIVES[0], {f"AAA_1{x}/{k}": v for x, t in (("a", a1), ("b", a2))
                                 for k, v in zip(("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"), t,
                                                 strict=True)})
    # stage 7: a file used whole holds as many libraries as its entry names samples (two folders, two samples)
    c1, c2 = _triple(rng, mu, fold, tag="a"), _triple(rng, mu, fold, tag="b")
    _tar(staging / ARCHIVES[1], {f"CCC_1{x}/{k}": v for x, t in (("a", c1), ("b", c2))
                                 for k, v in zip(("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"), t,
                                                 strict=True)})
    _tar(staging / ARCHIVES[2], {"DDD/matrix.mtx.gz": b"never read"})
    b = _join([_triple(rng, mu, fold, suffix="-1"), _triple(rng, mu, fold, suffix="-2")])
    for name, data in zip((LOOSE[2], LOOSE[1], LOOSE[0]), b, strict=True):
        (staging / name).write_bytes(data)
    for gsm in ("GSM9105", "GSM9201", "GSM9202"):
        for name, data in zip(("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"),
                              _triple(rng, mu, np.ones(len(PANEL))), strict=True):
            (staging / f"{gsm}_{name}").write_bytes(data)


def _fake_get(staging: Path, calls: list):
    def get(url, dst, max_bytes):
        calls.append(url)
        dst = cm.outside_repo(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(staging / Path(url).name, dst)
        return dst.stat().st_size, cm.sha256_file(dst)
    return get


def _clf(wt_path, C, present, genes, names, threads):
    calls = np.array([n.rsplit(":", 1)[1].split("/")[-1].split("_")[0] for n in names])
    return calls, np.full(len(names), 0.9), {"stub": True}


def _run_setup(tmp_path, monkeypatch, sister=False):
    rng = np.random.default_rng(5)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    wt = tmp_path / "wt.h5ad"
    wt.write_bytes(b"stand-in for the WT reference")
    refp = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats({"genes": PANEL, "groups": GROUPS, "shares": {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5},
                       "mu": mu, "count": {"CM": 20, "Endoderm": 30, "Neural": 50}}, refp, cm.sha256_file(wt))
    y = _yaml(tmp_path, named=sister)
    dest = tmp_path / "dest"
    assert _plan(y, dest, named=sister)["status"] == "ok"
    _pin_all(y, dest, sister=sister)
    assert _plan(y, dest, named=sister)["status"] == "ok"
    _stage(tmp_path / "staging", rng, mu)
    calls: list = []
    monkeypatch.setattr(fe, "_get", _fake_get(tmp_path / "staging", calls))
    return y, dest, wt, refp, calls


# ---------------------------------------------------------------------------------------------------------------------
# the declarations' form


def test_the_placeholder_control_series_is_inert_and_series_files_are_declared_per_screened_series(tmp_path):
    y = _yaml(tmp_path)
    ds = cm.dataset(ACC, y)
    assert cm.control_series(ds) is None and cm.all_screened(ds) == [ACC] and ds["d"] == 1.0
    assert ds["series_files"][0]["series"] == ACC
    bad = y.read_text().replace("      - series: GSE910001\n        patterns",
                                "      - series: GSE910009\n        patterns")
    y.write_text(bad)
    with pytest.raises(cm.Refused, match="series the dataset screens"):
        cm.load_datasets(y)
    y.write_text(DS.replace("CONTROL_SERIES\n", PLACEHOLDER_CS.replace("screen_series: []",
                                                                       "screen_series: [GSE910002]")))
    with pytest.raises(cm.Refused, match="placeholder"):
        cm.load_datasets(y)


@pytest.mark.parametrize("edit,match", [
    (lambda t: t.replace("reader_note: \"both readers", "reader_note: \"TO FILL both readers"), "reader_note"),
    (lambda t: t.replace("same_lab: \"same contact", "same_lab: \"TO FILL"), "relation.same_lab"),
    (lambda t: t.replace("same_publication:", "same_paper:"), "must state exactly"),
    (lambda t: t.replace("accession: GSE910002\n      screen_series: [GSE910002]",
                         "accession: GSE910001\n      screen_series: [GSE910001]"), "other than"),
    (lambda t: t.replace("screen_series: [GSE910002]", "screen_series: [GSE910003, GSE910002]"), "first"),
    (lambda t: t.replace("accession: GSE910002\n      screen_series: [GSE910002]",
                         "accession: GSE280587\n      screen_series: [GSE280587]"), "held or excluded"),
    (lambda t: t.replace("accession: GSE910002\n      screen_series: [GSE910002]",
                         "accession: GSE165300\n      screen_series: [GSE165300]"), "held or excluded"),
])
def test_a_named_control_series_needs_its_relation_note_and_no_held_series(tmp_path, edit, match):
    y = _yaml(tmp_path, named=True)
    ds = cm.dataset(ACC, y)
    assert cm.all_screened(ds) == [ACC, SIS] and cm.sister_series(ds) == [SIS]
    y.write_text(edit(y.read_text()))
    with pytest.raises(cm.Refused, match=match):
        cm.load_datasets(y)


def test_a_held_datasets_sister_series_is_held_too():
    allds = {"GSE1": {"status": "held_phenocopy", "screen_series": ["GSE1"],
                      "control_series": {"accession": "GSE2", "screen_series": ["GSE2", "GSE3"]}}}
    assert {"GSE1", "GSE2", "GSE3"} <= cm.held_series(allds)


def _pinned_yaml(tmp_path, named=False):
    y = _yaml(tmp_path, named=named)
    dest = tmp_path / "d"
    _plan(y, dest, named=named)
    _pin_all(y, dest, sister=named)
    return y, dest


UNIT_EDITS = [
    ("series: GSE910001\n        files: [\"GSE910001_AAA", "series: GSE910009\n        files: [\"GSE910001_AAA",
     "does not declare"),
    ("split: {by: none}", "split: {by: none, map: {\"x\": GSM9106}}", "takes no map"),
    ("split: {by: folder, map: {\"AAA_1a\": GSM9101, \"AAA_1b\": GSM9102}}",
     "split: {by: folder, map: {\"AAA_1a\": GSM9101}}", "every sample"),
    ("split: {by: folder, map: {\"AAA_1a\": GSM9101, \"AAA_1b\": GSM9102}}",
     "split: {by: folder, map: {\"../up\": GSM9101, \"AAA_1b\": GSM9102}}", "plain"),
    ("split: {by: barcode_suffix, map: {\"-1\": GSM9103, \"-2\": GSM9104}}",
     "split: {by: barcode_suffix, map: {\"1\": GSM9103, \"-2\": GSM9104}}", "plain"),
    ("split: {by: barcode_suffix, map: {\"-1\": GSM9103, \"-2\": GSM9104}}",
     "split: {by: folder, map: {\"a\": GSM9103, \"b\": GSM9104}}", "needs an archive"),
    ("split: {by: barcode_suffix", "split: {by: regex", "split.by"),
    ("samples: [GSM9106, GSM9107]", "samples: [GSM9106, GSM9107, GSM9101]", "two series file entries"),
    ("samples: [GSM9106, GSM9107]", "samples: [GSM9106, GSM9107, GSM9999]", "not a pinned sample"),
    ("files: [\"GSE910001_CCC_E8.5_1ab.tar.gz\"]", "files: [\"GSE910001_AAA_E8.5_1ab.tar.gz\"]",
     "two series file entries"),
    ("files: [\"GSE910001_CCC_E8.5_1ab.tar.gz\"]",
     "files: [\"GSE910001_CCC_E8.5_1ab.tar.gz\", \"GSE910001_DDD_E8.5_1.tar.gz\"]", "one archive"),
    ("files: [\"GSE910001_CCC_E8.5_1ab.tar.gz\"]", "files: [\"../GSE910001_CCC.tar.gz\"]", "plain file names"),
]


@pytest.mark.parametrize("old,new,match", UNIT_EDITS)
def test_the_series_file_pins_form_is_checked(tmp_path, old, new, match):
    y, _ = _pinned_yaml(tmp_path)
    assert cm.dataset(ACC, y)["pinned_series_files"]
    t = y.read_text()
    assert old in t
    y.write_text(t.replace(old, new, 1))
    with pytest.raises(cm.Refused, match=match):
        cm.load_datasets(y)


def test_a_unit_used_whole_holds_one_arm_gene_and_stage_and_a_sample_has_one_source(tmp_path):
    y, _ = _pinned_yaml(tmp_path)
    t = y.read_text()
    # GSM9107 re-pinned as gene Aaa: the CCC archive used as one unit would mix two genes
    i = t.index("- accession: GSM9107")
    y.write_text(t[:i] + t[i:].replace('gene: "Ccc"', 'gene: "Aaa"', 1))
    with pytest.raises(cm.Refused, match="ONE arm and stage, each pinned for the same genes"):
        cm.load_datasets(y)
    # a sample named by a series file but pinned with its own files too
    i = t.index("- accession: GSM9101")
    y.write_text(t[:i] + t[i:].replace("files: []", 'files: ["GSM9101_matrix.mtx.gz"]', 1))
    with pytest.raises(cm.Refused, match="files: \\[\\]"):
        cm.load_datasets(y)
    # a sample with no files that no series file names
    i = t.index("- accession: GSM9105")
    j = t.index("files: [", i)
    k = t.index("]", j)
    y.write_text(t[:j] + "files: []" + t[k + 1:])
    with pytest.raises(cm.Refused, match="non-empty list"):
        cm.load_datasets(y)


# ---------------------------------------------------------------------------------------------------------------------
# the series' own folder, and the plan table


@pytest.mark.parametrize("url,ok", [
    (f"{HOST}/series/GSE910nnn/GSE910001/suppl/GSE910001_AAA_E8.5_1ab.tar.gz", True),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE910nnn/GSE910001/suppl/GSE910001_x.tar.gz", True),
    (f"{HOST}/series/GSE910nnn/GSE910002/suppl/GSE910002_AAA.tar.gz", False),        # another series' folder
    (f"{HOST}/series/GSE910nnn/GSE910001/suppl/GSE910002_AAA.tar.gz", False),        # another series' name
    (f"{HOST}/series/GSE910nnn/GSE910001/suppl/AAA.tar.gz", False),                  # no accession prefix
    (f"{HOST}/series/GSE91nnn/GSE910001/suppl/GSE910001_AAA.tar.gz", False),         # a wrong stub
    (f"{HOST}/samples/GSM9101nnn/GSM9101/suppl/GSE910001_AAA.tar.gz", False),        # a sample folder
    (f"{HOST}/series/GSE910nnn/GSE910001/suppl/GSE910001_AAA.tar.gz?x=1", False),    # a query
    ("https://example.org/geo/series/GSE910nnn/GSE910001/suppl/GSE910001_AAA.tar.gz", False),
    ("ftp://ftp.ncbi.nlm.nih.gov:21/geo/series/GSE910nnn/GSE910001/suppl/GSE910001_AAA.tar.gz", False),
])
def test_a_series_file_is_taken_only_from_the_series_own_folder(url, ok):
    u, why = fe.own_series_file(ACC, url)
    assert (u is not None) == ok, why
    assert fe.series_stub("GSE123") == "GSEnnn" and fe.series_stub("GSE137337") == "GSE137nnn"


def test_the_plan_table_shows_the_series_files_offered_and_never_fetched_with_reasons(tmp_path):
    y = _yaml(tmp_path)
    files = ARCHIVES + LOOSE + OTHER_SERIES_FILES + ["GSE910001_ZZZ_E8.5_GSE165300.tar.gz"]
    pl = _plan(y, tmp_path / "d", files=files)
    assert pl["status"] == "ok"
    b = pl["series"][ACC]
    assert [f["name"] for f in b["files"]] == ARCHIVES + LOOSE
    reasons = {o["name"]: o["reason"] for o in b["other_files"]}
    assert reasons["GSE910001_RAW.tar"] == "a denied file type"
    assert "patterns" in reasons["GSE910001_notes.txt"]
    assert "held or excluded" in reasons["GSE910001_ZZZ_E8.5_GSE165300.tar.gz"]
    text = (tmp_path / "d" / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert ("series files allowed (this series' own folder, series_files patterns): "
            "['GSE910001_AAA_E8.5_1ab.tar.gz'") in text
    assert "GSE910001_RAW.tar (a denied file type)" in text


def test_pin_templates_for_series_files_print_blocks_that_load(tmp_path):
    y, dest = _pinned_yaml(tmp_path)
    ds = cm.dataset(ACC, y)
    assert len(ds["pinned_series_files"]) == 3
    assert {p["accession"]: p["files"] for p in ds["pinned"]}["GSM9101"] == []
    with pytest.raises(cm.Refused, match="offers"):
        fe.series_file_template(ACC, dest, ACC, ["GSE910001_RAW.tar"], ["GSM9101"], "none", datasets=y)
    with pytest.raises(cm.Refused, match="offers"):
        fe.series_file_template(ACC, dest, ACC, ["GSE910001_notes.txt"], ["GSM9101"], "none", datasets=y)
    with pytest.raises(cm.Refused, match="not in the plan table"):
        fe.series_file_template(ACC, dest, ACC, [ARCHIVES[2]], ["GSM9999"], "none", datasets=y)
    with pytest.raises(cm.Refused, match="map every named sample"):
        fe.series_file_template(ACC, dest, ACC, [ARCHIVES[2]], ["GSM9101", "GSM9102"], "folder", ["a=GSM9101"],
                                datasets=y)
    with pytest.raises(cm.Refused, match="no --key"):
        fe.series_file_template(ACC, dest, ACC, [ARCHIVES[2]], ["GSM9101"], "none", ["a=GSM9101"], datasets=y)
    with pytest.raises(cm.Refused, match="no per-sample file"):
        fe.pin_template(ACC, dest, "ko", "Aaa", ["GSM9101"], datasets=y)


# ---------------------------------------------------------------------------------------------------------------------
# fetch and process


def test_fetch_downloads_only_the_pinned_series_files_and_pinned_sample_files(tmp_path, monkeypatch):
    y, dest, wt, refp, calls = _run_setup(tmp_path, monkeypatch)
    man = fe.fetch(ACC, dest, y)
    names = sorted(Path(u).name for u in calls)
    assert names == sorted([ARCHIVES[0], ARCHIVES[1], *LOOSE, "GSM9105_barcodes.tsv.gz", "GSM9105_features.tsv.gz",
                            "GSM9105_matrix.mtx.gz"])
    assert all(re.search(r"/geo/(series/GSE910nnn/GSE910001|samples/GSM9nnn/GSM9105)/suppl/", u) for u in calls)
    assert ARCHIVES[2] not in names and "GSE910001_RAW.tar" not in names
    rec = man["files"][f"{ACC}/{ARCHIVES[0]}"]
    assert rec["series_metadata_sha256"] == cm.series_metadata_sha256(ACC, _metas()[ACC]["series"])
    assert rec["pins_sha256"] == cm.pins_sha256(cm.dataset(ACC, y))
    # nothing is unpacked: the series folder holds the pinned files only
    assert sorted(p.name for p in (dest / "raw" / ACC / ACC).iterdir()) == sorted([ARCHIVES[0], ARCHIVES[1], *LOOSE])
    # a stray file in the series folder goes; a second fetch reuses the rest
    (dest / "raw" / ACC / ACC / ARCHIVES[2]).write_bytes(b"stray")
    calls.clear()
    fe.fetch(ACC, dest, y)
    assert calls == [] and not (dest / "raw" / ACC / ACC / ARCHIVES[2]).exists()


def test_a_hand_edited_series_file_url_in_the_stored_plan_is_never_fetched(tmp_path, monkeypatch):
    y, dest, wt, refp, calls = _run_setup(tmp_path, monkeypatch)
    pp = dest / "meta" / ACC / "plan.json"
    pl = json.loads(pp.read_text())
    pl["series"][ACC]["files"][0]["url"] = pl["series"][ACC]["files"][0]["url"].replace(
        "GSE910001/suppl/GSE910001_AAA", "GSE910001/suppl/GSE910001_DDD")
    pp.write_text(json.dumps(pl))
    with pytest.raises(cm.Refused, match="series file list differs"):
        fe.fetch(ACC, dest, y)
    assert calls == []


def test_end_to_end_archives_split_by_folder_suffix_and_whole(tmp_path, monkeypatch):
    y, dest, wt, refp, calls = _run_setup(tmp_path, monkeypatch)
    fe.fetch(ACC, dest, y)
    before = sorted(p.relative_to(tmp_path) for p in tmp_path.rglob("*"))
    out = tmp_path / "processed"
    qc = pr.run(ACC, dest, wt, refp, out, y, classifier=_clf)
    assert qc["status"] == "ok"
    # nothing was unpacked anywhere: the only new files are the outputs
    after = sorted(p.relative_to(tmp_path) for p in tmp_path.rglob("*"))
    assert all(str(p).startswith("processed") for p in set(after) - set(before))
    units = {u["files"][0]: u for u in qc["series_file_units"]}
    assert [lib["label"] for lib in units[ARCHIVES[0]]["libraries"]] == ["GSM9101", "GSM9102"]
    assert [lib["label"] for lib in units[LOOSE[0]]["libraries"]] == ["GSM9103", "GSM9104"]
    whole = f"series:{ACC}/{ARCHIVES[1]}"
    assert [(lib["label"], lib["folder"]) for lib in units[ARCHIVES[1]]["libraries"]] == [(whole, "CCC_1a"),
                                                                                           (whole, "CCC_1b")]
    # the whole archive counts as ONE library: one confidence record, m = 1 for its knockout arm
    assert set(qc["confidence_by_sample"]) == {"GSM9101", "GSM9102", "GSM9103", "GSM9104", "GSM9105", whole}
    for gene in ("Aaa", "Bbb", "Ccc"):
        u = qc["units"][gene]
        assert u["status"] == "ok", (gene, u)
        assert u["cells_per_stage_and_arm"] == {"E8.5:control": 320, "E8.5:ko": 640}
    assert qc["units"]["Ccc"]["stages"]["E8.5"]["groups"]["CM"]["m"] == 1
    assert "control_series" not in qc and "d" not in qc["units"]["Aaa"]
    assert set(qc["count_files_sha256"]) >= {f"{ACC}/{ARCHIVES[0]}", f"{ACC}/{LOOSE[2]}",
                                             "GSM9105/GSM9105_matrix.mtx.gz"}
    r = co.combine(out, refp, tmp_path / "tables", y, dest=dest)
    assert "study_GSE910001" in r["tables"]


def test_process_refuses_a_stray_or_rebound_file_in_the_series_folder(tmp_path, monkeypatch):
    y, dest, wt, refp, calls = _run_setup(tmp_path, monkeypatch)
    fe.fetch(ACC, dest, y)
    stray = dest / "raw" / ACC / ACC / ARCHIVES[2]
    stray.write_bytes(b"stray")
    with pytest.raises(cm.Refused, match="not a pinned series file"):
        pr.run(ACC, dest, wt, refp, tmp_path / "o", y, classifier=_clf)
    stray.unlink()
    mp = dest / "raw" / ACC / "MANIFEST.json"
    man = json.loads(mp.read_text())
    man["files"][f"{ACC}/{ARCHIVES[0]}"]["series_metadata_sha256"] = "0" * 64
    mp.write_text(json.dumps(man))
    with pytest.raises(cm.Refused, match="another series text"):
        pr.run(ACC, dest, wt, refp, tmp_path / "o", y, classifier=_clf)
    assert not (tmp_path / "o" / ACC).exists()


# ---------------------------------------------------------------------------------------------------------------------
# what process refuses inside an archive


def _unit(split, samples=("GSM1", "GSM2")):
    return {"series": "GSE1", "files": ["GSE1_x.tar.gz"], "samples": list(samples), "split": split}


FOLDER = {"by": "folder", "map": {"A": "GSM1", "B": "GSM2"}}


def _good_members():
    rng = np.random.default_rng(1)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    t1, t2 = _triple(rng, mu, np.ones(len(PANEL)), n=3), _triple(rng, mu, np.ones(len(PANEL)), n=3)
    kinds = ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
    return {**{f"A/{k}": v for k, v in zip(kinds, t1, strict=True)},
            **{f"B/{k}": v for k, v in zip(kinds, t2, strict=True)}}


def _bad_tar(path, kind):
    path.parent.mkdir(parents=True, exist_ok=True)
    good = _good_members()
    with tarfile.open(path, "w:gz") as t:
        for name, data in good.items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            t.addfile(ti, io.BytesIO(data))
        ti = tarfile.TarInfo({"traversal": "A/../../evil/matrix.mtx.gz", "absolute": "/tmp/evil/matrix.mtx.gz",
                              "symlink": "A/link", "hardlink": "A/hard", "fifo": "A/pipe",
                              "unknown": "A/README.txt", "duplicate": "A/matrix.mtx.gz",
                              "screened": "A/Gata4_counts_matrix.mtx.gz", "nofolder": "matrix.mtx.gz",
                              "otherfolder": "C/matrix.mtx.gz", "backslash": "A\\x\\matrix.mtx.gz",
                              "dotdir": "A/./barcodes.tsv.gz", "drive": "C:/evil/matrix.mtx.gz",
                              "denied": "A/raw_feature_bc_matrix/matrix.mtx.gz"}[kind])
        data = b"x"
        if kind == "symlink":
            ti.type, ti.linkname = tarfile.SYMTYPE, "/etc/passwd"
        elif kind == "hardlink":
            ti.type, ti.linkname = tarfile.LNKTYPE, "A/matrix.mtx.gz"
        elif kind == "fifo":
            ti.type = tarfile.FIFOTYPE
        if ti.isfile():
            ti.size = len(data)
            t.addfile(ti, io.BytesIO(data))
        else:
            t.addfile(ti)
    return path


@pytest.mark.parametrize("kind,match", [
    ("traversal", "path traversal"), ("absolute", "absolute"), ("drive", "absolute"), ("backslash", "absolute"),
    ("symlink", "symlink"), ("hardlink", "hard link"), ("fifo", "not a regular file"),
    ("unknown", "no rule"), ("duplicate", "twice"), ("screened", "held-out screen"),
    ("nofolder", "no folder the split map"), ("otherfolder", "no folder the split map"), ("dotdir", "'.' path"),
    ("denied", "denied type"),
])
def test_an_archive_member_that_breaks_a_rule_refuses_the_file_and_nothing_is_written(tmp_path, kind, match):
    p = _bad_tar(tmp_path / "in" / "GSE1_x.tar.gz", kind)
    before = sorted(tmp_path.rglob("*"))
    with pytest.raises(cm.Refused, match=match):
        pr.read_series_unit([p], _unit(FOLDER))
    assert sorted(tmp_path.rglob("*")) == before
    assert not Path("/tmp/evil").exists()


def test_a_good_archive_splits_by_folder_and_reads_in_memory(tmp_path):
    p = _tar(tmp_path / "GSE1_x.tar.gz", _good_members())
    libs = pr.read_series_unit([p], _unit(FOLDER))
    assert [(g, f, len(b)) for g, f, b, _ in libs] == [("GSM1", "A", 6), ("GSM2", "B", 6)]
    assert sorted(x.name for x in tmp_path.iterdir()) == ["GSE1_x.tar.gz"]
    # a folder the map names but the archive lacks, and a library missing a member
    with pytest.raises(cm.Refused, match="hold no library"):
        pr.read_series_unit([p], _unit({"by": "folder", "map": {"A": "GSM1", "B": "GSM2", "Z": "GSM2"}}))
    m = {k: v for k, v in _good_members().items() if k != "B/barcodes.tsv.gz"}
    with pytest.raises(cm.Refused, match="exactly one matrix"):
        pr.read_series_unit([_tar(tmp_path / "GSE1_y.tar.gz", m)], _unit(FOLDER))
    # split none: the file is one unit; each folder one library, both under the entry's first sample
    whole = pr.read_series_unit([p], _unit({"by": "none"}))
    assert [(g, f, len(b)) for g, f, b, _ in whole] == [("GSM1", "A", 6), ("GSM1", "B", 6)]
    m = {**_good_members(), "B/x_matrix.mtx.gz": b"y"}
    with pytest.raises(cm.Refused, match="two 'matrix'"):
        pr.read_series_unit([_tar(tmp_path / "GSE1_z.tar.gz", m)], _unit({"by": "none"}))


def test_a_barcode_that_no_suffix_fits_refuses_and_a_suffix_with_no_barcode_refuses(tmp_path):
    rng = np.random.default_rng(2)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    one = np.ones(len(PANEL))
    t = _join([_triple(rng, mu, one, n=3, suffix="-1"), _triple(rng, mu, one, n=3, suffix="-2")])
    names = ["GSE1_m_matrix.mtx.gz", "GSE1_m_features.tsv.gz", "GSE1_m_barcodes.tsv.gz"]
    for n, d in zip(names, t, strict=True):
        (tmp_path / n).write_bytes(d)
    paths = [tmp_path / n for n in names]
    u = {"series": "GSE1", "files": names, "samples": ["GSM1", "GSM2"],
         "split": {"by": "barcode_suffix", "map": {"-1": "GSM1", "-2": "GSM2"}}}
    libs = pr.read_series_unit(paths, u)
    assert [(g, len(b)) for g, _, b, _ in libs] == [("GSM1", 6), ("GSM2", 6)]
    with pytest.raises(cm.Refused, match="no split key"):
        pr.read_series_unit(paths, {**u, "samples": ["GSM1"], "split": {"by": "barcode_suffix",
                                                                          "map": {"-1": "GSM1"}}})
    with pytest.raises(cm.Refused, match="no barcode carries"):
        pr.read_series_unit(paths, {**u, "samples": ["GSM1", "GSM2", "GSM3"],
                                    "split": {"by": "barcode_suffix", "map": {"-1": "GSM1", "-2": "GSM2",
                                                                              "-3": "GSM3"}}})


# ---------------------------------------------------------------------------------------------------------------------
# a control from a sister series


def test_a_sister_series_is_planned_shown_hashed_and_pinned_like_the_own_series(tmp_path):
    y = _yaml(tmp_path, named=True)
    with pytest.raises(cm.Refused, match="not exactly the screened series"):
        fe.plan(ACC, tmp_path / "d", y, metas=_metas(named=False))
    pl = _plan(y, tmp_path / "d", named=True)
    assert pl["status"] == "ok" and set(pl["series"]) == {ACC, SIS} and set(pl["tables"]) == {ACC, SIS}
    assert pl["series"][SIS]["role"] == "control series" and pl["control_series"]["accession"] == SIS
    own = (tmp_path / "d" / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    sis = (tmp_path / "d" / "meta" / ACC / f"plan_table_{SIS}.txt").read_text()
    for text in (own, sis):
        assert f"==== SERIES {SIS}" in text and f"==== SERIES {ACC}" in text
        assert "# CONTROL SERIES (datasets.yaml control_series)" in text and "reader note: both readers" in text
        assert "relation same_publication: both series cite the same publication record" in text
    assert "THIS IS THE CONTROL SERIES" in sis and "==== GSM9201" in sis and "title: WT E8.5 1" in sis
    assert "series membership: ok" in sis
    # the sister's samples are offered as controls only
    with pytest.raises(cm.Refused, match="not in the plan table"):
        fe.pin_template(ACC, tmp_path / "d", "ko", "Aaa", ["GSM9201"], datasets=y)
    txt = fe.pin_template(ACC, tmp_path / "d", "control", "Aaa", ["GSM9201"], datasets=y)
    assert "CONTROL SERIES GSE910002" in txt
    assert [b.count("- accession") for b in [fe.series_pin_template(ACC, tmp_path / "d", datasets=y)]] == [2]


def test_a_sister_sample_pinned_as_a_knockout_or_a_changed_sister_text_is_refused(tmp_path):
    y, dest = _pinned_yaml(tmp_path, named=True)
    assert _plan(y, dest, named=True)["status"] == "ok"
    ds = cm.dataset(ACC, y)
    pl = json.loads((dest / "meta" / ACC / "plan.json").read_text())
    got = {s["accession"]: s["sister"] for s in fe.verify_pins(ds, pl)}
    assert got["GSM9201"] is True and got["GSM9101"] is False
    bad = [{**p, "arm": "ko"} if p["accession"] == "GSM9201" and p["gene"] == "Aaa" else p for p in ds["pinned"]]
    with pytest.raises(cm.Refused, match="control only"):
        fe.verify_pins({**ds, "pinned": bad}, pl)
    live = _plan(y, tmp_path / "d2", named=True, sister_title="A synthetic series, corrected")
    assert live["status"] == "refused" and f"series {SIS}" in live["reason"] and "live text" in live["reason"]
    # the control series declaration is part of what the plan was made under
    with pytest.raises(cm.Refused, match="another control_series"):
        fe.verify_pins({**ds, "control_series": {**ds["control_series"], "reader_note": "other"}}, pl)


def test_a_unit_with_a_sister_control_takes_d_one_half_and_says_so(tmp_path, monkeypatch):
    y, dest, wt, refp, calls = _run_setup(tmp_path, monkeypatch, sister=True)
    assert cm.dataset(ACC, y)["d"] == 1.0          # germline against littermates: d would be 1
    fe.fetch(ACC, dest, y)
    assert any("/samples/GSM9nnn/GSM9201/suppl/" in u for u in calls)
    assert not any("GSM9202" in u for u in calls)
    qc = pr.run(ACC, dest, wt, refp, tmp_path / "processed", y, classifier=_clf)
    assert qc["status"] == "ok" and qc["sister_series_controls"] == ["GSM9201"]
    assert qc["control_series"]["accession"] == SIS
    for gene in ("Aaa", "Bbb", "Ccc"):
        u = qc["units"][gene]
        assert u["d"] == cm.SISTER_D == 0.5
        assert u["control_source"].startswith(f"sister series {SIS}")
        g = u["stages"]["E8.5"]["groups"]["CM"]
        assert g["k"] == pytest.approx(cm.shrink_k(list(g["cells"].values()), g["m"], 0.5))
        unit = json.loads(Path(u["unit_file"]).read_text())
        assert unit["d"] == 0.5 and unit["control_source"] == u["control_source"]
    assert cm.unit_design_factor({"d": 1.0}, True) == 0.5 and cm.unit_design_factor({"d": 1.0}, False) == 1.0
    assert cm.unit_design_factor({"d": 0.5}, True) == 0.5


# ---------------------------------------------------------------------------------------------------------------------
# the committed file, the README and the launcher


def test_the_committed_gse137337_declarations_offer_its_gene_archives_and_pin_nothing():
    ds = cm.dataset("GSE137337")
    assert ds["pinned"] == [] and ds["pinned_series"] == [] and ds["pinned_series_files"] == []
    assert cm.control_series(ds) is None and ds["control_series"]["accession"] is None
    assert cm.all_screened(ds) == ["GSE137337"] and ds["d"] == 0.5
    host = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137337/suppl/"
    names = ["GSE137337_RNF2_E8.5_1ab.tar.gz", "GSE137337_EED_E8.5_2.tar.gz", "GSE137337_KDM2B_barcodes.tsv.gz",
             "GSE137337_KDM2B_features.tsv.gz", "GSE137337_KDM2B_matrix.mtx.gz", "GSE137337_RAW.tar",
             "GSE137337_RNF2_E7.5_1.tar.gz", "GSE137337_scBS_methylation_E8.5_x.tar.gz",
             "GSE137337_WGBS_matrix.mtx.gz", "GSE137337_notes.txt"]
    ok, other = fe.series_files(ds, "GSE137337", {"supplementary_file": [host + n for n in names]})
    assert [f["name"] for f in ok] == names[:5]
    assert {o["name"] for o in other} == set(names[5:])
    # every other cleared dataset is as it was: no series files, no control series
    for acc, d in cm.load_datasets().items():
        if acc != "GSE137337":
            assert not d.get("series_files") and not d.get("control_series") and not d["pinned_series_files"]


def test_pins_hash_is_unchanged_for_a_dataset_without_stage_6_declarations():
    ds = {"pinned": [{"accession": "GSM1"}], "pinned_series": [{"accession": "GSE1"}]}
    want = __import__("hashlib").sha256(json.dumps({"samples": ds["pinned"], "series": ds["pinned_series"]},
                                                   sort_keys=True, ensure_ascii=False, separators=(",", ":"))
                                        .encode()).hexdigest()
    assert cm.pins_sha256(ds) == want
    assert cm.pins_sha256({**ds, "pinned_series_files": [{"series": "GSE1"}]}) != want


def test_the_readme_and_launcher_describe_series_files_and_the_sister_series():
    text = (EXT / "README.md").read_text()
    for s in ("## Series-level files", "## A control from a sister series", "series-file-template",
              "--from-series-file", "d = 0.5"):
        assert s in text, s
    sh = (EXT / "run_dgx.sh").read_text()
    assert "series files" in sh and "control series" in sh
    y = yaml.safe_load((EXT / "datasets.yaml").read_text())
    assert [d["accession"] for d in y["datasets"] if "series_files" in d or "control_series" in d] == ["GSE137337"]


# ---------------------------------------------------------------------------------------------------------------------
# stage 7 (campaign/t3-ext-3): the three blockers of the two reviews of 853b214, each replaying the reviewer's
# reproduction, and the lesser points


def _lib(rng, n=3):
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    return dict(zip(("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"), _triple(rng, mu, np.ones(len(PANEL)), n=n),
                    strict=True))


def test_blocker_1_split_none_refuses_a_file_with_more_libraries_than_samples_and_writes_nothing(tmp_path):
    """The first reviewer's reproduction: a tar.gz with folders KO_1a, WT_littermate and unlisted_3 pinned split none
    to one knockout sample GSM100 returned three libraries under GSM100 (another arm's and an unpinned sample's
    cells in the knockout arm). Now it is refused before any count is read."""
    rng = np.random.default_rng(11)
    members = {f"{d}/{k}": v for d in ("KO_1a", "WT_littermate", "unlisted_3") for k, v in _lib(rng).items()}
    p = _tar(tmp_path / "in" / "GSE1_x.tar.gz", members)
    before = sorted(tmp_path.rglob("*"))
    with pytest.raises(cm.Refused, match="split none: the file holds 3 libraries but the entry names 1 samples"):
        pr.read_series_unit([p], _unit({"by": "none"}, samples=("GSM100",)))
    assert sorted(tmp_path.rglob("*")) == before
    # one sample, two folders (the reviewer's suggested test) and a loose triple named for two samples: refused too
    two = _tar(tmp_path / "in" / "GSE1_y.tar.gz", {f"{d}/{k}": v for d in ("A", "B") for k, v in _lib(rng).items()})
    with pytest.raises(cm.Refused, match="holds 2 libraries but the entry names 1 samples"):
        pr.read_series_unit([two], _unit({"by": "none"}, samples=("GSM1",)))
    names = ["GSE1_m_matrix.mtx.gz", "GSE1_m_features.tsv.gz", "GSE1_m_barcodes.tsv.gz"]
    for n, v in zip(names, _lib(rng).values(), strict=True):
        (tmp_path / n).write_bytes(v)
    loose = {"series": "GSE1", "files": names, "samples": ["GSM1", "GSM2"], "split": {"by": "none"}}
    with pytest.raises(cm.Refused, match="holds 1 libraries but the entry names 2 samples"):
        pr.read_series_unit([tmp_path / n for n in names], loose)
    # equal counts pass, every library under the entry's first sample, its folder recorded
    libs = pr.read_series_unit([two], _unit({"by": "none"}, samples=("GSM1", "GSM2")))
    assert [(g, f) for g, f, _, _ in libs] == [("GSM1", "A"), ("GSM1", "B")]


def test_blocker_1_the_barcode_cap_applies_to_each_library_under_split_none(tmp_path):
    """The reviewers' second finding: under split none the cap was MAX_BARCODES times the number of samples per
    folder. One folder of 60,001 barcodes beside a small one, two samples: refused now (120,000 was allowed)."""
    rng = np.random.default_rng(12)
    n_bc, n_ft = cm.MAX_BARCODES + 1, 300
    big = sparse.coo_matrix((np.ones(10), (np.arange(10), np.arange(10))), shape=(n_ft, n_bc))
    b = io.BytesIO()
    spio.mmwrite(b, big)
    feats = "\n".join(f"E{i}\tx{i}\tGene Expression" for i in range(n_ft)).encode()
    members = {"A/matrix.mtx.gz": gzip.compress(b.getvalue()), "A/features.tsv.gz": gzip.compress(feats),
               "A/barcodes.tsv.gz": gzip.compress("\n".join(f"B{i}-1" for i in range(n_bc)).encode()),
               **{f"B/{k}": v for k, v in _lib(rng).items()}}
    p = _tar(tmp_path / "GSE1_x.tar.gz", members)
    with pytest.raises(cm.Refused, match=f"{n_bc} barcodes in one library"):
        pr.read_series_unit([p], _unit({"by": "none"}))


def _committed_with_sister(acc, screen=None):
    t = cm.DATASETS.read_text()
    old = ('    control_series:\n      accession: null\n      screen_series: []\n      relation:\n'
           '        same_lab: "TO FILL"\n        same_protocol: "TO FILL"\n        same_publication: "TO FILL"\n'
           '      reader_note: "TO FILL"\n')
    assert t.count(old) == 1
    listed = ", ".join(screen or [acc])
    return t.replace(old, f"    control_series:\n      accession: {acc}\n      screen_series: [{listed}]\n"
                          '      relation: {same_lab: "yes", same_protocol: "yes", same_publication: "yes"}\n'
                          '      reader_note: "yes"\n')


@pytest.mark.parametrize("acc,screen", [("GSE204906", None), ("GSE204908", None),
                                        ("GSE930001", ["GSE930001", "GSE167493"])])
def test_blocker_2_a_sister_series_registered_or_screened_by_another_dataset_is_refused(tmp_path, acc, screen):
    """The first reviewer's reproduction: control_series.accession = GSE204906 (the Dnmt3a/b study, its own cleared
    dataset) with relation strings 'yes' passed. Now load_datasets refuses it, its superseries, and a sister whose
    screen_series names another dataset's superseries."""
    y = tmp_path / "datasets.yaml"
    y.write_text(_committed_with_sister(acc, screen))
    with pytest.raises(cm.Refused, match="is registered as \\(or screened by\\) dataset GSE"):
        cm.load_datasets(y)
    # an unregistered sister passes this check (it is still matched against the records at plan time)
    y.write_text(_committed_with_sister("GSE930001"))
    assert cm.control_series(cm.load_datasets(y)["GSE137337"])["accession"] == "GSE930001"


def _sister_metas(own_extra=(), sister_extra=(), drop=("contact_institute",), extra_series=None):
    """Own and sister SOFT with chosen contact lines: the sister's listed fields dropped, then lines added."""
    own = _soft(ACC, OWN, ARCHIVES + LOOSE + OTHER_SERIES_FILES).replace(
        "!Series_contact_name = A,,Person\n", "!Series_contact_name = A,,Person\n" + "".join(own_extra))
    sis = _soft(SIS, SISTER, [])
    for k in drop:
        sis = "\n".join(x for x in sis.split("\n") if not x.startswith(f"!Series_{k} =")) + "\n"
    sis = sis.replace("!Series_overall_design = one embryo pool per sample\n",
                      "!Series_overall_design = one embryo pool per sample\n" + "".join(sister_extra))
    m = {ACC: fe.parse_soft(own), SIS: fe.parse_soft(sis)}
    m.update(extra_series or {})
    return m


def test_blocker_2_a_sister_series_needs_a_shared_institute_and_a_shared_contact_or_publication(tmp_path):
    y = _yaml(tmp_path, named=True)
    # a different institute: refused, the reason stated
    m = _sister_metas(sister_extra=["!Series_contact_institute = Another Institute\n"])
    pl = fe.plan(ACC, tmp_path / "a", y, metas=m)
    assert pl["status"] == "refused" and "do not share a contact_institute" in pl["reason"]
    # the same institute but no other shared contact, laboratory or publication value: refused
    m = _sister_metas(sister_extra=["!Series_contact_institute = A Synthetic Institute\n"],
                      drop=("contact_institute", "contact_name"))
    pl = fe.plan(ACC, tmp_path / "b", y, metas=m)
    assert pl["status"] == "refused" and "nothing else" in pl["reason"]
    # case is kept: 'a synthetic institute' is not the same value
    m = _sister_metas(sister_extra=["!Series_contact_institute = a synthetic institute\n"])
    assert fe.plan(ACC, tmp_path / "c", y, metas=m)["status"] == "refused"
    # a shared institute and a shared PubMed id: passes, and the plan tables print what the code matched
    m = _sister_metas(own_extra=["!Series_pubmed_id = 123\n"],
                      sister_extra=["!Series_contact_institute = A Synthetic Institute\n", "!Series_pubmed_id = 123\n"],
                      drop=("contact_institute", "contact_name"))
    pl = fe.plan(ACC, tmp_path / "d", y, metas=m)
    assert pl["status"] == "ok"
    assert pl["control_series_match"] == {"contact_institute": ["A Synthetic Institute"], "pubmed_id": ["123"]}
    for s in (ACC, SIS):
        text = (tmp_path / "d" / "meta" / ACC / f"plan_table_{s}.txt").read_text()
        assert ("the code matched (shared by the own series' and the control series' records): contact_institute = "
                "A Synthetic Institute") in text and "pubmed_id = 123" in text
    # the shared institute and contact name of the fixtures: passes
    assert _plan(y, tmp_path / "e", named=True)["control_series_match"] == {
        "contact_institute": ["A Synthetic Institute"], "contact_name": ["A,,Person"]}


def test_blocker_2_the_match_is_derived_again_when_pins_are_verified(tmp_path):
    y, dest = _pinned_yaml(tmp_path, named=True)
    assert _plan(y, dest, named=True)["status"] == "ok"
    ds = cm.dataset(ACC, y)
    pp = dest / "meta" / ACC / "plan.json"
    pl = json.loads(pp.read_text())
    fe.verify_pins(ds, pl)
    with pytest.raises(cm.Refused, match="control series match differs"):
        fe.verify_pins(ds, {**pl, "control_series_match": {"contact_institute": ["typed by hand"]}})
    # a stored sister block edited so that the records no longer share an institute: refused too
    bad = json.loads(pp.read_text())
    bad["series"][SIS]["fields"]["contact_institute"] = ["Another Institute"]
    with pytest.raises(cm.Refused, match="do not share a contact_institute"):
        fe.verify_pins(ds, bad)


def test_a_superseries_of_the_own_series_listed_under_the_sister_is_refused(tmp_path):
    """The second reviewer's lesser point 3: a superseries holding the knockout samples, listed under the sister
    alone, satisfied the knockouts' membership check through the sister's declaration."""
    y = _yaml(tmp_path, named=True)
    y.write_text(y.read_text().replace("screen_series: [GSE910002]", "screen_series: [GSE910002, GSE910003]"))
    sup = {"GSM9101": OWN["GSM9101"], "GSM9201": SISTER["GSM9201"]}
    m = _sister_metas(sister_extra=["!Series_contact_institute = A Synthetic Institute\n"],
                      extra_series={"GSE910003": fe.parse_soft(_soft("GSE910003", sup, []))})
    pl = fe.plan(ACC, tmp_path / "d", y, metas=m)
    assert pl["status"] == "refused" and "['GSE910003'] hold samples of the dataset's own series" in pl["reason"]


SIS_ARCH = "GSE910002_WT_E8.5_3.tar.gz"
SIS_PATTERN = "      - series: GSE910002\n        patterns: ['^GSE910002_WT_E8\\.5_[0-9]+\\.tar\\.gz$']\n"


def test_blocker_3_a_one_library_sister_file_controls_several_knockouts(tmp_path, monkeypatch):
    """The second reviewer's reproduction: a sister archive holding one wild-type library (one flat triple, one
    sample), pinned split none as control for two genes, was refused ('ONE arm, gene and stage'). Now it loads,
    fetches and processes, and every unit takes d = 0.5."""
    y = _yaml(tmp_path, named=True)
    t = y.read_text()
    i = t.index("    control_series:")
    y.write_text(t[:i] + SIS_PATTERN + t[i:])
    assert [e["series"] for e in cm.dataset(ACC, y)["series_files"]] == [ACC, SIS]
    sister = {**SISTER, "GSM9203": _sample("GSM9203", "WT E8.5 3", "wild type", SIS)}
    metas = {ACC: _metas()[ACC], SIS: fe.parse_soft(_soft(SIS, sister, [SIS_ARCH]))}
    dest = tmp_path / "dest"
    assert fe.plan(ACC, dest, y, metas=metas)["status"] == "ok"
    genes = ("Aaa", "Bbb", "Ccc")
    block = "".join([
        fe.pin_template_from_series_file(ACC, dest, "ko", "Aaa", ["GSM9101", "GSM9102"], datasets=y),
        fe.pin_template_from_series_file(ACC, dest, "ko", "Bbb", ["GSM9103", "GSM9104"], datasets=y),
        fe.pin_template_from_series_file(ACC, dest, "ko", "Ccc", ["GSM9106", "GSM9107"], datasets=y),
        *[fe.pin_template_from_series_file(ACC, dest, "control", g, ["GSM9203"], datasets=y) for g in genes]])
    units = "".join([
        fe.series_file_template(ACC, dest, ACC, [ARCHIVES[0]], ["GSM9101", "GSM9102"], "folder",
                                ["AAA_1a=GSM9101", "AAA_1b=GSM9102"], datasets=y),
        fe.series_file_template(ACC, dest, ACC, LOOSE, ["GSM9103", "GSM9104"], "barcode_suffix",
                                ["-1=GSM9103", "-2=GSM9104"], datasets=y),
        fe.series_file_template(ACC, dest, ACC, [ARCHIVES[1]], ["GSM9106", "GSM9107"], "none", datasets=y),
        fe.series_file_template(ACC, dest, SIS, [SIS_ARCH], ["GSM9203"], "none", datasets=y)])
    _paste(y, "    pinned_series:\n" + fe.series_pin_template(ACC, dest, datasets=y) + "    pinned_series_files:\n"
           + units + "    pinned:\n" + block)
    ds = cm.dataset(ACC, y)   # loads: one arm (control), one stage, the same three genes
    assert sorted(p["gene"] for p in ds["pinned"] if p["accession"] == "GSM9203") == sorted(genes)
    assert fe.plan(ACC, dest, y, metas=metas)["status"] == "ok"
    rng = np.random.default_rng(5)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    staging = tmp_path / "staging"
    _stage(staging, rng, mu)
    _tar(staging / SIS_ARCH, dict(zip(("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"),
                                      _triple(rng, mu, np.ones(len(PANEL))), strict=True)))
    calls: list = []
    monkeypatch.setattr(fe, "_get", _fake_get(staging, calls))
    fe.fetch(ACC, dest, y)
    assert any(u.endswith(f"/series/GSE910nnn/{SIS}/suppl/{SIS_ARCH}") for u in calls)
    wt = tmp_path / "wt.h5ad"
    wt.write_bytes(b"stand-in for the WT reference")
    refp = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats({"genes": PANEL, "groups": GROUPS, "shares": {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5},
                       "mu": mu, "count": {"CM": 20, "Endoderm": 30, "Neural": 50}}, refp, cm.sha256_file(wt))
    qc = pr.run(ACC, dest, wt, refp, tmp_path / "processed", y, classifier=_clf)
    assert qc["status"] == "ok" and qc["sister_series_controls"] == ["GSM9203"]
    for g in genes:
        u = qc["units"][g]
        assert u["status"] == "ok" and u["d"] == 0.5 and u["control_source"].endswith("['GSM9203']"), (g, u)
    # the gene condition still guards a file of SEVERAL samples: two samples pinned for different genes refuse
    t = y.read_text()
    i = t.index("- accession: GSM9107")
    y.write_text(t[:i] + t[i:].replace('gene: "Ccc"', 'gene: "Aaa"', 1))
    with pytest.raises(cm.Refused, match="each pinned for the same genes"):
        cm.load_datasets(y)


def test_the_per_sample_archive_route_checks_every_member_like_a_series_archive(tmp_path):
    """Both reviewers' lesser point: read_sample's tar route read members with no checks (and _pick took the first
    match). Now every member goes through archive_members and the archive holds one triple in one folder."""
    rng = np.random.default_rng(13)
    lib = _lib(rng)
    ok = _tar(tmp_path / "GSM1_ok.tar.gz", {f"GSM1/{k}": v for k, v in lib.items()})
    M, _, _, bcs = pr.read_sample([ok])
    assert M.shape[0] == len(bcs) == 6
    s_lib = {f"S/{k}": v for k, v in lib.items()}
    cases = {"extra": ({**s_lib, "S/README.txt": b"x"}, "no rule"),
             "two_folders": ({**s_lib, **{f"T/{k}": v for k, v in lib.items()}}, "more than one folder"),
             "two_triples": ({**s_lib, "S/x_matrix.mtx.gz": b"y"}, "two 'matrix'"),
             "traversal": ({**s_lib, "S/../../evil/matrix.mtx.gz": b"y"}, "path traversal"),
             "screened": ({f"S/Gata4_{k}": v for k, v in lib.items()}, "held-out screen")}
    for name, (members, match) in cases.items():
        p = _tar(tmp_path / f"GSM2_{name}.tar.gz", members)
        with pytest.raises(cm.Refused, match=match):
            pr.read_sample([p])
    sym = tmp_path / "GSM3_sym.tar.gz"
    with tarfile.open(sym, "w:gz") as t:
        for k, v in s_lib.items():
            ti = tarfile.TarInfo(k)
            ti.size = len(v)
            t.addfile(ti, io.BytesIO(v))
        ti = tarfile.TarInfo("S/link")
        ti.type, ti.linkname = tarfile.SYMTYPE, "/etc/passwd"
        t.addfile(ti)
    with pytest.raises(cm.Refused, match="symlink"):
        pr.read_sample([sym])


def test_a_gzip_member_is_decompressed_in_chunks_and_refused_past_the_cap(tmp_path, monkeypatch):
    rng = np.random.default_rng(14)
    for k, v in _lib(rng).items():
        (tmp_path / f"GSM1_{k}").write_bytes(v)
    paths = sorted(tmp_path.iterdir())
    assert pr.read_sample(paths)[0].shape[0] == 6
    monkeypatch.setattr(pr, "MEMBER_MAX_BYTES", 64)
    with pytest.raises(cm.Refused, match="decompresses to more than 64 bytes"):
        pr.read_sample(paths)
