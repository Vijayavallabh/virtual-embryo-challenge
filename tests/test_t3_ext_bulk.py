"""Stage 9, phase 1: bulk whole-embryo units from GEO (experiments/t3_ext/bulk.py), on synthetic data only: the dataset
keys and pin forms, the per-sample and series-table intakes (the header step, the exact tie rule, pins by column,
usecols), the counts-only rule, the id routes, the math and every guard, the bulk table family and its sign rule, and
that nothing of the single-cell path changes. No network: fetch._open and fetch._get are stubs."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
T3 = ROOT / "experiments" / "t3_decomp"
for p in (EXT, T3):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import bulk  # noqa: E402
import combine as co  # noqa: E402
import fetch as fe  # noqa: E402
import make_panel_ids as mpi  # noqa: E402
import process as pr  # noqa: E402
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402


@pytest.fixture(autouse=True)
def _geo_bulk_forms_open_for_these_tests(monkeypatch):
    """The GEO bulk forms are closed in the committed code (t3x_common.GEO_BULK_OPEN is False); the synthetic tests
    of those forms open them for their own run only."""
    monkeypatch.setattr(cm, "GEO_BULK_OPEN", True)

KO = "Kozz"
PANEL = [f"g{i}" for i in range(37)] + ["Xist", "Ddx3y", KO]
GROUPS = ["CM", "Endoderm", "Neural"]
SHARES = {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5}
FILLERS = [f"f{i}" for i in range(100)]
ID = {s: f"ENSMUSG{100000 + i:011d}" for i, s in enumerate(PANEL + FILLERS + list(cm.SEX_GENES))}

YAML = """version: 1
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
    modality: bulk
    count_form: per_sample
    id_type: symbol
    gene_columns: [gene]
    sex_field: "characteristics_ch1:sex"
    somite_field: "characteristics_ch1:somites"
    screen_series: [GSE910001]
    stages_used: [E8.5]
    sample_file_patterns: ['_counts\\.']
    pinned: []  # GSE910001
  - accession: GSE910002
    status: cleared
    tier: unrelated
    study: GSE910002
    design: germline
    control_relation: same_series
    repository: geo
    modality: bulk
    count_form: series_table
    id_type: ensembl
    gene_columns: [gene_id]
    stage_from_series: overall_design
    screen_series: [GSE910002]
    stages_used: [E8.5]
    pinned: []  # GSE910002
"""
MARK = {"GSE910001": "    pinned: []  # GSE910001\n", "GSE910002": "    pinned: []  # GSE910002\n"}


def stub(gsm: str) -> str:
    return "GSM" + gsm[3:-3] + "nnn"


def sample(gsm, geno, stage="E8.5", sex="female", som=12, files=("counts.txt.gz",), extra=None):
    ch = [f"genotype: {geno}"] + ([f"developmental stage: {stage}"] if stage else []) + [f"sex: {sex}",
                                                                                         f"somites: {som}"]
    return {"title": gsm, "source_name_ch1": "whole embryo", "characteristics_ch1": ch + (extra or []),
            **{f"supplementary_file_{i + 1}": f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{stub(gsm)}/{gsm}/suppl/"
               f"{gsm}_{f}" for i, f in enumerate(files)}}


def soft(acc, samples, design="whole embryos", series_files=(), relation=()):
    st = "GSE" + acc[3:-3] + "nnn"
    lines = [f"^SERIES = {acc}", f"!Series_title = A bulk series {acc}", f"!Series_overall_design = {design}",
             "!Series_summary = whole embryos, bulk RNA-seq", "!Series_contact_name = A,,Person",
             *[f"!Series_relation = {r}" for r in relation],
             *[f"!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/{st}/{acc}/suppl/{acc}_{f}"
               for f in series_files]]
    for gsm, f in samples.items():
        lines.append(f"^SAMPLE = {gsm}")
        for k, vals in {**f, "series_id": acc}.items():
            for v in (vals if isinstance(vals, list) else [vals]):
                lines.append(f"!Sample_{k} = {v}")
    return "\n".join(lines) + "\n"


@pytest.fixture()
def y(tmp_path):
    p = tmp_path / "datasets.yaml"
    p.write_text(YAML)
    return p


def make_ref(tmp_path, seed=5):
    rng = np.random.default_rng(seed)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    wt = tmp_path / "wt.h5ad"
    wt.write_bytes(b"stand-in WT")
    refp = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats({"genes": PANEL, "groups": GROUPS, "shares": SHARES, "mu": mu,
                       "count": {"CM": 20, "Endoderm": 30, "Neural": 50}}, refp, cm.sha256_file(wt))
    return wt, refp, tb.load_ref_stats(refp)


def lib(rng, ref, fold=None, sex="female", own=4000, panel_n=300_000, filler=17_000, profile=None):
    """{symbol: count} of one whole-embryo library: the panel drawn from the reference's pseudo-bulk (times fold),
    fillers, the sex genes, the knocked-out gene."""
    p = np.asarray(profile if profile is not None else tb.whole_den(ref), float).copy()
    if fold is not None:
        p = p * fold
    xp = rng.multinomial(panel_n, p / p.sum()).astype(float)
    out = dict(zip(PANEL, xp, strict=True))
    out.update({f: float(rng.poisson(filler)) for f in FILLERS})
    out["Xist"] = 3000.0 if sex == "female" else 0.0
    for g in cm.SEX_Y_GENES:
        out[g] = 0.0 if sex == "female" else 150.0
    out[KO] = float(own)
    out["__no_feature"] = 5000.0
    return out


def per_sample_file(counts: dict, col="counts") -> bytes:
    return gzip.compress(("gene\t" + col + "\n" + "".join(f"{g}\t{int(v)}\n" for g, v in counts.items())).encode())


def series_table(cols: dict[str, dict], text_cols=(), version=".3") -> bytes:
    genes = list(next(iter(cols.values())))
    head = ["gene_id", *cols, *text_cols]
    rows = ["\t".join(head)]
    for g in genes:
        gid = ID.get(g, g) + (version if g in ID else "")
        rows.append("\t".join([gid, *[str(int(c[g])) for c in cols.values()], *["some text" for _ in text_cols]]))
    return gzip.compress(("\n".join(rows) + "\n").encode())


def fake_get(staging: Path, calls: list):
    def get(url, dst, max_bytes):
        calls.append(url)
        dst = cm.outside_repo(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(staging / Path(url).name, dst)
        return dst.stat().st_size, cm.sha256_file(dst)
    return get


def paste(y: Path, acc: str, pins: str, series: str, table: str = "") -> None:
    block = (("    pinned_series_files:\n" + table if table else "") + "    pinned_series:\n" + series + "    pinned:\n"
             + pins)
    y.write_text(y.read_text().replace(MARK[acc], block))


def panel_file(tmp_path, ref) -> tuple[Path, str]:
    rows = [(ID[s], s) for s in PANEL + FILLERS + list(cm.SEX_GENES)]
    text = mpi.build(rows, ref["genes"], [KO], "GSE0/GSM0/GSM0_features.tsv.gz", "0" * 64)
    p = tmp_path / "panel_ids.tsv"
    p.write_text(text)
    return p, cm.sha256_file(p)


# ---------------------------------------------------------------------------------------------------------------------
# a synthetic per_sample dataset, walked from plan to tables


PS = {f"GSM91000{i}": ("Kozz -/-" if i <= 3 else "wild type") for i in range(1, 7)}


def run_per_sample(tmp_path, y, monkeypatch, folds=None, sexes=None, owns=None, extra_samples=None, cols=None):
    wt, refp, ref = make_ref(tmp_path)
    rng = np.random.default_rng(1)
    dest = tmp_path / "dest"
    smp = {g: sample(g, geno, sex=(sexes or {}).get(g, "female"), som=11 + (i % 3)) for i, (g, geno) in
           enumerate(PS.items())}
    smp.update(extra_samples or {})
    metas = {"GSE910001": fe.parse_soft(soft("GSE910001", smp))}
    pl = fe.plan("GSE910001", dest, y, metas=metas)
    assert pl["status"] == "ok", pl.get("reason")
    ko, ctl = list(PS)[:3], list(PS)[3:]
    pins = bulk.pin_template("GSE910001", dest, "ko", KO, ko, y) + bulk.pin_template("GSE910001", dest, "control",
                                                                                      KO, ctl, y)
    paste(y, "GSE910001", pins, fe.series_pin_template("GSE910001", dest, y))
    assert fe.plan("GSE910001", dest, y, metas=metas)["status"] == "ok"
    fold = np.ones(len(PANEL))
    fold[:4] = 2.0
    staging = tmp_path / "staging"
    staging.mkdir()
    for g in PS:
        c = lib(rng, ref, fold=(folds or {}).get(g, fold if g in ko else None), sex=(sexes or {}).get(g, "female"),
                own=(owns or {}).get(g, 200 if g in ko else 4000))
        (staging / f"{g}_counts.txt.gz").write_bytes(per_sample_file(c, col=(cols or {}).get(g, "counts")))
    calls: list = []
    monkeypatch.setattr(fe, "_get", fake_get(staging, calls))
    fe.fetch("GSE910001", dest, y)
    return dest, wt, refp, ref, calls


def test_a_per_sample_bulk_dataset_from_plan_to_tables(tmp_path, y, monkeypatch, capsys):
    dest, wt, refp, ref, calls = run_per_sample(tmp_path, y, monkeypatch)
    assert sorted(Path(u).name for u in calls) == sorted(f"{g}_counts.txt.gz" for g in PS)
    out = tmp_path / "processed"
    qc = pr.run("GSE910001", dest, wt, refp, out, y)
    assert qc["modality"] == "bulk" and qc["status"] == "ok", qc.get("reason")
    u = qc["units"][KO]
    assert u["status"] == "ok", u.get("reason")
    st = u["stages"]["E8.5"]
    assert st["m"] == 3 and st["dropped"] == {} and st["genotype_in_counts"].startswith("the knockout arm is lower")
    assert u["k_by_stage"]["E8.5"] == pytest.approx((3 / 4) ** 2 * 1.0)
    un = json.loads(Path(u["unit_file"]).read_text())
    assert un["modality"] == "bulk" and set(un["L"]) == {"__whole__"}
    L = np.array([np.nan if v is None else v for v in un["L"]["__whole__"]])
    for g in ("Xist", "Ddx3y", KO):   # always not measured, and the knocked-out gene's own row
        assert np.isnan(L[PANEL.index(g)])
    assert np.nanmin(L[:4]) > 0 and np.nanmax(np.abs(L[4:37])) < np.nanmin(L[:4])
    tdir = Path(u["tables"][f"GSE910001__{KO}"]["file"]).parents[1]
    assert not list(tdir.glob("*_per_group.csv"))
    r = co.combine(out, refp, tmp_path / "tables", y, dest=dest)
    assert r["tables"] == ["bulk_all", "bulk_study_GSE910001"]
    man = json.loads((tmp_path / "tables" / "bulk_all" / "bulk_all_manifest.json").read_text())
    assert man["per_group"] == "none (bulk)" and man["sign_rule"] == "single study" and man["modality"] == "bulk"
    assert man["checked"] == "by the practice read named in its TRIES row only"
    assert not (tmp_path / "tables" / "bulk_all" / "bulk_all_per_group.csv").exists()
    printed = capsys.readouterr().out
    assert not re.search(r"\d\.\d{3,}", printed)   # no gene-level value is printed


def test_the_one_unpinned_count_file_is_never_opened_and_g8_refuses_a_copied_file(tmp_path, y, monkeypatch):
    dest, wt, refp, ref, _ = run_per_sample(tmp_path, y, monkeypatch, extra_samples={
        "GSM910009": sample("GSM910009", "Kozz +/-")})
    assert not (dest / "raw" / "GSE910001" / "GSM910009").exists()
    # two pinned libraries with one file's bytes: G8
    a = dest / "raw" / "GSE910001" / "GSM910001" / "GSM910001_counts.txt.gz"
    b = dest / "raw" / "GSE910001" / "GSM910002" / "GSM910002_counts.txt.gz"
    shutil.copyfile(a, b)
    mp_ = dest / "raw" / "GSE910001" / "MANIFEST.json"
    man = json.loads(mp_.read_text())
    man["files"]["GSM910002/GSM910002_counts.txt.gz"]["sha256"] = cm.sha256_file(b)
    mp_.write_text(json.dumps(man))
    with pytest.raises(cm.Refused, match="G8"):
        pr.run("GSE910001", dest, wt, refp, tmp_path / "processed", y)


def test_g9_refuses_a_knockout_arm_that_is_not_lower_on_its_own_gene_and_prints_no_value(tmp_path, y, monkeypatch,
                                                                                            capsys):
    dest, wt, refp, ref, _ = run_per_sample(tmp_path, y, monkeypatch,
                                            owns={g: (6000 if g in list(PS)[:3] else 4000) for g in PS})
    qc = pr.run("GSE910001", dest, wt, refp, tmp_path / "processed", y)
    u = qc["units"][KO]
    assert u["status"] == "refused_genotype" and "G9" in u["reason"] and "unit_file" not in u
    assert not re.search(r"\d\.\d{3,}", json.dumps(u["stages"]["E8.5"]["genotype_in_counts"]))
    assert qc["status"] == "left_out"


def test_g3_drops_a_library_whose_clear_call_contradicts_its_recorded_sex(tmp_path, y, monkeypatch):
    # GSM910004 is recorded female but its counts call male: dropped; three controls become two (still m >= 2)
    dest, wt, refp, ref, _ = run_per_sample(tmp_path, y, monkeypatch)
    raw = dest / "raw" / "GSE910001" / "GSM910004" / "GSM910004_counts.txt.gz"
    t = gzip.decompress(raw.read_bytes()).decode()
    t = re.sub(r"\nXist\t\d+", "\nXist\t0", t)
    for g in cm.SEX_Y_GENES:
        t = re.sub(rf"\n{g}\t\d+", f"\n{g}\t300", t)
    raw.write_bytes(gzip.compress(t.encode()))
    mp_ = dest / "raw" / "GSE910001" / "MANIFEST.json"
    man = json.loads(mp_.read_text())
    man["files"]["GSM910004/GSM910004_counts.txt.gz"]["sha256"] = cm.sha256_file(raw)
    mp_.write_text(json.dumps(man))
    qc = pr.run("GSE910001", dest, wt, refp, tmp_path / "processed", y)
    st = qc["units"][KO]["stages"]["E8.5"]
    assert st["dropped"] == {"GSM910004": "G3 sex call contradicts the recorded sex"} and st["m"] == 2


# ---------------------------------------------------------------------------------------------------------------------
# a synthetic series_table dataset (Ensembl route): the header step, ties, pins by column, usecols


ST = {"GSM920001": ("KO_rep1", "Kozz null"), "GSM920002": ("KO_rep2", "Kozz null"), "GSM920003": ("WT_rep1", "wt"),
      "GSM920004": ("WT_rep2", "wt"), "GSM920005": ("HET_rep1", "Kozz het")}
TSHORT = "counts_table.txt.gz"
TNAME = f"GSE910002_{TSHORT}"   # the file name the series' own folder offers


def st_samples():
    return {g: {"title": t, "source_name_ch1": "whole embryo", "characteristics_ch1": [f"genotype: {geno}"]}
            for g, (t, geno) in ST.items()}


def st_metas(design="whole embryos at E8.5"):
    return {"GSE910002": fe.parse_soft(soft("GSE910002", st_samples(), design=design, series_files=(TSHORT,)))}


def opener(data: bytes, calls: list):
    class R(io.BytesIO):
        n_read = 0

        def read(self, n=-1):
            b = super().read(n)
            R.n_read += len(b)
            return b

    def op(url):
        calls.append(url)
        return R(data)
    return op, R


def run_series_table(tmp_path, y, monkeypatch, het_text=True):
    wt, refp, ref = make_ref(tmp_path)
    rng = np.random.default_rng(2)
    dest = tmp_path / "dest"
    assert fe.plan("GSE910002", dest, y, metas=st_metas())["status"] == "ok"
    fold = np.ones(len(PANEL))
    fold[:4] = 2.0
    cols = {t: lib(rng, ref, fold=fold if t.startswith("KO") else None, own=200 if t.startswith("KO") else 4000)
            for t, _ in ST.values() if not t.startswith("HET")}
    data = series_table(cols, text_cols=("HET_rep1",) if het_text else ())
    calls: list = []
    op, R = opener(data, calls)
    monkeypatch.setattr(fe, "_open", op)
    rec = bulk.header("GSE910002", dest, "GSE910002", TNAME, y)
    assert rec["status"] == "ok", rec.get("reason")
    assert calls == ["https://ftp.ncbi.nlm.nih.gov/geo/series/GSE910nnn/GSE910002/suppl/" + TNAME]
    pins = (bulk.pin_template("GSE910002", dest, "ko", KO, ["GSM920001", "GSM920002"], y)
            + bulk.pin_template("GSE910002", dest, "control", KO, ["GSM920003", "GSM920004"], y))
    paste(y, "GSE910002", pins, fe.series_pin_template("GSE910002", dest, y), bulk.table_template("GSE910002", dest, y))
    pl = fe.plan("GSE910002", dest, y, metas=st_metas())
    assert pl["status"] == "ok", pl.get("reason")
    staging = tmp_path / "staging2"
    staging.mkdir()
    (staging / TNAME).write_bytes(data)
    got: list = []
    monkeypatch.setattr(fe, "_get", fake_get(staging, got))
    return dest, wt, refp, ref, data, R, got


def test_the_header_step_reads_only_the_first_line_and_writes_only_it(tmp_path, y, monkeypatch):
    wt, refp, ref = make_ref(tmp_path)
    dest = tmp_path / "dest"
    fe.plan("GSE910002", dest, y, metas=st_metas())
    rng = np.random.default_rng(0)
    big = series_table({t: lib(rng, ref) for t in ("KO_rep1", "KO_rep2", "WT_rep1", "WT_rep2", "HET_rep1")})
    # a long table that does not compress well: the step must stop reading soon after the first newline
    pad = "".join(f"{rng.integers(1 << 62):x}\t{rng.integers(1 << 62):x}\n" for _ in range(60_000)).encode()
    big = gzip.compress(gzip.decompress(big) + pad)
    assert len(big) > 1 << 20
    op, R = opener(big, [])
    monkeypatch.setattr(fe, "_open", op)
    rec = bulk.header("GSE910002", dest, "GSE910002", TNAME, y)
    assert rec["status"] == "ok" and R.n_read < len(big)
    txt = (dest / "meta" / "GSE910002" / f"header_{TNAME}.txt").read_text()
    assert txt == "gene_id\tKO_rep1\tKO_rep2\tWT_rep1\tWT_rep2\tHET_rep1\n"
    assert [t["how"] for t in rec["ties"]] == ["gene column", *["title"] * 5]
    assert rec["header_sha256"] == hashlib.sha256(txt.rstrip("\n").encode()).hexdigest()
    # nothing else of the table is on disk
    assert sorted(p.name for p in (dest / "meta" / "GSE910002").iterdir() if "header" in p.name) == [
        f"header_{TNAME}.json", f"header_{TNAME}.txt"]


def test_a_first_line_over_the_cap_or_without_newline_is_refused():
    long = b"a" * (cm.HEADER_MAX_BYTES + 10) + b"\nrest"
    with pytest.raises(cm.Refused, match="longer than"):
        bulk.first_line(io.BytesIO(long).read, False, "t")
    with pytest.raises(cm.Refused, match="longer than"):
        bulk.first_line(io.BytesIO(gzip.compress(long)).read, True, "t")
    with pytest.raises(cm.Refused, match="before its first newline"):
        bulk.first_line(io.BytesIO(b"gene\tx").read, False, "t")
    assert bulk.first_line(io.BytesIO(b"a\tb\r\nc").read, False, "t") == b"a\tb"


@pytest.mark.parametrize("hdr,status,withheld", [("gene_id\tCtnnb1_KO1\tWT_rep1", "refused", True),
                                                 ("gene_id\tAxin2_KO\tWT_rep1", "refused", False)])
def test_a_header_with_a_held_out_term_is_withheld_and_a_review_term_refuses(tmp_path, y, monkeypatch, hdr, status,
                                                                             withheld):
    dest = tmp_path / "dest"
    fe.plan("GSE910002", dest, y, metas=st_metas())
    monkeypatch.setattr(fe, "_open", opener(gzip.compress((hdr + "\n1\t2\t3\n").encode()), [])[0])
    rec = bulk.header("GSE910002", dest, "GSE910002", TNAME, y)
    assert rec["status"] == status and rec["withheld"] is withheld
    stored = (dest / "meta" / "GSE910002" / f"header_{TNAME}.json").read_text()
    if withheld:
        assert "Ctnnb1" not in stored and not (dest / "meta" / "GSE910002" / f"header_{TNAME}.txt").exists()
        assert rec["n_columns"] == 3
    else:
        assert "review" in rec["reason"]


def test_a_review_term_in_a_header_is_cleared_only_by_an_entry_naming_the_series(tmp_path, y, monkeypatch):
    y.write_text(YAML.replace("    stage_from_series: overall_design\n", "    stage_from_series: overall_design\n"
                              "    review_cleared: [{sample: GSE910002, term: axin2, reason: a header read}]\n"))
    dest = tmp_path / "dest"
    fe.plan("GSE910002", dest, y, metas=st_metas())
    monkeypatch.setattr(fe, "_open", opener(gzip.compress(b"gene_id\tKO_rep1\tAxin2\n1\t2\t3\n"), [])[0])
    rec = bulk.header("GSE910002", dest, "GSE910002", TNAME, y)
    assert "review" not in (rec.get("reason") or "")   # cleared; it fails on the tie rule instead (Axin2 ties nothing)
    assert rec["status"] == "refused" and "tie rule" in rec["reason"]


def _rows(samples):
    return [{"accession": g, "fields": f} for g, f in samples.items()]


def test_the_exact_tie_rule():
    s = {"GSM1": {"title": ["KO 1"], "description": ["lib A"], "characteristics_ch1": ["id: S1", "genotype: ko"]},
         "GSM2": {"title": ["WT 1"], "characteristics_ch1": ["id: S2", "genotype: wt"]},
         "GSM3": {"title": ["WT 2"], "characteristics_ch1": ["id: S3", "genotype: wt"]}}
    ties, problems = bulk.tie_table(["gene", "GSM1", "WT 1", "S3"], _rows(s), ["gene"])
    assert problems == [] and [t["sample"] for t in ties] == [None, "GSM1", "GSM2", "GSM3"]
    assert [t["how"] for t in ties[1:]] == ["accession", "title", "characteristics value"]
    assert bulk.tie_table(["gene", "lib A"], _rows(s), ["gene"])[1] == []      # a whole description value
    # a token-only match does not tie
    assert bulk.tie_table(["gene", "KO"], _rows(s), ["gene"])[1]
    assert bulk.tie_table(["gene", "1"], _rows(s), ["gene"])[1]
    # a header tying to 2 samples; two headers tying to one sample; a missing or doubled gene column; a doubled header
    assert any("ties to 2" in p for p in bulk.tie_table(["gene", "wt"], _rows(s), ["gene"])[1])
    assert any("ties to 2 headers" in p for p in bulk.tie_table(["gene", "GSM1", "KO 1"], _rows(s), ["gene"])[1])
    assert bulk.tie_table(["GSM1", "GSM2"], _rows(s), ["gene"])[1]
    assert bulk.tie_table(["gene", "gene", "GSM1"], _rows(s), ["gene"])[1]
    assert bulk.tie_table(["gene", "GSM1", "GSM1"], _rows(s), ["gene"])[1]
    assert bulk.split_header(b'"gene","GSM1"', "x.csv") == ["gene", "GSM1"]
    with pytest.raises(cm.Refused, match="quote"):
        bulk.split_header(b'gene\ta"b', "x.tsv")


def test_a_series_table_dataset_from_header_to_tables_reads_only_pinned_columns(tmp_path, y, monkeypatch):
    dest, wt, refp, ref, data, R, got = run_series_table(tmp_path, y, monkeypatch)
    # the Ensembl route is refused while the id map is not pinned (the default until 2026-11-20; the committed pin is
    # tested in test_t3_ext_panel_ids.py)
    monkeypatch.setattr(cm, "PANEL_IDS_SHA256", None)
    with pytest.raises(cm.Refused, match="PANEL_IDS_SHA256 is None"):
        fe.fetch("GSE910002", dest, y)
    assert "PANEL_IDS_SHA256" in fe.check_gate(["GSE910002"], y)[0]
    pfile, psha = panel_file(tmp_path, ref)
    monkeypatch.setattr(cm, "PANEL_IDS", pfile)
    monkeypatch.setattr(cm, "PANEL_IDS_SHA256", psha)
    assert fe.check_gate(["GSE910002"], y) == []
    man = fe.fetch("GSE910002", dest, y)
    assert [Path(u).name for u in got] == [TNAME]
    assert man["files"][f"_series/{TNAME}"]["header_sha256"] == cm.dataset("GSE910002", y)["pinned_series_files"][0][
        "header_sha256"]
    # the HET column holds text: it would not parse as counts, and it is never parsed (usecols)
    qc = pr.run("GSE910002", dest, wt, refp, tmp_path / "processed", y)
    assert qc["status"] == "ok", qc.get("reason")
    assert qc["read"]["columns_read"] == 4 and qc["units"][KO]["status"] == "ok"
    assert qc["panel_match"]["GSM920001"]["route"] == "ensembl"
    # a changed first line refuses on header_sha256 (process reads the table again)
    p = dest / "raw" / "GSE910002" / "_series" / TNAME
    with pytest.raises(cm.Refused, match="header_sha256"):
        bulk.read_series_table(p, TNAME, ["gene_id"], ["KO_rep1"], "0" * 64)


def test_a_pinned_column_must_be_the_tie_and_the_stored_ties_must_derive_again(tmp_path, y, monkeypatch):
    dest, *_ = run_series_table(tmp_path, y, monkeypatch)
    ds = cm.dataset("GSE910002", y)
    pl = json.loads((dest / "meta" / "GSE910002" / "plan.json").read_text())
    bad = {**ds, "pinned": [{**ds["pinned"][0], "column": "KO_rep2"}, *ds["pinned"][1:]]}
    with pytest.raises(cm.Refused):
        fe.verify_pins(bad, pl)
    with pytest.raises(cm.Refused, match="is not the header the tie table ties"):
        bulk.pin_template("GSE910002", dest, "ko", KO, ["GSM920001"], y, column="KO_rep2")
    pl2 = json.loads(json.dumps(pl))
    pl2["count_table_headers"][TNAME]["ties"][1]["sample"] = "GSM920002"
    with pytest.raises(cm.Refused, match="derived again"):
        fe.verify_pins(ds, pl2)
    pl3 = json.loads(json.dumps(pl))
    pl3["count_table_headers"][TNAME]["text"] += "\tX"
    with pytest.raises(cm.Refused, match="hash"):
        fe.verify_pins(ds, pl3)


def test_fetch_refuses_a_table_whose_first_line_is_not_the_pinned_one(tmp_path, y, monkeypatch):
    dest, wt, refp, ref, data, R, got = run_series_table(tmp_path, y, monkeypatch)
    pfile, psha = panel_file(tmp_path, ref)
    monkeypatch.setattr(cm, "PANEL_IDS", pfile)
    monkeypatch.setattr(cm, "PANEL_IDS_SHA256", psha)
    staging = tmp_path / "staging3"
    staging.mkdir()
    (staging / TNAME).write_bytes(gzip.compress(b"gene_id\tKO_rep1\tOTHER\n1\t2\t3\n"))
    monkeypatch.setattr(fe, "_get", fake_get(staging, []))
    with pytest.raises(cm.Refused, match="first line does not hash"):
        fe.fetch("GSE910002", dest, y)
    assert not (dest / "raw" / "GSE910002" / "_series" / TNAME).exists()


def test_stage_from_series():
    ds = {"stage_from_series": "overall_design"}
    assert bulk.series_stage({"overall_design": ["whole embryos at E8.5"]}, "overall_design")[0] == 8.5
    assert bulk.series_stage({"overall_design": ["E8.5 and E9.5 embryos"]}, "overall_design")[0] is None
    assert bulk.series_stage({"overall_design": ["E8.5-E9.5 embryos"]}, "overall_design")[0] is None
    assert bulk.series_stage({"overall_design": ["E10.5 embryos"]}, "overall_design")[0] is None
    sf = {"overall_design": ["whole embryos at E8.5"]}
    assert bulk.row_stage(ds, {"title": ["KO 1"]}, sf)[0] == 8.5
    assert bulk.row_stage(ds, {"title": ["KO 1"], "characteristics_ch1": ["stage: E8.5"]}, sf)[0] == 8.5
    e, how = bulk.row_stage(ds, {"title": ["KO 1"], "characteristics_ch1": ["stage: E9.5"]}, sf)
    assert e is None and "E9.5" in how
    assert bulk.row_stage(ds, {"title": ["KO 1"], "characteristics_ch1": ["stage: E8.5-E9"]}, sf)[0] is None


def test_somite_only_staging_cannot_be_pinned(tmp_path, y):
    smp = {g: sample(g, geno) for g, geno in PS.items()}
    smp["GSM910001"]["characteristics_ch1"] = ["genotype: Kozz -/-", "developmental stage: 12 somites",
                                               "sex: female", "somites: 12"]
    metas = {"GSE910001": fe.parse_soft(soft("GSE910001", smp))}
    fe.plan("GSE910001", tmp_path / "d", y, metas=metas)
    with pytest.raises(cm.Refused, match="no single stage"):
        bulk.pin_template("GSE910001", tmp_path / "d", "ko", KO, ["GSM910001"], y)


# ---------------------------------------------------------------------------------------------------------------------
# per_sample file rules


def _rec(**kw):
    """A sample record as parse_soft gives it (every value a list)."""
    return {k: v if isinstance(v, list) else [v] for k, v in sample(**kw).items()}


def test_per_sample_files_own_folder_one_text_table_xlsx_rds_and_reads_refused(y):
    ds = cm.dataset("GSE910001", y)
    rec = _rec(gsm="GSM910001", geno="x", files=("counts.txt.gz",))
    assert [f["name"] for f in bulk.sample_files(ds, rec, "GSM910001")[0]] == ["GSM910001_counts.txt.gz"]
    for f, why in (("counts.xlsx", "not a plain text"), ("counts.rds", "not a plain text"),
                   ("R1.fastq.gz", "denied")):
        ok, other = bulk.sample_files(ds, _rec(gsm="GSM910001", geno="x", files=(f,)), "GSM910001")
        assert ok == [] and why in other[0]["reason"]
    rec = _rec(gsm="GSM910001", geno="x", files=("counts.txt.gz", "counts2_counts.tsv"))
    ok, other = bulk.sample_files(ds, rec, "GSM910001")
    assert ok == [] and all("one text table per sample" in o["reason"] for o in other)
    rec = _rec(gsm="GSM910001", geno="x", files=())
    rec["supplementary_file_1"] = ["ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM910nnn/GSM910002/suppl/"
                                   "GSM910002_counts.txt"]
    ok, other = bulk.sample_files(ds, rec, "GSM910001")
    assert ok == [] and "own GEO folder" in other[0]["reason"]
    ds2 = cm.dataset("GSE910002", y)
    ok, other = bulk.sample_files(ds2, _rec(gsm="GSM920001", geno="x"), "GSM920001")
    assert ok == [] and "count_form series_table" in other[0]["reason"]
    ok, other = bulk.series_files(ds2, "GSE910002", {"supplementary_file": [
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE910nnn/GSE910002/suppl/GSE910002_counts.xlsx",
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE910nnn/GSE910002/suppl/GSE910002_counts.tsv.gz"]})
    assert [f["name"] for f in ok] == ["GSE910002_counts.tsv.gz"] and "not a plain text" in other[0]["reason"]


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_bytes(gzip.compress(text.encode()) if name.endswith(".gz") else text.encode())
    return p


def test_read_count_file_rules(tmp_path):
    def yes(h):
        return ["GSM1"]

    def no(h):
        return []

    def other(h):
        return ["GSM2"]

    def screen(fields, text):
        return None

    def rd(p, name, count_column, hits):
        return bulk.read_count_file(p, name, ["gene"], count_column, hits, "GSM1", screen)

    p = _write(tmp_path, "a.txt", "gene\tcounts\ng1\t5\ng2\t7\n")
    ids, sym, x, col = rd(p, "a.txt", None, no)
    assert ids == ["g1", "g2"] and x.tolist() == [5, 7] and col == "counts"
    assert rd(p, "a.txt", None, yes)[2].tolist() == [5, 7]
    with pytest.raises(cm.Refused, match="another sample"):
        rd(p, "a.txt", None, other)    # round 2: a single count column tying to another sample
    p = _write(tmp_path, "b.txt", "gene\tcounts\tLength\ng1\t5\t900\n")
    with pytest.raises(cm.Refused, match="Length"):
        rd(p, "b.txt", None, yes)
    assert rd(p, "b.txt", "counts", yes)[2].tolist() == [5]
    with pytest.raises(cm.Refused, match="does not tie"):
        rd(p, "b.txt", "counts", no)    # count_column naming another sample's header
    p = _write(tmp_path, "c.txt", "g1\t5\ng2\t7\n")
    with pytest.raises(cm.Refused, match="headerless"):
        rd(p, "c.txt", None, yes)
    p = _write(tmp_path, "d_tpm.txt", "gene\tvalue\ng1\t5\n")
    with pytest.raises(cm.Refused, match="normalised"):
        rd(p, "d_tpm.txt", None, no)
    p = _write(tmp_path, "e.csv.gz", "gene,counts\ng1,5\n")
    assert rd(p, "e.csv.gz", None, no)[2].tolist() == [5]
    # the screen runs on the first line before any value is parsed (a value that would refuse is never reached)
    p = _write(tmp_path, "f.txt", "gene\tcounts\ng1\tnot a number\n")
    seen = []

    def stop(fields, text):
        seen.append(fields)
        raise cm.Refused("stopped by the screen")

    with pytest.raises(cm.Refused, match="stopped by the screen"):
        bulk.read_count_file(p, "f.txt", ["gene"], None, no, "GSM1", stop)
    assert seen == [["gene", "counts"]]
    # a tab-separated table is read without quote handling: a quote in a value cannot merge rows
    p = _write(tmp_path, "g.txt", 'gene\tcounts\n"g1\t5\ng2"\t7\n')
    ids, _, x, _ = rd(p, "g.txt", None, no)
    assert ids == ['"g1', 'g2"'] and x.tolist() == [5, 7]


def test_counts_only():
    with pytest.raises(cm.Refused, match="non-integral"):
        bulk.counts_column(["1.5", "2"], "x")
    with pytest.raises(cm.Refused, match="negative"):
        bulk.counts_column(["-1", "2"], "x")
    with pytest.raises(cm.Refused, match="not a number"):
        bulk.counts_column(["NA", "2"], "x")
    with pytest.raises(cm.Refused, match="non-finite"):
        bulk.counts_column(["inf", "2"], "x")
    assert bulk.counts_column(["3.0000000001", "2"], "x").tolist() == [3, 2]
    with pytest.raises(cm.Refused, match="CPM or TPM"):
        bulk.counts_only([1e6, 1e6 * 1.004, 999_100], "x")
    bulk.counts_only([1e6, 2.3e6], "x")
    assert bulk.library_total(["a", "__no_feature", "N_unmapped", "b"], np.array([1.0, 5, 7, 2])) == 3.0


# ---------------------------------------------------------------------------------------------------------------------
# the math


def test_mean_of_logs_scaling_median_floor_clip_and_masks():
    rng = np.random.default_rng(0)
    P = 30
    base = rng.uniform(50, 500, P)
    X_ctl = np.vstack([rng.poisson(base) for _ in range(3)]).astype(float)
    X_ko = X_ctl.copy()
    X_ko[2] *= 100                     # one deep library
    T_ctl = X_ctl.sum(1)
    T_ko = X_ko.sum(1)
    pres = np.ones(P, bool)
    L = bulk.bulk_log_ratio(X_ko, X_ctl, T_ko, T_ctl, pres)
    assert np.nanmax(np.abs(L)) < 0.5   # the deep library weighs one third, as any library
    # a per-library normalisation constant is a gene-independent offset: the finished phi does not change
    den = rng.dirichlet(np.ones(P))
    L2 = bulk.bulk_log_ratio(X_ko, X_ctl, T_ko * np.array([2.0, 3.0, 0.5]), T_ctl, pres)
    assert np.allclose(tb.finish(L, den), tb.finish(L2, den))
    # the median floor: a gene under 20 in the median of either arm is NaN
    X_ctl[:, 0] = [100, 5, 5]
    assert np.isnan(bulk.bulk_log_ratio(X_ko, X_ctl, T_ko, T_ctl, pres)[0])
    # the clip
    X_ko[:, 1] = 1e7
    assert bulk.bulk_log_ratio(X_ko, X_ctl, T_ko, T_ctl, pres)[1] == pytest.approx(np.log(4.0))
    # absent and masked genes
    pres[3] = False
    mask = np.zeros(P, bool)
    mask[4] = True
    L3 = bulk.bulk_log_ratio(X_ko, X_ctl, T_ko, T_ctl, pres, mask)
    assert np.isnan(L3[3]) and np.isnan(L3[4])


def test_shrink_k_bulk_equals_shrink_k_without_a_cell_term():
    for m in range(1, 21):
        for d in (0.5, 1.0):
            assert cm.shrink_k_bulk(m, d) == cm.shrink_k([1e9, 1e9], m, d)
    assert cm.shrink_k_bulk(0, 1.0) == 0.0


def test_sex_call_and_balance():
    assert bulk.sex_call(0, [10, 10, 0, 0], 1e6) == "male"
    assert bulk.sex_call(20, [0, 0, 0, 0], 1e6) == "female"
    assert bulk.sex_call(20, [20, 0, 0, 0], 1e6) == "unclear"
    assert bulk.sex_call(None, [1, 1, 1, 1], 1e6) == "not possible"
    assert bulk.sex_balance(["male", "female"], ["male", "female"])[0] == "balanced"
    assert bulk.sex_balance(["male", "male"], ["female", "female"])[0] == "imbalanced"
    assert bulk.sex_balance(["male", "not possible"], ["female"])[0] == "not possible"
    assert bulk.sex_balance(["unclear"], ["female"])[0] == "not possible"


# ---------------------------------------------------------------------------------------------------------------------
# the guards on a synthetic unit (bulk.bulk_unit)


def _unit_libs(ref, rng, n_ko=3, n_ctl=3, sexes=None, fold=None, profile_ctl=None, own_ko=100, totals=None,
               somites=None):
    libs, u = {}, {"gene": KO, "ko": [], "control": [], "stages": [8.5]}
    for arm, n in (("ko", n_ko), ("control", n_ctl)):
        for i in range(n):
            lab = f"{arm}{i}"
            sex = (sexes or {}).get(lab, "female")
            prof = profile_ctl if (arm == "control" and profile_ctl is not None) else None
            c = lib(rng, ref, fold=fold if arm == "ko" else None, sex=sex, own=own_ko if arm == "ko" else 4000,
                    profile=prof)
            ids = list(c)
            x = np.array([c[g] for g in ids])
            if totals and lab in totals:
                x = x * 0 + 1
                x[0] = totals[lab]
            idm = bulk.IdMap("symbol", ref["genes"], None)
            xp, pres, _ = idm.to_panel(ids, x)
            T = bulk.library_total(ids, x)
            sc_ = {g: idm.lookup(ids, x, g) for g in cm.SEX_GENES}
            libs[lab] = {"label": lab, "pin": {"stage": 8.5}, "ids": ids, "x": x, "xp": xp, "T": T,
                         "panel_total": float(xp.sum()), "sex_recorded": sex,
                         "sex_call": bulk.sex_call(sc_["Xist"], [sc_[g] for g in cm.SEX_Y_GENES], T),
                         "somites": (somites or {}).get(lab, 12)}
            u[arm].append(lab)
    return libs, u


def _run_unit(ref, libs, u, chrom=None, d=1.0):
    genes = ref["genes"]
    chrom = chrom or [None] * len(genes)
    present = np.ones(len(genes), bool)
    always = np.array([g in cm.ALWAYS_NAN_BULK or g in cm.SEX_Y_GENES or c == "Y" for g, c in
                       zip(genes, chrom, strict=True)])
    chrx = np.array([c == "X" for c in chrom])
    return bulk.bulk_unit({"d": d}, u, libs, ref, present, always, chrx, all(c is not None for c in chrom),
                          bulk.IdMap("symbol", genes, None))


def test_g1_drops_small_libraries_and_an_arm_under_two_is_left_out(tmp_path):
    _, _, ref = make_ref(tmp_path)
    libs, u = _unit_libs(ref, np.random.default_rng(3), n_ctl=2)
    libs["control1"]["T"] = 900_000
    out = _run_unit(ref, libs, u)
    assert out["stages"]["E8.5"]["dropped"] == {"control1": "G1 library size"}
    assert out["status"] == "left_out" and "_L" not in out


def test_g3_imbalance_sets_chrx_not_measured_or_refuses_without_chromosomes(tmp_path):
    _, _, ref = make_ref(tmp_path)
    sexes = {"ko0": "male", "ko1": "male", "ko2": "male"}
    libs, u = _unit_libs(ref, np.random.default_rng(4), sexes=sexes)
    out = _run_unit(ref, libs, u)
    assert out["status"] == "refused_sex_chromosomes"
    chrom = ["X" if g in ("g5", "g6") else ("Y" if g == "g7" else "1") for g in ref["genes"]]
    out = _run_unit(ref, libs, u, chrom=chrom)
    assert out["status"] == "ok" and out["stages"]["E8.5"]["sex"]["balance"] == "imbalanced"
    L = out["_L"]
    assert np.isnan(L[5]) and np.isnan(L[6]) and np.isnan(L[7]) and np.isfinite(L[8])
    libs, u = _unit_libs(ref, np.random.default_rng(4))
    out = _run_unit(ref, libs, u, chrom=chrom)
    assert out["stages"]["E8.5"]["sex"]["balance"] == "balanced"
    assert np.isfinite(out["_L"][5]) and np.isnan(out["_L"][7])   # chrY always, chrX only on imbalance


def test_g4_drops_an_outlier_only_when_both_conditions_hold_and_refuses_swapped_labels(tmp_path):
    rng = np.random.default_rng(0)
    base = rng.uniform(1, 10, 50)
    Y = np.vstack([base + rng.normal(0, 0.05, 50) for _ in range(5)])
    keep, recs = bulk.outlier_drop(Y)
    assert all(keep)
    Y2 = Y.copy()
    Y2[4] = rng.uniform(1, 10, 50)        # a library of another kind: far under the median
    keep, recs = bulk.outlier_drop(Y2)
    assert keep == [True, True, True, True, False]
    Y3 = Y.copy()
    Y3[4] = base + rng.normal(0, 0.2, 50)   # robust z very low (tiny MAD) but within 0.05 of the median: kept
    keep, recs = bulk.outlier_drop(Y3)
    assert recs[4]["rho"] > np.median([r["rho"] for r in recs]) - cm.OUTLIER_RHO_GAP and keep[4]
    assert bulk.outlier_drop(Y[:2])[0] == [True, True]
    other = rng.uniform(1, 10, 50)
    ko = np.vstack([base + rng.normal(0, 0.05, 50) for _ in range(2)])
    ctl_mixed = np.vstack([base + rng.normal(0, 0.05, 50), other + rng.normal(0, 0.05, 50)])
    assert bulk.control_grouping(ko, ctl_mixed)[0] is False
    assert bulk.control_grouping(ko, np.vstack([base + rng.normal(0, 0.05, 50) for _ in range(2)]))[0] is True


def test_g5_a_control_of_one_group_fails_a_pseudo_bulk_passes_and_knockouts_are_recorded_only(tmp_path):
    _, _, ref = make_ref(tmp_path)
    pres = np.ones(len(PANEL), bool)
    assert bulk.whole_embryo_guard(tb.whole_den(ref), ref, pres)[0]
    assert not bulk.whole_embryo_guard(ref["mu"]["Neural"], ref, pres)[0]
    # arms of two (no outlier test): a knockout of one group's profile is recorded and kept; such a control is dropped
    libs, u = _unit_libs(ref, np.random.default_rng(5), n_ko=2, n_ctl=2)
    libs["ko0"]["xp"] = ref["mu"]["CM"] * 300_000
    libs["control1"]["xp"] = ref["mu"]["Neural"] * 300_000
    out = _run_unit(ref, libs, u)
    st = out["stages"]["E8.5"]
    assert st["dropped"] == {"control1": "G5 dominated by one lineage (not more like pb than any one group)"}
    rec = st["whole_embryo"]["control1"]
    assert rec["threshold"] == cm.BULK_WHOLE_MARGIN == 0.0 and rec["margin"] <= 0.0 and rec["passed"] is False
    assert rec["margin"] == pytest.approx(rec["rho_whole"] - rec["rho_best_group"])
    assert st["whole_embryo"]["ko0"]["passed"] is False and "ko0" in st["libraries"]["ko"]
    assert out["status"] == "left_out"
    assert out["g5_rule"] == "bulk margin > 0 (BULK_WHOLE_MARGIN), decision 2026-11-21"


def test_g5_two_lineage_mix_fails_and_a_library_on_each_side_of_zero(tmp_path):
    """Bulk G5 (decision 2026-11-21): margin > 0, strict. A 50/50 mix of the two smallest-share groups fails; along
    the line from pb to one group's profile there is a profile with margin in (0, 0.05), which passes (it failed
    under the old 0.05), and one with margin under 0, which fails."""
    _, _, ref = make_ref(tmp_path)
    pres = np.ones(len(PANEL), bool)
    small = sorted(ref["groups"], key=lambda g: ref["shares"][g])[:2]
    assert small == ["CM", "Endoderm"]
    ok, rec = bulk.whole_embryo_guard(0.5 * ref["mu"]["CM"] + 0.5 * ref["mu"]["Endoderm"], ref, pres)
    assert ok is False and rec["margin"] <= 0.0
    pb = tb.whole_den(ref)
    above = below = None
    for t in np.linspace(0.0, 1.0, 2001):
        ok, rec = bulk.whole_embryo_guard((1 - t) * pb + t * ref["mu"]["Neural"], ref, pres)
        if above is None and 0.0 < rec["margin"] < cm.SPEARMAN_MARGIN:
            above = (ok, rec)
        if below is None and rec["margin"] < 0.0:
            below = (ok, rec)
    assert above is not None and below is not None
    assert above[0] is True and above[1]["threshold"] == 0.0
    assert below[0] is False


def test_g5_arm_mean_refuses_when_it_fails_and_passes_inside_0_0_05(tmp_path, monkeypatch):
    """The control arm mean under the new rule (the reviewer's lesser point 3). Two controls with one profile whose
    margin lies in (0, 0.05): both are kept (the old 0.05 would have dropped them) and the arm mean passes. An arm
    mean that fails refuses the unit with the new reason."""
    _, _, ref = make_ref(tmp_path)
    pres = np.ones(len(PANEL), bool)
    pb = tb.whole_den(ref)
    prof = None
    for t_ in np.linspace(0.0, 1.0, 2001):
        cand = (1 - t_) * pb + t_ * ref["mu"]["Neural"]
        ok, rec = bulk.whole_embryo_guard(cand, ref, pres)
        if 0.0 < rec["margin"] < cm.SPEARMAN_MARGIN:
            prof = cand
            break
    assert prof is not None
    libs, u = _unit_libs(ref, np.random.default_rng(5), n_ko=2, n_ctl=2)
    for lab in ("control0", "control1"):
        libs[lab]["xp"] = prof / prof.sum() * 300_000
    out = _run_unit(ref, libs, u)
    st = out["stages"]["E8.5"]
    assert not any("G5" in v for v in st["dropped"].values())
    assert st["whole_embryo"]["control0"]["passed"] and st["whole_embryo"]["control1"]["passed"]
    assert 0.0 < st["whole_embryo_control_mean"]["margin"] < cm.SPEARMAN_MARGIN
    assert st["whole_embryo_control_mean"]["threshold"] == 0.0 and out["status"] != "refused_whole_embryo"
    # an arm mean that fails: the guard's verdict on the fifth call (two knockouts, two controls, then the mean)
    real, calls = bulk.whole_embryo_guard, []

    def guard(profile, ref_, present):
        calls.append(1)
        ok, rec = real(profile, ref_, present)
        return (False, {**rec, "margin": -0.01}) if len(calls) == 5 else (ok, rec)

    monkeypatch.setattr(bulk, "whole_embryo_guard", guard)
    libs2, u2 = _unit_libs(ref, np.random.default_rng(5), n_ko=2, n_ctl=2)
    out2 = _run_unit(ref, libs2, u2)
    assert len(calls) == 5 and out2["status"] == "refused_whole_embryo"
    assert out2["reason"] == ("E8.5: the control arm's mean is dominated by one lineage (not more like pb than any one "
                              "group) (G5)")


def test_g5_tie_at_zero_fails_under_the_strict_rule(tmp_path):
    """A one-group reference: pb equals that group, so every profile's margin is exactly 0 and fails."""
    _, _, ref = make_ref(tmp_path)
    one = {**ref, "groups": ["CM"], "shares": {"CM": 1.0}, "mu": {"CM": ref["mu"]["CM"]}}
    pres = np.ones(len(PANEL), bool)
    for prof in (one["mu"]["CM"], ref["mu"]["Neural"], tb.whole_den(ref)):
        ok, rec = bulk.whole_embryo_guard(prof, one, pres)
        assert rec["margin"] == 0.0 and ok is False


def test_the_single_cell_group_check_keeps_spearman_margin_and_the_bulk_number_is_separate(tmp_path):
    """process.spearman_guard is unchanged: a control profile whose own group leads the next by a margin in
    (0, 0.05) fails it. A bulk G5 profile with a margin in the same interval passes. The two numbers are separate."""
    assert cm.SPEARMAN_MARGIN == 0.05 and cm.BULK_WHOLE_MARGIN == 0.0
    _, _, ref = make_ref(tmp_path)
    pres = np.ones(len(PANEL), bool)
    sc = None
    for t in np.linspace(0.0, 1.0, 2001):
        ok, rec = pr.spearman_guard((1 - t) * ref["mu"]["CM"] + t * ref["mu"]["Endoderm"], ref, "CM", pres)
        if rec.get("best") == "CM" and "margin" in rec and 0.0 < rec["margin"] < cm.SPEARMAN_MARGIN:
            sc = (ok, rec)
            break
    assert sc is not None and sc[0] is False
    pb = tb.whole_den(ref)
    bk = None
    for t in np.linspace(0.0, 1.0, 2001):
        ok, rec = bulk.whole_embryo_guard((1 - t) * pb + t * ref["mu"]["Neural"], ref, pres)
        if 0.0 < rec["margin"] < cm.SPEARMAN_MARGIN:
            bk = (ok, rec)
            break
    assert bk is not None and bk[0] is True


def test_g6_non_overlapping_somites_refuse(tmp_path):
    assert bulk.somite_record([10, 11], [14, 15])[0] is False
    ok, rec = bulk.somite_record([10, 14], [12, 15])
    assert ok and rec["median_difference"] == -1.5 and rec["ko_within_control_range"] == 1
    assert bulk.somite_record([None], [12])[1] == {"recorded": False}
    _, _, ref = make_ref(tmp_path)
    libs, u = _unit_libs(ref, np.random.default_rng(6), somites={"ko0": 8, "ko1": 8, "ko2": 9})
    assert _run_unit(ref, libs, u)["status"] == "refused_somites"


def test_g12_a_threefold_shift_of_a_main_group_refuses_a_small_group_passes(tmp_path):
    _, _, ref = make_ref(tmp_path)
    pres = np.ones(len(PANEL), bool)

    def mix(w):
        p = sum(w[g] * ref["mu"][g] for g in GROUPS)
        return p / p.sum()

    ctl = np.vstack([mix({"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5})] * 2)
    ko = np.vstack([mix({"CM": 0.6, "Endoderm": 0.15, "Neural": 0.25})] * 2)
    assert bulk.composition_guard(ko, ctl, ref, pres)[0] is False
    ctl = np.vstack([mix({"CM": 0.01, "Endoderm": 0.39, "Neural": 0.6})] * 2)
    ko = np.vstack([mix({"CM": 0.03, "Endoderm": 0.38, "Neural": 0.59})] * 2)
    ok, rec = bulk.composition_guard(ko, ctl, ref, pres)
    assert ok and rec["shares_control"]["CM"] < cm.COMPOSITION_MIN_SHARE


def test_g8_identical_vectors(tmp_path):
    a = {"label": "a", "vector_sha256": "x"}
    assert bulk.duplicate_vectors([a, {"label": "b", "vector_sha256": "x"}, {"label": "c", "vector_sha256": "y"}]) == [
        ["a", "b"]]


# ---------------------------------------------------------------------------------------------------------------------
# ids


def test_the_id_routes(tmp_path, monkeypatch):
    _, _, ref = make_ref(tmp_path)
    pfile, psha = panel_file(tmp_path, ref)
    panel = cm.load_panel_ids(ref["genes"], pfile, psha)
    idm = bulk.IdMap("ensembl", ref["genes"], panel)
    ids = [ID["g0"] + ".7", ID["g1"], "ENSMUSG99999999999", "__no_feature"]
    xp, pres, rec = idm.to_panel(ids, np.array([5.0, 6, 7, 8]))
    assert xp[0] == 5 and xp[1] == 6 and pres.sum() == 2 and rec["route"] == "ensembl"
    assert idm.lookup(ids, np.array([5.0, 6, 7, 8]), "g1") == 6
    assert idm.lookup(ids, np.array([5.0, 6, 7, 8]), "Uty") is None
    with pytest.raises(cm.Refused, match="repeat a stripped"):
        bulk.IdMap("ensembl", ref["genes"], panel).to_panel([ID["g0"] + ".1", ID["g0"] + ".2"], np.array([1.0, 2]))
    with pytest.raises(cm.Refused, match="not the reference"):
        cm.load_panel_ids(list(reversed(ref["genes"])), pfile, psha)
    with pytest.raises(cm.Refused, match="sha256"):
        cm.load_panel_ids(ref["genes"], pfile, "0" * 64)
    with monkeypatch.context() as mp:   # the unpinned state (the committed pin: test_t3_ext_panel_ids.py)
        mp.setattr(cm, "PANEL_IDS_SHA256", None)
        with pytest.raises(cm.Refused, match="None"):
            cm.load_panel_ids(ref["genes"])
    # an id naming two panel symbols: NA for both in make_panel_ids; a hand-made file with it refuses
    text = pfile.read_text().replace(f"g2\t{ID['g2']}", f"g2\t{ID['g1']}")
    bad = tmp_path / "bad.tsv"
    bad.write_text(text)
    with pytest.raises(cm.Refused, match="two panel symbols"):
        cm.load_panel_ids(ref["genes"], bad, cm.sha256_file(bad))
    built = mpi.build([(ID["g1"], "g1"), (ID["g1"], "g2"), ("E5", "g3"), ("E6", "g3")], ref["genes"], [], "s", "0" * 64)
    lines = {ln.split("\t")[0]: ln.split("\t") for ln in built.splitlines() if not ln.startswith("#")}
    assert lines["g1"][1] == lines["g2"][1] == "NA" and lines["g3"][3] == "2 ids in the source"
    assert lines["g4"][1:] == ["NA", "NA", "no id in the source"]
    # an NA in the map: the gene is not measured
    na = panel["ids"].index(ID["g0"])
    panel2 = {**panel, "ids": [None if i == na else v for i, v in enumerate(panel["ids"])]}
    assert not bulk.IdMap("ensembl", ref["genes"], panel2).to_panel([ID["g0"]], np.array([5.0]))[1][0]
    # the symbol route is process.to_panel's rule
    syms = ["g0", "g1", "g1", "T", "zz"]
    x = np.array([1.0, 2, 3, 4, 5])
    genes = ["g0", "g1", "Tbxt", "g9"]
    a = bulk.IdMap("symbol", genes, None).to_panel(syms, x)
    from scipy import sparse

    b = pr.to_panel(sparse.csr_matrix(x[None, :]), syms, ["Gene Expression"] * 5, genes)
    assert a[0].tolist() == b[0].toarray().ravel().tolist() and a[1].tolist() == b[1].tolist()


SC_YAML = """version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
  - accession: GSE930001
    status: cleared
    tier: unrelated
    study: GSE930001
    design: germline
    control_relation: littermate
    repository: geo
    screen_series: [GSE930001]
    stages_used: [E8.5]
    pinned_series: [{accession: GSE930001, series_metadata_sha256: "H0"}]
    pinned:
      - {accession: GSM1, arm: ko, gene: X, stage: E8.5, files: [GSM1_features.tsv.gz], metadata_sha256: "H1"}
      - {accession: GSM2, arm: control, gene: X, stage: E8.5, files: [GSM2_features.tsv.gz], metadata_sha256: "H2"}
""".replace("H0", "0" * 64).replace("H1", "1" * 64).replace("H2", "2" * 64)


def test_make_panel_ids_reads_the_features_of_a_pinned_fetched_sample_only(tmp_path):
    y2 = tmp_path / "sc.yaml"
    y2.write_text(SC_YAML)
    dest = tmp_path / "dest"
    f = dest / "raw" / "GSE930001" / "GSM1" / "GSM1_features.tsv.gz"
    f.parent.mkdir(parents=True)
    f.write_bytes(gzip.compress(b"E1\tg0\tGene Expression\nE2\tg1\tGene Expression\n"))
    rec = {"sha256": cm.sha256_file(f), "pins_sha256": cm.pins_sha256(cm.dataset("GSE930001", y2))}
    (dest / "raw" / "GSE930001" / "MANIFEST.json").write_text(json.dumps({"files": {"GSM1/GSM1_features.tsv.gz": rec}}))
    assert mpi.source_file(dest, "GSE930001", "GSM1", "GSM1_features.tsv.gz", y2) == f
    with pytest.raises(cm.Refused, match="not pinned"):
        mpi.source_file(dest, "GSE930001", "GSM1", "GSM1_barcodes.tsv.gz", y2)
    with pytest.raises(cm.Refused, match="current pins"):
        mpi.source_file(dest, "GSE930001", "GSM2", "GSM2_features.tsv.gz", y2)
    f.write_bytes(gzip.compress(b"E1\tg0\tGene Expression\n"))   # changed after the fetch
    with pytest.raises(cm.Refused, match="current pins"):
        mpi.source_file(dest, "GSE930001", "GSM1", "GSM1_features.tsv.gz", y2)
    with pytest.raises(cm.Refused, match="single-cell"):
        mpi.source_file(dest, "GSE910001", "GSM910001", "x", _yaml(tmp_path))


def _yaml(tmp_path):
    p = tmp_path / "b.yaml"
    p.write_text(YAML)
    return p


def test_make_panel_ids_takes_only_a_pinned_fetched_file(tmp_path):
    with pytest.raises(cm.Refused):
        mpi.source_file(tmp_path, "GSE168768", "GSM5169142", "x.tsv.gz")
    with pytest.raises(cm.Refused, match="not a features"):
        p = tmp_path / "GSM1_matrix.mtx.gz"
        p.write_bytes(b"")
        mpi.features_rows(p)
    p = tmp_path / "GSM1_features.tsv.gz"
    p.write_bytes(gzip.compress(b"E1\tg0\tGene Expression\nE2\tg1\tAntibody Capture\n"))
    assert mpi.features_rows(p) == ([("E1", "g0")], "")


# ---------------------------------------------------------------------------------------------------------------------
# the dataset keys and pins


@pytest.mark.parametrize("edit,match", [
    (("    count_form: per_sample\n", "    count_form: xlsx\n"), "count_form"),
    (("    id_type: symbol\n", "    id_type: entrez\n"), "id_type"),
    (("    modality: bulk\n    count_form: per_sample\n", "    modality: bulky\n    count_form: per_sample\n"),
     "modality"),
    (("    gene_columns: [gene]\n", "    gene_columns: []\n"), "gene_columns"),
    (("    gene_columns: [gene]\n", "    gene_columns: [gene]\n    cell_labels: true\n"), "single-cell"),
    (("    gene_columns: [gene]\n", "    gene_columns: [gene]\n    withheld_lines: []\n"), "arrayexpress only"),
    (("    gene_columns: [gene]\n", "    gene_columns: [gene]\n    stage_from_series: title\n"), "series_table"),
    (("    count_form: per_sample\n", "    count_form: arrayexpress\n"), "arrayexpress"),
])
def test_the_bulk_keys_are_checked(tmp_path, edit, match):
    p = tmp_path / "y.yaml"
    p.write_text(YAML.replace(*edit, 1))
    with pytest.raises(cm.Refused, match=match):
        cm.load_datasets(p)


def test_a_single_cell_dataset_holds_no_bulk_key_and_its_pins_hash_as_before(tmp_path):
    p = tmp_path / "y.yaml"
    p.write_text(YAML.replace("    modality: bulk\n    count_form: per_sample\n", "    count_form: per_sample\n", 1))
    with pytest.raises(cm.Refused, match="bulk keys"):
        cm.load_datasets(p)
    ds = {"pinned": [{"accession": "GSM1"}], "pinned_series": [{"accession": "GSE1"}]}
    want = hashlib.sha256(json.dumps({"samples": ds["pinned"], "series": ds["pinned_series"]}, sort_keys=True,
                                     ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
    assert cm.pins_sha256(ds) == want
    # the committed file: every dataset single-cell, as before, but the one many-line bulk resource added on 2026-11-20
    # (docs/PLAN.md decision log 2026-11-20 (eee)); no GEO bulk dataset (those forms are closed)
    bulk = {acc for acc, d in cm.load_datasets().items() if cm.is_bulk(d)}
    assert bulk == {"E-MTAB-7322"}
    for acc, d in cm.load_datasets().items():
        if acc not in bulk:
            assert cm.modality(d) == "single_cell"


def test_bulk_pins_need_two_libraries_per_arm_and_their_own_form(tmp_path, y):
    smp = {g: sample(g, geno) for g, geno in PS.items()}
    metas = {"GSE910001": fe.parse_soft(soft("GSE910001", smp))}
    dest = tmp_path / "d"
    fe.plan("GSE910001", dest, y, metas=metas)
    pins = bulk.pin_template("GSE910001", dest, "ko", KO, ["GSM910001"], y) + bulk.pin_template(
        "GSE910001", dest, "control", KO, ["GSM910004", "GSM910005"], y)
    paste(y, "GSE910001", pins, fe.series_pin_template("GSE910001", dest, y))
    with pytest.raises(cm.Refused, match="fewer than 2"):
        cm.load_datasets(y)
    y.write_text(YAML)
    pins = bulk.pin_template("GSE910001", dest, "ko", KO, ["GSM910001", "GSM910002"], y) + bulk.pin_template(
        "GSE910001", dest, "control", KO, ["GSM910004", "GSM910005"], y)
    paste(y, "GSE910001", pins.replace("        files:", "        column: x\n        files:", 1),
          fe.series_pin_template("GSE910001", dest, y))
    with pytest.raises(cm.Refused, match="exactly the keys"):
        cm.load_datasets(y)


def test_a_held_superseries_relation_refuses_the_plan(tmp_path, y):
    smp = {g: sample(g, geno) for g, geno in PS.items()}
    metas = {"GSE910001": fe.parse_soft(soft("GSE910001", smp, relation=("SubSeries of: GSE259342",)))}
    pl = fe.plan("GSE910001", tmp_path / "d", y, metas=metas)
    assert pl["status"] == "refused" and "GSE259342" in pl["reason"]


def test_the_header_step_and_ae_refuse_held_accessions_and_no_cli_takes_another_list(tmp_path):
    for acc in ("GSE165300", "E-MTAB-7325"):
        with pytest.raises(cm.Refused, match="NEVER_FETCH"):
            bulk.header(acc, tmp_path, "GSE1", "x.txt")
    with pytest.raises(SystemExit) as e:
        fe.main(["header", "GSE910002", "--dest", str(tmp_path), "--datasets", "x.yaml"])
    assert e.value.code == 2
    assert fe.main(["header", "GSE165300", "--dest", str(tmp_path)]) == 2
    assert fe.main(["table-template", "GSE165300", "--dest", str(tmp_path)]) == 2
    with pytest.raises(cm.Refused, match="inside the repository"):
        bulk.write_bulk_unit(EXT / "out", "A", "B", np.zeros(1), {"genes": ["a"], "file_sha256": "x", "shares": {},
                                                                   "groups": [], "mu": {}}, {})


def test_no_t1_or_t2_builder_imports_the_bulk_code():
    files = [*(ROOT / "src" / "vec").glob("*.py"), *(ROOT / "experiments" / "builds").glob("*.py")]
    for d in sorted((ROOT / "experiments").glob("t[12]_*")):
        files += [*d.glob("*.py")]
    for f in files:
        t = f.read_text(errors="replace")
        assert "import bulk" not in t and "import ae" not in t and "make_panel_ids" not in t, f


# ---------------------------------------------------------------------------------------------------------------------
# combine: the bulk family, its sign rule, and the single-cell family untouched


def _old_sign_rule_mean(V, n_table_studies):
    """sign_rule_mean as it was before stage 9 (copied), for the equality test."""
    V = np.atleast_2d(np.asarray(V, dtype=np.float64))
    held = np.isfinite(V)
    n = held.sum(0)
    m = np.where(n > 0, np.nansum(np.where(held, V, 0.0), 0) / np.maximum(n, 1), np.nan)
    sgn = np.sign(m)
    same = (held & (np.sign(np.where(held, V, 0.0)) == sgn) & (np.where(held, V, 0.0) != 0)).sum(0)
    if n_table_studies < 2:
        return m, same, n
    opp = np.where(held & (np.sign(np.where(held, V, 0.0)) == -sgn), np.abs(np.where(held, V, 0.0)), 0.0).max(0)
    keep = (same >= 2) & (opp <= np.abs(m)) & (sgn != 0)
    return np.where(n > 0, np.where(keep, m, 0.0), np.nan), same, n


def test_the_two_study_rule_is_unchanged_on_a_grid():
    rng = np.random.default_rng(9)
    for n in range(1, 8):
        V = rng.normal(0, 1, (n, 200))
        V[rng.uniform(size=V.shape) < 0.3] = np.nan
        V[rng.uniform(size=V.shape) < 0.1] = 0.0
        for a, b in zip(co.sign_rule_mean(V, n), _old_sign_rule_mean(V, n), strict=True):
            assert np.array_equal(a, b, equal_nan=True)


def test_the_two_thirds_rule():
    def keep(vals):
        V = np.array(vals, dtype=float)[:, None]
        return co.sign_rule_mean(V, len(vals), "two_thirds")[0][0]

    assert keep([0.3, 0.2, -0.1]) == pytest.approx(0.4 / 3)          # h = 3, 2 agree: kept
    assert keep([0.3, 0.2, -0.1, -0.05]) == 0.0                       # h = 4, 2 agree: 0
    assert keep([1] * 6 + [-0.1] * 4) == 0.0                          # h = 10, 6 agree: 0
    assert keep([1] * 7 + [-0.1] * 3) != 0.0                          # h = 10, 7 agree: kept
    assert keep([0.3, 0.0, 0.0]) == 0.0                               # exact zeros do not agree
    assert keep([0.3, 0.2, np.nan, np.nan]) == pytest.approx(0.25)   # h = 2 of a 4-study table: 2 agree
    with pytest.raises(cm.Refused):
        co.sign_rule_mean(np.zeros((3, 1)), 3, "median")


def test_write_table_with_per_group_is_byte_identical(tmp_path):
    genes = ["a", "b", "c"]
    phi = np.array([1.0, 2.0, 0.5])
    pg = {"CM": np.array([1.5, 1.0, 1.0])}
    files = tb.write_table(tmp_path / "t", "n", genes, phi, pg, {"x": 1})
    buf = io.StringIO()
    pd.DataFrame({"gene": genes, "factor": phi}).to_csv(buf, index=False)
    assert Path(files["n"]["file"]).read_text() == buf.getvalue()
    buf = io.StringIO()
    pd.DataFrame([("CM", g, float(v)) for g, v in zip(genes, pg["CM"], strict=True)],
                 columns=["group", "gene", "factor"]).to_csv(buf, index=False)
    assert Path(files["n_per_group"]["file"]).read_text() == buf.getvalue()
    man = json.loads((tmp_path / "t" / "n_manifest.json").read_text())
    assert "per_group" not in man and list(man["files"]) == ["n", "n_per_group"]
    files = tb.write_table(tmp_path / "u", "m", genes, phi, None, {"x": 1})
    assert list(files) == ["m"] and not (tmp_path / "u" / "m_per_group.csv").exists()


def _units(ref, rng, studies, mod="bulk", per=1):
    out = []
    for s in studies:
        for k in range(per):
            L = rng.normal(0, 0.2, len(PANEL))
            u = {"accession": f"A{s}", "gene": f"G{k}", "study": s, "tier": "unrelated", "sha256": "x",
                 "L": {"__whole__": L} if mod == "bulk" else {"CM": L, "Neural": L * 0.5}}
            if mod == "bulk":
                u.update(modality="bulk", k={"E8.5": 0.5}, libraries={"ko": ["a"], "control": ["b"]}, guards={})
            out.append(u)
    return out


def test_the_bulk_family_and_the_single_cell_bytes(tmp_path, monkeypatch):
    _, refp, ref = make_ref(tmp_path)
    rng = np.random.default_rng(1)
    sc_units = _units(ref, rng, ["S1", "S2", "S3"], mod="sc")
    b_units = _units(ref, rng, ["B1", "B2", "B3"])
    monkeypatch.setattr(co, "load_units", lambda *a, **k: (list(state["units"]), {"left_out": [], "not_cleared": [],
                                                                                   "inputs": {}}))
    state = {"units": sc_units}
    co.combine(tmp_path, refp, tmp_path / "sc_only", dest=tmp_path)
    state["units"] = sc_units + b_units
    r = co.combine(tmp_path, refp, tmp_path / "both", dest=tmp_path)
    sc_names = sorted(p.name for p in (tmp_path / "sc_only").iterdir() if p.is_dir())
    for name in sc_names:
        for f in (tmp_path / "sc_only" / name).rglob("*.csv"):
            assert f.read_bytes() == (tmp_path / "both" / f.relative_to(tmp_path / "sc_only")).read_bytes()
    bulk_names = {n for n in r["tables"] if n.startswith("bulk_")}
    assert bulk_names == {"bulk_all", "bulk_study_B1", "bulk_study_B2", "bulk_study_B3", "bulk_all__loso_B1",
                          "bulk_all__loso_B2", "bulk_all__loso_B3"}
    assert not any(n.startswith("both") for n in r["tables"]) and "all" in r["tables"]
    assert not [p for p in (tmp_path / "both").rglob("bulk_*_per_group.csv")]
    man = json.loads((tmp_path / "both" / "bulk_all" / "bulk_all_manifest.json").read_text())
    assert man["sign_rule"] == "two_thirds" and man["studies"] == ["B1", "B2", "B3"]
    man2 = json.loads((tmp_path / "both" / "bulk_all__loso_B1" / "bulk_all__loso_B1_manifest.json").read_text())
    assert man2["sign_rule"] == "two_study"
    man3 = json.loads((tmp_path / "both" / "all" / "all_manifest.json").read_text())
    assert man3["sign_rule"] == "applied" and man3["studies"] == ["S1", "S2", "S3"]
    # the bulk table is finished against the whole-embryo pseudo-bulk
    L = co.sign_rule_mean(np.vstack([u["L"]["__whole__"] for u in b_units]), 3, "two_thirds")[0]
    got = pd.read_csv(tmp_path / "both" / "bulk_all" / "factors" / "bulk_all.csv")["factor"].to_numpy()
    assert np.allclose(got, tb.finish(L, tb.whole_den(ref)))
    # a mixed-modality study refuses
    state["units"] = sc_units + _units(ref, rng, ["S1"])
    with pytest.raises(cm.Refused, match="span modalities"):
        co.combine(tmp_path, refp, tmp_path / "mixed", dest=tmp_path)


def test_a_fifty_unit_study_equals_the_mean_of_its_units(tmp_path):
    _, _, ref = make_ref(tmp_path)
    us = _units(ref, np.random.default_rng(2), ["MGP"], per=50)
    st = co.bulk_study_tables(us, len(PANEL))
    assert np.allclose(st["MGP"]["L"]["__whole__"], np.mean([u["L"]["__whole__"] for u in us], 0))
    assert len(st["MGP"]["line_agreement"]) == len(PANEL)


def test_load_units_refuses_a_whole_key_on_a_single_cell_unit_and_a_shared_count_file(tmp_path, y, monkeypatch):
    dest, wt, refp, ref, _ = run_per_sample(tmp_path, y, monkeypatch)
    out = tmp_path / "processed"
    pr.run("GSE910001", dest, wt, refp, out, y)
    qcp = out / "GSE910001" / "QC.json"
    qc = json.loads(qcp.read_text())
    up = Path(qc["units"][KO]["unit_file"])
    un = json.loads(up.read_text())
    un["L"] = {"CM": un["L"]["__whole__"]}
    up.write_text(json.dumps(un))
    qc["units"][KO]["unit_sha256"] = cm.sha256_file(up)
    qcp.write_text(json.dumps(qc))
    with pytest.raises(cm.Refused, match="__whole__"):
        co.load_units(out, y, ref, dest)
    un["modality"] = "single_cell"
    up.write_text(json.dumps(un))
    qc["units"][KO]["unit_sha256"] = cm.sha256_file(up)
    qcp.write_text(json.dumps(qc))
    with pytest.raises(cm.Refused, match="single_cell unit under a bulk dataset"):
        co.load_units(out, y, ref, dest)
    # G8 across datasets: the same count file under a second bulk dataset's QC
    monkeypatch.setattr(co, "check_plan", lambda *a, **k: None)
    (out / "GSE910002").mkdir()
    (out / "GSE910002" / "QC.json").write_text(json.dumps({
        "accession": "GSE910002", "status": "ok", "modality": "bulk", "pins_sha256": cm.pins_sha256(
            {**cm.dataset("GSE910002", y), "pinned": [{"accession": "GSM1"}]}), "plan_sha256": "x",
        "count_files_sha256": qc["count_files_sha256"], "units": {}}))
    un["modality"] = "bulk"
    un["L"] = {"__whole__": un["L"]["CM"]}
    up.write_text(json.dumps(un))
    qc["units"][KO]["unit_sha256"] = cm.sha256_file(up)
    qcp.write_text(json.dumps(qc))
    monkeypatch.setattr(cm, "load_datasets", _with_pins(cm.load_datasets))
    with pytest.raises(cm.Refused, match="G8"):
        co.load_units(out, y, ref, dest)


def _with_pins(f):
    def g(*a, **k):
        d = f(*a, **k)
        d["GSE910002"]["pinned"] = [{"accession": "GSM1"}]
        return d
    return g


# ---------------------------------------------------------------------------------------------------------------------
# stage 9 round 2: each blocker of the first review of 120af2f, replayed on the engineer's own synthetic studies


@pytest.mark.parametrize("hdr,term", [("Gata4_KO_1", None), ("Ctnnb1_null_rep2", None), ("Axin2_KO", "axin2")])
def test_r2_blocker2_a_per_sample_header_is_screened_before_any_value_is_read(tmp_path, y, monkeypatch, hdr, term):
    dest, wt, refp, ref, _ = run_per_sample(tmp_path, y, monkeypatch, cols={"GSM910002": hdr})
    qc = pr.run("GSE910001", dest, wt, refp, tmp_path / "processed", y)
    assert qc["status"] == "failed"
    stored = (tmp_path / "processed" / "GSE910001" / "QC.json").read_text()
    if term is None:   # a held-out hit: the text is withheld
        assert "held-out hit in a per-sample header (text withheld)" in qc["reason"] and hdr not in stored
        assert hdr.split("_")[0] not in stored
    else:
        assert "per-sample header screen" in qc["reason"] and "review" in qc["reason"]


def test_r2_blocker2_a_review_term_in_a_per_sample_header_is_cleared_only_by_an_entry_naming_that_sample(
        tmp_path, y, monkeypatch):
    for who, status in (("GSM910003", "failed"), ("GSM910002", "ok")):
        y.write_text(YAML.replace("    somite_field: \"characteristics_ch1:somites\"\n",
                                  "    somite_field: \"characteristics_ch1:somites\"\n    review_cleared: [{sample: "
                                  f"{who}, term: axin2, reason: a header read}}]\n"))
        d = tmp_path / who
        d.mkdir()
        dest, wt, refp, ref, _ = run_per_sample(d, y, monkeypatch, cols={"GSM910002": "Axin2_KO"})
        qc = pr.run("GSE910001", dest, wt, refp, d / "processed", y)
        assert qc["status"] == status, qc.get("reason")


def test_r2_blocker3_a_single_count_column_headed_by_another_sample_is_refused(tmp_path, y, monkeypatch):
    # the reviewer's probe: every file's one column headed 'GSM910006' (a control's accession)
    dest, wt, refp, ref, _ = run_per_sample(tmp_path, y, monkeypatch, cols={g: "GSM910006" for g in PS})
    qc = pr.run("GSE910001", dest, wt, refp, tmp_path / "processed", y)
    assert qc["status"] == "failed" and "ties by the exact rule to another sample" in qc["reason"]
    # a title or a characteristics value of another sample, likewise
    for col in ("GSM910005", "Kozz -/-"):
        d = tmp_path / col.replace(" ", "_").replace("/", "_")
        d.mkdir()
        y.write_text(YAML)
        dest, wt, refp, ref, _ = run_per_sample(d, y, monkeypatch, cols={"GSM910001": col})
        qc = pr.run("GSE910001", dest, wt, refp, d / "processed", y)
        assert qc["status"] == "failed" and "another sample" in qc["reason"], qc.get("reason")


@pytest.mark.parametrize("short", ["normalized_counts.txt.gz", "expected_count.tsv", "counts_TPM.txt"])
def test_r2_blocker4_a_series_table_named_as_normalised_is_refused_at_the_offer_header_step_and_read(
        tmp_path, y, monkeypatch, short):
    name = f"GSE910002_{short}"
    dest = tmp_path / "dest"
    metas = {"GSE910002": fe.parse_soft(soft("GSE910002", st_samples(), design="whole embryos at E8.5",
                                              series_files=(short,)))}
    pl = fe.plan("GSE910002", dest, y, metas=metas)
    assert pl["status"] == "ok", pl.get("reason")
    b = pl["series"]["GSE910002"]
    assert b["files"] == [] and "normalised or expected counts" in b["other_files"][0]["reason"]   # the offer
    calls: list = []
    monkeypatch.setattr(fe, "_open", opener(gzip.compress(b"gene_id\tKO_rep1\n1\t2\n"), calls)[0])
    with pytest.raises(cm.Refused, match="normalised or expected counts"):
        bulk.header("GSE910002", dest, "GSE910002", name, y)   # the header step
    assert calls == []
    rows = [r for rs in pl["tables"].values() for r in rs]
    with pytest.raises(cm.Refused, match="normalised or expected counts"):
        bulk.record_header(cm.dataset("GSE910002", y), dest, "GSE910002", name, "u", b"gene_id\tKO_rep1", rows, "x")
    p = tmp_path / name
    p.write_bytes(gzip.compress(b"gene_id\tKO_rep1\ng1\t5\n") if name.endswith(".gz") else b"gene_id\tKO_rep1\ng1\t5\n")
    with pytest.raises(cm.Refused, match="normalised or expected counts"):   # the read
        bulk.read_series_table(p, name, ["gene_id"], ["KO_rep1"], hashlib.sha256(b"gene_id\tKO_rep1").hexdigest())
    ds = {**cm.dataset("GSE910002", y), "pinned_series_files": [{"series": "GSE910002", "name": name,
                                                                 "header_sha256": "0" * 64}]}
    with pytest.raises(cm.Refused, match="normalised or expected counts"):
        bulk._header_check(ds, pl, rows)


def test_r2_blocker5_heterozygous_wording_refuses_a_pin_in_either_arm(tmp_path, y, monkeypatch):
    dest = tmp_path / "dest"
    smp = {g: sample(g, geno) for g, geno in PS.items()}
    smp["GSM910009"] = sample("GSM910009", "Kozz +/-")
    smp["GSM910010"] = sample("GSM910010", "wild type", extra=["note: littermate of a het embryo"])
    smp["GSM910011"] = {**sample("GSM910011", "Kozz -/-"), "title": "Kozz heterozygous embryo 2"}
    metas = {"GSE910001": fe.parse_soft(soft("GSE910001", smp))}
    pl = fe.plan("GSE910001", dest, y, metas=metas)
    assert pl["status"] == "ok", pl.get("reason")
    rows = {r["accession"]: r for r in pl["tables"]["GSE910001"]}
    for g in ("GSM910009", "GSM910010", "GSM910011"):
        assert "heterozygous wording" in rows[g]["not_pinnable"]
    assert rows["GSM910001"]["not_pinnable"] is None
    table = (dest / "meta" / "GSE910001" / "plan_table_GSE910001.txt").read_text()
    assert "heterozygous wording" in table
    # the reviewer's probe: a ko pin of the +/- sample beside a clean one; and the control arm
    with pytest.raises(cm.Refused, match="heterozygous"):
        bulk.pin_template("GSE910001", dest, "ko", KO, ["GSM910001", "GSM910009"], y)
    with pytest.raises(cm.Refused, match="heterozygous"):
        bulk.pin_template("GSE910001", dest, "control", KO, ["GSM910004", "GSM910010"], y)
    with pytest.raises(cm.Refused, match="heterozygous"):
        bulk.pin_template("GSE910001", dest, "ko", KO, ["GSM910011"], y)
    # the pin check: a pin of the het sample, forged past the template, is refused in either arm
    pins = (bulk.pin_template("GSE910001", dest, "ko", KO, list(PS)[:3], y)
            + bulk.pin_template("GSE910001", dest, "control", KO, list(PS)[3:], y))
    paste(y, "GSE910001", pins, fe.series_pin_template("GSE910001", dest, y))
    ds = cm.dataset("GSE910001", y)
    pl = fe.plan("GSE910001", dest, y, metas=metas)
    assert pl["status"] == "ok", pl.get("reason")
    for g, arm in (("GSM910009", "ko"), ("GSM910010", "control")):
        r = {x["accession"]: x for x in pl["tables"]["GSE910001"]}[g]
        forged = {"accession": g, "arm": arm, "gene": KO, "stage": "E8.5", "files": [f"{g}_counts.txt.gz"],
                  "metadata_sha256": r["metadata_sha256"]}
        with pytest.raises(cm.Refused, match="heterozygous"):
            fe.verify_pins({**ds, "pinned": [*ds["pinned"], forged]}, pl)
    # the series-table form: the 'Kozz het' sample titled HET_rep1
    st = tmp_path / "st"
    st.mkdir()
    dest2, *_ = run_series_table(st, y, monkeypatch)
    with pytest.raises(cm.Refused, match="heterozygous"):
        bulk.pin_template("GSE910002", dest2, "ko", KO, ["GSM920005"], y)


def test_r2_het_wording_rule():
    for v in ("Kozz +/-", "-/+", "het", "HET_rep1", "Pigf_het1", "heterozygous", "Heterozygote", "het embryo"):
        assert bulk.HET_WORDING.search(v), v
    for v in ("wild type", "Kozz -/-", "heterogeneous tissue", "methyl", "Chet1", "whole embryo", "null"):
        assert not bulk.HET_WORDING.search(v), v


def test_r2_the_header_step_reads_4096_byte_blocks():
    data = b"gene\tA\n" + b"x" * 100_000
    seen = []

    def read(n):
        seen.append(n)
        return io.BytesIO(data).read(n) if len(seen) == 1 else b""
    assert bulk.first_line(read, False, "t") == b"gene\tA" and seen == [bulk.READ_BLOCK] and bulk.READ_BLOCK == 4096
    src = io.BytesIO(gzip.compress(data))
    got = []
    d = bulk.first_line(lambda n: (got.append(n), src.read(n))[1], True, "t")
    assert d == b"gene\tA" and set(got) == {4096}


def test_r2_single_cell_manifests_and_combine_json_are_byte_identical_with_and_without_bulk_units(tmp_path,
                                                                                                  monkeypatch):
    _, refp, ref = make_ref(tmp_path)
    rng = np.random.default_rng(1)
    sc_units = _units(ref, rng, ["S1", "S2", "S3"], mod="sc")
    b_units = _units(ref, rng, ["B1", "B2", "B3"])
    sc_rec = {"left_out": ["AS1:Gz"], "not_cleared": ["GSE1"], "inputs": {
        f"AS{i}": {"qc": f"q{i}", "qc_sha256": "s", "plan_sha256": "p"} for i in (1, 2, 3)}}
    both_rec = {"left_out": ["AB1:Gz", "AS1:Gz"], "not_cleared": ["GSE1"], "inputs": {
        **{f"AB{i}": {"qc": f"b{i}", "qc_sha256": "s", "plan_sha256": "p"} for i in (1, 2, 3)}, **sc_rec["inputs"]}}
    state = {}
    monkeypatch.setattr(co, "load_units", lambda *a, **k: (list(state["units"]), json.loads(json.dumps(state["rec"]))))
    monkeypatch.setattr(co.time, "strftime", lambda *a, **k: "2026-11-21 00:00:00 UTC")
    out, a = tmp_path / "tables", tmp_path / "sc_only"   # one output folder (the files name their own paths)
    state.update(units=sc_units, rec=sc_rec)
    co.combine(tmp_path, refp, out, dest=tmp_path)
    shutil.copytree(out, a)
    shutil.rmtree(out)
    state.update(units=sc_units + b_units, rec=both_rec)
    co.combine(tmp_path, refp, out, dest=tmp_path)
    b = out
    files = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file())
    assert any(str(f).endswith("_manifest.json") for f in files) and Path("COMBINE.json") in files
    for f in files:
        assert (a / f).read_bytes() == (b / f).read_bytes(), f
    bulk_rec = json.loads((b / "COMBINE_bulk.json").read_text())["record"]
    assert bulk_rec == {"left_out": ["AB1:Gz"], "not_cleared": [], "inputs": {
        f"AB{i}": {"qc": f"b{i}", "qc_sha256": "s", "plan_sha256": "p"} for i in (1, 2, 3)}}
    man = json.loads((b / "bulk_all" / "bulk_all_manifest.json").read_text())
    assert man["inputs"] == bulk_rec and not (a / "COMBINE_bulk.json").exists()
