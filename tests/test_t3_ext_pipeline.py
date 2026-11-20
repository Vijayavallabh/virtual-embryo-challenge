"""fetch.py's plan (no network: synthetic SOFT text), process.py's readers and a synthetic end-to-end run, combine.py's
weights and refusals, and the rule that no step writes inside the repository or feeds a T1/T2 builder."""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import sys
import tarfile
from pathlib import Path

import numpy as np
import pytest
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
    sample_file_patterns: ['(matrix\\.mtx|barcodes\\.tsv|features\\.tsv)(\\.gz)?$']
  - accession: GSE900002
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
  - accession: GSE900003
    status: held_phenocopy
    tier: wnt
    study: GSE900003
    design: conditional
    control_relation: same_capture
    repository: geo
    screen_series: [GSE900003]
    stages_used: [E8.5]
"""


def _soft(samples: dict[str, dict], title="A series", design="embryos", summary="a summary") -> str:
    lines = ["^SERIES = GSE900001", f"!Series_title = {title}", f"!Series_overall_design = {design}",
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


def _plan(dsy, dest, samples, acc="GSE900001", **kw):
    metas = {acc: fe.parse_soft(_soft(samples, **kw))}
    return fe.plan(acc, dest, dsy, metas=metas)


# ---------------------------------------------------------------------------------------------------------------------
# fetch.py plan


def test_plan_resolves_arms_stages_and_files(dsy, tmp_path):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT"),
         "GSM3": _sample("GSM3", "Hoxz KO", stage="E9.5"), "GSM4": _sample("GSM4", "WT", stage="E9.5"),
         "GSM5": _sample("GSM5", "Hoxz KO", stage="E7.5"),
         "GSM6": _sample("GSM6", "WT", extra=["assay: ATAC"], files=("atac_fragments.tsv.gz",))}
    pl = _plan(dsy, tmp_path / "d", s)
    assert pl["status"] == "ok", pl.get("reason")
    assert sorted((x["gsm"], x["arm"], x["stage"]) for x in pl["samples"]) == [
        ("GSM1", "ko", 8.5), ("GSM2", "control", 8.5), ("GSM3", "ko", 9.5), ("GSM4", "control", 9.5)]
    assert pl["units"] == [{"gene": "Hoxz", "ko": ["GSM1", "GSM3"], "control": ["GSM2", "GSM4"]}]
    for x in pl["samples"]:
        assert all(f["url"].startswith("https://ftp.ncbi.nlm.nih.gov/") for f in x["files"])
    assert (tmp_path / "d" / "meta" / "GSE900001" / "plan.json").is_file()


@pytest.mark.parametrize("bad", ["Ctnnb1 fl/fl", "Gata4 KO", "beta-catenin shRNA", "XAV939 treated"])
def test_plan_refuses_the_whole_dataset_on_any_held_out_sample(dsy, tmp_path, bad):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT"), "GSM9": _sample("GSM9", bad)}
    pl = _plan(dsy, tmp_path / "d", s)
    assert pl["status"] == "refused" and "held_out" in pl["reason"]
    assert "samples" not in pl


def test_plan_refuses_review_terms_and_held_out_series_titles(dsy, tmp_path):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT; Gata6 F/+")}
    assert _plan(dsy, tmp_path / "d", s)["status"] == "refused"
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT")}
    assert _plan(dsy, tmp_path / "d", s, title="Gata4 and Hoxz in the heart")["status"] == "refused"


@pytest.mark.parametrize("key", ["title", "design", "summary"])
@pytest.mark.parametrize("text,kind", [("Gata6 and Hoxz", "review"), ("Wnt signalling inhibition", "review"),
                                       ("GATA binding protein 4", "held_out"), ("sgbcat embryos", "held_out")])
def test_series_title_design_and_summary_refuse_alike(dsy, tmp_path, key, text, kind):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT")}
    pl = _plan(dsy, tmp_path / "d", s, **{key: text})
    assert pl["status"] == "refused" and kind in pl["reason"]


def test_a_series_review_term_is_cleared_only_by_a_named_series_entry(dsy, tmp_path):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT")}
    ok = DS_YAML.replace("    exclude_match: ['atac']\n", "    exclude_match: ['atac']\n    review_cleared:\n"
                         "      - {sample: GSE900001, term: gata6, reason: names a neighbour only}\n", 1)
    dsy.write_text(ok)
    assert _plan(dsy, tmp_path / "d", s, summary="Gata6 is mentioned")["status"] == "ok"
    assert _plan(dsy, tmp_path / "d", s, summary="Gata6 and Gata4 are mentioned")["status"] == "refused"
    dsy.write_text(ok.replace("sample: GSE900001", "sample: GSM1"))
    assert _plan(dsy, tmp_path / "d", s, summary="Gata6 is mentioned")["status"] == "refused"


def test_exclude_match_on_title_ignores_the_other_fields(dsy):
    ds = {"exclude_match": ["(?<![a-z])atac(?![a-z])"], "exclude_match_on": "title",
          "arms": {"control": {"on": "text", "all": ["genotype: wt"]}, "ko": []}}
    rna = {**_sample("E8.5_WT_rep1", "WT", extra=["library: multiome RNA (paired ATAC)"]), "title": ["E8.5_WT_rep1"]}
    atac = {**_sample("E8.5_WT_rep1_ATAC", "WT"), "title": ["E8.5_WT_rep1_ATAC"]}
    assert fe.resolve_arm(ds, rna)[0] == "control"
    assert fe.resolve_arm(ds, atac)[0] is None
    assert fe.resolve_arm({**ds, "exclude_match_on": "text"}, rna)[0] is None


@pytest.mark.parametrize("stage", ["E10.5", "E9.5-E10.5", "unknown"])
def test_plan_fails_closed_on_a_late_range_or_missing_stage(dsy, tmp_path, stage):
    s = {"GSM1": _sample("GSM1", "Hoxz KO", stage=stage), "GSM2": _sample("GSM2", "WT")}
    pl = _plan(dsy, tmp_path / "d", s)
    assert pl["status"] == "refused" and "stage" in pl["reason"]


def test_plan_pairs_by_chemistry_and_takes_a_genes_own_control(dsy, tmp_path):
    def smp(t, chem):
        r = _sample(t, "x")
        r["extract_protocol_ch1"] = f"Chromium Single Cell 3' {chem} chemistry"
        return r
    s = {"AAA_E8.5": smp("AAA_E8.5", "v3"), "BBB_E8.5": smp("BBB_E8.5", "v2"), "WT_E8.5_a": smp("WT_E8.5_a", "v3"),
         "WT_AAA_E8.5": smp("WT_AAA_E8.5", "v3")}
    pl = _plan(dsy, tmp_path / "d", s, acc="GSE900002")
    assert pl["status"] == "refused" and "no matched control" in pl["reason"]   # Bbb has no v2 WT
    s["WT_E8.5_b"] = smp("WT_E8.5_b", "v2")
    pl = _plan(dsy, tmp_path / "d", s, acc="GSE900002")
    assert pl["status"] == "ok", pl.get("reason")
    u = {x["gene"]: x for x in pl["units"]}
    assert u["Aaa"]["control"] == ["WT_AAA_E8.5"] and u["Bbb"]["control"] == ["WT_E8.5_b"]


def test_fetch_refuses_unlisted_held_or_unplanned_accessions_and_paths_inside_the_repo(dsy, tmp_path):
    with pytest.raises(cm.Refused):
        fe.plan("GSE999999", tmp_path / "d", dsy, metas={})
    with pytest.raises(cm.Refused, match="held_phenocopy"):
        fe.plan("GSE900003", tmp_path / "d", dsy, metas={})
    with pytest.raises(cm.Refused, match="inside the repository"):
        fe.plan("GSE900001", ROOT / "scratch_dl", dsy, metas={"GSE900001": fe.parse_soft(_soft({}))})
    with pytest.raises(cm.Refused, match="no plan"):
        fe.fetch("GSE900001", tmp_path / "empty", dsy)
    s = {"GSM1": _sample("GSM1", "Ctnnb1 KO"), "GSM2": _sample("GSM2", "WT")}
    _plan(dsy, tmp_path / "d", s)
    with pytest.raises(cm.Refused, match="did not pass"):
        fe.fetch("GSE900001", tmp_path / "d", dsy)
    for url in ("https://example.org/x.mtx.gz", "http://ftp.ncbi.nlm.nih.gov/x", "file:///etc/passwd"):
        with pytest.raises(cm.Refused):
            fe.https_url(url)
    assert fe.DENY.search("GSM1_atac_fragments.tsv.gz") and fe.DENY.search("x_R1.fastq.gz")
    assert fe.DENY.search("GSE1_RAW_raw.tar") and not fe.DENY.search("GSM1_matrix.mtx.gz")


def test_fetch_refuses_when_the_dataset_list_changed_after_the_plan(dsy, tmp_path):
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT")}
    assert _plan(dsy, tmp_path / "d", s)["status"] == "ok"
    dsy.write_text(DS_YAML + "\n# edited\n")
    with pytest.raises(cm.Refused, match="changed"):
        fe.fetch("GSE900001", tmp_path / "d", dsy)


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
# a synthetic end-to-end run and combine


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


def _synthetic_run(tmp_path, dsy):
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
    s = {"GSM1": _sample("GSM1", "Hoxz KO"), "GSM2": _sample("GSM2", "WT")}
    pl = _plan(dsy, dest, s)
    assert pl["status"] == "ok"
    fold = np.ones(len(PANEL))
    fold[:4] = 3.0
    man = {"accession": "GSE900001", "files": {}}
    for gsm, f in (("GSM1", fold), ("GSM2", np.ones(len(PANEL)))):
        d = dest / "raw" / "GSE900001" / gsm
        _write_sample(d, gsm, {"CM": 160, "Neural": 160}, mu, f, rng)
        for p in d.iterdir():
            man["files"][f"{gsm}/{p.name}"] = {"sha256": cm.sha256_file(p), "size": p.stat().st_size}
    man["plan_sha256"] = cm.sha256_file(dest / "meta" / "GSE900001" / "plan.json")
    cm.dump_json(man, dest / "raw" / "GSE900001" / "MANIFEST.json")

    def classifier(wt_path, C, present, genes, names, threads):
        calls = np.array([n.split(":")[1].split("_")[0] for n in names])
        return calls, np.full(len(names), 0.9), {"stub": True}

    return dest, wt, refp, classifier


def test_end_to_end_synthetic_run_writes_only_outside_and_combine_takes_it(tmp_path, dsy, capsys):
    before = sorted(p.name for p in EXT.iterdir())
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy)
    out = tmp_path / "processed"
    qc = pr.run("GSE900001", dest, wt, refp, out, dsy, classifier=clf)
    assert qc["status"] == "ok" and qc["panel_genes_present"] == 40
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


def test_confident_shares_are_recorded_and_a_table_under_the_floor_is_refused(tmp_path, dsy):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy)

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
    dsy.write_text(DS_YAML.replace("qc_floors:\n  min_confident_share_per_sample: 0.5\n", ""))
    with pytest.raises(cm.Refused, match="qc_floors"):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o2", dsy, classifier=clf)


def test_process_refuses_changed_raw_files_and_a_sample_that_fails_the_screen(tmp_path, dsy):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy)
    p = next((dest / "raw" / "GSE900001" / "GSM1").glob("*barcodes*"))
    p.write_bytes(gzip.compress(b"changed"))
    with pytest.raises(cm.Refused, match="sha256"):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o", dsy, classifier=clf)
    dest, wt, refp, clf = _synthetic_run(tmp_path / "b", dsy)
    pp = dest / "meta" / "GSE900001" / "plan.json"
    pl = json.loads(pp.read_text())
    pl["samples"][0]["record"]["description"] = ["littermate of a Gata4 KO"]
    pp.write_text(json.dumps(pl))
    man = json.loads((dest / "raw" / "GSE900001" / "MANIFEST.json").read_text())
    man["plan_sha256"] = cm.sha256_file(pp)
    (dest / "raw" / "GSE900001" / "MANIFEST.json").write_text(json.dumps(man))
    with pytest.raises(cm.Refused, match="before reading counts"):
        pr.run("GSE900001", dest, wt, refp, tmp_path / "o2", dsy, classifier=lambda *a: pytest.fail("read"))


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


def test_combine_refuses_a_failed_qc_or_a_changed_unit(tmp_path, dsy):
    dest, wt, refp, clf = _synthetic_run(tmp_path, dsy)
    out = tmp_path / "processed"
    pr.run("GSE900001", dest, wt, refp, out, dsy, classifier=clf)
    unit = next((out / "GSE900001").glob("*.unit.json"))
    unit.write_text(unit.read_text().replace('"gene"', '"gene" ', 1))
    with pytest.raises(cm.Refused, match="changed"):
        co.combine(out, refp, tmp_path / "t", dsy)
    qcp = out / "GSE900001" / "QC.json"
    qc = json.loads(qcp.read_text())
    qc["status"] = "failed"
    qcp.write_text(json.dumps(qc))
    with pytest.raises(cm.Refused, match="QC failed"):
        co.combine(out, refp, tmp_path / "t", dsy)


def test_the_launcher_plans_by_default_and_writes_all_exit_on_a_kill():
    import subprocess

    sh = EXT / "run_dgx.sh"
    assert subprocess.run(["bash", "-n", str(sh)], capture_output=True).returncode == 0
    text = sh.read_text()
    assert 'MODE=${2:-plan}' in text
    assert """trap 'echo 99 > "$D/logs/ALL.exit"; exit 99' TERM INT HUP""" in text
    plan_part, all_part = text.split('if [ "$MODE" = plan ]; then', 1)[1].split("\nfi\n", 1)
    assert "fetch.py plan" in plan_part and "fetch.py fetch" not in plan_part and "process.py" not in plan_part
    assert "fetch.py plan" not in all_part and 'plan_$A.exit' in all_part


def test_outside_repo_guard(tmp_path):
    with pytest.raises(cm.Refused):
        cm.outside_repo(ROOT / "x" / "y")
    with pytest.raises(cm.Refused):
        cm.outside_repo(ROOT)
    with pytest.raises(cm.Refused):
        tb.write_table(ROOT / "t", "n", ["a"], np.ones(1), {}, {})
    assert cm.outside_repo(tmp_path / "x") == (tmp_path / "x").resolve()


def test_no_t1_or_t2_builder_reads_the_outside_knockout_code():
    files = [*(ROOT / "src" / "vec").glob("*.py"), *(ROOT / "experiments" / "builds").glob("*.py"),
             *(ROOT / "experiments" / "builds").glob("*.sh"), *(ROOT / "scripts").glob("*.py"),
             *(ROOT / "scripts").glob("*.sh")]
    for d in sorted((ROOT / "experiments").glob("t[12]_*")):
        files += [*d.glob("*.py"), *d.glob("*.sh")]
    assert files
    for f in files:
        assert "t3_ext" not in f.read_text(errors="replace"), f
