"""Stage 9, phase 2: the many-line ArrayExpress resource (experiments/t3_ext/ae.py) on a SYNTHETIC resource only:
several lines, technical-replicate runs (paired rows), a withheld line whose rows hold a held-out term, a duplicate
individual, a disagreeing record, a day-10 embryo, a het, the baseline line and an unannotated run column, walked from
plan to tables. No network: fetch._open and fetch._get_ae are stubs."""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
T3 = ROOT / "experiments" / "t3_decomp"
for p in (EXT, T3):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import ae  # noqa: E402
import bulk  # noqa: E402
import combine as co  # noqa: E402
import fetch as fe  # noqa: E402
import make_panel_ids as mpi  # noqa: E402
import process as pr  # noqa: E402
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402

ACC = "E-MTAB-7322"
LINES = ("Abcz1", "Defz2", "Ghiz3")
PANEL = [f"g{i}" for i in range(38)] + ["Xist", "Tsix"]
GROUPS = ["CM", "Endoderm", "Neural"]
FILLERS = [f"f{i}" for i in range(100)]
ID = {s: f"ENSMUSG{200000 + i:011d}" for i, s in enumerate(PANEL + FILLERS + list(cm.SEX_GENES[1:]) + list(LINES))}
PAT = r"(?P<gene>[^_]+)_(?P<genotype>wt|het|hom)\d+"   # the README example: <line>_<wt|het|hom><n>
WITHHELD = ["Trps1", "Anks6", "Atxn10", "B9d2", "Cc2d2a", "Ift140", "Kif3b", "Kifap3", "Rpgrip1l", "Mad2l2", "Tmem98",
            "Nek9", "Mtor", "Map2k7", "Traf6", "Traf2", "Rala", "Pigf", "Pigl", "Pgap2", "Hira", "Smc3",
            "1700007K13Rik", "4933434E20Rik", "Hmgxb3", "Zkscan14", "Zkscan17", "Dhx35", "Gpr107"]
YAML = f"""version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
  - accession: {ACC}
    status: cleared
    tier: unrelated
    study: MGP-DMDD
    design: germline
    control_relation: littermate
    repository: arrayexpress
    modality: bulk
    count_form: arrayexpress
    id_type: ensembl
    gene_columns: ["Gene ID"]
    source_name_pattern: '{PAT}'
    withheld_lines:
      - {{line: Trps1, reason: held}}
""" + "".join(f"      - {{line: {g}, reason: pending_phenocopy}}\n" for g in WITHHELD[1:]) + f"""\
    screen_series: [{ACC}]
    stages_used: [E9.5]
    pinned: []  # {ACC}
"""
MARK = f"    pinned: []  # {ACC}\n"
H = ["Source Name", "Comment[ENA_SAMPLE]", "Characteristics[organism]", "Characteristics[developmental stage]",
     "Characteristics[mutant line]", "Comment[allele]", "Characteristics[somite stage]", "Characteristics[sex]",
     "Characteristics[genotype]", "Characteristics[individual]", "Protocol REF", "Assay Name", "Comment[ENA_RUN]",
     "Comment[FASTQ_URI]", "Factor Value[somite stage]"]
GENO = {"hom": "homozygous knockout", "het": "heterozygous knockout", "wt": "wild type genotype"}


def embryos():
    """[(Source Name, line, token, stage, n runs, individual, overrides)]: the synthetic resource's embryos."""
    out, ind = [], 100
    plan = {"Abcz1": {"hom": 3, "wt": 3, "het": 1}, "Defz2": {"hom": 2, "wt": 3}, "Ghiz3": {"hom": 1, "wt": 3},
            "Trps1": {"hom": 2, "wt": 2}}
    runs = [1, 2, 4]
    for line, arms in plan.items():
        for tok, n in arms.items():
            for k in range(1, n + 1):
                ind += 1
                sn = f"{line}_{tok}{k}"   # the individual is the Source Name, as in the real resource
                out.append((sn, line, tok, "embryonic day 9.5", runs[(ind + k) % 3], sn, {}))
    ind += 1
    out.append(("Abcz1_hom9", "Abcz1", "hom", "embryonic day 10", 1, "Abcz1_hom9", {}))        # past the cap
    out.append(("Abcz1_hom8", "Abcz1", "hom", "embryonic day 8.5", 1, "Abcz1_hom8", {}))   # not E9.5
    out.append(("Abcz1_wt7", "Abcz1", "wt", "embryonic day 9.5", 2, "Abcz1_wt7", {"sex_row2": "male"}))   # (a)
    out.append(("Defz2_wt9", "Defz2", "wt", "embryonic day 9.5", 1, "Defz2_wt9", {"genotype": GENO["hom"]}))  # (b)
    out.append(("Defz2_wt8", "Defz2", "wt", "embryonic day 9.5", 1, "dup-1", {}))         # (e) duplicate individual
    out.append(("Ghiz3_wt8", "Ghiz3", "wt", "embryonic day 9.5", 1, "dup-1", {}))
    out.append(("baseline_wt1", "baseline", "wt", "embryonic day 9.5", 1, "baseline_wt1", {"allele": ""}))
    out.append(("baseline_wt2", "baseline", "wt", "embryonic day 9.5", 1, "baseline_wt2", {"allele": ""}))
    return out


def sdrf_and_runs():
    rows, runs_of, ers_of, err = ["\t".join(H)], {}, {}, 1000000
    for i, (sn, line, tok, stage, n, ind, ov) in enumerate(embryos()):
        ers = f"ERS{3000000 + i}"
        ers_of[sn] = ers
        runs_of[sn] = []
        allele = ov.get("allele", f"{line}<tm1a(EUCOMM)Wtsi>" if line != "Trps1" else "Trps1<Gata4-reporter>")
        for r in range(n):
            err += 1
            run = f"ERR{err}"
            runs_of[sn].append(run)
            for read in (1, 2):
                sex = ov.get("sex_row2", "female") if (r == 0 and read == 2) else "female"
                rows.append("\t".join([sn, ers, "Mus musculus", stage, line, allele, "12-somite stage", sex,
                                       ov.get("genotype", GENO[tok]), ind, "P-MTAB-1", f"{sn}_{r}", run,
                                       f"ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/{run}/{run}_{read}.fastq.gz",
                                       "12-somite stage"]))
    return ("\n".join(rows) + "\n").encode(), runs_of, ers_of


IDF = (b"Investigation Title\tWhole-embryo RNA-seq of knockout lines at E9.5\n"
       b"Experiment Description\tWhole embryos of knockout lines and wild-type littermates, one embryo a library\n"
       b"Protocol Description\tembryos dissected free of yolk sac\tpaired-end sequencing\n"
       b"Comment[ArrayExpressAccession]\tE-MTAB-7322\n")
METHODS = b"Analysis\tMethod\nRead filtering\tadapters trimmed\nMapping\tEnsembl release 98\nQuantification\tper run\n"


def texts(sdrf=None, idf=IDF):
    return {"idf": idf, "sdrf": sdrf if sdrf is not None else sdrf_and_runs()[0], "methods": METHODS}


@pytest.fixture()
def y(tmp_path):
    p = tmp_path / "datasets.yaml"
    p.write_text(YAML)
    return p


def make_ref(tmp_path):
    rng = np.random.default_rng(7)
    mu = {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}
    wt = tmp_path / "wt.h5ad"
    wt.write_bytes(b"stand-in WT")
    refp = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats({"genes": PANEL, "groups": GROUPS, "shares": {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5},
                       "mu": mu, "count": {"CM": 2, "Endoderm": 3, "Neural": 5}}, refp, cm.sha256_file(wt))
    return wt, refp, tb.load_ref_stats(refp)


def matrix(ref, runs_of, rng) -> bytes:
    """The count matrix: every annotated run (withheld lines' runs included, as the real matrix holds them) and one
    unannotated run column holding text (it must never be parsed)."""
    genes = PANEL + FILLERS + list(cm.SEX_GENES[1:]) + list(LINES)
    cols = {}
    for sn, runs in runs_of.items():
        line, tok = sn.split("_")[0], sn.split("_")[1].rstrip("0123456789")
        p = tb.whole_den(ref).copy()
        if tok == "hom" and line in LINES:
            k = LINES.index(line) * 4
            p[k:k + 4] *= 2.0
        xp = rng.multinomial(300_000, p / p.sum()).astype(float)
        lib = dict(zip(PANEL, xp, strict=True))
        lib["Xist"] = 3000.0
        lib.update({f: float(rng.poisson(17_000)) for f in FILLERS})
        lib.update({g: 0.0 for g in cm.SEX_GENES[1:]})
        lib.update({g: (100.0 if (g == line and tok == "hom") else 4000.0) for g in LINES})
        per_run = {r: {} for r in runs}
        for g in genes:
            parts = rng.multinomial(int(lib[g]), np.full(len(runs), 1 / len(runs)))
            for r, v in zip(runs, parts, strict=True):
                per_run[r][g] = int(v)
        cols.update(per_run)
    head = ["Gene ID", *cols, "ERR9999999"]
    lines = ["\t".join(head)] + ["\t".join([ID[g], *[str(cols[r][g]) for r in cols], "x"]) for g in genes]
    return ("\n".join(lines) + "\n").encode()


def opener(data: bytes, calls: list):
    import io

    def op(url):
        calls.append(url)
        return io.BytesIO(data)
    return op


def stage(tmp_path, y, monkeypatch):
    wt, refp, ref = make_ref(tmp_path)
    dest = tmp_path / "dest"
    sdrf, runs_of, ers_of = sdrf_and_runs()
    pl = fe.plan(ACC, dest, y, metas=texts(sdrf))
    assert pl["status"] == "ok", pl.get("reason")
    data = matrix(ref, runs_of, np.random.default_rng(3))
    calls: list = []
    monkeypatch.setattr(fe, "_open", opener(data, calls))
    rec = bulk.header(ACC, dest, datasets=y)
    assert rec["status"] == "ok", rec.get("reason")
    assert calls == ["https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7322/"
                     "E-MTAB-7322-raw-counts.tsv.undecorated"]
    assert pl["line_set_basis"]["final"] is False
    pl = fe.plan(ACC, dest, y, metas=texts(sdrf))   # plan again: the line set is final after the header step
    assert pl["status"] == "ok" and pl["line_set_basis"]["final"] is True, pl.get("reason")
    return dest, wt, refp, ref, pl, rec, data, runs_of, ers_of, sdrf


def pin_all(dest, y, pl):
    pins = "".join(ae.pin_template(ACC, dest, arm, line, v[arm], y) for line, v in pl["line_set"].items()
                   for arm in ("ko", "control"))
    block = ("    pinned_series_files:\n" + bulk.table_template(ACC, dest, y) + "    pinned_series:\n"
             + ae.series_pin_template(ACC, dest, y) + "    pinned:\n" + pins)
    y.write_text(y.read_text().replace(MARK, block))


def test_the_urls_are_derived_from_the_accession_on_one_host():
    u = ae.urls(ACC)
    assert u == {"idf": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/322/E-MTAB-7322/Files/E-MTAB-7322.idf.txt",
                 "sdrf": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/322/E-MTAB-7322/Files/E-MTAB-7322.sdrf.txt",
                 "methods": "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7322/"
                            "E-MTAB-7322-analysis-methods.tsv",
                 "matrix": "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-7322/"
                           "E-MTAB-7322-raw-counts.tsv.undecorated"}
    for bad in ("http://ftp.ebi.ac.uk/x", "https://ftp.ebi.ac.uk:21/x", "https://ftp.ebi.ac.uk/x?y=1",
                "https://www.ebi.ac.uk/x", "https://ftp.ncbi.nlm.nih.gov/x"):
        with pytest.raises(cm.Refused):
            fe.ae_url(bad)
    with pytest.raises(cm.Refused):
        ae.urls("E-MTAB-7325")
    with pytest.raises(cm.Refused):
        ae.urls("E-MTAB-5224")


def test_the_many_line_resource_from_plan_to_tables(tmp_path, y, monkeypatch):
    dest, wt, refp, ref, pl, rec, data, runs_of, ers_of, sdrf = stage(tmp_path, y, monkeypatch)
    # the withheld lines: dropped before the screen (a Trps1 row holds a held-out term), the views, the hashes, tables
    n_trps1 = sum(2 * e[4] for e in embryos() if e[1] == "Trps1")   # two rows (read 1, read 2) per run
    assert pl["withheld_rows"]["Trps1"] == {"rows": n_trps1, "reason": "held"}
    assert pl["withheld_rows"]["B9d2"] == {"rows": 0, "reason": "pending_phenocopy"}
    table = (dest / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert table.count("Trps1") == 1 and f"# Trps1: {n_trps1} rows withheld (held)" in table
    assert "Gata4" not in table and "Gata4" not in json.dumps(pl["tables"]) + json.dumps(pl["views"])
    assert all(r["line"] != "Trps1" for r in pl["tables"][ACC])
    # the comparisons
    c = pl["comparisons"]
    assert c["a"] == ["Abcz1_wt7"] and c["b"] == ["Defz2_wt9"] and c["e"] == ["Defz2_wt8", "Ghiz3_wt8"]
    assert c["c"] == [] and c["d"] == [] and c["f"] == []
    rows = {r["source_name"]: r for r in pl["tables"][ACC]}
    assert "after E9.5" in rows["Abcz1_hom9"]["not_pinnable"]
    assert "only 'embryonic day 9.5'" in rows["Abcz1_hom8"]["not_pinnable"]
    assert "heterozygous" in rows["Abcz1_het1"]["not_pinnable"]
    assert "baseline" in rows["baseline_wt1"]["not_pinnable"]
    # the record hash is the canonical serialisation of the Source Name's fields
    r = rows["Abcz1_hom1"]
    want = hashlib.sha256(json.dumps({"accession": "Abcz1_hom1", "fields": r["fields"]}, sort_keys=True,
                                     ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
    assert r["metadata_sha256"] == want and list(r["fields"])[0] == "001|Source Name"
    assert len(r["runs"]) == len(runs_of["Abcz1_hom1"]) and len(r["fields"]["013|Comment[ENA_RUN]"]) == 2 * len(
        r["runs"])
    # the committed rule fixes the line set
    assert pl["line_set"] == {"Abcz1": {"ko": ["Abcz1_hom1", "Abcz1_hom2", "Abcz1_hom3"],
                                        "control": ["Abcz1_wt1", "Abcz1_wt2", "Abcz1_wt3"]},
                              "Defz2": {"ko": ["Defz2_hom1", "Defz2_hom2"],
                                        "control": ["Defz2_wt1", "Defz2_wt2", "Defz2_wt3"]}}
    # the header: the withheld line's runs and the unannotated run are tied to nothing and never read
    trps_runs = sum(len(v) for k, v in runs_of.items() if k.startswith("Trps1"))
    assert rec["ties"]["untied"] == trps_runs + 1 and rec["ties"]["incomplete"] == []
    assert "ERR9999999" not in "\n".join(ae.header_lines(rec))
    # pins: exactly the rule's eligible Source Names, by line and arm
    with pytest.raises(cm.Refused, match="exactly the rule's eligible"):
        ae.pin_template(ACC, dest, "ko", "Abcz1", ["Abcz1_hom1", "Abcz1_hom2"], y)
    with pytest.raises(cm.Refused, match="not a line of the committed rule"):
        ae.pin_template(ACC, dest, "ko", "Trps1", ["Trps1_hom1", "Trps1_hom2"], y)
    with pytest.raises(cm.Refused, match="not a line of the committed rule"):
        ae.pin_template(ACC, dest, "ko", "Ghiz3", ["Ghiz3_hom1"], y)
    pin_all(dest, y, pl)
    ds = cm.dataset(ACC, y)
    assert {p["accession"] for p in ds["pinned"]} == {ers_of[s] for v in pl["line_set"].values() for a in v.values()
                                                      for s in a}
    pl2 = fe.plan(ACC, dest, y, metas=texts(sdrf))
    assert pl2["status"] == "ok", pl2.get("reason")
    # fetch: the one matrix, from its derived url, its first line checked
    staging = tmp_path / "staging"
    staging.mkdir()
    (staging / ae.matrix_name(ACC)).write_bytes(data)
    calls: list = []

    def get(url, dst, max_bytes):
        calls.append(url)
        dst = cm.outside_repo(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(staging / Path(url).name, dst)
        return dst.stat().st_size, cm.sha256_file(dst)

    monkeypatch.setattr(fe, "_get_ae", get)
    monkeypatch.setattr(fe, "_get", lambda *a, **k: (_ for _ in ()).throw(AssertionError("a GEO download")))
    with pytest.raises(cm.Refused, match="PANEL_IDS_SHA256 is None"):
        fe.fetch(ACC, dest, y)
    rows_ids = [(ID[s], s) for s in ID]
    text = mpi.build(rows_ids, ref["genes"], list(LINES), "GSE0/GSM0/GSM0_features.tsv.gz", "0" * 64)
    pf = tmp_path / "panel_ids.tsv"
    pf.write_text(text)
    monkeypatch.setattr(cm, "PANEL_IDS", pf)
    monkeypatch.setattr(cm, "PANEL_IDS_SHA256", cm.sha256_file(pf))
    man = fe.fetch(ACC, dest, y)
    assert calls == [ae.urls(ACC)["matrix"]] and list(man["files"]) == [f"_series/{ae.matrix_name(ACC)}"]
    # process: one library per Source Name, its runs summed; the text column is never parsed
    out = tmp_path / "processed"
    qc = pr.run(ACC, dest, wt, refp, out, y)
    assert qc["status"] == "ok", qc.get("reason")
    n_runs = sum(len(runs_of[s]) for v in pl["line_set"].values() for a in v.values() for s in a)
    assert qc["read"]["columns_read"] == n_runs
    assert set(qc["units"]) == {"Abcz1", "Defz2"} and all(u["status"] == "ok" for u in qc["units"].values()), {
        g: u.get("reason") for g, u in qc["units"].items()}
    lab = ers_of["Abcz1_hom1"]
    assert qc["libraries"][lab]["column"] == "+".join(runs_of["Abcz1_hom1"])
    for g in ("Abcz1", "Defz2"):
        L = json.loads(Path(qc["units"][g]["unit_file"]).read_text())["L"]["__whole__"]
        assert L[PANEL.index("Xist")] is None and L[PANEL.index("Tsix")] is None
    r = co.combine(out, refp, tmp_path / "tables", y, dest=dest)
    assert r["tables"] == ["bulk_all", "bulk_study_MGP-DMDD"]
    man = json.loads((tmp_path / "tables" / "bulk_all" / "bulk_all_manifest.json").read_text())
    assert man["line_agreement_by_study"]["MGP-DMDD"] and man["units_by_study"]["MGP-DMDD"] == [
        f"{ACC}:Abcz1", f"{ACC}:Defz2"]


def test_pins_that_are_not_the_rule_s_set_are_refused(tmp_path, y, monkeypatch):
    dest, wt, refp, ref, pl, rec, data, runs_of, ers_of, sdrf = stage(tmp_path, y, monkeypatch)
    pin_all(dest, y, pl)
    good = y.read_text()
    ds = cm.dataset(ACC, y)
    pl2 = fe.plan(ACC, dest, y, metas=texts(sdrf))
    # a hand-picked subset (one control dropped)
    drop = [p for p in ds["pinned"] if p["source_name"] != "Abcz1_wt3"]
    with pytest.raises(cm.Refused, match="not exactly the rule's eligible"):
        ae.verify_pins({**ds, "pinned": drop}, pl2)
    # a line dropped
    with pytest.raises(cm.Refused, match="line set"):
        ae.verify_pins({**ds, "pinned": [p for p in ds["pinned"] if p["gene"] != "Defz2"]}, pl2)
    # runs other than the record's; another arm; a stored record edited (comparison (g))
    p0 = ds["pinned"][0]
    with pytest.raises(cm.Refused, match="runs"):
        ae.verify_pins({**ds, "pinned": [{**p0, "runs": p0["runs"][:1] + ["ERR1999999"]}, *ds["pinned"][1:]]}, pl2)
    with pytest.raises(cm.Refused, match="genotype"):
        ae.verify_pins({**ds, "pinned": [{**p0, "arm": "control"}, *ds["pinned"][1:]]}, pl2)
    pl3 = json.loads(json.dumps(pl2))
    pl3["tables"][ACC][0]["fields"]["008|Characteristics[sex]"][0] = "male"
    with pytest.raises(cm.Refused, match="comparison"):
        ae.verify_pins(ds, pl3)
    # a plan made under another withheld_lines declaration
    y.write_text(good.replace("      - {line: Atxn10, reason: pending_phenocopy}\n", ""))
    with pytest.raises(cm.Refused, match="another declaration"):
        ae.verify_pins(cm.dataset(ACC, y), pl2)
    # the header record changed
    y.write_text(good)
    pl4 = json.loads(json.dumps(pl2))
    pl4["count_table_headers"][ae.matrix_name(ACC)]["ties"]["untied"] += 1
    with pytest.raises(cm.Refused, match="derived again"):
        ae.verify_pins(ds, pl4)


def test_grammar_and_text_refusals_withhold_the_text(tmp_path, y):
    sdrf, *_ = sdrf_and_runs()
    for bad in (sdrf.replace(b"Characteristics[sex]", b"Sex of embryo", 1),
                sdrf.replace(b"\n", b"\n\xff", 1),
                sdrf.replace(b"Mus musculus", b"Mus\x0cmusculus", 1),
                sdrf.replace(b"\tP-MTAB-1", b"", 1)):
        pl = fe.plan(ACC, tmp_path / "d", y, metas=texts(bad))
        assert pl["status"] == "refused" and pl["withheld"] is True
        t = (tmp_path / "d" / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
        assert "Abcz1" not in t and "WITHHELD" in t
    pl = fe.plan(ACC, tmp_path / "d", y, metas=texts(idf=b"Bad|Key\tx\n"))
    assert pl["status"] == "refused" and pl["withheld"] is True


def test_a_near_match_of_a_withheld_line_refuses_the_resource(tmp_path, y):
    sdrf, *_ = sdrf_and_runs()
    lines = sdrf.decode().split("\n")
    lines.append(lines[1].replace("Abcz1_hom1", "trps1_hom1").replace("\tAbcz1\t", "\tTRPS1\t"))
    pl = fe.plan(ACC, tmp_path / "d", y, metas=texts("\n".join(lines).encode()))
    assert pl["status"] == "refused" and "near the withheld" in pl["reason"]


def test_a_held_out_term_outside_the_withheld_rows_withholds_the_resource(tmp_path, y):
    sdrf, *_ = sdrf_and_runs()
    bad = sdrf.replace(b"Abcz1<tm1a(EUCOMM)Wtsi>", b"Abcz1<tm1a(EUCOMM)Wtsi>;Ctnnb1", 1)
    pl = fe.plan(ACC, tmp_path / "d", y, metas=texts(bad))
    assert pl["status"] == "refused" and pl["withheld"] is True
    t = (tmp_path / "d" / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert "Ctnnb1" not in t and "Abcz1" not in t
    pl = fe.plan(ACC, tmp_path / "d", y, metas=texts(idf=IDF + b"Comment[Related]\tE-MTAB-7325\n"))
    assert pl["status"] == "refused" and pl["withheld"] is True


def test_the_arrayexpress_keys_and_pins(tmp_path):
    p = tmp_path / "y.yaml"
    for old, new, match in (("reason: held}", "reason: maybe}", "withheld_lines"),
                            (f"source_name_pattern: '{PAT}'", "source_name_pattern: '(.*)'", "named groups"),
                            ("repository: arrayexpress", "repository: geo", "repository arrayexpress"),
                            ("    withheld_lines:\n", "    withheld_lines:\n      - {line: Trps1, reason: held}\n",
                             "withheld_lines")):
        p.write_text(YAML.replace(old, new, 1))
        with pytest.raises(cm.Refused, match=match):
            cm.load_datasets(p)
    pin = {"accession": "ERS1", "source_name": "a", "arm": "ko", "gene": "Abcz1", "stage": "E9.5",
           "runs": ["ERR1000001"], "metadata_sha256": "0" * 64}
    ds = cm.load_datasets(_y(tmp_path, YAML))[ACC]
    bad = [pin, {**pin, "accession": "ERS2", "runs": ["ERR1000002"]},
           {**pin, "accession": "ERS3", "arm": "control", "source_name": "b", "runs": ["ERR1000003"]},
           {**pin, "accession": "ERS4", "arm": "control", "source_name": "c", "runs": ["ERR1000004"]}]
    with pytest.raises(cm.Refused, match="two accessions"):
        cm.check_bulk_pins(ACC, {**ds, "pinned": bad})
    bad[1]["source_name"] = "d"
    bad[1]["runs"] = ["ERR1000001"]
    with pytest.raises(cm.Refused, match="run ERR1000001"):
        cm.check_bulk_pins(ACC, {**ds, "pinned": bad})
    with pytest.raises(cm.Refused, match="ERS or SAMEA"):
        cm.check_bulk_pins(ACC, {**ds, "pinned": [{**pin, "accession": "GSM1"}]})
    for acc in ("E-MTAB-7325",):
        with pytest.raises(cm.Refused, match="NEVER_FETCH"):
            fe.plan(acc, tmp_path / "d", None, metas=texts())


def test_the_run_ties_refuse_a_doubled_header_a_run_of_two_records_and_mark_partial_records(tmp_path, y):
    sdrf, runs_of, ers_of = sdrf_and_runs()
    pl = fe.plan(ACC, tmp_path / "d", y, metas=texts(sdrf))
    rows = pl["tables"][ACC]
    run = runs_of["Abcz1_hom1"][0]
    ties, problems = ae.run_ties(["Gene ID", run, run], rows, ["Gene ID"])
    assert any("twice" in p for p in problems)
    ties, problems = ae.run_ties(["ID", run], rows, ["Gene ID"])
    assert any("gene column" in p for p in problems)
    two = [dict(rows[0]), {**rows[1], "runs": [run]}]
    assert any("ties to 2" in p for p in ae.run_ties(["Gene ID", run], two, ["Gene ID"])[1])
    multi = [r for r in rows if len(r["runs"]) > 1][0]
    ties, problems = ae.run_ties(["Gene ID", multi["runs"][0]], rows, ["Gene ID"])
    assert problems == [] and multi["source_name"] in ties["incomplete"]
    ties, _ = ae.run_ties(["Gene ID", "Abcz1_hom1", *runs_of["Abcz1_hom1"]], rows, ["Gene ID"])
    assert "Abcz1_hom1" in ties["other_tie"]   # tied through a non-run header too: not pinnable (double count)
    dest = tmp_path / "d"
    rec = ae.record_header(cm.dataset(ACC, y), dest, pl, ae.matrix_name(ACC), ae.urls(ACC)["matrix"],
                           b"Gene ID\tCtnnb1_ERR1\t" + run.encode())
    assert rec["status"] == "refused" and rec["withheld"] and "Ctnnb1" not in json.dumps(rec)


def _y(tmp_path, text):
    p = tmp_path / "y2.yaml"
    p.write_text(text)
    return p


# ---------------------------------------------------------------------------------------------------------------------
# stage 9 round 2: each blocker of the two reviews of 120af2f, replayed


def drop_columns(data: bytes, drop: set[str]) -> bytes:
    lines = data.decode().rstrip("\n").split("\n")
    keep = [i for i, h in enumerate(lines[0].split("\t")) if h not in drop]
    return ("\n".join("\t".join(ln.split("\t")[i] for i in keep) for ln in lines) + "\n").encode()


def fetch_process_combine(tmp_path, dest, y, ref, wt, refp, data, monkeypatch):
    staging = tmp_path / "staging"
    staging.mkdir()
    (staging / ae.matrix_name(ACC)).write_bytes(data)

    def get(url, dst, max_bytes):
        dst = cm.outside_repo(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(staging / Path(url).name, dst)
        return dst.stat().st_size, cm.sha256_file(dst)

    monkeypatch.setattr(fe, "_get_ae", get)
    text = mpi.build([(ID[s], s) for s in ID], ref["genes"], list(LINES), "GSE0/GSM0/GSM0_features.tsv.gz", "0" * 64)
    pf = tmp_path / "panel_ids.tsv"
    pf.write_text(text)
    monkeypatch.setattr(cm, "PANEL_IDS", pf)
    monkeypatch.setattr(cm, "PANEL_IDS_SHA256", cm.sha256_file(pf))
    fe.fetch(ACC, dest, y)
    qc = pr.run(ACC, dest, wt, refp, tmp_path / "processed", y)
    r = co.combine(tmp_path / "processed", refp, tmp_path / "tables", y, dest=dest)
    return qc, r


def with_reused_run(sdrf: bytes, runs_of: dict) -> tuple[bytes, str]:
    """The first reviewer's probe: a kept Abcz1 wild-type row whose run is a withheld (Trps1) row's run."""
    trun = runs_of["Trps1_hom1"][0]
    lines = sdrf.decode().rstrip("\n").split("\n")
    f = next(ln for ln in lines[1:] if ln.startswith("Abcz1_wt1\t")).split("\t")
    f[0], f[1], f[9], f[11], f[12] = "Abcz1_wt9", "ERS3999999", "Abcz1_wt9", "Abcz1_wt9_0", trun
    f[13] = f"ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/{trun}/{trun}_1.fastq.gz"
    lines.append("\t".join(f))
    return ("\n".join(lines) + "\n").encode(), trun


def test_r2_blocker1_a_kept_row_reusing_a_withheld_run_is_withheld_and_the_run_is_never_pinned(tmp_path, y,
                                                                                                monkeypatch):
    sdrf, runs_of, ers_of = sdrf_and_runs()
    bad, trun = with_reused_run(sdrf, runs_of)
    pl = fe.plan(ACC, tmp_path / "d", y, metas=texts(bad))
    assert pl["status"] == "ok", pl.get("reason")
    rows = {r["source_name"]: r for r in pl["tables"][ACC]}
    # the reviewer's result was: not_pinnable None, in line_set['Abcz1']['control'], run_ties 'complete'
    assert "Abcz1_wt9" not in rows and "Abcz1_wt9" not in pl["line_set"]["Abcz1"]["control"]
    assert pl["withheld_linked_rows"] == 1 and ae.run_digest(trun) in pl["withheld_run_sha256"]
    text = (tmp_path / "d" / "meta" / ACC / "plan.json").read_text() + (
        tmp_path / "d" / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert trun not in text and "Abcz1_wt9" not in text and "Trps1_hom" not in text
    # an empty value links nothing (the second reviewer's lesser point 5)
    kept, _, linked, _ = ae.withhold_rows(cm.dataset(ACC, y), H, [
        ["Trps1_hom1", "", "x", "s", "Trps1", "a", "1", "f", "g", "", "p", "a", "ERR1", "u", "1"],
        ["Abcz1_wt1", "", "x", "s", "Abcz1", "a", "1", "f", "g", "", "p", "a", "ERR2", "u", "1"]])
    assert linked == 0 and len(kept) == 1
    # the pin check and the pin template refuse a withheld run
    dest, wt, refp, ref, pl, rec, data, runs_of, ers_of, sdrf = stage(tmp_path, y, monkeypatch)
    pin_all(dest, y, pl)
    ds = cm.dataset(ACC, y)
    pl2 = fe.plan(ACC, dest, y, metas=texts(sdrf))
    assert pl2["status"] == "ok", pl2.get("reason")
    trun = runs_of["Trps1_hom1"][0]
    p0 = ds["pinned"][0]
    with pytest.raises(cm.Refused, match="withheld line"):
        ae.verify_pins({**ds, "pinned": [{**p0, "runs": [*p0["runs"], trun]}, *ds["pinned"][1:]]}, pl2)
    pp = dest / "meta" / ACC / "plan.json"
    stored = json.loads(pp.read_text())
    for r in stored["tables"][ACC]:
        if r["source_name"] == "Abcz1_hom1":
            r["runs"].append(trun)
    pp.write_text(json.dumps(stored))
    with pytest.raises(cm.Refused, match="withheld line"):
        ae.pin_template(ACC, dest, "ko", "Abcz1", pl2["line_set"]["Abcz1"]["ko"], y)
    with pytest.raises(cm.Refused, match="withheld line"):
        ae.verify_pins(ds, stored)


def test_r2_blocker6_one_eligible_run_missing_from_the_matrix_the_rest_pins_fetches_processes_and_combines(
        tmp_path, y, monkeypatch):
    wt, refp, ref = make_ref(tmp_path)
    dest = tmp_path / "dest"
    sdrf, runs_of, ers_of = sdrf_and_runs()
    pl = fe.plan(ACC, dest, y, metas=texts(sdrf))
    table = (dest / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert pl["line_set_basis"]["final"] is False and "THE LINE SET IS PROVISIONAL" in table
    assert "Defz2_wt1" in pl["line_set"]["Defz2"]["control"]
    gone = runs_of["Defz2_wt1"][0]
    data = drop_columns(matrix(ref, runs_of, np.random.default_rng(3)), {gone})
    monkeypatch.setattr(fe, "_open", opener(data, []))
    rec = bulk.header(ACC, dest, datasets=y)
    assert rec["status"] == "ok" and rec["ties"]["incomplete"] == ["Defz2_wt1"]
    # the template refuses until the plan carries the passing header record (the line set is provisional)
    with pytest.raises(cm.Refused, match="provisional"):
        ae.pin_template(ACC, dest, "control", "Defz2", ["Defz2_wt2", "Defz2_wt3"], y)
    pl = fe.plan(ACC, dest, y, metas=texts(sdrf))
    assert pl["status"] == "ok" and pl["line_set_basis"]["final"] is True
    rows = {r["source_name"]: r for r in pl["tables"][ACC]}
    assert "(h) a run missing from the matrix header" in rows["Defz2_wt1"]["not_pinnable"]
    assert pl["comparisons"]["h"] == ["Defz2_wt1"]
    assert pl["line_set"]["Defz2"] == {"ko": ["Defz2_hom1", "Defz2_hom2"], "control": ["Defz2_wt2", "Defz2_wt3"]}
    table = (dest / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert "THE LINE SET IS FINAL" in table and "(h) 1 disagreements" in table
    # the reviewer's probe: the old set is refused, the remaining set pins
    with pytest.raises(cm.Refused, match="exactly the rule's eligible"):
        ae.pin_template(ACC, dest, "control", "Defz2", ["Defz2_wt1", "Defz2_wt2", "Defz2_wt3"], y)
    pin_all(dest, y, pl)
    pl2 = fe.plan(ACC, dest, y, metas=texts(sdrf))
    assert pl2["status"] == "ok", pl2.get("reason")
    qc, r = fetch_process_combine(tmp_path, dest, y, ref, wt, refp, data, monkeypatch)
    assert qc["status"] == "ok" and set(qc["units"]) == {"Abcz1", "Defz2"}
    assert all(u["status"] == "ok" for u in qc["units"].values()), {g: u.get("reason") for g, u in qc["units"].items()}
    assert ers_of["Defz2_wt1"] not in qc["libraries"] and ers_of["Defz2_wt2"] in qc["libraries"]
    assert r["tables"] == ["bulk_all", "bulk_study_MGP-DMDD"]
    # the stored ties edited after the re-plan: the pin check refuses
    ds = cm.dataset(ACC, y)
    pl3 = json.loads(json.dumps(pl2))
    pl3["line_set"]["Defz2"]["control"].insert(0, "Defz2_wt1")
    with pytest.raises(cm.Refused, match="line set"):
        ae.verify_pins(ds, pl3)


def test_r2_blocker6_a_line_that_drops_under_two_in_an_arm_leaves_the_set(tmp_path, y, monkeypatch):
    wt, refp, ref = make_ref(tmp_path)
    dest = tmp_path / "dest"
    sdrf, runs_of, ers_of = sdrf_and_runs()
    fe.plan(ACC, dest, y, metas=texts(sdrf))
    data = drop_columns(matrix(ref, runs_of, np.random.default_rng(3)), {runs_of["Defz2_hom1"][-1]})
    monkeypatch.setattr(fe, "_open", opener(data, []))
    assert bulk.header(ACC, dest, datasets=y)["status"] == "ok"
    pl = fe.plan(ACC, dest, y, metas=texts(sdrf))
    assert pl["status"] == "ok" and set(pl["line_set"]) == {"Abcz1"}
    with pytest.raises(cm.Refused, match="not a line of the committed rule"):
        ae.pin_template(ACC, dest, "ko", "Defz2", ["Defz2_hom2"], y)
    pin_all(dest, y, pl)
    assert fe.plan(ACC, dest, y, metas=texts(sdrf))["status"] == "ok"


def test_r2_blocker4_a_matrix_named_as_normalised_is_refused_at_the_header_step_and_the_read(tmp_path, y,
                                                                                             monkeypatch):
    assert not bulk.NORMALISED.search(ae.matrix_name(ACC))   # the one name the resource offers
    sdrf, runs_of, _ = sdrf_and_runs()
    pl = fe.plan(ACC, tmp_path / "d", y, metas=texts(sdrf))
    calls: list = []
    monkeypatch.setattr(fe, "_open", opener(b"Gene ID\tERR1\n", calls))
    for name in ("E-MTAB-7322-normalized-counts.tsv", "E-MTAB-7322-expected_counts.tsv"):
        monkeypatch.setattr(ae, "matrix_name", lambda acc, name=name: name)
        with pytest.raises(cm.Refused, match="normalised or expected counts"):
            ae.header(ACC, tmp_path / "d", y)
        with pytest.raises(cm.Refused, match="normalised or expected counts"):
            ae.record_header(cm.dataset(ACC, y), tmp_path / "d", pl, name, ae.urls(ACC)["matrix"], b"Gene ID\tERR1")
        with pytest.raises(cm.Refused, match="normalised or expected counts"):
            bulk.read_series_table(tmp_path / "absent", name, ["Gene ID"], ["ERR1"], "0" * 64)
    assert calls == []


# ---------------------------------------------------------------------------------------------------------------------
# what the resource's real description and sample tables showed (2026-11-20, metadata only), and two closures


def test_a_description_key_with_a_hyphen_is_read_and_other_punctuation_still_withholds(y):
    idf = IDF + b"MAGE-TAB Version\t1.1\n"
    pl = ae.plan_from_texts(cm.dataset(ACC, y), idf, sdrf_and_runs()[0], METHODS)
    assert pl["status"] == "ok" and pl["series"][ACC]["fields"]["idf|MAGE-TAB Version"] == ["1.1"]
    for bad in (b"Bad;Key\tx\n", b"Key/Other\tx\n", b"-Leading\tx\n", b"Key(1)\tx\n", b"Key--\tx\n", b"Key-\tx\n",
                b"A - B\tx\n"):
        got = ae.plan_from_texts(cm.dataset(ACC, y), IDF + bad, sdrf_and_runs()[0], METHODS)
        assert got["status"] == "refused" and "closed grammar" in got["reason"] and "text withheld" in got["reason"]


def test_an_allele_written_with_a_caret_passes_the_name_comparison_and_a_wrong_prefix_does_not(y):
    sdrf = sdrf_and_runs()[0].decode()
    caret = sdrf.replace("<tm1a(EUCOMM)Wtsi>", "^tm1a(EUCOMM)Wtsi")
    assert caret != sdrf and "Abcz1^tm1a(EUCOMM)Wtsi" in caret
    base = ae.plan_from_texts(cm.dataset(ACC, y), IDF, sdrf.encode(), METHODS)
    got = ae.plan_from_texts(cm.dataset(ACC, y), IDF, caret.encode(), METHODS)
    assert got["status"] == "ok" and got["line_set"] == base["line_set"]
    assert got["comparisons"]["c"] == base["comparisons"]["c"]
    wrong = caret.replace("Abcz1^tm1a(EUCOMM)Wtsi", "Abcz^tm1a(EUCOMM)Wtsi")
    bad = ae.plan_from_texts(cm.dataset(ACC, y), IDF, wrong.encode(), METHODS)
    assert len(bad["comparisons"]["c"]) > len(base["comparisons"]["c"]) and "Abcz1" not in bad["line_set"]


def test_the_named_resource_must_withhold_its_held_line(tmp_path):
    for text, ok in ((YAML, True),
                     (YAML.replace("      - {line: Trps1, reason: held}\n", ""), False),
                     (YAML.replace("{line: Trps1, reason: held}", "{line: Trps1, reason: pending_phenocopy}"), False)):
        p = tmp_path / f"d{ok}{len(text)}.yaml"
        p.write_text(text)
        if ok:
            assert cm.dataset(ACC, p)["withheld_lines"][0] == {"line": "Trps1", "reason": "held"}
        else:
            with pytest.raises(cm.Refused, match="withheld_lines must hold"):
                cm.load_datasets(p)


def test_the_geo_bulk_forms_are_closed_in_the_committed_code(tmp_path, monkeypatch):
    assert "GEO_BULK_OPEN = False" in (EXT / "t3x_common.py").read_text()
    monkeypatch.setattr(cm, "GEO_BULK_OPEN", False)
    geo = YAML.replace("repository: arrayexpress", "repository: geo")
    geo = geo.replace("count_form: arrayexpress", "count_form: per_sample")
    geo = geo.replace(f"accession: {ACC}", "accession: GSE910001")
    geo = geo.replace(f"screen_series: [{ACC}]", "screen_series: [GSE910001]")
    p = tmp_path / "geo.yaml"
    p.write_text(geo)
    with pytest.raises(cm.Refused, match="the GEO bulk forms are closed"):
        cm.load_datasets(p)
    # the arrayexpress form is not affected
    q = tmp_path / "ae.yaml"
    q.write_text(YAML)
    assert cm.dataset(ACC, q)["count_form"] == "arrayexpress"


def test_trailing_whitespace_after_the_allele_part_is_not_compared_but_two_prefixes_fail(y):
    sdrf = sdrf_and_runs()[0].decode().replace("<tm1a(EUCOMM)Wtsi>", "^tm1a(EUCOMM)Wtsi")
    base = ae.plan_from_texts(cm.dataset(ACC, y), IDF, sdrf.encode(), METHODS)
    lines = sdrf.split("\n")
    line = sorted(base["line_set"])[0]
    sn = base["line_set"][line]["ko"][0]                         # an eligible knockout embryo of the committed set
    rows_of = [i for i, ln in enumerate(lines) if ln.startswith(sn + "\t")]
    assert len(rows_of) >= 2 and sn not in base["comparisons"]["c"]
    k = rows_of[0]
    spaced = lines.copy()
    # a trailing space on one row of the embryo only
    spaced[k] = spaced[k].replace(f"{line}^tm1a(EUCOMM)Wtsi", f"{line}^tm1a(EUCOMM)Wtsi ", 1)
    assert spaced[k] != lines[k]
    got = ae.plan_from_texts(cm.dataset(ACC, y), IDF, "\n".join(spaced).encode(), METHODS)
    assert got["line_set"] == base["line_set"] and got["comparisons"]["c"] == base["comparisons"]["c"]
    two = lines.copy()
    two[k] = two[k].replace(f"{line}^tm1a(EUCOMM)Wtsi", f"{line}x^tm1a(EUCOMM)Wtsi", 1)
    bad = ae.plan_from_texts(cm.dataset(ACC, y), IDF, "\n".join(two).encode(), METHODS)
    assert sn in bad["comparisons"]["c"]


def test_a_row_naming_a_withheld_line_in_its_source_name_or_allele_is_withheld_too(y):
    sdrf, _, _ = sdrf_and_runs()
    lines = sdrf.decode().split("\n")
    k = next(i for i, ln in enumerate(lines) if ln.startswith("Abcz1_wt1\t"))
    row = lines[k].split("\t")
    by_name = row.copy()
    by_name[0] = "Trps1_hom9"                                   # its line column still says Abcz1
    by_allele = row.copy()
    by_allele[0], by_allele[H.index("Comment[allele]")] = "Abcz1_wt9", "Trps1<x>"
    for extra in (by_name, by_allele):
        e = extra.copy()
        e[H.index("Comment[ENA_RUN]")] = "ERR7777777"
        e[H.index("Comment[ENA_SAMPLE]")] = "ERS7777777"
        e[H.index("Characteristics[individual]")] = "new-individual"
        text = ("\n".join([*lines[:-1], "\t".join(e)]) + "\n").encode()
        pl = ae.plan_from_texts(cm.dataset(ACC, y), IDF, text, METHODS)
        assert pl["status"] == "ok" and e[0] not in {r["source_name"] for r in pl["tables"][ACC]}
        assert "ERR7777777" not in json.dumps(pl)

