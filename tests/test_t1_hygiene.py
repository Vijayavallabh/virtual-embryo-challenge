"""R2's gene modules (experiments/t1_hygiene/hygiene.py), their config (configs/t1_hygiene.yaml) and the extractor of
the external lists (experiments/external/hygiene_sources.py), on synthetic inputs and the committed extracts.

No gene symbol is typed here: the synthetic symbols are made up (and say so), and the committed extracts are checked
by their counts, their sha256 and Table 1's place inside Table 5, never against a typed list."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
HYG_DIR = ROOT / "experiments" / "t1_hygiene"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def hyg():
    return _load("t1_hygiene_hygiene", HYG_DIR / "hygiene.py")


@pytest.fixture(scope="module")
def src():
    return _load("external_hygiene_sources", ROOT / "experiments" / "external" / "hygiene_sources.py")


@pytest.fixture(scope="module")
def cfg(hyg):
    return hyg.load_config()


# ---------------------------------------------------------------------------------------------------------------------
# The config holds R2's and N12's registered numbers


def test_config_constants_are_r2s_registered_numbers(cfg):
    sys.path.insert(0, str(ROOT / "experiments" / "t1_m1"))
    import feng_dirs

    assert (cfg["gate"]["list_size"], cfg["gate"]["upper"], cfg["gate"]["floor"]) == (400, 133, 10)
    assert cfg["n12"]["reported_sizes"] == [100, 1000]
    null = cfg["n12"]["decile_null"]
    assert (null["draws"], null["seed"], null["bins"]) == (1000, 20261106, 10)
    d = cfg["direction"]
    assert d["fraction"] == 0.4 == pytest.approx(feng_dirs.bracket_fraction(d["from_day"], d["to_day"]))
    assert (d["from_day"], d["to_day"], d["pseudocount_cp10k"]) == (10.5, 12.5, 0.5)
    assert d["pseudocount_cp10k"] == feng_dirs.load_config()["metrics"]["pseudocount_cp10k"]
    f = cfg["n2_filter"]
    assert (f["type_of_gene"], f["min_nuclear_prop"], f["clusters_nuc"], f["clusters_cell"]) == ("ncRNA", 0.9, 11, 11)
    assert {k: v["n"] for k, v in cfg["extracts"].items()} == {"D": 32, "D_table5": 140, "N2": 19}
    assert cfg["patterns"] == {"A1": "^Hb[ab]-", "N1": "^mt-"} and cfg["lines"] == {"A1": True, "N1": True}
    assert cfg["n12"]["reported_patterns"] == {"ribosomal": "^Rp[ls][0-9]",
                                               "oxphos": "^(Nduf|Sdh[a-d]|Uqcr|Cox[4-8]|Atp5)", "mt": "^mt-"}
    union = ROOT / cfg["union"]["path"]
    assert cfg["union"]["path"] == "experiments/t1_hygiene/results/n12_hygiene_size/hygiene_genes.csv"
    pinned = cfg["union_sha256"]  # null until the build commit after N12; then the committed union's sha256
    assert pinned is None or (re.fullmatch(r"[0-9a-f]{64}", pinned)
                              and hashlib.sha256(union.read_bytes()).hexdigest() == pinned)


def test_config_sources_are_the_registered_files(cfg):
    """Each source's sha256 and size are the ones its EXTERNAL_DATA.md row records."""
    rows = {ln.split("|")[1].strip(): ln for ln in (ROOT / "EXTERNAL_DATA.md").read_text().splitlines()
            if ln.startswith("| X0")}
    for s in cfg["sources"].values():
        row = rows[s["register"]]
        assert s["sha256"] in row and f"{s['bytes']:,} bytes" in row and s["url"] in row


def test_committed_extracts_match_their_pinned_sha256_counts_and_subset(hyg, cfg):
    for ent in cfg["extracts"].values():
        assert hashlib.sha256((ROOT / ent["path"]).read_bytes()).hexdigest() == ent["sha256"]
    d, t5, n2 = (hyg.read_extract(m, cfg) for m in ("D", "D_table5", "N2"))
    assert (len(d), len(t5), len(n2)) == (32, 140, 19)
    assert set(d) <= set(t5)
    for ent in cfg["extracts"].values():  # each header names the source, its sha256, the extractor and the date
        head = dict(ln[2:].split(": ", 1) for ln in (ROOT / ent["path"]).read_text().splitlines() if ln[:2] == "# ")
        src = cfg["sources"][ent["source"]]
        assert head["source_sha256"] == src["sha256"] and head["source_bytes"] == str(src["bytes"])
        assert head["extractor"].startswith("experiments/external/hygiene_sources.py (sha256 ")
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", head["date"]) and head["rows"] == str(ent["n"])
        assert ("pdftotext version" in head["tool"]) == (ent["source"] == "vdb2017_si")


def test_extracts_hold_only_the_registered_columns(src, cfg):
    """X007 registers van den Brink's gene symbol columns only (Table 1, Table 5); X008 the S6 columns N2's filter
    reads. The committed extracts hold those and nothing else."""
    rows = {ln.split("|")[1].strip(): ln for ln in (ROOT / "EXTERNAL_DATA.md").read_text().splitlines()
            if ln.startswith("| X0")}
    columns = {m: src.table_rows((ROOT / e["path"]).read_text())[0] for m, e in cfg["extracts"].items()}
    assert columns["D"] == columns["D_table5"] == ["gene_symbol"]
    assert "column 'Gene symbol' (32 genes)" in rows["X007"] and "column 'Gene Symbol' (140 genes)" in rows["X007"]
    assert columns["N2"] == ["gene_symbol", *src.N2_COLUMNS[1:]]
    assert f"columns {', '.join(src.N2_COLUMNS)} |" in rows["X008"]


# ---------------------------------------------------------------------------------------------------------------------
# hygiene.py


def test_importing_hygiene_reads_no_file_and_changes_nothing():
    code = ("import sys\n"
            "seen = []\n"
            "def hook(ev, a):\n"
            "    if ev == 'open' and a and isinstance(a[0], str):\n"
            "        seen.append(a[0])\n"
            "sys.addaudithook(hook)\n"
            f"sys.path.insert(0, {str(HYG_DIR)!r})\n"
            "before = list(sys.path)\n"
            "import hygiene\n"
            "assert sys.path == before, 'sys.path changed'\n"
            "data = [p for p in seen if not p.endswith(('.py', '.pyc', '.so')) and '__pycache__' not in p\n"
            "        and not p.endswith(('/', 'site-packages')) and '.dist-info' not in p and '.pth' not in p]\n"
            "assert not data, data\n")
    r = subprocess.run([sys.executable, "-I", "-B", "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_modules_map_by_exact_case_sensitive_symbol_and_list_the_unmapped(hyg):
    t1 = ["SynA", "synB", "SynC", "Hba-syn1", "mt-Syn1"]
    m = hyg.map_symbols(["SynA", "SYNA", "SynB", "SynC ", "SynC", "Ghost-syn"], t1)
    assert m == {"mapped": ["SynA", "SynC"], "unmapped": ["SYNA", "SynB", "SynC ", "Ghost-syn"]}
    mods = hyg.map_modules({"D": ["SynC", "SynA", "Nope-syn"], "A1": ["Hba-syn1"]}, t1)
    assert mods["D"] == {"source": ["SynC", "SynA", "Nope-syn"], "mapped": ["SynC", "SynA"], "unmapped": ["Nope-syn"]}
    assert mods["A1"]["unmapped"] == []


def test_pattern_modules_hold_every_matching_feng_symbol_once(hyg, cfg):
    """R2: A1 = Feng symbols matching ^Hb[ab]-, N1 = Feng symbols beginning mt-, with no other rule: a symbol Feng
    uses more than once is listed once, not dropped, although dirlib.sym_index gives it no Feng value."""
    sys.path.insert(0, str(ROOT / "experiments" / "t1_direction"))
    import dirlib

    feng = ["Hba-syn1", "Hbb-syn2", "Hbb-syn2", "Hbq-syn3", "hba-syn4", "Hbaa-syn5", "mt-Syn1", "Mt-Syn2", "MT-Syn3",
            "mt-Syn4", "mt-Syn4", "Syn-mt-6", "Plain-syn"]
    a1, n1 = hyg.pattern_module(feng, cfg["patterns"]["A1"]), hyg.pattern_module(feng, cfg["patterns"]["N1"])
    assert a1 == ["Hba-syn1", "Hbb-syn2"] and n1 == ["mt-Syn1", "mt-Syn4"]  # each once; case and anchors respected
    assert dirlib.sym_index(feng, a1 + n1).tolist() == [0, -1, 6, -1]  # the repeated ones have no Feng value


def test_build_modules_reads_the_extracts_and_the_patterns(hyg, cfg):
    feng = ["Hba-syn1", "mt-Syn1", "Plain-syn"]
    mods = hyg.build_modules(feng, cfg, table5=True)
    assert list(mods) == ["D", "A1", "N1", "N2", "D_table5"]
    assert (mods["A1"], mods["N1"]) == (["Hba-syn1"], ["mt-Syn1"])
    assert [len(mods[m]) for m in ("D", "N2", "D_table5")] == [32, 19, 140]
    assert "D_table5" not in hyg.build_modules(feng, cfg)


def test_read_extract_refuses_a_file_whose_sha256_is_not_pinned(hyg, cfg, tmp_path):
    ent = cfg["extracts"]["N2"]
    text = (ROOT / ent["path"]).read_text()
    (tmp_path / Path(ent["path"]).parent).mkdir(parents=True)
    (tmp_path / ent["path"]).write_text(text + "Extra-syn\tncRNA\t2\t9\t2\t9\t1\n")  # one made-up row more
    with pytest.raises(hyg.Refused, match="sha256"):
        hyg.read_extract("N2", cfg, root=tmp_path)
    c = copy.deepcopy(cfg)
    c["extracts"]["N2"]["sha256"] = hashlib.sha256((tmp_path / ent["path"]).read_bytes()).hexdigest()
    with pytest.raises(hyg.Refused, match="registers 19"):  # pinned, but not the registered number of rows
        hyg.read_extract("N2", c, root=tmp_path)


def test_union_follows_the_lines_and_names_each_genes_modules(hyg, cfg):
    mapped = {"D": {"mapped": ["SynB", "SynA"]}, "A1": {"mapped": ["Hba-syn1"]}, "N1": {"mapped": ["mt-Syn1"]},
              "N2": {"mapped": ["SynA", "Nc-syn"]}}
    assert hyg.union_modules(cfg) == ["D", "A1", "N1", "N2"]
    u = hyg.union(mapped, hyg.union_modules(cfg))
    assert u == {"Hba-syn1": ["A1"], "Nc-syn": ["N2"], "SynA": ["D", "N2"], "SynB": ["D"], "mt-Syn1": ["N1"]}
    dropped = copy.deepcopy(cfg)
    dropped["lines"]["A1"] = False  # the captain drops a line
    assert hyg.union_modules(dropped) == ["D", "N1", "N2"]
    assert "Hba-syn1" not in hyg.union(mapped, hyg.union_modules(dropped))


def test_load_union_refuses_unless_its_sha256_is_pinned(hyg, cfg, tmp_path):
    members = {"Hba-syn1": ["A1"], "SynA": ["D", "N2"], "mt-Syn1": ["N1"]}
    path = tmp_path / "hygiene_genes.csv"
    path.write_text(hyg.union_text(members, ["synthetic union", "commit x"]))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(hyg.Refused, match="pins no union_sha256"):  # null until the build commit after N12
        hyg.load_union(path, dict(cfg, union_sha256=None))
    wrong = copy.deepcopy(cfg)
    wrong["union_sha256"] = "0" * 64
    with pytest.raises(hyg.Refused, match="is not the union_sha256"):
        hyg.load_union(path, wrong)
    pinned = copy.deepcopy(cfg)
    pinned["union_sha256"] = sha
    u = hyg.load_union(path, pinned, root=tmp_path)
    assert u == {"file": "hygiene_genes.csv", "sha256": sha, "genes": ["Hba-syn1", "SynA", "mt-Syn1"],
                 "modules": {"D": 1, "A1": 1, "N1": 1, "N2": 1},
                 "sources": {k: v["sha256"] for k, v in cfg["sources"].items()}}
    path.write_text(path.read_text() + "Late-syn,D\n")  # changed after pinning
    with pytest.raises(hyg.Refused, match="is not the union_sha256"):
        hyg.load_union(path, pinned)
    with pytest.raises(hyg.Refused, match="no union file"):
        hyg.load_union(tmp_path / "absent.csv", pinned)
    with pytest.raises(hyg.Refused, match="cannot hold"):
        hyg.union_text({"Bad,syn": ["D"]}, [])


def test_the_repository_union_is_refused_while_no_sha256_is_pinned(hyg, cfg):
    if cfg["union_sha256"] is None:
        with pytest.raises(hyg.Refused, match="pins no union_sha256"):
            hyg.load_union()


# ---------------------------------------------------------------------------------------------------------------------
# The extractor's parsers, on synthetic text and a synthetic xlsx


TEXT = """Body text that cites Supplementary Table 1. and Supplementary Table 5. in passing.
                 Supplementary table 1. List of all genes that are upregulated* in cluster 2-cells
                 *Differential gene expression analysis. Only the 3 genes with a p-value below 0.001 are shown.
                    Gene            Fold       P-value
                     Syn-a      3.09         1.89E-139
                   Syn-b2       0.85         8.15E-06
                    Syn-c       -1.10        2.00E-48
Nature Methods: doi:10.1038/nmeth.4437
                 Supplementary table 2. List of smFISH probes.
                 Supplementary table 5. List of genes of which expression levels were assumed to be
                 affected by the dissociation procedure.
                  Gene Symbol
                 Syn-a
                 Syn-x
Nature Methods: doi:10.1038/nmeth.4437
                 Syn-b2
                 Syn-c
                 *Differential gene expression analysis was performed with DESeq.
                 Syn-after-footnote
                 Supplementary table 6. List of genes differentially expressed
"""


def test_vdb_tables_are_read_between_their_captions(src):
    assert src.vdb_table1(TEXT) == [("Syn-a", "3.09", "1.89E-139"), ("Syn-b2", "0.85", "8.15E-06"),
                                    ("Syn-c", "-1.10", "2.00E-48")]
    assert src.vdb_table5(TEXT) == ["Syn-a", "Syn-x", "Syn-b2", "Syn-c"]
    with pytest.raises(SystemExit, match="caption states 4"):
        src.vdb_table1(TEXT.replace("Only the 3 genes", "Only the 4 genes"))


def _xlsx(path: Path, rows: list[list[str]]) -> Path:
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rel = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    strings = sorted({v for r in rows for v in r if not re.fullmatch(r"-?[0-9.]+", v)})
    sst = f'<sst {ns}>' + "".join(f"<si><t>{s}</t></si>" for s in strings) + "</sst>"

    def cell(col, i, v):
        ref = f"{chr(65 + col)}{i + 1}"
        if v == "":
            return ""  # an empty cell is left out, as Excel does
        if v in strings:
            return f'<c r="{ref}" t="s"><v>{strings.index(v)}</v></c>'
        return f'<c r="{ref}"><v>{v}</v></c>'

    sheet = (f"<worksheet {ns}><sheetData>" + "".join(
        f'<row r="{i + 1}">' + "".join(cell(c, i, v) for c, v in enumerate(r)) + "</row>" for i, r in enumerate(rows))
             + "</sheetData></worksheet>")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/workbook.xml", f'<workbook {ns} {rel}><sheets><sheet name="S" sheetId="1" r:id="rId1"/>'
                                      "</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/worksheets/sheet1.xml", sheet)
        z.writestr("xl/sharedStrings.xml", sst)
    return path


def test_bakken_rows_pass_n2s_filter_only_when_every_condition_holds(src, cfg, tmp_path):
    head = ["gene", "GeneID", "type_of_gene", "inh.clusters_nuc", "exc.clusters_nuc", "inh.clusters_cell",
            "exc.clusters_cell", "nuclear.prop"]
    rows = [head,
            ["Keep-syn1", "1", "ncRNA", "2", "9", "2", "9", "0.95"],
            ["Keep-syn2", "2", "ncRNA", "2", "9", "2", "9", "0.9"],  # at the threshold
            ["Low-syn", "3", "ncRNA", "2", "9", "2", "9", "0.89"],
            ["NA-syn", "4", "ncRNA", "2", "9", "2", "9", "NA"],  # never passes
            ["Coding-syn", "5", "protein-coding", "2", "9", "2", "9", "0.99"],
            ["Nuc10-syn", "6", "ncRNA", "1", "9", "2", "9", "0.99"],
            ["Cell10-syn", "7", "ncRNA", "2", "9", "2", "8", "0.99"],
            ["Empty-syn", "8", "ncRNA", "2", "9", "", "9", "0.99"]]
    columns, table = src.xlsx_rows(_xlsx(tmp_path / "s6.xlsx", rows))
    assert columns == head and table[0]["gene"] == "Keep-syn1" and "inh.clusters_cell" not in table[-1]
    kept = src.bakken_n2(columns, table, cfg["n2_filter"])
    assert [r[0] for r in kept] == ["Keep-syn1", "Keep-syn2"]
    assert kept[0] == ["Keep-syn1", "ncRNA", "2", "9", "2", "9", "0.95"]  # the extract keeps the values as stored
    with pytest.raises(SystemExit, match="missing"):
        src.bakken_n2(columns[:-1], table, cfg["n2_filter"])


def test_extractor_refuses_a_source_that_is_not_the_registered_file(src, cfg, tmp_path):
    fake = tmp_path / "si.pdf"
    fake.write_bytes(b"not the registered file")
    with pytest.raises(SystemExit, match="refused: not the registered file"):
        src.check_source(fake, cfg["sources"]["vdb2017_si"])


def test_extract_files_round_trip(src):
    text = src.extract_text([("source", "synthetic")], ["gene_symbol", "x"], [["Syn-a", "1"], ["Syn-b", "2"]])
    assert text.startswith("# source: synthetic\n")
    assert src.table_rows(text) == (["gene_symbol", "x"], [["Syn-a", "1"], ["Syn-b", "2"]])
    with pytest.raises(SystemExit, match="tab or a line break"):
        src.extract_text([], ["gene_symbol"], [["Syn\tbad"]])


def test_config_is_valid_yaml_with_no_symbol_lists():
    raw = (ROOT / "configs" / "t1_hygiene.yaml").read_text()
    cfg = yaml.safe_load(raw)
    assert set(cfg) == {"sources", "extracts", "n2_filter", "patterns", "lines", "direction", "gate", "n12", "union",
                        "union_sha256"}
    assert not any(isinstance(v, list) and len(v) > 3 for v in _walk(cfg))  # no list of symbols hides in it


def _walk(x):
    if isinstance(x, dict):
        for v in x.values():
            yield v
            yield from _walk(v)
    elif isinstance(x, list):
        for v in x:
            yield from _walk(v)
