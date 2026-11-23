"""experiments/t1_ext/load.py on synthetic files only: window samples removed FIRST (before any file is opened), a
stageless sample refuses the atlas, the QC record proves zero cells in the window, no label or embedding column is ever
read (synthetic tables whose forbidden columns raise when touched), the per-file gates (aggr suffixes, non-integer
counts, Multiome feature types), the dge reader, and that no step writes inside the repository. Third version: the
manifest stores each sample's full record, and the loader refuses any sample that is not pinned or whose record, end
or stage is not its pin's, before any file is opened."""

from __future__ import annotations

import copy
import gzip
import importlib.util
import io
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import io as spio
from scipy import sparse

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


ld = _load("load")
ACC = "GSE999002"


def config():
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["atlases"][ACC] = {"status": "use", "wave": 1, "series": [ACC], "form": "tenx_mtx",
                           "early": {"labels": ["E8.5"], "tissue": "embryo"},
                           "late": {"labels": ["E14.5"], "tissue": "heart"},
                           "assay_rule": {"library_strategy": "RNA-Seq"}, "early_tissue": "whole_embryo",
                           "groups_formed": ["CM"], "qc": {"min_umi": 5, "min_genes": 2, "max_mito": 0.5},
                           "protocol_terms_cleared": {}, "protocol_gate_cleared": {}}
    return cfg


class Forbidden:
    """A synthetic table: whitelisted columns answer, every other column raises when touched."""

    def __init__(self, cols: dict, forbidden=("cell_type", "cluster", "umap_1", "authors_label")):
        self.cols, self.forbidden = cols, forbidden

    def names(self):
        return (*self.cols, *self.forbidden)

    def column(self, name):
        if name in self.forbidden:
            raise AssertionError(f"forbidden column {name} was touched")
        return np.asarray(self.cols[name], dtype=object)


def counts(n_genes, n_cells, seed=0, mito=None):
    rng = np.random.default_rng(seed)
    M = rng.poisson(2.0, size=(n_genes, n_cells)).astype(np.int64)
    return sparse.csc_matrix(M)


def genes(n):
    return [f"Gene{i}" for i in range(n - 1)] + ["mt-Co1"]


def test_forbidden_columns_are_never_touched():
    rd = ld.Reader()
    g = genes(6)
    bt = Forbidden({"barcode": [f"AAAC{i}-1" for i in range(4)]})
    ft = Forbidden({"id": [f"ENSMUSG{i}" for i in range(6)], "symbol": g,
                    "feature_type": ["Gene Expression"] * 5 + ["Peaks"]})
    x = ld.tenx_from_tables(bt, ft, counts(6, 4), True, "-1", rd)
    assert x["C"].shape == (4, 5) and list(x["genes"]) == g[:5]
    assert set(rd.read) == {"barcodes:barcode", "features:id", "features:symbol", "features:feature_type"}
    with pytest.raises(cm.Refused):
        rd.column(ft, "features", "cell_type")
    with pytest.raises(cm.Refused):
        rd.column(ft, "features", "umap_1")


def test_dge_refuses_annotation_columns_without_reading_them():
    class Dge(Forbidden):
        def counts(self):
            return np.ones((3, 2))

    rd = ld.Reader()
    t = Dge({"gene": ["A", "B", "C"]}, forbidden=("cluster",))
    t.names = lambda: ("gene", "ACGTACGT", "cluster")
    with pytest.raises(cm.Refused):
        ld.dge_from_table(t, rd)
    ok = Dge({"gene": ["A", "B", "C"]}, forbidden=())
    ok.names = lambda: ("gene", "ACGTACGT", "TTGCAACG")
    x = ld.dge_from_table(ok, rd)
    assert x["C"].shape == (2, 3)


def test_gates_suffix_integers_and_multiome():
    rd = ld.Reader()
    g = genes(4)
    ft = Forbidden({"id": list("abcd"), "symbol": g, "feature_type": ["Gene Expression"] * 4})
    with pytest.raises(cm.Refused):  # an aggr matrix renumbers libraries -2, -3
        ld.tenx_from_tables(Forbidden({"barcode": ["A-1", "B-2"]}), ft, counts(4, 2), False, "-1", rd)
    M = counts(4, 2).astype(np.float64)
    M.data = M.data + 0.5
    with pytest.raises(cm.Refused):  # corrected or normalised values
        ld.tenx_from_tables(Forbidden({"barcode": ["A-1", "B-1"]}), ft, M, False, "-1", rd)

    class TwoCols(Forbidden):
        def names(self):
            return ("id", "symbol")

    with pytest.raises(cm.Refused):  # a Multiome features file must carry feature_type
        ld.tenx_from_tables(Forbidden({"barcode": ["A-1", "B-1"]}), TwoCols({"id": list("abcd"), "symbol": g}),
                            counts(4, 2), True, "-1", rd)


def test_dge_file_reader(tmp_path):
    p = tmp_path / "GSM1_E9.5DayEmbryo_dge.txt.gz"
    p.write_bytes(gzip.compress(b'"AAAC"\t"GGTT"\nMyl7\t3\t0\nTnnt2\t1\t2\nmt-Co1\t0\t1\n'))
    t = ld.DgeFile(p)
    x = ld.dge_from_table(t, ld.Reader())
    assert x["C"].toarray().tolist() == [[3, 1, 0], [0, 2, 1]]
    assert list(x["barcodes"]) == ["AAAC", "GGTT"]


E05 = ["The day of the vaginal plug was designated E0.5."]
SERIES_RECORDS = {ACC: [["title", ["synthetic"]], ["summary", ["Hearts."]],
                        ["overall_design", ["One library per sample."]]]}


def rec(title, label, files=(), extra=(), series_protocol=E05, **kw):
    """A stored sample entry as fetch.py writes it into the fetch manifest: the full record, in order, as pairs."""
    gsm = kw.pop("gsm", "GSM1")
    record = [["title", [title]], *[[k, list(v)] for k, v in extra],
              ["data_processing", ["Cell Ranger count per sample"]]]
    if files:
        record.append(["supplementary_file_1", [f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/x/{f}" for f in files]])
    return {"series": ACC, "end": "early" if label == "E8.5" else "late", "label": label, "record": record,
            "record_sha256": cm.record_sha256("SAMPLE", gsm, record), "series_protocol": list(series_protocol), **kw}


def test_window_filter_first_and_stageless_refuses():
    samples = {"GSM1": rec("E8.5 embryo", "E8.5"), "GSM2": rec("E12.5 heart", "E12.5")}
    kept, wrec = ld.window_filter(samples)
    assert sorted(kept) == ["GSM1"] and [r["gsm"] for r in wrec["samples_removed_window"]] == ["GSM2"]
    with pytest.raises(cm.Refused):
        ld.window_filter({**samples, "GSM3": rec("fetal heart", "E14.5")})
    with pytest.raises(cm.Refused):  # a manifest without the sample's own record
        ld.window_filter({"GSM4": {"label": "E14.5"}})
    with pytest.raises(cm.Refused):  # the label disagrees with the record
        ld.window_filter({"GSM5": rec("E16.5 heart", "E14.5")})
    with pytest.raises(cm.Refused):  # the loader re-reads the whole record: a window stage under another key
        ld.window_filter({"GSM6": rec("E14.5 heart", "E14.5", extra=[("characteristics_ch1", ["tissue: E12.5"])])})
    with pytest.raises(cm.Refused):  # ... and in the sample's own protocol text (the verifier's blocker 3)
        ld.window_filter({"GSM6": rec("E14.5 heart", "E14.5",
                                      extra=[("extract_protocol_ch1", ["Hearts were dissected from E12.5 embryos."])])})
    # and the convention: plug = E1 puts E14.5 at [13.5, 14.0], so the sample is removed before any file is opened
    kept, wrec = ld.window_filter({"GSM7": rec("E14.5 heart", "E14.5", series_protocol=["Plug day was E1."])})
    assert not kept and wrec["samples_removed_window"][0]["low"] == 13.5
    assert ld.drop_window_cells(np.array([8.5, 9.5, 10.5, 13.5, 14.5])).tolist() == [True, True, False, False,
                                                                                       True]
    with pytest.raises(cm.Refused):
        ld.drop_window_cells(np.array([8.5, np.nan]))


def write_10x(folder: Path, gsm: str, n_genes: int, n_cells: int, seed: int) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    g = genes(n_genes)
    files = {}
    bc = "\n".join(f"BC{seed}{i:04d}-1\tauthors_cluster_{i % 3}" for i in range(n_cells)) + "\n"
    ft = "\n".join(f"ENSMUSG{i}\t{s}\tGene Expression" for i, s in enumerate(g)) + "\n"
    buf = io.BytesIO()
    spio.mmwrite(buf, counts(n_genes, n_cells, seed))
    for role, name, data in (("barcodes", f"{gsm}_barcodes.tsv.gz", bc.encode()),
                             ("features", f"{gsm}_features.tsv.gz", ft.encode()),
                             ("matrix", f"{gsm}_matrix.mtx.gz", buf.getvalue())):
        (folder / name).write_bytes(gzip.compress(data))
        files[role] = {"name": name, "sha256": cm.sha256_file(folder / name)}
    return files


def raw_manifest(dest: Path, cfg, samples: dict) -> None:
    rdir = dest / "raw" / ACC
    rdir.mkdir(parents=True, exist_ok=True)
    man = {"accession": ACC, "atlas_entry_sha256": cm.entry_sha256(ACC, cfg), "chosen": {"early": "E8.5",
                                                                                         "late": "E14.5"},
           "samples": samples, "tasks": ["T1"], "series_records": SERIES_RECORDS}
    (rdir / "MANIFEST.json").write_text(json.dumps(man))


def pin_for(cfg, samples) -> dict:
    """The atlas entry with every manifest sample and the series pinned by their records' hashes."""
    cfg = copy.deepcopy(cfg)
    who = {"read_by": ["Reader A", "Reader B"], "date": "2026-11-21"}
    cfg["atlases"][ACC]["pinned"] = {
        "samples": [{"accession": g, "end": s["end"], "stage": s["label"],
                     "record_sha256": cm.record_sha256("SAMPLE", g, s["record"]), **who} for g, s in samples.items()],
        "series": [{"accession": ACC, "record_sha256": cm.record_sha256("SERIES", ACC, SERIES_RECORDS[ACC]), **who}]}
    cm.check_pins(ACC, cm.atlas(ACC, cfg))
    return cfg


def clear_for(cfg, samples) -> dict:
    """The atlas entry pinned (pin_for) and with processing_cleared set to the hash of the text the manifest's samples
    give."""
    cfg = pin_for(cfg, samples)
    text = cm.processing_text(ACC, [(ACC, SERIES_RECORDS[ACC])],
                              {g: {"series": s["series"], "text": cm.sample_text(s["record"])}
                               for g, s in samples.items()})
    cfg["atlases"][ACC]["processing_cleared"] = {"sha256": cm.text_sha256(text), "read_by": "test"}
    return cfg


def two_samples(tmp_path):
    samples = {}
    for gsm, title, label, seed in (("GSM11", "E8.5 embryo", "E8.5", 1), ("GSM12", "E14.5 heart", "E14.5", 2)):
        files = write_10x(tmp_path / "raw" / ACC / gsm, gsm, 30, 12, seed)
        samples[gsm] = {**rec(title, label, [f["name"] for f in files.values()], gsm=gsm), "files": files}
    return samples


def never(*a, **k):
    raise AssertionError("a file was opened")


def test_load_refuses_a_sample_that_is_not_pinned_or_not_its_pin(tmp_path):
    """The first line at load: an unpinned sample (a window sample whose files do not exist, or a permitted one), a
    record changed by one character, a pin for a stage the record does not state, or pins covering one end only refuse
    before any file is opened."""
    samples = two_samples(tmp_path)
    cfg = clear_for(config(), samples)
    raw_manifest(tmp_path, cfg, {**samples, "GSM13": {**rec("E12.5 heart", "E12.5", gsm="GSM13"),
                                                      "files": {"matrix": {"name": "absent.mtx.gz", "sha256": "0"}}}})
    with pytest.raises(cm.Refused, match="not pinned"):
        ld.load_atlas(ACC, tmp_path, cfg, check_frozen=False, reader=never)
    raw_manifest(tmp_path, cfg, {**samples, "GSM14": {**rec("E14.5 heart", "E14.5", gsm="GSM14"), "files": {}}})
    with pytest.raises(cm.Refused, match="not pinned"):
        ld.load_atlas(ACC, tmp_path, cfg, check_frozen=False, reader=never)
    changed = copy.deepcopy(samples)
    changed["GSM12"]["record"][0][1] = ["E14.5 hearT"]  # one character
    raw_manifest(tmp_path, cfg, changed)
    with pytest.raises(cm.Refused, match="not its pin"):
        ld.load_atlas(ACC, tmp_path, cfg, check_frozen=False, reader=never)
    other = copy.deepcopy(cfg)
    other["atlases"][ACC]["late"]["labels"] = ["E14.5", "E16.5"]
    other["atlases"][ACC]["pinned"]["samples"][1]["stage"] = "E16.5"  # a stage the record does not state
    raw_manifest(tmp_path, other, samples)
    with pytest.raises(cm.Refused, match="not its pin"):
        ld.load_atlas(ACC, tmp_path, other, check_frozen=False, reader=never)
    one_end = copy.deepcopy(cfg)
    one_end["atlases"][ACC]["pinned"]["samples"] = one_end["atlases"][ACC]["pinned"]["samples"][1:]
    raw_manifest(tmp_path, one_end, samples)
    with pytest.raises(cm.Refused, match="both ends"):
        ld.load_atlas(ACC, tmp_path, one_end, check_frozen=False, reader=never)


def test_load_atlas_end_to_end_on_files(tmp_path):
    samples = two_samples(tmp_path)
    cfg = clear_for(config(), samples)
    raw_manifest(tmp_path, cfg, samples)
    red = ld.load_atlas(ACC, tmp_path, cfg, check_frozen=False)
    w = red["qc"]["window"]
    assert w["cells_in_window_after"] == 0 and "recomputed per cell" in w["cells_in_window_after_how"]
    assert w["samples_removed_window"] == []
    assert set(w["cells_per_stage_before"]) == {"E8.5", "E14.5"} and w["cells_per_stage_before"] == \
        w["cells_per_stage_after"]
    assert red["qc"]["columns_read"] == ["barcodes:barcode", "features:id", "features:symbol"]
    assert red["C"].shape[1] == 30  # every symbol is unique; the mitochondrial gene stays as a gene
    # the reduction goes outside the repository only
    with pytest.raises(cm.Refused):
        ld.write_reduced(ACC, ROOT / "data" / "derived" / "t1_ext_test", red, cfg)
    p = ld.write_reduced(ACC, tmp_path, red, cfg)
    back = ld.read_reduced(ACC, tmp_path, cfg)
    assert back["C"].shape == red["C"].shape and json.loads(p.read_text())["tasks"] == ["T1"]


def test_processing_clearance_refuses_the_load(tmp_path):
    samples = two_samples(tmp_path)
    pinned = pin_for(config(), samples)
    raw_manifest(tmp_path, pinned, samples)
    with pytest.raises(cm.Refused, match="not cleared"):
        ld.load_atlas(ACC, tmp_path, pinned, check_frozen=False, reader=never)
    other = copy.deepcopy(pinned)
    other["atlases"][ACC]["processing_cleared"] = {"sha256": "0" * 64}
    with pytest.raises(cm.Refused, match="now hashes"):
        ld.load_atlas(ACC, tmp_path, other, check_frozen=False, reader=never)


def test_read_reduced_rechecks_every_row(tmp_path):
    samples = two_samples(tmp_path)
    cfg = clear_for(config(), samples)
    raw_manifest(tmp_path, cfg, samples)
    red = ld.load_atlas(ACC, tmp_path, cfg, check_frozen=False)
    ld.write_reduced(ACC, tmp_path, red, cfg)
    out = tmp_path / "reduced" / ACC
    man = json.loads((out / "MANIFEST.json").read_text())

    def tamper(cells_text=None, records=None):
        if cells_text is not None:
            (out / "cells.csv").write_text(cells_text)
        m = copy.deepcopy(man)
        m["outputs"]["cells.csv"]["sha256"] = cm.sha256_file(out / "cells.csv")
        if records is not None:
            m["stage_records"] = records
        (out / "MANIFEST.json").write_text(json.dumps(m))

    good = (out / "cells.csv").read_text()
    lines = good.splitlines()
    # a row whose day is in the window, with the manifest's hashes made to match: the re-check refuses
    bad = "\n".join([lines[0], lines[1].rsplit(",", 1)[0] + ",12.5", *lines[2:]]) + "\n"
    tamper(bad)
    with pytest.raises(cm.Refused):
        ld.read_reduced(ACC, tmp_path, cfg)
    # a row whose day is outside the window but is not its sample's stage
    bad = "\n".join([lines[0], lines[1].rsplit(",", 1)[0] + ",16.5", *lines[2:]]) + "\n"
    tamper(bad)
    with pytest.raises(cm.Refused):
        ld.read_reduced(ACC, tmp_path, cfg)
    # a stored record that reads in the window
    recs = copy.deepcopy(man["stage_records"])
    recs["GSM12"]["record"].append(["description", ["E12.5 littermates"]])
    tamper(good, recs)
    with pytest.raises(cm.Refused):
        ld.read_reduced(ACC, tmp_path, cfg)
    tamper(good, man["stage_records"])
    assert ld.read_reduced(ACC, tmp_path, cfg)["C"].shape == red["C"].shape


def test_stageless_sample_refuses_the_load(tmp_path):
    cfg = config()
    samples = {"GSM21": {**rec("heart", "E14.5", gsm="GSM21"), "files": {}}}
    raw_manifest(tmp_path, cfg, samples)
    with pytest.raises(cm.Refused):
        ld.load_atlas(ACC, tmp_path, cfg, check_frozen=False, reader=never)


# ------------------------------------------------------------------------------------------------ dense tables


class Dense:
    def __init__(self, rows: dict, barcodes=("AAAC", "GGTT", "CCAA")):
        self.rows, self.bc = rows, barcodes

    def names(self):
        return ("gene", *self.bc)

    def column(self, name):
        assert name == "gene"
        return np.asarray(list(self.rows), dtype=object)

    def counts(self):
        return np.array(list(self.rows.values()), dtype=np.float64)


GENES = {"Myl7": [3, 0, 2], "Tnnt2": [1, 2, 0], "Epcam": [0, 1, 1], "Snn": [0, 0, 1], "mt-Co1": [1, 0, 0]}


@pytest.mark.parametrize("row", ["Cluster", "nUMI", "nGene", "nCount_RNA", "nFeature_RNA", "percent.mt", "cell_type",
                                 "louvain", "leiden_res1", "seurat_clusters", "RNA_snn_res.0.8", "batch", "Sample",
                                 "stage", "orig.ident", "CellType"])
def test_dense_table_refuses_annotation_and_summary_rows(row):
    with pytest.raises(cm.Refused):
        ld.dge_from_table(Dense({**GENES, row: [1, 2, 1]}), ld.Reader())


def test_dense_table_refuses_a_column_sum_or_gene_count_row():
    tot = np.sum(list(GENES.values()), axis=0).tolist()
    with pytest.raises(cm.Refused):
        ld.dge_from_table(Dense({**GENES, "Xist": tot}), ld.Reader())
    nz = np.sum(np.array(list(GENES.values())) > 0, axis=0).tolist()
    with pytest.raises(cm.Refused):
        ld.dge_from_table(Dense({**GENES, "Xist": nz}), ld.Reader())


def test_dense_table_keeps_genes_whose_names_hold_annotation_words():
    rd = ld.Reader()
    x = ld.dge_from_table(Dense(GENES), rd)  # Epcam holds 'pca'; Snn is the stannin gene
    assert list(x["genes"]) == list(GENES) and x["C"].shape == (3, 5)
    assert "checked" in x["gates"]["row_names"] and "checked" in x["gates"]["summary_rows"]
    assert any(r.startswith("dge:row names") for r in rd.read)


def test_dge_file_with_a_text_annotation_row_refuses_before_parsing(tmp_path):
    p = tmp_path / "GSM1_dge.txt.gz"
    p.write_bytes(gzip.compress(b'"AAAC"\t"GGTT"\nMyl7\t3\t0\nCluster\tCM\tEndo\n'))
    with pytest.raises(cm.Refused, match="summaries or annotations"):
        ld.dge_from_table(ld.DgeFile(p), ld.Reader())


def test_tsv_table_reads_only_its_column(tmp_path):
    p = tmp_path / "GSM1_barcodes.tsv.gz"
    p.write_bytes(gzip.compress(b"AAAC-1\tauthor\tcluster 3\nGGTT-1\tauthor\n"))
    t = ld.TsvTable(p, ld.WHITELIST["barcodes"])
    assert t.names() == ("barcode",) and t.column("barcode").tolist() == ["AAAC-1", "GGTT-1"]


def test_draft_atlases_yaml_refuses_load():
    cfg = copy.deepcopy(cm.load_atlases())
    cfg["status"] = "draft"  # whatever the committed file says (the placement commit freezes it)
    with pytest.raises(cm.Refused):
        cm.require_frozen(cfg)


def test_cell_call_rule():
    C = sparse.csr_matrix(np.array([[10, 0, 0], [3, 3, 0], [1, 1, 8], [0, 0, 0]]))
    keep = ld.cell_call(C, np.array(["A", "B", "mt-Nd1"]), {"min_umi": 5, "min_genes": 2, "max_mito": 0.5})
    assert keep.tolist() == [False, True, False, False]
