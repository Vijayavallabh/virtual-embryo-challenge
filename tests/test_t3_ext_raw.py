"""Stage 7, piece 2 of experiments/t3_ext (campaign/t3-ext-3): the intake for per-sample UNFILTERED matrices (a
dataset's raw_matrices declaration, the one fixed cell call of t3x_common, the sparse two-pass reader). Synthetic SOFT
text, synthetic counts and a stub download only: nothing here touches the network or real data."""

from __future__ import annotations

import gzip
import importlib.util
import inspect
import io
import json
import re
import shutil
import sys
import tarfile
from pathlib import Path

import numpy as np
import pytest
from scipy import io as spio
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(f"t3_ext_raw_{name}", EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


fe, pr = _load("fetch"), _load("process")

ACC = "GSE920001"
HOST = "ftp://ftp.ncbi.nlm.nih.gov/geo"
PIN = "    pinned: []  # GSE920001\n"
RAW_DECL = """    raw_matrices:
      patterns: ['^GSM\\d+_[A-Za-z0-9_.+-]*raw_feature_bc_matrix\\.tar\\.gz$']
      reason: "the samples' only per-sample files are unfiltered matrices (synthetic)"
"""
DS = """version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
  - accession: GSE920001
    status: cleared
    tier: unrelated
    study: GSE920001
    design: germline
    control_relation: same_series
    repository: geo
    screen_series: [GSE920001]
    stages_used: [E8.5]
    sample_file_patterns: ['raw_feature_bc_matrix\\.tar\\.gz$', '(matrix\\.mtx|barcodes\\.tsv|features\\.tsv)(\\.gz)?$']
RAW
""" + PIN


def _yaml(tmp_path, declared=True) -> Path:
    p = tmp_path / "datasets.yaml"
    p.write_text(DS.replace("RAW\n", RAW_DECL if declared else ""))
    return p


def raw_name(gsm):
    return f"{gsm}_E8_5_lib_raw_feature_bc_matrix.tar.gz"


def _sample(gsm, title, geno, files):
    stub = "GSM" + gsm[3:-3] + "nnn"
    rec = {"title": [title], "source_name_ch1": ["whole embryo"],
           "characteristics_ch1": [f"genotype: {geno}", "developmental stage: E8.5"],
           "extract_protocol_ch1": ["Chromium Single Cell 3' v3"], "series_id": [ACC]}
    for i, f in enumerate(files):
        rec[f"supplementary_file_{i + 1}"] = [f"{HOST}/samples/{stub}/{gsm}/suppl/{f}"]
    return rec


SAMPLES = {"GSM9301": _sample("GSM9301", "Zzz KO E8.5 1", "Zzz -/-", [raw_name("GSM9301")]),
           "GSM9302": _sample("GSM9302", "WT E8.5 1", "wild type", [raw_name("GSM9302")]),
           "GSM9303": _sample("GSM9303", "Zzz KO E8.5 2", "Zzz -/-", [raw_name("GSM9303")]),
           # a sample with both forms: the filtered triple is used, its raw archive refused
           "GSM9304": _sample("GSM9304", "WT E8.5 2", "wild type",
                              [raw_name("GSM9304"), "GSM9304_matrix.mtx.gz", "GSM9304_barcodes.tsv.gz",
                               "GSM9304_features.tsv.gz"])}


def _soft(samples=None) -> str:
    lines = [f"^SERIES = {ACC}", f"!Series_title = A synthetic series {ACC}", "!Series_summary = embryos",
             "!Series_overall_design = one embryo per library", "!Series_contact_name = A,,Person"]
    for gsm, rec in (samples or SAMPLES).items():
        lines.append(f"^SAMPLE = {gsm}")
        for k, vals in rec.items():
            for v in (vals if isinstance(vals, list) else [vals]):
                lines.append(f"!Sample_{k} = {v}")
    return "\n".join(lines) + "\n"


def _plan(y, dest, samples=None):
    return fe.plan(ACC, dest, y, metas={ACC: fe.parse_soft(_soft(samples))})


# ---------------------------------------------------------------------------------------------------------------------
# a synthetic unfiltered library with an empty-droplet tail

PANEL = [f"g{i}" for i in range(40)]
EXTRA = 260
GROUPS = ["CM", "Endoderm", "Neural"]


def _library(rng, mu, fold=None, n_cells=600, border=40, n_empty=5000, n_ab=10):
    """(features x barcodes counts, barcode names, gene-expression totals). Cells: totals 6000 + 10 i (half CM, half
    Neural); `border` barcodes at 1170 and `border` at 1169; an empty-droplet tail of totals 0 to 20; n_ab barcodes
    with antibody counts only (no gene-expression UMI). The last feature is an Antibody Capture feature."""
    fold = np.ones(len(PANEL)) if fold is None else fold
    n_gex = len(PANEL) + EXTRA
    cols, names, totals = [], [], []

    def add(total, p, name):
        c = np.zeros(n_gex + 1, np.int64)
        c[:n_gex] = rng.multinomial(int(total), p) if total else 0
        cols.append(c)
        names.append(name)
        totals.append(int(total))

    uniform = np.full(n_gex, 1.0 / n_gex)
    for i in range(n_cells):
        g = "CM" if i % 2 == 0 else "Neural"
        q = mu[g] * fold
        p = np.concatenate([0.5 * q / q.sum(), np.full(EXTRA, 0.5 / EXTRA)])
        add(6000 + 10 * i, p, f"{g}_{i}-1")
    for i in range(border):
        add(1170, uniform, f"L_{i}-1")
    for i in range(border):
        add(1169, uniform, f"N_{i}-1")
    for i, t in enumerate(rng.integers(0, 21, size=n_empty)):
        add(t, uniform, f"E_{i}-1")
    for i in range(n_ab):
        c = np.zeros(n_gex + 1, np.int64)
        c[-1] = 5000
        cols.append(c)
        names.append(f"A_{i}-1")
        totals.append(0)
    perm = rng.permutation(len(names))   # cells are not first in the file
    X = np.vstack(cols).T[:, perm]
    return X, [names[j] for j in perm], np.array(totals)[perm]


def _triple_bytes(X, names):
    b = io.BytesIO()
    spio.mmwrite(b, sparse.coo_matrix(X))
    feats = [f"E{i}\t{s}\tGene Expression" for i, s in enumerate(PANEL + [f"x{j}" for j in range(EXTRA)])]
    feats.append(f"AB{0}\tab0\tAntibody Capture")
    return {"matrix.mtx.gz": gzip.compress(b.getvalue()), "features.tsv.gz": gzip.compress("\n".join(feats).encode()),
            "barcodes.tsv.gz": gzip.compress("\n".join(names).encode())}


def _tar(path: Path, members: dict[str, bytes]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as t:
        for name, data in members.items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            t.addfile(ti, io.BytesIO(data))
    return path


def _raw_archive(path, X, names, folder="raw_feature_bc_matrix"):
    return _tar(path, {f"{folder}/{k}": v for k, v in _triple_bytes(X, names).items()})


def _mu(rng):
    return {g: rng.dirichlet(np.full(len(PANEL), 0.5)) for g in GROUPS}


def _hand_call(totals):
    """The fixed cell call computed by hand: the 3000 largest totals in ascending order, the 99th percentile by
    linear interpolation between the two order statistics around position 0.99 (n - 1)."""
    top = sorted(totals, reverse=True)[:3000][::-1]
    pos = 0.99 * (len(top) - 1)
    lo = int(pos)
    m = top[lo] + (pos - lo) * (top[lo + 1] - top[lo])
    cutoff = max(500.0, m / 10.0)
    return m, cutoff, [t >= cutoff for t in totals]


# ---------------------------------------------------------------------------------------------------------------------
# the fixed cell call


def test_the_cell_call_equals_a_hand_computation_on_a_library_with_an_empty_droplet_tail(tmp_path):
    rng = np.random.default_rng(1)
    X, names, totals = _library(rng, _mu(rng))
    # by hand: ascending, the 3000 largest are 2320 empties, 40 at 1169, 40 at 1170, then the 600 cells (index 2400
    # on); position 0.99 x 2999 = 2969.01 falls between cell 569 (11690) and cell 570 (11700)
    m_hand = 11690 + (0.99 * 2999 - 2969) * (11700 - 11690)
    assert m_hand == pytest.approx(11690.1)
    cutoff_hand = max(500, m_hand / 10)               # 1169.01: 1170 is a cell, 1169 is not
    called, rec = cm.raw_cell_call(totals)
    assert rec["m"] == pytest.approx(m_hand) and rec["cutoff"] == pytest.approx(cutoff_hand)
    assert rec["cells_called"] == 640 and called.sum() == 640
    assert _hand_call(list(totals))[0] == pytest.approx(m_hand)
    want = sorted(n for n, t in zip(names, totals, strict=True) if t >= cutoff_hand)
    assert len(want) == 640 and not any(n.startswith(("N_", "E_", "A_")) for n in want)
    # the reader gives the same call on the archive, keeps only the called columns, and records the call
    p = _raw_archive(tmp_path / raw_name("GSM1"), X, names)
    M, symbols, types, bcs, r = pr.read_raw_sample([p], "GSM1")
    assert sorted(bcs) == want and M.shape == (640, len(PANEL) + EXTRA + 1)
    keep = [names.index(b) for b in bcs]
    assert (M.toarray() == X[:, keep].T).all()          # the called cells' counts, unchanged
    any_count = int((X.sum(0) > 0).sum())
    gex_tot = X[:-1].sum(0)
    assert r == {"barcodes_in_file": len(names), "barcodes_with_any_count": any_count, "m": pytest.approx(m_hand),
                 "cutoff": pytest.approx(cutoff_hand), "cells_called": 640,
                 "umi_share_called": pytest.approx(gex_tot[keep].sum() / gex_tot.sum()),
                 "features": len(PANEL) + EXTRA + 1, "entries": int((X > 0).sum()),
                 # stage 8 records, by hand: the called totals are 40 at 1170 and the 600 cells at 6000 + 10 i;
                 # positions 0.05, 0.25, 0.5 x 639 fall at 1170, cell 119.75 and cell 279.5; the 40 border barcodes
                 # are under twice the cutoff (2338.02)
                 "called_total_p5": 1170.0, "called_total_p25": pytest.approx(7197.5),
                 "called_total_p50": pytest.approx(8795.0), "called_under_twice_cutoff": 40}
    # a loose raw triple reads the same
    for k, v in _triple_bytes(X, names).items():
        (tmp_path / f"GSM2_raw_feature_bc_matrix_{k}").write_bytes(v)
    loose = sorted(tmp_path.glob("GSM2_*"))
    M2, _, _, bcs2, r2 = pr.read_raw_sample(loose, "GSM2")
    assert bcs2 == bcs and (M2 != M).nnz == 0 and r2 == r
    # nothing was unpacked
    assert sorted(x.name for x in tmp_path.iterdir() if x.is_dir()) == []


def test_the_cell_call_is_fixed_not_configurable():
    assert (cm.RAW_CALL_N, cm.RAW_CALL_PERCENTILE, cm.RAW_CALL_RATIO, cm.RAW_CALL_MIN_UMI, cm.RAW_MIN_CELLS) == (
        3000, 99.0, 10.0, 500, 500)
    assert list(inspect.signature(cm.raw_cell_call).parameters) == ["totals"]
    assert list(inspect.signature(pr.read_raw_sample).parameters) == ["paths", "label"]
    for f in ("t3x_common.py", "process.py", "fetch.py"):
        src = (EXT / f).read_text()
        assert "os.environ" not in src and "getenv" not in src and "--raw" not in src and "--cell" not in src, f
    # fewer than N barcodes: the percentile of all of them
    called, rec = cm.raw_cell_call([10000, 2000, 900, 400, 0])
    want = float(np.percentile([10000, 2000, 900, 400, 0], 99))
    assert rec["m"] == pytest.approx(want) and rec["cutoff"] == pytest.approx(want / 10)
    assert called.tolist() == [True, True, False, False, False]
    # the floor of 500 UMIs holds when m is small
    assert cm.raw_cell_call([3000, 600, 499])[1]["cutoff"] == 500.0


def test_a_library_outside_the_limits_is_refused(tmp_path, monkeypatch):
    rng = np.random.default_rng(2)
    mu = _mu(rng)
    X, names, _ = _library(rng, mu, n_cells=200, border=0, n_empty=300, n_ab=0)
    few = _raw_archive(tmp_path / "few" / raw_name("GSM1"), X, names)
    with pytest.raises(cm.Refused, match="keeps 200 cells, fewer than 500"):
        pr.read_raw_sample([few], "GSM1")
    X, names, _ = _library(rng, mu, n_empty=300, n_ab=0)
    ok = _raw_archive(tmp_path / "ok" / raw_name("GSM1"), X, names)
    # (with a short empty tail m is higher: the 600 cells are called, the 1170 border is not)
    for const, value, match in (("MAX_BARCODES", 599, "more than the cap 599"),
                                ("RAW_MAX_ENTRIES", 1000, "more than 1000"),
                                ("RAW_MAX_BARCODES", 100, "more than 100: refused"),
                                ("RAW_MAX_FILE_BYTES", 1000, "larger than 1000 bytes")):
        with monkeypatch.context() as mp:
            mp.setattr(cm, const, value)
            with pytest.raises(cm.Refused, match=match):
                pr.read_raw_sample([ok], "GSM1")
    h5 = tmp_path / "GSM1_raw_feature_bc_matrix.h5"
    h5.write_bytes(b"not read")
    with pytest.raises(cm.Refused, match="not an h5"):
        pr.read_raw_sample([h5], "GSM1")
    # a matrix whose header states more entries than it holds
    m = _triple_bytes(X, names)
    body = gzip.decompress(m["matrix.mtx.gz"]).decode().split("\n")
    i = next(k for k, line in enumerate(body) if not line.startswith("%"))
    f, b, n = body[i].split()
    body[i] = f"{f} {b} {int(n) + 5}"
    short = _tar(tmp_path / "short" / raw_name("GSM1"),
                 {**{f"raw_feature_bc_matrix/{k}": v for k, v in m.items()},
                  "raw_feature_bc_matrix/matrix.mtx.gz": gzip.compress("\n".join(body).encode())})
    with pytest.raises(cm.Refused, match="not the .* its header states"):
        pr.read_raw_sample([short], "GSM1")


@pytest.mark.parametrize("kind,match", [
    ("traversal", "path traversal"), ("symlink", "symlink"), ("extra", "no rule"), ("two_folders", "more than one"),
    ("fragments", "denied type"), ("absolute", "absolute"), ("screened", "held-out screen")])
def test_a_crafted_member_of_a_raw_archive_is_refused(tmp_path, kind, match):
    rng = np.random.default_rng(3)
    X, names, _ = _library(rng, _mu(rng), n_empty=50, n_ab=0)
    good = {f"raw_feature_bc_matrix/{k}": v for k, v in _triple_bytes(X, names).items()}
    extra = {"traversal": "raw_feature_bc_matrix/../../evil/matrix.mtx.gz", "extra": "raw_feature_bc_matrix/notes.txt",
             "two_folders": "other/matrix.mtx.gz", "fragments": "raw_feature_bc_matrix/atac_fragments.tsv.gz",
             "absolute": "/evil/barcodes.tsv.gz", "screened": "raw_feature_bc_matrix/Gata4_barcodes.tsv.gz",
             "symlink": "raw_feature_bc_matrix/link"}[kind]
    p = tmp_path / raw_name("GSM1")
    with tarfile.open(p, "w:gz") as t:
        for name, data in good.items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            t.addfile(ti, io.BytesIO(data))
        ti = tarfile.TarInfo(extra)
        if kind == "symlink":
            ti.type, ti.linkname = tarfile.SYMTYPE, "/etc/passwd"
            t.addfile(ti)
        else:
            ti.size = 1
            t.addfile(ti, io.BytesIO(b"x"))
    before = sorted(tmp_path.rglob("*"))
    with pytest.raises(cm.Refused, match=match):
        pr.read_raw_sample([p], "GSM1")
    assert sorted(tmp_path.rglob("*")) == before


# ---------------------------------------------------------------------------------------------------------------------
# the declaration


def test_without_the_declaration_a_raw_file_is_still_refused(tmp_path):
    rec = SAMPLES["GSM9301"]
    ds_no = cm.dataset(ACC, _yaml(tmp_path, declared=False))
    ok, other = fe.sample_files(ds_no, rec, "GSM9301")
    assert ok == [] and other == [{"name": raw_name("GSM9301"), "reason": "a denied file type"}]
    pl = _plan(_yaml(tmp_path, declared=False), tmp_path / "d")
    assert {r["accession"]: r["files"] for r in pl["tables"][ACC]}["GSM9301"] == []
    with pytest.raises(cm.Refused, match="no per-sample file"):
        fe.pin_template(ACC, tmp_path / "d", "ko", "Zzz", ["GSM9301"], datasets=_yaml(tmp_path, declared=False))
    # the filtered routes never read an unfiltered archive, declared or not
    rng = np.random.default_rng(4)
    X, names, _ = _library(rng, _mu(rng), n_empty=50, n_ab=0)
    p = _raw_archive(tmp_path / "a" / raw_name("GSM9301"), X, names)
    with pytest.raises(cm.Refused, match="denied type"):
        pr.read_sample([p])
    # declared: offered and marked raw; a name the declared pattern misses stays denied; every other term stays
    ds = cm.dataset(ACC, _yaml(tmp_path))
    ok, _ = fe.sample_files(ds, rec, "GSM9301")
    assert ok == [{"name": raw_name("GSM9301"), "url": ok[0]["url"], "raw": True}]
    for name in ("GSM9301_raw_feature_bc_matrix.h5", "GSM9301_RAW_FEATURE_BC_MATRIX.tar.gz",
                 "GSM9301_atac_raw_feature_bc_matrix.tar.gz", "GSM9301_fragments_raw_feature_bc_matrix.tar.gz"):
        r2 = _sample("GSM9301", "x", "y", [name])
        ok2, other2 = fe.sample_files(ds, r2, "GSM9301")
        assert ok2 == [] and other2[0]["reason"] in ("a denied file type",
                                                     "not matched by the dataset's sample_file_patterns"), name
    assert fe.denied("x_raw_feature_bc_matrix.tar.gz") and not fe.denied("x_raw_feature_bc_matrix.tar.gz", True)
    assert fe.denied("x_fragments.tsv.gz", True) and fe.denied("GSE1_RAW.tar", True)
    # the declaration's form
    y = _yaml(tmp_path)
    for old, new, match in (("reason: \"the samples'", "reason: \"TO FILL the samples'", "reason"),
                            ("patterns: ['^GSM", "patterns: ['(^GSM", "does not compile"),
                            ("    raw_matrices:\n", "    raw_matrices:\n      cutoff: 300\n", "exactly the keys")):
        y.write_text(DS.replace("RAW\n", RAW_DECL.replace(old, new)))
        with pytest.raises(cm.Refused, match=match):
            cm.load_datasets(y)


def test_a_sample_with_both_forms_uses_the_filtered_one(tmp_path):
    y = _yaml(tmp_path)
    ds = cm.dataset(ACC, y)
    ok, other = fe.sample_files(ds, SAMPLES["GSM9304"], "GSM9304")
    assert [f["name"] for f in ok] == ["GSM9304_matrix.mtx.gz", "GSM9304_barcodes.tsv.gz", "GSM9304_features.tsv.gz"]
    assert not any(f.get("raw") for f in ok)
    assert other == [{"name": raw_name("GSM9304"), "reason": "an unfiltered matrix of a sample that lists other files "
                      "of its own (GSM9304_matrix.mtx.gz, GSM9304_barcodes.tsv.gz, GSM9304_features.tsv.gz): taken "
                      "only when it is the sample's only file"}]
    pl = _plan(y, tmp_path / "d")
    text = (tmp_path / "d" / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert "# RAW MATRICES (datasets.yaml raw_matrices)" in text and "the one fixed cell call" in text
    assert f"UNFILTERED (raw_matrices; read through the fixed cell call; confirm it is ONE library of THIS sample): " \
           f"['{raw_name('GSM9301')}']" in text
    rows = {r["accession"]: r for r in pl["tables"][ACC]}
    assert [f["name"] for f in rows["GSM9304"]["files"]] == [f["name"] for f in ok]
    txt = fe.pin_template(ACC, tmp_path / "d", "control", "Zzz", ["GSM9304"], datasets=y)
    files_line = [x for x in txt.splitlines() if x.strip().startswith("files:")]
    assert len(files_line) == 1 and raw_name("GSM9304") not in files_line[0] and "GSM9304_matrix" in files_line[0]
    assert "UNFILTERED" not in txt
    assert "UNFILTERED" in fe.pin_template(ACC, tmp_path / "d", "ko", "Zzz", ["GSM9301"], datasets=y)


# ---------------------------------------------------------------------------------------------------------------------
# fetch and process


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


def _setup(tmp_path, monkeypatch):
    rng = np.random.default_rng(7)
    mu = _mu(rng)
    wt = tmp_path / "wt.h5ad"
    wt.write_bytes(b"stand-in for the WT reference")
    refp = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats({"genes": PANEL, "groups": GROUPS, "shares": {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5},
                       "mu": mu, "count": {"CM": 20, "Endoderm": 30, "Neural": 50}}, refp, cm.sha256_file(wt))
    y = _yaml(tmp_path)
    dest = tmp_path / "dest"
    assert _plan(y, dest)["status"] == "ok"
    block = (fe.pin_template(ACC, dest, "ko", "Zzz", ["GSM9301"], datasets=y)
             + fe.pin_template(ACC, dest, "control", "Zzz", ["GSM9302"], datasets=y))
    y.write_text(y.read_text().replace(PIN, "    pinned_series:\n" + fe.series_pin_template(ACC, dest, datasets=y)
                                       + "    pinned:\n" + block))
    assert _plan(y, dest)["status"] == "ok"
    staging = tmp_path / "staging"
    fold = np.ones(len(PANEL))
    fold[:4] = 3.0
    totals = {}
    for gsm, f in (("GSM9301", fold), ("GSM9302", None), ("GSM9303", fold), ("GSM9304", None)):
        X, names, totals[gsm] = _library(rng, mu, fold=f, n_empty=3000, n_ab=5)
        _raw_archive(staging / raw_name(gsm), X, names)
    calls: list = []
    monkeypatch.setattr(fe, "_get", _fake_get(staging, calls))
    return y, dest, wt, refp, calls, totals


def test_an_unpinned_samples_raw_file_is_never_fetched_or_read(tmp_path, monkeypatch):
    y, dest, wt, refp, calls, _ = _setup(tmp_path, monkeypatch)
    man = fe.fetch(ACC, dest, y)
    assert sorted(Path(u).name for u in calls) == [raw_name("GSM9301"), raw_name("GSM9302")]
    assert set(man["files"]) == {f"GSM9301/{raw_name('GSM9301')}", f"GSM9302/{raw_name('GSM9302')}"}
    assert not (dest / "raw" / ACC / "GSM9303").exists() and not (dest / "raw" / ACC / "GSM9304").exists()
    read: list = []
    real = pr.read_raw_sample

    def spy(paths, label):
        read.append((label, sorted(Path(p).name for p in paths)))
        return real(paths, label)

    monkeypatch.setattr(pr, "read_raw_sample", spy)
    monkeypatch.setattr(pr, "read_sample", lambda paths: pytest.fail("the filtered route read a raw file"))
    pr.run(ACC, dest, wt, refp, tmp_path / "processed", y, classifier=_clf)
    assert read == [("GSM9301", [raw_name("GSM9301")]), ("GSM9302", [raw_name("GSM9302")])]
    # an unpinned sample's archive placed by hand in the raw folder is never read (only pinned folders are opened)
    (dest / "raw" / ACC / "GSM9303").mkdir()
    shutil.copyfile(tmp_path / "staging" / raw_name("GSM9303"), dest / "raw" / ACC / "GSM9303" / raw_name("GSM9303"))
    read.clear()
    pr.run(ACC, dest, wt, refp, tmp_path / "processed", y, classifier=_clf)
    assert [r[0] for r in read] == ["GSM9301", "GSM9302"]


def test_the_same_cell_call_is_applied_to_both_arms_and_recorded(tmp_path, monkeypatch):
    y, dest, wt, refp, calls, totals = _setup(tmp_path, monkeypatch)
    fe.fetch(ACC, dest, y)
    qc = pr.run(ACC, dest, wt, refp, tmp_path / "processed", y, classifier=_clf)
    assert qc["status"] == "ok", qc
    rc = qc["raw_cell_call"]
    assert rc["rule"] == cm.raw_call_rule() and rc["raw_matrices"] == cm.dataset(ACC, y)["raw_matrices"]
    assert set(rc["samples"]) == {"GSM9301", "GSM9302"}           # the knockout and its control
    for gsm, rec in rc["samples"].items():
        m, cutoff, called = _hand_call(list(totals[gsm]))
        assert rec["m"] == pytest.approx(m) and rec["cutoff"] == pytest.approx(cutoff)
        assert rec["cells_called"] == sum(called) == 640
        assert set(rec) == {"barcodes_in_file", "barcodes_with_any_count", "m", "cutoff", "cells_called",
                            "umi_share_called", "features", "entries", "called_total_p5", "called_total_p25",
                            "called_total_p50", "called_under_twice_cutoff"}
        # stage 8 records (the second reviewer of ec86a1b, lesser point 1), by hand
        ct = np.array([t for t, c in zip(totals[gsm], called, strict=True) if c], float)
        for q in (5, 25, 50):
            assert rec[f"called_total_p{q}"] == pytest.approx(float(np.percentile(ct, q)))
        assert rec["called_under_twice_cutoff"] == int((ct < 2 * cutoff).sum())
    u = qc["units"]["Zzz"]
    assert u["status"] == "ok" and u["raw_called"] is True and u["raw_called_samples"] == ["GSM9301", "GSM9302"]
    # stage 8: the per-arm medians of the call and their knockout-to-control ratios, a record only
    arms = u["raw_call_arms"]
    for k in ("cutoff", "cells_called", "umi_share_called"):
        assert arms["median_ko"][k] == pytest.approx(rc["samples"]["GSM9301"][k])
        assert arms["median_control"][k] == pytest.approx(rc["samples"]["GSM9302"][k])
        assert arms["ko_to_control"][k] == pytest.approx(rc["samples"]["GSM9301"][k] / rc["samples"]["GSM9302"][k])
    assert "record only" in arms["note"]
    # downstream the called library is an ordinary one: the per-cell rule ran on the called cells only, the unit
    # file has no raw field, and the pins hash binds the declaration
    assert set(qc["cells_dropped_low_genes"]) == {"GSM9301", "GSM9302"}
    assert sum(u["cells_per_stage_and_arm"].values()) <= 2 * 640
    unit = json.loads(Path(u["unit_file"]).read_text())
    assert not any("raw" in k for k in unit)
    assert qc["pins_sha256"] == cm.pins_sha256(cm.dataset(ACC, y))
    ds = cm.dataset(ACC, y)
    assert cm.pins_sha256({k: v for k, v in ds.items() if k != "raw_matrices"}) != qc["pins_sha256"]


def test_a_plan_made_under_another_raw_declaration_is_refused(tmp_path, monkeypatch):
    y, dest, wt, refp, calls, _ = _setup(tmp_path, monkeypatch)
    ds = cm.dataset(ACC, y)
    pl = json.loads((dest / "meta" / ACC / "plan.json").read_text())
    assert pl["raw_matrices"] == ds["raw_matrices"]
    with pytest.raises(cm.Refused, match="another raw_matrices declaration"):
        fe.verify_pins({**ds, "raw_matrices": {**ds["raw_matrices"], "reason": "other"}}, pl)
    # a raw entry cannot be pinned beside a filtered file of the same sample, and a hand-edited stored list refuses
    bad = json.loads((dest / "meta" / ACC / "plan.json").read_text())
    for r in bad["tables"][ACC]:
        if r["accession"] == "GSM9301":
            r["files"][0].pop("raw")
    with pytest.raises(cm.Refused, match="stored file list differs"):
        fe.verify_pins(ds, bad)


# ---------------------------------------------------------------------------------------------------------------------
# the committed file, the README and the launcher


def test_the_committed_declaration_is_on_gse270969_only_with_no_pins():
    # The nine samples were pinned on 2026-11-20 and withdrawn the same day, before any table (docs/PLAN.md decision
    # log 2026-11-20 (bbb), (ddd)); the declaration stays on this dataset alone.
    allds = cm.load_datasets()
    assert [a for a, d in allds.items() if d.get("raw_matrices")] == ["GSE270969"]
    ds = allds["GSE270969"]
    assert ds["pinned"] == [] and ds["pinned_series"] == []
    host = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9999nnn/GSM9999001/suppl/"
    ok, _ = fe.sample_files(ds, {"supplementary_file_1": [host + "GSM9999001_KO_1_raw_feature_bc_matrix.tar.gz"]},
                            "GSM9999001")
    assert [f.get("raw") for f in ok] == [True]
    # another dataset still refuses the same name
    ok, other = fe.sample_files(allds["GSE277839"], {"supplementary_file_1": [
        host + "GSM9999001_KO_1_raw_feature_bc_matrix.tar.gz"]}, "GSM9999001")
    assert ok == [] and other[0]["reason"] == "a denied file type"


def test_the_readme_and_launcher_describe_the_unfiltered_intake():
    text = (EXT / "README.md").read_text()
    for s in ("## Unfiltered per-sample matrices", "raw_matrices", "max(500, m / 10)", "ONE library of ONE sample",
              "barcodes with any count"):
        assert s in text, s
    assert "raw_matrices" in (EXT / "run_dgx.sh").read_text()


# ---------------------------------------------------------------------------------------------------------------------
# stage 8 (campaign/t3-ext-4): the blockers of the first review of ec86a1b, each replaying the reviewer's reproduction,
# and the lesser points of both reviews


FILTERED_NAMES = ("GSM9301_E8_5_lib_filtered_feature_bc_matrix.h5",
                  "GSM9301_E8_5_lib_filtered_feature_bc_matrix.tar.gz", "GSM9301_filtered_matrix.tar.gz")


@pytest.mark.parametrize("filtered", FILTERED_NAMES)
def test_stage8_blocker_1_a_raw_file_beside_any_other_file_of_its_own_is_never_fetched(tmp_path, filtered):
    """The reviewer's reproduction: with GSE270969's own two patterns and its raw declaration, a sample listing its raw
    archive together with a filtered h5, a filtered archive or 'GSM9301_filtered_matrix.tar.gz' came back with ok =
    [the raw archive]: the filtered file, which the patterns do not offer, was put under 'other' first. Now the raw
    file is offered only when it is the sample's only file of its own."""
    rec = _sample("GSM9301", "Zzz KO E8.5 1", "Zzz -/-", [raw_name("GSM9301"), filtered])
    ds = cm.load_datasets()["GSE270969"]          # the committed record: its own patterns and declaration
    ok, other = fe.sample_files(ds, rec, "GSM9301")
    assert ok == []
    assert {o["name"]: o["reason"] for o in other}[raw_name("GSM9301")] == (
        f"an unfiltered matrix of a sample that lists other files of its own ({filtered}): taken only when it is the "
        "sample's only file")
    # the plan table lists the raw file as never fetched, with that reason, and nothing of the sample can be pinned
    y = _yaml(tmp_path)
    pl = _plan(y, tmp_path / "d", samples={**SAMPLES, "GSM9301": rec})
    assert pl["status"] == "ok"
    row = {r["accession"]: r for r in pl["tables"][ACC]}["GSM9301"]
    assert row["files"] == [] and raw_name("GSM9301") in {o["name"] for o in row["other_files"]}
    text = (tmp_path / "d" / "meta" / ACC / f"plan_table_{ACC}.txt").read_text()
    assert f"{raw_name('GSM9301')} (an unfiltered matrix of a sample that lists other files of its own" in text
    with pytest.raises(cm.Refused, match="no per-sample file"):
        fe.pin_template(ACC, tmp_path / "d", "ko", "Zzz", ["GSM9301"], datasets=y)


def test_stage8_blocker_1_two_raw_files_or_any_link_of_its_own_refuse_and_a_lone_raw_file_is_offered():
    ds = cm.load_datasets()["GSE270969"]
    two = _sample("GSM9301", "x", "y", [raw_name("GSM9301"), "GSM9301_second_raw_feature_bc_matrix.tar.gz"])
    ok, other = fe.sample_files(ds, two, "GSM9301")
    assert ok == [] and len(other) == 2 and all("other files of its own" in o["reason"] for o in other)
    # a file of its own whose name breaks the plain-name rule still counts as a file of its own
    odd = _sample("GSM9301", "x", "y", [raw_name("GSM9301"), "GSM9301_notes file.txt"])
    ok, other = fe.sample_files(ds, odd, "GSM9301")
    assert ok == [] and any("other files of its own" in o["reason"] for o in other)
    # a link into another folder is not a file of its own: the raw file stays offered (that link is never fetched)
    elsewhere = _sample("GSM9301", "x", "y", [raw_name("GSM9301")])
    elsewhere["supplementary_file_2"] = [f"{HOST}/series/GSE920nnn/{ACC}/suppl/{ACC}_RAW.tar"]
    ok, other = fe.sample_files(ds, elsewhere, "GSM9301")
    assert [f["name"] for f in ok] == [raw_name("GSM9301")] and ok[0]["raw"] is True
    assert other == [{"name": f"{ACC}_RAW.tar", "reason": "url not under the sample's own GEO folder"}]
    alone = _sample("GSM9301", "x", "y", [raw_name("GSM9301")])
    assert [f.get("raw") for f in fe.sample_files(ds, alone, "GSM9301")[0]] == [True]


def _mixed_setup(tmp_path, monkeypatch):
    """The reviewer's synthetic study: knockout GSM9301 pinned to its raw archive, control GSM9304 (which lists both
    forms) on its filtered triple."""
    rng = np.random.default_rng(17)
    mu = _mu(rng)
    wt = tmp_path / "wt.h5ad"
    wt.write_bytes(b"stand-in for the WT reference")
    refp = tmp_path / "ref" / "refstats.npz"
    tb.save_ref_stats({"genes": PANEL, "groups": GROUPS, "shares": {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5},
                       "mu": mu, "count": {"CM": 20, "Endoderm": 30, "Neural": 50}}, refp, cm.sha256_file(wt))
    y = _yaml(tmp_path)
    dest = tmp_path / "dest"
    assert _plan(y, dest)["status"] == "ok"
    before = json.loads((dest / "meta" / ACC / "plan.json").read_text())
    block = (fe.pin_template(ACC, dest, "ko", "Zzz", ["GSM9301"], datasets=y)
             + fe.pin_template(ACC, dest, "control", "Zzz", ["GSM9304"], datasets=y))
    y.write_text(y.read_text().replace(PIN, "    pinned_series:\n" + fe.series_pin_template(ACC, dest, datasets=y)
                                       + "    pinned:\n" + block))
    staging = tmp_path / "staging"
    fold = np.ones(len(PANEL))
    fold[:4] = 3.0
    X, names, _ = _library(rng, mu, fold=fold, n_empty=3000, n_ab=5)
    _raw_archive(staging / raw_name("GSM9301"), X, names)
    X, names, _ = _library(rng, mu, n_empty=0, n_ab=0, border=0)
    for k, v in _triple_bytes(X, names).items():
        (staging / f"GSM9304_{k}").write_bytes(v)
    calls: list = []
    monkeypatch.setattr(fe, "_get", _fake_get(staging, calls))
    return y, dest, wt, refp, calls, before


def test_stage8_blocker_2_a_unit_whose_arms_are_read_differently_is_refused(tmp_path, monkeypatch):
    """The reviewer's reproduction: knockout GSM9301 pinned to its raw archive, control GSM9304 moved to its filtered
    triple. The YAML loaded, the plan was ok, fetch fetched, and process returned unit Zzz ok with raw_called_samples
    ['GSM9301']. Now verify_pins (so plan with pins, fetch, process and combine) refuses the unit, and process
    refuses it again as a second guard."""
    y, dest, wt, refp, calls, before = _mixed_setup(tmp_path, monkeypatch)
    ds = cm.dataset(ACC, y)                          # the YAML still loads: the check needs the plan's file lists
    assert [p["accession"] for p in ds["pinned"]] == ["GSM9301", "GSM9304"]
    msg = ("Zzz: its knockout and control samples are not all read the same way (unfiltered matrices through the fixed "
           "cell call vs filtered matrices): refused")
    with pytest.raises(cm.Refused, match=re.escape(msg)):
        fe.verify_pins(ds, before)
    pl = _plan(y, dest)
    assert pl["status"] == "refused" and msg in pl["reason"]
    with pytest.raises(cm.Refused, match="did not pass"):
        fe.fetch(ACC, dest, y)
    assert calls == []
    # the second guard: with the first switched off in both module copies, the plan passes and fetch fetches, and
    # process refuses before any group call, writing nothing
    monkeypatch.setattr(fe, "check_unit_forms", lambda pins, per: None)
    monkeypatch.setattr(pr.fe, "check_unit_forms", lambda pins, per: None)
    assert _plan(y, dest)["status"] == "ok"
    fe.fetch(ACC, dest, y)
    assert sorted(Path(u).name for u in calls) == sorted([raw_name("GSM9301"), "GSM9304_matrix.mtx.gz",
                                                         "GSM9304_barcodes.tsv.gz", "GSM9304_features.tsv.gz"])
    with pytest.raises(cm.Refused, match="Zzz: an unfiltered and a filtered library in one unit: refused"):
        pr.run(ACC, dest, wt, refp, tmp_path / "processed", y,
               classifier=lambda *a: pytest.fail("a group call was made"))
    assert not (tmp_path / "processed" / ACC).exists()


def test_stage8_a_series_file_sample_counts_as_filtered_in_a_unit():
    """Stricter than the reviewer's line: a sample pinned with files [] (its cells in a pinned series file, always a
    filtered route) set against a raw-called sample of the same unit refuses too."""
    pins = [{"accession": "GSM1", "arm": "ko", "gene": "G", "stage": "E8.5"},
            {"accession": "GSM2", "arm": "control", "gene": "G", "stage": "E8.5"}]
    per = {"GSM1": {"files": [{"name": "a", "raw": True}]}, "GSM2": {"files": []}}
    with pytest.raises(cm.Refused, match="not all read the same way"):
        fe.check_unit_forms(pins, per)
    fe.check_unit_forms(pins, {"GSM1": {"files": [{"name": "a"}]}, "GSM2": {"files": []}})
    fe.check_unit_forms(pins, {"GSM1": {"files": [{"name": "a", "raw": True}]},
                               "GSM2": {"files": [{"name": "b", "raw": True}]}})


def test_stage8_a_matrix_line_that_is_not_three_fields_is_refused():
    """The first reviewer's lesser point 1: '1 2 3 9' was read with its first field as the index (r = [1, 2], c = [2,
    3], v = [9, 9]). Now a line of four fields, two fields or a fractional index is refused."""
    for block in (b"1 2 3 9\n2 3 4 9\n", b"1 2\n2 3\n", b"1.5 2 3\n", b"1 2 3\n2 3\n", b"1 2 x\n"):
        with pytest.raises(cm.Refused, match="not 'feature barcode value'"):
            pr._parse_entries(block, "L", (10, 10))
    r, c, v = pr._parse_entries(b"1 2 3\n2 3 4.5\n", "L", (10, 10))
    assert r.tolist() == [0, 1] and c.tolist() == [1, 2] and v.tolist() == [3.0, 4.5]


def test_stage8_a_long_comment_block_before_the_size_line_is_refused():
    head = b"%%MatrixMarket matrix coordinate integer general\n"
    ok = io.BytesIO(head + b"%c\n" * pr.RAW_MAX_COMMENT_LINES + b"3 4 0\n")
    assert pr._mtx_header(ok, "L") == (3, 4, 0)
    bad = io.BytesIO(head + b"%c\n" * (pr.RAW_MAX_COMMENT_LINES + 1) + b"3 4 0\n")
    with pytest.raises(cm.Refused, match="more than 1000 comment or empty lines"):
        pr._mtx_header(bad, "L")


def test_stage8_raw_call_arms_records_medians_and_ratios_only():
    rec = {"K1": {"cutoff": 1000.0, "cells_called": 3000, "umi_share_called": 0.8},
           "K2": {"cutoff": 1400.0, "cells_called": 5000, "umi_share_called": 0.6},
           "K3": {"cutoff": 1200.0, "cells_called": 4000, "umi_share_called": 0.7},
           "C1": {"cutoff": 600.0, "cells_called": 2000, "umi_share_called": 0.0}}
    got = pr.raw_call_arms(rec, ["K1", "K2", "K3"], ["C1"])
    assert got["median_ko"] == {"cutoff": 1200.0, "cells_called": 4000.0, "umi_share_called": pytest.approx(0.7)}
    assert got["ko_to_control"] == {"cutoff": 2.0, "cells_called": 2.0, "umi_share_called": None}
    # no guard reads it: the source holds no threshold on these keys
    src = (EXT / "process.py").read_text()
    assert "refused_call_mismatch" not in src


def test_stage8_the_readme_states_the_new_rules():
    text = (EXT / "README.md").read_text()
    for s in ("A raw file is taken only when it is the sample's only file of its", "never has its raw file fetched",
              "(the code refuses a unit whose arms are read differently)", "## Stage 8 (campaign/t3-ext-4)",
              "the barcode cap (60,000 per library) is the", 'plan it with `run_dgx.sh "$PWD" plan "GSE270969"`',
              "`raw_call_arms`", "HELD_EXTRA"):
        assert s in " ".join(text.split()) or s in text, s
    assert "a raw-called arm against a filtered arm is allowed by the code" not in text
    assert "must use the filtered one" not in text


# ---------------------------------------------------------------------------------------------------------------------
# macOS companion members (2026-11-20: a real archive of the pinned study held them)


def _tar_with_dirs(path: Path, members: dict[str, bytes | None]) -> Path:
    """A tar.gz whose members are written in order; a value of None makes a directory member."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as t:
        for name, data in members.items():
            ti = tarfile.TarInfo(name)
            if data is None:
                ti.type = tarfile.DIRTYPE
                t.addfile(ti)
            else:
                ti.size = len(data)
                t.addfile(ti, io.BytesIO(data))
    return path


def _mac_layout(X, names, folder="SIS1_raw_feature_bc_matrix", ad_size=176, ds_size=6148):
    """The member layout of an archive made on macOS: an AppleDouble companion beside the folder and beside every
    file, and a Finder .DS_Store with its own companion."""
    ad = b"\x00\x05\x16\x07" + b"\x00" * (ad_size - 4)
    out: dict[str, bytes | None] = {f"._{folder}": ad, folder: None, f"{folder}/._.DS_Store": ad,
                                    f"{folder}/.DS_Store": b"\x00" * ds_size}
    for k, v in _triple_bytes(X, names).items():
        out[f"{folder}/._{k}"] = ad
        out[f"{folder}/{k}"] = v
    return out


def test_macos_companion_members_are_skipped_unread_and_recorded(tmp_path):
    rng = np.random.default_rng(1)
    X, names, _ = _library(rng, _mu(rng))
    plain = pr.read_raw_sample([_raw_archive(tmp_path / "a" / raw_name("GSM1"), X, names)], "GSM1")
    mac = pr.read_raw_sample([_tar_with_dirs(tmp_path / "b" / raw_name("GSM1"), _mac_layout(X, names))], "GSM1")
    # the same cells, counts, symbols and call; '._matrix.mtx.gz' was never read as a matrix
    assert (mac[0] != plain[0]).nnz == 0 and mac[1:4] == plain[1:4]
    rec = dict(mac[4])
    assert rec.pop("mac_companions_skipped") == 6      # the folder's, .DS_Store and its own, and one per triple file
    assert rec == plain[4] and "mac_companions_skipped" not in plain[4]
    # the names reach a caller that asks for them, and only they are skipped
    with tarfile.open(tmp_path / "b" / raw_name("GSM1")) as t:
        skipped: list[str] = []
        items = pr.archive_members(t, "GSM1", raw_ok=True, skipped=skipped)
    assert sorted(i[1] for i in items) == ["barcodes", "features", "matrix"]
    assert all(s.rsplit("/", 1)[-1].startswith("._") or s.endswith("/.DS_Store") for s in skipped) and len(skipped) == 6


@pytest.mark.parametrize("change,match", [
    ("big_appledouble", "a '._' member larger than"),
    ("orphan_appledouble", "names no member beside it"),
    ("big_ds_store", "a '.DS_Store' member larger than"),
    ("other_dotfile", "not part of a 10x triple"),
    ("too_many", "hold more than"),
])
def test_a_macos_companion_outside_its_rule_is_refused(tmp_path, change, match):
    rng = np.random.default_rng(2)
    X, names, _ = _library(rng, _mu(rng))
    folder = "SIS1_raw_feature_bc_matrix"
    members = _mac_layout(X, names, folder)
    if change == "big_appledouble":
        members[f"{folder}/._matrix.mtx.gz"] = b"\x00" * (pr.MAC_APPLEDOUBLE_MAX + 1)
    elif change == "orphan_appledouble":
        members[f"{folder}/._notes.txt"] = b"\x00" * 176
    elif change == "big_ds_store":
        members[f"{folder}/.DS_Store"] = b"\x00" * (pr.MAC_DS_STORE_MAX + 1)
    elif change == "too_many":
        for i in range(20):       # 20 further folders, each with a 64 KiB .DS_Store: over 1 MiB skipped in total
            members[f"extra{i}"] = None
            members[f"extra{i}/.DS_Store"] = b"\x00" * pr.MAC_DS_STORE_MAX
    else:
        members[f"{folder}/.hidden"] = b"x"
    p = _tar_with_dirs(tmp_path / raw_name("GSM1"), members)
    with pytest.raises(pr.Refused, match=match):
        pr.read_raw_sample([p], "GSM1")
