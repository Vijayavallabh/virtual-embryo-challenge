"""Stage 7, piece 2 of experiments/t3_ext (campaign/t3-ext-3): the intake for per-sample UNFILTERED matrices (a
dataset's raw_matrices declaration, the one fixed cell call of t3x_common, the sparse two-pass reader). Synthetic SOFT
text, synthetic counts and a stub download only: nothing here touches the network or real data."""

from __future__ import annotations

import gzip
import importlib.util
import inspect
import io
import json
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


def _soft() -> str:
    lines = [f"^SERIES = {ACC}", f"!Series_title = A synthetic series {ACC}", "!Series_summary = embryos",
             "!Series_overall_design = one embryo per library", "!Series_contact_name = A,,Person"]
    for gsm, rec in SAMPLES.items():
        lines.append(f"^SAMPLE = {gsm}")
        for k, vals in rec.items():
            for v in (vals if isinstance(vals, list) else [vals]):
                lines.append(f"!Sample_{k} = {v}")
    return "\n".join(lines) + "\n"


def _plan(y, dest):
    return fe.plan(ACC, dest, y, metas={ACC: fe.parse_soft(_soft())})


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
                 "features": len(PANEL) + EXTRA + 1, "entries": int((X > 0).sum())}
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
    assert other == [{"name": raw_name("GSM9304"), "reason": "an unfiltered matrix of a sample that lists a filtered "
                      "file of its own: the filtered one is used"}]
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
                            "umi_share_called", "features", "entries"}
    u = qc["units"]["Zzz"]
    assert u["status"] == "ok" and u["raw_called"] is True and u["raw_called_samples"] == ["GSM9301", "GSM9302"]
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
