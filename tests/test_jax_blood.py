"""Row T66's reducer (experiments/external/jax_blood.py) on synthetic nuclei only: no JAX file, no released file and no
score is read. The committed yamls (feng_relabel.yaml and jax_blood.yaml) are read for their sets."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "external"))
import feng_relabel as fr  # noqa: E402
import jax_blood as jb  # noqa: E402

CFG = jb.load_config()
FILLER = [f"Gene{i}" for i in range(420)]
UP, DOWN = "Gene1", "Gene2"  # planted late change in the erythroid nuclei: x4 and x0.25


def _symbols():
    sets = [g for s in CFG["sets"].values() for g in s["genes"]]
    return [*sets, "Lyve1", "Hba-a1", "Hbb-bs", *FILLER]


SYMBOLS = _symbols()
GENE_IDS = [f"ENSMUSG{i:05d}" for i in range(len(SYMBOLS))]
COL = {g: i for i, g in enumerate(SYMBOLS)}


def _profile(kind: str, late: bool, rng) -> np.ndarray:
    p = rng.gamma(2.0, 1.0, len(SYMBOLS))
    p[COL[UP]] = p[COL[DOWN]] = 30.0
    if kind == "ery":
        for g in CFG["sets"]["blood"]["genes"]:
            p[COL[g]] = 150.0
        if late:
            p[COL[UP]], p[COL[DOWN]] = 120.0, 7.5
    else:
        for g in CFG["sets"]["blood"]["genes"]:
            p[COL[g]] = 0.01  # no globin outside the planted group
    if kind == "cm":
        for g in CFG["sets"]["cm"]["genes"]:
            p[COL[g]] = 40.0
    return p / p.sum()


def _embryo(n: int, n_ery: int, n_cm: int, n_other: int, umis: int = 900, seed: int = 0):
    rng = np.random.default_rng([seed, n])
    late = n in jb.ROW_EMBRYOS["late"]
    base = rng.gamma(2.0, 1.0, len(SYMBOLS))  # one draw of filler levels per embryo
    rows, names, kinds = [], [], []
    for kind, k in (("ery", n_ery), ("cm", n_cm), ("other", n_other)):
        for i in range(k):
            p = _profile(kind, late, np.random.default_rng([seed, n, len(rows)]))
            p = 0.9 * p + 0.1 * base / base.sum()
            rows.append(rng.multinomial(umis, p / p.sum()))
            names.append(f"e{n}_{kind}_{i}")
            kinds.append(kind)
    rows.append(rng.multinomial(umis, _profile("other", late, rng)))  # a nucleus df_cell does not list
    names.append(f"e{n}_unlisted")
    kinds.append("unlisted")
    rows.append(np.zeros(len(SYMBOLS), dtype=np.int64))  # a nucleus with no counts
    names.append(f"e{n}_empty")
    kinds.append("empty")
    return sparse.csr_matrix(np.vstack(rows)), names, kinds


def _write_inputs(tmp_path, n_ery=40, poison="Blood", seed=0):
    src = tmp_path / "jax"
    (src / "meta").mkdir(parents=True)
    with open(src / "meta" / "df_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "gene_id", "gene_short_name", "chr"])
        for i, (g, s) in enumerate(zip(GENE_IDS, SYMBOLS, strict=True)):
            w.writerow([i, g, s, "1"])
    store = {}
    with open(src / "meta" / "df_cell.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "cell_id", "day", "embryo_id", "celltype_update", "major_trajectory", "experimental_batch"])
        k = 0
        for n in [*jb.ROW_EMBRYOS["early"], *jb.ROW_EMBRYOS["late"], 30]:
            C, names, kinds = _embryo(n, n_ery, 20, 20, seed=seed)
            store[n] = (C, list(GENE_IDS), names)
            day = "E8.5" if n < 13 else "E16.0"
            for name, kind in zip(names, kinds, strict=True):
                if kind == "unlisted":
                    continue
                w.writerow([k, name, day, f"embryo_{n}", poison if kind == "ery" else "x", "y",
                            "run_4" if n < 13 else "run_20"])
                k += 1
    return src, store


def _loader(store, seen):
    def load(path):
        n = int(re.search(r"embryo_(\d+)_gene_count", Path(path).name).group(1))
        seen.append(n)
        return store[n]
    return load


GATE = {"row": "T66", "clean": True}


def _reduce(tmp_path, **kw):
    src, store = _write_inputs(tmp_path, **kw)
    seen = []
    res = jb.reduce(CFG, src, GATE, loader=_loader(store, seen), verify=lambda p: {"sha256": "s"})
    return res, seen


# ---------------------------------------------------------------------------------------------------------------------
# The group and the direction


def test_a_planted_erythroid_group_is_recovered_and_its_change_is_the_direction(tmp_path):
    res, seen = _reduce(tmp_path)
    rec, arr = res["record"], res["arrays"]
    assert seen == sorted(n for e in jb.ENDS for n in jb.ROW_EMBRYOS[e])  # embryo 30 is never loaded
    for n, r in rec["per_embryo"].items():
        assert r["erythroid_assigned"] == 40 and r["erythroid_entering"] == 40, n
        assert r["nuclei_read"] == 80 and r["empty_dropped"] == 1 and r["nuclei_in_matrix"] == 82
    assert rec["embryos_used"] == {e: list(jb.ROW_EMBRYOS[e]) for e in jb.ENDS} and rec["has_direction"]
    d = dict(zip(arr["symbol"].tolist(), arr["direction"].tolist(), strict=True))
    assert d[UP] == pytest.approx(np.log(4.0), abs=0.35) and d[DOWN] == pytest.approx(np.log(0.25), abs=0.35)
    assert abs(d["Gene300"]) < 0.35
    assert rec["authors_labels_read"] is False and rec["columns_read"] == list(jb.CELL_WHITELIST)
    assert rec["sets"]["blood"]["found"] == ["Hba-x", "Hbb-y", "Gata1"]
    named = fr.named_genes(jb.scoring_config(CFG)) | fr.panel_genes()
    assert not named & {g for s in rec["sets"].values() for g in s["controls"]}  # N54's exclusions hold


def test_the_direction_is_the_equal_embryo_mean_of_downsampled_nuclei(tmp_path):
    res, _ = _reduce(tmp_path)
    arr = res["arrays"]
    early = np.isin(arr["embryo"], jb.ROW_EMBRYOS["early"])
    m_e, m_l = arr["embryo_mean_cp10k"][early].mean(0), arr["embryo_mean_cp10k"][~early].mean(0)
    assert np.allclose(arr["m_early"], m_e) and np.allclose(arr["m_late"], m_l)
    assert np.allclose(arr["direction"], np.log((m_l + 0.5) / (m_e + 0.5)))
    assert np.allclose(arr["embryo_mean_cp10k"].sum(1), 1e4)


def test_author_columns_never_change_the_result(tmp_path):
    a, _ = _reduce(tmp_path / "a", poison="Blood")
    b, _ = _reduce(tmp_path / "b", poison="Neuron")
    assert all(np.array_equal(a["arrays"][k], b["arrays"][k]) for k in a["arrays"])


def test_an_embryo_under_the_floor_is_left_out_and_an_end_without_embryos_gives_no_direction():
    v = np.full(4, 2500.0)
    per = {n: (v, 40) for n in jb.ROW_EMBRYOS["early"]} | {n: (v * 2, 40) for n in jb.ROW_EMBRYOS["late"]}
    per[59] = (np.full(4, 1e6), 29)  # under min_embryo_nuclei: left out
    r = jb.end_direction(per, 30, 0.5)
    assert 59 not in r["embryos_used"]["late"] and np.allclose(r["m_late"], 5000.0)
    per = {n: (v, 40) for n in jb.ROW_EMBRYOS["early"]} | {n: (None, 0) for n in jb.ROW_EMBRYOS["late"]}
    assert jb.end_direction(per, 30, 0.5)["direction"] is None


def test_downsampling_gives_exactly_the_target_without_replacement():
    rng = np.random.default_rng(3)
    D = rng.poisson(3.0, (30, 200))
    D[0] = 0
    D[1, :] = 0
    D[1, 5] = 399  # one nucleus just under the target
    C = sparse.csr_matrix(D)
    ds, keep = jb.downsample(C, 400, np.random.default_rng(1))
    lib = np.asarray(C.sum(1)).ravel()
    assert np.array_equal(keep, lib >= 400) and ds.shape[0] == int(keep.sum())
    assert np.all(np.asarray(ds.sum(1)).ravel() == 400)
    assert np.all(ds.toarray() <= C.toarray()[keep])
    ds2, _ = jb.downsample(C, 400, np.random.default_rng(1))
    assert (ds != ds2).nnz == 0  # one seed, one result


def test_end_means_is_gene_means_of_feng_relabel():
    rng = np.random.default_rng(5)
    C = sparse.csr_matrix(rng.poisson(2.0, (50, 30)) + 1)
    stage = np.array(["early"] * 20 + ["late"] * 30)
    L = fr.log1p_cp10k(C)
    want = fr.gene_means(L, stage, ["early", "late"])
    sums = {e: np.asarray(L[np.flatnonzero(stage == e)].sum(0)).ravel() for e in jb.ENDS}
    assert np.allclose(jb.end_means(sums, {"early": 20, "late": 30}), want)
    with pytest.raises(jb.Refused):
        jb.end_means(sums, {"early": 20, "late": 0})


# ---------------------------------------------------------------------------------------------------------------------
# The whitelist, the embryo lists and the sets


@pytest.mark.parametrize("cols", [["celltype_update"], ["cell_id", "major_trajectory"], ["cluster"], ["umap_1"],
                                  ["cell_id", "somite_count"]])
def test_the_label_whitelist_refuses_before_any_file_is_opened(tmp_path, cols):
    with pytest.raises(jb.Refused):
        jb.check_columns(cols)
    with pytest.raises(jb.Refused):
        jb.read_cells(tmp_path / "does_not_exist.csv", [1], cols)  # refused, not FileNotFoundError
    assert {"celltype_update", "major_trajectory"} <= set(jb.REFUSED_COLUMNS)
    assert jb.check_columns(jb.CELL_WHITELIST) == jb.CELL_WHITELIST


def test_the_cell_reader_returns_only_the_whitelist_and_the_listed_embryos(tmp_path):
    src, _ = _write_inputs(tmp_path)
    got = jb.read_cells(src / "meta" / "df_cell.csv", [1, 59])
    assert set(got["columns"]) == set(jb.CELL_WHITELIST)
    assert set(got["columns"]["embryo_id"]) == {"embryo_1", "embryo_59"} and got["rows_other_embryos"] > 0


@pytest.mark.parametrize("lists", [{"early": [*range(1, 13), 13], "late": [59, 60, 61, 62, 63]},
                                   {"early": list(range(1, 13)), "late": [54, 59, 60, 61, 62, 63]},
                                   {"early": list(range(1, 12)), "late": [59, 60, 61, 62, 63]},
                                   {"early": list(range(1, 13)), "late": [59, 60, 61, 62, 63], "mid": [22]}])
def test_an_embryo_outside_the_lists_refuses(lists):
    with pytest.raises(jb.Refused):
        jb.check_embryos(lists)
    with pytest.raises(jb.Refused):
        jb.end_of(22)
    bad = {**CFG, "embryos": lists}
    with pytest.raises(jb.Refused):
        jb.reduce(bad, Path("/nonexistent"), GATE, loader=lambda p: pytest.fail("no file is opened"),
                  verify=lambda p: pytest.fail("no file is opened"))


def test_an_embryo_id_outside_the_lists_in_the_count_store_is_never_read(tmp_path):
    _, seen = _reduce(tmp_path)
    assert 30 not in seen and set(seen) <= set(jb.ROW_EMBRYOS["early"]) | set(jb.ROW_EMBRYOS["late"])


def test_the_sets_are_n54s_released_data_sets_and_the_exclusions_n54s():
    n54 = yaml.safe_load((ROOT / "experiments" / "external" / "feng_relabel.yaml").read_text())
    want = {k: {g: v for g, v in s["genes"].items() if v == "released_e95"} for k, s in n54["sets"].items()}
    assert {k: v for k, v in want.items() if v} == {k: s["genes"] for k, s in CFG["sets"].items()}
    mine = {g for s in CFG["sets"].values() for g in [*s["genes"], *(s.get("excluded") or {})]}
    assert set(CFG["dropped"]) == fr.named_genes(n54) - mine
    assert CFG["assign"] == n54["assign"] and CFG["score"] == n54["score"]
    assert CFG["sets"]["blood"]["genes"] == {"Hba-x": "released_e95", "Hbb-y": "released_e95", "Gata1": "released_e95"}
    assert CFG["erythroid_label"] == "r_blood" and CFG["umi_target"] == 400 and CFG["eps"] == 0.5


def test_the_yaml_hash_is_pinned():
    assert fr.file_sha256(jb.CONFIG) == jb.YAML_SHA256
    assert CFG["embryos"] == {e: list(jb.ROW_EMBRYOS[e]) for e in jb.ENDS} and CFG["status"] == "frozen"


# ---------------------------------------------------------------------------------------------------------------------
# The gate


@pytest.fixture
def gate_env(tmp_path, monkeypatch):
    tries, reg = tmp_path / "TRIES.md", tmp_path / "EXTERNAL_DATA.md"
    tries.write_text("| T66 (one row) | T1 | x |\n")
    reg.write_text(f"| X003 | Qiu | ... {AMENDED} |\n")
    monkeypatch.setattr(jb, "TRIES", tries)
    monkeypatch.setattr(jb, "REGISTER", reg)
    monkeypatch.setattr(jb.fr, "committed_at_head", lambda p: Path(p).read_bytes())
    monkeypatch.setattr(jb.provenance, "git_state", lambda repo: ("c" * 40, []))
    monkeypatch.setattr(jb, "check_inputs", lambda *a, **k: {"inputs_checked": []})
    return tries, reg


AMENDED = ("amended for row T66: model input for the held Blood direction; reads the raw counts and, from df_cell, "
           "`cell_id`, `embryo_id`, `day` and `experimental_batch`, and from df_gene `gene_id` and `gene_short_name`")


def test_the_gate_passes_and_each_refusal_fires(tmp_path, gate_env, monkeypatch):
    tries, reg = gate_env
    dst = tmp_path / "out"
    assert jb.preflight(jb.CONFIG, "c" * 40, dst)["row"] == "T66"
    with pytest.raises(jb.Refused, match="HEAD"):
        jb.preflight(jb.CONFIG, "d" * 40, dst)
    dst.mkdir()
    with pytest.raises(jb.Refused, match="exists"):
        jb.preflight(jb.CONFIG, "c" * 40, dst)
    dst = tmp_path / "out2"
    monkeypatch.setattr(jb.provenance, "git_state", lambda repo: ("c" * 40, ["x.py"]))
    with pytest.raises(jb.Refused, match="not clean"):
        jb.preflight(jb.CONFIG, "c" * 40, dst)
    monkeypatch.setattr(jb.fr, "committed_at_head", lambda p: b"other")
    with pytest.raises(jb.Refused, match="committed version"):
        jb.preflight(jb.CONFIG, "c" * 40, dst)
    monkeypatch.setattr(jb, "YAML_SHA256", "0" * 64)
    with pytest.raises(jb.Refused, match="sha256"):
        jb.preflight(jb.CONFIG, "c" * 40, dst)
    with pytest.raises(jb.Refused, match="only the committed"):
        jb.preflight(tmp_path / "other.yaml", "c" * 40, dst)
    reg.write_text("| X003 | Qiu | Diagnostics only |\n")
    with pytest.raises(jb.Refused, match="X003"):
        jb.preflight(jb.CONFIG, "c" * 40, dst)
    tries.write_text("| T64 (one row) | T1 | x |\n")
    with pytest.raises(jb.Refused, match="no row"):
        jb.preflight(jb.CONFIG, "c" * 40, dst)


def test_the_run_needs_the_gate():
    with pytest.raises(jb.Refused, match="gate"):
        jb.reduce(CFG, Path("/nonexistent"), {}, loader=lambda p: pytest.fail("opened"))
    with pytest.raises(jb.Refused, match="gate"):
        jb.reduce(CFG, Path("/nonexistent"), {"row": "T66", "clean": False}, loader=lambda p: pytest.fail("opened"))


def test_main_refuses_without_the_row_before_any_jax_file(monkeypatch, tmp_path):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T64 (one row) | T1 | x |\n")
    monkeypatch.setattr(jb, "TRIES", tries)
    monkeypatch.setattr(jb, "reduce", lambda *a, **k: pytest.fail("nothing is read before the gate"))
    with pytest.raises(SystemExit, match="refused"):
        jb.main(["--commit", "c" * 40])


def test_the_register_check_reads_x003s_row_only(tmp_path):
    reg = tmp_path / "r.md"
    reg.write_text("| X002 | model input for the held Blood direction |\n| X003 | diagnostics only |\n")
    assert jb.register_amended(reg) is False
    reg.write_text("| X003 | model input for the held\n Blood direction |\n")
    assert jb.register_amended(reg) is False  # one table row is one line
    reg.write_text("| X003 | now a model input for the held Blood direction of row T66 |\n")
    assert jb.register_amended(reg) is False  # the phrase alone does not name the columns read
    reg.write_text(f"| X003 | {AMENDED} |\n")
    assert jb.register_amended(reg) is True
    for col in jb.AMENDMENT_COLUMNS:  # each column read must be named
        reg.write_text(f"| X003 | {AMENDED.replace(f'`{col}`', col)} |\n")
        assert jb.register_amended(reg) is False, col
    assert set(jb.AMENDMENT_COLUMNS) == {"cell_id", "embryo_id", "day", "experimental_batch", "gene_id",
                                         "gene_short_name"}


# ---------------------------------------------------------------------------------------------------------------------
# The inputs' format, checked before the output folder is made


def _inputs_on_disk(tmp_path, cell_header=None, gene_header=None, skip=None):
    src = tmp_path / "jax"
    (src / "meta").mkdir(parents=True)
    (src / "embryo").mkdir()
    (src / "meta" / "df_cell.csv").write_text(",".join(cell_header or ["", "cell_id", "day", "embryo_id",
                                                                       "celltype_update", "experimental_batch"])
                                              + "\n0,a,E8.5,embryo_1,Blood,run_4\n")
    (src / "meta" / "df_gene.csv").write_text(",".join(gene_header or ["", "gene_id", "gene_short_name"]) + "\n")
    for e in jb.ENDS:
        for n in jb.ROW_EMBRYOS[e]:
            if n != skip:
                (src / "embryo" / f"embryo_{n}_gene_count.rds").write_bytes(b"x")
    return src


def _ok(name):
    return None


def test_the_inputs_check_passes_and_reads_header_lines_only(tmp_path, monkeypatch):
    src = _inputs_on_disk(tmp_path)
    seen = []
    got = jb.check_inputs(CFG, src, importer=_ok, verify=lambda p: seen.append(Path(p).name) or {"sha256": "s"})
    assert len(got["inputs_checked"]) == 2 + 17 and len(seen) == 19 and "df_cell.csv" in seen
    def header_only(f):
        yield next(iter(f)).rstrip("\n").split(",")
        pytest.fail("a row past the header was read")

    monkeypatch.setattr(jb.csv, "reader", header_only)
    jb.check_inputs(CFG, src, importer=_ok, verify=lambda p: {})  # only the header lines are taken


@pytest.mark.parametrize("case,part", [
    ("rdata", "cannot be imported"),
    ("missing_rds", "embryo_61"),
    ("manifest", "sha256"),
    ("batch", "experimental_batch"),
    ("cell_id", "cell_id"),
    ("gene", "gene_short_name"),
])
def test_the_inputs_check_refuses(tmp_path, case, part):
    importer, verify = _ok, (lambda p: {})
    kw = {}
    if case == "missing_rds":
        kw["skip"] = 61
    elif case == "batch":
        kw["cell_header"] = ["", "cell_id", "day", "embryo_id"]
    elif case == "cell_id":
        kw["cell_header"] = ["", "day", "embryo_id", "experimental_batch"]
    elif case == "gene":
        kw["gene_header"] = ["", "gene_id", "symbol"]
    src = _inputs_on_disk(tmp_path, **kw)
    if case == "rdata":
        def importer(name):
            raise ModuleNotFoundError(f"No module named {name!r}")
    elif case == "manifest":
        def verify(p):
            raise ValueError(f"{p} does not match the sha256 recorded")
    with pytest.raises(jb.Refused, match=part):
        jb.check_inputs(CFG, src, importer=importer, verify=verify)


def test_a_format_failure_leaves_no_output_folder(tmp_path, monkeypatch):
    """The checks run inside preflight, so main refuses before it makes the folder or writes STARTED.json."""
    src = _inputs_on_disk(tmp_path, cell_header=["", "cell_id", "day", "embryo_id"])
    tries, reg = tmp_path / "TRIES.md", tmp_path / "EXTERNAL_DATA.md"
    tries.write_text("| T66 (one row) | T1 | x |\n")
    reg.write_text(f"| X003 | Qiu | ... {AMENDED} |\n")
    monkeypatch.setattr(jb, "TRIES", tries)
    monkeypatch.setattr(jb, "REGISTER", reg)
    monkeypatch.setattr(jb.fr, "committed_at_head", lambda p: Path(p).read_bytes())
    monkeypatch.setattr(jb.provenance, "git_state", lambda repo: ("c" * 40, []))
    dst = tmp_path / "data" / "derived" / "jax_blood"
    monkeypatch.setattr(jb.fr, "output_dir_ok", lambda p: dst)
    real_check = jb.check_inputs
    monkeypatch.setattr(jb, "check_inputs", lambda cfg, s, **k: real_check(cfg, src, importer=_ok,
                                                                            verify=lambda p: {}))
    monkeypatch.setattr(jb, "reduce", lambda *a, **k: pytest.fail("nothing is reduced after a refusal"))
    with pytest.raises(SystemExit, match="experimental_batch"):
        jb.main(["--commit", "c" * 40])
    assert not dst.exists() and not (dst / jb.STARTED).exists()


def test_preflight_runs_the_inputs_check_on_the_yamls_input_folder(tmp_path, gate_env, monkeypatch):
    called = []
    monkeypatch.setattr(jb, "check_inputs", lambda cfg, src, **k: called.append(Path(src)) or {})
    jb.preflight(jb.CONFIG, "c" * 40, tmp_path / "out")
    assert called == [jb.fio.REPO / CFG["inputs"]["dir"]]
    monkeypatch.setattr(jb, "check_inputs", lambda *a, **k: (_ for _ in ()).throw(jb.Refused("format")))
    with pytest.raises(jb.Refused, match="format"):
        jb.preflight(jb.CONFIG, "c" * 40, tmp_path / "out")


def test_the_source_reads_no_label_keyed_file():
    text = Path(jb.__file__).read_text()
    assert "stats/" not in text.split('"""', 2)[2] and "load_jax" not in text and "jax_sel" not in text
    assert "celltype_update" in text and "major_trajectory" in text  # only in the refusal list
    for name in ("celltype_update", "major_trajectory"):
        assert text.count(f'"{name}"') == 1
