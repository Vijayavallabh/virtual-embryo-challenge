"""Row T71's reducer (experiments/external/jax_cm.py) on synthetic nuclei only: no JAX file, no released file and no
score is read. The committed yamls (feng_relabel.yaml and jax_cm.yaml) are read for their sets; the ledger, register
and plan checks run on temporary files."""

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
import jax_cm as jc  # noqa: E402

CFG = jc.load_config()
FILLER = [f"Gene{i}" for i in range(420)]
UP, DOWN, SKM_UP = "Gene1", "Gene2", "Gene3"  # planted late change in CM nuclei (x4, x0.25) and in myonuclei (x8)
SHARED = ("Actc1", "Tnnt2", "Myl4")  # sarcomere genes of the cm set that fetal skeletal muscle also expresses


def _symbols():
    sets = [g for s in CFG["sets"].values() for g in [*s["genes"], *(s.get("excluded") or {})]]
    return [*dict.fromkeys([*sets, *CFG["dropped"], *FILLER])]


SYMBOLS = _symbols()
GENE_IDS = [f"ENSMUSG{i:05d}" for i in range(len(SYMBOLS))]
COL = {g: i for i, g in enumerate(SYMBOLS)}
EXTRA = [f"ENSMUSG9{k:05d}" for k in range(11)] + [f"{GENE_IDS[0]}.1"]
DAYS = {**{n: "E8.5" for n in range(1, 13)}, **{n: "E9.25" for n in range(13, 23)}, 59: "E15.5", 60: "E15.5",
        61: "E16.0", 62: "E16.5", 63: "E16.5"}


def _profile(kind: str, late: bool, rng) -> np.ndarray:
    p = rng.gamma(2.0, 1.0, len(SYMBOLS))
    p[COL[UP]] = p[COL[DOWN]] = p[COL[SKM_UP]] = 30.0
    for g in [*CFG["sets"]["cm"]["genes"], *jc.SKM_GENES]:
        p[COL[g]] = 0.01
    if kind == "cm":
        for g in CFG["sets"]["cm"]["genes"]:
            p[COL[g]] = 40.0
        if late:
            p[COL[UP]], p[COL[DOWN]] = 120.0, 7.5
    if kind == "skm":
        for g in jc.SKM_GENES:
            p[COL[g]] = 60.0
        for g in SHARED:
            p[COL[g]] = 40.0
        p[COL[SKM_UP]] = 240.0
    return p / p.sum()


def _embryo(n: int, n_cm: int, n_skm: int, n_other: int, umis: int = 900, seed: int = 0):
    rng = np.random.default_rng([seed, n])
    late = n >= 59
    base = rng.gamma(2.0, 1.0, len(SYMBOLS))
    rows, names = [], []
    for kind, k in (("cm", n_cm), ("skm", n_skm if late else 0), ("other", n_other)):
        for i in range(k):
            p = 0.9 * _profile(kind, late, np.random.default_rng([seed, n, len(rows)])) + 0.1 * base / base.sum()
            rows.append(rng.multinomial(umis, p / p.sum()))
            names.append(f"e{n}_{kind}_{i}")
    rows.append(np.zeros(len(SYMBOLS), dtype=np.int64))  # a nucleus with no counts
    names.append(f"e{n}_empty")
    return sparse.csr_matrix(np.vstack(rows)), names


def _write_inputs(tmp_path, mode="practice", n_cm=40, n_skm=60, seed=0, days=None, extra=True, cm_of=None):
    src = tmp_path / "jax"
    (src / "meta").mkdir(parents=True)
    with open(src / "meta" / "df_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "gene_id", "gene_type", "gene_short_name", "chr"])
        for i, (g, s) in enumerate(zip(GENE_IDS, SYMBOLS, strict=True)):
            w.writerow([i, g, "protein_coding", s, "1"])
    store = {}
    days = {**DAYS, **(days or {})}
    with open(src / "meta" / "df_cell.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "cell_id", "day", "embryo_id", "celltype_update", "experimental_batch"])
        k = 0
        lists = jc.ROW_EMBRYOS[mode]
        for n in [*lists["early"], *lists["late"], 30]:
            C, names = _embryo(n, (cm_of or {}).get(n, n_cm), n_skm, 20, seed=seed)
            ids = list(GENE_IDS)
            if extra:  # matrix rows the gene table does not list: appended, with counts
                X = np.random.default_rng([seed, n, 7]).poisson(3.0, (C.shape[0], len(EXTRA)))
                C = sparse.csr_matrix(np.hstack([C.toarray(), X]))
                ids = [*ids, *EXTRA]
            store[n] = (C, ids, names)
            for name in names:
                w.writerow([k, name, days.get(n, "E12.5"), f"embryo_{n}", "Cardiomyocyte", "run_x"])
                k += 1
    return src, store


def _loader(store, seen):
    def load(path):
        n = int(re.search(r"embryo_(\d+)_gene_count", Path(path).name).group(1))
        seen.append(n)
        return store[n]
    return load


def _gate(mode="practice"):
    return {"row": "T71", "clean": True, "mode": mode}


def _reduce(tmp_path, mode="practice", cfg=None, **kw):
    src, store = _write_inputs(tmp_path, mode=mode, **kw)
    seen = []
    res = jc.reduce(cfg or CFG, src, _gate(mode), loader=_loader(store, seen), verify=lambda p: {"sha256": "s"})
    return res, seen


# ---------------------------------------------------------------------------------------------------------------------
# The group, the guard and the direction


def test_a_planted_cm_group_is_recovered_and_myonuclei_go_to_the_skm_class(tmp_path):
    res, seen = _reduce(tmp_path)
    rec, arr = res["record"], res["arrays"]
    assert seen == sorted(n for e in jc.ENDS for n in jc.ROW_EMBRYOS["practice"][e])  # embryo 30 never loaded
    for n, r in rec["per_embryo"].items():
        # a planted CM nucleus may miss the margin against a two-gene class's noise; no other nucleus enters r_cm
        assert 35 <= r["cm_assigned"] <= 40 and r["cm_entering"] == r["cm_assigned"] and r["cm_skm_high"] <= 2, n
        late = n >= 59
        assert rec["labels_per_embryo"][n]["r_skm"] == (60 if late else 0)
    d = dict(zip(arr["symbol"].tolist(), arr["log_ratio"].tolist(), strict=True))
    assert d[UP] == pytest.approx(np.log(4.0), abs=0.35) and d[DOWN] == pytest.approx(np.log(0.25), abs=0.35)
    assert abs(d[SKM_UP]) < 0.35 and abs(d["Gene300"]) < 0.35  # the myonuclei's own change is not in the direction
    assert rec["authors_labels_read"] is False and rec["columns_read"] == ["cell_id", "embryo_id", "day"]
    assert rec["sets"]["skm"]["found"] == list(jc.SKM_GENES) and rec["mode"] == "practice"
    controls = {g for s in rec["sets"].values() for g in s["controls"]}
    assert not controls & set(jc.SKM_GENES)  # a set gene is never a control
    assert not controls & (fr.named_genes(jb.scoring_config(CFG)) | fr.panel_genes())


def test_without_the_skm_class_myonuclei_would_enter_the_cm_group(tmp_path):
    """Why the guard: the cm set's shared sarcomere genes carry myonuclei into r_cm when no skm class competes."""
    no_skm = {**CFG, "sets": {k: v for k, v in CFG["sets"].items() if k != "skm"}}
    with pytest.raises(jb.Refused, match="skm"):
        jc.check_config(no_skm)  # a yaml without the guard class is refused
    # score the late myonuclei directly: their cm score clears min_score, and only the skm class beats it
    C, _ = _embryo(59, 0, 30, 0)
    C = C[:-1]
    L = fr.log1p_cp10k(C)
    idx = fr.unique_index(SYMBOLS)
    cm = fr.set_score(L, [idx[g] for g in CFG["sets"]["cm"]["genes"]], [])
    skm = fr.set_score(L, [idx[g] for g in jc.SKM_GENES], [])
    assert np.all(cm > 0.3) and np.all(skm > cm)


def test_the_direction_is_the_log_ratio_over_five_days_of_delta(tmp_path):
    res, _ = _reduce(tmp_path)
    arr, rec = res["arrays"], res["record"]
    early = np.isin(arr["embryo"], jc.ROW_EMBRYOS["practice"]["early"])
    m_e, m_l = arr["embryo_mean_cp10k"][early].mean(0), arr["embryo_mean_cp10k"][~early].mean(0)
    assert np.allclose(arr["m_early"], m_e) and np.allclose(arr["m_late"], m_l)
    assert np.allclose(arr["log_ratio"], np.log((m_l + 0.5) / (m_e + 0.5)))
    delta = np.mean([15.5, 15.5, 16.0, 16.5, 16.5]) - 8.5
    assert rec["delta_days"] == pytest.approx(delta)
    assert np.allclose(arr["direction"], arr["log_ratio"] * 5.0 / delta)
    assert np.allclose(arr["embryo_mean_cp10k"].sum(1), 1e4)


def test_an_end_needs_three_embryos_at_the_floor(tmp_path):
    under = {59: 20, 60: 20, 61: 20}  # three of five late embryos under 30 CM nuclei: two left
    res, _ = _reduce(tmp_path, cm_of=under)
    rec = res["record"]
    assert rec["embryos_used"]["late"] == [62, 63] and rec["has_direction"] is False
    assert np.all(np.isnan(res["arrays"]["direction"]))
    days = {n: "E8.5" for n in range(1, 13)} | {59: "E15.5", 60: "E16.5", 61: "E16.0", 62: "E15.5", 63: "E16.5"}
    v = np.full(4, 2500.0)
    per = {n: (v, 40) for n in jc.ROW_EMBRYOS["practice"]["early"]} | {n: (v * 2, 40) for n in (59, 60, 61)}
    r = jc.end_direction(per | {62: (None, 0), 63: (v, 29)}, days, "practice", 30, 3, 0.5, 5.0)
    assert r["embryos_used"]["late"] == [59, 60, 61] and np.allclose(r["m_late"], 5000.0)
    assert r["delta_days"] == pytest.approx(16.0 - 8.5)


def test_matrix_genes_outside_the_table_change_nothing(tmp_path):
    a, _ = _reduce(tmp_path / "a", extra=True)
    b, _ = _reduce(tmp_path / "b", extra=False)
    for k in a["arrays"]:
        assert np.array_equal(a["arrays"][k], b["arrays"][k], equal_nan=a["arrays"][k].dtype.kind == "f"), k
    ra, rb = a["record"], b["record"]
    assert ra["genes_left_out"] == len(EXTRA) and rb["genes_left_out"] == 0 and ra["sets"] == rb["sets"]
    assert all(0 < s < 1 for s in ra["end_left_out_umi_share"].values())
    assert all(s == 0 for s in rb["end_left_out_umi_share"].values())


def test_author_columns_never_enter(tmp_path):
    src, _ = _write_inputs(tmp_path)
    got = jb.read_cells(src / "meta" / "df_cell.csv", [1], jc.CELL_WHITELIST)
    assert set(got["columns"]) == {"cell_id", "embryo_id", "day"}
    assert "experimental_batch" not in jc.CELL_WHITELIST and "celltype_update" not in jc.CELL_WHITELIST


# ---------------------------------------------------------------------------------------------------------------------
# The lists and the day labels


@pytest.mark.parametrize("mode,lists", [
    ("practice", {"early": [*range(1, 13), 13], "late": [59, 60, 61, 62, 63]}),
    ("practice", {"early": list(range(1, 13)), "late": [54, 59, 60, 61, 62, 63]}),
    ("practice", {"early": list(range(13, 23)), "late": [59, 60, 61, 62, 63]}),
    ("test", {"early": list(range(1, 13)), "late": [59, 60, 61, 62, 63]}),
    ("test", {"early": [*range(13, 23), 23], "late": [59, 60, 61, 62, 63]}),
    ("other", {"early": list(range(1, 13)), "late": [59, 60, 61, 62, 63]})])
def test_any_other_list_refuses(mode, lists):
    with pytest.raises(jb.Refused):
        jc.check_embryos(lists, mode)
    assert all(n not in jc.ROW_EMBRYOS[m][e] for m in jc.MODES for e in jc.ENDS for n in range(23, 59))


def test_the_yaml_holds_both_modes_lists_and_refuses_another():
    assert CFG["embryos"] == {m: {e: list(jc.ROW_EMBRYOS[m][e]) for e in jc.ENDS} for m in jc.MODES}
    bad = {**CFG, "embryos": {**CFG["embryos"], "test": {"early": list(range(13, 22)), "late": [59, 60]}}}
    with pytest.raises(jb.Refused):
        jc.check_config(bad)
    with pytest.raises(jb.Refused):
        jc.reduce(CFG, Path("/nonexistent"), _gate("bogus"), loader=lambda p: pytest.fail("opened"))


def test_test_mode_reads_13_to_22_and_59_to_63(tmp_path):
    res, seen = _reduce(tmp_path, mode="test")
    assert seen == [*range(13, 23), 59, 60, 61, 62, 63] and res["record"]["mode"] == "test"
    assert res["record"]["delta_days"] == pytest.approx(np.mean([15.5, 15.5, 16.0, 16.5, 16.5]) - 9.25)


@pytest.mark.parametrize("mode,days,match", [
    ("practice", {3: "E9.0"}, "practice tasks' window"),
    ("practice", {60: "E11.0"}, "protected window|practice tasks' window"),
    ("test", {14: "E9.75"}, "protected window"),
    ("test", {62: "E13.5"}, "protected window"),
    ("test", {15: "E8.0"}, None)])
def test_the_day_labels_are_checked_before_any_count(tmp_path, mode, days, match):
    src, store = _write_inputs(tmp_path, mode=mode, days=days)
    if match is None:
        jc.reduce(CFG, src, _gate(mode), loader=_loader(store, []), verify=lambda p: {"sha256": "s"})
        return
    with pytest.raises(jb.Refused, match=match):
        jc.reduce(CFG, src, _gate(mode), loader=lambda p: pytest.fail("a count was read"),
                  verify=lambda p: {"sha256": "s"})


def test_an_embryo_with_two_day_labels_refuses():
    col = {"embryo_id": np.array(["embryo_1", "embryo_1"] + [f"embryo_{n}" for n in [*range(2, 13), *range(59, 64)]]),
           "day": np.array(["E8.5", "E8.25"] + ["E8.5"] * 11 + ["E16.0"] * 5)}
    with pytest.raises(jb.Refused, match="2 day labels"):
        jc.embryo_days(col, "practice")
    with pytest.raises(jb.Refused, match="test's bracket"):
        jc.embryo_days({"embryo_id": np.array([f"embryo_{n}" for n in [*range(13, 23), *range(59, 64)]]),
                        "day": np.array(["E9.5"] * 10 + ["E14.0"] * 4 + ["E8.5"])}, "test")


# ---------------------------------------------------------------------------------------------------------------------
# The yaml


def test_the_sets_are_n54s_plus_skm_and_the_settings_n54s():
    n54 = yaml.safe_load((ROOT / "experiments" / "external" / "feng_relabel.yaml").read_text())
    mine = {k: v for k, v in CFG["sets"].items() if k != "skm"}
    assert list(mine) == list(n54["sets"])  # the same order: each set keeps its control seed index
    for k, s in n54["sets"].items():
        assert mine[k]["genes"] == s["genes"] and mine[k]["role"] == s["role"], k
        assert set(mine[k].get("excluded") or {}) == set(s.get("excluded") or {}), k
        assert mine[k].get("label") == s.get("label"), k
    assert set(CFG["dropped"]) == set(n54["dropped"])
    assert CFG["assign"] == n54["assign"] and CFG["score"] == n54["score"]
    assert CFG["sets"]["skm"] == {"role": "competing", "label": "r_skm",
                                  "genes": {g: "X010_T71" for g in jc.SKM_GENES}}
    assert list(CFG["sets"])[-1] == "skm"
    srcs = {src for s in CFG["sets"].values() for src in s["genes"].values()}
    assert srcs <= set(CFG["sources"]) and srcs - {"X010_T71"} <= set(n54["sources"])


def test_the_settings_are_the_specifications():
    assert (CFG["umi_target"], CFG["downsample_seed"], CFG["min_embryo_nuclei"], CFG["min_end_embryos"]) == (
        400, 20261123, 30, 3)
    assert (CFG["eps"], CFG["bracket_days"], CFG["purity_score"], CFG["cm_label"]) == (0.5, 5.0, 0.3, "r_cm")
    assert CFG["row"] == jc.ROW == "T71" and CFG["status"] == "frozen" and CFG["gene_universe"] == "gene_table"
    assert CFG["output"]["practice"] == "data/derived/jax_cm_t71_practice"
    assert CFG["output"]["test"] == "data/derived/jax_cm_t71_test"


def test_the_yaml_hash_is_pinned():
    """A committed text file: git gives the same bytes on both machines."""
    assert fr.file_sha256(jc.CONFIG) == jc.YAML_SHA256


# ---------------------------------------------------------------------------------------------------------------------
# The gate


AMENDED = ("Amendment for row T71: T1 CM direction input, row T71; embryos 1-12 and 59-63 for practice; reads "
           "`cell_id`, `embryo_id` and `day`")
SKM = "| X010 | eight papers | ... skeletal-muscle competing class of row T71: `Tnnt3`, `Mylpf` (a source) |"


@pytest.fixture
def gate_env(tmp_path, monkeypatch):
    tries, reg, plan = tmp_path / "TRIES.md", tmp_path / "EXTERNAL_DATA.md", tmp_path / "PLAN.md"
    tries.write_text("| T71 (one row) | T1 | x |\n")
    reg.write_text(f"| X003 | Qiu | ... {AMENDED} |\n{SKM}\n")
    plan.write_text("## log\n\n- 2026-11-25 (a) the captain rules on row T71: a test-mode reduction.\n")
    monkeypatch.setattr(jc, "TRIES", tries)
    monkeypatch.setattr(jc, "REGISTER", reg)
    monkeypatch.setattr(jc, "PLAN", plan)
    monkeypatch.setattr(jc.fr, "committed_at_head", lambda p: Path(p).read_bytes())
    monkeypatch.setattr(jc.provenance, "git_state", lambda repo: ("c" * 40, []))
    monkeypatch.setattr(jc, "check_inputs", lambda *a, **k: {"inputs_checked": []})
    return tries, reg, plan


def test_the_gate_passes_and_each_refusal_fires(tmp_path, gate_env, monkeypatch):
    tries, reg, plan = gate_env
    dst = tmp_path / "out"
    assert jc.preflight(jc.CONFIG, "c" * 40, dst)["mode"] == "practice"
    with pytest.raises(jb.Refused, match="X003"):
        jc.preflight(jc.CONFIG, "c" * 40, dst, "test", "2026-11-25 (a)")  # no test phrase yet
    reg.write_text(f"| X003 | Qiu | ... {AMENDED}; T1 CM test early end, embryos 13-22 |\n{SKM}\n")
    for ruling in (None, "2026-11-25 (b)", "2026-11-25"):
        with pytest.raises(jb.Refused, match="captain's ruling"):
            jc.preflight(jc.CONFIG, "c" * 40, dst, "test", ruling)
    assert jc.preflight(jc.CONFIG, "c" * 40, dst, "test", "2026-11-25 (a)")["ruling"] == "2026-11-25 (a)"
    with pytest.raises(jb.Refused, match="HEAD"):
        jc.preflight(jc.CONFIG, "d" * 40, dst)
    dst.mkdir()
    with pytest.raises(jb.Refused, match="exists"):
        jc.preflight(jc.CONFIG, "c" * 40, dst)
    dst = tmp_path / "out2"
    monkeypatch.setattr(jc.provenance, "git_state", lambda repo: ("c" * 40, ["x.py"]))
    with pytest.raises(jb.Refused, match="not clean"):
        jc.preflight(jc.CONFIG, "c" * 40, dst)
    monkeypatch.setattr(jc.fr, "committed_at_head", lambda p: b"other")
    with pytest.raises(jb.Refused, match="committed version"):
        jc.preflight(jc.CONFIG, "c" * 40, dst)
    monkeypatch.setattr(jc, "YAML_SHA256", "0" * 64)
    with pytest.raises(jb.Refused, match="sha256"):
        jc.preflight(jc.CONFIG, "c" * 40, dst)
    with pytest.raises(jb.Refused, match="only the committed"):
        jc.preflight(tmp_path / "other.yaml", "c" * 40, dst)
    reg.write_text(f"| X003 | Qiu | ... {AMENDED} |\n| X010 | eight papers | `Tnnt3` `Mylpf` |\n")
    with pytest.raises(jb.Refused, match="X010"):
        jc.preflight(jc.CONFIG, "c" * 40, dst)
    reg.write_text(f"| X003 | Qiu | model input for the held Blood direction `cell_id` `embryo_id` `day` |\n{SKM}\n")
    with pytest.raises(jb.Refused, match="X003"):
        jc.preflight(jc.CONFIG, "c" * 40, dst)
    tries.write_text("| T68 (one row) | T1 | x |\n")
    with pytest.raises(jb.Refused, match="no row"):
        jc.preflight(jc.CONFIG, "c" * 40, dst)


def test_the_register_checks_read_their_rows_only(tmp_path):
    reg = tmp_path / "r.md"
    reg.write_text(f"| X002 | {AMENDED} |\n| X003 | Qiu | nothing |\n{SKM}\n")
    assert not jc.register_amended("practice", reg)
    reg.write_text(f"| X003 | Qiu | {AMENDED} |\n| X010 | skeletal-muscle competing class of row T71 `Tnnt3` |\n")
    assert jc.register_amended("practice", reg) and not jc.skm_registered(reg)  # one gene named is not two
    reg.write_text(f"| X003 | Qiu | {AMENDED} |\n| X011 | {SKM} |\n")
    assert not jc.skm_registered(reg)
    reg.write_text(f"{SKM}\n")
    assert jc.skm_registered(reg)


def test_the_run_needs_the_gate():
    with pytest.raises(jb.Refused, match="gate"):
        jc.reduce(CFG, Path("/nonexistent"), {}, loader=lambda p: pytest.fail("opened"))
    with pytest.raises(jb.Refused, match="gate"):
        jc.reduce(CFG, Path("/nonexistent"), {"row": "T68", "clean": True, "mode": "practice"},
                  loader=lambda p: pytest.fail("opened"))


def test_main_refuses_without_the_row_before_any_jax_file(monkeypatch, tmp_path):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| T68 (one row) | T1 | x |\n")
    monkeypatch.setattr(jc, "TRIES", tries)
    monkeypatch.setattr(jc.fr, "output_dir_ok", lambda d: tmp_path / "never")
    monkeypatch.setattr(jc, "reduce", lambda *a, **k: pytest.fail("no read"))
    with pytest.raises(SystemExit, match="refused"):
        jc.main(["--commit", "c" * 40, "--mode", "practice"])
    assert not (tmp_path / "never").exists()


def test_the_inputs_check_reads_header_lines_and_the_modes_files(tmp_path):
    src, _ = _write_inputs(tmp_path)
    for n in [*range(1, 13), *range(59, 64)]:
        p = src / "embryo" / f"embryo_{n}_gene_count.rds"
        p.parent.mkdir(exist_ok=True)
        p.write_text("x")
    checked = []
    got = jc.check_inputs(CFG, src, "practice", importer=lambda m: None, verify=lambda p: checked.append(p.name))
    assert len(got["inputs_checked"]) == 19 and len(checked) == 19
    with pytest.raises(jb.Refused, match="not in place"):
        jc.check_inputs(CFG, src, "test", importer=lambda m: None, verify=lambda p: None)


def test_the_source_reads_no_label_keyed_file():
    text = (ROOT / "experiments" / "external" / "jax_cm.py").read_text()
    body = text.split('"""', 2)[2]
    assert "stats/" not in body and "celltype" not in body and "experimental_batch\"" not in body
    assert "lpm" not in body.lower()
