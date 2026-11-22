"""T1 candidate rc18 (ledger rows N58 and T69): the outside heart atlas's direction over E9.5 to E18.

Held here, on synthetic arrays, hand-made files and throwaway repositories only (no Feng file, released file, board or
score is read; no launcher is run on real data):
  - row N58's labelling (experiments/external/feng_relabel_e18.py): the E18 cells scored against row N54's recorded
    control sets and labelled by the frozen yaml's rule; N54's E9.5 and E14.5 rows copied byte for byte; N54's file and
    record never written; its refusals in order;
  - feng_dirs' E9.5 -> E18 chord (`recomputed_e18_directions`, `e18_late_term`): its checks and refusals;
  - row T69's read (experiments/t1_m1/rc18_read.py): the arms, the size and permuted controls, the statement;
  - the read's launcher experiments/t1_m1/rc18.sh: its refusal paths and form.
The builder option (build_test.py --feng-bracket E18) and its test-build launcher are not on this branch: by the T1
order rule (decision log 2026-11-21 (y)) they are written, with their tests, only after an AHEAD reading of row T69 and
the captain's yes.
"""

from __future__ import annotations

import copy
import csv
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import feng2022_geo  # noqa: E402
import feng_dirs  # noqa: E402
import feng_relabel as fr  # noqa: E402
import feng_relabel_e18 as fre  # noqa: E402
import filterio as fio  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
from test_feng_relabel import GATE, GENES, SPEC, cells, counts, syn_cfg  # noqa: E402
from test_t1_scope import HELD  # noqa: E402
from test_t1_state_transport import _gmom_ctx, classified_pair  # noqa: E402

from vec import harness  # noqa: E402
from vec.core import counts_of  # noqa: E402
from vec.t1_scope import load_spec  # noqa: E402

READ = ROOT / "experiments" / "t1_m1" / "rc18_read.py"
LAUNCHER = ROOT / "experiments" / "t1_m1" / "rc18.sh"
GROUPS4 = ("CM", "Endocardium", "Endothelium", "Epicardium")
EPS = 0.5


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def rr18(monkeypatch, tmp_path):
    mod = _load("rc18_read_module", READ)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


# ---------------------------------------------------------------------------------------------------------------------
# Row N58: the E18 labelling


def _n54_like(cfg, types, seed=3):
    """fr.relabel on E9.5 and E14.5 cells: its labels and the record-style sets (genes found, controls)."""
    stage = np.array(["E9.5" if i % 2 == 0 else "E14.5" for i in range(len(types))])
    C = counts(types, seed=seed)
    return C, stage, fr.relabel(C, GENES, stage, cfg)


def test_the_fixed_controls_give_the_frozen_labellers_scores_and_labels():
    cfg = syn_cfg()
    types = cells(SPEC) + cells(SPEC)
    C, stage, res = _n54_like(cfg, types)
    got = fre.relabel_fixed(C, GENES, np.array(["E18"] * len(types)), cfg, res["sets"])
    assert got["labels"].tolist() == res["labels"].tolist()  # N54's controls: the same labels, cell for cell
    for name in cfg["sets"]:
        assert got["scores"][name].tobytes() == res["scores"][name].tobytes()
    assert set(got["unassigned_share"]) == {"E18"} and got["cap"]["passed"]
    assert got["sets"]["cm"]["controls"] == res["sets"]["cm"]["controls"]
    with pytest.raises(fre.Refused, match="outside"):
        fre.relabel_fixed(C, GENES, stage, cfg, res["sets"])  # only E18 cells are labelled here


def test_the_fixed_sets_refuse_anything_but_n54s_record():
    cfg = syn_cfg()
    _, _, res = _n54_like(cfg, cells(SPEC))
    sets = copy.deepcopy(res["sets"])
    sets["cm"]["found"] = sets["cm"]["found"][:-1]
    with pytest.raises(fre.Refused, match="differ from row N54"):
        fre.fixed_sets(GENES, cfg, sets)
    sets = copy.deepcopy(res["sets"])
    sets["epi"]["controls"] = [*sets["epi"]["controls"], "NotAGene"]
    with pytest.raises(fre.Refused, match="control genes are missing"):
        fre.fixed_sets(GENES, cfg, sets)
    sets = copy.deepcopy(res["sets"])
    sets["blood"]["controls"] = []
    with pytest.raises(fre.Refused, match="control genes are missing"):
        fre.fixed_sets(GENES, cfg, sets)
    sets = copy.deepcopy(res["sets"])
    del sets["vasc_split"]
    with pytest.raises(fre.Refused, match="vasc_split"):
        fre.fixed_sets(GENES, cfg, sets)
    set_idx, ctrl, shown = fre.fixed_sets(GENES, cfg, res["sets"])
    assert [GENES[i] for i in ctrl["cm"]] == res["sets"]["cm"]["controls"]


def _group_stats(stage_value, types, seed):
    C = counts(types, seed=seed)
    n = len(types)
    labels = np.array(["r_cm" if t == "cm" else "r_epicardial" if t == "epi" else "r_vascular_ec" for t in types])
    return feng2022_geo.group_stats(C, np.array([stage_value] * n), np.array(["LV"] * n), labels,
                                    np.array([f"{stage_value}_{i}" for i in range(n)]), GENES)


def _n54_stats():
    a = _group_stats("E9.5", cells(SPEC), 1)
    b = _group_stats("E14.5", cells(SPEC), 2)
    return {"keys": np.array([*a["keys"], *b["keys"]]), "genes": a["genes"],
            **{f: np.concatenate([a[f], b[f]]) for f in fre.ROW_FIELDS}}


def test_the_copied_rows_are_n54s_byte_for_byte_and_a_change_is_seen():
    n54 = _n54_stats()
    e18 = _group_stats("E18", cells(SPEC), 3)
    comb = fre.combine(n54, e18)
    assert [str(k) for k in comb["keys"][:len(n54["keys"])]] == [str(k) for k in n54["keys"]]
    ok = fre.copied_rows_identical(comb, n54)
    assert ok["identical"] and all(ok[k] for k in ("genes", "keys", "rows", "other_stages"))
    bad = {k: np.array(v, copy=True) for k, v in comb.items()}
    bad["S_cp"][1, 5] = np.nextafter(bad["S_cp"][1, 5], np.inf)  # one ulp on one copied row
    assert not fre.copied_rows_identical(bad, n54)["identical"]
    bad = {k: np.array(v, copy=True) for k, v in comb.items()}
    bad["keys"] = np.array([*bad["keys"][:-1], bad["keys"][0]])  # an E9.5 key twice
    assert not fre.copied_rows_identical(bad, n54)["identical"]
    bad = dict(comb, n=comb["n"].astype(np.float32))
    assert not fre.copied_rows_identical(bad, n54)["identical"]
    with pytest.raises(fre.Refused, match="one gene list"):
        fre.combine(n54, dict(e18, genes=np.array(list(reversed(e18["genes"])))))
    with pytest.raises(fre.Refused, match="other than E18"):
        fre.combine(n54, _group_stats("E16.5", cells(SPEC), 3))
    with pytest.raises(fre.Refused, match="outside"):
        fre.combine(e18, e18)


def _obs(path: Path, stages: list[str], author: list[str]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["", "orig.ident", "strain", "stage", "zone", "cell_type_spec", "seurat_clusters"])
        for i, s in enumerate(stages):
            w.writerow([f"c{i}", "s1", "CD1", s, ["LA", "RA", "LV", "RV"][i % 4], author[i], str(i % 7)])


def _src(root: Path, types: list[str], stages: list[str], author: list[str]) -> Path:
    src = root / "src"
    src.mkdir(parents=True, exist_ok=True)
    _obs(src / "cd1_legal_obs.csv", stages, author)
    sparse.save_npz(src / "cd1_legal_counts.npz", counts(types, seed=5))
    (src / "cd1_legal_genes.txt").write_text("\n".join(GENES) + "\n")
    fio.write_manifest(src, {"outputs": {f: {"sha256": fio.sha256_file(src / f)} for f in
                                         ("cd1_legal_obs.csv", "cd1_legal_counts.npz", "cd1_legal_genes.txt")}})
    return src


def _world(tmp_path, monkeypatch, author_mod=0):
    """A synthetic clean reduction holding E9.5, E14.5, E16.5 and E18 cells, row N54's run on it (fr.run, its record's
    bytes) and the module constants pointed at the synthetic N54 file and settings."""
    types = [t for t in cells(SPEC) for _ in range(3)]
    stages = [("E9.5", "E14.5", "E18", "E16.5")[i % 4] for i in range(len(types))]
    author = [f"x{(i * 7 + author_mod) % 5}" for i in range(len(types))]
    src = _src(tmp_path, types, stages, author)
    cfg = syn_cfg(status="frozen")
    cfg_path = tmp_path / "relabel.yaml"
    import yaml

    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    n54_dir = tmp_path / "n54"
    fr.run(cfg, cfg_path, src, n54_dir, GATE)
    n54_path = n54_dir / cfg["output"]["stats"]
    monkeypatch.setattr(fre, "YAML_SHA256", fio.sha256_file(cfg_path))
    monkeypatch.setattr(fre, "N54_STATS_SHA256", fio.sha256_file(n54_path))
    raw = (n54_dir / "MANIFEST.json").read_bytes()
    gate = {"row": "N70", "commit": "c" * 40, "clean": True, "config_sha256": "x", "n54_record_sha256": "y"}
    return cfg, src, n54_path, raw, gate, types, stages


def test_the_labelling_copies_n54s_rows_labels_e18_alone_and_never_writes_n54s_files(tmp_path, monkeypatch):
    cfg, src, n54_path, raw, gate, types, stages = _world(tmp_path, monkeypatch)
    before = n54_path.read_bytes()
    dst = tmp_path / "e18"
    man = fre.run(cfg, src, dst, gate, raw, n54_path)
    assert n54_path.read_bytes() == before  # N54's file untouched
    assert man["source"] == "feng_relabel_e18" and man["stages"] == ["E18"]
    assert man["stages_copied"] == ["E9.5", "E14.5"] and man["copied_rows"]["identical"] is True
    assert man["authors_labels_read"] is False and man["usable"] is True and man["gate"] == gate
    assert set(man["unassigned_share"]) == {"E18"}
    st = fre.load_raw(dst / fre.STATS)
    n54 = fre.load_raw(n54_path)
    assert fre.copied_rows_identical(st, n54)["identical"]
    assert {str(k).split("|")[0] for k in st["keys"]} == {"E9.5", "E14.5", "E18"}  # no E16.5 cell taken
    e18_rows = [i for i, k in enumerate(st["keys"]) if str(k).startswith("E18|")]
    assert len(e18_rows) == man["outputs"][fre.STATS]["e18_rows"] > 0
    # the E18 labels are those N54's controls give (the frozen labeller's own scores on the same cells)
    rec = json.loads(raw)
    keep = np.flatnonzero(np.array(stages) == "E18")
    C = counts(types, seed=5).tocsr()[keep]
    want = fre.relabel_fixed(C, GENES, np.array(["E18"] * len(keep)), cfg, rec["sets"])["labels"]
    table = list(csv.reader((dst / fre.CELLS).read_text().splitlines()))
    assert [r[3] for r in table[1:]] == want.tolist()
    assert fio.verify_manifest_output(dst / fre.STATS)["usable"] is True
    # the outputs do not depend on the authors' column
    other = tmp_path / "other"
    cfg2, src2, n54_2, raw2, gate2, *_ = _world(other, monkeypatch, author_mod=3)
    fre.run(cfg2, src2, other / "e18", gate2, raw2, n54_2)
    assert (other / "e18" / fre.STATS).read_bytes() == (dst / fre.STATS).read_bytes()
    assert (other / "e18" / fre.CELLS).read_bytes() == (dst / fre.CELLS).read_bytes()


def test_the_labelling_refuses_other_inputs_another_n54_file_or_no_gate(tmp_path, monkeypatch):
    cfg, src, n54_path, raw, gate, *_ = _world(tmp_path, monkeypatch)
    with pytest.raises(fre.Refused, match="no clean gate"):
        fre.run(cfg, src, tmp_path / "a", {"row": "T69", "clean": True}, raw, n54_path)
    monkeypatch.setattr(fre, "N54_STATS_SHA256", "0" * 64)
    with pytest.raises(fre.Refused, match="statistics sha256"):
        fre.run(cfg, src, tmp_path / "b", gate, raw, n54_path)
    monkeypatch.setattr(fre, "N54_STATS_SHA256", fio.sha256_file(n54_path))
    rec = json.loads(raw)
    rec["inputs"]["counts"]["sha256"] = "0" * 64
    with pytest.raises(fre.Refused, match="not the ones row N54 read"):
        fre.run(cfg, src, tmp_path / "c", gate, json.dumps(rec).encode(), n54_path)
    for change, msg in (({"route": "none"}, "route recomputed"), ({"usable": False}, "not usable"),
                        ({"stages": ["E9.5"]}, "stages"), ({"sets": {}}, "no control set")):
        with pytest.raises(fre.Refused, match=msg):
            fre.n54_record(json.dumps({**json.loads(raw), **change}).encode(), cfg)
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists() and not (tmp_path / "c").exists()


def test_the_labellers_paths_are_its_own():
    assert fre.RECORD != fr.RECORD and fre.RECORD.parent.name == "feng_relabel_e18"
    assert fre.E18_DIR == "data/derived/feng_relabel_e18" and fre.E18_DIR != fr.load_config()["output"]["dir"]
    assert fre.YAML_SHA256 == fio.sha256_file(fr.CONFIG)  # the frozen settings, unchanged
    code = (ROOT / "experiments" / "external" / "feng_relabel_e18.py").read_text().split('"""', 2)[2]
    assert "cell_type_spec" not in code and "feng_relabel_c" not in code
    assert "fr.relabel(" not in code and "control_sets(" not in code  # no control gene is drawn here
    assert "write_record(" in code and "fr.write_record" not in code and "fr.RECORD.write" not in code


def test_the_x002_check_and_the_row_line(tmp_path):
    reg = tmp_path / "EXTERNAL_DATA.md"
    reg.write_text("| X001 | feng_relabel_e18 |\n| X002 | Feng 2022 ... feng_relabel_late |\n")
    assert not fre.x002_names_this_use(reg)
    reg.write_text("| X002 | Feng 2022 ...; amended: data/derived/feng_relabel_e18 (E18) |\n")
    assert fre.x002_names_this_use(reg)
    tries = tmp_path / "TRIES.md"
    tries.write_text("| N70 (target-free) | T1 | experiments/external/feng_relabel_e18.py |\n| N7 | x |\n")
    assert "feng_relabel_e18.py" in fre.row_line("N70", tries) and fre.row_line("N7", tries) == "| N7 | x |"


def test_the_preflight_refuses_in_order(tmp_path, monkeypatch):
    dst = tmp_path / "data" / "derived" / "feng_relabel_e18"
    with pytest.raises(fre.Refused, match="not an N row"):
        fre.preflight("T69", dst)
    real = fre.YAML_SHA256
    monkeypatch.setattr(fre, "YAML_SHA256", "0" * 64)
    with pytest.raises(fre.Refused, match="frozen"):
        fre.preflight("N70", dst)
    monkeypatch.setattr(fre, "YAML_SHA256", real)
    stats = fr.load_config()["output"]["stats"]
    good = {"config": {"sha256": real}, "route": "recomputed", "usable": True, "stages": ["E9.5", "E14.5"],
            "outputs": {stats: {"sha256": fre.N54_STATS_SHA256}},
            "sets": {n: {"found": ["g"], "controls": ["h"]} for n in fr.load_config()["sets"]}}
    head = {fr.CONFIG: fr.CONFIG.read_bytes(), fr.RECORD: b"{}", fre.RECORD: None}
    monkeypatch.setattr(fr, "committed_at_head", lambda p: head.get(p))
    with pytest.raises(fre.Refused, match="route recomputed"):
        fre.preflight("N70", dst)
    head[fr.RECORD] = json.dumps(good).encode()
    monkeypatch.setattr(fre, "x002_names_this_use", lambda: False)
    with pytest.raises(fre.Refused, match="X002"):
        fre.preflight("N70", dst)
    monkeypatch.setattr(fre, "x002_names_this_use", lambda: True)
    monkeypatch.setattr(fre, "row_line", lambda row: "| N70 | T1 | another script |")
    with pytest.raises(fre.Refused, match="naming feng_relabel_e18.py"):
        fre.preflight("N70", dst)
    monkeypatch.setattr(fre, "row_line", lambda row: "| N70 | T1 | feng_relabel_e18.py |")
    from vec import provenance

    monkeypatch.setattr(provenance, "git_state", lambda repo: ("c" * 40, ["x.py"]))
    with pytest.raises(fre.Refused, match="not clean"):
        fre.preflight("N70", dst)
    monkeypatch.setattr(provenance, "git_state", lambda repo: ("c" * 40, []))
    monkeypatch.setattr(fre, "RECORD", tmp_path / "record" / "MANIFEST.json")
    head[fre.RECORD] = None
    dst.mkdir(parents=True)
    with pytest.raises(fre.Refused, match="already recorded or was started"):
        fre.preflight("N70", dst)
    gate = fre.preflight("N70", tmp_path / "fresh")
    assert gate["row"] == "N70" and gate["clean"] is True and gate["config_sha256"] == real
    assert gate["n54_record_sha256"] == fre.sha256_bytes(head[fr.RECORD])


# ---------------------------------------------------------------------------------------------------------------------
# feng_dirs: the E9.5 -> E18 chord

LABELS = {"r_cm": 30, "r_endocardial_ec": 20, "r_vascular_ec": 30, "r_epicardial": 12}
FGENES = [f"f{i}" for i in range(30)]
LEVEL = {s: np.random.default_rng(i).gamma(0.6, 20.0, (len(LABELS), len(FGENES)))
         for i, s in enumerate(("E9.5", "E14.5", "E18"))}


def _rows(stages, cells_of=lambda s, lab: LABELS[lab]):
    keys, n, cp = [], [], []
    for s in stages:
        for j, lab in enumerate(LABELS):
            k = cells_of(s, lab)
            if k == 0:
                continue
            for half in "AB":
                keys.append(f"{s}|LV|{lab}|{half}")
                n.append(float(k))
                cp.append(k * LEVEL[s][j])
    cp = np.array(cp)
    return {"keys": np.array(keys), "genes": np.array(FGENES), "n": np.array(n), "S_cp": cp,
            "S_log": np.zeros_like(cp), "S_nz": np.zeros_like(cp), "depth_sum": np.array(n) * 1e3}


def _chord_world(tmp_path, monkeypatch, e18_cells=lambda s, lab: LABELS[lab], tamper=None, **man_over):
    n54_dir, e18_dir = tmp_path / "n54", tmp_path / "e18"
    n54_dir.mkdir(parents=True)
    e18_dir.mkdir(parents=True)
    n54 = _rows(("E9.5", "E14.5"))
    n54_path = n54_dir / "cd1_stats_relabel.npz"
    np.savez(n54_path, **n54)
    fio.write_manifest(n54_dir, {"source": "feng_relabel", "config": {"sha256": fio.sha256_file(fr.CONFIG),
                                                                      "status": "frozen"},
                                 "unsourced": {}, "gate": dict(GATE), "provenance": {"git_uncommitted": []},
                                 "outputs": {n54_path.name: {"sha256": fio.sha256_file(n54_path), "usable": True}}})
    comb = fre.combine(fre.load_raw(n54_path), _rows(("E18",), e18_cells))
    if tamper is not None:
        tamper(comb)
    e18_path = e18_dir / fre.STATS
    np.savez(e18_path, **comb)
    man = {"source": fre.SOURCE, "config": {"sha256": fre.YAML_SHA256, "status": "frozen"}, "stages": ["E18"],
           "stages_copied": ["E9.5", "E14.5"], "n54_stats": {"sha256": fio.sha256_file(n54_path)},
           "copied_rows": {"identical": True}, "authors_labels_read": False, "usable": True,
           "gate": {"row": "N70", "clean": True}, "provenance": {"git_uncommitted": []},
           "outputs": {fre.STATS: {"sha256": fio.sha256_file(e18_path), "usable": True}}, **man_over}
    fio.write_manifest(e18_dir, man)
    monkeypatch.setattr(feng_dirs, "record_matches", lambda p: True)
    monkeypatch.setattr(feng_dirs, "e18_record_matches", lambda p: True)
    monkeypatch.setattr(feng_dirs, "recomputed_stats_path", lambda *a: n54_path)
    monkeypatch.setattr(feng_dirs, "e18_stats_path", lambda: e18_path)
    monkeypatch.setattr(fre, "N54_STATS_SHA256", fio.sha256_file(n54_path))
    return n54_path, e18_path


def _want(stage_a, stage_b, labels):
    a = sum(LABELS[lab] * 2 * LEVEL[stage_a][list(LABELS).index(lab)] for lab in labels) / sum(
        LABELS[lab] * 2 for lab in labels)
    b = sum(LABELS[lab] * 2 * LEVEL[stage_b][list(LABELS).index(lab)] for lab in labels) / sum(
        LABELS[lab] * 2 for lab in labels)
    return np.log((b + EPS) / (a + EPS))


def test_the_chord_is_the_e95_to_e18_log_ratio_on_the_recomputed_groups(tmp_path, monkeypatch):
    _chord_world(tmp_path, monkeypatch)
    genes = [*FGENES[:20], "absent"]
    d18, rec = feng_dirs.recomputed_e18_directions(genes)
    assert set(d18) == set(GROUPS4)
    assert np.allclose(d18["CM"][:20], _want("E9.5", "E18", ["r_cm"])[:20]) and np.isnan(d18["CM"][20])
    assert np.allclose(d18["Endothelium"][:20], _want("E9.5", "E18", ["r_endocardial_ec", "r_vascular_ec"])[:20])
    assert rec["bracket"] == ["E9.5", "E18"] and rec["bracket_days"] == [9.5, 18.0]
    assert rec["route"] == "recomputed" and rec["e18"]["copied_rows_identical"] is True
    assert rec["e18"]["groups"]["CM"]["E18"] == 2 * LABELS["r_cm"]
    late = feng_dirs.e18_late_term(genes)
    for g in GROUPS4:  # D18 = F + G in log ratios (to rounding)
        assert np.allclose(late["D18"][g], late["F"][g] + late["G"][g], equal_nan=True, rtol=0, atol=1e-12)
        assert late["D18"][g].tobytes() == d18[g].tobytes()
    assert np.allclose(late["basis"]["CM"][:20], ((LEVEL["E14.5"][0] + LEVEL["E18"][0]) / 2)[:20])
    assert feng_dirs.bracket_fraction(10.5, 12.5, feng_dirs.E18_BRACKET_DAYS) == pytest.approx(2 / 8.5, abs=1e-15)


@pytest.mark.parametrize("case, msg", [
    ("no_e18", "hold no stage"),
    ("small_group", "not row R13's"),
    ("tampered", "not row N54's byte for byte"),
    ("source", "manifest source"),
    ("unusable", "usable not true"),
    ("not_identical", "not identical"),
    ("gate", "no clean gate"),
    ("authors", "authors_labels_read"),
    ("stages", "stages"),
])
def test_the_chord_refuses(tmp_path, monkeypatch, case, msg):
    kw = {"no_e18": {"e18_cells": lambda s, lab: 0},
          "small_group": {"e18_cells": lambda s, lab: 3 if lab == "r_epicardial" else LABELS[lab]},
          "tampered": {"tamper": lambda c: c["S_cp"].__setitem__((0, 0), c["S_cp"][0, 0] + 1.0)},
          "source": {"source": "feng_relabel_late"}, "unusable": {"usable": False},
          "not_identical": {"copied_rows": {"identical": False}}, "gate": {"gate": {"row": "T69", "clean": True}},
          "authors": {"authors_labels_read": True}, "stages": {"stages": ["E17.5"]}}[case]
    if case == "no_e18":
        with pytest.raises(fre.Refused):  # no E18 row at all: the labeller's own combine refuses an empty table
            _chord_world(tmp_path, monkeypatch, **kw)
        n54_path, e18_path = _chord_world(tmp_path / "x", monkeypatch)
        st = fre.load_raw(n54_path)
        np.savez(e18_path, **st)  # the E18 file holding no E18 row
        man = json.loads((e18_path.parent / "MANIFEST.json").read_text())
        man["outputs"][fre.STATS]["sha256"] = fio.sha256_file(e18_path)
        fio.write_manifest(e18_path.parent, man)
    else:
        _chord_world(tmp_path, monkeypatch, **kw)
    with pytest.raises(ValueError, match=msg):
        feng_dirs.recomputed_e18_directions(FGENES)


def test_the_chord_refuses_a_group_that_fell_back_and_a_record_not_committed(tmp_path, monkeypatch):
    _chord_world(tmp_path, monkeypatch)
    real = feng_dirs.recomputed_directions

    def fell(*a, **k):
        d, ext = real(*a, **k)
        ext = {**ext, "fallbacks": {**ext["fallbacks"], "Endocardium": {"own": False, "fired": True,
                                                                         "takes": "Endothelium"}}}
        return d, ext

    monkeypatch.setattr(feng_dirs, "recomputed_directions", fell)
    with pytest.raises(ValueError, match="nothing stands in"):
        feng_dirs.recomputed_e18_directions(FGENES)
    monkeypatch.setattr(feng_dirs, "recomputed_directions", real)
    monkeypatch.setattr(feng_dirs, "e18_record_matches", lambda p: False)
    with pytest.raises(ValueError, match="RECORD"):
        feng_dirs.recomputed_e18_directions(FGENES)


# ---------------------------------------------------------------------------------------------------------------------
# Row T69's read: the directions and the arms


def _d(n=400, seed=1, nan_every=13):
    d = np.random.default_rng(seed).normal(0.0, 0.5, n)
    d[::nan_every] = np.nan
    return d


def test_the_size_control_has_rc18s_applied_size_on_the_moved_genes(rr18):
    genes = [f"g{i}" for i in range(400)]
    f, d18 = _d(seed=1), _d(seed=2) * 1.9
    held = ["g3", "g7", "g200"]
    mask = rr18.moved_mask(genes, held, f, d18)
    assert not mask[3] and not mask[0] and mask[1]  # held and NaN genes do not move
    sm = rr18.size_control(f, d18, mask)
    assert rr18.l2(sm / 8.5, mask) == pytest.approx(rr18.l2(d18 / 8.5, mask), rel=1e-12)
    assert np.isnan(sm[0]) and np.allclose(sm[mask] / f[mask], sm[1] / f[1])  # one factor per group
    with pytest.raises(ValueError, match="no size"):
        rr18.size_control(np.zeros(400), d18, mask)


def test_the_permuted_control_is_row_t15s_procedure(rr18):
    g, basis = _d(seed=3), np.random.default_rng(4).gamma(1.0, 5.0, 400)
    basis[5] = np.nan
    p = rr18.permuted_late(g, basis, "Endocardium")
    ok = np.isfinite(g) & np.isfinite(basis)
    want = harness.t3_decile_permutation(g[ok], basis[ok], [20261110, harness.LP_GROUPS.index("Endocardium")], 10)
    assert p[ok].tobytes() == want.tobytes() and p[5] == g[5] and np.isnan(p[0])
    assert sorted(p[ok]) == sorted(g[ok]) and not np.array_equal(p[ok], g[ok])
    assert not np.array_equal(rr18.permuted_late(g, basis, "CM")[ok], p[ok])  # each group its own seed
    with pytest.raises(ValueError, match="R2's guard"):
        rr18.permuted_late(g, basis, "Blood")


def _late(n, groups=GROUPS4):
    rng = np.random.default_rng(11)
    F = {g: _d(n, seed=20 + i, nan_every=17) for i, g in enumerate(groups)}
    G = {g: np.where(np.isfinite(F[g]), rng.normal(0, 0.4, n), np.nan) for g in groups}
    return {"F": F, "G": G, "D18": {g: F[g] + G[g] for g in groups},
            "basis": {g: rng.gamma(1.0, 5.0, n) for g in groups}}


def test_the_arm_directions(rr18):
    genes = [f"g{i}" for i in range(300)]
    late = _late(300)
    dirs, sizes = rr18.arm_directions(late, genes, ["g1"])
    assert set(dirs) == set(rr18.ARMS) and all(set(v) == set(GROUPS4) for v in dirs.values())
    assert dirs["m1rc_gmom"]["CM"].tobytes() == late["F"]["CM"].tobytes()
    assert dirs["rc18"]["CM"].tobytes() == late["D18"]["CM"].tobytes()
    for g in GROUPS4:
        assert sizes[g]["rc18_sm"] == pytest.approx(sizes[g]["rc18"], rel=1e-12)
        assert sizes[g]["m1rc_gmom"] == pytest.approx(0.2 * rr18.l2(late["F"][g], rr18.moved_mask(
            genes, ["g1"], late["F"][g], late["D18"][g])), rel=1e-12)
    assert rr18.FRACTIONS == {"m1rc_gmom": 0.2, "rc18": 1 / 8.5, "rc18_sm": 1 / 8.5, "rc18_perm": 1 / 8.5}
    with pytest.raises(ValueError, match="groups"):
        rr18.arm_directions({**late, "G": {"CM": late["G"]["CM"]}}, genes, [])


def _held_m1(genes):
    return sorted({*HELD, genes[3], genes[7]} & set(genes))


def test_the_cut_arms_are_t62s_candidate_build_on_each_arms_direction_and_fraction(rr18):
    start, ref = classified_pair(seed=5, n_scale=0.5)
    ctx = _gmom_ctx(ref, start, seed=1000003)
    fits = rr18.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    genes = list(ref.var_names)
    held_m1 = _held_m1(genes)
    dirs, _ = rr18.arm_directions(_late(ref.n_vars), genes, held_m1)
    arms = rr18.cut_arms(ctx, fits["file"], dirs, held_m1)
    assert tuple(arms) == rr18.ARMS
    for a in rr18.ARMS:
        alone = rr18.rc.cut_pair(ctx, fits["file"], dirs[a], rr18.FRACTIONS[a], held_m1)["m1rc_gmom"]
        assert np.asarray(arms[a].X).tobytes() == np.asarray(alone.X).tobytes(), a
    x = {a: np.asarray(v.X).tobytes() for a, v in arms.items()}
    assert len(set(x.values())) == 4


def test_the_library_arms(rr18):
    e95, e85 = classified_pair(seed=5, n_scale=0.5)[::-1]
    spec = load_spec()
    fb = {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}
    within = {"surrounding": ["Foregut"]}
    factors = rr18.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    shares = {"cm": 0.5, "endocardium": 0.2, "epicardium": 0.15, "blood": 0.05, "surrounding": 0.1}
    genes = list(e95.var_names)
    held_m1 = _held_m1(genes)
    dirs, _ = rr18.arm_directions(_late(e95.n_vars), genes, held_m1)
    arms = rr18.library_arms(e95, shares, factors, HELD, held_m1, 6, within, fb, dirs, n_cells=300)
    for a in rr18.ARMS:
        alone = rr18.rc.library_pair(e95, shares, factors, HELD, held_m1, 6, within, fb, dirs[a],
                                     rr18.FRACTIONS[a], n_cells=300)["m1rc_gmom"]
        got, want = counts_of(arms[a]), counts_of(alone)
        got = np.asarray(got.todense()) if hasattr(got, "todense") else np.asarray(got)
        want = np.asarray(want.todense()) if hasattr(want, "todense") else np.asarray(want)
        assert np.array_equal(got, want), a


# ---------------------------------------------------------------------------------------------------------------------
# Row T69's read: the statement on hand-made tables

WEIGHTS = {"de_recovery": 0.25, "direction": 0.25, "distribution": 0.30, "covariance": 0.20}
TR_YAML = "883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49"
TR_STATS = "36d1f7edcbbf3d7a80c68b51d5624a94bb753668079fe3d69aed720c71baa7e9"
E18_SHA = "e" * 64
COMMITTED = {k: {s: 50.0 + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}
E18_OK = {"stats_sha256": E18_SHA, "reasons": []}


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in WEIGHTS}}


def _task(key, gain=0.5, sm=0.1, perm=0.1, gate=0.0, seeds=range(10), commit="c" * 40, clean=True, identical=True,
          fractions=None, e18_sha=E18_SHA, inverted=False, per_seed=None):
    rows = []
    for s in seeds:
        base = COMMITTED[key][s] + gate
        g = gain if per_seed is None else per_seed[s]
        arms = {"m1rc_gmom": _arm(base), "rc18": _arm(base + g), "rc18_sm": _arm(base + sm),
                "rc18_perm": _arm(base + perm)}
        if key in ("1", "2", "3"):
            arms["copy_last"] = _arm(base - 1)
        row = {"seed": s, "arms": arms, "sha256": {a: "s" for a in arms}}
        if inverted:
            row["inverted_anchors"] = ["x"]
        rows.append(row)
    return {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean},
            "fractions": fractions or {"m1rc_gmom": 0.2, "rc18": 1 / 8.5, "rc18_sm": 1 / 8.5, "rc18_perm": 1 / 8.5},
            "kit": {"copied_rows_identical": identical, "e18_sha256": e18_sha, "sizes": {}},
            "relabel": {"sha256": TR_STATS, "stats_sha256_hashed": TR_STATS, "yaml_sha256": TR_YAML,
                        "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False},
            "e18": {"stats_sha256": e18_sha, "stats_sha256_hashed": e18_sha, "yaml_sha256": TR_YAML,
                    "authors_labels_read": False}}


def _loaded(**over):
    out = {k: _task(k) for k in ("chamber", "whole", "1", "2", "3")}
    out.update(over)
    return out


def _read(rr18, loaded, committed=COMMITTED, checks=None, e18=E18_OK, source=()):
    return rr18.t69_reading(loaded, committed, checks or {"reasons": []}, e18, source)


def test_the_constants(rr18):
    assert rr18.ROW == "T69" and rr18.TRIES_COUNT == 15 and rr18.GATE_TOL == 0.02
    assert rr18.COUNTED == ("rc18", "rc18_sm", "rc18_perm") and rr18.RC18_FRACTION == pytest.approx(1 / 8.5)
    assert rr18.PERM_SEED == 20261110 and rr18.PERM_BINS == 10 and rr18.MIN_SEEDS == 8
    assert rr18.CHAMBER_BAR == 0.3 and rr18.WHOLE_BAR == 0.0 and rr18.LIB_FLOOR == -0.3 and rr18.STOP_FLOOR == -0.3
    assert rr18.STOP_TASK == "3" and rr18.MIN_LIBRARIES == 2
    assert rr18.rc.STATS_SHA256 == TR_STATS and rr18.rc.YAML_SHA256 == TR_YAML == fre.YAML_SHA256
    assert fre.N54_STATS_SHA256 == TR_STATS


def test_the_default_table_is_ahead(rr18):
    r = _read(rr18, _loaded())
    assert r["outcome"] == "AHEAD" and r["ahead"] and not r["stop"] and all(r["parts"].values())
    assert r["j6_range"] == [0.0, pytest.approx(0.5)]
    assert r["chamber_beats_both"]["both"] == 10


@pytest.mark.parametrize("over, part", [
    ({"chamber": _task("chamber", gain=0.3)}, "chamber"),  # above +0.3, strictly
    ({"chamber": _task("chamber", per_seed=[1.0] * 7 + [-0.1] * 3)}, "chamber"),  # 7 of 10 ahead
    ({"whole": _task("whole", gain=0.0)}, "whole"),  # above 0, strictly
    ({"1": _task("1", inverted=True), "2": _task("2", inverted=True)}, "two_valid_libraries"),
    ({"2": _task("2", gain=-0.31)}, "no_valid_library_below_floor"),
    ({"chamber": _task("chamber", sm=0.6)}, "beats_both_controls"),
    ({"chamber": _task("chamber", perm=0.6)}, "beats_both_controls"),
])
def test_every_failed_part_reads_not_ahead(rr18, over, part):
    r = _read(rr18, _loaded(**over))
    assert r["outcome"].startswith("NOT AHEAD") and not r["parts"][part] and not r["ahead"]


def test_beating_both_controls_counts_the_same_seed(rr18):
    t = _task("chamber")
    for x in t["rows"][:3]:  # three seeds below rc18_sm, three other seeds below rc18_perm: 4 seeds beat both
        x["arms"]["rc18_sm"] = _arm(x["arms"]["rc18"]["score"] + 0.1)
    for x in t["rows"][3:6]:
        x["arms"]["rc18_perm"] = _arm(x["arms"]["rc18"]["score"] + 0.1)
    r = _read(rr18, _loaded(chamber=t))
    assert r["chamber_beats_both"] == {"both": 4, "over_sm": 7, "over_perm": 7, "of": 10}
    assert r["outcome"].startswith("NOT AHEAD")


def test_the_stop_ends_the_line_even_on_an_invalid_library_3(rr18):
    r = _read(rr18, _loaded(**{"3": _task("3", gain=-0.31, inverted=True)}))
    assert r["stop"] and r["outcome"].startswith("STOP") and not r["ahead"]
    assert all(r["parts"].values())  # every part of AHEAD holds; the stop alone ends it
    assert not _read(rr18, _loaded(**{"3": _task("3", gain=-0.3)}))["stop"]


@pytest.mark.parametrize("case", ["gate", "identical", "fraction", "seeds", "commits", "clean", "e18_sha", "kit_sha",
                                  "missing", "checks", "source", "stats", "route"])
def test_every_void_condition(rr18, case):
    loaded, committed, checks, e18, source = _loaded(), COMMITTED, {"reasons": []}, E18_OK, ()
    if case == "gate":
        loaded["whole"] = _task("whole", gate=0.03)
    elif case == "identical":
        loaded["1"] = _task("1", identical=False)
    elif case == "fraction":
        loaded["2"] = _task("2", fractions={"m1rc_gmom": 0.2, "rc18": 0.2, "rc18_sm": 1 / 8.5, "rc18_perm": 1 / 8.5})
    elif case == "seeds":
        loaded["3"] = _task("3", seeds=range(9))
    elif case == "commits":
        loaded["3"] = _task("3", commit="d" * 40)
    elif case == "clean":
        loaded["chamber"] = _task("chamber", clean=False)
    elif case == "e18_sha":
        loaded["chamber"] = _task("chamber", e18_sha="f" * 64)
    elif case == "kit_sha":
        loaded["chamber"]["kit"]["e18_sha256"] = "f" * 64
    elif case == "missing":
        loaded["2"] = None
    elif case == "checks":
        e18 = {"stats_sha256": E18_SHA, "reasons": ["the RECORD is not at HEAD"]}
    elif case == "source":
        source = ["the source names x"]
    elif case == "stats":
        loaded["1"]["relabel"]["sha256"] = "0" * 64
    elif case == "route":
        loaded["1"]["relabel"]["route"] = "none"
    r = _read(rr18, loaded, committed, checks, e18, source)
    assert r["outcome"] == "VOID" and r["reasons"], case


def test_a_gate_within_the_tolerance_passes(rr18):
    assert _read(rr18, _loaded(whole=_task("whole", gate=0.019)))["outcome"] == "AHEAD"


def test_the_read_writes_its_files_once_and_prints_outcome_lines_only(rr18, tmp_path, monkeypatch, capsys):
    out = rr18.RESULTS / "t69_rc18"
    (out / "tasks").mkdir(parents=True)
    for k, name in rr18.task_files().items():
        (out / "tasks" / name).write_text(json.dumps(_task(k)))
    monkeypatch.setattr(rr18.br, "tries_rows", lambda path: {"T69"})
    key_of = {v: k for k, v in rr18.task_files().items()}
    monkeypatch.setattr(rr18, "committed_t62", lambda name: COMMITTED[key_of[name]])
    monkeypatch.setattr(rr18.rc, "relabel_checks", lambda need_stats=True: {"reasons": []})
    monkeypatch.setattr(rr18, "e18_checks", lambda need_stats=True: E18_OK)
    monkeypatch.setattr(rr18.gv, "stamp", lambda s, inputs: {"commit": "c" * 40, "tree_clean": True, "utc": "now"})
    monkeypatch.setattr(rr18.gv, "tree_state", lambda: {"commit": "c" * 40, "clean": True})
    res = rr18.read(out)
    printed = capsys.readouterr().out.splitlines()
    assert res["reading"]["outcome"] == "AHEAD" and printed[0] == "Outcome: AHEAD"
    assert all(re.match(r"^(Outcome: |Parts: |What follows: |\(numbers in )", ln) for ln in printed)
    assert (out / "t69.json").exists() and "Tries counted: 15" in (out / "t69.txt").read_text()
    with pytest.raises(SystemExit, match="already holds"):
        rr18.read(out)


def test_the_runs_and_the_read_are_refused_without_the_row(rr18, monkeypatch, tmp_path):
    monkeypatch.setattr(rr18.br, "tries_rows", lambda path: {"T67"})
    for argv in (["run-cut", "--cut", "chamber", "--out", "x"], ["run-lib", "--library", "1", "--out", "x"],
                 ["read", "--out", "x"]):
        with pytest.raises(SystemExit, match="T69 is refused"):
            rr18.main(argv)


def test_a_run_is_refused_in_order(rr18, monkeypatch, tmp_path):
    out = tmp_path / "tasks"
    monkeypatch.setattr(rr18.br, "tries_rows", lambda path: {"T69"})
    monkeypatch.setattr(rr18.gv, "tree_state", lambda: {"commit": "c" * 40, "clean": False})
    with pytest.raises(SystemExit, match="not clean"):
        rr18.refuse(out, ("lib1.json",))
    monkeypatch.setattr(rr18.gv, "tree_state", lambda: {"commit": "c" * 40, "clean": True})
    monkeypatch.setattr(rr18, "source_guard", lambda: ["bad"])
    with pytest.raises(SystemExit, match="bad"):
        rr18.refuse(out, ("lib1.json",))
    monkeypatch.setattr(rr18, "source_guard", lambda: [])
    monkeypatch.setattr(rr18.rc, "source_guard", lambda: [])
    monkeypatch.setattr(rr18.rc, "relabel_checks", lambda: {"passes": False, "reasons": ["n54 off"]})
    with pytest.raises(SystemExit, match="n54 off"):
        rr18.refuse(out, ("lib1.json",))
    monkeypatch.setattr(rr18.rc, "relabel_checks", lambda: {"passes": True, "reasons": []})
    monkeypatch.setattr(rr18, "e18_checks", lambda: {"passes": False, "reasons": ["e18 off"]})
    with pytest.raises(SystemExit, match="e18 off"):
        rr18.refuse(out, ("lib1.json",))
    monkeypatch.setattr(rr18, "e18_checks", lambda: {"passes": True, "reasons": []})
    monkeypatch.setattr(rr18.fre, "x002_names_this_use", lambda: False)
    with pytest.raises(SystemExit, match="X002"):
        rr18.refuse(out, ("lib1.json",))
    monkeypatch.setattr(rr18.fre, "x002_names_this_use", lambda: True)
    assert rr18.refuse(out, ("lib1.json",))[0]["clean"] is True
    out.mkdir()
    (out / "lib1.json").write_text("{}")
    with pytest.raises(SystemExit, match="already holds"):
        rr18.refuse(out, ("lib1.json",))


def _e18_files(root, **change):
    root.mkdir(parents=True, exist_ok=True)
    stats = root / fre.STATS
    stats.write_bytes(b"e18 statistics")
    sha = fio.sha256_file(stats)
    man = {"source": "feng_relabel_e18", "config": {"sha256": TR_YAML}, "authors_labels_read": False,
           "stages": ["E18"], "stages_copied": ["E9.5", "E14.5"], "n54_stats": {"sha256": TR_STATS},
           "copied_rows": {"identical": True}, "usable": True, "gate": {"row": "N70", "clean": True},
           "outputs": {fre.STATS: {"sha256": sha}}, **change}
    path = root / "MANIFEST.json"
    path.write_text(json.dumps(man))
    return path, stats, sha


def test_the_e18_checks(rr18, tmp_path):
    path, stats, sha = _e18_files(tmp_path / "ok")

    def run(**kw):
        args = {"manifest": path, "stats_path": stats, "line": f"| T69 (x) | T1 | sha256 {sha} |",
                "record_at_head": path.read_bytes(), **kw}
        return rr18.e18_checks(**args)

    ok = run()
    assert ok["passes"] and ok["stats_sha256"] == sha == ok["stats_sha256_hashed"] and ok["e18_row"] == "N70"
    assert not run(line="| T69 (x) | T1 | sha256 <the N58 statistics sha256> |")["passes"]
    assert not run(record_at_head=b"{}")["passes"]
    assert not run(stats_path=tmp_path / "absent.npz")["passes"]
    assert run(stats_path=tmp_path / "absent.npz", need_stats=False)["passes"]
    for k, change in enumerate(({"source": "feng_relabel_late"}, {"config": {"sha256": "0" * 64}},
                                {"authors_labels_read": True}, {"stages": ["E17.5"]}, {"stages_copied": ["E9.5"]},
                                {"n54_stats": {"sha256": "0" * 64}}, {"copied_rows": {"identical": False}},
                                {"usable": False}, {"stopped": True}, {"gate": {"row": "T69", "clean": True}},
                                {"gate": {"row": "N70", "clean": False}})):
        p2, s2, sha2 = _e18_files(tmp_path / f"c{k}", **change)
        r = rr18.e18_checks(manifest=p2, stats_path=s2, line=f"| T69 | {sha2} |", record_at_head=p2.read_bytes())
        assert not r["passes"], change


def test_the_source_guard_passes_this_read_and_refuses_other_routes(rr18, tmp_path):
    assert rr18.source_guard() == [] and rr18.rc.source_guard() == []
    for text in ("import perm_read as pr\n", 'x = fd.feng_directions(list(g))\n', 'g = "auth' + 'ors"\n',
                 "x = __imp" + "ort__('eb_read')\n"):
        bad = tmp_path / "bad.py"
        bad.write_text(text)
        assert rr18.source_guard(bad), text


# ---------------------------------------------------------------------------------------------------------------------
# The read's launcher, experiments/t1_m1/rc18.sh


def test_the_launcher_deletes_nothing_and_holds_the_rules():
    text = LAUNCHER.read_text()
    code = [line.split("#")[0] for line in text.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    for msg in re.findall(r"\$\{[^}]*:\?([^}]*)\}", text):
        assert "'" not in msg and '"' not in msg, msg
    usage = text[text.index("usage() {"):text.index("}", text.index("usage() {"))]
    assert "'" not in usage.split('echo "', 1)[1].split('"', 1)[0]
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0" in text
    assert 'export UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in text and "python -I -B" in text
    assert "--basetemp=" in text and 'S="$HOME/v/scratch/t69"' in text
    assert "/v/vec-campaign" not in text.split("set -u", 1)[1]  # no checkout path in the code
    assert text.index('echo running > "$S/$NAME.exit"') < text.index('echo "$ROWC" > "$S/$NAME.head"')
    assert text.index("refuse 98") < text.index("echo running")  # every check before the running mark
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", text)
    assert "tests/test_t1_rc18.py" in listed and all((ROOT / f).exists() for f in listed)
    n54 = (ROOT / "experiments" / "external" / "n54_step.sh").read_text()
    assert "tests/test_t1_rc18.py" in n54.split("  tests)")[1].split(";;")[0]


def _env(tmp_path, fake_uv=7):
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    fake = tmp_path / "bin"
    fake.mkdir(exist_ok=True)
    (fake / "uv").write_text(f"#!/bin/sh\nexit {fake_uv}\n")  # a step that runs gets this exit; nothing real runs
    (fake / "uv").chmod(0o755)
    env["PATH"] = f"{fake}:{env.get('PATH', '')}"
    return env


def _repo(tmp_path, rows=True):
    repo = tmp_path / "repo"
    (repo / "experiments" / "t1_m1").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    shutil.copy(LAUNCHER, repo / "experiments" / "t1_m1" / "rc18.sh")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(
        "| T62 | T1 | x |\n" + ("| T69 (one row) | T1 | `experiments/t1_m1/rc18_read.py` |\n"
                                 "| N70 (target-free) | T1 | `experiments/external/feng_relabel_e18.py` |\n"
                                 if rows else ""))
    (repo / ".gitignore").write_text("data/\n")
    env = _env(tmp_path)
    for args in (["init", "-q"], ["add", "-A"], ["-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "x"]):
        subprocess.run(["git", *args], cwd=repo, env=env, check=True, capture_output=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, capture_output=True, text=True)
    return repo, head.stdout.strip()


def _launch(tmp_path, repo, *argv):
    return subprocess.run(["bash", "experiments/t1_m1/rc18.sh", *argv], cwd=repo, env=_env(tmp_path),
                          capture_output=True, text=True).returncode


def _s(tmp_path):
    return tmp_path / "home" / "v" / "scratch" / "t69"


def _read_exit(tmp_path, name):
    p = _s(tmp_path) / f"{name}.exit"
    return p.read_text().strip() if p.exists() else None


def _refused(tmp_path, name):
    return (_s(tmp_path) / f"{name}.refused.exit").read_text().strip()


def _mark(tmp_path, name, code, head):
    s = _s(tmp_path)
    s.mkdir(parents=True, exist_ok=True)
    (s / f"{name}.exit").write_text(f"{code}\n")
    (s / f"{name}.head").write_text(f"{head}\n")


def test_the_launchers_refusal_paths_never_touch_a_step_exit(tmp_path):
    repo, head = _repo(tmp_path)
    assert _launch(tmp_path, repo) == 2 and _refused(tmp_path, "usage") == "2"
    assert _launch(tmp_path, repo, "tests", head[:12]) == 2  # a prefix is never enough
    assert _launch(tmp_path, repo, "label", head, "T69") == 2  # the labelling takes an N row
    assert _launch(tmp_path, repo, "build", head) == 97 and _refused(tmp_path, "unknown") == "97"
    assert _launch(tmp_path, repo, "cut", "half", head) == 97
    assert _launch(tmp_path, repo, "tests", "0" * 40) == 95 and _refused(tmp_path, "tests") == "95"
    assert _read_exit(tmp_path, "tests") is None  # a refusal writes no step exit
    assert _launch(tmp_path, repo, "label", head, "N71") == 96  # no row N71 names the labeller
    assert _launch(tmp_path, repo, "label", head, "N70") == 92  # label-tests has not exited 0 here
    assert _launch(tmp_path, repo, "tests", head) == 92  # no labelling yet
    for step in (["cut", "chamber"], ["cut", "whole"], ["lib", "1"], ["lib", "2"], ["lib", "3"], ["read"]):
        assert _launch(tmp_path, repo, *step, head) == 92, step
    # a step that passes its checks runs: running, then its head, then its exit (the fake uv exits 7)
    assert _launch(tmp_path, repo, "label-tests", head, "N70") == 7
    assert _read_exit(tmp_path, "label_tests") == "7"
    assert (_s(tmp_path) / "label_tests.head").read_text().strip() == head
    _mark(tmp_path, "label_tests", 0, head)
    assert _launch(tmp_path, repo, "label-tests", head, "N70") == 94  # already exited 0 here
    assert _read_exit(tmp_path, "label_tests") == "0"  # the refusal left the finished exit alone
    assert _launch(tmp_path, repo, "label", head, "N70") == 93  # N54's statistics are not in place
    data = repo / "data" / "derived"
    (data / "feng_relabel_c").mkdir(parents=True)
    (data / "feng_relabel_c" / "cd1_stats_relabel.npz").write_text("x")
    (data / "feng_relabel_e18").mkdir()
    assert _launch(tmp_path, repo, "label", head, "N70") == 93  # its output folder exists
    _mark(tmp_path, "label", 0, "f" * 40)  # exited 0 at a commit this one does not descend from
    assert _launch(tmp_path, repo, "tests", head) == 92
    _mark(tmp_path, "label", 0, head)
    assert _launch(tmp_path, repo, "tests", head) == 93  # the E18 statistics are not in place
    (data / "feng_relabel_e18" / "cd1_stats_relabel_e18.npz").write_text("x")
    assert _launch(tmp_path, repo, "tests", head) == 7 and _read_exit(tmp_path, "tests") == "7"
    assert _launch(tmp_path, repo, "lib", "2", head) == 92  # tests did not exit 0
    _mark(tmp_path, "tests", 0, head)
    tasks = repo / "experiments" / "t1_m1" / "results" / "t69_rc18" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / "lib2.json").write_text("{}")  # untracked in the row's results folder: the tree stays clean
    assert _launch(tmp_path, repo, "lib", "2", head) == 91 and _refused(tmp_path, "lib_2") == "91"
    assert _read_exit(tmp_path, "lib_2") is None
    for name in ("cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3"):
        _mark(tmp_path, name, 0, head)
    assert _launch(tmp_path, repo, "read", head) == 90  # task files missing
    for f in ("T1_E9.5_E10.5_chamber.json", "T1_E9.5_E10.5.json", "lib1.json", "lib3.json"):
        (tasks / f).write_text("{}")
    (tasks.parent / "t69.txt").write_text("x")
    assert _launch(tmp_path, repo, "read", head) == 91
    (repo / "stray.txt").write_text("x")  # an untracked file outside the results folder
    assert _launch(tmp_path, repo, "read", head) == 98 and _refused(tmp_path, "read") == "98"


def test_the_launcher_refuses_without_the_rows(tmp_path):
    repo, head = _repo(tmp_path, rows=False)
    assert _launch(tmp_path, repo, "tests", head) == 96 and _refused(tmp_path, "tests") == "96"
    assert _launch(tmp_path, repo, "label-tests", head, "N70") == 96 and _refused(tmp_path, "label_tests") == "96"
