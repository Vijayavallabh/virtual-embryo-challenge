"""Rows N64 and T72 (experiments/harness/TRIES.md): experiments/t1_scope/late_mix.py (step 0 of row R54's late-mix
line) and experiments/t1_scope/t72_late_mix.py (its practice read) on synthetic stages, counts and scenario files only.
No released file, board or score is read; the committed relabelling manifests are read for their label names and their
frozen-record fields only (target-free), and the committed T30 scenario folder for nothing."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "t1_scope"))
import kind_gap as kg  # noqa: E402
import late_mix as lm  # noqa: E402
import t30_pricing as t30  # noqa: E402
import t72_late_mix as t72  # noqa: E402

from vec import harness, t1_scope  # noqa: E402
from vec.core import log1p_cp10k, quotas  # noqa: E402

C = "c" * 40
LABELS = {"V-CM": "cm", "OFT-CM": "cm", "Endocardium-1": "endocardium", "Early VIC": "endocardium",
          "Epicardium": "epicardium", "BEC": "vascular_endothelium", "Blood": "blood", "Macrophage": "blood",
          "Hepatocyte": "surrounding", "aSHF": "progenitor"}


def _stage(sizes: dict[str, dict[str, int]], seed: int = 3) -> ad.AnnData:
    """A synthetic E10.5: per library (its cell-name suffix) the cells of each label."""
    rng = np.random.default_rng(seed)
    labels, names = [], []
    for lib, per in sizes.items():
        for lab, n in per.items():
            labels += [lab] * n
            names += [f"c_{len(names) + i}_{lib}" for i in range(n)]
    counts = sparse.csr_matrix(rng.poisson(3.0, (len(labels), 12)).astype(np.int32) + 1)
    a = ad.AnnData(X=log1p_cp10k(counts), layers={"counts": counts})
    a.obs_names, a.var_names = names, [f"g{i}" for i in range(12)]
    a.obs["celltype"] = labels
    return a


SIZES = {"1": {"V-CM": 40, "Blood": 30, "Hepatocyte": 50, "aSHF": 20},
         "2": {"V-CM": 35, "Endocardium-1": 10, "Hepatocyte": 45},
         "3": {"V-CM": 61, "OFT-CM": 13, "Endocardium-1": 70, "Early VIC": 29, "Epicardium": 23, "BEC": 37,
               "Blood": 55, "Macrophage": 9, "Hepatocyte": 11}}
LATE = {"E15.5": {"r_cm": 120, "r_endocardial_ec": 9, "r_vascular_ec": 7, "r_ec_ambiguous": 0, "r_epicardial": 3,
                  "r_blood": 4, "r_neural_crest": 0, "r_hepatic": 0, "r_mesenchyme": 17, "r_epithelial": 0,
                  "r_mural": 2, "r_immune": 7, "r_lymphatic": 1, "r_unassigned": 9},
        "E16.5": {"r_cm": 93, "r_endocardial_ec": 9, "r_vascular_ec": 6, "r_ec_ambiguous": 1, "r_epicardial": 3,
                  "r_blood": 3, "r_neural_crest": 0, "r_hepatic": 0, "r_mesenchyme": 18, "r_epithelial": 0,
                  "r_mural": 3, "r_immune": 5, "r_lymphatic": 0, "r_unassigned": 6},
        "E17.5": {"r_cm": 91, "r_endocardial_ec": 12, "r_vascular_ec": 11, "r_ec_ambiguous": 1, "r_epicardial": 3,
                  "r_blood": 3, "r_neural_crest": 0, "r_hepatic": 0, "r_mesenchyme": 25, "r_epithelial": 1,
                  "r_mural": 4, "r_immune": 11, "r_lymphatic": 0, "r_unassigned": 8}}
EARLY = {"E9.5": {"r_cm": 70, "r_endocardial_ec": 5, "r_vascular_ec": 0, "r_ec_ambiguous": 1, "r_epicardial": 1,
                  "r_blood": 9, "r_neural_crest": 0, "r_hepatic": 0, "r_mesenchyme": 1, "r_epithelial": 0,
                  "r_mural": 0, "r_immune": 0, "r_lymphatic": 0, "r_unassigned": 2},
         "E14.5": {"r_cm": 177, "r_endocardial_ec": 8, "r_vascular_ec": 5, "r_ec_ambiguous": 0, "r_epicardial": 3,
                   "r_blood": 2, "r_neural_crest": 0, "r_hepatic": 0, "r_mesenchyme": 26, "r_epithelial": 0,
                   "r_mural": 1, "r_immune": 4, "r_lymphatic": 0, "r_unassigned": 9}}


def _logit(p):
    return math.log(p / (1 - p))


def _expit(x):
    return 1 / (1 + math.exp(-x))


# ------------------------------------------------------------------ step 0: the map, the rule, the recomputation


def test_the_map_places_every_label_of_the_committed_manifests_and_they_are_the_frozen_records():
    for path, stages, gate in ((lm.LATE_MANIFEST, tuple(lm.ANCHOR_DAYS), "N57"),
                               (lm.EARLY_MANIFEST, tuple(lm.EARLY_DAYS), "N54")):
        assert lm.manifest_reasons(json.loads(path.read_text()), stages, gate) == []
    placed = [*lm.HEART, *lm.BLOOD, *lm.OUT_LABELS, *lm.DROPPED]
    assert len(placed) == len(set(placed)) == 14
    assert set(lm.DESIGN_HEART) - set(lm.HEART) == {"r_mural", "r_mesenchyme"}  # the rule-3 fix, and only it
    assert {"r_mesenchyme", "r_mural"} <= set(lm.OUT_LABELS) and lm.CM in lm.HEART


def test_a_manifest_that_is_not_the_frozen_record_is_named():
    good = {"config": {"sha256": lm.RELABEL_YAML_SHA256}, "stages": list(lm.ANCHOR_DAYS), "gate": {"row": "N57"},
            "authors_labels_read": False, "usable": True, "label_counts": LATE}
    assert lm.manifest_reasons(good, tuple(lm.ANCHOR_DAYS), "N57") == []
    for change, word in (({"config": {"sha256": "x"}}, "settings"), ({"stages": ["E15.5"]}, "stages"),
                         ({"gate": {"row": "N58"}}, "gate row"), ({"authors_labels_read": True}, "authors"),
                         ({"usable": False}, "usable"),
                         ({"label_counts": {**LATE, "E16.5": {**LATE["E16.5"], "r_new": 3}}}, "does not place"),
                         ({"label_counts": {**LATE, "E17.5": {**LATE["E17.5"], "r_cm": 1.5}}}, "whole number")):
        assert any(word in r for r in lm.manifest_reasons({**good, **change}, tuple(lm.ANCHOR_DAYS), "N57")), word


def test_the_rule_by_hand():
    own = {"cm": 0.2, "endocardium": 0.3, "epicardium": 0.1, "vascular_endothelium": 0.2, "blood": 0.2}
    res = lm.compute(own, LATE, EARLY)
    heart = {st: sum(LATE[st][x] for x in lm.HEART) for st in LATE}
    blood = {st: LATE[st]["r_blood"] + LATE[st]["r_immune"] for st in LATE}
    cm_a = sum(LATE[st]["r_cm"] for st in LATE) / sum(heart.values())
    blood_a = sum(blood.values()) / (sum(heart.values()) + sum(blood.values()))
    day = sum({"E15.5": 15.5, "E16.5": 16.5, "E17.5": 17.5}[st] * (heart[st] + blood[st]) for st in LATE) / (
        sum(heart.values()) + sum(blood.values()))
    w = 2.0 / (day - 10.5)
    cm = _expit(_logit(0.25) + w * (_logit(cm_a) - _logit(0.25)))
    b = _expit(_logit(0.2) + w * (_logit(blood_a) - _logit(0.2)))
    assert res["own_cm_share_of_heart"] == pytest.approx(0.25, abs=1e-15)
    assert res["anchor"]["cm_share_of_heart"] == pytest.approx(cm_a, abs=1e-15)
    assert res["anchor"]["mean_day"] == pytest.approx(day, abs=1e-12) and res["weight"] == pytest.approx(w, abs=1e-12)
    assert res["cm_share_of_heart"] == pytest.approx(cm, abs=1e-12)
    assert res["blood_share_of_chamber"] == pytest.approx(b, abs=1e-12)
    s = res["shares"]
    heart_cats = list(t1_scope.load_spec().heart)
    assert t1_scope.share_of(s, "cm", heart_cats) == pytest.approx(cm, abs=1e-12) and s["blood"] == pytest.approx(b)
    assert s["endocardium"] / s["epicardium"] == pytest.approx(3.0) and sum(s.values()) == pytest.approx(1.0)
    mid = _expit((_logit(0.25) + _logit(cm)) / 2)
    assert res["t72_bar_c"]["threshold"] == pytest.approx(mid + 0.03, abs=1e-12)
    # the design's map (r_mesenchyme and r_mural in heart) lowers the anchor's cardiomyocyte share; reported only
    assert res["design_map"]["anchor"]["cm_share_of_heart"] < cm_a
    t = res["change_transfer"]["rule3_map"]
    h9 = sum(EARLY["E9.5"][x] for x in lm.HEART)
    h14 = sum(EARLY["E14.5"][x] for x in lm.HEART)
    want = _expit(_logit(0.25) + 0.4 * (_logit(EARLY["E14.5"]["r_cm"] / h14) - _logit(EARLY["E9.5"]["r_cm"] / h9)))
    assert t["cm_share_of_heart"] == pytest.approx(want, abs=1e-12)


def test_the_steps_own_shares_and_the_independent_path_agree():
    ref = _stage(SIZES)
    own = lm.own_shares(ref)
    assert own == pytest.approx(lm.own_shares_direct(ref), abs=1e-15)
    lib3 = SIZES["3"]
    chamber = {lab: n for lab, n in lib3.items() if LABELS[lab] != "surrounding"}
    assert own["cm"] == pytest.approx((61 + 13) / sum(chamber.values()))
    main = lm.compute(own, LATE, EARLY)
    again = lm.recompute(lm.own_shares_direct(ref), LATE)
    assert lm.largest_difference(main, again) < 1e-13
    off = {**again, "cm_share_of_heart": again["cm_share_of_heart"] + 2e-9}
    assert lm.largest_difference(main, off) > lm.TOL
    other = {**again, "shares": {**again["shares"], "pericardium": 0.0, "x": 0.0}}
    assert lm.largest_difference(main, other) == math.inf


def test_library_3_is_taken_by_its_suffix_alone():
    ref = _stage({"1": {"V-CM": 5}, "3": {"V-CM": 7, "Blood": 3}})
    assert lm.library3_pool(ref).n_obs == 10
    with pytest.raises(SystemExit, match="do not separate library 3"):
        lm.library3_pool(_stage({"3": {"V-CM": 7}}))


# ------------------------------------------------------------------ step 0: refusals and the outputs


def _ledger(tmp_path, first: str | None) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text("# ledger\n" + ("" if first is None else f"| {first} | T1 | step 0 |\n"))
    return p


def _register(tmp_path, amended: bool = True) -> Path:
    p = tmp_path / "EXTERNAL_DATA.md"
    tail = " Amendment of 2026-11-25 (rows N64, T72 and R54): the counts enter." if amended else ""
    p.write_text(f"| X001 | x |\n| X002 | Feng | data | stages | removed | T1 | No absolute external levels.{tail} |\n")
    return p


FIRST = "N64 (step 0, target-free: `experiments/t1_scope/late_mix.py` with `experiments/t1_scope/n64_step0.sh`)"


@pytest.fixture
def step0(tmp_path, monkeypatch):
    """Row N64's step with a temporary ledger, register, outputs and manifests, HEAD the commit C, a clean tree."""
    out = tmp_path / "out"
    late = tmp_path / "late.json"
    early = tmp_path / "early.json"
    base = {"config": {"sha256": lm.RELABEL_YAML_SHA256}, "authors_labels_read": False, "usable": True}
    late.write_text(json.dumps({**base, "stages": list(lm.ANCHOR_DAYS), "gate": {"row": "N57"}, "label_counts": LATE}))
    early.write_text(json.dumps({**base, "stages": list(lm.EARLY_DAYS), "gate": {"row": "N54"},
                                 "label_counts": EARLY}))
    monkeypatch.setattr(lm, "TRIES", _ledger(tmp_path, FIRST))
    monkeypatch.setattr(lm, "REGISTER", _register(tmp_path))
    monkeypatch.setattr(lm, "OUT", out)
    monkeypatch.setattr(lm, "LATE_MANIFEST", late)
    monkeypatch.setattr(lm, "EARLY_MANIFEST", early)
    git = {"rev-parse": C, "status": ""}
    monkeypatch.setattr(harness, "_git", lambda *a: git[a[0]])
    data = tmp_path / "data"
    (data / "T1").mkdir(parents=True)
    _stage(SIZES).write_h5ad(data / "T1" / "E10.5_RNA.h5ad")
    return {"out": out, "git": git, "data": data, "tmp": tmp_path}


def test_step0_writes_shares_the_builder_takes(step0, capsys):
    rec = lm.run(C, step0["data"])
    out = step0["out"]
    raw = (out / "shares.json").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == rec["shares_sha256"]
    assert capsys.readouterr().out.splitlines()[0] == f"shares.json sha256 {rec['shares_sha256']}"
    shares = json.loads(raw)
    assert shares["row"] == "N64" and shares["late_mix"] == "e15_17" and shares["commit"] == C
    assert shares["own_shares"] == lm.own_shares(_stage(SIZES))
    sys.path.insert(0, str(ROOT / "experiments" / "builds"))
    import build_test as bt

    values = bt._late_mix_values(raw, "e15_17")  # the builder's own reading of the file
    moved, _ = bt._t1_late_mix_shares(shares["own_shares"], {**values, "late_mix": "e15_17", "shares_file": "x",
                                                              "shares_sha256": rec["shares_sha256"]})
    assert moved == pytest.approx(shares["shares"], abs=0)
    assert rec["recomputation"]["largest_difference"] < 1e-13
    text = (out / "n64.txt").read_text()
    assert f"shares.json sha256 {rec['shares_sha256']}" in text and "row T72's bar (c)" in text
    assert json.loads((out / "n64.json").read_text())["t72_bar_c"]["threshold"] == rec["t72_bar_c"]["threshold"]


@pytest.mark.parametrize("first, match", [
    (None, "no single row N64"),
    ("N64 (step 0: `experiments/t1_scope/late_mix.pyc`)", "no single row N64"),
    ("N640 (step 0: `experiments/t1_scope/late_mix.py`)", "no single row N64"),
    ("N64 step 0: `experiments/t1_scope/late_mix.py`", "no single row N64"),
])
def test_step0_without_its_row_is_refused(step0, monkeypatch, first, match):
    monkeypatch.setattr(lm, "TRIES", _ledger(step0["tmp"], first))
    with pytest.raises(SystemExit, match=match):
        lm.run(C, step0["data"])
    assert not step0["out"].exists()


def test_step0_refusals_before_anything_loads(step0, monkeypatch):
    monkeypatch.setattr(lm, "load_stage", lambda *a, **k: pytest.fail("nothing loads before a refusal"))
    with pytest.raises(SystemExit, match="HEAD is not the commit given"):
        lm.run("d" * 40, step0["data"])
    with pytest.raises(SystemExit, match="HEAD is not the commit given"):
        lm.run(None, step0["data"])
    step0["git"]["status"] = " M x"
    with pytest.raises(SystemExit, match="not clean"):
        lm.run(C, step0["data"])
    step0["git"]["status"] = ""
    monkeypatch.setattr(lm, "REGISTER", _register(step0["tmp"], amended=False))
    with pytest.raises(SystemExit, match="does not hold \"Amendment of 2026-11-25"):
        lm.run(C, step0["data"])
    monkeypatch.setattr(lm, "REGISTER", _register(step0["tmp"]))
    bad = json.loads(lm.LATE_MANIFEST.read_text())
    bad["authors_labels_read"] = True
    lm.LATE_MANIFEST.write_text(json.dumps(bad))
    with pytest.raises(SystemExit, match="authors_labels_read false"):
        lm.run(C, step0["data"])
    step0["out"].mkdir()
    (step0["out"] / "n64.txt").write_text("x")
    with pytest.raises(SystemExit, match="never overwritten"):
        lm.run(C, step0["data"])


def test_step0_writes_nothing_when_the_recomputation_differs(step0, monkeypatch):
    inner = lm.recompute
    monkeypatch.setattr(lm, "recompute", lambda *a: {**inner(*a), "weight": inner(*a)["weight"] + 1e-6})
    with pytest.raises(SystemExit, match="the independent recomputation differs"):
        lm.run(C, step0["data"])
    assert not step0["out"].exists()


# ------------------------------------------------------------------ T72: truths, arms and the break-even


def _geometry(seed: int = 5) -> dict:
    """kind_gap's geometry on synthetic categories: library 3's halves A (truth) and B (the arms' pool)."""
    rng = np.random.default_rng(seed)
    chamber = ["cm", "endocardium", "epicardium", "vascular_endothelium", "blood"]
    cats = np.array(rng.choice([*chamber, "surrounding"], 900, p=[0.2, 0.3, 0.08, 0.12, 0.2, 0.1]), dtype=object)
    cat = cats.astype(str)
    half_a, half_b = np.arange(0, 450), np.arange(450, 900)
    pools = {"b": {c: half_b[cat[half_b] == c] for c in chamber}, "o": {c: np.array([], int) for c in chamber}}
    ch = np.isin(cat[half_b], chamber)
    mix = {c: float((cat[half_b][ch] == c).sum()) / ch.sum() for c in chamber}
    return {"cat": cat, "half_a": half_a, "half_b": half_b, "chamber": chamber, "pools": pools, "lib3_mix": mix}


def test_where_blood_is_own_the_truths_are_t30s():
    d = _geometry()
    assert np.array_equal(t72.truth_rows(d, None, None), d["half_a"])
    for cm in (0.30, 0.45, 0.60):
        assert np.array_equal(t72.truth_rows(d, cm, None), t30.cm_truth_rows(d, cm))


@pytest.mark.parametrize("cm, blood", [(None, 0.12), (0.30, 0.12), (0.44, 0.12), (0.60, 0.12)])
def test_where_blood_moves_the_truth_holds_both_shares(cm, blood):
    d = _geometry()
    rows = t72.truth_rows(d, cm, blood)
    assert set(rows.tolist()) <= set(d["half_a"].tolist()) and len(set(rows.tolist())) == len(rows)
    got_cm, got_blood = t72.shares_of(d, rows)
    own_cm, _ = t72.shares_of(d, d["half_a"])
    assert got_cm == pytest.approx(own_cm if cm is None else cm, abs=0.02)
    assert got_blood == pytest.approx(blood, abs=0.02)
    assert np.array_equal(rows, t72.truth_rows(d, cm, blood))  # its own stream: the same rows twice


def test_the_arm_at_half_bs_own_mix_is_lib3_bases_draw_and_never_a_truth_cell(monkeypatch):
    d = _geometry()
    heart = list(t1_scope.load_spec().heart)
    own = {"cm_share_of_heart": t1_scope.share_of(d["lib3_mix"], "cm", heart),
           "blood_share_of_chamber": d["lib3_mix"]["blood"]}
    d["lm16_mix"] = t72.lm16_mix(d["lib3_mix"], own)
    monkeypatch.setattr(kg, "_D", d)
    for n in (500, 777, 1000):
        assert {c: q for c, q in quotas(d["lm16_mix"], n).items() if q} == {
            c: q for c, q in quotas(d["lib3_mix"], n).items() if q}
        for seed in (1000001, 1000007):
            a, b = t72.arm_rows(t72.ARM, n, seed), t72.arm_rows(t72.BASE, n, seed)
            assert np.array_equal(a, b) and np.array_equal(b, kg.arm_rows("lib3_base", n, seed))
    d["lm16_mix"] = t72.lm16_mix(d["lib3_mix"], {"cm_share_of_heart": 0.45, "blood_share_of_chamber": 0.12})
    rows = t72.arm_rows(t72.ARM, 800, 1000003)
    assert not np.isin(rows, d["half_a"]).any() and not np.array_equal(rows, t72.arm_rows(t72.BASE, 800, 1000003))
    cats, n = np.unique(d["cat"][rows], return_counts=True)
    assert dict(zip(cats.tolist(), n.tolist(), strict=True)) == {
        c: q for c, q in quotas(d["lm16_mix"], 800).items() if q}


def test_the_smoke_test_builds_the_twelve_scenarios_without_scoring(monkeypatch):
    from vec.boards import Board

    d0 = _geometry()
    d0.update(X=sparse.csr_matrix(np.ones((900, 3), dtype=np.float32)), tot=np.arange(900, dtype=float) + 100.0,
              ref_rows=np.arange(450, 900))
    monkeypatch.setattr(t72, "load_boards", lambda: {kg.TWIN: Board(key=kg.TWIN, task="T1", genes=["g0"],
                                                                     needs_coords=False)})
    monkeypatch.setattr(t72, "score_board", lambda *a, **k: pytest.fail("the smoke test scores nothing"))
    shares = {"cm_share_of_heart": 0.45, "blood_share_of_chamber": 0.12}
    res = t72.smoke(d0, shares)
    assert [(r["cm"], r["blood"]) for r in res] == [(c, b) for c in t72.CMS for b in t72.BLOODS]
    for r in res:
        want_cm = {"own": t72.shares_of(d0, d0["half_a"])[0], "lm16": 0.45}.get(r["cm"], None)
        want_cm = float(r["cm"]) if want_cm is None else want_cm
        assert r["truth_cm_of_heart"] == pytest.approx(want_cm, abs=0.02)
        if r["blood"] == "lm16":
            assert r["truth_blood_of_chamber"] == pytest.approx(0.12, abs=0.02)
        assert r["lm16_cm_of_heart"] == pytest.approx(0.45, abs=0.01)
        assert r["lm16_blood_of_chamber"] == pytest.approx(0.12, abs=0.01)
        assert not r["draws_truth_cell"] and r["median_library_over_base"] > 0
    assert set(d0) == set(_geometry()) | {"X", "tot", "ref_rows"}  # the geometry is left as it was


@pytest.mark.parametrize("points, want", [
    ([(0.24, -5.0), (0.30, -2.0), (0.33, 0.5), (0.44, 3.0), (0.45, 3.5), (0.60, 6.0)], 0.30 + 2.0 / 2.5 * 0.03),
    ([(0.24, 1.0), (0.30, 2.0), (0.60, 6.0)], 0.24),
    ([(0.24, -1.0), (0.30, -2.0), (0.60, -0.1)], math.inf),
    ([(0.60, 6.0), (0.24, -4.0), (0.30, 0.0), (0.33, 1.0)], 0.30),
    ([(0.24, -1.0), (0.30, 1.0), (0.33, -1.0), (0.45, 1.0), (0.60, 2.0)], 0.39),
])
def test_the_break_even_share(points, want):
    assert t72.break_even(points) == pytest.approx(want)


# ------------------------------------------------------------------ T72: the reading


SHA = "e" * 64
SHARES = {"own_cm_share_of_heart": 0.243, "cm_share_of_heart": 0.44, "blood_share_of_chamber": 0.14}
TRUTH_CM = {"own": 0.237, "0.30": 0.30, "0.33": 0.33, "lm16": 0.44, "0.45": 0.45, "0.60": 0.60}


def _scenario(cm: str, blood: str, d: float | list[float], base: float = 80.0, **kw) -> dict:
    ds = [d] * 10 if not isinstance(d, list) else d
    rows = [{"seed": s, "scores": {t72.BASE: base + 0.01 * s, t72.ARM: base + 0.01 * s + ds[s]},
             "skills": {}, "median_library": {t72.BASE: 100.0, t72.ARM: 104.0},
             "repeat_share": {t72.BASE: 0.0, t72.ARM: 0.1}, "rows_sha256": {t72.BASE: f"b{s}", t72.ARM: f"a{s}"},
             "draws_truth_cell": {t72.BASE: False, t72.ARM: False}} for s in range(10)]
    return {"try": "T72", "commit": C, "tree_clean": True, "cm": cm, "blood": blood, "truth_cells": 3000,
            "truth_cm_of_heart": TRUTH_CM[cm], "truth_blood_of_chamber": 0.18 if blood == "own" else 0.14,
            "shares_sha256": SHA, "recomputation": {"largest_difference": 0.0, "tolerance": 1e-9}, "rows": rows, **kw}


OWN_BLOOD_D = {"own": -5.0, "0.30": -2.0, "0.33": 0.5, "lm16": 3.0, "0.45": 3.5, "0.60": 6.0}


def _loaded() -> dict:
    out = {}
    for cm in t72.CMS:
        for blood in t72.BLOODS:
            out[(cm, blood)] = _scenario(cm, blood, OWN_BLOOD_D[cm] if blood == "own" else OWN_BLOOD_D[cm] - 0.5)
    return out


def _stored(loaded: dict, shift: float = 0.0) -> dict:
    return {"rows": [{"seed": x["seed"], "scores": {t72.BASE: x["scores"][t72.BASE] + shift}}
                     for x in loaded[("own", "own")]["rows"]]}


def test_the_reading_passes_when_all_three_bars_hold():
    loaded = _loaded()
    r = t72.t72_reading(loaded, _stored(loaded, 0.019), SHARES, SHA, C)
    assert r["outcome"] == "PASS" and r["bars"] == {"a": True, "b": True, "c": True} and not r["failed"]
    assert r["break_even"]["share"] == pytest.approx(0.324)
    mid = _expit((_logit(0.243) + _logit(0.44)) / 2)
    assert r["break_even"]["threshold"] == pytest.approx(mid + 0.03)
    assert t72.outcome_lines(r)[0] == "Outcome: T72 PASS"
    assert r["by_scenario"]["cm lm16, blood lm16"]["median_library_over_base"] == pytest.approx(1.04)


@pytest.mark.parametrize("over, failed", [
    ({("lm16", "lm16"): _scenario("lm16", "lm16", 1.9)}, "(a)"),
    ({("lm16", "lm16"): _scenario("lm16", "lm16", [3.0] * 7 + [-1.0] * 3)}, "(a)"),
    ({("lm16", "own"): _scenario("lm16", "own", [4.0] * 7 + [-0.1] * 3)}, "(b)"),
    ({("0.33", "own"): _scenario("0.33", "own", -2.0), ("lm16", "own"): _scenario("lm16", "own", 2.1),
      ("0.30", "own"): _scenario("0.30", "own", -3.0)}, "(c)"),
])
def test_the_reading_stops_when_a_bar_fails(over, failed):
    loaded = _loaded()
    loaded.update(over)
    r = t72.t72_reading(loaded, _stored(loaded), SHARES, SHA, C)
    assert r["outcome"] == "STOP" and any(f.startswith(failed) for f in r["failed"]), r["failed"]
    assert t72.outcome_lines(r)[0].startswith("Outcome: T72 STOP (")


def test_an_arm_that_is_lib3_bases_draw_everywhere_reads_stop():
    loaded = _loaded()
    for r in loaded.values():
        for x in r["rows"]:
            x["rows_sha256"][t72.ARM] = x["rows_sha256"][t72.BASE]
    r = t72.t72_reading(loaded, _stored(loaded), SHARES, SHA, C)
    assert r["outcome"] == "STOP" and r["idle"] and any("lib3_base's draw" in f for f in r["failed"])


def _void(loaded, stored=None, shares=SHARES, sha=SHA, commit=C):
    r = t72.t72_reading(loaded, _stored(loaded) if stored is None else stored, shares, sha, commit)
    assert r["outcome"] == "VOID", r
    assert t72.outcome_lines(r)[0].startswith("Outcome: T72 VOID (")
    return " ".join(r["reasons"])


def test_each_void_condition():
    loaded = _loaded()
    del loaded[("0.45", "lm16")]
    assert "cm 0.45, blood lm16 is missing" in _void(loaded)
    loaded = _loaded()
    loaded[("0.30", "own")]["rows"].pop()
    assert "a seed or an arm is missing" in _void(loaded)
    loaded = _loaded()
    loaded[("own", "lm16")]["commit"] = "d" * 40
    assert "more than one commit" in _void(loaded)
    assert "did not run at the row's commit" in _void(_loaded(), commit="d" * 40)
    loaded = _loaded()
    loaded[("0.60", "own")]["tree_clean"] = False
    assert "dirty tree" in _void(loaded)
    loaded = _loaded()
    assert "does not reproduce T30's stored per-seed scores" in _void(loaded, _stored(loaded, 0.021))
    assert "T30's stored scores" in _void(_loaded(), stored={})  # (an empty file: no rows to compare)
    loaded = _loaded()
    loaded[("0.33", "lm16")]["recomputation"]["largest_difference"] = 2e-9
    assert "differ from the independent recomputation" in _void(loaded)
    loaded = _loaded()
    loaded[("lm16", "own")]["rows"][3]["draws_truth_cell"][t72.ARM] = True
    assert "draws a truth-half cell" in _void(loaded)
    loaded = _loaded()
    loaded[("own", "own")]["shares_sha256"] = "f" * 64
    assert "another shares file" in _void(loaded)
    assert "another shares file" in _void(_loaded(), sha=None)
    assert "shares file was not read" in _void(_loaded(), shares=None)


def test_the_two_rows_must_name_one_shares_file(tmp_path):
    p = tmp_path / "TRIES.md"
    n64 = f"| N64 (step 0) | T1 | written: `shares.json sha256 {SHA}` |"
    first = f"| T72 (the read; `experiments/t1_scope/t72_late_mix.py`; `shares.json sha256 {SHA}`) | T1 |"
    for text, want in ((f"{n64}\n{first}", SHA), (n64, None), (first, None),
                       (f"{n64}\n{first.replace(SHA, 'f' * 64)}", None),
                       (f"{n64}\n| T72 (the read) | T1 | `shares.json sha256 {SHA}` |", None),
                       (f"{n64}\n{first}\n{first}", None)):
        p.write_text(f"# ledger\n{text}\n")
        assert t72.shares_sha_named(p) == want


def test_the_shares_file_is_taken_only_as_committed_with_the_named_sha(tmp_path):
    rec = {"row": "N64", "late_mix": "e15_17", "cm_share_of_heart": 0.4}
    raw = json.dumps(rec).encode()
    path = tmp_path / "shares.json"
    path.write_bytes(raw)
    sha = hashlib.sha256(raw).hexdigest()
    assert t72.load_shares_file(sha, path, lambda p: raw) == rec
    with pytest.raises(SystemExit, match="not committed at HEAD as it is in the tree"):
        t72.load_shares_file(sha, path, lambda p: None)
    with pytest.raises(SystemExit, match="not committed at HEAD as it is in the tree"):
        t72.load_shares_file(sha, path, lambda p: raw + b" ")
    with pytest.raises(SystemExit, match="is not the one rows N64 and T72 name"):
        t72.load_shares_file("0" * 64, path, lambda p: raw)
    other = json.dumps({**rec, "late_mix": "e14"}).encode()
    path.write_bytes(other)
    with pytest.raises(SystemExit, match="is not row N64's output"):
        t72.load_shares_file(hashlib.sha256(other).hexdigest(), path, lambda p: other)


def test_a_run_is_refused_before_anything_loads(tmp_path, monkeypatch):
    p = tmp_path / "TRIES.md"
    first = f"| T72 (the read; `experiments/t1_scope/t72_late_mix.py`; `shares.json sha256 {SHA}`) | T1 |"
    p.write_text(f"# ledger\n| N64 (step 0) | T1 | `shares.json sha256 {SHA}` |\n{first}\n")
    monkeypatch.setattr(t72, "TRIES", p)
    git = {"rev-parse": C, "status": ""}
    monkeypatch.setattr(harness, "_git", lambda *a: git[a[0]])
    path = tmp_path / "t72_cmown_bloodown.json"
    assert t72.refuse_run(C, path) == SHA
    for commit in ("d" * 40, None):
        with pytest.raises(SystemExit, match="HEAD is not the commit given"):
            t72.refuse_run(commit, path)
    git["status"] = " M x"
    with pytest.raises(SystemExit, match="not clean"):
        t72.refuse_run(C, path)
    git["status"] = ""
    path.write_text("{}")
    with pytest.raises(SystemExit, match="never run twice"):
        t72.refuse_run(C, path)
    path.unlink()
    p.write_text(f"# ledger\n{first}\n")
    with pytest.raises(SystemExit, match="do not each name the one same"):
        t72.refuse_run(C, path)
    p.write_text("# ledger\n| T72 (the read) | T1 |\n")
    with pytest.raises(SystemExit, match="no single row T72"):
        t72.refuse_run(C, path)


def test_the_scripts_print_no_score():
    for name in ("t72_late_mix.py", "late_mix.py"):
        src = (ROOT / "experiments" / "t1_scope" / name).read_text()
        body = src.split('"""', 2)[2]
        prints = [ln.strip() for ln in body.splitlines() if "print(" in ln]
        for ln in prints:
            assert "score" not in ln.lower() or "no score" in ln, ln
    assert "upload" not in (ROOT / "experiments" / "t1_scope" / "t72_late_mix.py").read_text().split('"""', 2)[2]
