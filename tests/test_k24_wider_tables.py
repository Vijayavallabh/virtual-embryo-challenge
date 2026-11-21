"""Tests for experiments/t3_decomp/k24_wider_tables.py (try K24: the barred factor times two wider outside-knockout
tables and their equal-weight joint, each brought to the named form's table size, with permuted controls and size
twins, after a guard part) and its launcher run_k24.sh. Offline, synthetic and fast, in the manner of
tests/test_k23_outside_exponents.py: K18's plain fold fixture with stubbed construct, composition and assemble,
hand-made tables in temporary folders outside the repository, synthetic task files, a temporary git repository for the
launcher's refusals. No competition data and no outside data. K22's and K23's committed result files are read for their
pins and figures only."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3 = REPO / "experiments" / "t3_decomp"
T3X = REPO / "experiments" / "t3_ext"
for p in (T3, T3X):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402
import k22_outside_forms as k22  # noqa: E402
import k23_outside_exponents as k23  # noqa: E402
import k24_wider_tables as k24  # noqa: E402
import outside_forms as ext  # noqa: E402

K23_DIR = T3 / "results" / "k23_outside_exponents"
TWO_SHA = "9490f524e40e146480c775f85ae8829ec9c005607104ad8ace129c5c46ed6228"
G, N_POOL, N_DRAW = 14, 90, 150
GENES = [f"g{i}" for i in range(G)]
FILE_GENES = GENES[::-1] + ["extra"]
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES")


@pytest.fixture(autouse=True)
def clean_state():
    yield
    k24._K24.clear()
    k22._EXT.clear()
    k23._K23.clear()


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ (1) constants, pins


BYTE_IDENTICAL = {
    "rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13",
    "t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
    "build_submission.py": "1194a63993a27dde4dedfdb9180de035c5f8f2a57aaa7d97b38e7a15ea9e0c74",
    "outside_forms.py": "4837687d47f134ff690177ecaaec7c712952c238d0e7c806c5c8ef2cecb3ceb5",
    "k22_outside_forms.py": "6482585242ad6213eb550ceac54f3b9f911da6f8b0c2901216ea74ff24c7dc19",
    "k23_outside_exponents.py": "98bdac5b9f2413bc6d11e470fdbd2202df5ab3dc22548c8940812ac1f455eb98",
    "build_outside_form.py": "6eb64142ab58e743746b22371a3329347a68d159c37ffb0fdf4607d5fd91516f",
    "build_outside_form.sh": "d54429326d2e5cb1c7675840d4b2793859060dcffe5f8c528f864d7bda057a77",
    "build_two_knockout.py": "e6bff930c32d3377bd1b5dc8174c2d60814fa0867d160ab64f1d577f94a50c9b",
    "build_two_knockout.sh": "dc1f8cff8fc2eb9cfa8c6bf9b85f3414e5639f8dbf1ca4984fd56ed6ffd69add",
    "build_combination.py": "ef94fb068845ea3b23db69dc1d6b1d936fa8219f51b7cec07916b5f891ada915",
    "build_combination.sh": "496f371db726407106a776dc5dfac26f2c8b41badc02576147b317063d97c45e",
    "run_k22.sh": "fdadb4b4c07dc30c4823a97d6662d8f2fb17d0f537c96e778f53986fdae44001",
    "run_k23.sh": "14e3705edda1b4434bf76fab81d7d79fe28a5c05db597d05093f81a58358da48",
}


def test_the_files_k24_leaves_alone_are_byte_identical():
    got = {f: _sha(T3 / f) for f in BYTE_IDENTICAL}
    assert got == BYTE_IDENTICAL
    assert {k: BYTE_IDENTICAL[k] for k in ext.FROZEN_SHA256} == ext.FROZEN_SHA256


def test_arms_constants_and_the_shared_modules():
    assert k24.ARMS == ("barred_form", "x_e2", "x_sc3", "x_sc3_perm", "x_sc3_size", "x_bulk", "x_bulk_perm",
                        "x_bulk_size", "x_joint", "x_joint_perm", "x_joint_size")
    assert k24.GUARDS == ("barred_form", "x_e2") and k24.NAMEABLE == ("x_sc3", "x_bulk", "x_joint")
    assert k24.KEYS == ("bulk", "sc3", "two") and k24.MAX_TRIES == 9
    assert k24.SOURCES == {"two": ("VEC_T3_EXT_TABLES_TWO", "all", "COMBINE.json"),
                           "sc3": ("VEC_T3_EXT_TABLES", "all", "COMBINE.json"),
                           "bulk": ("VEC_T3_EXT_TABLES", "bulk_all", "COMBINE_bulk.json")}
    assert (k24.E_MAX, k24.E_MIN, k24.MIN_GENES, k24.MAX_ABS_LN_SCALED) == (8.0, 0.125, 100, 0.5)
    assert k24.JOINT_CORR == (0.0, 0.9) and k24.SC3_TWO_MAX_CORR == k22.ext.FILE2_MAX_CORR == ext.FILE2_MAX_CORR == 0.9
    assert (k24.BULK_MIN_LINES, k24.BULK_MIN_LIBS_UNIT, k24.BULK_MIN_LIBS_TOTAL) == (12, 2, 36)
    # (b4) as restated on 2026-11-20: AGREE_MULT = 2 gave way to an excess of 0.11 over chance, at n_g lines each
    assert (k24.AGREE, k24.AGREE_FLOOR, k24.AGREE_EXCESS, k24.AGREE_MIN_LINES) == (2 / 3, 50, 0.11, 12)
    assert not hasattr(k24, "AGREE_MULT") and not hasattr(k24, "expected_agreement")
    assert k24.BULK_SIGN_RULE == "single study" and k24.BULK_STUDY == "MGP-DMDD"
    assert (k24.MIN_CHANGED, k24.MAX_ABS_CORR, k24.MIN_STUDIES) == (10, 0.5, 2)
    assert (k24.MAX_AT_CLIP, k24.MAX_DEN_SHIFT, k24.MAX_SD_RATIO) == (k23.MAX_AT_CLIP, k23.MAX_DEN_SHIFT,
                                                                       k23.MAX_SD_RATIO) == (45, 0.05, 1.05)
    assert (k24.MIN_GAIN, k24.MIN_AHEAD, k24.MAX_EXTRA_GATED, k24.GUARD_MEAN, k24.GUARD_SEED) == (0.3, 8, 1, 0.1, 0.5)
    assert ext.PERM_SEED == 20261113 and ext.PERM_BINS == 10
    assert k24.K18_SHA256 == k23.K18_SHA256 == k21.K18_SHA256 and k24.STREAM == k18.STREAM == "t3_s2_form"
    assert k24.TASKS == k18.TASKS and k24.WORLDS == k18.WORLDS and k24.BASE == k18.BARRED
    assert k24.RESULTS == T3 / "results" / "k24_wider_tables" and k24.K23_RESULTS == K23_DIR
    assert k24.LISTING_FILE == T3X / "TABLES_K24.sha256"


def test_s_star_is_k23s_committed_table_figure_and_the_pins_wait_for_the_row():
    pre = json.loads((K23_DIR / "pregate.json").read_text())
    assert k24.S_STAR == pre["pregate"]["exponents"]["2"]["table"]["sd_ln"] == 0.04949423311987974
    # until TABLES_K24.sha256 is committed the pins are None; once it is (row K24's commit), they are its hash and lines
    if not k24.LISTING_FILE.exists():
        assert k24.K24_LISTING_SHA256 is None and k24.K24_TABLE_SHA256["sc3"] is None
        assert k24.K24_TABLE_SHA256["bulk"] is None
    else:
        raw = k24.LISTING_FILE.read_bytes()
        assert k24.K24_LISTING_SHA256 == hashlib.sha256(raw).hexdigest()
        assert k24.K24_TABLE_SHA256 == ext.parse_listing(raw) and ext.tracked(k24.LISTING_FILE)
    assert k24.K24_TABLE_SHA256["two"] == TWO_SHA == ext.K22_TABLE_SHA256["all"] == k23.K23_TABLE_SHA256
    r = json.loads((K23_DIR / "reading.json").read_text())
    assert r["ext"]["tables"] == {"all": TWO_SHA} and r["names"] == ["x_e2"]
    # K23's named form at its committed figures (guard 2's references)
    got = {(w, n): (round(r["worlds"][w][n]["arms"]["x_e2"]["mean"], 2), r["worlds"][w][n]["arms"]["x_e2"][
        "severity_gated"]) for w in k24.WORLDS for n in k24.TASKS}
    assert got == {("plain", "g1"): (60.13, 0), ("plain", "g2"): (62.97, 0), ("dm", "g1"): (55.83, 4),
                   ("dm", "g2"): (64.52, 0)}


def test_the_pins_refuse_until_filled():
    with pytest.raises(k24.Refused, match="not filled"):
        k24.check_pins(None, dict(k24.K24_TABLE_SHA256))
    with pytest.raises(k24.Refused, match="not filled"):
        k24.check_pins("a" * 64, {"two": TWO_SHA, "sc3": "b" * 64, "bulk": None})
    with pytest.raises(k24.Refused, match="not the `all` sha256"):
        k24.check_pins("a" * 64, {"two": "c" * 64, "sc3": "b" * 64, "bulk": "d" * 64})
    k24.check_pins("a" * 64, {"two": TWO_SHA, "sc3": "b" * 64, "bulk": "d" * 64})


# ------------------------------------------------------------------------------------------------ (2) the size rule


def test_the_size_rule_s_exponent_by_hand():
    s = 0.02
    T = np.exp(np.array([s, -s, 0.0, 0.0]))  # ln T: mean 0, population SD s / sqrt(2)
    e = k24.exponent(T)
    assert math.isclose(e, k24.S_STAR * math.sqrt(2) / s, rel_tol=1e-13)
    assert math.isclose(float(np.std(np.log(T ** e))), k24.S_STAR, rel_tol=1e-12)
    assert k24.exponent(np.ones(5)) is None
    # genes exactly 1 count in the SD
    T2 = np.exp(np.array([s, -s]))
    assert math.isclose(k24.exponent(T2), k24.S_STAR / s, rel_tol=1e-13)


def test_c_by_hand_and_the_joint_at_s_star():
    a = np.array([0.1, -0.1, 0.0, 0.0])
    b = np.array([0.0, 0.0, 0.2, -0.2])  # uncorrelated with a in this form
    sc3, bulk = np.exp(a), np.exp(b)
    e1, e2 = k24.exponent(sc3), k24.exponent(bulk)
    c = k24.joint_c(sc3, bulk, e1, e2)
    v = e1 * a + e2 * b
    assert math.isclose(c, k24.S_STAR / float(np.std(v)), rel_tol=1e-13)
    J = k24.k24_table("joint", {"sc3": sc3, "bulk": bulk}, {"e_sc3": e1, "e_bulk": e2, "c": c})
    assert math.isclose(float(np.std(np.log(J))), k24.S_STAR, rel_tol=1e-12)
    # orthogonal tables of equal size: c = 1 / sqrt(2)
    assert math.isclose(c, 1 / math.sqrt(2), rel_tol=1e-12)
    assert np.allclose(J, sc3 ** (c * e1) * bulk ** (c * e2), rtol=0, atol=0)
    with pytest.raises(RuntimeError):
        k24.k24_table("two", {"two": sc3}, {})


def test_each_arm_s_product_by_hand_for_one_gene():
    phi = np.array([2.0, 0.5, 1.5, 0.25, 4.0, 1.0])
    sc3 = np.array([1.1, 1.0, 0.9, 1.05, 0.8, 1.2])
    bulk = np.array([0.9, 1.0, 1.1, 0.95, 1.2, 1.0])
    two = np.array([1.2, 1.0, 0.95, 1.0, 1.1, 0.9])
    basis = np.array([0.1, 0.5, 0.2, 0.9, 0.3, 0.7])
    consts = {"e_sc3": 1.5, "e_bulk": 0.5, "c": 0.6,
              "twins": {"x_sc3": {"a": 1.25, "k": 0.875}, "x_bulk": {"a": 1.1, "k": 0.9},
                        "x_joint": {"a": 1.05, "k": 0.95}}}
    k24._K24["consts"] = consts
    F = {"phi_barred": phi, "ext": {"two": two, "sc3": sc3, "bulk": bulk}}
    F["scaled"] = {k: k24.k24_table(k, F["ext"], consts) for k in ("sc3", "bulk", "joint")}
    F["scaled_perm"] = {k: ext.permuted(v, basis) for k, v in F["scaled"].items()}
    assert math.isclose(k24.arm_factor("x_sc3", F)[0], 2.0 * 1.1 ** 1.5, rel_tol=1e-15)
    assert math.isclose(k24.arm_factor("x_bulk", F)[0], 2.0 * 0.9 ** 0.5, rel_tol=1e-15)
    assert math.isclose(k24.arm_factor("x_joint", F)[0], 2.0 * 1.1 ** 0.9 * 0.9 ** 0.3, rel_tol=1e-14)
    assert math.isclose(k24.arm_factor("x_e2", F)[0], 2.0 * 1.44, rel_tol=1e-15)
    assert math.isclose(k24.arm_factor("x_sc3_size", F)[0], 0.875 * 2.0 ** 1.25, rel_tol=1e-15)
    assert math.isclose(k24.arm_factor("x_bulk_size", F)[0], 0.9 * 2.0 ** 1.1, rel_tol=1e-15)
    assert math.isclose(k24.arm_factor("x_joint_size", F)[0], 0.95 * 2.0 ** 1.05, rel_tol=1e-15)
    assert np.array_equal(k24.arm_factor("barred_form", F), phi)
    # a gene absent from both tables keeps phi_M
    assert k24.arm_factor("x_joint", F)[1] == phi[1]
    # x_e2 is K23's x_e2 through K23's code path, value for value
    k23._K23.update({"2": {"a": 1.0, "k": 1.0}})
    assert np.array_equal(k24.arm_factor("x_e2", F), k23.arm_factor("x_e2", {"phi_barred": phi,
                                                                             "ext": {"all": two}}))
    # the controls: perm of the scaled table; perm(T ** e) = perm(T) ** e; perm(J) = perm(sc3 ** ce) perm(bulk ** ce)
    for kind in ("sc3", "bulk", "joint"):
        assert np.array_equal(k24.arm_factor(f"x_{kind}_perm", F), phi * ext.permuted(F["scaled"][kind], basis))
    assert np.array_equal(ext.permuted(sc3 ** 1.5, basis), ext.permuted(sc3, basis) ** 1.5)
    J = F["scaled"]["joint"]
    assert np.array_equal(ext.permuted(J, basis),
                          ext.permuted(sc3 ** (0.6 * 1.5), basis) * ext.permuted(bulk ** (0.6 * 0.5), basis))
    k24._K24.clear()
    with pytest.raises(RuntimeError, match="not loaded"):
        k24.arm_factor("x_sc3_size", F)


# ------------------------------------------------------------------------------------------------ (3) the pre-gate


def _bases(n, rng):
    return {"plain": {"genes": None, "basis": rng.random(n)}, "dm": {"genes": None, "basis": rng.random(n)}}


def _bulk_rec(n, shares, lines=25, libs=(4, 4), sha="b" * 64, **over):
    units = [f"E-MTAB-7322:L{i}" for i in range(lines)]
    rec = {"modality": "bulk", "studies": ["MGP-DMDD"], "tiers": ["unrelated"], "sign_rule": "single study",
           "genes_by_study": {"MGP-DMDD": [f"L{i}" for i in range(lines)]},
           "units_by_study": {"MGP-DMDD": units},
           "libraries_by_unit": {u: {"ko": [f"{u}k{j}" for j in range(libs[0])],
                                     "control": [f"{u}c{j}" for j in range(libs[1])]} for u in units},
           "k_by_unit": {u: {"E9.5": 0.64} for u in units}, "guards_by_unit": {u: {"E9.5": "passed"} for u in units},
           "line_agreement_by_study": {"MGP-DMDD": list(shares)}, "refstats_sha256": "r" * 64,
           "checked": "by the practice read named in its TRIES row only", "sha256": sha}
    rec.update(over)
    return rec


def _synth(seed=0, n=300, sc3=None, bulk=None, sc3_rec=None, bulk_rec=None, two=None):
    rng = np.random.default_rng(seed)
    phi = np.clip(np.exp(rng.normal(0, 0.6, n)), 0.25, 4.0)
    den = rng.random(n) + 0.05
    z1, z2, z3 = (rng.normal(size=n) for _ in range(3))
    m1, m2, m3 = (rng.random(n) < 0.85 for _ in range(3))
    two = np.where(m1, np.exp(0.025 * z1), 1.0) if two is None else two
    sc3 = np.where(m2, np.exp(0.03 * (0.6 * z1 + 0.8 * z2)), 1.0) if sc3 is None else sc3
    bulk = np.where(m3, np.exp(0.04 * (0.3 * z2 + 0.95 * z3)), 1.0) if bulk is None else bulk
    shares = [None if v == 1.0 else (0.8 if i % 2 == 0 else 0.5) for i, v in enumerate(bulk)]
    two_rec = {"studies": ["GSE204908", "GSE277839"], "tiers": ["unrelated"], "sign_rule": "applied",
               "sha256": "2" * 64, "refstats_sha256": "r" * 64}
    sc3_rec = sc3_rec or {"studies": ["GSE204908", "GSE277839", "GSE299068"], "tiers": ["unrelated"],
                          "sign_rule": "applied", "sha256": "3" * 64, "refstats_sha256": "r" * 64}
    bulk_rec = bulk_rec or _bulk_rec(n, shares)
    loaded = {"listing_sha256": "l" * 64,
              "tables": {"two": (two, two_rec), "sc3": (sc3, sc3_rec), "bulk": (bulk, bulk_rec)},
              "combine": {"two": {"all": "2" * 64}, "sc3": {"all": "3" * 64, "study_GSE204908": "x" * 64},
                          "bulk": {"bulk_all": bulk_rec["sha256"], "bulk_study_MGP-DMDD": bulk_rec["sha256"]}}}
    return SimpleNamespace(phi=phi, den=den, bases=_bases(n, rng), loaded=loaded, two=two, sc3=sc3, bulk=bulk,
                           n=n, shares=shares)


def test_the_pregate_figures_and_no_drop_for_well_formed_tables():
    s = _synth()
    pre = k24.pregate_k24(s.phi, s.loaded, s.bases, s.den)
    assert pre["dropped"] == {}, pre["dropped"]
    assert k24.active_arms(pre) == k24.ARMS and k24.try_count(k24.active_arms(pre)) == 9
    for key, T in (("sc3", s.sc3), ("bulk", s.bulk)):
        assert pre["e"][key] == k24.S_STAR / float(np.std(np.log(T)))
        assert math.isclose(pre["scaled"][key]["sd_ln"], k24.S_STAR, rel_tol=1e-12)
        v = pre["arms"][f"x_{key}"]
        prod = s.phi * T ** pre["e"][key]
        assert v["product"] == ext.figures(s.phi, prod, s.den)
        assert math.isclose(v["a"], np.std(np.log(prod)) / np.std(np.log(s.phi)), rel_tol=1e-15)
        twin = v["k"] * s.phi ** v["a"]
        assert math.isclose(s.den @ twin, s.den @ prod, rel_tol=1e-12)
        assert v["controls"]["plain"] == ext.figures(
            s.phi, s.phi * ext.permuted(T ** pre["e"][key], s.bases["plain"]["basis"]), s.den)
    c = pre["c"]
    assert c == k24.S_STAR / float(np.std(pre["e"]["sc3"] * np.log(s.sc3) + pre["e"]["bulk"] * np.log(s.bulk)))
    assert 0.527 < c <= 0.7072 and math.isclose(pre["scaled"]["joint"]["sd_ln"], k24.S_STAR, rel_tol=1e-12)
    assert 0 <= pre["corr"]["sc3_bulk"] < 0.9 and pre["corr"]["sc3_two"] < 0.9
    assert pre["bulk"]["lines"] == 25 and pre["bulk"]["libraries_ko"] == 100 and pre["bulk"]["libraries_wt"] == 100
    assert json.loads(json.dumps(pre)) == pre  # recorded to full precision
    co = k24.constants_of(pre)
    assert co["e_sc3"] == pre["e"]["sc3"] and co["c"] == c and set(co["twins"]) == set(k24.NAMEABLE)
    lines = "\n".join(k24.pregate_lines(pre))
    for x in ("S* 0.04949423311987974", "x_joint: a", "control (dm basis)", "dropped: none", "somites"):
        assert x in lines, x


def test_the_pregate_agrees_with_k23s_at_e_2_on_one_table():
    s = _synth(1)
    rec = s.loaded["tables"]["two"][1]
    k23p = k23.pregate_exponents(s.phi, {"tables": {"all": (s.two, rec)}, "combine": {}}, s.bases, s.den)
    got, why = k24.arm_figures(s.phi, s.two ** 2.0, s.bases, s.den)
    want = k23p["exponents"]["2"]
    for k in ("product", "controls", "size", "a", "k", "den_shift", "sd_ln_phi_m", "sd_ln_product"):
        assert got[k] == want[k], k
    assert why == []


def test_two_is_checked_against_k23s_committed_figures(monkeypatch):
    s = _synth(2)
    rec = s.loaded["tables"]["two"][1]
    k23p = json.loads(json.dumps(k23.pregate_exponents(s.phi, {"tables": {"all": (s.two, rec)}, "combine": {}},
                                                       s.bases, s.den)))
    monkeypatch.setattr(k24, "S_STAR", k23p["exponents"]["2"]["table"]["sd_ln"])
    pre = k24.pregate_k24(s.phi, s.loaded, s.bases, s.den, None, k23p)
    assert pre["two_checked"] and pre["two"]["a"] == k23p["exponents"]["2"]["a"]
    near = json.loads(json.dumps(k23p))
    near["exponents"]["2"]["k"] *= 1 + 1e-13
    k24.pregate_k24(s.phi, s.loaded, s.bases, s.den, None, near)
    for path, rel in ((("a",), 1e-10), (("controls", "dm", "sd_ln"), 1e-9), (("table", "sd_ln"), 1e-9)):
        bad = json.loads(json.dumps(k23p))
        d = bad["exponents"]["2"]
        for p in path[:-1]:
            d = d[p]
        d[path[-1]] *= 1 + rel
        with pytest.raises(k24.Refused, match="refused"):
            k24.pregate_k24(s.phi, s.loaded, s.bases, s.den, None, bad)
    monkeypatch.setattr(k24, "S_STAR", k23p["exponents"]["2"]["table"]["sd_ln"] * (1 + 1e-9))
    with pytest.raises(k24.Refused, match="S\\*"):
        k24.pregate_k24(s.phi, s.loaded, s.bases, s.den, None, k23p)


def _drops(**kw):
    s = _synth(**kw)
    return k24.pregate_k24(s.phi, s.loaded, s.bases, s.den)


def _arms_left(pre):
    return k24.active_arms(pre)


SC3_ONLY = ("barred_form", "x_e2", "x_sc3", "x_sc3_perm", "x_sc3_size")
BULK_ONLY = ("barred_form", "x_e2", "x_bulk", "x_bulk_perm", "x_bulk_size")


def _bulk_bad(case, n=300):
    shares = [0.8] * n
    if case == "b1_sign":
        return _bulk_rec(n, shares, sign_rule="two_study")
    if case == "b1_modality":
        return _bulk_rec(n, shares, modality="single_cell")
    if case == "b1_study":
        return _bulk_rec(n, shares, studies=["OTHER"])
    if case == "b1_tier":
        return _bulk_rec(n, shares, tiers=["cardiac"])
    if case == "b2":
        return _bulk_rec(n, shares, lines=11, libs=(4, 4))
    if case == "b3_unit":
        r = _bulk_rec(n, shares)
        r["libraries_by_unit"]["E-MTAB-7322:L3"]["ko"] = ["one"]
        return r
    if case == "b3_total":
        return _bulk_rec(n, shares, lines=12, libs=(2, 2))
    if case == "b4":
        return _bulk_rec(n, [0.5] * n)
    if case == "b5":
        r = _bulk_rec(n, shares)
        r["k_by_unit"] = None
        return r
    raise AssertionError(case)


@pytest.mark.parametrize("case", ["b1_sign", "b1_modality", "b1_study", "b1_tier", "b2", "b3_unit", "b3_total", "b4",
                                  "b5"])
def test_a_one_study_bulk_table_that_fails_its_conditions_drops_with_the_joint(case):
    pre = _drops(bulk_rec=_bulk_bad(case))
    assert "bulk" in pre["dropped"] and "joint" in pre["dropped"] and "sc3" not in pre["dropped"]
    tag = case.split("_")[0]
    assert f"({tag})" in pre["dropped"]["bulk"] or (case == "b1_tier" and "cardiac" in pre["dropped"]["bulk"])
    assert _arms_left(pre) == SC3_ONLY and k24.try_count(_arms_left(pre)) == 3


@pytest.mark.parametrize("sign", ["applied", "two_study", "two_thirds", "single study"])
def test_a_bulk_table_of_two_or_more_studies_is_dropped(sign):
    """(b1): row K24 reads the one-study MGP-DMDD bulk table only; a bulk table of two or more studies is dropped with
    the joint, whatever its sign rule, and sc3 is kept."""
    s0 = _synth()
    rec = _bulk_rec(s0.n, s0.shares, studies=["MGP-DMDD", "OTHER"], sign_rule=sign)
    rec["genes_by_study"]["OTHER"] = [f"M{i}" for i in range(25)]
    why, _ = k24.bulk_drops(s0.bulk, rec)
    assert [w for w in why if w.startswith("(b1)")] and "two or more studies is dropped" in why[0]
    pre = _drops(bulk_rec=rec)
    assert "(b1)" in pre["dropped"]["bulk"] and "joint" in pre["dropped"] and "sc3" not in pre["dropped"]
    assert _arms_left(pre) == SC3_ONLY
    # the well-formed one-study table passes (b1)
    why, _ = k24.bulk_drops(s0.bulk, _bulk_rec(s0.n, s0.shares))
    assert not [w for w in why if w.startswith("(b1)")]


def test_b4_by_hand_skips_none_shares_and_counts_lines_not_units():
    n = 300
    s = _synth()
    T = s.bulk
    ne1 = T != 1.0
    share = [None] * n
    idx = [i for i in range(n) if ne1[i]]
    for j, i in enumerate(idx):
        share[i] = 0.8 if j < 60 else (None if j < 120 else 0.4)
    rec = _bulk_rec(n, share, lines=25)
    rec["units_by_study"] = {"MGP-DMDD": [f"u{i}" for i in range(40)]}  # units are not lines
    why, figs = k24.bulk_drops(T, rec)
    G = len(idx) - 60
    assert figs["lines"] == 25 and figs["agree_G"] == G and figs["agree_count"] == 60
    p = sum(math.comb(25, j) for j in range(17, 26)) / 2 ** 25
    # restated (b4), no lines_held_by_study in this manifest: n_g = U = 25 for every gene, recorded as the fallback
    assert math.ceil(2 * 25 / 3) == 17 and math.isclose(figs["E0"], G * 2 * p, rel_tol=1e-14)
    assert figs["n_g_source"].startswith("fallback U") and figs["n_g_histogram"] == {"25": G}
    need = max(50, G * 2 * p + 0.11 * G)
    assert math.isclose(figs["agree_needed"], need, rel_tol=1e-14)
    assert (60 >= need) == (not any("(b4)" in w for w in why))
    assert not any(w.startswith("(b2)") or w.startswith("(b3)") for w in why)
    # a share of exactly 2/3 counts as agreeing
    share2 = [None if v is None else 2 / 3 for v in share]
    _, f2 = k24.bulk_drops(T, _bulk_rec(n, share2))
    assert f2["agree_count"] == f2["agree_G"]


def test_exponents_out_of_bounds_drop_their_table():
    s0 = _synth()
    tiny = np.where(s0.sc3 != 1.0, s0.sc3 ** 0.05, 1.0)  # e above 8
    pre = _drops(sc3=tiny)
    assert "(s1)" in pre["dropped"]["sc3"] and pre["e"]["sc3"] > 8 and _arms_left(pre) == BULK_ONLY
    big = np.where(s0.bulk != 1.0, np.clip(s0.bulk ** 60.0, 0.25, 4.0), 1.0)  # e below 1/8
    pre = _drops(bulk=big)
    assert "(s1)" in pre["dropped"]["bulk"] and pre["e"]["bulk"] < 0.125 and _arms_left(pre) == SC3_ONLY


def test_too_few_genes_and_an_outlier_drop_their_table():
    s0 = _synth()
    few = np.ones(s0.n)
    few[:99] = s0.sc3[:99] * 1.001
    pre = _drops(sc3=few)
    assert "(s2)" in pre["dropped"]["sc3"] and "joint" in pre["dropped"]
    out = np.where(s0.bulk != 1.0, np.exp(0.001 * np.log(s0.bulk)), 1.0)
    out[5] = 2.0  # one gene carries the SD: scaled to S*, its |ln| is far above 0.5
    pre = _drops(bulk=out)
    # (s3) is checked before the scaled table is formed (the fourth review): no scaled figures, no arm figures
    assert "(s3)" in pre["dropped"]["bulk"] and "bulk" not in pre["scaled"] and "x_bulk" not in pre["arms"]


def test_k23s_product_bounds_drop_a_table(monkeypatch):
    monkeypatch.setattr(k24, "MAX_DEN_SHIFT", 0.0)
    pre = _drops()
    assert all("den-weighted total" in pre["dropped"][k] for k in ("sc3", "bulk"))
    assert k24.nothing_runs(_arms_left(pre)) and _arms_left(pre) == ("barred_form", "x_e2")


def test_the_common_drops():
    s0 = _synth()
    rec = dict(s0.loaded["tables"]["sc3"][1])
    for over, why in (({"refstats_sha256": "q" * 64}, "reference statistics"), ({"tiers": ["wnt", "unrelated"]}, "wnt"),
                      ({"tiers": ["cardiac"]}, "cardiac"), ({"studies": ["GSE204908"]}, "1 study"),
                      ({"sign_rule": "single study"}, "sign rule"), ({"sha256": "x" * 64}, "bytes equal"),
                      ({"sha256": "2" * 64}, "`two`'s")):
        pre = _drops(sc3_rec={**rec, **over})
        assert why in pre["dropped"]["sc3"], (over, pre["dropped"])
        assert _arms_left(pre) == BULK_ONLY
    # |corr(ln phi_M, ln T)| at 0.5 or more
    s = _synth()
    loaded = s.loaded
    T = np.where(s.sc3 != 1.0, np.exp(0.05 * np.log(s.phi)), 1.0)
    loaded["tables"]["sc3"] = (T, loaded["tables"]["sc3"][1])
    pre = k24.pregate_k24(s.phi, loaded, s.bases, s.den)
    assert "corr(ln phi_M, ln T)" in pre["dropped"]["sc3"]


@pytest.mark.parametrize("case", ["at_0.9", "below_0.9", "undefined"])
def test_the_corr_with_two_drop_takes_the_joint(case, monkeypatch):
    """Restated 2026-11-20: a corr of 0.9 or more (the repeat clause) drops x_sc3's three arms only and the joint keeps
    its own drops (here sc3 = two ** 1.3 is uncorrelated with bulk, so the joint's window drops it on its own); an
    undefined corr still takes the joint."""
    s0 = _synth()
    if case == "at_0.9":
        sc3 = s0.two ** 1.3  # corr 1 with `two`
        pre = _drops(sc3=sc3)
        assert "repeat of the named form" in pre["dropped"]["sc3"] and pre["sc3_repeat_only"] is True
        assert "sc3 is dropped" not in pre.get("dropped", {}).get("joint", "")
        if "joint" in pre["dropped"]:
            assert _arms_left(pre) == BULK_ONLY
        else:
            assert _arms_left(pre) == BULK_ONLY + k24.ARMS_OF_KIND["joint"]
    elif case == "below_0.9":
        pre = _drops()
        assert pre["corr"]["sc3_two"] < 0.9 and "sc3" not in pre["dropped"] and "joint" not in pre["dropped"]
    else:
        original = ext.corr_between
        monkeypatch.setattr(ext, "corr_between", lambda a, b: None if b is s0.two else original(a, b))
        pre = k24.pregate_k24(s0.phi, s0.loaded, s0.bases, s0.den)
        assert "corr(ln sc3, ln two) undefined" in pre["dropped"]["sc3"] and "joint" in pre["dropped"]
        assert pre["sc3_repeat_only"] is False and "sc3 is dropped" in pre["dropped"]["joint"]


@pytest.mark.parametrize("case", ["negative", "at_0.9", "undefined"])
def test_the_joint_window(case, monkeypatch):
    s0 = _synth()
    if case == "negative":
        bulk = np.where(s0.bulk != 1.0, 1.0 / np.where(s0.sc3 != 1.0, s0.sc3, s0.bulk), 1.0)
        pre = _drops(bulk=bulk)
        assert pre["corr"]["sc3_bulk"] < 0
    elif case == "at_0.9":
        bulk = np.where(s0.sc3 != 1.0, s0.sc3 ** 1.7, 1.0)
        pre = _drops(bulk=bulk)
        assert pre["corr"]["sc3_bulk"] >= 0.9
    else:
        original = ext.corr_between
        monkeypatch.setattr(ext, "corr_between", lambda a, b: None if b is s0.bulk else original(a, b))
        pre = k24.pregate_k24(s0.phi, s0.loaded, s0.bases, s0.den)
    assert "corr(ln sc3, ln bulk)" in pre["dropped"]["joint"]
    assert "sc3" not in pre["dropped"] and "bulk" not in pre["dropped"]
    assert _arms_left(pre) == k24.ARMS[:8] and k24.try_count(_arms_left(pre)) == 6


def test_drops_as_arm_drops_and_the_count():
    assert k24.active_arms({"dropped": {}}) == k24.ARMS
    assert k24.active_arms({"dropped": {"joint": "x"}}) == k24.ARMS[:8]
    assert k24.active_arms({"dropped": {"bulk": "x"}}) == SC3_ONLY  # the joint goes with either
    assert k24.active_arms({"dropped": {"sc3": "x"}}) == BULK_ONLY
    assert k24.active_arms({"dropped": {"sc3": "x", "bulk": "y"}}) == ("barred_form", "x_e2")
    # the restated repeat clause: sc3 dropped by it alone leaves the joint in
    rep = ("barred_form", "x_e2", "x_bulk", "x_bulk_perm", "x_bulk_size", "x_joint", "x_joint_perm", "x_joint_size")
    assert k24.active_arms({"dropped": {"sc3": "x"}, "sc3_repeat_only": True}) == rep
    assert k24.active_arms({"dropped": {"sc3": "x", "bulk": "y"}, "sc3_repeat_only": True}) == ("barred_form", "x_e2")
    for r in range(len(k24.ARMS) + 1):
        for sub in itertools.combinations(k24.ARMS, r):
            assert k24.try_count(sub) <= 9


# ------------------------------------------------------------------------------------------------ (4) the tables


def _write_table(d: Path, name: str, genes, vals, man: dict) -> str:
    fd = d / name / "factors"
    fd.mkdir(parents=True, exist_ok=True)
    p = fd / f"{name}.csv"
    pd.DataFrame({"gene": list(genes), "factor": vals}).to_csv(p, index=False)
    sha = _sha(p)
    (fd / "factors.json").write_text(json.dumps({"files": {name: {"file": str(p), "sha256": sha}}}))
    (d / name / f"{name}_manifest.json").write_text(json.dumps({**man, "name": name,
                                                                 "files": {name: {"file": str(p), "sha256": sha}}}))
    return sha


def _git_env():
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                        "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, env=_git_env(),
                       text=True)
    return r.stdout.strip()


@pytest.fixture
def folders(tmp_path, monkeypatch):
    """Two tables folders in combine.py's layout, the listing committed in a throwaway repository, the pins patched."""
    s = _synth(n=30)
    genes = [f"G{i}" for i in range(30)]
    d1, d2 = tmp_path / "tables", tmp_path / "tables_two"
    sc3_man = {"studies": ["GSE204908", "GSE277839", "GSE299068"], "tiers": ["unrelated"], "sign_rule": "applied",
               "refstats_sha256": "r" * 64, "checked": "x"}
    bulk_man = {k: v for k, v in _bulk_rec(30, s.shares).items() if k != "sha256"}
    shas = {"sc3": _write_table(d1, "all", genes, s.sc3, sc3_man),
            "bulk": _write_table(d1, "bulk_all", genes, s.bulk, bulk_man),
            "two": _write_table(d2, "all", genes, s.two, {**sc3_man, "studies": ["GSE204908", "GSE277839"]})}
    (d1 / "COMBINE.json").write_text(json.dumps({"tables": {"all": {"all": {"sha256": shas["sc3"]}}}, "record": {}}))
    (d1 / "COMBINE_bulk.json").write_text(json.dumps({"tables": {"bulk_all": {"bulk_all": {"sha256": shas["bulk"]}}},
                                                      "record": {}}))
    (d2 / "COMBINE.json").write_text(json.dumps({"tables": {"all": {"all": {"sha256": shas["two"]}}}, "record": {}}))
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    lp = repo / "TABLES_K24.sha256"
    lp.write_text("".join(f"{shas[k]}  {k}\n" for k in k24.KEYS))
    _git(repo, "add", "TABLES_K24.sha256")
    _git(repo, "commit", "-q", "-m", "listing")
    monkeypatch.setattr(k24, "LISTING_FILE", lp)
    monkeypatch.setattr(k24, "K24_LISTING_SHA256", _sha(lp))
    monkeypatch.setattr(k24, "K24_TABLE_SHA256", dict(shas))
    monkeypatch.setattr(k24, "check_pins", lambda pin, pins: None if pin and all(pins.values()) else (
        _ for _ in ()).throw(k24.Refused("not filled")))
    env = {k24.ENV_TABLES: str(d1), k24.ENV_TWO: str(d2)}
    return SimpleNamespace(d1=d1, d2=d2, genes=genes, shas=shas, env=env, lp=lp, repo=repo, s=s)


def test_the_tables_load_by_key_with_their_records(folders):
    f = folders
    loaded = k24.load_tables(f.genes, env=f.env)
    assert set(loaded["tables"]) == set(k24.KEYS)
    T, rec = loaded["tables"]["bulk"]
    assert np.allclose(T, f.s.bulk, rtol=1e-15, atol=0) and rec["combine_file"] == "COMBINE_bulk.json"
    assert rec["file"] == "<VEC_T3_EXT_TABLES>/bulk_all/factors/bulk_all.csv" and rec["sha256"] == f.shas["bulk"]
    assert rec["modality"] == "bulk" and rec["sign_rule"] == "single study" and len(rec["libraries_by_unit"]) == 25
    assert loaded["tables"]["two"][1]["file"] == "<VEC_T3_EXT_TABLES_TWO>/all/factors/all.csv"
    assert loaded["tables"]["sc3"][1]["modality"] is None and loaded["tables"]["sc3"][1]["studies"][2] == "GSE299068"
    er = k24.ext_record(loaded)
    assert er["tables"] == f.shas and set(er["combine_sha256"]) == set(k24.KEYS)
    assert k24.listing_text(f.env) == f.lp.read_text()
    assert k24.pins_text(f.lp).startswith(f'K24_LISTING_SHA256: str | None = "{_sha(f.lp)}"')


def test_the_loader_refusals(folders, monkeypatch, tmp_path):
    f = folders
    # pins not filled
    monkeypatch.setattr(k24, "K24_TABLE_SHA256", {**f.shas, "bulk": None})
    with pytest.raises(k24.Refused, match="not filled"):
        k24.load_tables(f.genes, env=f.env)
    monkeypatch.setattr(k24, "K24_TABLE_SHA256", dict(f.shas))
    # unlisted (a listing of other keys)
    other = f.repo / "other.sha256"
    other.write_text(f"{f.shas['sc3']}  sc3\n{f.shas['two']}  two\n")
    with pytest.raises(k24.Refused, match="not tracked"):
        k24.read_listing(other, _sha(other), f.shas)
    with pytest.raises(k24.Refused, match="not exactly"):
        k24.read_listing(other, _sha(other), f.shas, check_tracked=False)
    with pytest.raises(k24.Refused, match="not one of K24's tables"):
        k24.load_k24_table("unrelated", f.genes, {}, "l", f.shas, f.env)
    # the listing off its pin
    with pytest.raises(k24.Refused, match="not the pinned K24_LISTING_SHA256"):
        k24.read_listing(f.lp, "0" * 64, f.shas)
    # bytes off the pins
    p = f.d1 / "bulk_all" / "factors" / "bulk_all.csv"
    good = p.read_bytes()
    p.write_bytes(good + b"\n")
    with pytest.raises(k24.Refused, match="not the listed"):
        k24.load_tables(f.genes, env=f.env)
    p.write_bytes(good)
    # a folder that is not a real directory: unset, missing, a symlink
    with pytest.raises(k24.Refused, match="VEC_T3_EXT_TABLES_TWO is not set"):
        k24.load_tables(f.genes, env={k24.ENV_TABLES: str(f.d1)})
    with pytest.raises(k24.Refused, match="VEC_T3_EXT_TABLES_TWO does not name a directory"):
        k24.load_tables(f.genes, env={**f.env, k24.ENV_TWO: str(tmp_path / "none")})
    link = tmp_path / "link"
    link.symlink_to(f.d2)
    with pytest.raises(k24.Refused, match="symlink"):
        k24.load_tables(f.genes, env={**f.env, k24.ENV_TWO: str(link)})
    # a COMBINE record of another sha256, and another form
    cb = f.d1 / "COMBINE_bulk.json"
    raw = cb.read_bytes()
    cb.write_text(json.dumps({"tables": {"bulk_all": {"bulk_all": {"sha256": "0" * 64}}}, "record": {}}))
    with pytest.raises(k24.Refused, match="records another sha256"):
        k24.load_tables(f.genes, env=f.env)
    cb.write_text(json.dumps({"bulk_all": {"sha256": f.shas["bulk"]}}))
    with pytest.raises(k24.Refused, match="combine.py's form"):
        k24.load_tables(f.genes, env=f.env)
    cb.write_bytes(raw)
    # genes off the barred factor's order
    with pytest.raises(k24.Refused, match="genes differ"):
        k24.load_tables(f.genes[::-1], env=f.env)
    k24.load_tables(f.genes, env=f.env)


def test_combine_bulk_json_and_libraries_by_unit_are_combine_pys_forms(tmp_path):
    """The bulk table written by main's combine.bulk_tables and its record written as combine() writes
    COMBINE_bulk.json load through K24's loader, and (b2) to (b4) read the manifest's forms."""
    import combine as cb
    import t3x_common as cm

    rng = np.random.default_rng(4)
    n = 60
    genes = [f"G{i}" for i in range(n)]
    ref = {"genes": genes, "groups": ["a", "b"], "shares": {"a": 0.6, "b": 0.4},
           "mu": {"a": rng.random(n) + 0.1, "b": rng.random(n) + 0.1}, "file_sha256": "r" * 64}
    common = rng.normal(0, 0.2, n)
    units = []
    for i in range(14):
        L = common + rng.normal(0, 0.1, n)
        L[:5] = np.nan
        units.append({"accession": "E-MTAB-7322", "gene": f"L{i}", "study": "MGP-DMDD", "tier": "unrelated",
                      "L": {"__whole__": L}, "sha256": "u" * 64, "modality": "bulk", "k": {"E9.5": 0.64}, "m": {},
                      "libraries": {"ko": [f"k{i}{j}" for j in range(3)], "control": [f"c{i}{j}" for j in range(3)]},
                      "guards": {"E9.5": "passed"}})
    out = tmp_path / "tables"
    rec = {"left_out": [], "not_cleared": [], "inputs": {}}
    written = cb.bulk_tables(units, ref, rec, out)
    cm.dump_json({"tables": written, "record": rec}, out / "COMBINE_bulk.json")
    comb, _ = k24.combine_record(out, "COMBINE_bulk.json")
    assert set(comb) == {"bulk_all", "bulk_study_MGP-DMDD"} and comb["bulk_all"] == comb["bulk_study_MGP-DMDD"]
    sha = comb["bulk_all"]
    T, r, _ = k24.load_k24_table("bulk", genes, {"bulk": sha}, "l" * 64, {"bulk": sha}, {k24.ENV_TABLES: str(out)})
    assert r["modality"] == "bulk" and r["studies"] == ["MGP-DMDD"] and r["sign_rule"] == "single study"
    assert r["libraries_by_unit"]["E-MTAB-7322:L0"] == {"ko": ["k00", "k01", "k02"], "control": ["c00", "c01", "c02"]}
    why, figs = k24.bulk_drops(T, r)
    assert figs["lines"] == 14 and figs["libraries_ko"] == 42 and figs["libraries_wt"] == 42
    assert figs["agree_G"] == n - 5 and figs["agree_count"] > 0
    assert any(w.startswith("(b3)") and "42" in w for w in why) is False  # 42 >= 36 in each arm
    assert not any(w.startswith(("(b1)", "(b2)", "(b5)")) for w in why), why


# ------------------------------------------------------------------------------------------------ (5) the fold


@pytest.fixture
def fold(monkeypatch):
    """K18's plain fold fixture (as tests/test_k23_outside_exponents.py), the barred factor and the three tables in a
    file gene order that differs from the reference's."""
    bs = harness._t3_head()
    rng = np.random.default_rng(5)
    C = rng.poisson(3.0, size=(N_POOL, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, 3, N_POOL)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    V = SimpleNamespace(p=np.array([0.5, 0.3, 0.2]))
    ref = bs.Reference(preset="gata4", rules={}, genes=GENES, V=V, pool=bs.t3d.Pool(C, lab, P, 3),
                       coords=rng.normal(size=(N_POOL, 3)))
    pool = ad.AnnData(X=bs.t3d.lognorm(C).astype(np.float32), var=pd.DataFrame(index=GENES))
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, N_POOL)]
    phi = np.exp(rng.normal(0, 0.5, G))
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {}}
    f_s2 = np.array([0.6, 0.25, 0.15])
    barred_file = np.clip(np.exp(rng.normal(0, 0.8, G + 1)), 0.25, 4.0)
    tabs = {k: np.where(rng.random(G + 1) < 0.7, np.exp(rng.normal(0, 0.3, G + 1)), 1.0) for k in k24.KEYS}
    spec, _, _ = k18.barred_spec()
    log = {"construct": [], "assemble": []}

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        assert r is ref and ratio is None and alphas is None
        rows = g.integers(0, N_POOL, n)
        log["construct"].append({"phi": None if phi is None else np.asarray(phi).copy(), "w": w, "row": rows.copy()})
        Cr = C[rows]
        if phi is not None:
            Cr = bs.t3d.apply_gene_factor(Cr, np.asarray(phi), P[rows], g)
        return bs.Built(X=bs.t3d.lognorm(Cr), lab=lab[rows], source=np.full(n, "wt", dtype=object), row=rows)

    def assemble(r, built, panel, g):
        log["assemble"].append({"row": built.row.copy(), "draw": g.random()})
        return built.X, None, {}

    monkeypatch.setattr(k21.k7, "_prepared", lambda ctx: prep)
    monkeypatch.setattr(bs, "construct", construct)
    monkeypatch.setattr(bs, "composition", lambda v, s, params: (f_s2, {"states_missing": []}))
    monkeypatch.setattr(bs, "assemble", assemble)
    monkeypatch.setattr(harness, "count_for", lambda board, n: N_DRAW)
    monkeypatch.setattr(harness, "finish", lambda ctx, p, X, coords: X)
    monkeypatch.setattr(k18, "_BARRED", {"phi": barred_file, "genes": FILE_GENES, "spec": spec, "record": {}})
    consts = {"e_sc3": 1.7, "e_bulk": 0.6, "c": 0.62,
              "twins": {"x_sc3": {"a": 1.02, "k": 0.97}, "x_bulk": {"a": 1.01, "k": 0.99},
                        "x_joint": {"a": 1.015, "k": 0.98}}}
    k18._FOLD.clear()
    k24._K24.clear()
    k24._K24.update(genes=list(FILE_GENES), T=dict(tabs), consts=consts)

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    order = lambda v: k18.in_gene_order(v, FILE_GENES, GENES)  # noqa: E731
    yield SimpleNamespace(bs=bs, ctx=ctx, log=log, pool=pool, phi_barred=order(barred_file),
                          T={k: order(v) for k, v in tabs.items()}, basis=harness.t3_mean_expression(pool),
                          consts=consts, tabs=tabs, spec=spec)
    k18._FOLD.clear()


def test_barred_form_is_k18s_arm_value_for_value_and_the_names_come_back(fold):
    for seed in (0, 7):
        k18._FOLD.clear()
        want = k18._build_plain(fold.ctx("g1", seed), k18.BARRED)
        with k24.k18_arms(k24.ARMS):
            got = k18._build_plain(fold.ctx("g1", seed), k24.BASE)
        assert np.array_equal(got, want)
    assert k18.ARMS == k24._K18_ORIGINAL["ARMS"] and k18.build_arm is k24._K18_ORIGINAL["build_arm"]
    assert k18.fold_factors is k24._K18_ORIGINAL["fold_factors"]
    with pytest.raises(ZeroDivisionError):
        with k24.k18_arms(k24.ARMS):
            assert k18.fold_factors is k24.fold_factors and k18.build_arm is k24.build_arm
            raise ZeroDivisionError
    assert k18.fold_factors is k24._K18_ORIGINAL["fold_factors"] and k18.ARMS == k24._K18_ORIGINAL["ARMS"]
    with pytest.raises(RuntimeError):
        with k24.k18_arms(("x_sc3",)):
            pass


def test_fold_factors_is_k22s_body_on_two(fold):
    k22._EXT.clear()
    k22._EXT.update(genes=list(FILE_GENES), T={"all": fold.tabs["two"]})
    with k22.k18_arms(("barred_form", "x_all")):
        k18._build_plain(fold.ctx("g1", 0), "barred_form")
        F22 = k18._FOLD["F"]
    with k24.k18_arms(k24.GUARDS):
        k18._build_plain(fold.ctx("g1", 0), "barred_form")
        F24 = k18._FOLD["F"]
    assert np.array_equal(F24["ext"]["two"], F22["ext"]["all"])
    assert np.array_equal(F24["ext_perm"]["two"], F22["ext_perm"]["all"])
    for k in ("phi_barred", "phi", "f_upload", "ones"):
        assert np.array_equal(F24[k], F22[k]), k
    assert set(F24["scaled"]) == {"sc3", "bulk", "joint"}


def test_every_arm_draws_the_barred_s1s_cells_and_x_e2_draws_k23s(fold):
    k22._EXT.clear()
    k22._EXT.update(genes=list(FILE_GENES), T={"all": fold.tabs["two"]}, record={})
    k23._K23.update({"2": {"a": 1.0, "k": 1.0}, "4": {"a": 1.0, "k": 1.0}})
    with k23.k18_arms(("barred_form", "x_e2")):
        fold.log["construct"].clear()
        x23 = k18._build_plain(fold.ctx("g2", 3), "x_e2")
    T, c = fold.T, fold.consts
    scaled = {"sc3": T["sc3"] ** c["e_sc3"], "bulk": T["bulk"] ** c["e_bulk"],
              "joint": T["sc3"] ** (c["c"] * c["e_sc3"]) * T["bulk"] ** (c["c"] * c["e_bulk"])}
    perm = {k: harness.t3_decile_permutation(v, fold.basis, 20261113, 10) for k, v in scaled.items()}
    pb = fold.phi_barred
    want = {"barred_form": pb, "x_e2": pb * T["two"] ** 2}
    for kind in ("sc3", "bulk", "joint"):
        tw = c["twins"][f"x_{kind}"]
        want |= {f"x_{kind}": pb * scaled[kind], f"x_{kind}_perm": pb * perm[kind],
                 f"x_{kind}_size": tw["k"] * pb ** tw["a"]}
    X = {}
    with k24.k18_arms(k24.ARMS):
        for arm in k24.ARMS:
            fold.log["construct"].clear()
            fold.log["assemble"].clear()
            X[arm] = k18._build_plain(fold.ctx("g2", 3), arm)
            (cc,) = fold.log["construct"]
            assert cc["w"] == 1.0 and np.allclose(cc["phi"], want[arm], rtol=1e-14, atol=0), arm
            fold.log[arm] = (cc["row"], fold.log["assemble"][0])
    assert np.array_equal(X["x_e2"], x23)
    base_rows, base_asm = fold.log[k24.BASE]
    for arm in k24.ARMS:
        rows, asm = fold.log[arm]
        assert np.array_equal(rows, base_rows) and asm["draw"] == base_asm["draw"]
        if arm != k24.BASE:
            assert not np.array_equal(X[arm], X[k24.BASE])


def test_check_world_compares_the_run_s_controls_with_the_pregate_s():
    s = _synth(3)
    pre = k24.pregate_k24(s.phi, s.loaded, s.bases, s.den)
    k24._K24["consts"] = k24.constants_of(pre)
    tabs = {k: v[0] for k, v in s.loaded["tables"].items()}
    F = {"phi_barred": s.phi, "ext": tabs}
    F["scaled"] = {k: k24.k24_table(k, tabs, k24._K24["consts"]) for k in ("sc3", "bulk", "joint")}
    F["scaled_perm"] = {k: ext.permuted(v, s.bases["plain"]["basis"]) for k, v in F["scaled"].items()}
    k24.check_world(pre, "plain", F, k24.ARMS)
    with pytest.raises(k24._Mismatch):
        k24.check_world(pre, "dm", F, k24.ARMS)


@pytest.mark.parametrize("flat", ["bulk", "sc3"])
def test_a_table_with_no_exponent_drops_with_the_joint_and_the_other_table_s_arms_are_still_made(fold, flat):
    """An all-ones table has zero log SD and no exponent: the pre-gate drops it and the joint and keeps the other
    table; fold_factors still makes the other table's scaled table and its permutation, so its try, control and twin
    build (review blocker of 7304c05: one missing exponent used to stop every kind)."""
    other = "sc3" if flat == "bulk" else "bulk"
    s0 = _synth()
    pre = _drops(**{flat: np.ones(s0.n)})
    assert pre["e"][flat] is None and pre["e"][other] is not None
    assert "no exponent" in pre["dropped"][flat] and "joint" in pre["dropped"] and other not in pre["dropped"]
    arms = k24.active_arms(pre)
    assert arms == (SC3_ONLY if other == "sc3" else BULK_ONLY) and k24.try_count(arms) == 3
    consts = k24.constants_of(pre)
    assert consts[f"e_{flat}"] is None and consts["c"] is None and set(consts["twins"]) == {f"x_{other}"}
    tabs = dict(fold.tabs)
    tabs[flat] = np.ones(G + 1)
    k24._K24.update(T=tabs, consts=consts)
    with k24.k18_arms(arms):
        k18._build_plain(fold.ctx("g1", 0), k24.BASE)
        F = k18._FOLD["F"]
        assert set(F["scaled"]) == {other} and set(F["scaled_perm"]) == {other}
        pb = fold.phi_barred
        scaled = fold.T[other] ** consts[f"e_{other}"]
        perm = harness.t3_decile_permutation(scaled, fold.basis, 20261113, 10)
        tw = consts["twins"][f"x_{other}"]
        want = {f"x_{other}": pb * scaled, f"x_{other}_perm": pb * perm, f"x_{other}_size": tw["k"] * pb ** tw["a"]}
        for arm, v in want.items():
            assert np.allclose(k24.arm_factor(arm, F), v, rtol=1e-14, atol=0), arm
            fold.log["construct"].clear()
            k18._build_plain(fold.ctx("g1", 0), arm)
            (cc,) = fold.log["construct"]
            assert np.allclose(cc["phi"], v, rtol=1e-14, atol=0), arm
        for arm in (f"x_{flat}", f"x_{flat}_perm", "x_joint", "x_joint_perm"):
            with pytest.raises(RuntimeError, match="not made"):
                k24.arm_factor(arm, F)


@pytest.mark.parametrize("flat", ["bulk", "sc3"])
@pytest.mark.parametrize("value", [1.2, 0.25, 4.0, 0.7])
def test_a_table_of_one_constant_value_has_no_exponent_and_the_other_table_is_still_read(flat, value):
    """A table holding one value that is not 1 has log SD exactly 0 by the rule, while np.std may return rounding
    noise for it (the second review's blocker: the exponent then came out near 1e15, the scaled table overflowed and
    the pre-gate raised, so nothing ran). It has no exponent: it drops with the joint, the other table keeps its
    three tries, and nothing raises."""
    other = "sc3" if flat == "bulk" else "bulk"
    s0 = _synth()
    const = np.full(s0.n, value)
    assert k24.exponent(const) is None and k24.exponent(np.full(300, value)) is None
    assert k24.joint_c(const, const, 1.0, 1.0) is None
    pre = _drops(**{flat: const})
    assert pre["e"][flat] is None and pre["e"][other] is not None
    assert "no exponent" in pre["dropped"][flat] and "joint" in pre["dropped"] and other not in pre["dropped"]
    arms = k24.active_arms(pre)
    assert arms == (SC3_ONLY if other == "sc3" else BULK_ONLY) and k24.try_count(arms) == 3
    consts = k24.constants_of(pre)
    assert consts[f"e_{flat}"] is None and consts["c"] is None and set(consts["twins"]) == {f"x_{other}"}


def test_the_exponent_of_a_table_that_is_not_constant_is_unchanged_by_the_constant_guard():
    rng = np.random.default_rng(7)
    for sd in (1e-9, 1e-6, 1e-3, 0.02, 0.3):
        T = np.exp(rng.normal(0.0, sd, 500))
        e = k24.exponent(T)
        assert e == k24.S_STAR / float(np.std(np.log(T)))
        assert math.isclose(float(np.std(np.log(T ** e))), k24.S_STAR, rel_tol=1e-12)
    two = np.exp(rng.normal(0.0, 0.02, 500))
    assert k24.joint_c(two, two[::-1].copy(), 1.5, 2.5) == k24.S_STAR / float(
        np.std(1.5 * np.log(two) + 2.5 * np.log(two[::-1])))


@pytest.mark.parametrize("flat", ["bulk", "sc3"])
@pytest.mark.parametrize("v,d", [(1.2, 1e-6), (4.0, 1e-3), (0.25, 1e-12), (1.2, 2.220446049250313e-16)])
def test_a_nearly_constant_table_drops_on_s1_before_its_scaled_table_is_formed(flat, v, d):
    """A table of one value with one gene a hair off it has a tiny log SD and a huge finite exponent: the rule drops
    it on (s1). The third review's blocker: the code formed T ** e_T first, which overflowed or underflowed, and the
    pre-gate raised instead of dropping. Now (s1) is checked before T ** e_T is formed."""
    other = "sc3" if flat == "bulk" else "bulk"
    s0 = _synth()
    T = np.full(s0.n, v)
    T[7] = v + d
    if float(np.ptp(np.log(T))) == 0.0:   # the hair vanished in rounding: the constant-table case, tested above
        pytest.skip("not representable")
    pre = _drops(**{flat: T})
    assert pre["e"][flat] > k24.E_MAX and "(s1)" in pre["dropped"][flat]
    assert "joint" in pre["dropped"] and other not in pre["dropped"]
    assert flat not in pre["scaled"] and f"x_{flat}" not in pre["arms"]
    arms = k24.active_arms(pre)
    assert arms == (SC3_ONLY if other == "sc3" else BULK_ONLY) and k24.try_count(arms) == 3


@pytest.mark.parametrize("flat", ["bulk", "sc3"])
@pytest.mark.parametrize("mean", [100.0, -100.0])
def test_a_table_far_from_one_drops_on_s3_before_its_scaled_table_is_formed(flat, mean):
    """Every gene's ln T near +100 or -100 with a small spread: the exponent is inside [1/8, 8], T ** e_T overflows
    or underflows, and the rule drops the table on (s3). The fourth review's blocker: the pre-gate raised instead."""
    other = "sc3" if flat == "bulk" else "bulk"
    s0 = _synth()
    T = np.exp(mean + np.random.default_rng(3).normal(0.0, k24.S_STAR / 7.9, s0.n))
    e = k24.exponent(T)
    assert k24.E_MIN <= e <= k24.E_MAX and e * float(np.max(np.abs(np.log(T)))) > 700.0
    pre = _drops(**{flat: T})
    assert "(s3)" in pre["dropped"][flat] and "(s1)" not in pre["dropped"][flat]
    assert flat not in pre["scaled"] and f"x_{flat}" not in pre["arms"]
    assert "joint" in pre["dropped"] and other not in pre["dropped"]
    assert k24.try_count(k24.active_arms(pre)) == 3


def test_the_joint_is_formed_through_logs_when_a_factor_alone_would_overflow():
    """A large c (the two size-matched tables nearly cancel): sc3 ** (c e_sc3) overflows and bulk ** (c e_bulk)
    underflows although J is small. The fourth review's blocker: the product was formed factor by factor and the
    pre-gate raised. For ordinary tables the form as written is used, value for value."""
    rng = np.random.default_rng(4)
    n = 300
    half = rng.random(n) < 0.5
    z = rng.normal(0.0, 1e-6, n)
    y = 0.6 * z + 0.8 * rng.normal(0.0, 1e-6, n)
    sc3 = np.where(half, np.exp(0.05 + z), 1.0)
    bulk = np.where(half, np.exp(-0.05 + y), 1.0)
    e_s, e_b = k24.exponent(sc3), k24.exponent(bulk)
    c = k24.joint_c(sc3, bulk, e_s, e_b)
    assert c > 1000.0
    consts = {"e_sc3": e_s, "e_bulk": e_b, "c": c}
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        naive = sc3 ** (c * e_s) * bulk ** (c * e_b)
    assert not np.all(np.isfinite(naive)) or np.any(naive == 0.0)
    J = k24.k24_table("joint", {"sc3": sc3, "bulk": bulk}, consts)
    assert np.all(np.isfinite(J)) and np.all(J > 0)
    assert np.allclose(np.log(J), c * e_s * np.log(sc3) + c * e_b * np.log(bulk), rtol=1e-12, atol=1e-15)
    assert math.isclose(float(np.std(np.log(J))), k24.S_STAR, rel_tol=1e-6)
    # ordinary tables: the form as written, value for value
    a = np.exp(rng.normal(0.0, 0.03, n))
    b = np.exp(rng.normal(0.0, 0.02, n))
    ea, eb = k24.exponent(a), k24.exponent(b)
    cc = k24.joint_c(a, b, ea, eb)
    got = k24.k24_table("joint", {"sc3": a, "bulk": b}, {"e_sc3": ea, "e_bulk": eb, "c": cc})
    assert np.array_equal(got, a ** (cc * ea) * b ** (cc * eb))


def test_a_share_that_is_not_a_number_drops_the_bulk_table_and_a_guards_list_reports_nothing():
    n, U = 200, 14
    for bad in ("x", {"a": 1}, [0.9], True):
        why, f = _b4(n, U, [bad] + [0.8] * (n - 1))
        assert _b4_fails(why) and any("not a number" in w for w in why), bad
        assert f["n_g_source"] == "malformed" and f["agree_G"] is None
    assert k24.somite_report([1, 2, 3]) == k24.somite_report(None)


# ------------------------------------------------------------------------------------------------ (6) the rule


def _rows(scores: dict, gated: dict | None = None, noise: dict | None = None) -> list[dict]:
    gated = gated or {}
    out = []
    for s in SEEDS:
        sc = {a: v + 0.1 * s + (noise[a][s] if noise and a in noise else 0.0) for a, v in scores.items()}
        out.append({"seed": s, "scores": sc,
                    "skills": {a: {"de_score": 0.5, "severity_slope": 0.6} for a in sc},
                    "t3": {a: {"severity": {"state": "r2_gate" if s < gated.get(a, 0) else "scored", "beta": 0.4,
                                            "r2": 0.05}} for a in sc}})
    return out


BASE_SCORE = 56.0
ALL_CELLS = [(w, n) for w in k24.WORLDS for n in k24.TASKS]


def _table(**arms) -> dict:
    return {harness.FLOOR: 50.0, **{a: BASE_SCORE for a in k24.ARMS}, **arms}


def _per(cells: dict, gated: dict | None = None, arms=k24.ARMS) -> dict:
    out: dict = {}
    for w in k24.WORLDS:
        a = (harness.FLOOR, *arms) if w == "plain" else arms
        out[w] = {n: k24._summarise(_rows({x: cells.get((w, n), _table())[x] for x in a},
                                          (gated or {}).get((w, n))), a) for n in k24.TASKS}
    return out


def _everywhere(**arms) -> dict:
    return {k: _table(**arms) for k in ALL_CELLS}


GOOD = {"x_sc3": 57.5, "x_e2": 57.0}  # +1.5 over barred_form, its control and its twin; +0.5 over x_e2


def test_the_summaries_carry_the_extra_pairs():
    p = _per(_everywhere(**GOOD))["plain"]["g1"]["paired"]
    for k in ("x_sc3 - barred_form", "x_sc3 - x_sc3_perm", "x_sc3 - x_sc3_size", "x_sc3 - x_e2", "x_bulk - x_e2",
              "x_joint - x_joint_size", "x_e2 - barred_form", "x_joint_perm - barred_form", "barred_form - floor"):
        assert k in p, k
    assert p["x_sc3 - x_e2"]["mean"] == pytest.approx(0.5) and p["x_bulk - x_e2"]["mean"] == pytest.approx(-1.0)


def test_one_qualified_arm_is_named():
    res = k24.decide(_per(_everywhere(**GOOD)), k24.ARMS)
    s = res["status"]["x_sc3"]
    assert s["qualified"] and all(s[c] for c in "abcdef") and res["names"] == ["x_sc3"]
    assert res["files"] == ["k24sc3-B (base B of decision log 2026-11-20 (q))"]
    assert "build_k24_form.py --form k24sc3-<B>" in res["outcome"] and res["rule"] == "one arm qualified"
    assert not res["status"]["x_bulk"]["qualified"] and not res["status"]["x_joint"]["qualified"]


def _fails_alone(per, letter, t="x_sc3"):
    res = k24.decide(per, k24.ARMS)
    s = res["status"][t]
    assert not s["qualified"] and not s[letter], s
    assert all(s[c] for c in "abcdef" if c != letter), s
    assert res["names"] == []
    return res


def test_each_clause_failing_alone():
    up = _table(**GOOD)
    # (a): AHEAD in plain g1 and dm g2 only, with one more gated seed than barred_form in dm g2 (x_e2 too, so (f) holds)
    _fails_alone(_per({("plain", "g1"): up, ("dm", "g2"): up}, {("dm", "g2"): {"x_sc3": 1, "x_e2": 1}}), "a")
    assert k24.decide(_per({("plain", "g1"): up, ("dm", "g2"): up}), k24.ARMS)["names"] == ["x_sc3"]
    # (b): behind in one cell
    _fails_alone(_per(_everywhere(**GOOD) | {("dm", "g1"): _table(x_sc3=55.5, x_e2=55.4)}), "b")
    # (c): two more gated seeds than barred_form in plain g1
    _fails_alone(_per(_everywhere(**GOOD), {("plain", "g1"): {"x_sc3": 2, "x_e2": 2}}), "c")
    # (d): its permuted control as high in one AHEAD cell
    _fails_alone(_per(_everywhere(**GOOD) | {("plain", "g2"): _table(**GOOD, x_sc3_perm=57.4)}), "d")
    # (d): no dm cell where it beats its control
    cells = _everywhere(**GOOD) | {("dm", n): _table(**GOOD, x_sc3_perm=57.4) for n in k24.TASKS}
    _fails_alone(_per(cells), "d")
    # (e): its size twin as high in one AHEAD cell
    _fails_alone(_per(_everywhere(**GOOD) | {("dm", "g2"): _table(**GOOD, x_sc3_size=57.3)}), "e")


def test_f_is_the_head_to_head_and_a_to_e_alone_is_reported_not_named():
    # (f) X - x_e2 under the bar everywhere
    res = _fails_alone(_per(_everywhere(x_sc3=57.5, x_e2=57.3)), "f")
    assert res["status"]["x_sc3"]["qualified_a_to_e"] and res["reported"] == {"x_sc3": k24.NOT_AHEAD}
    assert res["outcome"] == k24.NO_FILE == ("no file: no wider table beats the named form at a like size; the K23 "
                                             "file stands alone")
    # (f) ahead of x_e2 in plain cells only
    cells = _everywhere(**GOOD) | {("dm", n): _table(x_sc3=57.5, x_e2=57.4) for n in k24.TASKS}
    _fails_alone(_per(cells), "f")
    # (f) ahead of x_e2 in a dm cell only where the gated counts differ
    _fails_alone(_per(_everywhere(**GOOD), {("dm", "g1"): {"x_e2": 1}, ("dm", "g2"): {"x_e2": 1}}), "f")
    # ... and with one dm cell where they are equal it holds
    assert k24.decide(_per(_everywhere(**GOOD), {("dm", "g1"): {"x_e2": 1}}), k24.ARMS)["names"] == ["x_sc3"]
    # (f) below x_e2 by more than 0.3 in one cell
    _fails_alone(_per(_everywhere(**GOOD) | {("plain", "g2"): _table(x_sc3=57.5, x_e2=58.0)}), "f")


def test_several_qualified_named_by_the_gate_equal_cells():
    cells = _everywhere(x_sc3=57.5, x_bulk=57.8, x_e2=57.0)
    cells[("plain", "g1")] = _table(x_sc3=59.0, x_bulk=57.8, x_e2=57.0)
    # x_bulk has one more gated seed than barred_form in plain g1: that cell leaves the ranking
    res = k24.decide(_per(cells, {("plain", "g1"): {"x_bulk": 1}}), k24.ARMS)
    assert res["qualified"] == ["x_sc3", "x_bulk"] and res["names"] == ["x_bulk"]
    assert res["rank"]["gate_equal_cells"] and "plain g1" not in res["rank"]["cells"] and len(res["rank"]["cells"]) == 3
    # without the gate the four-cell mean would favour x_sc3
    res4 = k24.decide(_per(cells), k24.ARMS)
    assert res4["names"] == ["x_sc3"] and len(res4["rank"]["cells"]) == 4


def test_several_qualified_with_no_gate_equal_cell_use_the_four_cells():
    cells = _everywhere(x_sc3=57.5, x_bulk=57.8, x_e2=57.0)
    gated = {k: {"x_e2": 1} for k in ALL_CELLS}
    gated[("dm", "g1")] = {"x_e2": 1, "x_sc3": 1, "x_bulk": 1}
    res = k24.decide(_per(cells, gated), k24.ARMS)
    assert res["qualified"] == ["x_sc3", "x_bulk"] and res["names"] == ["x_bulk"]
    assert not res["rank"]["gate_equal_cells"] and len(res["rank"]["cells"]) == 4
    assert "no gate-equal cell" in res["rule"]


def test_exact_ties_follow_the_order_sc3_bulk_joint():
    res = k24.decide(_per(_everywhere(x_sc3=57.5, x_bulk=57.5, x_joint=57.5, x_e2=57.0)), k24.ARMS)
    assert res["qualified"] == ["x_sc3", "x_bulk", "x_joint"] and res["names"] == ["x_sc3"]
    res = k24.decide(_per(_everywhere(x_bulk=57.5, x_joint=57.5, x_e2=57.0)), k24.ARMS)
    assert res["names"] == ["x_bulk"]
    res = k24.decide(_per(_everywhere(x_joint=57.5, x_e2=57.0)), k24.ARMS)
    assert res["names"] == ["x_joint"] and res["files"] == ["k24joint-B (base B of decision log 2026-11-20 (q))"]


def test_nothing_qualifies_and_controls_twins_and_guards_are_never_named():
    res = k24.decide(_per(_everywhere(x_e2=60.0, x_sc3_perm=60.0, x_bulk_size=60.0, x_joint_perm=60.0)), k24.ARMS)
    assert res["names"] == [] and res["outcome"] == k24.NO_FILE and set(res["status"]) == set(k24.NAMEABLE)
    arms = SC3_ONLY
    res = k24.decide(_per(_everywhere(**GOOD), arms=arms), arms)
    assert res["status"]["x_bulk"] == {"ran": False, "qualified": False, "qualified_a_to_e": False}
    assert res["names"] == ["x_sc3"]
    rng = np.random.default_rng(21)
    for _ in range(100):
        cells = {k: _table(**{a: float(BASE_SCORE + rng.choice([-0.5, 0.0, 0.5, 1.0, 2.0]))
                              for a in k24.ARMS if a != k24.BASE}) for k in ALL_CELLS}
        res = k24.decide(_per(cells), k24.ARMS)
        assert len(res["names"]) <= 1 and set(res["names"]) <= set(k24.NAMEABLE)


def test_on_null_seeds_a_name_is_rare():
    """Every arm at x_e2's level, per-seed noise SD 0.25 (K23's paired SDs 0.17 to 0.30 around x_e2): the rule names a
    file in under 2% of 300 draws."""
    rng = np.random.default_rng(23)
    named, n = 0, 300
    for _ in range(n):
        per = {}
        for w in k24.WORLDS:
            a = (harness.FLOOR, *k24.ARMS) if w == "plain" else k24.ARMS
            per[w] = {}
            for t in k24.TASKS:
                noise = {x: rng.normal(0, 0.25, len(SEEDS)) for x in k24.ARMS if x != k24.BASE}
                per[w][t] = k24._summarise(_rows({x: _table(**{y: 57.0 for y in k24.ARMS[1:]})[x] for x in a},
                                                 noise=noise), a)
        named += bool(k24.decide(per, k24.ARMS)["names"])
    assert named / n < 0.02, named


# ------------------------------------------------------------------------------------------------ (7) the guard part


def _guard_refs(base=None, e2=None, e2_gated=None):
    base, e2 = base or {}, e2 or {}
    refs = {"k18": {}, "k23_rows": {}, "k23_cells": {}}
    for w in k24.WORLDS:
        for n in k24.TASKS:
            rows = _rows({k24.BASE: BASE_SCORE, "x_e2": 57.0})
            refs["k18"][(w, n)] = [{"seed": r["seed"], "scores": {k24.BASE: r["scores"][k24.BASE]
                                                                  + base.get((w, n), 0.0)}} for r in rows]
            refs["k23_rows"][(w, n)] = [{"seed": r["seed"], "scores": {"x_e2": r["scores"]["x_e2"]
                                                                       + e2.get((w, n), 0.0)}} for r in rows]
            refs["k23_cells"][(w, n)] = {"mean": float(np.mean([q["scores"]["x_e2"] for q in rows])),
                                         "severity_gated": (e2_gated or {}).get((w, n), 0)}
    return refs


def _guard_rows(shift=None, one_seed=None, gated=None):
    out = {}
    for w in k24.WORLDS:
        for n in k24.TASKS:
            rows = _rows({k24.BASE: BASE_SCORE, "x_e2": 57.0}, (gated or {}).get((w, n)))
            for arm, d in (shift or {}).get((w, n), {}).items():
                for r in rows:
                    r["scores"][arm] += d
            for arm, d in (one_seed or {}).get((w, n), {}).items():
                rows[3]["scores"][arm] += d
            out[(w, n)] = rows
    return out


@pytest.mark.parametrize("kw,failed", [
    ({}, []),
    ({"shift": {("plain", "g1"): {"barred_form": 0.11}}}, ["barred_form"]),
    ({"one_seed": {("dm", "g2"): {"barred_form": 0.6}}}, ["barred_form"]),
    ({"shift": {("dm", "g1"): {"x_e2": -0.11}}}, ["x_e2"]),
    ({"one_seed": {("plain", "g2"): {"x_e2": 0.6}}}, ["x_e2"]),
    ({"gated": {("dm", "g1"): {"x_e2": 1}}}, ["x_e2"]),
    ({"shift": {("plain", "g1"): {"barred_form": 0.2, "x_e2": 0.2}}}, ["barred_form", "x_e2"]),
])
def test_guard_1_and_guard_2(kw, failed):
    res = k24.guard_eval(_guard_rows(**kw), _guard_refs())
    assert res["failed"] == failed and res["passes"] == (not failed)
    # within the tolerance it holds
    ok = k24.guard_eval(_guard_rows(shift={("plain", "g1"): {"barred_form": 0.09, "x_e2": 0.09}}), _guard_refs())
    assert ok["passes"]
    # guard 2 against K23's committed mean too
    refs = _guard_refs()
    refs["k23_cells"][("dm", "g2")]["mean"] += 0.11
    assert k24.guard_eval(_guard_rows(), refs)["failed"] == ["x_e2"]


def test_k23s_per_seed_rows_are_read_at_their_sha256(tmp_path):
    runs = tmp_path / "k23run"
    files = {}
    for w in k24.WORLDS:
        for n in k24.TASKS:
            for s in SEEDS:
                f = f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json"
                row = {"seed": s, "scores": {"x_e2": 57.0 + s}, "t3": {"x_e2": {"severity": {"state": "scored"}}}}
                body = {"rows": [row]} if w == "plain" else row
                p = runs / f
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(json.dumps(body))
                files[f] = _sha(p)
    got = k24.k23_seed_rows({"files_sha256": files}, runs)
    assert [q["scores"]["x_e2"] for q in got[("dm", "g2")]] == [57.0 + s for s in SEEDS]
    (runs / "dm" / "g1_seed_3.json").write_text("{}")
    with pytest.raises(k24.Refused, match="not at the sha256"):
        k24.k23_seed_rows({"files_sha256": files}, runs)
    (runs / "dm" / "g1_seed_3.json").unlink()
    with pytest.raises(k24.Refused, match="missing"):
        k24.k23_seed_rows({"files_sha256": files}, runs)
    assert k24.k23_runs_default({"commit": "abcdef0123"}) == REPO / "runs" / "k23_abcdef0"


def test_k23s_run_folder_is_this_checkout_s_and_no_option_overrides_it(monkeypatch, capsys):
    import inspect

    assert "k23_runs" not in inspect.signature(k24.guard_refs).parameters
    assert "no option overrides it" in k24.__doc__
    for argv in (["pregate", "--out", "x", "--k23-runs", "y"],
                 ["run", "--part", "guard", "--world", "dm", "--out", "x", "--k23-runs", "y"],
                 ["guard", "--runs", "x", "--k23-runs", "y"]):
        with pytest.raises(SystemExit) as e:
            k24.main(argv)
        assert e.value.code == 2 and "unrecognized arguments: --k23-runs" in capsys.readouterr().err
    cells = {w: {n: {"arms": {"x_e2": {"mean": 57.0, "severity_gated": 0}}} for n in k24.TASKS} for w in k24.WORLDS}
    monkeypatch.setattr(k24, "k23_committed", lambda d=None: ({}, {"commit": "0123456789ab", "worlds": cells}))
    monkeypatch.setattr(k21, "k18_rows", lambda refs: ({}, []))
    seen = []
    monkeypatch.setattr(k24, "k23_seed_rows", lambda read, runs: seen.append(runs) or {})
    k24.guard_refs()
    assert seen == [REPO / "runs" / "k23_0123456"]


# ------------------------------------------------------------------------------------------------ (8) the reading


def _prov(commit="abc", clean=True, machine="m1"):
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "x"}


EREC = {"listing_sha256": "f" * 64, "tables": {"two": TWO_SHA, "sc3": "3" * 64, "bulk": "b" * 64},
        "combine_sha256": {"two": "c" * 64, "sc3": "d" * 64, "bulk": "e" * 64}}


@pytest.fixture
def pins(monkeypatch):
    monkeypatch.setattr(k24, "K24_LISTING_SHA256", EREC["listing_sha256"])
    monkeypatch.setattr(k24, "K24_TABLE_SHA256", dict(EREC["tables"]))


def _pre(dropped=None, repeat_only=False) -> dict:
    s = _synth(7, n=200)
    pre = json.loads(json.dumps(k24.pregate_k24(s.phi, s.loaded, s.bases, s.den)))
    assert pre["dropped"] == {}
    pre["dropped"] = dict(dropped or {})
    pre["sc3_repeat_only"] = repeat_only
    return pre


def _write(root: Path, cells=None, gcells=None, commit="abc", clean=True, machine="m1", dropped=None, drop=None,
           guard=None, ext_rec=None, bad_file=None, k23_ext=None, repeat_only=False):
    cells = cells or {}
    pre = _pre(dropped, repeat_only)
    arms = k24.active_arms(pre)
    erec = ext_rec or EREC
    root.mkdir(parents=True)
    rec = {"dropped": k24.dropped_arms(pre), "arms": list(arms), "ext": erec, "constants": k24.constants_of(pre)}
    k18.write_task_file(root / "pregate.json", {"pregate": pre, "arms": list(arms), "ext": erec,
                                                "constants": k24.constants_of(pre),
                                                "provenance": _prov(commit, clean, machine)})
    refs, gsha = {}, {}
    for part, prefix, part_arms in (("guard", "guard/", k24.GUARDS), ("score", "", arms)):
        for w in k24.WORLDS:
            for n in k24.TASKS:
                a = (harness.FLOOR, *part_arms) if w == "plain" else part_arms
                src = (gcells or cells) if part == "guard" else cells
                rows = _rows({x: src.get((w, n), _table())[x] for x in a})
                refs[(w, n)] = [{"seed": r["seed"], "scores": {k24.BASE: r["scores"][k24.BASE]}} for r in rows]
                for r in rows:
                    s = r["seed"]
                    rel = prefix + (f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json")
                    f = root / rel
                    f.parent.mkdir(parents=True, exist_ok=True)
                    if drop == rel:
                        continue
                    rr = rec if bad_file != rel else {**rec, "constants": {**rec["constants"], "c": 0.5}}
                    body = ({"task": {"key": k24.TASKS[n].key}, "rows": [r], "seed": s, "pregate": rr,
                             "provenance": _prov(commit, clean, machine)} if w == "plain" else
                            {**r, "check": {k18.DM_CHECK: 60.0 + s}, "target": k18.DM_TARGETS[n], "task": n,
                             "pregate": rr, "arms": list(part_arms), "provenance": _prov(commit, clean, machine)})
                    k18.write_task_file(f, body)
                    if part == "guard":
                        gsha[rel] = _sha(f)
    g = guard or {"passes": True, "failed": []}
    k18.write_task_file(root / "guard.json", {**g, "guard1": {}, "guard2": {}, "files_sha256": gsha,
                                              "provenance": _prov(commit, clean, machine)})
    k18_dir = root.parent / f"k18_{root.name}"
    k18_dir.mkdir()
    for n in k24.TASKS:
        (k18_dir / f"plain_{n}.json").write_text(json.dumps({"rows": refs[("plain", n)]}))
    (k18_dir / "dm.json").write_text(json.dumps({"rows": {n: refs[("dm", n)] for n in k24.TASKS}}))
    k18pins = {f: _sha(k18_dir / f) for f in k24.K18_SHA256}
    k23_dir = root.parent / f"k23_{root.name}"
    k23_dir.mkdir()
    (k23_dir / "pregate.json").write_text("{}")
    (k23_dir / "reading.json").write_text(json.dumps({"ext": k23_ext or {"tables": {"all": TWO_SHA},
                                                                         "combine_sha256": "c" * 64}}))
    return SimpleNamespace(root=root, refs=k18_dir, k18pins=k18pins, k23=k23_dir)


@pytest.fixture
def u3d(tmp_path):
    rows = [{"target": t, "seed": s, "world": "depth_matched", "variant": "dm:transfer:mab21l2_all", "t3": 60.0 + s}
            for t in ("rep1", "rep2") for s in SEEDS]
    p = tmp_path / "orth_records.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def _read(w, u3d, monkeypatch):
    monkeypatch.setattr(k24, "K18_SHA256", w.k18pins)
    return k24.reading(w.root, w.refs, u3d, w.k23)


def test_the_reading_of_a_run_folder(tmp_path, u3d, monkeypatch, pins):
    w = _write(tmp_path / "ok", _everywhere(**GOOD))
    res = _read(w, u3d, monkeypatch)
    assert not res["void"], res.get("why")
    assert res["names"] == ["x_sc3"] and res["tries_counted"] == 9
    assert all(k18.guard_passes(d) for d in res["guard_match"].values()) and len(res["guard_match"]) == 8
    assert len(res["files_sha256"]) == 1 + 1 + 40 + 40
    txt = k24.text(res)
    for x in ("WORLD plain", "x_sc3 - x_e2", "x_sc3 - x_sc3_size", "guard part: passes True", "pre-gate",
              "k24sc3-B", "barred_form - floor"):
        assert x in txt, x


def test_a_dropped_table_reads_with_its_count(tmp_path, u3d, monkeypatch, pins):
    w = _write(tmp_path / "d", _everywhere(**GOOD), dropped={"bulk": "test"})
    res = _read(w, u3d, monkeypatch)
    assert not res["void"] and res["tries_counted"] == 3 and res["names"] == ["x_sc3"]
    assert res["arms"] == list(SC3_ONLY) and res["status"]["x_joint"]["ran"] is False
    w = _write(tmp_path / "j", _everywhere(**GOOD), dropped={"joint": "test"})
    assert _read(w, u3d, monkeypatch)["tries_counted"] == 6


@pytest.mark.parametrize("failed,tries", [(["barred_form"], 1), (["x_e2"], 1), (["barred_form", "x_e2"], 2)])
def test_a_failed_guard_part_is_void_and_counts_only_its_failed_guard_arms(tmp_path, u3d, monkeypatch, pins, failed,
                                                                          tries):
    w = _write(tmp_path / "g", _everywhere(**GOOD), guard={"passes": False, "failed": failed})
    for f in [p for p in w.root.rglob("*.json") if p.relative_to(w.root).parts[0] != "guard"
              and p.name not in ("pregate.json", "guard.json") and not p.name.endswith(".prov.json")]:
        f.unlink()  # no scoring part ran
    res = _read(w, u3d, monkeypatch)
    assert res["void"] and res["tries_counted"] == tries and res["names"] == [] and "guard part failed" in res["why"][0]
    assert "VOID" in k24.text(res)


def test_a_guard_mismatch_in_a_scoring_part_voids_with_no_added_count(tmp_path, u3d, monkeypatch, pins):
    g = _everywhere(**GOOD) | {("dm", "g1"): _table(**GOOD, barred_form=BASE_SCORE + 0.2)}
    w = _write(tmp_path / "m", _everywhere(**GOOD), gcells=g)
    res = _read(w, u3d, monkeypatch)
    assert res["void"] and res["tries_counted"] == 9 and any("guard match: dm g1: barred_form" in x for x in res["why"])


@pytest.mark.parametrize("kw,why", [
    ({"clean": False}, "dirty tree"), ({"bad_file": "dm/g1_seed_3.json"}, "constants differ"),
    ({"ext_rec": {**EREC, "listing_sha256": "0" * 64}}, "not the pinned ones"),
    ({"k23_ext": {"tables": {"all": TWO_SHA}, "combine_sha256": "0" * 64}}, "K23's committed reading.json"),
])
def test_void_readings(tmp_path, u3d, monkeypatch, pins, kw, why):
    w = _write(tmp_path / "v", _everywhere(**GOOD), **kw)
    res = _read(w, u3d, monkeypatch)
    assert res["void"] and any(why in x for x in res["why"]), res["why"]
    assert res["names"] == [] and res["tries_counted"] == 9


def test_missing_files_are_not_a_reading(tmp_path, u3d, monkeypatch, pins, capsys):
    w = _write(tmp_path / "i", drop="plain_g2/seed_7.json")
    res = _read(w, u3d, monkeypatch)
    assert res["incomplete"] and res["missing"] == ["plain_g2/seed_7.json"]
    w2 = _write(tmp_path / "i2", drop="guard/dm/g1_seed_2.json")
    assert _read(w2, u3d, monkeypatch)["missing"] == ["guard/dm/g1_seed_2.json"]
    (w2.root / "guard.json").unlink()
    assert _read(w2, u3d, monkeypatch)["missing"] == ["guard.json"]
    monkeypatch.setattr(k24, "reading", lambda runs, f=k24.reading: f(runs, w.refs, u3d, w.k23))
    with pytest.raises(SystemExit) as e:
        k24.main(["read", "--runs", str(w.root), "--out", str(tmp_path / "reading")])
    assert e.value.code == 3 and not (tmp_path / "reading").exists() and "not a reading" in capsys.readouterr().out


def test_read_writes_once_and_records_its_code(tmp_path, u3d, monkeypatch, pins, capsys):
    w = _write(tmp_path / "r")
    monkeypatch.setattr(k24, "K18_SHA256", w.k18pins)
    monkeypatch.setattr(k24, "reading", lambda runs, f=k24.reading: f(runs, w.refs, u3d, w.k23))
    out = tmp_path / "reading"
    k24.main(["read", "--runs", str(w.root), "--out", str(out)])
    body = json.loads((out / "reading.json").read_text())
    rc = body["reading_code"]
    assert rc["script_sha256"] == _sha(T3 / "k24_wider_tables.py") and rc["k23_sha256"] == _sha(
        T3 / "k23_outside_exponents.py")
    assert body["outcome"] == k24.NO_FILE and body["tries_counted"] == 9
    first = (out / "reading.json").read_bytes()
    with pytest.raises(SystemExit) as e:
        k24.main(["read", "--runs", str(w.root), "--out", str(out)])
    assert e.value.code == 4 and (out / "reading.json").read_bytes() == first
    capsys.readouterr()


# ------------------------------------------------------------------------------------------------ (9) the runs


def test_the_run_refuses_without_its_row(tmp_path, capsys):
    p = tmp_path / "TRIES.md"
    p.write_text("| K23 (one row) | T3 | `experiments/t3_decomp/k23_outside_exponents.py` |\n")
    with pytest.raises(SystemExit) as e:
        k24.row_present(p)
    assert e.value.code == 2 and "no row" in capsys.readouterr().err
    p.write_text(p.read_text() + "| K24 (one row) | T3 | `experiments/t3_decomp/k24_wider_tables.py` with |\n")
    assert k24.row_present(p) == "K24"
    for argv in (["run", "--world", "plain", "--task", "g1", "--out", str(tmp_path / "o")],
                 ["run", "--part", "guard", "--world", "dm", "--out", str(tmp_path / "o")],
                 ["pregate", "--out", str(tmp_path / "o")], ["guard", "--runs", str(tmp_path / "o")]):
        with pytest.MonkeyPatch.context() as m:
            m.setattr(k24, "TRIES_MD", tmp_path / "none.md")
            (tmp_path / "none.md").write_text("no rows\n")
            with pytest.raises(SystemExit) as e:
                k24.main(argv)
            assert e.value.code == 2
    assert not (tmp_path / "o").exists()


def test_preconditions_refuse_unfilled_pins_a_dirty_tree_and_a_bad_folder(tmp_path, monkeypatch, capsys):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| K24 (one row) | T3 | `experiments/t3_decomp/k24_wider_tables.py` |\n")
    monkeypatch.setattr(k24, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"})
    argv = ["pregate", "--out", str(tmp_path / "o")]
    # the unfilled state (the committed pins are filled since row K24's commit: the test sets it itself)
    with monkeypatch.context() as mp:
        mp.setattr(k24, "K24_LISTING_SHA256", None)
        mp.setattr(k24, "K24_TABLE_SHA256", {**k24.K24_TABLE_SHA256, "sc3": None, "bulk": None})
        with pytest.raises(SystemExit) as e:
            k24.main(argv)
        assert e.value.code == 2 and "not filled" in capsys.readouterr().err
    monkeypatch.setattr(k24, "K24_LISTING_SHA256", "a" * 64)
    monkeypatch.setattr(k24, "K24_TABLE_SHA256", {**k24.K24_TABLE_SHA256, "sc3": "b" * 64, "bulk": "d" * 64})
    monkeypatch.setattr(k24, "read_listing", lambda *a, **k: ({}, "x"))
    monkeypatch.delenv(k24.ENV_TWO, raising=False)
    monkeypatch.setenv(k24.ENV_TABLES, str(tmp_path))
    with pytest.raises(SystemExit) as e:
        k24.main(argv)
    assert e.value.code == 2 and "VEC_T3_EXT_TABLES_TWO is not set" in capsys.readouterr().err
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": False, "machine": "m1", "utc": "x"})
    with pytest.raises(SystemExit) as e:
        k24.main(argv)
    assert e.value.code == 2 and "not clean" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()


@pytest.fixture
def fake_run(monkeypatch, tmp_path):
    calls = {"seeds": [], "names": []}
    lines = {}
    s = _synth(9, n=200)
    pre = k24.pregate_k24(s.phi, s.loaded, s.bases, s.den)
    tabs = {k: v[0] for k, v in s.loaded["tables"].items()}
    basis = s.bases["plain"]["basis"]

    def run_task(t, cands, seeds, data_dir, log=print):
        consts = k24._K24["consts"]
        F = {"phi_barred": s.phi, "ext": tabs, "basis": basis}
        F["scaled"] = {k: k24.k24_table(k, tabs, consts) for k in ("sc3", "bulk", "joint")}
        F["scaled_perm"] = {k: ext.permuted(v, basis) for k, v in F["scaled"].items()}
        k18._FOLD["F"] = F
        calls["seeds"] += list(seeds)
        calls["names"] = [c.name for c in cands]
        for sd in seeds:
            lines[sd] = f"{t.key} seed {sd}: floor 50.00, barred_form {60 - sd / 9:.2f}, x_sc3 {61 + sd / 7:.2f}"
            log(lines[sd])
        return {"task": {"key": t.key}, "rows": [{"seed": sd} for sd in seeds],
                "provenance": {"commit": "abc", "tree_clean": True}}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=s.phi, genes=[f"g{i}" for i in range(200)], spec={}, record={})
        return {}

    record = dict(EREC)

    def load_k24(genes, phi_m):
        k24._K24.clear()
        k24._K24.update(genes=list(genes), T=dict(tabs), record=dict(record))
        return {}

    out = tmp_path / "run"
    out.mkdir()
    (out / "pregate.json").write_text(json.dumps({"pregate": pre, "arms": list(k24.active_arms(pre)),
                                                  "ext": record}))
    st = {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"}
    monkeypatch.setattr(k24, "preconditions", lambda: ("K24", st))
    monkeypatch.setattr(k24, "guard_refs", lambda **k: {})
    monkeypatch.setattr(k24, "load_k24", load_k24)
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(k21.k7, "_PREP", {"fold": {"record": {}}})
    monkeypatch.setattr(k18, "stamp", lambda: st)
    yield SimpleNamespace(calls=calls, lines=lines, out=out, pre=pre, st=st)
    k18._FOLD.clear()


def test_the_guard_part_runs_the_two_guard_arms_and_logs_progress_only(fake_run, capsys):
    out = fake_run.out
    k24.main(["run", "--part", "guard", "--world", "plain", "--task", "g1", "--seeds", "0-1", "--out", str(out)])
    assert fake_run.calls["names"] == [harness.FLOOR, *k24.GUARDS] and fake_run.calls["seeds"] == [0, 1]
    printed = capsys.readouterr()
    for x in printed.out.splitlines():
        assert re.fullmatch(r"K24 guard plain g1 seed \d: scored|wrote \S+/guard/plain_g1/seed_\d\.json: commit abc, "
                            r"clean True, machine m1", x), x
    for line in fake_run.lines.values():
        for score in re.findall(r"\d+\.\d\d", line):
            assert score not in printed.out and score not in printed.err
    body = json.loads((out / "guard" / "plain_g1" / "seed_0.json").read_text())
    assert body["part"] == "guard" and body["pregate"]["arms"] == list(k24.ARMS)
    assert set(body["world_figures"]) == {"x_e2"}


def test_the_scoring_parts_refuse_without_a_passing_guard_of_this_commit_and_machine(fake_run, capsys):
    out = fake_run.out
    argv = ["run", "--world", "plain", "--task", "g2", "--seeds", "0", "--out", str(out)]
    with pytest.raises(SystemExit) as e:
        k24.main(argv)
    assert e.value.code == 2 and "guard part runs and passes first" in capsys.readouterr().err
    for body, why in (({"passes": False, "provenance": fake_run.st}, "did not pass"),
                      ({"passes": True, "provenance": {**fake_run.st, "commit": "def"}}, "another commit"),
                      ({"passes": True, "provenance": {**fake_run.st, "machine": "m2"}}, "another commit")):
        (out / "guard.json").write_text(json.dumps(body))
        with pytest.raises(SystemExit) as e:
            k24.main(argv)
        assert e.value.code == 2 and why in capsys.readouterr().err
    assert fake_run.calls["seeds"] == []
    (out / "guard.json").write_text(json.dumps({"passes": True, "provenance": fake_run.st}))
    k24.main(argv)
    assert fake_run.calls["names"] == [harness.FLOOR, *k24.ARMS] and fake_run.calls["seeds"] == [0]
    for x in capsys.readouterr().out.splitlines():
        assert re.fullmatch(r"K24 plain g2 seed 0: scored|wrote \S+/plain_g2/seed_0\.json: commit abc, clean True, "
                            r"machine m1", x), x
    body = json.loads((out / "plain_g2" / "seed_0.json").read_text())
    assert body["pregate"]["constants"] == k24.constants_of(fake_run.pre)
    assert set(body["world_figures"]) == set(k24.ARMS[1:])


def test_a_crash_is_masked_and_a_world_off_the_pregate_writes_nothing(fake_run, monkeypatch, capsys):
    (fake_run.out / "guard.json").write_text(json.dumps({"passes": True, "provenance": fake_run.st}))
    body = json.loads((fake_run.out / "pregate.json").read_text())
    body["pregate"]["arms"]["x_bulk"]["controls"]["plain"]["sd_ln"] += 1e-6
    (fake_run.out / "pregate.json").write_text(json.dumps(body))
    with pytest.raises(SystemExit) as e:
        k24.main(["run", "--world", "plain", "--task", "g2", "--seeds", "0", "--out", str(fake_run.out)])
    assert e.value.code == 1 and "not the pre-gate's" in capsys.readouterr().err
    assert not (fake_run.out / "plain_g2" / "seed_0.json").exists()

    def boom(*a, **k):
        raise RuntimeError("score 61.234 went wrong")

    monkeypatch.setattr(harness, "run_task", boom)
    with pytest.raises(SystemExit) as e:
        k24.main(["run", "--part", "guard", "--world", "plain", "--task", "g1", "--seeds", "0", "--out",
                  str(fake_run.out)])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "61.234" not in err and "<num>" in err
    assert k18.ARMS == k24._K18_ORIGINAL["ARMS"]


def test_the_guard_check_writes_guard_json_and_exits_6_on_a_failure(tmp_path, monkeypatch, capsys, pins):
    w = _write(tmp_path / "gc", _everywhere(**GOOD))
    (w.root / "guard.json").unlink()
    (w.root / "guard.json.prov.json").unlink()
    st = _prov()
    monkeypatch.setattr(k24, "preconditions", lambda: ("K24", st))
    refs = _guard_refs()
    monkeypatch.setattr(k24, "guard_refs", lambda **k: refs)
    k24.main(["guard", "--runs", str(w.root)])
    g = json.loads((w.root / "guard.json").read_text())
    assert g["passes"] and g["failed"] == [] and set(g["files_sha256"]) == set(k24.guard_files())
    with pytest.raises(SystemExit) as e:
        k24.main(["guard", "--runs", str(w.root)])
    assert e.value.code == 2  # made once
    (w.root / "guard.json").unlink()
    refs["k23_cells"][("dm", "g1")]["severity_gated"] = 4
    with pytest.raises(SystemExit) as e:
        k24.main(["guard", "--runs", str(w.root)])
    assert e.value.code == 6 and json.loads((w.root / "guard.json").read_text())["failed"] == ["x_e2"]
    assert "FAILS" in capsys.readouterr().out
    (w.root / "guard.json").unlink()
    (w.root / "guard" / "dm" / "g2_seed_9.json").unlink()
    with pytest.raises(SystemExit) as e:
        k24.main(["guard", "--runs", str(w.root)])
    assert e.value.code == 3 and not (w.root / "guard.json").exists()


def test_the_pregate_command_writes_its_file_and_exits_5_when_both_tables_drop(tmp_path, monkeypatch, capsys):
    import vec.prep

    n = 200
    genes = [f"g{i}" for i in range(n)]
    C = np.random.default_rng(0).poisson(3.0, size=(30, n)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    wt = ad.AnnData(X=np.log1p(C).astype(np.float32), var=pd.DataFrame(index=genes), layers={"counts": C})
    s = _synth(5, n=n)
    seen = {"loaded": s.loaded}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=s.phi, genes=genes, spec={}, record={"file": "f"})
        return {"file": "f"}

    def load_k24(g, phi_m):
        k24._K24.clear()
        k24._K24.update(genes=list(g), T={k: v[0] for k, v in seen["loaded"]["tables"].items()},
                        record=dict(EREC), table_records={})
        return seen["loaded"]

    monkeypatch.setattr(k24, "preconditions", lambda: ("K24", _prov()))
    monkeypatch.setattr(k24, "guard_refs", lambda **k: {})
    monkeypatch.setattr(k24, "k23_committed", lambda *a: ({"pregate": None}, {}))
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(vec.prep, "load_stage", lambda p, drop_duplicates=False, genes=None: wt)
    monkeypatch.setattr(k24, "load_k24", load_k24)
    monkeypatch.setattr(k24, "world_bases", lambda data, g, bs: {
        "plain": {"genes": list(g), "basis": np.arange(n, dtype=float)},
        "dm": {"genes": list(g)[::-1], "basis": np.arange(n, dtype=float)}})
    k24.main(["pregate", "--out", str(tmp_path / "a")])
    body = json.loads((tmp_path / "a" / "pregate.json").read_text())
    assert body["row"] == "K24" and body["ext"] == EREC and body["arms"] == list(k24.active_arms(body["pregate"]))
    assert body["tries_if_run"] == k24.try_count(body["arms"]) and body["constants"] == k24.constants_of(
        body["pregate"])
    assert "pre-gate (the tables, phi_M and the worlds' bases alone; no score)" in capsys.readouterr().out
    with pytest.raises(SystemExit) as e:
        k24.main(["pregate", "--out", str(tmp_path / "a")])
    assert e.value.code == 2
    few = np.ones(n)
    few[:3] = 1.1
    t = seen["loaded"]["tables"]
    seen["loaded"] = {**s.loaded, "tables": {**t, "sc3": (few, t["sc3"][1]), "bulk": (few.copy(), t["bulk"][1])}}
    with pytest.raises(SystemExit) as e:
        k24.main(["pregate", "--out", str(tmp_path / "b")])
    assert e.value.code == 5 and (tmp_path / "b" / "pregate.json").exists()
    assert "nothing runs" in capsys.readouterr().out


# ------------------------------------------------------------------------------------- (9b) the two restated clauses
# The cold challenger's restatement of 2026-11-20 (before any table value, agreement count or score was read): (b4) at
# n_g lines per gene with an excess of 0.11 over chance, and the repeat clause that drops x_sc3's arms only. Synthetic
# tables only. Numbers in the comments refer to the challenger's test list.


def _b4(n, U, shares, held=..., T=None, **over):
    """A one-study bulk table of n genes (none exactly 1 unless T says so) whose manifest records `held` as
    lines_held_by_study[MGP-DMDD] (an int for all genes, a list, or ... for no field)."""
    T = np.exp(np.linspace(-0.05, 0.05, n) + 1e-3) if T is None else T
    rec = _bulk_rec(n, shares, lines=U, libs=(4, 4), **over)
    if held is not ...:
        rec["lines_held_by_study"] = {"MGP-DMDD": [held] * n if isinstance(held, int) else list(held)}
    return k24.bulk_drops(T, rec)


def _b4_fails(why):
    return any(w.startswith("(b4)") for w in why)


def test_1_chance_agreement_is_the_exact_fraction():
    from fractions import Fraction

    want = {12: Fraction(397, 1024), 13: Fraction(1093, 4096), 14: Fraction(1471, 8192), 15: Fraction(309, 1024),
            16: Fraction(6885, 32768), 25: Fraction(1807781, 16777216)}
    for n, f in want.items():
        got = k24.chance_agreement(n)
        assert isinstance(got, Fraction) and got == f, n
        assert got == 2 * Fraction(sum(math.comb(n, j) for j in range(math.ceil(2 * n / 3), n + 1)), 2 ** n)


@pytest.mark.parametrize("U", [12, 16])
def test_2_a_count_of_ceil_need_passes_and_one_less_fails(U):
    from fractions import Fraction

    n = 300
    need = max(Fraction(50), n * k24.chance_agreement(U) + Fraction(11, 100) * n)
    k = math.ceil(need)
    assert need > 50 and k < n
    why, figs = _b4(n, U, [0.75] * k + [0.5] * (n - k), held=U)
    assert not _b4_fails(why), why
    assert figs["agree_G"] == n and figs["agree_count"] == k and figs["agree_needed"] == float(need)
    assert figs["E0"] == float(n * k24.chance_agreement(U)) and figs["n_g_source"].startswith("manifest")
    assert figs["n_g_histogram"] == {str(U): n} and figs["agree_left_out_few_lines"] == 0
    why, figs = _b4(n, U, [0.75] * (k - 1) + [0.5] * (n - k + 1), held=U)
    assert _b4_fails(why) and figs["agree_count"] == k - 1
    assert not any(w.startswith(("(b1)", "(b2)", "(b3)", "(b5)")) for w in why), why


def test_3_genes_held_in_fewer_than_12_lines_are_left_out_of_the_count_and_g():
    n, U = 300, 16
    held = [11] * 100 + [16] * 200
    shares = [0.9] * 100 + [0.5] * 200
    why, figs = _b4(n, U, shares, held=held)
    assert figs["agree_G"] == 200 and figs["agree_count"] == 0 and figs["agree_left_out_few_lines"] == 100
    assert figs["E0"] == float(200 * k24.chance_agreement(16)) and figs["n_g_histogram"] == {"11": 100, "16": 200}
    assert _b4_fails(why)
    # all held in 12 lines or more: the same agreeing genes now count
    _, f2 = _b4(n, U, shares, held=[12] * 100 + [16] * 200)
    assert f2["agree_G"] == 300 and f2["agree_count"] == 100


def test_4_the_fallback_n_g_u_is_used_and_recorded_and_a_malformed_field_drops_without_raising():
    n, U = 200, 14
    shares = [0.8] * n
    _, f = _b4(n, U, shares)  # no lines_held_by_study in the manifest
    assert f["n_g_source"] == "fallback U (the manifest records no lines_held_by_study)"
    assert f["n_g_histogram"] == {"14": n} and f["E0"] == float(n * k24.chance_agreement(14))
    bad = {"not a mapping": [U] * n, "one per gene": {"MGP-DMDD": [U] * (n - 1)},
           "none for MGP-DMDD": {"OTHER": [3] * n},   # the fourth review's point 1: a contradiction, not a fallback
           "not an integer": {"MGP-DMDD": [float(U)] * n}, "a bool": {"MGP-DMDD": [True] * n},
           "outside": {"MGP-DMDD": [U + 1] * n}, "negative": {"MGP-DMDD": [-1] + [U] * (n - 1)},
           "0 lines": {"MGP-DMDD": [0] + [U] * (n - 1)}}
    for case, field in bad.items():
        why, f = _b4(n, U, shares, lines_held_by_study=field)
        assert _b4_fails(why) and any("malformed" in w for w in why), case
        assert f["n_g_source"] == "malformed" and f["agree_G"] is None and f["agree_needed"] is None, case
    # 0 lines where the gene has no agreement value is well formed
    _, f = _b4(n, U, [None] + shares[1:], held=[0] + [U] * (n - 1))
    assert f["n_g_source"].startswith("manifest") and f["agree_G"] == n - 1
    # in the pre-gate a malformed field drops bulk (and the joint) and keeps sc3; nothing raises
    s0 = _synth()
    rec = _bulk_rec(s0.n, s0.shares, lines_held_by_study={"MGP-DMDD": [99] * s0.n})
    pre = _drops(bulk_rec=rec)
    assert "malformed" in pre["dropped"]["bulk"] and "joint" in pre["dropped"] and "sc3" not in pre["dropped"]
    assert pre["bulk"]["n_g_source"] == "malformed"


def test_5_a_share_of_exactly_two_thirds_counts():
    n, U = 120, 15
    held = [12] * 60 + [15] * 60
    shares = [8 / 12] * 60 + [10 / 15] * 60
    _, f = _b4(n, U, shares, held=held)
    assert f["agree_count"] == 120 == f["agree_G"]
    _, f = _b4(n, U, [7 / 12] * 60 + [2 / 3 - 1e-9] * 60, held=held)
    assert f["agree_count"] == 0


def test_6_genes_equal_to_1_and_none_shares_are_skipped():
    n, U = 200, 16
    T = np.exp(np.linspace(-0.05, 0.05, n) + 1e-3)
    T[:30] = 1.0
    shares = [0.9] * 30 + [None] * 20 + [0.9] * 150
    _, f = _b4(n, U, shares, held=U, T=T)
    assert f["agree_G"] == 150 and f["agree_count"] == 150 and f["n_g_histogram"] == {"16": 150}


def test_7_the_floor_of_50_binds_at_g_100_and_u_25():
    from fractions import Fraction

    n, U = 100, 25
    assert 100 * k24.chance_agreement(25) + Fraction(11) < 50
    why, f = _b4(n, U, [0.8] * 50 + [0.5] * 50, held=U)
    assert f["agree_needed"] == 50.0 and not _b4_fails(why)
    why, f = _b4(n, U, [0.8] * 49 + [0.5] * 51, held=U)
    assert _b4_fails(why) and f["agree_count"] == 49


def _sc3_like_two(corr: str) -> np.ndarray:
    """A synthetic sc3 in _synth()'s world, built from `two`, `bulk` and noise: corr(ln sc3, ln two) 0.953 ('0.95') or
    0.861 ('0.85'), corr(ln sc3, ln bulk) 0.26 or 0.41 (inside the joint's window)."""
    s0 = _synth()
    lt, lb = np.log(s0.two), np.log(s0.bulk)
    z = np.random.default_rng(99).normal(size=s0.n)
    b, g = {"0.95": (0.3, 0.1), "0.85": (0.5, 0.3)}[corr]
    ls = lt / lt.std() + b * lb / lb.std() + g * z
    return np.exp(0.03 * ls / ls.std())


REPEAT_ARMS = ("barred_form", "x_e2", "x_bulk", "x_bulk_perm", "x_bulk_size", "x_joint", "x_joint_perm", "x_joint_size")


def test_8_a_repeat_of_the_named_form_with_a_valid_bulk_runs_bulk_and_the_joint():
    pre = _drops(sc3=_sc3_like_two("0.95"))
    assert 0.95 <= pre["corr"]["sc3_two"] < 0.96 and 0 <= pre["corr"]["sc3_bulk"] < 0.9
    assert pre["sc3_repeat_only"] is True and "repeat of the named form" in pre["dropped"]["sc3"]
    assert set(pre["dropped"]) == {"sc3"}
    arms = k24.active_arms(pre)
    assert arms == REPEAT_ARMS and k24.try_count(arms) == 6
    assert not any(a.startswith("x_sc3") for a in arms) and not k24.nothing_runs(arms)
    assert set(k24.dropped_arms(pre)) == {"x_sc3", "x_sc3_perm", "x_sc3_size"}
    co = k24.constants_of(pre)
    assert co["e_sc3"] == pre["e"]["sc3"] is not None and co["c"] == pre["c"] is not None
    assert set(co["twins"]) == {"x_bulk", "x_joint"}
    # the joint is built from sc3 with the same formula, at S*
    assert math.isclose(pre["scaled"]["joint"]["sd_ln"], k24.S_STAR, rel_tol=1e-12)
    assert "x_joint" in pre["arms"] and "sc3 enters the joint unchanged" in "\n".join(k24.pregate_lines(pre))
    # fold_factors' constants: the joint's scaled table is made, and x_sc3's arms are not run
    F = {"phi_barred": np.ones(3), "ext": {"sc3": np.array([1.1, 0.9, 1.0]), "bulk": np.array([0.95, 1.05, 1.0]),
                                         "two": np.ones(3)}}
    F["scaled"] = {k: k24.k24_table(k, F["ext"], co) for k in ("bulk", "joint")}
    k24._K24["consts"] = co
    want = 1.1 ** (co["c"] * co["e_sc3"]) * 0.95 ** (co["c"] * co["e_bulk"])
    assert math.isclose(k24.arm_factor("x_joint", F)[0], want, rel_tol=1e-14)


def test_9_a_repeat_with_bulk_dropped_runs_nothing():
    pre = _drops(sc3=_sc3_like_two("0.95"), bulk_rec=_bulk_bad("b2"))
    assert pre["sc3_repeat_only"] is True and "bulk is dropped" in pre["dropped"]["joint"]
    assert "sc3 is dropped" not in pre["dropped"]["joint"]
    arms = k24.active_arms(pre)
    assert arms == ("barred_form", "x_e2") and k24.nothing_runs(arms) and k24.try_count(arms) == 0


def test_9_the_pregate_command_exits_5_on_a_repeat_with_bulk_dropped(tmp_path, monkeypatch, capsys):
    import vec.prep

    s = _synth()
    n = s.n
    genes = [f"g{i}" for i in range(n)]
    C = np.random.default_rng(0).poisson(3.0, size=(30, n)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    wt = ad.AnnData(X=np.log1p(C).astype(np.float32), var=pd.DataFrame(index=genes), layers={"counts": C})
    t = s.loaded["tables"]
    loaded = {**s.loaded, "tables": {**t, "sc3": (_sc3_like_two("0.95"), t["sc3"][1]),
                                     "bulk": (t["bulk"][0], _bulk_bad("b2"))}}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=s.phi, genes=genes, spec={}, record={"file": "f"})
        return {"file": "f"}

    def load_k24(g, phi_m):
        k24._K24.clear()
        k24._K24.update(genes=list(g), T={k: v[0] for k, v in loaded["tables"].items()}, record=dict(EREC),
                        table_records={})
        return loaded

    monkeypatch.setattr(k24, "preconditions", lambda: ("K24", _prov()))
    monkeypatch.setattr(k24, "guard_refs", lambda **k: {})
    monkeypatch.setattr(k24, "k23_committed", lambda *a: ({"pregate": None}, {}))
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(vec.prep, "load_stage", lambda p, drop_duplicates=False, genes=None: wt)
    monkeypatch.setattr(k24, "load_k24", load_k24)
    monkeypatch.setattr(k24, "world_bases", lambda data, g, bs: {
        "plain": {"genes": list(g), "basis": np.arange(n, dtype=float)},
        "dm": {"genes": list(g)[::-1], "basis": np.arange(n, dtype=float)}})
    with pytest.raises(SystemExit) as e:
        k24.main(["pregate", "--out", str(tmp_path / "r")])
    assert e.value.code == 5
    body = json.loads((tmp_path / "r" / "pregate.json").read_text())
    assert body["pregate"]["sc3_repeat_only"] is True and body["tries_if_run"] == 0 and body["arms"] == [
        "barred_form", "x_e2"]
    assert body["constants"]["e_sc3"] == body["pregate"]["e"]["sc3"] is not None
    assert "nothing runs" in capsys.readouterr().out


def test_10_a_repeat_with_another_sc3_failure_drops_the_joint(monkeypatch):
    sc3 = _sc3_like_two("0.95")
    tiny = np.where(sc3 != 1.0, sc3 ** 0.05, 1.0)  # e_sc3 above 8, the same corr with `two`
    pre = _drops(sc3=tiny)
    assert pre["corr"]["sc3_two"] >= 0.9 and "(s1)" in pre["dropped"]["sc3"]
    assert pre["sc3_repeat_only"] is False and "sc3 is dropped" in pre["dropped"]["joint"]
    assert k24.active_arms(pre) == BULK_ONLY and k24.try_count(k24.active_arms(pre)) == 3
    # a failed product bound of sc3 at e_sc3 (the first product the pre-gate makes is x_sc3's)
    original, calls = k24.arm_figures, []

    def bound(*a, **k):
        figs, why = original(*a, **k)
        calls.append(1)
        return figs, why + (["the product fails a bound (test)"] if len(calls) == 1 else [])

    monkeypatch.setattr(k24, "arm_figures", bound)
    pre = _drops(sc3=sc3)
    assert "fails a bound (test)" in pre["dropped"]["sc3"] and "repeat of the named form" in pre["dropped"]["sc3"]
    assert pre["sc3_repeat_only"] is False and "sc3 is dropped" in pre["dropped"]["joint"]
    assert k24.active_arms(pre) == BULK_ONLY


def test_11_bytes_equal_to_two_s_or_an_undefined_corr_with_two_drop_the_joint(monkeypatch):
    s0 = _synth()
    rec = {**s0.loaded["tables"]["sc3"][1], "sha256": "2" * 64}
    pre = _drops(sc3=_sc3_like_two("0.95"), sc3_rec=rec)
    assert "`two`'s" in pre["dropped"]["sc3"] and "repeat" in pre["dropped"]["sc3"]
    assert pre["sc3_repeat_only"] is False and "sc3 is dropped" in pre["dropped"]["joint"]
    assert k24.active_arms(pre) == BULK_ONLY
    original = ext.corr_between
    monkeypatch.setattr(ext, "corr_between", lambda a, b: None if b is s0.two else original(a, b))
    pre = k24.pregate_k24(s0.phi, s0.loaded, s0.bases, s0.den)
    assert "corr(ln sc3, ln two) undefined" in pre["dropped"]["sc3"] and "repeat" not in pre["dropped"]["sc3"]
    assert pre["sc3_repeat_only"] is False and "sc3 is dropped" in pre["dropped"]["joint"]
    assert k24.active_arms(pre) == BULK_ONLY


def test_12_a_corr_with_two_below_0_9_gives_nine_arms():
    pre = _drops(sc3=_sc3_like_two("0.85"))
    assert 0.85 <= pre["corr"]["sc3_two"] < 0.9 and pre["dropped"] == {} and pre["sc3_repeat_only"] is False
    assert k24.active_arms(pre) == k24.ARMS and k24.try_count(k24.ARMS) == 9


def test_the_try_counts_of_the_restated_clause():
    """The eight cases of the restatement, from the pre-gate's drop record."""
    def count(dropped, repeat=False):
        return k24.try_count(k24.active_arms({"dropped": dropped, "sc3_repeat_only": repeat}))

    assert count({}) == 9
    assert count({"sc3": "repeat"}, True) == 6
    assert count({"sc3": "repeat", "joint": "its own window"}, True) == 3
    assert count({"sc3": "(s1)", "joint": "sc3 is dropped"}) == 3
    assert count({"bulk": "(b2)", "joint": "bulk is dropped"}) == 3
    assert count({"bulk": "(b2)", "sc3": "repeat", "joint": "bulk is dropped"}, True) == 0
    assert count({"bulk": "(b2)", "sc3": "(s1)", "joint": "x"}) == 0
    assert count({"joint": "its own window"}) == 6
    # the joint is dropped with a table even when the pre-gate's record of the joint is missing
    assert count({"sc3": "(s1)"}) == 3 and count({"sc3": "repeat"}, True) == 6


def test_a_repeat_reads_with_six_tries_and_may_name_the_joint(tmp_path, u3d, monkeypatch, pins):
    w = _write(tmp_path / "rep", _everywhere(x_joint=57.5, x_e2=57.0),
               dropped={"sc3": "corr(ln sc3, ln two) 0.950 (>= 0.9: a repeat of the named form)"}, repeat_only=True)
    res = _read(w, u3d, monkeypatch)
    assert not res["void"], res.get("why")
    assert res["arms"] == list(REPEAT_ARMS) and res["tries_counted"] == 6 and res["names"] == ["x_joint"]
    assert res["status"]["x_sc3"]["ran"] is False and res["files"] == [
        "k24joint-B (base B of decision log 2026-11-20 (q))"]
    assert "x_sc3: not run (dropped by the pre-gate)" in k24.text(res)


def test_14_the_tie_order_and_f_are_unchanged_in_the_repeat_state():
    arms = REPEAT_ARMS
    res = k24.decide(_per(_everywhere(x_bulk=57.5, x_joint=57.5, x_e2=57.0), arms=arms), arms)
    assert res["qualified"] == ["x_bulk", "x_joint"] and res["names"] == ["x_bulk"]
    res = k24.decide(_per(_everywhere(x_joint=57.5, x_e2=57.3), arms=arms), arms)
    s = res["status"]["x_joint"]
    assert s["qualified_a_to_e"] and not s["f"] and res["names"] == [] and res["reported"] == {
        "x_joint": k24.NOT_AHEAD}
    assert k24.NAMEABLE == ("x_sc3", "x_bulk", "x_joint") and k24.MIN_GAIN == 0.3 and k24.MIN_AHEAD == 8


# ------------------------------------------------------------------------------------------------ (10) the launcher


def test_the_launcher_constants_and_its_parts_are_k23s():
    sh = (T3 / "run_k24.sh").read_text()
    k23sh = (T3 / "run_k23.sh").read_text()
    assert "SCRIPT=experiments/t3_decomp/k24_wider_tables.py" in sh
    assert "S=~/v/scratch/k24_${COMMIT:0:7}" in sh and "R=runs/k24_${COMMIT:0:7}" in sh
    assert "OUT=experiments/t3_decomp/results/k24_wider_tables" in sh
    assert "COMMIT=${1:?usage: run_k24.sh COMMIT}" in sh
    assert not re.search(r"\$\{[^}]*:\?[^}]*'[^}]*\}", sh)  # no apostrophe inside a ${var:?message}
    assert "TABLES=$HOME/v/vec-ext/t3/tables  # combine.py's --out (sc3, bulk); fixed, no override" in sh
    assert "TABLES_TWO=$HOME/v/vec-ext/t3/tables_two_study" in sh and "VEC_T3_EXT_TABLES_TWO=\"$TABLES_TWO\"" in sh
    for f in ("experiments/t3_ext/TABLES_K24.sha256", "experiments/t3_decomp/results/k23_outside_exponents/"
              "pregate.json", "experiments/t3_decomp/results/k23_outside_exponents/reading.json"):
        assert f in sh
    assert "for acc in E-MTAB-7322 GSE299068; do" in sh
    for t in ("tests/test_k24_wider_tables.py", "tests/test_build_k24_form.py", "tests/test_k23_outside_exponents.py",
              "tests/test_k22_outside_forms.py", "tests/test_t3_outside_forms.py", "tests/test_k21_factor_forms.py",
              "tests/test_k18_built_files.py"):
        assert t in sh, t
    for code in (10, 20, 22, 25, 26, 30, 41, 43, 44, 95, 96, 97, 98):
        assert f"finish {code}" in sh, code
    assert "if [ $rc = 6 ]; then" in sh and "finish 26" in sh and "run --part guard" in sh
    assert sh.index("run --part guard") < sh.index("PARTS=(") and sh.index("guard --runs") < sh.index("PARTS=(")
    # the scoring parts and the read are K23's, with K23 replaced by K24
    part = lambda s: s.split("PARTS=", 1)[1].split("rc=0\n$PY \"$SCRIPT\" read", 1)[0]  # noqa: E731
    assert part(sh) == part(k23sh).replace("K23", "K24")
    tail = lambda s: s.split("if [ $bad != 0 ]; then finish 30; fi\n", 1)[1]  # noqa: E731
    assert tail(sh) == tail(k23sh)


def _launcher_repo(tmp_path: Path) -> SimpleNamespace:
    repo = tmp_path / "repo"
    home = tmp_path / "home"
    for d in ("v/vec-ext/t3/tables", "v/vec-ext/t3/tables_two_study"):
        (home / d).mkdir(parents=True)
    files = {
        "experiments/t3_decomp/run_k24.sh": (T3 / "run_k24.sh").read_text(),
        "experiments/harness/TRIES.md": "| K24 (one row) | T3 | `experiments/t3_decomp/k24_wider_tables.py` |\n",
        "EXTERNAL_DATA.md": "| X015 | GEO GSE299068 ... |\n| X017 | ArrayExpress E-MTAB-7322 ... |\n",
        "experiments/t3_ext/TABLES_K24.sha256": "x\n",
        "experiments/t3_decomp/results/k23_outside_exponents/pregate.json": "{}\n",
        "experiments/t3_decomp/results/k23_outside_exponents/reading.json": "{}\n",
    }
    repo.mkdir()
    _git(repo, "init", "-q")
    for rel, t in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t)
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "init")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "uv").write_text("#!/bin/sh\nexit 1\n")
    (bindir / "uv").chmod(0o755)
    return SimpleNamespace(repo=repo, home=home, bin=bindir, head=_git(repo, "rev-parse", "HEAD"))


def _launch(L, *args) -> tuple[int, str]:
    env = {**_git_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}",
           "GIT_CEILING_DIRECTORIES": str(L.repo.parent)}
    r = subprocess.run(["bash", "experiments/t3_decomp/run_k24.sh", *args], cwd=L.repo, env=env, capture_output=True,
                       text=True, timeout=60)
    return r.returncode, r.stderr


def _all_exit(L, head) -> str:
    return (L.home / "v" / "scratch" / f"k24_{head[:7]}" / "all.exit").read_text().strip()


def test_the_launcher_s_refusal_paths_run_in_a_temporary_repository(tmp_path):
    L = _launcher_repo(tmp_path)
    code, err = _launch(L)
    assert code != 0 and "usage: run_k24.sh COMMIT" in err
    code, err = _launch(L, "0" * 40)
    assert code == 97 and "HEAD is not" in err and _all_exit(L, "0" * 40) == "97"
    # past every refusal, the launcher reaches the tests (a stub uv that fails: exit 10)
    code, err = _launch(L, L.head)
    assert code == 10 and _all_exit(L, L.head) == "10", err
    # a dirty tracked file
    (L.repo / "EXTERNAL_DATA.md").write_text("changed\n")
    code, err = _launch(L, L.head)
    assert code == 98 and "not clean" in err
    _git(L.repo, "checkout", "--", "EXTERNAL_DATA.md")
    # no row naming the script
    (L.repo / "experiments/harness/TRIES.md").write_text("| K23 | `experiments/t3_decomp/k23_outside_exponents.py` |\n")
    _git(L.repo, "commit", "-q", "-am", "no row")
    head = _git(L.repo, "rev-parse", "HEAD")
    code, err = _launch(L, head)
    assert code == 95 and "no TRIES.md row" in err
    _git(L.repo, "revert", "--no-edit", "HEAD")
    # the listing untracked
    _git(L.repo, "rm", "-q", "--cached", "experiments/t3_ext/TABLES_K24.sha256")
    _git(L.repo, "commit", "-q", "-m", "untrack")
    head = _git(L.repo, "rev-parse", "HEAD")
    code, err = _launch(L, head)
    assert code == 96 and "TABLES_K24.sha256 is not tracked" in err
    _git(L.repo, "add", "experiments/t3_ext/TABLES_K24.sha256")
    _git(L.repo, "commit", "-q", "-m", "track")
    head = _git(L.repo, "rev-parse", "HEAD")
    # a missing folder, then a symlinked one
    two = L.home / "v/vec-ext/t3/tables_two_study"
    shutil.rmtree(two)
    code, err = _launch(L, head)
    assert code == 96 and "tables_two_study is not a real directory" in err
    (L.home / "elsewhere").mkdir()
    two.symlink_to(L.home / "elsewhere")
    code, err = _launch(L, head)
    assert code == 96 and "not a real directory" in err
    two.unlink()
    two.mkdir()
    # a register without either accession's row
    for acc in ("E-MTAB-7322", "GSE299068"):
        text_ = (L.repo / "EXTERNAL_DATA.md").read_text()
        (L.repo / "EXTERNAL_DATA.md").write_text("".join(ln for ln in text_.splitlines(True) if acc not in ln))
        _git(L.repo, "commit", "-q", "-am", f"drop {acc}")
        h = _git(L.repo, "rev-parse", "HEAD")
        code, err = _launch(L, h)
        assert code == 96 and f"no row naming {acc}" in err
        _git(L.repo, "revert", "--no-edit", "HEAD")
    head = _git(L.repo, "rev-parse", "HEAD")
    code, _ = _launch(L, head)
    assert code == 10


UV_STUB = """#!/bin/sh
# uv run --frozen --no-sync python -I -B ARGS: the tests pass, a -c one-liner runs in python -I -B, and every call of
# the script is logged ($HOME/script_calls) and exits 1 (read: 0)
shift 6
case "$1" in
  -m) exit 0 ;;
  -c) exec "{py}" -I -B "$@" ;;
esac
echo "$*" >> "$HOME/script_calls"
[ "$2" = read ] && exit 0
exit 1
"""


def test_a_rerun_after_a_guard_failure_finishes_26_again_and_starts_no_scoring_part(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.bin / "uv").write_text(UV_STUB.replace("{py}", sys.executable))
    R = L.repo / "runs" / f"k24_{L.head[:7]}"
    R.mkdir(parents=True)
    (R / "pregate.json").write_text("{}")
    calls = L.home / "script_calls"
    script = "experiments/t3_decomp/k24_wider_tables.py"
    out = "experiments/t3_decomp/results/k24_wider_tables"
    called = lambda: calls.read_text().splitlines() if calls.exists() else []  # noqa: E731
    # guard.json says the guard part failed and OUT holds no reading: the VOID read once, then 26
    (R / "guard.json").write_text(json.dumps({"passes": False, "failed": ["x_e2"]}))
    code, err = _launch(L, L.head)
    assert code == 26 and _all_exit(L, L.head) == "26", err
    assert called() == [f"{script} read --runs runs/k24_{L.head[:7]} --out {out}"]
    # the reading is there (as the first finish 26 left it): 26 again, no call of the script at all
    (L.repo / out).mkdir(parents=True)
    (L.repo / out / "reading.json").write_text("{}")
    calls.unlink()
    for body in ({"passes": False, "failed": ["x_e2"]}, {"passes": "true"}, {"failed": []}):
        (R / "guard.json").write_text(json.dumps(body))
        code, err = _launch(L, L.head)
        assert code == 26 and _all_exit(L, L.head) == "26", (body, err)
        assert called() == [] and not (R / "exit" / "plain_g1_0-4").exists()
    # an unreadable guard.json finishes 20, no scoring part
    for raw in ("not json", "[]"):
        (R / "guard.json").write_text(raw)
        code, err = _launch(L, L.head)
        assert code == 20 and _all_exit(L, L.head) == "20", (raw, err)
        assert called() == []
    # a passing guard.json goes on to the scoring parts (the stub fails them: 30)
    (R / "guard.json").write_text(json.dumps({"passes": True, "failed": []}))
    code, err = _launch(L, L.head)
    assert code == 30 and _all_exit(L, L.head) == "30", err
    assert sorted(c.split()[1] for c in called()) == ["run"] * 5
