"""Tests for experiments/t3_decomp/k23_outside_exponents.py (try K23: the barred factor times K22's outside-knockout
table at exponents 2 and 4, each with its permuted control and its size twin). Offline, synthetic and fast, in the
manner of tests/test_k22_outside_forms.py: K18's plain fold fixture with stubbed construct, composition and assemble,
hand-made tables, synthetic task files. No competition data and no outside data. K22's committed result files are
read for their pins and code hashes only."""

from __future__ import annotations

import hashlib
import json
import math
import re
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
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402
import k22_outside_forms as k22  # noqa: E402
import k23_outside_exponents as k23  # noqa: E402
import outside_forms as ext  # noqa: E402

REAL_LISTING, REAL_TABLES = ext.K22_LISTING_SHA256, dict(ext.K22_TABLE_SHA256)  # the committed pins, before patching
K22_DIR = T3 / "results" / "k22_outside_forms"
G, N_POOL, N_DRAW = 14, 90, 150
GENES = [f"g{i}" for i in range(G)]
FILE_GENES = GENES[::-1] + ["extra"]
PINS = {"listing_sha256": k23.K23_LISTING_SHA256, "tables": {"all": k23.K23_TABLE_SHA256}}


@pytest.fixture(autouse=True)
def clean_state():
    yield
    k22._EXT.clear()
    k23._K23.clear()


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ (1) constants, pins


def test_arms_constants_and_the_shared_modules():
    assert k23.ARMS == ("barred_form", "x_e1", "x_e2", "x_e2_perm", "x_e2_size", "x_e4", "x_e4_perm", "x_e4_size")
    assert k23.NAMEABLE == ("x_e2", "x_e4") and k23.GUARDS == ("barred_form", "x_e1") and k23.EXPONENTS == (2.0, 4.0)
    assert k23.CONTROL_OF == {"x_e2": "x_e2_perm", "x_e4": "x_e4_perm"}
    assert k23.SIZE_OF == {"x_e2": "x_e2_size", "x_e4": "x_e4_size"}
    assert {a: k23.EXP_OF[a] for a in k23.ARMS if a != k23.BASE} == {
        "x_e1": 1.0, "x_e2": 2.0, "x_e2_perm": 2.0, "x_e2_size": 2.0, "x_e4": 4.0, "x_e4_perm": 4.0, "x_e4_size": 4.0}
    assert k23.TABLE == "all" and k23.BASE == k18.BARRED and k23.STREAM == k18.STREAM == "t3_s2_form"
    assert (k23.MAX_SD_LN_TABLE, k23.MAX_AT_CLIP, k23.MAX_DEN_SHIFT, k23.MAX_SD_RATIO) == (0.25, 45, 0.05, 1.05)
    assert k23.DOSE_MIN == 0.0 and k23.PIVOT_BEHIND == -k22.MIN_GAIN == -0.3
    assert (k23.MIN_GAIN, k23.MIN_AHEAD, k23.MAX_EXTRA_GATED) == (k22.MIN_GAIN, k22.MIN_AHEAD, k22.MAX_EXTRA_GATED)
    assert (k23.MIN_GAIN, k23.MIN_AHEAD, k23.MAX_EXTRA_GATED) == (0.3, 8, 1)
    assert k23.X_E1_TOL_MEAN == k18.GUARD_MEAN == k23.GUARD_MEAN == 0.1 and k23.GUARD_SEED == 0.5
    assert k23.ext is ext and k23.k22 is k22 and k23.TASKS == k18.TASKS and k23.WORLDS == k18.WORLDS
    assert k23.K18_SHA256 == k22.K18_SHA256 == k21.K18_SHA256
    assert k23.RESULTS == T3 / "results" / "k23_outside_exponents" and k23.K22_RESULTS == K22_DIR
    assert ext.PERM_SEED == 20261113 and ext.PERM_BINS == 10 and ext.EXPONENT == 1.0
    assert k23.TRIES_MD == REPO / "experiments" / "harness" / "TRIES.md"


def test_the_pins_are_k22s_and_the_shared_code_is_byte_identical():
    """No new pin: K22's listing and table sha256, as K22's committed reading records them; outside_forms.py and
    k22_outside_forms.py are the bytes K22's reading was made by."""
    assert k23.K23_LISTING_SHA256 == REAL_LISTING and k23.K23_TABLE_SHA256 == REAL_TABLES["all"]
    r = json.loads((K22_DIR / "reading.json").read_text())
    assert r["ext"]["listing_sha256"] == k23.K23_LISTING_SHA256 and r["ext"]["tables"] == {"all": k23.K23_TABLE_SHA256}
    pre = json.loads((K22_DIR / "pregate.json").read_text())
    assert pre["ext"]["tables"] == {"all": k23.K23_TABLE_SHA256}
    rc = r["reading_code"]
    assert _sha(T3 / "outside_forms.py") == rc["outside_forms_sha256"]
    assert _sha(T3 / "k22_outside_forms.py") == rc["script_sha256"]
    txt = (K22_DIR / "reading.txt").read_text()
    m = re.search(r"read by experiments/t3_decomp/k22_outside_forms\.py ([0-9a-f]{12}) \(outside_forms "
                  r"([0-9a-f]{12})\)", txt)
    assert m and _sha(T3 / "k22_outside_forms.py")[:12] == m.group(1) and _sha(T3 / "outside_forms.py")[:12] == \
        m.group(2)
    assert m.group(1) == "6482585242ad" and m.group(2) == "4837687d47f1"


def test_the_frozen_t3_files_are_at_their_pins():
    bs = harness._t3_head()
    got = {"t3_decompose.py": bs.sha256(bs.T3D_PATH), "rules.yaml": bs.sha256(bs.RULES_PATH),
           "build_submission.py": bs.sha256(bs.BUILDER_PATH)}
    assert got == ext.FROZEN_SHA256


# ------------------------------------------------------------------------------------------------ (2) refusals


def test_the_run_refuses_without_its_row(tmp_path, capsys):
    p = tmp_path / "TRIES.md"
    p.write_text("| K22 (one row) | T3 | `experiments/t3_decomp/k22_outside_forms.py` |\n")
    with pytest.raises(SystemExit) as e:
        k23.row_present(p)
    assert e.value.code == 2 and "no row" in capsys.readouterr().err
    p.write_text(p.read_text() + "| K23 (one row) | T3 | `experiments/t3_decomp/k23_outside_exponents.py` with |\n")
    assert k23.row_present(p) == "K23"
    p.write_text(p.read_text() + "| K24 (one row) | T3 | `experiments/t3_decomp/k23_outside_exponents.py` again |\n")
    with pytest.raises(SystemExit) as e:
        k23.row_present(p)
    assert e.value.code == 2
    for argv in (["run", "--world", "plain", "--task", "g1", "--out", str(tmp_path / "o")],
                 ["pregate", "--out", str(tmp_path / "o")]):
        with pytest.MonkeyPatch.context() as m:
            m.setattr(k23, "TRIES_MD", tmp_path / "none.md")
            (tmp_path / "none.md").write_text("no rows\n")
            with pytest.raises(SystemExit) as e:
                k23.main(argv)
            assert e.value.code == 2
    assert not (tmp_path / "o").exists()


@pytest.fixture
def row(tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| K23 (one row) | T3 | `experiments/t3_decomp/k23_outside_exponents.py` |\n")
    monkeypatch.setattr(k23, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"})
    return tries


def test_preconditions_refuse_pins_listing_tree_and_k22s_files(tmp_path, monkeypatch, capsys, row):
    argv = ["run", "--world", "plain", "--task", "g1", "--out", str(tmp_path / "o")]
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", None)
    with pytest.raises(SystemExit) as e:
        k23.main(argv)
    assert e.value.code == 2 and "not filled" in capsys.readouterr().err
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", "0" * 64)
    monkeypatch.setattr(ext, "TABLES_FILE", tmp_path / "TABLES.sha256")
    (tmp_path / "TABLES.sha256").write_text("")
    with pytest.raises(SystemExit) as e:
        k23.main(argv)
    assert e.value.code == 2 and "not tracked" in capsys.readouterr().err
    # the listing passes, the pins are not K22's
    monkeypatch.setattr(ext, "read_listing", lambda *a, **k: ({}, "x"))
    monkeypatch.setattr(ext, "K22_TABLE_SHA256", {"all": "1" * 64, "unrelated": None})
    with pytest.raises(SystemExit) as e:
        k23.main(argv)
    assert e.value.code == 2 and "not the listing and `all` sha256 K22 pinned" in capsys.readouterr().err
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", REAL_LISTING)
    monkeypatch.setattr(ext, "K22_TABLE_SHA256", dict(REAL_TABLES))
    monkeypatch.setattr(k22, "frozen_ok", lambda: {})
    monkeypatch.setattr(k23, "K22_RESULTS", tmp_path / "no_k22")
    with pytest.raises(SystemExit) as e:
        k23.main(argv)
    assert e.value.code == 2 and "K22's committed pregate.json" in capsys.readouterr().err
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": False, "machine": "m1", "utc": "x"})
    with pytest.raises(SystemExit) as e:
        k23.main(argv)
    assert e.value.code == 2 and "not clean" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()


def test_a_run_refuses_without_its_pregate_or_with_another_one(tmp_path, capsys):
    rec = {**PINS, "combine_sha256": "c" * 64}
    with pytest.raises(SystemExit) as e:
        k23.load_pregate(tmp_path, rec)
    assert e.value.code == 2 and "pre-gate runs first" in capsys.readouterr().err
    pre = _pre()
    (tmp_path / "pregate.json").write_text(json.dumps({"pregate": pre, "arms": list(k23.active_arms(pre)),
                                                       "ext": {**rec, "combine_sha256": "d"}}))
    with pytest.raises(SystemExit) as e:
        k23.load_pregate(tmp_path, rec)
    assert e.value.code == 2 and "other tables" in capsys.readouterr().err
    both = _pre(dropped={"2": "t", "4": "t"})
    (tmp_path / "pregate.json").write_text(json.dumps({"pregate": both, "arms": list(k23.active_arms(both)),
                                                       "ext": rec}))
    with pytest.raises(SystemExit) as e:
        k23.load_pregate(tmp_path, rec)
    assert e.value.code == 2 and "nothing runs" in capsys.readouterr().err
    (tmp_path / "pregate.json").write_text(json.dumps({"pregate": pre, "arms": list(k23.active_arms(pre)),
                                                       "ext": rec}))
    got, arms = k23.load_pregate(tmp_path, rec)
    assert arms == k23.ARMS and k23._K23 == k23.constants_of(pre) and set(k23._K23) == {"2", "4"}


def test_dropped_exponents_drop_their_arms():
    assert k23.active_arms({"dropped": {}}) == k23.ARMS
    assert k23.active_arms({"dropped": {"4": "x"}}) == ("barred_form", "x_e1", "x_e2", "x_e2_perm", "x_e2_size")
    assert k23.active_arms({"dropped": {"2": "x", "4": "y"}}) == ("barred_form",)
    assert k23.nothing_runs(("barred_form",)) and not k23.nothing_runs(("barred_form", "x_e1", "x_e2"))
    assert set(k23.dropped_arms({"dropped": {"4": "x"}})) == {"x_e4", "x_e4_perm", "x_e4_size"}


# ------------------------------------------------------------------------------------------------ (3) arm factors


def test_the_arm_factors_by_hand():
    phi = np.array([2.0, 0.5, 1.5, 0.25, 4.0, 1.0])
    T = np.array([1.1, 1.0, 0.9, 1.05, 0.8, 1.2])
    basis = np.array([0.1, 0.5, 0.2, 0.9, 0.3, 0.7])
    F = {"phi_barred": phi, "ext": {"all": T, "unrelated": T}, "ext_perm": {"all": ext.permuted(T, basis)}}
    k23._K23.update({"2": {"a": 1.25, "k": 0.875}, "4": {"a": 1.5, "k": 0.75}})
    # one gene by hand: phi 2.0, T 1.1
    assert math.isclose(k23.arm_factor("x_e1", F)[0], 2.2, rel_tol=1e-15)
    assert math.isclose(k23.arm_factor("x_e2", F)[0], 2.0 * 1.21, rel_tol=1e-15)
    assert math.isclose(k23.arm_factor("x_e4", F)[0], 2.0 * 1.4641, rel_tol=1e-14)
    assert math.isclose(k23.arm_factor("x_e2_size", F)[0], 0.875 * 2.0 ** 1.25, rel_tol=1e-15)
    assert math.isclose(k23.arm_factor("x_e4_size", F)[0], 0.75 * 2.0 ** 1.5, rel_tol=1e-15)
    # x_e1 is K22's x_all through the same code path, value for value
    assert np.array_equal(k23.arm_factor("x_e1", F), k22.arm_factor("x_all", F))
    for e in (2.0, 4.0):
        assert np.array_equal(k23.arm_factor(f"x_e{e:g}", F), phi * T ** e)
        # the permuted control: perm(T) ** E = perm(T ** E), K22's permutation
        assert np.array_equal(ext.permuted(T, basis) ** e, ext.permuted(T ** e, basis))
        assert np.array_equal(k23.arm_factor(f"x_e{e:g}_perm", F), phi * ext.permuted(T ** e, basis))
        # a gene absent from the table (exactly 1) keeps phi_M
        assert k23.arm_factor(f"x_e{e:g}", F)[1] == phi[1]
    assert np.array_equal(k23.arm_factor("barred_form", F), phi)
    k23._K23.clear()
    with pytest.raises(RuntimeError, match="not loaded"):
        k23.arm_factor("x_e2_size", F)


def test_the_size_twin_matches_the_product_s_sd_and_den_weighted_total():
    rng = np.random.default_rng(3)
    phi = np.clip(np.exp(rng.normal(0, 0.7, 200)), 0.25, 4.0)
    T = np.where(rng.random(200) < 0.8, np.exp(rng.normal(-0.01, 0.03, 200)), 1.0)
    den = rng.random(200) + 0.05
    loaded = _loaded(T)
    pre = k23.pregate_exponents(phi, loaded, _bases(200, rng), den)
    for e in (2.0, 4.0):
        v = pre["exponents"][f"{e:g}"]
        prod = phi * T ** e
        twin = v["k"] * phi ** v["a"]
        assert math.isclose(np.std(np.log(twin)), np.std(np.log(prod)), rel_tol=1e-12)
        assert math.isclose(den @ twin / den.sum(), den @ prod / den.sum(), rel_tol=1e-12)
        assert math.isclose(v["size"]["den_total"], v["product"]["den_total"], rel_tol=1e-12)
        assert math.isclose(v["a"], np.std(np.log(prod)) / np.std(np.log(phi)), rel_tol=1e-15)
        # recorded to full precision: a json round trip gives the same float
        assert json.loads(json.dumps(v["a"])) == v["a"] and json.loads(json.dumps(v["k"])) == v["k"]


# ------------------------------------------------------------------------------------------------ (4) the fold


@pytest.fixture
def fold(monkeypatch):
    """K18's plain fold fixture (as tests/test_k22_outside_forms.py), the barred factor and the table in a file gene
    order that differs from the reference's."""
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
    t_all = np.where(rng.random(G + 1) < 0.6, np.exp(rng.normal(0, 0.5, G + 1)), 1.0)
    spec, _, _ = k18.barred_spec()
    log = {"construct": [], "assemble": []}

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        assert r is ref and ratio is None and alphas is None
        rows = g.integers(0, N_POOL, n)
        log["construct"].append({"phi": None if phi is None else np.asarray(phi).copy(), "w": w, "row": rows.copy()})
        Cr = C[rows]
        if phi is not None:
            Cr = bs.t3d.apply_gene_factor(Cr, np.asarray(phi) ** w if w != 1.0 else np.asarray(phi), P[rows], g)
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
    k18._FOLD.clear()
    k22._EXT.clear()
    k22._EXT.update(genes=list(FILE_GENES), T={"all": t_all}, record=dict(PINS))
    k23._K23.update({"2": {"a": 1.02, "k": 0.97}, "4": {"a": 1.07, "k": 0.93}})

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    order = lambda v: k18.in_gene_order(v, FILE_GENES, GENES)  # noqa: E731
    yield SimpleNamespace(bs=bs, ctx=ctx, log=log, pool=pool, phi_barred=order(barred_file), t_all=order(t_all),
                          basis=harness.t3_mean_expression(pool))
    k18._FOLD.clear()


def test_barred_form_is_k18s_arm_value_for_value_and_the_names_come_back(fold):
    for seed in (0, 7):
        k18._FOLD.clear()
        want = k18._build_plain(fold.ctx("g1", seed), k18.BARRED)  # K18 as committed
        with k23.k18_arms(k23.ARMS):
            got = k18._build_plain(fold.ctx("g1", seed), k23.BASE)
        assert np.array_equal(got, want)
    assert k18.ARMS == k23._K18_ORIGINAL["ARMS"] and k18.build_arm is k23._K18_ORIGINAL["build_arm"]
    assert k18.fold_factors is k23._K18_ORIGINAL["fold_factors"]
    with pytest.raises(ZeroDivisionError):
        with k23.k18_arms(k23.ARMS):
            assert k18.fold_factors is k22.fold_factors and k18.build_arm is k23.build_arm
            raise ZeroDivisionError
    assert k18.ARMS == k23._K18_ORIGINAL["ARMS"] and k18.build_arm is k23._K18_ORIGINAL["build_arm"]
    assert k18.fold_factors is k23._K18_ORIGINAL["fold_factors"]
    with pytest.raises(RuntimeError):
        with k23.k18_arms(("x_e2",)):
            pass


def test_x_e1_draws_k22s_x_all_and_every_arm_draws_the_barred_s1s_cells(fold):
    k22_x = {}
    with k22.k18_arms(("barred_form", "x_all")):
        fold.log["construct"].clear()
        k22_x["X"] = k18._build_plain(fold.ctx("g2", 3), "x_all")
        k22_x["phi"] = fold.log["construct"][0]["phi"]
    perm = harness.t3_decile_permutation(fold.t_all, fold.basis, 20261113, 10)
    want = {"barred_form": fold.phi_barred, "x_e1": fold.phi_barred * fold.t_all,
            "x_e2": fold.phi_barred * fold.t_all ** 2, "x_e2_perm": fold.phi_barred * perm ** 2,
            "x_e2_size": 0.97 * fold.phi_barred ** 1.02, "x_e4": fold.phi_barred * fold.t_all ** 4,
            "x_e4_perm": fold.phi_barred * perm ** 4, "x_e4_size": 0.93 * fold.phi_barred ** 1.07}
    X = {}
    with k23.k18_arms(k23.ARMS):
        for arm in k23.ARMS:
            fold.log["construct"].clear()
            fold.log["assemble"].clear()
            X[arm] = k18._build_plain(fold.ctx("g2", 3), arm)
            (c,) = fold.log["construct"]
            assert c["w"] == 1.0 and np.allclose(c["phi"], want[arm], rtol=1e-14, atol=0), arm
            fold.log[arm] = (c["row"], fold.log["assemble"][0])
    assert np.array_equal(X["x_e1"], k22_x["X"]) and np.array_equal(fold.log["x_e1"][0], fold.log["x_e1"][0])
    base_rows, base_asm = fold.log[k23.BASE]
    for arm in k23.ARMS:
        rows, asm = fold.log[arm]
        assert np.array_equal(rows, base_rows) and asm["draw"] == base_asm["draw"]
        if arm != k23.BASE:
            assert not np.array_equal(X[arm], X[k23.BASE])


def test_check_world_compares_the_run_s_controls_with_the_pregate_s():
    rng = np.random.default_rng(8)
    n = 120
    phi = np.clip(np.exp(rng.normal(0, 0.6, n)), 0.25, 4.0)
    T = np.where(rng.random(n) < 0.7, np.exp(rng.normal(0, 0.05, n)), 1.0)
    bases = _bases(n, rng)
    pre = k23.pregate_exponents(phi, _loaded(T), bases, rng.random(n) + 0.1)
    k23._K23.update(k23.constants_of(pre))
    F = {"phi_barred": phi, "ext": {"all": T}, "ext_perm": {"all": ext.permuted(T, bases["plain"]["basis"])}}
    k23.check_world(pre, "plain", F, k23.ARMS)
    with pytest.raises(k23._Mismatch):
        k23.check_world(pre, "dm", F, k23.ARMS)  # the dm control was made on another basis


# ------------------------------------------------------------------------------------------------ (5) the pre-gate


def _loaded(T, studies=("GSE204908", "GSE277839"), tiers=("unrelated",)):
    rec = {"studies": list(studies), "tiers": list(tiers), "sign_rule": "applied", "sha256": "a" * 64}
    return {"listing_sha256": "l" * 64, "tables": {"all": (np.asarray(T, dtype=np.float64), rec)},
            "combine": {"all": "a" * 64}, "combine_sha256": "c" * 64}


def _bases(n, rng):
    return {"plain": {"genes": None, "basis": rng.random(n)}, "dm": {"genes": None, "basis": rng.random(n)}}


@pytest.fixture
def small():
    """A table of K22's size (log SD about 0.03, mostly down) on 300 genes, phi_M with a few genes at the clip."""
    rng = np.random.default_rng(11)
    n = 300
    phi = np.clip(np.exp(rng.normal(0, 0.6, n)), 0.25, 4.0)
    T = np.where(rng.random(n) < 0.85, np.exp(rng.normal(-0.01, 0.03, n)), 1.0)
    return SimpleNamespace(phi=phi, T=T, den=rng.random(n) + 0.05, bases=_bases(n, rng), loaded=_loaded(T))


def test_the_pregate_figures_and_no_drop_for_a_small_table(small):
    pre = k23.pregate_exponents(small.phi, small.loaded, small.bases, small.den)
    assert pre["dropped"] == {} and set(pre["exponents"]) == {"2", "4"}
    assert pre["k22"] == ext.pregate(small.phi, small.loaded, None, small.den)
    for e in (2.0, 4.0):
        v = pre["exponents"][f"{e:g}"]
        assert v["table"] == ext.table_figures(small.phi, small.T ** e, small.den)
        assert v["product"] == ext.figures(small.phi, small.phi * small.T ** e, small.den)
        assert set(v["controls"]) == {"plain", "dm"}
        assert v["controls"]["plain"] == ext.figures(
            small.phi, small.phi * ext.permuted(small.T, small.bases["plain"]["basis"]) ** e, small.den)
        assert math.isclose(v["table"]["corr_ln_phi_m"], pre["k22"]["tables"]["all"]["corr_ln_phi_m"], rel_tol=1e-12)
    assert pre["constants"]["exponents"] == [2.0, 4.0]
    assert k23.active_arms(pre) == k23.ARMS
    lines = "\n".join(k23.pregate_lines(pre))
    assert "e = 2: a_E" in lines and "control (dm basis)" in lines and "exponents dropped: none" in lines


def test_world_gene_orders_are_applied_before_the_permutation(small):
    n = small.phi.size
    genes = [f"g{i}" for i in range(n)]
    rev = genes[::-1]
    bases = {"plain": {"genes": rev, "basis": small.bases["plain"]["basis"]}}
    pre = k23.pregate_exponents(small.phi, small.loaded, bases, small.den, genes)
    phi_w, T_w, den_w = (np.ascontiguousarray(v[::-1]) for v in (small.phi, small.T, small.den))
    want = ext.figures(phi_w, phi_w * ext.permuted(T_w, bases["plain"]["basis"]) ** 2, den_w)
    assert pre["exponents"]["2"]["controls"]["plain"] == want


@pytest.mark.parametrize("bound", ["sd_table", "at_clip", "den", "sd_ratio"])
def test_each_bound_alone_drops_the_larger_exponent(small, monkeypatch, bound):
    relaxed = {"MAX_SD_LN_TABLE": 1e9, "MAX_AT_CLIP": 10 ** 9, "MAX_DEN_SHIFT": 1e9, "MAX_SD_RATIO": 1e9}
    for k, v in relaxed.items():
        monkeypatch.setattr(k23, k, v)
    if bound == "at_clip":  # five genes just under the clip that the table moves up: at e = 4 they pass it
        small.phi[:5], small.T[:5] = 3.7, 1.03
        small.loaded = _loaded(small.T)
    free = k23.pregate_exponents(small.phi, small.loaded, small.bases, small.den)
    assert free["dropped"] == {}
    e2, e4 = free["exponents"]["2"], free["exponents"]["4"]
    if bound == "sd_table":
        a, b = e2["table"]["sd_ln"], e4["table"]["sd_ln"]
        name, why = "MAX_SD_LN_TABLE", "SD(ln T ** 4)"
    elif bound == "at_clip":
        a = max(e2["product"]["genes_at_clip"], *(c["genes_at_clip"] for c in e2["controls"].values()))
        b = e4["product"]["genes_at_clip"]
        name, why = "MAX_AT_CLIP", "the product has"
    elif bound == "den":
        a, b = abs(e2["den_shift"]), abs(e4["den_shift"])
        name, why = "MAX_DEN_SHIFT", "den-weighted total moves"
    else:
        a, b = e2["sd_ln_product"] / e2["sd_ln_phi_m"], e4["sd_ln_product"] / e4["sd_ln_phi_m"]
        name, why = "MAX_SD_RATIO", "SD(ln product)"
    assert b > a, (bound, a, b)
    monkeypatch.setattr(k23, name, int(a) if name == "MAX_AT_CLIP" else (a + b) / 2)
    pre = k23.pregate_exponents(small.phi, small.loaded, small.bases, small.den)
    assert set(pre["dropped"]) == {"4"} and why in pre["dropped"]["4"], pre["dropped"]
    assert k23.active_arms(pre) == ("barred_form", "x_e1", "x_e2", "x_e2_perm", "x_e2_size")
    assert k23.constants_of(pre) == {"2": {"a": e2["a"], "k": e2["k"]}}


def test_the_clip_bound_tripped_by_the_control_only(monkeypatch):
    """phi_M at the clip on ten genes the table moves down: the product leaves the clip, its permutation does not."""
    rng = np.random.default_rng(4)
    n = 200
    phi = np.clip(np.exp(rng.normal(0, 0.4, n)), 0.5, 2.0)
    phi[:10], phi[10:20] = 4.0, 3.9
    T = np.ones(n)
    T[:10] = np.linspace(0.85, 0.94, 10)
    T[50:] = np.exp(rng.normal(0, 0.03, n - 50))
    loaded = _loaded(T)
    bases = {"plain": {"genes": None, "basis": np.arange(n, dtype=float)}}
    pre = k23.pregate_exponents(phi, loaded, bases, np.ones(n))
    assert "all" not in pre["k22"]["dropped"] and pre["dropped"] == {}
    assert pre["exponents"]["2"]["product"]["genes_at_clip"] == 0
    assert pre["exponents"]["2"]["controls"]["plain"]["genes_at_clip"] > 0
    monkeypatch.setattr(k23, "MAX_AT_CLIP", 0)
    pre = k23.pregate_exponents(phi, loaded, bases, np.ones(n))
    assert set(pre["dropped"]) == {"2", "4"}
    assert "the control in world plain" in pre["dropped"]["2"] and "the product has" not in pre["dropped"]["2"]
    assert k23.active_arms(pre) == ("barred_form",) and k23.nothing_runs(k23.active_arms(pre))


def test_a_drop_at_e2_drops_e4(small, monkeypatch):
    original = ext.product

    def odd(phi_m, x, T, e=1.0):  # the product at e = 2 alone moves the den-weighted total by a fifth
        out = original(phi_m, x, T, e)
        return out * 1.2 if e == 2.0 else out

    monkeypatch.setattr(ext, "product", odd)
    pre = k23.pregate_exponents(small.phi, small.loaded, small.bases, small.den)
    assert "den-weighted total" in pre["dropped"]["2"]
    assert pre["dropped"]["4"] == "a smaller exponent is dropped"
    assert k23.active_arms(pre) == ("barred_form",)


@pytest.mark.parametrize("case,why", [("few", "genes != 1"), ("one_study", "study"), ("cardiac", "cardiac"),
                                      ("corr", "corr")])
def test_k22s_drops_drop_both_exponents(case, why):
    rng = np.random.default_rng(2)
    n = 100
    phi = np.clip(np.exp(rng.normal(0, 0.6, n)), 0.25, 4.0)
    T = np.where(rng.random(n) < 0.8, np.exp(rng.normal(0, 0.03, n)), 1.0)
    kw = {}
    if case == "few":
        T = np.ones(n)
        T[:5] = 1.1
    elif case == "one_study":
        kw = {"studies": ("GSE204908",)}
    elif case == "cardiac":
        kw = {"tiers": ("cardiac",)}
    else:
        T = np.where(T != 1.0, phi ** 0.05, 1.0)
    loaded = _loaded(T, **kw)
    pre = k23.pregate_exponents(phi, loaded, _bases(n, rng), rng.random(n) + 0.1)
    assert "all" in pre["k22"]["dropped"] and why in pre["k22"]["dropped"]["all"]
    assert set(pre["dropped"]) == {"2", "4"} and all("K22's pre-gate drops" in v for v in pre["dropped"].values())
    assert k23.active_arms(pre) == ("barred_form",) and pre["exponents"] == {}


def test_the_pregate_command_writes_its_file_and_exits_5_when_both_exponents_drop(tmp_path, monkeypatch, capsys):
    import vec.prep

    genes = [f"g{i}" for i in range(G)]
    C = np.random.default_rng(0).poisson(3.0, size=(30, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    wt = ad.AnnData(X=np.log1p(C).astype(np.float32), var=pd.DataFrame(index=genes), layers={"counts": C})
    rng = np.random.default_rng(1)
    phi = np.clip(np.exp(rng.normal(0, 0.6, G)), 0.25, 4.0)
    T = np.where(rng.random(G) < 0.9, np.exp(rng.normal(0, 0.03, G)), 1.0)
    T[:11] = np.exp(rng.normal(0, 0.03, 11))
    seen: dict = {"T": T}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=phi, genes=genes, spec={}, record={"file": "f"})
        return {"file": "f"}

    def load_ext(g, phi_m, basis=None, den=None):
        loaded = _loaded(seen["T"])
        pre = ext.pregate(phi_m, loaded, basis, den)
        seen.update(den=den)
        k22._EXT.clear()
        k22._EXT.update(genes=list(g), T={"all": seen["T"]}, loaded=loaded, pre=pre,
                        record={**PINS, "combine_sha256": "c" * 64}, table_records={})
        return pre

    monkeypatch.setattr(k23, "preconditions", lambda: ("K23", {"commit": "abc", "tree_clean": True,
                                                               "machine": "m1", "utc": "x"}))
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(vec.prep, "load_stage", lambda p, drop_duplicates=False, genes=None: wt)
    monkeypatch.setattr(k22, "load_ext", load_ext)
    monkeypatch.setattr(k23, "world_bases", lambda data, g, bs: {
        "plain": {"genes": list(g), "basis": np.arange(G, dtype=float)},
        "dm": {"genes": list(g)[::-1], "basis": np.arange(G, dtype=float)}})
    k23.main(["pregate", "--out", str(tmp_path / "a")])
    body = json.loads((tmp_path / "a" / "pregate.json").read_text())
    assert body["arms"] == list(k23.ARMS) and body["row"] == "K23" and body["ext"]["listing_sha256"] == \
        k23.K23_LISTING_SHA256
    assert set(body["pregate"]["exponents"]) == {"2", "4"} and "k22" in body["pregate"]
    assert np.allclose(seen["den"], (C / C.sum(1, keepdims=True)).mean(0))
    out = capsys.readouterr().out
    assert "pre-gate (the table, phi_M and the worlds' bases alone; no score)" in out
    with pytest.raises(SystemExit) as e:
        k23.main(["pregate", "--out", str(tmp_path / "a")])
    assert e.value.code == 2
    few = np.ones(G)
    few[:3] = 1.1
    seen["T"] = few
    with pytest.raises(SystemExit) as e:
        k23.main(["pregate", "--out", str(tmp_path / "b")])
    assert e.value.code == 5 and (tmp_path / "b" / "pregate.json").exists()
    assert "nothing runs" in capsys.readouterr().out


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
ALL_CELLS = [(w, n) for w in k23.WORLDS for n in k23.TASKS]


def _table(**arms) -> dict:
    return {harness.FLOOR: 50.0, **{a: BASE_SCORE for a in k23.ARMS}, **arms}


def _per(cells: dict, gated: dict | None = None, arms=k23.ARMS, default=None) -> dict:
    out: dict = {}
    for w in k23.WORLDS:
        a = (harness.FLOOR, *arms) if w == "plain" else arms
        out[w] = {n: k23._summarise(_rows({x: cells.get((w, n), default or _table())[x] for x in a},
                                          (gated or {}).get((w, n))), a) for n in k23.TASKS}
    return out


def _everywhere(**arms) -> dict:
    return {k: _table(**arms) for k in ALL_CELLS}


GOOD2 = {"x_e2": 57.0}  # +1 over barred_form, its control, its twin and x_e1 in every cell


def test_the_summaries_carry_the_extra_pairs():
    per = _per(_everywhere(**GOOD2))
    p = per["plain"]["g1"]["paired"]
    for k in ("x_e2 - barred_form", "x_e2 - x_e2_perm", "x_e2 - x_e2_size", "x_e2 - x_e1", "x_e4 - x_e1",
              "x_e4 - x_e2", "x_e4 - x_e4_size", "x_e1 - barred_form", "x_e2_size - barred_form",
              "barred_form - floor"):
        assert k in p, k
    assert p["x_e2 - x_e1"]["mean"] == pytest.approx(1.0) and p["x_e4 - x_e2"]["mean"] == pytest.approx(-1.0)
    assert "x_e4 - x_e2" not in k23._summarise(_rows({a: 56.0 for a in ("barred_form", "x_e1", "x_e2")}),
                                               ("barred_form", "x_e1", "x_e2"))["paired"]


def test_only_e2_qualifies_and_is_named():
    res = k23.decide(_per(_everywhere(**GOOD2)), k23.ARMS)
    s = res["status"]["x_e2"]
    assert s["qualified"] and all(s[c] for c in "abcdef") and res["names"] == ["x_e2"]
    assert res["qualified"] == ["x_e2"] and res["rule"] == "only x_e2 qualified"
    assert res["files"] == ["xall-B --table-exponent 2.0 (base B of decision log 2026-11-20 (q))"]
    assert "--table-exponent 2.0" in res["outcome"] and "after (q)'s C list" in res["outcome"]
    assert not res["status"]["x_e4"]["qualified"]


def _fails_alone(per, letter):
    s = k23.decide(per, k23.ARMS)["status"]["x_e2"]
    assert not s["qualified"] and not s[letter], s
    assert all(s[c] for c in "abcdef" if c != letter), s


def test_each_qualification_clause_failing_alone():
    up = _table(**GOOD2)
    # (a) the gate clause: AHEAD in plain g1 and dm g2 only, with one more gated seed than barred_form in dm g2
    _fails_alone(_per({("plain", "g1"): up, ("dm", "g2"): up}, {("dm", "g2"): {"x_e2": 1}}), "a")
    # ... the same cells with equal gated counts qualify
    assert k23.decide(_per({("plain", "g1"): up, ("dm", "g2"): up}), k23.ARMS)["status"]["x_e2"]["qualified"]
    # (a) plain only (then (d)'s dm clause fails with it: no dm cell is AHEAD)
    s = k23.decide(_per({("plain", "g1"): up, ("plain", "g2"): up}), k23.ARMS)["status"]["x_e2"]
    assert not s["qualified"] and not s["a"] and not s["d"] and s["b"] and s["c"] and s["e"] and s["f"]
    # (b) behind in one cell
    _fails_alone(_per(_everywhere(**GOOD2) | {("dm", "g1"): _table(x_e2=55.5)}), "b")
    # (c) two more gated seeds than barred_form in plain g1
    _fails_alone(_per(_everywhere(**GOOD2), {("plain", "g1"): {"x_e2": 2}}), "c")
    # (d) its permuted control as high in one AHEAD cell
    _fails_alone(_per(_everywhere(**GOOD2) | {("plain", "g2"): _table(x_e2=57.0, x_e2_perm=56.9)}), "d")
    # (d) no dm cell where it beats its control
    cells = {("plain", "g1"): up, ("plain", "g2"): up, ("dm", "g1"): _table(x_e2=57.0, x_e2_perm=56.9),
             ("dm", "g2"): _table(x_e2=57.0, x_e2_perm=56.9)}
    _fails_alone(_per(cells), "d")
    # (e) its size twin as high in one AHEAD cell
    _fails_alone(_per(_everywhere(**GOOD2) | {("dm", "g2"): _table(x_e2=57.0, x_e2_size=56.8)}), "e")
    # (f) not above x_e1 in one AHEAD cell
    _fails_alone(_per(_everywhere(**GOOD2) | {("plain", "g1"): _table(x_e2=57.0, x_e1=57.0)}), "f")


def test_e4_needs_its_dose_step_over_e2():
    res = k23.decide(_per(_everywhere(x_e2=57.0, x_e4=57.0)), k23.ARMS)
    s4 = res["status"]["x_e4"]
    assert not s4["qualified"] and not s4["f"] and all(s4[c] for c in "abcde")
    assert res["names"] == ["x_e2"]


def test_both_qualify_and_the_direct_pair_names_e4():
    res = k23.decide(_per(_everywhere(x_e2=57.0, x_e4=58.0)), k23.ARMS)
    assert res["qualified"] == ["x_e2", "x_e4"] and res["names"] == ["x_e4"] and res["pivot"]["e4"]
    assert res["files"] == ["xall-B --table-exponent 4.0 (base B of decision log 2026-11-20 (q))"]


def test_an_e4_lead_from_dm_g1_alone_names_e2():
    cells = _everywhere(x_e2=57.0, x_e4=57.2) | {("dm", "g1"): _table(x_e2=57.0, x_e4=59.0)}
    res = k23.decide(_per(cells), k23.ARMS)
    assert res["qualified"] == ["x_e2", "x_e4"] and res["pivot"]["ahead_cells"] == ["dm g1"]
    assert not res["pivot"]["e4"] and res["names"] == ["x_e2"]


def test_an_e4_pair_below_the_bar_in_one_cell_names_e2():
    cells = _everywhere(x_e2=57.0, x_e4=58.0) | {("plain", "g2"): _table(x_e2=57.5, x_e4=56.2)}
    res = k23.decide(_per(cells), k23.ARMS)
    assert res["qualified"] == ["x_e2", "x_e4"] and res["pivot"]["below_cells"] == ["plain g2"]
    assert not res["pivot"]["e4"] and res["names"] == ["x_e2"]


def test_only_e4_with_e2_pointing_the_same_way_names_e4():
    cells = _everywhere(x_e2=57.0, x_e4=58.0) | {("plain", "g1"): _table(x_e2=57.0, x_e2_perm=57.0, x_e4=58.0)}
    res = k23.decide(_per(cells), k23.ARMS)
    assert res["qualified"] == ["x_e4"] and not res["status"]["x_e2"]["d"] and res["names"] == ["x_e4"]
    assert res["rule"] == "only x_e4 qualified; x_e2 points the same way"


@pytest.mark.parametrize("why", ["negative", "behind", "gated"])
def test_only_e4_without_e2_pointing_the_same_way_names_nothing(why):
    cells = _everywhere(x_e2=57.0, x_e4=58.0) | {("plain", "g1"): _table(x_e2=57.0, x_e2_perm=57.0, x_e4=58.0)}
    gated, arms = None, k23.ARMS
    if why == "negative":
        cells[("dm", "g2")] = _table(x_e2=55.9, x_e4=58.0)  # x_e2 below barred_form (not BEHIND) where x_e4 is AHEAD
    elif why == "behind":
        cells[("plain", "g2")] = _table(x_e2=55.5, x_e4=57.0)
    else:
        gated = {("plain", "g2"): {"x_e2": 2}}
    res = k23.decide(_per(cells, gated, arms=arms), arms)
    assert res["qualified"] == ["x_e4"] and res["names"] == [] and res["files"] == []
    assert res["outcome"] == k23.NOT_NAMED_E4 == ("e = 4 alone qualified; not named: the smaller exponent does not "
                                                  "point the same way")


def test_nothing_qualifies_and_controls_twins_and_x_e1_are_never_named():
    res = k23.decide(_per(_everywhere(x_e1=60.0, x_e2_perm=60.0, x_e4_perm=60.0, x_e2_size=60.0, x_e4_size=60.0)),
                     k23.ARMS)
    assert res["names"] == [] and res["outcome"] == "no file: nothing qualifies; T3's list is unchanged"
    assert set(res["status"]) == set(k23.NAMEABLE)
    arms = ("barred_form", "x_e1", "x_e2", "x_e2_perm", "x_e2_size")
    res = k23.decide(_per(_everywhere(**GOOD2), arms=arms), arms)
    assert res["status"]["x_e4"] == {"ran": False, "qualified": False} and res["names"] == ["x_e2"]


def test_never_two_names():
    rng = np.random.default_rng(21)
    for _ in range(150):
        cells = {k: _table(**{a: float(BASE_SCORE + rng.choice([-0.5, 0.0, 0.5, 1.0, 2.0]))
                              for a in k23.ARMS if a != k23.BASE}) for k in ALL_CELLS}
        res = k23.decide(_per(cells), k23.ARMS)
        assert len(res["names"]) <= 1 and set(res["names"]) <= set(k23.NAMEABLE)


# ------------------------------------------------------------------------------------------------ (7) a null


def test_on_null_seeds_with_paired_sd_growing_with_the_exponent_a_name_is_rare():
    """Every arm at barred_form's level, its per-seed noise SD 0.16 x E (x_e1 at E = 1): K22's control's paired SD at
    E = 1, scaled by the exponent. The rule names a file in under 2% of 500 draws."""
    rng = np.random.default_rng(23)
    named, n = 0, 500
    arms_noise = {a: 0.16 * k23.EXP_OF[a] for a in k23.ARMS if a != k23.BASE}
    for _ in range(n):
        per = {}
        for w in k23.WORLDS:
            a = (harness.FLOOR, *k23.ARMS) if w == "plain" else k23.ARMS
            per[w] = {}
            for t in k23.TASKS:
                noise = {x: rng.normal(0, sd, len(SEEDS)) for x, sd in arms_noise.items()}
                per[w][t] = k23._summarise(_rows({x: _table()[x] for x in a}, noise=noise), a)
        named += bool(k23.decide(per, k23.ARMS)["names"])
    assert named / n < 0.02, named


# ------------------------------------------------------------------------------------------------ (8) the reading


def _prov(commit="abc", clean=True, machine="m1"):
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "x"}


def _base_pre() -> dict:
    rng = np.random.default_rng(31)
    n = 60
    phi = np.clip(np.exp(rng.normal(0, 0.6, n)), 0.25, 4.0)
    T = np.where(rng.random(n) < 0.8, np.exp(rng.normal(-0.01, 0.03, n)), 1.0)
    return json.loads(json.dumps(k23.pregate_exponents(phi, _loaded(T), _bases(n, rng), rng.random(n) + 0.1)))


def _pre(dropped=None, k22_rel=0.0) -> dict:
    pre = _base_pre()
    pre["dropped"] = dict(dropped or {})
    pre["k22"]["products"]["all"]["sd_ln"] *= 1 + k22_rel
    return pre


def _write(root: Path, cells=None, commit="abc", clean=True, machine="m1", shift=0.0, drop=None, dropped=None,
           ext_rec=None, bad_ext_file=None, k22_shift=0.0, k22_gated=None, k22_rel=0.0, bad_exp_file=None):
    cells = cells or {}
    pre = _pre(dropped)
    arms = k23.active_arms(pre)
    erec = ext_rec or {**PINS, "combine_sha256": "c" * 64}
    root.mkdir(parents=True)
    rec = {"dropped": k23.dropped_arms(pre), "arms": list(arms), "ext": erec, "exponents": k23.constants_of(pre)}
    k18.write_task_file(root / "pregate.json", {"pregate": _pre(dropped, k22_rel), "arms": list(arms),
                                                "ext": erec, "provenance": _prov(commit, clean, machine)})
    refs, k22_worlds = {}, {"plain": {}, "dm": {}}
    for w in k23.WORLDS:
        for n in k23.TASKS:
            a = (harness.FLOOR, *arms) if w == "plain" else arms
            rows = _rows({x: cells.get((w, n), _table())[x] for x in a})
            refs[(w, n)] = [{"seed": r["seed"], "scores": {k23.BASE: r["scores"][k23.BASE] - shift}} for r in rows]
            if "x_e1" in arms:
                k22_worlds[w][n] = {"arms": {"x_all": {"mean": float(np.mean([r["scores"]["x_e1"] for r in rows]))
                                                       + k22_shift, "severity_gated": k22_gated or 0}}}
            for r in rows:
                s = r["seed"]
                f = root / (f"plain_{n}/seed_{s}.json" if w == "plain" else f"dm/{n}_seed_{s}.json")
                f.parent.mkdir(parents=True, exist_ok=True)
                rel = str(f.relative_to(root))
                if drop == rel:
                    continue
                rr = rec
                if bad_ext_file == rel:
                    rr = {**rec, "ext": {**erec, "combine_sha256": "d"}}
                if bad_exp_file == rel:
                    rr = {**rec, "exponents": {**rec["exponents"], "2": {"a": 1.0, "k": 1.0}}}
                body = ({"task": {"key": k23.TASKS[n].key}, "rows": [r], "seed": s, "pregate": rr,
                         "provenance": _prov(commit, clean, machine)} if w == "plain" else
                        {**r, "check": {k18.DM_CHECK: 60.0 + s}, "target": k18.DM_TARGETS[n], "task": n,
                         "pregate": rr, "arms": list(arms), "provenance": _prov(commit, clean, machine)})
                k18.write_task_file(f, body)
    k18_dir = root.parent / f"k18_{root.name}"
    k18_dir.mkdir(parents=True, exist_ok=True)
    for n in k23.TASKS:
        (k18_dir / f"plain_{n}.json").write_text(json.dumps({"rows": refs[("plain", n)]}))
    (k18_dir / "dm.json").write_text(json.dumps({"rows": {n: refs[("dm", n)] for n in k23.TASKS}}))
    pins = {f: hashlib.sha256((k18_dir / f).read_bytes()).hexdigest() for f in k23.K18_SHA256}
    k22_dir = root.parent / f"k22_{root.name}"
    k22_dir.mkdir()
    (k22_dir / "pregate.json").write_text(json.dumps({"pregate": _pre()["k22"]}))
    (k22_dir / "reading.json").write_text(json.dumps({"ext": {**PINS, "combine_sha256": "c" * 64},
                                                      "worlds": k22_worlds}))
    return SimpleNamespace(root=root, refs=k18_dir, pins=pins, k22=k22_dir)


@pytest.fixture
def u3d(tmp_path):
    rows = [{"target": t, "seed": s, "world": "depth_matched", "variant": "dm:transfer:mab21l2_all", "t3": 60.0 + s}
            for t in ("rep1", "rep2") for s in SEEDS]
    p = tmp_path / "orth_records.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def _read(w, u3d):
    return k23.reading(w.root, w.refs, u3d, w.k22)


def test_the_reading_of_a_run_folder(tmp_path, u3d, monkeypatch):
    w = _write(tmp_path / "ok", _everywhere(**GOOD2))
    monkeypatch.setattr(k23, "K18_SHA256", w.pins)
    res = _read(w, u3d)
    assert not res["void"] and res["names"] == ["x_e2"] and res["tries_counted"] == 6, res.get("why")
    assert all(d["max_abs_diff"] == 0.0 for d in res["guard"].values()) and len(res["files_sha256"]) == 41
    assert all(g["passes"] for g in res["guard2"].values()) and len(res["guard2"]) == 4
    txt = k23.text(res)
    for s in ("WORLD plain", "x_e2 - x_e2_size", "x_e2 - x_e1", "x_e4 - x_e2", "guard 2", "pre-gate",
              "--table-exponent 2.0", "barred_form - floor"):
        assert s in txt, s


def test_a_dropped_exponent_reads_with_three_tries(tmp_path, u3d, monkeypatch):
    w = _write(tmp_path / "d", _everywhere(**GOOD2), dropped={"4": "test"})
    monkeypatch.setattr(k23, "K18_SHA256", w.pins)
    res = _read(w, u3d)
    assert not res["void"] and res["tries_counted"] == 3 and res["names"] == ["x_e2"]
    assert res["arms"] == ["barred_form", "x_e1", "x_e2", "x_e2_perm", "x_e2_size"]
    assert res["status"]["x_e4"]["ran"] is False


@pytest.mark.parametrize("kw,why,tries", [
    ({"shift": 0.2}, "guard 1", 7), ({"k22_shift": 0.11}, "guard 2", 7), ({"k22_gated": 1}, "guard 2", 7),
    ({"clean": False}, "dirty tree", 6), ({"bad_ext_file": "dm/g1_seed_3.json"}, "tables differ", 6),
    ({"bad_exp_file": "plain_g2/seed_4.json"}, "a_E or k_E", 6),
    ({"ext_rec": {"listing_sha256": "0" * 64, "tables": PINS["tables"], "combine_sha256": "c" * 64}},
     "not the pinned ones", 6),
    ({"k22_rel": 1e-9}, "K22's product figures", 6)])
def test_void_readings(tmp_path, u3d, monkeypatch, kw, why, tries):
    w = _write(tmp_path / "v", _everywhere(**GOOD2), **kw)
    monkeypatch.setattr(k23, "K18_SHA256", w.pins)
    res = _read(w, u3d)
    assert res["void"] and any(why in x for x in res["why"]), res["why"]
    assert res["names"] == [] and "no file" in res["outcome"] and "VOID" in k23.text(res)
    assert res["tries_counted"] == tries


def test_k22s_figures_agree_within_the_relative_tolerance(tmp_path, u3d, monkeypatch):
    w = _write(tmp_path / "t", k22_rel=1e-13)
    monkeypatch.setattr(k23, "K18_SHA256", w.pins)
    assert not _read(w, u3d)["void"]


def test_a_missing_k22_record_voids(tmp_path, u3d, monkeypatch):
    w = _write(tmp_path / "m")
    (w.k22 / "reading.json").unlink()
    monkeypatch.setattr(k23, "K18_SHA256", w.pins)
    res = _read(w, u3d)
    assert res["void"] and any("K22's committed reading.json is missing" in x for x in res["why"])


def test_read_writes_once_records_its_code_and_never_overwrites(tmp_path, u3d, monkeypatch, capsys):
    w = _write(tmp_path / "r")
    monkeypatch.setattr(k23, "K18_SHA256", w.pins)
    monkeypatch.setattr(k23, "reading", lambda runs, f=k23.reading: f(runs, w.refs, u3d, w.k22))
    out = tmp_path / "reading"
    k23.main(["read", "--runs", str(w.root), "--out", str(out)])
    body = json.loads((out / "reading.json").read_text())
    rc = body["reading_code"]
    assert rc["script_sha256"] == _sha(T3 / "k23_outside_exponents.py")
    assert rc["outside_forms_sha256"] == _sha(T3 / "outside_forms.py") and rc["k22_sha256"] == _sha(
        T3 / "k22_outside_forms.py") and rc["k21_sha256"] == _sha(T3 / "k21_factor_forms.py")
    assert body["outcome"] == "no file: nothing qualifies; T3's list is unchanged"
    first = (out / "reading.json").read_bytes()
    with pytest.raises(SystemExit) as e:
        k23.main(["read", "--runs", str(w.root), "--out", str(out)])
    assert e.value.code == 4 and (out / "reading.json").read_bytes() == first
    capsys.readouterr()


def test_a_missing_file_is_not_a_reading(tmp_path, u3d, monkeypatch, capsys):
    w = _write(tmp_path / "i", drop="plain_g2/seed_7.json")
    monkeypatch.setattr(k23, "K18_SHA256", w.pins)
    (w.root / "dm" / "g1_seed_0.json").write_text("{ not JSON: never opened")
    res = _read(w, u3d)
    assert res["incomplete"] and res["missing"] == ["plain_g2/seed_7.json"]
    monkeypatch.setattr(k23, "reading", lambda runs, f=k23.reading: f(runs, w.refs, u3d, w.k22))
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k23.main(["read", "--runs", str(w.root), "--out", str(out)])
    assert e.value.code == 3 and not out.exists() and "not a reading" in capsys.readouterr().out


# ------------------------------------------------------------------------------------------------ (9) the plain run


@pytest.fixture
def fake_run(monkeypatch, tmp_path):
    calls = {"seeds": []}
    lines = {}
    rng = np.random.default_rng(1)
    n = 40
    phi = np.clip(np.exp(rng.normal(0, 0.6, n)), 0.25, 4.0)
    T = np.where(rng.random(n) < 0.8, np.exp(rng.normal(0, 0.03, n)), 1.0)
    basis = rng.random(n)
    pre = k23.pregate_exponents(phi, _loaded(T), {"plain": {"genes": None, "basis": basis},
                                                  "dm": {"genes": None, "basis": rng.random(n)}}, np.ones(n))
    record = {**PINS, "combine_sha256": "c" * 64}
    F = {"phi_barred": phi, "ext": {"all": T}, "ext_perm": {"all": ext.permuted(T, basis)}}

    def run_task(t, cands, seeds, data_dir, log=print):
        calls["seeds"] += list(seeds)
        calls["names"] = [c.name for c in cands]
        k18._FOLD["F"] = F
        for s in seeds:
            lines[s] = f"{t.key} seed {s}: floor 50.00, barred_form {60 - s / 9:.2f}, x_e2 {61 + s / 7:.2f}"
            log(lines[s])
        return {"task": {"key": t.key}, "rows": [{"seed": s} for s in seeds],
                "provenance": {"commit": "abc", "tree_clean": True}}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=phi, genes=[f"g{i}" for i in range(n)], spec={}, record={})
        return {}

    def load_ext(genes, phi_m, basis=None, den=None):
        k22._EXT.clear()
        k22._EXT.update(genes=list(genes), T={"all": T}, record=dict(record), pre=pre["k22"])
        return pre["k22"]

    out = tmp_path / "run"
    out.mkdir()
    (out / "pregate.json").write_text(json.dumps({"pregate": pre, "arms": list(k23.active_arms(pre)),
                                                  "ext": record}))
    monkeypatch.setattr(k23, "preconditions", lambda: ("K23", k18.stamp()))
    monkeypatch.setattr(k22, "load_ext", load_ext)
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(k21.k7, "_PREP", {"fold": {"record": {}}})
    monkeypatch.setattr(k18, "stamp", lambda: {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"})
    yield SimpleNamespace(calls=calls, lines=lines, out=out, pre=pre)
    k18._FOLD.clear()


def test_the_plain_run_is_seed_by_seed_resumable_and_logs_progress_only(fake_run, capsys):
    out = fake_run.out
    k23.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0-1", "--out", str(out)])
    assert fake_run.calls["seeds"] == [0, 1] and fake_run.calls["names"] == [harness.FLOOR, *k23.ARMS]
    printed = capsys.readouterr()
    for x in printed.out.splitlines():
        assert re.fullmatch(r"K23 plain g1 seed \d: scored|wrote \S+/seed_\d\.json: commit abc, clean True, "
                            r"machine m1", x), x
    for line in fake_run.lines.values():
        for score in re.findall(r"\d+\.\d\d", line):
            assert score not in printed.out and score not in printed.err
    body = json.loads((out / "plain_g1" / "seed_0.json").read_text())
    assert body["pregate"]["arms"] == list(k23.ARMS) and body["pregate"]["ext"]["listing_sha256"] == \
        k23.K23_LISTING_SHA256
    assert body["pregate"]["exponents"] == k23.constants_of(fake_run.pre) and set(body["world_figures"]) == set(
        k23.ARMS[1:])
    k23.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0-2", "--out", str(out)])
    assert fake_run.calls["seeds"] == [0, 1, 2] and "seed 0: present, skipped" in capsys.readouterr().out
    assert k18.ARMS == k23._K18_ORIGINAL["ARMS"] and k18.build_arm is k23._K18_ORIGINAL["build_arm"]


def test_a_crash_is_masked_and_writes_nothing(fake_run, monkeypatch, capsys):
    def boom(*a, **k):
        raise RuntimeError("score 61.234 went wrong")

    monkeypatch.setattr(harness, "run_task", boom)
    with pytest.raises(SystemExit) as e:
        k23.main(["run", "--world", "plain", "--task", "g1", "--seeds", "0", "--out", str(fake_run.out)])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "61.234" not in err and "<num>" in err
    assert not (fake_run.out / "plain_g1" / "seed_0.json").exists()
    assert k18.ARMS == k23._K18_ORIGINAL["ARMS"]


def test_a_world_whose_control_is_not_the_pregate_s_writes_nothing(fake_run, monkeypatch, capsys):
    body = json.loads((fake_run.out / "pregate.json").read_text())
    body["pregate"]["exponents"]["2"]["controls"]["plain"]["sd_ln"] += 1e-6
    (fake_run.out / "pregate.json").write_text(json.dumps(body))
    with pytest.raises(SystemExit) as e:
        k23.main(["run", "--world", "plain", "--task", "g2", "--seeds", "0", "--out", str(fake_run.out)])
    assert e.value.code == 1 and "not the pre-gate's" in capsys.readouterr().err
    assert not (fake_run.out / "plain_g2" / "seed_0.json").exists()


# ------------------------------------------------------------------------------------------------ (10) the launchers


def test_the_launcher_constants():
    sh = (T3 / "run_k23.sh").read_text()
    k22sh = (T3 / "run_k22.sh").read_text()
    assert "SCRIPT=experiments/t3_decomp/k23_outside_exponents.py" in sh and "k22_outside_forms.py`" not in sh
    assert "S=~/v/scratch/k23_${COMMIT:0:7}" in sh and "R=runs/k23_${COMMIT:0:7}" in sh
    assert "OUT=experiments/t3_decomp/results/k23_outside_exponents" in sh
    assert "grep -q '^| .*`experiments/t3_decomp/k23_outside_exponents.py`' experiments/harness/TRIES.md" in sh
    assert "TABLES=$HOME/v/vec-ext/t3/tables  # combine.py's --out; fixed, no override" in sh
    assert "VEC_T3_EXT_TABLES=\"$TABLES\"" in sh
    for f in ("experiments/t3_ext/TABLES.sha256", "experiments/t3_decomp/results/k22_outside_forms/pregate.json",
              "experiments/t3_decomp/results/k22_outside_forms/reading.json"):
        assert f in sh
    assert "tests/test_k23_outside_exponents.py tests/test_k22_outside_forms.py tests/test_t3_outside_forms.py" in sh
    assert "tests/test_k21_factor_forms.py tests/test_k18_built_files.py" in sh
    assert "^K23 plain g[12] seed [0-9]+: (scored|present, skipped)$" in sh
    for code in (10, 20, 22, 25, 30, 41, 43, 44, 95, 96, 97, 98):
        assert f"finish {code}" in sh, code
    assert "echo 99 > \"$S/all.exit\"; exit 99" in sh
    assert "if [ $rc = 2 ]; then finish 22; fi" in sh and "if [ $rc = 5 ]; then finish 25; fi" in sh
    assert "if [ $rc = 3 ]; then finish 43; fi" in sh and "if [ $rc = 4 ]; then finish 44; fi" in sh
    assert "export NUMPY_MADVISE_HUGEPAGE=0" in sh and "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/" \
        "gata4_20261028T122425Z  # read only; fixed, no override" in sh
    # the same parts and the same steps as K22's launcher
    part = lambda s: s.split("PARTS=", 1)[1].split("rc=0\n$PY \"$SCRIPT\" read", 1)[0]  # noqa: E731
    assert part(sh).replace("K23", "K22") == part(k22sh)
    bsh = (T3 / "build_outside_form.sh").read_text()
    assert "XFLAG=(--table-exponent \"$EXP\")" in bsh and "2.0|4.0) ;;" in bsh
    assert "if [ \"$KIND1\" != xall ] || [ -n \"$KIND2\" ]" in bsh
    assert "gate s1 20" in bsh and "one \"$form\" \"$f\" ${XFLAG[@]+\"${XFLAG[@]}\"} || rc=$?" in bsh
