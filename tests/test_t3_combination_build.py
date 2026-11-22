"""Tests for experiments/t3_decomp/build_combination.py: the C2 factor of row R17 times the square of row K23's
outside table (the combination file CK of row R18), with its identity gates.

Offline, synthetic and fast. It reuses the throwaway-repository fixture and the synthetic WT file of
tests/test_t3_barred_exponent_build.py, synthetic `transfer_mab21l2_all` and `gata4_pooled` tables in the run's
factors/, a hand-made outside table `all` in combine.py's layout, and a listing, a K23 reading and a ledger committed in
a throwaway repository. No competition data and no outside data.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
T3 = ROOT / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


tx = sys.modules.get("_t3_barred_exponent_helpers") or _load("_t3_barred_exponent_helpers",
                                                             ROOT / "tests" / "test_t3_barred_exponent_build.py")
bc = _load("build_combination", T3 / "build_combination.py")
tk, bo, ext, bs = bc.tk, bc.bo, bc.ext, bc.bs
bx, t3d = tx.bx, bs.t3d
own_repo = tx.own_repo  # builds are stamped, and the tree checked, against a throwaway repository
SEED, N, G, GENES = tx.SEED, tx.N, tx.G, tx.GENES
PHI_M = tx.PHI  # up 2 (0-9), down 0.5 (10-14), the clip ends 4 (15-19) and 0.25 (20-24), 1 elsewhere
PHI_G = np.ones(G)
PHI_G[0:5], PHI_G[5:10], PHI_G[10:13] = 2.0, 0.5, 4.0
PHI_G[15:18], PHI_G[18:20] = 0.25, 4.0
PHI_G[20:23], PHI_G[23:25] = 4.0, 0.25
PHI_G[25:30], PHI_G[30:33] = 2.0, 0.25
# the outside table: 17 genes off 1, |corr(ln phi_G, ln T)| about 0.13, SD(ln CK) within 1% of the base's
T_ALL = np.ones(G)
T_ALL[0:10] = [0.9, 1.1, 0.95, 1.05, 0.9, 1.1, 0.95, 1.05, 0.9, 1.1]
T_ALL[33:40] = [1.2, 0.85, 1.1, 0.9, 1.15, 0.8, 1.05]
REPO_MAIN = "main"


# ------------------------------------------------------------------------------------------------ helpers


def write_tables(d: Path, tables: dict) -> dict[str, str]:
    """Hand-made tables in combine.py's layout under d: {name: (values, studies, tiers)}; returns their sha256."""
    d.mkdir(parents=True, exist_ok=True)
    comb = {}
    for name, (vals, studies, tiers) in tables.items():
        fd = d / name / "factors"
        fd.mkdir(parents=True, exist_ok=True)
        p = fd / f"{name}.csv"
        pd.DataFrame({"gene": GENES, "factor": vals}).to_csv(p, index=False)
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        (fd / "factors.json").write_text(json.dumps({"files": {name: {"file": str(p), "sha256": sha}}}))
        files = {name: {"file": str(p), "sha256": sha}}
        (d / name / f"{name}_manifest.json").write_text(json.dumps({
            "name": name, "table": name, "studies": studies, "tiers": tiers, "sign_rule": "applied",
            "refstats_sha256": "r" * 64, "files": files, "checked": "test"}))
        comb[name] = files
    (d / "COMBINE.json").write_text(json.dumps({"tables": comb, "record": {}}))
    return {n: v[n]["sha256"] for n, v in comb.items()}


def commit_file(repo: Path, rel: str, text_: str) -> Path:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text_)
    tx._git(repo, "add", rel)
    tx._git(repo, "commit", "-q", "--allow-empty", "-m", f"add {rel}")
    return p


def write_run(run: Path, wt: Path, tables: dict) -> None:
    fdir = run / "factors"
    fdir.mkdir(parents=True, exist_ok=True)
    files = {}
    for name, phi in tables.items():
        p = fdir / f"{name}.csv"
        pd.DataFrame({"gene": GENES, "factor": phi}).to_csv(p, index=False)
        files[name] = {"file": str(p), "sha256": t3d.sha256(p)}
    (fdir / "factors.json").write_text(json.dumps({"preset": "gata4", "files": files, "parameters": {},
                                                   "inputs": {"wt": {"path": "wt", "sha256": t3d.sha256(wt)}}}))


ROWS = ("| R17 (2026-11-20) two-knockout list | ... |\n"
        "| R18 (2026-11-21) the combination CK = C2 x T ** 2 | ... |\n")


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    from vec.boards import Board

    d = tmp_path_factory.mktemp("t3ck")
    wt = d / "wt.h5ad"
    tx._write_wt(wt, np.random.default_rng(0))
    run = d / "run"
    write_run(run, wt, {tk.TABLE_M: PHI_M, tk.TABLE_G: PHI_G})
    board = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, min_cells=200, max_cells=3000)
    dec = tx._decision(d / "gata4_decisions.json")
    pins = {k: t3d.sha256(run / "factors" / f"{k}.csv") for k in (tk.TABLE_M, tk.TABLE_G)}
    return SimpleNamespace(d=d, wt=wt, run=run, board=board, dec=dec, pins=pins)


@pytest.fixture
def lst(tmp_path_factory):
    repo = tmp_path_factory.mktemp("listing_repo")
    tx._git(repo, "init", "-q")
    return repo


def k23_reading(lst, names=("x_e2",), void=False, tables=None, rel="k23_reading.json") -> Path:
    pins = {n: v for n, v in ext.K22_TABLE_SHA256.items() if v}
    body = {"incomplete": False, "void": void, "names": list(names), "outcome": "test",
            "ext": {"listing_sha256": ext.K22_LISTING_SHA256, "tables": tables or pins},
            "pregate": {"constants": {"exponents": [2.0, 4.0]}, "dropped": {}}}
    return commit_file(lst, rel, json.dumps(body))


def setup(tmp_path, lst, monkeypatch, world, T=T_ALL, rows=ROWS, tag="") -> SimpleNamespace:
    """Synthetic pins of phi_M and phi_G; the table `all` (listed alone), its listing, the K23 reading naming x_e2
    and the ledger, all committed in `lst`; the module paths patched to them."""
    monkeypatch.setattr(tk, "TABLE_SHA256", dict(world.pins))
    monkeypatch.setattr(bo, "TABLE_SHA256", world.pins[tk.TABLE_M])
    d = tmp_path / f"tables{tag}"
    shas = write_tables(d, {"all": (T, ["GSE204908", "GSE277839"], ["unrelated"])})
    p = commit_file(lst, f"TABLES{tag}.sha256", f"{shas['all']}  all\n")
    monkeypatch.setattr(ext, "TABLES_FILE", p)
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", hashlib.sha256(p.read_bytes()).hexdigest())
    monkeypatch.setattr(ext, "K22_TABLE_SHA256", {"all": shas["all"], "unrelated": None})
    r = k23_reading(lst, rel=f"k23_reading{tag}.json")
    monkeypatch.setattr(bo, "READING_K23", r)
    tries = commit_file(lst, f"TRIES{tag}.md", rows)
    monkeypatch.setattr(bc, "TRIES_PATH", tries)
    return SimpleNamespace(env={ext.ENV: str(d)}, d=d, sha=shas["all"], reading=r, tries=tries)


@pytest.fixture
def cw(world, tmp_path, lst, monkeypatch):
    return setup(tmp_path, lst, monkeypatch, world)


def build(world, cw, form, out: Path, **kw) -> dict:
    return bc.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=kw.pop("run_dir", world.run),
                         seed=SEED, threads=2, env=kw.pop("env", cw.env), **kw)


class Spy:
    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


def read(path: Path):
    import anndata as ad

    return ad.read_h5ad(path)


# ------------------------------------------------------------------------------------------ the modules, unchanged


def test_the_existing_builders_are_imported_unchanged_and_share_the_frozen_builder():
    assert tk.bs is bs and bo.bs is bs and bo.ext is ext
    assert bs.__file__ == str(T3 / "build_submission.py") and tk.bs.load_factor.__name__ == "load_factor"
    assert bc.TABLE_EXPONENT == 2.0 and bc.SD_BOUND == 1.05 and bc.TABLE_T == "all"
    assert bc.BASES == tk.BASES == bo.BASES
    assert bc.KINDS == {"identity": ("identity", None), "gate-c2": ("c2", 0.0), "gate-k23": ("identity", 2.0),
                        "ck": ("c2", 2.0)}
    assert bc.FORMS == tuple(f"{k}-{b}" for k in ("identity", "gate-c2", "gate-k23", "ck")
                             for b in ("s1", "s1x1.5", "s1x2.0", "s2"))
    assert bc.LIST_ENTRY == "docs/PLAN.md decision log 2026-11-21 (b); experiments/harness/TRIES.md row R18"
    assert bc.TRIES_PATH == ROOT / "experiments" / "harness" / "TRIES.md"


# ---------------------------------------------------------------------------------------------- the factor by hand


def test_the_product_by_hand_for_one_gene():
    m = np.array([2.0, 0.5, 1.0, 4.0])
    g = np.array([2.0, 0.5, 1.0, 0.25])
    T = np.array([1.5, 1.0, 0.8, 1.0])
    # C2 at x = 1: ln m = a[1, -1, 0, 2], ln g = a[1, -1, 0, -2]; c = SD(ln m) / SD(ln m + ln g) = sqrt(10) / 4
    c = math.sqrt(10.0) / 4.0
    phi, rec = bc.ck_factor("c2", 1.0, m, g, T, 2.0)
    # gene 0: phi_M ** c * phi_G ** c * T ** 2 = 2 ** c * 2 ** c * 1.5 ** 2 = 4 ** c * 2.25
    assert phi[0] == pytest.approx(4.0 ** c * 2.25, rel=1e-14, abs=0)
    assert np.allclose(phi, [4.0 ** c * 2.25, 0.25 ** c, 0.64, 1.0], rtol=1e-13, atol=0)
    assert rec["part"]["c"] == pytest.approx(c, rel=1e-14, abs=0) and rec["table_exponent"] == 2.0
    # at x = 2 the same gene: 2 ** (2 c) * 2 ** c * 2.25, with c = 2 sqrt(5 / 21)
    c2 = 2.0 * math.sqrt(5.0 / 21.0)
    phi, _ = bc.ck_factor("c2", 2.0, m, g, T, 2.0)
    assert phi[0] == pytest.approx(2.0 ** (3 * c2) * 2.25, rel=1e-13, abs=0)


@pytest.mark.parametrize("x", [1.0, 1.5, 2.0])
def test_the_factor_is_the_parts_in_their_order(x):
    m, g, T = np.asarray(PHI_M, dtype=np.float64), PHI_G, T_ALL
    part = tk.combine("c2", x, m, g)[0]
    ck, _ = bc.ck_factor("c2", x, m, g, T, 2.0)
    assert np.array_equal(ck, part * T ** 2.0)
    assert np.array_equal(bc.ck_factor("c2", x, m, g, T, 0.0)[0], part)          # gate-c2: C2's factor exactly
    assert np.array_equal(bc.ck_factor("identity", x, m, g, T, 2.0)[0], ext.product(m, x, T, 2.0))  # gate-k23
    assert np.array_equal(bc.ck_factor("identity", x, m, g, None, None)[0], tk.combine("identity", x, m, g)[0])
    # no rescale and no clip after the product: values beyond [0.25, 4] stay
    big = np.ones(G)
    big[15:20] = 1.2
    out = bc.ck_factor("c2", x, m, g, big, 2.0)[0]
    assert np.array_equal(out, part * big ** 2.0)


def test_c_is_never_recomputed_with_the_table():
    m = np.asarray(PHI_M, dtype=np.float64)
    for x in (1.0, 1.5, 2.0):
        c = tk.combine("c2", x, m, PHI_G)[1]["c"]
        for T in (T_ALL, np.ones(G), T_ALL ** 3):
            _, rec = bc.ck_factor("c2", x, m, PHI_G, T, 2.0)
            assert rec["part"]["c"] == c
        lm, lg = np.log(m), np.log(PHI_G)
        assert c == pytest.approx(np.std(x * lm) / np.std(x * lm + lg), rel=1e-14)


def test_a_table_needs_its_exponent_and_the_genes_must_agree():
    m = np.asarray(PHI_M, dtype=np.float64)
    with pytest.raises(bc.Refused):
        bc.ck_factor("c2", 1.0, m, PHI_G, T_ALL, None)
    with pytest.raises(bc.Refused):
        bc.ck_factor("c2", 1.0, m, PHI_G, None, 2.0)
    with pytest.raises(bc.Refused, match="same genes"):
        bc.ck_factor("c2", 1.0, m, PHI_G, T_ALL[:-1], 2.0)


# ------------------------------------------------------------------------------------------- identity gates (bytes)


@pytest.mark.parametrize("base", ["s1", "s1x1.5", "s1x2.0", "s2"])
def test_identity_gives_the_frozen_barred_file_on_every_base(world, cw, tmp_path, base):
    slot, x, _ = bc.BASES[base]
    if x == 1.0:
        ref = bs.build_submission(world.dec, "barred", slot, world.wt, world.board, N, tmp_path / "ref.h5ad",
                                  run_dir=world.run, seed=SEED, threads=2)
    else:
        ref = bx.build_exponent(world.dec, f"{x}", world.wt, world.board, N, tmp_path / "ref.h5ad",
                                run_dir=world.run, seed=SEED, threads=2)
    side = build(world, cw, f"identity-{base}", tmp_path / "id.h5ad", env={})  # no table is read
    assert (tmp_path / "ref.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()
    assert side["output"]["sha256"] == ref["output"]["sha256"]
    cb = json.loads((tmp_path / "id.json").read_text())["combination_build"]
    assert cb["form"] == f"identity-{base}" and cb["table_exponent"] is None and cb["tables"]["all"] is None
    assert cb["k23_reading"] is None and cb["ledger_rows"] is None
    assert bs.load_factor.__module__ == "build_submission"  # left as it was


@pytest.mark.parametrize("base", ["s1", "s1x1.5", "s1x2.0", "s2"])
def test_gate_c2_gives_the_two_knockout_builders_c2_file(world, cw, tmp_path, base):
    ref = tk.build_form(world.dec, f"c2-{base}", world.wt, world.board, N, tmp_path / "c2.h5ad", run_dir=world.run,
                        seed=SEED, threads=2)
    side = build(world, cw, f"gate-c2-{base}", tmp_path / "g.h5ad")
    assert (tmp_path / "c2.h5ad").read_bytes() == (tmp_path / "g.h5ad").read_bytes()
    assert side["output"]["sha256"] == ref["output"]["sha256"]
    cb = side["combination_build"]
    assert cb["c"] == ref["combined_build"]["combination"]["c"] and cb["table_exponent"] == 0.0
    assert cb["parts_sha256"] == {"c2": bc.C2_SHA256[base], "k23": bc.K23_SHA256[base]}


@pytest.mark.parametrize("base", ["s1", "s1x1.5", "s1x2.0", "s2"])
def test_gate_k23_gives_the_outside_builders_file_at_exponent_2(world, cw, tmp_path, base):
    ref = bo.build_form(world.dec, f"xall-{base}", world.wt, world.board, N, tmp_path / "k.h5ad", run_dir=world.run,
                        seed=SEED, threads=2, env=cw.env, table_exponent=2.0, k23_reading=cw.reading)
    side = build(world, cw, f"gate-k23-{base}", tmp_path / "g.h5ad")
    assert (tmp_path / "k.h5ad").read_bytes() == (tmp_path / "g.h5ad").read_bytes()
    assert side["output"]["sha256"] == ref["output"]["sha256"]
    cb = side["combination_build"]
    assert cb["gata4_weight"] == 0.0 and cb["c"] == 1.0 and cb["table_exponent"] == 2.0
    assert cb["tables"]["all"]["sha256"] == cw.sha == ref["outside_build"]["table"]["sha256"]


# ---------------------------------------------------------------------------------------------------- CK built


@pytest.mark.parametrize("base", ["s1", "s1x1.5", "s1x2.0", "s2"])
def test_ck_applies_c2_times_the_square_once_and_records_it(world, cw, tmp_path, monkeypatch, base):
    spy = Spy(monkeypatch)
    side = build(world, cw, f"ck-{base}", tmp_path / "ck.h5ad")
    (phi,) = spy.seen
    x = bc.BASES[base][1]
    m = np.asarray(PHI_M, dtype=np.float64)
    part = tk.combine("c2", x, m, PHI_G)[0]
    assert np.array_equal(phi, part * T_ALL ** 2.0)
    assert side["construction"]["factor_exponent"] == 1.0 and side["spec"]["factor"] == bc.FROZEN_FACTOR
    assert side["inputs"]["factor_note"] == bc.FACTOR_NOTE
    on_disk = json.loads((tmp_path / "ck.json").read_text())
    cb = on_disk["combination_build"]
    assert cb["form"] == f"ck-{base}" and cb["x"] == x and cb["table_exponent"] == 2.0
    assert cb["c"] == tk.combine("c2", x, m, PHI_G)[1]["c"]
    assert cb["mab21l2_exponent"] == pytest.approx(x * cb["c"]) and cb["gata4_weight"] == cb["c"]
    assert set(cb["tables"]) == {tk.TABLE_M, tk.TABLE_G, "all"}
    assert cb["tables"]["all"]["sha256"] == cw.sha and cb["listing_sha256"] == ext.K22_LISTING_SHA256
    assert all(cb["tables"][k]["sha256"] == cb["tables"][k]["pinned_sha256"] for k in cb["tables"])
    assert cb["k23_reading"]["names"] == ["x_e2"] and cb["ledger_rows"]["rows"] == ["R17", "R18"]
    assert cb["parts_sha256"] == {"c2": bc.C2_SHA256[base], "k23": bc.K23_SHA256[base]}
    f = cb["figures"]
    lm = np.log(m)
    assert f["sd_ln_base"] == pytest.approx(np.std(x * lm), rel=1e-12)
    assert f["sd_ln_c2"] == pytest.approx(np.std(x * lm), rel=1e-12)
    assert f["sd_ln_factor"] == pytest.approx(np.std(np.log(phi)), rel=1e-12)
    assert f["genes_changed_against_c2"] == int((T_ALL != 1).sum()) == 17
    mask = T_ALL != 1
    assert f["corr_ln_t_ln_phi_g"] == pytest.approx(np.corrcoef(np.log(PHI_G[mask]), np.log(T_ALL[mask]))[0, 1])
    assert f["corr_ln_t_ln_phi_m"] == pytest.approx(np.corrcoef(lm[mask], np.log(T_ALL[mask]))[0, 1])
    assert f["factor"]["genes_at_clip"] >= 0 and f["c2"]["min"] == pytest.approx(part.min())
    assert cb["frozen_files_sha256"] == tk.FROZEN_SHA256 and cb["tree_at_start"]["clean"] is True
    assert cb["script_sha256"] == t3d.sha256(T3 / "build_combination.py") and cb["list_entry"] == bc.LIST_ENTRY
    assert set(on_disk) == set(side) and cb["depth_step"] is None

    def keys(o):
        if isinstance(o, dict):
            for k, v in o.items():
                yield k
                yield from keys(v)
    assert not any("cell" in k.lower() for k in keys(cb))
    assert bs.load_factor.__module__ == "build_submission"


@pytest.mark.parametrize("base", ["s1", "s2"])
def test_ck_draws_the_same_cells_as_its_base(world, cw, tmp_path, base):
    build(world, cw, f"identity-{base}", tmp_path / "id.h5ad")
    s = build(world, cw, f"ck-{base}", tmp_path / "ck.h5ad")
    a, b = read(tmp_path / "id.h5ad"), read(tmp_path / "ck.h5ad")
    assert list(a.obs_names) == list(b.obs_names) and a.obs.equals(b.obs) and a.shape == b.shape
    assert np.array_equal(a.obsm["spatial_3D"], b.obsm["spatial_3D"])
    assert not np.array_equal(np.asarray(a.X), np.asarray(b.X))
    assert s["inputs"]["rng_streams"] == [f"submission:barred:{base}", f"submission:barred:{base}:assemble"]
    assert s["inputs"]["seed"] == SEED and s["output"]["check"]["ok"]


def test_ck_twice_gives_the_same_bytes(world, cw, tmp_path):
    build(world, cw, "ck-s1x1.5", tmp_path / "a.h5ad")
    build(world, cw, "ck-s1x1.5", tmp_path / "b.h5ad")
    assert (tmp_path / "a.h5ad").read_bytes() == (tmp_path / "b.h5ad").read_bytes()


def test_the_existing_builders_default_outputs_are_unchanged(world, cw, tmp_path):
    """The two builders' own files are the same before and after a CK build in this process, and the frozen loader
    is restored."""
    original = bs.load_factor
    tk.build_form(world.dec, "c2-s1", world.wt, world.board, N, tmp_path / "c2a.h5ad", run_dir=world.run, seed=SEED,
                  threads=2)
    bo.build_form(world.dec, "identity-s1", world.wt, world.board, N, tmp_path / "ida.h5ad", run_dir=world.run,
                  seed=SEED, threads=2)
    build(world, cw, "ck-s1", tmp_path / "ck.h5ad")
    assert bs.load_factor is original
    tk.build_form(world.dec, "c2-s1", world.wt, world.board, N, tmp_path / "c2b.h5ad", run_dir=world.run, seed=SEED,
                  threads=2)
    bo.build_form(world.dec, "identity-s1", world.wt, world.board, N, tmp_path / "idb.h5ad", run_dir=world.run,
                  seed=SEED, threads=2)
    assert (tmp_path / "c2a.h5ad").read_bytes() == (tmp_path / "c2b.h5ad").read_bytes()
    assert (tmp_path / "ida.h5ad").read_bytes() == (tmp_path / "idb.h5ad").read_bytes()


# ------------------------------------------------------------------------------------------- refusals (exit 2)


@pytest.mark.parametrize("form", ["c1-s1", "c2-s1", "c2-s1x2.0", "xall-s1", "xunrel-s2", "x_e2-s1", "x_e4-s1",
                                  "x_e2_perm-s1", "x_e2_size-s1", "xall_perm-s1", "rep_all-s1", "barred_form-s1",
                                  "x_e1-s1", "ck", "ck-s3", "CK-s1", "ck-s1 ", " ck-s1", "ck-s1\n", "ck_t4-s1",
                                  "c1xT-s1", "gate-c2", "gate-ck-s1", "identity", "", None])
def test_other_forms_are_refused_by_name(world, cw, tmp_path, form):
    with pytest.raises(bc.Refused, match="form"):
        build(world, cw, form, tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


def test_the_other_builders_forms_point_to_their_builders():
    with pytest.raises(bc.Refused, match="build_two_knockout"):
        bc.parse_form("c2-s1")
    with pytest.raises(bc.Refused, match="build_outside_form"):
        bc.parse_form("xall-s1")
    with pytest.raises(bc.Refused, match="build_outside_form"):
        bc.parse_form("x_e2-s1")
    assert bc.parse_form("gate-c2-s1x1.5") == ("gate-c2", "s1x1.5") and bc.parse_form("ck-s2") == ("ck", "s2")


def test_an_unclean_tree_and_changed_frozen_files_are_refused(world, cw, tmp_path, monkeypatch, own_repo):
    for name, attr in (("t3_decompose.py", "T3D_PATH"), ("rules.yaml", "RULES_PATH"),
                       ("build_submission.py", "BUILDER_PATH")):
        changed = tmp_path / name
        shutil.copy(getattr(bs, attr), changed)
        with changed.open("a") as fh:
            fh.write("\n# changed\n")
        with monkeypatch.context() as m:
            m.setattr(bs, attr, changed)
            with pytest.raises(bc.Refused, match="frozen file"):
                build(world, cw, "ck-s1", tmp_path / f"{name}.h5ad")
        assert not (tmp_path / f"{name}.h5ad").exists()
    (own_repo / "stray.py").write_text("x = 1\n")
    with pytest.raises(bc.Refused, match="not clean"):
        build(world, cw, "ck-s1", tmp_path / "u.h5ad")
    assert not (tmp_path / "u.h5ad").exists()


def test_tables_off_their_pins_are_refused(world, cw, tmp_path, monkeypatch):
    for name in (tk.TABLE_M, tk.TABLE_G):
        with monkeypatch.context() as m:
            m.setattr(tk, "TABLE_SHA256", {**world.pins, name: "0" * 64})
            with pytest.raises(bc.Refused, match="pinned table"):
                build(world, cw, "ck-s1", tmp_path / f"{name}.h5ad")
        assert not (tmp_path / f"{name}.h5ad").exists()
    # T off its pin: the file changed after the listing
    p = cw.d / "all" / "factors" / "all.csv"
    good = p.read_bytes()
    p.write_bytes(good + b"\n")
    with pytest.raises(bc.Refused, match="not the listed"):
        build(world, cw, "ck-s1", tmp_path / "t.h5ad")
    p.write_bytes(good)
    # the table pin itself other than the listing's
    with monkeypatch.context() as m:
        m.setattr(ext, "K22_TABLE_SHA256", {"all": "0" * 64, "unrelated": None})
        with pytest.raises(bc.Refused):
            build(world, cw, "ck-s1", tmp_path / "pin.h5ad")
    with monkeypatch.context() as m:
        m.setattr(ext, "K22_LISTING_SHA256", None)
        with pytest.raises(bc.Refused, match="not filled"):
            build(world, cw, "gate-k23-s1", tmp_path / "pins.h5ad")
    with pytest.raises(bc.Refused, match="not set"):
        build(world, cw, "ck-s1", tmp_path / "env.h5ad", env={})
    assert not any(tmp_path.glob("*.h5ad"))


def test_the_k23_reading_must_name_x_e2_alone(world, cw, lst, tmp_path, monkeypatch):
    cases = ((k23_reading(lst, names=("x_e4",), rel="r1.json"), "does not name x_e2"),
             (k23_reading(lst, names=("x_e2", "x_e4"), rel="r2.json"), "does not name x_e2"),
             (k23_reading(lst, void=True, rel="r3.json"), "VOID"),
             (k23_reading(lst, tables={"all": "0" * 64}, rel="r4.json"), "other tables"),
             (tmp_path / "none.json", "no K23 reading"))
    for r, why in cases:
        for form in ("ck-s1", "gate-k23-s1", "gate-c2-s1"):
            with monkeypatch.context() as m:
                m.setattr(bo, "READING_K23", r)
                with pytest.raises(bc.Refused, match=why):
                    build(world, cw, form, tmp_path / "x.h5ad")
    loose = tmp_path / "loose.json"
    loose.write_bytes(cw.reading.read_bytes())
    monkeypatch.setattr(bo, "READING_K23", loose)
    with pytest.raises(bc.Refused, match="not committed"):
        build(world, cw, "ck-s1", tmp_path / "x.h5ad")
    # identity reads no table and no reading
    build(world, cw, "identity-s1", tmp_path / "id.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


def test_the_ledger_rows(world, cw, lst, tmp_path, monkeypatch):
    no18 = commit_file(lst, "no18.md", "| R17 (2026-11-20) the list |\n| R16 (x) CK |\n")
    no_ck = commit_file(lst, "nock.md", "| R17 (2026-11-20) the list |\n| R18 (x) something else |\n")
    no17 = commit_file(lst, "no17.md", "| R18 (x) CK |\n")
    indented = commit_file(lst, "indented.md", " | R17 (a) |\n | R18 (b) CK |\n")
    for path, form, why in ((no18, "ck-s1", "no row R18"), (no_ck, "ck-s1", "no row R18"),
                            (no17, "ck-s1", "no row R17"), (no17, "gate-c2-s1", "no row R17"),
                            (indented, "ck-s1", "no row R17"), (tmp_path / "none.md", "ck-s1", "no ledger")):
        monkeypatch.setattr(bc, "TRIES_PATH", path)
        with pytest.raises(bc.Refused, match=why):
            build(world, cw, form, tmp_path / "x.h5ad")
    loose = tmp_path / "loose.md"
    loose.write_text(ROWS)
    monkeypatch.setattr(bc, "TRIES_PATH", loose)
    with pytest.raises(bc.Refused, match="not committed"):
        build(world, cw, "ck-s1", tmp_path / "x.h5ad")
    # gate-c2 needs R17 only; gate-k23 and identity read no ledger
    monkeypatch.setattr(bc, "TRIES_PATH", no18)
    build(world, cw, "gate-c2-s1", tmp_path / "g.h5ad")
    monkeypatch.setattr(bc, "TRIES_PATH", tmp_path / "none.md")
    build(world, cw, "gate-k23-s1", tmp_path / "k.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


def test_frozen_builder_refusals_still_apply(world, cw, tmp_path):
    out = tmp_path / "x.h5ad"
    build(world, cw, "ck-s2", out)
    with pytest.raises(bc.Refused, match="exists"):
        build(world, cw, "ck-s2", out)
    with pytest.raises(bc.Refused, match="seed"):
        bc.build_form(world.dec, "ck-s1", world.wt, world.board, N, tmp_path / "s.h5ad", run_dir=world.run, seed=0,
                      threads=2, env=cw.env)
    dec = tx._decision(tmp_path / "dec.json", script_sha256="0" * 64)
    with pytest.raises(bc.Refused, match="script"):
        bc.build_form(dec, "ck-s1", world.wt, world.board, N, tmp_path / "d.h5ad", run_dir=world.run, seed=SEED,
                      threads=2, env=cw.env)


# ------------------------------------------------------------------------------------- the form's own refusals (4)


def _sd_table():
    T = np.ones(G)
    T[0:10], T[10:15] = 1.6, 0.6  # SD(ln CK) far above 1.05 SD(x ln phi_M); |corr| with ln phi_G about 0.495
    return T


def _corr_table():
    T = np.ones(G)
    T[0:5], T[5:10], T[25:30], T[30:33] = 1.1, 0.9, 1.1, 0.9  # ln T follows ln phi_G (corr about 0.97)
    return T


def _few_table():
    T = np.ones(G)
    T[[0, 5, 33]] = 1.2, 1.2, 0.8  # three genes; corr with ln phi_G 0
    return T


@pytest.mark.parametrize("T,why", [(_sd_table(), "SD"), (_corr_table(), "corr"), (_few_table(), "3 genes")])
def test_the_form_s_own_refusals(world, tmp_path, lst, monkeypatch, T, why):
    setup(tmp_path, lst, monkeypatch, world, T=T, tag="own")
    with pytest.raises(bc.FormRefused, match=why):
        bc.build_form(world.dec, "ck-s1", world.wt, world.board, N, tmp_path / "x.h5ad", run_dir=world.run,
                      seed=SEED, threads=2, env={ext.ENV: str(tmp_path / "tablesown")})
    assert not (tmp_path / "x.h5ad").exists() and not (tmp_path / "x.json").exists()
    # the gates are not the form: they build on the same table
    bc.build_form(world.dec, "gate-k23-s1", world.wt, world.board, N, tmp_path / "g.h5ad", run_dir=world.run,
                  seed=SEED, threads=2, env={ext.ENV: str(tmp_path / "tablesown")})


def test_the_own_refusals_by_hand():
    m = np.asarray(PHI_M, dtype=np.float64)
    part = tk.combine("c2", 1.0, m, PHI_G)[0]
    phi = part * T_ALL ** 2.0
    bc.own_refusals(1.0, m, PHI_G, T_ALL, part, phi)  # passes
    bad = phi.copy()
    bad[3] = np.inf
    with pytest.raises(bc.FormRefused, match="finite and positive"):
        bc.own_refusals(1.0, m, PHI_G, T_ALL, part, bad)
    bad[3] = 0.0
    with pytest.raises(bc.FormRefused, match="finite and positive"):
        bc.own_refusals(1.0, m, PHI_G, T_ALL, part, bad)
    flat = np.ones(G)
    flat[33:40] = [1.1, 0.9, 1.05, 0.95, 1.1, 0.9, 1.02]  # ln phi_G is 0 wherever T != 1: the correlation is undefined
    with pytest.raises(bc.FormRefused, match="undefined"):
        bc.own_refusals(1.0, m, PHI_G, flat, part, part * flat ** 2.0)
    # the SD bound at exactly the edge: a factor whose SD(ln) is 1.05 SD(ln phi_M) passes, above it is refused
    lm = np.log(m)
    sd = np.std(lm)
    for scale, ok in ((1.0499, True), (1.0501, False)):
        f = np.exp(lm * scale)
        assert np.std(np.log(f)) == pytest.approx(scale * sd, rel=1e-12)
        if ok:
            bc.own_refusals(1.0, m, PHI_G, T_ALL, part, f)
        else:
            with pytest.raises(bc.FormRefused, match="SD"):
                bc.own_refusals(1.0, m, PHI_G, T_ALL, part, f)
    assert ext.MAX_ABS_CORR == 0.5 and ext.MIN_CHANGED == 10


# --------------------------------------------------------------------------------------------------------------- CLI


def _boards(tmp_path: Path) -> Path:
    boards = tmp_path / "boards"
    boards.mkdir(exist_ok=True)
    (boards / "T3__synthetic.genes.txt").write_text("\n".join(GENES) + "\n")
    (boards / "index.json").write_text(json.dumps({"T3:synthetic": {
        "task": "T3", "genes_file": "T3__synthetic.genes.txt", "needs_coords": True, "min_cells": 200,
        "max_cells": N}}))
    return boards


def test_cli_takes_names_only_and_its_exit_codes(world, cw, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    from vec.io import check_submission

    monkeypatch.setenv(ext.ENV, cw.env[ext.ENV])
    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(_boards(tmp_path)), "--run-dir", str(world.run), "--threads", "2"]
    runner = CliRunner()
    res = runner.invoke(bc.app, [*base, "--form", "c2-s1", "--out", str(tmp_path / "a.h5ad")])
    assert res.exit_code == 2 and "refused" in res.output and not (tmp_path / "a.h5ad").exists()
    res = runner.invoke(bc.app, [*base, "--form", "ck-s1", "--out", str(tmp_path / "b.h5ad")])
    assert res.exit_code == 0, res.output
    assert check_submission(tmp_path / "b.h5ad", world.board).ok
    side = json.loads((tmp_path / "b.json").read_text())
    assert side["combination_build"]["form"] == "ck-s1" and side["branch"] == "barred" and side["submission"] == "s1"
    for extra in (["--table-exponent", "2.0"], ["--exponent", "1.5"], ["--c", "1"], ["--gata4-weight", "0.5"],
                  ["--table", "all"], ["--tables-dir", str(cw.d)], ["--branch", "allowed"], ["--depth"],
                  ["--seed", "1000002"], ["--n-cells", str(N)], ["--n-cells", "MAX"]):
        out = tmp_path / f"{extra[0].strip('-')}.h5ad"
        res = runner.invoke(bc.app, [*base, "--form", "ck-s1", *extra, "--out", str(out)])
        assert res.exit_code != 0 and not out.exists(), extra
    res = runner.invoke(bc.app, [*base, "--form", "ck-s1", "--seed", "1000002", "--out", str(tmp_path / "s.h5ad")])
    assert res.exit_code == 2 and "the list fixes seed 1000001 and max_cells" in res.output


def test_cli_exit_code_of_the_form_s_own_refusal(world, tmp_path, lst, monkeypatch):
    from typer.testing import CliRunner

    s = setup(tmp_path, lst, monkeypatch, world, T=_few_table(), tag="cli")
    monkeypatch.setenv(ext.ENV, s.env[ext.ENV])
    res = CliRunner().invoke(bc.app, ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic",
                                      "--boards-dir", str(_boards(tmp_path)), "--run-dir", str(world.run),
                                      "--threads", "2", "--form", "ck-s1", "--out", str(tmp_path / "c.h5ad")])
    assert res.exit_code == 4 and "own refusal" in res.output and not (tmp_path / "c.h5ad").exists()


# ------------------------------------------------------------------------------------------ the launcher's constants


def _sums(name: str) -> dict[str, str]:
    out = {}
    for ln in (ROOT / "experiments" / "builds" / "results" / name).read_text().splitlines():
        sha, f = ln.split("  ")
        out[f] = sha
    return out


def _array(sh: str, name: str) -> dict[str, str]:
    m = re.search(rf"declare -A {name}=\(\n(.*?)\n\)", sh, re.S)
    assert m, name
    return dict(re.fullmatch(r"  \[([^\]]+)\]=([0-9a-f]{64})", ln).groups() for ln in m.group(1).split("\n"))


def test_the_launchers_constants():
    sh = (T3 / "build_combination.sh").read_text()
    c2, k23 = _sums("t3two_fb01487.SHA256SUMS"), _sums("t3k23_49f3fd5.SHA256SUMS")
    assert _array(sh, "REF") == {b: v[2] for b, v in bc.BASES.items()}
    assert _array(sh, "C2REF") == {b: c2[f"c2-{b}.h5ad"] for b in bc.BASES} == bc.C2_SHA256
    assert _array(sh, "K23REF") == {b: k23[f"xall-{b}.h5ad"] for b in bc.BASES} == bc.K23_SHA256
    assert "TABLES=$HOME/v/vec-ext/t3/tables" in sh and "VEC_T3_EXT_TABLES=\"$TABLES\"" in sh
    assert "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z" in sh
    assert "S=~/v/scratch/t3ck_${HEAD_SHA:0:7}_" in sh and "O=runs/t3_combination_${HEAD_SHA:0:7}_$BASE" in sh
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and "TMPDIR=\"$S/tmp\"" in sh and "--seed 1000001" in sh
    assert "--table-exponent" not in sh
    for t in ("tests/test_t3_combination_build.py", "tests/test_t3_two_knockout_build.py",
              "tests/test_t3_outside_forms.py", "tests/test_k23_outside_exponents.py",
              "tests/test_k22_outside_forms.py", "tests/test_t3_barred_exponent_build.py",
              "experiments/t3_decomp/test_build_submission.py"):
        assert t in sh
    for code in ("finish 2", "finish 97", "finish 98", "finish 96", "finish 95", "finish 10", "finish 46",
                 "finish 45", "finish 40", "finish 41", "finish 42", "finish 43", "finish 44"):
        assert code in sh, code
    assert 'gate identity-s1 "${REF[s1]}" 20' in sh and 'gate "identity-$BASE" "${REF[$BASE]}" 30' in sh
    assert 'gate "gate-c2-$BASE" "${C2REF[$BASE]}" 50' in sh and 'gate "gate-k23-$BASE" "${K23REF[$BASE]}" 60' in sh
    assert sh.count("> /dev/null 2>&1") >= 3


def test_the_ck_hashes_cannot_equal_a_part_s():
    """The launcher compares CK only with itself; the parts' hashes are pinned so the record can show CK differs."""
    assert len(set(bc.C2_SHA256.values()) | set(bc.K23_SHA256.values())) == 8


# ------------------------------------------------------------------------------------------ byte-identity pins


PINNED = {
    "experiments/t3_decomp/rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13",
    "experiments/t3_decomp/t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
    "experiments/t3_decomp/build_submission.py": "1194a63993a27dde4dedfdb9180de035c5f8f2a57aaa7d97b38e7a15ea9e0c74",
    "experiments/t3_decomp/outside_forms.py": "4837687d47f134ff690177ecaaec7c712952c238d0e7c806c5c8ef2cecb3ceb5",
    "experiments/t3_decomp/build_two_knockout.py": "e6bff930c32d3377bd1b5dc8174c2d60814fa0867d160ab64f1d577f94a50c9b",
    "experiments/t3_decomp/build_two_knockout.sh": "dc1f8cff8fc2eb9cfa8c6bf9b85f3414e5639f8dbf1ca4984fd56ed6ffd69add",
    # build_outside_form.py and .sh: no longer pinned here; they gained row R22's T3 line (--captain-ruling R22, off by
    # default), and their unchanged paths are held by tests/test_t3_outside_forms.py and
    # tests/test_t3_exponent_line_build.py on the synthetic fixture
    "experiments/t3_decomp/k23_outside_exponents.py":
        "98bdac5b9f2413bc6d11e470fdbd2202df5ab3dc22548c8940812ac1f455eb98",
    "experiments/builds/results/t3two_fb01487.SHA256SUMS":
        "cad9d132221f811c3ae686311d827c46e75e385d86f26fda4830b464cc72aa64",
    "experiments/builds/results/t3two_fb01487.txt": "caf79cedb808583ccec11c8badd6d7229800fec4fe0e0efb3f7294a57ca08d67",
    "experiments/builds/results/t3k23_49f3fd5.SHA256SUMS":
        "a8333bed33e670159a592fbece1d6c96a940720cccc4f3f13f643dd97d8f3673",
    "experiments/builds/results/t3k23_49f3fd5.txt": "417837d73a066ea54fcad0045fbcbd17a42c6e56b9b945291d6cf5391c9caece",
}


def test_the_existing_files_are_byte_identical():
    for rel, sha in PINNED.items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == sha, rel
    for name in ("t3_decompose.py", "rules.yaml", "build_submission.py"):
        assert PINNED[f"experiments/t3_decomp/{name}"] == tk.FROZEN_SHA256[name] == ext.FROZEN_SHA256[name]


def test_the_existing_files_are_mains_bytes():
    env = {k: v for k, v in os.environ.items() if k not in tx.GIT_VARS}
    env["GIT_OPTIONAL_LOCKS"] = "0"
    for rel, sha in PINNED.items():
        res = subprocess.run(["git", "-C", str(ROOT), "show", f"{REPO_MAIN}:{rel}"], capture_output=True, env=env)
        if res.returncode != 0:
            pytest.skip("branch main is not available in this checkout")
        assert hashlib.sha256(res.stdout).hexdigest() == sha, rel
