"""Tests for experiments/t3_decomp/outside_forms.py (the outside tables by name, the product, the permuted control, the
pre-gate) and experiments/t3_decomp/build_outside_form.py (the thin builder). Offline, synthetic and fast: hand-made
tables in a temporary folder outside the repository, a temporary git repository for the listing and the reading, and
the helpers and throwaway-repository fixture of tests/test_t3_barred_exponent_build.py. No competition data and no
outside data."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
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
bo = sys.modules.get("build_outside_form") or _load("build_outside_form", T3 / "build_outside_form.py")
bk = sys.modules.get("build_two_knockout") or _load("build_two_knockout", T3 / "build_two_knockout.py")
bf = sys.modules.get("build_factor_form") or _load("build_factor_form", T3 / "build_factor_form.py")
import factor_forms as ff  # noqa: E402
import k7_wt_revert as k7  # noqa: E402
import outside_forms as ext  # noqa: E402

from vec import harness  # noqa: E402

bx, bs = tx.bx, bo.bs
t3d = bs.t3d
own_repo = tx.own_repo  # builds are stamped against a throwaway repository
SEED, N, G, GENES, PHI = tx.SEED, tx.N, tx.G, tx.GENES, tx.PHI


# ------------------------------------------------------------------------------------------------ helpers


def table_values(seed: int, n_moved: int = 20, genes: int = G) -> np.ndarray:
    """A synthetic finished table: n_moved genes off 1 in [0.4, 2.5], the rest exactly 1 (absent genes)."""
    rng = np.random.default_rng(seed)
    T = np.ones(genes)
    idx = rng.choice(genes, n_moved, replace=False)
    T[idx] = np.exp(rng.uniform(-0.9, 0.9, n_moved))
    return T


def write_tables(d: Path, genes: list[str], tables: dict) -> dict[str, str]:
    """Hand-made tables in combine.py's layout under d: tables = {name: (values, studies, tiers)}; returns sha256s."""
    d.mkdir(parents=True, exist_ok=True)
    comb = {}
    for name, (vals, studies, tiers) in tables.items():
        fd = d / name / "factors"
        fd.mkdir(parents=True, exist_ok=True)
        p = fd / f"{name}.csv"
        pd.DataFrame({"gene": list(genes), "factor": vals}).to_csv(p, index=False)
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        (fd / "factors.json").write_text(json.dumps({"files": {name: {"file": str(p), "sha256": sha}}}))
        pg = d / name / f"{name}_per_group.csv"
        pg.write_text("group,gene,factor\n")
        files = {name: {"file": str(p), "sha256": sha},
                 f"{name}_per_group": {"file": str(pg), "sha256": hashlib.sha256(pg.read_bytes()).hexdigest()}}
        (d / name / f"{name}_manifest.json").write_text(json.dumps({
            "name": name, "table": name, "studies": studies, "tiers": tiers,
            "sign_rule": "applied" if len(studies) >= 2 else "single study", "refstats_sha256": "r" * 64,
            "files": files, "checked": "by try K21 only, after its row"}))
        comb[name] = files
    (d / "COMBINE.json").write_text(json.dumps({"tables": comb, "record": {}}))
    return {n: v[n]["sha256"] for n, v in comb.items()}


def listing_text(shas: dict[str, str]) -> str:
    return "".join(f"{shas[n]}  {n}\n" for n in sorted(shas))


def commit_file(repo: Path, rel: str, text_: str) -> Path:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text_)
    tx._git(repo, "add", rel)
    tx._git(repo, "commit", "-q", "--allow-empty", "-m", f"add {rel}")
    return p


@pytest.fixture
def lst(tmp_path_factory):
    """A throwaway git repository holding a committed listing."""
    repo = tmp_path_factory.mktemp("listing_repo")
    tx._git(repo, "init", "-q")
    return repo


TABLES = {"all": (table_values(1), ["GSE137337", "GSE167493"], ["cardiac", "unrelated"]),
          "unrelated": (table_values(2), ["GSE137337", "GSE204908"], ["unrelated"]),
          "cardiac": (table_values(3), ["GSE167493"], ["cardiac"]),
          "study_GSE167493": (table_values(3), ["GSE167493"], ["cardiac"])}


def setup_tables(tmp_path: Path, repo: Path, monkeypatch, tables=None, listed=None, genes=None):
    """Tables in tmp_path/tables, their listing committed in repo, the pins and TABLES_FILE patched; returns env."""
    tables = TABLES if tables is None else tables
    d = tmp_path / "tables"
    shas = write_tables(d, genes or GENES, tables)
    listed = sorted(shas) if listed is None else listed
    p = commit_file(repo, "TABLES.sha256", listing_text({n: shas[n] for n in listed}))
    monkeypatch.setattr(ext, "TABLES_FILE", p)
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", hashlib.sha256(p.read_bytes()).hexdigest())
    monkeypatch.setattr(ext, "K22_TABLE_SHA256", {n: shas[n] if n in listed else None for n in ext.TRIES})
    return SimpleNamespace(env={ext.ENV: str(d)}, d=d, shas=shas, listing=p)


# ------------------------------------------------------------------------------------------------ constants


def test_constants_are_the_record_s():
    assert ext.TRIES == ("all", "unrelated") and ext.EXPONENT == 1.0 and ext.ENV == "VEC_T3_EXT_TABLES"
    assert ext.PERM_SEED == k7.PERM_SEED == 20261113 and ext.PERM_BINS == harness.K2_PERM_BINS == 10
    assert ext.MIN_CHANGED == ff.MIN_CHANGED == 10 and ext.CLIP == ff.CLIP
    assert (ext.MIN_STUDIES, ext.MAX_ABS_CORR, ext.FILE2_MAX_CORR) == (2, 0.5, 0.9)
    assert ext.TABLES_FILE == ROOT / "experiments" / "t3_ext" / "TABLES.sha256"


def test_frozen_hashes_are_pinned_and_the_files_match():
    assert ext.FROZEN_SHA256 == bk.FROZEN_SHA256 == bf.FROZEN_SHA256 == bo.FROZEN_SHA256
    for name, sha in ext.FROZEN_SHA256.items():
        assert hashlib.sha256((T3 / name).read_bytes()).hexdigest() == sha
    assert bo.BASES == bk.BASES and bo.TABLE_SHA256 == bk.TABLE_SHA256[bo.TABLE_M]


def test_the_pins_are_unset_or_the_committed_listing_s():
    """Until TABLES.sha256 is committed the pins are None; once it is, they must be its hash and lines."""
    if not ext.TABLES_FILE.exists():
        assert ext.K22_LISTING_SHA256 is None and ext.K22_TABLE_SHA256 == {"all": None, "unrelated": None}
        return
    raw = ext.TABLES_FILE.read_bytes()
    lines = ext.parse_listing(raw)
    assert ext.K22_LISTING_SHA256 == hashlib.sha256(raw).hexdigest()
    assert ext.K22_TABLE_SHA256 == {n: lines.get(n) for n in ext.TRIES}
    assert ext.tracked(ext.TABLES_FILE)


# ------------------------------------------------------------------------------------------------ the forms


def test_the_product_by_hand_with_a_gene_absent_from_the_table():
    phi_m = np.array([2.0, 0.5, 4.0, 1.0, 0.25])
    T = np.array([1.5, 1.0, 0.5, 2.0, 4.0])  # gene 1 absent from the outside data: exactly 1
    assert np.array_equal(ext.product(phi_m, 1.0, T), np.array([3.0, 0.5, 2.0, 2.0, 1.0]))
    assert np.allclose(ext.product(phi_m, 2.0, T), [6.0, 0.25, 8.0, 2.0, 0.25], rtol=0, atol=1e-15)
    assert np.allclose(ext.product(phi_m, 1.5, T), [2 ** 1.5 * 1.5, 0.5 ** 1.5, 8.0 * 0.5, 2.0, 0.125 * 4], atol=1e-15)
    assert ext.product(phi_m, 1.0, T)[1] == phi_m[1]  # an absent gene leaves phi_M as it is
    # nothing rescaled or clipped after the product: 4 * 4 = 16 stays 16
    assert ext.product(np.array([4.0, 1.0]), 1.0, np.array([4.0, 1.0]))[0] == 16.0
    for bad in (np.array([1.0, 0.0, 1, 1, 1]), np.array([1.0, np.nan, 1, 1, 1]), np.ones(4)):
        with pytest.raises(ext.FormRefused):
            ext.product(phi_m, 1.0, bad)


def test_the_permuted_control_uses_the_fixed_seed_and_keeps_each_decile():
    rng = np.random.default_rng(3)
    T, basis = table_values(4, 30, 60), rng.normal(size=60)
    p = ext.permuted(T, basis)
    assert np.array_equal(p, harness.t3_decile_permutation(T, basis, 20261113, 10))
    assert np.array_equal(p, ext.permuted(T, basis))  # one permutation, whatever the call
    assert not np.array_equal(p, harness.t3_decile_permutation(T, basis, 20261114, 10))
    assert not np.array_equal(p, T) and sorted(p) == sorted(T) and (p != 1).sum() == (T != 1).sum()
    for genes in np.array_split(np.argsort(basis, kind="stable"), 10):
        assert sorted(p[genes]) == sorted(T[genes])


# ------------------------------------------------------------------------------------------------ the listing


def test_the_listing_form(tmp_path):
    a, b = "a" * 64, "b" * 64
    assert ext.parse_listing(f"{a}  all\n{b}  unrelated\n".encode()) == {"all": a, "unrelated": b}
    for bad in (f"{a}  all", f"{a} all\n", f"{a.upper()}  all\n", f"{a}  all\r\n", f"{a}  all\n\n",
                f"# c\n{a}  all\n", f"{b}  unrelated\n{a}  all\n", f"{a}  all\n{a}  all\n",
                f"{a}  all_per_group\n", f"{a}  ../x\n", f"{a[:-1]}  all\n"):
        with pytest.raises(ext.TableRefused):
            ext.parse_listing(bad.encode())


def test_listing_refusals(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    listing, sha = ext.read_listing()
    assert listing == s.shas and sha == ext.K22_LISTING_SHA256
    with pytest.raises(ext.TableRefused, match="not filled"):
        ext.read_listing(pin=None)
    with pytest.raises(ext.TableRefused, match="pinned K22_LISTING_SHA256"):
        ext.read_listing(pin="0" * 64)
    with pytest.raises(ext.TableRefused, match="the pin"):
        ext.read_listing(pins={**ext.K22_TABLE_SHA256, "all": "0" * 64})
    with pytest.raises(ext.TableRefused, match="the pin"):  # listed but not pinned
        ext.read_listing(pins={**ext.K22_TABLE_SHA256, "unrelated": None})
    untracked = lst / "OTHER.sha256"
    untracked.write_bytes(s.listing.read_bytes())
    with pytest.raises(ext.TableRefused, match="not tracked"):
        ext.read_listing(untracked)
    link = lst / "LINK.sha256"
    link.symlink_to(s.listing)
    with pytest.raises(ext.TableRefused, match="symlink"):
        ext.read_listing(link)


def test_listing_and_pins_from_combine(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    assert ext.listing_from_combine(s.env) == s.listing.read_text()
    out = ext.pins_text(s.listing)
    assert f'K22_LISTING_SHA256: str | None = "{ext.K22_LISTING_SHA256}"' in out
    assert f'"all": "{s.shas["all"]}", "unrelated": "{s.shas["unrelated"]}"' in out
    p = s.d / "cardiac" / "factors" / "cardiac.csv"
    p.write_bytes(p.read_bytes() + b"\n")
    with pytest.raises(ext.TableRefused, match="COMBINE.json says"):
        ext.listing_from_combine(s.env)


# ------------------------------------------------------------------------------------------------ the loader


def _load_all(s, name="all", genes=None):
    listing, sha = ext.read_listing()
    return ext.load_table(name, genes or GENES, listing, sha, env=s.env)


def test_a_table_loads_by_name_with_its_record(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    T, rec = _load_all(s)
    assert np.array_equal(T, TABLES["all"][0])
    assert rec["sha256"] == s.shas["all"] == rec["pinned_sha256"] and rec["listing_sha256"] == ext.K22_LISTING_SHA256
    assert rec["studies"] == ["GSE137337", "GSE167493"] and rec["tiers"] == ["cardiac", "unrelated"]
    assert rec["file"] == "<VEC_T3_EXT_TABLES>/all/factors/all.csv" and rec["dir_realpath"] == str(s.d.resolve())


def test_loader_refusals_names_and_folder(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    listing, sha = ext.read_listing()
    for name, why in (("cardiac", "not one of K22's tables"), ("all_per_group", "not a pooled table name"),
                      ("../all", "not a pooled table name"), ("study_GSE167493", "not one of K22's tables")):
        with pytest.raises(ext.TableRefused, match=why):
            ext.load_table(name, GENES, listing, sha, env=s.env)
    with pytest.raises(ext.TableRefused, match="not set"):
        ext.load_table("all", GENES, listing, sha, env={})
    with pytest.raises(ext.TableRefused, match="absolute"):
        ext.load_table("all", GENES, listing, sha, env={ext.ENV: "tables"})
    link = tmp_path / "linked"
    link.symlink_to(s.d)
    with pytest.raises(ext.TableRefused, match="symlink"):
        ext.load_table("all", GENES, listing, sha, env={ext.ENV: str(link)})
    with pytest.raises(ext.TableRefused, match="not name a directory"):
        ext.load_table("all", GENES, listing, sha, env={ext.ENV: str(tmp_path / "none")})
    with pytest.raises(ext.TableRefused, match="inside the repository"):
        ext.load_table("all", GENES, listing, sha, env={ext.ENV: str(ROOT / "experiments")})
    with pytest.raises(ext.TableRefused, match="not listed"):
        ext.load_table("unrelated", GENES, {"all": listing["all"]}, sha, env=s.env)


def test_loader_refusals_bytes_and_records(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    p = s.d / "all" / "factors" / "all.csv"
    good = p.read_bytes()
    p.write_bytes(good + b"\n")
    with pytest.raises(ext.TableRefused, match="not the listed"):
        _load_all(s)
    p.write_bytes(good)
    for f, mutate, why in (
            (s.d / "COMBINE.json", lambda j: j["tables"]["all"]["all"].update(sha256="0" * 64), "COMBINE.json"),
            (s.d / "all" / "factors" / "factors.json", lambda j: j["files"]["all"].update(sha256="0" * 64),
             "factors.json"),
            (s.d / "all" / "all_manifest.json", lambda j: j["files"]["all"].update(sha256="0" * 64), "manifest")):
        raw = f.read_text()
        j = json.loads(raw)
        mutate(j)
        f.write_text(json.dumps(j))
        with pytest.raises(ext.TableRefused, match=why):
            _load_all(s)
        f.unlink()
        with pytest.raises(ext.TableRefused, match="missing or a symlink"):
            _load_all(s)
        f.write_text(raw)
    _load_all(s)
    # a symlinked csv, factors folder or table folder
    real = tmp_path / "real.csv"
    real.write_bytes(good)
    p.unlink()
    p.symlink_to(real)
    with pytest.raises(ext.TableRefused, match="symlink"):
        _load_all(s)
    p.unlink()
    p.write_bytes(good)


def test_loader_refusals_genes_and_values(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    with pytest.raises(ext.TableRefused, match="order included"):
        _load_all(s, genes=GENES[::-1])
    for vals, genes, why in ((np.r_[0.0, np.ones(G - 1)], GENES, "finite and positive"),
                             (np.r_[np.inf, np.ones(G - 1)], GENES, "finite and positive"),
                             (np.ones(G), [*GENES[:-1], GENES[0]], "duplicated")):
        d = tmp_path / f"t{why[:3]}{len(set(genes))}{vals[0]}"
        r = tmp_path / f"repo_{d.name}"
        r.mkdir()
        tx._git(r, "init", "-q")
        s = setup_tables(d, r, monkeypatch, tables={**TABLES, "all": (vals, ["a", "b"], ["unrelated"])}, genes=genes)
        with pytest.raises(ext.TableRefused, match=why):
            _load_all(s, genes=genes)
    # columns other than gene, factor
    s = setup_tables(tmp_path / "cols", lst, monkeypatch)
    p = s.d / "all" / "factors" / "all.csv"
    p.write_text(p.read_text().replace("gene,factor", "gene,phi", 1))
    sha = hashlib.sha256(p.read_bytes()).hexdigest()
    for f in (s.d / "COMBINE.json", s.d / "all" / "factors" / "factors.json", s.d / "all" / "all_manifest.json"):
        f.write_text(f.read_text().replace(s.shas["all"], sha))
    shas = {**s.shas, "all": sha}
    lp = commit_file(lst, "TABLES.sha256", listing_text(shas))
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", hashlib.sha256(lp.read_bytes()).hexdigest())
    monkeypatch.setattr(ext, "K22_TABLE_SHA256", {n: shas[n] for n in ext.TRIES})
    with pytest.raises(ext.TableRefused, match="columns"):
        _load_all(s)


# ------------------------------------------------------------------------------------------------ the pre-gate


def _phi_m(seed=9) -> np.ndarray:
    return np.clip(np.exp(np.random.default_rng(seed).normal(0, 0.6, G)), 0.25, 4.0)


def test_pregate_figures_and_no_drop_for_two_studies(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    loaded = ext.load_listed(GENES, env=s.env)
    phi_m = _phi_m()
    pre = ext.pregate(phi_m, loaded, basis=np.arange(G, dtype=float), den=np.full(G, 1.0 / G))
    for n in ext.TRIES:
        c = pre["tables"][n]["corr_ln_phi_m"]
        assert c is not None and abs(c) < 0.5, c
    assert pre["dropped"] == {}
    t = pre["tables"]["all"]
    assert t["genes_ne_1"] == 20 and t["studies"] == ["GSE137337", "GSE167493"] and "den_total" in t
    assert pre["products"]["all"]["genes_changed"] == 20 and pre["controls"]["all"]["genes_changed"] == 20
    assert pre["corr_ln_all_ln_unrelated"] == ext.corr_between(TABLES["all"][0], TABLES["unrelated"][0])
    # the drops do not depend on the basis or den
    assert ext.pregate(phi_m, loaded)["dropped"] == pre["dropped"]


@pytest.mark.parametrize("case,want", [
    ("one_study", {"all": "1 study", "unrelated": "1 study"}),
    ("cardiac_only", {"all": "only cardiac-tier"}),
    ("all_is_cardiac", {"all": "its bytes equal ['cardiac', 'study_GSE167493']"}),
    ("few_genes", {"all": "genes != 1 (< 10)"}),
    ("aligned", {"all": "|corr(ln phi_M, ln T)|"}),
    ("same_bytes", {"unrelated": "its bytes equal all's"}),
    ("unlisted", {"unrelated": "not listed"}),
])
def test_pregate_drops(tmp_path, lst, monkeypatch, case, want):
    phi_m = _phi_m()
    tabs = dict(TABLES)
    listed = None
    if case == "one_study":
        tabs["all"] = (tabs["all"][0], ["GSE167493"], ["cardiac", "unrelated"])
        tabs["unrelated"] = (tabs["unrelated"][0], ["GSE204908"], ["unrelated"])
    elif case == "cardiac_only":
        tabs["all"] = (tabs["all"][0], ["GSE167493", "GSE1"], ["cardiac"])
    elif case == "all_is_cardiac":
        tabs["all"] = (tabs["cardiac"][0], ["GSE167493", "GSE1"], ["cardiac", "unrelated"])
    elif case == "few_genes":
        tabs["all"] = (table_values(5, 9), tabs["all"][1], tabs["all"][2])
    elif case == "aligned":
        T = np.ones(G)
        T[:20] = phi_m[:20] ** 0.5
        tabs["all"] = (T, tabs["all"][1], tabs["all"][2])
    elif case == "same_bytes":
        tabs["unrelated"] = (tabs["all"][0], tabs["unrelated"][1], tabs["unrelated"][2])
    elif case == "unlisted":
        listed = ["all", "cardiac", "study_GSE167493"]
    s = setup_tables(tmp_path, lst, monkeypatch, tables=tabs, listed=listed)
    pre = ext.pregate(phi_m, ext.load_listed(GENES, env=s.env))
    assert set(pre["dropped"]) == set(want), pre["dropped"]
    for n, why in want.items():
        assert why in pre["dropped"][n], (n, pre["dropped"][n])


def test_a_listed_table_of_other_bytes_refuses_the_run(tmp_path, lst, monkeypatch):
    s = setup_tables(tmp_path, lst, monkeypatch)
    p = s.d / "unrelated" / "factors" / "unrelated.csv"
    p.write_bytes(p.read_bytes() + b"\n")
    with pytest.raises(ext.TableRefused, match="not the listed"):
        ext.load_listed(GENES, env=s.env)


# ------------------------------------------------------------------------------------------------ the builder


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    from vec.boards import Board

    d = tmp_path_factory.mktemp("t3outside")
    wt = d / "wt.h5ad"
    tx._write_wt(wt, np.random.default_rng(0))
    run = d / "run"
    fdir = run / "factors"
    fdir.mkdir(parents=True)
    p = fdir / f"{bo.TABLE_M}.csv"
    pd.DataFrame({"gene": GENES, "factor": PHI}).to_csv(p, index=False)
    (fdir / "factors.json").write_text(json.dumps({
        "preset": "gata4", "files": {bo.TABLE_M: {"file": str(p), "sha256": t3d.sha256(p)}}, "parameters": {},
        "inputs": {"wt": {"path": "wt", "sha256": t3d.sha256(wt)}}}))
    board = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, min_cells=200, max_cells=3000)
    dec = tx._decision(d / "gata4_decisions.json")
    return SimpleNamespace(d=d, wt=wt, run=run, board=board, dec=dec, pin=t3d.sha256(p))


@pytest.fixture
def bw(world, tmp_path, lst, monkeypatch):
    """The synthetic phi_M's pin, tables and listing, and a committed K22 reading naming x_all and x_unrel."""
    monkeypatch.setattr(bo, "TABLE_SHA256", world.pin)
    s = setup_tables(tmp_path, lst, monkeypatch)
    pins = {n: v for n, v in ext.K22_TABLE_SHA256.items() if v}

    def reading(names=("x_all", "x_unrel"), void=False, tables=None, rel="reading.json"):
        body = {"incomplete": False, "void": void, "names": list(names), "outcome": "test",
                "ext": {"listing_sha256": ext.K22_LISTING_SHA256, "tables": tables or pins}}
        return commit_file(lst, rel, json.dumps(body))

    s.reading = reading(rel="reading.json")
    s.make_reading = reading
    return s


def build(world, bw, form, out: Path, **kw) -> dict:
    return bo.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=kw.pop("run_dir", world.run),
                         seed=SEED, threads=2, env=kw.pop("env", bw.env), reading=kw.pop("reading", bw.reading), **kw)


class Spy:
    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


@pytest.mark.parametrize("base", ["s1", "s2"])
def test_identity_gives_the_frozen_barred_file_byte_for_byte(world, bw, tmp_path, base):
    frozen = bs.build_submission(world.dec, "barred", base, world.wt, world.board, N, tmp_path / "frozen.h5ad",
                                 run_dir=world.run, seed=SEED, threads=2)
    side = build(world, bw, f"identity-{base}", tmp_path / "id.h5ad", env={}, reading=tmp_path / "none.json")
    assert (tmp_path / "frozen.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()
    assert side["output"]["sha256"] == frozen["output"]["sha256"] and side["construction"] == frozen["construction"]
    ob = json.loads((tmp_path / "id.json").read_text())["outside_build"]
    assert ob["form"] == f"identity-{base}" and ob["table"] is None and ob["k22_reading"] is None
    assert ob["script_sha256"] == t3d.sha256(T3 / "build_outside_form.py")
    assert ob["outside_forms_sha256"] == t3d.sha256(T3 / "outside_forms.py")
    assert ob["frozen_files_sha256"] == bo.FROZEN_SHA256 and ob["factor"]["applied"]["genes_changed"] == 0
    assert not any("cell" in k.lower() for k in ob)
    assert bs.load_factor.__module__ == "build_submission"  # left as it was


def test_identity_gives_the_exponent_builders_file_byte_for_byte(world, bw, tmp_path):
    ref = bx.build_exponent(world.dec, "2.0", world.wt, world.board, N, tmp_path / "bx.h5ad", run_dir=world.run,
                            seed=SEED, threads=2)
    side = build(world, bw, "identity-s1x2.0", tmp_path / "id.h5ad")
    assert (tmp_path / "bx.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()
    assert side["output"]["sha256"] == ref["output"]["sha256"]


@pytest.mark.parametrize("form", ["xall-s1", "xunrel-s1", "xall-s1x1.5", "xunrel-s1x2.0", "xall-s2"])
def test_the_factor_applied_is_the_shared_product_at_the_base_s_exponent(world, bw, tmp_path, monkeypatch, form):
    spy = Spy(monkeypatch)
    calls = []
    original = ext.product
    monkeypatch.setattr(ext, "product", lambda *a, **k: calls.append(a[1]) or original(*a, **k))
    side = build(world, bw, form, tmp_path / "f.h5ad")
    (phi,) = spy.seen
    kind, base = form.split("-", 1)
    x = bo.BASES[base][1]
    T = TABLES[bo.KINDS[kind][0]][0]
    assert calls == [x]  # the one implementation the practice arm applies
    assert np.array_equal(phi, (PHI if x == 1.0 else PHI ** x) * T)
    ob = side["outside_build"]
    assert ob["form"] == form and ob["base"]["x"] == x and ob["exponent"] == 1.0
    assert ob["table"]["name"] == bo.KINDS[kind][0] and ob["table"]["sha256"] == bw.shas[bo.KINDS[kind][0]]
    assert ob["table"]["listing_sha256"] == ext.K22_LISTING_SHA256 and ob["k22_reading"]["names"] == ["x_all",
                                                                                                     "x_unrel"]
    assert ob["table"]["register_accessions"]["GSE204908" if kind == "xunrel" else "GSE167493"] == (
        ["GSE204906"] if kind == "xunrel" else ["GSE168768", "GSE186934"])
    assert side["inputs"]["factor_note"] == bo.FACTOR_NOTE and side["construction"]["factor_exponent"] == 1.0


def test_an_outside_form_draws_the_same_cells_as_its_base(world, bw, tmp_path):
    import anndata as ad

    build(world, bw, "identity-s1", tmp_path / "id.h5ad")
    build(world, bw, "xall-s1", tmp_path / "x.h5ad")
    a, b = ad.read_h5ad(tmp_path / "id.h5ad"), ad.read_h5ad(tmp_path / "x.h5ad")
    assert list(a.obs_names) == list(b.obs_names) and a.obs.equals(b.obs)
    assert np.array_equal(a.obsm["spatial_3D"], b.obsm["spatial_3D"])
    assert not np.array_equal(np.asarray(a.X), np.asarray(b.X))


@pytest.mark.parametrize("form", ["xall_perm-s1", "xunrel_perm-s2", "rep_all-s1", "x_all-s1", "barred_form-s1",
                                  "xall", "xall-s3", "XALL-s1", "xall-s1 ", "cardiac-s1", "wnt-s1", "study_x-s1"])
def test_controls_and_other_forms_are_refused_by_name(world, bw, tmp_path, form):
    with pytest.raises(bo.Refused):
        build(world, bw, form, tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


def test_the_reading_must_name_the_form(world, bw, tmp_path, monkeypatch):
    for reading, why in ((bw.make_reading(names=["x_unrel"], rel="r1.json"), "does not name x_all"),
                         (bw.make_reading(void=True, rel="r2.json"), "VOID"),
                         (bw.make_reading(tables={"all": "0" * 64}, rel="r3.json"), "other tables"),
                         (tmp_path / "none.json", "no K22 reading")):
        with pytest.raises(bo.Refused, match=why):
            build(world, bw, "xall-s1", tmp_path / "x.h5ad", reading=reading)
    loose = tmp_path / "loose.json"
    loose.write_bytes(bw.reading.read_bytes())
    with pytest.raises(bo.Refused, match="not committed"):
        build(world, bw, "xall-s1", tmp_path / "x.h5ad", reading=loose)
    monkeypatch.setattr(ext, "K22_LISTING_SHA256", None)
    with pytest.raises(bo.Refused, match="not filled"):
        build(world, bw, "xall-s1", tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


def test_the_builder_refuses_what_the_loader_refuses(world, bw, tmp_path):
    with pytest.raises(bo.Refused, match="not set"):
        build(world, bw, "xall-s1", tmp_path / "x.h5ad", env={})
    p = bw.d / "all" / "factors" / "all.csv"
    good = p.read_bytes()
    p.write_bytes(good + b"\n")
    with pytest.raises(bo.Refused, match="not the listed"):
        build(world, bw, "xall-s1", tmp_path / "x.h5ad")
    p.write_bytes(good)
    assert not (tmp_path / "x.h5ad").exists()


def test_the_form_s_own_refusal_is_distinct(world, tmp_path, lst, monkeypatch):
    monkeypatch.setattr(bo, "TABLE_SHA256", world.pin)
    few = np.ones(G)
    few[30:33] = 2.0  # three genes: the product changes fewer than 10 genes
    s = setup_tables(tmp_path, lst, monkeypatch, tables={**TABLES, "all": (few, ["a", "b"], ["unrelated"])})
    pins = {n: v for n, v in ext.K22_TABLE_SHA256.items() if v}
    r = commit_file(lst, "reading.json", json.dumps({"incomplete": False, "void": False, "names": ["x_all"],
                                                     "ext": {"listing_sha256": ext.K22_LISTING_SHA256,
                                                             "tables": pins}}))
    with pytest.raises(bo.FormRefused, match="fewer than 10|< 10"):
        bo.build_form(world.dec, "xall-s1", world.wt, world.board, N, tmp_path / "x.h5ad", run_dir=world.run,
                      seed=SEED, threads=2, env=s.env, reading=r)
    assert not (tmp_path / "x.h5ad").exists()


def test_changed_frozen_files_and_an_unclean_tree_are_refused(world, bw, tmp_path, monkeypatch, own_repo):
    import shutil

    changed = tmp_path / "t3_decompose.py"
    shutil.copy(bs.T3D_PATH, changed)
    with changed.open("a") as fh:
        fh.write("\n# changed\n")
    with monkeypatch.context() as m:
        m.setattr(bs, "T3D_PATH", changed)
        with pytest.raises(bo.Refused, match="frozen file"):
            build(world, bw, "xall-s1", tmp_path / "a.h5ad")
    (own_repo / "stray.py").write_text("x = 1\n")
    with pytest.raises(bo.Refused, match="not clean"):
        build(world, bw, "xall-s1", tmp_path / "b.h5ad")
    assert not (tmp_path / "a.h5ad").exists() and not (tmp_path / "b.h5ad").exists()


def _boards(tmp_path: Path) -> Path:
    boards = tmp_path / "boards"
    boards.mkdir()
    (boards / "T3__synthetic.genes.txt").write_text("\n".join(GENES) + "\n")
    (boards / "index.json").write_text(json.dumps({"T3:synthetic": {
        "task": "T3", "genes_file": "T3__synthetic.genes.txt", "needs_coords": True, "min_cells": 200,
        "max_cells": N}}))
    return boards


def test_cli_takes_names_only_and_exit_codes(world, bw, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bo, "READING", bw.reading)
    monkeypatch.setenv(ext.ENV, bw.env[ext.ENV])
    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(_boards(tmp_path)), "--run-dir", str(world.run), "--threads", "2"]
    runner = CliRunner()
    res = runner.invoke(bo.app, [*base, "--form", "xall_perm-s1", "--out", str(tmp_path / "a.h5ad")])
    assert res.exit_code == 2 and "refused" in res.output and not (tmp_path / "a.h5ad").exists()
    res = runner.invoke(bo.app, [*base, "--form", "xunrel-s1", "--out", str(tmp_path / "b.h5ad")])
    assert res.exit_code == 0, res.output
    for extra in (["--seed", "1000002"], ["--n-cells", str(N)], ["--table", "all"], ["--exponent", "0.5"],
                  ["--tables-dir", str(bw.d)], ["--path", str(bw.d / "all" / "factors" / "all.csv")]):
        out = tmp_path / f"{extra[0].strip('-')}.h5ad"
        res = runner.invoke(bo.app, [*base, "--form", "xall-s1", *extra, "--out", str(out)])
        assert res.exit_code != 0 and not out.exists(), extra


def test_the_cli_s_form_refusal_exit_code(world, tmp_path, lst, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bo, "TABLE_SHA256", world.pin)
    few = np.ones(G)
    few[30:33] = 2.0
    s = setup_tables(tmp_path, lst, monkeypatch, tables={**TABLES, "all": (few, ["a", "b"], ["unrelated"])})
    pins = {n: v for n, v in ext.K22_TABLE_SHA256.items() if v}
    r = commit_file(lst, "reading.json", json.dumps({"incomplete": False, "void": False, "names": ["x_all"],
                                                     "ext": {"listing_sha256": ext.K22_LISTING_SHA256,
                                                             "tables": pins}}))
    monkeypatch.setattr(bo, "READING", r)
    monkeypatch.setenv(ext.ENV, s.env[ext.ENV])
    res = CliRunner().invoke(bo.app, ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic",
                                      "--boards-dir", str(_boards(tmp_path)), "--run-dir", str(world.run),
                                      "--threads", "2", "--form", "xall-s1", "--out", str(tmp_path / "c.h5ad")])
    assert res.exit_code == 4 and "own refusal" in res.output and not (tmp_path / "c.h5ad").exists()


def test_the_constants_of_the_launchers():
    sh = (T3 / "build_outside_form.sh").read_text()
    assert "VEC_T3_EXT_TABLES=\"$TABLES\"" in sh and "TABLES=$HOME/v/vec-ext/t3/tables" in sh
    assert "if [ $b = 4 ]; then return 3; fi" in sh and "finish 46" in sh
    for b, (_, _, sha) in bo.BASES.items():
        assert f"[{b}]={sha}" in sh
    run = (T3 / "run_k22.sh").read_text()
    assert "VEC_T3_EXT_TABLES=\"$TABLES\"" in run and "k22_outside_forms.py" in run
    assert "if [ $rc = 5 ]; then finish 25; fi" in run and "if [ $rc = 4 ]; then finish 44; fi" in run
    assert math.isclose(ext.EXPONENT, 1.0) and os.sep == "/"
