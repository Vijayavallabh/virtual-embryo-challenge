"""Tests for experiments/t3_decomp/build_k24_form.py (the thin builder of the one form row K24's reading may name) and
its launcher build_k24_form.sh. Offline, synthetic and fast: the synthetic world and throwaway-repository fixture of
tests/test_t3_barred_exponent_build.py (as tests/test_t3_outside_forms.py uses them), hand-made tables in a temporary
folder outside the repository, a committed synthetic K24 reading. No competition data and no outside data."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
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
bo = sys.modules.get("build_outside_form") or _load("build_outside_form", T3 / "build_outside_form.py")
bk24 = sys.modules.get("build_k24_form") or _load("build_k24_form", T3 / "build_k24_form.py")
import k24_wider_tables as k24  # noqa: E402
import outside_forms as ext  # noqa: E402

bs = bo.bs
t3d = bs.t3d
own_repo = tx.own_repo  # builds are stamped against a throwaway repository
SEED, N, G, GENES, PHI = tx.SEED, tx.N, tx.G, tx.GENES, tx.PHI
TWO_SHA = "9490f524e40e146480c775f85ae8829ec9c005607104ad8ace129c5c46ed6228"
REGISTER = ("| X013 | GEO GSE204906, a SubSeries of GSE204908 |\n| X014 | GEO GSE277839 |\n"
            "| X015 | GEO GSE299068 |\n| X017 | ArrayExpress E-MTAB-7322 |\n")


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _values(seed: int, moved: int = 30, scale: float = 0.05) -> np.ndarray:
    rng = np.random.default_rng(seed)
    T = np.ones(G)
    idx = rng.choice(G, moved, replace=False)
    T[idx] = np.exp(rng.normal(0, scale, moved))
    return T


def _write_table(d: Path, name: str, vals, man: dict) -> str:
    fd = d / name / "factors"
    fd.mkdir(parents=True, exist_ok=True)
    p = fd / f"{name}.csv"
    pd.DataFrame({"gene": GENES, "factor": vals}).to_csv(p, index=False)
    sha = _sha(p)
    (fd / "factors.json").write_text(json.dumps({"files": {name: {"file": str(p), "sha256": sha}}}))
    (d / name / f"{name}_manifest.json").write_text(json.dumps({**man, "name": name,
                                                                 "files": {name: {"file": str(p), "sha256": sha}}}))
    return sha


def commit_file(repo: Path, rel: str, text_: str) -> Path:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text_)
    tx._git(repo, "add", rel)
    tx._git(repo, "commit", "-q", "--allow-empty", "-m", f"add {rel}")
    return p


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    from vec.boards import Board

    d = tmp_path_factory.mktemp("t3k24")
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


def setup_tables(tmp_path: Path, monkeypatch, sc3=None, bulk=None):
    """sc3 and bulk in combine.py's layout, the listing (with `two`'s pinned sha256) committed in a throwaway
    repository, the pins patched, a register file; returns what the builds need."""
    repo = tmp_path / "lst"
    repo.mkdir()
    tx._git(repo, "init", "-q")
    d = tmp_path / "tables"
    sc3 = _values(1) if sc3 is None else sc3
    bulk = _values(2, scale=0.08) if bulk is None else bulk
    shas = {"sc3": _write_table(d, "all", sc3, {"studies": ["GSE204908", "GSE277839", "GSE299068"],
                                               "tiers": ["unrelated"], "sign_rule": "applied",
                                               "refstats_sha256": "r" * 64}),
            "bulk": _write_table(d, "bulk_all", bulk, {"studies": ["MGP-DMDD"], "tiers": ["unrelated"],
                                                      "modality": "bulk", "sign_rule": "single study",
                                                      "refstats_sha256": "r" * 64}),
            "two": TWO_SHA}
    (d / "COMBINE.json").write_text(json.dumps({"tables": {"all": {"all": {"sha256": shas["sc3"]}}}, "record": {}}))
    (d / "COMBINE_bulk.json").write_text(json.dumps({"tables": {"bulk_all": {"bulk_all": {"sha256": shas["bulk"]}}},
                                                     "record": {}}))
    lp = commit_file(repo, "TABLES_K24.sha256", "".join(f"{shas[k]}  {k}\n" for k in k24.KEYS))
    monkeypatch.setattr(k24, "LISTING_FILE", lp)
    monkeypatch.setattr(k24, "K24_LISTING_SHA256", _sha(lp))
    monkeypatch.setattr(k24, "K24_TABLE_SHA256", {"two": TWO_SHA, "sc3": shas["sc3"], "bulk": shas["bulk"]})
    e = {"sc3": k24.exponent(sc3), "bulk": k24.exponent(bulk)}
    c = k24.joint_c(sc3, bulk, e["sc3"], e["bulk"])
    reg = tmp_path / "EXTERNAL_DATA.md"
    reg.write_text(REGISTER)
    s = SimpleNamespace(repo=repo, d=d, env={ext.ENV: str(d)}, shas=shas, sc3=sc3, bulk=bulk, e=e, c=c, reg=reg)

    def reading(names=("x_sc3",), void=False, incomplete=False, tables=None, dropped=None, e_=None, c_=None,
                rel="reading.json"):
        pre = {"dropped": dict(dropped or {}), "e": dict(e_ or s.e), "c": s.c if c_ is None else c_}
        body = {"incomplete": incomplete, "void": void, "names": list(names), "outcome": "test",
                "arms": list(k24.active_arms(pre)), "pregate": pre,
                "ext": {"listing_sha256": k24.K24_LISTING_SHA256, "tables": tables or dict(k24.K24_TABLE_SHA256)}}
        return commit_file(repo, rel, json.dumps(body))

    s.make_reading = reading
    return s


@pytest.fixture
def tb(world, tmp_path, monkeypatch):
    monkeypatch.setattr(bo, "TABLE_SHA256", world.pin)
    monkeypatch.setattr(bk24, "TABLE_SHA256", world.pin)
    return setup_tables(tmp_path, monkeypatch)


def build(world, tb, form, out: Path, reading=None, **kw) -> dict:
    return bk24.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=world.run, seed=SEED, threads=2,
                           env=kw.pop("env", tb.env), reading=reading, register=kw.pop("register", tb.reg), **kw)


class Spy:
    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


def test_constants_and_the_reused_builder():
    assert bk24.FORMS == tuple(f"{k}-{b}" for k in ("k24sc3", "k24bulk", "k24joint")
                               for b in ("s1", "s1x1.5", "s1x2.0", "s2"))
    assert bk24.BASES is bo.BASES and bk24.TABLE_M == bo.TABLE_M == "transfer_mab21l2_all"
    assert bk24.READING == T3 / "results" / "k24_wider_tables" / "reading.json"
    assert bk24.PRICED_AT == "x = 1 only" and "row K24" in bk24.LIST_ENTRY and "2026-11-21 (c)" in bk24.LIST_ENTRY
    assert {k: v[0] for k, v in bk24.KINDS.items()} == {"k24sc3": "x_sc3", "k24bulk": "x_bulk", "k24joint": "x_joint"}


@pytest.mark.parametrize("form", ["k24sc3-s1", "k24bulk-s1x1.5", "k24joint-s1x2.0", "k24joint-s2"])
def test_the_factor_applied_is_k24s_scaled_table_on_the_base(world, tb, tmp_path, monkeypatch, form):
    kind, base = form.split("-", 1)
    arm = bk24.KINDS[kind][0]
    r = tb.make_reading(names=(arm,), rel=f"{form}.json")
    spy = Spy(monkeypatch)
    side = build(world, tb, form, tmp_path / "f.h5ad", reading=r)
    (phi,) = spy.seen
    x = bo.BASES[base][1]
    phi_x = PHI if x == 1.0 else PHI ** x
    e, c = tb.e, tb.c
    want = {"k24sc3": tb.sc3 ** e["sc3"], "k24bulk": tb.bulk ** e["bulk"],
            "k24joint": tb.sc3 ** (c * e["sc3"]) * tb.bulk ** (c * e["bulk"])}[kind]
    assert np.allclose(phi, phi_x * want, rtol=1e-14, atol=0)
    assert np.array_equal(phi, ext.product(PHI, x, k24.k24_table(bk24.KINDS[kind][1], {"sc3": tb.sc3, "bulk": tb.bulk},
                                                                   {"e_sc3": e["sc3"], "e_bulk": e["bulk"], "c": c})))
    ob = side["outside_build"]
    assert ob["form"] == form and ob["arm"] == arm and ob["base"]["x"] == x and ob["priced_at"] == "x = 1 only"
    assert ob["k24_reading"]["names"] == [arm] and ob["k24_reading"]["sha256"] == _sha(r)
    assert ob["script_sha256"] == t3d.sha256(T3 / "build_k24_form.py") and ob["s_star"] == k24.S_STAR
    assert ob["c"] == (c if kind == "k24joint" else None)
    keys = bk24.KINDS[kind][2]
    assert set(ob["tables"]) == set(keys) and set(ob["e"]) == {f"e_{k}" for k in keys}
    if "sc3" in keys:
        reg = ob["tables"]["sc3"]["register"]
        assert reg["GSE204908"]["rows"] == ["X013"] and reg["GSE299068"]["rows"] == ["X015"]
        assert ob["tables"]["sc3"]["sha256"] == tb.shas["sc3"]
    if "bulk" in keys:
        assert ob["tables"]["bulk"]["register"]["MGP-DMDD"] == {"accessions": ["E-MTAB-7322"], "rows": ["X017"]}
        assert ob["tables"]["bulk"]["combine_file"] == "COMBINE_bulk.json"
    assert side["inputs"]["factor_note"] == bo.FACTOR_NOTE and side["construction"]["factor_exponent"] == 1.0
    assert not any("cell" in k.lower() for k in ob)
    assert bs.load_factor.__module__ == "build_submission"


def test_a_k24_form_draws_the_same_cells_as_its_base(world, tb, tmp_path):
    import anndata as ad

    bo.build_form(world.dec, "identity-s1", world.wt, world.board, N, tmp_path / "id.h5ad", run_dir=world.run,
                  seed=SEED, threads=2, env={}, reading=tmp_path / "none.json")
    build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=tb.make_reading())
    a, b = ad.read_h5ad(tmp_path / "id.h5ad"), ad.read_h5ad(tmp_path / "x.h5ad")
    assert list(a.obs_names) == list(b.obs_names) and a.obs.equals(b.obs)
    assert not np.array_equal(np.asarray(a.X), np.asarray(b.X))


@pytest.mark.parametrize("form", ["x_sc3-s1", "x_sc3_perm-s1", "x_sc3_size-s1", "x_e2-s1", "barred_form-s1",
                                  "identity-s1", "xall-s1", "k24sc3-s3", "k24two-s1", "k24sc3", "K24SC3-s1",
                                  "k24joint_perm-s1", "k24sc3-s1 "])
def test_other_forms_are_refused_by_name(world, tb, tmp_path, form):
    with pytest.raises(bk24.Refused):
        build(world, tb, form, tmp_path / "x.h5ad", reading=tb.make_reading(rel="r.json"))
    assert not (tmp_path / "x.h5ad").exists()


def test_the_reading_must_name_exactly_the_arm_undropped_with_its_constants(world, tb, tmp_path):
    cases = ((dict(names=("x_bulk",)), "does not name exactly x_sc3"),
             (dict(names=("x_sc3", "x_bulk")), "does not name exactly"),
             (dict(void=True), "VOID"), (dict(incomplete=True), "VOID or incomplete"),
             (dict(tables={"two": TWO_SHA, "sc3": "0" * 64, "bulk": tb.shas["bulk"]}), "other tables"),
             (dict(dropped={"sc3": "test"}), "dropped"),
             (dict(e_={"sc3": tb.e["sc3"] * (1 + 1e-9), "bulk": tb.e["bulk"]}), "e_sc3"))
    for i, (kw, why) in enumerate(cases):
        with pytest.raises(bk24.Refused, match=why):
            build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=tb.make_reading(rel=f"r{i}.json", **kw))
    with pytest.raises(bk24.Refused, match="^c: "):
        build(world, tb, "k24joint-s1", tmp_path / "x.h5ad",
              reading=tb.make_reading(names=("x_joint",), c_=tb.c * (1 + 1e-9), rel="rc.json"))
    with pytest.raises(bk24.Refused, match="no K24 reading"):
        build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=tmp_path / "none.json")
    loose = tmp_path / "loose.json"
    loose.write_bytes(tb.make_reading(rel="ok.json").read_bytes())
    with pytest.raises(bk24.Refused, match="not committed"):
        build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=loose)
    # e within the relative 1e-12 is the reading's
    build(world, tb, "k24sc3-s1", tmp_path / "ok.h5ad",
          reading=tb.make_reading(e_={"sc3": tb.e["sc3"] * (1 + 1e-14), "bulk": tb.e["bulk"]}, rel="near.json"))
    assert not (tmp_path / "x.h5ad").exists()


def test_pins_register_and_loader_refusals(world, tb, tmp_path, monkeypatch):
    r = tb.make_reading()
    reg = tmp_path / "reg.md"
    reg.write_text("| X013 | GSE204906 |\n")  # no row for GSE277839 or GSE299068
    with pytest.raises(bk24.Refused, match="no row naming"):
        build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=r, register=reg)
    assert not (tmp_path / "x.h5ad").exists() and not (tmp_path / "x.json").exists()
    with pytest.raises(bk24.Refused, match="not set"):
        build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=r, env={})
    p = tb.d / "all" / "factors" / "all.csv"
    good = p.read_bytes()
    p.write_bytes(good + b"\n")
    with pytest.raises(bk24.Refused, match="not the listed"):
        build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=r)
    p.write_bytes(good)
    monkeypatch.setattr(k24, "K24_LISTING_SHA256", None)
    with pytest.raises(bk24.Refused, match="not filled"):
        build(world, tb, "k24sc3-s1", tmp_path / "x.h5ad", reading=r)
    assert not (tmp_path / "x.h5ad").exists()


def test_the_form_s_own_refusal_is_distinct(world, tmp_path, monkeypatch):
    monkeypatch.setattr(bk24, "TABLE_SHA256", world.pin)
    few = np.ones(G)
    few[30:33] = [1.1, 0.9, 1.05]  # three genes: the product changes fewer than 10 genes
    s = setup_tables(tmp_path, monkeypatch, sc3=few)
    with pytest.raises(bk24.FormRefused, match="fewer than 10|< 10"):
        build(world, s, "k24sc3-s1", tmp_path / "x.h5ad", reading=s.make_reading())
    assert not (tmp_path / "x.h5ad").exists()


def test_changed_frozen_files_and_an_unclean_tree_are_refused(world, tb, tmp_path, monkeypatch, own_repo):
    import shutil

    changed = tmp_path / "t3_decompose.py"
    shutil.copy(bs.T3D_PATH, changed)
    with changed.open("a") as fh:
        fh.write("\n# changed\n")
    r = tb.make_reading()
    with monkeypatch.context() as m:
        m.setattr(bs, "T3D_PATH", changed)
        with pytest.raises(bk24.Refused, match="frozen file"):
            build(world, tb, "k24sc3-s1", tmp_path / "a.h5ad", reading=r)
    (own_repo / "stray.py").write_text("x = 1\n")
    with pytest.raises(bk24.Refused, match="not clean"):
        build(world, tb, "k24sc3-s1", tmp_path / "b.h5ad", reading=r)
    assert not (tmp_path / "a.h5ad").exists() and not (tmp_path / "b.h5ad").exists()


def _boards(tmp_path: Path) -> Path:
    boards = tmp_path / "boards"
    boards.mkdir()
    (boards / "T3__synthetic.genes.txt").write_text("\n".join(GENES) + "\n")
    (boards / "index.json").write_text(json.dumps({"T3:synthetic": {
        "task": "T3", "genes_file": "T3__synthetic.genes.txt", "needs_coords": True, "min_cells": 200,
        "max_cells": N}}))
    return boards


def test_cli_takes_names_only_and_exit_codes(world, tb, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk24, "READING", tb.make_reading())
    monkeypatch.setattr(bk24, "REGISTER", tb.reg)
    monkeypatch.setenv(ext.ENV, tb.env[ext.ENV])
    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(_boards(tmp_path)), "--run-dir", str(world.run), "--threads", "2"]
    runner = CliRunner()
    res = runner.invoke(bk24.app, [*base, "--form", "x_sc3_perm-s1", "--out", str(tmp_path / "a.h5ad")])
    assert res.exit_code == 2 and "refused" in res.output and not (tmp_path / "a.h5ad").exists()
    res = runner.invoke(bk24.app, [*base, "--form", "k24sc3-s1", "--out", str(tmp_path / "b.h5ad")])
    assert res.exit_code == 0, res.output
    for extra in (["--seed", "1000002"], ["--n-cells", str(N)], ["--table", "sc3"], ["--exponent", "2.0"],
                  ["--table-exponent", "2.0"]):
        out = tmp_path / f"{extra[0].strip('-')}.h5ad"
        res = runner.invoke(bk24.app, [*base, "--form", "k24sc3-s1", *extra, "--out", str(out)])
        assert res.exit_code != 0 and not out.exists(), extra


def test_the_cli_s_form_refusal_exit_code(world, tmp_path, monkeypatch):
    from typer.testing import CliRunner

    monkeypatch.setattr(bk24, "TABLE_SHA256", world.pin)
    few = np.ones(G)
    few[30:33] = [1.1, 0.9, 1.05]
    s = setup_tables(tmp_path, monkeypatch, sc3=few)
    monkeypatch.setattr(bk24, "READING", s.make_reading())
    monkeypatch.setattr(bk24, "REGISTER", s.reg)
    monkeypatch.setenv(ext.ENV, s.env[ext.ENV])
    res = CliRunner().invoke(bk24.app, ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic",
                                        "--boards-dir", str(_boards(tmp_path)), "--run-dir", str(world.run),
                                        "--threads", "2", "--form", "k24sc3-s1", "--out", str(tmp_path / "c.h5ad")])
    assert res.exit_code == 4 and "own refusal" in res.output and not (tmp_path / "c.h5ad").exists()


# ------------------------------------------------------------------------------------------------ the launcher


def test_the_launcher_constants():
    sh = (T3 / "build_k24_form.sh").read_text()
    bsh = (T3 / "build_outside_form.sh").read_text()
    assert "HEAD_SHA=${1:?usage: build_k24_form.sh HEAD_SHA BASE KIND}" in sh
    assert "sc3|bulk|joint) ;;" in sh and "FORM=\"k24$KIND-$BASE\"" in sh
    assert "one build_outside_form.py \"identity-$1\" \"$f\" || rc=$?" in sh
    assert "one build_k24_form.py \"$FORM\" \"$f\" || rc=$?" in sh
    for b, h in bo.BASES.items():
        assert f"[{b}]={h[2]}" in sh and f"[{b}]={h[2]}" in bsh
    assert "gate s1 20" in sh and "TABLES=$HOME/v/vec-ext/t3/tables  # combine.py's --out; fixed, no override" in sh
    for code in (2, 10, 40, 41, 42, 43, 44, 45, 46, 95, 96, 97, 98):
        assert f"finish {code}" in sh, code
    assert "tests/test_build_k24_form.py tests/test_k24_wider_tables.py" in sh
    import re

    assert not re.search(r"\$\{[^}]*:\?[^}]*'[^}]*\}", sh)  # no apostrophe inside a ${var:?message}


def test_the_launcher_s_refusal_paths_run(tmp_path):
    repo = tmp_path / "repo"
    (repo / "experiments" / "t3_decomp").mkdir(parents=True)
    (repo / "experiments" / "t3_decomp" / "build_k24_form.sh").write_text((T3 / "build_k24_form.sh").read_text())
    tx._git(repo, "init", "-q")
    tx._git(repo, "add", ".")
    tx._git(repo, "commit", "-q", "-m", "init")
    home = tmp_path / "home"
    (home / "v" / "vec-ext" / "t3" / "tables").mkdir(parents=True)
    env = {**tx._git_env(), "HOME": str(home), "GIT_CEILING_DIRECTORIES": str(tmp_path)}

    def run(*args):
        r = subprocess.run(["bash", "experiments/t3_decomp/build_k24_form.sh", *args], cwd=repo, env=env,
                           capture_output=True, text=True, timeout=60)
        return r.returncode, r.stderr

    code, err = run()
    assert code != 0 and "usage: build_k24_form.sh HEAD_SHA BASE KIND" in err
    code, err = run("0" * 40, "s3", "sc3")
    assert code == 2 and "unknown base" in err
    code, err = run("0" * 40, "s1", "two")
    assert code == 2 and "unknown kind" in err
    code, err = run("0" * 40, "s1", "sc3")
    assert code == 97 and "HEAD is not" in err
    assert (home / "v" / "scratch" / "t3k24_0000000_s1" / "all.exit").read_text().strip() == "97"
    head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, env=env,
                          check=True).stdout.strip()
    (home / "v" / "vec-ext" / "t3" / "tables").rmdir()
    code, err = run(head, "s1", "joint")
    assert code == 96 and "not a real directory" in err
    assert os.environ.get("GIT_DIR") is None
