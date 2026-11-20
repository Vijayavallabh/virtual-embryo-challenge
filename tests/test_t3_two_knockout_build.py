"""Tests for experiments/t3_decomp/build_two_knockout.py: the barred file with its Mab21l2 factor combined with the
released Gata4 knockout's factor, in the fixed forms of the two-knockout list.

Offline, synthetic and fast. It reuses the helpers and the throwaway-repository fixture of
tests/test_t3_barred_exponent_build.py: a small synthetic WT file for WT E8.75, synthetic `transfer_mab21l2_all` and
`gata4_pooled` tables in the run's factors/, and a decision file with D7m and the frozen hashes. No competition data.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
T3 = ROOT / "experiments" / "t3_decomp"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


tx = sys.modules.get("_t3_barred_exponent_helpers") or _load("_t3_barred_exponent_helpers",
                                                             ROOT / "tests" / "test_t3_barred_exponent_build.py")
bk = _load("build_two_knockout", T3 / "build_two_knockout.py")
bx, bs = tx.bx, bk.bs
t3d = bs.t3d
own_repo = tx.own_repo  # the exponent builder's autouse fixture: builds are stamped against a throwaway repository
SEED, N, G, GENES = tx.SEED, tx.N, tx.G, tx.GENES
PHI_M = tx.PHI  # up 2 (0-9), down 0.5 (10-14), the clip ends 4 (15-19) and 0.25 (20-24), 1 elsewhere
PHI_G = np.ones(G)
PHI_G[0:5], PHI_G[5:10], PHI_G[10:13] = 2.0, 0.5, 4.0       # agree, cancel, oppose
PHI_G[15:18], PHI_G[18:20] = 0.25, 4.0                      # cancel the top clip end; double it (16, beyond the clip)
PHI_G[20:23], PHI_G[23:25] = 4.0, 0.25                      # cancel the bottom clip end; double it (1/16)
PHI_G[25:30], PHI_G[30:33] = 2.0, 0.25                      # moved by the Gata4 table only
# C1 = phi_M * phi_G by hand, gene by gene
C1_HAND = np.ones(G)
C1_HAND[0:5], C1_HAND[5:10], C1_HAND[10:13], C1_HAND[13:15] = 4.0, 1.0, 2.0, 0.5
C1_HAND[15:18], C1_HAND[18:20], C1_HAND[20:23], C1_HAND[23:25] = 1.0, 16.0, 1.0, 0.0625
C1_HAND[25:30], C1_HAND[30:33] = 2.0, 0.25


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    from types import SimpleNamespace

    from vec.boards import Board

    d = tmp_path_factory.mktemp("t3two")
    wt = d / "wt.h5ad"
    tx._write_wt(wt, np.random.default_rng(0))
    run = d / "run"
    write_run(run, wt, {bk.TABLE_M: PHI_M, bk.TABLE_G: PHI_G})
    board = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, min_cells=200, max_cells=3000)
    boards_dir = d / "boards"
    boards_dir.mkdir()
    (boards_dir / "T3__synthetic.genes.txt").write_text("\n".join(GENES) + "\n")
    (boards_dir / "index.json").write_text(json.dumps({"T3:synthetic": {
        "task": "T3", "genes_file": "T3__synthetic.genes.txt", "needs_coords": True, "min_cells": 200,
        "max_cells": 3000}}))
    dec = tx._decision(d / "gata4_decisions.json")
    pins = {k: t3d.sha256(run / "factors" / f"{k}.csv") for k in (bk.TABLE_M, bk.TABLE_G)}
    return SimpleNamespace(d=d, wt=wt, run=run, board=board, dec=dec, pins=pins)


def write_run(run: Path, wt: Path, tables: dict, genes: dict | None = None) -> None:
    fdir = run / "factors"
    fdir.mkdir(parents=True, exist_ok=True)
    files = {}
    for name, phi in tables.items():
        p = fdir / f"{name}.csv"
        pd.DataFrame({"gene": (genes or {}).get(name, GENES), "factor": phi}).to_csv(p, index=False)
        files[name] = {"file": str(p), "sha256": t3d.sha256(p)}
    (fdir / "factors.json").write_text(json.dumps({"preset": "gata4", "files": files, "parameters": {},
                                                   "inputs": {"wt": {"path": "wt", "sha256": t3d.sha256(wt)}}}))


@pytest.fixture(autouse=True)
def synthetic_pins(world, monkeypatch):
    """The builder pins the committed tables' sha256; here the synthetic tables' (the pin itself is tested below)."""
    monkeypatch.setattr(bk, "TABLE_SHA256", dict(world.pins))


def build(world, form, out: Path, **kw) -> dict:
    return bk.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=kw.pop("run_dir", world.run),
                         seed=SEED, threads=2, **kw)


def read(path: Path):
    import anndata as ad

    return ad.read_h5ad(path)


class Spy:
    """Records the factor apply_gene_factor receives."""

    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


# ---------------------------------------------------------------------------------------------------- identity


@pytest.mark.parametrize("base", ["s1", "s2"])
def test_identity_forms_give_the_frozen_barred_files_byte_for_byte(world, tmp_path, base):
    frozen = bs.build_submission(world.dec, "barred", base, world.wt, world.board, N, tmp_path / "frozen.h5ad",
                                 run_dir=world.run, seed=SEED, threads=2)
    side = build(world, f"identity-{base}", tmp_path / "id.h5ad")
    assert (tmp_path / "frozen.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()
    assert side["output"]["sha256"] == frozen["output"]["sha256"]
    assert side["spec"] == frozen["spec"] and side["construction"] == frozen["construction"]
    assert side["inputs"]["rng_streams"] == frozen["inputs"]["rng_streams"] == [f"submission:barred:{base}",
                                                                                f"submission:barred:{base}:assemble"]
    assert side["inputs"]["factor"] == frozen["inputs"]["factor"]
    on_disk = json.loads((tmp_path / "id.json").read_text())
    assert set(on_disk) == set(frozen) | {"combined_build"}
    cb = on_disk["combined_build"]
    assert cb["form"] == f"identity-{base}" and cb["base"]["reference_sha256"] == bk.BASES[base][2]
    assert cb["script"] == "experiments/t3_decomp/build_two_knockout.py"
    assert cb["script_sha256"] == t3d.sha256(T3 / "build_two_knockout.py")
    assert cb["frozen_files_sha256"] == bk.FROZEN_SHA256 and cb["tree_at_start"]["clean"] is True
    assert cb["combination"]["gata4_weight"] == 0.0 and cb["combination"]["c"] == 1.0
    assert set(cb["tables"]) == {bk.TABLE_M, bk.TABLE_G}
    assert all(t["sha256"] == t["pinned_sha256"] for t in cb["tables"].values())
    assert cb["depth_step"] is None and "list_entry" in cb
    assert not any("cell" in k.lower() for k in cb) and not any("cell" in k.lower() for k in cb["combination"])
    stamp = json.loads((tmp_path / "id.h5ad.stamp.json").read_text())
    assert stamp["sha256"] == side["output"]["sha256"] and stamp["clean"] is True
    # the frozen loader is left as it was
    assert bs.load_factor.__module__ == "build_submission" and bs.load_factor.__name__ == "load_factor"


@pytest.mark.parametrize("x", ["1.5", "2.0"])
def test_identity_forms_give_the_exponent_builders_files_byte_for_byte(world, tmp_path, x):
    ref = bx.build_exponent(world.dec, x, world.wt, world.board, N, tmp_path / "bx.h5ad", run_dir=world.run,
                            seed=SEED, threads=2)
    side = build(world, f"identity-s1x{x}", tmp_path / "id.h5ad")
    assert (tmp_path / "bx.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()
    assert side["output"]["sha256"] == ref["output"]["sha256"]
    assert side["inputs"]["rng_streams"] == ref["inputs"]["rng_streams"]


# ---------------------------------------------------------------------------------------- the factor of each form


def test_combine_by_hand_on_a_small_case():
    a = math.log(2.0)
    m = np.array([2.0, 0.5, 1.0, 4.0])
    g = np.array([2.0, 0.5, 1.0, 0.25])
    # identity: the base table, at x = 1 the very array
    phi, rec = bk.combine("identity", 1.0, m, g)
    assert phi is not None and np.array_equal(phi, m) and rec["c"] == 1.0 and rec["gata4_weight"] == 0.0
    phi, _ = bk.combine("identity", 2.0, m, g)
    assert np.array_equal(phi, [4.0, 0.25, 1.0, 16.0])
    # C1: the product, beyond the clip where both move the same way, never rescaled or clipped again
    phi, rec = bk.combine("c1", 1.0, m, g)
    assert np.array_equal(phi, [4.0, 0.25, 1.0, 1.0]) and rec["c"] == 1.0 and rec["gata4_weight"] == 1.0
    # C2 at x = 1: ln m = a[1, -1, 0, 2], ln g = a[1, -1, 0, -2]; SD(ln m) = a sqrt5 / 2, SD(ln m + ln g) = a sqrt2
    c = math.sqrt(10.0) / 4.0
    phi, rec = bk.combine("c2", 1.0, m, g)
    assert rec["c"] == pytest.approx(c, rel=1e-14, abs=0)
    assert np.allclose(phi, [4.0 ** c, 0.25 ** c, 1.0, 1.0], rtol=1e-13, atol=0)
    assert rec["sd_ln_base"] == pytest.approx(a * math.sqrt(5.0) / 2.0, rel=1e-14)
    assert rec["sd_ln_base_plus_g"] == pytest.approx(a * math.sqrt(2.0), rel=1e-14)
    assert rec["factor"]["sd_ln"] == pytest.approx(rec["sd_ln_base"], rel=1e-12, abs=0)
    # C2 at x = 2: 2 ln m + ln g = a[3, -3, 0, 2]; SD(2 ln m) = a sqrt5, SD(sum) = a sqrt21 / 2
    c = 2.0 * math.sqrt(5.0 / 21.0)
    phi, rec = bk.combine("c2", 2.0, m, g)
    assert rec["c"] == pytest.approx(c, rel=1e-14, abs=0) and rec["mab21l2_exponent"] == pytest.approx(2.0 * c)
    assert np.allclose(phi, [2.0 ** (3 * c), 2.0 ** (-3 * c), 1.0, 2.0 ** (2 * c)], rtol=1e-13, atol=0)
    assert rec["factor"]["sd_ln"] == pytest.approx(a * math.sqrt(5.0), rel=1e-12, abs=0)
    # the common shift, from the factor alone: C1's ln is a[2, -2, 0, 0], mean 0, cos^2 0; the base's is a[1,-1,0,2]
    _, rec = bk.combine("c1", 1.0, m, g)
    assert rec["factor"]["gene_mean_ln"] == pytest.approx(0.0, abs=1e-15)
    assert rec["factor"]["cos2_ln_with_ones"] == pytest.approx(0.0, abs=1e-15)
    assert rec["base_factor"]["gene_mean_ln"] == pytest.approx(a / 2.0, rel=1e-14)
    assert rec["base_factor"]["cos2_ln_with_ones"] == pytest.approx(4.0 / (4 * 6.0), rel=1e-14)
    assert rec["corr_ln_m_ln_g"] == pytest.approx(np.corrcoef(np.log(m), np.log(g))[0, 1], rel=1e-14)


def test_combine_refuses_tables_of_other_genes_and_bad_values():
    m = np.array([2.0, 0.5, 1.0])
    with pytest.raises(bk.Refused, match="same genes"):
        bk.combine("c1", 1.0, m, np.array([2.0, 0.5]))
    with pytest.raises(bk.Refused, match="finite and positive"):
        bk.combine("c1", 1.0, m, np.array([2.0, 0.0, 1.0]))
    with pytest.raises(bk.Refused, match="c is undefined"):
        bk.combine("c2", 1.0, np.ones(3), np.ones(3))
    with pytest.raises(bk.Refused, match="kind"):
        bk.combine("c3", 1.0, m, m)


@pytest.mark.parametrize("form", ["c1-s1", "c1-s2", "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2", "identity-s1x1.5"])
def test_the_factor_applied_is_the_forms_combination(world, tmp_path, monkeypatch, form):
    spy = Spy(monkeypatch)
    side = build(world, form, tmp_path / "f.h5ad")
    assert len(spy.seen) == 1
    phi = spy.seen[0]
    kind, base = form.split("-", 1)
    x = bk.BASES[base][1]
    lm, lg = np.log(PHI_M), np.log(PHI_G)
    if kind == "c1":
        assert np.array_equal(phi, C1_HAND)                     # by hand; no rescale, no clip after the product
    elif kind == "identity":
        assert np.allclose(phi, PHI_M ** x, rtol=1e-15, atol=0)
    else:
        c = np.std(x * lm) / np.std(x * lm + lg)                 # the rule, written out again
        assert np.allclose(phi, np.exp(c * (x * lm + lg)), rtol=1e-12, atol=0)
        assert np.std(np.log(phi)) == pytest.approx(np.std(x * lm), rel=1e-12, abs=0)   # held at the base's size
        assert side["combined_build"]["combination"]["c"] == pytest.approx(c, rel=1e-13)
    cb = side["combined_build"]
    assert cb["form"] == form and cb["base"]["x"] == x and side["construction"]["factor_exponent"] == 1.0
    assert side["spec"]["factor"] == bk.FROZEN_FACTOR and side["construction"]["depth_ratios"] is None
    assert side["construction"]["stage_alphas"] is None and not side["spec"].get("depth")


# ------------------------------------------------------------------------------------------- same cells, new expression


@pytest.mark.parametrize("base,forms", [("s1", ["c1-s1", "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "identity-s1x2.0"]),
                                        ("s2", ["c1-s2", "c2-s2"])])
def test_every_form_draws_the_same_cells_as_its_base(world, tmp_path, base, forms):
    ref_side = build(world, f"identity-{base}", tmp_path / "id.h5ad")
    ref = read(tmp_path / "id.h5ad")
    for f in forms:
        s = build(world, f, tmp_path / f"{f}.h5ad")
        a = read(tmp_path / f"{f}.h5ad")
        assert list(a.obs_names) == list(ref.obs_names) and list(a.var_names) == list(ref.var_names)
        assert a.obs.equals(ref.obs) and a.shape == ref.shape
        assert np.array_equal(a.obsm["spatial_3D"], ref.obsm["spatial_3D"])
        assert not np.array_equal(np.asarray(a.X), np.asarray(ref.X))
        assert s["spec"] == ref_side["spec"] and s["inputs"]["rng_streams"] == ref_side["inputs"]["rng_streams"]
        assert s["construction"]["composition"] == ref_side["construction"]["composition"]
        assert s["inputs"]["seed"] == SEED and s["output"]["check"]["ok"]


# ---------------------------------------------------------------------------------------------------------- refusals


@pytest.mark.parametrize("form", ["c1-s1x1.5", "c1-s1x2.0"])
def test_c1_on_the_larger_exponent_bases_is_refused(world, tmp_path, form):
    with pytest.raises(bk.Refused, match="x = 1 bases only"):
        build(world, form, tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


@pytest.mark.parametrize("form", ["c3-s1", "C1-s1", "c1", "identity", "c2-s1x1.25", "c2-s3", "g-s1", "m-s1", "s-s1",
                                  "p-s1", "", None, "1.5", "c2-s1x1.5,c1-s1"])
def test_other_forms_are_refused(world, tmp_path, form):
    with pytest.raises(bk.Refused, match="form"):
        build(world, form, tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists() and not (tmp_path / "x.json").exists()


def test_the_fixed_forms():
    assert bk.FORMS == ("identity-s1", "identity-s1x1.5", "identity-s1x2.0", "identity-s2", "c1-s1", "c1-s2",
                        "c2-s1", "c2-s1x1.5", "c2-s1x2.0", "c2-s2")
    assert {k: v[:2] for k, v in bk.BASES.items()} == {"s1": ("s1", 1.0), "s1x1.5": ("s1", 1.5),
                                                       "s1x2.0": ("s1", 2.0), "s2": ("s2", 1.0)}
    assert bk.BASES["s1"][2] == "3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2"
    assert bk.KINDS == {"identity": 0.0, "c1": 1.0, "c2": 1.0}


def test_the_table_pins_are_the_committed_decision_files_record(monkeypatch):
    monkeypatch.undo()  # the real pins, not the synthetic ones
    dec = json.loads((T3 / "results" / "gata4" / "gata4_decisions.json").read_text())
    rec = dec["submissions"]["parameters"]["factor_files"]
    assert {k: rec[k]["sha256"] for k in (bk.TABLE_M, bk.TABLE_G)} == bk.TABLE_SHA256


def test_a_table_other_than_the_pinned_one_is_refused(world, tmp_path, monkeypatch):
    for name in (bk.TABLE_M, bk.TABLE_G):
        with monkeypatch.context() as m:
            m.setattr(bk, "TABLE_SHA256", {**world.pins, name: "0" * 64})
            with pytest.raises(bk.Refused, match="pinned table"):
                build(world, "c1-s1", tmp_path / f"{name}.h5ad")
        assert not (tmp_path / f"{name}.h5ad").exists()
    # a table changed after the run recorded it: the frozen loader refuses it first
    run = tmp_path / "run"
    shutil.copytree(world.run, run)
    p = run / "factors" / f"{bk.TABLE_G}.csv"
    p.write_text(p.read_text().replace("G25,2.0", "G25,3.0"))
    with pytest.raises(bk.Refused, match="differs from the one recorded"):
        build(world, "c2-s1", tmp_path / "changed.h5ad", run_dir=run)
    assert not (tmp_path / "changed.h5ad").exists()


def test_a_gene_absent_from_one_table_is_refused(world, tmp_path, monkeypatch):
    """No form fills a gene missing from one table: the frozen loader refuses a table whose genes differ from the
    reference's, so a combination is only ever made over the same genes in the same order."""
    run = tmp_path / "run"
    genes = {bk.TABLE_G: GENES[:-1]}
    write_run(run, world.wt, {bk.TABLE_M: PHI_M, bk.TABLE_G: PHI_G[:-1]}, genes)
    monkeypatch.setattr(bk, "TABLE_SHA256", {k: t3d.sha256(run / "factors" / f"{k}.csv") for k in bk.TABLE_SHA256})
    with pytest.raises(bk.Refused, match="genes differ"):
        build(world, "c1-s1", tmp_path / "absent.h5ad", run_dir=run)
    assert not (tmp_path / "absent.h5ad").exists()


def test_changed_frozen_files_are_refused(world, tmp_path, monkeypatch):
    for name, attr in (("t3_decompose.py", "T3D_PATH"), ("rules.yaml", "RULES_PATH"),
                       ("build_submission.py", "BUILDER_PATH")):
        changed = tmp_path / name
        shutil.copy(getattr(bs, attr), changed)
        with changed.open("a") as fh:
            fh.write("\n# changed\n")
        with pytest.raises(bk.Refused, match=name.replace(".", r"\.")):
            bk.check_frozen({**bk.frozen_paths(), name: changed})
        with monkeypatch.context() as m:
            m.setattr(bs, attr, changed)
            with pytest.raises(bk.Refused, match="frozen file"):
                build(world, "c1-s1", tmp_path / f"{name}.h5ad")
        assert not (tmp_path / f"{name}.h5ad").exists()
    for field, what in (("script_sha256", "script"), ("rules_sha256", "rules")):
        dec = tx._decision(tmp_path / f"dec_{what}.json", **{field: "0" * 64})
        with pytest.raises(bk.Refused, match=what):
            bk.build_form(dec, "c1-s1", world.wt, world.board, N, tmp_path / f"{what}.h5ad", run_dir=world.run,
                          seed=SEED, threads=2)
        assert not (tmp_path / f"{what}.h5ad").exists()


def test_unclean_tree_is_refused(world, tmp_path, own_repo):
    (own_repo / "stray.py").write_text("x = 1\n")
    with pytest.raises(bk.Refused, match="not clean"):
        build(world, "c2-s2", tmp_path / "a.h5ad")
    assert not (tmp_path / "a.h5ad").exists()


def test_frozen_builder_refusals_still_apply(world, tmp_path):
    out = tmp_path / "x.h5ad"
    build(world, "c1-s1", out)
    with pytest.raises(bk.Refused, match="exists"):
        build(world, "c1-s1", out)
    with pytest.raises(bk.Refused, match="seed"):
        bk.build_form(world.dec, "c1-s1", world.wt, world.board, N, tmp_path / "s.h5ad", run_dir=world.run, seed=0,
                      threads=2)


# ------------------------------------------------------------------------------------------------ the frozen bytes


def test_the_frozen_files_are_the_pinned_bytes_and_the_exponent_builders():
    assert bk.FROZEN_SHA256 == bx.FROZEN_SHA256
    assert bk.check_frozen() == bk.FROZEN_SHA256
    for name, sha in bk.FROZEN_SHA256.items():
        assert hashlib.sha256((T3 / name).read_bytes()).hexdigest() == sha


# --------------------------------------------------------------------------------------------------------------- CLI


def test_cli_takes_a_fixed_form_and_no_free_number(world, tmp_path):
    from typer.testing import CliRunner

    from vec.io import check_submission

    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(world.d / "boards"), "--n-cells", str(N), "--run-dir", str(world.run), "--seed", str(SEED),
            "--threads", "2"]
    runner = CliRunner()
    res = runner.invoke(bk.app, [*base, "--form", "c1-s1x2.0", "--out", str(tmp_path / "a.h5ad")])
    assert res.exit_code == 2 and "refused" in res.output and not (tmp_path / "a.h5ad").exists()
    res = runner.invoke(bk.app, [*base, "--form", "c2-s1", "--out", str(tmp_path / "b.h5ad")])
    assert res.exit_code == 0, res.output
    assert check_submission(tmp_path / "b.h5ad", world.board).ok
    side = json.loads((tmp_path / "b.json").read_text())
    assert side["combined_build"]["form"] == "c2-s1" and side["branch"] == "barred" and side["submission"] == "s1"
    for extra in (["--exponent", "1.5"], ["--weight", "0.5"], ["--gata4-weight", "0.5"], ["--c", "1"],
                  ["--table", str(world.run / "factors" / "gata4_pooled.csv")], ["--branch", "allowed"],
                  ["--submission", "s2"], ["--depth"]):
        out = tmp_path / f"{extra[0].strip('-')}.h5ad"
        res = runner.invoke(bk.app, [*base, "--form", "c1-s1", *extra, "--out", str(out)])
        assert res.exit_code != 0 and not out.exists(), extra
