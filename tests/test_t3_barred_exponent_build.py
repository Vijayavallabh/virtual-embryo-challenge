"""Tests for experiments/t3_decomp/build_barred_exponent.py: the barred s1 with its Mab21l2 factor at a named exponent.

Offline, synthetic and fast, in the style of experiments/t3_decomp/test_build_submission.py. No competition data: a
small synthetic WT file stands in for WT E8.75, a synthetic `transfer_mab21l2_all` factor for the run's factors/, and
the decision file carries only what the barred branch reads (D7m) plus the frozen hashes.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
T3 = ROOT / "experiments" / "t3_decomp"
_spec = importlib.util.spec_from_file_location("build_barred_exponent", T3 / "build_barred_exponent.py")
bx = importlib.util.module_from_spec(_spec)
sys.modules["build_barred_exponent"] = bx
_spec.loader.exec_module(bx)
bs = bx.bs
t3d = bs.t3d
RULES = yaml.safe_load((T3 / "rules.yaml").read_text())
SEED = 1_000_001
N = 2500
G = 40
GENES = [f"G{j}" for j in range(G)]
# (celltype, cm_celltype, share), as the frozen builder's tests: every depleted scenario state is in the vocabulary
CELLS = [("V-CM", "OFT/RV-CM", 0.10), ("V-CM", "vCM1", 0.14), ("aSHF", "Unknown", 0.08), ("pSHF", "Unknown", 0.06),
         ("aPHM", "Unknown", 0.08), ("pPHM", "Unknown", 0.08), ("NCC", "Unknown", 0.12), ("SE", "Unknown", 0.12),
         ("Foregut", "Unknown", 0.08), ("Peri", "Unknown", 0.07), ("Endo", "Unknown", 0.07)]
M7M = 0.7
# the synthetic barred factor: up 2, down 0.5, and both clip ends (4 and 0.25) as transfer_factors writes them
PHI = np.ones(G)
PHI[0:10], PHI[10:15], PHI[15:20], PHI[20:25] = 2.0, 0.5, 4.0, 0.25
# phi ** x by hand: 2^1.5 = 2 sqrt 2, 0.5^1.5 = 1 / (2 sqrt 2), 4^1.5 = 8, 0.25^1.5 = 1/8; at 2: 4, 1/4, 16, 1/16
HAND = {1.5: (2.8284271247461903, 0.3535533905932738, 8.0, 0.125), 2.0: (4.0, 0.25, 16.0, 0.0625)}
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES", "GIT_CEILING_DIRECTORIES")


def _git_env() -> dict:
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}


def _git(repo: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                    "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, env=_git_env())


@pytest.fixture(autouse=True)
def own_repo(tmp_path_factory, monkeypatch):
    """Built files are stamped, and the tree checked, against a clean throwaway repository, never the real checkout."""
    from vec import provenance

    repo = tmp_path_factory.mktemp("repo")
    _git(repo, "init", "-q")
    _git(repo, "commit", "-q", "--allow-empty", "-m", "init")
    monkeypatch.setattr(provenance, "REPO_ROOT", repo)
    provenance._record_start(repo)
    return repo


def _write_wt(path: Path, rng, n: int = 900) -> None:
    import anndata as ad
    from scipy import sparse

    from vec.prep import log1p_cp10k

    prof = rng.lognormal(0, 1.0, size=(len(CELLS), G))
    w = np.array([c[2] for c in CELLS])
    lab = rng.choice(len(CELLS), n, p=w / w.sum())
    C = np.vstack([rng.poisson(prof[i] / prof[i].sum() * rng.lognormal(np.log(150), 0.3)) for i in lab])
    C = C.astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [CELLS[i][0] for i in lab], "cm_celltype": [CELLS[i][1] for i in lab]},
                       index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 400, n), rng.uniform(0, 400, n), rng.integers(0, 12, n) * 10.0]
    a.write_h5ad(path)


def _decision(path: Path, **over) -> Path:
    rules = {"D7m": {"outcome": True, "m": M7M, "n_counted": 1}}
    dec = {"preset": "gata4", "rules_sha256": t3d.sha256(bs.RULES_PATH), "script_sha256": t3d.sha256(bs.T3D_PATH),
           "protocol_as_registered": True, "seeds": list(range(10)), "rules": rules, **over}
    path.write_text(json.dumps(dec))
    return path


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    from types import SimpleNamespace

    from vec.boards import Board

    d = tmp_path_factory.mktemp("t3exp")
    wt = d / "wt.h5ad"
    _write_wt(wt, np.random.default_rng(0))
    run = d / "run"
    fdir = run / "factors"
    fdir.mkdir(parents=True)
    p = fdir / "transfer_mab21l2_all.csv"
    pd.DataFrame({"gene": GENES, "factor": PHI}).to_csv(p, index=False)
    (fdir / "factors.json").write_text(json.dumps({
        "preset": "gata4", "files": {"transfer_mab21l2_all": {"file": str(p), "sha256": t3d.sha256(p)}},
        "parameters": {}, "inputs": {"wt": {"path": "wt", "sha256": t3d.sha256(wt)}}}))
    board = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, min_cells=200, max_cells=3000)
    boards_dir = d / "boards"
    boards_dir.mkdir()
    (boards_dir / "T3__synthetic.genes.txt").write_text("\n".join(GENES) + "\n")
    (boards_dir / "index.json").write_text(json.dumps({"T3:synthetic": {
        "task": "T3", "genes_file": "T3__synthetic.genes.txt", "needs_coords": True, "min_cells": 200,
        "max_cells": 3000}}))
    dec = _decision(d / "gata4_decisions.json")
    return SimpleNamespace(d=d, wt=wt, run=run, board=board, boards_dir=boards_dir, dec=dec)


def build_x(world, x, out: Path, **kw) -> dict:
    return bx.build_exponent(world.dec, x, world.wt, world.board, N, out, run_dir=world.run, seed=SEED, threads=2,
                             **kw)


def build_frozen(world, out: Path) -> dict:
    return bs.build_submission(world.dec, "barred", "s1", world.wt, world.board, N, out, run_dir=world.run, seed=SEED,
                               threads=2)


def read(path: Path):
    import anndata as ad

    return ad.read_h5ad(path)


# ---------------------------------------------------------------------------------------------------- identity


def test_exponent_one_is_the_frozen_barred_s1_byte_for_byte(world, tmp_path):
    frozen = build_frozen(world, tmp_path / "frozen.h5ad")
    side = build_x(world, "1.0", tmp_path / "x1.h5ad")
    assert (tmp_path / "frozen.h5ad").read_bytes() == (tmp_path / "x1.h5ad").read_bytes()
    assert side["output"]["sha256"] == frozen["output"]["sha256"] == t3d.sha256(tmp_path / "x1.h5ad")
    assert side["spec"] == frozen["spec"] and side["spec"]["factor"] == {"source": "mab21l2", "exponent": 1.0}
    assert side["construction"] == frozen["construction"]
    assert side["inputs"]["rng_streams"] == frozen["inputs"]["rng_streams"] == ["submission:barred:s1",
                                                                                 "submission:barred:s1:assemble"]
    # the record is the frozen builder's plus one block; the stamp is written for the file as the frozen builder does
    on_disk = json.loads((tmp_path / "x1.json").read_text())
    assert set(on_disk) == set(frozen) | {"exponent_build"}
    eb = on_disk["exponent_build"]
    assert eb["script"] == "experiments/t3_decomp/build_barred_exponent.py" and eb["exponent"] == 1.0
    assert eb["script_sha256"] == t3d.sha256(T3 / "build_barred_exponent.py")
    assert eb["frozen_files_sha256"] == bx.FROZEN_SHA256 and eb["tree_at_start"]["clean"] is True
    assert not any("cell" in k.lower() for k in eb)
    stamp = json.loads((tmp_path / "x1.h5ad.stamp.json").read_text())
    assert stamp["sha256"] == side["output"]["sha256"] and stamp["clean"] is True
    # the frozen builder is left as it was
    assert bs.barred_table.__module__ == "build_submission" and bs.barred_table.__name__ == "barred_table"


# ------------------------------------------------------------------------------------------- same cells, new expression


def test_other_exponents_keep_the_cells_and_composition_and_change_only_expression(world, tmp_path):
    sides = {x: build_x(world, x, tmp_path / f"x{x}.h5ad") for x in (1.0, 1.5, 2.0)}
    files = {x: read(tmp_path / f"x{x}.h5ad") for x in sides}
    base, base_side = files[1.0], sides[1.0]
    for x in (1.5, 2.0):
        a, s = files[x], sides[x]
        assert list(a.obs_names) == list(base.obs_names) and list(a.var_names) == list(base.var_names)
        assert np.array_equal(a.obsm["spatial_3D"], base.obsm["spatial_3D"])
        assert a.obs.equals(base.obs) and a.shape == base.shape
        assert not np.array_equal(np.asarray(a.X), np.asarray(base.X))
        assert s["spec"]["composition"] == base_side["spec"]["composition"]
        assert s["construction"]["composition"] == base_side["construction"]["composition"]
        assert s["spec"]["factor"] == {"source": "mab21l2", "exponent": x}
        assert s["construction"]["factor_exponent"] == x and s["exponent_build"]["exponent"] == x
        assert s["inputs"]["seed"] == SEED and s["inputs"]["rng_streams"] == base_side["inputs"]["rng_streams"]
        assert s["inputs"]["factor"]["sha256"] == base_side["inputs"]["factor"]["sha256"]
        assert s["output"]["check"]["ok"]
        # the spec differs from the frozen one in the exponent alone
        assert {**s["spec"], "factor": base_side["spec"]["factor"]} == base_side["spec"]
    assert not np.array_equal(np.asarray(files[1.5].X), np.asarray(files[2.0].X))


def test_the_factor_applied_is_the_frozen_factor_to_the_power(world, tmp_path, monkeypatch):
    seen: list = []
    original = t3d.apply_gene_factor

    def spy(C, phi, P, rng):
        seen.append(np.array(phi, dtype=np.float64))
        return original(C, phi, P, rng)

    monkeypatch.setattr(t3d, "apply_gene_factor", spy)
    for x in (1.0, 1.5, 2.0):
        seen.clear()
        build_x(world, x, tmp_path / f"p{x}.h5ad")
        assert len(seen) == 1
        phi = seen[0]
        if x == 1.0:
            assert np.array_equal(phi, PHI)
            continue
        up, down, top, bottom = HAND[x]
        assert np.allclose(phi[0:10], up, rtol=1e-15, atol=0) and np.allclose(phi[10:15], down, rtol=1e-15, atol=0)
        assert np.allclose(phi[15:20], top, rtol=1e-15, atol=0)        # beyond the clip: not clipped again
        assert np.allclose(phi[20:25], bottom, rtol=1e-15, atol=0)
        assert np.array_equal(phi[25:], np.ones(G - 25))
        assert np.allclose(phi, PHI ** x, rtol=1e-15, atol=0)


# ---------------------------------------------------------------------------------------------------------- refusals


@pytest.mark.parametrize("x", ["0.5", "0.75", "1.25", "1.4999", "2.5", "3", "0", "-1", "nan", "inf", "abc", None])
def test_other_exponents_are_refused(world, tmp_path, x):
    with pytest.raises(bx.Refused, match="exponent"):
        build_x(world, x, tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists() and not (tmp_path / "x.json").exists()


def test_accepted_exponents_parse():
    assert [bx.parse_exponent(s) for s in ("1", "1.0", "1.5", "2", "2.0")] == [1.0, 1.0, 1.5, 2.0, 2.0]
    assert bx.EXPONENTS == (1.0, 1.5, 2.0)


def test_changed_frozen_files_are_refused(world, tmp_path, monkeypatch):
    # a changed copy of each frozen file, named where the frozen builder names it
    for name, attr in (("t3_decompose.py", "T3D_PATH"), ("rules.yaml", "RULES_PATH"),
                       ("build_submission.py", "BUILDER_PATH")):
        changed = tmp_path / name
        shutil.copy(getattr(bs, attr), changed)
        with changed.open("a") as fh:
            fh.write("\n# changed\n")
        with pytest.raises(bx.Refused, match=name.replace(".", r"\.")):
            bx.check_frozen({**bx.frozen_paths(), name: changed})
        with monkeypatch.context() as m:
            m.setattr(bs, attr, changed)
            with pytest.raises(bx.Refused, match="frozen file"):
                build_x(world, "1.5", tmp_path / f"{name}.h5ad")
        assert not (tmp_path / f"{name}.h5ad").exists()
    # a decision file made by other code or rules: the frozen builder's own check refuses, and nothing overrides it
    for field, what in (("script_sha256", "script"), ("rules_sha256", "rules")):
        dec = _decision(tmp_path / f"dec_{what}.json", **{field: "0" * 64})
        with pytest.raises(bx.Refused, match=what):
            bx.build_exponent(dec, "1.5", world.wt, world.board, N, tmp_path / f"{what}.h5ad", run_dir=world.run,
                              seed=SEED, threads=2)
        assert not (tmp_path / f"{what}.h5ad").exists()


def test_unclean_tree_is_refused(world, tmp_path, own_repo):
    (own_repo / "runs").mkdir()
    (own_repo / "runs" / "output.txt").write_text("an output folder file is allowed\n")
    assert bx.check_clean()["clean"] is True
    (own_repo / "stray.py").write_text("x = 1\n")                     # untracked outside results/ and runs/
    with pytest.raises(bx.Refused, match="not clean"):
        build_x(world, "1.5", tmp_path / "a.h5ad")
    _git(own_repo, "add", "stray.py")
    _git(own_repo, "commit", "-q", "-m", "stray")
    assert bx.check_clean()["clean"] is True
    (own_repo / "stray.py").write_text("x = 2\n")                     # a tracked change
    with pytest.raises(bx.Refused, match="not clean"):
        build_x(world, "1.5", tmp_path / "b.h5ad")
    assert not (tmp_path / "a.h5ad").exists() and not (tmp_path / "b.h5ad").exists()


def test_frozen_builder_refusals_still_apply(world, tmp_path):
    out = tmp_path / "x.h5ad"
    build_x(world, "2.0", out)
    with pytest.raises(bx.Refused, match="exists"):
        build_x(world, "2.0", out)
    with pytest.raises(bx.Refused, match="seed"):
        bx.build_exponent(world.dec, "2.0", world.wt, world.board, N, tmp_path / "s.h5ad", run_dir=world.run,
                          seed=0, threads=2)
    unreg = _decision(tmp_path / "unreg.json", protocol_as_registered=False)
    with pytest.raises(bx.Refused, match="registered protocol"):
        bx.build_exponent(unreg, "2.0", world.wt, world.board, N, tmp_path / "u.h5ad", run_dir=world.run, seed=SEED,
                          threads=2)


# ------------------------------------------------------------------------------------------------ the frozen bytes


def test_the_frozen_files_are_the_pinned_bytes():
    assert bx.check_frozen() == bx.FROZEN_SHA256
    for name, sha in bx.FROZEN_SHA256.items():
        assert hashlib.sha256((T3 / name).read_bytes()).hexdigest() == sha


def test_the_frozen_files_are_mains_bytes():
    for name, sha in bx.FROZEN_SHA256.items():
        res = subprocess.run(["git", "-C", str(ROOT), "show", f"main:experiments/t3_decomp/{name}"],
                             capture_output=True, env=_git_env())
        if res.returncode != 0:
            pytest.skip("branch main is not available in this checkout")
        assert hashlib.sha256(res.stdout).hexdigest() == sha, name


# --------------------------------------------------------------------------------------------------------------- CLI


def test_cli(world, tmp_path):
    from typer.testing import CliRunner

    from vec.io import check_submission

    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(world.boards_dir), "--n-cells", str(N), "--run-dir", str(world.run), "--seed", str(SEED),
            "--threads", "2"]
    runner = CliRunner()
    res = runner.invoke(bx.app, [*base, "--exponent", "0.5", "--out", str(tmp_path / "a.h5ad")])
    assert res.exit_code == 2 and "refused" in res.output and not (tmp_path / "a.h5ad").exists()
    res = runner.invoke(bx.app, [*base, "--exponent", "1.5", "--out", str(tmp_path / "b.h5ad")])
    assert res.exit_code == 0, res.output
    assert check_submission(tmp_path / "b.h5ad", world.board).ok
    side = json.loads((tmp_path / "b.json").read_text())
    assert side["exponent_build"]["exponent"] == 1.5 and side["branch"] == "barred" and side["submission"] == "s1"
    res = runner.invoke(bx.app, [*base, "--exponent", "1.5", "--branch", "refused", "--out", str(tmp_path / "c.h5ad")])
    assert res.exit_code != 0 and not (tmp_path / "c.h5ad").exists()   # the branch and slot are fixed
