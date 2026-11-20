"""Tests for experiments/t3_decomp/build_factor_form.py: the barred file with its Mab21l2 factor in one of try K21's
fixed forms. Offline, synthetic and fast: the helpers and throwaway-repository fixture of
tests/test_t3_barred_exponent_build.py, a synthetic `transfer_mab21l2_all` table and synthetic wild-type-only factors.
No competition data."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
from pathlib import Path

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
bf = _load("build_factor_form", T3 / "build_factor_form.py")
bk = sys.modules.get("build_two_knockout") or _load("build_two_knockout", T3 / "build_two_knockout.py")
import factor_forms as ff  # noqa: E402
import k21_factor_forms as k21  # noqa: E402

bx, bs = tx.bx, bf.bs
t3d = bs.t3d
own_repo = tx.own_repo  # builds are stamped against a throwaway repository
SEED, N, G, GENES = tx.SEED, tx.N, tx.G, tx.GENES
# a synthetic barred factor: four genes at the clip, twenty distinct moves, sixteen at 1
PHI_M = np.ones(G)
PHI_M[0:3], PHI_M[3] = 4.0, 0.25
_mag = np.linspace(0.1, 1.2, 20)
PHI_M[4:24] = np.exp(np.where(np.arange(20) % 2 == 0, _mag, -_mag))
K_TEST = 8


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    from types import SimpleNamespace

    from vec.boards import Board

    d = tmp_path_factory.mktemp("t3form")
    wt = d / "wt.h5ad"
    tx._write_wt(wt, np.random.default_rng(0))
    run = d / "run"
    fdir = run / "factors"
    fdir.mkdir(parents=True)
    p = fdir / f"{bf.TABLE_M}.csv"
    pd.DataFrame({"gene": GENES, "factor": PHI_M}).to_csv(p, index=False)
    (fdir / "factors.json").write_text(json.dumps({
        "preset": "gata4", "files": {bf.TABLE_M: {"file": str(p), "sha256": t3d.sha256(p)}}, "parameters": {},
        "inputs": {"wt": {"path": "wt", "sha256": t3d.sha256(wt)}}}))
    board = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, min_cells=200, max_cells=3000)
    dec = tx._decision(d / "gata4_decisions.json")
    return SimpleNamespace(d=d, wt=wt, run=run, board=board, dec=dec, pin=t3d.sha256(p))


@pytest.fixture(autouse=True)
def synthetic(world, monkeypatch):
    """The synthetic table's pin (the real pin is tested below) and a top-K for 40 genes."""
    monkeypatch.setattr(bf, "TABLE_SHA256", world.pin)
    monkeypatch.setattr(ff, "TOP_K", K_TEST)


@pytest.fixture
def wt_factors(monkeypatch):
    """The frozen wt_revert_series and wt_revert_factors stubbed with synthetic wild-type-only factors fitted to the
    pool's rows (two groups and one without a factor of its own)."""
    rng = np.random.default_rng(7)
    per = {"A": np.exp(rng.normal(0, 0.4, G)), "B": np.exp(rng.normal(0, 0.4, G))}
    pooled = np.exp(rng.normal(0, 0.4, G))
    seen = {}

    def series(wt, given=None):
        return {"E8.5": Path("E8.5.h5ad")}

    def factors(wt, ser, genes, pool_counts=None):
        assert list(genes) == GENES and pool_counts is not None
        group = np.array(["A", "B", "Unknown"])[np.arange(len(pool_counts)) % 3]
        seen["group"] = group
        return pooled, per, group, {"rule": "synthetic", "series": {"E8.5": {"path": "E8.5.h5ad", "sha256": "0" * 64}}}

    monkeypatch.setattr(bs, "wt_revert_series", series)
    monkeypatch.setattr(bs, "wt_revert_factors", factors)
    return {"per": per, "pooled": pooled, "seen": seen}


def build(world, form, out: Path, **kw) -> dict:
    return bf.build_form(world.dec, form, world.wt, world.board, N, out, run_dir=kw.pop("run_dir", world.run),
                         seed=SEED, threads=2, **kw)


def read(path: Path):
    import anndata as ad

    return ad.read_h5ad(path)


class Spy:
    def __init__(self, monkeypatch):
        self.seen: list = []
        original = t3d.apply_gene_factor

        def spy(C, phi, P, rng):
            self.seen.append(np.array(phi, dtype=np.float64))
            return original(C, phi, P, rng)

        monkeypatch.setattr(t3d, "apply_gene_factor", spy)


# ---------------------------------------------------------------------------------------------------- identity


@pytest.mark.parametrize("base", ["s1", "s2"])
def test_identity_gives_the_frozen_barred_file_byte_for_byte(world, tmp_path, base):
    frozen = bs.build_submission(world.dec, "barred", base, world.wt, world.board, N, tmp_path / "frozen.h5ad",
                                 run_dir=world.run, seed=SEED, threads=2)
    side = build(world, f"identity-{base}", tmp_path / "id.h5ad")
    assert (tmp_path / "frozen.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()
    assert side["output"]["sha256"] == frozen["output"]["sha256"] and side["construction"] == frozen["construction"]
    assert {k: v for k, v in side["inputs"].items() if k != "factor_note"} == frozen["inputs"]
    fb = json.loads((tmp_path / "id.json").read_text())["form_build"]
    assert fb["form"] == f"identity-{base}" and fb["base"]["reference_sha256"] == bf.BASES[base][2]
    assert fb["script_sha256"] == t3d.sha256(T3 / "build_factor_form.py")
    assert fb["factor_forms_sha256"] == t3d.sha256(T3 / "factor_forms.py")
    assert fb["frozen_files_sha256"] == bf.FROZEN_SHA256 and fb["tree_at_start"]["clean"] is True
    assert fb["table"]["sha256"] == fb["table"]["pinned_sha256"] and fb["factor"]["applied"]["genes_changed"] == 0
    assert not any("cell" in k.lower() for k in fb)
    assert bs.load_factor.__module__ == bs.construct.__module__ == "build_submission"  # left as they were


@pytest.mark.parametrize("x", ["1.5", "2.0"])
def test_identity_gives_the_exponent_builders_file_byte_for_byte(world, tmp_path, x):
    ref = bx.build_exponent(world.dec, x, world.wt, world.board, N, tmp_path / "bx.h5ad", run_dir=world.run,
                            seed=SEED, threads=2)
    side = build(world, f"identity-s1x{x}", tmp_path / "id.h5ad")
    assert (tmp_path / "bx.h5ad").read_bytes() == (tmp_path / "id.h5ad").read_bytes()
    assert side["output"]["sha256"] == ref["output"]["sha256"]


# ------------------------------------------------------------------------- the form applied is the practice arm's


@pytest.mark.parametrize("form", ["top100-s1", "cap90-s1", "top100-s1x2.0", "cap90-s1x1.5", "top100-s2", "cap90-s2"])
def test_the_factor_applied_is_the_practice_arms_at_the_bases_exponent(world, tmp_path, monkeypatch, form):
    spy = Spy(monkeypatch)
    side = build(world, form, tmp_path / "f.h5ad")
    (phi,) = spy.seen
    kind, base = form.split("-", 1)
    x = bf.BASES[base][1]
    arm = k21.arm_factor(kind, PHI_M)  # the very function the practice arm applies
    assert np.array_equal(phi, arm if x == 1.0 else arm ** x)
    # by hand: top-K keeps the four clip genes and the four largest moves (genes 23, 22, 21, 20), 1 elsewhere
    if kind == "top100":
        keep = [0, 1, 2, 3, 20, 21, 22, 23]
        want = np.ones(G)
        want[keep] = PHI_M[keep]
        assert np.array_equal(arm, want)
    else:  # the cap: the 0.9 quantile of |ln phi| off the clip
        off = np.abs(np.log(PHI_M))[4:]
        t = float(np.quantile(off, 0.9))
        assert side["form_build"]["factor"]["transform"]["t"] == t
        assert np.array_equal(arm, np.clip(PHI_M, np.exp(-t), np.exp(t)))
    fb = side["form_build"]
    assert fb["form"] == form and fb["base"]["x"] == x and side["construction"]["factor_exponent"] == 1.0
    assert side["inputs"]["factor_note"] == bf.FACTOR_NOTE and fb["constants"]["top_k"] == K_TEST


@pytest.mark.parametrize("base,forms", [("s1", ["top100-s1", "cap90-s1", "wtprod-s1"]),
                                        ("s2", ["top100-s2", "wtprod-s2"])])
def test_every_form_draws_the_same_cells_as_its_base(world, tmp_path, wt_factors, base, forms):
    ref_side = build(world, f"identity-{base}", tmp_path / "id.h5ad")
    ref = read(tmp_path / "id.h5ad")
    for f in forms:
        s = build(world, f, tmp_path / f"{f}.h5ad")
        a = read(tmp_path / f"{f}.h5ad")
        assert list(a.obs_names) == list(ref.obs_names) and a.obs.equals(ref.obs) and a.shape == ref.shape
        assert np.array_equal(a.obsm["spatial_3D"], ref.obsm["spatial_3D"])
        assert not np.array_equal(np.asarray(a.X), np.asarray(ref.X))
        assert s["spec"] == ref_side["spec"] and s["output"]["check"]["ok"]


@pytest.mark.parametrize("form", ["wtprod-s1", "wtprod-s1x2.0"])
def test_wtprod_moves_by_group_with_phi_m_squared_on_the_frozen_stream(world, tmp_path, monkeypatch, wt_factors, form):
    seen = {}
    original = bs.apply_wt_revert

    def spy(ref, built, phi, per, group, rng):
        seen.update(phi=phi, per=per, group=group, state=rng.bit_generator.state)
        return original(ref, built, phi, per, group, rng)

    monkeypatch.setattr(bs, "apply_wt_revert", spy)
    side = build(world, form, tmp_path / "w.h5ad")
    x = bf.BASES[form.split("-", 1)[1]][1]
    sq = (PHI_M if x == 1.0 else PHI_M ** x) ** 2
    assert np.array_equal(seen["phi"], sq * wt_factors["pooled"])
    assert set(seen["per"]) == {"A", "B"} and all(np.array_equal(seen["per"][g], sq * wt_factors["per"][g])
                                                  for g in seen["per"])
    assert seen["state"] == t3d.stable_rng(SEED, "submission:barred:s1:wt_revert").bit_generator.state
    fb = side["form_build"]
    assert fb["factor"]["wtprod"]["stream"] == "submission:barred:s1:wt_revert"
    assert side["inputs"]["rng_streams"][-1] == "submission:barred:s1:wt_revert"
    assert fb["factor"]["wtprod"]["factor_by_group"] == {"A": "its own factor", "B": "its own factor",
                                                         "Unknown": "the pooled factor"}
    want_max = max(float(np.max((sq * v) ** 0.5)) for v in [*wt_factors["per"].values(), wt_factors["pooled"]])
    assert fb["factor"]["wtprod"]["applied_range"]["max"] == want_max
    assert fb["factor"]["wtprod"]["wt_revert"]["series"]["E8.5"]["sha256"] == "0" * 64
    assert side["construction"]["wt_revert"] is None  # the frozen record: its own --wt-revert was not used


def test_wtprod_refuses_a_base_with_a_step(world, tmp_path, monkeypatch, wt_factors):
    real = bs.barred_table

    def with_depth(dec, rules):
        t = dict(real(dec, rules))
        t["s1"] = {**t["s1"], "depth": True}
        return t

    monkeypatch.setattr(bs, "barred_table", with_depth)
    with pytest.raises(bf.Refused, match="carries a step"):
        build(world, "wtprod-s1", tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


# ---------------------------------------------------------------------------------------------------------- refusals


@pytest.mark.parametrize("form", ["uni_top-s1", "uni_cap-s2", "wtprod_perm-s1", "soft50s-s1", "rand100-s1",
                                  "barred_form-s1"])
def test_controls_and_dropped_forms_are_refused_by_name(world, tmp_path, form):
    with pytest.raises(bf.Refused, match="refused"):
        build(world, form, tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


@pytest.mark.parametrize("form", ["top100", "top100-s3", "TOP100-s1", "top100-s1 ", " cap90-s1", "cap90-s1\n",
                                  "c1-s1", "identity", "", None, "top50-s1", "cap80-s1", "wtprod-s1x1.25"])
def test_other_forms_are_refused(world, tmp_path, form):
    with pytest.raises(bf.Refused, match="form"):
        build(world, form, tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


def test_the_fixed_forms():
    assert bf.FORMS == tuple(f"{k}-{b}" for k in ("identity", "top100", "cap90", "wtprod")
                             for b in ("s1", "s1x1.5", "s1x2.0", "s2"))
    assert bf.parse_form("cap90-s1x1.5") == ("cap90", "s1x1.5")


def test_top100_is_refused_when_the_clip_fills_the_kept_set(world, tmp_path, monkeypatch):
    monkeypatch.setattr(ff, "TOP_K", 4)  # four genes sit at the clip
    with pytest.raises(bf.Refused, match="at the clip"):
        build(world, "top100-s1", tmp_path / "x.h5ad")
    assert not (tmp_path / "x.h5ad").exists()


def test_the_pins_are_the_other_builders_and_the_decision_files(monkeypatch):
    monkeypatch.undo()
    assert bf.FROZEN_SHA256 == bk.FROZEN_SHA256 == bx.FROZEN_SHA256
    assert bf.BASES == bk.BASES and bf.TABLE_SHA256 == bk.TABLE_SHA256[bf.TABLE_M] and bf.TABLE_M == bk.TABLE_M
    dec = json.loads((T3 / "results" / "gata4" / "gata4_decisions.json").read_text())
    assert dec["submissions"]["parameters"]["factor_files"][bf.TABLE_M]["sha256"] == bf.TABLE_SHA256
    assert bf.check_frozen() == bf.FROZEN_SHA256
    for name, sha in bf.FROZEN_SHA256.items():
        assert hashlib.sha256((T3 / name).read_bytes()).hexdigest() == sha
    assert bf.WT_REVERT_STREAM == "wt_revert" and ff.WT_EXPONENT == 0.5


def test_a_table_other_than_the_pinned_one_is_refused(world, tmp_path, monkeypatch):
    monkeypatch.setattr(bf, "TABLE_SHA256", "0" * 64)
    with pytest.raises(bf.Refused, match="pinned table"):
        build(world, "cap90-s1", tmp_path / "x.h5ad")
    run = tmp_path / "run"
    shutil.copytree(world.run, run)
    p = run / "factors" / f"{bf.TABLE_M}.csv"
    p.write_text(p.read_text().replace("G0,4.0", "G0,3.0"))
    monkeypatch.setattr(bf, "TABLE_SHA256", world.pin)
    with pytest.raises(bf.Refused, match="differs from the one recorded"):
        build(world, "cap90-s1", tmp_path / "y.h5ad", run_dir=run)
    assert not (tmp_path / "x.h5ad").exists() and not (tmp_path / "y.h5ad").exists()


def test_changed_frozen_files_are_refused(world, tmp_path, monkeypatch):
    for name, attr in (("t3_decompose.py", "T3D_PATH"), ("rules.yaml", "RULES_PATH"),
                       ("build_submission.py", "BUILDER_PATH")):
        changed = tmp_path / name
        shutil.copy(getattr(bs, attr), changed)
        with changed.open("a") as fh:
            fh.write("\n# changed\n")
        with monkeypatch.context() as m:
            m.setattr(bs, attr, changed)
            with pytest.raises(bf.Refused, match="frozen file"):
                build(world, "top100-s1", tmp_path / f"{name}.h5ad")
        assert not (tmp_path / f"{name}.h5ad").exists()


def test_unclean_tree_is_refused(world, tmp_path, own_repo):
    (own_repo / "stray.py").write_text("x = 1\n")
    with pytest.raises(bf.Refused, match="not clean"):
        build(world, "cap90-s2", tmp_path / "a.h5ad")
    assert not (tmp_path / "a.h5ad").exists()


# --------------------------------------------------------------------------------------------------------------- CLI


def _boards(tmp_path: Path) -> Path:
    boards = tmp_path / "boards"
    boards.mkdir()
    (boards / "T3__synthetic.genes.txt").write_text("\n".join(GENES) + "\n")
    (boards / "index.json").write_text(json.dumps({"T3:synthetic": {
        "task": "T3", "genes_file": "T3__synthetic.genes.txt", "needs_coords": True, "min_cells": 200,
        "max_cells": N}}))
    return boards


def test_cli_takes_a_fixed_form_and_no_free_number(world, tmp_path):
    from typer.testing import CliRunner

    from vec.io import check_submission

    base = ["build", str(world.dec), "--wt", str(world.wt), "--board", "T3:synthetic", "--boards-dir",
            str(_boards(tmp_path)), "--run-dir", str(world.run), "--threads", "2"]
    runner = CliRunner()
    res = runner.invoke(bf.app, [*base, "--form", "uni_cap-s1", "--out", str(tmp_path / "a.h5ad")])
    assert res.exit_code == 2 and "refused" in res.output and not (tmp_path / "a.h5ad").exists()
    res = runner.invoke(bf.app, [*base, "--form", "cap90-s1", "--out", str(tmp_path / "b.h5ad")])
    assert res.exit_code == 0, res.output
    assert check_submission(tmp_path / "b.h5ad", world.board).ok
    for extra in (["--seed", "1000002"], ["--n-cells", str(N)], ["--wt-series", "E8.5=x.h5ad"], ["--k", "50"],
                  ["--q", "0.8"], ["--exponent", "1.5"], ["--branch", "allowed"]):
        out = tmp_path / f"{extra[0].strip('-')}.h5ad"
        res = runner.invoke(bf.app, [*base, "--form", "cap90-s1", *extra, "--out", str(out)])
        assert res.exit_code != 0 and not out.exists(), extra
