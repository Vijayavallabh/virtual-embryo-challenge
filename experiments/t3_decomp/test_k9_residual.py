"""Tests for tries K9a and K9b (experiments/t3_decomp/k9_residual.py), on synthetic data only: no competition file is
opened, nothing is scored, and no prediction is compared with any target."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness
from vec.adoption import SEEDS
from vec.boards import Board
from vec.prep import load_stage, log1p_cp10k

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402
import k9_residual as k9  # noqa: E402

SOURCE = HERE / "k9_residual.py"
G = 16
CLIP = (0.25, 4.0)


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


# ----------------------------------------------------------------------------------------------- constants and tasks
def test_the_fixed_constants():
    assert k9.RES_EXPONENT == 0.5 == k7.EXPONENT
    assert k9.PERM_SEED == 20261113 and k9.STREAM == harness.K4_STREAM
    assert (k9.MIN_GAIN, k9.MIN_AHEAD, k9.X2_FLOOR) == (0.3, 8, -0.3)
    assert k9.ARMS == ("raw", "raw_res", "raw_res_perm")
    assert tuple(harness._t3_head().load_rules()["components"]["factor"]["clip"]) == CLIP


def test_the_tasks_are_k4s_folds_with_k7s_series_added():
    assert set(k9.TASKS) == {"x1", "x2", "r12"} and set(k9.FOLDS.values()) <= set(harness.K4_FOLDS)
    series = {ref: {f"T2/{s}.h5ad" for s in k7.FULL_PANEL if f"T2/{s}.h5ad" != ref} for ref in
              ("T2/E8.75.h5ad", "T2/E9.5.h5ad")}
    x1, x2, r12 = (k9.TASKS[n] for n in ("x1", "x2", "r12"))
    assert x1.reference == "T2/E8.75.h5ad" and x1.target == "T3/E8.75_gata4_ko_2.h5ad"
    assert (x1.inputs["ko_train"], x1.inputs["wt_train"]) == ("T3/E9.5_mab21l2_ko.h5ad", "T2/E9.5.h5ad")
    assert "T3/E8.75_gata4_ko_1.h5ad" in x1.forbidden
    assert x2.reference == "T2/E9.5.h5ad" and x2.target == "T3/E9.5_mab21l2_ko.h5ad"
    assert (x2.inputs["ko_train"], x2.inputs["ko_train_2"], x2.inputs["wt_train"]) == (
        "T3/E8.75_gata4_ko_1.h5ad", "T3/E8.75_gata4_ko_2.h5ad", "T2/E8.75.h5ad")
    assert r12.reference == "T2/E8.75.h5ad" and r12.target == "T3/E8.75_gata4_ko_2.h5ad" and r12.replicate
    assert (r12.inputs["ko_train"], r12.inputs["wt_train"]) == ("T3/E8.75_gata4_ko_1.h5ad", "T2/E8.75.h5ad")
    for name, t in k9.TASKS.items():
        assert harness.check_task(t) is t
        fold = harness.K4_FOLDS[k9.FOLDS[name]]
        assert tuple(t.inputs[r] for r in fold["ko_roles"]) == fold["source_files"]
        assert t.inputs["wt_train"] == fold["wt_file"]
        s = {f for r, f in t.inputs.items() if r.startswith(k9.SERIES_PREFIX)}
        assert s == series[t.reference]
        assert {f for r, f in t.inputs.items() if not r.startswith(k9.SERIES_PREFIX)} == set(
            harness.K4_FOLDS[k9.FOLDS[name]]["source_files"]) | {fold["wt_file"], t.reference}
        assert not t.counted and t.test_board == "T3:bcatenin" and t.target not in t.inputs.values()
        assert harness.TaskView.of(t).reference_role == "wt"
    assert all(c.applies(t) for c in k9.candidates()[1:] for t in k9.TASKS.values())
    assert not any(c.applies(t) for c in k9.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)


def test_the_source_opens_no_file_beyond_the_tasks_inputs():
    src = SOURCE.read_text()
    for call in ("read_h5ad", "load_stage", "read_stage", "h5py", "np.load", "open(", "boards/", "runs/", "data/"):
        assert call not in src, call
    # the only .h5ad path the script spells is the series one, built from k7.FULL_PANEL
    assert re.findall(r'"[^"]*\.h5ad"', src) == ['"T2/{s}.h5ad"']
    assert "k7.FULL_PANEL" in src


# ----------------------------------------------------------------------------------------------- the factors
def test_the_combined_factor_is_phi_raw_times_rho_clipped_per_group():
    rng = np.random.default_rng(0)
    phi_raw = np.clip(np.exp(rng.normal(0, 0.8, G)), *CLIP)
    pool = np.exp(rng.normal(0, 0.3, G))
    per = {"A": pool * np.exp(rng.normal(0, 1.5, G)), "B": pool * np.exp(rng.normal(0, 0.2, G))}
    rho = k9.residual_factors(per, pool)
    for g in per:
        assert np.allclose(rho[g], np.sqrt(per[g] / pool))
        want = np.clip(phi_raw * rho[g], *CLIP)
        assert np.array_equal(k9.combined_factor(phi_raw, rho[g], CLIP), want)
        assert want.min() >= CLIP[0] and want.max() <= CLIP[1]
    assert (np.clip(phi_raw * rho["A"], *CLIP) != phi_raw * rho["A"]).any()  # the clip binds somewhere
    assert np.array_equal(k9.combined_factor(phi_raw, None, CLIP), phi_raw)  # no per-group factor: rho = 1
    drawn = np.array(["B", "A", "Unknown", "C", "A", "B", "C"])
    classes = k9.factor_classes(phi_raw, rho, drawn, CLIP)
    got = {tuple(rows): fac for fac, rows in classes}
    assert np.array_equal(got[(1, 4)], np.clip(phi_raw * rho["A"], *CLIP))
    assert np.array_equal(got[(0, 5)], np.clip(phi_raw * rho["B"], *CLIP))
    assert np.array_equal(got[(2, 3, 6)], phi_raw)  # "C" and "Unknown": no per-group factor, rho = 1, one class


def test_the_residual_is_one_without_exponent_or_without_a_group_part(monkeypatch):
    rng = np.random.default_rng(1)
    pool = np.exp(rng.normal(0, 0.3, G))
    per = {"A": pool * 2.0, "B": pool.copy()}
    assert np.all(k9.residual_factors(per, pool)["B"] == 1.0)
    monkeypatch.setattr(k9, "RES_EXPONENT", 0.0)
    assert all(np.all(v == 1.0) for v in k9.residual_factors(per, pool).values())


def test_the_control_uses_one_permutation_for_every_group_and_keeps_each_deciles_values():
    rng = np.random.default_rng(2)
    n = 50
    basis = rng.gamma(2.0, size=n)
    rho = {g: np.exp(rng.normal(0, 0.4, n)) for g in "ABC"}
    perm = k9.permuted_residuals(rho, basis)
    idx = harness.t3_decile_permutation(np.arange(n, dtype=np.float64), basis, k9.PERM_SEED).astype(int)
    assert not np.array_equal(idx, np.arange(n))
    for g in rho:
        assert np.array_equal(perm[g], rho[g][idx])  # the same permutation for every group
        for genes in np.array_split(np.argsort(basis, kind="stable"), 10):
            assert np.array_equal(np.sort(perm[g][genes]), np.sort(rho[g][genes]))  # each decile keeps its values


def test_the_cosine_is_scale_free_and_weighted_by_the_reference_cells():
    rng = np.random.default_rng(3)
    phi_raw = np.exp(rng.normal(0, 0.5, G))
    rho = {"A": phi_raw ** 0.3, "B": phi_raw ** 0.1}
    assert k9.residual_cosine(phi_raw, rho, {"A": 300, "B": 900}) == pytest.approx(1.0)
    opp = {"A": phi_raw ** -0.5}
    assert k9.residual_cosine(phi_raw, opp, {"A": 300}) == pytest.approx(-1.0)
    mixed = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    m = 0.25 * np.log(mixed["A"]) + 0.75 * np.log(mixed["B"])
    a = np.log(phi_raw)
    assert k9.residual_cosine(phi_raw, mixed, {"A": 100, "B": 300}) == pytest.approx(
        a @ m / np.linalg.norm(a) / np.linalg.norm(m))
    assert k9.residual_cosine(phi_raw, {}, {}) is None


# ----------------------------------------------------------------------------------------------- the arms
def _reference(bs, rng, n=120):
    C = rng.poisson(4.0, size=(n, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, 3, n)
    P = ((C + 0.5) / (C + 0.5).sum(1, keepdims=True)).astype(np.float32)
    pool = bs.t3d.Pool(C, lab, P, 3)
    V = SimpleNamespace(L=3, p=np.array([0.5, 0.3, 0.2]))
    return bs.Reference(preset="gata4", rules={"components": {"stage": {"min_cells": 20}}}, genes=list(range(G)), V=V,
                        pool=pool, coords=rng.normal(size=(n, 3)))


def _spy(bs, monkeypatch):
    calls, real = [], bs.t3d.apply_gene_factor

    def spy(C, phi, P, rng):
        calls.append((C.copy(), np.array(phi)))
        return real(C, phi, P, rng)

    monkeypatch.setattr(bs.t3d, "apply_gene_factor", spy)
    return calls


def test_every_arm_draws_the_same_cells_and_applies_each_groups_factor(bs, monkeypatch):
    rng = np.random.default_rng(4)
    ref = _reference(bs, rng)
    n_pool = ref.pool.C.shape[0]
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, n_pool)]
    phi_raw = np.clip(np.exp(rng.normal(0, 0.6, G)), *CLIP)
    pool = np.exp(rng.normal(0, 0.3, G))
    rho = k9.residual_factors({"A": pool * np.exp(rng.normal(0, 0.8, G)), "B": pool * np.exp(rng.normal(0, 0.8, G))},
                              pool)
    rho_perm = k9.permuted_residuals(rho, rng.gamma(2.0, size=G))
    calls = _spy(bs, monkeypatch)
    built = {arm: k9.arm_built(bs, ref, 200, 7, ref.V.p, arm, phi_raw, r, group, CLIP)
             for arm, r in (("raw", None), ("raw_res", rho), ("raw_res_perm", rho_perm))}
    rows = built["raw"].row
    assert all(np.array_equal(b.row, rows) and np.array_equal(b.lab, built["raw"].lab) for b in built.values())
    assert not np.array_equal(built["raw"].X, built["raw_res"].X)
    assert not np.array_equal(built["raw_res"].X, built["raw_res_perm"].X)
    # raw: one call on all drawn cells with phi_raw; the residual arms: one call per group on the same cells
    assert len(calls) == 1 + 3 + 3
    assert np.array_equal(calls[0][0], ref.pool.C[rows]) and np.array_equal(calls[0][1], phi_raw)
    drawn = group[rows]
    for k, r in ((1, rho), (4, rho_perm)):
        for (C, phi), g in zip(calls[k:k + 3], ("A", "B", "Unknown"), strict=True):
            assert np.array_equal(C, ref.pool.C[rows][drawn == g])
            assert np.array_equal(phi, k9.combined_factor(phi_raw, r.get(g), CLIP))
    again = k9.arm_built(bs, ref, 200, 7, ref.V.p, "raw_res", phi_raw, rho, group, CLIP)
    assert np.array_equal(again.X, built["raw_res"].X)  # a seed's cells are reproducible


@pytest.mark.parametrize("case", ["exponent_zero", "per_equals_pool"])
def test_without_a_residual_raw_res_is_raw_bit_for_bit(bs, monkeypatch, case):
    rng = np.random.default_rng(5)
    ref = _reference(bs, rng)
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, ref.pool.C.shape[0])]
    phi_raw = np.clip(np.exp(rng.normal(0, 0.6, G)), *CLIP)
    pool = np.exp(rng.normal(0, 0.3, G))
    if case == "exponent_zero":
        monkeypatch.setattr(k9, "RES_EXPONENT", 0.0)
        per = {"A": pool * np.exp(rng.normal(0, 0.8, G)), "B": pool * 3.0}
    else:
        per = {"A": pool.copy(), "B": pool.copy()}
    rho = k9.residual_factors(per, pool)
    rho_perm = k9.permuted_residuals(rho, rng.gamma(2.0, size=G))
    raw = k9.arm_built(bs, ref, 200, 9, ref.V.p, "raw", phi_raw, None, group, CLIP)
    for r in (rho, rho_perm):
        res = k9.arm_built(bs, ref, 200, 9, ref.V.p, "raw_res", phi_raw, r, group, CLIP)
        assert np.array_equal(res.X, raw.X) and np.array_equal(res.row, raw.row) and np.array_equal(res.lab, raw.lab)


# ----------------------------------------------------------------------------------------------- end to end, synthetic
E875 = [("V-CM", "OFT/RV-CM", 0.15), ("V-CM", "vCM1", 0.15), ("aSHF", "Unknown", 0.10), ("Endo", "Unknown", 0.08),
        ("NCC", "Unknown", 0.14), ("SE", "Unknown", 0.13), ("Foregut", "Unknown", 0.12),
        ("Neural Tube", "Unknown", 0.13)]
E95 = [("V-CM", None, 0.30), ("A-CM", None, 0.10), ("SE", None, 0.15), ("Foregut", None, 0.15), ("NCC", None, 0.20),
       ("Endo", None, 0.10)]
GENES = [f"G{j}" for j in range(40)]


def _label(cell) -> str:
    return cell[1] if cell[1] not in (None, "Unknown") else cell[0]


def _cells(cells, prof, n, rng, depth, factor=None, name="c") -> ad.AnnData:
    w = np.array([c[2] for c in cells])
    lab = rng.choice(len(cells), n, p=w / w.sum())
    C = np.asarray([rng.poisson(prof[_label(cells[i])] * (1.0 if factor is None else factor)
                                / (prof[_label(cells[i])] * (1.0 if factor is None else factor)).sum()
                                * rng.lognormal(np.log(depth), 0.3)) for i in lab], dtype=np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [cells[i][0] for i in lab]}, index=[f"{name}_{i}" for i in range(n)])
    if any(c[1] is not None for c in cells):
        obs["cm_celltype"] = [cells[i][1] or "Unknown" for i in lab]
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 500, n), rng.uniform(0, 500, n), rng.integers(0, 20, n) * 10.0]
    return a


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Synthetic files in task x1's roles (the Mab21l2 -> Gata4rep2 fold plus the wild-type series), written and loaded
    as the harness loads them."""
    d = tmp_path_factory.mktemp("k9")
    rng = np.random.default_rng(0)
    prof = {k: rng.lognormal(0, 1.0, len(GENES)) for k in sorted({_label(c) for c in E875 + E95})}
    t = k9.TASKS["x1"]
    files = {"wt": _cells(E875, prof, 700, rng, 150, rng.lognormal(0, 0.3, len(GENES)), name="ref"),
             "wt_train": _cells(E95, prof, 900, rng, 80, name="e95"),
             "ko_train": _cells(E95, prof, 700, rng, 70, rng.lognormal(0, 0.4, len(GENES)), name="mab")}
    for r in t.inputs:
        if r.startswith(k9.SERIES_PREFIX):
            files[r] = _cells(E875, prof, 500, rng, 120, rng.lognormal(0, 0.1, len(GENES)), name=r)
    files["s_E9.5"] = files["wt_train"]
    loaded = {}
    for r, a in files.items():
        p = d / f"{r}.h5ad"
        a.write_h5ad(p)
        loaded[r] = load_stage(p, genes=GENES, drop_duplicates=True)
    return loaded


def _ctx(task, inputs, seed=harness.BUILD_SEED_OFFSET + 3):
    board = Board(key=f"pseudo:{task.key}", task="T3", genes=list(GENES), needs_coords=True, split="pseudo",
                  target_ko="gata4", min_cells=100, max_cells=900, truth_cells=120, ref_cells=150)
    return harness.BuildContext(task=harness.TaskView.of(task), board=harness._candidate_board(board), inputs=inputs,
                                seed=seed)


def test_end_to_end_the_raw_arm_is_t3_factor_raw_and_the_residual_moves_it(world, monkeypatch):
    monkeypatch.setattr(k7, "MIN_CELLS", 20)
    t = k9.TASKS["x1"]
    k9._PREP.clear()
    harness._K4_PREP.clear()
    ctx = _ctx(t, world)
    with harness.guarded_reads() as g:
        out = {arm: k9._build(ctx, arm) for arm in k9.ARMS}
    assert not g["refused"]
    prep = k9._PREP["fold"]
    rec = prep["record"]
    assert rec["groups_with_a_factor"] >= 2 and -1.0 <= rec["cosine_log_phi_raw_mean_log_rho"] <= 1.0
    assert rec["series"] == sorted(s for s in k7.FULL_PANEL if s != "E8.75")
    assert json.loads(json.dumps(rec)) == rec  # summaries only, plain JSON
    fold = {x.key: x for x in harness.load_tasks()["T3:bcatenin"]}[k9.FOLDS["x1"]]
    k4ctx = _ctx(fold, {r: world[r] for r in fold.inputs})
    with harness.guarded_reads():
        want = harness.CANDIDATES["t3_factor_raw"].build(k4ctx)
    assert np.array_equal(harness._K4_PREP["fold"]["raw"], prep["phi_raw"])
    assert np.array_equal(out["raw"].X, want.X) and np.array_equal(out["raw"].coords, want.coords)
    assert not np.array_equal(out["raw_res"].X, out["raw"].X)
    k9._PREP.clear()
    monkeypatch.setattr(k9, "RES_EXPONENT", 0.0)
    with harness.guarded_reads():
        none = k9._build(ctx, "raw_res")
    assert np.array_equal(none.X, out["raw"].X) and np.array_equal(none.coords, out["raw"].coords)
    k9._PREP.clear()


# ----------------------------------------------------------------------------------------------- the reading
def _result(key: str, scores: dict[str, list[float]], clean: bool = True, commit: str = "abc") -> dict:
    rows = [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5 + (0.01 if a == k9.MAIN else 0.0), "severity": 0.6} for a in scores},
             "t3": {a: {"severity": {"state": "scored" if a == k9.MAIN else "r2_gate"}} for a in scores}}
            for i, s in enumerate(SEEDS)]
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
            "factor_record": {"groups_with_a_factor": 9, "cosine_log_phi_raw_mean_log_rho": 0.4}}


def _scores(main: float, control: float = 50.0, raw: float = 52.0) -> dict[str, list[float]]:
    n = len(SEEDS)
    return {harness.FLOOR: [50.0] * n, k9.BASE: [raw] * n, k9.MAIN: [main] * n, k9.CONTROL: [control] * n}


def _read(d: Path, **per_task) -> dict:
    d.mkdir()
    for name, sc in per_task.items():
        (d / f"{name}.json").write_text(json.dumps(_result(k9.TASKS[name].key, sc)))
    return k9.reading(d)


def test_the_reading_needs_x1_and_r12_against_both_and_x2_not_below(tmp_path):
    yes = _read(tmp_path / "yes", x1=_scores(54.0), r12=_scores(53.0), x2=_scores(51.8))
    assert not yes["void"] and yes["k9a_supported"] and yes["verdict"].startswith("SUPPORTED")
    assert yes["tasks"]["x1"]["paired"]["raw_res - raw"]["mean"] == pytest.approx(2.0)
    assert yes["tasks"]["x1"]["skill_diff_raw_res_minus_raw"] == pytest.approx({"de_score": 0.01, "severity": 0.0})
    assert yes["tasks"]["x2"]["severity_states"][k9.MAIN] == {"scored": len(SEEDS)}
    assert yes["tasks"]["x2"]["counted_fold"] and not yes["tasks"]["r12"]["counted_fold"]
    txt = k9.text(yes)
    assert "NOT the adoption rule" in txt and "never counts" in txt and "WT E9.5" in txt and "flattered" in txt
    assert "groups with a per-group factor 9" in txt
    # x2 at exactly -0.3 is not below it; just under is
    assert _read(tmp_path / "edge", x1=_scores(54.0), r12=_scores(53.0), x2=_scores(51.7))["k9a_supported"]
    assert not _read(tmp_path / "x2", x1=_scores(54.0), r12=_scores(53.0), x2=_scores(51.6))["k9a_supported"]
    # r12 failing against raw, or x1 failing against the control, closes it
    assert not _read(tmp_path / "r12", x1=_scores(54.0), r12=_scores(52.2), x2=_scores(53.0))["k9a_supported"]
    no = _read(tmp_path / "ctrl", x1=_scores(54.0, control=53.9), r12=_scores(53.0), x2=_scores(53.0))
    assert not no["k9a_supported"] and no["verdict"].startswith("CLOSED")
    assert not _read(tmp_path / "r12c", x1=_scores(54.0), r12=_scores(53.0, control=53.0),
                     x2=_scores(53.0))["k9a_supported"]
    # x2 gaining does not stand in for r12
    assert not _read(tmp_path / "swap", x1=_scores(54.0), r12=_scores(52.0), x2=_scores(60.0))["k9a_supported"]


def test_ahead_in_eight_seeds_is_required(tmp_path):
    n = len(SEEDS)
    sc = _scores(54.0)
    sc[k9.MAIN] = [70.0] * 7 + [40.0] * (n - 7)
    res = _read(tmp_path / "a", x1=sc, r12=_scores(53.0), x2=_scores(53.0))
    assert not res["k9a_supported"] and res["tasks"]["x1"]["paired"]["raw_res - raw"]["ahead"] == 7


def test_the_reading_is_void_on_missing_dirty_or_mixed_results(tmp_path):
    missing = _read(tmp_path / "m", x1=_scores(54.0), r12=_scores(53.0))
    assert missing["void"] and any("missing task x2" in w for w in missing["why"])
    d = tmp_path / "dirty"
    d.mkdir()
    for name in k9.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k9.TASKS[name].key, _scores(54.0), clean=name != "x2")))
    assert k9.reading(d)["void"]
    d = tmp_path / "commits"
    d.mkdir()
    for i, name in enumerate(k9.TASKS):
        (d / f"{name}.json").write_text(json.dumps(_result(k9.TASKS[name].key, _scores(54.0), commit=f"c{i % 2}")))
    assert k9.reading(d)["void"]


def test_run_and_read_refuse_to_overwrite(tmp_path):
    (tmp_path / "x1.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k9.main(["run", "--task", "x1", "--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])
    runs = tmp_path / "runs"
    _read(runs, x1=_scores(54.0), r12=_scores(53.0), x2=_scores(53.0))
    out = tmp_path / "out"
    k9.main(["read", "--runs", str(runs), "--out", str(out)])
    assert json.loads((out / "reading.json").read_text())["k9a_supported"]
    with pytest.raises(SystemExit, match="exists"):
        k9.main(["read", "--runs", str(runs), "--out", str(out)])
