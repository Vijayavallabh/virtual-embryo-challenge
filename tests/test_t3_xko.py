"""Tries K2a and K2b (experiments/harness/TRIES.md): the harness candidates t3_factor_xko and t3_factor_xko_perm and the
T3 parts their rows' fit checks read, on synthetic data only (no competition file, no score of the fold). The factor is
the one the pre-registered T3 run writes for s1 on the same files (t3_decompose.build_context, submission_parameters,
build_submission.load_factor); the control's permutation stays within deciles, is fixed by its seed and keeps each
decile's values; the build applies the factor through t3_decompose.apply_gene_factor with the kNN-smoothed top-up
profile, at the reference's own composition, never through vec.core.apply_factors."""

from __future__ import annotations

import json
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import spearmanr

from vec import core, harness
from vec.boards import GROUPS, Board
from vec.evaluate import Stage, panel_metrics, t3_response_detail
from vec.prep import load_stage, log1p_cp10k

G = 40
GENES = [f"G{j}" for j in range(G)]
# WT E8.75-like, in the gata4 preset's combined labels (celltype, cm_celltype, share): strict non-Mesp1 labels (NCC,
# SE, Foregut, Neural Tube) and Mesp1-derived ones
E875 = [("V-CM", "OFT/RV-CM", 0.11), ("V-CM", "vCM1", 0.11), ("aSHF", "Unknown", 0.08), ("pSHF", "Unknown", 0.07),
        ("aPHM", "Unknown", 0.08), ("Peri", "Unknown", 0.07), ("Endo", "Unknown", 0.06), ("NCC", "Unknown", 0.12),
        ("SE", "Unknown", 0.10), ("Foregut", "Unknown", 0.10), ("Neural Tube", "Unknown", 0.10)]
# WT E9.5-like, the fold's reference (the mab21l2 preset: celltype labels)
E95 = [("V-CM", None, 0.25), ("A-CM", None, 0.10), ("Peri", None, 0.10), ("SE", None, 0.10), ("Foregut", None, 0.10),
       ("NCC", None, 0.15), ("ST", None, 0.10), ("Endo", None, 0.10)]
# sha256 recorded by the T3 run that wrote s1's factor (experiments/t3_decomp/results/gata4/gata4_decisions.json:
# script_sha256 and rules_sha256)
S1_RUN = {"t3_decompose.py": "3e671c3bb5ff57bba5ceec80ea36fcc739b1f111399a1e2854400b5c9f5133d2",
          "rules.yaml": "cd2d94aa6e202529a2eead2a7ad60b312d9a2fd9ff0bd989d9320662de750c13"}
FOLD = "T3/Gata4->Mab21l2"


def _label(cell) -> str:
    return cell[1] if cell[1] not in (None, "Unknown") else cell[0]


def _cells(cells, prof, n, rng, depth, factor=None, condition=None, name="c") -> ad.AnnData:
    """A MERFISH-like AnnData: labels drawn from the shares, Poisson counts around per-label profiles (times a per-gene
    factor), log1p CP10k X, a counts layer and section-like positions."""
    w = np.array([c[2] for c in cells])
    lab = rng.choice(len(cells), n, p=w / w.sum())
    rows = []
    for i in lab:
        p = prof[_label(cells[i])] * (1.0 if factor is None else factor)
        rows.append(rng.poisson(p / p.sum() * rng.lognormal(np.log(depth), 0.3)))
    C = np.asarray(rows, dtype=np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [cells[i][0] for i in lab]}, index=[f"{name}_{i}" for i in range(n)])
    if any(c[1] is not None for c in cells):
        obs["cm_celltype"] = [cells[i][1] or "Unknown" for i in lab]
    if condition is not None:
        obs["condition"] = condition
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 500, n), rng.uniform(0, 500, n), rng.integers(0, 20, n) * 10.0]
    return a


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Synthetic files in the fold's roles, and the inputs as the harness loads them (load_stage, duplicates
    dropped). The knockout replicates are shallower than their WT (so the harmoniser refits at their depths), share a
    per-gene processing factor, carry fewer V-CM cells, and replicate 1 holds littermate wild types."""
    d = tmp_path_factory.mktemp("k2")
    rng = np.random.default_rng(0)
    prof = {k: rng.lognormal(0, 1.0, G) for k in sorted({_label(c) for c in E875 + E95})}
    ko_factor = rng.lognormal(0, 0.5, G)
    ko = [(c, m, s * (0.5 if c == "V-CM" else 1.0)) for c, m, s in E875]
    files = {"wt_train": _cells(E875, prof, 900, rng, 150, name="e875"),
             "ko_train": _cells(ko, prof, 600, rng, 90, ko_factor, name="rep1",
                                condition=np.where(rng.random(600) < 0.2, "WT", "KO")),
             "ko_train_2": _cells(ko, prof, 600, rng, 95, ko_factor, name="rep2", condition=np.full(600, "KO")),
             "wt": _cells(E95, prof, 1500, rng, 75, name="e95")}
    paths = {role: d / f"{role}.h5ad" for role in files}
    for role, a in files.items():
        a.write_h5ad(paths[role])
    inputs = {role: load_stage(p, genes=GENES, drop_duplicates=True) for role, p in paths.items()}
    # a Mab21l2-like knockout of the E9.5 reference, the scored target of the synthetic fold
    target = _cells([(c, m, s * (2.0 if c == "ST" else 1.0)) for c, m, s in E95], prof, 2400, rng, 90,
                    rng.lognormal(0, 0.4, G), name="ko95")
    return SimpleNamespace(paths=paths, inputs=inputs, ko_factor=ko_factor,
                           target=Stage(X=np.asarray(target.X.toarray(), np.float32),
                                        coords=np.asarray(target.obsm["spatial_3D"], np.float32)))


@pytest.fixture(scope="module")
def t3():
    bs = harness._t3_head()
    return SimpleNamespace(bs=bs, t3d=bs.t3d, rules=bs.load_rules())


def _fold() -> harness.PseudoTask:
    return {t.key: t for t in harness.load_tasks()["T3:bcatenin"]}[FOLD]


def _board(max_cells: int = 900) -> Board:
    return Board(key=f"pseudo:{FOLD}", task="T3", genes=list(GENES), needs_coords=True, split="pseudo",
                 target_ko="mab21l2", min_cells=100, max_cells=max_cells, truth_cells=120, ref_cells=150)


# ------------------------------------------------------------------------------------------------ the factor
def test_k2_uses_the_t3_code_that_wrote_s1s_factor(t3):
    for name, sha in S1_RUN.items():
        assert t3.t3d.sha256(harness.T3_DIR / name) == sha, name
    fc = t3.rules["components"]["factor"]
    assert (fc["eps"], list(fc["clip"]), fc["min_cells"]) == (1e-5, [0.25, 4.0], 20)  # the rows' eps and clip
    assert t3.rules["presets"]["gata4"]["label_mode"] == "combined"
    assert t3.rules["presets"]["mab21l2"]["label_mode"] == "celltype"


@pytest.mark.filterwarnings("ignore:.*obs names of the model's training cells")  # the stage pools are the WT file
def test_k2a_factor_is_the_factor_the_t3_run_writes_for_s1(world, t3, tmp_path):
    p = world.paths
    paths = {"wt": str(p["wt_train"]), "targets": [("rep1", str(p["ko_train"])), ("rep2", str(p["ko_train_2"]))],
             "stage:E9.5": str(p["wt_train"]), "stage:E8.25_late": str(p["wt_train"])}
    ctx, _ = t3.t3d.build_context("gata4", paths, t3.rules, n_jobs=1, threads=1, out=tmp_path)
    params = t3.t3d.submission_parameters(ctx, tmp_path)
    s1_phi, _ = t3.bs.load_factor(tmp_path, "gata4_pooled", ctx.genes, {"submissions": {"parameters": params}})
    assert ctx.targets[0].littermate is not None and ctx.targets[1].littermate is None
    inp = world.inputs
    before = {r: (list(a.obs.columns), sorted(a.uns)) for r, a in inp.items()}
    phi = harness.t3_gata4_factor(inp["wt_train"], [inp["ko_train"], inp["ko_train_2"]], t3.rules)
    assert ctx.genes == GENES and np.array_equal(phi, s1_phi)
    assert before == {r: (list(a.obs.columns), sorted(a.uns)) for r, a in inp.items()}  # inputs left as they were
    # not a trivial factor: it recovers the replicates' processing factor, and pooling both replicates matters
    assert spearmanr(np.log(phi), np.log(world.ko_factor)).correlation > 0.7
    assert not np.array_equal(phi, harness.t3_gata4_factor(inp["wt_train"], [inp["ko_train"]], t3.rules))


# ------------------------------------------------------------------------------------------------ the control
def test_k2b_permutation_stays_within_deciles_and_is_fixed_by_its_seed():
    assert (harness.K2_PERM_SEED, harness.K2_PERM_BINS) == (20261104, 10)  # fixed in K2b's row
    rng = np.random.default_rng(3)
    for n in (500, 103):
        phi, basis = rng.lognormal(0, 0.5, n), rng.gamma(1.0, 1.0, n)
        out = harness.t3_decile_permutation(phi, basis)
        assert np.array_equal(out, harness.t3_decile_permutation(phi, basis, seed=20261104))
        assert not np.array_equal(out, harness.t3_decile_permutation(phi, basis, seed=20261105))
        bins = np.array_split(np.argsort(basis, kind="stable"), 10)
        assert sorted(map(len, bins))[0] >= n // 10 and sorted(map(len, bins))[-1] <= -(-n // 10)
        decile = np.empty(n, int)
        for i, b in enumerate(bins):
            decile[b] = i
            assert np.array_equal(np.sort(out[b]), np.sort(phi[b]))  # each decile keeps its values
        source = {v: g for g, v in enumerate(phi)}
        assert all(decile[source[out[g]]] == decile[g] for g in range(n))  # every gene's value is from its decile
        assert (out != phi).mean() > 0.7
    # ties keep panel order: one flat basis is cut into runs of consecutive genes
    flat = harness.t3_decile_permutation(np.arange(20.0), np.zeros(20), bins=2)
    assert set(flat[:10]) == set(range(10)) and set(flat[10:]) == set(range(10, 20))


def test_k2b_basis_is_mean_released_expression(world):
    a = world.inputs["wt_train"]
    want = np.asarray(a.X.toarray(), np.float64).mean(0)
    assert np.allclose(harness.t3_mean_expression(a), want, rtol=0, atol=1e-12)


# ------------------------------------------------------------------------------------------------ the build
def test_k2_builds_through_t3_decompose_with_the_knn_smoothed_profile(world, t3, monkeypatch):
    assert t3.rules["components"]["smoothing"] == {"n_pcs": 30, "k": 15}  # rules note 5's top-up profile
    calls, real = [], t3.t3d.apply_gene_factor

    def spy(C, phi, P, rng):
        calls.append((C.copy(), np.array(phi), np.array(P)))
        return real(C, phi, P, rng)

    def never(*a, **k):
        raise AssertionError("vec.core.apply_factors used (own-count top-up)")

    monkeypatch.setattr(t3.t3d, "apply_gene_factor", spy)
    monkeypatch.setattr(core, "apply_factors", never)
    harness._K2_PREP.clear()  # the factor and the reference are made inside the guard below
    board = _board()
    ctx = harness.BuildContext(task=harness.TaskView.of(_fold()), board=harness._candidate_board(board),
                               inputs=world.inputs, seed=harness.BUILD_SEED_OFFSET + 3)
    with harness.guarded_reads() as g:
        a = harness.CANDIDATES["t3_factor_xko"].build(ctx)
        b = harness.CANDIDATES["t3_factor_xko_perm"].build(ctx)
    assert not g["refused"] and "experiments/t3_decomp/rules.yaml" in g["log"]
    prep = harness._K2_PREP["fold"]
    (Ca, phia, Pa), (Cb, phib, Pb) = calls
    assert np.array_equal(phia, prep["phi"]) and np.array_equal(phib, prep["phi_perm"])  # exponent 1: phi as is
    assert np.array_equal(Ca, Cb) and np.array_equal(Pa, Pb)  # both arms factor the same drawn cells
    pool, V = prep["ref"].pool, prep["ref"].V
    where = {r.tobytes(): i for i, r in enumerate(pool.C)}
    rows = np.array([where[c.tobytes()] for c in Ca])
    assert np.array_equal(Pa, pool.P[rows])  # the top-up profile is each drawn cell's kNN-smoothed proportions
    assert np.allclose(pool.P, t3.t3d.smoothed_props(pool.C, 30, 15), rtol=0, atol=1e-6)
    assert np.array_equal(np.bincount(pool.lab[rows], minlength=V.L), t3.t3d.alloc(board.max_cells, V.p))
    for s in (a, b):
        assert s.X.shape == (board.max_cells, G) and np.unique(s.coords, axis=0).shape[0] == board.max_cells
    assert not np.array_equal(a.X, b.X)


def test_k2_reference_is_build_submissions_reference(world, t3):
    with t3.bs.thread_limits(1):
        mine = harness.t3_reference(world.inputs["wt"], "mab21l2", t3.rules)
    theirs = t3.bs.load_reference(world.paths["wt"], "mab21l2", t3.rules, threads=1)
    assert mine.genes == theirs.genes and mine.V.labels == theirs.V.labels
    for x, y in ((mine.V.p, theirs.V.p), (mine.V.nm, theirs.V.nm), (mine.pool.C, theirs.pool.C),
                 (mine.pool.lab, theirs.pool.lab), (mine.pool.P, theirs.pool.P), (mine.coords, theirs.coords)):
        assert np.array_equal(x, y)


def test_k2_runs_on_the_gata4_to_mab21l2_fold_only():
    tasks = harness.load_tasks()
    for name in ("t3_factor_xko", "t3_factor_xko_perm"):
        c = harness.CANDIDATES[name]
        assert [t.key for ts in tasks.values() for t in ts if c.applies(t)] == [FOLD]


def _k2_result(gains, seeds=harness.SEEDS):
    cands = ("floor", "copy", "shift_transfer", "t3_factor_xko", "t3_factor_xko_perm")
    rows = [{"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
             "scores": {"floor": 50.0, "copy": 50.5, "shift_transfer": 47.0, "t3_factor_xko_perm": 49.0,
                        "t3_factor_xko": 49.0 + g}} for s, g in zip(seeds, gains, strict=True)]
    prov = {"commit": "abc", "tree_clean": True, "config_sha256": harness._sha256(harness.CONFIG),
            "files_sha256": {"T2/E9.5.h5ad": "f"}}
    return {"task": {"key": FOLD, "test_board": "T3:bcatenin"}, "board": {"key": f"pseudo:{FOLD}"},
            "candidates": list(cands), "rows": rows, "protocol": "bake", "build": "main", "anchors": "computed",
            "frac": 0.1, "seeds": list(seeds), "provenance": prov}


def test_decide_counts_only_the_fold_for_k2a_against_k2b_and_the_baselines():
    configured = harness.load_tasks()["T3:bcatenin"]
    d = harness.decide([_k2_result([0.5] * 8 + [-0.1] * 2)], configured, "t3_factor_xko_perm", "t3_factor_xko")
    assert d["keep"] and d["counted"] == [FOLD] and d["need"] == 1
    assert d["excluded"] == {"T3/Gata4rep1->rep2": "replicate fold", "T3/Mab21l2->Gata4rep2": "change does not apply"}
    assert not harness.decide([_k2_result([0.5] * 7 + [-0.1] * 3)], configured, "t3_factor_xko_perm",
                              "t3_factor_xko")["keep"]  # ahead in 7 of 10 seeds
    for base in ("floor", "copy", "shift_transfer"):
        assert harness.decide([_k2_result([0.5] * 10)], configured, base, "t3_factor_xko")["counted"] == [FOLD]


# ------------------------------------------------------------------------------------------------ the fit checks
def _ko_stages(world):
    ref = world.inputs["wt"]
    return world.target, Stage(X=np.asarray(ref.X.toarray(), np.float32),
                               coords=np.asarray(ref.obsm["spatial_3D"], np.float32))


def test_t3_response_detail_reproduces_the_panel(world):
    target, ref = _ko_stages(world)
    truth, refw = target.take(np.arange(0, 2400, 2)), ref.take(np.arange(0, 1500, 3))
    rng = np.random.default_rng(1)
    dt = truth.X.mean(0) - refw.X.mean(0)
    preds = {"no_response": refw, "anti": Stage(X=np.clip(refw.X - 2.0 * dt, 0, None)),
             "partial": Stage(X=np.clip(refw.X + 0.5 * dt + rng.normal(0, 0.05, (refw.n, G)), 0, None).astype(
                 np.float32)), "noise": Stage(X=np.clip(refw.X + rng.normal(0, 0.3, (refw.n, G)), 0, None).astype(
                     np.float32))}
    states = {}
    for name, pred in preds.items():
        m = panel_metrics("T3", pred, truth, refw)
        det = t3_response_detail(pred.X, truth.X, refw.X, metrics=m)
        assert all(det["matches_panel"].values()), name
        assert len(det["de_direction_terms"]) == G
        assert abs(np.clip(det["de_direction"], -1, 1) - m["de_direction"]) <= 1e-4
        if det["severity"]["state"] == "scored":
            assert m["severity_slope"] == round(float(np.log(det["severity"]["beta"])), 4)
        states[name] = det["severity"]["state"]
    assert states["no_response"] == "no_response" and states["anti"] == "beta_gate" and states["partial"] == "scored"
    assert not any(t3_response_detail(refw.X, truth.X, refw.X)["de_direction_terms"])


def test_score_seeds_rows_carry_groups_and_t3_parts_for_the_fit_checks(world):
    target, ref = _ko_stages(world)
    cands = [harness.CANDIDATES[n] for n in ("floor", "copy", "t3_factor_xko", "t3_factor_xko_perm")]
    rows = harness.score_seeds(_fold(), _board(), target, ref, world.inputs, cands, seeds=(0, 1))
    for row in rows:
        assert set(row["t3"]) == set(row["groups"]) == set(row["scores"])
        assert all(all(det["matches_panel"].values()) for det in row["t3"].values())
        assert row["t3"]["floor"]["severity"]["state"] == "no_response"
        assert row["n_cells"]["t3_factor_xko"] == row["n_cells"]["t3_factor_xko_perm"] == 900
        for c in row["scores"]:
            total = sum(100.0 * w * row["groups"][c][g] for g, w, _ in GROUPS["T3"])
            assert abs(total - row["scores"][c]) < 1e-9
    fc = harness.t3_fit_checks({"rows": rows}, "t3_factor_xko_perm", "t3_factor_xko")
    margin = np.mean([r["scores"]["t3_factor_xko"] - r["scores"]["t3_factor_xko_perm"] for r in rows])
    ddir = np.mean([r["metrics"]["t3_factor_xko"]["de_direction"] - r["metrics"]["t3_factor_xko_perm"]["de_direction"]
                    for r in rows])
    assert fc["seeds"] == [0, 1] and fc["matches_panel"] and abs(fc["score_margin"] - margin) < 1e-9
    assert abs(fc["de_direction_margin"] - ddir) <= 2e-4 and set(fc["group_margins"]) == {g for g, _, _ in GROUPS["T3"]}
    assert harness.t3_fit_checks({"rows": rows}, "floor", "shift_transfer") is None


def test_fit_checks_read_the_margin_by_group_the_top_genes_and_the_gate():
    def row(seed, base_terms, change_terms, states):
        det = {n: {"de_direction_terms": list(t), "severity": {"state": s, "gated": s != "scored"},
                   "matches_panel": {"de_direction": True, "severity": True}}
               for n, t, s in (("b", base_terms, states[0]), ("c", change_terms, states[1]))}
        groups = {"b": {g: 0.5 for g, _, _ in GROUPS["T3"]},
                  "c": {g: 0.5 + (0.1 if g == "response_direction" else 0.0) for g, _, _ in GROUPS["T3"]}}
        return {"seed": seed, "t3": det, "groups": groups}

    zero, change = np.zeros(100), np.r_[np.full(10, 0.05), np.full(90, 0.001)]
    fc = harness.t3_fit_checks({"rows": [row(0, zero, change, ("scored", "scored")),
                                         row(1, zero, change, ("r2_gate", "scored"))]}, "b", "c")
    assert fc["group_margins"]["response_direction"] == pytest.approx(2.5) and fc["score_margin"] == pytest.approx(2.5)
    assert fc["de_direction_margin"] == pytest.approx(0.59) and fc["top_share"] == pytest.approx(0.5 / 0.59)
    assert fc["severity_gate_differs"] == [1] and fc["severity_states"][1] == {"b": "r2_gate", "c": "scored"}
    neg = harness.t3_fit_checks({"rows": [row(0, change, zero, ("scored", "scored"))]}, "b", "c")
    assert neg["de_direction_margin"] < 0 and neg["top_share"] is None


def test_run_decide_writes_the_fit_checks(tmp_path):
    import importlib.util

    spec = importlib.util.spec_from_file_location("harness_run_for_k2", harness.REPO_ROOT / "experiments" / "harness"
                                                  / "run.py")
    run = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run)
    res = _k2_result([0.5] * 10)
    for r in res["rows"]:
        r["groups"] = {c: {g: 0.5 for g, _, _ in GROUPS["T3"]} for c in res["candidates"]}
        r["t3"] = {c: {"de_direction_terms": [0.0] * 5, "severity": {"state": "scored", "gated": False},
                       "matches_panel": {"de_direction": True, "severity": True}} for c in res["candidates"]}
    (tmp_path / "runs").mkdir()
    (tmp_path / "runs" / "T3_Gata4_Mab21l2.json").write_text(json.dumps(res))
    run.main(["decide", "--runs", str(tmp_path / "runs"), "--board", "T3:bcatenin", "--base", "t3_factor_xko_perm",
              "--change", "t3_factor_xko", "--out", str(tmp_path / "out")])
    dec = json.loads((tmp_path / "out" / "decisions_T3_bcatenin.json").read_text())["decisions"][0]
    assert dec["keep"] and dec["fit_checks"][FOLD]["severity_gate_differs"] == []
    assert "fit checks" in (tmp_path / "out" / "decisions_T3_bcatenin.txt").read_text()
