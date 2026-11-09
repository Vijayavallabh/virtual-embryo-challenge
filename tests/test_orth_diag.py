"""Tries U3 and U3d (experiments/harness/TRIES.md, experiments/t3_decomp/orth_diag.py) on synthetic data only: the
removal of the stage axis and gene_factor's rescale; the decision rule at its boundaries; the depth matching; the N11
rerun's equivalence with experiments/gate/t3_backup_sensitivity.py where the measured cells equal N11's; the gate before
any scoring; and a synthetic end-to-end run through the real commands, whose rows pair with the T3 script's own scored
variant. No competition file, score of a project target or board file is read."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.boards import RANKED_METRICS
from vec.evaluate import Stage, aggregate, panel_metrics, skills_from_anchors, split_half
from vec.prep import load_stage, log1p_cp10k
from vec.provenance import read_stamp, stamp_path

ROOT = Path(__file__).resolve().parents[1]
G = 40
GENES = [f"G{j}" for j in range(G)]
# WT E8.75-like labels in the gata4 preset's combined form (celltype, cm_celltype, share); NCC, SE, Foregut and Neural
# Tube are the strict non-Mesp1 labels the axis is measured on
CELLS = [("V-CM", "OFT/RV-CM", 0.11), ("V-CM", "vCM1", 0.11), ("aSHF", "Unknown", 0.08), ("pSHF", "Unknown", 0.07),
         ("aPHM", "Unknown", 0.08), ("Peri", "Unknown", 0.07), ("Endo", "Unknown", 0.06), ("NCC", "Unknown", 0.13),
         ("SE", "Unknown", 0.11), ("Foregut", "Unknown", 0.09), ("Neural Tube", "Unknown", 0.09)]
FAST_DRAWS = 3000


@pytest.fixture(scope="module")
def od():
    spec = importlib.util.spec_from_file_location("t3_orth_diag", ROOT / "experiments" / "t3_decomp" / "orth_diag.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _label(cell) -> str:
    return cell[1] if cell[1] != "Unknown" else cell[0]


def _write(path, prof, n, rng, depth, factor=None):
    """A MERFISH-like file: labels from CELLS, Poisson counts around per-label profiles times a per-gene factor, log1p
    CP10k X, a counts layer and section-like positions."""
    w = np.array([c[2] for c in CELLS], dtype=np.float64)
    lab = rng.choice(len(CELLS), n, p=w / w.sum())
    f = np.ones(G) if factor is None else factor
    C = np.array([rng.poisson(prof[_label(CELLS[i])] * f / (prof[_label(CELLS[i])] * f).sum()
                              * rng.lognormal(np.log(depth), 0.25)) for i in lab], dtype=np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [CELLS[i][0] for i in lab], "cm_celltype": [CELLS[i][1] for i in lab]},
                       index=[f"{path.stem}_{i}" for i in range(n)])
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 500, n), rng.uniform(0, 500, n), rng.integers(0, 20, n) * 10.0]
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)


@pytest.fixture(scope="module")
def world(tmp_path_factory, od):
    """Synthetic released files in the run's roles, and what the T3 run committed about them (the hashed factors, the
    decision file's hashes, the meta's label count, D5/D6's correlation, the records rows). The factors are
    t3_decompose.transfer_factors', so the gate can reproduce them. D5/D6's values are computed here by t3_checks.py
    d5d6's formulas, not through orth_diag's helpers."""
    d = tmp_path_factory.mktemp("orth")
    data = d / "raw"
    rng = np.random.default_rng(0)
    prof = {k: rng.lognormal(0, 1.0, G) for k in sorted({_label(c) for c in CELLS})}
    stage = rng.normal(0, 0.5, G)                                      # a per-gene developmental direction
    mab = np.exp(0.7 * stage + rng.normal(0, 0.4, G))                  # the Mab21l2 factor leans on the stage direction
    g4 = np.exp(0.2 * stage + rng.normal(0, 0.4, G))                   # the Gata4 factor less
    _write(data / "T2" / "E8.75.h5ad", prof, 2400, rng, 130)
    _write(data / "T2" / "E9.5.h5ad", prof, 2400, rng, 80, np.exp(0.6 * stage))
    _write(data / "T2" / "E8.25_late.h5ad", prof, 2400, rng, 200, np.exp(-0.4 * stage))
    _write(data / "T2" / "E8.5.h5ad", prof, 2000, rng, 170, np.exp(-0.1 * stage))
    _write(data / "T3" / "E9.5_mab21l2_ko.h5ad", prof, 2400, rng, 90, mab)
    _write(data / "T3" / "E8.75_gata4_ko_1.h5ad", prof, 1600, rng, 65, g4)
    _write(data / "T3" / "E8.75_gata4_ko_2.h5ad", prof, 1600, rng, 75, g4)
    t3d, rules = od.t3d, od.bs.load_rules()
    tmp = d / "tmp"
    tmp.mkdir()
    with od.bs.thread_limits(1):
        ctx, _ = t3d.build_context("gata4", od.context_paths(data, with_targets=False), rules, n_jobs=1, threads=1,
                                   out=tmp)
        meta_rec: dict = {}
        factors = t3d.transfer_factors(ctx, str(data / "T3" / "E9.5_mab21l2_ko.h5ad"), meta_rec)
    run = d / "run"
    (run / "factors").mkdir(parents=True)
    files = {}
    for k, v in (("transfer_mab21l2_all", factors["mab21l2_all"]), ("transfer_mab21l2_vs95", factors["mab21l2_vs95"])):
        p = run / "factors" / f"{k}.csv"
        pd.DataFrame({"gene": ctx.genes, "factor": v}).to_csv(p, index=False)
        files[k] = {"file": str(p), "sha256": t3d.sha256(p)}
    decision = d / "gata4_decisions.json"
    decision.write_text(json.dumps({"preset": "gata4", "rules_sha256": t3d.sha256(t3d.RULES_PATH),
                                    "script_sha256": t3d.sha256(Path(t3d.__file__)), "rules": {},
                                    "submissions": {"parameters": {"factor_files": files}}}))
    meta = d / "gata4_meta.json"
    meta.write_text(json.dumps({"transfer_ko": {"labels_used_all": meta_rec["transfer_ko"]["labels_used_all"]}}))
    # D5's correlation, by t3_checks.py's formulas on this context: the axis is gene_factor(E9.5 against E8.25_late) on
    # the strict non-Mesp1 labels with >= 100 cells in both stages, centred, and the factor the centred log
    fc, V = rules["components"]["factor"], ctx.V
    mu = {sp.name: t3d.props_type_means(sp.pool.C, sp.pool.lab, V.L) for sp in ctx.stages}
    cnt = {sp.name: np.asarray(sp.count) for sp in ctx.stages}
    mask = np.asarray(V.nm, dtype=bool) & (cnt["E9.5"] >= 100) & (cnt["E8.25_late"] >= 100)
    phi, _ = t3d.gene_factor(mu["E9.5"], mu["E8.25_late"], V.p, mask, fc["eps"], tuple(fc["clip"]))
    d1 = np.log(phi) - np.log(phi).mean()
    all_ = np.log(factors["mab21l2_all"]) - np.log(factors["mab21l2_all"]).mean()
    v95 = np.log(factors["mab21l2_vs95"]) - np.log(factors["mab21l2_vs95"]).mean()
    corr = float(np.corrcoef(all_, d1)[0, 1])
    dd = d / "d5d6.json"
    dd.write_text(json.dumps({"D5": {"primary": {"corr_fM_d": corr}, "secondary": {"corr_fM_d": corr},
                                     "context": {"corr_d_primary_with_log_all_over_vs95":
                                                 float(np.corrcoef(d1, all_ - v95)[0, 1])}}}))
    # the committed rows of B and vs95 for the two replicates, ten seeds (values on the 0-100 scale, B's the higher)
    r2 = np.random.default_rng(1)
    rec = [{"target": t, "seed": s, "variant": v, "t3": float(loc + r2.normal(0, 1)), "sev_floor": False}
           for v, loc in (("transfer:mab21l2_all", 60.0), ("transfer:mab21l2_vs95", 50.0))
           for t in ("rep1", "rep2") for s in range(10)]
    records = d / "gata4_records.csv"
    pd.DataFrame(rec).to_csv(records, index=False)
    return SimpleNamespace(root=d, data=data, run=run, decision=decision, meta=meta, d5d6=dd, records=records, ctx=ctx,
                           d1=d1, corr=corr, rules=rules, phi_all=factors["mab21l2_all"],
                           labels_used=meta_rec["transfer_ko"]["labels_used_all"])


@pytest.fixture(scope="module")
def ctx_targets(world, od):
    """The T3 script's context with the two synthetic Gata4 replicates, as `score` builds it."""
    tmp = world.root / "tmp_targets"
    tmp.mkdir()
    with od.bs.thread_limits(1):
        ctx, _ = od.t3d.build_context("gata4", od.context_paths(world.data, with_targets=True), world.rules, n_jobs=1,
                                      threads=1, out=tmp)
    ref = od.reference_counts(ctx, world.data / "T2" / "E8.75.h5ad", ctx.preset["label_mode"])
    return SimpleNamespace(ctx=ctx, ref=ref)


@pytest.fixture
def run_args(world, tmp_path):
    out, fdir = tmp_path / "out", tmp_path / "factors"
    common = ["--data", str(world.data), "--run-dir", str(world.run), "--decision", str(world.decision)]
    return SimpleNamespace(out=out, fdir=fdir, axis=["axis", "--out", str(out), *common, "--factor-dir", str(fdir),
                                                     "--d5d6", str(world.d5d6), "--meta", str(world.meta),
                                                     "--threads", "1"],
                           score=["score", "--out", str(out), *common, "--records", str(world.records),
                                  "--workers", "1", "--threads", "1", "--seeds", "0,1", "--code-check"])


# ---------------------------------------------------------------------------------------------------------------------
# The removal of the axis and the rescale


def test_remove_axis_is_the_ols_projection_then_gene_factors_rescale_and_clip(od):
    rng = np.random.default_rng(0)
    d = rng.normal(0, 1, 300)
    d -= d.mean()
    r = rng.normal(0, 0.4, 300)
    r -= r.mean()
    r -= d * (r @ d) / (d @ d)
    v = 0.1 + 0.35 * d + r                                     # log phi_B: a mean, the axis part and an orthogonal part
    den = rng.dirichlet(np.ones(300)) * 0.4
    out = od.remove_axis(np.exp(v), d, den, (1e-6, 1e6))
    assert out["b"] == pytest.approx(0.35)
    lp = np.log(out["phi_perp"])
    want = v - 0.35 * d
    np.testing.assert_allclose(lp - lp.mean(), want - want.mean(), atol=1e-9)
    assert den @ out["phi_perp"] == pytest.approx(den.sum())  # gene_factor's rescale: the weighted total is unchanged
    assert abs(np.corrcoef(lp, d)[0, 1]) < 1e-6
    assert out["corr_before_clip_with_axis"] == pytest.approx(np.corrcoef(lp, d)[0, 1], abs=1e-9)
    assert out["share_of_genes_clipped"] == 0.0
    # the clip is gene_factor's: after the rescale, and only on phi
    tight = od.remove_axis(np.exp(v), d, den, (0.9, 1.1))
    assert tight["phi_perp"].min() >= 0.9 and tight["phi_perp"].max() <= 1.1 and tight["share_of_genes_clipped"] > 0
    scaled = np.exp(want) / (den @ np.exp(want) / den.sum())
    np.testing.assert_allclose(tight["phi_perp"], np.clip(scaled, 0.9, 1.1), atol=1e-12)


def test_remove_axis_equals_gene_factor_rescale_of_the_same_raw_factor(od):
    """The rescale is gene_factor's: dividing by the den-weighted mean of the raw factor, then the clip."""
    t3d = od.t3d
    rng = np.random.default_rng(1)
    mu_wt = rng.dirichlet(np.ones(50), size=4)
    raw = np.exp(rng.normal(0, 0.7, 50))
    w = np.array([0.4, 0.3, 0.2, 0.1])
    mask = np.array([1, 1, 0, 1], dtype=bool)
    # gene_factor on a knockout whose ratio to WT is exactly `raw` gives the rescaled, clipped raw factor
    phi, den = t3d.gene_factor(mu_wt * raw, mu_wt, w, mask, eps=0.0, clip=(0.25, 4.0))
    want = np.clip(raw / (den @ raw / den.sum()), 0.25, 4.0)
    np.testing.assert_allclose(phi, want, atol=1e-12)
    d = rng.normal(0, 1, 50)
    d -= d.mean()
    b = float((np.log(phi) - np.log(phi).mean()) @ d / (d @ d))
    out = od.remove_axis(phi, d, den, (0.25, 4.0))
    assert out["b"] == pytest.approx(b)
    expect = np.clip(np.exp(np.log(phi) - b * d) / (den @ np.exp(np.log(phi) - b * d) / den.sum()), 0.25, 4.0)
    np.testing.assert_allclose(out["phi_perp"], expect, atol=1e-12)


def test_barred_denominator_keeps_the_labels_with_100_cells_on_both_sides(od, tmp_path):
    """The mask of the rescale is t3_decompose.transfer_factors': labels with >= 100 Mab21l2 cells and >= 100 WT cells;
    the factor and the denominator are gene_factor's on that mask."""
    t3d = od.t3d
    rng = np.random.default_rng(9)
    genes = [f"G{j}" for j in range(12)]
    labels = np.array(["A"] * 150 + ["B"] * 100 + ["C"] * 99 + ["D"] * 300 + ["nope"] * 20)
    C = rng.poisson(30, size=(labels.size, 12)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": ["V-CM"] * labels.size}, index=[f"c{i}" for i in range(labels.size)])
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=genes))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 9, labels.size), rng.uniform(0, 9, labels.size),
                                 rng.integers(0, 5, labels.size) * 1.0]
    p = tmp_path / "mab.h5ad"
    a.write_h5ad(p)
    V = SimpleNamespace(labels=["A", "B", "C", "D"], L=4, p=np.array([0.25, 0.25, 0.25, 0.25]))
    wt_mu = rng.dirichlet(np.ones(12), size=4)
    ctx = SimpleNamespace(genes=genes, V=V,
                          wt=SimpleNamespace(count=np.array([500, 500, 500, 90]), mu=wt_mu),
                          rules={"harmonisation": {"seed": 0},
                                 "components": {"factor": {"eps": 1e-5, "clip": [0.25, 4.0]}}},
                          harm=SimpleNamespace(calls=lambda q, seed: (pd.DataFrame({"label_raw": labels}), 7.0)))
    both, phi, den, key = od.barred_denominator(ctx, p)
    assert key == 7.0 and both.tolist() == [True, True, False, False]  # C has 99 cells; D has 90 WT cells
    lab = t3d.calls_to_index(V, pd.DataFrame({"label_raw": labels}))
    mu = t3d.props_type_means(C, lab, 4)
    want_phi, want_den = t3d.gene_factor(mu, wt_mu, V.p, both, 1e-5, (0.25, 4.0))
    np.testing.assert_array_equal(phi, want_phi)
    np.testing.assert_array_equal(den, want_den)


def test_stage_axis_is_t3_checks_formula(od):
    t3d = od.t3d
    rng = np.random.default_rng(2)
    mu_late = rng.dirichlet(np.ones(60), size=4)
    mu_early = rng.dirichlet(np.ones(60), size=4)
    V = SimpleNamespace(nm=np.array([1, 1, 0, 1], dtype=bool), p=np.array([0.25, 0.25, 0.25, 0.25]), L=4)
    cnt = {"a": np.array([100, 99, 500, 300]), "b": np.array([100, 500, 500, 100])}
    fc = {"eps": 1e-5, "clip": [0.25, 4.0]}
    got = od.stage_axis({"a": mu_late, "b": mu_early}, cnt, V, fc, "a", "b")
    mask = V.nm & (cnt["a"] >= 100) & (cnt["b"] >= 100)
    assert mask.tolist() == [True, False, False, True]  # label 1 has 99 cells in one stage; label 2 is not non-Mesp1
    phi, _ = t3d.gene_factor(mu_late, mu_early, V.p, mask, 1e-5, (0.25, 4.0))
    np.testing.assert_allclose(got, np.log(phi) - np.log(phi).mean(), atol=1e-12)
    assert got.mean() == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------------------------------------------------
# The gate and the decision rule


@pytest.mark.parametrize("corr,repro,labels,sd_perp,expected", [
    (0.42, 1e-5, 28, 0.60, (True, True, False, True)),
    (0.4225 + 0.0104, 1e-5, 28, 0.60, (False, True, False, False)),
    (0.4225 - 0.0099, 1e-5, 28, 0.60, (True, True, False, True)),
    (float("nan"), 1e-5, 28, 0.60, (False, True, False, False)),
    (0.42, 0.02, 28, 0.60, (True, False, False, False)),
    (0.42, 1e-5, 27, 0.60, (True, False, False, False)),
    (0.42, 1e-5, 28, 0.645, (True, True, True, True)),   # a no-op: U3 inconclusive; `proceed` for U3d
    (0.42, 1e-5, 28, 0.630, (True, True, False, True)),  # SD change exactly 0.02 is not below it
])
def test_gate_outcome_stops_on_the_checks_and_leaves_a_noop_inconclusive(od, corr, repro, labels, sd_perp, expected):
    out = od.gate_outcome(corr, 0.4225, repro, labels, 28, 0.65, sd_perp)
    assert (out["self_check_passed"], out["barred_reproduced"], out["noop"], out["proceed"]) == expected


def _gate(od, corr=0.42, repro=1e-6, sd_perp=0.5):
    return {"corr_fM_d": corr, "corr_fM_d_committed": 0.4225, "tolerance": 0.01, "reported_only": {
        "corr_all_minus_vs95_d": 0.5, "corr_all_minus_vs95_d_committed": 0.59,
        "secondary_axis": {"axis": "E9.5 vs E8.5", "corr_fM_d": 0.3, "corr_fM_d_committed": 0.3, "scored": False}},
        "barred": {"factor": "transfer_mab21l2_all", "sha256": "x", "labels_in_mask": 28,
                   "labels_in_mask_committed": 28, "max_abs_diff_log_phi": repro, "harmoniser_depth_target": None},
        "removal": {"b": 0.1, "corr_before_clip_with_axis": 0.0, "share_of_genes_clipped": 0.0, "sd_log_phi_B": 0.6,
                    "sd_log_phi_perp": sd_perp, "sd_change": sd_perp - 0.6, "corr_log_phi_perp_with_axis": 0.0,
                    "corr_log_phi_perp_with_barred": 0.9},
        **od.gate_outcome(corr, 0.4225, repro, 28, 28, 0.6, sd_perp)}


def test_axis_text_says_each_verdict(od):
    text = od.axis_text
    res = {"factor_perp": {"file": "f", "sha256": "s"}, "stamp": {"commit": "c", "tree_clean": True, "utc": "u"}}
    assert "STOP and fix the code" in text({"gate": _gate(od, corr=0.5), **res})
    assert "NOT REPRODUCED" in text({"gate": _gate(od, repro=0.5), **res})
    noop = "U3 is inconclusive and B-perp is not scored; U3d still runs"
    assert noop in text({"gate": _gate(od, sd_perp=0.595), **res})
    assert "proceed to score U3 and U3d" in text({"gate": _gate(od), **res})


def test_load_gate_refuses_a_missing_or_failed_gate_or_another_script(od, tmp_path):
    with pytest.raises(SystemExit, match="run `orth_diag.py axis` first"):
        od.load_gate(tmp_path, "sha")
    p = tmp_path / od.AXIS_JSON
    p.write_text(json.dumps({"gate": {"proceed": False, "noop": False}, "script_sha256": "sha"}))
    with pytest.raises(SystemExit, match="did not pass"):
        od.load_gate(tmp_path, "sha")
    p.write_text(json.dumps({"gate": {"proceed": True, "noop": False}, "script_sha256": "sha"}))
    with pytest.raises(SystemExit, match="another version of this script"):
        od.load_gate(tmp_path, "other")
    assert od.load_gate(tmp_path, "sha")["gate"]["proceed"]
    p.write_text(json.dumps({"gate": {"proceed": True, "noop": True}, "script_sha256": "sha", "factor_perp": {}}))
    assert od.load_phi_perp(od.load_gate(tmp_path, "sha"), ["G0"]) is None  # a no-op: B-perp is not loaded or scored


def test_load_phi_perp_checks_the_hash_and_the_gene_order(od, tmp_path):
    p = tmp_path / "perp.csv"
    pd.DataFrame({"gene": ["a", "b"], "factor": [1.5, 0.5]}).to_csv(p, index=False)
    h = od.t3d.sha256(p)
    gate = {"gate": {"noop": False}, "factor_perp": {"file": str(p), "sha256": h}}
    assert od.load_phi_perp(gate, ["a", "b"]).tolist() == [1.5, 0.5]
    with pytest.raises(SystemExit, match="genes differ"):
        od.load_phi_perp(gate, ["b", "a"])
    gate["factor_perp"]["sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="differs from the gate"):
        od.load_phi_perp(gate, ["a", "b"])


@pytest.mark.parametrize("q,floored,outcome", [
    (0.75, 2, "i"), (0.7499, 2, "iii"), (0.75, 3, "iii"), (0.9, 0, "i"), (1.4, 1, "i"),
    (0.5, 0, "ii"), (0.5001, 0, "iii"), (0.9, 5, "ii"), (0.2, 3, "ii"), (-0.2, 0, "ii"),
    (0.75, 4, "iii"), (0.76, 4, "iii"), (float("nan"), 0, "iii"),
])
def test_u3_rule_at_its_boundaries(od, q, floored, outcome):
    assert od.u3_outcome(q, floored) == outcome


def _recs(rows):
    return pd.DataFrame([{"target": t, "seed": s, "world": w, "variant": v, "t3": x, "sev_floor": f}
                         for t, s, w, v, x, f in rows])


def _committed(od, b=(60.0, 62.0), v95=(48.0, 50.0)):
    def table(xs):
        return pd.DataFrame([{"target": t, "seed": s, "t3": x, "sev_floor": False}
                             for (t, s), x in zip([("rep1", 0), ("rep1", 1)], xs, strict=True)])

    return {od.V_B: {"rows": table(b), "mean": float(np.mean(b))},
            "transfer:mab21l2_vs95": {"rows": table(v95), "mean": float(np.mean(v95))}}


def test_summarise_u3_q_hat_and_paired_gain(od):
    committed = _committed(od, b=(60.0, 62.0), v95=(48.0, 50.0))  # means 61 and 49
    recs = _recs([("rep1", 0, "as_released", od.V_PERP, 58.0, "True"),
                  ("rep1", 1, "as_released", od.V_PERP, 60.0, False),
                  ("rep1", 0, "as_released", od.V_B, 61.0, False), ("rep1", 1, "as_released", od.V_B, 62.5, False)])
    s = od.summarise_u3(recs, committed, expect_rows=2)
    assert s["mean_t3"] == 59.0 and s["rows_at_severity_floor"] == 1
    assert s["q_hat"] == pytest.approx((59.0 - 49.0) / (61.0 - 49.0))
    assert s["paired_gain_over_B_rerun"] == pytest.approx(np.mean([58.0 - 61.0, 60.0 - 62.5]))
    assert s["outcome"] == "i" and s["per_replicate"]["rep1"] == 59.0 and "no-answer and zero-answer" in s["reading"]
    with pytest.raises(SystemExit, match="expected 20"):
        od.summarise_u3(recs, committed, expect_rows=20)


def test_summarise_u3_after_a_noop_removal_reads_inconclusive(od):
    committed = _committed(od)
    recs = _recs([("rep1", 0, "as_released", od.V_B, 61.0, False), ("rep1", 1, "as_released", od.V_B, 62.0, False)])
    s = od.summarise_u3(recs, committed, expect_rows=0)
    assert s["outcome"] == "no-op" and s["q_hat"] is None and s["n_rows"] == 0
    assert "inconclusive" in s["reading"] and "B-perp was not scored" in s["reading"]


def test_replay_check_matches_within_tolerance_and_reads_floor_flags_as_strings(od):
    committed = _committed(od, b=(60.0, 62.0))
    good = _recs([("rep1", 0, "as_released", od.V_B, 60.005, False), ("rep1", 1, "as_released", od.V_B, 62.01, False)])
    assert od.replay_check(good, committed)["matched"]
    far = _recs([("rep1", 0, "as_released", od.V_B, 60.5, False), ("rep1", 1, "as_released", od.V_B, 62.0, False)])
    assert not od.replay_check(far, committed)["matched"]
    flag = _recs([("rep1", 0, "as_released", od.V_B, 60.0, True), ("rep1", 1, "as_released", od.V_B, 62.0, False)])
    out = od.replay_check(flag, committed)
    assert not out["floor_flags_equal"] and not out["matched"]
    # flags as read back from a CSV ("True", "False") compare with booleans
    committed[od.V_B]["rows"]["sev_floor"] = ["False", "False"]
    assert od.replay_check(good, committed)["matched"]
    assert not od.replay_check(_recs([("rep2", 0, "as_released", od.V_B, 60.0, False)]), committed)["matched"]


def test_committed_standins_enforce_the_row_numbers(od, tmp_path):
    rows = [{"target": t, "seed": s, "variant": v, "t3": loc, "sev_floor": False}
            for v, loc in (("transfer:mab21l2_all", 60.7745), ("transfer:mab21l2_vs95", 49.1371))
            for t in ("rep1", "rep2") for s in range(10)]
    p = tmp_path / "r.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    assert od.committed_standins(p)["transfer:mab21l2_all"]["mean"] == pytest.approx(60.7745)
    pd.DataFrame(rows[:-1]).to_csv(p, index=False)
    with pytest.raises(SystemExit, match="expected 20"):
        od.committed_standins(p)
    rows[0]["t3"] += 2.0
    pd.DataFrame(rows).to_csv(p, index=False)
    with pytest.raises(SystemExit, match="differs from the 60.7745"):
        od.committed_standins(p)
    assert len(od.committed_standins(p, check=False)["transfer:mab21l2_all"]["rows"]) == 20  # a code check: not held


def test_refuse_existing_needs_a_registered_row_other_than_its_own(od, tmp_path):
    existing = tmp_path / "out.json"
    existing.write_text("{}")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        od.refuse_existing([existing], None)
    od.refuse_existing([tmp_path / "absent.json"], None)  # nothing to overwrite
    rows = od.tries_rows()
    assert "N12" in rows and "E11" in rows and "K2b" in rows
    od.refuse_existing([existing], "N12")  # a row of the ledger other than U3 and U3d
    for bad in ("U3", "U3d", "not-a-row", "Z99"):
        with pytest.raises(SystemExit, match="--rerun-row"):
            od.refuse_existing([existing], bad)


# ---------------------------------------------------------------------------------------------------------------------
# The depth matching


def test_matched_counts_apply_each_labels_ratio_through_rescale_depth(od):
    t3d = od.t3d
    rng = np.random.default_rng(3)
    C = rng.poisson(30, size=(6000, 8)).astype(np.int32)
    props = np.full(C.shape, 1 / 8, dtype=np.float32)
    lab = np.tile(np.arange(-1, 3), 1500)  # labels -1, 0, 1 and 2 in turn
    ratio = np.array([0.5, 1.0, 1.6])
    out = od.matched_counts(C, lab, ratio, props, np.random.default_rng(4))
    r = np.array([1.0, 0.5, 1.0, 1.6])[lab + 1]  # a label's ratio per cell; label -1 keeps its depth
    want = t3d.rescale_depth(C, r, props, np.random.default_rng(4))
    assert np.array_equal(out, want)
    for k, kept in ((-1, 1.0), (0, 0.5), (1, 1.0), (2, 1.6)):
        got = out[lab == k].sum() / C[lab == k].sum()
        assert abs(got - kept) < 0.03, (k, got)


def test_reference_counts_reject_a_reference_that_does_not_reproduce_its_expression(od, tmp_path):
    rng = np.random.default_rng(5)
    C = rng.poisson(20, size=(300, 12)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    ct = np.array(["V-CM", "NCC", "Foregut"])[rng.integers(0, 3, 300)]
    obs = pd.DataFrame({"celltype": ct, "cm_celltype": ["Unknown"] * 300}, index=[f"c{i}" for i in range(300)])
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=[f"G{j}" for j in range(12)]))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 9, 300), rng.uniform(0, 9, 300), rng.integers(0, 5, 300) * 1.0]
    p = tmp_path / "wt.h5ad"
    a.write_h5ad(p)
    V = SimpleNamespace(labels=["V-CM", "NCC"])  # Foregut is outside this vocabulary: label -1
    ctx = SimpleNamespace(genes=[f"G{j}" for j in range(12)], ref_X=np.asarray(a.X.toarray(), np.float32), V=V,
                          rules=od.bs.load_rules())
    C2, lab, props, gap = od.reference_counts(ctx, p, "combined")
    assert np.array_equal(C2, C) and gap < 1e-3
    assert set(lab[ct == "Foregut"]) == {-1} and set(lab[ct == "NCC"]) == {1}
    assert props.shape == C.shape and np.allclose(props.sum(1), 1.0, atol=1e-5)
    ctx.ref_X = ctx.ref_X + 0.05  # the released expression does not come from these counts: stop and look
    with pytest.raises(SystemExit, match="do not reproduce"):
        od.reference_counts(ctx, p, "combined")
    ctx.ref_X = np.asarray(a.X.toarray(), np.float32)[:-1]
    with pytest.raises(SystemExit, match="do not line up"):
        od.reference_counts(ctx, p, "combined")


# ---------------------------------------------------------------------------------------------------------------------
# N11 with the measured cells


def _n11_rows(rng):
    v = {"xrep:F+depth": 74.0, "fac:depth": 55.0, "blind:G2_bcat_central_d50": 50.5, "transfer:mab21l2_all": 60.8,
         "transfer:mab21l2_vs95": 49.1}
    keys = [(t, s) for t in ("rep1", "rep2") for s in range(10)]
    return {name: {"key": keys, "t3": rng.normal(m, 2.0, 20), "floor": rng.random(20) < 0.3} for name, m in v.items()}


def _others(rng):
    def o():
        return {"a": rng.normal(60, 3, 7), "b": rng.normal(55, 3, 7), "w": np.ones(7) / 7}

    return {"T1": o(), "heart": o(), "embryo": o()}


@pytest.mark.parametrize("pf", [0.0, 1.0])
def test_simulate_measured_equals_n11_when_the_cells_equal_n11s(od, monkeypatch, pf):
    """Where the measured cells equal N11's (the judgement cell made a point: SD 0, floor probability 0 or 1; the
    barred file as committed), the rerun is N11's simulate number for number: same worlds, same random numbers."""
    S = od.S
    monkeypatch.setattr(od.P, "SD_J", 0.0)
    monkeypatch.setattr(od.E, "N_DRAWS", FAST_DRAWS)
    rng = np.random.default_rng(6)
    rows, others = _n11_rows(rng), _others(rng)
    m = 51.0
    measured = {od.V_S1_DM: {"t3": np.full(20, m), "floor": np.full(20, bool(pf))},
                od.V_B_DM: rows["transfer:mab21l2_all"]}
    for pu, pa, q in ((0.7, 0.65, 1.0), (0.7, 0.5, 0.75), (0.7, 0.65, 0.25)):
        a = S.simulate(rows, others, pu, pa, (m, pf), q, np.random.default_rng(od.P.SEED))
        b = od.simulate_measured(rows, measured, others, pu, pa, q, np.random.default_rng(od.P.SEED))
        for reading in ("half", "zero"):
            assert b[reading]["dE"] == pytest.approx(a[reading]["dE"], abs=1e-12), reading
            assert b[reading]["dP"]["lambda1.0_q90"] == pytest.approx(a[reading]["dP_lambda1"][90], abs=1e-12)
            assert b[reading]["dP"]["lambda1.0_q95"] == pytest.approx(a[reading]["dP_lambda1"][95], abs=1e-12)
            assert b[reading]["wc_gain"] == pytest.approx(a[reading]["wc_gain"], abs=1e-12)
            assert b[reading]["p_wc"] == pytest.approx(a[reading]["p_wc"], abs=1e-12)
        assert b["mixture"]["dE"] == pytest.approx(a["mixture"]["dE"], abs=1e-12)
        assert set(b["half"]["dP"]) == {"lambda0.5_q90", "lambda0.5_q95", "lambda1.0_q90", "lambda1.0_q95"}


def test_simulate_measured_changes_only_the_not_shallow_world(od, monkeypatch):
    """The measured cells enter only in U1 & anchor & not U2; a lower measured barred score there lowers the swap's
    gain in that world, and the gain outside it does not depend on the measured cells."""
    monkeypatch.setattr(od.E, "N_DRAWS", FAST_DRAWS)
    rng = np.random.default_rng(7)
    rows, others = _n11_rows(rng), _others(rng)
    keep = {od.V_S1_DM: {"t3": np.full(20, 52.0), "floor": np.zeros(20, bool)},
            od.V_B_DM: {"t3": np.full(20, 60.0), "floor": np.zeros(20, bool)}}
    low = {od.V_S1_DM: keep[od.V_S1_DM], od.V_B_DM: {"t3": np.full(20, 40.0), "floor": np.zeros(20, bool)}}
    a = od.simulate_measured(rows, keep, others, 0.7, 0.65, 1.0, np.random.default_rng(1))
    b = od.simulate_measured(rows, low, others, 0.7, 0.65, 1.0, np.random.default_rng(1))
    assert b["half"]["wc_gain"] < a["half"]["wc_gain"] and b["half"]["dE"] < a["half"]["dE"]
    assert a["half"]["p_wc"] == b["half"]["p_wc"]
    # the difference is exactly the world's gain difference, times that world's probability
    assert (a["half"]["dE"] - b["half"]["dE"]) == pytest.approx(
        (a["half"]["wc_gain"] - b["half"]["wc_gain"]) * a["half"]["p_wc"], abs=1e-9)


def test_measured_rows_follow_n11s_row_order_and_refuse_a_mismatch(od):
    keys = [("rep1", 1), ("rep1", 0), ("rep2", 0)]
    rec = _recs([("rep2", 0, "depth_matched", od.V_B_DM, 3.0, "True"),
                 ("rep1", 1, "depth_matched", od.V_B_DM, 1.0, False),
                 ("rep1", 0, "depth_matched", od.V_B_DM, 2.0, "False"),
                 ("rep1", 1, "depth_matched", od.V_S1_DM, 10.0, False),
                 ("rep1", 0, "depth_matched", od.V_S1_DM, 20.0, False),
                 ("rep2", 0, "depth_matched", od.V_S1_DM, 30.0, True)])
    out = od.measured_rows(rec, keys)
    assert out[od.V_B_DM]["t3"].tolist() == [1.0, 2.0, 3.0]
    assert out[od.V_S1_DM]["t3"].tolist() == [10.0, 20.0, 30.0]
    for v in (od.V_B_DM, od.V_S1_DM):
        assert out[v]["floor"].tolist() == [False, False, True]
    with pytest.raises(SystemExit, match="do not pair"):
        od.measured_rows(rec, [("rep1", 1), ("rep1", 0)])


def _cell(mix_de, dp_ok=True):
    keys = [f"lambda{lam}_q{q}" for lam in (0.5, 1.0) for q in (90, 95)]
    return {"mixture": {"dE": mix_de, "dP": dict.fromkeys(keys, 0.1 if dp_ok else -0.01)}}


def test_no_answer_rule_needs_rule_i_and_n10s_bar_at_both_bands(od):
    good = [_cell(0.6), _cell(0.5)]
    assert od.no_answer_rule("i", good)["change_no_answer_default_to_swap"]
    assert not od.no_answer_rule("iii", good)["change_no_answer_default_to_swap"]
    assert not od.no_answer_rule("no-op", [])["change_no_answer_default_to_swap"]  # no q-hat, no cells at q-hat
    assert not od.no_answer_rule("i", [_cell(0.6), _cell(0.49)])["n10_bar_met_at_both_bands"]
    assert not od.no_answer_rule("i", [_cell(0.9), _cell(0.8, dp_ok=False)])["change_no_answer_default_to_swap"]
    assert od.n10_bar_met(_cell(0.5)) and not od.n10_bar_met(_cell(0.49))


def test_n11_command_rejects_a_missing_run_and_writes_stamped_outputs(od, tmp_path, monkeypatch, own_repo):
    """`n11` on synthetic records of the two measured variants; the rest of N11 (the other boards, the committed rows)
    are the committed results. q-hat comes from the decision file; a no-op run has no q-hat."""
    with pytest.raises(SystemExit, match="run `orth_diag.py score` first"):
        od.main(["n11", "--out", str(tmp_path / "none")])
    monkeypatch.setattr(od.E, "N_DRAWS", FAST_DRAWS)
    rng = np.random.default_rng(8)
    rows = [{"target": t, "seed": s, "world": "depth_matched", "variant": v, "t3": float(rng.normal(m, 2)),
             "sev_floor": bool(rng.random() < 0.2)}
            for v, m in ((od.V_B_DM, 57.0), (od.V_S1_DM, 58.0)) for t in ("rep1", "rep2") for s in range(10)]
    for name, q in (("full", 0.82), ("noop", None)):
        out = tmp_path / name
        out.mkdir()
        pd.DataFrame(rows).to_csv(out / od.RECORDS_CSV, index=False)
        (out / od.DECISION_JSON).write_text(json.dumps({"u3": {"q_hat": q, "outcome": "i" if q else "no-op"}}))
        od.main(["n11", "--out", str(out)])
        res = json.loads((out / od.N11_JSON).read_text())
        assert res["q_hat"] == q and res["stamp"]["commit"]
        grid = {c["q_barred"] for c in res["cells"]}
        assert grid == ({0.82, 1.0, 0.75, 0.5, 0.25} if q else {1.0, 0.75, 0.5, 0.25})
        assert len(res["cells"]) == 2 * len(grid)  # both N11 bands
        assert read_stamp(stamp_path(out / od.N11_JSON))["sha256"]
        with pytest.raises(SystemExit, match="refusing to overwrite"):
            od.main(["n11", "--out", str(out)])
        if q is None:
            assert not res["no_answer_rule"]["change_no_answer_default_to_swap"]
            assert "no q-hat" in (out / od.N11_TXT).read_text()


# ---------------------------------------------------------------------------------------------------------------------
# The commands, on synthetic files


@pytest.mark.filterwarnings("ignore:.*obs names of the model's training cells")
def test_axis_gate_is_t3_checks_formula_on_the_hashed_factor(od, world, run_args, own_repo, capsys):
    od.main(run_args.axis)
    res = json.loads((run_args.out / od.AXIS_JSON).read_text())
    g = res["gate"]
    assert g["self_check_passed"] and g["barred_reproduced"] and g["proceed"] and not g["noop"]
    assert g["corr_fM_d"] == pytest.approx(world.corr, abs=1e-6)  # d5d6's formula, against the hashed factor
    assert g["barred"]["labels_in_mask"] == world.labels_used and g["barred"]["max_abs_diff_log_phi"] < 1e-6
    assert res["script_sha256"] == od.t3d.sha256(Path(od.__file__))
    assert (run_args.out / od.AXIS_TXT).read_text() in capsys.readouterr().out
    # B-perp is exactly the hashed factor minus its stage part, by the OLS slope, rescaled by the mask's denominator
    perp = pd.read_csv(run_args.fdir / od.PERP_CSV, float_precision="round_trip")
    assert res["factor_perp"]["sha256"] == od.t3d.sha256(run_args.fdir / od.PERP_CSV)
    ctx = world.ctx
    assert list(perp["gene"]) == list(ctx.genes)
    v = np.log(world.phi_all)
    b = float((v - v.mean()) @ world.d1 / (world.d1 @ world.d1))
    assert g["removal"]["b"] == pytest.approx(b, abs=1e-6) and abs(b) > 1e-3
    mab = load_stage(world.data / "T3" / "E9.5_mab21l2_ko.h5ad", genes=ctx.genes, drop_duplicates=True)
    calls, _ = ctx.harm.calls(mab, seed=world.rules["harmonisation"]["seed"])
    lab = od.t3d.calls_to_index(ctx.V, calls)
    mu = od.t3d.props_type_means(od.t3d.dense_counts(mab), lab, ctx.V.L)
    cnt = np.bincount(lab[lab >= 0], minlength=ctx.V.L)
    both = (cnt >= 100) & (ctx.wt.count >= 100)
    fc = world.rules["components"]["factor"]
    _, den = od.t3d.gene_factor(mu, ctx.wt.mu, ctx.V.p, both, fc["eps"], tuple(fc["clip"]))
    raw = np.exp(v - b * world.d1)
    want = np.clip(raw / (den @ raw / den.sum()), *tuple(fc["clip"]))
    np.testing.assert_allclose(perp["factor"].to_numpy(np.float64), want, rtol=0, atol=1e-5)
    assert read_stamp(stamp_path(run_args.fdir / od.PERP_CSV))["sha256"]
    assert read_stamp(stamp_path(run_args.out / od.AXIS_JSON))["sha256"]
    assert res["stamp"]["tree_clean"] is True and res["stamp"]["commit"] not in ("", "unknown")
    assert set(res["inputs_sha256"]) >= {str(world.decision), str(world.d5d6), str(world.meta)}
    assert g["reported_only"]["secondary_axis"]["scored"] is False  # the secondary axis is reported, not scored
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        od.main(run_args.axis)


@pytest.mark.filterwarnings("ignore:.*obs names of the model's training cells")
def test_axis_exits_nonzero_when_the_self_check_fails(od, world, run_args, own_repo, tmp_path):
    bad = tmp_path / "d5d6_bad.json"
    d = json.loads(world.d5d6.read_text())
    d["D5"]["primary"]["corr_fM_d"] = world.corr + 0.05
    bad.write_text(json.dumps(d))
    args = list(run_args.axis)
    args[args.index("--d5d6") + 1] = str(bad)
    args[args.index("--out") + 1] = str(tmp_path / "out_bad")
    with pytest.raises(SystemExit, match="the gate failed"):
        od.main(args)
    res = json.loads((tmp_path / "out_bad" / od.AXIS_JSON).read_text())
    assert not res["gate"]["self_check_passed"] and not res["gate"]["proceed"]
    assert "STOP" in (tmp_path / "out_bad" / od.AXIS_TXT).read_text()
    with pytest.raises(SystemExit, match="did not pass"):  # scoring refuses a failed gate
        od.main(["score", "--out", str(tmp_path / "out_bad"), "--data", str(world.data), "--run-dir", str(world.run),
                 "--decision", str(world.decision), "--records", str(world.records), "--workers", "1", "--threads", "1",
                 "--seeds", "0", "--code-check"])


@pytest.mark.filterwarnings("ignore:.*obs names of the model's training cells")
def test_score_flow_rows_and_outputs(od, world, run_args, own_repo, capsys):
    od.main(run_args.axis)
    gate = json.loads((run_args.out / od.AXIS_JSON).read_text())
    capsys.readouterr()
    od.main(run_args.score)
    printed = capsys.readouterr().out
    recs = pd.read_csv(run_args.out / od.RECORDS_CSV)
    dec = json.loads((run_args.out / od.DECISION_JSON).read_text())
    assert (run_args.out / od.DECISION_TXT).read_text() in printed
    assert dec["seeds"] == [0, 1] and dec["protocol_as_registered"] is False
    want = {(od.V_B, "as_released"), (od.V_PERP, "as_released"), (od.V_B_DM, "depth_matched"),
            (od.V_S1_DM, "depth_matched")}
    assert set(zip(recs["variant"], recs["world"], strict=True)) == want
    assert len(recs) == 2 * 2 * 4 and set(recs["target"]) == {"rep1", "rep2"} and recs["t3"].notna().all()
    assert dec["tries"]["U3"] == 1 and dec["tries"]["U3d"] == 2
    # B-perp's rows are the rule's: the committed fake means give q-hat, and the rerun of B cannot match the fake rows
    perp = recs[recs["variant"] == od.V_PERP]
    ref = pd.read_csv(world.records)
    m_b, m_95 = (ref[ref["variant"] == v]["t3"].mean() for v in ("transfer:mab21l2_all", "transfer:mab21l2_vs95"))
    assert dec["u3"]["q_hat"] == pytest.approx((perp["t3"].mean() - m_95) / (m_b - m_95))
    assert not dec["replay_of_B"]["matched"] and dec["tries"]["rerun_of_B"] == 1
    assert dec["u3d"][od.V_S1_DM]["n_rows"] == 4
    # the depth step really ran: depth kept below 1 where the replicates are shallower than WT (65-75 against 130)
    assert 0.2 < dec["u3d"][od.V_B_DM]["mean_dm_pool_depth_kept"] < 1.0
    assert 0.2 < dec["u3d"][od.V_B_DM]["mean_dm_ref_depth_kept"] < 1.0
    # the stamp, the inputs, the refusal to overwrite
    assert read_stamp(stamp_path(run_args.out / od.RECORDS_CSV))["sha256"]
    assert dec["stamp"]["tree_clean"] is True
    assert set(dec["inputs_sha256"]) >= {str(run_args.out / od.AXIS_JSON), str(world.records)}
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        od.main(run_args.score)
    # a gate from another script version is refused before anything is loaded (the rerun row lets it past the refusal)
    (run_args.out / od.AXIS_JSON).write_text(json.dumps({**gate, "script_sha256": "0" * 64}))
    with pytest.raises(SystemExit, match="another version of this script"):
        od.main([*run_args.score, "--rerun-row", "N12"])


@pytest.mark.filterwarnings("ignore:.*obs names of the model's training cells")
def test_rows_pair_with_the_t3_scripts_own_scored_variant(od, world, ctx_targets, monkeypatch):
    """The B rerun's row equals t3_decompose.score_seed's `transfer:mab21l2_all` row on the same context and streams, so
    B-perp's rows pair with the committed B rows; B-perp is built on B's cells; and U3d's depth-matched row is the
    scorer's own on the matched reference with its own anchors."""
    t3d, ctx = od.t3d, ctx_targets.ctx
    C, lab, props, _ = ctx_targets.ref
    ctx.transfer = {"mab21l2_all": world.phi_all}
    phi_perp = np.clip(np.exp(np.log(world.phi_all) - 0.2 * world.d1), 0.25, 4.0)
    st = {"ctx": ctx, "phi_b": world.phi_all, "phi_perp": phi_perp, "ref_counts": C, "ref_lab": lab, "ref_props": props}
    seed, ti = 3, 0
    # the T3 script's row of B for the same (target, seed), through score_seed with only that variant left in
    real = t3d.variant_specs
    monkeypatch.setattr(t3d, "variant_specs", lambda c, i, e: [s for s in real(c, i, e) if s[0] == od.V_B])
    theirs = {r["variant"]: r for r in t3d.score_seed(ctx, ti, seed)}[od.V_B]
    monkeypatch.setattr(t3d, "variant_specs", real)
    captured, rescales = [], []
    real_apply, real_rescale = t3d.apply_gene_factor, t3d.rescale_depth

    def spy(C_, phi, P_, rng):
        captured.append((C_.copy(), np.array(phi)))
        return real_apply(C_, phi, P_, rng)

    def spy_rescale(C_, ratio, P_, rng):
        rescales.append((C_.shape[0], np.array(ratio)))
        return real_rescale(C_, ratio, P_, rng)

    monkeypatch.setattr(t3d, "apply_gene_factor", spy)
    monkeypatch.setattr(t3d, "rescale_depth", spy_rescale)
    ours = {(r["variant"], r["world"]): r for r in od.score_row(st, ti, seed)}
    monkeypatch.undo()
    b = ours[(od.V_B, "as_released")]
    for c in ("t3", "de_score", "de_direction", "severity_slope", "mmd_u"):
        assert b[c] == pytest.approx(theirs[c], abs=1e-12, nan_ok=True), c
    assert b["sev_floor"] == theirs["sev_floor"]
    # B and B-perp factor the same drawn cells (the stream's first draws are the cells), through different factors
    (Cb, phib), (Cp, phip) = captured[0], captured[1]
    assert np.array_equal(Cb, Cp) and np.array_equal(phib, world.phi_all) and np.array_equal(phip, phi_perp)
    # U3d: the matched pool gives B's cells at the replicate's depth; s1's stand-in draws from it too, then its factor
    C_dm, phi_dm = captured[2]
    assert np.array_equal(phi_dm, world.phi_all)
    assert C_dm.shape == Cb.shape and (C_dm.sum(1) <= Cb.sum(1)).all() and C_dm.sum() < 0.8 * Cb.sum()
    C_s1, phi_s1 = captured[3]
    S1 = ctx.targets[1 - ti].full
    assert np.array_equal(phi_s1, S1.phi_nm) and C_s1.sum() < 0.8 * ctx.wt.pool.C.sum(1).mean() * C_s1.shape[0]
    # the three depth steps of the matched world (reference copy, pool, s1's own), after the estimate's
    pr = world.rules["protocol"]
    T = ctx.targets[ti]
    rows_r = t3d.working_copy_rows(len(ctx.ref_X), pr["working_copy_fraction"], seed + 1)
    rows_t = t3d.working_copy_rows(len(T.X_rel), pr["working_copy_fraction"], seed)
    truth, other = split_half(Stage(X=T.X_rel[rows_t]), seed=pr["split_seed"])
    ref = Stage(X=ctx.ref_X[rows_r])
    n_pred = pr["prediction_cells_multiplier"] * max(truth.n, ref.n)
    assert [n for n, _ in rescales[1:]] == [len(rows_r), ctx.wt.pool.C.shape[0], n_pred]
    # the last depth step is s1's own: the other replicate's per-label ratios, in the build
    assert set(np.unique(rescales[-1][1])) <= set(np.unique(S1.ratio)) and not set(np.unique(rescales[-1][1])) <= {1.0}
    # the matched world's rows are the scorer's on the matched reference, with anchors recomputed on it
    pool_rows = np.flatnonzero(~np.isin(T.hash_pool, T.hash_rel[rows_t]))
    est = t3d.estimate(ctx, T, pool_rows, t3d.stable_rng(seed, f"estimate:{T.name}"))
    C_r = od.matched_counts(C[rows_r], lab[rows_r], est.ratio, props[rows_r], t3d.stable_rng(seed, f"dm_ref:{T.name}"))
    C_p = od.matched_counts(ctx.wt.pool.C, ctx.wt.pool.lab, est.ratio, ctx.wt.pool.P,
                            t3d.stable_rng(seed, f"dm_pool:{T.name}"))
    ref_m = Stage(X=t3d.lognorm(C_r))
    src = t3d.Pool(C_p, ctx.wt.pool.lab, ctx.wt.pool.P, ctx.V.L)
    X_bm = t3d.build(ctx, n_pred, t3d.stable_rng(seed, od.B_STREAM), f=ctx.V.p, src=src, phi=world.phi_all)
    X_s1 = t3d.build(ctx, n_pred, t3d.stable_rng(seed, od.S1_STREAM), f=ctx.V.p, src=src, phi=S1.phi_nm, ratio=S1.ratio)
    floor = panel_metrics("T3", ref_m, truth, ref_m, seed=seed)
    ceil = panel_metrics("T3", other, truth, ref_m, seed=seed)
    anchors = {m: {"floor": floor.get(m), "ceiling": ceil.get(m)} for m in RANKED_METRICS["T3"]
               if floor.get(m) is not None and ceil.get(m) is not None}
    for variant, X in ((od.V_B_DM, X_bm), (od.V_S1_DM, X_s1)):
        skills = skills_from_anchors("T3", panel_metrics("T3", Stage(X=X), truth, ref_m, seed=seed), anchors)
        assert ours[(variant, "depth_matched")]["t3"] == pytest.approx(aggregate("T3", skills)[0], abs=1e-12), variant
    assert ours[(od.V_B_DM, "depth_matched")]["dm_pool_depth_kept"] == pytest.approx(C_p.sum() / ctx.wt.pool.C.sum())


@pytest.mark.filterwarnings("ignore:.*obs names of the model's training cells")
def test_noop_removal_scores_u3d_only(od, world, ctx_targets):
    """After a no-op removal `score_row` skips B-perp and keeps U3d's two rows and the rerun of B."""
    C, lab, props, _ = ctx_targets.ref
    st = {"ctx": ctx_targets.ctx, "phi_b": world.phi_all, "phi_perp": None, "ref_counts": C, "ref_lab": lab,
          "ref_props": props}
    recs = od.score_row(st, 1, 0)
    assert {r["variant"] for r in recs} == {od.V_B, od.V_B_DM, od.V_S1_DM}


@pytest.mark.filterwarnings("ignore:.*obs names of the model's training cells")
def test_noop_flow_through_the_commands(od, world, run_args, own_repo, monkeypatch, capsys):
    """A no-op removal: `axis` reports it and exits normally; `score` then runs U3d only, with no try for U3."""
    monkeypatch.setattr(od, "NOOP_SD", 100.0)  # any removal is a no-op against this threshold
    od.main(run_args.axis)
    assert json.loads((run_args.out / od.AXIS_JSON).read_text())["gate"]["noop"]
    assert "U3 is inconclusive" in (run_args.out / od.AXIS_TXT).read_text()
    od.main(run_args.score)
    capsys.readouterr()
    recs = pd.read_csv(run_args.out / od.RECORDS_CSV)
    dec = json.loads((run_args.out / od.DECISION_JSON).read_text())
    assert od.V_PERP not in set(recs["variant"]) and {od.V_B, od.V_B_DM, od.V_S1_DM} <= set(recs["variant"])
    assert dec["u3"]["outcome"] == "no-op" and dec["u3"]["q_hat"] is None and dec["tries"]["U3"] == 0
    assert "B-perp was not scored" in (run_args.out / od.DECISION_TXT).read_text()
    assert str(run_args.fdir / od.PERP_CSV) not in dec["inputs_sha256"]  # the factor file is not an input of this run
