"""Synthetic checks of the T1 direction readout (no competition or external data).

Run: uv run pytest experiments/t1_direction -q
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import yaml
from common.core_metrics import de_direction, de_score
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import dirlib as dl  # noqa: E402
import readout  # noqa: E402
import t1_stats  # noqa: E402

CFG = yaml.safe_load((HERE / "config.yaml").read_text())


# ------------------------------------------------------------------ metrics agree with veckit's matrix versions
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_vector_metrics_match_veckit(seed):
    rng = np.random.default_rng(seed)
    G, n = 300, 400
    mu = rng.uniform(1.0, 3.0, G)
    d_true = np.where(rng.random(G) < 0.15, rng.choice([-1.0, 1.0], G), 0.0)  # clear DE genes, BH step inert
    d_pred = 0.6 * d_true + rng.normal(0, 0.4, G) * (rng.random(G) < 0.5)
    ref_X = mu + rng.normal(0, 0.3, (n, G))
    true_X = mu + d_true + rng.normal(0, 0.3, (n, G))
    pred_X = mu + d_pred + rng.normal(0, 0.3, (n, G))
    pb = lambda X: X.mean(0)  # noqa: E731
    dp, dt, ref = pb(pred_X) - pb(ref_X), pb(true_X) - pb(ref_X), pb(ref_X)
    assert dl.de_score_vec(dp, dt, ref)[0] == pytest.approx(de_score(pred_X, true_X, ref_X)["score"], abs=1e-12)
    assert dl.de_direction_vec(dp, dt, ref) == pytest.approx(de_direction(pred_X, true_X, ref_X), abs=1e-12)


def test_metric_guards():
    rng = np.random.default_rng(0)
    dt, ref = rng.normal(0, 0.5, 200), rng.uniform(0, 2, 200)
    assert dl.de_score_vec(np.zeros(200), dt, ref)[0] == 0.0
    assert dl.de_direction_vec(np.zeros(200), dt, ref) == 0.0
    assert np.isnan(dl.de_score_vec(dt, np.zeros(200), ref)[0])


# ------------------------------------------------------------------ decision rule
RULE = CFG["decision"]


@pytest.mark.parametrize("p,b,cm,verdict,failed", [
    ((0.30, 0.30), (0.20, 0.20), 0.80, "build", []),
    ((0.25, 0.30), (0.20, 0.20), 0.80, "build", []),                     # exactly 0.05 passes
    ((0.24, 0.30), (0.20, 0.20), 0.80, "do not build", ["delta_de_score"]),
    ((0.30, 0.24), (0.20, 0.20), 0.80, "do not build", ["delta_de_direction"]),
    ((0.30, 0.30), (0.20, 0.20), 0.69, "do not build", ["cm_sign_agreement"]),
    ((float("nan"), 0.30), (0.20, 0.20), 0.80, "do not build", ["delta_de_score"]),
    ((0.10, 0.10), (0.20, 0.20), 0.90, "do not build", ["delta_de_score", "delta_de_direction"]),
])
def test_decision_rule(p, b, cm, verdict, failed):
    readings = {"primary": {"de_score": p[0], "de_direction": p[1]},
                "baseline": {"de_score": b[0], "de_direction": b[1]}, "cm_sign_agreement": cm}
    out = dl.decide(RULE, readings)
    assert out["verdict"] == verdict and out["failed"] == failed


def test_decision_fails_closed_without_readings():
    out = dl.decide(RULE, {})
    assert out["verdict"] == "do not build" and len(out["failed"]) == 3


def test_preregistered_values():
    assert RULE["primary_predictor"] == "feng_E14.5" and RULE["baseline"] == "mom_xfit"
    assert RULE["min_delta_de_score"] == RULE["min_delta_de_direction"] == 0.05
    assert RULE["min_cm_sign_agreement"] == 0.70 and RULE["population_design"] == "matched_coverage"


# ------------------------------------------------------------------ harmonisation
@pytest.mark.parametrize("label,group", [
    ("OFT/RV-CM", "CM"), ("IFT-CM", "CM"), ("V-CM", "ventricularCM"), ("LV-CM", "ventricularCM"),
    ("Trabecular CM", "ventricularCM"), ("Compact-CM", "ventricularCM"), ("Endocardium-1", "Endocardium"),
    ("BEC", "Endothelium"), ("JCF", "SHF"), ("aSHF", "SHF"), ("Proepicardium", "Epicardium"), ("EPDC", "Epicardium"),
    ("ST", "ST"), ("Hepatocyte", "Hepatocyte"), ("Epithelium", "Foregut"), ("NCC-derived", "NCC"), ("aPHM", "PHM"),
    ("Surface Ectoderm", "SurfaceEcto"), ("Cushion mesenchyme", "Mesenchyme"), ("EXEM", "other:EXEM"),
])
def test_harmonise(label, group):
    assert dl.harmonise([label], CFG["harmonise"])[label] == group


def test_override_and_fallback():
    m = dl.harmonise(["EXEM", "Endothelium", "RV-CM", "OFT/RV-CM"], CFG["harmonise"], {"EXEM": "Blood"})
    assert m["EXEM"] == "Blood"
    assert dl.group_labels(m, "CM", CFG["cm_includes"]) == ["OFT/RV-CM", "RV-CM"]
    assert dl.antecedent_labels(m, "Endocardium", CFG["cm_includes"], CFG["fallback"]) == ["Endothelium"]
    assert dl.antecedent_labels(m, "Epicardium", CFG["cm_includes"], CFG["fallback"]) == []


# ------------------------------------------------------------------ synthetic statistics
GENES = np.array([f"g{i}" for i in range(298)] + ["Xist", "Cwc22"], dtype=object)


def t1_like(means: dict[str, np.ndarray], n: int = 60, det: dict[str, np.ndarray] | None = None) -> dict:
    keys, N, SL, SC, SZ = [], [], [], [], []
    for lab, m in means.items():
        for h in "AB":
            keys.append(f"{lab}|1|{h}")
            N.append(n / 2)
            SL.append(n / 2 * m)
            SC.append(n / 2 * np.expm1(m))
            SZ.append(n / 2 * (det[lab] if det else (m > 0).astype(float)))
    return {"genes": GENES, "keys": np.array(keys, dtype=object), "n": np.array(N), "S_log": np.vstack(SL),
            "S_cp": np.vstack(SC), "S_nz": np.vstack(SZ), "depth_sum": np.array(N) * 1e4}


def feng_like(cps: dict[tuple[str, str], np.ndarray], n: int = 40) -> dict:
    keys, N, SL, SC, SZ = [], [], [], [], []
    for (stage, typ), cp in cps.items():
        for h in "AB":
            keys.append(f"{stage}|LV|{typ}|{h}")
            N.append(n / 2)
            SC.append(n / 2 * cp)
            SL.append(n / 2 * np.log1p(cp))
            SZ.append(np.full(len(cp), n / 2))
    return {"genes": GENES.copy(), "keys": np.array(keys, dtype=object), "n": np.array(N), "S_log": np.vstack(SL),
            "S_cp": np.vstack(SC), "S_nz": np.vstack(SZ), "depth_sum": np.array(N) * 1e4}


def synthetic_readout(external_right: bool, seed: int = 0):
    rng = np.random.default_rng(seed)
    G = len(GENES)
    fam = {"CM": ["OFT/RV-CM", "V-CM"], "Endo": ["Endocardium"], "BEC": ["BEC"], "Blood": ["Blood"], "SHF": ["pSHF"]}
    true = {f: np.clip(rng.normal(0, 0.5, G), -1, 1) for f in fam}
    noise = {f: np.clip(rng.normal(0, 0.5, G), -1, 1) for f in fam}
    mom = noise if external_right else true
    ext = true if external_right else noise
    base = {lab: rng.uniform(1.5, 2.5, G) for labs in fam.values() for lab in labs}
    fam_of = {lab: f for f, labs in fam.items() for lab in labs}
    e95 = {lab: base[lab] for lab in base}
    e105 = {lab: base[lab] + true[fam_of[lab]] for lab in base}
    e85 = {"OFT/RV-CM": base["OFT/RV-CM"] - mom["CM"], "RV-CM": base["V-CM"] - mom["CM"],
           "Endothelium": base["BEC"] - mom["BEC"], "Blood": base["Blood"] - mom["Blood"],
           "pSHF": base["pSHF"] - mom["SHF"]}
    cp0 = {t: rng.uniform(2, 20, G) for t in ("ventricular_cm", "atrial_cm", "endocardial_ec", "vascular_ec", "blood")}
    efam = {"ventricular_cm": "CM", "atrial_cm": "CM", "endocardial_ec": "Endo", "vascular_ec": "BEC", "blood": "Blood"}
    cps = {}
    for t, c in cp0.items():
        cps[("E9.5", t)] = c
        cps[("E14.5", t)] = c * np.exp(2 * ext[efam[t]])
        cps[("E15.5", t)] = c * np.exp(2 * ext[efam[t]])
    st = {"prev": t1_like(e85), "start": t1_like(e95), "target": t1_like(e105)}
    H = {r: dl.harmonise({k.split("|")[0] for k in s["keys"]}, CFG["harmonise"]) for r, s in st.items()}
    F = feng_like(cps)
    exts, idx = {"feng": F}, {"feng": dl.sym_index(F["genes"], GENES)}
    flat = dl.flat_mask(GENES, ["Xist", "Cwc22"])
    groups = readout.build_groups(CFG, st, H, exts, idx, smoke=False)
    type_rows = readout.per_type_rows(CFG, groups, flat, "yes")
    pop_rows = readout.population_rows(CFG, st, groups, flat, "yes")
    readings = readout.readings_for_decision(CFG, type_rows, pop_rows)
    return groups, type_rows, pop_rows, readings, dl.decide(RULE, readings)


def test_synthetic_external_right_builds():
    groups, type_rows, pop_rows, readings, dec = synthetic_readout(external_right=True)
    assert set(groups) == {"CM", "ventricularCM", "Endocardium", "Endothelium", "Blood", "SHF"}
    assert groups["Endocardium"]["Lp"] == ["Endothelium"]  # E8.5 fallback antecedent
    assert "feng_E14.5" not in groups["SHF"]["ext"]
    assert readings["primary"]["groups"] == "Blood,CM,Endocardium,Endothelium"  # matched coverage, no SHF
    assert readings["cm_sign_agreement"] == pytest.approx(1.0)
    assert dec["verdict"] == "build", dec


def test_synthetic_momentum_right_does_not_build():
    *_, readings, dec = synthetic_readout(external_right=False)
    assert dec["verdict"] == "do not build"
    assert readings["baseline"]["de_score"] > readings["primary"]["de_score"]


def test_cross_fitting_uses_disjoint_halves():
    """With halves that differ only at the start stage, fold A->B and B->A must see different targets."""
    groups, type_rows, *_ = synthetic_readout(external_right=True)
    r = groups["CM"]
    assert r["S_half"]["A"]["n"] + r["S_half"]["B"]["n"] == r["S_all"]["n"]
    folds = {x["target"] for x in type_rows if x["predictor"] == "mom_xfit"}
    assert folds == {"xfit_A", "xfit_B", "xfit_mean"}


def test_held_flat_genes_are_zero_in_predictions():
    flat = dl.flat_mask(GENES, ["Xist", "Cwc22"])
    vec = {"CM": np.ones(len(GENES)), "Blood": np.full(len(GENES), np.nan)}
    pred = dl.population_prediction(vec, {"CM": 0.5, "Blood": 0.5}, ["CM", "Blood"], flat)
    assert pred[flat].tolist() == [0.0, 0.0] and np.allclose(pred[~flat], 0.5)
    v = readout.held(np.ones(len(GENES)), flat)
    assert np.isnan(v[flat]).all() and np.isfinite(v[~flat]).all()


def test_pool_rule_and_held_flat_list():
    G = len(GENES)
    rng = np.random.default_rng(3)
    labs = [f"type{i}" for i in range(10)]
    det_a = {lab: rng.uniform(0.3, 0.6, G) for lab in labs}
    det_b = {lab: det_a[lab].copy() for lab in labs}
    for i, lab in enumerate(labs):
        det_b[lab][0] -= 0.3                      # g0: every type -> flagged
        det_b[lab][2] += 0.2 if i < 9 else 0.0    # g2: 9 of 10 types -> flagged at 0.9
        det_b[lab][3] += 0.1                      # g3: below the 0.15 floor -> kept
    det_b["type0"][1] -= 0.5                      # g1: one type only -> kept
    means = {lab: np.ones(G) for lab in labs}
    a, b = t1_like(means, n=100, det=det_a), t1_like(means, n=100, det=det_b)
    mask, info = dl.pool_effect_genes(a, b, 50, 0.15, 0.9)
    assert GENES[mask].tolist() == ["g0", "g2"] and info["min_types"] == 9
    cfg = {**CFG["held_flat"], "pool_rule": {**CFG["held_flat"]["pool_rule"], "stages": ["x", "y"]}}
    df = dl.held_flat_genes(cfg, GENES, {"x": a, "y": b})
    assert set(df.gene) == {"g0", "g2", "Xist", "Cwc22"}
    assert df.set_index("gene").loc["Xist", "reasons"] == "sex-linked"


def test_align_genes():
    rng = np.random.default_rng(1)
    st = t1_like({"V-CM": rng.uniform(1, 2, len(GENES))})
    perm = rng.permutation(len(GENES))
    shuffled = {**st, "genes": GENES[perm], "S_log": st["S_log"][:, perm], "S_cp": st["S_cp"][:, perm],
                "S_nz": st["S_nz"][:, perm]}
    back, how = readout.align_genes(shuffled, GENES)
    assert how.startswith("reordered") and np.array_equal(back["S_log"], st["S_log"])
    assert readout.align_genes(st, GENES)[1] == "same order"
    with pytest.raises(SystemExit):
        readout.align_genes({**shuffled, "genes": shuffled["genes"][1:], "S_log": shuffled["S_log"][:, 1:],
                             "S_cp": shuffled["S_cp"][:, 1:], "S_nz": shuffled["S_nz"][:, 1:]}, GENES)


# ------------------------------------------------------------------ t1_stats and the command
def test_t1_stats_on_tiny_file(tmp_path):
    rng = np.random.default_rng(0)
    n, G = 40, 12
    counts = rng.poisson(3, (n, G))
    counts[:, 0] = 1  # every cell has a count of 1 somewhere, so depth is recoverable
    lib = counts.sum(1, keepdims=True)
    X = sparse.csr_matrix(np.log1p(counts / lib * 1e4).astype(np.float32))
    obs = {"celltype": np.repeat(["V-CM", "Blood"], n // 2)}
    a = ad.AnnData(X=X, obs=obs)
    a.obs_names = [f"cell{i}-1_{1 + i % 2}" for i in range(n)]
    a.var_names = [f"g{i}" for i in range(G)]
    a.write_h5ad(tmp_path / "t.h5ad")
    st = t1_stats.compute(tmp_path / "t.h5ad", "E0", "celltype", 0)
    assert int(st["n"].sum()) == n and set(st["keys"]) <= {f"{c}|{li}|{h}" for c in ("V-CM", "Blood")
                                                           for li in "12" for h in "AB"}
    assert np.allclose(st["S_log"].sum(0), X.toarray().sum(0), rtol=1e-5)
    assert np.allclose(st["S_cp"].sum(0), (counts / lib * 1e4).sum(0), rtol=1e-4)
    assert np.allclose(st["depth_sum"].sum(), lib.sum(), rtol=1e-4)
    halves = {k.split("|")[2] for k in st["keys"]}
    assert halves == {"A", "B"}
    with pytest.raises(KeyError):
        t1_stats.compute(tmp_path / "t.h5ad", "E0", "cell_type", 0)
    tags, note = t1_stats.library_tags([f"c{i}" for i in range(20)])
    assert set(tags) == {"0"} and "not a library tag" in note


def test_missing_e105_exits_cleanly(tmp_path):
    p = subprocess.run([sys.executable, str(HERE / "readout.py"), str(tmp_path / "E10.5_RNA.h5ad")],
                       capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, p.stderr
    assert "not found" in p.stdout


def test_label_map_outside_the_repo_is_refused(tmp_path):
    """A label map is part of the procedure: a real run refuses one that is not a committed repo file."""
    e105 = tmp_path / "E10.5_RNA.h5ad"
    e105.write_bytes(b"")
    label_map = tmp_path / "label_map.json"
    label_map.write_text("{}")
    p = subprocess.run([sys.executable, str(HERE / "readout.py"), str(e105), "--label-map", str(label_map),
                        "--out", str(tmp_path / "out")], capture_output=True, text=True, timeout=120)
    assert p.returncode == 3, p.stderr
    assert "label_map.json" in p.stderr and not (tmp_path / "out").exists()
