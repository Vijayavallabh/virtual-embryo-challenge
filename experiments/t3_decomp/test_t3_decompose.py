"""Tests for the pre-registered T3 decomposition (experiments/t3_decomp).

Offline and fast: `uv run pytest experiments/t3_decomp -q`. The end-to-end test on synthetic files runs both commands
(mab21l2 and gata4) with one seed and is skipped unless RUN_T3_E2E=1 (run it on the DGX; it takes a few minutes).
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
spec = importlib.util.spec_from_file_location("t3_decompose", HERE / "t3_decompose.py")
t3d = importlib.util.module_from_spec(spec)
sys.modules["t3_decompose"] = t3d
spec.loader.exec_module(t3d)
RULES = yaml.safe_load((HERE / "rules.yaml").read_text())


def _vocab(p, nm, mesp1=None, model_idx=None):
    p = np.asarray(p, float)
    L = len(p)
    labels = [f"L{i}" for i in range(L)]
    mesp1 = np.array(mesp1 if mesp1 is not None else ["yes"] * L)
    return t3d.Vocab(labels=labels, group=np.array(["CM"] * L), mesp1=mesp1, nm=np.asarray(nm, bool),
                     broad=mesp1 != "yes", model_idx=np.arange(L) if model_idx is None else np.asarray(model_idx),
                     p=p / p.sum(), count=(p * 1000).astype(int))


# ------------------------------------------------------------------------------------------------- rules and tables

def test_eval_condition():
    names = {"D2": True, "D10": False, "D7": None}
    assert t3d.eval_condition("D2 and not D10", names)
    assert not t3d.eval_condition("D2 and D10", names)
    assert not t3d.eval_condition("D7", names)          # not evaluated counts as false
    assert t3d.eval_condition("true", names) and t3d.eval_condition(True, names)
    with pytest.raises(KeyError):
        t3d.eval_condition("D99", names)
    with pytest.raises(ValueError):
        t3d.eval_condition("__import__('os')", names)
    with pytest.raises(ValueError):
        t3d.eval_condition("D2 == 1", names)


def _tables():
    sub = RULES["submissions"]
    for branch in ("allowed", "refused"):
        for key, tab in sub[branch].items():
            if isinstance(tab, list):
                yield branch, key, tab


def test_decision_tables_complete_and_named():
    known = set(RULES["rules"]) | {"m_raw_positive"}
    for branch, key, tab in _tables():
        assert str(tab[-1]["when"]).strip().lower() == "true", (branch, key)
        for row in tab:
            for name in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", str(row["when"])):
                assert name in known | {"and", "or", "not", "true", "false"}, (branch, key, name)
    # every combination of outcomes picks a row in every table
    ids = sorted(known)
    for combo in itertools.product([True, False, None], repeat=4):
        names = dict.fromkeys(ids, False)
        names.update(dict(zip(["D2", "D2u", "D10", "D7"], combo, strict=True)))
        for _, _, tab in _tables():
            assert t3d.first_match(tab, names) is not None


def test_factor_sources_and_compositions_are_defined():
    sub = RULES["submissions"]
    for branch in ("allowed", "refused"):
        b = sub[branch]
        for row in b["factor_table"]:
            for s in ("s1", "s2"):
                assert row[s] in b["factor_exponent"], (branch, row[s])
        for key, tab in b.items():
            if key.startswith("composition_s") and isinstance(tab, list):
                for row in tab:
                    assert row["value"] in ("none", *RULES["weights"]["composition"]), row


def test_rules_name_no_panel_gene():
    """Provenance: the pre-registration names cell states and lineage groups, never a gene on the T3 panel."""
    genes_file = REPO / "boards" / "T3__gata4.genes.txt"
    if not genes_file.exists():
        pytest.skip("panel gene list not present")
    genes = {g.strip() for g in genes_file.read_text().splitlines() if g.strip()}
    genotypes = {"Gata4", "Gata6", "Mab21l2", "Ctnnb1"}       # the knockouts themselves (presets, genotype text)
    text = (HERE / "rules.yaml").read_text()
    tokens = set(re.findall(r"[A-Za-z][A-Za-z0-9-]*", text))
    assert not (tokens & (genes - genotypes)), sorted(tokens & (genes - genotypes))


def test_scenario_matches_config():
    out = t3d.check_scenario(RULES)
    assert out["depleted_states"] == RULES["presets"]["gata4"]["scenario"]["depleted_states"]


# ------------------------------------------------------------------------------------------------- count space

def test_gene_factor_normalised_and_clipped():
    rng = np.random.default_rng(0)
    mu_wt = rng.dirichlet(np.ones(50), size=4)
    mu_ko = mu_wt * rng.lognormal(0, 0.3, size=50)
    mu_ko /= mu_ko.sum(1, keepdims=True)
    w = np.array([0.4, 0.3, 0.2, 0.1])
    phi, den = t3d.gene_factor(mu_ko, mu_wt, w, np.array([1, 1, 0, 1]), eps=0.0, clip=(1e-9, 1e9))
    assert np.isclose(den @ phi / den.sum(), 1.0)
    phi_c, _ = t3d.gene_factor(mu_ko * np.r_[100.0, np.ones(49)], mu_wt, w, np.ones(4))
    assert phi_c.max() <= 4.0 and phi_c.min() >= 0.25
    # no label in the mask: the factor is neutral
    assert np.all(t3d.gene_factor(mu_ko, mu_wt, w, np.zeros(4))[0] == 1.0)


def test_apply_gene_factor_expectation():
    rng = np.random.default_rng(1)
    C = rng.poisson(20, size=(4000, 3)).astype(np.int32)
    P = np.full(C.shape, 1 / 3, dtype=np.float32)
    phi = np.array([0.5, 1.0, 2.0])
    out = t3d.apply_gene_factor(C, phi, P, rng)
    ratio = out.mean(0) / C.mean(0)
    assert abs(ratio[0] - 0.5) < 0.02 and abs(ratio[1] - 1.0) < 1e-9
    assert abs(ratio[2] - (1 + 1.0 * (C.sum(1).mean() / 3) / C[:, 2].mean())) < 0.05  # top-up from P


def test_rescale_depth_and_thin_never_empty():
    rng = np.random.default_rng(2)
    C = rng.poisson(1.0, size=(500, 10)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    P = np.full(C.shape, 0.1, dtype=np.float32)
    out = t3d.rescale_depth(C, np.full(500, 0.01), P, rng)
    assert (out.sum(1) > 0).all()
    up = t3d.rescale_depth(C, np.full(500, 3.0), P, rng)
    assert up.sum() > 2.0 * C.sum()


def test_alloc_and_lognorm():
    k = t3d.alloc(1001, [0.5, 0.3, 0.2])
    assert k.sum() == 1001 and list(k) == [501, 300, 200]
    from vec.prep import log1p_cp10k

    C = np.random.default_rng(3).poisson(3, size=(20, 7))
    C[C.sum(1) == 0, 0] = 1
    assert np.allclose(t3d.lognorm(C), log1p_cp10k(C).toarray(), atol=1e-6)


# ------------------------------------------------------------------------------------------------- compositions

def test_bias_matched_identity_and_shift():
    V = _vocab([0.5, 0.3, 0.2], nm=[False, False, True])
    q = np.array([0.4, 0.35, 0.25])
    assert np.allclose(t3d.bias_matched(V, q, q), V.p)           # no change in calls -> WT composition
    f = t3d.bias_matched(V, q, np.array([0.5, 0.3, 0.2]))
    assert np.isclose(f.sum(), 1.0) and np.allclose(f, q)


def test_bias_matched_keeps_labels_outside_the_classifier():
    V = _vocab([0.4, 0.4, 0.2], nm=[False, True, False], model_idx=[0, 1])
    f = t3d.bias_matched(V, np.array([0.75, 0.25]), np.array([0.5, 0.5]))
    assert np.isclose(f[2], 0.2) and np.isclose(f.sum(), 1.0)


def test_coverage_and_mesp1_compose_to_full():
    V = _vocab([0.3, 0.3, 0.2, 0.2], nm=[False, False, True, True])
    f = np.array([0.4, 0.2, 0.1, 0.3])
    cov, m1 = t3d.coverage_comp(V, f), t3d.mesp1_comp(V, f)
    assert np.allclose(cov[2:], f[2:]) and np.allclose(cov[:2] / cov[:2].sum(), V.p[:2] / V.p[:2].sum())
    assert np.allclose(m1[2:], V.p[2:]) and np.allclose(m1[:2] / m1[:2].sum(), f[:2] / f[:2].sum())
    assert np.isclose(cov.sum(), 1) and np.isclose(m1.sum(), 1)


def test_rule_comp_depletion_and_magnitude():
    V = _vocab([0.5, 0.25, 0.25], nm=[False, False, False])
    V.labels = ["A", "B", "C"]
    f, missing = t3d.rule_comp(V, {"deplete": {"A": 0.5, "Z": 0.5}}, {})
    assert missing == ["Z"] and np.allclose(f, np.array([0.25, 0.25, 0.25]) / 0.75)
    g, _ = t3d.rule_comp(V, {"deplete": {"A": 0.0}}, {}, log2_magnitude=1.0)
    assert np.allclose(g, f)
    h, _ = t3d.rule_comp(V, {"deplete_list": "cm", "fraction": 0.5, "enrich": {"C": 2.0}}, {"cm": ["A"]})
    assert np.allclose(h, np.array([0.25, 0.25, 0.5]) / 1.0)


# ------------------------------------------------------------------------------------------------- emulator, factorial

def test_working_copy_rows_match_vec_evaluate():
    from vec.evaluate import Stage, subsample

    X = np.arange(5029 * 2, dtype=np.float32).reshape(-1, 2)
    for seed in (0, 1, 7):
        rows = t3d.working_copy_rows(len(X), 0.1, seed)
        assert np.array_equal(subsample(Stage(X=X), 0.1, seed).X, X[rows])


def test_fac_names_roundtrip():
    for e in ("top", "em"):
        for on in [(), ("factor",), ("cov", "depth"), ("mesp1", "stage", "factor"), t3d.COMPONENTS]:
            name = t3d.fac_name(e, on)
            est, comps = t3d.fac_parse(name)
            assert set(comps) == set(on)
            assert (est == e) == bool({"cov", "mesp1"} & set(on))
    cells = t3d.factorial_cells(t3d.COMPONENTS)
    assert len({t3d.fac_name(e or "top", on) for e, on in cells}) == 8 + 2 * 24


def test_shapley_additive_game():
    comps = t3d.COMPONENTS
    contrib = dict(zip(comps, [1.0, 2.0, -0.5, 0.25, 10.0], strict=True))
    rows = []
    for seed in (0, 1):
        for e, on in t3d.factorial_cells(comps):
            val = sum(contrib[c] for c in on) + (0.1 if e == "em" else 0.0) * ({"cov", "mesp1"} <= set(on))
            rows.append({"target": "x", "seed": seed, "variant": t3d.fac_name(e or "top", on), "t3": val,
                         "de_direction": val, "de_score": val, "skill_severity_slope": val, "skill_mmd_u": val,
                         "skill_variogram": val})
    sh = t3d.shapley_table(pd.DataFrame(rows), comps)
    top = sh[(sh["estimator"] == "top") & (sh["metric"] == "t3")].set_index("component")
    for c in comps:
        assert np.isclose(top.at[c, "shapley"], contrib[c])
        assert np.isclose(top.at[c, "alone"], contrib[c]) and np.isclose(top.at[c, "last"], contrib[c])
    # an interaction is split equally between its two players (the em game adds 0.1 when cov and mesp1 are both on)
    em = sh[(sh["estimator"] == "em") & (sh["metric"] == "t3")].set_index("component")
    assert np.isclose(em.at["cov", "shapley"], contrib["cov"] + 0.05)
    assert np.isclose(em["shapley"].sum(), sum(contrib.values()) + 0.1)


def test_row_hashes_match_across_storage():
    """An as-released file may store counts as float32 (interpolated records); its integer rows must hash like the
    deduplicated pool's int32 rows, so copies of working-copy cells are found."""
    import anndata as ad
    from scipy import sparse

    C = np.array([[0, 3, 1], [2, 0, 0], [0, 3, 1]], dtype=np.int32)
    a_int = ad.AnnData(X=np.zeros((3, 3), dtype=np.float32))
    a_int.layers["counts"] = sparse.csr_matrix(C)
    a_flt = ad.AnnData(X=np.zeros((4, 3), dtype=np.float32))
    a_flt.layers["counts"] = sparse.csr_matrix(np.vstack([C, [[0.5, 1.5, 0.0]]]).astype(np.float32))
    h_int, h_flt = t3d.row_hashes_layer(a_int), t3d.row_hashes_layer(a_flt)
    assert np.array_equal(h_int, h_flt[:3]) and h_int[0] == h_int[2] and h_flt[3] not in set(h_int)
    assert np.array_equal(t3d.dense_counts(a_flt)[:3], C)


def test_stable_rng_independent_of_order():
    a = t3d.stable_rng(3, "fac:none").integers(0, 10**9, 5)
    _ = t3d.stable_rng(3, "fac:factor").integers(0, 10**9, 5)
    b = t3d.stable_rng(3, "fac:none").integers(0, 10**9, 5)
    assert np.array_equal(a, b)


# ------------------------------------------------------------------------------------------------- command line

def test_gata4_exits_cleanly_without_files(tmp_path, capsys):
    rc = t3d.main(["gata4", str(tmp_path / "rep1.h5ad"), str(tmp_path / "rep2.h5ad"), str(tmp_path / "wt.h5ad"),
                   "--out", str(tmp_path / "out")])
    err = capsys.readouterr().err
    assert rc == 0 and "nothing was run" in err and "rep1.h5ad" in err
    assert not (tmp_path / "out").exists()


# ------------------------------------------------------------------------------------------------- synthetic end to end

def _write_stage(path, shares, prof, n, rng, depth=120, cm=None, factor=None, n_levels=20, counts=True,
                 condition=None):
    """A synthetic MERFISH-like file: labels drawn from `shares`, Poisson counts around per-label profiles."""
    import anndata as ad
    from scipy import sparse

    from vec.prep import log1p_cp10k

    names = list(shares)
    w = np.array([shares[k] for k in names], float)
    lab = rng.choice(len(names), n, p=w / w.sum())
    rows = []
    for i in lab:
        p = prof[names[i]] * (1.0 if factor is None else factor)
        rows.append(rng.poisson(p / p.sum() * rng.lognormal(np.log(depth), 0.3)))
    C = np.asarray(rows, dtype=np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [names[i] for i in lab]}, index=[f"{path.stem}_{i}" for i in range(n)])
    if cm is not None:
        obs["cm_celltype"] = [cm.get(names[i], "Unknown") for i in lab]
    if condition is not None:
        obs["condition"] = condition(n, rng)
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=[f"G{j}" for j in range(C.shape[1])]))
    if counts:
        a.layers["counts"] = sparse.csr_matrix(C)
    z = rng.integers(0, n_levels, n) * 10.0
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 500, n), rng.uniform(0, 500, n), z]
    a.write_h5ad(path)


@pytest.mark.skipif(os.environ.get("RUN_T3_E2E") != "1", reason="set RUN_T3_E2E=1 (run on the DGX)")
def test_end_to_end_synthetic(tmp_path):
    rng = np.random.default_rng(0)
    G = 60
    wt95 = ["V-CM", "A-CM", "Peri", "SE", "Foregut", "NCC", "Branch Arch", "ST", "Proepi", "Hepatocytes",
            "Dorsal-PAM", "Lateral-PAM", "Chamber-Endo", "pPHM"]
    e875_ct = ["V-CM", "aPHM", "pPHM", "JCF", "Neural Tube", "NCC", "d-CSE", "D-FG", "a-FG", "Endo", "PAM-1",
               "ExEM-1", "Peri", "Intra-Endoth-1"]
    allnames = sorted(set(wt95) | set(e875_ct) | {"Unknown"})
    profiles = pd.DataFrame(rng.lognormal(0, 1.2, size=(len(allnames), G)), index=allnames,
                            columns=[f"G{j}" for j in range(G)])
    prof = {k: profiles.loc[k].to_numpy() for k in allnames}
    d = tmp_path / "raw"
    (d / "T2").mkdir(parents=True)
    (d / "T3").mkdir()
    share95 = dict.fromkeys(wt95, 1.0)
    _write_stage(d / "T2" / "E9.5.h5ad", share95, prof, 2400, rng, depth=80)
    share875 = dict.fromkeys(e875_ct, 1.0)
    cm = {"V-CM": "vCM1", "aPHM": "aSHF", "pPHM": "pSHF"}
    _write_stage(d / "T2" / "E8.75.h5ad", share875, prof, 2400, rng, depth=130, cm=cm)
    _write_stage(d / "T2" / "E8.25_late.h5ad", share875, prof, 2000, rng, depth=220, cm=cm)
    ko_factor = rng.lognormal(0, 0.4, G)
    ko_share = {**share95, "ST": 2.0, "V-CM": 1.8}
    _write_stage(d / "T3" / "E9.5_mab21l2_ko.h5ad", ko_share, prof, 2400, rng, depth=95, factor=ko_factor)
    g4_share = {**share875, "V-CM": 0.5}
    _write_stage(d / "rep1.h5ad", g4_share, prof, 2000, rng, depth=110, cm=cm, factor=ko_factor,
                 condition=lambda n, r: np.where(r.random(n) < 0.3, "WT", "Gata4 cKO"))
    _write_stage(d / "rep2.h5ad", g4_share, prof, 2000, rng, depth=110, cm=cm, factor=ko_factor, counts=False)
    out1 = tmp_path / "out_mab"
    rc = t3d.main(["mab21l2", "--wt", str(d / "T2" / "E9.5.h5ad"), "--ko", str(d / "T3" / "E9.5_mab21l2_ko.h5ad"),
                   "--e875", str(d / "T2" / "E8.75.h5ad"), "--e825", str(d / "T2" / "E8.25_late.h5ad"),
                   "--out", str(out1), "--seeds", "0", "--workers", "1", "--threads", "2"])
    assert rc == 0
    dec = json.loads((out1 / "mab21l2_decisions.json").read_text())
    assert dec["protocol_as_registered"] is False and "D3" in dec["rules"]
    summ = pd.read_csv(out1 / "mab21l2_summary.csv")
    assert {"fac:none", "fac:factor", "fac[em]:cov+mesp1", "R"} <= set(summ["variant"])
    out2 = tmp_path / "out_g4"
    rc = t3d.main(["gata4", str(d / "rep1.h5ad"), str(d / "rep2.h5ad"), str(d / "T2" / "E8.75.h5ad"),
                   "--e95", str(d / "T2" / "E9.5.h5ad"), "--e825", str(d / "T2" / "E8.25_late.h5ad"),
                   "--mab21l2", str(d / "T3" / "E9.5_mab21l2_ko.h5ad"), "--out", str(out2), "--seeds", "0",
                   "--workers", "2", "--threads", "1"])
    assert rc == 0
    dec = json.loads((out2 / "gata4_decisions.json").read_text())
    for branch in ("allowed", "refused"):
        assert {"s1", "s2"} <= set(dec["submissions"][branch])
    assert dec["rules"]["D2"]["outcome"] in (True, False)
    meta = json.loads((out2 / "gata4_meta.json").read_text())
    assert meta["targets"]["rep1"]["littermate_wt_cells"] > 0
