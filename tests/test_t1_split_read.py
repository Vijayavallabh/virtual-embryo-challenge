"""experiments/t1_m1/split_read.py (try T51 of experiments/harness/TRIES.md) on synthetic Feng statistics and synthetic
AnnData stages built here. No released file, statistics file, committed result or score is read."""

from __future__ import annotations

import importlib.util
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec.boards import REPO_ROOT

SCRIPT = REPO_ROOT / "experiments" / "t1_m1" / "split_read.py"
EPS = 0.5
N_CELLS = 21  # per label and library: above the config's min_cells_t1
FENG_LABELS = ("ventricular_cm", "endocardial_ec", "vascular_ec", "epicardial")
T1_LABEL = {"CM": "V-CM", "Endocardium": "Endocardium", "Endothelium": "BEC", "Epicardium": "Epicardium"}
G = 120


@pytest.fixture
def sp(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("split_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    return mod


@pytest.fixture
def cfg(sp):
    return sp.feng_dirs.load_config()


# ---------------------------------------------------------------------------------------------------------------------
# Synthetic inputs


def split_dirs(seed=0, n=G - 1, orth=True, sd_c=0.3, sd_s=0.3):
    """A common part c0 and four specific parts S (rows, summing to zero over the groups), every one centred over the
    genes; with orth each row of S is orthogonal to c0 (so d_g = c0 + S_g has its share along c and s_g alike)."""
    rng = np.random.default_rng(seed)
    c0 = rng.normal(0, sd_c, n)
    c0 -= c0.mean()
    S = rng.normal(0, sd_s, (len(FENG_LABELS), n))
    S -= S.mean(0)
    S -= S.mean(1, keepdims=True)
    if orth:
        S -= np.outer(S @ c0 / (c0 @ c0), c0)
    return c0, S


def feng_world(c0, S, seed=0, half_noise=0.0):
    """Feng-format statistics at E9.5 and E14.5 on genes g1..g{G-1} (g0 absent) whose group directions are exactly
    D_g = c0 + S_g (groups in T1_LABEL's order). ventricular_cm (CM), endocardial_ec (Endocardium) and epicardial
    (Epicardium): E9.5 means a, E14.5 means (a + EPS) exp(D_g) - EPS. Endothelium pools endocardial_ec and vascular_ec
    (the config's labels): vascular_ec has a hundred times endocardial_ec's cells, the same E9.5 means, and E14.5
    means chosen so that the pool's are (a + EPS) exp(D_Endothelium) - EPS. With half_noise, each group's two E14.5
    halves carry independent noise (sd half_noise) in the exponent, so every group's direction has the same noise; the
    pooled direction is then the log of the halves' average."""
    rng = np.random.default_rng(seed + 100)
    genes = [f"g{i}" for i in range(G)]
    D = {g: c0 + S[i] for i, g in enumerate(T1_LABEL)}
    n_e, n_v = 30, 3000  # cells per half: endocardial_ec and every other label, vascular_ec

    def noise():
        return rng.normal(0, half_noise, G - 1) if half_noise else 0.0

    a = {lab: rng.uniform(2.0, 6.0, G - 1) for lab in ("ventricular_cm", "endocardial_ec", "epicardial")}
    a["vascular_ec"] = a["endocardial_ec"]
    e14 = {}
    for lab, g in (("ventricular_cm", "CM"), ("endocardial_ec", "Endocardium"), ("epicardial", "Epicardium")):
        e14[lab] = {h: (a[lab] + EPS) * np.exp(D[g] + noise()) - EPS for h in "AB"}
    e14["vascular_ec"] = {}
    for h in "AB":
        pool = (a["endocardial_ec"] + EPS) * np.exp(D["Endothelium"] + noise()) - EPS
        e14["vascular_ec"][h] = ((n_e + n_v) * pool - n_e * e14["endocardial_ec"][h]) / n_v
        assert np.all(e14["vascular_ec"][h] > 0)
    keys, ns, Ssum = [], [], []
    for st in ("E9.5", "E14.5"):
        for lab in FENG_LABELS:
            n = n_v if lab == "vascular_ec" else n_e
            for h in "AB":
                keys.append(f"{st}|z|{lab}|{h}")
                ns.append(n)
                Ssum.append((a[lab] if st == "E9.5" else e14[lab][h]) * n)
    Ssum = np.array(Ssum, dtype=np.float64)
    ns = np.array(ns, float)
    fst = {"genes": np.array(genes[1:], dtype=object), "keys": np.array(keys, dtype=object), "n": ns,
           "S_log": np.log1p(Ssum / ns[:, None]) * ns[:, None], "S_cp": Ssum, "S_nz": np.repeat(ns[:, None], G - 1, 1),
           "depth_sum": 1e4 * ns}
    return genes, fst, a


def stage(blocks: dict, genes) -> ad.AnnData:
    """blocks[(label, library)] = mean CP10k vector (on all genes); N_CELLS identical cells each."""
    rows, labs, names = [], [], []
    for (lab, lib), cp in blocks.items():
        for _ in range(N_CELLS):
            rows.append(np.log1p(cp))
            labs.append(lab)
            names.append(f"c{len(names)}_{lib}")
    return ad.AnnData(X=sparse.csr_matrix(np.array(rows, dtype=np.float32)),
                      obs=pd.DataFrame({"celltype": labs}, index=names), var=pd.DataFrame(index=list(genes)))


def counts_stage(blocks: dict, genes, rng, n=60, depth=200.0) -> ad.AnnData:
    """blocks[(label, library)] = mean vector; n cells of Poisson counts (mean x depth) each, X = log1p CP10k of them
    (counts recoverable: gene g0 has a count of one in every cell)."""
    from vec import t1_scope

    rows, labs, names = [], [], []
    for (lab, lib), mean in blocks.items():
        rows.append(rng.poisson(np.maximum(mean, 0.0) * depth, size=(n, len(genes))))
        labs += [lab] * n
        names += [f"c{len(names) + i}_{lib}" for i in range(n)]
    C = np.vstack(rows).astype(np.int64)
    C[:, 0] = 1
    return ad.AnnData(X=t1_scope.log1p_cp10k(sparse.csr_matrix(C)), obs=pd.DataFrame({"celltype": labs}, index=names),
                      var=pd.DataFrame(index=list(genes)))


def t1_world(sp, cfg, genes, fst, change, counts=False, seed=0):
    """Three T1 stages: E9.5 at Feng's E9.5 means (libraries 1 and 2); E10.5 libraries 1, 2 and 3 moved by
    change[group] x the library's factor (1, 1, 2); E8.5 moved by -0.5 change. change[group] is a vector on the T1 genes
    (0 on g0). Identical cells: moved exactly in log(cp + EPS). Counts: moved multiplicatively in count space (their
    CP10k are about twenty times the means, so the pseudocount damps a change by 1-2%)."""
    _, _, idx = sp.feng_dirs.directions_from_stats(fst, genes, cfg, sp.sr.BRACKET)
    e95 = {}
    for g in T1_LABEL:
        v = sp.dl.on_t1(sp.dl.feng_sel(fst, "E9.5", cfg["external_labels"][g][0])["cp"], idx)
        v[~np.isfinite(v)] = 5.0
        e95[g] = v

    def moved(g, k):  # exact in log(cp + EPS) for identical cells; multiplicative (count space) for counts
        return e95[g] * np.exp(k * change[g]) if counts else (e95[g] + EPS) * np.exp(k * change[g]) - EPS

    E95 = {(T1_LABEL[g], li): e95[g] for g in T1_LABEL for li in ("1", "2")}
    E105 = {(T1_LABEL[g], li): moved(g, k) for g in T1_LABEL for li, k in (("1", 1.0), ("2", 1.0), ("3", 2.0))}
    E85 = {(T1_LABEL[g], "1"): moved(g, -0.5) for g in T1_LABEL}
    if not counts:
        return {"E8.5": stage(E85, genes), "E9.5": stage(E95, genes), "E10.5": stage(E105, genes)}
    rng = np.random.default_rng(seed + 7)
    return {s: counts_stage(b, genes, rng) for s, b in (("E8.5", E85), ("E9.5", E95), ("E10.5", E105))}


def on_t1(vec):
    return np.concatenate([[0.0], vec])


def group_change(c0, S, a_c, a_s):
    """Per T1 group (Feng label order), the change a_c c0 + a_s S_label on the T1 genes."""
    return {g: on_t1(a_c * c0 + a_s * S[i]) for i, g in enumerate(T1_LABEL)}


def t47_of(sp, stages, fst, cfg, held=()):
    """T47's record of the same world, from share_read.reading (as t47_record reads the committed json)."""
    r = sp.sr.reading(stages, fst, list(held), {}, cfg, n_boot=1)
    return {"A": r["reading"]["A"]["value"], "L": r["reading"]["L"]["value"],
            "L_corr": r["reading"]["L_corr"]["value"], "support": r["support"]["A_K_support"]}


# ---------------------------------------------------------------------------------------------------------------------
# The split


def test_c_and_s_add_up_to_d_on_M_and_the_specific_parts_sum_to_zero(sp):
    rng = np.random.default_rng(1)
    dirs = {g: rng.normal(size=50) for g in "abcd"}
    dirs["b"][3] = np.nan
    masks = {g: np.isfinite(v) for g, v in dirs.items()}
    masks["c"][7] = False  # a held gene in one group's mask
    M = sp.common_mask(masks)
    assert not M[3] and not M[7] and M.sum() == 48
    sv = sp.split_vectors(dirs, M, {"a": 3.0, "b": 1.0, "c": 0.0, "d": 2.0})
    for g in dirs:
        assert np.allclose(sv["c"][M] + sv["s"][g][M], dirs[g][M])
        assert np.allclose(sv["c_w"][M] + sv["s_w"][g][M], dirs[g][M])
        assert np.all(np.isnan(sv["s"][g][~M]))
    assert np.allclose(sum(sv["s"][g][M] for g in dirs), 0.0, atol=1e-12)
    assert np.allclose(sv["c_w"][M], (3 * dirs["a"][M] + dirs["b"][M] + 2 * dirs["d"][M]) / 6)  # zero weight out
    e = sp.energy_parts(dirs["a"], sv["c"], sv["s"]["a"], M)
    assert e["e_c"] + e["e_s"] + e["cross"] == pytest.approx(1.0)


def test_l_eff_by_hand(sp):
    assert sp.l_eff(1.2, 0.5, 0.6, 0.4) == pytest.approx(1.2 * (0.6 + 0.5 * 0.4) / 1.0)
    assert sp.l_eff(1.2, 1.7, 0.6, 0.4) == pytest.approx(1.2)  # min(1, R_corr)
    assert np.isnan(sp.l_eff(1.2, float("nan"), 0.6, 0.4))


def test_reliabilities_by_hand_from_the_halves(sp):
    rng = np.random.default_rng(2)
    groups = ["a", "b", "c"]
    dt = {g: rng.normal(size=80) for g in groups}
    halves = {g: {"d_a": dt[g] + rng.normal(0, 0.3, 80), "d_b": dt[g] + rng.normal(0, 0.3, 80)} for g in groups}
    d = {g: (halves[g]["d_a"] + halves[g]["d_b"]) / 2 for g in groups}
    M = np.ones(80, bool)
    w = dict.fromkeys(groups, 1.0)
    sv = sp.split_vectors(d, M, w)
    vecs = {"M": d, "c": dict.fromkeys(groups, sv["c"]), "s": sv["s"], "c_w": dict.fromkeys(groups, sv["c_w"]),
            "s_w": sv["s_w"]}
    rel = sp.reliabilities(vecs, sp.half_vectors(halves, groups, M, w), M)
    ca = np.mean([halves[g]["d_a"] for g in groups], 0)
    cb = np.mean([halves[g]["d_b"] for g in groups], 0)
    c = sv["c"]
    hand = (np.mean(c ** 2) - np.mean(((ca - cb) / 2) ** 2)) / np.mean(c ** 2)
    assert rel["c"]["a"] == pytest.approx(hand) and rel["c"]["b"] == rel["c"]["a"]
    s, sa, sb = sv["s"]["b"], halves["b"]["d_a"] - ca, halves["b"]["d_b"] - cb
    assert rel["s"]["b"] == pytest.approx((np.mean(s ** 2) - np.mean(((sa - sb) / 2) ** 2)) / np.mean(s ** 2))
    # a group without halves leaves c's and every s_g's reliability undefined, and its own r_d
    rel2 = sp.reliabilities(vecs, sp.half_vectors({g: halves[g] for g in ("a", "b")}, groups, M, w), M)
    assert np.isnan(rel2["c"]["a"]) and np.isnan(rel2["s"]["a"]) and np.isnan(rel2["M"]["c"])
    assert rel2["M"]["a"] == rel["M"]["a"]


def test_c_is_more_reliable_than_s_with_equal_noise_per_group(sp):
    """The bias the correction removes: with four groups and equal noise in every group's halves, c (a mean of four)
    carries a quarter of a direction's noise variance and s_g three quarters; with c and s_g of equal signal per group
    (0.09 and 0.0675 per gene) and noise 0.045 per gene in d, r_c is 0.89 and r_s 0.67 (seeds 3-7: 0.888-0.893 and
    0.661-0.677), so the raw R reads about 1.33 where the truth is 1."""
    rng = np.random.default_rng(1003)
    groups = ["a", "b", "c", "d"]
    c0, S = split_dirs(seed=3, n=4000, sd_c=0.3, sd_s=0.3)
    dt = {g: c0 + S[i] for i, g in enumerate(groups)}
    halves = {g: {"d_a": dt[g] + rng.normal(0, 0.3, 4000), "d_b": dt[g] + rng.normal(0, 0.3, 4000)} for g in groups}
    d = {g: (halves[g]["d_a"] + halves[g]["d_b"]) / 2 for g in groups}
    M = np.ones(4000, bool)
    w = dict.fromkeys(groups, 1.0)
    sv = sp.split_vectors(d, M, w)
    vecs = {"M": d, "c": dict.fromkeys(groups, sv["c"]), "s": sv["s"], "c_w": dict.fromkeys(groups, sv["c_w"]),
            "s_w": sv["s_w"]}
    rel = sp.reliabilities(vecs, sp.half_vectors(halves, groups, M, w), M)
    r_c, r_s = rel["c"]["a"], np.mean([rel["s"][g] for g in groups])
    # noise variance per gene in d: 0.3^2 / 2 = 0.045; in c: 0.045 / 4; in s_g: 0.045 * 3 / 4
    assert r_c == pytest.approx(0.09 / (0.09 + 0.045 / 4), abs=0.02)
    assert r_s == pytest.approx(0.0675 / (0.0675 + 0.045 * 0.75), abs=0.02)
    assert r_c / r_s == pytest.approx(1.33, abs=0.05)


# ---------------------------------------------------------------------------------------------------------------------
# Worlds with a known answer, through the whole read


def _inputs(sp, cfg, c0, S, seed=0):
    genes, fst, _ = feng_world(c0, S, seed=seed)
    stages = t1_world(sp, cfg, genes, fst, group_change(c0, S, 0.25, 0.25), seed=seed)
    return stages, fst, []


def _read(sp, cfg, c0, S, a_c, a_s, counts=False, n_boot=2, self_check=False, half_noise=0.0, seed=0, **kw):
    genes, fst, _ = feng_world(c0, S, seed=seed, half_noise=half_noise)
    stages = t1_world(sp, cfg, genes, fst, group_change(c0, S, a_c, a_s), counts=counts, seed=seed)
    t47 = t47_of(sp, stages, fst, cfg)
    return sp.split_reading(stages, fst, [], cfg, t47, n_boot=n_boot, self_check=self_check, **kw), t47


def test_a_world_moved_along_d_reads_its_share_along_c_and_s(sp, cfg):
    c0, S = split_dirs(seed=4)
    r, t47 = _read(sp, cfg, c0, S, 0.25, 0.25)
    rd = r["reading"]
    for k in ("A_M", "A_c", "A_s", "A_M_T47_mask"):
        assert rd[k]["value"] == pytest.approx(0.25, abs=1e-6)
    assert rd["R"]["value"] == pytest.approx(1.0, abs=1e-5)
    assert rd["A_M"]["value"] == pytest.approx(t47["A"], abs=1e-9)  # the mask is T47's here (no held gene)
    setup, _ = sp.build_setup(*_inputs(sp, cfg, c0, S), cfg, t47)
    for i, g in enumerate(T1_LABEL):  # the world's directions are c0 + S_g in T1_LABEL's order
        assert np.allclose(setup["dirs"][g][1:], c0 + S[i], atol=1e-9)
        assert np.allclose(setup["vecs"]["c"][g][1:], c0, atol=1e-9)
    assert rd["K_c"]["value"] == pytest.approx(0.25, abs=1e-6)  # library 3 moves twice as far
    assert rd["P_c"]["value"] == pytest.approx(0.125, abs=1e-6)  # E8.5 -> E9.5 is half a step
    # identical Feng halves: every r_v is 1 and the corrected values are the raw ones
    assert all(v == pytest.approx(1.0) for v in r["reliability"]["c"].values())
    assert rd["R_corr"]["value"] == pytest.approx(rd["R"]["value"])
    # energy: c and s_g orthogonal, so the cross part is 0
    assert r["energy"]["unweighted"]["cross"] == pytest.approx(0.0, abs=1e-9)
    assert r["energy"]["unweighted"]["e_c"] + r["energy"]["unweighted"]["e_s"] == pytest.approx(1.0)
    # L_eff = L where R_corr = 1; the courses from A_M and L_eff
    assert rd["L_eff"]["value"] == pytest.approx(t47["L"], rel=1e-5)
    c = sp.sr.courses(rd["A_M"]["value"], rd["L_eff"]["value"])
    assert rd["f_even"]["value"] == pytest.approx(c["f_even"]) and rd["f_decay"]["value"] == pytest.approx(c["f_decay"])
    # without the self-check the reading is void
    assert r["statement"]["statement"].startswith("void") and "the self-check did not run" in r["statement"]["void"]


def test_a_world_moved_along_c_only_reads_A_s_near_zero(sp, cfg):
    c0, S = split_dirs(seed=5)
    r, _ = _read(sp, cfg, c0, S, 0.25, 0.0, counts=True, n_boot=20, self_check=True)
    rd = r["reading"]
    assert rd["A_c"]["value"] == pytest.approx(0.25, abs=0.02)
    assert abs(rd["A_s"]["value"]) < 0.02
    st = r["statement"]
    # the rule's outcome: R is undefined wherever A_s <= 0 (here at the point and in about half the draws)
    assert rd["R_corr"]["nonfinite_draws"] > 0
    assert st["void"] == [] and st["statement"].startswith(sp.NOT_COMPUTABLE)
    assert "of 20 draws undefined" in st["statement"]


def test_a_world_whose_common_part_moves_at_half_reads_less(sp, cfg):
    c0, S = split_dirs(seed=6)
    r, t47 = _read(sp, cfg, c0, S, 0.15, 0.3, counts=True, n_boot=20, self_check=True)
    rd, st = r["reading"], r["statement"]
    assert rd["R"]["value"] == pytest.approx(0.5, abs=0.05)
    assert rd["R_corr"]["value"] == pytest.approx(0.5, abs=0.05)
    assert st["void"] == [] and st["statement"] == sp.LESS and st["R"] == st["R_corr"] == sp.LESS
    sc = r["self_check"]
    assert sc["known"]["R"] == 1.0 and sc["measured"]["R"] == pytest.approx(1.0, abs=0.1)
    assert sc["bias"]["R"] == sc["measured"]["R"] - 1.0
    assert sc["measured"]["R_corr"] == sc["measured"]["R"]  # every r_v = 1 in the world
    assert sc["known_by_definition"]["R"] == pytest.approx(1.0, abs=0.05)  # orthogonal parts: near the brief's 1
    txt = sp.text({"label": sp.LABEL_T51, "read": r, "stamp": {"commit": "abc", "tree_clean": True, "utc": "t"}})
    assert f"Reading: {sp.LESS}" in txt and all(f"Limit: {n}" in txt for n in sp.LEDGER_NOTES)
    json.dumps(r)  # the record serialises


def test_non_orthogonal_parts_move_the_known_R_from_one(sp, cfg):
    """A change made exactly along d_g reads 1 + <s_g, c> / <c, c> times its share along c: the self-check's
    definitional known values differ from the brief's 1 by the cross terms."""
    c0, S = split_dirs(seed=7, orth=False)
    genes, fst, _ = feng_world(c0, S, seed=7)
    stages = t1_world(sp, cfg, genes, fst, group_change(c0, S, 0.2, 0.2))
    t47 = t47_of(sp, stages, fst, cfg)
    r = sp.split_reading(stages, fst, [], cfg, t47, n_boot=1)
    setup, _ = sp.build_setup(stages, fst, [], cfg, t47)
    M, w = setup["M"], setup["weights"]
    want = {}
    for v in ("c", "s"):
        vals = {g: 0.2 * float(setup["dirs"][g][M] @ setup["vecs"][v][g][M]) / float(
            setup["vecs"][v][g][M] @ setup["vecs"][v][g][M]) for g in t47["support"]}
        want[v] = sp.sr.combine(vals, w)[0]
    assert r["reading"]["A_c"]["value"] == pytest.approx(want["c"], abs=1e-6)
    assert r["reading"]["A_s"]["value"] == pytest.approx(want["s"], abs=1e-6)
    assert r["reading"]["A_M"]["value"] == pytest.approx(0.2, abs=1e-6)


def test_noise_in_the_halves_lowers_the_raw_share_along_s_and_the_correction_recovers_it(sp, cfg):
    """Cells move along the true direction by 0.3; Feng's E14.5 halves carry independent noise (sd 0.25 in the log
    ratio, equal per group). The raw share along s_g is the more attenuated (raw R 1.14-1.33 where the truth is 1); the
    corrected shares recover 0.3 within the tolerance measured on seeds 10-19 with 119 genes: largest error 0.022 on
    A_c_corr, 0.041 on A_s_corr and 0.13 on R_corr (the reliabilities rest on few genes); asserted at 0.03, 0.045 and
    0.15."""
    for seed in range(10, 20):
        c0, S = split_dirs(seed=seed, sd_c=0.3, sd_s=0.3)
        r, _ = _read(sp, cfg, c0, S, 0.3, 0.3, half_noise=0.25, seed=seed, n_boot=1)
        rd = r["reading"]
        assert rd["A_s"]["value"] < rd["A_c"]["value"] < 0.3  # raw: both attenuated, s_g the more
        assert rd["A_s"]["value"] < 0.3 - 0.05
        assert rd["R"]["value"] > 1.1  # the raw ratio is biased upwards by the unequal attenuation
        assert rd["A_c_corr"]["value"] == pytest.approx(0.3, abs=0.03)
        assert rd["A_s_corr"]["value"] == pytest.approx(0.3, abs=0.045)
        assert rd["R_corr"]["value"] == pytest.approx(1.0, abs=0.15)
        assert abs(rd["R_corr"]["value"] - 1.0) < abs(rd["R"]["value"] - 1.0)
        assert np.mean(list(r["reliability"]["c"].values())) > np.mean(list(r["reliability"]["s"].values()))


# ---------------------------------------------------------------------------------------------------------------------
# The statement and the void cases


def _rec(v, lo, hi, nf=0):
    return {"value": v, "p5": lo, "p95": hi, "nonfinite_draws": nf}


def test_each_statement_value_is_reached(sp):
    nb = 10
    less = {"R": _rec(0.5, 0.45, 0.55), "R_corr": _rec(0.5, 0.45, 0.55)}
    much = {"R": _rec(1.0, 0.95, 1.05), "R_corr": _rec(1.0, 0.95, 1.05)}
    bias0 = {"R": 0.0, "R_corr": 0.0}
    assert sp.statement(less, bias0, nb)["statement"] == sp.LESS
    assert sp.statement(much, bias0, nb)["statement"] == sp.AS_MUCH
    mid = {"R": _rec(0.85, 0.82, 0.88), "R_corr": _rec(0.85, 0.82, 0.88)}
    assert sp.statement(mid, bias0, nb)["statement"] == sp.NOT_SEPARATED
    # the bias widens both sides: a bias of 0.3 (either sign) on R_corr takes 'less' to 'not separated'
    st = sp.statement(less, {"R": 0.0, "R_corr": -0.3}, nb)
    assert st["statement"] == sp.NOT_SEPARATED and st["R_corr"] == sp.NOT_SEPARATED and st["R"] == sp.LESS
    assert st["turned_on"].startswith("R_corr")
    st = sp.statement(less, {"R": 0.3, "R_corr": 0.0}, nb)
    assert st["statement"] == sp.NOT_SEPARATED and st["turned_on"].startswith("R (raw)")
    st = sp.statement({"R": much["R"], "R_corr": less["R_corr"]}, bias0, nb)
    assert st["statement"] == sp.NOT_SEPARATED and "opposite" in st["turned_on"]
    # not computable where R_corr is undefined in a draw or at the point
    st = sp.statement({"R": less["R"], "R_corr": _rec(0.5, 0.45, 0.55, nf=3)}, bias0, nb)
    assert st["statement"] == "not computable (3 of 10 draws undefined)"
    st = sp.statement({"R": less["R"], "R_corr": _rec(float("nan"), 0.45, 0.55)}, bias0, nb)
    assert st["statement"] == "not computable (undefined at the point; 0 of 10 draws undefined)"
    # R raw undefined: the two differ, so 'not separated', turned on R
    st = sp.statement({"R": _rec(0.5, 0.45, 0.55, nf=1), "R_corr": less["R_corr"]}, bias0, nb)
    assert st["statement"] == sp.NOT_SEPARATED and st["turned_on"].startswith("R (raw): not computable")
    # the thresholds are strict
    assert sp.rule_statement(0.5, 0.8) == sp.NOT_SEPARATED and sp.rule_statement(0.9, 1.2) == sp.NOT_SEPARATED
    # an undefined bias of R_corr: not computable
    assert sp.statement(less, {"R": 0.0}, nb)["statement"].startswith("not computable (the self-check's bias of R_corr")
    # void beats everything
    assert sp.statement(less, bias0, nb, ["x"])["statement"] == "void (x)"


def test_the_void_cases(sp, cfg):
    c0, S = split_dirs(seed=8)
    genes, fst, _ = feng_world(c0, S, seed=8)
    stages = t1_world(sp, cfg, genes, fst, group_change(c0, S, 0.2, 0.2))
    t47 = t47_of(sp, stages, fst, cfg)
    setup, extra = sp.build_setup(stages, fst, [], cfg, t47)
    point = sp.measure_split(extra["tabs"], {s: sp.sr.block_sums(t) for s, t in extra["tabs"].items()}, setup)
    good_sc = {"measured": {"R": 1.0}}
    assert sp.void_reasons(point, setup, t47["A"], good_sc) == []
    assert sp.void_reasons(point, setup, t47["A"] + 0.019, good_sc) == []
    (why,) = sp.void_reasons(point, setup, t47["A"] + 0.021, good_sc)
    assert why.startswith("A_M") and "within 0.02" in why
    (why,) = sp.void_reasons(point, setup, t47["A"], {"measured": {"R": float("nan")}})
    assert why == "the self-check's R is undefined"
    (why,) = sp.void_reasons(point, setup, t47["A"], None)
    assert why == "the self-check did not run"
    # a group of T47's A set without a value (its library-3 cells gone)
    e105 = stages["E10.5"]
    drop = (e105.obs["celltype"] == "Epicardium") & np.array([n.endswith("_3") for n in e105.obs_names])
    stages2 = {**stages, "E10.5": e105[~drop.to_numpy()].copy()}
    setup2, extra2 = sp.build_setup(stages2, fst, [], cfg, t47)
    point2 = sp.measure_split(extra2["tabs"], {s: sp.sr.block_sums(t) for s, t in extra2["tabs"].items()}, setup2)
    reasons = sp.void_reasons(point2, setup2, t47["A"], good_sc)
    assert any("['Epicardium']" in x for x in reasons)


def test_the_held_genes_leave_M_and_A_M_stays_close_to_t47(sp, cfg):
    c0, S = split_dirs(seed=9)
    genes, fst, _ = feng_world(c0, S, seed=9)
    stages = t1_world(sp, cfg, genes, fst, group_change(c0, S, 0.2, 0.2))
    t47 = t47_of(sp, stages, fst, cfg, held=[genes[5]])
    setup, _ = sp.build_setup(stages, fst, [genes[5]], cfg, t47)
    assert not setup["M"][5] and not setup["M"][0] and setup["M"].sum() == G - 2
    r = sp.split_reading(stages, fst, [genes[5]], cfg, t47, n_boot=1)
    assert abs(r["reading"]["A_M"]["value"] - t47["A"]) <= 0.02


def test_the_calibration_reads_t47s_own_step_along_the_split(sp, cfg):
    c0, S = split_dirs(seed=11)
    genes, fst, _ = feng_world(c0, S, seed=11)
    stages = t1_world(sp, cfg, genes, fst, group_change(c0, S, 0.2, 0.2), counts=True, seed=11)
    t47 = t47_of(sp, stages, fst, cfg)
    r = sp.split_reading(stages, fst, [], cfg, t47, n_boot=1, calibrate=True)
    cal = r["calibration"]
    t47_cal = sp.sr.reading(stages, fst, [], {}, cfg, n_boot=1, calibrate=True)["calibration"]
    assert cal["combined_weighted"]["d_T47"] == pytest.approx(t47_cal["combined_weighted"], rel=1e-12)  # T47's arm
    assert cal["stream"] == [47, 1] and cal["fraction"] == 0.2
    for v in ("M", "c", "s"):
        assert cal["combined_weighted"][v] == pytest.approx(0.2, abs=0.02)
    assert cal["R"] == pytest.approx(1.0, abs=0.1)


def test_the_bootstrap_is_seeded_on_its_own_stream(sp, cfg):
    assert sp.BOOT_SEED == (51, 0) and sp.N_BOOT == 200
    c0, S = split_dirs(seed=12)
    genes, fst, _ = feng_world(c0, S, seed=12)
    stages = t1_world(sp, cfg, genes, fst, group_change(c0, S, 0.2, 0.2), counts=True, seed=12)
    t47 = t47_of(sp, stages, fst, cfg)
    r1 = sp.split_reading(stages, fst, [], cfg, t47, n_boot=3)
    r2 = sp.split_reading(stages, fst, [], cfg, t47, n_boot=3)
    r3 = sp.split_reading(stages, fst, [], cfg, t47, n_boot=3, seed=(51, 1))
    assert json.dumps(r1["values"], sort_keys=True) == json.dumps(r2["values"], sort_keys=True)
    assert r1["values"]["reading/R"]["p5"] != r3["values"]["reading/R"]["p5"]


def test_the_t47_record_is_read_from_its_committed_layout(sp, tmp_path):
    p = tmp_path / "t47.json"
    p.write_text(json.dumps({"read": {"reading": {"A": {"value": 0.1}, "L": {"value": 0.9}, "L_corr": {"value": 1.1}},
                                      "support": {"A_K_support": ["CM"]}}}))
    assert sp.t47_record(p) == {"A": 0.1, "L": 0.9, "L_corr": 1.1, "support": ["CM"]}
    assert sp.T47_JSON.parts[-2:] == ("t47_share_read", "t47.json")


def test_the_limits_are_in_the_docstring_and_the_output_holds_no_cell_count(sp):
    starts = ("A lower share along the common part can mean", "The challenge's own sample differences",
              "The specific part is a difference of noisy directions", "L_eff assumes the specific part's scale is L",
              "No statement says how far the hidden test's days move")
    assert all(any(n.startswith(s) for n in sp.LEDGER_NOTES) for s in starts)
    src = SCRIPT.read_text()
    assert '"weights"' not in src.split("res = {", 1)[1].split("}", 1)[0]  # the record carries no weights


# ---------------------------------------------------------------------------------------------------------------------
# Refusals


def test_refused_without_the_row(sp, monkeypatch):
    monkeypatch.setattr(sp.n12, "tries_rows", lambda *a: {"T47"})
    with pytest.raises(SystemExit, match="no row"):
        sp.run("x")


def test_refused_on_an_existing_output(sp, monkeypatch):
    monkeypatch.setattr(sp.n12, "tries_rows", lambda *a: {"T51"})
    out = sp.RESULTS / "x"
    out.mkdir(parents=True)
    (out / "t51.txt").write_text("")
    with pytest.raises(SystemExit, match="already holds"):
        sp.run("x")


def test_refused_on_a_dirty_tree(sp, monkeypatch):
    monkeypatch.setattr(sp.n12, "tries_rows", lambda *a: {"T51"})
    monkeypatch.setattr(sp.gv, "tree_state", lambda: {"commit": "abc", "clean": False})
    monkeypatch.setattr(sp.feng_dirs, "load_config", lambda *a: pytest.fail("read past the refusal"))
    with pytest.raises(SystemExit, match="not clean"):
        sp.run("x")
