"""Tries K4 and K4c and diagnostic N19 (experiments/harness/TRIES.md), on synthetic data only: the empirical-Bayes
shrinkage of the T3 source factor (experiments/t3_decomp/eb_factor.py), the harness candidates t3_factor_raw,
t3_factor_eb and t3_factor_eb_ctrl, N19's pre-gate and K4's keep-or-stop rule (no competition file, no score of a fold).
The raw factor is K2a's estimator exactly (`t3_source_frame`, pinned against `t3_gata4_factor`) and the raw arm K2a's
build on the Gata4 -> Mab21l2 fold; the jackknife is the definition pinned in K4's row; the control permutes the weights
within the factor's own deciles (K2b's logic)."""

from __future__ import annotations

import json
import shutil
import sys
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from vec import harness, provenance
from vec.boards import GROUPS, Board
from vec.evaluate import Stage
from vec.prep import load_stage, log1p_cp10k

sys.path.insert(0, str(harness.T3_DIR))
import eb_factor as eb  # noqa: E402
import k4_verdict  # noqa: E402
import n19_eb_pregate as n19  # noqa: E402

G = 40
GENES = [f"G{j}" for j in range(G)]
# WT E8.75-like, in the gata4 preset's combined labels (celltype, cm_celltype, share): strict non-Mesp1 labels (NCC,
# SE, Foregut, Neural Tube) and Mesp1-derived ones
E875 = [("V-CM", "OFT/RV-CM", 0.11), ("V-CM", "vCM1", 0.11), ("aSHF", "Unknown", 0.08), ("pSHF", "Unknown", 0.07),
        ("aPHM", "Unknown", 0.08), ("Peri", "Unknown", 0.07), ("Endo", "Unknown", 0.06), ("NCC", "Unknown", 0.12),
        ("SE", "Unknown", 0.10), ("Foregut", "Unknown", 0.10), ("Neural Tube", "Unknown", 0.10)]
# WT E9.5-like (the mab21l2 preset: celltype labels)
E95 = [("V-CM", None, 0.25), ("A-CM", None, 0.10), ("Peri", None, 0.10), ("SE", None, 0.10), ("Foregut", None, 0.10),
       ("NCC", None, 0.15), ("ST", None, 0.10), ("Endo", None, 0.10)]
G2M, M2G, REP = "T3/Gata4->Mab21l2", "T3/Mab21l2->Gata4rep2", "T3/Gata4rep1->rep2"
RELEASED = {"wt_e875": "T2/E8.75.h5ad", "wt_e95": "T2/E9.5.h5ad", "gata4_1": "T3/E8.75_gata4_ko_1.h5ad",
            "gata4_2": "T3/E8.75_gata4_ko_2.h5ad", "mab_ko": "T3/E9.5_mab21l2_ko.h5ad"}


def _label(cell) -> str:
    return cell[1] if cell[1] not in (None, "Unknown") else cell[0]


def _cells(cells, prof, n, rng, depth, factor=None, condition=None, name="c") -> ad.AnnData:
    """A MERFISH-like AnnData: labels drawn from the shares, Poisson counts around per-label profiles (times a per-gene
    factor), log1p CP10k X, a counts layer and section-like positions on 20 z-levels."""
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
    """Synthetic files in the T3 folds' roles, loaded as the harness loads them (load_stage, duplicates dropped).
    Gata4-like replicates: shallower than their WT, a shared per-gene factor, fewer V-CM cells, littermate WT cells in
    replicate 1. A Mab21l2-like knockout of the WT E9.5-like file with its own per-gene factor."""
    d = tmp_path_factory.mktemp("k4")
    rng = np.random.default_rng(0)
    prof = {k: rng.lognormal(0, 1.0, G) for k in sorted({_label(c) for c in E875 + E95})}
    gata4 = [(c, m, s * (0.5 if c == "V-CM" else 1.0)) for c, m, s in E875]
    gfac = rng.lognormal(0, 0.5, G)
    files = {"wt_e875": _cells(E875, prof, 900, rng, 150, name="e875"),
             "gata4_1": _cells(gata4, prof, 600, rng, 90, gfac, name="rep1",
                               condition=np.where(rng.random(600) < 0.2, "WT", "KO")),
             "gata4_2": _cells(gata4, prof, 600, rng, 95, gfac, name="rep2", condition=np.full(600, "KO")),
             "wt_e95": _cells(E95, prof, 1200, rng, 75, name="e95"),
             "mab_ko": _cells([(c, m, s * (2.0 if c == "ST" else 1.0)) for c, m, s in E95], prof, 800, rng, 60,
                              rng.lognormal(0, 0.4, G), name="mab")}
    paths = {role: d / f"{role}.h5ad" for role in files}
    for role, a in files.items():
        a.write_h5ad(paths[role])
    loaded = {role: load_stage(p, genes=GENES, drop_duplicates=True) for role, p in paths.items()}
    target = _cells([(c, m, s * (2.0 if c == "ST" else 1.0)) for c, m, s in E95], prof, 2000, rng, 90,
                    rng.lognormal(0, 0.4, G), name="ko95")
    return SimpleNamespace(paths=paths, loaded=loaded,
                           target=Stage(X=np.asarray(target.X.toarray(), np.float32),
                                        coords=np.asarray(target.obsm["spatial_3D"], np.float32)),
                           inputs={G2M: {"wt": loaded["wt_e95"], "ko_train": loaded["gata4_1"],
                                         "ko_train_2": loaded["gata4_2"], "wt_train": loaded["wt_e875"]},
                                   M2G: {"wt": loaded["wt_e875"], "ko_train": loaded["mab_ko"],
                                         "wt_train": loaded["wt_e95"]},
                                   REP: {"wt": loaded["wt_e875"], "ko_train": loaded["gata4_1"],
                                         "wt_train": loaded["wt_e875"]}})


@pytest.fixture(scope="module")
def t3():
    bs = harness._t3_head()
    return SimpleNamespace(bs=bs, t3d=bs.t3d, rules=bs.load_rules())


@pytest.fixture(scope="module")
def frames(world, t3):
    return {G2M: harness.t3_source_frame(world.inputs[G2M]["wt_train"],
                                         [world.inputs[G2M]["ko_train"], world.inputs[G2M]["ko_train_2"]], t3.rules,
                                         "gata4"),
            M2G: harness.t3_source_frame(world.inputs[M2G]["wt_train"], [world.inputs[M2G]["ko_train"]], t3.rules,
                                         "mab21l2")}


def _fold(key) -> harness.PseudoTask:
    return {t.key: t for t in harness.load_tasks()["T3:bcatenin"]}[key]


def _board(key=G2M, max_cells: int = 900) -> Board:
    return Board(key=f"pseudo:{key}", task="T3", genes=list(GENES), needs_coords=True, split="pseudo",
                 target_ko="mab21l2", min_cells=100, max_cells=max_cells, truth_cells=120, ref_cells=150)


# ----------------------------------------------------------------------------------------------- the jackknife, by hand
def _toy(seed=0, eps=1e-5):
    """Cells with labels (one label outside the mask, one label living in a single section, one cell unlabelled) and
    section blocks over two 'replicates'."""
    rng = np.random.default_rng(seed)
    n, g, L = 300, 12, 5
    C = rng.poisson(rng.lognormal(0, 1, (1, g)), (n, g)).astype(np.int32)
    C[:, 3] *= rng.integers(0, 2, n)  # a sparser gene
    lab = rng.integers(0, 4, n)
    block = rng.integers(0, 6, n)
    lab[:7] = 4  # label 4 lives only in section 5 (a deletion empties it; it stays in the mask)
    block[:7] = 5
    lab[7] = -1  # an unlabelled cell
    mask = np.array([True, True, False, True, True])
    w = np.array([0.3, 0.2, 0.2, 0.2, 0.1])
    mu_wt = rng.dirichlet(np.ones(g), L)
    return SimpleNamespace(C=C, lab=lab, block=block, mu_wt=mu_wt, w=w, mask=mask, eps=eps, L=L, G=g)


def _manual_delete_one(t):
    """The definition in K4's row, spelled out: numerator recomputed without each section; denominator, mask and label
    weights fixed; a label a deletion empties keeps its full-sample mean; K = the sections holding estimator cells."""
    props = t.C / np.maximum(t.C.sum(1, keepdims=True), 1)
    keep = (t.lab >= 0) & t.mask[np.maximum(t.lab, 0)]
    keep[t.lab < 0] = False
    ww = t.w * t.mask
    den = (ww[:, None] * t.mu_wt).sum(0)

    def num(sel):
        mus = []
        for label in range(t.L):
            cells = keep & (t.lab == label) & sel
            full = keep & (t.lab == label)
            mus.append(props[cells].mean(0) if cells.any() else (props[full].mean(0) if full.any() else np.zeros(t.G)))
        return (ww[:, None] * np.stack(mus)).sum(0)

    y_full = np.log((num(np.ones(t.lab.size, bool)) + t.eps) / (den + t.eps))
    sections = np.unique(t.block[keep])
    ys = np.stack([np.log((num(t.block != s) + t.eps) / (den + t.eps)) for s in sections])
    K = sections.size
    v = (K - 1) / K * ((ys - ys.mean(0)) ** 2).sum(0)
    return y_full, v, den, K


def test_jackknife_is_the_rows_definition():
    t = _toy()
    jack = eb.jackknife(t.C, t.lab, t.block, t.mu_wt, t.w, t.mask, t.eps)
    y, v, den, K = _manual_delete_one(t)
    assert jack.n_blocks == K == 6
    assert np.array_equal(jack.den, den)  # the denominator is gene_factor's, untouched
    assert np.allclose(jack.y, y, rtol=0, atol=1e-12) and np.allclose(jack.v, v, rtol=0, atol=1e-12)
    # sections holding none of the estimator's cells are not blocks; unlabelled and unmasked cells never count
    b = t.block.copy()
    b[(t.lab == 2) | (t.lab < 0)] = 99
    assert eb.jackknife(t.C, t.lab, b, t.mu_wt, t.w, t.mask, t.eps).n_blocks == 6
    with pytest.raises(ValueError, match="at least two sections"):
        eb.jackknife(t.C, t.lab, np.zeros_like(t.block), t.mu_wt, t.w, t.mask, t.eps)
    with pytest.raises(ValueError, match="no WT mass"):
        eb.jackknife(t.C, t.lab, t.block, t.mu_wt, t.w, np.zeros(t.L, bool), t.eps)


def test_the_denominator_is_never_jackknifed():
    # the numerator side only: the WT means enter y through a fixed denominator, so with eps negligible, multiplying
    # them by a constant moves y by the log of that constant and leaves the jackknife variance unchanged
    t = _toy(3, eps=1e-12)
    a = eb.jackknife(t.C, t.lab, t.block, t.mu_wt, t.w, t.mask, t.eps)
    b = eb.jackknife(t.C, t.lab, t.block, t.mu_wt * 3.0, t.w, t.mask, t.eps)
    assert np.allclose(a.y - b.y, np.log(3.0), atol=1e-6) and np.allclose(a.v, b.v, atol=1e-6)


def test_jackknife_y_is_gene_factors_ratio_and_finalize_is_its_rescale_and_clip(t3):
    t = _toy(1)
    jack = eb.jackknife(t.C, t.lab, t.block, t.mu_wt, t.w, t.mask, t.eps)
    mu_ko = t3.t3d.props_type_means(t.C, t.lab, t.L)
    phi, den = t3.t3d.gene_factor(mu_ko, t.mu_wt, t.w, t.mask, t.eps, (0.25, 4.0))
    assert np.array_equal(jack.den, den)
    ww = (t.w * t.mask)[:, None]
    assert np.allclose(np.exp(jack.y), ((ww * mu_ko).sum(0) + t.eps) / (den + t.eps))
    assert np.allclose(eb.finalize(jack.y, den), phi)  # with B = 1 and c = 1 the EB path is gene_factor's phi
    assert np.array_equal(np.clip(eb.finalize(np.log(1e3 * np.ones(t.G)), den), 0.25, 4.0),
                          eb.finalize(np.log(1e3 * np.ones(t.G)), den))


def test_block_ids_are_replicate_by_level():
    z1, z2 = np.array([0.0, 0.0, 10.0, 20.0, 20.0]), np.array([5.0, 15.0, 15.0])
    block, n = eb.block_ids([z1, z2])
    assert n == 5 and list(block) == [0, 0, 1, 2, 2, 3, 4, 4]  # levels of the second file start after the first's
    assert eb.block_ids([z1])[1] == 3
    # float noise within a level stays one level (snap_z)
    assert eb.block_ids([np.array([0.0, 1e-9, 10.0])])[1] == 2


# -------------------------------------------------------------------------------------------- the prior and the weights
def test_eb_weights_are_dersimonian_laird_with_prior_mean_zero():
    # two deciles of the basis: in the second, tau2 > 0; in the first, tau2 = 0
    den = np.array([0.1, 0.2, 0.3, 0.4, 1.0, 2.0, 3.0, 4.0])
    y = np.array([0.0, 0.1, -0.1, 0.05, 2.0, -1.5, 1.0, 0.0])
    v = np.array([0.5, 0.5, 0.5, 0.5, 0.1, 0.2, 0.05, 0.1])
    B, tau2 = eb.eb_weights(y, v, den, bins=2)
    w = 1.0 / v[4:]
    t2 = ((w * y[4:] ** 2).sum() - 4.0) / w.sum()
    assert tau2[1] == pytest.approx(t2) and t2 > 0 and tau2[0] == 0.0  # low spread, tau2 floored at zero
    assert np.array_equal(B[:4], np.zeros(4))  # tau2 = 0: full shrinkage to f = 1
    assert np.allclose(B[4:], t2 / (t2 + v[4:]))
    # prior mean zero, not a fitted mean: a shifted cloud with little spread still gives tau2 > 0
    shifted, t3_ = eb.eb_weights(np.full(10, 2.0), np.full(10, 0.01), np.arange(10.0), bins=1)
    assert t3_[0] == pytest.approx(((1.0 / 0.01 * 4.0 * 10) - 10) / (10 / 0.01)) and (shifted > 0.99).all()


def test_eb_weights_guards_a_gene_with_no_measured_noise_and_a_decile_with_none():
    den = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    y = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    v = np.array([0.0, 0.1, 0.3, 0.0, 0.0, 0.0])  # second decile: no noise measured anywhere
    B, tau2 = eb.eb_weights(y, v, den, bins=2)
    v0 = np.array([0.2, 0.1, 0.3])  # the zero takes the decile's median positive variance, 0.2
    w = 1.0 / v0
    t2 = max(0.0, (w.sum() - 3.0) / w.sum())
    assert tau2[0] == pytest.approx(t2) and np.allclose(B[:3], t2 / (t2 + v0))
    assert np.array_equal(B[3:], np.ones(3)) and tau2[1] == 0.0  # none measured: no shrinkage


def test_robust_sd_and_amplitude_constant():
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 100.0])
    assert eb.robust_sd(x) == pytest.approx(1.4826 * np.median(np.abs(x - np.median(x))))
    y = np.linspace(-1.0, 1.0, 40)
    den = np.ones(40)
    jack = eb.Jackknife(y=y, v=np.full(40, 0.05), den=den, n_blocks=30)  # a small variance: tau2 > 0, B positive
    phi_raw = np.exp(y * 0.5)
    fit = eb.eb_factors(phi_raw, jack, bins=1)
    m = fit.weights * y
    assert (fit.weights > 0).all()
    assert fit.constant == pytest.approx(eb.robust_sd(np.log(phi_raw)) / eb.robust_sd(m))
    # the EB arm matches the raw factor's robust SD (the clip does not bind); the final rescale preserves the WT total
    assert eb.robust_sd(np.log(fit.eb)) == pytest.approx(eb.robust_sd(np.log(phi_raw)), rel=1e-6)
    assert fit.raw is not phi_raw and np.array_equal(fit.raw, phi_raw)


def test_a_fully_shrunk_factor_is_flat_and_needs_no_constant():
    y = np.array([0.5, -0.5, 0.2, 0.1])
    jack = eb.Jackknife(y=y, v=np.full(4, 100.0), den=np.ones(4), n_blocks=10)
    fit = eb.eb_factors(np.exp(y), jack, bins=1)
    assert np.array_equal(fit.weights, np.zeros(4)) and fit.constant == 1.0
    assert np.array_equal(fit.eb, np.ones(4)) and np.array_equal(fit.ctrl, np.ones(4))


def test_a_section_unstable_gene_is_shrunk_harder_than_a_stable_one():
    # gene 0: the same effect in every section; gene 1: the whole effect comes from one section of the ten
    rng = np.random.default_rng(5)
    n, sections = 4000, 10
    block = rng.integers(0, sections, n)
    C = rng.poisson(5.0, (n, 8)).astype(np.int32)
    C[:, 0] = rng.poisson(10.0, n)              # twice the WT level everywhere
    C[:, 1] = rng.poisson(2.5, n) + 15 * rng.poisson(2.0, n) * (block == 0)   # strong only in one section
    lab, L = np.zeros(n, int), 2
    mu_wt = np.ones((L, 8)) / 8
    jack = eb.jackknife(C, lab, block, mu_wt, np.array([1.0, 0.0]), np.array([True, False]), 1e-5)
    assert jack.v[1] > 10 * jack.v[0] > 0
    B, _ = eb.eb_weights(jack.y, jack.v, jack.den, bins=1)
    assert B[1] < B[0]


def test_the_shrinkage_recovers_a_known_factor_better_on_section_unstable_genes():
    # a known per-gene log factor; half the genes carry a large section-to-section random effect, independent of it
    rng = np.random.default_rng(7)
    genes, sections, cells = 400, 30, 8000
    p = rng.lognormal(0, 1, genes)
    p /= p.sum()
    theta = rng.normal(0, 0.35, genes)
    unstable = rng.random(genes) < 0.5
    noise = rng.normal(0, 1, (sections, genes)) * np.where(unstable, 0.9, 0.05)
    block = rng.integers(0, sections, cells)
    rate = p * np.exp(theta + noise[block])
    C = rng.poisson(200 * rate / rate.sum(1, keepdims=True))
    jack = eb.jackknife(C, np.zeros(cells, int), block, p[None, :], np.array([1.0]), np.array([True]), 1e-5)
    B, _ = eb.eb_weights(jack.y, jack.v, jack.den)
    assert np.median(B[~unstable]) > 0.95 and np.median(B[~unstable]) > np.median(B[unstable]) + 0.1
    err_raw = ((jack.y - theta)[unstable] ** 2).mean()
    err_eb = (((B * jack.y) - theta)[unstable] ** 2).mean()
    assert err_eb < 0.8 * err_raw


# ------------------------------------------------------------------------------------------------------- the control
def test_the_control_is_k2bs_permutation_of_the_weights_within_the_deciles():
    assert (harness.K4_CTRL_SEED, harness.K4_BINS) == (20261109, 10)  # fixed in K4c's row
    rng = np.random.default_rng(3)
    for n in (500, 103):
        x, basis = rng.random(n), rng.gamma(1.0, 1.0, n)
        out = eb.decile_permutation(x, basis, harness.K4_CTRL_SEED)
        assert np.array_equal(out, harness.t3_decile_permutation(x, basis, seed=harness.K4_CTRL_SEED, bins=10))
        for b in np.array_split(np.argsort(basis, kind="stable"), 10):
            assert np.array_equal(np.sort(out[b]), np.sort(x[b]))  # each decile keeps its values
        assert not np.array_equal(out, eb.decile_permutation(x, basis, 20261110))


def test_eb_factors_are_deterministic_and_the_control_shares_c(frames):
    frame = frames[G2M]
    jack = eb.jackknife(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask, frame.eps)
    fit = eb.eb_factors(frame.phi, jack, clip=frame.clip, seed=harness.K4_CTRL_SEED, bins=harness.K4_BINS)
    again = harness.t3_k4_fit(frame)
    for name in ("raw", "eb", "ctrl", "weights", "weights_ctrl", "tau2"):
        assert np.array_equal(getattr(fit, name), getattr(again, name)), name
    assert fit.constant == again.constant
    # the raw arm is K2a's factor as is; EB and control are finite, inside the clip, and differ from it
    assert np.array_equal(fit.raw, frame.phi)
    for phi in (fit.eb, fit.ctrl):
        assert np.isfinite(phi).all() and (phi >= 0.25).all() and (phi <= 4.0).all()
        assert not np.array_equal(phi, frame.phi)
    # the same constant on the control: its log factor is c x the permuted weights x y, finalized as the EB arm
    c = eb.robust_sd(np.log(frame.phi)) / eb.robust_sd(fit.weights * jack.y)
    assert fit.constant == pytest.approx(c)
    assert np.array_equal(fit.ctrl, eb.finalize(fit.constant * (fit.weights_ctrl * jack.y), jack.den, frame.clip))
    assert np.array_equal(fit.eb, eb.finalize(fit.constant * (fit.weights * jack.y), jack.den, frame.clip))
    assert np.array_equal(np.sort(fit.weights), np.sort(fit.weights_ctrl))


# ------------------------------------------------------------------------------------------------- the source frame
def test_raw_is_k2as_factor_on_the_gata4_source(world, frames, t3):
    inp = world.inputs[G2M]
    k2a = harness.t3_gata4_factor(inp["wt_train"], [inp["ko_train"], inp["ko_train_2"]], t3.rules)
    assert np.array_equal(frames[G2M].phi, k2a)
    # the frame's jackknife denominator is gene_factor's, and the frame holds the pooled kept cells only (littermates
    # off)
    frame = frames[G2M]
    jack = eb.jackknife(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask, frame.eps)
    assert np.array_equal(jack.den, frame.den)
    a1, a2 = inp["ko_train"], inp["ko_train_2"]
    lit = np.asarray(a1.obs["condition"].astype(str) == "WT").sum()
    assert frame.C.shape[0] == a1.n_obs + a2.n_obs - lit
    # replicate-by-level blocks: the two files' sections are distinct
    n_lev1 = eb.block_ids([np.asarray(a1.obsm["spatial_3D"])[:, 2]])[1]
    assert frame.n_levels > n_lev1


def test_the_mab21l2_source_uses_its_own_wt_vocabulary_and_preset(world, frames, t3):
    inp = world.inputs[M2G]
    frame = frames[M2G]
    t3d = t3.t3d
    # the same steps, spelled out for the mab21l2 preset: celltype labels of WT E9.5
    w = harness._shallow(inp["wt_train"])
    t3d.add_label_column(w, t3.rules["presets"]["mab21l2"]["label_mode"])
    V, wlab = harness._t3_vocab(t3d, w, t3.rules)
    assert "vCM1" not in V.labels and "V-CM" in V.labels  # celltype labels, not the combined ones
    Cw = t3d.dense_counts(w)
    harm = t3d.Harmoniser(w, "t3_label", t3d.label_groups(w, "t3_label"), t3.rules["harmonisation"], 1)
    a = harness._shallow(inp["ko_train"])
    calls, _ = harm.calls(a, seed=t3.rules["harmonisation"]["seed"])
    C, lab = t3d.dense_counts(a), t3d.calls_to_index(V, calls)
    fc = t3.rules["components"]["factor"]
    ok = ((np.bincount(lab[lab >= 0], minlength=V.L) >= fc["min_cells"])
          & (np.bincount(wlab, minlength=V.L) >= fc["min_cells"]))
    phi, _ = t3d.gene_factor(t3d.props_type_means(C, lab, V.L), t3d.props_type_means(Cw, wlab, V.L), V.p, ok & V.nm,
                             fc["eps"], tuple(fc["clip"]))
    assert np.array_equal(frame.phi, phi) and (frame.block.min(), frame.n_levels) == (0, len(np.unique(frame.block)))


# ----------------------------------------------------------------------------------------------------- the candidates
def test_k4_arms_apply_on_the_counted_folds_and_the_reported_replicate_fold():
    tasks = harness.load_tasks()
    assert sorted(k for k, f in harness.K4_FOLDS.items() if f["counted"]) == sorted([G2M, M2G])
    assert [k for k, f in harness.K4_FOLDS.items() if not f["counted"]] == [REP]
    for name in ("t3_factor_raw", "t3_factor_eb", "t3_factor_eb_ctrl"):
        c = harness.CANDIDATES[name]
        assert sorted(t.key for ts in tasks.values() for t in ts if c.applies(t)) == sorted([G2M, M2G, REP])


def _build_context(world, key) -> harness.BuildContext:
    return harness.BuildContext(task=harness.TaskView.of(_fold(key)), board=harness._candidate_board(_board(key)),
                                inputs=world.inputs[key], seed=harness.BUILD_SEED_OFFSET + 2)


@pytest.mark.parametrize("key", [G2M, M2G, REP])
def test_k4_arms_build_the_same_cells_with_their_own_factors(world, t3, key, monkeypatch):
    calls, real = [], t3.t3d.apply_gene_factor

    def spy(C, phi, P, rng):
        calls.append((C.copy(), np.array(phi)))
        return real(C, phi, P, rng)

    monkeypatch.setattr(t3.t3d, "apply_gene_factor", spy)
    harness._K4_PREP.clear()
    ctx = _build_context(world, key)
    with harness.guarded_reads() as g:
        built = {name: harness.CANDIDATES[f"t3_factor_{name}"].build(ctx) for name in ("raw", "eb", "eb_ctrl")}
    assert not g["refused"] and "experiments/t3_decomp/rules.yaml" in g["log"]
    prep = harness._K4_PREP["fold"]
    (Cr, pr), (Ce, pe), (Cc, pc) = calls
    assert np.array_equal(Cr, Ce) and np.array_equal(Cr, Cc)  # the same drawn cells in all three arms
    assert np.array_equal(pr, prep["raw"]) and np.array_equal(pe, prep["eb"]) and np.array_equal(pc, prep["ctrl"])
    assert not np.array_equal(built["raw"].X, built["eb"].X)
    assert all(s.X.shape == (900, G) for s in built.values())
    # the reference is the fold's own WT in the fold's reference preset's labels (combined labels for WT E8.75, the
    # reference of the Mab21l2 -> Gata4rep2 fold and of the replicate fold)
    assert ("vCM1" in prep["ref"].V.labels) == (key != G2M)


def test_the_prepared_factors_follow_the_fold_across_folds(world, frames):
    # the folds hold different numbers of input roles, and one harness process runs them in turn
    harness._K4_PREP.clear()
    for key in (G2M, M2G, REP, G2M):
        with harness.guarded_reads():
            harness.CANDIDATES["t3_factor_raw"].build(_build_context(world, key))
        hit = harness._K4_PREP["fold"]
        assert hit["key"] == key and hit["inputs"][0] is world.inputs[key]["wt_train"]
    assert np.array_equal(harness._K4_PREP["fold"]["raw"], frames[G2M].phi)
    # the second G2M call had to rebuild it after REP; a fold handed new objects rebuilds too
    before = harness._K4_PREP["fold"]["inputs"]
    harness.CANDIDATES["t3_factor_raw"].build(_build_context(world, G2M))
    assert harness._K4_PREP["fold"]["inputs"] is before  # same objects: the cached factors are reused


def test_the_raw_arm_is_k2as_build_on_the_gata4_to_mab21l2_fold(world):
    harness._K2_PREP.clear()
    harness._K4_PREP.clear()
    ctx = _build_context(world, G2M)
    with harness.guarded_reads() as g:
        k2a = harness.CANDIDATES["t3_factor_xko"].build(ctx)
        raw = harness.CANDIDATES["t3_factor_raw"].build(ctx)
    assert not g["refused"]
    assert np.array_equal(raw.X, k2a.X) and np.array_equal(raw.coords, k2a.coords)


def _k4_result(key, gains, seeds=harness.SEEDS):
    cands = ("floor", "copy", "t3_factor_raw", "t3_factor_eb", "t3_factor_eb_ctrl")
    rows = [{"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
             "scores": {"floor": 50.0, "copy": 50.5, "t3_factor_raw": 49.0, "t3_factor_eb_ctrl": 49.0,
                        "t3_factor_eb": 49.0 + g}} for s, g in zip(seeds, gains, strict=True)]
    prov = {"commit": "abc", "tree_clean": True, "config_sha256": harness._sha256(harness.CONFIG),
            "files_sha256": {"T2/E9.5.h5ad": "f"}}
    return {"task": {"key": key, "test_board": "T3:bcatenin"}, "board": {"key": f"pseudo:{key}"},
            "candidates": list(cands), "rows": rows, "protocol": "bake", "build": "main", "anchors": "computed",
            "frac": 0.1, "seeds": list(seeds), "provenance": prov}


@pytest.mark.parametrize("base", ["t3_factor_raw", "t3_factor_eb_ctrl"])
def test_decide_needs_both_counted_folds_against_each_base(base):
    configured = harness.load_tasks()["T3:bcatenin"]
    both = [_k4_result(G2M, [0.5] * 10), _k4_result(M2G, [0.5] * 10)]
    d = harness.decide(both, configured, base, "t3_factor_eb")
    assert d["keep"] and sorted(d["counted"]) == sorted([G2M, M2G]) and d["need"] == 2
    assert d["excluded"] == {REP: "replicate fold"}
    one_short = [_k4_result(G2M, [0.5] * 10), _k4_result(M2G, [0.1] * 10)]
    assert not harness.decide(one_short, configured, base, "t3_factor_eb")["keep"]
    assert not harness.decide([_k4_result(G2M, [0.5] * 10)], configured, base, "t3_factor_eb")["keep"]  # incomplete
    # the replicate fold's result is reported and never counted: a bad replicate fold changes nothing
    with_rep = [*both, _k4_result(REP, [-5.0] * 10)]
    d2 = harness.decide(with_rep, configured, base, "t3_factor_eb")
    assert d2["keep"] and d2["counted"] == d["counted"] and REP not in d2["tasks"]


def test_score_seeds_runs_the_k4_arms_and_the_fit_checks_on_them(world):
    key = G2M
    ref = world.inputs[key]["wt"]
    target, refstage = world.target, Stage(X=np.asarray(ref.X.toarray(), np.float32),
                                           coords=np.asarray(ref.obsm["spatial_3D"], np.float32))
    cands = [harness.CANDIDATES[n] for n in ("floor", "t3_factor_raw", "t3_factor_eb", "t3_factor_eb_ctrl")]
    rows = harness.score_seeds(_fold(key), _board(key), target, refstage, world.inputs[key], cands, seeds=(0, 1))
    for row in rows:
        assert row["n_cells"]["t3_factor_eb"] == 900
        assert all(all(d["matches_panel"].values()) for d in row["t3"].values())
    for base in ("t3_factor_raw", "t3_factor_eb_ctrl"):
        fc = harness.t3_fit_checks({"rows": rows}, base, "t3_factor_eb")
        assert fc["seeds"] == [0, 1] and set(fc["group_margins"]) == {g for g, _, _ in GROUPS["T3"]}


# ------------------------------------------------------------------------------------------------------ N19's pre-gate
def test_predicted_change_and_ranks():
    den = np.array([0.5, 0.5])
    assert np.allclose(n19.predicted_change([2.0, 1.0], den), [np.log1p(1e4) - np.log1p(5e3), 0.0], atol=1e-12)
    assert list(n19.ranks([1.0, 3.0, 3.0, 2.0])) == [3, 0, 1, 2]  # ties in panel order


def test_source_gate_flags_a_no_op_only_when_weights_and_ranks_both_say_so():
    rng = np.random.default_rng(0)
    den, raw = rng.gamma(1.0, 1.0, 300), np.exp(rng.normal(0, 0.5, 300))
    assert n19.source_gate(raw, raw, np.full(300, 1.0), den) == {
        "top": 100, "median_weight": 1.0, "n_moved": 0, "noop": True}
    # weights low: not a no-op even when the ranks do not move
    assert not n19.source_gate(raw, raw, np.full(300, 0.5), den)["noop"]
    # weights high but ranks moved: not a no-op
    shuffled = rng.permutation(raw)
    gate = n19.source_gate(raw, shuffled, np.full(300, 1.0), den)
    assert gate["n_moved"] >= 10 and not gate["noop"]
    assert n19.source_gate(raw[:40], raw[:40], np.ones(40), den[:40])["top"] == 40  # a smaller panel


def test_verdict_prunes_only_when_both_sources_are_no_ops():
    noop = {"gate": {"noop": True}}
    live = {"gate": {"noop": False}}
    assert n19.verdict({"a": noop, "b": noop}) == "prune"
    assert n19.verdict({"a": noop, "b": live}) == "run"


def _lists(obj) -> list:
    """Every list inside a nested record."""
    if isinstance(obj, dict):
        return [x for v in obj.values() for x in _lists(v)]
    return [obj] if isinstance(obj, list) else []


def test_source_record_reports_summaries_only(world, frames):
    frame = frames[M2G]
    fit = harness.t3_k4_fit(frame)
    rec = n19.source_record(frame, fit)
    assert rec["sections_used"] >= 2 and set(rec["gate"]) == {"top", "median_weight", "n_moved", "noop"}
    assert rec["robust_sd_ratio_to_raw"]["eb"] == pytest.approx(1.0, abs=0.2)
    assert "n_cells" not in json.dumps(rec)
    # no per-gene arrays in the report: the only list is the ten deciles' tau2
    assert [len(x) for x in _lists(rec)] == [10]


def test_n19_refuses_to_overwrite_and_reads_nothing_before_that(tmp_path):
    out = tmp_path / "n19"
    out.mkdir()
    (out / "something").write_text("x")
    with pytest.raises(SystemExit, match="already written"):
        n19.main(["--out", str(out), "--data", str(tmp_path / "no_data_here")])


def test_n19_runs_end_to_end_on_the_source_files_only(world, tmp_path, monkeypatch, capsys):
    # only the two counted folds' sources exist: a read of anything else (a reference, a target file) would raise
    data = tmp_path / "raw"
    for role, rel in RELEASED.items():
        (data / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(world.paths[role], data / rel)
    board = Board(key="T3:gata4", task="T3", genes=list(GENES), needs_coords=True, split="val", target_ko="gata4",
                  min_cells=100, max_cells=900, truth_cells=120, ref_cells=150)
    monkeypatch.setattr("vec.boards.load_boards", lambda *a, **k: {"T3:gata4": board})
    out = tmp_path / "out"
    n19.main(["--data", str(data), "--out", str(out)])
    res = json.loads((out / "n19_eb_pregate.json").read_text())
    assert set(res["sources"]) == {G2M, M2G} and res["verdict"] in ("prune", "run")
    assert set(res["inputs"]) == set(RELEASED.values())  # the five source files, each with its sha256
    assert all(len(sha) == 64 for sha in res["inputs"].values())
    assert res["stamp"]["commit"] and "tree_clean" in res["stamp"] and res["stamp"]["utc"].endswith("Z")
    for name in ("n19_eb_pregate.json", "n19_eb_pregate.txt"):
        assert (out / name).exists() and provenance.stamp_path(out / name).exists()
    text = capsys.readouterr().out
    assert "Verdict:" in text and "Inputs (sha256)" in text
    assert "cells" not in " ".join(k for s in res["sources"].values() for k in s)  # no cell count among the fields
    with pytest.raises(SystemExit, match="already written"):
        n19.main(["--data", str(data), "--out", str(out)])


# -------------------------------------------------------------------------------------------- K4's keep-or-stop rule
def _decision(base, keep=True, margins=None, gate=(), matches_panel=True):
    m = margins or {"response_identity": 0.5, "response_direction": 0.3, "response_magnitude": 0.1, "mutant_state": 0.0}
    fc = {"group_margins": m, "score_margin": sum(m.values()), "severity_gate_differs": list(gate),
          "matches_panel": matches_panel}
    return {"base": base, "change": "t3_factor_eb", "keep": keep, "reason": "fixture", "counted": [G2M, M2G], "need": 2,
            "fit_checks": {G2M: fc, M2G: fc}}


def test_k4_verdict_keeps_only_when_every_condition_holds():
    good = (_decision("t3_factor_raw"), _decision("t3_factor_eb_ctrl"))
    assert k4_verdict.verdict(*good)["kept"]
    # the adoption rule fails on one comparison
    v = k4_verdict.verdict(_decision("t3_factor_raw"), _decision("t3_factor_eb_ctrl", keep=False))
    assert not v["kept"] and v["parts"]["vs_ctrl"]["problems"] and not v["parts"]["vs_raw"]["problems"]
    # more than 40% of the margin on magnitude and mutant state (0.45 / 1.0)
    outside = {"response_identity": 0.3, "response_direction": 0.25, "response_magnitude": 0.3, "mutant_state": 0.15}
    assert not k4_verdict.verdict(_decision("t3_factor_raw", margins=outside), good[1])["kept"]
    edge = {"response_identity": 0.35, "response_direction": 0.25, "response_magnitude": 0.3, "mutant_state": 0.1}
    assert k4_verdict.verdict(_decision("t3_factor_raw", margins=edge), good[1])["kept"]  # exactly 0.40 passes
    # the severity gate differs in 3 of the seeds
    assert not k4_verdict.verdict(_decision("t3_factor_raw", gate=(1, 2, 3)), good[1])["kept"]
    assert k4_verdict.verdict(_decision("t3_factor_raw", gate=(1, 2)), good[1])["kept"]
    # a decide file with another base or a missing fold is refused as evidence
    assert not k4_verdict.verdict(_decision("t3_factor_eb_ctrl"), good[1])["kept"]
    missing = {**_decision("t3_factor_raw"), "counted": [G2M], "need": 1}
    assert not k4_verdict.verdict(missing, good[1])["kept"]
    # parts that do not reproduce the panel, or a fold without fit checks: the fit checks cannot be read, so not kept
    assert not k4_verdict.verdict(_decision("t3_factor_raw", matches_panel=False), good[1])["kept"]
    no_fc = _decision("t3_factor_raw")
    no_fc["fit_checks"] = {G2M: no_fc["fit_checks"][G2M]}
    assert not k4_verdict.verdict(no_fc, good[1])["kept"]


def _results_with_parts(key, gains, skills):
    """_k4_result plus the group skills and T3 parts that `run.py decide`'s fit checks read."""
    res = _k4_result(key, gains)
    for row in res["rows"]:
        row["groups"] = {c: dict(skills.get(c, skills["base"])) for c in res["candidates"]}
        row["t3"] = {c: {"de_direction_terms": [0.0] * 5, "severity": {"state": "scored", "gated": False},
                         "matches_panel": {"de_direction": True, "severity": True}} for c in res["candidates"]}
    return res


def _run_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location("harness_run_for_k4",
                                                  harness.REPO_ROOT / "experiments" / "harness" / "run.py")
    run = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run)
    return run


@pytest.mark.parametrize("on_target, expected_kept", [(True, True), (False, False)])
def test_decide_twice_then_the_verdict_through_the_real_cli(tmp_path, on_target, expected_kept):
    # the path of K4's row: `run.py decide` against the raw base and against the control, then k4_verdict on both files
    base = {g: 0.5 for g, _, _ in GROUPS["T3"]}
    moved = dict(base)
    moved.update({"response_identity": 0.51, "response_direction": 0.51} if on_target else {"response_magnitude": 0.55})
    skills = {"base": base, "t3_factor_eb": moved}
    runs = tmp_path / "runs"
    runs.mkdir()
    for key, name in ((G2M, "T3_Gata4_Mab21l2"), (M2G, "T3_Mab21l2_Gata4rep2"), (REP, "T3_Gata4rep1_rep2")):
        (runs / f"{name}.json").write_text(json.dumps(_results_with_parts(key, [0.5] * 10, skills)))
    run = _run_module()
    for stem, b in (("vs_raw", "t3_factor_raw"), ("vs_ctrl", "t3_factor_eb_ctrl")):
        run.main(["decide", "--runs", str(runs), "--board", "T3:bcatenin", "--base", b, "--change", "t3_factor_eb",
                  "--out", str(tmp_path / stem)])
    verdict_out = tmp_path / "verdict"
    k4_verdict.main(["--vs-raw", str(tmp_path / "vs_raw" / "decisions_T3_bcatenin.json"),
                     "--vs-ctrl", str(tmp_path / "vs_ctrl" / "decisions_T3_bcatenin.json"), "--out", str(verdict_out)])
    res = json.loads((verdict_out / "k4_verdict.json").read_text())
    assert res["kept"] is expected_kept
    for part in res["parts"].values():
        assert part["keep_by_adoption_rule"] is True and sorted(part["folds"]) == sorted([G2M, M2G])
        if not on_target:
            assert any("outside identity and direction" in p for p in part["problems"])
    # the replicate fold is reported by decide but never enters the verdict
    assert "replicate fold" in json.dumps(json.loads((tmp_path / "vs_raw" / "decisions_T3_bcatenin.json").read_text()))


def test_k4_verdict_main_writes_and_refuses_a_rerun(tmp_path):
    def write(name, dec, runs="runs/harness/t3_k4_v1"):
        p = tmp_path / name
        p.write_text(json.dumps({"board": "T3:bcatenin", "runs": runs, "decisions": [dec]}))
        return str(p)

    raw = write("raw.json", _decision("t3_factor_raw"))
    ctrl = write("ctrl.json", _decision("t3_factor_eb_ctrl"))
    out = tmp_path / "verdict"
    k4_verdict.main(["--vs-raw", raw, "--vs-ctrl", ctrl, "--out", str(out)])
    res = json.loads((out / "k4_verdict.json").read_text())
    assert res["kept"] is True and set(res["inputs"]) == {"vs_raw", "vs_ctrl"}
    assert res["stamp"]["commit"] and res["stamp"]["utc"].endswith("Z") and "tree_clean" in res["stamp"]
    with pytest.raises(SystemExit, match="already written"):
        k4_verdict.main(["--vs-raw", raw, "--vs-ctrl", ctrl, "--out", str(out)])
    other = write("ctrl2.json", _decision("t3_factor_eb_ctrl"), runs="runs/harness/other")
    out2 = tmp_path / "verdict2"
    k4_verdict.main(["--vs-raw", raw, "--vs-ctrl", other, "--out", str(out2)])
    assert json.loads((out2 / "k4_verdict.json").read_text())["kept"] is False  # two different run folders
