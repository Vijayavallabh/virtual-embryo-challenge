"""Tries K4 and K4c and diagnostic N19 (experiments/harness/TRIES.md): a gene-wise empirical-Bayes shrinkage of a T3
source factor. The estimator only: pure array code, no file is read here.

The source factor is K2a's estimator (`vec.harness.t3_source_frame`, t3_decompose.gene_factor): a plug-in ratio of
w-weighted mean within-cell proportions over the strict non-Mesp1 labels, knockout file(s) against their matched wild
type (eps 1e-5, clip 0.25-4), used at exponent 1. Replicates of one file agree only at Spearman 0.72 (rule D11a), and
that within-file instability enters the factor as if it were signal. The estimator shrinks each gene's log ratio
towards f = 1 by how unstable the gene's numerator is across the knockout file's own sections. Pinned in K4's row
before any fold is read:

  y_g    the log ratio log((num_g + eps) / (den_g + eps)) of gene_factor, before its rescale and clip.
  v_g    the delete-one-section jackknife variance of y_g, (K - 1) / K x sum_k (y_g^(-k) - mean_k y_g^(-k))^2, numerator
         side only: the knockout cells are deleted by section level; the denominator (at test the WT side is the shared
         anchor and is not shrunk), the label mask, the labels' weights and the harmoniser calls are those of the full
         estimate. A section is a z-level of a knockout file (vec.spatial.snap_z per file); pooled replicates give
         replicate-by-level ids. K counts the sections holding at least one of the estimator's cells. A label a deletion
         empties keeps its full-sample mean. A gene with v_g = 0 (no measured noise: no count in the numerator) takes
         its decile's median positive variance; a decile with none takes B_g = 1.
  tau2   per decile of the factor's own denominator den (the w-weighted WT mean proportion; ties in panel order),
         DerSimonian-Laird with the prior mean fixed at 0 (f = 1): tau2_b = max(0, (Q_b - k_b) / sum_b w_g),
         Q_b = sum_b w_g y_g^2, w_g = 1 / v_g.
  B_g    tau2_b / (tau2_b + v_g), the posterior mean being B_g y_g (B_g = 0 where tau2_b = 0).
  c      one constant: the raw factor's robust SD over the panel genes divided by the robust SD of B y (robust SD =
         1.4826 x the median absolute deviation from the median; c = 1 if the robust SD of B y is 0). It is computed on
         the EB arm and applied to both arms, so shrinkage changes which genes move and the raw factor's amplitude is
         kept (severity's beta otherwise falls with it).
  phi    exp(c B y), then the last two steps of gene_factor: rescaled so the WT-weighted total proportion is unchanged,
         then clipped at the rules' [0.25, 4]. The raw factor stays gene_factor's own output (t3_source_frame).

The control (K4c) permutes the weights B across genes within the same deciles of den (K2b's logic: one generator,
default_rng(seed), the deciles from the lowest up, `decile_permutation`; seed fixed in K4c's row), with the same c. It
answers "both folds reward any shrinkage of a factor that lacks the shared WT reference": it moves as many genes by the
same amounts, at random within those deciles. No second prior, variance model or gene filter is tried (K4's row).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy import sparse

from vec.spatial import snap_z

BINS = 10  # deciles of the factor's denominator (tau2 and the control's permutation)
MAD_SCALE = 1.4826  # robust SD = MAD_SCALE x median absolute deviation from the median
DEFAULT_CTRL_SEED = 20261109  # K4c: the one permutation of the weights (its row; vec.harness.K4_CTRL_SEED passes it)


@dataclass(frozen=True)
class Jackknife:
    y: np.ndarray  # the log ratio of gene_factor, before its rescale and clip, on all of the estimator's cells
    v: np.ndarray  # its delete-one-section jackknife variance (numerator side)
    den: np.ndarray  # gene_factor's denominator: the w-weighted WT mean proportion over the masked labels
    n_blocks: int  # K: the sections holding at least one of the estimator's cells


@dataclass(frozen=True)
class EbFit:
    raw: np.ndarray  # K2a's estimator's factor, as t3_decompose.gene_factor returned it
    eb: np.ndarray  # K4: exp(c B y), rescaled and clipped as gene_factor
    ctrl: np.ndarray  # K4c: the same with B permuted within the deciles of den (same c)
    weights: np.ndarray  # B, the gene-wise shrinkage weights
    weights_ctrl: np.ndarray  # B permuted (K4c)
    constant: float  # c, the amplitude constant (computed on the EB arm)
    tau2: np.ndarray  # one per decile of den, lowest first
    jack: Jackknife


def block_ids(z_by_file: Sequence[np.ndarray]) -> tuple[np.ndarray, int]:
    """Each cell's section block over the knockout files, in file order, and the number of blocks over the files. One
    block per file and snapped z-level (`snap_z` per file, its levels 0..n-1), so pooled replicates get
    replicate-by-level ids."""
    blocks, offset = [], 0
    for z in z_by_file:
        lev = snap_z(np.asarray(z, dtype=np.float64)).index
        blocks.append(lev + offset)
        offset += int(lev.max()) + 1 if lev.size else 0
    return (np.concatenate(blocks) if blocks else np.zeros(0, dtype=np.int64)), offset


def decile_bins(basis, bins: int = BINS) -> list[np.ndarray]:
    """Gene positions cut into `bins` runs of consecutive ranks of `basis` (a stable sort, ties in panel order), sizes
    differing by at most one: the deciles of K4's tau2 and of its control's permutation."""
    return np.array_split(np.argsort(np.asarray(basis, dtype=np.float64), kind="stable"), bins)


def decile_permutation(x, basis, seed: int = DEFAULT_CTRL_SEED, bins: int = BINS) -> np.ndarray:
    """`x` permuted across genes within the bins of `basis`: K2b's logic (vec.harness.t3_decile_permutation), one
    generator, default_rng(seed), permuting each bin's values in turn, from the lowest basis up."""
    x = np.asarray(x, dtype=np.float64)
    rng = np.random.default_rng(seed)
    out = x.copy()
    for genes in decile_bins(basis, bins):
        out[genes] = x[genes[rng.permutation(len(genes))]]
    return out


def jackknife(C, lab, block, mu_wt, w, mask, eps) -> Jackknife:
    """gene_factor's log ratio y on the estimator's cells and its delete-one-section jackknife variance v, numerator
    side (the module docstring). C: knockout cells x genes counts; lab: each cell's label index (-1: none); block: each
    cell's section block; mu_wt, w, mask, eps: gene_factor's mu_wt (labels x genes), label weights, label mask and
    eps."""
    mu_wt = np.asarray(mu_wt, dtype=np.float64)
    ww = (np.asarray(w, dtype=np.float64) * np.asarray(mask, dtype=np.float64))[:, None]
    den = (ww * mu_wt).sum(0)
    if den.sum() <= 0:
        raise ValueError("the factor's mask holds no WT mass, so there is no factor to shrink")
    L, G = mu_wt.shape
    lab, block = np.asarray(lab), np.asarray(block)
    labelled = lab >= 0
    keep = labelled.copy()
    keep[labelled] = np.asarray(mask, dtype=bool)[lab[labelled]]
    sections, bi = np.unique(block[keep], return_inverse=True)
    K = sections.size
    if K < 2:
        raise ValueError(f"the jackknife needs at least two sections holding the estimator's cells, found {K}")
    Ck = np.asarray(C)[keep]
    props = Ck / np.maximum(Ck.sum(1, keepdims=True), 1)  # the proportions of t3_decompose.props_type_means
    idx = lab[keep].astype(np.int64) * K + bi
    n = np.bincount(idx, minlength=L * K).reshape(L, K)
    sums = np.asarray(sparse.csr_matrix((np.ones(idx.size), (idx, np.arange(idx.size))), shape=(L * K, idx.size))
                      @ props).reshape(L, K, G)
    n_t, s_t = n.sum(1), sums.sum(1)
    mu_full = np.divide(s_t, n_t[:, None], out=np.zeros_like(s_t), where=n_t[:, None] > 0)
    ys = np.empty((K, G))
    for k in range(K):
        rem = n_t - n[:, k]
        mu_k = np.divide(s_t - sums[:, k], rem[:, None], out=mu_full.copy(), where=rem[:, None] > 0)
        ys[k] = np.log(((ww * mu_k).sum(0) + eps) / (den + eps))
    v = (K - 1) / K * ((ys - ys.mean(0)) ** 2).sum(0)
    y = np.log(((ww * mu_full).sum(0) + eps) / (den + eps))
    return Jackknife(y=y, v=v, den=den, n_blocks=K)


def eb_weights(y, v, den, bins: int = BINS) -> tuple[np.ndarray, np.ndarray]:
    """The shrinkage weights B and the per-decile tau2 (the module docstring): DerSimonian-Laird with the prior mean
    fixed at 0, per decile of den, lowest first."""
    y, v = np.asarray(y, dtype=np.float64), np.asarray(v, dtype=np.float64)
    B = np.empty_like(y)
    tau2 = []
    for genes in decile_bins(den, bins):
        vb = v[genes].copy()
        pos = vb > 0
        if not pos.any():  # no noise measured in this decile: no shrinkage
            B[genes] = 1.0
            tau2.append(0.0)
            continue
        vb[~pos] = np.median(vb[pos])  # a gene with no measured noise takes the decile's median variance
        w = 1.0 / vb
        t2 = max(0.0, float((w * y[genes] ** 2).sum() - genes.size) / float(w.sum()))
        B[genes] = t2 / (t2 + vb) if t2 > 0 else 0.0
        tau2.append(t2)
    return B, np.asarray(tau2)


def robust_sd(x) -> float:
    """1.4826 x the median absolute deviation from the median."""
    x = np.asarray(x, dtype=np.float64)
    return MAD_SCALE * float(np.median(np.abs(x - np.median(x))))


def finalize(log_phi, den, clip=(0.25, 4.0)) -> np.ndarray:
    """exp(log_phi), then gene_factor's last two steps: rescaled so the WT-weighted total proportion is unchanged, then
    clipped."""
    den = np.asarray(den, dtype=np.float64)
    phi = np.exp(np.asarray(log_phi, dtype=np.float64))
    phi = phi / (den @ phi / max(den.sum(), 1e-12))
    return np.clip(phi, *clip)


def eb_factors(phi_raw, jack: Jackknife, clip=(0.25, 4.0), seed: int = DEFAULT_CTRL_SEED, bins: int = BINS) -> EbFit:
    """K4's and K4c's factors from K2a's raw factor and the jackknife: B, c, then exp(c B y) rescaled and clipped as
    gene_factor; the control permutes B within the deciles of den (same c)."""
    raw = np.array(phi_raw, dtype=np.float64)
    weights, tau2 = eb_weights(jack.y, jack.v, jack.den, bins)
    weights_ctrl = decile_permutation(weights, jack.den, seed, bins)
    m, m_ctrl = weights * jack.y, weights_ctrl * jack.y
    sd_m = robust_sd(m)
    constant = robust_sd(np.log(raw)) / sd_m if sd_m > 0 else 1.0
    return EbFit(raw=raw, eb=finalize(constant * m, jack.den, clip), ctrl=finalize(constant * m_ctrl, jack.den, clip),
                 weights=weights, weights_ctrl=weights_ctrl, constant=float(constant), tau2=tau2, jack=jack)
