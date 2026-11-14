"""Count-model latent interpolation for the T2 embryo board: embryo try E14 (`fgw_latent`, experiments/harness/TRIES.md)
and the parts its target-free fidelity gate N18 shares with it (experiments/harness/n18_latent_fidelity.py).

fgw's mixture is real cells of the two brackets, and a whole cell cannot be an intermediate state; averaging counts is
too narrow a population (try E8). This module trains a negative-binomial variational autoencoder (scvi-tools `SCVI`:
n_latent 20, gene-wise dispersion, NB likelihood, observed library, no batch key) once per task on the panel counts of
the two brackets, encodes every cell (encoder means), and for each output cell of fgw's coupling-paired draw takes

    z_t = (1 - t) z_b + t z_a + sqrt(2 t (1 - t) (1 - rho)) eps,   eps ~ N(0, S),

z_b and z_a the before and after cell of the pair (the output cell and its coupled partner), S = (1 - t) S_b + t S_a the
brackets' within-bracket covariance of the encoder means, and rho the coupling-weighted correlation of the two brackets'
latent means over the coupling's anchors (`coupling_correlation`). The noise term restores the latent variance that
averaging two coupled cells loses: a pair average has variance (1 + rho) / 2 of V at t = 1/2 for correlation rho, the
failure E8 measured. One form, no sweep. The decoder's NB mean is evaluated at the drawn real cell's panel library (the
library of the real cell at fgw's row, so the output population's depth distribution is fgw's own and the closed depth
line is not touched); counts are sampled from the NB, never decoded means (`docs/PLAN.md`, principle 1), and normalised
as the released files are, log(1 + 1e4 c / N) with N the sampled panel counts plus the real cell's counts outside the
panel. Positions, rows and partners are fgw's (`vec.t2_build.interpolate_clouds`).

The model is trained once per task and per input hash (`task_model`, cached per process), so the harness's 10 seeds vary
the pairs and the NB draws only and the training randomness is unsampled; the training seed is fixed (TRAIN_SEED). CPU
(16 threads): only a CPU build can pass the determinism check of two builds at the freeze commit.

Everything except `train_nbvae` is numpy: torch and scvi-tools (the `gpu` dependency group) are imported only inside it,
so this module and the rest of the package import without them. Tests give `train_nbvae` a stand-in (`NBVAE`).
"""

from __future__ import annotations

import hashlib
import tempfile
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
from scipy import sparse

from vec import t2_build

# Fixed in rows N18 and E14 of experiments/harness/TRIES.md before any result
N_LATENT = 20
TRAIN_SEED = 20261109  # one model per task; the build seeds vary the pairs and NB draws only
TRAIN_THREADS = 16
MAX_EPOCHS = 400  # scvi-tools' own default cap (its heuristic gives 400 at about 20,000 cells)
BATCH_SIZE = 128  # scvi-tools' default
EXPRESSION_STREAM = 20  # random stream of the expression step, apart from the interpolation's and the other steps'
COV_RIDGE = 1e-6  # of the mean diagonal: keeps a near-singular latent covariance factorisable
DECODE_CHUNK = 4096  # cells per decoder call (memory only)
PROBE_CELLS = 64  # train_nbvae's check that the decoder returns NB means at the library it was given
LATENT_RULE = ("latent: z_t = (1 - t) z_b + t z_a + sqrt(2 t (1 - t) (1 - rho)) eps, eps ~ N(0, S), rho the "
               "coupling-weighted correlation of the anchors' latent means, S = (1 - t) S_b + t S_a the within-bracket "
               "latent covariance; NB rate decoded at the drawn real cell's panel library, counts sampled, normalised "
               "with the real cell's counts outside the panel")


@dataclass
class NBVAE:
    """What the interpolation needs from a trained count model: `latent`, the encoder means of its training cells in
    training order (the two brackets' cells, before then after), and `decode(z, library)`, which returns the NB means
    (cells x genes) at latent z and panel library sizes `library`, and the NB inverse dispersions (genes, or cells x
    genes). `train_nbvae` builds it from an scvi-tools model; tests give a stand-in."""

    latent: np.ndarray
    decode: Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]
    record: dict = field(default_factory=dict)
    encode: Callable[[object], np.ndarray] | None = None  # counts (cells x the training genes) -> encoder means


def _batch_size(n_cells: int, batch_size: int) -> int:
    """scvi-tools' batch norm fails on a last batch of 1-2 cells (its data splitter warns about it): one more cell per
    batch moves such a remainder, deterministically."""
    b = int(batch_size)
    while n_cells % b in (1, 2):
        b += 1
    return b


def train_nbvae(counts, seed: int = TRAIN_SEED, *, n_latent: int = N_LATENT, max_epochs: int = MAX_EPOCHS,
                batch_size: int = BATCH_SIZE, threads: int = TRAIN_THREADS, device: str = "cpu",
                history: dict | None = None) -> NBVAE:
    """Train the negative-binomial VAE on a cells x genes count matrix (the two brackets' panel counts, before then
    after): scvi-tools `SCVI` with `n_latent`, NB likelihood, gene-wise dispersion, observed library, no batch key, on
    every cell (no validation split, no early stopping, the fixed seed; `device` "cpu" or "cuda:<i>", the CPU path being
    the one the determinism check can pass). Returns its encoder means and decoder as an `NBVAE`.

    torch and scvi-tools are imported only here: the `gpu` dependency group, not installed in the project's main
    environment (docs/COMPUTE.md). history (default None: nothing): a dict that receives the last epoch's value of each
    training curve scvi-tools logs (e.g. "train_loss") and the epochs run; the model and its record are the same
    with or without it (vec.latent_extrap reads the final training loss through it)."""
    import scvi
    import torch

    if device != "cpu" and not device.startswith("cuda"):
        raise ValueError(f"device must be 'cpu' or 'cuda[:<i>]', got {device!r}")
    C = sparse.csr_matrix(counts, dtype=np.float32)
    C.eliminate_zeros()
    if (np.asarray(C.sum(1)).ravel() <= 0).any():
        raise ValueError("a cell has no panel counts: the observed library would be zero")
    torch.set_num_threads(int(threads))
    scvi.settings.seed = int(seed)
    adata = ad.AnnData(X=C)
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata, n_latent=int(n_latent), gene_likelihood="nb", dispersion="gene")
    accelerator, devices = (("cpu", 1) if device == "cpu"
                            else ("gpu", [int(device.split(":")[1])] if ":" in device else 1))
    with tempfile.TemporaryDirectory() as tmp:  # anything the trainer might write stays out of the repository
        model.train(max_epochs=int(max_epochs), accelerator=accelerator, devices=devices, train_size=1.0,
                    batch_size=_batch_size(C.shape[0], batch_size), early_stopping=False, enable_progress_bar=False,
                    enable_checkpointing=False, default_root_dir=tmp)
    if history is not None:
        for name, curve in model.history.items():
            if len(curve):
                history[name] = float(np.asarray(curve, dtype=np.float64).ravel()[-1])
        history["epochs_run"] = int(max((len(c) for c in model.history.values()), default=0))
    module = model.module
    module.eval()
    dev = next(module.parameters()).device
    latent = np.asarray(model.get_latent_representation(adata, give_mean=True, batch_size=DECODE_CHUNK),
                        dtype=np.float64)

    def decode(z, library):
        zt, lib = np.asarray(z, dtype=np.float32), np.asarray(library, dtype=np.float64)
        if (lib <= 0).any():
            raise ValueError("the decoder needs positive libraries")
        mus, theta = [], None
        with torch.no_grad():
            for lo in range(0, len(zt), DECODE_CHUNK):
                hi = lo + DECODE_CHUNK
                z_t = torch.as_tensor(zt[lo:hi], device=dev)
                log_lib = torch.as_tensor(np.log(lib[lo:hi]), dtype=torch.float32, device=dev).unsqueeze(1)
                batch = torch.zeros((len(z_t), 1), dtype=torch.int64, device=dev)  # one batch: no batch key
                px = module.generative(z=z_t, library=log_lib, batch_index=batch)["px"]  # NB(mean mu, theta)
                mus.append(px.mu.detach().cpu().numpy())
                if theta is None:
                    theta = px.theta.detach().cpu().numpy()
        mu = np.concatenate(mus, axis=0).astype(np.float64)
        th = np.asarray(theta, dtype=np.float64)
        return mu, (th[0] if th.ndim == 2 else th)

    # The decoder's scale sums to one per cell, so its NB means must sum to the library they were decoded at: this
    # catches a misread of the module's library convention before any cell is generated.
    n_probe = min(PROBE_CELLS, latent.shape[0])
    probe_lib = np.asarray(C[:n_probe].sum(1)).ravel().astype(np.float64)
    probe_mu, _ = decode(latent[:n_probe], probe_lib)
    if not np.allclose(probe_mu.sum(1), probe_lib, rtol=1e-3):
        raise RuntimeError("the decoder's NB means do not sum to the given library; the module's convention changed")
    record = {"n_latent": int(n_latent), "gene_likelihood": "nb", "dispersion": "gene", "max_epochs": int(max_epochs),
              "batch_size": _batch_size(C.shape[0], batch_size), "seed": int(seed), "threads": int(threads),
              "device": device, "scvi_tools": str(getattr(scvi, "__version__", "unknown")),
              "torch": str(torch.__version__)}
    def encode(counts) -> np.ndarray:
        """Encoder means of other cells' counts on the training genes, in the training columns' order (vec.latent_extrap
        encodes the stages before the fold's `start` with it)."""
        X = sparse.csr_matrix(counts, dtype=np.float32)
        if X.shape[1] != C.shape[1]:
            raise ValueError("encode needs counts on the model's training genes")
        return np.asarray(model.get_latent_representation(ad.AnnData(X=X), give_mean=True, batch_size=DECODE_CHUNK),
                          dtype=np.float64)

    return NBVAE(latent=latent, decode=decode, record=record, encode=encode)


def input_key(count_matrices: Sequence, genes: Sequence[str]) -> str:
    """sha256 of the panel count matrices (in order) and the gene names: what `task_model` keys its cache on."""
    h = hashlib.sha256()
    for C in count_matrices:
        C = sparse.csr_matrix(C, dtype=np.int64)
        C.eliminate_zeros()
        C.sort_indices()
        for part in (np.asarray(C.shape, dtype=np.int64), C.indptr.astype(np.int64), C.indices.astype(np.int64),
                     C.data.astype(np.int64)):
            h.update(part.tobytes())
    h.update("\n".join(str(g) for g in genes).encode())
    return h.hexdigest()


_MODELS: dict[str, NBVAE] = {}  # task_model's cache: one trained model per set of inputs, per process


def task_model(before: ad.AnnData, after: ad.AnnData, genes: Sequence[str], *, seed: int = TRAIN_SEED,
               threads: int = TRAIN_THREADS, device: str = "cpu") -> NBVAE:
    """The task's model: trained once per process on the two brackets' panel counts (before then after, the harness's
    loaded inputs, no target), cached on their sha256 (`input_key`), the seed and the settings. The 10 seeds of a
    harness run therefore share one model and vary only the pairs and the NB draws; training randomness is unsampled
    (row E14 says so)."""
    Cb = sparse.csr_matrix(before[:, list(genes)].layers["counts"], dtype=np.int64)
    Ca = sparse.csr_matrix(after[:, list(genes)].layers["counts"], dtype=np.int64)
    key = input_key((Cb, Ca), genes) + f"|{seed}|{threads}|{device}"
    if key not in _MODELS:
        model = train_nbvae(sparse.vstack([Cb, Ca]).tocsr(), seed, threads=threads, device=device)
        if model.latent.shape[0] != Cb.shape[0] + Ca.shape[0]:
            raise ValueError("the model's latent means do not match its training cells")
        _MODELS[key] = model
    return _MODELS[key]


def coupling_correlation(z_before, z_after, coupling) -> float:
    """rho of the variance term: the correlation of the two brackets' latent means under the anchors' coupling, i.e.
    tr(C) / sqrt(tr(V_b) tr(V_a)) where C is the coupling-weighted cross-covariance and V_b, V_a the covariances under
    the coupling's marginals (z_before and z_after: the anchors' latent means; coupling: anchors x anchors mass, the
    one that set the moves and partners). |rho| <= 1 by construction (Cauchy-Schwarz on the coupling)."""
    Zb, Za = np.asarray(z_before, dtype=np.float64), np.asarray(z_after, dtype=np.float64)
    T = np.asarray(coupling, dtype=np.float64)
    if T.shape != (Zb.shape[0], Za.shape[0]) or (T < 0).any():
        raise ValueError("the coupling must be a non-negative anchors x anchors matrix")
    total = float(T.sum())
    if not total > 0:
        raise ValueError("the coupling carries no mass")
    T = T / total
    mb, ma = T.sum(1), T.sum(0)
    Xb, Xa = Zb - mb @ Zb, Za - ma @ Za
    cross = float(((Xb.T @ T) * Xa.T).sum())
    vb, va = float((mb * (Xb * Xb).sum(1)).sum()), float((ma * (Xa * Xa).sum(1)).sum())
    if not (vb > 0 and va > 0):
        raise ValueError("the anchors' latent means have no variance under the coupling")
    return float(np.clip(cross / np.sqrt(vb * va), -1.0, 1.0))


def latent_covariance(z) -> np.ndarray:
    """The covariance of a bracket's latent means over its cells (ddof 0, rows = cells)."""
    Z = np.asarray(z, dtype=np.float64)
    X = Z - Z.mean(0)
    return X.T @ X / max(len(X), 1)


def _covariance_root(cov) -> np.ndarray:
    """A lower-triangular L with L L^T = cov plus a small ridge (COV_RIDGE of the mean diagonal): exact and
    deterministic for a positive definite cov, stable for a near-singular one."""
    S = (np.asarray(cov, dtype=np.float64) + np.asarray(cov, dtype=np.float64).T) / 2.0
    ridge = COV_RIDGE * max(float(np.trace(S)) / S.shape[0], 1e-12)
    return np.linalg.cholesky(S + ridge * np.eye(S.shape[0]))


def interpolate_latent(z_before, z_after, t: float, rho: float, cov, rng: np.random.Generator) -> np.ndarray:
    """z_t = (1 - t) z_b + t z_a + sqrt(2 t (1 - t) (1 - rho)) eps, eps ~ N(0, cov), one independent eps per row. The
    linear term alone has variance (1 - t)^2 V_b + t^2 V_a + 2 t (1 - t) rho sqrt(V_b V_a), so with V_b = V_a = V the
    noise term restores V exactly. At t = 0 or 1, or for rho = 1, the noise term is zero: the cells' own latents."""
    Zb, Za = np.asarray(z_before, dtype=np.float64), np.asarray(z_after, dtype=np.float64)
    if Zb.shape != Za.shape or Zb.ndim != 2:
        raise ValueError("z_before and z_after must be the same cells x latent matrix shape")
    if not 0.0 <= t <= 1.0:
        raise ValueError(f"t must lie in [0, 1], got {t}")
    if not -1.0 <= rho <= 1.0:
        raise ValueError(f"rho must lie in [-1, 1], got {rho}")
    z = (1.0 - t) * Zb + t * Za
    w = 2.0 * t * (1.0 - t) * (1.0 - rho)
    if w <= 0:
        return z
    L = _covariance_root(cov)
    return z + np.sqrt(w) * rng.standard_normal(z.shape) @ L.T


def sample_nb(mu, theta, rng: np.random.Generator) -> np.ndarray:
    """Counts from NB(mean mu, inverse dispersion theta), by its Gamma-Poisson form: lambda ~ Gamma(theta, mu / theta),
    c ~ Poisson(lambda). Numpy's generator, so the draw is deterministic given the seed. mu: cells x genes; theta: genes
    or cells x genes, positive."""
    m = np.asarray(mu, dtype=np.float64)
    th = np.asarray(theta, dtype=np.float64)
    if (m < 0).any() or not np.isfinite(m).all():
        raise ValueError("NB means must be finite and non-negative")
    if not (th > 0).all() or not np.isfinite(th).all():
        raise ValueError("NB inverse dispersions must be finite and positive")
    lam = rng.gamma(shape=th, scale=m / th)
    return rng.poisson(lam)


def _streams(seed) -> tuple[np.random.Generator, np.random.Generator]:
    """Two independent generators from `seed` (an int or a sequence of ints): the latent noise and the NB sampling."""
    s = [int(x) for x in seed] if isinstance(seed, (list, tuple, np.ndarray)) else [int(seed)]
    return np.random.default_rng([*s, 0]), np.random.default_rng([*s, 1])


def decode_and_sample(model: NBVAE, z, library, other, rng: np.random.Generator) -> tuple[np.ndarray, int]:
    """The last step of the latent candidate and of its gate: NB means decoded at latent z and the real cells' panel
    libraries (`library`), counts sampled (never a decoded mean), and the cells normalised as the released files are,
    log(1 + 1e4 c / N), N = sampled panel counts + `other` (the real cells' counts outside the panel). A sampled cell
    with no panel count is sampled again, as `vec.t2_build._keep_some_counts` does for thinned cells. Returns
    expression (cells x genes, float32) and how many cells were sampled more than once."""
    lib = np.asarray(library, dtype=np.float64)
    extra = np.asarray(other, dtype=np.float64)
    if lib.ndim != 1 or lib.shape != extra.shape or len(lib) != len(z):
        raise ValueError("z, library and other must describe the same cells")
    if (lib <= 0).any():
        raise ValueError("a cell has no panel library to decode at")
    mu, theta = model.decode(np.asarray(z, dtype=np.float64), lib)
    mu = np.asarray(mu, dtype=np.float64)
    if mu.shape[0] != len(lib):
        raise ValueError("the decoder returned the wrong number of cells")
    th = np.asarray(theta, dtype=np.float64)

    def make(idx, g):
        C = sparse.csr_matrix(sample_nb(mu[idx], th if th.ndim == 1 else th[idx], g), dtype=np.int64)
        return C, np.asarray(C.sum(1)).ravel().astype(np.float64) + extra[idx]

    C, N, retried = t2_build._keep_some_counts(make, len(lib), rng)
    return t2_build._log_normalise(C, N), retried


def latent_midpoint_counts(model: NBVAE, before: ad.AnnData, after: ad.AnnData, interp: t2_build.Interpolation, genes,
                           t: float, seed) -> tuple[np.ndarray, dict]:
    """Try E14's expression: for every output cell of `interp` (`vec.t2_build.interpolate_clouds` with
    `with_partners=True` and `with_coupling=True`; output cell i is row `rows[i]` of stage `stage[i]`, its partner a
    row of the other stage), the pair's latent means interpolated with the variance-preserving term (LATENT_RULE) and
    NB counts sampled at the drawn real cell's panel library. `model`: the task's model (`task_model`), trained on
    `before` then `after`. Returns expression (cells x genes, float32, log(1 + 1e4 c / N)) and a record."""
    if interp.partners is None or interp.coupling is None:
        raise ValueError("latent_midpoint_counts needs interpolate_clouds(..., with_partners=True, with_coupling=True)")
    if not 0.0 <= t <= 1.0:
        raise ValueError(f"t must lie in [0, 1], got {t}")
    (_, Cb, lb, ob), (_, Ca, la, oa) = (t2_build._stage_counts(a, genes) for a in (before, after))
    n_b, n_a = Cb.shape[0], Ca.shape[0]
    Z = np.asarray(model.latent, dtype=np.float64)
    if Z.shape[0] != n_b + n_a:
        raise ValueError("the model's latent means are not this task's two brackets (before then after)")
    Zb, Za = Z[:n_b], Z[n_b:]
    c = interp.coupling
    rho = coupling_correlation(Zb[c["anchors_before"]], Za[c["anchors_after"]], c["moves"])
    cov = (1.0 - t) * latent_covariance(Zb) + t * latent_covariance(Za)
    rng_eps, rng_nb = _streams(seed)
    is_b = np.asarray(interp.stage) == 0
    rows, partners = np.asarray(interp.rows, dtype=np.int64), np.asarray(interp.partners, dtype=np.int64)
    rb, ra = np.where(is_b, rows, partners), np.where(is_b, partners, rows)
    zt = interpolate_latent(Zb[rb], Za[ra], t, rho, cov, rng_eps)
    n = len(rows)
    panel_b, panel_a = lb - ob, la - oa  # the real cells' panel libraries
    library, other = np.empty(n), np.empty(n)
    library[is_b], other[is_b] = panel_b[rows[is_b]], ob[rows[is_b]]
    library[~is_b], other[~is_b] = panel_a[rows[~is_b]], oa[rows[~is_b]]
    X, retried = decode_and_sample(model, zt, library, other, rng_nb)
    return X, {"rule": LATENT_RULE, "t": float(t), "rho": float(rho), "retried": int(retried),
               "library": "each output cell's panel library is the drawn real cell's (fgw's row)"}
