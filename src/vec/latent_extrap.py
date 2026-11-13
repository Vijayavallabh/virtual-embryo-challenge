"""B0 in a count model's latent space: the extrapolation counterpart of embryo row E21's learned generator, on the
two extrapolation boards (T1 and the heart), candidates `b0_lat_dec` and `b0_lat_ratio` of `vec.harness`.

B0 (`vec.harness._b0_t1`, `vec.harness._b0_t2`) draws real reference-stage cells at the count rule and moves each in
count space by B0_DAMPING of its category's mean change from the `start` stage to the reference. Here the same cells
(B0's own draw, without B0's count-space shift) are moved by the same logic in the latent space of a count model:

    z_new = z_own + damping * (mean z(reference cells of the category) - mean z(start cells of the category)),

z the encoder means. The categories, and which of them move, are B0's (`b0_moved_types`, the rule of
`vec.t1_scope.momentum_shift`: a category moves only when it has at least B0_MIN_CELLS cells at both stages and is not
represented by stand-in labels at either; the cells of a category at a stage are those momentum_shift averages). On T1
the categories are B0's chamber categories (obs["category"], with B0's per-file stand-ins); on the heart the harmonised
lineage obs["merged"], "Untracked" tagged per stage so that it never matches and never moves (as `_b0_t2_shift` does).
A cell of a category that does not move keeps its own latent exactly.

The count model is row E14's (`vec.t2_latent.train_nbvae`: scvi-tools SCVI, n_latent 20, gene-wise dispersion, NB
likelihood, observed library, no batch key, 400 epochs, batch size 128, seed 20261109), trained once per fold on the
counts of all cells of the `start` stage and of the reference stage (T1: every gene the harness hands over; heart: the
panel), start cells first. It is cached per fold in-process on the identity of the two input objects, as
`vec.harness._g_mom_factors` caches (the harness hands every seed the same objects), so the seeds vary the draw and the
sampling only. Device: the first visible CUDA device when torch sees one, else the CPU (`fit_device`; choose the GPU
with CUDA_VISIBLE_DEVICES). A GPU fit is not bit-reproducible across fits; the CPU path is (row E14).

Two arms, one form each, no option:

- decoded (`b0_lat_dec`, the class of row E21): counts sampled from the negative binomial at the decoded rate of the
  moved latent, each cell decoded at its real cell's own library (T1: the cell's total counts, since every gene is
  modelled; heart: its panel counts), then normalised as the released files are, log(1 + 1e4 c / N), N the sampled
  counts plus the real cell's counts outside the modelled genes (zero on T1), as `vec.t2_latent.decode_and_sample`
  does (a sampled cell without counts is sampled again, `vec.t2_build._keep_some_counts`). Stream [seed,
  DECODE_STREAM].
- ratio (`b0_lat_ratio`, the real cell kept): the drawn real cell's counts edited through `vec.core.apply_factors` by
  the per-gene factor decoder-normalised-mean(z_new) / decoder-normalised-mean(z_own) (both decoded at library 1, so
  the library cancels), expected depth kept, zeros kept; T1's held-flat genes (B0's list) get factor 1 and are left
  exactly as they are. The log factor is clipped at +-log 3 (CLIP, T11's `vec.t1_state_transport.CLIP_AT`); there is
  no zeroing threshold and no other clip. A cell whose latent did not move gets factor 1 without being decoded (its
  two means are the same by definition); a gene whose decoded mean is zero at both latents (a float underflow, never
  seen in practice) gets factor 1. Normalised as B0 normalises its moved cells, log1p CP10k of the edited counts
  (`vec.prep.log1p_cp10k`), so with every factor 1 the arm is B0 at damping 0 to the last bit. Streams [seed,
  RATIO_STREAM, chunk], one per CELL_CHUNK cells in draw order.

Coordinates (heart), the count, the draw, the composition and the scope are B0's.

torch and scvi-tools are imported only when a model is trained (`vec.t2_latent.train_nbvae`) or a device is chosen, so
this module imports in the project's main environment.
"""

from __future__ import annotations

import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import anndata as ad
import numpy as np
from scipy import sparse

from vec import t2_build, t2_latent
from vec.core import apply_factors, counts_of
from vec.prep import log1p_cp10k

# Constants this module adds to row E14's (all of whose constants are vec.t2_latent's, unchanged)
B0_MIN_CELLS = 20  # vec.t1_scope.momentum_shift's default min_cells, which B0 uses (a test pins the two equal)
CLIP = float(np.log(3.0))  # the ratio arm's clip of |log factor|: T11's vec.t1_state_transport.CLIP_AT (a test pins it)
RATIO_STREAM = 21  # apply_factors' stream of the ratio arm, apart from B0's 1 and 3, g_mom's 5 and E14's 20
DECODE_STREAM = 22  # the decoded arm's NB sampling stream
CELL_CHUNK = 1024  # cells per decoder call and per apply_factors call (memory; it also cuts the ratio arm's streams)
ARMS = ("decoded", "ratio")
T2_LINEAGE_KEY = "merged"
T2_UNTRACKED = "Untracked"
T2_TYPE_KEY = "b0_type"  # as vec.harness._b0_t2_shift names it
RULE = ("latent B0: z_new = z_own + damping * (mean z(ref cells of the category) - mean z(start cells of the "
        "category)) for B0's moved categories, z the encoder means of row E14's count model fitted once per fold on "
        "the start and reference stages; decoded: NB counts at the decoded rate of z_new and the real cell's own "
        "library, log1p CP10k with the counts outside the modelled genes; ratio: the real cell's counts times "
        "decoded-mean(z_new) / decoded-mean(z_own) via apply_factors (expected depth kept, zeros kept, held-flat "
        "genes 1), log factor clipped at +-log 3")


@dataclass
class FoldModel:
    """A fold's count model: `model` trained on the start stage's cells then the reference's (`n_start` of the
    former), on `genes`, and the fit's record (device, wall time, final training loss)."""

    model: t2_latent.NBVAE
    n_start: int
    genes: list[str]
    record: dict = field(default_factory=dict)

    @property
    def z_start(self) -> np.ndarray:
        return np.asarray(self.model.latent[:self.n_start], dtype=np.float64)

    @property
    def z_ref(self) -> np.ndarray:
        return np.asarray(self.model.latent[self.n_start:], dtype=np.float64)


_FOLD: dict = {}  # fold_model's cache: one fold, keyed on the identity of its two input objects and the genes


def fit_device() -> str:
    """'cuda:0' (the first device CUDA_VISIBLE_DEVICES exposes) when torch sees a GPU, else 'cpu'."""
    import torch

    return "cuda:0" if torch.cuda.is_available() else "cpu"


def fold_model(ref: ad.AnnData, start: ad.AnnData, genes: Sequence[str]) -> FoldModel:
    """The fold's count model, trained once per pair of input objects (`ref` and `start` as the harness hands them,
    never copies) on all their cells' counts on `genes`, start cells first, with row E14's settings
    (`vec.t2_latent.train_nbvae` at TRAIN_SEED, its defaults otherwise) on `fit_device()`."""
    genes = list(map(str, genes))
    hit = _FOLD.get("fold")
    if hit is None or hit["inputs"][0] is not ref or hit["inputs"][1] is not start or hit["genes"] != genes:
        S, R = counts_of(start, genes), counts_of(ref, genes)
        device = fit_device()
        history: dict = {}
        t0 = time.time()
        model = t2_latent.train_nbvae(sparse.vstack([S, R]).tocsr(), t2_latent.TRAIN_SEED, device=device,
                                      history=history)
        wall = time.time() - t0
        if model.latent.shape[0] != S.shape[0] + R.shape[0]:
            raise ValueError("the model's latent means do not match its training cells")
        record = {**dict(model.record), "device": device, "fit_seconds": round(wall, 1), "history": history,
                  "final_train_loss": history.get("train_loss"), "n_genes": len(genes)}
        _FOLD.clear()
        hit = _FOLD["fold"] = {"inputs": (ref, start), "genes": genes,
                               "fold": FoldModel(model=model, n_start=S.shape[0], genes=genes, record=record)}
    return hit["fold"]


# ------------------------------------------------------------------------------------------------ B0's categories
def b0_moved_types(start: ad.AnnData, ref: ad.AnnData, types: Sequence[str], key: str, *,
                   label_key: str = "celltype", start_fallback: Mapping[str, Sequence[str]] | None = None,
                   ref_fallback: Mapping[str, Sequence[str]] | None = None,
                   min_cells: int = B0_MIN_CELLS) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """The categories among `types` that B0's momentum moves, each with its cells at the start stage and at the
    reference (row indices of `start` and `ref`): exactly `vec.t1_scope.momentum_shift`'s rule with B0's arguments (no
    `within`, its stand-ins, standin_momentum False): the rows are `_type_rows`', and a category moves only with at
    least min_cells rows at both stages, none of them stand-ins."""
    from vec.t1_scope import _type_rows

    out = {}
    for t in sorted(set(map(str, types))):
        s_rows, s_how = _type_rows(start, t, key, label_key, None, start_fallback, min_cells)
        r_rows, r_how = _type_rows(ref, t, key, label_key, None, ref_fallback, min_cells)
        if len(s_rows) < min_cells or len(r_rows) < min_cells:
            continue
        if s_how.startswith("stand-in") or r_how.startswith("stand-in"):
            continue
        out[t] = (np.asarray(s_rows), np.asarray(r_rows))
    return out


def moved_latent(fold: FoldModel, rows: np.ndarray, cell_types: np.ndarray,
                 moved: Mapping[str, tuple[np.ndarray, np.ndarray]], damping: float) -> tuple[np.ndarray, np.ndarray,
                                                                                               dict]:
    """(z_own, z_new, record) for the drawn reference rows `rows` (their categories `cell_types`): each moved
    category's cells shifted by damping times its mean latent displacement from the start stage to the reference;
    every other cell keeps its own latent exactly."""
    z_ref, z_start = fold.z_ref, fold.z_start
    z_own = z_ref[np.asarray(rows, dtype=np.int64)]
    z_new = z_own.copy()
    cell_types = np.asarray(cell_types).astype(str)
    disp = {}
    for t, (s_rows, r_rows) in moved.items():
        d = z_ref[r_rows].mean(0) - z_start[s_rows].mean(0)
        disp[t] = d
        m = cell_types == t
        z_new[m] = z_own[m] + float(damping) * d
    step = np.abs(z_new - z_own)
    is_moved = np.isin(cell_types, sorted(moved))
    record = {"damping": float(damping), "categories_moved": sorted(moved), "n_categories_moved": len(moved),
              "share_cells_moved": float(is_moved.mean()) if len(is_moved) else 0.0,
              "mean_abs_latent_step": float(step.mean()) if step.size else 0.0,
              "mean_abs_latent_step_moved_cells": float(step[is_moved].mean()) if is_moved.any() else 0.0,
              "displacement_norm": {t: float(np.linalg.norm(d)) for t, d in disp.items()}}
    return z_own, z_new, record


# ------------------------------------------------------------------------------------------------ the two arms
def decoded_counts(model: t2_latent.NBVAE, z: np.ndarray, library: np.ndarray, other: np.ndarray,
                   rng: np.random.Generator, chunk: int = CELL_CHUNK) -> tuple[np.ndarray, sparse.csr_matrix, int]:
    """The decoded arm: NB counts sampled (never a decoded mean) at the decoder's rate for latent z and the real
    cells' libraries, CELL_CHUNK cells per decoder call (`vec.t2_latent.decode_and_sample` decodes all cells at once,
    too large on T1's genes), normalised as log(1 + 1e4 c / N), N = sampled counts + `other`. Returns (expression,
    dense float32; the sampled counts; how many cells were sampled again for having no counts)."""
    z = np.asarray(z, dtype=np.float64)
    lib = np.asarray(library, dtype=np.float64)
    extra = np.asarray(other, dtype=np.float64)
    if lib.ndim != 1 or lib.shape != extra.shape or len(lib) != len(z):
        raise ValueError("z, library and other must describe the same cells")
    if (lib <= 0).any():
        raise ValueError("a cell has no library to decode at")

    def make(idx, g):
        parts = []
        for lo in range(0, len(idx), chunk):
            sel = idx[lo:lo + chunk]
            mu, theta = model.decode(z[sel], lib[sel])
            mu, th = np.asarray(mu, dtype=np.float64), np.asarray(theta, dtype=np.float64)
            if mu.shape[0] != len(sel):
                raise ValueError("the decoder returned the wrong number of cells")
            parts.append(sparse.csr_matrix(t2_latent.sample_nb(mu, th, g), dtype=np.int64))
        C = sparse.vstack(parts).tocsr() if parts else sparse.csr_matrix((0, 0), dtype=np.int64)
        return C, np.asarray(C.sum(1)).ravel().astype(np.float64) + extra[idx]

    C, N, retried = t2_build._keep_some_counts(make, len(lib), rng)
    return t2_build._log_normalise(C, N), C, retried


def ratio_counts(model: t2_latent.NBVAE, C, z_own: np.ndarray, z_new: np.ndarray, held: np.ndarray, seed: Sequence[int],
                 chunk: int = CELL_CHUNK) -> tuple[sparse.csr_matrix, dict]:
    """The ratio arm: each real cell's counts `C` (cells x genes, the model's genes) edited by
    `vec.core.apply_factors` (keep_depth) with the per-gene factor decoded-mean(z_new) / decoded-mean(z_own), the log
    clipped at +-CLIP, held genes 1; cells whose latent is unchanged get factor 1 undecoded. Stream [*seed, chunk
    index] per CELL_CHUNK cells. Returns the edited counts and a record (the log factor's spread over the entries it
    acts on: the non-held genes with a count in a moved cell)."""
    C = sparse.csr_matrix(C, dtype=np.int64)
    z_own, z_new = np.asarray(z_own, dtype=np.float64), np.asarray(z_new, dtype=np.float64)
    held = np.asarray(held, dtype=bool)
    n, g = C.shape
    if z_own.shape != z_new.shape or len(z_own) != n or held.shape != (g,):
        raise ValueError("counts, latents and the held mask must describe the same cells and genes")
    moved = np.any(z_new != z_own, axis=1)
    parts, acting, both_zero, restored = [], [], 0, 0
    for k, lo in enumerate(range(0, n, chunk)):
        hi = min(lo + chunk, n)
        Ck = C[lo:hi]
        mk = np.flatnonzero(moved[lo:hi])
        logf = np.zeros((len(mk) + 1, g))  # one row per moved cell, then a row of ones (log 0) for the rest
        if len(mk):
            ones = np.ones(len(mk))
            mu_new, _ = model.decode(z_new[lo + mk], ones)
            mu_own, _ = model.decode(z_own[lo + mk], ones)
            mu_new, mu_own = np.asarray(mu_new, dtype=np.float64), np.asarray(mu_own, dtype=np.float64)
            with np.errstate(divide="ignore", invalid="ignore"):
                lf = np.log(mu_new) - np.log(mu_own)
            zero = (mu_new == 0) & (mu_own == 0)
            both_zero += int(zero.sum())
            lf[zero] = 0.0
            lf = np.clip(lf, -CLIP, CLIP)
            lf[:, held] = 0.0
            logf[:len(mk)] = lf
            sub = Ck[mk].tocoo()
            acting.append(lf[sub.row, sub.col][~held[sub.col]])
        row = np.full(hi - lo, len(mk), dtype=np.int64)
        row[mk] = np.arange(len(mk))
        new, _, dead = apply_factors(Ck, np.exp(logf), row, held, keep_depth=True, seed=[*map(int, seed), k])
        restored += int(dead.sum())
        parts.append(new)
    out = sparse.vstack(parts).tocsr() if parts else C.copy()
    L = np.concatenate(acting).astype(np.float64) if acting else np.zeros(0)
    record = {"n_cells_moved": int(moved.sum()), "n_acting_entries": int(L.size),
              "sd_log_factor": float(L.std()) if L.size else 0.0,
              "p1_log_factor": float(np.percentile(L, 1)) if L.size else 0.0,
              "p99_log_factor": float(np.percentile(L, 99)) if L.size else 0.0,
              "n_at_clip": int((np.abs(L) >= CLIP - 1e-12).sum()), "n_both_zero_means": both_zero,
              "restored_cells": restored, "clip": CLIP}
    return out, record


def _libraries_record(C_real, C_out) -> dict:
    lr = np.asarray(sparse.csr_matrix(C_real).sum(1)).ravel()
    lo = np.asarray(sparse.csr_matrix(C_out).sum(1)).ravel()
    size = float(np.prod(C_real.shape)) or 1.0
    return {"median_library_out": float(np.median(lo)) if len(lo) else 0.0,
            "median_library_real": float(np.median(lr)) if len(lr) else 0.0,
            "zero_share_out": 1.0 - sparse.csr_matrix(C_out).nnz / size,
            "zero_share_real": 1.0 - sparse.csr_matrix(C_real).nnz / size}


def _arm(arm: str, fold: FoldModel, C: sparse.csr_matrix, z_own, z_new, held, other, seed: int):
    """Expression (cells x the fold's genes) of one arm for the drawn cells' real counts C, and its record."""
    if arm == "decoded":
        lib = np.asarray(C.sum(1)).ravel().astype(np.float64)
        X, Cs, retried = decoded_counts(fold.model, z_new, lib, other, np.random.default_rng([int(seed),
                                                                                             DECODE_STREAM]))
        return X, {"arm": arm, "retried": int(retried), **_libraries_record(C, Cs)}
    if arm == "ratio":
        new, rec = ratio_counts(fold.model, C, z_own, z_new, held, [int(seed), RATIO_STREAM])
        X = log1p_cp10k(new)
        return X, {"arm": arm, **rec, **_libraries_record(C, new)}
    raise ValueError(f"arm must be one of {ARMS}, got {arm!r}")


# ------------------------------------------------------------------------------------------------ the two boards
def t1_population(ref_in: ad.AnnData, start_in: ad.AnnData, ref: ad.AnnData, start: ad.AnnData,
                  shares: Mapping[str, float], n_cells: int, seed: int, held_flat: Sequence[str], *, ref_source: str,
                  start_source: str, damping: float, arm: str) -> tuple[sparse.spmatrix | np.ndarray, list[str], dict]:
    """T1: B0's draw (`vec.t1_scope.reweight` exactly as `vec.t1_state_transport.b0_population` and B0's build call
    it, the same seed and stand-ins; no shift), the fold's model on the inputs `ref_in` and `start_in` (the harness's
    objects, the cache key), B0's moved chamber categories, and one arm. ref and start: shallow copies of the inputs,
    classified in place as B0 classifies them. Returns (expression, cells x genes; the genes; record)."""
    from vec import t1_scope

    spec = t1_scope.load_spec()
    t1_scope.classify(ref, t1_scope.LABEL_KEY, ref_source, spec)
    fallback = spec.fallbacks(ref_source)
    pop, _ = t1_scope.reweight(ref, shares, n_cells, seed, label_key=t1_scope.LABEL_KEY, within=None,
                               fallback=fallback, spec=spec)
    t1_scope.classify(start, t1_scope.LABEL_KEY, start_source, spec)
    if not ref.obs_names.is_unique:
        raise ValueError("the reference's cell names must be unique to find the drawn rows")
    rows = ref.obs_names.get_indexer(pop.obs["source_obs"].astype(str))
    if (rows < 0).any():
        raise ValueError("a drawn cell is not a reference cell")
    cats = pop.obs[t1_scope.CATEGORY_KEY].astype(str).to_numpy()
    moved = b0_moved_types(start, ref, cats, t1_scope.CATEGORY_KEY, label_key=t1_scope.LABEL_KEY,
                           start_fallback=spec.fallbacks(start_source), ref_fallback=fallback)
    genes = list(map(str, ref.var_names))
    fold = fold_model(ref_in, start_in, genes)
    z_own, z_new, lat = moved_latent(fold, rows, cats, moved, damping)
    C = counts_of(pop, genes)
    held_set = set(map(str, held_flat))
    held = np.array([g in held_set for g in genes], dtype=bool)
    X, rec = _arm(arm, fold, C, z_own, z_new, held, np.zeros(C.shape[0]), seed)
    return X, genes, {"rule": RULE, "latent": lat, "expression": rec, "n_held_genes": int(held.sum()),
                      "fold": fold.record}


def _typed(a: ad.AnnData, tag: str) -> ad.AnnData:
    """A shallow copy with B0's T2 type (obs["merged"], "Untracked" tagged by stage so it never matches)."""
    m = np.asarray(a.obs[T2_LINEAGE_KEY]).astype(str)
    obs = a.obs.copy()
    obs[T2_TYPE_KEY] = np.where(m == T2_UNTRACKED, f"{T2_UNTRACKED}@{tag}", m)
    return ad.AnnData(X=a.X, obs=obs, var=a.var, layers=dict(a.layers), obsm=dict(a.obsm))


def t2_population(pool: ad.AnnData, start: ad.AnnData, rows: np.ndarray, seed: int, *, damping: float,
                  arm: str) -> tuple[np.ndarray | sparse.spmatrix, list[str], dict]:
    """The heart (a T2 extrapolation board): B0's drawn rows of the reference pool (`vec.harness.draw_copy`, passed
    in), the fold's model on `pool` and `start` (the harness's objects) on the pool's genes (the panel), B0's moved
    lineages (obs["merged"], `_b0_t2_shift`'s typing, no gene held), and one arm. The decoded arm's library is the real
    cell's panel counts and its normaliser adds the real cell's counts outside the panel (`vec.t2_build._stage_counts`,
    as `vec.t2_latent`). Returns (expression, cells x genes; the genes; record)."""
    genes = list(map(str, pool.var_names))
    tp, ts_ = _typed(pool, "end"), _typed(start, "start")
    rows = np.asarray(rows, dtype=np.int64)
    types = np.asarray(tp.obs[T2_TYPE_KEY]).astype(str)[rows]
    moved = b0_moved_types(ts_, tp, types, T2_TYPE_KEY)
    fold = fold_model(pool, start, genes)
    z_own, z_new, lat = moved_latent(fold, rows, types, moved, damping)
    _, Cpool, _, other = t2_build._stage_counts(pool, genes)
    C = Cpool[rows]
    X, rec = _arm(arm, fold, C, z_own, z_new, np.zeros(len(genes), dtype=bool), other[rows], seed)
    rec["any_counts_outside_panel"] = bool((other[rows] > 0).any())
    return X, genes, {"rule": RULE, "latent": lat, "expression": rec, "n_held_genes": 0, "fold": fold.record}


__all__ = ["ARMS", "B0_MIN_CELLS", "CELL_CHUNK", "CLIP", "DECODE_STREAM", "RATIO_STREAM", "RULE", "FoldModel",
           "b0_moved_types", "decoded_counts", "fit_device", "fold_model", "moved_latent", "ratio_counts",
           "t1_population", "t2_population"]
