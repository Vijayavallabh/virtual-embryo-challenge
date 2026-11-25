"""The world blur: the heart's world-blur hedge, file WBH of row R69 (experiments/harness/TRIES.md), and the arms of its
gate, row H37 (experiments/harness/h37_heart_world_blur.py). HEART_MOON2_design.md's P1 with every fix of its cold
challenge (HEART_MOON2_challenge.md, section 5): the prior's unexpressed genes take the median spread (fix 3), the
scorer's pseudobulk is held gene by gene (fix 2), G00's kernel-scale check (fix 4a) and the KO-against-replicate check
(fix 4c).

It re-measures a built file's cells in many plausible measurement worlds at once. Only expression moves: rows, obs
names, variables, coordinates and the cell number are the input's, and each cell's expected depth is kept.

  prior   per file pair (E8.75 WT-G4r1, WT-G4r2, G4r1-G4r2; E9.5 WT-Mab21l2 KO; the captain's words in row R69 allow
          these three T3 files, with their matched wild types, as inputs to the world prior only) and per merged
          lineage group that is not Mesp1-derived (configs/lineage_map.yaml, mesp1 "no", the released labels) with at
          least 200 cells on both sides: d_g = log(mean CP10k_A + 0.5) - log(mean CP10k_B + 0.5); a pair's profile is
          the median over its groups; a pair is used when it has 2 or more groups. s_g^2 is the mean over the used
          pairs of the squared profile. A gene with mean CP10k under 1 on both sides in every used group of every used
          pair is flagged (the Mesp1 exclusion would leave the heart's own genes almost unblurred) and takes the median
          s of the other genes; every other gene's log s is shrunk halfway to that median.
  counts  from X = log1p CP10k: v = expm1(x); per cell N = 1e4 m / min positive v for the first m of 1, 2, 3 that makes
          c = v N / 1e4 whole to 1e-3 in every entry, else refusal (CountRecoveryError).
  worlds  blocks (file WBH and arms W, P): KMeans(32, n_init 10, random_state = the seed) on the coordinates; blocks
          paired by size rank; pair j draws z_j ~ N(0, I) on stream [seed, 63]; its blocks take +z_j and -z_j.
          cells (arm U): random antithetic cell pairs and their z on stream [seed, 65]; an odd cell out keeps its
          counts. P: s permuted across genes on stream [seed, 64]. Factor f_g = exp(s_g z_g - s_g^2 / 2).
  depth   a per-cell scalar k (bisection) makes the expected count total of f' = k f equal the cell's own total T.
  re-measure  f' < 1: Binomial(c, f'). f' > 1: c + Poisson((f' - 1) lambda), lambda = (c + T pbar) / 2, pbar the mean
          proportions of the cell's 10 nearest cells in 20 PCs of its X (so a gene can switch on).
  hold    (fix 2) the scorer's pseudobulk is the mean of X = log1p(1e4 c / N) with each cell's N fixed. Per gene, one
          scalar a in [0, 8] on the top-up rate of every f' > 1 cell brings the gene's mean X to the input's within
          1e-3; failing that, one scalar b in [0, 1] on the thinning depth (1 - f') of every f' < 1 cell, the top-up at
          its own rate; failing both, the gene keeps the input's counts. The draws are coupled (a Poisson process for
          the top-up, one uniform per count for the thinning), so each mean moves monotonically in its scalar and the
          scalar's plateau nearest the input's mean is taken. A cell left empty keeps its counts (the genes are then
          held again with that cell fixed). Refused if any gene ends more than 1e-3 from the input's mean X.
  X       entries whose count changed take float32(log1p(1e4 c' / N)); every other entry keeps its input bytes.
  depth   the mean library of the re-measure before the hold (a = b = 1) is what the depth scalar keeps, and G0's and
          the read's 2% bound reads it (library_ratio_design). The hold itself raises the count library by the Jensen
          gap of log1p (a thinned 1 -> 0 costs about five times what a top-up 1 -> 2 gives back): 4 to 10% at s 0.2 to
          0.5 with the thinning scalar first, 10 to 38% with the top-up scalar first, on synthetic counts at MERFISH
          depth. That figure (library_ratio) is recorded and decides nothing: the scorer sees X, whose per-gene mean
          the hold keeps within 1e-3, and a 2% bound on it would stop G0 for the fix's own reason.

Every draw comes from a stream named by the seed, so the same input, prior and seed give the same bytes, and a prior
of s = 0 gives the input's X exactly. Nothing here reads a target or a score.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import sparse

PAIRS = (("T2/E8.75.h5ad", "T3/E8.75_gata4_ko_1.h5ad"), ("T2/E8.75.h5ad", "T3/E8.75_gata4_ko_2.h5ad"),
         ("T3/E8.75_gata4_ko_1.h5ad", "T3/E8.75_gata4_ko_2.h5ad"), ("T2/E9.5.h5ad", "T3/E9.5_mab21l2_ko.h5ad"))
PAIR_NAMES = ("E8.75 WT-G4r1", "E8.75 WT-G4r2", "E8.75 G4r1-G4r2", "E9.5 WT-Mab21l2")
REPLICATE = PAIR_NAMES[2]  # the one genotype-matched pair
KO_WT = (PAIR_NAMES[0], PAIR_NAMES[1], PAIR_NAMES[3])
PRIOR_FILES = ("T2/E8.75.h5ad", "T3/E8.75_gata4_ko_1.h5ad", "T3/E8.75_gata4_ko_2.h5ad", "T2/E9.5.h5ad",
               "T3/E9.5_mab21l2_ko.h5ad")  # the only files the prior reads (row R69); no E10.5 cell enters it
NON_MESP1 = "no"  # configs/lineage_map.yaml's mesp1 flag of the groups the prior reads
UNTRACKED = "Untracked"  # vec.prep.add_lineage's merged value for a group without a stable cross-stage meaning
MIN_GROUP_CELLS = 200
MIN_GROUPS = 2  # a pair is used when it has at least this many groups
PSEUDO = 0.5  # the profile's pseudo-count, in CP10k
FLAG_CP10K = 1.0  # fix 3: a gene under this mean CP10k on both sides in every used group is flagged
SHRINK = 0.5  # log s is shrunk halfway to the median
LIBRARY = 1e4
N_BLOCKS = 32
KMEANS_N_INIT = 10
KNN_PCS, KNN_K = 20, 10
TOPUP_MAX = 8.0  # fix 2: the top-up scalar's range is [0, 8]
PB_TOL = 1e-3  # fix 2: each gene's mean X within this of the input's
RECOVERY_TOL = 1e-3
MIN_COUNTS = (1, 2, 3)  # the smallest count a cell is taken to hold, tried in this order
STREAM_WORLDS, STREAM_PERM, STREAM_CELLS, STREAM_G00, STREAM_COUNTS = 63, 64, 65, 66, 67
BLOCKS, CELLS = "blocks", "cells"
MODES = (BLOCKS, CELLS)
# G00 (fix 4a), in the scorer's geometry (common.core_metrics.mmd_unbiased: PCA 30 on the truth, gamma0 = 1 / the
# median squared distance, scales 0.25 to 4)
G00_CELLS, G00_PCS, G00_DRAWS = 2000, 30, 16
G00_SCALES = (0.25, 0.5, 1.0, 2.0, 4.0)
G00_MIN_CUT = 0.10


class CountRecoveryError(ValueError):
    """A cell's counts cannot be recovered from X exactly."""


class HoldError(ValueError):
    """A gene's mean X could not be held to the input's."""


# ------------------------------------------------------------------------------------------------ the prior
@dataclass
class Prior:
    """The per-gene world spread s and how it was formed (every figure from the released prior files)."""

    genes: tuple
    s: np.ndarray  # the spread each gene's world factor is drawn with
    raw: np.ndarray  # sqrt of the mean over the used pairs of the squared profile (nan when no pair is used)
    flagged: np.ndarray  # fix 3's unexpressed genes, which take the median s
    profiles: dict = field(default_factory=dict)  # used pair -> its profile over the genes
    groups: dict = field(default_factory=dict)  # every pair -> the groups it has (used or not)

    @property
    def used(self) -> list[str]:
        return [p for p in PAIR_NAMES if p in self.profiles]

    def median_s(self) -> float:
        return float(np.median(self.s))

    def spread_of(self, pairs) -> np.ndarray:
        """Per gene, sqrt of the mean over `pairs` (used ones only) of the squared profile."""
        ps = [self.profiles[p] for p in pairs if p in self.profiles]
        if not ps:
            return np.full(len(self.genes), np.nan)
        return np.sqrt(np.mean(np.square(ps), axis=0))

    def summary(self) -> dict:
        """The prior's figures (released files only; no target, no test file)."""
        keep = ~self.flagged
        ko = self.spread_of(KO_WT)[keep]
        rep = np.abs(self.profiles[REPLICATE])[keep] if REPLICATE in self.profiles else np.array([np.nan])
        return {"pairs_used": self.used, "groups": {p: list(g) for p, g in self.groups.items()},
                "n_genes": len(self.genes), "flagged": int(self.flagged.sum()), "median_s": self.median_s(),
                "ko_wt_median_spread": float(np.median(ko)) if ko.size else float("nan"),
                "replicate_median_spread": float(np.median(rep)) if rep.size else float("nan")}


def group_means(counts, merged, mesp1) -> dict[str, tuple[int, np.ndarray]]:
    """Per merged group that is not Mesp1-derived: its cell number and its mean CP10k per gene (cells without a count
    left out)."""
    C = sparse.csr_matrix(counts, dtype=np.float64)
    tot = np.asarray(C.sum(1)).ravel()
    merged, mesp1 = np.asarray(merged).astype(str), np.asarray(mesp1).astype(str)
    ok = (mesp1 == NON_MESP1) & (merged != UNTRACKED) & (tot > 0)
    out = {}
    for g in sorted(set(merged[ok])):
        rows = np.flatnonzero(ok & (merged == g))
        cp = sparse.diags(LIBRARY / tot[rows]) @ C[rows]
        out[g] = (len(rows), np.asarray(cp.mean(0)).ravel())
    return out


def prior_from_tables(tables: dict[str, dict], genes) -> Prior:
    """The prior from each file's group means (group_means), keyed by the file names of PAIRS."""
    genes = tuple(genes)
    ng = len(genes)
    profiles, groups, low = {}, {}, []
    for name, (fa, fb) in zip(PAIR_NAMES, PAIRS, strict=True):
        ta, tb = tables[fa], tables[fb]
        both = [g for g in sorted(set(ta) & set(tb)) if ta[g][0] >= MIN_GROUP_CELLS and tb[g][0] >= MIN_GROUP_CELLS]
        groups[name] = both
        if len(both) < MIN_GROUPS:
            continue
        d = np.stack([np.log(ta[g][1] + PSEUDO) - np.log(tb[g][1] + PSEUDO) for g in both])
        profiles[name] = np.median(d, axis=0)
        low += [(ta[g][1] < FLAG_CP10K) & (tb[g][1] < FLAG_CP10K) for g in both]
    if not profiles:
        nan = np.full(ng, np.nan)
        return Prior(genes, nan.copy(), nan, np.zeros(ng, dtype=bool), profiles, groups)
    flagged = np.all(np.stack(low), axis=0)
    raw = np.sqrt(np.mean(np.square(list(profiles.values())), axis=0))
    keep = ~flagged
    logs = np.log(np.maximum(raw, 1e-6))
    med = float(np.median(logs[keep])) if keep.any() else 0.0
    s = np.exp((1 - SHRINK) * logs + SHRINK * med)
    s[flagged] = np.exp(med)
    return Prior(genes, s, raw, flagged, profiles, groups)


def load_prior(data_dir: str | Path, genes) -> Prior:
    """The prior from PRIOR_FILES under data_dir, each loaded as the heart builder loads its inputs (duplicate records
    dropped, the lineage map's groups from the released labels), on the board's genes."""
    from vec.core import counts_of
    from vec.prep import load_stage

    tables = {}
    for f in PRIOR_FILES:
        a = load_stage(Path(data_dir) / f, genes=list(genes), drop_duplicates=True)
        tables[f] = group_means(counts_of(a, list(genes)), a.obs["merged"], a.obs["mesp1"])
        del a
    return prior_from_tables(tables, genes)


# ------------------------------------------------------------------------------------------------ counts from X
def _dense(X) -> np.ndarray:
    return np.asarray(X.toarray() if sparse.issparse(X) else X)


def recover_counts(X) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """From X = log1p CP10k: the counts (float64, whole), each cell's N and the smallest count it was taken to hold (0
    for a cell without a count). CountRecoveryError unless every cell recovers to RECOVERY_TOL with one of
    MIN_COUNTS."""
    x = _dense(X).astype(np.float64)
    if not np.isfinite(x).all() or (x < 0).any():
        raise CountRecoveryError("X holds a negative or non-finite value")
    v = np.expm1(x)
    vmin = np.where(v > 0, v, np.inf).min(1)
    n = x.shape[0]
    counts, N, used = np.zeros_like(v), np.zeros(n), np.zeros(n, dtype=np.int64)
    todo = np.isfinite(vmin)
    for m in MIN_COUNTS:
        if not todo.any():
            break
        rows = np.flatnonzero(todo)
        Nm = LIBRARY * m / vmin[rows]
        c = v[rows] * (Nm / LIBRARY)[:, None]
        ok = np.abs(c - np.rint(c)).max(1) <= RECOVERY_TOL
        good = rows[ok]
        counts[good], N[good], used[good] = np.rint(c[ok]), Nm[ok], m
        todo[good] = False
    if todo.any():
        raise CountRecoveryError(f"{int(todo.sum())} cells do not recover to whole counts with a smallest count of "
                                 f"{', '.join(map(str, MIN_COUNTS))}")
    return counts, N, used


# ------------------------------------------------------------------------------------------------ the worlds
@dataclass
class Worlds:
    """Each cell's world (a row of z, -1 for none) and its sign; z holds one row per antithetic pair."""

    unit: np.ndarray
    sign: np.ndarray
    z: np.ndarray
    mode: str
    pairs: list  # blocks: (block +z, block -z) per pair; cells: the number of pairs

    def factors(self, s: np.ndarray) -> np.ndarray:
        """f = exp(sign s z - s^2 / 2) per cell and gene; 1 for a cell without a world."""
        s = np.asarray(s, dtype=np.float64)
        F = np.ones((len(self.unit), len(s)))
        w = self.unit >= 0
        F[w] = np.exp(self.sign[w, None] * s[None, :] * self.z[self.unit[w]] - 0.5 * s[None, :] ** 2)
        return F

    def antithetic(self) -> bool:
        """Every pair holds exactly its two members, one +z and one -z (blocks: whole blocks)."""
        u, sg = self.unit[self.unit >= 0], self.sign[self.unit >= 0]
        if self.mode == CELLS:
            return bool(np.array_equal(np.bincount(u[sg > 0], minlength=len(self.z)), np.ones(len(self.z)))
                        and np.array_equal(np.bincount(u[sg < 0], minlength=len(self.z)), np.ones(len(self.z))))
        return bool(len(self.pairs) == len(self.z) and all(a != b for a, b in self.pairs)
                    and len({b for p in self.pairs for b in p}) == 2 * len(self.pairs))


def block_worlds(coords, seed: int, n_genes: int) -> Worlds:
    from sklearn.cluster import KMeans

    C = np.asarray(coords, dtype=np.float64)
    if C.shape[0] < N_BLOCKS:
        raise ValueError(f"{C.shape[0]} cells cannot form {N_BLOCKS} blocks")
    lab = KMeans(N_BLOCKS, n_init=KMEANS_N_INIT, random_state=seed).fit(C).labels_
    size = np.bincount(lab, minlength=N_BLOCKS)
    order = sorted(range(N_BLOCKS), key=lambda b: (-int(size[b]), b))
    pairs = [(order[2 * j], order[2 * j + 1]) for j in range(N_BLOCKS // 2)]
    z = np.random.default_rng([seed, STREAM_WORLDS]).standard_normal((len(pairs), n_genes))
    unit, sign = np.full(len(lab), -1), np.zeros(len(lab))
    for j, (bp, bm) in enumerate(pairs):
        unit[lab == bp], sign[lab == bp] = j, 1.0
        unit[lab == bm], sign[lab == bm] = j, -1.0
    return Worlds(unit, sign, z, BLOCKS, pairs)


def cell_worlds(n: int, seed: int, n_genes: int) -> Worlds:
    rng = np.random.default_rng([seed, STREAM_CELLS])
    perm = rng.permutation(n)
    k = n // 2
    z = rng.standard_normal((k, n_genes))
    unit, sign = np.full(n, -1), np.zeros(n)
    unit[perm[0:2 * k:2]], sign[perm[0:2 * k:2]] = np.arange(k), 1.0
    unit[perm[1:2 * k:2]], sign[perm[1:2 * k:2]] = np.arange(k), -1.0
    return Worlds(unit, sign, z, CELLS, [k])


def permuted(s: np.ndarray, seed: int) -> np.ndarray:
    """Arm P's spread: s permuted across genes on stream [seed, 64] (the same sizes, the wrong pattern)."""
    s = np.asarray(s)
    return s[np.random.default_rng([seed, STREAM_PERM]).permutation(len(s))]


# ------------------------------------------------------------------------------------------------ the re-measure
def neighbour_proportions(x: np.ndarray, counts: np.ndarray, totals: np.ndarray) -> np.ndarray:
    """pbar: the mean count proportions of each cell's KNN_K nearest other cells in KNN_PCS PCs of its X."""
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    n = x.shape[0]
    k = min(KNN_K, n - 1)
    npc = min(KNN_PCS, x.shape[1], n - 1)
    emb = PCA(npc, svd_solver="full", random_state=0).fit_transform(x)
    idx = NearestNeighbors(n_neighbors=k + 1).fit(emb).kneighbors(emb, return_distance=False)
    own = idx == np.arange(n)[:, None]
    keep = ~own
    keep[~own.any(1), -1] = False
    nb = idx[keep].reshape(n, k)
    props = np.divide(counts, totals[:, None], out=np.zeros_like(counts), where=totals[:, None] > 0)
    return props[nb].mean(1)


def _expected_total(k, c, F, lam):
    f = k[:, None] * F
    return np.where(f < 1, c * f, c + (f - 1) * lam).sum(1)


def depth_scalars(c: np.ndarray, F: np.ndarray, lam: np.ndarray, T: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """Per cell of `rows`, the scalar k whose factors k F keep the expected count total at T (bisection); 1
    elsewhere."""
    k = np.ones(len(T))
    if not len(rows):
        return k
    cc, FF, ll, tt = c[rows], F[rows], lam[rows], T[rows]
    lo, hi = np.zeros(len(rows)), np.ones(len(rows))
    for _ in range(80):
        short = _expected_total(hi, cc, FF, ll) < tt
        if not short.any():
            break
        hi = np.where(short, hi * 2, hi)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        under = _expected_total(mid, cc, FF, ll) < tt
        lo, hi = np.where(under, mid, lo), np.where(under, hi, mid)
    k[rows] = 0.5 * (lo + hi)
    return k


def _ranks(cell: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Each event's rank (1, 2, ...) among its own cell's events in order of position."""
    if not len(cell):
        return np.zeros(0, dtype=np.int64)
    o = np.lexsort((pos, cell))
    cs = cell[o]
    start = np.r_[0, np.flatnonzero(np.diff(cs)) + 1]
    first = np.repeat(start, np.diff(np.r_[start, len(cs)]))
    r = np.empty(len(cell), dtype=np.int64)
    r[o] = np.arange(len(cs)) - first + 1
    return r


def _plateau(base: float, pos: np.ndarray, delta: np.ndarray, target: float, n: int) -> tuple[int, float, np.ndarray]:
    """The number of events (taken in order of position) whose mean is nearest the target: (m, |mean - target|, the
    order of the events)."""
    o = np.argsort(pos, kind="stable")
    means = (base + np.r_[0.0, np.cumsum(delta[o])]) / n
    m = int(np.argmin(np.abs(means - target)))
    return m, float(abs(means[m] - target)), o


def _val(cnt, N):
    return np.log1p(LIBRARY * cnt / N)


def hold_gene(c: np.ndarray, fp: np.ndarray, lam: np.ndarray, N: np.ndarray, live: np.ndarray,
              rng: np.random.Generator) -> tuple[np.ndarray, str, np.ndarray]:
    """One gene's re-measured counts with its mean X held to the input's (fix 2): (counts, how, design), how one of
    "topup" (the top-up scalar a in [0, 8], the thinning whole), "thin" (the thinning scalar b in [0, 1], the top-up at
    its own rate) and "kept" (neither holds: the input's counts); design is the re-measure before the hold (a = b = 1),
    whose count total the depth scalar keeps. Only `live` cells move."""
    n = len(c)
    target = float(_val(c, N).sum()) / n
    up = np.flatnonzero(live & (fp > 1) & (lam > 0))
    dn = np.flatnonzero(live & (fp < 1) & (c > 0))
    K = rng.poisson(TOPUP_MAX * (fp[up] - 1) * lam[up])
    ucell = np.repeat(up, K)
    upos = TOPUP_MAX * rng.random(len(ucell))
    ci = c[dn].astype(np.int64)
    dcell = np.repeat(dn, ci)
    dpos = (1 - rng.random(len(dcell))) / (1 - fp[dcell])
    urank, drank = _ranks(ucell, upos), _ranks(dcell, dpos)
    dropped = dpos <= 1
    thin = c.copy()
    np.subtract.at(thin, dcell[dropped], 1)
    design = thin.copy()
    np.add.at(design, ucell[upos <= 1], 1)
    # the top-up scalar, the thinning whole (b = 1)
    base = float(_val(thin, N).sum())
    udelta = _val(c[ucell] + urank, N[ucell]) - _val(c[ucell] + urank - 1, N[ucell])
    m, gap, o = _plateau(base, upos, udelta, target, n)
    if gap <= PB_TOL:
        out = thin.copy()
        np.add.at(out, ucell[o[:m]], 1)
        return out, "topup", design
    # the thinning scalar, the top-up at its own rate (a = 1)
    topped = c.copy()
    np.add.at(topped, ucell[upos <= 1], 1)
    base = float(_val(topped, N).sum())
    dc, dr, dp = dcell[dropped], drank[dropped], dpos[dropped]
    ddelta = _val(c[dc] - dr, N[dc]) - _val(c[dc] - dr + 1, N[dc])
    m, gap, o = _plateau(base, dp, ddelta, target, n)
    if gap <= PB_TOL:
        out = topped.copy()
        np.subtract.at(out, dc[o[:m]], 1)
        return out, "thin", design
    return c.copy(), "kept", design


@dataclass
class Blurred:
    X: object  # the blurred X, the input's type and dtype
    record: dict


def blur(X, coords, s, seed: int, mode: str = BLOCKS, permute: bool = False) -> Blurred:
    """The world blur of one file's X (log1p CP10k) at its coordinates: `mode` blocks (WBH, W, P) or cells (U);
    `permute` gives arm P's spread. Raises CountRecoveryError or HoldError (a refusal)."""
    if mode not in MODES:
        raise ValueError(f"mode is one of {MODES}")
    dense_in = _dense(X)
    x32 = dense_in.astype(np.float32, copy=False)
    n, ng = x32.shape
    s = np.asarray(s, dtype=np.float64)
    if s.shape != (ng,) or not np.isfinite(s).all() or (s < 0).any():
        raise ValueError(f"the spread must be {ng} finite values at least 0")
    if permute:
        s = permuted(s, seed)
    c, N, used = recover_counts(x32)
    T = c.sum(1)
    worlds = block_worlds(coords, seed, ng) if mode == BLOCKS else cell_worlds(n, seed, ng)
    F = worlds.factors(s)
    moved = (worlds.unit >= 0) & (T > 0) & ~np.all(F == 1, axis=1)
    pbar = neighbour_proportions(x32.astype(np.float64), c, T)
    lam = 0.5 * (c + T[:, None] * pbar)
    k = depth_scalars(c, F, lam, T, np.flatnonzero(moved))
    FP = np.where(moved[:, None], k[:, None] * F, 1.0)
    frozen = ~moved
    for _ in range(6):
        out, design, hows = np.empty_like(c), np.empty_like(c), []
        for g in range(ng):
            out[:, g], how, design[:, g] = hold_gene(c[:, g], FP[:, g], lam[:, g], np.where(N > 0, N, 1.0), ~frozen,
                                                     np.random.default_rng([seed, STREAM_COUNTS, g]))
            hows.append(how)
        empty = (out.sum(1) == 0) & (T > 0) & ~frozen
        if not empty.any():
            break
        frozen |= empty
    else:
        raise HoldError("cells kept emptying after six rounds")
    changed = out != c
    y = x32.copy()
    Nc = np.broadcast_to(N[:, None], c.shape)[changed]
    y[changed] = np.log1p(LIBRARY * out[changed] / Nc).astype(np.float32)
    diff = np.abs(y.astype(np.float64).mean(0) - x32.astype(np.float64).mean(0))
    if diff.max(initial=0.0) > PB_TOL:
        raise HoldError(f"{int((diff > PB_TOL).sum())} genes end more than {PB_TOL} from the input's mean X")
    if sparse.issparse(X):
        Y = sparse.csr_matrix(y, dtype=X.dtype)
        Y.eliminate_zeros()
    else:
        Y = y.astype(np.asarray(X).dtype, copy=False)
    T2 = out.sum(1)
    record = {"mode": mode, "seed": int(seed), "permuted": bool(permute), "cells": int(n),
              "recovery": {"exact": True, "smallest_count": {str(m): int((used == m).sum()) for m in MIN_COUNTS}},
              "worlds": len(worlds.z), "antithetic": worlds.antithetic(), "moved_cells": int(moved.sum()),
              "frozen_cells": int((frozen & moved).sum()),
              "genes": {h: hows.count(h) for h in ("topup", "thin", "kept")},
              "library_ratio_design": float(design.sum(1).mean() / T.mean()) if T.mean() > 0 else float("nan"),
              "library_ratio": float(T2.mean() / T.mean()) if T.mean() > 0 else float("nan"),
              "changed_cells": float(changed.any(1).mean()), "pb_max_abs_diff": float(diff.max(initial=0.0)),
              "mean_depth_scalar": float(k[moved].mean()) if moved.any() else 1.0}
    return Blurred(Y, record)


# ------------------------------------------------------------------------------------------------ checks
def pseudobulk(X) -> np.ndarray:
    return np.asarray(_dense(X), dtype=np.float64).mean(0)


def spearman(a, b) -> float:
    from scipy.stats import spearmanr

    r = spearmanr(np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)).statistic
    return float(r) if np.isfinite(r) else float("nan")


OUTPUT_CHECKS = ("count recovery", "pseudobulk rank", "pseudobulk gap", "dp rank", "library", "changed cells",
                 "antithetic")
PB_RANK_MIN, PB_GAP_MAX, DP_RANK_MIN, LIBRARY_TOL, CHANGED_MIN = 0.99, 0.15, 0.99, 0.02, 0.5


def output_checks(x_in, out: Blurred | None, refs: dict[str, np.ndarray], error: str | None = None) -> dict:
    """G0's checks of one blurred output against its input (fix 4b's dp included): each check's name -> True when it
    holds. `refs` maps a reference stage's name to its pseudobulk. A refused blur (error) fails count recovery or the
    hold and every check after it."""
    if out is None:
        rec = {k: False for k in OUTPUT_CHECKS}
        if error and "recover" not in error.lower():
            rec["count recovery"] = True
        return rec
    pb_c, pb_w = pseudobulk(x_in), pseudobulk(out.X)
    dp = [spearman(pb_w - r, pb_c - r) >= DP_RANK_MIN for r in refs.values()]
    r = out.record
    return {"count recovery": bool(r["recovery"]["exact"]),
            "pseudobulk rank": bool(spearman(pb_w, pb_c) >= PB_RANK_MIN),
            "pseudobulk gap": bool(np.abs(pb_w - pb_c).max(initial=0.0) <= PB_GAP_MAX),
            "dp rank": bool(dp and all(dp)),
            "library": bool(abs(r["library_ratio_design"] - 1) <= LIBRARY_TOL),
            "changed cells": bool(r["changed_cells"] >= CHANGED_MIN),
            "antithetic": bool(r["antithetic"])}


def g00(x, s, seed: int) -> dict:
    """Fix 4a's kernel-scale check on a released stage's X (log1p CP10k), target-free: c, the mean over the scorer's
    scales of 2 scale gamma0 E||Pi30 Delta||^2, Delta one world draw's shift of the log1p of the expressed entries
    (G00_DRAWS draws on stream [seed, 66]); the cut (1 - e^(-c/2)) / 2; PASS when the cut is at least 0.10."""
    from sklearn.decomposition import PCA

    x = _dense(x).astype(np.float64)
    rng = np.random.default_rng(0)
    sub = x[rng.choice(x.shape[0], min(G00_CELLS, x.shape[0]), replace=False)]
    pca = PCA(n_components=min(G00_PCS, sub.shape[1]), random_state=0).fit(sub)
    B = pca.transform(sub)
    sq = (B ** 2).sum(1)
    d2 = np.maximum(sq[:, None] + sq[None, :] - 2 * B @ B.T, 0.0)
    gamma0 = 1.0 / (float(np.median(d2[d2 > 0])) + 1e-9)
    s = np.asarray(s, dtype=np.float64)
    z = np.random.default_rng([seed, STREAM_G00]).standard_normal((G00_DRAWS, len(s)))
    v = np.expm1(sub)
    e = 0.0
    for zi in z:
        f = np.exp(s * zi - 0.5 * s ** 2)
        delta = np.where(v > 0, np.log1p(v * f[None, :]) - sub, 0.0)
        e += float(((delta @ pca.components_.T) ** 2).sum(1).mean())
    e /= G00_DRAWS
    c = float(np.mean([2 * sc * gamma0 * e for sc in G00_SCALES]))
    cut = (1 - np.exp(-c / 2)) / 2
    return {"c": c, "cut": float(cut), "gamma0": gamma0, "shift_sq": e, "pass": bool(cut >= G00_MIN_CUT)}


PRIOR_CHECKS = ("pairs", "median s", "ko against replicate")
MIN_USED_PAIRS, S_LO, S_HI, KO_REPLICATE_MAX = 3, 0.10, 1.00, 2.0


def prior_checks(prior: Prior) -> dict:
    """G0's checks of the prior: at least 3 pairs with 2 or more groups; median s in [0.10, 1.00]; the KO-against-WT
    pairs' median spread at most twice the replicate pair's (fix 4c). The rank-correlation check is dropped (fix 3)."""
    sm = prior.summary()
    ko, rep = sm["ko_wt_median_spread"], sm["replicate_median_spread"]
    return {"pairs": len(prior.used) >= MIN_USED_PAIRS,
            "median s": bool(np.isfinite(sm["median_s"]) and S_LO <= sm["median_s"] <= S_HI),
            "ko against replicate": bool(np.isfinite(ko) and np.isfinite(rep) and ko <= KO_REPLICATE_MAX * rep)}
