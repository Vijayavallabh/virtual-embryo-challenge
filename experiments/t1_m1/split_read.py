#!/usr/bin/env python
"""Try T51 (experiments/harness/TRIES.md): the released day's share along the part of Feng's direction that all cell
groups share, and along the rest. A labelled diagnostic on try T47's route; it scores nothing and calls no scorer.

Try T47 (share_read.py) measured the share A of the released day E9.5 -> E10.5 along Feng's E9.5 -> E14.5 direction d
and a scale L (between-type contrasts at E9.5, challenge against Feng). L comes from contrasts inside one stage, in
which anything that differs between whole samples cancels; d contrasts samples five days apart, in which it does not.
The part of d most exposed to such a sample difference is the part all cell groups share. This read measures the
released day's share along that common part and along the rest.

Definitions (T47's wherever one exists; share_read.py's functions are imported, not copied):
  groups, masks, cells, change, weights, bootstrap   T47's (share_read: group_labels_by_stage, stage_table,
           block_sums, resample_weights, cell_set, share_stats, combine, summarise; dirlib.lfc_ratio).
  M        the genes in every group's mask (Feng's direction finite in every group, the gene not held).
  c        per gene of M, the unweighted mean over every direction group of d_g; s_g = d_g - c on M. Beside it,
           reported only: c_w, the mean weighted by the groups' cells in E10.5's library 3, and s_g,w = d_g - c_w; c_n,
           the unweighted mean over one group per overlapping Feng label set (groups in sorted order, a group kept
           unless its Feng labels meet those of a group already kept: Endothelium pools endocardial_ec with
           vascular_ec, so it and Endocardium do not both enter), and s_g,n = d_g - c_n.
  energy   per group over M: e_c = sum(c^2) / sum(d_g^2), e_s = sum(s_g^2) / sum(d_g^2), cross = 2 sum(c s_g) /
           sum(d_g^2); combined with T47's weights over the read set.
  share    of a change ch along v over M: sum(ch v) / sum(v v) (share_read.share_stats' share). Read along d_g on M
           ('M'), c, s_g, c_w, s_g,w, c_n and s_g,n for every comparison of T47; T47's own share along d_g on its own
           mask beside it.
  r_v      from Feng's two halves (feng_dirs.direction_halves, called per group as T47 does) with T47's formula
           (share_read.contrast_reliability) over the genes of M where v, v_a and v_b are finite; c_a and c_b are the
           unweighted means of the groups' d_a and d_b, s_g,a = d_g,a - c_a (likewise b). A corrected share is the
           share over r_v where r_v > 0 (share_read.share_stats' share_corrected), per group, then combined.
  read set the groups of T47's committed A set (t47.json's A_K_support) with a raw and a corrected share along d_g on
           M, c and s_g in libraries 1, 2 and 3 of (b) at the point. A_v, A_v,corr, K_v, R and R_corr all rest on this
           one set (fixed at the point for every draw; a draw in which a group of the set has no value reads undefined),
           so the common and the specific part are never combined over different groups.
  A_v, K_v T47's A and K along v: A = the library-3-weighted combined share of (b) libraries 1 and 2, averaged; K =
           library 3's combined share - A; over the read set. A_M is A along d_g on M.
  R        A_c / A_s; R_corr = A_c,corr / A_s,corr; each NaN where its denominator is not above 0.
  calibration  T47's calibration arm unchanged (t1_scope.external_shift, fraction 0.2, stream (47, 1), on E10.5
           library 3's cells), its before -> after change read along d_g on M, c and s_g (raw), combined over the read
           set (T47's own combined share along d_g on its masks beside it, over every group with a value, as T47).
  self-check   T47's self-check world unchanged (E9.5's cells split on stream (47, 3), H1 moved along d by 0.2 on stream
           (47, 2), H2 the world's E9.5), read along M, c and s_g at the point over the groups of the read set with a
           value in libraries 1 and 2 there (the set it reports is the set its values rest on); as in T47 the world's
           direction is d itself, so every r_v is 1 there and the corrected values equal the raw ones. Known: A_M = A_c
           = A_s = 0.2, R = R_corr = 1 (the brief's); beside it the known values of the definitions themselves,
           0.2 <d_g, v> / <v, v> combined (they differ from 0.2 by the cross terms), reported only. Bias = measured
           minus the brief's known. Because every r_v is 1 there, the self-check does not test the noise correction.
  correction check (reported, read by nothing)   the same world with H1 moved along half A's direction d_g,a by 0.2
           (stream (47, 2)), read along the vectors built from half B alone (d_g,b, c_b, s_g,b), each corrected by its
           single-half reliability max(0, mean(v_b^2) - 2 mean(((v_a - v_b) / 2)^2)) / mean(v_b^2) (one half carries
           twice the noise variance of their mean); half A's noise is independent of half B's, so the corrected shares
           should read 0.2 and R_corr 1 up to the cross terms: its bias is the correction's own error, at twice the
           noise of the read. Beyond it the correction's error is known only from synthetic data
           (tests/test_t1_split_read.py): through the whole read at 119 genes, R_corr within 0.13 of the truth (seeds
           10-19, r_s 0.4-0.9); on vectors alone (four groups, unequal weights, Student-t(3) genes, half noise sd 0.5,
           truth 0.5, seeds 0-19), largest error 0.116 at 119 genes and 0.017 at 2000 genes (r_s 0.61 and above). A
           group whose r_s is near 0 is not covered: its corrected share is its share over r_s, unbounded, and the
           bootstrap (r fixed) cannot show it; r_s <= 0 voids the read, a small positive r_s shows only in the
           leave-one-out values (e.g. r_s 0.033 gave R_corr 0.833 against a truth of 0.5, 0.495 without that group).
  Feng half gauge (reported, read by nothing)   at the point, R and R_corr with every vector built from half A alone
           and from half B alone (each corrected by its single-half reliability), and their difference: Feng's own
           noise in the ratio, which T47's bootstrap (Feng fixed) does not hold.
  leave one out (reported, read by nothing)   R_corr with each group of the read set left out in turn, with its
           bootstrap interval: whether one group (a small r_s) carries R_corr.
  L_eff    L and L_corr from the committed t47_share_read/t47.json (not recomputed); L_eff = L (e_s + min(1,
           R_corr) e_c) / (e_s + e_c), L_eff_corr the same from L_corr; the courses of T47 (share_read.courses):
           f_even = (L_eff - A_M) / 2 and f_decay from (A_M, L_eff), and from (A_M,corr, L_eff_corr). L_eff is
           conditional on the sample-difference reading (under the pace reading L is unchanged and the course's shape
           differs instead) and one-sided by construction (min(1, .) can only lower it); beside it, reported only,
           L_eff_uncapped with R_corr in place of min(1, R_corr). The energy parts come from the noisy d, whose noise
           inflates e_s more than e_c and so pulls L_eff toward L.
Every value has T47's bootstrap over the challenge's cells (200 draws, stream (51, 0); Feng's statistics, the
reliabilities, the energy parts, L, the calibration, the self-check, the correction check and the gauge fixed), its 5th
and 95th percentiles and its number of undefined draws.

The reading (fixed in row T51 before the run):
  void unless A_M is within 0.02 of T47's committed A, A on T47's own masks equals T47's committed A within 1e-9 (the
  inputs are T47's), every group of T47's A set is in the read set (a raw and a corrected share along M, c and s_g in
  libraries 1, 2 and 3 of (b)), every group of T47's A set has r_M, r_c and r_s above 0, and the self-check's R is
  defined.
  On R_corr with its bootstrap interval widened on both sides by the absolute self-check bias of R_corr: "the common
  part moves less" if the widened upper end is below 0.8; "the common part moves as much" if its lower end is above 0.9;
  otherwise "not separated". The same on R (raw, widened by R's bias); where the two differ the reading is "not
  separated" and the text says which of the two it turned on. Raw R is not an independent check: the specific part's
  share is the more attenuated one, so R sits above R_corr by about r_c / r_s; where R reads the higher statement and
  sits above R_corr, the text says that this "not separated" is that attenuation bias of raw R, not cell-sampling
  uncertainty, and prints R / R_corr and r_c / r_s beside it ("less" on both needs a true ratio below about 0.8 r_s /
  r_c). "not computable" where R_corr is undefined at the point or in any draw (with the number of draws).
  Reported and read by nothing: library 3 (K's parts), (a) E8.5 -> E9.5, per group, everything with c_w and with c_n
  for c, the leave-one-out R_corr, the Feng half gauge, the correction check and L_eff_uncapped.

Ledger sentences (the limits; printed in t51.txt as notes):
- A lower share along the common part can mean that part of d is a sample difference, or that what all cell types share
  runs at another pace than what separates them; the read cannot tell these apart.
- The read does not test the sample-difference hypothesis; it tests only whether the common part's share is lower, and
  no control in it (the (a) E8.5 -> E9.5 ratio included) separates a Feng sample difference from another pace.
- The challenge's own sample differences between E9.5 and E10.5 also lie along the common part and can raise or lower
  A_c.
- The specific part is a difference of noisy directions and its share is the more attenuated one, which the correction
  addresses only as far as Feng's two halves measure the noise.
- The bootstrap intervals hold only the challenge's cell sampling: Feng's noise, the error of the reliabilities and the
  choice of group set are outside them (the gauge, the correction check and the leave-one-out values are beside them).
- L_eff assumes the specific part's scale is L, which between-type contrasts support, and gives the common part R_corr
  of it; it holds only under the sample-difference reading and can only lower L.
- No statement says how far the hidden test's days move.
- c and s_g are not orthogonal within a group: a share along each is a projection, and a change made exactly along d_g
  reads 1 + <s_g, c> / <c, c> times its share along c and 1 + <c, s_g> / <s_g, s_g> times it along s_g; the self-check
  reports these definitional values beside the brief's known R = 1.
- T51 puts nothing forward by itself; its result can only go to the team captain.

  uv run --frozen --no-sync python -I -B experiments/t1_m1/split_read.py run --out t51_split_read

Writes experiments/t1_m1/results/<out>/t51.{json,txt}; refused without row T51 in TRIES.md, on an existing output and
on a dirty tree.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import share_read as sr  # noqa: E402  (it puts the experiment folders T47 imports from on sys.path)

dl, feng_dirs, gv, gc, hygiene, n12 = sr.dl, sr.feng_dirs, sr.gv, sr.gc, sr.hygiene, sr.n12
t1_scope = sr.t1_scope

ROW = "T51"
STEM = "t51"
RESULTS = HERE / "results"
T47_JSON = RESULTS / "t47_share_read" / "t47.json"
N_BOOT = sr.N_BOOT  # 200
BOOT_SEED = (51, 0)  # row T51's stream: no existing step draws on it
A_TOL = 0.02  # A_M within this of T47's committed A, else void
MASK_TOL = 1e-9  # A on T47's own masks equals T47's committed A within this, else void (the inputs are not T47's)
LESS_AT, AS_MUCH_AT = 0.8, 0.9
SELF_KNOWN = 0.2  # the self-check world's share along d (T47's), and so along c and s_g by the brief
VECS = ("M", "c", "s", "c_w", "s_w", "c_n", "s_n")  # read along: d_g on M, c, s_g, c_w, s_g,w, c_n, s_g,n
READ_VECS = ("M", "c", "s")  # the read set needs a raw and a corrected share along each
FIELDS = ("share", "share_corrected")
PAIRS = {"": ("c", "s"), "_w": ("c_w", "s_w"), "_n": ("c_n", "s_n")}  # R from (common, specific); _w, _n reported
LIBS12, LIB3 = sr.AK_COMPARISONS[:2], sr.AK_COMPARISONS[2]
PREV = "a_E8.5_E9.5"
SELF_KEYS = ("A_M", "A_c", "A_s", "R", "R_corr", "A_c_w", "A_s_w", "R_w", "R_w_corr")
GAUGE_KEYS = ("A_M", "A_c", "A_s", "R", "A_M_corr", "A_c_corr", "A_s_corr", "R_corr")

LESS = "the common part moves less"
AS_MUCH = "the common part moves as much"
NOT_SEPARATED = sr.NOT_SEPARATED
NOT_COMPUTABLE = sr.NOT_COMPUTABLE
VOID = "void"
ORDER = {LESS: 0, NOT_SEPARATED: 1, AS_MUCH: 2}
ATTENUATION = ("raw R reads the higher statement and sits above R_corr by the attenuation of the noisier specific "
               "part: this 'not separated' is that bias of raw R, not cell-sampling uncertainty")
LABEL_T51 = ("try T51, a labelled diagnostic on try T47's route, counted, never an adoption result: the released day's "
             "share along the part of Feng's E9.5 -> E14.5 direction that all cell groups share and along the rest; it "
             "compares a predictor with the cuts' target stage, scores nothing and calls no scorer; the read puts "
             "nothing forward by itself and its result can only go to the team captain")
LEDGER_NOTES = sr._ledger_notes(__doc__)


# ---------------------------------------------------------------------------------------------------------------------
# The split of d, its energy parts and reliabilities


def common_mask(masks: Mapping[str, np.ndarray]) -> np.ndarray:
    """M: the genes in every group's mask."""
    return np.logical_and.reduce([np.asarray(m, bool) for m in masks.values()])


def non_overlapping_groups(groups: Sequence[str], feng_labels: Mapping[str, Sequence[str] | None]) -> list[str]:
    """One group per overlapping Feng label set: groups in sorted order, each kept unless its Feng labels meet those of
    a group already kept (a group without Feng labels is kept)."""
    kept, used = [], set()
    for g in sorted(groups):
        labs = set(feng_labels.get(g) or [])
        if labs & used:
            continue
        kept.append(g)
        used |= labs
    return kept


def split_vectors(dirs: Mapping[str, np.ndarray], M: np.ndarray, weights: Mapping[str, float] | None = None,
                  among: Sequence[str] | None = None) -> dict:
    """c (the unweighted mean of the groups' d_g on M, NaN off M) and s_g = d_g - c; with weights, c_w (the mean
    weighted by the positive weights; NaN everywhere where they sum to 0) and s_g,w = d_g - c_w; with among, c_n (the
    unweighted mean over the groups of `among`; NaN where none) and s_g,n = d_g - c_n for every group."""
    groups = sorted(dirs)
    D = np.vstack([np.asarray(dirs[g], np.float64) for g in groups])
    c = np.full(D.shape[1], np.nan)
    c[M] = D[:, M].mean(0)
    out = {"c": c, "s": {g: np.where(M, D[i] - c, np.nan) for i, g in enumerate(groups)}}
    if weights is not None:
        w = np.array([max(float(weights.get(g, 0.0)), 0.0) for g in groups])
        cw = np.full(D.shape[1], np.nan)
        if w.sum() > 0:
            cw[M] = (w[:, None] * D[:, M]).sum(0) / w.sum()
        out["c_w"] = cw
        out["s_w"] = {g: np.where(M, D[i] - cw, np.nan) for i, g in enumerate(groups)}
    if among is not None:
        idx = [i for i, g in enumerate(groups) if g in set(among)]
        cn = np.full(D.shape[1], np.nan)
        if idx:
            cn[M] = D[idx][:, M].mean(0)
        out["c_n"] = cn
        out["s_n"] = {g: np.where(M, D[i] - cn, np.nan) for i, g in enumerate(groups)}
    return out


def vector_set(dirs: Mapping[str, np.ndarray], sv: Mapping) -> dict:
    """The read vectors per group from split_vectors' output: d_g ('M'), c, s_g and whichever of c_w, c_n it holds."""
    groups = sorted(dirs)
    out = {"M": {g: dirs[g] for g in groups}}
    for vc, vs in PAIRS.values():
        if vc in sv:
            out[vc], out[vs] = dict.fromkeys(groups, sv[vc]), dict(sv[vs])
    return out


def energy_parts(d: np.ndarray, c: np.ndarray, s: np.ndarray, M: np.ndarray) -> dict:
    """e_c = sum(c^2) / sum(d^2), e_s = sum(s^2) / sum(d^2), cross = 2 sum(c s) / sum(d^2), over M (NaN where sum(d^2)
    = 0); the three add up to 1 where s = d - c."""
    dd = float(np.sum(np.asarray(d)[M] ** 2))
    if dd == 0 or not np.isfinite(dd):
        return dict.fromkeys(("e_c", "e_s", "cross"), float("nan"))
    cm, sm = np.asarray(c)[M], np.asarray(s)[M]
    return {"e_c": float(cm @ cm) / dd, "e_s": float(sm @ sm) / dd, "cross": 2.0 * float(cm @ sm) / dd}


def half_vectors(halves: Mapping[str, Mapping[str, np.ndarray]], groups: Sequence[str], M: np.ndarray,
                 weights: Mapping[str, float], among: Sequence[str] | None = None) -> dict:
    """The half versions of every read vector: per half h in (a, b), d_g,h, c_h, s_g,h = d_g,h - c_h and the same with
    c_w (and c_n with among); a vector whose half needs a group without halves is all NaN."""
    nan = np.full(len(M), np.nan)
    out: dict = {}
    complete = all(g in halves for g in groups)
    for h, key in (("a", "d_a"), ("b", "d_b")):
        dh = {g: np.asarray(halves[g][key], np.float64) if g in halves else nan for g in groups}
        if complete:
            sv = split_vectors(dh, M, weights, among)
        else:
            sv = {"c": nan, "s": dict.fromkeys(groups, nan), "c_w": nan, "s_w": dict.fromkeys(groups, nan)}
            if among is not None:
                sv.update({"c_n": nan, "s_n": dict.fromkeys(groups, nan)})
        out[h] = vector_set(dh, sv)
    return out


def reliabilities(vecs: Mapping[str, Mapping[str, np.ndarray]], hv: Mapping, M: np.ndarray) -> dict:
    """r_v per read vector and group: T47's formula (share_read.contrast_reliability) over the genes of M where v, v_a
    and v_b are finite."""
    return {v: {g: sr.contrast_reliability(vecs[v][g], hv["a"][v][g], hv["b"][v][g], M) for g in vecs[v]}
            for v in vecs}


def half_reliability(v_h: np.ndarray, v_a: np.ndarray, v_b: np.ndarray, M: np.ndarray) -> float:
    """The reliability of one half's vector v_h: max(0, mean(v_h^2) - 2 mean(((v_a - v_b) / 2)^2)) / mean(v_h^2) over
    the genes of M where the three are finite (one half carries twice the noise variance of the halves' mean); NaN
    where no gene or mean(v_h^2) = 0."""
    v_h, v_a, v_b = (np.asarray(x, np.float64) for x in (v_h, v_a, v_b))
    ok = np.asarray(M, bool) & np.isfinite(v_h) & np.isfinite(v_a) & np.isfinite(v_b)
    if not ok.any():
        return float("nan")
    msq = float(np.mean(v_h[ok] ** 2))
    if msq == 0:
        return float("nan")
    return max(msq - 2.0 * float(np.mean(((v_a[ok] - v_b[ok]) / 2.0) ** 2)), 0.0) / msq


def half_reliabilities(hv: Mapping, M: np.ndarray) -> dict:
    """Per half h, read vector and group: half_reliability of the half's own vector."""
    return {h: {v: {g: half_reliability(hv[h][v][g], hv["a"][v][g], hv["b"][v][g], M) for g in hv[h][v]}
                for v in hv[h]} for h in ("a", "b")}


# ---------------------------------------------------------------------------------------------------------------------
# One measurement


def _ratio(a: float, b: float) -> float:
    return float(a / b) if np.isfinite(a) and np.isfinite(b) and b > 0 else float("nan")


def l_eff(L: float, R_corr: float, e_s: float, e_c: float, cap: bool = True) -> float:
    """L (e_s + min(1, R_corr) e_c) / (e_s + e_c) (cap=False: R_corr in place of min(1, R_corr)); NaN where any input
    is not finite or e_s + e_c = 0."""
    vals = (L, R_corr, e_s, e_c)
    if not all(np.isfinite(x) for x in vals) or e_s + e_c == 0:
        return float("nan")
    r = min(1.0, R_corr) if cap else R_corr
    return float(L * (e_s + r * e_c) / (e_s + e_c))


def change_shares(ch: np.ndarray, g: str, setup: Mapping) -> dict:
    """One change's shares along each read vector over M (raw and corrected by its r_v) and along d_g on T47's mask."""
    out = {"d_T47/share": sr.share_stats(ch, setup["dirs"][g], setup["masks"][g], None)["share"]}
    for v in VECS:
        st = sr.share_stats(ch, setup["vecs"][v][g], setup["M"], setup["rel"][v][g])
        out[f"{v}/share"], out[f"{v}/share_corrected"] = st["share"], st["share_corrected"]
    return out


def strict_combine(values: Mapping[str, float], among: Sequence[str], weights: Mapping[str, float]) -> float:
    """share_read.combine's weighted mean over the groups of `among`; NaN where `among` is empty or any of its groups
    has no value (the set is never shrunk silently)."""
    if not among or not all(np.isfinite(values[g]) for g in among):
        return float("nan")
    return sr.combine({g: values[g] for g in among}, weights)[0]


def read_set(out: Mapping[str, float], candidates: Sequence[str], comparisons: Sequence[str] = sr.AK_COMPARISONS,
             vecs: Sequence[str] = READ_VECS) -> list[str]:
    """The groups of `candidates` with a raw and a corrected share along every vector of `vecs` in every comparison."""
    return [g for g in candidates if not missing_shares(out, [g], comparisons, vecs)]


def missing_shares(out: Mapping[str, float], candidates: Sequence[str], comparisons: Sequence[str] = sr.AK_COMPARISONS,
                   vecs: Sequence[str] = READ_VECS) -> dict[str, list[str]]:
    """Per group of `candidates` lacking one: the comparison/vector/field shares that are not finite."""
    miss = {}
    for g in candidates:
        bad = [f"{n}/{v}/{f}" for n in comparisons for v in vecs for f in FIELDS
               if not np.isfinite(out.get(f"{n}/{g}/{v}/{f}", float("nan")))]
        if bad:
            miss[g] = bad
    return miss


def split_values(out: Mapping[str, float], setup: Mapping) -> dict:
    """A, K, the library-3 and (a) shares along every read vector over the read set (setup['read_set']; (a) over
    setup['p_set']), per group, R and its parts, R_corr with each group left out, L_eff and the courses, from one
    measurement's per-group shares."""
    w, read, pset = setup["weights"], list(setup["read_set"]), list(setup["p_set"])

    def comb(name, v, field, among=read):
        return strict_combine({g: out[f"{name}/{g}/{v}/{field}"] for g in among}, among, w)

    rv: dict[str, float] = {}
    for v in VECS:
        for field, tag in (("share", ""), ("share_corrected", "_corr")):
            A = (comb(LIBS12[0], v, field) + comb(LIBS12[1], v, field)) / 2.0
            rv[f"A_{v}{tag}"] = A
            rv[f"lib3_{v}{tag}"] = comb(LIB3, v, field)
            rv[f"K_{v}{tag}"] = rv[f"lib3_{v}{tag}"] - A
            rv[f"P_{v}{tag}"] = comb(PREV, v, field, pset)
    rv["A_M_T47_mask"] = (comb(LIBS12[0], "d_T47", "share") + comb(LIBS12[1], "d_T47", "share")) / 2.0
    for g in read:  # per group, A's two libraries averaged
        for v in READ_VECS:
            for field, tag in (("share", ""), ("share_corrected", "_corr")):
                rv[f"group/{g}/A_{v}{tag}"] = (out[f"{LIBS12[0]}/{g}/{v}/{field}"]
                                               + out[f"{LIBS12[1]}/{g}/{v}/{field}"]) / 2.0
    for wtag, (vc, vs) in PAIRS.items():
        for tag in ("", "_corr"):
            rv[f"R{wtag}{tag}"] = _ratio(rv[f"A_{vc}{tag}"], rv[f"A_{vs}{tag}"])
            rv[f"R_lib3{wtag}{tag}"] = _ratio(rv[f"lib3_{vc}{tag}"], rv[f"lib3_{vs}{tag}"])
            rv[f"R_K{wtag}{tag}"] = _ratio(rv[f"K_{vc}{tag}"], rv[f"K_{vs}{tag}"])
            rv[f"R_a{wtag}{tag}"] = _ratio(rv[f"P_{vc}{tag}"], rv[f"P_{vs}{tag}"])
    if len(read) > 1:
        for g in read:  # leave one out: reported only
            rest = [h for h in read if h != g]
            a_c = (comb(LIBS12[0], "c", "share_corrected", rest) + comb(LIBS12[1], "c", "share_corrected", rest)) / 2
            a_s = (comb(LIBS12[0], "s", "share_corrected", rest) + comb(LIBS12[1], "s", "share_corrected", rest)) / 2
            rv[f"R_corr_without/{g}"] = _ratio(a_c, a_s)
    for wtag in PAIRS:
        E = setup["energy"][wtag or "unweighted"]
        R_corr = rv[f"R{wtag}_corr"]
        for ltag, L, A in (("", setup["L"], rv["A_M"]), ("_corr", setup["L_corr"], rv["A_M_corr"])):
            le = l_eff(L, R_corr, E["e_s"], E["e_c"])
            rv[f"L_eff{wtag}{ltag}"] = le
            rv[f"L_eff_uncapped{wtag}{ltag}"] = l_eff(L, R_corr, E["e_s"], E["e_c"], cap=False)
            cc = sr.courses(A, le)
            rv[f"A_M_over_L_eff{wtag}{ltag}"] = cc["A_over_L"]
            for k in sr.COURSE_KEYS:
                rv[f"{k}{wtag}{ltag}"] = cc[k]
    return rv


def _mean_of(tabs, sums, setup, spec, g):
    stage, libs = spec
    labs = setup["labels"].get(stage, {}).get(g)
    if not labs or stage not in tabs:
        return None
    return sr.cell_set(tabs[stage], sums[stage], labs, libs, setup["min_cells"])


def measure_shares(tabs: Mapping[str, dict], sums: Mapping[str, np.ndarray], setup: Mapping) -> dict:
    """Per comparison and group, the shares along each read vector from one set of block sums."""
    out: dict[str, float] = {}
    keys = ["d_T47/share", *(f"{v}/{f}" for v in VECS for f in FIELDS)]
    for name, (sa, sb) in sr.COMPARISONS.items():
        for g in setup["groups"]:
            a, b = _mean_of(tabs, sums, setup, sa, g), _mean_of(tabs, sums, setup, sb, g)
            st = (dict.fromkeys(keys, float("nan")) if a is None or b is None
                  else change_shares(dl.lfc_ratio(a, b, setup["eps"]), g, setup))
            for k, v in st.items():
                out[f"{name}/{g}/{k}"] = v
    return out


def measure_split(tabs: Mapping[str, dict], sums: Mapping[str, np.ndarray], setup: Mapping) -> dict:
    """Every scalar of the read from one set of block sums: per comparison and group the shares along each read vector,
    and the reading's values (over setup's read set) under 'reading/'."""
    out = measure_shares(tabs, sums, setup)
    for k, v in split_values(out, setup).items():
        out[f"reading/{k}"] = v
    return out


def finish_setup(setup: dict, shares: Mapping[str, float]) -> dict:
    """From the point's per-group shares: the read set (T47's A set with every share along M, c and s_g in libraries
    1, 2 and 3), the (a) set (every group with every share along M, c and s_g in (a)) and the energy parts combined
    over the read set."""
    setup["read_set"] = read_set(shares, setup["support"])
    setup["p_set"] = read_set(shares, setup["groups"], (PREV,))
    setup["energy"] = {k: {p: sr.combine({g: e[g][p] for g in setup["read_set"]}, setup["weights"])[0]
                           for p in ("e_c", "e_s", "cross")} for k, e in setup["energy_g"].items()}
    return setup


# ---------------------------------------------------------------------------------------------------------------------
# The statement


def rule_statement(lo: float, hi: float) -> str:
    """'less' if the (widened) upper end is below 0.8, 'as much' if the lower end is above 0.9, else 'not separated'."""
    if hi < LESS_AT:
        return LESS
    if lo > AS_MUCH_AT:
        return AS_MUCH
    return NOT_SEPARATED


def side_statement(rec: Mapping, bias: float | None, n_boot: int, name: str) -> str:
    """One ratio's statement on its bootstrap interval widened by |bias|: 'not computable' where it is undefined at the
    point or in a draw (share_read.undefined) or where its bias is undefined."""
    u = sr.undefined([rec], n_boot)
    if u is not None:
        return u
    if bias is None or not np.isfinite(bias):
        return f"{NOT_COMPUTABLE} (the self-check's bias of {name} is undefined)"
    w = sr.widen(rec, bias)
    return rule_statement(w["p5"], w["p95"])


def statement(r: Mapping[str, Mapping], bias: Mapping[str, float] | None, n_boot: int,
              void: Sequence[str] = (), r_c_over_r_s: float | None = None) -> dict:
    """The reading: void where any void condition holds; 'not computable' where R_corr is; the shared statement of
    R_corr and R where they agree, else 'not separated' with the side it turned on. 'why_not_separated' says when the
    disagreement is the attenuation bias of raw R (R reads the higher statement and sits above R_corr) and gives R /
    R_corr and r_c / r_s beside it."""
    bias = bias or {}
    corr = side_statement(r["R_corr"], bias.get("R_corr"), n_boot, "R_corr")
    raw = side_statement(r["R"], bias.get("R"), n_boot, "R")
    out = {"R_corr": corr, "R": raw, "turned_on": None, "void": list(void), "why_not_separated": None}
    if void:
        out["statement"] = f"{VOID} ({'; '.join(void)})"
    elif corr.startswith(NOT_COMPUTABLE):
        out["statement"] = corr
    elif raw == corr:
        out["statement"] = corr
        if corr == NOT_SEPARATED:
            out["why_not_separated"] = "the widened interval of R_corr (and of R) does not separate"
    else:
        out["statement"] = NOT_SEPARATED
        if raw.startswith(NOT_COMPUTABLE):
            out["turned_on"] = f"R (raw): {raw}"
        elif corr == NOT_SEPARATED:
            out["turned_on"] = f"R_corr (R_corr: {corr}; R: {raw})"
        elif raw == NOT_SEPARATED:
            out["turned_on"] = f"R (raw) (R_corr: {corr}; R: {raw})"
        else:
            out["turned_on"] = f"R_corr and R read opposite statements (R_corr: {corr}; R: {raw})"
        q = _ratio(r["R"]["value"], r["R_corr"]["value"])
        if raw in ORDER and ORDER[raw] > ORDER[corr] and np.isfinite(q) and q > 1:
            rs = f"{r_c_over_r_s:.3f}" if r_c_over_r_s is not None and np.isfinite(r_c_over_r_s) else "undefined"
            out["why_not_separated"] = f"{ATTENUATION} (R / R_corr {q:.3f} at the point; r_c / r_s {rs})"
        else:
            out["why_not_separated"] = "the readings on R and R_corr differ"
    return out


def void_reasons(point: Mapping[str, float], setup: Mapping, t47_A: float, self_check: Mapping | None) -> list[str]:
    """The void conditions, each as the sentence the text prints."""
    out = []
    a_m = point["reading/A_M"]
    if not (np.isfinite(a_m) and np.isfinite(t47_A) and abs(a_m - t47_A) <= A_TOL):
        out.append(f"A_M {a_m:+.4f} is not within {A_TOL} of T47's committed A {t47_A:+.4f}")
    a_t = point["reading/A_M_T47_mask"]
    if not (np.isfinite(a_t) and np.isfinite(t47_A) and abs(a_t - t47_A) <= MASK_TOL):
        out.append(f"A on T47's own masks {a_t:+.10f} is not T47's committed A {t47_A:+.10f} within {MASK_TOL:g} "
                   "(the inputs are not T47's)")
    sup = list(setup["support"])
    if not sup:
        out.append("T47's A set is empty")
    missing = missing_shares(point, sup)
    if missing:
        out.append("groups of T47's A set without a raw and a corrected share along M, c and s_g in libraries 1, 2 "
                   f"and 3: {sorted(missing)} (first missing: "
                   + "; ".join(f"{g} {missing[g][0]}" for g in sorted(missing)) + ")")
    rel = setup["rel"]
    low = {g: {v: rel[v].get(g, float("nan")) for v in READ_VECS} for g in sup
           if not all(np.isfinite(rel[v].get(g, float("nan"))) and rel[v].get(g, float("nan")) > 0
                      for v in READ_VECS)}
    if low:
        out.append("groups of T47's A set whose reliability along M, c or s_g is not above 0: "
                   + "; ".join(f"{g} " + ", ".join(f"r_{v} {x:.4f}" for v, x in low[g].items()) for g in sorted(low)))
    if self_check is None:
        out.append("the self-check did not run")
    elif not np.isfinite(self_check["measured"]["R"]):
        out.append("the self-check's R is undefined")
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The calibration, the self-check and the correction check (T47's arms unchanged, read along the split)


def _full(dirs: Mapping[str, np.ndarray], pos: np.ndarray, n_vars: int) -> dict[str, np.ndarray]:
    full = {}
    for g, v in dirs.items():
        w = np.full(n_vars, np.nan)
        w[pos] = v
        full[g] = w
    return full


def calibration_split(adata: ad.AnnData, genes: Sequence[str], setup: Mapping, held: Sequence[str],
                      label_key: str = "celltype") -> dict:
    """T47's calibration arm (share_read.calibration's shift: external_shift at CAL_FRACTION on CAL_STREAM, library 3's
    cells of the direction groups at E10.5) read along d_g on M, c and s_g (raw; c_w, s_g,w, c_n and s_g,n beside them):
    per group and combined (library-3 weighted) over the read set (T47's mask over every group with a value, as T47),
    and the ratio of the combined shares along c and s_g."""
    labels = setup["labels"]["E10.5"]
    dirs, w, read = setup["dirs"], setup["weights"], setup["read_set"]
    libs, _ = sr.t1_stats.library_tags(adata.obs_names)
    lab = adata.obs[label_key].astype(str).to_numpy()
    group_of = {lb: g for g in sorted(dirs) for lb in labels.get(g) or []}
    cells = np.flatnonzero((np.asarray(libs) == sr.LIB_HEART) & np.isin(lab, list(group_of)))
    pos = adata.var_names.get_indexer(list(genes))
    per = {g: dict.fromkeys(("d_T47/share", *(f"{v}/share" for v in VECS)), float("nan")) for g in dirs}
    if cells.size:
        sub = adata[cells].copy()
        moved, _ = t1_scope.external_shift(sub, _full(dirs, pos, adata.n_vars), sr.CAL_FRACTION, group_of, held,
                                           label_key=label_key, seed=list(sr.CAL_STREAM))
        sl = lab[cells]

        def cp(a, m):
            X = a.X[np.flatnonzero(m)]
            X = X.toarray() if sparse.issparse(X) else np.asarray(X)
            return {"cp": np.expm1(np.asarray(X, np.float64)[:, pos]).mean(0)}

        for g in dirs:
            m = np.isin(sl, list(labels.get(g) or []))
            if m.sum() == 0 or m.sum() < setup["min_cells"]:
                continue
            st = change_shares(dl.lfc_ratio(cp(sub, m), cp(moved, m), setup["eps"]), g, setup)
            per[g] = {k: v for k, v in st.items() if not k.endswith("share_corrected")}
    comb = {v: strict_combine({g: per[g][f"{v}/share"] for g in read}, read, w) for v in VECS}
    comb["d_T47"] = sr.combine({g: per[g]["d_T47/share"] for g in per}, w)[0]
    return {"fraction": sr.CAL_FRACTION, "stream": list(sr.CAL_STREAM), "groups": per, "combined_weighted": comb,
            "combined_over": list(read), "R": _ratio(comb["c"], comb["s"]), "R_w": _ratio(comb["c_w"], comb["s_w"]),
            "R_n": _ratio(comb["c_n"], comb["s_n"])}


def world_read(stage95: ad.AnnData, genes: Sequence[str], setup: Mapping, held: Sequence[str], label_key: str,
               shift: Mapping[str, np.ndarray], vecs: Mapping, rel: Mapping) -> tuple[dict, list[str]]:
    """T47's self-check world (share_read.split_cells on SPLIT_STREAM; H1 moved by external_shift along `shift` at
    CAL_FRACTION on SELF_STREAM; H2 the world's E9.5) read at the point with the vectors `vecs` and reliabilities
    `rel`, over the groups of the read set with a raw and a corrected share along M, c and s_g in libraries 1 and 2
    there: (the reading's values, that set)."""
    labels95 = setup["labels"]["E9.5"]
    parts = sr.split_cells(stage95, label_key)
    h1 = np.sort(np.concatenate([parts["H1a"], parts["H1b"]]))
    pos = stage95.var_names.get_indexer(list(genes))
    group_of = {lb: g for g in sorted(shift) for lb in labels95.get(g) or []}
    sub = stage95[h1].copy()
    moved, _ = t1_scope.external_shift(sub, _full(shift, pos, stage95.n_vars), sr.CAL_FRACTION, group_of, held,
                                       label_key=label_key, seed=list(sr.SELF_STREAM))
    keep = {lab for g in shift for lab in labels95.get(g) or []}
    tabs = {"E9.5": sr.stage_table(stage95[parts["H2"]], genes, keep, label_key),
            "E10.5": sr.stage_table(moved, genes, keep, label_key)}
    world = {**setup, "labels": {"E9.5": labels95, "E10.5": labels95}, "vecs": vecs, "rel": rel}
    out = measure_shares(tabs, {s: sr.block_sums(t) for s, t in tabs.items()}, world)
    sup = read_set(out, setup["read_set"], LIBS12)
    return split_values(out, {**world, "read_set": sup, "p_set": []}), sup


def definitional_known(setup: Mapping, sup: Sequence[str]) -> dict:
    """The self-check's known values of the definitions themselves: 0.2 <d_g, v> / <v, v> combined over `sup`, and
    the ratios of them."""
    M, w, dirs = setup["M"], setup["weights"], setup["dirs"]

    def proj(v):
        vals = {}
        for g in sup:
            x, d = setup["vecs"][v][g][M], np.asarray(dirs[g])[M]
            vals[g] = SELF_KNOWN * float(d @ x) / float(x @ x) if float(x @ x) > 0 else float("nan")
        return strict_combine(vals, sup, w)

    defin = {f"A_{v}": proj(v) for v in VECS}
    for wtag, (vc, vs) in PAIRS.items():
        defin[f"R{wtag}"] = defin[f"R{wtag}_corr"] = _ratio(defin[f"A_{vc}"], defin[f"A_{vs}"])
    return {k: defin[k] for k in SELF_KEYS}


def self_check_split(stage95: ad.AnnData, genes: Sequence[str], setup: Mapping, held: Sequence[str],
                     label_key: str = "celltype") -> dict:
    """T47's self-check world measured at the point with every r_v = 1 (the world's direction is d itself): the
    measured A_M, A_c, A_s, R, R_corr (and with c_w) over the reported set, the brief's known values (0.2 and R = 1),
    the bias (measured minus that), and beside it the definitions' own known values 0.2 <d_g, v> / <v, v> combined over
    the same set."""
    dirs = setup["dirs"]
    vals, sup = world_read(stage95, genes, setup, held, label_key, dirs, setup["vecs"],
                           {v: dict.fromkeys(dirs, 1.0) for v in VECS})
    measured = {k: vals[k] for k in SELF_KEYS}
    known = {k: (1.0 if k.startswith("R") else SELF_KNOWN) for k in SELF_KEYS}
    defin = definitional_known(setup, sup)
    return {"fraction": sr.CAL_FRACTION, "stream": list(sr.SELF_STREAM), "split_stream": list(sr.SPLIT_STREAM),
            "support": sup, "measured": measured, "known": known,
            "bias": {k: measured[k] - known[k] for k in SELF_KEYS},
            "known_by_definition": defin,
            "bias_from_definition": {k: measured[k] - defin[k] for k in SELF_KEYS}}


def correction_check(stage95: ad.AnnData, genes: Sequence[str], setup: Mapping, held: Sequence[str],
                     label_key: str = "celltype") -> dict | None:
    """Reported, read by nothing: the self-check world with H1 moved along half A's d_g,a, read along half B's vectors
    each corrected by its single-half reliability; the measured values, the brief's known values (0.2, R = 1), the
    bias, and the bias from the definitional known values of the full vectors. None where a group lacks halves."""
    halves = setup["halves"]
    if any(g not in halves for g in setup["dirs"]):
        return None
    shift = {g: halves[g]["d_a"] for g in setup["dirs"]}
    vals, sup = world_read(stage95, genes, setup, held, label_key, shift, setup["half_vecs"]["b"],
                           setup["half_rel"]["b"])
    measured = {k: vals[k] for k in SELF_KEYS}
    known = {k: (1.0 if k.startswith("R") else SELF_KNOWN) for k in SELF_KEYS}
    defin = definitional_known(setup, sup)
    return {"fraction": sr.CAL_FRACTION, "stream": list(sr.SELF_STREAM), "support": sup, "measured": measured,
            "known": known, "bias": {k: measured[k] - known[k] for k in SELF_KEYS},
            "bias_from_definition": {k: measured[k] - defin[k] for k in SELF_KEYS}}


def half_gauge(tabs: Mapping[str, dict], sums: Mapping[str, np.ndarray], setup: Mapping) -> dict:
    """Reported, read by nothing: at the point, A and R (raw and corrected) with every vector built from half A alone
    and from half B alone (each corrected by its single-half reliability), over the read set, and their difference."""
    out = {}
    for h in ("a", "b"):
        s2 = {**setup, "vecs": setup["half_vecs"][h], "rel": setup["half_rel"][h]}
        vals = split_values(measure_shares(tabs, sums, s2), s2)
        out[h] = {k: vals[k] for k in GAUGE_KEYS}
    out["a_minus_b"] = {k: out["a"][k] - out["b"][k] for k in GAUGE_KEYS}
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The read


def t47_record(path: Path = T47_JSON) -> dict:
    """From T47's committed record: A, L, L_corr (the point values) and its A set."""
    rep = json.loads(Path(path).read_text())
    r = rep["read"]
    return {"A": r["reading"]["A"]["value"], "L": r["reading"]["L"]["value"],
            "L_corr": r["reading"]["L_corr"]["value"], "support": list(r["support"]["A_K_support"])}


def build_setup(stages: Mapping[str, ad.AnnData], feng: dict, held: Sequence[str], cfg: Mapping,
                t47: Mapping) -> tuple[dict, dict]:
    """T47's groups, masks, labels, tables, weights (as share_read.reading builds them), and the split: M, the read
    vectors, their reliabilities from Feng's two halves, the half vectors and their single-half reliabilities, and the
    energy parts per group (combined over the read set by finish_setup)."""
    label_key = cfg.get("label_key", "celltype")
    genes = list(stages[sr.GENE_STAGE].var_names)
    dirs, _, _ = feng_dirs.directions_from_stats(feng, genes, cfg, sr.BRACKET)
    groups = sorted(dirs)
    dirs = {g: dirs[g] for g in groups}
    held_set = set(map(str, held))
    not_held = np.array([g not in held_set for g in genes], dtype=bool)
    masks = {g: np.isfinite(dirs[g]) & not_held for g in groups}
    labels, source = sr.group_labels_by_stage(stages, groups, cfg, label_key)
    keep = {s: {lab for g in groups for lab in labels[s][g]} for s in stages}
    tabs = {s: sr.stage_table(stages[s], genes, keep[s], label_key) for s in stages}
    heart = tabs.get("E10.5")
    weights = {}
    for g in groups:  # share_read.reading's weights: the group's cells in E10.5's library 3
        labs = set(labels.get("E10.5", {}).get(g) or [])
        weights[g] = (float(sum(n for (lab, li), n in zip(heart["keys"], heart["n_block"], strict=True)
                                if lab in labs and li == sr.LIB_HEART)) if heart is not None else 0.0)
    M = common_mask(masks)
    among = non_overlapping_groups(groups, {g: (cfg["external_labels"].get(g) or [None])[0] for g in groups})
    sv = split_vectors(dirs, M, weights, among)
    vecs = vector_set(dirs, sv)
    halves = {}
    for g in groups:  # per group, as share_read.direction_reliability calls it
        try:
            halves.update(feng_dirs.direction_halves(feng, genes, dict(cfg), sr.BRACKET, groups=[g]))
        except ValueError:
            pass
    hv = half_vectors(halves, groups, M, weights, among)
    rel = reliabilities(vecs, hv, M)
    energy_g = {"unweighted": {g: energy_parts(dirs[g], sv["c"], sv["s"][g], M) for g in groups},
                "_w": {g: energy_parts(dirs[g], sv["c_w"], sv["s_w"][g], M) for g in groups},
                "_n": {g: energy_parts(dirs[g], sv["c_n"], sv["s_n"][g], M) for g in groups}}
    setup = {"eps": cfg["metrics"]["pseudocount_cp10k"], "min_cells": cfg["metrics"]["min_cells_t1"],
             "groups": groups, "dirs": dirs, "masks": masks, "labels": labels, "weights": weights, "M": M,
             "vecs": vecs, "rel": rel, "support": list(t47["support"]), "energy_g": energy_g, "L": t47["L"],
             "L_corr": t47["L_corr"], "c_groups": list(groups), "c_n_groups": among, "halves": halves,
             "half_vecs": hv, "half_rel": half_reliabilities(hv, M)}
    extra = {"genes": genes, "tabs": tabs, "source": source,
             "halves_missing": [g for g in groups if g not in halves], "label_key": label_key}
    return setup, extra


def prepare(stages: Mapping[str, ad.AnnData], feng: dict, held: Sequence[str], cfg: Mapping,
            t47: Mapping) -> tuple[dict, dict, dict]:
    """build_setup, then the point's shares, the read set (finish_setup) and the point's values: (setup, extra,
    point)."""
    setup, extra = build_setup(stages, feng, held, cfg, t47)
    tabs = extra["tabs"]
    extra["point_sums"] = {s: sr.block_sums(t) for s, t in tabs.items()}
    shares = measure_shares(tabs, extra["point_sums"], setup)
    finish_setup(setup, shares)
    point = {**shares, **{f"reading/{k}": v for k, v in split_values(shares, setup).items()}}
    return setup, extra, point


def r_c_over_r_s(setup: Mapping) -> float:
    """The weighted r_c over the weighted r_s, both combined over the read set (NaN where either is undefined)."""
    read, w, rel = setup["read_set"], setup["weights"], setup["rel"]
    rc = strict_combine({g: rel["c"][g] for g in read}, read, w)
    rs = strict_combine({g: rel["s"][g] for g in read}, read, w)
    return _ratio(rc, rs)


def weight_shares(setup: Mapping) -> dict[str, float]:
    """Each read-set group's fraction of the read set's positive library-3 weight (a fraction, no cell count)."""
    w = {g: max(float(setup["weights"].get(g, 0.0)), 0.0) for g in setup["read_set"]}
    tot = sum(w.values())
    return {g: (x / tot if tot > 0 else float("nan")) for g, x in w.items()}


def split_reading(stages: Mapping[str, ad.AnnData], feng: dict, held: Sequence[str], cfg: Mapping, t47: Mapping,
                  n_boot: int = N_BOOT, seed=BOOT_SEED, calibrate: bool = False, self_check: bool = False) -> dict:
    """The whole read on loaded inputs: T1 stages, Feng's bracket statistics, M1's held list, T47's committed record
    (t47_record). calibrate / self_check: run T47's two arms read along the split (both need counts recoverable from
    .X; the self-check also runs the correction check); without the self-check the reading is void."""
    setup, extra, point = prepare(stages, feng, held, cfg, t47)
    tabs = extra["tabs"]
    rng = np.random.default_rng(list(seed))
    draws = {k: np.empty(n_boot) for k in point}
    t0 = time.time()
    for i in range(n_boot):
        sums = {s: sr.block_sums(t, sr.resample_weights(t, rng)) for s, t in tabs.items()}
        m = measure_split(tabs, sums, setup)
        for k in draws:
            draws[k][i] = m[k]
    seconds = time.time() - t0
    values = sr.summarise(point, draws, n_boot)
    r = {k.split("/", 1)[1]: values[k] for k in values if k.startswith("reading/")}
    sc = self_check_split(stages["E9.5"], extra["genes"], setup, held, extra["label_key"]) if self_check else None
    void = void_reasons(point, setup, t47["A"], sc)
    ratio = r_c_over_r_s(setup)
    st = statement(r, sc["bias"] if sc else None, n_boot, void, ratio)
    st["R_corr_without_each_group"] = {  # reported, read by nothing: the R_corr side with one group left out
        g: side_statement(r[f"R_corr_without/{g}"], (sc["bias"] if sc else {}).get("R_corr"), n_boot, "R_corr")
        for g in setup["read_set"] if f"R_corr_without/{g}" in r}
    res = {"groups": setup["groups"], "t47": {k: t47[k] for k in ("A", "L", "L_corr", "support")},
           "read_set": setup["read_set"], "p_set": setup["p_set"], "c_groups": setup["c_groups"],
           "c_n_groups": setup["c_n_groups"], "reliability": setup["rel"], "half_reliability": setup["half_rel"],
           "r_c_over_r_s": ratio, "weight_share": weight_shares(setup), "energy": setup["energy"],
           "energy_per_group": setup["energy_g"],
           "halves_missing": extra["halves_missing"], "values": values, "reading": r, "statement": st,
           "gauge": half_gauge(tabs, extra["point_sums"], setup),
           "n_boot": n_boot, "seed": list(seed), "percentiles": list(sr.PERCENTILES),
           "bootstrap_seconds": round(seconds, 1)}
    if calibrate:
        res["calibration"] = calibration_split(stages["E10.5"], extra["genes"], setup, held, extra["label_key"])
    if sc is not None:
        res["self_check"] = sc
        res["correction_check"] = correction_check(stages["E9.5"], extra["genes"], setup, held, extra["label_key"])
    return res


# ---------------------------------------------------------------------------------------------------------------------
# Text and the run


def _fmt(rec: Mapping, n_boot: int) -> str:
    return f"{rec['value']:+.3f} [{rec['p5']:+.3f}, {rec['p95']:+.3f}] ({rec['nonfinite_draws']} of {n_boot} undefined)"


def _num(x) -> str:
    return f"{x:+.4f}" if x is not None and np.isfinite(x) else "undefined"


def text(res: Mapping) -> str:
    r = res["read"]
    rd, st, nb = r["reading"], r["statement"], r["n_boot"]

    def f(k):
        return _fmt(rd[k], nb)

    rel, read, w = r["reliability"], r["read_set"], r["t47"]["support"]
    lines = [res["label"],
             f"Reading: {st['statement']}" + (f" [turned on {st['turned_on']}]" if st["turned_on"] else ""),
             f"    on R_corr (widened by the self-check's bias): {st['R_corr']}; on R (raw, widened): {st['R']}; "
             f"against {LESS_AT} and {AS_MUCH_AT}",
             f"    why not separated: {st['why_not_separated']}" if st.get("why_not_separated") else "",
             f"R {f('R')}; R_corr {f('R_corr')}; R / R_corr at the point "
             f"{_num(_ratio(rd['R']['value'], rd['R_corr']['value']))}; r_c / r_s (weighted) {_num(r['r_c_over_r_s'])}",
             f"A_M {f('A_M')} (T47's committed A {r['t47']['A']:+.4f}; A on T47's own masks {f('A_M_T47_mask')}); "
             f"A_M_corr {f('A_M_corr')}",
             f"A_c {f('A_c')}; A_s {f('A_s')}; A_c_corr {f('A_c_corr')}; A_s_corr {f('A_s_corr')}",
             f"    A, A_corr, K, R and R_corr all rest on the read set {read} (T47's A set {w}); (a) on {r['p_set']}",
             f"    c averages {r['c_groups']}; c_n averages {r['c_n_groups']}",
             f"Reliabilities (Feng's two halves; groups without halves: {r['halves_missing']}): r_c "
             f"{_num(next(iter(rel['c'].values()), float('nan')))}",
             "Per group of the read set (library-3 weight share; r_M, r_s; A's libraries 1 and 2 averaged, raw / "
             "corrected, along c and s_g):"]
    wt = r.get("weight_share", {})
    for g in read:
        lines.append(f"    {g}: weight share {_num(wt.get(g))}; r_M {_num(rel['M'][g])}, r_s {_num(rel['s'][g])}; "
                     f"along c {f(f'group/{g}/A_c')} / {f(f'group/{g}/A_c_corr')}; "
                     f"along s_g {f(f'group/{g}/A_s')} / {f(f'group/{g}/A_s_corr')}")
    lines += ["Energy parts (combined over the read set, library-3 weighted): "
              + "; ".join(f"{k} {v:+.4f}" for k, v in r["energy"]["unweighted"].items())
              + " | with c_w: " + "; ".join(f"{k} {v:+.4f}" for k, v in r["energy"]["_w"].items())
              + " | with c_n: " + "; ".join(f"{k} {v:+.4f}" for k, v in r["energy"]["_n"].items()),
              f"L {r['t47']['L']:+.4f}, L_corr {r['t47']['L_corr']:+.4f} (T47's committed); L_eff {f('L_eff')}; "
              f"L_eff_corr {f('L_eff_corr')} (conditional on the sample-difference reading; one-sided: it can only "
              f"lower L); uncapped, reported only: {f('L_eff_uncapped')}, {f('L_eff_uncapped_corr')}",
              f"    courses with A_M and L_eff: f_even {f('f_even')}; f_decay {f('f_decay')}; f_max {f('f_max')}",
              f"    courses with A_M_corr and L_eff_corr: f_even {f('f_even_corr')}; f_decay {f('f_decay_corr')}; "
              f"f_max {f('f_max_corr')}",
              "Reported and read by nothing:",
              "    R_corr with one group left out (its R_corr-side statement): "
              + "; ".join(f"without {g} {f(f'R_corr_without/{g}')} ({st.get('R_corr_without_each_group', {}).get(g)})"
                          for g in read if f"R_corr_without/{g}" in rd),
              f"    library 3: along c {f('lib3_c')}, along s_g {f('lib3_s')}, ratio {f('R_lib3')} (corrected "
              f"{f('R_lib3_corr')}); K_c {f('K_c')}, K_s {f('K_s')}, K_M {f('K_M')}",
              f"    (a) E8.5 -> E9.5: along c {f('P_c')}, along s_g {f('P_s')}, ratio {f('R_a')} (corrected "
              f"{f('R_a_corr')})",
              f"    with c_w: A_c_w {f('A_c_w')}; A_s_w {f('A_s_w')}; R_w {f('R_w')}; R_w_corr {f('R_w_corr')}; "
              f"L_eff_w {f('L_eff_w')}; f_even_w {f('f_even_w')}; f_decay_w {f('f_decay_w')}",
              f"    with c_n: A_c_n {f('A_c_n')}; A_s_n {f('A_s_n')}; R_n {f('R_n')}; R_n_corr {f('R_n_corr')}"]
    g_ = r["gauge"]
    lines.append("    Feng half gauge (vectors from one half alone, at the point): "
                 + "; ".join(f"{k} a {_num(g_['a'][k])}, b {_num(g_['b'][k])}, a - b {_num(g_['a_minus_b'][k])}"
                             for k in ("R", "R_corr", "A_c_corr", "A_s_corr")))
    lines.append(f"(value [p{sr.PERCENTILES[0]:g}, p{sr.PERCENTILES[1]:g}], {nb} draws, stream {r['seed']})")
    cal = r.get("calibration")
    if cal is not None:
        cw = cal["combined_weighted"]
        lines.append(f"Calibration: M1's step of {cal['fraction']} on E10.5 library 3 reads along d_g on M "
                     f"{_num(cw['M'])} (T47's mask {_num(cw['d_T47'])}), along c {_num(cw['c'])}, along s_g "
                     f"{_num(cw['s'])}, ratio {_num(cal['R'])} (with c_w {_num(cal['R_w'])}, with c_n "
                     f"{_num(cal['R_n'])}; over {cal['combined_over']})")
    sc = r.get("self_check")
    if sc is not None:
        lines.append(f"Self-check (T47's world: H2 against H1 moved along d by {sc['fraction']} on stream "
                     f"{sc['stream']}, split on stream {sc['split_stream']}; every r_v = 1, so it does not test the "
                     f"noise correction; at the point; over {sc['support']}): "
                     + "; ".join(f"{k} {_num(sc['measured'][k])} (known {sc['known'][k]:+.3f}, bias "
                                 f"{_num(sc['bias'][k])}; by definition {_num(sc['known_by_definition'][k])})"
                                 for k in SELF_KEYS))
    cc = r.get("correction_check")
    if cc is not None:
        lines.append(f"Correction check, read by nothing (H1 moved along half A's direction, read along half B's "
                     f"vectors with single-half reliabilities; over {cc['support']}): "
                     + "; ".join(f"{k} {_num(cc['measured'][k])} (bias {_num(cc['bias'][k])}, from the definition "
                                 f"{_num(cc['bias_from_definition'][k])})" for k in ("A_c", "A_s", "R_corr")))
    v = r["values"]
    for name in sr.COMPARISONS:
        lines.append(f"{name}:")
        for g in r["groups"]:
            lines.append(f"    {g}: " + "; ".join(f"{vv} {_fmt(v[f'{name}/{g}/{vv}/share'], nb)}"
                                                 for vv in ("d_T47", *VECS)))
    lines.extend(f"Limit: {s}" for s in LEDGER_NOTES)
    lines.append(f"Commit {res['stamp']['commit']}, tree clean {res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}")
    return "\n".join(x for x in lines if x) + "\n"


def run(out_name: str, files: Mapping[str, Path] = sr.T1_FILES) -> dict:
    if ROW not in n12.tries_rows():
        sys.exit(f"{ROW} is refused: it has no row in experiments/harness/TRIES.md (a row comes before its run)")
    out = RESULTS / out_name
    gv.refuse_existing(out, (f"{STEM}.json", f"{STEM}.txt"))
    started = gv.tree_state()
    if started["clean"] is not True:
        sys.exit(f"{ROW} is refused: the tree is not clean (a counted read runs on a committed tree)")
    cfg = feng_dirs.load_config()
    union = hygiene.load_union()
    held = gc.held_feng(union, t1_scope.held_flat_genes())
    feng, feng_rec = sr.feng_inputs()
    t47 = t47_record()
    stages = {s: ad.read_h5ad(files[s]) for s in sr.T1_STAGES}
    read = split_reading(stages, feng, held, cfg, t47, calibrate=True, self_check=True)
    res = {"label": LABEL_T51, "read": read, "feng": feng_rec}
    res["stamp"] = gv.stamp(started, [*files.values(), feng_dirs.clean_stats_path(), t1_scope.HELD_FLAT_GENES,
                                      sr.REPO / union["file"], T47_JSON])
    gv.write_outputs(out, STEM, res, text(res))
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="try T51's read")
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    run(a.out)


if __name__ == "__main__":
    main()
