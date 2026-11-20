"""The fixed regularised forms of the barred factor (the Mab21l2 knockout against WT E8.75), one implementation shared
by the practice read (k21_factor_forms.py) and the thin builder (build_factor_form.py), so that a file built in a form
applies exactly the factor its practice arm read.

Pure functions of a factor vector and module constants; nothing here loads a file, draws a cell or reads a score.

Notation: phi_M is the table `transfer_mab21l2_all` as the frozen loader returns it (t3_decompose.gene_factor already
rescaled it and clipped it to [0.25, 4]); l = ln phi_M; SD is the population SD (ddof 0) over the genes. No form is
rescaled or clipped after its transform, as the exponent builder and the two-knockout builder treat every factor.

  identity   phi_M
  top100     phi_M on the TOP_K genes of largest |l| (stable sort, ties in the vector's order), 1 elsewhere. Refused
             when TOP_K or more genes sit at the clip, or when the |l| at rank TOP_K equals the one at rank TOP_K + 1
             (the kept set would then follow the gene order, not the factor)
  uni_top    phi_M ** s_t, s_t = SD(ln top100) / SD(l): the control of top100 at its SD (an exponent below 1, a control
             only, never built or named)
  cap90      t = the CAP_Q quantile (numpy's linear) of |l| over the genes off the clip (|l| < ln 4 - CLIP_TOL);
             phi = clip(phi_M, e^-t, e^t). Refused when no gene is off the clip
  uni_cap    phi_M ** s_c, s_c = SD(ln cap90) / SD(l): the control of cap90 at its SD (a control only)
  wtprod     not a factor vector: the cells are drawn with no factor, then each cell is moved per lineage group at the
             wild-type-only factor's exponent WT_EXPONENT (0.5) with per'_g = phi_M**2 * per_g and the fallback
             phi_M**2 * phi_pooled, which applies phi_M * per_g ** 0.5 in one step (wtprod_factors)
At a base of exponent x the builder applies T(phi_M) ** x, the transform's constants (t, the kept genes, s) fixed from
phi_M itself; wtprod at x uses phi_M ** x in place of phi_M.
"""

from __future__ import annotations

import math

import numpy as np

CLIP = (0.25, 4.0)            # t3_decompose.gene_factor's clip, as the table was written
LN_CLIP = math.log(CLIP[1])   # = -ln CLIP[0]
CLIP_TOL = 1e-9               # |l| >= ln 4 - CLIP_TOL counts as at the clip
TOP_K = 100                   # of the panel's 500 genes
CAP_Q = 0.9                   # the quantile of |l| off the clip
WT_EXPONENT = 0.5             # k7_wt_revert.EXPONENT, the wild-type-only factor's exponent (a test pins the equality)
MIN_CHANGED = 10              # pre-gate: an arm that changes fewer genes than this against phi_M is not run

TRIES = ("top100", "cap90", "wtprod")                       # the forms a reading may name
CONTROL_OF = {"top100": "uni_top", "cap90": "uni_cap", "wtprod": "wtprod_perm"}
VECTOR_FORMS = ("identity", "top100", "uni_top", "cap90", "uni_cap")


class FormRefused(ValueError):
    """A form that cannot be made from this factor (a pre-gate drop, or a builder refusal)."""


def _ln(phi) -> np.ndarray:
    phi = np.asarray(phi, dtype=np.float64)
    if phi.ndim != 1 or phi.size == 0 or not np.all(np.isfinite(phi) & (phi > 0)):
        raise FormRefused("the factor must be a non-empty vector of finite positive values")
    return np.log(phi)


def sd(v) -> float:
    return float(np.std(np.asarray(v, dtype=np.float64), ddof=0))


def at_clip(phi) -> np.ndarray:
    """The genes whose factor sits at either clip end."""
    return np.abs(_ln(phi)) >= LN_CLIP - CLIP_TOL


def top_k(phi, k: int | None = None) -> tuple[np.ndarray, dict]:
    """phi on the k genes of largest |ln phi| (np.argsort(-|l|, kind='stable')), 1.0 elsewhere."""
    k = TOP_K if k is None else int(k)
    phi = np.asarray(phi, dtype=np.float64)
    a = np.abs(_ln(phi))
    n = a.size
    n_clip = int(at_clip(phi).sum())
    if k >= n:
        return phi.copy(), {"k": k, "n_genes": n, "genes_at_clip": n_clip, "cut": None}
    if n_clip >= k:
        raise FormRefused(f"{n_clip} genes sit at the clip, at least the {k} kept: the kept set would follow the "
                          "gene order, not the factor")
    order = np.argsort(-a, kind="stable")
    if a[order[k - 1]] == a[order[k]]:
        raise FormRefused(f"|ln phi| ties across the cut at rank {k}: the kept set would follow the gene order")
    out = np.ones(n)
    keep = order[:k]
    out[keep] = phi[keep]
    return out, {"k": k, "n_genes": n, "genes_at_clip": n_clip, "cut": float(a[order[k - 1]]),
                 "next": float(a[order[k]])}


def cap_at(phi, t: float) -> np.ndarray:
    """phi clipped to [e^-t, e^t]."""
    return np.clip(np.asarray(phi, dtype=np.float64), math.exp(-t), math.exp(t))


def cap(phi, q: float | None = None) -> tuple[np.ndarray, dict]:
    """t = the q quantile (numpy's linear) of |ln phi| over the genes off the clip; phi clipped to [e^-t, e^t]."""
    q = CAP_Q if q is None else float(q)
    a = np.abs(_ln(phi))
    off = ~at_clip(phi)
    if not off.any():
        raise FormRefused("no gene lies off the clip: the cap's quantile is undefined")
    t = float(np.quantile(a[off], q))
    return cap_at(phi, t), {"q": q, "t": t, "genes_off_clip": int(off.sum()), "genes_at_clip": int((~off).sum())}


def uniform(phi, target) -> tuple[np.ndarray, dict]:
    """phi ** s with s = SD(ln target) / SD(ln phi): the base's shape at the target's SD."""
    s_phi = sd(_ln(phi))
    if not s_phi > 0:
        raise FormRefused("SD(ln phi) is 0: the size-matched exponent is undefined")
    s = sd(_ln(target)) / s_phi
    return np.asarray(phi, dtype=np.float64) ** s, {"s": float(s)}


def vector_form(name: str, phi) -> tuple[np.ndarray, dict]:
    """The factor of a vector form, from phi (= phi_M) as loaded; returns (factor, record)."""
    phi = np.asarray(phi, dtype=np.float64)
    _ln(phi)
    if name == "identity":
        return phi, {}
    if name == "top100":
        return top_k(phi)
    if name == "uni_top":
        t, rt = top_k(phi)
        out, ru = uniform(phi, t)
        return out, {**rt, **ru}
    if name == "cap90":
        return cap(phi)
    if name == "uni_cap":
        c, rc = cap(phi)
        out, ru = uniform(phi, c)
        return out, {**rc, **ru}
    raise FormRefused(f"{name!r} is not a vector form ({list(VECTOR_FORMS)})")


def at_exponent(phi_t: np.ndarray, x: float) -> np.ndarray:
    """T(phi_M) ** x, the frozen build's own operation for w = x; at x = 1 the very array (no pow)."""
    return phi_t if x == 1.0 else np.asarray(phi_t, dtype=np.float64) ** x


def wtprod_factors(phi_mx, per: dict, pooled) -> tuple[dict, np.ndarray]:
    """(per'_g, fallback') = (phi_mx**2 * per_g, phi_mx**2 * pooled): raised to WT_EXPONENT (0.5) by the per-group
    step, each group's applied factor is phi_mx * per_g ** 0.5 (up to rounding). phi_mx is phi_M ** x."""
    sq = np.asarray(phi_mx, dtype=np.float64) ** 2
    _ln(sq)
    return {g: sq * np.asarray(v, dtype=np.float64) for g, v in per.items()}, sq * np.asarray(pooled, dtype=np.float64)


def describe(phi_m, phi) -> dict:
    """Target-free figures of a factor against phi_M: genes changed, genes at the clip, SD and gene mean of ln phi
    (the common shift), the min and max."""
    phi_m, phi = np.asarray(phi_m, dtype=np.float64), np.asarray(phi, dtype=np.float64)
    lp = _ln(phi)
    return {"genes_changed": int((phi != phi_m).sum()), "genes_set_to_1": int(((phi == 1.0) & (phi_m != 1.0)).sum()),
            "genes_at_clip": int(at_clip(phi).sum()), "sd_ln": sd(lp), "gene_mean_ln": float(lp.mean()),
            "min": float(phi.min()), "max": float(phi.max())}


def pregate(phi_m, per: dict | None = None, pooled=None) -> dict:
    """The pre-gate, from phi_M alone (and, for wtprod's printed range only, the wild-type-only factors): each arm's
    figures and the arms dropped. A form that cannot be made, or a try or its control changing fewer than MIN_CHANGED
    genes against phi_M, drops the try with its control; a dropped arm is not run and not counted. wtprod is never
    dropped here: its applied factor differs from phi_M wherever a group's own factor is not 1."""
    phi_m = np.asarray(phi_m, dtype=np.float64)
    out: dict = {"n_genes": int(phi_m.size), "genes_at_clip_phi_m": int(at_clip(phi_m).sum()),
                 "phi_m": describe(phi_m, phi_m), "arms": {}, "dropped": {}}
    for name in VECTOR_FORMS[1:]:
        try:
            phi, rec = vector_form(name, phi_m)
            out["arms"][name] = {**rec, **describe(phi_m, phi)}
        except FormRefused as e:
            out["arms"][name] = {"refused": str(e)}
    for t in ("top100", "cap90"):
        c = CONTROL_OF[t]
        why = [f"{a}: {out['arms'][a]['refused']}" for a in (t, c) if "refused" in out["arms"][a]]
        why += [f"{a} changes {out['arms'][a]['genes_changed']} genes (< {MIN_CHANGED})" for a in (t, c)
                if "refused" not in out["arms"][a] and out["arms"][a]["genes_changed"] < MIN_CHANGED]
        if why:
            out["dropped"][t] = out["dropped"][c] = "; ".join(why)
    w: dict = {"applied": "phi_M * per_g ** 0.5 per lineage group (fallback phi_M * phi_pooled ** 0.5)"}
    if per is not None and pooled is not None:
        p2, f2 = wtprod_factors(phi_m, per, pooled)
        applied = {g: v ** WT_EXPONENT for g, v in {**p2, "pooled_fallback": f2}.items()}
        w |= {"min": float(min(v.min() for v in applied.values())),
              "max": float(max(v.max() for v in applied.values())), "groups": sorted(per),
              "sd_ln_by_group": {g: sd(np.log(v)) for g, v in applied.items()}}
    out["arms"]["wtprod"] = w
    return out


def active_arms(pre: dict, all_arms: tuple[str, ...]) -> tuple[str, ...]:
    """The arms that run: every arm but the pre-gate's drops, in the fixed order."""
    return tuple(a for a in all_arms if a not in pre["dropped"])
