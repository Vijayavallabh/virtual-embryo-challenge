"""Per-gene timing of T1's recomputed Feng direction from Feng's post-window slope (ledger row T67; off everywhere).

The recomputed-groups file moves each recomputed Feng group by a fraction of the E9.5 -> E14.5 log chord D (natural
log of mean CP10k plus the readout's pseudocount, per gene), as if each gene changed at a steady rate. This module holds
the arithmetic of the moonshot that replaces that steady rate with a cubic Hermite through the bracket:

  - m0, the early slope, is fixed at the linear value D / 5 per day (the bracket spans five days);
  - m1, the late slope, is the per-day ordinary least-squares slope of log(mean CP10k + eps) over Feng's E14.5, E15.5,
    E16.5 and E17.5 cells of the same recomputed group (`ols_slope`);
  - the Hermite's increment over the test span tau 0.2 -> 0.6 is w = 0.512 D - 0.56 m1 (`hermite_increment`); a
    steady-rate gene (m1 = D / 5) gets exactly 0.4 D;
  - the monotone limiter (`limit`): m1 whose sign is not D's is set to 0; |m1| is capped at 3 |D| / 5;
  - the reliability fallback (`timing_m1`): m1 := D / 5 for a gene whose two per-cell split halves give slopes of
    different sign or more than |m1| apart, and for every gene of a group with fewer than MIN_LATE_CELLS cells at any
    late stage;
  - per group, w is rescaled to the entry's own L2 size (`norm_match`), so only the genes' order and relative sizes
    change, never the step's size. As a direction (the object the builder multiplies by the bracket fraction) that is
    `timing_direction` = shape x ||D|| / ||shape||, shape = w / 0.4 = D - 1.4 (m1 - D / 5). The two forms are equal in
    exact arithmetic; the second makes m1 = D / 5 give D back bit for bit (the row's void condition).

The control (`permuted_m1`): the deviations m1 - D / 5 permuted among genes within each group, inside PERM_D_BINS
equal-count bins of |D| crossed with PERM_LEVEL_BINS equal-count bins of the group's Feng E9.5 level, seed
[PERM_SEED, the group's index in GROUP_ORDER]; then the same limiter and norm match.

The idle pre-gate (`idle_ratio`): ||direction - D|| / ||D|| per group; below IDLE_BAR in every group the arm is idle.

Arrays only: nothing here reads a file, a cell or a target. NaN (an unmapped gene) stays NaN throughout.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

BRACKET_DAYS = 5.0  # E9.5 -> E14.5
SLOPE_STAGES = ("E14.5", "E15.5", "E16.5", "E17.5")
SLOPE_DAYS = (14.5, 15.5, 16.5, 17.5)
LATE_STAGES = SLOPE_STAGES[1:]  # the stages the late labelling covers; E14.5 is the bracket's own end
HERMITE_D = 0.512  # H01(0.6) - H01(0.2) + (H10(0.6) - H10(0.2)) x 5 x (1/5)
HERMITE_M1 = 0.56  # -(H11(0.6) - H11(0.2)) x 5
TEST_FRACTION = 0.4  # the test span's share of the bracket (tau 0.2 -> 0.6)
SHAPE_M1 = HERMITE_M1 / TEST_FRACTION  # 1.4: w / 0.4 = D - 1.4 (m1 - D / 5)
CAP = 3.0 / 5.0  # |m1| <= 3 |D| / 5
MIN_LATE_CELLS = 50
PERM_SEED = 20261122
PERM_D_BINS = 20
PERM_LEVEL_BINS = 5
IDLE_BAR = 0.05
GROUP_ORDER = ("CM", "Endocardium", "Endothelium", "Epicardium")


def steady_m1(d) -> np.ndarray:
    """The linear (steady-rate) slope per day, D / 5."""
    return np.asarray(d, dtype=np.float64) / BRACKET_DAYS


def ols_slope(levels, days: Sequence[float] = SLOPE_DAYS) -> np.ndarray:
    """Per gene (columns), the ordinary least-squares slope of `levels` (one row per day of `days`) on the day; NaN for
    a gene with a non-finite level on any day."""
    y = np.asarray(levels, dtype=np.float64)
    x = np.asarray(days, dtype=np.float64)
    if y.ndim != 2 or y.shape[0] != x.size or x.size < 2:
        raise ValueError(f"levels must hold one row per day ({x.size}); got shape {y.shape}")
    xc = x - x.mean()
    out = (xc[:, None] * (y - y.mean(axis=0))).sum(axis=0) / float((xc * xc).sum())
    out[~np.isfinite(y).all(axis=0)] = np.nan
    return out


def log_level(cp, eps: float) -> np.ndarray:
    """log(mean CP10k + eps), the unit of D (dirlib.lfc_ratio)."""
    return np.log(np.asarray(cp, dtype=np.float64) + float(eps))


def limit(m1, d) -> np.ndarray:
    """The monotone limiter: m1 whose sign is not D's becomes 0; |m1| is capped at 3 |D| / 5. NaN stays NaN."""
    m1, d = np.asarray(m1, dtype=np.float64), np.asarray(d, dtype=np.float64)
    out = np.where(np.sign(m1) != np.sign(d), 0.0, m1)
    cap = CAP * np.abs(d)
    out = np.clip(out, -cap, cap)
    out[~(np.isfinite(m1) & np.isfinite(d))] = np.nan
    return out


def unreliable(m1, m1_a, m1_b) -> np.ndarray:
    """The split-half check: the halves' slopes differ in sign or by more than |m1|; a non-finite value fails."""
    m1, a, b = (np.asarray(x, dtype=np.float64) for x in (m1, m1_a, m1_b))
    finite = np.isfinite(m1) & np.isfinite(a) & np.isfinite(b)
    with np.errstate(invalid="ignore"):
        bad = (np.sign(a) != np.sign(b)) | (np.abs(a - b) > np.abs(m1))
    return ~finite | bad


def timing_m1(d, m1_raw, m1_a, m1_b, group_ok: bool) -> tuple[np.ndarray, np.ndarray]:
    """The arm's m1 per gene and the genes that fell back to D / 5: every gene when the group fails the late-cell
    minimum (`group_ok` False); otherwise the genes the split-half check marks (`unreliable`, on the measured m1), the
    rest limited (`limit`). NaN where D is NaN."""
    d = np.asarray(d, dtype=np.float64)
    steady = steady_m1(d)
    finite = np.isfinite(d)
    fell = finite.copy() if not group_ok else finite & unreliable(m1_raw, m1_a, m1_b)
    out = np.where(fell, steady, limit(m1_raw, d))
    out[~finite] = np.nan
    return out, fell


def hermite_increment(d, m1) -> np.ndarray:
    """The cubic Hermite's increment over tau 0.2 -> 0.6 (m0 = D / 5): w = 0.512 D - 0.56 m1 (reported)."""
    return HERMITE_D * np.asarray(d, dtype=np.float64) - HERMITE_M1 * np.asarray(m1, dtype=np.float64)


def timing_shape(d, m1) -> np.ndarray:
    """w / 0.4 written as D - 1.4 (m1 - D / 5), so that m1 = D / 5 gives D itself."""
    d = np.asarray(d, dtype=np.float64)
    return d - SHAPE_M1 * (np.asarray(m1, dtype=np.float64) - steady_m1(d))


def l2(v) -> float:
    v = np.asarray(v, dtype=np.float64)
    ok = np.isfinite(v)
    return float(np.sqrt(np.sum(v[ok] * v[ok])))


def norm_match(shape, d) -> np.ndarray:
    """`shape` rescaled to D's L2 size over the finite genes (NaN stays NaN). A shape of size 0 is refused."""
    shape = np.asarray(shape, dtype=np.float64)
    size = l2(shape)
    if not size > 0:
        raise ValueError("the timing shape has no size: it cannot be matched to the entry's")
    return shape * (l2(d) / size)


def timing_direction(d, m1) -> np.ndarray:
    """The direction the builder applies in place of D (at its own bracket fraction): the Hermite shape, norm-matched
    to D. Its test step at 2/5 is w rescaled to ||0.4 D||; at the practice 1/5 it is the same shape at ||0.2 D||."""
    d = np.asarray(d, dtype=np.float64)
    out = norm_match(timing_shape(d, m1), d)
    gone = ~np.isfinite(d)
    out[gone] = d[gone]  # an unmapped gene keeps D's own NaN, bit for bit
    return out


def idle_ratio(direction, d) -> float:
    """||direction - D|| / ||D|| over D's finite genes (= ||w_scaled - 0.4 D|| / ||0.4 D||)."""
    d = np.asarray(d, dtype=np.float64)
    ok = np.isfinite(d)
    size = l2(d[ok])
    return l2(np.asarray(direction, dtype=np.float64)[ok] - d[ok]) / size if size > 0 else float("nan")


def equal_count_bins(values, bins: int) -> np.ndarray:
    """Each entry's bin among `bins` equal-count runs of the stable ranking of `values` (all finite)."""
    v = np.asarray(values, dtype=np.float64)
    out = np.empty(v.size, dtype=np.int64)
    for b, run in enumerate(np.array_split(np.argsort(v, kind="stable"), bins)):
        out[run] = b
    return out


def permuted_m1(d, m1, level, group: str) -> np.ndarray:
    """The control's m1: the deviations m1 - D / 5 of the genes with finite D, m1 and level, permuted within the cells
    PERM_D_BINS bins of |D| crossed with PERM_LEVEL_BINS bins of `level` (the group's Feng E9.5 mean CP10k), seed
    [PERM_SEED, GROUP_ORDER.index(group)], added back to each gene's own D / 5 and limited (`limit`)."""
    if group not in GROUP_ORDER:
        raise ValueError(f"{group}: not one of {GROUP_ORDER}: it has no permutation seed")
    d, m1, level = (np.asarray(x, dtype=np.float64) for x in (d, m1, level))
    delta = m1 - steady_m1(d)
    ok = np.flatnonzero(np.isfinite(d) & np.isfinite(m1) & np.isfinite(level))
    cell = (equal_count_bins(np.abs(d[ok]), PERM_D_BINS) * PERM_LEVEL_BINS
            + equal_count_bins(level[ok], PERM_LEVEL_BINS))
    rng = np.random.default_rng([PERM_SEED, GROUP_ORDER.index(group)])
    moved = delta.copy()
    for c in np.unique(cell):
        members = ok[cell == c]
        moved[members] = delta[members[rng.permutation(members.size)]]
    return limit(steady_m1(d) + moved, d)


def group_timing(group: str, d, m1_raw, m1_a, m1_b, level, late_cells: Mapping[str, int]) -> dict:
    """One group's arm and control: the late-cell check over LATE_STAGES, the arm's m1 and its fallback genes, the
    arm's and the control's directions, both idle ratios and the reported Hermite increment."""
    d = np.asarray(d, dtype=np.float64)
    missing = [s for s in LATE_STAGES if s not in late_cells]
    if missing:
        raise ValueError(f"{group}: no late cell count for {missing}")
    group_ok = all(int(late_cells[s]) >= MIN_LATE_CELLS for s in LATE_STAGES)
    m1, fell = timing_m1(d, m1_raw, m1_a, m1_b, group_ok)
    m1_perm = permuted_m1(d, m1, level, group)
    direction, perm_direction = timing_direction(d, m1), timing_direction(d, m1_perm)
    finite = np.isfinite(d)
    return {"m1": m1, "fell_back": fell, "m1_perm": m1_perm, "direction": direction, "perm_direction": perm_direction,
            "w": hermite_increment(d, m1), "group_ok": bool(group_ok),
            "fallback_share": float(fell[finite].mean()) if finite.any() else float("nan"),
            "idle_ratio": idle_ratio(direction, d), "perm_idle_ratio": idle_ratio(perm_direction, d),
            "genes": int(finite.sum())}


def is_idle(ratios: Mapping[str, float], bar: float = IDLE_BAR) -> bool:
    """The pre-gate: idle when every group's ratio is below the bar (a NaN ratio is not below it)."""
    if not ratios:
        raise ValueError("no group to read")
    return all(r < bar for r in ratios.values())


def steady_identity(d) -> bool:
    """m1 = D / 5 for every gene gives D back bit for bit (NaN positions included)."""
    d = np.asarray(d, dtype=np.float64)
    return np.ascontiguousarray(timing_direction(d, steady_m1(d))).tobytes() == np.ascontiguousarray(d).tobytes()
