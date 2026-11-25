"""Row E57's step, one run per side (experiments/harness/TRIES.md, row E57; the score-blind ERROR designer's proposal P1
of the embryo design round of 2026-11-24/25, as its cold challenger fixed it): on O's final cloud (fgw_shape0x_ov, the
overlay of the two brackets' real rows), each bracket's expression rows are moved to its own side of the earlier
bracket's midline plane, within fixed coarse-lineage pools. Positions never move: only which of O's expression rows
sits at which of O's positions changes. The coordinates (their bytes and their order), the count, the expression rows
as a multiset and each position's pool stay O's, and every row keeps its own expression. Nothing here is random: the
normals are a fixed Fibonacci hemisphere and the pairings are scipy's linear_sum_assignment.

The step (sort_rows(..., mode="run_sorted")), on the final coordinates (after the shape form, the size rule, the
duplicate jitter and the cloud guard: the place --row-order takes), in the row's words:
  (1) u = the normal, among N_NORMALS fixed Fibonacci hemisphere normals, of the plane through the centroid that best
      mirrors the earlier bracket's drawn rows onto rows of the same released label of that file (score: the median
      over labels with at least MIN_LABEL_ROWS drawn rows of the median distance from a mirrored row to its label's
      nearest row);
  (2) earlier-bracket rows to side (x - c).u < 0, later-bracket rows to the other;
  (3) within each pool of POOLS, wrong-side rows of the two brackets are paired by scipy linear_sum_assignment on
      |x_j - mirror(x_i)|^2 and each pair swaps expression rows; leftovers stay.
The placebo (mode="run_sorted_same") makes the same number of swaps per pool and side, each between two rows of one
bracket across the plane, by the same cost. mode="as_built" changes nothing (O itself).

Readings where the row left a choice, fixed with this module before any data were read:
  * c is the centroid of the final cloud (every row, float64). A row exactly on the plane, (x - c).u = 0, is on the
    later bracket's side.
  * The normals: the Fibonacci lattice on the upper unit hemisphere, z_i = 1 - (i + 1/2)/n and azimuth i times the
    golden angle (a plane's two normals are one plane, so a hemisphere covers every plane once). The smallest score
    wins, the first normal on a tie (numpy's argmin).
  * A label's nearest row may be the mirrored row's own position (a row on the plane is its own mirror image).
  * The contrast is the best score over the median score of all normals. The later bracket's contrast is its score at
    the earlier bracket's u, on its own file's labels, over the median of its own scores (the cold challenge's fix 3).
  * A row's pool is the pool of POOLS that holds its lineage-map group (vec.prep's obs["group"]); a group in no pool
    refuses the step (RunSortRefused).
  * The placebo, per pool with RS's k swaps: ceil(k/2) swaps between rows of the earlier bracket and floor(k/2) between
    rows of the later bracket, so each bracket has about as many rows displaced as under RS; where one bracket lacks
    pairs across the plane, the other takes the rest, and any swap still missing is recorded (placebo_shortfall).
    A bracket's swaps are the cheapest pairs of its own linear_sum_assignment across the plane (cost
    |x_j - mirror(x_i)|^2, i on the side (x - c).u >= 0, j on the other), ties in the order the assignment returns them.
    Every swap crosses the plane, so each side gives up as many rows as under RS.
  * A swap's homology error is |x_j - mirror(x_i)|, the same for both of its rows.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

ROW = "E57"
MODES = ("as_built", "run_sorted", "run_sorted_same")  # O; the arm E57b (RS); its placebo E57c (rs_same)
BUILDABLE = MODES[:2]  # the placebo is harness-only and never built
N_NORMALS = 2000  # fixed Fibonacci hemisphere normals (the row's "2,000 fixed Fibonacci hemisphere normals")
MIN_LABEL_ROWS = 20  # a label scores the plane only with at least this many drawn rows (the row's "at least 20")
K = 15  # the scorer's neighbourhood: a cell and its 14 nearest others (N15's f, the 15-nearest-cell reach)
CHUNK = 64  # normals mirrored at once (memory only; the result does not depend on it)
# The pools the cold challenger fixed (row E57, cell 3): every group of configs/lineage_map.yaml in exactly one pool
POOLS = (("EpiblastStreak", "Neural", "Ectoderm", "NCC"),
         ("CM", "CardiacProgenitor", "Endocardium", "Endothelium", "EpicardiumST", "Pericardium", "PharyngealMesoderm",
          "ParaxialMesoderm", "LPM"),
         ("ExtraEmbryonicMeso",),
         ("Endoderm", "Hepatocyte"),
         ("Blood",),
         ("Unknown",))
RULE = ("vec.run_sort (row E57): on the final coordinates, u = the Fibonacci-hemisphere normal (2000) of the plane "
        "through the cloud's centroid that best mirrors the earlier bracket's rows onto rows of their own released "
        "label (median over labels with >= 20 rows of the median distance to the label's nearest row); earlier rows "
        "to (x - c).u < 0, later rows to the other side; within each fixed pool the wrong-side rows of the two "
        "brackets paired by linear_sum_assignment on |x_j - mirror(x_i)|^2 swap expression rows; leftovers stay; "
        "positions, count and the expression-row multiset untouched")


class RunSortRefused(ValueError):
    """The step cannot be taken as the row defines it (a group in no pool, no label to score the plane, bad input)."""


@dataclass
class RunSort:
    """The step's result: `src` such that X[src] puts at position i the expression row the step assigns to it (O's own
    order for mode="as_built"); `pairs` (k, 2) the swapped positions (i on the side (x - c).u >= 0, j on the other);
    `errors` each swap's homology error; `stage_after` the bracket of the row at each position afterwards; and a
    record of flat values (no cell number but the per-pool tallies in `pools`)."""

    src: np.ndarray
    pairs: np.ndarray
    errors: np.ndarray
    stage_after: np.ndarray
    normal: np.ndarray
    centre: np.ndarray
    record: dict
    pools: dict = field(default_factory=dict)


def fibonacci_hemisphere(n: int = N_NORMALS) -> np.ndarray:
    """n unit normals on the upper hemisphere: z_i = 1 - (i + 1/2)/n, azimuth i times the golden angle."""
    i = np.arange(int(n), dtype=np.float64)
    z = 1.0 - (i + 0.5) / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    phi = i * np.pi * (3.0 - np.sqrt(5.0))
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def mirror(P, c, u) -> np.ndarray:
    """P reflected through the plane through c with unit normal u."""
    P = np.asarray(P, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    return P - 2.0 * ((P - c) @ u)[:, None] * u[None, :]


def pool_index(groups) -> np.ndarray:
    """Each row's pool (an index into POOLS) from its lineage-map group; a group in no pool refuses."""
    groups = np.asarray(groups).astype(str)
    of = {g: p for p, members in enumerate(POOLS) for g in members}
    unknown = sorted(set(groups) - set(of))
    if unknown:
        raise RunSortRefused(f"groups in no pool of vec.run_sort.POOLS: {unknown}")
    return np.array([of[g] for g in groups], dtype=np.int64)


def plane_scores(P, labels, c, normals, min_rows: int = MIN_LABEL_ROWS) -> np.ndarray:
    """Per normal, the row's score of the plane through c: the median over labels with at least `min_rows` rows of the
    median distance from a mirrored row to the nearest row of its own label (unmirrored, the row itself included)."""
    P = np.asarray(P, dtype=np.float64)
    labels = np.asarray(labels).astype(str)
    normals = np.asarray(normals, dtype=np.float64)
    keep = [lab for lab in sorted(set(labels)) if int((labels == lab).sum()) >= min_rows]
    if not keep:
        raise RunSortRefused(f"no label holds {min_rows} rows: the plane cannot be scored")
    per = np.empty((len(keep), len(normals)))
    for li, lab in enumerate(keep):
        Q = P[labels == lab]
        tree = cKDTree(Q)
        proj = (Q - c) @ normals.T  # (rows, normals)
        for a in range(0, len(normals), CHUNK):
            U = normals[a:a + CHUNK]
            M = Q[None, :, :] - 2.0 * proj[:, a:a + CHUNK].T[:, :, None] * U[:, None, :]
            d = tree.query(M.reshape(-1, 3), k=1)[0].reshape(len(U), len(Q))
            per[li, a:a + CHUNK] = np.median(d, axis=1)
    return np.median(per, axis=0)


def midline(P, labels, c, normals=None) -> tuple[np.ndarray, float, np.ndarray]:
    """(u, contrast, scores): the best plane's normal, the best score over the median score, every normal's score."""
    normals = fibonacci_hemisphere() if normals is None else np.asarray(normals, dtype=np.float64)
    scores = plane_scores(P, labels, c, normals)
    best = int(np.argmin(scores))
    med = float(np.median(scores))
    return normals[best].copy(), (float(scores[best]) / med if med > 0 else float("nan")), scores


def contrast_at(P, labels, c, u, normals=None) -> float:
    """Another bracket's contrast at a given plane: its score at u over the median of its own scores."""
    normals = fibonacci_hemisphere() if normals is None else np.asarray(normals, dtype=np.float64)
    scores = plane_scores(P, labels, c, normals)
    at = float(plane_scores(P, labels, c, np.asarray(u, dtype=np.float64)[None, :])[0])
    med = float(np.median(scores))
    return at / med if med > 0 else float("nan")


def _assign(C, plus, minus, c, u) -> tuple[np.ndarray, np.ndarray]:
    """linear_sum_assignment of the rows `plus` against the rows `minus` on |x_j - mirror(x_i)|^2: the pairs (as
    positions) and each pair's cost, in the order the assignment returns them."""
    if len(plus) == 0 or len(minus) == 0:
        return np.empty((0, 2), dtype=np.int64), np.empty(0)
    cost = cdist(mirror(C[plus], c, u), C[minus], "sqeuclidean")
    r, k = linear_sum_assignment(cost)
    return np.column_stack([plus[r], minus[k]]).astype(np.int64), cost[r, k]


def _checked(coords, stage, labels, groups) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    C = np.asarray(coords, dtype=np.float64)
    if C.ndim != 2 or C.shape[1] != 3 or len(C) == 0 or not np.isfinite(C).all():
        raise RunSortRefused("coordinates must be finite and (n, 3)")
    stage = np.asarray(stage)
    if stage.shape != (len(C),) or not np.isin(stage, (0, 1)).all():
        raise RunSortRefused("stage must hold 0 (earlier bracket) or 1 (later bracket) for every row")
    labels = np.asarray(labels).astype(str)
    if labels.shape != (len(C),):
        raise RunSortRefused("one released label per row")
    pools = pool_index(groups)
    if pools.shape != (len(C),):
        raise RunSortRefused("one lineage group per row")
    if not (stage == 0).any() or not (stage == 1).any():
        raise RunSortRefused("both brackets must hold rows")
    return C, stage.astype(np.int64), labels, pools


def sort_rows(coords, stage, labels, groups, mode: str = MODES[1], later_contrast: bool = False,
              normals=None) -> RunSort:
    """The step on one cloud. coords: the final coordinates, one row per position; stage: 0 for a row of the earlier
    bracket, 1 for the later; labels: each row's released label in its own file; groups: each row's lineage-map group.
    later_contrast also scores the later bracket's rows at the chosen plane (step 0's check (i)). Returns a RunSort;
    mode="as_built" returns O's order (src = 0..n-1) without computing anything."""
    if mode not in MODES:
        raise RunSortRefused(f"mode must be one of {MODES}, got {mode!r}")
    C, stage, labels, pools = _checked(coords, stage, labels, groups)
    n = len(C)
    if mode == MODES[0]:
        return RunSort(src=np.arange(n, dtype=np.int64), pairs=np.empty((0, 2), dtype=np.int64), errors=np.empty(0),
                       stage_after=stage.copy(), normal=np.full(3, np.nan), centre=C.mean(0),
                       record={"mode": mode, "rule": RULE})
    normals = fibonacci_hemisphere() if normals is None else np.asarray(normals, dtype=np.float64)
    c = C.mean(0)
    early = stage == 0
    u, contrast_e, scores = midline(C[early], labels[early], c, normals)
    contrast_l = contrast_at(C[~early], labels[~early], c, u, normals) if later_contrast else None
    s = (C - c) @ u
    plus = s >= 0  # the later bracket's side
    tallies, chunks, costs = {}, [], []
    for p in range(len(POOLS)):
        here = pools == p
        wrong_e = np.flatnonzero(here & early & plus)
        wrong_l = np.flatnonzero(here & ~early & ~plus)
        k = min(len(wrong_e), len(wrong_l))
        tally = {"wrong_earlier": int(len(wrong_e)), "wrong_later": int(len(wrong_l)), "rs_swaps": int(k)}
        if mode == MODES[1]:
            pr, cost = _assign(C, wrong_e, wrong_l, c, u)
        else:
            pr, cost, extra = _same_bracket(C, early, plus, here, k, c, u)
            tally.update(extra)
        tally["swaps"] = int(len(pr))
        tallies[POOLS[p][0]] = tally
        chunks.append(pr)
        costs.append(cost)
    pairs = np.vstack(chunks) if chunks else np.empty((0, 2), dtype=np.int64)
    errors = np.sqrt(np.concatenate(costs)) if costs else np.empty(0)
    src = np.arange(n, dtype=np.int64)
    if len(pairs):
        if len(np.unique(pairs.ravel())) != pairs.size:
            raise RunSortRefused("a row would be swapped twice")  # cannot happen: the sets are disjoint
        src[pairs[:, 0]], src[pairs[:, 1]] = pairs[:, 1], pairs[:, 0]
    after = stage[src]
    rs = sum(t["rs_swaps"] for t in tallies.values())
    rec = {"mode": mode, "rule": RULE, "normal_x": float(u[0]), "normal_y": float(u[1]), "normal_z": float(u[2]),
           "contrast_earlier": contrast_e, "best_score_earlier": float(scores.min()),
           "median_score_earlier": float(np.median(scores)),
           "contrast_later": contrast_l,
           "moved_share": float(2 * len(pairs) / n),
           "median_homology_error": float(np.median(errors)) if len(errors) else None,
           "placebo_matches_swaps": (None if mode == MODES[1] else bool(len(pairs) == rs)),
           "placebo_shortfall": (None if mode == MODES[1]
                                 else int(sum(t.get("placebo_shortfall", 0) for t in tallies.values())))}
    return RunSort(src=src, pairs=pairs, errors=errors, stage_after=after, normal=u, centre=c, record=rec,
                   pools=tallies)


def _same_bracket(C, early, plus, here, k: int, c, u) -> tuple[np.ndarray, np.ndarray, dict]:
    """The placebo in one pool: k swaps, ceil(k/2) within the earlier bracket and floor(k/2) within the later one across
    the plane (topped up from the other bracket where one lacks pairs), each bracket's the cheapest pairs of its own
    linear_sum_assignment across the plane."""
    sides = {}
    for b, rows in ((0, early), (1, ~early)):
        sides[b] = (np.flatnonzero(here & rows & plus), np.flatnonzero(here & rows & ~plus))
    cap = {b: min(len(p), len(m)) for b, (p, m) in sides.items()}
    q_e = min((k + 1) // 2, cap[0])
    q_l = min(k - q_e, cap[1])
    q_e = min(k - q_l, cap[0])
    out, cost = [], []
    for b, q in ((0, q_e), (1, q_l)):
        if q == 0:
            continue
        pr, cs = _assign(C, *sides[b], c, u)
        pick = np.argsort(cs, kind="stable")[:q]
        out.append(pr[pick])
        cost.append(cs[pick])
    pairs = np.vstack(out) if out else np.empty((0, 2), dtype=np.int64)
    costs = np.concatenate(cost) if cost else np.empty(0)
    return pairs, costs, {"placebo_earlier": int(q_e), "placebo_later": int(q_l),
                          "placebo_shortfall": int(k - q_e - q_l)}


# ------------------------------------------------------------------------------------------------ step 0's figures
def later_share(coords, later, k: int = K) -> np.ndarray:
    """N15's f: per position, the share of later-bracket rows among its k nearest positions, itself included."""
    C = np.asarray(coords, dtype=np.float64)
    nb = cKDTree(C).query(C, k=k)[1]
    return np.asarray(later, dtype=np.float64)[nb].mean(1)


def bracket_mix(coords, later, k: int = K) -> dict:
    """SD(f) (population SD, as N15) and the mean own-run share: for each position, the share of its k nearest
    positions (itself included) holding a row of its own bracket."""
    later = np.asarray(later, dtype=np.int64)
    f = later_share(coords, later, k)
    own = np.where(later == 1, f, 1.0 - f)
    return {"sd_f": float(f.std()), "own_run_share": float(own.mean())}


def reach(coords, k: int = K) -> float:
    """The median 15-nearest-cell reach: the median distance from a position to the farthest of its k nearest
    positions, itself included (vec.t2_build.neighbourhood_radius's definition)."""
    C = np.asarray(coords, dtype=np.float64)
    d = np.asarray(cKDTree(C).query(C, k=k)[0]).reshape(len(C), k)
    return float(np.median(d[:, k - 1]))
