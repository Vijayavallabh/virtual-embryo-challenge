#!/usr/bin/env python
"""Tries K12 (the T3 arm raw_r8; specification fixed by the main session before any code): recipe R8 of
experiments/research_2026_11_runs/REPORT.md (line 96) on K9's knockout-derived factor, in the simplest form the heart's
try H20 built (campaign/heart-bfix, vec.harness.B0_TREND_R8_RULE), with the recipe's own gate. A labelled reading for
a ruling, NOT the adoption rule.

R8: stack control contrasts, take a low-rank basis, build on what is left after removing it; "kill unless a left-out
control contrast lies inside the basis more than inside a random subspace". K11 left it out as "design not on record";
this arm takes H20's form: rank 1, a least-squares projection, the interior-stage residuals of the wild-type series
with recipe R2's weights, each first cleaned of its part monotone in reference level. It uses the released wild-type
series only: no replicate difference, so no fold truth.

Contrasts (per fold, from the fold's inputs only). The series is the fold's reference and the released wild-type
stages K9 adds (roles s_<stage>), in stage-day order. For every stage with a released stage on each side and every
lineage group K9's revert uses for that reference (K7b's groups: held by the reference and by at least k7.MIN_SAMPLES
other samples), with at least k7.MIN_CELLS cells at all three stages: r = y_k - (w y_(k-1) + (1 - w) y_(k+1)),
y = log(group mean proportion + k7.EPS) (k7.group_means), w = (d_(k+1) - d_k) / (d_(k+1) - d_(k-1)) (R2's weights:
0.75 and 0.25 at E8.75, 0.571 and 0.429 at E9.5). Each contrast loses its part monotone in the reference's mean
expression (harness.t3_mean_expression, the mean over its cells of log1p CP10k: H20's reference pb) by H20's
level_monotone_fit (copied below unchanged). Basis u: the rank-1 uncentred SVD direction of the cleaned contrasts.

Gate (target-free, before any arm; amended by fix brief D2 before any gate output existed): one lineage GROUP out.
For a contrast c_j of group g, u_(-g) is the rank-1 uncentred SVD direction of the cleaned contrasts of all other
groups (groups hold disjoint cells, so nothing of g enters its own basis); s_j = cos^2(c_j, u_(-g)); G = mean s_j.
With fewer than two groups that have contrasts, or a zero basis for every group, "not inside" with the reason. Null:
the same G with every left-out contrast permuted across genes within deciles of the reference's mean expression
(harness.t3_decile_permutation; draw p uses the seed [k9.PERM_SEED, GATE_STREAM, p], one permutation for every
contrast of a draw, the bases unpermuted), GATE_PERMS draws. "inside" if G is above the null's 95th percentile.
Computed for WT E8.75 (x1, r12) and WT E9.5 (x2); the arms are run only if it says "inside" for WT E8.75; otherwise the
line is closed by the recipe's own gate. A run records the gate file's commit and tree state; the reading is void if
that commit differs from the runs' or the gate's tree was dirty.

Arms (K9's stream, the same drawn cells in every arm; all through k9.arm_built's raw path):
  raw           K9's raw arm exactly.
  raw_r8        eb_factor.finalize(y - (y . u) u), y the factor's log ratio as raw_c3 takes it
                (k11.section_estimates); with no contrast or a zero basis, raw's factor bit for bit.
  raw_r8_perm   the control: the same with u permuted across genes within deciles of the reference's mean expression
                (harness.t3_decile_permutation, k9.PERM_SEED).

Reading, fixed (K9's form, the adoption rule's constants): raw_r8 is SUPPORTED if on x1 AND on r12 it beats raw by a
paired mean above 0.3 with at least 8 of 10 seeds ahead and beats raw_r8_perm the same way, and on x2 its paired mean
against raw is not below -0.3; anything else CLOSED. Void if K9's committed raw scores (RESULTS_K9) are absent or raw
does not reproduce them to 1e-6 (k11.k9_reproduction), on a missing seed or arm, or a dirty tree.

  gate  --out DIR                     the gate; writes gate.json and gate.txt (refuses to overwrite)
  run   --task x1|x2|r12 --out DIR    one task per process (refuses an existing task file; refuses unless gate.json
                                      in DIR's parent, or --gate, says "inside" for WT E8.75)
  read  --runs DIR --out DIR          apply the reading; write reading.json and reading.txt (refuses to overwrite)
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (group means, trend groups, MIN_CELLS, EPS)
import k9_residual as k9  # noqa: E402  (K9's folds, base arm, control seed and reading form)
import k11_factor_arms as k11  # noqa: E402  (raw_c3's log ratio y and the reproduction check)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.prep import load_stage  # noqa: E402

PERM_SEED = k9.PERM_SEED
STREAM = k9.STREAM
GATE_PERMS = 200  # the brief's number of null draws
GATE_STREAM = 1  # the null draws' stream: seed [PERM_SEED, GATE_STREAM, p] (fix brief D2, item 2)
GATE_REFERENCE = "E8.75"  # the gate that decides: the reference of x1 and r12, and the test's own
GATE_TASKS = {"E8.75": "x1", "E9.5": "x2"}  # a task that carries each reference and its series
INSIDE, NOT_INSIDE = "inside", "not inside"
BASE = k9.BASE  # "raw"
R8, R8_PERM = "raw_r8", "raw_r8_perm"
ARMS = (BASE, R8, R8_PERM)
CONTROLS = {R8: R8_PERM}
MIN_GAIN, MIN_AHEAD, X2_FLOOR = k9.MIN_GAIN, k9.MIN_AHEAD, k9.X2_FLOOR
TASKS = k9.TASKS
KEYS = k9.KEYS
FOLD_OF = k9.FOLD_OF
COUNTED = k9.COUNTED
RESULTS_K9 = k11.RESULTS_K9
R8_RULE = ("y - (y . u) u then eb_factor.finalize, y gene_factor's log ratio as raw_c3 takes it, u the rank-1 "
           "uncentred SVD direction of the wild-type series' interior-stage residuals per K7b group (R2's weights), "
           "each first cleaned of its part monotone in the reference's mean expression (H20's level_monotone_fit)")
R8_PERM_RULE = ("as raw_r8 with u permuted across genes within deciles of the reference's t3_mean_expression "
                "(harness.t3_decile_permutation, k9.PERM_SEED)")
_PREP: dict = {}


# ------------------------------------------------------------------------------------------------ H20's cleaning step
# pava_blocks and level_monotone_fit are copied unchanged from campaign/heart-bfix:src/vec/harness.py (try H20).
def pava_blocks(y: Sequence[float], counts: Sequence[float] | None = None) -> list[list[float]]:
    """Pool-adjacent-violators for a non-decreasing fit (a tie pools): the blocks [sum, count, first index] left
    after pooling. y holds the values (counts None: one value per position, count 1) or, with counts, each position's
    sum of values and its number of values; a block's fitted value is sum / count."""
    blocks: list[list[float]] = []  # [sum, count, first index]
    for i, v in enumerate(y):
        blocks.append([float(v), 1.0 if counts is None else float(counts[i]), float(i)])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] >= blocks[-1][0] / blocks[-1][1]:
            s, c, _ = blocks.pop()
            blocks[-1][0] += s
            blocks[-1][1] += c
    return blocks


def level_monotone_fit(v: np.ndarray, level: np.ndarray) -> tuple[np.ndarray, str]:
    """The part of v (one value per gene) monotone in `level`: genes ordered by level (genes of equal level form one
    position, so they share one fitted value), v fitted non-decreasing and non-increasing by pool-adjacent-violators
    (pava_blocks; equal gene weights), and the fit with the smaller residual sum of squares returned (the
    non-decreasing one on a tie) with its direction, "up" or "down"."""
    v = np.asarray(v, dtype=np.float64)
    level = np.asarray(level, dtype=np.float64)
    if v.shape != level.shape or v.ndim != 1:
        raise ValueError(f"v and level must be 1-d of one length; got {v.shape} and {level.shape}")
    if v.size == 0:
        return v.copy(), "up"
    order = np.argsort(level, kind="stable")
    lv, vv = level[order], v[order]
    starts = np.flatnonzero(np.r_[True, lv[1:] != lv[:-1]])
    sums = np.add.reduceat(vv, starts)
    counts = np.diff(np.r_[starts, len(vv)])

    def fit(sign: float) -> np.ndarray:
        blocks = pava_blocks(sign * sums, counts)
        edges = [int(b[2]) for b in blocks] + [len(sums)]
        g = np.empty(len(sums))
        for b, lo, hi in zip(blocks, edges[:-1], edges[1:], strict=True):
            g[lo:hi] = sign * (b[0] / b[1])
        return np.repeat(g, counts)

    up, down = fit(1.0), fit(-1.0)
    use_up = float(((vv - up) ** 2).sum()) <= float(((vv - down) ** 2).sum())
    out = np.empty_like(v)
    out[order] = up if use_up else down
    return out, "up" if use_up else "down"


# ------------------------------------------------------------------------------------------------ the contrasts
def interior_weight(d_prev: float, d: float, d_next: float) -> float:
    """R2's weight on the earlier neighbour: (d_next - d) / (d_next - d_prev); the later one gets 1 minus it."""
    return (d_next - d) / (d_next - d_prev)


def trend_groups(ref_means: dict, series_means: list[dict], series_days: list[float], ref_day: float) -> list[str]:
    """K7b's groups for a reference: held by the reference (k7.MIN_CELLS) and by enough other samples for K7's trend
    (k7.trend_at's MIN_SAMPLES), as k7.revert_factors takes them."""
    trend = k7.trend_at(series_means, series_days, ref_day)
    return [g for g in sorted(ref_means) if g in trend]


def run_contrasts(stage_means: list[dict], days: list[float], groups: Sequence[str]) -> tuple[np.ndarray, list]:
    """(contrasts as rows x genes, their (stage index, group) labels): for every interior stage of the day-ordered
    series and every group in `groups` held (k7.group_means) at that stage and both neighbours, the residual of
    log(mean + k7.EPS) against the linear interpolation in stage days of its two neighbours."""
    rows, labels = [], []
    for k in range(1, len(days) - 1):
        w = interior_weight(days[k - 1], days[k], days[k + 1])
        for g in groups:
            if not all(g in stage_means[i] for i in (k - 1, k, k + 1)):
                continue
            y = [np.log(np.asarray(stage_means[i][g], dtype=np.float64) + k7.EPS) for i in (k - 1, k, k + 1)]
            rows.append(y[1] - (w * y[0] + (1.0 - w) * y[2]))
            labels.append((k, g))
    return (np.vstack(rows) if rows else np.zeros((0, 0))), labels


def clean_contrasts(M: np.ndarray, level: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """Each contrast less its part monotone in `level` (level_monotone_fit), with the direction and the share of the
    contrast's squared norm removed."""
    M = np.array(M, dtype=np.float64)
    rec = []
    for i in range(M.shape[0]):
        fit, direction = level_monotone_fit(M[i], level)
        e0 = float((M[i] ** 2).sum())
        rec.append({"direction": direction, "removed_energy_share": float((fit ** 2).sum()) / e0 if e0 > 0 else None})
        M[i] = M[i] - fit
    return M, rec


def rank1_basis(M: np.ndarray) -> np.ndarray | None:
    """The rank-1 uncentred SVD direction of the rows of M (a unit vector over genes); None with no row or a zero
    matrix."""
    M = np.asarray(M, dtype=np.float64)
    if M.ndim != 2 or M.shape[0] == 0 or not np.any(M):
        return None
    _, sv, Vt = np.linalg.svd(M, full_matrices=False)
    return Vt[0].copy() if sv[0] > 0 else None


def _cos2(c: np.ndarray, u: np.ndarray | None) -> float:
    nc = float(np.linalg.norm(c))
    if u is None or nc == 0:
        return 0.0
    return float((c @ u) ** 2 / (nc ** 2 * float(u @ u)))


def group_out_bases(M: np.ndarray, groups: Sequence) -> list[np.ndarray | None]:
    """u_(-g) for every contrast j of group g: the rank-1 basis of the contrasts of all OTHER groups (groups hold
    disjoint cells, so nothing of group g enters its own basis)."""
    groups = list(groups)
    if len(groups) != M.shape[0]:
        raise ValueError(f"{len(groups)} group labels for {M.shape[0]} contrasts")
    lab = np.array([str(g) for g in groups], dtype=object)
    by_group = {g: rank1_basis(M[lab != g]) for g in sorted(set(lab))}
    return [by_group[g] for g in lab]


def gate_statistic(M: np.ndarray, groups: Sequence, bases: list | None = None) -> float | None:
    """G = mean over contrasts of cos^2(c_j, u_(-g)), g the group of contrast j; None with fewer than two groups."""
    if len(set(map(str, groups))) < 2:
        return None
    bases = group_out_bases(M, groups) if bases is None else bases
    return float(np.mean([_cos2(M[j], bases[j]) for j in range(M.shape[0])]))


def null_seed(p: int, seed: int = PERM_SEED) -> list[int]:
    """The seed of the null's draw p: [seed, GATE_STREAM, p], a stream no other step uses ([seed, 0] would be the
    generator of `seed` alone, the control arm's own permutation)."""
    return [int(seed), GATE_STREAM, int(p)]


def gate(M: np.ndarray, basis: np.ndarray, groups: Sequence, perms: int | None = None, seed: int = PERM_SEED
         ) -> dict:
    """The recipe's gate on cleaned contrasts M, one group out: G (gate_statistic), its null (every left-out contrast
    permuted within deciles of `basis`, the bases unpermuted; draw p with seed null_seed(p), one permutation for every
    contrast of the draw), and the verdict, "inside" if G is above the null's 95th percentile."""
    perms = GATE_PERMS if perms is None else int(perms)
    n = int(M.shape[0])
    n_groups = len(set(map(str, groups)))
    rec = {"G": None, "null_median": None, "null_p95": None, "n_contrasts": n, "n_groups_with_contrasts": n_groups,
           "n_permutations": perms, "gate_stream": GATE_STREAM}
    if n_groups < 2:
        return rec | {"verdict": NOT_INSIDE, "why": "fewer than two lineage groups with contrasts"}
    bases = group_out_bases(M, groups)
    if all(b is None for b in bases):
        return rec | {"verdict": NOT_INSIDE, "why": "a zero basis for every group"}
    G = gate_statistic(M, groups, bases)
    ar = np.arange(M.shape[1], dtype=np.float64)
    null = []
    for p in range(perms):
        idx = harness.t3_decile_permutation(ar, basis, null_seed(p, seed)).astype(np.int64)
        null.append(float(np.mean([_cos2(M[j][idx], bases[j]) for j in range(n)])))
    p95 = float(np.percentile(null, 95))
    return rec | {"G": G, "null_median": float(np.median(null)), "null_p95": p95,
                  "verdict": INSIDE if G > p95 else NOT_INSIDE}


def series_contrasts(pool, series: list, series_days: list[float], ref_day: float, genes: Sequence[str]) -> dict:
    """The fold's cleaned contrasts and basis from the reference `pool` and the other wild-type stages: the series is
    the reference plus the others, in day order; groups are K7b's for the reference."""
    for a in series:
        if list(a.var_names) != list(genes):
            raise ValueError("a wild-type series file holds other genes than the factor")
    if list(pool.var_names) != list(genes):
        raise ValueError("the reference holds other genes than the factor")
    ref_means, _ = k7.group_means(pool)
    other_means = [k7.group_means(a)[0] for a in series]
    groups = trend_groups(ref_means, other_means, list(series_days), ref_day)
    order = np.argsort(np.r_[ref_day, np.asarray(series_days, dtype=np.float64)], kind="stable")
    all_means, all_days = [ref_means, *other_means], [ref_day, *series_days]
    means, days = [all_means[i] for i in order], [float(all_days[i]) for i in order]
    if len(set(days)) != len(days):
        raise ValueError("two series stages share a day")
    M, labels = run_contrasts(means, days, groups)
    level = harness.t3_mean_expression(pool)
    if M.shape[0]:
        Mc, clean = clean_contrasts(M, level)
    else:
        Mc, clean = np.zeros((0, len(genes))), []
    return {"M": Mc, "u": rank1_basis(Mc), "level": level, "groups": groups, "days": days,
            "labels": [(days[k], g) for k, g in labels], "cleaning": clean}


# ------------------------------------------------------------------------------------------------ the factors
def remove_basis(y: np.ndarray, u: np.ndarray | None) -> np.ndarray:
    """y - (y . u) u (u a unit vector); y itself with no basis."""
    y = np.asarray(y, dtype=np.float64)
    if u is None or not np.any(u):
        return y.copy()
    return y - (y @ u) * np.asarray(u, dtype=np.float64)


def r8_factor(y: np.ndarray, den: np.ndarray, clip: tuple[float, float], u: np.ndarray | None,
              phi_raw: np.ndarray) -> np.ndarray:
    """raw_r8's factor: eb_factor.finalize(y - (y . u) u, den, clip); with no basis or a zero one, phi_raw itself (the
    arm is raw bit for bit)."""
    if u is None or not np.any(u):
        return np.array(phi_raw, dtype=np.float64)
    return harness._t3_eb().finalize(remove_basis(y, u), den, clip)


def permuted_basis(u: np.ndarray | None, basis: np.ndarray, seed: int = PERM_SEED) -> np.ndarray | None:
    """The control's basis: u permuted across genes within deciles of `basis` (harness.t3_decile_permutation)."""
    return None if u is None else harness.t3_decile_permutation(u, basis, seed)


def removed_share(y: np.ndarray, u: np.ndarray | None) -> float | None:
    """The share of y's squared norm the projection removes."""
    y = np.asarray(y, dtype=np.float64)
    e = float(y @ y)
    if u is None or e == 0:
        return None if e == 0 else 0.0
    return float((y @ u) ** 2 / e)


# ------------------------------------------------------------------------------------------------ the arms
def _series_roles(task) -> list[str]:
    return sorted(r for r in task.inputs if r.startswith(k9.SERIES_PREFIX))


def _prepared(ctx) -> dict:
    """The fold's raw factor (K9's), the basis from its wild-type series, raw_r8's factor and its control's, and the
    reference in t3_decompose's form, once per loaded fold."""
    roles = sorted(ctx.inputs)
    objs = tuple(ctx.inputs[r] for r in roles)
    hit = _PREP.get("fold")
    if hit is None or hit["key"] != ctx.task.key or len(hit["inputs"]) != len(objs) or any(
            a is not b for a, b in zip(hit["inputs"], objs, strict=True)):
        fold = harness.K4_FOLDS[FOLD_OF[ctx.task.key]]
        rules = harness._t3_head().load_rules()
        ref_role = ctx.task.reference_role
        wt, pool = ctx.inputs["wt_train"], ctx.inputs[ref_role]
        if list(wt.var_names) != list(pool.var_names):
            raise ValueError(f"{ctx.task.key}: the factor's WT and the reference hold different genes")
        frame = harness.t3_source_frame(wt, [ctx.inputs[r] for r in fold["ko_roles"]], rules, fold["source_preset"])
        phi_raw = np.array(frame.phi, dtype=np.float64)
        series = [r for r in roles if r.startswith(k9.SERIES_PREFIX)]
        sc = series_contrasts(pool, [ctx.inputs[r] for r in series], [ctx.task.role_days[r] for r in series],
                              ctx.task.role_days[ref_role], list(pool.var_names))
        u, basis = sc["u"], sc["level"]
        u_perm = permuted_basis(u, basis, PERM_SEED)
        y, _, den, _ = k11.section_estimates(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask,
                                             frame.eps)
        clip = tuple(float(c) for c in rules["components"]["factor"]["clip"])
        _PREP.clear()
        hit = _PREP["fold"] = {
            "key": ctx.task.key, "inputs": objs, "phi_raw": phi_raw, "u": u, "u_perm": u_perm,
            "phi_r8": r8_factor(y, den, clip, u, phi_raw), "phi_r8_perm": r8_factor(y, den, clip, u_perm, phi_raw),
            "clip": clip, "ref": harness.t3_reference(pool, fold["reference_preset"], rules),
            "record": {"fold": FOLD_OF[ctx.task.key], "reference": ctx.task.role_stems[ref_role],
                       "series_days": sc["days"], "groups": sc["groups"], "n_contrasts": int(sc["M"].shape[0]),
                       "basis": u is not None,
                       "level_monotone": sc["cleaning"],
                       "removed_share_y_r8": removed_share(y, u),
                       "removed_share_y_r8_perm": removed_share(y, u_perm)}}
    return hit


def arm_built(bs, ref, n: int, seed: int, f, arm: str, prep: dict):
    """One arm's cells through k9.arm_built's raw path, each with its own factor: phi_raw, phi_r8 or phi_r8_perm."""
    phi = {BASE: "phi_raw", R8: "phi_r8", R8_PERM: "phi_r8_perm"}[arm]
    return k9.arm_built(bs, ref, n, seed, f, BASE, prep[phi], None, None, prep["clip"])


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = _prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    built = arm_built(bs, ref, harness.count_for(ctx.board, pool.n_obs), ctx.seed, ref.V.p, arm, prep)
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params=harness.K4_PARAMS | {
                                     "r8": R8_RULE if arm in (R8, R8_PERM) else None,
                                     "perm_seed": PERM_SEED if arm == R8_PERM else None} | (
                                     {"r8_perm": R8_PERM_RULE} if arm == R8_PERM else {}),
                                 requires=lambda t: t.key in KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in ARMS)]


# ------------------------------------------------------------------------------------------------ the gate command
def gate_files(name: str) -> dict[str, str]:
    """The files the gate loads for a task: its reference and the wild-type series K9 adds; no knockout, no target."""
    t = TASKS[name]
    ref_role = harness.base_role(t)
    return {r: f for r, f in t.inputs.items() if r == ref_role or r.startswith(k9.SERIES_PREFIX)}


def gate_reference(inputs: dict, task) -> dict:
    """The gate for one reference from its loaded files (roles as gate_files names them)."""
    ref_role = harness.base_role(task)
    pool = inputs[ref_role]
    series = sorted(r for r in inputs if r.startswith(k9.SERIES_PREFIX))
    days = {r: harness.STAGE_DAYS[f] for r, f in task.inputs.items()}
    sc = series_contrasts(pool, [inputs[r] for r in series], [days[r] for r in series], days[ref_role],
                          list(pool.var_names))
    return gate(sc["M"], sc["level"], [g for _, g in sc["labels"]]) | {"groups": len(sc["groups"])}


def cmd_gate(a) -> None:
    out = Path(a.out)
    for p in (out / "gate.json", out / "gate.txt"):
        if p.exists():
            sys.exit(f"{p} exists")
    boards = load_boards()
    res = {}
    for stem, name in GATE_TASKS.items():
        t = TASKS[name]
        genes = boards[t.panel_board].genes
        inputs = {r: load_stage(Path(a.data) / f, genes=genes, drop_duplicates=True)
                  for r, f in gate_files(name).items()}
        res[stem] = gate_reference(inputs, t)
    doc = {"references": res, "decides": GATE_REFERENCE, "perm_seed": PERM_SEED, "gate_stream": GATE_STREAM,
           "leave_out": "one lineage group",
           "commit": harness._git("rev-parse", "HEAD"),
           "tree_clean": harness._git("status", "--porcelain", "--untracked-files=no") == ""}
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate.json").write_text(json.dumps(doc, indent=1))
    (out / "gate.txt").write_text(gate_text(doc))
    print(gate_text(doc), end="")


def gate_text(doc: dict) -> str:
    lines = [f"Tries K12 gate (R8's own, target-free; commit {doc['commit']}): the arms run only if WT "
             f"{doc['decides']} says \"{INSIDE}\""]
    for stem, g in doc["references"].items():
        lines.append(f"WT {stem}: {g['verdict']}; G {_fmt(g['G'])}, null median {_fmt(g['null_median'])}, null 95th "
                     f"percentile {_fmt(g['null_p95'])}, contrasts {g['n_contrasts']}, permutations "
                     f"{g['n_permutations']}" + (f"; {g['why']}" if g.get("why") else ""))
    return "\n".join(lines) + "\n"


def gate_says_inside(path: Path) -> tuple[bool, dict | None, dict | None]:
    """(True, the deciding reference's entry, the gate file's commit and tree state) when the gate file says "inside"
    for GATE_REFERENCE."""
    if not path.exists():
        return False, None, None
    doc = json.loads(path.read_text())
    entry = (doc.get("references") or {}).get(GATE_REFERENCE) or {}
    return entry.get("verdict") == INSIDE, entry, {"commit": doc.get("commit"), "tree_clean": doc.get("tree_clean")}


def cmd_run(a) -> None:
    out = Path(a.out)
    gpath = Path(a.gate) if a.gate else out.resolve().parent / "gate.json"
    ok, entry, gate_file = gate_says_inside(gpath)
    if not ok:
        sys.exit(f"{gpath} does not say \"{INSIDE}\" for WT {GATE_REFERENCE}: the line is closed by the recipe's own "
                 f"gate, not run and not counted")
    t = TASKS[a.task]
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{a.task}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, candidates(), seeds=SEEDS, data_dir=Path(a.data))
    res["task_source"] = "experiments/t3_decomp/k12_r8_arm.py"
    res["factor_record"] = _PREP["fold"]["record"]
    res["gate"] = entry
    res["gate_file"] = gate_file
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def decide(per: dict) -> tuple[bool, dict]:
    """raw_r8 passes only when it beats raw and raw_r8_perm on x1 and on r12, and its paired mean against raw on x2 is
    not below X2_FLOOR."""
    checks = {f"{n}: {R8} - {y}": per[n]["paired"][f"{R8} - {y}"]["passes"]
              for n in ("x1", "r12") for y in (BASE, R8_PERM)}
    checks[f"x2: {R8} - {BASE} not below {X2_FLOOR}"] = per["x2"]["paired"][f"{R8} - {BASE}"]["mean"] >= X2_FLOOR
    return all(checks.values()), checks


def reading(runs: Path, k9_dir: Path = RESULTS_K9) -> dict:
    why, per, commits, all_rows = [], {}, set(), {}
    for name, t in TASKS.items():
        p = runs / f"{name}.json"
        if not p.exists():
            why.append(f"missing task {name}")
            continue
        r = json.loads(p.read_text())
        rows = sorted(r["rows"], key=lambda x: x["seed"])
        if r["task"]["key"] != t.key or [x["seed"] for x in rows] != list(SEEDS):
            why.append(f"{name}: another task or other seeds")
            continue
        if not r["provenance"]["tree_clean"]:
            why.append(f"{name}: a dirty tree")
        gf = r.get("gate_file") or {}
        if gf.get("commit") != r["provenance"]["commit"]:
            why.append(f"{name}: the gate's commit {gf.get('commit')} is not the run's {r['provenance']['commit']}")
        if gf.get("tree_clean") is not True:
            why.append(f"{name}: the gate's tree was dirty or is not recorded")
        missing = sorted(set(ARMS) - set.intersection(*(set(q["scores"]) for q in rows)))
        if missing:
            why.append(f"{name}: missing arms {missing}")
            continue
        commits.add(r["provenance"]["commit"])
        all_rows[name] = rows
        arms = (harness.FLOOR, *ARMS)
        metrics = list(rows[0]["skills"][BASE])
        rec = r.get("factor_record") or {}
        per[name] = {
            "counted_fold": name in COUNTED,
            "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in arms},
            "paired": {f"{x} - {y}": k9._paired(rows, x, y)
                       for x, y in ((BASE, harness.FLOOR), (R8, BASE), (R8, R8_PERM), (R8_PERM, BASE))},
            "severity_states": {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                                    for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in arms},
            "skill_diff_minus_raw": {a: {m: float(np.mean([q["skills"][a][m] - q["skills"][BASE][m] for q in rows]))
                                         for m in metrics} for a in (R8, R8_PERM)},
            "removed_share_y_r8": rec.get("removed_share_y_r8"),
            "removed_share_y_r8_perm": rec.get("removed_share_y_r8_perm"),
            "gate": r.get("gate"), "record": rec}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if not Path(k9_dir).is_dir():
        why.append(f"K9's committed raw scores ({k9_dir}) are absent from this checkout")
    status, repro_why = k11.k9_reproduction(all_rows, k9_dir) if not why else ("not checked", [])
    why += repro_why
    if why:
        return {"void": True, "why": why}
    ok, checks = decide(per)
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "k9_reproduction": status,
            "checks": checks, "supported": bool(ok),
            "verdict": f"SUPPORTED: {NEXT}" if ok else "CLOSED: no request", "label": LABEL}


NEXT = ("read next on the first-choice file's own form (as K10 did) before any request; no request goes on this "
        "reading alone")
LABEL = ("A labelled reading for a ruling, NOT the adoption rule: only x1 and x2 are counted cross-knockout folds; r12 "
         "is a replicate fold that never counts; x2's reference is WT E9.5, where K7 found reverts flattered.")


def _fmt(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K12 (commit {res['commit']}): {R8} {res['verdict']}", res["label"],
             f"raw against K9's committed raw scores: {res['k9_reproduction']}",
             "Checks: " + "; ".join(f"{k} {'yes' if x else 'no'}" for k, x in res["checks"].items()),
             "Reported, deciding nothing:"]
    for n, t in TASKS.items():
        v = res["tasks"][n]
        kind = "counted cross-knockout fold" if v["counted_fold"] else "replicate fold, never counted"
        lines.append(f"{n} {FOLD_OF[t.key]} ({kind}): " + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})"
                         + ("  passes" if d["passes"] else ""))
        for a in (harness.FLOOR, *ARMS):
            lines.append(f"    {a:<12} severity " + ", ".join(f"{s} {c}" for s, c in v["severity_states"][a].items()))
        for a, sk in v["skill_diff_minus_raw"].items():
            lines.append(f"    skills {a} - raw: " + ", ".join(f"{m} {x:+.3f}" for m, x in sk.items()))
        lines.append(f"    share of y's squared norm removed: raw_r8 {_fmt(v['removed_share_y_r8'])}, raw_r8_perm "
                     f"{_fmt(v['removed_share_y_r8_perm'])}; contrasts {v['record'].get('n_contrasts')}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    for p in (out / "reading.json", out / "reading.txt"):
        if p.exists():
            sys.exit(f"{p} exists")
    res = reading(Path(a.runs))
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate")
    g.add_argument("--out", required=True)
    g.add_argument("--data", default=str(harness.DATA_DIR))
    r = sub.add_parser("run")
    r.add_argument("--task", choices=sorted(TASKS), required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--gate", default=None)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"gate": cmd_gate, "run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
