#!/usr/bin/env python
"""Tries K14 (four T3 arms: raw_half, raw_half_ctrl, raw_xt, raw_mn; specification fixed by the main session before any
code): the all-boards forms of the ledger's "Pruned before running" (within-upload mixtures over shift sizes,
exact-total kernels, negative-binomial-like resampling) read on T3's knockout folds, in the forms read on the heart
(tries H16, H18, H19). A labelled reading for a ruling, NOT the adoption rule.

Why. The round's record closed these forms on T3 with "T3's kernel is fixed by its frozen rules"; the frozen rules fix
the upload, not what may be read on the knockout folds, where factor arms were read (tries K9 to K13). So the three
forms are read on T3 as they were read elsewhere, unchanged.

Arms (K9's stream, the same drawn cells in every arm; build seed = the harness's ctx.seed):
  raw            K9's raw arm exactly (k9.arm_built's raw path with K9's phi_raw).
  raw_half       try H16's form: the residual-arm path (construct without a factor, the same draws), then phi_raw
                 applied through t3_decompose.apply_gene_factor, continuing the build's generator, to the counts of a
                 uniformly random subset of floor(n/2) of the n drawn cells (half_mask: not stratified, drawn on the
                 stream t3d.stable_rng(build seed, HALF_STREAM)); the other half keep the reference's counts as drawn.
                 With the mask of all cells it is raw bit for bit.
  raw_half_ctrl  its control, mean-matched (fix brief J2): every drawn cell moved through K9's raw path (one construct
                 call, as raw) by phi_ctrl = (1 + phi_raw) / 2 elementwise, phi_raw the fold's raw factor exactly as raw
                 uses it (K9's, after its own rescale and clip), with no finalize. raw_half's expected count of an entry
                 is exactly that of (1 + phi_raw) / 2 applied to every cell (thinning: c (1 + phi) / 2; top-up:
                 c + P n ((1 + phi) / 2 - 1)), so the control moves each gene's expected count exactly as far as the
                 mixture does, in every cell, without two kinds of cells in the file. It asks whether a gain comes from
                 having two kinds of cells in one file or only from a smaller shift. (9ec3aba's control,
                 eb_factor.finalize(y / 2), matched the mixture's mean only to first order.)
  raw_xt         try H18's form: the residual-arm path with phi_raw applied to all drawn cells in one call (raw's build
                 exactly), then per cell the total over the entries the step moves brought to its rounded expected total
                 by H18's repair (campaign/heart-cfix vec.core._repair_totals, its logic unchanged) on the stream
                 t3d.stable_rng(build seed, XT_STREAM). A cell already at its total is left as drawn. With the repair
                 switched off it is raw bit for bit.
  raw_mn         try H19's form: the same path; every moved entry of a cell redrawn together as one
                 Multinomial(rint(S), e / S), e the entries' expected counts and S their sum (campaign/heart-cfix
                 vec.core.apply_factors kernel="multinomial"), on the stream t3d.stable_rng(build seed, MN_STREAM).

Expected counts (expected_counts, derived from t3_decompose.apply_gene_factor as it stands). Per gene it draws
Binomial(c, phi) where phi < 1, c + Poisson(P n (phi - 1)) where phi > 1 (n the cell's total before the step, P the
cell's smoothed proportions) and leaves c where phi == 1; then a cell whose output is all zero keeps its counts c. On
that event every draw is 0, so E[out] = E[draw] + c Pr(all zero), Pr(all zero) the product over genes of (1 - phi)^c
(phi < 1), [c = 0] (phi == 1) and [c = 0] exp(-P n (phi - 1)) (phi > 1). The step moves an entry when its factor is not
exactly 1 and its expected count is above 0 (a top-up can move a zero count; a zero thinned stays zero); entries the
step does not move, and entries with expected count 0, are left as drawn (an entry with expected count 0 is drawn 0).
After either kernel a cell left with no counts keeps its counts c, apply_gene_factor's own rule (H18's and H19's too).

Each run prints one check line per task and arm (numbers only): share = the share of drawn cells whose counts the arm
changes relative to raw's build, mean_abs_off = for raw_xt the mean over the cells of |drawn - rint(expected)| total
over the moved entries, in counts (0 otherwise); both the mean over the seeds, and the per-seed values in the task file.

Reading, fixed (K9's form with the adoption rule's constants, as k12 and k13):
  raw_half SUPPORTED if on x1 AND on r12 it beats raw by a paired mean above 0.3 with at least 8 of 10 seeds ahead and
  beats raw_half_ctrl the same way, and on x2 its paired mean against raw is not below -0.3; anything else CLOSED.
  raw_xt, raw_mn (no gene's factor changes, so no control): SUPPORTED if on x1 AND on r12 it beats raw that way, and on
  x2 its paired mean against raw is not below -0.3; anything else CLOSED.
  Void if K9's committed raw scores (RESULTS_K9) are absent or raw does not reproduce them to 1e-6, on a missing seed,
  arm or check record, or a dirty tree. Reported, deciding nothing: raw_half_ctrl against raw, per-metric skill
  differences, the mean board score of each arm on each fold.

  run   --task x1|x2|r12 --out DIR    one task per process (refuses an existing task file)
  read  --runs DIR --out DIR          apply the reading; write reading.json and reading.txt (refuses to overwrite)
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k9_residual as k9  # noqa: E402  (K9's folds, base arm, stream and reading form)
import k11_factor_arms as k11  # noqa: E402  (raw_c3's log ratio y and the reproduction check)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

PERM_SEED = k9.PERM_SEED
STREAM = k9.STREAM
BASE = k9.BASE  # "raw"
HALF, HALF_CTRL, XT, MN = "raw_half", "raw_half_ctrl", "raw_xt", "raw_mn"
ARMS = (BASE, HALF, HALF_CTRL, XT, MN)
READ = (HALF, XT, MN)  # the arms the reading decides
CONTROLS = {HALF: HALF_CTRL, XT: None, MN: None}
# the new draws: t3d.stable_rng(build seed, name), names no other step uses (the build uses STREAM and
# f"{STREAM}:assemble"; K12's gate [PERM_SEED, 1, p]; K13's control [PERM_SEED, 2]; the decile permutations PERM_SEED)
HALF_STREAM, XT_STREAM, MN_STREAM = f"{STREAM}:k14_half", f"{STREAM}:k14_xt", f"{STREAM}:k14_mn"
MIN_GAIN, MIN_AHEAD, X2_FLOOR = k9.MIN_GAIN, k9.MIN_AHEAD, k9.X2_FLOOR
TASKS = k9.TASKS
KEYS = k9.KEYS
FOLD_OF = k9.FOLD_OF
COUNTED = k9.COUNTED
DECIDING = ("x1", "r12")
RESULTS_K9 = k11.RESULTS_K9
CHECK_KEYS = ("share", "mean_abs_off")
HALF_RULE = ("phi_raw applied through t3_decompose.apply_gene_factor (continuing K9's build generator) to a uniformly "
             "random floor(n/2) of the n drawn cells, drawn by t3d.stable_rng(build seed, HALF_STREAM); the other "
             "cells keep the reference's counts as drawn (try H16's form)")
HALF_CTRL_RULE = ("every drawn cell through K9's raw path with phi_ctrl = (1 + phi_raw) / 2 elementwise (no finalize): "
                  "it moves each gene's expected count exactly as far as raw_half's mixture does, in every cell, "
                  "without two kinds of cells in the file")
XT_RULE = ("raw's build, then per cell the total over the moved entries brought to rint of its expected total by try "
           "H18's repair (hypergeometric removal over the drawn units, multinomial addition at the expected profile) "
           "on t3d.stable_rng(build seed, XT_STREAM)")
MN_RULE = ("raw's build, then per cell the moved entries redrawn as one Multinomial(rint(S), e / S) (try H19's form) "
           "on t3d.stable_rng(build seed, MN_STREAM)")
_PREP: dict = {}
_CHECK: dict = {}  # arm -> build seed -> {"share", "mean_abs_off"}, filled by the builds of one task


# ------------------------------------------------------------------------------------------------ the half
def half_mask(n: int, seed: int) -> np.ndarray:
    """raw_half's cells: a uniformly random subset of floor(n/2) of n cells (not stratified), drawn on
    t3d.stable_rng(seed, HALF_STREAM), seed the build seed."""
    rng = harness._t3_head().t3d.stable_rng(seed, HALF_STREAM)
    out = np.zeros(n, dtype=bool)
    out[rng.choice(n, n // 2, replace=False)] = True
    return out


def half_ctrl_factor(phi_raw: np.ndarray) -> np.ndarray:
    """raw_half_ctrl's factor: (1 + phi_raw) / 2 elementwise, phi_raw the fold's raw factor exactly as raw uses it
    (K9's, after its own rescale and clip). Its expected count per entry is raw_half's: thinning gives c (1 + phi) / 2
    and top-up c + P n ((1 + phi) / 2 - 1), the mean of moving half the cells by phi and leaving the other half as
    drawn. No finalize is applied: a second rescale would break that match, and the factor lies inside the clip range
    whenever phi_raw does (the clip range holds 1, and (1 + phi) / 2 lies between 1 and phi)."""
    return (1.0 + np.asarray(phi_raw, dtype=np.float64)) / 2.0


# ------------------------------------------------------------------------------------------------ the kernels
def expected_counts(C, phi, P) -> np.ndarray:
    """E[t3_decompose.apply_gene_factor(C, phi, P, rng)] per entry (cells x genes), derived from its code (module
    docstring), the restore of an emptied cell included."""
    c = np.asarray(C, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    n = c.sum(1, keepdims=True)
    down, up = phi < 1.0, phi > 1.0
    mean = c.copy()
    zero = (c == 0.0).astype(np.float64)
    mean[:, down] = c[:, down] * phi[down]
    zero[:, down] = (1.0 - phi[down]) ** c[:, down]
    lam = np.asarray(P, dtype=np.float64)[:, up] * n * (phi[up] - 1.0)
    mean[:, up] = c[:, up] + lam
    zero[:, up] = np.where(c[:, up] == 0.0, np.exp(-lam), 0.0)
    return mean + zero.prod(1, keepdims=True) * c


def moved_entries(phi, expected: np.ndarray) -> np.ndarray:
    """The entries the step moves: factor not exactly 1 and expected count above 0."""
    return (np.asarray(phi, dtype=np.float64) != 1.0)[None, :] & (np.asarray(expected) > 0.0)


def repair_totals(out: np.ndarray, expected: np.ndarray, moved: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Try H18's repair (campaign/heart-cfix vec.core._repair_totals, its logic unchanged, on dense counts), in place on
    `out`: per cell, d = (sum of out over its moved entries) - rint(sum of expected over them); d > 0 removes d units by
    one multivariate hypergeometric draw over the moved entries' drawn units, d < 0 adds -d units by one multinomial
    draw with probabilities expected / sum(expected) over them; d == 0 draws nothing; cells in row order. Returns d."""
    w = np.where(moved, expected, 0.0)
    target = np.rint(w.sum(1)).astype(np.int64)
    have = np.where(moved, out, 0).sum(1).astype(np.int64)
    d = have - target
    for i in np.flatnonzero(d):
        idx = np.flatnonzero(moved[i])
        if d[i] > 0:
            out[i, idx] -= rng.multivariate_hypergeometric(out[i, idx].astype(np.int64), int(d[i])).astype(out.dtype)
        else:
            wi = w[i, idx]
            out[i, idx] += rng.multinomial(int(-d[i]), wi / wi.sum()).astype(out.dtype)
    return d


def redraw_multinomial(out: np.ndarray, expected: np.ndarray, moved: np.ndarray, rng: np.random.Generator) -> None:
    """Try H19's form (campaign/heart-cfix vec.core.apply_factors kernel="multinomial"), in place on `out`: per cell
    with a moved entry, in row order, its moved entries redrawn together as Multinomial(rint(S), e / S), e their
    expected counts and S their sum."""
    for i in np.flatnonzero(moved.any(1)):
        idx = np.flatnonzero(moved[i])
        e = expected[i, idx]
        s = float(e.sum())
        out[i, idx] = rng.multinomial(int(np.rint(s)), e / s).astype(out.dtype)


def restore_empty(out: np.ndarray, C: np.ndarray) -> None:
    """apply_gene_factor's rule, in place: a cell left with no counts keeps its counts C."""
    empty = out.sum(1) == 0
    out[empty] = C[empty]


def kernel_counts(t3d, C, phi, P, rng: np.random.Generator, arm: str, seed: int, repair: bool = True
                  ) -> tuple[np.ndarray, np.ndarray, dict]:
    """raw_xt's or raw_mn's counts from the drawn counts C: (the arm's counts, raw's build of the same cells, the
    check's numbers). raw's build is apply_gene_factor on all cells continuing `rng` (K9's build generator); the kernel
    draws on its own stream; repair=False returns raw's build (the idea switched off)."""
    raw = t3d.apply_gene_factor(C, phi, P, rng)
    out = raw.copy()
    off = 0.0
    if repair:
        e = expected_counts(C, phi, P)
        moved = moved_entries(phi, e)
        if arm == XT:
            d = repair_totals(out, e, moved, t3d.stable_rng(seed, XT_STREAM))
            off = float(np.abs(d).mean()) if len(d) else 0.0
        elif arm == MN:
            redraw_multinomial(out, e, moved, t3d.stable_rng(seed, MN_STREAM))
        else:
            raise ValueError(f"no kernel for arm {arm!r}")
        restore_empty(out, C)
    return out, raw, {"share": share_changed(out, raw), "mean_abs_off": off}


def share_changed(a: np.ndarray, b: np.ndarray) -> float:
    """The share of cells whose counts differ between a and b."""
    return float((np.asarray(a) != np.asarray(b)).any(1).mean()) if len(a) else 0.0


def check_line(task: str, arm: str, chk: dict) -> str:
    return f"k14 check: task={task} arm={arm} share={chk['share']:.4f} mean_abs_off={chk['mean_abs_off']:.4f}"


# ------------------------------------------------------------------------------------------------ the arms
def _drawn(bs, ref, n: int, seed: int, f):
    """The residual-arm path's draws: K9's generator, construct without a factor; (generator, drawn, C, P)."""
    rng = bs.t3d.stable_rng(seed, STREAM)
    drawn = bs.construct(ref, n, rng, f=f)
    if not np.all(drawn.source == "wt"):
        raise RuntimeError("the K14 arms expect every drawn cell to come from the reference pool")
    return rng, drawn, ref.pool.C[drawn.row].copy(), ref.pool.P[drawn.row]


def arm_built(bs, ref, n: int, seed: int, f, arm: str, prep: dict, mask: np.ndarray | None = None,
              repair: bool = True, check: dict | None = None):
    """One arm's cells before assembly. raw and raw_half_ctrl through k9.arm_built's raw path with phi_raw or
    phi_half_ctrl; raw_half, raw_xt and raw_mn through the residual-arm path (the same draws). `mask` overrides
    raw_half's cells (all True: raw bit for bit); repair=False switches raw_xt's and raw_mn's kernel off (raw bit for
    bit). When `check` is a dict the check's numbers for this build go into it."""
    t3d = bs.t3d
    if arm in (BASE, HALF_CTRL):
        phi = prep["phi_raw" if arm == BASE else "phi_half_ctrl"]
        built = k9.arm_built(bs, ref, n, seed, f, BASE, phi, None, None, prep["clip"])
        if check is not None:
            if arm == BASE:
                check.update(share=0.0, mean_abs_off=0.0)
            else:  # the same cells on the residual path, raw's factor and the control's from one generator state
                rng, _, C, P = _drawn(bs, ref, n, seed, f)
                ctrl = t3d.apply_gene_factor(C, phi, P, copy.deepcopy(rng))
                if not np.array_equal(t3d.lognorm(ctrl), built.X):
                    raise RuntimeError("the control's build on the residual path differs from K9's raw path")
                check.update(share=share_changed(ctrl, t3d.apply_gene_factor(C, prep["phi_raw"], P, rng)),
                             mean_abs_off=0.0)
        return built
    rng, drawn, C, P = _drawn(bs, ref, n, seed, f)
    phi = prep["phi_raw"]
    if arm == HALF:
        m = half_mask(len(C), seed) if mask is None else np.asarray(mask, dtype=bool)
        if m.shape != (len(C),):
            raise ValueError("the mask must hold one entry per drawn cell")
        raw = t3d.apply_gene_factor(C, phi, P, copy.deepcopy(rng)) if check is not None else None
        out = C.copy()
        out[m] = t3d.apply_gene_factor(C[m], phi, P[m], rng)
        chk = None if raw is None else {"share": share_changed(out, raw), "mean_abs_off": 0.0}
    elif arm in (XT, MN):
        out, _, chk = kernel_counts(t3d, C, phi, P, rng, arm, seed, repair)
    else:
        raise ValueError(f"unknown arm {arm!r}")
    if check is not None:
        check.update(chk)
    return bs.Built(X=t3d.lognorm(out), lab=drawn.lab, source=drawn.source, row=drawn.row)


def _prepared(ctx) -> dict:
    """The fold's raw factor (K9's), raw_half_ctrl's factor ((1 + phi_raw) / 2) and the reference in t3_decompose's
    form, once per loaded fold."""
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
        y, _, den, _ = k11.section_estimates(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask,
                                             frame.eps)
        clip = tuple(float(c) for c in rules["components"]["factor"]["clip"])
        phi_ctrl = half_ctrl_factor(phi_raw)
        full = harness._t3_eb().finalize(y, den, clip)
        _PREP.clear()
        hit = _PREP["fold"] = {
            "key": ctx.task.key, "inputs": objs, "phi_raw": phi_raw, "phi_half_ctrl": phi_ctrl, "clip": clip,
            "ref": harness.t3_reference(pool, fold["reference_preset"], rules),
            "record": {"fold": FOLD_OF[ctx.task.key], "reference": ctx.task.role_stems[ref_role],
                       # how far finalize(y) is from phi_raw (kept from 9ec3aba as a record; the control no longer
                       # uses y)
                       "max_abs_log_finalize_y_over_phi_raw": float(np.max(np.abs(np.log(full / phi_raw)))),
                       "streams": {HALF: HALF_STREAM, XT: XT_STREAM, MN: MN_STREAM}}}
    return hit


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = _prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    chk: dict = {}
    built = arm_built(bs, ref, harness.count_for(ctx.board, pool.n_obs), ctx.seed, ref.V.p, arm, prep, check=chk)
    _CHECK.setdefault(arm, {})[int(ctx.seed)] = chk
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    rules = {HALF: HALF_RULE, HALF_CTRL: HALF_CTRL_RULE, XT: XT_RULE, MN: MN_RULE}
    streams = {HALF: HALF_STREAM, XT: XT_STREAM, MN: MN_STREAM}

    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params=harness.K4_PARAMS | {"k14_form": rules.get(arm),
                                                             "k14_stream": streams.get(arm)},
                                 requires=lambda t: t.key in KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in ARMS)]


def task_check(per_seed: dict) -> dict:
    """The task file's check record: per arm the mean over the seeds of each number, and the per-seed values."""
    return {arm: {**{k: float(np.mean([v[k] for v in per_seed[arm].values()])) for k in CHECK_KEYS},
                  "per_seed": {str(s): v for s, v in sorted(per_seed[arm].items())}}
            for arm in ARMS if per_seed.get(arm)}


def cmd_run(a) -> None:
    t = TASKS[a.task]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{a.task}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    _CHECK.clear()
    res = harness.run_task(t, candidates(), seeds=SEEDS, data_dir=Path(a.data))
    chk = task_check(_CHECK)
    res["task_source"] = "experiments/t3_decomp/k14_form_arms.py"
    res["factor_record"] = _PREP["fold"]["record"]
    res["k14_check"] = chk
    path.write_text(json.dumps(res, indent=1, default=float))
    for arm in ARMS:
        if arm in chk:
            print(check_line(a.task, arm, chk[arm]), flush=True)
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def decide(per: dict) -> dict[str, tuple[bool, dict]]:
    """Per read arm: passes only when it beats raw (and, for raw_half, raw_half_ctrl) on x1 and on r12, and its paired
    mean against raw on x2 is not below X2_FLOOR."""
    out = {}
    for arm in READ:
        against = (BASE,) if CONTROLS[arm] is None else (BASE, CONTROLS[arm])
        checks = {f"{n}: {arm} - {y}": per[n]["paired"][f"{arm} - {y}"]["passes"] for n in DECIDING for y in against}
        checks[f"x2: {arm} - {BASE} not below {X2_FLOOR}"] = per["x2"]["paired"][f"{arm} - {BASE}"]["mean"] >= X2_FLOOR
        out[arm] = (all(checks.values()), checks)
    return out


PAIRS = ((BASE, harness.FLOOR), (HALF, BASE), (HALF, HALF_CTRL), (HALF_CTRL, BASE), (XT, BASE), (MN, BASE))


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
        missing = sorted(set(ARMS) - set.intersection(*(set(q["scores"]) for q in rows)))
        if missing:
            why.append(f"{name}: missing arms {missing}")
            continue
        chk = r.get("k14_check")
        if not isinstance(chk, dict) or any(not isinstance(chk.get(a), dict) or not set(CHECK_KEYS) <= set(chk[a])
                                            for a in ARMS):
            why.append(f"{name}: the k14 check record is missing")
            continue
        commits.add(r["provenance"]["commit"])
        all_rows[name] = rows
        arms = (harness.FLOOR, *ARMS)
        metrics = list(rows[0]["skills"][BASE])
        per[name] = {
            "counted_fold": name in COUNTED,
            "k14_check": {a: {k: chk[a][k] for k in CHECK_KEYS} for a in ARMS},
            "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in arms},
            "paired": {f"{x} - {y}": k9._paired(rows, x, y) for x, y in PAIRS},
            "severity_states": {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                                    for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in arms},
            "skill_diff_minus_raw": {a: {m: float(np.mean([q["skills"][a][m] - q["skills"][BASE][m] for q in rows]))
                                         for m in metrics} for a in ARMS if a != BASE},
            "record": r.get("factor_record") or {}}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if not Path(k9_dir).is_dir():
        why.append(f"K9's committed raw scores ({k9_dir}) are absent from this checkout")
    status, repro_why = k11.k9_reproduction(all_rows, k9_dir) if not why else ("not checked", [])
    why += repro_why
    if why:
        return {"void": True, "why": why}
    decided = decide(per)
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "k9_reproduction": status,
            "arms": {a: {"supported": bool(ok), "checks": checks,
                         "verdict": f"SUPPORTED: {NEXT}" if ok else "CLOSED: no request"}
                     for a, (ok, checks) in decided.items()},
            "label": LABEL, "control": {HALF_CTRL: HALF_CTRL_RULE}}


NEXT = ("read next on the first-choice file's own form (as K10 did) before any request; no request goes on this "
        "reading alone")
LABEL = ("A labelled reading for a ruling, NOT the adoption rule: only x1 and x2 are counted cross-knockout folds; r12 "
         "is a replicate fold that never counts; x2's reference is WT E9.5, where K7 found reverts flattered.")


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K14 (commit {res['commit']}):"]
    lines += [f"  {a} {v['verdict']}" for a, v in res["arms"].items()]
    lines += [res["label"], f"{HALF_CTRL}, the control of {HALF}: {HALF_CTRL_RULE}.",
              f"raw against K9's committed raw scores: {res['k9_reproduction']}"]
    for a, v in res["arms"].items():
        lines.append(f"Checks {a}: " + "; ".join(f"{k} {'yes' if x else 'no'}" for k, x in v["checks"].items()))
    lines.append("Reported, deciding nothing:")
    for n, t in TASKS.items():
        v = res["tasks"][n]
        kind = "counted cross-knockout fold" if v["counted_fold"] else "replicate fold, never counted"
        lines.append(f"{n} {FOLD_OF[t.key]} ({kind}): " + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for a in ARMS:
            lines.append("    " + check_line(n, a, v["k14_check"][a]))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})"
                         + ("  passes" if d["passes"] else ""))
        for a in (harness.FLOOR, *ARMS):
            lines.append(f"    {a:<14} severity " + ", ".join(f"{s} {c}" for s, c in v["severity_states"][a].items()))
        for a, sk in v["skill_diff_minus_raw"].items():
            lines.append(f"    skills {a} - raw: " + ", ".join(f"{m} {x:+.3f}" for m, x in sk.items()))
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
    r = sub.add_parser("run")
    r.add_argument("--task", choices=sorted(TASKS), required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
