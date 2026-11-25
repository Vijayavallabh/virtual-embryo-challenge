"""The heart's trend-form line: four directions of the banked heart step (row R55 of experiments/harness/TRIES.md, the
captain's order of 2026-11-24; the HEART PORTFOLIO design of the step-change round as its cold challenge fixes it,
OKS125 struck; the screen is row H35, experiments/harness/h35_heart_trend_forms.py).

The banked step is row H28's strength-scaled b0_trend step as experiments/builds/build_test.py's `_trend_strength_step`
runs it on file D's context (OTS125 at strength 1.25): the draw (`vec.harness.draw_copy`), the B0 typing (obs["merged"],
Untracked tagged by stage), the trend factors (`vec.harness.b0_trend_factors`) times the strength, momentum_shift's
fallback damping B0_DAMPING times the strength, the stream [seed, 3]. Each form changes ONE thing, how each lineage's
log factor is formed, and keeps the draw, the rows, the coordinates, the stream and the damping:

  category  OTK125  the merged Endothelial lineage split by its chamber category (configs/t2_chamber_classes.yaml, the
                    cells' own labels), the B0 types Endothelial/endocardium and Endothelial/vascular_endothelium at
                    every stage, each with its own slope; every other cell keeps its merged type, so CM, epicardium
                    (merged null) and every Untracked type are as in OTS125
  contrast  OTX125  each lineage's trend factor minus s-bar, the unweighted per-gene mean of the trend factors fitted on
                    the context whose reference input is the whole reference stage, over that stage's lineages (row
                    H12's contrast); each moved lineage's L2 norm held at the banked factor's
  whole     OTW125  the trend factors fitted on the context whose reference input is the whole reference stage; the
                    draw, the typing and momentum_shift stay on D's pool
  portion   OTP125  each lineage's slope over its usable input stages before the reference (the reference's own point
                    out), at the banked budget; its L2 norm held at the banked factor's; a lineage with fewer than two
                    such stages keeps the banked factor

A lineage the trend mapping leaves out (B0's own factor through momentum_shift's damping) is as in OTS125 in every
form. `zero=True` takes each form's change out (category: no split; contrast: s-bar = 0; whole: the context's own
reference in place of the whole stage; portion: the reference's point kept), and the step is then the banked step bit
for bit (tests/test_heart_trend_line_h.py). Nothing here reads a file, a target or a score; the records hold no cell
number (the momentum_shift report's counts are left out)."""

from __future__ import annotations

from dataclasses import dataclass, replace

import anndata as ad
import numpy as np

from vec import harness, t1_scope, t2_build
from vec.core import counts_of, mean_cp10k
from vec.io import dense

FORMS = ("category", "contrast", "whole", "portion")
BANKED = "banked"  # the banked step itself (OTS125's), the forms' control
FILES = {"category": "OTK125", "contrast": "OTX125", "whole": "OTW125", "portion": "OTP125"}
TAGS = {"category": "tk", "contrast": "tx", "whole": "tw", "portion": "tp"}  # files T2heart_organ_ts125_<tag>
SPLIT_LINEAGE = "Endothelial"  # the one merged lineage the category form splits (the challenge's fix F3)
SPLIT_KEY = "b0_split"  # the obs column the category form's trend factors are keyed by at every stage
NORM_HELD = ("contrast", "portion")  # forms whose moved lineages keep the banked factor's L2 norm
WHOLE_FORMS = ("contrast", "whole")  # forms that read the whole reference stage
RULES = {
    "category": "the merged Endothelial lineage typed Endothelial/<chamber category> at every stage "
                "(configs/t2_chamber_classes.yaml, the cells' own labels), each part with its own trend slope; every "
                "other cell keeps its merged type",
    "contrast": "each lineage's trend factor minus the unweighted per-gene mean of the trend factors fitted on the "
                "whole reference stage's context over that stage's lineages, rescaled to the banked factor's L2 norm",
    "whole": "the trend factors fitted on the context whose reference input is the whole reference stage; the draw, "
             "the typing and momentum_shift on D's pool",
    "portion": "each lineage's least-squares slope over its usable input stages before the reference, at the banked "
               "budget, rescaled to the banked factor's L2 norm; fewer than two such stages: the banked factor",
}


@dataclass
class StepArgs:
    """Everything a form hands vec.t1_scope.momentum_shift, and the draw it came from."""

    pool: ad.AnnData  # the context's reference input (D's pool), the draw's source
    rows: np.ndarray  # the drawn rows of the pool
    coords: np.ndarray | None  # the drawn cells' coordinates (before any size step)
    start: ad.AnnData  # `start`, typed
    end: ad.AnnData  # the pool, typed
    like: ad.AnnData  # the drawn cells, typed
    damping: float  # momentum_shift's fallback damping: B0_DAMPING x the strength
    mapping: dict[str, np.ndarray]  # lineage -> log factor at the strength (logf_by_type)
    record: dict  # target-free, no cell number


def merged_types(a: ad.AnnData, tag: str | None) -> np.ndarray:
    """B0's typing: obs["merged"], Untracked tagged by stage (`tag`; None leaves it plain, as the factors' key)."""
    m = np.asarray(a.obs[harness.T2_LINEAGE_KEY]).astype(str)
    return m if tag is None else np.where(m == harness.T2_UNTRACKED, f"{harness.T2_UNTRACKED}@{tag}", m)


def split_types(a: ad.AnnData, tag: str | None, spec: t2_build.T2ChamberSpec | None = None) -> np.ndarray:
    """The category form's typing: merged_types with each Endothelial cell typed Endothelial/<its chamber category>
    (vec.t2_build.t2_chamber_categories of its label: an unmapped label refuses, UnmappedT2LabelError). Every other cell
    keeps its merged type."""
    t = merged_types(a, tag).astype(object)
    endo = np.flatnonzero(t == SPLIT_LINEAGE)
    if len(endo):
        if t2_build.T2_LABEL_KEY not in a.obs:
            raise ValueError(f"a stage has no obs[{t2_build.T2_LABEL_KEY!r}]: the category form needs the labels")
        labels = np.asarray(a.obs[t2_build.T2_LABEL_KEY]).astype(str)[endo]
        cats = t2_build.t2_chamber_categories(labels, spec)
        t[endo] = [f"{SPLIT_LINEAGE}/{c}" for c in cats]
    return t.astype(str)


def _typed(a: ad.AnnData, types: np.ndarray) -> ad.AnnData:
    return harness._shallow(a, b0_type=types)


def _types_of(like: ad.AnnData) -> list[str]:
    return sorted(set(np.asarray(like.obs["b0_type"]).astype(str)))


def _norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def _cosine(x: np.ndarray, y: np.ndarray) -> float:
    return harness._cosine(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))


def _held(v: np.ndarray, c: np.ndarray, what: str) -> np.ndarray:
    """`c` rescaled to v's L2 norm (the banked factor's), in this order of operations so that c == v gives v back bit
    for bit: c x (||v|| / ||c||). A direction of zero or non-finite length cannot carry the norm: refused."""
    nv, nc = _norm(v), _norm(c)
    if not (np.isfinite(nc) and nc > 0 and np.isfinite(nv)):
        raise ValueError(f"{what}: the form's direction has no finite, non-zero length; its norm cannot be held")
    return c * (nv / nc)


def _banked(ctx: harness.BuildContext, strength: float) -> tuple:
    """The banked step's own arguments, line by line as build_test.py's `_trend_strength_step` forms them."""
    harness._b0_trend_task_check(ctx.task)
    s = float(strength)
    pool, rows, _, coords = harness.draw_copy(ctx)
    start, end = ctx.inputs["start"], pool
    like = _typed(pool, merged_types(pool, "end"))[rows]
    genes = list(like.var_names)
    logf, rec = harness.b0_trend_factors(ctx, genes, _types_of(like))
    return s, pool, rows, coords, start, end, like, genes, logf, rec


def _whole_context(ctx: harness.BuildContext, whole: ad.AnnData | None, form: str) -> harness.BuildContext:
    if whole is None:
        raise ValueError(f"the {form} form needs the whole reference stage")
    return replace(ctx, inputs={**ctx.inputs, ctx.task.reference_role: whole})


def form_arguments(ctx: harness.BuildContext, strength: float, form: str | None = None,
                   whole: ad.AnnData | None = None, zero: bool = False,
                   spec: t2_build.T2ChamberSpec | None = None) -> StepArgs:
    """What the step hands momentum_shift under `form` (None or BANKED: the banked step) at `strength`, on the context's
    draw. `whole` is the whole reference stage (contrast and whole only; the context's reference input is D's pool).
    `zero` takes the form's change out (the banked step bit for bit). Raises on a form that cannot be formed (no
    whole-stage trend lineage for the contrast, a direction of zero length where a norm is held, an unmapped label)."""
    form = BANKED if form is None else form
    if form != BANKED and form not in FORMS:
        raise ValueError(f"unknown trend form {form!r}: one of {FORMS}")
    s, pool, rows, coords, start, end, like, genes, logf, rec = _banked(ctx, strength)
    banked_start, banked_end = _typed(start, merged_types(start, "start")), _typed(end, merged_types(end, "end"))
    record: dict = {"form": form, "strength": s, "zero": bool(zero), "budget_days": float(rec["budget_days"]),
                    "rule": RULES.get(form, "the banked step"), "trend": rec}
    ref_role = ctx.task.reference_role
    if form == BANKED:
        mapping = {k: s * v for k, v in logf.items()}
        return StepArgs(pool, rows, coords, banked_start, banked_end, like, harness.B0_DAMPING * s, mapping, record)
    if form == "category":
        spec = None if zero else (spec or t2_build.load_t2_chamber_spec())

        def kind(a, tag):
            return merged_types(a, tag) if zero else split_types(a, tag, spec)
        like_k = _typed(pool, kind(pool, "end"))[rows]
        keyed = replace(ctx, inputs={r: harness._shallow(a, **{SPLIT_KEY: kind(a, None)})
                                     for r, a in ctx.inputs.items()})
        logf_k, rec_k = harness.b0_trend_factors(keyed, genes, _types_of(like_k), key=SPLIT_KEY)
        mapping = {k: s * v for k, v in logf_k.items()}
        record.update(trend=rec_k, split_lineage=SPLIT_LINEAGE,
                      split_types=sorted(t for t in _types_of(like_k) if t.startswith(SPLIT_LINEAGE + "/")),
                      collapsed=collapsed(ctx, rows, kind))
        return StepArgs(pool, rows, coords, _typed(start, kind(start, "start")), _typed(end, kind(end, "end")), like_k,
                        harness.B0_DAMPING * s, mapping, record)
    mapping, per = {}, {}
    if form == "whole":
        wctx = ctx if zero else _whole_context(ctx, whole, form)
        logf_w, rec_w = harness.b0_trend_factors(wctx, genes, _types_of(like))
        mapping = {k: s * v for k, v in logf_w.items()}
        for k in sorted(set(logf) | set(logf_w)):
            per[k] = ({"cosine": _cosine(logf_w[k], logf[k]), "norm_ratio": _norm(logf_w[k]) / _norm(logf[k])}
                      if k in logf and k in logf_w and _norm(logf[k]) > 0
                      else {"in_banked": k in logf, "in_whole": k in logf_w})
        record.update(trend_whole=rec_w, lineages=per)
    elif form == "contrast":
        if zero:
            sbar, shared = np.zeros(len(genes)), []
        else:
            wctx = _whole_context(ctx, whole, form)
            wlike = _typed(whole, merged_types(whole, "end"))
            logf_w, _ = harness.b0_trend_factors(wctx, genes, _types_of(wlike))
            if not logf_w:
                raise ValueError("the contrast form: no lineage of the whole reference stage has a trend factor")
            shared = sorted(logf_w)
            sbar = np.mean(np.vstack([np.asarray(logf_w[k], dtype=np.float64) for k in shared]), axis=0)
        for k in sorted(logf):
            v = logf[k]
            c = v - sbar
            per[k] = {"cosine": _cosine(c, v), "norm_ratio_before_hold": _norm(c) / _norm(v) if _norm(v) > 0
                      else float("nan")}
            mapping[k] = s * _held(v, c, f"the contrast form, lineage {k}")
        record.update(shared_lineages=shared, shared_norm=_norm(sbar), lineages=per)
    else:  # portion
        eps = float(harness._momentum_default("eps"))
        budget = float(rec["budget_days"])
        by_lineage = {r["lineage"]: r for r in rec["lineages"]}
        merged = {r: np.asarray(a.obs[harness.T2_LINEAGE_KEY]).astype(str) for r, a in ctx.inputs.items()}
        counts: dict = {}
        kept = []
        for k in sorted(logf):
            usable = list(zip(by_lineage[k]["stages"], by_lineage[k]["days"], strict=True))
            used = usable if zero else [(r, d) for r, d in usable if r != ref_role]
            if len(used) < 2:
                mapping[k] = s * logf[k]
                kept.append(k)
                continue
            for r, _ in used:
                if r not in counts:
                    counts[r] = counts_of(ctx.inputs[r], genes)
            y = [np.log(mean_cp10k(counts[r][np.flatnonzero(merged[r] == k)]) + eps) for r, _ in used]
            p = harness.trend_log_factor(y, [d for _, d in used], budget)
            v = logf[k]
            per[k] = {"stages": [r for r, _ in used], "cosine": _cosine(p, v),
                      "norm_ratio_before_hold": _norm(p) / _norm(v) if _norm(v) > 0 else float("nan")}
            mapping[k] = s * _held(v, p, f"the portion form, lineage {k}")
        record.update(lineages=per, banked_factor_kept=kept)
    return StepArgs(pool, rows, coords, banked_start, banked_end, like, harness.B0_DAMPING * s, mapping, record)


def collapsed(ctx: harness.BuildContext, rows: np.ndarray, kind=None) -> dict[str, bool]:
    """Per input stage, and for the draw: whether the split lineage's cells there carry at most one chamber category
    (the category form then types them as the merged lineage does, up to the name). True everywhere is the form's
    fall-through (G0 of row H35). Booleans only: no cell number."""
    kind = kind or (lambda a, tag: split_types(a, tag))
    out = {}
    for r, a in ctx.inputs.items():
        t = kind(a, None)
        out[r] = len({x for x in t if x.startswith(SPLIT_LINEAGE + "/")}) <= 1
    t = kind(ctx.inputs[ctx.task.reference_role], None)[np.asarray(rows)]
    out["draw"] = len({x for x in t if x.startswith(SPLIT_LINEAGE + "/")}) <= 1
    return out


def shift(ctx: harness.BuildContext, args: StepArgs):
    """momentum_shift on the form's arguments, as the banked step calls it: key b0_type, no held-flat gene, the
    stream [seed, 3]."""
    return t1_scope.momentum_shift(args.start, args.end, args.like, args.damping, (), key="b0_type",
                                   seed=[ctx.seed, 3], logf_by_type=args.mapping)


def statuses(report) -> list[dict]:
    """momentum_shift's per-type outcome, the type and its status only (its counts are left out)."""
    return [{"type": str(r["type"]), "status": str(r["status"])} for r in report.rows]


def form_step(ctx: harness.BuildContext, strength: float, form: str | None = None, whole: ad.AnnData | None = None,
              zero: bool = False, spec: t2_build.T2ChamberSpec | None = None) -> tuple:
    """The step under `form` through `harness.finish`, as the banked step ends (build_test.py `_trend_strength_step`):
    returns the stage, the shifted counts, the draw's rows and the record (target-free; per type its status, shifted
    or skipped, and no cell number)."""
    args = form_arguments(ctx, strength, form, whole, zero, spec)
    shifted, report = shift(ctx, args)
    stage = harness.finish(ctx, args.pool, dense(shifted[:, ctx.board.genes].X), args.coords)
    record = {**args.record, "types": statuses(report),
              "shifted": sorted(r["type"] for r in report.rows if r["status"] == "shifted")}
    return stage, counts_of(shifted, list(ctx.board.genes)), args.rows, record


# ------------------------------------------------------------------------------------------------ G0 of row H35
def applied(ctx: harness.BuildContext, args: StepArgs) -> dict:
    """The log factors momentum_shift applies under `args`, taken at its call to apply_factors, which hands the counts
    back as they are (no count is drawn and the stream is not touched): the types in momentum_shift's order, their
    log-factor rows (before each cell's depth constant), each drawn cell's row, and each type's status from
    momentum_shift's own report. vec.t1_scope is as it was afterwards, whatever happens in between."""
    from scipy import sparse

    seen: list = []
    own = t1_scope.apply_factors

    def take(C, F, tidx, held=None, keep_depth=True, seed=0):
        seen.append((np.log(np.asarray(F, dtype=np.float64)), np.asarray(tidx).copy()))
        n = C.shape[0]
        return sparse.csr_matrix(C, dtype=np.int64, copy=True), np.ones(n), np.zeros(n, dtype=bool)

    t1_scope.apply_factors = take
    try:
        _, report = shift(ctx, args)
    finally:
        t1_scope.apply_factors = own
    if len(seen) != 1:
        raise RuntimeError("momentum_shift did not reach apply_factors exactly once")
    logF, tidx = seen[0]
    types = sorted(set(np.asarray(args.like.obs["b0_type"]).astype(str)))
    return {"types": types, "logF": logF, "tidx": tidx, "statuses": statuses(report)}


def g0_statistics(form: dict, banked: dict) -> dict:
    """Row H35's G0 figures for one form against the banked step, over the drawn cells (the same cells in the same
    order in both): the changed share of the step, sum over cells of ||lf_form - lf_banked||^2 over the sum of
    ||lf_banked||^2 (lf a cell's applied log-factor vector), and the mean |log f| over cells and genes of each. Ratios
    and means only: no cell number."""
    tf, tb = np.asarray(form["tidx"]), np.asarray(banked["tidx"])
    if tf.shape != tb.shape:
        raise ValueError("the form and the banked step drew different cells")
    Lf, Lb = form["logF"], banked["logF"]
    pairs, n = np.unique(np.stack([tf, tb], axis=1), axis=0, return_counts=True)
    num = float(sum(c * float(((Lf[i] - Lb[j]) ** 2).sum()) for (i, j), c in zip(pairs, n, strict=True)))
    den = float(sum(c * float((Lb[j] ** 2).sum()) for (_, j), c in zip(pairs, n, strict=True)))
    cells = float(n.sum())
    mean_f = float(sum(c * float(np.abs(Lf[i]).mean()) for (i, _), c in zip(pairs, n, strict=True)) / cells)
    mean_b = float(sum(c * float(np.abs(Lb[j]).mean()) for (_, j), c in zip(pairs, n, strict=True)) / cells)
    # a banked step that moves nothing: a form that moves nothing either changes none of it, any other all of it
    share = num / den if den > 0 else (0.0 if num == 0 else float("inf"))
    wild = mean_f / mean_b if mean_b > 0 else (0.0 if mean_f == 0 else float("inf"))
    return {"changed_share": share, "mean_abs_logf": mean_f, "mean_abs_logf_banked": mean_b, "wild_ratio": wild}
