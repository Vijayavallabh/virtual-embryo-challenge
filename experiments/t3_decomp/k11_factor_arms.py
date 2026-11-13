#!/usr/bin/env python
"""Tries K11 (the T3 batch of the campaign; specification fixed by the main session, before any code): four recorded
lines read on K9's folds as arms on K9's knockout-derived factor. A labelled reading for a ruling, NOT the adoption
rule.

Folds, base, stream, permutation and reading are K9's (k9_residual.py, imported, unchanged): the tasks x1, x2 and r12
are k9.TASKS under their own keys, so the base arm `raw` is K9's own raw arm (k9.arm_built with the fold's factor
phi_raw, t3_factor_raw's build) and must reproduce K9's committed raw scores. Every arm is drawn from K9's stream
(vec.harness.K4_STREAM), so all arms hold the same drawn cells, coordinates and cell number; each differs from `raw` in
one idea only.

Arms built:
  raw           the base (not a try): K9's raw arm exactly.
  raw_did       the difference-in-differences lineage factor (experiments/research_2026_11_runs/REPORT.md R5, line 90):
                from the fold's SOURCE knockout alone (t3_source_frame's cells, labels and WT means, the frame's own
                vocabulary), D(g) = [KO - WT in L] - [KO - WT in N], each bracket gene_factor's unrescaled log ratio of
                w-weighted mean proportions (eps rules.yaml components.factor.eps), L the Mesp1-lineage labels
                (lineage map mesp1 "yes") and N the strict non-Mesp1 labels of the factor (frame.mask), each with at
                least min_cells cells on both sides. G = DID_C x D, added in log to the factor of the reference's
                Mesp1-lineage cells (reference vocabulary, mesp1 "yes"): clip(phi_raw x exp(G)); every other cell keeps
                phi_raw. Applied once per class through K9's factor_classes and arm_built. DID_C = 1: R5's
                c = K rho / (1 + (K - 1) rho) needs the correlation rho between two knockouts' contrasts, which one
                source knockout cannot give without the fold's truth; with no estimate on record the simplest
                parameter-free form, the unshrunk contrast, is used. Nothing from the fold's truth enters.
  raw_did_perm  its control: exp(G) permuted across genes within deciles of the reference's mean expression (K9's
                k9.permuted_residuals, seed k9.PERM_SEED).
  raw_c3        the section-robust source factor (experiments/research_2026_11_sota/REPORT.md lines 543-551; TRIES.md
                "section-artefact"): the same estimator's log ratio y re-estimated robustly across the source's
                sections: the per-gene median over sections of the delete-one-section jackknife pseudo-values
                K y - (K - 1) y^(-k) (eb_factor.jackknife's sections, deletions and conventions, recomputed here because
                that function returns only their variance), then gene_factor's last two steps (eb_factor.finalize:
                rescale to the WT-weighted total, clip). The median of the pseudo-values is the simplest
                parameter-free robust form; the record names no estimator. No vector is added, so it has no control.

Blocked (design not on record; BLOCKED below): raw_lvc (R4) and raw_r8 (R8).

Reading, fixed (K9's form, per arm, the adoption rule's constants: a paired mean gain above 0.3 with at least 8 of 10
seeds ahead): an arm is SUPPORTED if on x1 AND on r12 it beats raw that way, and beats its permuted control that way
where it has one, and on x2 its paired mean against raw is not below -0.3; anything else CLOSED. A supported arm is next
read on the first-choice file's own form (as K10 did) before any request. The read also checks that `raw` reproduces
K9's committed raw scores (RESULTS_K9, if that folder exists in the checkout) to 1e-6 per seed, and is void if not.

  run   --task x1|x2|r12 --out DIR    one task per process (refuses an existing task file)
  read  --runs DIR --out DIR          apply the reading; write reading.json and reading.txt (refuses to overwrite)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k9_residual as k9  # noqa: E402  (K9's folds, base arm, per-group application, control and reading form)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

DID_C = 1.0  # R5's c with one source knockout: no rho is estimable from it, so the unshrunk contrast
C3_STRENGTH = 1.0  # 1: the robust estimate; 0: the raw factor itself (the idea switched off)
PERM_SEED = k9.PERM_SEED
STREAM = k9.STREAM
MESP1 = "yes"  # the lineage map's flag for the Mesp1 lineage (L)
BASE = k9.BASE  # "raw"
DID, DID_PERM, C3 = "raw_did", "raw_did_perm", "raw_c3"
ARMS = (BASE, DID, DID_PERM, C3)
CONTROLS = {DID: DID_PERM, C3: None}  # each read arm and its permuted control (None: no vector is added)
MIN_GAIN, MIN_AHEAD, X2_FLOOR = k9.MIN_GAIN, k9.MIN_AHEAD, k9.X2_FLOOR
TASKS = k9.TASKS
KEYS = k9.KEYS
FOLD_OF = k9.FOLD_OF
COUNTED = k9.COUNTED
RESULTS_K9 = HERE / "results" / "k9_residual" / "tasks"
REPRO_TOL = 1e-6
DID_RULE = ("G = DID_C x ([KO - WT in Mesp1-lineage labels] - [KO - WT in strict non-Mesp1 labels]) from the fold's "
            "source knockout alone (gene_factor's unrescaled log ratios); clip(phi_raw x exp(G)) on the reference's "
            "Mesp1-lineage cells, phi_raw elsewhere; applied once per class with t3_decompose.apply_gene_factor")
C3_RULE = ("the per-gene median over the source's sections of the delete-one-section jackknife pseudo-values of "
           "gene_factor's log ratio (eb_factor.jackknife's sections and conventions), then eb_factor.finalize")
BLOCKED = {
    "raw_lvc": "design not on record: R4's variance is a replicate difference of the two Gata4 runs, which holds the "
               "fold's truth on x1 and r12 and does not exist for x1's one-file Mab21l2 source; R4 names no estimator "
               "for its signal and run terms on one file; the only source-side level-windowed shrinkage on record is "
               "K4's (eb_factor.py, t3_factor_eb), already a counted try.",
    "raw_r8": "design not on record: R8 gives its rank as 1-3, no method for the sparse remainder, the interior-stage "
              "residual or the removal of the part monotone in level it asks for first, and no damping strength (it "
              "removes the in-basis run part; the brief's damping of the outside part is not in the report).",
}
_PREP: dict = {}


# ------------------------------------------------------------------------------------------------ raw_did
def masked_log_ratio(mu_ko: np.ndarray, mu_wt: np.ndarray, w: np.ndarray, mask: np.ndarray, eps: float
                     ) -> np.ndarray | None:
    """gene_factor's log ratio before its rescale and clip, log((num + eps) / (den + eps)) over the labels in `mask`
    weighted by w; None when the mask holds no WT mass."""
    ww = (np.asarray(w, dtype=np.float64) * np.asarray(mask, dtype=np.float64))[:, None]
    num, den = (ww * mu_ko).sum(0), (ww * mu_wt).sum(0)
    if den.sum() <= 0:
        return None
    return np.log((num + eps) / (den + eps))


def did_contrast(mu_ko: np.ndarray, mu_wt: np.ndarray, w: np.ndarray, mask_l: np.ndarray, mask_n: np.ndarray,
                 eps: float) -> np.ndarray:
    """R5's D(g) = [KO - WT in L] - [KO - WT in N] for one knockout source; zeros where L or N holds no WT mass."""
    yl = masked_log_ratio(mu_ko, mu_wt, w, mask_l, eps)
    yn = masked_log_ratio(mu_ko, mu_wt, w, mask_n, eps)
    if yl is None or yn is None:
        return np.zeros(np.asarray(mu_wt).shape[1])
    return yl - yn


def did_factors(D: np.ndarray, c: float | None = None) -> dict[str, np.ndarray]:
    """The per-group multiplier on phi_raw: {"mesp1": exp(c D)} (c default DID_C), in K9's rho form."""
    c = DID_C if c is None else c
    return {"mesp1": np.exp(c * np.asarray(D, dtype=np.float64))}


def lineage_classes(V, lab: np.ndarray) -> np.ndarray:
    """Each pool row's class for raw_did: "mesp1" where its label's lineage flag is Mesp1 "yes", else "other"."""
    flag = np.asarray(V.mesp1).astype(str)[np.asarray(lab)] == MESP1
    return np.where(flag, "mesp1", "other")


def source_vocab_masks(wt, frame, rules, preset: str) -> tuple[np.ndarray, np.ndarray]:
    """(mask_l, mask_n) in the source frame's vocabulary: the labels with at least min_cells cells on both sides
    (t3_source_frame's `ok`) that are Mesp1 "yes" (L) and the frame's own strict non-Mesp1 mask (N), recomputed by the
    calls t3_source_frame makes (harness._shallow, add_label_column, harness._t3_vocab), since the frame does not
    return its vocabulary."""
    t3d = harness._t3_head().t3d
    w = harness._shallow(wt)
    t3d.add_label_column(w, rules["presets"][preset]["label_mode"])
    V, wlab = harness._t3_vocab(t3d, w, rules)
    mc = int(rules["components"]["factor"]["min_cells"])
    ok = (np.bincount(frame.lab[frame.lab >= 0], minlength=V.L) >= mc) & (np.bincount(wlab, minlength=V.L) >= mc)
    if not np.array_equal(ok & V.nm, np.asarray(frame.mask)):
        raise RuntimeError("the recomputed vocabulary does not give the frame's own non-Mesp1 mask")
    return ok & (np.asarray(V.mesp1).astype(str) == MESP1), np.asarray(frame.mask, dtype=bool)


# ------------------------------------------------------------------------------------------------ raw_c3
def section_estimates(C, lab, block, mu_wt, w, mask, eps):
    """eb_factor.jackknife's lines, returning (y, the delete-one-section estimates y^(-k) as sections x genes, den, K)
    instead of their variance: the same kept cells, sections, deletions and conventions (a label a deletion empties
    keeps its full-sample mean)."""
    mu_wt = np.asarray(mu_wt, dtype=np.float64)
    ww = (np.asarray(w, dtype=np.float64) * np.asarray(mask, dtype=np.float64))[:, None]
    den = (ww * mu_wt).sum(0)
    if den.sum() <= 0:
        raise ValueError("the factor's mask holds no WT mass, so there is no factor to re-estimate")
    L, G = mu_wt.shape
    lab, block = np.asarray(lab), np.asarray(block)
    labelled = lab >= 0
    keep = labelled.copy()
    keep[labelled] = np.asarray(mask, dtype=bool)[lab[labelled]]
    sections, bi = np.unique(block[keep], return_inverse=True)
    K = sections.size
    if K < 2:
        raise ValueError(f"the section-robust factor needs at least two sections holding the estimator's cells, "
                         f"found {K}")
    Ck = np.asarray(C)[keep]
    props = Ck / np.maximum(Ck.sum(1, keepdims=True), 1)
    idx = lab[keep].astype(np.int64) * K + bi
    n = np.bincount(idx, minlength=L * K).reshape(L, K)
    sums = np.asarray(sparse.csr_matrix((np.ones(idx.size), (idx, np.arange(idx.size))), shape=(L * K, idx.size))
                      @ props).reshape(L, K, G)
    n_t, s_t = n.sum(1), sums.sum(1)
    mu_full = np.divide(s_t, n_t[:, None], out=np.zeros_like(s_t), where=n_t[:, None] > 0)
    ys = np.empty((K, G))
    for k in range(K):
        rem = n_t - n[:, k]
        mu_k = np.divide(s_t - sums[:, k], rem[:, None], out=mu_full.copy(), where=rem[:, None] > 0)
        ys[k] = np.log(((ww * mu_k).sum(0) + eps) / (den + eps))
    y = np.log(((ww * mu_full).sum(0) + eps) / (den + eps))
    return y, ys, den, K


def robust_log_ratio(y: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """The per-gene median over sections of the jackknife pseudo-values K y - (K - 1) y^(-k)."""
    K = ys.shape[0]
    return np.median(K * np.asarray(y)[None, :] - (K - 1) * np.asarray(ys), axis=0)


def c3_factor(frame, strength: float | None = None) -> np.ndarray:
    """raw_c3's factor from a fold's source frame: eb_factor.finalize of the robust log ratio; at strength 0 the frame's
    own raw factor (the idea switched off); between, y moved that share of the way to the robust value."""
    s = C3_STRENGTH if strength is None else strength
    if s == 0:
        return np.array(frame.phi, dtype=np.float64)
    y, ys, den, _ = section_estimates(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask,
                                      frame.eps)
    m = robust_log_ratio(y, ys)
    return harness._t3_eb().finalize(y + s * (m - y), den, frame.clip)


# ------------------------------------------------------------------------------------------------ the arms
def _cosine(a: np.ndarray, b: np.ndarray) -> float | None:
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    return None if na == 0 or nb == 0 else float(a @ b / (na * nb))


def _prepared(ctx) -> dict:
    """The fold's raw factor (t3_source_frame's phi, K9's), raw_did's multiplier and its permutation, raw_c3's factor,
    the reference in t3_decompose's form and its rows' lineage classes, once per loaded fold."""
    roles = sorted(ctx.inputs)
    objs = tuple(ctx.inputs[r] for r in roles)
    hit = _PREP.get("fold")
    if hit is None or hit["key"] != ctx.task.key or len(hit["inputs"]) != len(objs) or any(
            a is not b for a, b in zip(hit["inputs"], objs, strict=True)):
        fold = harness.K4_FOLDS[FOLD_OF[ctx.task.key]]
        rules = harness._t3_head().load_rules()
        t3d = harness._t3_head().t3d
        wt, pool = ctx.inputs["wt_train"], ctx.inputs[ctx.task.reference_role]
        if list(wt.var_names) != list(pool.var_names):
            raise ValueError(f"{ctx.task.key}: the factor's WT and the reference hold different genes")
        frame = harness.t3_source_frame(wt, [ctx.inputs[r] for r in fold["ko_roles"]], rules, fold["source_preset"])
        phi_raw = np.array(frame.phi, dtype=np.float64)
        mask_l, mask_n = source_vocab_masks(wt, frame, rules, fold["source_preset"])
        D = did_contrast(t3d.props_type_means(frame.C, frame.lab, frame.mu_wt.shape[0]), frame.mu_wt, frame.weights,
                         mask_l, mask_n, frame.eps)
        rho = did_factors(D)
        ref = harness.t3_reference(pool, fold["reference_preset"], rules)
        phi_c3 = c3_factor(frame)
        _PREP.clear()
        hit = _PREP["fold"] = {
            "key": ctx.task.key, "inputs": objs, "phi_raw": phi_raw, "rho_did": rho,
            "rho_did_perm": k9.permuted_residuals(rho, harness.t3_mean_expression(pool), PERM_SEED),
            "phi_c3": phi_c3, "clip": tuple(float(c) for c in rules["components"]["factor"]["clip"]),
            "classes": lineage_classes(ref.V, ref.pool.lab), "ref": ref,
            "record": {"fold": FOLD_OF[ctx.task.key], "reference": ctx.task.role_stems[ctx.task.reference_role],
                       "source_mesp1_labels_in_contrast": int(mask_l.sum()),
                       "source_nonmesp1_labels": int(mask_n.sum()),
                       "cosine_log_phi_raw_did": _cosine(np.log(phi_raw), np.asarray(D)),
                       "cosine_log_phi_raw_log_phi_c3": _cosine(np.log(phi_raw), np.log(phi_c3)),
                       "sections": int(frame.n_levels)}}
    return hit


def arm_built(bs, ref, n: int, seed: int, f, arm: str, prep: dict):
    """One arm's cells, before assembly, all through k9.arm_built on K9's stream: raw with phi_raw; raw_did and its
    control with phi_raw and the per-class multiplier (the same construct draws, then one apply_gene_factor per class);
    raw_c3 as raw with its own factor in place of phi_raw."""
    if arm == BASE:
        return k9.arm_built(bs, ref, n, seed, f, BASE, prep["phi_raw"], None, None, prep["clip"])
    if arm == C3:
        return k9.arm_built(bs, ref, n, seed, f, BASE, prep["phi_c3"], None, None, prep["clip"])
    rho = prep["rho_did" if arm == DID else "rho_did_perm"]
    return k9.arm_built(bs, ref, n, seed, f, arm, prep["phi_raw"], rho, prep["classes"], prep["clip"])


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
                                     "did": DID_RULE if arm in (DID, DID_PERM) else None,
                                     "did_c": DID_C if arm in (DID, DID_PERM) else None,
                                     "perm_seed": PERM_SEED if arm == DID_PERM else None,
                                     "c3": C3_RULE if arm == C3 else None},
                                 requires=lambda t: t.key in KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in ARMS)]


def cmd_run(a) -> None:
    t = TASKS[a.task]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{a.task}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, candidates(), seeds=SEEDS, data_dir=Path(a.data))
    res["task_source"] = "experiments/t3_decomp/k11_factor_arms.py"
    res["factor_record"] = _PREP["fold"]["record"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def decide(per: dict) -> dict[str, tuple[bool, dict]]:
    """The fixed reading per arm: supported only when the arm passes against raw (and against its control where it
    has one) on x1 and on r12, and its paired mean against raw on x2 is not below X2_FLOOR."""
    out = {}
    for arm, ctrl in CONTROLS.items():
        ys = (BASE,) if ctrl is None else (BASE, ctrl)
        checks = {f"{n}: {arm} - {y}": per[n]["paired"][f"{arm} - {y}"]["passes"] for n in ("x1", "r12") for y in ys}
        checks[f"x2: {arm} - {BASE} not below {X2_FLOOR}"] = per["x2"]["paired"][f"{arm} - {BASE}"]["mean"] >= X2_FLOOR
        out[arm] = (all(checks.values()), checks)
    return out


def k9_reproduction(runs_rows: dict[str, list[dict]], k9_dir: Path = RESULTS_K9) -> tuple[str, list[str]]:
    """Checks the raw arm's per-seed scores against K9's committed raw scores (k9_dir/<task>.json) to REPRO_TOL:
    (status, reasons to void)."""
    if not k9_dir.is_dir():
        return f"not checked: {k9_dir.name} folder absent in this checkout", []
    why = []
    for name, rows in runs_rows.items():
        p = k9_dir / f"{name}.json"
        if not p.exists():
            why.append(f"{name}: K9's committed result is missing from {k9_dir}")
            continue
        theirs = {r["seed"]: r["scores"][BASE] for r in json.loads(p.read_text())["rows"]}
        ours = {r["seed"]: r["scores"][BASE] for r in rows}
        if set(theirs) != set(ours) or any(abs(ours[s] - theirs[s]) > REPRO_TOL for s in ours):
            why.append(f"{name}: raw does not reproduce K9's committed raw scores")
    return ("reproduced" if not why else "not reproduced"), why


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
        if set(ARMS) - set(rows[0]["scores"]):
            why.append(f"{name}: missing arms {sorted(set(ARMS) - set(rows[0]['scores']))}")
            continue
        commits.add(r["provenance"]["commit"])
        all_rows[name] = rows
        arms = (harness.FLOOR, *ARMS)
        metrics = list(rows[0]["skills"][BASE])
        pairs = [(BASE, harness.FLOOR)] + [(a, BASE) for a in CONTROLS] + [
            (a, c) for a, c in CONTROLS.items() if c] + [(c, BASE) for c in CONTROLS.values() if c]
        per[name] = {
            "counted_fold": name in COUNTED,
            "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in arms},
            "paired": {f"{x} - {y}": k9._paired(rows, x, y) for x, y in pairs},
            "severity_states": {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                                    for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in arms},
            "skill_diff_minus_raw": {a: {m: float(np.mean([q["skills"][a][m] - q["skills"][BASE][m] for q in rows]))
                                         for m in metrics} for a in CONTROLS},
            "record": r.get("factor_record") or {}}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    status, repro_why = k9_reproduction(all_rows, k9_dir) if not why else ("not checked", [])
    why += repro_why
    if why:
        return {"void": True, "why": why}
    decided = decide(per)
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "k9_reproduction": status,
            "arms": {a: {"supported": bool(ok), "checks": checks,
                         "verdict": (f"SUPPORTED: {NEXT}" if ok else "CLOSED: no request")}
                     for a, (ok, checks) in decided.items()},
            "blocked": BLOCKED, "label": LABEL}


NEXT = ("read next on the first-choice file's own form (as K10 did) before any request; no request goes on this "
        "reading alone")
LABEL = ("A labelled reading for a ruling, NOT the adoption rule: only x1 and x2 are counted cross-knockout folds; r12 "
         "is a replicate fold that never counts; x2's reference is WT E9.5, where K7 found reverts flattered.")


def _fmt(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K11 (commit {res['commit']}): " + "; ".join(f"{a} {v['verdict']}" for a, v in res["arms"].items()),
             res["label"], f"raw against K9's committed raw scores: {res['k9_reproduction']}"]
    for a, v in res["arms"].items():
        lines.append(f"Checks {a}: " + "; ".join(f"{k} {'yes' if x else 'no'}" for k, x in v["checks"].items()))
    lines.append("Blocked, not built: " + "; ".join(f"{a}: {r}" for a, r in res["blocked"].items()))
    lines.append("Reported, deciding nothing:")
    for n, t in TASKS.items():
        v = res["tasks"][n]
        kind = "counted cross-knockout fold" if v["counted_fold"] else "replicate fold, never counted"
        lines.append(f"{n} {FOLD_OF[t.key]} ({kind}): " + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})"
                         + ("  passes" if d["passes"] else ""))
        for a in (harness.FLOOR, *ARMS):
            lines.append(f"    {a:<13} severity " + ", ".join(f"{s} {c}" for s, c in v["severity_states"][a].items()))
        for a, sk in v["skill_diff_minus_raw"].items():
            lines.append(f"    skills {a} - raw: " + ", ".join(f"{m} {x:+.3f}" for m, x in sk.items()))
        rec = v["record"]
        lines.append(f"    cosine of log(phi_raw) with D {_fmt(rec.get('cosine_log_phi_raw_did'))}; with log(phi_c3) "
                     f"{_fmt(rec.get('cosine_log_phi_raw_log_phi_c3'))}; Mesp1-lineage labels in the contrast "
                     f"{rec.get('source_mesp1_labels_in_contrast')}")
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
