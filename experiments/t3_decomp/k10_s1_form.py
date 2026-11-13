#!/usr/bin/env python
"""Tries K10a and K10b (specification fixed by the main session on 2026-11-13, before any code): does K9's per-group
residual of the wild-type revert carry to the form of T3's first-choice file? The reading a request on that file
would rest on; NOT the adoption rule.

Why. Tries K9 (k9_residual.py) read the residual on top of a knockout factor ALONE (the harness's t3_factor_raw: no
depth, no composition). T3's first-choice file moves WT E8.75's cells by the pooled non-Mesp1 Gata4 factor at exponent
1 WITH Gata4's per-label depth step (and a composition term K10 leaves out). The frozen script's stand-in for that file
on a held-out replicate is its variant "xrep:F+depth" (t3_decompose.variant_specs: f=V.p, phi=S.phi_nm, ratio=S.ratio,
S the other replicate's full-pool estimate); build_submission.construct takes the same f, phi and ratio. K10 reads
the residual on that form.

Tasks (diagnostic pseudo-tasks, checked by vec.harness.check_task; configs/pseudo_tasks.yaml unchanged; adoption
False, never counted). Both are replicate folds on one slide:
  r12  source Gata4 replicate 1 against WT E8.75; reference WT E8.75; truth Gata4 replicate 2: K9's r12 (the
       configured fold "T3/Gata4rep1->rep2" with the released wild-type series added), under its own key
  r21  the mirror: r12 with the two replicate files (source and target) swapped, every other field kept; its
       forbidden list has the same swap applied and its truth added

Arms, all drawn from one random stream (K9's, vec.harness.K4_STREAM), so every arm edits the same drawn cells:
  s1x           the base (not a try): build_submission.construct(ref, n, rng, f=ref.V.p, phi=phi_src,
                ratio=ratio_src), then assemble. phi_src: the source replicate's non-Mesp1 factor at exponent 1
                (vec.harness.t3_source_frame's phi, as K9's raw arm takes it). ratio_src: the source replicate's
                per-label depth ratio, knockout / WT, by the lines t3_decompose.estimate computes S.ratio with
                (label_medians; rules.yaml components.depth.min_cells), on the source cells t3_source_frame calls into
                WT E8.75's vocabulary and on that WT's cells in the same vocabulary. No composition term, no stage
                mixing.
  s1x_res       K10a: the same drawn cells, each cell's factor multiplied per gene by K9's residual
                rho_g = (per_g / phi_pool) ** RES_EXPONENT of the reference's wild-type revert, phi_src x rho_g clipped
                to rules.yaml's factor clip and applied once per class with t3_decompose.apply_gene_factor (K9's
                factor_classes and combined_factor), then s1x's depth step unchanged (t3_decompose.rescale_depth with
                ratio_src per cell label), all continuing the build's own generator.
  s1x_res_perm  K10b, the control: as s1x_res with K9's decile-permuted residual (seed PERM_SEED).

Reading, fixed in the specification (the adoption rule's constants: a paired mean gain above 0.3 with at least 8 of 10
seeds ahead): CONFIRMED if on r12 AND on r21 s1x_res beats s1x that way and beats s1x_res_perm that way; then K9 and
K10 together go to the captain as one priced request. Anything else: NOT CONFIRMED, no request, and K9's reading is
reported as not carried to the first-choice file's form.

  run   --task r12|r21 --out DIR    one task per process (refuses an existing task file)
  read  --runs DIR --out DIR        apply the reading; write reading.json and reading.txt (refuses to overwrite)
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402  (the wild-type revert factors of tries K7, unchanged)
import k9_residual as k9  # noqa: E402  (the residual, the control, the per-group application, the task format)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

RES_EXPONENT = k9.RES_EXPONENT  # 0.5, K9's
PERM_SEED = k9.PERM_SEED
STREAM = k9.STREAM
BASE, MAIN, CONTROL = "s1x", "s1x_res", "s1x_res_perm"
ARMS = (BASE, MAIN, CONTROL)
MIN_GAIN, MIN_AHEAD = k9.MIN_GAIN, k9.MIN_AHEAD
SOURCE = harness.K4_FOLDS[k9.FOLDS["r12"]]  # knockout role, presets: the same for both replicate folds
DEPTH_RULE = ("t3_decompose.estimate's per-label depth ratio of the source replicate: median cell depth per label, "
              "knockout / WT, 1.0 where either side has fewer than rules.yaml components.depth.min_cells cells")
_PREP: dict = {}


# ------------------------------------------------------------------------------------------------ the tasks
def _mirror(t: harness.PseudoTask) -> harness.PseudoTask:
    """r12 with its two replicate files (source and target) swapped; the forbidden list gets the same swap and the
    truth; every other field is kept, except the board: r12's is the published validation board, whose counts are
    replicate 2's, so the mirror is scored on a synthetic board with that board's gene panel and its own files'
    counts (vec.harness.task_board's twin form, as k7_wt_revert.py's replicate-1 task)."""
    src, truth = t.inputs["ko_train"], t.target
    swap = {src: truth, truth: src}
    forbidden = tuple(dict.fromkeys([swap.get(f, f) for f in t.forbidden] + [src]))
    return harness.check_task(replace(
        t, key="T3k10/r21", target=src, inputs={r: swap.get(f, f) for r, f in t.inputs.items()}, forbidden=forbidden,
        validation_board=None, twin_board=t.validation_board or t.twin_board,
        notes="Tries K10: the mirror of fold T3/Gata4rep1->rep2 (replicate 2 against WT E8.75, truth replicate 1), "
              "with the released wild-type series added for the reference's revert."))


R12 = harness.check_task(replace(k9.TASKS["r12"], key="T3k10/r12", notes=(
    "Tries K10: fold T3/Gata4rep1->rep2 with the released wild-type series added for the reference's revert "
    "(K9's r12).")))
TASKS = {"r12": R12, "r21": _mirror(R12)}
KEYS = frozenset(t.key for t in TASKS.values())


# ------------------------------------------------------------------------------------------------ the source side
def depth_ratio(C: np.ndarray, lab: np.ndarray, wt_count: np.ndarray, wt_med: np.ndarray, L: int, dmin: int
                ) -> np.ndarray:
    """The per-label depth ratio knockout / WT exactly as t3_decompose.estimate computes Estimate.ratio (its lines,
    unchanged): the knockout's per-label cell counts and median depths (t3_decompose.label_medians) against the WT's
    (WT.count, WT.med), 1.0 where either side has fewer than `dmin` cells or the WT median is not positive."""
    t3d = harness._t3_head().t3d
    count = np.bincount(lab[lab >= 0], minlength=L)
    med = t3d.label_medians(lab, C.sum(1), L)
    ok_d = (count >= dmin) & (wt_count >= dmin) & np.isfinite(med) & (wt_med > 0)
    return np.where(ok_d, med / np.where(wt_med > 0, wt_med, 1.0), 1.0)


def source_side(wt, kos: list, rules, preset: str):
    """(phi_src, ratio_src, the WT vocabulary's labels) for a fold's source: phi_src is t3_source_frame's phi (K9's
    raw factor); ratio_src is depth_ratio on t3_source_frame's pooled source counts and labels against the WT's cells
    in the same vocabulary, prepared by the calls t3_source_frame makes (harness._shallow, add_label_column,
    harness._t3_vocab, dense_counts), since the frame does not return the WT's counts and labels."""
    t3d = harness._t3_head().t3d
    frame = harness.t3_source_frame(wt, kos, rules, preset)
    w = harness._shallow(wt)
    t3d.add_label_column(w, rules["presets"][preset]["label_mode"])
    V, wlab = harness._t3_vocab(t3d, w, rules)
    Cw = t3d.dense_counts(w)
    ratio = depth_ratio(frame.C, frame.lab, np.bincount(wlab, minlength=V.L), t3d.label_medians(wlab, Cw.sum(1), V.L),
                        V.L, int(rules["components"]["depth"]["min_cells"]))
    return np.array(frame.phi, dtype=np.float64), ratio, list(V.labels)


# ------------------------------------------------------------------------------------------------ the arms
def arm_built(bs, ref, n: int, seed: int, f, arm: str, phi_src: np.ndarray, ratio_src: np.ndarray,
              rho: dict[str, np.ndarray] | None, group: np.ndarray | None, clip: tuple[float, float]):
    """One arm's cells, before assembly. s1x: build_submission.construct with phi_src and ratio_src (the frozen
    "xrep:F+depth" build). The residual arms: the same construct call without factor or ratio (the same draws from the
    same generator), each cell's counts moved once by its combined factor (K9's classes, t3_decompose.apply_gene_factor)
    and then all cells once by ratio_src at their labels (t3_decompose.rescale_depth), continuing that generator, in
    t3_decompose.build's order; when every drawn cell's combined factor is phi_src this is s1x's build exactly."""
    rng = bs.t3d.stable_rng(seed, STREAM)
    if arm == BASE:
        return bs.construct(ref, n, rng, f=f, phi=phi_src, ratio=ratio_src)
    drawn = bs.construct(ref, n, rng, f=f)
    if not np.all(drawn.source == "wt"):
        raise RuntimeError("the residual arms expect every drawn cell to come from the reference pool")
    if len(group) != ref.pool.C.shape[0]:
        raise RuntimeError("the lineage groups do not line up with the reference pool's rows")
    C, P = ref.pool.C[drawn.row].copy(), ref.pool.P[drawn.row]
    for fac, rows in k9.factor_classes(phi_src, rho, group[drawn.row], clip):
        C[rows] = bs.t3d.apply_gene_factor(C[rows], fac, P[rows], rng)
    C = bs.t3d.rescale_depth(C, np.asarray(ratio_src)[drawn.lab], P, rng)
    return bs.Built(X=bs.t3d.lognorm(C), lab=drawn.lab, source=drawn.source, row=drawn.row)


def _prepared(ctx) -> dict:
    """The source replicate's factor and depth ratio, K7's factors of the reference, K9's residuals and the reference in
    t3_decompose's form, once per loaded fold."""
    roles = sorted(ctx.inputs)
    objs = tuple(ctx.inputs[r] for r in roles)
    hit = _PREP.get("fold")
    if hit is None or hit["key"] != ctx.task.key or len(hit["inputs"]) != len(objs) or any(
            a is not b for a, b in zip(hit["inputs"], objs, strict=True)):
        rules = harness._t3_head().load_rules()
        ref_role = ctx.task.reference_role
        wt, pool = ctx.inputs["wt_train"], ctx.inputs[ref_role]
        if list(wt.var_names) != list(pool.var_names):
            raise ValueError(f"{ctx.task.key}: the factor's WT and the reference hold different genes")
        phi_src, ratio_src, labels = source_side(wt, [ctx.inputs[r] for r in SOURCE["ko_roles"]], rules,
                                                 SOURCE["source_preset"])
        ref = harness.t3_reference(pool, SOURCE["reference_preset"], rules)
        if list(ref.V.labels) != labels:
            raise ValueError(f"{ctx.task.key}: the depth ratio's vocabulary is not the reference's")
        series = [r for r in roles if r.startswith(k9.SERIES_PREFIX)]
        phi_pool, per, sizes = k9.revert_parts(pool, [ctx.inputs[r] for r in series],
                                               [ctx.task.role_days[r] for r in series], ctx.task.role_days[ref_role])
        rho = k9.residual_factors(per, phi_pool, RES_EXPONENT)
        clip = tuple(float(c) for c in rules["components"]["factor"]["clip"])
        _PREP.clear()
        hit = _PREP["fold"] = {
            "key": ctx.task.key, "inputs": objs, "phi_src": phi_src, "ratio_src": ratio_src, "rho": rho,
            "rho_perm": k9.permuted_residuals(rho, harness.t3_mean_expression(pool), PERM_SEED), "clip": clip,
            "group": np.asarray(pool.obs["group"]).astype(str), "ref": ref,
            "record": {"source": ctx.task.role_stems[SOURCE["ko_roles"][0]],
                       "reference": ctx.task.role_stems[ref_role],
                       "series": [ctx.task.role_stems[r] for r in series], "groups_with_a_factor": len(per),
                       "labels_with_a_depth_ratio": int((ratio_src != 1.0).sum()),
                       "cosine_log_phi_raw_mean_log_rho": k9.residual_cosine(phi_src, rho, sizes)}}
    return hit


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = _prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    rho = None if arm == BASE else prep["rho" if arm == MAIN else "rho_perm"]
    built = arm_built(bs, ref, harness.count_for(ctx.board, pool.n_obs), ctx.seed, ref.V.p, arm, prep["phi_src"],
                      prep["ratio_src"], rho, prep["group"], prep["clip"])
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


PARAMS = {"factor": "the source replicate's factor: K2a's estimator against WT E8.75 (t3_source_frame): strict "
                    "non-Mesp1 labels, gene_factor (rules.yaml components.factor), exponent 1",
          "depth": DEPTH_RULE, "stream": STREAM,
          "build": "build_submission.construct and assemble with f=V.p, phi and ratio (the frozen stand-in "
                   "xrep:F+depth); no composition, no stage mixing"}


def candidates() -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params=PARAMS | {
                                     "residual": k9.RESIDUAL_RULE if arm != BASE else None,
                                     "revert_rule": k7.WT_REVERT_RULE if arm != BASE else None,
                                     "res_exponent": RES_EXPONENT if arm != BASE else None,
                                     "perm_seed": PERM_SEED if arm == CONTROL else None},
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
    res["task_source"] = "experiments/t3_decomp/k10_s1_form.py"
    res["factor_record"] = _PREP["fold"]["record"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def decide(per: dict) -> tuple[bool, dict]:
    """The fixed reading: confirmed only when s1x_res passes against s1x and against s1x_res_perm on r12 and on r21."""
    checks = {f"{n}: {MAIN} - {y}": per[n]["paired"][f"{MAIN} - {y}"]["passes"] for n in TASKS for y in (BASE, CONTROL)}
    return all(checks.values()), checks


def reading(runs: Path) -> dict:
    why, per, commits = [], {}, set()
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
        arms = (harness.FLOOR, *ARMS)
        metrics = list(rows[0]["skills"][MAIN])
        rec = r.get("factor_record") or {}
        per[name] = {
            "counted_fold": False,
            "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in arms},
            "paired": {f"{x} - {y}": k9._paired(rows, x, y) for x, y in (
                (MAIN, BASE), (MAIN, CONTROL), (CONTROL, BASE), (BASE, harness.FLOOR), (MAIN, harness.FLOOR))},
            "severity_states": {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                                    for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in arms},
            "skill_diff_s1x_res_minus_s1x": {m: float(np.mean([q["skills"][MAIN][m] - q["skills"][BASE][m]
                                                               for q in rows])) for m in metrics},
            "cosine_log_phi_raw_mean_log_rho": rec.get("cosine_log_phi_raw_mean_log_rho"),
            "groups_with_a_factor": rec.get("groups_with_a_factor"),
            "labels_with_a_depth_ratio": rec.get("labels_with_a_depth_ratio")}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if why:
        return {"void": True, "why": why}
    confirmed, checks = decide(per)
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "checks": checks,
            "k10_confirmed": bool(confirmed),
            "verdict": ("CONFIRMED: K9 and K10 together go to the captain as one priced request" if confirmed
                        else "NOT CONFIRMED: no request; K9's reading is not carried to the first-choice file's form"),
            "label": LABEL}


LABEL = ("A reading for a request, NOT the adoption rule: both folds (r12 and its mirror r21) are replicate folds on "
         "one slide, never counted by the adoption rule. The base s1x is the frozen stand-in for the first-choice file "
         "(xrep:F+depth: factor and depth from the other replicate), without the file's composition term, which is "
         "fixed by the frozen rules and whose magnitude was set from both replicates.")


def _fmt(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K10a-K10b (commit {res['commit']}): {res['verdict']}", res["label"],
             "Checks: " + "; ".join(f"{k} {'yes' if v else 'no'}" for k, v in res["checks"].items()),
             "Reported, deciding nothing:"]
    for n in TASKS:
        v = res["tasks"][n]
        lines.append(f"{n} (replicate fold, never counted): "
                     + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})"
                         + ("  passes" if d["passes"] else ""))
        for a in (harness.FLOOR, *ARMS):
            lines.append(f"    {a:<13} severity " + ", ".join(f"{s} {c}" for s, c in v["severity_states"][a].items()))
        lines.append("    skills s1x_res - s1x: " + ", ".join(f"{m} {x:+.3f}"
                                                            for m, x in v["skill_diff_s1x_res_minus_s1x"].items()))
        lines.append(f"    cosine of log(phi_src) with the weighted mean log(rho): "
                     f"{_fmt(v['cosine_log_phi_raw_mean_log_rho'])}; groups with a per-group factor "
                     f"{v['groups_with_a_factor']}; labels with a depth ratio {v['labels_with_a_depth_ratio']}")
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
