#!/usr/bin/env python
"""Tries K15a and K15b (brief AO; specification fixed by the main session before any code): Task 3's knockout factor at
half size, (1 + factor) / 2 per gene (try K14b's form), read in the first-choice file's own form, with its per-label
depth step (try K10's stand-in s1x). A labelled reading for a ruling, NOT the adoption rule. It changes no file: the
first-choice file changes only by the captain's own ruling (decision log 2026-11-18 (c) fixes the reading's use).

Why. Try K14b (a control) read the half factor on the factor alone, with no depth step; row K14 fixed that a supported
arm is next read on the first-choice file's own form (as tries K10 did) before any request. This is that read.

Folds. r12 and r21 exactly as try K10 defines them (k10_s1_form.TASKS, the same keys): the two Gata4 replicate
directions against WT E8.75. The cross-knockout folds x1 and x2 of try K9 are NOT built: design not on record. Their
source's wild type is another stage than the fold's reference (x1: Mab21l2 against WT E9.5 in the mab21l2 preset,
reference WT E8.75 in the gata4 preset; x2: Gata4 against WT E8.75, reference WT E9.5 in the mab21l2 preset), so try
K10's per-label depth ratio comes out in the source WT's label vocabulary, not the reference's, and k10_s1_form refuses
that case. The record fixes no rule for taking it across: which labels correspond between the two stages' vocabularies
(by name or otherwise), what ratio a reference label without a counterpart takes, and whether the ratio is the source's
knockout / its own WT (a within-stage depth change) or is also to carry the depth difference between the source's WT
stage and the reference's. `run --task x1|x2` refuses with that note.

Arms per fold and seed, all on try K10's drawn cells and stream (k9.STREAM; build seed = the harness's ctx.seed):
  s1x            the base (not a try): k10_s1_form.arm_built's s1x path exactly, build_submission.construct(ref, n,
                 rng, f=V.p, phi=phi_src, ratio=ratio_src), phi_src and ratio_src by k10_s1_form.source_side.
  s1x_half       K15a: the same call with phi_src replaced per gene by (1 + phi_src) / 2
                 (k14_form_arms.half_ctrl_factor, try K14b's form; no finalize), the depth step unchanged.
  s1x_half_perm  K15b, the control: the same call with the half factor's deviation from one permuted over genes within
                 deciles of the reference's mean expression (harness.t3_decile_permutation, basis
                 harness.t3_mean_expression of the reference, seed k9.PERM_SEED: try K9's decile permutation and seed).
                 Permuting phi_half itself is the same permutation of phi_half - 1, and keeps every value bit for bit.
No new random stream is drawn.

Reading, fixed (row K14's form, the adoption rule's constants: a paired mean above 0.3 with at least 8 of 10 seeds
ahead):
  Void unless one commit, a clean tree, every built fold present (r12 and r21), every seed and arm present, and s1x
  reproduces try K10's committed s1x scores on r12 and r21 to 1e-6 per seed (RESULTS_K10). A fold where severity is not
  scored in every seed of every arm is reported as not readable.
  Per fold: s1x_half - s1x and s1x_half - s1x_half_perm, each with its mean, seeds ahead and per-seed range.
  "Supported in the file's own form" if on x1 AND on r12 s1x_half beats s1x that way and beats s1x_half_perm that way,
  and on x2 and on r21 its paired mean against s1x is not below -0.3; anything else "not supported", with the failing
  folds named (a not-readable fold fails its checks). Where x1 and x2 could not be built (as here) the reading says
  "not readable on the cross-knockout folds" and gives the replicate folds' figures only.
  Reported, read by nothing: each arm's mean board score per fold, each difference's parts by metric, the share of genes
  the half factor moves by more than 1.5-fold against the full factor's.

  run   --task r12|r21|x1|x2 --out DIR    one task per process (refuses an existing task file; x1, x2: not on record)
  read  --runs DIR --out DIR [--k10 DIR]  apply the reading; write reading.json and reading.txt (refuses to overwrite)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k9_residual as k9  # noqa: E402  (the decile permutation's seed, the stream, the paired comparison)
import k10_s1_form as k10  # noqa: E402  (the replicate folds, the source side, the s1x build)
import k14_form_arms as k14  # noqa: E402  (try K14b's half factor)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

PERM_SEED = k9.PERM_SEED
STREAM = k10.STREAM
BASE, HALF, PERM = "s1x", "s1x_half", "s1x_half_perm"
ARMS = (BASE, HALF, PERM)
MIN_GAIN, MIN_AHEAD, FLOOR_MEAN = k9.MIN_GAIN, k9.MIN_AHEAD, k9.X2_FLOOR
REPLICATE = ("r12", "r21")
CROSS = ("x1", "x2")
DECIDING = ("x1", "r12")  # must beat s1x and the control
GUARDED = ("x2", "r21")  # paired mean against s1x not below FLOOR_MEAN
TASKS = {n: k10.TASKS[n] for n in REPLICATE}  # the folds built; x1 and x2 are not on record
KEYS = frozenset(t.key for t in TASKS.values())
SOURCE = k10.SOURCE
TASK_SOURCE = "experiments/t3_decomp/k15_file_form_half.py"
RESULTS_K10 = HERE / "results" / "k10_s1_form" / "tasks"
REPRO_TOL = 1e-6
FOLD_RATIO = 1.5  # the reported share's threshold, fixed in the brief
NOT_ON_RECORD = (
    "design not on record: the cross-knockout source's wild type is another stage than the fold's reference, so try "
    "K10's per-label depth ratio is in the source WT's label vocabulary, not the reference's; the record fixes no rule "
    "for carrying it across (which labels correspond, the ratio of a reference label without a counterpart, and "
    "whether the ratio also carries the depth difference between the two wild-type stages)")
HALF_RULE = "phi_half = (1 + phi_src) / 2 per gene (k14_form_arms.half_ctrl_factor, try K14b's form; no finalize)"
PERM_RULE = ("phi_half's deviation from one permuted over genes within deciles of the reference's t3_mean_expression "
             "(harness.t3_decile_permutation, seed k9.PERM_SEED: try K9's decile permutation and seed)")
_PREP: dict = {}


# ------------------------------------------------------------------------------------------------ the factors
def half_factor(phi_src) -> np.ndarray:
    """K15a's factor: try K14b's (1 + phi) / 2 per gene."""
    return k14.half_ctrl_factor(phi_src)


def permuted_half(phi_half, basis, seed: int = PERM_SEED) -> np.ndarray:
    """K15b's factor: phi_half permuted over genes within deciles of `basis` (try K9's permutation and seed); the same
    permutation of phi_half - 1, with every value kept bit for bit."""
    return harness.t3_decile_permutation(phi_half, basis, seed)


def share_beyond(phi, ratio: float = FOLD_RATIO) -> float:
    """The share of genes a factor moves by more than `ratio`-fold, either way."""
    a = np.abs(np.log(np.asarray(phi, dtype=np.float64)))
    return float((a > np.log(ratio)).mean()) if a.size else 0.0


# ------------------------------------------------------------------------------------------------ the arms
def arm_built(bs, ref, n: int, seed: int, f, arm: str, phi_src, ratio_src, phi_perm=None):
    """One arm's cells before assembly: try K10's s1x call (k10_s1_form.arm_built's s1x path), with the factor phi_src
    (s1x), its half (s1x_half) or the permuted half (s1x_half_perm, `phi_perm`); nothing else changes."""
    if arm == BASE:
        phi = np.asarray(phi_src, dtype=np.float64)
    elif arm == HALF:
        phi = half_factor(phi_src)
    elif arm == PERM:
        if phi_perm is None:
            raise ValueError("the control needs its permuted half factor")
        phi = np.asarray(phi_perm, dtype=np.float64)
    else:
        raise ValueError(f"unknown arm {arm!r}")
    return k10.arm_built(bs, ref, n, seed, f, k10.BASE, phi, ratio_src, None, None, None)


def _prepared(ctx) -> dict:
    """The source replicate's factor and depth ratio (k10_s1_form.source_side), the half and its permutation, and the
    reference in t3_decompose's form, once per loaded fold."""
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
        phi_src, ratio_src, labels = k10.source_side(wt, [ctx.inputs[r] for r in SOURCE["ko_roles"]], rules,
                                                     SOURCE["source_preset"])
        ref = harness.t3_reference(pool, SOURCE["reference_preset"], rules)
        if list(ref.V.labels) != labels:
            raise ValueError(f"{ctx.task.key}: the depth ratio's vocabulary is not the reference's")
        phi_half = half_factor(phi_src)
        _PREP.clear()
        hit = _PREP["fold"] = {
            "key": ctx.task.key, "inputs": objs, "phi_src": phi_src, "ratio_src": ratio_src, "ref": ref,
            "phi_perm": permuted_half(phi_half, harness.t3_mean_expression(pool)),
            "record": {"source": ctx.task.role_stems[SOURCE["ko_roles"][0]],
                       "reference": ctx.task.role_stems[ref_role],
                       "labels_with_a_depth_ratio": int((ratio_src != 1.0).sum()),
                       "share_beyond_1.5_fold": {"full": share_beyond(phi_src), "half": share_beyond(phi_half)}}}
    return hit


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = _prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    built = arm_built(bs, ref, harness.count_for(ctx.board, pool.n_obs), ctx.seed, ref.V.p, arm, prep["phi_src"],
                      prep["ratio_src"], prep["phi_perm"])
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    rules = {BASE: None, HALF: HALF_RULE, PERM: f"{HALF_RULE}; {PERM_RULE}"}

    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params=k10.PARAMS | {"k15_form": rules[arm],
                                                      "perm_seed": PERM_SEED if arm == PERM else None},
                                 requires=lambda t: t.key in KEYS)

    return [harness.CANDIDATES[harness.FLOOR], *(one(a) for a in ARMS)]


def cmd_run(a) -> None:
    if a.task in CROSS:
        sys.exit(f"{a.task}: not built, {NOT_ON_RECORD}")
    t = TASKS[a.task]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{a.task}.json"
    if path.exists():
        sys.exit(f"{path} exists")
    res = harness.run_task(t, candidates(), seeds=SEEDS, data_dir=Path(a.data))
    res["task_source"] = TASK_SOURCE
    res["factor_record"] = _PREP["fold"]["record"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def paired(rows: list[dict], x: str, y: str) -> dict:
    """k9._paired (mean, sd, seeds ahead, n, passes at the adoption rule's constants) and the per-seed range."""
    d = np.array([r["scores"][x] - r["scores"][y] for r in rows])
    return k9._paired(rows, x, y) | {"min": float(d.min()), "max": float(d.max())}


def k10_reproduction(rows_by_fold: dict[str, list[dict]], k10_dir: Path = RESULTS_K10) -> list[str]:
    """Reasons to void: s1x's per-seed scores against try K10's committed s1x scores (k10_dir/<fold>.json), to
    REPRO_TOL."""
    if not Path(k10_dir).is_dir():
        return [f"try K10's committed scores ({k10_dir}) are absent from this checkout"]
    why = []
    for name, rows in rows_by_fold.items():
        p = Path(k10_dir) / f"{name}.json"
        if not p.exists():
            why.append(f"{name}: try K10's committed result is missing from {k10_dir}")
            continue
        theirs = {r["seed"]: r["scores"][k10.BASE] for r in json.loads(p.read_text())["rows"]}
        ours = {r["seed"]: r["scores"][BASE] for r in rows}
        if set(theirs) != set(ours) or any(abs(ours[s] - theirs[s]) > REPRO_TOL for s in ours):
            why.append(f"{name}: s1x does not reproduce try K10's committed s1x scores")
    return why


def fold_figures(rows: list[dict], counted: bool) -> dict:
    """One fold's figures: readability (severity scored in every seed of every arm), the paired differences, each arm's
    mean board score, each difference's parts by metric and the severity states."""
    present = (harness.FLOOR, *ARMS) if all(harness.FLOOR in q["scores"] for q in rows) else ARMS
    states = {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                  for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in present}
    metrics = list(rows[0]["skills"][HALF])
    return {
        "counted_fold": counted,
        "readable": all(states[x] == {"scored": len(rows)} for x in ARMS),
        "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in present},
        "paired": {f"{x} - {y}": paired(rows, x, y) for x, y in ((HALF, BASE), (HALF, PERM), (PERM, BASE))},
        "severity_states": states,
        "skill_diff": {f"{x} - {y}": {m: float(np.mean([q["skills"][x][m] - q["skills"][y][m] for q in rows]))
                                      for m in metrics} for x, y in ((HALF, BASE), (HALF, PERM))}}


def decide(per: dict) -> tuple[bool, dict, list[str]]:
    """The fixed reading over the four folds: (supported, checks, the failing folds). A not-readable fold fails."""
    checks: dict[str, bool] = {}
    fails: list[str] = []
    for n in DECIDING:
        for y in (BASE, PERM):
            ok = per[n]["readable"] and per[n]["paired"][f"{HALF} - {y}"]["passes"]
            checks[f"{n}: {HALF} - {y} passes"] = bool(ok)
            if not ok and n not in fails:
                fails.append(n)
    for n in GUARDED:
        ok = per[n]["readable"] and per[n]["paired"][f"{HALF} - {BASE}"]["mean"] >= FLOOR_MEAN
        checks[f"{n}: {HALF} - {BASE} not below {FLOOR_MEAN}"] = bool(ok)
        if not ok:
            fails.append(n)
    return not fails, checks, fails


def verdict(per: dict) -> dict:
    """The verdict on the folds read: not readable on the cross-knockout folds when x1 or x2 is absent, else decide."""
    missing = [n for n in CROSS if n not in per]
    if missing:
        return {"supported": None, "checks": {}, "failing_folds": [],
                "verdict": "NOT READABLE ON THE CROSS-KNOCKOUT FOLDS (x1 and x2 not built: "
                           f"{NOT_ON_RECORD}); the replicate folds' figures only"}
    ok, checks, fails = decide(per)
    return {"supported": bool(ok), "checks": checks, "failing_folds": fails,
            "verdict": ("SUPPORTED IN THE FILE'S OWN FORM" if ok
                        else "NOT SUPPORTED (failing folds: " + ", ".join(fails) + ")")}


def reading(runs: Path, k10_dir: Path = RESULTS_K10) -> dict:
    why, per, commits, all_rows = [], {}, set(), {}
    for name, t in TASKS.items():
        p = Path(runs) / f"{name}.json"
        if not p.exists():
            why.append(f"missing task {name}")
            continue
        r = json.loads(p.read_text())
        rows = sorted(r["rows"], key=lambda x: x["seed"])
        if r["task"]["key"] != t.key or r.get("task_source") != TASK_SOURCE or [x["seed"] for x in rows] != list(SEEDS):
            why.append(f"{name}: another task, another script or other seeds")
            continue
        if not r["provenance"]["tree_clean"]:
            why.append(f"{name}: a dirty tree")
        missing = sorted(set(ARMS) - set.intersection(*(set(q["scores"]) for q in rows)))
        if missing:
            why.append(f"{name}: missing arms {missing}")
            continue
        commits.add(r["provenance"]["commit"])
        all_rows[name] = rows
        per[name] = fold_figures(rows, counted=False) | {"record": r.get("factor_record") or {}}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if not why:
        why += k10_reproduction(all_rows, k10_dir)
    if why:
        return {"void": True, "why": why}
    return {"void": False, "why": [], "commit": sorted(commits)[0], "k10_reproduction": "reproduced",
            "tasks": per, "not_readable_folds": [n for n in per if not per[n]["readable"]],
            "cross_knockout": NOT_ON_RECORD, **verdict(per), "label": LABEL, "limits": LIMITS}


LABEL = ("A labelled reading for a ruling, NOT the adoption rule; it changes no file: the first-choice file changes "
         "only by the captain's own ruling.")
LIMITS = ("Limits: the base is the frozen stand-in without the file's composition term; replicate folds are two "
          "samples on one slide; no released fold has the test's structure (a factor from one knockout applied to "
          "another at the same stage), and which fold the test resembles is not known; the folds' targets differ in "
          "depth against their references, which moved other reads on these folds; one factor, ten seeds of draw "
          "noise.")


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K15a-K15b (commit {res['commit']}): {res['verdict']}", res["label"], res["limits"],
             f"s1x against try K10's committed scores: {res['k10_reproduction']}"]
    if res["checks"]:
        lines.append("Checks: " + "; ".join(f"{k} {'yes' if v else 'no'}" for k, v in res["checks"].items()))
    for n, v in res["tasks"].items():
        kind = "counted cross-knockout fold" if v["counted_fold"] else "replicate fold, never counted"
        lines.append(f"{n} ({kind}){'' if v['readable'] else ', NOT READABLE (severity not scored throughout)'}")
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (ahead {d['ahead']} of {d['n']}, per-seed range "
                         f"{d['min']:+.2f} to {d['max']:+.2f}, sd {d['sd']:.2f})" + ("  passes" if d["passes"] else ""))
    lines.append("Reported, read by nothing:")
    for n, v in res["tasks"].items():
        lines.append(f"{n} mean board scores: " + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for k, sk in v["skill_diff"].items():
            lines.append(f"    skills {k}: " + ", ".join(f"{m} {x:+.3f}" for m, x in sk.items()))
        for a, st in v["severity_states"].items():
            lines.append(f"    {a:<14} severity " + ", ".join(f"{s} {c}" for s, c in st.items()))
        sh = v["record"].get("share_beyond_1.5_fold") or {}
        if sh:
            lines.append(f"    share of genes moved more than 1.5-fold: half {sh['half']:.3f}, full {sh['full']:.3f}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    for p in (out / "reading.json", out / "reading.txt"):
        if p.exists():
            sys.exit(f"{p} exists")
    res = reading(Path(a.runs), Path(a.k10))
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", choices=sorted(REPLICATE + CROSS), required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--out", required=True)
    d.add_argument("--k10", default=str(RESULTS_K10))
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
