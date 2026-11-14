#!/usr/bin/env python
"""Tries K13 (the T3 arm raw_tail; specification fixed by the main session before any code): the tail decoder on
K9's knockout-derived factor, in the form of the heart's try H17b (campaign/heart-b, vec.harness.tail_infeasible,
tail_filter, B0_TREND_TAIL_RULE). A labelled reading for a ruling, NOT the adoption rule.

Why. The ledger's "Pruned before running" holds "T3, heart and embryo: tail decoders"; T3's was never read. A truth
down-gene needs a reference level of at least 0.25 (the scorer's threshold on the reference's mean log1p
expression), so a gene below that level with a predicted decrease is a sure miss in the predicted down-list. The T3
arms change expression by a per-gene factor, as the heart's arm does, so the heart's form applies as it stands: no
kernel, no bound, no solve.

Arms (K9's stream, the same drawn cells in every arm; all through k9.arm_built's raw path):
  raw             K9's raw arm exactly.
  raw_tail        eb_factor.finalize(y'), y the factor's log ratio as raw_c3 takes it (k11.section_estimates) and y'
                  y with every negative entry on a down-infeasible gene set to 0; a gene is down-infeasible if the
                  fold's reference level (harness.t3_mean_expression of the reference: the mean over its cells of
                  log1p CP10k, the level k12 uses) is below TAIL_PB_MIN. With no such gene, or none with a negative
                  entry, raw's factor bit for bit.
  raw_tail_perm   the control: the same filter with the down-infeasible mask permuted across genes within deciles of
                  the reference's mean expression (harness.t3_decile_permutation, k9.PERM_SEED).

Each run prints one check line per task (numbers of genes only) and writes the same numbers into the task file.

Reading, fixed (K9's form, the adoption rule's constants, as k12): raw_tail is SUPPORTED if on x1 AND on r12 it beats
raw by a paired mean above 0.3 with at least 8 of 10 seeds ahead and beats raw_tail_perm the same way, and on x2 its
paired mean against raw is not below -0.3; anything else CLOSED. A fold where the filter changes no gene is IDLE (the
arm is raw there); idle on x1 or on r12 means CLOSED with the word "idle" in the verdict. Void if K9's committed raw
scores (RESULTS_K9) are absent or raw does not reproduce them to 1e-6 (k11.k9_reproduction), on a missing seed, arm or
check record, or a dirty tree.

  run   --task x1|x2|r12 --out DIR    one task per process (refuses an existing task file)
  read  --runs DIR --out DIR          apply the reading; write reading.json and reading.txt (refuses to overwrite)
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
import k9_residual as k9  # noqa: E402  (K9's folds, base arm, control seed and reading form)
import k11_factor_arms as k11  # noqa: E402  (raw_c3's log ratio y and the reproduction check)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

# a truth down-gene needs pb(reference) >= 0.25 (the scorer's threshold; the heart's vec.harness.TAIL_PB_MIN on
# campaign/heart-b, experiments/research_2026_11/diagnostics.py, N3/N5); this branch does not hold it, so defined here
TAIL_PB_MIN = 0.25
PERM_SEED = k9.PERM_SEED
STREAM = k9.STREAM
BASE = k9.BASE  # "raw"
TAIL, TAIL_PERM = "raw_tail", "raw_tail_perm"
ARMS = (BASE, TAIL, TAIL_PERM)
CONTROLS = {TAIL: TAIL_PERM}
MIN_GAIN, MIN_AHEAD, X2_FLOOR = k9.MIN_GAIN, k9.MIN_AHEAD, k9.X2_FLOOR
TASKS = k9.TASKS
KEYS = k9.KEYS
FOLD_OF = k9.FOLD_OF
COUNTED = k9.COUNTED
DECIDING = ("x1", "r12")  # the folds the arm must win on; idle on either closes the line
RESULTS_K9 = k11.RESULTS_K9
TAIL_RULE = ("eb_factor.finalize(y'), y gene_factor's log ratio as raw_c3 takes it and y' y with every negative entry "
             "on a down-infeasible gene set to 0: a gene whose reference level (harness.t3_mean_expression, the mean "
             "over the reference's cells of log1p CP10k) is below TAIL_PB_MIN cannot be a truth down-gene (the "
             "heart's B0_TREND_TAIL_RULE form, campaign/heart-b); no gene changed: raw's factor bit for bit")
TAIL_PERM_RULE = ("as raw_tail with the down-infeasible mask permuted across genes within deciles of the reference's "
                  "t3_mean_expression (harness.t3_decile_permutation, k9.PERM_SEED)")
_PREP: dict = {}


# ------------------------------------------------------------------------------------------------ the filter
def tail_infeasible(level: np.ndarray, pb_min: float | None = None) -> np.ndarray:
    """The down-infeasible mask: the reference level (harness.t3_mean_expression) below TAIL_PB_MIN."""
    return np.asarray(level, dtype=np.float64) < (TAIL_PB_MIN if pb_min is None else pb_min)


def tail_filter(y: np.ndarray, infeasible: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(y with every negative entry on an infeasible gene set to 0, the mask of the genes it changed); every other
    entry is returned as it came."""
    out = np.array(y, dtype=np.float64)
    hit = np.asarray(infeasible, dtype=bool) & (out < 0)
    out[hit] = 0.0
    return out, hit


def tail_factor(y: np.ndarray, den: np.ndarray, clip: tuple[float, float], infeasible: np.ndarray,
                phi_raw: np.ndarray) -> np.ndarray:
    """raw_tail's factor: eb_factor.finalize of the filtered log ratio; with no gene changed, phi_raw itself (the arm
    is raw bit for bit)."""
    yf, hit = tail_filter(y, infeasible)
    if not hit.any():
        return np.array(phi_raw, dtype=np.float64)
    return harness._t3_eb().finalize(yf, den, clip)


def permuted_mask(infeasible: np.ndarray, level: np.ndarray, seed: int = PERM_SEED) -> np.ndarray:
    """The control's mask: the infeasible mask permuted across genes within deciles of `level`
    (harness.t3_decile_permutation), so each decile keeps its number of masked genes."""
    m = np.asarray(infeasible, dtype=bool).astype(np.float64)
    return harness.t3_decile_permutation(m, level, seed) > 0.5


def tail_check(y: np.ndarray, infeasible: np.ndarray, control: np.ndarray) -> dict[str, int]:
    """The check's numbers of genes: below the level, of those with a negative entry (the genes the filter changes),
    and the genes the control changes."""
    neg = np.asarray(y, dtype=np.float64) < 0
    infeasible, control = np.asarray(infeasible, dtype=bool), np.asarray(control, dtype=bool)
    return {"infeasible": int(infeasible.sum()), "negative": int((infeasible & neg).sum()),
            "control_changed": int((control & neg).sum())}


def check_line(task: str, chk: dict) -> str:
    return (f"raw_tail check: task={task} infeasible={chk['infeasible']} negative={chk['negative']} "
            f"control_changed={chk['control_changed']}")


# ------------------------------------------------------------------------------------------------ the arms
def _prepared(ctx) -> dict:
    """The fold's raw factor (K9's), the down-infeasible mask from its reference, raw_tail's factor and its control's,
    and the reference in t3_decompose's form, once per loaded fold."""
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
        level = harness.t3_mean_expression(pool)
        infeasible = tail_infeasible(level)
        control = permuted_mask(infeasible, level, PERM_SEED)
        y, _, den, _ = k11.section_estimates(frame.C, frame.lab, frame.block, frame.mu_wt, frame.weights, frame.mask,
                                             frame.eps)
        clip = tuple(float(c) for c in rules["components"]["factor"]["clip"])
        chk = tail_check(y, infeasible, control)
        _PREP.clear()
        hit = _PREP["fold"] = {
            "key": ctx.task.key, "inputs": objs, "phi_raw": phi_raw, "infeasible": infeasible, "control": control,
            "phi_tail": tail_factor(y, den, clip, infeasible, phi_raw),
            "phi_tail_perm": tail_factor(y, den, clip, control, phi_raw),
            "clip": clip, "ref": harness.t3_reference(pool, fold["reference_preset"], rules),
            "record": {"fold": FOLD_OF[ctx.task.key], "reference": ctx.task.role_stems[ref_role],
                       "pb_min": TAIL_PB_MIN, "tail_check": chk,
                       "both_changed": int((infeasible & control & (np.asarray(y) < 0)).sum())}}
    return hit


def arm_built(bs, ref, n: int, seed: int, f, arm: str, prep: dict):
    """One arm's cells through k9.arm_built's raw path, each with its own factor: phi_raw, phi_tail or phi_tail_perm."""
    phi = {BASE: "phi_raw", TAIL: "phi_tail", TAIL_PERM: "phi_tail_perm"}[arm]
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
                                     "tail": TAIL_RULE if arm in (TAIL, TAIL_PERM) else None,
                                     "pb_min": TAIL_PB_MIN if arm in (TAIL, TAIL_PERM) else None,
                                     "perm_seed": PERM_SEED if arm == TAIL_PERM else None} | (
                                     {"tail_perm": TAIL_PERM_RULE} if arm == TAIL_PERM else {}),
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
    rec = _PREP["fold"]["record"]
    res["task_source"] = "experiments/t3_decomp/k13_tail_arm.py"
    res["factor_record"] = rec
    res["tail_check"] = rec["tail_check"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(check_line(a.task, rec["tail_check"]), flush=True)
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def decide(per: dict) -> tuple[bool, dict]:
    """raw_tail passes only when it beats raw and raw_tail_perm on x1 and on r12, and its paired mean against raw on
    x2 is not below X2_FLOOR."""
    checks = {f"{n}: {TAIL} - {y}": per[n]["paired"][f"{TAIL} - {y}"]["passes"]
              for n in DECIDING for y in (BASE, TAIL_PERM)}
    checks[f"x2: {TAIL} - {BASE} not below {X2_FLOOR}"] = per["x2"]["paired"][f"{TAIL} - {BASE}"]["mean"] >= X2_FLOOR
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
        missing = sorted(set(ARMS) - set.intersection(*(set(q["scores"]) for q in rows)))
        if missing:
            why.append(f"{name}: missing arms {missing}")
            continue
        chk = r.get("tail_check")
        if not isinstance(chk, dict) or not {"infeasible", "negative", "control_changed"} <= set(chk):
            why.append(f"{name}: the raw_tail check record is missing")
            continue
        commits.add(r["provenance"]["commit"])
        all_rows[name] = rows
        arms = (harness.FLOOR, *ARMS)
        metrics = list(rows[0]["skills"][BASE])
        per[name] = {
            "counted_fold": name in COUNTED,
            "idle": int(chk["negative"]) == 0,
            "tail_check": chk,
            "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in arms},
            "paired": {f"{x} - {y}": k9._paired(rows, x, y)
                       for x, y in ((BASE, harness.FLOOR), (TAIL, BASE), (TAIL, TAIL_PERM), (TAIL_PERM, BASE))},
            "severity_states": {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                                    for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in arms},
            "skill_diff_minus_raw": {a: {m: float(np.mean([q["skills"][a][m] - q["skills"][BASE][m] for q in rows]))
                                         for m in metrics} for a in (TAIL, TAIL_PERM)},
            "record": r.get("factor_record") or {}}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if not Path(k9_dir).is_dir():
        why.append(f"K9's committed raw scores ({k9_dir}) are absent from this checkout")
    status, repro_why = k11.k9_reproduction(all_rows, k9_dir) if not why else ("not checked", [])
    why += repro_why
    if why:
        return {"void": True, "why": why}
    ok, checks = decide(per)
    idle = [n for n in DECIDING if per[n]["idle"]]
    supported = bool(ok) and not idle
    if idle:
        verdict = (f"CLOSED: idle on {' and '.join(idle)} (the filter changes no gene there, so the arm is raw); "
                   f"no request")
    else:
        verdict = f"SUPPORTED: {NEXT}" if supported else "CLOSED: no request"
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "k9_reproduction": status,
            "checks": checks, "idle_folds": [n for n in TASKS if per[n]["idle"]], "supported": supported,
            "verdict": verdict, "label": LABEL}


NEXT = ("read next on the first-choice file's own form (as K10 did) before any request; no request goes on this "
        "reading alone")
LABEL = ("A labelled reading for a ruling, NOT the adoption rule: only x1 and x2 are counted cross-knockout folds; r12 "
         "is a replicate fold that never counts; x2's reference is WT E9.5, where K7 found reverts flattered.")


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K13 (commit {res['commit']}): {TAIL} {res['verdict']}", res["label"],
             f"raw against K9's committed raw scores: {res['k9_reproduction']}",
             "Checks: " + "; ".join(f"{k} {'yes' if x else 'no'}" for k, x in res["checks"].items()),
             "Reported, deciding nothing:"]
    for n, t in TASKS.items():
        v = res["tasks"][n]
        kind = "counted cross-knockout fold" if v["counted_fold"] else "replicate fold, never counted"
        lines.append(f"{n} {FOLD_OF[t.key]} ({kind}){' IDLE' if v['idle'] else ''}: "
                     + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        lines.append("    " + check_line(n, v["tail_check"]))
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
