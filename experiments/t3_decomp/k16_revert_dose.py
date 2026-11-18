#!/usr/bin/env python
"""Try K16 (experiments/harness/TRIES.md, registered before the run): the exponent of T3's wild-type-only factor, read
on the practice folds. A labelled reading for a ruling, not the adoption rule: g1, g2 and m are diagnostic folds that
tries K7 and K8 already read. Whatever it shows, the try changes no entry; at most a priced question goes to the
captain, who alone can change T3's frozen second-upload rule.

Why. T3's second upload (the refused branch's s1 built with --wt-revert; docs/S2_RULES.md) moves each drawn cell by
its lineage group's wild-type-only factor at exponent 0.5 (k7_wt_revert.EXPONENT, the arm of tries K7b and K8b). The
exponent was fixed in row K7 without a read of any other value. The scorer's severity is log beta, beta = dt.dp / dt.dt
over the truth's DE genes, gated to the worst value when beta <= 0 or r2 <= 0.01: r2 does not depend on the factor's
size and beta is proportional to it, de_score and de_direction are rank-based, and mmd_u and the variogram pay for
every move away from the truth. So the exponent trades response magnitude against mutant state.

Arms, all from one stream per task (the same drawn cells in every arm; each arm's per-group step starts from the same
generator), every arm the per-group factor of tries K7b and K8b (k7_wt_revert._prepared, unchanged):
  grp_e050       exponent 0.5: the base, not a try (on g1 and g2 it is K8b's arm, on m K7b's, value for value)
  grp_e075       exponent 0.75
  grp_e100       exponent 1.0
  grp_e150       exponent 1.5
  perm_e050      the control: each group's factor permuted across genes within deciles of the reference's mean
                 expression (K7c's permutation: vec.harness.t3_decile_permutation, seed k7_wt_revert.PERM_SEED; one
                 gene permutation, applied to every group's factor), exponent 0.5
  perm_e100      the same control at exponent 1.0
Forms: on g1 and g2 the upload's form (k8_s2_form: the central composition of the refused branch's s1, stream
"t3_s2_form"); on m K7b's form (the reference's cells at their own composition, stream "t3_wt_revert"), since the
central scenario is defined on WT E8.75's vocabulary.

The reading (fixed in the row before any run; `decide`):
  CARRIES      an exponent e in {0.75, 1.0, 1.5} carries if on both g1 and g2 it beats grp_e050 by a paired mean above
               0.3 with at least 8 of 10 seeds ahead and severity is scored (not gated) in at least 9 of 10 seeds, and
               on m its paired mean over grp_e050 is not below -0.3.
  SPECIFIC     at exponent 1.0, on both g1 and g2, (grp_e100 - grp_e050) exceeds (perm_e100 - perm_e050) by more than
               0.3 (paired means).
  STEP         every step up to e from the exponent tried before it beats that exponent by a paired mean above 0.3 on
               both g1 and g2: 0.75 (its step from 0.5 is part of carrying); 1.0 over 0.75, waived when 0.75 does not
               carry; 1.5 over 1.0.
  OUTCOME      the largest e that carries and whose every step up to it holds, and only if SPECIFIC holds: a priced
               question to the captain. Anything else: no request; the second upload stays at 0.5.

  run   --task g1|g2|m --out DIR    one task per process; start the three together (one commit, a clean tree)
  read  --runs DIR --out DIR        apply the reading; write reading.json and reading.txt
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
import k7_wt_revert as k7  # noqa: E402  (the factors, the tasks, the permutation seed, the paired comparison)
import k8_s2_form as k8  # noqa: E402  (the upload's form: the central composition and its stream)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

BASE_EXPONENT = k7.EXPONENT  # 0.5, the second upload's
EXPONENTS = (0.5, 0.75, 1.0, 1.5)
PERM_EXPONENTS = (0.5, 1.0)
MIN_GAIN, MIN_AHEAD = k7.MIN_GAIN, k7.MIN_AHEAD  # 0.3 and 8 of 10
MIN_SCORED = 9
M_FLOOR = -0.3


def arm(e: float, perm: bool = False) -> str:
    return f"{'perm' if perm else 'grp'}_e{round(e * 100):03d}"


BASE = arm(BASE_EXPONENT)
TRIED = tuple(e for e in EXPONENTS if e != BASE_EXPONENT)
ARM_SPEC = {**{arm(e): (e, False) for e in EXPONENTS}, **{arm(e, True): (e, True) for e in PERM_EXPONENTS}}
ARMS = tuple(ARM_SPEC)
TASKS = dict(k7.TASKS)  # g1, g2 and m
S2_TASKS = ("g1", "g2")  # read in the upload's form (k8); m in K7b's
S2_KEYS = frozenset(TASKS[n].key for n in S2_TASKS)
KEYS = frozenset(t.key for t in TASKS.values())
_PERM: dict = {}


def permuted_factors(prep: dict, pool) -> tuple[dict, np.ndarray]:
    """K7c's permutation applied to each group's factor and to the pooled fallback: one gene permutation (deciles of
    the reference's mean expression, seed k7.PERM_SEED), the same for every group. Made once per prepared fold."""
    if _PERM.get("prep") is not prep:
        basis = harness.t3_mean_expression(pool)
        _PERM.clear()
        _PERM.update(prep=prep, per={g: harness.t3_decile_permutation(v, basis, k7.PERM_SEED)
                                     for g, v in prep["per"].items()},
                     phi=harness.t3_decile_permutation(prep["phi"], basis, k7.PERM_SEED))
    return _PERM["per"], _PERM["phi"]


def move_by_group(bs, ref, built, group: np.ndarray, per: dict, fallback: np.ndarray, exponent: float, rng):
    """Each drawn cell's counts moved by its lineage group's factor raised to `exponent` (a group without its own
    factor takes `fallback`), through t3_decompose.apply_gene_factor with the pool's smoothed profile, the groups in
    sorted order on one generator: the step of tries K7b and K8b with the exponent as a parameter."""
    if not np.all(built.source == "wt"):
        raise RuntimeError("the grouped arm expects every drawn cell to come from the reference pool")
    C, P = ref.pool.C[built.row].copy(), ref.pool.P[built.row]
    grp = group[built.row]
    for g in sorted(set(grp)):
        rows = np.flatnonzero(grp == g)
        C[rows] = bs.t3d.apply_gene_factor(C[rows], per.get(g, fallback) ** exponent, P[rows], rng)
    return bs.Built(X=bs.t3d.lognorm(C), lab=built.lab, source=built.source, row=built.row)


def _build(ctx, name: str):
    bs = harness._t3_head()
    prep = k7._prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    exponent, perm = ARM_SPEC[name]
    if ctx.task.key in S2_KEYS:
        f, info = bs.composition(ref.V, k8.spec(), {})
        if info["states_missing"]:
            raise RuntimeError(f"states missing from the reference's vocabulary: {info['states_missing']}")
        stream = k8.STREAM
    else:
        f, stream = ref.V.p, k7.STREAM
    n = harness.count_for(ctx.board, pool.n_obs)
    built = bs.construct(ref, n, bs.t3d.stable_rng(ctx.seed, stream), f=f)
    per, fallback = permuted_factors(prep, pool) if perm else (prep["per"], prep["phi"])
    built = move_by_group(bs, ref, built, prep["group"], per, fallback, exponent,
                          bs.t3d.stable_rng(ctx.seed, f"{stream}:grouped"))
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{stream}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    def one(name):
        e, perm = ARM_SPEC[name]
        return harness.Candidate(name, frozenset({"T3"}), (), lambda ctx, name=name: _build(ctx, name),
                                 params={"factor": k7.WT_REVERT_RULE + " (per lineage group)", "exponent": e,
                                         "perm_seed": k7.PERM_SEED if perm else None,
                                         "form": {"g1/g2": {"composition": k8.spec()["composition"],
                                                            "stream": k8.STREAM},
                                                  "m": {"composition": "the reference's own", "stream": k7.STREAM}}},
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
    res["task_source"] = "experiments/t3_decomp/k16_revert_dose.py"
    res["factor_record"] = k7._PREP["fold"]["record"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def _median(v: list) -> float | None:
    v = [float(x) for x in v if x is not None]
    return float(np.median(v)) if v else None


def summarise(rows: list[dict]) -> dict:
    """Per arm: the board score's mean and sd, each metric's mean skill, the severity states and the number of seeds
    scored, the median beta and r2; the paired differences named by the reading."""
    arms = (harness.FLOOR, *ARMS)
    out = {"arms": {}, "paired": {}}
    for x in arms:
        s = np.array([q["scores"][x] for q in rows], dtype=float)
        sev = [q["t3"][x]["severity"] for q in rows]
        out["arms"][x] = {
            "mean": float(s.mean()), "sd": float(s.std(ddof=1)),
            "skills": {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]},
            "severity_states": {st: sum(v["state"] == st for v in sev) for st in sorted({v["state"] for v in sev})},
            "severity_scored": int(sum(v["state"] == "scored" for v in sev)),
            "median_beta": _median([v.get("beta") for v in sev]), "median_r2": _median([v.get("r2") for v in sev])}
    pairs = [(x, BASE) for x in ARMS if x != BASE] + [(arm(1.0, True), arm(0.5, True)), (arm(1.0), arm(0.75)),
                                                       (arm(1.5), arm(1.0)), (BASE, harness.FLOOR)]
    out["paired"] = {f"{x} - {y}": k7._paired(rows, x, y) for x, y in pairs}
    return out


def decide(per: dict) -> dict:
    """The reading fixed in row K16, on summarise()'s output per task."""
    def pm(n, x, y):
        return per[n]["paired"][f"{x} - {y}"]

    def beats(x, y):
        return all(pm(n, x, y)["mean"] > MIN_GAIN for n in S2_TASKS)

    carries = {}
    for e in TRIED:
        a = arm(e)
        carries[a] = bool(all(pm(n, a, BASE)["mean"] > MIN_GAIN and pm(n, a, BASE)["ahead"] >= MIN_AHEAD
                              and per[n]["arms"][a]["severity_scored"] >= MIN_SCORED for n in S2_TASKS)
                          and pm("m", a, BASE)["mean"] >= M_FLOOR)
    gain = {n: pm(n, arm(1.0), BASE)["mean"] - pm(n, arm(1.0, True), arm(0.5, True))["mean"] for n in S2_TASKS}
    specific = bool(all(v > MIN_GAIN for v in gain.values()))
    step = {arm(0.75): True,
            arm(1.0): bool(beats(arm(1.0), arm(0.75)) or not carries[arm(0.75)]),
            arm(1.5): bool(beats(arm(1.5), arm(1.0)))}
    named = None
    for i, e in enumerate(TRIED):
        if carries[arm(e)] and all(step[arm(x)] for x in TRIED[:i + 1]):
            named = e
    request = named is not None and specific
    return {"carries": carries, "specific_gain_over_control": gain, "factor_specific": specific, "step": step,
            "named_exponent": named, "request": bool(request),
            "outcome": (f"exponent {named:g}: a priced question to the captain (T3's frozen second-upload rule; a "
                        "builder option and a rehearsal build); the try changes no entry" if request else
                        "no request; the second upload stays at exponent 0.5")}


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
        per[name] = summarise(rows) | {"factor_record": r.get("factor_record")}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if why:
        return {"void": True, "why": why}
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per} | decide(per)


def _f(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Try K16 (commit {res['commit']}): {res['outcome']}.",
             "carries: " + ", ".join(f"{a} {'yes' if v else 'no'}" for a, v in res["carries"].items())
             + "; steps: " + ", ".join(f"{a} {'yes' if v else 'no'}" for a, v in res["step"].items())
             + "; factor-specific at 1.0: " + ("yes" if res["factor_specific"] else "no") + " ("
             + ", ".join(f"{n} {v:+.2f}" for n, v in res["specific_gain_over_control"].items()) + ")"]
    for n, t in TASKS.items():
        v = res["tasks"][n]
        lines.append(f"{n} {t.key} ({'the upload form' if n in S2_TASKS else 'K7b form'}):")
        for a, s in v["arms"].items():
            lines.append(f"    {a:<10} {s['mean']:.2f} (sd {s['sd']:.2f}); severity scored {s['severity_scored']} of "
                         f"{len(SEEDS)}, median beta {_f(s['median_beta'])}, r2 {_f(s['median_r2'])}; skills "
                         + ", ".join(f"{m} {x:.3f}" for m, x in s["skills"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    res = reading(Path(a.runs))
    out = Path(a.out)
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
