#!/usr/bin/env python
"""Tries K8a and K8b (experiments/harness/TRIES.md, registered before the run): the wild-type-only factor of tries K7
on T3's second upload as it is built today, on the two Gata4 replicates. Information for the captain, who alone can
reopen the frozen T3 rules: no request follows from any outcome, and K7's verdict (not supported by its reading)
stands whatever this shows.

Why. Tries K7 applied the factor to the reference's cells at their own composition. T3's second upload is not that
file: it is the refused branch's s1 of the frozen decision file, wild-type cells with the central composition (five
states depleted by 2^-0.7178, the Mab21l2 magnitude) and no factor. A composition change moves the pseudobulk on its
own, and the factor's gain on Gata4 was the scorer's severity gate cleared by a thin margin (r2 0.01 to 0.07), so the
gain may not carry to the upload's own form. Decision log 2026-11-13 (b) lists that as an open risk.

Arms, all drawn from one stream (the same cells in every arm), built by the T3 builder's own functions
(build_submission.composition, construct and assemble):
  s2_form             today's second upload: the central composition, no factor (the base; not a try)
  s2_form_revert      K8a: the same with K7a's pooled factor at exponent 0.5
  s2_form_revert_grp  K8b: the same with K7b's factor, applied per lineage group
The factors are K7's, unchanged (k7_wt_revert._prepared). Tasks: K7's g1 and g2 (Gata4 replicate 1 and 2 against WT
E8.75, the other wild-type stage files as the only inputs). The Mab21l2 task is left out: its reference has another
vocabulary, and the central scenario is defined on WT E8.75's.

Reading, fixed before the run: the gain "carries to the upload's form" for an arm if it beats s2_form by a paired mean
above 0.3 with at least 8 of 10 seeds ahead on both replicates. Reported with it: each arm against the floor, the
severity states and the skills per metric. Nothing is requested on any outcome.

  run   --task g1|g2 --out DIR    one task per process; start the two together (one commit, a clean tree)
  read  --runs DIR --out DIR      apply the reading; write reading.json and reading.txt
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
import k7_wt_revert as k7  # noqa: E402  (the factors and the tasks of tries K7, unchanged)

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

# experiments/t3_decomp/results/gata4/gata4_decisions.json, submissions.refused.s1 (a test pins both to that file)
DEPLETED_STATES = ("OFT/RV-CM", "aSHF", "pSHF", "aPHM", "pPHM")
LOG2_MAGNITUDE = 0.717813870986706
STREAM = "t3_s2_form"
BASE, POOLED, GROUPED = "s2_form", "s2_form_revert", "s2_form_revert_grp"
ARMS = (BASE, POOLED, GROUPED)
TASKS = {n: k7.TASKS[n] for n in ("g1", "g2")}
KEYS = frozenset(t.key for t in TASKS.values())


def spec() -> dict:
    """The refused branch's s1 as the frozen decision file states it: no factor, the central composition."""
    return {"factor": {"source": "none", "exponent": 0.0},
            "composition": {"kind": "central", "depleted_states": list(DEPLETED_STATES),
                            "log2_magnitude": LOG2_MAGNITUDE}}


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = k7._prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    f, info = bs.composition(ref.V, spec(), {})
    if info["states_missing"]:
        raise RuntimeError(f"states missing from the reference's vocabulary: {info['states_missing']}")
    n = harness.count_for(ctx.board, pool.n_obs)
    rng = bs.t3d.stable_rng(ctx.seed, STREAM)
    if arm == POOLED:
        built = bs.construct(ref, n, rng, f=f, phi=prep["phi"], w=k7.EXPONENT)
    else:
        built = bs.construct(ref, n, rng, f=f)
        if arm == GROUPED:
            if not np.all(built.source == "wt"):
                raise RuntimeError("the grouped arm expects every drawn cell to come from the reference pool")
            C, P = ref.pool.C[built.row].copy(), ref.pool.P[built.row]
            grp = prep["group"][built.row]
            rng_g = bs.t3d.stable_rng(ctx.seed, f"{STREAM}:grouped")
            for g in sorted(set(grp)):
                rows = np.flatnonzero(grp == g)
                C[rows] = bs.t3d.apply_gene_factor(C[rows], prep["per"].get(g, prep["phi"]) ** k7.EXPONENT, P[rows],
                                                   rng_g)
            built = bs.Built(X=bs.t3d.lognorm(C), lab=built.lab, source=built.source, row=built.row)
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params={"composition": spec()["composition"], "stream": STREAM,
                                         "factor": "none" if arm == BASE else k7.WT_REVERT_RULE},
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
    res["task_source"] = "experiments/t3_decomp/k8_s2_form.py"
    res["factor_record"] = k7._PREP["fold"]["record"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


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
        per[name] = {
            "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in arms},
            "skills": {x: {m: float(np.mean([q["skills"][x][m] for q in rows])) for m in rows[0]["skills"][x]}
                       for x in arms},
            "severity_states": {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                                    for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in arms},
            "paired": {f"{x} - {y}": k7._paired(rows, x, y) for x, y in (
                (POOLED, BASE), (GROUPED, BASE), (BASE, harness.FLOOR), (POOLED, harness.FLOOR),
                (GROUPED, harness.FLOOR))}}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if why:
        return {"void": True, "why": why}
    carries = {x: bool(all(per[n]["paired"][f"{x} - {BASE}"]["passes"] for n in TASKS)) for x in (POOLED, GROUPED)}
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "carries": carries}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K8a-K8b (commit {res['commit']}): the gain carries to the upload's form: "
             + ", ".join(f"{a} {'yes' if v else 'no'}" for a, v in res["carries"].items())
             + ". No request follows; K7's verdict stands."]
    for n, t in TASKS.items():
        v = res["tasks"][n]
        lines.append(f"{n} {t.key}: " + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})"
                         + ("  passes" if d["passes"] else ""))
        for a in (harness.FLOOR, *ARMS):
            lines.append(f"    {a:<19} severity " + ", ".join(f"{s} {c}" for s, c in v["severity_states"][a].items())
                         + "; skills " + ", ".join(f"{m} {x:.3f}" for m, x in v["skills"][a].items()))
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
