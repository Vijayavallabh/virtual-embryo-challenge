#!/usr/bin/env python
"""Tries K7a to K7c (experiments/harness/TRIES.md, registered before the run): the wild-type reference's own per-gene
deviation, estimated from the other released wild-type samples alone, as the factor of a T3 prediction that holds no
knockout data. A labelled diagnostic for a ruling: T3's uploads are fixed by frozen rules and change only by the
captain's own say-so, whatever this shows.

Why. The campaign of 2026-11-12 read that released MERFISH samples differ in a per-gene profile of their own, larger
than the change between neighbouring stages (docs/DATA_NOTES.md, "a series is not one run", "successive changes
reverse at short gaps"). A knockout is scored against a wild-type reference of the same stage, so what the scorer
calls its response is (the knockout sample's own deviation) minus (the reference sample's own deviation) plus the
genotype's effect. The second term is a property of the reference alone, and it is the same for every knockout scored
against that reference: for the test, WT E8.75. Today's second T3 upload carries no per-gene factor, so its severity
sits at the gate. The pre-registered script's component R tried this term once, from two or three samples at
exponent 1, and scored at the floor (experiments/t3_decomp/results/gata4: R 50.49 and 48.17). The committed rows also
show the mechanism from the other side: the Mab21l2 factor measured against WT E8.75 scores 59.40 and 62.15 on Gata4,
and the same knockout's factor against its own WT E9.5 scores 52.04 and 46.24.

The factor (WT_REVERT_RULE, every constant fixed here). Samples: the task's reference and the six other released
wild-type stage files that carry the full panel (E7.5 to E10.5). Per sample and lineage group (obs "group" of
vec.prep.load_stage, "Unknown" left out, at least MIN_CELLS cells with counts): the mean within-cell panel proportion
of each gene. Per group held by at least MIN_SAMPLES other samples: log(mean + EPS) fitted on the samples' stage days
by least squares and read at the reference's day (rules.yaml's series_factor form, with every sample instead of two or
three). The factor is t3_decompose.gene_factor of that trend against the reference's own group means over the groups
both hold, weights the reference's group sizes, eps 1e-5, clip 0.25 to 4, applied at exponent EXPONENT = 0.5 by the
T3 build itself (build_submission.construct: binomial thinning, Poisson top-up from the kNN-smoothed profile).

Arms, all drawn from one stream (the same wild-type cells in every arm):
  wt_draw          the reference's cells at its own composition, no factor (the base; not a try)
  wt_revert        K7a: the pooled factor above
  wt_revert_grp    K7b: each cell moved by its own group's factor (the same trend, one gene_factor per group; a group
                   without a trend takes the pooled factor)
  wt_revert_perm   K7c, the control: K7a's factor permuted across genes within deciles of the reference's mean
                   expression (vec.harness.t3_decile_permutation, seed PERM_SEED)

Tasks (diagnostic pseudo-tasks checked by vec.harness.check_task; configs/pseudo_tasks.yaml unchanged). No knockout
file is an input of any of them:
  g1   Gata4 replicate 1 against WT E8.75      g2   Gata4 replicate 2 against WT E8.75
  m    the Mab21l2 knockout against WT E9.5
The two Gata4 replicates share a slide and the reference, so they are one knockout run read twice; m is a second run
with another reference. On g1 and g2 the factor is the very vector a beta-catenin build would apply.

Reading, fixed before the run (the adoption rule's constants: a paired mean gain above 0.3 with at least 8 of 10
seeds ahead). K7a is supported if on each of g1, g2 and m it beats wt_draw, the floor and wt_revert_perm that way. Then
a priced request goes to the captain to add the factor to T3's second upload, which would still hold no number from
the Gata4 files. K7b replaces K7a in that request only if it beats K7a that way on all three tasks. Anything else:
closed, no request. Reported with it: each arm's severity states over the seeds and its skills per metric.

  run   --task g1|g2|m --out DIR    one task per process; start the three together (one commit, a clean tree)
  read  --runs DIR --out DIR        apply the reading; write reading.json and reading.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

from vec import harness
from vec.adoption import SEEDS

EPS = 1e-5
MIN_CELLS = 200
MIN_SAMPLES = 3
EXPONENT = 0.5
CLIP = (0.25, 4.0)
PERM_SEED = 20261113
STREAM = "t3_wt_revert"
FULL_PANEL = ("E7.5", "E8.0", "E8.25_late", "E8.5", "E8.75", "E9.5", "E10.5")
PRESET = {"E8.75": "gata4", "E9.5": "mab21l2"}  # rules.yaml preset whose reference (hence label mode) is that stage
BASE, MAIN, GROUPED, CONTROL = "wt_draw", "wt_revert", "wt_revert_grp", "wt_revert_perm"
ARMS = (BASE, MAIN, GROUPED, CONTROL)
MIN_GAIN, MIN_AHEAD = 0.3, 8
WT_REVERT_RULE = ("per lineage group with >= 200 cells: mean within-cell panel proportions; per group held by >= 3 "
                  "other wild-type samples: log(mean + 1e-5) fitted on stage day by least squares, read at the "
                  "reference's day; t3_decompose.gene_factor of that trend against the reference's group means "
                  "(weights: the reference's group sizes; eps 1e-5; clip 0.25-4); exponent 0.5")
_PREP: dict = {}


def _task(key: str, target: str, reference: str, validation: str | None = None, forbidden: tuple[str, ...] = ()):
    inputs = {"wt": f"T2/{reference}.h5ad"} | {f"s_{s}": f"T2/{s}.h5ad" for s in FULL_PANEL if s != reference}
    return harness.check_task(harness.PseudoTask(
        key=key, test_board="T3:bcatenin", task="T3", kind="ko", target=f"T3/{target}.h5ad",
        reference=f"T2/{reference}.h5ad", inputs=inputs, validation_board=validation,
        twin_board=None if validation else "T3:gata4", forbidden=forbidden, adoption=False,
        notes="Diagnostic K7: a knockout against its wild type, with the other wild-type stages as the only inputs."))


TASKS = {
    "g1": _task("T3w/WT->Gata4rep1", "E8.75_gata4_ko_1", "E8.75", forbidden=("T3/E8.75_gata4_ko_2.h5ad",)),
    "g2": _task("T3w/WT->Gata4rep2", "E8.75_gata4_ko_2", "E8.75", validation="T3:gata4",
                forbidden=("T3/E8.75_gata4_ko_1.h5ad",)),
    "m": _task("T3w/WT->Mab21l2", "E9.5_mab21l2_ko", "E9.5"),
}
KEYS = frozenset(t.key for t in TASKS.values())


def group_means(a) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    """Per lineage group with at least MIN_CELLS cells that hold counts: the mean within-cell proportion of each gene,
    and the number of those cells. "Unknown" is left out."""
    C = sparse.csr_matrix(a.layers["counts"], dtype=np.float64)
    tot = np.asarray(C.sum(1)).ravel()
    P = sparse.diags(1.0 / np.maximum(tot, 1.0)) @ C
    grp = np.asarray(a.obs["group"]).astype(str)
    means, sizes = {}, {}
    for g in sorted(set(grp) - {"Unknown"}):
        rows = np.flatnonzero((grp == g) & (tot > 0))
        if len(rows) >= MIN_CELLS:
            means[str(g)], sizes[str(g)] = np.asarray(P[rows].mean(0)).ravel(), int(len(rows))
    return means, sizes


def trend_at(means: list[dict[str, np.ndarray]], days: list[float], day: float) -> dict[str, np.ndarray]:
    """Per group held by at least MIN_SAMPLES of the samples: log(mean + EPS) fitted on the samples' days by least
    squares and read at `day`, returned as a proportion (never below 0). Samples of one day give their geometric
    mean."""
    out = {}
    for g in sorted(set().union(*means)):
        have = [(d, m[g]) for d, m in zip(days, means, strict=True) if g in m]
        if len(have) < MIN_SAMPLES:
            continue
        x = np.array([d for d, _ in have], dtype=np.float64)
        Y = np.log(np.stack([m for _, m in have]) + EPS)
        xc = x - x.mean()
        sxx = float(xc @ xc)
        slope = (xc[:, None] * (Y - Y.mean(0))).sum(0) / sxx if sxx > 0 else np.zeros(Y.shape[1])
        out[g] = np.clip(np.exp(Y.mean(0) + slope * (day - x.mean())) - EPS, 0.0, None)
    return out


def revert_factors(ref_means: dict, ref_sizes: dict, trend: dict) -> tuple[np.ndarray, dict[str, np.ndarray], dict]:
    """(the pooled factor, one factor per group, a record). The pooled factor is t3_decompose.gene_factor of the trend
    against the reference's group means over the groups both hold, weighted by the reference's group sizes; a group's
    own factor is the same function on that group alone."""
    t3d = harness._t3_head().t3d
    groups = [g for g in sorted(ref_means) if g in trend]
    if not groups:
        raise ValueError("no lineage group is held by the reference and by enough other samples")
    mu_o, mu_r = np.stack([trend[g] for g in groups]), np.stack([ref_means[g] for g in groups])
    w = np.array([ref_sizes[g] for g in groups], dtype=np.float64)
    phi, _ = t3d.gene_factor(mu_o, mu_r, w / w.sum(), np.ones(len(groups), bool), EPS, CLIP)
    one = np.ones(1)
    per = {g: t3d.gene_factor(mu_o[i:i + 1], mu_r[i:i + 1], one, one.astype(bool), EPS, CLIP)[0]
           for i, g in enumerate(groups)}
    rec = {"groups": groups, "sd_log_phi": float(np.log(phi).std()),
           "genes_at_clip": int(((phi <= CLIP[0]) | (phi >= CLIP[1])).sum()),
           "genes_beyond_twofold": int((np.abs(np.log2(phi)) > 1).sum()),
           "sd_log_phi_per_group": {g: float(np.log(v).std()) for g, v in per.items()}}
    return phi, per, rec


def _prepared(ctx) -> dict:
    """The fold's factors and its reference in t3_decompose's form, made once per loaded fold (none depends on the
    seed, and the harness hands every seed the same input objects)."""
    roles = sorted(ctx.inputs)
    objs = tuple(ctx.inputs[r] for r in roles)
    hit = _PREP.get("fold")
    if hit is None or len(hit["inputs"]) != len(objs) or any(
            a is not b for a, b in zip(hit["inputs"], objs, strict=True)):
        ref_role = ctx.task.reference_role
        pool = ctx.inputs[ref_role]
        others = [r for r in roles if r != ref_role]
        ref_means, ref_sizes = group_means(pool)
        trend = trend_at([group_means(ctx.inputs[r])[0] for r in others], [ctx.task.role_days[r] for r in others],
                         ctx.task.role_days[ref_role])
        phi, per, rec = revert_factors(ref_means, ref_sizes, trend)
        rules = harness._t3_head().load_rules()
        rec |= {"other_samples": [ctx.task.role_stems[r] for r in others], "reference": ctx.task.role_stems[ref_role]}
        _PREP.clear()
        hit = _PREP["fold"] = {
            "inputs": objs, "phi": phi, "per": per, "record": rec,
            "phi_perm": harness.t3_decile_permutation(phi, harness.t3_mean_expression(pool), PERM_SEED),
            "group": np.asarray(pool.obs["group"]).astype(str),
            "ref": harness.t3_reference(pool, PRESET[ctx.task.role_stems[ref_role]], rules)}
    return hit


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = _prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    n = harness.count_for(ctx.board, pool.n_obs)
    rng = bs.t3d.stable_rng(ctx.seed, STREAM)
    if arm == MAIN:
        built = bs.construct(ref, n, rng, f=ref.V.p, phi=prep["phi"], w=EXPONENT)
    elif arm == CONTROL:
        built = bs.construct(ref, n, rng, f=ref.V.p, phi=prep["phi_perm"], w=EXPONENT)
    else:
        built = bs.construct(ref, n, rng, f=ref.V.p)
        if arm == GROUPED:
            if not np.all(built.source == "wt"):
                raise RuntimeError("the grouped arm expects every drawn cell to come from the reference pool")
            C, P = ref.pool.C[built.row].copy(), ref.pool.P[built.row]
            grp = prep["group"][built.row]
            rng_g = bs.t3d.stable_rng(ctx.seed, f"{STREAM}:grouped")
            for g in sorted(set(grp)):
                rows = np.flatnonzero(grp == g)
                phi = prep["per"].get(g, prep["phi"]) ** EXPONENT
                C[rows] = bs.t3d.apply_gene_factor(C[rows], phi, P[rows], rng_g)
            built = bs.Built(X=bs.t3d.lognorm(C), lab=built.lab, source=built.source, row=built.row)
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params={"rule": WT_REVERT_RULE if arm != BASE else "no factor", "stream": STREAM,
                                         "exponent": EXPONENT if arm != BASE else None,
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
    res["task_source"] = "experiments/t3_decomp/k7_wt_revert.py"
    res["factor_record"] = _PREP["fold"]["record"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


def _paired(rows: list[dict], x: str, y: str) -> dict:
    d = np.array([r["scores"][x] - r["scores"][y] for r in rows])
    return {"mean": float(d.mean()), "sd": float(d.std(ddof=1)), "ahead": int((d > 0).sum()), "n": int(len(d)),
            "passes": bool(d.mean() > MIN_GAIN and (d > 0).sum() >= MIN_AHEAD)}


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
            "paired": {f"{x} - {y}": _paired(rows, x, y) for x, y in (
                (MAIN, BASE), (MAIN, harness.FLOOR), (MAIN, CONTROL), (GROUPED, BASE), (GROUPED, MAIN),
                (CONTROL, BASE))},
            "factor_record": r.get("factor_record")}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if why:
        return {"void": True, "why": why}
    need = [f"{MAIN} - {y}" for y in (BASE, harness.FLOOR, CONTROL)]
    supported = all(per[n]["paired"][k]["passes"] for n in TASKS for k in need)
    grouped = supported and all(per[n]["paired"][f"{GROUPED} - {MAIN}"]["passes"] for n in TASKS)
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "k7a_supported": bool(supported),
            "k7b_replaces_k7a": bool(grouped),
            "verdict": ("supported: a priced request goes to the captain" if supported
                        else "not supported: closed, no request")}


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K7a-K7c (commit {res['commit']}): {res['verdict']}"
             + ("; the grouped factor replaces the pooled one" if res["k7b_replaces_k7a"] else "")]
    for n, t in TASKS.items():
        v = res["tasks"][n]
        lines.append(f"{n} {t.key}: " + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})"
                         + ("  passes" if d["passes"] else ""))
        for a in (harness.FLOOR, *ARMS):
            lines.append(f"    {a:<15} severity " + ", ".join(f"{s} {c}" for s, c in v["severity_states"][a].items())
                         + "; skills " + ", ".join(f"{m} {x:.3f}" for m, x in v["skills"][a].items()))
        f = v["factor_record"] or {}
        lines.append(f"    factor: groups {len(f.get('groups', []))}, SD of log factor "
                     f"{f.get('sd_log_phi', float('nan')):.3f}, genes beyond twofold "
                     f"{f.get('genes_beyond_twofold')}, at the clip {f.get('genes_at_clip')}")
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
