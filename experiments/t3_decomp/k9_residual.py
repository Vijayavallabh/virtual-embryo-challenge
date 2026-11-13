#!/usr/bin/env python
"""Tries K9a and K9b (experiments/harness/TRIES.md; specification fixed by the main session on 2026-11-13, before any
code): does adding to a knockout-derived factor only each lineage group's part of the wild-type reference's own
deviation beyond its pooled part ("the residual") help? A labelled reading for a ruling, NOT the adoption rule.

Why. T3's first-choice file applies one per-gene factor measured from the Gata4 knockout against WT E8.75 in the
non-Mesp1 compartments, so it already holds the pooled part of WT E8.75's own deviation. Tries K7
(k7_wt_revert.py) built that deviation from the released wild-type series alone, pooled (K7a) or per lineage group
(K7b); diagnostic N29 found the per-group factor far from the pooled one (F = 0.89). K9 adds only the per-group part
beyond the pooled one to a knockout-derived factor.

Tasks (diagnostic pseudo-tasks, checked by vec.harness.check_task; configs/pseudo_tasks.yaml unchanged). Each is one
of vec.harness.K4_FOLDS' folds, its definition copied from the configured fold, with the released wild-type series
added as inputs (roles s_<stage>, every stage of k7_wt_revert.FULL_PANEL but the reference's, as K7 adds them):
  x1   source the Mab21l2 knockout against WT E9.5; reference WT E8.75; truth Gata4 replicate 2
       (fold "T3/Mab21l2->Gata4rep2"; a counted cross-knockout fold)
  x2   source both Gata4 replicates pooled against WT E8.75; reference WT E9.5; truth the Mab21l2 knockout
       (fold "T3/Gata4->Mab21l2"; a counted cross-knockout fold)
  r12  source Gata4 replicate 1 against WT E8.75; reference WT E8.75; truth Gata4 replicate 2
       (fold "T3/Gata4rep1->rep2"; a replicate fold: reported, never counted)

Arms, all drawn from one random stream (K4's, vec.harness.K4_STREAM), so every arm edits the same drawn cells:
  raw           the base: the fold's knockout-derived factor exactly as the harness candidate t3_factor_raw computes
                and applies it (exponent 1; build_submission.construct and assemble; no depth, no composition, no
                stage mixing). Not a try.
  raw_res       K9a: the same cells and factor, each cell's factor multiplied per gene by
                rho_g = (per_g / phi_pool) ** RES_EXPONENT, per_g K7b's per-group revert factor of the reference for
                the cell's lineage group and phi_pool K7a's pooled one (k7.group_means, trend_at, revert_factors,
                unchanged); rho = 1 for a group without a per-group factor. phi_raw x rho_g is clipped to rules.yaml's
                factor clip and applied once, with t3_decompose.apply_gene_factor (thinning, Poisson top-up from the
                kNN-smoothed profile), continuing the build's own generator.
  raw_res_perm  K9b, the control: as raw_res with each rho_g permuted across genes within deciles of the reference's
                mean expression (vec.harness.t3_decile_permutation, basis t3_mean_expression of the reference, seed
                PERM_SEED; the same permutation for every group).

Reading, fixed in the specification (the adoption rule's constants: a paired mean gain above 0.3 with at least 8 of 10
seeds ahead). K9a is SUPPORTED if on x1 AND on r12 raw_res beats raw that way and beats raw_res_perm that way, and on x2
raw_res's paired mean against raw is not below -0.3. Then a priced request may go to the captain; anything else is
CLOSED with no request. Only x1 and x2 are counted cross-knockout folds; r12 is a replicate fold that never counts;
x2's reference is WT E9.5, where K7 found reverts flattered.

  run   --task x1|x2|r12 --out DIR    one task per process; start the three together (one commit, a clean tree)
  read  --runs DIR --out DIR          apply the reading; write reading.json and reading.txt
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

from vec import harness  # noqa: E402
from vec.adoption import SEEDS  # noqa: E402

RES_EXPONENT = 0.5  # K7's and ruling R5's exponent
PERM_SEED = 20261113
STREAM = harness.K4_STREAM  # the stream t3_factor_raw draws from
BASE, MAIN, CONTROL = "raw", "raw_res", "raw_res_perm"
ARMS = (BASE, MAIN, CONTROL)
MIN_GAIN, MIN_AHEAD = 0.3, 8
X2_FLOOR = -0.3  # x2: raw_res's paired mean against raw must not be below this
SERIES_PREFIX = "s_"
FOLDS = {"x1": "T3/Mab21l2->Gata4rep2", "x2": "T3/Gata4->Mab21l2", "r12": "T3/Gata4rep1->rep2"}
COUNTED = ("x1", "x2")  # the counted cross-knockout folds; r12 is a replicate fold that never counts
RESIDUAL_RULE = ("rho_g = (K7b's per-group revert factor / K7a's pooled revert factor) ** 0.5 for the cell's lineage "
                 "group (1 for a group without a per-group factor); phi_raw x rho_g clipped to rules.yaml's factor "
                 "clip; applied once with t3_decompose.apply_gene_factor")
_PREP: dict = {}


def series_files(reference: str) -> dict[str, str]:
    """The released wild-type series K7 adds for a reference file: role s_<stage> -> file, every stage of
    k7.FULL_PANEL but the reference's own."""
    stem = harness.stem(reference)
    if stem not in k7.FULL_PANEL:
        raise ValueError(f"the reference {reference} is not one of the stages {list(k7.FULL_PANEL)}")
    return {f"{SERIES_PREFIX}{s}": f"T2/{s}.h5ad" for s in k7.FULL_PANEL if s != stem}


def _task(name: str) -> harness.PseudoTask:
    """K4's fold, as configured, with the wild-type series added as inputs; a diagnostic (never counted)."""
    fold = next(t for ts in harness.load_tasks().values() for t in ts if t.key == FOLDS[name])
    clash = set(fold.inputs) & set(series_files(fold.reference))
    if clash:
        raise ValueError(f"{fold.key}: roles {sorted(clash)} are already inputs")
    return harness.check_task(replace(
        fold, key=f"T3k9/{name}", inputs=dict(fold.inputs) | series_files(fold.reference), adoption=False,
        notes=f"Tries K9: fold {fold.key} with the released wild-type series added for the reference's revert."))


TASKS = {n: _task(n) for n in FOLDS}
KEYS = frozenset(t.key for t in TASKS.values())
FOLD_OF = {t.key: FOLDS[n] for n, t in TASKS.items()}


# ------------------------------------------------------------------------------------------------ the factors
def residual_factors(per: dict[str, np.ndarray], phi_pool: np.ndarray, exponent: float | None = None
                     ) -> dict[str, np.ndarray]:
    """rho_g = (per_g / phi_pool) ** exponent (default RES_EXPONENT) for every group with a per-group factor."""
    e = RES_EXPONENT if exponent is None else exponent
    pool = np.asarray(phi_pool, dtype=np.float64)
    return {g: (np.asarray(v, dtype=np.float64) / pool) ** e for g, v in per.items()}


def combined_factor(phi_raw: np.ndarray, rho: np.ndarray | None, clip: tuple[float, float]) -> np.ndarray:
    """phi_raw x rho clipped to `clip`; rho None (a group without a per-group factor) is rho = 1."""
    phi_raw = np.asarray(phi_raw, dtype=np.float64)
    return np.clip(phi_raw if rho is None else phi_raw * rho, *clip)


def permuted_residuals(rho: dict[str, np.ndarray], basis: np.ndarray, seed: int = PERM_SEED
                       ) -> dict[str, np.ndarray]:
    """K9b: each group's rho permuted across genes within deciles of `basis`, by t3_decile_permutation at one seed,
    so every group gets the same permutation."""
    return {g: harness.t3_decile_permutation(v, basis, seed) for g, v in rho.items()}


def factor_classes(phi_raw, rho: dict[str, np.ndarray], drawn_groups: np.ndarray,
                   clip: tuple[float, float]) -> list[tuple[np.ndarray, np.ndarray]]:
    """(factor, cell positions) per distinct combined factor among the drawn cells: groups whose combined factors are
    equal share one class, cells kept in drawn order; classes ordered by their first group name (sorted)."""
    classes: list[tuple[np.ndarray, list[str]]] = []
    for g in sorted(set(drawn_groups.tolist())):
        phi = combined_factor(phi_raw, rho.get(g), clip)
        for fac, names in classes:
            if np.array_equal(fac, phi):
                names.append(g)
                break
        else:
            classes.append((phi, [g]))
    return [(fac, np.flatnonzero(np.isin(drawn_groups, names))) for fac, names in classes]


def arm_built(bs, ref, n: int, seed: int, f, arm: str, phi_raw: np.ndarray, rho: dict[str, np.ndarray] | None,
              group: np.ndarray | None, clip: tuple[float, float]):
    """One arm's cells, before assembly. raw: build_submission.construct with phi_raw (t3_factor_raw's build). The
    residual arms: the same construct call without a factor (the same draws from the same generator), then each cell's
    counts moved once by its combined factor through t3_decompose.apply_gene_factor, continuing that generator; when
    every drawn cell's combined factor is the same vector this is one call on all cells, i.e. raw's build exactly."""
    rng = bs.t3d.stable_rng(seed, STREAM)
    if arm == BASE:
        return bs.construct(ref, n, rng, f=f, phi=phi_raw)
    drawn = bs.construct(ref, n, rng, f=f)
    if not np.all(drawn.source == "wt"):
        raise RuntimeError("the residual arms expect every drawn cell to come from the reference pool")
    if len(group) != ref.pool.C.shape[0]:
        raise RuntimeError("the lineage groups do not line up with the reference pool's rows")
    C, P = ref.pool.C[drawn.row].copy(), ref.pool.P[drawn.row]
    for fac, rows in factor_classes(phi_raw, rho, group[drawn.row], clip):
        C[rows] = bs.t3d.apply_gene_factor(C[rows], fac, P[rows], rng)
    return bs.Built(X=bs.t3d.lognorm(C), lab=drawn.lab, source=drawn.source, row=drawn.row)


def residual_cosine(phi_raw: np.ndarray, rho: dict[str, np.ndarray], sizes: dict[str, int]) -> float | None:
    """The cosine between log(phi_raw) and the reference-cell-weighted mean over groups of log(rho_g) (weights the
    reference's group sizes as K7 counts them; a group without a factor adds zero and leaves the cosine unchanged)."""
    if not rho:
        return None
    w = np.array([float(sizes[g]) for g in rho])
    m = sum(wi * np.log(rho[g]) for wi, g in zip(w / w.sum(), rho, strict=True))
    a = np.log(np.asarray(phi_raw, dtype=np.float64))
    na, nm = float(np.linalg.norm(a)), float(np.linalg.norm(m))
    return None if na == 0 or nm == 0 else float(a @ m / (na * nm))


def revert_parts(pool, series: list, series_days: list[float], ref_day: float):
    """K7's factors for the reference `pool` from the series alone, through K7's own functions, unchanged:
    (phi_pool, per-group factors, the reference's group sizes)."""
    ref_means, ref_sizes = k7.group_means(pool)
    trend = k7.trend_at([k7.group_means(a)[0] for a in series], series_days, ref_day)
    phi_pool, per, _ = k7.revert_factors(ref_means, ref_sizes, trend)
    return phi_pool, per, ref_sizes


def _prepared(ctx) -> dict:
    """The fold's raw factor (t3_factor_raw's: harness.t3_source_frame's phi, which t3_k4_fit returns as `raw`), K7's
    factors of its reference, the residuals and the fold's reference in t3_decompose's form, once per loaded fold."""
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
        series = [r for r in roles if r.startswith(SERIES_PREFIX)]
        phi_pool, per, sizes = revert_parts(pool, [ctx.inputs[r] for r in series],
                                            [ctx.task.role_days[r] for r in series], ctx.task.role_days[ref_role])
        rho = residual_factors(per, phi_pool)
        clip = tuple(float(c) for c in rules["components"]["factor"]["clip"])
        _PREP.clear()
        hit = _PREP["fold"] = {
            "key": ctx.task.key, "inputs": objs, "phi_raw": phi_raw, "rho": rho,
            "rho_perm": permuted_residuals(rho, harness.t3_mean_expression(pool)), "clip": clip,
            "group": np.asarray(pool.obs["group"]).astype(str),
            "ref": harness.t3_reference(pool, fold["reference_preset"], rules),
            "record": {"fold": FOLD_OF[ctx.task.key], "reference": ctx.task.role_stems[ref_role],
                       "series": [ctx.task.role_stems[r] for r in series], "groups_with_a_factor": len(per),
                       "cosine_log_phi_raw_mean_log_rho": residual_cosine(phi_raw, rho, sizes)}}
    return hit


def _build(ctx, arm: str):
    bs = harness._t3_head()
    prep = _prepared(ctx)
    ref, pool = prep["ref"], ctx.inputs[ctx.task.reference_role]
    rho = None if arm == BASE else prep["rho" if arm == MAIN else "rho_perm"]
    built = arm_built(bs, ref, harness.count_for(ctx.board, pool.n_obs), ctx.seed, ref.V.p, arm, prep["phi_raw"], rho,
                      prep["group"], prep["clip"])
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{STREAM}:assemble"))
    return harness.finish(ctx, pool, X, coords if ctx.board.needs_coords else None)


def candidates() -> list:
    def one(arm):
        return harness.Candidate(arm, frozenset({"T3"}), (), lambda ctx, arm=arm: _build(ctx, arm),
                                 params=harness.K4_PARAMS | {
                                     "residual": RESIDUAL_RULE if arm != BASE else None,
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
    res["task_source"] = "experiments/t3_decomp/k9_residual.py"
    res["factor_record"] = _PREP["fold"]["record"]
    path.write_text(json.dumps(res, indent=1, default=float))
    print(f"wrote {path}", flush=True)


# ------------------------------------------------------------------------------------------------ the reading
def _paired(rows: list[dict], x: str, y: str) -> dict:
    d = np.array([r["scores"][x] - r["scores"][y] for r in rows])
    return {"mean": float(d.mean()), "sd": float(d.std(ddof=1)), "ahead": int((d > 0).sum()), "n": int(len(d)),
            "passes": bool(d.mean() > MIN_GAIN and (d > 0).sum() >= MIN_AHEAD)}


def decide(per: dict) -> tuple[bool, dict]:
    """The fixed reading on per-task paired comparisons: supported only when raw_res passes against raw and against
    raw_res_perm on x1 and on r12, and its paired mean against raw on x2 is not below X2_FLOOR."""
    checks = {f"{n}: {MAIN} - {y}": per[n]["paired"][f"{MAIN} - {y}"]["passes"]
              for n in ("x1", "r12") for y in (BASE, CONTROL)}
    checks[f"x2: {MAIN} - {BASE} not below {X2_FLOOR}"] = per["x2"]["paired"][f"{MAIN} - {BASE}"]["mean"] >= X2_FLOOR
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
            "counted_fold": name in COUNTED,
            "mean_score": {x: float(np.mean([q["scores"][x] for q in rows])) for x in arms},
            "paired": {f"{x} - {y}": _paired(rows, x, y) for x, y in (
                (MAIN, BASE), (MAIN, CONTROL), (CONTROL, BASE), (BASE, harness.FLOOR), (MAIN, harness.FLOOR))},
            "severity_states": {x: {s: sum(q["t3"][x]["severity"]["state"] == s for q in rows)
                                    for s in sorted({q["t3"][x]["severity"]["state"] for q in rows})} for x in arms},
            "skill_diff_raw_res_minus_raw": {m: float(np.mean([q["skills"][MAIN][m] - q["skills"][BASE][m]
                                                                for q in rows])) for m in metrics},
            "cosine_log_phi_raw_mean_log_rho": rec.get("cosine_log_phi_raw_mean_log_rho"),
            "groups_with_a_factor": rec.get("groups_with_a_factor")}
    if len(commits) > 1:
        why.append(f"results from different commits {sorted(commits)}")
    if why:
        return {"void": True, "why": why}
    supported, checks = decide(per)
    return {"void": False, "why": [], "commit": sorted(commits)[0], "tasks": per, "checks": checks,
            "k9a_supported": bool(supported),
            "verdict": ("SUPPORTED: a priced request may go to the captain" if supported
                        else "CLOSED: no request"),
            "label": LABEL}


LABEL = ("A labelled reading for a ruling, NOT the adoption rule: only x1 and x2 are counted cross-knockout folds; r12 "
         "is a replicate fold that never counts; x2's reference is WT E9.5, where K7 found reverts flattered.")


def _fmt(x) -> str:
    return "n/a" if x is None else f"{x:.3f}"


def text(res: dict) -> str:
    if res["void"]:
        return "VOID: " + "; ".join(res["why"]) + "\n"
    lines = [f"Tries K9a-K9b (commit {res['commit']}): {res['verdict']}", res["label"],
             "Checks: " + "; ".join(f"{k} {'yes' if v else 'no'}" for k, v in res["checks"].items()),
             "Reported, deciding nothing:"]
    for n, t in TASKS.items():
        v = res["tasks"][n]
        kind = "counted cross-knockout fold" if v["counted_fold"] else "replicate fold, never counted"
        lines.append(f"{n} {FOLD_OF[t.key]} ({kind}): " + ", ".join(f"{a} {s:.2f}" for a, s in v["mean_score"].items()))
        for k, d in v["paired"].items():
            lines.append(f"    {k}: {d['mean']:+.2f} (sd {d['sd']:.2f}, ahead {d['ahead']} of {d['n']})"
                         + ("  passes" if d["passes"] else ""))
        for a in (harness.FLOOR, *ARMS):
            lines.append(f"    {a:<13} severity " + ", ".join(f"{s} {c}" for s, c in v["severity_states"][a].items()))
        lines.append("    skills raw_res - raw: " + ", ".join(f"{m} {x:+.3f}"
                                                            for m, x in v["skill_diff_raw_res_minus_raw"].items()))
        lines.append(f"    cosine of log(phi_raw) with the weighted mean log(rho): "
                     f"{_fmt(v['cosine_log_phi_raw_mean_log_rho'])} (above 0.8: the residual is more of what the "
                     f"factor already holds); groups with a per-group factor {v['groups_with_a_factor']}")
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    res = reading(Path(a.runs))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for p in (out / "reading.json", out / "reading.txt"):
        if p.exists():
            sys.exit(f"{p} exists")
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
