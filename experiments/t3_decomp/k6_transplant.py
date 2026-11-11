#!/usr/bin/env python
"""Tries K6a, K6b and K6c (experiments/harness/TRIES.md, registered before the run): real cells of the sibling
replicate's run inside a T3 prediction, on the Gata4 replicate fold. A labelled diagnostic for a ruling: the replicate
fold never counts for adoption, and no upload changes whatever it shows.

The pre-registered script (t3_decompose.py, hash-pinned, not edited) builds every borrowed variant from wild-type
E8.75 cells moved by the sibling replicate's factor and depth ("xrep:F+depth", the stand-in of T3's s1). Its
"xrep:draw" row (the sibling's own cells at the sibling's composition) scores about 8 points higher, and its "xrep:sys"
row (the sibling's non-Mesp1 cells beside untouched wild-type cells, at the sibling's coverage) about 8 lower. No row
separates what the sibling's real cells add where the knockout cannot act from what they add as the same genotype.
This script scores three more rows with the frozen script's own functions, context, protocol and scorer:

  xrep:hyb     K6a. Strict non-Mesp1 labels: the sibling replicate's real cells. Every other label: wild-type cells
               moved by the sibling's factor and per-label depth, as xrep:F+depth moves them. Wild-type composition.
  xrep:within  K6b. Every label: the sibling replicate's real cells, at the wild-type composition (the sibling-run
               form of oracle:ko_within). It carries the shared genotype, so it bounds K6a from above and says nothing
               about another knockout.
  xrep:sys@p   K6c, the control of K6a: K6a without the factor and depth on the wild-type part.

It also rebuilds xrep:F+depth, xrep:draw and xrep:sys at the same seeds: the reproduction gate compares them with the
committed experiments/t3_decomp/results/gata4/gata4_records.csv before any new row is read.

Run from the repository root (about 10 minutes on 8 workers):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
    uv run --frozen --no-sync python -I -B experiments/t3_decomp/k6_transplant.py --out runs/t3_k6/<stamp>
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))

import t3_decompose as t3d  # noqa: E402  (the pre-registered script, imported as tests import it)

DATA = REPO / "data" / "raw"
COMMITTED = HERE / "results" / "gata4" / "gata4_records.csv"
REPRO = ("xrep:F+depth", "xrep:draw", "xrep:sys")
NEW = ("xrep:hyb", "xrep:within", "xrep:sys@p")
REPRO_TOL = 0.05  # board points per (target, seed) row, as the reproduction gates of TRIES rows T20-T25
_CTX = None


def variants(ctx, ti: int):
    """(name, builder(n, seed)) for the three reproduced rows and the three new ones, target ti's sibling as source."""
    V, wt = ctx.V, ctx.wt
    S = ctx.targets[1 - ti].full
    has = np.array([len(S.pool.idx[t]) > 0 for t in range(V.L)])
    nm = V.nm & has  # strict non-Mesp1 labels the sibling's pool holds; any other label is built from wild-type cells
    specs = []

    def add(name, fn):
        specs.append((name, fn))

    add("xrep:F+depth", lambda n, seed: t3d.build(ctx, n, t3d.stable_rng(seed, "xrep:F+depth"), f=V.p, phi=S.phi_nm,
                                                  ratio=S.ratio))
    add("xrep:draw", lambda n, seed: t3d.build(ctx, n, t3d.stable_rng(seed, "xrep:draw"), f=S.q_hard, src=S.pool))
    add("xrep:sys", lambda n, seed: t3d.build(ctx, n, t3d.stable_rng(seed, "xrep:sys"), f=S.f["top"]["cov"],
                                              src_by_label=lambda t: S.pool if V.nm[t] else wt.pool))

    f_nm, f_rest = np.where(nm, V.p, 0.0), np.where(nm, 0.0, V.p)
    share_nm = float(f_nm.sum() / (f_nm.sum() + f_rest.sum()))

    def hyb(n, seed, moved: bool):
        n_nm = int(round(n * share_nm))
        name = "xrep:hyb" if moved else "xrep:sys@p"
        parts = []
        if n_nm > 0:
            parts.append(t3d.build(ctx, n_nm, t3d.stable_rng(seed, f"{name}:nm"), f=f_nm, src=S.pool))
        if n - n_nm > 0:
            kw = {"phi": S.phi_nm, "ratio": S.ratio} if moved else {}
            parts.append(t3d.build(ctx, n - n_nm, t3d.stable_rng(seed, f"{name}:rest"), f=f_rest, **kw))
        return np.vstack(parts)

    add("xrep:hyb", lambda n, seed: hyb(n, seed, True))
    pm = np.where(has, V.p, 0.0)
    add("xrep:within", lambda n, seed: t3d.build(ctx, n, t3d.stable_rng(seed, "xrep:within"), f=pm / pm.sum(),
                                                 src=S.pool))
    add("xrep:sys@p", lambda n, seed: hyb(n, seed, False))
    return specs, {"share_nm": share_nm, "nm_labels": [lab for lab, b in zip(V.labels, nm, strict=True) if b],
                   "labels_absent_from_sibling": [lab for lab, b in zip(V.labels, has, strict=True) if not b]}


def score_seed(ctx, ti: int, seed: int) -> list[dict]:
    """t3_decompose.score_seed's protocol, anchors and record, for this script's variants only."""
    from common.core_metrics import SEVERITY_LOG_WORST, de_genes

    from vec.boards import RANKED_METRICS
    from vec.evaluate import Stage, aggregate, panel_metrics, skills_from_anchors, split_half

    t0 = time.time()
    T, pr = ctx.targets[ti], ctx.rules["protocol"]
    rows_t = t3d.working_copy_rows(len(T.X_rel), pr["working_copy_fraction"], seed)
    rows_r = t3d.working_copy_rows(len(ctx.ref_X), pr["working_copy_fraction"], seed + 1)
    truth, other = split_half(Stage(X=T.X_rel[rows_t]), seed=pr["split_seed"])
    ref = Stage(X=ctx.ref_X[rows_r])
    floor_v = panel_metrics("T3", ref, truth, ref, seed=seed)
    ceil_v = panel_metrics("T3", other, truth, ref, seed=seed)
    anchors = {m: {"floor": floor_v.get(m), "ceiling": ceil_v.get(m)} for m in RANKED_METRICS["T3"]}
    anchors = {m: a for m, a in anchors.items() if a["floor"] is not None and a["ceiling"] is not None}
    up, dn, _ = de_genes(truth.X, ref.X)
    n_pred = pr["prediction_cells_multiplier"] * max(truth.n, ref.n)
    floor_cut = SEVERITY_LOG_WORST + pr["severity_floor_tolerance"]
    recs = []
    specs, _ = variants(ctx, ti)
    for name, builder in specs:
        X = builder(n_pred, seed)
        vals = panel_metrics("T3", Stage(X=X), truth, ref, seed=seed)
        skills = skills_from_anchors("T3", vals, anchors)
        score, groups = aggregate("T3", skills)
        sev = vals.get("severity_slope")
        recs.append({"target": T.name, "seed": seed, "variant": name, "n_pred": int(X.shape[0]), "t3": score,
                     **{f"skill_{m}": s for m, s in skills.items()}, **{f"group_{g}": s for g, s in groups.items()},
                     **{k: vals.get(k) for k in ("de_score", "de_direction", "severity_slope", "mmd_u", "variogram")},
                     "sev_floor": bool(sev is None or sev <= floor_cut), **t3d._diagnostics(X, truth.X, ref.X, up, dn)})
    print(f"[{T.name} seed {seed}] {len(recs)} variants, {time.time() - t0:.0f}s", flush=True)
    return recs


def _task(task):
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=_CTX.threads):
        return score_seed(_CTX, *task)


def reproduction(rec: pd.DataFrame) -> dict:
    """The gate: every rebuilt (target, seed) row of the three frozen variants against the committed records."""
    old = pd.read_csv(COMMITTED)
    old = old[old["variant"].isin(REPRO)][["target", "seed", "variant", "t3"]].rename(columns={"t3": "t3_committed"})
    new = rec[rec["variant"].isin(REPRO)][["target", "seed", "variant", "t3"]]
    m = new.merge(old, on=["target", "seed", "variant"], how="left")
    m["diff"] = (m["t3"] - m["t3_committed"]).abs()
    return {"rows": int(len(m)), "missing": int(m["t3_committed"].isna().sum()), "max_abs_diff": float(m["diff"].max()),
            "tolerance": REPRO_TOL, "passed": bool(m["t3_committed"].notna().all() and m["diff"].max() <= REPRO_TOL),
            "by_variant": {v: float(g["diff"].max()) for v, g in m.groupby("variant")}}


def paired(rec: pd.DataFrame, change: str, base: str) -> dict:
    out = {}
    skills = [c for c in rec.columns if c.startswith("skill_")]
    for target, g in rec.groupby("target"):
        a = g[g["variant"] == change].set_index("seed").sort_index()
        b = g[g["variant"] == base].set_index("seed").sort_index()
        d = a["t3"] - b["t3"]
        out[target] = {"mean_change": float(a["t3"].mean()), "mean_base": float(b["t3"].mean()),
                       "paired_mean": float(d.mean()), "paired_sd": float(d.std(ddof=1)), "ahead": int((d > 0).sum()),
                       "seeds": int(len(d)), "floored_change": int(a["sev_floor"].sum()),
                       "skill_gain": {s[6:]: float((a[s] - b[s]).mean()) for s in skills}}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", default=None, help="default: the registered seeds of rules.yaml (0-9)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args(argv)
    from threadpoolctl import threadpool_limits

    rules = yaml.safe_load(t3d.RULES_PATH.read_text())
    seeds = t3d.parse_seeds(a.seeds) if a.seeds else list(rules["protocol"]["seeds"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    paths = {"wt": str(DATA / "T2" / "E8.75.h5ad"),
             "targets": [("rep1", str(DATA / "T3" / "E8.75_gata4_ko_1.h5ad")),
                         ("rep2", str(DATA / "T3" / "E8.75_gata4_ko_2.h5ad"))],
             "stage:E9.5": str(DATA / "T2" / "E9.5.h5ad"), "stage:E8.25_late": str(DATA / "T2" / "E8.25_late.h5ad"),
             "series:E8.25_late": str(DATA / "T2" / "E8.25_late.h5ad"), "series:E9.5": str(DATA / "T2" / "E9.5.h5ad"),
             "series:E8.5": str(DATA / "T2" / "E8.5.h5ad"), "transfer_ko": None}
    n_par = a.workers * a.threads
    global _CTX
    with threadpool_limits(limits=n_par, user_api="blas"), threadpool_limits(limits=1, user_api="openmp"):
        ctx, _ = t3d.build_context("gata4", paths, rules, n_jobs=n_par, threads=a.threads, out=out)
        _CTX = ctx
        info = {T.name: variants(ctx, ti)[1] for ti, T in enumerate(ctx.targets)}
        tasks = [(ti, s) for ti in range(len(ctx.targets)) for s in seeds]
        recs = []
        if a.workers > 1:
            with mp.get_context("fork").Pool(min(a.workers, len(tasks))) as pool:
                for r in pool.imap_unordered(_task, tasks):
                    recs.extend(r)
        else:
            for t in tasks:
                recs.extend(_task(t))
    rec = pd.DataFrame(recs).sort_values(["target", "variant", "seed"])
    rec.to_csv(out / "k6_records.csv", index=False)
    gate = reproduction(rec)
    summary = {"seeds": seeds, "protocol_as_registered": seeds == list(rules["protocol"]["seeds"]),
               "git": t3d.git_state(), "sha256": {"rules": t3d.sha256(t3d.RULES_PATH),
                                                  "t3_decompose": t3d.sha256(HERE / "t3_decompose.py"),
                                                  "script": t3d.sha256(Path(__file__))},
               "variants": info, "reproduction_gate": gate,
               "means": {v: {t: float(g["t3"].mean()) for t, g in gv.groupby("target")}
                         for v, gv in rec.groupby("variant")},
               "K6a_xrep:hyb_vs_xrep:F+depth": paired(rec, "xrep:hyb", "xrep:F+depth"),
               "K6b_xrep:within_vs_xrep:F+depth": paired(rec, "xrep:within", "xrep:F+depth"),
               "K6c_xrep:sys@p_vs_xrep:F+depth": paired(rec, "xrep:sys@p", "xrep:F+depth"),
               "K6a_minus_K6c": paired(rec, "xrep:hyb", "xrep:sys@p"),
               "K6b_minus_K6a": paired(rec, "xrep:within", "xrep:hyb")}
    (out / "k6_summary.json").write_text(json.dumps(t3d.json_safe(summary), indent=1))
    print(json.dumps(t3d.json_safe({k: summary[k] for k in ("reproduction_gate", "means")}), indent=1))
    for k in ("K6a_xrep:hyb_vs_xrep:F+depth", "K6b_xrep:within_vs_xrep:F+depth", "K6c_xrep:sys@p_vs_xrep:F+depth"):
        for t, d in summary[k].items():
            print(f"{k} {t}: {d['paired_mean']:+.2f} (sd {d['paired_sd']:.2f}, ahead {d['ahead']} of {d['seeds']})")
    return 0 if gate["passed"] else 3


if __name__ == "__main__":
    sys.exit(main())
