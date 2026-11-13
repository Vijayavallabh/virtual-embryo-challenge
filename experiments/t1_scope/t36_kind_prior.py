#!/usr/bin/env python
"""Try T36 (experiments/harness/TRIES.md, registered before the run): the strong-null down-tail prior of tries T12a and
T12b on the library-3 base, in try T27's zero-gap geometry. A labelled diagnostic for a ruling; no upload changes.

Why. The captain left the prior out on 2026-11-10 (decision log (l)). Diagnostic N28 found that the library-3 base's own
change against E10.5 as released does not carry the prior's genes, and its row says a counted try in T27's geometry may
follow. T36 is that try.

Geometry: T27's, imported unchanged from kind_gap.py (its task, SPLIT_SEED, reference, truth half, build seeds and its
lib3_base arm; the T1 board emulator, bake protocol, build main, computed anchors, the floor checked at 50 per seed).

Arms (n cells by the count rule, real cells, none from the truth half):
  lib3_base     kind_gap.py's own arm, draw for draw (the base, not a try): it must reproduce T27's committed per-seed
                scores (experiments/t1_scope/results/t27_kind_gap/kind_gap.json) to REPRO_TOL, or the reading is void
  lib3_base_sn  T36: lib3_base's drawn cells, the same draw, with the prior applied as the harness candidate
                b0_strong_null applies it on a T1 pseudo-task (vec.t1_scope.build_submission's down_prior step):
                vec.t1_scope.apply_down_prior with harness._t12_prior (DownPriorTable.level_only(), every DP_* constant
                unchanged, covered cells: those whose label harness._t12_labels takes from the task's reference),
                pb_ref the mean of the reference's log1p CP10k over the board's panel (every E10.5 cell outside the
                truth half), held: the pre-registered held-flat list b0_strong_null is handed (harness.DERIVED
                "t1_held_flat"), seed [build seed, 4] as build_submission passes it

The held-flat list. b0_strong_null holds the list its momentum step holds (computed on E8.5 and E9.5, committed before
E10.5 reached us). The harness hands that list only to a pseudo-task whose inputs include E8.5 and E9.5
(harness.derived_ok); T27's task has none, so the list is carried here as a fixed set of names that says nothing of
E10.5 (`held_flat` refuses one computed from E10.5 or later). The main session confirms this before the run.

Reading (fixed before the run, the adoption rule's constants vec.adoption.MIN_GAIN and MIN_AHEAD): "the prior adds in
the kind geometry" if lib3_base_sn beats lib3_base by a paired mean gain above MIN_GAIN with at least MIN_AHEAD of the
10 seeds ahead; otherwise "closed". A zero-gap, in-kind task is not adoption evidence either way: a pass leads only to a
priced question to the captain about his ruling of 2026-11-10 (l), and "closed" ends the line.
Reported: both arms' means, the paired difference, per-metric skill differences, and how many genes the prior promoted
per seed (a number of genes, never of cells), with the prior's own target-free checks.
Void if a seed or an arm is missing or non-finite, the run's tree was dirty or its commit unknown, T27's stored result
is missing a seed, the split seed differs, or lib3_base differs from T27's stored per-seed score by more than REPRO_TOL.

  run   --out DIR             the ten seeds; writes t36.json (refuses a folder that already holds anything)
  read  --runs DIR --out DIR  applies the reading; writes reading.json and reading.txt (refuses likewise)

Usage (repo root of the campaign checkout; about 15 minutes):
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/t1_scope/t36_kind_prior.py run --out runs/t1_t36/<stamp>
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kind_gap as kg  # noqa: E402  (T27's geometry, draws and scoring, unchanged)

from vec import harness, t1_scope  # noqa: E402
from vec.adoption import MIN_AHEAD, MIN_GAIN, SEEDS  # noqa: E402
from vec.boards import RANKED_METRICS  # noqa: E402
from vec.evaluate import SUBSAMPLE_FRAC, Stage, score_board, subsample  # noqa: E402
from vec.prep import log1p_cp10k  # noqa: E402

TRY = "T36"
BASE = "lib3_base"
PRIOR_ARM = "lib3_base_sn"
ARMS = (BASE, PRIOR_ARM)
PRIOR_SEED_STREAM = 4  # build_submission's seed=[seed, 4] for apply_down_prior
REPRO_TOL = 1e-6
VERDICTS = ("the prior adds in the kind geometry", "closed")
HERE = Path(__file__).resolve().parent
T27_RESULT = HERE / "results" / "t27_kind_gap" / "kind_gap.json"
REF_STAGE = "T1/E10.5_RNA.h5ad"  # the stage of T27's task (truth and reference both inside it)
NOTE = ("A zero-gap, in-kind task is not adoption evidence either way: a pass leads only to a priced question to the "
        "captain about his ruling of 2026-11-10 (l), and \"closed\" ends the line.")


def held_flat() -> tuple[str, ...]:
    """The held-flat genes b0_strong_null is handed (harness.DERIVED["t1_held_flat"], loaded as the harness loads it).
    Refuses unless every stage the list was computed from is earlier than E10.5, the stage that holds both T27's truth
    and its reference (the list then carries nothing about the truth half)."""
    d = harness.DERIVED["t1_held_flat"]
    late = [s for s in d.sources() if harness.STAGE_DAYS[s] >= harness.STAGE_DAYS[REF_STAGE]]
    if late:
        raise ValueError(f"the held-flat list was computed from stages not earlier than E10.5: {late}")
    return tuple(str(g) for g in d.load(d.path))


def covered_prior(ref_labels) -> t1_scope.DownPrior:
    """b0_strong_null's prior (harness._t12_prior) for a task whose reference cells carry `ref_labels`."""
    ctx = SimpleNamespace(task=SimpleNamespace(reference_role="ref"),
                          inputs={"ref": SimpleNamespace(obs=pd.DataFrame({"celltype": np.asarray(ref_labels)}))})
    return harness._t12_prior(ctx)


def prepare(data_dir: Path) -> dict:
    """kind_gap.prepare at T27's settings (truth_depth 1), plus what the prior needs from the same loaded file: the
    gene names, the cell labels, the prior, the held genes and pb_ref."""
    seen = {}
    real = kg.load_stage

    def capture(*args, **kw):
        seen["a"] = real(*args, **kw)
        return seen["a"]

    kg.load_stage = capture
    try:
        d = kg.prepare(data_dir)
    finally:
        kg.load_stage = real
    a = seen["a"]
    d.update(extras(d, list(map(str, a.var_names)), np.asarray(a.obs["celltype"]).astype(str)))
    return d


def extras(d: dict, var_names: list[str], labels: np.ndarray) -> dict:
    board = d["board"]
    genes = [str(g) for g in board.genes]
    cols = pd.Index(var_names).get_indexer(genes)
    if (cols < 0).any():
        raise KeyError(f"{int((cols < 0).sum())} panel genes are missing from the file (names not shown)")
    Xref = d["X"][d["ref_rows"]]
    # pb_ref as build_submission takes it from its `ref`: the reference's mean log1p CP10k per panel gene
    pb_ref = np.asarray(sparse.csr_matrix(Xref[:, cols], dtype=np.float64).mean(axis=0)).ravel()
    prior = covered_prior(labels[d["ref_rows"]])
    return {"t36_var_names": list(var_names), "t36_labels": labels, "t36_genes": genes, "t36_pb_ref": pb_ref,
            "t36_prior": prior, "t36_held": held_flat()}


def population(d: dict, rows: np.ndarray) -> ad.AnnData:
    """The drawn cells as the prior's step reads them: the file's log1p CP10k in X, its counts in layers, the labels."""
    rows = np.asarray(rows)
    pop = ad.AnnData(X=sparse.csr_matrix(d["X"][rows]),
                     obs=pd.DataFrame({t1_scope.LABEL_KEY: d["t36_labels"][rows]},
                                      index=[f"c{i}" for i in range(len(rows))]),
                     var=pd.DataFrame(index=pd.Index(d["t36_var_names"])))
    pop.layers["counts"] = sparse.csr_matrix(d["C"][rows]).astype(np.int32)
    return pop


def prior_arm(d: dict, rows: np.ndarray, build_seed: int, prior: t1_scope.DownPrior | None = None):
    """lib3_base's rows with the prior applied as build_submission applies it. Returns (the edited population, in the
    file's gene order, the prior's report, the covered mask). `prior` defaults to b0_strong_null's."""
    prior = d["t36_prior"] if prior is None else prior
    pop = population(d, rows)
    labels = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy()
    if prior.covered_labels is None:
        raise ValueError("b0_strong_null's prior names its covered labels")
    cov = np.isin(labels, sorted(prior.covered_labels))
    out, rep = t1_scope.apply_down_prior(pop, d["t36_pb_ref"], list(d["t36_genes"]), cov, list(d["t36_held"]), prior,
                                         seed=[build_seed, PRIOR_SEED_STREAM])
    if list(out.var_names) != list(d["t36_var_names"]):
        raise ValueError("the prior's output changed the gene order")
    return out, rep, cov


def arm_x(arm: str, n: int, seed: int) -> tuple[sparse.csr_matrix, dict]:
    """The arm's prediction (the file's gene order) for evaluation seed `seed`, and what is recorded about it."""
    d = kg._D
    build_seed = harness.BUILD_SEED_OFFSET + seed
    rows = kg.arm_rows(BASE, n, build_seed)
    if arm == BASE:
        return d["X"][rows], {}
    if arm != PRIOR_ARM:
        raise ValueError(f"unknown arm {arm}")
    out, rep, cov = prior_arm(d, rows, build_seed)
    cov_rows = rows[cov]
    gap = sparse.csr_matrix(log1p_cp10k(d["C"][cov_rows]) - d["X"][cov_rows])
    info = {"n_promoted_genes": int(rep["n_promoted"]), "budget": int(rep["budget"]),
            "covered_share": float(cov.mean()), "checks": rep["checks"],
            # reported, deciding nothing: how far the file's X of the covered cells sits from log1p CP10k of their
            # counts, which the prior's step writes for those cells (0 when the file stores exactly that)
            "recompute_max_abs": float(abs(gap).max()) if gap.nnz else 0.0}
    return sparse.csr_matrix(out.X, dtype=np.float32), info


def score_seed(seed: int) -> dict:
    d = kg._D
    t0 = time.time()
    X, board = d["X"], d["board"]
    target, ref = Stage(X=d["x_truth"]), Stage(X=X[d["ref_rows"]])
    f = score_board(board, subsample(ref, SUBSAMPLE_FRAC, seed), target, ref, seed=seed, anchors="computed",
                    build=harness.BUILD, protocol=harness.PROTOCOL)
    if abs(f.score - 50.0) > 1e-9 or set(RANKED_METRICS["T1"]) - set(f.anchors):
        raise ValueError(f"seed {seed}: the floor scores {f.score}")
    n = harness.count_for(board, len(d["ref_rows"]))
    row = {"seed": seed, "scores": {}, "skills": {}}
    for arm in ARMS:
        pred, info = arm_x(arm, n, seed)
        r = score_board(board, Stage(X=pred), target, ref, seed=seed, anchors=f.anchors, build=harness.BUILD,
                        protocol=harness.PROTOCOL)
        row["scores"][arm], row["skills"][arm] = float(r.score), {m: float(v) for m, v in r.skills.items()}
        if info:
            row["prior"] = info
    row["seconds"] = round(time.time() - t0, 1)
    print(f"seed {seed}: " + ", ".join(f"{k} {v:.2f}" for k, v in row["scores"].items())
          + f"; genes promoted {row['prior']['n_promoted_genes']}", flush=True)
    return row


def _state() -> tuple[str, bool]:
    commit = harness._git("rev-parse", "HEAD")
    return commit, commit != "unknown" and harness._git("status", "--porcelain", "--untracked-files=no") == ""


def _refuse_full(out: Path) -> None:
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise SystemExit(f"{out} already holds something: T36 never overwrites (a rerun needs a new folder and row)")


def cmd_run(a) -> None:
    out = Path(a.out)
    _refuse_full(out)
    commit, clean = _state()
    data = Path(a.data)
    kg._D.update(prepare(data))
    d = kg._D
    if a.jobs > 1:
        with mp.get_context("fork").Pool(a.jobs) as pool:
            rows = pool.map(score_seed, list(SEEDS))
    else:
        rows = [score_seed(s) for s in SEEDS]
    rows = sorted(rows, key=lambda r: r["seed"])
    commit_end, clean_end = _state()
    res = {"try": TRY, "commit": commit, "tree_clean": bool(clean and clean_end and commit == commit_end),
           "file": kg.FILE, "split_seed": kg.SPLIT_SEED, "seeds": list(SEEDS), "arms": list(ARMS),
           "protocol": harness.PROTOCOL, "build": harness.BUILD, "frac": SUBSAMPLE_FRAC,
           "inputs": {kg.FILE: harness._sha256(data / kg.FILE),
                      "experiments/t1_scope/kind_gap.py": harness._sha256(Path(kg.__file__)),
                      str(harness.DERIVED["t1_held_flat"].path.relative_to(harness.REPO_ROOT)):
                          harness._sha256(harness.DERIVED["t1_held_flat"].path)},
           "prior": {"table": "DownPriorTable.level_only()", "budget": d["t36_prior"].budget,
                     "band": d["t36_prior"].band, "f_min": d["t36_prior"].f_min, "groups": list(harness.T12_GROUPS),
                     "held_flat_sources": list(harness.DERIVED["t1_held_flat"].sources()),
                     "seed": f"[BUILD_SEED_OFFSET + seed, {PRIOR_SEED_STREAM}]"},
           "rows": rows, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out.mkdir(parents=True, exist_ok=True)
    path = out / "t36.json"
    path.write_text(json.dumps(res, indent=1))
    print(f"wrote {path}", flush=True)


def _finite(x) -> bool:
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def reading(run: dict | None, t27: dict | None) -> dict:
    """T36's reading of one run against T27's stored result (module docstring)."""
    why = []
    if run is None:
        why.append("the run's t36.json is missing")
    if t27 is None:
        why.append("T27's stored result is missing")
    if why:
        return {"void": True, "why": why, "verdict": "none: the reading is void", "note": NOTE}
    by_seed = {r.get("seed"): r for r in run.get("rows", [])}
    if run.get("try") != TRY:
        why.append("the result is not T36's")
    if sorted(by_seed) != list(SEEDS) or len(run.get("rows", [])) != len(SEEDS):
        why.append("a seed is missing or repeated")
    for s, r in sorted(by_seed.items(), key=lambda kv: str(kv[0])):
        if any(a not in r.get("scores", {}) or a not in r.get("skills", {}) for a in ARMS) or "prior" not in r:
            why.append(f"seed {s}: an arm is missing")
        elif not all(_finite(r["scores"][a]) for a in ARMS) or not all(
                _finite(v) for a in ARMS for v in r["skills"][a].values()):
            why.append(f"seed {s}: a non-finite score or skill")
    if run.get("commit") in (None, "", "unknown"):
        why.append("the run's commit is unknown")
    if run.get("tree_clean") is not True:
        why.append("the run comes from a dirty tree")
    if run.get("split_seed") != kg.SPLIT_SEED or t27.get("split_seed") != kg.SPLIT_SEED:
        why.append("the split seed differs from T27's")
    old = {r.get("seed"): r.get("scores", {}).get(BASE) for r in t27.get("rows", [])}
    if any(not _finite(old.get(s)) for s in SEEDS):
        why.append("T27's stored result lacks a seed's lib3_base score")
    worst = None
    if not why:
        worst = max(abs(float(by_seed[s]["scores"][BASE]) - float(old[s])) for s in SEEDS)
        if worst > REPRO_TOL:
            why.append(f"lib3_base does not reproduce T27's stored scores (largest difference {worst:.3g})")
    if why:
        return {"void": True, "why": why, "verdict": "none: the reading is void", "note": NOTE,
                "reproduction_largest_difference": worst}
    rows = [by_seed[s] for s in SEEDS]
    base = np.array([r["scores"][BASE] for r in rows], dtype=float)
    sn = np.array([r["scores"][PRIOR_ARM] for r in rows], dtype=float)
    diff = sn - base
    mean, ahead = float(diff.mean()), int((diff > 0).sum())
    passes = bool(mean > MIN_GAIN and ahead >= MIN_AHEAD)
    metrics = list(rows[0]["skills"][BASE])
    return {"void": False, "why": [], "verdict": VERDICTS[0] if passes else VERDICTS[1], "passes": passes,
            "note": NOTE, "rule": {"min_gain": MIN_GAIN, "min_ahead": MIN_AHEAD, "seeds": len(SEEDS)},
            "reproduction_largest_difference": worst, "reproduction_tolerance": REPRO_TOL,
            "mean_score": {BASE: float(base.mean()), PRIOR_ARM: float(sn.mean())},
            "paired": {"mean": mean, "sd": float(diff.std(ddof=1)), "ahead": ahead, "behind": int((diff < 0).sum()),
                       "per_seed": {str(s): float(x) for s, x in zip(SEEDS, diff, strict=True)}},
            "skill_difference": {m: float(np.mean([r["skills"][PRIOR_ARM][m] - r["skills"][BASE][m] for r in rows]))
                                 for m in metrics},
            "genes_promoted_per_seed": {str(r["seed"]): int(r["prior"]["n_promoted_genes"]) for r in rows},
            "prior_checks_passed_per_seed": {str(r["seed"]): bool(r["prior"]["checks"].get("passed")) for r in rows},
            "commit": run["commit"]}


def text(res: dict) -> str:
    if res["void"]:
        return "Try T36: VOID: " + "; ".join(res["why"]) + "\n" + res["note"] + "\n"
    p = res["paired"]
    lines = [f"Try T36 (commit {res['commit']}): {res['verdict']} (rule: mean gain > {res['rule']['min_gain']:g} and "
             f"at least {res['rule']['min_ahead']} of {res['rule']['seeds']} seeds ahead)",
             res["note"],
             f"lib3_base against T27's stored scores: largest difference {res['reproduction_largest_difference']:.2g} "
             f"(tolerance {res['reproduction_tolerance']:g})",
             "mean board score: " + ", ".join(f"{k} {v:.2f}" for k, v in res["mean_score"].items()),
             f"lib3_base_sn - lib3_base: {p['mean']:+.3f} (sd {p['sd']:.3f}, ahead {p['ahead']}, behind {p['behind']})",
             "  per seed: " + ", ".join(f"{s} {v:+.3f}" for s, v in p["per_seed"].items()),
             "per-metric skill difference (mean over seeds): "
             + ", ".join(f"{m} {v:+.4f}" for m, v in res["skill_difference"].items()),
             "genes promoted per seed: " + ", ".join(f"{s} {v}" for s, v in res["genes_promoted_per_seed"].items()),
             "prior's target-free checks passed per seed: "
             + ", ".join(f"{s} {v}" for s, v in res["prior_checks_passed_per_seed"].items())]
    return "\n".join(lines) + "\n"


def cmd_read(a) -> None:
    out = Path(a.out)
    _refuse_full(out)
    run_path, t27_path = Path(a.runs) / "t36.json", Path(a.t27)
    run = json.loads(run_path.read_text()) if run_path.is_file() else None
    t27 = json.loads(t27_path.read_text()) if t27_path.is_file() else None
    res = reading(run, t27)
    out.mkdir(parents=True, exist_ok=True)
    (out / "reading.json").write_text(json.dumps(res, indent=1))
    (out / "reading.txt").write_text(text(res))
    print(text(res), end="")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", required=True)
    r.add_argument("--data", default=str(harness.DATA_DIR))
    r.add_argument("--jobs", type=int, default=5)
    d = sub.add_parser("read")
    d.add_argument("--runs", required=True)
    d.add_argument("--t27", default=str(T27_RESULT))
    d.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "read": cmd_read}[a.cmd](a)


if __name__ == "__main__":
    main()
