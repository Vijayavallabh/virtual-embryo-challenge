#!/usr/bin/env python
"""Row N18 of experiments/harness/TRIES.md (T2 embryo; a diagnostic, not a try): the target-free fidelity gate of try
E14 (`fgw_latent`, vec.t2_latent). It runs before any counted run, and on a fail the try stops.

Can the count model reproduce the cells of the stages it was trained on, at all? On the two brackets of the E7.5
backtest (T2/E7.25.h5ad and T2/E8.0.h5ad, loaded as the harness loads them; the backtest's own panel), 20% of each
bracket's cells (seeded, one split per bracket) are held out of training. The negative-binomial VAE (vec.t2_latent,
trained once on the other 80% of both brackets, as try E14 will train on its inputs) encodes the training cells; NB
counts are sampled at t = 0 (the encoder means of the before-bracket's training cells) and at t = 1 (those of the
after-bracket's), each at the training cell's own panel library, and the sampled cells are compared with the held-out
real cells of that bracket on the scorer's two cell-state metrics (the pinned scorer's
`common.core_metrics.mmd_unbiased` and `variogram_score`, called as its T2 panel calls them: seed only, here the repeat
number), with a real-versus-real control (the held-out cells' two halves) and the distance between the brackets (real
cells of the other bracket's training set against the same held-out half).

For each bracket and metric, over REPEATS independent splits of the held-out cells into halves Ha and Hb: g =
metric(generated cells, Hb), the generated cells as many as Hb; c = metric(Ha, Hb); o = metric(the other bracket's
real cells, Hb). g, c and o are the means over the repeats. The rule, fixed in row N18 before this runs:
  variogram: g <= 1.3 c (MAX_VARIOGRAM_RATIO, the ratio to the real-versus-real control the design review proposed);
  mmd_u:     (g - c) <= 0.1 (o - c) (MAX_MMD_SHARE). mmd_u is an unbiased estimate of MMD^2: its real-versus-real
             control is centred on 0 and can be negative, so a ratio to it is not defined. The generated cells' excess
             over the control is measured instead as a share of the real distance between the brackets.
N18 passes if all four checks (two brackets, two metrics) hold; on a fail, try E14 does not run.

Target-free: it reads the backtest's two input stages only (E7.25 and E8.0, the brackets of the E7.5 backtest), not its
target (E7.5), and no score, board anchor or built file. E7.25 is also the target of the other counted embryo backtest
(E7.25 from E6.75 and E7.5); N18 uses it only to measure the generator at a bracket, and nothing it computes enters a
build. It prints and writes no cell count and no test-board figure: only metric values, ratios, settings and hashes.

Compute: CPU by default (--device cpu), 16 BLAS and torch threads; it trains one model (20-40 minutes at 16 threads)
and then scores (minutes). --device auto uses a GPU only if one is free (at least 90% of its memory free) and records
which device it used. The CPU path is the one the candidate's determinism check can pass. Needs the `gpu` dependency
group's environment (a separate venv outside the repository; docs/COMPUTE.md, UV_PROJECT_ENVIRONMENT); nothing here
changes the project's main environment.

It refuses to run while its output folder holds any file (a rerun needs a new TRIES row). Outputs, in
experiments/harness/results/n18_latent_fidelity/ (no target-derived content, no counts): n18_latent_fidelity.json and
n18_latent_fidelity.txt, each stamped with vec.provenance (commit, clean-tree flag, UTC, input hashes).

Usage (DGX, repo root, the gpu venv, CPU):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 CUDA_VISIBLE_DEVICES= \
      UV_PROJECT_ENVIRONMENT=$HOME/v/venvs/vec-gpu uv run --frozen --no-sync python -I -B \
      experiments/harness/n18_latent_fidelity.py --device cpu
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import textwrap
import time
from pathlib import Path

import numpy as np
from scipy import sparse

from vec import harness, provenance, t2_build, t2_latent
from vec.boards import load_boards
from vec.prep import load_stage

OUT = Path(__file__).resolve().parent / "results" / "n18_latent_fidelity"
TASK_KEY = "T2e/E7.25+E8.0->E7.5"  # the E7.5 backtest, the task whose two brackets the gate reads
EXPECTED_INPUTS = {"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"}  # never its target (T2/E7.5.h5ad)
BRACKETS = (("before", 0.0), ("after", 1.0))  # (role, t)
METRICS = ("mmd_u", "variogram")
HOLD_FRAC = 0.2
REPEATS = 5
HOLD_SEED = 20261109  # the one 80/20 split of each bracket
REPEAT_SEED = 7  # the halves, the generated rows and the other bracket's rows of each repeat
GEN_SEED = 8  # the NB draws of each repeat
MAX_VARIOGRAM_RATIO = 1.3  # row N18: variogram g / c at most this
MAX_MMD_SHARE = 0.1  # row N18: (g - c) / (o - c) of mmd_u at most this (see the docstring: the ratio to c is undefined)
GPU_FREE_FRACTION = 0.9  # --device auto: a GPU counts as free with at least this much of its memory free
NAMES = {"json": "n18_latent_fidelity.json", "text": "n18_latent_fidelity.txt"}  # write order
LABEL = "N18: the target-free fidelity gate of try E14 (fgw_latent), a diagnostic, not a try"


def metric_values(pred, truth, seed: int) -> dict[str, float]:
    """The pinned scorer's two cell-state metrics, called as its T2 panel calls them (`seed=` only; the panel rounds
    its values, this gate keeps them unrounded)."""
    from common.core_metrics import mmd_unbiased, variogram_score

    return {"mmd_u": float(mmd_unbiased(pred, truth, seed=seed)),
            "variogram": float(variogram_score(pred, truth, seed=seed))}


def holdout_split(n: int, seed) -> tuple[np.ndarray, np.ndarray]:
    """(train rows, held rows), each sorted: a seeded permutation of n rows, its first round(HOLD_FRAC x n) held out."""
    perm = np.random.default_rng(seed).permutation(int(n))
    k = int(round(HOLD_FRAC * n))
    return np.sort(perm[k:]), np.sort(perm[:k])


def repeat_rows(rng: np.random.Generator, held: np.ndarray, n_train: int, n_other_train: int):
    """One repeat's disjoint halves of the held rows (equal sizes), and rows into the bracket's training cells and the
    other bracket's training cells, as many as one half, drawn without replacement."""
    perm = rng.permutation(held)
    half = len(perm) // 2
    if half < 2:
        raise ValueError("too few held-out cells for a split-half control")
    if half > n_train or half > n_other_train:
        raise ValueError("the training cells are fewer than half the held-out cells")
    return (perm[:half], perm[half:2 * half], rng.choice(n_train, half, replace=False),
            rng.choice(n_other_train, half, replace=False))


def check_metric(metric: str, g: float, c: float, o: float) -> dict:
    """One bracket-metric check of N18's rule (module docstring; thresholds fixed in the row)."""
    base = {"generated": g, "control": c, "other_bracket": o}
    if metric == "variogram":
        ratio = g / c if c > 0 else float("inf")
        share = (g - c) / (o - c) if o > c else float("nan")
        return {**base, "ratio_to_control": ratio, "share_of_bracket_gap": share, "threshold": MAX_VARIOGRAM_RATIO,
                "statistic": "ratio_to_control", "pass": bool(ratio <= MAX_VARIOGRAM_RATIO)}
    if metric == "mmd_u":
        share = (g - c) / (o - c) if o > c else float("inf")
        return {**base, "share_of_bracket_gap": share, "threshold": MAX_MMD_SHARE,
                "statistic": "share_of_bracket_gap", "pass": bool(share <= MAX_MMD_SHARE)}
    raise ValueError(f"unknown metric {metric!r}")


def verdict(values: dict[str, dict[str, dict[str, float]]]) -> tuple[bool, dict]:
    """N18's verdict from each bracket's and metric's generated, control and other-bracket means (g, c, o): the checks
    (`check_metric`), and whether all of them pass."""
    checks = {stage: {m: check_metric(m, v["generated"], v["control"], v["other_bracket"]) for m, v in per.items()}
              for stage, per in values.items()}
    return all(ck["pass"] for per in checks.values() for ck in per.values()), checks


def pick_device(requested: str) -> str:
    """'cpu', 'cuda:<i>', or 'auto': the first GPU with at least GPU_FREE_FRACTION of its memory free, else the CPU."""
    if requested != "auto":
        return requested
    import torch

    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            free, total = torch.cuda.mem_get_info(i)
            if free >= GPU_FREE_FRACTION * total:
                return f"cuda:{i}"
    return "cpu"


def tree_state(exempt: list[Path]) -> dict:
    """HEAD and whether the tree is clean by vec.provenance's rule, the outputs exempt (None when git cannot tell)."""
    commit, unclean = provenance.git_state(provenance.REPO_ROOT, ignore=exempt)
    return {"commit": commit, "clean": None if unclean is None else not unclean}


def report_text(res: dict) -> str:
    lines = [f"{LABEL}: {res['verdict'].upper()}"]
    lines += textwrap.wrap(
        f"On the E7.5 backtest's two brackets ({EXPECTED_INPUTS['before']} at t = 0, {EXPECTED_INPUTS['after']} at t = "
        f"1), {HOLD_FRAC:.0%} of each bracket's cells held out of training; means over {REPEATS} repeats: generated "
        "cells, the real-versus-real control and the other bracket's real cells, each against the same held-out half. "
        f"Rule: variogram generated / control at most {MAX_VARIOGRAM_RATIO}; mmd_u (generated - control) / (other "
        f"bracket - control) at most {MAX_MMD_SHARE} (mmd_u's real-versus-real control is centred on 0, so no ratio to "
        "it is defined); all four checks must pass.", 118)
    lines.append("")
    for stage, per in res["checks"].items():
        for m, ck in per.items():
            lines.append(f"  {stage:<8s}{m:<11s}generated {ck['generated']:<12.6g}control {ck['control']:<12.6g}other "
                         f"bracket {ck['other_bracket']:<12.6g}{ck['statistic']} {ck[ck['statistic']]:.4g} "
                         f"(threshold {ck['threshold']:g}): {'pass' if ck['pass'] else 'FAIL'}")
    s = res["provenance"]
    lines += ["", f"Device {res['model']['device']}, {res['model']['threads']} threads; commit {s['commit']}, tree "
                  f"clean {s['tree_clean']}, UTC {s['utc']}",
              "Inputs: " + "; ".join(f"{k} {v['path']} sha256 {v['sha256']}" for k, v in res["inputs"].items()),
              f"verdict: {res['verdict']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(OUT), help="output folder (default: the registered one)")
    ap.add_argument("--data", default=str(harness.DATA_DIR), help="the released data (default: data/raw)")
    ap.add_argument("--device", default="cpu", help="cpu (default), cuda:<i>, or auto (a GPU only if one is free)")
    ap.add_argument("--threads", type=int, default=t2_latent.TRAIN_THREADS, help="torch CPU threads (default 16)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N18 has already written its outputs; a rerun needs a new TRIES row")
    paths = {k: out / v for k, v in NAMES.items()}
    exempt = [*paths.values(), *(provenance.stamp_path(p) for p in paths.values())]
    started = tree_state(exempt)
    t = next(x for x in harness.load_tasks()["T2:embryo:test_interp"] if x.key == TASK_KEY)
    if {r: t.inputs[r] for r in EXPECTED_INPUTS} != EXPECTED_INPUTS or t.target in EXPECTED_INPUTS.values():
        sys.exit(f"{TASK_KEY}: its brackets are not {EXPECTED_INPUTS}; N18 reads no other stage")
    genes = list(load_boards()[t.panel_board].genes)  # the backtest's own panel; no test-board figure is read
    stages = {r: load_stage(Path(a.data) / f, genes=genes, drop_duplicates=True) for r, f in EXPECTED_INPUTS.items()}
    device = pick_device(a.device)
    parts = {r: t2_build._stage_counts(stages[r], genes) for r in ("before", "after")}
    split = {r: holdout_split(parts[r][1].shape[0], [HOLD_SEED, i]) for i, (r, _) in enumerate(BRACKETS)}
    train_counts = sparse.vstack([parts[r][1][split[r][0]] for r, _ in BRACKETS]).tocsr()
    model = t2_latent.train_nbvae(train_counts, t2_latent.TRAIN_SEED, threads=a.threads, device=device)
    n_b = len(split["before"][0])
    if model.latent.shape[0] != train_counts.shape[0]:
        raise SystemExit("the model's latent means are not the training cells")
    latent = {"before": model.latent[:n_b], "after": model.latent[n_b:]}
    values: dict[str, dict] = {}
    repeats_log: dict[str, dict] = {}
    for i, (role, _) in enumerate(BRACKETS):
        other = "after" if role == "before" else "before"
        E, _, lib, other_counts = parts[role]
        Eo = parts[other][0]
        train, held = split[role]
        other_train = split[other][0]
        panel = lib - other_counts
        collected = {m: {"generated": [], "control": [], "other_bracket": []} for m in METRICS}
        for r in range(REPEATS):
            rng = np.random.default_rng([REPEAT_SEED, i, r])
            ha, hb, rows, orows = repeat_rows(rng, held, len(train), len(other_train))
            X_gen, _ = t2_latent.decode_and_sample(model, latent[role][rows], panel[train[rows]],
                                                   other_counts[train[rows]], np.random.default_rng([GEN_SEED, i, r]))
            truth, real, other_real = E[hb], E[ha], Eo[other_train[orows]]
            for name, pred in (("generated", X_gen), ("control", real), ("other_bracket", other_real)):
                v = metric_values(pred, truth, seed=r)
                for m in METRICS:
                    collected[m][name].append(v[m])
        values[role] = {m: {k: float(np.mean(v)) for k, v in collected[m].items()} for m in METRICS}
        repeats_log[role] = {m: {k: [float(x) for x in v] for k, v in collected[m].items()} for m in METRICS}
        print(f"{role}: " + ", ".join(f"{m} generated {values[role][m]['generated']:.5g} control "
                                      f"{values[role][m]['control']:.5g} other {values[role][m]['other_bracket']:.5g}"
                                      for m in METRICS), flush=True)
    passed, checks = verdict(values)
    ended = tree_state(exempt)
    clean = (None if started["clean"] is None or ended["clean"] is None
             else started["clean"] and ended["clean"] and started["commit"] == ended["commit"])
    inputs = {r: {"path": EXPECTED_INPUTS[r], "sha256": harness._sha256(Path(a.data) / EXPECTED_INPUTS[r])}
              for r, _ in BRACKETS}
    inputs["panel"] = {"path": f"board {t.panel_board} genes",
                       "sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest()}
    res = {"label": LABEL, "task": TASK_KEY, "panel_genes": len(genes),
           "settings": {"hold_frac": HOLD_FRAC, "repeats": REPEATS, "hold_seed": HOLD_SEED, "repeat_seed": REPEAT_SEED,
                        "gen_seed": GEN_SEED, "max_variogram_ratio": MAX_VARIOGRAM_RATIO,
                        "max_mmd_share": MAX_MMD_SHARE, "n_latent": t2_latent.N_LATENT,
                        "max_epochs": t2_latent.MAX_EPOCHS, "train_seed": t2_latent.TRAIN_SEED},
           "checks": checks, "repeats": repeats_log,
           "model": dict(model.record) | {"device": device, "threads": int(a.threads)}, "inputs": inputs,
           "verdict": "pass" if passed else "fail",
           "provenance": {"commit": ended["commit"], "tree_clean": clean,
                          "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    paths["json"].write_text(json.dumps(res, indent=1) + "\n")
    text = report_text(res)
    paths["text"].write_text(text)
    for p in paths.values():
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
