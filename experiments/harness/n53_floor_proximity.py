#!/usr/bin/env python
"""Row N53 of experiments/harness/TRIES.md (a diagnostic, not a try; target-free): at the test pair, which embryo
floor-form arms build a cloud in the floor's form?

The floor of T2:embryo:test_interp is E7.5 as released, positions included. This read builds, at the test pair (E7.5
and E8.0 as the brackets, t = 1/2) and at the build seeds 1,000,001 to 1,000,010, the positions of three arms with the
harness's own builders (`vec.harness.CANDIDATES[arm].build`, guard and checks included):
  fgw_shape0   the record's floor-form arm (E7.5's rows in place, E8.0's rows at their partners' mean position)
  fgw_shape0x  try E47a's arm (every cell on a real, otherwise unused E7.5 position)
  fgw_aligned  try E10's aligned moves at t = 1/2
and compares each with an independent 10 percent draw of E7.5's own cells (the floor's form) on the two shape metrics
the first upload lost: occupancy_dice (scorer build "main", as vec.evaluate calls it) and d2_shape (veckit's
d2_distance, the raw value; lower is closer). The noise floor is two disjoint 10 percent draws of E7.5 against each
other at the same seeds. The first upload's own geometry (fgw) is not an arm of this read and is never built here.

The panel and the count rule are the embryo validation board's (T2:embryo:val_interp), never the test board's: the
script reads no test board and writes no cell count, pool size or file size. Both metrics ignore the size rule's scale
(each canonicalises or divides by the RMS radius), so the arms' size is the harness's own, read on the two brackets.

Rule, fixed before the read: an arm is "in the floor's form" if both its means (over the ten seeds) lie within
NOISE_SDS noise standard deviations (the seed SD of the E7.5-against-E7.5 value) of the noise floor's means. The read
reports which arms are, and decides nothing else.

Refused: a run without row N53 in TRIES.md, a dirty tree, an existing output.

Usage (DGX campaign checkout, repo root):
  OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 \
    uv run --frozen --no-sync python -I -B experiments/harness/n53_floor_proximity.py run \
    --out experiments/harness/results/n53_floor_proximity
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))
import n12_hygiene_size as n12  # noqa: E402

from vec import harness  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.evaluate import _shape_metrics, scorer_build  # noqa: E402
from vec.io import COORD_KEY  # noqa: E402
from vec.prep import load_stage  # noqa: E402
from vec.provenance import repo_state  # noqa: E402

ROW = "N53"
TRIES = HERE / "TRIES.md"
OUT = HERE / "results" / "n53_floor_proximity"
BOARD = "T2:embryo:val_interp"  # the panel and the count rule; never the test board
FILES = {"before": ("T2/E7.5.h5ad", 7.5), "after": ("T2/E8.0.h5ad", 8.0)}
TARGET_DAY = 7.75
SEEDS = tuple(range(harness.BUILD_SEED_OFFSET + 1, harness.BUILD_SEED_OFFSET + 11))  # 1,000,001 to 1,000,010
ARMS = ("fgw_shape0", "fgw_shape0x", "fgw_aligned")
NOT_AN_ARM = "fgw"  # the first upload's moved geometry: never built in this read
REF_FRAC = 0.1  # each reference draw: this share of E7.5's cells
REF_STREAM = 153  # [build seed, 153]: the two disjoint reference draws, a stream no build step uses
METRIC_SEED = 0  # the scorer's subsampling seed for both metrics
NOISE_SDS = 2.0
METRICS = ("occupancy_dice", "d2_shape")
assert NOT_AN_ARM not in ARMS


def view() -> harness.TaskView:
    """The test pair as the harness shows it to a builder (build_test.py's view), named by the validation board."""
    return harness.TaskView(key="N53:E7.5+E8.0->E7.75", task="T2", kind="interp", setting="embryo", test_board=BOARD,
                            roles=tuple(FILES), role_days={r: d for r, (_, d) in FILES.items()},
                            role_stems={r: harness.stem(f) for r, (f, _) in FILES.items()}, reference_role="before",
                            target_day=TARGET_DAY)


def context(board, inputs: dict, seed: int) -> harness.BuildContext:
    """A builder's context at the test pair with the validation board's panel and limits (as a candidate sees them)."""
    return harness.BuildContext(task=view(), board=harness._candidate_board(board), inputs=inputs, seed=seed)


def positions(ctx: harness.BuildContext, arm: str) -> np.ndarray:
    """One arm's positions as its file holds them (float32), built by the harness's own candidate."""
    if arm not in ARMS:
        raise ValueError(f"{arm!r} is not an arm of row N53 (arms: {ARMS}; the first upload's geometry is never built)")
    return np.asarray(harness.CANDIDATES[arm].build(ctx).coords, dtype=np.float32)


def reference_draws(n_pool: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Two disjoint uniform draws of round(REF_FRAC x n_pool) rows each (stream [seed, REF_STREAM]), rows sorted."""
    m = int(round(REF_FRAC * n_pool))
    if m < 1 or 2 * m > n_pool:
        raise ValueError("the pool is too small for two disjoint reference draws")
    p = np.random.default_rng([seed, REF_STREAM]).permutation(n_pool)
    return np.sort(p[:m]), np.sort(p[m:2 * m])


def _value(v) -> float:
    return float(v[0] if isinstance(v, tuple | list) else v)


def shape_pair(pred: np.ndarray, true: np.ndarray) -> dict[str, float]:
    """occupancy_dice and d2_shape (raw) of a cloud against a reference, under scorer build "main"."""
    sm = _shape_metrics()
    with scorer_build("main"):
        return {"occupancy_dice": _value(sm.occupancy_dice(pred, true, seed=METRIC_SEED)),
                "d2_shape": _value(sm.d2_distance(pred, true, seed=METRIC_SEED))}


def measure(board, inputs: dict, seeds=SEEDS, log=print) -> dict:
    """Per seed, each arm against reference draw A and draw A against draw B."""
    P = np.asarray(inputs["before"].obsm[COORD_KEY], dtype=np.float32)[:, :3]
    rows = []
    for s in seeds:
        t0 = time.time()
        a, b = reference_draws(len(P), s)
        ctx = context(board, inputs, s)
        row = {"seed": int(s), "noise": shape_pair(P[b], P[a])}
        for arm in ARMS:
            row[arm] = shape_pair(positions(ctx, arm), P[a])
        rows.append(row)
        log(f"seed {s}: " + "; ".join(f"{k} " + ", ".join(f"{m} {v[m]:.4f}" for m in METRICS)
                                     for k, v in row.items() if k != "seed") + f" ({time.time() - t0:.1f} s)",
            flush=True)
    return {"rows": rows, **summary(rows)}


def summary(rows: list[dict]) -> dict:
    """Mean and seed SD per arm and metric, the noise floor's, and the rule."""
    def stats(key):
        return {m: {"mean": float(np.mean([r[key][m] for r in rows])),
                    "sd": float(np.std([r[key][m] for r in rows], ddof=1))} for m in METRICS}

    noise = stats("noise")
    arms = {arm: stats(arm) for arm in ARMS}
    return {"noise": noise, "arms": arms, "in_form": {arm: in_form(arms[arm], noise) for arm in ARMS}}


def in_form(arm: dict, noise: dict) -> dict:
    """The rule: both metrics' means within NOISE_SDS noise SDs of the noise floor's means."""
    out = {}
    for m in METRICS:
        dist = abs(arm[m]["mean"] - noise[m]["mean"])
        out[m] = {"distance": float(dist), "bound": float(NOISE_SDS * noise[m]["sd"]),
                  "inside": bool(dist <= NOISE_SDS * noise[m]["sd"])}
    out["in_form"] = all(out[m]["inside"] for m in METRICS)
    return out


def _state() -> tuple[str, bool, str]:
    s = repo_state(REPO)
    return s.commit, s.unclean == (), s.utc


def cmd_run(a) -> None:
    if ROW not in n12.tries_rows(TRIES):
        sys.exit(f"{ROW} has no row in experiments/harness/TRIES.md: a row is committed before its run")
    out = Path(a.out)
    files = [out / "n53.json", out / "n53.txt"]
    existing = [str(f) for f in files if f.exists()]
    if existing:
        sys.exit(f"refusing to overwrite: {', '.join(existing)}")
    commit, clean, utc = _state()
    if not clean:
        sys.exit(f"the tree is not clean: {ROW} runs on a committed tree")
    board = load_boards()[BOARD]
    genes = list(board.genes)
    inputs = {r: load_stage(Path(a.data) / f, genes=genes, drop_duplicates=True) for r, (f, _) in FILES.items()}
    res = measure(board, inputs, seeds=SEEDS)
    _, clean_end, utc_end = _state()
    report = {"what": __doc__.split("\n")[0], "row": ROW, "commit": commit, "clean": bool(clean and clean_end),
              "utc": [utc, utc_end], "board": BOARD, "brackets": [f for f, _ in FILES.values()], "t": 0.5,
              "seeds": list(SEEDS), "arms": list(ARMS), "ref_frac": REF_FRAC, "ref_stream": REF_STREAM,
              "metric_seed": METRIC_SEED, "scorer_build": "main", "noise_sds": NOISE_SDS, **res}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(report, indent=1) + "\n")
    files[1].write_text(text(report))
    print(text(report), end="")


def text(report: dict) -> str:
    noise = report["noise"]
    lines = [f"Row {report['row']} (a diagnostic, not a try; target-free): the floor-form arms at the test pair "
             f"against {report['ref_frac']:.0%} draws of E7.5's own cells, build seeds {report['seeds'][0]:,} to "
             f"{report['seeds'][-1]:,}, panel and count rule of {report['board']}, scorer build main. Commit "
             f"{report['commit']}, clean tree: {report['clean']}.",
             "  (mean, seed SD)       occupancy_dice (higher closer)   d2_shape (lower closer)   in the floor's form",
             "  noise (E7.5 vs E7.5)  " + "   ".join(f"{noise[m]['mean']:.4f} ({noise[m]['sd']:.4f})" for m in METRICS)]
    for arm in report["arms"]:
        s, r = report["arms"][arm], report["in_form"][arm]
        lines.append(f"  {arm:<20}  " + "   ".join(f"{s[m]['mean']:.4f} ({s[m]['sd']:.4f})" for m in METRICS)
                     + f"   {r['in_form']} (" + ", ".join(f"{m} {r[m]['distance']:.4f} vs {r[m]['bound']:.4f}"
                                                         for m in METRICS) + ")")
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", default=str(OUT))
    r.add_argument("--data", default=str(harness.DATA_DIR))
    a = ap.parse_args(argv)
    cmd_run(a)


if __name__ == "__main__":
    main()
