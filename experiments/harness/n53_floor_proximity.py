#!/usr/bin/env python
"""Row N53 of experiments/harness/TRIES.md (a diagnostic, not a try; target-free): at the test pair, which embryo
floor-form arms build a cloud in the floor's form?

The floor of T2:embryo:test_interp is E7.5 as released, positions included. This read builds, at the test pair (E7.5
and E8.0 as the brackets, t = 1/2) and at the build seeds 1,000,001 to 1,000,010, the positions of three arms with the
harness's own builders (`vec.harness.CANDIDATES[arm].build`, guard and checks included):
  fgw_shape0   the record's floor-form arm (E7.5's rows in place, E8.0's rows at their partners' mean position)
  fgw_shape0x  try E47a's arm (every cell on a real, otherwise unused E7.5 position)
  fgw_aligned  try E10's aligned moves at t = 1/2
and compares each with an independent 10 percent draw of E7.5's own cells (draw A, the floor's form) on the two
shape metrics the first upload lost: occupancy_dice (scorer build "main", as vec.evaluate calls it) and d2_shape
(veckit's d2_distance, the raw value; lower is closer). The first upload's own geometry (fgw) is not an arm of this read
and is never built here.

The noise floor is draw B against draw A at the same seed: B is a uniform draw of E7.5's own cells at the arms' count
(the count rule below; `vec.core.draw_balanced`, stream [seed, COPY_STREAM]), independent of A. It is therefore a
copy of E7.5 at the same count as the arms, it reaches the scorer at the arms' size and goes through the same
subsampling, and after that subsampling it shares about the same share of cells with A (about 10 percent) as
fgw_shape0x does. fgw_shape0 shares about half that (only its E7.5 half sits on E7.5 positions) and fgw_aligned none;
the row states this bias. A second reference value, a 10 percent draw disjoint from A against A (the first version of
this read), is reported beside it and decides nothing.

The panel and the count rule are the embryo validation board's (T2:embryo:val_interp), never the test board's: the
script reads no test board and writes no cell count, pool size or file size. Whether fgw_shape0x's floor-exact
fallback fired under that count rule is written as one boolean. Both metrics ignore the size rule's scale (each
canonicalises or divides by the RMS radius), so the arms' size is the harness's own, read on the two brackets.

Rule, fixed before the read: an arm is "in the floor's form" if, for both metrics, the mean over the ten seeds of the
seed-paired difference (arm against A minus B against A) lies within NOISE_SES standard errors of that mean (its seed
SD over sqrt(10)). "In form" means no difference was detected at that bar, not that the arm equals the floor's form.
Beside it the read reports the wider bar of the first version (the arm's mean within NOISE_SDS single-draw SDs of the
noise floor's mean); that bar decides nothing. The read reports which arms are in form, and decides nothing else.

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
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "experiments" / "t1_hygiene"))
import e47_floor_exact as e47  # noqa: E402
import n12_hygiene_size as n12  # noqa: E402

from vec import harness  # noqa: E402
from vec.boards import load_boards  # noqa: E402
from vec.core import draw_balanced  # noqa: E402
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
REF_FRAC = 0.1  # draw A (and the disjoint reference draw): this share of E7.5's cells
REF_STREAM = 153  # [build seed, 153]: draw A and the disjoint draw, a stream no build step uses
COPY_STREAM = 154  # [build seed, 154]: draw B, the copy of E7.5 at the arms' count
METRIC_SEED = 0  # the scorer's subsampling seed for both metrics
NOISE_SES = 2.0  # the rule's bar: standard errors of the ten-seed mean of the paired difference
NOISE_SDS = 2.0  # the first version's bar, reported only: single-draw SDs of the noise floor
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
    """Draw A and a draw disjoint from it, round(REF_FRAC x n_pool) rows each (stream [seed, REF_STREAM]), sorted."""
    m = int(round(REF_FRAC * n_pool))
    if m < 1 or 2 * m > n_pool:
        raise ValueError("the pool is too small for two disjoint reference draws")
    p = np.random.default_rng([seed, REF_STREAM]).permutation(n_pool)
    return np.sort(p[:m]), np.sort(p[m:2 * m])


def copy_draw(n_pool: int, n: int, seed: int) -> np.ndarray:
    """Draw B: n rows of the pool, uniformly (without replacement when n <= n_pool; `draw_balanced`), on stream
    [seed, COPY_STREAM], independent of draw A."""
    return np.sort(draw_balanced(n_pool, n, np.random.default_rng([seed, COPY_STREAM])))


def _value(v) -> float:
    return float(v[0] if isinstance(v, tuple | list) else v)


def shape_pair(pred: np.ndarray, true: np.ndarray) -> dict[str, float]:
    """occupancy_dice and d2_shape (raw) of a cloud against a reference, under scorer build "main"."""
    sm = _shape_metrics()
    with scorer_build("main"):
        return {"occupancy_dice": _value(sm.occupancy_dice(pred, true, seed=METRIC_SEED)),
                "d2_shape": _value(sm.d2_distance(pred, true, seed=METRIC_SEED))}


def measure(board, inputs: dict, seeds=SEEDS, log=print) -> dict:
    """Per seed, each arm against draw A, draw B (a copy at the arms' count) against A, and the disjoint draw against
    A; and whether fgw_shape0x's fallback fired at any seed (a boolean)."""
    P = np.asarray(inputs["before"].obsm[COORD_KEY], dtype=np.float32)[:, :3]
    n = harness.count_for(harness._candidate_board(board), inputs["before"].n_obs)  # the arms' count; never written
    rows = []
    with e47.floor_exact_records() as seen:
        for s in seeds:
            t0 = time.time()
            a, c = reference_draws(len(P), s)
            ctx = context(board, inputs, s)
            row = {"seed": int(s), "noise": shape_pair(P[copy_draw(len(P), n, s)], P[a]),
                   "noise_disjoint": shape_pair(P[c], P[a])}
            for arm in ARMS:
                row[arm] = shape_pair(positions(ctx, arm), P[a])
            rows.append(row)
            log(f"seed {s}: " + "; ".join(f"{k} " + ", ".join(f"{m} {v[m]:.4f}" for m in METRICS)
                                         for k, v in row.items() if k != "seed") + f" ({time.time() - t0:.1f} s)",
                flush=True)
    return {"rows": rows, "floor_exact_fallback": bool(any(seen.values())), **summary(rows)}


def summary(rows: list[dict]) -> dict:
    """Mean and seed SD per arm and metric, the noise floor's and the disjoint draw's, and the rule."""
    def stats(key):
        return {m: {"mean": float(np.mean([r[key][m] for r in rows])),
                    "sd": float(np.std([r[key][m] for r in rows], ddof=1))} for m in METRICS}

    return {"noise": stats("noise"), "noise_disjoint": stats("noise_disjoint"),
            "arms": {arm: stats(arm) for arm in ARMS}, "in_form": {arm: in_form(rows, arm) for arm in ARMS}}


def in_form(rows: list[dict], arm: str) -> dict:
    """The rule, per metric: the seed-paired difference arm - noise (each against draw A); inside if the absolute
    ten-seed mean is at most NOISE_SES standard errors of that mean. Also the first version's wider bar (the mean
    distance against NOISE_SDS single-draw SDs of the noise floor), reported only."""
    out = {}
    for m in METRICS:
        d = np.array([r[arm][m] - r["noise"][m] for r in rows])
        se = float(d.std(ddof=1) / np.sqrt(len(d)))
        sd_noise = float(np.std([r["noise"][m] for r in rows], ddof=1))
        dist = abs(float(d.mean()))
        out[m] = {"difference": float(d.mean()), "se": se, "bound": float(NOISE_SES * se),
                  "inside": bool(dist <= NOISE_SES * se),
                  "in_noise_sds": float(dist / sd_noise) if sd_noise > 0 else float("inf"),
                  "inside_wide_bar": bool(dist <= NOISE_SDS * sd_noise)}
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
              "copy_stream": COPY_STREAM, "metric_seed": METRIC_SEED, "scorer_build": "main", "noise_ses": NOISE_SES,
              "noise_sds": NOISE_SDS, **res}
    out.mkdir(parents=True, exist_ok=True)
    files[0].write_text(json.dumps(report, indent=1) + "\n")
    files[1].write_text(text(report))
    print(text(report), end="")


def text(report: dict) -> str:
    noise, disjoint = report["noise"], report["noise_disjoint"]
    lines = [f"Row {report['row']} (a diagnostic, not a try; target-free): the floor-form arms at the test pair "
             f"against {report['ref_frac']:.0%} draws of E7.5's own cells, build seeds {report['seeds'][0]:,} to "
             f"{report['seeds'][-1]:,}, panel and count rule of {report['board']}, scorer build main. Commit "
             f"{report['commit']}, clean tree: {report['clean']}. fgw_shape0x's floor-exact fallback under this count "
             f"rule: {report['floor_exact_fallback']}.",
             "  (mean, seed SD)            occupancy_dice (higher closer)   d2_shape (lower closer)",
             "  noise (copy of E7.5 vs A)  " + "   ".join(f"{noise[m]['mean']:.4f} ({noise[m]['sd']:.4f})"
                                                       for m in METRICS),
             "  disjoint draw vs A         " + "   ".join(f"{disjoint[m]['mean']:.4f} ({disjoint[m]['sd']:.4f})"
                                                       for m in METRICS)]
    for arm in report["arms"]:
        s = report["arms"][arm]
        lines.append(f"  {arm:<25}  " + "   ".join(f"{s[m]['mean']:.4f} ({s[m]['sd']:.4f})" for m in METRICS))
    lines.append(f"  rule: the paired difference from the noise floor within {report['noise_ses']:g} SE of its "
                 f"ten-seed mean; beside it, the distance in single-draw noise SDs (the first version's bar, "
                 f"{report['noise_sds']:g} SD, reported only):")
    for arm in report["arms"]:
        r = report["in_form"][arm]
        lines.append(f"  {arm:<25}  in the floor's form: {r['in_form']} (" + "; ".join(
            f"{m} {r[m]['difference']:+.4f} vs {r[m]['bound']:.4f}, {r[m]['in_noise_sds']:.2f} SD" for m in METRICS)
            + ")")
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
