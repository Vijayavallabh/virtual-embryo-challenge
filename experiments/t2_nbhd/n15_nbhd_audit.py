#!/usr/bin/env python
"""N15 (experiments/harness/TRIES.md; a diagnostic, not a try): do the embryo entry's two brackets segregate in space,
beyond what coarse cell type explains? The gate of try E15 (fgw_perm). Target-free: no target expression and no score.

For each bracket pair of the embryo board (the two counted backtests' and the test's) and each build seed
1,000,000-1,000,009, the entry's own interpolation (vec.t2_build.interpolate_clouds, called as vec.harness and the s1
build call it, at the board's count; the size step is a uniform scaling, so it is left out and the neighbourhoods are
the same) gives the drawn cells: which bracket each comes from, its coarse type (the lineage map's `group`) and its
position. For every drawn cell, f = the share of after-bracket cells among its 15 nearest cells, itself included (the
scorer's neighbourhoods). Per draw, as ratios only:
  ratio_null          SD(f) / the mean SD(f) under 100 random permutations of the after-bracket label within coarse
                      type, positions fixed (decides: segregation that coarse type does not explain);
  ratio_binomial      SD(f) / sqrt(t (1 - t) / 15), t the pair's bracket fraction (reported);
  null_over_binomial  the null's SD over the binomial SD (reported: what coarse type alone explains);
  after_over_built    SD(f) after the candidate's own re-arrangement (vec.harness.nbhd_rearrangement, E15's fixed radius
                      and stream) over SD(f) as built (reported: what the candidate can remove at its radius).

The pass rule, fixed in the row before this runs: the median of ratio_null over the 10 build seeds exceeds 1.5 on the
test pair and on at least one backtest pair; otherwise the line closes and E15 is not registered. A ratio near 1 closes
it; the other two ratios decide nothing.

No cell count of any stage, of a draw or of a board is printed or written (a count would reveal a test board's limit);
the counts stay in local variables. Reads the released E6.75, E7.25, E7.5 and E8.0 files; E7.25 and E7.5 only as
bracket inputs of the pairs that use them. Their released row counts, which the two backtests' boards need, are the
loader's n_input for those same files (as vec.harness counts rows), not read separately.

Usage (DGX, repo root; about 20 minutes, at most 16 BLAS threads, no GPU; it refuses to overwrite its outputs):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \
    uv run --frozen --no-sync python -I -B experiments/t2_nbhd/n15_nbhd_audit.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import textwrap
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from vec import harness, provenance, t2_build
from vec.adoption import SEEDS
from vec.boards import REPO_ROOT, load_boards
from vec.io import COORD_KEY, dense
from vec.prep import LINEAGE_MAP, load_stage

HERE = Path(__file__).resolve().parent
DATA = REPO_ROOT / "data" / "raw" / "T2"
OUT = HERE / "results" / "n15_nbhd_audit"
EMBRYO_BOARD = "T2:embryo:test_interp"
VAL_BOARD = "T2:embryo:val_interp"
# (key, before, after, fraction): the two counted embryo backtests (the pseudo-tasks' own keys) and the test pair
PAIRS = (("T2e/E6.75+E7.5->E7.25", "E6.75", "E7.5", 2 / 3),
         ("T2e/E7.25+E8.0->E7.5", "E7.25", "E8.0", 1 / 3),
         ("test/E7.5+E8.0->E7.75", "E7.5", "E8.0", 1 / 2))
TEST_PAIR = PAIRS[2][0]
BACKTEST_PAIRS = (PAIRS[0][0], PAIRS[1][0])
K = harness.NBHD_K  # neighbours per cell, self included, as the scorer's
NULL_PERMUTATIONS = 100
NULL_STREAM = 21  # the null's random stream, apart from the interpolation's and the candidate's (harness.NBHD_STREAM)
BAR = 1.5  # the pass rule's bar on the median of ratio_null, fixed in the row before this runs
METRICS = ("ratio_null", "ratio_binomial", "null_over_binomial", "after_over_built")
RULE = (f"Rule (fixed in row N15 before this ran): the median of ratio_null over the 10 build seeds exceeds {BAR} on "
        "the test pair and on at least one backtest pair; otherwise the line closes and E15 is not registered.")


def after_fraction(flags: np.ndarray, neighbours: np.ndarray) -> np.ndarray:
    """f per position: the share of 1s (after-bracket) in `flags` (one value per position) over its neighbours."""
    return np.asarray(flags, dtype=np.float64)[neighbours].mean(1)


def null_sd(flags: np.ndarray, codes: np.ndarray, neighbours: np.ndarray, rng: np.random.Generator,
            permutations: int = NULL_PERMUTATIONS) -> float:
    """The mean of SD(f) over random permutations of `flags` within each value of `codes`, positions fixed."""
    groups = [np.flatnonzero(codes == c) for c in np.unique(codes)]
    sds = []
    for _ in range(permutations):
        shuffled = flags.copy()
        for idx in groups:
            shuffled[idx] = rng.permutation(flags[idx])
        sds.append(after_fraction(shuffled, neighbours).std())
    return float(np.mean(sds))


def _ratio(a: float, b: float) -> float:
    return float(a) / float(b) if b else float("nan")


def audit_draw(coords, after_flags, labels, t: float, seed: int) -> dict[str, float]:
    """The four ratios of one draw (see the module docstring): `coords`, one row per drawn cell; `after_flags`, whether
    it comes from the after bracket; `labels`, its coarse type; `t`, the pair's bracket fraction; `seed`, the build seed
    of the random streams (the candidate's re-arrangement uses the same stream as its build)."""
    C = np.asarray(coords, dtype=np.float64)
    flags = np.asarray(after_flags, dtype=np.float64)
    neighbours = cKDTree(C).query(C, k=K)[1]
    built = float(after_fraction(flags, neighbours).std())
    binomial = float(np.sqrt(t * (1.0 - t) / K))
    codes = np.unique(np.asarray(labels).astype(str), return_inverse=True)[1].reshape(-1)
    null = null_sd(flags, codes, neighbours, np.random.default_rng([seed, NULL_STREAM]))
    at = harness.nbhd_rearrangement(C, labels, seed)
    after = float(after_fraction(flags[np.argsort(at)], neighbours).std())
    return {"ratio_null": _ratio(built, null), "ratio_binomial": _ratio(built, binomial),
            "null_over_binomial": _ratio(null, binomial), "after_over_built": _ratio(after, built)}


def verdict(medians: Mapping[str, float]) -> str:
    """The row's pass rule on the medians of ratio_null per pair (strictly above the bar)."""
    test = medians[TEST_PAIR]
    backs = [k for k in BACKTEST_PAIRS if medians[k] > BAR]
    if test > BAR and backs:
        return f"pass: the median ratio_null exceeds {BAR} on the test pair and on {backs[0]}: register E15 (fgw_perm)"
    failed = []
    if not test > BAR:
        failed.append("the test pair")
    if not backs:
        failed.append("both backtest pairs")
    return f"close: the median ratio_null does not exceed {BAR} on {' or '.join(failed)}: E15 is not registered"


def report_text(res: dict) -> str:
    """The text report: ratios only. `res` is what main writes to the JSON."""
    lines = textwrap.wrap(
        f"N15 (experiments/harness/TRIES.md; a diagnostic, not a try): bracket segregation of the embryo entry's "
        f"mixture, per bracket pair at build seeds {res['build_seeds']}. Target-free: no target and no score; ratios "
        f"only, no cell count of any stage or board. Commit {res['stamp']['commit']}, clean tree: "
        f"{res['stamp']['tree_clean']}, UTC {res['stamp']['utc']}.", 118)
    lines += ["", "ratio_null: SD(f) over its within-coarse-type permutation null (median and range over the draws).",
              "ratio_binomial: SD(f) over sqrt(t(1-t)/15). null_over_binomial: the null's SD over the binomial's.",
              "after_over_built: SD(f) after the candidate's own re-arrangement over SD(f) as built.", ""]
    for key, p in res["pairs"].items():
        lines.append(f"{key} (t = {p['t']:.3f}; board {p['board']}; exceeds the bar: {p['exceeds_bar']}):")
        for m in METRICS:
            v = p["summary"][m]
            lines.append(f"  {m:<20} median {v['median']:.3f} (range {v['min']:.3f} to {v['max']:.3f})")
    lines += ["", *textwrap.wrap(RULE, 118), f"verdict: {res['verdict']}"]
    return "\n".join(lines) + "\n"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _genes_sha256(genes: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(genes).encode()).hexdigest()


def draw_count(key: str, before, boards: Mapping, tasks: Mapping, rows_released: Mapping[str, int]) -> int:
    """The number of cells the candidate's build draws on this pair's board: the harness's count rule on the harness's
    board (on the test pair the test board's own limit, read by code only). Kept in a local variable, never written."""
    if key == TEST_PAIR:
        board = boards[EMBRYO_BOARD]
    else:
        t = tasks[key]
        board = harness.task_board(t, boards, rows_released[harness.stem(t.target)],
                                   rows_released[harness.stem(t.reference)])
    return harness.count_for(harness._candidate_board(board), before.n_obs)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--data", default=str(DATA))
    ap.add_argument("--seeds", type=int, nargs="*", default=list(SEEDS),
                    help="build seeds minus the harness offset; any other than 0-9 gives no verdict")
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N15 has already written its outputs; a rerun needs a new TRIES row")
    if not a.seeds:
        sys.exit("--seeds: give at least one seed")
    commit, unclean = provenance.git_state(provenance.REPO_ROOT)
    data = Path(a.data)
    boards = load_boards()
    harness.check_index(boards)  # no test-board anchors in memory, as the harness's runs refuse
    tasks = {t.key: t for t in harness.load_tasks()[EMBRYO_BOARD]}
    panels = {"val": list(boards[VAL_BOARD].genes), "test": list(boards[EMBRYO_BOARD].genes)}
    stems = {"val": ("E6.75", "E7.25", "E7.5", "E8.0"), "test": ("E7.5", "E8.0")}
    stages = {(panel, s): load_stage(data / f"{s}.h5ad", genes=panels[panel], drop_duplicates=True)
              for panel, ss in stems.items() for s in ss}
    rows_released = {s: int(stages[("val", s)].uns["qc"]["n_input"]) for s in stems["val"]}  # as the harness's rows
    res = {"what": "row N15 of experiments/harness/TRIES.md (a diagnostic, not a try)", "seeds": list(a.seeds),
           "build_seeds": ", ".join(f"{harness.BUILD_SEED_OFFSET + s:,}" for s in a.seeds),
           "constants": {"k": K, "bar": BAR, "null_permutations": NULL_PERMUTATIONS,
                         "radius_factor": harness.NBHD_RADIUS_FACTOR, "sweeps": harness.NBHD_SWEEPS,
                         "labels": "obs['group'] (configs/lineage_map.yaml)"},
           "pairs": {}}
    t0 = time.time()
    for key, sb, sa, t in PAIRS:
        panel = "test" if key == TEST_PAIR else "val"
        genes = panels[panel]
        before, after = stages[(panel, sb)], stages[(panel, sa)]
        n = draw_count(key, before, boards, tasks, rows_released)  # a test-board limit on the test pair: local only
        Cb, Ca = (np.asarray(s.obsm[COORD_KEY], float)[:, :3] for s in (before, after))
        Eb, Ea = (dense(s[:, genes].X).astype(np.float32) for s in (before, after))
        draws = []
        for s in a.seeds:
            seed = harness.BUILD_SEED_OFFSET + s
            res_i = t2_build.interpolate_clouds(Cb, Ca, Eb, Ea, t, n, seed, target_rms=1.0)
            labels = t2_build.drawn_labels(res_i, before.obs["group"], after.obs["group"])
            r = audit_draw(res_i.coords, res_i.stage == 1, labels, t, seed)
            draws.append({"seed": s, "build_seed": seed, **{m: round(r[m], 4) for m in METRICS}})
            print(f"{key} seed {s}: " + " ".join(f"{m} {r[m]:.3f}" for m in METRICS), flush=True)
        board = EMBRYO_BOARD if key == TEST_PAIR else tasks[key].panel_board
        summary = {m: {"median": round(float(np.median([d[m] for d in draws])), 4),
                       "min": round(min(d[m] for d in draws), 4), "max": round(max(d[m] for d in draws), 4)}
                   for m in METRICS}
        res["pairs"][key] = {"before": sb, "after": sa, "t": t, "board": board, "summary": summary,
                             "exceeds_bar": bool(summary["ratio_null"]["median"] > BAR), "draws": draws}
    medians = {k: p["summary"]["ratio_null"]["median"] for k, p in res["pairs"].items()}
    res["verdict"] = verdict(medians) if sorted(a.seeds) == list(SEEDS) else "partial run: no verdict"
    res["inputs"] = {"files_sha256": {f"T2/{s}.h5ad": _sha256(data / f"{s}.h5ad") for s in stems["val"]},
                     "panels_sha256": {p: _genes_sha256(g) for p, g in panels.items()},
                     "pseudo_tasks_sha256": _sha256(harness.CONFIG), "lineage_map_sha256": _sha256(LINEAGE_MAP)}
    res["stamp"] = {"commit": commit, "tree_clean": unclean == [],
                    "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "seconds": round(time.time() - t0, 1)}
    out.mkdir(parents=True, exist_ok=True)
    (out / "n15_nbhd_audit.json").write_text(json.dumps(res, indent=1) + "\n")
    text = report_text(res)
    (out / "n15_nbhd_audit.txt").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
