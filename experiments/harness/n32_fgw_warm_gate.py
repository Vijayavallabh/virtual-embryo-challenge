#!/usr/bin/env python
"""Row N32 of experiments/harness/TRIES.md (T2 embryo; a diagnostic, not a try): the target-free gate of the candidate
fgw_warm, the embryo's FGW coupling solved a second time from the aligned coupling of try E10 as its start
(vec.t2_build.interpolate_clouds(move_coupling="fgw_warm"), vec.harness.FGW_WARM_RULE).

What it reads: the inputs of the three embryo bracket pairs that experiments/harness/diag_emb_coupling.py reads (its
PAIRS: the two counted embryo backtests' brackets and the test's), loaded as that script loads them
(`load_stage(drop_duplicates=True)`), each pair on its own board's panel: the backtests' on T2:embryo:val_interp,
the test pair's on T2:embryo:test_interp (the gene list only, through vec.boards.load_boards). Each stage is read as an
input of the pairs it brackets and of no other: no target stage of any pair is opened, nothing is compared with a
target, no score is computed and no test-board file is built.

For each pair and each build seed 1,000,000 to 1,000,009 it runs the entry's interpolation (move_coupling="fgw") and
the candidate's ("fgw_warm"), which share the anchors, the first FGW solution, the rotation, the draws and every random
stream, and reports for each, with diag_emb_coupling's own statistics (`_measure`):
- xcorr_min: the smallest eigenvalue of the symmetrised whitened cross-correlation of the anchors' positions and their
  partners' (about 1 when partners sit where their anchors do, below 0 when they are mirrored);
- thin_ratio: the moved pool's thinnest principal SD over its RMS radius at the pair's fraction, over the thinner
  bracket's;
- mean_anchor_shift: the mean length of the anchors' moves (the record's two values, before and after anchors, and
  their mean, which the gate reads);
- label_agreement: the share of coupled mass whose two cells carry the same merged lineage label (untracked aside);
- the FGW objective (vec.core.fgw_objective, POT's fgw_dist) of the first solution, of the aligned start and of the
  warm solution.

The gate, fixed in N32's row before the run (module constants below): a draw passes if fgw_warm's xcorr_min >= 0.6 AND
its thin_ratio >= 0.9 AND its mean anchor shift >= 0.9 x fgw's in the same draw. Verdict "open" (a counted try of
fgw_warm on the two embryo backtests may be registered) only if at least 8 of the 10 draws pass on EACH of the three
pairs; otherwise "prune". A run on fewer seeds gives no verdict.

Outputs, in experiments/harness/results/n32_fgw_warm_gate/: n32_fgw_warm_gate.json and n32_fgw_warm_gate.txt, each
with the commit, the clean-tree flag, the UTC time and the inputs' sha256, and a vec.provenance stamp each. Summaries
per pair and draw only: no cell number of any stage or draw, no coupling count, no per-gene value, no test-board limit.
N32 runs once: it refuses when its output folder already holds anything (a rerun needs a new TRIES row).

Usage (DGX campaign checkout, repo root; after N32's row is committed):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/harness/n32_fgw_warm_gate.py \\
      --out experiments/harness/results/n32_fgw_warm_gate
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import textwrap
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "n32_fgw_warm_gate"

COUPLINGS = ("fgw", "fgw_warm")  # the entry's moves and the candidate's
# The gate, fixed in N32's row (experiments/harness/TRIES.md) before the run.
XCORR_MIN = 0.6  # fgw_warm's smallest whitened anchor-partner cross-correlation at least this ...
THIN_RATIO_MIN = 0.9  # ... and its moved pool's thinnest axis at least this fraction of the thinner bracket's ...
SHIFT_RATIO_MIN = 0.9  # ... and its mean anchor shift at least this fraction of fgw's in the same draw
DRAWS_MIN = 8  # passing draws (of the 10) needed on each of the three pairs for "open"
BOARD_OF_PAIR = {"T2e/E6.75+E7.5->E7.25": "T2:embryo:val_interp", "T2e/E7.25+E8.0->E7.5": "T2:embryo:val_interp",
                 "test/E7.5+E8.0->E7.75": "T2:embryo:test_interp"}  # each pair on its own board's panel
REPORTED = ("xcorr_min", "xcorr_eigs", "thin_ratio", "label_agreement", "mean_anchor_shift", "rotation_deg",
            "shape_at_t", "nn15_kept")  # what each draw keeps of diag_emb_coupling._measure (no count field)


def _diag():
    """experiments/harness/diag_emb_coupling.py, loaded by path (its statistics, pairs and seeds)."""
    spec = importlib.util.spec_from_file_location("n32_diag_emb_coupling", HERE / "diag_emb_coupling.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


dc = _diag()


def shift(m: dict) -> float:
    """A solution's mean anchor shift: the mean of the record's two values (the before anchors' and the after's)."""
    return float(np.mean(m["mean_anchor_shift"]))


def draw_gate(fgw: dict, warm: dict) -> dict:
    """N32's gate on one draw, from the two measurements."""
    ratio = shift(warm) / shift(fgw) if shift(fgw) > 0 else float("inf")
    checks = {"xcorr": bool(warm["xcorr_min"] >= XCORR_MIN), "thin": bool(warm["thin_ratio"] >= THIN_RATIO_MIN),
              "shift": bool(ratio >= SHIFT_RATIO_MIN)}
    return {"shift_ratio": ratio, "checks": checks, "passes": all(checks.values())}


def verdict(pairs: dict[str, dict], full: bool) -> str:
    """"open" only if at least DRAWS_MIN draws pass on each pair of a full run, "prune" otherwise."""
    if not full:
        return "partial run: no verdict"
    return "open" if all(p["n_pass"] >= DRAWS_MIN for p in pairs.values()) else "prune"


def _span(v) -> str:
    v = [x for x in v if x is not None]
    if not v:
        return "-"
    lo, hi = f"{min(v):.2f}", f"{max(v):.2f}"
    return lo if lo == hi else f"{lo} to {hi}"


def report_text(res: dict) -> str:
    g = res["gate"]
    lines = textwrap.wrap(
        f"N32 (experiments/harness/TRIES.md; a diagnostic, not a try): the embryo coupling started from the aligned "
        f"coupling (fgw_warm) against the entry's (fgw), on the inputs of the three embryo bracket pairs, build seeds "
        f"{res['seeds_label']}. Target-free: no target stage opened, no score, no test-board file.", 118)
    lines += textwrap.wrap(
        f"Gate (fixed before the run): a draw passes if fgw_warm's xcorr_min >= {g['XCORR_MIN']}, its thin_ratio >= "
        f"{g['THIN_RATIO_MIN']} and its mean anchor shift >= {g['SHIFT_RATIO_MIN']} x fgw's; open only if at least "
        f"{g['DRAWS_MIN']} of 10 draws pass on each pair.", 118)
    for key, p in res["pairs"].items():
        lines += ["", f"{key} ({p['genes_board']}'s panel): {p['n_pass']} of {len(p['draws'])} draws pass",
                  "  seed  fgw: xcorr  thin  shift  label | fgw_warm: xcorr  thin  shift  label | shift ratio  "
                  "objective fgw / start / warm  pass"]
        for d in p["draws"]:
            f, w, o = d["fgw"], d["fgw_warm"], d["objective"]

            def lab(m):
                return "  -  " if m["label_agreement"] is None else f"{m['label_agreement']:.3f}"

            lines.append(f"  {d['seed']:>4}  {f['xcorr_min']:+.2f}  {f['thin_ratio']:.2f}  {shift(f):.3f}  {lab(f)} | "
                         f"{w['xcorr_min']:+.2f}  {w['thin_ratio']:.2f}  {shift(w):.3f}  {lab(w)} | "
                         f"{d['gate']['shift_ratio']:.2f}  {o['fgw']:.5f} / {o['aligned_start']:.5f} / "
                         f"{o['fgw_warm']:.5f}  {d['gate']['passes']}")
        ds = p["draws"]
        def rng(c, k, ds=ds):
            return _span([d[c][k] for d in ds])

        lines += textwrap.wrap(
            f"  ranges: fgw_warm xcorr {rng('fgw_warm', 'xcorr_min')} (fgw {rng('fgw', 'xcorr_min')}); thin ratio "
            f"{rng('fgw_warm', 'thin_ratio')} (fgw {rng('fgw', 'thin_ratio')}); shift ratio "
            f"{rng('gate', 'shift_ratio')}.", 118, subsequent_indent="  ")
    s = res["stamp"]
    lines += ["", f"Inputs (sha256): {res['inputs']}",
              f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}", f"verdict: {res['verdict']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="N32: the target-free gate of fgw_warm")
    ap.add_argument("--data", default=None, help="the released T2 stages (default: data/raw/T2)")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--seeds", type=int, nargs="*", default=list(dc.SEEDS))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N32 has already written its outputs; a rerun needs a new TRIES row")
    if not a.seeds:
        sys.exit("--seeds: give at least one seed")

    from vec import harness, provenance
    from vec.t2_build import interpolate_clouds

    started = provenance.git_state(provenance.REPO_ROOT)
    data = Path(a.data) if a.data else dc.DATA
    boards = dc.load_boards()
    pairs, inputs = {}, {}
    for key, sb_, sa_, t in dc.PAIRS:
        board = BOARD_OF_PAIR[key]
        genes = list(boards[board].genes)
        before, after = (dc.load_stage(data / f"{s}.h5ad", genes=genes, drop_duplicates=True) for s in (sb_, sa_))
        for s in (sb_, sa_):
            if f"{s}.h5ad" not in inputs:
                inputs[f"{s}.h5ad"] = harness._sha256(data / f"{s}.h5ad")
        Cb, Ca, Eb, Ea, lab_b, lab_a = dc._inputs(before, after, genes)
        draws = []
        for s in a.seeds:
            seed, t0 = harness.BUILD_SEED_OFFSET + s, time.time()
            nn_rows = dc._nn_rows(seed, len(Cb), len(Ca))
            row, first = {"seed": s, "build_seed": seed}, None
            for mc in COUPLINGS:
                res = interpolate_clouds(Cb, Ca, Eb, Ea, t, dc.N_DRAW, seed, target_rms=1.0, move_coupling=mc,
                                         with_coupling=True)
                if first is None:
                    first = res.coupling["fgw"]
                elif not np.array_equal(first, res.coupling["fgw"]):
                    raise RuntimeError("the first FGW solution is not reproducible within a seed")
                m = dc._measure(res, t, lab_b, lab_a, nn_rows)
                row[mc] = {k: m[k] for k in REPORTED}
                if mc == "fgw_warm":
                    row["objective"] = dict(res.record["fgw_objective"])
            row["gate"] = draw_gate(row["fgw"], row["fgw_warm"])
            row["seconds"] = round(time.time() - t0, 1)
            draws.append(row)
            f, w = row["fgw"], row["fgw_warm"]
            print(f"{key} seed {s}: fgw xcorr {f['xcorr_min']:+.2f} thin {f['thin_ratio']:.2f} shift {shift(f):.3f} | "
                  f"fgw_warm xcorr {w['xcorr_min']:+.2f} thin {w['thin_ratio']:.2f} shift {shift(w):.3f} | pass "
                  f"{row['gate']['passes']} ({row['seconds']} s)", flush=True)
        pairs[key] = {"before": sb_, "after": sa_, "t": t, "genes_board": board,
                      "n_pass": sum(d["gate"]["passes"] for d in draws), "draws": draws}
    full = sorted(a.seeds) == list(dc.SEEDS)
    ended = provenance.git_state(provenance.REPO_ROOT)
    clean = None if started[1] is None or ended[1] is None else (not started[1] and not ended[1]
                                                                 and started[0] == ended[0])
    seeds_label = (f"{harness.BUILD_SEED_OFFSET + dc.SEEDS[0]:,} to {harness.BUILD_SEED_OFFSET + dc.SEEDS[-1]:,}"
                   if full else ", ".join(f"{harness.BUILD_SEED_OFFSET + s:,}" for s in a.seeds))
    res = {"label": "N32 (experiments/harness/TRIES.md): the target-free gate of fgw_warm; bracket inputs only",
           "verdict": verdict(pairs, full), "seeds": list(a.seeds), "seeds_label": seeds_label,
           "gate": {"XCORR_MIN": XCORR_MIN, "THIN_RATIO_MIN": THIN_RATIO_MIN, "SHIFT_RATIO_MIN": SHIFT_RATIO_MIN,
                    "DRAWS_MIN": DRAWS_MIN, "shift": "the mean of the record's two mean anchor shifts"},
           "pairs": pairs, "inputs": inputs,
           "stamp": {"commit": ended[0], "tree_clean": clean,
                     "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    text = report_text(res)
    paths = {out / "n32_fgw_warm_gate.json": json.dumps(res, indent=1) + "\n", out / "n32_fgw_warm_gate.txt": text}
    for p, body in paths.items():
        p.write_text(body)
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
