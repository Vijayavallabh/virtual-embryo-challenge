#!/usr/bin/env python
"""Row N34 of experiments/harness/TRIES.md (T2 embryo; a diagnostic, not a try): the target-free gate of the candidate
fgw_reflect (try E30), the embryo's FGW interpolation with the later bracket brought into the earlier one's frame by
the orthogonal matrix that best fits the coupling, a reflection allowed (vec.t2_build.interpolate_clouds(reflect=True),
vec.core.coupling_rotation(allow_reflection=True), vec.harness.FGW_REFLECT_RULE).

Why: diagnostic N32 found the coupling's partners at the test pair mirrored along one axis in every draw, even from an
unmirrored start and at a lower objective. One reading is that the coupling is right and the two samples' coordinate
frames are mirror images, which no proper rotation aligns; the best orthogonal fit would then be a reflection, the
partners would sit together and the moved cloud would not flatten. This measures whether that happens, on inputs only.

What it reads: the inputs of the three embryo bracket pairs that experiments/harness/diag_emb_coupling.py reads (its
PAIRS: the two counted embryo backtests' brackets and the test's), loaded as that script loads them
(`load_stage(drop_duplicates=True)`), each pair on its own board's panel: the backtests' on T2:embryo:val_interp, the
test pair's on T2:embryo:test_interp (the gene list only, through vec.boards.load_boards). Each stage is read as an
input of the pairs it brackets and of no other: no target stage of any pair is opened, nothing is compared with a
target, no score is computed and no test-board file is built.

For each pair and each build seed 1,000,000 to 1,000,009 it runs the entry's interpolation (fgw) and the candidate's
(fgw_reflect), which share the anchors, the FGW coupling, the draws and every random stream and differ only in the
frame alignment, and reports for each, with diag_emb_coupling's own statistics (`_measure`):
- det: the sign of the frame alignment's determinant (+1 a proper rotation, -1 a reflection);
- xcorr_min: the smallest eigenvalue of the symmetrised whitened cross-correlation of the anchors' positions and their
  partners' (about 1 when partners sit where their anchors do, below 0 when they are mirrored);
- thin_ratio: the moved pool's thinnest principal SD over its RMS radius at the pair's fraction, over the thinner
  bracket's;
- mean_anchor_shift: the mean length of the anchors' moves (the record's two values, before and after anchors);
- label_agreement: the share of coupled mass whose two cells carry the same merged lineage label (untracked aside;
  the same coupling in both, so the same value: shown as a check).

The gate, fixed by the main session before any code (module constants below; row E10's own fix criteria): a draw
passes if fgw_reflect's xcorr_min >= 0.6 AND its thin_ratio >= 0.9. Verdict "open" (the counted try E30 on the two
embryo backtests may run) only if at least 8 of the 10 draws pass on the TEST pair; otherwise "prune". The two backtest
pairs' numbers are reported and decide nothing. A run on fewer seeds gives no verdict.

Outputs, in experiments/harness/results/n34_reflect_gate/: n34_reflect_gate.json and n34_reflect_gate.txt, each with
the commit, the clean-tree flag, the UTC time and the inputs' sha256, and a vec.provenance stamp each. Summaries per
pair and draw only: no cell number of any stage or draw, no coupling count, no per-gene value, no test-board limit.
N34 runs once: it refuses when its output folder already holds anything (a rerun needs a new TRIES row).

Usage (DGX campaign checkout, repo root; after N34's row is committed):
  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMPY_MADVISE_HUGEPAGE=0 \\
      uv run --frozen --no-sync python -I -B experiments/harness/n34_reflect_gate.py \\
      --out experiments/harness/results/n34_reflect_gate
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
OUT = HERE / "results" / "n34_reflect_gate"

ARMS = ("fgw", "fgw_reflect")  # the entry's frame alignment (a proper rotation) and the candidate's
# The gate, fixed by the main session before any code (row E10's own fix criteria).
XCORR_MIN = 0.6  # fgw_reflect's smallest whitened anchor-partner cross-correlation at least this ...
THIN_RATIO_MIN = 0.9  # ... and its moved pool's thinnest axis at least this fraction of the thinner bracket's
DRAWS_MIN = 8  # passing draws (of the 10) needed on the test pair for "open"
BOARD_OF_PAIR = {"T2e/E6.75+E7.5->E7.25": "T2:embryo:val_interp", "T2e/E7.25+E8.0->E7.5": "T2:embryo:val_interp",
                 "test/E7.5+E8.0->E7.75": "T2:embryo:test_interp"}  # each pair on its own board's panel
REPORTED = ("xcorr_min", "xcorr_eigs", "thin_ratio", "label_agreement", "mean_anchor_shift", "rotation_deg",
            "shape_at_t", "nn15_kept")  # what each draw keeps of diag_emb_coupling._measure (no count field)


def _diag():
    """experiments/harness/diag_emb_coupling.py, loaded by path (its statistics, pairs and seeds)."""
    spec = importlib.util.spec_from_file_location("n34_diag_emb_coupling", HERE / "diag_emb_coupling.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


dc = _diag()


def shift(m: dict) -> float:
    """A solution's mean anchor shift: the mean of the record's two values (the before anchors' and the after's)."""
    return float(np.mean(m["mean_anchor_shift"]))


def det_sign(R) -> int:
    """The sign of a frame alignment's determinant: +1 a proper rotation, -1 a reflection."""
    return 1 if np.linalg.det(np.asarray(R, dtype=np.float64)) > 0 else -1


def draw_gate(reflect: dict) -> dict:
    """N34's gate on one draw, from fgw_reflect's measurement."""
    checks = {"xcorr": bool(reflect["xcorr_min"] >= XCORR_MIN),
              "thin": bool(reflect["thin_ratio"] >= THIN_RATIO_MIN)}
    return {"checks": checks, "passes": all(checks.values())}


def verdict(pairs: dict[str, dict], full: bool, test_pair: str) -> str:
    """"open" only if at least DRAWS_MIN draws pass on the test pair of a full run, "prune" otherwise; the backtest
    pairs decide nothing."""
    if not full:
        return "partial run: no verdict"
    return "open" if pairs[test_pair]["n_pass"] >= DRAWS_MIN else "prune"


def _span(v) -> str:
    v = [x for x in v if x is not None]
    if not v:
        return "-"
    lo, hi = f"{min(v):.2f}", f"{max(v):.2f}"
    return lo if lo == hi else f"{lo} to {hi}"


def report_text(res: dict) -> str:
    g = res["gate"]
    lines = textwrap.wrap(
        f"N34 (experiments/harness/TRIES.md; a diagnostic, not a try): the embryo frame alignment with a reflection "
        f"allowed (fgw_reflect) against the entry's proper rotation (fgw), on the inputs of the three embryo bracket "
        f"pairs, build seeds {res['seeds_label']}. Target-free: no target stage opened, no score, no test-board file.",
        118)
    lines += textwrap.wrap(
        f"Gate (fixed before any code): a draw passes if fgw_reflect's xcorr_min >= {g['XCORR_MIN']} and its "
        f"thin_ratio >= {g['THIN_RATIO_MIN']}; open only if at least {g['DRAWS_MIN']} of 10 draws pass on the test "
        f"pair ({g['TEST_PAIR']}); the backtest pairs are reported and decide nothing.", 118)
    for key, p in res["pairs"].items():
        role = "decides" if key == g["TEST_PAIR"] else "reported only"
        lines += ["", f"{key} ({p['genes_board']}'s panel; {role}): {p['n_pass']} of {len(p['draws'])} draws pass",
                  "  seed  fgw: det  xcorr  thin  shift  label | fgw_reflect: det  xcorr  thin  shift  label | pass"]
        for d in p["draws"]:
            f, r = d["fgw"], d["fgw_reflect"]

            def lab(m):
                return "  -  " if m["label_agreement"] is None else f"{m['label_agreement']:.3f}"

            lines.append(f"  {d['seed']:>4}  {f['det']:+d}  {f['xcorr_min']:+.2f}  {f['thin_ratio']:.2f}  "
                         f"{shift(f):.3f}  {lab(f)} | {r['det']:+d}  {r['xcorr_min']:+.2f}  {r['thin_ratio']:.2f}  "
                         f"{shift(r):.3f}  {lab(r)} | {d['gate']['passes']}")
        ds = p["draws"]

        def rng(c, k, ds=ds):
            return _span([d[c][k] for d in ds])

        lines += textwrap.wrap(
            f"  ranges: fgw_reflect xcorr {rng('fgw_reflect', 'xcorr_min')} (fgw {rng('fgw', 'xcorr_min')}); thin "
            f"ratio {rng('fgw_reflect', 'thin_ratio')} (fgw {rng('fgw', 'thin_ratio')}); reflections chosen in "
            f"{sum(d['fgw_reflect']['det'] < 0 for d in ds)} of {len(ds)} draws.", 118, subsequent_indent="  ")
    s = res["stamp"]
    lines += ["", f"Inputs (sha256): {res['inputs']}",
              f"Commit {s['commit']}, tree clean {s['tree_clean']}, UTC {s['utc']}", f"verdict: {res['verdict']}"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="N34: the target-free gate of fgw_reflect")
    ap.add_argument("--data", default=None, help="the released T2 stages (default: data/raw/T2)")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--seeds", type=int, nargs="*", default=list(dc.SEEDS))
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        sys.exit("N34 has already written its outputs; a rerun needs a new TRIES row")
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
            for arm in ARMS:
                res = interpolate_clouds(Cb, Ca, Eb, Ea, t, dc.N_DRAW, seed, target_rms=1.0,
                                         reflect=arm == "fgw_reflect", with_coupling=True)
                if first is None:
                    first = res.coupling["fgw"]
                elif not np.array_equal(first, res.coupling["fgw"]):
                    raise RuntimeError("the FGW coupling is not reproducible within a seed")
                m = dc._measure(res, t, lab_b, lab_a, nn_rows)
                row[arm] = {**{k: m[k] for k in REPORTED}, "det": det_sign(res.coupling["rotation"])}
            row["gate"] = draw_gate(row["fgw_reflect"])
            row["seconds"] = round(time.time() - t0, 1)
            draws.append(row)
            f, r = row["fgw"], row["fgw_reflect"]
            print(f"{key} seed {s}: fgw det {f['det']:+d} xcorr {f['xcorr_min']:+.2f} thin {f['thin_ratio']:.2f} | "
                  f"fgw_reflect det {r['det']:+d} xcorr {r['xcorr_min']:+.2f} thin {r['thin_ratio']:.2f} | pass "
                  f"{row['gate']['passes']} ({row['seconds']} s)", flush=True)
        pairs[key] = {"before": sb_, "after": sa_, "t": t, "genes_board": board,
                      "n_pass": sum(d["gate"]["passes"] for d in draws), "draws": draws}
    full = sorted(a.seeds) == list(dc.SEEDS)
    ended = provenance.git_state(provenance.REPO_ROOT)
    clean = None if started[1] is None or ended[1] is None else (not started[1] and not ended[1]
                                                                 and started[0] == ended[0])
    seeds_label = (f"{harness.BUILD_SEED_OFFSET + dc.SEEDS[0]:,} to {harness.BUILD_SEED_OFFSET + dc.SEEDS[-1]:,}"
                   if full else ", ".join(f"{harness.BUILD_SEED_OFFSET + s:,}" for s in a.seeds))
    res = {"label": "N34 (experiments/harness/TRIES.md): the target-free gate of fgw_reflect; bracket inputs only",
           "verdict": verdict(pairs, full, dc.TEST_PAIR), "seeds": list(a.seeds), "seeds_label": seeds_label,
           "gate": {"XCORR_MIN": XCORR_MIN, "THIN_RATIO_MIN": THIN_RATIO_MIN, "DRAWS_MIN": DRAWS_MIN,
                    "TEST_PAIR": dc.TEST_PAIR, "decides": "the test pair only; the backtest pairs are reported"},
           "pairs": pairs, "inputs": inputs,
           "stamp": {"commit": ended[0], "tree_clean": clean,
                     "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}}
    out.mkdir(parents=True, exist_ok=True)
    text = report_text(res)
    paths = {out / "n34_reflect_gate.json": json.dumps(res, indent=1) + "\n", out / "n34_reflect_gate.txt": text}
    for p, body in paths.items():
        p.write_text(body)
        provenance.write_stamp(p)
    print(text, end="")


if __name__ == "__main__":
    main()
