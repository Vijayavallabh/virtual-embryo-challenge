#!/usr/bin/env python
"""Target-free base checks of row R75's five files (the embryo's refinement round on O: T2embryo_s2x_ovt712.h5ad,
T2embryo_s2x_ovt56.h5ad, T2embryo_s2x_ovt512.h5ad, T2embryo_s2xs_ovt23.h5ad and T2embryo_s2x_ovd.h5ad, each O's command
with one option string of build_test.EMBRYO_R75_TABLE). experiments/builds/embryo_r75_line.sh runs it on the test board
after its identity rebuilds and compares the file it writes with a known file's sha256: on the file's own path (row
R75's route through build_test._r75_refusals and build_embryo) with the file's one change put back, the build must give
that known file byte for byte. It prints nothing and exits 0 when the route ran as named below, else 1; a refusal of
the builder, or a command that is not one of row R75's, exits as the builder does or with this script's refusal.

  half ARGS...     build_test.py's own main on ARGS (`embryo ... --captain-ruling R75 --out <file>`, a file with a lean:
                   the three leans and T2embryo_s2xs_ovt23.h5ad) with the file's t_mix put back at O's 1/2, passed to
                   vec.t2_build.interpolate_clouds as 0.5: a lean of 7/12, 5/6 or 5/12 must then give O's file, series
                   at 2/3 SO's. Exits 1 unless the command took row R75's route once, interpolate_clouds ran once with
                   t_mix 0.5, and the file's own expression ran as the file has it: series_rendered_expression once at
                   t = 1/2 with stream [seed, 7] for the series file, no expression function for the mixture files.
  mixture ARGS...  the same with the file's expression put back at O's mixture (a file with --expression:
                   T2embryo_s2xs_ovt23.h5ad, which must then give row R59's file Ot23, and T2embryo_s2x_ovd.h5ad, which
                   must then give O's, coordinates bit for bit): the file's own expression function
                   (series_rendered_expression or depth_matched_expression) runs once on the route, at t = 1/2 with
                   stream [seed, 7], and must leave the interpolation's rows, brackets and coordinates as they were, bit
                   for bit; then every row takes its bracket's released expression, as O's mixture gives it. Exits 1
                   unless the command took row R75's route once, interpolate_clouds ran once with the file's own t_mix
                   (2/3 for the series file, unset for the depth file) and the expression function ran once as named.

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/embryo_r75_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

HALF = Fraction(1, 2)  # O's t_mix: interpolate_clouds then draws as O's call, which leaves t_mix unset (t = 1/2)
# the file's expression -> its function in vec.t2_build, the one build_embryo calls for it
EXPRESSIONS = {"series": "series_rendered_expression", "depth": "depth_matched_expression"}
MODES = ("half", "mixture")
UNSET = "unset"  # interpolate_clouds called without t_mix, as O's call is
USAGE = "usage: embryo_r75_check.py half|mixture embryo ARGS... --out FILE"


def mixture_expression(before, after, interp, genes) -> np.ndarray:
    """O's expression of the interpolation's rows: every row its bracket's released cell, as build_embryo's mixture
    branch gives it (the panel's dense float32 expression of each bracket, read at the row's index)."""
    Eb = bt.dense(before[:, genes].X).astype(np.float32)
    Ea = bt.dense(after[:, genes].X).astype(np.float32)
    X = np.empty((len(interp.rows), len(genes)), dtype=np.float32)
    X[interp.stage == 0], X[interp.stage == 1] = Eb[interp.rows[interp.stage == 0]], Ea[interp.rows[interp.stage == 1]]
    return X


def _bits(interp) -> tuple:
    """The interpolation's rows, brackets and coordinates, each as its dtype, shape and bytes."""
    return tuple((str(np.asarray(x).dtype), np.shape(x), np.ascontiguousarray(x).tobytes())
                 for x in (interp.rows, interp.stage, interp.coords))


def at_base(mode: str, argv: list[str]) -> int:
    """build_test.main(argv) on row R75's route with the file's one change put back (module docstring). Returns 0 when
    the build ran through with every call as named there, else 1."""
    seen = {"files": [], "t_mix": [], "expr": {kind: [] for kind in EXPRESSIONS}, "kept": []}
    inner_r75, inner_interp = bt._r75_refusals, bt.t2_build.interpolate_clouds
    inner_expr = {kind: getattr(bt.t2_build, fn) for kind, fn in EXPRESSIONS.items()}

    def r75(a):
        values = inner_r75(a)
        if values is None or values.get("file") not in bt.EMBRYO_R75_FILES:
            sys.exit("embryo_r75_check.py checks row R75's files only: the command must be one of them, with "
                     "--captain-ruling R75")
        seen["files"].append(dict(values))
        if mode == "half":
            if values["t_mix"] is None:
                sys.exit(f"embryo_r75_check.py half takes a file with a lean, not {values['file']}")
            return {**values, "t_mix": HALF}
        if values["expression"] not in EXPRESSIONS:
            sys.exit(f"embryo_r75_check.py mixture takes a file with --expression, not {values['file']}")
        return values

    def interp(*a, **k):
        seen["t_mix"].append(k["t_mix"] if "t_mix" in k else UNSET)
        return inner_interp(*a, **k)

    def expression(kind):
        def run(before, after, res, genes, t, seed, *a, **k):
            seen["expr"][kind].append((t, list(seed)))
            if mode == "half":
                return inner_expr[kind](before, after, res, genes, t, seed, *a, **k)
            was = _bits(res)
            _, rec = inner_expr[kind](before, after, res, genes, t, seed, *a, **k)
            seen["kept"].append(_bits(res) == was)
            return mixture_expression(before, after, res, genes), rec
        return run

    bt._r75_refusals, bt.t2_build.interpolate_clouds = r75, interp
    for kind, fn in EXPRESSIONS.items():
        setattr(bt.t2_build, fn, expression(kind))
    try:
        bt.main(argv)
    finally:
        bt._r75_refusals, bt.t2_build.interpolate_clouds = inner_r75, inner_interp
        for kind, fn in EXPRESSIONS.items():
            setattr(bt.t2_build, fn, inner_expr[kind])
    if len(seen["files"]) != 1:
        return 1
    values = seen["files"][0]
    own = values["expression"]
    stream = [bt.harness.BUILD_SEED_OFFSET + 1, 7]  # the default seed's expression stream, as build_embryo's
    want_expr = {kind: ([(0.5, stream)] if kind == own else []) for kind in EXPRESSIONS}
    if mode == "half":
        ok = seen["t_mix"] == [float(HALF)] and seen["kept"] == []
    else:
        ok = (seen["t_mix"] == [UNSET if values["t_mix"] is None else float(values["t_mix"])]
              and seen["kept"] == [True])
    return 0 if ok and seen["expr"] == want_expr else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] not in MODES or argv[1] != "embryo":
        sys.exit(USAGE)
    return at_base(argv[0], argv[1:])


if __name__ == "__main__":
    sys.exit(main())
