#!/usr/bin/env python
"""Target-free base check of row R78's two files (the embryo's fourth refinement round on O: T2embryo_s2x_ovt1112.h5ad
and T2embryo_s2x_ovt13.h5ad, each O's command with one option string of build_test.EMBRYO_R78_TABLE), after
experiments/builds/embryo_r75_check.py's half mode. experiments/builds/embryo_r78_line.sh runs it on the test board
after O's identity rebuild and compares the file it writes with O's sha256: on the file's own path (row R78's route
through build_test._r78_refusals and build_embryo) with the file's one change, its t_mix, put back at O's 1/2, the build
must give O's file byte for byte. It prints nothing and exits 0 when the route ran as named below, else 1; a refusal of
the builder, or a command that is not one of row R78's, exits as the builder does or with this script's refusal.

  half ARGS...  build_test.py's own main on ARGS (`embryo ... --captain-ruling R78 --out <file>`, one of row R78's two
                commands) with the file's t_mix put back at O's 1/2, passed to vec.t2_build.interpolate_clouds as 0.5:
                it must then give O's file. Exits 1 unless the command took row R78's route once, interpolate_clouds
                ran once and with t_mix 0.5, and no expression function of vec.t2_build ran (both files are O's
                mixture: every row its bracket's released cell, depth moved only through the shares).

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/embryo_r78_check.py half embryo ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

HALF = Fraction(1, 2)  # O's t_mix: interpolate_clouds then draws as O's call, which leaves t_mix unset (t = 1/2)
# every expression function of vec.t2_build that build_embryo calls for an --expression: none runs for a mixture file
EXPRESSIONS = ("series_rendered_expression", "depth_matched_expression", "depth_rendered_expression",
               "series_back_expression", "blend_expression")
MODES = ("half",)
UNSET = "unset"  # interpolate_clouds called without t_mix, as O's call is
USAGE = "usage: embryo_r78_check.py half embryo ARGS... --out FILE"


def at_half(argv: list[str]) -> int:
    """build_test.main(argv) on row R78's route with the file's t_mix put back at 1/2 (module docstring). Returns 0 when
    the build ran through with every call as named there, else 1."""
    seen = {"files": [], "t_mix": [], "expr": []}
    inner_r78, inner_interp = bt._r78_refusals, bt.t2_build.interpolate_clouds
    inner_expr = {fn: getattr(bt.t2_build, fn) for fn in EXPRESSIONS}

    def r78(a):
        values = inner_r78(a)
        if values is None or values.get("file") not in bt.EMBRYO_R78_FILES:
            sys.exit("embryo_r78_check.py checks row R78's files only: the command must be one of them, with "
                     "--captain-ruling R78")
        seen["files"].append(dict(values))
        return {**values, "t_mix": HALF}

    def interp(*a, **k):
        seen["t_mix"].append(k["t_mix"] if "t_mix" in k else UNSET)
        return inner_interp(*a, **k)

    def expression(fn):
        def run(*a, **k):
            seen["expr"].append(fn)
            return inner_expr[fn](*a, **k)
        return run

    bt._r78_refusals, bt.t2_build.interpolate_clouds = r78, interp
    for fn in EXPRESSIONS:
        setattr(bt.t2_build, fn, expression(fn))
    try:
        bt.main(argv)
    finally:
        bt._r78_refusals, bt.t2_build.interpolate_clouds = inner_r78, inner_interp
        for fn in EXPRESSIONS:
            setattr(bt.t2_build, fn, inner_expr[fn])
    ok = (len(seen["files"]) == 1 and seen["files"][0]["expression"] == "mixture"
          and seen["t_mix"] == [float(HALF)] and seen["expr"] == [])
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] not in MODES or argv[1] != "embryo":
        sys.exit(USAGE)
    return at_half(argv[1:])


if __name__ == "__main__":
    sys.exit(main())
