#!/usr/bin/env python
"""Target-free base check of row R67's three files (T1_s2rcMSBm125, T1_s2rcMSBm15 and T1_s2rcMSBm075, each
T1_s2rcMSB's command with MS's own damping at 1.25, 1.5 or 0.75; the table is build_test.T1_R67_TABLE).
experiments/builds/t1_r67_line.sh runs it on the test board after its identity rebuilds and compares the file with
T1_s2rcMSB's sha256: with MS's damping put back at T1_s2rcMSB's value, each file's own path must give MSB's file byte
for byte. It prints nothing and exits 0 when it passes, else 1; a refusal of the builder exits as the builder does.

  base ARGS...    build_test.py's own main on ARGS (one of the three commands: `t1 ... --out <file>`) on the file's own
                  path (row R67's route, MS's step with the file's own note, g_mom's hook as the file has them), with
                  MS's damping put back at harness.B0_DAMPING (T1_s2rcMSB's). Exits 1 unless the command took row
                  R67's route once, MS's step ran once at 1.5 x the calendar fraction with MS's damping at B0's (and
                  the file's own note, which a damping at B0's leaves out of the record), MS's momentum ran once at
                  B0's damping, and g_mom's hook was made once at damping 0.5.

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/t1_r67_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

from vec import t1_scope  # noqa: E402

USAGE = "usage: t1_r67_check.py base t1 ARGS... --out FILE"


def at_base(argv: list[str]) -> int:
    """build_test.main(argv) on row R67's route with MS's damping at B0's (module docstring). Returns 0 when the build
    ran through with every call as named there, else 1."""
    seen = {"files": [], "ms": [], "g_mom": [], "shift": []}
    inner_r67, inner_ms, inner_g, inner_shift = bt._t1_r67, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift
    b0 = float(bt.harness.B0_DAMPING)

    def r67(a):
        values = inner_r67(a)
        if values is None:
            return None
        seen["files"].append(values["file"])
        return {**values, "direction_momentum_damping": b0}

    def ms_step(ref, library, directions, fraction, groups, held, seed, **k):
        seen["ms"].append((fraction, dict(k)))
        return inner_ms(ref, library, directions, fraction, groups, held, seed, **k)

    def g_mom(*a, **k):
        seen["g_mom"].append(k.get("damping"))
        return inner_g(*a, **k)

    def shift(start, ref, target_like, damping, *a, **k):
        seen["shift"].append(damping)
        return inner_shift(start, ref, target_like, damping, *a, **k)

    bt._t1_r67, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift = r67, ms_step, g_mom, shift
    try:
        bt.main(argv)
    finally:
        bt._t1_r67, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift = inner_r67, inner_ms, inner_g, inner_shift
    if len(seen["files"]) != 1 or seen["files"][0] not in bt.T1_R67_FILES:
        return 1
    fraction = bt.T1_MSB_MULT * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    want_kw = {"damping": b0, "option_note": bt.T1_R67_NOTES[seen["files"][0]]}
    ok = (seen["ms"] == [(fraction, want_kw)] and seen["shift"] == [b0]
          and seen["g_mom"] == [bt.G_MOM_DAMPINGS[0]])
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] != "base" or argv[1] != "t1":
        sys.exit(USAGE)
    return at_base(argv[1:])


if __name__ == "__main__":
    sys.exit(main())
