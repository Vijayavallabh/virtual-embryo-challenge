#!/usr/bin/env python
"""Target-free check of row R60's seven files (T1's line 12, each one change to T1_s2rcMS20's command: T1_s2rcMS20d10,
T1_s2rcMS20H05, T1_s2rcMS20d20, T1_s2rcMS20m10, T1_s2rcMS20d15, T1_s2rcMS20H10, T1_s2rcMS20m20) and of row R46's pair
T1_s2rcB25MSH05 (row R44's Blood step on T1_s2rcB25MS's command). experiments/builds/t1_l12_line.sh runs it on the
test board after its identity builds. It prints nothing and exits 0 when it passes, else 1; a refusal of the builder
exits as the builder does.

  base ARGS...    build_test.py's own main on ARGS (one of the eight commands: `t1 ... --out <file>`) on the file's own
                  path (row R51's route through _t1_l12, MS's step, Blood's step and g_mom's hook as the file has them),
                  with the one setting the file changes put at its base file's value: g_mom's damping at 0.5
                  (T1_s2rcMS20d10, d15, d20), MS's damping at harness.B0_DAMPING (T1_s2rcMS20m10, m20), or Blood's
                  factors at 1, the damping of Blood's momentum_shift call forced to 0 (T1_s2rcMS20H05, H10 and
                  T1_s2rcB25MSH05). Exits 1 unless the command took row R60's or R46's route once, MS's step ran once at
                  the file's multiplier x the calendar fraction and MS's momentum once at B0's damping, g_mom's hook was
                  made once at damping 0.5, and Blood's step ran once at the file's damping in the three H files and
                  never in the others. The launcher compares the file with its base file's sha256 (T1_s2rcMS20's for
                  row R60's files, T1_s2rcB25MS's for row R46's): with its one setting at the base's value, each new
                  path must give the base's file byte for byte.

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/t1_l12_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

from vec import t1_scope  # noqa: E402


def _setting(values: tuple) -> str:
    """The one setting a file changes from its base (a key of _t1_l12's values), from its T1_R60_FILES entry."""
    _, _, damping, held, ms_damping = values
    if damping is not None:
        return "g_mom_damping"
    return "held_momentum" if held is not None else "direction_momentum_damping"


# file -> the one setting it changes (a key of _t1_l12's values)
SETTING = {**{name: _setting(v) for name, v in bt.T1_R60_FILES.items()}, bt.T1_R46_FILE: "held_momentum"}
# the base files' value of each setting but Blood's (the bases have no Blood step: the factors are put at 1 instead)
BASE_VALUES = {"g_mom_damping": bt.G_MOM_DAMPINGS[0], "direction_momentum_damping": float(bt.harness.B0_DAMPING)}
BLOOD_WITHIN = {bt.T1_HM_CATEGORY: [bt.T1_HM_LABEL]}  # Blood's momentum_shift call, and no other, restricts to it
USAGE = "usage: t1_l12_check.py base t1 ARGS... --out FILE"


def at_base(argv: list[str]) -> int:
    """build_test.main(argv) on row R60's or R46's route with the file's one setting at its base's value (module
    docstring). Returns 0 when the build ran through with every call as named there, else 1."""
    seen = {"files": [], "ms": [], "g_mom": [], "blood": [], "other": []}
    inner_l12, inner_ms, inner_g, inner_shift = bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift

    def l12(a):
        values = inner_l12(a)
        if values is None:
            return None
        seen["files"].append((values["file"], values["feng_fraction_mult"], values["held_momentum"]))
        key = SETTING[values["file"]]
        return values if key == "held_momentum" else {**values, key: BASE_VALUES[key]}

    def ms_step(ref, library, directions, fraction, groups, held, seed, **k):
        seen["ms"].append((fraction, k.get("damping")))
        return inner_ms(ref, library, directions, fraction, groups, held, seed, **k)

    def g_mom(*a, **k):
        seen["g_mom"].append(k.get("damping"))
        return inner_g(*a, **k)

    def shift(start, ref, target_like, damping, *a, **k):
        if k.get("within") == BLOOD_WITHIN:  # Blood's step: its factors at 1
            seen["blood"].append(damping)
            damping = 0.0
        else:
            seen["other"].append(damping)
        return inner_shift(start, ref, target_like, damping, *a, **k)

    bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift = l12, ms_step, g_mom, shift
    try:
        bt.main(argv)
    finally:
        bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift = inner_l12, inner_ms, inner_g, inner_shift
    if len(seen["files"]) != 1:
        return 1
    name, mult, held = seen["files"][0]
    fraction = mult * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    b0 = BASE_VALUES["direction_momentum_damping"]
    ok = (seen["ms"] == [(fraction, b0)] and seen["other"] == [b0] and seen["g_mom"] == [BASE_VALUES["g_mom_damping"]]
          and seen["blood"] == ([held] if SETTING[name] == "held_momentum" else []))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] != "base" or argv[1] != "t1":
        sys.exit(USAGE)
    return at_base(argv[1:])


if __name__ == "__main__":
    sys.exit(main())
