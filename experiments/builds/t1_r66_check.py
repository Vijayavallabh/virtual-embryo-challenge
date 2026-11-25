#!/usr/bin/env python
"""Target-free check of row R66's thirty files (T1's third-round step-on slots on T1_s2rcMS20's command: for K1, K2
and K3, two or three of line 12's steps together; for C1, g_mom's damping 2.5; for C3, Blood's step at 2.0; the table
is build_test.T1_R66_FILES). experiments/builds/t1_r66_line.sh runs it on the test board after its identity builds.
It prints nothing and exits 0 when it passes, else 1; a refusal of the builder exits as the builder does.

  base ARGS...    build_test.py's own main on ARGS (one of the thirty commands: `t1 ... --out <file>`) on the file's
                  own path (row R66's route through _t1_l12, MS's step, Blood's step and g_mom's hook as the file has
                  them), with every setting the file changes put at T1_s2rcMS20's value: g_mom's damping at 0.5 (the d
                  files), MS's damping at harness.B0_DAMPING (the m files), and Blood's factors at 1, the damping of
                  Blood's momentum_shift call forced to 0 (the H files). Exits 1 unless the command took row R66's
                  route once, MS's step ran once at 2.0 x the calendar fraction and MS's momentum once at B0's damping,
                  g_mom's hook was made once at damping 0.5, and Blood's step ran once at the file's damping in the
                  files with Blood's step and never in the others. The launcher compares the file with T1_s2rcMS20's
                  sha256: with every one of its steps at T1_s2rcMS20's value, each new path must give T1_s2rcMS20's
                  file byte for byte.

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/t1_r66_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

from vec import t1_scope  # noqa: E402

# the record key each of a file's option flags sets in _t1_l12's values
KEYS = {"--g-mom-damping": "g_mom_damping", "--held-momentum": "held_momentum",
        "--direction-momentum-damping": "direction_momentum_damping"}
# file -> the settings it changes from T1_s2rcMS20 (keys of _t1_l12's values), read off its option string
SETTINGS = {name: tuple(KEYS[w] for w in option.split() if w in KEYS) for name, (option, *_) in bt.T1_R66_FILES.items()}
# T1_s2rcMS20's value of each setting but Blood's (T1_s2rcMS20 has no Blood step: the factors are put at 1 instead)
BASE_VALUES = {"g_mom_damping": bt.G_MOM_DAMPINGS[0], "direction_momentum_damping": float(bt.harness.B0_DAMPING)}
BLOOD_WITHIN = {bt.T1_HM_CATEGORY: [bt.T1_HM_LABEL]}  # Blood's momentum_shift call, and no other, restricts to it
USAGE = "usage: t1_r66_check.py base t1 ARGS... --out FILE"


def at_base(argv: list[str]) -> int:
    """build_test.main(argv) on row R66's route with every setting the file changes at T1_s2rcMS20's value (module
    docstring). Returns 0 when the build ran through with every call as named there, else 1."""
    seen = {"files": [], "ms": [], "g_mom": [], "blood": [], "other": []}
    inner_l12, inner_ms, inner_g, inner_shift = bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift

    def l12(a):
        values = inner_l12(a)
        if values is None:
            return None
        seen["files"].append((values["captain_ruling"], values["file"], values["feng_fraction_mult"],
                              values["held_momentum"]))
        if values["captain_ruling"] != bt.T1_R66_ROW or values["file"] not in SETTINGS:
            sys.exit(f"t1_r66_check.py checks row {bt.T1_R66_ROW}'s files only, not row {values['captain_ruling']}'s "
                     f"{values['file']}")
        return {**values, **{key: BASE_VALUES[key] for key in SETTINGS[values["file"]] if key != "held_momentum"}}

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
    _, name, mult, held = seen["files"][0]
    fraction = mult * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    b0 = BASE_VALUES["direction_momentum_damping"]
    ok = (mult == 2.0 and seen["ms"] == [(fraction, b0)] and seen["other"] == [b0]
          and seen["g_mom"] == [BASE_VALUES["g_mom_damping"]]
          and seen["blood"] == ([held] if "held_momentum" in SETTINGS[name] else []))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] != "base" or argv[1] != "t1":
        sys.exit(USAGE)
    return at_base(argv[1:])


if __name__ == "__main__":
    sys.exit(main())
