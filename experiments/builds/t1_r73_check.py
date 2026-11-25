#!/usr/bin/env python
"""Target-free base check of row R73's four files (T1's second refinement round: T1_s2rcMSBm10d10, T1_s2rcMSBm10H05,
T1_s2rcMSBm10d15 and T1_s2rcMSBm10H10, each T1_s2rcMSBm10's command with one further change, g_mom's within-type
damping or Blood's own step; the table is build_test.T1_R73_TABLE). experiments/builds/t1_r73_line.sh runs it on the
test board after its identity rebuilds and compares the file with T1_s2rcMSBm10's sha256: with the one further step put
back at T1_s2rcMSBm10's value, each file's own path must give T1_s2rcMSBm10's file byte for byte. It prints nothing and
exits 0 when it passes, else 1; a refusal of the builder exits as the builder does.

  base ARGS...    build_test.py's own main on ARGS (one of the four commands: `t1 ... --out <file>`) on the file's own
                  path (row R73's route through _t1_l12, MS's step, Blood's step and g_mom's hook as the file has
                  them), with the file's one further step at T1_s2rcMSBm10's value: g_mom's damping at 0.5 (the d
                  files), Blood's factors at 1, the damping of Blood's momentum_shift call forced to 0 (the H files);
                  MS's own damping is kept at the file's 1.0, T1_s2rcMSBm10's. Exits 1 unless the command took row
                  R73's route once, MS's step ran once at 1.5 x the calendar fraction with MS's damping 1.0 and MS's
                  momentum once at 1.0, g_mom's hook was made once at damping 0.5, and Blood's step ran once at the
                  file's damping in the H files and never in the d files.

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/t1_r73_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

from vec import t1_scope  # noqa: E402

# the record key each file's one further flag sets in _t1_l12's values
KEYS = {"--g-mom-damping": "g_mom_damping", "--held-momentum": "held_momentum"}
# file -> the setting its one further step changes from T1_s2rcMSBm10, read off its option string
SETTINGS = {name: tuple(KEYS[w] for w in option.split() if w in KEYS) for name, (option, *_) in bt.T1_R73_FILES.items()}
# T1_s2rcMSBm10's value of each setting but Blood's (it has no Blood step: Blood's factors are put at 1 instead), and
# MS's own damping, which every file of row R73 keeps at T1_s2rcMSBm10's 1.0
BASE_VALUES = {"g_mom_damping": bt.G_MOM_DAMPINGS[0], "direction_momentum_damping": 1.0}
BLOOD_WITHIN = {bt.T1_HM_CATEGORY: [bt.T1_HM_LABEL]}  # Blood's momentum_shift call, and no other, restricts to it
USAGE = "usage: t1_r73_check.py base t1 ARGS... --out FILE"


def at_base(argv: list[str]) -> int:
    """build_test.main(argv) on row R73's route with the file's one further step at T1_s2rcMSBm10's value (module
    docstring). Returns 0 when the build ran through with every call as named there, else 1."""
    seen = {"files": [], "ms": [], "g_mom": [], "blood": [], "other": []}
    inner_l12, inner_ms, inner_g, inner_shift = bt._t1_l12, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift

    def l12(a):
        values = inner_l12(a)
        if values is None:
            return None
        seen["files"].append((values["captain_ruling"], values["file"], values["feng_fraction_mult"],
                              values["held_momentum"], values["direction_momentum_damping"]))
        if values["captain_ruling"] != bt.T1_R73_ROW or values["file"] not in SETTINGS:
            sys.exit(f"t1_r73_check.py checks row {bt.T1_R73_ROW}'s files only, not row {values['captain_ruling']}'s "
                     f"{values['file']}")
        return {**values, **{key: BASE_VALUES[key] for key in SETTINGS[values["file"]] if key != "held_momentum"}}

    def ms_step(ref, library, directions, fraction, groups, held, seed, **k):
        seen["ms"].append((fraction, dict(k)))
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
    _, name, mult, held, ms_damping = seen["files"][0]
    fraction = bt.T1_MSB_MULT * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    m10 = BASE_VALUES["direction_momentum_damping"]
    ok = (mult == bt.T1_MSB_MULT and ms_damping == m10 and seen["ms"] == [(fraction, {"damping": m10})]
          and seen["other"] == [m10] and seen["g_mom"] == [BASE_VALUES["g_mom_damping"]]
          and seen["blood"] == ([held] if "held_momentum" in SETTINGS[name] else []))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] != "base" or argv[1] != "t1":
        sys.exit(USAGE)
    return at_base(argv[1:])


if __name__ == "__main__":
    sys.exit(main())
