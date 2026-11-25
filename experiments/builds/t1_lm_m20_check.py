#!/usr/bin/env python
"""Target-free identity builds of rows R54's and R58's files (T1_s2rcMSB_LM16 and T1_s2rcMSBm20, each one change to
T1_s2rcMSB's command). experiments/builds/t1_lm_m20_line.sh runs them on the test board after its identity rebuilds
and compares each file with T1_s2rcMSB's sha256: with its one change put back at T1_s2rcMSB's value, each new path must
give MSB's file byte for byte. Each mode exits 0 when the build ran through the file's own route with every call as
named below, else 1; a refusal of the builder exits as the builder does. It prints nothing.

  m20msb ARGS...   build_test.py's own main on row R58's command (`t1 ... --out <file>`) on the file's own path (row
                   R58's route, MS's step, g_mom's hook as the file has them), with MS's damping put back at
                   harness.B0_DAMPING (T1_s2rcMSB's). Exits 1 unless the command took row R58's route once, MS's step
                   ran once at 1.5 x the calendar fraction with MS's damping at B0's, and g_mom's hook was made once at
                   damping 0.5.
  lm16own ARGS...  build_test.py's own main on row R54's command on the file's own path (row R54's route, its shares
                   file and guards, MS's step, g_mom's hook), with the moved shares put back at the pool's own chamber
                   mix (the late-mix step returns the mix it is given, so build_submission draws at T1_s2rcMSB's
                   shares). Exits 1 unless the command took row R54's route once, the late-mix step ran once and was
                   given the pool's own mix, MS's step ran once at 1.5 x the calendar fraction with no damping of its
                   own, and g_mom's hook was made once at damping 0.5. (Identity at the moved shares' own values is a
                   quota question, not a byte one: with_share rescales floats; the synthetic tests hold it on quotas.)

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/t1_lm_m20_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

USAGE = "usage: t1_lm_m20_check.py m20msb|lm16own t1 ARGS... --out FILE"


def _ms_calls(seen: dict):
    inner_ms, inner_g = bt._t1_ms_step, bt._t1_g_mom

    def ms_step(ref, library, directions, fraction, groups, held, seed, **k):
        seen["ms"].append((fraction, dict(k)))
        return inner_ms(ref, library, directions, fraction, groups, held, seed, **k)

    def g_mom(*a, **k):
        seen["g_mom"].append(k.get("damping"))
        return inner_g(*a, **k)

    return inner_ms, inner_g, ms_step, g_mom


def m20msb(argv: list[str]) -> int:
    """Row R58's file with MS's damping at B0's (module docstring); 0 when every call was as named, else 1."""
    seen = {"route": [], "ms": [], "g_mom": []}
    inner_r58 = bt._t1_r58
    inner_ms, inner_g, ms_step, g_mom = _ms_calls(seen)

    def r58(a):
        values = inner_r58(a)
        if values is None:
            return None
        seen["route"].append(values["file"])
        return {**values, "direction_momentum_damping": float(bt.harness.B0_DAMPING)}

    bt._t1_r58, bt._t1_ms_step, bt._t1_g_mom = r58, ms_step, g_mom
    try:
        bt.main(argv)
    finally:
        bt._t1_r58, bt._t1_ms_step, bt._t1_g_mom = inner_r58, inner_ms, inner_g
    fraction = bt.T1_MSB_MULT * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    ok = (seen["route"] == [bt.T1_R58_FILE] and len(seen["ms"]) == 1 and seen["ms"][0][0] == fraction
          and seen["ms"][0][1].get("damping") == float(bt.harness.B0_DAMPING)
          and seen["g_mom"] == [bt.G_MOM_DAMPINGS[0]])
    return 0 if ok else 1


def lm16own(argv: list[str]) -> int:
    """Row R54's file with the moved shares put back at the pool's own mix (module docstring); 0 or 1."""
    seen = {"route": [], "late": [], "ms": [], "g_mom": []}
    inner_r54, inner_late = bt._t1_r54, bt._t1_late_mix_shares
    inner_ms, inner_g, ms_step, g_mom = _ms_calls(seen)

    def r54(a):
        values = inner_r54(a)
        if values is not None:
            seen["route"].append(values["file"])
        return values

    def late(own, values):
        seen["late"].append(dict(own) == dict(values["own_shares"]))
        return own, {"identity_check": "the pool's own chamber mix, unmoved (experiments/builds/t1_lm_m20_check.py)"}

    bt._t1_r54, bt._t1_late_mix_shares, bt._t1_ms_step, bt._t1_g_mom = r54, late, ms_step, g_mom
    try:
        bt.main(argv)
    finally:
        bt._t1_r54, bt._t1_late_mix_shares, bt._t1_ms_step, bt._t1_g_mom = inner_r54, inner_late, inner_ms, inner_g
    fraction = bt.T1_MSB_MULT * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    ok = (seen["route"] == [bt.T1_R54_FILE] and seen["late"] == [True] and len(seen["ms"]) == 1
          and seen["ms"][0] == (fraction, {}) and seen["g_mom"] == [bt.G_MOM_DAMPINGS[0]])
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] not in ("m20msb", "lm16own") or argv[1] != "t1":
        sys.exit(USAGE)
    return {"m20msb": m20msb, "lm16own": lm16own}[argv[0]](argv[1:])


if __name__ == "__main__":
    sys.exit(main())
