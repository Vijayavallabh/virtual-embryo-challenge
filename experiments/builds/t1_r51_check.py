#!/usr/bin/env python
"""Target-free checks of row R51's seven files (T1's line 10, each one change to T1_s2rcMSB's command: T1_s2rcMS20,
T1_s2rcMS30, T1_s2rcMSBd10, T1_s2rcMSBd15, T1_s2rcMSBH05, T1_s2rcMSBH10, T1_s2rcMSBm10). experiments/builds/
t1_r51_line.sh runs them on the test board after its identity builds. Each check prints booleans only (a name and True
or False per line, never a value, count or shape) or nothing, and exits 0 when it passes, else 1; a refusal of the
builder exits as the builder does.

  msb ARGS...     build_test.py's own main on ARGS (one of row R51's seven commands: `t1 ... --out <file>`) on the
                  file's own path (row R51's route, its MS step, Blood's step and g_mom's hook as the file has them),
                  with the one setting the file changes put at T1_s2rcMSB's value: the multiplier at 1.5 (T1_s2rcMS20,
                  T1_s2rcMS30), g_mom's damping at 0.5 (T1_s2rcMSBd10, T1_s2rcMSBd15), MS's damping at
                  harness.B0_DAMPING (T1_s2rcMSBm10), or Blood's factors at 1, the damping of Blood's momentum_shift
                  call forced to 0 (T1_s2rcMSBH05, T1_s2rcMSBH10). Exits 1 unless the command took row R51's route
                  once, MS's step ran once at 1.5 x the calendar fraction and MS's momentum once at B0's damping,
                  g_mom's hook was made once at damping 0.5, and Blood's step ran once at the file's damping in the two
                  H files and never in the others. The launcher compares the file with T1_s2rcMSB's sha256: with its one
                  setting at MSB's value, each new path must give MSB's file byte for byte.
  rows MSB H...   T1_s2rcMSB's draw, by the builder's own calls (build_test.py's main on MSB_ARGV at the default seed on
                  the test board, stopped inside build_submission at the population its post_cells would get: nothing
                  is written), then: MSB's file holds the drawn Blood cells at their drawn counts (so the draw's rows
                  are the files' rows), and for each H file every row but the drawn Blood cells' is MSB's byte for
                  byte and at least one Blood row moved.

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/t1_r51_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I puts no script folder on sys.path
import build_test as bt  # noqa: E402

from vec import t1_scope  # noqa: E402
from vec.core import log1p_cp10k  # noqa: E402
from vec.io import dense  # noqa: E402

# T1_s2rcMSB's command (row R26's; t1_r51_line.sh's MSB_ARGS after "t1"), at the builder's default seed and board
MSB_ARGV = ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "recomputed", "--feng-fraction-mult", "1.5",
            "--direction-momentum", "all", "--captain-ruling", bt.T1_R26_ROW]
# file -> the one setting it changes (a key of _t1_r51's values)
SETTING = {"T1_s2rcMS20": "feng_fraction_mult", "T1_s2rcMS30": "feng_fraction_mult",
           "T1_s2rcMSBd10": "g_mom_damping", "T1_s2rcMSBd15": "g_mom_damping",
           "T1_s2rcMSBH05": "held_momentum", "T1_s2rcMSBH10": "held_momentum",
           "T1_s2rcMSBm10": "direction_momentum_damping"}
# MSB's value of each setting but Blood's (MSB has no Blood step: its factors are put at 1 instead)
MSB_VALUES = {"feng_fraction_mult": bt.T1_MSB_MULT, "g_mom_damping": bt.G_MOM_DAMPINGS[0],
              "direction_momentum_damping": float(bt.harness.B0_DAMPING)}
BLOOD_WITHIN = {bt.T1_HM_CATEGORY: [bt.T1_HM_LABEL]}  # Blood's momentum_shift call, and no other, restricts to it
USAGE = "usage: t1_r51_check.py msb t1 ARGS... --out FILE | rows MSB_FILE H_FILE..."


class _Drawn(Exception):
    """Raised inside build_submission once the drawn population is captured, so that nothing is written."""


def at_msb(argv: list[str]) -> int:
    """build_test.main(argv) on row R51's route with the file's one setting at MSB's value (module docstring). Returns
    0 when the build ran through with every call as named there, else 1."""
    seen = {"files": [], "ms": [], "g_mom": [], "blood": [], "other": []}
    inner_r51, inner_ms, inner_g, inner_shift = bt._t1_r51, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift

    def r51(a):
        values = inner_r51(a)
        if values is None:
            return None
        seen["files"].append((values["file"], values["held_momentum"]))
        key = SETTING[values["file"]]
        return values if key == "held_momentum" else {**values, key: MSB_VALUES[key]}

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

    bt._t1_r51, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift = r51, ms_step, g_mom, shift
    try:
        bt.main(argv)
    finally:
        bt._t1_r51, bt._t1_ms_step, bt._t1_g_mom, t1_scope.momentum_shift = inner_r51, inner_ms, inner_g, inner_shift
    if len(seen["files"]) != 1:
        return 1
    name, held = seen["files"][0]
    msb_fraction = MSB_VALUES["feng_fraction_mult"] * bt.feng_dirs.bracket_fraction(10.5, 12.5)
    ok = (seen["ms"] == [(msb_fraction, MSB_VALUES["direction_momentum_damping"])]
          and seen["other"] == [MSB_VALUES["direction_momentum_damping"]]
          and seen["g_mom"] == [MSB_VALUES["g_mom_damping"]]
          and seen["blood"] == ([held] if SETTING[name] == "held_momentum" else []))
    return 0 if ok else 1


def msb_population() -> ad.AnnData:
    """The population MSB's post_cells gets (the builder's draw: library 3's cells at its own chamber mix, the count
    rule, the default seed, the test board; before MS's step), from build_test.main on MSB_ARGV, stopped there: no file
    is written (the output path is in a fresh temporary folder)."""
    inner = t1_scope.build_submission
    seen = {}

    def capture(*a, **k):
        def stop(pop):
            seen["pop"] = pop.copy()
            raise _Drawn

        return inner(*a, **{**k, "post_cells": stop})

    out = Path(tempfile.mkdtemp(prefix="t1r51_draw_")) / "T1_draw.h5ad"
    t1_scope.build_submission = capture
    try:
        bt.main([*MSB_ARGV, "--out", str(out)])
    except _Drawn:
        pass
    finally:
        t1_scope.build_submission = inner
    if "pop" not in seen or out.exists():
        raise SystemExit(1)
    return seen["pop"]


def _x(a: ad.AnnData) -> np.ndarray:
    return np.ascontiguousarray(dense(a.X))


def _same(p: np.ndarray, q: np.ndarray) -> bool:
    p, q = np.ascontiguousarray(p), np.ascontiguousarray(q)
    return p.dtype == q.dtype and p.shape == q.shape and p.tobytes() == q.tobytes()


def rows(msb_path: str, h_paths: list[str]) -> dict[str, bool]:
    """The rows check (module docstring). The drawn Blood cells' expected rows in MSB's file are their drawn counts as
    the file writes every row: log1p CP10k over all genes (vec.core.log1p_cp10k, row by row), the board's genes,
    float32, clipped at 0 (vec.io.make_submission)."""
    pop = msb_population()
    board = bt.load_boards()[bt.TEST_BOARDS["t1"]]
    blood = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy() == bt.T1_HM_LABEL
    cols = pop.var_names.get_indexer(list(board.genes))
    counts = sparse.csr_matrix(pop.layers["counts"])[np.flatnonzero(blood)]
    want = np.clip(dense(log1p_cp10k(counts)[:, cols]).astype(np.float32, copy=False), 0.0, None)
    msb = ad.read_h5ad(msb_path)
    xm = _x(msb)
    res = {"msb_file_holds_the_draw": (xm.shape[0] == pop.n_obs and (cols >= 0).all()
                                       and list(msb.var_names) == list(board.genes)),
           "blood_drawn": bool(blood.any() and not blood.all())}
    res["msb_blood_rows_at_drawn_counts"] = res["msb_file_holds_the_draw"] and _same(xm[blood], want)
    for p in h_paths:
        h = ad.read_h5ad(p)
        xh = _x(h)
        name = Path(p).stem
        frame = (xh.shape == xm.shape and list(h.var_names) == list(msb.var_names)
                 and list(h.obs_names) == list(msb.obs_names))
        res[f"{name}_frame_is_msbs"] = frame
        res[f"{name}_other_rows_are_msbs"] = frame and _same(xh[~blood], xm[~blood])
        res[f"{name}_a_blood_row_moved"] = frame and bool((xh[blood] != xm[blood]).any())
    return res


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in ("msb", "rows"):
        sys.exit(USAGE)
    cmd, rest = argv[0], argv[1:]
    if cmd == "msb":
        if not rest or rest[0] != "t1":
            sys.exit(USAGE)
        return at_msb(rest)
    if len(rest) < 2:
        sys.exit(USAGE)
    res = rows(rest[0], rest[1:])
    for k, v in res.items():
        print(k, bool(v))
    return 0 if res and all(res.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
