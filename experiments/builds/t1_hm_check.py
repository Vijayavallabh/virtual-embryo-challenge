#!/usr/bin/env python
"""Target-free checks of the files the T1 design round of 2026-11-24 puts on the banked command (T1_s2rcB25H05, H10 and
H20: Blood's released momentum, --held-momentum; T1_s2rcB25MS: MS at 2.5). experiments/builds/t1_hm_line.sh runs them
on the test board after its builds. Each check prints booleans only (a name and True or False per line, never a value,
count or shape) and exits 0 when every boolean is True, else 1; a refusal of the builder exits as the builder does.

  one ARGS...           build_test.py's own main on ARGS (a build: `t1 ... --out <file>`), with the damping of every
                        vec.t1_scope.momentum_shift call forced to 0, so the momentum's factors are exactly 1 (Blood's
                        step in a --held-momentum command, MS's step in the MS command); exits 1 unless that step ran
                        exactly once. The launcher compares the file with the banked file's sha256: with the factors
                        at 1 the file must be the banked file byte for byte.
  rows BANKED H...      the banked command's draw, by the builder's own calls (build_test.py's main on BANKED_ARGV at
                        the default seed on the test board, stopped inside build_submission at the population g_mom's
                        hook would get: nothing is written), then: the banked file holds the drawn Blood cells at their
                        drawn counts (so the draw's rows are the files' rows), and for each H file every row but the
                        drawn Blood cells' is the banked file's byte for byte and at least one Blood row moved.
  record H V [H V ...]  each H file's build record (beside it, <file>.build.json): only its t1_candidate.momentum fields
                        and its value are read: the fit runs from E9.5 to E10.5's libraries 1 and 2 on Blood-labelled
                        cells at both ends (library 3 out), Blood shifted (not skipped), at damping V, on stream 6.

Run from the repository root as the launcher does (uv run --frozen --no-sync python -I -B
experiments/builds/t1_hm_check.py ...); nothing here uploads or deletes anything.
"""

from __future__ import annotations

import json
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

# The banked file's command (T1_s2rcB25; t1_hm_line.sh's X25_ARGS after "t1"), at the builder's default seed and board
BANKED_ARGV = ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "recomputed", "--feng-fraction-mult", "2.5",
               "--captain-ruling", bt.T1_R22_ROW]
USAGE = ("usage: t1_hm_check.py one t1 ARGS... --out FILE | rows BANKED_FILE H_FILE... | record H_FILE VALUE "
         "[H_FILE VALUE ...]")


class _Drawn(Exception):
    """Raised inside build_submission once the drawn population is captured, so that nothing is written."""


def factors_one(argv: list[str]) -> int:
    """build_test.main(argv) with every t1_scope.momentum_shift call's damping forced to 0 (log f = 0: every factor
    exactly 1). Returns 0 when the build ran through and the momentum step was called exactly once, else 1."""
    inner = t1_scope.momentum_shift
    calls = []

    def one(start, ref, target_like, damping, *a, **k):
        calls.append(1)
        return inner(start, ref, target_like, 0.0, *a, **k)

    t1_scope.momentum_shift = one
    try:
        bt.main(argv)
    finally:
        t1_scope.momentum_shift = inner
    return 0 if len(calls) == 1 else 1


def banked_population() -> ad.AnnData:
    """The population the banked command's post_cells hook gets (the builder's draw: library 3's cells at its own
    chamber mix, the count rule, the default seed, the test board; then the atlas step), from build_test.main on
    BANKED_ARGV, stopped there: no file is written (the output path is in a fresh temporary folder)."""
    inner = t1_scope.build_submission
    seen = {}

    def capture(*a, **k):
        def stop(pop):
            seen["pop"] = pop.copy()
            raise _Drawn

        return inner(*a, **{**k, "post_cells": stop})

    out = Path(tempfile.mkdtemp(prefix="t1hm_draw_")) / "T1_draw.h5ad"
    t1_scope.build_submission = capture
    try:
        bt.main([*BANKED_ARGV, "--out", str(out)])
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


def rows(banked_path: str, h_paths: list[str]) -> dict[str, bool]:
    """The rows check (module docstring). The drawn Blood cells' expected rows in the banked file are their drawn
    counts as the file writes every row: log1p CP10k over all genes (vec.core.log1p_cp10k, row by row), the board's
    genes, float32, clipped at 0 (vec.io.make_submission)."""
    pop = banked_population()
    board = bt.load_boards()[bt.TEST_BOARDS["t1"]]
    blood = pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy() == bt.T1_HM_LABEL
    cols = pop.var_names.get_indexer(list(board.genes))
    counts = sparse.csr_matrix(pop.layers["counts"])[np.flatnonzero(blood)]
    want = np.clip(dense(log1p_cp10k(counts)[:, cols]).astype(np.float32, copy=False), 0.0, None)
    banked = ad.read_h5ad(banked_path)
    xb = _x(banked)
    res = {"banked_file_holds_the_draw": (xb.shape[0] == pop.n_obs and (cols >= 0).all()
                                          and list(banked.var_names) == list(board.genes)),
           "blood_drawn": bool(blood.any() and not blood.all())}
    res["banked_blood_rows_at_drawn_counts"] = res["banked_file_holds_the_draw"] and _same(xb[blood], want)
    for p in h_paths:
        h = ad.read_h5ad(p)
        xh = _x(h)
        name = Path(p).stem
        frame = (xh.shape == xb.shape and list(h.var_names) == list(banked.var_names)
                 and list(h.obs_names) == list(banked.obs_names))
        res[f"{name}_frame_is_the_banked_files"] = frame
        res[f"{name}_other_rows_are_the_banked_files"] = frame and _same(xh[~blood], xb[~blood])
        res[f"{name}_a_blood_row_moved"] = frame and bool((xh[blood] != xb[blood]).any())
    return res


def record(pairs: list[tuple[str, str]]) -> dict[str, bool]:
    """The record check (module docstring): for each (H file, value), only its record's t1_candidate.held_momentum and
    t1_candidate.momentum fields are read."""
    res = {}
    for path, value in pairs:
        name = Path(path).stem
        cand = json.loads(Path(path).with_suffix(".build.json").read_text()).get("t1_candidate") or {}
        m = cand.get("momentum") or {}
        fit = m.get("fit") or {}
        res[f"{name}_fit_e9_5_to_e10_5_libraries_1_and_2"] = (
            fit.get("start") == "E9.5_RNA" and fit.get("later") == "E10.5_RNA"
            and fit.get("later_libraries") == list(bt.T1_MS_FIT_LIBRARIES))
        res[f"{name}_library_3_out"] = fit.get("library_3_in_fit") is False
        res[f"{name}_blood_cells_only"] = (
            fit.get("within") == {bt.T1_HM_CATEGORY: [bt.T1_HM_LABEL]} and m.get("label") == bt.T1_HM_LABEL
            and fit.get("start_cells") == "own, restricted" and fit.get("later_cells") == "own, restricted")
        res[f"{name}_blood_shifted"] = m.get("status") == "shifted"
        res[f"{name}_its_damping"] = (value in bt.T1_HM_CHOICES and cand.get("held_momentum") == float(value)
                                      and m.get("damping") == float(value))
        res[f"{name}_its_own_stream"] = list(m.get("stream") or [])[1:] == [bt.T1_HM_STREAM]
    return res


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in ("one", "rows", "record"):
        sys.exit(USAGE)
    cmd, rest = argv[0], argv[1:]
    if cmd == "one":
        if not rest or rest[0] != "t1":
            sys.exit(USAGE)
        return factors_one(rest)
    if cmd == "rows":
        if len(rest) < 2:
            sys.exit(USAGE)
        res = rows(rest[0], rest[1:])
    else:
        if len(rest) < 2 or len(rest) % 2:
            sys.exit(USAGE)
        res = record(list(zip(rest[::2], rest[1::2], strict=True)))
    for k, v in res.items():
        print(k, bool(v))
    return 0 if res and all(res.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
