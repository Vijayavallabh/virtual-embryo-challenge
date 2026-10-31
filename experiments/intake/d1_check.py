#!/usr/bin/env python
"""Decision point D1: what each final-phase test board publishes, and assumption A3 (docs/P3_RUNBOOK.md section 5).

Run after `vec boards update` has written the test boards to boards/ (anchors stripped in memory at intake). Reads
only submission-spec fields from boards/index.json: never anchors (the intake drops them, and the script refuses an
index in which a non-validation board still carries any) and never the public board file itself.

A3 (ASSUMPTIONS.md): each test board uses the last released stage as copy_last and DE reference (T1 and T2 heart:
E10.5; T2 embryo: E7.5; T3: the matched wild type E8.75) and publishes truth_cells. The server's working copies are
10% of a stage (runbook step 4.6: round(N * 0.1) rows), so a board's ref_cells names its reference stage when it
equals 10% of exactly one released stage's row count. Rule (fixed 2026-10-31, before any test spec was read): A3
passes for a board if its ref_cells matches the expected stage and truth_cells is published; it fails if ref_cells
matches another released stage (the builders' deltas must then be rebuilt against that stage); it is "unresolved" if
ref_cells is missing or matches no released stage. The test boards' cell limits, gene lists and required keys are
listed for the builders, and each gene list is compared with the validation board of the same task and setting.

A test board's truth_cells, min_cells and max_cells can encode the held-out stage's cell count (on the validation
boards max_cells is 3 x the larger of truth_cells and ref_cells, and min_cells follows truth_cells below 1,000), so
this script never prints or stores them: it reports only whether each is published. The builders read them from
boards/ in code; nobody quotes them in docs, notes or agent output (red-team review, 2026-10-31). ref_cells is 10% of
a released stage and is reported.

Usage (local or DGX, repo root):  uv run python experiments/intake/d1_check.py --out experiments/intake/results
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from vec.boards import REPO_ROOT, VALIDATION_BOARDS, load_boards

TEST_BOARDS = ("T1:test", "T2:heart:test_extrap", "T2:embryo:test_interp", "T3:bcatenin")
VALIDATION_TWIN = {"T1:test": "T1:val", "T2:heart:test_extrap": "T2:heart:val_extrap",
                   "T2:embryo:test_interp": "T2:embryo:val_interp", "T3:bcatenin": "T3:gata4"}
# Released row counts (docs/DATA_NOTES.md), by task setting.
RELEASED = {
    "T1": {"E8.5": 16787, "E9.5": 17057, "E10.5": 24979},
    "T2:heart": {"E8.25_late": 58716, "E8.5": 23710, "E8.75": 24826, "E9.5": 53742, "E10.5": 167864},
    "T2:embryo": {"E6.75": 7093, "E7.25": 13295, "E7.5": 11665, "E8.0": 31671},
    "T3": {"WT E8.75": 24826, "Mab21l2 KO E9.5": 50294, "Gata4 KO 1 E8.75": 42257, "Gata4 KO 2 E8.75": 46396,
           "WT E9.5": 53742},
}
EXPECTED_REF = {"T1:test": "E10.5", "T2:heart:test_extrap": "E10.5", "T2:embryo:test_interp": "E7.5",
                "T3:bcatenin": "WT E8.75"}


def ten_percent(n: int) -> int:
    return max(2, int(round(n * 0.1)))


def family(board) -> str:
    return board.task if board.task in ("T1", "T3") else f"{board.task}:{board.setting}"


def check(boards: dict) -> dict:
    leaking = [k for k, b in boards.items() if k not in VALIDATION_BOARDS and b.anchors]
    if leaking:  # test-board anchors encode hidden targets and are never read
        sys.exit(f"refused: {len(leaking)} non-validation board(s) still carry anchors; re-run the intake")
    out = {}
    for key in TEST_BOARDS:
        b = boards.get(key)
        if b is None:
            out[key] = {"present": False, "A3": "unresolved", "why": "board not in boards/index.json"}
            continue
        stages = RELEASED[family(b)]
        hits = sorted(s for s, n in stages.items() if b.ref_cells is not None and ten_percent(n) == b.ref_cells)
        if b.ref_cells is None or not hits:
            a3 = "unresolved"
            why = "ref_cells missing" if b.ref_cells is None else "ref_cells matches no released stage"
        elif hits == [EXPECTED_REF[key]] and b.truth_cells is not None:
            a3, why = "pass", f"reference stage {hits[0]}; truth_cells published"
        elif EXPECTED_REF[key] in hits and b.truth_cells is not None:
            a3, why = "pass", f"ref_cells fits {hits} (same count); truth_cells published"
        elif EXPECTED_REF[key] in hits:
            a3, why = "fail", "truth_cells not published"
        else:
            a3, why = "fail", f"reference stage looks like {hits}, not {EXPECTED_REF[key]}"
        twin = boards.get(VALIDATION_TWIN[key])
        out[key] = {
            "present": True, "task": b.task, "setting": b.setting, "split": b.split, "target_ko": b.target_ko,
            "n_genes": b.n_genes, "needs_coords": b.needs_coords, "ref_cells": b.ref_cells,
            "floor_model": b.floor_model,
            "published": {"min_cells": b.min_cells is not None, "max_cells": b.max_cells is not None,
                          "truth_cells": b.truth_cells is not None},
            "ref_stage_hits": hits, "expected_ref_stage": EXPECTED_REF[key], "A3": a3, "why": why,
            "genes_same_as_validation_twin": None if twin is None else b.genes == twin.genes,
            "genes_same_set_as_validation_twin": None if twin is None else set(b.genes) == set(twin.genes),
        }
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--boards-dir", type=Path, default=REPO_ROOT / "boards")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "experiments" / "intake" / "results")
    a = ap.parse_args(argv)
    res = check(load_boards(a.boards_dir))
    verdicts = {k: v["A3"] for k, v in res.items()}
    overall = ("pass" if all(v == "pass" for v in verdicts.values())
               else "fail" if "fail" in verdicts.values() else "unresolved")
    summary = {"utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), "A3": overall, "boards": res}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "d1_check.json").write_text(json.dumps(summary, indent=1))
    lines = [f"D1 / A3: {overall}"] + [
        f"  {k}: A3 {v['A3']} ({v['why']})" + (
            f"; cell limits and truth_cells published {v['published']} (values not shown), ref_cells {v['ref_cells']}, "
            f"{v['n_genes']} genes (same list as {VALIDATION_TWIN[k]}: {v['genes_same_as_validation_twin']}), "
            f"coords {v['needs_coords']}, floor {v['floor_model']}" if v["present"] else "")
        for k, v in res.items()]
    (a.out / "d1_check.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
