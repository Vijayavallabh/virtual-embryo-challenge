#!/usr/bin/env python3
"""Rows H37 and R69 of experiments/harness/TRIES.md, read by code (the standard library only, so the launchers run it
with python3 -I -B before any virtual environment is touched).

  h37_row.py check TRIES     exit 0 when the ledger holds exactly one row H37 of eleven cells whose first cell is
                             FIRST_CELL and whose outcome cell holds no outcome yet; 3 when the row is missing or
                             malformed; 4 when its outcome cell already holds an outcome
  h37_row.py passed TRIES    exit 0 and print the PASS label ("world pattern" or "spread only") when row H37 is
                             formed
                             and its outcome cell holds exactly one outcome statement and that statement is the read's
                             PASS line; 3 malformed; 5 VOID; 6 STOP; 4 anything else
  h37_row.py ruling TRIES    exit 0 when the ledger holds exactly one row R69 whose first cell is the captain's order
                             ("ruling, no run") naming file WBH, its gate row H37 and `--captain-ruling R69`; 3
                             otherwise

experiments/harness/h37_heart_world_blur.py, experiments/builds/world_blur.py and both launchers read the rows
through this module and nothing else.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROW = "H37"
RULING_ROW = "R69"
# Row H37's first cell, exactly: every step of the gate script, its launcher (run_h37.sh), the file's builder
# (experiments/builds/world_blur.py) and the build launcher (heart_wbh_build.sh) refuse a row H37 whose first cell is
# not this text
FIRST_CELL = ("H37 (the gate of row R69's file WBH, the heart's world-blur hedge: six counted tries, arms W, U and P "
              "on two transitions; G0 first, target-free, not a try)")
N_CELLS = 11
OUTCOME_CELL = N_CELLS - 1
MARK = f"OUTCOME {ROW}:"
PASS_LABELS = ("world pattern", "spread only")
PASS_RE = re.compile(rf"OUTCOME {ROW}: PASS \(({'|'.join(PASS_LABELS)})\)\.")
STOP_RE = re.compile(rf"OUTCOME {ROW}: STOP(:| at) ")
# what row R69's first cell must hold (the captain's order of 2026-11-25)
RULING_NAMES = ("File WBH", f"row {ROW}", f"`--captain-ruling {RULING_ROW}`")


def rows(text: str, row: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln.startswith(f"| {row} ")]


def cells(text: str, row: str = ROW) -> list[str] | None:
    """The row's cells when the text holds exactly one such row of N_CELLS cells; None otherwise."""
    found = rows(text, row)
    if len(found) != 1:
        return None
    ln = found[0].rstrip()
    if not ln.endswith("|"):
        return None
    out = [c.strip() for c in ln[1:-1].split("|")]
    return out if len(out) == N_CELLS else None


def formed(text: str) -> bool:
    c = cells(text)
    return c is not None and c[0] == FIRST_CELL


def is_open(text: str) -> bool:
    return formed(text) and MARK not in cells(text)[OUTCOME_CELL]


def outcome(text: str) -> tuple[str, str | None] | None:
    """(word, label) from row H37's outcome cell: ("PASS", label), ("STOP", None), ("VOID", None) or ("OTHER", None);
    None when the row is not formed."""
    if not formed(text):
        return None
    cell = cells(text)[OUTCOME_CELL]
    if cell.count(MARK) != 1:
        return ("OTHER", None)
    if f"{MARK} VOID" in cell:
        return ("VOID", None)
    if STOP_RE.search(cell):
        return ("STOP", None)
    m = PASS_RE.search(cell)
    return ("PASS", m.group(1)) if m else ("OTHER", None)


def ruling_ok(text: str) -> bool:
    """Exactly one row R69 whose first cell is the captain's order, ruling, no run, naming RULING_NAMES."""
    found = rows(text, RULING_ROW)
    if len(found) != 1:
        return False
    first = found[0][2:].split("|", 1)[0]
    return (first.startswith(f"{RULING_ROW} (the captain's order") and "ruling, no run" in first
            and all(s in first for s in RULING_NAMES))


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2 or argv[0] not in ("check", "passed", "ruling"):
        print("usage: h37_row.py check|passed|ruling TRIES", file=sys.stderr)
        return 2
    try:
        text = Path(argv[1]).read_text(encoding="utf-8")
    except OSError:
        text = ""
    if argv[0] == "ruling":
        return 0 if ruling_ok(text) else 3
    if not formed(text):
        print(f"row {ROW} is missing or malformed (one row of {N_CELLS} cells, its first cell the fixed text)",
              file=sys.stderr)
        return 3
    if argv[0] == "check":
        return 0 if is_open(text) else 4
    word, label = outcome(text)
    if word == "PASS":
        print(label)
        return 0
    print(f"row {ROW}'s outcome cell holds {word}, not the read's PASS line", file=sys.stderr)
    return {"VOID": 5, "STOP": 6}.get(word, 4)


if __name__ == "__main__":
    sys.exit(main())
