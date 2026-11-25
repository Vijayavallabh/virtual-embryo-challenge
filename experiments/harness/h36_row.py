#!/usr/bin/env python3
"""Row H36 of experiments/harness/TRIES.md as code reads it (standard library only, so that the launchers can ask it
with the system python3 before any environment is touched): whether the row is there and well formed, and what its
outcome cell (the eleventh) says. Used by experiments/harness/h36_heart_run_profile.py (every step's refusals),
experiments/harness/run_h36.sh (its check, exit 93) and experiments/builds/heart_xpc_build.sh (its gate: OXP125 is
built only when the outcome cell holds the test profile's PASS statement, read here).

  python3 -I -B experiments/harness/h36_row.py check <TRIES.md>    exit 0 when row H36 is there and well formed
  python3 -I -B experiments/harness/h36_row.py passed <TRIES.md>   print the profile's sha256 and exit 0 when the
                                                                   outcome cell holds exactly one final PASS statement
                                                                   and no statement that ends the line

Well formed: exactly one table line opens with "| H36 "; it closes with "|"; it has CELLS cells (no "|" inside a
cell); its first cell is FIRST_CELL exactly and its second "heart". Every refusal names its reason on stderr and shows
no figure.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROW = "H36"
CELLS = 11  # the ledger's columns: the row's ten cells and its outcome cell
BOARD_CELL = "heart"
FIRST_CELL = ("H36 (one row, one script; eleven counted heart arms on one backtest, plus the readings the eighth cell "
              "counts; a labelled diagnostic on one truth under rule 9 J5; runs only after ruling row R56 records the "
              "captain's named A15 exception for file OXP125; from the score-blind cold challenge of round_1125_gap)")
FILE = "OXP125"
# The read's statements (experiments/harness/h36_heart_run_profile.py writes them; the main session copies them into
# the outcome cell word for word)
PASSED_AT_PRACTICE = (f"Try {ROW}: PASSED AT PRACTICE: the test profile step (T) may follow once this outcome is "
                      "committed.")
PASS_FORM = re.compile(rf"Try {ROW}: PASS: {FILE} may be built with --run-profile ([0-9a-f]{{64}})\.")
# Openings of the statements that end the line (VOID, the stops and NOT PASSED): none may stand in the outcome cell of
# a file that is built
ENDINGS = (f"Try {ROW}: VOID", f"Try {ROW}: STOPPED", f"Try {ROW}: NOT PASSED", f"Try {ROW}: T STOPPED")


def pass_statement(sha256: str) -> str:
    return f"Try {ROW}: PASS: {FILE} may be built with --run-profile {sha256}."


def row_lines(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln.startswith(f"| {ROW} ")]


def row_cells(text: str) -> list[str] | None:
    """The cells of the one row H36 (None unless exactly one line opens it and closes with "|")."""
    lines = row_lines(text)
    if len(lines) != 1 or not lines[0].rstrip().endswith("|"):
        return None
    return [c.strip() for c in lines[0].strip()[1:-1].split("|")]


def row_reasons(text: str) -> list[str]:
    """Why row H36 is missing or malformed (empty when it is there and well formed)."""
    lines = row_lines(text)
    if not lines:
        return [f"no row {ROW} in TRIES.md"]
    if len(lines) > 1:
        return [f"{len(lines)} rows open with {ROW}, not one"]
    cells = row_cells(text)
    if cells is None:
        return [f"row {ROW} does not close with |"]
    why = []
    if len(cells) != CELLS:
        why.append(f"row {ROW} has {len(cells)} cells, not {CELLS} (no | may stand inside a cell)")
    if cells[0] != FIRST_CELL:
        why.append(f"row {ROW}'s first cell is not the wording fixed in experiments/harness/h36_row.py")
    if len(cells) < 2 or cells[1] != BOARD_CELL:
        why.append(f"row {ROW}'s second cell is not {BOARD_CELL!r}")
    return why


def outcome_cell(text: str) -> str | None:
    """The eleventh cell of a well-formed row H36 (None otherwise)."""
    return None if row_reasons(text) else row_cells(text)[CELLS - 1]


def passed_sha256(text: str) -> tuple[str | None, list[str]]:
    """The profile's sha256 when the outcome cell holds exactly one PASS statement, the practice read's PASSED AT
    PRACTICE before it, and no statement that ends the line; else None and the reasons."""
    why = row_reasons(text)
    if why:
        return None, why
    cell = outcome_cell(text)
    found = PASS_FORM.findall(cell)
    why = [f"the outcome cell holds {e!r}" for e in ENDINGS if e in cell]
    if len(found) != 1:
        why.append(f"the outcome cell holds {len(found)} PASS statements, not one")
    if PASSED_AT_PRACTICE not in cell:
        why.append("the outcome cell does not hold the practice read's PASSED AT PRACTICE statement")
    elif found and cell.index(PASSED_AT_PRACTICE) > cell.index(pass_statement(found[0])):
        why.append("the PASS statement stands before the practice read's statement")
    return (found[0] if not why else None), why


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2 or argv[0] not in ("check", "passed"):
        print("usage: h36_row.py check|passed <TRIES.md>", file=sys.stderr)
        return 2
    try:
        text = Path(argv[1]).read_text(encoding="utf-8")
    except OSError as e:
        print(f"refused: {argv[1]} cannot be read ({type(e).__name__})", file=sys.stderr)
        return 1
    if argv[0] == "check":
        why = row_reasons(text)
        for w in why:
            print(f"refused: {w}", file=sys.stderr)
        return 1 if why else 0
    sha, why = passed_sha256(text)
    for w in why:
        print(f"refused: {w}", file=sys.stderr)
    if sha is None:
        return 1
    print(sha)
    return 0


if __name__ == "__main__":
    sys.exit(main())
