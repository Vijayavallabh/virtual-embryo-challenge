#!/usr/bin/env python3
"""Row H38 of experiments/harness/TRIES.md as code reads it (standard library only, so the launchers can ask it with the
system python3): whether the row is there and well formed, and what its outcome cell (the eleventh) says. Used by
experiments/heart_ext/run_oic_g0.sh (its check), experiments/harness/h38_heart_insitu.py and run_h38.sh (every step's
refusals) and experiments/builds/heart_oic_build.sh (its gate: OIC125 is built only when the outcome cell holds G0's
PASS statement and, after it, the read's PASS statement, and nothing that ends the line).

  python3 -I -B experiments/harness/h38_row.py check <TRIES.md>   exit 0 when row H38 is there, well formed, and its
                                                                  outcome cell holds no statement yet (before G0)
  python3 -I -B experiments/harness/h38_row.py g0 <TRIES.md>      print the three chords' sha256 (E10.5, E9.5, E8.75)
                                                                  and exit 0 when the outcome cell holds exactly G0's
                                                                  PASS statement and no read outcome (before the read)
  python3 -I -B experiments/harness/h38_row.py passed <TRIES.md>  print E10.5's chord sha256 and exit 0 when the cell
                                                                  holds G0's PASS statement, then "OUTCOME H38: PASS."
                                                                  and no statement that ends the line (the build)

Well formed: exactly one table line opens with "| H38 "; it closes with "|"; it has CELLS cells (no "|" inside a
cell); its first cell is FIRST_CELL exactly and its second "heart". Every refusal names its reason on stderr and shows
no figure.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROW = "H38"
CELLS = 11
BOARD_CELL = "heart"
FILE = "OIC125"
RULING = "R70"
FIRST_CELL = ("H38 (one row, one script; the counted read of heart file OIC125 under ruling row R70, ten counted heart "
              "tries on one backtest truth, a labelled diagnostic under rule 9 J5; G0 first, target-free, not a try)")
SHA = "[0-9a-f]{64}"
G0_FORM = re.compile(rf"G0 {ROW}: PASS\. chords sha256: E10\.5 ({SHA}); E9\.5 ({SHA}); E8\.75 ({SHA})\.")
G0_MARK = f"G0 {ROW}:"
OUTCOME_MARK = f"OUTCOME {ROW}:"
OUTCOMES = ("PASS", "STOP (IDLE)", "STOP", "VOID")
OUTCOME_FORM = re.compile(rf"OUTCOME {ROW}: (PASS|STOP \(IDLE\)|STOP|VOID)\.")


def g0_statement(e10: str, e95: str, e875: str) -> str:
    return f"G0 {ROW}: PASS. chords sha256: E10.5 {e10}; E9.5 {e95}; E8.75 {e875}."


def outcome_statement(word: str) -> str:
    if word not in OUTCOMES:
        raise ValueError(f"unknown outcome {word!r}")
    return f"OUTCOME {ROW}: {word}."


def row_lines(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if ln.startswith(f"| {ROW} ")]


def row_cells(text: str) -> list[str] | None:
    lines = row_lines(text)
    if len(lines) != 1 or not lines[0].rstrip().endswith("|"):
        return None
    return [c.strip() for c in lines[0].strip()[1:-1].split("|")]


def row_reasons(text: str) -> list[str]:
    """Why row H38 is missing or malformed (empty when it is there and well formed)."""
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
    if not cells or cells[0] != FIRST_CELL:
        why.append(f"row {ROW}'s first cell is not the wording fixed in experiments/harness/h38_row.py")
    if len(cells) < 2 or cells[1] != BOARD_CELL:
        why.append(f"row {ROW}'s second cell is not {BOARD_CELL!r}")
    return why


def outcome_cell(text: str) -> str | None:
    return None if row_reasons(text) else row_cells(text)[CELLS - 1]


def before_g0(text: str) -> list[str]:
    """Reasons the row is not ready for G0: malformed, or its outcome cell already holds a statement."""
    why = row_reasons(text)
    if why:
        return why
    cell = outcome_cell(text)
    if G0_MARK in cell or OUTCOME_MARK in cell:
        return ["the outcome cell already holds a statement (G0 is run once)"]
    return []


def g0_shas(text: str) -> tuple[tuple[str, str, str] | None, list[str]]:
    """The three chords' sha256 when the outcome cell holds exactly one G0 PASS statement, no other G0 statement and
    no read outcome yet."""
    why = row_reasons(text)
    if why:
        return None, why
    cell = outcome_cell(text)
    found = G0_FORM.findall(cell)
    if len(found) != 1 or cell.count(G0_MARK) != 1:
        why.append("the outcome cell does not hold exactly one G0 statement, and that one PASS")
    if OUTCOME_MARK in cell:
        why.append("the outcome cell already holds the read's outcome (the read is made once)")
    return (tuple(found[0]) if not why else None), why


def passed_sha256(text: str) -> tuple[str | None, list[str]]:
    """E10.5's chord sha256 when the cell holds G0's PASS statement, then exactly one read outcome, PASS."""
    why = row_reasons(text)
    if why:
        return None, why
    cell = outcome_cell(text)
    found = G0_FORM.findall(cell)
    outs = OUTCOME_FORM.findall(cell)
    if len(found) != 1 or cell.count(G0_MARK) != 1:
        why.append("the outcome cell does not hold exactly one G0 statement, and that one PASS")
    if cell.count(OUTCOME_MARK) != 1 or outs != ["PASS"]:
        why.append("the outcome cell does not hold exactly one read outcome, and that one PASS")
    if not why and cell.index(OUTCOME_MARK) < cell.index(G0_MARK):
        why.append("the read's outcome stands before G0's statement")
    return (found[0][0] if not why else None), why


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2 or argv[0] not in ("check", "g0", "passed"):
        print("usage: h38_row.py check|g0|passed <TRIES.md>", file=sys.stderr)
        return 2
    try:
        text = Path(argv[1]).read_text(encoding="utf-8")
    except OSError as e:
        print(f"refused: {argv[1]} cannot be read ({type(e).__name__})", file=sys.stderr)
        return 1
    if argv[0] == "check":
        why = before_g0(text)
        out = None
    elif argv[0] == "g0":
        shas, why = g0_shas(text)
        out = " ".join(shas) if shas else None
    else:
        out, why = passed_sha256(text)
    for w in why:
        print(f"refused: {w}", file=sys.stderr)
    if why:
        return 1
    if out:
        print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
