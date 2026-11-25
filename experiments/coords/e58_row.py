#!/usr/bin/env python3
"""Row E58 of experiments/harness/TRIES.md (gate G1 of row R59's file T2embryo_s2x_ovt23.h5ad), read by code: the
row's form, the read's fixed outcome words and the committed reading. Standard library only, so the launchers run it
with the system python3 before anything else (`python3 -I -B experiments/coords/e58_row.py ...`); the read script
(experiments/coords/e58_mix_lean_read.py) imports the same constants.

  row TRIES.md                exit 0 when TRIES.md holds exactly one row E58, of ROW_CELLS cells, whose first cell
                              starts with FIRST_CELL_OPENING (the row says in its first line that it moves the overlay's
                              form); exit 1 otherwise. Prints nothing on stdout.
  passed TRIES.md             as `row`, and the row's outcome cell (its last, cell OUTCOME_CELL) holds PASS_TEXT word
                              for word and neither "G1: STOP" nor "G1: VOID"; exit 1 otherwise
  reading READING.json TXT    exit 0 and print the run commit when the committed reading says PASS: its outcome is
                              "PASS" and its outcome_text PASS_TEXT, it is not void, its run commit is a full hash, its
                              reading code ran on a clean tree, and the text holds PASS_TEXT and no STOP or VOID line;
                              exit 1 otherwise
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROW_ID = "E58"
ROW_PREFIX = f"| {ROW_ID} "  # the row's line in TRIES.md starts with this
ROW_CELLS = 11  # every row of the ledger's tables has eleven cells
OUTCOME_CELL = 11  # 1-based: the last cell holds the outcome, in the read's fixed words
FIRST_CELL_OPENING = ("E58 (gate G1 of row R59's file T2embryo_s2x_ovt23.h5ad: it moves the overlay's form; one "
                      "counted try, the arm `fgw_shape0x_ov_t23`")
PASS_TEXT = ("G1: PASS. Row R59's file T2embryo_s2x_ovt23.h5ad is built as named, twice at one clean commit; a "
             "diagnostic (J5, noise floor 6.54).")
STOP_TEXT = "G1: STOP. No file; the handed-over embryo files stand."
IDLE_TEXT = ("G1: STOP (the idle case: the merged gain is below 0.3 in size). No file; the handed-over embryo files "
             "stand.")
VOID_TEXT = ("G1: VOID. Read as nothing and no file. No task or seed is re-run selectively, and nothing further "
             "follows under this row. A rerun needs a new row placed before it, and its try counts again. A try "
             "counted at the run step's start stays counted (rule 4).")
NOT_PASS = ("G1: STOP", "G1: VOID")
HEX40 = re.compile(r"[0-9a-f]{40}")


def cells(line: str) -> list[str] | None:
    """The cells of one table line, split on pipes that are not escaped; None if the line is not a table row."""
    parts = re.split(r"(?<!\\)\|", line)
    if len(parts) < 3 or parts[0].strip() or parts[-1].strip():
        return None
    return [p.strip() for p in parts[1:-1]]


def row_lines(path) -> list[str]:
    p = Path(path)
    if not p.is_file():
        return []
    return [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln.startswith(ROW_PREFIX)]


def row_problems(path) -> list[str]:
    """Why TRIES.md at `path` does not hold row E58 in the form the gate needs (empty: it does)."""
    lines = row_lines(path)
    if len(lines) != 1:
        return [f"TRIES.md holds {len(lines)} rows {ROW_ID}, not exactly one"]
    c = cells(lines[0])
    if c is None or len(c) != ROW_CELLS:
        return [f"row {ROW_ID} is not a table row of {ROW_CELLS} cells"]
    if not c[0].startswith(FIRST_CELL_OPENING):
        return [f"row {ROW_ID}'s first cell does not open with the gate's fixed words"]
    return []


def outcome_cell(path) -> str | None:
    """Row E58's outcome cell when the row is well formed, else None."""
    if row_problems(path):
        return None
    return cells(row_lines(path)[0])[OUTCOME_CELL - 1]


def passed(path) -> bool:
    """Whether row E58 is well formed and its outcome cell holds PASS_TEXT word for word and no STOP or VOID words."""
    cell = outcome_cell(path)
    return cell is not None and PASS_TEXT in cell and not any(w in cell for w in NOT_PASS)


def reading_commit(json_path, txt_path) -> str | None:
    """The run commit of a committed reading that says PASS (see the module's `reading`), else None."""
    try:
        r = json.loads(Path(json_path).read_text(encoding="utf-8"))
        text = Path(txt_path).read_text(encoding="utf-8")
    except (OSError, ValueError):
        return None
    if not isinstance(r, dict):
        return None
    code = r.get("reading_code") if isinstance(r.get("reading_code"), dict) else {}
    commit = str(r.get("commit"))
    ok = (r.get("void") is False and r.get("outcome") == "PASS" and r.get("outcome_text") == PASS_TEXT
          and r.get("passes") is True and HEX40.fullmatch(commit) is not None and code.get("clean") is True
          and code.get("head") == commit and PASS_TEXT in text.splitlines() and not any(w in text for w in NOT_PASS))
    return commit if ok else None


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) == 2 and argv[0] == "row":
        return 0 if not row_problems(argv[1]) else 1
    if len(argv) == 2 and argv[0] == "passed":
        return 0 if passed(argv[1]) else 1
    if len(argv) == 3 and argv[0] == "reading":
        commit = reading_commit(argv[1], argv[2])
        if commit is None:
            return 1
        print(commit)
        return 0
    print("usage: e58_row.py row TRIES.md | passed TRIES.md | reading READING.json READING.txt", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
