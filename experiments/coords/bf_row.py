#!/usr/bin/env python3
"""Row E59 of experiments/harness/TRIES.md (the gate of row R71's files T2embryo_s2x_ovbf.h5ad and
T2embryo_s2xs_ovbf.h5ad, the embryo's consensus-form files BF on O and on SO), read by code: the row's form, the read's
fixed outcome words, the files a PASS names and the committed reading. Standard library only, so the launchers run it
with the system python3 before anything else (`python3 -I -B experiments/coords/bf_row.py ...`); the read script
(experiments/coords/bf_read.py) imports the same constants.

  row TRIES.md                exit 0 when TRIES.md holds exactly one row E59 of ROW_CELLS cells whose first cell opens
                              with FIRST_CELL_OPENING; exit 1 otherwise. Prints nothing.
  passed TRIES.md             as `row`, and the outcome cell (the last) holds exactly one of the PASS texts word for
                              word and no STOP or VOID words: prints the files that PASS names, one per line, exit 0;
                              exit 1 otherwise
  reading READING.json TXT    when the committed reading says PASS (its outcome a PASS key and its outcome_text that
                              key's words, not void, a full-hash run commit, its reading code at that commit on a clean
                              tree, the text holding the words and no STOP or VOID line): prints the run commit, then
                              the files, one per line, exit 0; exit 1 otherwise
  outcome READING.json        prints the reading's outcome line alone (what the main session reads); exit 1 if none
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROW_ID = "E59"
ROW_PREFIX = f"| {ROW_ID} "  # the row's line in TRIES.md starts with this
ROW_CELLS = 11  # every row of the ledger's tables has eleven cells
OUTCOME_CELL = 11  # 1-based: the last cell holds the outcome, in the read's fixed words
RULING_ROW = "R71"
FILE_O = "T2embryo_s2x_ovbf.h5ad"
FILE_SO = "T2embryo_s2xs_ovbf.h5ad"
FIRST_CELL_OPENING = (f"{ROW_ID} (the gate of row {RULING_ROW}'s files {FILE_O} (BF_O) and {FILE_SO} (BF_SO): it moves "
                      "the overlay's form, the occupancy term and the z extent; four counted tries")
TAG = "GBF"
DIAGNOSTIC = "a diagnostic (J5, noise floor 6.54)"
PASS_TEXTS = {
    "PASS_BOTH": (f"{TAG}: PASS for both files. Row {RULING_ROW}'s files {FILE_O} (BF_O) and {FILE_SO} (BF_SO) are "
                  f"built as named, each twice at one clean commit; {DIAGNOSTIC}."),
    "PASS_O": (f"{TAG}: PASS for {FILE_O} (BF_O) alone. Row {RULING_ROW}'s file {FILE_O} is built as named, twice at "
               f"one clean commit; {FILE_SO} is not built; {DIAGNOSTIC}."),
    "PASS_SO": (f"{TAG}: PASS for {FILE_SO} (BF_SO) alone. Row {RULING_ROW}'s file {FILE_SO} is built as named, twice "
                f"at one clean commit; {FILE_O} is not built; {DIAGNOSTIC}."),
}
PASS_FILES = {"PASS_BOTH": (FILE_O, FILE_SO), "PASS_O": (FILE_O,), "PASS_SO": (FILE_SO,)}
STOP_TEXT = f"{TAG}: STOP. No file; the handed-over embryo files stand."
IDLE_TEXT = (f"{TAG}: STOP (the idle case: neither arm's merged board gain reaches 0.3 in size). No file; the "
             "handed-over embryo files stand.")
G0_STOP_TEXT = (f"{TAG}: STOP at G0 (target-free, before any try). No file; the handed-over embryo files stand; "
                "nothing is counted.")
VOID_TEXT = (f"{TAG}: VOID. Read as nothing and no file. No task or seed is re-run selectively, and nothing further "
             "follows under this row. A rerun needs a new row placed before it, and its tries count again. Tries "
             "counted at the run step's start stay counted (rule 4).")
OUTCOME_TEXTS = {**PASS_TEXTS, "STOP": STOP_TEXT, "IDLE": IDLE_TEXT, "G0_STOP": G0_STOP_TEXT, "VOID": VOID_TEXT}
NOT_PASS = (f"{TAG}: STOP", f"{TAG}: VOID")
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
    """Why TRIES.md at `path` does not hold row E59 in the form the gate needs (empty: it does)."""
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
    """Row E59's outcome cell when the row is well formed, else None."""
    if row_problems(path):
        return None
    return cells(row_lines(path)[0])[OUTCOME_CELL - 1]


def pass_key(text: str) -> str | None:
    """The PASS key whose words the text holds, when it holds exactly one PASS text and no STOP or VOID words."""
    keys = [k for k, words in PASS_TEXTS.items() if words in text]
    if len(keys) != 1 or any(w in text for w in NOT_PASS):
        return None
    return keys[0]


def passed_files(path) -> tuple[str, ...]:
    """The files row E59's outcome cell names as passing (empty: none, or the row is not well formed)."""
    cell = outcome_cell(path)
    key = None if cell is None else pass_key(cell)
    return () if key is None else PASS_FILES[key]


def reading_pass(json_path, txt_path) -> tuple[str, tuple[str, ...]] | None:
    """(run commit, files) of a committed reading that says PASS (see the module's `reading`), else None."""
    try:
        r = json.loads(Path(json_path).read_text(encoding="utf-8"))
        text = Path(txt_path).read_text(encoding="utf-8")
    except (OSError, ValueError):
        return None
    if not isinstance(r, dict):
        return None
    code = r.get("reading_code") if isinstance(r.get("reading_code"), dict) else {}
    commit, key = str(r.get("commit")), r.get("outcome")
    ok = (key in PASS_TEXTS and r.get("outcome_text") == PASS_TEXTS[key] and r.get("void") is False
          and HEX40.fullmatch(commit) is not None and code.get("clean") is True and code.get("head") == commit
          and PASS_TEXTS[key] in text.splitlines() and not any(w in text for w in NOT_PASS)
          and sorted(r.get("files_passing") or []) == sorted(PASS_FILES[key]))
    return (commit, PASS_FILES[key]) if ok else None


def outcome_line(json_path) -> str | None:
    """The reading's outcome line, one of OUTCOME_TEXTS, else None."""
    try:
        r = json.loads(Path(json_path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    text = r.get("outcome_text") if isinstance(r, dict) else None
    return text if text in OUTCOME_TEXTS.values() else None


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) == 2 and argv[0] == "row":
        return 0 if not row_problems(argv[1]) else 1
    if len(argv) == 2 and argv[0] == "passed":
        files = passed_files(argv[1])
        if not files:
            return 1
        print("\n".join(files))
        return 0
    if len(argv) == 3 and argv[0] == "reading":
        got = reading_pass(argv[1], argv[2])
        if got is None:
            return 1
        print("\n".join([got[0], *got[1]]))
        return 0
    if len(argv) == 2 and argv[0] == "outcome":
        line = outcome_line(argv[1])
        if line is None:
            return 1
        print(line)
        return 0
    print("usage: bf_row.py row TRIES.md | passed TRIES.md | reading READING.json READING.txt | outcome READING.json",
          file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
