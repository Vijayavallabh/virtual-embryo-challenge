#!/usr/bin/env python
"""Row K28 (experiments/harness/TRIES.md, T3 candidate CHIPQ-C): for a query gene named to code, code picks three other
antigens of a pinned generic ChIP database (register row X018, as reopened for this row) whose binding near
transcription start sites overlaps the query's most specifically, from kept experiments only. This file is where the
row's steps live. For now it holds one: the first step of the database-only gate G0.

  g0a --commit C    G0a, STRUCTURE (the specification's file 02, steps 1 to 3, in experiments/t3_ext/chipq_partners.py;
                    file 04). For Ctnnb1 and for Gata4:
                      a1  2 or more usable experiments in the query's class;
                      a2  11 or more candidates in that class.
                    It prints the number of classes that class (vi) excludes and, per query, its usable experiments in
                    its class and its candidates: three lines of counts, nothing else.

G0 reads the pinned database files only: no fold, no barred factor, no panel, no score. It is uncounted. By the
specification the code of G0b, G0c and the read is written only after G0a has passed on the real database; those steps
are added here as further subcommands and change nothing that G0a computes.

A step refuses (exit 2, before any database file is opened) unless: the row's results folder holds no reading; row K28
is in TRIES.md; register row X018 in EXTERNAL_DATA.md carries the amendment's opening words (REGISTER_WORDS); the
tracked tree is clean; and HEAD is the commit given. Then every file of the fetch manifest is hashed again (a file off
its sha256 refuses, exit 2) and the committed kept list is loaded, before anything is counted.

Exit codes: 0 pass; 3 a bar not met; 2 refused. A bar that cannot be computed is not met. At exit 3 the step makes the
row's reading once, with the sentence S_G0, into experiments/t3_decomp/results/k28_chip_partners/ (reading.json and
reading.txt): nothing is counted, no fold is read, no file is named. A reading is never overwritten, and once one is
there the candidate is closed: every step refuses. A refusal is a fault to repair in code, never a reading.

No line this script prints, logs or writes carries an antigen name other than the two queries', a query's class, a
gene symbol, a database score or an experiment id: counts only. A crash shows its type, never its message.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import re
import socket
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
T3X = REPO / "experiments" / "t3_ext"
if str(T3X) not in sys.path:
    sys.path.insert(0, str(T3X))
import chipq_partners as cp  # noqa: E402  (steps 1 to 3 and the gate's counts)

cq, cf = cp.cq, cp.cf
Refused = cp.Refused

SCRIPT = "experiments/t3_decomp/k28_chip_partners.py"
ROW = "K28"
CANDIDATE = "CHIPQ-C"
TRIES_MD = REPO / "experiments" / "harness" / "TRIES.md"  # looked up at call time
REGISTER = REPO / "EXTERNAL_DATA.md"  # the disclosure register, looked up at call time
REGISTER_ROW = "X018"
REGISTER_WORDS = "Reopened by the captain's order of 2026-11-24"  # the amendment's opening words
RESULTS = HERE / "results" / "k28_chip_partners"  # the row's results folder, named as row K27's
READING_FILES = ("reading.json", "reading.txt")
QUERIES = cp.QUERIES  # ("Ctnnb1", "Gata4"): the test query, then the practice query
EXIT_PASS, EXIT_REFUSED, EXIT_BAR = 0, 2, 3

# THE SENTENCE THE READING PRINTS when a step of G0 ends at exit 3 (the specification's file 04, fixed)
S_G0 = "Try K28: STOP at G0: the database bars are not met; no fold read, no try, no file."

# What G0a prints: these three lines, counts only
LINE_CLASSES = "G0a: classes excluded by class (vi): {n}"
LINE_QUERY = "G0a {q}: usable experiments in its class: {u}; candidates: {c}"

DIRTY = ("the tree is not clean (git status --porcelain --untracked-files=no prints a tracked file, or git failed): "
         "not started, nothing read or written; commit or restore first: refused")
READING_EXISTS = ("the row's results folder holds a reading already: it is never overwritten, and the candidate is "
                  "closed (no step runs again): refused")


# ---------------------------------------------------------------------------------------------------------------------
# The ledger row, the register's wording, the tree and the commit


def _lines(path: Path) -> list[str]:
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise Refused(f"{path.name} is missing or a symlink: refused")
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def row_present(path: Path | None = None) -> str:
    """Row K28 is in the ledger: exactly one line of TRIES.md begins `| K28 ` (the row's id cell). Else Refused."""
    path = TRIES_MD if path is None else Path(path)
    n = sum(1 for ln in _lines(path) if ln.startswith(f"| {ROW} "))
    if n != 1:
        raise Refused(f"{'no' if n == 0 else 'more than one'} row {ROW} in {path.name}: the row is committed before "
                      "any step: refused")
    return ROW


def register_reopened(path: Path | None = None) -> None:
    """Register row X018 carries the amendment's opening words: exactly one line of EXTERNAL_DATA.md begins `| X018 `
    and it holds REGISTER_WORDS. Else Refused: no step opens a database file before the amendment is committed."""
    path = REGISTER if path is None else Path(path)
    rows = [ln for ln in _lines(path) if ln.startswith(f"| {REGISTER_ROW} ")]
    if len(rows) != 1:
        raise Refused(f"{'no' if not rows else 'more than one'} row {REGISTER_ROW} in {path.name}: refused")
    if REGISTER_WORDS not in rows[0]:
        raise Refused(f"row {REGISTER_ROW} of {path.name} does not carry the amendment's opening words "
                      f"({REGISTER_WORDS!r}): refused")


def _git(*args: str) -> tuple[int, str]:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["GIT_OPTIONAL_LOCKS"] = "0"
    try:
        r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, env=env, check=False)
    except OSError:
        return 1, ""
    return r.returncode, r.stdout.strip()


def stamp() -> dict:
    """HEAD, the clean flag (tracked files only, the harness's rule; False when git fails) and the machine."""
    rc_head, head = _git("rev-parse", "HEAD")
    rc_status, status = _git("status", "--porcelain", "--untracked-files=no")
    return {"commit": head if rc_head == 0 and head else None, "tree_clean": rc_status == 0 and status == "",
            "machine": socket.gethostname()}


def reading_absent(out: Path | None = None) -> None:
    """Refused when the row's results folder holds a reading (either of its two files, or a link under its name)."""
    out = RESULTS if out is None else Path(out)
    if any((out / f).exists() or (out / f).is_symlink() for f in READING_FILES):
        raise Refused(READING_EXISTS)


def preconditions(commit: str) -> dict:
    """Before any database file is opened: no reading yet, the row, the register's wording, a clean tracked tree, and
    HEAD at the commit given (its full hash). Returns the stamp."""
    reading_absent()
    row_present()
    register_reopened()
    st = stamp()
    if not st["tree_clean"]:
        raise Refused(DIRTY)
    if not re.fullmatch(r"[0-9a-f]{40}", str(commit)) or st["commit"] != commit:
        raise Refused("HEAD is not the commit given (its full hash): refused")
    return st


# ---------------------------------------------------------------------------------------------------------------------
# The reading (made once, at exit 3)


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reading_code(st: dict) -> dict:
    return {"script": SCRIPT, "script_sha256": _sha(__file__), "chipq_partners_sha256": _sha(T3X / "chipq_partners.py"),
            "chipq_sha256": _sha(T3X / "chipq.py"), "chipq_fetch_sha256": _sha(T3X / "chipq_fetch.py"),
            "commit": st["commit"], "tree_clean": st["tree_clean"], "machine": st["machine"]}


def g0a_lines(res: dict) -> list[str]:
    """G0a's printed block: the classes line, then one line per query in QUERIES' order. Counts only."""
    return [LINE_CLASSES.format(n=int(res["classes_excluded"])),
            *(LINE_QUERY.format(q=q, u=int(res["queries"][q]["usable_in_class"]),
                                c=int(res["queries"][q]["candidates"])) for q in QUERIES)]


def reading_record(step: str, counts: dict, st: dict, kept_sha256: str, fetch_manifest_sha256: str) -> dict:
    """The STOP reading's record: the fixed sentence, the step that stopped, its counts and bars, what was read (by
    sha256) and which code read it. Counts, the two queries' names and sha256 values only."""
    return {"outcome": S_G0, "stop": "G0", "step": step, "passes": False, "tries_counted": 0, "row": ROW,
            "candidate": CANDIDATE, "register_row": REGISTER_ROW, step: counts,
            "bars": {"a1_min_usable_in_class": cp.QUERY_MIN_USABLE, "a2_min_candidates": cp.MIN_CANDIDATES},
            "readings": cp.readings(), "kept_sha256": kept_sha256, "fetch_manifest_sha256": fetch_manifest_sha256,
            "reading_code": reading_code(st)}


def reading_text(rec: dict) -> str:
    """reading.txt: the outcome sentence on a line of its own beginning OUTCOME, then the counts."""
    lines = [f"OUTCOME {rec['outcome']}", f"tries counted: {rec['tries_counted']}",
             f"stopped at step {rec['step']} of gate G0 (database files only: no fold, no barred factor, no score)"]
    if rec.get("g0a") is not None:
        lines += g0a_lines(rec["g0a"])
        b = rec["bars"]
        lines.append(f"bars: a1, {b['a1_min_usable_in_class']} or more usable experiments in the query's class; a2, "
                     f"{b['a2_min_candidates']} or more candidates in that class")
        lines.append("bars met: " + "; ".join(
            f"{q} a1 {'yes' if v['a1'] else 'no'}, a2 {'yes' if v['a2'] else 'no'}"
            for q, v in ((q, rec["g0a"]["queries"][q]) for q in QUERIES)))
    lines.append(f"kept list {rec['kept_sha256']}; fetch manifest {rec['fetch_manifest_sha256']}")
    rc = rec["reading_code"]
    lines.append(f"read by {rc['script']} {rc['script_sha256'][:12]} (chipq_partners {rc['chipq_partners_sha256'][:12]}"
                 f", chipq {rc['chipq_sha256'][:12]}, chipq_fetch {rc['chipq_fetch_sha256'][:12]}) at commit "
                 f"{rc['commit']}, clean {rc['tree_clean']}, machine {rc['machine']}")
    return "\n".join(lines) + "\n"


def write_reading(rec: dict, out: Path | None = None) -> None:
    """The reading's two files, each created new (never over an existing file). Refused when a reading is there."""
    out = RESULTS if out is None else Path(out)
    reading_absent(out)
    bodies = {"reading.json": json.dumps(rec, indent=1) + "\n", "reading.txt": reading_text(rec)}
    out.mkdir(parents=True, exist_ok=True)
    for name in READING_FILES:
        try:
            with open(out / name, "x", encoding="utf-8") as f:
                f.write(bodies[name])
        except FileExistsError:
            raise Refused(READING_EXISTS) from None


# ---------------------------------------------------------------------------------------------------------------------
# G0a


def cmd_g0a(a) -> int:
    st = preconditions(a.commit)
    with contextlib.redirect_stdout(sys.stderr):  # anything a loader prints goes to stderr: stdout is the block's
        root = cf.root_of(None)
        if not cf.step_verify(root):
            raise Refused("a file of the pinned database is off its sha256 in MANIFEST.json, or the manifest is "
                          "empty: refused")
        kept = cq.load_kept()
        res = cp.g0a(root, kept)
        fetch_sha = cf.sha256_file(root / "MANIFEST.json")
    block = "".join(f"{ln}\n" for ln in g0a_lines(res))
    if not res["passes"]:
        write_reading(reading_record("g0a", res, st, kept.sha256, fetch_sha))
    sys.stdout.write(block)
    sys.stdout.flush()
    return EXIT_PASS if res["passes"] else EXIT_BAR


STEPS = {"g0a": cmd_g0a}  # G0b, G0c and the read are added here after G0a has passed on the real database


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("g0a")
    g.add_argument("--commit", required=True, help="the full hash of the commit the step runs at (HEAD)")
    a = ap.parse_args(sys.argv[1:] if argv is None else list(argv))
    try:
        return STEPS[a.cmd](a)
    except Refused as e:
        print(f"refused: {e}", file=sys.stderr, flush=True)
        return EXIT_REFUSED
    except Exception as e:  # noqa: BLE001  (a message could carry a name from a file: the type only)
        print(f"refused: {type(e).__name__} in step {a.cmd}", file=sys.stderr, flush=True)
        return EXIT_REFUSED


if __name__ == "__main__":
    sys.exit(main())
