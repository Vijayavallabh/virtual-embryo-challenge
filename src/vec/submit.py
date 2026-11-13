"""Submission client (feat-016): free format checks by default, guarded official uploads, one ledger row per attempt.

The challenge API treats an upload as OFFICIAL unless the request sends format_only=true, and the final phase allows 2
scored submissions per board for the whole phase, which cannot be withdrawn. So this module:

- sends format_only=true unless `official` is set AND `confirm` equals the first 12 hex characters of the file's sha256;
- refuses an official upload unless every guard passes: the local format check, a board whose upload fields and split
  are known, a clean working tree (vec.provenance's clean-tree rule: no uncommitted change to a tracked file and no
  untracked, non-ignored file outside results/, runs/ and experiments/**/results/; the ledger and the file's own build
  stamp are exempt), a build stamp that ties the file to a clean build at the local HEAD or an ancestor of it (below),
  at least one known disclosure id from EXTERNAL_DATA.md, a passing server format check of the same file for the same
  board in the ledger, a phase that takes submissions for the board's split, a readable quota and submission list, and
  room in the board's budget (below); it also refuses a file the board already holds through a format check the server
  stored, unless an official attempt of that file is still unanswered (its identical retry cannot take a second slot);
- refuses an official upload while the captain's upload gate stands (guard "upload gate"; 2026-11-05, his words: "Until
  our estimated score in leaderboard crosses 209, don't submit our 1st submission"; docs/PLAN.md section 8): the
  committed record UPLOAD_GATE must hold an estimate above UPLOAD_GATE_MIN with a named independent review, or the
  captain's lift, or its deadline exception must have opened (2026-11-27 00:00 UTC; the captain: "I am fine with
  submitting our 1st submission at Nov 27 anyway"); an unreadable or inconsistent record refuses too. Format
  checks never read it;
- refuses an official upload of a file that is not in the current reference list (guard "reference list"; AGENTS.md:
  an official upload sends only a file whose sha256 is in it), so that a superseded file left in an older archive
  folder with a clean build stamp cannot go: the committed pointer (results/reference_list.json of the repository
  the client runs in, rewritten after each rehearsal) names the list (sha256sum output, a path inside the repository),
  its sha256 and, per board, the list's names meant for that board (and, as "held", names no board may send, with the
  reason). The pointer must be a regular file and parse (at most REFERENCE_LIST_MAX_BYTES), the list must lie inside
  the repository, exist and have that sha256, the board must have names, none held, each of which occurs in the list
  exactly once, and the plan's sha256 (for a pulled file, its copy's) must be the list's hash for one of the board's
  names. A format check never refuses on it, and warns when an official upload of the file would fail it;
- ties an official upload to the code that built the file (guard "build stamp"): the stamp `<file>.stamp.json` next to
  the file's real path (written with the file by vec.io.write_submission; see vec.provenance) must exist and parse (at
  most STAMP_MAX_BYTES), its sha256 must be the plan's (the bytes that will be sent), it must say the build was clean
  (the tree passed the clean-tree rule when the build started, at the builder's first import of vec.io, and when the
  file was written, at one HEAD that did not move in between), and its commit must be a full hex commit that is the
  local HEAD or an ancestor of it (`git merge-base --is-ancestor`, run in the local repo). So a file built before a
  later fix is recorded under the commit that built it, and a file whose build saw a pull, checkout, sync or stash is
  refused and its build commit marked +dirty in the ledger (an uncommitted edit made and undone while a build runs is
  the one change git cannot see). A format check records a matching stamp's commit and warns when the stamp is missing
  or would fail the guard; it never requires one;
- refuses a format check too when the phase (read without the token) cannot be read or its split is not the board's:
  the server routes T1 and T2 uploads to the phase's split by itself, so a format check it stored would land elsewhere;
- refuses every POST, format check or official, to every board while a format check that came back as a submission
  (UNEXPECTED_SUBMISSION: the server ignored format_only, which would not be limited to one board) or unconfirmed
  (format_unconfirmed) waits for a person's `vec submit resolve`. A format check answered with a submission id is
  recorded, counted and blocked under the board the server named (the answer's bundle_key) and, if different, under
  the board the request named; each record is resolved on its own board;
- returns the earlier record instead of uploading when the server lists the same file on the board, or when the ledger
  holds the server's answer to an official upload of it, but never a format check's submission (made without the key:
  the identical request goes out under the key instead); a replay row (mode replay) never counts or settles anything;
- takes a file built on the DGX (`from_dgx`, `vec submit PATH --from-dgx`) by pulling it for the call only: PATH is
  relative to the DGX repo root (DGX_REPO on DGX_HOST) or absolute inside it, and one that resolves outside it (through
  .., a symlink or an absolute path) or holds a newline or NUL is refused. One remote call reports the resolved path,
  its sha256 and size, the build stamp next to that real path (its sha256 and size; only its size when it is larger than
  STAMP_MAX_BYTES; or that there is none, that it is a link or not a regular file, or that it cannot be read), the DGX
  repo's HEAD and its `git status --porcelain -z`; `cat` then streams the file into a 0600 file in a 0700 temporary
  directory while it is hashed, then a stamp reported at no more than STAMP_MAX_BYTES into a 0600 file next to it
  (stopping one byte past that size), and the pull is refused ("transfer mismatch") unless each copy has the reported
  sha256 and size. The copy then goes through everything here as a local file would, its stamp with it, and the
  directory is removed when the call ends, also on errors, Ctrl-C, SIGHUP (a closed terminal or dropped session) and
  SIGTERM, which unwind the call; a copy left by a call killed outright (SIGKILL, power loss) is removed by the next
  call once it is STALE_COPY_HOURS old. A refused pull, like a missing local file, writes no ledger row. An official
  upload also needs the guard "dgx provenance": the DGX's HEAD is the local HEAD and its tree passes the clean-tree rule
  above (the same parser reads its status; the ledger and the file's stamp are exempt); a format check records the DGX
  state without requiring it. Every remote command goes through dgx_run (`ssh -o BatchMode=yes DGX_HOST -- COMMAND`,
  paths quoted with shlex.quote), whose environment leaves out VEC_TOKEN and every variable whose name contains TOKEN in
  any case;
- sends a private snapshot: the file is opened before the snapshot exists, then copied once into a 0600 file in a 0700
  temporary directory and hashed from the copy (failing to open or read the file is "cannot read FILE"; only failing to
  write the copy points at the temporary directory), nothing is sent unless that hash equals the plan's ("file changed
  before upload"), every attempt sends the snapshot's bytes with its Content-Length, and the snapshot is deleted
  afterwards;
- keeps the filename, form fields and Idempotency-Key (derived from the board and the file hash) identical on every
  retry of an official upload, because the server answers a changed payload under a used key with 409;
- appends a row to results/submissions.csv for every attempt, including refusals, replays and resolutions, and writes a
  "sending" row before each official request so a crash mid-upload still leaves a trace. Its git_commit column is the
  local HEAD at the call (plus "+dirty" when the tree is not clean). Its last two columns are origin, "local" for a
  local file and "PrakashDGX_H2:<path>@<DGX HEAD>" (plus "+dirty" when the DGX tree is not clean, or "+unknown" when
  its git status cannot be read) for a pulled one, whose file column holds that repo-relative DGX path (symlinks
  resolved) and whose upload filename is the basename of PATH as given, as for a local file; and build_commit, the
  commit in the file's build stamp when the stamp matches the file (plus "+dirty" when it says the build was not
  clean), else empty.

The budget of a board is the larger of the quota endpoint's count (the source of truth) and the slots below, against
the quota's limit (at most 2); if the quota or the submission list cannot be read or parsed, official uploads are
refused. Every server entry on the board that is not rejected takes one slot (counted by id). The ledger gives each file
slots of its own: one per server id from an answer with a status other than rejected or from a format check that came
back as a submission, one for an ambiguous official outcome (no response, HTTP 3xx or 5xx, a 2xx without both an id and
a status, or a "sending" row with no outcome after it), and one per format check answered format_unconfirmed. Each
server entry of the file, rejected or not, settles one of these: the slot of its id if the ledger has that id, else
(matched by sha256) one without an id; the others count on top of the server's entries. An ambiguous file keeps
counting, however old and even after a later 4xx, until a definite answer (a 2xx with id and status) or a server
listing settles it, or a person who has checked the portal records `vec submit resolve SHA12 --board KEY --as
not-stored --note TEXT` (a later listing still counts it). The other resolution, `--as acknowledged`, is for a format
check that came back as a submission: it keeps counting, and POSTs resume once no such check waits on any board. The
plan's own unanswered official attempts are left out of its room check, since an identical retry under the same key
cannot take a second slot.

The token is loaded with vec.data.load_token at call time (it refuses a token with whitespace, control or non-ASCII
characters), sent only to the API base (no redirects are followed; never to the DGX, whose commands, arguments and
environment never hold it), and never printed or logged: _call redacts it from everything the server sends back (the
parsed JSON's keys and values at any depth, and the raw error text before it is cut) and from error messages, so ids,
statuses, guard details and `vec submit status` only ever see "<token>".

Endpoints and fields, from public material only (unauthenticated GETs, 2026-10-19): the API's OpenAPI document
(https://api.virtualembryo.ai/openapi.json) and the portal's API guide and upload code in its public JavaScript bundle
(https://virtualembryo.ai/_nuxt/CFNGMlCM.js, component ChallengeApiGuide; upload function in _nuxt/DxvW8jU0.js).
- POST /challenge/submissions, multipart: task (T1|T2|T3), setting and mode (T2), target (T3; the portal derives it from
  the board key and the OpenAPI form lists it), optional model and architecture, format_only (default false), file;
  optional Idempotency-Key header (1-128 of A-Za-z0-9._:-). A format check answers {ok, problems, format_only} and
  creates nothing; an upload answers {id, status, task, bundle_key, phase, bytes, sha256, problems, used_credit}, and a
  rejected one can still come back with HTTP 200.
- GET /challenge/integration/quota: daily counts and, in the final phase, `official_per_board`, a phase-wide allowance
  per board. The shape of its entries is not published; this module reads {board: {"used": n, "limit": m}} (or
  "remaining" in place of "used", or both when they agree; any other field only with a string or null value) and
  treats anything else as unknown, which refuses official uploads.
- GET /challenge/submissions (latest 200, each with bundle_key, status and the server's sha256), and
  GET /challenge/submissions/{id}.
- GET /challenge/phase (public): split and accepts_submissions.
"""

from __future__ import annotations

import contextlib
import csv
import fcntl
import hashlib
import http.client
import json
import os
import re
import shlex
import shutil
import signal
import stat
import subprocess
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import BinaryIO

from vec import data as vec_data
from vec.boards import BOARDS_DIR, REPO_ROOT, Board, get_board, load_boards
from vec.io import CheckReport, check_submission
from vec.provenance import (
    LEDGER_PATH,
    STAMP_MAX_BYTES,
    STAMP_SUFFIX,
    _repo_paths,
    git_state,
    is_ancestor,
    read_stamp,
    stamp_path,
    unclean_paths,
)

# the clean-tree rule moved to vec.provenance: the names it had here, used here or not, stay importable from here
from vec.provenance import OUTPUT_DIRS as OUTPUT_DIRS
from vec.provenance import _output_path as _output_path

API_BASE = vec_data.API_BASE
LEDGER = REPO_ROOT / LEDGER_PATH
EXTERNAL_DATA = REPO_ROOT / "EXTERNAL_DATA.md"
UPLOAD_GATE = REPO_ROOT / "results" / "upload_gate.json"  # the captain's upload gate record (docs/PLAN.md section 8)
UPLOAD_GATE_MIN = 209.0  # the captain's figure (2026-11-05); only he changes it
UPLOAD_GATE_MAX_BYTES = 64 * 1024
REFERENCE_POINTER = "results/reference_list.json"  # the pointer, relative to the repository root
REFERENCE_LIST = REPO_ROOT / REFERENCE_POINTER  # the pointer read without a repo (guard below)
REFERENCE_LIST_MAX_BYTES = 64 * 1024  # the pointer file
REFERENCE_SUMS_MAX_BYTES = 1024 * 1024  # the list file it names (sha256sum output)
SUMS_LINE = re.compile(r"^([0-9a-f]{64}) [ *](.+)$")  # one line of sha256sum output: hash, mode character, name
FINAL_PHASE_LIMIT = 2  # scored submissions per board for the whole final phase
LEDGER_FIELDS = ["utc_time", "board", "mode", "file", "upload_name", "sha256", "server_status", "server_id",
                 "http_status", "git_commit", "disclosure", "idempotency_key", "fields", "note", "detail", "origin",
                 "build_commit"]
REJECTED = "rejected"  # the only server status that spends no official attempt (cancelled ones are not refunded)
CLIENT_STATUSES = frozenset({"refused", "sending", "no_response", "http_error", "unconfirmed", "unknown", "format_ok",
                             "format_failed", "format_unconfirmed", "UNEXPECTED_SUBMISSION"})  # never a server status
FORMAT_SUBMISSION = ("UNEXPECTED_SUBMISSION", "format_unconfirmed")  # format checks that may have made a submission
RESOLUTIONS = ("not-stored", "acknowledged")  # vec submit resolve --as
QUOTA_FIELDS = ("used", "limit", "remaining")
RETRY_DELAYS = (2.0, 4.0, 8.0)  # the API guide's back-off for timeouts and 5xx; 4xx are never retried automatically
GET_TIMEOUT = 60.0
UPLOAD_TIMEOUT = 900.0
IDEMPOTENCY_KEY = re.compile(r"^[A-Za-z0-9._:-]{1,128}$")
DISCLOSURE_ID = re.compile(r"^X\d{3}$")
LOOPBACK = ("http://127.0.0.1:", "http://localhost:")
DGX_HOST = "PrakashDGX_H2"  # where predictions are built (vec submit PATH --from-dgx); the token never goes there
DGX_REPO = "~/v/virtual-embryo-challenge"  # the repo's clone on DGX_HOST: a pulled PATH must resolve inside it
STALE_COPY_HOURS = 6.0  # a private copy this old is left over from a call that was killed outright (SIGKILL, power)
COPY_PREFIXES = ("vec-dgx-", "vec-submit-")  # the temporary directories of pulled copies and of snapshots
# One call on the DGX: refuse a path ($2) that is missing, not a regular file or outside the repo root ($1), then print
# the root, the file's real path, sha256, size, the repo's HEAD and the build stamp next to the real path (the real path
# plus $3): "<sha256> <size>", "large <size>" above $4 bytes (never hashed), "none", "special" (a link or not a regular
# file: never read) or "unreadable"; NUL-separated, and, last, the repo's git status. The path is an argument and never
# part of the script; `printf x` keeps a real path's trailing newline for the check here.
DGX_REPORT = r"""cd -- "$1" && root=$(pwd -P) || exit 3
f=$(realpath -e -- "$2" && printf x) || exit 4
f=${f%?x}
[ -f "$f" ] || { echo "not a regular file: $f" >&2; exit 4; }
case $f in "$root"/*) ;; *) echo "resolves to $f, outside the repository $root" >&2; exit 4 ;; esac
sum=$(sha256sum < "$f") && size=$(stat -c %s -- "$f") || exit 5
s=$f$3
if [ -L "$s" ] || { [ -e "$s" ] && [ ! -f "$s" ]; }; then stamp=special
elif [ ! -e "$s" ]; then stamp=none
elif ! ssize=$(stat -c %s -- "$s"); then stamp=unreadable
elif [ "$ssize" -gt "$4" ]; then stamp="large $ssize"
elif ssum=$(sha256sum < "$s"); then stamp="${ssum%% *} $ssize"
else stamp=unreadable; fi
head=$(git rev-parse --verify -q HEAD) || head=unknown
printf 'vec-dgx\0%s\0%s\0%s\0%s\0%s\0%s\0' "$root" "$f" "${sum%% *}" "$size" "$head" "$stamp"
GIT_OPTIONAL_LOCKS=0 exec git status --porcelain -z --untracked-files=all"""


class SubmitRefused(RuntimeError):
    """A guard refused the request; nothing was sent."""

    def __init__(self, message: str, plan: Plan | None = None):
        super().__init__(message)
        self.plan = plan


@dataclass
class Guard:
    name: str
    ok: bool
    detail: str = ""


@dataclass
class Plan:
    """What one attempt will send, and the guard results; built and shown before anything is sent."""

    file: Path
    board: str
    official: bool
    sha256: str
    upload_name: str
    fields: dict[str, str]
    idempotency_key: str
    disclosure: list[str]
    note: str
    git_commit: str
    guards: list[Guard] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    replay: dict | None = None  # earlier record when this file was already sent officially to this board
    remote: str = ""  # the file's path in the DGX repo when it was pulled from there (the ledger's file column)
    origin: str = "local"  # the ledger's origin column: local, or DGX_HOST:<remote>@<DGX HEAD>[+dirty|+unknown]
    build_commit: str = ""  # the ledger's build_commit column: a matching build stamp's commit[+dirty], else ""

    @property
    def ok(self) -> bool:
        return all(g.ok for g in self.guards)

    @property
    def failed(self) -> list[str]:
        return [g.name for g in self.guards if not g.ok]


@dataclass
class Result:
    plan: Plan
    server_status: str  # format_ok | format_failed | format_unconfirmed | UNEXPECTED_SUBMISSION | a server status
    # (queued, scored, rejected, ...) | unconfirmed (a 2xx without both an id and a status) | http_error | no_response
    server_id: str = ""
    http_status: int | None = None
    detail: str = ""
    replayed: bool = False
    boards: tuple[str, ...] = ()  # an UNEXPECTED_SUBMISSION's ledger boards: the server's bundle_key, then the plan's

    @property
    def ok(self) -> bool:
        if self.plan.official:
            return self.server_status not in (REJECTED, "cancelled", "http_error", "no_response", "unconfirmed",
                                              "unknown")
        return self.server_status == "format_ok"


@dataclass
class Budget:
    """Official uploads already spent on one board: the server's quota and entries, and files only the ledger knows."""

    board: str
    server_used: int | None  # from the quota endpoint (source of truth); None if unknown
    server_limit: int | None
    listed: list[str]  # ids of the server's entries on the board that are not rejected: one slot each
    unlisted: dict[str, int]  # sha256 -> slots of a file the server does not list (yet), from the ledger
    ambiguous: set[str]  # files whose outcome is unknown: they count until listed or resolved by a person
    retryable: set[str]  # files with unanswered official attempts (one slot each; an identical retry is free)
    unresolved: list[dict] = field(default_factory=list)  # format checks that came back as a submission: no POST
    reason: str = ""  # why the server budget is unknown

    @property
    def limit(self) -> int:
        return FINAL_PHASE_LIMIT if self.server_limit is None else min(self.server_limit, FINAL_PHASE_LIMIT)

    def spent_besides(self, sha: str | None = None) -> int:
        """Slots spent, leaving out the slot of `sha`'s unanswered official attempts: an identical retry under the same
        Idempotency-Key cannot take a second one."""
        return max(self.server_used or 0, len(self.listed) + sum(self.unlisted.values()) - (sha in self.retryable))

    @property
    def spent(self) -> int:
        return self.spent_besides()


@dataclass
class Reply:
    status: int  # HTTP status; 0 = no response (timeout or connection failure)
    body: dict | None
    detail: str


# ---------------------------------------------------------------------------------------------------------------------
# Files, boards, disclosures


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def board_fields(board: Board) -> dict[str, str]:
    """The form fields that name a board, derived from its key exactly as the portal does.

    T1 -> task; T2 -> task, setting, mode (the part of the split after '_'); T3 -> task, target. The server picks the
    phase's validation or test split by itself, so the board's split is checked against the phase before an official
    upload.
    """
    parts = board.key.split(":")
    task = parts[0]
    if task == "T1" and len(parts) == 2:
        out = {"task": task}
    elif task == "T2" and len(parts) == 3 and "_" in parts[2]:
        out = {"task": task, "setting": parts[1], "mode": parts[2].split("_", 1)[1]}
        if out["setting"] not in ("heart", "embryo") or out["mode"] not in ("interp", "extrap"):
            raise ValueError(f"unexpected T2 setting or mode in board key {board.key!r}")
    elif task == "T3" and len(parts) == 2 and re.fullmatch(r"[A-Za-z0-9._-]+", parts[1]):
        out = {"task": task, "target": parts[1]}
    else:
        raise ValueError(f"cannot derive upload fields from board key {board.key!r}")
    if board.task != task:
        raise ValueError(f"board {board.key!r} has task {board.task!r}")
    if board.setting and out.get("setting", board.setting) != board.setting:
        raise ValueError(f"board {board.key!r} has setting {board.setting!r}")
    if board.target_ko and out.get("target", board.target_ko) != board.target_ko:
        raise ValueError(f"board {board.key!r} has target {board.target_ko!r}")
    return out


def subtask_key(board: Board) -> str:
    """The portal's per-subtask quota key: T2/<setting>/<mode> for T2, the task otherwise."""
    f = board_fields(board)
    return f"T2/{f['setting']}/{f['mode']}" if f["task"] == "T2" else f["task"]


def board_split(board: Board) -> str | None:
    """'val' or 'test', from the board's published split (or its key); None if neither."""
    split = board.split or (board.key.split(":")[1] if board.task == "T1" else board.key.split(":")[-1])
    split = split.lower()
    return "test" if split.startswith("test") else "val" if split.startswith("val") else None


def idempotency_key(board: str, sha256: str) -> str:
    key = f"vec.{re.sub(r'[^A-Za-z0-9.:-]', '-', board)}.{sha256}"
    if not IDEMPOTENCY_KEY.match(key):
        raise ValueError(f"cannot build an Idempotency-Key for board {board!r}")
    return key


def upload_name(path: Path) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", path.name)


def check_local(file: str | Path, board: str | Board, boards_dir: Path | None = None) -> CheckReport:
    """Free local format check (vec.io.check_submission plus the portal's file checks); refuses on any error."""
    b = board if isinstance(board, Board) else get_board(board, boards_dir or BOARDS_DIR)
    path = Path(file)
    if not path.name.lower().endswith(".h5ad"):
        raise SubmitRefused("local format check: the file must be an AnnData .h5ad")
    if path.stat().st_size == 0:
        raise SubmitRefused("local format check: the file is empty")
    try:
        rep = check_submission(path, b)
    except Exception as e:  # unreadable or malformed file
        raise SubmitRefused(f"local format check: cannot read {path.name}: {type(e).__name__}: {e}") from e
    if not rep.ok:
        raise SubmitRefused("local format check: " + "; ".join(rep.errors))
    return rep


def disclosure_register(path: Path | None = None) -> dict[str, str]:
    """Disclosure ids in EXTERNAL_DATA.md's table, each with its source cell."""
    try:
        text = Path(path or EXTERNAL_DATA).read_text()
    except OSError:
        return {}
    reg = {}
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) > 1 and DISCLOSURE_ID.match(cells[0]):
            reg[cells[0]] = cells[1]
    return reg


def parse_ids(disclose: str | Sequence[str]) -> list[str]:
    items = disclose.split(",") if isinstance(disclose, str) else list(disclose)
    return list(dict.fromkeys(i.strip().upper() for i in items if i.strip()))


# ---------------------------------------------------------------------------------------------------------------------
# Ledger


class Ledger:
    """results/submissions.csv, locked for the whole call so that two clients cannot both pass the budget guard."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._fh = None
        self._rows: list[dict] = []

    def __enter__(self) -> Ledger:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fh = self.path.open("a+", newline="")
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            fh.close()
            raise SubmitRefused(f"{self.path} is locked by another vec submit; try again when it finishes") from None
        fh.seek(0)
        reader = csv.DictReader(fh)
        rows = list(reader)
        if reader.fieldnames is not None and list(reader.fieldnames) != LEDGER_FIELDS:
            fh.close()
            raise SubmitRefused(f"{self.path} has unexpected columns; refusing to append to it")
        self._fh, self._rows = fh, rows
        return self

    def __exit__(self, *exc) -> None:
        if self._fh is not None:
            fcntl.flock(self._fh, fcntl.LOCK_UN)
            self._fh.close()
            self._fh = None

    def rows(self) -> list[dict]:
        return list(self._rows)

    def append(self, row: dict) -> None:
        assert self._fh is not None, "ledger is not open"
        self._fh.seek(0, os.SEEK_END)
        w = csv.DictWriter(self._fh, fieldnames=LEDGER_FIELDS)
        if self._fh.tell() == 0:
            w.writeheader()
        w.writerow(row)
        self._fh.flush()
        os.fsync(self._fh.fileno())
        self._rows.append(row)


def read_ledger(path: Path | None = None) -> list[dict]:
    path = Path(path or LEDGER)
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _row(plan: Plan, status: str, server_id: str = "", http_status: int | None = None, detail: str = "") -> dict:
    return {"utc_time": _now(), "board": plan.board,
            "mode": "official" if plan.official else "format", "file": plan.remote or _display_path(plan.file),
            "upload_name": plan.upload_name, "sha256": plan.sha256, "server_status": status, "server_id": server_id,
            "http_status": "" if http_status is None else str(http_status), "git_commit": plan.git_commit,
            "disclosure": ";".join(plan.disclosure), "idempotency_key": plan.idempotency_key,
            "fields": json.dumps(plan.fields), "note": plan.note, "detail": detail, "origin": plan.origin,
            "build_commit": plan.build_commit}


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


@dataclass
class _File:
    """What the ledger says about one file (sha256) on one board, read from its rows in order."""

    answer: dict | None = None  # latest answer to an official upload: a row with the server's submission id and status
    ambiguous: bool = False  # an official attempt may have stored it, and no answer or resolution has settled it since
    sending: bool = False  # a "sending" row still waiting for its outcome row (a crash mid-upload leaves one)
    format_ids: dict[str, dict] = field(default_factory=dict)  # format checks that came back as a submission, by id
    unconfirmed: int = 0  # format checks answered format_unconfirmed since the last not-stored resolution
    unresolved: list[dict] = field(default_factory=list)  # rows of those two kinds that no person has resolved yet

    def ids(self) -> set[str]:
        return set(self.format_ids) | ({self.answer["server_id"]} if self.answer else set())

    def held(self) -> set[str]:
        """The server ids that take a slot: format checks that came back as a submission, and the answer unless it was
        rejected. The file's other slots have no id: one if ambiguous, and one per format_unconfirmed check."""
        ids = set(self.format_ids)
        if self.answer is not None:
            ids.discard(self.answer["server_id"])
            if self.answer["server_status"] != REJECTED:
                ids.add(self.answer["server_id"])
        return ids


def _attempt(row: dict) -> str:
    """What one official outcome row says: "answer" (the server's submission id and status), "ambiguous" (the upload may
    have been stored), "not_stored" (HTTP 4xx) or "" (not an outcome row this client writes)."""
    status, code = row.get("server_status", ""), row.get("http_status", "")
    if status == "http_error":
        return "not_stored" if code.startswith("4") else "ambiguous"  # 3xx and 5xx may have been processed
    if row.get("server_id") and status and status not in CLIENT_STATUSES:
        return "answer"
    return "ambiguous" if status in ("no_response", "unconfirmed") or code else ""


def _ledger_files(rows: list[dict], board: str) -> dict[str, _File]:
    """What the ledger says about each file sent to `board` (the rules are in the module docstring). Replay rows sent
    nothing and are ignored."""
    files: dict[str, _File] = {}
    for r in rows:
        if r.get("board") != board or not r.get("sha256"):
            continue
        f, mode, status = files.setdefault(r["sha256"], _File()), r.get("mode"), r.get("server_status", "")
        if mode == "official" and status == "sending":
            f.ambiguous |= f.sending and f.answer is None  # the attempt before this one never got an outcome row
            f.sending = True
        elif mode == "official" and status != "refused":
            f.sending, kind = False, _attempt(r)
            if kind == "answer":
                f.answer, f.ambiguous = r, False
            elif kind == "ambiguous" and f.answer is None:
                f.ambiguous = True  # a later 4xx does not undo this
        elif mode == "format" and status in FORMAT_SUBMISSION:
            if status == "UNEXPECTED_SUBMISSION":
                f.format_ids[r.get("server_id") or r.get("utc_time", "")] = r
            else:
                f.unconfirmed += 1
            f.unresolved.append(r)
        elif mode == "resolve" and status == "not-stored":
            f.ambiguous = f.sending = False
            f.unconfirmed = 0
            f.unresolved = [x for x in f.unresolved if x.get("server_status") == "UNEXPECTED_SUBMISSION"]
        elif mode == "resolve" and status == "acknowledged":
            f.unresolved = []
    for f in files.values():
        f.ambiguous |= f.sending and f.answer is None
    return files


def _unresolved(rows: list[dict], board: str | None = None) -> list[dict]:
    """Format-check rows that came back as a submission or unconfirmed and await a person (ledger order): on `board`,
    or on every board."""
    boards = [board] if board is not None else list(dict.fromkeys(r.get("board", "") for r in rows))
    pending = {id(r) for b in boards for f in _ledger_files(rows, b).values() for r in f.unresolved}
    return [r for r in rows if id(r) in pending]


def resolve_commands(rows: list[dict]) -> dict[tuple[str, str], tuple[str, str]]:
    """For each (board, sha256) among format-check rows that await a person: what came back, and the `vec submit
    resolve` command (--as acknowledged once the server gave a submission id; not-stored can only settle the rest)."""
    kinds: dict[tuple[str, str], list[str]] = {}
    for r in rows:
        k = kinds.setdefault((r["board"], r["sha256"]), [])
        k += [r["server_status"]] if r["server_status"] not in k else []
    return {(board, sha): (" and ".join(k), f"vec submit resolve {sha[:12]} --board {shlex.quote(board)} --as "
                           f"{'acknowledged' if 'UNEXPECTED_SUBMISSION' in k else 'acknowledged (or not-stored)'} "
                           "--note TEXT") for (board, sha), k in kinds.items()}


# ---------------------------------------------------------------------------------------------------------------------
# Server state


def official_quota(quota: dict | None, board: Board) -> tuple[int, int] | str:
    """(used, limit) for `board` from the quota endpoint's official_per_board, or why it is unknown."""
    if not isinstance(quota, dict):
        return "the quota endpoint could not be read"
    per_board = quota.get("official_per_board")
    if per_board is None:
        return "the quota response has no official_per_board (not the final phase?)"
    if not isinstance(per_board, dict):
        return f"official_per_board has an unexpected shape ({type(per_board).__name__})"
    keys = list(dict.fromkeys([board.key, subtask_key(board)]))
    hits = [k for k in keys if k in per_board]
    if len(hits) != 1:
        return (f"official_per_board has no single entry for {' or '.join(keys)} "
                f"(keys: {', '.join(map(str, per_board))})")
    entry = per_board[hits[0]]
    # any other field that is not a string may be a count this parser does not understand
    if isinstance(entry, dict) and all(k in QUOTA_FIELDS or v is None or isinstance(v, str) for k, v in entry.items()):
        used, limit, remaining = entry.get("used"), entry.get("limit"), entry.get("remaining")
        if _count(used) and _count(limit) and (remaining is None or (_count(remaining) and used + remaining == limit)):
            return used, limit
        if used is None and _count(remaining) and _count(limit) and remaining <= limit:
            return limit - remaining, limit
    return f"official_per_board[{hits[0]!r}] has an unexpected shape: {json.dumps(entry)[:160]}"


def _count(x) -> bool:
    return isinstance(x, int) and not isinstance(x, bool) and x >= 0


def server_entries(subs: list[dict], board: str, ident: dict[str, str] | None = None) -> list[dict]:
    """Server submissions on `board`: by bundle_key, or (entries without one) by task and the board's fields."""
    out = []
    for e in subs:
        if not isinstance(e, dict):
            continue
        bundle = e.get("bundle_key")
        if bundle is not None:
            if bundle == board:
                out.append(e)
        elif ident and e.get("task") == ident.get("task") and all(
                e.get(k) in (None, v) for k, v in ident.items() if k != "task"):
            out.append(e)  # conservative: may over-count, never under-count
    return out


def board_budget(board: str, rows: list[dict], quota: tuple[int, int] | str, subs: list[dict] | None,
                 ident: dict[str, str] | None = None) -> Budget:
    """The board's budget: server entries by id, plus the ledger's slots that no server entry has settled."""
    entries = server_entries(subs or [], board, ident)
    listed = list(dict.fromkeys(str(e.get("id") or f"#{i}") for i, e in enumerate(entries)
                                if e.get("status") != REJECTED))
    unlisted, ambiguous, retryable, unresolved = {}, set(), set(), []
    for sha, f in _ledger_files(rows, board).items():
        unresolved += f.unresolved
        ids, held, anonymous = f.ids(), f.held(), int(f.ambiguous) + f.unconfirmed
        for e in entries:  # each entry of the file settles one slot: its id's, else (same sha256) one without an id
            if (eid := str(e.get("id") or "")) and eid in ids:
                held.discard(eid)
            elif e.get("sha256") == sha:
                anonymous = max(anonymous - 1, 0)
        if n := len(held) + anonymous:
            unlisted[sha] = n
        if anonymous:
            ambiguous.add(sha)
        if anonymous and f.ambiguous:
            retryable.add(sha)
    used, limit, reason = (quota[0], quota[1], "") if isinstance(quota, tuple) else (None, None, quota)
    return Budget(board, used, limit, listed, unlisted, ambiguous, retryable, unresolved, reason)


def _subs_list(reply: Reply) -> list[dict] | None:
    if reply.status == 200 and isinstance(reply.body, dict) and isinstance(reply.body.get("submissions"), list):
        return reply.body["submissions"]
    return None


# ---------------------------------------------------------------------------------------------------------------------
# HTTP


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: D102
        return None  # never follow a redirect: the token goes to the API base and nowhere else


def _opener(url: str) -> urllib.request.OpenerDirector:
    if url.startswith(LOOPBACK):  # local test servers: never through a proxy
        return urllib.request.build_opener(_NoRedirect, urllib.request.ProxyHandler({}))
    return urllib.request.build_opener(_NoRedirect)


def _check_base(base: str) -> str:
    base = base.rstrip("/")
    if not (base.startswith("https://") or base.startswith(LOOPBACK)):
        raise SubmitRefused(f"refusing to send the token to {base!r}: https (or a loopback test server) only")
    return base


def _call(req: urllib.request.Request, timeout: float, token: str | None) -> Reply:
    """One request. The token is redacted here, once, from everything that comes back (the parsed JSON's keys and
    values at any depth, the raw error text before it is cut) and from error messages, so no caller ever sees it."""
    try:
        with _opener(req.full_url).open(req, timeout=timeout) as r:
            status, raw = r.status, r.read()
    except urllib.error.HTTPError as e:
        status, raw = e.code, e.read()
    except (OSError, ValueError, http.client.HTTPException) as e:  # ValueError: a bad header value, a garbled answer
        return Reply(0, None, vec_data.redact(f"no response ({type(e).__name__}: {getattr(e, 'reason', e)})", token))
    try:
        parsed = vec_data.redact(json.loads(raw.decode()), token) if raw else None
    except (UnicodeDecodeError, json.JSONDecodeError):
        parsed = None
    body = parsed if isinstance(parsed, dict) else None
    detail = ""
    if not 200 <= status < 300:
        d = body.get("detail") if body is not None else parsed if parsed is not None else \
            vec_data.redact(raw.decode(errors="replace"), token)
        detail = f"HTTP {status}: {(d if isinstance(d, str) else json.dumps(d))[:300]}"
    return Reply(status, body, detail)


def _get(base: str, path: str, token: str | None = None) -> Reply:
    headers = {"Accept": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return _call(urllib.request.Request(f"{base}{path}", headers=headers, method="GET"), GET_TIMEOUT, token)


class _Multipart:
    """multipart/form-data body streamed from the snapshot. The boundary comes from the file hash, so retries are
    identical."""

    def __init__(self, fields: dict[str, str], path: Path, size: int, name: str, boundary: str):
        parts = [f'--{boundary}\r\nContent-Disposition: form-data; name="{k}"\r\n\r\n{v}\r\n'
                 for k, v in fields.items()]
        parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="file"; filename="{name}"\r\n'
                     "Content-Type: application/octet-stream\r\n\r\n")
        self.head = "".join(parts).encode()
        self.tail = f"\r\n--{boundary}--\r\n".encode()
        self.path = path
        self.length = len(self.head) + size + len(self.tail)

    def __iter__(self) -> Iterator[bytes]:
        yield self.head
        with self.path.open("rb") as f:
            while chunk := f.read(1 << 20):
                yield chunk
        yield self.tail


def _post(base: str, token: str, plan: Plan, snapshot: Path, size: int) -> Reply:
    """POST the plan with the snapshot's bytes (`size` of them) as the file."""
    boundary = f"vec-{plan.sha256}"
    body = _Multipart(plan.fields, snapshot, size, plan.upload_name, boundary)
    headers = {"Authorization": f"Bearer {token}", "Accept": "application/json",
               "Content-Type": f"multipart/form-data; boundary={boundary}", "Content-Length": str(body.length)}
    if plan.official:
        headers["Idempotency-Key"] = plan.idempotency_key
    req = urllib.request.Request(f"{base}/challenge/submissions", data=body, headers=headers, method="POST")
    return _call(req, UPLOAD_TIMEOUT, token)


class _Unreadable(Exception):
    """Reading the file to upload failed (its cause is the OSError); writing the snapshot did not."""


def _private_copy(src: BinaryIO, dst: Path, cap: int | None = None) -> tuple[str, int]:
    """Copy the open file `src` into a new file `dst` (mode 0600); returns the sha256 and the size of the bytes written.
    An error reading `src` is raised as _Unreadable, so that it is never taken for a full disk. With `cap`, the copy
    stops after cap + 1 bytes, more than any file it may hold."""
    h, size = hashlib.sha256(), 0
    with open(os.open(dst, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600), "wb") as fo:
        while cap is None or size <= cap:
            try:
                chunk = src.read(1 << 20 if cap is None else min(1 << 20, cap + 1 - size))
            except OSError as e:
                raise _Unreadable from e
            if not chunk:
                break
            fo.write(chunk)
            h.update(chunk)
            size += len(chunk)
    return h.hexdigest(), size


# ---------------------------------------------------------------------------------------------------------------------
# Files built on the DGX: pulled into a private copy for one call; the token stays here


def dgx_run(command: str, stderr: BinaryIO) -> subprocess.Popen:
    """Start `ssh -o BatchMode=yes DGX_HOST -- COMMAND` with its standard output piped and its standard error into
    `stderr` (a file, so that neither can fill up while the other is read). Every remote command goes through here, and
    tests replace it. ssh forwards the variables its SendEnv names, so VEC_TOKEN and every variable whose name contains
    TOKEN, in any case, are left out of its environment."""
    env = {k: v for k, v in os.environ.items() if "TOKEN" not in k.upper()}
    return subprocess.Popen(["ssh", "-o", "BatchMode=yes", DGX_HOST, "--", command], env=env,
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=stderr)


@dataclass
class DgxFile:
    """A file pulled from the DGX for one call: its private copy, and what the DGX reported just before the pull."""

    path: str  # relative to the DGX repo root, symlinks resolved
    copy: Path  # a 0600 file in a 0700 temporary directory, removed when the call ends
    sha256: str
    size: int
    head: str  # the DGX repo's HEAD ("unknown" when git cannot tell)
    status: str | None  # its `git status --porcelain -z --untracked-files=all` (None when git cannot tell)
    stamp: Path | None = None  # the private copy of the build stamp next to the file's real path, next to `copy`
    stamp_missing: str = ""  # why there is no such copy: no stamp, not a regular file, unreadable or too large


def _dgx_word(path: str) -> str:
    """A DGX path as one shell word; a leading ~/ becomes "$HOME"/, since a quoted ~ is never expanded."""
    return f'"$HOME"/{shlex.quote(path[2:])}' if path.startswith("~/") else shlex.quote(path)


def _dgx_start(command: str, stderr: BinaryIO) -> subprocess.Popen:
    try:
        return dgx_run(command, stderr)
    except OSError as e:  # no ssh on this machine
        raise SubmitRefused(f"cannot run ssh to {DGX_HOST} ({type(e).__name__}: {e})") from None


def _dgx_error(stderr: BinaryIO, code: int | None) -> str:
    stderr.seek(0)
    return f"{stderr.read().decode(errors='replace').strip() or 'no message'} (exit {code})"


def pull_from_dgx(path: str, stack: contextlib.ExitStack) -> DgxFile:
    """Pull `path` (relative to the DGX repo root, or absolute inside it) into a private temporary copy, which `stack`
    deletes, with the build stamp next to its real path. One remote call reports the file, the stamp and the repo's
    state, then `cat` streams the file, and the stamp unless it is larger than STAMP_MAX_BYTES, into the copies while
    they are hashed. Refuses a path outside the repo or with a newline or NUL, a file the DGX cannot report, and a pull
    whose bytes are not the reported ones ("transfer mismatch"). Reads no token and sends nothing to the challenge
    API."""
    if not path or "\n" in path or "\0" in path:
        raise SubmitRefused(f"refusing the DGX path {path!r}: it is empty or holds a newline or NUL")
    report = (f"bash -c {shlex.quote(DGX_REPORT)} vec-dgx {_dgx_word(DGX_REPO)} {shlex.quote(path)} "
              f"{shlex.quote(STAMP_SUFFIX)} {STAMP_MAX_BYTES}")
    with tempfile.TemporaryFile() as err:
        with _dgx_start(report, err) as proc:
            out = proc.stdout.read()
        parts = out.split(b"\0", 7)
        if len(parts) < 8 or parts[0] != b"vec-dgx":
            raise SubmitRefused(f"{DGX_HOST} cannot report {path}: {_dgx_error(err, proc.returncode)}")
    try:
        root, real, sha, size, head, stamp = (p.decode() for p in parts[1:7])
    except UnicodeDecodeError:
        raise SubmitRefused(f"{DGX_HOST} reports a path that is not UTF-8 for {path}") from None
    if "\n" in real:
        raise SubmitRefused(f"refusing the DGX path {path!r}: it resolves to {real!r}, which holds a newline")
    rel = real[len(root) + 1:] if root and real.startswith(f"{root}/") else ""
    if not (rel and re.fullmatch(r"[0-9a-f]{64}", sha) and re.fullmatch(r"[0-9]+", size)
            and re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}|unknown", head)
            and re.fullmatch(r"(?:[0-9a-f]{64}|large) [0-9]+|none|special|unreadable", stamp)):
        raise SubmitRefused(f"unexpected report from {DGX_HOST} for {path}: {out[:300]!r}")
    status = parts[7].decode(errors="replace") if proc.returncode == 0 else None  # git status failed
    try:
        tmp = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix="vec-dgx-", ignore_cleanup_errors=True)))
    except OSError as e:
        raise SubmitRefused(f"cannot make a temporary directory in {tempfile.gettempdir()} ({type(e).__name__}: {e}); "
                            "free some space there or point TMPDIR elsewhere") from None
    copy = tmp / path.rstrip("/").rsplit("/", 1)[-1]  # the name as given (as for a local file), not a link's target
    n = _pull(real, rel, copy, sha, int(size))
    kind, _, stamp_size = stamp.partition(" ")
    if kind in ("none", "special", "unreadable"):
        return DgxFile(rel, copy, sha, n, head, status, stamp_missing={
            "none": f"no build stamp next to {rel} on {DGX_HOST}",
            "special": f"{rel}{STAMP_SUFFIX} on {DGX_HOST} is a link or not a regular file",
            "unreadable": f"{DGX_HOST} cannot read {rel}{STAMP_SUFFIX}"}[kind])
    if kind == "large" or int(stamp_size) > STAMP_MAX_BYTES:
        return DgxFile(rel, copy, sha, n, head, status, stamp_missing=f"{rel}{STAMP_SUFFIX} on {DGX_HOST} is "
                       f"{stamp_size} bytes, more than a build stamp's {STAMP_MAX_BYTES // 1024} KiB")
    _pull(real + STAMP_SUFFIX, rel + STAMP_SUFFIX, stamp_path(copy), kind, int(stamp_size), cap=STAMP_MAX_BYTES)
    return DgxFile(rel, copy, sha, n, head, status, stamp=stamp_path(copy))


def _pull(real: str, rel: str, copy: Path, sha: str, size: int, cap: int | None = None) -> int:
    """`cat` the DGX file `real` (shown as `rel`) into the private file `copy` while hashing it (with `cap`, reading at
    most cap + 1 bytes); refuses unless the bytes have the reported `sha` and `size`, which it returns."""
    with tempfile.TemporaryFile() as err:
        with _dgx_start(f"cat -- {shlex.quote(real)}", err) as proc:
            try:
                got, n = _private_copy(proc.stdout, copy, cap)
            except _Unreadable as e:
                proc.kill()
                raise SubmitRefused(f"the pull of {rel} from {DGX_HOST} broke off ({e.__cause__})") from None
            except OSError as e:
                proc.kill()
                raise SubmitRefused(f"cannot write the copy into {copy.parent} ({type(e).__name__}: {e}); free some "
                                    "space there or point TMPDIR elsewhere") from None
            over = cap is not None and n > cap
            if over:
                proc.kill()  # more than a valid file: refused below, as a transfer mismatch
        if proc.returncode != 0 and not over:
            raise SubmitRefused(f"cannot pull {rel} from {DGX_HOST}: {_dgx_error(err, proc.returncode)}")
    if (got, n) != (sha, size):
        raise SubmitRefused(f"transfer mismatch: the pull of {rel} gave {n} bytes with sha256 {got[:12]}, but "
                            f"{DGX_HOST} reported {size} bytes with sha256 {sha[:12]}; check the file and run again")
    return n


@contextlib.contextmanager
def _unwind_on_hangup() -> Iterator[None]:
    """While a call runs in the main thread, SIGHUP and SIGTERM raise SystemExit instead of killing the process, so
    the call unwinds and its private copies are removed (a closed terminal or a dropped ssh session sends SIGHUP)."""
    if threading.current_thread() is not threading.main_thread():
        yield
        return

    def stop(signum, _frame):
        raise SystemExit(128 + signum)

    before = {s: signal.signal(s, stop) for s in (signal.SIGHUP, signal.SIGTERM)}
    try:
        yield
    finally:
        for s, handler in before.items():
            signal.signal(s, handler if handler is not None else signal.SIG_DFL)


def remove_stale_copies(now: float | None = None) -> list[str]:
    """Remove this user's private copies left in the temporary directory by a call that was killed outright (SIGKILL,
    power loss), which could not remove its own: directories named with COPY_PREFIXES, older than STALE_COPY_HOURS
    (no call runs that long). Returns their names."""
    now, gone = time.time() if now is None else now, []
    root = Path(tempfile.gettempdir())
    try:
        entries = list(root.iterdir())
    except OSError:
        return gone
    for d in entries:
        if not d.name.startswith(COPY_PREFIXES):
            continue
        try:
            st = d.lstat()  # a link is never followed
        except OSError:
            continue
        if stat.S_ISDIR(st.st_mode) and st.st_uid == os.getuid() and now - st.st_mtime > STALE_COPY_HOURS * 3600:
            shutil.rmtree(d, ignore_errors=True)
            gone.append(d.name)
    return gone



def upload_gate_guard(path: str | Path | None = None, now: datetime | None = None) -> Guard:
    """The guard "upload gate": the captain's rule of 2026-11-05, "Until our estimated score in leaderboard crosses 209,
    don't submit our 1st submission" (docs/PLAN.md section 8). It reads the committed record (UPLOAD_GATE, JSON) and
    passes when the captain has lifted the gate there (`lifted_by_captain` true) or when its `estimate`, our estimated
    leaderboard total from our own evidence, is above UPLOAD_GATE_MIN and an independent `review` is named. It also
    passes from the record's `deadline_exception_utc` on (an ISO time with its zone; the captain, 2026-11-05: "I am
    fine with submitting our 1st submission at Nov 27 anyway"). An unreadable record or exception time, a missing
    or non-finite estimate, or a `threshold` other than UPLOAD_GATE_MIN refuses. `now` is for tests."""
    p = Path(path) if path is not None else UPLOAD_GATE
    name = "upload gate"
    try:
        if p.stat().st_size > UPLOAD_GATE_MAX_BYTES:
            return Guard(name, False, f"{p.name} is larger than {UPLOAD_GATE_MAX_BYTES} bytes")
        rec = json.loads(p.read_text())
    except (OSError, ValueError) as e:
        return Guard(name, False, f"cannot read the captain's upload gate record {p.name} ({type(e).__name__})")
    if not isinstance(rec, dict):
        return Guard(name, False, f"{p.name} is not a JSON object")
    if rec.get("lifted_by_captain") is True:
        return Guard(name, True, f"lifted by the captain ({str(rec.get('lifted_note') or 'no note')[:120]})")
    opens = ""
    if rec.get("deadline_exception_utc") is not None:
        raw = str(rec["deadline_exception_utc"])
        try:
            when = datetime.fromisoformat(raw.replace("Z", "+00:00"))
        except ValueError:
            when = None
        if when is None or when.tzinfo is None:
            return Guard(name, False, f"{p.name}: cannot read deadline_exception_utc {raw[:40]!r} (an ISO time "
                         "with its zone)")
        if (now or datetime.now(UTC)) >= when:
            return Guard(name, True, f"the captain's deadline exception, open since {when.isoformat()}")
        opens = f"; it opens at the captain's deadline exception, {when.isoformat()}"
    try:
        threshold = float(rec.get("threshold", UPLOAD_GATE_MIN))
        estimate = float(rec["estimate"])
    except (KeyError, TypeError, ValueError):
        return Guard(name, False, f"{p.name} holds no numeric estimate")
    if threshold != UPLOAD_GATE_MIN:
        return Guard(name, False, f"{p.name} says threshold {threshold:g}; the captain's figure is {UPLOAD_GATE_MIN:g}")
    review = str(rec.get("review") or "").strip()
    if not (estimate == estimate and abs(estimate) != float("inf")):
        return Guard(name, False, f"{p.name} holds no finite estimate")
    if estimate > UPLOAD_GATE_MIN and review:
        return Guard(name, True, f"estimate {estimate:g} above {UPLOAD_GATE_MIN:g}, reviewed: {review[:120]}")
    why = (f"estimate {estimate:g} is not above {UPLOAD_GATE_MIN:g}" if estimate <= UPLOAD_GATE_MIN
           else "no independent review of the estimate is recorded")
    return Guard(name, False, f"the captain's upload gate stands: {why} (docs/PLAN.md section 8){opens}; only the "
                 "captain lifts it earlier")


def _inside(path: Path, root: Path) -> str | None:
    """`path` relative to `root` (both resolved, links followed) as a POSIX string, or None when it lies outside."""
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except (OSError, RuntimeError, ValueError):
        return None


def reference_list_guard(sha: str, board: str, pointer: str | Path | None = None,
                         repo: str | Path | None = None) -> Guard:
    """The guard "reference list" (AGENTS.md: an official upload sends only a file whose sha256 is in the current
    reference list). The pointer (REFERENCE_LIST, JSON: {"list", "list_sha256", "boards": {board: [names]}, "held":
    {name: reason}, ...}) names the list, sha256sum output at a path inside `repo` (default REPO_ROOT); without
    `pointer` it is `repo`/results/reference_list.json (REFERENCE_LIST when `repo` is None too). "held"
    (optional) names list entries no board may send, such as files a ruling replaced before a new list exists. Passes
    only when the pointer is a regular file and parses (at most REFERENCE_LIST_MAX_BYTES), the list lies inside the
    repository, exists and has the sha256 `list_sha256`, `board` has names, none of them held, each of which occurs in
    the list exactly once, and `sha` is the list's hash for one of them. Every refusal says which condition failed and
    names the pointer file; no message holds a path outside `repo`."""
    name = "reference list"
    root = Path(repo) if repo is not None else REPO_ROOT
    p = Path(pointer) if pointer is not None else REFERENCE_LIST if repo is None else root / REFERENCE_POINTER
    shown = _inside(p, root) or p.name  # never a path outside the repository

    def refuse(why: str) -> Guard:
        return Guard(name, False, f"{shown}: {why}")

    try:
        st = p.stat()
        if not (stat.S_ISREG(st.st_mode) or stat.S_ISDIR(st.st_mode)):  # a directory: IsADirectoryError below
            return refuse("the pointer is not a regular file")  # a FIFO would block the read, /dev/zero never end it
        if st.st_size > REFERENCE_LIST_MAX_BYTES:
            return refuse(f"the pointer is larger than {REFERENCE_LIST_MAX_BYTES} bytes")
        with p.open() as fh:  # read no more than the limit, whatever the file became after the stat
            text = fh.read(REFERENCE_LIST_MAX_BYTES + 1)
        if len(text) > REFERENCE_LIST_MAX_BYTES:
            return refuse(f"the pointer is larger than {REFERENCE_LIST_MAX_BYTES} bytes")
        rec = json.loads(text)
        json.dumps(rec, ensure_ascii=False).encode()  # a lone surrogate ("\ud800") in a string: UnicodeEncodeError
    except (OSError, ValueError, RecursionError) as e:  # RecursionError: JSON nested too deep to parse
        return refuse(f"cannot read the reference list pointer ({type(e).__name__})")
    if not isinstance(rec, dict):
        return refuse("the pointer is not a JSON object")
    rel, want, boards = rec.get("list"), rec.get("list_sha256"), rec.get("boards")
    if not isinstance(rel, str) or not rel.strip():
        return refuse('the pointer names no list ("list" is not a non-empty string)')
    if not isinstance(want, str) or not re.fullmatch(r"[0-9a-f]{64}", want):
        return refuse('"list_sha256" is not 64 lowercase hex characters')
    if not isinstance(boards, dict) or not all(
            isinstance(v, list) and all(isinstance(n, str) and n for n in v) for v in boards.values()):
        return refuse('"boards" is not an object of board keys with lists of names')
    held = rec.get("held", {})
    if not isinstance(held, dict) or not all(isinstance(n, str) and n and isinstance(r, str) and r.strip()
                                             for n, r in held.items()):
        return refuse('"held" is not an object of names with reasons')
    parts = Path(rel).parts
    if Path(rel).is_absolute() or rel.startswith("~") or ".." in parts or re.search(r"[\x00-\x1f]", rel):
        return refuse("the list it names lies outside the repository (an absolute path, ~, .. or a control "
                      "character)")
    list_path = root / rel
    inside = _inside(list_path, root)
    if inside is None:
        return refuse("the list it names lies outside the repository (through a link)")
    try:
        st = list_path.stat()
        if not stat.S_ISREG(st.st_mode):
            return refuse(f"the list {inside} is not a regular file")
        if st.st_size > REFERENCE_SUMS_MAX_BYTES:
            return refuse(f"the list {inside} is larger than {REFERENCE_SUMS_MAX_BYTES} bytes")
        raw = list_path.read_bytes()
    except FileNotFoundError:
        return refuse(f"the list {inside} does not exist")
    except OSError as e:
        return refuse(f"cannot read the list {inside} ({type(e).__name__})")
    got = hashlib.sha256(raw).hexdigest()
    if got != want:
        return refuse(f"the list {inside} has sha256 {got[:12]}, not the pointer's list_sha256 {want[:12]}")
    entries: list[tuple[str, str]] = []
    try:
        lines = raw.decode().splitlines()
    except UnicodeDecodeError:
        return refuse(f"the list {inside} is not UTF-8")
    for i, line in enumerate(lines, 1):
        if not line.strip():
            continue
        m = SUMS_LINE.match(line)
        if m is None:
            return refuse(f"line {i} of the list {inside} is not sha256sum output")
        entries.append((m.group(1), m.group(2)))
    names = boards.get(board)
    if not names:
        return refuse(f"the pointer lists no names for board {board}")
    counts: dict[str, int] = {}
    for _, n in entries:
        counts[n] = counts.get(n, 0) + 1
    for n in names:
        if n in held:
            return refuse(f"{board}'s name {n} is held by the pointer ({held[n]})")
        if counts.get(n, 0) != 1:
            return refuse(f"{board}'s name {n} occurs {counts.get(n, 0)} times in the list {inside}, not once")
    hits = [n for h, n in entries if h == sha and n in names]
    list_name = Path(inside).name
    if hits:
        return Guard(name, True, f"{list_name}: {hits[0]}")
    others = [n for h, n in entries if h == sha]
    if others:
        return refuse(f"sha256 {sha[:12]} is in the list {list_name} only as {', '.join(others[:5])}, none of them "
                      f"one of {board}'s names ({', '.join(names)})")
    return refuse(f"sha256 {sha[:12]} is not in the list {list_name} ({board}'s names: {', '.join(names)})")


def _dgx_provenance(dgx: DgxFile, changed: list[str] | None, local_head: str) -> Guard:
    """An official upload of a pulled file needs the DGX repo at the local HEAD and a tree that passes the clean-tree
    rule."""
    if dgx.head == "unknown" or changed is None:
        return Guard("dgx provenance", False, f"cannot read the git state of {DGX_HOST}:{DGX_REPO}")
    problems = [] if dgx.head == local_head else [f"DGX HEAD {dgx.head[:12]} is not the local HEAD {local_head[:12]}"]
    if changed:
        more = " ..." if len(changed) > 5 else ""
        problems.append(f"uncommitted or untracked on the DGX: {', '.join(changed[:5])}{more}")
    return Guard("dgx provenance", not problems, "; ".join(problems) or f"DGX HEAD {dgx.head[:12]} is the local HEAD "
                 "and the DGX tree is clean")


def _build_stamp(stamp: Path | None, sha: str, repo: Path, missing: str = "") -> tuple[Guard, str]:
    """The guard "build stamp" and the ledger's build_commit. An official upload needs the file's stamp (`stamp`; None
    with the reason `missing`), parsed, for the bytes to be sent (`sha`), from a clean build (vec.provenance: a clean
    tree when the build started and when the file was written, at one HEAD that did not move), at a commit that is the
    local HEAD or an ancestor of it in `repo`. build_commit is the stamp's commit (+dirty if it says the build was not
    clean) when the stamp matches the file, else ""."""
    def refuse(detail: str, build: str = "") -> tuple[Guard, str]:
        return Guard("build stamp", False, detail), build

    if stamp is None:
        return refuse(missing or "no build stamp")
    try:
        s = read_stamp(stamp)
    except FileNotFoundError:
        return refuse(f"no build stamp next to the file ({stamp.name}): write it with vec.io.write_submission, as "
                      "every builder does")
    except (OSError, ValueError) as e:
        return refuse(f"cannot read the build stamp {stamp.name}: {e}")
    if s["sha256"] != sha:
        return refuse(f"the build stamp is for other bytes (sha256 {s['sha256'][:12]}, the file's is {sha[:12]}): "
                      "the file was changed or rebuilt without a new stamp")
    commit = s["commit"]
    build = commit if s["clean"] or commit == "unknown" else f"{commit}+dirty"
    if not s["clean"]:  # the stamp says why; a stamp without a reason is described from its paths
        paths, n = s["unclean"], s["unclean_count"]
        why = s.get("why") or (
            f"uncommitted or untracked: {', '.join(paths[:5])}{' ...' if (n or 0) > 5 else ''}" if paths else
            "its git state could not be read" if commit == "unknown" or paths is None else "no reason recorded")
        return refuse(f"the build was not clean ({why}): rebuild the file from a clean tree, in a new process", build)
    if commit == "unknown":
        return refuse("the build stamp has no commit", build)
    ancestor = is_ancestor(repo, commit)
    if ancestor is None:
        return refuse(f"build commit {commit[:12]} is not in the local repository (fetch it, or rebuild at a local "
                      "commit)", build)
    if not ancestor:
        return refuse(f"build commit {commit[:12]} is neither the local HEAD nor an ancestor of it (built on a branch "
                      "that is not merged?)", build)
    return Guard("build stamp", True, f"built from a clean tree at {commit[:12]}, the local HEAD or an ancestor of "
                 "it"), build


# ---------------------------------------------------------------------------------------------------------------------
# Plan, send, resolve, status


def prepare(file: str | Path, board: str, official: bool, confirm: str | None, *, disclose: str | Sequence[str],
            note: str, model: str | None, architecture: str | None, token: str, base: str, rows: list[dict],
            repo: Path, ledger_path: Path, boards_dir: Path, external_data: Path | None,
            dgx: DgxFile | None = None) -> Plan:
    """Hash the file, derive the request and evaluate every guard. Official plans read phase, quota and submissions.
    `dgx` describes a file pulled from the DGX (`file` is then its copy)."""
    path = Path(file)
    sha = sha256_file(path)
    guards: list[Guard] = []
    warnings: list[str] = []
    b, ident = None, {}
    try:
        b = get_board(board, boards_dir)
    except KeyError as e:
        guards.append(Guard("board", False, str(e).strip("\"'")))
    else:
        try:
            ident = board_fields(b)
            guards.append(Guard("board", True, ", ".join(f"{k}={v}" for k, v in ident.items())))
        except ValueError as e:
            guards.append(Guard("board", False, str(e)))
    if b is not None:
        try:
            rep = check_local(path, b)
            warnings += rep.warnings
            count = (("max_cells" if rep.n_cells == b.max_cells else "count not shown") if b.counts_hidden
                     else rep.n_cells)
            guards.append(Guard("local format check", True, f"{count} cells, {len(rep.warnings)} warning(s)"))
        except SubmitRefused as e:
            guards.append(Guard("local format check", False, str(e)))

    extra = {k: v for k, v in (("model", model), ("architecture", architecture)) if v}
    bad = [k for k, v in extra.items() if len(v) > 200 or re.search(r"[\x00-\x1f\"]", v)]
    if bad:
        guards.append(Guard("fields", False,
                            f"{', '.join(bad)}: at most 200 characters, no quotes or control characters"))
    fields = {**ident, **extra}
    if not official:
        fields["format_only"] = "true"

    ids = parse_ids(disclose)
    register = disclosure_register(external_data)
    unknown = [i for i in ids if i not in register]
    if unknown:
        guards.append(Guard("disclosure", False, f"not in EXTERNAL_DATA.md: {', '.join(unknown)}"))
    elif official and not ids:
        guards.append(Guard("disclosure", False, "an official upload needs --disclose with EXTERNAL_DATA.md ids"))
    elif ids:
        guards.append(Guard("disclosure", True, "; ".join(f"{i} {register[i][:60]}" for i in ids)))

    stray = resolve_commands(_unresolved(rows))
    if stray:  # the server may be ignoring format_only, and not on one board only: no POST anywhere until someone looks
        seen = ", ".join(f"{sha[:12]} on {b} ({kind})" for (b, sha), (kind, _) in stray.items())
        guards.append(Guard("unresolved upload", False, f"possible submissions from format checks ({seen}): nothing "
                            "is sent to any board until a person has checked the portal and run "
                            + "; ".join(cmd for _, cmd in stray.values())))

    # the build stamp next to a local file's real path, or the copy pulled with a DGX file; it and the ledger are exempt
    # from the clean-tree rule
    stamp = dgx.stamp if dgx is not None else stamp_path(path.resolve())
    built, build_commit = _build_stamp(stamp, sha, repo, dgx.stamp_missing if dgx is not None else "")
    commit, changed = git_state(repo, ignore=[ledger_path] if dgx is not None else [ledger_path, stamp])
    git_commit = commit if changed is None or not changed else f"{commit}+dirty"
    name, key = upload_name(path), idempotency_key(board, sha) if official and ident else ""
    plan = Plan(path, board, official, sha, name, fields, key, ids, note, git_commit, guards, warnings,
                build_commit=build_commit)
    provenance = None
    if dgx is not None:  # the ledger names the DGX file and the DGX repo's state; only official uploads require it
        skip = _repo_paths(repo, [ledger_path]) | {dgx.path + STAMP_SUFFIX}
        dgx_changed = None if dgx.status is None else unclean_paths(dgx.status, skip)
        state = "+unknown" if dgx_changed is None else "+dirty" if dgx_changed else ""  # unknown: git status failed
        plan.remote, plan.origin = dgx.path, f"{DGX_HOST}:{dgx.path}@{dgx.head}{state}"
        provenance = _dgx_provenance(dgx, dgx_changed, commit)
    listed = reference_list_guard(sha, board, repo=repo)  # the plan's sha256: for a pulled file, its copy's
    if not official:
        if confirm:
            warnings.append("--confirm is ignored without --official: this is a format check")
        if not built.ok:
            warnings.append(f"an official upload of this file would fail build stamp: {built.detail}")
        if provenance is not None and not provenance.ok:
            warnings.append(f"an official upload of this file would fail dgx provenance: {provenance.detail}")
        if not listed.ok:
            warnings.append(f"an official upload of this file would fail reference list: {listed.detail}")
        if b is not None and ident:  # a format check the server stored anyway would land on the phase's split
            guards.append(_phase_guard(base, b, official=False))
        return plan

    ok = (confirm or "").strip().lower() == sha[:12]
    guards.insert(0, Guard("confirm", ok, "matches sha256[:12]" if ok else
                           "--confirm must equal the first 12 hex characters of the file's sha256"))
    if changed is None:
        guards.append(Guard("clean tree", False, f"cannot read the git state of {repo}"))
    else:
        guards.append(Guard("clean tree", not changed, "nothing uncommitted (untracked files only under results/, "
                            "runs/ and experiments/**/results/)" if not changed else
                            f"uncommitted or untracked: {', '.join(changed[:5])}{' ...' if len(changed) > 5 else ''}"))
    guards.append(upload_gate_guard())
    guards.append(built)
    if provenance is not None:
        guards.append(provenance)
    guards.append(listed)  # always, so that a refused plan shows every reason
    _reuse_earlier_request(plan, rows)
    fmt_ok = any(r.get("mode") == "format" and r.get("board") == board and r.get("sha256") == sha
                 and r.get("server_status") == "format_ok" for r in rows)
    again = f"{shlex.quote(dgx.path)} --board {board} --from-dgx" if dgx is not None else f"{path} --board {board}"
    guards.append(Guard("server format check", fmt_ok, "passed earlier for this file and board" if fmt_ok else
                        f"run a format check first: vec submit {again}"))
    if b is None or not ident:  # the board guard has failed; no server reads
        return plan

    guards.append(_phase_guard(base, b, official=True))
    q = _get(base, "/challenge/integration/quota", token)
    quota = official_quota(q.body if q.status == 200 else None, b)
    if q.status != 200:
        quota = f"the quota endpoint could not be read ({q.detail or 'no response'})"
    guards.append(Guard("quota", isinstance(quota, tuple),
                        f"server: {quota[0]} of {quota[1]} used on {board}" if isinstance(quota, tuple) else quota))
    s = _get(base, "/challenge/submissions", token)
    subs = _subs_list(s)
    guards.append(Guard("submission list", subs is not None, f"{len(subs)} server entries read" if subs is not None
                        else f"cannot read /challenge/submissions ({s.detail or s.status})"))

    plan.replay = _find_replay(board, sha, rows, subs, ident)
    if plan.replay is not None:
        guards.append(Guard("budget", True, f"already sent ({plan.replay['source']} record {plan.replay['id'] or '?'}, "
                            f"status {plan.replay['status']}): replayed, nothing new is spent"))
        return plan
    held = _ledger_files(rows, board).get(sha)
    if held is not None and held.format_ids and not held.ambiguous:  # an unanswered official attempt may still retry
        guards.append(Guard("already held", False, f"the board already holds this file as submission "
                            f"{', '.join(sorted(held.format_ids))} (a format check the server stored): an official "
                            "upload would spend another slot on the same bytes"))
    budget = board_budget(board, rows, quota, subs, ident)
    spent = budget.spent_besides(sha)
    room = isinstance(quota, tuple) and subs is not None and spent < budget.limit
    retry = "; this file's unanswered attempt is not counted again" if sha in budget.retryable else ""
    guards.append(Guard("budget", room, f"{spent} of {budget.limit} spent (server quota "
                        f"{budget.server_used if budget.server_used is not None else 'unknown'}; server entries not "
                        f"rejected: {len(budget.listed)}; unlisted in the ledger: {sum(budget.unlisted.values())}; "
                        f"ambiguous: {len(budget.ambiguous)}{retry})"))
    return plan


def _phase_guard(base: str, b: Board, official: bool) -> Guard:
    """The phase must take the board's split, since the server routes T1 and T2 uploads to the phase's split by itself;
    an official upload also needs accepts_submissions. The phase is public: no token."""
    phase, split = _get(base, "/challenge/phase"), board_split(b)
    if phase.status != 200 or phase.body is None:
        return Guard("phase", False, f"cannot read /challenge/phase ({phase.detail or phase.status})")
    p_split, accepts = str(phase.body.get("split", "")), phase.body.get("accepts_submissions")
    ok = split is not None and p_split.lower().startswith(split) and (accepts is True or not official)
    return Guard("phase", ok, f"phase {phase.body.get('phase')} split {p_split!r}, accepts submissions {accepts}; "
                 f"board split {b.split or split!r}")


def _reuse_earlier_request(plan: Plan, rows: list[dict]) -> None:
    """A retry of an official upload reuses the earlier filename, fields and Idempotency-Key byte for byte."""
    earlier = [r for r in rows if r.get("mode") == "official" and r.get("board") == plan.board
               and r.get("sha256") == plan.sha256 and r.get("server_status") != "refused" and r.get("idempotency_key")]
    if not earlier:
        return
    first = earlier[0]
    try:
        fields = json.loads(first["fields"])
    except ValueError:
        fields = None
    if not isinstance(fields, dict) or not first.get("upload_name"):
        plan.guards.append(Guard("retry", False, f"the earlier attempt's ledger row ({first['utc_time']}) is "
                                 "unreadable"))
        return
    if fields != plan.fields or first["upload_name"] != plan.upload_name:
        plan.warnings.append(f"retry of an earlier official attempt ({first['utc_time']}): reusing its filename "
                             f"{first['upload_name']!r} and fields {fields}")
    plan.fields, plan.upload_name, plan.idempotency_key = fields, first["upload_name"], first["idempotency_key"]


def _find_replay(board: str, sha: str, rows: list[dict], subs: list[dict] | None, ident: dict) -> dict | None:
    """The earlier record of this file on this board: the server's entry, or the ledger's answer to an official upload
    (a 2xx with the server's submission id and status). Never a format check that came back as a submission, which was
    made without the key, and never a ledger row without an id: the identical request goes out again then, and under
    the same key the server answers it with the original if it holds one."""
    f = _ledger_files(rows, board).get(sha) or _File()
    if subs is not None:
        same = [e for e in server_entries(subs, board, ident)
                if e.get("sha256") == sha and str(e.get("id") or "") not in f.format_ids]
        if same:
            e = next((x for x in same if x.get("status") != REJECTED), same[0])
            return {"source": "server", "id": str(e.get("id") or ""), "status": str(e.get("status") or "unknown")}
    if f.answer is not None:
        return {"source": "ledger", "id": f.answer["server_id"], "status": f.answer["server_status"]}
    return None


def submit(file: str | Path, board: str, official: bool = False, confirm: str | None = None, *,
           disclose: str | Sequence[str] = (), note: str = "", model: str | None = None,
           architecture: str | None = None, from_dgx: bool = False, on_plan: Callable[[Plan], None] | None = None,
           echo: Callable[[str], None] | None = None, sleep: Callable[[float], None] = time.sleep,
           api_base: str | None = None, token_file: Path | None = None, ledger: Path | None = None,
           repo: Path | None = None, boards_dir: Path | None = None, external_data: Path | None = None) -> Result:
    """Format check (default) or official upload of one prediction file, with every attempt in the ledger.

    Sends format_only=true unless `official` is True and `confirm` equals sha256[:12] of the file. Raises SubmitRefused
    (after `on_plan` has seen the guard results) when a guard fails; nothing is sent then. With `from_dgx`, `file` is a
    path in the DGX repo: pull_from_dgx copies it here for this call only, and the copy is what is checked and sent.
    """
    base = _check_base(api_base or API_BASE)
    if not from_dgx and not Path(file).is_file():
        raise SubmitRefused(f"no such file: {file}")
    token = vec_data.load_token(token_file) if token_file else vec_data.load_token()
    ledger_path = Path(ledger or LEDGER)
    remove_stale_copies()
    # a pulled copy, deleted on the way out (SIGHUP and SIGTERM unwind too), and the ledger lock
    with _unwind_on_hangup(), contextlib.ExitStack() as stack:
        dgx = pull_from_dgx(str(file), stack) if from_dgx else None
        led = stack.enter_context(Ledger(ledger_path))
        plan = prepare(dgx.copy if dgx else file, board, official, confirm, disclose=disclose, note=note, model=model,
                       architecture=architecture, token=token, base=base, rows=led.rows(), repo=Path(repo or REPO_ROOT),
                       ledger_path=ledger_path, boards_dir=Path(boards_dir or BOARDS_DIR), external_data=external_data,
                       dgx=dgx)
        if on_plan is not None:
            on_plan(plan)
        if not plan.ok:
            led.append(_row(plan, "refused", detail="failed: " + ", ".join(plan.failed)))
            raise SubmitRefused("guards failed: " + ", ".join(plan.failed), plan)
        if plan.replay is not None:
            return _replay(plan, led, base, token)
        return _send(plan, led, base, token, echo or (lambda _msg: None), sleep)


def _replay(plan: Plan, led: Ledger, base: str, token: str) -> Result:
    status, sid = plan.replay["status"], plan.replay["id"]
    if sid:  # refresh the status (read-only)
        r = _get(base, f"/challenge/submissions/{urllib.parse.quote(sid, safe='')}", token)
        if r.status == 200 and r.body and r.body.get("status"):
            status = str(r.body["status"])
    detail = f"replay of the earlier {plan.replay['source']} record: nothing was sent"
    led.append({**_row(plan, status, sid, detail=detail), "mode": "replay"})  # never counts or settles a slot
    return Result(plan, status, sid, None, detail, replayed=True)


def _send(plan: Plan, led: Ledger, base: str, token: str, echo: Callable[[str], None],
          sleep: Callable[[float], None]) -> Result:
    if plan.official == ("format_only" in plan.fields) or (not plan.official and plan.fields["format_only"] != "true"):
        raise SubmitRefused("internal error: format_only does not match the mode", plan)  # never reached
    if plan.official and not any(g.name == "confirm" and g.ok for g in plan.guards):
        raise SubmitRefused("internal error: official upload without a matching --confirm", plan)  # never reached
    with contextlib.ExitStack() as stack:  # every attempt streams the snapshot; it is deleted on the way out
        snapshot, size = _snapshot(plan, led, stack)
        attempts = len(RETRY_DELAYS) + 1
        for i in range(attempts):
            tag = f"attempt {i + 1} of {attempts}"
            if plan.official:
                led.append(_row(plan, "sending", detail=tag))
            reply = _post(base, token, plan, snapshot, size)
            res = _outcome(plan, reply)
            led.append(_row(plan, res.server_status, res.server_id, res.http_status,
                            f"{tag}: {res.detail}" if res.detail else tag))
            for board in res.boards:  # a format check's submission counts and blocks where the server put it, too
                if board != plan.board:
                    led.append({**_row(plan, res.server_status, res.server_id, res.http_status,
                                       f"{tag}: the server named this board for submission {res.server_id}; the "
                                       f"request named {plan.board}"), "board": board})
            if not (reply.status == 0 or reply.status >= 500) or i == attempts - 1:
                return res
            echo(f"{tag}: {res.detail or reply.status}; retrying the identical request in {RETRY_DELAYS[i]:.0f} s")
            sleep(RETRY_DELAYS[i])
    raise AssertionError("unreachable")


def _snapshot(plan: Plan, led: Ledger, stack: contextlib.ExitStack) -> tuple[Path, int]:
    """Copy the plan's file into a 0600 file in a 0700 temporary directory (mkdtemp), which `stack` deletes, hashing the
    copy; returns its path and size. The file is opened before the directory exists, so that one that has gone is not
    blamed on the disk. Refuses, with a ledger row, if the file cannot be read, the copy cannot be written, or the
    copy's sha256 is not the plan's."""
    def refuse(detail: str, message: str) -> SubmitRefused:
        led.append(_row(plan, "refused", detail=detail))
        return SubmitRefused(message, plan)

    def unreadable(e: BaseException | None) -> SubmitRefused:
        return refuse(f"cannot read the file ({type(e).__name__})",
                      f"cannot read {plan.file} ({type(e).__name__}: {e}); check it and run again")

    try:
        src = plan.file.open("rb")
    except OSError as e:
        raise unreadable(e) from None
    with src:
        try:
            tmp = stack.enter_context(tempfile.TemporaryDirectory(prefix="vec-submit-", ignore_cleanup_errors=True))
        except OSError as e:
            where = tempfile.gettempdir()
            raise refuse(f"cannot make a temporary directory ({type(e).__name__})", f"cannot make a temporary "
                         f"directory in {where} ({type(e).__name__}: {e}); free some space there or point TMPDIR "
                         "elsewhere") from None
        snapshot = Path(tmp) / "snapshot"
        try:
            sha, size = _private_copy(src, snapshot)
        except _Unreadable as e:
            raise unreadable(e.__cause__) from None
        except OSError as e:
            raise refuse(f"cannot copy the file ({type(e).__name__})", f"cannot copy the file into {tmp} "
                         f"({type(e).__name__}: {e}); free some space there or point TMPDIR elsewhere") from None
    if sha != plan.sha256:
        raise refuse("file changed before upload", "file changed before upload: its sha256 no longer matches the plan; "
                     "check it and run again")
    return snapshot, size


def _outcome(plan: Plan, reply: Reply) -> Result:
    if reply.status == 0:
        return Result(plan, "no_response", detail=reply.detail)
    if not 200 <= reply.status < 300:
        return Result(plan, "http_error", http_status=reply.status, detail=reply.detail)
    body = reply.body or {}
    problems = "; ".join(map(str, body.get("problems") or []))
    if not plan.official:
        if body.get("id"):  # recorded under the board the server named, and under the plan's
            sid, bundle = str(body["id"]), body.get("bundle_key")
            boards = tuple(dict.fromkeys([bundle, plan.board] if isinstance(bundle, str) and bundle else [plan.board]))
            cmds = "; ".join(f"vec submit resolve {plan.sha256[:12]} --board {shlex.quote(b)} --as acknowledged --note "
                             "TEXT" for b in boards)
            return Result(plan, "UNEXPECTED_SUBMISSION", sid, reply.status,
                          f"a format check came back as submission {sid} on {' and '.join(boards)}: it counts there, "
                          f"and nothing is sent to any board until a person has checked the portal and run {cmds}",
                          boards=boards)
        if body.get("format_only") is False:
            return Result(plan, "format_unconfirmed", "", reply.status,
                          "the server says format_only is false, so this may have been stored as a submission: it "
                          f"counts on {plan.board}, and nothing is sent to any board until a person has checked the "
                          f"portal and run vec submit resolve {plan.sha256[:12]} --board {shlex.quote(plan.board)} "
                          "--as acknowledged (it is there) or --as not-stored --note TEXT")
        return Result(plan, "format_ok" if body.get("ok") is True else "format_failed", "", reply.status, problems)
    notes = [f"problems: {problems}"] if problems else []
    if body.get("reject_reason"):
        notes.append(f"reject_reason: {body['reject_reason']}")
    if body.get("sha256") and body["sha256"] != plan.sha256:
        notes.append(f"SERVER SHA256 DIFFERS: {body['sha256']}")
    if body.get("bundle_key") and body["bundle_key"] != plan.board:
        notes.append(f"SERVER BOARD DIFFERS: {body['bundle_key']}")
    sid, status = str(body.get("id") or ""), str(body.get("status") or "")
    if not (sid and status):  # stored or not: only a listing, the identical retry or a person can tell
        notes.insert(0, f"the answer has no submission id and status (id {sid or '-'}, status {status or '-'}): the "
                        "upload may or may not be stored, and it counts until the server lists it")
        return Result(plan, "unconfirmed", sid, reply.status, "; ".join(notes))
    return Result(plan, status, sid, reply.status, "; ".join(notes))


def resolve(sha: str, board: str, decision: str, note: str, *, ledger: Path | None = None) -> dict:
    """Record a person's decision about a file whose outcome on `board` the ledger cannot settle. Local: nothing is
    sent and no token is read.

    not-stored: the person checked the portal and the file is not stored; its unanswered official attempts and its
    format checks answered format_unconfirmed stop counting and blocking (a later server listing still counts it).
    acknowledged: a format check of the file came back as a submission (or unconfirmed); it keeps counting on `board`
    and stops blocking. POSTs resume once no format check waits on any board; a submission recorded under two boards is
    resolved on each. Returns the ledger row; raises SubmitRefused, writing nothing, when there is nothing to resolve.
    """
    prefix = sha.strip().lower()
    if decision not in RESOLUTIONS:
        raise SubmitRefused(f"--as must be one of: {', '.join(RESOLUTIONS)}")
    if not re.fullmatch(r"[0-9a-f]{12,64}", prefix):
        raise SubmitRefused("give at least the first 12 hex characters of the file's sha256")
    if not note.strip():
        raise SubmitRefused("a note is required: say what was checked (for example the portal's submission list)")
    with Ledger(Path(ledger or LEDGER)) as led:
        rows = led.rows()
        files = {s: f for s, f in _ledger_files(rows, board).items() if s.startswith(prefix)}
        if len(files) != 1:
            raise SubmitRefused(f"{prefix} matches {len(files)} files on {board} in the ledger; give more characters"
                                if files else f"no ledger row for {prefix} on {board}")
        [(full, f)] = files.items()
        if decision == "not-stored" and not (f.ambiguous or f.unconfirmed):
            why = (f"the server answered with submission {f.answer['server_id']} ({f.answer['server_status']})"
                   if f.answer else "a format check of it came back as a submission: use --as acknowledged"
                   if f.format_ids else "no attempt of it is ambiguous")
            raise SubmitRefused(f"nothing to resolve as not-stored for {full[:12]} on {board}: {why}")
        if decision == "acknowledged" and not f.unresolved:
            raise SubmitRefused(f"nothing to acknowledge for {full[:12]} on {board}: no format check of it is waiting "
                                "for a person")
        last = next(r for r in reversed(rows) if r.get("board") == board and r.get("sha256") == full)
        settled = ("its unanswered attempts and format_unconfirmed checks stop counting; a server listing still counts "
                   "it" if decision == "not-stored" else
                   "the format check's submission keeps counting; POSTs resume once no format check waits on any board")
        row = {**dict.fromkeys(LEDGER_FIELDS, ""), "utc_time": _now(), "board": board, "mode": "resolve",
               "file": last.get("file", ""), "upload_name": last.get("upload_name", ""), "sha256": full,
               "server_status": decision, "note": note.strip(), "detail": f"resolved by a person: {settled}",
               "origin": last.get("origin", ""),
               "server_id": ";".join(r["server_id"] for r in f.unresolved if r.get("server_id"))
               if decision == "acknowledged" else ""}
        led.append(row)
    return row


@dataclass
class BoardStatus:
    board: str
    budget: Budget
    ledger_rows: int
    format_ok: int  # distinct files with a passing server format check
    entries: list[dict]  # server submissions on this board
    uploads: list[dict] = field(default_factory=list)  # the ledger's latest official row per file sent: its server
    # status and id, git_commit (the local HEAD when it was sent) and build_commit (the commit that built it)


@dataclass
class StatusReport:
    phase: dict | None
    quota: dict | None
    notes: list[str]
    boards: list[BoardStatus]


def status(board: str | None = None, *, api_base: str | None = None, token_file: Path | None = None,
           ledger: Path | None = None, boards_dir: Path | None = None) -> StatusReport:
    """Read-only view: phase, quota and submission list (GETs), and the ledger, per board."""
    base = _check_base(api_base or API_BASE)
    rows = read_ledger(ledger)
    notes: list[str] = []
    phase = _get(base, "/challenge/phase")
    if phase.status != 200:
        notes.append(f"phase: {phase.detail or 'no response'}")
    quota, subs = None, None
    try:
        token = vec_data.load_token(token_file) if token_file else vec_data.load_token()
    except vec_data.AuthError as e:
        notes.append(f"auth: {e} (ledger only)")
    else:
        q = _get(base, "/challenge/integration/quota", token)
        quota = q.body if q.status == 200 else None
        if quota is None:
            notes.append(f"quota: {q.detail or 'no response'}")
        s = _get(base, "/challenge/submissions", token)
        subs = _subs_list(s)
        if subs is None:
            notes.append(f"submissions: {s.detail or 'unreadable response'}")
    try:
        known = load_boards(Path(boards_dir or BOARDS_DIR))
    except (OSError, ValueError, KeyError) as e:
        known = {}
        notes.append(f"boards: {e}")
    keys = [board] if board else sorted(set(known) | {r.get("board", "") for r in rows} - {""} |
                                        {str(e["bundle_key"]) for e in subs or [] if isinstance(e, dict)
                                         and e.get("bundle_key")})
    out = []
    for key in keys:
        b = known.get(key)
        ident: dict[str, str] = {}
        if b is not None:
            try:
                ident = board_fields(b)
            except ValueError:
                pass
        q_entry = official_quota(quota, b) if b is not None and ident else "board not in boards/index.json"
        if quota is None:
            q_entry = "the quota endpoint could not be read"
        budget = board_budget(key, rows, q_entry, subs, ident)
        fmt = {r.get("sha256") for r in rows if r.get("board") == key and r.get("mode") == "format"
               and r.get("server_status") == "format_ok"}
        sent = {r.get("sha256"): r for r in rows if r.get("board") == key and r.get("mode") == "official"
                and r.get("server_status") != "refused"}
        out.append(BoardStatus(key, budget, sum(r.get("board") == key for r in rows), len(fmt),
                               server_entries(subs or [], key, ident), list(sent.values())))
    return StatusReport(phase.body if phase.status == 200 else None, quota, notes, out)
