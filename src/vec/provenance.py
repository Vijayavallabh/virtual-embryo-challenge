"""Provenance of built files: the clean-tree rule, and the build stamp that ties a prediction file to the code that
built it.

The clean-tree rule (its single home: write_stamp applies it when a build starts and when it writes its file, vec.submit
to this tree and to the DGX's when a file is uploaded): a tree is clean when no tracked file has an uncommitted change
and no untracked, non-ignored file lies outside results/, runs/ and experiments/**/results/. A caller names the files it
exempts (`ignore`): vec.submit the ledger and the uploaded file's own stamp, write_stamp the file, its stamp and the
ledger.

A build stamp is `<file>.stamp.json` (STAMP_SUFFIX) next to the file's real path, written right after the file by
vec.io.write_submission (atomically: a temporary file in the same directory, then a rename). It holds the stamp's
version, the file's name, the sha256 and size of the file as written, the repository's HEAD when the file was written
(full hex, or "unknown"), whether the build was clean (below), the paths that kept it from being clean (at most
STAMP_MAX_PATHS of them, and their count; null when git cannot tell), why it was not clean (null when it was), the HEAD
when the build started and when that was (start_commit, started_utc), and when, where and how the file was built
(created_utc, hostname, argv).

A build is clean when its tree passed the rule above both when the build started and when the file was written, at the
same HEAD, and HEAD did not move in between: its reflog has the same entries (a checkout, reset, stash, commit or pull
adds one, also when HEAD ends where it was). The start is this process's: the state of this package's repository at the
first import of this module, which vec.io imports, so every builder records it before its heavy work. A pull, checkout,
rsync or stash while a build runs, or in a notebook between importing vec and writing a file, therefore leaves a stamp
that is not clean; only an uncommitted edit made and undone while a build runs leaves no trace git can see. A process
that stamps against another repository (tests) records its start with _record_start. Reading the git state never fails
the write: the stamp then says commit "unknown" (start_commit "unknown" when only the start could not be read) and clean
false. vec.submit's guard "build stamp" reads it back (read_stamp, is_ancestor), so that an official upload names the
code that built the file, not the code checked out when it is sent.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import secrets
import socket
import subprocess
import sys
from collections.abc import Collection, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from vec.boards import REPO_ROOT

# untracked files may stay here (and under experiments/**/results/); anything else, code or not, must be committed so
# that a recorded commit describes what produced a file
OUTPUT_DIRS = ("results", "runs")
LEDGER_PATH = "results/submissions.csv"  # vec.submit's ledger, relative to the repository root
STAMP_SUFFIX = ".stamp.json"
STAMP_VERSION = 1
STAMP_MAX_BYTES = 64 * 1024  # a stamp is read (or pulled from the DGX) only up to this size
STAMP_MAX_PATHS = 20  # unclean paths listed in a stamp; all of them are counted
FULL_COMMIT = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")  # a full SHA-1 or SHA-256 object name
SHA256 = re.compile(r"[0-9a-f]{64}")


# ---------------------------------------------------------------------------------------------------------------------
# The clean-tree rule


def _output_path(path: str) -> bool:
    """An untracked path that may stay uncommitted: under results/, runs/ or experiments/**/results/."""
    parts = path.split("/")
    return (parts[0] in OUTPUT_DIRS and len(parts) > 1) or (parts[0] == "experiments" and "results" in parts[1:-1])


def _git(repo: str | Path, *args: str) -> str:
    """The standard output of git in `repo`; raises OSError or CalledProcessError. GIT_OPTIONAL_LOCKS=0: status never
    takes the index lock, so reading the state during a build cannot make a git command in the same tree fail."""
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, check=True,
                          env={**os.environ, "GIT_OPTIONAL_LOCKS": "0"}).stdout.decode(errors="replace")


def git_state(repo: Path, ignore: Sequence[Path] = ()) -> tuple[str, list[str] | None]:
    """HEAD commit and the paths that keep the tree from being clean (None when git cannot tell; see unclean_paths),
    leaving out `ignore` (the files the caller exempts)."""
    try:
        commit = _git(repo, "rev-parse", "HEAD").strip()
        status = _git(repo, "status", "--porcelain", "-z", "--untracked-files=all")
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None
    return commit, unclean_paths(status, _repo_paths(repo, ignore))


def _repo_paths(repo: Path, paths: Sequence[Path]) -> set[str]:
    """`paths` relative to `repo`, as git status prints them; paths outside it are left out."""
    out = set()
    for p in paths:
        try:
            out.add(Path(p).resolve().relative_to(Path(repo).resolve()).as_posix())
        except ValueError:
            pass
    return out


def unclean_paths(status: str, skip: Collection[str] = ()) -> list[str]:
    """The paths in `git status --porcelain -z --untracked-files=all` output (of this tree or the DGX's) that keep the
    tree from being clean: tracked files with uncommitted changes, and untracked, non-ignored files outside results/,
    runs/ and experiments/**/results/. -z never quotes the paths, so `skip` (the exempt files) compares exactly."""
    changed, entries = [], iter(status.split("\0"))
    for entry in entries:
        code, paths = entry[:2], [entry[3:]]
        if "R" in code or "C" in code:
            paths.append(next(entries, ""))  # a rename or copy: the new path, then the old one
        changed += [p for p in paths if p and p not in skip and not (code == "??" and _output_path(p))]
    return changed


# ---------------------------------------------------------------------------------------------------------------------
# Build stamps


def stamp_path(path: str | Path) -> Path:
    """The build stamp of `path`: `<path>.stamp.json` next to it (pass the real path: a stamp sits next to the file a
    link points at)."""
    path = Path(path)
    return path.with_name(path.name + STAMP_SUFFIX)


@dataclass(frozen=True)
class RepoState:
    """A repository's state, as a build stamp compares it between the start of a build and the write of its file."""

    commit: str  # HEAD, full hex; "unknown" when git cannot tell (and then unclean is None)
    unclean: tuple[str, ...] | None  # the paths that keep the tree from being clean, before any exemption
    # HEAD's reflog: its number of entries and the newest, "<commit> HEAD@{<unix time>} <subject>" ((0, "") for none).
    # Both: entries made within the second of the newest one can repeat it exactly (away and back to the same commit)
    reflog: tuple[int, str]
    utc: str  # when it was read


def repo_state(repo: str | Path) -> RepoState:
    """The state of `repo` now. Never raises: whatever keeps git from telling (no git, not a repository, no commit, a
    failure of any kind) gives commit "unknown"."""
    utc = _utc()
    try:
        commit, changed = git_state(Path(repo))
        log = _git(repo, "reflog", "--date=unix", "--format=%H %gd %gs").splitlines() if changed is not None else []
    except Exception:  # whatever goes wrong reading the state, it is recorded as unknown
        commit, changed, log = "unknown", None, []
    if changed is None or not FULL_COMMIT.fullmatch(commit):
        return RepoState("unknown", None, (0, ""), utc)
    return RepoState(commit, tuple(changed), (len(log), log[0] if log else ""), utc)


_STARTS: dict[Path, RepoState] = {}  # each repository's state when this process's build started (_record_start)


def _record_start(repo: str | Path | None = None) -> RepoState:
    """Record the state of `repo` (default: this package's repository) now as the start of this process's build, which
    write_stamp compares with the state when it writes a stamp. This module records it for this package's repository
    at its first import. Recording it again later in the same process would hide a change made after the code was
    loaded, so only tests, and a process that stamps against another repository, call it."""
    root = Path(REPO_ROOT if repo is None else repo).resolve()
    _STARTS[root] = state = repo_state(root)
    return state


def write_stamp(path: str | Path, repo: str | Path | None = None) -> Path:
    """Stamp the file at `path`, just written, with the state of `repo` (default: this package's repository) now and
    when this process's build started (_record_start); returns the stamp's path. The file, its stamp and the ledger are
    exempt from the clean-tree rule here, at both times, so that a rebuild or a committed ledger with new rows never
    reads as unclean code. Never raises over the git state."""
    real = Path(path).resolve()
    sha, size = _sha256_size(real)
    out = stamp_path(real)
    try:
        root = Path(REPO_ROOT if repo is None else repo).resolve()
        start, now = _STARTS.get(root), repo_state(root)
        exempt = _repo_paths(root, [real, out, root / LEDGER_PATH])
    except Exception:  # whatever goes wrong reading the state, the stamp records it as unknown and the write stands
        start, now, exempt = None, RepoState("unknown", None, (0, ""), _utc()), set()
    paths = None if now.unclean is None else [p for p in now.unclean if p not in exempt]
    before = [] if start is None or start.unclean is None else [
        p for p in start.unclean if p not in exempt and p not in (paths or ())]  # unclean only when the build started
    why = _why(start, now, paths, before)
    listed = None if paths is None else paths + before
    stamp = {"version": STAMP_VERSION, "file": real.name, "sha256": sha, "size": size, "commit": now.commit,
             "clean": why is None, "why": why, "unclean": None if listed is None else listed[:STAMP_MAX_PATHS],
             "unclean_count": None if listed is None else len(listed),
             "start_commit": "unknown" if start is None else start.commit,
             "started_utc": None if start is None else start.utc, "created_utc": _utc(),
             "hostname": socket.gethostname(), "argv": list(sys.argv)}
    _write_atomic(out, json.dumps(stamp, indent=1) + "\n")
    return out


def _why(start: RepoState | None, now: RepoState, paths: list[str] | None, before: list[str]) -> str | None:
    """Why a build was not clean (None when it was), from the repository's state when it started (`start`; None when
    none was recorded) and when its file was written (`now`), and the paths that kept the tree from being clean when
    the file was written (`paths`) and only when the build started (`before`), the exempt files left out."""
    if paths is None:
        return "git cannot tell the repository's state when the file was written"
    if start is None:
        return "no record of the repository's state when the build started"
    if start.unclean is None:
        return "git could not tell the repository's state when the build started"
    if start.commit != now.commit:
        return f"HEAD moved during the build, from {start.commit[:12]} to {now.commit[:12]}"
    if paths:
        return f"uncommitted or untracked: {_some(paths)}"
    if before:
        return f"uncommitted or untracked when the build started: {_some(before)}"
    if start.reflog != now.reflog:  # left and came back: a checkout, reset or stash, as the newest entry says
        subject = now.reflog[1].split(" ", 2)[2:]
        return f"HEAD moved during the build and back to {now.commit[:12]}" + (f" ({subject[0][:100]})" if subject
                                                                                 else "")
    return None


def _some(paths: Sequence[str]) -> str:
    return f"{', '.join(paths[:5])}{' ...' if len(paths) > 5 else ''}"


def _utc() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _sha256_size(path: Path) -> tuple[str, int]:
    h, size = hashlib.sha256(), 0
    with path.open("rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
            size += len(chunk)
    return h.hexdigest(), size


def _write_atomic(path: Path, text: str) -> None:
    """Write `text` to a temporary file next to `path`, then rename it over `path`: a reader never sees a partial stamp,
    and a failed write leaves an earlier one as it was."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.{secrets.token_hex(4)}.tmp")
    try:
        with open(os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666), "w") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def read_stamp(path: str | Path) -> dict:
    """The build stamp at `path`, checked: at most STAMP_MAX_BYTES of JSON, this version, and the fields an upload reads
    with their types (sha256, size, commit: full hex or "unknown", clean, why, unclean, unclean_count). Raises OSError
    (FileNotFoundError when there is none) or ValueError, saying what is wrong."""
    with Path(path).open("rb") as f:
        raw = f.read(STAMP_MAX_BYTES + 1)
    if len(raw) > STAMP_MAX_BYTES:
        raise ValueError(f"it is larger than {STAMP_MAX_BYTES // 1024} KiB")
    try:
        s = json.loads(raw)
    except ValueError as e:  # not JSON, or not text
        raise ValueError(f"it is not JSON ({e})") from None
    if not isinstance(s, dict) or type(s.get("version")) is not int:
        raise ValueError("it is not a build stamp (no version)")
    if s["version"] != STAMP_VERSION:
        raise ValueError(f"its version is {s['version']}, not {STAMP_VERSION}")
    sha, commit, unclean, count = s.get("sha256"), s.get("commit"), s.get("unclean"), s.get("unclean_count")
    bad = [k for k, ok in (
        ("sha256", isinstance(sha, str) and SHA256.fullmatch(sha)),
        ("size", type(s.get("size")) is int and s["size"] >= 0),
        ("commit", isinstance(commit, str) and (commit == "unknown" or FULL_COMMIT.fullmatch(commit))),
        ("clean", isinstance(s.get("clean"), bool)),
        ("why", s.get("why") is None or isinstance(s["why"], str)),
        ("unclean", unclean is None or (isinstance(unclean, list) and all(isinstance(p, str) for p in unclean))),
        ("unclean_count", count is None or (type(count) is int and count >= 0))) if not ok]
    if bad:
        raise ValueError(f"malformed {', '.join(bad)}")
    return s


def is_ancestor(repo: Path, commit: str) -> bool | None:
    """Whether `commit` is the HEAD of `repo` or an ancestor of it (`git merge-base --is-ancestor`); None when git
    cannot tell: not a full hex object name, a commit the repository does not have, no HEAD, no git."""
    if not FULL_COMMIT.fullmatch(commit):
        return None
    try:
        r = subprocess.run(["git", "-C", str(repo), "merge-base", "--is-ancestor", commit, "HEAD"], capture_output=True)
    except OSError:
        return None
    return {0: True, 1: False}.get(r.returncode)


# The start of this process's build: the state of this package's repository at the first import of this module, which
# vec.io imports at its top (and the T3 builder imports vec.io first thing in its build_submission)
_record_start()
