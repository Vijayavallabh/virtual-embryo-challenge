"""The guards of row T72's data steps (candidate MX2s; experiments/harness/TRIES.md, its cell 8, "FIXED NOW"): fetch,
load and tables each call guards() before any file is downloaded or opened, and refuse unless ALL of these hold:

  - row T72 is in experiments/harness/TRIES.md;
  - the register (EXTERNAL_DATA.md): X002's line holds 'consensus member of T1 candidate MX2s' and 'row T72'; X019's
    and X020's lines hold 'Candidate MX2s' and 'row T72'; X010's line holds 'skeletal-muscle competing class of row
    T71';
  - experiments/t1_ext/atlases.yaml has the sha256 pinned here (FROZEN_YAML_SHA256, filled at placement), status
    frozen, and equals its version at HEAD; the derived classes' file has the sha256 the yaml names (cm.require_frozen);
  - HEAD is the row's commit (the commit given on the command line, full 40 hex digits, which is the commit that adds
    row T72 to TRIES.md: its parent's TRIES.md has no row T72) on a clean tree (vec.provenance.git_state);
  - the step's output folder is absent and inside the checkout's own data/ (a real directory, not a link): fetch
    creates data/external/t1_ext_t72/ (and data/derived/t1_ext_t72/ must be absent too), load creates
    data/external/t1_ext_t72/reduced/, tables creates data/derived/t1_ext_t72/; each is created exclusively, so a
    folder left by a stopped data step blocks a rerun (a rerun is a new row);
  - no accession of the never-fetched list and no RAW.tar is requested (t1x_common.NEVER_FETCHED, a literal list,
    checked on every request by fetch.py).
Every manifest says tasks: [T1]. Nothing here reads a score, a truth or a test board.
"""

from __future__ import annotations

import re
import subprocess
from collections.abc import Callable, Mapping
from pathlib import Path

import t1x_common as cm

Refused = cm.Refused
ROW = "T72"
TRIES = cm.REPO / "experiments" / "harness" / "TRIES.md"
REGISTER = cm.REPO / "EXTERNAL_DATA.md"
ATLASES_REL = "experiments/t1_ext/atlases.yaml"
TRIES_REL = "experiments/harness/TRIES.md"
# Filled at placement, in the row's commit, from `sha256sum experiments/t1_ext/atlases.yaml` of the frozen yaml (row
# T72, cell 3: "atlases.yaml frozen, sha256 <filled at placement>, the same value pinned in the code").
FROZEN_YAML_SHA256 = "PENDING"
REGISTER_PHRASES: dict[str, tuple[str, ...]] = {
    "X002": ("consensus member of T1 candidate MX2s", "row T72"),
    "X019": ("Candidate MX2s", "row T72"),
    "X020": ("Candidate MX2s", "row T72"),
    "X010": ("skeletal-muscle competing class of row T71",),
}
ROW_EXTERNAL_REL = "data/external/t1_ext_t72"
ROW_DERIVED_REL = "data/derived/t1_ext_t72"
STEPS = ("fetch", "load", "tables")
_ROW_ID = re.compile(r"[A-Z][0-9]+[a-z]?")


def row_external(repo: Path | None = None) -> Path:
    return Path(repo or cm.REPO) / ROW_EXTERNAL_REL


def row_derived(repo: Path | None = None) -> Path:
    return Path(repo or cm.REPO) / ROW_DERIVED_REL


def row_in(text: str, row: str = ROW) -> bool:
    """`row` is a row ID of a TRIES.md text (a table row's first word; feng_relabel.row_in_tries's rule)."""
    for ln in text.splitlines():
        if ln[:2] == "| ":
            w = ln[2:].split("|", 1)[0].split()
            if w and _ROW_ID.fullmatch(w[0]) and w[0] == row:
                return True
    return False


def register_line(text: str, row: str) -> str | None:
    for line in text.splitlines():
        if line.startswith(f"| {row} |"):
            return " ".join(line.split())
    return None


def register_refusals(text: str) -> list[str]:
    out = []
    for row, phrases in REGISTER_PHRASES.items():
        line = register_line(text, row)
        if line is None:
            out.append(f"the register has no row {row}")
            continue
        out += [f"register row {row} does not hold {p!r}" for p in phrases if p not in line]
    return out


def default_git(repo: Path) -> Callable[..., str | None]:
    def git(*args: str) -> str | None:
        try:
            return subprocess.run(["git", *args], cwd=repo, capture_output=True, check=True).stdout.decode("utf-8")
        except (OSError, subprocess.CalledProcessError):
            return None
    return git


def default_state(repo: Path) -> tuple[str, list[str] | None]:
    from vec import provenance

    return provenance.git_state(repo)


def output_dir_ok(dst: Path, repo: Path) -> Path:
    """The folder resolves inside the checkout's own data/ (a real directory, not a link), so a campaign checkout whose
    data/external links into the reference checkout never writes there (feng_relabel.output_dir_ok's rule)."""
    data = Path(repo) / "data"
    if data.is_symlink() or not data.is_dir():
        raise Refused(f"{data} is a link or missing: the row's outputs go inside the checkout's own data directory")
    real = data.resolve()
    out = Path(dst).resolve()
    if real not in out.parents:
        raise Refused(f"the output folder must lie inside {real} (the checkout's own data directory); got {out}")
    return out


def guards(step: str, commit: str, repo: Path | None = None, git: Callable[..., str | None] | None = None,
           state: Callable[[Path], tuple[str, list[str] | None]] | None = None,
           yaml_sha256: str | None = None, cfg: Mapping | None = None) -> dict:
    """Every guard of cell 8 for one data step (before any download or file open); returns the record every manifest
    of the step carries. `git`, `state` and `yaml_sha256` default to the real ones (tests pass synthetic ones)."""
    if step not in STEPS:
        raise ValueError(step)
    repo = Path(repo or cm.REPO)
    git = git or default_git(repo)
    want = FROZEN_YAML_SHA256 if yaml_sha256 is None else yaml_sha256
    tries, register = repo / TRIES_REL, repo / "EXTERNAL_DATA.md"
    if not tries.is_file() or not row_in(tries.read_text(encoding="utf-8")):
        raise Refused(f"row {ROW} is not in {TRIES_REL}: a row is committed before its run")
    bad = register_refusals(register.read_text(encoding="utf-8") if register.is_file() else "")
    if bad:
        raise Refused("; ".join(bad))
    yp = repo / ATLASES_REL
    if want in (None, "", cm.PENDING) or not re.fullmatch(r"[0-9a-f]{64}", str(want)):
        raise Refused("row_t72.FROZEN_YAML_SHA256 is not filled: the placement commit pins the frozen yaml's sha256")
    if cm.sha256_file(yp) != want:
        raise Refused(f"{ATLASES_REL} does not have the sha256 the code pins ({want})")
    cfg = cfg if cfg is not None else cm.load_atlases(yp)
    if cfg.get("status") != "frozen":
        raise Refused(f"{ATLASES_REL} is not frozen")
    at_head = git("show", f"HEAD:{ATLASES_REL}")
    if at_head is None or at_head.encode("utf-8") != yp.read_bytes():
        raise Refused(f"{ATLASES_REL} differs from its version at HEAD")
    cm.require_frozen(cfg, repo)
    if not re.fullmatch(r"[0-9a-f]{40}", commit or ""):
        raise Refused("give the row's commit as its full 40-digit sha")
    head, unclean = (state or default_state)(repo)
    if head != commit:
        raise Refused(f"HEAD is {head}, not the row's commit {commit}")
    if unclean is None or unclean:
        raise Refused(f"the tree is not clean ({unclean}): the row's steps run on its committed tree")
    t_head, t_parent = git("show", f"{commit}:{TRIES_REL}"), git("show", f"{commit}^:{TRIES_REL}")
    if t_head is None or not row_in(t_head) or (t_parent is not None and row_in(t_parent)):
        raise Refused(f"{commit} is not the commit that adds row {ROW} to {TRIES_REL}")
    ext, der = row_external(repo), row_derived(repo)
    if step == "fetch" and (ext.exists() or der.exists()):
        raise Refused(f"{ext} or {der} exists: a folder left by a stopped data step blocks a rerun (a new row)")
    if step == "load" and ((ext / "reduced").exists() or der.exists()):
        raise Refused(f"{ext / 'reduced'} or {der} exists: a rerun is a new row")
    if step == "tables" and der.exists():
        raise Refused(f"{der} exists: a rerun is a new row")
    for d in (ext, der):
        output_dir_ok(d, repo)
    return {"row": ROW, "step": step, "commit": commit, "clean": True, "atlases_yaml_sha256": want,
            "derived_sha256": cfg["extra_competing"]["derived_sha256"], "tasks": list(cm.TASKS)}


def make_exclusive(d: Path) -> Path:
    """Create the step's output folder; refused if it exists (created exclusively)."""
    try:
        d.mkdir(parents=False, exist_ok=False)
    except FileExistsError:
        raise Refused(f"{d} exists: a rerun is a new row") from None
    except FileNotFoundError:
        d.parent.mkdir(parents=True, exist_ok=True)
        d.mkdir(exist_ok=False)
    return d
