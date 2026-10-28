"""Board specifications and the scoring weights of the official metric panel.

A board is one scored question (task x setting x held-out stage or knockout). The public specs come from
`boards/index.json` on the challenge site; test boards appear there at the start of the final phase. The group
weights are Table 1 of the accepted competition proposal (docs/reference/neurips-2026-competition-proposal.pdf).

Two intakes write `boards/`, and both strip test-board anchors in memory before anything is written: `update_boards`
fetches the public panel files, and `update_boards_from_file` reads a board index the user saved from the site with
the browser ("Save link as", never opened; added 2026-10-28, when the panel files answered HTTP 403 to scripts). Their
messages name boards, fields and files, never values from the index (docs/DATA_POLICY.md, do-not-use list).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
from collections.abc import Callable
from dataclasses import dataclass, field
from http.client import HTTPException
from pathlib import Path
from urllib.request import urlopen

from vec.data import API_BASE

REPO_ROOT = Path(__file__).resolve().parents[2]
BOARDS_DIR = REPO_ROOT / "boards"
PANELS_URL = "https://virtualembryo.ai/challenge/panels"
PHASE_URL = f"{API_BASE}/challenge/phase"  # public, no token: phase, split, accepts_submissions

# Only these boards keep their published anchors. Any other board (the final-phase test boards) keeps only the fields
# below: its anchors are distances between the reference stage and the hidden truth, so they encode the hidden target
# and must never be read (docs/DATA_POLICY.md, do-not-use list).
VALIDATION_BOARDS = frozenset({"T1:val", "T2:embryo:val_interp", "T2:heart:val_extrap", "T2:heart:val_interp",
                               "T3:gata4"})
VALIDATION_SPLITS = frozenset({"val", "val_interp", "val_extrap"})
SPEC_FIELDS = frozenset({"key", "label", "task", "setting", "split", "target_ko", "public", "n_genes", "genes_file",
                         "genes_sha256", "needs_coords", "obs_required", "obsm_required", "min_cells", "max_cells",
                         "truth_cells", "ref_cells", "floor_model"})
# Stages whose data are released, or released with the final phase, and may appear in the published T1 composition file.
# E10.5, the T1 validation answer, comes with the final phase (moved by the organisers from 2026-10-20 to 2026-10-27):
# a composition file that describes it is held back until the public phase endpoint reports the test split.
RELEASED_T1_STAGES = frozenset({"E8.5", "E9.5", "E10.5"})
FINAL_PHASE_T1_STAGES = frozenset({"E10.5"})

# Metrics where smaller is better. severity_slope and scale_log_ratio are scored on their absolute value.
LOWER_IS_BETTER = {"mmd_u", "variogram", "d2_shape", "neighborhood_mmd", "severity_slope", "scale_log_ratio"}
ABS_METRICS = {"severity_slope", "scale_log_ratio"}

# task -> [(group, group_weight, [(metric, within_group_weight), ...]), ...]
GROUPS: dict[str, list[tuple[str, float, list[tuple[str, float]]]]] = {
    "T1": [
        ("de_recovery", 0.25, [("de_score", 1.0)]),
        ("direction", 0.25, [("de_direction", 1.0)]),
        ("distribution", 0.30, [("mmd_u", 1.0)]),
        ("covariance", 0.20, [("variogram", 1.0)]),
    ],
    "T2": [
        ("expression_change", 0.25, [("de_score", 0.5), ("de_direction", 0.5)]),
        ("cell_state", 0.25, [("mmd_u", 0.6), ("variogram", 0.4)]),
        ("shape_growth", 0.25, [("d2_shape", 1 / 3), ("occupancy_dice", 1 / 3), ("scale_log_ratio", 1 / 3)]),
        ("local_organisation", 0.25, [("neighborhood_mmd", 1.0)]),
    ],
    "T3": [
        ("response_identity", 0.30, [("de_score", 1.0)]),
        ("response_direction", 0.25, [("de_direction", 1.0)]),
        ("response_magnitude", 0.25, [("severity_slope", 1.0)]),
        ("mutant_state", 0.20, [("mmd_u", 0.6), ("variogram", 0.4)]),
    ],
}

RANKED_METRICS = {task: [m for _, _, ms in groups for m, _ in ms] for task, groups in GROUPS.items()}


@dataclass
class Board:
    key: str
    task: str
    genes: list[str]
    needs_coords: bool
    setting: str | None = None
    split: str | None = None
    target_ko: str | None = None
    min_cells: int = 1000
    max_cells: int | None = None
    truth_cells: int | None = None
    ref_cells: int | None = None
    floor_model: str = "copy_last"
    anchors: dict[str, dict[str, float]] = field(default_factory=dict)

    @property
    def n_genes(self) -> int:
        return len(self.genes)


def read_genes(path: Path) -> list[str]:
    return [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]


def load_boards(boards_dir: Path = BOARDS_DIR) -> dict[str, Board]:
    index = json.loads((boards_dir / "index.json").read_text())
    boards = {}
    for key, spec in index.items():
        boards[key] = Board(
            key=key,
            task=spec["task"],
            genes=read_genes(boards_dir / spec["genes_file"]),
            needs_coords=bool(spec["needs_coords"]),
            setting=spec.get("setting"),
            split=spec.get("split"),
            target_ko=spec.get("target_ko"),
            min_cells=int(spec.get("min_cells", 1000)),
            max_cells=spec.get("max_cells"),
            truth_cells=spec.get("truth_cells"),
            ref_cells=spec.get("ref_cells"),
            floor_model=spec.get("floor_model", "copy_last"),
            anchors=spec.get("anchors", {}),
        )
    return boards


def get_board(key: str, boards_dir: Path = BOARDS_DIR) -> Board:
    boards = load_boards(boards_dir)
    if key not in boards:
        raise KeyError(f"unknown board {key!r}; known: {sorted(boards)}")
    return boards[key]


def is_validation_board(key: str, spec: dict) -> bool:
    return key in VALIDATION_BOARDS and spec.get("split") in VALIDATION_SPLITS


class IntakeError(ValueError):
    """The intake refused a board index or one of its files. Messages name boards, fields and files, never values."""


GENES_FILE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+-]*\.genes\.txt")  # a plain gene-list name next to index.json
REQUIRED_SPEC = ("genes_file", "needs_coords", "task")  # load_boards needs these
# The types a reduced (test) board's kept fields may have: those the validation boards publish (2026-10-28), and null
# where load_boards has a default. Any other type is refused, not dropped, so nothing nested under a kept name gets
# through and no field the builders need goes missing unnoticed.
_TEXT, _COUNT = (str, type(None)), (int, type(None))
SPEC_TYPES = {"key": _TEXT, "label": _TEXT, "task": (str,), "setting": _TEXT, "split": _TEXT, "target_ko": _TEXT,
              "public": (bool, type(None)), "n_genes": _COUNT, "genes_file": (str,), "genes_sha256": _TEXT,
              "needs_coords": (bool,), "obs_required": (list,), "obsm_required": (list,), "min_cells": _COUNT,
              "max_cells": _COUNT, "truth_cells": _COUNT, "ref_cells": _COUNT, "floor_model": _TEXT}
assert set(SPEC_TYPES) == SPEC_FIELDS
# The top-level fields of the published T1 composition file (as of 2026-10-28) and their types. Any other field is
# refused: it could describe a stage the released-stage rule never sees.
COMPOSITION_TYPES = {"note": (str,), "stages": (dict,), "shared_types": (list,), "n_shared_types": (int,)}


def _typed(value, kinds: tuple) -> bool:
    """value has one of `kinds`; true/false is never a count, and a list holds text only."""
    if isinstance(value, bool) and bool not in kinds:
        return False
    if isinstance(value, list):
        return list in kinds and all(isinstance(v, str) for v in value)
    return isinstance(value, kinds)


def _json(data: bytes, what: str):
    try:
        return json.loads(data)
    except ValueError:  # neither JSONDecodeError's nor UnicodeDecodeError's message is shown: the file stays unseen
        raise IntakeError(f"{what} is not JSON (perhaps a page the site served instead of the file; do not open it "
                          "to check)") from None


def _checked_index(data: bytes) -> dict:
    """The parsed board index, refused unless it is an object of board specs load_boards can read."""
    index = _json(data, "the board index")
    if not isinstance(index, dict) or not index:
        raise IntakeError("the board index is not a JSON object of boards")
    for key, spec in index.items():
        if not isinstance(spec, dict):
            raise IntakeError(f"board {key!r}: its spec is not a JSON object")
        missing = [f for f in REQUIRED_SPEC if f not in spec]
        if missing:
            raise IntakeError(f"board {key!r}: no {', '.join(missing)}")
        if not isinstance(spec["genes_file"], str) or not GENES_FILE.fullmatch(spec["genes_file"]):
            raise IntakeError(f"board {key!r}: genes_file is not a plain gene-list name (NAME.genes.txt)")
        if not isinstance(spec.get("genes_sha256", ""), str | None):
            raise IntakeError(f"board {key!r}: genes_sha256 is not a string")
    return index


def _gene_list_ok(data: bytes) -> bool:
    """Gene symbols, one per line: not empty, and none of the characters a JSON document or a web page needs."""
    return bool(data.strip()) and not any(c in data for c in b"{}[]<>")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sanitise_index(index: dict) -> tuple[dict, list[str]]:
    """Keep validation boards as published; reduce every other board to its submission spec (no anchors).

    A reduced board keeps only SPEC_FIELDS, each of the type SPEC_TYPES allows (otherwise IntakeError, naming the board
    and the fields only). Returns the cleaned index and the keys that were reduced. Nothing is printed, so stripped
    values are never seen.
    """
    clean, reduced = {}, []
    for key, spec in index.items():
        if is_validation_board(key, spec):
            clean[key] = spec
        else:
            clean[key] = {k: v for k, v in spec.items() if k in SPEC_FIELDS}
            odd = sorted(f"{k} ({type(v).__name__})" for k, v in clean[key].items() if not _typed(v, SPEC_TYPES[k]))
            if odd:
                raise IntakeError(f"board {key!r}: unexpected value type for {', '.join(odd)} (values not shown)")
            reduced.append(key)
    return clean, reduced


def _fetch(url: str) -> bytes:
    with urlopen(url, timeout=60) as r:  # noqa: S310 (fixed https URL)
        return r.read()


def final_phase_state(fetch: Callable[[str], bytes] = _fetch) -> str:
    """"open" when the public phase endpoint reports the test split, "closed" for any other split, "unreadable" when
    it gives no JSON object."""
    try:
        phase = json.loads(fetch(PHASE_URL))
    except (OSError, ValueError, HTTPException):
        return "unreadable"
    if not isinstance(phase, dict):
        return "unreadable"
    return "open" if str(phase.get("split", "")).lower().startswith("test") else "closed"


def final_phase_open(fetch: Callable[[str], bytes] = _fetch) -> bool:
    """True only when the public phase endpoint reports the test split; an unreadable answer counts as not open."""
    return final_phase_state(fetch) == "open"


def update_boards(boards_dir: Path = BOARDS_DIR, fetch: Callable[[str], bytes] = _fetch) -> dict:
    """Refresh boards/ from the public panel files, stripping test-board anchors in memory before anything is written.

    Gene lists are checked against their published sha256. The T1 composition file is written only if it describes
    released stages alone (E10.5 only once the final phase is open); otherwise it is left untouched and reported
    (without its contents). Everything is fetched and checked before the first file is written.
    """

    def gene_list(key: str, spec: dict) -> bytes:
        data = fetch(f"{PANELS_URL}/{spec['genes_file']}")
        if spec.get("genes_sha256") and _sha256(data) != spec["genes_sha256"]:
            raise IntakeError(f"{key}: gene list does not match its published sha256")
        if not _gene_list_ok(data):
            raise IntakeError(f"{key}: {spec['genes_file']} is not a plain gene list")
        return data

    return _intake(fetch(f"{PANELS_URL}/index.json"), gene_list, lambda: fetch(f"{PANELS_URL}/t1_composition.json"),
                   boards_dir, fetch)


def _inside_git_tree(path: Path) -> bool:
    return any((p / ".git").exists() for p in path.parents)


def update_boards_from_file(index_path: Path, boards_dir: Path = BOARDS_DIR, fetch: Callable[[str], bytes] = _fetch,
                            delete_source: bool = False) -> dict:
    """Offline intake: the same as update_boards, from a board index the user saved with the browser.

    The user saves the site's `panels/index.json` with "Save link as" (never opening it) into a new folder outside
    every git work tree; a file inside this repository or any other work tree (on this machine the repository's
    parent folder is one) is refused before it is read, so neither git nor a search from there reaches it. Each
    board's gene list comes from that folder, or from a file already in boards/ whose sha256 equals the one the index
    publishes (a list without a published sha256 must be in the folder). A `t1_composition.json` saved in the folder
    goes through update_boards' released-stage rule; without one, boards/t1_composition.json is left as it is. Only
    the public phase endpoint is fetched (for that rule). With delete_source the saved index, and a saved composition
    file, are deleted after a successful intake, unless the composition was held back only because the phase
    endpoint could not be read (both are kept for the re-run). A refusal keeps them and says so.
    """
    src = Path(index_path).expanduser().resolve()
    if src == REPO_ROOT or REPO_ROOT in src.parents or _inside_git_tree(src):
        raise IntakeError("the saved board index is inside a git work tree (this repository or a folder around it), "
                          "where git and searches reach it; it was not read. Do not open, diff or git add it: move it "
                          "into a new folder outside every git work tree (e.g. ~/vec-panels/), and if it replaced a "
                          "tracked file, restore that file with git checkout -- <path>")
    if not src.is_file():
        raise IntakeError(f"no saved board index at {src}")
    folder, comp_file, reused = src.parent, src.parent / "t1_composition.json", []
    supplied = comp_file.is_file()

    def gene_list(key: str, spec: dict) -> bytes:
        name, sha = spec["genes_file"], spec.get("genes_sha256")
        saved = folder / name
        if saved.is_file():
            if saved.resolve() == src:
                raise IntakeError(f"{key}: the gene list {name} is the saved index itself")
            data = saved.read_bytes()
            if not _gene_list_ok(data):
                raise IntakeError(f"{key}: the saved {name} is not a plain gene list (do not open it; delete it and "
                                  "save that link again)")
            if not sha or _sha256(data) == sha:
                return data
        if sha:
            for local in sorted(Path(boards_dir).glob("*.genes.txt")):
                data = local.read_bytes()
                if _sha256(data) == sha and _gene_list_ok(data):
                    reused.append(name if local.name == name else f"{name} = {local.name}")
                    return data
        why = "does not match its published sha256" if saved.is_file() else "is not in the folder"
        hint = ", and no file in boards/ has its published sha256" if sha else " (no sha256 published)"
        raise IntakeError(f"{key}: gene list {name} {why}{hint}: save it next to the index")

    def composition() -> bytes | None:
        if not supplied:
            return None
        if comp_file.resolve() == src:
            raise IntakeError("t1_composition.json in the folder is the saved index itself")
        return comp_file.read_bytes()

    try:
        summary = _intake(src.read_bytes(), gene_list, composition, boards_dir, fetch)
    except IntakeError as e:
        raise IntakeError(f"{e}. The saved index was kept for a re-run and still holds the test anchors: do not open "
                          f"it to check; delete it (rm -- {shlex.quote(str(src))}) if you will not re-run") from None
    kept = supplied and summary["final_phase"] == "unreadable"  # the re-run needs both files
    deleted, undeletable = [], []
    for f in (src, comp_file) if delete_source and not kept else ():
        if f.is_file():
            try:
                f.unlink()
                deleted.append(str(f))
            except OSError as e:
                undeletable.append(f"{f} ({e.strerror})")
    copies = sorted(p.name for p in folder.iterdir() if p.is_file() and p not in (src, comp_file)
                    and (p.suffix.lower() in {".json", ".crdownload", ".part", ".download"}
                         or p.name.lower().startswith("index")))
    return {**summary, "source": str(src), "t1_composition_supplied": supplied, "genes_reused": reused,
            "deleted": deleted, "undeletable": undeletable, "kept_for_rerun": kept, "possible_copies": copies}


def _write_all(boards_dir: Path, files: dict[str, bytes]) -> None:
    """Write every file under a temporary name first, then move each into place (index.json last), so a failed write
    replaces nothing."""
    tmp: dict[str, Path] = {}
    try:
        for name, data in files.items():
            tmp[name] = boards_dir / f".{name}.intake-tmp"
            tmp[name].write_bytes(data)
    except OSError:
        for t in tmp.values():
            try:
                t.unlink()
            except OSError:
                pass
        raise
    for name, t in tmp.items():
        os.replace(t, boards_dir / name)


def _intake(index_data: bytes, gene_list: Callable[[str, dict], bytes], composition: Callable[[], bytes | None],
            boards_dir: Path, fetch: Callable[[str], bytes]) -> dict:
    """Shared by both intakes: check and strip the index, collect every file, then write. Nothing is printed.

    Reads in update_boards' order (index, gene lists, composition, then the phase only if the composition names
    E10.5), and writes only after every read and check has passed. The summary names the fields dropped from each
    reduced board (names only, for decision point D1), never their values.
    """
    published = _checked_index(index_data)
    index, reduced = sanitise_index(published)
    dropped = {k: sorted(set(published[k]) - set(index[k])) for k in reduced if set(published[k]) - set(index[k])}
    files: dict[str, bytes] = {}
    for key, spec in index.items():
        data = gene_list(key, spec)
        if files.setdefault(spec["genes_file"], data) != data:
            raise IntakeError(f"{key}: another board names {spec['genes_file']} with a different gene list")
    extra, phase, comp_data = set(), None, composition()
    if comp_data is not None:
        comp = _json(comp_data, "t1_composition.json")
        if not isinstance(comp, dict) or not isinstance(comp.get("stages"), dict):
            raise IntakeError("t1_composition.json has no object of stages")
        odd = sorted(k for k, v in comp.items() if k not in COMPOSITION_TYPES or not _typed(v, COMPOSITION_TYPES[k]))
        if odd:
            raise IntakeError(f"t1_composition.json: unexpected top-level field(s) {', '.join(odd)} (values not "
                              "shown)")
        stages = set(comp["stages"])
        extra = stages - RELEASED_T1_STAGES
        phase = final_phase_state(fetch) if stages & FINAL_PHASE_T1_STAGES else None
        if phase not in (None, "open"):
            extra |= stages & FINAL_PHASE_T1_STAGES
        if not extra:
            files["t1_composition.json"] = (json.dumps(comp, indent=2) + "\n").encode()
    files["index.json"] = (json.dumps(index, indent=2) + "\n").encode()  # last, after the files it names
    _write_all(Path(boards_dir), files)
    return {"boards": sorted(index), "reduced": reduced, "dropped_fields": dropped,
            "t1_composition_written": comp_data is not None and not extra,
            "t1_composition_unreleased_stages": len(extra), "final_phase": phase}
