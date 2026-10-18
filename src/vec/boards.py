"""Board specifications and the scoring weights of the official metric panel.

A board is one scored question (task x setting x held-out stage or knockout). The public specs come from
`boards/index.json` on the challenge site; test boards appear there at the start of the final phase. The group
weights are Table 1 of the accepted competition proposal (docs/reference/neurips-2026-competition-proposal.pdf).
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from urllib.request import urlopen

REPO_ROOT = Path(__file__).resolve().parents[2]
BOARDS_DIR = REPO_ROOT / "boards"
PANELS_URL = "https://virtualembryo.ai/challenge/panels"

# Only these boards keep their published anchors. Any other board (the final-phase test boards) keeps only the fields
# below: its anchors are distances between the reference stage and the hidden truth, so they encode the hidden target
# and must never be read (docs/DATA_POLICY.md, do-not-use list).
VALIDATION_BOARDS = frozenset({"T1:val", "T2:embryo:val_interp", "T2:heart:val_extrap", "T2:heart:val_interp",
                               "T3:gata4"})
VALIDATION_SPLITS = frozenset({"val", "val_interp", "val_extrap"})
SPEC_FIELDS = frozenset({"key", "label", "task", "setting", "split", "target_ko", "public", "n_genes", "genes_file",
                         "genes_sha256", "needs_coords", "obs_required", "obsm_required", "min_cells", "max_cells",
                         "truth_cells", "ref_cells", "floor_model"})
# Stages whose data are released (or released on 2026-10-20) and may appear in the published T1 composition file.
RELEASED_T1_STAGES = frozenset({"E8.5", "E9.5", "E10.5"})

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


def sanitise_index(index: dict) -> tuple[dict, list[str]]:
    """Keep validation boards as published; reduce every other board to its submission spec (no anchors).

    Returns the cleaned index and the keys that were reduced. Nothing is printed, so stripped values are never seen.
    """
    clean, reduced = {}, []
    for key, spec in index.items():
        if is_validation_board(key, spec):
            clean[key] = spec
        else:
            clean[key] = {k: v for k, v in spec.items() if k in SPEC_FIELDS}
            reduced.append(key)
    return clean, reduced


def _fetch(url: str) -> bytes:
    with urlopen(url, timeout=60) as r:  # noqa: S310 (fixed https URL)
        return r.read()


def update_boards(boards_dir: Path = BOARDS_DIR, fetch: Callable[[str], bytes] = _fetch) -> dict:
    """Refresh boards/ from the public panel files, stripping test-board anchors in memory before anything is written.

    Gene lists are checked against their published sha256. The T1 composition file is written only if it describes
    released stages alone; otherwise it is left untouched and reported (without its contents).
    """
    index, reduced = sanitise_index(json.loads(fetch(f"{PANELS_URL}/index.json")))
    for key, spec in index.items():
        if Path(spec["genes_file"]).name != spec["genes_file"]:
            raise ValueError(f"{key}: unexpected genes_file path {spec['genes_file']!r}")
        data = fetch(f"{PANELS_URL}/{spec['genes_file']}")
        if spec.get("genes_sha256") and hashlib.sha256(data).hexdigest() != spec["genes_sha256"]:
            raise ValueError(f"{key}: gene list does not match its published sha256")
        (boards_dir / spec["genes_file"]).write_bytes(data)
    (boards_dir / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    comp = json.loads(fetch(f"{PANELS_URL}/t1_composition.json"))
    extra = sorted(set(comp.get("stages", {})) - RELEASED_T1_STAGES)
    if not extra:
        (boards_dir / "t1_composition.json").write_text(json.dumps(comp, indent=2) + "\n")
    return {"boards": sorted(index), "reduced": reduced, "t1_composition_written": not extra,
            "t1_composition_unreleased_stages": len(extra)}
