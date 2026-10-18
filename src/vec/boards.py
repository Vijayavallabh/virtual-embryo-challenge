"""Board specifications and the scoring weights of the official metric panel.

A board is one scored question (task x setting x held-out stage or knockout). The public specs come from
`boards/index.json` on the challenge site; test boards appear there at the start of the final phase. The group
weights are Table 1 of the accepted competition proposal (docs/reference/neurips-2026-competition-proposal.pdf).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BOARDS_DIR = REPO_ROOT / "boards"

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
