"""The run registry, results/registry.csv: one row per scored candidate (AGENTS.md, "Evidence over claims").

Rows record how a score was made (scoring protocol, scorer build, anchors, gate, working-copy fraction), so a number
from the uncalibrated emulator protocol can never pass for a calibrated one. `vec eval --register` and the pseudo-task
harness (vec.harness, experiments/harness/run.py) both write through `append`.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping
from pathlib import Path

from vec.boards import REPO_ROOT

REGISTRY = REPO_ROOT / "results" / "registry.csv"
HEADER = ("utc_time", "commit", "tree_clean", "source", "run_id", "test_board", "pseudo_task", "scoring_board",
          "candidate", "params", "seeds", "score_mean", "score_sd", "protocol", "build", "anchors", "gate", "frac",
          "scores_json")


def append(rows: Iterable[Mapping], path: str | Path = REGISTRY) -> int:
    """Append rows (every HEADER field, nothing else) and return how many were written. A new file gets the header;
    an existing one must already carry exactly this header."""
    rows = [dict(r) for r in rows]
    for r in rows:
        if set(r) != set(HEADER):
            raise ValueError(f"registry row fields differ from the header: missing {sorted(set(HEADER) - set(r))}, "
                             f"extra {sorted(set(r) - set(HEADER))}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    new = not path.exists() or path.stat().st_size == 0
    if not new:
        with path.open(newline="") as f:
            first = next(csv.reader(f), None)
        if tuple(first or ()) != HEADER:
            raise ValueError(f"{path} has another header; move it aside before appending")
    with path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER)
        if new:
            w.writeheader()
        w.writerows(rows)
    return len(rows)
