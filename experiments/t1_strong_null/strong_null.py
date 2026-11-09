"""The approval switch of the strong-null down-tail prior (experiments/harness/TRIES.md rows T12a and T12b;
configs/t1_strong_null.yaml): `build_test.py t1 --strong-null` builds only when the config says `approved: true`. The
main session sets it at the build commit after the rows' stop rule and the captain's ruling (docs/PLAN.md decision log
2026-11-10 (d)), so that nothing can flip the build on by accident.

Importing this module reads no file and changes nothing.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "configs" / "t1_strong_null.yaml"


class Refused(RuntimeError):
    """The prior is not approved for a build: nothing is built with it."""


def load_config(path: str | Path = CONFIG) -> dict:
    return yaml.safe_load(Path(path).read_text())


def check_approved(cfg: Mapping | None = None) -> None:
    """Refuses unless configs/t1_strong_null.yaml says `approved: true` (exactly the boolean)."""
    cfg = load_config() if cfg is None else cfg
    if cfg.get("approved") is not True:
        raise Refused("configs/t1_strong_null.yaml does not say approved: true yet (the build commit after the stop "
                      "rule of TRIES rows T12a and T12b and the captain's ruling sets it)")
