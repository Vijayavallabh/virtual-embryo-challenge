"""T12's pinned table (experiments/harness/TRIES.md row T12; configs/t1_strong_null.yaml): the production table of the
strong-null down-tail prior, given to `build_test.py t1 --strong-null` only when the config pins its sha256 (null until
the build commit after T12's calibration and stop rule; R2's `hygiene.load_union` pattern), so that nothing can flip the
build on by accident and the table the build reads is the one the row's run wrote.

Importing this module reads no file and changes nothing.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

import yaml

from vec import t1_scope

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "configs" / "t1_strong_null.yaml"


class Refused(RuntimeError):
    """The table is not the one configs/t1_strong_null.yaml pins: nothing is built from it."""


def load_config(path: str | Path = CONFIG) -> dict:
    return yaml.safe_load(Path(path).read_text())


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(1 << 22):
            h.update(block)
    return h.hexdigest()


def table_record(cfg: Mapping | None = None, root: str | Path = REPO) -> dict:
    """The production table's file and pinned sha256, refused unless the config pins it and the file is the pinned
    one."""
    cfg = load_config() if cfg is None else cfg
    rel, pinned = cfg["table_file"], cfg.get("table_sha256")
    path = Path(root) / rel
    if not pinned:
        raise Refused(f"configs/t1_strong_null.yaml pins no table_sha256 for {rel} yet (the build commit after T12's "
                      "stop rule pins it)")
    if not path.exists():
        raise Refused(f"{rel}: no such file (experiments/t1_strong_null/calibrate.py writes it)")
    sha = sha256_file(path)
    if sha != pinned:
        raise Refused(f"{rel}: sha256 {sha} is not the one configs/t1_strong_null.yaml pins ({pinned})")
    return {"file": rel, "sha256": sha}


def load_table(cfg: Mapping | None = None, root: str | Path = REPO) -> t1_scope.DownPriorTable:
    """The production table (`DownPriorTable`), refused unless pinned (`table_record`)."""
    cfg = load_config() if cfg is None else cfg
    rec = table_record(cfg, root)
    return t1_scope.DownPriorTable.from_dict(json.loads((Path(root) / rec["file"]).read_text()))
