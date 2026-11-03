"""docs/S2_RULES.md was frozen on 2026-11-03: its sha256 must equal docs/S2_RULES.sha256. A change after the freeze
needs the captain's own say-so, in a commit that updates both files and logs his words, the reason and the new sha256
in the docs/PLAN.md decision log (docs/S2_RULES.md, "Timing and freeze")."""

from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_s2_rules_match_their_frozen_hash():
    want, name = (ROOT / "docs" / "S2_RULES.sha256").read_text().split()
    assert name == "docs/S2_RULES.md"  # `sha256sum -c docs/S2_RULES.sha256` runs from the repo root
    got = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    assert got == want, ("docs/S2_RULES.md differs from its frozen sha256: a rule changes only with the captain's "
                         "say-so, logged with his words, the reason and the new sha256")
