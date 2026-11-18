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


def test_the_heart_commands_of_build_all_are_the_ones_the_rules_give():
    """Each heart build command of experiments/builds/build_all.sh (from `build_test.py heart` up to its --out)
    appears word for word in docs/S2_RULES.md, so the script and the frozen rules cannot drift apart unnoticed (a
    review finding of 2026-11-13: after the heart ruling the script carried --trend before the rules said so)."""
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
    script = (ROOT / "experiments" / "builds" / "build_all.sh").read_text().splitlines()
    commands = [line[line.index("build_test.py heart"):line.index(" --out")].strip()
                for line in script if "build_test.py heart" in line and " --out" in line]
    # the first upload's line, the second upload's as it was frozen, and the second's since the amendment of 2026-11-18
    assert len(commands) == 3, commands
    for c in commands:
        assert c in rules, f"build_all.sh runs `{c}`, which docs/S2_RULES.md does not give"
