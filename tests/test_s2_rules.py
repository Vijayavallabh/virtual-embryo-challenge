"""docs/S2_RULES.md was frozen on 2026-11-03: its sha256 must equal docs/S2_RULES.sha256. A change after the freeze
needs the captain's own say-so, in a commit that updates both files and logs his words, the reason and the new sha256
in the docs/PLAN.md decision log (docs/S2_RULES.md, "Timing and freeze")."""

from __future__ import annotations

import hashlib
import re
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


def test_the_t1_commands_of_build_all_are_the_ones_the_rules_give():
    """Each T1 build command of experiments/builds/build_all.sh (from `build_test.py t1` up to its --out) appears word
    for word in docs/S2_RULES.md, as the heart's do: the first upload's (row R6), the second as ruled on 2026-11-13
    (row R7), and T1's last upload with no outside data, by the amendment of 2026-11-19 (set t1none, row R8)."""
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
    script = (ROOT / "experiments" / "builds" / "build_all.sh").read_text().splitlines()
    commands = [line[line.index("build_test.py t1"):line.index(" --out")].strip()
                for line in script if "build_test.py t1" in line and " --out" in line]
    assert commands == ["build_test.py t1 --g-mom --base lib3",
                        "build_test.py t1 --g-mom --base mid --cm-share 0.60",
                        "build_test.py t1 --g-mom --base lib3 --feng-groups none"], commands
    for c in commands:
        # with its --out after it: the first upload's command is a prefix of the t1none command, so a bare `in` check
        # would pass on the t1none section alone
        assert c + " --out" in rules, f"build_all.sh runs `{c}`, which docs/S2_RULES.md does not give"


def test_the_t1none_set_checks_the_rules_hash_and_builds_one_file():
    """The set t1none checks docs/S2_RULES.sha256 before it builds, refuses if the rules hold the draft heading
    ("DRAFT, not in force"; gone since the amendment of 2026-11-19), builds one file, T1_s2none.h5ad, by the first
    upload's command plus --feng-groups none, and checks it locally for T1:test."""
    text = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    block = text[text.index("\n  t1none)\n"):]
    block = block[:block.index("\n    ;;\n")]
    lines = [ln.strip() for ln in block.strip().splitlines()]
    assert lines == [
        "t1none)",
        "sha256sum -c docs/S2_RULES.sha256",
        "if grep -q 'DRAFT, not in force' docs/S2_RULES.md; then "
        'echo "refused: the t1none amendment is a draft, not in force" >&2; exit 3; fi',
        "$U run python experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups none "
        "--out $O/builds/T1_s2none.h5ad",
        'files="$O/builds/T1_s2none.h5ad:T1:test"',
    ], lines
    assert "--out runs/builds/T1_s2none.h5ad" in " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())


def test_the_t1none_guard_string_marks_only_the_draft_heading():
    """The guard of set t1none greps docs/S2_RULES.md for "DRAFT, not in force". While the amendment is a draft the
    string is in its heading, so the set refuses; once the captain's words are in, the heading is rewritten and the
    string must be gone from the whole file, or the set would refuse a ruled amendment."""
    text = (ROOT / "docs" / "S2_RULES.md").read_text()
    if "[THE CAPTAIN'S WORDS]" in text:
        assert "\n## DRAFT, not in force: amendment for T1's second file" in text
    else:
        assert "DRAFT, not in force" not in text


def test_the_embryo_commands_of_build_all_are_the_ones_the_rules_give():
    """Each embryo build command of experiments/builds/build_all.sh (from `build_test.py embryo` up to its --out, its
    continuation lines joined and the script's EMBRYO_EXPRESSION put in) appears word for word in docs/S2_RULES.md.
    The amendment of 2026-11-19 (row R9 of experiments/harness/TRIES.md) made the embryo's last upload the first
    upload's command with --method fgw_shape0 and --positions floor_exact: the same expression step as the set s1."""
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
    text = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    expression = re.search(r"^EMBRYO_EXPRESSION=(\S+)", text, re.M).group(1)
    assert expression == "mixture"  # the first upload's expression step, which the last upload shares
    script = " ".join(text.replace("\\\n", " ").replace('"$EMBRYO_EXPRESSION"', expression).split())
    commands = [c.strip() for c in re.findall(r"build_test\.py embryo .*? (?=--out )", script)]
    # the first upload's line, the second upload's of ruling R4 (not uploaded), and the last upload's since 2026-11-19
    assert len(commands) == 3, commands
    for c in commands:
        assert c in rules, f"build_all.sh runs `{c}`, which docs/S2_RULES.md does not give"
    assert ("build_test.py embryo --method fgw_shape0 --positions floor_exact --expression mixture "
            "--out $O/builds/T2embryo_s2x.h5ad") in script


def test_the_embryo2_set_checks_the_rules_hash_and_builds_one_file():
    """The set embryo2 checks docs/S2_RULES.sha256 before it builds, builds one file, T2embryo_s2x.h5ad, by the one
    command of the amendment of 2026-11-19 for the embryo (row R9), and checks it locally for T2:embryo:test_interp.
    The rules give the same command with runs/ as its output folder."""
    text = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    block = text[text.index("\n  embryo2)\n"):]
    block = block[:block.index("\n    ;;\n")]
    lines = [ln.strip() for ln in block.strip().splitlines()]
    assert lines == [
        "embryo2)",
        "sha256sum -c docs/S2_RULES.sha256",
        "$U run python experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact \\",
        '--expression "$EMBRYO_EXPRESSION" --out $O/builds/T2embryo_s2x.h5ad',
        'files="$O/builds/T2embryo_s2x.h5ad:T2:embryo:test_interp"',
    ], lines
    assert ("build_test.py embryo --method fgw_shape0 --positions floor_exact --expression mixture "
            "--out runs/builds/T2embryo_s2x.h5ad") in " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
