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
    # the first upload's line, the second upload's as it was frozen, the second's since the amendment of 2026-11-18,
    # and the last upload's since the amendment of 2026-11-19 for the heart (the full trim, set heart3)
    assert commands == ["build_test.py heart --entry s1 --trend",
                        "build_test.py heart --entry s2_trim --trend",
                        "build_test.py heart --entry s1 --expression copy --z-mode raw",
                        "build_test.py heart --entry core_trim --expression copy --z-mode raw"], commands
    for c in commands:
        assert c in rules, f"build_all.sh runs `{c}`, which docs/S2_RULES.md does not give"


def test_the_heart3_set_checks_the_rules_hash_and_builds_the_full_trim_and_its_half():
    """The set heart3 checks docs/S2_RULES.sha256 before it builds, builds T2heart_s2core.h5ad by the full trim's
    command of the amendment of 2026-11-19 for the heart (the captain's own named exception), halves it with halve.py
    into T2heart_s2core_half.h5ad (the uploaded file) and checks both locally for T2:heart:test_extrap. The rules give
    both commands word for word with runs/ as the output folder."""
    text = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    block = text[text.index("\n  heart3)\n"):]
    block = block[:block.index("\n    ;;\n")]
    lines = [ln.strip() for ln in block.strip().splitlines()]
    assert lines == [
        "heart3)",
        "sha256sum -c docs/S2_RULES.sha256",
        "$U run python experiments/builds/build_test.py heart --entry core_trim --expression copy --z-mode raw "
        "--out $O/builds/T2heart_s2core.h5ad",
        "$U run python experiments/builds/halve.py $O/builds/T2heart_s2core.h5ad $O/builds/T2heart_s2core_half.h5ad",
        'files="$O/builds/T2heart_s2core.h5ad:T2:heart:test_extrap '
        '$O/builds/T2heart_s2core_half.h5ad:T2:heart:test_extrap"',
    ], lines
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
    assert ("uv run python experiments/builds/build_test.py heart --entry core_trim --expression copy --z-mode raw "
            "--out runs/builds/T2heart_s2core.h5ad") in rules
    assert ("uv run python experiments/builds/halve.py runs/builds/T2heart_s2core.h5ad "
            "runs/builds/T2heart_s2core_half.h5ad") in rules


def test_the_usage_line_and_the_unknown_set_message_name_every_set():
    text = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    assert "usage: build_all.sh s1|s2|heart2|t3alt|t1none|embryo2|heart3|embryo3|t1rc}" in text
    assert "unknown set $set_name (s1, s2, heart2, t3alt, t1none, embryo2, heart3, embryo3 or t1rc)" in text
    for name in ("s1", "s2", "heart2", "t3alt", "t1none", "embryo2", "heart3", "embryo3", "t1rc"):
        assert f"\n  {name})\n" in text, name


def test_the_t1_commands_of_build_all_are_the_ones_the_rules_give():
    """Each T1 build command of experiments/builds/build_all.sh (from `build_test.py t1` up to its --out) appears word
    for word in docs/S2_RULES.md, as the heart's do: the first upload's (row R6), the second as ruled on 2026-11-13
    (row R7), T1's last upload with no outside data, by the amendment of 2026-11-19 (set t1none, row R8), and the
    recomputed-groups file that replaced it by the later amendment of 2026-11-19 for T1 (set t1rc)."""
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
    script = (ROOT / "experiments" / "builds" / "build_all.sh").read_text().splitlines()
    commands = [line[line.index("build_test.py t1"):line.index(" --out")].strip()
                for line in script if "build_test.py t1" in line and " --out" in line]
    assert commands == ["build_test.py t1 --g-mom --base lib3",
                        "build_test.py t1 --g-mom --base mid --cm-share 0.60",
                        "build_test.py t1 --g-mom --base lib3 --feng-groups none",
                        "build_test.py t1 --g-mom --base lib3 --feng-groups recomputed"], commands
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


def test_the_t1rc_set_checks_the_rules_hash_and_builds_the_recomputed_groups_file():
    """The set t1rc checks docs/S2_RULES.sha256 before it builds, builds one file, T1_s2rc.h5ad, by the first upload's
    command plus --feng-groups recomputed (the later amendment of 2026-11-19 for T1; the captain's own say-so), and
    checks it locally for T1:test. The rules' section for it, the last of the file, gives the same command with runs/
    as its output folder, the captain's answer, the three sha256 it rests on and the condition on row T63's read."""
    text = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    block = text[text.index("\n  t1rc)\n"):]
    block = block[:block.index("\n    ;;\n")]
    lines = [ln.strip() for ln in block.strip().splitlines()]
    assert lines == [
        "t1rc)",
        "sha256sum -c docs/S2_RULES.sha256",
        "$U run python experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups recomputed "
        "--out $O/builds/T1_s2rc.h5ad",
        'files="$O/builds/T1_s2rc.h5ad:T1:test"',
    ], lines
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().split())
    heading = ("## Amendment of 2026-11-19: T1's last upload is the recomputed-groups file (the captain's own "
               "say-so)")
    assert rules.count(heading) == 1
    section = rules[rules.index(heading):]
    assert " ## " not in section  # the last section of the file
    assert ("uv run python experiments/builds/build_test.py t1 --g-mom --base lib3 --feng-groups recomputed "
            "--out runs/builds/T1_s2rc.h5ad") in section
    assert '"proceed with your recommendation"' in section  # the captain's answer, decision log 2026-11-19 (ii)
    for sha in ("883e74c2c1ef5f01d5244555f8b3dcae24255568c731a6d53ee2788b7b677b49",  # the labeller's frozen settings
                "36d1f7edcbbf3d7a80c68b51d5624a94bb753668079fe3d69aed720c71baa7e9",  # the recomputed statistics
                "ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f"):  # the check builds at 8138fa2
        assert sha in section, sha
    assert "**Built after row T63's read.**" in section


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
    # the first upload's line, the second upload's of ruling R4 (not uploaded), the last upload's of the first
    # amendment of 2026-11-19 for the embryo (set embryo2, not uploaded), and the series twin of the later one (embryo3)
    assert len(commands) == 4, commands
    for c in commands:
        assert c in rules, f"build_all.sh runs `{c}`, which docs/S2_RULES.md does not give"
    assert ("build_test.py embryo --method fgw_shape0 --positions floor_exact --expression mixture "
            "--out $O/builds/T2embryo_s2x.h5ad") in script
    assert ("build_test.py embryo --method fgw_shape0 --positions floor_exact --expression series "
            "--out $O/builds/T2embryo_s2xs.h5ad") in script


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


def test_the_embryo3_set_checks_the_rules_hash_and_builds_the_series_twin():
    """The set embryo3 checks docs/S2_RULES.sha256 before it builds, builds one file, T2embryo_s2xs.h5ad, by the one
    command of the later amendment of 2026-11-19 for the embryo (the series twin, try E47b; the captain's own say-so),
    and checks it locally for T2:embryo:test_interp. The rules give the same command with runs/ as its output folder."""
    text = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    block = text[text.index("\n  embryo3)\n"):]
    block = block[:block.index("\n    ;;\n")]
    lines = [ln.strip() for ln in block.strip().splitlines()]
    assert lines == [
        "embryo3)",
        "sha256sum -c docs/S2_RULES.sha256",
        "$U run python experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact \\",
        "--expression series --out $O/builds/T2embryo_s2xs.h5ad",
        'files="$O/builds/T2embryo_s2xs.h5ad:T2:embryo:test_interp"',
    ], lines
    assert ("uv run python experiments/builds/build_test.py embryo --method fgw_shape0 --positions floor_exact "
            "--expression series --out runs/builds/T2embryo_s2xs.h5ad") in " ".join(
        (ROOT / "docs" / "S2_RULES.md").read_text().split())


def test_t3s_last_file_is_t3alts_barred_s1_on_the_reference_list():
    """The amendment of 2026-11-19 for T3 names the barred branch's s1 with its sha256; set t3alt builds it and the
    reference list of the third rehearsal holds that sha256 under t3_build/barred_s1.h5ad. No new set builds it."""
    # the command's continuation lines joined, as for the embryo's commands of build_all.sh
    rules = " ".join((ROOT / "docs" / "S2_RULES.md").read_text().replace("\\\n", " ").split())
    sha = "3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2"
    section = rules[rules.index("## Amendment of 2026-11-19: T3's last upload is the barred branch's s1"):]
    assert sha in section
    assert ("--branch barred --submission s1 --wt data/raw/T2/E8.75.h5ad --run-dir "
            "runs/t3_decomp/gata4_20261028T122425Z --board T3:bcatenin --seed 1000001 "
            "--out runs/t3_build/barred_s1.h5ad") in section
    listing = (ROOT / "experiments" / "builds" / "results" / "builds_dc7b6ae.SHA256SUMS").read_text().split("\n")
    assert f"{sha}  t3_build/barred_s1.h5ad" in listing
    script = (ROOT / "experiments" / "builds" / "build_all.sh").read_text()
    block = script[script.index("\n  t3alt)\n"):]
    assert "t3 barred s1 $O/t3_build/barred_s1.h5ad\n" in block[:block.index("\n    ;;\n")]
