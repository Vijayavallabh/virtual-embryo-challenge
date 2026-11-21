"""Row R19's embryo overlay file (experiments/harness/TRIES.md; docs/S2_RULES.md, the amendment of 2026-11-21).

The file is the exact-position file's command (set embryo2 of experiments/builds/build_all.sh: `build_test.py embryo
--method fgw_shape0 --positions floor_exact --expression mixture`) with one named change, `--shape-interp overlay`,
and the test-board guard's flag `--shape-interp-test-ok`, which is the captain's yes of 2026-11-21. These tests hold the
launcher experiments/builds/embryo_overlay.sh to exactly that: the command it runs is the one the rules note gives, the
positions file must come out first with its known bytes, the overlay file is built twice, no output of a test build is
kept, and nothing in the launcher deletes anything. One test runs the launcher's refusal path for real (a wrong commit),
with its scratch folder under a temporary home, so the set-up lines and the exit file are exercised, not only read."""

import os
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "embryo_overlay.sh"
BASE = "--method fgw_shape0 --positions floor_exact --expression mixture"
OVERLAY = BASE + " --shape-interp overlay --shape-interp-test-ok"
POSITIONS_SHA = "a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9"


def _flat(path: Path) -> str:
    return " ".join(path.read_text().split())


def test_the_launcher_builds_the_positions_file_then_the_overlay_file_twice():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^(BASE|OV)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"BASE": BASE, "OV": OVERLAY}
    assert "build D T2embryo_s2x $BASE_ARGS" in sh
    assert 'build "$STEP" T2embryo_s2x_ov $OV_ARGS' in sh and "for STEP in OVa OVb" in sh
    assert f"POSITIONS_SHA={POSITIONS_SHA}" in sh and '[ "$got" = "$POSITIONS_SHA" ]' in sh
    assert 'cmp -s "$O/OVa/T2embryo_s2x_ov.h5ad" "$O/OVb/T2embryo_s2x_ov.h5ad"' in sh
    assert "vec check" in sh and "BOARD=T2:embryo:test_interp" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "git status --porcelain" in sh and "S2_RULES.sha256" in sh and "'^| R19 '" in sh
    assert "PATH=$HOME/.local/bin:$PATH" in sh
    # every builder and check call sends its output to /dev/null: a test build is read by its exit code only
    for line in sh.splitlines():
        if re.search(r'build_test\.py embryo|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert "tests/test_embryo_overlay.py" in listed and all((REPO / f).exists() for f in listed)


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in LAUNCHER.read_text().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", LAUNCHER.read_text()):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


def test_the_command_is_the_one_the_rules_note_gives_and_the_positions_command_is_set_embryo2s():
    rules = _flat(REPO / "docs" / "S2_RULES.md")
    assert f"experiments/builds/build_test.py embryo {OVERLAY} --out runs/builds/T2embryo_s2x_ov.h5ad" in rules
    assert "R19" in rules and "experiments/builds/embryo_overlay.sh" in rules
    build_all = _flat(REPO / "experiments" / "builds" / "build_all.sh").replace(" \\ ", " ")  # continuations joined
    build_all = build_all.replace('"$EMBRYO_EXPRESSION"', "mixture")
    assert f"build_test.py embryo {BASE} --out $O/builds/T2embryo_s2x.h5ad" in build_all
    ledger = (REPO / "experiments" / "harness" / "TRIES.md").read_text().splitlines()
    row = [x for x in ledger if x.startswith("| R19 ")]
    assert len(row) == 1 and OVERLAY in row[0] and '"Build it"' in row[0]


def _run(tmp_path: Path, *argv: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO, env=env, capture_output=True, text=True)


def test_a_wrong_commit_is_refused_at_the_gate_with_an_exit_file_and_nothing_built(tmp_path):
    out = _run(tmp_path, "0" * 40, "t1")
    scratch = tmp_path / "v" / "scratch" / "embov" / "t1"
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    assert (scratch / "gate.exit").read_text().strip() == "96"
    assert sorted(p.name for p in scratch.iterdir()) == ["gate.exit", "tmp"]  # no later step left a mark
    again = _run(tmp_path, "0" * 40, "t1")  # the same run name: refused before any step, the first record untouched
    assert again.returncode == 91 and "never reused" in again.stderr
    assert (scratch / "gate.exit").read_text().strip() == "96"


def test_a_missing_or_odd_run_name_is_refused_before_any_folder_is_made(tmp_path):
    assert _run(tmp_path, "0" * 40).returncode != 0
    odd = _run(tmp_path, "0" * 40, "a/b")
    assert odd.returncode == 90 and not (tmp_path / "v").exists()
