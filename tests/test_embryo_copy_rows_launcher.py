"""The test-board build launcher of the embryo's copy-rows file (practice row E50; docs/PLAN.md decision log 2026-11-21
(l) and (m); the specification experiments/coords/E50_COPY_ROWS_SPEC.md).

The launcher experiments/builds/embryo_copy_rows.sh runs only what the specification's start conditions name: the
default command must first give the exact-position file's known bytes, then the copy-rows file is built twice from one
clean commit and the two are compared. Its gate refuses unless row E50's committed reading is not void and says the
test file is built (reading (i) holds, mixture arm). These tests hold the launcher to that, check that it deletes
nothing and keeps no output of a test build, and run its refusal paths for real under a temporary home (a wrong commit,
a reused run name, an odd run name). The launcher is never run here with the real commit: that would build files."""

import json
import os
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "embryo_copy_rows.sh"
BASE = "--method fgw_shape0 --positions floor_exact --expression mixture"
COPY_ROWS = BASE + " --copy-rows --copy-rows-test-ok"
POSITIONS_SHA = "a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9"


def test_the_launcher_builds_the_positions_file_then_the_copy_rows_file_twice():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^(BASE|CR)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"BASE": BASE, "CR": COPY_ROWS}
    assert "build D T2embryo_s2x $BASE_ARGS" in sh
    assert 'build "$STEP" T2embryo_s2x_cr $CR_ARGS' in sh and "for STEP in CRa CRb" in sh
    assert f"POSITIONS_SHA={POSITIONS_SHA}" in sh and '[ "$got" = "$POSITIONS_SHA" ]' in sh
    assert 'cmp -s "$O/CRa/T2embryo_s2x_cr.h5ad" "$O/CRb/T2embryo_s2x_cr.h5ad"' in sh
    assert "vec check" in sh and "BOARD=T2:embryo:test_interp" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "git status --porcelain" in sh and "S2_RULES.sha256" in sh and "'^| E50 '" in sh
    assert "READING=experiments/coords/results/e50_copy_rows/reading.json" in sh
    assert "PATH=$HOME/.local/bin:$PATH" in sh
    # every builder and check call sends its output to /dev/null: a test build is read by its exit code only
    for line in sh.splitlines():
        if re.search(r'build_test\.py embryo|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert "tests/test_embryo_copy_rows_launcher.py" in listed and all((REPO / f).exists() for f in listed)


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in LAUNCHER.read_text().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", LAUNCHER.read_text()):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


def test_the_command_is_the_specifications_and_the_builder_knows_both_options():
    spec = " ".join((REPO / "experiments" / "coords" / "E50_COPY_ROWS_SPEC.md").read_text().split())
    assert "--copy-rows-test-ok" in spec and "--expression mixture" in spec and POSITIONS_SHA in spec
    builder = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert '"--copy-rows"' in builder and '"--copy-rows-test-ok"' in builder
    build_all = " ".join((REPO / "experiments" / "builds" / "build_all.sh").read_text().split()).replace(" \\ ", " ")
    build_all = build_all.replace('"$EMBRYO_EXPRESSION"', "mixture")
    assert f"build_test.py embryo {BASE} --out $O/builds/T2embryo_s2x.h5ad" in build_all


def _reading_check() -> str:
    """The gate's one-line check of row E50's reading, as the launcher holds it."""
    m = re.search(r"python3 -I -B -c '([^']+)' \"\$READING\"", LAUNCHER.read_text())
    assert m, "the gate's reading check is not where the test expects it"
    return m.group(1)


def test_the_gate_takes_only_a_reading_that_is_not_void_holds_and_is_the_mixture_arms(tmp_path):
    good = {"void": False, "build_test_file": True, "i": {"holds": True}, "handed_over_expression": "mixture"}
    cases = [(good, 0), ({**good, "void": True}, 1), ({**good, "build_test_file": False}, 1),
             ({**good, "i": {"holds": False}}, 1), ({**good, "handed_over_expression": "series"}, 1),
             ({"void": True, "why": ["a gate miss"]}, 1), ({**good, "i": {}}, 1), ({**good, "void": None}, 1)]
    for k, (reading, want) in enumerate(cases):
        p = tmp_path / f"reading{k}.json"
        p.write_text(json.dumps(reading))
        got = subprocess.run(["python3", "-I", "-B", "-c", _reading_check(), str(p)], capture_output=True)
        assert got.returncode == want, (reading, got.returncode)


def _run(tmp_path: Path, *argv: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO, env=env, capture_output=True, text=True)


def test_a_wrong_commit_is_refused_at_the_gate_with_an_exit_file_and_nothing_built(tmp_path):
    out = _run(tmp_path, "0" * 40, "t1")
    scratch = tmp_path / "v" / "scratch" / "embcr" / "t1"
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
