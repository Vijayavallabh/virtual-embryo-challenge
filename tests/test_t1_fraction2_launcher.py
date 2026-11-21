"""The test-board build launcher of row R20's T1 file (experiments/harness/TRIES.md; the captain's order of 2026-11-21,
docs/PLAN.md decision log 2026-11-21 (x); the last amendment of docs/S2_RULES.md).

The launcher experiments/builds/t1_fraction2.sh first rebuilds R13's file by the default command and candidate B by
the 1.5 command and stops unless each gives its known bytes; then it builds the 2.0 file (build_test.py with
--captain-ruling R20) twice from one clean commit and compares the two. These tests hold its command strings, hashes
and gate lines to that, check that it deletes nothing and reads no output of a test build, and run its refusal paths
for real under a temporary home (a wrong commit, a shortened commit, a reused run name, an odd run name). The launcher
is never run here with the real commit: that would build files."""

import os
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_fraction2.sh"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
R13_SHA = "ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f"  # set t1rc (row R13)
B_SHA = "c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062"  # candidate B (set t1cand_1420c70)


def test_the_launcher_builds_r13s_file_and_b_then_the_r20_file_twice():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^(BASE|B|X2)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"BASE": BASE, "B": f"{BASE} --feng-fraction-mult 1.5",
                    "X2": f"{BASE} --feng-fraction-mult 2.0 --captain-ruling R20"}
    assert "build D T1_s2rc $BASE_ARGS" in sh and "build B T1_s2rcB $B_ARGS" in sh
    assert 'build "$STEP" T1_s2rcB2 $X2_ARGS' in sh and "for STEP in X2a X2b" in sh
    assert f"R13_SHA={R13_SHA}" in sh and '[ "$got" = "$R13_SHA" ]' in sh and "finish 20" in sh
    assert f"B_SHA={B_SHA}" in sh and '[ "$got" = "$B_SHA" ]' in sh and "finish 21" in sh
    assert sh.index("STEP=D") < sh.index("STEP=B") < sh.index("for STEP in X2a X2b") < sh.index("STEP=done")
    assert 'cmp -s "$O/X2a/T1_s2rcB2.h5ad" "$O/X2b/T1_s2rcB2.h5ad" || finish 31' in sh
    assert "sha256sum D/*.h5ad B/*.h5ad X2a/*.h5ad X2b/*.h5ad" in sh
    assert "vec check" in sh and "BOARD=T1:test" in sh
    assert "--seed" not in sh and "--board T1" not in sh  # seed and board: the builder's defaults, as the references
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh  # the builds, as the references
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh  # the tests
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=D")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert '[ "$head" = "$COMMIT" ]' in sh  # the full hash, never a prefix
    assert '[ "$st" = "$head True" ] || return 15' in sh  # a build counts only with a stamp of this commit, clean
    assert "PATH=$HOME/.local/bin:$PATH" in sh
    # every builder and check call sends its output to /dev/null: a test build is read by its exit code only
    for line in sh.splitlines():
        if re.search(r'build_test\.py t1|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    assert ".build.json" not in sh and "tests.out" in sh  # no sidecar record is read


def test_the_gate_lines():
    sh = LAUNCHER.read_text()
    assert "git status --porcelain --untracked-files=no" in sh and "finish 98" in sh
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in sh
    assert "grep -q '^| R20 ' experiments/harness/TRIES.md" in sh and "finish 95" in sh
    assert "grep -q 'row R20' docs/S2_RULES.md" in sh and "finish 94" in sh
    assert "[ -d data/derived/feng_relabel_c ]" in sh and "finish 92" in sh
    assert "O=runs/t1x2_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1x2/$RUN" in sh
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index("finish 96") < gate.index("finish 98") < gate.index("finish 97") < gate.index("finish 95")
    assert gate.index("finish 95") < gate.index("finish 94") < gate.index("finish 92") < gate.index("finish 93")
    # what the gate looks for is where it looks for it in this commit
    tries = (REPO / "experiments" / "harness" / "TRIES.md").read_text().splitlines()
    assert sum(ln.startswith("| R20 ") for ln in tries) == 1
    assert "row R20" in (REPO / "docs" / "S2_RULES.md").read_text()


def test_the_tests_step_runs_the_builders_targeted_tests():
    sh = LAUNCHER.read_text()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    assert set(listed) == {"tests/test_t1_fraction2_launcher.py", "tests/test_t1_fraction2_build.py",
                           "tests/test_t1_candidates_build.py", "tests/test_t1_bets.py", "tests/test_hidden_counts.py",
                           "tests/test_s2_rules.py"}
    assert all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in LAUNCHER.read_text().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", LAUNCHER.read_text()):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


def test_the_builder_holds_the_commands_the_launcher_runs():
    rules = " ".join((REPO / "docs" / "S2_RULES.md").read_text().split())
    want = f"build_test.py t1 {BASE} --feng-fraction-mult 2.0 --captain-ruling R20 --out runs/builds/T1_s2rcB2.h5ad"
    assert want in rules
    record = (REPO / "experiments" / "builds" / "results" / "t1cand_1420c70.txt").read_text()
    assert R13_SHA in record and B_SHA in record
    builder = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert "FENG_FRACTION_MULT_R20 = 2.0" in builder and 'T1_R20_ROW = "R20"' in builder


def _run(tmp_path: Path, *argv: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO, env=env, capture_output=True, text=True)


def test_a_wrong_commit_is_refused_at_the_gate_with_an_exit_file_and_nothing_built(tmp_path):
    out = _run(tmp_path, "0" * 40, "t1")
    scratch = tmp_path / "v" / "scratch" / "t1x2" / "t1"
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    assert (scratch / "gate.exit").read_text().strip() == "96"
    assert sorted(p.name for p in scratch.iterdir()) == ["gate.exit", "tmp"]  # no later step left a mark
    assert not list(REPO.glob("runs/t1x2_0000000_*"))  # the gate refused before any output folder
    again = _run(tmp_path, "0" * 40, "t1")  # the same run name: refused before any step, the first record untouched
    assert again.returncode == 91 and "never reused" in again.stderr
    assert (scratch / "gate.exit").read_text().strip() == "96"


def test_a_commit_prefix_is_refused_like_a_wrong_commit(tmp_path):
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    out = _run(tmp_path, head[:12], "t2")  # the right commit, shortened: never enough
    gate = tmp_path / "v" / "scratch" / "t1x2" / "t2" / "gate.exit"
    assert out.returncode == 96 and gate.read_text().strip() == "96"


def test_a_missing_or_odd_run_name_is_refused_before_any_folder_is_made(tmp_path):
    assert _run(tmp_path, "0" * 40).returncode != 0
    odd = _run(tmp_path, "0" * 40, "a/b")
    assert odd.returncode == 90 and not (tmp_path / "v").exists()
