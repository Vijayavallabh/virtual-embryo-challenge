"""The test-board build launcher of row R22's first search line (experiments/harness/TRIES.md; the captain's order of
2026-11-22, docs/PLAN.md decision log 2026-11-22 (j); the amendment of docs/S2_RULES.md "a score-guided search of round
steps may run").

experiments/builds/t1_fraction_line.sh first rebuilds R13's file by the default command, candidate B by the 1.5 command
and row R20's file by the 2.0 command, and stops unless each gives its known bytes; then it builds the 2.5 and the 3.0
file (build_test.py with --captain-ruling R22) twice each from one clean commit and compares each pair. These tests hold
its command strings, hashes and step order to that, check that it deletes nothing and reads no output of a test build,
and run its gate's refusal paths for real in throwaway repositories under a temporary home. The launcher is never run
here with a commit that passes its gate: that would build files. Nothing here reads the repository's own ledger or
rules file for its state of today."""

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_fraction_line.sh"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
R13_SHA = "ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f"  # set t1rc (row R13)
B_SHA = "c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062"  # candidate B (set t1cand_1420c70)
X2_SHA = "52897fe375b8a81afffa788592047c30361824148b00a8176514385060518837"  # row R20's file (t1x2_3d18841)
STEPS = ["gate", "tests", "D", "B", "X2", "X25a", "X25b", "X30a", "X30b", "done"]


def _sh() -> str:
    return LAUNCHER.read_text()


# ------------------------------------------------------------------ what the launcher runs


def test_the_launcher_reproduces_three_files_then_builds_each_line_file_twice():
    sh = _sh()
    args = dict(re.findall(r'^(BASE|B|X2|X25|X30)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"BASE": BASE, "B": f"{BASE} --feng-fraction-mult 1.5",
                    "X2": f"{BASE} --feng-fraction-mult 2.0 --captain-ruling R20",
                    "X25": f"{BASE} --feng-fraction-mult 2.5 --captain-ruling R22",
                    "X30": f"{BASE} --feng-fraction-mult 3.0 --captain-ruling R22"}
    assert "build D T1_s2rc $BASE_ARGS" in sh and "build B T1_s2rcB $B_ARGS" in sh
    assert "build X2 T1_s2rcB2 $X2_ARGS" in sh
    assert 'build "$STEP" T1_s2rcB25 $X25_ARGS' in sh and "for STEP in X25a X25b" in sh
    assert 'build "$STEP" T1_s2rcB30 $X30_ARGS' in sh and "for STEP in X30a X30b" in sh
    assert f"R13_SHA={R13_SHA}" in sh and '[ "$got" = "$R13_SHA" ]' in sh and "finish 20" in sh
    assert f"B_SHA={B_SHA}" in sh and '[ "$got" = "$B_SHA" ]' in sh and "finish 21" in sh
    assert f"X2_SHA={X2_SHA}" in sh and '[ "$got" = "$X2_SHA" ]' in sh and "finish 22" in sh
    order = [sh.index(s) for s in ("STEP=D", "STEP=B", "STEP=X2", "for STEP in X25a X25b", "for STEP in X30a X30b",
                                   "STEP=done")]
    assert order == sorted(order)
    assert 'cmp -s "$O/X25a/T1_s2rcB25.h5ad" "$O/X25b/T1_s2rcB25.h5ad" || finish 31' in sh
    assert 'cmp -s "$O/X30a/T1_s2rcB30.h5ad" "$O/X30b/T1_s2rcB30.h5ad" || finish 32' in sh
    assert ("sha256sum D/*.h5ad B/*.h5ad X2/*.h5ad X25a/*.h5ad X25b/*.h5ad X30a/*.h5ad X30b/*.h5ad" in sh
            and 'cp "$S/SHA256SUMS" "$O/SHA256SUMS"' in sh)
    assert "vec check" in sh and "BOARD=T1:test" in sh
    assert "--seed" not in sh and "--board T1" not in sh  # seed and board: the builder's defaults, as the references
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh  # the builds, as the references
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh  # the tests
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=D")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and "PATH=$HOME/.local/bin:$PATH" in sh
    assert '[ "$head" = "$COMMIT" ]' in sh  # the full hash, never a prefix
    assert '[ "$st" = "$head True" ] || return 15' in sh  # a build counts only with a stamp of this commit, clean
    # every builder and check call sends its output to /dev/null: a test build is read by its exit code only
    for line in sh.splitlines():
        if re.search(r'build_test\.py t1|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    assert ".build.json" not in sh and "tests.out" in sh  # no sidecar record is read
    assert "upload" not in sh.replace("Nothing here uploads anything", "")


def test_the_known_hashes_are_the_recorded_ones():
    """The three identity hashes are the ones the earlier launcher and its run record hold (no build here)."""
    x2 = (REPO / "experiments" / "builds" / "results" / "t1x2_3d18841.SHA256SUMS").read_text().split()
    assert dict(zip(x2[1::2], x2[0::2])) == {"D/T1_s2rc.h5ad": R13_SHA, "B/T1_s2rcB.h5ad": B_SHA,
                                             "X2a/T1_s2rcB2.h5ad": X2_SHA, "X2b/T1_s2rcB2.h5ad": X2_SHA}
    old = (REPO / "experiments" / "builds" / "t1_fraction2.sh").read_text()
    assert f"R13_SHA={R13_SHA}" in old and f"B_SHA={B_SHA}" in old


def test_the_builder_takes_the_commands_the_launcher_runs():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    assert bt.FENG_FRACTION_MULTS_R22 == (2.5, 3.0) and bt.T1_R22_ROW == "R22"
    assert bt.FENG_FRACTION_MULT_R20 == 2.0 and bt.T1_R20_ROW == "R20"
    sh = _sh()
    args = dict(re.findall(r'^(X2|X25|X30)_ARGS="([^"]+)"$', sh, flags=re.M))
    for key, (mult, row) in {"X2": (2.0, "R20"), "X25": (2.5, "R22"), "X30": (3.0, "R22")}.items():
        a = bt.build_parser().parse_args(["t1", *args[key].split(), "--out", "x.h5ad"])
        assert a.feng_fraction_mult == mult and a.captain_ruling == row
        assert bt.FENG_FRACTION_MULTS_RULED[mult][0] == row


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    assert "grep -q '^| R22 ' experiments/harness/TRIES.md" in gate
    assert "grep -q '^| R20 ' experiments/harness/TRIES.md" in gate
    assert "grep -q 'row R22' docs/S2_RULES.md" in gate
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "O=runs/t1line_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1line/$RUN" in sh
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 94, 92, 93]


def test_the_tests_step_runs_the_builders_and_launchers_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    assert set(listed) == {"tests/test_t1_fraction_line_launcher.py", "tests/test_t1_fraction_line_build.py",
                           "tests/test_t1_fraction2_build.py", "tests/test_t1_fraction2_launcher.py",
                           "tests/test_t1_candidates_build.py", "tests/test_t1_bets.py", "tests/test_hidden_counts.py",
                           "tests/test_s2_rules.py"}
    assert all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in _sh().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 2
    for msg in usages:  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


# ------------------------------------------------------------------ the gate's refusals, run for real


ROW_R20 = "| R20 (the captain's order of a day; ruling, no run) | T1 | a file |\n"
ROW_R22 = "| R22 (the captain's order of a later day; ruling, no run) | T1 | a line |\n"


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R20 + ROW_R22, rules: str = "the amendment names row R22\n",
          relabel: bool = True, good_sha: bool = True) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher and the files its gate reads; returns it, the env and its HEAD."""
    env = _env(tmp_path)
    (tmp_path / "home").mkdir(exist_ok=True)
    repo = tmp_path / "repo"
    (repo / "experiments" / "builds").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "docs").mkdir()
    (repo / "experiments" / "builds" / "t1_fraction_line.sh").write_text(_sh())
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / "docs" / "S2_RULES.md").write_text(rules)
    sha = hashlib.sha256(rules.encode()).hexdigest() if good_sha else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    (repo / ".gitignore").write_text("data/\nruns/\n")
    if relabel:
        (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "experiments/builds/t1_fraction_line.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1line" / run


def _stopped_at_the_gate(tmp_path: Path, repo: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no later step left a mark
    assert not list((s / "tmp").iterdir())


def test_a_wrong_commit_and_a_prefix_are_refused_and_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1")
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, repo, "r1", 96)
    assert not (repo / "runs").exists()  # refused before any output folder
    again = _run(repo, env, head, "r1")  # the same run name, now with the right commit: refused before any step
    assert again.returncode == 91 and "never reused" in again.stderr
    _stopped_at_the_gate(tmp_path, repo, "r1", 96)
    short = _run(repo, env, head[:12], "r2")  # the right commit, shortened: never enough
    assert short.returncode == 96
    _stopped_at_the_gate(tmp_path, repo, "r2", 96)


def test_a_missing_or_odd_run_name_is_refused_before_any_folder_is_made(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head).returncode != 0
    assert _run(repo, env).returncode != 0
    odd = _run(repo, env, head, "a/b")
    assert odd.returncode == 90 and not (tmp_path / "home" / "v").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R20}{ROW_R22}more\n")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, repo, "r1", 98)


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ROW_R20}, 95, "row R22 is not in TRIES.md"),
    ({"tries": ROW_R22}, 95, "row R20 is not in TRIES.md"),
    ({"tries": ROW_R20 + "| R22b (a row with another ID) |\n"}, 95, "row R22 is not in TRIES.md"),
    ({"rules": "the amendment names row R20 only\n"}, 94, "does not name row R22"),
    ({"relabel": False}, 92, "feng_relabel_c is not in this checkout"),
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1")
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, repo, "r1", code)
    assert not (repo / "runs" / f"t1line_{head[:7]}_r1").exists()


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t1line_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, repo, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]
