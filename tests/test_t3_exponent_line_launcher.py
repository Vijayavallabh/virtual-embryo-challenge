"""The launcher of row R22's T3 search line: experiments/t3_decomp/build_outside_form.sh with --table-exponent 4.0
--captain-ruling R22 HEAD s2 xall. It gates on the barred s1 and s2 by identity and on row K23's file on s2 by its own
command (2.0, no ruling), then builds the 4.0 file twice. These tests hold its strings and codes to that, check that it
deletes nothing and reads no output of a test build, and run its refusal paths for real in throwaway repositories
under a temporary home with a stand-in uv. No test lets a call reach a build. Nothing here reads the repository's own
ledger for its state of today."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "t3_decomp" / "build_outside_form.sh"
K23_S2 = "77ff0b947f7f409c9ce171bdc59a44de3d3c2ef2f3e5a7332011a459ce9ee9c3"
ROW = ("| R22 (the captain's order of a day; ruling, no run, not a try; the builder is to take `--feng-fraction-mult "
       "2.5` only with `--captain-ruling R22`; and the T3 outside-form builder `--table-exponent 4.0`, each only with "
       "`--captain-ruling R22`) | T1 | x |\n")


def _sh() -> str:
    return LAUNCHER.read_text()


# ------------------------------------------------------------------ what the launcher runs


def test_the_ruled_mode_s_strings():
    sh = _sh()
    assert 'XFLAG=(--table-exponent "$EXP" --captain-ruling "$RULING")' in sh
    assert 'if [ "$RULING" != R22 ] || [ "$EXP" != 4.0 ] || [ "$BASE" != s2 ]; then' in sh
    assert f"K23_S2={K23_S2}" in sh
    assert 'one xall-s2 "$f" --table-exponent 2.0 || rc=$?' in sh and 'if [ "$h" != "$K23_S2" ]; then finish 52' in sh
    assert "if [ $b = 5 ]; then return 5; fi" in sh and "if [ $rc = 5 ]; then finish 47; fi" in sh
    assert "XTESTS+=(tests/test_t3_exponent_line_build.py tests/test_t3_exponent_line_launcher.py)" in sh
    assert all((REPO / f).exists() for f in ("tests/test_t3_exponent_line_build.py",
                                             "tests/test_t3_exponent_line_launcher.py"))
    order = [sh.index(s) for s in ("gate s1 20", 'gate "$BASE" 30', "one xall-s2", "for form in $CANDS")]
    assert order == sorted(order)
    # the checks, then the word running, then the head file
    assert sh.index("r22_row ||") < sh.index("echo running > \"$S/all.exit\"") < sh.index("git rev-parse HEAD > ")
    assert '[ "$(git rev-parse HEAD)" = "$HEAD_SHA" ]' in sh  # the full hash
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and "threads 16" in sh


def test_the_known_hash_is_the_build_record_s():
    rec = (REPO / "experiments" / "builds" / "results" / "t3k23_49f3fd5.SHA256SUMS").read_text().split()
    assert dict(zip(rec[1::2], rec[0::2], strict=True))["xall-s2.h5ad"] == K23_S2
    assert "[s2]=1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f" in _sh()
    assert "[s1]=3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2" in _sh()


def test_build_and_check_output_goes_to_dev_null():
    lines = _sh().splitlines()
    for i, line in enumerate(lines):
        if "build_outside_form.py build" in line:
            assert "> /dev/null 2>&1" in lines[i + 1], line
        if re.search(r"vec check \"", line):
            assert "> /dev/null 2>&1" in line, line
    assert ".json\"" not in _sh().replace('"$1.stamp.json"', "")  # no sidecar record is read, only the stamp


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in _sh().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 5
    for msg in usages:
        assert "'" not in msg and '"' not in msg, msg


# ------------------------------------------------------------------ refusals, run for real


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    stand = tmp_path / "bin"
    stand.mkdir(exist_ok=True)
    uv = stand / "uv"  # a stand-in: records what the run's exit and head files hold when it is first called, fails
    uv.write_text('#!/usr/bin/env bash\nfor d in "$HOME"/v/scratch/t3outside_*; do\n'
                  '  { cat "$d/all.exit"; cat "$d/head"; } > "$HOME/seen" 2> /dev/null\ndone\nexit 1\n')
    uv.chmod(0o755)
    env["PATH"] = f"{stand}:{env.get('PATH', '/usr/bin:/bin')}"
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW, tables: bool = True) -> tuple[Path, dict, str]:
    env = _env(tmp_path)
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    if tables:
        (home / "v" / "vec-ext" / "t3" / "tables").mkdir(parents=True)
    repo = tmp_path / "repo"
    (repo / "experiments" / "t3_decomp").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "experiments" / "t3_decomp" / "build_outside_form.sh").write_text(_sh())
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / ".gitignore").write_text("runs/\n")
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str, cwd: Path | None = None,
         script: str = "experiments/t3_decomp/build_outside_form.sh") -> subprocess.CompletedProcess:
    return subprocess.run(["bash", script, *argv], cwd=cwd or repo, env=env, capture_output=True, text=True)


def _ruled(head: str, *, exp="4.0", ruling="R22", base="s2", kind="xall") -> list[str]:
    return ["--table-exponent", exp, "--captain-ruling", ruling, head, base, kind]


def _s(tmp_path: Path, head: str, base: str = "s2", exp: str = "4.0") -> Path:
    return tmp_path / "home" / "v" / "scratch" / f"t3outside_{head[:7]}_{base}_e{exp}"


def _refused(tmp_path: Path, head: str, code: int, **kw) -> None:
    s = _s(tmp_path, head, **kw)
    assert (s / "refused.exit").read_text().strip() == str(code)
    assert not (s / "all.exit").exists() and not (s / "head").exists()
    assert not (tmp_path / "home" / "seen").exists()  # nothing ran


def test_a_wrong_commit_or_a_prefix_is_refused_without_touching_all_exit(tmp_path):
    repo, env, head = _repo(tmp_path)
    for given in ("0" * 40, head[:12]):
        out = _run(repo, env, *_ruled(given))
        assert out.returncode == 97 and "HEAD is not" in out.stderr
    _refused(tmp_path, "0" * 40, 97)
    _refused(tmp_path, head[:12], 97)
    assert not (repo / "runs").exists()


@pytest.mark.parametrize("argv", [
    dict(ruling="R20"), dict(ruling="R23"), dict(exp="2.0"), dict(base="s1"), dict(base="s1x2.0"), dict(kind="xunrel"),
])
def test_other_rulings_exponents_bases_and_kinds_are_refused(tmp_path, argv):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, *_ruled(head, **argv))
    assert out.returncode == 2
    s = _s(tmp_path, head, base=argv.get("base", "s2"), exp=argv.get("exp", "4.0"))
    assert (s / "refused.exit").read_text().strip() == "2" and not (s / "all.exit").exists()
    assert not (repo / "runs").exists()


def test_two_kinds_are_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, *_ruled(head), "xunrel").returncode == 2
    _refused(tmp_path, head, 2)


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW}more\n")
    out = _run(repo, env, *_ruled(head))
    assert out.returncode == 98
    _refused(tmp_path, head, 98)


def test_a_missing_tables_folder_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path, tables=False)
    assert _run(repo, env, *_ruled(head)).returncode == 96
    _refused(tmp_path, head, 96)


@pytest.mark.parametrize("tries", [
    "",
    ROW + ROW,
    ROW.replace("`--table-exponent 4.0`", "--table-exponent 4.0"),
    ROW.replace("(the captain's order", "(the captain's ruling"),
    ROW.replace("ruling, no run", "a run"),
    ROW.replace("`--captain-ruling R22`", "the flag"),
    ROW.replace("; and the T3 outside-form builder `--table-exponent 4.0`, each only with `--captain-ruling R22`) | T1 | x",
                ") | T1 | `--table-exponent 4.0`"),
    ROW.replace("| R22 (", "| R22b ("),
])
def test_a_ledger_row_that_does_not_name_the_option_is_refused(tmp_path, tries):
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, *_ruled(head))
    assert out.returncode == 94 and "row R22" in out.stderr
    _refused(tmp_path, head, 94)
    assert not (repo / "runs").exists()


def test_a_call_away_from_the_root_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, *_ruled(head), cwd=repo / "experiments", script="t3_decomp/build_outside_form.sh")
    assert out.returncode == 93 and "root" in out.stderr
    _refused(tmp_path, head, 93)


def test_a_finished_run_s_exit_is_never_touched_by_a_refused_call(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t3_outside_form_{head[:7]}_s2_e4.0"
    o.mkdir(parents=True)
    s = _s(tmp_path, head)
    s.mkdir(parents=True)
    (s / "all.exit").write_text("0\n")
    (s / "head").write_text(f"{head}\n")
    out = _run(repo, env, *_ruled(head))
    assert out.returncode == 95 and "exists" in out.stderr
    assert (s / "all.exit").read_text() == "0\n" and (s / "refused.exit").read_text().strip() == "95"
    assert sorted(p.name for p in o.iterdir()) == []
    assert _run(repo, env, *_ruled("0" * 40)).returncode == 97 and (s / "all.exit").read_text() == "0\n"


def test_a_passing_gate_writes_running_then_the_head_then_runs(tmp_path):
    """Every check passes: all.exit holds the word running and the head file the commit when the first step (the
    targeted tests, through the stand-in uv, which fails) runs; the run then stops with 10, nothing built."""
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, *_ruled(head))
    assert out.returncode == 10, out.stderr
    assert (tmp_path / "home" / "seen").read_text().split() == ["running", head]
    s = _s(tmp_path, head)
    assert (s / "all.exit").read_text().strip() == "10" and (s / "tests.exit").read_text().strip() != "0"
    assert not (s / "refused.exit").exists()
    o = repo / "runs" / f"t3_outside_form_{head[:7]}_s2_e4.0"
    assert sorted(p.name for p in o.iterdir()) == ["a", "b", "gate"] and not any(o.rglob("*.h5ad"))


def test_without_the_ruling_a_refusal_writes_all_exit_as_before(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "--table-exponent", "4.0", "0" * 40, "s2", "xall")
    assert out.returncode == 97
    s = _s(tmp_path, "0" * 40)
    assert (s / "all.exit").read_text().strip() == "97" and not (s / "refused.exit").exists()
