"""Tests for experiments/t3_decomp/run_k26.sh (row K26's launcher): every refusal path and every step's exit mapping,
traced in a temporary git repository under a temporary HOME with a stub `uv` whose exit status the test sets. No
Python step of the try runs; no data of any kind is read."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]
SH = REPO / "experiments" / "t3_decomp" / "run_k26.sh"
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES", "GIT_CEILING_DIRECTORIES")
ROW = "| K26 (one row) | T3 | `experiments/t3_decomp/k26_query_map.py` |\n"
REGISTER = "".join(f"| {x} | source | ... | K26: per-gene units as map labels |\n" for x in ("X013", "X014", "X015",
                                                                                           "X017"))


# a stub uv: the gate check and the read take their own exit status; the guard.json one-liner prints STUB_PRINT
STUB = """#!/bin/sh
case " $* " in
  *" guard --runs "*) exit "${STUB_GUARD_RC:-${STUB_RC:-1}}" ;;
  *" read --runs "*) exit "${STUB_READ_RC:-${STUB_RC:-1}}" ;;
  *" -c "*) [ -n "${STUB_PRINT:-}" ] && echo "$STUB_PRINT"; exit 0 ;;
esac
exit "${STUB_RC:-1}"
"""


def _env():
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}


def _git(repo: Path, *args) -> str:
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                           "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, text=True,
                          env=_env()).stdout.strip()


def _setup(tmp_path: Path) -> SimpleNamespace:
    repo, home = tmp_path / "repo", tmp_path / "home"
    for d in ("v/vec-ext/t3/tables", "v/vec-ext/t3/tables_two_study", "v/vec-ext/t3/processed",
              "v/vec-ext/t3/reference"):
        (home / d).mkdir(parents=True)
    (home / "v/vec-ext/t3/reference/refstats.npz").write_text("x")
    files = {"experiments/t3_decomp/run_k26.sh": SH.read_text(), "experiments/harness/TRIES.md": ROW,
             "EXTERNAL_DATA.md": REGISTER, "experiments/t3_ext/TABLES_K24.sha256": "x\n"}
    repo.mkdir()
    _git(repo, "init", "-q")
    for rel, t in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t)
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "init")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "uv").write_text(STUB)
    (bindir / "uv").chmod(0o755)
    return SimpleNamespace(repo=repo, home=home, bin=bindir, head=_git(repo, "rev-parse", "HEAD"))


def _run(L, *args, rc=1, trace=False, **stub) -> tuple[int, str]:
    env = {**_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}", "STUB_RC": str(rc),
           "GIT_CEILING_DIRECTORIES": str(L.repo.parent), **{k: str(v) for k, v in stub.items()}}
    cmd = ["bash", *(["-x"] if trace else []), "experiments/t3_decomp/run_k26.sh", *args]
    r = subprocess.run(cmd, cwd=L.repo, env=env, capture_output=True, text=True, timeout=60)
    return r.returncode, r.stderr


def _exit(L, step, head=None) -> str:
    return (L.home / "v" / "scratch" / f"k26_{(head or L.head)[:7]}" / f"{step}.exit").read_text().strip()


def test_the_launcher_parses_and_never_removes_anything():
    sh = SH.read_text()
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    assert not re.search(r"(^|[\s;|&(])rm(\s|$)", sh) and "unlink" not in sh and "rmdir" not in sh
    assert "set -u" in sh and "finish()" in sh


def test_usage_and_unknown_steps(tmp_path):
    L = _setup(tmp_path)
    code, err = _run(L, trace=True)
    assert code == 2 and "usage: run_k26.sh" in err
    assert (L.home / "v/scratch/k26/usage.exit").read_text().strip() == "2"
    code, _ = _run(L, "bogus", L.head)
    assert code == 97 and _exit(L, "bogus") == "97"


def test_each_refusal_before_the_tests(tmp_path):
    L = _setup(tmp_path)
    code, err = _run(L, "tests", "0" * 40)
    assert code == 99 and "HEAD is not" in err and _exit(L, "tests", "0" * 40) == "99"
    (L.repo / "EXTERNAL_DATA.md").write_text("changed\n")
    code, err = _run(L, "tests", L.head)
    assert code == 98 and "not clean" in err and _exit(L, "tests") == "98"
    _git(L.repo, "checkout", "--", "EXTERNAL_DATA.md")
    # no row naming the script, and a row of another label
    for text in ("| K25 | `experiments/t3_decomp/k26_query_map.py` |\n", "| K26 | another script |\n"):
        (L.repo / "experiments/harness/TRIES.md").write_text(text)
        _git(L.repo, "commit", "-q", "-am", "row")
        h = _git(L.repo, "rev-parse", "HEAD")
        code, err = _run(L, "tests", h)
        assert code == 95 and "no TRIES.md row K26" in err and _exit(L, "tests", h) == "95"
        _git(L.repo, "revert", "--no-edit", "HEAD")
    h = _git(L.repo, "rev-parse", "HEAD")
    # an outside folder missing, then a symlink, then the reference statistics missing
    proc = L.home / "v/vec-ext/t3/processed"
    proc.rmdir()
    code, err = _run(L, "tests", h)
    assert code == 94 and "processed is not a real directory" in err
    (L.home / "elsewhere").mkdir()
    proc.symlink_to(L.home / "elsewhere")
    code, err = _run(L, "tests", h)
    assert code == 94 and "not a real directory" in err
    proc.unlink()
    proc.mkdir()
    rs = L.home / "v/vec-ext/t3/reference/refstats.npz"
    rs.unlink()
    code, err = _run(L, "tests", h)
    assert code == 94 and "refstats.npz is not a plain file" in err
    rs.write_text("x")
    # a register row that does not name K26
    (L.repo / "EXTERNAL_DATA.md").write_text(REGISTER.replace("| X014 | source | ... | K26: per-gene units as map "
                                                              "labels |", "| X014 | source | ... | K24 only |"))
    _git(L.repo, "commit", "-q", "-am", "register")
    h2 = _git(L.repo, "rev-parse", "HEAD")
    code, err = _run(L, "tests", h2)
    assert code == 94 and "row X014 does not name K26" in err
    _git(L.repo, "revert", "--no-edit", "HEAD")
    # the listing untracked
    _git(L.repo, "rm", "-q", "--cached", "experiments/t3_ext/TABLES_K24.sha256")
    _git(L.repo, "commit", "-q", "-m", "untrack")
    h3 = _git(L.repo, "rev-parse", "HEAD")
    code, err = _run(L, "tests", h3)
    assert code == 94 and "TABLES_K24.sha256 is not tracked" in err


def test_each_step_maps_the_scripts_exit_and_refuses_without_the_step_before(tmp_path):
    L = _setup(tmp_path)
    code, err = _run(L, "features", L.head, rc=0)
    assert code == 96 and "step tests has not exited 0" in err and _exit(L, "features") == "96"
    assert _run(L, "tests", L.head, rc=1)[0] == 10 and _exit(L, "tests") == "10"
    assert _run(L, "tests", L.head, rc=0)[0] == 0
    assert (L.home / "v/scratch" / f"k26_{L.head[:7]}" / "tests.head").read_text().strip() == L.head
    for rc, want in ((2, 22), (4, 24), (1, 20), (0, 0)):
        assert _run(L, "features", L.head, rc=rc)[0] == want and _exit(L, "features") == str(want)
    for rc, want in ((7, 27), (2, 22), (1, 20)):
        assert _run(L, "stage0", L.head, rc=rc)[0] == want
        assert _run(L, "pregate", L.head, rc=0)[0] == 96  # stage 0 did not exit 0: nothing after it runs
    assert _run(L, "stage0", L.head, rc=0)[0] == 0
    for rc, want in ((5, 25), (7, 27), (2, 22), (4, 24), (1, 20)):
        assert _run(L, "pregate", L.head, rc=rc)[0] == want
        assert _run(L, "guard", L.head, rc=0)[0] == 96
    assert _run(L, "pregate", L.head, rc=0)[0] == 0
    assert _run(L, "guard", L.head, rc=1)[0] == 30 and _exit(L, "guard") == "30"
    R = L.repo / f"runs/k26_{L.head[:7]}"
    parts = sorted(p.name for p in (R / "exit").iterdir() if p.name.startswith("guard_") and "." not in p.name)
    assert parts == sorted(f"guard_{p}" for p in ("plain_g1_0-4", "plain_g1_5-9", "plain_g2_0-4", "plain_g2_5-9",
                                                  "dm", "x2_x2_0-4", "x2_x2_5-9"))
    assert _run(L, "run", L.head, rc=0)[0] == 96
    assert _run(L, "guard", L.head, rc=0)[0] == 0  # every part and the gate exit 0
    assert _run(L, "read", L.head, rc=0)[0] == 96  # run has not exited 0
    assert _run(L, "run", L.head, rc=1)[0] == 30
    assert _run(L, "run", L.head, rc=0)[0] == 0
    for rc, want in ((3, 43), (4, 44), (1, 41), (0, 0)):
        assert _run(L, "read", L.head, rc=rc)[0] == want and _exit(L, "read") == str(want)
    # a step exited 0 at another commit does not open the next one
    (L.repo / "x.txt").write_text("x")
    _git(L.repo, "add", "x.txt")
    _git(L.repo, "commit", "-q", "-m", "next")
    h = _git(L.repo, "rev-parse", "HEAD")
    assert _run(L, "features", h, rc=0)[0] == 96


def test_a_failed_gate_writes_the_void_reading_once(tmp_path):
    L = _setup(tmp_path)
    for step in ("tests", "features", "stage0", "pregate"):
        assert _run(L, step, L.head, rc=0)[0] == 0
    S = L.home / "v/scratch" / f"k26_{L.head[:7]}"
    # the parts exit 0 and the gate exits 6: the VOID read runs (exit 0) and the step finishes 26
    assert _run(L, "guard", L.head, rc=0, STUB_GUARD_RC=6, STUB_READ_RC=0)[0] == 26
    assert (S / "void_read.exit").read_text().strip() == "0" and _exit(L, "guard") == "26"
    assert _run(L, "run", L.head, rc=0)[0] == 96  # no scoring part after a failed gate
    # a failing gate whose void read crashed or was incomplete
    assert _run(L, "guard", L.head, rc=0, STUB_GUARD_RC=6, STUB_READ_RC=3)[0] == 43
    assert _run(L, "guard", L.head, rc=0, STUB_GUARD_RC=6, STUB_READ_RC=1)[0] == 41
    # a rerun over a recorded guard.json: pass -> 0; fail with a reading already there -> 26 without a second read;
    # unreadable -> 20
    R = L.repo / f"runs/k26_{L.head[:7]}"
    (R / "guard.json").write_text("{}")
    assert _run(L, "guard", L.head, rc=1, STUB_PRINT="pass")[0] == 0
    out = L.repo / "experiments/t3_decomp/results/k26_query_map"
    out.mkdir(parents=True)
    (out / "reading.json").write_text("{}")
    assert _run(L, "guard", L.head, rc=1, STUB_PRINT="fail", STUB_READ_RC=1)[0] == 26
    assert _run(L, "guard", L.head, rc=1)[0] == 20
