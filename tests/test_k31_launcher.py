"""The launcher of row R72's G0 and try K31's read (experiments/t3_decomp/run_k31.sh): its text, and its refusal paths
executed in a throwaway git repository under the test's own temporary folder, with a hand-made ledger (rows R72 and
K31), a hand-made manifest, dummy K18 task files and a `uv` stub that stands in for every Python call (nothing is
built, scored or read). The identity build can never give C2's recorded bytes here, so a stubbed call that passes
every gate stops at the identity step. No competition file, no committed ledger, no remote machine."""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
SH = REPO / "experiments" / "t3_decomp" / "run_k31.sh"
SH_REL = "experiments/t3_decomp/run_k31.sh"
MAN_REL = "experiments/t3_decomp/results/t3_cube/MANIFEST.json"
K18_REL = "experiments/t3_decomp/results/k18_built_files/tasks"
RUN_DIR = "v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
VERTS = ("ppn", "pnp", "pnn", "npp", "npn", "nnp", "nnn")
MAN_TEXT = '{"tables": "a hand-made manifest"}\n'
MAN_SHA = hashlib.sha256(MAN_TEXT.encode()).hexdigest()
R72 = ("| R72 (the captain's order of 2026-11-25, ruling, no run: \"Yes ... each only with `--captain-ruling R72` ...: "
       + ", ".join(f"T3_cube_{v}.h5ad `--form c2-s2 --vertex {v}`" for v in VERTS) + "\"; not a try) | T3 | x |\n")
K31 = (f"| K31 (one row, one script `experiments/t3_decomp/k31_cube_read.py`; cube tables manifest sha256 `{MAN_SHA}`;"
       " four tries) | T3 | the row's other cells |\n")
LEDGER = "| K30 (another row) | T3 | x |\n" + R72 + K31
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES")

UV_STUB = """#!{py} -I
# uv run --frozen --no-sync python -I -B ARGS, or uv run --frozen --no-sync vec check ...: the tests exit
# $HOME/tests_rc; a -c one-liner runs in python -I -B; `k31_cube_read.py row` prints the sha256 the row records (read
# from the ledger as the script does) and exits 2 when it holds none; the identity build exits $HOME/build_rc or
# writes a file that is not C2's bytes and a stamp of $HOME/stamp (default: HEAD clean); `vec check` exits
# $HOME/check_rc; every call is logged in $HOME/calls
import json, os, re, subprocess, sys
from pathlib import Path
a = sys.argv[1:]
home = Path(os.environ["HOME"])
def code(name):
    p = home / name
    return int(p.read_text()) if p.exists() else 0
with open(home / "calls", "a") as fh:
    fh.write(" ".join(a) + "\\n")
if a[:4] == ["run", "--frozen", "--no-sync", "vec"]:
    sys.exit(code("check_rc"))
if a[:6] != ["run", "--frozen", "--no-sync", "python", "-I", "-B"]:
    sys.exit(97)
a = a[6:]
if a[0] == "-m":
    sys.exit(code("tests_rc"))
if a[0] == "-c":
    os.execv(sys.executable, [sys.executable, "-I", "-B", *a])
if a[0].endswith("k31_cube_read.py") and a[1] == "row":
    lines = [ln for ln in Path("experiments/harness/TRIES.md").read_text().splitlines() if ln.startswith("| K31 ")]
    shas = set(re.findall(r"cube tables manifest sha256 `([0-9a-f]{64})`", lines[0])) if len(lines) == 1 else set()
    if len(shas) != 1:
        sys.exit(2)
    print(shas.pop())
    sys.exit(0)
if a[0].endswith("build_two_knockout.py"):
    if code("build_rc"):
        sys.exit(code("build_rc"))
    out = Path(a[a.index("--out") + 1])
    out.write_bytes(b"not C2's bytes")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    st = (home / "stamp").read_text().split() if (home / "stamp").exists() else [head, "true"]
    Path(str(out) + ".stamp.json").write_text(json.dumps({"commit": st[0], "clean": st[1] == "true"}))
    sys.exit(0)
sys.exit(code(Path(a[0]).stem + "_" + a[1] + "_rc"))
"""


def _git_env() -> dict:
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                        "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, env=_git_env(),
                       text=True)
    return r.stdout.strip()


def _repo(tmp_path: Path, ledger: str = LEDGER, man: str = MAN_TEXT) -> SimpleNamespace:
    repo, home, bindir = tmp_path / "repo", tmp_path / "home", tmp_path / "bin"
    (home / RUN_DIR).mkdir(parents=True)
    files = {SH_REL: SH.read_text(), "experiments/harness/TRIES.md": ledger, ".gitignore": "runs/\n", MAN_REL: man,
             **{f"{K18_REL}/{f}": "{}\n" for f in ("plain_g1.json", "plain_g2.json")}}
    repo.mkdir()
    _git(repo, "init", "-q")
    for rel, t in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t)
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "init")
    bindir.mkdir()
    (bindir / "uv").write_text(UV_STUB.replace("{py}", sys.executable))
    (bindir / "uv").chmod(0o755)
    return SimpleNamespace(repo=repo, home=home, bin=bindir, head=_git(repo, "rev-parse", "HEAD"))


def _launch(L, *args) -> tuple[int, str]:
    env = {**_git_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}",
           "GIT_CEILING_DIRECTORIES": str(L.repo.parent)}
    r = subprocess.run(["bash", SH_REL, *args], cwd=L.repo, env=env, capture_output=True, text=True, timeout=120)
    return r.returncode, r.stderr


def _exit(L, run: str, step: str) -> str:
    p = L.home / "v" / "scratch" / "k31" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else "absent"


def _calls(L) -> list[str]:
    p = L.home / "calls"
    return p.read_text().splitlines() if p.exists() else []


def test_the_launcher_text():
    sh = SH.read_text()
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    assert not re.search(r"\$\{[^}]*:\?[^}]*'[^}]*\}", sh)  # no apostrophe inside ${var:?message}
    assert not re.search(r"(^|[\s;|&(])rm(\s|$)", sh, flags=re.M)  # nothing deleted
    threads = [int(x) for x in re.findall(r"threads (\d+)", sh)]
    assert threads and max(threads) <= 16
    order = [sh.index(f"begin {s} ") for s in ("preflight", "tests", "identity", "g0", "dryrun", "run", "read")]
    assert order == sorted(order)
    for step, prev in (("tests", "preflight"), ("identity", "tests"), ("g0", "identity"), ("dryrun", "g0"),
                       ("run", "dryrun"), ("read", "run")):
        assert f"begin {step} {prev}" in sh
    for t in ("tests/test_sign_cube.py", "tests/test_cube_tables.py", "tests/test_k31_cube_read.py",
              "tests/test_k31_launcher.py", "tests/test_k18_built_files.py"):
        assert t in sh and (REPO / t).exists()
    c2 = (REPO / "experiments/builds/results/t3two_fb01487.SHA256SUMS").read_text()
    sha = re.search(r"^C2_SHA=([0-9a-f]{64})", sh, flags=re.M).group(1)
    assert f"{sha}  c2-s2.h5ad" in c2
    assert 'S="$HOME/v/scratch/k31/$RUN"' in sh and 'R="runs/k31_$RUN"' in sh and f"RUN_DIR=$HOME/{RUN_DIR}" in sh
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh and "set -u" in sh
    assert "build_two_knockout.py build \"$DEC\" --form c2-s2 " in sh  # C2's own command
    assert 'names "$first" "--form c2-s2 --vertex $v"' in sh
    assert "build_two_knockout.py build" in sh and "THE TRIES COUNT FROM THE START OF THIS STEP" in sh


def test_usage_and_the_run_name(tmp_path):
    L = _repo(tmp_path)
    for args in ((), (L.head,), (L.head, "r", "x"), ("", "r"), (L.head, "")):
        assert _launch(L, *args)[0] == 2, args
    for bad in (".hidden", "a/b", "x y"):
        assert _launch(L, L.head, bad)[0] == 90
    assert _calls(L) == [] and not (L.home / "v" / "scratch").exists()


def test_head_the_run_folder_and_the_tree(tmp_path):
    L = _repo(tmp_path)
    code, err = _launch(L, "0" * 40, "r1")
    assert code == 97 and "HEAD is not" in err and _exit(L, "r1", "preflight") == "97"
    assert _exit(L, "r1", "done") == "97" and _exit(L, "r1", "tests") == "absent"
    code, err = _launch(L, L.head, "r1")
    assert code == 91 and "another commit" in err
    (L.repo / "experiments/harness/TRIES.md").write_text(LEDGER + "an edit\n")
    code, err = _launch(L, L.head, "r2")
    assert code == 98 and "not clean" in err and _exit(L, "r2", "preflight") == "98"
    assert _calls(L) == []


@pytest.mark.parametrize("ledger,code,why", [
    (K31, 95, "exactly one row R72"),
    (R72 * 2 + K31, 95, "exactly one row R72"),
    (R72.replace("R72 (the captain's order", "R72 (a ruling") + K31, 95, "not the captain's order"),
    (R72.replace("ruling, no run", "ruled") + K31, 95, "ruling, no run"),
    (R72.replace("`--captain-ruling R72`", "--captain-ruling R72") + K31, 95, "--captain-ruling R72"),
    (R72.replace("`--form c2-s2 --vertex nnp`", "--form c2-s2 --vertex nnp") + K31, 95, "vertex nnp"),
    (R72.replace("`--form c2-s2 --vertex pnn`", "`--form c2-s2 --vertex pnn --x`") + K31, 95, "vertex pnn"),
    (R72, 94, "exactly one row K31"),
    (R72 + K31.replace("k31_cube_read.py", "k29_cd7_read.py"), 94, "exactly one row K31"),
    (R72 + K31.replace(MAN_SHA, "b" * 64), 94, "manifest sha256"),
    (R72 + K31.replace(f"`{MAN_SHA}`", MAN_SHA), 94, "manifest sha256"),
])
def test_the_preflight_refuses_rows_r72_and_k31_not_as_required(tmp_path, ledger, code, why):
    L = _repo(tmp_path, ledger)
    got, err = _launch(L, L.head, "r")
    assert got == code and why in err and _exit(L, "r", "preflight") == str(code), err
    assert not any("build_two_knockout" in c or "cube_tables.py g0" in c for c in _calls(L))


def test_a_failing_test_step_stops_before_any_build(tmp_path):
    L = _repo(tmp_path)
    (L.home / "tests_rc").write_text("1")
    code, _ = _launch(L, L.head, "t")
    assert code == 10 and _exit(L, "t", "tests") == "10" and _exit(L, "t", "identity") == "absent"
    assert not any("build_two_knockout" in c for c in _calls(L))


@pytest.mark.parametrize("setup,code", [
    ({"build_rc": "2"}, 12), ({"build_rc": "1"}, 11), ({"check_rc": "1"}, 13), ({"stamp": "f" * 40 + " true"}, 15),
    ({"stamp": "HEAD false"}, 15), ({}, 16),
])
def test_the_identity_rebuild_stops_on_anything_but_c2s_bytes(tmp_path, setup, code):
    L = _repo(tmp_path)
    for k, v in setup.items():
        (L.home / k).write_text(v.replace("HEAD", L.head))
    got, err = _launch(L, L.head, "i")
    assert got == code and _exit(L, "i", "identity") == str(code) and _exit(L, "i", "g0") == "absent", err
    assert not any("cube_tables.py g0" in c or "k31_cube_read.py run" in c for c in _calls(L))
