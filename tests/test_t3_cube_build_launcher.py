"""The launcher of row R72's cube builds (experiments/t3_decomp/t3_cube_build.sh): its text, and its refusal paths
executed in a throwaway git repository under the test's own temporary folder, with a hand-made ledger (rows R72 and
K31), a hand-made manifest and readings, and a `uv` stub that stands in for every Python call (nothing is built,
scored or read). The identity gate can never give C2's recorded bytes here, so a stubbed call that passes every
gate stops at the identity gate. No competition file, no committed ledger, no remote machine."""

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
SH = REPO / "experiments" / "t3_decomp" / "t3_cube_build.sh"
SH_REL = "experiments/t3_decomp/t3_cube_build.sh"
MAN_REL = "experiments/t3_decomp/results/t3_cube/MANIFEST.json"
READ_REL = "experiments/t3_decomp/results/k31_cube_read"
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
# uv run --frozen --no-sync python -I -B ARGS or uv run --frozen --no-sync vec check ...: the tests exit
# $HOME/tests_rc; a -c one-liner runs in python -I -B; `k31_cube_read.py row` prints the sha256 row K31 records;
# `build_sign_cube.py vertices` prints $HOME/vertices and exits $HOME/vertices_rc; `build_sign_cube.py build` exits
# $HOME/build_rc or writes a file that is not C2's bytes and a stamp ($HOME/stamp, default HEAD clean); `vec check`
# exits $HOME/check_rc; every call is logged in $HOME/calls
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
if a[0].endswith("build_sign_cube.py") and a[1] == "vertices":
    if code("vertices_rc"):
        sys.exit(code("vertices_rc"))
    print((home / "vertices").read_text().strip() if (home / "vertices").exists() else "ppn npp")
    sys.exit(0)
if a[0].endswith("build_sign_cube.py") and a[1] == "build":
    if code("build_rc"):
        sys.exit(code("build_rc"))
    out = Path(a[a.index("--out") + 1])
    out.write_bytes(b"not C2's bytes")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    st = (home / "stamp").read_text().split() if (home / "stamp").exists() else [head, "true"]
    Path(str(out) + ".stamp.json").write_text(json.dumps({"commit": st[0], "clean": st[1] == "true"}))
    sys.exit(0)
sys.exit(96)
"""


def _git_env() -> dict:
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                        "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, env=_git_env(),
                       text=True)
    return r.stdout.strip()


def _repo(tmp_path: Path, ledger: str = LEDGER, readings: bool = True) -> SimpleNamespace:
    repo, home, bindir = tmp_path / "repo", tmp_path / "home", tmp_path / "bin"
    (home / RUN_DIR).mkdir(parents=True)
    files = {SH_REL: SH.read_text(), "experiments/harness/TRIES.md": ledger, ".gitignore": "runs/\n", MAN_REL: MAN_TEXT}
    if readings:
        files |= {f"{READ_REL}/reading.json": "{}\n", f"{READ_REL}/g0.json": "{}\n"}
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
    p = L.home / "v" / "scratch" / "t3cube" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else "absent"


def _calls(L) -> list[str]:
    p = L.home / "calls"
    return p.read_text().splitlines() if p.exists() else []


def _builds(L) -> list[str]:
    return [c for c in _calls(L) if "build_sign_cube.py build" in c]


def test_the_launcher_text():
    sh = SH.read_text()
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    assert not re.search(r"\$\{[^}]*:\?[^}]*'[^}]*\}", sh)  # no apostrophe inside ${var:?message}
    assert not re.search(r"(^|[\s;|&(])rm(\s|$)", sh, flags=re.M)  # nothing deleted
    assert "set -o pipefail" in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    c2 = (REPO / "experiments/builds/results/t3two_fb01487.SHA256SUMS").read_text()
    assert f"{re.search(r'^C2_SHA=([0-9a-f]{64})', sh, flags=re.M).group(1)}  c2-s2.h5ad" in c2
    assert "--captain-ruling R72" in sh and '--cube-tables "$MANIFEST"' in sh and "build gate ppp" in sh
    order = [sh.index(f"STEP={s}\n") for s in ("tests", "vertices", "ppp", "pair", "done")]
    assert order == sorted(order) and sh.index("STEP=$v$side") > sh.index("STEP=ppp")
    for t in ("tests/test_t3_cube_build.py", "tests/test_t3_cube_build_launcher.py", "tests/test_sign_cube.py"):
        assert t in sh and (REPO / t).exists()
    assert "for side in a b" in sh and "cmp -s" in sh  # each vertex twice, one sha256


def test_usage_and_the_run_name(tmp_path):
    L = _repo(tmp_path)
    for args in ((), (L.head,), (L.head, "r", "x"), ("", "r"), (L.head, "")):
        assert _launch(L, *args)[0] == 89, args
    for bad in (".hidden", "a/b", "x y", "a.b"):
        assert _launch(L, L.head, bad)[0] == 90
    assert _calls(L) == []


def test_head_the_tree_and_the_folders(tmp_path):
    L = _repo(tmp_path)
    code, err = _launch(L, "0" * 40, "r1")
    assert code == 96 and "HEAD is not" in err and _exit(L, "r1", "gate") == "96" and _exit(L, "r1", "all") == "96"
    code, err = _launch(L, "0" * 40, "r1")
    assert code == 91 and "never reused" in err
    (L.repo / "stray.txt").write_text("x")
    code, err = _launch(L, L.head, "r2")
    assert code == 98 and "not clean" in err
    (L.repo / "stray.txt").rename(L.home / "stray.txt")
    (L.repo / "runs" / f"t3_cube_{L.head[:7]}_r3").mkdir(parents=True)
    code, err = _launch(L, L.head, "r3")
    assert code == 92 and "output folder" in err
    assert _builds(L) == []


@pytest.mark.parametrize("ledger,code,why", [
    (K31, 95, "exactly one row R72"),
    (R72.replace("R72 (the captain's order", "R72 (a ruling") + K31, 95, "does not open"),
    (R72.replace("ruling, no run", "ruled") + K31, 95, "ruling, no run"),
    (R72.replace("`--captain-ruling R72`", "--captain-ruling R72") + K31, 95, "--captain-ruling R72"),
    (R72.replace("`--form c2-s2 --vertex npn`", "`--form c2-s2 --vertex npn --x`") + K31, 95, "vertex npn"),
    (R72, 94, "exactly one row K31"),
    (R72 + K31.replace(MAN_SHA, "c" * 64), 94, "manifest sha256"),
])
def test_the_gate_refuses_rows_r72_and_k31_not_as_required(tmp_path, ledger, code, why):
    L = _repo(tmp_path, ledger)
    got, err = _launch(L, L.head, "r")
    assert got == code and why in err and _exit(L, "r", "gate") == str(code), err
    assert _builds(L) == []


def test_the_readings_must_be_committed(tmp_path):
    L = _repo(tmp_path, readings=False)
    code, err = _launch(L, L.head, "r")
    assert code == 93 and "not committed" in err and _builds(L) == []


def test_failing_tests_and_an_empty_reading_build_nothing(tmp_path):
    L = _repo(tmp_path)
    (L.home / "tests_rc").write_text("1")
    assert _launch(L, L.head, "t1")[0] == 10 and _exit(L, "t1", "vertices") == "absent"
    (L.home / "tests_rc").write_text("0")
    (L.home / "vertices_rc").write_text("3")
    assert _launch(L, L.head, "t2")[0] == 30 and _exit(L, "t2", "vertices") == "30"
    (L.home / "vertices_rc").write_text("2")
    assert _launch(L, L.head, "t3")[0] == 31
    (L.home / "vertices_rc").write_text("0")
    (L.home / "vertices").write_text("ppn ppp")
    assert _launch(L, L.head, "t4")[0] == 31  # ppp is never a file of the order
    assert _builds(L) == []


@pytest.mark.parametrize("setup,code", [
    ({"build_rc": "2"}, 12), ({"build_rc": "6"}, 16), ({"build_rc": "5"}, 17), ({"build_rc": "1"}, 11),
    ({"check_rc": "1"}, 13), ({"stamp": "f" * 40 + " true"}, 15), ({"stamp": "HEAD false"}, 15), ({}, 20),
])
def test_the_identity_gate_stops_on_anything_but_c2s_bytes(tmp_path, setup, code):
    L = _repo(tmp_path)
    for k, v in setup.items():
        (L.home / k).write_text(v.replace("HEAD", L.head))
    got, err = _launch(L, L.head, "i")
    assert got == code and _exit(L, "i", "ppp") == str(code) and _exit(L, "i", "ppna") == "absent", err
    assert len(_builds(L)) == 1 and "--vertex ppp" in _builds(L)[0] and "--captain-ruling R72" in _builds(L)[0]
