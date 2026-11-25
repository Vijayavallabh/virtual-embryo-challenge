"""The launchers of heart file OIC125's line: experiments/heart_ext/run_oic_g0.sh (G0), experiments/harness/run_h38.sh
(the counted read) and experiments/builds/heart_oic_build.sh (the build). Every call runs in a throwaway git repository
under the test's own folder, with a temporary home and a stand-in uv first on PATH: no script step runs, nothing reads a
competition or outside file, nothing is built or scored. Each refusal path is executed and its exit code checked; the
success paths make their calls in order and stop at the first step that fails; the build's identity step bites before
any OIC125 build."""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
G0_SH = ROOT / "experiments" / "heart_ext" / "run_oic_g0.sh"
H38_SH = ROOT / "experiments" / "harness" / "run_h38.sh"
BUILD_SH = ROOT / "experiments" / "builds" / "heart_oic_build.sh"
sys.path.insert(0, str(ROOT / "experiments" / "harness"))
import h38_row  # noqa: E402

SHA = {"E10.5": "a" * 64, "E9.5": "b" * 64, "E8.75": "c" * 64}
R70 = ("R70 (the captain's order of 2026-11-25, ruling, no run: \"a synthetic stand-in for OIC125\"; built only with "
       "`--captain-ruling R70`; not a try)")

STUB_UV = r'''#!/usr/bin/env python3
import json, os, shutil, subprocess, sys
from pathlib import Path
home = Path(os.environ["HOME"])
args = sys.argv[1:]
with open(home / "uvcalls", "a") as fh:
    fh.write(json.dumps(args) + "\n")
if not args or args[0] != "run":
    sys.exit(90)
a = args[1:]
while a and a[0] in ("--frozen", "--no-sync"):
    a = a[1:]
if a[0] == "vec":
    sys.exit(int(os.environ.get("STUB_CHECK_EXIT", "0")))
if a[0] != "python":
    sys.exit(91)
a = a[1:]
while a and a[0] in ("-I", "-B"):
    a = a[1:]
if a[0] == "-m":
    sys.exit(int(os.environ.get("STUB_TESTS_EXIT", "0")))
if a[0] == "-c":
    sys.exit(subprocess.run([sys.executable, "-I", "-B", *a]).returncode)
script, rest = Path(a[0]).name, a[1:]
def opt(n):
    return rest[rest.index(n) + 1] if n in rest else None
if script == "h38_heart_insitu.py":
    step = rest[0]
    if step == "refusals":
        sys.exit(1 if os.environ.get("STUB_REFUSE") == opt("--step") else 0)
    out = Path(opt("--out"))
    out.mkdir(parents=True, exist_ok=True)
    if step == "gate":
        (out / "gate.json").write_text("{}")
        sys.exit(int(os.environ.get("STUB_GATE_EXIT", "0")))
    if step == "run":
        if os.environ.get("STUB_RUN_BAD") == opt("--seeds"):
            sys.exit(7)
        (out / "h38.json").write_text("{}")
        sys.exit(0)
    if step == "read":
        (out / "reading.json").write_text("{}")
        (out / "reading.txt").write_text("x")
        sys.exit(0)
if script == "oic_g0.py":
    step = rest[0]
    code = int(os.environ.get("STUB_" + step.upper() + "_EXIT", "0"))
    if step in ("ranks", "units", "regions"):
        Path(opt("--out"), step + ".json").write_text("{}")
    sys.exit(code)
if script == "build_test.py":
    out = Path(opt("--out"))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("file " + str(opt("--trend-form")) + " " + str(opt("--insitu-chords")))
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    Path(str(out) + ".stamp.json").write_text(json.dumps({"commit": head, "clean": True}))
    sys.exit(int(os.environ.get("STUB_BUILD_EXIT", "0")))
if script == "halve.py":
    Path(rest[1]).write_text(Path(rest[0]).read_text() + " half")
    sys.exit(0)
sys.exit(92)
'''


def h38_line(cell="pending", first=h38_row.FIRST_CELL):
    return f"| {first} | heart | " + " | ".join(["x"] * 8) + f" | {cell} |\n"


G0_CELL = h38_row.g0_statement(SHA["E10.5"], SHA["E9.5"], SHA["E8.75"])


def _env(tmp, extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GIT_", "OIC_", "H38_", "STUB_"))}
    env.update(HOME=str(tmp / "home"), GIT_CEILING_DIRECTORIES=str(tmp),
               PATH=f"{tmp / 'bin'}{os.pathsep}{os.environ['PATH']}")
    env.update(extra or {})
    return env


@pytest.fixture()
def repo(tmp_path):
    root = tmp_path / "repo"
    for d in ("experiments/harness", "experiments/heart_ext/results", "experiments/builds", "docs", "data/raw/T2",
              ".venv/bin"):
        (root / d).mkdir(parents=True, exist_ok=True)
    (tmp_path / "home").mkdir()
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "uv").write_text(STUB_UV)
    (tmp_path / "bin" / "uv").chmod(0o755)
    (root / ".venv/bin/python").write_text('#!/bin/sh\ncase "$*" in *pyarrow*) exit 1 ;; esac\nexit 0\n')
    (root / ".venv/bin/python").chmod(0o755)
    (root / "experiments/harness/h38_row.py").write_text((ROOT / "experiments/harness/h38_row.py").read_text())
    for s in ("experiments/harness/h38_heart_insitu.py", "experiments/heart_ext/oic_g0.py",
              "experiments/builds/build_test.py", "experiments/builds/halve.py"):
        (root / s).write_text("# a stand-in\n")
    (root / "docs/S2_RULES.md").write_text("rules\n")
    (root / "docs/S2_RULES.sha256").write_text(hashlib.sha256(b"rules\n").hexdigest() + "  docs/S2_RULES.md\n")
    (root / ".gitignore").write_text("runs/\n.venv/\ndata/\n")
    (root / "experiments/harness/TRIES.md").write_text(f"| {R70} | heart | x |\n" + h38_line())
    fetched = tmp_path / "fetched"
    fetched.mkdir()
    (fetched / "MANIFEST.json").write_text(json.dumps({"files": [{"name": "E15.5_E1S1_GEM_bin1.tsv.gz",
                                                                  "study": "X022"}]}))

    def git(*a):
        return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "init.defaultBranch=x", *a],
                              cwd=root, env=_env(tmp_path), capture_output=True, text=True, check=True).stdout.strip()
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "one")
    return SimpleNamespace(root=root, git=git, head=git("rev-parse", "HEAD"), tmp=tmp_path, fetched=fetched,
                           home=tmp_path / "home")


def commit_change(repo, path: str, text: str) -> str:
    p = repo.root / path
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    repo.git("add", "-f", path)
    repo.git("commit", "-q", "-m", "change")
    repo.head = repo.git("rev-parse", "HEAD")
    return repo.head


def launch(repo, script, *args, **extra):
    r = subprocess.run(["bash", str(script), *args], cwd=repo.root, env=_env(repo.tmp, extra), capture_output=True,
                       text=True, timeout=120)
    return r.returncode, r.stderr


def calls(repo) -> list[list[str]]:
    f = repo.home / "uvcalls"
    return [json.loads(x) for x in f.read_text().splitlines()] if f.exists() else []


def script_steps(repo, name) -> list[str]:
    out = []
    for c in calls(repo):
        words = list(c)
        while words and words[0] in ("run", "--frozen", "--no-sync", "python", "-I", "-B"):
            words = words[1:]
        if words and words[0].endswith(name):
            out.append(words[1] if words[1] != "refusals" else "ask " + words[words.index("--step") + 1])
        elif words and words[0] == "-m":
            out.append("tests")
    return out


def test_launchers_parse_and_carry_no_apostrophe_in_a_required_parameter():
    for sh in (G0_SH, H38_SH, BUILD_SH):
        assert subprocess.run(["bash", "-n", str(sh)], capture_output=True).returncode == 0
        text = sh.read_text()
        for line in text.splitlines():
            if ":?" in line:
                assert "'" not in line.split(":?", 1)[1].split("}", 1)[0], line
        assert not re.search(r"(^|[;&|(\s])rm\s", text) and "unlink" not in text
        assert "export PATH=$HOME/.local/bin:$PATH" in text


# ---- run_oic_g0.sh ------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("args", [(), ("x",), ("a" * 40, "r1", "prep"), ("a" * 40, "r1", "other", "f"),
                                  ("a" * 40, "r 1", "prep", "f"), ("abc", "r1", "prep", "f"),
                                  ("a" * 40, "r1", "g0", "f"), ("a" * 40, "r1", "prep", "f", "extra")])
def test_g0_odd_arguments_refuse(repo, args):
    code, _ = launch(repo, G0_SH, *args)
    assert code in (1, 90, 2) and not calls(repo)


def test_g0_refusals(repo):
    f = str(repo.fetched)
    assert launch(repo, G0_SH, "f" * 40, "r1", "prep", f)[0] == 96
    assert launch(repo, G0_SH, repo.head, "r1", "prep", f)[0] == 91  # the run folder exists now
    (repo.root / "experiments/harness/TRIES.md").write_text(f"| {R70} | heart | x |\n" + h38_line() + "dirty\n")
    assert launch(repo, G0_SH, repo.head, "r2", "prep", f)[0] == 98
    repo.git("checkout", "--", "experiments/harness/TRIES.md")
    commit_change(repo, "experiments/harness/TRIES.md", f"| {R70} | heart | x |\n" + h38_line(G0_CELL))
    assert launch(repo, G0_SH, repo.head, "r3", "prep", f)[0] == 93
    commit_change(repo, "experiments/harness/TRIES.md", f"| {R70} | heart | x |\n" + h38_line())
    assert launch(repo, G0_SH, repo.head, "r4", "prep", str(repo.tmp / "nowhere"))[0] == 92
    (repo.fetched / "MANIFEST.json").write_text(json.dumps({"files": [{"study": "X021",
                                                                        "name": "x_transcripts.parquet"}]}))
    assert launch(repo, G0_SH, repo.head, "r5", "prep", f)[0] == 87  # the stand-in python has no pyarrow
    assert not calls(repo)


def test_g0_prep_runs_its_steps_in_order(repo):
    code, err = launch(repo, G0_SH, repo.head, "p1", "prep", str(repo.fetched))
    assert code == 0, err
    assert script_steps(repo, "oic_g0.py") == ["tests", "ranks", "units", "regions"]
    s = repo.home / "v/scratch/oic_g0/p1"
    assert {p.name: p.read_text().strip() for p in s.glob("*.exit")} == {
        "check.exit": "0", "tests.exit": "0", "ranks.exit": "0", "units.exit": "0", "regions.exit": "0",
        "done.exit": "0"}
    assert not (repo.home / "v/scratch/oic_g0/g0.exit").exists()


def test_g0_prep_stops_at_the_first_failing_step(repo):
    code, _ = launch(repo, G0_SH, repo.head, "p2", "prep", str(repo.fetched), STUB_UNITS_EXIT="4")
    assert code == 4 and script_steps(repo, "oic_g0.py") == ["tests", "ranks", "units"]
    assert (repo.home / "v/scratch/oic_g0/p2/done.exit").read_text().strip() == "4"


def test_g0_phase_refusals_and_once_only(repo):
    f = str(repo.fetched)
    assert launch(repo, G0_SH, repo.head, "g1", "g0", f, "nothere")[0] == 84
    assert launch(repo, G0_SH, repo.head, "p1", "prep", f)[0] == 0
    assert launch(repo, G0_SH, repo.head, "g2", "g0", f, "p1")[0] == 85  # the reader's file is not committed
    commit_change(repo, "experiments/heart_ext/oic_reader.json", "{}")
    (repo.root / "experiments/heart_ext/results/g0.json").write_text("{}")
    assert launch(repo, G0_SH, repo.head, "g3", "g0", f, "p1")[0] == 86
    commit_change(repo, "experiments/heart_ext/results/g0.json", "{}")
    repo.git("rm", "-q", "experiments/heart_ext/results/g0.json")
    repo.git("commit", "-q", "-m", "away")
    repo.head = repo.git("rev-parse", "HEAD")
    n = len(calls(repo))
    code, err = launch(repo, G0_SH, repo.head, "g4", "g0", f, "p1", STUB_G0_EXIT="3")
    assert code == 3, err
    assert script_steps(repo, "oic_g0.py")[-1] == "g0" and len(calls(repo)) == n + 2
    assert (repo.home / "v/scratch/oic_g0/g0.exit").read_text().strip() == "3"
    assert launch(repo, G0_SH, repo.head, "g5", "g0", f, "p1")[0] == 88  # G0 is started once, ever


# ---- run_h38.sh ---------------------------------------------------------------------------------------------------
@pytest.fixture()
def read_repo(repo):
    commit_change(repo, "experiments/harness/TRIES.md", f"| {R70} | heart | x |\n" + h38_line(G0_CELL))
    return repo


@pytest.mark.parametrize("args", [(), ("x",), ("abc", "r"), ("a" * 40, "r/1"), ("a" * 40, "r", "extra")])
def test_h38_odd_arguments_refuse(read_repo, args):
    assert launch(read_repo, H38_SH, *args)[0] in (1, 90) and not calls(read_repo)


def test_h38_refusals(read_repo):
    r = read_repo
    assert launch(r, H38_SH, "f" * 40, "a")[0] == 96
    assert launch(r, H38_SH, r.head, "a")[0] == 91
    commit_change(r, "experiments/harness/TRIES.md", f"| {R70} | heart | x |\n" + h38_line())
    assert launch(r, H38_SH, r.head, "b")[0] == 93  # no G0 statement
    commit_change(r, "experiments/harness/TRIES.md",
                  f"| {R70} | heart | x |\n" + h38_line(G0_CELL + " " + h38_row.outcome_statement("STOP")))
    assert launch(r, H38_SH, r.head, "c")[0] == 93  # a read outcome stands already
    commit_change(r, "experiments/harness/TRIES.md", f"| {R70} | heart | x |\n" + h38_line(G0_CELL))
    out = r.root / "experiments/harness/results/h38_heart_insitu"
    out.mkdir(parents=True)
    (out / "reading.json").write_text("{}")
    assert launch(r, H38_SH, r.head, "d")[0] == 95
    (out / "reading.json").rename(out / "kept.json")
    assert launch(r, H38_SH, r.head, "e", H38_VENV=str(r.tmp / "novenv"))[0] == 99
    assert launch(r, H38_SH, r.head, "f", H38_DATA=str(r.tmp / "nodata"))[0] == 97
    assert not calls(r)
    assert launch(r, H38_SH, r.head, "g", STUB_REFUSE="gate")[0] == 89
    assert script_steps(r, "h38_heart_insitu.py") == ["tests", "ask gate"]
    assert not (r.home / "v/scratch/h38/gate.exit").exists()


def test_h38_runs_its_steps_once_each_in_order(read_repo):
    r = read_repo
    code, err = launch(r, H38_SH, r.head, "go")
    assert code == 0, err
    assert script_steps(r, "h38_heart_insitu.py") == ["tests", "ask gate", "gate", "ask run", *["run"] * 10,
                                                       "ask read", "read"]
    base = r.home / "v/scratch/h38"
    assert {p.name: p.read_text().strip() for p in base.glob("*.exit")} == {"gate.exit": "0", "run.exit": "0",
                                                                            "read.exit": "0"}
    assert launch(r, H38_SH, r.head, "again")[0] in (88, 95)


def test_h38_stops_at_a_failing_gate_and_never_restarts_it(read_repo):
    r = read_repo
    assert launch(r, H38_SH, r.head, "x1", STUB_GATE_EXIT="3")[0] == 3
    assert script_steps(r, "h38_heart_insitu.py")[-1] == "gate"
    assert (r.home / "v/scratch/h38/gate.exit").read_text().strip() == "3"
    (r.root / "experiments/harness/results/h38_heart_insitu/gate/gate.json").rename(r.tmp / "gate.json")
    assert launch(r, H38_SH, r.head, "x2")[0] == 88


def test_h38_a_failing_part_fails_the_run(read_repo):
    r = read_repo
    assert launch(r, H38_SH, r.head, "y1", STUB_RUN_BAD="4-5")[0] == 1
    assert "read" not in script_steps(r, "h38_heart_insitu.py")


# ---- heart_oic_build.sh -------------------------------------------------------------------------------------------
@pytest.fixture()
def build_repo(repo):
    chords = '{"row": "H38"}\n'
    sha = hashlib.sha256(chords.encode()).hexdigest()
    cell = h38_row.g0_statement(sha, SHA["E9.5"], SHA["E8.75"]) + " " + h38_row.outcome_statement("PASS")
    (repo.root / "experiments/heart_ext/results/chords_E10.5.json").write_text(chords)
    repo.git("add", "-f", "experiments/heart_ext/results/chords_E10.5.json")
    commit_change(repo, "experiments/harness/TRIES.md", f"| {R70} | heart | x |\n" + h38_line(cell))
    repo.sha = sha
    return repo


@pytest.mark.parametrize("args", [(), ("x",), ("abc", "r"), ("a" * 40, "r 1"), ("a" * 40, "r", "extra")])
def test_build_odd_arguments_refuse(build_repo, args):
    assert launch(build_repo, BUILD_SH, *args)[0] in (1, 90) and not calls(build_repo)


def test_build_gate_refusals(build_repo):
    r = build_repo
    assert launch(r, BUILD_SH, "f" * 40, "a")[0] == 96
    assert launch(r, BUILD_SH, r.head, "a")[0] == 91
    (r.root / "docs/S2_RULES.md").write_text("changed\n")
    assert launch(r, BUILD_SH, r.head, "b")[0] == 98
    commit_change(r, "docs/S2_RULES.md", "changed\n")
    assert launch(r, BUILD_SH, r.head, "c")[0] == 97
    commit_change(r, "docs/S2_RULES.md", "rules\n")
    text = (r.root / "experiments/harness/TRIES.md").read_text()
    commit_change(r, "experiments/harness/TRIES.md", text.replace("| R70 (", "| R71 ("))
    assert launch(r, BUILD_SH, r.head, "d")[0] == 95
    commit_change(r, "experiments/harness/TRIES.md", text.replace("OUTCOME H38: PASS.", "OUTCOME H38: STOP."))
    assert launch(r, BUILD_SH, r.head, "e")[0] == 94
    commit_change(r, "experiments/harness/TRIES.md", text)
    commit_change(r, "experiments/heart_ext/results/chords_E10.5.json", '{"row": "other"}\n')
    assert launch(r, BUILD_SH, r.head, "f")[0] == 92
    commit_change(r, "experiments/heart_ext/results/chords_E10.5.json", '{"row": "H38"}\n')
    (r.root / "runs").mkdir(exist_ok=True)
    (r.root / f"runs/heartoic_{r.head[:7]}_g").mkdir()
    assert launch(r, BUILD_SH, r.head, "g")[0] == 93
    assert not calls(r)


def test_build_identity_step_bites_before_any_oic_build(build_repo):
    r = build_repo
    code, err = launch(r, BUILD_SH, r.head, "id")
    assert code == 20, err
    builds = [c for c in calls(r) if any(w.endswith("build_test.py") for w in c)]
    assert len(builds) == 1 and "category" in builds[0] and "insitu" not in builds[0]


def test_build_builds_oic125_twice_after_the_identity(build_repo, tmp_path):
    r = build_repo
    otk = hashlib.sha256(b"file category None").hexdigest()
    half = hashlib.sha256(b"file category None half").hexdigest()
    text = BUILD_SH.read_text()
    text = text.replace("d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893", otk)
    text = text.replace("9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54", half)
    patched = tmp_path / "heart_oic_build_patched.sh"
    patched.write_text(text)
    code, err = launch(r, patched, r.head, "ok")
    assert code == 0, err
    builds = [c for c in calls(r) if any(w.endswith("build_test.py") for w in c)]
    assert len(builds) == 3
    for c in builds[1:]:
        i = c.index("--insitu-chords")
        assert c[i + 1] == r.sha and c[c.index("--captain-ruling") + 1] == "R70" and "insitu" in c
    s = r.home / "v/scratch/heartoic/ok"
    assert (s / "done.exit").read_text().strip() == "0"
    assert "OIC125a/T2heart_organ_ts125_ti.h5ad" in (s / "SHA256SUMS").read_text()
    code, _ = launch(r, patched, r.head, "bad", STUB_CHECK_EXIT="1")
    assert code == 13
