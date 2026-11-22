"""The launchers of the heart's organ-cut file D: experiments/harness/run_h29b.sh (row H29b's practice read) and
experiments/builds/heart_organ.sh (its test-board builds). Every call runs in a throwaway git repository under the
test's own folder, with a temporary home and a stand-in uv (and, where named, a stand-in sha256sum) first on PATH:
no launcher runs with a real commit, nothing reads a competition file and nothing is built."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
READ_LAUNCHER = REPO_ROOT / "experiments" / "harness" / "run_h29b.sh"
BUILD_LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_organ.sh"
HEART3 = REPO_ROOT / "experiments" / "builds" / "results" / "heart3_3d52b83.SHA256SUMS"
C_SHA = "341058344dbd73cd4267f6512bce0f185431ab16b6691d87b28c379f8094c299"
C_HALF_SHA = "04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35"

# A stand-in for uv: it records each call and imitates each kind of call the two launchers make.
STUB_UV = r"""#!/usr/bin/env bash
echo "$*" >> "$HOME/uvcalls"
args=("$@")
[ "${args[0]}" = run ] || exit 90
args=("${args[@]:1}")
while [ "${args[0]}" = --frozen ] || [ "${args[0]}" = --no-sync ]; do args=("${args[@]:1}"); done
if [ "${args[0]}" = vec ]; then exit "${STUB_CHECK_EXIT:-0}"; fi
[ "${args[0]}" = python ] || exit 91
args=("${args[@]:1}")
while [ "${args[0]}" = -I ] || [ "${args[0]}" = -B ]; do args=("${args[@]:1}"); done
out=""; prev=""
for a in "${args[@]}"; do [ "$prev" = "--out" ] && out="$a"; prev="$a"; done
case "${args[0]}" in
  -m) exit "${STUB_PYTEST_EXIT:-0}" ;;
  -c) exec python3 -I -B "${args[@]}" ;;
  experiments/harness/h29b_heart_organ.py)
    case "${args[1]}" in
      gate) mkdir -p "$out"; echo "GATE PASS"; exit "${STUB_GATE_EXIT:-0}" ;;
      run) mkdir -p "$out"; echo '{"stub": 1}' > "$out/h29b.json"; echo "T2h/x seed 0 (1.0 s) done"
           echo "a warning line"; echo "wrote $out/h29b.json"; exit 0 ;;
      read) echo "OUTCOME H29b: READ"; echo "H29b READ: follows"; exit 0 ;;
    esac ;;
  experiments/builds/build_test.py)
    [ "${STUB_BUILD_EXIT:-0}" = 0 ] || exit "$STUB_BUILD_EXIT"
    noise=""; [ -n "${STUB_NONDET:-}" ] && noise="$RANDOM$RANDOM$(date +%N)"
    spec=""; for a in "${args[@]:1}"; do [ "$a" = --out ] && break; spec="$spec $a"; done
    echo "built$spec $noise" > "$out"
    clean=True; [ -n "${STUB_DIRTY:-}" ] && clean=False
    python3 -I -B -c 'import json,sys; json.dump({"commit": sys.argv[2], "clean": sys.argv[3] == "True"}, open(sys.argv[1], "w"))' "$out.stamp.json" "$(git rev-parse HEAD)" "$clean"
    exit 0 ;;
  experiments/builds/halve.py)
    { echo "half of"; cat "${args[1]}"; } > "${args[2]}"; exit 0 ;;
esac
exit 92
"""

# A stand-in sha256sum that gives set heart3's hashes for the full trim's two files and the real hash otherwise.
STUB_SHA = r"""#!/usr/bin/env bash
REAL=__REAL__
if [ "${1:-}" = -c ]; then exec "$REAL" "$@"; fi
for f in "$@"; do
  case "$(basename "$f")" in
    T2heart_s2core.h5ad) echo "__C__  $f" ;;
    T2heart_s2core_half.h5ad) echo "__CH__  $f" ;;
    *) "$REAL" "$f" ;;
  esac
done
"""


def _env(tmp_path, bin_dir=None, extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path))
    if bin_dir is not None:
        env["PATH"] = f"{bin_dir}{os.pathsep}{os.environ['PATH']}"
    env.update(extra or {})
    return env


def _git_fn(repo, tmp_path):
    def git(*a):
        return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "init.defaultBranch=x", *a],
                              cwd=repo, env=_env(tmp_path), capture_output=True, text=True, check=True).stdout.strip()
    return git


def _bin(tmp_path, sha_stub=False):
    b = tmp_path / ("bin_sha" if sha_stub else "bin")
    b.mkdir(exist_ok=True)
    (b / "uv").write_text(STUB_UV)
    (b / "uv").chmod(0o755)
    if sha_stub:
        real = shutil.which("sha256sum")
        (b / "sha256sum").write_text(STUB_SHA.replace("__REAL__", real).replace("__CH__", C_HALF_SHA)
                                     .replace("__C__", C_SHA))
        (b / "sha256sum").chmod(0o755)
    return b


def _calls(tmp_path):
    f = tmp_path / "home" / "uvcalls"
    return f.read_text().splitlines() if f.exists() else []


# ================================================================================================ run_h29b.sh
@pytest.fixture()
def read_repo(tmp_path):
    repo = tmp_path / "repo"
    (repo / "experiments" / "harness").mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| H29a | x |\n| H29b (the row) | heart | x |\n")
    (repo / ".gitignore").write_text("runs/\n")
    git = _git_fn(repo, tmp_path)
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "row")
    return repo, git, git("rev-parse", "HEAD"), tmp_path / "home" / "v" / "scratch" / "h29b"


def _rlaunch(tmp_path, *args, extra=None, stub=True):
    r = subprocess.run(["bash", str(READ_LAUNCHER), *args], cwd=tmp_path / "repo", capture_output=True, text=True,
                       env=_env(tmp_path, _bin(tmp_path) if stub else None, extra), timeout=120)
    return r.returncode, r.stdout, r.stderr


def test_the_read_launcher_refusals_never_touch_a_finished_exit(tmp_path, read_repo):
    repo, git, head, S = read_repo

    def f(name):
        p = S / name
        return p.read_text().strip() if p.exists() else None
    assert _rlaunch(tmp_path)[0] == 2 and f("usage.refused.exit") == "2"
    assert _rlaunch(tmp_path, "tests")[0] == 2
    assert _rlaunch(tmp_path, "build", head)[0] == 97 and f("unknown.refused.exit") == "97"
    assert _rlaunch(tmp_path, "a/b", head)[0] == 97 and not (S / "a").exists()
    for step in ("gate", "run", "read"):  # before the step before it: 96, its own exit file never written
        assert _rlaunch(tmp_path, step, head)[0] == 96 and f(f"{step}.refused.exit") == "96" and f(f"{step}.exit") is None
    code, out, _ = _rlaunch(tmp_path, "tests", "d" * 40)
    assert code == 99 and f("tests.refused.exit") == "99" and f("tests.exit") is None and out == ""
    assert _rlaunch(tmp_path, "tests", head[:12])[0] == 99
    (repo / "g.txt").write_text("y\n")
    assert _rlaunch(tmp_path, "tests", head)[0] == 98 and f("tests.refused.exit") == "98"
    (repo / "g.txt").rename(tmp_path / "g.txt")
    # the tests step passes: running, then its head, then its exit
    assert _rlaunch(tmp_path, "tests", head)[0] == 0 and f("tests.exit") == "0" and f("tests.head") == head
    assert any("-m pytest" in c and "tests/test_h29b_heart_organ.py" in c for c in _calls(tmp_path))
    # once a step has exited 0 at this commit it is refused (94) and its exit stays 0
    assert _rlaunch(tmp_path, "tests", head)[0] == 94 and f("tests.refused.exit") == "94" and f("tests.exit") == "0"
    # a failing gate (STOP 3): the run is refused
    assert _rlaunch(tmp_path, "gate", head, extra={"STUB_GATE_EXIT": "3"})[0] == 3 and f("gate.exit") == "3"
    assert _rlaunch(tmp_path, "run", head)[0] == 96 and f("run.exit") is None
    # a step killed part way leaves "running": the next step refuses
    (S / "gate.exit").write_text("running\n")
    assert _rlaunch(tmp_path, "run", head)[0] == 96


def test_the_read_launcher_in_order_with_a_stand_in_uv(tmp_path, read_repo):
    repo, git, head, S = read_repo
    for step in ("tests", "gate", "run"):
        assert _rlaunch(tmp_path, step, head)[0] == 0, step
        assert (S / f"{step}.exit").read_text().strip() == "0" and (S / f"{step}.head").read_text().strip() == head
    R = repo / "runs" / "h29b" / head
    assert len(list(R.rglob("h29b.json"))) == 15
    runs = [c for c in _calls(tmp_path) if " run --task " in c]
    assert len(runs) == 15 and all(f"--gate runs/h29b/{head}/gate --commit {head}" in c for c in runs)
    assert {c.split("--task ")[1].split(" ")[0] for c in runs} == {"T2h/<=E9.5->E10.5", "T2h/E8.25_late->E8.5",
                                                                   "T2h/<=E8.75->E10.5"}
    # the stand-in prints one other line per part: counted, never printed
    assert (R / "exit" / "other_lines").read_text().strip() == "15"
    gate = [c for c in _calls(tmp_path) if " gate --out " in c]
    assert gate == [f"run --frozen --no-sync python -I -B experiments/harness/h29b_heart_organ.py gate --out "
                    f"runs/h29b/{head}/gate --commit {head}"]
    code, out, _ = _rlaunch(tmp_path, "read", head)
    assert code == 0 and out == "" and (S / "read.out").read_text() == "OUTCOME H29b: READ\nH29b READ: follows\n"
    read = [c for c in _calls(tmp_path) if " read --runs " in c]
    assert read == [f"run --frozen --no-sync python -I -B experiments/harness/h29b_heart_organ.py read --runs "
                    f"runs/h29b/{head} --gate runs/h29b/{head}/gate --commit {head} --out "
                    "experiments/harness/results/h29b_heart_organ"]


def test_the_read_launcher_refuses_without_the_row_and_over_a_reading(tmp_path, read_repo):
    repo, git, head, S = read_repo
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| H29a | x |\n")
    git("commit", "-q", "-am", "no row")
    h2 = git("rev-parse", "HEAD")
    assert _rlaunch(tmp_path, "tests", h2)[0] == 93 and (S / "tests.refused.exit").read_text().strip() == "93"
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| H29b x | heart |\n")
    out = repo / "experiments" / "harness" / "results" / "h29b_heart_organ"
    out.mkdir(parents=True)
    (out / "reading.txt").write_text("the one reading\n")
    git("add", "-A")
    git("commit", "-q", "-m", "a reading")
    h3 = git("rev-parse", "HEAD")
    for step in ("tests", "gate", "run"):
        (S / f"{step}.exit").write_text("0\n")
        (S / f"{step}.head").write_text(h3 + "\n")
    assert _rlaunch(tmp_path, "read", h3)[0] == 95 and (out / "reading.txt").read_text() == "the one reading\n"
    assert not (S / "read.exit").exists()


# ================================================================================================ heart_organ.sh
@pytest.fixture()
def build_repo(tmp_path):
    """A throwaway repository: the rules file and its pin, the ledger with rows R21 and H29b, a stand-in src/, then
    (the run commit's descendant) row H29b's committed reading made at the run commit."""
    repo = tmp_path / "repo"
    for d in ("docs", "experiments/harness/results/h29b_heart_organ", "src"):
        (repo / d).mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (repo / "docs" / "S2_RULES.md").write_text("rules\n")
    sha = subprocess.run(["sha256sum", "docs/S2_RULES.md"], cwd=repo, capture_output=True, text=True).stdout
    (repo / "docs" / "S2_RULES.sha256").write_text(sha)
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| R21 (ruling) | heart |\n| H29b (row) | heart |\n")
    (repo / "src" / "a.py").write_text("x = 1\n")
    (repo / ".gitignore").write_text("runs/\n")
    git = _git_fn(repo, tmp_path)
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "row")
    runc = git("rev-parse", "HEAD")
    reading = repo / "experiments" / "harness" / "results" / "h29b_heart_organ" / "reading.json"

    def commit_reading(**over):
        r = {"void": False, "outcome": "READ", "commit": runc, "gate": {"verdict": "PASS", "commit": runc}}
        r.update(over)
        reading.write_text(json.dumps(r))
        git("add", "-A")
        git("commit", "-q", "-m", "reading")
        return git("rev-parse", "HEAD")
    head = commit_reading()
    return repo, git, runc, head, commit_reading


def _blaunch(tmp_path, *args, sha_stub=True, extra=None):
    r = subprocess.run(["bash", str(BUILD_LAUNCHER), *args], cwd=tmp_path / "repo", capture_output=True, text=True,
                       env=_env(tmp_path, _bin(tmp_path, sha_stub), extra), timeout=120)
    return r.returncode, r.stdout, r.stderr


def _bexit(tmp_path, run, step):
    p = tmp_path / "home" / "v" / "scratch" / "heartorgan" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


def test_the_build_launcher_end_to_end_with_stand_ins(tmp_path, build_repo):
    repo, git, runc, head, _ = build_repo
    code, out, err = _blaunch(tmp_path, head, "r1")
    assert code == 0, err
    assert out == ""
    for step in ("gate", "tests", "C", "Da", "Db", "done"):
        assert _bexit(tmp_path, "r1", step) == "0", step
    S = tmp_path / "home" / "v" / "scratch" / "heartorgan" / "r1"
    O = repo / "runs" / f"heartorgan_{head[:7]}_r1"
    assert (S / "head").read_text().strip() == head
    assert (S / "stamps").read_text().splitlines() == [f"C/T2heart_s2core {head} True",
                                                       f"Da/T2heart_organ {head} True",
                                                       f"Db/T2heart_organ {head} True"]
    sums = (S / "SHA256SUMS").read_text()
    assert sums == (O / "SHA256SUMS").read_text()
    assert [ln.split()[1] for ln in sums.splitlines()] == [
        "C/T2heart_s2core.h5ad", "C/T2heart_s2core_half.h5ad", "Da/T2heart_organ.h5ad", "Da/T2heart_organ_half.h5ad",
        "Db/T2heart_organ.h5ad", "Db/T2heart_organ_half.h5ad"]
    calls = _calls(tmp_path)
    builds = [c for c in calls if "build_test.py" in c]
    assert builds == [
        f"run --frozen --no-sync python -I -B experiments/builds/build_test.py heart --entry core_trim --expression "
        f"copy --z-mode raw --out {O.relative_to(repo)}/C/T2heart_s2core.h5ad",
        *(f"run --frozen --no-sync python -I -B experiments/builds/build_test.py heart --entry organ_scope "
          f"--expression copy --z-mode raw --out {O.relative_to(repo)}/{d}/T2heart_organ.h5ad" for d in ("Da", "Db"))]
    halves = [c for c in calls if "experiments/builds/halve.py" in c]
    assert len(halves) == 3 and all("--seed" not in c for c in halves)  # halve.py's default seed, 1000001
    checks = [c for c in calls if " vec check " in c]
    assert len(checks) == 6 and all(c.endswith("--board T2:heart:test_extrap") for c in checks)
    # a run folder is never reused
    assert _blaunch(tmp_path, head, "r1")[0] == 91


def test_the_build_launcher_stops_unless_the_full_trim_gives_set_heart3(tmp_path, build_repo):
    _, _, _, head, _ = build_repo
    assert _blaunch(tmp_path, head, "r2", sha_stub=False)[0] == 20
    assert _bexit(tmp_path, "r2", "C") == "20" and _bexit(tmp_path, "r2", "Da") is None
    assert not any("--entry organ_scope" in c for c in _calls(tmp_path))


def test_the_build_launcher_fails_on_a_build_a_stamp_or_two_different_files(tmp_path, build_repo):
    _, _, _, head, _ = build_repo
    assert _blaunch(tmp_path, head, "b1", extra={"STUB_BUILD_EXIT": "1"})[0] == 11 and _bexit(tmp_path, "b1", "C") == "11"
    assert _blaunch(tmp_path, head, "b2", extra={"STUB_CHECK_EXIT": "1"})[0] == 13
    assert _blaunch(tmp_path, head, "b3", extra={"STUB_DIRTY": "1"})[0] == 15
    assert _blaunch(tmp_path, head, "b4", extra={"STUB_PYTEST_EXIT": "1"})[0] == 10
    assert _bexit(tmp_path, "b4", "tests") == "10" and _bexit(tmp_path, "b4", "C") is None
    code = _blaunch(tmp_path, head, "b5", extra={"STUB_NONDET": "1"})[0]
    assert code == 31 and _bexit(tmp_path, "b5", "done") == "31"


def test_the_build_launcher_gate_refusals(tmp_path, build_repo):
    repo, git, runc, head, commit_reading = build_repo
    assert _blaunch(tmp_path, head[:12], "g0")[0] == 90
    assert _blaunch(tmp_path, head, "a/b")[0] == 90
    assert _blaunch(tmp_path, head, "it's")[0] == 90
    assert _blaunch(tmp_path, runc, "g1")[0] == 96 and _bexit(tmp_path, "g1", "gate") == "96"
    (repo / "src" / "a.py").write_text("x = 2\n")
    assert _blaunch(tmp_path, head, "g2")[0] == 98
    git("checkout", "-q", "--", "src/a.py")
    (repo / "runs").mkdir()
    (repo / "runs" / f"heartorgan_{head[:7]}_g3").mkdir()
    assert _blaunch(tmp_path, head, "g3")[0] == 93
    for i, over in enumerate(({"void": True}, {"outcome": "VOID"}, {"gate": {"verdict": "STOP", "commit": runc}},
                              {"gate": {"verdict": "PASS", "commit": "d" * 40}}, {"commit": runc[:12]},
                              {"gate": None})):
        h = commit_reading(**over)
        assert _blaunch(tmp_path, h, f"v{i}")[0] == 92, over
    # the reading's run commit must be an ancestor, and the builder's code unchanged since it
    other = "e" * 40
    h = commit_reading(commit=other, gate={"verdict": "PASS", "commit": other})
    assert _blaunch(tmp_path, h, "anc")[0] == 88
    commit_reading()
    (repo / "src" / "a.py").write_text("x = 3\n")
    git("commit", "-q", "-am", "code changed after the read")
    assert _blaunch(tmp_path, git("rev-parse", "HEAD"), "code")[0] == 89
    git("revert", "--no-edit", "HEAD")
    # the rules pin, the ledger rows, a tracked reading
    (repo / "docs" / "S2_RULES.md").write_text("rules, changed\n")
    git("commit", "-q", "-am", "rules changed without a new pin")
    assert _blaunch(tmp_path, git("rev-parse", "HEAD"), "pin")[0] == 97
    git("revert", "--no-edit", "HEAD")
    for row in ("R21", "H29b"):
        t = repo / "experiments" / "harness" / "TRIES.md"
        t.write_text("".join(ln for ln in t.read_text().splitlines(keepends=True) if not ln.startswith(f"| {row} ")))
        git("commit", "-q", "-am", f"no row {row}")
        assert _blaunch(tmp_path, git("rev-parse", "HEAD"), f"row{row}")[0] == 95
        git("revert", "--no-edit", "HEAD")
    git("rm", "-q", "--cached", "experiments/harness/results/h29b_heart_organ/reading.json")
    git("commit", "-q", "-m", "reading untracked")
    assert _blaunch(tmp_path, git("rev-parse", "HEAD"), "untracked")[0] == 94
    assert not any("build_test.py" in c for c in _calls(tmp_path))


# ================================================================================================ both, as text
@pytest.mark.parametrize("path", [READ_LAUNCHER, BUILD_LAUNCHER])
def test_the_launchers_delete_nothing_and_parse(path):
    sh = path.read_text()
    body = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert not re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|truncate|shred)\s", body)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    assert "python -I -B" in sh and 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'TMPDIR="$S/tmp"' in sh and "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh
    assert "/mnt/" not in sh and "vec-campaign/" not in body
    assert subprocess.run(["bash", "-n", str(path)], capture_output=True).returncode == 0


def test_the_build_launcher_reads_only_exit_sha_and_stamp(tmp_path):
    sh = BUILD_LAUNCHER.read_text()
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    for line in (ln for ln in sh.splitlines() if not ln.lstrip().startswith("#")):
        if line.lstrip().startswith(("$PY experiments/builds/", "uv run --frozen --no-sync vec check")):
            assert "> /dev/null 2>&1" in line, line
    body = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert 'print(s["commit"], s["clean"])' in body and '"size"' not in body
    assert not re.search(r"(^|[\s;&|(])(stat|du|wc|ls|cat|head|tail|less|h5ls)\s", body.replace('cat "$S', ""))
    # the identity step's hashes are set heart3's committed DGX list
    listed = dict(reversed(ln.split()) for ln in HEART3.read_text().splitlines())
    assert listed["builds/T2heart_s2core.h5ad"] == C_SHA and listed["builds/T2heart_s2core_half.h5ad"] == C_HALF_SHA
    assert f"C_SHA={C_SHA}" in sh and f"C_HALF_SHA={C_HALF_SHA}" in sh
    assert 'C_ARGS="--entry core_trim --expression copy --z-mode raw"' in sh
    assert 'D_ARGS="--entry organ_scope --expression copy --z-mode raw"' in sh
    assert "BOARD=T2:heart:test_extrap" in sh
