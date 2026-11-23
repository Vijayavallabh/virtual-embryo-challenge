"""experiments/builds/heart_trend_line_d.sh, the test-board builds of the heart search line's step on file D, files
OTS05 and OTS10 (row R29 of experiments/harness/TRIES.md). Every call runs in a throwaway git repository under the
test's own folder, with a temporary home and a stand-in uv (and, where named, a stand-in sha256sum) first on PATH: no
launcher runs with a real commit, nothing reads a competition file and nothing is built. The uv stand-in is
tests/test_heart_organ_launchers' own, so the launchers are held to the same calls; the sha256sum stand-in gives the
two identity steps' committed hashes."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_heart_organ_launchers import STUB_UV, _calls, _env, _git_fn  # noqa: E402

LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line_d.sh"
LINE = REPO_ROOT / "experiments" / "builds" / "heart_trend_line.sh"
ORGAN = REPO_ROOT / "experiments" / "builds" / "heart_organ.sh"
D_LIST = REPO_ROOT / "experiments" / "builds" / "results" / "heartorgan_72b9662_d1.SHA256SUMS"
TS_LIST = REPO_ROOT / "experiments" / "builds" / "results" / "heartline_d41cd14_ts1.SHA256SUMS"
D_SHA = "76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b"
D_HALF_SHA = "1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00"
TS05_SHA = "a6edffa9ab27bd8740db54db652832cbf56278f09154d7ce33d72b884db7d741"
TS05_HALF_SHA = "b6b8155be9b5469a2dbc8943a59df567f52c7cbe013bbf0b55e4da706fda27b7"
R29 = ("| R29 (the captain's order, \"i want it too\": a synthetic stand-in; ruling, no run, not a try; the heart "
       "builder is to take `--entry organ_scope` with `--trend-strength 0.5` and `--trend-strength 1.0` only with "
       "`--captain-ruling R29`) | heart | x |\n")
STEPS = ("gate", "tests", "D", "T", "OTS05a", "OTS05b", "OTS10a", "OTS10b", "done")

# A stand-in sha256sum that gives the committed hashes for file D's and TS05's two files and the real hash otherwise.
STUB_SHA = r"""#!/usr/bin/env bash
REAL=__REAL__
if [ "${1:-}" = -c ]; then exec "$REAL" "$@"; fi
for f in "$@"; do
  case "$(basename "$f")" in
    T2heart_organ.h5ad) echo "__D__  $f" ;;
    T2heart_organ_half.h5ad) echo "__DH__  $f" ;;
    T2heart_core_ts05.h5ad) echo "${STUB_TS05:-__T__}  $f" ;;
    T2heart_core_ts05_half.h5ad) echo "__TH__  $f" ;;
    *) "$REAL" "$f" ;;
  esac
done
"""


def _bin(tmp_path, sha_stub=False):
    b = tmp_path / ("bin_sha" if sha_stub else "bin")
    b.mkdir(exist_ok=True)
    (b / "uv").write_text(STUB_UV)
    (b / "uv").chmod(0o755)
    if sha_stub:
        real = shutil.which("sha256sum")
        (b / "sha256sum").write_text(STUB_SHA.replace("__REAL__", real).replace("__DH__", D_HALF_SHA)
                                     .replace("__D__", D_SHA).replace("__TH__", TS05_HALF_SHA)
                                     .replace("__T__", TS05_SHA))
        (b / "sha256sum").chmod(0o755)
    return b


@pytest.fixture()
def repo(tmp_path):
    """A throwaway repository: the rules file and its pin, the ledger with row R29 as the captain's order, a stand-in
    src/."""
    repo = tmp_path / "repo"
    for d in ("docs", "experiments/harness", "src"):
        (repo / d).mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (repo / "docs" / "S2_RULES.md").write_text("rules\n")
    sha = subprocess.run(["sha256sum", "docs/S2_RULES.md"], cwd=repo, capture_output=True, text=True).stdout
    (repo / "docs" / "S2_RULES.sha256").write_text(sha)
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| R28 (other) | heart | x |\n" + R29
                                                               + "| H28 (row) | heart | x |\n")
    (repo / "src" / "a.py").write_text("x = 1\n")
    (repo / ".gitignore").write_text("runs/\n")
    git = _git_fn(repo, tmp_path)
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "row")
    return repo, git, git("rev-parse", "HEAD")


def _launch(tmp_path, *args, sha_stub=True, extra=None):
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=tmp_path / "repo", capture_output=True, text=True,
                       env=_env(tmp_path, _bin(tmp_path, sha_stub), extra), timeout=120)
    return r.returncode, r.stdout, r.stderr


def _exit(tmp_path, run, step):
    p = tmp_path / "home" / "v" / "scratch" / "heartlined" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


def test_the_launcher_end_to_end_with_stand_ins(tmp_path, repo):
    repo, git, head = repo
    code, out, err = _launch(tmp_path, head, "r1")
    assert code == 0, err
    assert out == ""
    for step in STEPS:
        assert _exit(tmp_path, "r1", step) == "0", step
    S = tmp_path / "home" / "v" / "scratch" / "heartlined" / "r1"
    rel = Path("runs") / f"heartlined_{head[:7]}_r1"
    assert (S / "head").read_text().strip() == head
    assert (S / "stamps").read_text().splitlines() == [f"D/T2heart_organ {head} True",
                                                       f"T/T2heart_core_ts05 {head} True",
                                                       f"OTS05a/T2heart_organ_ts05 {head} True",
                                                       f"OTS05b/T2heart_organ_ts05 {head} True",
                                                       f"OTS10a/T2heart_organ_ts10 {head} True",
                                                       f"OTS10b/T2heart_organ_ts10 {head} True"]
    sums = (S / "SHA256SUMS").read_text()
    assert sums == (repo / rel / "SHA256SUMS").read_text()
    assert [ln.split()[1] for ln in sums.splitlines()] == [
        "D/T2heart_organ.h5ad", "D/T2heart_organ_half.h5ad",
        "T/T2heart_core_ts05.h5ad", "T/T2heart_core_ts05_half.h5ad",
        "OTS05a/T2heart_organ_ts05.h5ad", "OTS05a/T2heart_organ_ts05_half.h5ad",
        "OTS05b/T2heart_organ_ts05.h5ad", "OTS05b/T2heart_organ_ts05_half.h5ad",
        "OTS10a/T2heart_organ_ts10.h5ad", "OTS10a/T2heart_organ_ts10_half.h5ad",
        "OTS10b/T2heart_organ_ts10.h5ad", "OTS10b/T2heart_organ_ts10_half.h5ad"]
    calls = _calls(tmp_path)
    py = "run --frozen --no-sync python -I -B experiments/builds/build_test.py heart"
    assert [c for c in calls if "build_test.py" in c] == [
        f"{py} --entry organ_scope --expression copy --z-mode raw --out {rel}/D/T2heart_organ.h5ad",
        f"{py} --entry core_trim --trend-strength 0.5 --z-mode raw --captain-ruling R28 "
        f"--out {rel}/T/T2heart_core_ts05.h5ad",
        *(f"{py} --entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29 "
          f"--out {rel}/{d}/T2heart_organ_ts05.h5ad" for d in ("OTS05a", "OTS05b")),
        *(f"{py} --entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29 "
          f"--out {rel}/{d}/T2heart_organ_ts10.h5ad" for d in ("OTS10a", "OTS10b"))]
    halves = [c for c in calls if "experiments/builds/halve.py" in c]
    assert len(halves) == 6 and all("--seed" not in c for c in halves)  # halve.py's default seed, 1000001
    assert halves[2].endswith(f"{rel}/OTS05a/T2heart_organ_ts05.h5ad {rel}/OTS05a/T2heart_organ_ts05_half.h5ad")
    checks = [c for c in calls if " vec check " in c]
    assert len(checks) == 12 and all(c.endswith("--board T2:heart:test_extrap") for c in checks)
    tests = [c for c in calls if " -m pytest " in c]
    assert len(tests) == 1
    for f in ("test_heart_trend_line_d.py", "test_heart_trend_line_d_launcher.py", "test_heart_trend_line.py",
              "test_heart_trend_line_launcher.py", "test_build_heart_organ_scope.py", "test_heart_organ_launchers.py",
              "test_h28_heart_trend_strength.py", "test_halve.py", "test_hidden_counts.py"):
        assert f"tests/{f}" in tests[0], f
    # a run folder is never reused
    assert _launch(tmp_path, head, "r1")[0] == 91


def test_the_launcher_stops_unless_file_d_and_ts05_give_their_committed_bytes(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "r2", sha_stub=False)[0] == 20
    assert _exit(tmp_path, "r2", "D") == "20" and _exit(tmp_path, "r2", "T") is None
    assert not any("--trend-strength" in c for c in _calls(tmp_path))
    assert _launch(tmp_path, head, "r3", extra={"STUB_TS05": "0" * 64})[0] == 22
    assert _exit(tmp_path, "r3", "D") == "0" and _exit(tmp_path, "r3", "T") == "22"
    assert _exit(tmp_path, "r3", "OTS05a") is None
    assert not any("--captain-ruling R29" in c for c in _calls(tmp_path))


def test_the_launcher_fails_on_a_build_a_check_a_stamp_the_tests_or_two_different_files(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "b1", extra={"STUB_BUILD_EXIT": "1"})[0] == 11
    assert _exit(tmp_path, "b1", "D") == "11"
    assert _launch(tmp_path, head, "b2", extra={"STUB_CHECK_EXIT": "1"})[0] == 13
    assert _launch(tmp_path, head, "b3", extra={"STUB_DIRTY": "1"})[0] == 15
    assert _launch(tmp_path, head, "b4", extra={"STUB_PYTEST_EXIT": "1"})[0] == 10
    assert _exit(tmp_path, "b4", "tests") == "10" and _exit(tmp_path, "b4", "D") is None
    code = _launch(tmp_path, head, "b5", extra={"STUB_NONDET": "1"})[0]
    assert code == 31 and _exit(tmp_path, "b5", "done") == "31"


def test_the_launcher_gate_refusals(tmp_path, repo):
    repo, git, head = repo
    assert _launch(tmp_path, head[:12], "g0")[0] == 90
    assert _launch(tmp_path, head, "a/b")[0] == 90
    assert _launch(tmp_path, head, "it's")[0] == 90
    assert _launch(tmp_path, "d" * 40, "g1")[0] == 96 and _exit(tmp_path, "g1", "gate") == "96"
    (repo / "src" / "a.py").write_text("x = 2\n")
    assert _launch(tmp_path, head, "g2")[0] == 98
    git("checkout", "-q", "--", "src/a.py")
    (repo / "runs").mkdir()
    (repo / "runs" / f"heartlined_{head[:7]}_g3").mkdir()
    assert _launch(tmp_path, head, "g3")[0] == 93
    (repo / "docs" / "S2_RULES.md").write_text("rules, changed\n")
    git("commit", "-q", "-am", "rules changed without a new pin")
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "pin")[0] == 97
    git("revert", "--no-edit", "HEAD")
    t = repo / "experiments" / "harness" / "TRIES.md"
    good = t.read_text()
    bad = {"none": good.replace(R29, ""), "two": good + R29, "other": good.replace("| R29 (", "| R30 ("),
           **{f"drop{i}": good.replace(s, "x") for i, s in enumerate(
               ("`--entry organ_scope`", "`--trend-strength 0.5`", "`--trend-strength 1.0`", "`--captain-ruling R29`",
                "ruling, no run", "the captain's order"))},
           "later_cell": good.replace("`--trend-strength 0.5` and ", "").replace(
               "| heart | x |\n| H28", "| `--trend-strength 0.5` | x |\n| H28")}
    for name, text in bad.items():
        t.write_text(text)
        git("commit", "-q", "-am", f"ledger {name}")
        assert _launch(tmp_path, git("rev-parse", "HEAD"), f"row_{name}")[0] == 95, name
        assert _exit(tmp_path, f"row_{name}", "gate") == "95"
        t.write_text(good)
        git("commit", "-q", "-am", "ledger back")
    assert not any("build_test.py" in c or " -m pytest " in c for c in _calls(tmp_path))
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "ok", extra={"STUB_PYTEST_EXIT": "1"})[0] == 10


# ================================================================================================ as text
def test_the_launcher_deletes_nothing_parses_and_keeps_the_first_lines_conventions():
    sh = LAUNCHER.read_text()
    body = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert not re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|truncate|shred)\s", body)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    assert "python -I -B" in sh and 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'TMPDIR="$S/tmp"' in sh and "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "/mnt/" not in sh and "vec-campaign/" not in body
    assert "S=$HOME/v/scratch/heartlined/$RUN" in sh and "O=runs/heartlined_${head:0:7}_$RUN" in sh
    # PATH is set before the first uv call
    assert body.index("export PATH=$HOME/.local/bin:$PATH") < body.index("uv run")
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    # the shared parts are heart_trend_line.sh's (and so heart_organ.sh's), line for line
    for other in (LINE.read_text(), ORGAN.read_text()):
        for start, end in (("finish() {", "step_ok() {"), ("stamp() {", ': > "$S/stamps"')):
            assert sh[sh.index(start):sh.index(end)] == other[other.index(start):other.index(end)]
    line = LINE.read_text()
    for start, end in (("set -u\n", "export PATH="), ("head=$(git rev-parse HEAD)", "python3 -I -B -c")):
        assert sh[sh.index(start):sh.index(end)] == line[line.index(start):line.index(end)].replace(
            "heart_trend_line.sh", "heart_trend_line_d.sh")


def test_the_launcher_reads_only_exit_sha_and_stamp():
    sh = LAUNCHER.read_text()
    for line in (ln for ln in sh.splitlines() if not ln.lstrip().startswith("#")):
        if line.lstrip().startswith(("$PY experiments/builds/", "uv run --frozen --no-sync vec check")):
            assert "> /dev/null 2>&1" in line, line
    body = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert 'print(s["commit"], s["clean"])' in body and '"size"' not in body
    assert not re.search(r"(^|[\s;&|(])(stat|du|wc|ls|cat|head|tail|less|h5ls)\s", body.replace('cat "$S', ""))
    assert 'D_ARGS="--entry organ_scope --expression copy --z-mode raw"' in sh
    assert 'TS05_ARGS="--entry core_trim --trend-strength 0.5 --z-mode raw --captain-ruling R28"' in sh
    assert 'OTS05_ARGS="--entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29"' in sh
    assert 'OTS10_ARGS="--entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29"' in sh
    assert "BOARD=T2:heart:test_extrap" in sh


def test_the_identity_hashes_are_the_committed_lists():
    """D_SHA and D_HALF_SHA are file D's in experiments/builds/results/heartorgan_72b9662_d1.SHA256SUMS (both of its
    builds), TS05_SHA and TS05_HALF_SHA row R28's TS05 in experiments/builds/results/heartline_d41cd14_ts1.SHA256SUMS
    (both of its builds); the launcher holds exactly these."""
    sh = LAUNCHER.read_text()

    def listed(path):
        out: dict[str, set] = {}
        for ln in path.read_text().splitlines():
            digest, name = ln.split()
            out.setdefault(name.split("/", 1)[1], set()).add(digest)
        return out

    d, ts = listed(D_LIST), listed(TS_LIST)
    assert d["T2heart_organ.h5ad"] == {D_SHA} and d["T2heart_organ_half.h5ad"] == {D_HALF_SHA}
    assert ts["T2heart_core_ts05.h5ad"] == {TS05_SHA} and ts["T2heart_core_ts05_half.h5ad"] == {TS05_HALF_SHA}
    for var, digest in (("D_SHA", D_SHA), ("D_HALF_SHA", D_HALF_SHA), ("TS05_SHA", TS05_SHA),
                        ("TS05_HALF_SHA", TS05_HALF_SHA)):
        assert re.search(rf"^{var}={digest}\b", sh, re.M), var


def test_the_launchers_commands_are_the_builders():
    """The commands the launcher builds are the builder's own: HEART_ORGAN_SCOPE, row R28's HEART_TS_COMMANDS[0.5] and
    row R29's HEART_OTS_COMMANDS."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    sh = LAUNCHER.read_text()
    assert f'D_ARGS="{bt.HEART_ORGAN_SCOPE[len("heart "):]}"' in sh
    assert f'TS05_ARGS="{bt.HEART_TS_COMMANDS[0.5][len("heart "):]}"' in sh
    for s, var in ((0.5, "OTS05_ARGS"), (1.0, "OTS10_ARGS")):
        assert f'{var}="{bt.HEART_OTS_COMMANDS[s][len("heart "):]}"' in sh
