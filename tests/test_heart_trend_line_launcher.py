"""experiments/builds/heart_trend_line.sh, the test-board builds of the heart search line's two files TS05 and TS10
(row R28 of experiments/harness/TRIES.md). Every call runs in a throwaway git repository under the test's own folder,
with a temporary home and a stand-in uv (and, where named, a stand-in sha256sum) first on PATH: no launcher runs with a
real commit, nothing reads a competition file and nothing is built. The stand-ins are tests/test_heart_organ_launchers'
own, so the two launchers are held to the same calls."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_heart_organ_launchers import C_HALF_SHA, C_SHA, HEART3, _bin, _calls, _env, _git_fn  # noqa: E402

LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line.sh"
ORGAN = REPO_ROOT / "experiments" / "builds" / "heart_organ.sh"
R28 = ("| R28 (the captain's order of 2026-11-23, \"start the heart search\": a synthetic stand-in; ruling, no run, "
       "not a try; the heart builder is to take `--trend-strength 0.5` and `--trend-strength 1.0` only with "
       "`--captain-ruling R28`) | heart | x |\n")
STEPS = ("gate", "tests", "C", "TS05a", "TS05b", "TS10a", "TS10b", "done")


@pytest.fixture()
def repo(tmp_path):
    """A throwaway repository: the rules file and its pin, the ledger with row R28 as the captain's order, a stand-in
    src/."""
    repo = tmp_path / "repo"
    for d in ("docs", "experiments/harness", "src"):
        (repo / d).mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (repo / "docs" / "S2_RULES.md").write_text("rules\n")
    sha = subprocess.run(["sha256sum", "docs/S2_RULES.md"], cwd=repo, capture_output=True, text=True).stdout
    (repo / "docs" / "S2_RULES.sha256").write_text(sha)
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| R27 (other) | T1 | x |\n" + R28
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
    p = tmp_path / "home" / "v" / "scratch" / "heartline" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


def test_the_launcher_end_to_end_with_stand_ins(tmp_path, repo):
    repo, git, head = repo
    code, out, err = _launch(tmp_path, head, "r1")
    assert code == 0, err
    assert out == ""
    for step in STEPS:
        assert _exit(tmp_path, "r1", step) == "0", step
    S = tmp_path / "home" / "v" / "scratch" / "heartline" / "r1"
    rel = Path("runs") / f"heartline_{head[:7]}_r1"
    assert (S / "head").read_text().strip() == head
    assert (S / "stamps").read_text().splitlines() == [f"C/T2heart_s2core {head} True",
                                                       f"TS05a/T2heart_core_ts05 {head} True",
                                                       f"TS05b/T2heart_core_ts05 {head} True",
                                                       f"TS10a/T2heart_core_ts10 {head} True",
                                                       f"TS10b/T2heart_core_ts10 {head} True"]
    sums = (S / "SHA256SUMS").read_text()
    assert sums == (repo / rel / "SHA256SUMS").read_text()
    assert [ln.split()[1] for ln in sums.splitlines()] == [
        "C/T2heart_s2core.h5ad", "C/T2heart_s2core_half.h5ad",
        "TS05a/T2heart_core_ts05.h5ad", "TS05a/T2heart_core_ts05_half.h5ad",
        "TS05b/T2heart_core_ts05.h5ad", "TS05b/T2heart_core_ts05_half.h5ad",
        "TS10a/T2heart_core_ts10.h5ad", "TS10a/T2heart_core_ts10_half.h5ad",
        "TS10b/T2heart_core_ts10.h5ad", "TS10b/T2heart_core_ts10_half.h5ad"]
    calls = _calls(tmp_path)
    py = "run --frozen --no-sync python -I -B experiments/builds/build_test.py heart"
    assert [c for c in calls if "build_test.py" in c] == [
        f"{py} --entry core_trim --expression copy --z-mode raw --out {rel}/C/T2heart_s2core.h5ad",
        *(f"{py} --entry core_trim --trend-strength 0.5 --z-mode raw --captain-ruling R28 "
          f"--out {rel}/{d}/T2heart_core_ts05.h5ad" for d in ("TS05a", "TS05b")),
        *(f"{py} --entry core_trim --trend-strength 1.0 --z-mode raw --captain-ruling R28 "
          f"--out {rel}/{d}/T2heart_core_ts10.h5ad" for d in ("TS10a", "TS10b"))]
    halves = [c for c in calls if "experiments/builds/halve.py" in c]
    assert len(halves) == 5 and all("--seed" not in c for c in halves)  # halve.py's default seed, 1000001
    assert halves[1].endswith(f"{rel}/TS05a/T2heart_core_ts05.h5ad {rel}/TS05a/T2heart_core_ts05_half.h5ad")
    checks = [c for c in calls if " vec check " in c]
    assert len(checks) == 10 and all(c.endswith("--board T2:heart:test_extrap") for c in checks)
    tests = [c for c in calls if " -m pytest " in c]
    assert len(tests) == 1 and "tests/test_heart_trend_line.py" in tests[0]
    assert "tests/test_heart_trend_line_launcher.py" in tests[0]
    assert "tests/test_h28_heart_trend_strength.py" in tests[0]
    # a run folder is never reused
    assert _launch(tmp_path, head, "r1")[0] == 91


def test_the_launcher_stops_unless_the_full_trim_gives_set_heart3(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "r2", sha_stub=False)[0] == 20
    assert _exit(tmp_path, "r2", "C") == "20" and _exit(tmp_path, "r2", "TS05a") is None
    assert not any("--trend-strength" in c for c in _calls(tmp_path))


def test_the_launcher_fails_on_a_build_a_check_a_stamp_the_tests_or_two_different_files(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "b1", extra={"STUB_BUILD_EXIT": "1"})[0] == 11
    assert _exit(tmp_path, "b1", "C") == "11"
    assert _launch(tmp_path, head, "b2", extra={"STUB_CHECK_EXIT": "1"})[0] == 13
    assert _launch(tmp_path, head, "b3", extra={"STUB_DIRTY": "1"})[0] == 15
    assert _launch(tmp_path, head, "b4", extra={"STUB_PYTEST_EXIT": "1"})[0] == 10
    assert _exit(tmp_path, "b4", "tests") == "10" and _exit(tmp_path, "b4", "C") is None
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
    (repo / "runs" / f"heartline_{head[:7]}_g3").mkdir()
    assert _launch(tmp_path, head, "g3")[0] == 93
    (repo / "docs" / "S2_RULES.md").write_text("rules, changed\n")
    git("commit", "-q", "-am", "rules changed without a new pin")
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "pin")[0] == 97
    git("revert", "--no-edit", "HEAD")
    t = repo / "experiments" / "harness" / "TRIES.md"
    good = t.read_text()
    bad = {"none": good.replace(R28, ""), "two": good + R28, "other": good.replace("| R28 (", "| R29 ("),
           **{f"drop{i}": good.replace(s, "x") for i, s in enumerate(
               ("`--trend-strength 0.5`", "`--trend-strength 1.0`", "`--captain-ruling R28`", "ruling, no run",
                "the captain's order"))},
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
def test_the_launcher_deletes_nothing_parses_and_keeps_heart_organs_conventions():
    sh = LAUNCHER.read_text()
    body = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert not re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|truncate|shred)\s", body)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    assert "python -I -B" in sh and 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'TMPDIR="$S/tmp"' in sh and "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "/mnt/" not in sh and "vec-campaign/" not in body
    assert "S=$HOME/v/scratch/heartline/$RUN" in sh and "O=runs/heartline_${head:0:7}_$RUN" in sh
    # PATH is set before the first uv call
    assert body.index("export PATH=$HOME/.local/bin:$PATH") < body.index("uv run")
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    # the shared parts are heart_organ.sh's, line for line
    organ = ORGAN.read_text()
    for start, end in (("finish() {", "step_ok() {"), ("stamp() {", ': > "$S/stamps"')):
        assert sh[sh.index(start):sh.index(end)] == organ[organ.index(start):organ.index(end)]


def test_the_launcher_reads_only_exit_sha_and_stamp():
    sh = LAUNCHER.read_text()
    for line in (ln for ln in sh.splitlines() if not ln.lstrip().startswith("#")):
        if line.lstrip().startswith(("$PY experiments/builds/", "uv run --frozen --no-sync vec check")):
            assert "> /dev/null 2>&1" in line, line
    body = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert 'print(s["commit"], s["clean"])' in body and '"size"' not in body
    assert not re.search(r"(^|[\s;&|(])(stat|du|wc|ls|cat|head|tail|less|h5ls)\s", body.replace('cat "$S', ""))
    # the identity step's hashes are set heart3's committed DGX list, as heart_organ.sh holds them
    listed = dict(reversed(ln.split()) for ln in HEART3.read_text().splitlines())
    assert listed["builds/T2heart_s2core.h5ad"] == C_SHA and listed["builds/T2heart_s2core_half.h5ad"] == C_HALF_SHA
    assert f"C_SHA={C_SHA}" in sh and f"C_HALF_SHA={C_HALF_SHA}" in sh
    assert f"C_SHA={C_SHA}" in ORGAN.read_text()
    assert 'C_ARGS="--entry core_trim --expression copy --z-mode raw"' in sh
    assert 'TS05_ARGS="--entry core_trim --trend-strength 0.5 --z-mode raw --captain-ruling R28"' in sh
    assert 'TS10_ARGS="--entry core_trim --trend-strength 1.0 --z-mode raw --captain-ruling R28"' in sh
    assert "BOARD=T2:heart:test_extrap" in sh


def test_the_launchers_commands_are_the_builders(tmp_path):
    """The two commands the launcher builds are the builder's own HEART_TS_COMMANDS and C's HEART_CORE_TRIM."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    sh = LAUNCHER.read_text()
    for s, var in ((0.5, "TS05_ARGS"), (1.0, "TS10_ARGS")):
        assert f'{var}="{bt.HEART_TS_COMMANDS[s][len("heart "):]}"' in sh
    assert f'C_ARGS="{bt.HEART_CORE_TRIM[len("heart "):]}"' in sh
