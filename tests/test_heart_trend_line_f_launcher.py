"""experiments/builds/heart_trend_line_f.sh, the test-board builds of file D's three further strengths, files OTS025,
OTS075 and OTS125 (the ruling row given as its third argument). Every call runs in a throwaway git repository under the
test's own folder, with a temporary home and a stand-in uv (and, where named, a stand-in sha256sum) first on PATH: no
launcher runs with a real commit, nothing reads a competition file and nothing is built. The uv stand-in is
tests/test_heart_organ_launchers' own, so the launchers are held to the same calls; the sha256sum stand-in gives the
three identity steps' committed hashes. The gate's ledger check is also held to the builder's own on the same
ledgers."""

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

LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line_f.sh"
D_LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line_d.sh"
D_LIST = REPO_ROOT / "experiments" / "builds" / "results" / "heartorgan_72b9662_d1.SHA256SUMS"
OD_LIST = REPO_ROOT / "experiments" / "builds" / "results" / "heartlined_0d47a67_od1.SHA256SUMS"
D_SHA = "76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b"
D_HALF_SHA = "1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00"
OTS05_SHA = "c3b5fc9c7b47d0d0aafc4e5876bac7f06988cb2d628a0030b4b643806da8bceb"
OTS05_HALF_SHA = "f918e76731aac18eda3b46c8a35c18f26ae00a1a6a1db36d33e7515704cc43ea"
OTS10_SHA = "efaaabe0b744d82dc40135ed9269a6baf8938c9cdd83f85fb3917c0b102e7a2b"
OTS10_HALF_SHA = "ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9"
ROW = "R44"  # a synthetic stand-in for the ruling row
F_ROW = ("| R44 (the captain's order, a synthetic stand-in: ruling, no run, not a try; the heart builder is to take "
         "`--entry organ_scope` with `--trend-strength 0.25`, `--trend-strength 0.75` and `--trend-strength 1.25` only "
         "with `--captain-ruling R44`) | heart | x |\n")
F_STRINGS = ("`--entry organ_scope`", "`--trend-strength 0.25`", "`--trend-strength 0.75`", "`--trend-strength 1.25`",
             "`--captain-ruling R44`")
STEPS = ("gate", "tests", "D", "OTS05", "OTS10", "OTS025a", "OTS025b", "OTS075a", "OTS075b", "OTS125a", "OTS125b",
         "done")
TRIO = {"025": "0.25", "075": "0.75", "125": "1.25"}

# A stand-in sha256sum that gives the committed hashes for file D's, OTS05's and OTS10's files and halves (each one
# replaceable through the environment) and the real hash otherwise.
STUB_SHA = r"""#!/usr/bin/env bash
REAL=__REAL__
if [ "${1:-}" = -c ]; then exec "$REAL" "$@"; fi
for f in "$@"; do
  case "$(basename "$f")" in
    T2heart_organ.h5ad) echo "${STUB_SHA_D:-__SHA_D__}  $f" ;;
    T2heart_organ_half.h5ad) echo "${STUB_SHA_DH:-__SHA_DH__}  $f" ;;
    T2heart_organ_ts05.h5ad) echo "${STUB_SHA_O5:-__SHA_O5__}  $f" ;;
    T2heart_organ_ts05_half.h5ad) echo "${STUB_SHA_O5H:-__SHA_O5H__}  $f" ;;
    T2heart_organ_ts10.h5ad) echo "${STUB_SHA_O10:-__SHA_O10__}  $f" ;;
    T2heart_organ_ts10_half.h5ad) echo "${STUB_SHA_O10H:-__SHA_O10H__}  $f" ;;
    *) "$REAL" "$f" ;;
  esac
done
"""
# The shared uv stand-in, with its non-determinism narrowed to the builds whose arguments hold STUB_NONDET_ON (all of
# them when it is unset, as the stand-in itself has it).
NONDET = 'noise=""; [ -n "${STUB_NONDET:-}" ] && noise="$RANDOM$RANDOM$(date +%N)"'
STUB_UV_F = STUB_UV.replace(NONDET, 'noise=""; [ -n "${STUB_NONDET:-}" ] && case "$*" in *"${STUB_NONDET_ON:-}"*) '
                                    'noise="$RANDOM$RANDOM$(date +%N)" ;; esac')


def _bin(tmp_path, sha_stub=False):
    b = tmp_path / ("bin_sha" if sha_stub else "bin")
    b.mkdir(exist_ok=True)
    (b / "uv").write_text(STUB_UV_F)
    (b / "uv").chmod(0o755)
    if sha_stub:
        text = STUB_SHA.replace("__REAL__", shutil.which("sha256sum"))
        for key, digest in (("__SHA_D__", D_SHA), ("__SHA_DH__", D_HALF_SHA), ("__SHA_O5__", OTS05_SHA),
                            ("__SHA_O5H__", OTS05_HALF_SHA), ("__SHA_O10__", OTS10_SHA),
                            ("__SHA_O10H__", OTS10_HALF_SHA)):
            text = text.replace(key, digest)
        assert "__SHA_" not in text
        (b / "sha256sum").write_text(text)
        (b / "sha256sum").chmod(0o755)
    return b


@pytest.fixture()
def repo(tmp_path):
    """A throwaway repository: the rules file and its pin, the ledger with row R44 as the captain's order for the three
    strengths, a stand-in src/."""
    repo = tmp_path / "repo"
    for d in ("docs", "experiments/harness", "src"):
        (repo / d).mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (repo / "docs" / "S2_RULES.md").write_text("rules\n")
    sha = subprocess.run(["sha256sum", "docs/S2_RULES.md"], cwd=repo, capture_output=True, text=True).stdout
    (repo / "docs" / "S2_RULES.sha256").write_text(sha)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(
        "| R28 (other) | heart | x |\n| R29 (other) | heart | x |\n" + F_ROW + "| H28 (row) | heart | x |\n")
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


def _scratch(tmp_path, run):
    return tmp_path / "home" / "v" / "scratch" / "heartf" / run


def _exit(tmp_path, run, step):
    p = _scratch(tmp_path, run) / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


def test_the_stand_in_uv_is_the_shared_one_with_a_narrower_noise():
    assert NONDET in STUB_UV and NONDET not in STUB_UV_F and STUB_UV_F.count("STUB_NONDET_ON") == 1


def test_the_launcher_end_to_end_with_stand_ins(tmp_path, repo):
    repo, git, head = repo
    code, out, err = _launch(tmp_path, head, "r1", ROW)
    assert code == 0, err
    assert out == ""
    for step in STEPS:
        assert _exit(tmp_path, "r1", step) == "0", step
    S = _scratch(tmp_path, "r1")
    rel = Path("runs") / f"heartf_{head[:7]}_r1"
    assert (S / "head").read_text().strip() == head and (S / "row").read_text().strip() == ROW
    assert (S / "stamps").read_text().splitlines() == [
        f"D/T2heart_organ {head} True", f"OTS05/T2heart_organ_ts05 {head} True",
        f"OTS10/T2heart_organ_ts10 {head} True",
        *(f"OTS{k}{x}/T2heart_organ_ts{k} {head} True" for k in TRIO for x in "ab")]
    sums = (S / "SHA256SUMS").read_text()
    assert sums == (repo / rel / "SHA256SUMS").read_text()
    assert [ln.split()[1] for ln in sums.splitlines()] == [
        "D/T2heart_organ.h5ad", "D/T2heart_organ_half.h5ad",
        "OTS05/T2heart_organ_ts05.h5ad", "OTS05/T2heart_organ_ts05_half.h5ad",
        "OTS10/T2heart_organ_ts10.h5ad", "OTS10/T2heart_organ_ts10_half.h5ad",
        *(f"OTS{k}{x}/T2heart_organ_ts{k}{h}.h5ad" for k in TRIO for x in "ab" for h in ("", "_half"))]
    digests = {ln.split()[1]: ln.split()[0] for ln in sums.splitlines()}
    assert (digests["D/T2heart_organ.h5ad"], digests["OTS05/T2heart_organ_ts05.h5ad"],
            digests["OTS10/T2heart_organ_ts10_half.h5ad"]) == (D_SHA, OTS05_SHA, OTS10_HALF_SHA)
    for k in TRIO:
        assert digests[f"OTS{k}a/T2heart_organ_ts{k}.h5ad"] == digests[f"OTS{k}b/T2heart_organ_ts{k}.h5ad"]
    calls = _calls(tmp_path)
    py = "run --frozen --no-sync python -I -B experiments/builds/build_test.py heart"
    assert [c for c in calls if "build_test.py" in c] == [
        f"{py} --entry organ_scope --expression copy --z-mode raw --out {rel}/D/T2heart_organ.h5ad",
        f"{py} --entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29 "
        f"--out {rel}/OTS05/T2heart_organ_ts05.h5ad",
        f"{py} --entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29 "
        f"--out {rel}/OTS10/T2heart_organ_ts10.h5ad",
        *(f"{py} --entry organ_scope --trend-strength {v} --z-mode raw --captain-ruling {ROW} "
          f"--out {rel}/OTS{k}{x}/T2heart_organ_ts{k}.h5ad" for k, v in TRIO.items() for x in "ab")]
    halves = [c for c in calls if "experiments/builds/halve.py" in c]
    assert len(halves) == 9 and all("--seed" not in c for c in halves)  # halve.py's default seed, 1000001
    assert halves[3].endswith(f"{rel}/OTS025a/T2heart_organ_ts025.h5ad {rel}/OTS025a/T2heart_organ_ts025_half.h5ad")
    checks = [c for c in calls if " vec check " in c]
    assert len(checks) == 18 and all(c.endswith("--board T2:heart:test_extrap") for c in checks)
    tests = [c for c in calls if " -m pytest " in c]
    assert len(tests) == 1
    for f in ("test_heart_trend_line_f.py", "test_heart_trend_line_f_launcher.py", "test_heart_trend_line_d.py",
              "test_heart_trend_line_d_launcher.py", "test_heart_trend_line.py", "test_heart_trend_line_launcher.py",
              "test_build_heart_organ_scope.py", "test_heart_organ_launchers.py", "test_h28_heart_trend_strength.py",
              "test_halve.py", "test_hidden_counts.py", "test_s2_rules.py"):
        assert f"tests/{f}" in tests[0], f
    # a run folder is never reused
    assert _launch(tmp_path, head, "r1", ROW)[0] == 91


@pytest.mark.parametrize("env, code, step", [
    ({"STUB_SHA_D": "0" * 64}, 20, "D"), ({"STUB_SHA_DH": "0" * 64}, 21, "D"),
    ({"STUB_SHA_O5": "0" * 64}, 22, "OTS05"), ({"STUB_SHA_O5H": "0" * 64}, 23, "OTS05"),
    ({"STUB_SHA_O10": "0" * 64}, 24, "OTS10"), ({"STUB_SHA_O10H": "0" * 64}, 25, "OTS10")])
def test_the_launcher_stops_unless_d_ots05_and_ots10_give_their_committed_bytes(tmp_path, repo, env, code, step):
    _, _, head = repo
    assert _launch(tmp_path, head, "s", ROW, extra=env)[0] == code
    assert _exit(tmp_path, "s", step) == str(code)
    later = STEPS[STEPS.index(step) + 1:]
    assert all(_exit(tmp_path, "s", s) is None for s in later)
    assert not any(f"--captain-ruling {ROW}" in c for c in _calls(tmp_path))  # none of the three is built


def test_the_launcher_stops_at_file_d_on_its_real_hash(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "r2", ROW, sha_stub=False)[0] == 20
    assert _exit(tmp_path, "r2", "D") == "20" and _exit(tmp_path, "r2", "OTS05") is None
    assert not any("--trend-strength" in c for c in _calls(tmp_path))


def test_the_launcher_fails_on_a_build_a_check_a_stamp_the_tests_or_two_different_files(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "b1", ROW, extra={"STUB_BUILD_EXIT": "1"})[0] == 11
    assert _exit(tmp_path, "b1", "D") == "11"
    assert _launch(tmp_path, head, "b2", ROW, extra={"STUB_CHECK_EXIT": "1"})[0] == 13
    assert _launch(tmp_path, head, "b3", ROW, extra={"STUB_DIRTY": "1"})[0] == 15
    assert _launch(tmp_path, head, "b4", ROW, extra={"STUB_PYTEST_EXIT": "1"})[0] == 10
    assert _exit(tmp_path, "b4", "tests") == "10" and _exit(tmp_path, "b4", "D") is None
    for run, on, code in (("b5", "", 31), ("b6", "--trend-strength 0.75 ", 32), ("b7", "--trend-strength 1.25 ", 33)):
        assert _launch(tmp_path, head, run, ROW, extra={"STUB_NONDET": "1", "STUB_NONDET_ON": on})[0] == code, run
        assert _exit(tmp_path, run, "done") == str(code)
        assert all(_exit(tmp_path, run, s) == "0" for s in STEPS[:-1])


def test_the_launcher_gate_refusals(tmp_path, repo):
    repo, git, head = repo
    assert _launch(tmp_path, head[:12], "g0", ROW)[0] == 90
    assert _launch(tmp_path, head, "a/b", ROW)[0] == 90
    assert _launch(tmp_path, head, "it's", ROW)[0] == 90
    for i, bad in enumerate(("r44", "R-F", "44", "R44 ", " R44", "R44;x", "R", "R44ab", "R44'", "R4.4", "R44`",
                             "R44 --board T2:heart:val_extrap", "$(id)", "*")):
        assert _launch(tmp_path, head, f"w{i}", bad)[0] == 90, bad
        assert not _scratch(tmp_path, f"w{i}").exists(), bad
    for args in ((), (head,), (head, "m2"), (head, "m3", "")):
        code = _launch(tmp_path, *args)[0]
        assert code != 0 and not (tmp_path / "home" / "v" / "scratch" / "heartf").exists(), args
    assert _launch(tmp_path, "d" * 40, "g1", ROW)[0] == 96 and _exit(tmp_path, "g1", "gate") == "96"
    (repo / "src" / "a.py").write_text("x = 2\n")
    assert _launch(tmp_path, head, "g2", ROW)[0] == 98
    git("checkout", "-q", "--", "src/a.py")
    (repo / "runs").mkdir()
    (repo / "runs" / f"heartf_{head[:7]}_g3").mkdir()
    assert _launch(tmp_path, head, "g3", ROW)[0] == 93
    (repo / "docs" / "S2_RULES.md").write_text("rules, changed\n")
    git("commit", "-q", "-am", "rules changed without a new pin")
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "pin", ROW)[0] == 97
    git("revert", "--no-edit", "HEAD")
    t = repo / "experiments" / "harness" / "TRIES.md"
    good = t.read_text()
    bad = {"none": good.replace(F_ROW, ""), "two": good + F_ROW, "other": good.replace("| R44 (", "| R45 ("),
           "longer": good.replace("| R44 (", "| R440 ("),
           "ruled": good.replace("`--captain-ruling R44`", "`--captain-ruling R45`"),
           **{f"drop{i}": good.replace(s, "x") for i, s in enumerate((*F_STRINGS, "ruling, no run",
                                                                       "the captain's order"))},
           "later_cell": good.replace("`--trend-strength 0.75` and ", "").replace(
               "| heart | x |\n| H28", "| `--trend-strength 0.75` | x |\n| H28"),
           "wider": good.replace("`--trend-strength 1.25`", "`--trend-strength 1.250`")}
    for name, text in bad.items():
        assert text != good, name
        t.write_text(text)
        git("commit", "-q", "-am", f"ledger {name}")
        assert _launch(tmp_path, git("rev-parse", "HEAD"), f"row_{name}", ROW)[0] == 95, name
        assert _exit(tmp_path, f"row_{name}", "gate") == "95"
        t.write_text(good)
        git("commit", "-q", "-am", "ledger back")
    # the good ledger with another row given, R29's or one it does not hold
    for i, other in enumerate(("R45", "R29", "R4")):
        assert _launch(tmp_path, git("rev-parse", "HEAD"), f"given{i}", other)[0] == 95, other
    assert not any("build_test.py" in c or " -m pytest " in c for c in _calls(tmp_path))
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "ok", ROW, extra={"STUB_PYTEST_EXIT": "1"})[0] == 10


def test_the_gates_row_check_is_the_builders(tmp_path, monkeypatch):
    """The gate's inline check and the builder's _heart_ftrio_ordered take and refuse the same ledgers."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    sh = LAUNCHER.read_text()
    code = sh[sh.index("python3 -I -B -c '") + len("python3 -I -B -c '"):]
    code = code[:code.index("' \\\n")]
    base = "| R28 (other) | heart | x |\n" + F_ROW + "| H28 (row) | heart | x |\n"
    ledgers = {"ok": base, "none": "| R27 (x) | T1 | x |\n", "two": base + F_ROW,
               "other": base.replace("| R44 (", "| R45 ("), "longer": base.replace("| R44 (", "| R440 ("),
               "ruled": base.replace("`--captain-ruling R44`", "`--captain-ruling R45`"),
               "wider": base.replace("`--trend-strength 0.25`", "`--trend-strength 0.250`"),
               **{f"drop{i}": base.replace(s, "x") for i, s in enumerate((*F_STRINGS, "ruling, no run",
                                                                           "the captain's order"))}}
    seen = set()
    for name, text in ledgers.items():
        path = tmp_path / f"{name}.md"
        path.write_text(text)
        monkeypatch.setattr(bt, "TRIES_FILE", path)
        for row in (ROW, "R45", "R440"):
            gate = subprocess.run(["python3", "-I", "-B", "-c", code, str(path), row]).returncode == 0
            assert gate == bt._heart_ftrio_ordered(row), (name, row)
            seen.add(gate)
    assert seen == {True, False}


# ================================================================================================ as text
def _body(sh: str) -> str:
    return "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))


def test_the_launcher_deletes_nothing_parses_and_keeps_the_d_launchers_conventions():
    sh, d = LAUNCHER.read_text(), D_LAUNCHER.read_text()
    body = _body(sh)
    assert not re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|truncate|shred)\s", body)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    assert len(re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh)) == 3
    assert "python -I -B" in sh and 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'TMPDIR="$S/tmp"' in sh and "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "/mnt/" not in sh and "vec-campaign/" not in body
    assert "S=$HOME/v/scratch/heartf/$RUN" in sh and "O=runs/heartf_${head:0:7}_$RUN" in sh
    # PATH is set before the first uv call, and ROW is checked before it is used
    assert body.index("export PATH=$HOME/.local/bin:$PATH") < body.index("uv run")
    assert body.index('[[ "$ROW" =~ ^R[0-9]+[a-z]?$ ]]') < body.index('mkdir -p "$HOME')
    assert body.index('[[ "$ROW" =~ ^R[0-9]+[a-z]?$ ]]') < body.index("$ROW\"\n")
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    # the shared parts are heart_trend_line_d.sh's, line for line: (start, this launcher's end, the D launcher's end)
    for start, mine_end, their_end in (("finish() {", "step_ok() {", "step_ok() {"),
                                       ("stamp() {", ': > "$S/stamps"', ': > "$S/stamps"'),
                                       ("head=$(git rev-parse HEAD)", "python3 -I -B -c", "python3 -I -B -c"),
                                       ("mkdir -p runs\n", 'echo "$head"', 'echo "$head"'),
                                       ("STEP=D\n", "STEP=OTS05\n", "STEP=T\n"),
                                       ("STEP=done\n", "(cd ", "(cd ")):
        i, j = sh.index(start), d.index(start)
        assert sh[i:sh.index(mine_end, i)] == d[j:d.index(their_end, j)], start
    assert sh.count("build() {") == d.count("build() {") == 1
    # the arguments' block is the D launcher's with ROW added
    mine = sh[sh.index("set -u\n"):sh.index("export PATH=")]
    theirs = d[d.index("set -u\n"):d.index("export PATH=")].replace("heart_trend_line_d.sh COMMIT RUN",
                                                                     "heart_trend_line_f.sh COMMIT RUN ROW")
    added = [ln for ln in mine.splitlines() if ln not in theirs.splitlines()]
    assert [ln for ln in mine.splitlines() if ln in theirs.splitlines()] == theirs.splitlines()
    assert added == ["ROW=${3:?usage: heart_trend_line_f.sh COMMIT RUN ROW}",
                     '[[ "$ROW" =~ ^R[0-9]+[a-z]?$ ]] || { echo "refused: ROW is a ruling row id, R and its number '
                     '(e.g. R44)" >&2; exit 90; }']


def test_the_launcher_reads_only_exit_sha_and_stamp():
    sh = LAUNCHER.read_text()
    for line in (ln for ln in sh.splitlines() if not ln.lstrip().startswith("#")):
        if line.lstrip().startswith(("$PY experiments/builds/", "uv run --frozen --no-sync vec check")):
            assert "> /dev/null 2>&1" in line, line
    body = _body(sh)
    assert 'print(s["commit"], s["clean"])' in body and '"size"' not in body
    assert not re.search(r"(^|[\s;&|(])(stat|du|wc|ls|cat|head|tail|less|h5ls)\s", body)
    assert 'D_ARGS="--entry organ_scope --expression copy --z-mode raw"' in sh
    assert 'OTS05_ARGS="--entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29"' in sh
    assert 'OTS10_ARGS="--entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29"' in sh
    for k, v in TRIO.items():
        assert f'OTS{k}_ARGS="--entry organ_scope --trend-strength {v} --z-mode raw --captain-ruling $ROW"' in sh
    assert "BOARD=T2:heart:test_extrap" in sh


def test_the_identity_hashes_are_the_committed_lists():
    """D_SHA and D_HALF_SHA are file D's in experiments/builds/results/heartorgan_72b9662_d1.SHA256SUMS (both of its
    builds) and in heartlined_0d47a67_od1.SHA256SUMS; OTS05's and OTS10's are row R29's in
    heartlined_0d47a67_od1.SHA256SUMS (both of their builds); the launcher holds exactly these."""
    sh = LAUNCHER.read_text()

    def listed(path):
        out: dict[str, set] = {}
        for ln in path.read_text().splitlines():
            digest, name = ln.split()
            out.setdefault(name.split("/", 1)[1], set()).add(digest)
        return out

    d, od = listed(D_LIST), listed(OD_LIST)
    assert d["T2heart_organ.h5ad"] == od["T2heart_organ.h5ad"] == {D_SHA}
    assert d["T2heart_organ_half.h5ad"] == od["T2heart_organ_half.h5ad"] == {D_HALF_SHA}
    assert od["T2heart_organ_ts05.h5ad"] == {OTS05_SHA} and od["T2heart_organ_ts05_half.h5ad"] == {OTS05_HALF_SHA}
    assert od["T2heart_organ_ts10.h5ad"] == {OTS10_SHA} and od["T2heart_organ_ts10_half.h5ad"] == {OTS10_HALF_SHA}
    for var, digest in (("D_SHA", D_SHA), ("D_HALF_SHA", D_HALF_SHA), ("OTS05_SHA", OTS05_SHA),
                        ("OTS05_HALF_SHA", OTS05_HALF_SHA), ("OTS10_SHA", OTS10_SHA),
                        ("OTS10_HALF_SHA", OTS10_HALF_SHA)):
        assert re.search(rf"^{var}={digest}\b", sh, re.M), var


def test_the_launchers_commands_are_the_builders():
    """The commands the launcher builds are the builder's own: HEART_ORGAN_SCOPE, row R29's HEART_OTS_COMMANDS and the
    three further strengths' HEART_FTRIO_COMMANDS with the row given."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    sh = LAUNCHER.read_text()
    assert f'D_ARGS="{bt.HEART_ORGAN_SCOPE[len("heart "):]}"' in sh
    for s, var in ((0.5, "OTS05_ARGS"), (1.0, "OTS10_ARGS")):
        assert f'{var}="{bt.HEART_OTS_COMMANDS[s][len("heart "):]}"' in sh
    for s, var in ((0.25, "OTS025_ARGS"), (0.75, "OTS075_ARGS"), (1.25, "OTS125_ARGS")):
        assert f'{var}="{bt.HEART_FTRIO_COMMANDS[s].format(row="$ROW")[len("heart "):]}"' in sh
    assert bt.HEART_FTRIO_ROW_FORM.pattern == "R[0-9]+[a-z]?" and "^R[0-9]+[a-z]?$" in sh
