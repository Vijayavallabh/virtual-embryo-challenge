"""experiments/builds/heart_trend_line_g.sh, the test-board builds of file D's strength 1.5, file OTS15 (row R52, the
captain's order of 2026-11-24; the row fixed, no ROW argument). Every call runs in a throwaway git repository under
the test's own folder, with a temporary home and a stand-in uv (and, where named, a stand-in sha256sum) first on PATH:
no launcher runs with a real commit, nothing reads a competition file and nothing is built. The uv stand-in is
tests/test_heart_trend_line_f_launcher's (tests/test_heart_organ_launchers' own with its noise narrowed), so the
launchers are held to the same calls; the sha256sum stand-in gives the five identity steps' committed hashes. The
gate's ledger check is also held to the builder's own on the same ledgers."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_heart_organ_launchers import _calls, _env, _git_fn  # noqa: E402
from test_heart_trend_line_f_launcher import STUB_UV_F  # noqa: E402

LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line_g.sh"
F_LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line_f.sh"
RESULTS = REPO_ROOT / "experiments" / "builds" / "results"
D_LIST = RESULTS / "heartorgan_72b9662_d1.SHA256SUMS"
OD_LIST = RESULTS / "heartlined_0d47a67_od1.SHA256SUMS"
F_LIST = RESULTS / "heartf_ee53e77_f1.SHA256SUMS"
# the five identity steps: step -> (file name, its sha256, its half's sha256), as the committed lists give them
IDENTITY = {
    "D": ("T2heart_organ", "76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b",
          "1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00"),
    "OTS05": ("T2heart_organ_ts05", "c3b5fc9c7b47d0d0aafc4e5876bac7f06988cb2d628a0030b4b643806da8bceb",
              "f918e76731aac18eda3b46c8a35c18f26ae00a1a6a1db36d33e7515704cc43ea"),
    "OTS10": ("T2heart_organ_ts10", "efaaabe0b744d82dc40135ed9269a6baf8938c9cdd83f85fb3917c0b102e7a2b",
              "ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9"),
    "OTS075": ("T2heart_organ_ts075", "109f7f26966c4a3cd0cd18807032ea8eb8f7786df833cd7b3064a21298b1cdbb",
               "1a66b5a5a9d5388d806f548d5c51570345f02e340a6ef2d4ce3490d1e2673f5d"),
    "OTS125": ("T2heart_organ_ts125", "79d0b813db758ce88bda22ed1f8c8b58a560214360bfbfc2795c7eee372419da",
               "adc3243665949f9802056216d86443e2ea0d45d07534013f99e26afb7e57ade0"),
}
STOP_CODES = {"D": (20, 21), "OTS05": (22, 23), "OTS10": (24, 25), "OTS075": (26, 27), "OTS125": (28, 29)}
ARGS = {"D": "--entry organ_scope --expression copy --z-mode raw",
        "OTS05": "--entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29",
        "OTS10": "--entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29",
        "OTS075": "--entry organ_scope --trend-strength 0.75 --z-mode raw --captain-ruling R49",
        "OTS125": "--entry organ_scope --trend-strength 1.25 --z-mode raw --captain-ruling R49",
        "OTS15": "--entry organ_scope --trend-strength 1.5 --z-mode raw --captain-ruling R52"}
R52_ROW = ("| R52 (the captain's order of 2026-11-24, ruling, no run: a synthetic stand-in; the heart builder is to "
           "take `--entry organ_scope` with `--trend-strength 1.5` only with `--captain-ruling R52`; not a try) "
           "| heart | x |\n")
R52_STRINGS = ("`--entry organ_scope`", "`--trend-strength 1.5`", "`--captain-ruling R52`")
STEPS = ("gate", "tests", "D", "OTS05", "OTS10", "OTS075", "OTS125", "OTS15a", "OTS15b", "done")

# A stand-in sha256sum that gives the committed hashes for the five identity steps' files and halves (each one
# replaceable through the environment: STUB_SHA_<STEP> and STUB_SHA_<STEP>_HALF) and the real hash otherwise.
STUB_SHA_HEAD = r"""#!/usr/bin/env bash
REAL=__REAL__
if [ "${1:-}" = -c ]; then exec "$REAL" "$@"; fi
for f in "$@"; do
  case "$(basename "$f")" in
"""
STUB_SHA_TAIL = r"""    *) "$REAL" "$f" ;;
  esac
done
"""


def _stub_sha() -> str:
    cases = "".join(f'    {name}.h5ad) echo "${{STUB_SHA_{step}:-{full}}}  $f" ;;\n'
                    f'    {name}_half.h5ad) echo "${{STUB_SHA_{step}_HALF:-{half}}}  $f" ;;\n'
                    for step, (name, full, half) in IDENTITY.items())
    return STUB_SHA_HEAD.replace("__REAL__", shutil.which("sha256sum")) + cases + STUB_SHA_TAIL


def _bin(tmp_path, sha_stub=False):
    b = tmp_path / ("bin_sha" if sha_stub else "bin")
    b.mkdir(exist_ok=True)
    (b / "uv").write_text(STUB_UV_F)
    (b / "uv").chmod(0o755)
    if sha_stub:
        (b / "sha256sum").write_text(_stub_sha())
        (b / "sha256sum").chmod(0o755)
    return b


@pytest.fixture()
def repo(tmp_path):
    """A throwaway repository: the rules file and its pin, the ledger with row R52 as the captain's order for the file,
    a stand-in src/."""
    repo = tmp_path / "repo"
    for d in ("docs", "experiments/harness", "src"):
        (repo / d).mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (repo / "docs" / "S2_RULES.md").write_text("rules\n")
    sha = subprocess.run(["sha256sum", "docs/S2_RULES.md"], cwd=repo, capture_output=True, text=True).stdout
    (repo / "docs" / "S2_RULES.sha256").write_text(sha)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(
        "| R29 (other) | heart | x |\n| R49 (other) | heart | x |\n" + R52_ROW + "| H34 (row) | heart | x |\n")
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
    return tmp_path / "home" / "v" / "scratch" / "heartg" / run


def _exit(tmp_path, run, step):
    p = _scratch(tmp_path, run) / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


def test_the_stand_in_sha256sum_gives_the_committed_hashes(tmp_path):
    b = _bin(tmp_path, sha_stub=True)
    for step, (name, full, half) in IDENTITY.items():
        for f, digest in ((f"x/{name}.h5ad", full), (f"x/{name}_half.h5ad", half)):
            out = subprocess.run([str(b / "sha256sum"), f], capture_output=True, text=True).stdout
            assert out == f"{digest}  {f}\n", step


def test_the_launcher_end_to_end_with_stand_ins(tmp_path, repo):
    repo, git, head = repo
    code, out, err = _launch(tmp_path, head, "r1")
    assert code == 0, err
    assert out == "" and err == ""
    for step in STEPS:
        assert _exit(tmp_path, "r1", step) == "0", step
    S = _scratch(tmp_path, "r1")
    rel = Path("runs") / f"heartg_{head[:7]}_r1"
    assert (S / "head").read_text().strip() == head and not (S / "row").exists()
    assert (S / "stamps").read_text().splitlines() == [
        *(f"{step}/{name} {head} True" for step, (name, _, _) in IDENTITY.items()),
        f"OTS15a/T2heart_organ_ts15 {head} True", f"OTS15b/T2heart_organ_ts15 {head} True"]
    sums = (S / "SHA256SUMS").read_text()
    assert sums == (repo / rel / "SHA256SUMS").read_text()
    assert [ln.split()[1] for ln in sums.splitlines()] == [
        *(f"{step}/{name}{h}.h5ad" for step, (name, _, _) in IDENTITY.items() for h in ("", "_half")),
        *(f"OTS15{x}/T2heart_organ_ts15{h}.h5ad" for x in "ab" for h in ("", "_half"))]
    digests = {ln.split()[1]: ln.split()[0] for ln in sums.splitlines()}
    for step, (name, full, half) in IDENTITY.items():
        assert digests[f"{step}/{name}.h5ad"] == full and digests[f"{step}/{name}_half.h5ad"] == half, step
    for h in ("", "_half"):
        assert digests[f"OTS15a/T2heart_organ_ts15{h}.h5ad"] == digests[f"OTS15b/T2heart_organ_ts15{h}.h5ad"]
    calls = _calls(tmp_path)
    py = "run --frozen --no-sync python -I -B experiments/builds/build_test.py heart"
    assert [c for c in calls if "build_test.py" in c] == [
        *(f"{py} {ARGS[step]} --out {rel}/{step}/{name}.h5ad" for step, (name, _, _) in IDENTITY.items()),
        *(f"{py} {ARGS['OTS15']} --out {rel}/OTS15{x}/T2heart_organ_ts15.h5ad" for x in "ab")]
    halves = [c for c in calls if "experiments/builds/halve.py" in c]
    assert len(halves) == 7 and all("--seed" not in c for c in halves)  # halve.py's default seed, 1000001
    assert halves[5].endswith(f"{rel}/OTS15a/T2heart_organ_ts15.h5ad {rel}/OTS15a/T2heart_organ_ts15_half.h5ad")
    checks = [c for c in calls if " vec check " in c]
    assert len(checks) == 14 and all(c.endswith("--board T2:heart:test_extrap") for c in checks)
    tests = [c for c in calls if " -m pytest " in c]
    assert len(tests) == 1
    for f in ("test_heart_trend_line_g.py", "test_heart_trend_line_g_launcher.py", "test_heart_trend_line_f.py",
              "test_heart_trend_line_f_launcher.py", "test_heart_trend_line_d.py",
              "test_heart_trend_line_d_launcher.py", "test_heart_trend_line.py", "test_heart_trend_line_launcher.py",
              "test_build_heart_organ_scope.py",
              "test_heart_organ_launchers.py", "test_h28_heart_trend_strength.py", "test_halve.py",
              "test_hidden_counts.py", "test_s2_rules.py"):
        assert f"tests/{f}" in tests[0], f
    # a run folder is never reused
    assert _launch(tmp_path, head, "r1")[0] == 91


@pytest.mark.parametrize("step", list(IDENTITY))
@pytest.mark.parametrize("half", [False, True])
def test_the_launcher_stops_unless_each_identity_step_gives_its_committed_bytes(tmp_path, repo, step, half):
    _, _, head = repo
    var = f"STUB_SHA_{step}" + ("_HALF" if half else "")
    code = STOP_CODES[step][int(half)]
    assert _launch(tmp_path, head, "s", extra={var: "0" * 64})[0] == code
    assert _exit(tmp_path, "s", step) == str(code)
    later = STEPS[STEPS.index(step) + 1:]
    assert all(_exit(tmp_path, "s", s) is None for s in later)
    assert not any("--trend-strength 1.5" in c for c in _calls(tmp_path))  # OTS15 is never built


def test_the_launcher_stops_at_file_d_on_its_real_hash(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "r2", sha_stub=False)[0] == 20
    assert _exit(tmp_path, "r2", "D") == "20" and _exit(tmp_path, "r2", "OTS05") is None
    assert not any("--trend-strength" in c for c in _calls(tmp_path))


def test_the_launcher_fails_on_a_build_a_check_a_stamp_the_tests_or_two_different_files(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "b1", extra={"STUB_BUILD_EXIT": "1"})[0] == 11
    assert _exit(tmp_path, "b1", "D") == "11"
    assert _launch(tmp_path, head, "b2", extra={"STUB_CHECK_EXIT": "1"})[0] == 13
    assert _launch(tmp_path, head, "b3", extra={"STUB_DIRTY": "1"})[0] == 15
    assert _launch(tmp_path, head, "b4", extra={"STUB_PYTEST_EXIT": "1"})[0] == 10
    assert _exit(tmp_path, "b4", "tests") == "10" and _exit(tmp_path, "b4", "D") is None
    # OTS15a and OTS15b differ (their builds only): done fails, every step before it passed
    assert _launch(tmp_path, head, "b5", extra={"STUB_NONDET": "1", "STUB_NONDET_ON": "--trend-strength 1.5 "})[0] == 31
    assert _exit(tmp_path, "b5", "done") == "31"
    assert all(_exit(tmp_path, "b5", s) == "0" for s in STEPS[:-1])


def test_the_launcher_gate_refusals(tmp_path, repo):
    repo, git, head = repo
    assert _launch(tmp_path, head[:12], "g0")[0] == 90
    assert _launch(tmp_path, head, "a/b")[0] == 90
    assert _launch(tmp_path, head, "it's")[0] == 90
    assert _launch(tmp_path, head, "$(id)")[0] == 90
    # the row is fixed: a third argument (a ROW given by habit, or anything else) is refused before any folder
    for i, extra in enumerate(("R52", "R49", "", "--board T2:heart:val_extrap")):
        assert _launch(tmp_path, head, f"x{i}", extra)[0] == 90, extra
        assert not _scratch(tmp_path, f"x{i}").exists(), extra
    for args in ((), (head,), (head, "")):
        code = _launch(tmp_path, *args)[0]
        assert code != 0 and not (tmp_path / "home" / "v" / "scratch" / "heartg").exists(), args
    assert _launch(tmp_path, "d" * 40, "g1")[0] == 96 and _exit(tmp_path, "g1", "gate") == "96"
    (repo / "src" / "a.py").write_text("x = 2\n")
    assert _launch(tmp_path, head, "g2")[0] == 98
    git("checkout", "-q", "--", "src/a.py")
    (repo / "runs").mkdir()
    (repo / "runs" / f"heartg_{head[:7]}_g3").mkdir()
    assert _launch(tmp_path, head, "g3")[0] == 93
    (repo / "docs" / "S2_RULES.md").write_text("rules, changed\n")
    git("commit", "-q", "-am", "rules changed without a new pin")
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "pin")[0] == 97
    git("revert", "--no-edit", "HEAD")
    t = repo / "experiments" / "harness" / "TRIES.md"
    good = t.read_text()
    bad = {"none": good.replace(R52_ROW, ""), "two": good + R52_ROW, "other": good.replace("| R52 (", "| R53 ("),
           "longer": good.replace("| R52 (", "| R520 ("),
           "ruled": good.replace("`--captain-ruling R52`", "`--captain-ruling R49`"),
           **{f"drop{i}": good.replace(s, "x") for i, s in enumerate((*R52_STRINGS, "ruling, no run",
                                                                       "the captain's order"))},
           "later_cell": good.replace("with `--trend-strength 1.5` only", "only").replace(
               "| heart | x |\n| H34", "| `--trend-strength 1.5` | x |\n| H34"),
           "wider": good.replace("`--trend-strength 1.5`", "`--trend-strength 1.50`"),
           "kind": good.replace("R52 (the captain's order", "R52 (the captain's named exception")}
    for name, text in bad.items():
        assert text != good, name
        t.write_text(text)
        git("commit", "-q", "-am", f"ledger {name}")
        assert _launch(tmp_path, git("rev-parse", "HEAD"), f"row_{name}")[0] == 95, name
        assert _exit(tmp_path, f"row_{name}", "gate") == "95"
        t.write_text(good)
        git("commit", "-q", "-am", "ledger back")
    assert not any("build_test.py" in c or " -m pytest " in c for c in _calls(tmp_path))
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "ok", extra={"STUB_PYTEST_EXIT": "1"})[0] == 10


def _gate_code() -> str:
    sh = LAUNCHER.read_text()
    code = sh[sh.index("python3 -I -B -c '") + len("python3 -I -B -c '"):]
    return code[:code.index("' \\\n")]


def test_the_gates_row_check_is_the_builders(tmp_path, monkeypatch):
    """The gate's inline check and the builder's _heart_ots15_ordered take and refuse the same ledgers."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    code = _gate_code()
    base = "| R29 (other) | heart | x |\n" + R52_ROW + "| H34 (row) | heart | x |\n"
    ledgers = {"ok": base, "only": R52_ROW, "none": "| R27 (x) | T1 | x |\n", "two": base + R52_ROW,
               "other": base.replace("| R52 (", "| R53 ("), "longer": base.replace("| R52 (", "| R520 ("),
               "ruled": base.replace("`--captain-ruling R52`", "`--captain-ruling R49`"),
               "wider": base.replace("`--trend-strength 1.5`", "`--trend-strength 1.50`"),
               "exp": base.replace("`--trend-strength 1.5`", "`--trend-strength 1.5e0`"),
               "more": base.replace("with `--trend-strength 1.5`", "with `--trend-strength 0.75` and "
                                                                   "`--trend-strength 1.5`"),
               "kind": base.replace("R52 (the captain's order", "R52 (the captain's named exception"),
               "unpiped": base.replace("| R52 (", "|R52 ("),
               **{f"drop{i}": base.replace(s, "x") for i, s in enumerate((*R52_STRINGS, "ruling, no run",
                                                                           "the captain's order"))}}
    seen = set()
    for name, text in ledgers.items():
        path = tmp_path / f"{name}.md"
        path.write_text(text)
        monkeypatch.setattr(bt, "TRIES_FILE", path)
        gate = subprocess.run(["python3", "-I", "-B", "-c", code, str(path)]).returncode == 0
        assert gate == bt._heart_ots15_ordered(), name
        seen.add(gate)
    assert seen == {True, False}


def test_the_gate_takes_the_committed_ledger():
    """The gate's inline check passes on experiments/harness/TRIES.md as committed (read by code only), as the
    builder's own does (tests/test_heart_trend_line_g.py)."""
    r = subprocess.run(["python3", "-I", "-B", "-c", _gate_code(), str(REPO_ROOT / "experiments" / "harness" /
                                                                         "TRIES.md")])
    assert r.returncode == 0


# ================================================================================================ as text
def _body(sh: str) -> str:
    return "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))


def test_the_launcher_deletes_nothing_parses_and_keeps_the_f_launchers_conventions():
    sh, f = LAUNCHER.read_text(), F_LAUNCHER.read_text()
    body = _body(sh)
    assert not re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|truncate|shred)\s", body)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    assert len(re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh)) == 2 and "${3" not in sh and "$ROW" not in sh
    assert "python -I -B" in sh and 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'TMPDIR="$S/tmp"' in sh and "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "/mnt/" not in sh and "vec-campaign/" not in body
    assert "S=$HOME/v/scratch/heartg/$RUN" in sh and "O=runs/heartg_${head:0:7}_$RUN" in sh
    # PATH is set before the first uv call, and the argument count is checked before any folder is made
    assert body.index("export PATH=$HOME/.local/bin:$PATH") < body.index("uv run")
    assert body.index('[ "$#" = 2 ]') < body.index('mkdir -p "$HOME')
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    # the shared parts are heart_trend_line_f.sh's, line for line: (start, this launcher's end, the F launcher's end)
    for start, mine_end, their_end in (("finish() {", "step_ok() {", "step_ok() {"),
                                       ("stamp() {", ': > "$S/stamps"', ': > "$S/stamps"'),
                                       ("head=$(git rev-parse HEAD)", "python3 -I -B -c", "python3 -I -B -c"),
                                       ("mkdir -p runs\n", 'echo "$head"', 'echo "$head"'),
                                       ("STEP=D\n", "STEP=OTS075\n", "for STEP in OTS075a")):
        i, j = sh.index(start), f.index(start)
        assert sh[i:sh.index(mine_end, i)] == f[j:f.index(their_end, j)], start
    assert sh.count("build() {") == f.count("build() {") == 1
    # the arguments' block is the F launcher's without ROW, with the argument count added
    mine = sh[sh.index("set -u\n"):sh.index("export PATH=")]
    theirs = f[f.index("set -u\n"):f.index("export PATH=")].replace("heart_trend_line_f.sh COMMIT RUN ROW",
                                                                     "heart_trend_line_g.sh COMMIT RUN")
    kept = [ln for ln in theirs.splitlines() if "ROW" not in ln]
    assert [ln for ln in mine.splitlines() if ln in theirs.splitlines()] == kept
    assert [ln for ln in mine.splitlines() if ln not in theirs.splitlines()] == [
        '[ "$#" = 2 ] || { echo "refused: COMMIT and RUN only (the row, R52, is fixed)" >&2; exit 90; }']
    # the two new identity steps are the OTS10 step's, line for line, with their own names and codes
    ots10 = sh[sh.index("STEP=OTS10\n"):sh.index("STEP=OTS075\n")]
    for step, row, codes in (("OTS075", "R49", (26, 27)), ("OTS125", "R49", (28, 29))):
        name = IDENTITY[step][0]
        want = (ots10.replace("T2heart_organ_ts10", name).replace("OTS10", step).replace("row R29's", f"row {row}'s")
                .replace("finish 24", f"finish {codes[0]}").replace("finish 25", f"finish {codes[1]}"))
        nxt = "STEP=OTS125\n" if step == "OTS075" else "for STEP in OTS15a"
        assert sh[sh.index(f"STEP={step}\n"):sh.index(nxt)] == want, step


def test_the_launcher_reads_only_exit_sha_and_stamp():
    sh = LAUNCHER.read_text()
    for line in (ln for ln in sh.splitlines() if not ln.lstrip().startswith("#")):
        if line.lstrip().startswith(("$PY experiments/builds/", "uv run --frozen --no-sync vec check")):
            assert "> /dev/null 2>&1" in line, line
    body = _body(sh)
    assert 'print(s["commit"], s["clean"])' in body and '"size"' not in body
    assert not re.search(r"(^|[\s;&|(])(stat|du|wc|ls|cat|head|tail|less|h5ls)\s", body)
    for step, args in ARGS.items():
        assert f'{step}_ARGS="{args}"' in sh, step
    assert "BOARD=T2:heart:test_extrap" in sh


def test_the_identity_hashes_are_the_committed_lists():
    """D's hashes are file D's in experiments/builds/results/heartorgan_72b9662_d1.SHA256SUMS,
    heartlined_0d47a67_od1.SHA256SUMS and heartf_ee53e77_f1.SHA256SUMS; OTS05's and OTS10's are row R29's in the last
    two (every build of them); OTS075's and OTS125's are row R49's in heartf_ee53e77_f1.SHA256SUMS (both of its builds);
    the launcher holds exactly these."""
    sh = LAUNCHER.read_text()

    def listed(path):
        out: dict[str, set] = {}
        for ln in path.read_text().splitlines():
            digest, name = ln.split()
            out.setdefault(name.split("/", 1)[1], set()).add(digest)
        return out

    d, od, f1 = listed(D_LIST), listed(OD_LIST), listed(F_LIST)
    for step, (name, full, half) in IDENTITY.items():
        lists = (d, od, f1) if step == "D" else (od, f1) if step in ("OTS05", "OTS10") else (f1,)
        for lst in lists:
            assert lst[f"{name}.h5ad"] == {full} and lst[f"{name}_half.h5ad"] == {half}, step
        assert re.search(rf"^{step}_SHA={full}\b", sh, re.M), step
        assert re.search(rf"^{step}_HALF_SHA={half}\b", sh, re.M), step
    assert sum(1 for ln in F_LIST.read_text().splitlines() if "ts075" in ln or "ts125" in ln) == 8  # a and b


def test_the_launchers_commands_are_the_builders():
    """The commands the launcher builds are the builder's own: HEART_ORGAN_SCOPE, row R29's HEART_OTS_COMMANDS, row
    R49's HEART_FURTHER_COMMANDS with the pinned row and row R52's HEART_OTS15_COMMAND."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    sh = LAUNCHER.read_text()
    assert f'D_ARGS="{bt.HEART_ORGAN_SCOPE[len("heart "):]}"' in sh
    for s, var in ((0.5, "OTS05_ARGS"), (1.0, "OTS10_ARGS")):
        assert f'{var}="{bt.HEART_OTS_COMMANDS[s][len("heart "):]}"' in sh
    for s, var in ((0.75, "OTS075_ARGS"), (1.25, "OTS125_ARGS")):
        assert f'{var}="{bt.HEART_FURTHER_COMMANDS[s].format(row=bt.HEART_FURTHER_ROW)[len("heart "):]}"' in sh
    assert f'OTS15_ARGS="{bt.HEART_OTS15_COMMAND[len("heart "):]}"' in sh
    assert bt.HEART_FURTHER_ROW == "R49" and bt.HEART_OTS15_ROW == "R52"
    # the gate's strings are the builder's
    code = _gate_code()
    for s in (bt.HEART_OTS_ENTRY, bt.HEART_OTS15_OPTION, f"`--captain-ruling {bt.HEART_OTS15_ROW}`"):
        assert s in code, s
