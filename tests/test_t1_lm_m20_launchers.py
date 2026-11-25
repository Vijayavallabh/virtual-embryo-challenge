"""The four launchers of rows N64, T72, T73, R54 and R58 (the step-change round's T1 lines):
  experiments/t1_scope/n64_step0.sh       row N64's step 0 (late_mix.py)
  experiments/t1_scope/t72_late_mix.sh    row T72's read (t72_late_mix.py)
  experiments/t1_m1/t73_m20.sh            row T73's read (m20_read.py), gate G0 of row R58's file
  experiments/builds/t1_lm_m20_line.sh    the test-board builds of rows R54's and R58's files
These tests parse each with bash, hold their refusal codes, PATH and thread settings to their headers, check that none
deletes anything or reads a test build's log, size or shape, and run every gate path for real in throwaway repositories
under a temporary home, with a stand-in `uv` in ~/.local/bin (each launcher puts that folder first on PATH): one that
fails, so a call passing the gate stops at its tests step, and, for the build launcher, one that stands in for every
build so a whole run goes through every step in order. Nothing here reads the repository's own ledger or rules file."""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
from pathlib import Path

import pytest
from test_t1_lm_m20_build import ROW_R54, ROW_R58

REPO = Path(__file__).resolve().parents[1]
N64 = REPO / "experiments" / "t1_scope" / "n64_step0.sh"
T72 = REPO / "experiments" / "t1_scope" / "t72_late_mix.sh"
T73 = REPO / "experiments" / "t1_m1" / "t73_m20.sh"
LINE = REPO / "experiments" / "builds" / "t1_lm_m20_line.sh"
LAUNCHERS = (N64, T72, T73, LINE)
SHARES_PATH = "experiments/t1_scope/results/n64_late_mix/shares.json"
SHARES = '{"row": "N64", "late_mix": "e15_17"}\n'
SHA = hashlib.sha256(SHARES.encode()).hexdigest()
RULES = "the amendment names the rows\n"
FAILING_UV = "#!/bin/sh\nexit 1\n"
PASSING_UV = "#!/bin/sh\necho \"$*\" >> \"$HOME/uv_calls.txt\"\nexit 0\n"


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _cells(*cells: str) -> str:
    return "| " + " | ".join(cells) + " |"


def _row(first: str, last: str = "pending") -> str:
    """A ledger row of eleven cells."""
    return _cells(first, "T1", *[f"c{i}" for i in range(3, 11)], last)


N64_ROW = _row("N64 (step 0 of row R54's late-mix line, target-free, not a try: `experiments/t1_scope/late_mix.py`, "
               "launcher `experiments/t1_scope/n64_step0.sh`)", f"Run: `shares.json sha256 {SHA}`")
T72_FIRST = ("T72 (labelled diagnostic, the practice read that binds row R54's file: "
             "`experiments/t1_scope/t72_late_mix.py`, launcher `experiments/t1_scope/t72_late_mix.sh`, at row N64's "
             f"shares (`shares.json sha256 {SHA}`))")
T72_ROW = _row(T72_FIRST, "Run; in the read's words: \"Outcome: T72 PASS\" ...")
T73_FIRST = ("T73 (gate G0 of row R58's file: `experiments/t1_m1/m20_read.py`, launcher "
             "`experiments/t1_m1/t73_m20.sh`; fifteen tries)")
T73_ROW = _row(T73_FIRST, "Run; in the read's words: \"Outcome: T73 PASS (m20's step is AHEAD)\"")
REGISTER = ("| X001 | x |\n| X002 | Feng | n/a | data | stages | removed | T1 | No absolute external levels. "
            "Amendment of 2026-11-25 (rows N64, T72 and R54): the counts enter. | added |\n")


def _env(tmp_path: Path, uv: str = FAILING_UV) -> dict:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    env["UV_PROJECT_ENVIRONMENT"] = str(tmp_path / "elsewhere")
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text(uv)
    (fake / "uv").chmod(0o755)
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: list[str], files: dict[str, str] | None = None, ignored: tuple[str, ...] = (),
          venv: str | None = "own", uv: str = FAILING_UV, good_rules: bool = True) -> tuple[Path, dict, str]:
    """A throwaway repository holding the four launchers, a ledger of `tries`, the rules file and its hash, `files`
    (tracked) and `ignored` (untracked, ignored), a .venv whose vec.pth names its src (or another's, or none), and a
    home whose ~/.local/bin holds the stand-in `uv`. Returns the repository, the env and its HEAD."""
    env = _env(tmp_path, uv)
    repo = tmp_path / "repo"
    for path in LAUNCHERS:
        dst = repo / path.relative_to(REPO)
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(path.read_text())
    (repo / "experiments" / "harness").mkdir(parents=True, exist_ok=True)
    (repo / "experiments" / "harness" / "TRIES.md").write_text("# ledger\n" + "".join(f"{t}\n" for t in tries))
    (repo / "docs").mkdir(exist_ok=True)
    (repo / "docs" / "S2_RULES.md").write_text(RULES)
    sha = hashlib.sha256(RULES.encode()).hexdigest() if good_rules else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    (repo / "src").mkdir(exist_ok=True)
    (repo / "src" / "keep.txt").write_text("src\n")
    (repo / ".gitignore").write_text("data/\nruns/\n.venv/\n")
    for rel, text in (files or {}).items():
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_text(text)
    for rel in ignored:
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_text("x")
    if venv is not None:
        site = repo / ".venv" / "lib" / "python3.12" / "site-packages"
        site.mkdir(parents=True)
        (site / "vec.pth").write_text(f"{repo / 'src' if venv == 'own' else tmp_path / 'other' / 'src'}\n")
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, launcher: Path, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(launcher.relative_to(REPO)), *argv], cwd=repo, env=env, capture_output=True,
                          text=True)


# ------------------------------------------------------------------ what every launcher holds


@pytest.mark.parametrize("launcher", LAUNCHERS)
def test_bash_parses_it_and_it_puts_the_local_bin_first_before_any_uv(launcher):
    out = subprocess.run(["bash", "-n", str(launcher)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    sh = launcher.read_text()
    code = _code(sh)  # the header's usage lines name `uv run` too
    assert "export PATH=$HOME/.local/bin:$PATH" in code
    assert code.index("export PATH=$HOME/.local/bin:$PATH") < code.index("uv run")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert "-p no:cacheprovider" in sh and "--basetemp=" in sh


@pytest.mark.parametrize("launcher", LAUNCHERS)
def test_it_deletes_nothing_uploads_nothing_and_fetches_nothing(launcher):
    code = _code(launcher.read_text())
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code.splitlines())
    assert not re.search(r"find .*-delete|git (clean|reset|checkout|stash|commit|push|fetch|pull)", code)
    for word in ("upload", "curl", "wget", "scp ", "ssh "):
        assert word not in code, word


@pytest.mark.parametrize("launcher", LAUNCHERS)
def test_the_documented_codes_are_the_ones_used(launcher):
    sh = launcher.read_text()
    head = sh[:sh.index("set -u")]
    code = _code(sh)
    used = {int(c) for c in re.findall(r"(?:finish|exit|return|refuse) (\d+)", code)}
    used |= {int(c) for c in re.findall(r'same [A-Za-z0-9]+ [A-Za-z0-9_]+ "\$[A-Z0-9_]+" (\d+)', code)}
    used -= {0, 1}  # the Python steps' own codes, documented as such
    documented = {int(c) for c in re.findall(r"(?<![\d.])(\d{1,2})(?=[ ),;:])", head)}
    assert used <= documented, (sorted(used - documented), launcher.name)


def test_the_build_launcher_reads_no_build_output_but_codes_hashes_and_stamps():
    sh = LINE.read_text()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for line in lines:
        if re.search(r'build_test\.py t1|t1_lm_m20_check\.py "\$mode"|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    code = _code(sh)
    for word in (".build.json", "du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -",
                 "tail ", "cat "):
        assert word not in code, word
    assert 'print(s["commit"], s["clean"])' in sh and "OMP_NUM_THREADS=16" in sh
    assert "unset UV_PROJECT_ENVIRONMENT" in sh
    base = "--g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult 1.5 --direction-momentum all"
    assert f'LM16_ARGS="{base} --late-mix e15_17 --captain-ruling R54"' in sh
    assert f'M20_ARGS="{base} --direction-momentum-damping 2.0 --captain-ruling R58"' in sh
    assert "MSB_SHA=78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2" in sh


def test_the_build_launchers_identity_hashes_are_the_r51_launchers():
    r51 = (REPO / "experiments" / "builds" / "t1_r51_line.sh").read_text()
    sh = LINE.read_text()
    for var in ("R13_SHA", "B_SHA", "MS_SHA", "MSB_SHA"):
        mine = re.search(rf"^{var}=([0-9a-f]{{64}})", sh, re.M).group(1)
        assert f"{var}={mine}" in r51, var
    for var in ("BASE_ARGS", "B_ARGS", "MS_ARGS", "MSB_ARGS"):
        mine = re.search(rf'^{var}="([^"]+)"', sh, re.M).group(1)
        assert f'{var}="{mine}"' in r51, var


# ------------------------------------------------------------------ the build launcher, run for real


def _line_repo(tmp_path, tries=None, shares=SHARES, track_shares=True, relabel=True, **kw):
    tries = [ROW_R54, ROW_R58, N64_ROW, T72_ROW, T73_ROW] if tries is None else tries
    files = {SHARES_PATH: shares} if track_shares else {}
    ignored = ("data/derived/feng_relabel_c/x",) if relabel else ()
    repo, env, head = _repo(tmp_path, tries, files, ignored, **kw)
    if not track_shares:
        (repo / SHARES_PATH).parent.mkdir(parents=True, exist_ok=True)
        (repo / SHARES_PATH).write_text(shares)
    return repo, env, head


def _scratch(tmp_path, run):
    return tmp_path / "home" / "v" / "scratch" / "t1lmm20" / run


def test_a_call_that_passes_the_gate_stops_at_the_failing_tests_step_and_never_reuses_a_run(tmp_path):
    repo, env, head = _line_repo(tmp_path)
    out = _run(repo, env, LINE, head, "r1", "R54", "R58")
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head and (s / "rows").read_text().strip() == "R54 R58"
    assert not list((repo / "runs" / f"t1lmm20_{head[:7]}_r1").iterdir())
    again = _run(repo, env, LINE, head, "r1", "R58")
    assert again.returncode == 91 and "never reused" in again.stderr
    for rows in (("R54",), ("R58",)):
        name = "r" + "".join(rows)
        assert _run(repo, env, LINE, head, name, *rows).returncode == 10


@pytest.mark.parametrize("argv, code", [
    ([], None), (["r1"], None), (["a/b", "R54"], 90), (["r 1", "R54"], 90),
    *[(["r1", bad], 89) for bad in ("R51", "R26", "r54", "R540", "R5", "54", "R54 ", "$(id)")],
    (["r1", "R54", "R54"], 89), (["r1", "R58", "R54", "R58"], 89),
])
def test_odd_arguments_are_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _line_repo(tmp_path)
    out = _run(repo, env, LINE, head, *argv)
    assert out.returncode != 0 if code is None else out.returncode == code, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


def _replace(rows, old, new):
    return [r.replace(old, new) for r in rows]


BASE_ROWS = [ROW_R54, ROW_R58, N64_ROW, T72_ROW, T73_ROW]


@pytest.mark.parametrize("kw, rows, code, says", [
    ({"good_rules": False}, ("R54",), 97, None),
    ({"tries": [r for r in BASE_ROWS if r is not ROW_R54]}, ("R54",), 95, "exactly one row R54"),
    ({"tries": [*BASE_ROWS, ROW_R58]}, ("R58",), 95, "exactly one row R58"),
    ({"tries": _replace(BASE_ROWS, "R54 (the captain's", "R54 (a draft of the captain's")}, ("R54",), 95,
     "not the captain's ruling"),
    ({"tries": _replace(BASE_ROWS, "`--late-mix e15_17`", "--late-mix e15_17")}, ("R54",), 95,
     "does not name `--late-mix e15_17`"),
    ({"tries": _replace(BASE_ROWS, "damping 2.0`", "damping 2.0 x`")}, ("R58",), 95, "does not name"),
    ({"tries": _replace(BASE_ROWS, "`--captain-ruling R58`", "`--captain-ruling R51`")}, ("R58",), 95,
     "does not name `--captain-ruling R58`"),
    ({"tries": [r for r in BASE_ROWS if r is not T72_ROW]}, ("R54",), 87, "exactly one row T72"),
    ({"tries": _replace(BASE_ROWS, "Outcome: T72 PASS", "Outcome: T72 STOP (x)")}, ("R54",), 87, "does not say"),
    ({"tries": _replace(BASE_ROWS, "Outcome: T73 PASS", "Outcome: T73 VOID (x)")}, ("R58",), 87, "does not say"),
    ({"tries": _replace(BASE_ROWS, "\"Outcome: T72 PASS\" ...", "\"Outcome: T72 PASS\", once \"Outcome: T72 STOP\"")},
     ("R54",), 87, "also says STOP or VOID"),
    ({"tries": _replace(BASE_ROWS, "Outcome: T72 PASS", "Outcome: T72 PASSES")}, ("R54",), 87, "does not say"),
    ({"tries": _replace(BASE_ROWS, "Outcome: T73 PASS", "Outcome: T73 PASS_")}, ("R58",), 87, "does not say"),
    ({"tries": _replace(BASE_ROWS, "| c5 | c6 |", "| c5 |")}, ("R54", "R58"), 87, "eleven cells"),
    ({"tries": _replace(BASE_ROWS, "\"Outcome: T72 PASS\" ...", "PASS")}, ("R54",), 87, "does not say"),
    ({"tries": [r for r in BASE_ROWS if r is not N64_ROW]}, ("R54",), 86, "exactly one row N64"),
    ({"tries": _replace(BASE_ROWS, f"Run: `shares.json sha256 {SHA}`", "Run: written")}, ("R54",), 86,
     "does not name exactly one"),
    ({"tries": _replace(BASE_ROWS, f"Run: `shares.json sha256 {SHA}`",
                        f"Run: `shares.json sha256 {SHA}` `shares.json sha256 {'f' * 64}`")}, ("R54",), 86,
     "does not name exactly one"),
    ({"tries": _replace(BASE_ROWS, f"(`shares.json sha256 {SHA}`))", "())")}, ("R54",), 86,
     "row T72's first cell does not name"),
    ({"track_shares": False}, ("R54",), 86, "is not tracked"),
    ({"shares": SHARES + " "}, ("R54",), 86, "is not the file row N64 names"),
    ({"relabel": False}, ("R58",), 92, "feng_relabel_c is not in this checkout"),
    ({"venv": None}, ("R58",), 88, "no .venv of its own"),
    ({"venv": "other"}, ("R54",), 88, "no .venv of its own"),
])
def test_each_gate_refusal(tmp_path, kw, rows, code, says):
    repo, env, head = _line_repo(tmp_path, **kw)
    out = _run(repo, env, LINE, head, "r1", *rows)
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]
    assert not (repo / "runs" / f"t1lmm20_{head[:7]}_r1").exists()


def test_a_wrong_commit_a_prefix_a_dirty_tree_and_an_existing_output_are_refused(tmp_path):
    repo, env, head = _line_repo(tmp_path)
    assert _run(repo, env, LINE, "0" * 40, "r1", "R54").returncode == 96
    assert _run(repo, env, LINE, head[:12], "r2", "R54").returncode == 96
    o = repo / "runs" / f"t1lmm20_{head[:7]}_r3"
    o.mkdir(parents=True)
    out = _run(repo, env, LINE, head, "r3", "R58")
    assert out.returncode == 93 and "output folder" in out.stderr
    (repo / "experiments" / "harness" / "TRIES.md").write_text("# ledger\nchanged\n")
    out = _run(repo, env, LINE, head, "r4", "R54")
    assert out.returncode == 98 and "tracked files differ" in out.stderr


FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/build_test.py|experiments/builds/t1_lm_m20_check.py)
        prog=$1; shift
        [ "$prog" = experiments/builds/t1_lm_m20_check.py ] && { echo "$1 $2" >> calls.txt; shift; }
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        [ "$prog" = experiments/builds/t1_lm_m20_check.py ] && name=T1_s2rcMSB
        echo "file $name" > "$out"
        printf '{"commit": "%s", "clean": true}' "$(git rev-parse HEAD)" > "$out.stamp.json"
        exit 0 ;;
    esac ;;
esac
exit 1
"""


@pytest.mark.parametrize("rows, folders", [
    (("R54", "R58"), ["D", "B", "MS", "MSB", "LM16own", "m20msb", "LM16a", "LM16b", "m20a", "m20b"]),
    (("R58",), ["D", "B", "MS", "MSB", "m20msb", "m20a", "m20b"]),
    (("R54",), ["D", "B", "MS", "MSB", "LM16own", "LM16a", "LM16b"]),
])
def test_a_whole_stand_in_run_passes_every_step_in_order(tmp_path, rows, folders):
    """The launcher with its four identity hashes set to the stand-in files' (in a copy) and the stand-in `uv` above:
    every step exits 0, the stamps list every build, the checker gets its modes, the hash list holds every built file
    once and each pair is compared."""
    names = {"R13_SHA": "T1_s2rc", "B_SHA": "T1_s2rcB", "MS_SHA": "T1_s2rc_MS", "MSB_SHA": "T1_s2rcMSB"}
    repo, env, head = _line_repo(tmp_path, uv=FULL_UV)
    sh = LINE.read_text()
    for var, name in names.items():
        old = re.search(rf"^{var}=([0-9a-f]{{64}})", sh, re.M).group(1)
        sh = sh.replace(f"{var}={old}", f"{var}={hashlib.sha256(f'file {name}'.encode() + b'\n').hexdigest()}")
    (repo / LINE.relative_to(REPO)).write_text(sh)
    _git(repo, env, "commit", "-q", "-am", "stand-in hashes")
    head = _git(repo, env, "rev-parse", "HEAD")
    out = _run(repo, env, LINE, head, "r1", *rows)
    assert out.returncode == 0, out.stderr
    s = _scratch(tmp_path, "r1")
    codes = {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")}
    assert codes == {step: "0" for step in ["gate", "tests", *folders, "done"]}
    stamps = (s / "stamps").read_text().splitlines()
    assert [line.split()[0].split("/")[0] for line in stamps] == folders
    assert all(line.endswith(f" {head} True") for line in stamps)
    listed = (repo / "runs" / f"t1lmm20_{head[:7]}_r1" / "SHA256SUMS").read_text().splitlines()
    assert [line.split()[1].split("/")[0] for line in listed] == folders
    calls = (repo / "calls.txt").read_text().splitlines()
    assert calls == [f"{m} t1" for m, r in (("lm16own", "R54"), ("m20msb", "R58")) if r in rows]


def test_a_stand_in_run_stops_at_the_first_unknown_bytes(tmp_path):
    repo, env, head = _line_repo(tmp_path, uv=FULL_UV)
    out = _run(repo, env, LINE, head, "r1", "R58")
    assert out.returncode == 20 and "step D does not give the known bytes of T1_s2rc" in out.stderr
    codes = {p.name: p.read_text().strip() for p in _scratch(tmp_path, "r1").glob("*.exit")}
    assert codes == {"gate.exit": "0", "tests.exit": "0", "D.exit": "20"}


# ------------------------------------------------------------------ the read launchers, run for real


def _read_scratch(tmp_path, row):
    return tmp_path / "home" / "v" / "scratch" / row


N64_FILES = {"EXTERNAL_DATA.md": REGISTER, "experiments/external/results/feng_relabel_late/MANIFEST.json": "{}",
             "experiments/external/results/feng_relabel/MANIFEST.json": "{}"}
T72_FILES = {SHARES_PATH: SHARES, "experiments/t1_scope/results/t30_pricing/scenarios/t30_d1.0_cmown.json": "{}"}


def _n64_repo(tmp_path, tries=None, files=None, data=True, **kw):
    tries = [N64_ROW] if tries is None else tries
    ignored = ("data/raw/T1/E10.5_RNA.h5ad",) if data else ()
    return _repo(tmp_path, tries, N64_FILES if files is None else files, ignored, **kw)


def test_step0_runs_its_two_steps_in_order(tmp_path):
    repo, env, head = _n64_repo(tmp_path, uv=PASSING_UV)
    s = _read_scratch(tmp_path, "n64")
    assert _run(repo, env, N64, "step0", head).returncode == 94  # the tests first
    assert (s / "step0.refused").read_text().strip() == "94" and not (s / "step0.exit").exists()
    assert _run(repo, env, N64, "tests", head).returncode == 0
    assert (s / "tests.exit").read_text().strip() == "0" and (s / "tests.head").read_text().strip() == head
    assert _run(repo, env, N64, "step0", head).returncode == 0
    calls = (tmp_path / "home" / "uv_calls.txt").read_text().splitlines()
    assert calls[-1] == f"run --frozen --no-sync python -I -B experiments/t1_scope/late_mix.py --commit {head}"
    assert "tests/test_t1_late_mix.py" in calls[0]


@pytest.mark.parametrize("case, code, says", [
    ("usage", 2, "usage"), ("unknown", 97, "unknown step"), ("dirty", 98, "not clean"),
    ("no_row", 96, "exactly one row N64"), ("two_rows", 96, "exactly one row N64"),
    ("not_first", 96, "does not start 'N64 ('"), ("no_launcher", 96, "n64_step0"),
    ("longer_path", 96, "late_mix"), ("short", 95, ""), ("other_head", 95, "HEAD is not"),
    ("venv", 88, "no .venv of its own"), ("no_data", 93, "E10.5_RNA.h5ad is not in place"),
    ("untracked", 93, "is not tracked"), ("no_amendment", 92, "does not hold the amendment"),
    ("output", 91, "exists"),
])
def test_each_step0_refusal(tmp_path, case, code, says):
    tries, files, data, venv = None, None, True, "own"
    if case == "no_row":
        tries = []
    elif case == "two_rows":
        tries = [N64_ROW, N64_ROW]
    elif case == "not_first":
        tries = [N64_ROW.replace("| N64 (", "| N64 step (")]
    elif case == "no_launcher":
        tries = [N64_ROW.replace("`experiments/t1_scope/n64_step0.sh`", "the launcher")]
    elif case == "longer_path":
        tries = [N64_ROW.replace("late_mix.py`", "late_mix.pyc`")]
    elif case == "venv":
        venv = "other"
    elif case == "no_data":
        data = False
    elif case == "untracked":
        files = {k: v for k, v in N64_FILES.items() if "feng_relabel/" not in k}
    elif case == "no_amendment":
        files = {**N64_FILES, "EXTERNAL_DATA.md": REGISTER.replace("Amendment of 2026-11-25", "Amendment")}
    repo, env, head = _n64_repo(tmp_path, tries, files, data, venv=venv, uv=PASSING_UV)
    if case == "untracked":
        p = repo / "experiments/external/results/feng_relabel/MANIFEST.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        (repo / ".gitignore").write_text("data/\nruns/\n.venv/\nexperiments/external/results/feng_relabel/\n")
        p.write_text("{}")
        _git(repo, env, "commit", "-q", "-am", "ignore")
        head = _git(repo, env, "rev-parse", "HEAD")
    s = _read_scratch(tmp_path, "n64")
    if case in ("no_data", "untracked", "no_amendment", "output"):
        assert _run(repo, env, N64, "tests", head).returncode == 0
    if case == "output":
        out = repo / "experiments/t1_scope/results/n64_late_mix"
        out.mkdir(parents=True)
        (out / "n64.txt").write_text("x")
    argv = {"usage": [], "unknown": ["../x", head], "dirty": ["tests", head], "short": ["tests", head[:12]],
            "other_head": ["tests", "d" * 40]}.get(case, ["step0" if code in (93, 92, 91) else "tests", head])
    if case == "dirty":
        (repo / "experiments" / "harness" / "TRIES.md").write_text("changed\n")
    out = _run(repo, env, N64, *argv)
    assert out.returncode == code, out.stderr
    assert says in out.stderr
    name = {"usage": "usage", "unknown": "unknown"}.get(case, argv[0] if argv else "usage")
    assert (s / f"{name}.refused").read_text().strip() == str(code)
    assert not (s / "x.refused").exists() and not (tmp_path / "home" / "v" / "scratch" / "x.refused").exists()


def _t72_repo(tmp_path, tries=None, files=None, data=True, **kw):
    tries = [N64_ROW, T72_ROW] if tries is None else tries
    ignored = ("data/raw/T1/E10.5_RNA.h5ad",) if data else ()
    return _repo(tmp_path, tries, T72_FILES if files is None else files, ignored, **kw)


def test_t72_runs_its_steps_in_order_and_starts_the_twelve_scenarios_together(tmp_path):
    repo, env, head = _t72_repo(tmp_path, uv=PASSING_UV)
    s = _read_scratch(tmp_path, "t72")
    assert _run(repo, env, T72, "runs", head).returncode == 94
    assert _run(repo, env, T72, "tests", head).returncode == 0
    assert _run(repo, env, T72, "read", head).returncode == 94  # the runs first
    assert _run(repo, env, T72, "runs", head).returncode == 0
    assert sorted(p.name for p in (s / "runs").glob("*.exit")) == sorted(
        f"{cm}_{b}.exit" for cm in ("own", "0.30", "0.33", "lm16", "0.45", "0.60") for b in ("own", "lm16"))
    calls = (tmp_path / "home" / "uv_calls.txt").read_text().splitlines()
    runs = sorted(c for c in calls if " run --cm " in c)
    assert len(runs) == 12 and all(c.endswith(f"--commit {head} --jobs 2") for c in runs)
    assert any("run --cm lm16 --blood own " in c for c in runs)
    assert any("run --cm 0.30 --blood lm16 " in c for c in runs)
    assert _run(repo, env, T72, "runs", head).returncode == 91  # never twice
    assert _run(repo, env, T72, "read", head).returncode == 0
    assert calls == (tmp_path / "home" / "uv_calls.txt").read_text().splitlines()[:len(calls)]
    assert (tmp_path / "home" / "uv_calls.txt").read_text().splitlines()[-1].endswith(
        f"experiments/t1_scope/t72_late_mix.py read --commit {head}")


@pytest.mark.parametrize("case, code, says", [
    ("usage", 2, "usage"), ("unknown", 97, "unknown step"), ("no_row", 96, "exactly one row T72"),
    ("not_first", 96, "does not start 'T72 ('"), ("no_launcher", 96, "t72_late_mix"),
    ("no_sha", 96, "do not each name the one same"), ("other_sha", 96, "do not each name the one same"),
    ("no_n64", 96, "do not each name the one same"), ("short", 95, ""), ("other_head", 95, "HEAD is not"),
    ("venv", 88, "no .venv"), ("no_data", 93, "is not in place"), ("bad_shares", 93, "is not the file rows"),
    ("no_t30", 93, "is not tracked"), ("dirty", 98, "not clean"),
])
def test_each_t72_refusal(tmp_path, case, code, says):
    tries, files, data, venv = None, None, True, "own"
    if case == "no_row":
        tries = [N64_ROW]
    elif case == "not_first":
        tries = [N64_ROW, T72_ROW.replace("| T72 (", "| T72 read (")]
    elif case == "no_launcher":
        tries = [N64_ROW, T72_ROW.replace("`experiments/t1_scope/t72_late_mix.sh`", "its launcher")]
    elif case == "no_sha":
        tries = [N64_ROW, T72_ROW.replace(f"(`shares.json sha256 {SHA}`)", "")]
    elif case == "other_sha":
        tries = [N64_ROW, T72_ROW.replace(SHA, "f" * 64)]
    elif case == "no_n64":
        tries = [T72_ROW]
    elif case == "venv":
        venv = None
    elif case == "no_data":
        data = False
    elif case == "bad_shares":
        files = {**T72_FILES, SHARES_PATH: SHARES + "x"}
    elif case == "no_t30":
        files = {SHARES_PATH: SHARES}
    repo, env, head = _t72_repo(tmp_path, tries, files, data, venv=venv, uv=PASSING_UV)
    if code == 93:
        assert _run(repo, env, T72, "tests", head).returncode == 0
    if case == "dirty":
        (repo / SHARES_PATH).write_text("changed\n")
    argv = {"usage": ["runs"], "unknown": ["x", head], "short": ["tests", head[:12]],
            "other_head": ["tests", "d" * 40]}.get(case, ["runs" if code == 93 else "tests", head])
    out = _run(repo, env, T72, *argv)
    assert out.returncode == code, out.stderr
    assert says in out.stderr
    name = {"usage": "usage", "unknown": "unknown"}.get(case, argv[0])
    assert (_read_scratch(tmp_path, "t72") / f"{name}.refused").read_text().strip() == str(code)


def _t73_repo(tmp_path, tries=None, relabel=True, **kw):
    tries = [T73_ROW] if tries is None else tries
    ignored = ("data/derived/feng_relabel_c/cd1_stats_relabel.npz",) if relabel else ()
    return _repo(tmp_path, tries, {}, ignored, **kw)


def test_t73_runs_its_steps_in_order(tmp_path):
    repo, env, head = _t73_repo(tmp_path, uv=PASSING_UV)
    assert _run(repo, env, T73, "cut", "chamber", head).returncode == 94
    assert _run(repo, env, T73, "tests", head).returncode == 0
    for step in (("cut", "chamber"), ("cut", "whole"), ("lib", "1"), ("lib", "2")):
        assert _run(repo, env, T73, *step, head).returncode == 0
    assert _run(repo, env, T73, "read", head).returncode == 94  # library 3 has not run
    assert _run(repo, env, T73, "lib", "3", head).returncode == 0
    assert _run(repo, env, T73, "read", head).returncode == 0
    calls = (tmp_path / "home" / "uv_calls.txt").read_text().splitlines()
    assert calls[1].endswith(f"experiments/t1_m1/m20_read.py run-cut --cut chamber --out t73_m20 --commit {head}")
    assert calls[-1].endswith(f"experiments/t1_m1/m20_read.py read --out t73_m20 --commit {head}")
    task = repo / "experiments/t1_m1/results/t73_m20/tasks/lib1.json"
    task.parent.mkdir(parents=True)
    task.write_text("{}")
    assert _run(repo, env, T73, "lib", "1", head).returncode == 91  # never overwritten


@pytest.mark.parametrize("case, code, says", [
    ("usage", 2, "usage"), ("usage_cut", 2, "usage"), ("unknown", 97, "unknown step"),
    ("unknown_lib", 97, "unknown step"), ("no_row", 96, "exactly one row T73"),
    ("not_first", 96, "does not start 'T73 ('"), ("no_launcher", 96, "t73_m20"), ("short", 95, ""),
    ("other_head", 95, "HEAD is not"), ("venv", 88, "no .venv"), ("no_relabel", 93, "relabelling's output"),
    ("dirty", 98, "not clean"),
])
def test_each_t73_refusal(tmp_path, case, code, says):
    tries, relabel, venv = None, True, "own"
    if case == "no_row":
        tries = []
    elif case == "not_first":
        tries = [T73_ROW.replace("| T73 (", "| T73 gate (")]
    elif case == "no_launcher":
        tries = [T73_ROW.replace("`experiments/t1_m1/t73_m20.sh`", "its launcher")]
    elif case == "venv":
        venv = "other"
    elif case == "no_relabel":
        relabel = False
    repo, env, head = _t73_repo(tmp_path, tries, relabel, venv=venv, uv=PASSING_UV)
    if case == "dirty":
        (repo / "docs" / "S2_RULES.md").write_text("changed\n")
    argv = {"usage": ["tests"], "usage_cut": ["cut", head], "unknown": ["x", head], "unknown_lib": ["lib", "4", head],
            "short": ["tests", head[:12]], "other_head": ["tests", "d" * 40]}.get(case, ["tests", head])
    out = _run(repo, env, T73, *argv)
    assert out.returncode == code, out.stderr
    assert says in out.stderr
    name = {"usage": "usage", "usage_cut": "usage", "unknown": "unknown", "unknown_lib": "unknown"}.get(case, "tests")
    assert (_read_scratch(tmp_path, "t73") / f"{name}.refused").read_text().strip() == str(code)
