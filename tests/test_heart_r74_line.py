"""The test-board build launcher of row R74's four files (OTK125's command with the category form's strength from the
builder's HEART_R74_TABLE), one file per call: experiments/builds/heart_r74_line.sh <COMMIT> <RUN> <FILE_NAME>.

It refuses any name but the four of its literal table, odd arguments, a reused run folder and an odd commit before it
makes a folder (the near-miss names OTK125, OTK1, OTK100, OTK0.75 and otk10 among them); its gate checks the commit,
the clean tree, the rules file's hash, row R74's first cell (`--entry organ_scope`, `--z-mode raw`, the file's option
string and `--captain-ruling R74` as whole backticked strings, the option string paired there with the file's name),
the five released stages and the venv; then it runs the targeted tests, rebuilds OTK125 and its half (halve.py --seed
1000001) and stops unless both give their recorded sha256; then it builds the file twice from one clean commit, halves
each at seed 1000001, checks each half, and compares the two builds and refuses a file with OTK125's bytes. These tests
hold its table to the builder's, its hashes to the committed list, its commands to the ones the builder takes, check
that it deletes nothing and reads no output of a test build but exit codes, sha256 and the stamps' commit and clean
fields, and run it for real in throwaway repositories under a temporary home: with a stand-in `uv` that fails (every
refusal path, and a call passing the gate stops at the tests step and builds nothing) and with stand-ins that pass (a
whole run of each of the four names goes through every step in order; a stand-in of other bytes, a failing build,
halve or check, an unclean stamp and a file with OTK125's bytes stop the run where their step says). Nothing here reads
the repository's own ledger or rules file; build_test is imported only to check its table and that it takes the commands
the launcher runs, against a temporary ledger."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import GIT_LOCATION_VARS
from test_heart_r74_build import NAMES, R74_FIRST, TABLE
from test_heart_trend_line_g import LEDGER
from test_heart_trend_line_h import R55_ROW

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "heart_r74_line.sh"
WBH = REPO / "experiments" / "builds" / "heart_wbh_build.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
OTK125_SHA = "d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893"
OTK125_HALF_SHA = "9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54"
OTK125_ARGS = "--entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55"
STAGES = ("E10.5", "E9.5", "E8.75", "E8.5", "E8.25_late")
ROW = f"| {R74_FIRST} | heart | the refinement round's four files |\n"
FIRST = "OTK15"  # the file the gate's refusals are run on
FIRST_OPT = TABLE[FIRST][0]
RULES = "the second-upload rules\n"


def _sh() -> str:
    return LAUNCHER.read_text()


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _table(sh: str) -> dict[str, tuple[str, str]]:
    """The launcher's FILE_NAME -> (option string, file stem) table (its case lines)."""
    rows = re.findall(r'^  (OTK[0-9]+)\) OPT="([^"]+)"; STEM=([A-Za-z0-9_]+) ;;$', sh, flags=re.M)
    return {name: (opt, stem) for name, opt, stem in rows}


def _args(sh: str) -> str:
    """The launcher's ARGS line, the file's command with $OPT for its option string."""
    found = re.findall(r'^ARGS="([^"]+)"$', sh, flags=re.M)
    assert len(found) == 1
    return found[0]


def _steps(name: str) -> list[str]:
    return ["gate", "tests", "OTK125", f"{name}a", f"{name}b", "done"]


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launchers_literal_table_is_the_builders_one_table():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    table = _table(sh)
    assert list(table) == NAMES == [n for n, _ in bt.HEART_R74_TABLE]
    assert [(n, v[0]) for n, v in table.items()] == [tuple(p) for p in bt.HEART_R74_TABLE]
    assert {n: v[1] for n, v in table.items()} == bt.HEART_R74_STEMS == {n: v[2] for n, v in TABLE.items()}
    assert _args(sh) == "--entry organ_scope $OPT --z-mode raw --captain-ruling R74"
    assert "heart " + _args(sh) == bt.HEART_R74_COMMAND.format(option="$OPT")
    assert '  *) echo "refused: FILE_NAME is not one of row R74\'s four files" >&2; exit 89 ;;' in sh
    code = _code(sh)
    assert code.count('OPT="--trend-strength ') == 4 and code.count("--trend-form category") == 5  # four and OTK125's
    # a literal table: no name is read by a pattern
    case = sh[sh.index('case "$FILE" in'):sh.index("esac", sh.index('case "$FILE" in'))]
    assert [line.split(")")[0].strip() for line in case.splitlines()[1:]] == [*NAMES, "*"]
    assert "BASH_REMATCH" not in sh and "=~" not in sh and "${FILE#" not in sh and "${FILE%" not in sh
    assert "HALF_SEED=1000001" in code and code.count("--seed") == 1  # halve.py's, at 1000001; no build seed


def test_the_identity_hashes_are_the_committed_list_and_r55s_command():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    hearth = _sums("hearth_8a8bc3a_h2.SHA256SUMS")
    assert hearth["tka/T2heart_organ_ts125_tk.h5ad"] == hearth["tkb/T2heart_organ_ts125_tk.h5ad"] == OTK125_SHA
    assert (hearth["tka/T2heart_organ_ts125_tk_half.h5ad"] == hearth["tkb/T2heart_organ_ts125_tk_half.h5ad"]
            == OTK125_HALF_SHA)
    assert f"OTK125_SHA={OTK125_SHA}" in sh and sh.count(OTK125_SHA) == 1
    assert f"OTK125_HALF_SHA={OTK125_HALF_SHA}" in sh and sh.count(OTK125_HALF_SHA) == 1
    assert f'OTK125_ARGS="{OTK125_ARGS}"' in sh
    assert "heart " + OTK125_ARGS == bt.HEART_TF_COMMANDS["category"]  # R55's command for OTK125
    wbh = WBH.read_text()  # the earlier heart launcher's own identity step, the same three
    assert f"OTK125_SHA={OTK125_SHA}" in wbh and f"OTK125_HALF_SHA={OTK125_HALF_SHA}" in wbh
    assert f'OTK125_ARGS="{OTK125_ARGS}"' in wbh


def test_the_steps_run_in_order_and_each_build_is_halved_checked_and_stamped():
    sh = _sh()
    marks = ("STEP=gate", "STEP=tests", "STEP=OTK125", "build OTK125 T2heart_organ_ts125_tk $OTK125_ARGS || finish $?",
             'for STEP in "${FILE}a" "${FILE}b"; do begin; build "$STEP" "$STEM" $ARGS || finish $?; step_ok; done',
             "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    build = sh[sh.index("build() {"):sh.index(': > "$S/stamps"')]
    lines = [ln.strip() for ln in build.splitlines()]
    want = ['$PY experiments/builds/build_test.py heart "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11',
            '$PY experiments/builds/halve.py "$d/$name.h5ad" "$d/${name}_half.h5ad" --seed "$HALF_SEED" > /dev/null '
            '2>&1 \\',
            "|| return 12",
            'uv run --frozen --no-sync vec check "$d/${name}_half.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13',
            'st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14',
            '[ "$st" = "$head True" ] || return 15']
    at = [lines.index(w) for w in want]
    assert at == sorted(at)
    ident = sh[sh.index("STEP=OTK125"):sh.index('for STEP in "${FILE}a"')]
    assert '[ "$(sha "$O/OTK125/T2heart_organ_ts125_tk.h5ad")" = "$OTK125_SHA" ]' in ident and "finish 20" in ident
    assert '[ "$(sha "$O/OTK125/T2heart_organ_ts125_tk_half.h5ad")" = "$OTK125_HALF_SHA" ]' in ident
    assert "finish 21" in ident
    done = sh[sh.index("STEP=done"):]
    assert '(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 30' in done
    assert 'for f in "${STEM}_half" "$STEM"; do' in done
    assert 'cmp -s "$O/${FILE}a/$f.h5ad" "$O/${FILE}b/$f.h5ad"' in done and "finish 31" in done
    assert '[ "$(sha "$O/${FILE}a/$STEM.h5ad")" != "$OTK125_SHA" ]' in done
    assert '[ "$(sha "$O/${FILE}a/${STEM}_half.h5ad")" != "$OTK125_HALF_SHA" ]' in done and "finish 32" in done
    assert "sha() { sha256sum \"$1\" | cut -d' ' -f1; }" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=OTK125")
    assert "export NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    first_uv = sh.index("uv run")
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < first_uv
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh and "BOARD=T2:heart:test_extrap" in sh
    assert '[ "$head" = "$COMMIT" ]' in sh
    assert "upload" not in sh.replace("Nothing here uploads", "").replace("goes up", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")
    assert "O=runs/heartr74_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/heartr74/$RUN" in sh


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for i, line in enumerate(lines):  # every builder, halve and check call: output to /dev/null
        if re.search(r'build_test\.py heart|halve\.py "|vec check "', line):
            assert "> /dev/null 2>&1" in line + lines[i + 1], line
    code = _code(sh)
    for word in (".build.json", "du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -",
                 "tail ", "cat "):
        assert word not in code, word
    assert "tests.out" in sh  # the tests' log: synthetic stages only
    assert 'print(s["commit"], s["clean"])' in sh and "stamp.json" in sh  # of the stamp, its commit and clean only
    assert "curl" not in sh and "wget" not in sh and "git fetch" not in sh and "git pull" not in sh
    assert "ssh" not in code and "scp" not in code


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    for step in ("tests", "OTK125", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    for code in ("exit 89", "exit 90", "exit 91"):
        assert sh.index(code) < sh.index("STEP=gate"), code
    assert sh.index("exit 89 ;;") < sh.index('mkdir -p "$HOME/v/scratch/heartr74"')  # every refusal before any folder
    assert sh.index("exit 90") < sh.index('mkdir -p "$HOME/v/scratch/heartr74"')
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    row = gate[gate.index("first_cell() {"):gate.index("\n}\n")]
    assert "[ \"$(grep -c \"^| $id \" experiments/harness/TRIES.md)\" = 1 ]" in row
    assert "first=$(grep \"^| $id \" experiments/harness/TRIES.md | cut -d'|' -f2)" in row
    assert "grep -q \"^ $id (the captain's \"" in row and 'grep -qF -- "ruling, no run"' in row
    assert "r74=$(first_cell R74 95) || exit $?" in gate
    assert 'for opt in "--entry organ_scope" "--z-mode raw" "$OPT" "--captain-ruling R74"; do' in gate
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in gate  # whole backticked strings, line by line
    assert "grep -qF -- \"\\`$OPT\\` ($FILE)\"" in gate  # the string paired with the file's name
    assert 'STAGES="E10.5 E9.5 E8.75 E8.5 E8.25_late"' in sh and '[ -f "data/raw/T2/$stage.h5ad" ]' in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 92, 88, 93]


def test_the_stages_are_the_ones_the_builder_reads():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    assert bt.DATA == REPO / "data" / "raw"
    want = {"T2/E10.5.h5ad", "T2/E9.5.h5ad", *(f for f, _ in bt.HEART_TREND_STAGES.values())}
    assert want == {f"T2/{s}.h5ad" for s in STAGES}


def test_the_tests_step_runs_the_new_tests_and_the_earlier_heart_lines():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    assert {"tests/test_heart_r74_line.py", "tests/test_heart_r74_build.py", "tests/test_heart_trend_line_h.py",
            "tests/test_heart_wbh_build.py", "tests/test_build_heart_organ_scope.py", "tests/test_halve.py",
            "tests/test_hidden_counts.py", "tests/test_s2_rules.py"} <= set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step and '--basetemp="$S/tmp/pytest"' in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = _code(_sh()).splitlines()
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    assert not re.search(r"(^|[^A-Za-z0-9_])rm\s", _code(_sh()), flags=re.M)
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 3
    for msg in usages:  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg
    assert not re.search(r"\$\{[A-Za-z_][A-Za-z0-9_]*:\?[^}]*'", _sh())  # nor in any other ${var:?message}
    assert "export PATH=$HOME/.local/bin:$PATH" in _sh()


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", sh)}
    documented = {int(c) for c in re.findall(r"(?:^|[\s,#])(\d\d) ", doc, flags=re.M)}
    assert used == documented, (used, documented)


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with rows R55's and R74's orders, not the repository's
    ledger.write_text(LEDGER + R55_ROW + ROW)
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)

    def parse(args: str):
        a = bt.build_parser().parse_args(["heart", *args.split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        return a

    a = parse(OTK125_ARGS)
    assert bt._trend_strength_refusals(a) == 1.25 and bt._r74_file(a) is None
    sh = _sh()
    args = _args(sh)
    for name, (opt, _) in _table(sh).items():
        a = parse(args.replace("$OPT", opt))
        assert bt._trend_strength_refusals(a) == float(TABLE[name][1]) != 1.25  # the table's strength
        assert bt._r74_file(a) == name


# ------------------------------------------------------------------ the launcher, run for real


def _env(tmp_path: Path, **extra) -> dict:
    import os

    env = {k: v for k, v in os.environ.items() if k not in GIT_LOCATION_VARS}
    env.update(HOME=str(tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path), GIT_CONFIG_GLOBAL=os.devnull,
               GIT_CONFIG_NOSYSTEM="1", UV_PROJECT_ENVIRONMENT=str(tmp_path / "elsewhere"))  # the launcher unsets it
    env.update(extra)
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW, data: tuple[str, ...] = STAGES, good_sha: bool = True,
          venv: str | None = "own", launcher: str | None = None) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (or the given copy of it) and the files its gate reads (the given
    released stages as empty stand-ins, and a .venv whose vec.pth names its src, or another checkout's src, or none),
    and a home whose ~/.local/bin holds a `uv` that fails (the launcher puts that folder first on PATH): returns the
    repository, the env and its HEAD."""
    env = _env(tmp_path)
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text("#!/bin/sh\nexit 1\n")
    (fake / "uv").chmod(0o755)
    repo = tmp_path / "repo"
    (repo / "experiments" / "builds").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "docs").mkdir()
    (repo / "src").mkdir()
    (repo / "src" / "keep.txt").write_text("src\n")
    (repo / "experiments" / "builds" / "heart_r74_line.sh").write_text(_sh() if launcher is None else launcher)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / "docs" / "S2_RULES.md").write_text(RULES)
    sha = hashlib.sha256(RULES.encode()).hexdigest() if good_sha else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    (repo / ".gitignore").write_text("data/\nruns/\n.venv/\ncalls.txt\n")
    (repo / "data" / "raw" / "T2").mkdir(parents=True)
    for stage in data:
        (repo / "data" / "raw" / "T2" / f"{stage}.h5ad").write_text("stand-in\n")
    if venv is not None:
        site = repo / ".venv" / "lib" / "python3.12" / "site-packages"
        site.mkdir(parents=True)
        src = repo / "src" if venv == "own" else tmp_path / "other" / "src"
        (site / "vec.pth").write_text(f"{src}\n")
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "experiments/builds/heart_r74_line.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True, timeout=120)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "heartr74" / run


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no head file, no later step's mark
    assert not list((s / "tmp").iterdir())


@pytest.mark.parametrize("name", NAMES)
def test_each_name_passes_the_gate_on_the_row_and_stops_at_the_failing_tests_step(tmp_path, name):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1", name)
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"heartr74_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", FIRST).returncode == 10
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "OTK10")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", FIRST)
    assert other.returncode == 10 and (repo / "runs" / f"heartr74_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_is_refused_at_the_gate(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1", FIRST)
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()


@pytest.mark.parametrize("argv, code", [
    ([], None), (["r1"], None), (["r1", ""], None),  # missing arguments (an empty name is none)
    (["a/b", FIRST], 90), (["r 1", FIRST], 90), (["$(id)", FIRST], 90),
    (["r1", FIRST, "extra"], 89),
    *[(["r1", bad], 89) for bad in (
        "OTK125", "OTK1", "OTK100", "OTK0.75", "OTK0.750", "OTK075 ", " OTK075", "otk10", "Otk10", "OTK10.h5ad",
        "OTK10_half", "OTK15a", "OTK17", "OTK175b", "OTK20", "OTK05", "OTK125h2", "OTK125h3", "OTS10", "OTP10",
        "OTKP125", "OTX10", "T2heart_organ_ts10_tk", "T2heart_organ_ts10_tk_half.h5ad", "10", "1.0", "ts10_tk",
        "$(id)", "*", "OTK*", "OTK1?")],
])
def test_odd_arguments_and_unknown_names_are_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, *argv)
    assert out.returncode != 0 if code is None else out.returncode == code, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


@pytest.mark.parametrize("commit", ["HEAD12", "", "g" * 40, "A" * 40])
def test_a_short_or_odd_commit_is_refused_before_any_folder_is_made(tmp_path, commit):
    repo, env, head = _repo(tmp_path)
    commit = head[:12] if commit == "HEAD12" else commit
    out = _run(repo, env, commit, "r1", FIRST)
    assert out.returncode == 90 if commit else out.returncode != 0, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW}more\n")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ""}, 95, "exactly one row R74"),
    ({"tries": ROW + ROW}, 95, "exactly one row R74"),  # two rows with one ID
    ({"tries": ROW.replace("| R74 (", "| R74b (")}, 95, "exactly one row R74"),  # another ID
    ({"tries": ROW.replace("| R74 (", "| R740 (")}, 95, "exactly one row R74"),
    ({"tries": ROW.replace("R74 (the captain's order", "R74 (a proposal")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("R74 (the captain's", "R74 awaits (the captain's")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("ruling, no run", "a try")}, 95, "does not say ruling, no run"),
    # the file's string inside a longer backticked string, or with no backticks, is not named
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", f"`{FIRST_OPT} --z-mode raw`")}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", FIRST_OPT)}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace("`--entry organ_scope`", "`--entry organ_scope --z-mode raw`")}, 95,
     "does not name `--entry organ_scope`"),
    ({"tries": ROW.replace("`--z-mode raw`", "--z-mode raw")}, 95, "does not name `--z-mode raw`"),
    ({"tries": ROW.replace("`--captain-ruling R74`", "`--captain-ruling R55`")}, 95,
     "does not name `--captain-ruling R74`"),
    ({"tries": ROW.replace("`--captain-ruling R74`", "`--captain-ruling R74 now`")}, 95,
     "does not name `--captain-ruling R74`"),
    # the string named, but paired with another file's name
    ({"tries": ROW.replace(f"`{FIRST_OPT}` ({FIRST})", f"`{FIRST_OPT}` (OTK175)")}, 95,
     f"does not pair `{FIRST_OPT}` with ({FIRST})"),
    ({"tries": (f"| R74 (the captain's order of X; ruling, no run) | heart | `--entry organ_scope`, `--z-mode raw`, "
                f"`{FIRST_OPT}` ({FIRST}) `--captain-ruling R74` |\n")}, 95, "does not name `--entry organ_scope`"),
    *[({"data": tuple(s for s in STAGES if s != gone)}, 92, f"data/raw/T2/{gone}.h5ad is not in this checkout")
      for gone in STAGES],
    ({"venv": None}, 88, "no .venv of its own"),
    ({"venv": "other"}, 88, "no .venv of its own"),  # a venv whose vec.pth names another checkout's src
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"heartr74_{head[:7]}_r1").exists()


@pytest.mark.parametrize("name", NAMES)
def test_the_gate_takes_a_short_row_in_its_own_words(tmp_path, name):
    tries = (f"| R74 (the captain's order of X, ruling, no run: `--entry organ_scope`, `--z-mode raw`, "
             f"`{TABLE[name][0]}` ({name}), `--captain-ruling R74`) | heart |\n")
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, head, "r1", name)
    assert out.returncode == 10, out.stderr  # past the gate: the stand-in uv fails the tests step
    assert (_scratch(tmp_path, "r1") / "gate.exit").read_text().strip() == "0"
    other = [f for f in NAMES if f != name][0]  # a file the short row does not name
    out = _run(repo, env, head, "r2", other)
    assert out.returncode == 95 and f"does not name `{TABLE[other][0]}`" in out.stderr


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"heartr74_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


# A stand-in `uv` for runs past the tests step: the tests pass; a build writes "file <its output's name>" (with
# SAME_AS_OTK125 set, OTK125's) and a stamp naming HEAD and clean (or CLEAN's word); halve.py writes "half <seed>" and
# its input's bytes; `vec check` passes; BUILD_FAILS, HALVE_FAILS and CHECK_FAILS fail their calls, and with B_DIFFERS
# set the second build of the file writes other bytes. Every builder, halve and check call is listed in calls.txt; the
# stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec)
    echo "$*" >> calls.txt
    [ -n "${CHECK_FAILS:-}" ] && exit 1
    exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/build_test.py)
        echo "$*" >> calls.txt
        [ -n "${BUILD_FAILS:-}" ] && exit 1
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        [ -n "${SAME_AS_OTK125:-}" ] && name=T2heart_organ_ts125_tk
        case "$out" in *b/*) [ -n "${B_DIFFERS:-}" ] && name="$name other" ;; esac
        echo "file $name" > "$out"
        printf '{"commit": "%s", "clean": %s}' "$(git rev-parse HEAD)" "${CLEAN:-true}" > "$out.stamp.json"
        exit 0 ;;
      experiments/builds/halve.py)
        echo "$*" >> calls.txt
        [ -n "${HALVE_FAILS:-}" ] && exit 1
        seed=""; prev=""
        for a in "$@"; do [ "$prev" = "--seed" ] && seed=$a; prev=$a; done
        { echo "half $seed"; cat "$2"; } > "$3"
        exit 0 ;;
    esac ;;
esac
exit 1
"""
OTK125_STAND_IN = hashlib.sha256(b"file T2heart_organ_ts125_tk\n").hexdigest()
OTK125_HALF_STAND_IN = hashlib.sha256(b"half 1000001\nfile T2heart_organ_ts125_tk\n").hexdigest()


def _stand_in(sh: str, keep: tuple[str, ...] = ()) -> str:
    """The launcher with its identity hashes set to the stand-in files' (those of the variables in keep left as they
    are)."""
    for var, sha, stand_in in (("OTK125_SHA", OTK125_SHA, OTK125_STAND_IN),
                               ("OTK125_HALF_SHA", OTK125_HALF_SHA, OTK125_HALF_STAND_IN)):
        if var not in keep:
            sh = sh.replace(f"{var}={sha}", f"{var}={stand_in}")
    return sh


def _stand_in_run(tmp_path: Path, name: str, keep: tuple[str, ...] = (), **env_extra):
    sh = _stand_in(_sh(), keep)
    assert sh != _sh() or keep
    repo, env, head = _repo(tmp_path, launcher=sh)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(FULL_UV)
    return repo, head, _run(repo, {**env, **env_extra}, head, "r1", name)


def _exits(tmp_path: Path) -> dict[str, str]:
    return {p.name[:-5]: p.read_text().strip() for p in _scratch(tmp_path, "r1").glob("*.exit")}


@pytest.mark.parametrize("name", NAMES)
def test_a_whole_stand_in_run_passes_every_step_in_order(tmp_path, name):
    """With only its identity hashes set to the stand-in files' (in a copy) and the stand-in `uv` above: a call runs
    every step in order to done.exit 0, each step's exit file says 0, the stamps list every build, the builder, halve.py
    and `vec check` get the commands their steps name, the hash list holds every built file once and the two builds of
    the file are one file."""
    repo, head, out = _stand_in_run(tmp_path, name)
    assert out.returncode == 0, out.stderr
    s = _scratch(tmp_path, "r1")
    assert _exits(tmp_path) == {step: "0" for step in _steps(name)}
    r = f"runs/heartr74_{head[:7]}_r1"
    o = repo / r
    folders = _steps(name)[2:-1]
    assert sorted(p.name for p in o.iterdir()) == sorted([*folders, "SHA256SUMS"])
    stamps = (s / "stamps").read_text().splitlines()
    assert [line.split()[0].split("/")[0] for line in stamps] == folders
    assert all(line.endswith(f" {head} True") for line in stamps)
    listed = (o / "SHA256SUMS").read_text().splitlines()
    stem = TABLE[name][2]
    files = ["OTK125/T2heart_organ_ts125_tk.h5ad", "OTK125/T2heart_organ_ts125_tk_half.h5ad",
             *(f"{d}/{stem}{h}.h5ad" for d in folders[1:] for h in ("", "_half"))]
    assert [line.split()[1] for line in listed] == files
    assert (s / "SHA256SUMS").read_text() == (o / "SHA256SUMS").read_text()
    sums = dict(line.split()[::-1] for line in listed)
    assert sums["OTK125/T2heart_organ_ts125_tk.h5ad"] == OTK125_STAND_IN
    assert sums[f"{name}a/{stem}.h5ad"] == sums[f"{name}b/{stem}.h5ad"] != OTK125_STAND_IN
    args = f"--entry organ_scope {TABLE[name][0]} --z-mode raw --captain-ruling R74"
    want = []
    for d, n, a in (("OTK125", "T2heart_organ_ts125_tk", OTK125_ARGS), (f"{name}a", stem, args),
                    (f"{name}b", stem, args)):
        want += [f"experiments/builds/build_test.py heart {a} --out {r}/{d}/{n}.h5ad",
                 f"experiments/builds/halve.py {r}/{d}/{n}.h5ad {r}/{d}/{n}_half.h5ad --seed 1000001",
                 f"vec check {r}/{d}/{n}_half.h5ad --board T2:heart:test_extrap"]
    assert (repo / "calls.txt").read_text().splitlines() == want
    for f in (f"{stem}.h5ad", f"{stem}_half.h5ad"):
        assert (o / f"{name}a" / f).read_bytes() == (o / f"{name}b" / f).read_bytes()


def test_a_call_past_the_tests_step_stops_at_otk125s_recorded_bytes(tmp_path):
    """The real hashes: no stand-in file has OTK125's recorded bytes, so the run stops at the identity step."""
    repo, head, out = _stand_in_run(tmp_path, FIRST, keep=("OTK125_SHA", "OTK125_HALF_SHA"))
    assert out.returncode == 20, out.stderr
    assert "stopped: OTK125's command does not give its recorded file" in out.stderr
    assert _exits(tmp_path) == {"gate": "0", "tests": "0", "OTK125": "20"}
    assert (_scratch(tmp_path, "r1") / "stamps").read_text() == f"OTK125/T2heart_organ_ts125_tk {head} True\n"
    assert sorted(p.name for p in (repo / "runs" / f"heartr74_{head[:7]}_r1").iterdir()) == ["OTK125"]


def test_a_half_of_other_bytes_stops_the_run_at_the_identity_step(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, FIRST, keep=("OTK125_HALF_SHA",))
    assert out.returncode == 21, out.stderr
    assert "stopped: OTK125's half is not its recorded half" in out.stderr
    assert _exits(tmp_path) == {"gate": "0", "tests": "0", "OTK125": "21"}


@pytest.mark.parametrize("env_extra, code", [
    ({"BUILD_FAILS": "1"}, 11), ({"HALVE_FAILS": "1"}, 12), ({"CHECK_FAILS": "1"}, 13), ({"CLEAN": "false"}, 15),
])
def test_a_failing_build_halve_check_or_stamp_stops_the_identity_step(tmp_path, env_extra, code):
    repo, head, out = _stand_in_run(tmp_path, FIRST, **env_extra)
    assert out.returncode == code, out.stderr
    assert _exits(tmp_path) == {"gate": "0", "tests": "0", "OTK125": str(code)}
    assert not (repo / "runs" / f"heartr74_{head[:7]}_r1" / f"{FIRST}a").exists()


def test_two_builds_of_other_bytes_stop_the_run_at_done(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, FIRST, B_DIFFERS="1")
    assert out.returncode == 31, out.stderr
    assert _exits(tmp_path) == {**{step: "0" for step in _steps(FIRST)[:-1]}, "done": "31"}
    assert (_scratch(tmp_path, "r1") / "SHA256SUMS").exists()


@pytest.mark.parametrize("name", NAMES)
def test_a_file_with_otk125s_bytes_stops_the_run_at_done(tmp_path, name):
    """The trap at run time: a build that took R55's strength would give OTK125's bytes; the run says so (32)."""
    repo, head, out = _stand_in_run(tmp_path, name, SAME_AS_OTK125="1")
    assert out.returncode == 32, out.stderr
    assert f"stopped: {name} or its half has OTK125's bytes" in out.stderr
    assert _exits(tmp_path) == {**{step: "0" for step in _steps(name)[:-1]}, "done": "32"}
