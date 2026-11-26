"""The test-board build launcher of row R77's two files (OTK125's command with the category form's strength from the
builder's HEART_R77_TABLE), one file per call: experiments/builds/heart_r77_line.sh <COMMIT> <RUN> <FILE_NAME>.

It refuses any name but the two of its literal table, odd arguments, a reused run folder and an odd commit before it
makes a folder (OTK125, row R74's four names, OTK5, OTK200 and otk05 among them); its gate checks the commit, the clean
tree, the rules file's hash, row R77's first cell (`--entry organ_scope`, `--z-mode raw`, the file's option string and
`--captain-ruling R77` as whole backticked strings, the option string paired there with the file's name), the five
released stages and the venv; then it runs the targeted tests, rebuilds OTK125 and its half (halve.py --seed 1000001)
and stops unless both give their recorded sha256; then it builds the file twice from one clean commit, halves each at
seed 1000001, checks each half, compares the two builds, and refuses a file or half with OTK125's bytes, with the
recorded bytes of a file of row R74, or with those of row R77's other file as any list of this line records them (a run
at the same time included). These tests hold its table to the builder's, its hashes to the committed lists, its
commands to the ones the builder takes and its gate, identity and build steps to row R74's launcher's; check that it
deletes nothing and reads no output of a test build but exit codes, sha256 and the stamps' commit and clean fields; and
run it for real in throwaway repositories under a temporary home: with a stand-in `uv` that fails (every refusal path,
and a call passing the gate stops at the tests step and builds nothing) and with stand-ins that pass (a whole run of
each name goes through every step in order; a stand-in of other bytes, a failing build, halve or check, an unclean
stamp, and a file with OTK125's, an R74 file's or the other R77 file's recorded bytes stop the run where their step
says). Nothing here reads the repository's own ledger or rules file; build_test is imported only to check its table and
that it takes the commands the launcher runs, against a temporary ledger."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import GIT_LOCATION_VARS
from test_heart_r74_build import R74_FIRST
from test_heart_r74_build import TABLE as R74_TABLE
from test_heart_r77_build import NAMES, R77_FIRST, TABLE
from test_heart_trend_line_g import LEDGER
from test_heart_trend_line_h import R55_ROW

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "heart_r77_line.sh"
R74_LAUNCHER = REPO / "experiments" / "builds" / "heart_r74_line.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
OTK125_SHA = "d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893"
OTK125_HALF_SHA = "9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54"
OTK125_ARGS = "--entry organ_scope --trend-strength 1.25 --trend-form category --z-mode raw --captain-ruling R55"
R74_LISTS = {"OTK10": "heartr74_4b3e7a9_r74k10.SHA256SUMS", "OTK15": "heartr74_4b3e7a9_r74k15.SHA256SUMS",
             "OTK075": "heartr74_4b3e7a9_r74k075.SHA256SUMS", "OTK175": "heartr74_4b3e7a9_r74k175.SHA256SUMS"}
STAGES = ("E10.5", "E9.5", "E8.75", "E8.5", "E8.25_late")
ROW = f"| {R77_FIRST} | heart | refinement round 4's two files |\n"
FIRST = "OTK20"  # the file the gate's refusals are run on
FIRST_OPT = TABLE[FIRST][0]
OTHER = {"OTK05": "OTK20", "OTK20": "OTK05"}
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


def _r74_shas(sh: str) -> list[str]:
    """The launcher's R74_SHAS list, in its order."""
    found = re.findall(r'^R74_SHAS="([0-9a-f\n]+)"$', sh, flags=re.M)
    assert len(found) == 1
    return found[0].split()


def _steps(name: str) -> list[str]:
    return ["gate", "tests", "OTK125", f"{name}a", f"{name}b", "done"]


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


def _between(sh: str, start: str, end: str) -> str:
    return sh[sh.index(start):sh.index(end, sh.index(start))]


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launchers_literal_table_is_the_builders_one_table():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    table = _table(sh)
    assert list(table) == NAMES == [n for n, _ in bt.HEART_R77_TABLE]
    assert [(n, v[0]) for n, v in table.items()] == [tuple(p) for p in bt.HEART_R77_TABLE]
    assert {n: v[1] for n, v in table.items()} == bt.HEART_R77_STEMS == {n: v[2] for n, v in TABLE.items()}
    assert _args(sh) == "--entry organ_scope $OPT --z-mode raw --captain-ruling R77"
    assert "heart " + _args(sh) == bt.HEART_R77_COMMAND.format(option="$OPT")
    assert '  *) echo "refused: FILE_NAME is not one of row R77\'s two files" >&2; exit 89 ;;' in sh
    assert f'R77_STEMS="{" ".join(bt.HEART_R77_STEMS.values())}"' in sh  # both stems, the table's order
    code = _code(sh)
    assert code.count('OPT="--trend-strength ') == 2 and code.count("--trend-form category") == 3  # two and OTK125's
    # a literal table: no name is read by a pattern
    case = sh[sh.index('case "$FILE" in'):sh.index("esac", sh.index('case "$FILE" in'))]
    assert [line.split(")")[0].strip() for line in case.splitlines()[1:]] == [*NAMES, "*"]
    assert "BASH_REMATCH" not in sh and "=~" not in sh and "${FILE#" not in sh and "${FILE%" not in sh
    assert "HALF_SEED=1000001" in code and code.count("--seed") == 1  # halve.py's, at 1000001; no build seed
    assert not set(table) & set(R74_TABLE) and not {v[0] for v in table.values()} & {v[0] for v in R74_TABLE.values()}


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
    r74 = R74_LAUNCHER.read_text()  # row R74's launcher's identity step, the same three
    assert f"OTK125_SHA={OTK125_SHA}" in r74 and f"OTK125_HALF_SHA={OTK125_HALF_SHA}" in r74
    assert f'OTK125_ARGS="{OTK125_ARGS}"' in r74


def test_the_r74_hashes_are_its_four_committed_lists():
    """R74_SHAS holds exactly row R74's four files and their halves as their committed lists record them (a and b
    alike), file then half in row R74's order; none is OTK125's."""
    want = []
    for name, (_, _, stem) in R74_TABLE.items():
        sums = _sums(R74_LISTS[name])
        for suffix in ("", "_half"):
            a, b = sums[f"{name}a/{stem}{suffix}.h5ad"], sums[f"{name}b/{stem}{suffix}.h5ad"]
            assert a == b
            want.append(a)
        assert sums["OTK125/T2heart_organ_ts125_tk.h5ad"] == OTK125_SHA
        assert sums["OTK125/T2heart_organ_ts125_tk_half.h5ad"] == OTK125_HALF_SHA
    assert sorted(R74_LISTS.values()) == sorted(p.name for p in RESULTS.glob("heartr74_*.SHA256SUMS"))
    assert _r74_shas(_sh()) == want and len(set(want)) == 8 and not {OTK125_SHA, OTK125_HALF_SHA} & set(want)


def test_the_gate_identity_and_build_steps_are_row_r74s_launchers():
    """Copied from row R74's launcher, not refactored: the gate (row R77 for R74), the identity step, the build and
    stamp functions and the two builds are R74's lines; the tests step lists R74's files and this line's two."""
    sh, r74 = _sh(), R74_LAUNCHER.read_text()
    assert _between(sh, "STEP=gate", "STEP=tests") == _between(r74, "STEP=gate", "STEP=tests").replace(
        "R74", "R77").replace("r74", "r77")
    assert _between(sh, "set -u\n", 'case "$FILE" in') == _between(r74, "set -u\n", 'case "$FILE" in').replace(
        "r74_line", "r77_line").replace("R74", "R77")
    assert _between(sh, "HALF_SEED=", "# row R74's four files") == _between(r74, "HALF_SEED=", "STAGES=").replace(
        "heartr74", "heartr77")
    assert _between(sh, "export OMP_NUM_THREADS=16", "recorded_other() {") == _between(
        r74, "export OMP_NUM_THREADS=16", ': > "$S/stamps"')
    assert _between(sh, ': > "$S/stamps"', "STEP=done") == _between(r74, ': > "$S/stamps"', "STEP=done")
    assert _between(sh, "STEP=done", "file_sha=") == _between(r74, "STEP=done", '[ "$(sha "$O/${FILE}a/$STEM.h5ad")"')
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", _between(sh, "STEP=tests", "step_ok"))
    r74_listed = re.findall(r"tests/test_[a-z0-9_]+\.py", _between(r74, "STEP=tests", "step_ok"))
    assert listed == ["tests/test_heart_r77_line.py", "tests/test_heart_r77_build.py", *r74_listed]


def test_the_steps_run_in_order_and_each_build_is_halved_checked_and_stamped():
    sh = _sh()
    marks = ("STEP=gate", "STEP=tests", "STEP=OTK125", "build OTK125 T2heart_organ_ts125_tk $OTK125_ARGS || finish $?",
             'for STEP in "${FILE}a" "${FILE}b"; do begin; build "$STEP" "$STEM" $ARGS || finish $?; step_ok; done',
             "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    build = sh[sh.index("build() {"):sh.index("recorded_other() {")]
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
    steps = ['(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 30',
             'cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 30',
             'for f in "${STEM}_half" "$STEM"; do',
             'cmp -s "$O/${FILE}a/$f.h5ad" "$O/${FILE}b/$f.h5ad"',
             'file_sha=$(sha "$O/${FILE}a/$STEM.h5ad") || finish 30',
             'half_sha=$(sha "$O/${FILE}a/${STEM}_half.h5ad") || finish 30',
             '[ "$file_sha" != "$OTK125_SHA" ] && [ "$half_sha" != "$OTK125_HALF_SHA" ]',
             "finish 32", "for x in $R74_SHAS; do", '[ "$file_sha" != "$x" ] && [ "$half_sha" != "$x" ]', "finish 33",
             "other=$(recorded_other)", "finish 34", "for x in $other; do",
             '[ "$file_sha" != "$x" ] && [ "$half_sha" != "$x" ]', "finish 34", "step_ok"]
    at, pos = [], 0
    for s in steps:  # each in this order, the list written before any comparison
        pos = done.index(s, pos)
        at.append(pos)
    assert "finish 31" in done
    other = sh[sh.index("recorded_other() {"):sh.index(': > "$S/stamps"')]
    assert ('for f in "$HOME"/v/scratch/heartr77/*/SHA256SUMS experiments/builds/results/heartr77_*.SHA256SUMS; '
            'do') in other
    assert '[ -f "$f" ] || continue' in other and "for s in $R77_STEMS; do" in other
    assert '[ "$s" = "$STEM" ] && continue' in other  # this file's own lines (a run of it again) are not compared
    assert ("awk -v a=\"$s.h5ad\" -v b=\"${s}_half.h5ad\" '{ n = $2; sub(/.*\\//, \"\", n); if (n == a || n == b) "
            "print $1 }' \"$f\"") in other  # the other file's name and its half's, whole, after the folder
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
    assert "O=runs/heartr77_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/heartr77/$RUN" in sh


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
    # awk reads sha256 lists only (this line's own and the committed ones), never a built file
    assert code.count("awk ") == 1 and '"$f"' in code[code.index("awk "):].split("\n")[0]


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    for step in ("tests", "OTK125", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    for code in ("exit 89", "exit 90", "exit 91"):
        assert sh.index(code) < sh.index("STEP=gate"), code
    assert sh.index("exit 89 ;;") < sh.index('mkdir -p "$HOME/v/scratch/heartr77"')  # every refusal before any folder
    assert sh.index("exit 90") < sh.index('mkdir -p "$HOME/v/scratch/heartr77"')
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
    assert "r77=$(first_cell R77 95) || exit $?" in gate
    assert 'for opt in "--entry organ_scope" "--z-mode raw" "$OPT" "--captain-ruling R77"; do' in gate
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in gate  # whole backticked strings, line by line
    assert "grep -qF -- \"\\`$OPT\\` ($FILE)\"" in gate  # the string paired with the file's name
    assert 'STAGES="E10.5 E9.5 E8.75 E8.5 E8.25_late"' in sh and '[ -f "data/raw/T2/$stage.h5ad" ]' in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 92, 88, 93]
    assert "R74" not in gate


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
    assert {"tests/test_heart_r77_line.py", "tests/test_heart_r77_build.py", "tests/test_heart_r74_line.py",
            "tests/test_heart_r74_build.py", "tests/test_heart_trend_line_h.py", "tests/test_heart_wbh_build.py",
            "tests/test_build_heart_organ_scope.py", "tests/test_halve.py", "tests/test_hidden_counts.py",
            "tests/test_s2_rules.py"} <= set(listed)
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

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with rows R55's, R74's and R77's orders, not the repository's
    ledger.write_text(LEDGER + R55_ROW + f"| {R74_FIRST} | heart | x |\n" + ROW)
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)

    def parse(args: str):
        a = bt.build_parser().parse_args(["heart", *args.split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        return a

    a = parse(OTK125_ARGS)
    assert bt._trend_strength_refusals(a) == 1.25 and bt._r74_file(a) is None and bt._r77_file(a) is None
    sh = _sh()
    args = _args(sh)
    for name, (opt, _) in _table(sh).items():
        a = parse(args.replace("$OPT", opt))
        assert bt._trend_strength_refusals(a) == float(TABLE[name][1]) != 1.25  # the table's strength
        assert bt._r77_file(a) == name and bt._r74_file(a) is None


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
          venv: str | None = "own", launcher: str | None = None,
          extra: dict[str, str] | None = None) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (or the given copy of it), the files its gate reads (the given
    released stages as empty stand-ins, and a .venv whose vec.pth names its src, or another checkout's src, or none)
    and any extra committed files, and a home whose ~/.local/bin holds a `uv` that fails (the launcher puts that
    folder first on PATH): returns the repository, the env and its HEAD."""
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
    (repo / "experiments" / "builds" / "heart_r77_line.sh").write_text(_sh() if launcher is None else launcher)
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
    for rel, text in (extra or {}).items():
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_text(text)
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "experiments/builds/heart_r77_line.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True, timeout=120)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "heartr77" / run


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
    o = repo / "runs" / f"heartr77_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", FIRST).returncode == 10
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "OTK05")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", FIRST)
    assert other.returncode == 10 and (repo / "runs" / f"heartr77_{head[:7]}_r2").is_dir()


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
        "OTK125", "OTK5", "OTK2", "OTK050", "OTK200", "OTK0.5", "OTK2.0", "OTK05 ", " OTK05", "otk05", "Otk20",
        "OTK05.h5ad", "OTK05_half", "OTK05a", "OTK20b", *R74_TABLE, "OTK0875", "OTK1125", "OTK1375", "OTK1625",
        "OTS05", "OTS20", "OTP10", "OTKP125", "T2heart_organ_ts05_tk", "T2heart_organ_ts20_tk_half.h5ad", "05", "0.5",
        "2.0", "ts20_tk", "$(id)", "*", "OTK*", "OTK0?", "OTK2?")],
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
    ({"tries": ""}, 95, "exactly one row R77"),
    ({"tries": ROW + ROW}, 95, "exactly one row R77"),  # two rows with one ID
    ({"tries": ROW.replace("| R77 (", "| R77b (")}, 95, "exactly one row R77"),  # another ID
    ({"tries": ROW.replace("| R77 (", "| R770 (")}, 95, "exactly one row R77"),
    ({"tries": f"| {R74_FIRST} | heart | x |\n"}, 95, "exactly one row R77"),  # row R74 alone
    ({"tries": ROW.replace("R77 (the captain's order", "R77 (a proposal")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("R77 (the captain's", "R77 awaits (the captain's")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("ruling, no run", "a try")}, 95, "does not say ruling, no run"),
    # the file's string inside a longer backticked string, or with no backticks, is not named
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", f"`{FIRST_OPT} --z-mode raw`")}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", FIRST_OPT)}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", f"`{FIRST_OPT.replace('2.0', '2.00')}`")}, 95,
     f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace("`--entry organ_scope`", "`--entry organ_scope --z-mode raw`")}, 95,
     "does not name `--entry organ_scope`"),
    ({"tries": ROW.replace("`--z-mode raw`", "--z-mode raw")}, 95, "does not name `--z-mode raw`"),
    ({"tries": ROW.replace("`--captain-ruling R77`", "`--captain-ruling R74`")}, 95,
     "does not name `--captain-ruling R77`"),
    ({"tries": ROW.replace("`--captain-ruling R77`", "`--captain-ruling R77 now`")}, 95,
     "does not name `--captain-ruling R77`"),
    # the string named, but paired with another file's name
    ({"tries": ROW.replace(f"`{FIRST_OPT}` ({FIRST})", f"`{FIRST_OPT}` (OTK05)")}, 95,
     f"does not pair `{FIRST_OPT}` with ({FIRST})"),
    ({"tries": (f"| R77 (the captain's order of X; ruling, no run) | heart | `--entry organ_scope`, `--z-mode raw`, "
                f"`{FIRST_OPT}` ({FIRST}) `--captain-ruling R77` |\n")}, 95, "does not name `--entry organ_scope`"),
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
    assert not (repo / "runs" / f"heartr77_{head[:7]}_r1").exists()


@pytest.mark.parametrize("name", NAMES)
def test_the_gate_takes_a_short_row_in_its_own_words(tmp_path, name):
    tries = (f"| R77 (the captain's order of X, ruling, no run: `--entry organ_scope`, `--z-mode raw`, "
             f"`{TABLE[name][0]}` ({name}), `--captain-ruling R77`) | heart |\n")
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, head, "r1", name)
    assert out.returncode == 10, out.stderr  # past the gate: the stand-in uv fails the tests step
    assert (_scratch(tmp_path, "r1") / "gate.exit").read_text().strip() == "0"
    out = _run(repo, env, head, "r2", OTHER[name])  # a file the short row does not name
    assert out.returncode == 95 and f"does not name `{TABLE[OTHER[name]][0]}`" in out.stderr


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"heartr77_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


# A stand-in `uv` for runs past the tests step: the tests pass; a build writes "file <its output's name>" (with SAME_AS
# set, every build but OTK125's writes "file <SAME_AS>", another file's stand-in bytes) and a stamp naming HEAD and
# clean (or CLEAN's word); halve.py writes "half <seed>" and its input's bytes; `vec check` passes; BUILD_FAILS,
# HALVE_FAILS and CHECK_FAILS fail their calls, and with B_DIFFERS set the second build of the file writes other bytes.
# Every builder, halve and check call is listed in calls.txt; the stamp reader is python itself.
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
        case "$out" in */OTK125/*) ;; *) [ -n "${SAME_AS:-}" ] && name=$SAME_AS ;; esac
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


def _file_sha(stem: str) -> str:
    """The stand-in build's bytes for a file of this stem, and its half's at seed 1000001."""
    return hashlib.sha256(f"file {stem}\n".encode()).hexdigest()


def _half_sha(stem: str) -> str:
    return hashlib.sha256(f"half 1000001\nfile {stem}\n".encode()).hexdigest()


OTK125_STAND_IN = _file_sha("T2heart_organ_ts125_tk")
OTK125_HALF_STAND_IN = _half_sha("T2heart_organ_ts125_tk")


def _stand_in(sh: str, keep: tuple[str, ...] = (), r74: tuple[str, ...] = ()) -> str:
    """The launcher with its identity hashes set to the stand-in files' (those of the variables in keep left as they
    are), and with R74_SHAS holding the stand-in hashes of row R74's files (file and half, in its order) where r74
    names "file" or "half"."""
    for var, sha, stand_in in (("OTK125_SHA", OTK125_SHA, OTK125_STAND_IN),
                               ("OTK125_HALF_SHA", OTK125_HALF_SHA, OTK125_HALF_STAND_IN)):
        if var not in keep:
            sh = sh.replace(f"{var}={sha}", f"{var}={stand_in}")
    if r74:
        old = _r74_shas(sh)
        new = []
        for i, (_, _, stem) in enumerate(R74_TABLE.values()):
            new += [_file_sha(stem) if "file" in r74 else old[2 * i],
                    _half_sha(stem) if "half" in r74 else old[2 * i + 1]]
        sh = sh.replace('R74_SHAS="' + "\n".join(old) + '"', 'R74_SHAS="' + "\n".join(new) + '"')
        assert _r74_shas(sh) == new
    return sh


def _stand_in_repo(tmp_path: Path, keep: tuple[str, ...] = (), r74: tuple[str, ...] = (), **kw):
    sh = _stand_in(_sh(), keep, r74)
    assert sh != _sh() or keep
    repo, env, head = _repo(tmp_path, launcher=sh, **kw)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(FULL_UV)
    return repo, env, head


def _stand_in_run(tmp_path: Path, name: str, keep: tuple[str, ...] = (), r74: tuple[str, ...] = (), **env_extra):
    repo, env, head = _stand_in_repo(tmp_path, keep, r74)
    return repo, head, _run(repo, {**env, **env_extra}, head, "r1", name)


def _exits(tmp_path: Path, run: str = "r1") -> dict[str, str]:
    return {p.name[:-5]: p.read_text().strip() for p in _scratch(tmp_path, run).glob("*.exit")}


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
    r = f"runs/heartr77_{head[:7]}_r1"
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
    assert sums[f"{name}a/{stem}.h5ad"] == sums[f"{name}b/{stem}.h5ad"] == _file_sha(stem) != OTK125_STAND_IN
    args = f"--entry organ_scope {TABLE[name][0]} --z-mode raw --captain-ruling R77"
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
    assert sorted(p.name for p in (repo / "runs" / f"heartr77_{head[:7]}_r1").iterdir()) == ["OTK125"]


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
    assert not (repo / "runs" / f"heartr77_{head[:7]}_r1" / f"{FIRST}a").exists()


def test_two_builds_of_other_bytes_stop_the_run_at_done(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, FIRST, B_DIFFERS="1")
    assert out.returncode == 31, out.stderr
    assert _exits(tmp_path) == {**{step: "0" for step in _steps(FIRST)[:-1]}, "done": "31"}
    assert (_scratch(tmp_path, "r1") / "SHA256SUMS").exists()


@pytest.mark.parametrize("name", NAMES)
def test_a_file_with_otk125s_bytes_stops_the_run_at_done(tmp_path, name):
    """The trap at run time: a build that took R55's strength would give OTK125's bytes; the run says so (32)."""
    repo, head, out = _stand_in_run(tmp_path, name, SAME_AS="T2heart_organ_ts125_tk")
    assert out.returncode == 32, out.stderr
    assert f"stopped: {name} or its half has OTK125's bytes" in out.stderr
    assert _exits(tmp_path) == {**{step: "0" for step in _steps(name)[:-1]}, "done": "32"}


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("r74_name, which", [*((n, ("file", "half")) for n in R74_TABLE), ("OTK10", ("half",)),
                                             ("OTK175", ("file",))])
def test_a_file_with_an_r74_files_recorded_bytes_stops_the_run_at_done(tmp_path, name, r74_name, which):
    """R74_SHAS (in a copy) holding the stand-in hashes of row R74's files, of the files and halves, of the halves
    only, or of the files only: a build with an R74 file's bytes stops at done (33), its list written first."""
    repo, head, out = _stand_in_run(tmp_path, name, r74=which, SAME_AS=R74_TABLE[r74_name][2])
    assert out.returncode == 33, out.stderr
    assert f"stopped: {name} or its half has the recorded bytes of a file of row R74" in out.stderr
    assert _exits(tmp_path) == {**{step: "0" for step in _steps(name)[:-1]}, "done": "33"}
    assert (_scratch(tmp_path, "r1") / "SHA256SUMS").exists()


def test_a_whole_run_passes_with_row_r74s_stand_in_hashes_of_other_files(tmp_path):
    """The same copy, a file of its own bytes: R74's hashes do not stop it."""
    repo, head, out = _stand_in_run(tmp_path, FIRST, r74=("file", "half"))
    assert out.returncode == 0, out.stderr


def _sibling_list(stem: str, name: str, half_only: bool = False) -> str:
    """A list as this line writes it for a run of file `name` (stem `stem`) whose builds have the stand-in bytes."""
    lines = [f"{OTK125_STAND_IN}  OTK125/T2heart_organ_ts125_tk.h5ad",
             f"{OTK125_HALF_STAND_IN}  OTK125/T2heart_organ_ts125_tk_half.h5ad"]
    for d in (f"{name}a", f"{name}b"):
        if not half_only:
            lines.append(f"{_file_sha(stem)}  {d}/{stem}.h5ad")
        lines.append(f"{_half_sha(stem)}  {d}/{stem}_half.h5ad")
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("half_only", [False, True])
def test_the_other_files_bytes_in_another_runs_list_stop_the_run_at_done(tmp_path, name, half_only):
    """Another run's scratch list (a run of the other file, at the same time or before) records the other file's
    bytes; a build of this file with those bytes stops at done (34)."""
    other = OTHER[name]
    repo, env, head = _stand_in_repo(tmp_path)
    s0 = _scratch(tmp_path, "r0")
    s0.mkdir(parents=True)
    (s0 / "SHA256SUMS").write_text(_sibling_list(TABLE[other][2], other, half_only))
    out = _run(repo, {**env, "SAME_AS": TABLE[other][2]}, head, "r1", name)
    assert out.returncode == 34, out.stderr
    assert f"stopped: {name} or its half has the recorded bytes of row R77's other file" in out.stderr
    assert _exits(tmp_path) == {**{step: "0" for step in _steps(name)[:-1]}, "done": "34"}


@pytest.mark.parametrize("name", NAMES)
def test_the_other_files_bytes_in_a_committed_list_stop_the_run_at_done(tmp_path, name):
    other = OTHER[name]
    listed = {f"experiments/builds/results/heartr77_0000000_{other.lower()}.SHA256SUMS":
              _sibling_list(TABLE[other][2], other)}
    repo, env, head = _stand_in_repo(tmp_path, extra=listed)
    out = _run(repo, {**env, "SAME_AS": TABLE[other][2]}, head, "r1", name)
    assert out.returncode == 34, out.stderr
    assert _exits(tmp_path) == {**{step: "0" for step in _steps(name)[:-1]}, "done": "34"}
    out = _run(repo, env, head, "r2", name)  # its own bytes: the committed list does not stop it
    assert out.returncode == 0, out.stderr


def test_both_files_in_turn_and_a_file_again_pass(tmp_path):
    """Each file run once, then one again under a new run name: every list holds the files' own bytes, so nothing
    stops (a list's lines of this file itself are never compared)."""
    repo, env, head = _stand_in_repo(tmp_path)
    for run, name in (("r1", "OTK05"), ("r2", "OTK20"), ("r3", "OTK05")):
        out = _run(repo, env, head, run, name)
        assert out.returncode == 0, (run, out.stderr)
        assert _exits(tmp_path, run)["done"] == "0"


def test_two_runs_at_the_same_time_with_one_files_bytes_stop_at_done(tmp_path):
    """Both files at once, the OTK05 call building OTK20's stand-in bytes: each run writes its list before it compares,
    so at least one of the two sees the other's list and stops (34)."""
    repo, env, head = _stand_in_repo(tmp_path)
    calls = [subprocess.Popen(["bash", "experiments/builds/heart_r77_line.sh", head, run, name], cwd=repo,
                              env={**env, **extra}, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for run, name, extra in (("p1", "OTK05", {"SAME_AS": TABLE["OTK20"][2]}), ("p2", "OTK20", {}))]
    outs = [c.communicate(timeout=120) for c in calls]
    codes = [c.returncode for c in calls]
    assert 34 in codes and set(codes) <= {0, 34}, (codes, outs)
