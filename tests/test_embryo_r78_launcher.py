"""The test-board build launcher of row R78's two files (the embryo's fourth refinement round on O, each O's command
with one option string of the builder's EMBRYO_R78_TABLE and --captain-ruling R78), one file per call:
experiments/builds/embryo_r78_line.sh.

It refuses any name but the two of its literal table, an odd run name or commit, a fourth argument, or a run started
after 10:00 UTC on 2026-11-26, before it makes a folder; its gate checks the commit, the clean tree, the rules file's
hash, row R78's first cell (the file's option string and `--captain-ruling R78` as whole backticked strings, the string
paired there with the file's name), the released embryo stages, the venv and the pins against their committed SHA256SUMS
lists; then it runs the targeted tests (the embryo launchers' lists and its own), rebuilds O and stops unless it gives
its recorded bytes, and in the first file's run only rebuilds row R75's Ot56 and Ot512 the same way; then it builds the
file once on its own path with its t_mix put back at 1/2 (experiments/builds/embryo_r78_check.py half) and stops unless
that gives O's bytes; then it builds the file twice from one clean commit, compares the two, and stops if the file has
O's, Ot23's or one of row R75's files' bytes. These tests hold its table to the builder's, its hashes to the committed
lists, its commands to the ones the builder takes and to the earlier launchers', its clock to the deadline, check that
it deletes nothing and reads no output of a test build but exit codes, sha256 and the stamps' commit and clean fields,
and run it for real in throwaway repositories under a temporary home: with a stand-in `uv` that fails (every refusal
path, traced with bash -x, and a call passing the gate stops at the tests step and builds nothing) and with stand-ins
that pass (a whole run of each name goes through every step in order; a stand-in file of other bytes stops the run where
its step says). Runs use a copy whose deadline lies far ahead, so that they do not depend on the day they run; the real
deadline is checked as written and by a copy whose deadline has passed. Nothing here reads the repository's own ledger
or rules file but one check, by code, that its row R78 passes the launcher's own gate lines; build_test is imported only
to check its table and that it takes the commands the launcher runs, against a temporary ledger."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import pytest
from conftest import git, git_env
from test_embryo_r75_build import ROW_R75
from test_embryo_r78_build import ARGS, NAMES, OPTS, ROW_R78, TABLE

REPO = Path(__file__).resolve().parents[1]
REL = "experiments/builds/embryo_r78_line.sh"
LAUNCHER = REPO / REL
RESULTS = "experiments/builds/results"
FLOOR = "--method fgw_shape0 --positions floor_exact"
OV_TAIL = "--shape-interp overlay --shape-interp-test-ok"
DEADLINE = int(datetime(2026, 11, 26, 10, 0, tzinfo=UTC).timestamp())
REAL_CLOCK = f"NOT_AFTER={DEADLINE}  #"
FAR = 4102444800  # 2100-01-01 00:00 UTC: the run copies' deadline
# every pin: variable prefix -> (its file's stem, sha256, its committed list, its entry there)
PINS = {
    "O": ("T2embryo_s2x_ov", "1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea",
          "embov_3b813b5.SHA256SUMS", "OVa/T2embryo_s2x_ov.h5ad"),
    "OT23": ("T2embryo_s2x_ovt23", "753900fab42738a8695d8048645cccd953ded07cd15539101964c26259aad22f",
             "embot23_045a92e_o2.SHA256SUMS", "OTa/T2embryo_s2x_ovt23.h5ad"),
    "OT712": ("T2embryo_s2x_ovt712", "ab14288faad3ad10c990073bb370c25dbb2cb3258d605cea28d5774d2dd35a42",
              "embryor75_ee22b75_r712.SHA256SUMS", "ovt712a/T2embryo_s2x_ovt712.h5ad"),
    "OT56": ("T2embryo_s2x_ovt56", "c58fe13bbca497fd8b7f11abab89094b6583b4fc6beebe0fa1a2c1db386d213c",
             "embryor75_ee22b75_r56.SHA256SUMS", "ovt56a/T2embryo_s2x_ovt56.h5ad"),
    "OT512": ("T2embryo_s2x_ovt512", "598db8bb7cd125dabe45af22da7496a0d20b072f967a4ad4c87dc26bbeced7b5",
              "embryor75_ee22b75_r512.SHA256SUMS", "ovt512a/T2embryo_s2x_ovt512.h5ad"),
    "SOT23": ("T2embryo_s2xs_ovt23", "baaeff58340a0ab60fa921f11fdead7effd1e9df85bf61330caec673cb257d06",
              "embryor75_ee22b75_rs23.SHA256SUMS", "sovt23a/T2embryo_s2xs_ovt23.h5ad"),
    "OVD": ("T2embryo_s2x_ovd", "96a13734d618eea3329c0a843c6e0e5b9ec2369d17a6c37aff2f7994733fbda8",
            "embryor75_ee22b75_rd.SHA256SUMS", "ovda/T2embryo_s2x_ovd.h5ad"),
}
# identity step -> (pin prefix, args variable, the args, exit code)
KNOWN = {
    "O": ("O", "O_ARGS", f"{FLOOR} --expression mixture {OV_TAIL}", 21),
    "ovt56": ("OT56", "OT56_ARGS", f"{FLOOR} {OV_TAIL} --mix-lean later2 --captain-ruling R75", 25),
    "ovt512": ("OT512", "OT512_ARGS", f"{FLOOR} {OV_TAIL} --mix-lean earlier_half --captain-ruling R75", 26),
}
STEMS = {"T2embryo_s2x_ovt1112.h5ad": "ovt1112", "T2embryo_s2x_ovt13.h5ad": "ovt13"}
R75_RUN = {"T2embryo_s2x_ovt1112.h5ad": "yes", "T2embryo_s2x_ovt13.h5ad": "no"}  # Ot56 and Ot512 rebuilt: one run
ROWS = f"{ROW_R75}\n{ROW_R78}\n"
FIRST = "T2embryo_s2x_ovt13.h5ad"  # the file the gate's refusals are run on
NEAR_MISS_NAMES = ["T2embryo_s2x_ovt1112", "T2embryo_s2x_ovt13", "t2embryo_s2x_ovt13.h5ad", "T2embryo_s2x_ovt13.H5AD",
                   "T2embryo_s2x_ovt13.h5ad ", " T2embryo_s2x_ovt13.h5ad", "T2embryo_s2x_ovt2224.h5ad",
                   "T2embryo_s2x_ovt26.h5ad", "T2embryo_s2xs_ovt13.h5ad", "T2embryo_s2xs_ovt56.h5ad",
                   "T2embryo_s2x_ovt56.h5ad", "T2embryo_s2x_ovt512.h5ad", "T2embryo_s2x_ovt712.h5ad",
                   "T2embryo_s2x_ovt23.h5ad", "T2embryo_s2x_ov.h5ad", "T2embryo_s2x_ovt11.h5ad",
                   "T2embryo_s2x_ovt1.h5ad", "T2embryo_s2x_ovt13.h5ad.h5ad", "runs/T2embryo_s2x_ovt13.h5ad",
                   "ovt1112", "ovt13", "*", "T2embryo_s2x_ovt1*", "$(id)"]
RULES = "the amendment names the rows\n"


def _sh() -> str:
    return LAUNCHER.read_text()


def _runnable(sh: str | None = None, not_after: int = FAR) -> str:
    """The launcher (or the given copy) with its deadline moved, the one change a run copy makes."""
    sh = _sh() if sh is None else sh
    assert sh.count(REAL_CLOCK) == 1
    return sh.replace(REAL_CLOCK, f"NOT_AFTER={not_after}  #")


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _table(sh: str) -> dict[str, tuple[str, str, str]]:
    """The launcher's FILE_NAME -> (option string, STEM, R75) table (its case lines)."""
    rows = re.findall(r'^  (T2embryo_[A-Za-z0-9_]+\.h5ad)\) OPT="([^"]+)"; STEM=([a-z0-9]+); R75=(yes|no) ;;$', sh,
                      flags=re.M)
    return {name: (opt, stem, r75) for name, opt, stem, r75 in rows}


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (REPO / RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


def _steps(name: str) -> list[str]:
    stem = STEMS[name]
    r75 = ["ovt56", "ovt512"] if R75_RUN[name] == "yes" else []
    return ["gate", "tests", "O", *r75, f"{stem}half", f"{stem}a", f"{stem}b", "done"]


def _pytest_list(text: str) -> list[str]:
    """The test files of a launcher's one pytest call: from `-m pytest` to the first line that no backslash
    continues."""
    lines = text[text.index("-m pytest"):].splitlines()
    end = next(i for i, line in enumerate(lines) if not line.rstrip().endswith("\\"))
    return re.findall(r"tests/test_[a-z0-9_]+\.py", "\n".join(lines[:end + 1]))


# ----------------------------------------------------------------------------------------------- what the launcher runs
def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launchers_literal_table_is_the_builders_one_table():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    table = _table(sh)
    assert [(f, v[0]) for f, v in table.items()] == [tuple(p) for p in bt.EMBRYO_R78_TABLE] == TABLE
    assert {f: v[1] for f, v in table.items()} == STEMS and {f: v[2] for f, v in table.items()} == R75_RUN
    assert len(set(STEMS.values())) == 2 and not set(STEMS.values()) & set(KNOWN)
    assert list(R75_RUN.values()).count("yes") == 1 and R75_RUN[NAMES[0]] == "yes"  # the row's first file's run only
    args = 'ARGS="--method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok $OPT ' \
           '--captain-ruling R78"'
    assert args in sh and bt.EMBRYO_R78_COMMAND == "embryo " + args[6:-1].replace("$OPT", "{option}")
    assert '  *) echo "refused: FILE_NAME is not one of row R78\'s two files" >&2; exit 89 ;;' in sh
    assert "NAME=${FILE%.h5ad}" in sh
    # a literal table: no name is read by a pattern (no regex, glob or slice of FILE_NAME but the suffix)
    case = sh[sh.index('case "$FILE" in'):sh.index("esac", sh.index('case "$FILE" in'))]
    assert [line.split(")")[0].strip() for line in case.splitlines()[1:] if not line.strip().startswith("*)")] == NAMES
    assert "BASH_REMATCH" not in sh and "=~" not in sh and sh.count('case "$FILE" in') == 1


def test_the_pins_are_the_committed_lists_and_the_args_the_earlier_launchers():
    sh = _sh()
    for prefix, (stem, sha, sums, entry) in PINS.items():
        assert f"{prefix}_SHA={sha}" in sh and sh.count(sha) == 1, prefix  # the pin, written once
        assert _sums(sums)[entry] == sha and entry.endswith(f"/{stem}.h5ad"), prefix  # the committed list's entry
        assert f"{prefix}_SUMS=$R/{sums}" in sh and f"{prefix}_ENTRY={entry}\n" in sh, prefix
    assert "R=experiments/builds/results\n" in sh
    assert 'KNOWN="O_SHA OT23_SHA OT712_SHA OT56_SHA OT512_SHA SOT23_SHA OVD_SHA"' in sh
    assert set(re.findall(r"^([A-Z0-9]+)_SHA=", sh, flags=re.M)) == set(PINS)
    assert not re.search(r"results/[a-z0-9_]+\.txt", sh)  # the pins come from SHA256SUMS lists, never a .txt record
    # each list's twin build agrees with its pin
    assert _sums(PINS["O"][2])["OVb/T2embryo_s2x_ov.h5ad"] == PINS["O"][1]
    for prefix in ("OT712", "OT56", "OT512", "SOT23", "OVD"):
        stem, sha, sums, entry = PINS[prefix]
        assert _sums(sums)[entry.replace("a/", "b/", 1)] == sha, prefix
    args = dict(re.findall(r'^([A-Z0-9]+_ARGS)="([^"]+)"$', sh, flags=re.M))
    assert args == {v[1]: v[2] for v in KNOWN.values()}
    assert f'OV_ARGS="{KNOWN["O"][2]}"' in (REPO / "experiments" / "builds" / "embryo_overlay.sh").read_text()
    r75 = (REPO / "experiments" / "builds" / "embryo_r75_line.sh").read_text()
    r75_args = re.search(r'^ARGS="([^"]+)"$', r75, flags=re.M).group(1)
    for step, opt in (("ovt56", "--mix-lean later2"), ("ovt512", "--mix-lean earlier_half")):
        assert KNOWN[step][2] == r75_args.replace("$OPT", opt), step  # row R75's own command for the file
        assert f'OPT="{opt}"; STEM={step}; BASES="half:O" ;;' in r75


def test_the_steps_run_in_order_and_each_build_is_checked_and_stamped():
    sh = _sh()
    for step, (prefix, key, _, code) in KNOWN.items():
        assert f'identity {step} {PINS[prefix][0]} "${prefix}_SHA" {code} ${key}\n' in sh, step
    marks = ("STEP=gate", "STEP=tests", "identity O ", 'if [ "$R75" = yes ]; then', "identity ovt56 ",
             "identity ovt512 ", "STEP=${STEM}half",
             'for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$NAME" $ARGS || finish $?; step_ok; done',
             "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    r75 = sh[sh.index('if [ "$R75" = yes ]; then'):sh.index("STEP=${STEM}half")]
    assert r75.count("identity ") == 2 and r75.rstrip().endswith("fi")
    half = sh[sh.index("STEP=${STEM}half"):sh.index('for STEP in "${STEM}a"')]
    assert half.splitlines()[1] == "begin"
    assert 'build_base "$STEP" "${NAME}_half" half $ARGS || finish $?' in half
    assert 'same "$STEP" "${NAME}_half" "$O_SHA" 51 T2embryo_s2x_ov' in half and "step_ok" in half
    build_base = sh[sh.index("build_base() {"):sh.index("same() {")]
    assert ('$PY experiments/builds/embryo_r78_check.py "$mode" embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 '
            "|| return 12") in build_base and 'checked "$d" "$name"' in build_base
    build = sh[sh.index("build() {"):sh.index("build_base() {")]
    assert ('$PY experiments/builds/build_test.py embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11'
            in build) and 'checked "$d" "$name"' in build
    done = sh[sh.index("STEP=done"):]
    assert '(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40' in done
    assert 'cmp -s "$O/${STEM}a/$FILE" "$O/${STEM}b/$FILE" || finish 41' in done
    assert 'new=$(sha256sum "$O/${STEM}a/$FILE" | cut -d\' \' -f1)' in done
    assert 'for var in $KNOWN; do\n  [ "$new" != "${!var}" ] || {' in done and "finish 42" in done
    assert done.index("finish 41") < done.index("finish 42") < done.index("step_ok")
    assert 'got=$(sha256sum "$O/$1/$2.h5ad" | cut -d\' \' -f1)' in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("identity O ")
    assert "export NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    first_uv = sh.index("uv run")
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < first_uv and ".light" not in sh
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert "--seed" not in sh and "--board T2" not in sh and "vec check" in sh and "BOARD=T2:embryo:test_interp" in sh
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")
    assert "O=runs/embryor78_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/embryor78/$RUN" in sh


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for line in lines:  # every builder, checker build and check call: output to /dev/null
        if re.search(r'build_test\.py embryo|embryo_r78_check\.py "|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    code = _code(sh)
    for word in (".build.json", "du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -",
                 "tail ", "cat ", "python3 -c", "anndata"):
        assert word not in code, word
    assert "tests.out" in sh  # the tests' log: synthetic stages only
    assert 'print(s["commit"], s["clean"])' in sh and "stamp.json" in sh  # of the stamp, its commit and clean only
    assert "curl" not in sh and "wget" not in sh and "git fetch" not in sh and "git pull" not in sh


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    for step in ("tests", "${STEM}half", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    ident = sh[sh.index("identity() {"):sh.index("identity O ")]
    assert ident.index("STEP=$1") < ident.index("begin") < ident.index("build ") < ident.index("step_ok")
    folder = sh.index('mkdir -p "$HOME/v/scratch/embryor78"')
    for code in ("exit 90", "exit 89", "exit 87"):  # every refusal of the arguments and the clock before any folder
        assert sh.rindex(code) < folder, code
    assert sh.index("exit 89 ;;") < sh.index("exit 87") < folder < sh.index("exit 91") < sh.index("STEP=gate")
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def test_the_deadline_is_ten_utc_on_the_day_read_from_bashs_own_clock():
    sh = _sh()
    assert DEADLINE == 1795687200 and sh.count(REAL_CLOCK) == 1
    assert "printf -v now '%(%s)T' -1" in sh  # bash's builtin clock: no `date` on PATH can stand in for it
    refusal = '[ "$now" -le "$NOT_AFTER" ] || { echo "refused: no run starts after 10:00 UTC on 2026-11-26" >&2; '
    assert refusal + "exit 87; }" in sh
    assert not re.search(r"(^|[\s;&|($])date(\s|$)", _code(sh), flags=re.M)
    assert "export PATH" not in sh[:sh.index("printf -v now")]  # the clock is read before PATH changes


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    assert '[ "$(grep -c "^| R78 " experiments/harness/TRIES.md)" = 1 ]' in gate
    assert "r78=$(grep \"^| R78 \" experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert "grep -q \"^ R78 (the captain's order\"" in gate and 'grep -qF -- "ruling, no run"' in gate
    assert 'for opt in "$OPT" "--captain-ruling R78"; do' in gate
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in gate  # whole backticked strings, line by line
    assert "grep -qF -- \"\\`$OPT\\` ($FILE)\"" in gate  # the string paired with the file's name
    assert 'for stage in E7.5 E8.0; do' in gate and '[ -f "data/raw/T2/$stage.h5ad" ]' in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    assert "for var in $KNOWN; do" in gate and "awk -v e=\"${!entry}\" '$2 == e {print $1}' \"${!sums}\"" in gate
    assert '[ "$got" = "${!var}" ] || { echo "refused: ${!sums} does not give $var for ${!entry}" >&2; finish 85; }' \
        in gate
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 95, 95, 95, 92, 88, 85, 93]


def test_the_gate_lines_pass_the_repositorys_row_r78_for_each_file():
    """The launcher's own row lines (bash, grep) on the repository's ledger, printing nothing: each file passes."""
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    lines = gate[gate.index('[ "$(grep -c "^| R78 "'):gate.index("for stage in E7.5 E8.0")]
    for name, option in TABLE:
        script = f'set -u\nfinish() {{ exit "$1"; }}\nOPT="{option}"\nFILE="{name}"\n{lines}\nexit 0\n'
        out = subprocess.run(["bash", "-c", script], cwd=REPO, capture_output=True, text=True)
        assert out.returncode == 0, name  # nothing of the row printed


def test_the_tests_step_runs_the_embryo_launchers_lists_and_its_own():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = _pytest_list(sh)
    assert {"tests/test_embryo_r78_launcher.py", "tests/test_embryo_r78_build.py"} <= set(listed)
    for rel in ("experiments/builds/embryo_r75_line.sh", "experiments/builds/embryo_ot23.sh",
                "experiments/builds/embryo_bf.sh", "experiments/coords/run_e58.sh"):
        theirs = _pytest_list((REPO / rel).read_text())
        assert theirs and set(theirs) <= set(listed), rel
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step and '--basetemp="$S/tmp/pytest"' in step and "-p no:cacheprovider" in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = _code(_sh()).splitlines()
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    assert "rm " not in _sh()
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 3
    for msg in usages:  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg
    assert not re.search(r"\$\{[A-Za-z_][A-Za-z0-9_]*:\?[^}]*'", _sh())  # nor in any other ${var:?message}


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", sh)}
    used |= {int(c) for c in re.findall(r'identity \S+ \S+ "\$\w+" (\d+)', sh)}
    used |= {int(c) for c in re.findall(r'same "\$STEP" "\$\{NAME\}_half" "\$O_SHA" (\d+)', sh)}
    documented = {int(c) for c in re.findall(r"(?:^|[\s,#])(\d\d) ", doc, flags=re.M)}
    assert used == documented, (used, documented)


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt
    import embryo_r78_check as chk

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with rows R75's and R78's first cells, not the repository's
    ledger.write_text(f"# ledger\n{ROWS}")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)

    def parse(args: str, out: str = "x.h5ad"):
        a = bt.build_parser().parse_args(["embryo", *args.split(), "--out", out])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        return a

    for name, (opt, stem, _) in _table(_sh()).items():
        a = parse(f"{FLOOR} {OV_TAIL} {opt} --captain-ruling R78", f"runs/x/{stem}a/{name}")  # as the a build's
        assert bt._r78_refusals(a)["file"] == name
        assert " ".join(ARGS[name]) == f"{FLOOR} {OV_TAIL} {opt} --captain-ruling R78"
    assert chk.MODES == ("half",)
    a = parse(KNOWN["O"][2])
    assert bt._r78_refusals(a) is None and bt._r75_refusals(a) is None  # O goes its own route, as before
    for step, file in (("ovt56", "T2embryo_s2x_ovt56.h5ad"), ("ovt512", "T2embryo_s2x_ovt512.h5ad")):
        a = parse(KNOWN[step][2])
        assert bt._r78_refusals(a) is None and bt._r75_refusals(a)["file"] == file, step  # row R75's route


# ----------------------------------------------------------------------------------------------- the launcher, for real
def _env(tmp_path: Path) -> dict:
    env = git_env()
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    env["UV_PROJECT_ENVIRONMENT"] = str(tmp_path / "elsewhere")  # the launcher unsets it
    return env


def _repo(tmp_path: Path, tries: str = ROWS, stages: tuple[str, ...] = ("E7.5", "E8.0"), venv: str | None = "own",
          good_sha: bool = True, records: dict | None = None, launcher: str | None = None) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (a run copy, or the given text) and the files its gate reads (the
    committed SHA256SUMS lists of every pin, or the given texts of them; stand-in released stages; a .venv whose vec.pth
    names its src, or another checkout's src, or none), and a home whose ~/.local/bin holds a `uv` that fails and lists
    its calls (the launcher puts that folder first on PATH): returns the repository, the env and its HEAD."""
    env = _env(tmp_path)
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text("#!/bin/sh\necho \"$*\" >> \"$HOME/uvcalls\"\nexit 1\n")
    (fake / "uv").chmod(0o755)
    repo = tmp_path / "repo"
    for d in ("experiments/builds/results", "experiments/harness", "docs", "src"):
        (repo / d).mkdir(parents=True, exist_ok=True)
    (repo / "src" / "keep.txt").write_text("src\n")
    (repo / REL).write_text(_runnable() if launcher is None else launcher)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / "docs" / "S2_RULES.md").write_text(RULES)
    sha = hashlib.sha256(RULES.encode()).hexdigest() if good_sha else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    for _, _, sums, _ in PINS.values():
        text = (records or {}).get(sums, (REPO / RESULTS / sums).read_text())
        if text is not None:
            (repo / RESULTS / sums).write_text(text)
    (repo / ".gitignore").write_text("data/\nruns/\n.venv/\n")
    for stage in stages:
        (repo / "data" / "raw" / "T2").mkdir(parents=True, exist_ok=True)
        (repo / "data" / "raw" / "T2" / f"{stage}.h5ad").write_text("a stand-in, never read\n")
    if venv is not None:
        site = repo / ".venv" / "lib" / "python3.12" / "site-packages"
        site.mkdir(parents=True)
        src = repo / "src" if venv == "own" else tmp_path / "other" / "src"
        (site / "vec.pth").write_text(f"{src}\n")
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "throwaway")
    return repo, env, git(repo, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str, trace: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", *(["-x"] if trace else []), REL, *argv], cwd=repo, env=env, capture_output=True,
                          text=True, timeout=120)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "embryor78" / run


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no head file, no later step's mark
    assert not list((s / "tmp").iterdir())
    assert not (tmp_path / "home" / "uvcalls").exists()  # nothing ran through uv


@pytest.mark.parametrize("name", NAMES)
def test_each_name_passes_the_gate_on_the_row_and_stops_at_the_failing_tests_step(tmp_path, name):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1", name)
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"embryor78_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built
    calls = (tmp_path / "home" / "uvcalls").read_text().splitlines()
    assert len(calls) == 1 and "-m pytest" in calls[0] and "tests/test_embryo_r78_launcher.py" in calls[0]


def test_a_traced_refusal_and_a_traced_gate_run_nothing_but_git_and_the_shell(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "t1", FIRST, trace=True)
    assert out.returncode == 96 and "HEAD is not" in out.stderr and "finish 96" in out.stderr
    assert not re.search(r"^\++ (uv|python3?|date) ", out.stderr, flags=re.M)
    out = _run(repo, env, head, "t2", FIRST, trace=True)
    assert out.returncode == 10 and "finish 10" in out.stderr
    traced = re.findall(r"^\++ (uv|python3?) .*$", out.stderr, flags=re.M)
    assert traced == ["uv"]  # the tests step's one call, which the stand-in fails
    assert not re.search(r"build_test\.py embryo|embryo_r78_check\.py", out.stderr)


def test_a_run_started_after_the_deadline_is_refused_before_any_folder(tmp_path):
    now = int(time.time())
    repo, env, head = _repo(tmp_path, launcher=_runnable(not_after=now - 3600))
    out = _run(repo, env, head, "late", FIRST, trace=True)
    assert out.returncode == 87 and "no run starts after 10:00 UTC on 2026-11-26" in out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()
    assert not (tmp_path / "home" / "uvcalls").exists()
    assert not re.search(r"^\++ date ", out.stderr, flags=re.M)  # bash's own clock, no command
    # a deadline an hour ahead lets the same call through to the gate (here, past it to the failing tests step)
    repo2, env2, head2 = _repo(tmp_path / "early", launcher=_runnable(not_after=now + 3600))
    assert _run(repo2, env2, head2, "early", FIRST).returncode == 10


def test_a_stand_in_date_on_the_path_does_not_move_the_clock(tmp_path):
    """A `date` in ~/.local/bin (first on PATH once the launcher sets it) cannot make a late run look early."""
    now = int(time.time())
    repo, env, head = _repo(tmp_path, launcher=_runnable(not_after=now - 3600))
    fake = tmp_path / "home" / ".local" / "bin" / "date"
    fake.write_text("#!/bin/sh\necho 1000000000\n")
    fake.chmod(0o755)
    env = {**env, "PATH": f"{fake.parent}:{env['PATH']}"}
    out = _run(repo, env, head, "late", FIRST)
    assert out.returncode == 87 and not (tmp_path / "home" / "v").exists()


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", FIRST).returncode == 10
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", NAMES[0])
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", FIRST)
    assert other.returncode == 10 and (repo / "runs" / f"embryor78_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_is_refused_at_the_gate_and_an_odd_one_before_any_folder(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1", FIRST)
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()
    for odd in (head[:12], head.upper(), head + "0", f"{head[:39]}g"):
        out = _run(repo, env, odd, "r2", FIRST)
        assert out.returncode == 90 and "COMMIT must be a full hash" in out.stderr, odd
        assert not _scratch(tmp_path, "r2").exists()
    out = _run(repo, env, "", "r2", FIRST)
    assert out.returncode != 0 and not _scratch(tmp_path, "r2").exists()


@pytest.mark.parametrize("argv, code", [
    ([], None), (["r1"], None), (["r1", ""], None),  # missing arguments (an empty name is none)
    (["a/b", FIRST], 90), (["r 1", FIRST], 90), (["$(id)", FIRST], 90), (["r1;x", FIRST], 90),
    (["r1", FIRST, "extra"], 89), (["r1", FIRST, ""], 89),
    *[(["r1", bad], 89) for bad in NEAR_MISS_NAMES],
])
def test_odd_arguments_and_unknown_names_are_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, *argv)
    assert out.returncode != 0 if code is None else out.returncode == code, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROWS}more\n")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


FIRST_OPT = OPTS[FIRST]


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": f"{ROW_R75}\n"}, 95, "exactly one row R78"),
    ({"tries": ROWS + f"{ROW_R78}\n"}, 95, "exactly one row R78"),  # two rows with one ID
    ({"tries": ROWS.replace("| R78 (", "| R78b (")}, 95, "exactly one row R78"),  # another ID
    ({"tries": ROWS.replace("R78 (the captain's order", "R78 (a proposal")}, 95, "not the captain's order"),
    ({"tries": ROWS.replace("R78 (the captain's order", "R78 (the captain's ruling")}, 95, "not the captain's order"),
    ({"tries": ROWS.replace("R78 (the captain's", "R78 awaits (the captain's")}, 95, "not the captain's order"),
    ({"tries": ROWS.replace("ruling, no run: \"Yes to T2embryo_s2x_ovt1112", "a try: \"Yes to T2embryo_s2x_ovt1112")},
     95, "does not say ruling, no run"),
    # the file's string inside a longer backticked string, or with no backticks, is not named
    ({"tries": ROWS.replace(f"`{FIRST_OPT}`", f"`{FIRST_OPT} --herd`")}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROWS.replace(f"`{FIRST_OPT}`", FIRST_OPT)}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROWS.replace("`--captain-ruling R78`", "`--captain-ruling R75`")}, 95,
     "does not name `--captain-ruling R78`"),
    ({"tries": ROWS.replace("`--captain-ruling R78`", "`--captain-ruling R78 now`")}, 95,
     "does not name `--captain-ruling R78`"),
    # the string named, but paired with another file's name, or not followed by " (FILE_NAME)"
    ({"tries": ROWS.replace(f"`{FIRST_OPT}` ({FIRST})", f"`{FIRST_OPT}` (T2embryo_s2x_ovt1112.h5ad)")}, 95,
     f"does not pair `{FIRST_OPT}` with ({FIRST})"),
    ({"tries": ROWS.replace(f"`{FIRST_OPT}` ({FIRST})", f"`{FIRST_OPT}`, ({FIRST})")}, 95,
     f"does not pair `{FIRST_OPT}` with ({FIRST})"),
    ({"tries": (f"| R78 (the captain's order of X; ruling, no run) | T2 embryo | `{FIRST_OPT}` ({FIRST}) "
                "`--captain-ruling R78` |\n")}, 95, f"does not name `{FIRST_OPT}`"),  # named outside the first cell
    ({"stages": ()}, 92, "data/raw/T2/E7.5.h5ad is not in this checkout"),
    ({"stages": ("E7.5",)}, 92, "data/raw/T2/E8.0.h5ad is not in this checkout"),
    ({"stages": ("E8.0",)}, 92, "data/raw/T2/E7.5.h5ad is not in this checkout"),
    ({"venv": None}, 88, "no .venv of its own"),
    ({"venv": "other"}, 88, "no .venv of its own"),  # a venv whose vec.pth names another checkout's src
    *[({"records": {PINS[k][2]: (REPO / RESULTS / PINS[k][2]).read_text().replace(PINS[k][1], "0" * 64)}}, 85,
       f"does not give {k}_SHA for {PINS[k][3]}") for k in PINS],
    *[({"records": {PINS[k][2]: None}}, 85, f"does not give {k}_SHA") for k in PINS],  # no list at all
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"embryor78_{head[:7]}_r1").exists()


@pytest.mark.parametrize("name", NAMES)
def test_the_gate_takes_a_short_row_in_its_own_words(tmp_path, name):
    tries = (f"| R78 (the captain's order of X, ruling, no run: `{OPTS[name]}` ({name}), `--captain-ruling R78`) "
             "| T2 embryo |\n")
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, head, "r1", name)
    assert out.returncode == 10, out.stderr  # past the gate: the stand-in uv fails the tests step
    assert (_scratch(tmp_path, "r1") / "gate.exit").read_text().strip() == "0"
    other = [f for f in NAMES if f != name][0]  # a file the short row does not name
    out = _run(repo, env, head, "r2", other)
    assert out.returncode == 95 and f"does not name `{OPTS[other]}`" in out.stderr


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"embryor78_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


# A stand-in `uv` for runs past the tests step: the tests pass; a build writes "file <its output's name>" with a stamp
# naming HEAD and clean (the checker writes O's stand-in, or BASE_GIVES's); FILE_GIVES sets what both builds of the file
# write; with B_DIFFERS set the second build of the file writes other bytes; `vec check` passes; every builder and
# checker call is listed in calls.txt; the stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec "$REAL_PY" -I -B -c "$@" ;;
      experiments/builds/build_test.py|experiments/builds/embryo_r78_check.py)
        prog=$1; shift
        echo "$prog $*" >> calls.txt
        mode=""
        [ "$prog" = experiments/builds/embryo_r78_check.py ] && { mode=$1; shift; }
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        if [ -n "$mode" ]; then
          name=T2embryo_s2x_ov
          [ -n "${BASE_GIVES:-}" ] && name=$BASE_GIVES
        fi
        case "$out" in
          *a/T2embryo_s2x_ovt1112.h5ad|*b/T2embryo_s2x_ovt1112.h5ad|*a/T2embryo_s2x_ovt13.h5ad|*b/T2embryo_s2x_ovt13.h5ad)
            [ -n "${FILE_GIVES:-}" ] && name=$FILE_GIVES ;;
        esac
        case "$out" in *b/*) [ -n "${B_DIFFERS:-}" ] && name="$name other" ;; esac
        echo "file $name" > "$out"
        printf '{"commit": "%s", "clean": true}' "$(git rev-parse HEAD)" > "$out.stamp.json"
        exit 0 ;;
    esac ;;
esac
exit 1
"""


def _stand_in_sha(name: str) -> str:
    return hashlib.sha256(f"file {name}\n".encode()).hexdigest()


def _stand_in_run(tmp_path: Path, name: str, keep: tuple[str, ...] = (), **env_extra):
    """A run copy of the launcher and the committed lists with every pin set to its stand-in file's (those named in
    keep left as they are), the stand-in uv above, one call for `name`."""
    sh, records = _runnable(), {}
    for prefix, (stem, sha, sums, entry) in PINS.items():
        if prefix in keep:
            continue
        stand_in = _stand_in_sha(stem)
        sh = sh.replace(f"{prefix}_SHA={sha}", f"{prefix}_SHA={stand_in}")
        text = records.get(sums, (REPO / RESULTS / sums).read_text())
        assert text.count(f"{sha}  {entry}\n") == 1
        records[sums] = text.replace(f"{sha}  {entry}\n", f"{stand_in}  {entry}\n")
    assert sh != _runnable()
    repo, env, head = _repo(tmp_path, launcher=sh, records=records)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(FULL_UV)
    env = {**env, "REAL_PY": sys.executable, **env_extra}
    return repo, head, _run(repo, env, head, "r1", name)


@pytest.mark.parametrize("name", NAMES)
def test_a_whole_stand_in_run_passes_every_step_in_order(tmp_path, name):
    """With only its pins set to the stand-in files' (in a copy, with the lists to match) and the stand-in `uv` above:
    a call runs every step in order to done.exit 0, each step's exit file says 0, the stamps list every build, the
    builder and the checker get the commands their steps name (row R75's two rebuilds in the first file's run only),
    the hash list holds every built file once and the two builds of the file are one file."""
    repo, head, out = _stand_in_run(tmp_path, name)
    assert out.returncode == 0, out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(name)
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {step: "0" for step in steps}
    r = f"runs/embryor78_{head[:7]}_r1"
    o = repo / r
    folders = steps[2:-1]
    assert sorted(p.name for p in o.iterdir()) == sorted([*folders, "SHA256SUMS"])
    stamps = (s / "stamps").read_text().splitlines()
    assert [line.split()[0].split("/")[0] for line in stamps] == folders
    assert all(line.endswith(f" {head} True") for line in stamps)
    listed = (o / "SHA256SUMS").read_text().splitlines()
    assert [line.split()[1].split("/")[0] for line in listed] == folders and (s / "SHA256SUMS").exists()
    stem, base = STEMS[name], name.removesuffix(".h5ad")
    sums = dict(line.split()[::-1] for line in listed)
    assert sums[f"{stem}half/{base}_half.h5ad"] == _stand_in_sha("T2embryo_s2x_ov")
    args = " ".join(ARGS[name])
    rebuilt = [k for k in KNOWN if k in steps]
    want = [*(f"experiments/builds/build_test.py embryo {KNOWN[k][2]} --out {r}/{k}/{PINS[KNOWN[k][0]][0]}.h5ad"
              for k in rebuilt),
            f"experiments/builds/embryo_r78_check.py half embryo {args} --out {r}/{stem}half/{base}_half.h5ad",
            f"experiments/builds/build_test.py embryo {args} --out {r}/{stem}a/{name}",
            f"experiments/builds/build_test.py embryo {args} --out {r}/{stem}b/{name}"]
    assert (repo / "calls.txt").read_text().splitlines() == want
    assert rebuilt == (["O", "ovt56", "ovt512"] if R75_RUN[name] == "yes" else ["O"])
    assert (o / f"{stem}a" / name).read_bytes() == (o / f"{stem}b" / name).read_bytes()


@pytest.mark.parametrize("name, step", [(n, k) for n in NAMES for k in KNOWN if k in _steps(n)])
def test_an_identity_file_of_other_bytes_stops_the_run_at_its_step(tmp_path, name, step):
    repo, head, out = _stand_in_run(tmp_path, name, (KNOWN[step][0],))
    code = KNOWN[step][3]
    assert out.returncode == code, out.stderr
    assert f"stopped: step {step} does not give the known bytes of {PINS[KNOWN[step][0]][0]}" in out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(name)
    ran = steps[:steps.index(step) + 1]
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {**{x: "0" for x in ran[:-1]},
                                                                              step: str(code)}
    assert not (repo / "runs" / f"embryor78_{head[:7]}_r1" / f"{STEMS[name]}a").exists()


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("gives", ["T2embryo_s2x_ovt56", "T2embryo_s2x_ovt23", "own"])
def test_a_half_check_of_other_bytes_stops_the_run_before_the_file_is_built(tmp_path, name, gives):
    gives = name.removesuffix(".h5ad") if gives == "own" else gives
    repo, head, out = _stand_in_run(tmp_path, name, BASE_GIVES=gives)
    stem = STEMS[name]
    assert out.returncode == 51, out.stderr
    assert f"stopped: step {stem}half does not give the known bytes of T2embryo_s2x_ov" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / f"{stem}half.exit").read_text().strip() == "51" and not (s / f"{stem}a.exit").exists()
    assert not (repo / "runs" / f"embryor78_{head[:7]}_r1" / f"{stem}a").exists()


def test_two_builds_of_other_bytes_stop_the_run_at_done(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, FIRST, B_DIFFERS="1")
    assert out.returncode == 41, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "done.exit").read_text().strip() == "41" and (s / "SHA256SUMS").exists()
    o = repo / "runs" / f"embryor78_{head[:7]}_r1"
    assert (o / "ovt13a" / FIRST).exists() and (o / "ovt13b" / FIRST).exists()  # both kept: nothing deleted


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("prefix", list(PINS))
def test_a_file_with_a_known_files_bytes_stops_the_run_at_done(tmp_path, name, prefix):
    """Both builds of the file agree, but with the bytes of O, Ot23 or one of row R75's files: the run stops at done
    (42), the hash list written and both builds kept."""
    repo, head, out = _stand_in_run(tmp_path, name, FILE_GIVES=PINS[prefix][0])
    assert out.returncode == 42, out.stderr
    assert f"stopped: {name} has the bytes of a known file ({prefix}_SHA)" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "done.exit").read_text().strip() == "42" and (s / "SHA256SUMS").exists()
    o = repo / "runs" / f"embryor78_{head[:7]}_r1"
    assert (o / f"{STEMS[name]}a" / name).exists() and (o / f"{STEMS[name]}b" / name).exists()


def test_a_build_that_fails_or_misses_its_stamp_stops_its_step(tmp_path):
    """A build that exits non-zero stops its step with 11 (the checker's with 12); a stamp of another commit with 15."""
    name = FIRST
    repo, head, _ = _stand_in_run(tmp_path, name)
    env = {**_env(tmp_path), "REAL_PY": sys.executable}
    uv = tmp_path / "home" / ".local" / "bin" / "uv"
    failing_build = FULL_UV.replace('        prog=$1; shift\n', '        prog=$1; shift\n'
                                    '        [ "$prog" = experiments/builds/build_test.py ] && exit 1\n')
    uv.write_text(failing_build)
    out = _run(repo, env, head, "b1", name)
    assert out.returncode == 11 and (_scratch(tmp_path, "b1") / "O.exit").read_text().strip() == "11"
    failing_check = FULL_UV.replace('        prog=$1; shift\n', '        prog=$1; shift\n'
                                    '        [ "$prog" = experiments/builds/embryo_r78_check.py ] && exit 1\n')
    uv.write_text(failing_check)
    out = _run(repo, env, head, "b2", name)
    assert out.returncode == 12 and (_scratch(tmp_path, "b2") / "ovt13half.exit").read_text().strip() == "12"
    other_commit = FULL_UV.replace('"$(git rev-parse HEAD)"', '"0000000000000000000000000000000000000000"')
    uv.write_text(other_commit)
    out = _run(repo, env, head, "b3", name)
    assert out.returncode == 15 and (_scratch(tmp_path, "b3") / "O.exit").read_text().strip() == "15"
