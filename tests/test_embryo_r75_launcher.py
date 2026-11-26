"""The test-board build launcher of row R75's five files (the embryo's refinement round on O, each O's command with one
option string of the builder's EMBRYO_R75_TABLE and --captain-ruling R75), one file per call:
experiments/builds/embryo_r75_line.sh.

It refuses any name but the five of its literal table, an odd run name or commit, or a fourth argument, before it makes
a folder; its gate checks the commit, the clean tree, the rules file's hash, row R75's first cell (the file's option
string and `--captain-ruling R75` as whole backticked strings, the string paired there with the file's name), the
released embryo stages, the venv and the identity pins against their committed SHA256SUMS lists; then it runs the
targeted tests, rebuilds O, row R59's Ot23, SO and row R48's W and stops unless each gives its recorded bytes; then it
builds the file once per base check on its own path with its one change put back
(experiments/builds/embryo_r75_check.py: a lean's t_mix at 1/2, an expression at O's mixture) and stops unless that
gives the known file's bytes; then it builds the file twice from one clean commit and compares the two. These tests
hold its table to the builder's, its hashes to the committed lists, its commands to the ones the builder takes and to
the earlier launchers', check that it deletes nothing and reads no output of a test build but exit codes, sha256 and
the stamps' commit and clean fields, and run it for real in throwaway repositories under a temporary home: with a
stand-in `uv` that fails (every refusal path, traced with bash -x, and a call passing the gate stops at the tests step
and builds nothing) and with stand-ins that pass (a whole run of each of the five names goes through every step in
order; a stand-in file of other bytes stops the run where its step says). Nothing here reads the repository's own
ledger or rules file but one check, by code, that its row R75 passes the launcher's own gate lines; build_test is
imported only to check its table and that it takes the commands the launcher runs, against a temporary ledger."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import git, git_env
from test_embryo_r75_build import ARGS, NAMES, OPTS, ROW_R75, T_MIX, TABLE

REPO = Path(__file__).resolve().parents[1]
REL = "experiments/builds/embryo_r75_line.sh"
LAUNCHER = REPO / REL
RESULTS = "experiments/builds/results"
FLOOR = "--method fgw_shape0 --positions floor_exact"
OV_TAIL = "--shape-interp overlay --shape-interp-test-ok"
# identity step -> (file name, args variable, the args, sha variable, sha256, exit code, its committed list and entry,
# the earlier launcher whose own variable holds the same args)
KNOWN = {
    "O": ("T2embryo_s2x_ov", "O_ARGS", f"{FLOOR} --expression mixture {OV_TAIL}", "O_SHA",
          "1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea", 21,
          "embov_3b813b5.SHA256SUMS", "OVa/T2embryo_s2x_ov.h5ad", ("embryo_overlay.sh", "OV_ARGS")),
    "OT": ("T2embryo_s2x_ovt23", "OT_ARGS",
           f"{FLOOR} --expression mixture {OV_TAIL} --mix-lean later --captain-ruling R59",
           "OT_SHA", "753900fab42738a8695d8048645cccd953ded07cd15539101964c26259aad22f", 22,
           "embot23_045a92e_o2.SHA256SUMS", "OTa/T2embryo_s2x_ovt23.h5ad", ("embryo_ot23.sh", "OT_ARGS")),
    "SO": ("T2embryo_s2xs_ov", "SO_ARGS", f"{FLOOR} --expression series {OV_TAIL}", "SO_SHA",
           "77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165", 23,
           "embsohi_93f832e_r1.SHA256SUMS", "SOa/T2embryo_s2xs_ov.h5ad", ("embryo_so_ohi.sh", "SO_ARGS")),
    "W": ("T2embryo_s2x_ovsb", "W_ARGS", f"{FLOOR} --expression series_back {OV_TAIL} --captain-ruling R48", "W_SHA",
          "f576b64b8b93d79b1269b580df31bd5c6e58246cd4f3eda12cba7d77387087ab", 24,
          "embovsb_20b2381_w1.SHA256SUMS", "Wa/T2embryo_s2x_ovsb.h5ad", ("embryo_ovsb.sh", "W_ARGS")),
}
STEMS = {"T2embryo_s2x_ovt712.h5ad": "ovt712", "T2embryo_s2x_ovt56.h5ad": "ovt56", "T2embryo_s2x_ovt512.h5ad": "ovt512",
         "T2embryo_s2xs_ovt23.h5ad": "sovt23", "T2embryo_s2x_ovd.h5ad": "ovd"}
# each file's base checks in order: (mode, the identity step whose bytes the file's path must then give)
BASES = {"T2embryo_s2x_ovt712.h5ad": [("half", "O")], "T2embryo_s2x_ovt56.h5ad": [("half", "O")],
         "T2embryo_s2x_ovt512.h5ad": [("half", "O")], "T2embryo_s2xs_ovt23.h5ad": [("half", "SO"), ("mixture", "OT")],
         "T2embryo_s2x_ovd.h5ad": [("mixture", "O")]}
BASE_CODES = {"half": 51, "mixture": 52}
ROW = f"{ROW_R75}\n"
FIRST = "T2embryo_s2x_ovt712.h5ad"  # the file the gate's refusals are run on
NEAR_MISS_NAMES = ["T2embryo_s2x_ovt13.h5ad", "T2embryo_s2x_ovt23.h5ad", "T2embryo_s2x_ov.h5ad",
                   "T2embryo_s2xs_ov.h5ad",
                   "T2embryo_s2x_ovsb.h5ad", "T2embryo_s2x_ovhi.h5ad", "T2embryo_s2xd_ov.h5ad", "T2embryo_s2x_ovt712",
                   "t2embryo_s2x_ovt712.h5ad", "T2embryo_s2x_ovt712.h5ad ", " T2embryo_s2x_ovt712.h5ad",
                   "T2embryo_s2x_ovt712.H5AD", "T2embryo_s2x_ovt1012.h5ad", "T2embryo_s2x_ovt56d.h5ad",
                   "T2embryo_s2xs_ovt712.h5ad", "T2embryo_s2x_ovt23d.h5ad", "T2embryo_s2x_ovd.h5ad.h5ad",
                   "runs/T2embryo_s2x_ovd.h5ad", "ovt712", "ovd", "*", "T2embryo_s2x_ovt7*", "$(id)"]
RULES = "the amendment names the rows\n"


def _sh() -> str:
    return LAUNCHER.read_text()


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _table(sh: str) -> dict[str, tuple[str, str, list[tuple[str, str]]]]:
    """The launcher's FILE_NAME -> (option string, STEM, base checks) table (its case lines)."""
    rows = re.findall(r'^  (T2embryo_[A-Za-z0-9_]+\.h5ad)\) OPT="([^"]+)"; STEM=([a-z0-9]+); BASES="([^"]+)" ;;$', sh,
                      flags=re.M)
    return {name: (opt, stem, [tuple(b.split(":")) for b in bases.split()]) for name, opt, stem, bases in rows}


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (REPO / RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


def _steps(name: str) -> list[str]:
    stem = STEMS[name]
    return ["gate", "tests", *KNOWN, *(f"{stem}{mode}" for mode, _ in BASES[name]), f"{stem}a", f"{stem}b", "done"]


def _listed_tests(sh: str) -> list[str]:
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    return re.findall(r"tests/test_[a-z0-9_]+\.py", step)


# ----------------------------------------------------------------------------------------------- what the launcher runs
def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launchers_literal_table_is_the_builders_one_table():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    table = _table(sh)
    assert [(f, v[0]) for f, v in table.items()] == [tuple(p) for p in bt.EMBRYO_R75_TABLE] == TABLE
    assert {f: v[1] for f, v in table.items()} == STEMS and {f: v[2] for f, v in table.items()} == BASES
    assert len(set(STEMS.values())) == 5 and not set(STEMS.values()) & set(KNOWN)
    args = 'ARGS="--method fgw_shape0 --positions floor_exact --shape-interp overlay --shape-interp-test-ok $OPT ' \
           '--captain-ruling R75"'
    assert args in sh and bt.EMBRYO_R75_COMMAND == "embryo " + args[6:-1].replace("$OPT", "{option}")
    assert '  *) echo "refused: FILE_NAME is not one of row R75\'s five files" >&2; exit 89 ;;' in sh
    assert "NAME=${FILE%.h5ad}" in sh
    # a literal table: no name is read by a pattern (no regex, glob or slice of FILE_NAME but the suffix)
    case = sh[sh.index('case "$FILE" in'):sh.index("esac", sh.index('case "$FILE" in'))]
    assert [line.split(")")[0].strip() for line in case.splitlines()[1:] if not line.strip().startswith("*)")] == NAMES
    assert "BASH_REMATCH" not in sh and "=~" not in sh and sh.count('case "$FILE" in') == 1
    # the base checks are each file's own changes: half for a lean, mixture for an expression
    for name, bases in BASES.items():
        modes = [m for m, _ in bases]
        lean, expression = T_MIX[name] is not None, "--expression" in OPTS[name]
        assert modes == (["half"] if lean else []) + (["mixture"] if expression else []), name


def test_the_known_hashes_are_the_committed_lists_and_the_args_the_earlier_launchers():
    sh = _sh()
    args = dict(re.findall(r'^([A-Z]+_ARGS)="([^"]+)"$', sh, flags=re.M))
    assert set(args) == {v[1] for v in KNOWN.values()}
    for step, (name, key, full, var, sha, _, sums, entry, (earlier, earlier_key)) in KNOWN.items():
        assert f"{var}={sha}" in sh and sh.count(sha) == 1, step  # the pin, written once
        assert _sums(sums)[entry] == sha and entry.endswith(f"/{name}.h5ad"), step  # the committed list's entry
        prefix = var.removesuffix("_SHA")
        assert f"{prefix}_SUMS={RESULTS}/{sums}" in sh and f"{prefix}_ENTRY={entry}" in sh
        assert args[key] == full, step
        assert f'{earlier_key}="{full}"' in (REPO / "experiments" / "builds" / earlier).read_text(), step
    assert not re.search(r"results/[a-z0-9_]+\.txt", sh)  # the pins come from SHA256SUMS lists, never a .txt record
    # the twin builds of the earlier runs agree with the pins
    assert _sums(KNOWN["O"][6])["OVb/T2embryo_s2x_ov.h5ad"] == KNOWN["O"][4]
    assert _sums(KNOWN["OT"][6])["OTb/T2embryo_s2x_ovt23.h5ad"] == KNOWN["OT"][4]
    assert _sums(KNOWN["SO"][6])["SOb/T2embryo_s2xs_ov.h5ad"] == KNOWN["SO"][4]
    assert _sums(KNOWN["W"][6])["Wb/T2embryo_s2x_ovsb.h5ad"] == KNOWN["W"][4]


def test_the_steps_run_in_order_and_each_build_is_checked_and_stamped():
    sh = _sh()
    for step, (name, key, _, var, _, code, *_) in KNOWN.items():
        assert f'identity {step} {name} "${var}" {code} ${key}\n' in sh, step
    marks = ("STEP=gate", "STEP=tests", "identity O ", "identity OT ", "identity SO ", "identity W ",
             "for base in $BASES; do", 'for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$NAME" $ARGS '
             "|| finish $?; step_ok; done", "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    loop = sh[sh.index("for base in $BASES; do"):sh.index('for STEP in "${STEM}a"')]
    assert 'build_base "$STEP" "${NAME}_$mode" "$mode" $ARGS || finish $?' in loop
    assert 'same "$STEP" "${NAME}_$mode" "${!sha}" "$code" "$known"' in loop
    assert "half) code=51 ;;" in loop and "mixture) code=52 ;;" in loop and "STEP=$STEM$mode" in loop
    build_base = sh[sh.index("build_base() {"):sh.index("same() {")]
    assert ('$PY experiments/builds/embryo_r75_check.py "$mode" embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 '
            "|| return 12") in build_base and 'checked "$d" "$name"' in build_base
    build = sh[sh.index("build() {"):sh.index("build_base() {")]
    assert ('$PY experiments/builds/build_test.py embryo "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11'
            in build) and 'checked "$d" "$name"' in build
    done = sh[sh.index("STEP=done"):]
    assert '(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40' in done
    assert 'cmp -s "$O/${STEM}a/$FILE" "$O/${STEM}b/$FILE" || finish 41' in done
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
    assert "O=runs/embryor75_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/embryor75/$RUN" in sh


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for line in lines:  # every builder, checker build and check call: output to /dev/null
        if re.search(r'build_test\.py embryo|embryo_r75_check\.py "|vec check "', line):
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
    for step in ("tests", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    loop = sh[sh.index("for base in $BASES; do"):sh.index('for STEP in "${STEM}a"')]
    assert loop.index("STEP=$STEM$mode") < loop.index("begin") < loop.index("build_base ") < loop.index("step_ok")
    ident = sh[sh.index("identity() {"):sh.index("identity O ")]
    assert ident.index("STEP=$1") < ident.index("begin") < ident.index("build ") < ident.index("step_ok")
    folder = sh.index('mkdir -p "$HOME/v/scratch/embryor75"')
    for code in ("exit 90", "exit 89"):  # every refusal of the arguments before any folder
        assert sh.rindex(code) < folder, code
    assert sh.index("exit 89 ;;") < folder < sh.index("exit 91") < sh.index("STEP=gate")
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    assert '[ "$(grep -c "^| R75 " experiments/harness/TRIES.md)" = 1 ]' in gate
    assert "r75=$(grep \"^| R75 \" experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert "grep -q \"^ R75 (the captain's order\"" in gate and 'grep -qF -- "ruling, no run"' in gate
    assert 'for opt in "$OPT" "--captain-ruling R75"; do' in gate
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in gate  # whole backticked strings, line by line
    assert "grep -qF -- \"\\`$OPT\\` ($FILE)\"" in gate  # the string paired with the file's name
    assert 'for stage in E7.5 E8.0; do' in gate and '[ -f "data/raw/T2/$stage.h5ad" ]' in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    assert '[ "$got" = "${!var}" ] || { echo "refused: $sums does not give $var for $entry" >&2; finish 85; }' in gate
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 95, 95, 95, 92, 88, 85, 93]


def test_the_gate_lines_pass_the_repositorys_row_r75_for_each_file():
    """The launcher's own row lines (bash, grep) on the repository's ledger, printing nothing: each file passes."""
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    lines = gate[gate.index('[ "$(grep -c "^| R75 "'):gate.index("for stage in E7.5 E8.0")]
    for name, option in TABLE:
        script = f'set -u\nfinish() {{ exit "$1"; }}\nOPT="{option}"\nFILE="{name}"\n{lines}\nexit 0\n'
        out = subprocess.run(["bash", "-c", script], cwd=REPO, capture_output=True, text=True)
        assert out.returncode == 0, name  # nothing of the row printed


def test_the_tests_step_runs_the_new_tests_and_the_earlier_embryo_routes():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = _listed_tests(sh)
    assert {"tests/test_embryo_r75_launcher.py", "tests/test_embryo_r75_build.py", "tests/test_mix_lean_build.py",
            "tests/test_embryo_ot23_launcher.py", "tests/test_series_back_build.py",
            "tests/test_embryo_ovsb_launcher.py", "tests/test_bf_build.py", "tests/test_s2_rules.py",
            "tests/test_hidden_counts.py", "tests/test_t2_build.py"} <= set(listed)
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
    used |= {int(c) for c in re.findall(r"(?:half|mixture)\) code=(\d+)", sh)}
    documented = {int(c) for c in re.findall(r"(?:^|[\s,#])(\d\d) ", doc, flags=re.M)}
    assert used == documented, (used, documented)


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt
    import embryo_r75_check as chk

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with row R75's first cell, not the repository's
    ledger.write_text(f"# ledger\n{ROW}")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)

    def parse(args: str):
        a = bt.build_parser().parse_args(["embryo", *args.split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        return a

    for name, (opt, _, bases) in _table(_sh()).items():
        a = parse(f"{FLOOR} {OV_TAIL} {opt} --captain-ruling R75")
        assert bt._r75_refusals(a)["file"] == name
        assert " ".join(ARGS[name]) == f"{FLOOR} {OV_TAIL} {opt} --captain-ruling R75"
        assert all(mode in chk.MODES for mode, _ in bases)
    for step, (_, _, args, *_) in KNOWN.items():  # the identity commands go their own routes, as before
        assert bt._r75_refusals(parse(args)) is None, step


# ----------------------------------------------------------------------------------------------- the launcher, for real
def _env(tmp_path: Path) -> dict:
    env = git_env()
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    env["UV_PROJECT_ENVIRONMENT"] = str(tmp_path / "elsewhere")  # the launcher unsets it
    return env


def _repo(tmp_path: Path, tries: str = ROW, stages: tuple[str, ...] = ("E7.5", "E8.0"), venv: str | None = "own",
          good_sha: bool = True, records: dict | None = None, launcher: str | None = None) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (or the given copy of it) and the files its gate reads (the four
    committed SHA256SUMS lists, or the given texts of them; stand-in released stages; a .venv whose vec.pth names its
    src, or another checkout's src, or none), and a home whose ~/.local/bin holds a `uv` that fails and lists its calls
    (the launcher puts that folder first on PATH): returns the repository, the env and its HEAD."""
    env = _env(tmp_path)
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text("#!/bin/sh\necho \"$*\" >> \"$HOME/uvcalls\"\nexit 1\n")
    (fake / "uv").chmod(0o755)
    repo = tmp_path / "repo"
    for d in ("experiments/builds/results", "experiments/harness", "docs", "src"):
        (repo / d).mkdir(parents=True, exist_ok=True)
    (repo / "src" / "keep.txt").write_text("src\n")
    (repo / REL).write_text(_sh() if launcher is None else launcher)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / "docs" / "S2_RULES.md").write_text(RULES)
    sha = hashlib.sha256(RULES.encode()).hexdigest() if good_sha else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    for _, _, _, _, _, _, sums, _, _ in KNOWN.values():
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
    return tmp_path / "home" / "v" / "scratch" / "embryor75" / run


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
    o = repo / "runs" / f"embryor75_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built
    calls = (tmp_path / "home" / "uvcalls").read_text().splitlines()
    assert len(calls) == 1 and "-m pytest" in calls[0] and "tests/test_embryo_r75_launcher.py" in calls[0]


def test_a_traced_refusal_and_a_traced_gate_run_nothing_but_git_and_the_shell(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "t1", FIRST, trace=True)
    assert out.returncode == 96 and "HEAD is not" in out.stderr and "finish 96" in out.stderr
    assert not re.search(r"^\++ (uv|python3?) ", out.stderr, flags=re.M)
    out = _run(repo, env, head, "t2", FIRST, trace=True)
    assert out.returncode == 10 and "finish 10" in out.stderr
    traced = re.findall(r"^\++ (uv|python3?) .*$", out.stderr, flags=re.M)
    assert traced == ["uv"]  # the tests step's one call, which the stand-in fails
    assert not re.search(r"build_test\.py embryo|embryo_r75_check\.py", out.stderr)


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", FIRST).returncode == 10
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "T2embryo_s2x_ovd.h5ad")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", FIRST)
    assert other.returncode == 10 and (repo / "runs" / f"embryor75_{head[:7]}_r2").is_dir()


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
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW}more\n")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


FIRST_OPT = OPTS[FIRST]


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ""}, 95, "exactly one row R75"),
    ({"tries": ROW + ROW}, 95, "exactly one row R75"),  # two rows with one ID
    ({"tries": ROW.replace("| R75 (", "| R75b (")}, 95, "exactly one row R75"),  # another ID
    ({"tries": ROW.replace("R75 (the captain's order", "R75 (a proposal")}, 95, "not the captain's order"),
    ({"tries": ROW.replace("R75 (the captain's order", "R75 (the captain's ruling")}, 95, "not the captain's order"),
    ({"tries": ROW.replace("R75 (the captain's", "R75 awaits (the captain's")}, 95, "not the captain's order"),
    ({"tries": ROW.replace("ruling, no run", "a try")}, 95, "does not say ruling, no run"),
    # the file's string inside a longer backticked string, or with no backticks, is not named
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", f"`{FIRST_OPT} --herd`")}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", FIRST_OPT)}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace("`--captain-ruling R75`", "`--captain-ruling R59`")}, 95,
     "does not name `--captain-ruling R75`"),
    ({"tries": ROW.replace("`--captain-ruling R75`", "`--captain-ruling R75 now`")}, 95,
     "does not name `--captain-ruling R75`"),
    # the string named, but paired with another file's name, or not followed by " (FILE_NAME)"
    ({"tries": ROW.replace(f"`{FIRST_OPT}` ({FIRST})", f"`{FIRST_OPT}` (T2embryo_s2x_ovt56.h5ad)")}, 95,
     f"does not pair `{FIRST_OPT}` with ({FIRST})"),
    ({"tries": ROW.replace(f"`{FIRST_OPT}` ({FIRST})", f"`{FIRST_OPT}`, ({FIRST})")}, 95,
     f"does not pair `{FIRST_OPT}` with ({FIRST})"),
    ({"tries": (f"| R75 (the captain's order of X; ruling, no run) | T2 embryo | `{FIRST_OPT}` ({FIRST}) "
                "`--captain-ruling R75` |\n")}, 95, f"does not name `{FIRST_OPT}`"),  # named outside the first cell
    ({"stages": ()}, 92, "data/raw/T2/E7.5.h5ad is not in this checkout"),
    ({"stages": ("E7.5",)}, 92, "data/raw/T2/E8.0.h5ad is not in this checkout"),
    ({"stages": ("E8.0",)}, 92, "data/raw/T2/E7.5.h5ad is not in this checkout"),
    ({"venv": None}, 88, "no .venv of its own"),
    ({"venv": "other"}, 88, "no .venv of its own"),  # a venv whose vec.pth names another checkout's src
    *[({"records": {KNOWN[k][6]: (REPO / RESULTS / KNOWN[k][6]).read_text().replace(KNOWN[k][4], "0" * 64)}}, 85,
       f"does not give {KNOWN[k][3]} for {KNOWN[k][7]}") for k in KNOWN],
    *[({"records": {KNOWN[k][6]: None}}, 85, f"does not give {KNOWN[k][3]}") for k in KNOWN],  # no list at all
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"embryor75_{head[:7]}_r1").exists()


@pytest.mark.parametrize("name", NAMES)
def test_the_gate_takes_a_short_row_in_its_own_words(tmp_path, name):
    tries = (f"| R75 (the captain's order of X, ruling, no run: `{OPTS[name]}` ({name}), `--captain-ruling R75`) "
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
    o = repo / "runs" / f"embryor75_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


# A stand-in `uv` for runs past the tests step: the tests pass; a build writes "file <its output's name>" with a stamp
# naming HEAD and clean (the checker writes the stand-in of the known file its base check must give for that file and
# mode, or BASE_GIVES's); with B_DIFFERS set the second build of the file writes other bytes; `vec check` passes; every
# builder and checker call is listed in calls.txt; the stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec "$REAL_PY" -I -B -c "$@" ;;
      experiments/builds/build_test.py|experiments/builds/embryo_r75_check.py)
        prog=$1; shift
        echo "$prog $*" >> calls.txt
        mode=""
        [ "$prog" = experiments/builds/embryo_r75_check.py ] && { mode=$1; shift; }
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        if [ -n "$mode" ]; then
          case "$mode:$name" in
            half:T2embryo_s2xs_ovt23_half) name=T2embryo_s2xs_ov ;;
            half:*) name=T2embryo_s2x_ov ;;
            mixture:T2embryo_s2xs_ovt23_mixture) name=T2embryo_s2x_ovt23 ;;
            mixture:*) name=T2embryo_s2x_ov ;;
          esac
          [ -n "${BASE_GIVES:-}" ] && name=$BASE_GIVES
        fi
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
    """The launcher and the committed lists with the identity pins set to the stand-in files' (those of the steps in
    keep left as they are), the stand-in uv above, one call for `name`."""
    sh, records = _sh(), {}
    for step, (fname, _, _, var, sha, _, sums, entry, _) in KNOWN.items():
        if step in keep:
            continue
        stand_in = _stand_in_sha(fname)
        sh = sh.replace(f"{var}={sha}", f"{var}={stand_in}")
        text = records.get(sums, (REPO / RESULTS / sums).read_text())
        assert text.count(f"{sha}  {entry}\n") == 1
        records[sums] = text.replace(f"{sha}  {entry}\n", f"{stand_in}  {entry}\n")
    assert sh != _sh()
    repo, env, head = _repo(tmp_path, launcher=sh, records=records)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(FULL_UV)
    env = {**env, "REAL_PY": sys.executable, **env_extra}
    return repo, head, _run(repo, env, head, "r1", name)


@pytest.mark.parametrize("name", NAMES)
def test_a_whole_stand_in_run_passes_every_step_in_order(tmp_path, name):
    """With only its identity pins set to the stand-in files' (in a copy, with the lists to match) and the stand-in
    `uv` above: a call runs every step in order to done.exit 0, each step's exit file says 0, the stamps list every
    build, the builder and the checker get the commands their steps name, the hash list holds every built file once
    and the two builds of the file are one file."""
    repo, head, out = _stand_in_run(tmp_path, name)
    assert out.returncode == 0, out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(name)
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {step: "0" for step in steps}
    r = f"runs/embryor75_{head[:7]}_r1"
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
    for mode, known in BASES[name]:
        assert sums[f"{stem}{mode}/{base}_{mode}.h5ad"] == _stand_in_sha(KNOWN[known][0])
    args = " ".join(ARGS[name])
    want = [*(f"experiments/builds/build_test.py embryo {KNOWN[k][2]} --out {r}/{k}/{KNOWN[k][0]}.h5ad"
              for k in KNOWN),
            *(f"experiments/builds/embryo_r75_check.py {mode} embryo {args} --out {r}/{stem}{mode}/{base}_{mode}.h5ad"
              for mode, _ in BASES[name]),
            f"experiments/builds/build_test.py embryo {args} --out {r}/{stem}a/{name}",
            f"experiments/builds/build_test.py embryo {args} --out {r}/{stem}b/{name}"]
    assert (repo / "calls.txt").read_text().splitlines() == want
    assert (o / f"{stem}a" / name).read_bytes() == (o / f"{stem}b" / name).read_bytes()


@pytest.mark.parametrize("name, step", [(n, k) for n in (FIRST, "T2embryo_s2xs_ovt23.h5ad") for k in KNOWN])
def test_an_identity_file_of_other_bytes_stops_the_run_at_its_step(tmp_path, name, step):
    repo, head, out = _stand_in_run(tmp_path, name, (step,))
    code = KNOWN[step][5]
    assert out.returncode == code, out.stderr
    assert f"stopped: step {step} does not give the known bytes of {KNOWN[step][0]}" in out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(name)
    ran = steps[:steps.index(step) + 1]
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {**{x: "0" for x in ran[:-1]},
                                                                              step: str(code)}
    assert not (repo / "runs" / f"embryor75_{head[:7]}_r1" / f"{STEMS[name]}a").exists()


@pytest.mark.parametrize("name, gives", [
    (n, g) for n in NAMES for g in [*(v[0] for v in KNOWN.values()), n.removesuffix(".h5ad")]
    if any(KNOWN[k][0] != g for _, k in BASES[n])])
def test_a_base_check_of_other_bytes_stops_the_run_before_the_file_is_built(tmp_path, name, gives):
    """A checker build giving another file's bytes (BASE_GIVES: an identity file's stand-in or the file's own) stops
    the run at the first base check whose known file it is not, before the file is built; T2embryo_s2xs_ovt23.h5ad's
    mixture check stops on SO's bytes, which its half check takes."""
    repo, head, out = _stand_in_run(tmp_path, name, BASE_GIVES=gives)
    stem = STEMS[name]
    mode, known = next((m, k) for m, k in BASES[name] if KNOWN[k][0] != gives)
    assert out.returncode == BASE_CODES[mode], out.stderr
    assert f"stopped: step {stem}{mode} does not give the known bytes of {known}" in out.stderr
    s = _scratch(tmp_path, "r1")
    passed = [m for m, _ in BASES[name]][:[m for m, _ in BASES[name]].index(mode)]
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob(f"{stem}*.exit")} == {
        **{f"{stem}{m}": "0" for m in passed}, f"{stem}{mode}": str(BASE_CODES[mode])}
    assert not (repo / "runs" / f"embryor75_{head[:7]}_r1" / f"{stem}a").exists()


def test_two_builds_of_other_bytes_stop_the_run_at_done(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, FIRST, B_DIFFERS="1")
    assert out.returncode == 41, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "done.exit").read_text().strip() == "41" and (s / "SHA256SUMS").exists()
    o = repo / "runs" / f"embryor75_{head[:7]}_r1"
    assert (o / "ovt712a" / FIRST).exists() and (o / "ovt712b" / FIRST).exists()  # both kept: nothing deleted


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
                                    '        [ "$prog" = experiments/builds/embryo_r75_check.py ] && exit 1\n')
    uv.write_text(failing_check)
    out = _run(repo, env, head, "b2", name)
    assert out.returncode == 12 and (_scratch(tmp_path, "b2") / "ovt712half.exit").read_text().strip() == "12"
    other_commit = FULL_UV.replace('"$(git rev-parse HEAD)"', '"0000000000000000000000000000000000000000"')
    uv.write_text(other_commit)
    out = _run(repo, env, head, "b3", name)
    assert out.returncode == 15 and (_scratch(tmp_path, "b3") / "O.exit").read_text().strip() == "15"
