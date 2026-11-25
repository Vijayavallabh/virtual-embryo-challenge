"""The test-board build launcher of row R73's four files (T1_s2rcMSBm10's command with one further change, g_mom's
damping or Blood's step), one file per call.

experiments/builds/t1_r73_line.sh refuses any name but the four of its literal table before it makes a folder (the
near-miss names T1_s2rcMSBd10m10, T1_s2rcMSBm10m10 and T1_s2rcMSBm10d10H05 among them); its gate checks the commit,
the clean tree, the rules file's hash, row R73's first cell (the file's option string and `--captain-ruling R73` as
whole backticked strings, the string paired there with the file's name), the relabelling's output and the venv; then it
runs the targeted tests, rebuilds T1_s2rcMSB, T1_s2rcMSBm10 and row R51's single-step file at the file's value and
stops unless each gives its known bytes; then it builds the file once on its own path with its one further step put
back at T1_s2rcMSBm10's value (experiments/builds/t1_r73_check.py base) and stops unless that is T1_s2rcMSBm10's bytes;
then it builds the file twice from one clean commit and compares the two. These tests hold its table to the builder's,
its hashes to the committed lists, its commands to the ones the builder takes and to the earlier launchers', check that
it deletes nothing and reads no output of a test build but exit codes, sha256 and the stamps' commit and clean fields,
and run it for real in throwaway repositories under a temporary home: with a stand-in `uv` that fails (every refusal
path, and a call passing the gate stops at the tests step and builds nothing) and with stand-ins that pass (a whole run
of each of the four names goes through every step in order; a stand-in file of other bytes stops the run where its
step says). Nothing here reads the repository's own ledger or rules file; build_test is imported only to check its
table and that it takes the commands the launcher runs, against a temporary ledger."""

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
import test_t1_r51_launcher as r51l
from test_t1_l67pb_launcher import PASSING_UV, _sums
from test_t1_lm_m20_build import ROW_R58
from test_t1_r73_build import NEAR_MISS_NAMES, OPTS, ROW_R73, SINGLE

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_r73_line.sh"
R67_LINE = REPO / "experiments" / "builds" / "t1_r67_line.sh"
R66_LINE = REPO / "experiments" / "builds" / "t1_r66_line.sh"
R51_LINE = REPO / "experiments" / "builds" / "t1_r51_line.sh"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
MSB_OPT = "--feng-fraction-mult 1.5 --direction-momentum all"
# identity step -> (file name, args variable, the args, hash variable, sha256, exit code, the earlier launcher that
# builds it with the same args); MSB and m10 for every file, then the one single-step file at the file's value
KNOWN = {
    "MSB": ("T1_s2rcMSB", "MSB", f"{BASE} {MSB_OPT} --captain-ruling R26", "MSB_SHA",
            "78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2", 23, R67_LINE),
    "m10": ("T1_s2rcMSBm10", "M10", f"{BASE} {MSB_OPT} --direction-momentum-damping 1.0 --captain-ruling R51",
            "M10_SHA", "a0c46efac87deb3756b3ca489a215694640707ea8a8053e254b8a86ee9d7d964", 24, R67_LINE),
    "d10": ("T1_s2rcMSBd10", "D10", f"{BASE} {MSB_OPT} --g-mom-damping 1.0 --captain-ruling R51", "D10_SHA",
            "835e68fb9c390e6271e0ea8ca2ddfb5431f2f6d333b739fe17ef0102a0989ae8", 71, R51_LINE),
    "d15": ("T1_s2rcMSBd15", "D15", f"{BASE} {MSB_OPT} --g-mom-damping 1.5 --captain-ruling R51", "D15_SHA",
            "7d2745d884dba02d94db6dc01133d9971c52aa63fe5191aca60473a40b7a8c75", 72, R51_LINE),
    "H05": ("T1_s2rcMSBH05", "H05", f"{BASE} {MSB_OPT} --held-momentum 0.5 --captain-ruling R51", "H05_SHA",
            "f8bbec6990354d42ed1b754d795a9ff617b2ee6c40a2d5376a274486e2f77c07", 74, R51_LINE),
    "H10": ("T1_s2rcMSBH10", "H10", f"{BASE} {MSB_OPT} --held-momentum 1.0 --captain-ruling R51", "H10_SHA",
            "1cfdfc54d838e03294c444dd284c319724519eb731349667e5bf7f4e79f79250", 75, R51_LINE),
}
NAMES = list(OPTS)
ONE = {name: SINGLE[name].removeprefix("T1_s2rcMSB") for name in NAMES}  # d10, H05, d15, H10
ROW = f"{ROW_R73}\n"
FIRST = "T1_s2rcMSBm10H05"  # the file the gate's refusals are run on
FIRST_OPT = OPTS[FIRST]


def _sh() -> str:
    return LAUNCHER.read_text()


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _table(sh: str) -> dict[str, tuple[str, str]]:
    """The launcher's FILE_NAME -> (option string, ONE) table (its case lines), $MSB_OPT written out."""
    rows = re.findall(r'^  (T1_s2rcMSBm10[A-Za-z0-9]+)\) OPT="\$MSB_OPT ([^"]+)"; ONE=([A-Za-z0-9]+) ;;$', sh,
                      flags=re.M)
    return {name: (f"{MSB_OPT} {rest}", one) for name, rest, one in rows}


def _stem(name: str) -> str:
    return name.removeprefix("T1_s2rcMSB")


def _steps(name: str) -> list[str]:
    stem = _stem(name)
    return ["gate", "tests", "MSB", "m10", ONE[name], f"{stem}base", f"{stem}a", f"{stem}b", "done"]


def _listed_tests(path: Path) -> list[str]:
    sh = path.read_text()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    return re.findall(r"tests/test_[a-z0-9_]+\.py", step)


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launchers_literal_table_is_the_builders_one_table():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    table = _table(sh)
    assert [(f, v[0]) for f, v in table.items()] == [(f, v[0]) for f, v in bt.T1_R73_FILES.items()]
    assert [tuple(p) for p in bt.T1_R73_TABLE] == [(f, v[0]) for f, v in table.items()] == list(OPTS.items())
    assert {f: v[1] for f, v in table.items()} == ONE
    assert {f: f"T1_s2rcMSB{one}" for f, one in ONE.items()} == bt.T1_R73_SINGLE
    assert f'MSB_OPT="{MSB_OPT}"' in sh
    assert 'ARGS="--g-mom --base lib3 --feng-groups recomputed $OPT --captain-ruling R73"' in sh
    assert '  *) echo "refused: FILE_NAME is not one of row R73\'s four files" >&2; exit 89 ;;' in sh
    assert sh.count('OPT="$MSB_OPT ') == 4 and "STEM=${FILE#T1_s2rcMSB}" in sh
    # a literal table: no name is read by a pattern (no regex, glob or slice of FILE_NAME but the stem)
    case = sh[sh.index('case "$FILE" in'):sh.index("esac", sh.index('case "$FILE" in'))]
    assert all(re.fullmatch(r"T1_s2rcMSBm10(d10|d15|H05|H10)", line.split(")")[0].strip())
               for line in case.splitlines()[1:] if not line.strip().startswith("*)"))
    assert "BASH_REMATCH" not in sh and "=~" not in sh
    assert {_stem(n) for n in NAMES} == {"m10d10", "m10H05", "m10d15", "m10H10"}
    assert not {_stem(n) for n in NAMES} & set(KNOWN) and set(ONE.values()) <= set(KNOWN)


def test_the_known_hashes_are_the_committed_lists_and_the_args_the_earlier_launchers():
    sh = _sh()
    msc = _sums("t1msc_a3e879a_m1.SHA256SUMS")
    r51 = _sums("t1r51_0cdcc6f_r1.SHA256SUMS")
    assert msc["MSBa/T1_s2rcMSB.h5ad"] == msc["MSBb/T1_s2rcMSB.h5ad"] == r51["MSB/T1_s2rcMSB.h5ad"] == KNOWN["MSB"][4]
    for step in ("m10", "d10", "d15", "H05", "H10"):
        name, sha = KNOWN[step][0], KNOWN[step][4]
        assert r51[f"{step}a/{name}.h5ad"] == r51[f"{step}b/{name}.h5ad"] == sha, step
        assert r51[f"{step}msb/{name}_msb.h5ad"] == KNOWN["MSB"][4], step  # each one's own path at MSB's value
    args = dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert set(args) == {v[1] for v in KNOWN.values()}
    for step, (_, key, full, var, sha, _, earlier) in KNOWN.items():
        assert f"{var}={sha}" in sh and sh.count(sha) == 1, step
        assert args[key] == full, step
        assert f'{key}_ARGS="{full}"' in earlier.read_text(), step  # the earlier launcher's own command


def test_the_steps_run_in_order_and_each_build_is_checked_and_stamped():
    sh = _sh()
    for step in ("MSB", "m10"):
        name, key, _, var, _, code, _ = KNOWN[step]
        assert f'identity {step} {name} "${var}" {code} ${key}_ARGS\n' in sh
    singles = sh[sh.index('case "$ONE" in'):sh.index("esac", sh.index('case "$ONE" in'))]
    for step in ("d10", "d15", "H05", "H10"):
        name, key, _, var, _, code, _ = KNOWN[step]
        assert f'  {step}) identity {step} {name} "${var}" {code} ${key}_ARGS ;;\n' in singles
    assert singles.count("identity ") == 4
    marks = ("STEP=gate", "STEP=tests", "identity MSB ", "identity m10 ", 'case "$ONE" in', "STEP=${STEM}base",
             'for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$FILE" $ARGS || finish $?; step_ok; done',
             "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    base = sh[sh.index("STEP=${STEM}base"):sh.index('for STEP in "${STEM}a"')]
    assert 'build_base "$STEP" "${FILE}_base" $ARGS || finish $?' in base
    assert 'same "$STEP" "${FILE}_base" "$M10_SHA" 51' in base
    build_base = sh[sh.index("build_base() {"):sh.index("same() {")]
    assert ('$PY experiments/builds/t1_r73_check.py base t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12'
            in build_base) and 'checked "$d" "$name"' in build_base
    build = sh[sh.index("build() {"):sh.index("build_base() {")]
    assert ('$PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11'
            in build) and 'checked "$d" "$name"' in build
    done = sh[sh.index("STEP=done"):]
    assert '(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40' in done
    assert 'cmp -s "$O/${STEM}a/$FILE.h5ad" "$O/${STEM}b/$FILE.h5ad" || finish 41' in done
    assert 'got=$(sha256sum "$O/$1/$2.h5ad" | cut -d\' \' -f1)' in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("identity MSB ")
    assert "export NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert "export PATH=$HOME/.local/bin:$PATH" in sh
    first_uv = sh.index("uv run")
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < first_uv and ".light" not in sh
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert "--seed" not in sh and "--board T1" not in sh and "vec check" in sh and "BOARD=T1:test" in sh
    assert "upload" not in sh.replace("Nothing here uploads", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")
    assert "O=runs/t1r73_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1r73/$RUN" in sh


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for line in lines:  # every builder, checker build and check call: output to /dev/null
        if re.search(r'build_test\.py t1|t1_r73_check\.py base|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    code = _code(sh)
    for word in (".build.json", "du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -",
                 "tail ", "cat "):
        assert word not in code, word
    assert "tests.out" in sh  # the tests' log: synthetic stages only
    assert 'print(s["commit"], s["clean"])' in sh and "stamp.json" in sh  # of the stamp, its commit and clean only
    assert "curl" not in sh and "wget" not in sh and "git fetch" not in sh and "git pull" not in sh


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    for step in ("tests", "${STEM}base", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    ident = sh[sh.index("identity() {"):sh.index("identity MSB ")]
    assert ident.index("STEP=$1") < ident.index("begin") < ident.index("build ") < ident.index("step_ok")
    for code in ("exit 90", "exit 89", "exit 91"):
        assert sh.index(code) < sh.index("STEP=gate"), code
    assert sh.index("exit 89 ;;") < sh.index('mkdir -p "$HOME/v/scratch/t1r73"')  # every refusal before any folder
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
    assert 'r73=$(first_cell R73 95) || exit $?' in gate
    assert 'for opt in "$OPT" "--captain-ruling R73"; do' in gate
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in gate  # whole backticked strings, line by line
    assert "grep -qF -- \"\\`$OPT\\` ($FILE)\"" in gate  # the string paired with the file's name
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 92, 88, 93]


def test_the_tests_step_runs_the_new_tests_and_r67s_and_r66s():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = _listed_tests(LAUNCHER)
    assert {"tests/test_t1_r73_launcher.py", "tests/test_t1_r73_build.py", *_listed_tests(R67_LINE),
            *_listed_tests(R66_LINE)} == set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step and '--basetemp="$S/tmp/pytest"' in step


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
    used |= {int(c) for c in re.findall(r'same "\$STEP" \S+ "\$M10_SHA" (\d+)', sh)}
    documented = {int(c) for c in re.findall(r"(?:^|[\s,#])(\d\d) ", doc, flags=re.M)}
    assert used == documented, (used, documented)


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    sys.path.insert(0, str(REPO / "experiments" / "external"))
    import build_test as bt
    import feng_relabel as fr

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with the ruling rows' orders, not the repository's
    ledger.write_text(r51l.ROWS_ALL + r51l.ROW_R51 + f"{ROW_R58}\n{ROW}")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())

    def parse(args: str):
        a = bt.build_parser().parse_args(["t1", *args.split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        assert bt._t1_r54(a) is None and bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None
        return a

    a = parse(KNOWN["MSB"][2])  # nothing read but the ledger
    assert bt._msb_ordered(a) and bt._t1_r58(a) is None and bt._t1_r67(a) is None and bt._t1_r51(a) is None
    for step in ("m10", "d10", "d15", "H05", "H10"):
        a = parse(KNOWN[step][2])
        assert bt._t1_l12(a) is None and bt._t1_r51(a)["file"] == KNOWN[step][0], step
    for name, (opt, _) in _table(_sh()).items():
        a = parse(f"{BASE} {opt} --captain-ruling R73")
        got = bt._t1_r51(a)
        assert got == bt._t1_l12(a) and got["file"] == name and got["option"] == f"`{opt}`"
        assert got["captain_ruling"] == "R73" and got["base_file"] == "T1_s2rcMSBm10"
        assert bt._t1_r58(a) is None and bt._t1_r67(a) is None and bt._t1_pooled(a) is None and bt._l12_ordered(a)


# ------------------------------------------------------------------ the launcher, run for real


def _repo(tmp_path: Path, tries: str = ROW, launcher: str | None = None, **kw):
    """tests/test_t1_r51_launcher.py's throwaway repository (its launcher path holds this launcher, or the given copy
    of it) with a stand-in `uv` that fails."""
    return r51l._repo(tmp_path, tries=tries, launcher=_sh() if launcher is None else launcher, **kw)


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return r51l._run(repo, env, *argv)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1r73" / run


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
    o = repo / "runs" / f"t1r73_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", FIRST).returncode == 10
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "T1_s2rcMSBm10d15")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", FIRST)
    assert other.returncode == 10 and (repo / "runs" / f"t1r73_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_and_a_prefix_are_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1", FIRST)
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()
    short = _run(repo, env, head[:12], "r2", FIRST)
    assert short.returncode == 96
    _stopped_at_the_gate(tmp_path, "r2", 96)


@pytest.mark.parametrize("argv, code", [
    ([], None), (["r1"], None), (["r1", ""], None),  # missing arguments (an empty name is none)
    (["a/b", FIRST], 90), (["r 1", FIRST], 90), (["$(id)", FIRST], 90),
    (["r1", FIRST, "extra"], 89),
    *[(["r1", bad], 89) for bad in (
        *NEAR_MISS_NAMES,  # T1_s2rcMSBd10m10 (row R66's order), T1_s2rcMSBm10m10, T1_s2rcMSBm10d10H05
        "T1_s2rcMSBm10", "T1_s2rcMSB", "T1_s2rcMSBd10", "T1_s2rcMSBH05", "T1_s2rcMSBm10d20", "T1_s2rcMSBm10H20",
        "T1_s2rcMSBm10d05", "T1_s2rcMSBm10H5", "T1_s2rcMSBm10d1.0", "T1_s2rcMSBm10d10 ", " T1_s2rcMSBm10d10",
        "t1_s2rcMSBm10d10", "T1_s2rcMSBM10d10", "T1_s2rcMS20d10m10", "T1_s2rcMS20m10", "T1_s2rcMSBm125",
        "T1_s2rcMSBm10d10.h5ad", "m10d10", "d10", "$(id)", "*", "T1_s2rcMSBm10*")],
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


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ""}, 95, "exactly one row R73"),
    ({"tries": ROW + ROW}, 95, "exactly one row R73"),  # two rows with one ID
    ({"tries": ROW.replace("| R73 (", "| R73b (")}, 95, "exactly one row R73"),  # another ID
    ({"tries": ROW.replace("R73 (the captain's order", "R73 (a proposal")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("R73 (the captain's", "R73 awaits (the captain's")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("ruling, no run", "a try")}, 95, "does not say ruling, no run"),
    # the file's string inside a longer backticked string, or with no backticks, is not named
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", f"`{FIRST_OPT} --herd`")}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace(f"`{FIRST_OPT}`", FIRST_OPT)}, 95, f"does not name `{FIRST_OPT}`"),
    ({"tries": ROW.replace("`--captain-ruling R73`", "`--captain-ruling R67`")}, 95,
     "does not name `--captain-ruling R73`"),
    ({"tries": ROW.replace("`--captain-ruling R73`", "`--captain-ruling R73 now`")}, 95,
     "does not name `--captain-ruling R73`"),
    # the string named, but paired with another file's name
    ({"tries": ROW.replace(f"`{FIRST_OPT}` ({FIRST})", f"`{FIRST_OPT}` (T1_s2rcMSBm10H10)")}, 95,
     f"does not pair `{FIRST_OPT}` with ({FIRST})"),
    ({"tries": (f"| R73 (the captain's order of X; ruling, no run) | T1 | `{FIRST_OPT}` ({FIRST}) "
                "`--captain-ruling R73` |\n")}, 95, f"does not name `{FIRST_OPT}`"),  # named outside the first cell
    ({"relabel": False}, 92, "feng_relabel_c is not in this checkout"),
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
    assert not (repo / "runs" / f"t1r73_{head[:7]}_r1").exists()


@pytest.mark.parametrize("name", NAMES)
def test_the_gate_takes_a_short_row_in_its_own_words(tmp_path, name):
    tries = (f"| R73 (the captain's order of X, ruling, no run: `{OPTS[name]}` ({name}), `--captain-ruling R73`) "
             "| T1 |\n")
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, head, "r1", name)
    assert out.returncode == 10, out.stderr  # past the gate: the stand-in uv fails the tests step
    assert (_scratch(tmp_path, "r1") / "gate.exit").read_text().strip() == "0"
    other = [f for f in NAMES if f != name][0]  # a file the short row does not name
    out = _run(repo, env, head, "r2", other)
    assert out.returncode == 95 and f"does not name `{OPTS[other]}`" in out.stderr


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t1r73_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_past_the_tests_step_stops_at_the_first_unknown_bytes(tmp_path):
    """tests/test_t1_l67pb_launcher.py's stand-in `uv`: the tests pass, a build writes a stand-in file whose stamp
    names HEAD and clean, `vec check` passes; no stand-in file has a known file's bytes, so the run stops at MSB."""
    repo, env, head = _repo(tmp_path)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(PASSING_UV)
    out = _run(repo, env, head, "r1", FIRST)
    assert out.returncode == 23, out.stderr
    assert "stopped: step MSB does not give the known bytes of T1_s2rcMSB" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert {p.name: p.read_text().strip() for p in s.glob("*.exit")} == {"gate.exit": "0", "tests.exit": "0",
                                                                          "MSB.exit": "23"}
    assert (s / "stamps").read_text() == f"MSB/T1_s2rcMSB {head} True\n"
    o = repo / "runs" / f"t1r73_{head[:7]}_r1"
    assert sorted(p.name for p in o.iterdir()) == ["MSB"]


# A stand-in `uv` for whole runs: the tests pass; a build writes "file <its output's name>" with a stamp naming HEAD
# and clean (the checker's base writes T1_s2rcMSBm10's stand-in, as each file's path with its one further step at
# T1_s2rcMSBm10's value must, or BASE_GIVES's); with B_DIFFERS set the second build of the file writes other bytes;
# `vec check` passes; every builder and checker call is listed in calls.txt; the stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/build_test.py|experiments/builds/t1_r73_check.py)
        prog=$1; shift
        echo "$prog $*" >> calls.txt
        [ "$prog" = experiments/builds/t1_r73_check.py ] && shift
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        [ "$prog" = experiments/builds/t1_r73_check.py ] && name=${BASE_GIVES:-T1_s2rcMSBm10}
        case "$out" in *b/*) [ -n "${B_DIFFERS:-}" ] && name="$name other" ;; esac
        echo "file $name" > "$out"
        printf '{"commit": "%s", "clean": true}' "$(git rev-parse HEAD)" > "$out.stamp.json"
        exit 0 ;;
    esac ;;
esac
exit 1
"""


def _stand_in(sh: str, keep: tuple[str, ...] = ()) -> str:
    """The launcher with its known hashes set to the stand-in files' (those of the variables in keep left as they
    are)."""
    for name, _, _, var, sha, _, _ in KNOWN.values():
        if var not in keep:
            stand_in = hashlib.sha256(f"file {name}\n".encode()).hexdigest()
            sh = sh.replace(f"{var}={sha}", f"{var}={stand_in}")
    return sh


def _stand_in_run(tmp_path: Path, name: str, keep: tuple[str, ...] = (), **env_extra):
    sh = _stand_in(_sh(), keep)
    assert sh != _sh()
    repo, env, head = _repo(tmp_path, launcher=sh)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(FULL_UV)
    return repo, head, _run(repo, {**env, **env_extra}, head, "r1", name)


@pytest.mark.parametrize("name", NAMES)
def test_a_whole_stand_in_run_passes_every_step_in_order(tmp_path, name):
    """With only its known hashes set to the stand-in files' (in a copy) and the stand-in `uv` above: a call runs every
    step in order to done.exit 0, each step's exit file says 0, the stamps list every build, the builder and the
    checker get the commands their steps name, the hash list holds every built file once and the two builds of the
    file are one file."""
    repo, head, out = _stand_in_run(tmp_path, name)
    assert out.returncode == 0, out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(name)
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {step: "0" for step in steps}
    r = f"runs/t1r73_{head[:7]}_r1"
    o = repo / r
    folders = steps[2:-1]
    assert sorted(p.name for p in o.iterdir()) == sorted([*folders, "SHA256SUMS"])
    stamps = (s / "stamps").read_text().splitlines()
    assert [line.split()[0].split("/")[0] for line in stamps] == folders
    assert all(line.endswith(f" {head} True") for line in stamps)
    listed = (o / "SHA256SUMS").read_text().splitlines()
    assert [line.split()[1].split("/")[0] for line in listed] == folders and (s / "SHA256SUMS").exists()
    stem, m10 = _stem(name), hashlib.sha256(b"file T1_s2rcMSBm10\n").hexdigest()
    assert dict(line.split()[::-1] for line in listed)[f"{stem}base/{name}_base.h5ad"] == m10
    args = f"{BASE} {OPTS[name]} --captain-ruling R73"
    want = [*(f"experiments/builds/build_test.py t1 {KNOWN[step][2]} --out {r}/{step}/{KNOWN[step][0]}.h5ad"
              for step in ("MSB", "m10", ONE[name])),
            f"experiments/builds/t1_r73_check.py base t1 {args} --out {r}/{stem}base/{name}_base.h5ad",
            f"experiments/builds/build_test.py t1 {args} --out {r}/{stem}a/{name}.h5ad",
            f"experiments/builds/build_test.py t1 {args} --out {r}/{stem}b/{name}.h5ad"]
    assert (repo / "calls.txt").read_text().splitlines() == want
    assert (o / f"{stem}a" / f"{name}.h5ad").read_bytes() == (o / f"{stem}b" / f"{name}.h5ad").read_bytes()


@pytest.mark.parametrize("name, step", [(n, s) for n in NAMES for s in ("MSB", "m10", ONE[n])])
def test_an_identity_file_of_other_bytes_stops_the_run_at_its_step(tmp_path, name, step):
    repo, head, out = _stand_in_run(tmp_path, name, (KNOWN[step][3],))
    code = KNOWN[step][5]
    assert out.returncode == code, out.stderr
    assert f"stopped: step {step} does not give the known bytes of {KNOWN[step][0]}" in out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(name)
    ran = steps[:steps.index(step) + 1]
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {
        **{x: "0" for x in ran[:-1]}, step: str(code)}
    assert not (repo / "runs" / f"t1r73_{head[:7]}_r1" / f"{_stem(name)}a").exists()


@pytest.mark.parametrize("gives", [FIRST, "T1_s2rcMSB", "T1_s2rcMSBH05"])
def test_a_base_check_of_other_bytes_stops_the_run_before_the_file_is_built(tmp_path, gives):
    repo, head, out = _stand_in_run(tmp_path, FIRST, BASE_GIVES=gives)
    assert out.returncode == 51, out.stderr
    stem = _stem(FIRST)
    assert f"stopped: step {stem}base does not give the known bytes of {FIRST}_base" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / f"{stem}base.exit").read_text().strip() == "51" and not (s / f"{stem}a.exit").exists()
    assert not (repo / "runs" / f"t1r73_{head[:7]}_r1" / f"{stem}a").exists()


def test_two_builds_of_other_bytes_stop_the_run_at_done(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, FIRST, B_DIFFERS="1")
    assert out.returncode == 41, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "done.exit").read_text().strip() == "41" and (s / "SHA256SUMS").exists()
