"""The test-board build launcher of row R67's three files (T1_s2rcMSB's command with MS's own damping at 1.25, 1.5 or
0.75), one file per call.

experiments/builds/t1_r67_line.sh refuses any name but the three before it makes a folder; its gate checks the commit,
the clean tree, the rules file's hash, row R67's first cell (the file's option string and `--captain-ruling R67` as
whole backticked strings, the string paired there with the file's name), the relabelling's output and the venv; then it
runs the targeted tests, rebuilds T1_s2rcMSB, T1_s2rcMSBm10 and T1_s2rcMSBm20 and stops unless each gives its known
bytes; then it builds the file once on its own path with MS's damping put back at B0's
(experiments/builds/t1_r67_check.py base) and stops unless that is T1_s2rcMSB's bytes; then it builds the file twice
from one clean commit and compares the two. These tests hold its table to the builder's, its hashes to the committed
lists, its commands to the ones the builder takes and to the earlier launchers', check that it deletes nothing and reads
no output of a test build but exit codes, sha256 and the stamps' commit and clean fields, and run it for real in
throwaway repositories under a temporary home: with a stand-in `uv` that fails (every refusal path, and a call passing
the gate stops at the tests step and builds nothing) and with stand-ins that pass (a whole run of each of the three
names goes through every step in order; a stand-in file of other bytes stops the run where its step says). Nothing here
reads the repository's own ledger or rules file for its state of today; build_test is imported only to check its table
and that it takes the commands the launcher runs, against a temporary ledger."""

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
from test_t1_r67_build import OPTS, ROW_R67

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_r67_line.sh"
LM_M20 = REPO / "experiments" / "builds" / "t1_lm_m20_line.sh"
R51_LINE = REPO / "experiments" / "builds" / "t1_r51_line.sh"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
MSB_OPT = "--feng-fraction-mult 1.5 --direction-momentum all"
# identity step -> (file name, args variable, the args, hash variable, sha256, exit code, the earlier launcher that
# builds it with the same args)
KNOWN = {
    "MSB": ("T1_s2rcMSB", "MSB", f"{BASE} {MSB_OPT} --captain-ruling R26", "MSB_SHA",
            "78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2", 23, LM_M20),
    "m10": ("T1_s2rcMSBm10", "M10", f"{BASE} {MSB_OPT} --direction-momentum-damping 1.0 --captain-ruling R51",
            "M10_SHA", "a0c46efac87deb3756b3ca489a215694640707ea8a8053e254b8a86ee9d7d964", 24, R51_LINE),
    "m20": ("T1_s2rcMSBm20", "M20", f"{BASE} {MSB_OPT} --direction-momentum-damping 2.0 --captain-ruling R58",
            "M20_SHA", "d82c3db2e5f43a3abeefe32a04e4d58d2da41cefa8fa979eb642de884c435401", 25, LM_M20),
}
NAMES = list(OPTS)
ROW = f"{ROW_R67}\n"
ONE = "T1_s2rcMSBm125"  # the file the gate's refusals are run on
ONE_OPT = OPTS[ONE]


def _sh() -> str:
    return LAUNCHER.read_text()


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _table(sh: str) -> dict[str, str]:
    """The launcher's FILE_NAME -> option string table (its case lines), $MSB_OPT written out."""
    rows = re.findall(r'^  (T1_s2rcMSBm[0-9]+)\) OPT="\$MSB_OPT ([^"]+)" ;;$', sh, flags=re.M)
    return {name: f"{MSB_OPT} {rest}" for name, rest in rows}


def _stem(name: str) -> str:
    return name.removeprefix("T1_s2rcMSB")


def _steps(name: str) -> list[str]:
    stem = _stem(name)
    return ["gate", "tests", *KNOWN, f"{stem}base", f"{stem}a", f"{stem}b", "done"]


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launchers_table_is_the_builders_one_table():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    table = _table(sh)
    assert list(table.items()) == [(f, v[0]) for f, v in bt.T1_R67_FILES.items()] == list(OPTS.items())
    assert [tuple(p) for p in bt.T1_R67_TABLE] == list(OPTS.items())
    assert f'MSB_OPT="{MSB_OPT}"' in sh
    assert 'ARGS="--g-mom --base lib3 --feng-groups recomputed $OPT --captain-ruling R67"' in sh
    assert '*) echo "refused: FILE_NAME is not one of row R67\'s three files" >&2; exit 89 ;;' in sh
    assert sh.count('OPT="$MSB_OPT ') == 3 and "STEM=${FILE#T1_s2rcMSB}" in sh
    assert {_stem(n) for n in NAMES} == {"m125", "m15", "m075"} and not {_stem(n) for n in NAMES} & set(KNOWN)


def test_the_known_hashes_are_the_committed_lists_and_the_args_the_earlier_launchers():
    sh = _sh()
    msc = _sums("t1msc_a3e879a_m1.SHA256SUMS")
    r51 = _sums("t1r51_0cdcc6f_r1.SHA256SUMS")
    lm = _sums("t1lmm20_89700ef_lm1.SHA256SUMS")
    msb, m10, m20 = (KNOWN[k][4] for k in ("MSB", "m10", "m20"))
    assert msc["MSBa/T1_s2rcMSB.h5ad"] == msc["MSBb/T1_s2rcMSB.h5ad"] == r51["MSB/T1_s2rcMSB.h5ad"] == msb
    assert lm["MSB/T1_s2rcMSB.h5ad"] == lm["m20msb/T1_s2rcMSBm20_msb.h5ad"] == msb
    assert r51["m10msb/T1_s2rcMSBm10_msb.h5ad"] == msb
    assert r51["m10a/T1_s2rcMSBm10.h5ad"] == r51["m10b/T1_s2rcMSBm10.h5ad"] == m10
    assert lm["m20a/T1_s2rcMSBm20.h5ad"] == lm["m20b/T1_s2rcMSBm20.h5ad"] == m20
    args = dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert set(args) == {v[1] for v in KNOWN.values()}
    for step, (_, key, full, var, sha, _, earlier) in KNOWN.items():
        assert f"{var}={sha}" in sh and sh.count(sha) == 1, step
        assert args[key] == full, step
        assert f'{key}_ARGS="{full}"' in earlier.read_text(), step  # the earlier launcher's own command


def test_the_steps_run_in_order_and_each_build_is_checked_and_stamped():
    sh = _sh()
    for step, (name, key, _, var, _, code, _) in KNOWN.items():
        assert f'identity {step} {name} "${var}" {code} ${key}_ARGS\n' in sh
    marks = ("STEP=gate", "STEP=tests", "identity MSB ", "identity m10 ", "identity m20 ", "STEP=${STEM}base",
             'for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$FILE" $ARGS || finish $?; step_ok; done',
             "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    base = sh[sh.index("STEP=${STEM}base"):sh.index('for STEP in "${STEM}a"')]
    assert 'build_base "$STEP" "${FILE}_base" $ARGS || finish $?' in base
    assert 'same "$STEP" "${FILE}_base" "$MSB_SHA" 51' in base
    build_base = sh[sh.index("build_base() {"):sh.index("same() {")]
    assert ('$PY experiments/builds/t1_r67_check.py base t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12'
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
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert "export NUMPY_MADVISE_HUGEPAGE=0" in sh and "export PATH=$HOME/.local/bin:$PATH" in sh
    first_uv = sh.index("uv run")
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < first_uv and ".light" not in sh
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert "--seed" not in sh and "--board T1" not in sh and "vec check" in sh and "BOARD=T1:test" in sh
    assert "upload" not in sh.replace("Nothing here uploads", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")
    assert "O=runs/t1r67_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1r67/$RUN" in sh


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for line in lines:  # every builder, checker build and check call: output to /dev/null
        if re.search(r'build_test\.py t1|t1_r67_check\.py base|vec check "', line):
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
    assert 'r67=$(first_cell R67 95) || exit $?' in gate
    assert 'for opt in "$OPT" "--captain-ruling R67"; do' in gate
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in gate  # whole backticked strings, line by line
    assert "grep -qF -- \"\\`$OPT\\` ($FILE)\"" in gate  # the string paired with the file's name
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 92, 88, 93]


def test_the_tests_step_runs_the_new_tests_and_t1_lm_m20_line_shs():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    lm = LM_M20.read_text()
    lm_step = lm[lm.index("STEP=tests"):lm.index("step_ok", lm.index("STEP=tests"))]
    assert {"tests/test_t1_r67_launcher.py", "tests/test_t1_r67_build.py",
            *re.findall(r"tests/test_[a-z0-9_]+\.py", lm_step)} == set(listed)
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


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", sh)}
    used |= {int(c) for c in re.findall(r'identity \S+ \S+ "\$\w+" (\d+)', sh)}
    used |= {int(c) for c in re.findall(r'same "\$STEP" \S+ "\$MSB_SHA" (\d+)', sh)}
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
    a = parse(KNOWN["m10"][2])
    assert bt._t1_r58(a) is None and bt._t1_r67(a) is None and bt._t1_r51(a)["file"] == "T1_s2rcMSBm10"
    a = parse(KNOWN["m20"][2])
    assert bt._t1_r58(a)["file"] == "T1_s2rcMSBm20" and bt._t1_r67(a) is None
    for name, opt in _table(_sh()).items():
        a = parse(f"{BASE} {opt} --captain-ruling R67")
        got = bt._t1_r67(a)
        assert got["file"] == name and got["option"] == f"`{opt}`" and got["captain_ruling"] == "R67"
        assert bt._t1_r58(a) is None and bt._t1_pooled(a) is None and bt._r67_ordered(a)


# ------------------------------------------------------------------ the launcher, run for real


def _repo(tmp_path: Path, tries: str = ROW, launcher: str | None = None, **kw):
    """tests/test_t1_r51_launcher.py's throwaway repository (its launcher path holds this launcher, or the given copy
    of it) with a stand-in `uv` that fails."""
    return r51l._repo(tmp_path, tries=tries, launcher=_sh() if launcher is None else launcher, **kw)


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return r51l._run(repo, env, *argv)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1r67" / run


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
    o = repo / "runs" / f"t1r67_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", ONE).returncode == 10
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "T1_s2rcMSBm075")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", ONE)
    assert other.returncode == 10 and (repo / "runs" / f"t1r67_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_and_a_prefix_are_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1", ONE)
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()
    short = _run(repo, env, head[:12], "r2", ONE)
    assert short.returncode == 96
    _stopped_at_the_gate(tmp_path, "r2", 96)


@pytest.mark.parametrize("argv, code", [
    ([], None), (["r1"], None), (["r1", ""], None),  # missing arguments (an empty name is none)
    (["a/b", ONE], 90), (["r 1", ONE], 90), (["$(id)", ONE], 90),
    (["r1", ONE, "extra"], 89),
    *[(["r1", bad], 89) for bad in (
        "T1_s2rcMSB", "T1_s2rcMSBm10", "T1_s2rcMSBm20", "T1_s2rcMSBm12", "T1_s2rcMSBm1.25", "T1_s2rcMSBm075 ",
        "T1_s2rcMS20m125", "T1_s2rcMSBm175", "T1_s2rcMS175m10", "t1_s2rcMSBm125", "T1_s2rcMSBd125", "$(id)", "*")],
])
def test_odd_arguments_and_unknown_names_are_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, *argv)
    assert out.returncode != 0 if code is None else out.returncode == code, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW}more\n")
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ""}, 95, "exactly one row R67"),
    ({"tries": ROW + ROW}, 95, "exactly one row R67"),  # two rows with one ID
    ({"tries": ROW.replace("| R67 (", "| R67b (")}, 95, "exactly one row R67"),  # another ID
    ({"tries": ROW.replace("R67 (the captain's order", "R67 (a proposal")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("R67 (the captain's", "R67 awaits (the captain's")}, 95, "not the captain's ruling"),
    ({"tries": ROW.replace("ruling, no run", "a try")}, 95, "does not say ruling, no run"),
    # the file's string inside a longer backticked string, or with no backticks, is not named
    ({"tries": ROW.replace(f"`{ONE_OPT}`", f"`{ONE_OPT} --herd`")}, 95, f"does not name `{ONE_OPT}`"),
    ({"tries": ROW.replace(f"`{ONE_OPT}`", ONE_OPT)}, 95, f"does not name `{ONE_OPT}`"),
    ({"tries": ROW.replace("`--captain-ruling R67`", "`--captain-ruling R58`")}, 95,
     "does not name `--captain-ruling R67`"),
    ({"tries": ROW.replace("`--captain-ruling R67`", "`--captain-ruling R67 now`")}, 95,
     "does not name `--captain-ruling R67`"),
    # the string named, but paired with another file's name
    ({"tries": ROW.replace(f"`{ONE_OPT}` ({ONE})", f"`{ONE_OPT}` (T1_s2rcMSBm15)")}, 95,
     f"does not pair `{ONE_OPT}` with ({ONE})"),
    ({"tries": (f"| R67 (the captain's order of X; ruling, no run) | T1 | `{ONE_OPT}` ({ONE}) "
                "`--captain-ruling R67` |\n")}, 95, f"does not name `{ONE_OPT}`"),  # named outside the first cell only
    ({"relabel": False}, 92, "feng_relabel_c is not in this checkout"),
    ({"venv": None}, 88, "no .venv of its own"),
    ({"venv": "other"}, 88, "no .venv of its own"),  # a venv whose vec.pth names another checkout's src
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"t1r67_{head[:7]}_r1").exists()


@pytest.mark.parametrize("name", NAMES)
def test_the_gate_takes_a_short_row_in_its_own_words(tmp_path, name):
    tries = (f"| R67 (the captain's order of X, ruling, no run: `{OPTS[name]}` ({name}), `--captain-ruling R67`) "
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
    o = repo / "runs" / f"t1r67_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_past_the_tests_step_stops_at_the_first_unknown_bytes(tmp_path):
    """tests/test_t1_l67pb_launcher.py's stand-in `uv`: the tests pass, a build writes a stand-in file whose stamp
    names HEAD and clean, `vec check` passes; no stand-in file has a known file's bytes, so the run stops at MSB."""
    repo, env, head = _repo(tmp_path)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(PASSING_UV)
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == 23, out.stderr
    assert "stopped: step MSB does not give the known bytes of T1_s2rcMSB" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert {p.name: p.read_text().strip() for p in s.glob("*.exit")} == {"gate.exit": "0", "tests.exit": "0",
                                                                          "MSB.exit": "23"}
    assert (s / "stamps").read_text() == f"MSB/T1_s2rcMSB {head} True\n"
    o = repo / "runs" / f"t1r67_{head[:7]}_r1"
    assert sorted(p.name for p in o.iterdir()) == ["MSB"]


# A stand-in `uv` for whole runs: the tests pass; a build writes "file <its output's name>" with a stamp naming HEAD
# and clean (the checker's base writes T1_s2rcMSB's stand-in, as each file's path with MS's damping at B0's must, or
# BASE_GIVES's); with B_DIFFERS set the second build of the file writes other bytes; `vec check` passes; every builder
# and checker call is listed in calls.txt; the stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/build_test.py|experiments/builds/t1_r67_check.py)
        prog=$1; shift
        echo "$prog $*" >> calls.txt
        [ "$prog" = experiments/builds/t1_r67_check.py ] && shift
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        [ "$prog" = experiments/builds/t1_r67_check.py ] && name=${BASE_GIVES:-T1_s2rcMSB}
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
    r = f"runs/t1r67_{head[:7]}_r1"
    o = repo / r
    folders = steps[2:-1]
    assert sorted(p.name for p in o.iterdir()) == sorted([*folders, "SHA256SUMS"])
    stamps = (s / "stamps").read_text().splitlines()
    assert [line.split()[0].split("/")[0] for line in stamps] == folders
    assert all(line.endswith(f" {head} True") for line in stamps)
    listed = (o / "SHA256SUMS").read_text().splitlines()
    assert [line.split()[1].split("/")[0] for line in listed] == folders and (s / "SHA256SUMS").exists()
    stem, msb = _stem(name), hashlib.sha256(b"file T1_s2rcMSB\n").hexdigest()
    assert dict(line.split()[::-1] for line in listed)[f"{stem}base/{name}_base.h5ad"] == msb
    args = f"{BASE} {OPTS[name]} --captain-ruling R67"
    want = [*(f"experiments/builds/build_test.py t1 {v[2]} --out {r}/{step}/{v[0]}.h5ad" for step, v in KNOWN.items()),
            f"experiments/builds/t1_r67_check.py base t1 {args} --out {r}/{stem}base/{name}_base.h5ad",
            f"experiments/builds/build_test.py t1 {args} --out {r}/{stem}a/{name}.h5ad",
            f"experiments/builds/build_test.py t1 {args} --out {r}/{stem}b/{name}.h5ad"]
    assert (repo / "calls.txt").read_text().splitlines() == want
    assert (o / f"{stem}a" / f"{name}.h5ad").read_bytes() == (o / f"{stem}b" / f"{name}.h5ad").read_bytes()


@pytest.mark.parametrize("keep, step, code", [(("MSB_SHA",), "MSB", 23), (("M10_SHA",), "m10", 24),
                                              (("M20_SHA",), "m20", 25)])
def test_an_identity_file_of_other_bytes_stops_the_run_at_its_step(tmp_path, keep, step, code):
    repo, head, out = _stand_in_run(tmp_path, ONE, keep)
    assert out.returncode == code, out.stderr
    assert f"stopped: step {step} does not give the known bytes of {KNOWN[step][0]}" in out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(ONE)
    ran = steps[:steps.index(step) + 1]
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {
        **{x: "0" for x in ran[:-1]}, step: str(code)}
    assert not (repo / "runs" / f"t1r67_{head[:7]}_r1" / f"{_stem(ONE)}a").exists()


def test_a_base_check_of_other_bytes_stops_the_run_before_the_file_is_built(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, ONE, BASE_GIVES=ONE)
    assert out.returncode == 51, out.stderr
    stem = _stem(ONE)
    assert f"stopped: step {stem}base does not give the known bytes of {ONE}_base" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / f"{stem}base.exit").read_text().strip() == "51" and not (s / f"{stem}a.exit").exists()
    assert not (repo / "runs" / f"t1r67_{head[:7]}_r1" / f"{stem}a").exists()


def test_two_builds_of_other_bytes_stop_the_run_at_done(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, ONE, B_DIFFERS="1")
    assert out.returncode == 41, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "done.exit").read_text().strip() == "41" and (s / "SHA256SUMS").exists()
