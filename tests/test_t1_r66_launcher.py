"""The test-board build launcher of row R66's thirty files (T1's third-round step-on slots on T1_s2rcMS20's command),
one file per call.

experiments/builds/t1_r66_line.sh refuses any name but the thirty before it makes a folder; its gate checks the commit,
the clean tree, the rules file's hash, row R66's first cell (the file's option string and `--captain-ruling R66` as
whole backticked strings, the string paired there with the file's name), the relabelling's output and the venv; then it
runs the targeted tests, rebuilds T1_s2rcMS20 and the line-12 files whose step the file uses and stops unless each
gives its known bytes; then it builds the file once on its own path with every step at T1_s2rcMS20's value
(experiments/builds/t1_r66_check.py base) and stops unless that is T1_s2rcMS20's bytes; then it builds the file twice
from one clean commit and compares the two. These tests hold its table to the builder's, its hashes to the committed
lists, its commands to the ones the builder takes, check that it deletes nothing and reads no output of a test build
but exit codes, sha256 and the stamps' commit and clean fields, and run it for real in throwaway repositories under a
temporary home: with a stand-in `uv` that fails (every refusal path, and a call passing the gate stops at the tests
step and builds nothing) and with stand-ins that pass (a whole run of each of the thirty names goes through every step
in order; a stand-in file of other bytes stops the run where its step says). Nothing here reads the repository's own
ledger or rules file for its state of today; build_test is imported only to check its table and that it takes the
commands the launcher runs, against a temporary ledger."""

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
import test_t1_l12_launcher as l12l
import test_t1_r51_launcher as r51l
from test_t1_l12_build import REAL_R46, REAL_R60
from test_t1_l67pb_launcher import PASSING_UV, _sums
from test_t1_r66_build import OPTS, REAL_R66

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_r66_line.sh"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
MS20_OPT = "--feng-fraction-mult 2.0 --direction-momentum all"
L12 = f"{BASE} {MS20_OPT}"
# line-12 step -> (file name, args variable, the args, hash variable, sha256, exit code), as
# experiments/builds/t1_l12_line.sh builds them and t1l12_128dfa1_b1.SHA256SUMS holds them
LINE12 = {
    "d10": ("T1_s2rcMS20d10", "D10", f"{L12} --g-mom-damping 1.0 --captain-ruling R60", "D10_SHA",
            "16f35132bb96919f1dcc89ace1dd12f5a3ced88749e889086904e189ab9fd2a8", 71),
    "d15": ("T1_s2rcMS20d15", "D15", f"{L12} --g-mom-damping 1.5 --captain-ruling R60", "D15_SHA",
            "cfea3b689696bd2e08abe04ad607ab5a826b300ad6d841c3f2d8fc475f536f63", 72),
    "d20": ("T1_s2rcMS20d20", "D20", f"{L12} --g-mom-damping 2.0 --captain-ruling R60", "D20_SHA",
            "c5340b1505e4e032e764465d69973910b7071cedbd266487866606c0dafb8cd9", 73),
    "H05": ("T1_s2rcMS20H05", "H05", f"{L12} --held-momentum 0.5 --captain-ruling R60", "H05_SHA",
            "97589938a95ef957d2aed1117f84147cd027f2a0c3d1c48fbcaf3c67570a87cc", 74),
    "H10": ("T1_s2rcMS20H10", "H10", f"{L12} --held-momentum 1.0 --captain-ruling R60", "H10_SHA",
            "5888b0f4da0124e7561d17a002006d79b26efe48c14e346cf3e45faf1222be3e", 75),
    "m10": ("T1_s2rcMS20m10", "M10", f"{L12} --direction-momentum-damping 1.0 --captain-ruling R60", "M10_SHA",
            "316b67b44df2501db6f5d43ea099bec107ba9422498ed90f1861fc32d70cc1c5", 76),
    "m20": ("T1_s2rcMS20m20", "M20", f"{L12} --direction-momentum-damping 2.0 --captain-ruling R60", "M20_SHA",
            "961df1c7b4b8cce40006729f8439701313cb5d2b41b159763e20d1dbd5097000", 77),
}
MS20 = ("T1_s2rcMS20", "MS20", f"{BASE} {MS20_OPT} --captain-ruling R51", "MS20_SHA",
        "629d3a1cffa8151f0345a169ba97420e4ee3b43e379bc3bf7d33ab947dcb62ef", 39)
NAMES = list(OPTS)
ROW_R66 = f"| {REAL_R66} | T1 | step-on |\n"
ONE = "T1_s2rcMS20d15H05m20"  # the file the gate's refusals are run on


def _sh() -> str:
    return LAUNCHER.read_text()


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _table(sh: str) -> dict[str, str]:
    """The launcher's FILE_NAME -> option string table (its case lines), $MS20_OPT written out."""
    rows = re.findall(r'^  (T1_s2rcMS20[A-Za-z0-9]+)\) OPT="\$MS20_OPT ([^"]+)" ;;$', sh, flags=re.M)
    return {name: f"{MS20_OPT} {rest}" for name, rest in rows}


def _stem(name: str) -> str:
    return name.removeprefix("T1_s2rcMS20")


def _ids(name: str) -> list[str]:
    """The line-12 files whose step the file uses, as the header fixes them: on each axis, the line-12 file at the
    file's value, or every line-12 file on that axis where line 12 has none at that value (d25, H20)."""
    out = []
    for step in re.findall(r"[dHm]\d\d", _stem(name)):
        out += [step] if step in LINE12 else [s for s in LINE12 if s[0] == step[0]]
    return out


def _steps(name: str) -> list[str]:
    stem = _stem(name)
    return ["gate", "tests", "MS20", *_ids(name), f"{stem}base", f"{stem}a", f"{stem}b", "done"]


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launchers_table_is_the_builders_one_table():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    sh = _sh()
    table = _table(sh)
    assert list(table.items()) == [(f, v[0]) for f, v in bt.T1_R66_FILES.items()] == list(OPTS.items())
    assert f'MS20_OPT="{MS20_OPT}"' in sh and f'L12="{L12}"' in sh
    assert 'ARGS="--g-mom --base lib3 --feng-groups recomputed $OPT --captain-ruling R66"' in sh
    assert '*) echo "refused: FILE_NAME is not one of row R66\'s thirty files" >&2; exit 89 ;;' in sh
    assert sh.count('OPT="$MS20_OPT ') == 30


def test_the_known_hashes_are_the_committed_lists_and_the_args_line_12s():
    sh = _sh()
    l12 = _sums("t1l12_128dfa1_b1.SHA256SUMS")
    r51 = _sums("t1r51_0cdcc6f_r1.SHA256SUMS")
    l12_sh = (REPO / "experiments" / "builds" / "t1_l12_line.sh").read_text()
    args = dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    for step, (name, key, full, var, sha, _) in LINE12.items():
        assert l12[f"{step}a/{name}.h5ad"] == l12[f"{step}b/{name}.h5ad"] == sha, step
        assert f"{var}={sha}" in sh and sh.count(sha) == 1
        assert args[key].replace("$L12", L12) == full
        assert f'{key}_ARGS="{args[key]}"' in l12_sh  # the line-12 launcher's own command
    name, key, full, var, sha, _ = MS20
    assert r51["MS20a/T1_s2rcMS20.h5ad"] == r51["MS20b/T1_s2rcMS20.h5ad"] == l12["MS20/T1_s2rcMS20.h5ad"] == sha
    assert f"{var}={sha}" in sh and sh.count(sha) == 1 and args[key] == full and f'MS20_ARGS="{full}"' in l12_sh
    assert set(args) == {"MS20", *(v[1] for v in LINE12.values())}


def test_the_steps_run_in_order_and_each_build_is_checked_and_stamped():
    sh = _sh()
    assert "identity MS20 T1_s2rcMS20 \"$MS20_SHA\" 39 $MS20_ARGS" in sh
    for step, (name, key, _, var, _, code) in LINE12.items():
        assert f'    {step}) identity {step} {name} "${var}" {code} ${key}_ARGS ;;' in sh
    marks = ("STEP=gate", "STEP=tests", "identity MS20 ", "for id in $IDS; do", "STEP=${STEM}base",
             'for STEP in "${STEM}a" "${STEM}b"; do begin; build "$STEP" "$FILE" $ARGS || finish $?; step_ok; done',
             "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    base = sh[sh.index("STEP=${STEM}base"):sh.index('for STEP in "${STEM}a"')]
    assert 'build_base "$STEP" "${FILE}_base" $ARGS || finish $?' in base
    assert 'same "$STEP" "${FILE}_base" "$MS20_SHA" 51' in base
    build_base = sh[sh.index("build_base() {"):sh.index("same() {")]
    assert ('$PY experiments/builds/t1_r66_check.py base t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12'
            in build_base) and 'checked "$d" "$name"' in build_base
    build = sh[sh.index("build() {"):sh.index("build_base() {")]
    assert ('$PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11'
            in build) and 'checked "$d" "$name"' in build
    done = sh[sh.index("STEP=done"):]
    assert '(cd "$O" && sha256sum "${BUILT[@]}") > "$S/SHA256SUMS" || finish 40' in done
    assert 'cmp -s "$O/${STEM}a/$FILE.h5ad" "$O/${STEM}b/$FILE.h5ad" || finish 41' in done
    assert 'got=$(sha256sum "$O/$1/$2.h5ad" | cut -d\' \' -f1)' in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("identity MS20 ")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    first_uv = sh.index("uv run")
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < first_uv and ".light" not in sh
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert "--seed" not in sh and "--board T1" not in sh and "vec check" in sh and "BOARD=T1:test" in sh
    assert "upload" not in sh.replace("Nothing here uploads", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")
    assert "O=runs/t1r66_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1r66/$RUN" in sh


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for line in lines:  # every builder, checker build and check call: output to /dev/null
        if re.search(r'build_test\.py t1|t1_r66_check\.py base|vec check "', line):
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
    ident = sh[sh.index("identity() {"):sh.index("identity MS20 ")]
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
    assert 'r66=$(first_cell R66 95) || exit $?' in gate
    assert 'for opt in "$OPT" "--captain-ruling R66"; do' in gate
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in gate  # whole backticked strings, line by line
    assert "grep -qF -- \"\\`$OPT\\` ($FILE)\"" in gate  # the string paired with the file's name
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 92, 88, 93]


def test_the_tests_step_runs_the_new_tests_and_t1_l12_line_shs():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    l12 = l12l.LAUNCHER.read_text()
    l12_step = l12[l12.index("STEP=tests"):l12.index("step_ok", l12.index("STEP=tests"))]
    assert {"tests/test_t1_r66_launcher.py", "tests/test_t1_r66_build.py",
            *re.findall(r"tests/test_[a-z0-9_]+\.py", l12_step)} == set(listed)
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
    used |= {int(c) for c in re.findall(r'same "\$STEP" \S+ "\$MS20_SHA" (\d+)', sh)}
    documented = {int(c) for c in re.findall(r"(?:^|[\s,#])(\d\d) ", doc, flags=re.M)}
    assert used == documented, (used, documented)


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    sys.path.insert(0, str(REPO / "experiments" / "external"))
    import build_test as bt
    import feng_relabel as fr

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with the ruling rows' orders, not the repository's
    ledger.write_text(r51l.ROWS_ALL + r51l.ROW_R51 + f"| {REAL_R60} | T1 |\n| {REAL_R46} | T1 |\n{ROW_R66}")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    known = [MS20, *LINE12.values()]
    for name, _, full, *_ in known:
        a = bt.build_parser().parse_args(["t1", *full.split(), "--out", "x.h5ad"])
        got = bt._t1_r51(a)  # nothing read but the ledger
        assert got["file"] == name and got["captain_ruling"] in ("R51", "R60")
    for name, opt in _table(_sh()).items():
        a = bt.build_parser().parse_args(["t1", *f"{BASE} {opt} --captain-ruling R66".split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        got = bt._t1_r51(a)
        assert got["file"] == name and got["option"] == f"`{opt}`" and got["captain_ruling"] == "R66"
        assert bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None and bt._t1_pooled(a) is None


def test_the_line_12_files_a_name_rebuilds():
    assert _ids("T1_s2rcMS20d15H05m20") == ["d15", "H05", "m20"]
    assert _ids("T1_s2rcMS20d10H10") == ["d10", "H10"] and _ids("T1_s2rcMS20H05m10") == ["H05", "m10"]
    assert _ids("T1_s2rcMS20d25") == ["d10", "d15", "d20"] and _ids("T1_s2rcMS20H20") == ["H05", "H10"]
    assert all(1 <= len(_ids(n)) <= 3 for n in NAMES)


# ------------------------------------------------------------------ the launcher, run for real


def _repo(tmp_path: Path, tries: str = ROW_R66, launcher: str | None = None, **kw):
    """tests/test_t1_r51_launcher.py's throwaway repository (its launcher path holds this launcher, or the given copy
    of it) with a stand-in `uv` that fails."""
    return r51l._repo(tmp_path, tries=tries, launcher=_sh() if launcher is None else launcher, **kw)


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return r51l._run(repo, env, *argv)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1r66" / run


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no head file, no later step's mark
    assert not list((s / "tmp").iterdir())


@pytest.mark.parametrize("name", NAMES)
def test_each_name_passes_the_gate_on_the_real_row_and_stops_at_the_failing_tests_step(tmp_path, name):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1", name)
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t1r66_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", ONE).returncode == 10
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "T1_s2rcMS20d25")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", ONE)
    assert other.returncode == 10 and (repo / "runs" / f"t1r66_{head[:7]}_r2").is_dir()


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
        "T1_s2rcMS20", "T1_s2rcMS20d10", "T1_s2rcMS20H05", "T1_s2rcMS20m20",  # MS20 and line 12's files
        "T1_s2rcMS20H05d10", "T1_s2rcMS20m10H05",  # suffixes out of d, H, m order
        "T1_s2rcMS20d25H05", "T1_s2rcMS20d10H20", "T1_s2rcMS20m30", "T1_s2rcMS20d10H05m10x", "t1_s2rcMS20d10H05",
        "T1_s2rcMSBd10H05", "T1_s2rcB25MSH05", "T1_s2rcMS20d10H05 ", "$(id)", "*")],
])
def test_odd_arguments_and_unknown_names_are_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, *argv)
    assert out.returncode != 0 if code is None else out.returncode == code, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R66}more\n")
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


ONE_OPT = OPTS[ONE]


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ""}, 95, "exactly one row R66"),
    ({"tries": ROW_R66 + ROW_R66}, 95, "exactly one row R66"),  # two rows with one ID
    ({"tries": ROW_R66.replace("| R66 (", "| R66b (")}, 95, "exactly one row R66"),  # another ID
    ({"tries": ROW_R66.replace("R66 (the captain's order", "R66 (a proposal")}, 95, "not the captain's ruling"),
    ({"tries": ROW_R66.replace("R66 (the captain's", "R66 awaits (the captain's")}, 95, "not the captain's ruling"),
    ({"tries": ROW_R66.replace("ruling, no run", "a try")}, 95, "does not say ruling, no run"),
    # the file's string inside a longer backticked string, or with no backticks, is not named
    ({"tries": ROW_R66.replace(f"`{ONE_OPT}`", f"`{ONE_OPT} --herd`")}, 95, f"does not name `{ONE_OPT}`"),
    ({"tries": ROW_R66.replace(f"`{ONE_OPT}`", ONE_OPT)}, 95, f"does not name `{ONE_OPT}`"),
    ({"tries": ROW_R66.replace("`--captain-ruling R66`", "`--captain-ruling R60`")}, 95,
     "does not name `--captain-ruling R66`"),
    ({"tries": ROW_R66.replace("`--captain-ruling R66`", "`--captain-ruling R66 now`")}, 95,
     "does not name `--captain-ruling R66`"),
    # the string named, but paired with another file's name
    ({"tries": ROW_R66.replace(f"`{ONE_OPT}` ({ONE})", f"`{ONE_OPT}` (T1_s2rcMS20d15H05m10)")}, 95,
     f"does not pair `{ONE_OPT}` with ({ONE})"),
    ({"tries": (f"| R66 (the captain's order of X; ruling, no run) | T1 | `{ONE_OPT}` ({ONE}) "
                "`--captain-ruling R66` |\n")}, 95, f"does not name `{ONE_OPT}`"),  # named outside the first cell only
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
    assert not (repo / "runs" / f"t1r66_{head[:7]}_r1").exists()


def test_the_gate_takes_a_short_row_in_its_own_words(tmp_path):
    tries = f"| R66 (the captain's order of X, ruling, no run: `{ONE_OPT}` ({ONE}), `--captain-ruling R66`) | T1 |\n"
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == 10, out.stderr  # past the gate: the stand-in uv fails the tests step
    assert (_scratch(tmp_path, "r1") / "gate.exit").read_text().strip() == "0"


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t1r66_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_past_the_tests_step_stops_at_the_first_unknown_bytes(tmp_path):
    """tests/test_t1_l67pb_launcher.py's stand-in `uv`: the tests pass, a build writes a stand-in file whose stamp
    names HEAD and clean, `vec check` passes; no stand-in file has a known file's bytes, so the run stops at MS20."""
    repo, env, head = _repo(tmp_path)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(PASSING_UV)
    out = _run(repo, env, head, "r1", ONE)
    assert out.returncode == 39, out.stderr
    assert "stopped: step MS20 does not give the known bytes of T1_s2rcMS20" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert {p.name: p.read_text().strip() for p in s.glob("*.exit")} == {"gate.exit": "0", "tests.exit": "0",
                                                                          "MS20.exit": "39"}
    assert (s / "stamps").read_text() == f"MS20/T1_s2rcMS20 {head} True\n"
    o = repo / "runs" / f"t1r66_{head[:7]}_r1"
    assert sorted(p.name for p in o.iterdir()) == ["MS20"]


# A stand-in `uv` for whole runs: the tests pass; a build writes "file <its output's name>" with a stamp naming HEAD
# and clean (the checker's base writes T1_s2rcMS20's stand-in, as each new path with its steps at MS20's values must,
# or BASE_GIVES's); with B_DIFFERS set the second build of the file writes other bytes; `vec check` passes; every
# builder and checker call is listed in calls.txt; the stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/build_test.py|experiments/builds/t1_r66_check.py)
        prog=$1; shift
        echo "$prog $*" >> calls.txt
        [ "$prog" = experiments/builds/t1_r66_check.py ] && shift
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        [ "$prog" = experiments/builds/t1_r66_check.py ] && name=${BASE_GIVES:-T1_s2rcMS20}
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
    for name, _, _, var, sha, _ in (MS20, *LINE12.values()):
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
    r = f"runs/t1r66_{head[:7]}_r1"
    o = repo / r
    folders = steps[2:-1]
    assert sorted(p.name for p in o.iterdir()) == sorted([*folders, "SHA256SUMS"])
    stamps = (s / "stamps").read_text().splitlines()
    assert [line.split()[0].split("/")[0] for line in stamps] == folders
    assert all(line.endswith(f" {head} True") for line in stamps)
    listed = (o / "SHA256SUMS").read_text().splitlines()
    assert [line.split()[1].split("/")[0] for line in listed] == folders and (s / "SHA256SUMS").exists()
    stem, ms20 = _stem(name), hashlib.sha256(b"file T1_s2rcMS20\n").hexdigest()
    assert dict(line.split()[::-1] for line in listed)[f"{stem}base/{name}_base.h5ad"] == ms20
    args = f"{BASE} {OPTS[name]} --captain-ruling R66"
    want = [f"experiments/builds/build_test.py t1 {MS20[2]} --out {r}/MS20/T1_s2rcMS20.h5ad",
            *(f"experiments/builds/build_test.py t1 {LINE12[i][2]} --out {r}/{i}/{LINE12[i][0]}.h5ad"
              for i in _ids(name)),
            f"experiments/builds/t1_r66_check.py base t1 {args} --out {r}/{stem}base/{name}_base.h5ad",
            f"experiments/builds/build_test.py t1 {args} --out {r}/{stem}a/{name}.h5ad",
            f"experiments/builds/build_test.py t1 {args} --out {r}/{stem}b/{name}.h5ad"]
    assert (repo / "calls.txt").read_text().splitlines() == want
    assert (o / f"{stem}a" / f"{name}.h5ad").read_bytes() == (o / f"{stem}b" / f"{name}.h5ad").read_bytes()


@pytest.mark.parametrize("name, keep, step, code", [
    ("T1_s2rcMS20d15H05m20", ("D15_SHA",), "d15", 72),
    ("T1_s2rcMS20d15H05m20", ("M20_SHA",), "m20", 77),
    ("T1_s2rcMS20d25", ("D20_SHA",), "d20", 73),
    ("T1_s2rcMS20H20", ("H05_SHA",), "H05", 74),
    ("T1_s2rcMS20H10m10", ("M10_SHA",), "m10", 76),
])
def test_a_line_12_file_of_other_bytes_stops_the_run_at_its_step(tmp_path, name, keep, step, code):
    repo, head, out = _stand_in_run(tmp_path, name, keep)
    assert out.returncode == code, out.stderr
    assert f"stopped: step {step} does not give the known bytes of T1_s2rcMS20{step}" in out.stderr
    s = _scratch(tmp_path, "r1")
    steps = _steps(name)
    ran = steps[:steps.index(step) + 1]
    assert {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")} == {
        **{x: "0" for x in ran[:-1]}, step: str(code)}


def test_a_base_check_of_other_bytes_stops_the_run_before_the_file_is_built(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, ONE, BASE_GIVES=ONE)
    assert out.returncode == 51, out.stderr
    stem = _stem(ONE)
    assert f"stopped: step {stem}base does not give the known bytes of {ONE}_base" in out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / f"{stem}base.exit").read_text().strip() == "51" and not (s / f"{stem}a.exit").exists()
    assert not (repo / "runs" / f"t1r66_{head[:7]}_r1" / f"{stem}a").exists()


def test_two_builds_of_other_bytes_stop_the_run_at_done(tmp_path):
    repo, head, out = _stand_in_run(tmp_path, ONE, B_DIFFERS="1")
    assert out.returncode == 41, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "done.exit").read_text().strip() == "41" and (s / "SHA256SUMS").exists()
