"""The test-board build launcher of row R51's seven files (T1's line 10, each one change to T1_s2rcMSB's command:
T1_s2rcMS20, T1_s2rcMS30, T1_s2rcMSBd10, T1_s2rcMSBd15, T1_s2rcMSBH05, T1_s2rcMSBH10, T1_s2rcMSBm10), each only with
--captain-ruling R51.

experiments/builds/t1_r51_line.sh first rebuilds nineteen known files (t1_hm_line.sh's fifteen and the four
t1_l67pb_line.sh built twice with one sha256) and stops unless each gives its known bytes; then it builds each new
command once on its own path with its one setting at T1_s2rcMSB's value (experiments/builds/t1_r51_check.py msb) and
stops unless each gives MSB's bytes; then it builds the seven files twice each from one clean commit, runs the rows
check on MSB's file and the two H files and compares each pair. These tests hold its command strings, hashes and step
order to that, check that it deletes nothing and reads no output of a test build but exit codes, sha256, the stamps'
commit and clean fields and the rows check's booleans, and run its gate's paths for real in throwaway repositories under
a temporary home, with a stand-in `uv` that fails, so that a call passing the gate stops at the tests step and builds
nothing (and, with stand-ins that pass, that a whole run goes through every step in order). Nothing here reads the
repository's own ledger or rules file for its state of today; build_test is imported only to check that it takes the
commands the launcher runs, against a temporary ledger."""

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
from test_t1_l67pb_launcher import KNOWN as KNOWN15
from test_t1_l67pb_launcher import PASSING_UV, _sums

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_r51_line.sh"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
# step -> (file name, args variable, the args, hash variable, sha256, exit code): the fifteen of t1_hm_line.sh and the
# four t1_l67pb_line.sh built twice with one sha256
KNOWN = {step: (name, key, f"{BASE}{tail}", var, sha, code) for step, (name, key, tail, var, sha, code) in
         KNOWN15.items()}
KNOWN |= {
    "Cm30": ("T1_s2rcCm30", "CM30", f"{BASE} --g-mom-damping 1.0 --feng-fraction-mult 3.0 --captain-ruling R22",
             "CM30_SHA", "807fe097c88ee779430896002a34635aeffbc99b433c68ebdb7bd41d8af2cbbe", 35),
    "B25d15": ("T1_s2rcB25d15", "B25D15", f"{BASE} --g-mom-damping 1.5 --feng-fraction-mult 2.5 --captain-ruling R22",
               "B25D15_SHA", "5bc2b32e4688439f29c54c1b1203325595d996c37df30ba2b9ff26f27a60f4c6", 36),
    "B25d20": ("T1_s2rcB25d20", "B25D20", f"{BASE} --g-mom-damping 2.0 --feng-fraction-mult 2.5 --captain-ruling R22",
               "B25D20_SHA", "e6520dbee0fdc3f1336bc3f9de1e99e8a501de2d440b3f9085c88b5c7ecfb3f6", 37),
    "pooled": ("T1_s2rc_pooled", "POOLED", "--g-mom --base pooled --feng-groups recomputed --captain-ruling R37",
               "POOLED_SHA", "91341a5f8c71dcc103a47f9c19a5d7794b1358c1c642c5c9251624eb6feb8bdc", 38),
}
# file -> its option string, as row R51 names it
OPTS = {"T1_s2rcMS20": "--feng-fraction-mult 2.0 --direction-momentum all",
        "T1_s2rcMS30": "--feng-fraction-mult 3.0 --direction-momentum all",
        "T1_s2rcMSBd10": "--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.0",
        "T1_s2rcMSBd15": "--feng-fraction-mult 1.5 --direction-momentum all --g-mom-damping 1.5",
        "T1_s2rcMSBH05": "--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 0.5",
        "T1_s2rcMSBH10": "--feng-fraction-mult 1.5 --direction-momentum all --held-momentum 1.0",
        "T1_s2rcMSBm10": "--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 1.0"}
KEYS = {"T1_s2rcMS20": ("MS20", "MS20"), "T1_s2rcMS30": ("MS30", "MS30"), "T1_s2rcMSBd10": ("D10", "d10"),
        "T1_s2rcMSBd15": ("D15", "d15"), "T1_s2rcMSBH05": ("H05", "H05"), "T1_s2rcMSBH10": ("H10", "H10"),
        "T1_s2rcMSBm10": ("M10", "m10")}  # file -> (args variable, step stem)
NEW = {KEYS[f][0]: f"{BASE} {opt} --captain-ruling R51" for f, opt in OPTS.items()}
# new step -> (file name, args variable, exit code when it is not MSB's bytes)
MSBS = {f"{KEYS[f][1]}msb": (f"{f}_msb", KEYS[f][0], 51 + i) for i, f in enumerate(OPTS)}
# new step pair -> (file name, args variable, exit code of the pair's comparison)
PAIRS = {(f"{KEYS[f][1]}a", f"{KEYS[f][1]}b"): (f, KEYS[f][0], 41 + i) for i, f in enumerate(OPTS)}
STEPS = ["gate", "tests", *KNOWN, *MSBS, *(s for pair in PAIRS for s in pair), "rows", "done"]
MSB_SHA = KNOWN["MSB"][4]


def _sh() -> str:
    return LAUNCHER.read_text()


def _code(sh: str) -> str:
    return "\n".join(line.split("#")[0] for line in sh.splitlines())


def _args(sh: str) -> dict[str, str]:
    return dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launcher_reproduces_nineteen_files_checks_msbs_bytes_then_builds_each_new_file_twice():
    sh = _sh()
    want = {key: args for _, key, args, _, _, _ in KNOWN.values()} | NEW
    assert _args(sh) == want
    for step, (name, key, _, var, sha, code) in KNOWN.items():
        assert f"build {step} {name} ${key}_ARGS || finish $?" in sh
        assert f'same {step} {name} "${var}" {code}' in sh
        assert f"{var}={sha}" in sh and sh.count(sha) == 1
    for step, (name, key, code) in MSBS.items():
        assert f"build_msb {step} {name} ${key}_ARGS || finish $?" in sh
        assert f'same {step} {name} "$MSB_SHA" {code}' in sh  # T1_s2rcMSB's bytes, rebuilt in step MSB
    assert 'got=$(sha256sum "$O/$1/$2.h5ad" | cut -d\' \' -f1)' in sh
    assert '[ "$got" = "$3" ] || {' in sh and 'finish "$4"; }' in sh
    for (a, b), (name, key, code) in PAIRS.items():
        assert f'for STEP in {a} {b}; do begin; build "$STEP" {name} ${key}_ARGS || finish $?; step_ok; done' in sh
        assert f'cmp -s "$O/{a}/{name}.h5ad" "$O/{b}/{name}.h5ad" || finish {code}' in sh
    order = [sh.index(s) for s in ("STEP=gate", "STEP=tests", *(f"STEP={k}\n" for k in KNOWN),
                                   *(f"STEP={k}\n" for k in MSBS), *(f"for STEP in {a} {b}" for a, b in PAIRS),
                                   "STEP=rows", "STEP=done")]
    assert order == sorted(order)
    listing = sh[sh.index("(cd \"$O\" && sha256sum"):sh.index('> "$S/SHA256SUMS"')]
    assert re.findall(r"([A-Za-z0-9]+)/\*\.h5ad", listing) == STEPS[2:-2]  # every built folder, once, in order
    assert 'cp "$S/SHA256SUMS" "$O/SHA256SUMS"' in sh
    assert "vec check" in sh and "BOARD=T1:test" in sh
    assert "--seed" not in sh and "--board T1" not in sh  # seed and board: the builder's defaults, as the references
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh  # the builds, as the references
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh  # the tests
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=D\n")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    first_uv = sh.index("uv run")
    assert sh.index("PATH=$HOME/.local/bin:$PATH") < first_uv  # PATH set before the first uv call
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv  # uv runs this checkout's own .venv
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh  # the full hash, never a prefix
    assert '[ "$st" = "$head True" ] || return 15' in sh  # a build counts only with a stamp of this commit, clean
    assert 'echo "$(basename "$d")/$name $st" >> "$S/stamps"' in sh
    assert "upload" not in sh.replace("Nothing here uploads", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")  # runs from the folder it is called in
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")  # the rules file is checked by its hash only


def test_the_msb_builds_go_through_the_checker_and_every_build_is_checked_and_stamped():
    sh = _sh()
    msb = sh[sh.index("build_msb() {"):sh.index("same() {")]
    assert ('$PY experiments/builds/t1_r51_check.py msb t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12'
            in msb)
    assert 'checked "$d" "$name"' in msb
    build = sh[sh.index("build() {"):sh.index("build_msb() {")]
    assert ('$PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11'
            in build)
    assert 'checked "$d" "$name"' in build
    checked = sh[sh.index("checked() {"):sh.index("build() {")]
    assert 'vec check "$d/$name.h5ad" --board "$BOARD" > /dev/null 2>&1 || return 13' in checked
    assert 'st=$(stamp "$d/$name.h5ad" 2> /dev/null) || return 14' in checked


def test_the_rows_step_checks_the_files_the_launcher_built():
    sh = _sh()
    rows = sh[sh.index("STEP=rows"):sh.index("STEP=done")]
    files = re.findall(r'"\$O/([A-Za-z0-9]+)/([A-Za-z0-9_]+)\.h5ad"', rows)
    assert files == [("MSB", "T1_s2rcMSB"), ("H05a", "T1_s2rcMSBH05"), ("H10a", "T1_s2rcMSBH10")]  # MSB's first
    assert '> "$S/rows.out" 2> /dev/null || finish 48' in rows and "t1_r51_check.py rows" in rows
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import t1_r51_check

    assert t1_r51_check.MSB_ARGV == ["t1", *_args(sh)["MSB"].split()]  # the rows check's draw is MSB's
    assert {f"{f}_msb" for f in t1_r51_check.SETTING} == {name for name, _, _ in MSBS.values()}


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    lines = [line for line in sh.splitlines() if not line.lstrip().startswith("#")]
    for line in lines:  # every builder, checker build and check call: output to /dev/null
        if re.search(r'build_test\.py t1|t1_r51_check\.py msb|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    checks = [i for i, line in enumerate(lines) if re.search(r"t1_r51_check\.py rows", line)]
    assert len(checks) == 1
    for i in checks:  # the rows check prints booleans only, to the scratch folder; its errors are dropped
        assert lines[i].endswith("\\") and re.search(r'> "\$S/rows\.out" 2> /dev/null \|\| finish 48$', lines[i + 1])
    code = _code(sh)
    for word in (".build.json", "du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -",
                 "tail ", "cat "):
        assert word not in code, word
    assert "tests.out" in sh  # the tests' log: synthetic stages only
    assert 'print(s["commit"], s["clean"])' in sh and "stamp.json" in sh  # of the stamp, its commit and clean only
    assert "curl" not in sh and "wget" not in sh and "git fetch" not in sh and "git pull" not in sh  # fetches nothing


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    for step in ("tests", *KNOWN, *MSBS, "rows", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    for a, b in PAIRS:
        assert f"for STEP in {a} {b}; do begin; build" in sh
    for code in ("exit 90", "exit 89", "exit 91"):
        assert sh.index(code) < sh.index("STEP=gate"), code
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", sh)}
    used |= {code for *_, code in KNOWN.values()} | {code for *_, code in MSBS.values()}
    used |= {int(c) for c in re.findall(r'ruling_row "\$ROW" (\d+)', sh)}
    documented = {int(c) for c in re.findall(r"(?<![\d.])(\d\d)(?= )", doc)}
    assert used == documented, (used, documented)


def test_the_nineteen_known_hashes_are_the_recorded_ones():
    """The nineteen identity hashes are the ones the earlier runs' committed lists hold (no build here): every one in
    t1_l67pb_line.sh's run list, the pairs there built twice with one sha256, and the fifteen in t1_hm_line.sh."""
    sh, hm = _sh(), (REPO / "experiments" / "builds" / "t1_hm_line.sh").read_text()
    l67 = _sums("t1l67pb_ce1fb91_l1.SHA256SUMS")
    for step, (name, key, args, var, sha, _) in KNOWN.items():
        folders = [step] if f"{step}/{name}.h5ad" in l67 else [f"{step}a", f"{step}b"]
        assert all(l67[f"{f}/{name}.h5ad"] == sha for f in folders), step
        assert f'{key}_ARGS="{args}"' in sh
        if step in KNOWN15:
            assert f"{var}={sha}" in hm and f'{key}_ARGS="{args}"' in hm, step
    assert _sums("t1msc_a3e879a_m1.SHA256SUMS")["MSBa/T1_s2rcMSB.h5ad"] == MSB_SHA
    assert MSB_SHA == "78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2"


ROWS_ALL = ("# ledger\n"
            "| R20 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.0` `--captain-ruling R20`) |\n"
            "| R22 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` "
            "`--g-mom-damping 1.5` `--g-mom-damping 2.0` `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
            "`--g-mom-damping 1.0 --feng-fraction-mult 2.5` `--g-mom-damping 1.5 --feng-fraction-mult 2.0` "
            "`--g-mom-damping 2.0 --feng-fraction-mult 2.0` `--g-mom-damping 1.0 --feng-fraction-mult 3.0`, "
            "`--g-mom-damping 1.5 --feng-fraction-mult 2.5` and `--g-mom-damping 2.0 --feng-fraction-mult 2.5`, each "
            "only with `--captain-ruling R22`) |\n"
            "| R23 (the captain's order of X; ruling, no run; `--direction-momentum all` `--captain-ruling R23`) |\n"
            "| R26 (the captain's order of X; ruling, no run; `--feng-fraction-mult 1.5 --direction-momentum all` "
            "`--captain-ruling R26`) |\n"
            "| R37 (the captain's named exception of X; ruling, no run; `--base pooled` with `--feng-groups "
            "recomputed` and `--captain-ruling R37`) |\n")
ROW_R51 = ("| R51 (the captain's order of a later day, ruling, no run: \"Yes, the seven files\"; T1's line 10; the T1 "
           "builder is to take " + ", ".join(f"`{o}`" for o in OPTS.values()) + f" on `t1 {BASE}`, each only with "
           "`--captain-ruling R51`; not a try) | T1 | line 10 |\n")


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    sys.path.insert(0, str(REPO / "experiments" / "external"))
    import build_test as bt
    import feng_relabel as fr

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with the ruling rows' orders, not the repository's
    ledger.write_text(ROWS_ALL + ROW_R51)
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    args = _args(_sh())
    want = {"BASE": None, "A": (1.0, 1.0, None), "B": (0.5, 1.5, None), "C": (1.0, 1.5, None),
            "X2": (0.5, 2.0, "R20"), "X25": (0.5, 2.5, "R22"), "A15": (1.5, 1.0, "R22"), "A20": (2.0, 1.0, "R22"),
            "CM20": (1.0, 2.0, "R22"), "MS": None, "MSB": (0.5, 1.5, None), "X30": (0.5, 3.0, "R22"),
            "CM25": (1.0, 2.5, "R22"), "B2D15": (1.5, 2.0, "R22"), "B2D20": (2.0, 2.0, "R22"),
            "CM30": (1.0, 3.0, "R22"), "B25D15": (1.5, 2.5, "R22"), "B25D20": (2.0, 2.5, "R22"), "POOLED": None}
    assert set(want) | set(NEW) == set(args)
    for key, values in want.items():
        a = bt.build_parser().parse_args(["t1", *args[key].split(), "--out", "x.h5ad"])
        assert bt._t1_r51(a) is None  # no earlier file takes row R51's route
        cand = bt._t1_candidate(a)
        assert (None if cand is None else (cand["g_mom_damping"], cand["feng_fraction_mult"],
                                           cand.get("captain_ruling"))) == values, key
    for f, (key, _) in KEYS.items():
        a = bt.build_parser().parse_args(["t1", *args[key].split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        r51 = bt._t1_r51(a)  # nothing read but the ledger
        assert r51["file"] == f and r51["option"] == f"`{OPTS[f]}`" and r51["captain_ruling"] == "R51"
        assert bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None and bt._t1_pooled(a) is None


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    row = gate[gate.index("ruling_row() {"):gate.index("\n}\n")]
    assert "[ \"$(grep -c \"^| $id \" experiments/harness/TRIES.md)\" = 1 ]" in row
    assert "first=$(grep \"^| $id \" experiments/harness/TRIES.md | cut -d'|' -f2)" in row
    assert "grep -q \"^ $id (the captain's \"" in row and 'grep -qF -- "ruling, no run"' in row
    assert 'for opt in "$@" "--captain-ruling $id"; do' in row
    assert "grep -o '`[^`]*`' | grep -qxF -- \"\\`$opt\\`\"" in row  # whole backticked strings, line by line
    assert row.count('finish "$code"') == 4
    named = re.findall(r'"(--feng-fraction-mult [^"]+)"', gate[gate.index('ruling_row "$ROW" 95'):])
    assert named == list(OPTS.values())
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    assert "O=runs/t1r51_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1r51/$RUN" in sh
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 92, 88, 93]
    assert '[ "$ROW" = R51 ] || {' in sh


def test_the_tests_step_runs_the_builders_and_launchers_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    hm = (REPO / "experiments" / "builds" / "t1_hm_line.sh").read_text()
    hm_step = hm[hm.index("STEP=tests"):hm.index("step_ok", hm.index("STEP=tests"))]
    assert {"tests/test_t1_r51_launcher.py", "tests/test_t1_r51_build.py",
            *re.findall(r"tests/test_[a-z0-9_]+\.py", hm_step)} == set(listed)
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


# ------------------------------------------------------------------ the gate, run for real


RULES = "the amendment names the rows\n"


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    env["UV_PROJECT_ENVIRONMENT"] = str(tmp_path / "elsewhere")  # the launcher unsets it
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R51, relabel: bool = True, good_sha: bool = True,
          venv: str | None = "own", launcher: str | None = None) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (or the given copy of it) and the files its gate reads (with a
    .venv whose vec.pth names its src, or another checkout's src, or none), and a home whose ~/.local/bin holds a `uv`
    that fails (the launcher puts that folder first on PATH): returns the repository, the env and its HEAD."""
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
    (repo / "experiments" / "builds" / "t1_r51_line.sh").write_text(_sh() if launcher is None else launcher)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / "docs" / "S2_RULES.md").write_text(RULES)
    sha = hashlib.sha256(RULES.encode()).hexdigest() if good_sha else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    (repo / ".gitignore").write_text("data/\nruns/\n.venv/\n")
    if relabel:
        (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
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
    return subprocess.run(["bash", "experiments/builds/t1_r51_line.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1r51" / run


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no head file, no later step's mark
    assert not list((s / "tmp").iterdir())


def test_a_call_that_passes_the_gate_writes_the_head_and_stops_at_the_failing_tests_step(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1", "R51")
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t1r51_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "R51")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2", "R51")
    assert other.returncode == 10 and (repo / "runs" / f"t1r51_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_and_a_prefix_are_refused_and_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1", "R51")
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()
    again = _run(repo, env, head, "r1", "R51")
    assert again.returncode == 91 and "never reused" in again.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    short = _run(repo, env, head[:12], "r2", "R51")
    assert short.returncode == 96
    _stopped_at_the_gate(tmp_path, "r2", 96)


@pytest.mark.parametrize("argv, code", [
    ([], None), (["r1"], None),  # missing arguments
    (["a/b", "R51"], 90), (["r 1", "R51"], 90),
    *[(["r1", bad], 89) for bad in ("R44", "R26", "r51", "R510", "R51a", "R5", "51", "R51 ", "$(id)", "R51`")],
])
def test_a_missing_or_odd_run_name_or_row_id_is_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, *argv)
    assert out.returncode != 0 if code is None else out.returncode == code, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R51}more\n")
    out = _run(repo, env, head, "r1", "R51")
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


FIRST = "| R51 (the captain's order of X; ruling, no run; {opts} `--captain-ruling R51`) | T1 |\n"
ALL = " ".join(f"`{o}`" for o in OPTS.values())


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ""}, 95, "exactly one row R51"),
    ({"tries": ROW_R51 + ROW_R51}, 95, "exactly one row R51"),  # two rows with one ID
    ({"tries": ROW_R51.replace("| R51 (", "| R51b (")}, 95, "exactly one row R51"),  # another ID
    ({"tries": ROW_R51.replace("R51 (the captain's order", "R51 (a proposal")}, 95, "not the captain's ruling"),
    ({"tries": ROW_R51.replace("R51 (the captain's", "R51 awaits (the captain's")}, 95, "not the captain's ruling"),
    ({"tries": ROW_R51.replace("ruling, no run", "a try")}, 95, "does not say ruling, no run"),
    *[({"tries": FIRST.format(opts=ALL.replace(f"`{o}`", ""))}, 95, f"does not name `{o}`") for o in OPTS.values()],
    # an option inside a longer backticked string, or with no backticks, is not named
    *[({"tries": FIRST.format(opts=ALL.replace(f"`{o}`", f"`{o} --g-mom-damping 2.0`"))}, 95, f"does not name `{o}`")
      for o in OPTS.values()],
    ({"tries": FIRST.format(opts=ALL.replace("`", ""))}, 95, f"does not name `{OPTS['T1_s2rcMS20']}`"),
    ({"tries": FIRST.format(opts=ALL).replace("`--captain-ruling R51`", "`--captain-ruling R22`")}, 95,
     "does not name `--captain-ruling R51`"),
    ({"tries": FIRST.format(opts=ALL).replace("`--captain-ruling R51`", "`--captain-ruling R51 now`")}, 95,
     "does not name `--captain-ruling R51`"),
    ({"tries": f"| R51 (the captain's order of X; ruling, no run) | T1 | {ALL} `--captain-ruling R51` |\n"}, 95,
     f"does not name `{OPTS['T1_s2rcMS20']}`"),  # named outside the first cell only
    ({"relabel": False}, 92, "feng_relabel_c is not in this checkout"),
    ({"venv": None}, 88, "no .venv of its own"),
    ({"venv": "other"}, 88, "no .venv of its own"),  # a venv whose vec.pth names another checkout's src
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1", "R51")
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"t1r51_{head[:7]}_r1").exists()


def test_the_gate_takes_the_row_in_its_own_words_and_in_any_order(tmp_path):
    tries = FIRST.format(opts=" and ".join(f"`{o}`" for o in reversed(OPTS.values())))  # the names in another order
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, head, "r1", "R51")
    assert out.returncode == 10, out.stderr  # past the gate: the stand-in uv fails the tests step
    assert (_scratch(tmp_path, "r1") / "gate.exit").read_text().strip() == "0"


def test_a_call_past_the_tests_step_stops_at_the_first_unknown_bytes(tmp_path):
    """A stand-in `uv` lets the call run past the tests step (t1_l67pb_line.sh's tests' stand-in: the tests pass, a
    build writes a stand-in file whose stamp names HEAD and clean, `vec check` passes); no stand-in file has a known
    file's bytes, so the run stops at step D."""
    repo, env, head = _repo(tmp_path)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(PASSING_UV)
    out = _run(repo, env, head, "r1", "R51")
    assert out.returncode == 20, out.stderr
    assert "stopped: step D does not give the known bytes of T1_s2rc" in out.stderr
    s = _scratch(tmp_path, "r1")
    codes = {p.name: p.read_text().strip() for p in s.glob("*.exit")}
    assert codes == {"gate.exit": "0", "tests.exit": "0", "D.exit": "20"}
    assert (s / "stamps").read_text() == f"D/T1_s2rc {head} True\n"
    o = repo / "runs" / f"t1r51_{head[:7]}_r1"
    assert sorted(p.name for p in o.iterdir()) == ["D"]
    assert (o / "D" / "T1_s2rc.h5ad").read_text() == (f"stand-in experiments/builds/build_test.py t1 {BASE} --out "
                                                      f"runs/t1r51_{head[:7]}_r1/D/T1_s2rc.h5ad\n")


# A stand-in `uv` for a whole run: the tests pass; a build writes a file named after its output (the checker's msb
# writes MSB's, as each new path at MSB's values must) with a stamp naming HEAD and clean; `vec check` passes; the rows
# check prints one True; the stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/build_test.py|experiments/builds/t1_r51_check.py)
        prog=$1; shift
        if [ "$prog" = experiments/builds/t1_r51_check.py ]; then
          case "$1" in
            rows) echo "$* stand_in True" >> calls.txt; echo "stand_in True"; exit 0 ;;
          esac
          echo "$1 $2" >> calls.txt
          shift
        fi
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        name=$(basename "$out" .h5ad)
        [ "$prog" = experiments/builds/t1_r51_check.py ] && name=T1_s2rcMSB
        echo "file $name" > "$out"
        printf '{"commit": "%s", "clean": true}' "$(git rev-parse HEAD)" > "$out.stamp.json"
        exit 0 ;;
    esac ;;
esac
exit 1
"""


def test_a_whole_stand_in_run_passes_every_step_in_order(tmp_path):
    """The launcher with only its nineteen known hashes set to the stand-in files' (in a copy) and the stand-in `uv`
    above: a call runs every step in order to done.exit 0, each step's exit file says 0, the stamps list every build,
    the checker gets the files its steps name, the hash list holds every built file once and each pair is compared."""
    sh = _sh()
    for name, _, _, var, sha, _ in KNOWN.values():
        content = f"file {name}\n".encode()
        sh = sh.replace(f"{var}={sha}", f"{var}={hashlib.sha256(content).hexdigest()}")
    assert sh != _sh()
    repo, env, head = _repo(tmp_path, launcher=sh)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(FULL_UV)
    out = _run(repo, env, head, "r1", "R51")
    assert out.returncode == 0, out.stderr
    s = _scratch(tmp_path, "r1")
    codes = {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")}
    assert codes == {step: "0" for step in STEPS}
    o = repo / "runs" / f"t1r51_{head[:7]}_r1"
    folders = STEPS[2:-2]
    assert sorted(p.name for p in o.iterdir()) == sorted([*folders, "SHA256SUMS"])
    stamps = (s / "stamps").read_text().splitlines()
    assert [line.split()[0].split("/")[0] for line in stamps] == folders
    assert all(line.endswith(f" {head} True") for line in stamps)
    listed = (o / "SHA256SUMS").read_text().splitlines()
    assert [line.split()[1].split("/")[0] for line in listed] == folders and (s / "SHA256SUMS").exists()
    msb = hashlib.sha256(b"file T1_s2rcMSB\n").hexdigest()
    for step in MSBS:
        assert next(line for line in listed if line.split()[1].startswith(f"{step}/")).split()[0] == msb
    calls = (repo / "calls.txt").read_text().splitlines()
    r = f"runs/t1r51_{head[:7]}_r1"
    assert calls == [*(["msb t1"] * len(MSBS)),
                     f"rows {r}/MSB/T1_s2rcMSB.h5ad {r}/H05a/T1_s2rcMSBH05.h5ad {r}/H10a/T1_s2rcMSBH10.h5ad "
                     "stand_in True"]
    assert (s / "rows.out").read_text() == "stand_in True\n"
    m10 = (o / "m10a" / "T1_s2rcMSBm10.h5ad").read_text()
    assert m10 == "file T1_s2rcMSBm10\n" == (o / "m10b" / "T1_s2rcMSBm10.h5ad").read_text()


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t1r51_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", "R51")
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]
