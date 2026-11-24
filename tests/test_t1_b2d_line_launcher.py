"""The test-board build launcher of row R22's fifth T1 search line (experiments/harness/TRIES.md; T1 line 5 of
docs/PLAN.md's decision log): g_mom's within-type damping at 1.5 and at 2.0 on the 2.0 file's command (the outside
atlas's direction step at 2.0 x the calendar fraction), nothing else changed.

experiments/builds/t1_b2d_line.sh first rebuilds eleven known files (R13's, A, B, C, the 2.0 file, T1_s2rcB25,
T1_s2rcA15, T1_s2rcA20, T1_s2rcCm20, MS and T1_s2rcMSB) and stops unless each gives its known bytes; then it builds
T1_s2rcB2d15 and T1_s2rcB2d20 (build_test.py with --feng-fraction-mult 2.0 and --captain-ruling R22) twice each from
one clean commit and compares each pair. These tests hold its command strings, hashes and step order to that, check the
eleven hashes against the committed build records, check that it deletes nothing and reads no output of a test build,
and run its gate's paths for real in throwaway repositories under a temporary home, with a stand-in `uv` that fails, so
that a call passing the gate stops at the tests step and builds nothing. Nothing here reads the repository's own ledger
or rules file for its state of today."""

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_b2d_line.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
# step -> (file name, args variable, args after BASE, hash variable, sha256, exit code)
KNOWN = {
    "D": ("T1_s2rc", "BASE", "", "R13_SHA", "ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f", 20),
    "A": ("T1_s2rcA", "A", " --g-mom-damping 1.0", "A_SHA",
          "76865826f5a3741c7585f81c4241048ca1a95b5791ebb3c8feb1f2b536678ce4", 21),
    "B": ("T1_s2rcB", "B", " --feng-fraction-mult 1.5", "B_SHA",
          "c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062", 22),
    "C": ("T1_s2rcC", "C", " --g-mom-damping 1.0 --feng-fraction-mult 1.5", "C_SHA",
          "31477ef7df5d5374a74a8249c570217020c593916ab4d3780049b7b3f87b0620", 23),
    "X2": ("T1_s2rcB2", "X2", " --feng-fraction-mult 2.0 --captain-ruling R20", "X2_SHA",
           "52897fe375b8a81afffa788592047c30361824148b00a8176514385060518837", 24),
    "X25": ("T1_s2rcB25", "X25", " --feng-fraction-mult 2.5 --captain-ruling R22", "X25_SHA",
            "ccb7c457c0e8ac8ac668bcea17bc61ed22f414f58ab5a63e946cb2f08c4f4e77", 25),
    "A15": ("T1_s2rcA15", "A15", " --g-mom-damping 1.5 --captain-ruling R22", "A15_SHA",
            "e2334137b54d98b04e9dece90213987e61fe9aa15a3fb8d85e6c262d8fb5ba75", 26),
    "A20": ("T1_s2rcA20", "A20", " --g-mom-damping 2.0 --captain-ruling R22", "A20_SHA",
            "372713686782be273a52427140eab652517e35b679f54fc2a969484dee730218", 27),
    "Cm20": ("T1_s2rcCm20", "CM20", " --g-mom-damping 1.0 --feng-fraction-mult 2.0 --captain-ruling R22", "CM20_SHA",
             "acce5f5f7097269bffebbfd583833e4c9c1c5ae77f1d5294dde6d13e197f863e", 28),
    "MS": ("T1_s2rc_MS", "MS", " --direction-momentum all --captain-ruling R23", "MS_SHA",
           "6103d10535665634ee0cad3390590a5df132d1805682fb1d309bb19581fdaa5e", 29),
    "MSB": ("T1_s2rcMSB", "MSB", " --feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26", "MSB_SHA",
            "78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2", 30),
}
STEPS = ["gate", "tests", *KNOWN, "B2d15a", "B2d15b", "B2d20a", "B2d20b", "done"]
OPT15 = "`--g-mom-damping 1.5 --feng-fraction-mult 2.0`"
OPT20 = "`--g-mom-damping 2.0 --feng-fraction-mult 2.0`"
NEW = {"B2D15": f"{BASE} --g-mom-damping 1.5 --feng-fraction-mult 2.0 --captain-ruling R22",
       "B2D20": f"{BASE} --g-mom-damping 2.0 --feng-fraction-mult 2.0 --captain-ruling R22"}


def _sh() -> str:
    return LAUNCHER.read_text()


def _args(sh: str) -> dict[str, str]:
    return dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launcher_reproduces_eleven_files_then_builds_each_line_file_twice():
    sh = _sh()
    want = {key: f"{BASE}{tail}" for _, key, tail, _, _, _ in KNOWN.values()} | NEW
    assert _args(sh) == want
    for step, (name, key, _, var, sha, code) in KNOWN.items():
        assert f"build {step} {name} ${key}_ARGS || finish $?" in sh
        assert f'same {step} {name} "${var}" {code}' in sh
        assert f"{var}={sha}" in sh and sh.count(sha) == 1
    assert 'got=$(sha256sum "$O/$1/$2.h5ad" | cut -d\' \' -f1)' in sh
    assert '[ "$got" = "$3" ] || {' in sh and 'finish "$4"; }' in sh
    assert 'build "$STEP" T1_s2rcB2d15 $B2D15_ARGS' in sh and "for STEP in B2d15a B2d15b" in sh
    assert 'build "$STEP" T1_s2rcB2d20 $B2D20_ARGS' in sh and "for STEP in B2d20a B2d20b" in sh
    order = [sh.index(s) for s in ("STEP=gate", "STEP=tests", *(f"STEP={k}\n" for k in KNOWN),
                                   "for STEP in B2d15a B2d15b", "for STEP in B2d20a B2d20b", "STEP=done")]
    assert order == sorted(order)
    assert 'cmp -s "$O/B2d15a/T1_s2rcB2d15.h5ad" "$O/B2d15b/T1_s2rcB2d15.h5ad" || finish 41' in sh
    assert 'cmp -s "$O/B2d20a/T1_s2rcB2d20.h5ad" "$O/B2d20b/T1_s2rcB2d20.h5ad" || finish 42' in sh
    listing = sh[sh.index("(cd \"$O\" && sha256sum"):sh.index('> "$S/SHA256SUMS"')]
    assert re.findall(r"([A-Za-z0-9]+)/\*\.h5ad", listing) == STEPS[2:-1]  # every built folder, once, in order
    assert 'cp "$S/SHA256SUMS" "$O/SHA256SUMS"' in sh
    assert "vec check" in sh and "BOARD=T1:test" in sh
    assert "--seed" not in sh and "--board T1" not in sh  # seed and board: the builder's defaults, as the references
    assert "--feng-fraction-mult 3" not in sh and "--g-mom-damping 2.5" not in sh and "--g-mom-damping 0" not in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh  # the builds, as the references
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh  # the tests
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=D\n")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert sh.index("PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")  # PATH set before the first uv call
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh  # the full hash, never a prefix
    assert '[ "$st" = "$head True" ] || return 15' in sh  # a build counts only with a stamp of this commit, clean
    for line in sh.splitlines():  # every builder and check call: output to /dev/null, read by its exit code only
        if re.search(r'build_test\.py t1|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    assert ".build.json" not in sh and "tests.out" in sh  # no sidecar record is read
    assert "upload" not in sh.replace("Nothing here uploads anything", "")
    assert "/v/vec-campaign" not in sh.replace("(~/v/vec-campaign-d, never the reference checkout)", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")  # runs from the folder it is called in
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")  # the rules file is checked by its hash only


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    for step in ("tests", *KNOWN, "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    assert "for STEP in B2d15a B2d15b; do begin; build" in sh and "for STEP in B2d20a B2d20b; do begin; build" in sh
    assert sh.index("exit 91") < sh.index("STEP=gate") and sh.index("exit 90") < sh.index("STEP=gate")
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", sh)}
    used |= {code for *_, code in KNOWN.values()}
    documented = {int(c) for c in re.findall(r"(?<![\d.])(\d\d)(?= )", doc)}
    assert used == documented, (used, documented)


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


def test_the_eleven_known_hashes_are_the_recorded_ones():
    """The eleven identity hashes are the ones the earlier runs' committed lists hold (no build here)."""
    sha = {step: v[4] for step, v in KNOWN.items()}
    cand = _sums("t1cand_1420c70.SHA256SUMS")
    assert (cand["T1_s2rcA.h5ad"], cand["T1_s2rcB.h5ad"], cand["T1_s2rcC.h5ad"]) == (sha["A"], sha["B"], sha["C"])
    x2 = _sums("t1x2_3d18841.SHA256SUMS")
    assert x2["D/T1_s2rc.h5ad"] == sha["D"] and x2["X2a/T1_s2rcB2.h5ad"] == x2["X2b/T1_s2rcB2.h5ad"] == sha["X2"]
    line1 = _sums("t1line_e3cc0e4_r1.SHA256SUMS")
    assert line1["X25a/T1_s2rcB25.h5ad"] == line1["X25b/T1_s2rcB25.h5ad"] == sha["X25"]
    damp = _sums("t1damp_9a10157_r1.SHA256SUMS")
    assert damp["X15a/T1_s2rcA15.h5ad"] == damp["X15b/T1_s2rcA15.h5ad"] == sha["A15"]
    assert damp["X20a/T1_s2rcA20.h5ad"] == damp["X20b/T1_s2rcA20.h5ad"] == sha["A20"]
    cline = _sums("t1cline_2a1c248_c1.SHA256SUMS")
    assert cline["Cm20a/T1_s2rcCm20.h5ad"] == cline["Cm20b/T1_s2rcCm20.h5ad"] == sha["Cm20"]
    ms = _sums("t1ms_6bc36d9_r1.SHA256SUMS")
    assert ms["MSa/T1_s2rc_MS.h5ad"] == ms["MSb/T1_s2rc_MS.h5ad"] == sha["MS"]
    msc = _sums("t1msc_a3e879a_m1.SHA256SUMS")
    assert msc["MSBa/T1_s2rcMSB.h5ad"] == msc["MSBb/T1_s2rcMSB.h5ad"] == sha["MSB"]


def test_the_builder_takes_the_line_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with the ruling rows' orders, not the repository's
    ledger.write_text("# ledger\n| R20 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.0` "
                      "`--captain-ruling R20`) |\n"
                      "| R22 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.5` "
                      "`--g-mom-damping 1.5` `--g-mom-damping 2.0` `--g-mom-damping 1.0 --feng-fraction-mult 2.0` "
                      f"{OPT15} and {OPT20}, each only with `--captain-ruling R22`) |\n")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    assert bt.T1_R22_LINE5_OPTIONS == {1.5: OPT15, 2.0: OPT20}
    args = _args(_sh())
    want = {"A": (1.0, 1.0, None), "B": (0.5, 1.5, None), "C": (1.0, 1.5, None), "X2": (0.5, 2.0, "R20"),
            "X25": (0.5, 2.5, "R22"), "A15": (1.5, 1.0, "R22"), "A20": (2.0, 1.0, "R22"), "CM20": (1.0, 2.0, "R22"),
            "B2D15": (1.5, 2.0, "R22"), "B2D20": (2.0, 2.0, "R22")}
    for key, (damping, mult, row) in want.items():
        a = bt.build_parser().parse_args(["t1", *args[key].split(), "--out", "x.h5ad"])
        assert a.captain_ruling == row and a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        cand = bt._t1_candidate(a)  # the refusals pass: nothing is read but the ledger
        assert (cand["g_mom_damping"], cand["feng_fraction_mult"], cand.get("captain_ruling")) == (damping, mult, row)
        if key.startswith("B2D"):
            assert "fifth T1 search line" in cand["order"]


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    assert "[ \"$(grep -c '^| R22 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "first=$(grep '^| R22 ' experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert f"for opt in '{OPT15}' '{OPT20}' \\\n  '`--captain-ruling R22`'; do" in gate
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "O=runs/t1b2dline_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1b2dline/$RUN" in sh
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 92, 93]


def test_the_tests_step_runs_the_builders_and_launchers_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    assert {"tests/test_t1_b2d_line_launcher.py", "tests/test_t1_line5_build.py", "tests/test_t1_line3_build.py",
            "tests/test_t1_damping_line_build.py", "tests/test_t1_fraction2_build.py", "tests/test_t1_ms_build.py",
            "tests/test_t1_msb_build.py", "tests/test_t1_candidates_build.py", "tests/test_s2_rules.py"} <= set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in _sh().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    assert "rm " not in _sh()
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 2
    for msg in usages:  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


# ------------------------------------------------------------------ the gate, run for real


ROW_R22 = (f"| R22 (the captain's order of a later day; ruling, no run; {OPT15} and {OPT20}, each only with "
           "`--captain-ruling R22`) | T1 | the lines |\n")


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R22, rules: str = "the amendment names row R22\n",
          relabel: bool = True, good_sha: bool = True) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher and the files its gate reads, and a home whose ~/.local/bin holds a
    `uv` that fails (the launcher puts that folder first on PATH): returns the repository, the env and its HEAD."""
    env = _env(tmp_path)
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text("#!/bin/sh\nexit 1\n")
    (fake / "uv").chmod(0o755)
    repo = tmp_path / "repo"
    (repo / "experiments" / "builds").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "docs").mkdir()
    (repo / "experiments" / "builds" / "t1_b2d_line.sh").write_text(_sh())
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / "docs" / "S2_RULES.md").write_text(rules)
    sha = hashlib.sha256(rules.encode()).hexdigest() if good_sha else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    (repo / ".gitignore").write_text("data/\nruns/\n")
    if relabel:
        (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "experiments/builds/t1_b2d_line.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1b2dline" / run


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no head file, no later step's mark
    assert not list((s / "tmp").iterdir())


def test_a_call_that_passes_the_gate_writes_the_head_and_stops_at_the_failing_tests_step(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1")
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t1b2dline_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2")
    assert other.returncode == 10 and (repo / "runs" / f"t1b2dline_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_and_a_prefix_are_refused_and_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1")
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()
    again = _run(repo, env, head, "r1")
    assert again.returncode == 91 and "never reused" in again.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    short = _run(repo, env, head[:12], "r2")
    assert short.returncode == 96
    _stopped_at_the_gate(tmp_path, "r2", 96)


def test_a_missing_or_odd_run_name_is_refused_before_any_folder_is_made(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head).returncode != 0
    assert _run(repo, env).returncode != 0
    odd = _run(repo, env, head, "a/b")
    assert odd.returncode == 90 and not (tmp_path / "home" / "v").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R22}more\n")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": "| R20 (another row) |\n"}, 95, "exactly one row R22"),
    ({"tries": f"| R22b (a row with another ID; {OPT15} {OPT20} `--captain-ruling R22`) |\n"}, 95,
     "exactly one row R22"),
    ({"tries": ROW_R22 + ROW_R22}, 95, "exactly one row R22"),  # two rows with one ID
    ({"tries": "| R22 (the captain's order; `--g-mom-damping 1.5` `--g-mom-damping 2.0` `--feng-fraction-mult 2.0` "
               "`--captain-ruling R22`) |\n"}, 95, f"first cell does not name {OPT15}"),  # the bare options only
    ({"tries": f"| R22 (the captain's order; {OPT15} only, with `--captain-ruling R22`) |\n"}, 95,
     f"first cell does not name {OPT20}"),
    ({"tries": f"| R22 (the captain's order; {OPT15} and {OPT20}) |\n"}, 95,
     "first cell does not name `--captain-ruling R22`"),
    ({"tries": f"| R22 (the captain's order; {OPT15} and {OPT20}, with `--captain-ruling R20`) |\n"}, 95,
     "first cell does not name `--captain-ruling R22`"),
    ({"tries": f"| R22 (the captain's order) | T1 | {OPT15} {OPT20} `--captain-ruling R22` |\n"}, 95,
     f"first cell does not name {OPT15}"),  # named outside the first cell only
    ({"tries": f"| R22 (the captain's order) |\n| R21 ({OPT15} {OPT20} `--captain-ruling R22`) |\n"}, 95,
     f"first cell does not name {OPT15}"),  # named in another row only
    ({"relabel": False}, 92, "feng_relabel_c is not in this checkout"),
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1")
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"t1b2dline_{head[:7]}_r1").exists()


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t1b2dline_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]
