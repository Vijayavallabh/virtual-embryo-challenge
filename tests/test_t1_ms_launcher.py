"""The test-board build launcher of candidate MS under row R23 (experiments/harness/TRIES.md; the captain's order of
2026-11-23, docs/PLAN.md decision log 2026-11-23 (c)): R13's file with B0's released per-category damped mean momentum
added on the cells whose label's group has a recomputed direction, before the atlas step and g_mom.

experiments/builds/t1_ms.sh first rebuilds R13's file (the default command) and stops unless it gives its known bytes;
then it builds the MS file (build_test.py with --direction-momentum all --captain-ruling R23) twice from one clean
commit and compares the pair. These tests hold its command strings, hash and step order to that, check the hash against
the committed build records, check that it deletes nothing and reads no output of a test build, and run its gate's
paths for real in throwaway repositories under a temporary home, with a stand-in `uv` that fails, so that a call
passing the gate stops at the tests step and builds nothing. Nothing here reads the repository's own ledger or rules
file for its state of today."""

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_ms.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
R13_SHA = "ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f"  # set t1rc (row R13)
STEPS = ["gate", "tests", "D", "MSa", "MSb", "done"]
OPT = "`--direction-momentum all`"
RULING = "`--captain-ruling R23`"


def _sh() -> str:
    return LAUNCHER.read_text()


# ------------------------------------------------------------------ what the launcher runs


def test_the_launcher_reproduces_r13s_file_then_builds_the_ms_file_twice():
    sh = _sh()
    args = dict(re.findall(r'^(BASE|MS)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"BASE": BASE, "MS": f"{BASE} --direction-momentum all --captain-ruling R23"}
    assert "build D T1_s2rc $BASE_ARGS || finish $?" in sh
    assert 'got=$(sha256sum "$O/D/T1_s2rc.h5ad" | cut -d\' \' -f1)' in sh
    assert f"R13_SHA={R13_SHA}" in sh
    line = next(ln for ln in sh.splitlines() if '[ "$got" = "$R13_SHA" ]' in ln)
    assert line.endswith("finish 20; }"), line
    assert "for STEP in MSa MSb; do begin; build \"$STEP\" T1_s2rc_MS $MS_ARGS || finish $?; step_ok; done" in sh
    order = [sh.index(s) for s in ("STEP=gate", "STEP=tests", "STEP=D\n", "for STEP in MSa MSb", "STEP=done")]
    assert order == sorted(order)
    assert 'cmp -s "$O/MSa/T1_s2rc_MS.h5ad" "$O/MSb/T1_s2rc_MS.h5ad" || finish 31' in sh
    listing = sh[sh.index("(cd \"$O\" && sha256sum"):sh.index('> "$S/SHA256SUMS"')]
    assert re.findall(r"([A-Za-z0-9]+)/\*\.h5ad", listing) == STEPS[2:-1]  # every built folder, once, in order
    assert 'cp "$S/SHA256SUMS" "$O/SHA256SUMS"' in sh
    assert "vec check" in sh and "BOARD=T1:test" in sh
    assert "--seed" not in sh and "--board T1" not in sh  # seed and board: the builder's defaults, as the references
    for other in ("--g-mom-damping", "--feng-fraction-mult", "--released-bet", "--bet-strength",
                  "--no-hold-unbracketed"):
        assert other not in sh  # MS is never stacked on A, B or C
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh  # the builds, as the references
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh  # the tests
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=D\n")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert sh.index("PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")  # PATH set before the first uv call
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
    for step in ("tests", "D", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    assert "for STEP in MSa MSb; do begin; build" in sh
    assert sh.index("exit 91") < sh.index("STEP=gate") and sh.index("exit 90") < sh.index("STEP=gate")
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


def test_the_identity_hash_is_the_recorded_one():
    """R13's hash is the one the earlier runs' committed lists hold (no build here)."""
    assert _sums("t1rc_5028cab.SHA256SUMS")["builds/T1_s2rc.h5ad"] == R13_SHA
    assert _sums("t1x2_3d18841.SHA256SUMS")["D/T1_s2rc.h5ad"] == R13_SHA
    assert _sh().count(R13_SHA) == 1


def test_the_builder_takes_the_command_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test as bt

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with row R23's order, not the repository's
    ledger.write_text(f"# ledger\n| R23 (the captain's order of X; ruling, no run; {OPT} only with {RULING}) |\n")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    assert bt.T1_R23_OPTION == OPT
    args = dict(re.findall(r'^(BASE|MS)_ARGS="([^"]+)"$', _sh(), flags=re.M))
    a = bt.build_parser().parse_args(["t1", *args["MS"].split(), "--out", "x.h5ad"])
    assert a.captain_ruling == "R23" and a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
    ms = bt._t1_ms(a)  # the refusals pass: nothing is read but the ledger
    assert (ms["direction_momentum"], ms["captain_ruling"]) == ("all", "R23")
    assert bt._t1_candidate(a) is None and bt._bet_refusals(a) is None
    assert " ".join(["t1", *args["MS"].split()]) == bt.T1_MS_COMMAND
    d = bt.build_parser().parse_args(["t1", *args["BASE"].split(), "--out", "x.h5ad"])
    assert bt._t1_ms(d) is None and bt._t1_candidate(d) is None


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    assert "[ \"$(grep -c '^| R23 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "first=$(grep '^| R23 ' experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert f"for opt in '{OPT}' '{RULING}'; do" in gate
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "O=runs/t1ms_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1ms/$RUN" in sh
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 92, 93]


def test_the_tests_step_runs_the_builders_and_launchers_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    assert set(listed) == {"tests/test_t1_ms_launcher.py", "tests/test_t1_ms_build.py", "tests/test_t65_ms_read.py",
                           "tests/test_t1_line3_build.py", "tests/test_t1_damping_line_build.py",
                           "tests/test_t1_fraction_line_build.py", "tests/test_t1_fraction2_build.py",
                           "tests/test_t1_candidates_build.py", "tests/test_t1_bets.py", "tests/test_hidden_counts.py",
                           "tests/test_s2_rules.py"}
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in _sh().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 2
    for msg in usages:  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


# ------------------------------------------------------------------ the gate, run for real


ROW_R23 = f"| R23 (the captain's order of a later day; ruling, no run; {OPT} only with {RULING}) | T1 | MS |\n"


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R23, rules: str = "the amendment names row R23\n",
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
    (repo / "experiments" / "builds" / "t1_ms.sh").write_text(_sh())
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
    return subprocess.run(["bash", "experiments/builds/t1_ms.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1ms" / run


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
    o = repo / "runs" / f"t1ms_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2")
    assert other.returncode == 10 and (repo / "runs" / f"t1ms_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_and_a_prefix_are_refused_and_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1")
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()  # refused before any output folder
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
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R23}more\n")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": "| R22 (another row) |\n"}, 95, "exactly one row R23"),
    ({"tries": f"| R23b (a row with another ID; {OPT} {RULING}) |\n"}, 95, "exactly one row R23"),
    ({"tries": ROW_R23 + ROW_R23}, 95, "exactly one row R23"),  # two rows with one ID
    ({"tries": f"| R23 (the captain's order; ruling, no run; only with {RULING}) |\n"}, 95,
     f"first cell does not name {OPT}"),
    ({"tries": "| R23 (the captain's order; --direction-momentum all with `--captain-ruling R23`) |\n"}, 95,
     f"first cell does not name {OPT}"),  # not in backticks
    ({"tries": f"| R23 (the captain's order; {OPT}) |\n"}, 95, f"first cell does not name {RULING}"),
    ({"tries": f"| R23 (the captain's order; {OPT} with `--captain-ruling R22`) |\n"}, 95,
     f"first cell does not name {RULING}"),
    ({"tries": f"| R23 (the captain's order) | T1 | {OPT} {RULING} |\n"}, 95,
     f"first cell does not name {OPT}"),  # named outside the first cell only
    ({"tries": f"| R23 (the captain's order) |\n| R24 ({OPT} {RULING}) |\n"}, 95,
     f"first cell does not name {OPT}"),  # named in another row only
    ({"relabel": False}, 92, "feng_relabel_c is not in this checkout"),
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1")
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"t1ms_{head[:7]}_r1").exists()


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t1ms_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]
