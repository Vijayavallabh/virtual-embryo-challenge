"""The test-board build launcher of row R22's third T3 line's file:
experiments/t3_decomp/t3_c2g_build.sh <COMMIT> <RUN>.

It gates on a clean tree at COMMIT and on rows R22 and R41 of the ledger, runs the targeted tests, rebuilds five known
files (identity-s1, identity-s2, c1-s2, c2-s2 and line 2's c1-s2 file) and stops unless each gives its recorded bytes,
then builds T3_c2_s2_g2.h5ad twice into separate folders (the builder's guard and factor check run inside each build),
compares the pair and writes one SHA256SUMS. These tests hold its strings, hashes and step order to that, check the
hashes against the committed SHA256SUMS records, check that it deletes, fetches and uploads nothing and reads no output
of a test build, and run it for real in throwaway repositories under a temporary home with a stand-in `uv` (it logs
each call and fakes the tests, the builder and `vec check`): every gate refusal, every stop, and a whole run through a
copy of the launcher whose five recorded hashes are replaced by the stand-in's. Nothing here reads the repository's own
ledger or builds anything from released data.
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "t3_decomp" / "t3_c2g_build.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
REL = "experiments/t3_decomp/t3_c2g_build.sh"
DEC = "experiments/t3_decomp/results/gata4/gata4_decisions.json"
WT = "data/raw/T2/E8.75.h5ad"
FILE = "T3_c2_s2_g2.h5ad"
C2G = "--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22"
IDENTITY = [("ids1", "identity-s1.h5ad", "IDS1_SHA", 20, "--form identity-s1"),
            ("ids2", "identity-s2.h5ad", "IDS2_SHA", 21, "--form identity-s2"),
            ("c1s2", "c1-s2.h5ad", "C1S2_SHA", 22, "--form c1-s2"),
            ("c2s2", "c2-s2.h5ad", "C2S2_SHA", 23, "--form c2-s2"),
            ("c1g2", "T3_c1_s2_g2.h5ad", "C1G2_SHA", 24, "--form c1-s2 --gata4-exponent 2.0 --captain-ruling R22")]
STEPS = ["gate", "tests", *[i[0] for i in IDENTITY], "a", "b", "pair", "done"]
OPT22 = ["`--form c2-s2 --gata4-exponent 2.0`", "`--gata4-exponent 2.0`", "`--captain-ruling R22`"]
OPT41 = ["`T3_c2_s2_g2.h5ad`", "`--form c2-s2 --gata4-exponent 2.0`"]
ROW_R22 = ("| R22 (the captain's order of a day; ruling, no run, not a try; by T3's second search line the builder is "
           "also to take `--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`; by T3's third "
           "search line also `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`) | T3 | x |\n")
ROW_R41 = ("| R41 (the captain's yes of a day, \"yes, T3_c2_s2_g2\": one file, `T3_c2_s2_g2.h5ad`; ruling, no run, "
           "not a try; the builder is to take `--form c2-s2 --gata4-exponent 2.0` only with `--captain-ruling R22`) "
           "| T3 | y |\n")

# a stand-in for uv: it logs each call and fakes pytest, the builder (its exit set per output folder) and vec check
FAKE_UV = r"""#!/bin/bash
F="$HOME/fake"
mkdir -p "$F"
printf '%s\n' "$*" >> "$F/calls.log"
ctl() { if [ -e "$F/$1" ]; then cat "$F/$1"; else echo 0; fi; }
[ "${1:-}" = run ] && shift
[ "${1:-}" = --frozen ] && shift
[ "${1:-}" = --no-sync ] && shift
if [ "${1:-}" = vec ]; then exit "$(ctl check_exit)"; fi
[ "${1:-}" = python ] && shift
while [ "${1:-}" = -I ] || [ "${1:-}" = -B ]; do shift; done
if [ "${1:-}" = -c ]; then exec "$FAKE_PY" -I -B "$@"; fi
if [ "${1:-}" = -m ] && [ "${2:-}" = pytest ]; then exit "$(ctl tests_exit)"; fi
shift
[ "${1:-}" = build ] || exit 64
shift 2
out=""
rest=""
while [ $# -gt 0 ]; do
  case "$1" in
    --out) out=$2; shift 2 ;;
    --wt|--run-dir|--board|--seed) shift 2 ;;
    *) rest="$rest $1"; shift ;;
  esac
done
if [ -e "$F/fail" ]; then read -r pat code < "$F/fail"; case "$out" in $pat) exit "$code" ;; esac; fi
extra=""
if [ -e "$F/differ" ]; then pat=$(cat "$F/differ"); case "$out" in $pat) extra=" $out" ;; esac; fi
printf 'fake build%s%s\n' "$rest" "$extra" > "$out"
clean=true
if [ -e "$F/unclean" ]; then pat=$(cat "$F/unclean"); case "$out" in $pat) clean=false ;; esac; fi
printf '{"commit": "%s", "clean": %s}\n' "$(git rev-parse HEAD)" "$clean" > "$out.stamp.json"
"""


def _sh() -> str:
    return LAUNCHER.read_text()


def _code() -> list[str]:
    return [line.split("#")[0] for line in _sh().splitlines()]


def _sums(name: str) -> dict[str, str]:
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


def _builder():
    """experiments/t3_decomp/build_two_knockout.py as loaded in this session (loaded here if no test module has)."""
    mod = sys.modules.get("build_two_knockout")
    if mod is None:
        spec = importlib.util.spec_from_file_location("build_two_knockout", REPO / "experiments" / "t3_decomp" /
                                                      "build_two_knockout.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["build_two_knockout"] = mod
        spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ what the launcher runs


def test_the_launcher_parses_and_deletes_fetches_and_uploads_nothing():
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], check=False).returncode == 0
    code = _code()
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)(\s|$)", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash|fetch|pull|push|clone)", line)
                   for line in code)
    assert not any(re.search(r"(^|[\s;&|(])(curl|wget|scp|rsync|ssh|sftp|ftp)(\s|$)", line) for line in code)
    assert not any(re.search(r"upload|submit|format_only|token", line, re.I) for line in code)
    assert re.findall(r"\$\{\d+:\?([^}]*)\}", _sh()) == ["usage: t3_c2g_build.sh COMMIT RUN"] * 2
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", _sh()):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


def test_the_recorded_hashes_are_the_committed_build_records():
    sh = _sh()
    dc7, two, g2 = _sums("builds_dc7b6ae.SHA256SUMS"), _sums("t3two_fb01487.SHA256SUMS"), _sums(
        "t3g2_a3c511c.SHA256SUMS")
    want = {"IDS1_SHA": dc7["t3_build/barred_s1.h5ad"], "IDS2_SHA": dc7["t3_build/barred_s2.h5ad"],
            "C1S2_SHA": two["c1-s2.h5ad"], "C2S2_SHA": two["c2-s2.h5ad"], "C1G2_SHA": g2["a/T3_c1_s2_g2.h5ad"]}
    assert dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", sh, flags=re.M)) == want
    for sha in want.values():
        assert sh.count(sha) == 1
    # the same five as row R36's launcher
    cd7 = (REPO / "experiments" / "t3_decomp" / "t3_cd7_build.sh").read_text()
    assert dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", cd7, flags=re.M)) == want


def test_the_commands_and_the_step_order():
    sh = _sh()
    assert f'C2G_ARGS="{C2G}"' in sh and f"FILE={FILE}" in sh
    assert f"DEC={DEC}" in sh and f"WT={WT}" in sh and "BOARD=T3:bcatenin" in sh
    assert "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only" in sh
    assert ('$PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board '
            '"$BOARD" \\\n    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?') in sh
    assert ("  if [ $rc = 2 ]; then return 12; fi\n  if [ $rc = 6 ]; then return 16; fi\n"
            "  if [ $rc = 5 ]; then return 17; fi\n  if [ $rc != 0 ]; then return 11; fi\n") in sh
    assert 'uv run --frozen --no-sync vec check "$d/$name" --board "$BOARD" > /dev/null 2>&1 || return 13' in sh
    for step, name, var, code, args in IDENTITY:
        assert f'STEP={step}\nbegin\nidentity {name} "${var}" {code} {args}\nstep_ok\n' in sh
    assert 'for STEP in a b; do begin; build "$STEP" "$FILE" $C2G_ARGS || finish $?; step_ok; done' in sh
    assert 'cmp -s "$O/a/$FILE" "$O/b/$FILE" ||' in sh
    assert '(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 33' in sh
    marks = ["STEP=gate", "STEP=tests", *[f"STEP={i[0]}\n" for i in IDENTITY], "for STEP in a b", "STEP=pair",
             "STEP=done"]
    order = [sh.index(m) for m in marks]
    assert order == sorted(order)
    assert '[ "$st" = "$COMMIT True" ] || return 15' in sh
    assert "O=runs/t3c2g_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t3c2g/$RUN" in sh
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert sh.index("PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=ids1")
    for line in sh.splitlines():
        if re.search(r"build_two_knockout\.py build|vec check \"", line):
            assert "> /dev/null 2>&1" in line or line.rstrip().endswith("\\"), line
    assert ".json" not in sh.replace("$1.stamp.json", "").replace("gata4_decisions.json", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")
    assert " g1 " not in sh and " g2 " not in sh  # row R36's gates are not this line's


def test_the_tests_step_lists_existing_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"(?:tests|experiments/t3_decomp)/test_[a-z0-9_]+\.py", step)
    assert {"tests/test_t3_c2g_build.py", "tests/test_t3_c2g_old_new.py", "tests/test_t3_c2g_launcher.py",
            "tests/test_t3_two_knockout_build.py", "tests/test_t3_gata4_weight.py",
            "experiments/t3_decomp/test_build_submission.py"} <= set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert '> "$S/tests.out" 2>&1 || finish 10' in step


def test_the_gate_lines_and_the_builder_s_strings():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert '[ -z "$(git status --porcelain 2> /dev/null)" ]' in gate
    assert "[ \"$(grep -c '^| R22 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "for opt in " + " ".join(f"'{o}'" for o in OPT22) + "; do" in gate
    assert "[ \"$(grep -c '^| R41 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "for opt in " + " ".join(f"'{o}'" for o in OPT41) + "; do" in gate
    assert [int(c) for c in re.findall(r"finish (\d+)", gate)] == [96, 98, 95, 95, 94, 94, 93]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")
    bk = _builder()
    assert bk.C2G_OPTION == OPT22[0] and bk.G2_OPTION == OPT22[1] and bk.C2G_FILE == FILE
    assert (bk.EXIT_GUARD_FAILED, bk.EXIT_FACTOR_MISMATCH, bk.EXIT_REFUSED) == (6, 5, 2)
    assert C2G.split() == ["--form", bk.C2G_FORM, "--gata4-exponent", "2.0", "--captain-ruling", bk.G2_ROW]


def test_the_builder_takes_the_command_the_launcher_runs(tmp_path):
    bk = _builder()
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R22}{ROW_R41}")
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", p)["option"] == OPT22[0]
    assert bk.check_gata4("c1-s2", 2.0, "R22", p)["option"] == OPT22[1]  # the last identity gate's line


# ------------------------------------------------------------------ the launcher, run for real


def _content(args: str) -> bytes:
    return f"fake build {args}\n".encode()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _with_stand_in_hashes(sh: str) -> str:
    out = sh
    for _, _, var, _, args in IDENTITY:
        out, n = re.subn(rf"^{var}=[0-9a-f]{{64}}", f"{var}={_sha(_content(args))}", out, flags=re.M)
        assert n == 1, var
    assert [a for a, b in zip(sh.splitlines(), out.splitlines(), strict=True) if a != b] == [
        ln for ln in sh.splitlines() if re.match(r"^[A-Z0-9]+_SHA=", ln)]
    return out


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env |= {"HOME": str(tmp_path / "home"), "GIT_CEILING_DIRECTORIES": str(tmp_path), "FAKE_PY": sys.executable,
            "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R22 + ROW_R41, sh: str | None = None, ctl: dict | None = None
          ) -> tuple[Path, dict, str]:
    env = _env(tmp_path)
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text(FAKE_UV)
    (fake / "uv").chmod(0o755)
    (tmp_path / "home" / "fake").mkdir(exist_ok=True)
    for name, value in (ctl or {}).items():
        (tmp_path / "home" / "fake" / name).write_text(f"{value}\n")
    repo = tmp_path / "repo"
    (repo / "experiments" / "t3_decomp").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / REL).write_text(sh if sh is not None else _sh())
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / ".gitignore").write_text("data/\nruns/\n")
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", REL, *argv], cwd=repo, env=env, capture_output=True, text=True, timeout=120,
                          check=False)


def _scratch(tmp_path: Path, run: str = "r1") -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t3c2g" / run


def _exits(s: Path) -> dict[str, str]:
    return {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")}


def _calls(tmp_path: Path) -> list[str]:
    log = tmp_path / "home" / "fake" / "calls.log"
    return log.read_text().splitlines() if log.exists() else []


def test_a_whole_run_with_the_stand_in_hashes(tmp_path):
    repo, env, head = _repo(tmp_path, sh=_with_stand_in_hashes(_sh()))
    res = _run(repo, env, head, "r1")
    assert res.returncode == 0, res.stderr
    s = _scratch(tmp_path)
    assert _exits(s) == {**dict.fromkeys(STEPS, "0"), "all": "0"}
    assert (s / "head").read_text().strip() == head
    o = repo / "runs" / f"t3c2g_{head[:7]}_r1"
    sums = (s / "SHA256SUMS").read_text()
    assert (o / "SHA256SUMS").read_text() == sums
    words = sums.split()
    want = {f"gate/{name}": _sha(_content(args)) for _, name, _, _, args in IDENTITY}
    want |= {f"a/{FILE}": _sha(_content(C2G)), f"b/{FILE}": _sha(_content(C2G))}
    assert dict(zip(words[1::2], words[0::2], strict=True)) == want
    assert (s / "stamps").read_text().splitlines() == [f"gate/{name} {head} True" for _, name, _, _, _ in IDENTITY] + [
        f"a/{FILE} {head} True", f"b/{FILE} {head} True"]
    run_dir = f"{tmp_path / 'home'}/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
    calls = _calls(tmp_path)
    assert calls[0].startswith("run --frozen --no-sync python -I -B -m pytest -q -p no:cacheprovider")
    fixed = f"--wt {WT} --run-dir {run_dir} --board T3:bcatenin --seed 1000001"
    tool = [c for c in calls if "build_two_knockout.py" in c]
    want_tool = [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} {fixed} "
                 f"{args} --out runs/t3c2g_{head[:7]}_r1/gate/{name}" for _, name, _, _, args in IDENTITY]
    want_tool += [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} "
                  f"{fixed} {C2G} --out runs/t3c2g_{head[:7]}_r1/{side}/{FILE}" for side in "ab"]
    assert tool == want_tool
    assert len([c for c in calls if c.startswith("run --frozen --no-sync vec check")]) == 7
    assert sorted(p.name for p in o.iterdir()) == ["SHA256SUMS", "a", "b", "gate"]


@pytest.mark.parametrize("ctl,step,code", [
    ({"tests_exit": 1}, "tests", 10),
    ({"fail": "*/gate/identity-s1.h5ad 1"}, "ids1", 11),
    ({"fail": "*/gate/identity-s1.h5ad 2"}, "ids1", 12),
    ({"fail": "*/gate/T3_c1_s2_g2.h5ad 5"}, "c1g2", 17),  # line 2's own factor check
    ({"check_exit": 1}, "ids1", 13),
    ({"unclean": "*/gate/c2-s2.h5ad"}, "c2s2", 15),
    ({"fail": "*/a/* 6"}, "a", 16),  # the guard stops the file: nothing written
    ({"fail": "*/a/* 5"}, "a", 17),  # the factor is not the guard's
    ({"fail": "*/b/* 2"}, "b", 12),
    ({"unclean": "*/a/*"}, "a", 15),
    ({"differ": "*/b/*"}, "pair", 31),
])
def test_each_stop_after_the_gate(tmp_path, ctl, step, code):
    repo, env, head = _repo(tmp_path, sh=_with_stand_in_hashes(_sh()), ctl=ctl)
    res = _run(repo, env, head, "r1")
    assert res.returncode == code, res.stderr
    s = _scratch(tmp_path)
    ex = _exits(s)
    assert ex.pop("all") == str(code) and ex.pop(step) == str(code)
    assert ex == dict.fromkeys(STEPS[:STEPS.index(step)], "0")  # every earlier step passed, no later one began
    assert not (s / "SHA256SUMS").exists()
    o = repo / "runs" / f"t3c2g_{head[:7]}_r1"
    if code in (16, 17) and step == "a":
        assert not (o / "a" / FILE).exists() and not (o / "b").exists()


def test_the_launcher_s_own_hashes_stop_the_stand_in_at_the_first_identity_gate(tmp_path):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, head, "r1")
    assert res.returncode == 20 and "gate/identity-s1.h5ad is not its recorded bytes" in res.stderr
    assert _exits(_scratch(tmp_path)) == {"gate": "0", "tests": "0", "ids1": "20", "all": "20"}


def _stopped_at_the_gate(tmp_path: Path, code: int, run: str = "r1") -> None:
    s = _scratch(tmp_path, run)
    assert _exits(s) == {"gate": str(code), "all": str(code)}
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "tmp"]
    assert not list((s / "tmp").iterdir()) and not _calls(tmp_path)


@pytest.mark.parametrize("tries,code,says", [
    (ROW_R41, 95, "exactly one row R22"),
    (ROW_R22 + ROW_R22 + ROW_R41, 95, "exactly one row R22"),
    (ROW_R22.replace("| R22 (", "| R22b (") + ROW_R41, 95, "exactly one row R22"),
    (ROW_R22.replace("`--form c2-s2 --gata4-exponent 2.0`", "--form c2-s2 --gata4-exponent 2.0") + ROW_R41, 95,
     "does not name `--form c2-s2 --gata4-exponent 2.0`"),  # line 3 not in backticks
    (ROW_R22.replace("`--gata4-exponent 2.0` on", "the weight on") + ROW_R41, 95,
     "does not name `--gata4-exponent 2.0`"),  # line 2 gone (the last identity gate needs it)
    (ROW_R22.replace("`--captain-ruling R22`", "R22") + ROW_R41, 95, "does not name `--captain-ruling R22`"),
    (ROW_R22.replace("; by T3's third search line also `--form c2-s2 --gata4-exponent 2.0`, only with "
                     "`--captain-ruling R22`) | T3 | x |", ") | T3 | `--form c2-s2 --gata4-exponent 2.0` |") + ROW_R41,
     95, "does not name `--form c2-s2 --gata4-exponent 2.0`"),  # named outside the first cell only
    (ROW_R22, 94, "exactly one row R41"),
    (ROW_R22 + ROW_R41 + ROW_R41, 94, "exactly one row R41"),
    (ROW_R22 + ROW_R41.replace("`T3_c2_s2_g2.h5ad`", "T3_c2_s2_g2.h5ad"), 94, "does not name `T3_c2_s2_g2.h5ad`"),
    (ROW_R22 + ROW_R41.replace("`--form c2-s2 --gata4-exponent 2.0`", "the option"), 94,
     "does not name `--form c2-s2 --gata4-exponent 2.0`"),
])
def test_each_ledger_refusal(tmp_path, tries, code, says):
    repo, env, head = _repo(tmp_path, tries=tries)
    res = _run(repo, env, head, "r1")
    assert res.returncode == code and says in res.stderr, res.stderr
    _stopped_at_the_gate(tmp_path, code)
    assert not (repo / "runs").exists()


def test_a_wrong_commit_a_prefix_and_an_unclean_tree_are_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, "0" * 40, "r1")
    assert res.returncode == 96 and "HEAD is not" in res.stderr
    _stopped_at_the_gate(tmp_path, 96)
    assert _run(repo, env, head[:12], "r2").returncode == 96
    _stopped_at_the_gate(tmp_path, 96, "r2")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R22}{ROW_R41}more\n")
    res = _run(repo, env, head, "r3")
    assert res.returncode == 98 and "not clean" in res.stderr
    _stopped_at_the_gate(tmp_path, 98, "r3")
    _git(repo, env, "checkout", "-q", "--", "experiments/harness/TRIES.md")
    (repo / "stray.py").write_text("x = 1\n")
    assert _run(repo, env, head, "r4").returncode == 98
    _stopped_at_the_gate(tmp_path, 98, "r4")
    assert not (repo / "runs").exists()


def test_a_run_name_is_never_reused_and_an_odd_one_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1").returncode == 10
    s = _scratch(tmp_path)
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    odd = _run(repo, env, head, "a/b")
    assert odd.returncode == 90 and not _scratch(tmp_path, "a").exists()
    assert _run(repo, env, head).returncode != 0 and _run(repo, env).returncode != 0


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t3c2g_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    res = _run(repo, env, head, "r1")
    assert res.returncode == 93 and "output folder" in res.stderr
    _stopped_at_the_gate(tmp_path, 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_that_passes_the_gate_writes_the_head_and_its_folder(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1").returncode == 10
    s = _scratch(tmp_path)
    assert _exits(s) == {"gate": "0", "tests": "10", "all": "10"} and (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t3c2g_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())
