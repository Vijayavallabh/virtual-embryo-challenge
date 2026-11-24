"""The test-board build launcher of row R36's T3 file cD7, experiments/t3_decomp/t3_cd7_build.sh <COMMIT> <RUN>.

It gates on a clean tree at COMMIT and on rows R36 and R22 of the ledger, runs the targeted tests, rebuilds five known
files (identity-s1, identity-s2, c1-s2, c2-s2 and the c1-s2 file of row R22's line) and stops unless each gives its
recorded bytes, runs G1, builds T3_c2_s2_cd7.h5ad twice into separate folders and compares the pair, runs G2, and
writes one SHA256SUMS. These tests hold its strings, hashes and step order to that, check the hashes against the
committed SHA256SUMS records, check that it deletes, fetches and uploads nothing and reads no output of a test build,
and run it for real in throwaway repositories under a temporary home, with a stand-in `uv` (it logs each call and fakes
the tests, the builder, its gates and `vec check`): every gate refusal, every stop, and a whole run through a copy of
the launcher whose five recorded hashes are replaced by the stand-in's. Nothing here reads the repository's own ledger
or builds anything from released data.
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "t3_decomp" / "t3_cd7_build.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
REL = "experiments/t3_decomp/t3_cd7_build.sh"
DEC = "experiments/t3_decomp/results/gata4/gata4_decisions.json"
WT = "data/raw/T2/E8.75.h5ad"
FILE = "T3_c2_s2_cd7.h5ad"
CD7 = "--form c2-s2 --composition-from allowed-s1 --captain-ruling R36"
# step, folder file, sha256 variable, its exit code, the form arguments
IDENTITY = [("ids1", "identity-s1.h5ad", "IDS1_SHA", 20, "--form identity-s1"),
            ("ids2", "identity-s2.h5ad", "IDS2_SHA", 21, "--form identity-s2"),
            ("c1s2", "c1-s2.h5ad", "C1S2_SHA", 22, "--form c1-s2"),
            ("c2s2", "c2-s2.h5ad", "C2S2_SHA", 23, "--form c2-s2"),
            ("c1g2", "T3_c1_s2_g2.h5ad", "C1G2_SHA", 24, "--form c1-s2 --gata4-exponent 2.0 --captain-ruling R22")]
STEPS = ["gate", "tests", *[i[0] for i in IDENTITY], "g1", "a", "b", "pair", "g2", "done"]
OPT36 = ["`--composition-from allowed-s1`", "`--form c2-s2`", "`--captain-ruling R36`"]
OPT22 = ["`--gata4-exponent 2.0`", "`--captain-ruling R22`"]
ROW_R36 = ("| R36 (the captain's named exception of a day; ruling, no run, not a try; the T3 two-knockout builder is "
           "to take `--composition-from allowed-s1` with `--form c2-s2` only with `--captain-ruling R36`) | T3 | x |\n")
ROW_R22 = ("| R22 (the captain's order of a day; ruling, no run; the T3 two-knockout builder is also to take "
           "`--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`) | T3 | y |\n")


def _sh() -> str:
    return LAUNCHER.read_text()


def _code() -> list[str]:
    return [line.split("#")[0] for line in _sh().splitlines()]


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


# ------------------------------------------------------------------ what the launcher runs


def test_the_launcher_parses_and_deletes_fetches_and_uploads_nothing():
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], check=False).returncode == 0
    code = _code()
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)(\s|$)", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash|fetch|pull|push|clone)", line)
                   for line in code)
    assert not any(re.search(r"(^|[\s;&|(])(curl|wget|scp|rsync|ssh|sftp|ftp)(\s|$)", line) for line in code)
    assert not any(re.search(r"upload|submit|format_only|token", line.replace("G2 failed: no upload; C2 stays.", ""),
                             re.I) for line in code)  # the fixed line of a G2 failure says no upload
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert usages == ["usage: t3_cd7_build.sh COMMIT RUN"] * 2
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", _sh()):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


def test_the_recorded_hashes_are_the_committed_build_records():
    sh = _sh()
    dc7, two, g2 = _sums("builds_dc7b6ae.SHA256SUMS"), _sums("t3two_fb01487.SHA256SUMS"), _sums(
        "t3g2_a3c511c.SHA256SUMS")
    want = {"IDS1_SHA": dc7["t3_build/barred_s1.h5ad"], "IDS2_SHA": dc7["t3_build/barred_s2.h5ad"],
            "C1S2_SHA": two["c1-s2.h5ad"], "C2S2_SHA": two["c2-s2.h5ad"], "C1G2_SHA": g2["a/T3_c1_s2_g2.h5ad"]}
    assert g2["a/T3_c1_s2_g2.h5ad"] == g2["b/T3_c1_s2_g2.h5ad"]
    assert g2["gate/c1-s2.h5ad"] == want["C1S2_SHA"] and g2["gate/c2-s2.h5ad"] == want["C2S2_SHA"]
    assert g2["gate/identity-s1.h5ad"] == want["IDS1_SHA"] and g2["gate/identity-s2.h5ad"] == want["IDS2_SHA"]
    got = dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", sh, flags=re.M))
    assert got == want
    assert want == {"IDS1_SHA": "3712c3f9988d0f8aa4198b22c37b866a84ef8757fd500a9163944f0663a520a2",
                    "IDS2_SHA": "1b69521e3c28bc84185410e838efc6fcbcff898434e73f9a7547208a4c82711f",
                    "C1S2_SHA": "8718330866429dc089eaffd2f96f66936fe54134eba2d50f8126266e6a3e54b8",
                    "C2S2_SHA": "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7",
                    "C1G2_SHA": "e01e24facb6841a18a503655bcc5e9816304a993fe9a0dcfeca458c42ead1e16"}
    for sha in want.values():
        assert sh.count(sha) == 1


def test_the_commands_and_the_step_order():
    sh = _sh()
    assert f'CD7_ARGS="{CD7}"' in sh and f"FILE={FILE}" in sh
    assert f"DEC={DEC}" in sh and f"WT={WT}" in sh and "BOARD=T3:bcatenin" in sh
    assert "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only" in sh
    assert ('$PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board '
            '"$BOARD" \\\n    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?') in sh
    assert 'uv run --frozen --no-sync vec check "$d/$name" --board "$BOARD" > /dev/null 2>&1 || return 13' in sh
    for step, name, var, code, args in IDENTITY:
        assert f'STEP={step}\nbegin\nidentity {name} "${var}" {code} {args}\nstep_ok\n' in sh
    assert ('$PY experiments/t3_decomp/build_two_knockout.py g1 "$DEC" --wt "$WT" --composition-from allowed-s1 \\\n'
            '  > "$S/g1.out" 2> /dev/null || rc=$?') in sh
    assert ('$PY experiments/t3_decomp/build_two_knockout.py g2 "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board '
            '"$BOARD" \\\n  --seed 1000001 --composition-from allowed-s1 > "$S/g2.out" 2> /dev/null || rc=$?') in sh
    assert 'for STEP in a b; do begin; build "$STEP" "$FILE" $CD7_ARGS || finish $?; step_ok; done' in sh
    assert 'cmp -s "$O/a/$FILE" "$O/b/$FILE" ||' in sh
    assert '(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 33' in sh
    assert 'cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 33' in sh
    marks = ["STEP=gate", "STEP=tests", *[f"STEP={i[0]}\n" for i in IDENTITY], "STEP=g1", "for STEP in a b",
             "STEP=pair", "STEP=g2", "STEP=done"]
    order = [sh.index(m) for m in marks]
    assert order == sorted(order)
    assert '[ "$st" = "$COMMIT True" ] || return 15' in sh  # a build counts only with a stamp of this commit, clean
    assert "O=runs/t3cd7_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t3cd7/$RUN" in sh
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert sh.index("PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("STEP=ids1")
    for line in sh.splitlines():  # every builder and check call: output to /dev/null, read by its exit code only
        if re.search(r"build_two_knockout\.py build|vec check \"", line):
            assert "> /dev/null 2>&1" in line or line.rstrip().endswith("\\"), line
    # no sidecar record is read: the only JSON files named are the stamp and the decision file the builder takes
    assert ".json" not in sh.replace("$1.stamp.json", "").replace("gata4_decisions.json", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")  # runs from the folder it is called in


def test_the_tests_step_lists_existing_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"(?:tests|experiments/t3_decomp)/test_[a-z0-9_]+\.py", step)
    assert set(listed) == {"tests/test_t3_cd7_build.py", "tests/test_t3_cd7_old_new.py",
                           "tests/test_t3_cd7_launcher.py",
                           "tests/test_t3_two_knockout_build.py", "tests/test_t3_gata4_weight.py",
                           "tests/test_t3_barred_exponent_build.py", "tests/test_t3_combination_build.py",
                           "tests/test_t3_outside_forms.py", "experiments/t3_decomp/test_build_submission.py"}
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert '> "$S/tests.out" 2>&1 || finish 10' in step


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert '[ -z "$(git status --porcelain 2> /dev/null)" ]' in gate
    assert "[ \"$(grep -c '^| R36 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "first=$(grep '^| R36 ' experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert "for opt in " + " ".join(f"'{o}'" for o in OPT36) + "; do" in gate
    assert "[ \"$(grep -c '^| R22 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "for opt in " + " ".join(f"'{o}'" for o in OPT22) + "; do" in gate
    assert [int(c) for c in re.findall(r"finish (\d+)", gate)] == [96, 98, 95, 95, 94, 94, 93]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")


def _builder():
    """experiments/t3_decomp/build_two_knockout.py as loaded in this session (loaded here if no test module has)."""
    import importlib.util

    mod = sys.modules.get("build_two_knockout")
    if mod is None:
        spec = importlib.util.spec_from_file_location("build_two_knockout", REPO / "experiments" / "t3_decomp" /
                                                      "build_two_knockout.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["build_two_knockout"] = mod
        spec.loader.exec_module(mod)
    return mod


def test_the_builder_takes_the_options_the_launcher_names(tmp_path):
    """The option strings the gate looks for are the builder's own (R36_NAMES, G2_OPTION)."""
    mod = _builder()
    assert list(mod.R36_NAMES) == OPT36 and mod.G2_OPTION == OPT22[0] and mod.R36_FILE == FILE
    assert CD7.split() == ["--form", mod.R36_FORM, "--composition-from", "allowed-s1", "--captain-ruling",
                           mod.R36_ROW]
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R36}{ROW_R22}")
    assert mod.check_composition_from("c2-s2", "allowed-s1", "R36", None, p)["row"] == "R36"


# ------------------------------------------------------------------ the launcher, run for real


FAKE_UV = r"""#!/bin/bash
# a stand-in for uv: it logs each call and fakes pytest, the builder, its two gates and vec check
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
case "${1:-}" in
  g1) for i in 1 2 3 4 5; do echo 0.3835; done; exit "$(ctl g1_exit)" ;;
  g2) for b in 0.02 0.03 0.04 0.025 0.035 0.045; do echo "$b"; done; exit "$(ctl g2_exit)" ;;
  build) shift 2 ;;
  *) exit 64 ;;
esac
out=""
rest=""
while [ $# -gt 0 ]; do
  case "$1" in
    --out) out=$2; shift 2 ;;
    --wt|--run-dir|--board|--seed) shift 2 ;;
    *) rest="$rest $1"; shift ;;
  esac
done
code=$(ctl build_exit)
[ "$code" = 0 ] || exit "$code"
extra=""
if [ -e "$F/differ" ]; then pat=$(cat "$F/differ"); case "$out" in $pat) extra=" $out" ;; esac; fi
printf 'fake build%s%s\n' "$rest" "$extra" > "$out"
clean=true
if [ -e "$F/unclean" ]; then pat=$(cat "$F/unclean"); case "$out" in $pat) clean=false ;; esac; fi
printf '{"commit": "%s", "clean": %s}\n' "$(git rev-parse HEAD)" "$clean" > "$out.stamp.json"
"""


def _content(args: str) -> bytes:
    """What the stand-in builder writes for these form arguments."""
    return f"fake build {args}\n".encode()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _with_stand_in_hashes(sh: str) -> str:
    """The launcher with its five recorded hashes replaced by those of the stand-in's files, nothing else changed."""
    out = sh
    for _, _, var, _, args in IDENTITY:
        out, n = re.subn(rf"^{var}=[0-9a-f]{{64}}", f"{var}={_sha(_content(args))}", out, flags=re.M)
        assert n == 1, var
    assert [a for a, b in zip(sh.splitlines(), out.splitlines(), strict=True) if a != b] == [
        ln for ln in sh.splitlines() if re.match(r"^[A-Z0-9]+_SHA=", ln)]
    return out


def _env(tmp_path: Path) -> dict:
    """The caller's environment with a temporary home, git kept inside tmp_path, and the python running these tests
    for the stand-in's stamp reads (the launcher puts the home's ~/.local/bin, the stand-in uv, first on PATH)."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env |= {"HOME": str(tmp_path / "home"), "GIT_CEILING_DIRECTORIES": str(tmp_path), "FAKE_PY": sys.executable,
            "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R36 + ROW_R22, sh: str | None = None, ctl: dict | None = None
          ) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (or `sh`) and the ledger its gate reads, and a home whose
    ~/.local/bin holds the stand-in uv (the launcher puts that folder first on PATH), with its control files."""
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
    return tmp_path / "home" / "v" / "scratch" / "t3cd7" / run


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
    assert (s / "g1.out").read_text() == "0.3835\n" * 5
    assert (s / "g2.out").read_text() == "0.02\n0.03\n0.04\n0.025\n0.035\n0.045\n"
    assert not (s / "g2.verdict").exists() and not (s / "g2.rc").exists()
    assert "G2 failed" not in res.stderr
    o = repo / "runs" / f"t3cd7_{head[:7]}_r1"
    sums = (s / "SHA256SUMS").read_text()
    assert (o / "SHA256SUMS").read_text() == sums
    words = sums.split()
    got = dict(zip(words[1::2], words[0::2], strict=True))
    want = {f"gate/{name}": _sha(_content(args)) for _, name, _, _, args in IDENTITY}
    want |= {f"a/{FILE}": _sha(_content(CD7)), f"b/{FILE}": _sha(_content(CD7))}
    assert got == want
    stamps = (s / "stamps").read_text().splitlines()
    assert stamps == [f"gate/{name} {head} True" for _, name, _, _, _ in IDENTITY] + [
        f"a/{FILE} {head} True", f"b/{FILE} {head} True"]
    # the calls, in order: the tests, five identity builds, G1, the pair, G2; each build checked, then its stamp read
    run_dir = f"{tmp_path / 'home'}/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
    calls = _calls(tmp_path)
    assert calls[0].startswith("run --frozen --no-sync python -I -B -m pytest -q -p no:cacheprovider")
    tool = [c for c in calls if "build_two_knockout.py" in c]
    fixed = f"--wt {WT} --run-dir {run_dir} --board T3:bcatenin --seed 1000001"
    want_tool = [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} {fixed} "
                 f"{args} --out runs/t3cd7_{head[:7]}_r1/gate/{name}" for _, name, _, _, args in IDENTITY]
    want_tool.append(f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py g1 {DEC} --wt "
                     f"{WT} --composition-from allowed-s1")
    want_tool += [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} "
                  f"{fixed} {CD7} --out runs/t3cd7_{head[:7]}_r1/{side}/{FILE}" for side in "ab"]
    want_tool.append(f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py g2 {DEC} --wt "
                     f"{WT} --run-dir {run_dir} --board T3:bcatenin --seed 1000001 --composition-from allowed-s1")
    assert tool == want_tool
    checks = [c for c in calls if c.startswith("run --frozen --no-sync vec check")]
    assert len(checks) == 7 and all(c.endswith("--board T3:bcatenin") for c in checks)
    assert sorted(p.name for p in o.iterdir()) == ["SHA256SUMS", "a", "b", "gate"]


@pytest.mark.parametrize("ctl,step,code", [
    ({"tests_exit": 1}, "tests", 10),
    ({"build_exit": 1}, "ids1", 11),
    ({"build_exit": 2}, "ids1", 12),
    ({"check_exit": 1}, "ids1", 13),
    ({"unclean": "*/gate/c2-s2.h5ad"}, "c2s2", 15),
    ({"g1_exit": 1}, "g1", 30),
    ({"unclean": "*/a/*"}, "a", 15),
    ({"differ": "*/b/*"}, "pair", 31),
    ({"g2_exit": 1}, "g2", 32),
])
def test_each_stop_after_the_gate(tmp_path, ctl, step, code):
    repo, env, head = _repo(tmp_path, sh=_with_stand_in_hashes(_sh()), ctl=ctl)
    res = _run(repo, env, head, "r1")
    assert res.returncode == code, res.stderr
    s = _scratch(tmp_path)
    ex = _exits(s)
    assert ex.pop("all") == str(code) and ex.pop(step) == str(code)
    done = STEPS[:STEPS.index(step)]
    assert ex == dict.fromkeys(done, "0")  # every earlier step passed, no later one began
    assert not (s / "SHA256SUMS").exists()
    if STEPS.index(step) <= STEPS.index("g1"):
        assert not any(f"/a/{FILE}" in c for c in _calls(tmp_path))  # the file is never built before G1 passes


@pytest.mark.parametrize("g2_exit", [1, 2, 3])
def test_a_g2_failure_writes_the_fixed_line_after_the_two_builds(tmp_path, g2_exit):
    """G2 not passing (exit 1), refused (2) or not computed (3): the fixed line, no SHA256SUMS, exit 32, and the two
    builds of the file already made and compared."""
    repo, env, head = _repo(tmp_path, sh=_with_stand_in_hashes(_sh()), ctl={"g2_exit": g2_exit})
    res = _run(repo, env, head, "r1")
    assert res.returncode == 32
    s = _scratch(tmp_path)
    assert (s / "g2.verdict").read_text() == "G2 failed: no upload; C2 stays.\n"
    assert res.stderr.splitlines()[-1] == "G2 failed: no upload; C2 stays."
    assert (s / "g2.rc").read_text().strip() == str(g2_exit)
    assert _exits(s) == {**dict.fromkeys(STEPS[:STEPS.index("g2")], "0"), "g2": "32", "all": "32"}
    o = repo / "runs" / f"t3cd7_{head[:7]}_r1"
    assert (o / "a" / FILE).read_bytes() == (o / "b" / FILE).read_bytes() == _content(CD7)
    assert not (s / "SHA256SUMS").exists() and not (o / "SHA256SUMS").exists()


def test_the_fixed_line_is_the_builder_s():
    mod = _builder()
    sh = _sh()
    assert mod.G2_FAILED_LINE == "G2 failed: no upload; C2 stays."
    assert f'echo "{mod.G2_FAILED_LINE}" > "$S/g2.verdict"' in sh and f'echo "{mod.G2_FAILED_LINE}" >&2' in sh
    g2 = sh[sh.index("STEP=g2"):sh.index("STEP=done")]
    assert g2.index('echo "$rc" > "$S/g2.rc"') < g2.index("finish 32") and g2.count("finish 32") == 1


def test_the_launcher_s_own_hashes_stop_the_stand_in_at_the_first_identity_gate(tmp_path):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, head, "r1")
    assert res.returncode == 20 and "gate/identity-s1.h5ad is not its recorded bytes" in res.stderr
    ex = _exits(_scratch(tmp_path))
    assert ex == {"gate": "0", "tests": "0", "ids1": "20", "all": "20"}
    assert len([c for c in _calls(tmp_path) if "build_two_knockout.py" in c]) == 1


def _stopped_at_the_gate(tmp_path: Path, code: int, run: str = "r1") -> None:
    s = _scratch(tmp_path, run)
    assert _exits(s) == {"gate": str(code), "all": str(code)}
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "tmp"]  # no head file
    assert not list((s / "tmp").iterdir()) and not _calls(tmp_path)  # nothing ran


@pytest.mark.parametrize("tries,code,says", [
    (ROW_R22, 95, "exactly one row R36"),
    (ROW_R36 + ROW_R36 + ROW_R22, 95, "exactly one row R36"),
    (ROW_R36.replace("| R36 (", "| R36b (") + ROW_R22, 95, "exactly one row R36"),
    (ROW_R36.replace("`--composition-from allowed-s1`", "--composition-from allowed-s1") + ROW_R22, 95,
     "does not name `--composition-from allowed-s1`"),
    (ROW_R36.replace("`--form c2-s2`", "`--form c1-s2`") + ROW_R22, 95, "does not name `--form c2-s2`"),
    (ROW_R36.replace("`--captain-ruling R36`", "R36") + ROW_R22, 95, "does not name `--captain-ruling R36`"),
    (ROW_R36.replace(" with `--form c2-s2` only with `--captain-ruling R36`) | T3 | x |",
                     ") | T3 | `--form c2-s2` `--captain-ruling R36` |") + ROW_R22, 95,
     "does not name `--form c2-s2`"),  # named outside the first cell only
    (ROW_R36, 94, "exactly one row R22"),
    (ROW_R36 + ROW_R22 + ROW_R22, 94, "exactly one row R22"),
    (ROW_R36 + ROW_R22.replace("`--gata4-exponent 2.0`", "`--gata4-exponent 3.0`"), 94,
     "does not name `--gata4-exponent 2.0`"),
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
    res = _run(repo, env, head[:12], "r2")
    assert res.returncode == 96
    _stopped_at_the_gate(tmp_path, 96, "r2")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R36}{ROW_R22}more\n")
    res = _run(repo, env, head, "r3")
    assert res.returncode == 98 and "not clean" in res.stderr
    _stopped_at_the_gate(tmp_path, 98, "r3")
    _git(repo, env, "checkout", "-q", "--", "experiments/harness/TRIES.md")
    (repo / "stray.py").write_text("x = 1\n")  # an untracked file
    res = _run(repo, env, head, "r4")
    assert res.returncode == 98
    _stopped_at_the_gate(tmp_path, 98, "r4")
    assert not (repo / "runs").exists()


def test_a_run_name_is_never_reused_and_an_odd_one_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    res = _run(repo, env, head, "r1")
    assert res.returncode == 10
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
    o = repo / "runs" / f"t3cd7_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    res = _run(repo, env, head, "r1")
    assert res.returncode == 93 and "output folder" in res.stderr
    _stopped_at_the_gate(tmp_path, 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_that_passes_the_gate_writes_the_head_and_its_folder(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    res = _run(repo, env, head, "r1")
    assert res.returncode == 10
    s = _scratch(tmp_path)
    assert _exits(s) == {"gate": "0", "tests": "10", "all": "10"} and (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t3cd7_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built
