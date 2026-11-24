"""The test-board build launcher of the options on C2's command:
experiments/t3_decomp/t3_gonly_build.sh <COMMIT> <RUN> [LIFT_ROW].

It gates on a clean tree at COMMIT and on row R22 of the ledger (each file's option string, the two old lines' strings
its last identity gates need and `--captain-ruling R22`, each as one whole backtick-delimited string), and, given
LIFT_ROW, on that row (a ruling saying "the 5% check gates nothing for it" and naming at least one file's string; the
files it names are built with `--den-ruling LIFT_ROW`). It runs the targeted tests, rebuilds six known files
(identity-s1, identity-s2, c1-s2, c2-s2, line 2's T3_c1_s2_g2.h5ad and line 3's T3_c2_s2_g2.h5ad) and stops unless
each gives its recorded bytes, then builds F1 T3_c2_s2_gonly.h5ad, F3 T3_c2_s2_gonly_s067.h5ad and F2
T3_c2_s2_s067.h5ad each twice into separate folders (the builder's range guard, 5% den check and factor check run
inside each build), compares each pair and writes one SHA256SUMS. These tests hold its strings, hashes and step order
to that, check the hashes against the committed SHA256SUMS records, check that it deletes, fetches and uploads nothing
and reads no output of a test build, and run it for real in throwaway repositories under a temporary home with a
stand-in `uv` (it logs each call and fakes the tests, the builder and `vec check`): every gate refusal, every stop, and
whole runs through a copy of the launcher whose six recorded hashes are replaced by the stand-in's. Nothing here reads
the repository's own ledger or builds anything from released data.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "t3_decomp" / "t3_gonly_build.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
REL = "experiments/t3_decomp/t3_gonly_build.sh"
DEC = "experiments/t3_decomp/results/gata4/gata4_decisions.json"
WT = "data/raw/T2/E8.75.h5ad"
USAGE = "usage: t3_gonly_build.sh COMMIT RUN [LIFT_ROW]"
FILES = {"F1": "T3_c2_s2_gonly.h5ad", "F2": "T3_c2_s2_s067.h5ad", "F3": "T3_c2_s2_gonly_s067.h5ad"}
OPT = {"F1": "--form c2-s2 --gata4-exponent inf", "F2": "--form c2-s2 --size-scale 0.67",
       "F3": "--form c2-s2 --gata4-exponent inf --size-scale 0.67"}
ARGS = {f: f"{OPT[f]} --captain-ruling R22" for f in OPT}
IDENTITY = [("ids1", "identity-s1.h5ad", "IDS1_SHA", 20, "--form identity-s1"),
            ("ids2", "identity-s2.h5ad", "IDS2_SHA", 21, "--form identity-s2"),
            ("c1s2", "c1-s2.h5ad", "C1S2_SHA", 22, "--form c1-s2"),
            ("c2s2", "c2-s2.h5ad", "C2S2_SHA", 23, "--form c2-s2"),
            ("c1g2", "T3_c1_s2_g2.h5ad", "C1G2_SHA", 24, "--form c1-s2 --gata4-exponent 2.0 --captain-ruling R22"),
            ("c2g2", "T3_c2_s2_g2.h5ad", "C2G2_SHA", 25, "--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22")]
BUILDS = [("f1a", "a", "F1"), ("f1b", "b", "F1"), ("f3a", "a", "F3"), ("f3b", "b", "F3"), ("f2a", "a", "F2"),
          ("f2b", "b", "F2")]
STEPS = ["gate", "tests", *[i[0] for i in IDENTITY], *[b[0] for b in BUILDS], "pair", "done"]
GATE22 = [OPT["F1"], OPT["F3"], OPT["F2"], "--form c2-s2 --gata4-exponent 2.0", "--gata4-exponent 2.0",
          "--captain-ruling R22"]
ROW_R22 = ("| R22 (the captain's order of a day; ruling, no run, not a try; by T3's second search line the builder is "
           "also to take `--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`; by T3's third "
           "search line also `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`; by the options on "
           f"C2's command also `{OPT['F1']}`, `{OPT['F3']}` and `{OPT['F2']}`, each only with `--captain-ruling R22`) "
           "| T3 | x |\n")
ROW_R41 = "| R41 (the captain's yes of a day; ruling, no run, not a try) | T3 | y |\n"


def lift_row(*names: str, row: str = "R44", words: str = "the 5% check gates nothing for it",
             ruling: str = "ruling, no run, not a try") -> str:
    files = " and ".join(f"`{OPT[n]}`" for n in names)
    return f"| {row} (the captain's yes of a day; {ruling}; for {files}: {words}) | T3 | y |\n"


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
    assert re.findall(r"\$\{\d+:\?([^}]*)\}", _sh()) == [USAGE] * 2
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", _sh()):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg
    assert _sh().startswith("#!/usr/bin/env bash\n") and "\nset -u\n" in _sh() and "LIFT=${3:-}" in _sh()


def test_the_recorded_hashes_are_the_committed_build_records():
    sh = _sh()
    dc7, two, g2 = _sums("builds_dc7b6ae.SHA256SUMS"), _sums("t3two_fb01487.SHA256SUMS"), _sums(
        "t3g2_a3c511c.SHA256SUMS")
    c2g = _sums("t3c2g_fc57ca9_g1.SHA256SUMS")
    assert c2g["a/T3_c2_s2_g2.h5ad"] == c2g["b/T3_c2_s2_g2.h5ad"]
    want = {"IDS1_SHA": dc7["t3_build/barred_s1.h5ad"], "IDS2_SHA": dc7["t3_build/barred_s2.h5ad"],
            "C1S2_SHA": two["c1-s2.h5ad"], "C2S2_SHA": two["c2-s2.h5ad"], "C1G2_SHA": g2["a/T3_c1_s2_g2.h5ad"],
            "C2G2_SHA": c2g["a/T3_c2_s2_g2.h5ad"]}
    assert dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", sh, flags=re.M)) == want
    assert want["C2G2_SHA"] == "014d2bf9d8ba7fdc56dcfba2a3b8aa77f0ea20bb806ebc2c1193e483d62712ac"
    for sha in want.values():
        assert sh.count(sha) == 1
    # the first five as t3_c2g_build.sh's, whose own identity gates t3c2g_fc57ca9_g1 records
    c2g_sh = (REPO / "experiments" / "t3_decomp" / "t3_c2g_build.sh").read_text()
    five = dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", c2g_sh, flags=re.M))
    assert five == {k: v for k, v in want.items() if k != "C2G2_SHA"}
    assert {c2g[f"gate/{name}"] for _, name, _, _, _ in IDENTITY[:5]} == set(five.values())


def test_the_commands_and_the_step_order():
    sh = _sh()
    for name in FILES:
        assert f'{name}_OPT="{OPT[name]}"' in sh and f'{name}_ARGS="${name}_OPT --captain-ruling R22"' in sh
        assert f"{name}={FILES[name]}" in sh
    assert f"DEC={DEC}" in sh and f"WT={WT}" in sh and "BOARD=T3:bcatenin" in sh
    assert "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only" in sh
    assert ('$PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board '
            '"$BOARD" \\\n    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?') in sh
    assert ("  if [ $rc = 2 ]; then return 12; fi\n  if [ $rc = 6 ]; then return 16; fi\n"
            "  if [ $rc = 5 ]; then return 17; fi\n  if [ $rc != 0 ]; then return 11; fi\n") in sh
    assert 'uv run --frozen --no-sync vec check "$d/$name" --board "$BOARD" > /dev/null 2>&1 || return 13' in sh
    for step, name, var, code, args in IDENTITY:
        assert f'STEP={step}\nbegin\nidentity {name} "${var}" {code} {args}\nstep_ok\n' in sh
    for n in ("1", "3", "2"):
        assert (f'for STEP in f{n}a f{n}b; do begin; build "${{STEP:2}}" "$F{n}" $F{n}_ARGS || finish $?; step_ok; '
                'done') in sh
    for n, code in (("1", 31), ("3", 33), ("2", 32)):
        assert (f'cmp -s "$O/a/$F{n}" "$O/b/$F{n}" || {{ echo "stopped: the two builds of $F{n} differ" >&2; '
                f'finish {code}; }}') in sh
    assert '(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 39' in sh
    marks = ["STEP=gate", "STEP=tests", *[f"STEP={i[0]}\n" for i in IDENTITY], "for STEP in f1a f1b",
             "for STEP in f3a f3b", "for STEP in f2a f2b", "STEP=pair", "STEP=done"]
    order = [sh.index(m) for m in marks]
    assert order == sorted(order)
    assert '[ "$st" = "$COMMIT True" ] || return 15' in sh
    assert "O=runs/t3gonly_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t3gonly/$RUN" in sh
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
    assert {"tests/test_t3_gonly_build.py", "tests/test_t3_gonly_old_new.py", "tests/test_t3_gonly_launcher.py",
            "tests/test_t3_c2g_build.py", "tests/test_t3_c2g_old_new.py", "tests/test_t3_c2g_launcher.py",
            "tests/test_t3_cd7_old_new.py", "tests/test_t3_two_knockout_build.py", "tests/test_t3_gata4_weight.py",
            "experiments/t3_decomp/test_build_submission.py"} <= set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert '> "$S/tests.out" 2>&1 || finish 10' in step


def test_the_gate_lines_and_the_builder_s_strings():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert '[ -z "$(git status --porcelain 2> /dev/null)" ]' in gate
    assert "[ \"$(grep -c '^| R22 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert ("for opt in \"$F1_OPT\" \"$F3_OPT\" \"$F2_OPT\" '--form c2-s2 --gata4-exponent 2.0' '--gata4-exponent "
            "2.0' \\\n    '--captain-ruling R22'; do\n  if ! names \"$first\" \"$opt\"; then") in gate
    assert "printf '%s\\n' \"$1\" | grep -o '`[^`]*`' | grep -qxF -- \"\\`$2\\`\"" in sh  # one whole string
    assert "for words in 'ruling, no run' 'the 5% check gates nothing for it'; do" in gate
    for n in ("1", "3", "2"):
        assert (f'if names "$lift" "$F{n}_OPT"; then F{n}_ARGS="$F{n}_ARGS --den-ruling $LIFT"; '
                f'LIFTED="$LIFTED F{n}"; fi') in gate
    assert [int(c) for c in re.findall(r"finish (\d+)", gate)] == [96, 98, 95, 95, 94, 94, 94, 93]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")
    bk = _builder()
    assert OPT == bk.GONLY_OPTIONS and FILES == bk.GONLY_FILES
    assert GATE22[3:] == [bk.C2G_OPTION.strip("`"), bk.G2_OPTION.strip("`"), f"--captain-ruling {bk.G2_ROW}"]
    assert bk.GONLY_LIFT_WORDS == "the 5% check gates nothing for it"
    assert (bk.EXIT_GUARD_FAILED, bk.EXIT_FACTOR_MISMATCH, bk.EXIT_REFUSED) == (6, 5, 2)


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path):
    bk = _builder()
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R22}{ROW_R41}{lift_row('F1', 'F2', 'F3')}")
    for name, args in ARGS.items():
        a = args.split()
        g = math.inf if "--gata4-exponent" in a else None
        s = float(a[a.index("--size-scale") + 1]) if "--size-scale" in a else None
        assert bk.is_gonly_command(a[1], g, s)
        rec = bk.check_c2_gonly(a[1], g, s, a[a.index("--captain-ruling") + 1], p)
        assert rec["name"] == name and rec["file"] == FILES[name] and rec["den_ruling"] is None
        lifted = bk.check_c2_gonly(a[1], g, s, "R22", p, den_ruling="R44")
        assert lifted["den_ruling"]["row"] == "R44" and lifted["den_ruling"]["option"] == f"`{OPT[name]}`"
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", p)["option"] == f"`{GATE22[3]}`"  # the last identity gate's line
    assert bk.check_gata4("c1-s2", 2.0, "R22", p)["option"] == f"`{GATE22[4]}`"  # the one before


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
    return tmp_path / "home" / "v" / "scratch" / "t3gonly" / run


def _exits(s: Path) -> dict[str, str]:
    return {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")}


def _calls(tmp_path: Path) -> list[str]:
    log = tmp_path / "home" / "fake" / "calls.log"
    return log.read_text().splitlines() if log.exists() else []


@pytest.mark.parametrize("lift,lifted", [(None, set()), ("F1 F3", {"F1", "F3"}), ("F2", {"F2"}),
                                         ("F1 F2 F3", {"F1", "F2", "F3"})])
def test_a_whole_run_with_the_stand_in_hashes(tmp_path, lift, lifted):
    tries = ROW_R22 + ROW_R41 + (lift_row(*lift.split()) if lift else "")
    repo, env, head = _repo(tmp_path, tries=tries, sh=_with_stand_in_hashes(_sh()))
    res = _run(repo, env, head, "r1", *(["R44"] if lift else []))
    assert res.returncode == 0, res.stderr
    s = _scratch(tmp_path)
    assert _exits(s) == {**dict.fromkeys(STEPS, "0"), "all": "0"}
    assert (s / "head").read_text().strip() == head
    assert (s / "lift").read_text().strip() == ("R44 " + " ".join(f for f in ("F1", "F3", "F2") if f in lifted)
                                                if lift else "none")
    args = {f: ARGS[f] + (" --den-ruling R44" if f in lifted else "") for f in FILES}
    o = repo / "runs" / f"t3gonly_{head[:7]}_r1"
    sums = (s / "SHA256SUMS").read_text()
    assert (o / "SHA256SUMS").read_text() == sums
    words = sums.split()
    want = {f"gate/{name}": _sha(_content(a)) for _, name, _, _, a in IDENTITY}
    want |= {f"{side}/{FILES[f]}": _sha(_content(args[f])) for _, side, f in BUILDS}
    assert dict(zip(words[1::2], words[0::2], strict=True)) == want and len(words) == 2 * 12
    assert len({want[f"a/{FILES[f]}"] for f in FILES}) == 3  # three different files, each pair one sha256
    assert (s / "stamps").read_text().splitlines() == [f"gate/{name} {head} True" for _, name, _, _, _ in IDENTITY] + [
        f"{side}/{FILES[f]} {head} True" for _, side, f in BUILDS]
    run_dir = f"{tmp_path / 'home'}/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
    calls = _calls(tmp_path)
    assert calls[0].startswith("run --frozen --no-sync python -I -B -m pytest -q -p no:cacheprovider")
    fixed = f"--wt {WT} --run-dir {run_dir} --board T3:bcatenin --seed 1000001"
    tool = [c for c in calls if "build_two_knockout.py" in c]
    want_tool = [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} {fixed} "
                 f"{a} --out runs/t3gonly_{head[:7]}_r1/gate/{name}" for _, name, _, _, a in IDENTITY]
    want_tool += [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} "
                  f"{fixed} {args[f]} --out runs/t3gonly_{head[:7]}_r1/{side}/{FILES[f]}" for _, side, f in BUILDS]
    assert tool == want_tool
    assert not any("--den-ruling" in c for c in tool[:6])  # the identity gates never take the lifting row
    assert len([c for c in calls if c.startswith("run --frozen --no-sync vec check")]) == 12
    assert sorted(p.name for p in o.iterdir()) == ["SHA256SUMS", "a", "b", "gate"]
    assert sorted(p.name for p in (o / "a").iterdir() if p.suffix == ".h5ad") == sorted(FILES.values())


@pytest.mark.parametrize("ctl,step,code", [
    ({"tests_exit": 1}, "tests", 10),
    ({"fail": "*/gate/identity-s1.h5ad 1"}, "ids1", 11),
    ({"fail": "*/gate/identity-s1.h5ad 2"}, "ids1", 12),
    ({"fail": "*/gate/T3_c1_s2_g2.h5ad 5"}, "c1g2", 17),  # line 2's own factor check
    ({"fail": "*/gate/T3_c2_s2_g2.h5ad 6"}, "c2g2", 16),  # line 3's own 5% guard
    ({"fail": "*/gate/T3_c2_s2_g2.h5ad 5"}, "c2g2", 17),
    ({"check_exit": 1}, "ids1", 13),
    ({"unclean": "*/gate/c2-s2.h5ad"}, "c2s2", 15),
    ({"unclean": "*/gate/T3_c2_s2_g2.h5ad"}, "c2g2", 15),
    ({"fail": "*/a/T3_c2_s2_gonly.h5ad 6"}, "f1a", 16),  # a check stops F1: a stop, nothing written
    ({"fail": "*/a/T3_c2_s2_gonly.h5ad 5"}, "f1a", 17),  # the factor is not the checks'
    ({"fail": "*/b/T3_c2_s2_gonly.h5ad 2"}, "f1b", 12),
    ({"fail": "*/a/T3_c2_s2_gonly_s067.h5ad 6"}, "f3a", 16),
    ({"unclean": "*/b/T3_c2_s2_gonly_s067.h5ad"}, "f3b", 15),
    ({"fail": "*/a/T3_c2_s2_s067.h5ad 6"}, "f2a", 16),
    ({"fail": "*/b/T3_c2_s2_s067.h5ad 1"}, "f2b", 11),
    ({"differ": "*/b/T3_c2_s2_gonly.h5ad"}, "pair", 31),
    ({"differ": "*/b/T3_c2_s2_gonly_s067.h5ad"}, "pair", 33),
    ({"differ": "*/b/T3_c2_s2_s067.h5ad"}, "pair", 32),
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
    o = repo / "runs" / f"t3gonly_{head[:7]}_r1"
    if code in (16, 17) and step.startswith("f"):
        side, f = {b[0]: (b[1], b[2]) for b in BUILDS}[step]
        assert not (o / side / FILES[f]).exists()


def test_the_launcher_s_own_hashes_stop_the_stand_in_at_the_first_identity_gate(tmp_path):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, head, "r1")
    assert res.returncode == 20 and "gate/identity-s1.h5ad is not its recorded bytes" in res.stderr
    assert _exits(_scratch(tmp_path)) == {"gate": "0", "tests": "0", "ids1": "20", "all": "20"}


def test_the_last_identity_gate_stops_on_bytes_other_than_line_3_s_record(tmp_path):
    """The stand-in's hashes for the first five gates only: the sixth (T3_c2_s2_g2.h5ad) keeps its recorded sha256."""
    sh = _with_stand_in_hashes(_sh())
    sh = re.sub(r"^C2G2_SHA=[0-9a-f]{64}", "C2G2_SHA=014d2bf9d8ba7fdc56dcfba2a3b8aa77f0ea20bb806ebc2c1193e483d62712ac",
                sh, flags=re.M)
    repo, env, head = _repo(tmp_path, sh=sh)
    res = _run(repo, env, head, "r1")
    assert res.returncode == 25 and "gate/T3_c2_s2_g2.h5ad is not its recorded bytes" in res.stderr
    ex = _exits(_scratch(tmp_path))
    assert ex.pop("c2g2") == "25" and ex.pop("all") == "25" and ex == dict.fromkeys(STEPS[:STEPS.index("c2g2")], "0")


def _stopped_at_the_gate(tmp_path: Path, code: int, run: str = "r1") -> None:
    s = _scratch(tmp_path, run)
    assert _exits(s) == {"gate": str(code), "all": str(code)}
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "tmp"]
    assert not list((s / "tmp").iterdir()) and not _calls(tmp_path)


R22_REFUSALS = [
    (ROW_R41, "exactly one row R22"),
    (ROW_R22 + ROW_R22 + ROW_R41, "exactly one row R22"),
    (ROW_R22.replace("| R22 (", "| R22b (") + ROW_R41, "exactly one row R22"),
    *[(ROW_R22.replace(f"`{o}`", o) + ROW_R41, f"does not name `{o}` as one whole string") for o in GATE22],
    (ROW_R22.replace(f"`{OPT['F1']}`", "`--gata4-exponent inf`") + ROW_R41,
     f"does not name `{OPT['F1']}` as one whole string"),  # without its form
    (ROW_R22.replace(f"`{OPT['F1']}`", f"`a`{OPT['F1']}`b`") + ROW_R41,
     f"does not name `{OPT['F1']}` as one whole string"),  # in the text, not as one whole string
    (ROW_R22.replace(f"`{OPT['F2']}`", "`--size-scale 0.67`") + ROW_R41,
     f"does not name `{OPT['F2']}` as one whole string"),
    (ROW_R22.replace(f"; by the options on C2's command also `{OPT['F1']}`, `{OPT['F3']}` and `{OPT['F2']}`, each "
                     "only with `--captain-ruling R22`) | T3 | x |",
                     f") | T3 | `{OPT['F1']}` `{OPT['F3']}` `{OPT['F2']}` |") + ROW_R41,
     f"does not name `{OPT['F1']}` as one whole string"),  # named outside the first cell only
]


@pytest.mark.parametrize("tries,says", R22_REFUSALS, ids=[f"r22-{i}" for i in range(len(R22_REFUSALS))])
def test_each_ledger_refusal(tmp_path, tries, says):
    repo, env, head = _repo(tmp_path, tries=tries)
    res = _run(repo, env, head, "r1")
    assert res.returncode == 95 and says in res.stderr, res.stderr
    _stopped_at_the_gate(tmp_path, 95)
    assert not (repo / "runs").exists()


LIFT_REFUSALS = [
    ("", "exactly one row R44"),
    (lift_row("F1") * 2, "exactly one row R44"),
    (lift_row("F1", row="R45"), "exactly one row R44"),
    (lift_row("F1", ruling="a note"), "does not say: ruling, no run"),
    (lift_row("F1", words="the 5% check gates nothing for them"), "does not say: the 5% check gates nothing for it"),
    (lift_row("F1", "F2", "F3").replace("`", ""), "names none of the three option strings"),
    (lift_row("F1").replace(f"`{OPT['F1']}`", "`--gata4-exponent inf`"), "names none of the three option strings"),
]


@pytest.mark.parametrize("lift,says", LIFT_REFUSALS, ids=[f"lift-{i}" for i in range(len(LIFT_REFUSALS))])
def test_each_lifting_row_refusal(tmp_path, lift, says):
    repo, env, head = _repo(tmp_path, tries=ROW_R22 + ROW_R41 + lift)
    res = _run(repo, env, head, "r1", "R44")
    assert res.returncode == 94 and says in res.stderr, res.stderr
    _stopped_at_the_gate(tmp_path, 94)
    assert not (repo / "runs").exists()


def test_an_odd_lifting_row_or_a_fourth_argument_is_refused_before_anything_is_made(tmp_path):
    repo, env, head = _repo(tmp_path, tries=ROW_R22 + ROW_R41 + lift_row("F1"))
    for i, bad in enumerate(["r44", "R", "R44a", "44", "R4-4", "R 44", "R44;x"]):
        res = _run(repo, env, head, f"r{i}", bad)
        assert res.returncode == 89 and "LIFT_ROW takes R and digits only" in res.stderr, (bad, res.stderr)
        assert not _scratch(tmp_path, f"r{i}").exists()
    res = _run(repo, env, head, "r9", "R44", "extra")
    assert res.returncode == 89 and "at most three arguments" in res.stderr and not _scratch(tmp_path, "r9").exists()
    assert not _calls(tmp_path) and not (repo / "runs").exists()


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
    for argv in ((head,), ()):
        res = _run(repo, env, *argv)
        assert res.returncode != 0 and USAGE in res.stderr


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t3gonly_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    res = _run(repo, env, head, "r1")
    assert res.returncode == 93 and "output folder" in res.stderr
    _stopped_at_the_gate(tmp_path, 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_that_passes_the_gate_writes_the_head_the_lift_and_its_folder(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1").returncode == 10
    s = _scratch(tmp_path)
    assert _exits(s) == {"gate": "0", "tests": "10", "all": "10"} and (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "head", "lift", "tests.exit", "tests.out",
                                                   "tmp"]
    assert (s / "lift").read_text().strip() == "none"
    o = repo / "runs" / f"t3gonly_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())
