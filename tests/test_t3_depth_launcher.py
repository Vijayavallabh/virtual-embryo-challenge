"""The test-board build launcher of the depth line: experiments/t3_decomp/t3_depth_build.sh <COMMIT> <RUN> <ROW> <ARM>.

ROW is R53 (row R53, the captain's order of 2026-11-24: his named exception on T3's depth for T3_c2_s2_dd, T3_c2_s2_ddh
and T3_c2_s2_dg, bound to row K30's read), the only row the builder takes for this line; ARM is dd, ddh, dg or all. The
launcher gates on a clean tree at COMMIT, on row R22 of the ledger (the two old lines' strings its last identity gates
need and `--captain-ruling R22`, each as one whole backtick-delimited string) and on row ROW (it opens "ROW (the
captain's order", says "ruling, no run" and names each named file's option string, `--form c2-s2` and
`--captain-ruling ROW`, each as one whole string). It runs the targeted tests, rebuilds six known files (identity-s1,
identity-s2, c1-s2, c2-s2, line 2's T3_c1_s2_g2.h5ad and line 3's T3_c2_s2_g2.h5ad) and stops unless each gives its
recorded bytes, then builds each named file twice into separate folders (the builder's own checks run inside each
build: exit 7, nothing written), compares each pair and writes one SHA256SUMS. These tests hold its strings, hashes and
step order to that, check the hashes against the committed SHA256SUMS records, check that it deletes, fetches and
uploads nothing and reads no output of a test build, and run it for real in throwaway repositories under a temporary
home with a stand-in `uv` (it logs each call and fakes the tests, the builder and `vec check`): every refusal, every
stop, whole runs through a copy of the launcher whose six recorded hashes are replaced by the stand-in's, and the gate
on row R53's first cell as the repository's ledger holds it (read by code, never printed). Nothing here builds anything
from released data.
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
LAUNCHER = REPO / "experiments" / "t3_decomp" / "t3_depth_build.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
REPO_LEDGER = REPO / "experiments" / "harness" / "TRIES.md"
REL = "experiments/t3_decomp/t3_depth_build.sh"
DEC = "experiments/t3_decomp/results/gata4/gata4_decisions.json"
WT = "data/raw/T2/E8.75.h5ad"
USAGE = "usage: t3_depth_build.sh COMMIT RUN ROW ARM"
ARMS = ("dd", "ddh", "dg")
FILES = {"dd": "T3_c2_s2_dd.h5ad", "ddh": "T3_c2_s2_ddh.h5ad", "dg": "T3_c2_s2_dg.h5ad"}
OPT = {"dd": "--depth-from allowed-s1", "ddh": "--depth-from allowed-s1 --depth-exponent 0.5",
       "dg": "--depth-from allowed-s1 --depth-form global"}
FORM = "--form c2-s2"
ARGS = {a: f"{FORM} {OPT[a]} --captain-ruling R53" for a in ARMS}
PAIR = {"dd": 31, "ddh": 32, "dg": 33}
IDENTITY = [("ids1", "identity-s1.h5ad", "IDS1_SHA", 20, "--form identity-s1"),
            ("ids2", "identity-s2.h5ad", "IDS2_SHA", 21, "--form identity-s2"),
            ("c1s2", "c1-s2.h5ad", "C1S2_SHA", 22, "--form c1-s2"),
            ("c2s2", "c2-s2.h5ad", "C2S2_SHA", 23, "--form c2-s2"),
            ("c1g2", "T3_c1_s2_g2.h5ad", "C1G2_SHA", 24, "--form c1-s2 --gata4-exponent 2.0 --captain-ruling R22"),
            ("c2g2", "T3_c2_s2_g2.h5ad", "C2G2_SHA", 25, "--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22")]
GATE22 = ["--form c2-s2 --gata4-exponent 2.0", "--gata4-exponent 2.0", "--captain-ruling R22"]
ROW_R22 = ("| R22 (the captain's order of a day; ruling, no run, not a try; by T3's second search line the builder is "
           "also to take `--gata4-exponent 2.0` on `--form c1-s2`, only with `--captain-ruling R22`; by T3's third "
           "search line also `--form c2-s2 --gata4-exponent 2.0`, only with `--captain-ruling R22`) | T3 | x |\n")
ROW_R53 = ("| R53 (the captain's order of a day, ruling, no run: \"Yes: my named exception on T3's depth\"; the T3 "
           f"two-knockout builder is to take `{OPT['dd']}`, `{OPT['ddh']}` and `{OPT['dg']}` (the options "
           "`--depth-exponent 0.5` and `--depth-form global` only beside `--depth-from allowed-s1`) on `--form c2-s2` "
           "only with `--captain-ruling R53`; not a try) | T3 | y |\n")


def steps(arms) -> list[str]:
    return ["gate", "tests", *[i[0] for i in IDENTITY], *[f"{a}{s}" for a in arms for s in "ab"], "pair", "done"]


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
    assert re.findall(r"\$\{\d+:\?([^}]*)\}", _sh()) == [USAGE] * 4
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", _sh()):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg
    sh = _sh()
    assert sh.startswith("#!/usr/bin/env bash\n") and "\nset -u\nset -o pipefail\n" in sh
    assert "grep -q" not in sh  # no pipe ends in an early-closing grep under pipefail


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
    assert want["C2S2_SHA"] == "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"  # C2's own bytes
    for sha in want.values():
        assert sh.count(sha) == 1
    # the first five as t3_c2g_build.sh's, all six as t3_gonly_build.sh's
    c2g_sh = (REPO / "experiments" / "t3_decomp" / "t3_c2g_build.sh").read_text()
    five = dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", c2g_sh, flags=re.M))
    assert five == {k: v for k, v in want.items() if k != "C2G2_SHA"}
    gonly_sh = (REPO / "experiments" / "t3_decomp" / "t3_gonly_build.sh").read_text()
    assert dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", gonly_sh, flags=re.M)) == want


def test_the_commands_and_the_step_order():
    sh = _sh()
    assert f'FORM_OPT="{FORM}"' in sh and 'DD_OPT="--depth-from allowed-s1"' in sh
    assert 'DDH_OPT="--depth-from allowed-s1 --depth-exponent 0.5"' in sh
    assert 'DG_OPT="--depth-from allowed-s1 --depth-form global"' in sh
    for a in ARMS:
        assert f"{a.upper()}={FILES[a]}" in sh
    assert f"DEC={DEC}" in sh and f"WT={WT}" in sh and "BOARD=T3:bcatenin" in sh
    assert "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only" in sh
    assert ('$PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board '
            '"$BOARD" \\\n    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?') in sh
    assert ("  if [ $rc = 2 ]; then return 12; fi\n  if [ $rc = 6 ]; then return 16; fi\n"
            "  if [ $rc = 5 ]; then return 17; fi\n  if [ $rc = 7 ]; then return 18; fi\n"
            "  if [ $rc != 0 ]; then return 11; fi\n") in sh
    assert 'uv run --frozen --no-sync vec check "$d/$name" --board "$BOARD" > /dev/null 2>&1 || return 13' in sh
    for step, name, var, code, args in IDENTITY:
        assert f'STEP={step}\nbegin\nidentity {name} "${var}" {code} {args}\nstep_ok\n' in sh
    assert ('    build "$side" "$(arm_file "$a")" $FORM_OPT $(arm_opt "$a") --captain-ruling "$ROW" || finish $?\n'
            in sh)
    assert ('  cmp -s "$O/a/$f" "$O/b/$f" || { echo "stopped: the two builds of $f differ" >&2; finish '
            '"$(arm_pair "$a")"; }') in sh
    assert "arm_pair() { case \"$1\" in dd) echo 31 ;; ddh) echo 32 ;; dg) echo 33 ;; esac; }" in sh
    assert '(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 39' in sh
    marks = ["STEP=gate", "STEP=tests", *[f"STEP={i[0]}\n" for i in IDENTITY], "for a in $ARMS; do\n  for side",
             "STEP=pair", "STEP=done"]
    order = [sh.index(m) for m in marks]
    assert order == sorted(order)
    assert '[ "$st" = "$COMMIT True" ] || return 15' in sh
    assert "O=runs/t3depth_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t3depth/$RUN" in sh
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
    assert {"tests/test_t3_depth_build.py", "tests/test_t3_depth_old_new.py", "tests/test_t3_depth_launcher.py",
            "tests/test_t3_gonly_build.py", "tests/test_t3_gonly_old_new.py", "tests/test_t3_gonly_launcher.py",
            "tests/test_t3_c2g_build.py", "tests/test_t3_c2g_old_new.py", "tests/test_t3_cd7_old_new.py",
            "tests/test_t3_two_knockout_build.py", "experiments/t3_decomp/test_build_submission.py"} <= set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert '> "$S/tests.out" 2>&1 || finish 10' in step


def test_the_gate_lines_and_the_builder_s_strings():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert '[ -z "$(git status --porcelain 2> /dev/null)" ]' in gate
    assert "[ \"$(grep -c '^| R22 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert ("for opt in '--form c2-s2 --gata4-exponent 2.0' '--gata4-exponent 2.0' '--captain-ruling R22'; do\n"
            "  if ! names \"$first\" \"$opt\"; then") in gate
    assert '[ "$(grep -c "^| $ROW " experiments/harness/TRIES.md)" = 1 ]' in gate
    assert "[[ $first == \" $ROW (the captain's order\"* ]]" in gate and '[[ $first == *"ruling, no run"* ]]' in gate
    assert ('  for opt in "$(arm_opt "$a")" "$FORM_OPT" "--captain-ruling $ROW"; do\n'
            '    if ! names "$first" "$opt"; then') in gate
    assert [int(c) for c in re.findall(r"finish (\d+)", gate)] == [96, 98, 95, 95, 94, 94, 94, 94, 93]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")
    bk = _builder()
    assert OPT == bk.DEPTH_OPTIONS and FILES == bk.DEPTH_FILES and FORM == bk.DEPTH_FORM_OPTION
    assert bk.DEPTH_ROW == "R53" and '[ "$ROW" = R53 ]' in sh
    assert GATE22 == [bk.C2G_OPTION.strip("`"), bk.G2_OPTION.strip("`"), f"--captain-ruling {bk.G2_ROW}"]
    assert (bk.EXIT_DEPTH_CHECK, bk.EXIT_GUARD_FAILED, bk.EXIT_FACTOR_MISMATCH, bk.EXIT_REFUSED) == (7, 6, 5, 2)


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path):
    bk = _builder()
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R22}{ROW_R53}")
    for arm, args in ARGS.items():
        a = args.split()
        kw = {"depth_from": a[a.index("--depth-from") + 1]}
        if "--depth-exponent" in a:
            kw["depth_exponent"] = float(a[a.index("--depth-exponent") + 1])
        if "--depth-form" in a:
            kw["depth_form"] = a[a.index("--depth-form") + 1]
        assert bk.is_depth_command(**kw)
        rec = bk.check_depth_line(a[a.index("--form") + 1], captain_ruling=a[a.index("--captain-ruling") + 1],
                                  tries=p, **kw)
        assert rec["arm"] == arm and rec["file"] == FILES[arm] and rec["row"] == "R53"
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", p)["option"] == f"`{GATE22[0]}`"  # the last identity gate's line
    assert bk.check_gata4("c1-s2", 2.0, "R22", p)["option"] == f"`{GATE22[1]}`"  # the one before


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


def _repo(tmp_path: Path, tries: str = ROW_R22 + ROW_R53, sh: str | None = None, ctl: dict | None = None
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
    return tmp_path / "home" / "v" / "scratch" / "t3depth" / run


def _exits(s: Path) -> dict[str, str]:
    return {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")}


def _calls(tmp_path: Path) -> list[str]:
    log = tmp_path / "home" / "fake" / "calls.log"
    return log.read_text().splitlines() if log.exists() else []


@pytest.mark.parametrize("arm", [*ARMS, "all"])
def test_a_whole_run_with_the_stand_in_hashes(tmp_path, arm):
    arms = ARMS if arm == "all" else (arm,)
    repo, env, head = _repo(tmp_path, sh=_with_stand_in_hashes(_sh()))
    res = _run(repo, env, head, "r1", "R53", arm)
    assert res.returncode == 0, res.stderr
    s = _scratch(tmp_path)
    assert _exits(s) == {**dict.fromkeys(steps(arms), "0"), "all": "0"}
    assert (s / "head").read_text().strip() == head
    assert (s / "arm").read_text().strip() == f"R53 {arm}: " + " ".join(FILES[a] for a in arms)
    o = repo / "runs" / f"t3depth_{head[:7]}_r1"
    sums = (s / "SHA256SUMS").read_text()
    assert (o / "SHA256SUMS").read_text() == sums
    words = sums.split()
    want = {f"gate/{name}": _sha(_content(a)) for _, name, _, _, a in IDENTITY}
    want |= {f"{side}/{FILES[a]}": _sha(_content(ARGS[a])) for a in arms for side in "ab"}
    assert dict(zip(words[1::2], words[0::2], strict=True)) == want and len(words) == 2 * (6 + 2 * len(arms))
    assert len({want[f"a/{FILES[a]}"] for a in arms}) == len(arms)  # different files, each pair one sha256
    assert (s / "stamps").read_text().splitlines() == [f"gate/{name} {head} True" for _, name, _, _, _ in IDENTITY] + [
        f"{side}/{FILES[a]} {head} True" for a in arms for side in "ab"]
    run_dir = f"{tmp_path / 'home'}/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
    calls = _calls(tmp_path)
    assert calls[0].startswith("run --frozen --no-sync python -I -B -m pytest -q -p no:cacheprovider")
    fixed = f"--wt {WT} --run-dir {run_dir} --board T3:bcatenin --seed 1000001"
    tool = [c for c in calls if "build_two_knockout.py" in c]
    want_tool = [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} {fixed} "
                 f"{a} --out runs/t3depth_{head[:7]}_r1/gate/{name}" for _, name, _, _, a in IDENTITY]
    want_tool += [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} "
                  f"{fixed} {ARGS[a]} --out runs/t3depth_{head[:7]}_r1/{side}/{FILES[a]}"
                  for a in arms for side in "ab"]
    assert tool == want_tool
    assert not any("--depth" in c for c in tool[:6])  # the identity gates never take the line's options
    assert len([c for c in calls if c.startswith("run --frozen --no-sync vec check")]) == 6 + 2 * len(arms)
    assert sorted(p.name for p in o.iterdir()) == ["SHA256SUMS", "a", "b", "gate"]
    assert sorted(p.name for p in (o / "a").iterdir() if p.suffix == ".h5ad") == sorted(FILES[a] for a in arms)


@pytest.mark.parametrize("ctl,step,code", [
    ({"tests_exit": 1}, "tests", 10),
    ({"fail": "*/gate/identity-s1.h5ad 1"}, "ids1", 11),
    ({"fail": "*/gate/identity-s1.h5ad 2"}, "ids1", 12),
    ({"fail": "*/gate/T3_c1_s2_g2.h5ad 5"}, "c1g2", 17),  # line 2's own factor check
    ({"fail": "*/gate/T3_c2_s2_g2.h5ad 6"}, "c2g2", 16),  # line 3's own 5% guard
    ({"check_exit": 1}, "ids1", 13),
    ({"unclean": "*/gate/c2-s2.h5ad"}, "c2s2", 15),
    ({"fail": "*/a/T3_c2_s2_dd.h5ad 7"}, "dda", 18),  # the depth line's checks: nothing written
    ({"fail": "*/b/T3_c2_s2_dd.h5ad 2"}, "ddb", 12),
    ({"fail": "*/a/T3_c2_s2_ddh.h5ad 7"}, "ddha", 18),
    ({"fail": "*/b/T3_c2_s2_ddh.h5ad 1"}, "ddhb", 11),
    ({"unclean": "*/a/T3_c2_s2_dg.h5ad"}, "dga", 15),
    ({"fail": "*/b/T3_c2_s2_dg.h5ad 7"}, "dgb", 18),
    ({"differ": "*/b/T3_c2_s2_dd.h5ad"}, "pair", 31),
    ({"differ": "*/b/T3_c2_s2_ddh.h5ad"}, "pair", 32),
    ({"differ": "*/b/T3_c2_s2_dg.h5ad"}, "pair", 33),
])
def test_each_stop_after_the_gate(tmp_path, ctl, step, code):
    repo, env, head = _repo(tmp_path, sh=_with_stand_in_hashes(_sh()), ctl=ctl)
    res = _run(repo, env, head, "r1", "R53", "all")
    assert res.returncode == code, res.stderr
    s = _scratch(tmp_path)
    ex = _exits(s)
    all_steps = steps(ARMS)
    assert ex.pop("all") == str(code) and ex.pop(step) == str(code)
    assert ex == dict.fromkeys(all_steps[:all_steps.index(step)], "0")  # every earlier step passed, no later began
    assert not (s / "SHA256SUMS").exists()
    if code == 18:
        side, f = step[-1], FILES[step[:-1]]
        assert not (repo / "runs" / f"t3depth_{head[:7]}_r1" / side / f).exists()


def test_one_arm_stops_only_on_its_own_builds(tmp_path):
    """ARM dg: dd's and ddh's stand-in failures never fire, and the run builds dg alone."""
    ctl = {"fail": "*/a/T3_c2_s2_dd.h5ad 7"}
    repo, env, head = _repo(tmp_path, sh=_with_stand_in_hashes(_sh()), ctl=ctl)
    res = _run(repo, env, head, "r1", "R53", "dg")
    assert res.returncode == 0, res.stderr
    assert _exits(_scratch(tmp_path)) == {**dict.fromkeys(steps(("dg",)), "0"), "all": "0"}


def test_the_launcher_s_own_hashes_stop_the_stand_in_at_the_first_identity_gate(tmp_path):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, head, "r1", "R53", "dd")
    assert res.returncode == 20 and "gate/identity-s1.h5ad is not its recorded bytes" in res.stderr
    assert _exits(_scratch(tmp_path)) == {"gate": "0", "tests": "0", "ids1": "20", "all": "20"}


def test_the_c2_identity_gate_stops_on_bytes_other_than_c2_s_record(tmp_path):
    """The stand-in's hashes for every gate but C2's: C2 (c2-s2.h5ad) keeps its recorded sha256 47111ac1d07e..."""
    sh = _with_stand_in_hashes(_sh())
    sh = re.sub(r"^C2S2_SHA=[0-9a-f]{64}", "C2S2_SHA=47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7",
                sh, flags=re.M)
    repo, env, head = _repo(tmp_path, sh=sh)
    res = _run(repo, env, head, "r1", "R53", "all")
    assert res.returncode == 23 and "gate/c2-s2.h5ad is not its recorded bytes" in res.stderr
    ex = _exits(_scratch(tmp_path))
    all_steps = steps(ARMS)
    assert ex.pop("c2s2") == "23" and ex.pop("all") == "23"
    assert ex == dict.fromkeys(all_steps[:all_steps.index("c2s2")], "0")


def _stopped_at_the_gate(tmp_path: Path, code: int, run: str = "r1") -> None:
    s = _scratch(tmp_path, run)
    assert _exits(s) == {"gate": str(code), "all": str(code)}
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "tmp"]
    assert not list((s / "tmp").iterdir()) and not _calls(tmp_path)


R22_REFUSALS = [
    (ROW_R53, "exactly one row R22"),
    (ROW_R22 + ROW_R22 + ROW_R53, "exactly one row R22"),
    (ROW_R22.replace("| R22 (", "| R22b (") + ROW_R53, "exactly one row R22"),
    *[(ROW_R22.replace(f"`{o}`", o) + ROW_R53, f"does not name `{o}` as one whole string") for o in GATE22],
]


@pytest.mark.parametrize("tries,says", R22_REFUSALS, ids=[f"r22-{i}" for i in range(len(R22_REFUSALS))])
def test_each_row_r22_refusal(tmp_path, tries, says):
    repo, env, head = _repo(tmp_path, tries=tries)
    res = _run(repo, env, head, "r1", "R53", "all")
    assert res.returncode == 95 and says in res.stderr, res.stderr
    _stopped_at_the_gate(tmp_path, 95)
    assert not (repo / "runs").exists()


R53_REFUSALS = [
    ("", "all", "exactly one row R53"),
    (ROW_R53, "all", None),  # two rows R53 (the ledger's own and this one)
    (ROW_R53.replace("| R53 (", "| R53b ("), "all", "exactly one row R53"),
    (ROW_R53.replace("the captain's order", "the captain's named exception"), "all", "does not open"),
    (ROW_R53.replace("the captain's order", "the captain's ruling"), "all", "does not open"),
    (ROW_R53.replace("ruling, no run", "a run"), "all", "does not say: ruling, no run"),
    (ROW_R53.replace("`--form c2-s2`", "--form c2-s2"), "dd", "does not name `--form c2-s2` as one whole string"),
    (ROW_R53.replace("`--captain-ruling R53`", "R53"), "dg",
     "does not name `--captain-ruling R53` as one whole string"),
    (ROW_R53.replace(f"`{OPT['ddh']}`", OPT["ddh"]), "ddh", f"does not name `{OPT['ddh']}` as one whole string"),
    (ROW_R53.replace(f"`{OPT['ddh']}`", OPT["ddh"]), "all", f"does not name `{OPT['ddh']}` as one whole string"),
    (ROW_R53.replace(f"`{OPT['dg']}`", f"`a`{OPT['dg']}`b`"), "dg", f"does not name `{OPT['dg']}` as one whole"),
    (ROW_R53.replace(f"`{OPT['dg']}`", "`--depth-form global`"), "dg", f"does not name `{OPT['dg']}`"),
    ("| R53 (the captain's order of a day, ruling, no run; not a try) | T3 | " + " ".join(f"`{o}`" for o in
                                                                                      OPT.values())
     + f" `{FORM}` `--captain-ruling R53` |\n", "dd", f"does not name `{OPT['dd']}`"),  # outside the first cell
]


@pytest.mark.parametrize("r53,arm,says", R53_REFUSALS, ids=[f"r53-{i}" for i in range(len(R53_REFUSALS))])
def test_each_row_r53_refusal(tmp_path, r53, arm, says):
    tries = ROW_R22 + (ROW_R53 + r53 if says is None else r53)
    repo, env, head = _repo(tmp_path, tries=tries)
    res = _run(repo, env, head, "r1", "R53", arm)
    assert res.returncode == 94 and (says or "exactly one row R53") in res.stderr, res.stderr
    _stopped_at_the_gate(tmp_path, 94)
    assert not (repo / "runs").exists()


def test_a_row_naming_one_file_opens_that_file_only(tmp_path):
    one = ("| R53 (the captain's order of a day, ruling, no run; the builder is to take `--depth-from allowed-s1 "
           "--depth-form global` on `--form c2-s2` only with `--captain-ruling R53`) | T3 | y |\n")
    repo, env, head = _repo(tmp_path, tries=ROW_R22 + one, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1", "R53", "dg").returncode == 10  # past the gate
    for i, arm in enumerate(("dd", "ddh", "all")):
        res = _run(repo, env, head, f"x{i}", "R53", arm)
        assert res.returncode == 94 and "does not name" in res.stderr, (arm, res.stderr)


def test_the_gate_passes_on_row_r53_as_the_repository_s_ledger_holds_it(tmp_path):
    """Row R53's first cell read by code from the repository's ledger (never printed), set in a throwaway ledger
    beside a synthetic row R22: the gate passes for every ARM and the run goes on to the tests step."""
    lines = [ln for ln in REPO_LEDGER.read_text(encoding="utf-8").splitlines() if ln.startswith("| R53 ")]
    assert len(lines) == 1, "the repository's ledger holds no single row R53"
    first = lines[0][2:].split("|", 1)[0]
    repo, env, head = _repo(tmp_path, tries=ROW_R22 + f"| {first}| T3 | x |\n", ctl={"tests_exit": 1})
    for i, arm in enumerate([*ARMS, "all"]):
        res = _run(repo, env, head, f"r{i}", "R53", arm)
        assert res.returncode == 10, (arm, res.stderr)
        assert _exits(_scratch(tmp_path, f"r{i}")) == {"gate": "0", "tests": "10", "all": "10"}


def test_an_odd_arm_row_or_argument_count_is_refused_before_anything_is_made(tmp_path):
    repo, env, head = _repo(tmp_path)
    for i, bad in enumerate(["", "DD", "d", "dd ", "ddh,dg", "both", "dd;x", "a"]):
        res = _run(repo, env, head, f"a{i}", "R53", bad)
        if bad == "":
            assert res.returncode != 0 and USAGE in res.stderr
        else:
            assert res.returncode == 88 and "ARM takes dd, ddh, dg or all" in res.stderr, (bad, res.stderr)
        assert not _scratch(tmp_path, f"a{i}").exists()
    for i, bad in enumerate(["R22", "R36", "R54", "r53", "R53 ", "53", "R530"]):
        res = _run(repo, env, head, f"b{i}", bad, "dd")
        assert res.returncode == 89 and "ROW is R53" in res.stderr, (bad, res.stderr)
        assert not _scratch(tmp_path, f"b{i}").exists()
    res = _run(repo, env, head, "c1", "R53", "dd", "extra")
    assert res.returncode == 89 and "four arguments" in res.stderr and not _scratch(tmp_path, "c1").exists()
    for argv in ((head, "r9", "R53"), (head, "r9"), (head,), ()):
        res = _run(repo, env, *argv)
        assert res.returncode != 0 and USAGE in res.stderr and not _scratch(tmp_path, "r9").exists()
    assert not _calls(tmp_path) and not (repo / "runs").exists()


def test_a_wrong_commit_a_prefix_and_an_unclean_tree_are_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, "0" * 40, "r1", "R53", "all")
    assert res.returncode == 96 and "HEAD is not" in res.stderr
    _stopped_at_the_gate(tmp_path, 96)
    assert _run(repo, env, head[:12], "r2", "R53", "all").returncode == 96
    _stopped_at_the_gate(tmp_path, 96, "r2")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R22}{ROW_R53}more\n")
    res = _run(repo, env, head, "r3", "R53", "all")
    assert res.returncode == 98 and "not clean" in res.stderr
    _stopped_at_the_gate(tmp_path, 98, "r3")
    _git(repo, env, "checkout", "-q", "--", "experiments/harness/TRIES.md")
    (repo / "stray.py").write_text("x = 1\n")
    assert _run(repo, env, head, "r4", "R53", "dd").returncode == 98
    _stopped_at_the_gate(tmp_path, 98, "r4")
    assert not (repo / "runs").exists()


def test_a_run_name_is_never_reused_and_an_odd_one_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1", "R53", "all").returncode == 10
    s = _scratch(tmp_path)
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "R53", "all")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    odd = _run(repo, env, head, "a/b", "R53", "all")
    assert odd.returncode == 90 and not _scratch(tmp_path, "a").exists()


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t3depth_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    res = _run(repo, env, head, "r1", "R53", "all")
    assert res.returncode == 93 and "output folder" in res.stderr
    _stopped_at_the_gate(tmp_path, 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_that_passes_the_gate_writes_the_head_the_arm_and_its_folder(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1", "R53", "ddh").returncode == 10
    s = _scratch(tmp_path)
    assert _exits(s) == {"gate": "0", "tests": "10", "all": "10"} and (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "arm", "gate.exit", "head", "tests.exit", "tests.out",
                                                   "tmp"]
    assert (s / "arm").read_text().strip() == "R53 ddh: T3_c2_s2_ddh.h5ad"
    o = repo / "runs" / f"t3depth_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())
