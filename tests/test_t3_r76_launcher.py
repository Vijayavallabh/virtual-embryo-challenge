"""The test-board build launcher of row R76's five files:
experiments/t3_decomp/t3_r76_line.sh <COMMIT> <RUN> <FILE_NAME>.

It refuses any name but the five of its literal table, and any odd argument, before it makes a folder; its gate checks
the commit, the clean tracked tree, row R76's first cell (the file's option string and `--captain-ruling R76`, and
`--den-ruling R76` for the three lifted strings, as whole backticked strings, the string paired there with the file's
name), the data and the decision file, the gata4 run folder (its two tables at the builder's pins) and the venv; then it
runs the targeted tests, rebuilds C2, T3_c2_s2_g2, T3_c2_s2_s067 and T3_c2_s2_gonly and stops unless each gives its
recorded bytes, then builds the file twice into separate folders from one clean commit, writes one sha256 list and
compares the two. These tests hold its table to the builder's, its hashes to the committed SHA256SUMS records and its
identity commands to the earlier launchers', check that it deletes, fetches and uploads nothing and reads no output of a
test build but exit codes, sha256 and the stamps' commit and clean fields, and run it for real in throwaway
repositories under a temporary home with a stand-in `uv` (it logs each call and fakes the tests, the builder and
`vec check`): every refusal, every stop, whole runs of each of the five names through a copy of the launcher whose pins
are the stand-ins', and the gate on row R76's first cell as the repository's ledger holds it (read by code, never
printed). Nothing here builds anything from released data.
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
LAUNCHER = REPO / "experiments" / "t3_decomp" / "t3_r76_line.sh"
GONLY_LAUNCHER = REPO / "experiments" / "t3_decomp" / "t3_gonly_build.sh"
DEPTH_LAUNCHER = REPO / "experiments" / "t3_decomp" / "t3_depth_build.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
REPO_LEDGER = REPO / "experiments" / "harness" / "TRIES.md"
REL = "experiments/t3_decomp/t3_r76_line.sh"
DEC = "experiments/t3_decomp/results/gata4/gata4_decisions.json"
WT = "data/raw/T2/E8.75.h5ad"
RUN_DIR_REL = "v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
USAGE = "usage: t3_r76_line.sh COMMIT RUN FILE_NAME"
TABLE = (
    ("T3_c2_s2_s085", "--form c2-s2 --size-scale 0.85", True),
    ("T3_c2_s2_s05", "--form c2-s2 --size-scale 0.5", True),
    ("T3_c2_s2_g15", "--form c2-s2 --gata4-exponent 1.5", False),
    ("T3_c2_s2_g05", "--form c2-s2 --gata4-exponent 0.5", False),
    ("T3_c2_s2_g2_s067", "--form c2-s2 --gata4-exponent 2.0 --size-scale 0.67", True),
)
NAMES = [n for n, _, _ in TABLE]
OPT = {n: o for n, o, _ in TABLE}
LIFTED = [n for n, _, lifted in TABLE if lifted]
WEIGHT_ONLY = [n for n, _, lifted in TABLE if not lifted]
ARGS = {n: f"{OPT[n]} --captain-ruling R76" + (" --den-ruling R76" if n in LIFTED else "") for n in NAMES}
# identity step -> (file name in gate/, hash variable, exit code, the builder's options)
IDENTITY = [("c2", "c2-s2.h5ad", "C2S2_SHA", 21, "--form c2-s2"),
            ("c2g2", "T3_c2_s2_g2.h5ad", "C2G2_SHA", 22, "--form c2-s2 --gata4-exponent 2.0 --captain-ruling R22"),
            ("s067", "T3_c2_s2_s067.h5ad", "S067_SHA", 23,
             "--form c2-s2 --size-scale 0.67 --captain-ruling R22 --den-ruling R47"),
            ("gonly", "T3_c2_s2_gonly.h5ad", "GONLY_SHA", 24,
             "--form c2-s2 --gata4-exponent inf --captain-ruling R22 --den-ruling R47")]
TABLES = {"TABLE_M_SHA": ("transfer_mab21l2_all.csv", "gene,factor\nG0,2.0\n"),
          "TABLE_G_SHA": ("gata4_pooled.csv", "gene,factor\nG0,0.5\n")}
FIRST = "T3_c2_s2_s085"  # the file the gate's refusals are run on


def r76_row(names=tuple(NAMES), lifted=tuple(LIFTED), opener="R76 (the captain's order of a day",
            ruling="ruling, no run", captain="`--captain-ruling R76`", den="`--den-ruling R76`", pairs=None) -> str:
    pairs = pairs if pairs is not None else ", ".join(f"`{OPT[n]}` ({n})" for n in names)
    return (f"| {opener}, {ruling}: \"Yes to the five files under R76\"; T3's refinement round on C2's command; the "
            f"T3 two-knockout builder is to take exactly these strings, each only with {captain}: {pairs}; for each "
            "of " + ", ".join(f"`{OPT[n]}`" for n in lifted) + f", the 5% check gates nothing for it ({den}); not a "
            "try) | T3 | z |\n")


ROW = r76_row()


def steps() -> list[str]:
    return ["gate", "tests", *[i[0] for i in IDENTITY], "a", "b", "done"]


# a stand-in for uv: it logs each call and fakes pytest, the builder (its exit set per output path) and vec check
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


def _code(sh: str | None = None) -> list[str]:
    return [line.split("#")[0] for line in (sh or _sh()).splitlines()]


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


def _pins(sh: str) -> dict[str, str]:
    return dict(re.findall(r"^([A-Z0-9_]+_SHA)=([0-9a-f]{64})  #", sh, flags=re.M))


# ------------------------------------------------------------------ what the launcher runs


def test_the_launcher_parses_and_deletes_fetches_and_uploads_nothing():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True, check=False)
    assert out.returncode == 0, out.stderr
    code = _code()
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)(\s|$)", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash|fetch|pull|push|clone)", line)
                   for line in code)
    assert not any(re.search(r"(^|[\s;&|(])(curl|wget|scp|rsync|ssh|sftp|ftp)(\s|$)", line) for line in code)
    assert not any(re.search(r"upload|submit|format_only|token", line, re.I) for line in code)
    sh = _sh()
    assert re.findall(r"\$\{\d+:\?([^}]*)\}", sh) == [USAGE] * 3
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg
    assert sh.startswith("#!/usr/bin/env bash\n") and "\nset -u\nset -o pipefail\n" in sh
    assert "grep -q" not in sh  # no pipe ends in an early-closing grep under pipefail


def test_the_launchers_literal_table_is_the_builders_one_table():
    sh = _sh()
    rows = re.findall(r'^  (T3_c2_s2_[A-Za-z0-9_]+)\) OPT="([^"]+)"; LIFT=(yes|no) ;;$', sh, flags=re.M)
    table = tuple((n, o, lift == "yes") for n, o, lift in rows)
    assert table == TABLE == _builder().R76_TABLE
    case = sh[sh.index('case "$FILE" in'):sh.index("esac", sh.index('case "$FILE" in'))]
    patterns = [line.split(")")[0].strip() for line in case.splitlines()[1:]]
    assert patterns == [*NAMES, "*"]  # a literal table: no name is read by a pattern, the default refuses
    assert '  *) echo "refused: FILE_NAME is not one of row R76\'s five files" >&2; exit 89 ;;' in case
    assert 'ARGS="$OPT --captain-ruling R76"\nif [ "$LIFT" = yes ]; then ARGS="$ARGS --den-ruling R76"; fi\n' in sh
    assert "BASH_REMATCH" not in sh and "=~" not in sh
    assert sh.index("exit 89 ;;") < sh.index('mkdir -p "$HOME/v/scratch/t3r76"')  # every refusal before any folder


def test_the_recorded_hashes_are_the_committed_build_records_and_the_table_pins_the_builder_s():
    sh = _sh()
    two, c2g, gonly = (_sums("t3two_fb01487.SHA256SUMS"), _sums("t3c2g_fc57ca9_g1.SHA256SUMS"),
                       _sums("t3gonly_466e94d_g1.SHA256SUMS"))
    for d, name in ((c2g, "T3_c2_s2_g2.h5ad"), (gonly, "T3_c2_s2_s067.h5ad"), (gonly, "T3_c2_s2_gonly.h5ad")):
        assert d[f"a/{name}"] == d[f"b/{name}"], name  # each a pair of one sha256
    bk = _builder()
    want = {"C2S2_SHA": two["c2-s2.h5ad"], "C2G2_SHA": c2g["a/T3_c2_s2_g2.h5ad"],
            "S067_SHA": gonly["a/T3_c2_s2_s067.h5ad"], "GONLY_SHA": gonly["a/T3_c2_s2_gonly.h5ad"],
            "TABLE_M_SHA": bk.TABLE_SHA256[bk.TABLE_M], "TABLE_G_SHA": bk.TABLE_SHA256[bk.TABLE_G]}
    assert _pins(sh) == want
    assert want["C2S2_SHA"] == bk.C2_S2_SHA256 == "47111ac1d07e73c0c26e990a474dff0de94ae22d6a4474a569bcb674ae2b85a7"
    for sha in want.values():
        assert sh.count(sha) == 1
    for other in (GONLY_LAUNCHER, DEPTH_LAUNCHER):  # C2 and T3_c2_s2_g2 as the earlier launchers gate them
        pins = _pins(other.read_text())
        assert pins["C2S2_SHA"] == want["C2S2_SHA"] and pins["C2G2_SHA"] == want["C2G2_SHA"]


def test_the_identity_commands_are_the_earlier_launchers_and_the_builder_takes_every_command(tmp_path):
    sh = _sh()
    gonly_sh = GONLY_LAUNCHER.read_text()
    assert 'F2_OPT="--form c2-s2 --size-scale 0.67"' in gonly_sh and 'F1_OPT="--form c2-s2 --gata4-exponent inf"' in (
        gonly_sh)
    assert 'F2_ARGS="$F2_OPT --captain-ruling R22"' in gonly_sh and '"$F2_ARGS --den-ruling $LIFT"' in gonly_sh
    assert "identity T3_c2_s2_g2.h5ad \"$C2G2_SHA\" 25 --form c2-s2 --gata4-exponent 2.0 --captain-ruling R22" in (
        gonly_sh)
    args = dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"C2G2": IDENTITY[1][4], "S067": IDENTITY[2][4], "GONLY": IDENTITY[3][4]}
    bk = _builder()
    rows = ("| R22 (the captain's order of a day; ruling, no run; also `--form c2-s2 --gata4-exponent 2.0`, "
            "`--form c2-s2 --gata4-exponent inf` and `--form c2-s2 --size-scale 0.67`, only with "
            "`--captain-ruling R22`) | T3 | x |\n"
            "| R47 (the captain's order of a day, ruling, no run; for each of `--form c2-s2 --gata4-exponent inf` "
            "and `--form c2-s2 --size-scale 0.67`, the 5% check gates nothing for it) | T3 | y |\n" + ROW)
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{rows}")
    assert bk.check_c2_gata4("c2-s2", 2.0, "R22", p)["option"] == "`--form c2-s2 --gata4-exponent 2.0`"
    assert bk.check_c2_gonly("c2-s2", None, 0.67, "R22", p, den_ruling="R47")["name"] == "F2"
    assert bk.check_c2_gonly("c2-s2", float("inf"), None, "R22", p, den_ruling="R47")["name"] == "F1"
    for name in NAMES:  # each file's command as the launcher passes it, read back from its words
        a = ARGS[name].split()
        kw = {"gata4_exponent": float(a[a.index("--gata4-exponent") + 1]) if "--gata4-exponent" in a else None,
              "size_scale": float(a[a.index("--size-scale") + 1]) if "--size-scale" in a else None,
              "den_ruling": a[a.index("--den-ruling") + 1] if "--den-ruling" in a else None}
        rec = bk.check_r76(a[a.index("--form") + 1], captain_ruling=a[a.index("--captain-ruling") + 1], tries=p, **kw)
        assert rec["name"] == name and rec["file"] == f"{name}.h5ad"


def test_the_commands_and_the_step_order():
    sh = _sh()
    assert f"DEC={DEC}" in sh and f"WT={WT}" in sh and "BOARD=T3:bcatenin" in sh
    assert "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z  # read only" in sh
    assert ('$PY experiments/t3_decomp/build_two_knockout.py build "$DEC" --wt "$WT" --run-dir "$RUN_DIR" --board '
            '"$BOARD" \\\n    --seed 1000001 "$@" --out "$d/$name" > /dev/null 2>&1 || rc=$?') in sh
    assert ("  if [ $rc = 2 ]; then return 12; fi\n  if [ $rc = 6 ]; then return 16; fi\n"
            "  if [ $rc = 5 ]; then return 17; fi\n  if [ $rc != 0 ]; then return 11; fi\n") in sh
    assert 'uv run --frozen --no-sync vec check "$d/$name" --board "$BOARD" > /dev/null 2>&1 || return 13' in sh
    assert '[ "$st" = "$COMMIT True" ] || return 15' in sh
    for step, name, var, code, _ in IDENTITY:
        tail = "--form c2-s2" if step == "c2" else f"${step.upper()}_ARGS"
        assert f'identity {step} {name} "${var}" {code} {tail}\n' in sh
    assert 'for STEP in a b; do begin; build "$STEP" "$FILE.h5ad" $ARGS || finish $?; step_ok; done' in sh
    assert '(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 40' in sh
    assert 'cp "$S/SHA256SUMS" "$O/SHA256SUMS" || finish 40' in sh
    assert 'cmp -s "$O/a/$FILE.h5ad" "$O/b/$FILE.h5ad" || {' in sh and "finish 41; }" in sh
    marks = ["STEP=gate", "STEP=tests", "identity c2 ", "identity c2g2 ", "identity s067 ", "identity gonly ",
             "for STEP in a b; do", "STEP=done"]
    order = [sh.index(m) for m in marks]
    assert order == sorted(order)
    assert "O=runs/t3r76_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t3r76/$RUN" in sh
    assert "export NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < sh.index("uv run")
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("identity c2 ")


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    for line in sh.splitlines():
        if re.search(r"build_two_knockout\.py build|vec check \"", line) and not line.lstrip().startswith("#"):
            assert "> /dev/null 2>&1" in line or line.rstrip().endswith("\\"), line
    code = "\n".join(_code())
    for word in (".build.json", "du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -",
                 "tail ", "cat ", "jq "):
        assert word not in code, word
    assert ".json" not in sh.replace("$1.stamp.json", "").replace("gata4_decisions.json", "")
    assert 'print(s["commit"], s["clean"])' in sh
    cds = [line for line in _code() if re.search(r"(^|[\s;&|(])cd\s", line)]
    assert cds == ['(cd "$O" && sha256sum gate/*.h5ad a/*.h5ad b/*.h5ad) > "$S/SHA256SUMS" || finish 40']
    assert " g1 " not in sh and " g2 " not in sh  # row R36's gates are not this line's


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh and 'step_ok() { echo 0 > "$S/$STEP.exit"; }' in sh
    assert 'echo running > "$S/all.exit"\ntrap abort ERR TERM INT HUP EXIT\nbegin\n' in sh
    for step in ("tests", "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step
    ident = sh[sh.index("identity() {"):sh.index("identity c2 ")]
    assert ident.index("STEP=$1") < ident.index("begin") < ident.index("build gate") < ident.index("step_ok")
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("begin\n\nhead="):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert 'n76=$(grep -c "^| R76 " experiments/harness/TRIES.md 2> /dev/null) || n76=0' in gate
    assert "first=$(grep \"^| R76 \" experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert "[[ $first == \" R76 (the captain's \"* ]]" in gate and '[[ $first == *"ruling, no run"* ]]' in gate
    assert 'NAMED=("$OPT" "--captain-ruling R76")\nif [ "$LIFT" = yes ]; then NAMED+=("--den-ruling R76"); fi' in gate
    assert '[[ $first == *"$BT$OPT$BT ($FILE)"* ]]' in gate  # the string paired with the file's name
    assert '[ -f "$WT" ] && [ -f "$DEC" ]' in gate
    assert ('[ -d "$RUN_DIR" ] && pinned "$RUN_DIR/factors/transfer_mab21l2_all.csv" "$TABLE_M_SHA" \\\n'
            '  && pinned "$RUN_DIR/factors/gata4_pooled.csv" "$TABLE_G_SHA" \\\n') in gate
    assert "for pth in .venv/lib/python3*/site-packages/vec.pth; do" in gate
    assert '[ "$(readlink -f "$(< "$pth")")" = "$(readlink -f src)" ]' in gate
    assert [int(c) for c in re.findall(r"finish (\d+)", gate)] == [96, 98, 95, 95, 95, 95, 95, 92, 87, 88, 93]


def test_the_tests_step_lists_existing_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"(?:tests|experiments/t3_decomp)/test_[a-z0-9_]+\.py", step)
    depth = DEPTH_LAUNCHER.read_text()
    depth_step = depth[depth.index("STEP=tests"):depth.index("step_ok", depth.index("STEP=tests"))]
    earlier = re.findall(r"(?:tests|experiments/t3_decomp)/test_[a-z0-9_]+\.py", depth_step)
    assert {"tests/test_t3_r76_build.py", "tests/test_t3_r76_old_new.py", "tests/test_t3_r76_launcher.py",
            *earlier} == set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert '> "$S/tests.out" 2>&1 || finish 10' in step and '--basetemp="$S/tmp/pytest"' in step


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    code = "\n".join(_code())  # comments dropped
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", code)}
    used |= {int(c) for c in re.findall(r'identity \S+ \S+ "\$\w+" (\d+)', code)}
    used -= {0, 1}  # finish 0 at the end; a helper's own false
    documented = {int(c) for c in re.findall(r"(?<![\w.])(\d\d)(?= )", doc)}
    assert used == documented, (used, documented)
    assert 99 in used and {21, 22, 23, 24} <= used


# ------------------------------------------------------------------ the launcher, run for real


def _content(args: str) -> bytes:
    return f"fake build {args}\n".encode()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _stand_in(sh: str, tables: bool = True, ids: bool = True, keep: tuple[str, ...] = ()) -> str:
    """The launcher with its table pins and/or its identity hashes set to the stand-ins' (those in keep left)."""
    out = sh
    subs = {}
    if tables:
        subs |= {var: _sha(text.encode()) for var, (_, text) in TABLES.items()}
    if ids:
        subs |= {var: _sha(_content(args)) for _, _, var, _, args in IDENTITY}
    for var, sha in subs.items():
        if var in keep:
            continue
        out, n = re.subn(rf"^{var}=[0-9a-f]{{64}}", f"{var}={sha}", out, flags=re.M)
        assert n == 1, var
    changed = [a for a, b in zip(sh.splitlines(), out.splitlines(), strict=True) if a != b]
    assert all(re.match(r"^[A-Z0-9_]+_SHA=", ln) for ln in changed)
    return out


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env |= {"HOME": str(tmp_path / "home"), "GIT_CEILING_DIRECTORIES": str(tmp_path), "FAKE_PY": sys.executable,
            "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW, sh: str | None = None, ctl: dict | None = None, data: bool = True,
          dec: bool = True, run_folder: str | None = "stand-in", venv: str | None = "own"
          ) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (or the given copy of it, by default the one with stand-in pins)
    and the files its gate reads, and a home whose ~/.local/bin holds the stand-in `uv` and whose ~/v holds the gata4
    run folder with stand-in tables; returns the repository, the env and its HEAD."""
    env = _env(tmp_path)
    home = tmp_path / "home"
    fake = home / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text(FAKE_UV)
    (fake / "uv").chmod(0o755)
    (home / "fake").mkdir(exist_ok=True)
    for name, value in (ctl or {}).items():
        (home / "fake" / name).write_text(f"{value}\n")
    if run_folder is not None:
        factors = home / RUN_DIR_REL / "factors"
        factors.mkdir(parents=True, exist_ok=True)
        for _, (fname, text) in TABLES.items():
            (factors / fname).write_text(text if run_folder == "stand-in" else text + "changed\n")
    repo = tmp_path / "repo"
    (repo / "experiments" / "t3_decomp").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "src").mkdir()
    (repo / "src" / "keep.txt").write_text("src\n")
    (repo / REL).write_text(_stand_in(_sh()) if sh is None else sh)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    if dec:
        (repo / DEC).parent.mkdir(parents=True)
        (repo / DEC).write_text("{}\n")
    if data:
        (repo / WT).parent.mkdir(parents=True)
        (repo / WT).write_text("stand-in\n")
    if venv is not None:
        site = repo / ".venv" / "lib" / "python3.12" / "site-packages"
        site.mkdir(parents=True)
        (site / "vec.pth").write_text(f"{repo / 'src' if venv == 'own' else tmp_path / 'other' / 'src'}\n")
    (repo / ".gitignore").write_text("data/\nruns/\n.venv/\n")
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", REL, *argv], cwd=repo, env=env, capture_output=True, text=True, timeout=120,
                          check=False)


def _scratch(tmp_path: Path, run: str = "r1") -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t3r76" / run


def _exits(s: Path) -> dict[str, str]:
    return {p.name[:-5]: p.read_text().strip() for p in s.glob("*.exit")}


def _calls(tmp_path: Path) -> list[str]:
    log = tmp_path / "home" / "fake" / "calls.log"
    return log.read_text().splitlines() if log.exists() else []


@pytest.mark.parametrize("name", NAMES)
def test_a_whole_run_with_the_stand_in_pins(tmp_path, name):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, head, "r1", name)
    assert res.returncode == 0, res.stderr
    s = _scratch(tmp_path)
    assert _exits(s) == {**dict.fromkeys(steps(), "0"), "all": "0"}
    assert (s / "head").read_text().strip() == head
    o = repo / "runs" / f"t3r76_{head[:7]}_r1"
    sums = (s / "SHA256SUMS").read_text()
    assert (o / "SHA256SUMS").read_text() == sums
    words = sums.split()
    want = {f"gate/{n}": _sha(_content(a)) for _, n, _, _, a in IDENTITY}
    want |= {f"{side}/{name}.h5ad": _sha(_content(ARGS[name])) for side in "ab"}
    assert dict(zip(words[1::2], words[0::2], strict=True)) == want and len(words) == 2 * (4 + 2)
    assert (s / "stamps").read_text().splitlines() == [f"gate/{n} {head} True" for _, n, _, _, _ in IDENTITY] + [
        f"{side}/{name}.h5ad {head} True" for side in "ab"]
    run_dir = f"{tmp_path / 'home'}/{RUN_DIR_REL}"
    calls = _calls(tmp_path)
    assert calls[0].startswith("run --frozen --no-sync python -I -B -m pytest -q -p no:cacheprovider")
    fixed = f"--wt {WT} --run-dir {run_dir} --board T3:bcatenin --seed 1000001"
    tool = [c for c in calls if "build_two_knockout.py" in c]
    r = f"runs/t3r76_{head[:7]}_r1"
    want_tool = [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} {fixed} "
                 f"{a} --out {r}/gate/{n}" for _, n, _, _, a in IDENTITY]
    want_tool += [f"run --frozen --no-sync python -I -B experiments/t3_decomp/build_two_knockout.py build {DEC} "
                  f"{fixed} {ARGS[name]} --out {r}/{side}/{name}.h5ad" for side in "ab"]
    assert tool == want_tool
    assert not any("R76" in c for c in tool[:4])  # the identity rebuilds never take row R76's ruling
    assert ("--den-ruling R76" in tool[4]) is (name in LIFTED)
    assert len([c for c in calls if c.startswith("run --frozen --no-sync vec check")]) == 6
    assert sorted(p.name for p in o.iterdir()) == ["SHA256SUMS", "a", "b", "gate"]
    assert sorted(p.name for p in (o / "a").iterdir()) == [f"{name}.h5ad", f"{name}.h5ad.stamp.json"]


@pytest.mark.parametrize("ctl,step,code", [
    ({"tests_exit": 1}, "tests", 10),
    ({"fail": "*/gate/c2-s2.h5ad 1"}, "c2", 11),
    ({"fail": "*/gate/c2-s2.h5ad 2"}, "c2", 12),
    ({"fail": "*/gate/T3_c2_s2_g2.h5ad 6"}, "c2g2", 16),  # line 3's own 5% guard
    ({"fail": "*/gate/T3_c2_s2_s067.h5ad 5"}, "s067", 17),
    ({"fail": "*/gate/T3_c2_s2_gonly.h5ad 2"}, "gonly", 12),
    ({"check_exit": 1}, "c2", 13),
    ({"unclean": "*/gate/T3_c2_s2_s067.h5ad"}, "s067", 15),
    ({"fail": "*/a/T3_c2_s2_g15.h5ad 6"}, "a", 16),  # the 5% den check stops the file: nothing written
    ({"fail": "*/b/T3_c2_s2_g15.h5ad 5"}, "b", 17),
    ({"fail": "*/a/T3_c2_s2_g15.h5ad 2"}, "a", 12),
    ({"fail": "*/b/T3_c2_s2_g15.h5ad 1"}, "b", 11),
    ({"unclean": "*/b/T3_c2_s2_g15.h5ad"}, "b", 15),
    ({"differ": "*/b/T3_c2_s2_g15.h5ad"}, "done", 41),
])
def test_each_stop_after_the_gate(tmp_path, ctl, step, code):
    repo, env, head = _repo(tmp_path, ctl=ctl)
    res = _run(repo, env, head, "r1", "T3_c2_s2_g15")
    assert res.returncode == code, res.stderr
    s = _scratch(tmp_path)
    ex = _exits(s)
    assert ex.pop("all") == str(code) and ex.pop(step) == str(code)
    all_steps = steps()
    assert ex == dict.fromkeys(all_steps[:all_steps.index(step)], "0")  # every earlier step passed, no later began
    assert (s / "SHA256SUMS").exists() is (step == "done")
    if code == 41:
        assert "the two builds of T3_c2_s2_g15.h5ad differ" in res.stderr


@pytest.mark.parametrize("step,name,var,code,args", IDENTITY, ids=[i[0] for i in IDENTITY])
def test_an_identity_rebuild_of_other_bytes_stops_the_run_at_its_step(tmp_path, step, name, var, code, args):
    """The stand-in's hashes for every pin but this rebuild's, which keeps its recorded sha256."""
    repo, env, head = _repo(tmp_path, sh=_stand_in(_sh(), keep=(var,)))
    res = _run(repo, env, head, "r1", "T3_c2_s2_g2_s067")
    assert res.returncode == code and f"gate/{name} is not its recorded bytes" in res.stderr, res.stderr
    ex = _exits(_scratch(tmp_path))
    assert ex.pop(step) == str(code) and ex.pop("all") == str(code)
    all_steps = steps()
    assert ex == dict.fromkeys(all_steps[:all_steps.index(step)], "0")
    assert not (repo / "runs" / f"t3r76_{head[:7]}_r1" / "a").exists()


def _stopped_at_the_gate(tmp_path: Path, code: int, run: str = "r1") -> None:
    s = _scratch(tmp_path, run)
    assert _exits(s) == {"gate": str(code), "all": str(code)}
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "tmp"]
    assert not list((s / "tmp").iterdir()) and not _calls(tmp_path)


ROW_REFUSALS = [
    ("", FIRST, "exactly one row R76"),
    (ROW + ROW, FIRST, "exactly one row R76"),
    (ROW.replace("| R76 (", "| R76b ("), FIRST, "exactly one row R76"),
    (r76_row(opener="R76 (a proposal of a day"), FIRST, "not the captain's ruling"),
    (r76_row(opener="R76 awaits (the captain's order"), FIRST, "not the captain's ruling"),
    (r76_row(ruling="a try"), FIRST, "does not say ruling, no run"),
    (r76_row(captain="--captain-ruling R76"), "T3_c2_s2_g15", "does not name `--captain-ruling R76`"),
    (r76_row(captain="`--captain-ruling R76 now`"), FIRST, "does not name `--captain-ruling R76`"),
    (r76_row(den="--den-ruling R76"), FIRST, "does not name `--den-ruling R76`"),
    (r76_row(den="`--den-ruling R47`"), "T3_c2_s2_g2_s067", "does not name `--den-ruling R76`"),
    (r76_row(pairs=", ".join(f"`{OPT[n]} --x` ({n})" for n in NAMES), lifted=()), FIRST,
     f"does not name `{OPT[FIRST]}`"),
    (r76_row(pairs=", ".join(f"{OPT[n]} ({n})" for n in NAMES), lifted=()), FIRST, f"does not name `{OPT[FIRST]}`"),
    (r76_row(pairs=", ".join(f"`a`{OPT[n]}`b` ({n})" for n in NAMES), lifted=()), "T3_c2_s2_g05",
     f"does not name `{OPT['T3_c2_s2_g05']}`"),
    (r76_row(pairs=", ".join(f"`{OPT[n]}` ({m})" for n, m in zip(NAMES, [*NAMES[1:], NAMES[0]], strict=True))),
     FIRST, f"does not pair `{OPT[FIRST]}` with ({FIRST})"),
    (r76_row(pairs=", ".join(f"`{OPT[n]}`" for n in NAMES)), "T3_c2_s2_g15", "does not pair"),
    (r76_row(pairs=", ".join(f"`{OPT[n]}` ({n}x)" for n in NAMES)), "T3_c2_s2_g05", "does not pair"),
    ("| R76 (the captain's order of a day, ruling, no run; not a try) | T3 | "
     + ", ".join(f"`{OPT[n]}` ({n})" for n in NAMES) + " `--captain-ruling R76` `--den-ruling R76` |\n", FIRST,
     f"does not name `{OPT[FIRST]}`"),  # named outside the first cell
]


@pytest.mark.parametrize("tries,name,says", ROW_REFUSALS, ids=[f"r76-{i}" for i in range(len(ROW_REFUSALS))])
def test_each_row_r76_refusal(tmp_path, tries, name, says):
    repo, env, head = _repo(tmp_path, tries=tries)
    res = _run(repo, env, head, "r1", name)
    assert res.returncode == 95 and says in res.stderr, res.stderr
    _stopped_at_the_gate(tmp_path, 95)
    assert not (repo / "runs").exists()


def test_a_row_without_the_lift_opens_the_weight_only_files_only(tmp_path):
    """The gate asks `--den-ruling R76` of the three lifted files only (the builder asks the row's lifting part of
    every file, so such a row stops each build there: its exit 2, here the stand-in's)."""
    rows = r76_row(den="(no lift)", lifted=())
    repo, env, head = _repo(tmp_path, tries=rows, ctl={"tests_exit": 1})
    for i, name in enumerate(NAMES):
        res = _run(repo, env, head, f"r{i}", name)
        if name in WEIGHT_ONLY:
            assert res.returncode == 10, (name, res.stderr)
        else:
            assert res.returncode == 95 and "does not name `--den-ruling R76`" in res.stderr, (name, res.stderr)


@pytest.mark.parametrize("kw,code,says", [
    ({"data": False}, 92, "the data data/raw/T2/E8.75.h5ad or the decision file"),
    ({"dec": False}, 92, "the data data/raw/T2/E8.75.h5ad or the decision file"),
    ({"run_folder": None}, 87, "does not hold the two tables at their pins"),
    ({"run_folder": "changed"}, 87, "does not hold the two tables at their pins"),
    ({"sh": "as is"}, 87, "does not hold the two tables at their pins"),  # the real pins against stand-in tables
    ({"venv": None}, 88, "no .venv of its own"),
    ({"venv": "other"}, 88, "no .venv of its own"),  # a venv whose vec.pth names another checkout's src
])
def test_each_later_gate_refusal(tmp_path, kw, code, says):
    kw = dict(kw)
    if kw.get("sh") == "as is":
        kw["sh"] = _sh()
    repo, env, head = _repo(tmp_path, **kw)
    res = _run(repo, env, head, "r1", FIRST)
    assert res.returncode == code and says in res.stderr, res.stderr
    _stopped_at_the_gate(tmp_path, code)
    assert not (repo / "runs").exists()


def test_the_gate_passes_on_row_r76_as_the_repository_s_ledger_holds_it(tmp_path):
    """Row R76's first cell read by code from the repository's ledger (never printed), set in a throwaway ledger: the
    gate passes for each of the five names and the run goes on to the tests step."""
    lines = [ln for ln in REPO_LEDGER.read_text(encoding="utf-8").splitlines() if ln.startswith("| R76 ")]
    assert len(lines) == 1, "the repository's ledger holds no single row R76"
    first = lines[0][2:].split("|", 1)[0]
    repo, env, head = _repo(tmp_path, tries=f"| {first}| T3 | x |\n", ctl={"tests_exit": 1})
    for i, name in enumerate(NAMES):
        res = _run(repo, env, head, f"r{i}", name)
        assert res.returncode == 10, (name, res.stderr)
        assert _exits(_scratch(tmp_path, f"r{i}")) == {"gate": "0", "tests": "10", "all": "10"}


@pytest.mark.parametrize("argv,code", [
    ([], None), (["r1"], None), (["r1", ""], None),  # missing arguments (an empty name is none)
    (["a/b", FIRST], 90), (["r 1", FIRST], 90), (["$(id)", FIRST], 90), (["r.1", FIRST], 90),
    (["r1", FIRST, "extra"], 89),
    *[(["r1", bad], 89) for bad in (
        "T3_c2_s2_s11", "T3_c2_s2_s067", "T3_c2_s2_g2", "T3_c2_s2_gonly", "T3_c2_s2_gonly_s067", "T3_c2_s2",
        "T3_c2_s2_s085.h5ad", "t3_c2_s2_s085", "T3_C2_S2_S085", "T3_c2_s2_s085 ", " T3_c2_s2_s085", "T3_c2_s2_s08",
        "T3_c2_s2_s5", "T3_c2_s2_s050", "T3_c2_s2_g1.5", "T3_c2_s2_g150", "T3_c2_s2_g2_s067x", "T3_c2_s2_g20_s067",
        "T3_c2_s2_s067_g2", "T3_c2_s2_g15_s085", "T3_c2_s2_s035", "T3_c2_s2_dd", "T3_c2_s2_cd7", "T3_c2_s2_mxe2",
        "s085", "c2-s2", "*", "T3_c2_s2_*", "T3_c2_s2_s085*", "$(id)")],
])
def test_odd_arguments_and_unknown_names_are_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _repo(tmp_path)
    res = _run(repo, env, head, *argv)
    assert (res.returncode != 0 and USAGE in res.stderr) if code is None else res.returncode == code, res.stderr
    assert not (tmp_path / "home" / "v" / "scratch").exists() and not (repo / "runs").exists()
    assert not _calls(tmp_path)


def test_a_wrong_commit_a_prefix_and_an_unclean_tree_are_refused_and_an_untracked_file_is_not(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    res = _run(repo, env, "0" * 40, "r1", FIRST)
    assert res.returncode == 96 and "HEAD is not" in res.stderr
    _stopped_at_the_gate(tmp_path, 96)
    assert _run(repo, env, head[:12], "r2", FIRST).returncode == 96
    _stopped_at_the_gate(tmp_path, 96, "r2")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW}more\n")
    res = _run(repo, env, head, "r3", FIRST)
    assert res.returncode == 98 and "tracked files differ" in res.stderr
    _stopped_at_the_gate(tmp_path, 98, "r3")
    _git(repo, env, "checkout", "-q", "--", "experiments/harness/TRIES.md")
    (repo / "stray.py").write_text("x = 1\n")  # untracked: the stamps' rule is the builder's own check
    assert _run(repo, env, head, "r4", FIRST).returncode == 10
    assert not (repo / "runs" / f"t3r76_{head[:7]}_r1").exists()


def test_a_run_name_is_never_reused_and_an_existing_output_folder_is_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1", FIRST).returncode == 10
    s = _scratch(tmp_path)
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", "T3_c2_s2_g05")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    o = repo / "runs" / f"t3r76_{head[:7]}_r2"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    res = _run(repo, env, head, "r2", FIRST)
    assert res.returncode == 93 and "output folder" in res.stderr
    assert _exits(_scratch(tmp_path, "r2")) == {"gate": "93", "all": "93"}
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


def test_a_call_that_passes_the_gate_writes_the_head_and_its_folder(tmp_path):
    repo, env, head = _repo(tmp_path, ctl={"tests_exit": 1})
    assert _run(repo, env, head, "r1", "T3_c2_s2_g2_s067").returncode == 10
    s = _scratch(tmp_path)
    assert _exits(s) == {"gate": "0", "tests": "10", "all": "10"} and (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["all.exit", "gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t3r76_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())
