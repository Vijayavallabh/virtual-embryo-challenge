"""experiments/harness/run_h35.sh, row H35's launcher (G0, the gate, the run and the read of
experiments/harness/h35_heart_trend_forms.py). Every call runs in a throwaway git repository under the test's own
folder, with a temporary home and a stand-in uv first on PATH: no step of the script runs, nothing reads a competition
file and nothing is built or scored. The launcher's row check is held to the script's own on the same ledgers, and
each call it makes of the script to the script's own command line."""

from __future__ import annotations

import importlib.util
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = REPO_ROOT / "experiments" / "harness" / "run_h35.sh"
SCRIPT = "experiments/harness/h35_heart_trend_forms.py"
OUT_DIR = "experiments/harness/results/h35_heart_trend_forms"
KEY1, KEY2 = "T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"
FIRST = "H35 (the counted screen of row R55's trend forms: ten tries; G0 first, target-free, not a try)"
STEPS = ("check", "tests", "g0", "gate", "run", "read")
ONE_CALL = ["tests", "ask g0", "g0", "ask gate", "gate", "ask run", *["run"] * 10, "ask read", "read"]


def h35_row(cell="pending", first=FIRST, cells=11):
    return f"| {first} | heart | " + " | ".join(["x"] * (cells - 3)) + f" | {cell} |\n"


# A stand-in for uv: it records each call with the environment the launcher gave it and, at that moment, the step's
# exit file in the call's folder, the call's done.exit and the step's once-only exit file; it imitates each call.
STUB_UV = r"""#!/usr/bin/env bash
BASE="$HOME/v/scratch/h35"
S="$BASE/$STUB_RUN"
threads="$OMP_NUM_THREADS|$OPENBLAS_NUM_THREADS|$MKL_NUM_THREADS"
echo "$UV_PROJECT_ENVIRONMENT|$threads|$NUMPY_MADVISE_HUGEPAGE|$TMPDIR|$*" >> "$HOME/uvcalls"
quoted=$(printf '%q ' "$@")
echo "$quoted" >> "$HOME/uvargs"
args=("$@")
[ "${args[0]}" = run ] || exit 90
args=("${args[@]:1}")
while [ "${args[0]}" = --frozen ] || [ "${args[0]}" = --no-sync ]; do args=("${args[@]:1}"); done
[ "${args[0]}" = python ] || exit 91
args=("${args[@]:1}")
while [ "${args[0]}" = -I ] || [ "${args[0]}" = -B ]; do args=("${args[@]:1}"); done
out=""; seeds=""; step=""; prev=""
for a in "${args[@]}"; do
  [ "$prev" = "--out" ] && out="$a"
  [ "$prev" = "--seeds" ] && seeds="$a"
  [ "$prev" = "--step" ] && step="$a"
  prev="$a"
done
state() {
  once=$(cat "$BASE/$2.exit" 2>/dev/null || echo none)
  echo "$1 $2 $(cat "$S/$2.exit" 2>/dev/null) $(cat "$S/done.exit" 2>/dev/null) $once" >> "$HOME/states"
}
case "${args[0]}" in
  -m) state call tests; echo "some test output"; exit "${STUB_TESTS_EXIT:-0}" ;;
  experiments/harness/h35_heart_trend_forms.py)
    case "${args[1]}" in
      refusals)
        state ask "$step"
        if [ "${STUB_REFUSE:-}" = "$step" ]; then echo "H35 is refused: a reason"; exit 1; fi
        if [ "${STUB_RACE:-}" = "$step" ]; then echo running > "$BASE/$step.exit"; fi
        exit 0 ;;
      g0)
        state call g0
        if [ "${STUB_G0_EXIT:-0}" = kill ]; then kill -9 "$PPID"; exit 137; fi
        if [ "${STUB_G0_EXIT:-0}" = 1 ]; then echo "Traceback"; exit 1; fi
        mkdir -p "$out"; echo '{"stub": "g0"}' > "$out/g0.json"
        echo "G0 category PASS"; echo "G0 contrast IDLE"; echo "G0 whole DROP"; echo "G0 portion PASS"
        echo "a library warning"
        case "${STUB_G0_EXIT:-0}" in
          3) echo "G0 STOP"; echo "OUTCOME H35: STOP: category IDLE; contrast IDLE; whole DROP; portion IDLE."
             echo "A screen, not evidence of development. What follows: nothing is built." ;;
          5) echo "G0 VOID"; echo "OUTCOME H35: VOID: G0: x; nothing is built."
             echo "A screen, not evidence of development. What follows: nothing is built." ;;
          *) echo "G0 PASS" ;;
        esac
        exit "${STUB_G0_EXIT:-0}" ;;
      gate)
        state call gate
        if [ "${STUB_GATE_EXIT:-0}" = kill ]; then kill -9 "$PPID"; exit 137; fi
        if [ "${STUB_GATE_EXIT:-0}" = 1 ]; then echo "Traceback"; exit 1; fi
        mkdir -p "$out"; echo '{"stub": "gate"}' > "$out/gate.json"
        echo "T2h/<=E9.5->E10.5 seed 0 (1.0 s) done"; echo "a library warning"
        if [ "${STUB_GATE_EXIT:-0}" = 3 ]; then
          echo "GATE FAIL"; echo "OUTCOME H35: VOID: v3: the gate fails (G1 (M): x); nothing is built."
          echo "A screen, not evidence of development. What follows: nothing is built."
        else
          echo "GATE PASS"
        fi
        exit "${STUB_GATE_EXIT:-0}" ;;
      run)
        state call run
        if [ "${STUB_RUN_BAD:-}" = "$seeds" ]; then echo "boom"; exit 7; fi
        mkdir -p "$out"
        [ "${STUB_RUN_NOFILE:-}" = "$seeds" ] || echo '{"stub": "run"}' > "$out/h35.json"
        echo "T2h/x seed 0 (1.0 s) done"; echo "a library warning"; echo "wrote $out/h35.json"; exit 0 ;;
      read)
        state call read
        if [ "${STUB_READ_EXIT:-0}" = kill ]; then kill -9 "$PPID"; exit 137; fi
        [ "${STUB_READ_EXIT:-0}" = 0 ] || exit "$STUB_READ_EXIT"
        mkdir -p "$out"; echo "the reading" > "$out/reading.txt"; echo "{}" > "$out/reading.json"
        echo "OUTCOME H35: category KEPT; contrast IDLE; whole DROP; portion IDLE."
        echo "A screen, not evidence of development. What follows: x"; exit 0 ;;
    esac ;;
esac
exit 92
"""


def _env(tmp_path, extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GIT_", "H35_", "STUB_"))}
    env.update(HOME=str(tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path),
               PATH=f"{tmp_path / 'bin'}{os.pathsep}{os.environ['PATH']}")
    env.update(extra or {})
    return env


def _stand_in_env(path):
    (path / "bin").mkdir(parents=True)
    (path / "bin" / "python").write_text("#!/bin/sh\nexit 0\n")
    (path / "bin" / "python").chmod(0o755)
    return path


@pytest.fixture()
def repo(tmp_path):
    """A throwaway git repository with a ledger that holds row H35 in its fixed form and no outcome, a stand-in for the
    script, for the virtual environment and for the data folder, a temporary home and the stand-in uv."""
    root = tmp_path / "repo"
    (root / "experiments" / "harness").mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "uv").write_text(STUB_UV)
    (tmp_path / "bin" / "uv").chmod(0o755)
    (root / "experiments" / "harness" / "TRIES.md").write_text("| H34 (a row) | heart | x |\n" + h35_row())
    (root / SCRIPT).write_text("# a stand-in\n")
    (root / ".gitignore").write_text("runs/\n.venv/\ndata/\n")
    _stand_in_env(root / ".venv")
    (root / "data" / "raw" / "T2").mkdir(parents=True)

    def git(*a):
        return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "init.defaultBranch=x", *a],
                              cwd=root, env=_env(tmp_path), capture_output=True, text=True, check=True).stdout.strip()
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "one")
    return SimpleNamespace(root=root, git=git, head=git("rev-parse", "HEAD"), tmp=tmp_path,
                           scratch=tmp_path / "home" / "v" / "scratch" / "h35")


def _launch(repo, *args, cwd=None, **extra):
    run = args[1] if len(args) > 1 else ""
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=cwd or repo.root,
                       env=_env(repo.tmp, {"STUB_RUN": run, **extra}), capture_output=True, text=True, timeout=120)
    return r.returncode, r.stderr


def _uv_calls(repo):
    f = repo.tmp / "home" / "uvcalls"
    return f.read_text().splitlines() if f.exists() else []


def _uv_args(repo):
    f = repo.tmp / "home" / "uvargs"
    return [shlex.split(ln) for ln in f.read_text().splitlines()] if f.exists() else []


def _kinds(repo):
    out = []
    for c in _uv_calls(repo):
        words = c.split("python -I -B ", 1)[1].split()
        out.append("tests" if words[0] == "-m" else f"ask {words[3]}" if words[1] == "refusals" else words[1])
    return out


def _exits(repo, run):
    return {p.name.removesuffix(".exit"): p.read_text().strip() for p in sorted((repo.scratch / run).glob("*.exit"))}


def _once(repo):
    return {p.name.removesuffix(".exit"): p.read_text().strip() for p in sorted(repo.scratch.glob("*.exit"))}


def _next_commit(repo):
    (repo.root / "note.txt").write_text("x\n")
    repo.git("add", "note.txt")
    repo.git("commit", "-q", "-m", "two")
    return repo.git("rev-parse", "HEAD")


@pytest.fixture(scope="module")
def h35():
    spec = importlib.util.spec_from_file_location("h35_for_run_h35_launcher", REPO_ROOT / SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_the_launcher_text():
    sh = LAUNCHER.read_text()
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    assert os.access(LAUNCHER, os.X_OK)
    code = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert not re.search(r"\brm\b|\brmdir\b|\bunlink\b|\bmv\b|-delete\b|\bshred\b", sh)  # it deletes nothing
    assert not re.search(r"\bscp\b|\brsync\b|\bcurl\b|\bwget\b|\bssh\b|vec submit|\bupload", sh)  # and sends nothing
    for m in re.finditer(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in m.group(1) and '"' not in m.group(1), m.group(0)
    assert len(re.findall(r"\$\{[A-Za-z0-9_]+:\?", sh)) == 2 and "set -o pipefail" in code
    assert code.index("export PATH=$HOME/.local/bin:$PATH") < code.index("uv run")
    assert 'PY="uv run --frozen --no-sync python -I -B"' in code
    assert "VENV=${H35_VENV:-$PWD/.venv}\n" in code and 'export UV_PROJECT_ENVIRONMENT="$VENV"\n' in code
    assert "DATA=${H35_DATA:-data/raw}\n" in code
    assert max(len(ln) for ln in sh.splitlines()) <= 120
    assert "export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0 " \
           'TMPDIR="$S/tmp"' in code
    at = {s: code.index(f"STEP={s}\n") for s in STEPS}
    assert list(at.values()) == sorted(at.values()) and code.count("STEP=") == 6
    assert f"OUT={OUT_DIR}\n" in code
    assert all(">&2" in ln for ln in code.splitlines() if re.search(r"\becho \"(refused|usage)", ln))
    for test in ("tests/test_h35_heart_trend_forms.py", "tests/test_run_h35_launcher.py",
                 "tests/test_heart_trend_line_h.py", "tests/test_h34_heart_organ_trend.py",
                 "tests/test_h29b_heart_organ.py", "tests/test_h27_heart_forecast_size.py"):
        assert test in code, test
    # G0, the gate, the run and the read are each started through start_once, before the step's own call
    assert re.findall(r"^start_once (.+)$", code, flags=re.M) == [
        '--out "$OUT/g0"', '--out "$OUT/gate" --g0 "$OUT/g0"', '--g0 "$OUT/g0" --gate "$OUT/gate"', '--out "$OUT"']
    assert code.count("set -o noclobber; echo running > \"$BASE/$STEP.exit\"") == 1
    for step in ("g0", "gate", "run", "read"):
        assert at[step] < code.index("\nstart_once --", at[step]) < code.index(f'$PY "$SCRIPT" {step} ', at[step])
    body = code[code.index("start_once() {"):code.index("count_other() {")]
    order = [body.index(s) for s in ('if [ -e "$BASE/$STEP.exit" ]; then', "finish 88", "finish 92",
                                     '"$SCRIPT" refusals --step ', "finish 89", "set -o noclobber", "OWNED=$STEP",
                                     '> "$BASE/$STEP.head"')]
    assert order == sorted(order) and "cat " not in body and "return" not in body


def _check_code() -> str:
    sh = LAUNCHER.read_text()
    i = sh.index("python3 -I -B -c '\nimport sys") + len("python3 -I -B -c '")
    return sh[i:sh.index("' experiments/harness/TRIES.md\ncode=$?", i)]


def test_the_launchers_row_check_is_the_scripts(h35, tmp_path):
    """The check step's inline reader of row H35 passes exactly the ledgers the script takes (row_formed and
    row_open): 0 open, 4 formed with an outcome, 3 missing or malformed."""
    led = "| H34 (a row) | heart | x |\n"
    cases = [led + h35_row(), led + h35_row("OUTCOME H35: VOID: x; nothing is built."), led, led + h35_row() * 2,
             led + h35_row(cells=10), led + h35_row(cells=12), led + h35_row(first=FIRST.replace("ten", "eight")),
             led + h35_row(first="**" + FIRST + "**"), led + h35_row().replace("| H35 (", "| H350 ("),
             led + h35_row("OUTCOME H35: category KEPT; contrast IDLE; whole DROP; portion IDLE."),
             led + h35_row("a ||pipe|| in a cell"), led + h35_row().rstrip(" |\n") + "\n"]
    seen = set()
    for i, text in enumerate(cases):
        p = tmp_path / f"t{i}.md"
        p.write_text(text)
        got = subprocess.run(["python3", "-I", "-B", "-c", _check_code(), str(p)]).returncode
        want = 0 if h35.row_open(text) else 4 if h35.row_formed(text) else 3
        assert got == want, (i, got, want)
        seen.add(got)
    assert seen == {0, 3, 4}
    assert h35.FIRST_CELL == FIRST


def test_the_launcher_runs_the_five_steps_in_order(repo):
    code, err = _launch(repo, repo.head, "r1")
    assert (code, err) == (0, "")
    S = repo.scratch / "r1"
    assert _exits(repo, "r1") == {s: "0" for s in (*STEPS, "done")}
    assert (S / "head").read_text().strip() == repo.head
    assert _once(repo) == {"g0": "0", "gate": "0", "run": "0", "read": "0"}
    assert all((repo.scratch / f"{s}.head").read_text() == f"{repo.head} r1\n" for s in ("g0", "gate", "run", "read"))
    states = (repo.tmp / "home" / "states").read_text().splitlines()
    assert states == ["call tests running running none", "ask g0 running running none",
                      "call g0 running running running", "ask gate running running none",
                      "call gate running running running", "ask run running running none",
                      *["call run running running running"] * 10, "ask read running running none",
                      "call read running running running"]
    calls = _uv_calls(repo)
    env = f"{repo.root}/.venv|8|8|8|0|{S}/tmp|run --frozen --no-sync python -I -B "
    assert len(calls) == 18 and all(c.startswith(env) for c in calls) and _kinds(repo) == ONE_CALL
    g0, gate = f"{OUT_DIR}/g0", f"{OUT_DIR}/gate"
    assert calls[0].startswith(env + f"-m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp={S}/tmp/pytest ")
    assert calls[1] == env + f"{SCRIPT} refusals --step g0 --commit {repo.head} --out {g0}"
    assert calls[2] == env + f"{SCRIPT} g0 --out {g0} --commit {repo.head}"
    assert calls[3] == env + f"{SCRIPT} refusals --step gate --commit {repo.head} --out {gate} --g0 {g0}"
    assert calls[4] == env + f"{SCRIPT} gate --out {gate} --g0 {g0} --commit {repo.head}"
    assert calls[5] == env + f"{SCRIPT} refusals --step run --commit {repo.head} --g0 {g0} --gate {gate}"
    parts = ("0-1", "2-3", "4-5", "6-7", "8-9")
    for block, (key, short) in zip((calls[6:11], calls[11:16]), ((KEY1, "e95"), (KEY2, "e875")), strict=True):
        assert sorted(block) == sorted(env + f"{SCRIPT} run --task {key} --seeds {p} --g0 {g0} --gate {gate} "
                                       f"--commit {repo.head} --out {OUT_DIR}/parts/{short}/{p}" for p in parts)
    assert calls[16] == env + f"{SCRIPT} refusals --step read --commit {repo.head} --out {OUT_DIR}"
    assert calls[17] == (env + f"{SCRIPT} read --runs {OUT_DIR}/parts --g0 {g0} --gate {gate} --commit {repo.head} "
                         f"--out {OUT_DIR}")
    assert len(sorted((repo.root / OUT_DIR / "parts").rglob("h35.json"))) == 10
    for short in ("e95", "e875"):
        for p in parts:
            assert (S / "exit" / f"{short}-{p}").read_text().strip() == "0"
            assert (S / "exit" / f"{short}-{p}.other_lines").read_text().strip() == "1"
    assert (S / "exit" / "other_lines").read_text().strip() == "10"
    assert (S / "g0.other_lines").read_text().strip() == "1" and (S / "gate.other_lines").read_text().strip() == "1"
    assert (S / "read.out").read_text().splitlines() == [
        "OUTCOME H35: category KEPT; contrast IDLE; whole DROP; portion IDLE.",
        "A screen, not evidence of development. What follows: x"]
    assert repo.git("status", "--porcelain", "--untracked-files=no") == ""
    # the run folder is never reused, and a later call at the same commit finds the reading: nothing is touched
    code, err = _launch(repo, repo.head, "r1")
    assert code == 91 and "the run folder" in err
    code, err = _launch(repo, repo.head, "r2")
    assert code == 95 and "already holds a reading" in err and len(_uv_calls(repo)) == 18
    # with the reading out of the way a later call gets as far as G0, which was started once already
    (repo.root / OUT_DIR / "reading.txt").rename(repo.tmp / "aside.txt")
    (repo.root / OUT_DIR / "reading.json").rename(repo.tmp / "aside.json")
    for i in range(2):
        code, err = _launch(repo, _next_commit(repo) if i else repo.head, f"r{3 + i}")
        assert code == 88 and "the g0 step was started before" in err
        assert _exits(repo, f"r{3 + i}") == {"check": "0", "tests": "0", "g0": "88", "done": "88"}
    assert _once(repo) == {"g0": "0", "gate": "0", "run": "0", "read": "0"}


def test_the_launchers_calls_are_calls_the_script_takes(repo, h35, monkeypatch):
    """Every call the launcher makes of the script, word for word, parses under the script's own command line and
    reaches the step it names; the gate is given no reading and no parts of its own choosing (the committed ones)."""
    data = repo.tmp / "data dir"  # a space in the path: each argument stays one argument
    (data / "T2").mkdir(parents=True)
    assert _launch(repo, repo.head, "x1", H35_DATA=str(data)) == (0, "")
    reached = []
    for name in ("g0", "gate", "run", "read", "refusals"):
        monkeypatch.setattr(h35, f"cmd_{name}", lambda a, name=name: reached.append((name, vars(a))))
    calls = _uv_args(repo)
    assert len(calls) == len(ONE_CALL)
    for args in calls:
        assert args[:6] == ["run", "--frozen", "--no-sync", "python", "-I", "-B"]
        if SCRIPT in args:
            h35.main(args[args.index(SCRIPT) + 1:])
    assert [name for name, _ in reached] == [k.split()[0].replace("ask", "refusals") for k in ONE_CALL[1:]]
    by: dict = {}
    for name, a in reached:
        by.setdefault(name, []).append(a)
    g0, gate = f"{OUT_DIR}/g0", f"{OUT_DIR}/gate"
    assert by["g0"] == [{"cmd": "g0", "out": g0, "commit": repo.head, "data": str(data)}]
    assert by["gate"] == [{"cmd": "gate", "out": gate, "g0": g0, "commit": repo.head, "data": str(data),
                           "h29b": str(h35.h34.H29B_READING), "h27": str(h35.h34.H27_PARTS)}]
    parts = ("0-1", "2-3", "4-5", "6-7", "8-9")
    want = {(k, p, f"{OUT_DIR}/parts/{short}/{p}") for k, short in ((KEY1, "e95"), (KEY2, "e875")) for p in parts}
    assert {(a["task"], a["seeds"], a["out"]) for a in by["run"]} == want and len(by["run"]) == 10
    assert all(a["g0"] == g0 and a["gate"] == gate and a["commit"] == repo.head and a["data"] == str(data)
               for a in by["run"])
    assert by["read"] == [{"cmd": "read", "runs": f"{OUT_DIR}/parts", "g0": g0, "gate": gate, "out": OUT_DIR,
                           "commit": repo.head}]
    assert [(a["step"], a["out"], a["g0"], a["gate"]) for a in by["refusals"]] == [
        ("g0", g0, None, None), ("gate", gate, g0, None), ("run", None, g0, gate), ("read", OUT_DIR, None, None)]
    assert h35.RESULTS == h35.HERE / OUT_DIR and h35.PART == "h35.json"


@pytest.mark.parametrize("extra, step, code, n_calls, once", [
    ({"STUB_TESTS_EXIT": "1"}, "tests", 1, 1, {}),
    ({"STUB_REFUSE": "g0"}, "g0", 89, 2, {}),
    ({"STUB_G0_EXIT": "3"}, "g0", 3, 3, {"g0": "3"}),
    ({"STUB_G0_EXIT": "5"}, "g0", 5, 3, {"g0": "5"}),
    ({"STUB_G0_EXIT": "1"}, "g0", 1, 3, {"g0": "1"}),
    ({"STUB_REFUSE": "gate"}, "gate", 89, 4, {"g0": "0"}),
    ({"STUB_GATE_EXIT": "3"}, "gate", 3, 5, {"g0": "0", "gate": "3"}),
    ({"STUB_REFUSE": "run"}, "run", 89, 6, {"g0": "0", "gate": "0"}),
    ({"STUB_RUN_BAD": "4-5"}, "run", 1, 11, {"g0": "0", "gate": "0", "run": "1"}),
    ({"STUB_RUN_NOFILE": "8-9"}, "run", 1, 11, {"g0": "0", "gate": "0", "run": "1"}),
    ({"STUB_REFUSE": "read"}, "read", 89, 17, {"g0": "0", "gate": "0", "run": "0"}),
    ({"STUB_READ_EXIT": "1"}, "read", 1, 18, {"g0": "0", "gate": "0", "run": "0", "read": "1"}),
])
def test_the_launcher_stops_at_the_first_step_that_does_not_exit_0(repo, extra, step, code, n_calls, once):
    got, err = _launch(repo, repo.head, "s1", **extra)
    assert got == code
    want = {s: "0" for s in STEPS[:STEPS.index(step)]} | {step: str(code), "done": str(code)}
    assert _exits(repo, "s1") == want
    assert _kinds(repo) == ONE_CALL[:n_calls]
    assert _once(repo) == once
    S = repo.scratch / "s1"
    if code == 89:
        assert f"the script refuses the {step} step" in err and step not in _kinds(repo)
        assert not (repo.scratch / f"{step}.head").exists()
    if step == "g0" and code in (3, 5):  # G0's own outcome lines are kept and counted as expected ones
        lines = (S / "g0.out").read_text().splitlines()
        assert lines[-3] == ("G0 STOP" if code == 3 else "G0 VOID") and lines[-2].startswith("OUTCOME H35: ")
        assert (S / "g0.other_lines").read_text().strip() == "1"
    if step != "read" or code != 0:
        assert not (repo.root / OUT_DIR / "reading.txt").exists()


@pytest.mark.parametrize("extra, first, once", [
    ({"STUB_G0_EXIT": "3"}, 3, {"g0": "3"}),  # G0 said STOP
    ({"STUB_G0_EXIT": "kill"}, -9, {"g0": "running"}),  # the call was killed in G0: no verdict
    ({"STUB_GATE_EXIT": "3"}, 3, {"g0": "0", "gate": "3"}),
    ({"STUB_RUN_BAD": "4-5"}, 1, {"g0": "0", "gate": "0", "run": "1"}),
    ({"STUB_READ_EXIT": "kill"}, -9, {"g0": "0", "gate": "0", "run": "0", "read": "running"}),
])
def test_no_step_is_started_a_second_time(repo, extra, first, once):
    assert _launch(repo, repo.head, "a1", **extra)[0] == first
    assert _once(repo) == once
    n = len(_uv_calls(repo))
    for i in range(2):  # a later call at the same commit, then at another
        got, err = _launch(repo, _next_commit(repo) if i else repo.head, f"a{2 + i}")
        assert got == 88 and "the g0 step was started before" in err
        assert _exits(repo, f"a{2 + i}") == {"check": "0", "tests": "0", "g0": "88", "done": "88"}
        assert _kinds(repo)[n:] == ["tests"] * (1 + i)
    assert _once(repo) == once


@pytest.mark.parametrize("step", ["g0", "gate", "run", "read"])
@pytest.mark.parametrize("holds", ["running", "0", ""])
def test_each_of_the_four_steps_has_its_own_once_only_file(repo, step, holds):
    repo.scratch.mkdir(parents=True)
    (repo.scratch / f"{step}.exit").write_text(holds + "\n")
    got, err = _launch(repo, repo.head, "o1")
    assert got == 88 and f"the {step} step was started before" in err
    assert _exits(repo, "o1") == {s: "0" for s in STEPS[:STEPS.index(step)]} | {step: "88", "done": "88"}
    assert _kinds(repo) == ONE_CALL[:ONE_CALL.index(f"ask {step}")]
    assert (repo.scratch / f"{step}.exit").read_text() == holds + "\n"


def test_parts_made_outside_the_launcher_and_a_race_start_nothing(repo):
    parts = repo.root / OUT_DIR / "parts"
    parts.mkdir(parents=True)
    got, err = _launch(repo, repo.head, "p1")
    assert got == 92 and "made outside this launcher" in err
    assert _kinds(repo) == ["tests", "ask g0", "g0", "ask gate", "gate"] and _once(repo) == {"g0": "0", "gate": "0"}
    got, err = _launch(repo, repo.head, "c1", STUB_RACE="g0")
    assert got == 88 and _exits(repo, "c1")["g0"] == "88"


def test_the_launchers_own_refusals(repo):
    def refused(run, code, words, *args, **kw):
        got, err = _launch(repo, *(args or (repo.head, run)), **kw)
        assert got == code and words in err, (run, got, err)
        assert _exits(repo, run) == {"check": str(code), "done": str(code)}
    for args in ((), (repo.head,), (repo.head, ""), ("", "x")):
        got, err = _launch(repo, *args)
        assert got == 1 and "usage: run_h35.sh COMMIT RUN" in err
    for odd in ("a/b", "..", "a b", "r.1", "x;y"):
        got, err = _launch(repo, repo.head, odd)
        assert got == 90 and "letters, digits" in err
    got, err = _launch(repo, repo.head, "x", "y")
    assert got == 90 and "COMMIT and RUN only" in err
    assert not repo.scratch.exists() or not list(repo.scratch.iterdir())
    refused("h1", 96, "HEAD is not", "d" * 40, "h1")
    refused("h2", 96, "HEAD is not", repo.head[:12], "h2")
    (repo.root / ".gitignore").write_text("runs/\n# edited\n")
    refused("c1", 98, "tracked files differ from the commit")
    repo.git("checkout", "-q", ".gitignore")
    refused("n1", 94, "run from the root of the checkout", cwd=repo.root / "experiments")
    ledger = repo.root / "experiments" / "harness" / "TRIES.md"
    for i, text in enumerate(("| H34 (a row) | heart | x |\n", h35_row() * 2, h35_row(cells=10),
                              h35_row(first=FIRST.replace("ten tries", "eight tries")),
                              h35_row().replace("| H35 (", "| H36 ("))):
        ledger.write_text(text)
        repo.git("commit", "-q", "-am", f"ledger {i}")
        refused(f"t{i}", 93, "must hold exactly one row H35 of eleven cells", repo.git("rev-parse", "HEAD"), f"t{i}")
    ledger.write_text(h35_row("OUTCOME H35: VOID: v3: no gate record; nothing is built."))
    repo.git("commit", "-q", "-am", "an outcome")
    refused("u1", 95, "already holds an outcome", repo.git("rev-parse", "HEAD"), "u1")
    ledger.write_text(h35_row())
    repo.git("commit", "-q", "-am", "the row")
    head = repo.git("rev-parse", "HEAD")
    (repo.root / OUT_DIR).mkdir(parents=True)
    (repo.root / OUT_DIR / "reading.json").write_text("the one reading\n")
    refused("k1", 95, "already holds a reading", head, "k1")
    (repo.root / OUT_DIR / "reading.json").rename(repo.tmp / "aside.json")
    refused("v0", 99, "no virtual environment at", head, "v0", H35_VENV=str(repo.tmp / "none"))
    refused("d0", 97, "no T2 stages under", head, "d0", H35_DATA=str(repo.tmp / "none"))
    assert _uv_calls(repo) == [] and _once(repo) == {}  # none of these refusals started anything
    assert _launch(repo, head, "ok")[0] == 0
    assert _once(repo) == {"g0": "0", "gate": "0", "run": "0", "read": "0"}
