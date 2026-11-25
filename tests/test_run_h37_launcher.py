"""experiments/harness/run_h37.sh, row H37's launcher (the gate of row R69's file WBH). Every call runs in a throwaway
git repository under the test's own folder, with a temporary home and stand-ins for uv and sha256sum first on PATH: no
launcher runs with a real commit, nothing reads a competition file, nothing is built and nothing is scored. Each of
the launcher's own refusals is executed, then a whole call through the stand-ins, then the steps it may start only
once (88), the script's refusal (89), parts made outside it (92), a G0 STOP and a failing tests step. The banked
identity step (build_test.py, then halve.py) stops on another sha256 and on a stamp that is not this commit's."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = REPO_ROOT / "experiments" / "harness" / "run_h37.sh"
ROWS = REPO_ROOT / "experiments" / "harness" / "h37_row.py"
FIRST = ("H37 (the gate of row R69's file WBH, the heart's world-blur hedge: six counted tries, arms W, U and P "
         "on two transitions; G0 first, target-free, not a try)")
R69 = ("| R69 (the captain's order of 2026-11-25, ruling, no run: \"File WBH may use ...\", a synthetic stand-in; "
       "built only after its gate row H37 reads PASS, only with `--captain-ruling R69`; not a try) | heart | x |\n")
BANKED = {"T2heart_organ_ts125_tk.h5ad": "d4d8d75f8bed48813bf143f5d7726828e5c30cdb27bd708b6fc2d9a455d21893",
          "T2heart_organ_ts125_tk_half.h5ad": "9835fa5d47d7021760beb2056ce44c2c846ded6fc31dacf64bbb2d9f0f7d1b54",
          "T2heart_organ_ts125_tk_wb_half.h5ad": "1" * 64}

# uv: logs its arguments; makes the files a builder would write (build_test.py --out, halve.py and world_blur.py OUT)
# and a run step's part;
# prints STUB_STAMP for a stamp read; exits STUB_PYTEST for pytest, STUB_REFUSALS for `refusals`, STUB_<STEP> for a
# step of the script (g0, gate, run, read) and 0 otherwise
STUB_UV = r"""#!/usr/bin/env bash
echo "$*" >> "$HOME/uvcalls"
args=("$@")
for i in "${!args[@]}"; do
  a=${args[$i]}
  case "$a" in
    *.stamp.json) echo "${STUB_STAMP:-none}"; exit 0 ;;
    --out) o=${args[$((i + 1))]}
           case " $* " in
             *" refusals "*) ;;
             *"h37_heart_world_blur.py run "*) mkdir -p "$o"; echo "{}" > "$o/h37.json" ;;
             *) case "$o" in *.h5ad) mkdir -p "$(dirname "$o")"; echo built > "$o" ;; esac ;;
           esac ;;
    *halve.py|*world_blur.py) out=${args[$((i + 2))]}; mkdir -p "$(dirname "$out")"; echo built > "$out" ;;
  esac
done
case " $* " in
  *" pytest "*) exit "${STUB_PYTEST:-0}" ;;
  *" refusals "*) exit "${STUB_REFUSALS:-0}" ;;
  *"h37_heart_world_blur.py g0 "*) exit "${STUB_G0:-0}" ;;
  *"h37_heart_world_blur.py gate "*) exit "${STUB_GATE:-0}" ;;
esac
exit 0
"""
# sha256sum: the recorded hashes for the banked file's names (STUB_SHA_FILE overrides the file's), the real one with -c
STUB_SHA = r"""#!/usr/bin/env bash
if [ "${1:-}" = -c ]; then exec __REAL__ "$@"; fi
for f in "$@"; do
  case "$(basename "$f")" in
    T2heart_organ_ts125_tk.h5ad) echo "${STUB_SHA_FILE:-__FILE__}  $f" ;;
    T2heart_organ_ts125_tk_half.h5ad) echo "${STUB_SHA_HALF:-__HALF__}  $f" ;;
    *) echo "__WB__  $f" ;;
  esac
done
"""


def _env(tmp_path, extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    b = tmp_path / "bin"
    if not b.exists():
        b.mkdir()
        (b / "uv").write_text(STUB_UV)
        sha = (STUB_SHA.replace("__REAL__", shutil.which("sha256sum"))
               .replace("__FILE__", BANKED["T2heart_organ_ts125_tk.h5ad"])
               .replace("__HALF__", BANKED["T2heart_organ_ts125_tk_half.h5ad"])
               .replace("__WB__", BANKED["T2heart_organ_ts125_tk_wb_half.h5ad"]))
        (b / "sha256sum").write_text(sha)
        for f in ("uv", "sha256sum"):
            (b / f).chmod(0o755)
    env.update(HOME=str(tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path),
               PATH=f"{b}{os.pathsep}{os.environ['PATH']}")
    env.update(extra or {})
    return env


def git_fn(repo, tmp_path):
    def git(*a):
        return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "init.defaultBranch=x", *a],
                              cwd=repo, env=_env(tmp_path), capture_output=True, text=True, check=True).stdout.strip()
    return git


def h37_line(outcome="pending", first=FIRST) -> str:
    return f"| {first} | heart | x | x | x | x | x | x | x | x | {outcome} |\n"


def make_repo(tmp_path, ledger=None, files=()) -> tuple[Path, object, str]:
    """A throwaway repository holding the row module, a stand-in script and a ledger (default: row R69 and an open
    row H37), a stand-in virtual environment and data folder (both ignored), committed; returns it, git and HEAD."""
    repo = tmp_path / "repo"
    (repo / "experiments" / "harness").mkdir(parents=True)
    (tmp_path / "home").mkdir(exist_ok=True)
    shutil.copy(ROWS, repo / "experiments" / "harness" / "h37_row.py")
    (repo / "experiments" / "harness" / "h37_heart_world_blur.py").write_text("x = 1\n")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(R69 + h37_line() if ledger is None else ledger)
    for f in files:
        (repo / f).parent.mkdir(parents=True, exist_ok=True)
        (repo / f).write_text("x = 1\n")
    (repo / ".gitignore").write_text("runs/\n.venv/\ndata/\nexperiments/harness/results/\n")
    (repo / ".venv" / "bin").mkdir(parents=True)
    (repo / ".venv" / "bin" / "python").write_text("#!/bin/sh\n")
    (repo / ".venv" / "bin" / "python").chmod(0o755)
    for d in ("T2", "T3"):
        (repo / "data" / "raw" / d).mkdir(parents=True)
    git = git_fn(repo, tmp_path)
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "row")
    return repo, git, git("rev-parse", "HEAD")


def launch(tmp_path, *args, extra=None):
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=tmp_path / "repo", capture_output=True, text=True,
                       env=_env(tmp_path, extra), timeout=180)
    return r.returncode, r.stderr


def scratch(tmp_path, run=None) -> Path:
    base = tmp_path / "home" / "v" / "scratch" / "h37"
    return base / run if run else base


def exit_of(tmp_path, run, step):
    p = scratch(tmp_path, run) / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


# ------------------------------------------------------------------------------------------------ arguments
def test_odd_arguments_are_refused(tmp_path):
    _, _, head = make_repo(tmp_path)
    assert launch(tmp_path, head)[0] != 0
    assert launch(tmp_path, head, "r1", "x")[0] == 90
    assert launch(tmp_path, head, "r 1")[0] == 90
    assert launch(tmp_path, head[:12], "r1")[0] == 90
    assert launch(tmp_path, head.upper(), "r1")[0] == 90
    assert not scratch(tmp_path).exists()


def test_a_run_folder_is_never_reused(tmp_path):
    _, _, head = make_repo(tmp_path)
    launch(tmp_path, head, "r1", extra={"STUB_PYTEST": "1"})
    assert launch(tmp_path, head, "r1")[0] == 91


# ------------------------------------------------------------------------------------------------ the check step
def test_head_and_a_clean_tree(tmp_path):
    repo, git, head = make_repo(tmp_path)
    assert launch(tmp_path, "b" * 40, "r1")[0] == 96 and exit_of(tmp_path, "r1", "check") == "96"
    (repo / "experiments" / "harness" / "TRIES.md").write_text(R69 + h37_line() + "\n")
    assert launch(tmp_path, head, "r2")[0] == 98
    assert exit_of(tmp_path, "r2", "done") == "98"


@pytest.mark.parametrize("ledger,code", [
    (R69, 93),                                                   # no row H37
    (R69 + h37_line(first="H37 (another first cell)"), 93),      # not the fixed first cell
    (R69 + h37_line() + h37_line(), 93),                         # two rows H37
    (R69 + h37_line("OUTCOME H37: STOP: IDLE."), 95),           # the gate was read
    (h37_line(), 93),                                            # no row R69
    (R69.replace("ruling, no run", "ruling") + h37_line(), 93),  # R69 not the captain's order
])
def test_the_rows(tmp_path, ledger, code):
    _, _, head = make_repo(tmp_path, ledger)
    assert launch(tmp_path, head, "r1")[0] == code


def test_root_reading_venv_and_data(tmp_path):
    repo, git, head = make_repo(tmp_path)
    out = repo / "experiments" / "harness" / "results" / "h37_heart_world_blur"
    out.mkdir(parents=True)
    (out / "reading.json").write_text("{}")
    assert launch(tmp_path, head, "r1")[0] == 95
    (out / "reading.json").unlink()
    assert launch(tmp_path, head, "r2", extra={"H37_VENV": str(tmp_path / "none")})[0] == 99
    (tmp_path / "d" / "T2").mkdir(parents=True)
    assert launch(tmp_path, head, "r3", extra={"H37_DATA": str(tmp_path / "d")})[0] == 97
    git("rm", "-q", "experiments/harness/h37_heart_world_blur.py")
    git("commit", "-q", "-m", "no script")
    assert launch(tmp_path, git("rev-parse", "HEAD"), "r4")[0] == 94


# ------------------------------------------------------------------------------------------------ the steps
def test_a_failing_tests_step_stops_the_call(tmp_path):
    _, _, head = make_repo(tmp_path)
    assert launch(tmp_path, head, "r1", extra={"STUB_PYTEST": "1"})[0] == 1
    assert exit_of(tmp_path, "r1", "check") == "0" and exit_of(tmp_path, "r1", "tests") == "1"
    assert exit_of(tmp_path, "r1", "banked") is None


def test_the_banked_identity_step_stops_on_other_bytes_and_stamps(tmp_path):
    _, _, head = make_repo(tmp_path)
    stamp = {"STUB_STAMP": f"{head} True"}
    assert launch(tmp_path, head, "r1", extra={**stamp, "STUB_SHA_FILE": "0" * 64})[0] == 20
    assert launch(tmp_path, head, "r2", extra={**stamp, "STUB_SHA_HALF": "0" * 64})[0] == 21
    assert launch(tmp_path, head, "r3", extra={"STUB_STAMP": f"{head} False"})[0] == 15
    assert exit_of(tmp_path, "r3", "g0") is None and not (scratch(tmp_path) / "g0.exit").exists()


def test_a_whole_call_then_every_step_only_once(tmp_path):
    repo, _, head = make_repo(tmp_path)
    stamp = {"STUB_STAMP": f"{head} True"}
    assert launch(tmp_path, head, "r1", extra=stamp)[0] == 0
    for step in ("check", "tests", "banked", "g0", "gate", "run", "read", "done"):
        assert exit_of(tmp_path, "r1", step) == "0", step
    for step in ("g0", "gate", "run", "read"):
        assert (scratch(tmp_path) / f"{step}.exit").read_text().strip() == "0"
        assert (scratch(tmp_path) / f"{step}.head").read_text().split() == [head, "r1"]
    calls = (tmp_path / "home" / "uvcalls").read_text()
    assert "experiments/builds/build_test.py heart --entry organ_scope --trend-strength 1.25 --trend-form category " \
           "--z-mode raw --captain-ruling R55" in calls
    assert calls.count("h37_heart_world_blur.py run --task") == 10 and calls.count("refusals --step") == 4
    assert launch(tmp_path, head, "r2", extra=stamp)[0] == 88
    assert exit_of(tmp_path, "r2", "g0") == "88"


def test_the_script_refusal_parts_and_a_g0_stop(tmp_path):
    repo, _, head = make_repo(tmp_path)
    stamp = {"STUB_STAMP": f"{head} True"}
    assert launch(tmp_path, head, "r1", extra={**stamp, "STUB_REFUSALS": "1"})[0] == 89
    assert not (scratch(tmp_path) / "g0.exit").exists()
    assert launch(tmp_path, head, "r2", extra={**stamp, "STUB_G0": "3"})[0] == 3
    assert (scratch(tmp_path) / "g0.exit").read_text().strip() == "3" and exit_of(tmp_path, "r2", "gate") is None


def test_parts_made_outside_the_launcher(tmp_path):
    repo, _, head = make_repo(tmp_path)
    (repo / "experiments" / "harness" / "results" / "h37_heart_world_blur" / "parts").mkdir(parents=True)
    assert launch(tmp_path, head, "r1", extra={"STUB_STAMP": f"{head} True"})[0] == 92
    assert exit_of(tmp_path, "r1", "gate") == "0" and not (scratch(tmp_path) / "run.exit").exists()
    assert launch(tmp_path, head, "r2", extra={"STUB_STAMP": f"{head} True", "STUB_GATE": "4"})[0] == 88
