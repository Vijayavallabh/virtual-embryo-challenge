"""The launcher of try K30 (experiments/t3_decomp/run_k30.sh): its gate's refusals, its steps' order and exit records,
the S3 weight read from row K30 and checked against ruling row R53, the VOID reading after a failed guard or past the
cutoff. Each test runs the launcher in a throwaway git repository under the test's own temporary folder, with a
hand-made ledger, a hand-made docs/S2_RULES.md and its pin, dummy committed inputs and a `uv` stub that stands in for
every Python call (nothing is scored, read or built). No competition file, no committed ledger, no remote machine."""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
SH = REPO / "experiments" / "t3_decomp" / "run_k30.sh"
SH_REL = "experiments/t3_decomp/run_k30.sh"
SCRIPT_REL = "experiments/t3_decomp/k30_depth_read.py"
OUT_REL = "experiments/t3_decomp/results/k30_depth_read"
TABLE_REL = "experiments/t3_decomp/results/gata4/gata4_composition.csv"
K18_REL = "experiments/t3_decomp/results/k18_built_files/tasks"
RUN_DIR = "v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
ROW = ("| K30 (one row, one script `experiments/t3_decomp/k30_depth_read.py` and its launcher "
       "`experiments/t3_decomp/run_k30.sh <COMMIT> <RUN>`; S3 weight `--s3-weight plain`; three tries) | T3 | the "
       "row's other cells |\n")
R53 = ("| R53 (the captain's order of 2026-11-24, ruling, no run: \"Yes: ...\"; row K30's world weight is `plain`; not "
       "a try) | T3 | the ruling's other cells |\n")
LEDGER = "| K29 (another row) | T3 | x |\n" + R53 + ROW
S2_TEXT = "synthetic second-upload rules for the launcher's test\n"
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES")


def test_the_launcher_text():
    sh = SH.read_text()
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    assert ":?" not in sh and not re.search(r"\$\{[^}]*:\?[^}]*'[^}]*\}", sh)  # no ${var:?message}, so no apostrophe
    assert not re.search(r"(^|[\s;|&(])rm(\s|$)", sh, flags=re.M)  # nothing deleted
    threads = [int(x) for x in re.findall(r"^\s*threads (\d+)", sh, flags=re.M)]
    assert threads and max(threads) <= 16 and "--workers 8 --threads 2" in sh
    assert f"SCRIPT={SCRIPT_REL}" in sh and f"OUT={OUT_REL}" in sh and f"TABLE={TABLE_REL}" in sh
    assert 'S="$HOME/v/scratch/k30/$RUN"' in sh and 'R="runs/k30_$RUN"' in sh
    assert f"RUN_DIR=$HOME/{RUN_DIR}" in sh
    for step, prev in (("gate", '""'), ("tests", "gate"), ("guard", "tests"), ("run", "guard"), ("read", "run")):
        assert f"begin {step} {prev}" in sh
    order = [sh.index(f"begin {s} ") for s in ("gate", "tests", "guard", "run", "read")]
    assert order == sorted(order)
    for t in ("tests/test_k30_depth_read.py", "tests/test_run_k30.py", "tests/test_k18_built_files.py"):
        assert t in sh
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh and "set -u" in sh
    assert sh.count('--s3-weight "$WEIGHT"') == 4  # two kinds of part, the void reading and the reading
    gate = sh.split("begin gate", 1)[1].split("begin tests", 1)[0]
    codes = [int(x) for x in re.findall(r"finish (\d+)", gate)]
    assert codes[:3] == [96, 98, 97] and set(codes) == {96, 98, 97, 95, 89, 94}  # HEAD first, then the tree, then S2


def _git_env() -> dict:
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                        "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, env=_git_env(),
                       text=True)
    return r.stdout.strip()


UV_STUB = """#!{py} -I
# uv run --frozen --no-sync python -I -B ARGS: the tests exit $HOME/tests_rc; a -c one-liner runs in python -I -B; every
# call of the script is logged in $HOME/script_calls; `guard` exits $HOME/guard_rc when it is there (writing nothing),
# else writes guard.json (failing, exit 6, when $HOME/guard_fails is there); `run` exits $HOME/run_rc and `read`
# $HOME/read_rc (0 when absent)
import json, os, sys
from pathlib import Path
a = sys.argv[1:]
if a[:6] != ["run", "--frozen", "--no-sync", "python", "-I", "-B"]:
    sys.exit(97)
a = a[6:]
home = Path(os.environ["HOME"])
def code(name):
    p = home / name
    return int(p.read_text()) if p.exists() else 0
if a[0] == "-m":
    sys.exit(code("tests_rc"))
if a[0] == "-c":
    os.execv(sys.executable, [sys.executable, "-I", "-B", *a])
with open(home / "script_calls", "a") as fh:
    fh.write(" ".join(a) + "\\n")
opt = dict(zip(a[2::2], a[3::2]))
if a[1] == "guard" and code("guard_rc"):
    sys.exit(code("guard_rc"))
if a[1] == "guard":
    ok = not (home / "guard_fails").exists()
    Path(opt["--runs"]).mkdir(parents=True, exist_ok=True)
    Path(opt["--runs"], "guard.json").write_text(json.dumps({"passes": ok, "provenance": {"commit": opt["--commit"]}}))
    sys.exit(0 if ok else 6)
sys.exit(code(a[1] + "_rc"))
"""


def _launcher_repo(tmp_path: Path, ledger: str = LEDGER) -> SimpleNamespace:
    repo, home, bindir = tmp_path / "repo", tmp_path / "home", tmp_path / "bin"
    (home / RUN_DIR).mkdir(parents=True)
    s2_sha = hashlib.sha256(S2_TEXT.encode()).hexdigest()
    files = {SH_REL: SH.read_text(), "experiments/harness/TRIES.md": ledger, ".gitignore": "runs/\n",
             "docs/S2_RULES.md": S2_TEXT, "docs/S2_RULES.sha256": f"{s2_sha}  docs/S2_RULES.md\n",
             TABLE_REL: "label,wt_cells,rep1_depth_ratio,rep2_depth_ratio\n",
             **{f"{K18_REL}/{f}": "{}\n" for f in ("plain_g1.json", "plain_g2.json", "dm.json")}}
    repo.mkdir()
    _git(repo, "init", "-q")
    for rel, t in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t)
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "init")
    bindir.mkdir()
    (bindir / "uv").write_text(UV_STUB.replace("{py}", sys.executable))
    (bindir / "uv").chmod(0o755)
    return SimpleNamespace(repo=repo, home=home, bin=bindir, head=_git(repo, "rev-parse", "HEAD"))


def _commit(L, rel: str, text_: str) -> str:
    (L.repo / rel).write_text(text_)
    _git(L.repo, "commit", "-q", "-am", f"edit {rel}")
    return _git(L.repo, "rev-parse", "HEAD")


def _launch(L, *args) -> tuple[int, str]:
    env = {**_git_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}",
           "GIT_CEILING_DIRECTORIES": str(L.repo.parent)}
    r = subprocess.run(["bash", SH_REL, *args], cwd=L.repo, env=env, capture_output=True, text=True, timeout=120)
    return r.returncode, r.stderr


def _exit(L, run: str, step: str) -> str:
    p = L.home / "v" / "scratch" / "k30" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else "absent"


def _calls(L) -> list[str]:
    p = L.home / "script_calls"
    return p.read_text().splitlines() if p.exists() else []


def _kinds(L) -> list[str]:
    return [c.split()[1] for c in _calls(L)]


def test_usage_and_the_run_name(tmp_path):
    L = _launcher_repo(tmp_path)
    for args in ((), (L.head,), (L.head, "r", "plain", "x"), (L.head, "r", "0.7"), (L.head, "r", ""), ("", "r")):
        code, err = _launch(L, *args)
        assert code == 2 and "usage: run_k30.sh COMMIT RUN [WEIGHT]" in err, args
    for bad in (".hidden", "a/b", "x y"):
        assert _launch(L, L.head, bad)[0] == 90
    assert _calls(L) == [] and not (L.home / "v" / "scratch").exists()


def test_the_gates_first_refusals_in_order(tmp_path):
    L = _launcher_repo(tmp_path)
    code, err = _launch(L, "0" * 40, "r1")
    assert code == 96 and "HEAD is not" in err
    assert _exit(L, "r1", "done") == "96" and _exit(L, "r1", "gate") == "96" and _exit(L, "r1", "tests") == "absent"
    code, err = _launch(L, L.head, "r1")  # r1 was made for another commit
    assert code == 91 and "another commit" in err
    (L.repo / "experiments/harness/TRIES.md").write_text(LEDGER + "an edit\n")
    code, err = _launch(L, L.head, "r2")
    assert code == 98 and "not clean" in err and _exit(L, "r2", "gate") == "98"
    _git(L.repo, "checkout", "--", "experiments/harness/TRIES.md")
    head = _commit(L, "docs/S2_RULES.md", S2_TEXT + "an amendment without its pin\n")
    code, err = _launch(L, head, "r3")
    assert code == 97 and "docs/S2_RULES.sha256 does not check" in err and _exit(L, "r3", "gate") == "97"
    assert _calls(L) == []  # no refusal called the script


SCRIPT_P, LAUNCHER_P = "`experiments/t3_decomp/k30_depth_read.py`", "`experiments/t3_decomp/run_k30.sh"


@pytest.mark.parametrize("ledger,code,why", [
    (R53, 95, "exactly one row K30"),
    (R53 + ROW * 2, 95, "exactly one row K30"),
    (R53 + ROW.replace("| K30 (", "| K30a ("), 95, "exactly one row K30"),
    (R53 + ROW.replace("| K30 (", "| K30 [", 1), 95, "does not start 'K30 ('"),
    (R53 + ROW.replace(SCRIPT_P, "`experiments/t3_decomp/k29_cd7_read.py`"), 95, "k30_depth_read"),
    (R53 + ROW.replace(LAUNCHER_P, "`experiments/t3_decomp/run_k29.sh"), 95, "run_k30"),
    (R53 + ROW.replace("k30_depth_read.py`", "k30_depth_read.py.bak`"), 95, "k30_depth_read"),
    (R53 + ROW.replace(LAUNCHER_P, "`old/experiments/t3_decomp/run_k30.sh"), 95, "run_k30"),
    (R53 + f"| K30 (one row; script and launcher in cell 3; `--s3-weight plain`) | {SCRIPT_P} {LAUNCHER_P}` |\n",
     95, "k30_depth_read"),
    (R53 + ROW.replace("S3 weight `--s3-weight plain`; ", ""), 95, "exactly one S3 weight"),
    (R53 + ROW.replace("`--s3-weight plain`", "--s3-weight plain"), 95, "exactly one S3 weight"),
    (R53 + ROW.replace("`--s3-weight plain`", "`--s3-weight 0.7`"), 95, "exactly one S3 weight"),
    (R53 + ROW.replace("`--s3-weight plain`", "`--s3-weight plain` and `--s3-weight 0.5`"), 95, "exactly one S3"),
    (R53 + ROW.replace("`--s3-weight plain`", "`--s3-weight plain --s3-weight 0.5`"), 95, "exactly one S3 weight"),
    (ROW, 89, "exactly one row R53"),
    (R53 * 2 + ROW, 89, "exactly one row R53"),
    (R53.replace("R53 (the captain's ", "R53 (a ruling of ") + ROW, 89, "does not start"),
    (R53.replace("row K30's world weight is `plain`; ", "") + ROW, 89, "does not name row K30's world weight"),
    (R53.replace("`plain`", "`0.5`") + ROW, 89, "does not name row K30's world weight as plain"),
    (R53.replace("is `plain`", "is `plain`, row K30's world weight is `0.5`") + ROW, 89, "world weight as plain"),
])
def test_the_gate_refuses_rows_k30_and_r53_not_as_required(tmp_path, ledger, code, why):
    L = _launcher_repo(tmp_path, ledger)
    got, err = _launch(L, L.head, "r")
    assert got == code and why in err and _exit(L, "r", "gate") == str(code) and _calls(L) == [], err


def test_the_weight_is_read_from_the_row_and_a_given_one_must_be_it(tmp_path):
    L = _launcher_repo(tmp_path)
    code, err = _launch(L, L.head, "w1", "0.5")
    assert code == 95 and "the weight given, 0.5, is not the one row K30 names (plain)" in err and _calls(L) == []
    assert _launch(L, L.head, "w2", "plain")[0] == 0
    assert (L.home / "v/scratch/k30/w2/weight").read_text().strip() == "plain"
    assert all("--s3-weight plain" in c for c in _calls(L) if c.split()[1] in ("run", "read"))
    L2 = _launcher_repo(tmp_path / "half", R53.replace("`plain`", "`0.5`") + ROW.replace("`--s3-weight plain`",
                                                                                          "`--s3-weight 0.5`"))
    assert _launch(L2, L2.head, "h")[0] == 0
    assert all("--s3-weight 0.5" in c for c in _calls(L2) if c.split()[1] in ("run", "read"))


@pytest.mark.parametrize("what", ["k18", "table", "run_dir"])
def test_the_gate_refuses_missing_inputs(tmp_path, what):
    L = _launcher_repo(tmp_path)
    head = L.head
    if what == "run_dir":
        (L.home / RUN_DIR).rename(L.home / RUN_DIR.replace("gata4_", "elsewhere_"))
        why = "gata4 run folder"
    else:
        rel = f"{K18_REL}/dm.json" if what == "k18" else TABLE_REL
        _git(L.repo, "rm", "-q", "--cached", rel)
        _git(L.repo, "commit", "-q", "-m", f"untrack {rel}")
        head, why = _git(L.repo, "rev-parse", "HEAD"), f"{rel} is not tracked"
    code, err = _launch(L, head, "i")
    assert code == 94 and why in err and _exit(L, "i", "gate") == "94" and _calls(L) == []


def test_failing_tests_stop_the_call(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.home / "tests_rc").write_text("1")
    assert _launch(L, L.head, "t")[0] == 10 and _exit(L, "t", "tests") == "10" and _exit(L, "t", "gate") == "0"
    assert _calls(L) == []


def test_the_launcher_runs_gate_tests_guard_run_and_read_in_order_with_exit_records(tmp_path):
    L = _launcher_repo(tmp_path)
    code, err = _launch(L, L.head, "k30a")
    assert code == 0, err
    for step in ("gate", "tests", "guard", "run", "read", "done"):
        assert _exit(L, "k30a", step) == "0", step
        if step != "done":
            assert (L.home / "v/scratch/k30/k30a" / f"{step}.head").read_text().strip() == L.head
    calls = _calls(L)
    assert _kinds(L) == ["run"] * 5 + ["guard"] + ["run"] * 5 + ["read"]
    assert all(c.startswith(f"{SCRIPT_REL} ") and f"--commit {L.head}" in c for c in calls)
    assert all("--part guard" in c for c in calls[0:5]) and all("--part score" in c for c in calls[6:11])
    assert sorted(c.split("--world ")[1].split()[0] for c in calls[0:5]) == ["dm", "plain", "plain", "plain", "plain"]
    plain_parts = [re.search(r"--task (g\d) --seeds (\S+)", c).groups() for c in calls[6:11] if "--world plain" in c]
    assert sorted(plain_parts) == [("g1", "0-4"), ("g1", "5-9"), ("g2", "0-4"), ("g2", "5-9")]
    assert sum("--world dm --workers 8 --threads 2 --s3-weight plain" in c for c in calls) == 2
    assert all(f"--run-dir {L.home / RUN_DIR}" in c for c in calls if c.split()[1] == "run")
    assert calls[5] == f"{SCRIPT_REL} guard --commit {L.head} --runs runs/k30_k30a"
    assert calls[11] == f"{SCRIPT_REL} read --commit {L.head} --runs runs/k30_k30a --s3-weight plain --out {OUT_REL}"
    (L.home / "script_calls").unlink()
    code, err = _launch(L, L.head, "k30a")  # the same RUN at the same commit: guard.json read again
    assert code == 0, err
    assert _kinds(L) == ["run"] * 5 + ["read"]


def test_a_failed_guard_makes_the_void_reading_once_and_starts_no_scoring_part(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.home / "guard_fails").write_text("")
    code, err = _launch(L, L.head, "g")
    assert code == 26 and _exit(L, "g", "done") == "26" and _exit(L, "g", "guard") == "26", err
    assert _exit(L, "g", "run") == "absent" and _exit(L, "g", "read") == "absent" and _exit(L, "g", "void_read") == "0"
    assert _kinds(L) == ["run"] * 5 + ["guard", "read"]
    assert _calls(L)[-1] == f"{SCRIPT_REL} read --commit {L.head} --runs runs/k30_g --s3-weight plain --out {OUT_REL}"
    (L.repo / OUT_REL).mkdir(parents=True)
    (L.repo / OUT_REL / "reading.json").write_text("{}")
    (L.home / "script_calls").unlink()
    assert _launch(L, L.head, "g")[0] == 26 and _calls(L) == []  # the reading is there: 26 again, nothing called
    R = L.repo / "runs" / "k30_g"
    (R / "guard.json").write_text(json.dumps({"passes": True, "provenance": {"commit": "f" * 40}}))
    assert _launch(L, L.head, "g")[0] == 21 and _calls(L) == []
    (R / "guard.json").write_text("not json")
    assert _launch(L, L.head, "g")[0] == 20 and _calls(L) == []


@pytest.mark.parametrize("rc,want", [(2, 22), (3, 30), (1, 20)])
def test_a_refused_or_crashed_guard_check_ends_the_call_without_a_reading(tmp_path, rc, want):
    L = _launcher_repo(tmp_path)
    (L.home / "guard_rc").write_text(str(rc))
    assert _launch(L, L.head, "c")[0] == want and _exit(L, "c", "guard") == str(want)
    assert _kinds(L) == ["run"] * 5 + ["guard"]


def test_a_failed_part_stops_the_call_before_the_next_step(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.home / "run_rc").write_text("1")
    assert _launch(L, L.head, "a")[0] == 30 and _exit(L, "a", "guard") == "30"
    assert _kinds(L) == ["run"] * 5  # the guard parts only: no check, no scoring part
    (L.home / "script_calls").unlink()
    R = L.repo / "runs" / "k30_b"
    R.mkdir(parents=True)
    (R / "guard.json").write_text(json.dumps({"passes": True, "provenance": {"commit": L.head}}))
    assert _launch(L, L.head, "b")[0] == 30 and _exit(L, "b", "run") == "30" and _exit(L, "b", "read") == "absent"
    assert _kinds(L) == ["run"] * 5


@pytest.mark.parametrize("rc,want", [(3, 43), (4, 44), (2, 42), (1, 41)])
def test_the_read_steps_codes(tmp_path, rc, want):
    L = _launcher_repo(tmp_path)
    (L.home / "read_rc").write_text(str(rc))
    assert _launch(L, L.head, "r")[0] == want and _exit(L, "r", "read") == str(want) and _exit(L, "r", "run") == "0"


def _passing(L, run: str) -> None:
    R = L.repo / "runs" / f"k30_{run}"
    R.mkdir(parents=True)
    (R / "guard.json").write_text(json.dumps({"passes": True, "provenance": {"commit": L.head}}))


@pytest.mark.parametrize("stub,passing,kinds,step", [
    ({"run_rc": "7"}, False, ["run"] * 5 + ["read"], "guard"),  # the guard parts
    ({"guard_rc": "7"}, False, ["run"] * 5 + ["guard", "read"], "guard"),  # the guard check
    ({"run_rc": "7"}, True, ["run"] * 5 + ["read"], "run"),  # the scoring parts
])
def test_a_refusal_past_the_cutoff_makes_the_reading_once_and_finishes_27(tmp_path, stub, passing, kinds, step):
    L = _launcher_repo(tmp_path)
    for name, value in stub.items():
        (L.home / name).write_text(value)
    if passing:
        _passing(L, "late")
    code, err = _launch(L, L.head, "late")
    assert code == 27 and _exit(L, "late", step) == "27" and _exit(L, "late", "void_read") == "0", err
    assert _kinds(L) == kinds
    want = f"{SCRIPT_REL} read --commit {L.head} --runs runs/k30_late --s3-weight plain --out {OUT_REL}"
    assert _calls(L)[-1] == want


def test_a_second_run_after_a_reading_scores_no_arm(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.repo / OUT_REL).mkdir(parents=True)
    (L.repo / OUT_REL / "reading.json").write_text("{}")
    code, err = _launch(L, L.head, "second")  # a fresh RUN: no guard part starts
    assert code == 44 and _exit(L, "second", "guard") == "44" and _calls(L) == [], err
    _passing(L, "third")  # a RUN whose guard passed: no scoring part starts
    code, err = _launch(L, L.head, "third")
    assert code == 44 and _exit(L, "third", "run") == "44" and _calls(L) == [], err
