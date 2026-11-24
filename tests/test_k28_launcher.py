"""Tests for experiments/t3_decomp/run_k28.sh (row K28's launcher, step g0a): every gate, each exit mapping, the
once-only stop and what it prints, traced in a temporary git repository under a temporary HOME with a stub `uv` whose
exit status per step the test sets. No Python step of the row runs; no data of any kind is read."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]
SH = REPO / "experiments" / "t3_decomp" / "run_k28.sh"
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES", "GIT_CEILING_DIRECTORIES")
ROW = "| K28 (one row, one script) | T3 | `experiments/t3_decomp/k28_chip_partners.py` |\n"
WORDS = (  # the amendment's whole text
    "Reopened by the captain's order of 2026-11-24 for the database-only gate of one practice row (K28, T3 candidate "
    "CHIPQ-C); a fold is read only under the captain's ruling R33 lifting row K27's closing sentence for that row "
    "alone. No new fetch, no new kind of file: the files pinned on 2026-11-23, verified. New use: for a query gene "
    "named to code, code picks three other antigens whose binding within 5 kb of transcription start sites overlaps "
    "the query's most specifically in one cell-type class, by a statistic computed from kept experiments only; panel "
    "genes those antigens bind form the table, under row K27's one sign convention. Experiments are wild type and "
    "untreated as far as the public metadata say: the K27 kept list for the two query genes; for other antigens the "
    "committed keyword filter, two further classes in code (tumour-derived or transformed lines; cell-type classes "
    "where pathway agents are common) and two screening readers on every row that can enter a table, antigen names "
    "masked. No experiment perturbing either held-out gene and no expression measurement enters. Not checked: "
    "genotypes and culture conditions the metadata do not state. No person or agent reads a partner name or a bound "
    "gene. Practice read only. Licence unchanged (CC BY-SA 4.0): share-alike for a derived table in any distributed "
    "file or report. No file without the captain's rulings on provenance (PLAN 4.4), A15 and the team's pathway-tier "
    "bar."
)
OPENING = "Reopened by the captain's order of 2026-11-24"  # the amendment's opening words alone
ONE_WORD = WORDS.replace(" bar.", " rule.")  # the whole text with one word, its last, changed
REGISTER = f"| X018 | a generic database of public ChIP-seq | fetched 2026-11-23. {WORDS} |\n"
FETCH_DIR = "chip" "atlas"  # the fetch folder the launcher checks (spelt apart: no file here reaches the database)
OUT = "experiments/t3_decomp/results/k28_chip_partners"
BLOCK = ("G0a: classes excluded by class (vi): 2\nG0a Ctnnb1: usable experiments in its class: 4; candidates: 9\n"
         "G0a Gata4: usable experiments in its class: 6; candidates: 13\n")

# The stub `uv`: a python -c one-liner runs for real (the launcher's json read); the tests step and the gate exit with
# the code the test sets, both with noise on stdout or stderr; the gate prints its block at exit 0 and 3 and, at exit 3,
# leaves the reading the real step would leave (reading.json naming this HEAD, and reading.txt).
STUB = """#!{py}
import json, os, subprocess, sys
a = sys.argv[1:]
log = os.environ.get("STUB_ARGS")
if log:
    open(log, "a").write(" ".join(a) + "\\n")
if "-c" in a:
    os.execv(sys.executable, [sys.executable, *a[4:]])
rc = lambda k: int(os.environ.get(k, "0"))
if "pytest" in a:
    print("noise from the tests: Ant01 gs00017")
    print("more noise from the tests", file=sys.stderr)
    sys.exit(rc("STUB_TESTS_RC"))
script = next((x for x in a if x.endswith(".py")), "")
step = a[a.index(script) + 1] if script else ""
if script.endswith("k28_chip_partners.py") and step == "g0a":
    r = rc("STUB_G0A_RC")
    print("refused: the words of a refusal", file=sys.stderr)
    if r in (0, 3):
        sys.stdout.write({block!r})
    else:
        print("words on stdout before a refusal or a crash")
    if r == 3 and not os.environ.get("STUB_NO_READING"):
        head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        os.makedirs({out!r}, exist_ok=True)
        open(os.path.join({out!r}, "reading.json"), "w").write(json.dumps({{"reading_code": {{"commit": head}}}}))
        open(os.path.join({out!r}, "reading.txt"), "w").write("OUTCOME a stop\\n")
    sys.exit(r)
sys.exit(99)
"""


def _env():
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}


def _git(repo: Path, *args) -> str:
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                           "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, text=True,
                          env=_env()).stdout.strip()


def _setup(tmp_path: Path) -> SimpleNamespace:
    repo, home = tmp_path / "repo", tmp_path / "home"
    (home / "v" / "vec-ext" / "t3" / FETCH_DIR).mkdir(parents=True)
    files = {"experiments/t3_decomp/run_k28.sh": SH.read_text(), "experiments/harness/TRIES.md": ROW,
             "EXTERNAL_DATA.md": REGISTER, "experiments/t3_decomp/k28_chip_partners.py": "# the stub uv stands in\n",
             "experiments/t3_ext/chipq_kept.yaml": "Ctnnb1: []\n", "experiments/t3_ext/CHIPQ_KEPT.sha256": "x\n"}
    repo.mkdir(parents=True)
    _git(repo, "init", "-q")
    for rel, text in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "registration")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "uv").write_text(STUB.format(py=sys.executable, block=BLOCK, out=OUT))
    (bindir / "uv").chmod(0o755)
    return SimpleNamespace(repo=repo, home=home, bin=bindir, head=_git(repo, "rev-parse", "HEAD"), tmp=tmp_path)


def _run(L, *args, trace=False, **stub) -> subprocess.CompletedProcess:
    env = {**_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}", "STUB_ARGS": str(L.tmp / "args"),
           "GIT_CEILING_DIRECTORIES": str(L.repo.parent), **{k: str(v) for k, v in stub.items()}}
    cmd = ["bash", *(["-x"] if trace else []), "experiments/t3_decomp/run_k28.sh", *args]
    return subprocess.run(cmd, cwd=L.repo, env=env, capture_output=True, text=True, timeout=120)


def _S(L, commit: str | None = None) -> Path:
    return L.home / "v" / "scratch" / "k28" / (commit or L.head)[:7] / "g0a"


def _all(L, commit: str | None = None) -> str:
    return (_S(L, commit) / "all.exit").read_text().strip()


def _args(L) -> str:
    p = L.tmp / "args"
    return p.read_text() if p.exists() else ""


def _calls(L) -> dict[str, int]:
    """How often the stub ran the tests step, the gate and the launcher's json read."""
    lines = _args(L).splitlines()
    return {"tests": sum(1 for ln in lines if " -m pytest " in ln),
            "gate": sum(1 for ln in lines if " experiments/t3_decomp/k28_chip_partners.py g0a " in ln),
            "json": sum(1 for ln in lines if " -c " in ln)}


def _commit(L, rel: str, text: str | None) -> str:
    """A further commit that writes (or, with None, stops tracking) one file; its hash."""
    if text is None:
        _git(L.repo, "rm", "-q", "--cached", rel)
    else:
        (L.repo / rel).write_text(text)
        _git(L.repo, "add", rel)
    _git(L.repo, "commit", "-q", "-m", "a change")
    return _git(L.repo, "rev-parse", "HEAD")


# ------------------------------------------------------------------------------------------------------------- text


def test_the_launcher_parses_deletes_fetches_and_uploads_nothing_and_sets_its_environment_first():
    sh = SH.read_text()
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    assert not re.search(r"(^|[\s;|&(])rm(\s|$)", sh) and "unlink" not in sh and "rmdir" not in sh
    code = "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))
    assert not re.search(r"\b(curl|wget|scp|sftp|rsync|ssh|mv)\b", code)
    assert "chipq_fetch.py" not in sh and "submit" not in sh and "format_only" not in sh
    assert ":?" not in sh  # no ${var:?message}: an apostrophe inside one stopped an earlier launcher at its first use
    assert "set -u" in sh and "finish()" in sh and "trap abort ERR TERM INT HUP EXIT" in sh
    body = sh.split("\nset -u\n", 1)[1]
    first_uv = min(body.index(x) for x in ('$PY ', '"$PY"') if x in body)
    assert body.index('export PATH="$HOME/.local/bin:$PATH"') < first_uv
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh and 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in sh
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert [int(x) for x in re.findall(r"^threads (\d+)", sh, flags=re.M)] == [8]
    assert f'CHIPATLAS=$HOME/v/vec-ext/t3/{FETCH_DIR} ' in sh and 'S="$HOME/v/scratch/k28/${COMMIT:0:7}/$STEP"' in sh
    assert f"OUT={OUT}\n" in sh and f'WORDS="{WORDS}"' in sh
    # the tests come first, then the gate; every process's output goes to /dev/null but the gate's printed block
    assert body.index("-m pytest") < body.index('$PY "$SCRIPT" g0a --commit "$COMMIT"')
    tests_cmd = body.split("-m pytest", 1)[1].split("|| rc=$?", 1)[0]
    assert tests_cmd.rstrip().endswith("> /dev/null 2>&1")
    for f in ("tests/test_chipq_partners.py", "tests/test_k28_chip_partners.py", "tests/test_k28_launcher.py",
              "tests/test_chipq.py"):
        assert f in tests_cmd
    assert all((REPO / f).is_file() for f in re.findall(r"tests/test_[a-z0-9_]+\.py", tests_cmd))
    assert '$PY "$SCRIPT" g0a --commit "$COMMIT" > "$S/g0a.out" 2> /dev/null || rc=$?' in body
    assert body.count("$PY ") == 3 and body.count("2> /dev/null") >= 3  # the json read, the tests, the gate
    assert "> \"$S/tests.log\"" not in sh and ".log" not in code


def test_usage(tmp_path):
    L = _setup(tmp_path)
    for args in ((), ("g0a",), ("g0a", L.head[:12]), ("g0a", L.head, L.head), ("g0b", L.head), ("read", L.head),
                 ("g0a", "../../../elsewhere"), ("g0a", L.head.upper()), ("g0a", "")):
        r = _run(L, *args, trace=True)
        assert r.returncode == 2 and "usage: run_k28.sh g0a COMMIT" in r.stderr and r.stdout == ""
        assert (L.home / "v/scratch/k28/usage.exit").read_text().strip() == "2"
    assert sorted(p.name for p in (L.home / "v/scratch/k28").iterdir()) == ["usage.exit"]
    assert _args(L) == ""  # nothing ran


# ------------------------------------------------------------------------------------------------------------- gates


def test_each_gate(tmp_path):
    L = _setup(tmp_path)
    r = _run(L, "g0a", "0" * 40)
    assert r.returncode == 97 and "HEAD is not" in r.stderr and _all(L, "0" * 40) == "97"
    (L.repo / "EXTERNAL_DATA.md").write_text("changed\n")
    assert _run(L, "g0a", L.head).returncode == 98 and _all(L) == "98"
    _git(L.repo, "checkout", "--", "EXTERNAL_DATA.md")

    def committed(rel, text, want, msg):
        old = (L.repo / rel).read_text()
        h = _commit(L, rel, text)
        r = _run(L, "g0a", h)
        assert r.returncode == want and msg in r.stderr and _all(L, h) == str(want), (rel, r.stderr)
        assert r.stdout == ""
        return _commit(L, rel, old)

    committed("experiments/harness/TRIES.md", "| K27 (one row) | T3 |\n| K280 | T3 |\n", 95, "no TRIES.md row K28")
    committed("EXTERNAL_DATA.md", "| X017 | other |\n", 96, "X018")
    committed("EXTERNAL_DATA.md", REGISTER.replace(WORDS, "Opened on 2026-11-23"), 96, "whole text")
    # the amendment's opening words alone, and its whole text with one word changed
    whole, changed = WORDS.split(" "), ONE_WORD.split(" ")
    assert WORDS.startswith(OPENING) and sum(x != y for x, y in zip(whole, changed, strict=True)) == 1
    committed("EXTERNAL_DATA.md", REGISTER.replace(WORDS, OPENING), 96, "whole text")
    committed("EXTERNAL_DATA.md", REGISTER.replace(WORDS, ONE_WORD), 96, "whole text")
    committed("EXTERNAL_DATA.md", REGISTER.replace("| X018 ", "| X019 "), 96, "whole text")
    committed("EXTERNAL_DATA.md", REGISTER.replace("2026-11-24", "2026-11-25"), 96, "whole text")
    for rel in ("experiments/t3_ext/chipq_kept.yaml", "experiments/t3_ext/CHIPQ_KEPT.sha256"):
        h = _commit(L, rel, None)
        r = _run(L, "g0a", h)
        assert r.returncode == 96 and f"{rel} is not tracked" in r.stderr
        _git(L.repo, "add", rel)
        _git(L.repo, "commit", "-q", "-m", "track")
    h = _git(L.repo, "rev-parse", "HEAD")
    d = L.home / "v/vec-ext/t3" / FETCH_DIR
    d.rmdir()
    r = _run(L, "g0a", h)
    assert r.returncode == 96 and "is not a real directory" in r.stderr
    (L.home / "elsewhere").mkdir()
    d.symlink_to(L.home / "elsewhere")
    r = _run(L, "g0a", h)
    assert r.returncode == 96 and "is not a real directory" in r.stderr and _all(L, h) == "96"
    assert _args(L) == ""  # no gate let a step start: no test ran and the row's script was never called
    d.unlink()
    d.mkdir()
    r = _run(L, "g0a", h)
    assert r.returncode == 0 and r.stdout == BLOCK and (_S(L, h) / "head").read_text().strip() == h


# ------------------------------------------------------------------------------------------------------------- steps


def test_the_step_maps_each_exit_and_prints_the_gates_block_only(tmp_path):
    L = _setup(tmp_path)
    r = _run(L, "g0a", L.head, STUB_TESTS_RC=1)
    assert r.returncode == 10 and _all(L) == "10" and r.stdout == "" and r.stderr == ""
    assert _calls(L) == {"tests": 1, "gate": 0, "json": 0}  # the tests come first: the gate did not run
    for rc, want in ((2, 22), (1, 20), (139, 20)):
        r = _run(L, "g0a", L.head, STUB_G0A_RC=rc)
        assert r.returncode == want and _all(L) == str(want) and r.stdout == "" and r.stderr == ""
        assert (_S(L) / "g0a.exit").read_text().strip() == str(rc)
        assert not (L.repo / OUT).exists()
    r = _run(L, "g0a", L.head)
    assert r.returncode == 0 and _all(L) == "0"
    assert r.stdout == BLOCK and r.stderr == ""  # the three lines of counts, and no other process's output
    assert (_S(L) / "g0a.out").read_text() == BLOCK
    for step in ("tests", "g0a"):
        assert (_S(L) / f"{step}.exit").read_text().strip() == "0", step
    assert (_S(L) / "head").read_text().strip() == L.head and (_S(L) / "tmp").is_dir()
    args = _args(L).splitlines()
    assert args[-1].endswith(f"experiments/t3_decomp/k28_chip_partners.py g0a --commit {L.head}")
    assert "-m pytest -q -p no:cacheprovider -o pythonpath=src --basetemp=" in args[-2]
    assert f"--basetemp={_S(L)}/tmp/pytest" in args[-2]
    assert not (L.repo / OUT).exists()  # a pass leaves no reading
    # nothing was written outside the scratch folder, the stub's own log and (never here) the results folder
    assert _git(L.repo, "status", "--porcelain") == ""
    # a traced run of the same path, as the first real launch is traced
    r = _run(L, "g0a", L.head, trace=True)
    assert r.returncode == 0 and r.stdout == BLOCK and "+ finish 0" in r.stderr


def test_a_stop_is_made_once_and_closes_the_step(tmp_path):
    L = _setup(tmp_path)
    r = _run(L, "g0a", L.head, STUB_G0A_RC=3)
    assert r.returncode == 28 and _all(L) == "28" and r.stdout == BLOCK and r.stderr == ""
    reading = L.repo / OUT / "reading.json"
    assert json.loads(reading.read_text())["reading_code"]["commit"] == L.head
    assert _calls(L) == {"tests": 1, "gate": 1, "json": 0}
    # the untracked reading leaves the tracked tree clean; a rerun at this commit finishes 28 again and runs nothing
    r = _run(L, "g0a", L.head)
    assert r.returncode == 28 and r.stdout == "" and _all(L) == "28"
    assert _calls(L) == {"tests": 1, "gate": 1, "json": 1}
    # a reading made at another commit, naming none, unreadable, or half there: 44, and nothing runs
    for body in (json.dumps({"reading_code": {"commit": "f" * 40}}), json.dumps({"outcome": "x"}), "not json", None):
        if body is None:
            reading.rename(reading.with_name("reading.json.moved"))
        else:
            reading.write_text(body)
        r = _run(L, "g0a", L.head)
        assert r.returncode == 44 and _all(L) == "44" and r.stdout == ""
    assert _calls(L) == {"tests": 1, "gate": 1, "json": 5}
    # exit 3 with no reading left behind is not a stop
    L2 = _setup(tmp_path / "two")
    r = _run(L2, "g0a", L2.head, STUB_G0A_RC=3, STUB_NO_READING=1)
    assert r.returncode == 41 and _all(L2) == "41" and r.stdout == BLOCK
