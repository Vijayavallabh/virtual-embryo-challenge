"""Tests for experiments/t3_decomp/run_k27.sh (row K27's launcher, two invocations): every gate, every step's exit
mapping and the once-only stop readings, traced in a temporary git repository under a temporary HOME with a stub `uv`
whose exit status per step the test sets. No Python step of the try runs; no data of any kind is read."""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]
SH = REPO / "experiments" / "t3_decomp" / "run_k27.sh"
SCRIPT = REPO / "experiments" / "t3_decomp" / "k27_chip_query.py"
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES", "GIT_CEILING_DIRECTORIES")
ROW = "| K27 (one row) | T3 | `experiments/t3_decomp/k27_chip_query.py`, cloned from k24 |\n"
REGISTER = "| X018 | ChIP-Atlas, a generic database of public ChIP-seq | ... |\n"
NAMES = ["chipq_Gata4", *(f"chipq_rand_seed{s}" for s in range(10))]

# The stub `uv`: a python -c one-liner runs for real (the launcher's json reads); every step exits with the code the
# test sets and writes what the real step would leave for the launcher to read (p2.json, pregate.json, guard.json,
# reading.json), so the launcher's own logic is what runs.
STUB = """#!{py}
import json, os, subprocess, sys
a = sys.argv[1:]
log = os.environ.get("STUB_ARGS")
if log:
    open(log, "a").write(" ".join(a) + "\\n")
if "-c" in a:
    os.execv(sys.executable, [sys.executable, *a[4:]])
rc = lambda k: int(os.environ.get(k, "0"))
opt = lambda n: a[a.index(n) + 1] if n in a else None
if "pytest" in a:
    sys.exit(rc("STUB_TESTS_RC"))
script = next((x for x in a if x.endswith(".py")), "")
step = a[a.index(script) + 1] if script else ""
if script.endswith("chipq_fetch.py"):
    sys.exit(rc("STUB_VERIFY_RC"))
if script.endswith("chipq.py"):
    print("4 22")
    sys.exit(rc("STUB_P1_RC"))
def put(p, body):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(json.dumps(body))
if step == "p2":
    r = rc("STUB_P2_RC")
    if r in (0, 3):
        put(os.path.join(opt("--out"), "p2.json"), {{"passes": r == 0}})
    sys.exit(r)
if step == "pins":
    print("# lines")
    sys.exit(rc("STUB_PINS_RC"))
if step == "check":
    sys.exit(rc("STUB_CHECK_RC"))
if step == "pregate":
    r = rc("STUB_PREGATE_RC")
    if r in (0, 5):
        put(os.path.join(opt("--out"), "pregate.json"), {{"pregate": {{"dropped": "corr" if r == 5 else None}}}})
    sys.exit(r)
if step == "run":
    sys.exit(rc("STUB_GPART_RC" if opt("--part") == "guard" else "STUB_SPART_RC"))
if step == "guard":
    r = rc("STUB_GUARD_RC")
    if r in (0, 6):
        put(os.path.join(opt("--runs"), "guard.json"), {{"passes": r == 0}})
    sys.exit(r)
if step == "read":
    open(os.environ["STUB_READS"], "a").write("read\\n")
    r = rc("STUB_READ_RC")
    if r == 0:
        head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        put(os.path.join(opt("--out"), "reading.json"), {{"reading_code": {{"commit": head}}}})
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


def _filled_script(text_: str) -> tuple[str, str]:
    """The script's text with the 11 pins and the listing pin filled, and the matching listing."""
    tabs = {n: hashlib.sha256(n.encode()).hexdigest() for n in NAMES}
    listing = "".join(f"{tabs[n]}  {n}\n" for n in NAMES)
    lsha = hashlib.sha256(listing.encode()).hexdigest()
    t = text_.replace("K27_LISTING_SHA256: str | None = None", f'K27_LISTING_SHA256: str | None = "{lsha}"')
    for n in NAMES:
        t = t.replace(f'    "{n}": None,', f'    "{n}": "{tabs[n]}",')
    return t, listing


def _setup(tmp_path: Path) -> SimpleNamespace:
    repo, home = tmp_path / "repo", tmp_path / "home"
    for d in ("v/vec-ext/t3/tables_k27", "v/vec-ext/t3/tables_two_study", "v/vec-ext/t3/chipatlas"):
        (home / d).mkdir(parents=True)
    files = {"experiments/t3_decomp/run_k27.sh": SH.read_text(), "experiments/harness/TRIES.md": ROW,
             "EXTERNAL_DATA.md": REGISTER, "experiments/t3_ext/TABLES_K27.sha256": "",
             "experiments/t3_decomp/k27_chip_query.py": SCRIPT.read_text(),
             "experiments/t3_ext/chipq_kept.yaml": "Ctnnb1: []\n", "experiments/t3_ext/CHIPQ_KEPT.sha256": "x\n",
             "experiments/t3_decomp/results/k23_outside_exponents/pregate.json": "{}",
             "experiments/t3_decomp/results/k23_outside_exponents/reading.json": "{}"}
    repo.mkdir()
    _git(repo, "init", "-q")
    for rel, t in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t)
    (repo / ".gitignore").write_text("runs/\nexperiments/t3_decomp/results/k27_chip_query/\n")
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "registration")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "uv").write_text(STUB.format(py=sys.executable))
    (bindir / "uv").chmod(0o755)
    return SimpleNamespace(repo=repo, home=home, bin=bindir, reg=_git(repo, "rev-parse", "HEAD"), tmp=tmp_path)


def _second(L) -> str:
    """The row's second commit: the listing's 11 lines and the pins, nothing else."""
    t, listing = _filled_script(SCRIPT.read_text())
    (L.repo / "experiments/t3_decomp/k27_chip_query.py").write_text(t)
    (L.repo / "experiments/t3_ext/TABLES_K27.sha256").write_text(listing)
    _git(L.repo, "commit", "-q", "-am", "second")
    return _git(L.repo, "rev-parse", "HEAD")


def _run(L, *args, trace=False, **stub) -> tuple[int, str]:
    env = {**_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}", "STUB_ARGS": str(L.tmp / "args"),
           "STUB_READS": str(L.tmp / "reads"), "GIT_CEILING_DIRECTORIES": str(L.repo.parent),
           **{k: str(v) for k, v in stub.items()}}
    cmd = ["bash", *(["-x"] if trace else []), "experiments/t3_decomp/run_k27.sh", *args]
    r = subprocess.run(cmd, cwd=L.repo, env=env, capture_output=True, text=True, timeout=120)
    return r.returncode, r.stderr


def _S(L, mode) -> Path:
    return L.home / "v" / "scratch" / f"k27_{L.reg[:7]}" / mode


def _all(L, mode) -> str:
    return (_S(L, mode) / "all.exit").read_text().strip()


def _reads(L) -> int:
    p = L.tmp / "reads"
    return len(p.read_text().splitlines()) if p.exists() else 0


def _db_ok(L) -> None:
    code, err = _run(L, "db", L.reg)
    assert code == 0 and _all(L, "db") == "0", err


# ------------------------------------------------------------------------------------------------------------- text


def test_the_launcher_parses_removes_nothing_and_sets_its_environment_first():
    sh = SH.read_text()
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    assert not re.search(r"(^|[\s;|&(])rm(\s|$)", sh) and "unlink" not in sh and "rmdir" not in sh
    assert "set -u" in sh and "finish()" in sh and "trap abort ERR TERM INT HUP EXIT" in sh
    body = sh.split("\nset -u\n", 1)[1]
    first_uv = min(body.index(x) for x in ('$PY ', '"$PY"') if x in body)
    assert body.index('export PATH="$HOME/.local/bin:$PATH"') < first_uv
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert max(int(x) for x in re.findall(r"threads (\d+)", sh)) == 16
    tests_line = sh.split("-m pytest", 1)[1].split("> \"$S/tests.log\"", 1)[0]
    for f in ("tests/test_k27_chip_query.py", "tests/test_k27_launcher.py", "tests/test_chipq.py",
              "tests/test_k24_wider_tables.py", "tests/test_k23_outside_exponents.py"):
        assert f in tests_line
    assert all((REPO / f).is_file() for f in re.findall(r"tests/test_[a-z0-9_]+\.py", tests_line))


def test_usage(tmp_path):
    L = _setup(tmp_path)
    for args in ((), ("db",), ("run", L.reg), ("bogus", L.reg), ("db", L.reg, L.reg)):
        code, err = _run(L, *args, trace=True)
        assert code == 2 and "usage: run_k27.sh" in err
        assert (L.home / "v/scratch/k27/usage.exit").read_text().strip() == "2"


# ------------------------------------------------------------------------------------------------------------- db


def test_each_gate_of_the_database_part(tmp_path):
    L = _setup(tmp_path)
    code, err = _run(L, "db", "0" * 40)
    assert code == 97 and "HEAD is not" in err
    assert (L.home / "v/scratch/k27_0000000/db/all.exit").read_text().strip() == "97"
    (L.repo / "EXTERNAL_DATA.md").write_text("changed\n")
    assert _run(L, "db", L.reg)[0] == 98 and _all(L, "db") == "98"
    _git(L.repo, "checkout", "--", "EXTERNAL_DATA.md")

    def committed(rel, text_, want, msg):
        p = L.repo / rel
        old = p.read_text()
        p.write_text(text_)
        _git(L.repo, "commit", "-q", "-am", "x")
        h = _git(L.repo, "rev-parse", "HEAD")
        code, err = _run(L, "db", h)
        assert code == want and msg in err, (rel, err)
        p.write_text(old)
        _git(L.repo, "commit", "-q", "-am", "back")
        return _git(L.repo, "rev-parse", "HEAD")

    committed("experiments/harness/TRIES.md", "| K27 | another script |\n", 95, "no TRIES.md row")
    committed("experiments/t3_ext/TABLES_K27.sha256", "x\n", 22, "not empty")
    filled, _ = _filled_script(SCRIPT.read_text())
    committed("experiments/t3_decomp/k27_chip_query.py", filled, 22, "not empty")
    half = SCRIPT.read_text().replace('    "chipq_rand_seed4": None,', f'    "chipq_rand_seed4": "{"a" * 64}",')
    committed("experiments/t3_decomp/k27_chip_query.py", half, 22, "not empty")
    committed("EXTERNAL_DATA.md", "| X017 | other |\n", 96, "X018")
    _git(L.repo, "rm", "-q", "--cached", "experiments/t3_ext/chipq_kept.yaml")
    _git(L.repo, "commit", "-q", "-m", "untrack")
    h = _git(L.repo, "rev-parse", "HEAD")
    code, err = _run(L, "db", h)
    assert code == 96 and "chipq_kept.yaml is not tracked" in err
    _git(L.repo, "add", "experiments/t3_ext/chipq_kept.yaml")
    _git(L.repo, "commit", "-q", "-m", "track")
    h = _git(L.repo, "rev-parse", "HEAD")
    d = L.home / "v/vec-ext/t3/tables_k27"
    d.rmdir()
    (L.home / "elsewhere").mkdir()
    d.symlink_to(L.home / "elsewhere")
    code, err = _run(L, "db", h)
    assert code == 96 and "tables_k27 is not a real directory" in err


def test_the_database_part_maps_each_step(tmp_path):
    L = _setup(tmp_path)
    assert _run(L, "db", L.reg, STUB_TESTS_RC=1)[0] == 10 and _all(L, "db") == "10"
    assert _run(L, "db", L.reg, STUB_VERIFY_RC=2)[0] == 22
    assert _run(L, "db", L.reg, STUB_P1_RC=2)[0] == 22
    assert _run(L, "db", L.reg, STUB_P1_RC=1)[0] == 20
    assert _run(L, "db", L.reg, STUB_P2_RC=2)[0] == 22 and not (L.repo / f"runs/k27_{L.reg[:7]}/p2.json").exists()
    assert _run(L, "db", L.reg, STUB_P2_RC=1)[0] == 20
    assert _run(L, "db", L.reg, STUB_PINS_RC=2)[0] == 22
    S = _S(L, "db")
    R = L.repo / f"runs/k27_{L.reg[:7]}"
    assert (R / "p1.out").read_text() == "4 22\n" and (R / "p1.exit").read_text().strip() == "0"
    _db_ok(L)
    assert (S / "second_commit.txt").read_text() == "# lines\n"
    for step in ("tests", "verify", "p1", "pins"):
        assert (S / f"{step}.exit").read_text().strip() == "0", step
    args = (L.tmp / "args").read_text()
    assert f"p1 --data {L.repo}/data/raw --run-dir" in args and f"p2 --out runs/k27_{L.reg[:7]}" in args
    assert " pregate " not in args and " guard " not in args and _reads(L) == 0   # nothing counted runs


def test_a_p1_or_p2_stop_is_read_once(tmp_path):
    L = _setup(tmp_path)
    assert _run(L, "db", L.reg, STUB_P1_RC=3, STUB_READ_RC=1)[0] == 41
    assert _run(L, "db", L.reg, STUB_P1_RC=3, STUB_READ_RC=3)[0] == 43
    assert _run(L, "db", L.reg, STUB_P1_RC=3)[0] == 28 and _reads(L) == 3 and _all(L, "db") == "28"
    assert _run(L, "db", L.reg, STUB_P1_RC=3)[0] == 28 and _reads(L) == 3   # this commit's reading: not read again
    out = L.repo / "experiments/t3_decomp/results/k27_chip_query/reading.json"
    out.write_text(json.dumps({"reading_code": {"commit": "f" * 40}}))
    assert _run(L, "db", L.reg, STUB_P1_RC=3)[0] == 44 and _reads(L) == 3
    L2 = _setup(tmp_path / "two")
    assert _run(L2, "db", L2.reg, STUB_P2_RC=3)[0] == 29 and _reads(L2) == 1
    # a rerun finds P2's record and its failure: 29 again, no second reading, no second P2
    assert _run(L2, "db", L2.reg, STUB_P2_RC=0)[0] == 29 and _reads(L2) == 1
    assert (L2.tmp / "args").read_text().count(" p2 --out ") == 1


# ------------------------------------------------------------------------------------------------------------- run


def test_each_gate_of_the_counted_part(tmp_path):
    L = _setup(tmp_path)
    # at the registration commit (pins empty): HEAD is not a second commit of REG
    code, err = _run(L, "run", L.reg, L.reg)
    assert code == 97 and "strict ancestor" in err and _all(L, "run") == "97"
    _db_ok(L)
    (L.repo / "x.txt").write_text("x")
    _git(L.repo, "add", "x.txt")
    _git(L.repo, "commit", "-q", "-m", "not the pins")
    h = _git(L.repo, "rev-parse", "HEAD")
    code, err = _run(L, "run", h, L.reg)
    assert code == 22 and "not filled" in err
    second = _second(L)
    code, err = _run(L, "run", second, L.reg)
    assert code == 22 and "more than the listing and the pins" in err   # x.txt changed between the two commits
    L = _setup(tmp_path / "b")
    second = _second(L)
    code, err = _run(L, "run", second, L.reg)
    assert code == 96 and "has not finished 0" in err
    code, err = _run(L, "run", second, "1" * 40)
    assert code == 97
    code, err = _run(L, "run", L.reg, second)
    assert code == 97
    # the database part at a commit whose pins and listing are filled
    code, err = _run(L, "db", second)
    assert code == 22 and "not empty" in err


def _ready(tmp_path) -> tuple[SimpleNamespace, str]:
    L = _setup(tmp_path)
    second = _second(L)
    S = L.home / "v" / "scratch" / f"k27_{L.reg[:7]}" / "db"
    S.mkdir(parents=True)
    (S / "all.exit").write_text("0\n")
    R = L.repo / f"runs/k27_{L.reg[:7]}"
    R.mkdir(parents=True)
    (R / "p2.json").write_text(json.dumps({"passes": True}))
    return L, second


def test_the_counted_part_maps_each_step(tmp_path):
    L, h = _ready(tmp_path)
    assert _run(L, "run", h, L.reg, STUB_TESTS_RC=1)[0] == 10
    assert _run(L, "run", h, L.reg, STUB_VERIFY_RC=2)[0] == 22
    assert _run(L, "run", h, L.reg, STUB_CHECK_RC=2)[0] == 22
    assert _run(L, "run", h, L.reg, STUB_CHECK_RC=1)[0] == 20
    assert _run(L, "run", h, L.reg, STUB_PREGATE_RC=2)[0] == 22
    assert _run(L, "run", h, L.reg, STUB_PREGATE_RC=1)[0] == 20
    assert _run(L, "run", h, L.reg, STUB_GPART_RC=1)[0] == 30
    R = L.repo / f"runs/k27_{L.reg[:7]}"
    gparts = sorted(p.name for p in (R / "exit").iterdir() if p.name.startswith("guard_") and "." not in p.name)
    assert gparts == sorted(f"guard_{p}" for p in ("plain_g1_0-4", "plain_g2_0-4", "plain_g1_5-9", "plain_g2_5-9",
                                                   "dm"))
    assert _run(L, "run", h, L.reg, STUB_GUARD_RC=2)[0] == 22
    assert _run(L, "run", h, L.reg, STUB_GUARD_RC=3)[0] == 30
    assert _run(L, "run", h, L.reg, STUB_SPART_RC=1)[0] == 30   # the guard passed and is recorded; a part fails
    for rc, want in ((3, 43), (4, 44), (1, 41), (0, 0)):
        assert _run(L, "run", h, L.reg, STUB_READ_RC=rc)[0] == want and _all(L, "run") == str(want)
    args = (L.tmp / "args").read_text()
    assert args.count(" pregate --out ") == 3   # refused (2), crashed (1), then made once: never again
    n_guard = args.count(" run --part guard ")
    assert n_guard % 5 == 0 and n_guard >= 10 and args.count(" run --part score --world dm ") >= 2
    assert all(f"--data {L.repo}/data/raw" in ln for ln in args.splitlines() if " run --part " in ln)
    for step in ("tests", "verify", "check", "guard", "read"):
        assert (_S(L, "run") / f"{step}.exit").exists(), step


def test_a_drop_at_p3_is_read_once_and_counts(tmp_path):
    L, h = _ready(tmp_path)
    assert _run(L, "run", h, L.reg, STUB_PREGATE_RC=5)[0] == 25 and _reads(L) == 1
    assert _run(L, "run", h, L.reg)[0] == 25 and _reads(L) == 1   # pregate.json records the drop: no second read
    args = (L.tmp / "args").read_text()
    assert " run --part " not in args


def test_a_failed_guard_is_read_once_and_starts_no_scoring_part(tmp_path):
    L, h = _ready(tmp_path)
    assert _run(L, "run", h, L.reg, STUB_GUARD_RC=6)[0] == 26 and _reads(L) == 1 and _all(L, "run") == "26"
    assert _run(L, "run", h, L.reg)[0] == 26 and _reads(L) == 1
    assert " run --part score " not in (L.tmp / "args").read_text()
    out = L.repo / "experiments/t3_decomp/results/k27_chip_query/reading.json"
    out.write_text(json.dumps({"reading_code": {"commit": "f" * 40}}))
    assert _run(L, "run", h, L.reg)[0] == 44
    out.write_text("not json")
    assert _run(L, "run", h, L.reg)[0] == 44
    L2, h2 = _ready(tmp_path / "two")
    assert _run(L2, "run", h2, L2.reg, STUB_GUARD_RC=6, STUB_READ_RC=1)[0] == 41
    assert _run(L2, "run", h2, L2.reg, STUB_READ_RC=3)[0] == 43


def test_the_launchers_logs_are_never_printed(tmp_path):
    L, h = _ready(tmp_path)
    r = subprocess.run(["bash", "experiments/t3_decomp/run_k27.sh", "run", h, L.reg], cwd=L.repo,
                       env={**_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}",
                            "STUB_READS": str(L.tmp / "reads"), "GIT_CEILING_DIRECTORIES": str(L.repo.parent)},
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0 and r.stdout == "" and r.stderr == ""
