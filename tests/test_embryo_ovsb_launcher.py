"""The test-board build launcher of row R48's file T2embryo_s2x_ovsb.h5ad: experiments/builds/embryo_ovsb.sh.

Its gate builds the file only when row R48 is the captain's order naming the option and the ruling (the builder's own
check), and step 0's reading (row E56, experiments/coords/ovsb_step0.py) is committed with the outcome "Step 0: PASS",
made on a clean tree at an ancestor of HEAD that holds row E56; O's, SO's and O_hi's commands must give their committed
sha256 (read by code from their committed build records) before the file is built; the file is built twice and
compared; the identity step prints booleans only. These tests hold the launcher's strings and step order, run every
gate refusal for real in throwaway git repositories under a temporary home with a stand-in uv that fails (so nothing is
built and no competition file is read), walk every step once with a stand-in uv whose builds copy fixtures, and run
its gate's two checks and its identity step on synthetic ledgers, readings and files. The launcher is never run here
with a real commit. Built through experiments/builds/build_test.py."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from conftest import git, git_env

REPO = Path(__file__).resolve().parents[1]
REL = "experiments/builds/embryo_ovsb.sh"
LAUNCHER = REPO / REL
FLOOR = "--method fgw_shape0 --positions floor_exact"
OV = f"{FLOOR} --expression mixture --shape-interp overlay --shape-interp-test-ok"
SO = f"{FLOOR} --expression series --shape-interp overlay --shape-interp-test-ok"
OHI = f"{FLOOR} --expression render_hi --shape-interp overlay --shape-interp-test-ok"
W = f"{FLOOR} --expression series_back --shape-interp overlay --shape-interp-test-ok --captain-ruling R48"
SHAS = {"O": "1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea",
        "SO": "77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165",
        "OH": "7e4d1eaafc40080d9d9c461d4b9396f0cee6353686538e79fb92ca31d1d2deeb"}
RECORDS = {"O": ("experiments/builds/results/embov_3b813b5.SHA256SUMS", "OVa/T2embryo_s2x_ov.h5ad"),
           "SO": ("experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS", "SOa/T2embryo_s2xs_ov.h5ad"),
           "OH": ("experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS", "OHa/T2embryo_s2x_ovhi.h5ad")}
READING_DIR = "experiments/coords/results/ovsb_step0"
E56 = "| E56 (step 0 of the file of row R48; synthetic) | T2 embryo |"
R48 = ("| R48 (the captain's order of 2026-11-24, ruling, no run: \"Build T2embryo_s2x_ovsb.h5ad as named\"; the "
       "embryo builder is to take `--expression series_back` on O's command only with `--captain-ruling R48`; not a "
       "try) | T2 embryo | synthetic |")


# ------------------------------------------------------------------------------------------------ the launcher's text
def test_the_commands_are_the_committed_ones_and_the_files():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^([A-Z_]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"OV": OV, "SO": SO, "OH": OHI, "W": W}
    # O's command is embryo_overlay.sh's; SO's and O_hi's are embryo_so_ohi.sh's; W's is the builder's named command
    assert f'OV_ARGS="{OV}"' in (REPO / "experiments/builds/embryo_overlay.sh").read_text()
    sohi = (REPO / "experiments/builds/embryo_so_ohi.sh").read_text()
    assert f'SO_ARGS="{SO}"' in sohi and f'OH_ARGS="{OHI}"' in sohi
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    assert build_test.EMBRYO_OVSB_COMMAND == f"embryo {W}"
    for k, sha in SHAS.items():
        assert f"{k}_SHA={sha}" in sh and f'[ "$got" = "${k}_SHA" ]' in sh
        sums, entry = RECORDS[k]
        assert f"{k}_SUMS={sums}" in sh and f"{k}_ENTRY={entry}" in sh
    assert "BOARD=T2:embryo:test_interp" in sh and f"READING_DIR={READING_DIR}" in sh
    assert "build O T2embryo_s2x_ov $OV_ARGS" in sh and "build SO T2embryo_s2xs_ov $SO_ARGS" in sh
    assert "build OH T2embryo_s2x_ovhi $OH_ARGS" in sh and 'build "$STEP" T2embryo_s2x_ovsb $W_ARGS' in sh
    assert 'cmp -s "$O/Wa/T2embryo_s2x_ovsb.h5ad" "$O/Wb/T2embryo_s2x_ovsb.h5ad" || finish 31' in sh
    builds = [ln for ln in sh.splitlines() if "build_test.py embryo" in ln]
    assert builds and not any("--board" in ln for ln in builds)  # no dry run: the option refuses one


def test_the_committed_records_hold_the_three_sha256():
    for k, (sums, entry) in RECORDS.items():
        rows = dict(reversed(line.split()) for line in (REPO / sums).read_text().splitlines() if line.strip())
        assert rows[entry] == SHAS[k]


def test_the_steps_run_in_order_each_writing_running_first():
    sh = LAUNCHER.read_text()
    order = []
    for m in re.finditer(r"^STEP=(\w+)$|^for STEP in ([\w ]+); do$", sh, flags=re.M):
        order += [m.group(1)] if m.group(1) else m.group(2).split()
    assert order == ["gate", "tests", "O", "SO", "OH", "Wa", "Wb", "identity", "done"]
    assert sh.count("echo running > \"$S/$STEP.exit\"") == 8  # gate, tests, O, SO, OH, the loop, identity, done
    assert sh.index('echo running > "$S/$STEP.exit"') < sh.index("head=$(git rev-parse HEAD)")
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh
    i8 = sh.index("OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8")
    i16 = sh.index("OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16")
    assert i8 < sh.index("STEP=tests") < i16 < sh.index("STEP=O\n")
    so_ohi = (REPO / "experiments/builds/embryo_so_ohi.sh").read_text()  # the existing launcher's pins
    assert so_ohi.index("OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16") < so_ohi.index("STEP=D\n")
    assert "S=$HOME/v/scratch/embovsb/$RUN" in sh and "O=runs/embovsb_${head:0:7}_$RUN" in sh
    assert '> "$S/identity.out" 2> /dev/null || finish 40' in sh
    assert "sha256sum O/*.h5ad SO/*.h5ad OH/*.h5ad Wa/*.h5ad Wb/*.h5ad) > \"$S/SHA256SUMS\"" in sh


def test_every_builder_and_check_call_goes_to_dev_null_and_nothing_is_deleted():
    text = LAUNCHER.read_text()
    joined = text.replace("\\\n", " ")
    for line in joined.splitlines():
        if re.search(r"build_test\.py embryo|vec check \"", line):
            assert "> /dev/null 2>&1" in line, line
    code = [line.split("#")[0] for line in text.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", text):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg
    assert "upload" not in " ".join(code)
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", text)
    assert {"tests/test_embryo_ovsb_launcher.py", "tests/test_series_back_build.py",
            "tests/test_ovsb_step0.py"} <= set(listed) and all((REPO / f).exists() for f in listed)
    assert "tests/test_embryo_ovsb_launcher.py" in (REPO / "experiments/external/n54_step.sh").read_text()


# ------------------------------------------------------------------------------------------- the gate's two checks
def _snippet(start: str) -> str:
    m = re.search(re.escape(start) + r"'([^']+)'", LAUNCHER.read_text())
    assert m, f"the check after {start!r} is not where the test expects it"
    return m.group(1)


def _python(code: str, *paths) -> subprocess.CompletedProcess:
    return subprocess.run(["python3", "-I", "-B", "-c", code, *map(str, paths)], capture_output=True, text=True)


def test_the_ledger_check_wants_exactly_one_row_r48_naming_the_option_and_the_ruling(tmp_path):
    code = _snippet("python3 -I -B -c ")
    cases = [(E56 + "\n" + R48, 0), (R48, 0), (E56, 1), (R48 + "\n" + R48, 1),
             (R48.replace("`--expression series_back`", "--expression series_back"), 1),
             (R48.replace("`--captain-ruling R48`", "--captain-ruling R48"), 1),
             (R48.replace("ruling, no run", "ruling"), 1), (R48.replace("the captain's order", "the order"), 1),
             ("| R48 (the captain's order, ruling, no run) | `--expression series_back` `--captain-ruling R48` |", 1),
             (R48.replace("| R48 ", "| R48a "), 1)]
    for k, (rows, want) in enumerate(cases):
        p = tmp_path / f"TRIES{k}.md"
        p.write_text("| ID | what |\n|---|---|\n" + rows + "\n")
        assert _python(code, p).returncode == want, rows
    # the builder's own check agrees on every case
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    for k, (_, want) in enumerate(cases):
        got = build_test._captain_order_names("R48", build_test.EMBRYO_OVSB_OPTION, tmp_path / f"TRIES{k}.md")
        assert got == (want == 0), k


def _reading(commit: str, outcome: str = "Step 0: PASS", json_outcome: str | None = None, clean=True) -> dict:
    return {"txt": f"{outcome}\n(i) reported only: synthetic\n",
            "json": {"outcome": outcome if json_outcome is None else json_outcome, "row": "E56",
                     "provenance": {"commit": commit, "clean": clean}}}


def test_the_reading_check_wants_a_pass_on_a_clean_tree_at_a_named_commit(tmp_path):
    code = _snippet("RUNC=$(python3 -I -B -c ")
    run = "a" * 40
    cases = [(_reading(run), run), (_reading(run, "Step 0: STOP (ii)"), None),
             (_reading(run, "Step 0: STOP (iii)"), None), (_reading(run, "Step 0: STOP (ii, iii)"), None),
             (_reading(run, json_outcome="Step 0: STOP (ii)"), None), (_reading(run, clean=False), None),
             (_reading(run, clean="True"), None), (_reading(run[:7]), None), (_reading("x" * 40), None),
             (_reading(run, "Step 0: PASS (ii)"), None), (_reading(run, " Step 0: PASS"), None)]
    for k, (r, want) in enumerate(cases):
        txt, js = tmp_path / f"r{k}.txt", tmp_path / f"r{k}.json"
        txt.write_text(r["txt"])
        js.write_text(json.dumps(r["json"]))
        got = _python(code, txt, js)
        if want is None:
            assert got.returncode == 1, r
        else:
            assert got.returncode == 0 and got.stdout.strip() == want, (r, got.stderr)
    js = tmp_path / "list.json"
    js.write_text("[]")
    assert _python(code, txt, js).returncode == 1  # not a reading at all


# ------------------------------------------------------------------------------------------------ the gate, for real
def _checkout(tmp_path: Path, at_run=(E56, R48), at_head=None, reading=lambda c: _reading(c), track=("txt", "json"),
              records=True) -> tuple[Path, str]:
    """A throwaway repository: the launcher, the rules file and its hash, the committed build records and a ledger with
    `at_run`'s rows, committed (the step 0 commit); then the ledger with `at_head`'s rows (default unchanged) and step
    0's reading made at that commit (unless `reading` is None), the files in `track` committed."""
    repo = tmp_path / "checkout"
    for rel in (REL, "docs/S2_RULES.md", "docs/S2_RULES.sha256", *sorted({r[0] for r in RECORDS.values()})):
        if not records and rel.endswith(".SHA256SUMS"):
            continue
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, repo / rel)
    tries = repo / "experiments/harness/TRIES.md"
    tries.parent.mkdir(parents=True)
    tries.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in at_run))
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "rows")
    run = git(repo, "rev-parse", "HEAD")
    if at_head is not None:
        tries.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in at_head))
    if reading is not None:
        r = reading(run)
        d = repo / READING_DIR
        d.mkdir(parents=True)
        (d / "reading.txt").write_text(r["txt"])
        (d / "reading.json").write_text(json.dumps(r["json"]))
        for kind in track:
            git(repo, "add", f"{READING_DIR}/reading.{kind}")
    git(repo, "add", "experiments/harness/TRIES.md")
    git(repo, "commit", "-qm", "step 0 read", "--allow-empty")
    return repo, run


def _home(tmp_path: Path, name: str = "home") -> Path:
    """A home whose ~/.local/bin holds a `uv` that fails (the launcher puts that folder first on PATH)."""
    home = tmp_path / name
    fake = home / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text("#!/bin/sh\necho \"$*\" >> \"$HOME/uvcalls\"\nexit 1\n")
    (fake / "uv").chmod(0o755)
    return home


def _run(repo: Path, home: Path, *argv: str, trace: bool = False) -> subprocess.CompletedProcess:
    env = {**git_env(), "HOME": str(home)}
    return subprocess.run(["bash", *(["-x"] if trace else []), REL, *argv], cwd=repo, env=env, capture_output=True,
                          text=True, timeout=60)


def _scratch(home: Path, run: str) -> Path:
    return home / "v" / "scratch" / "embovsb" / run


def _refused(repo, home, run, code, words):
    out = _run(repo, home, git(repo, "rev-parse", "HEAD"), run)
    assert out.returncode == code and words in out.stderr, (out.returncode, out.stderr)
    s = _scratch(home, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]
    assert not (repo / "runs").exists() and not (home / "uvcalls").exists()


def test_a_wrong_commit_is_refused_and_its_trace_writes_the_exit_record_before_any_python(tmp_path):
    repo, _ = _checkout(tmp_path)
    home = _home(tmp_path)
    out = _run(repo, home, "0" * 40, "t1", trace=True)
    assert out.returncode == 96 and "HEAD is not" in out.stderr and "finish 96" in out.stderr
    assert not re.search(r"^\++ (uv|python3?) ", out.stderr, flags=re.M)  # nothing ran but git and the shell
    assert (_scratch(home, "t1") / "gate.exit").read_text().strip() == "96"


def test_bad_arguments_and_a_reused_run_folder_are_refused(tmp_path):
    repo, _ = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    assert _run(repo, home, head, "a/b").returncode == 90
    assert _run(repo, home, head[:12], "r0").returncode == 90
    assert _run(repo, home, head.upper(), "r0").returncode == 90
    assert _run(repo, home, head).returncode != 0 and not (home / "v").exists()
    assert _run(repo, home, "0" * 40, "r0").returncode == 96
    again = _run(repo, home, head, "r0")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert (_scratch(home, "r0") / "gate.exit").read_text().strip() == "96"


def test_an_unclean_tree_or_a_changed_rules_file_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a")
    (repo / "experiments/harness/TRIES.md").write_text(E56 + "\n" + R48 + "\n| E57 changed |\n")
    _refused(repo, _home(tmp_path, "h1"), "u1", 98, "tracked files differ")
    repo, _ = _checkout(tmp_path / "b")
    with open(repo / "docs/S2_RULES.md", "a") as f:
        f.write("\nchanged\n")
    git(repo, "commit", "-qam", "rules changed")
    _refused(repo, _home(tmp_path, "h2"), "u2", 97, "")


def test_without_row_r48_as_the_builder_needs_it_the_gate_refuses(tmp_path):
    cases = [(E56,), (E56, R48, R48), (E56, R48.replace("`--expression series_back`", "series_back")),
             (E56, R48.replace("ruling, no run", "ruling"))]
    for k, rows in enumerate(cases):
        repo, _ = _checkout(tmp_path / f"c{k}", at_run=rows)
        _refused(repo, _home(tmp_path, f"h{k}"), f"r{k}", 95, "row R48")


def test_without_a_committed_reading_the_gate_refuses(tmp_path):
    repo, _ = _checkout(tmp_path / "a", reading=None)
    _refused(repo, _home(tmp_path, "h1"), "n1", 94, "no committed reading")
    repo, _ = _checkout(tmp_path / "b", track=("json",))  # the text is there but not tracked
    _refused(repo, _home(tmp_path, "h2"), "n2", 94, "reading.txt must be tracked")
    repo, _ = _checkout(tmp_path / "c", track=("txt",))
    _refused(repo, _home(tmp_path, "h3"), "n3", 94, "reading.json must be tracked")


def test_a_reading_that_is_not_a_pass_is_refused(tmp_path):
    for k, make in enumerate((lambda c: _reading(c, "Step 0: STOP (ii)"), lambda c: _reading(c, "Step 0: STOP (iii)"),
                              lambda c: _reading(c, "Step 0: STOP (ii, iii)"),
                              lambda c: _reading(c, json_outcome="Step 0: STOP (ii)"),
                              lambda c: _reading(c, clean=False), lambda c: _reading(c[:7]))):
        repo, _ = _checkout(tmp_path / f"c{k}", reading=make)
        _refused(repo, _home(tmp_path, f"h{k}"), f"p{k}", 92, "step 0 did not pass")


def test_a_foreign_step0_commit_or_one_without_row_e56_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a", reading=lambda c: _reading("f" * 40))
    _refused(repo, _home(tmp_path, "h1"), "c1", 88, "not an ancestor")
    repo, _ = _checkout(tmp_path / "b", at_run=(R48,), at_head=(E56, R48))  # row E56 came after step 0 ran
    _refused(repo, _home(tmp_path, "h2"), "c2", 86, "row E56 is not in TRIES.md at step 0")


def test_a_missing_or_different_build_record_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a", records=False)
    _refused(repo, _home(tmp_path, "h1"), "s1", 85, "does not give O_SHA")
    repo, _ = _checkout(tmp_path / "b")
    sums = repo / RECORDS["OH"][0]
    sums.write_text(sums.read_text().replace(SHAS["OH"], "0" * 64))
    git(repo, "commit", "-qam", "record changed")
    _refused(repo, _home(tmp_path, "h2"), "s2", 85, "does not give OH_SHA")


def test_an_existing_output_folder_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    (repo / "runs" / f"embovsb_{head[:7]}_o1").mkdir(parents=True)
    out = _run(repo, _home(tmp_path), head, "o1")
    assert out.returncode == 93 and "output folder" in out.stderr


def _exit_files(s: Path) -> dict:
    return {p.name: (p.read_text(), p.stat().st_mtime_ns) for p in s.glob("*.exit")}


def test_a_passing_reading_passes_the_gate_and_a_failing_uv_stops_the_tests_step(tmp_path):
    repo, run = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    out = _run(repo, home, head, "g1")
    s = _scratch(home, "g1")
    assert out.returncode == 10, out.stderr
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert sorted(p.name for p in s.glob("*.exit")) == ["gate.exit", "tests.exit"]
    assert (s / "head").read_text().strip() == head and (s / "step0_commit").read_text().strip() == run
    assert (repo / "runs" / f"embovsb_{head[:7]}_g1").is_dir()
    calls = (home / "uvcalls").read_text().splitlines()
    assert len(calls) == 1 and "pytest" in calls[0] and "tests/test_embryo_ovsb_launcher.py" in calls[0]
    before = _exit_files(s)
    again = _run(repo, home, head, "g1")  # a refused second call leaves the finished steps' exit files as they were
    assert again.returncode == 91 and _exit_files(s) == before
    assert len((home / "uvcalls").read_text().splitlines()) == 1


# ------------------------------------------------------------------------------------------------ every step, traced
FAKE_UV = """#!/bin/sh
# a stand-in uv: pytest passes, vec check passes, a python -c runs in the test's own python, and a build copies the
# fixture of its output's name and stamps it with HEAD, clean
echo "$*" >> "$HOME/uvcalls"
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) exec "$REAL_PY" -I -B "$@" ;;
      experiments/builds/build_test.py)
        out=""
        prev=""
        for a in "$@"; do
          if [ "$prev" = "--out" ]; then out="$a"; fi
          prev="$a"
        done
        cp "$FIXTURES/$(basename "$out")" "$out" || exit 1
        printf '{"commit": "%s", "clean": true}' "$(git rev-parse HEAD)" > "$out.stamp.json"
        exit 0 ;;
    esac ;;
esac
exit 1
"""


def test_every_step_runs_in_order_with_a_stand_in_uv(tmp_path):
    """The launcher as committed but for its three sha256 constants, put to the fixtures' (the builds are stand-ins
    that copy fixtures): every step writes 0, the stamps say HEAD and clean, the identity booleans are all True, one
    sha256 list covers the five files and the twin builds compare equal."""
    import hashlib

    rng = np.random.default_rng(5)
    fx = tmp_path / "fixtures"
    fx.mkdir()
    C = rng.normal(size=(50, 3)).astype(np.float32)
    names = {"O": "T2embryo_s2x_ov.h5ad", "SO": "T2embryo_s2xs_ov.h5ad", "OH": "T2embryo_s2x_ovhi.h5ad",
             "W": "T2embryo_s2x_ovsb.h5ad"}
    for name in names.values():
        _h5ad(fx / name, rng.gamma(2.0, 1.0, size=(50, 4)), C)
    sha = {k: hashlib.sha256((fx / n).read_bytes()).hexdigest() for k, n in names.items()}
    repo, _ = _checkout(tmp_path)
    sh = (repo / REL).read_text()
    for k in ("O", "SO", "OH"):
        sh = sh.replace(f"{k}_SHA={SHAS[k]}", f"{k}_SHA={sha[k]}")
        sums, entry = RECORDS[k]
        p = repo / sums
        p.write_text(p.read_text().replace(f"{SHAS[k]}  {entry}", f"{sha[k]}  {entry}"))
    (repo / REL).write_text(sh)
    git(repo, "commit", "-qam", "the fixtures' sha256")
    home = _home(tmp_path)
    (home / ".local" / "bin" / "uv").write_text(FAKE_UV)
    head = git(repo, "rev-parse", "HEAD")
    env = {**git_env(), "HOME": str(home), "REAL_PY": sys.executable, "FIXTURES": str(fx)}
    out = subprocess.run(["bash", REL, head, "f1"], cwd=repo, env=env, capture_output=True, text=True, timeout=300)
    s = _scratch(home, "f1")
    assert out.returncode == 0, (out.stderr, {p.name: p.read_text() for p in s.glob("*.exit")})
    steps = ["gate", "tests", "O", "SO", "OH", "Wa", "Wb", "identity", "done"]
    assert {p.name: p.read_text().strip() for p in s.glob("*.exit")} == {f"{k}.exit": "0" for k in steps}
    assert (s / "stamps").read_text().splitlines() == [f"{d}/{n.removesuffix('.h5ad')} {head} True" for d, n in
                                                       (("O", names["O"]), ("SO", names["SO"]), ("OH", names["OH"]),
                                                        ("Wa", names["W"]), ("Wb", names["W"]))]
    assert (s / "identity.out").read_text().split() == ["w_coords_equal_o", "True", "w_obs_names_equal_o", "True",
                                                         "w_var_names_equal_o", "True", "w_x_not_o", "True"]
    listed = [ln.split() for ln in (s / "SHA256SUMS").read_text().splitlines()]
    assert [p for _, p in listed] == ["O/T2embryo_s2x_ov.h5ad", "SO/T2embryo_s2xs_ov.h5ad", "OH/T2embryo_s2x_ovhi.h5ad",
                                      "Wa/T2embryo_s2x_ovsb.h5ad", "Wb/T2embryo_s2x_ovsb.h5ad"]
    assert (repo / "runs" / f"embovsb_{head[:7]}_f1" / "SHA256SUMS").read_text() == (s / "SHA256SUMS").read_text()
    calls = (home / "uvcalls").read_text().splitlines()
    builds = [c for c in calls if "build_test.py embryo" in c]
    assert [c[c.index("embryo") + 7:c.index(" --out")] for c in builds] == [OV, SO, OHI, W, W]


# ------------------------------------------------------------------------------------------------ the identity step
def _identity_code() -> str:
    m = re.search(r"STEP=identity\n.*?\$PY -c '(.+?)' \"\$O/O/", LAUNCHER.read_text(), flags=re.S)
    assert m, "the identity step is not where the test expects it"
    return m.group(1)


def _h5ad(path: Path, X, coords, names=None, genes=None) -> str:
    a = anndata.AnnData(X=np.asarray(X, np.float32),
                        obs=pd.DataFrame(index=names or [f"cell_{i}" for i in range(len(X))]),
                        var=pd.DataFrame(index=genes or [f"g{j}" for j in range(X.shape[1])]))
    a.obsm["spatial_3D"] = np.asarray(coords, np.float32)
    a.write_h5ad(path)
    return str(path)


def _identity(*paths) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-I", "-B", "-c", _identity_code(), *paths], capture_output=True, text=True,
                          timeout=120, env={**os.environ})


def test_the_identity_step_prints_booleans_only_and_holds_on_os_rows_with_other_expression(tmp_path):
    rng = np.random.default_rng(3)
    X_o, X_w = (rng.gamma(2.0, 1.0, size=(200, 5)).astype(np.float32) for _ in range(2))
    C = rng.normal(size=(200, 3)).astype(np.float32)
    o = _h5ad(tmp_path / "o.h5ad", X_o, C)
    out = _identity(o, _h5ad(tmp_path / "w.h5ad", X_w, C))
    assert out.returncode == 0, out.stderr
    lines = out.stdout.strip().splitlines()
    assert [ln.split()[0] for ln in lines] == ["w_coords_equal_o", "w_obs_names_equal_o", "w_var_names_equal_o",
                                                "w_x_not_o"]
    assert all(re.fullmatch(r"[a-z0-9_]+ True", ln) for ln in lines)
    nudged = C.copy()
    nudged[7, 2] = np.nextafter(nudged[7, 2], np.float32(np.inf))
    bad = _identity(o, _h5ad(tmp_path / "w1.h5ad", X_w, nudged))
    assert bad.returncode == 1 and "w_coords_equal_o False" in bad.stdout
    names = [f"cell_{i}" for i in range(200)]
    names[3], names[4] = names[4], names[3]
    bad = _identity(o, _h5ad(tmp_path / "w2.h5ad", X_w, C, names=names))
    assert bad.returncode == 1 and "w_obs_names_equal_o False" in bad.stdout
    bad = _identity(o, _h5ad(tmp_path / "w3.h5ad", X_o, C))  # O itself is not W
    assert bad.returncode == 1 and "w_x_not_o False" in bad.stdout
    assert all(re.fullmatch(r"[a-z0-9_]+ (True|False)", ln) for ln in bad.stdout.strip().splitlines())
