"""The test-board build launcher of row R59's file T2embryo_s2x_ovt23.h5ad: experiments/builds/embryo_ot23.sh.

Its gate builds the file only when row R59 is the captain's order naming `--mix-lean later` and `--captain-ruling R59`
(the builder's own check), row E58 (gate G1) is well formed and its outcome cell, read by code, holds the read's PASS
words, and G1's reading is committed, says PASS on a clean tree at its run commit, an ancestor of HEAD whose ledger
holds row E58; the identity pins of O and W are read from their committed SHA256SUMS lists (never from a .txt record),
and both files are rebuilt to those pins before the file is built twice and compared. These tests hold the launcher's
strings and step order, run every gate refusal for real in throwaway git repositories under a temporary home with a
stand-in uv that fails (so nothing is built and no competition file is read), trace the refusals with bash -x, and walk
every step once with a stand-in uv whose builds copy fixtures. The launcher is never run here with a real commit. Built
through experiments/builds/build_test.py."""

from __future__ import annotations

import hashlib
import importlib.util
import json
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
REL = "experiments/builds/embryo_ot23.sh"
LAUNCHER = REPO / REL
ROW_REL = "experiments/coords/e58_row.py"
FLOOR = "--method fgw_shape0 --positions floor_exact"
OV = f"{FLOOR} --expression mixture --shape-interp overlay --shape-interp-test-ok"
W = f"{FLOOR} --expression series_back --shape-interp overlay --shape-interp-test-ok --captain-ruling R48"
OT = f"{OV} --mix-lean later --captain-ruling R59"
SHAS = {"O": "1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea",
        "W": "f576b64b8b93d79b1269b580df31bd5c6e58246cd4f3eda12cba7d77387087ab"}
RECORDS = {"O": ("experiments/builds/results/embov_3b813b5.SHA256SUMS", "OVa/T2embryo_s2x_ov.h5ad"),
           "W": ("experiments/builds/results/embovsb_20b2381_w1.SHA256SUMS", "Wa/T2embryo_s2x_ovsb.h5ad")}
READING_DIR = "experiments/coords/results/e58_mix_lean"
R59 = ("| R59 (the captain's order of 2026-11-24, ruling, no run: \"Yes to T2embryo_s2x_ovt23.h5ad (Ot23): my named "
       "exception on the embryo's overlay form for this file, bound to its gate G1.\"; the embryo builder is to take "
       "`--mix-lean later` on O's command only with `--captain-ruling R59`; not a try) | T2 embryo | synthetic |")


def _er():
    spec = importlib.util.spec_from_file_location("e58_row_for_launcher", REPO / ROW_REL)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


ER = _er()


def _e58(outcome="not yet run (placed before its run)", cells=11, first=None) -> str:
    first = ER.FIRST_CELL_OPENING + " (synthetic))" if first is None else first
    return "| " + " | ".join([first, *[f"c{k}" for k in range(2, cells)], outcome]) + " |"


PASS_CELL = f"Run at the run commit; in the read's fixed words: \"{ER.PASS_TEXT}\""


# ------------------------------------------------------------------------------------------------ the launcher's text
def test_the_commands_are_the_committed_ones_and_the_files():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^([A-Z]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"OV": OV, "W": W, "OT": OT}
    assert f'OV_ARGS="{OV}"' in (REPO / "experiments/builds/embryo_overlay.sh").read_text()
    assert f'W_ARGS="{W}"' in (REPO / "experiments/builds/embryo_ovsb.sh").read_text()
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    # the builder's named command is O's without the default --expression mixture (tests/test_mix_lean_build.py
    # shows both spellings give the same bytes)
    assert build_test.EMBRYO_OT23_COMMAND == "embryo " + OT.replace(" --expression mixture", "")
    for k, sha in SHAS.items():
        assert f"{k}_SHA={sha}" in sh and f'[ "$got" = "${k}_SHA" ]' in sh
        sums, entry = RECORDS[k]
        assert f"{k}_SUMS={sums}" in sh and f"{k}_ENTRY={entry}" in sh
        assert sums.endswith(".SHA256SUMS")  # the pins come from SHA256SUMS lists, never from a .txt record
    assert not re.search(r"results/[a-z0-9_]+\.txt", sh)
    assert "BOARD=T2:embryo:test_interp" in sh and f"READING_DIR={READING_DIR}" in sh
    assert "build O T2embryo_s2x_ov $OV_ARGS" in sh and "build W T2embryo_s2x_ovsb $W_ARGS" in sh
    assert 'build "$STEP" T2embryo_s2x_ovt23 $OT_ARGS' in sh
    assert 'cmp -s "$O/OTa/T2embryo_s2x_ovt23.h5ad" "$O/OTb/T2embryo_s2x_ovt23.h5ad" || finish 31' in sh
    builds = [ln for ln in sh.splitlines() if "build_test.py embryo" in ln]
    assert builds and not any("--board" in ln for ln in builds)  # no dry run: the option refuses one


def test_the_committed_records_hold_the_two_sha256():
    for k, (sums, entry) in RECORDS.items():
        rows = dict(reversed(line.split()) for line in (REPO / sums).read_text().splitlines() if line.strip())
        assert rows[entry] == SHAS[k]


def test_the_steps_run_in_order_each_writing_running_first():
    sh = LAUNCHER.read_text()
    order = []
    for m in re.finditer(r"^STEP=(\w+)$|^for STEP in ([\w ]+); do$", sh, flags=re.M):
        order += [m.group(1)] if m.group(1) else m.group(2).split()
    assert order == ["gate", "tests", "O", "W", "OTa", "OTb", "done"]
    assert sh.count("echo running > \"$S/$STEP.exit\"") == 6  # gate, tests, O, W, the loop, done
    assert sh.index('echo running > "$S/$STEP.exit"') < sh.index("head=$(git rev-parse HEAD)")
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh
    i8 = sh.index("OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8")
    i16 = sh.index("OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16")
    assert i8 < sh.index("STEP=tests") < i16 < sh.index("STEP=O\n")
    assert "S=$HOME/v/scratch/embot23/$RUN" in sh and "O=runs/embot23_${head:0:7}_$RUN" in sh
    assert "sha256sum O/*.h5ad W/*.h5ad OTa/*.h5ad OTb/*.h5ad) > \"$S/SHA256SUMS\"" in sh
    # the identity files are rebuilt before the file
    assert sh.index("STEP=O\n") < sh.index("STEP=W\n") < sh.index("for STEP in OTa OTb; do")


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
    assert {"tests/test_embryo_ot23_launcher.py", "tests/test_mix_lean_build.py",
            "tests/test_e58_mix_lean_read.py"} <= set(listed) and all((REPO / f).exists() for f in listed)
    assert "tests/test_embryo_ot23_launcher.py" in (REPO / "experiments/external/n54_step.sh").read_text()


# ------------------------------------------------------------------------------------------- the R59 check
def _snippet(start: str) -> str:
    m = re.search(re.escape(start) + r"'([^']+)'", LAUNCHER.read_text())
    assert m, f"the check after {start!r} is not where the test expects it"
    return m.group(1)


def _python(code: str, *paths) -> subprocess.CompletedProcess:
    return subprocess.run(["python3", "-I", "-B", "-c", code, *map(str, paths)], capture_output=True, text=True)


def test_the_r59_check_agrees_with_the_builders_own(tmp_path):
    code = _snippet("python3 -I -B -c ")
    cases = [(R59, 0), (_e58() + "\n" + R59, 0), (_e58(), 1), (R59 + "\n" + R59, 1),
             (R59.replace("`--mix-lean later`", "--mix-lean later"), 1),
             (R59.replace("`--captain-ruling R59`", "--captain-ruling R59"), 1),
             (R59.replace("ruling, no run", "ruling"), 1), (R59.replace("the captain's order", "the order"), 1),
             ("| R59 (the captain's order, ruling, no run) | `--mix-lean later` `--captain-ruling R59` |", 1),
             (R59.replace("| R59 ", "| R59a "), 1)]
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    for k, (rows, want) in enumerate(cases):
        p = tmp_path / f"TRIES{k}.md"
        p.write_text("| ID | what |\n|---|---|\n" + rows + "\n")
        assert _python(code, p).returncode == want, rows
        got = build_test._captain_order_names("R59", build_test.EMBRYO_OT23_OPTION, p)
        assert got == (want == 0), k
    real = _python(code, REPO / "experiments/harness/TRIES.md")
    assert real.returncode == 0  # R59's real first cell, by this launcher's check


# ------------------------------------------------------------------------------------------------ the gate, for real
def _reading(commit: str, **over) -> dict:
    r = {"void": False, "outcome": "PASS", "outcome_text": ER.PASS_TEXT, "passes": True, "commit": commit,
         "reading_code": {"head": commit, "clean": True}}
    txt = over.pop("_txt", f"Row E58 (run commit {commit}): synthetic\n{ER.PASS_TEXT}\n")
    r.update(over)
    return {"txt": txt, "json": r}


def _checkout(tmp_path: Path, at_run=None, at_head=None, reading=lambda c: _reading(c), track=("txt", "json"),
              records=True) -> tuple[Path, str]:
    """A throwaway repository: the launcher, the row check, the rules file and its hash, the committed build records
    and a ledger with `at_run`'s rows (default R59 and a well-formed E58 before its run), committed (G1's run commit);
    then the ledger with `at_head`'s rows (default: R59 and E58 with the PASS words in its outcome cell) and G1's
    reading made at that commit (unless `reading` is None), the files in `track` committed."""
    at_run = (R59, _e58()) if at_run is None else at_run
    at_head = (R59, _e58(PASS_CELL)) if at_head is None else at_head
    repo = tmp_path / "checkout"
    for rel in (REL, ROW_REL, "docs/S2_RULES.md", "docs/S2_RULES.sha256", *sorted({r[0] for r in RECORDS.values()})):
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
    git(repo, "commit", "-qm", "G1 read", "--allow-empty")
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
    return home / "v" / "scratch" / "embot23" / run


def _refused(repo, home, run, code, words, trace=True):
    out = _run(repo, home, git(repo, "rev-parse", "HEAD"), run, trace=trace)
    assert out.returncode == code and words in out.stderr, (out.returncode, out.stderr[-3000:])
    assert f"finish {code}" in out.stderr  # the trace shows the refusal's own line
    s = _scratch(home, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert not (repo / "runs").exists() and not (home / "uvcalls").exists()
    return out


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
    (repo / "experiments/harness/TRIES.md").write_text(R59 + "\n" + _e58(PASS_CELL) + "\n| E59 changed |\n")
    _refused(repo, _home(tmp_path, "h1"), "u1", 98, "tracked files differ")
    repo, _ = _checkout(tmp_path / "b")
    with open(repo / "docs/S2_RULES.md", "a") as f:
        f.write("\nchanged\n")
    git(repo, "commit", "-qam", "rules changed")
    _refused(repo, _home(tmp_path, "h2"), "u2", 97, "")


def test_without_row_r59_as_the_builder_needs_it_the_gate_refuses(tmp_path):
    e58 = _e58(PASS_CELL)
    cases = [(e58,), (e58, R59, R59), (e58, R59.replace("`--mix-lean later`", "--mix-lean later")),
             (e58, R59.replace("ruling, no run", "ruling"))]
    for k, rows in enumerate(cases):
        repo, _ = _checkout(tmp_path / f"c{k}", at_head=rows)
        _refused(repo, _home(tmp_path, f"h{k}"), f"r{k}", 95, "row R59")


def test_without_a_well_formed_row_e58_the_gate_refuses(tmp_path):
    cases = [(R59,), (R59, _e58(PASS_CELL), _e58(PASS_CELL)), (R59, _e58(PASS_CELL, cells=10)),
             (R59, _e58(PASS_CELL, first="E58 (gate G1; it moves the mix)"))]
    for k, rows in enumerate(cases):
        repo, _ = _checkout(tmp_path / f"c{k}", at_head=rows)
        _refused(repo, _home(tmp_path, f"h{k}"), f"e{k}", 94, "well-formed row E58")


def test_an_outcome_cell_without_the_pass_words_is_refused(tmp_path):
    cells = ["not yet run", ER.STOP_TEXT, ER.IDLE_TEXT, ER.VOID_TEXT, PASS_CELL + " " + ER.STOP_TEXT,
             ER.PASS_TEXT.replace("PASS", "Pass")]
    for k, cell in enumerate(cells):
        repo, _ = _checkout(tmp_path / f"c{k}", at_head=(R59, _e58(cell)))
        _refused(repo, _home(tmp_path, f"h{k}"), f"o{k}", 92, "does not say that G1 passed")


def test_without_a_committed_reading_the_gate_refuses(tmp_path):
    repo, _ = _checkout(tmp_path / "a", reading=None)
    _refused(repo, _home(tmp_path, "h1"), "n1", 89, "no committed reading")
    repo, _ = _checkout(tmp_path / "b", track=("json",))  # the text is there but not tracked
    _refused(repo, _home(tmp_path, "h2"), "n2", 89, "reading.txt must be tracked")
    repo, _ = _checkout(tmp_path / "c", track=("txt",))
    _refused(repo, _home(tmp_path, "h3"), "n3", 89, "reading.json must be tracked")


def test_a_reading_that_is_not_a_pass_is_refused(tmp_path):
    bad = [dict(outcome="STOP"), dict(outcome_text=ER.STOP_TEXT), dict(void=True), dict(passes=False),
           dict(reading_code={"head": "x", "clean": True}), dict(_txt=f"{ER.STOP_TEXT}\n"),
           dict(_txt=f"{ER.PASS_TEXT}\n{ER.VOID_TEXT}\n")]
    for k, over in enumerate(bad):
        repo, _ = _checkout(tmp_path / f"c{k}", reading=lambda c, over=over: _reading(c, **dict(over)))
        _refused(repo, _home(tmp_path, f"h{k}"), f"p{k}", 88, "does not say PASS")
    repo, _ = _checkout(tmp_path / "short", reading=lambda c: _reading(c[:7]))
    _refused(repo, _home(tmp_path, "hs"), "ps", 88, "does not say PASS")


def test_a_foreign_run_commit_or_one_without_row_e58_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a", reading=lambda c: _reading("f" * 40))
    _refused(repo, _home(tmp_path, "h1"), "c1", 87, "not an ancestor")
    repo, _ = _checkout(tmp_path / "b", at_run=(R59,))  # row E58 came after G1 ran
    _refused(repo, _home(tmp_path, "h2"), "c2", 86, "row E58 is not in TRIES.md at G1")


def test_a_missing_or_different_build_record_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a", records=False)
    _refused(repo, _home(tmp_path, "h1"), "s1", 85, "does not give O_SHA")
    repo, _ = _checkout(tmp_path / "b")
    sums = repo / RECORDS["W"][0]
    sums.write_text(sums.read_text().replace(SHAS["W"], "0" * 64))
    git(repo, "commit", "-qam", "record changed")
    _refused(repo, _home(tmp_path, "h2"), "s2", 85, "does not give W_SHA")


def test_an_existing_output_folder_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    (repo / "runs" / f"embot23_{head[:7]}_o1").mkdir(parents=True)
    out = _run(repo, _home(tmp_path), head, "o1")
    assert out.returncode == 93 and "output folder" in out.stderr


def test_a_passing_gate_and_a_failing_uv_stops_the_tests_step(tmp_path):
    repo, run = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    out = _run(repo, home, head, "g1")
    s = _scratch(home, "g1")
    assert out.returncode == 10, out.stderr
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert sorted(p.name for p in s.glob("*.exit")) == ["gate.exit", "tests.exit"]
    assert (s / "head").read_text().strip() == head and (s / "g1_run_commit").read_text().strip() == run
    assert (repo / "runs" / f"embot23_{head[:7]}_g1").is_dir()
    calls = (home / "uvcalls").read_text().splitlines()
    assert len(calls) == 1 and "pytest" in calls[0] and "tests/test_embryo_ot23_launcher.py" in calls[0]


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


def _h5ad(path: Path, X, coords) -> str:
    a = anndata.AnnData(X=np.asarray(X, np.float32), obs=pd.DataFrame(index=[f"cell_{i}" for i in range(len(X))]),
                        var=pd.DataFrame(index=[f"g{j}" for j in range(X.shape[1])]))
    a.obsm["spatial_3D"] = np.asarray(coords, np.float32)
    a.write_h5ad(path)
    return str(path)


def test_every_step_runs_in_order_with_a_stand_in_uv(tmp_path):
    """The launcher as committed but for its two sha256 constants, put to the fixtures' (the builds are stand-ins that
    copy fixtures): every step writes 0, the stamps say HEAD and clean, one sha256 list covers the four files and the
    twin builds compare equal; the build calls are O's, W's and the file's command twice, in that order."""
    rng = np.random.default_rng(5)
    fx = tmp_path / "fixtures"
    fx.mkdir()
    names = {"O": "T2embryo_s2x_ov.h5ad", "W": "T2embryo_s2x_ovsb.h5ad", "OT": "T2embryo_s2x_ovt23.h5ad"}
    for name in names.values():
        _h5ad(fx / name, rng.gamma(2.0, 1.0, size=(50, 4)), rng.normal(size=(50, 3)))
    sha = {k: hashlib.sha256((fx / n).read_bytes()).hexdigest() for k, n in names.items()}
    repo, _ = _checkout(tmp_path)
    sh = (repo / REL).read_text()
    for k in ("O", "W"):
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
    steps = ["gate", "tests", "O", "W", "OTa", "OTb", "done"]
    assert {p.name: p.read_text().strip() for p in s.glob("*.exit")} == {f"{k}.exit": "0" for k in steps}
    assert (s / "stamps").read_text().splitlines() == [f"{d}/{n.removesuffix('.h5ad')} {head} True" for d, n in
                                                       (("O", names["O"]), ("W", names["W"]), ("OTa", names["OT"]),
                                                        ("OTb", names["OT"]))]
    listed = [ln.split() for ln in (s / "SHA256SUMS").read_text().splitlines()]
    assert [p for _, p in listed] == ["O/T2embryo_s2x_ov.h5ad", "W/T2embryo_s2x_ovsb.h5ad",
                                      "OTa/T2embryo_s2x_ovt23.h5ad", "OTb/T2embryo_s2x_ovt23.h5ad"]
    assert listed[2][0] == listed[3][0] == sha["OT"]  # one sha256 for the twin builds
    assert (repo / "runs" / f"embot23_{head[:7]}_f1" / "SHA256SUMS").read_text() == (s / "SHA256SUMS").read_text()
    calls = (home / "uvcalls").read_text().splitlines()
    builds = [c for c in calls if "build_test.py embryo" in c]
    assert [c[c.index("embryo") + 7:c.index(" --out")] for c in builds] == [OV, W, OT, OT]


def test_a_build_that_misses_its_pin_stops_before_the_file(tmp_path):
    """O's rebuild giving another sha256 than its pin stops the run at step O (21): W and the file are never built."""
    rng = np.random.default_rng(6)
    fx = tmp_path / "fixtures"
    fx.mkdir()
    for name in ("T2embryo_s2x_ov.h5ad", "T2embryo_s2x_ovsb.h5ad", "T2embryo_s2x_ovt23.h5ad"):
        _h5ad(fx / name, rng.gamma(2.0, 1.0, size=(20, 3)), rng.normal(size=(20, 3)))
    repo, _ = _checkout(tmp_path)
    home = _home(tmp_path)
    (home / ".local" / "bin" / "uv").write_text(FAKE_UV)
    head = git(repo, "rev-parse", "HEAD")
    env = {**git_env(), "HOME": str(home), "REAL_PY": sys.executable, "FIXTURES": str(fx)}
    out = subprocess.run(["bash", REL, head, "m1"], cwd=repo, env=env, capture_output=True, text=True, timeout=300)
    s = _scratch(home, "m1")
    assert out.returncode == 21 and (s / "O.exit").read_text().strip() == "21"
    assert not (s / "W.exit").exists() and not (s / "OTa.exit").exists()
    builds = [c for c in (home / "uvcalls").read_text().splitlines() if "build_test.py embryo" in c]
    assert len(builds) == 1 and OV in builds[0]
