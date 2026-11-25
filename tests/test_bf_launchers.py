"""The DGX launchers of row R71's files and their gate row E59: experiments/builds/embryo_bf.sh (the test-board builds)
and experiments/coords/run_bf.sh (the gate). Their strings and step order are held here, and every refusal path is run
for real in throwaway git repositories under a temporary home, with a stand-in uv first on PATH (so nothing is built,
no competition file is read and no score is made); refusals are traced with bash -x, and each launcher is walked once
through its later steps with the stand-in. Neither launcher is run here with a real commit. Built through
experiments/builds/build_test.py."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from conftest import git, git_env

REPO = Path(__file__).resolve().parents[1]
BUILD_REL = "experiments/builds/embryo_bf.sh"
GATE_REL = "experiments/coords/run_bf.sh"
ROW_REL = "experiments/coords/bf_row.py"
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import bf_row as br  # noqa: E402

SHAS = {"O": "1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea",
        "SO": "77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165"}
RECORDS = {"O": ("experiments/builds/results/embov_3b813b5.SHA256SUMS", "OVa/T2embryo_s2x_ov.h5ad"),
           "SO": ("experiments/builds/results/embsohi_93f832e_r1.SHA256SUMS", "SOa/T2embryo_s2xs_ov.h5ad")}
READING_DIR = "experiments/coords/results/bf_gate"
R71 = ("| R71 (the captain's order of 2026-11-25, ruling, no run: \"My named exception on the embryo's overlay form "
       "and the z-extent frozen rule for each BF file, as R59, with the challenger's fixes.\"; the embryo's "
       "consensus-form files BF on O and on SO only; built only after their gate row reads PASS, only with "
       "`--captain-ruling R71`; not a try) | T2 embryo | synthetic |")
STUB = r'''#!/usr/bin/env python3
import json, os, subprocess, sys
args = sys.argv[1:]
with open(os.path.join(os.environ["HOME"], "uvcalls"), "a") as f:
    f.write(" ".join(args) + "\n")
env = lambda k: int(os.environ.get(k, "0"))
if os.environ.get("STUB_FAIL_ALL") == "1":
    sys.exit(1)
if "pytest" in args:
    sys.exit(env("STUB_PYTEST"))
if any(a.endswith("build_test.py") for a in args):
    out = args[args.index("--out") + 1]
    with open(out, "wb") as f:
        f.write(b"stand-in " + out.encode())
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    with open(out + ".stamp.json", "w") as f:
        json.dump({"commit": head, "clean": True}, f)
    sys.exit(0)
if "vec" in args and "check" in args:
    sys.exit(0)
if any(a.endswith("bf_read.py") for a in args):
    i = [k for k, a in enumerate(args) if a.endswith("bf_read.py")][0]
    sys.exit(env("STUB_" + args[i + 1].upper()))
if "-c" in args:
    i = args.index("-c")
    os.execvp("python3", ["python3", "-I", "-B", "-c", *args[i + 1:]])
sys.exit(1)
'''


def _e59(outcome="not yet run (placed before its run)", cells=11, first=None) -> str:
    first = br.FIRST_CELL_OPENING + " (synthetic))" if first is None else first
    return "| " + " | ".join([first, *[f"c{k}" for k in range(2, cells)], outcome]) + " |"


def _cell(key: str) -> str:
    return f"Run at the run commit; in the read's fixed words: \"{br.PASS_TEXTS[key]}\""


def _home(tmp_path: Path, name: str = "home") -> Path:
    home = tmp_path / name
    fake = home / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text(STUB)
    (fake / "uv").chmod(0o755)
    return home


def _run(rel: str, repo: Path, home: Path, *argv: str, trace=False, **stub) -> subprocess.CompletedProcess:
    env = {**git_env(), "HOME": str(home), **{k: str(v) for k, v in stub.items()}}
    return subprocess.run(["bash", *(["-x"] if trace else []), rel, *argv], cwd=repo, env=env, capture_output=True,
                          text=True, timeout=120)


def _calls(home: Path) -> list[str]:
    p = home / "uvcalls"
    return p.read_text().splitlines() if p.exists() else []


# ------------------------------------------------------------------------------------------------ the texts
def test_the_build_launchers_commands_pins_and_steps():
    sh = (REPO / BUILD_REL).read_text()
    args = dict(re.findall(r'^([A-Z]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    floor = "--method fgw_shape0 --positions floor_exact"
    assert args["OV"] == f"{floor} --expression mixture --shape-interp overlay --shape-interp-test-ok"
    assert args["SO"] == f"{floor} --expression series --shape-interp overlay --shape-interp-test-ok"
    assert f'OV_ARGS="{args["OV"]}"' in (REPO / "experiments/builds/embryo_overlay.sh").read_text()
    assert f'SO_ARGS="{args["SO"]}"' in (REPO / "experiments/builds/embryo_so_ohi.sh").read_text()
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    assert "embryo " + args["BFO"] == build_test.EMBRYO_BF_COMMANDS[br.FILE_O]
    assert "embryo " + args["BFSO"] == build_test.EMBRYO_BF_COMMANDS[br.FILE_SO]
    for k, sha in SHAS.items():
        sums, entry = RECORDS[k]
        assert f"{k}_SHA={sha}" in sh and f"{k}_SUMS={sums}" in sh and f"{k}_ENTRY={entry}" in sh
        rows = dict(reversed(line.split()) for line in (REPO / sums).read_text().splitlines() if line.strip())
        assert rows[entry] == sha
    assert not re.search(r"results/[a-z0-9_]+\.txt", sh) and f"READING_DIR={READING_DIR}" in sh
    order = []
    for m in re.finditer(r"^STEP=(\w+)$|^  for STEP in ([\w ]+); do$", sh, flags=re.M):
        order += [m.group(1)] if m.group(1) else m.group(2).split()
    assert order == ["gate", "tests", "O", "SO", "BFOa", "BFOb", "BFSOa", "BFSOb", "done"]
    assert '[ "$st" = "$head True" ] || return 15' in sh and "vec check" in sh
    assert 'cmp -s "$O/BFOa/$BFO_FILE.h5ad" "$O/BFOb/$BFO_FILE.h5ad" || finish 31' in sh


def test_both_launchers_send_output_to_dev_null_delete_nothing_and_parse():
    for rel in (BUILD_REL, GATE_REL):
        text = (REPO / rel).read_text()
        joined = text.replace("\\\n", " ")
        for line in joined.splitlines():
            if re.search(r"build_test\.py embryo|vec check \"|\$SCRIPT\" (g0|run|read|g0read) ", line):
                assert "/dev/null" in line, line
        code = [line.split("#")[0] for line in text.splitlines()]
        assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
        for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", text):
            assert "'" not in msg and '"' not in msg, msg
        assert "upload" not in " ".join(code)
        assert subprocess.run(["bash", "-n", str(REPO / rel)], capture_output=True).returncode == 0
        listed = re.findall(r"tests/test_[a-z0-9_]+\.py", text)
        assert "tests/test_bf_launchers.py" in listed and all((REPO / f).exists() for f in listed)


def test_the_r71_checks_agree_with_the_builders_own(tmp_path):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    codes = [re.search(r"python3 -I -B -c '([^']+)'", (REPO / rel).read_text()).group(1)
             for rel in (BUILD_REL, GATE_REL)]
    assert codes[0] == codes[1]
    cases = [(R71, 0), (R71 + "\n" + R71, 1), (_e59(), 1), (R71.replace("BF on O and on SO only", "BF"), 1),
             (R71.replace("`--captain-ruling R71`", "--captain-ruling R71"), 1),
             (R71.replace("ruling, no run", "x"), 1),
             (R71.replace("| R71 ", "| R71a "), 1)]
    for k, (rows, want) in enumerate(cases):
        p = tmp_path / f"TRIES{k}.md"
        p.write_text("| ID | what |\n|---|---|\n" + rows + "\n")
        got = subprocess.run(["python3", "-I", "-B", "-c", codes[0], str(p)], capture_output=True)
        assert got.returncode == want, rows
        assert build_test._captain_order_names("R71", build_test.EMBRYO_BF_WORDS, p) == (want == 0)
    real = subprocess.run(["python3", "-I", "-B", "-c", codes[0], str(REPO / "experiments/harness/TRIES.md")])
    assert real.returncode == 0  # R71's real first cell


# ------------------------------------------------------------------------------------------------ the build launcher
def _reading(commit: str, key="PASS_BOTH", **over) -> dict:
    r = {"void": False, "outcome": key, "outcome_text": br.PASS_TEXTS.get(key, "x"), "commit": commit,
         "files_passing": list(br.PASS_FILES.get(key, ())), "reading_code": {"head": commit, "clean": True}}
    txt = over.pop("_txt", f"Row E59 (run commit {commit}): synthetic\n{br.PASS_TEXTS.get(key, 'x')}\n")
    r.update(over)
    return {"txt": txt, "json": r}


def _checkout(tmp_path: Path, at_run=None, at_head=None, reading=lambda c: _reading(c), track=("txt", "json"),
              records=True) -> tuple[Path, str]:
    """A throwaway repository: the launcher, the row check, the rules file and its hash, the committed build records,
    a ledger with `at_run`'s rows committed (the gate's run commit), then `at_head`'s rows and the gate's reading."""
    at_run = (R71, _e59()) if at_run is None else at_run
    at_head = (R71, _e59(_cell("PASS_BOTH"))) if at_head is None else at_head
    repo = tmp_path / "checkout"
    sums = sorted({r[0] for r in RECORDS.values()})
    for rel in (BUILD_REL, ROW_REL, "docs/S2_RULES.md", "docs/S2_RULES.sha256", *sums):
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
    git(repo, "commit", "-qm", "gate read", "--allow-empty")
    return repo, run


def _scratch(home: Path, name: str, run: str) -> Path:
    return home / "v" / "scratch" / name / run


def _refused(repo, home, run, code, words):
    out = _run(BUILD_REL, repo, home, git(repo, "rev-parse", "HEAD"), run, trace=True, STUB_FAIL_ALL=1)
    assert out.returncode == code and words in out.stderr, (out.returncode, out.stderr[-2000:])
    assert f"finish {code}" in out.stderr
    assert (_scratch(home, "embbf", run) / "gate.exit").read_text().strip() == str(code)
    assert not (repo / "runs").exists() and not _calls(home)


def test_bad_arguments_a_wrong_commit_and_a_reused_folder_are_refused(tmp_path):
    repo, _ = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    assert _run(BUILD_REL, repo, home, head, "a/b").returncode == 90
    assert _run(BUILD_REL, repo, home, head[:12], "r0").returncode == 90
    assert _run(BUILD_REL, repo, home, head).returncode != 0 and not (home / "v").exists()
    out = _run(BUILD_REL, repo, home, "0" * 40, "r0", trace=True)
    assert out.returncode == 96 and "finish 96" in out.stderr
    assert not re.search(r"^\++ (uv|python3?) ", out.stderr, flags=re.M)
    again = _run(BUILD_REL, repo, home, head, "r0")
    assert again.returncode == 91 and "never reused" in again.stderr


def test_an_unclean_tree_or_a_changed_rules_file_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a")
    (repo / "experiments/harness/TRIES.md").write_text(R71 + "\n")
    _refused(repo, _home(tmp_path, "h1"), "u1", 98, "tracked files differ")
    repo, _ = _checkout(tmp_path / "b")
    with open(repo / "docs/S2_RULES.md", "a") as f:
        f.write("\nchanged\n")
    git(repo, "commit", "-qam", "rules changed")
    _refused(repo, _home(tmp_path, "h2"), "u2", 97, "")


def test_the_rows_and_the_outcome_cell_gate_the_build(tmp_path):
    cases = [((_e59(_cell("PASS_BOTH")),), 95, "row R71"),
             ((_e59(_cell("PASS_BOTH")), R71, R71), 95, "row R71"),
             ((R71,), 94, "well-formed row E59"),
             ((R71, _e59(_cell("PASS_BOTH"), cells=10)), 94, "well-formed row E59"),
             ((R71, _e59()), 92, "does not say that the gate passed"),
             ((R71, _e59(br.STOP_TEXT)), 92, "does not say"),
             ((R71, _e59(br.IDLE_TEXT)), 92, "does not say"),
             ((R71, _e59(br.G0_STOP_TEXT)), 92, "does not say"),
             ((R71, _e59(br.VOID_TEXT)), 92, "does not say"),
             ((R71, _e59(_cell("PASS_O") + " " + br.STOP_TEXT)), 92, "does not say")]
    for k, (rows, code, words) in enumerate(cases):
        repo, _ = _checkout(tmp_path / f"c{k}", at_head=rows)
        _refused(repo, _home(tmp_path, f"h{k}"), f"r{k}", code, words)


def test_the_committed_reading_gates_the_build(tmp_path):
    repo, _ = _checkout(tmp_path / "a", reading=None)
    _refused(repo, _home(tmp_path, "h1"), "n1", 89, "no committed reading")
    repo, _ = _checkout(tmp_path / "b", track=("json",))
    _refused(repo, _home(tmp_path, "h2"), "n2", 89, "reading.txt must be tracked")
    bad = [dict(key="STOP"), dict(void=True), dict(reading_code={"head": "x", "clean": True}),
           dict(_txt=f"{br.STOP_TEXT}\n"), dict(files_passing=[br.FILE_O])]
    for k, over in enumerate(bad):
        repo, _ = _checkout(tmp_path / f"p{k}", reading=lambda c, over=over: _reading(c, **dict(over)))
        _refused(repo, _home(tmp_path, f"hp{k}"), f"p{k}", 88, "does not say PASS")
    repo, _ = _checkout(tmp_path / "d", reading=lambda c: _reading(c, "PASS_O"))  # the cell says both files
    _refused(repo, _home(tmp_path, "hd"), "d1", 88, "pass different files")
    repo, _ = _checkout(tmp_path / "e", reading=lambda c: _reading("f" * 40))
    _refused(repo, _home(tmp_path, "he"), "e1", 87, "not an ancestor")
    repo, _ = _checkout(tmp_path / "f", at_run=(R71,))
    _refused(repo, _home(tmp_path, "hf"), "f1", 86, "row E59 is not in TRIES.md")
    repo, _ = _checkout(tmp_path / "g", records=False)
    _refused(repo, _home(tmp_path, "hg"), "g1", 85, "does not give")


def test_a_walk_through_the_build_steps_stops_at_the_identity_rebuild(tmp_path):
    repo, _ = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    out = _run(BUILD_REL, repo, home, head, "w1")
    s = _scratch(home, "embbf", "w1")
    assert out.returncode == 21 and "does not give O's committed file" in out.stderr
    assert [(s / f"{k}.exit").read_text().strip() for k in ("gate", "tests", "O")] == ["0", "0", "21"]
    assert (s / "files").read_text().split() == [br.FILE_O, br.FILE_SO]
    assert (s / "stamps").read_text().split() == ["O/T2embryo_s2x_ov", head, "True"]
    calls = _calls(home)
    assert any("pytest" in c and "tests/test_bf_build.py" in c for c in calls)
    assert any("build_test.py embryo" in c and "--shape-interp overlay " in c for c in calls)
    assert not any("overlay_bayes" in c for c in calls)  # the identity file first, and it stopped everything


# ------------------------------------------------------------------------------------------------ the gate launcher
def _gate_checkout(tmp_path: Path, rows=None) -> Path:
    rows = (R71, _e59()) if rows is None else rows
    repo = tmp_path / "gate"
    for rel in (GATE_REL, ROW_REL):
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, repo / rel)
    tries = repo / "experiments/harness/TRIES.md"
    tries.parent.mkdir(parents=True)
    tries.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in rows))
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "row E59")
    return repo


def _exits(home: Path, run: str) -> dict[str, str]:
    s = _scratch(home, "bf_gate", run)
    return {p.stem: p.read_text().strip() for p in s.glob("*.exit")}


def test_the_gate_launchers_preflight_refusals(tmp_path):
    repo = _gate_checkout(tmp_path / "a")
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    assert _run(GATE_REL, repo, home).returncode == 2
    assert _run(GATE_REL, repo, home, head, ".x").returncode == 90
    out = _run(GATE_REL, repo, home, "0" * 40, "p1", trace=True)
    assert out.returncode == 96 and "finish 96" in out.stderr and _exits(home, "p1") == {"preflight": "96",
                                                                                          "done": "96"}
    assert _run(GATE_REL, repo, home, head, "p1").returncode == 91
    (repo / "experiments/harness/TRIES.md").write_text("changed\n")
    assert _run(GATE_REL, repo, home, head, "p2").returncode == 98
    for k, (rows, code) in enumerate((((R71,), 95), ((R71, _e59(), _e59()), 95), ((_e59(),), 94),
                                      ((R71.replace("BF on O and on SO only", "BF"), _e59()), 94))):
        repo = _gate_checkout(tmp_path / f"r{k}", rows)
        h = _home(tmp_path, f"h{k}")
        assert _run(GATE_REL, repo, h, git(repo, "rev-parse", "HEAD"), "q").returncode == code
        assert _exits(h, "q")["preflight"] == str(code) and not _calls(h)


def test_the_gate_launchers_steps_with_a_stand_in(tmp_path):
    repo = _gate_checkout(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    h = _home(tmp_path, "tests")  # the tests fail: nothing after them runs
    assert _run(GATE_REL, repo, h, head, "t1", STUB_PYTEST=1).returncode == 1
    assert _exits(h, "t1") == {"preflight": "0", "tests": "1", "done": "1"}
    h = _home(tmp_path, "g0stop")  # G0 says STOP: the run never starts, the read writes the G0 STOP reading
    assert _run(GATE_REL, repo, h, head, "t2", STUB_G0READ=4).returncode == 4
    ex = _exits(h, "t2")
    assert ex["g0"] == "0" and ex["g0read"] == "4" and ex["read"] == "0" and ex["done"] == "4" and "run" not in ex
    calls = _calls(h)
    assert sum(" g0 --seed " in c for c in calls) == 10 and not any(" run --task " in c for c in calls)
    assert sum(" read --commit " in c for c in calls) == 1
    h = _home(tmp_path, "runfail")  # a task fails: the read is still called once (a VOID reading)
    assert _run(GATE_REL, repo, h, head, "t3", STUB_RUN=1).returncode == 1
    ex = _exits(h, "t3")
    assert ex["run"] == "1" and ex["read"] == "0" and ex["done"] == "1"
    h = _home(tmp_path, "walk")  # every step exits 0
    assert _run(GATE_REL, repo, h, head, "t4").returncode == 0
    ex = _exits(h, "t4")
    assert ex == {k: "0" for k in ("preflight", "tests", "g0", "g0read", "run", "read", "done")}
    calls = _calls(h)
    assert sum(" run --task " in c for c in calls) == 9 and sum(" g0 --seed " in c for c in calls) == 10
    assert all(f"--commit {head}" in c for c in calls if "experiments/coords/bf_read.py" in c)
