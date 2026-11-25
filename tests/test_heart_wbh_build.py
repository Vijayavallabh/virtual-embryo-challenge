"""experiments/builds/heart_wbh_build.sh, the test-board builds of file WBH (row R69, the captain's order of
2026-11-25; the row fixed, no ROW argument), only when row H37's outcome cell, read by code, holds the read's PASS line.
Every call runs in a throwaway git repository under the test's own folder, with a temporary home and the stand-ins of
tests/test_run_h37_launcher.py first on PATH: no launcher runs with a real commit, nothing reads a competition file and
nothing is built (build_test.py, halve.py and world_blur.py are never run). Each refusal of the gate step is
executed, then the identity step's stops (bytes, stamps), a failing tests step, and a whole call through the
stand-ins, which lists every file and compares the two builds."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_run_h37_launcher import R69, _env, h37_line, make_repo  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_wbh_build.sh"
RECORDS = Path("experiments") / "harness" / "results" / "h37_heart_world_blur"
PASS_LINE = "OUTCOME H37: PASS (world pattern)."
GATE_FILES = ("src/a.py", "configs/c.yaml", "experiments/builds/build_test.py", "experiments/builds/halve.py",
              "experiments/builds/world_blur.py", "experiments/harness/h35_heart_trend_forms.py")


def records(repo, run_commit, label="world pattern", **over):
    g0 = {"row": "H37", "axis": "h37-heart-world-blur:g0", "commit": run_commit, "verdict": "PASS",
          "provenance": {"commit": run_commit, "tree_clean": True}, **over.get("g0", {})}
    read = {"row": "H37", "commit": run_commit, "void": False, "outcome": "PASS", "label": label,
            **over.get("read", {})}
    (repo / RECORDS / "g0").mkdir(parents=True, exist_ok=True)
    (repo / RECORDS / "g0" / "g0.json").write_text(json.dumps(g0))
    (repo / RECORDS / "reading.json").write_text(json.dumps(read))


@pytest.fixture()
def repo(tmp_path):
    """The gate's run commit (row H37 open), then a commit with its PASS line and records; returns HEAD's."""
    repo, git, run_commit = make_repo(tmp_path, files=GATE_FILES)
    (repo / "docs").mkdir()
    (repo / "docs" / "S2_RULES.md").write_text("rules\n")
    sha = subprocess.run(["sha256sum", "docs/S2_RULES.md"], cwd=repo, capture_output=True, text=True).stdout
    (repo / "docs" / "S2_RULES.sha256").write_text(sha)
    (repo / ".gitignore").write_text("runs/\n.venv/\ndata/\n")
    git("add", ".")
    git("commit", "-q", "-m", "rules")
    run_commit = git("rev-parse", "HEAD")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(R69 + h37_line(f'In the read\'s words: "{PASS_LINE}"'))
    records(repo, run_commit)
    git("add", ".")
    git("commit", "-q", "-m", "the outcome and the records")
    return repo, git, run_commit, git("rev-parse", "HEAD")


def launch(tmp_path, *args, extra=None):
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=tmp_path / "repo", capture_output=True, text=True,
                       env=_env(tmp_path, extra), timeout=180)
    return r.returncode, r.stderr


def exit_of(tmp_path, run, step):
    p = tmp_path / "home" / "v" / "scratch" / "heartwbh" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


def recommit(repo, git, ledger=None, **rec):
    if ledger is not None:
        (repo / "experiments" / "harness" / "TRIES.md").write_text(ledger)
    if rec:
        records(repo, rec.pop("run_commit"), **rec)
    git("commit", "-q", "-am", "change")
    return git("rev-parse", "HEAD")


def test_odd_arguments_and_a_reused_run_folder(tmp_path, repo):
    _, _, _, head = repo
    assert launch(tmp_path, head, "r1", "R69")[0] == 90
    assert launch(tmp_path, head, "r/1")[0] == 90
    assert launch(tmp_path, head[:7], "r1")[0] == 90
    assert launch(tmp_path, "b" * 40, "r1")[0] == 96
    assert launch(tmp_path, head, "r1")[0] == 91


def test_a_clean_tree_and_the_rules_pin(tmp_path, repo):
    r, git, _, head = repo
    (r / "docs" / "S2_RULES.md").write_text("other\n")
    assert launch(tmp_path, head, "r1")[0] == 98
    head = recommit(r, git)
    assert launch(tmp_path, head, "r2")[0] == 97


@pytest.mark.parametrize("ledger,code", [
    (h37_line(f'"{PASS_LINE}"'), 95),                                          # no row R69
    (R69, 92),                                                                 # no row H37
    (R69 + h37_line("pending"), 92),                                           # no outcome yet
    (R69 + h37_line(f'"{PASS_LINE}"', first="H37 (another first cell)"), 92),  # not the fixed first cell
    (R69 + h37_line("OUTCOME H37: VOID: v4: x; nothing is built."), 87),
    (R69 + h37_line("OUTCOME H37: STOP: IDLE."), 86),
    (R69 + h37_line("OUTCOME H37: STOP at G0: G00's cut (1 - e^(-c/2)) / 2 on E9.5 is under 0.10."), 86),
])
def test_the_rows(tmp_path, repo, ledger, code):
    r, git, _, _ = repo
    head = recommit(r, git, ledger)
    assert launch(tmp_path, head, "r1")[0] == code and exit_of(tmp_path, "r1", "gate") == str(code)


def test_the_gate_records(tmp_path, repo):
    r, git, run_commit, _ = repo
    head = recommit(r, git, run_commit=run_commit, label="spread only")  # the cell says world pattern
    assert launch(tmp_path, head, "r1")[0] == 84
    head = recommit(r, git, run_commit=run_commit, read={"void": True})
    assert launch(tmp_path, head, "r2")[0] == 84
    head = recommit(r, git, run_commit=run_commit, g0={"provenance": {"commit": run_commit, "tree_clean": False}})
    assert launch(tmp_path, head, "r3")[0] == 84
    head = recommit(r, git, run_commit="c" * 40)  # not in the history
    assert launch(tmp_path, head, "r4")[0] == 84
    head = recommit(r, git, run_commit=run_commit)
    (r / "experiments" / "builds" / "world_blur.py").write_text("x = 2\n")
    head = recommit(r, git)
    assert launch(tmp_path, head, "r5")[0] == 85


def test_the_output_folder_and_the_tests(tmp_path, repo):
    r, _, _, head = repo
    (r / "runs" / f"heartwbh_{head[:7]}_r1").mkdir(parents=True)
    assert launch(tmp_path, head, "r1")[0] == 93
    assert launch(tmp_path, head, "r2", extra={"STUB_PYTEST": "1"})[0] == 10
    assert exit_of(tmp_path, "r2", "gate") == "0" and exit_of(tmp_path, "r2", "OTK125") is None


def test_the_identity_step_and_the_stamps(tmp_path, repo):
    _, _, _, head = repo
    stamp = {"STUB_STAMP": f"{head} True"}
    assert launch(tmp_path, head, "r1", extra={**stamp, "STUB_SHA_FILE": "0" * 64})[0] == 20
    assert launch(tmp_path, head, "r2", extra={**stamp, "STUB_SHA_HALF": "0" * 64})[0] == 21
    assert launch(tmp_path, head, "r3", extra={"STUB_STAMP": f"{head} False"})[0] == 15
    assert exit_of(tmp_path, "r3", "WBHa") is None


def test_a_whole_call_builds_the_file_twice(tmp_path, repo):
    r, _, _, head = repo
    assert launch(tmp_path, head, "r1", extra={"STUB_STAMP": f"{head} True"})[0] == 0
    for step in ("gate", "tests", "OTK125", "WBHa", "WBHb", "done"):
        assert exit_of(tmp_path, "r1", step) == "0", step
    calls = (tmp_path / "home" / "uvcalls").read_text()
    wbh = ("experiments/builds/world_blur.py runs/heartwbh_{h}_r1/OTK125/T2heart_organ_ts125_tk_half.h5ad "
           "runs/heartwbh_{h}_r1/{d}/T2heart_organ_ts125_tk_wb_half.h5ad --prior same_stage --blocks 32 "
           "--seed 1000001 --captain-ruling R69")
    for d in ("WBHa", "WBHb"):
        assert wbh.format(h=head[:7], d=d) in calls
    assert calls.count("vec check") == 2
    sums = (r / "runs" / f"heartwbh_{head[:7]}_r1" / "SHA256SUMS").read_text()
    assert "WBHa/T2heart_organ_ts125_tk_wb_half.h5ad" in sums and "OTK125/T2heart_organ_ts125_tk.h5ad" in sums
