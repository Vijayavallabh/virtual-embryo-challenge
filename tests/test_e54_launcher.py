"""Row E54's launcher experiments/coords/run_e54.sh: its steps and order, that it removes nothing, that its fold and
task lists and its allowed log lines are the script's, and its refusal paths run for real under a temporary home (no
arguments, an unknown step, a step whose predecessor has not exited 0 at the commit, a wrong commit, an unclean tree
and a fit without a GPU). No step that runs anything is reached: every call here stops at a refusal."""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "coords" / "run_e54.sh"
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import e54_mmflow_read as e54  # noqa: E402


def _run(home: Path, *argv: str, cuda: str | None = None) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items()
           if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "CUDA_VISIBLE_DEVICES")}
    env["HOME"] = str(home)
    if cuda is not None:
        env["CUDA_VISIBLE_DEVICES"] = cuda
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO, env=env, capture_output=True, text=True)


def _head() -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()


def _clean() -> bool:
    return subprocess.run(["git", "status", "--porcelain"], cwd=REPO, capture_output=True, text=True).stdout == ""


def _exit(home: Path, step: str) -> str:
    return (home / "v" / "scratch" / "e54" / f"{step}.exit").read_text().strip()


def test_the_steps_their_order_and_their_lists_are_the_scripts():
    sh = LAUNCHER.read_text()
    order = dict(re.findall(r"^  (\w+)\) PREV=\"(\w*)\" ;;$", sh, flags=re.M))
    assert order == {"tests": "", "fit": "tests", "repro": "fit", "gate": "repro", "run": "gate", "read": "run"}
    tasks = re.search(r"^TASKS=\(([^)]*)\)$", sh, flags=re.M).group(1).split()
    folds = re.search(r"^FOLDS=\(([^)]*)\)$", sh, flags=re.M).group(1).split()
    assert sorted(tasks) == sorted(e54.nt.all_tasks()) and sorted(folds) == sorted(e54.folds())
    assert "SCRIPT=experiments/coords/e54_mmflow_read.py" in sh and "OUT=experiments/coords/results/e54_mmflow" in sh
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert "CUBLAS_WORKSPACE_CONFIG=:4096:8" in sh
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh.split("  tests)\n")[1].split(";;")[0])
    assert {"tests/test_t2_mmflow.py", "tests/test_e54_mmflow_read.py", "tests/test_e54_launcher.py",
            "tests/test_harness.py"} <= set(listed)
    assert all((REPO / f).exists() for f in listed)
    assert 'import torch, torchdiffeq' in sh and "finish 94" in sh  # the flow's tests cannot be skipped there


def test_the_allowed_log_lines_are_the_scripts_prints():
    ok = re.compile(re.search(r"^OK_LINE='([^']+)'$", LAUNCHER.read_text(), flags=re.M).group(1))
    for line in ("fit x80_825_late done", "repro e725 pass", "repro x1 fail", "gate pass", "gate fail",
                 "new2 seed 9 done", "wrote runs/e54/run/half"):
        assert ok.fullmatch(line), line
    for line in ("half seed 3: floor 50.00, fgw_shape0x 51.20", "n_cells 1234", "gate maybe", "Traceback"):
        assert not ok.fullmatch(line), line
    src = (REPO / "experiments" / "coords" / "e54_mmflow_read.py").read_text()
    prints = set(re.findall(r"print\((f\"[^\"]*\")", src))
    assert prints == {'f"fit {fold} done"', 'f"repro {a.task} {\'pass\' if rec[\'pass\'] else \'fail\'}"',
                      'f"gate {\'pass\' if g[\'pass\'] else \'fail\'}"', 'f"{a.task} seed {s} done"',
                      'f"wrote {folder}"', 'f"wrote {files[1]}"'}, prints
    assert src.count("print(") == 7  # those six, the fit's twice


def test_the_launcher_removes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in LAUNCHER.read_text().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", LAUNCHER.read_text()):
        assert "'" not in msg and '"' not in msg, msg


def test_no_arguments_and_an_unknown_step(tmp_path):
    out = _run(tmp_path)
    assert out.returncode == 2 and (tmp_path / "v" / "scratch" / "e54" / "usage.exit").read_text().strip() == "2"
    out = _run(tmp_path, "bogus", "0" * 40)
    assert out.returncode == 97 and _exit(tmp_path, "bogus") == "97"


def test_each_step_needs_its_predecessor_at_the_same_commit(tmp_path):
    for step in ("fit", "repro", "gate", "run", "read"):
        out = _run(tmp_path, step, _head())
        assert out.returncode == 96 and _exit(tmp_path, step) == "96" and "has not exited 0" in out.stderr
    s = tmp_path / "v" / "scratch" / "e54"
    (s / "repro.exit").write_text("0\n")
    (s / "repro.head").write_text("0" * 40 + "\n")  # exited 0, but at another commit
    assert _run(tmp_path, "gate", _head()).returncode == 96
    (s / "repro.head").write_text(_head() + "\n")
    (s / "repro.exit").write_text("3\n")  # at this commit, but a miss
    assert _run(tmp_path, "gate", _head()).returncode == 96


def test_a_wrong_commit_or_an_unclean_tree_is_refused(tmp_path):
    out = _run(tmp_path, "tests", "0" * 40)
    want = 99 if _clean() else 98
    assert out.returncode == want and _exit(tmp_path, "tests") == str(want)
    assert not (tmp_path / "v" / "scratch" / "e54" / "tests.head").exists()  # nothing ran


def test_the_fits_refuse_without_a_gpu(tmp_path):
    s = tmp_path / "v" / "scratch" / "e54"
    s.mkdir(parents=True)
    (s / "tests.exit").write_text("0\n")
    (s / "tests.head").write_text(_head() + "\n")
    for cuda in (None, ""):
        out = _run(tmp_path, "fit", _head(), cuda=cuda)
        want = 95 if _clean() else 98
        assert out.returncode == want and _exit(tmp_path, "fit") == str(want)


def _parallel(tmp_path: Path, step: str, fake_exit: dict[str, int], cuda: str = "4,6") -> int:
    """The launcher's own parallel() text, run with a fake interpreter that records its arguments and the GPU it was
    given and exits as told; the launcher's lists and allowed log lines, a temporary runs folder. No git, no python."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    sh = LAUNCHER.read_text()
    func = re.search(r"^parallel\(\) \{.*?^\}$", sh, flags=re.M | re.S).group(0)
    lists = "\n".join(re.findall(r"^(?:TASKS|FOLDS|OK_LINE)=.*$", sh, flags=re.M))
    fake = tmp_path / "fakepy"
    codes = " ".join(f"{k}) code={v} ;;" for k, v in fake_exit.items())
    fake.write_text('echo "$CUDA_VISIBLE_DEVICES $2 $4" >> "$(dirname "$0")/calls"\n'
                    f'case "$4" in {codes} *) code=0 ;; esac\n'
                    'if [ "$2" = fit ]; then echo "fit $4 done"; else echo "repro $4 pass"; echo "n_cells 12"; fi\n'
                    'exit $code\n')
    script = (f'set -u\nR="{tmp_path}/runs"\nPY="bash {fake}"\nSCRIPT=x.py\nCOMMIT=c\nDATA=d\nOUT=o\n'
              f'export CUDA_VISIBLE_DEVICES={cuda}\n{lists}\n{func}\nparallel {step}\nexit $?\n')
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True).returncode


def test_parallel_spreads_the_fits_over_the_listed_gpus(tmp_path):
    assert _parallel(tmp_path, "fit", {}) == 0
    calls = [c.split() for c in (tmp_path / "calls").read_text().splitlines()]
    assert sorted(c[2] for c in calls) == sorted(e54.folds()) and all(c[1] == "fit" for c in calls)
    assert {c[0] for c in calls} == {"4", "6"}  # one GPU id per process, from the list
    ex = tmp_path / "runs" / "exit"
    assert all((ex / f"fit-{f}").read_text().strip() == "0" for f in e54.folds())
    assert all((ex / f"fit-{f}.other_lines").read_text().strip() == "0" for f in e54.folds())


def test_parallel_returns_3_on_a_miss_and_1_on_an_error_and_counts_other_lines(tmp_path):
    assert _parallel(tmp_path / "a", "repro", {"x2": 3}) == 3
    ex = tmp_path / "a" / "runs" / "exit"
    assert (ex / "repro-x2").read_text().strip() == "3"
    assert (ex / "repro-half.other_lines").read_text().strip() == "1"  # the count line is counted, never printed
    calls = [c.split() for c in (tmp_path / "a" / "calls").read_text().splitlines()]
    assert sorted(c[2] for c in calls) == sorted(e54.nt.all_tasks())
    assert _parallel(tmp_path / "b", "repro", {"x2": 3, "e75": 1}) == 1
