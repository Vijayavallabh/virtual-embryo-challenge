"""The test-board build launcher of the embryo's copy-rows file at the team's pre-registered size, C_s (practice row
E52; the captain's order of 2026-11-21, docs/PLAN.md decision log 2026-11-21 (v)).

experiments/builds/embryo_copy_rows_sized.sh: the gate refuses unless row E52's committed reading is not void, holds
and was made at the given run commit, with the builder's code unchanged since; P's command must give P's known bytes,
C's command C's committed sha256 (read by code from C's committed SHA256SUMS); C_s is built twice and compared; the
identity step prints booleans only. These tests hold the launcher to that, run its gate's refusal paths for real in a
temporary repository under a temporary home, and run its identity step on synthetic files. The launcher is never run
here with the real commit (that would build files). Built through experiments/builds/build_test.py."""

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

from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "embryo_copy_rows_sized.sh"
BASE = "--method fgw_shape0 --positions floor_exact --expression mixture"
COPY_ROWS = BASE + " --copy-rows --copy-rows-test-ok"
SIZED = BASE + " --copy-rows --copy-rows-size --copy-rows-test-ok --copy-rows-size-test-ok"
POSITIONS_SHA = "a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9"
CR_SUMS = "experiments/builds/results/embcr_59ad7fa.SHA256SUMS"


def test_the_launcher_builds_p_then_c_then_c_s_twice_and_checks_identity():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^(BASE|CR|CS)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"BASE": BASE, "CR": COPY_ROWS, "CS": SIZED}
    assert "build D T2embryo_s2x $BASE_ARGS" in sh and f"POSITIONS_SHA={POSITIONS_SHA}" in sh
    assert '[ "$got" = "$POSITIONS_SHA" ]' in sh
    assert "build CR T2embryo_s2x_cr $CR_ARGS" in sh and '[ "$got" = "$CR_SHA" ]' in sh
    assert f"CR_SUMS={CR_SUMS}" in sh and "CR_ENTRY=CRa/T2embryo_s2x_cr.h5ad" in sh
    assert 'build "$STEP" T2embryo_s2x_crs $CS_ARGS' in sh and "for STEP in CSa CSb" in sh
    assert 'cmp -s "$O/CSa/T2embryo_s2x_crs.h5ad" "$O/CSb/T2embryo_s2x_crs.h5ad"' in sh
    assert "STEP=identity" in sh and '> "$S/identity.out" 2> /dev/null || finish 40' in sh
    assert "vec check" in sh and "BOARD=T2:embryo:test_interp" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "git status --porcelain" in sh and "S2_RULES.sha256" in sh and "'^| E52 '" in sh
    assert "READING=experiments/coords/results/e52_copy_rows_sized/reading.json" in sh
    assert 'git ls-files --error-unmatch "$READING"' in sh
    assert 'git diff --quiet "$E52_RUN_COMMIT" HEAD -- $CODE_PATHS' in sh
    assert 'git merge-base --is-ancestor "$E52_RUN_COMMIT" HEAD' in sh
    paths = re.search(r'^CODE_PATHS="([^"]+)"$', sh, flags=re.M).group(1).split()
    assert {"src", "experiments/builds/build_test.py", "uv.lock", "pyproject.toml", "boards", "configs",
            "experiments/t1_m1", "experiments/t1_hygiene", "experiments/t1_strong_null",
            "experiments/density/results/release/count_decision.json", "experiments/coords/results/p3"} == set(paths)
    assert all((REPO / q).exists() for q in paths)
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    for line in sh.splitlines():  # a test build is read by its exit code only
        if re.search(r'build_test\.py embryo|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert "tests/test_embryo_copy_rows_sized_launcher.py" in listed and all((REPO / f).exists() for f in listed)


def test_cs_committed_sha_is_there_for_the_cr_step():
    rows = dict(reversed(line.split()) for line in (REPO / CR_SUMS).read_text().splitlines() if line.strip())
    assert re.fullmatch(r"[0-9a-f]{64}", rows["CRa/T2embryo_s2x_cr.h5ad"])
    assert rows["CRa/T2embryo_s2x_cr.h5ad"] == rows["CRb/T2embryo_s2x_cr.h5ad"]


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    text = LAUNCHER.read_text()
    code = [line.split("#")[0] for line in text.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    for msg in re.findall(r"\$\{\d+:\?([^}]*)\}", text):  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0


def _reading_check() -> str:
    m = re.search(r"python3 -I -B -c '([^']+)' \"\$READING\" \"\$E52_RUN_COMMIT\"", LAUNCHER.read_text())
    assert m, "the gate's reading check is not where the test expects it"
    return m.group(1)


def test_the_gate_takes_only_a_reading_that_is_not_void_holds_and_is_the_arms(tmp_path):
    run = "a" * 40
    good = {"void": False, "build_test_file": True, "i": {"holds": True}, "arm": "fgw_copyrows_sized", "commit": run}
    cases = [(good, 0), ({**good, "void": True}, 1), ({**good, "build_test_file": False}, 1),
             ({**good, "i": {"holds": False}}, 1), ({**good, "arm": "fgw_copyrows"}, 1), ({**good, "i": {}}, 1),
             ({**good, "commit": run[:7]}, 1), ({k: v for k, v in good.items() if k != "commit"}, 1)]
    for k, (reading, want) in enumerate(cases):
        p = tmp_path / f"reading{k}.json"
        p.write_text(json.dumps(reading))
        got = subprocess.run(["python3", "-I", "-B", "-c", _reading_check(), str(p), run], capture_output=True)
        assert got.returncode == want, (reading, got.returncode)


# ------------------------------------------------------------------------------------------------ the gate, for real
def _fake_checkout(tmp_path: Path, reading: dict | None, sums: bool = False) -> tuple[Path, str]:
    """A temporary repository holding the launcher, the rules file and its hash, a ledger with row E52 and (if given)
    row E52's reading, committed; returns it and its HEAD."""
    repo = tmp_path / "checkout"
    for rel in ("experiments/builds/embryo_copy_rows_sized.sh", "docs/S2_RULES.md", "docs/S2_RULES.sha256"):
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, repo / rel)
    (repo / "experiments/harness").mkdir(parents=True)
    (repo / "experiments/harness/TRIES.md").write_text("| E52 (one try; synthetic) | embryo |\n")
    if sums:
        (repo / CR_SUMS).write_text((REPO / CR_SUMS).read_text())
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "row E52")
    head = git(repo, "rev-parse", "HEAD")
    if reading is not None:
        p = repo / "experiments/coords/results/e52_copy_rows_sized/reading.json"
        p.parent.mkdir(parents=True)
        p.write_text(json.dumps({**reading, "commit": reading.get("commit", head)}))
        git(repo, "add", "-A")
        git(repo, "commit", "-qm", "row E52 read")
    return repo, head


def _run(repo: Path, home: Path, *argv: str, trace: bool = False) -> subprocess.CompletedProcess:
    env = {**git_env(), "HOME": str(home)}
    return subprocess.run(["bash", *(["-x"] if trace else []), "experiments/builds/embryo_copy_rows_sized.sh", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, timeout=60)


GOOD = {"void": False, "build_test_file": True, "i": {"holds": True}, "arm": "fgw_copyrows_sized"}


def test_without_e52s_committed_reading_the_gate_refuses_and_builds_nothing(tmp_path):
    repo, run = _fake_checkout(tmp_path, None)
    head = git(repo, "rev-parse", "HEAD")
    out = _run(repo, tmp_path / "home", head, "g1", run)
    scratch = tmp_path / "home" / "v" / "scratch" / "embcrs" / "g1"
    assert out.returncode == 94 and "no committed reading" in out.stderr
    assert (scratch / "gate.exit").read_text().strip() == "94"
    assert sorted(p.name for p in scratch.iterdir()) == ["gate.exit", "tmp"] and not (repo / "runs").exists()


def test_a_void_or_foreign_reading_or_a_short_run_commit_is_refused(tmp_path):
    repo, run = _fake_checkout(tmp_path, {**GOOD, "void": True})
    head = git(repo, "rev-parse", "HEAD")
    assert _run(repo, tmp_path / "h1", head, "g2", run).returncode == 92
    assert _run(repo, tmp_path / "h2", head, "g3", run[:12]).returncode == 88
    assert _run(repo, tmp_path / "h3", head, "g4", "f" * 40).returncode == 87  # not an ancestor
    repo2, run2 = _fake_checkout(tmp_path / "b", {**GOOD, "commit": "0" * 40})
    assert _run(repo2, tmp_path / "h4", git(repo2, "rev-parse", "HEAD"), "g5", run2).returncode == 92


def test_a_good_reading_passes_to_the_c_sums_check_and_stops_there_without_them(tmp_path):
    repo, run = _fake_checkout(tmp_path, GOOD)
    head = git(repo, "rev-parse", "HEAD")
    out = _run(repo, tmp_path / "home", head, "g6", run)
    assert out.returncode == 85 and "no committed sha256 of C" in out.stderr and not (repo / "runs").exists()


def test_a_wrong_commit_is_refused_and_the_trace_of_that_path_stops_before_python(tmp_path):
    repo, run = _fake_checkout(tmp_path, GOOD)
    out = _run(repo, tmp_path / "home", "0" * 40, "t1", run, trace=True)
    assert out.returncode == 96 and "HEAD is not" in out.stderr and "finish 96" in out.stderr
    assert not re.search(r"^\++ (uv|python3?) ", out.stderr, flags=re.M)  # nothing ran but git and the shell
    again = _run(repo, tmp_path / "home", "0" * 40, "t1", run)
    assert again.returncode == 91 and "never reused" in again.stderr
    odd = _run(repo, tmp_path / "home2", "0" * 40, "a/b", run)
    assert odd.returncode == 90 and not (tmp_path / "home2" / "v").exists()


# ------------------------------------------------------------------------------------------------ the identity step
def _identity_code() -> str:
    m = re.search(r"STEP=identity\n\$PY -c '(.+?)' \"\$O/CR", LAUNCHER.read_text(), flags=re.S)
    assert m, "the identity step is not where the test expects it"
    return m.group(1)


def _h5ad(path: Path, X, coords, names=None):
    a = anndata.AnnData(X=np.asarray(X, np.float32),
                        obs=pd.DataFrame(index=names or [f"cell_{i}" for i in range(len(X))]),
                        var=pd.DataFrame(index=[f"g{j}" for j in range(X.shape[1])]))
    a.obsm["spatial_3D"] = np.asarray(coords, np.float32)
    a.write_h5ad(path)
    return path


def _identity(tmp_path, cs_coords=None, cs_X=None, size=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(5)
    X = rng.gamma(2.0, 1.0, size=(300, 5)).astype(np.float32)
    C = (rng.normal(size=(300, 3)) * [3.0, 1.5, 0.5] + 9.0).astype(np.float32)
    target = 2.75
    want = tb.rescale_to_rms(C.astype(np.float64), target).astype(np.float32)
    cr = _h5ad(tmp_path / "cr.h5ad", X, C)
    cs = _h5ad(tmp_path / "cs.h5ad", X if cs_X is None else cs_X, want if cs_coords is None else cs_coords)
    sj = tmp_path / "size.json"
    sj.write_text(json.dumps(size or {"board": "embryo", "valid": True,
                                      "preregistered": {"z_mode": "raw", "pred_rms_full": target}}))
    out = subprocess.run([sys.executable, "-I", "-B", "-c", _identity_code(), str(cr), str(cs), str(sj)],
                         capture_output=True, text=True, timeout=120, env={**os.environ})
    return out, want, X


def test_the_identity_step_prints_booleans_only_and_holds_on_the_rescaled_copy(tmp_path):
    out, _, _ = _identity(tmp_path)
    assert out.returncode == 0, out.stderr
    lines = out.stdout.strip().splitlines()
    assert [ln.split()[0] for ln in lines] == ["size_json_valid_raw", "x_equal", "obs_names_equal", "var_names_equal",
                                                "coords_equal_rescaled_c"]
    assert all(re.fullmatch(r"[a-z_]+ (True|False)", ln) for ln in lines) and all(ln.endswith("True") for ln in lines)


def test_the_identity_step_fails_on_one_bit_or_one_x_value(tmp_path):
    _, want, X = _identity(tmp_path / "a")
    nudged = want.copy()
    nudged[4, 0] = np.nextafter(nudged[4, 0], np.float32(np.inf))
    (tmp_path / "b").mkdir()
    out, _, _ = _identity(tmp_path / "b", cs_coords=nudged)
    assert out.returncode == 1 and "coords_equal_rescaled_c False" in out.stdout
    X2 = X.copy()
    X2[0, 0] += 1
    (tmp_path / "c").mkdir()
    out, _, _ = _identity(tmp_path / "c", cs_X=X2)
    assert out.returncode == 1 and "x_equal False" in out.stdout
    assert all(re.fullmatch(r"[a-z_]+ (True|False)", ln) for ln in out.stdout.strip().splitlines())


def test_the_identity_step_refuses_a_size_json_that_is_not_valid_raw(tmp_path):
    out, _, _ = _identity(tmp_path, size={"board": "embryo", "valid": False,
                                          "preregistered": {"z_mode": "raw", "pred_rms_full": 2.75}})
    assert out.returncode == 1 and "size_json_valid_raw False" in out.stdout
