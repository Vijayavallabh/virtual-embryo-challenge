"""The test-board build launcher of row E51's two embryo candidates, SO (T2embryo_s2xs_ov) and O_hi (T2embryo_s2x_ovhi):
experiments/builds/embryo_so_ohi.sh.

Its gate builds a candidate only when row E51's committed reading states, in that arm's own fields, that it passes;
P's, O's and s2xs's commands must give their committed sha256 (read by code from their committed build records) before
either candidate is built; each candidate is built twice and compared; the identity step prints booleans only. These
tests hold the launcher's strings and step order, run every gate refusal for real in throwaway git repositories under a
temporary home with a stand-in uv that fails (so nothing is built and no competition file is read), and run its gate's
reading check and its identity step on synthetic records and files. The launcher is never run here with a real commit.
Built through experiments/builds/build_test.py."""

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
REL = "experiments/builds/embryo_so_ohi.sh"
LAUNCHER = REPO / REL
BASE = "--method fgw_shape0 --positions floor_exact --expression mixture"
OVERLAY = BASE + " --shape-interp overlay --shape-interp-test-ok"
SERIES = "--method fgw_shape0 --positions floor_exact --expression series"
SO = SERIES + " --shape-interp overlay --shape-interp-test-ok"
OHI = "--method fgw_shape0 --positions floor_exact --expression render_hi --shape-interp overlay --shape-interp-test-ok"
SO_DRY = SERIES + " --shape-interp overlay"
OHI_DRY = "--method fgw_shape0 --positions floor_exact --expression render_hi --shape-interp overlay"
SHAS = {"P": "a03611a343faf1d37edf06ce0e2bccc86cce6f14d3580b823864e5faf3ba03f9",
        "O": "1cfc91be3a20bfdd1a4f7bcfebd7316e1143c6d4048e578c1f30bcfbe1e0e1ea",
        "S": "bc50d012cd2c9481cd18c1f5bcf5b01093e22fd16ac65582a5c9f34321f13115"}
RECORDS = {"P": ("experiments/builds/results/embryo2_2dccfee.SHA256SUMS", "builds/T2embryo_s2x.h5ad"),
           "O": ("experiments/builds/results/embov_3b813b5.SHA256SUMS", "OVa/T2embryo_s2x_ov.h5ad"),
           "S": ("experiments/builds/results/embryo3_3d52b83.SHA256SUMS", "builds/T2embryo_s2xs.h5ad")}
READING = "experiments/coords/results/e51_so_ohi/reading.json"
SO_ARM, OH_ARM = "fgw_shape0x_series_ov", "fgw_shape0x_ov_hi"


# ------------------------------------------------------------------------------------------------ the launcher's text
def test_the_commands_are_the_committed_ones_and_the_specifications():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^([A-Z_]+)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert args == {"BASE": BASE, "OV": OVERLAY, "S": SERIES, "SO": SO, "OH": OHI, "SO_DRY": SO_DRY, "OH_DRY": OHI_DRY}
    # P's and O's commands are embryo_overlay.sh's; s2xs's is set embryo3's of build_all.sh
    ov = (REPO / "experiments/builds/embryo_overlay.sh").read_text()
    assert f'BASE_ARGS="{BASE}"' in ov and f'OV_ARGS="{OVERLAY}"' in ov
    ba = " ".join((REPO / "experiments/builds/build_all.sh").read_text().split())
    assert f"build_test.py embryo {SERIES} --out $O/builds/T2embryo_s2xs.h5ad" in ba
    for k, sha in SHAS.items():
        assert f"{k}_SHA={sha}" in sh and f'[ "$got" = "${k}_SHA" ]' in sh
        sums, entry = RECORDS[k]
        assert f"{k}_SUMS={sums}" in sh and f"{k}_ENTRY={entry}" in sh
    assert "BOARD=T2:embryo:test_interp" in sh and "VAL_BOARD=T2:embryo:val_interp" in sh
    assert '--board "$VAL_BOARD" $SO_DRY_ARGS' in sh and '--board "$VAL_BOARD" $OH_DRY_ARGS' in sh
    assert '"$S/V/T2embryo_s2xs_ov_dry.h5ad"' in sh and '"$S/V/T2embryo_s2x_ovhi_dry.h5ad"' in sh
    assert "build D T2embryo_s2x $BASE_ARGS" in sh and "build O T2embryo_s2x_ov $OV_ARGS" in sh
    assert "build S T2embryo_s2xs $S_ARGS" in sh
    assert 'build "$STEP" T2embryo_s2xs_ov $SO_ARGS' in sh and 'build "$STEP" T2embryo_s2x_ovhi $OH_ARGS' in sh
    assert 'cmp -s "$O/SOa/T2embryo_s2xs_ov.h5ad" "$O/SOb/T2embryo_s2xs_ov.h5ad" || finish 31' in sh
    assert 'cmp -s "$O/OHa/T2embryo_s2x_ovhi.h5ad" "$O/OHb/T2embryo_s2x_ovhi.h5ad" || finish 32' in sh


def test_the_committed_records_hold_the_three_sha256():
    for k, (sums, entry) in RECORDS.items():
        rows = dict(reversed(line.split()) for line in (REPO / sums).read_text().splitlines() if line.strip())
        assert rows[entry] == SHAS[k]


def test_the_steps_run_in_order_each_writing_running_first():
    sh = LAUNCHER.read_text()
    order = []
    for m in re.finditer(r"^STEP=(\w+)$|^for STEP in ([\w ]+); do$", sh, flags=re.M):
        order += [m.group(1)] if m.group(1) else m.group(2).split()
    assert order == ["gate", "tests", "V", "D", "O", "S", "SOa", "SOb", "OHa", "OHb", "identity", "done"]
    assert sh.count("echo running > \"$S/$STEP.exit\"") == 10  # gate, tests, V, D, O, S, the two loops, identity, done
    assert sh.index('echo running > "$S/$STEP.exit"') < sh.index("head=$(git rev-parse HEAD)")
    assert sh.count("echo skipped > \"$S/$STEP.exit\"") == 2
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh
    i8 = sh.index("OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8")
    i16 = sh.index("OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16")
    assert i8 < sh.index("STEP=tests") < i16 < sh.index("STEP=V")
    assert "S=$HOME/v/scratch/embsohi/$RUN" in sh and "O=runs/embsohi_${head:0:7}_$RUN" in sh
    assert '> "$S/identity.out" 2> /dev/null' in sh and "|| finish 40" in sh


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
    assert "tests/test_embryo_so_ohi_launcher.py" in listed and all((REPO / f).exists() for f in listed)


# ------------------------------------------------------------------------------------------- the gate's reading check
def _arm(label, control, passes=True, void=False, clauses=None):
    return {"void": void, "label": label, "control": control, "passes": passes,
            "clauses": clauses if clauses is not None else {"a": True, "b": True, "c": True}}


def _reading(commit, so=None, oh=None, void=False):
    return {"void": void, "commit": commit,
            "arms": {SO_ARM: so or _arm("SO", "fgw_shape0x_series"), OH_ARM: oh or _arm("O_hi", "fgw_shape0x_ov")}}


def _reading_check() -> str:
    m = re.search(r"PASSES=\$\(python3 -I -B -c '([^']+)' \"\$READING\"", LAUNCHER.read_text())
    assert m, "the gate's reading check is not where the test expects it"
    return m.group(1)


def test_the_reading_check_takes_each_arms_own_fields(tmp_path):
    run = "a" * 40
    so_fail = _arm("SO", "fgw_shape0x_series", passes=False)
    cases = [(_reading(run), f"{run} 1 1"),
             (_reading(run, so=so_fail), f"{run} 0 1"),
             (_reading(run, oh=_arm("O_hi", "fgw_shape0x_ov", void=True)), f"{run} 1 0"),
             (_reading(run, oh=_arm("O_hi", "fgw_shape0x_ov", clauses={"a": True, "b": False, "c": True})),
              f"{run} 1 0"),
             (_reading(run, oh=_arm("O_hi", "fgw_shape0x_ov", clauses={"a": True, "b": True})), f"{run} 1 0"),
             (_reading(run, oh=_arm("O_hi", "fgw_shape0x_ov", passes="True")), f"{run} 1 0"),
             (_reading(run, so=_arm("SO", "fgw_shape0x_ov")), f"{run} 0 1"),  # not the row's control
             (_reading(run, so=so_fail, oh=_arm("O_hi", "fgw_shape0x_ov", passes=False)), None),
             (_reading(run, void=True), None), (_reading(run[:7]), None),
             ({k: v for k, v in _reading(run).items() if k != "commit"}, None),
             ({"void": False, "commit": run, "outcome": "SO: passes; O_hi: passes"}, None)]
    for k, (reading, want) in enumerate(cases):
        p = tmp_path / f"reading{k}.json"
        p.write_text(json.dumps(reading))
        got = subprocess.run(["python3", "-I", "-B", "-c", _reading_check(), str(p)], capture_output=True, text=True)
        if want is None:
            assert got.returncode == 1, (reading, got.stdout)
        else:
            assert got.returncode == 0 and got.stdout.strip() == want, (reading, got.stdout, got.stderr)


# ------------------------------------------------------------------------------------------------ the gate, for real
def _checkout(tmp_path: Path, reading=lambda c: _reading(c), row_at_run=True, records=True) -> tuple[Path, str]:
    """A throwaway repository: the launcher, the rules file and its hash, the three committed build records and a ledger
    with row E51 at the reading's run commit; then (unless `reading` is None) row E51's reading made at that commit,
    committed. Returns the repository and the run commit."""
    repo = tmp_path / "checkout"
    for rel in (REL, "docs/S2_RULES.md", "docs/S2_RULES.sha256", *(r[0] for r in RECORDS.values())):
        if not records and rel.endswith(".SHA256SUMS"):
            continue
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, repo / rel)
    (repo / "experiments/harness").mkdir(parents=True)
    tries = repo / "experiments/harness/TRIES.md"
    tries.write_text("| E51 (one row; synthetic) | embryo |\n" if row_at_run else "| E50 (synthetic) | embryo |\n")
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "row E51")
    run = git(repo, "rev-parse", "HEAD")
    if not row_at_run:
        tries.write_text("| E50 (synthetic) | embryo |\n| E51 (one row; synthetic) | embryo |\n")
    if reading is not None:
        p = repo / READING
        p.parent.mkdir(parents=True)
        p.write_text(json.dumps(reading(run)))
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "row E51 read", "--allow-empty")
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
    return home / "v" / "scratch" / "embsohi" / run


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
    assert not (home / "v").exists()
    assert _run(repo, home, "0" * 40, "r0").returncode == 96
    again = _run(repo, home, head, "r0")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert (_scratch(home, "r0") / "gate.exit").read_text().strip() == "96"


def test_an_unclean_tree_or_a_changed_rules_file_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a")
    (repo / "experiments/harness/TRIES.md").write_text("| E51 changed |\n")
    _refused(repo, _home(tmp_path, "h1"), "u1", 98, "tracked files differ")
    repo, _ = _checkout(tmp_path / "b")
    with open(repo / "docs/S2_RULES.md", "a") as f:
        f.write("\nchanged\n")
    git(repo, "commit", "-qam", "rules changed")
    _refused(repo, _home(tmp_path, "h2"), "u2", 97, "")


def test_without_row_e51_or_its_committed_reading_the_gate_refuses(tmp_path):
    repo, _ = _checkout(tmp_path / "a")
    (repo / "experiments/harness/TRIES.md").write_text("| E50 (synthetic) | embryo |\n")
    git(repo, "commit", "-qam", "no row")
    _refused(repo, _home(tmp_path, "h1"), "r1", 95, "row E51 is not in TRIES.md")
    repo, _ = _checkout(tmp_path / "b", reading=None)
    _refused(repo, _home(tmp_path, "h2"), "r2", 94, "no committed reading")
    repo, _ = _checkout(tmp_path / "c")  # the reading exists but is not tracked
    git(repo, "rm", "-q", "--cached", READING)
    git(repo, "commit", "-qm", "untrack the reading")
    _refused(repo, _home(tmp_path, "h3"), "r3", 94, "no committed reading")


def test_a_void_reading_or_one_that_passes_neither_candidate_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a", reading=lambda c: _reading(c, void=True))
    _refused(repo, _home(tmp_path, "h1"), "v1", 92, "states no candidate")
    neither = lambda c: _reading(c, so=_arm("SO", "fgw_shape0x_series", passes=False),  # noqa: E731
                                 oh=_arm("O_hi", "fgw_shape0x_ov", void=True))
    repo, _ = _checkout(tmp_path / "b", reading=neither)
    _refused(repo, _home(tmp_path, "h2"), "v2", 92, "states no candidate")
    repo, _ = _checkout(tmp_path / "c", reading=lambda c: {"void": False, "commit": c, "outcome": "SO: passes"})
    _refused(repo, _home(tmp_path, "h3"), "v3", 92, "states no candidate")


def test_a_foreign_run_commit_or_one_without_row_e51_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a", reading=lambda c: _reading("f" * 40))
    _refused(repo, _home(tmp_path, "h1"), "c1", 88, "not an ancestor")
    repo, _ = _checkout(tmp_path / "b", row_at_run=False)
    _refused(repo, _home(tmp_path, "h2"), "c2", 86, "not in TRIES.md at the reading")


def test_a_missing_or_different_build_record_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path / "a", records=False)
    _refused(repo, _home(tmp_path, "h1"), "s1", 85, "does not give P_SHA")
    repo, _ = _checkout(tmp_path / "b")
    sums = repo / RECORDS["S"][0]
    sums.write_text(sums.read_text().replace(SHAS["S"], "0" * 64))
    git(repo, "commit", "-qam", "record changed")
    _refused(repo, _home(tmp_path, "h2"), "s2", 85, "does not give S_SHA")


def test_an_existing_output_folder_is_refused(tmp_path):
    repo, _ = _checkout(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    (repo / "runs" / f"embsohi_{head[:7]}_o1").mkdir(parents=True)
    out = _run(repo, _home(tmp_path), head, "o1")
    assert out.returncode == 93 and "output folder" in out.stderr


def _exit_files(s: Path) -> dict:
    return {p.name: (p.read_text(), p.stat().st_mtime_ns) for p in s.glob("*.exit")}


def test_a_good_reading_passes_the_gate_and_a_failing_uv_stops_the_tests_step(tmp_path):
    repo, run = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    out = _run(repo, home, head, "g1")
    s = _scratch(home, "g1")
    assert out.returncode == 10, out.stderr
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert sorted(p.name for p in s.glob("*.exit")) == ["gate.exit", "tests.exit"]
    assert (s / "head").read_text().strip() == head
    assert (s / "candidates").read_text() == "SO 1\nO_hi 1\n"
    assert (repo / "runs" / f"embsohi_{head[:7]}_g1").is_dir()
    calls = (home / "uvcalls").read_text().splitlines()
    assert len(calls) == 1 and "pytest" in calls[0] and "tests/test_embryo_so_ohi_launcher.py" in calls[0]
    before = _exit_files(s)
    again = _run(repo, home, head, "g1")  # a refused second call leaves the finished steps' exit files as they were
    assert again.returncode == 91 and _exit_files(s) == before
    assert len((home / "uvcalls").read_text().splitlines()) == 1


def test_a_candidate_whose_pass_the_reading_does_not_state_is_marked_not_built(tmp_path):
    only_so = lambda c: _reading(c, oh=_arm("O_hi", "fgw_shape0x_ov", passes=False))  # noqa: E731
    repo, _ = _checkout(tmp_path, reading=only_so)
    home = _home(tmp_path)
    out = _run(repo, home, git(repo, "rev-parse", "HEAD"), "p1")
    assert out.returncode == 10
    assert (_scratch(home, "p1") / "candidates").read_text() == "SO 1\nO_hi 0\n"


# ------------------------------------------------------------------------------------------------ the identity step
def _identity_code() -> str:
    m = re.search(r"STEP=identity\n.*?\$PY -c '(.+?)' \"\$O/O/", LAUNCHER.read_text(), flags=re.S)
    assert m, "the identity step is not where the test expects it"
    return m.group(1)


def _h5ad(path: Path, X, coords, names=None, genes=None):
    a = anndata.AnnData(X=np.asarray(X, np.float32),
                        obs=pd.DataFrame(index=names or [f"cell_{i}" for i in range(len(X))]),
                        var=pd.DataFrame(index=genes or [f"g{j}" for j in range(X.shape[1])]))
    a.obsm["spatial_3D"] = np.asarray(coords, np.float32)
    a.write_h5ad(path)
    return str(path)


def _files(tmp_path: Path, so_coords=None, so_X=None, oh_coords=None, oh_names=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(3)
    X_mix = rng.gamma(2.0, 1.0, size=(200, 5)).astype(np.float32)
    X_ser = rng.gamma(2.0, 1.0, size=(200, 5)).astype(np.float32)
    X_hi = rng.gamma(3.0, 1.0, size=(200, 5)).astype(np.float32)
    C_ov = rng.normal(size=(200, 3)).astype(np.float32)
    C_ex = rng.normal(size=(200, 3)).astype(np.float32)
    o = _h5ad(tmp_path / "o.h5ad", X_mix, C_ov)
    s = _h5ad(tmp_path / "s.h5ad", X_ser, C_ex)
    so = _h5ad(tmp_path / "so.h5ad", X_ser if so_X is None else so_X, C_ov if so_coords is None else so_coords)
    oh = _h5ad(tmp_path / "oh.h5ad", X_hi, C_ov if oh_coords is None else oh_coords, names=oh_names)
    return o, s, so, oh, C_ov, X_ser


def _identity(*paths):
    return subprocess.run([sys.executable, "-I", "-B", "-c", _identity_code(), *paths], capture_output=True, text=True,
                          timeout=120, env={**os.environ})


def test_the_identity_step_prints_booleans_only_and_holds_on_the_recombined_files(tmp_path):
    o, s, so, oh, _, _ = _files(tmp_path)
    out = _identity(o, s, so, oh)
    assert out.returncode == 0, out.stderr
    lines = out.stdout.strip().splitlines()
    assert [ln.split()[0] for ln in lines] == ["so_coords_equal_o", "so_x_equal_s2xs", "so_obs_names_equal_s2xs",
                                                "so_var_names_equal_s2xs", "ohi_coords_equal_o",
                                                "ohi_obs_names_equal_o", "ohi_var_names_equal_o"]
    assert all(re.fullmatch(r"[a-z0-9_]+ True", ln) for ln in lines)
    only_oh = _identity(o, s, "-", oh)
    assert only_oh.returncode == 0 and [ln.split()[0] for ln in only_oh.stdout.split("\n") if ln][:1] == [
        "ohi_coords_equal_o"]
    assert _identity(o, s, "-", "-").returncode == 1  # nothing checked is not a pass


def test_the_identity_step_fails_on_one_bit_one_value_or_one_name(tmp_path):
    _, _, _, _, C, X = _files(tmp_path / "a")
    nudged = C.copy()
    nudged[7, 2] = np.nextafter(nudged[7, 2], np.float32(np.inf))
    out = _identity(*_files(tmp_path / "b", so_coords=nudged)[:4])
    assert out.returncode == 1 and "so_coords_equal_o False" in out.stdout
    X2 = X.copy()
    X2[0, 0] += 1
    out = _identity(*_files(tmp_path / "c", so_X=X2)[:4])
    assert out.returncode == 1 and "so_x_equal_s2xs False" in out.stdout
    out = _identity(*_files(tmp_path / "d", oh_coords=nudged)[:4])
    assert out.returncode == 1 and "ohi_coords_equal_o False" in out.stdout
    names = [f"cell_{i}" for i in range(200)]
    names[3], names[4] = names[4], names[3]
    out = _identity(*_files(tmp_path / "e", oh_names=names)[:4])
    assert out.returncode == 1 and "ohi_obs_names_equal_o False" in out.stdout
    assert all(re.fullmatch(r"[a-z0-9_]+ (True|False)", ln) for ln in out.stdout.strip().splitlines())
