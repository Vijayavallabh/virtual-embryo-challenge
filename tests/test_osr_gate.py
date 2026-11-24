"""The OSR gate (experiments/coords/osr_gate.py) and its launcher: the per-seed computation on synthetic brackets
against a hand computation from vec.shape_interp, the stop clauses S1 to S5 and PASS on hand-made pair files, the pair
command on a task with its target absent (the loaders swapped for synthetic stages that record every path asked for),
the refusals, and the launcher's syntax and refusal exits. No released data, no DGX."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "coords" / "osr_gate.sh"
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import osr_gate as og  # noqa: E402

from vec import shape_interp as si  # noqa: E402
from vec import t2_build as tb  # noqa: E402
from vec.boards import Board  # noqa: E402
from vec.io import COORD_KEY  # noqa: E402

GENES = [f"g{i}" for i in range(12)]


def _cloud(sd, n, seed):
    return np.random.default_rng(seed).normal(size=(n, 3)) * np.asarray(sd, float) + 5.0


def _expr(n, seed):
    return np.random.default_rng(seed).gamma(1.0, 1.0, (n, len(GENES))).astype(np.float32)


def _hand(Cb, Ca, Eb, Ea, t, n, seed):
    res = tb.interpolate_clouds(Cb, Ca, Eb, Ea, t, n, seed, target_rms=1.0, t_pos=0.0, positions="floor_exact",
                                shape_form="overlay")
    U = np.asarray(res.coords, np.float64)
    Zb = (Cb - Cb.mean(0)) / tb.rms_radius(Cb)
    Za = (Ca - Ca.mean(0)) / tb.rms_radius(Ca)
    blocks = si.bracket_blocks(Zb, Za)
    s_o, Vt = si.principal(U)
    f = si.stretch_factors(s_o, si.log_interp_sds(si.normalised_sds(Zb), si.normalised_sds(Za), t), blocks)
    fa = si.stretch_factors(s_o, si.log_interp_sds(si.normalised_sds(Zb), si.normalised_sds(Za), t, anti=True),
                            blocks)
    return U, blocks, f, fa, Vt


# ------------------------------------------------------------------------------------------------ the computation
def test_g_on_a_1_23_bracket_is_log_f1_over_f2():
    Cb, Ca, Eb, Ea = _cloud((4, 1, 1), 600, 1), _cloud((5, 1.2, 1.2), 500, 2), _expr(600, 3), _expr(500, 4)
    q = og.seed_quantities(Cb, Ca, Eb, Ea, 0.5, 300, 1_000_001)
    U, blocks, f, fa, Vt = _hand(Cb, Ca, Eb, Ea, 0.5, 300, 1_000_001)
    assert si.blocks_text(blocks) == "1|2,3" == q["blocks"]
    assert q["g"] == pytest.approx(np.log(f[0] / f[1]), abs=1e-12)
    assert q["g_anti"] == pytest.approx(np.log(fa[0] / fa[1]), abs=1e-12)
    assert q["thin_sd_union"] == pytest.approx(si.normalised_sds(U)[-1], abs=1e-12)
    m = U.mean(0)
    assert q["thin_sd_stretched"] == pytest.approx(
        si.normalised_sds(m + (U - m) @ si.stretch_matrix(Vt, f))[-1], abs=1e-12)
    assert q["angle_axis1_z_deg"] == pytest.approx(np.degrees(np.arccos(abs(Vt[0][2]))), abs=1e-9)
    assert 0.0 <= q["angle_axis1_z_deg"] <= 90.0
    assert "n" not in q  # no count is recorded


def test_other_blocks_take_the_max_log_ratio_with_the_sign_of_f1_over_f3():
    Cb, Ca, Eb, Ea = _cloud((6, 3, 1), 600, 1), _cloud((7, 3, 1), 500, 2), _expr(600, 3), _expr(500, 4)
    q = og.seed_quantities(Cb, Ca, Eb, Ea, 0.5, 300, 1_000_002)
    _, blocks, f, fa, _ = _hand(Cb, Ca, Eb, Ea, 0.5, 300, 1_000_002)
    assert si.blocks_text(blocks) == q["blocks"] != "1|2,3"
    lf = np.log(f)
    want = np.sign(lf[0] - lf[2]) * max(abs(lf[j] - lf[k]) for j in range(3) for k in range(3))
    assert q["g"] == pytest.approx(want, abs=1e-12)
    assert q["g"] != pytest.approx(lf[0] - lf[1], abs=1e-6)  # the other branch, not log(f1/f2)


def test_g_of_by_hand():
    f = [1.2, 1.0, 0.8]
    assert og.g_of(f, ((0,), (1, 2))) == pytest.approx(np.log(1.2))
    assert og.g_of(f, ((0, 1), (2,))) == pytest.approx(np.log(1.2 / 0.8))
    assert og.g_of([0.8, 1.0, 1.2], ((0,), (1,), (2,))) == pytest.approx(-np.log(1.5))
    assert og.g_of([1.0, 1.0, 1.0], ((0, 1, 2),)) == 0.0


# ------------------------------------------------------------------------------------------------ the reading
def _pair(out: Path, name: str, gs, blocks="1|2,3", errors=()):
    per = {}
    for s, g in zip(og.SEEDS, gs, strict=True):
        per[str(s)] = {"error": "ShapeRefused: x"} if s in errors else {"blocks": blocks, "g": float(g), "g_anti": 0.0}
    ok = [v["g"] for v in per.values() if "g" in v]
    (out / f"{name}.json").write_text(json.dumps({"name": name, "commit": "c" * 40, "per_seed": per,
                                                  "median_g": float(np.median(ok)) if ok else None}))


def _world(out: Path, test=(0.2,) * 10, like=(0.1, 0.2, 0.3, 0.25), test_blocks="1|2,3", skip=()):
    if "test" not in skip:
        _pair(out, "test", test, blocks=test_blocks)
    for n, m in zip(og.TEST_LIKE, like, strict=True):
        _pair(out, n, [m] * 10)
    for n in ("e75", "e725", "x1", "x2", "x3"):
        if n not in skip:
            _pair(out, n, [-0.4] * 10, blocks="1|2|3")


def test_pass(tmp_path):
    _world(tmp_path)
    r = og.reading(tmp_path)
    assert r["verdict"] == "PASS" and r["stopped_by"] == []
    assert r["median_g"]["test"] == pytest.approx(0.2)
    assert r["reported_only"]["test_abs_median_g_inside_test_like_range"] is True
    assert r["reported_only"]["statement"] == og.REPORTED_ONLY


def test_reported_only_outside_the_range_decides_nothing(tmp_path):
    _world(tmp_path, test=(0.9,) * 10)
    r = og.reading(tmp_path)
    assert r["verdict"] == "PASS" and r["reported_only"]["test_abs_median_g_inside_test_like_range"] is False


def test_s1(tmp_path):
    _world(tmp_path, test=(0.01,) * 8 + (0.3, 0.3))
    assert og.reading(tmp_path)["stopped_by"] == ["S1"]
    _world(tmp_path, test=(0.01,) * 7 + (0.3,) * 3)  # seven small seeds: not S1
    assert "S1" not in og.reading(tmp_path)["stopped_by"]


def test_s2(tmp_path):
    _world(tmp_path, test=(0.2,) * 7 + (-0.2,) * 3)
    assert og.reading(tmp_path)["stopped_by"] == ["S2"]
    _world(tmp_path, test=(0.2,) * 8 + (-0.2,) * 2)  # two against the majority: not S2
    assert og.reading(tmp_path)["verdict"] == "PASS"


def test_s3(tmp_path):
    _world(tmp_path, test_blocks="1|2|3")
    assert og.reading(tmp_path)["stopped_by"] == ["S3"]


def test_s4_opposite_sign(tmp_path):
    _world(tmp_path, like=(-0.1, -0.2, -0.3, 0.25))
    assert og.reading(tmp_path)["stopped_by"] == ["S4"]
    _world(tmp_path, like=(-0.1, -0.2, 0.3, 0.25))  # two opposite: not S4
    assert og.reading(tmp_path)["verdict"] == "PASS"


def test_s4_practice_cannot_host_the_step(tmp_path):
    _world(tmp_path, like=(0.01, -0.02, 0.04, 0.0))
    assert og.reading(tmp_path)["stopped_by"] == ["S4"]


def test_s5_seed_error_missing_file_and_pair_error(tmp_path):
    _world(tmp_path)
    _pair(tmp_path, "x2", [0.1] * 10, errors=(og.SEEDS[3],))
    assert og.reading(tmp_path)["stopped_by"] == ["S5"]
    _world(tmp_path)
    (tmp_path / "x3.json").unlink()
    assert og.reading(tmp_path)["stopped_by"] == ["S5"]
    _world(tmp_path)
    (tmp_path / "e75.json").write_text(json.dumps({"name": "e75", "error": "ValueError: x"}))
    assert og.reading(tmp_path)["stopped_by"] == ["S5"]


def test_read_command_writes_both_files(tmp_path):
    _world(tmp_path)
    assert og.main(["read", "--out", str(tmp_path)]) == 0
    assert json.loads((tmp_path / "gate.json").read_text())["verdict"] == "PASS"
    assert (tmp_path / "gate.txt").read_text().startswith("OSR gate: PASS")


# ------------------------------------------------------------------------------------------------ the pair command
def _stage(sd, n, seed):
    return ad.AnnData(X=_expr(n, seed + 100), obsm={COORD_KEY: _cloud(sd, n, seed)}, var=pd.DataFrame(index=GENES))


@pytest.mark.parametrize("name", ["half", "e75"])
def test_pair_for_a_task_never_asks_for_its_target(tmp_path, monkeypatch, name):
    task = og.nt.all_tasks()[name]
    data = tmp_path / "data"
    data.mkdir()
    target = data / task.target
    assert not target.exists()  # absent: any read or stat of it would fail
    asked = []

    def load(path, genes):
        asked.append(Path(path))
        assert Path(path) != target
        return _stage((4, 1, 1), 400, len(asked))

    def rows(path, genes):
        asked.append(Path(path))
        assert Path(path) != target
        return 400

    def no_read(*a, **k):
        raise AssertionError("a file was read outside the gate's loaders")

    board = Board(key="T2:embryo:val_interp", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=10, max_cells=300)
    monkeypatch.setattr(og, "_load", load)
    monkeypatch.setattr(og, "_rows_as_released", rows)
    monkeypatch.setattr(og, "_boards", lambda: {"T2:embryo:val_interp": board})
    monkeypatch.setattr(og, "commit_problems", lambda commit, repo=og.REPO: [])
    monkeypatch.setattr(og, "SEEDS", og.SEEDS[:2])
    monkeypatch.setattr(ad, "read_h5ad", no_read)
    out = tmp_path / "out"
    assert og.main(["pair", "--name", name, "--data", str(data), "--out", str(out), "--commit", "c" * 40]) == 0
    assert (out / f"{name}.exit").read_text().strip() == "0"
    rec = json.loads((out / f"{name}.json").read_text())
    assert set(rec["per_seed"]) == {str(s) for s in og.SEEDS[:2]} and rec["median_g"] is not None
    assert all("g" in v for v in rec["per_seed"].values())
    assert target not in asked and {p.relative_to(data).as_posix() for p in asked} == {
        task.inputs["before"], task.inputs["after"]}
    days = og.harness.TaskView.of(task)
    lo, hi = days.role_days["before"], days.role_days["after"]
    assert rec["t"] == pytest.approx((days.target_day - lo) / (hi - lo))
    assert not any("n" in v for v in rec["per_seed"].values()) and "n" not in rec  # no count is recorded


def test_a_seed_that_raises_is_recorded_and_the_files_are_written(tmp_path, monkeypatch):
    board = Board(key="T2:embryo:val_interp", task="T2", genes=list(GENES), needs_coords=True, min_cells=10,
                  max_cells=300)
    monkeypatch.setattr(og, "_load", lambda path, genes: _stage((4, 1, 1), 400, 1))
    monkeypatch.setattr(og, "_boards", lambda: {"T2:embryo:val_interp": board})
    monkeypatch.setattr(og, "commit_problems", lambda commit, repo=og.REPO: [])
    monkeypatch.setattr(og, "SEEDS", og.SEEDS[:2])

    def boom(*a, **k):
        raise si.ShapeRefused("refused here")

    monkeypatch.setattr(og, "seed_quantities", boom)
    out = tmp_path / "out"
    assert og.main(["pair", "--name", "test", "--data", str(tmp_path), "--out", str(out), "--commit", "c" * 40]) == 3
    assert (out / "test.exit").read_text().strip() == "3"
    rec = json.loads((out / "test.json").read_text())
    assert all("ShapeRefused" in v["error"] for v in rec["per_seed"].values())


def test_refused_at_another_commit(tmp_path):
    out = tmp_path / "out"
    assert og.main(["pair", "--name", "test", "--data", str(tmp_path), "--out", str(out), "--commit", "0" * 40]) == 2
    assert (out / "test.exit").read_text().strip() == "2"
    assert not (out / "test.json").exists()


def test_refused_when_the_pair_file_exists(tmp_path):
    (tmp_path / "half.json").write_text("{}")
    (tmp_path / "half.exit").write_text("0\n")
    assert og.main(["pair", "--name", "half", "--out", str(tmp_path), "--commit", "0" * 40]) == 2
    assert (tmp_path / "half.exit").read_text() == "0\n" and (tmp_path / "half.json").read_text() == "{}"


def test_unknown_name_is_refused(tmp_path):
    with pytest.raises(SystemExit) as e:
        og.main(["pair", "--name", "E8.0", "--out", str(tmp_path), "--commit", "0" * 40])
    assert e.value.code == 2


def test_commit_problems_on_a_scratch_repository(tmp_path, monkeypatch):
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path.parent))
    for k in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE"):
        monkeypatch.delenv(k, raising=False)
    repo = tmp_path / "r"
    repo.mkdir()
    git = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false"]
    subprocess.run([*git, "init", "-q"], check=True)
    (repo / "a").write_text("1")
    subprocess.run([*git, "add", "a"], check=True)
    subprocess.run([*git, "commit", "-q", "-m", "a"], check=True)
    head = subprocess.run([*git, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (repo / "untracked").write_text("x")  # untracked files do not count
    assert og.commit_problems(head, repo) == []
    assert any("HEAD" in w for w in og.commit_problems("0" * 40, repo))
    (repo / "a").write_text("2")
    assert og.commit_problems(head, repo) == ["the tracked tree is not clean"]


# ------------------------------------------------------------------------------------------------ the launcher
def _launch(home: Path, *argv: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(home)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=REPO, env=env, capture_output=True, text=True)


def test_launcher_syntax_and_contents():
    assert subprocess.run(["bash", "-n", str(LAUNCHER)]).returncode == 0
    sh = LAUNCHER.read_text()
    assert not re.search(r"(^|[\s;|&(])rm\s", sh)
    assert re.search(r"^NAMES=\(([^)]*)\)$", sh, flags=re.M).group(1).split() == list(og.NAMES)
    assert "OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh and "> /dev/null 2>&1" in sh
    assert 'PATH="$HOME/.local/bin:$PATH"' in sh and 'TMPDIR="$S/tmp"' in sh
    assert "upload" not in sh.split("\nset -u", 1)[1]


@pytest.mark.parametrize("run", ["a/b", "..", ".x", "a b", "x;y", "$(id)"])
def test_launcher_refuses_an_odd_run_name(tmp_path, run):
    r = _launch(tmp_path, "c" * 40, run)
    assert r.returncode == 90
    assert not (tmp_path / "v").exists()


def test_launcher_refuses_an_existing_run_folder(tmp_path):
    s = tmp_path / "v" / "scratch" / "osr_gate" / "r1"
    s.mkdir(parents=True)
    r = _launch(tmp_path, "c" * 40, "r1")
    assert r.returncode == 91 and list(s.iterdir()) == []


def test_launcher_usage(tmp_path):
    assert _launch(tmp_path, "c" * 40).returncode == 2


FAKE_UV = """#!/usr/bin/env bash
# stands in for uv: records its arguments and plays the pair or read command's exit files
echo "$*" >> "$HOME/calls"
args=("$@")
for i in "${!args[@]}"; do
  case "${args[$i]}" in
    --name) name="${args[$((i + 1))]}" ;;
    --out) out="${args[$((i + 1))]}" ;;
  esac
done
if [[ " $* " == *" pair "* ]]; then
  code=0
  if [ "$name" = "${FAIL_PAIR:-}" ]; then code=3; fi
  echo "$code" > "$out/$name.exit"
  exit "$code"
fi
exit 0
"""


@pytest.mark.parametrize("fail, want", [("", "0"), ("x2", "4")])
def test_launcher_runs_ten_pairs_then_read(tmp_path, fail, want):
    bin_ = tmp_path / ".local" / "bin"
    bin_.mkdir(parents=True)
    (bin_ / "uv").write_text(FAKE_UV)
    (bin_ / "uv").chmod(0o755)
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env.update(HOME=str(tmp_path), FAIL_PAIR=fail)
    r = subprocess.run(["bash", str(LAUNCHER), "c" * 40, "r1"], cwd=REPO, env=env, capture_output=True, text=True)
    s = tmp_path / "v" / "scratch" / "osr_gate" / "r1"
    assert (s / "done.exit").read_text().strip() == want and str(r.returncode) == want
    calls = (tmp_path / "calls").read_text().splitlines()
    pairs = [c for c in calls if " pair " in c]
    assert len(pairs) == 10 and calls[-1].endswith(f"read --out {s}")
    assert all(c.startswith("run --frozen --no-sync python -I -B experiments/coords/osr_gate.py pair --name ")
               and c.endswith(f"--out {s} --commit {'c' * 40}") for c in pairs)
    assert sorted(c.split("--name ")[1].split()[0] for c in pairs) == sorted(og.NAMES)
    assert all((s / f"{n}.rc").read_text().strip() == ("3" if n == fail else "0") for n in og.NAMES)
