"""Try E48's run script (experiments/coords/e48_row_order_read.py) and its launcher, on synthetic arrays, hand-made
result records and the committed E47 score files only. Nothing is read from data/, nothing is built from a stage and
nothing is scored.

Pinned: the frozen-code check ignores docstrings and sees code; the wrapped arm re-orders rows only, against the floor
copy score_seeds draws, at the scoring seed, and leaves the file as built when the option refuses; the pre-run gate
passes E47's own scores and fails a seed moved by more than 0.5; the reading applies row E48's void conditions and
readings (i) to (iii), with a refused seed counted as the control in the primary figure and left out of the secondary.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec import evaluate as ev
from vec import harness
from vec import t2_build as tb
from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
SEEDS = list(range(10))
COMMIT = "a" * 40


@pytest.fixture(scope="module")
def e48():
    spec = importlib.util.spec_from_file_location("e48_row_order_read",
                                                  REPO / "experiments" / "coords" / "e48_row_order_read.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _cloud(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n, 3)) * [3.0, 1.6, 0.7]
    Z[: n // 5] = [3.5, 2.2, 1.0] + rng.normal(size=(n // 5, 3)) * 0.4
    return (Z[rng.permutation(n)] + [10.0, -4.0, 2.0]).astype(np.float32)


# ------------------------------------------------------------------------------------------------ checks
def test_the_frozen_code_check_ignores_docstrings_and_sees_code(e48):
    old = 'def f(x):\n    """Old words."""\n    return x + 1\n'
    new = 'def f(x):\n    """New words, longer."""\n    return x + 1\n'
    assert e48.strip_docstrings(old) == e48.strip_docstrings(new)
    assert e48.strip_docstrings(old) != e48.strip_docstrings('def f(x):\n    """Old words."""\n    return x + 2\n')
    same = e48.frozen_changes(lambda p: old, lambda p: old.replace("Old", "New"))
    assert same == []
    assert e48.frozen_changes(lambda p: old, lambda p: old.replace("+ 1", "+ 2")) == list(e48.FROZEN)
    assert e48.frozen_changes(lambda p: "", lambda p: old) == list(e48.FROZEN)  # git show failed: not equal


def test_the_frozen_files_equal_616984b_apart_from_docstrings_here(e48):
    """At this branch the frozen files differ from 616984b in docstrings only (the reviewer's sentences)."""
    if e48._git("cat-file", "-e", f"{e48.FROZEN_AT}^{{commit}}").returncode != 0:
        pytest.skip("616984b is not in this repository")
    if (REPO / "src" / "vec" / "shape_interp.py").exists():
        # a branch that carries row E49's shape option (campaign/embryo-shape-interp) changes the code of t2_build.py
        # and harness.py, so E48's run refuses there by design: E48 runs at its own commit (9cd4b6e)
        assert e48.frozen_changes() == ["src/vec/t2_build.py", "src/vec/harness.py"]
        return
    assert e48.frozen_changes() == []


def test_a_run_is_refused_without_row_e48(e48, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| E47a and E47b (one row) | embryo |\n")
    monkeypatch.setattr(e48, "TRIES", tries)
    with pytest.raises(SystemExit, match="E48 has no row"):
        e48._refuse_unless_ready(COMMIT)


def test_the_constants_name_the_row_its_arms_and_its_bases(e48):
    assert e48.ROW == "E48" and e48.ARMS == {"fgw_shape0x_series_cp": "fgw_shape0x_series",
                                             "fgw_shape0x_cp": "fgw_shape0x"}
    assert e48.BUILDS == ("main", "prefix") and e48.BASES == ("1420c70", "616984b")
    assert (e48.GUARD_MEAN, e48.GUARD_SEED) == (0.1, 0.5) and e48.E47_COMMIT == "9967cf1"
    assert sorted(e48.nt.all_tasks()) == sorted(["half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"])


# ------------------------------------------------------------------------------------------------ the arms
def _setup(e48, n_ref: int = 6000, n_file: int = 1800, truth: int = 250, s: int = 3):
    ref = ev.Stage(X=np.zeros((n_ref, 2), np.float32), coords=_cloud(n_ref, 1))
    board = Board(key="pseudo:T2e/x", task="T2", genes=["a", "b"], needs_coords=True, split="pseudo", min_cells=10,
                  max_cells=n_file, truth_cells=truth)
    rng = np.random.default_rng(9)
    rows = rng.permutation(n_ref)[:n_file]
    X = rng.gamma(2.0, 1.0, size=(n_file, 2)).astype(np.float32)
    base = harness.Candidate("fgw_shape0x_series", frozenset({"T2"}), (),
                             lambda ctx: ev.Stage(X=X, coords=ref.coords[rows]), {"what": "synthetic"})
    ctx = SimpleNamespace(seed=harness.BUILD_SEED_OFFSET + s)
    return ref, board, base, ctx, s


def test_the_arm_reorders_rows_only_against_the_scored_floor_copy(e48):
    ref, board, base, ctx, s = _setup(e48)
    with e48.captured() as seen:
        seen.update(board=board, ref=ref)
        floor = harness.subsample(ref, ev.SUBSAMPLE_FRAC, s)  # what score_seeds draws as the floor at seed s
        ctl = e48.observed(base, seen).build(ctx)
        arm = e48.reordered("fgw_shape0x_series_cp", base, seen)
        out = arm.build(ctx)
    assert harness.subsample is ev.subsample  # the spy is removed afterwards
    rec = seen["arm"][("fgw_shape0x_series_cp", s)]
    assert rec["copy_is_floor"] and rec["scoring_seed_known"] and rec["metric_seed"] == s
    assert not rec["refused"] and rec["permutation"] and rec["n_true_handed"] == board.truth_cells
    assert seen["digest"][("fgw_shape0x_series_cp", s)] == seen["digest"][("fgw_shape0x_series", s)]
    assert out.n == ctl.n and not np.array_equal(out.coords, ctl.coords)  # the identity is never tried
    assert ev.prefix_parity(out.coords, board.truth_cells, s) == ev.prefix_parity(floor.coords, board.truth_cells, s)
    assert arm.name == "fgw_shape0x_series_cp" and arm.params["row_order"] == "copy_parity"
    assert arm.requires == base.requires and arm.tasks == base.tasks


def test_the_arm_ships_the_file_as_built_when_the_option_refuses(e48):
    ref, board, base, ctx, s = _setup(e48, n_ref=2000, truth=400)  # the copy (200 rows) is below the truth count
    with e48.captured() as seen:
        seen.update(board=board, ref=ref)
        harness.subsample(ref, ev.SUBSAMPLE_FRAC, s)
        ctl = base.build(ctx)
        out = e48.reordered("fgw_shape0x_series_cp", base, seen).build(ctx)
    assert seen["arm"][("fgw_shape0x_series_cp", s)]["refused"]
    np.testing.assert_array_equal(out.coords, ctl.coords)


def test_the_arm_flags_a_copy_that_is_not_the_floor(e48):
    ref, board, base, ctx, s = _setup(e48)
    with e48.captured() as seen:
        seen.update(board=board, ref=ref)
        harness.subsample(ref, ev.SUBSAMPLE_FRAC, s + 1)  # the floor of another seed only
        e48.reordered("fgw_shape0x_series_cp", base, seen).build(ctx)
    rec = seen["arm"][("fgw_shape0x_series_cp", s)]
    assert not rec["copy_is_floor"] and not rec["scoring_seed_known"]


def test_the_arm_refuses_a_test_board(e48):
    ref, board, base, ctx, s = _setup(e48)
    test = Board(**{**board.__dict__, "key": "T2:embryo:test_interp", "split": "test_interp"})
    with e48.captured() as seen:
        seen.update(board=test, ref=ref)
        with pytest.raises(RuntimeError, match="test board"):
            e48.reordered("fgw_shape0x_series_cp", base, seen).build(ctx)


def test_one_seed_runs_through_the_harness_on_a_synthetic_task(e48, monkeypatch, tmp_path):
    """The whole per-seed path through harness.run_task (both builds), on test_harness's synthetic T2 task with the
    harness's copy candidate standing in for both base candidates: the spy hands the arms the scored floor copy, the
    arms change row order only, and where the option acted and the truth rows are shared the arm's prefix mode against
    the fold truth is the floor's."""
    sys.path.insert(0, str(REPO / "tests"))
    import test_harness as th

    from vec.harness import CANDIDATES

    rng = np.random.default_rng(11)
    base = rng.gamma(2.0, 2.0, size=(3, len(th.GENES)))
    later = base.copy()
    later[:, :6] *= 3.0
    start = th._pool(3000, rng, base, [0.5, 0.3, 0.2])
    ref = th._pool(12000, rng, base * 1.2, [0.4, 0.35, 0.25], scale=1.1)
    target = th._pool(6000, rng, later, [0.3, 0.4, 0.3], scale=1.3, stretch=(1.0, 1.0, 3.0))
    board = Board(key="pseudo:synthetic", task="T2", genes=list(th.GENES), needs_coords=True, setting="heart",
                  split="pseudo", min_cells=100, max_cells=3600, truth_cells=300, ref_cells=1200)
    task = th._task(key="syn")
    th._fake_io(monkeypatch, tmp_path, task, th._stage(target), th._stage(ref), {"ref": ref, "start": start},
                {"scored": [], "inputs": []})
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    monkeypatch.setattr(harness, "load_boards", lambda: {"T2:heart:val_extrap": replace(board, key="T2:heart:val_"
                                                                                                      "extrap")})
    for c in e48.CONTROLS:
        monkeypatch.setitem(harness.CANDIDATES, c, replace(CANDIDATES["copy"], name=c))
    out = e48.one_seed(task, 0, tmp_path)
    assert set(out["builds"]) == {"main", "prefix"} and out["facts"]["files_identical_across_builds"]
    f = out["facts"]
    for a, c in e48.ARMS.items():
        x = f["arms"][a]
        assert x["copy_is_floor"] and x["scoring_seed_known"] and x["metric_seed"] == 0
        assert x["permutation"] and x["pairs_equal_control"] and x["n_true_is_truth_half"]
        assert not x["refused"] and f["truth_rows_shared"][a]
        # at this synthetic seed the file as built falls in the other mode than the floor; the re-ordered file does not
        assert f["hands"][c]["prefix_mode"] != f["hands"]["floor"]["prefix_mode"]
        assert f["hands"][a]["prefix_mode"] == f["hands"]["floor"]["prefix_mode"]
        for b in e48.BUILDS:
            row = out["builds"][b]["rows"][0]
            assert a in row["scores"] and c in row["scores"] and row["seed"] == 0
            assert out["builds"][b]["build"] == b
        # under main the mode is not read; under prefix the arm's occupancy is the floor-mode value
        hp = f["hands"][a]
        prefix_occ = out["builds"]["prefix"]["rows"][0]["metrics"][a]["occupancy_dice"]
        assert prefix_occ == (hp["dice_proper"] if hp["prefix_mode"] == "proper" else hp["dice_mirror"])


# ------------------------------------------------------------------------------------------------ the gate
def test_the_gate_passes_e47s_own_scores_and_fails_a_moved_seed(e48):
    key = e48.nt.all_tasks()["e75"].key
    old, why = e48.e47_scores(key)
    assert why == [] and set(e48.CONTROLS) <= set(old)
    g = e48.gate_check("e75", key, {c: old[c].copy() for c in e48.CONTROLS})
    assert g["pass"] and g["worst"] == {"mean": 0.0, "seed": 0.0}
    moved = {c: old[c].copy() for c in e48.CONTROLS}
    moved["fgw_shape0x"][4] += 0.6
    g = e48.gate_check("e75", key, moved)
    assert not g["pass"] and any("in a seed" in w for w in g["why"])
    assert not e48.gate_check("e75", "T2e/none", moved)["pass"]


# ------------------------------------------------------------------------------------------------ the reading
METRICS = ("occupancy_dice", "d2_shape", "de_score")


def _row(seed: int, arm_shift: dict[str, float]) -> dict:
    base = {"floor": 50.0, "fgw_shape0x_series": 48.0 + 0.1 * seed, "fgw_shape0x": 47.0}
    scores = {n: base.get(n, base[n.removesuffix("_cp")] + arm_shift.get(n, 0.0))
              for n in ("floor", "fgw_shape0x_series", "fgw_shape0x", "fgw_shape0x_series_cp", "fgw_shape0x_cp")}
    skills = {n: {m: 0.5 + (v - 50.0) / 100.0 for m in METRICS} for n, v in scores.items()}
    return {"seed": seed, "scores": scores, "skills": skills}


def _facts(shared=True, modes=None, refused=False) -> dict:
    modes = modes or {}
    hands = {n: {"prefix_mode": modes.get(n, "proper"), "dice_proper": 0.8, "dice_mirror": 0.75, "dice_gap": 0.05,
                 "parity_pred": 1} for n in ("floor", *("fgw_shape0x_series", "fgw_shape0x"),
                                             "fgw_shape0x_series_cp", "fgw_shape0x_cp")}
    arm = {"metric_seed": None, "scoring_seed_known": True, "copy_is_floor": True, "refused": refused,
           "permutation": True, "n_true_is_truth_half": True, "pairs_equal_control": True}
    return {"board_split": "pseudo", "hands": hands,
            "truth_rows_shared": {n: shared for n in ("fgw_shape0x_series", "fgw_shape0x", "fgw_shape0x_series_cp",
                                                      "fgw_shape0x_cp")},
            "arms": {a: dict(arm) for a in ("fgw_shape0x_series_cp", "fgw_shape0x_cp")}}


def _recs(e48, shift=None, refused_seeds=(), modes_on=None) -> dict:
    out = {}
    for name in e48.nt.all_tasks():
        rs = []
        for s in SEEDS:
            ref = s in refused_seeds
            sh = {} if ref else (shift or {})
            f = _facts(modes=(modes_on or {}).get((name, s)), refused=ref)
            for a in f["arms"].values():
                a["metric_seed"] = s
            builds = {b: {"provenance": {"commit": COMMIT, "tree_clean": True}, "build": b, "protocol": "bake",
                          "board": {"split": "pseudo"}, "rows": [_row(s, sh)]} for b in e48.BUILDS}
            rs.append({"seed": s, "run_commit": COMMIT, "builds": builds, "facts": f})
        out[name] = rs
    return out


def test_the_reading_counts_a_refused_seed_as_the_control_in_the_primary_figure(e48):
    recs = _recs(e48, shift={"fgw_shape0x_series_cp": 1.0}, refused_seeds=(0, 1))
    res = e48.reading(recs, [], COMMIT)
    assert not res["void"]
    p = res["tasks"]["e75"]["paired"]["fgw_shape0x_series_cp"]["main"]
    assert p["board"]["mean"] == pytest.approx(0.8) and p["board"]["n"] == 10
    assert p["secondary"]["board"]["mean"] == pytest.approx(1.0) and p["secondary"]["refused_left_out"] == 2
    assert p["occupancy_skill"]["mean"] == pytest.approx(0.008)
    assert set(p["skills"]) == set(METRICS)
    assert res["refused_task_seeds"]["fgw_shape0x_series_cp"] == 18
    r = res["rule_i"]["fgw_shape0x_series_cp"]
    assert r["can_price"] and not r["construction_void"] and r["shared_task_seeds"] == 90
    assert res["tasks"]["e75"]["gap"]["floor"]["gap"]["mean"] == pytest.approx(0.05)
    assert res["parity_agree_as_built"]["fgw_shape0x"]["e75"] == 1.0
    text = e48.text(res)
    assert "equals the floor's on all 90 shared task-seeds" in text and "VOID" not in text


def test_a_mode_differing_on_a_shared_seed_voids_the_construction_not_the_run(e48):
    recs = _recs(e48, modes_on={("x1", 4): {"fgw_shape0x_series_cp": "mirror"}})
    res = e48.reading(recs, [], COMMIT)
    assert not res["void"]
    r = res["rule_i"]["fgw_shape0x_series_cp"]
    assert r["construction_void"] and r["mode_differs_on"] == ["x1 seed 4"]
    assert not res["rule_i"]["fgw_shape0x_cp"]["construction_void"]
    assert "the construction is void" in e48.text(res)


def test_no_shared_seed_means_the_read_cannot_price_the_option(e48):
    recs = _recs(e48)
    for rs in recs.values():
        for r in rs:
            r["facts"]["truth_rows_shared"] = dict.fromkeys(r["facts"]["truth_rows_shared"], False)
    res = e48.reading(recs, [], COMMIT)
    assert not res["rule_i"]["fgw_shape0x_series_cp"]["can_price"]
    assert "cannot price the option" in e48.text(res)


@pytest.mark.parametrize("breach", ["seed", "commit", "copy", "pairs", "n_true", "gate", "test_board", "unclean"])
def test_the_void_conditions(e48, breach):
    recs = _recs(e48)
    gates = []
    r = recs["late"][3]
    if breach == "seed":
        recs["late"].pop(3)
    elif breach == "commit":
        r["builds"]["prefix"]["provenance"]["commit"] = "b" * 40
    elif breach == "copy":
        r["facts"]["arms"]["fgw_shape0x_cp"]["copy_is_floor"] = False
    elif breach == "pairs":
        r["facts"]["arms"]["fgw_shape0x_series_cp"]["pairs_equal_control"] = False
    elif breach == "n_true":
        r["facts"]["arms"]["fgw_shape0x_series_cp"]["n_true_is_truth_half"] = False
    elif breach == "gate":
        gates = ["the gate did not pass on late at the run commit"]
    elif breach == "test_board":
        r["facts"]["board_split"] = "test_interp"
    else:
        r["builds"]["main"]["provenance"]["tree_clean"] = False
    res = e48.reading(recs, gates, COMMIT)
    assert res["void"] and res["why"]
    assert e48.text(res).startswith(f"Try E48 (run commit {COMMIT}): VOID")


def test_gate_files_are_required_at_the_run_commit(e48, tmp_path):
    assert len(e48.gate_problems(tmp_path, COMMIT)) == 9
    for name in e48.nt.all_tasks():
        p = e48._gate_path(tmp_path, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"pass": True, "commit": COMMIT}))
    assert e48.gate_problems(tmp_path, COMMIT) == []
    assert len(e48.gate_problems(tmp_path, "c" * 40)) == 9


def test_the_read_refuses_to_overwrite(e48, tmp_path):
    (tmp_path / "reading.txt").write_text("x")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        e48.main(["read", "--commit", COMMIT, "--runs", str(tmp_path), "--out", str(tmp_path)])


# ------------------------------------------------------------------------------------------------ the launcher
def test_the_launcher_caps_threads_runs_nine_processes_and_writes_exit_files():
    sh = (REPO / "experiments" / "coords" / "run_e48.sh").read_text()
    for need in ("OMP_NUM_THREADS=8", "OPENBLAS_NUM_THREADS=8", "MKL_NUM_THREADS=8", "NUMPY_MADVISE_HUGEPAGE=0",
                 'finish 98', 'finish 99', '"$S/$STEP.exit"', "uv run --frozen --no-sync python -I -B",
                 "tests/test_e48_row_order_read.py", "worktree"):
        assert need in sh
    tasks = sh.split("TASKS=(")[1].split(")")[0].split()
    assert sorted(tasks) == sorted(["half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"])
    assert len(tasks) <= 10 and "gate)" in sh and "run)" in sh and "read)" in sh


def test_the_option_record_and_constants_are_the_builders(e48):
    assert tb.ROW_ORDERS == ("as_built", "copy_parity") and tb.COPY_PARITY_STREAM == 48
