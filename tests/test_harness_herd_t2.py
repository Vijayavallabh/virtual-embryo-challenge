"""Tries H24 and E46 (experiments/harness/TRIES.md; decision log 2026-11-18 (c)): herded draws on each T2 board's
first-upload entry. The heart's new candidate `b0_trend_herd` and the embryo's twin `fgw_herd` against their entries,
and the reading `experiments/harness/n53_herd_t2_readings.py` on synthetic score tables. Synthetic data only; nothing
is read from data/ or from a committed run, and nothing is scored."""

from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from test_b0_trend import BOARD, _ctx, three  # noqa: F401
from test_harness import GENES, INTERP, _pool, _task

from vec import harness
from vec.boards import Board
from vec.harness import CANDIDATES, load_tasks

HARNESS = Path(__file__).resolve().parents[1] / "experiments" / "harness"


def _bound(fn, args, kwargs) -> dict:
    b = inspect.signature(fn).bind(*args, **kwargs)
    b.apply_defaults()
    return dict(b.arguments)


def _rejecting(monkeypatch, log: list):
    """Replace t2_build's herd_rows by the real one aimed at the drawn rows' own mean: every proposal is drawn from its
    stream as before and every one is rejected (no swap brings the mean closer to where it already is)."""
    real = harness.t2_build.herd_rows

    def rejecting(X, pool, picked, rng, proposals=None, target=None):
        Xs = sparse.csr_matrix(X, dtype=np.float64)
        own = np.asarray(Xs[np.asarray(picked)].sum(0)).ravel() / len(picked)
        rows, rec = real(X, pool, picked, rng, proposals=proposals, target=own)
        log.append(rec)
        return rows, rec

    monkeypatch.setattr(harness.t2_build, "herd_rows", rejecting)


# ---------------------------------------------------------------------------------------------------- the heart
def test_b0_trend_herd_is_registered_beside_b0_trend():
    tasks = {t.key: t for ts in load_tasks().values() for t in ts}
    c, base = CANDIDATES["b0_trend_herd"], CANDIDATES["b0_trend"]
    assert {k for k, t in tasks.items() if c.applies(t)} == {k for k, t in tasks.items() if base.applies(t)} == {
        "T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"}
    assert c.tasks == base.tasks and c.requires is base.requires and c.tuned_on == base.tuned_on == frozenset()
    assert {k: v for k, v in c.params.items() if k != "what"} == {k: v for k, v in base.params.items() if k != "what"}


def test_b0_trend_herd_takes_b0_trends_route_with_herd_set(three, monkeypatch):  # noqa: F811
    """The same steps with the same arguments, herd set; and the copy's draw is the herded draw_copy (b0_herd's)."""
    calls, draws = [], []
    steps, draw = harness._b0_t2_steps, harness.draw_copy

    def spy_steps(*a, **kw):
        calls.append(_bound(steps, a, kw))
        return steps(*a, **kw)

    def spy_draw(ctx, herd=False):
        draws.append(herd)
        return draw(ctx, herd=herd)

    monkeypatch.setattr(harness, "_b0_t2_steps", spy_steps)
    monkeypatch.setattr(harness, "draw_copy", spy_draw)
    ctx = _ctx(three)
    plain, herded = CANDIDATES["b0_trend"].build(ctx), CANDIDATES["b0_trend_herd"].build(ctx)
    assert len(calls) == 2 and draws == [False, True]
    assert calls[0]["ctx"] is calls[1]["ctx"] is ctx
    assert {**calls[0], "herd": True} == calls[1]
    assert calls[0] == {"ctx": ctx, "depth": False, "herd": False, "cfm": False, "trend": True, "trend_record": None}
    assert plain.n == herded.n and not np.array_equal(plain.X, herded.X)  # herding moves the draw here
    # the route that b0_trend and the build use still refuses the trend with herd
    with pytest.raises(ValueError, match="depth, herd or cfm"):
        harness._b0_t2(ctx, trend=True, herd=True)


def test_b0_trend_herd_draw_moves_towards_the_pool_mean_and_keeps_positions(three):  # noqa: F811
    ctx = _ctx(three)
    pool, rows, X0, _ = harness.draw_copy(ctx)
    _, rows_h, Xh, coords_h = harness.draw_copy(ctx, herd=True)
    assert len(rows_h) == len(rows) and len(np.unique(rows_h)) == len(rows_h)
    np.testing.assert_allclose(coords_h, np.asarray(pool.obsm["spatial_3D"], float)[rows_h], atol=1e-5)
    goal = np.asarray(pool[:, BOARD.genes].X.mean(0)).ravel()
    assert ((Xh.mean(0) - goal) ** 2).sum() < ((X0.mean(0) - goal) ** 2).sum()


def test_b0_trend_herd_with_every_proposal_rejected_is_b0_trend_bit_for_bit(three, monkeypatch):  # noqa: F811
    ctx = _ctx(three)
    plain = CANDIDATES["b0_trend"].build(ctx)
    log: list = []
    _rejecting(monkeypatch, log)
    herded = CANDIDATES["b0_trend_herd"].build(ctx)
    assert len(log) == 1 and log[0]["proposals"] > 0 and log[0]["accepted"] == 0
    np.testing.assert_array_equal(herded.X, plain.X)
    np.testing.assert_array_equal(herded.coords, plain.coords)


# ---------------------------------------------------------------------------------------------------- the embryo
@pytest.fixture(scope="module")
def embryo():
    rng = np.random.default_rng(31)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    before = _pool(900, rng, base, [0.5, 0.3, 0.2])
    after = _pool(1200, rng, base * 1.5, [0.3, 0.4, 0.3], scale=1.8, stretch=(1.0, 1.4, 1.0))
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=600)
    task = _task(**INTERP, key="syn_embryo_herd", inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad"})
    return harness.BuildContext(task=harness.TaskView.of(task), board=board, inputs={"before": before, "after": after},
                                seed=harness.BUILD_SEED_OFFSET + 2)


def test_the_embryos_first_upload_is_fgw_and_its_twin_is_fgw_herd():
    """build_all.sh's first round builds the embryo with `--method fgw --expression mixture` and the move coupling and
    draw at their defaults: the harness's fgw (_fgw's defaults), whose herded twin is the registered fgw_herd."""
    root = HARNESS.parents[1]
    script = (root / "experiments" / "builds" / "build_all.sh").read_text()
    assert "EMBRYO_EXPRESSION=mixture" in script
    s1 = script.split("s1)")[1].split(";;")[0]
    assert 'embryo --method fgw --expression "$EMBRYO_EXPRESSION"' in s1 and "--move-coupling" not in s1
    sig = inspect.signature(harness._fgw).parameters
    assert sig["expression"].default == "mixture" and sig["move_coupling"].default == harness.t2_build.MOVE_COUPLINGS[0]
    assert sig["draw"].default == harness.t2_build.DRAWS[0] and sig["t_pos"].default is None
    assert "fgw_trend_herd" not in CANDIDATES and "fgw_herd" in CANDIDATES


def test_fgw_herd_takes_fgws_route_with_herd_set(embryo, monkeypatch):
    calls = []
    real = harness.t2_build.interpolate_clouds

    def spy(*a, **kw):
        calls.append(_bound(real, a, kw))
        return real(*a, **kw)

    monkeypatch.setattr(harness.t2_build, "interpolate_clouds", spy)
    plain, herded = CANDIDATES["fgw"].build(embryo), CANDIDATES["fgw_herd"].build(embryo)
    assert len(calls) == 2
    a, b = calls
    assert set(a) == set(b) and a["herd"] is False and b["herd"] is True
    for k in a:
        if k != "herd":
            v, w = a[k], b[k]
            assert (np.array_equal(v, w) if isinstance(v, np.ndarray) else v == w), k
    assert plain.n == herded.n and not np.array_equal(plain.X, herded.X)


def test_fgw_herd_with_every_proposal_rejected_is_fgw_bit_for_bit(embryo, monkeypatch):
    plain = CANDIDATES["fgw"].build(embryo)
    log: list = []
    _rejecting(monkeypatch, log)
    herded = CANDIDATES["fgw_herd"].build(embryo)
    assert len(log) == 2 and all(r["proposals"] > 0 and r["accepted"] == 0 for r in log)
    np.testing.assert_array_equal(herded.X, plain.X)
    np.testing.assert_array_equal(herded.coords, plain.coords)


# ---------------------------------------------------------------------------------------------------- the reading
@pytest.fixture(scope="module")
def n53():
    sys.path.insert(0, str(HARNESS))
    spec = importlib.util.spec_from_file_location("n53_herd_t2_readings", HARNESS / "n53_herd_t2_readings.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # its dataclass resolves its annotations through the module
    spec.loader.exec_module(mod)
    return mod


COMMIT = "a" * 40


def _file(key: str, board: str, scores: dict[str, np.ndarray], commit: str = COMMIT, clean: bool = True,
          shape: float = 0.0) -> dict:
    """A harness task file as run.py writes it, with the fields the reading and check_results look at."""
    arms = list(scores)
    rows = []
    for s in range(10):
        metrics = {a: {"d2_shape": 0.1 + shape * (a.endswith("herd")), "occupancy_dice": 0.5, "scale_log_ratio": 0.02}
                   for a in arms}
        skills = {a: {"d2_shape": 0.4, "occupancy_dice": 0.5, "scale_log_ratio": 0.6} for a in arms}
        rows.append({"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
                     "scores": {a: float(scores[a][s]) for a in arms}, "skills": skills, "metrics": metrics})
    return {"task": {"key": key, "test_board": board}, "candidates": arms, "rows": rows, "protocol": harness.PROTOCOL,
            "build": harness.BUILD, "anchors": "computed", "frac": harness.SUBSAMPLE_FRAC, "seeds": list(range(10)),
            "provenance": {"commit": commit, "tree_clean": clean, "config_sha256": harness._sha256(harness.CONFIG),
                           "files_sha256": {}, "derived_sha256": {}}}


def _write(folder: Path, files: list[dict]) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    for f in files:
        name = "".join(ch if ch.isalnum() or ch == "." else "_" for ch in f["task"]["key"])
        (folder / f"{name}.json").write_text(json.dumps(f))
    return folder


def _heart(n53, gain: float | list[float] = 0.0, pair_gain: float = 0.0, **kw) -> list[dict]:
    spec = n53.BOARDS["heart"]
    ramp = np.linspace(0.0, 0.9, 10)
    out = []
    for i, key in enumerate(("T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5")):
        g = gain[i] if isinstance(gain, list) else gain
        entry = 60.0 + ramp
        out.append(_file(key, spec.test_board, {"floor": np.full(10, 50.0), "b0": 58.0 + ramp,
                                                "b0_herd": 58.0 + ramp + pair_gain, "b0_trend": entry,
                                                "b0_trend_herd": entry + g}, **kw))
    return out


def test_the_reading_counts_the_boards_backtests_and_names_the_arms(n53):
    heart, embryo = n53.BOARDS["heart"], n53.BOARDS["embryo"]
    assert heart.arms == ("b0_trend", "b0_trend_herd", "b0", "b0_herd")
    assert embryo.arms == ("fgw", "fgw_herd")  # the registered pair is the entry and its twin
    tasks = load_tasks()
    assert n53.counted(heart, tasks[heart.test_board]) == ["T2h/<=E9.5->E10.5", "T2h/<=E8.75->E10.5"]
    assert n53.counted(embryo, tasks[embryo.test_board]) == ["T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25"]
    assert n53.POSITION_ONLY == ("d2_shape", "occupancy_dice", "scale_log_ratio")


def test_a_herded_entry_ahead_on_both_backtests_passes(n53, tmp_path):
    gain = np.r_[np.full(8, 0.6), -0.1, 0.2]  # mean 0.49, 9 of 10 ahead, largest single loss 0.1
    files = _heart(n53)
    for f in files:
        for row, g in zip(f["rows"], gain, strict=True):
            row["scores"]["b0_trend_herd"] += g
    res = n53.reading("heart", _write(tmp_path / "run", files), committed_dir=tmp_path / "none")
    assert not res["void"] and res["passes"] is True and res["statement"]["keep"] is True
    t = res["statement"]["tasks"]["T2h/<=E9.5->E10.5"]
    assert t["mean_gain"] == pytest.approx(0.49) and t["seeds_ahead"] == 9 and t["largest_loss"] == pytest.approx(0.1)
    assert res["registered_pair"]["keep"] is False and not res["pair_is_the_statement"]
    assert all(not r["checked"] for r in res["reproduction"].values())
    txt = n53.text(res)
    assert "the herded entry passes the adoption rule" in txt and "could not be checked" in txt
    assert "decides nothing" in txt or "deciding nothing" in txt
    assert "no biological content" in txt and "working copy" in txt


@pytest.mark.parametrize("gain, why", [(0.2, "passes rule 7's test"), ([0.5, -0.5], "loses more than")])
def test_a_small_or_mixed_gain_does_not_pass(n53, tmp_path, gain, why):
    files = _heart(n53, gain=gain)  # every seed ahead (or behind) by the same gain
    res = n53.reading("heart", _write(tmp_path / "run", files), committed_dir=tmp_path / "none")
    assert not res["void"] and res["passes"] is False and why in res["statement"]["reason"]
    assert "the herded entry does not pass the adoption rule" in n53.text(res)


def test_the_entry_must_reproduce_its_committed_scores(n53, tmp_path):
    files = _heart(n53)
    committed = _heart(n53, commit="b" * 40)
    ok = n53.reading("heart", _write(tmp_path / "run", files), committed_dir=_write(tmp_path / "c1", committed))
    assert not ok["void"] and all(r["checked"] and r["max_abs"] == 0.0 for r in ok["reproduction"].values())
    assert "reproduces its committed scores" in n53.text(ok)
    committed[0]["rows"][3]["scores"]["b0_trend"] += 1e-6
    bad = n53.reading("heart", tmp_path / "run", committed_dir=_write(tmp_path / "c2", committed))
    assert bad["void"] and any("differs from its committed scores" in w for w in bad["why"])
    committed[0]["rows"][3]["scores"]["b0_trend"] -= 1e-6 - 5e-10  # within 1e-9
    near = n53.reading("heart", tmp_path / "run", committed_dir=_write(tmp_path / "c3", committed))
    assert not near["void"]
    one = n53.reading("heart", tmp_path / "run", committed_dir=_write(tmp_path / "c4", committed[1:]))
    assert not one["void"] and [r["checked"] for r in one["reproduction"].values()] == [False, True]


def test_void_cases(n53, tmp_path):
    def void(files, extra=None):
        d = tmp_path / f"run{len(list(tmp_path.iterdir()))}"
        _write(d, files)
        if extra:
            extra(d)
        res = n53.reading("heart", d, committed_dir=tmp_path / "none")
        assert res["void"], res
        assert n53.text(res).startswith("heart (T2:heart:test_extrap): VOID")
        return " ".join(res["why"])

    assert "no task file" in void([])
    two = _heart(n53)
    assert "no task file in the run" in void(two[:1])
    assert "commits" in void([two[0], _heart(n53, commit="c" * 40)[1]])
    assert "dirty tree" in void(_heart(n53, clean=False))
    missing = _heart(n53)
    for row in missing[1]["rows"]:
        del row["scores"]["b0_herd"]
    assert "b0_herd is not scored" in void(missing)
    nine = _heart(n53)
    nine[0]["rows"] = nine[0]["rows"][:9]
    assert "b0_trend is not scored on each of the ten seeds" in void(nine)
    nan = _heart(n53)
    nan[1]["rows"][2]["scores"]["b0_trend_herd"] = float("nan")
    assert "b0_trend_herd is not scored" in void(nan)
    assert "a backtest failed" in void(_heart(n53), lambda d: (d / "x.failed.txt").write_text("boom"))
    seeds = _heart(n53)
    seeds[0]["rows"][0]["build_seed"] += 1
    assert "build seeds" in void(seeds)


def test_the_embryo_reports_its_pair_once_and_position_metrics_beside_the_entry(n53, tmp_path):
    spec = n53.BOARDS["embryo"]
    ramp = np.linspace(0.0, 0.9, 10)
    files = [_file(k, spec.test_board, {"floor": np.full(10, 50.0), "fgw": 55.0 + ramp, "fgw_herd": 55.0 + ramp - 0.2},
                   shape=0.05) for k in ("T2e/E7.25+E8.0->E7.5", "T2e/E6.75+E7.5->E7.25")]
    res = n53.reading("embryo", _write(tmp_path / "run", files), committed_dir=tmp_path / "none")
    assert not res["void"] and res["pair_is_the_statement"] and res["passes"] is False
    assert res["registered_pair"]["tasks"] == res["statement"]["tasks"]
    pos = res["position_only"]["T2e/E7.25+E8.0->E7.5"]
    assert pos["fgw_herd"]["d2_shape"] - pos["fgw"]["d2_shape"] == pytest.approx(0.05)
    assert pos["fgw"]["shape_growth"] == pytest.approx(0.5)
    txt = n53.text(res)
    assert "is the entry and its herded twin themselves" in txt and "(minus the entry's: +0.0500" in txt


def test_the_cli_writes_both_files_and_refuses_to_overwrite(n53, tmp_path, monkeypatch, capsys):
    run = _write(tmp_path / "run", _heart(n53))
    monkeypatch.setitem(n53.BOARDS, "heart", n53.BOARDS["heart"].__class__(
        **{**n53.BOARDS["heart"].__dict__, "committed": tmp_path / "none"}))
    out = tmp_path / "out"
    assert n53.main(["read", "--board", "heart", "--runs", str(run), "--out", str(out)]) == 0
    res = json.loads((out / "n53_heart.json").read_text())
    assert res["void"] is False and "provenance" in res and (out / "n53_heart.txt").read_text() == n53.text(
        {k: v for k, v in res.items() if k != "provenance"})
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        n53.main(["read", "--board", "heart", "--runs", str(run), "--out", str(out)])
