"""vec.t2_build on synthetic data only (no competition data, no board anchors).

Size JSONs are produced by experiments/coords/size_procedure.py itself on synthetic stage geometries with random
growth, so the reader is tested against the procedure's real output schema; none of the numbers means anything.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import re
import sys
import tokenize
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from conftest import GENES, simulate_stage

from vec import io as vec_io
from vec import t2_build as tb
from vec.boards import Board
from vec.evaluate import SCORER_BUILDS, Stage, score_board
from vec.io import COORD_KEY, check_submission, read_stage
from vec.provenance import read_stamp, stamp_path

REPO = Path(__file__).resolve().parents[1]


def t2_board(min_cells=10, max_cells=3000, truth_cells=None, key="T2:synthetic", task="T2", needs_coords=True):
    return Board(key=key, task=task, genes=GENES, needs_coords=needs_coords, min_cells=min_cells,
                 max_cells=max_cells, truth_cells=truth_cells)


def stack(n, seed=0, extent=(120.0, 80.0, 60.0), step=15.0):
    """A uniformly filled box cut into z sections, rows stored section by section in descending z (like the released
    MERFISH files)."""
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 1.0, size=(n, 3)) * np.asarray(extent)
    C[:, 2] = np.round(C[:, 2] / step) * step
    return C[np.argsort(-C[:, 2], kind="stable")]


def random_rotation(seed):
    Q, R = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


@pytest.fixture(scope="module")
def sp():
    """experiments/coords/size_procedure.py, loaded by path."""
    spec = importlib.util.spec_from_file_location("coords_size_procedure",
                                                  REPO / "experiments" / "coords" / "size_procedure.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------------------------------- cell count
def test_count_is_max_cells_with_a_record():
    n, rec = tb.n_cells(t2_board(min_cells=100, max_cells=3000, truth_cells=800), pool_size=5000)
    assert n == 3000
    assert rec["rule"] == "max_cells" and rec["floor"] == 800 and rec["floor_source"] == "truth_cells"
    assert rec["exceeds_pool"] is False and rec["repeated_rows"] == 0 and rec["fraction_of_max_cells"] == 1.0


def test_count_never_below_truth_cells_and_overrides_are_refused_not_clipped():
    board = t2_board(min_cells=100, max_cells=3000, truth_cells=800)
    assert tb.n_cells(board, 5000, override="0.5max")[0] == 1500
    assert tb.n_cells(board, 5000, override="1max")[0] == 3000
    assert tb.n_cells(board, 5000, override=800)[0] == 800  # exactly truth_cells is allowed
    for bad in (799, "0.25max", 3001):  # below truth_cells (750 cells for 0.25max) or above max_cells
        with pytest.raises(ValueError, match="refused"):
            tb.n_cells(board, 5000, override=bad)
    for bad in ("1.5max", "half", "0max", "2x"):
        with pytest.raises(ValueError):
            tb.n_cells(board, 5000, override=bad)
    for bad in (0.5, True):
        with pytest.raises(TypeError):
            tb.n_cells(board, 5000, override=bad)


def test_count_floor_is_min_cells_without_truth_cells_and_bad_boards_are_refused():
    n, rec = tb.n_cells(t2_board(min_cells=100, max_cells=3000), pool_size=5000, override="0.25max")
    assert n == 750 and rec["floor"] == 100 and rec["floor_source"] == "min_cells" and rec["truth_cells"] is None
    with pytest.raises(ValueError, match="no max_cells"):
        tb.n_cells(t2_board(max_cells=None), 5000)
    with pytest.raises(ValueError, match="inconsistent"):
        tb.n_cells(t2_board(min_cells=10, max_cells=500, truth_cells=800), 5000)
    with pytest.raises(ValueError):
        tb.n_cells(t2_board(), 0)


def test_count_above_the_pool_is_annotated():
    n, rec = tb.n_cells(t2_board(min_cells=100, max_cells=3000), pool_size=1200)
    assert n == 3000 and rec["exceeds_pool"] and rec["repeated_rows"] == 1800
    assert rec["x_pool"] == pytest.approx(2.5) and "2 or 3 times" in rec["note"]


# ---------------------------------------------------------------------------------------------------- count rule
def _rerun(gains=None, seeds=range(10), base=60.0, seed=0):
    """Density-rerun scores {level: {seed: score}}: "1max" around `base`, and each level in `gains` the 1max score of
    the same seed plus its gain there (one number for every seed, or one per seed). None of the numbers means
    anything."""
    rng = np.random.default_rng(seed)
    one = {s: base + float(rng.normal()) for s in seeds}
    out = {"1max": one}
    for level, g in (gains or {}).items():
        g = np.broadcast_to(np.asarray(g, dtype=float), (len(one),))
        out[level] = {s: one[s] + float(d) for s, d in zip(one, g, strict=True)}
    return out


LOSERS = {"0.25max": -2.0, "0.5max": -1.0}  # levels behind max_cells in every seed


def test_count_rule_keeps_max_cells_when_no_level_passes():
    for gains in ({**LOSERS, "0.75max": 0.2},  # ahead in all 10 seeds, mean gain 0.2
                  {**LOSERS, "0.75max": [1.0] * 7 + [-0.5] * 3},  # mean gain 0.55, 7 seeds ahead
                  {**LOSERS, "0.75max": [0.3] * 8 + [-0.2] * 2}):  # 8 seeds ahead, mean gain 0.2
        override, rec = tb.count_override(_rerun(gains))
        assert override is None and rec["override"] is None and rec["complete"] and rec["passing"] == []
        assert rec["decision"].startswith("keep max_cells: no lower level")
        assert [rec["per_level"][lv]["passes"] for lv in tb.COUNT_RULE_LEVELS] == [False] * 3
    assert rec["per_level"]["0.75max"]["seeds_ahead"] == 8
    assert rec["per_level"]["0.75max"]["mean_gain"] == pytest.approx(0.2)


def test_count_rule_needs_8_seeds_ahead_besides_the_mean_gain():
    override, rec = tb.count_override(_rerun({**LOSERS, "0.75max": [1.0] * 7 + [-0.5] * 3}))
    lv = rec["per_level"]["0.75max"]
    assert override is None and lv["mean_gain"] == pytest.approx(0.55) and lv["seeds_ahead"] == 7
    assert lv["passes"] is False
    override, rec = tb.count_override(_rerun({**LOSERS, "0.75max": [0.6] * 8 + [-0.2] * 2}))  # 8 ahead, mean 0.44
    assert override == "0.75max" and rec["per_level"]["0.75max"]["seeds_ahead"] == 8
    override, rec = tb.count_override(_rerun({**LOSERS, "0.75max": [1.0] * 7 + [0.0] * 3}))  # a tie is not ahead
    assert override is None and rec["per_level"]["0.75max"]["seeds_ahead"] == 7


def test_count_rule_one_level_passes_and_n_cells_takes_it():
    override, rec = tb.count_override(_rerun({**LOSERS, "0.75max": 0.5}))
    assert override == "0.75max" == rec["override"] and rec["passing"] == ["0.75max"] and rec["complete"]
    assert rec["decision"].startswith("use 0.75max: paired mean gain +0.500 points over max_cells, 10 of the 10 seeds")
    lv = rec["per_level"]["0.75max"]
    assert lv["paired_seeds"] == list(range(10)) and lv["differences"] == pytest.approx([0.5] * 10)
    assert lv["mean_gain"] == pytest.approx(0.5) and lv["seeds_ahead"] == 10 and lv["passes"] is True
    assert rec["per_level"]["0.25max"]["mean_gain"] == pytest.approx(-2.0)
    assert rec["approved"] == "2026-10-21" and rec["rule"].startswith("approved rule 7")
    assert rec["condition"] == "shuffle50" and rec["seeds"] == list(range(10)) and rec["base_seeds"] == list(range(10))
    n, count = tb.n_cells(t2_board(min_cells=100, max_cells=3000, truth_cells=800), 5000, override=override)
    assert n == 2250 and count["rule"] == "override '0.75max'"


def test_count_rule_largest_gain_wins_and_a_tie_goes_to_the_higher_count():
    override, rec = tb.count_override(_rerun({"0.25max": -2.0, "0.5max": 0.8, "0.75max": 0.5}))
    assert override == "0.5max" and rec["passing"] == ["0.5max", "0.75max"]
    assert "the largest mean gain of the passing levels ['0.5max', '0.75max']" in rec["decision"]
    override, rec = tb.count_override(_rerun({"0.25max": 0.6, "0.5max": -1.0, "0.75max": 0.6}))
    assert rec["per_level"]["0.25max"]["mean_gain"] == rec["per_level"]["0.75max"]["mean_gain"]  # an exact tie
    assert override == "0.75max"
    assert tb.count_override(_rerun({"0.25max": 0.9, "0.5max": 0.9, "0.75max": 0.4}))[0] == "0.5max"


def test_count_rule_keeps_max_cells_when_the_rerun_is_incomplete():
    gains = {**LOSERS, "0.75max": 0.5}  # with every seed, 0.75max passes
    override, rec = tb.count_override(_rerun(gains, seeds=range(9)))  # 9 seeds
    assert override is None and rec["complete"] is False
    assert rec["incomplete"] == [f"{lv} lacks seeds [9]" for lv in ("1max", *tb.COUNT_RULE_LEVELS)]
    assert "the rerun is incomplete" in rec["decision"] and "keeps max_cells" in rec["decision"]
    assert all(lv["passes"] is None for lv in rec["per_level"].values()) and rec["passing"] == []
    scores = _rerun(gains)
    del scores["0.25max"][4]  # one seed of a level that loses anyway: the rule reads every level's 10 seeds
    override, rec = tb.count_override(scores)
    assert override is None and rec["incomplete"] == ["0.25max lacks seeds [4]"]
    scores = _rerun(gains)
    del scores["1max"]
    override, rec = tb.count_override(scores)
    assert override is None and rec["incomplete"] == ["no 1max scores"] and rec["base_seeds"] == []
    assert rec["per_level"]["0.75max"]["paired_seeds"] == [] and rec["per_level"]["0.75max"]["mean_gain"] is None
    scores = _rerun(gains)
    del scores["0.5max"]
    assert tb.count_override(scores)[1]["incomplete"] == ["no 0.5max scores"]
    override, rec = tb.count_override({})  # nothing ran
    assert override is None and len(rec["incomplete"]) == 4


def test_count_rule_treats_nan_scores_as_missing_seeds():
    gains = {**LOSERS, "0.75max": 0.5}
    for missing in (float("nan"), np.nan, None, float("inf"), np.float32("nan")):
        scores = _rerun(gains)
        scores["0.75max"][2] = missing
        override, rec = tb.count_override(scores)
        assert override is None and rec["incomplete"] == ["0.75max lacks seeds [2]"]
        assert rec["per_level"]["0.75max"]["paired_seeds"] == [0, 1, 3, 4, 5, 6, 7, 8, 9]
    scores = _rerun(gains)
    scores["1max"][7] = float("nan")
    override, rec = tb.count_override(scores)
    assert override is None and rec["incomplete"] == ["1max lacks seeds [7]"]


def test_count_rule_ignores_other_seeds_and_levels_and_refuses_malformed_input():
    scores = _rerun({**LOSERS, "0.75max": -0.5})
    for level in ("1max", *tb.COUNT_RULE_LEVELS):
        scores[level][10] = 0.0 if level == "1max" else 99.0  # a seed beyond the release seeds is never counted
    scores["0.9max"] = {s: 99.0 for s in range(10)}  # a level the rerun does not try
    override, rec = tb.count_override(scores)
    assert override is None and rec["complete"] and rec["ignored_seeds"] == [10]
    assert rec["ignored_levels"] == ["0.9max"]
    typed = {lv: {np.int64(s): np.float64(v) for s, v in d.items()} for lv, d in
             _rerun({**LOSERS, "0.75max": 0.5}).items()}
    assert tb.count_override(typed)[0] == "0.75max"  # numpy seeds and scores
    for bad_level in ("half", "75%", 0.75):
        with pytest.raises(ValueError, match="<f>max"):
            tb.count_override({bad_level: {0: 60.0}})
    for bad_score in ("60", True):
        with pytest.raises(TypeError, match="score"):
            tb.count_override({"1max": {0: bad_score}})
    for bad_seed in (True, "0", 1.0):
        with pytest.raises(TypeError, match="seeds"):
            tb.count_override({"1max": {bad_seed: 60.0}})
    with pytest.raises(TypeError):
        tb.count_override([("1max", {0: 60.0})])
    with pytest.raises(TypeError):
        tb.count_override({"1max": [60.0]})


def _load_by_path(name, *parts):
    spec = importlib.util.spec_from_file_location(name, REPO.joinpath(*parts))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ds():
    """experiments/density/summarise.py (its decide mode is the count rule's release-day entry point)."""
    return _load_by_path("density_summarise", "experiments", "density", "summarise.py")


def test_count_rule_names_are_the_density_reruns(ds):
    rd = _load_by_path("density_run_density", "experiments", "density", "run_density.py")
    assert [f"{f:g}max" for f in rd.MAX_FRACS] == [*tb.COUNT_RULE_LEVELS, tb.COUNT_RULE_BASE]
    labels = {lab for lab, *_ in rd.conditions(np.arange(4000), np.arange(3000), 100, 1000, 0)}
    assert {f"{tb.COUNT_RULE_CONDITION}_{lv}" for lv in (*tb.COUNT_RULE_LEVELS, tb.COUNT_RULE_BASE)} <= labels
    assert rd.SHUFFLE_FRAC == 0.5  # the intermediate-quality condition: half the expression rows shuffled
    assert tb.COUNT_RULE_TARGETS == {"heart": "heart_E10.5", "embryo": "embryo_E7.5"}
    for kind, target in tb.COUNT_RULE_TARGETS.items():  # the validation answers, each on its own kind's board
        assert rd.TARGETS[target].get("release") and f":{kind}:" in rd.TARGETS[target]["board"]
    assert all((REPO / p).is_file() for p in ds.RULE_FILES)
    # the release decision's record: the --out of the release command in summarise.py's docstring, in this repo
    assert ds.release_record() == REPO / "experiments" / "density" / "results" / "release" / "count_decision.json"
    assert f"--out {ds.RELEASE_OUT} --decide" in " ".join(ds.__doc__.split())


def _raw_row(target, seed, condition, score, basis="max", n_pred=None, max_cells=4000, pool_cells=36000):
    """A raw row as experiments/density/run_density.py writes it, with the fields the decide mode reads."""
    pattern = condition.split("_")[0]
    return {"target": target, "board": "T2:synthetic", "seed": seed, "condition": condition, "source": "target",
            "pattern": pattern, "basis": basis, "n_pred": n_pred, "max_cells": max_cells,
            "max_cells_source": "board", "pool_cells": pool_cells, "score": score}


CLEAN_CODE = {"git_commit": "synthetic", "git_unclean": [], "rule_files_uncommitted": [],
              "rule_files_sha256": dict.fromkeys(("src/vec/t2_build.py", "experiments/density/summarise.py"),
                                                 "synthetic")}


@pytest.fixture
def release_json(ds, tmp_path, monkeypatch):
    """The release decision's canonical record (experiments/density/results/release/count_decision.json in the repo),
    moved under tmp_path and absent until a test writes it: the tests do not depend on whether this tree holds the
    release decision."""
    path = tmp_path / "canonical" / "results" / "release" / ds.DECISION_FILE
    monkeypatch.setattr(ds, "release_record", lambda: path)
    return path


@pytest.fixture
def clean_code(ds, monkeypatch, release_json):
    """The decide mode's view of the rule code, as committed: the tests do not depend on this tree's git state (nor,
    through release_json, on its release decision)."""
    monkeypatch.setattr(ds, "rule_code_state", lambda repo: dict(CLEAN_CODE))


def _write_raw(raw, target, scores, seeds=range(10)):
    """One <target>_s<seed>.jsonl per seed with the shuffle50 rows at the max_cells levels of `scores` and, around
    them, rows the rule must not read (other conditions, a truth-multiple level)."""
    n_of = {"0.25max": 1000, "0.5max": 2000, "0.75max": 3000, "1max": 4000}
    for seed in seeds:
        rows = []
        for level, per_seed in scores.items():
            if seed in per_seed:
                rows.append(_raw_row(target, seed, f"shuffle50_{level}", per_seed[seed], n_pred=n_of[level]))
                rows.append(_raw_row(target, seed, f"random_{level}", per_seed[seed] + 9.0, n_pred=n_of[level]))
        rows.append(_raw_row(target, seed, "shuffle50_1x", 99.0, basis="truth", n_pred=500))
        (raw / f"{target}_s{seed}.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))


def test_decide_mode_writes_the_count_decision_of_each_test_board(ds, tmp_path, capsys, clean_code):
    raw, out = tmp_path / "raw_release", tmp_path / "release"
    raw.mkdir()
    _write_raw(raw, "heart_E10.5", _rerun({**LOSERS, "0.75max": 0.5}, seed=1))  # 0.75max passes
    embryo = _rerun({"0.25max": 1.0, "0.5max": 1.0, "0.75max": 1.0}, seed=2)
    embryo["0.5max"][3] = float("nan")  # run_density writes NaN for a score it could not compute
    _write_raw(raw, "embryo_E7.5", embryo, seeds=range(9))  # seed 9 never ran
    _write_raw(raw, "heart_E9.5", _rerun({"0.25max": 5.0, "0.5max": 5.0, "0.75max": 5.0}))  # a pseudo-target
    first = raw / "heart_E10.5_s0.jsonl"  # a stale row, replaced by the later one for the same seed and condition
    stale = _raw_row("heart_E10.5", 0, "shuffle50_0.75max", -50.0, n_pred=3000)
    first.write_text(json.dumps(stale) + "\n" + first.read_text())

    res = ds.decide(raw, out)  # out does not exist yet: decide makes it
    text = (out / "count_decision.json").read_text()
    assert json.loads(text) == res and "NaN" not in text
    heart, embryo_b = res["boards"]["heart"], res["boards"]["embryo"]
    assert heart["validation_target"] == "heart_E10.5" and heart["override"] == "0.75max"
    assert heart["record"]["per_level"]["0.75max"]["differences"] == pytest.approx([0.5] * 10)
    assert heart["decision"] == heart["record"]["decision"] and heart["row_problems"] == []
    assert heart["max_cells"] == [4000] and heart["max_cells_source"] == ["board"] and heart["rows"] == 40
    assert heart["n_pred"] == {"0.25max": [1000], "0.5max": [2000], "0.75max": [3000], "1max": [4000]}
    assert embryo_b["validation_target"] == "embryo_E7.5" and embryo_b["override"] is None
    assert embryo_b["record"]["incomplete"] == ["1max lacks seeds [9]", "0.25max lacks seeds [9]",
                                                "0.5max lacks seeds [3, 9]", "0.75max lacks seeds [9]"]
    assert res["rows_replaced"] == 1 and len(res["inputs_sha256"]) == 29
    assert res["approved"] == "2026-10-21" and res["rule"] == tb.COUNT_RULE
    assert res["allow_uncommitted"] is False and res["rule_files_uncommitted"] == []
    assert res["git_commit"] == "synthetic" and set(res["rule_files_sha256"]) == set(ds.RULE_FILES)
    assert res["earlier_decisions"] == []
    printed = capsys.readouterr().out
    assert "heart board (heart_E10.5): override '0.75max'; use 0.75max" in printed
    assert "embryo board (embryo_E7.5): override None; keep max_cells" in printed
    assert "WARNING" not in printed
    with pytest.raises(SystemExit, match="exists"):  # the first decision is the record
        ds.decide(raw, out)
    second = ds.decide(raw, out, overwrite=True)
    assert second["boards"]["heart"]["override"] == "0.75max"
    kept = sorted(out.glob("count_decision.replaced-*.json"))  # --overwrite keeps the replaced record
    assert len(kept) == 1 and json.loads(kept[0].read_text()) == res
    assert second["earlier_decisions"] == [{"path": str(out / "count_decision.json"),
                                            "sha256": hashlib.sha256(text.encode()).hexdigest(),
                                            "kept_as": str(kept[0])}]
    ds.decide(raw, out, overwrite=True)  # twice in a second: a second name, nothing lost
    assert len(list(out.glob("count_decision.replaced-*.json"))) == 2


def test_decide_mode_keeps_max_cells_when_the_rows_do_not_measure_the_rule_levels(ds, tmp_path, clean_code):
    """Rows that disagree on max_cells, a max_cells not from the validation board, or a level capped at the pool: the
    rows do not measure the rule's levels on the validation board, so max_cells stays (whatever the scores say)."""
    def decide_heart(name, edit):
        raw = tmp_path / name
        raw.mkdir()
        _write_raw(raw, "heart_E10.5", _rerun({**LOSERS, "0.75max": 0.5}))  # 0.75max passes on clean rows
        for p in raw.glob("*.jsonl"):
            rows = [json.loads(line) for line in p.read_text().splitlines()]
            p.write_text("".join(json.dumps(edit(r)) + "\n" for r in rows))
        return ds.decide(raw, tmp_path / f"{name}_out")["boards"]["heart"]

    def capped(r):  # a pool of 3000 cells: 1max draws 3000, short of max_cells
        return {**r, "pool_cells": 3000, "n_pred": min(r["n_pred"], 3000)}

    def mixed(r):  # seeds 5-9 ran with --max-cells 6000
        return {**r, "max_cells": 6000, "max_cells_source": "cli"} if r["seed"] >= 5 else r

    def ruled(r):  # max_cells from the documented cap rule, not the board
        return {**r, "max_cells_source": "rule"}

    def short(r):  # 0.75max drew as many cells as 1max: not a lower count
        return {**r, "n_pred": 4000} if r["condition"].endswith("_0.75max") else r

    for name, edit, problem in (
            ("capped", capped, "1max drew [3000] cells, not round(1 x 4000) = 4000 (pool [3000])"),
            ("mixed", mixed, "the rows disagree on max_cells or lack it ([4000, 6000])"),
            ("ruled", ruled, "max_cells came from ['rule'], not only from the validation board's published limit"),
            ("short", short, "0.75max drew [4000] cells, not round(0.75 x 4000) = 3000")):
        heart = decide_heart(name, edit)
        assert heart["override"] is None and heart["record"]["override"] == "0.75max", name
        assert any(p.startswith(problem) for p in heart["row_problems"]), (name, heart["row_problems"])
        assert heart["decision"].startswith("keep max_cells: the rows do not measure approved rule 7's levels")
        assert heart["decision"].endswith("(count_override on these rows alone gives '0.75max')")
    assert decide_heart("clean", lambda r: r)["override"] == "0.75max"


def test_decide_mode_refuses_uncommitted_rule_code_unless_allowed(ds, tmp_path, monkeypatch, capsys, release_json):
    real_state = ds.rule_code_state
    raw, out = tmp_path / "raw", tmp_path / "out"
    raw.mkdir()
    _write_raw(raw, "heart_E10.5", _rerun({**LOSERS, "0.75max": 0.5}))
    for dirty in (["src/vec/t2_build.py"], None):  # None: git cannot tell
        monkeypatch.setattr(ds, "rule_code_state", lambda repo, d=dirty: {**CLEAN_CODE, "rule_files_uncommitted": d})
        with pytest.raises(SystemExit, match="uncommitted changes to the files that apply approved rule 7"):
            ds.decide(raw, out)
        assert not out.exists()  # nothing written
    monkeypatch.setattr(ds, "rule_code_state",
                        lambda repo: {**CLEAN_CODE, "rule_files_uncommitted": ["src/vec/t2_build.py"]})
    res = ds.decide(raw, out, allow_uncommitted=True)
    assert res["allow_uncommitted"] is True and res["rule_files_uncommitted"] == ["src/vec/t2_build.py"]
    assert res["boards"]["heart"]["override"] == "0.75max"
    assert "WARNING (allowed, recorded): uncommitted rule code" in capsys.readouterr().out
    # the real reader: this tree's git state, whatever it is, and the sha256 of the rule files
    state = real_state(REPO)
    assert state["rule_files_sha256"] == {p: hashlib.sha256((REPO / p).read_bytes()).hexdigest()
                                          for p in ds.RULE_FILES}
    assert state["rule_files_uncommitted"] is None or set(state["rule_files_uncommitted"]) <= set(ds.RULE_FILES)


def test_decide_mode_another_out_does_not_sidestep_the_release_decision(ds, tmp_path, capsys, clean_code,
                                                                         release_json):
    """The first decision is the record: once the release decision exists, a decide with another --out (the rerun run
    again, say) is refused too, as experiments/t1_direction/readout.py refuses; with --overwrite the release decision
    stays where it is, unchanged, and the new decision names it."""
    raw = tmp_path / "raw"
    raw.mkdir()
    _write_raw(raw, "heart_E10.5", _rerun({**LOSERS, "0.75max": 0.5}))
    first = ds.decide(raw, release_json.parent)  # the release command's --out
    text = release_json.read_text()
    sha = hashlib.sha256(text.encode()).hexdigest()
    assert json.loads(text) == first and first["earlier_decisions"] == []
    again = tmp_path / "rerun"
    with pytest.raises(SystemExit, match="the release decision .* exists .* another --out does not sidestep it"):
        ds.decide(raw, again)
    assert not again.exists()  # nothing written
    capsys.readouterr()
    second = ds.decide(raw, again, overwrite=True)
    assert release_json.read_text() == text and not list(release_json.parent.glob("*.replaced-*"))  # left in place
    assert json.loads((again / ds.DECISION_FILE).read_text()) == second
    assert second["earlier_decisions"] == [{"path": str(release_json), "sha256": sha, "kept_as": str(release_json)}]
    assert f"WARNING (--overwrite): the release decision {release_json} stays where it is" in capsys.readouterr().out
    # the release directory named another way is the same record: refused as the same --out, and replaced (kept next
    # to the new one) with --overwrite, listed once
    alias = release_json.parent / ".." / release_json.parent.name
    with pytest.raises(SystemExit, match=r"exists \(the first decision is the record\); pass --overwrite"):
        ds.decide(raw, alias)
    third = ds.decide(raw, alias, overwrite=True)
    kept = list(release_json.parent.glob("count_decision.replaced-*.json"))
    assert len(kept) == 1 and kept[0].read_text() == text and json.loads(release_json.read_text()) == third
    assert len(third["earlier_decisions"]) == 1 and third["earlier_decisions"][0]["sha256"] == sha
    assert Path(third["earlier_decisions"][0]["kept_as"]).resolve() == kept[0].resolve()


def test_decide_mode_refuses_partial_rows_and_a_directory_without_rule_rows(ds, tmp_path, clean_code):
    raw, out = tmp_path / "raw", tmp_path / "out"
    raw.mkdir()
    _write_raw(raw, "heart_E10.5", _rerun({**LOSERS, "0.75max": 0.5}))
    p = raw / "heart_E10.5_s3.jsonl"
    text = p.read_text()
    p.write_text(text[:-100])  # a killed job: the last row is cut short
    n_lines = len(text.splitlines())
    with pytest.raises(SystemExit, match=rf"heart_E10.5_s3.jsonl line {n_lines} is not a raw row"):
        ds.decide(raw, out)
    p.write_text(text + "42\n")  # a line that is JSON but not a row
    with pytest.raises(SystemExit, match=rf"line {n_lines + 1} is not a raw row \(not a JSON object\)"):
        ds.decide(raw, out)
    assert not out.exists()
    v2 = tmp_path / "raw_v2"  # the pre-release directory: pseudo-target rows only
    v2.mkdir()
    _write_raw(v2, "heart_E9.5", _rerun({**LOSERS, "0.75max": 0.5}))
    with pytest.raises(SystemExit, match="no shuffle50 rows at the max_cells levels for heart_E10.5 or embryo_E7.5"):
        ds.decide(v2, out)
    assert not out.exists()


def test_decide_mode_from_the_command_line_writes_only_the_decision(ds, tmp_path, monkeypatch, capsys, clean_code):
    raw, out = tmp_path / "raw", tmp_path / "out"
    raw.mkdir()
    _write_raw(raw, "heart_E10.5", _rerun({**LOSERS, "0.75max": 0.5}), seeds=[0])
    monkeypatch.setattr(sys, "argv", ["summarise.py", "--raw", str(raw), "--out", str(out), "--decide"])
    ds.main()
    assert [p.name for p in out.iterdir()] == ["count_decision.json"]  # no summary tables or plots
    boards = json.loads((out / "count_decision.json").read_text())["boards"]
    assert boards["heart"]["override"] is None and boards["embryo"]["override"] is None  # the rerun did not run
    assert boards["embryo"]["record"]["incomplete"] == ["no 1max scores", "no 0.25max scores", "no 0.5max scores",
                                                        "no 0.75max scores"]
    assert boards["embryo"]["rows"] == 0 and "WARNING: no rows for embryo_E7.5" in capsys.readouterr().out
    fresh = tmp_path / "fresh"
    monkeypatch.setattr(sys, "argv", ["summarise.py", "--raw", str(tmp_path / "nope"), "--out", str(fresh),
                                      "--decide"])
    with pytest.raises(SystemExit, match="no raw directory"):
        ds.main()
    for argv in (["--out", str(fresh)], ["--raw", str(raw)], []):  # the defaults are the pre-release run's
        monkeypatch.setattr(sys, "argv", ["summarise.py", *argv, "--decide"])
        with pytest.raises(SystemExit) as e:
            ds.main()
        assert e.value.code == 2 and "--decide needs" in capsys.readouterr().err
    assert not fresh.exists()  # a refused decide leaves no empty directory behind


def _decision(heart="0.75max", embryo=None, gains=None):
    """A count decision as decide mode writes it, reduced to what n_cells_from_decision reads. `gains` gives a board's
    passing levels and their mean gains ({kind: {level: gain}}); by default the board's override alone passes."""
    boards = {}
    for (kind, target), lv in zip(tb.COUNT_RULE_TARGETS.items(), (heart, embryo), strict=True):
        g = (gains or {}).get(kind, {} if lv is None else {lv: 0.5})
        per_level = {x: {"mean_gain": g.get(x, -1.0), "seeds_ahead": 10 if x in g else 0, "passes": x in g}
                     for x in tb.COUNT_RULE_LEVELS}
        boards[kind] = {"validation_target": target, "override": lv, "decision": f"decision for {kind}",
                        "record": {"complete": True, "passing": [x for x in tb.COUNT_RULE_LEVELS if x in g],
                                   "per_level": per_level}}
    return {"rule": tb.COUNT_RULE, "approved": tb.RULES_APPROVED, "git_commit": "synthetic", "allow_uncommitted": False,
            "boards": boards}


HEART_TEST = {"key": "T2:heart:test_extrap", "min_cells": 100, "max_cells": 3000, "truth_cells": 800}
EMBRYO_TEST = {"key": "T2:embryo:test_interp", "min_cells": 100, "max_cells": 4000, "truth_cells": 600}


def test_n_cells_from_decision_applies_the_decision_file_by_code(ds, tmp_path, clean_code):
    raw, out = tmp_path / "raw", tmp_path / "release"
    raw.mkdir()
    _write_raw(raw, "heart_E10.5", _rerun({**LOSERS, "0.75max": 0.5}))  # the heart board: 0.75max
    _write_raw(raw, "embryo_E7.5", _rerun({**LOSERS, "0.75max": -0.5}))  # the embryo board: max_cells
    ds.decide(raw, out)
    path = out / "count_decision.json"
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    n, rec = tb.n_cells_from_decision(t2_board(**HEART_TEST), 5000, path)
    assert n == 2250 and rec["rule"] == "override '0.75max'"
    assert rec["count_decision"] == {
        "decision_file": str(path), "decision_sha256": sha, "kind": "heart", "validation_target": "heart_E10.5",
        "override": "0.75max", "ranking": ["0.75max"],
        "decision": json.loads(path.read_text())["boards"]["heart"]["decision"],
        "decision_git_commit": "synthetic", "decision_allowed_uncommitted": False}
    assert "count_floor" not in rec  # 2250 cells: within the test board's limits, used as decided
    assert (n, rec) == tb.n_cells_from_decision(t2_board(**HEART_TEST), 5000, str(path))
    n, rec = tb.n_cells_from_decision(t2_board(**EMBRYO_TEST), 5000, path)
    assert n == 4000 and rec["rule"] == "max_cells" and rec["count_decision"]["validation_target"] == "embryo_E7.5"
    assert rec["count_decision"]["override"] is None and rec["count_decision"]["ranking"] == []
    # the same as n_cells with the override copied by hand, plus the record of where it came from
    want_n, want = tb.n_cells(t2_board(**HEART_TEST), 5000, override="0.75max")
    n, rec = tb.n_cells_from_decision(t2_board(**HEART_TEST), 5000, _decision())
    assert n == want_n and {k: v for k, v in rec.items() if k != "count_decision"} == want
    assert rec["count_decision"]["decision_file"] is None and rec["count_decision"]["decision_sha256"] is None


def test_read_count_decision_refuses_what_is_not_the_boards_own_decision():
    heart, unnamed = t2_board(**HEART_TEST), t2_board(**{**HEART_TEST, "key": "T2:test"})
    assert tb.read_count_decision(_decision(), heart)[0] == "0.75max"
    assert tb.read_count_decision(_decision(), unnamed, kind="heart")[0] == "0.75max"  # a key naming neither kind
    for board, kind, match in ((unnamed, None, "names neither heart nor embryo"),
                               (heart, "embryo", "does not agree with kind"),  # a heart board never takes the embryo's
                               (t2_board(key="T2:heart:embryo"), None, "does not agree"),
                               (heart, "lung", "kind must be one of"),
                               (t2_board(key="T3:heart", task="T3"), None, "T2 boards only")):
        with pytest.raises(ValueError, match=match):
            tb.read_count_decision(_decision(), board, kind)
    swapped = _decision()
    swapped["boards"]["heart"]["validation_target"] = "embryo_E7.5"
    for bad, match in ((swapped, "decided by its own validation target"),
                       (_decision(heart="0.57max"), "neither None"),  # a typo
                       (_decision(heart="1max"), "neither None"),  # not a lower level
                       ({**_decision(), "rule": "keep max_cells unless ..."}, "another rule text"),
                       ({**_decision(), "approved": "2026-10-20"}, "another rule text"),
                       ({"boards": []}, "not a count decision"),
                       ({**_decision(), "boards": {"embryo": {}}}, "no decision for the heart")):
        with pytest.raises(ValueError, match=match):
            tb.read_count_decision(bad, heart)


def test_read_count_decision_refuses_a_record_its_override_does_not_head():
    heart = t2_board(**HEART_TEST)
    two = {"heart": {"0.5max": 0.9, "0.75max": 0.5}}
    assert tb.read_count_decision(_decision(heart="0.5max", gains=two), heart)[1]["ranking"] == ["0.5max", "0.75max"]
    no_record = _decision()
    del no_record["boards"]["heart"]["record"]
    mislisted = _decision(gains={"heart": {"0.5max": 0.4, "0.75max": 0.5}})
    mislisted["boards"]["heart"]["record"]["passing"] = ["0.75max"]  # 0.5max passes too
    incomplete = _decision()
    incomplete["boards"]["heart"]["record"]["complete"] = False

    def with_numbers(level, **numbers):  # a decision whose record carries other numbers for one level
        d = _decision(gains={"heart": {"0.75max": 0.5}})
        d["boards"]["heart"]["record"]["per_level"][level].update(numbers)
        return d

    for bad, match in ((_decision(heart="0.75max", gains=two), "not the first of its passing levels"),
                       (_decision(heart="0.75max", gains={"heart": {"0.5max": 0.9}}), "not the first"),  # not passing
                       (no_record, "no record of its passing levels"),
                       (incomplete, "the rerun is incomplete"),
                       (mislisted, "per-level verdicts pass"),
                       (with_numbers("0.75max", seeds_ahead=7), "disagrees with its numbers"),  # marked passing
                       (with_numbers("0.75max", mean_gain=0.3), "disagrees with its numbers"),  # not above 0.3
                       (with_numbers("0.5max", mean_gain=0.9, seeds_ahead=10), "disagrees"),  # marked failing, passes
                       (_decision(gains={"heart": {"0.75max": 0.2}}), "disagrees with its numbers"),
                       (_decision(gains={"heart": {"0.75max": float("nan")}}), "no mean gain, seed count"),
                       (_decision(gains={"heart": {"0.75max": True}}), "no mean gain, seed count"),
                       (_decision(gains={"heart": {"0.75max": 10 ** 400}}), "no mean gain, seed count"),
                       (with_numbers("0.75max", seeds_ahead=11), "seed count \\(0-10\\)"),  # more than the 10 seeds
                       (with_numbers("0.75max", seeds_ahead=10.0), "no mean gain, seed count"),
                       (with_numbers("0.75max", seeds_ahead=True), "no mean gain, seed count"),
                       (with_numbers("0.25max", seeds_ahead=None), "no mean gain, seed count"),
                       (with_numbers("0.25max", passes=None), "and verdict for 0.25max: -1.0, 0, None")):
        with pytest.raises(ValueError, match=match):
            tb.read_count_decision(bad, heart)
    assert tb.read_count_decision(with_numbers("0.75max", seeds_ahead=8), heart)[1]["ranking"] == ["0.75max"]
    typed = with_numbers("0.75max", mean_gain=np.float32(0.5), seeds_ahead=np.int64(9), passes=np.True_)
    typed["boards"]["heart"]["record"]["complete"] = np.True_  # numpy numbers and bools, as in count_override's inputs
    assert tb.read_count_decision(typed, heart)[1]["ranking"] == ["0.75max"]
    keep = _decision(heart=None)
    del keep["boards"]["heart"]["record"]  # keeping max_cells reads no passing levels
    assert tb.read_count_decision(keep, heart)[1]["ranking"] == []


def test_n_cells_from_decision_takes_the_next_passing_level_below_the_test_boards_lower_limit():
    """COUNT_FLOOR_RULE (the user approved "the next passing count the board allows, or max_cells if none does" on
    2026-10-21 and left the reading of "next" to the AI assistant, which chose the next count up): a winning level that
    gives the test board fewer cells than its lower limit gives way to the passing level with the fewest cells that the
    board allows, whatever its gain, else max_cells. No level is clipped to the limit."""
    board = t2_board(key="T2:heart:test_extrap", min_cells=100, max_cells=6000, truth_cells=2000)
    # 0.25max wins with 1500 cells, below truth_cells: the next count up is 0.5max, although 0.75max measured higher
    gains = {"heart": {"0.25max": 0.9, "0.5max": 0.4, "0.75max": 0.6}}
    n, rec = tb.n_cells_from_decision(board, 9000, _decision(heart="0.25max", gains=gains))
    assert n == 3000 and rec["rule"] == "override '0.5max'" and rec["count_decision"]["override"] == "0.25max"
    assert rec["count_decision"]["ranking"] == ["0.25max", "0.75max", "0.5max"]
    assert rec["count_floor"] == {
        "rule": tb.COUNT_FLOOR_RULE, "floor": 2000, "floor_source": "truth_cells",
        "below_floor": [{"level": "0.25max", "n": 1500}], "level": "0.5max",
        "decision": "0.25max gives 1500 cells, below truth_cells 2000: use 0.5max, the passing level with the fewest "
                    "cells the board allows"}
    gains = {"heart": {"0.25max": 0.9, "0.5max": 0.6, "0.75max": 0.4}}  # the other levels' gains play no part
    assert tb.n_cells_from_decision(board, 9000, _decision(heart="0.25max", gains=gains))[0] == 3000
    gains = {"heart": {"0.25max": 0.9, "0.75max": 0.4}}  # 0.5max did not pass: the next passing count up is 0.75max
    n, rec = tb.n_cells_from_decision(board, 9000, _decision(heart="0.25max", gains=gains))
    assert n == 4500 and rec["count_floor"]["level"] == "0.75max"
    # every passing level below the limit is recorded, in rule order: 1500 and 3000 cells are below 3500, 4500 is not
    high = t2_board(key="T2:heart:test_extrap", min_cells=100, max_cells=6000, truth_cells=3500)
    gains = {"heart": {"0.25max": 0.9, "0.5max": 0.8, "0.75max": 0.4}}
    n, rec = tb.n_cells_from_decision(high, 9000, _decision(heart="0.25max", gains=gains))
    assert n == 4500 and rec["count_floor"]["below_floor"] == [{"level": "0.25max", "n": 1500},
                                                               {"level": "0.5max", "n": 3000}]
    # no passing level is allowed: max_cells
    n, rec = tb.n_cells_from_decision(board, 9000, _decision(heart="0.25max", gains={"heart": {"0.25max": 0.9}}))
    assert n == 6000 and rec["rule"] == "max_cells" and rec["count_floor"]["level"] is None
    assert rec["count_floor"]["decision"].endswith("no passing level is allowed, so keep max_cells")
    # without truth_cells the lower limit is min_cells
    no_truth = t2_board(key="T2:heart:test_extrap", min_cells=2000, max_cells=6000)
    gains = {"heart": {"0.25max": 0.9, "0.5max": 0.4}}
    n, rec = tb.n_cells_from_decision(no_truth, 9000, _decision(heart="0.25max", gains=gains))
    assert n == 3000 and rec["count_floor"]["floor"] == 2000 and rec["count_floor"]["floor_source"] == "min_cells"
    # a winning level the board allows is used as decided
    gains = {"heart": {"0.25max": 0.4, "0.5max": 0.9}}
    n, rec = tb.n_cells_from_decision(board, 9000, _decision(heart="0.5max", gains=gains))
    assert n == 3000 and rec["rule"] == "override '0.5max'" and "count_floor" not in rec
    # ... even when a passing level with fewer cells is allowed too: the count-up step acts only below the limit
    gains = {"heart": {"0.5max": 0.4, "0.75max": 0.9}}
    n, rec = tb.n_cells_from_decision(board, 9000, _decision(heart="0.75max", gains=gains))
    assert n == 4500 and rec["rule"] == "override '0.75max'" and "count_floor" not in rec
    # a winner below the limit is recorded first, then the other passing levels below it, in rule order
    high = t2_board(key="T2:heart:test_extrap", min_cells=100, max_cells=6000, truth_cells=3500)
    gains = {"heart": {"0.25max": 0.6, "0.5max": 0.9, "0.75max": 0.4}}
    n, rec = tb.n_cells_from_decision(high, 9000, _decision(heart="0.5max", gains=gains))
    assert n == 4500 and rec["count_floor"]["below_floor"] == [{"level": "0.5max", "n": 3000},
                                                               {"level": "0.25max", "n": 1500}]
    # two allowed levels that round to the same count: the lower fraction is used (the count is the same)
    tiny = t2_board(key="T2:heart:test_extrap", min_cells=2, max_cells=3)
    gains = {"heart": {"0.25max": 0.9, "0.5max": 0.5, "0.75max": 0.6}}
    n, rec = tb.n_cells_from_decision(tiny, 9000, _decision(heart="0.25max", gains=gains))
    assert n == 2 and rec["count_floor"]["level"] == "0.5max"
    # the limits are inclusive: a level giving exactly the lower limit is allowed, as decided or as the next one
    exact = t2_board(key="T2:heart:test_extrap", min_cells=100, max_cells=3000, truth_cells=750)
    n, rec = tb.n_cells_from_decision(exact, 9000, _decision(heart="0.25max", gains={"heart": {"0.25max": 0.9}}))
    assert n == 750 and "count_floor" not in rec
    exact = t2_board(key="T2:heart:test_extrap", min_cells=100, max_cells=3000, truth_cells=2250)
    gains = {"heart": {"0.25max": 0.9, "0.75max": 0.4}}
    n, rec = tb.n_cells_from_decision(exact, 9000, _decision(heart="0.25max", gains=gains))
    assert n == 2250 and rec["count_floor"]["level"] == "0.75max"
    # a board decide mode kept at max_cells (rows that do not measure the rule's levels) never falls back to a level
    gains = {"heart": {"0.25max": 0.9, "0.5max": 0.4}}
    n, rec = tb.n_cells_from_decision(board, 9000, _decision(heart=None, gains=gains))
    assert n == 6000 and rec["rule"] == "max_cells" and "count_floor" not in rec
    # other refusals of n_cells are passed on unchanged
    with pytest.raises(ValueError, match="no max_cells"):
        tb.n_cells_from_decision(t2_board(key="T2:heart:test_extrap", max_cells=None), 9000, _decision())
    with pytest.raises(ValueError, match="inconsistent limits"):
        tb.n_cells_from_decision(t2_board(key="T2:heart:test_extrap", min_cells=10, max_cells=500, truth_cells=800),
                                 9000, _decision())


def test_the_floor_rule_reads_the_decision_decide_mode_writes(ds, tmp_path, clean_code):
    raw, out = tmp_path / "raw", tmp_path / "release"
    raw.mkdir()
    # all three pass; 0.75max with exactly 8 of the 10 seeds ahead (mean gain 0.44)
    _write_raw(raw, "heart_E10.5", _rerun({"0.25max": 0.9, "0.5max": 0.4, "0.75max": [0.6] * 8 + [-0.2] * 2}))
    _write_raw(raw, "embryo_E7.5", _rerun({**LOSERS, "0.75max": -0.5}))
    ds.decide(raw, out)
    board = t2_board(key="T2:heart:test_extrap", min_cells=100, max_cells=6000, truth_cells=2000)
    n, rec = tb.n_cells_from_decision(board, 9000, out / "count_decision.json")
    assert n == 3000 and rec["count_decision"]["override"] == "0.25max"  # the next count up, not the next best gain
    assert rec["count_decision"]["ranking"] == ["0.25max", "0.75max", "0.5max"]
    assert rec["count_floor"]["below_floor"] == [{"level": "0.25max", "n": 1500}]
    assert rec["count_floor"]["level"] == "0.5max"
    assert json.loads((out / "count_decision.json").read_text())["boards"]["heart"]["record"]["per_level"]["0.75max"][
        "seeds_ahead"] == 8


# ---------------------------------------------------------------------------------------------------- sampling
def test_sampling_without_replacement_is_random_rows_in_random_order():
    pool = stack(2000)
    s = tb.sample_cells(pool, 500, seed=3)
    assert s.rows.shape == (500,) and np.unique(s.rows).size == 500 and s.rows.min() >= 0 and s.rows.max() < 2000
    assert not (np.all(np.diff(s.rows) > 0) or np.all(np.diff(s.rows) < 0))  # neither a block nor file order
    np.testing.assert_array_equal(s.coords, pool[s.rows])  # nothing to jitter
    assert s.record["with_replacement"] is False and s.record["n_jittered"] == 0 and s.record["max_shift"] == 0.0
    np.testing.assert_array_equal(tb.sample_cells(pool, 500, seed=3).rows, s.rows)  # deterministic per seed
    assert not np.array_equal(tb.sample_cells(pool, 500, seed=4).rows, s.rows)


def test_sampling_above_the_pool_repeats_rows_evenly_and_jitter_breaks_every_duplicate_within_the_bound():
    pool = stack(2000, seed=1)
    n = 5300
    s = tb.sample_cells(pool, n, seed=0)
    counts = np.bincount(s.rows, minlength=2000)
    assert set(np.unique(counts)) == {2, 3} and counts.sum() == n
    assert s.record["copies_per_row"] == [2, 3] and s.record["n_repeated_rows"] == n - 2000
    assert s.record["n_jittered"] == n - 2000 == s.record["dup_from_repeated_rows"]
    assert np.unique(s.coords, axis=0).shape[0] == n
    assert np.unique(s.coords.astype(np.float32), axis=0).shape[0] == n  # as the scorer reads them
    nn = tb.median_nn_distance(pool)
    assert s.record["median_nn"] == pytest.approx(nn)
    shift = np.linalg.norm(s.coords - pool[s.rows], axis=1)
    moved = shift > 0
    assert moved.sum() == n - 2000
    bound = tb.JITTER_NN_FRAC * nn
    assert shift.max() <= bound * (1 + 1e-12) and shift[moved].min() >= tb.JITTER_MIN_FRAC * bound * (1 - 1e-12)
    # the first occurrence of each row keeps its exact position
    _, first = np.unique(s.rows, return_index=True)
    np.testing.assert_array_equal(s.coords[first], pool[s.rows[first]])
    # isotropic: no axis is preferred
    d = (s.coords - pool[s.rows])[moved]
    assert np.all(np.abs(d.mean(0)) < 0.1 * bound)
    assert np.all(np.abs(np.abs(d).mean(0) / np.abs(d).mean() - 1) < 0.15)


def test_sampling_breaks_duplicate_positions_already_in_the_pool():
    pool = stack(1000, seed=2)
    pool[11] = pool[10]
    pool[500] = pool[10]
    pool[700] = pool[20] + np.array([0.0, 0.0, 1e-12])  # distinct in float64, one position as the scorer reads it
    assert np.unique(pool, axis=0).shape[0] == 998
    s = tb.sample_cells(pool, 1000, seed=0)
    assert s.record["dup_from_pool_positions"] == 3 and s.record["dup_from_repeated_rows"] == 0
    assert np.unique(s.coords.astype(np.float32), axis=0).shape[0] == 1000


def test_sampling_inputs_and_refusals():
    pool = stack(800, seed=3)
    a = ad.AnnData(X=np.zeros((800, len(GENES)), np.float32))
    a.obsm[COORD_KEY] = pool.astype(np.float32)
    st = Stage(X=np.zeros((800, len(GENES)), np.float32), coords=pool.astype(np.float32))
    rows = [tb.sample_cells(x, 300, seed=5).rows for x in (pool, a, st)]
    np.testing.assert_array_equal(rows[0], rows[1])
    np.testing.assert_array_equal(rows[0], rows[2])
    s = tb.sample_cells(pool, 2000, seed=1, jitter_frac=0.1)
    assert s.record["max_shift"] <= 0.1 * tb.median_nn_distance(pool) * (1 + 1e-12)
    with pytest.raises(ValueError, match="free budget"):
        tb.sample_cells(pool, 900, seed=0, jitter_frac=0.3)
    with pytest.raises(ValueError):
        tb.sample_cells(pool, 0, seed=0)
    with pytest.raises(ValueError):
        tb.sample_cells(pool, 10, seed=-1)
    with pytest.raises(ValueError, match="fewer than two distinct"):
        tb.sample_cells(np.zeros((1, 3)), 2, seed=0)
    with pytest.raises(KeyError):
        tb.sample_cells(ad.AnnData(X=np.zeros((5, 2), np.float32)), 2, seed=0)


# ---------------------------------------------------------------------------------------------------- coordinates
def test_centre_and_rms_radius():
    C = stack(500, seed=4) + np.array([1000.0, -50.0, 7.0])
    np.testing.assert_allclose(tb.centre(C).mean(0), 0.0, atol=1e-9)
    X = C - C.mean(0)
    assert tb.rms_radius(C) == pytest.approx(np.sqrt((X ** 2).sum(1).mean()))


def test_proper_rotation_rotates_about_the_centroid_and_refuses_improper_matrices():
    C = stack(400, seed=5) + 300.0
    R = random_rotation(0)
    D = tb.proper_rotation(C, R)
    np.testing.assert_allclose(D.mean(0), C.mean(0), atol=1e-9)
    np.testing.assert_allclose(np.linalg.norm(D[1:] - D[:-1], axis=1), np.linalg.norm(C[1:] - C[:-1], axis=1),
                               rtol=1e-12)
    assert tb.rms_radius(D) == pytest.approx(tb.rms_radius(C))
    np.testing.assert_allclose(tb.proper_rotation(C, R, about=np.zeros(3)), C @ R.T, atol=1e-9)
    mirror = np.diag([1.0, 1.0, -1.0])
    for improper in (mirror, R @ mirror, -np.eye(3)):
        with pytest.raises(ValueError, match="improper"):
            tb.proper_rotation(C, improper)
    with pytest.raises(ValueError, match="not orthonormal"):
        tb.proper_rotation(C, 2.0 * np.eye(3))
    with pytest.raises(ValueError):
        tb.check_rotation(np.eye(2))


def test_rescale_to_rms_is_proper_isotropic_and_matches_the_size_procedure(sp):
    C = stack(600, seed=6) + 40.0
    D = tb.rescale_to_rms(C, 321.0)
    assert tb.rms_radius(D) == pytest.approx(321.0)
    np.testing.assert_allclose(D.mean(0), C.mean(0), atol=1e-9)
    f = 321.0 / tb.rms_radius(C)
    np.testing.assert_allclose(D - D.mean(0), f * (C - C.mean(0)), atol=1e-9)
    np.testing.assert_allclose(D, sp.rescale_to_rms(C, 321.0), rtol=0, atol=1e-9)
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            tb.rescale_to_rms(C, bad)
    with pytest.raises(ValueError, match="no spread"):
        tb.rescale_to_rms(np.ones((5, 3)), 1.0)


# ---------------------------------------------------------------------------------------------------- size JSON
def _stage_cloud(scale, sd, seed):
    """A synthetic stage's cloud and core mask: a core blob inside a wider shell, scaled by `scale`, z on 15-unit
    sections."""
    rng = np.random.default_rng(seed)
    core = rng.normal(size=(1500, 3)) * np.asarray(sd) * 0.6
    shell = rng.normal(size=(1500, 3)) * np.asarray(sd)
    C = np.vstack([core, shell]) * scale
    C[:, 2] = np.round(C[:, 2] / 15.0) * 15.0
    return C, np.r_[np.ones(1500, bool), np.zeros(1500, bool)]


def _geom(sp, name, t, scale, sd, seed):
    """A synthetic stage's geometry, as the size procedure measures it (`_stage_cloud`)."""
    return sp.geometry_from_arrays(name, t, *_stage_cloud(scale, sd, seed))


def _size_json(sp, board, geoms, target_t, path):
    """What size_procedure.run writes, built from synthetic geometries (no files read)."""
    rep = sp.heart_report(geoms, target_t) if board == "heart" else sp.embryo_report(geoms, target_t)
    res = {"board": board, "target_t": float(target_t), "timestamp": "2026-10-19T00:00:00Z", "valid": rep["valid"],
           "preregistered_settings": dict(sp.PREREGISTERED),
           "provenance": {"git_commit": "synthetic", "git_uncommitted": [], "size_procedure_sha256": "synthetic"},
           "inputs": [], **rep}
    path.write_text(sp.dumps(res))
    return res


def _heart_stages(sd, seed=0):
    """The synthetic heart stages `_heart` measures, at t = 8.75, 9.5 and 10.5 with random growth, as
    (name, t, cloud, core mask): the clouds are the reference clouds a caller of size_or_hold passes."""
    rng = np.random.default_rng(seed)
    scales = np.cumprod(1.0 + rng.uniform(0.05, 0.6, size=3))  # random growth
    return [(f"s{i}", t, *_stage_cloud(s, sd, seed + i))
            for i, (t, s) in enumerate(zip((8.75, 9.5, 10.5), scales, strict=True))]


def _heart(sp, tmp_path, sd, name="size_heart.json", seed=0):
    geoms = [sp.geometry_from_arrays(*stage) for stage in _heart_stages(sd, seed)]
    path = tmp_path / name
    return _size_json(sp, "heart", geoms, 12.5, path), path


@pytest.mark.parametrize("hypothesis", ["H1", "H0", "H2"])
def test_size_raw_mode_is_an_isotropic_rescale_to_the_hypothesis_rms(sp, tmp_path, hypothesis):
    res, path = _heart(sp, tmp_path, sd=(30.0, 80.0, 100.0))  # core thin axis along x: auto resolves to raw
    assert res["preregistered"]["z_mode"] == "raw" and res["valid"] is True
    C = stack(900, seed=7)
    out, rec = tb.apply_size_procedure(C, path, hypothesis, expect_board="heart")
    target = res["preregistered"]["scope"][hypothesis]["pred_rms_full"]
    assert tb.rms_radius(out) == pytest.approx(target) and rec["target_rms_full"] == pytest.approx(target)
    np.testing.assert_allclose(out.mean(0), C.mean(0), atol=1e-9)
    np.testing.assert_allclose(out - out.mean(0), rec["scale"] * (C - C.mean(0)), atol=1e-9)
    assert rec["z_mode_recorded"] == rec["z_mode_used"] == "raw" and rec["z_held"] is False
    assert rec["size_json_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert rec["git_commit"] == "synthetic" and rec["warnings"] == []


def test_size_z_hold_mode_leaves_z_native(sp, tmp_path):
    res, path = _heart(sp, tmp_path, sd=(100.0, 80.0, 30.0))  # core thin axis along z: auto resolves to z hold
    pre = res["preregistered"]
    assert pre["z_mode"] == "xy_grow_z_hold"
    C = stack(900, seed=8)
    out, rec = tb.apply_size_procedure(C, path, "H2")
    np.testing.assert_array_equal(out[:, 2], C[:, 2])  # z untouched, bit for bit
    assert tb.rms_radius(out) == pytest.approx(pre["scope"]["H2"]["pred_rms_full"])
    assert rec["z_held"] is True and rec["target_rms_xy_full"] == pytest.approx(pre["scope"]["H2"]["pred_rms_xy_full"])
    out_xy, rec_xy = tb.apply_size_procedure(C, path, "H2", match="xy")
    X = out_xy - out_xy.mean(0)
    assert np.sqrt((X[:, :2] ** 2).sum(1).mean()) == pytest.approx(pre["scope"]["H2"]["pred_rms_xy_full"])
    np.testing.assert_array_equal(out_xy[:, 2], C[:, 2])
    # a logged alternative can be named for diagnosis; the record says which mode was used
    out_raw, rec_raw = tb.apply_size_procedure(C, path, "H1", z_mode="raw")
    assert rec_raw["z_mode_recorded"] == "xy_grow_z_hold" and rec_raw["z_mode_used"] == "raw"
    assert tb.rms_radius(out_raw) == pytest.approx(res["z_modes"]["raw"]["scope"]["H1"]["pred_rms_full"])


def test_size_embryo_has_no_hypothesis(sp, tmp_path):
    rng = np.random.default_rng(3)
    scales = np.cumprod(1.0 + rng.uniform(0.05, 0.6, size=4))
    geoms = [_geom(sp, f"e{i}", t, s, (60.0, 50.0, 90.0), 10 + i)
             for i, (t, s) in enumerate(zip((6.75, 7.25, 7.5, 8.0), scales, strict=True))]
    res = _size_json(sp, "embryo", geoms, 7.75, tmp_path / "size_embryo.json")
    C = stack(700, seed=9)
    out, rec = tb.apply_size_procedure(C, tmp_path / "size_embryo.json", expect_board="embryo")
    assert tb.rms_radius(out) == pytest.approx(res["preregistered"]["pred_rms_full"])
    assert rec["z_mode_used"] == "raw" and rec["reference"] == res["preregistered"]["bracket"]
    with pytest.raises(ValueError, match="no scope hypotheses"):
        tb.apply_size_procedure(C, tmp_path / "size_embryo.json", "H1")
    out_hold, _ = tb.apply_size_procedure(C, res, z_mode="xy_grow_z_hold")  # a dict works too
    np.testing.assert_array_equal(out_hold[:, 2], C[:, 2])
    assert tb.rms_radius(out_hold) == pytest.approx(res["z_modes"]["xy_grow_z_hold"]["pred_rms_full"])


def test_size_refusals(sp, tmp_path):
    res, path = _heart(sp, tmp_path, sd=(30.0, 80.0, 100.0))
    C = stack(300, seed=10)
    with pytest.raises(ValueError, match="scope hypothesis"):
        tb.apply_size_procedure(C, path)
    with pytest.raises(ValueError, match="scope hypothesis"):
        tb.apply_size_procedure(C, path, "H3")
    with pytest.raises(ValueError, match="diagnosis only"):
        tb.apply_size_procedure(C, path, "H1", z_mode="split")
    with pytest.raises(ValueError, match="cannot be applied"):
        tb.apply_size_procedure(C, path, "H1", z_mode="auto")
    with pytest.raises(ValueError, match="expected embryo"):
        tb.apply_size_procedure(C, path, "H1", expect_board="embryo")
    with pytest.raises(ValueError, match="only to the xy_grow_z_hold"):
        tb.apply_size_procedure(C, path, "H1", match="xy")
    with pytest.raises(ValueError, match="not a size procedure output"):
        tb.apply_size_procedure(C, res["preregistered"], "H1")
    # dry runs carry redacted forecasts: refused
    with pytest.raises(ValueError, match="redacted"):
        tb.apply_size_procedure(C, sp.redact_forecasts(res), "H1")
    # no growth rule applicable (no rate interval starts at or after E9.0): valid=false, refused
    early = [_geom(sp, "a", 8.25, 1.0, (30.0, 80.0, 100.0), 1), _geom(sp, "b", 8.75, 1.2, (30.0, 80.0, 100.0), 2)]
    invalid = _size_json(sp, "heart", early, 12.5, tmp_path / "invalid.json")
    assert invalid["valid"] is False
    with pytest.raises(ValueError, match="valid=false"):
        tb.apply_size_procedure(C, tmp_path / "invalid.json", "H1")
    # uncommitted procedure files are reported
    res["provenance"]["git_uncommitted"] = [" M experiments/coords/size_procedure.py"]
    assert "uncommitted" in tb.apply_size_procedure(C, res, "H1")[1]["warnings"][0]


# ---------------------------------------------------------------------------------------------------- size fallback
def _invalid_heart(sp, tmp_path, with_source=True, name="invalid_heart.json"):
    """A heart size JSON with valid=false (no rate interval starts at or after E9.0: there is no E9.5), from the
    procedure's own report. with_source adds a stage at E10.5 that the procedure measured. Returns the JSON, its path,
    and the E10.5 stage's full cloud and core mask (the fallback's reference cloud and a scope-filtered part of it)."""
    sd = (30.0, 80.0, 100.0)
    C, core = _stage_cloud(1.6, sd, seed=5)
    geoms = [_geom(sp, "a", 8.25, 1.0, sd, 1), _geom(sp, "b", 8.75, 1.2, sd, 2)]
    if with_source:
        geoms.append(sp.geometry_from_arrays("E10.5", 10.5, C, core))
    path = tmp_path / name
    res = _size_json(sp, "heart", geoms, 12.5, path)
    assert res["valid"] is False and res["preregistered"]["valid"] is False
    return res, path, C, core


def _embryo(sp, tmp_path):
    """A valid embryo size JSON (as in test_size_embryo_has_no_hypothesis) and its path."""
    rng = np.random.default_rng(3)
    scales = np.cumprod(1.0 + rng.uniform(0.05, 0.6, size=4))
    geoms = [_geom(sp, f"e{i}", t, s, (60.0, 50.0, 90.0), 10 + i)
             for i, (t, s) in enumerate(zip((6.75, 7.25, 7.5, 8.0), scales, strict=True))]
    path = tmp_path / "size_embryo.json"
    return _size_json(sp, "embryo", geoms, 7.75, path), path


@pytest.mark.parametrize("sd", [(30.0, 80.0, 100.0), (100.0, 80.0, 30.0)], ids=["raw", "z_hold"])
def test_size_or_hold_returns_what_apply_size_procedure_returns_for_a_valid_json(sp, tmp_path, sd):
    res, path = _heart(sp, tmp_path, sd=sd)
    assert res["valid"] is True
    name, t, ref, _ = _heart_stages(sd)[-1]  # the full E10.5 cloud, as the JSON measured it
    measured = next(g["rms_full"] for g in res["preregistered"]["stages"] if g["name"] == name)
    assert t == 10.5 and tb.rms_radius(ref) == pytest.approx(measured, rel=1e-12)
    C = stack(900, seed=7)
    for hypothesis in ("H1", "H2"):
        for size in (path, res):
            want = tb.apply_size_procedure(C, size, hypothesis, expect_board="heart")
            got = tb.size_or_hold(C, size, hypothesis, source_stage="E10.5", reference_cloud=ref, expect_board="heart")
            np.testing.assert_array_equal(got[0], want[0])
            assert got[1] == want[1] and "fallback" not in got[1]
    for kw in ({"z_mode": "raw"}, {"z_mode": "xy_grow_z_hold", "match": "xy"}):  # logged alternatives pass through
        want = tb.apply_size_procedure(C, path, "H1", **kw)
        got = tb.size_or_hold(C, path, "H1", source_stage="E10.5", reference_cloud=ref, **kw)
        np.testing.assert_array_equal(got[0], want[0])
        assert got[1] == want[1]


@pytest.mark.parametrize("valid", [True, False], ids=["valid", "invalid"])
@pytest.mark.parametrize("source", ["E9.5", "E8.5", "E10.5_val"])
def test_size_or_hold_returns_heart_cells_from_other_stages_to_the_captain(sp, tmp_path, valid, source):
    """Approved rule 3: if E10.5 is barred (A20 fails) the heart plan is rebuilt from E9.5 and the size rule comes back
    to the captain, so cells from another stage get no size, from a valid JSON (it measured E10.5) or an invalid one."""
    if valid:
        res, path = _heart(sp, tmp_path, sd=(30.0, 80.0, 100.0))
        ref = stack(3000, seed=8)
        assert res["valid"] is True and 10.5 in [g["t"] for g in res["preregistered"]["stages"]]
    else:
        _, path, ref, _ = _invalid_heart(sp, tmp_path)
    for h in ("H1", "H2"):
        with pytest.raises(ValueError, match=r"approved rule 3 \(2026-10-21\) sizes heart cells from E10.5 only") as e:
            tb.size_or_hold(stack(300, seed=4), path, h, source_stage=source, reference_cloud=ref,
                            expect_board="heart")
        msg = str(e.value)
        assert repr(source) in msg and "back to the captain" in msg and "A20" in msg and "No size is applied" in msg
        # a spelling of E10.5 (a file name) gets the hint; another stage does not
        assert ("pass source_stage='E10.5'" in msg) == (source == "E10.5_val")


def test_size_or_hold_passes_a_valid_embryo_json_through(sp, tmp_path):
    res, path = _embryo(sp, tmp_path)
    C = stack(700, seed=9)
    for size in (path, res):
        want = tb.apply_size_procedure(C, size, expect_board="embryo")
        got = tb.size_or_hold(C, size, expect_board="embryo")  # the embryo needs no source stage or reference cloud
        np.testing.assert_array_equal(got[0], want[0])
        assert got[1] == want[1]


@pytest.mark.parametrize("hypothesis", ["H1", "H0", "H2"])
def test_size_or_hold_keeps_the_e105_size_for_an_invalid_heart_json(sp, tmp_path, hypothesis):
    res, path, ref, core = _invalid_heart(sp, tmp_path)
    target = tb.rms_radius(ref)
    clouds = {"every lineage": tb.sample_cells(ref, 2000, seed=1).coords,  # an s1-like scope
              "core only": tb.sample_cells(ref[core], 900, seed=2).coords}  # an s2-like scope: a smaller cloud
    assert tb.rms_radius(clouds["core only"]) < 0.8 * target
    forecasts = [repr(v["pred_rms_full"]) for v in res["preregistered"]["scope"].values()]
    for name, C in clouds.items():
        out, rec = tb.size_or_hold(C, path, hypothesis, source_stage="E10.5", reference_cloud=ref,
                                   expect_board="heart")
        assert tb.rms_radius(out) == pytest.approx(target, rel=1e-12), name  # E10.5's size, whatever the scope
        np.testing.assert_allclose(out.mean(0), C.mean(0), atol=1e-9)
        np.testing.assert_allclose(out - out.mean(0), rec["scale"] * (C - C.mean(0)), atol=1e-9)  # isotropic
        assert rec["fallback"] is True and rec["fallback_rule"].startswith("approved rule 3")
        assert rec["fallback_approved"] == "2026-10-21" == tb.RULES_APPROVED
        assert rec["source_stage"] == "E10.5" and rec["hypothesis"] == hypothesis
        assert rec["hypothesis_changes_size"] is False and rec["z_held"] is False and rec["board"] == "heart"
        assert rec["target_rms_full"] == pytest.approx(target) and rec["rms_after"] == pytest.approx(target)
        assert rec["rms_before"] == pytest.approx(tb.rms_radius(C))
        assert rec["scale"] == pytest.approx(target / tb.rms_radius(C))
        assert rec["log_rms_after_vs_target"] == pytest.approx(0.0, abs=1e-12)
        why = rec["fallback_reason"]
        assert why["invalid_growth"] == {"core": {"applicable": []}} and why["min_rate_start_t"] == 9.0
        assert why["stages"] == {"a": 8.25, "b": 8.75, "E10.5": 10.5}
        assert rec["size_json_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert rec["git_commit"] == "synthetic" and rec["target_t"] == 12.5
        assert rec["size_json_stage_of_source"] == "E10.5"
        assert rec["size_json_rms_full_of_source"] == pytest.approx(target) and rec["warnings"] == []
        assert not any(f in json.dumps(rec) for f in forecasts)  # no forecast of the invalid JSON is copied


def test_size_or_hold_fallback_inputs_and_its_check_of_the_reference_cloud(sp, tmp_path):
    res, path, ref, core = _invalid_heart(sp, tmp_path)
    C = stack(800, seed=3)
    out, rec = tb.size_or_hold(C, path, "H1", source_stage="E10.5", reference_cloud=ref)
    a = ad.AnnData(X=np.zeros((len(ref), 2), np.float32))
    a.obsm[COORD_KEY] = ref
    np.testing.assert_array_equal(tb.size_or_hold(C, path, "H1", source_stage="E10.5", reference_cloud=a)[0], out)
    assert tb.size_or_hold(C, res, "H2", source_stage="E10.5", reference_cloud=ref)[1]["size_json_sha256"] is None
    # a reference cloud that is not the one the procedure measured as E10.5 (scope-filtered, rescaled, another
    # stage's) does not have E10.5's size: refused, never applied with a warning
    other, _ = _stage_cloud(1.2, (30.0, 80.0, 100.0), seed=2)
    for wrong in (ref[core], ref * 1.1, other):
        with pytest.raises(ValueError, match="is not the cloud the size procedure measured as E10.5") as e:
            tb.size_or_hold(C, res, "H2", source_stage="E10.5", reference_cloud=wrong)
        assert f"({len(wrong)} cells, RMS radius" in str(e.value) and "no size is applied" in str(e.value)
    # the procedure never measured E10.5: the layout of E10.5 barred (A20) or not the release-day JSON, refused
    bad2, path2, ref2, _ = _invalid_heart(sp, tmp_path, with_source=False, name="no_e105.json")
    assert bad2["preregistered"]["reference"] == "b"
    for reference in (ref2, ref, other):
        with pytest.raises(ValueError, match="the procedure did not measure E10.5") as e:
            tb.size_or_hold(C, path2, "H1", source_stage="E10.5", reference_cloud=reference)
        assert "{'a': 8.25, 'b': 8.75}" in str(e.value) and "back to the captain" in str(e.value)
    # provenance warnings, as apply_size_procedure gives them
    res["provenance"]["git_uncommitted"] = [" M experiments/coords/size_procedure.py"]
    assert "uncommitted" in tb.size_or_hold(C, res, "H1", source_stage="E10.5", reference_cloud=ref)[1]["warnings"][0]
    del res["provenance"]
    assert "no provenance" in tb.size_or_hold(C, res, "H1", source_stage="E10.5", reference_cloud=ref)[1]["warnings"][0]


def test_size_or_hold_refuses_the_a20_barred_layout(sp, tmp_path):
    """The release-day heart stages without the validation answers E8.5 and E10.5 (A20 fails): valid=false, reference
    E9.5. A stale source_stage='E10.5' with the E9.5 cloud gets no size: the size rule comes back to the captain."""
    sd = (30.0, 80.0, 100.0)
    e95, e95_core = _stage_cloud(1.5, sd, seed=4)
    geoms = [_geom(sp, "E8.25_late", 8.25, 1.0, sd, 1), _geom(sp, "E8.75", 8.75, 1.2, sd, 2),
             sp.geometry_from_arrays("E9.5", 9.5, e95, e95_core)]
    res = _size_json(sp, "heart", geoms, 12.5, tmp_path / "barred.json")
    assert res["valid"] is False and res["preregistered"]["reference"] == "E9.5"
    with pytest.raises(ValueError, match="the procedure did not measure E10.5"):
        tb.size_or_hold(stack(500, seed=6), res, "H1", source_stage="E10.5", reference_cloud=e95,
                        expect_board="heart")
    with pytest.raises(ValueError, match="sizes heart cells from E10.5 only"):
        tb.size_or_hold(stack(500, seed=6), res, "H1", source_stage="E9.5", reference_cloud=e95,
                        expect_board="heart")


def test_size_or_hold_checks_the_reference_cloud_of_a_valid_json_too(sp, tmp_path):
    """Every valid heart JSON measured E10.5, so the reference cloud is checked against that measurement as the
    fallback checks it: a stale source_stage='E10.5' with the E9.5 cloud (the heart plan rebuilt from E9.5 after A20
    fails), or a scope-filtered or rescaled E10.5 cloud, gets no size from a valid JSON either."""
    sd = (30.0, 80.0, 100.0)
    res, path = _heart(sp, tmp_path, sd=sd)
    (_, t95, e95, _), (_, t105, e105, core) = _heart_stages(sd)[1:]
    assert (t95, t105) == (9.5, 10.5) and res["valid"] is True
    C = stack(500, seed=6)
    for wrong in (e95, e105[core], e105 * 1.1):
        for size in (path, res):
            with pytest.raises(ValueError, match="is not the cloud the size procedure measured as E10.5") as e:
                tb.size_or_hold(C, size, "H1", source_stage="E10.5", reference_cloud=wrong, expect_board="heart")
            msg = str(e.value)
            assert f"({len(wrong)} cells, RMS radius" in msg and "not shown to come from E10.5" in msg
            assert "back to the captain" in msg and msg.endswith("no size is applied")
    for hypothesis in ("H1", "H2"):  # the E10.5 cloud itself: exactly what apply_size_procedure returns
        want = tb.apply_size_procedure(C, path, hypothesis, expect_board="heart")
        got = tb.size_or_hold(C, path, hypothesis, source_stage="E10.5", reference_cloud=e105, expect_board="heart")
        np.testing.assert_array_equal(got[0], want[0])
        assert got[1] == want[1]
    # apply_size_procedure's refusals still come first, whatever the reference cloud
    with pytest.raises(ValueError, match="redacted"):
        tb.size_or_hold(C, sp.redact_forecasts(res), "H1", source_stage="E10.5", reference_cloud=e95)
    # a valid JSON that never measured E10.5 is not the release-day JSON: refused, whatever the reference cloud
    geoms = [_geom(sp, f"x{i}", t, s, sd, 20 + i) for i, (t, s) in enumerate(((8.75, 1.0), (9.25, 1.3), (10.25, 1.6)))]
    other = _size_json(sp, "heart", geoms, 12.5, tmp_path / "valid_without_e105.json")
    assert other["valid"] is True
    for reference in (e105, e95):
        with pytest.raises(ValueError, match="is valid but the procedure did not measure E10.5") as e:
            tb.size_or_hold(C, other, "H1", source_stage="E10.5", reference_cloud=reference)
        assert "{'x0': 8.75, 'x1': 9.25, 'x2': 10.25}" in str(e.value) and "back to the captain" in str(e.value)


def test_size_or_hold_has_no_embryo_fallback(sp, tmp_path):
    res, path = _embryo(sp, tmp_path)
    res["valid"] = res["preregistered"]["valid"] = False
    path.write_text(sp.dumps(res))
    C = stack(700, seed=9)
    for size in (res, path):
        with pytest.raises(ValueError, match="no fallback was approved for the embryo"):
            tb.size_or_hold(C, size, expect_board="embryo")
        with pytest.raises(ValueError, match="no fallback was approved for the embryo"):
            tb.size_or_hold(C, size, source_stage="E10.5", reference_cloud=stack(500, seed=1))
    with pytest.raises(ValueError, match="valid=false"):
        tb.apply_size_procedure(C, path)  # unchanged


def test_size_or_hold_refuses_redacted_and_unknown_inputs_as_before(sp, tmp_path):
    res, path = _heart(sp, tmp_path, sd=(30.0, 80.0, 100.0))
    bad, path_bad, ref, _ = _invalid_heart(sp, tmp_path)
    C = stack(300, seed=10)
    kw = {"source_stage": "E10.5", "reference_cloud": ref}
    # a valid JSON: apply_size_procedure's refusals
    with pytest.raises(ValueError, match="redacted"):
        tb.size_or_hold(C, sp.redact_forecasts(res), "H1", **kw)
    with pytest.raises(ValueError, match="scope hypothesis"):
        tb.size_or_hold(C, path, "H3", **kw)
    with pytest.raises(ValueError, match="cannot be applied"):
        tb.size_or_hold(C, path, "H1", z_mode="auto", **kw)
    with pytest.raises(ValueError, match="not a size procedure output"):
        tb.size_or_hold(C, res["preregistered"], "H1", **kw)
    # an invalid JSON: the same refusals before any fallback
    with pytest.raises(ValueError, match="redacted"):
        tb.size_or_hold(C, sp.redact_forecasts(bad), "H1", **kw)
    for hypothesis in ("H3", None):
        with pytest.raises(ValueError, match="scope hypothesis"):
            tb.size_or_hold(C, path_bad, hypothesis, **kw)
    with pytest.raises(ValueError, match="expected embryo"):
        tb.size_or_hold(C, path_bad, "H1", expect_board="embryo", **kw)
    with pytest.raises(ValueError, match="not a size procedure output"):
        tb.size_or_hold(C, {"board": "heart", "valid": False}, "H1", **kw)
    with pytest.raises(ValueError, match="unknown size board"):
        tb.size_or_hold(C, {**bad, "board": "lung"}, "H1", **kw)
    for flag in (None, "false", 0):  # a valid flag that is neither true nor false is refused, never a fallback
        with pytest.raises(ValueError, match="valid=false"):
            tb.size_or_hold(C, {**bad, "valid": flag}, "H1", **kw)
    for extra in ({"z_mode": "xy_grow_z_hold"}, {"z_mode": "raw"}, {"match": "xy"}):
        with pytest.raises(ValueError, match="isotropic"):
            tb.size_or_hold(C, path_bad, "H1", **extra, **kw)
    # a heart call states the source stage and passes the reference cloud, valid JSON or not
    for p in (path, path_bad):
        with pytest.raises(ValueError, match="source_stage"):
            tb.size_or_hold(C, p, "H1", reference_cloud=ref)
        with pytest.raises(ValueError, match="reference_cloud"):
            tb.size_or_hold(C, p, "H1", source_stage="E10.5")
        with pytest.raises(ValueError, match=r"\(n, 3\)"):
            tb.size_or_hold(C, p, "H1", source_stage="E10.5", reference_cloud=np.zeros((5, 2)))
    with pytest.raises(ValueError, match="is not the cloud the size procedure measured"):  # no spread: not E10.5
        tb.size_or_hold(C, path_bad, "H1", source_stage="E10.5", reference_cloud=np.ones((5, 3)))


# ---------------------------------------------------------------------------------------------------- guard
def test_guard_passes_random_subsets_and_flags_contiguous_blocks():
    pool = stack(20000, seed=11)
    for m in (1000, 4000):
        g = tb.cloud_guard(tb.sample_cells(pool, m, seed=1).coords, pool, seed=0)
        assert g.ok, g.flags
        assert 0.9 < g.cover_ratio < 1.1 and abs(np.log(g.rms_ratio)) < 0.03
        mid = (len(pool) - m) // 2
        for block in (pool[mid:mid + m], pool[:m]):  # mid-file slab and the first rows (what a[:n] submits)
            gb = tb.cloud_guard(block, pool, seed=0)
            assert not gb.ok and any(f.startswith("block") for f in gb.flags)
            assert gb.cover_ratio < 0.6
    assert tb.cloud_guard(pool[:1000], pool, seed=0).as_dict()["ok"] is False


def test_guard_handles_translation_rescaling_rotation_and_oversampling():
    pool = stack(20000, seed=12)
    sub = tb.sample_cells(pool, 3000, seed=2).coords
    assert tb.cloud_guard(tb.centre(sub) + 500.0, pool).ok  # translation drops out
    big = tb.rescale_to_rms(sub, 1.8 * tb.rms_radius(sub))
    assert any(f.startswith("rms") for f in tb.cloud_guard(big, pool).flags)
    g = tb.cloud_guard(big, pool, expected_rms_ratio=tb.rms_radius(big) / tb.rms_radius(pool))
    assert g.ok and g.expected_rms_source == "argument"
    rot = tb.proper_rotation(sub, random_rotation(1))
    assert not tb.cloud_guard(rot, pool).ok  # frame="source" assumes no rotation
    assert tb.cloud_guard(rot, pool, frame="own").ok
    assert not tb.cloud_guard(pool[:3000], pool, frame="own").ok
    over = tb.sample_cells(pool, 26000, seed=3).coords  # copy-like above the pool size
    g = tb.cloud_guard(over, pool)
    assert g.ok and g.n_pred == 26000
    with pytest.raises(ValueError, match="frame"):
        tb.cloud_guard(sub, pool, frame="pca")
    with pytest.raises(ValueError, match="at least"):
        tb.cloud_guard(sub[:3], pool)
    with pytest.raises(ValueError):
        tb.cloud_guard(sub, pool, expected_rms_ratio=0.0)


def test_guard_weighted_reference_for_composition_weighted_draws():
    left = stack(10000, seed=15, extent=(60.0, 80.0, 60.0)) + np.array([-70.0, 0.0, 0.0])
    right = stack(10000, seed=16, extent=(60.0, 80.0, 60.0)) + np.array([70.0, 0.0, 0.0])
    pool = np.vstack([left, right])  # two cell types in two regions
    w = np.r_[np.ones(10000), np.full(10000, 0.03)]  # the right-hand type strongly down-weighted
    pred = pool[np.random.default_rng(17).choice(20000, 2000, replace=False, p=w / w.sum())]
    assert not tb.cloud_guard(pred, pool).ok  # uniform draws are the wrong reference for a weighted draw
    g = tb.cloud_guard(pred, pool, source_weights=w)
    assert g.ok and g.weighted_reference, g.flags
    mid = 10000 - 1000
    assert not tb.cloud_guard(pool[mid:mid + 2000], pool, source_weights=np.ones(20000)).ok
    with pytest.raises(ValueError, match="source_weights"):
        tb.cloud_guard(pred, pool, source_weights=w[:-1])
    with pytest.raises(ValueError, match="source_weights"):
        tb.cloud_guard(pred, pool, source_weights=-w)


# ---------------------------------------------------------------------------------------------------- builder
def _cells(n, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(rng.gamma(2.0, 1.0, size=(1, len(GENES))), size=(n, len(GENES)))
    lib = counts.sum(1, keepdims=True).clip(min=1)
    return np.log1p(1e4 * counts / lib).astype(np.float32)


def _counts(n, n_genes, seed=0):
    """MERFISH-like raw counts: a few per gene, a median library of about 150, no empty cell."""
    rng = np.random.default_rng(seed)
    mu = rng.gamma(0.5, 300.0 / n_genes, size=(1, n_genes))
    counts = rng.poisson(mu * rng.lognormal(0.0, 0.4, size=(n, 1)), size=(n, n_genes)).astype(np.float64)
    counts[counts.sum(1) == 0, 0] = 1.0
    return counts


def _cp10k(counts):
    return 1e4 * counts / counts.sum(1, keepdims=True)


def test_builder_orders_genes_and_passes_the_format_check(tmp_path, own_repo):
    board = t2_board(min_cells=10, max_cells=1000)
    X = _cells(200)
    order = np.random.default_rng(1).permutation(len(GENES))
    extra = np.ones((200, 2), np.float32)
    a = ad.AnnData(X=np.hstack([X[:, order], extra]))
    a.var_names = [GENES[i] for i in order] + ["Extra1", "Extra2"]
    coords = stack(200, seed=13)
    sub, rep = tb.build_submission(board, a, coords, out=tmp_path / "pred.h5ad")
    assert rep.ok and list(sub.var_names) == GENES
    np.testing.assert_array_equal(sub.X, X)
    assert sub.obsm[COORD_KEY].dtype == np.float32 and sub.obsm[COORD_KEY].shape == (200, 3)
    assert any("outside the board panel dropped" in w for w in rep.warnings)
    assert check_submission(tmp_path / "pred.h5ad", board).ok
    assert read_stamp(stamp_path(tmp_path / "pred.h5ad"))["clean"]  # stamped by default, against the test's repository
    back = read_stage(tmp_path / "pred.h5ad", GENES)
    np.testing.assert_array_equal(np.asarray(back.X), X)
    # matrix with gene names, Stage with gene names, matrix already in panel order
    assert tb.build_submission(board, X[:, order], coords, genes=[GENES[i] for i in order])[1].ok
    assert tb.build_submission(board, Stage(X=X[:, order]), coords, genes=[GENES[i] for i in order])[1].ok
    _, rep = tb.build_submission(board, X, coords)
    assert rep.ok and any("panel order" in w for w in rep.warnings)


def test_builder_refusals(tmp_path):
    board = t2_board(min_cells=50, max_cells=300)
    X, coords = _cells(100), stack(100, seed=14)
    with pytest.raises(ValueError, match="no gene names"):
        tb.build_submission(board, X[:, :-1], coords)
    with pytest.raises(KeyError):
        tb.build_submission(board, X[:, :-1], coords, genes=GENES[:-1])
    with pytest.raises(ValueError, match="do not pass genes"):
        tb.build_submission(board, ad.AnnData(X=X), coords, genes=GENES)
    with pytest.raises(ValueError, match="negative"):
        tb.build_submission(board, X - 1.0, coords)
    bad = X.copy()
    bad[0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        tb.build_submission(board, bad, coords)
    with pytest.raises(ValueError, match="coordinate rows"):
        tb.build_submission(board, X, coords[:99])
    with pytest.raises(ValueError, match="outside"):
        tb.build_submission(board, X[:40], coords[:40])
    with pytest.raises(ValueError, match="outside"):
        tb.build_submission(t2_board(min_cells=10, max_cells=99), X, coords)
    dup = coords.copy()
    dup[1] = dup[0]
    with pytest.raises(ValueError, match="repeat the exact"):
        tb.build_submission(board, X, dup)
    _, rep = tb.build_submission(board, X, dup, allow_duplicate_coords=True)
    assert rep.ok and any("exact coordinates" in w for w in rep.warnings)
    with pytest.raises(ValueError, match="takes no coordinates"):
        tb.build_submission(t2_board(task="T1", needs_coords=False), X, coords)
    counts = np.random.default_rng(0).poisson(3.0, size=(100, len(GENES))).astype(np.float32)
    _, rep = tb.build_submission(board, counts, coords)
    assert any("log1p CP10k" in w for w in rep.warnings)
    assert not list(tmp_path.iterdir())  # nothing written without `out`


def test_builder_stamps_its_file_but_never_the_temporary_one(tmp_path, monkeypatch):
    """The file written to `out` gets its build stamp (vec.io.write_submission); the temporary file written only to be
    checked, without `out`, gets none."""
    stamped = []
    monkeypatch.setattr(vec_io, "write_stamp", lambda path, repo=None: stamped.append(Path(path)))
    board, X, coords = t2_board(min_cells=10, max_cells=1000), _cells(200), stack(200, seed=13)
    assert tb.build_submission(board, X, coords, genes=GENES)[1].ok and stamped == []
    tb.build_submission(board, X, coords, genes=GENES, out=tmp_path / "pred.h5ad")
    assert stamped == [tmp_path / "pred.h5ad"]


# ---------------------------------------------------------------------------------------------------- end to end
def test_end_to_end_copy_like_submission_scores_under_both_builds(tmp_path, profiles, own_repo):
    _, later = profiles
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float)

    def draw(seed, n=4000):
        return simulate_stage(n, later, [0.3, 0.3, 0.4], np.random.default_rng(seed), coords_scale=3.0,
                              coords_shift=1.5 * centres)

    target, pool = draw(21), draw(22)
    base, _ = profiles
    ref = simulate_stage(4000, base, [0.5, 0.3, 0.2], np.random.default_rng(23), coords_scale=2.0,
                         coords_shift=centres)
    board = t2_board(min_cells=100, max_cells=1500, truth_cells=200)
    n, _ = tb.n_cells(board, pool.n)
    s = tb.sample_cells(pool, n, seed=0)
    assert tb.cloud_guard(s.coords, pool.coords, seed=0).ok
    tb.build_submission(board, pool.X[s.rows], s.coords, genes=GENES, out=tmp_path / "copy.h5ad")
    pred = Stage.from_anndata(read_stage(tmp_path / "copy.h5ad", GENES), GENES, need_coords=True)
    assert pred.n == 1500
    for build in SCORER_BUILDS:
        r = score_board(board, pred, target, ref, seed=0, build=build)
        assert r.extras["build"] == build and np.isfinite(r.score)
        assert r.score > 75.0, (build, r.group_skill)  # cells from the target's own distribution


# ---------------------------------------------------------------------------------------------------- no constants
# Every numeric literal allowed in vec/t2_build.py: indices and generic factors, the scorer's grid (16, 3.0), the
# jitter budget (0.25) and its lower end (0.5), guard thresholds (0.9, 0.03, 0.1), reference draws (5), the warning
# factor (2.0), the implied library band (0.1 and 2.0 x 1e4), the smallest cloud (4) and tolerances. A size, a count or
# a composition share of any stage fails.
ALLOWED_LITERALS = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 16.0, 0.25, 0.5, 0.9, 0.03, 0.1, 1e-6, 1e-9, 1e-12}
# The count rule's own thresholds as the user approved them on 2026-10-21 (docs/RELEASE_RULES_SUMMARY.md, answer 7: a
# paired mean gain above 0.3 points, at least 8 of the 10 release seeds ahead), each allowed only on the module-level
# line that defines its constant, so the same value anywhere else (a 0.3 composition share, say) still fails.
# The interpolation's method settings, fixed on 2026-11-01 before any result (feat-010; sizes of the computation, not
# values of any stage), are allowed the same way.
RULE_LITERALS = {"COUNT_RULE_SEEDS": {10.0}, "COUNT_RULE_MIN_GAIN": {0.3}, "COUNT_RULE_MIN_AHEAD": {8.0},
                 "INTERP_N_ANCHOR": {2000.0}, "INTERP_N_PC": {30.0}, "INTERP_K": {8.0}, "INTERP_MAX_ITER": {200.0},
                 "INTERP_EMD_ITER": {10_000_000.0},
                 # 2026-11-03, tries E7 and E8 (experiments/harness/TRIES.md): sizes of the computation, set before
                 # their run (a retry budget for thinned cells left without counts; the partner-sampling chunk)
                 "DEPTH_ZERO_RETRIES": {20.0}, "PARTNER_CHUNK": {4096.0},
                 "HERD_STREAMS": {6.0, 7.0},  # random-stream ids of the herded draws (try E9), not values of a stage
                 # 2026-11-09, try E13 (the midpoint draw): sizes of the computation and its fixed k-means seed, set
                 # before its run; MIDPOINT_STREAMS are random-stream ids, not values of a stage
                 "MIDPOINT_K_STRATA": {20.0}, "MIDPOINT_K_ANCHORS": {15.0}, "MIDPOINT_KMEANS_SEED": {20261109.0},
                 "MIDPOINT_STREAMS": {8.0, 9.0},
                 # 2026-11-09, try H6 (experiments/harness/TRIES.md): the smallest lineage the covariance filter acts
                 # on, a fixed setting, not a value of a stage
                 "CFM_MIN_CELLS": {200.0},
                 # 2026-11-12, tries E19a-E19c (experiments/harness/TRIES.md): sizes of the depth rendering's
                 # neighbour profile (principal components and neighbours), set before their run, not values of a stage
                 "RENDER_N_PC": {30.0}, "RENDER_K": {15.0},
                 # 2026-11-12, try E24 (the series rendering): the T3 script's gene_factor settings (rules.yaml,
                 # components.factor: eps 1e-5; its clip, 0.25-4, is made of allowed literals) and the smallest
                 # matched label, set before its run
                 "SERIES_MIN_CELLS": {100.0}, "SERIES_EPS": {1e-5}}


def _stray_literals(src: str) -> tuple[list[tuple[int, str]], set[str]]:
    """The numeric literals of `src` that neither ALLOWED_LITERALS nor RULE_LITERALS (on its constant's own line)
    allows, as (line, text), and the rule constants found on their lines."""
    bad, used = [], set()
    for t in tokenize.generate_tokens(io.StringIO(src).readline):
        if t.type != tokenize.NUMBER or float(t.string) in ALLOWED_LITERALS:
            continue
        m = re.match(r"([A-Z][A-Z0-9_]*) = ", t.line)  # a module-level constant's own line
        if m and float(t.string) in RULE_LITERALS.get(m[1], ()):
            used.add(m[1])
        else:
            bad.append((t.start[0], t.string))
    return bad, used


def test_module_holds_no_stage_constants():
    bad, used = _stray_literals((REPO / "src" / "vec" / "t2_build.py").read_text())
    assert not bad, bad
    assert used == set(RULE_LITERALS)  # each threshold is still defined where this test expects it


def test_the_literal_guard_catches_a_rule_value_used_elsewhere():
    src = (REPO / "src" / "vec" / "t2_build.py").read_text() + "\nSHARE = 0.3\n\n\ndef f(x):\n    return 8 * x\n"
    assert [text for _, text in _stray_literals(src)[0]] == ["0.3", "8"]  # a 0.3 share, an 8 outside the rule


def test_builder_enforces_the_truth_floor_unless_a_control_asks():
    board = t2_board(min_cells=10, max_cells=1000, truth_cells=300)
    X, coords = _cells(200), stack(200, seed=21)
    with pytest.raises(ValueError, match="below the board's truth_cells"):
        tb.build_submission(board, X, coords)
    assert tb.build_submission(board, X, coords, allow_below_truth=True)[1].ok


def test_no_raw_counts_warning_on_a_float32_single_gene_cell():
    board = t2_board(min_cells=10, max_cells=1000)
    X = _cells(50)
    X[0] = 0.0
    X[0, 0] = np.float32(np.log1p(1e4))  # a valid cell with all its counts in one gene, stored as float32
    _, rep = tb.build_submission(board, X, stack(50, seed=22))
    assert not any("raw counts" in w for w in rep.warnings)


def test_builder_warns_on_log1p_of_counts_never_scaled_to_1e4():
    board = t2_board(min_cells=10, max_cells=1000)
    counts, coords = _counts(300, len(GENES)), stack(300, seed=23)
    _, rep = tb.build_submission(board, np.log1p(counts), coords)  # fractional and below the ceiling: only size shows
    assert any("implied library" in w for w in rep.warnings), rep.warnings
    _, rep = tb.build_submission(board, np.log1p(_cp10k(counts)), coords)  # the same cells, scaled to 1e4 first
    assert not any("implied library" in w for w in rep.warnings), rep.warnings


def test_no_library_warning_on_released_like_cells_over_part_of_the_measured_panel():
    measured = [f"m{i}" for i in range(500)]
    board = Board(key="T2:synthetic", task="T2", genes=measured[:249] + measured[251:], needs_coords=True,
                  min_cells=10, max_cells=1000)  # 498 of the 500 measured genes
    X = np.log1p(_cp10k(_counts(300, 500, seed=1)))  # scaled over all 500 measured genes, as released
    _, rep = tb.build_submission(board, X, stack(300, seed=24), genes=measured)
    assert not any("implied library" in w or "raw counts" in w for w in rep.warnings), rep.warnings
    assert any("2 genes outside the board panel dropped" in w for w in rep.warnings)


def test_builder_warns_on_non_integer_values_above_the_log1p_ceiling():
    board = t2_board(min_cells=10, max_cells=1000)
    cp10k = _cp10k(_counts(300, len(GENES), seed=2))
    for X in (cp10k, np.log1p(100.0 * cp10k)):  # the log forgotten; log1p CPM
        nz = X[X > 0]
        assert not np.all(nz == np.rint(nz)) and X.max() > np.log1p(1e4)  # the max branch, not the whole-number one
        _, rep = tb.build_submission(board, X, stack(300, seed=25))
        raw = [w for w in rep.warnings if "raw counts" in w]
        assert raw and "whole numbers only: False" in raw[0], rep.warnings


# ---------------------------------------------------------------------------------------------------------------------
# Interpolation between two staged clouds (feat-010)
def test_pchip_rms_reads_the_embryo_size_rule_strictly_inside_the_stages():
    from scipy.interpolate import PchipInterpolator

    d, r = [6.75, 7.25, 8.0], [100.0, 150.0, 300.0]
    v = tb.pchip_rms(d, r, 7.5)
    assert 150.0 < v < 300.0 and v == pytest.approx(float(np.exp(PchipInterpolator(d, np.log(r))(7.5))))
    assert tb.pchip_rms(d[::-1], r[::-1], 7.5) == pytest.approx(v)  # stage order does not matter
    for bad in (6.75, 8.0, 9.0):
        with pytest.raises(ValueError, match="strictly inside"):
            tb.pchip_rms(d, r, bad)
    with pytest.raises(ValueError, match="share a day"):
        tb.pchip_rms([7.0, 7.0, 8.0], r, 7.5)


def _rot_z(a):
    return np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])


def test_interpolation_between_a_cloud_and_its_moved_copy_keeps_the_shape_and_lands_on_the_size():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(600, 3)) * [4.0, 2.0, 1.0] + 20.0
    E = rng.gamma(2.0, 1.0, size=(600, 12))
    perm = rng.permutation(600)
    B, EB = (2.0 * (A - A.mean(0)) @ _rot_z(0.6).T - 50.0)[perm], E[perm]  # rotated, doubled, moved, shuffled
    for t, stage in ((0.0, 0), (1.0, 1)):
        res = tb.interpolate_clouds(A, B, E, EB, t, 300, seed=1, target_rms=10.0, n_anchor=200, k=4)
        assert res.coords.shape == (300, 3) and (res.stage == stage).all()
        assert tb.rms_radius(res.coords) == pytest.approx(10.0)
    res = tb.interpolate_clouds(A, B, E, EB, 0.5, 400, seed=2, target_rms=10.0, n_anchor=600, k=4)  # all as anchors
    assert (res.stage == 1).sum() == 200 and res.record["rotation_deg"] == pytest.approx(np.degrees(0.6), abs=1.0)
    assert res.pool_coords.shape == (1200, 3) and res.pool_weights.sum() == pytest.approx(1.0)
    # same shape at both ends: every after cell lands where its own copy among the before cells lands
    pb, pa = res.pool_coords[:600], res.pool_coords[600:]
    assert np.abs(pa - pb[perm]).max() < 0.05 * tb.rms_radius(res.pool_coords)
    # placed at the before stage's centroid (up to the jitter of cells drawn twice: a before cell and its own copy)
    assert res.coords.mean(0) == pytest.approx(A.mean(0), abs=0.01 * tb.rms_radius(res.coords))
    assert np.unique(res.coords.astype(np.float32), axis=0).shape[0] == 400
    with pytest.raises(ValueError, match="t must"):
        tb.interpolate_clouds(A, B, E, EB, 1.5, 10, seed=0, target_rms=1.0)


def test_interpolation_moves_cells_part_of_the_way_between_two_shapes():
    rng = np.random.default_rng(4)
    A = rng.normal(size=(500, 3)) * [1.0, 1.0, 1.0]
    E = rng.gamma(2.0, 1.0, size=(500, 10))
    B = A * [3.0, 1.0, 1.0]  # the after shape is stretched along x; the same cells in the same order
    long_axis = []
    for t in (0.0, 0.5, 1.0):
        res = tb.interpolate_clouds(A, B, E, E, t, 500, seed=0, target_rms=1.0, n_anchor=250, k=4)
        C = res.coords - res.coords.mean(0)
        ev = np.sort(np.linalg.eigvalsh(np.cov(C.T)))[::-1]
        long_axis.append(ev[0] / ev[1])
    assert long_axis[0] < long_axis[1] < long_axis[2]  # elongation grows with t



# ------------------------------------------------------------------ heart s2 scope trim (docs/S2_RULES.md, 2026-11-02)
def _shelled_cloud(seed: int = 0, n_core: int = 300, n_shell: int = 900):
    rng = np.random.default_rng(seed)
    core = rng.normal(0, 1.0, (n_core, 3))
    shell = rng.normal(0, 3.0, (n_shell, 3)) + np.array([1.5, -0.5, 0.0])
    C = np.vstack([core, shell])
    return C, np.r_[np.ones(n_core, bool), np.zeros(n_shell, bool)]


@pytest.mark.parametrize("frac", [0.0, 0.25, 0.5, 1.0])
def test_radial_trim_mask_reaches_the_scope_with_fewest_removed(frac):
    C, core = _shelled_cloud()
    keep = tb.radial_trim_mask(C, core, frac)
    s0 = tb.xy_scope(C, core)
    assert keep[core].all()  # core cells always kept
    assert tb.xy_scope(C[keep], core[keep]) <= frac * s0 + 1e-12
    removed = np.flatnonzero(~keep)
    if len(removed):  # the farthest shell cells go first, and one fewer removal would miss the target
        d = np.linalg.norm(C[:, :2] - C[core, :2].mean(0), axis=1)
        assert d[removed].min() >= d[keep & ~core].max() - 1e-12
        back = keep.copy()
        back[removed[np.argmin(d[removed])]] = True
        assert tb.xy_scope(C[back], core[back]) > frac * s0
    if frac == 1.0:
        assert keep.all()


def test_radial_trim_mask_matches_the_payoff_study():
    spec = importlib.util.spec_from_file_location(
        "s2_payoff_for_trim", Path(__file__).resolve().parents[1] / "experiments" / "s2_hedge" / "payoff.py")
    payoff = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = payoff
    spec.loader.exec_module(payoff)
    for seed in range(3):
        C, core = _shelled_cloud(seed)
        for frac in (0.0, 0.5):
            ours = tb.radial_trim_mask(C, core, frac)
            theirs = payoff.trim_mask(C, core, frac, "radial", np.random.default_rng(0))
            assert np.array_equal(ours, theirs)


def test_radial_trim_mask_refuses_bad_input():
    C, core = _shelled_cloud()
    with pytest.raises(ValueError):
        tb.radial_trim_mask(C, np.zeros(len(C), bool), 0.5)
    with pytest.raises(ValueError):
        tb.radial_trim_mask(C, core, 1.5)


# ------------------------------------------------------------------ embryo s2 shape variants (t_pos), 2026-11-02
def test_interpolation_t_pos_sets_the_geometry_and_leaves_the_default_unchanged():
    rng = np.random.default_rng(4)
    A = rng.normal(size=(500, 3))
    E = rng.gamma(2.0, 1.0, size=(500, 10))
    B = A * [3.0, 1.0, 1.0]  # stretched along x; same cells, same order
    kw = dict(seed=0, target_rms=1.0, n_anchor=250, k=4)
    base = tb.interpolate_clouds(A, B, E, E, 0.5, 500, **kw)
    same = tb.interpolate_clouds(A, B, E, E, 0.5, 500, t_pos=0.5, **kw)
    assert np.array_equal(base.coords, same.coords) and np.array_equal(base.rows, same.rows)
    assert base.record["t_pos"] == 0.5
    elong = {}
    for tp in (0.0, 1.0):
        res = tb.interpolate_clouds(A, B, E, E, 0.5, 500, t_pos=tp, **kw)
        assert np.array_equal(res.rows, base.rows) and np.array_equal(res.stage, base.stage)  # the same mixture
        assert (res.stage == 1).sum() == 250 and tb.rms_radius(res.coords) == pytest.approx(1.0)
        C = res.coords - res.coords.mean(0)
        ev = np.sort(np.linalg.eigvalsh(np.cov(C.T)))[::-1]
        elong[tp] = ev[0] / ev[1]
    C = base.coords - base.coords.mean(0)
    ev = np.sort(np.linalg.eigvalsh(np.cov(C.T)))[::-1]
    assert elong[0.0] < ev[0] / ev[1] < elong[1.0]  # before's round shape, the middle, after's stretched shape
    with pytest.raises(ValueError, match="t_pos"):
        tb.interpolate_clouds(A, B, E, E, 0.5, 10, t_pos=1.2, **kw)


# ------------------------------------------------------------------- the moves' coupling (embryo try E10), 2026-11-03
def _left_right_pair(rng, n=500):
    """Two clouds whose expression is left-right symmetric (each gene a bump in |x|, the relative y and z, never in
    the sign of x), with x the thinnest axis as E8.0's left-right axis is; the after cloud is longer along y."""
    A = rng.normal(size=(n, 3)) * [0.5, 2.0, 1.0]
    B = rng.normal(size=(n, 3)) * [0.5, 3.0, 1.0]
    centres = rng.normal(size=(12, 3)) * [0.5, 1.0, 1.0]
    centres[:, 0] = np.abs(centres[:, 0])

    def expr(C, length):
        f = np.c_[np.abs(C[:, 0]), C[:, 1] / length, C[:, 2]]
        return np.exp(-((f[:, None, :] - centres[None]) ** 2).sum(-1))

    return A, B, expr(A, 2.0), expr(B, 3.0)


def _thinnest(C, w=None):
    """A cloud's smallest principal SD over its RMS radius (weights optional)."""
    w = np.full(len(C), 1.0 / len(C)) if w is None else np.asarray(w, float) / np.sum(w)
    X = C - w @ C
    ev = np.linalg.eigvalsh((X * w[:, None]).T @ X)
    return float(np.sqrt(ev.min() / ev.sum()))


def _xcorr_min(X, Y):
    """The smallest eigenvalue of the symmetrised whitened cross-correlation of paired rows of X and Y."""
    X, Y = X - X.mean(0), Y - Y.mean(0)

    def isqrt(S):
        ev, V = np.linalg.eigh(S)
        return V @ np.diag(ev ** -0.5) @ V.T

    K = isqrt(X.T @ X / len(X)) @ (X.T @ Y / len(X)) @ isqrt(Y.T @ Y / len(Y))
    return float(np.linalg.eigvalsh((K + K.T) / 2).min())


def test_interpolation_move_coupling_defaults_to_fgw_and_aligned_changes_only_the_moves():
    rng = np.random.default_rng(5)
    A, B, Eb, Ea = _left_right_pair(rng)
    kw = dict(seed=1, target_rms=1.0, n_anchor=200, k=4)
    base = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 400, **kw)
    same = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 400, move_coupling="fgw", **kw)
    assert np.array_equal(base.coords, same.coords) and base.record == same.record
    assert "move_coupling" not in base.record and base.coupling is None
    al = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 400, move_coupling="aligned", with_coupling=True, **kw)
    assert np.array_equal(al.rows, base.rows) and np.array_equal(al.stage, base.stage)  # the same mixture
    assert al.record["rotation_deg"] == base.record["rotation_deg"]  # FGW's rotation is kept
    assert al.record["coupling_nonzero"] == base.record["coupling_nonzero"]
    assert al.record["move_coupling"] == "aligned" and al.record["coupling_nonzero_moves"] == 200  # one to one
    c = al.coupling
    assert c["fgw"].shape == c["moves"].shape == (200, 200)
    np.testing.assert_allclose(c["after_frame"], ((B - B.mean(0)) / tb.rms_radius(B)) @ c["rotation"])
    P = np.vstack([c["before_frame"] + 0.5 * c["shift_before"], c["after_frame"] + 0.5 * c["shift_after"]])
    scale = al.record["scale_from_canonical"]
    np.testing.assert_allclose(al.pool_coords - al.pool_coords[0], (P - P[0]) * scale, atol=1e-9)
    with pytest.raises(ValueError, match="move_coupling"):
        tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 10, move_coupling="tissue", **kw)


def test_aligned_moves_keep_a_left_right_symmetric_interpolation_as_thick_as_its_brackets():
    rng = np.random.default_rng(6)
    A, B, Eb, Ea = _left_right_pair(rng)
    out = {}
    for mc in tb.MOVE_COUPLINGS:
        res = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 500, seed=2, target_rms=1.0, n_anchor=300, k=4,
                                    move_coupling=mc, with_coupling=True)
        c = res.coupling
        anchors_b = c["before_frame"][c["anchors_before"]]
        partners = (c["moves"] @ c["after_frame"][c["anchors_after"]]) / c["moves"].sum(1, keepdims=True)
        out[mc] = _xcorr_min(anchors_b, partners), _thinnest(res.pool_coords, res.pool_weights)
    # FGW cannot tell the sides apart: its partners are uncorrelated along the symmetric direction (the defect)
    assert out["fgw"][0] < 0.2
    # in the common frame position tells the side, and the moved cloud keeps its brackets' thickness
    assert out["aligned"][0] >= out["fgw"][0] + 0.25
    assert out["aligned"][1] > out["fgw"][1]
    assert out["aligned"][1] >= 0.9 * min(_thinnest(A), _thinnest(B))


# ------------------------------------------------------------------- detection depth of an interpolated stage (E7, E8)
def _counted(rng, n, depth, panel=12, extra=3, profile=None):
    """An AnnData as the released files store it: Poisson counts at a given depth on `panel` genes plus `extra` genes
    outside the panel, X = log(1 + 1e4 c / N) with N over every gene. Returns it, the panel gene names and N."""
    import pandas as pd
    from scipy import sparse as sp

    prof = rng.gamma(2.0, 1.0, size=panel + extra) if profile is None else np.asarray(profile, float)
    prof = prof / prof.sum()
    C = rng.poisson(prof * depth, size=(n, panel + extra))
    C[C[:, :panel].sum(1) == 0, 0] += 1  # every cell has a panel count, as every released cell does
    N = C.sum(1)
    genes = [f"g{i}" for i in range(panel + extra)]
    a = ad.AnnData(X=sp.csr_matrix(np.log1p(1e4 * C / N[:, None]).astype(np.float32)),
                   var=pd.DataFrame(index=genes))
    a.layers["counts"] = sp.csr_matrix(C.astype(np.int32))
    a.obsm[COORD_KEY] = rng.normal(size=(n, 3)).astype(np.float32)
    return a, genes[:panel], N


def test_library_sizes_recover_the_whole_library_from_the_panel_genes():
    rng = np.random.default_rng(30)
    a, panel, N = _counted(rng, 300, 80.0)
    sub = a[:, panel]
    lib = tb.library_sizes(sub.X, sub.layers["counts"])
    np.testing.assert_array_equal(lib, N)  # genes outside the panel are counted, as the released N counts them
    assert (lib > np.asarray(sub.layers["counts"].sum(1)).ravel()).any()
    with pytest.raises(ValueError, match="non-zero"):
        tb.library_sizes(sub.X, sub.layers["counts"][::-1])


def test_depth_target_is_log_linear_between_the_brackets():
    assert tb.depth_target(62.0, 212.0, 0.5) == pytest.approx(np.sqrt(62.0 * 212.0))
    assert tb.depth_target(35.0, 212.0, 0.0) == pytest.approx(35.0)
    assert tb.depth_target(35.0, 212.0, 1.0) == pytest.approx(212.0)
    with pytest.raises(ValueError):
        tb.depth_target(0.0, 10.0, 0.5)


def _interp_pair(rng, d_before, d_after, n=600, t=0.5, partners=False, **kw):
    before, panel, _ = _counted(rng, n, d_before, **kw)
    after, _, _ = _counted(rng, n, d_after, **kw)
    Eb, Ea = (np.asarray(x[:, panel].X.toarray(), np.float32) for x in (before, after))
    res = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, t, 500, seed=3,
                                target_rms=1.0, n_anchor=150, k=4, with_partners=partners)
    mix = np.empty((len(res.rows), len(panel)), np.float32)
    mix[res.stage == 0], mix[res.stage == 1] = Eb[res.rows[res.stage == 0]], Ea[res.rows[res.stage == 1]]
    return before, after, panel, res, mix


def test_depth_matching_thins_only_the_deeper_bracket_to_the_predicted_depth():
    rng = np.random.default_rng(31)
    before, after, panel, res, mix = _interp_pair(rng, 40.0, 200.0)
    X, rec = tb.depth_matched_expression(before, after, res, panel, 0.5, seed=[1, 7])
    b, a = res.stage == 0, res.stage == 1
    np.testing.assert_array_equal(X[b], mix[b])  # the shallower bracket keeps its released rows exactly
    assert rec["factor"]["before"] == 1.0 and rec["factor"]["after"] == pytest.approx(rec["target"] / 200, rel=0.05)
    assert rec["target"] == pytest.approx(np.sqrt(40.0 * 200.0), rel=0.08)
    assert not np.array_equal(X[a], mix[a])
    # the thinned cells' counts: as many detected genes as cells near the target depth, fewer than at 200
    assert (X[a] > 0).sum(1).mean() < (mix[a] > 0).sum(1).mean()
    assert (X[a] > 0).sum(1).min() >= 1 and np.isfinite(X).all() and X.max() <= np.log1p(1e4) + 1e-4
    # equal depths: nothing is thinned, so the result is the mixture
    before2, after2, panel2, res2, mix2 = _interp_pair(np.random.default_rng(32), 60.0, 60.0)
    X2, rec2 = tb.depth_matched_expression(before2, after2, res2, panel2, 0.5, seed=0)
    if rec2["factor"]["before"] == 1.0 and rec2["factor"]["after"] == 1.0:
        np.testing.assert_array_equal(X2, mix2)
    assert min(rec2["factor"].values()) > 0.9


def test_t_mix_changes_only_the_brackets_shares_and_defaults_to_t():
    """Embryo try E23: t_mix sets the later bracket's share of the drawn cells apart from t, which still sets the
    moves. Unset, the result is the entry's, bit for bit; at 1 every cell comes from the after stage."""
    rng = np.random.default_rng(51)
    before, after, panel, res, _ = _interp_pair(rng, 40.0, 200.0)
    Eb, Ea = (np.asarray(x[:, panel].X.toarray(), np.float32) for x in (before, after))
    kw = dict(seed=3, target_rms=1.0, n_anchor=150, k=4)
    same = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, t_mix=0.5, **kw)
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
        np.testing.assert_array_equal(getattr(res, f), getattr(same, f))
    late = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, t_mix=1.0, **kw)
    assert (late.stage == 1).all() and late.record["n_before"] == 0 and late.record["t_mix"] == 1.0
    assert late.record["t"] == 0.5 and late.record["t_pos"] == 0.5  # the moves are still the half-way moves
    assert late.pool_weights[: before.n_obs].sum() == 0.0 and late.pool_weights.sum() == pytest.approx(1.0)
    assert tb.rms_radius(late.coords) == pytest.approx(1.0, rel=1e-6)
    with pytest.raises(ValueError, match="t_mix"):
        tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, t_mix=1.5, **kw)


def _implied_library(X, C):
    """Median library size the released form implies for rows whose first panel gene is detected: 1e4 c / expm1(x)."""
    x, c = np.asarray(X, np.float64)[:, 0], np.asarray(C, np.float64)[:, 0]
    m = (x > 0) & (c > 0)
    return float(np.median(1e4 * c[m] / np.expm1(x[m])))


@pytest.mark.parametrize("mode", tb.RENDER_MODES)
def test_depth_rendering_moves_every_bracket_to_one_depth(mode):
    """Embryo tries E19a-E19c: mid renders both brackets to the log-linear depth, lo thins the deeper bracket to the
    shallower one's depth and hi tops the shallower one up to the deeper one's; a bracket already there is untouched."""
    rng = np.random.default_rng(41)
    before, after, panel, res, mix = _interp_pair(rng, 40.0, 200.0)
    X, rec = tb.depth_rendered_expression(before, after, res, panel, 0.5, seed=[1, 7], mode=mode)
    b, a = res.stage == 0, res.stage == 1
    d = rec["median_depth"]
    want = {"mid": np.sqrt(d["before"] * d["after"]), "lo": d["before"], "hi": d["after"]}[mode]
    assert rec["target"] == pytest.approx(want)
    assert np.isfinite(X).all() and (X >= 0).all() and X.max() <= np.log1p(1e4) + 1e-4
    assert ((X > 0).sum(1) >= 1).all()  # no cell is left without a panel count
    # each bracket's rendered median library size sits at the target (Poisson and binomial noise only)
    for name in ("before", "after"):
        assert rec["rendered_median"][name] == pytest.approx(want, rel=0.12)
    if mode == "lo":  # the shallower bracket keeps its released rows; the deeper one loses detected genes
        np.testing.assert_array_equal(X[b], mix[b])
        assert rec["topped_up"] == {"before": 0, "after": 0} and rec["thinned"]["after"] == int(a.sum())
        assert (X[a] > 0).sum(1).mean() < (mix[a] > 0).sum(1).mean()
    if mode == "hi":  # the deeper bracket keeps its released rows; the shallower one gains detected genes
        np.testing.assert_array_equal(X[a], mix[a])
        assert rec["thinned"] == {"before": 0, "after": 0} and rec["topped_up"]["before"] == int(b.sum())
        assert (X[b] > 0).sum(1).mean() > (mix[b] > 0).sum(1).mean()
    if mode == "mid":  # both move, in opposite directions
        assert rec["topped_up"]["before"] == int(b.sum()) and rec["thinned"]["after"] == int(a.sum())
        assert (X[b] > 0).sum(1).mean() > (mix[b] > 0).sum(1).mean()
        assert (X[a] > 0).sum(1).mean() < (mix[a] > 0).sum(1).mean()
    # the draws and positions are fgw's: only the expression step differs, and it is reproducible from its seed
    X2, _ = tb.depth_rendered_expression(before, after, res, panel, 0.5, seed=[1, 7], mode=mode)
    np.testing.assert_array_equal(X, X2)


def test_depth_rendering_top_up_follows_the_cells_own_profile_and_refuses_an_unknown_mode():
    """A topped-up cell's added counts follow its neighbours' proportions: genes its state does not express stay near
    zero. Two states with disjoint panel genes must not leak counts into each other."""
    rng = np.random.default_rng(42)
    pa = np.r_[np.ones(6), np.full(6, 1e-9), np.ones(3)]
    pb = np.r_[np.full(6, 1e-9), np.ones(6), np.ones(3)]
    s1, panel, _ = _counted(rng, 300, 40.0, profile=pa)
    s2, _, _ = _counted(rng, 300, 40.0, profile=pb)
    before = ad.concat([s1, s2])
    before.obs_names_make_unique()
    after, _, _ = _counted(rng, 600, 200.0, profile=(pa + pb) / 2)
    Eb, Ea = (np.asarray(x[:, panel].X.toarray(), np.float32) for x in (before, after))
    res = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, seed=3,
                                target_rms=1.0, n_anchor=150, k=4)
    X, rec = tb.depth_rendered_expression(before, after, res, panel, 0.5, seed=5, mode="hi")
    b = res.stage == 0
    rows = res.rows[b]
    first = rows < 300  # cells of the first state: genes 6-11 must stay almost empty after the top-up
    c = np.expm1(X[b].astype(np.float64))
    off = np.where(first, c[:, 6:12].sum(1), c[:, :6].sum(1)) / c.sum(1)
    assert rec["factor"]["before"] > 3 and off.mean() < 0.02
    with pytest.raises(ValueError, match="mode"):
        tb.depth_rendered_expression(before, after, res, panel, 0.5, seed=5, mode="deep")


def _labelled_pair(rng):
    """Two brackets whose cells are of two labelled kinds: an extraembryonic label both stages hold, and an embryonic
    one. Between the stages every cell's genes 0-3 double and genes 4-7 halve (a run's per-gene factor); the embryonic
    label also switches genes 8-11 on (development), which the series factor must not pick up."""
    import pandas as pd

    base = np.r_[np.ones(4), np.ones(4), np.full(4, 0.2), np.ones(3)]
    run = np.r_[np.full(4, 2.0), np.full(4, 0.5), np.ones(4), np.ones(3)]
    dev = np.r_[np.ones(8), np.full(4, 6.0), np.ones(3)]
    parts_b, parts_a = [], []
    for lab, grp, d in (("EXE-Endoderm", "Endoderm", np.ones(15)), ("SOM", "ParaxialMesoderm", dev)):
        b, panel, _ = _counted(rng, 300, 40.0, profile=base)
        a, _, _ = _counted(rng, 300, 200.0, profile=base * run * d)
        for x in (b, a):
            x.obs = pd.DataFrame({"celltype": [lab] * x.n_obs, "group": [grp] * x.n_obs}, index=x.obs_names)
        parts_b.append(b)
        parts_a.append(a)
    before, after = ad.concat(parts_b), ad.concat(parts_a)
    for x in (before, after):
        x.obs_names_make_unique()
    return before, after, panel


def test_series_rendering_moves_the_shallower_bracket_by_the_extraembryonic_factor():
    """Embryo try E24: the per-gene factor comes from the extraembryonic label both brackets hold, so it carries the
    run's factor and not the embryonic label's own development; the deeper bracket's cells are untouched."""
    rng = np.random.default_rng(61)
    before, after, panel = _labelled_pair(rng)
    assert tb.series_labels(before.obs["celltype"], before.obs["group"], after.obs["celltype"]) == ["EXE-Endoderm"]
    Eb, Ea = (np.asarray(x[:, panel].X.toarray(), np.float32) for x in (before, after))
    res = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, seed=3,
                                target_rms=1.0, n_anchor=150, k=4)
    X, rec = tb.series_rendered_expression(before, after, res, panel, 0.5, seed=[1, 7])
    b, a = res.stage == 0, res.stage == 1
    np.testing.assert_array_equal(X[a], Ea[res.rows[a]])  # the deeper bracket keeps its released rows
    assert rec["moved"] == "before" and rec["labels"] == ["EXE-Endoderm"] and rec["ratio"] > 3
    assert rec["rendered_median"] == pytest.approx(rec["median_depth"]["after"], rel=0.15)
    c = np.expm1(X[b].astype(np.float64))
    c0 = np.expm1(Eb[res.rows[b]].astype(np.float64))
    up = (c[:, :4].sum(1) / c.sum(1)).mean() / (c0[:, :4].sum(1) / c0.sum(1)).mean()
    down = (c[:, 4:8].sum(1) / c.sum(1)).mean() / (c0[:, 4:8].sum(1) / c0.sum(1)).mean()
    assert up > 1.2 and down < 0.8  # the run's factor is applied (shares of the doubled genes rise, of the halved fall)
    som = np.asarray(before.obs["celltype"])[res.rows[b]] == "SOM"
    dev_share = (c[som][:, 8:12].sum(1) / c[som].sum(1)).mean() / (c0[som][:, 8:12].sum(1) / c0[som].sum(1)).mean()
    assert dev_share < 1.5  # the embryonic label's sixfold switch is not carried over
    assert np.isfinite(X).all() and (X >= 0).all() and ((X > 0).sum(1) >= 1).all()
    X2, _ = tb.series_rendered_expression(before, after, res, panel, 0.5, seed=[1, 7])
    np.testing.assert_array_equal(X, X2)
    with pytest.raises(ValueError, match="no matched extraembryonic label"):
        tb.series_factor(sparse_counts(before, panel), before.obs["celltype"], sparse_counts(after, panel),
                         after.obs["celltype"], [])


def sparse_counts(a, panel):
    from scipy import sparse as sp

    return sp.csr_matrix(a[:, panel].layers["counts"], dtype=np.int64)


def test_thinning_never_adds_counts_and_keeps_integers():
    rng = np.random.default_rng(33)
    from scipy import sparse as sp

    C = sp.csr_matrix(rng.poisson(5.0, size=(200, 10)))
    other = rng.poisson(3.0, size=200)
    T, o = tb._binomial_rows(C, other, np.full(200, 0.4), rng)
    assert (T.toarray() <= C.toarray()).all() and (o <= other).all() and T.dtype.kind == "i"
    assert T.sum() == pytest.approx(0.4 * C.sum(), rel=0.05)
    same, o1 = tb._binomial_rows(C, other, np.ones(200), rng)
    assert (same.toarray() == C.toarray()).all() and (o1 == other).all()


def test_blend_splits_the_expected_counts_by_the_bracket_fraction():
    rng = np.random.default_rng(34)
    pb = np.r_[np.ones(6), np.full(6, 1e-9), np.ones(3)]  # the earlier stage expresses genes 0-5, the later 6-11
    pa = np.r_[np.full(6, 1e-9), np.ones(6), np.ones(3)]
    before, panel, _ = _counted(rng, 600, 40.0, profile=pb)
    after, _, _ = _counted(rng, 600, 200.0, profile=pa)
    Eb, Ea = (np.asarray(x[:, panel].X.toarray(), np.float32) for x in (before, after))
    res = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, seed=3,
                                target_rms=1.0, n_anchor=150, k=4, with_partners=True)
    assert res.partners is not None and len(res.partners) == 500
    assert (res.partners[res.stage == 0] < after.n_obs).all() and (res.partners[res.stage == 1] < before.n_obs).all()
    X, rec = tb.blend_expression(before, after, res, panel, 0.5, seed=5)
    c = np.expm1(X.astype(np.float64))  # proportional to each cell's counts (the library size cancels in a share)
    share = c[:, :6].sum(1) / c.sum(1)
    want = rec["expected_share_before"]
    assert rec["weights"]["before"] == 1.0 and 0.0 < rec["weights"]["after"] < 1.0  # the shallow side is capped
    assert share.mean() == pytest.approx(want, abs=0.05) and 0.3 < want < 0.6
    assert (share > 0).all() and (share < 1).all()  # every cell carries counts of both its cells
    no_partners = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, seed=3,
                                        target_rms=1.0, n_anchor=150, k=4)
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):  # partners come from their own stream
        np.testing.assert_array_equal(getattr(res, f), getattr(no_partners, f))
    with pytest.raises(ValueError, match="with_partners"):
        tb.blend_expression(before, after, no_partners, panel, 0.5, seed=5)


def test_extrapolated_depth_carries_half_the_latest_log_rate_over_the_horizon():
    # the heart series at test: E9.5 75, E10.5 56, two days ahead at B0's damping 0.5
    assert tb.extrapolated_depth(75.0, 56.0, 9.5, 10.5, 12.5, 0.5) == pytest.approx(56.0 * (56.0 / 75.0))
    assert tb.extrapolated_depth(130.0, 75.0, 8.75, 9.5, 10.5, 0.5) == pytest.approx(
        75.0 * np.exp(0.5 * np.log(75.0 / 130.0) / 0.75))
    assert tb.extrapolated_depth(50.0, 50.0, 1.0, 2.0, 4.0, 0.5) == pytest.approx(50.0)
    with pytest.raises(ValueError):
        tb.extrapolated_depth(75.0, 56.0, 10.5, 9.5, 12.5, 0.5)


def test_thin_to_depth_thins_by_the_ratio_and_leaves_a_deeper_target_alone():
    rng = np.random.default_rng(35)
    a, panel, N = _counted(rng, 400, 120.0)
    E = np.asarray(a[:, panel].X.toarray(), np.float32)
    same, rec = tb.thin_to_depth(a, panel, 120.0, 150.0, seed=1)
    np.testing.assert_array_equal(same, E)
    assert rec["factor"] == 1.0 and rec["thinned"] == 0
    X, rec = tb.thin_to_depth(a, panel, 120.0, 60.0, seed=1)
    assert rec["factor"] == pytest.approx(0.5) and rec["thinned"] == 400
    assert (X > 0).sum(1).mean() < (E > 0).sum(1).mean() and ((X > 0).sum(1) >= 1).all()
    assert tb.median_depth(a, panel) == pytest.approx(float(np.median(N)))



def test_interpolation_herd_moves_each_draw_mean_and_keeps_the_rest():
    rng = np.random.default_rng(36)
    before, after, panel, res, mix = _interp_pair(rng, 60.0, 60.0)
    Eb, Ea = (np.asarray(x[:, panel].X.toarray(), np.float32) for x in (before, after))
    kw = dict(seed=3, target_rms=1.0, n_anchor=150, k=4)
    h = tb.interpolate_clouds(before.obsm[COORD_KEY], after.obsm[COORD_KEY], Eb, Ea, 0.5, 500, herd=True, **kw)
    assert np.array_equal(h.stage, res.stage) and not np.array_equal(h.rows, res.rows)
    for s_, E in ((0, Eb), (1, Ea)):
        m0, m1 = E[res.rows[res.stage == s_]].mean(0), E[h.rows[h.stage == s_]].mean(0)
        assert ((m1 - E.mean(0)) ** 2).sum() < ((m0 - E.mean(0)) ** 2).sum()
    assert set(h.record["herd"]) == {"before", "after"}


# ------------------------------------------------------------------- the midpoint draw (embryo try E13), 2026-11-09
def _midpoint_pair(seed: int = 0, n: int = 800, genes: int = 30):
    """Two brackets of the same four states placed in four spatial blocks; each state's after cells are moved along
    its own expression direction (the local transport). Returns coordinates, expression, log depths and the states."""
    rng = np.random.default_rng(seed)
    states_b, states_a = rng.integers(0, 4, n), rng.integers(0, 4, n)
    blocks = np.array([[-8.0, 0.0, 0.0], [8.0, 0.0, 0.0], [0.0, -8.0, 0.0], [0.0, 8.0, 0.0]])
    A = blocks[states_b] + rng.normal(size=(n, 3))
    B = 1.1 * blocks[states_a] + rng.normal(size=(n, 3))
    base = 4.0 * rng.normal(size=(4, genes))
    shift = rng.normal(size=(4, genes))
    shift = 3.0 * shift / np.linalg.norm(shift, axis=1, keepdims=True)
    Eb = base[states_b] + 0.7 * rng.normal(size=(n, genes))
    Ea = base[states_a] + shift[states_a] + 0.7 * rng.normal(size=(n, genes))
    return A, B, Eb, Ea, rng.normal(3.0, 0.4, n), rng.normal(4.0, 0.4, n)


@pytest.fixture(scope="module")
def midpoint_case():
    A, B, Eb, Ea, Lb, La = _midpoint_pair()
    kw = dict(target_rms=1.0, n_anchor=300, k=4, seed=11)
    mid = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 600, draw="midpoint", log_depth=(Lb, La), with_coupling=True, **kw)
    return A, B, Eb, Ea, Lb, La, kw, mid, tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 600, **kw)


def test_midpoint_weights_move_each_stratum_half_way_to_zero_and_halve_its_variance():
    rng = np.random.default_rng(40)
    mu, sd = 3.0, 2.0
    s = np.r_[rng.normal(mu, sd, 400_000), rng.normal(-1.0, 0.5, 400_000)]  # two strata, each its own variance
    strata = np.r_[np.zeros(400_000, int), np.ones(400_000, int)]
    w = tb._midpoint_weights(s, strata)
    for st, (m0, v0) in {0: (mu, sd**2), 1: (-1.0, 0.25)}.items():
        x, p = s[strata == st], w[strata == st] / w[strata == st].sum()
        mean = float((p * x).sum())
        assert mean == pytest.approx(m0 / 2.0, abs=0.02)  # the Gaussian product rule: half way to the midpoint
        assert float((p * (x - mean) ** 2).sum()) == pytest.approx(v0 / 2.0, rel=0.02)
    # a stratum without spread stays flat, and weights never reach zero
    assert np.all(tb._midpoint_weights(np.full(5, 2.0), np.zeros(5, int)) == 1.0)
    far = tb._midpoint_weights(np.r_[0.0, 100.0], np.array([0, 0]))
    assert far[1] >= 1e-12 > 0.0


def test_draw_weighted_keeps_exact_stratum_quotas_and_repeats_only_past_the_pool():
    strata = np.r_[np.zeros(50, int), np.ones(30, int), np.full(20, 2)]
    w = np.random.default_rng(41).uniform(0.01, 1.0, 100)
    rows, expected = tb._draw_weighted(strata, w, 25, np.random.default_rng(1))
    assert np.bincount(strata[rows], minlength=3).tolist() == [13, 7, 5]  # largest remainders of 12.5, 7.5, 5
    assert len(set(rows.tolist())) == 25 and expected.sum() == pytest.approx(25)
    for st, k in ((0, 13), (1, 7), (2, 5)):
        assert expected[strata == st].sum() == pytest.approx(k)
    rows2, _ = tb._draw_weighted(strata, w, 25, np.random.default_rng(1))
    assert np.array_equal(rows, rows2)  # the stream alone decides the draw
    many, exp2 = tb._draw_weighted(strata, w, 140, np.random.default_rng(2))
    assert np.bincount(strata[many], minlength=3).tolist() == [70, 42, 28]  # the pool's shares, past each pool's size
    assert np.bincount(many, minlength=100).min() >= 1 and exp2.sum() == pytest.approx(140)
    # within a stratum, rows are drawn in proportion to the weights, the rest of the pool almost never first
    ws = np.ones(40)
    ws[:4] = 1000.0
    hits = np.zeros(40)
    for s_ in range(300):
        r, _ = tb._draw_weighted(np.zeros(40, int), ws, 4, np.random.default_rng([s_, 5]))
        hits[r] += 1
    assert hits[:4].min() >= 0.95 * 300 and hits[4:].max() < 0.05 * 300
    rows0, exp0 = tb._draw_weighted(strata, w, 0, np.random.default_rng(3))
    assert rows0.size == 0 and exp0.sum() == 0.0


def test_joint_pc_fit_reproduces_the_entrys_features_and_depth_is_regressed_out():
    rng = np.random.default_rng(42)
    Eb, Ea = rng.normal(size=(200, 15)), rng.normal(1.0, 2.0, size=(150, 15))
    Fb, Fa = tb._joint_pcs(Eb, Ea, 8)  # the entry's function, untouched by the midpoint draw (a pinned requirement)
    c, L = tb._joint_pc_loadings(Eb, Ea, 8)
    assert L.shape == (15, 8) and c.shape == (15,)
    np.testing.assert_allclose((Eb - c) @ L, Fb, atol=1e-8)
    np.testing.assert_allclose((Ea - c) @ L, Fa, atol=1e-8)
    # the fit is the one `_joint_pcs` makes, also when fewer components than requested exist
    X, y = rng.normal(size=(6, 15)), rng.normal(size=(2, 15))
    assert tb._joint_pc_loadings(X, y, 30)[1].shape[1] == tb._joint_pcs(X, y, 30)[0].shape[1]
    P = rng.normal(size=(500, 4)) + np.array([1.0, -2.0, 0.5, 3.0])
    depth = rng.normal(2.0, 0.6, 500)
    Q = tb._residualise_depth(P + depth[:, None] @ np.array([[2.0, -1.0, 0.0, 0.5]]), depth)
    np.testing.assert_allclose(Q.mean(0), P.mean(0) + depth.mean() * np.array([2.0, -1.0, 0.0, 0.5]), atol=1e-9)
    for j in range(4):
        assert abs(np.corrcoef(Q[:, j], depth)[0, 1]) < 1e-6
    np.testing.assert_array_equal(tb._residualise_depth(P, np.full(500, 2.0)), P)  # no depth variation: nothing moves


def test_midpoint_draw_keeps_the_nominal_mixture_and_leaves_the_default_path_identical(midpoint_case):
    A, B, Eb, Ea, Lb, La, kw, mid, base = midpoint_case
    # the default path: identical to an explicit uniform draw (log_depth, when given, ignored), no new record field
    uni = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 600, draw="uniform", log_depth=(Lb, La), **kw)
    for f in ("stage", "rows", "coords", "pool_coords", "pool_weights"):
        np.testing.assert_array_equal(getattr(uni, f), getattr(base, f))
    assert uni.record == base.record and "draw" not in base.record and "midpoint" not in base.record
    assert base.pool_weights[0] == pytest.approx(0.5 / 800)
    # exact bracket counts, the same shuffle of stages as the uniform draw, the same moves
    assert (mid.stage == 0).sum() == (base.stage == 0).sum() == 300 and (mid.stage == 1).sum() == 300
    np.testing.assert_array_equal(mid.stage, base.stage)

    def unit(C):  # the pool up to the similarity the draw sets (centre and scale follow the drawn cells)
        C = C - C.mean(0)
        return C / np.sqrt((C**2).sum(1).mean())

    np.testing.assert_allclose(unit(mid.pool_coords), unit(base.pool_coords), atol=1e-9)  # the same moved positions
    assert not np.array_equal(mid.rows, base.rows)  # only which rows are drawn changes
    assert mid.record["draw"] == "midpoint" and mid.record["midpoint"]["k_strata"] == tb.MIDPOINT_K_STRATA
    assert mid.record["midpoint"]["k_anchors"] == tb.MIDPOINT_K_ANCHORS
    m = mid.coupling["midpoint"]
    # every drawn row is a real cell of its bracket, at its own moved position (the pool's), the draw's mass kept
    Xb, Xa = mid.pool_coords[:800], mid.pool_coords[800:]
    rows_b, rows_a = mid.rows[mid.stage == 0], mid.rows[mid.stage == 1]
    assert rows_b.max() < 800 and rows_a.max() < 800
    assert len(set(rows_b.tolist())) == 300 and len(set(rows_a.tolist())) == 300  # no repeats below the pool's size
    for rows, X, sel in ((rows_b, Xb, mid.stage == 0), (rows_a, Xa, mid.stage == 1)):
        np.testing.assert_allclose(np.sort(mid.coords[sel], axis=0), np.sort(X[rows], axis=0), atol=1e-6)
    assert mid.pool_weights.sum() == pytest.approx(1.0)
    assert mid.pool_weights[:800].sum() == pytest.approx(0.5) and mid.pool_weights[800:].sum() == pytest.approx(0.5)
    # each bracket keeps each stratum's share of the pool, up to the rounding of one cell
    for rows, strata in ((rows_b, m["strata_before"]), (rows_a, m["strata_after"])):
        drawn = np.bincount(strata[rows], minlength=tb.MIDPOINT_K_STRATA) / len(rows)
        pool = np.bincount(strata, minlength=tb.MIDPOINT_K_STRATA) / len(strata)
        assert np.abs(drawn - pool).max() <= 1.0 / len(rows) + 1e-9
    # the uniform draw of the same seed is the one the entry would draw (for the gate's comparison)
    np.testing.assert_array_equal(np.sort(m["uniform_before"]), np.sort(base.rows[base.stage == 0]))
    np.testing.assert_array_equal(np.sort(m["uniform_after"]), np.sort(base.rows[base.stage == 1]))
    # the same call twice gives the same draw (deterministic k-means, fixed streams, none of them the herd's)
    again = tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 600, draw="midpoint", log_depth=(Lb, La), **kw)
    np.testing.assert_array_equal(again.rows, mid.rows)
    assert not set(tb.MIDPOINT_STREAMS) & {0, 1, 2, 3, 4, 5, *tb.HERD_STREAMS}


def test_midpoint_draw_tilts_each_bracket_towards_the_local_midpoint(midpoint_case):
    *_, mid, base = midpoint_case
    m = mid.coupling["midpoint"]
    s_b, s_a = m["s_before"], m["s_after"]
    # the anchor pairs' local midpoint sits between the brackets: before cells below it, after cells above it
    assert s_b.mean() < -0.5 and s_a.mean() > 0.5
    # the tilt moves each stratum's draw about half way to the midpoint (the expected draw shares, noise-free)
    moved = []
    for s, w, strata in ((s_b, m["expected_before"], m["strata_before"]),
                         (s_a, m["expected_after"], m["strata_after"])):
        for st in np.unique(strata):
            mask = strata == st
            if mask.sum() >= 40 and abs(s[mask].mean()) > s[mask].std():  # strata clearly off the midpoint
                moved.append((w[mask] / w[mask].sum() * s[mask]).sum() / s[mask].mean())
    assert len(moved) >= 4 and 0.35 < float(np.median(moved)) < 0.65
    # the draw itself follows: the drawn cells sit closer to the midpoint than the uniform draw's
    for s, rows, uni in ((s_b, mid.rows[mid.stage == 0], m["uniform_before"]),
                         (s_a, mid.rows[mid.stage == 1], m["uniform_after"])):
        assert abs(s[rows].mean()) < 0.8 * abs(s[uni].mean())
    # the two brackets stay on their own sides: a tilt, not an exchange of cells
    assert s_b[mid.rows[mid.stage == 0]].mean() < 0 < s_a[mid.rows[mid.stage == 1]].mean()
    # the effect is visible against two uniform draws of the same size (the gate's measure)
    assert m["w_before"].std() > 0.0 and m["w_after"].std() > 0.0
    ess = mid.record["midpoint"]["ess_frac"]
    assert 0.0 < ess["before"] < 1.0 and 0.0 < ess["after"] < 1.0
    assert base.record.get("midpoint") is None


def test_midpoint_draw_refuses_bad_arguments(midpoint_case):
    A, B, Eb, Ea, Lb, La, kw, *_ = midpoint_case
    with pytest.raises(ValueError, match="draw must be one of"):
        tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 100, draw="balanced", **kw)
    with pytest.raises(ValueError, match="needs log_depth"):
        tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 100, draw="midpoint", **kw)
    with pytest.raises(ValueError, match="log_depth"):
        tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 100, draw="midpoint", log_depth=(Lb[:-1], La), **kw)
    with pytest.raises(ValueError, match="both choose"):
        tb.interpolate_clouds(A, B, Eb, Ea, 0.5, 100, draw="midpoint", log_depth=(Lb, La), herd=True, **kw)


def test_log_panel_counts_sums_the_panels_counts_whatever_the_sparsity():
    rng = np.random.default_rng(43)
    a, panel, _ = _counted(rng, 50, 60.0)
    want = np.log1p(np.asarray(a[:, panel].layers["counts"].sum(1)).ravel())
    np.testing.assert_allclose(tb.log_panel_counts(a, panel), want)
    from scipy import sparse as sp

    dense_counts = ad.AnnData(X=a.X.copy(), var=a.var, layers={"counts": sp.csr_matrix(a.layers["counts"]).toarray()})
    np.testing.assert_allclose(tb.log_panel_counts(dense_counts, panel), want)
