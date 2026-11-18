"""experiments/gate/second_file_value.py: its rules on hand-made tables and synthetic draws, its reader's refusals, and
its committed output against the files it read. Nothing here reads a competition file or a target."""

import hashlib
import importlib.util
import json
import random
import statistics
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "experiments" / "gate" / "second_file_value.py"
RES = REPO / "experiments" / "gate" / "results"
OUT_JSON = RES / "second_file_value_2026-11-18.json"
OUT_TXT = RES / "second_file_value_2026-11-18.txt"


def _load():
    spec = importlib.util.spec_from_file_location("second_file_value", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load()


# ------------------------------------------------------------------------------------------------- the arithmetic
def test_paired_is_candidate_minus_first_seed_by_seed():
    assert S.paired([50.0, 51.0, 52.0], [51.0, 51.0, 50.0]) == [1.0, 0.0, -2.0]
    with pytest.raises(AssertionError):
        S.paired([1.0, 2.0], [1.0])


def test_seed_summary():
    s = S.seed_summary([1.0, -1.0, 3.0])
    assert s["mean"] == pytest.approx(1.0) and s["ahead"] == 2 and s["seeds"] == 3
    assert s["seed_sd"] == pytest.approx(statistics.stdev([1.0, -1.0, 3.0]))


def test_tasks_sharing_a_truth_are_one_unit():
    units = S.truth_units({"a": "T/E8.0", "b": "T/E8.0", "c": "T/E7.5"})
    assert units == {"T/E8.0": ["a", "b"], "T/E7.5": ["c"]}
    v = S.unit_values({"a": [2.0, 2.0], "b": [-2.0, -2.0], "c": [1.0, 3.0]}, units)
    assert v["T/E8.0"]["mean_d"] == pytest.approx(0.0)  # the mean of its tasks' means
    assert v["T/E8.0"]["seed_v"] == pytest.approx(1.0)  # mean over tasks and seeds of max(D, 0): (2 + 2 + 0 + 0) / 4
    assert v["T/E7.5"]["mean_d"] == pytest.approx(2.0) and v["T/E7.5"]["seed_v"] == pytest.approx(2.0)


def test_merging_changes_the_value():
    """Two tasks on one truth that disagree count once: merged, their gain and loss cancel before max(., 0)."""
    d = {"a": [2.0] * 10, "b": [-2.0] * 10}
    merged = S.summarise(S.unit_values(d, S.truth_units({"a": "T", "b": "T"})))
    apart = S.summarise(S.unit_values(d, S.truth_units({"a": "T1", "b": "T2"})))
    assert merged["v_task"] == pytest.approx(0.0) and merged["share_positive"] == 0.0
    assert apart["v_task"] == pytest.approx(1.0) and apart["share_positive"] == pytest.approx(0.5)


def test_summarise_equal_weights():
    units = {"u1": {"mean_d": 2.0, "seed_v": 2.5}, "u2": {"mean_d": -1.0, "seed_v": 0.5},
             "u3": {"mean_d": 0.5, "seed_v": 1.0}}
    s = S.summarise(units)
    assert s["mean_d"] == pytest.approx(0.5) and s["sd_d"] == pytest.approx(statistics.stdev([2.0, -1.0, 0.5]))
    assert s["v_task"] == pytest.approx((2.0 + 0.0 + 0.5) / 3) and s["v_seed"] == pytest.approx(4.0 / 3)
    assert s["share_positive"] == pytest.approx(2 / 3) and s["sd_kind"] == "n-1" and s["units"] == 3


def test_summarise_under_a_prior():
    units = {"u1": {"mean_d": 2.0, "seed_v": 2.5}, "u2": {"mean_d": -1.0, "seed_v": 0.5}}
    s = S.summarise(units, {"u1": 1.0, "u2": 3.0})  # normalised to 0.25 / 0.75
    assert s["mean_d"] == pytest.approx(0.25 * 2 - 0.75) and s["v_task"] == pytest.approx(0.5)
    assert s["v_seed"] == pytest.approx(0.25 * 2.5 + 0.75 * 0.5) and s["share_positive"] == pytest.approx(0.25)
    assert s["sd_d"] == pytest.approx((0.25 * (2 + 0.25) ** 2 + 0.75 * (-1 + 0.25) ** 2) ** 0.5)
    one = S.summarise({"u1": {"mean_d": 2.0, "seed_v": 2.5}})
    assert one["sd_d"] is None and one["v_task"] == pytest.approx(2.0)  # one unit: no spread
    zero = S.summarise(units, {"u1": 0.0, "u2": 1.0})  # a zero-weight unit counts for nothing
    assert zero["v_task"] == 0.0 and zero["mean_d"] == pytest.approx(-1.0)


def test_identical_candidate_is_worth_nothing_and_a_different_one_is_not():
    """The control (the first file's analogue) against itself gives D = 0 and V = 0; the same table with a candidate
    that differs gives a different value, so the comparison can tell them apart."""
    rng = random.Random(1)
    first = {f"t{i}": [50 + rng.gauss(0, 2) for _ in range(10)] for i in range(4)}
    targets = {t: t for t in first}
    same = S.candidate_table(first, first, targets)["summary"]
    assert same["v_task"] == 0.0 and same["v_seed"] == 0.0 and same["mean_d"] == 0.0 and same["share_positive"] == 0.0
    better = {t: [x + 1.0 for x in v] for t, v in first.items()}
    up = S.candidate_table(first, better, targets)["summary"]
    assert up["mean_d"] == pytest.approx(1.0) and up["v_task"] == pytest.approx(1.0) and up["share_positive"] == 1.0
    worse = {t: [x - 1.0 for x in v] for t, v in first.items()}
    assert S.candidate_table(first, worse, targets)["summary"]["v_task"] == 0.0


def test_pure_noise_gives_a_positive_value():
    """On pure noise (D ~ N(0, 1) per seed, no gain anywhere) V is not zero: at seed level it is E[max(Z, 0)] = 0.40,
    at task level about 0.40 / sqrt(10) = 0.13 per ten-seed task. V is an option value of the better of two, so a
    positive V is no evidence of a gain; the mean D and the share of units ahead are what show one."""
    rng = random.Random(20261118)
    n_tasks = 400
    d = {f"t{i}": [rng.gauss(0, 1) for _ in range(10)] for i in range(n_tasks)}
    s = S.summarise(S.unit_values(d, {t: [t] for t in d}))
    assert s["v_seed"] == pytest.approx(0.3989, abs=0.03)
    assert s["v_task"] == pytest.approx(0.3989 / 10 ** 0.5, abs=0.03)
    assert abs(s["mean_d"]) < 0.06 and s["share_positive"] == pytest.approx(0.5, abs=0.07)


# --------------------------------------------------------------------------------------- constants of a module
def test_module_constants_reads_arithmetic_literals_without_importing():
    src = ("import numpy as np\nDEPTHS = (0.35, 1.37)\nDEPTH_PRIORS = {'census': (1 / 7, 6 / 7), 'flat': (0.5,) * 2}\n"
           "MAIN = ('census', 'flat')\nOTHER = np.zeros(3)\n")
    c = S.module_constants(src, ("DEPTHS", "DEPTH_PRIORS", "MAIN"))
    assert c["DEPTHS"] == (0.35, 1.37) and c["MAIN"] == ("census", "flat")
    assert c["DEPTH_PRIORS"]["census"] == pytest.approx((1 / 7, 6 / 7)) and c["DEPTH_PRIORS"]["flat"] == (0.5, 0.5)
    with pytest.raises(ValueError):
        S.module_constants(src, ("OTHER",))  # a call is not a constant
    with pytest.raises(ValueError):
        S.module_constants(src, ("MISSING",))


def test_t30_priors_are_the_committed_ones():
    c = S.module_constants((REPO / S.T30_SRC).read_text(), ("DEPTH_PRIORS", "CM_PRIORS", "MAIN", "CHECKS"))
    assert set(c["DEPTH_PRIORS"]) == {"census", "latest", "flat"} and set(c["CM_PRIORS"]) == {"flat", "like", "rich"}
    for table in (c["DEPTH_PRIORS"], c["CM_PRIORS"]):
        assert all(sum(w) == pytest.approx(1.0) for w in table.values())
    assert c["MAIN"] == ("census", "flat")


# ---------------------------------------------------------------------------------------------------- the reader
@pytest.mark.parametrize("rel", ["boards/index.json", "data/raw/x.h5ad", "runs/builds/T1_s1.h5ad",
                                 "docs/research/digest_2026-10-19.md", "docs/PLAN.md",
                                 "experiments/x/notes/a.json", "experiments/../boards/index.json",
                                 "experiments/harness/results/index.json", "/etc/passwd", ".claude/x.json"])
def test_reader_refuses(rel, tmp_path):
    with pytest.raises(PermissionError):
        S.Reader(tmp_path).bytes(rel)


def test_reader_records_each_file(tmp_path):
    p = tmp_path / "experiments" / "a" / "results" / "x.json"
    p.parent.mkdir(parents=True)
    p.write_text('{"v": 1}')
    rd = S.Reader(tmp_path)
    assert rd.json("experiments/a/results/x.json") == {"v": 1}
    assert rd.read == {"experiments/a/results/x.json": hashlib.sha256(b'{"v": 1}').hexdigest()}


def _task(tmp_path, rows, clean=True):
    p = tmp_path / "experiments" / "t" / "results" / "task.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"task": {"key": "K", "target": "T/E8.0.h5ad", "reference": "T/E7.5.h5ad"},
                             "build": "main", "provenance": {"tree_clean": clean}, "rows": rows}))
    return "experiments/t/results/task.json"


def test_task_scores_needs_ten_seeds_and_a_clean_tree(tmp_path):
    rows = [{"seed": s, "scores": {"floor": 50.0, "a": 50.0 + s}} for s in reversed(range(10))]
    t = S.task_scores(S.Reader(tmp_path), _task(tmp_path, rows), [])
    assert t["arms"]["a"] == [50.0 + s for s in range(10)] and t["target"] == "T/E8.0.h5ad"
    with pytest.raises(AssertionError):
        S.task_scores(S.Reader(tmp_path), _task(tmp_path, rows[:9]), [])
    with pytest.raises(AssertionError):
        S.task_scores(S.Reader(tmp_path), _task(tmp_path, rows, clean=False), [])


# ------------------------------------------------------------------------------------------- the committed output
@pytest.fixture(scope="module")
def out():
    return json.loads(OUT_JSON.read_text())


def test_inputs_are_committed_result_files(out):
    assert out["inputs_sha256"]
    for rel, sha in out["inputs_sha256"].items():
        assert rel.startswith(("experiments/", "results/registry.csv")), rel
        assert not any(p in rel.split("/") for p in ("boards", "data", "runs", "notes", "research")), rel
        assert hashlib.sha256((REPO / rel).read_bytes()).hexdigest() == sha, f"{rel} changed since the run"


def _same(x, y, where=""):
    """Equal structure and strings; floats to 1e-9 (statistics' last digits differ between Python versions)."""
    assert type(x) is type(y) or {type(x), type(y)} <= {int, float}, where
    if isinstance(x, dict):
        assert x.keys() == y.keys(), where
        for k in x:
            _same(x[k], y[k], f"{where}/{k}")
    elif isinstance(x, list):
        assert len(x) == len(y), where
        for i, (u, v) in enumerate(zip(x, y, strict=True)):
            _same(u, v, f"{where}[{i}]")
    elif isinstance(x, float):
        assert x == pytest.approx(y, abs=1e-9), where
    else:
        assert x == y, where


def test_rerun_gives_the_committed_output(out):
    res = S.build()
    _same(json.loads(json.dumps(res)), out)
    assert S.text(res) == OUT_TXT.read_text()


def test_output_states_the_structure_asked_for(out):
    e = out["embryo"]
    assert e["first"] == "fgw" and set(e["candidates"]) == set(S.EMB_CANDS)
    assert {u: sorted(v) for u, v in e["units"].items()} == {
        "T2/E8.0.h5ad": ["half", "late"], "T2/E8.25_late.h5ad": ["new1", "x1"], "T2/E8.5.h5ad": ["new2", "x2", "x3"],
        "T2/E7.5.h5ad": ["e75"], "T2/E7.25.h5ad": ["e725"]}
    t1 = out["T1"]
    assert t1["first"] == "lib3_base" and len(t1["priors"]) == 9 and len(t1["scenarios"]) == 12
    assert all(set(v["by_prior"]) == set(t1["priors"]) for v in t1["candidates"].values())
    h = out["heart"]
    assert h["first"] == "b0_trend" and h["candidate"] == "copy" and len(h["units"]) == 1
    assert h["s2_hedge_heart"]["holds_first_or_candidate"] == []
    t3 = out["T3"]
    assert t3["paired_analogue"] is None and "no seed-paired D exists" in t3["why_none"]
    for v in t3["candidates"].values():  # T3 is against the floor: no value is computed
        for tab in (v, v["by_knockout"]):
            assert "v_task" not in tab["summary"] and "v_seed" not in tab["summary"]
            assert all("seed_v" not in u for u in tab["units"].values())
    for board in ("embryo", "T1", "heart", "T3"):
        assert out[board]["limits"], f"{board}: its limitations are in the output"
