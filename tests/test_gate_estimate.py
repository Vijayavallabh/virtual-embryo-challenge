"""experiments/gate/estimate.py: our estimated leaderboard total from committed results (aggregation, not a try).

The procedure on the committed files (its self-checks, reading A at least reading B on every board, K2a's place), on
small synthetic inputs, missing or inconsistent inputs failing loudly, the one door refusing board files, and no
leaderboard figure in the module."""

from __future__ import annotations

import ast
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "experiments" / "gate" / "estimate.py"

# Every number the module may hold, each a procedure or format constant with its source beside it in the module: the
# published floor (50), our own committed outlook of decision log 2026-11-04 (c) (63, 52, 60, 65, 184, 170-198), the
# band's percentiles, weights, tolerances, the draws and the RNG seed. A leaderboard figure is never added here, and the
# captain's gate threshold is read from results/upload_gate.json, so it never appears in the module.
ALLOWED_NUMBERS = {0, 1e-9, 0.25, 0.5, 1, 2, 10, 50, 52, 60, 63, 65, 90, 100, 118, 170, 184, 198, 200_000, 20261105}
ALLOWED_TOKENS = {"0", "0.0", "0.25", "0.5", "0.60", "1", "1.0", "2", "8", "8.2", "9", "9.2", "10", "50", "50.0", "52",
                  "60", "63", "65", "90", "100", "118", "170", "184", "198", "20261105"}
# a number standing alone in the text: not a file:line reference, a range's end, a date, a stage, a label or a path
STANDALONE = re.compile(r"(?<![\w.:\-/])\d+(?:\.\d+)?(?![\w.\-]*\d)")


@pytest.fixture(scope="module")
def est():
    spec = importlib.util.spec_from_file_location("gate_estimate", MODULE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def res(est):
    return est.estimate()


# ---------------------------------------------------------------------------------------------- committed files
def test_readings_on_the_committed_files(est, res):
    r = res["readings"]
    for name in ("A", "A_as_of_c", "B", "B_as_of_c"):
        d = r[name]
        assert d["Overall"] == pytest.approx(d["T1"] + (d["heart"] + d["embryo"]) / 2 + d["T3"])
        for b in (*est.BOARDS, "Overall"):
            lo, hi = d["band"][b]
            assert lo <= hi
        lo, hi = d["band"]["Overall"]
        assert lo < d["Overall"] < hi
    # a board shows the better of its two uploads, so A never reads below the first round, board by board
    for name, first in (("A", "B"), ("A_as_of_c", "B_as_of_c")):
        for b in (*est.BOARDS, "Overall"):
            assert r[name][b] >= r[first][b]
    for b in (*est.BOARDS, "Overall"):  # A and B read the same drawn scenarios, so each percentile keeps the order
        assert all(x >= y for x, y in zip(r["A"]["band"][b], r["B"]["band"][b], strict=True))
    # K2a moves T3 only, and it points down; reading B moves only by K2a's fold standing in for a quarter of s2's
    assert r["A"]["T3"] < r["A_as_of_c"]["T3"]
    assert all(r["A"][b] == r["A_as_of_c"][b] and r["B"][b] == r["B_as_of_c"][b] for b in ("T1", "heart", "embryo"))
    assert abs(r["B"]["T3"] - r["B_as_of_c"]["T3"]) < 0.1
    assert res["first_round"] == {"T1": "s1", "heart": "s1", "embryo": "s1", "T3": "s2"}
    # every self-check ran, and every input came through the door from committed results
    assert len(res["checks"]) == 20
    assert set(res["inputs_sha256"]) >= {est.T1_GRID_PAYOFF, est.T1_GRID_HEDGE, est.T1_DIAG, est.T1_D2,
                                         est.T3_RECORDS, est.T3_SUMMARY, est.T3_K2, est.GATE}
    assert all(k.startswith(est.READABLE) for k in res["inputs_sha256"])
    assert not any(k.startswith(est.REFUSED) or "index.json" in k for k in res["inputs_sha256"])


def test_committed_anchors(est, res):
    """The readings rest on the committed summaries the brief names: T1's grid (the pair 68.28 and s1 alone 63.54 under
    the uniform prior), the embryo's s2 rule (60.87 and fgw alone 59.95) and T3's stand-ins (TRIES.md)."""
    hedge = json.loads((ROOT / est.T1_GRID_HEDGE).read_text())
    t1 = res["steps"]["T1"]
    assert t1["grid_A"] == pytest.approx(hedge["candidates"][est.T1_S2]["uniform"]["mean"], abs=1e-9)
    assert t1["grid_B"] == pytest.approx(hedge["candidates"][est.T1_S1]["uniform"]["mean"], abs=1e-9)
    assert (round(t1["grid_A"], 2), round(t1["grid_B"], 2)) == (68.28, 63.54)
    r = res["readings"]
    assert r["A"]["T1"] == pytest.approx(t1["grid_A"] + t1["m1_minus_b0"])
    assert r["B"]["T1"] == pytest.approx(t1["grid_B"] + t1["m1_minus_b0"])
    assert (round(r["A"]["embryo"], 2), round(r["B"]["embryo"], 2)) == (60.87, 59.95)
    means = res["steps"]["T3_stand_in_means"]
    assert [round(means[k], 2) for k in ("A", "B", "s2")] == [73.33, 60.77, 50.87]
    assert r["B_as_of_c"]["T3"] == pytest.approx(means["s2"])
    assert res["steps"]["prior"] == "uniform"


def test_outlook_comparison(est, res):
    assert res["outlook"] == {"T1": 63, "heart": 52, "embryo": 60, "T3": 65, "Overall": 184, "range": [170, 198],
                              "source": "docs/PLAN.md decision log 2026-11-04 (c)"}
    a = res["readings"]["A_as_of_c"]
    for b in est.BOARDS:
        assert res["A_as_of_c_minus_outlook"][b] == pytest.approx(a[b] - est.OUTLOOK[b])
    assert res["reproduces_outlook"] == all(abs(res["A_as_of_c_minus_outlook"][b]) <= 1 for b in est.BOARDS)


def test_k2a_is_split_symmetrically(est):
    c, k = est.T3_WEIGHTS_AS_OF_C, est.T3_WEIGHTS_K2A
    assert sum(c.values()) == sum(k.values()) == 1
    assert k["A"] == c["A"] and k["B"] == k["K2a"] == c["B"] / 2


def test_cli_prints_and_writes_json(est, res, tmp_path, capsys):
    out = tmp_path / "gate" / "estimate.json"
    assert est.main(["--json", str(out)]) == 0
    printed = capsys.readouterr().out
    assert "not a try" in printed and "reproduces the outlook" in printed
    assert all(len(line) <= 120 for line in printed.splitlines())
    d = json.loads(out.read_text())
    out.unlink()
    assert d["readings"] == json.loads(json.dumps(res["readings"]))  # the fixed seed reproduces the bands too
    assert d["gate"]["threshold"] == json.loads((ROOT / est.GATE).read_text())["threshold"]


# ------------------------------------------------------------------------------------------------ synthetic inputs
def test_grid_outcomes_synthetic(est):
    rng = np.random.default_rng(0)
    s1 = {"t1": rng.normal(60, 2, 5), "t2": rng.normal(50, 2, 5)}
    s2 = {"t1": rng.normal(55, 2, 5), "t2": rng.normal(58, 2, 5)}
    o = est.grid_outcomes(s1, s2, {"t1": 1.0, "t2": 3.0})  # an unnormalised prior is normalised
    assert o["w"].sum() == pytest.approx(1.0)
    assert np.all(o["a"] >= o["b"])
    assert est.reading(o, "b") == pytest.approx(0.25 * s1["t1"].mean() + 0.75 * s1["t2"].mean())
    better = {t: np.maximum(s1[t], s2[t]).mean() for t in s1}
    assert est.reading(o, "a") == pytest.approx(0.25 * better["t1"] + 0.75 * better["t2"])
    twin = est.grid_outcomes(s1, s1, {"t1": 1.0, "t2": 1.0}, shift=1.5)  # an s2 equal to s1 adds nothing
    assert est.reading(twin, "a") == pytest.approx(est.reading(twin, "b"))
    assert est.reading(twin, "b") == pytest.approx((s1["t1"].mean() + s1["t2"].mean()) / 2 + 1.5)
    with pytest.raises(ValueError):
        est.grid_outcomes({"t": np.ones(3)}, {"t": np.ones(4)}, {"t": 1.0})


def test_t3_outcomes_synthetic(est):
    s2 = {"rep1": np.array([50.0, 52.0]), "rep2": np.array([48.0, 49.0])}
    standins = {"A": {"rep1": np.array([70.0, 71.0]), "rep2": np.array([72.0, 73.0])},
                "B": {"rep1": np.array([51.0, 51.0]), "rep2": np.array([60.0, 61.0])},
                "K2a": (np.array([49.0, 53.0]), np.array([50.0, 52.0]))}
    c = est.t3_outcomes(standins, s2, {"A": 0.5, "B": 0.5})
    k = est.t3_outcomes(standins, s2, {"A": 0.5, "B": 0.25, "K2a": 0.25})
    # reading B is s2 alone in every scenario: on the replicates s2's stand-in, on K2a's fold our WT copy there
    assert est.reading(c, "b") == pytest.approx(np.mean([*s2["rep1"], *s2["rep2"]]))
    assert est.reading(k, "b") == pytest.approx(0.75 * est.reading(c, "b") + 0.25 * np.mean([50.0, 52.0]))
    # a stand-in on the replicates is paired with s2's by seed: on rep1, B's second seed loses to s2's 52
    a_a = (np.mean([70, 71]) + np.mean([72, 73])) / 2
    a_b = (np.mean([51, 52]) + np.mean([60, 61])) / 2
    assert est.reading(c, "a") == pytest.approx(0.5 * a_a + 0.5 * a_b)
    # K2a, on another fold, is paired by seed with the s2 stand-in of that fold (our WT copy there): the same truth
    a_k = np.mean([max(49, 50), max(53, 52)])
    assert est.reading(k, "a") == pytest.approx(0.5 * a_a + 0.25 * a_b + 0.25 * a_k)
    with pytest.raises(TypeError):
        est.t3_outcomes({"K2a": np.array([49.0, 53.0])}, s2, {"K2a": 1.0})
    for o in (c, k):
        assert o["w"].sum() == pytest.approx(1.0) and np.all(o["a"] >= o["b"])


def test_mixture_and_bands_synthetic(est):
    def point(v):
        return {"a": np.array([v + 1.0]), "b": np.array([v]), "w": np.array([1.0])}

    boards = {"T1": point(60.0), "heart": point(50.0), "embryo": point(60.0), "T3": point(50.0)}
    got = est.bands(boards, n=1000, seed=1)
    assert got["A"]["Overall"] == pytest.approx([61 + (51 + 61) / 2 + 51] * 2)
    assert got["B"]["Overall"] == pytest.approx([60 + (50 + 60) / 2 + 50] * 2)
    assert est.readings(boards)["A"]["Overall"] == pytest.approx(168.0)
    m = est.mixture([(1.0, point(40.0)), (3.0, point(80.0))])
    assert list(m["w"]) == pytest.approx([0.25, 0.75]) and est.reading(m, "b") == pytest.approx(70.0)
    spread = {**boards, "T1": m}
    assert est.bands(spread, n=5000, seed=7) == est.bands(spread, n=5000, seed=7)  # a fixed seed reproduces
    lo, hi = est.bands(spread, n=5000, seed=7)["B"]["T1"]
    assert (lo, hi) == (40.0, 80.0)


# ----------------------------------------------------------------------------------------- failing loudly, the door
@pytest.fixture
def copy(est, tmp_path):
    """A copy of exactly the inputs `required()` lists, removed after the test."""
    dst = tmp_path / "copy"
    for rel in est.required():
        if rel.endswith("/"):
            shutil.copytree(ROOT / rel, dst / rel)
        else:
            (dst / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / rel, dst / rel)
    yield dst
    shutil.rmtree(dst, ignore_errors=True)


def test_missing_inputs_fail_loudly(est, tmp_path, copy):
    assert est.missing_inputs(ROOT) == []
    with pytest.raises(SystemExit, match="missing committed inputs") as err:
        est.estimate(root=tmp_path / "empty")
    for rel in (est.T1_GRID_HEDGE, est.T1_DIAG, est.T3_RECORDS, est.T3_K2, est.GATE, est.T1_GRID_RUNS):
        assert rel in str(err.value)
    assert est.estimate(root=copy, n_draws=1000)["readings"]["B"]["T1"] == pytest.approx(
        est.estimate(root=ROOT, n_draws=1000)["readings"]["B"]["T1"])  # the listed inputs are all it needs
    (copy / est.T3_K2).unlink()
    assert est.missing_inputs(copy) == [est.T3_K2]
    with pytest.raises(SystemExit, match=re.escape(est.T3_K2)):
        est.estimate(root=copy)


def test_a_missing_part_or_field_names_its_file(est, copy):
    parts = sorted(p for p in (copy / est.HEART[0][0]).iterdir() if p.is_dir())
    shutil.rmtree(parts[-1])  # one per-seed part of the heart payoff: the folder still holds others
    with pytest.raises(SystemExit, match=re.escape(est.HEART[0][1])):
        est.estimate(root=copy, n_draws=1000)
    shutil.copytree(ROOT / est.HEART[0][0] / parts[-1].name, parts[-1])
    path = copy / est.T1_DIAG
    d = json.loads(path.read_text())
    for r in d["rows"]:
        r["scores"].pop(est.T1_M1)
    path.write_text(json.dumps(d))
    with pytest.raises(SystemExit, match=re.escape(est.T1_DIAG) + ".*missing or malformed field"):
        est.estimate(root=copy, n_draws=1000)


def test_weighted_percentile_is_exact(est):
    x, w = np.array([2.0, 1.0, 4.0, 3.0]), np.full(4, 0.25)
    assert [est.weighted_percentile(x, w, q) for q in (10, 25, 90)] == [1.0, 1.5, 4.0]
    assert est.weighted_percentile(np.array([1.0, 2.0]), np.array([0.1, 0.9]), 10) == 1.5
    assert est.weighted_percentile(np.array([1.0, 2.0]), np.array([0.1, 0.9]), 90) == 2.0


def test_inconsistent_inputs_fail_loudly(est, copy):
    path = copy / est.T1_GRID_HEDGE
    hedge = json.loads(path.read_text())
    hedge["candidates"][est.T1_S2]["uniform"]["mean"] += 0.01
    path.write_text(json.dumps(hedge))
    with pytest.raises(SystemExit, match="check failed: T1 grid"):
        est.estimate(root=copy)


@pytest.mark.parametrize("rel", ["boards/index.json", "panels/index.json", "docs/research/digest_2026-10-19.md",
                                 "data/raw/T1/E10.5.h5ad", "runs/harness/d2_v1/scores.json",
                                 "experiments/harness/results/../../../boards/index.json",
                                 "experiments/harness/results/index.json", "docs/PLAN.md", "/etc/hostname"])
def test_door_refuses_before_opening(est, rel):
    door = est.Door(ROOT)
    with pytest.raises(SystemExit, match="refusing"):
        door.path(rel)
    assert door.read == {}


def test_door_refuses_a_link_out_of_the_readable_folders(est, tmp_path):
    (tmp_path / "boards").mkdir()
    (tmp_path / "boards" / "x.json").write_text("{}")
    inside = tmp_path / "experiments" / "harness" / "results"
    inside.mkdir(parents=True)
    (inside / "x.json").symlink_to(tmp_path / "boards" / "x.json")
    door = est.Door(tmp_path)
    with pytest.raises(SystemExit, match="refusing .* it resolves to boards/x.json"):
        door.path("experiments/harness/results/x.json")
    assert door.read == {}


# ---------------------------------------------------------------------------------------- no leaderboard figure
def test_no_leaderboard_figure_in_the_module(est):
    text = MODULE.read_text()
    tree = ast.parse(text)
    numbers = {n.value for n in ast.walk(tree)
               if isinstance(n, ast.Constant) and isinstance(n.value, int | float) and not isinstance(n.value, bool)}
    assert numbers <= ALLOWED_NUMBERS, sorted(numbers - ALLOWED_NUMBERS)
    tokens = set(STANDALONE.findall(text))
    assert tokens <= ALLOWED_TOKENS, sorted(tokens - ALLOWED_TOKENS)
    threshold = json.loads((ROOT / est.GATE).read_text())["threshold"]
    assert threshold not in numbers and str(threshold) not in tokens  # read from the record, never written here
    assert "boards/index.json" not in text and "panels/index.json" not in text
    assert not any(p.startswith(("boards", "panels", "docs/research", "data", "runs")) for p in est.READABLE)
