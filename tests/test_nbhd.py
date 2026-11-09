"""N15 and E15 (experiments/harness/TRIES.md), on synthetic data only: the local within-type re-arrangement
(vec.t2_build.local_type_permutation), the audit's ratios and pass rule (experiments/t2_nbhd/n15_nbhd_audit.py), its
outputs (ratios only, never a cell count, and they refuse to overwrite), and the embryo candidate fgw_perm, which moves
only positions among the drawn cells so that only the neighbourhood metric of the scorer moves."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import GENES, simulate_stage
from scipy.spatial import cKDTree

from vec import harness, t2_build
from vec.boards import Board
from vec.evaluate import Stage, panel_metrics
from vec.harness import CANDIDATES, PseudoTask
from vec.io import COORD_KEY

REPO = Path(__file__).resolve().parents[1]
K = 15


@pytest.fixture(scope="module")
def n15():
    path = REPO / "experiments" / "t2_nbhd" / "n15_nbhd_audit.py"
    spec = importlib.util.spec_from_file_location("n15_nbhd_audit", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------------- the local within-type chain
def _cloud(n, rng):
    """A rectangular 3-D cloud with one random label per cell."""
    return rng.uniform(size=(n, 3)) * np.array([1.0, 1.0, 0.4]), rng.choice(["a", "b", "c"], size=n)


def test_the_permutation_keeps_every_row_in_its_label_and_within_the_radius():
    rng = np.random.default_rng(3)
    C, labels = _cloud(1200, rng)
    labels[7] = "singleton"
    radius = 2.0 * t2_build.neighbourhood_radius(C, K)
    at = t2_build.local_type_permutation(C, labels, radius, 30, np.random.default_rng(1))
    assert np.array_equal(np.sort(at), np.arange(len(C)))  # a permutation of the rows
    assert np.array_equal(labels[at], labels)  # each row's new position is a position of its label
    assert np.linalg.norm(C[at] - C, axis=1).max() <= radius  # each row ends within the radius of its start
    assert (at != np.arange(len(C))).mean() > 0.5  # it does move most rows
    assert at[7] == 7  # the only cell of its label has no partner and stays
    # reproducible from the generator alone; another generator gives another arrangement
    again = t2_build.local_type_permutation(C, labels, radius, 30, np.random.default_rng(1))
    other = t2_build.local_type_permutation(C, labels, radius, 30, np.random.default_rng(2))
    assert np.array_equal(at, again) and not np.array_equal(at, other)
    # no rounds, or a radius too small for any partner, leave every row in place
    assert np.array_equal(t2_build.local_type_permutation(C, labels, radius, 0, np.random.default_rng(1)),
                          np.arange(len(C)))
    assert np.array_equal(t2_build.local_type_permutation(C, labels, 1e-9, 30, np.random.default_rng(1)),
                          np.arange(len(C)))


def test_the_permutation_refuses_bad_arguments():
    rng = np.random.default_rng(4)
    C, labels = _cloud(30, rng)
    for radius in (0.0, -1.0, np.inf, np.nan):
        with pytest.raises(ValueError, match="radius"):
            t2_build.local_type_permutation(C, labels, radius, 3, np.random.default_rng(0))
    with pytest.raises(ValueError, match="sweeps"):
        t2_build.local_type_permutation(C, labels, 0.1, -1, np.random.default_rng(0))
    with pytest.raises(ValueError, match="one entry per cell"):
        t2_build.local_type_permutation(C, labels[:-1], 0.1, 3, np.random.default_rng(0))
    with pytest.raises(ValueError, match="fewer than k"):
        t2_build.neighbourhood_radius(C[:5], K)


def test_the_partner_lists_do_not_depend_on_the_tree():
    """The chain's only source of randomness is the generator: its partner lists are sorted, not in the tree's order."""
    rng = np.random.default_rng(5)
    C, labels = _cloud(400, rng)
    codes = np.unique(labels, return_inverse=True)[1]
    reordered = rng.permutation(len(C))
    q0, s0, l0 = t2_build._within_label_partners(C, codes, 0.2)
    q1, s1, l1 = t2_build._within_label_partners(C[reordered], codes[reordered], 0.2)
    back = np.argsort(reordered)
    for p in range(len(C)):
        a = q0[s0[p]:s0[p] + l0[p]]
        assert np.all(np.diff(a) > 0)  # each position's list is sorted, whatever order the tree found them in
        b = np.sort(reordered[q1[s1[back[p]]:s1[back[p]] + l1[back[p]]]])
        assert np.array_equal(a, b)


def _patchy_after(n, t, rng, wave=0.3):
    """Coordinates and an after-bracket flag thresholded from a smooth random field (patches of about `wave`), a share t
    of the cells after, whatever each cell's label: segregation that coarse type does not explain."""
    C = rng.uniform(size=(n, 3)) * np.array([1.0, 1.0, 0.4])
    w = rng.normal(size=(20, 3)) * 2 * np.pi / wave
    field = np.cos(C @ w.T + rng.uniform(0, 2 * np.pi, 20)).sum(1)
    return C, field > np.quantile(field, 1.0 - t)


def test_the_chain_reaches_its_plateau_within_its_fixed_rounds(n15):
    """harness.NBHD_SWEEPS is a size of the computation: more rounds change the bracket mixing by no more than noise."""
    rng = np.random.default_rng(11)
    n = 5000
    C, after = _patchy_after(n, 0.5, rng)
    labels = rng.choice(["a", "b", "c"], size=n)
    neighbours = cKDTree(C).query(C, k=K)[1]
    built = n15.after_fraction(after, neighbours).std()
    radius = harness.NBHD_RADIUS_FACTOR * t2_build.neighbourhood_radius(C, K)
    sds = {}
    for sweeps in (5, harness.NBHD_SWEEPS, 150):
        at = t2_build.local_type_permutation(C, labels, radius, sweeps, np.random.default_rng([7, harness.NBHD_STREAM]))
        sds[sweeps] = n15.after_fraction(after[np.argsort(at)], neighbours).std()
    assert sds[harness.NBHD_SWEEPS] < 0.8 * built  # the chain removes most of the patchiness at this radius
    assert abs(sds[harness.NBHD_SWEEPS] - sds[150]) < 0.02  # 30 rounds are on the plateau 150 rounds reach
    assert sds[5] >= sds[150] - 0.02  # and the plateau is reached from above: more rounds never hurt the mixing


# --------------------------------------------------------------------------------- the audit's ratios and rule
def test_the_audit_ratios_tell_chance_from_within_type_segregation_from_coarse_type(n15):
    rng = np.random.default_rng(12)
    n = 4000
    # chance: flags independent of position, type and each other
    C = rng.uniform(size=(n, 3)) * np.array([1.0, 1.0, 0.4])
    labels = rng.choice(["a", "b", "c"], size=n)
    r = n15.audit_draw(C, rng.random(n) < 0.5, labels, 0.5, 1_000_000)
    assert 0.95 < r["ratio_null"] < 1.05 and 0.95 < r["ratio_binomial"] < 1.1 and 0.95 < r["null_over_binomial"] < 1.1
    assert 0.9 < r["after_over_built"] < 1.1
    # patches of one bracket within every type: the ratio to the null exceeds the bar and the candidate removes them
    C, after = _patchy_after(n, 0.5, rng)
    labels = rng.choice(["a", "b", "c"], size=n)
    r = n15.audit_draw(C, after, labels, 0.5, 1_000_000)
    assert r["ratio_null"] > n15.BAR and r["ratio_binomial"] > n15.BAR and r["after_over_built"] < 0.85
    # coarse type alone explains the segregation (two compartments, each wholly from one bracket): the ratio to the
    # binomial is large but the within-type null reproduces it exactly, so the ratio that decides stays at 1
    C = rng.uniform(size=(n, 3)) * np.array([1.0, 1.0, 0.4])
    left = C[:, 0] < 0.5
    labels = np.where(left, "a", "b")
    r = n15.audit_draw(C, left, labels, float(left.mean()), 1_000_000)
    assert r["ratio_binomial"] > n15.BAR and r["ratio_null"] == pytest.approx(1.0)
    assert r["after_over_built"] == pytest.approx(1.0)


def test_after_over_built_reads_each_row_at_its_new_position(n15, monkeypatch):
    """`at[r]` is the position row r moves to, so the flags at position p are the flags of the row whose `at` is p."""
    rng = np.random.default_rng(13)
    C, labels = _cloud(120, rng)
    flags = rng.random(120) < 0.4
    at = rng.permutation(120)  # a known arrangement in place of the chain's
    monkeypatch.setattr(n15.harness, "nbhd_rearrangement", lambda coords, labs, seed: at)
    r = n15.audit_draw(C, flags, labels, 0.4, 1_000_000)
    neighbours = cKDTree(C).query(C, k=K)[1]
    at_position = np.argsort(at)  # the row at each position
    want = flags[at_position][neighbours].mean(1).std() / flags[neighbours].mean(1).std()
    assert r["after_over_built"] == pytest.approx(want)


@pytest.mark.parametrize("test, a, b, passes", [
    (1.6, 1.6, 1.0, True),
    (1.6, 1.0, 1.6, True),
    (1.6, 1.6, 1.6, True),
    (1.5, 1.6, 1.6, False),  # the bar is exceeded strictly
    (1.6, 1.5, 1.4, False),
    (1.2, 1.9, 1.9, False),  # never without the test pair
    (1.0, 1.0, 1.0, False),
])
def test_the_pass_rule_needs_the_test_pair_and_one_backtest(n15, test, a, b, passes):
    verdict = n15.verdict({n15.TEST_PAIR: test, n15.BACKTEST_PAIRS[0]: a, n15.BACKTEST_PAIRS[1]: b})
    assert verdict.startswith("pass" if passes else "close")
    assert "E15" in verdict


def test_the_pairs_are_the_harness_tasks_and_fractions(n15):
    tasks = {t.key: t for t in harness.load_tasks()[n15.EMBRYO_BOARD]}
    for key, sb, sa, t in n15.PAIRS:
        if key == n15.TEST_PAIR:
            assert (sb, sa, t) == ("E7.5", "E8.0", 0.5)
            assert t == (7.75 - 7.5) / (8.0 - 7.5)
            continue
        task = tasks[key]
        assert harness.stem(task.inputs["before"]) == sb and harness.stem(task.inputs["after"]) == sa
        days = harness.STAGE_DAYS
        lo, hi, target = days[task.inputs["before"]], days[task.inputs["after"]], days[task.target]
        assert t == (target - lo) / (hi - lo)  # exactly the fraction vec.harness computes for the same draw
    assert set(n15.BACKTEST_PAIRS) == set(tasks)
    assert n15.K == harness.NBHD_K == 15


# --------------------------------------------------------------------------------- the script's outputs, on fakes
ALL_GENES = [f"g{i}" for i in range(12)]
VAL_GENES = ALL_GENES[:8]
TEST_GENES = ALL_GENES[2:10]
SIZES = {"E6.75": 1229, "E7.25": 1251, "E7.5": 1263, "E8.0": 1241}  # rows as released of the synthetic stages
SECRET = {"min_cells": 1777, "max_cells": 7321, "truth_cells": 2468, "ref_cells": 1357}  # stand for a test board's


def _boards():
    val = harness.pseudo_counts(SIZES["E7.5"], SIZES["E7.25"])  # the validation board's published counts
    return {"T2:embryo:val_interp": Board(key="T2:embryo:val_interp", task="T2", genes=VAL_GENES, needs_coords=True,
                                          split="val_interp", **val),
            "T2:embryo:test_interp": Board(key="T2:embryo:test_interp", task="T2", genes=TEST_GENES, needs_coords=True,
                                           split="test_interp", **SECRET)}


def _stage(name, genes):
    n = SIZES[name]
    rng = np.random.default_rng(n)
    a = ad.AnnData(X=rng.poisson(1.0, size=(n, len(genes))).astype(np.float32),
                   obs=pd.DataFrame({"group": rng.choice(["Neural", "Mesoderm", "Blood"], size=n)},
                                    index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=genes))
    a.obsm[COORD_KEY] = rng.normal(size=(n, 3))
    a.uns["qc"] = {"n_input": n, "drop_duplicates": True}
    return a


def _fake_load(calls):
    def load_stage(path, genes=None, drop_duplicates=False):
        name = Path(path).name.removesuffix(".h5ad")
        calls.append((name, list(genes), drop_duplicates))
        return _stage(name, list(genes))

    return load_stage


def _fake_interpolate(patchy_t, seen):
    """In place of the FGW interpolation: drawn cells in a box whose bracket flags are patchy for the fractions in
    `patchy_t` and random otherwise; records the cell count each call asks for, which the outputs must not contain."""
    def interpolate_clouds(Cb, Ca, Eb, Ea, t, n, seed, target_rms=1.0):
        seen.append(n)
        rng = np.random.default_rng([seed, int(round(t * 1000))])
        C = rng.uniform(size=(n, 3)) * np.array([1.0, 1.0, 0.4])
        if any(abs(t - p) < 1e-9 for p in patchy_t):
            w = rng.normal(size=(20, 3)) * 2 * np.pi / 0.3
            field = np.cos(C @ w.T + rng.uniform(0, 2 * np.pi, 20)).sum(1)
            flag = field > np.quantile(field, 1.0 - t)
        else:
            flag = rng.random(n) < t
        stage = flag.astype(np.int64)
        rows = np.where(flag, rng.integers(0, len(Ca), n), rng.integers(0, len(Cb), n))
        return t2_build.Interpolation(stage=stage, rows=rows, coords=C, pool_coords=C, pool_weights=np.ones(n) / n,
                                      record={"n": n})

    return interpolate_clouds


def _patch(n15, monkeypatch, patchy_t, tmp_path):
    calls, seen = [], []
    monkeypatch.setattr(n15, "load_boards", _boards)
    monkeypatch.setattr(n15, "load_stage", _fake_load(calls))
    monkeypatch.setattr(n15.provenance, "git_state", lambda repo: ("cafebabe", []))
    monkeypatch.setattr(n15.t2_build, "interpolate_clouds", _fake_interpolate(patchy_t, seen))
    data = tmp_path / "data"
    data.mkdir()
    for s in SIZES:
        (data / f"{s}.h5ad").write_bytes(s.encode())
    return data, tmp_path / "out", calls, seen


def _numbers(text: str) -> set[int]:
    """Every whole number in a text, thousands separators dropped."""
    return {int(d) for d in re.findall(r"\d+", re.sub(r"(?<=\d),(?=\d{3})", "", text))}


def _leaves(x, ints: set, texts: list) -> None:
    """Every whole number (booleans aside) and every string (dict keys included) in a parsed JSON document."""
    if isinstance(x, dict):
        for k, v in x.items():
            texts.append(str(k))
            _leaves(v, ints, texts)
    elif isinstance(x, list):
        for v in x:
            _leaves(v, ints, texts)
    elif isinstance(x, int) and not isinstance(x, bool):
        ints.add(x)
    elif isinstance(x, str):
        texts.append(x)


def test_the_script_writes_ratios_only_and_decides_by_the_row(n15, tmp_path, monkeypatch, capsys):
    # patchy on the test pair (t = 1/2) and on one backtest (t = 2/3): the rule says pass
    data, out, calls, seen = _patch(n15, monkeypatch, (0.5, 2 / 3), tmp_path)
    monkeypatch.setattr(n15, "SEEDS", (0, 1))  # a full run in two draws
    n15.main(["--out", str(out), "--data", str(data), "--seeds", "0", "1"])
    printed = capsys.readouterr().out
    # the test pair on the test board's panel, the backtests on the validation board's, deduplicated, once per panel
    want = [(s, VAL_GENES, True) for s in SIZES] + [(s, TEST_GENES, True) for s in ("E7.5", "E8.0")]
    assert sorted(calls) == sorted(want)
    # the draws asked for the count rule's numbers: the test board's limit on the test pair, read by code only
    assert seen == [5000, 5000] * 2 + [SECRET["max_cells"]] * 2
    rep = json.loads((out / "n15_nbhd_audit.json").read_text())
    txt = (out / "n15_nbhd_audit.txt").read_text()
    assert rep["verdict"].startswith("pass") and txt.splitlines()[-1] == f"verdict: {rep['verdict']}"
    assert printed.endswith(txt)
    medians = {k: p["summary"]["ratio_null"]["median"] for k, p in rep["pairs"].items()}
    assert medians[n15.TEST_PAIR] > n15.BAR and medians[n15.BACKTEST_PAIRS[0]] > n15.BAR
    assert medians[n15.BACKTEST_PAIRS[1]] < 1.2
    assert rep["pairs"][n15.TEST_PAIR]["exceeds_bar"] and rep["pairs"][n15.BACKTEST_PAIRS[0]]["exceeds_bar"]
    assert rep["stamp"]["commit"] == "cafebabe" and rep["stamp"]["tree_clean"] is True and rep["stamp"]["utc"]
    assert set(rep["inputs"]["files_sha256"]) == {f"T2/{s}.h5ad" for s in SIZES}
    assert all(len(h) == 64 for h in rep["inputs"]["files_sha256"].values())
    assert rep["inputs"]["pseudo_tasks_sha256"] and rep["inputs"]["lineage_map_sha256"]
    assert len(rep["inputs"]["panels_sha256"]) == 2
    # nothing written or printed is a cell count: the only whole numbers are seeds and the row's constants
    ints: set[int] = set()
    texts: list[str] = []
    _leaves(rep, ints, texts)
    assert ints <= {0, 1, 1_000_000, 1_000_001, K, n15.NULL_PERMUTATIONS, harness.NBHD_SWEEPS}
    forbidden = {*SECRET.values(), *SIZES.values(), 5000, *seen}
    for text in (*texts, txt, printed):
        assert not _numbers(text) & forbidden, text
    # the text gives every pair's ratios and the rule
    for key in rep["pairs"]:
        assert key in txt
    assert "ratio_null" in txt and n15.RULE.splitlines()[0][:40] in " ".join(txt.split())
    # and a rerun refuses, before it reads anything
    calls.clear()
    with pytest.raises(SystemExit, match="already written"):
        n15.main(["--out", str(out), "--data", str(data), "--seeds", "0", "1"])
    assert calls == []


def test_the_script_closes_the_line_when_the_test_pair_does_not_segregate(n15, tmp_path, monkeypatch, capsys):
    data, out, _, _ = _patch(n15, monkeypatch, (2 / 3,), tmp_path)  # only a backtest patchy: not enough
    monkeypatch.setattr(n15, "SEEDS", (0,))
    n15.main(["--out", str(out), "--data", str(data), "--seeds", "0"])
    capsys.readouterr()
    rep = json.loads((out / "n15_nbhd_audit.json").read_text())
    assert rep["verdict"].startswith("close") and "the test pair" in rep["verdict"]


def test_a_partial_run_gives_no_verdict(n15, tmp_path, monkeypatch, capsys):
    data, out, _, _ = _patch(n15, monkeypatch, (0.5,), tmp_path)
    n15.main(["--out", str(out), "--data", str(data), "--seeds", "0"])  # one of the rule's ten seeds
    capsys.readouterr()
    rep = json.loads((out / "n15_nbhd_audit.json").read_text())
    assert rep["verdict"] == "partial run: no verdict"


SMALL = {"E6.75": 229, "E7.25": 251, "E7.5": 263, "E8.0": 241}  # tiny stages for the real interpolation's test


def test_the_script_runs_end_to_end_on_the_real_interpolation(n15, tmp_path, monkeypatch, capsys):
    """The script's own call of vec.t2_build.interpolate_clouds and the candidate's re-arrangement, on tiny synthetic
    stages (one seed, so no verdict; the fakes above decide the rule's paths)."""
    calls = []
    monkeypatch.setitem(SIZES, "E6.75", SMALL["E6.75"])
    monkeypatch.setitem(SIZES, "E7.25", SMALL["E7.25"])
    monkeypatch.setitem(SIZES, "E7.5", SMALL["E7.5"])
    monkeypatch.setitem(SIZES, "E8.0", SMALL["E8.0"])  # _stage reads the sizes through this dictionary
    monkeypatch.setattr(n15, "load_boards", _boards)
    monkeypatch.setattr(n15, "load_stage", _fake_load(calls))
    monkeypatch.setattr(n15.provenance, "git_state", lambda repo: ("cafebabe", []))
    data = tmp_path / "data"
    data.mkdir()
    for s in SIZES:
        (data / f"{s}.h5ad").write_bytes(s.encode())
    out = tmp_path / "out"
    n15.main(["--out", str(out), "--data", str(data), "--seeds", "0"])
    capsys.readouterr()
    rep = json.loads((out / "n15_nbhd_audit.json").read_text())
    assert rep["verdict"] == "partial run: no verdict" and set(rep["pairs"]) == {p[0] for p in n15.PAIRS}
    for p in rep["pairs"].values():
        assert len(p["draws"]) == 1 and np.isfinite(p["draws"][0]["ratio_null"])


def test_draw_count_uses_the_harness_rule_on_each_board(n15):
    boards = _boards()
    tasks = {t.key: t for t in harness.load_tasks()[n15.EMBRYO_BOARD]}
    before = _stage("E7.5", TEST_GENES)
    assert n15.draw_count(n15.TEST_PAIR, before, boards, tasks, {}) == SECRET["max_cells"]
    rows = {s: SIZES[s] for s in SIZES}
    # the validation board's published limits, and the twin board's counts from the files (the harness's rule)
    want = boards["T2:embryo:val_interp"].max_cells
    assert n15.draw_count("T2e/E7.25+E8.0->E7.5", before, boards, tasks, rows) == want
    want = harness.pseudo_counts(SIZES["E7.25"], SIZES["E6.75"])["max_cells"]
    assert n15.draw_count("T2e/E6.75+E7.5->E7.25", before, boards, tasks, rows) == want


# --------------------------------------------------------------------------------- the candidate fgw_perm (E15)
EMBRYO_TASK = PseudoTask(key="syn_embryo", test_board="T2:embryo:test_interp", task="T2", setting="embryo",
                         kind="interp", target="T2/E7.5.h5ad", reference="T2/E7.25.h5ad",
                         inputs={"before": "T2/E7.25.h5ad", "after": "T2/E8.0.h5ad", "earlier": "T2/E6.75.h5ad"},
                         twin_board="T2:embryo:val_interp")


def _pool(n, rng, means, props, scale=1.0):
    centres = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0]], float) * scale
    st = simulate_stage(n, means, props, rng, coords_scale=1.5 * scale, coords_shift=centres)
    a = ad.AnnData(X=st.X, obs=pd.DataFrame({"celltype": st.celltype}, index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.obsm[COORD_KEY] = st.coords.astype(np.float32)
    a.obs["group"] = np.where(st.celltype == "type0", "Neural", "Mesoderm")  # coarse type: the lineage map's group
    return a


@pytest.fixture(scope="module")
def built():
    rng = np.random.default_rng(31)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    earlier = _pool(600, rng, base * 0.9, [0.6, 0.3, 0.1], scale=0.7)
    before = _pool(800, rng, base, [0.5, 0.3, 0.2])
    after = _pool(1100, rng, base * 1.6, [0.3, 0.4, 0.3], scale=1.7)
    target = _pool(1500, rng, base * 1.3, [0.4, 0.35, 0.25], scale=1.2)
    board = Board(key="pseudo:synthetic_embryo", task="T2", genes=list(GENES), needs_coords=True, setting="embryo",
                  min_cells=100, max_cells=700)
    ctx = harness.BuildContext(task=harness.TaskView.of(EMBRYO_TASK), board=board,
                               inputs={"before": before, "after": after, "earlier": earlier},
                               seed=harness.BUILD_SEED_OFFSET + 3)
    fgw = CANDIDATES["fgw"].build(ctx)
    perm = CANDIDATES["fgw_perm"].build(ctx)
    return fgw, perm, target, before


def test_fgw_perm_is_fgw_with_only_its_positions_rearranged_among_the_drawn_cells(built):
    fgw, perm, *_ = built
    assert CANDIDATES["fgw_perm"].applies(EMBRYO_TASK)
    tasks = {t.key: t for ts in harness.load_tasks().values() for t in ts}
    assert CANDIDATES["fgw_perm"].applies(tasks["T2e/E7.25+E8.0->E7.5"])
    assert CANDIDATES["fgw_perm"].applies(tasks["T2e/E6.75+E7.5->E7.25"])
    assert not CANDIDATES["fgw_perm"].applies(tasks["T2h/<=E9.5->E10.5"])
    np.testing.assert_array_equal(perm.X, fgw.X)  # the same cells and expression rows, bit for bit, in the same order
    # the same positions as a multiset; some rows really moved
    assert perm.n == fgw.n and not np.array_equal(perm.coords, fgw.coords)
    a, b = fgw.coords, perm.coords
    np.testing.assert_array_equal(a[np.lexsort(a.T)], b[np.lexsort(b.T)])
    # every row stays within the fixed radius of where fgw put it
    radius = harness.NBHD_RADIUS_FACTOR * t2_build.neighbourhood_radius(np.asarray(a, float), harness.NBHD_K)
    assert np.linalg.norm(np.asarray(b, float) - np.asarray(a, float), axis=1).max() <= radius * 1.001


def _stage_of(a: ad.AnnData) -> Stage:
    return Stage(X=np.asarray(a.X, np.float32), coords=np.asarray(a.obsm[COORD_KEY], np.float32),
                 celltype=np.asarray(a.obs["celltype"]).astype(str))


def test_only_the_neighbourhood_metric_moves(built):
    """TRIES row E15: the permutation leaves X untouched, so the expression metrics are bit-identical, and the whole
    multiset of positions is unchanged, so scale_log_ratio is identical up to the scorer's rounding. The scorer's
    d2_shape and occupancy_dice read fixed-seed row subsets of size min(n_pred, n_true) (common.shape_metrics._match_n):
    when the truth has at least as many cells as the prediction, every position is read and they are identical too;
    when it has fewer, the fixed subset of rows carries other positions after the permutation (at most the radius away),
    a small per-seed wobble with no expected direction."""
    fgw, perm, target, ref = built
    truth_large, ref_stage = _stage_of(target), _stage_of(ref)
    truth_small = truth_large.take(np.arange(300))  # fewer truth cells than the prediction's, as on the real boards
    assert truth_small.n < fgw.n < truth_large.n
    for truth, every_position_read in ((truth_large, True), (truth_small, False)):
        m0 = panel_metrics("T2", fgw, truth, ref_stage, seed=0)
        m1 = panel_metrics("T2", perm, truth, ref_stage, seed=0)
        for name in ("de_score", "de_direction", "mmd_u", "energy_distance", "variogram", "pb_rel_err"):
            assert m0[name] == m1[name], name  # bit-identical: the rows of X did not change
        assert m1["scale_log_ratio"] == pytest.approx(m0["scale_log_ratio"], abs=1e-3)
        assert m1["neighborhood_mmd"] != m0["neighborhood_mmd"]
        if every_position_read:
            assert m1["occupancy_dice"] == pytest.approx(m0["occupancy_dice"], abs=1e-3)
            assert m1["sliced_wasserstein"] == pytest.approx(m0["sliced_wasserstein"], abs=1e-3)
        else:
            assert abs(m1["occupancy_dice"] - m0["occupancy_dice"]) < 0.08
        assert abs(m1["d2_shape"] - m0["d2_shape"]) < 0.02  # 200,000 random pairs: at most a small drift
