"""Step 0 of row R48's file T2embryo_s2x_ovsb.h5ad (row E56): experiments/coords/ovsb_step0.py, on synthetic
brackets written as released files (read back through vec.prep.load_stage, as the script reads the real ones) and a
throwaway git repository standing in for the ledger and the clean-tree rule. Nothing is read from data/ or the real
ledger, and the script never runs here on released data.

Pinned: the factor applied in expectation is apply_gene_factor's mean; the carry check holds when the embryonic labels
share the run's factor and stops when they do not, or when no embryonic label qualifies; the depth bound is 10% either
side, inclusive; the outcome line is one of the four fixed strings; reading.txt starts with it and reading.json holds
the same outcome, the statement, the method and the commit; the script prints only "wrote <path>"; and it refuses
without row E56, on an unclean tree and over an existing reading, before any data are read."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from conftest import git
from test_series_back_build import PANEL, RUN, _stage

from vec import core
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "experiments" / "coords" / "ovsb_step0.py"
# the lineage map's own labels (configs/lineage_map.yaml), so the groups come from vec.prep.load_stage as on the files
WORLD = {"EXE-Endoderm": ("Endoderm", 300, 250, 0), "ExEM-1": ("ExtraEmbryonicMeso", 300, 250, 1),
         "Anterior Epiblast": ("EpiblastStreak", 300, 200, 2), "Primitive Streak": ("EpiblastStreak", 300, 200, 2),
         "Foregut": ("Endoderm", 0, 100, 1)}
BUMP = {"EXE-Endoderm": [0, 4], "ExEM-1": [1, 5], "Anterior Epiblast": [8, 10], "Primitive Streak": [2, 9],
        "Foregut": [3, 11]}
DEV = {"Primitive Streak": np.r_[np.ones(8), np.full(2, 4.0), np.ones(5)]}
EMBRYONIC = ["Anterior Epiblast", "Primitive Streak"]
NO_RUN = {lab: (RUN if lab in ("EXE-Endoderm", "ExEM-1") else np.ones(15)) for lab in WORLD}  # embryonic: no run


@pytest.fixture(scope="module")
def s0():
    spec = importlib.util.spec_from_file_location("ovsb_step0", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _pair(run=None, world=WORLD, seed=0):
    kw = dict(world=world, bump=BUMP, dev=DEV)
    return _stage("E7.5", 101 + seed, **kw), _stage("E8.0", 102 + seed, run=run, **kw)


def _ledger_repo(tmp_path: Path, row: bool = True) -> Path:
    repo = tmp_path / "checkout"
    (repo / "experiments" / "harness").mkdir(parents=True)
    rows = "| E55 (synthetic) | T2 embryo |\n" + ("| E56 (step 0 of row R48's file) | T2 embryo |\n" if row else "")
    (repo / "experiments" / "harness" / "TRIES.md").write_text("| ID | what |\n|---|---|\n" + rows)
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "ledger")
    return repo


@pytest.fixture
def world(s0, tmp_path, monkeypatch):
    """A clean throwaway checkout with row E56, the brackets written under data/T2/, the panel; returns a setter that
    writes another pair of brackets."""
    repo = _ledger_repo(tmp_path)
    monkeypatch.setattr(s0, "REPO", repo)
    monkeypatch.setattr(s0, "TRIES", repo / "experiments" / "harness" / "TRIES.md")
    monkeypatch.setattr(s0, "load_boards", lambda: {s0.BOARD: SimpleNamespace(genes=list(PANEL))})
    data = tmp_path / "data"
    (data / "T2").mkdir(parents=True)

    def write(run=None, seed=0):
        before, after = _pair(run=run, seed=seed)
        before.write_h5ad(data / "T2" / "E7.5.h5ad")
        after.write_h5ad(data / "T2" / "E8.0.h5ad")
        return data

    write()
    return SimpleNamespace(repo=repo, data=data, write=write)


def _run(s0, tmp_path, data, capsys, *extra):
    out = tmp_path / "results" / "ovsb_step0"
    s0.main(["--out", str(out), "--data", str(data), *extra])
    printed = capsys.readouterr().out.splitlines()
    return out, printed


# ------------------------------------------------------------------------------------------------ the checks
def test_the_factor_in_expectation_is_apply_gene_factors_mean(s0):
    rng = np.random.default_rng(1)
    C = rng.poisson(rng.gamma(2.0, 3.0, size=(6, 12)), size=(6, 12)).astype(np.int32)
    C[:, 0] += 5  # no cell left empty: the rule's "an empty cell keeps its counts" never fires
    phi = np.r_[0.25, 0.5, 0.8, 1.0, 1.0, 1.3, 2.0, 4.0, 0.6, 1.0, 3.0, 0.9]
    P = rng.dirichlet(np.ones(12), size=6)
    draws = np.mean([core.apply_gene_factor(C, phi, P, np.random.default_rng([7, k])) for k in range(4000)], axis=0)
    want = s0.expected_counts(C, phi, P)
    np.testing.assert_allclose(draws, want, rtol=0.06, atol=0.06)
    np.testing.assert_array_equal(want[:, phi == 1.0], C[:, phi == 1.0])


def test_the_carry_check_holds_when_the_embryonic_labels_share_the_run(s0):
    before, after = _pair()
    phi, matched = tb.series_back_factor(before, after, PANEL)
    assert matched == ["EXE-Endoderm", "ExEM-1"]
    ii = s0.carry_check(before, after, PANEL, phi, matched)
    assert ii["labels"] == EMBRYONIC and ii["raised"] is False and ii["sd_with"] < ii["sd_without"]
    assert sum(r["weight"] for r in ii["per_label"]) == pytest.approx(1.0)
    assert [r["weight"] for r in ii["per_label"]] == pytest.approx([0.5, 0.5])  # E8.0's 200 and 200 cells
    # the figure is the E8.0-cell-weighted mean of each label's SD over genes of the log ratio of mean proportions
    eps = tb.SERIES_EPS
    C75, C80 = tb._stage_counts(before, PANEL)[1], tb._stage_counts(after, PANEL)[1]
    n75, n80 = (np.asarray(x.obs["celltype"]).astype(str) for x in (before, after))
    sds = [np.log((s0.mean_props(C80, np.flatnonzero(n80 == lab)) + eps)
                  / (s0.mean_props(C75, np.flatnonzero(n75 == lab)) + eps)).std() for lab in EMBRYONIC]
    assert ii["sd_without"] == pytest.approx(np.mean(sds), rel=1e-12)


def test_the_carry_check_stops_when_the_embryonic_labels_do_not_share_it(s0):
    before, after = _pair(run=NO_RUN)
    phi, matched = tb.series_back_factor(before, after, PANEL)
    ii = s0.carry_check(before, after, PANEL, phi, matched)
    assert ii["labels"] == EMBRYONIC and ii["raised"] is True and ii["sd_with"] > ii["sd_without"]


def test_a_factor_of_one_changes_no_figure_and_does_not_stop(s0):
    before, after = _pair(run=NO_RUN)
    _, matched = tb.series_back_factor(before, after, PANEL)
    ii = s0.carry_check(before, after, PANEL, np.ones(len(PANEL)), matched)
    assert ii["sd_with"] == ii["sd_without"] and ii["raised"] is False
    assert all(r["sd_with"] == r["sd_without"] for r in ii["per_label"])


def test_the_carry_check_stops_without_an_embryonic_label(s0):
    world = {k: v for k, v in WORLD.items() if k in ("EXE-Endoderm", "ExEM-1")}
    world["Anterior Epiblast"] = ("EpiblastStreak", 300, 99, 2)  # one cell short at E8.0
    before, after = _pair(world=world)
    phi, matched = tb.series_back_factor(before, after, PANEL)
    ii = s0.carry_check(before, after, PANEL, phi, matched)
    assert ii["labels"] == [] and ii["raised"] is True and "cannot pass" in ii["reason"]
    assert s0.embryonic_labels(before.obs["celltype"], after.obs["celltype"], matched) == []


def test_the_depth_bound_and_the_outcome_lines(s0):
    assert s0.depth_within(110.0, 100.0) and s0.depth_within(90.0, 100.0) and s0.depth_within(100.0, 100.0)
    assert not s0.depth_within(110.01, 100.0) and not s0.depth_within(89.99, 100.0)
    lines = {(a, b): s0.outcome(a, b) for a in (False, True) for b in (False, True)}
    assert lines == {(False, False): "Step 0: PASS", (True, False): "Step 0: STOP (ii)",
                     (False, True): "Step 0: STOP (iii)", (True, True): "Step 0: STOP (ii, iii)"}
    assert set(lines.values()) == set(s0.OUTCOMES)


# ------------------------------------------------------------------------------------------------ the run
def test_a_passing_run_writes_the_outcome_line_then_the_figures(s0, world, tmp_path, capsys):
    out, printed = _run(s0, tmp_path, world.data, capsys)
    assert printed == [f"wrote {out / 'reading.txt'}", f"wrote {out / 'reading.json'}"]  # nothing else is printed
    text = (out / "reading.txt").read_text()
    res = json.loads((out / "reading.json").read_text())
    lines = text.splitlines()
    assert lines[0] == "Step 0: PASS" == res["outcome"]
    assert lines[1].startswith("(i) reported only:") and lines[2].startswith("(ii) carry check over 2 embryonic")
    assert any(ln.startswith("(iii) depth:") for ln in lines)
    assert s0.STATEMENT in text and res["statement"] == s0.STATEMENT and res["method"] in text
    assert "STOP if it raises the cell-weighted SD over genes of log(E8.0 mean proportions / E7.5's)" in s0.STATEMENT
    assert "within 10% of E7.5's" in s0.STATEMENT
    assert res["row"] == "E56" and res["ruling"] == "R48" and res["file"] == "T2embryo_s2x_ovsb.h5ad"
    p = res["provenance"]
    assert p["commit"] == git(world.repo, "rev-parse", "HEAD") and p["clean"] is True
    assert p["inputs"]["after"]["file"] == "T2/E8.0.h5ad" and re.fullmatch("[0-9a-f]{64}", p["script_sha256"])
    # the figures are the build's: its factor and labels, the rendering of every E8.0 cell at the build's seed
    assert res["i"]["matched_labels"] == ["EXE-Endoderm", "ExEM-1"] and res["i"]["n_genes"] == len(PANEL)
    assert res["iii"]["within"] is True and res["iii"]["seed"] == [1_000_001, 7]
    assert res["iii"]["ratio_to_e75"] == pytest.approx(res["iii"]["rendered_median_e80"] / res["iii"]["median_e75"])
    assert res["ii"]["raised"] is False
    # no cell number of a bracket is written: shares and medians only
    assert not re.search(r"(?<![\d.])(1200|1100|1000|600|500|400|250|200|100|99)(?![\d.%])",
                         text.split("Statement")[0].replace("10%", ""))


def _rendered_at(s0, monkeypatch, scale: float) -> None:
    """The rendering as the script calls it, its rendered median put at `scale` x E7.5's (to set check (iii))."""
    inner = s0.tb.series_back_expression

    def at(*a, **kw):
        X, rec = inner(*a, **kw)
        return X, {**rec, "rendered_median": rec["median_depth"]["before"] * scale}

    monkeypatch.setattr(s0.tb, "series_back_expression", at)


def test_a_run_whose_embryonic_labels_lack_the_run_stops_at_ii(s0, world, tmp_path, capsys, monkeypatch):
    data = world.write(run=NO_RUN)
    _rendered_at(s0, monkeypatch, 1.0)
    out, _ = _run(s0, tmp_path, data, capsys)
    assert (out / "reading.txt").read_text().splitlines()[0] == "Step 0: STOP (ii)"
    res = json.loads((out / "reading.json").read_text())
    assert res["outcome"] == "Step 0: STOP (ii)" and res["ii"]["raised"] is True and res["iii"]["within"] is True
    assert "raised (STOP)" in (out / "reading.txt").read_text()


@pytest.mark.parametrize("run, scale, want", [(None, 1.11, "Step 0: STOP (iii)"), (None, 0.89, "Step 0: STOP (iii)"),
                                              (NO_RUN, 1.11, "Step 0: STOP (ii, iii)")])
def test_a_rendered_depth_off_by_more_than_ten_percent_stops_at_iii(s0, world, tmp_path, capsys, monkeypatch, run,
                                                                    scale, want):
    data = world.write(run=run)
    _rendered_at(s0, monkeypatch, scale)
    out, _ = _run(s0, tmp_path, data, capsys)
    assert (out / "reading.txt").read_text().splitlines()[0] == want
    assert "outside (STOP)" in (out / "reading.txt").read_text()


def _no_load(s0, monkeypatch):
    def refuse(*a, **kw):
        raise AssertionError("a bracket was read")

    monkeypatch.setattr(s0, "load_stage", refuse)


def test_without_row_e56_nothing_is_read_or_written(s0, world, tmp_path, monkeypatch, capsys):
    repo = _ledger_repo(tmp_path / "other", row=False)
    monkeypatch.setattr(s0, "REPO", repo)
    monkeypatch.setattr(s0, "TRIES", repo / "experiments" / "harness" / "TRIES.md")
    _no_load(s0, monkeypatch)
    with pytest.raises(SystemExit, match="row E56 is not in"):
        _run(s0, tmp_path, world.data, capsys)
    assert not (tmp_path / "results").exists()


def test_an_unclean_tree_is_refused(s0, world, tmp_path, monkeypatch, capsys):
    (world.repo / "experiments" / "harness" / "TRIES.md").write_text(
        (world.repo / "experiments" / "harness" / "TRIES.md").read_text() + "| E57 (uncommitted) | x |\n")
    _no_load(s0, monkeypatch)
    with pytest.raises(SystemExit, match="not clean"):
        _run(s0, tmp_path, world.data, capsys)
    assert not (tmp_path / "results").exists()


def test_a_reading_is_never_overwritten(s0, world, tmp_path, monkeypatch, capsys):
    out = tmp_path / "results" / "ovsb_step0"
    out.mkdir(parents=True)
    (out / "reading.json").write_text("{}")
    _no_load(s0, monkeypatch)
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        _run(s0, tmp_path, world.data, capsys)
    assert (out / "reading.json").read_text() == "{}" and not (out / "reading.txt").exists()


def test_the_script_names_the_launchers_path_and_is_in_the_step_lists():
    src = SCRIPT.read_text()
    assert 'RESULTS = "experiments/coords/results/ovsb_step0"' in src
    assert "READING_DIR=experiments/coords/results/ovsb_step0" in (
        REPO / "experiments" / "builds" / "embryo_ovsb.sh").read_text()
    assert "tests/test_ovsb_step0.py" in (REPO / "experiments" / "external" / "n54_step.sh").read_text()
