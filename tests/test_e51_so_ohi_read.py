"""Row E51 as specified (job_embryo_so_ohi): the harness arms fgw_shape0x_series_ov (SO) and fgw_shape0x_ov_hi (O_hi),
the read script experiments/coords/e51_so_ohi_read.py (the gate, each pass clause, the void conditions, the tables
that decide nothing) on hand-made per-seed records, and its launcher's refusal paths under a temporary home. Nothing
is read from data/, nothing is built from a stage and nothing is scored. Two structural checks open committed
practice result files for their keys only (no value is compared with anything)."""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments" / "coords"))
import e51_so_ohi_read as e51  # noqa: E402

C = "c" * 40
TASKS = e51.nt.all_tasks()
NAMES = (harness.FLOOR, *e51.RUN_NAMES)
METRICS = (*e51.EXPR_METRICS, *e51.POS_METRICS, e51.NBHD)
HI_EXPR = 0.47
SKILLS = {  # expression skills and position skills per file, as the construction makes them equal
    harness.FLOOR: (0.30, 0.30), e51.FGW: (0.40, 0.40), e51.SERIES: (0.45, 0.42), e51.OV: (0.41, 0.55),
    e51.SO: (0.45, 0.55), e51.OHI: (HI_EXPR, 0.55), e51.HI: (HI_EXPR, 0.50)}
SCORES = {harness.FLOOR: 50.0, e51.FGW: 51.0, e51.SERIES: 52.0, e51.OV: 53.0, e51.SO: 54.0, e51.OHI: 54.0,
          e51.HI: 52.5}


def _skills(name):
    e, p = SKILLS[name]
    return {m: (e if m in e51.EXPR_METRICS else p if m in e51.POS_METRICS else 0.6) for m in METRICS}


def _rec(task, seed, scores=None, groups=None, skills=None, commit=C, drop=()):
    builds = {}
    for b in e51.BUILDS:
        sc = {n: SCORES[n] for n in NAMES} | (scores or {})
        gr = {n: {g: 0.5 for g in e51.GROUP_NAMES} for n in NAMES}
        for n, d in (groups or {}).items():
            gr[n] = {**gr[n], **d}
        sk = {n: _skills(n) for n in NAMES}
        for (bb, n), d in (skills or {}).items():
            if bb == b:
                sk[n] = {**sk[n], **d}
        row = {"seed": seed, "scores": sc, "groups": gr, "skills": sk}
        for n in drop:
            for k in ("scores", "groups", "skills"):
                row[k].pop(n, None)
        builds[b] = {"provenance": {"commit": commit, "tree_clean": True}, "build": b, "protocol": "bake",
                     "board": {"split": "pseudo"}, "rows": [row]}
    return {"task": task, "seed": seed, "run_commit": commit, "builds": builds,
            "facts": {f"board_split_{b}": "pseudo" for b in e51.BUILDS}}


def _recs(fn=None):
    """Records for every task and seed; fn(task, seed) -> keyword arguments of _rec."""
    return {t: [_rec(t, s, **((fn or (lambda *_: {}))(t, s))) for s in e51.SEEDS] for t in TASKS}


def _hi(value=HI_EXPR, board=SCORES[e51.HI]):
    """The committed fgw_shape0_hi figures, synthetic: expression skills and board scores per seed."""
    return {t.key: {**{m: np.full(10, value) for m in e51.EXPR_METRICS}, "board": np.full(10, board)}
            for t in TASKS.values()}


def _read(recs, why=(), hi=None, why_hi=()):
    return e51.reading(recs, list(why), C, _hi() if hi is None else hi, list(why_hi))


# ------------------------------------------------------------------------------------------------ the harness arms
def test_the_two_arms_call_the_pipeline_with_the_specifications_arguments(monkeypatch):
    calls = []
    monkeypatch.setattr(harness, "_fgw", lambda ctx, **kw: calls.append(kw) or "stage")
    assert harness.CANDIDATES["fgw_shape0x_series_ov"].build("ctx") == "stage"
    assert harness.CANDIDATES["fgw_shape0x_ov_hi"].build("ctx") == "stage"
    assert calls == [dict(t_pos=0.0, expression="series", positions="floor_exact", shape_form="overlay"),
                     dict(t_pos=0.0, expression="render_hi", positions="floor_exact", shape_form="overlay")]
    assert "floor_exact" in harness.t2_build.POSITIONS and "render_hi" in harness.RENDER_EXPRESSIONS
    assert "overlay" in harness.SHAPE_CANDIDATES.values()


def test_the_registry_entries():
    so, hi = harness.CANDIDATES["fgw_shape0x_series_ov"], harness.CANDIDATES["fgw_shape0x_ov_hi"]
    assert so.name == "fgw_shape0x_series_ov" and hi.name == "fgw_shape0x_ov_hi"
    assert so.params["shape_form"] == hi.params["shape_form"] == harness.CANDIDATES["fgw_shape0x_ov"].params[
        "shape_form"] == "overlay"
    assert so.params["positions"] == hi.params["positions"] == harness.t2_build.POSITIONS[1]
    assert so.requires is harness.CANDIDATES["fgw_shape0x_series"].requires  # series needs the shared labels
    assert hi.requires is harness.CANDIDATES["fgw_shape0x_ov"].requires
    assert so.tasks == hi.tasks == frozenset({"T2"}) and not so.needs_reference_cells and not hi.needs_reference_cells
    assert harness.SO_OHI_CANDIDATES == {"fgw_shape0x_series_ov": "series", "fgw_shape0x_ov_hi": "render_hi"}
    assert not set(harness.SO_OHI_CANDIDATES) & set(harness.SHAPE_CANDIDATES)  # E49's seven stay seven


def _canon(x):
    """A repr that does not depend on the hash seed (sets and mappings sorted, recursively)."""
    if isinstance(x, (set, frozenset)):
        return "{" + ",".join(sorted(_canon(v) for v in x)) + "}"
    if isinstance(x, dict) or hasattr(x, "items"):
        return "{" + ",".join(sorted(f"{_canon(k)}:{_canon(v)}" for k, v in x.items())) + "}"
    if isinstance(x, (list, tuple)):
        return "[" + ",".join(_canon(v) for v in x) + "]"
    return repr(x)


LATER_ROWS = ("GEODESIC_CANDIDATES",)  # harness registries of rows merged after this pin was made


def _old_registry_digest() -> tuple[int, str]:
    """Every candidate but the two new ones: name, roles, params, the builder's and the condition's qualified names
    and the builder's bytecode. Pinned on main before the arms were added. Candidates of rows merged into main beside
    this one (LATER_ROWS: row E53's geodesic pair, which its own tests hold) are left out too, so the pin still says
    that no candidate older than this row changed."""
    h = hashlib.sha256()
    later = set().union(*(set(getattr(harness, k, {})) for k in ("SO_OHI_CANDIDATES", *LATER_ROWS)))
    names = sorted(set(harness.CANDIDATES) - later)
    for n in names:
        c = harness.CANDIDATES[n]
        b = c.build
        h.update(_canon((n, c.tasks, c.needs, c.params,
                         None if b is None else (b.__qualname__, b.__code__.co_code.hex()),
                         None if c.requires is None else c.requires.__qualname__, c.tuned_on, c.externals,
                         c.derived, c.derived_optional, c.needs_reference_cells)).encode())
    return len(names), h.hexdigest()


# computed with _old_registry_digest on main's src/ (5b22c73) before the two arms were added
OLD_REGISTRY_PINNED = (63, "2a0f61f21da6d826c56569acfcb4e08e61811b9ea915163c4122093530da856b")


def test_no_existing_candidate_changed():
    n, d = _old_registry_digest()
    assert (n, d) == OLD_REGISTRY_PINNED


# ------------------------------------------------------------------------------------------------ constants, row
def test_the_constants():
    assert e51.ARMS == ("fgw_shape0x_series_ov", "fgw_shape0x_ov_hi")
    assert e51.CONTROL == {e51.SO: "fgw_shape0x_series", e51.OHI: "fgw_shape0x_ov"}
    assert e51.BEST_OF == ("fgw", "fgw_shape0x_series", "fgw_shape0x_ov")
    assert (e51.BAR, e51.AHEAD, e51.EXPR_TOL, e51.POS_TOL, e51.GATE_TOL) == (0.3, 8, 1e-9, 1e-6, 0.01)
    assert e51.TEST_LIKE == ("half", "late", "new1", "new2") and e51.BUILDS == ("main", "prefix")
    assert e51.POS_METRICS == ("d2_shape", "occupancy_dice", "scale_log_ratio")
    assert e51.EXPR_METRICS == ("de_score", "de_direction", "mmd_u", "variogram")
    assert set(e51.GROUP_NAMES) == {"expression_change", "cell_state", "shape_growth", "local_organisation"}
    assert e51.HI == "fgw_shape0_hi" and e51.E49_RUN_COMMIT == "d76ccd5"
    assert e51.RUN_NAMES == ("fgw", "fgw_shape0x_series", "fgw_shape0x_ov", "fgw_shape0_hi", "fgw_shape0x_series_ov",
                             "fgw_shape0x_ov_hi")
    assert e51.HI_REF == ("fgw_shape0_hi",) and e51.CONTROLS == ("fgw_shape0x_series", "fgw_shape0x_ov")
    assert e51.REPORTED == ("fgw",)


def test_the_run_includes_fgw_shape0_hi_uncounted():
    assert e51.HI in e51.RUN_NAMES and e51.HI not in e51.ARMS and e51.HI not in e51.GATE_RUN
    assert e51.HI not in e51.CONTROLS and e51.HI not in e51.REPORTED
    assert e51.HI_TIE_TOL == e51.GATE_TOL == 0.01
    # the run step refuses a task where an arm does not apply: fgw_shape0_hi applies wherever O_hi does (E24, E28)
    hi, ohi = harness.CANDIDATES[e51.HI], harness.CANDIDATES[e51.OHI]
    assert hi.requires is ohi.requires and hi.tasks == ohi.tasks and hi.derived == ohi.derived == ()
    assert all(hi.applies(t) for t in TASKS.values())
    assert [c.name for c in e51.candidates({})] == [harness.FLOOR, *e51.RUN_NAMES]


def test_has_row(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text("| E50 (one try) | embryo |\n")
    assert not e51.has_row(p)
    p.write_text("| E51 (one row) | embryo |\n")
    assert e51.has_row(p)


def test_the_script_prints_only_step_lines():
    src = (REPO / e51.SCRIPT).read_text()
    prints = re.findall(r"print\((f?\"[^\"]*\")", src)
    assert prints == ['f"gate {a.task} {\'pass\' if rec[\'pass\'] else \'fail\'}"', 'f"{a.task} seed {s} done"',
                      'f"wrote {folder}"', 'f"wrote {files[1]}"']


# ------------------------------------------------------------------------------------------------ pass clauses
def test_both_pass_and_both_go_to_the_captain():
    res = _read(_recs())
    assert not res["void"]
    so, hi = res["arms"][e51.SO], res["arms"][e51.OHI]
    assert so["passes"] and hi["passes"] and not so["void"] and not hi["void"]
    assert so["merged_board"]["mean"] == pytest.approx(2.0) and hi["merged_board"]["mean"] == pytest.approx(1.0)
    assert "both go to the captain" in res["outcome"] and "PASSES" in e51.text(res)


def _so_gain(gain_of):
    return lambda t, s: {"scores": {e51.SO: SCORES[e51.SERIES] + gain_of(t, s)}}


def _ohi_gain(gain_of):
    return lambda t, s: {"scores": {e51.OHI: SCORES[e51.OV] + gain_of(t, s)}}


def test_so_a_needs_the_test_like_board_above_the_bar():
    res = _read(_recs(_so_gain(lambda t, s: 0.3 if t in e51.TEST_LIKE else 2.0)))
    so = res["arms"][e51.SO]
    assert so["clauses"] == {"a": False, "b": True, "c": True} and not so["passes"]
    assert so["test_like_board"]["mean"] == pytest.approx(0.3)  # at the bar is not above it
    assert "SO: does not pass; no file" in res["outcome"] and "both" not in res["outcome"]
    res = _read(_recs(_so_gain(lambda t, s: 0.31 if t in e51.TEST_LIKE else 2.0)))
    assert res["arms"][e51.SO]["clauses"]["a"]


@pytest.mark.parametrize("ahead, want", [(8, True), (7, False)])
def test_so_a_and_b_need_eight_seeds_ahead(ahead, want):
    res = _read(_recs(_so_gain(lambda t, s: 1.0 if s < ahead else -0.1)))
    so = res["arms"][e51.SO]
    assert so["test_like_board"]["ahead"] == so["merged_board"]["ahead"] == ahead
    assert so["clauses"]["a"] is want and so["clauses"]["b"] is want


def test_so_b_needs_the_five_truth_board_above_zero():
    res = _read(_recs(_so_gain(lambda t, s: 2.0 if t in e51.TEST_LIKE else -5.0)))
    so = res["arms"][e51.SO]
    assert so["test_like_board"]["mean"] == pytest.approx(2.0) and so["merged_board"]["mean"] < 0
    assert so["clauses"] == {"a": True, "b": False, "c": True}


@pytest.mark.parametrize("group", ["expression_change", "cell_state", "shape_growth", "local_organisation"])
def test_c_no_group_behind_by_e49s_p4(group):
    def fn(t, s, n_behind=8):
        return {"groups": {e51.SO: {group: 0.49 if s < n_behind else 0.5},
                           e51.OHI: {group: 0.49 if s < n_behind else 0.5}}}

    res = _read(_recs(fn))
    for a in e51.ARMS:
        r = res["arms"][a]
        assert not r["clauses"]["c"] and r["clauses"]["a"] and r["clauses"]["b"] and not r["passes"]
        assert r["groups"][group]["is_behind"] and r["groups"][group]["behind"] == 8
        listed = r["reported_deciding_nothing"]["behind_on_e49_p4_groups"]
        assert listed == ([group] if group in e51.P4_GROUPS_E49 else [])
    res = _read(_recs(lambda t, s: fn(t, s, n_behind=7)))  # seven seeds behind: not behind
    assert all(res["arms"][a]["clauses"]["c"] for a in e51.ARMS)


def test_ohi_a_needs_the_merged_board_above_the_bar_and_b_the_test_like_above_zero():
    res = _read(_recs(_ohi_gain(lambda t, s: 0.3)))
    hi = res["arms"][e51.OHI]
    assert hi["clauses"] == {"a": False, "b": True, "c": True}
    res = _read(_recs(_ohi_gain(lambda t, s: 0.31)))
    assert res["arms"][e51.OHI]["passes"]
    res = _read(_recs(_ohi_gain(lambda t, s: -1.0 if t in e51.TEST_LIKE else 3.0)))
    hi = res["arms"][e51.OHI]
    assert hi["merged_board"]["mean"] > 0.3 and hi["test_like_board"]["mean"] == pytest.approx(-1.0)
    assert hi["clauses"] == {"a": True, "b": False, "c": True}
    res = _read(_recs(_ohi_gain(lambda t, s: 1.0 if s < 7 else -0.1)))
    assert res["arms"][e51.OHI]["clauses"] == {"a": False, "b": True, "c": True}  # 7 of 10 ahead


def test_so_reads_against_the_series_twin_and_ohi_against_the_overlay():
    # SO equal to the overlay and far above the series twin passes; O_hi equal to the series twin fails
    res = _read(_recs(lambda t, s: {"scores": {e51.SO: SCORES[e51.OV], e51.OHI: SCORES[e51.SERIES]}}))
    assert res["arms"][e51.SO]["passes"] and not res["arms"][e51.OHI]["passes"]


# ------------------------------------------------------------------------------------------------ reported tables
def test_the_reported_tables_by_hand():
    gain = {"half": 2.0, "late": 1.0, "new1": -1.0, "new2": 0.5, "e75": -3.0, "e725": 1.0, "x1": 0.0, "x2": 0.5,
            "x3": -0.5}
    res = _read(_recs(lambda t, s: {"scores": {e51.SO: 52.0 + gain[t], e51.FGW: 55.0 if t == "e75" else 51.0}}))
    so = res["arms"][e51.SO]
    rep = so["reported_deciding_nothing"]
    truths = res["truths"]
    want = {tr: float(np.mean([gain[n] for n in ns])) for tr, ns in truths.items()}
    assert rep["against"][e51.SERIES]["per_truth"] == pytest.approx(want)
    assert rep["against"][e51.OV]["per_truth"] == pytest.approx({k: v - 1.0 for k, v in want.items()})
    assert so["merged_board"]["mean"] == pytest.approx(np.mean(list(want.values())))
    inc = rep["increment_over_best"]
    for tr, ns in truths.items():
        best = 55.0 if ns == ["e75"] else 53.0  # fgw on e75's truth, else the overlay
        assert inc[tr]["increment"] == pytest.approx(52.0 + want[tr] - best)
        assert inc[tr]["best"] == (e51.FGW if ns == ["e75"] else e51.OV)
    deep = np.mean([gain[t] for t in e51.nt.DEEP])
    other = np.mean([gain[t] for t in e51.nt.OTHER])
    assert rep["deep_other"] == pytest.approx({"DEEP": deep, "OTHER": other})
    assert rep["per_task"]["half"][e51.NBHD] == pytest.approx(0.6)
    assert rep["per_task"]["half"][f"{e51.NBHD}_minus_{e51.SERIES}"] == pytest.approx(0.0)
    assert rep["against"][e51.SERIES]["prefix_board"]["mean"] == pytest.approx(so["merged_board"]["mean"])
    hi = res["arms"][e51.OHI]["reported_deciding_nothing"]
    assert set(hi["increment_over_best"]) == set(truths) and set(hi["deep_other"]) == {"DEEP", "OTHER"}


def test_the_test_like_merge_is_three_truths():
    res = _read(_recs())
    assert {k: sorted(v) for k, v in res["test_like_truths"].items()} == {
        "T2/E8.0.h5ad": ["half", "late"], "T2/E8.25_late.h5ad": ["new1"], "T2/E8.5.h5ad": ["new2"]}
    assert len(res["truths"]) == 5


# ------------------------------------------------------------------------------------------------ void
def _skill_off(build, name, metrics, by):
    return lambda t, s: ({"skills": {(build, name): {m: _skills(name)[m] + by for m in metrics}}}
                         if (t, s) == ("x2", 4) else {})


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_so_void_on_expression_off_the_series_twin_or_positions_off_the_overlay(build):
    res = _read(_recs(_skill_off(build, e51.SO, ["mmd_u"], 2e-9)))
    assert not res["void"] and res["arms"][e51.SO]["void"] and not res["arms"][e51.OHI]["void"]
    assert "expression skills differ" in res["arms"][e51.SO]["why"][0] and "SO: void" in res["outcome"]
    assert "VOID" in e51.text(res) and "O_hi" in res["outcome"]
    assert not _read(_recs(_skill_off(build, e51.SO, ["mmd_u"], 5e-10)))["arms"][e51.SO]["void"]
    for m in e51.POS_METRICS:
        res = _read(_recs(_skill_off(build, e51.SO, [m], 2e-6)))
        assert res["arms"][e51.SO]["void"] and "position skills" in res["arms"][e51.SO]["why"][0]
        assert not _read(_recs(_skill_off(build, e51.SO, [m], 5e-7)))["arms"][e51.SO]["void"]


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_ohi_void_on_positions_off_the_overlay(build):
    for m in e51.POS_METRICS:
        res = _read(_recs(_skill_off(build, e51.OHI, [m], 2e-6)))
        assert res["arms"][e51.OHI]["void"] and not res["arms"][e51.SO]["void"]
        assert not _read(_recs(_skill_off(build, e51.OHI, [m], 5e-7)))["arms"][e51.OHI]["void"]


def _only_ohi_void(res):
    assert not res["void"] and res["why"] == [] and res["arms"][e51.OHI]["void"]
    assert not res["arms"][e51.SO]["void"] and res["arms"][e51.SO]["passes"]


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_ohi_void_on_expression_off_the_same_run_fgw_shape0_hi(build):
    assert not _read(_recs())["arms"][e51.OHI]["void"]  # O_hi's skills equal the same-run fgw_shape0_hi's
    for m in e51.EXPR_METRICS:
        res = _read(_recs(_skill_off(build, e51.OHI, [m], 2e-9)))
        _only_ohi_void(res)
        assert res["arms"][e51.OHI]["why"] == [
            f"x2 seed 4: fgw_shape0x_ov_hi's expression skills differ from the same-run fgw_shape0_hi's ({build})"]
        assert not _read(_recs(_skill_off(build, e51.OHI, [m], 5e-10)))["arms"][e51.OHI]["void"]
    # the same-run reference moving voids O_hi just the same
    _only_ohi_void(_read(_recs(_skill_off(build, e51.HI, ["variogram"], 2e-9))))


def test_ohi_not_void_when_only_the_committed_expression_differs():
    res = _read(_recs(), hi=_hi(value=HI_EXPR + 1e-3))
    hi = res["arms"][e51.OHI]
    assert not res["void"] and not hi["void"] and hi["passes"]
    rep = hi["reported_deciding_nothing"]
    assert rep["expr_vs_committed_hi"]["overall"] == pytest.approx(1e-3)
    assert rep["expr_vs_committed_hi"]["per_task"] == pytest.approx({t: 1e-3 for t in TASKS})
    assert rep["hi_board_vs_committed"] == {t: {"seed": 0.0, "mean": 0.0} for t in TASKS}
    assert "expr_vs_committed_hi" not in res["arms"][e51.SO]["reported_deciding_nothing"]
    assert "committed fgw_shape0_hi" in e51.text(res)
    # the comparison with the committed skills is under main only: the prefix build's skills do not enter it
    res = _read(_recs(lambda t, s: {"skills": {("prefix", n): {m: 0.9 for m in e51.EXPR_METRICS}
                                               for n in (e51.OHI, e51.HI)}}))
    assert not res["arms"][e51.OHI]["void"]
    assert res["arms"][e51.OHI]["reported_deciding_nothing"]["expr_vs_committed_hi"]["overall"] == 0.0


def _hi_board(by_seed):
    """The run's fgw_shape0_hi board (both builds) moved by by_seed(seed) from the committed score."""
    return lambda t, s: {"scores": {e51.HI: SCORES[e51.HI] + by_seed(s)}}


def test_the_hi_board_tie():
    for d in (0.009, -0.009):
        res = _read(_recs(_hi_board(lambda s, d=d: d)))
        assert not res["arms"][e51.OHI]["void"]
        assert res["arms"][e51.OHI]["reported_deciding_nothing"]["hi_board_vs_committed"]["x3"] == pytest.approx(
            {"seed": 0.009, "mean": 0.009})
    # one seed 0.011 off, the ten-seed mean within 0.01
    res = _read(_recs(_hi_board(lambda s: 0.011 if s == 3 else 0.0)))
    _only_ohi_void(res)
    assert len(res["arms"][e51.OHI]["why"]) == 9 and all("in a seed" in w for w in res["arms"][e51.OHI]["why"])
    # every seed 0.011 off the same way: void on the mean as well
    res = _read(_recs(_hi_board(lambda s: -0.011)))
    _only_ohi_void(res)
    assert "half: the same-run fgw_shape0_hi is 0.0110 off its committed score in the ten-seed mean" in res["arms"][
        e51.OHI]["why"]
    # one task directly (a mean off by more than 0.01 always has a seed off by more than 0.01, so both are named)
    rs = [_rec("x1", s) for s in e51.SEEDS]
    assert e51.hi_tie_problems("x1", "k", rs, {"k": {"board": np.full(10, SCORES[e51.HI] - 0.0105)}}) == [
        "x1: the same-run fgw_shape0_hi is 0.0105 off its committed score in a seed",
        "x1: the same-run fgw_shape0_hi is 0.0105 off its committed score in the ten-seed mean"]
    off = np.full(10, SCORES[e51.HI]) + np.array([0.0] * 9 + [0.0099])
    assert e51.hi_tie_problems("x1", "k", rs, {"k": {"board": off}}) == []
    # a committed task missing, its board missing, of another length, or the committed files refused
    hi = _hi()
    hi.pop(TASKS["x3"].key)
    res = _read(_recs(), hi=hi)
    _only_ohi_void(res)
    assert res["arms"][e51.OHI]["why"] == ["x3: no committed fgw_shape0_hi"]
    hi = _hi()
    hi[TASKS["half"].key].pop("board")
    _only_ohi_void(_read(_recs(), hi=hi))
    hi = _hi()
    hi[TASKS["e75"].key]["board"] = np.full(9, SCORES[e51.HI])
    res = _read(_recs(), hi=hi)
    _only_ohi_void(res)
    assert res["arms"][e51.OHI]["why"] == ["e75: the committed fgw_shape0_hi is of another length"]
    _only_ohi_void(_read(_recs(), why_hi=["a dirty committed file"]))


def test_the_hi_board_tie_reads_build_main_only():
    recs = _recs()
    for rs in recs.values():
        for r in rs:
            r["builds"]["prefix"]["rows"][0]["scores"][e51.HI] = 60.0
    assert not _read(recs)["arms"][e51.OHI]["void"]


@pytest.mark.parametrize("build", ["main", "prefix"])
def test_a_missing_fgw_shape0_hi_voids_ohi_only(build):
    recs = _recs()
    row = recs["new1"][6]["builds"][build]["rows"][0]
    for k in ("scores", "groups", "skills"):
        row[k].pop(e51.HI)
    res = _read(recs)
    _only_ohi_void(res)
    assert res["arms"][e51.OHI]["why"] == [f"new1 seed 6: build {build} lacks fgw_shape0_hi"]
    assert "O_hi: void" in res["outcome"] and "SO: passes" in res["outcome"]


def test_a_missing_arm_voids_that_arm_only():
    res = _read(_recs(lambda t, s: {"drop": (e51.OHI,)} if (t, s) == ("late", 9) else {}))
    assert res["arms"][e51.OHI]["void"] and not res["arms"][e51.SO]["void"]
    res = _read(_recs(lambda t, s: {"drop": (e51.SO,)} if (t, s) == ("late", 9) else {}))
    assert res["arms"][e51.SO]["void"] and not res["arms"][e51.OHI]["void"]


@pytest.mark.parametrize("breaker, needle", [
    (lambda r: r["x3"].pop(), "seeds are not"),
    (lambda r: r.pop("e725"), "e725: seeds are not"),
    (lambda r: r["x3"][0].update(run_commit="d" * 40), "another run commit"),
    (lambda r: r["x3"][0]["builds"].pop("prefix"), "build prefix missing"),
    (lambda r: r["x3"][0]["builds"]["main"]["provenance"].update(tree_clean=False), "unclean"),
    (lambda r: r["x3"][0]["builds"]["main"]["provenance"].update(commit="e" * 40), "another commit"),
    (lambda r: r["x3"][0]["builds"]["main"].update(protocol="emulator"), "bake protocol"),
    (lambda r: r["x3"][0]["builds"]["main"]["board"].update(split="test"), "a test board"),
    (lambda r: r["x3"][0]["facts"].update(board_split_prefix="test_interp"), "a test board"),
    (lambda r: r["x3"][0]["builds"]["main"]["rows"][0]["scores"].pop(e51.OV), "lacks fgw_shape0x_ov"),
    (lambda r: r["x3"][0]["builds"]["main"]["rows"][0]["skills"].pop(e51.SERIES), "lacks fgw_shape0x_series"),
    (lambda r: r["x3"][0]["builds"]["prefix"]["rows"][0]["groups"].pop(e51.FGW), "lacks fgw"),
])
def test_void_for_both(breaker, needle):
    recs = _recs()
    breaker(recs)
    res = _read(recs)
    assert res["void"] and any(needle in w for w in res["why"]) and "VOID" in e51.text(res)


def test_a_gate_or_commit_problem_voids_both():
    res = _read(_recs(), why=["the gate did not pass on x1 at the run commit"])
    assert res["void"] and res["why"] == ["the gate did not pass on x1 at the run commit"]


# ------------------------------------------------------------------------------------------------ gate
def _gate(series=0.0, seed=0.0, ov=0.0, floor=50.0, c0=0.0):
    old = {e51.SERIES: np.full(10, 55.0), e51.C0: np.full(10, 54.0)}
    new = {harness.FLOOR: np.full(10, floor), e51.C0: np.full(10, 54.0 + c0), e51.SERIES: np.full(10, 55.0 + series),
           e51.OV: np.full(10, 54.0 + c0 + 2.0 + ov)}
    new[e51.SERIES][3] += seed
    return e51.gate_check("x1", "key", new, old, [], 2.0, [])


def test_gate_branches():
    assert _gate()["pass"]
    assert _gate(series=0.009, ov=-0.009)["pass"]
    assert not _gate(series=0.011)["pass"]
    assert not _gate(seed=0.02)["pass"]  # within 0.01 in the mean (0.002), off by 0.02 in one seed
    assert _gate(seed=0.0099)["pass"]
    assert not _gate(ov=0.011)["pass"] and not _gate(ov=-0.011)["pass"]
    g = _gate(floor=49.0, c0=3.0)  # the floor and fgw_shape0x are reported, deciding nothing
    assert g["pass"] and g["reported_deciding_nothing"]["floor_off_50"] == pytest.approx(1.0)
    assert g["reported_deciding_nothing"]["c0_vs_e47"]["mean"] == pytest.approx(3.0)
    assert not e51.gate_check("x1", "key", {harness.FLOOR: np.full(10, 50.0)}, {e51.SERIES: np.full(10, 1.0)}, [],
                              2.0, [])["pass"]
    assert not e51.gate_check("x1", "key", {}, None, ["no committed E47 file"], 2.0, [])["pass"]
    assert not e51.gate_check("x1", "key", {}, None, [], None, ["no committed E49 reading"])["pass"]


def test_the_e49_figure_and_its_refusals(tmp_path):
    p = tmp_path / "reading.json"
    good = {"void": False, "commit": "d76ccd5" + "0" * 33, "tasks": {"half": {e51.OV: {"board": {"mean": 1.5}}}}}
    p.write_text(json.dumps(good))
    assert e51.e49_ov_mean("half", p) == (1.5, [])
    v, why = e51.e49_ov_mean("x1", p)
    assert v is None and why == ["x1: no committed E49 figure for fgw_shape0x_ov"]
    p.write_text(json.dumps({**good, "void": True}))
    assert e51.e49_ov_mean("half", p)[1]
    p.write_text(json.dumps({**good, "commit": "f" * 40}))
    assert e51.e49_ov_mean("half", p)[1]
    assert e51.e49_ov_mean("half", tmp_path / "none.json") == (None, ["no committed E49 reading"])


def _task_file(d: Path, key: str, value: float, clean=True, seeds=e51.SEEDS, arm=e51.HI):
    d.mkdir(parents=True, exist_ok=True)
    rows = [{"seed": s, "skills": {arm: {m: value + s for m in METRICS}}, "scores": {arm: 50.0 + value + s}}
            for s in reversed(seeds)]  # stored out of seed order: read back in seed order
    (d / f"T2e_{key.replace('/', '_')}.json").write_text(json.dumps(
        {"task": {"key": key}, "build": "main", "protocol": "bake", "provenance": {"tree_clean": clean},
         "rows": rows}))


def test_committed_hi_reads_main_bake_runs_and_refuses_others(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    _task_file(a, "k1", 0.1)
    _task_file(a, "k0", 0.1, arm="fgw")  # a file without the arm is passed over
    _task_file(b, "k2", 0.2)
    hi, why = e51.committed_hi(((a, "x"), (b, "y")))
    assert why == [] and set(hi) == {"k1", "k2"}
    assert hi["k2"]["mmd_u"] == pytest.approx(0.2 + np.arange(10))
    assert hi["k2"]["board"] == pytest.approx(50.2 + np.arange(10)) and set(hi["k1"]) == {*e51.EXPR_METRICS, "board"}
    _task_file(b, "k3", 0.3, clean=False)
    _task_file(tmp_path / "c", "k1", 0.1)
    hi, why = e51.committed_hi(((a, "x"), (b, "y"), (tmp_path / "c", "z")))
    assert "k3" not in hi and any("not a clean" in w for w in why) and any("committed twice" in w for w in why)


def test_the_committed_figures_cover_the_nine_tasks():
    """Keys and provenance only: no committed value is compared with anything here."""
    hi, why = e51.committed_hi()
    assert why == [] and set(hi) == {t.key for t in TASKS.values()}
    assert all(set(v) == {*e51.EXPR_METRICS, "board"} and len(v["board"]) == 10 for v in hi.values())
    r = json.loads((e51.E49_DIR / "reading.json").read_text())
    assert r["void"] is False and r["commit"].startswith(e51.E49_RUN_COMMIT)
    assert all(e51.OV in r["tasks"][n] and "mean" in r["tasks"][n][e51.OV]["board"] for n in TASKS)


# ------------------------------------------------------------------------------------------------ the launcher
SH = REPO / "experiments" / "coords" / "run_e51.sh"


def _env(tmp_path):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(tmp_path), GIT_CEILING_DIRECTORIES=str(tmp_path.parent))
    return env


def _launch(tmp_path, *args, x=False):
    return subprocess.run(["bash", *(["-x"] if x else []), str(SH), *args], cwd=tmp_path, env=_env(tmp_path),
                          capture_output=True, text=True, timeout=60)


def test_the_launcher_refuses_out_of_order_and_writes_an_exit_file(tmp_path):
    S = tmp_path / "v" / "scratch" / "e51"
    r = _launch(tmp_path)
    assert r.returncode == 2 and "usage" in r.stderr and (S / "usage.exit").read_text().strip() == "2"
    r = _launch(tmp_path, "gate", C)
    assert r.returncode == 96 and "step tests has not exited 0" in r.stderr
    assert (S / "gate.exit").read_text().strip() == "96"
    (S / "gate.exit").write_text("0\n")
    (S / "gate.head").write_text("e" * 40 + "\n")  # passed, but at another commit
    r = _launch(tmp_path, "run", C)
    assert r.returncode == 96 and (S / "run.exit").read_text().strip() == "96"
    r = _launch(tmp_path, "bogus", C)
    assert r.returncode == 97 and (S / "bogus.exit").read_text().strip() == "97"
    trace = _launch(tmp_path, "read", C, x=True)  # a harmless path traced: it stops before git or Python
    assert trace.returncode == 96 and "finish 96" in trace.stderr and "uv run" not in trace.stderr


def test_the_launcher_refuses_another_head_outside_a_repository(tmp_path):
    S = tmp_path / "v" / "scratch" / "e51"
    r = _launch(tmp_path, "tests", C)
    assert r.returncode == 99 and (S / "tests.exit").read_text().strip() == "99" and not (S / "tests.head").exists()
    assert not (tmp_path / "runs").exists()


def test_the_launcher_never_overwrites_a_reading(tmp_path):
    S = tmp_path / "v" / "scratch" / "e51"
    S.mkdir(parents=True)
    (S / "run.exit").write_text("0\n")
    (S / "run.head").write_text(C + "\n")
    out = tmp_path / "experiments" / "coords" / "results" / "e51_so_ohi"
    out.mkdir(parents=True)
    (out / "reading.txt").write_text("kept")
    r = _launch(tmp_path, "read", C)
    assert r.returncode == 95 and "never overwritten" in r.stderr and (S / "read.exit").read_text().strip() == "95"
    assert (out / "reading.txt").read_text() == "kept"
    (out / "reading.txt").rename(out / "kept.txt")  # without a reading the step goes on to the tree checks
    r = _launch(tmp_path, "read", C)
    assert r.returncode == 99 and (S / "read.exit").read_text().strip() == "99"


def test_the_launcher_deletes_nothing_parses_and_lists_existing_tests():
    text = SH.read_text()
    code = [line.split("#")[0] for line in text.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    assert ":?" not in text
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", text)
    assert "tests/test_e51_so_ohi_read.py" in listed and all((REPO / f).exists() for f in listed)
    assert "SCRIPT=experiments/coords/e51_so_ohi_read.py" in text and "R=runs/e51" in text
    assert "TASKS=(half late new1 new2 e75 e725 x1 x2 x3)" in text
    assert set(re.search(r"TASKS=\(([^)]+)\)", text).group(1).split()) == set(TASKS)
