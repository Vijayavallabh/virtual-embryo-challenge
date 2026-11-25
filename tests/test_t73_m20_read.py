"""Row T73 (experiments/harness/TRIES.md; gate G0 of row R58's file T1_s2rcMSBm20): experiments/t1_m1/m20_read.py on
synthetic stages and hand-made task files only (the fixtures of tests/test_t65_ms_read.py). No released file,
statistics file, board file or score is read. The script reaches the builder through ms_read and recomputed_read."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_scope import HELD
from test_t1_state_transport import _gmom_ctx
from test_t65_ms_read import HELD_M1, SHARES, WITHIN, _arm, _classified, _directions, _fb

from vec import t1_scope
from vec.io import dense

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "t1_m1"))
import m20_read as m20  # noqa: E402

C = "c" * 40
GATE, MS025, MS05, MS1 = m20.GATE, m20.MS025, m20.MS05, m20.MS1


def _x(stage) -> bytes:
    return np.ascontiguousarray(dense(stage.X), dtype=np.float32).tobytes()


def _spy(monkeypatch) -> dict:
    seen = {"damping": [], "fraction": []}
    inner_m, inner_e = t1_scope.momentum_shift, t1_scope.external_shift

    def mom(start, ref, target_like, damping, *a, **k):
        seen["damping"].append(float(damping))
        return inner_m(start, ref, target_like, damping, *a, **k)

    def ext(target_like, directions, fraction, *a, **k):
        seen["fraction"].append(float(fraction))
        return inner_e(target_like, directions, fraction, *a, **k)

    monkeypatch.setattr(t1_scope, "momentum_shift", mom)
    monkeypatch.setattr(t1_scope, "external_shift", ext)
    return seen


# ------------------------------------------------------------------ the arms


def test_at_damping_0_5_the_population_is_ms_reads():
    e95, e85 = _classified()
    fb = _fb()
    dirs = _directions(e95.n_vars)
    base0, _ = t1_scope.reweight(e95, SHARES, 300, 4, within=WITHIN, fallback=fb["E9.5"])
    groups = m20.rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), dirs)[0]
    kw = {"within": WITHIN, "start_fallback": fb["E8.5"], "ref_fallback": fb["E9.5"], "min_cells": 20, "eps": 0.5}
    a, ra = m20.ms_population_at(base0, e85, e95, dirs, 0.3, groups, HELD_M1, 4, kw, 0.5)
    b, rb = m20.ms.ms_population(base0, e85, e95, dirs, 0.3, groups, HELD_M1, 4, kw)
    assert _x(a) == _x(b) and ra.as_dict() == rb.as_dict()
    c, _ = m20.ms_population_at(base0, e85, e95, dirs, 0.3, groups, HELD_M1, 4, kw, 1.0)
    assert _x(c) != _x(a)


def test_on_the_library_route_the_gate_is_t65s_arm_and_each_arm_takes_its_settings(monkeypatch):
    e95, e85 = _classified()
    factors = m20.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    dirs, fb = _directions(e95.n_vars), _fb()
    t65, _ = m20.ms.library_arms(e85, e95, SHARES, factors, HELD, HELD_M1, 4, WITHIN, fb, dirs, 0.2, n_cells=300)
    seen = _spy(monkeypatch)
    arms, recs = m20.library_arms(e85, e95, SHARES, factors, HELD, HELD_M1, 4, WITHIN, fb, dirs, 0.2, n_cells=300)
    assert tuple(arms) == m20.ARMS and _x(arms[GATE]) == _x(t65[m20.ms.MOM])
    assert seen["damping"] == [0.5, 0.25, 0.5, 1.0]
    assert seen["fraction"] == pytest.approx([0.2, 0.3, 0.3, 0.3], abs=1e-15)
    assert len({_x(v) for v in arms.values()}) == 4
    for arm in m20.COUNTED:
        assert recs[arm]["non_covered_equal"] is True and recs[arm]["covered_cells_changed"] > 0
        assert recs[arm]["covered_set"] == recs[arm]["builder_covered_set"]


def test_on_the_cut_route_the_gate_is_t65s_arm(monkeypatch):
    ref, start, _ = _small_inputs()
    fits = m20.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    t65, _ = m20.ms.cut_arms(_gmom_ctx(ref, start, seed=1000002), fits["file"], dirs, 0.2, HELD_M1)
    seen = _spy(monkeypatch)
    arms, recs = m20.cut_arms(_gmom_ctx(ref, start, seed=1000002), fits["file"], dirs, 0.2, HELD_M1)
    assert _x(arms[GATE]) == _x(t65[m20.ms.MOM])
    assert seen["damping"] == [0.5, 0.25, 0.5, 1.0]
    assert seen["fraction"] == pytest.approx([0.2, 0.3, 0.3, 0.3], abs=1e-15)
    assert all(recs[a]["non_covered_equal"] for a in m20.COUNTED)


def test_a_counted_arm_that_moves_an_uncovered_cell_stops_the_run():
    e95, _ = _classified()
    dirs = _directions(e95.n_vars)
    groups = m20.rc.direction_kit(e95.obs[t1_scope.LABEL_KEY].astype(str).unique(), dirs)[0]
    pop = e95.copy()
    x = np.asarray(dense(pop.X), dtype=np.float32)
    y = x.copy()
    cov = m20.ms.covered_mask(pop.obs[t1_scope.LABEL_KEY].astype(str).to_numpy(), groups, dirs)
    y[np.flatnonzero(~cov)[0], 0] += 1.0
    with pytest.raises(SystemExit, match="a cell outside the covered set differs"):
        m20.covered_record(pop, groups, dirs, x, y, MS1, None)


# ------------------------------------------------------------------ the reading


def _relabel(**change):
    return {"sha256": m20.rc.STATS_SHA256, "stats_sha256_hashed": m20.rc.STATS_SHA256,
            "yaml_sha256": m20.rc.YAML_SHA256, "route": "recomputed", "record_route": "recomputed",
            "authors_labels_read": False, **change}


def _covered():
    return {"covered_set": ["V-CM"], "builder_covered_set": ["V-CM"], "covered_cells": 30, "covered_cells_changed": 30,
            "non_covered_equal": True, "shift": []}


SETTINGS = {GATE: {"damping": 0.5, "fraction": 0.2}, MS025: {"damping": 0.25, "fraction": 1.5 * 0.2},
            MS05: {"damping": 0.5, "fraction": 1.5 * 0.2}, MS1: {"damping": 1.0, "fraction": 1.5 * 0.2}}


def _task(ms1: float | list[float], ms05: float = 0.3, ms025: float = 0.0, lib=False, base=50.0, commit=C,
          clean=True, seeds=range(10), inverted=False, same=(), **top):
    rows = []
    for s in seeds:
        g1 = ms1[s] if isinstance(ms1, list) else ms1
        arms = {GATE: _arm(base + s / 100), MS025: _arm(base + s / 100 + ms025), MS05: _arm(base + s / 100 + ms05),
                MS1: _arm(base + s / 100 + g1)}
        sha = {GATE: "a" * 64, MS025: "b" * 64, MS05: "d" * 64, MS1: "e" * 64}
        for x, y in same:
            sha[y] = sha[x]
        row = {"seed": s, "arms": arms, "sha256": sha, "covered": {a: _covered() for a in m20.COUNTED}}
        if lib:
            arms["copy_last"], sha["copy_last"] = _arm(40.0), "f" * 64
            row["inverted_anchors"] = ["x"] if inverted and s == 3 else []
        rows.append(row)
    out = {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "fraction": 0.2,
           "settings": json.loads(json.dumps(SETTINGS)), "momentum_stream": 1, "direction_stream": 2,
           "relabel": _relabel()}
    out.update(top)
    return out


def _loaded(chamber=1.0, whole=0.5, libs=(0.2, 0.2, 0.2), **kw):
    out = {"chamber": _task(chamber), "whole": _task(whole)}
    out.update({lib: _task(g, lib=True) for lib, g in zip(("1", "2", "3"), libs, strict=True)})
    out.update(kw)
    return out


def _committed(base=50.0):
    return {k: {s: base + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}


OK = {"reasons": []}


def test_pass():
    r = m20.t73_reading(_loaded(), _committed(), OK, (), C)
    assert r["outcome"] == "PASS" and r["ahead"] and not r["failed"]
    assert r["ms1_minus_ms05"]["chamber"]["mean"] == pytest.approx(0.7)
    assert r["ms1_minus_ms025"]["chamber"]["mean"] == pytest.approx(1.0)
    assert r["valid_libraries"] == ["1", "2", "3"] and r["j6_range"] == pytest.approx([0.0, 0.7])
    lines = m20.outcome_lines(r)
    assert lines[0] == "Outcome: T73 PASS (m20's step is AHEAD)" and "t1_lm_m20_line.sh" in lines[-1]


def _with(key, task):
    d = _loaded()
    d[key] = task
    return d


def test_the_boundaries_that_pass():
    p = {"mean": 0.3, "ahead": 8}
    tasks1 = {"chamber": p, "whole": {"mean": 1e-9}, "1": {"mean": -0.3}, "2": {"mean": -0.3}}
    assert m20.statement(tasks1, {"chamber": p}, ["1", "2"], {}) == (True, [])
    for t1, t2 in (({**tasks1, "chamber": {"mean": 0.3 - 1e-9, "ahead": 8}}, {"chamber": p}),
                   ({**tasks1, "chamber": {"mean": 0.3, "ahead": 7}}, {"chamber": p}),
                   ({**tasks1, "whole": {"mean": 0.0}}, {"chamber": p}),
                   ({**tasks1, "2": {"mean": -0.3 - 1e-9}}, {"chamber": p}),
                   (tasks1, {"chamber": {"mean": 0.3, "ahead": 7}})):
        assert m20.statement(t1, t2, ["1", "2"], {})[0] is False


@pytest.mark.parametrize("loaded, part", [
    (_loaded(chamber=0.59), "on the chamber cut ms1_x15 - ms05_x15 is not +0.3"),
    (_loaded(chamber=[1.0] * 7 + [0.2] * 3), "on the chamber cut ms1_x15 - ms05_x15 is not +0.3"),
    (_loaded(whole=0.3), "on the whole cut the mean of ms1_x15 - ms05_x15 is not above 0"),
    (_loaded(libs=(0.2, -0.01, 0.2)), "a valid library task's mean of ms1_x15 - ms05_x15 is below -0.3 (library 2)"),
    (_with("chamber", _task(1.0, ms025=0.71)), "on the chamber cut ms1_x15 - ms025_x15 is not +0.3"),
    (_with("whole", _task(0.5, same=((MS05, MS1),))), "ms05_x15 and ms1_x15 are one file on every seed of whole"),
    (_with("1", _task(0.2, lib=True, same=((MS025, MS05),))), "ms025_x15 and ms05_x15 are one file on every seed of 1"),
])
def test_stop(loaded, part):
    r = m20.t73_reading(loaded, _committed(), OK, (), C)
    assert r["outcome"] == "STOP" and any(part in f for f in r["failed"]), r["failed"]
    assert m20.outcome_lines(r)[0].startswith("Outcome: T73 STOP (")


def test_an_invalid_library_is_not_read():
    loaded = _loaded()
    loaded["2"] = _task(-2.0, lib=True, inverted=True)
    r = m20.t73_reading(loaded, _committed(), OK, (), C)
    assert r["outcome"] == "PASS" and r["valid_libraries"] == ["1", "3"]


VOIDS = {
    "no task file": lambda d: {**d, "1": None},
    "seeds": lambda d: {**d, "whole": _task(0.5, seeds=range(9))},
    "arms": lambda d: {**d, "chamber": _without_arm(_task(1.0), MS025)},
    "one commit": lambda d: {**d, "whole": _task(0.5, commit="d" * 40)},
    "not clean": lambda d: {**d, "whole": _task(0.5, clean=False)},
    "the gate failed": lambda d: {**d, "chamber": _task(1.0, base=50.021)},
    "damping or fraction": lambda d: {**d, "3": _task(0.2, lib=True, settings={**SETTINGS, MS1: {"damping": 2.0,
                                                                                               "fraction": 0.3}})},
    "practice fraction": lambda d: {**d, "3": _task(0.2, lib=True, fraction=0.4)},
    "a stream": lambda d: {**d, "3": _task(0.2, lib=True, momentum_stream=3)},
    "statistics file": lambda d: {**d, "1": _task(0.2, lib=True, relabel=_relabel(sha256="x"))},
    "settings file": lambda d: {**d, "1": _task(0.2, lib=True, relabel=_relabel(yaml_sha256="x"))},
    "route": lambda d: {**d, "1": _task(0.2, lib=True, relabel=_relabel(route="authors_x"))},
    "authors_labels_read": lambda d: {**d, "1": _task(0.2, lib=True, relabel=_relabel(authors_labels_read=None))},
    "non_covered_equal": lambda d: {**d, "chamber": _with_covered(_task(1.0), non_covered_equal=False)},
    "the covered set": lambda d: {**d, "chamber": _with_covered(_task(1.0), builder_covered_set=["V-CM", "x"])},
}


def _without_arm(t, arm):
    for x in t["rows"]:
        x["arms"].pop(arm)
    return t


def _with_covered(t, **change):
    t["rows"][2]["covered"][MS1].update(change)
    return t


@pytest.mark.parametrize("word", list(VOIDS))
def test_every_void_reason_fires(word):
    loaded = VOIDS[word](_loaded())
    r = m20.t73_reading(loaded, _committed(), OK, (), C)
    assert r["outcome"] == "VOID" and any(word in x for x in r["reasons"]), r["reasons"]
    assert m20.outcome_lines(r)[0].startswith("Outcome: T73 VOID (")


def test_the_gate_passes_within_0_02_and_void_reasons_from_the_read_and_the_source():
    assert m20.t73_reading(_with("chamber", _task(1.0, base=50.0199)), _committed(), OK, (), C)["outcome"] == "PASS"
    r = m20.t73_reading(_loaded(), _committed(), {"reasons": ["x"]}, ("y",), C)
    assert r["outcome"] == "VOID" and "at the read: x" in r["reasons"] and "source: y" in r["reasons"]
    r = m20.t73_reading(_loaded(), _committed(), OK, (), "d" * 40)
    assert r["outcome"] == "VOID" and any("row's commit" in x for x in r["reasons"])


def test_the_counts_and_constants():
    assert m20.TRIES_COUNT == 15 and m20.COUNTED == (MS025, MS05, MS1) and m20.GATE == "m1rc_gmom_mom"
    assert m20.SETTINGS == {GATE: (0.5, 1.0), MS025: (0.25, 1.5), MS05: (0.5, 1.5), MS1: (1.0, 1.5)}
    assert m20.PRACTICE_FRACTION == 0.2 and m20.REPRO_TOL == 0.02 and m20.BAR == 0.3 and m20.LIB_FLOOR == -0.3
    assert m20.ms.source_guard(Path(m20.__file__)) == []


# ------------------------------------------------------------------ J7, refusals and the read


def _task_file(folder: Path, name: str, arm: str, scores: dict[int, float]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / name).write_text(json.dumps({"rows": [{"seed": s, "arms": {arm: {"score": v}}}
                                                    for s, v in scores.items()]}))


def test_j7_compares_committed_t65_and_t62_files(tmp_path):
    for k, f in m20.rc.task_files().items():
        _task_file(tmp_path / "t65", f, "m1rc_gmom", {s: 50.0 + s / 100 for s in range(10)})
        _task_file(tmp_path / "t62", f, "m1rc_gmom", {s: 50.0 + s / 100 + (0.01 if k == "2" else 0.0)
                                                      for s in range(10)})
    res = m20.j7(tmp_path / "t65", tmp_path / "t62")
    assert res["holds"] and res["largest_difference"] == pytest.approx(0.01)
    assert res["tasks"]["2"]["seeds_compared"] == 10
    _task_file(tmp_path / "t62", m20.rc.task_files()["1"], "m1rc_gmom", {s: 50.03 for s in range(9)})
    res = m20.j7(tmp_path / "t65", tmp_path / "t62")
    assert not res["holds"] and res["tasks"]["1"]["seeds_compared"] == 9


def _ledger(tmp_path, first: str | None) -> Path:
    p = tmp_path / "TRIES.md"
    p.write_text("# ledger\n" + ("" if first is None else f"| {first} | T1 | G0 |\n"))
    return p


def test_the_row_must_name_the_script_whole(tmp_path):
    for first, ok in (("T73 (G0: `experiments/t1_m1/m20_read.py` and `experiments/t1_m1/t73_m20.sh`)", True),
                      ("T73 (G0: experiments/t1_m1/m20_read.py)", True), (None, False),
                      ("T73 (G0: `experiments/t1_m1/m20_read.pyc`)", False),
                      ("T73 (G0: `xexperiments/t1_m1/m20_read.py`)", False), ("T73 G0 m20_read.py", False)):
        assert m20.row_names_script(_ledger(tmp_path, first)) is ok, first


def test_the_runs_are_refused_before_anything_loads(tmp_path, monkeypatch):
    monkeypatch.setattr(m20, "TRIES", _ledger(tmp_path, "T73 (G0: `experiments/t1_m1/m20_read.py`)"))
    state = {"commit": C, "clean": True}
    monkeypatch.setattr(m20.gv, "tree_state", lambda: dict(state))
    monkeypatch.setattr(m20.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    out = tmp_path / "tasks"
    assert m20.refuse(out, "lib1.json", C)[0]["commit"] == C
    with pytest.raises(SystemExit, match="HEAD is not the run commit"):
        m20.refuse(out, "lib1.json", "d" * 40)
    state["clean"] = False
    with pytest.raises(SystemExit, match="not clean"):
        m20.refuse(out, "lib1.json", C)
    state["clean"] = True
    out.mkdir()
    (out / "lib1.json").write_text("{}")
    with pytest.raises(SystemExit):
        m20.refuse(out, "lib1.json", C)
    monkeypatch.setattr(m20.rc, "relabel_checks", lambda *a, **k: {"passes": False, "reasons": ["no record"]})
    with pytest.raises(SystemExit, match="no record"):
        m20.refuse(tmp_path / "other", "lib1.json", C)
    monkeypatch.setattr(m20, "TRIES", _ledger(tmp_path, None))
    with pytest.raises(SystemExit, match="no single row T73"):
        m20.refuse(tmp_path / "other", "lib1.json", C)


def test_the_read_prints_the_outcome_words_only_and_never_overwrites(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(m20, "TRIES", _ledger(tmp_path, "T73 (G0: `experiments/t1_m1/m20_read.py`)"))
    monkeypatch.setattr(m20.gv, "tree_state", lambda: {"commit": C, "clean": True, "utc": "now"})
    monkeypatch.setattr(m20.gv, "stamp", lambda started, paths: {"commit": C, "tree_clean": True, "utc": "now"})
    monkeypatch.setattr(m20.rc, "relabel_checks", lambda *a, **k: OK)
    monkeypatch.setattr(m20, "committed_scores", lambda f, arm, tasks: {s: 50.0 + s / 100 for s in range(10)})
    out = tmp_path / "t73_m20"
    for k, f in m20.rc.task_files().items():
        (out / "tasks").mkdir(parents=True, exist_ok=True)
        task = _task(1.0 if k == "chamber" else 0.5 if k == "whole" else 0.2, lib=k in ("1", "2", "3"))
        (out / "tasks" / f).write_text(json.dumps(task))
    res = m20.read(out, C)
    printed = capsys.readouterr().out
    assert res["reading"]["outcome"] == "PASS" and printed.splitlines()[0] == "Outcome: T73 PASS (m20's step is AHEAD)"
    assert "50." not in printed and (out / "t73.txt").exists() and (out / "t73.json").exists()
    with pytest.raises(SystemExit):
        m20.read(out, C)
