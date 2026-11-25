"""Row E57's read script (experiments/coords/e57_rs_read.py) and its launcher (experiments/coords/run_e57.sh), on
hand-made tables, synthetic pair files, synthetic stages and throwaway git repositories only. Nothing is read from data/
or boards/, no target is opened and nothing is scored against a real truth.

Pinned: step 0's three clauses decide on the test pair alone, each failing alone reads STOP, and an incomplete test pair
reads STOP; the gate holds O to E55's committed per-task figure within 0.01; the void conditions are the row's; each of
the clauses (a) to (e) works on hand-made tables, alone and at its boundary, and the idle case reads STOP; the outcome's
first line is the launcher's and the build launcher's word; and the launcher refuses a wrong commit (96), an unclean
tree (98), a changed rules file (97), a missing or malformed row E57 (94) and an existing reading (93) before any python
of the repository runs, then calls its steps in order and stops where the row says."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pytest
from conftest import git, git_env

from vec import evaluate as ev
from vec import harness
from vec import t2_build as tb
from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
REL = "experiments/coords/run_e57.sh"
LAUNCHER = REPO / REL
COMMIT = "c" * 40
O, ARM, PL = "fgw_shape0x_ov", "fgw_shape0x_ov_rs", "fgw_shape0x_ov_rs_same"  # noqa: E741 (the row names O)
TASKS = ("half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3")
TRUTH = {"half": "T2/E8.0.h5ad", "late": "T2/E8.0.h5ad", "new1": "T2/E8.25_late.h5ad", "x1": "T2/E8.25_late.h5ad",
         "new2": "T2/E8.5.h5ad", "x2": "T2/E8.5.h5ad", "x3": "T2/E8.5.h5ad", "e75": "T2/E7.5.h5ad",
         "e725": "T2/E7.25.h5ad"}
E57_ROW = ("| E57 (E57a to E57c; one row, one script: E57a step 0, a counted reading; if it passes, two tries, E57b RS "
           "and E57c its placebo; `experiments/coords/e57_rs_read.py` through `experiments/coords/run_e57.sh`; void "
           "if raw de_score and de_direction more than 1e-4 off O's; synthetic) | embryo | synthetic |")


@pytest.fixture(scope="module")
def e57():
    sys.path.insert(0, str(REPO / "experiments" / "coords"))
    import e57_rs_read

    return e57_rs_read


# ------------------------------------------------------------------------------------------------ hand-made tables
BASE_SKILLS = {"de_score": 0.55, "de_direction": 0.6, "mmd_u": 0.62, "variogram": 0.58, "d2_shape": 0.52,
               "occupancy_dice": 0.5, "scale_log_ratio": 0.7, "neighborhood_mmd": 0.6}


def _scored(skills: dict) -> dict:
    score, groups = ev.aggregate("T2", skills)
    return {"score": score, "groups": groups, "skills": dict(skills)}


def table(arm=None, placebo=None, o_nb=None, noise: float = 1e-4, seeds=range(10), tasks=TASKS, commit=COMMIT,
          facts=None) -> dict:
    """Per task, ten per-seed records as the run step writes them. `arm` and `placebo`: (task, seed) -> {metric: skill
    change against O}; `o_nb`: task -> O's neighborhood_mmd skill. Board and groups follow from the skills by
    vec.evaluate.aggregate, as the scorer's do."""
    rng = np.random.default_rng(0)
    out = {}
    for t in tasks:
        recs = []
        for s in seeds:
            o = dict(BASE_SKILLS, neighborhood_mmd=(o_nb or {}).get(t, BASE_SKILLS["neighborhood_mmd"]))
            o["mmd_u"] += rng.normal(scale=0.01)
            arms = {harness.FLOOR: {m: 0.5 for m in BASE_SKILLS}, O: o}
            for name, change in ((ARM, arm), (PL, placebo)):
                sk = dict(o)
                for m, v in ((change or (lambda t, s: {}))(t, s)).items():
                    sk[m] = sk[m] + v
                sk["mmd_u"] += rng.normal(scale=noise)  # the scorer's fixed draws read other rows
                arms[name] = sk
            sc = {n: _scored(sk) for n, sk in arms.items()}
            row = {"seed": s, "scores": {n: v["score"] for n, v in sc.items()},
                   "skills": {n: v["skills"] for n, v in sc.items()}, "groups": {n: v["groups"] for n, v in sc.items()}}
            builds = {b: {"provenance": {"commit": commit, "tree_clean": True}, "build": b, "protocol": "bake",
                          "rows": [json.loads(json.dumps(row))], "board": {"split": "pseudo"}}
                      for b in ("main", "prefix")}
            f = {"same_coords": {ARM: True, PL: True}, "same_count": {ARM: True, PL: True},
                 "same_sorted_rows": {ARM: True, PL: True}, "rows_moved": {ARM: True, PL: True},
                 "files_identical_across_builds": True,
                 "steps": {"run_sorted": {"contrast_earlier": 0.5, "moved_share": 0.4, "median_homology_error": 0.1,
                                          "own_run_share_o": 0.55, "own_run_share": 0.93, "sd_f_o": 0.15, "sd_f": 0.45},
                           "run_sorted_same": {"contrast_earlier": 0.5, "moved_share": 0.4,
                                               "median_homology_error": 0.08, "placebo_matches_swaps": True,
                                               "own_run_share_o": 0.55, "own_run_share": 0.55}}}
            f.update((facts or (lambda t, s: {}))(t, s))
            recs.append({"task": t, "seed": s, "run_commit": commit, "builds": builds, "facts": f})
        out[t] = recs
    return out


def _nb(v):
    return lambda t, s: {"neighborhood_mmd": v(t, s) if callable(v) else v}


def test_a_clean_gain_passes_every_clause_and_names_the_question(e57):
    res = e57.reading(table(arm=_nb(0.03)), [], COMMIT)
    assert res["void"] is False and res["clauses"] == dict.fromkeys("abcde", True) and not res["idle"]
    assert res["outcome"] == "E57: PASS" and res["passes"] is True
    # 100 x 0.25 x 0.03, give or take the scorer's draw noise on cell_state (mmd_u at 1e-4)
    assert res["merged"]["board"]["mean"] == pytest.approx(0.75, abs=1e-3) and res["merged"]["board"]["ahead"] == 10
    q = res["question"]
    assert q["file"] == "T2embryo_s2x_ovrs.h5ad" and "--neighbourhood run_sorted --captain-ruling <ROW>" in q["command"]
    assert set(q["per_truth"]) == set(TRUTH.values()) and set(q["half_and_x3"]) == {"half", "x3"}
    assert set(q["per_class"]) == {"class 1 (half, late, new1)", "class 2 (x1, x2, x3)", "own profile (new2, e75)",
                                   "e725"}
    text = e57.text(res)
    assert text.split("\n", 1)[0] == "E57: PASS" and "no file without his ruling" in text


def test_the_merge_is_by_truth_at_equal_weight(e57):
    """A gain on half alone: E8.0's truth holds half and late, so the five-truth board moves by gain / 2 / 5."""
    res = e57.reading(table(arm=lambda t, s: {"neighborhood_mmd": 0.08 if t == "half" else 0.0}, noise=0.0), [],
                      COMMIT)
    assert res["merged"]["board"]["mean"] == pytest.approx(25 * 0.08 / 2 / 5, abs=1e-9)
    assert res["per_truth"]["T2/E8.0.h5ad"] == pytest.approx(25 * 0.08 / 2, abs=1e-9)
    assert res["half_and_x3"]["half"]["board"] == pytest.approx(2.0, abs=1e-9)
    assert res["merged"]["board"]["ahead"] == 10 and res["half_and_x3"]["x3"]["board"] == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("case, failing, idle", [
    ("a_small", {"a"}, True),  # +0.29: below the bar, and idle
    ("a_seeds", {"a"}, False),  # +0.5 but only 7 seeds ahead
    ("b_carried", {"b"}, False),  # half the gain on cell_state
    ("b_cell", {"b"}, False),  # cell_state moves by more than 3 seed SEs
    ("c", {"c"}, False),  # the placebo gains as much
    ("d", {"d"}, False),  # new2 and e75 lose 0.5
    ("e", {"e"}, False),  # class 2 loses on local organisation
])
def test_each_clause_fails_alone_and_reads_stop(e57, case, failing, idle):
    lo = 0.03
    arm, placebo = _nb(lo), None
    if case == "a_small":
        arm = _nb(0.29 / 25)
    elif case == "a_seeds":
        arm = _nb(lambda t, s: 0.04 if s < 7 else -0.02)
        placebo = _nb(lambda t, s: 0.0 if s < 7 else -0.05)  # so that the arm still leads its placebo every seed
    elif case == "b_carried":
        arm = lambda t, s: {"neighborhood_mmd": 0.02, "mmd_u": 0.02 / 0.6, "variogram": 0.0}  # noqa: E731
    elif case == "b_cell":
        arm = lambda t, s: {"neighborhood_mmd": lo, "mmd_u": 0.0004}  # noqa: E731
    elif case == "c":
        placebo = _nb(lo)
    elif case == "d":
        arm = _nb(lambda t, s: -0.02 if t in ("new2", "e75") else 0.07)
    elif case == "e":
        arm = _nb(lambda t, s: -0.02 if t in ("x1", "x2", "x3") else 0.06)
    res = e57.reading(table(arm=arm, placebo=placebo), [], COMMIT)
    assert res["void"] is False
    assert {k for k, v in res["clauses"].items() if not v} == failing, (res["clauses"], res["merged"])
    assert res["idle"] is idle and res["outcome"] == "E57: STOP" and res["question"] is None
    assert e57.text(res).split("\n", 1)[0] == "E57: STOP" and e57.STOP_TEXT in e57.text(res)


def test_the_clauses_at_their_boundaries(e57):
    eps = 1e-6
    D = np.full(10, 0.30 + eps)
    base = {"board": D, "carried": D * (0.8 + eps), "cell": np.r_[0.001, -0.001].repeat(5),
            "lo_pl": np.full(10, 0.01 + eps), "own": np.full(10, -0.30 + eps),
            "classes": {k: np.full(10, eps) for k in e57.CLASSES}}
    assert e57.decide(base) == {"clauses": dict.fromkeys("abcde", True), "idle": False, "passes": True}
    seven = np.r_[np.full(7, 0.6), np.full(3, -0.2)]  # mean 0.36, but 7 of 10 seeds ahead
    for change, clause in (({"board": seven, "carried": seven}, "a"),
                           ({"carried": D * (0.8 - eps)}, "b"), ({"lo_pl": np.full(10, 0.01 - eps)}, "c"),
                           ({"lo_pl": np.r_[np.full(7, 0.05), np.full(3, -0.01)]}, "c"),
                           ({"own": np.full(10, -0.30 - eps)}, "d"),
                           ({"classes": {k: np.full(10, 0.0 if "2" in k else eps) for k in e57.CLASSES}}, "e")):
        dec = e57.decide({**base, **change})
        assert [k for k, v in dec["clauses"].items() if not v] == [clause] and not dec["passes"], change
    cell = np.r_[np.full(9, 0.01), 0.011]  # a mean far beyond 3 seed SEs
    assert e57.decide({**base, "cell": cell})["clauses"]["b"] is False
    assert e57.decide({**base, "cell": np.zeros(10)})["clauses"]["b"] is True  # no change at all
    idle = e57.decide({**base, "board": np.full(10, 0.30 - eps)})
    assert idle["idle"] is True and idle["clauses"]["a"] is False and not idle["passes"]
    assert e57.decide({**base, "board": np.full(10, -0.30 - eps)})["idle"] is False  # a loss is not idle


def test_the_gated_shape_part_counts_only_where_os_neighbourhood_skill_is_below_half(e57):
    """x1 with O's neighbourhood skill at 0.4: the gate min(1, skill/0.5) binds, so the arm's neighbourhood gain also
    raises gated shape_growth there, and (b) counts that part as carried by local organisation."""
    t = table(arm=_nb(0.04), o_nb={"x1": 0.4})
    r = t["x1"][0]["builds"]["main"]["rows"][0]
    sg = r["groups"][ARM]["shape_growth"] - r["groups"][O]["shape_growth"]
    assert sg > 0  # raw shape skills are O's; only the gate moved
    car = e57.carried(t["x1"])
    lo = r["groups"][ARM]["local_organisation"] - r["groups"][O]["local_organisation"]
    assert car[0] == pytest.approx(25 * lo + 25 * sg, rel=1e-12)
    car_half = e57.carried(t["half"])  # O at 0.6: the gate does not bind, no shape part
    assert car_half[0] == pytest.approx(25 * 0.04, rel=1e-9)
    res = e57.reading(t, [], COMMIT)
    assert res["clauses"]["b"] is True and res["merged"]["carried"]["mean"] >= 0.8 * res["merged"]["board"]["mean"]


@pytest.mark.parametrize("kind", ["seed", "task", "commit", "build_commit", "unclean", "arm", "coords", "count", "rows",
                                  "d2_shape", "occupancy_dice", "scale_log_ratio", "de_score", "de_direction"])
def test_every_void_condition_voids(e57, kind):
    t = table(arm=_nb(0.03))
    if kind == "seed":
        t["new1"] = t["new1"][:9]
    elif kind == "task":
        del t["x2"]
    elif kind == "commit":
        t["late"][3]["run_commit"] = "d" * 40
    elif kind == "build_commit":
        t["late"][3]["builds"]["prefix"]["provenance"]["commit"] = "d" * 40
    elif kind == "unclean":
        t["late"][3]["builds"]["main"]["provenance"]["tree_clean"] = False
    elif kind == "arm":
        del t["e75"][0]["builds"]["main"]["rows"][0]["scores"][PL]
    elif kind in ("coords", "count", "rows"):
        key = {"coords": "same_coords", "count": "same_count", "rows": "same_sorted_rows"}[kind]
        t["x3"][5]["facts"][key][ARM] = False
    else:  # the shape terms at 1e-9, de_score and de_direction at 1e-4
        t["x3"][5]["builds"]["prefix"]["rows"][0]["skills"][PL][kind] += 2e-4 if kind.startswith("de_") else 2e-9
    res = e57.reading(t, [], COMMIT)
    assert res["void"] is True and res["outcome"] == "E57: VOID" and res["why"]
    assert e57.text(res).split("\n", 1)[0] == "E57: VOID" and e57.VOID_TEXT in e57.text(res)
    if kind.startswith("de_"):
        assert "more than 1e-4 off O's (prefix)" in " ".join(res["why"])


def test_a_deviation_inside_the_tolerances_does_not_void(e57):
    """The shape terms read O's coordinate bytes, so 1e-9; de_score and de_direction read a float32 pseudobulk whose
    rounding follows the row order, so 1e-4 (the row's clause): 1.7e-5, the largest move measured on synthetic
    permutations, is inside it, and so is anything up to 1e-4 itself."""
    t = table(arm=_nb(0.03))
    row = t["x3"][5]["builds"]["main"]["rows"][0]["skills"]
    row[ARM]["d2_shape"] += 5e-10
    row[ARM]["de_direction"] += 1.7e-5
    row[PL]["de_score"] -= 1e-4
    assert e57.reading(t, [], COMMIT)["void"] is False
    assert (e57.SHAPE_TOL, e57.DE_TOL) == (1e-9, 1e-4)
    t = table(arm=_nb(0.03))
    t["x3"][5]["builds"]["main"]["rows"][0]["skills"][ARM]["d2_shape"] += 2e-9  # the shape terms stay at 1e-9
    assert e57.reading(t, [], COMMIT)["void"] is True


def test_the_tolerance_the_row_names_is_the_one_the_read_applies(e57):
    m = re.fullmatch(r"de_score and de_direction more than (\S+) off O's", e57.DE_TOL_TEXT)
    assert m and float(m.group(1)) == e57.DE_TOL and e57.DE_TOL_TEXT in e57.ROW_NAMES
    assert e57.DE_TOL_TEXT in E57_ROW and e57.DE_TOL_TEXT.replace("O's", "O\\x27s") in _snippet()
    assert '"the absolute difference in de_score and in de_direction against O\'s in the same task and seed, which ' \
           'covers float32 summation-order noise"' in " ".join(e57.__doc__.split())


def test_reasons_from_outside_the_run_void_it(e57):
    res = e57.reading(table(arm=_nb(0.03)), ["step 0 did not pass"], COMMIT)
    assert res["void"] is True and res["why"][0] == "step 0 did not pass"


def test_float32_pseudobulk_follows_row_order():
    """Why the module docstring's NOTE: the scorer's pseudobulk is np.mean(X, 0) of float32 rows, so one swap of two
    rows changes its bits, and a rank-based metric then moves wherever two genes' changes sit within that noise."""
    rng = np.random.default_rng(0)
    X = np.log1p(rng.gamma(0.3, 20, size=(3000, 200))).astype(np.float32)
    Y = X.copy()
    Y[[10, 2000]] = Y[[2000, 10]]
    assert not np.array_equal(X.mean(0), Y.mean(0))
    assert np.array_equal(np.sort(X, axis=0), np.sort(Y, axis=0))  # the same rows


def test_the_multiplicative_gate_and_the_prefix_build_are_reported_only(e57):
    t = table(arm=_nb(0.03))
    for recs in t.values():
        for r in recs:
            r["builds"]["prefix"]["rows"][0]["scores"][ARM] -= 5.0  # prefix disagrees: decides nothing
    res = e57.reading(t, [], COMMIT)
    rep = res["reported_deciding_nothing"]
    assert res["outcome"] == "E57: PASS" and rep["prefix_board"]["mean"] < 0
    assert rep["multiplicative_gate_board"]["mean"] > 0 and "adoption_rule_e75_e725" in rep


# ------------------------------------------------------------------------------------------------ the ledger row
def test_the_row_must_be_one_row_e57_naming_e57a_the_script_the_launcher_and_the_tolerance(e57, tmp_path):
    cases = [(E57_ROW, True), (E57_ROW + "\n" + E57_ROW, False), ("| E56 (x) | y |", False),
             (E57_ROW.replace("E57a", "E57x"), False),
             (E57_ROW.replace("experiments/coords/e57_rs_read.py", "e57.py"), False),
             (E57_ROW.replace("experiments/coords/run_e57.sh", "run.sh"), False),
             (E57_ROW.replace("more than 1e-4 off", "more than 1e-9 off"), False),
             (E57_ROW.replace("de_score and de_direction more", "de_direction more"), False),
             (E57_ROW.replace("| E57 (", "| E57a ("), False), (E57_ROW.replace("| E57 (", "| E57 ["), False),
             ("| E57 (E57a) | `experiments/coords/e57_rs_read.py` `experiments/coords/run_e57.sh` de_score and "
              "de_direction more than 1e-4 off O's |", False)]
    code = _snippet()
    for k, (rows, want) in enumerate(cases):
        p = tmp_path / f"TRIES{k}.md"
        p.write_text("| ID | what |\n|---|---|\n" + rows + "\n")
        assert e57.row_ok(p) is want, rows
        got = subprocess.run(["python3", "-I", "-B", "-c", code, str(p)], capture_output=True, text=True)
        assert (got.returncode == 0) is want, rows  # the launcher's own check agrees


def _snippet() -> str:
    m = re.search(r"python3 -I -B -c '([^']+)' experiments/harness/TRIES.md", LAUNCHER.read_text())
    assert m, "the launcher's row check is not where the test expects it"
    return m.group(1)


# ------------------------------------------------------------------------------------------------ step 0 (E57a)
def _pair(name, **kw) -> dict:
    v = {"contrast_earlier": 0.5, "contrast_later": 0.6, "normal": [1.0, 0.0, 0.0],
         "sd_f": {"O": 0.15, "RS": 0.45, "placebo": 0.15}, "own_run_share": {"O": 0.55, "RS": 0.93, "placebo": 0.55},
         "homology_error": {"RS": 0.10, "placebo": 0.08}, "reach": 0.2, "placebo_matches_swaps": True,
         "placebo_short": False}
    for k, x in kw.items():
        if isinstance(x, dict):
            v[k] = {**v[k], **x}
        else:
            v[k] = x
    return v


def _pairs(folder: Path, test=None, skip=(), commit=COMMIT, errors=None):
    folder.mkdir(parents=True, exist_ok=True)
    for name in ("test", *TASKS):
        if name in skip:
            continue
        per = {str(1_000_000 + s): (_pair(name, **(test or {})) if name == "test" else _pair(name)) for s in range(10)}
        if name == "test" and errors:
            per[str(1_000_003)] = {"error": errors}
        (folder / f"{name}.json").write_text(json.dumps({"name": name, "commit": commit, "per_seed": per}))
    return folder


@pytest.mark.parametrize("test, fails", [
    ({}, []),
    ({"contrast_earlier": 0.71}, ["i"]),
    ({"contrast_later": 0.86}, ["i"]),
    ({"sd_f": {"RS": 0.2249}}, ["ii"]),
    ({"own_run_share": {"RS": 0.849}}, ["ii"]),
    ({"homology_error": {"RS": 0.21, "placebo": 0.08}}, ["iii"]),
    ({"homology_error": {"RS": 0.10, "placebo": 0.11}}, ["iii"]),
    ({"contrast_earlier": 0.9, "sd_f": {"RS": 0.15}, "homology_error": {"RS": 0.3}}, ["i", "ii", "iii"]),
    ({"contrast_earlier": 0.70, "contrast_later": 0.85, "sd_f": {"RS": 0.225}, "own_run_share": {"RS": 0.85},
      "homology_error": {"RS": 0.2, "placebo": 0.2}}, []),  # every bar met with equality
])
def test_step0_decides_on_the_test_pair_and_each_clause_alone_stops(e57, tmp_path, test, fails):
    res = e57.step0_reading(_pairs(tmp_path / "p", test=test))
    want = "E57a step 0: PASS" if not fails else f"E57a step 0: STOP ({', '.join(fails)})"
    assert res["outcome"] == want and res["pass"] is (not fails)
    assert {k for k, v in res["clauses"].items() if not v} == set(fails)
    assert e57.step0_text({**res, "provenance": {}}).split("\n", 1)[0] == want


def test_step0_ignores_the_practice_pairs_and_stops_on_an_incomplete_test_pair(e57, tmp_path):
    res = e57.step0_reading(_pairs(tmp_path / "a", skip=("half", "x3")))
    assert res["pass"] is True and res["practice_pairs_reported"]["half"] is None
    res = e57.step0_reading(_pairs(tmp_path / "b", skip=("test",)))
    assert res["outcome"] == "E57a step 0: STOP (incomplete)" and res["pass"] is False
    res = e57.step0_reading(_pairs(tmp_path / "c", errors="RunSortRefused: no label"))
    assert res["outcome"] == "E57a step 0: STOP (incomplete)" and "no label" in " ".join(res["notes"])


def test_the_run_starts_only_after_a_pass_at_the_run_commit(e57, tmp_path):
    def reading(outcome, commit=COMMIT, clean=True, pairs=(COMMIT,)):
        d = tmp_path / f"r{len(list(tmp_path.iterdir()))}"
        d.mkdir()
        (d / "reading.json").write_text(json.dumps({"outcome": outcome, "pair_commits": list(pairs),
                                                    "provenance": {"commit": commit, "clean": clean}}))
        (d / "reading.txt").write_text(outcome + "\nmore\n")
        return d

    assert e57.step0_passed(reading("E57a step 0: PASS"), COMMIT) == []
    assert e57.step0_passed(reading("E57a step 0: STOP (ii)"), COMMIT) == ["step 0 did not pass"]
    assert e57.step0_passed(reading("E57a step 0: PASS", commit="d" * 40), COMMIT)
    assert e57.step0_passed(reading("E57a step 0: PASS", clean=False), COMMIT)
    assert e57.step0_passed(reading("E57a step 0: PASS", pairs=(COMMIT, "d" * 40)), COMMIT)
    assert e57.step0_passed(tmp_path / "nothing", COMMIT) == [f"no step-0 reading in {tmp_path / 'nothing'}"]


def test_step0s_figures_on_a_synthetic_pair_are_the_steps_own(e57, tmp_path, monkeypatch):
    """The step-0 computation end to end on synthetic brackets (the test pair's slot, a synthetic validation board):
    ten seeds of figures, no count, and the decision drawn from them."""
    import test_run_sorted as trs

    stages = {"E7.5.h5ad": trs.stage_file("E7.5", 31), "E8.0.h5ad": trs.stage_file("E8.0", 32)}
    board = Board(key="T2:embryo:val_interp", task="T2", genes=list(trs.PANEL), needs_coords=True,
                  split="val_interp", setting="embryo", min_cells=50, max_cells=trs.SECRET_MAX)
    monkeypatch.setattr(e57, "_boards", lambda: {board.key: board})
    monkeypatch.setattr(e57, "_load", lambda path, genes: stages[Path(path).name][:, genes].copy())
    monkeypatch.setattr(tb, "interpolate_clouds", partial(tb.interpolate_clouds, n_anchor=150))
    out = tmp_path / "step0"
    out.mkdir()
    assert e57.run_pair("test", tmp_path, out, COMMIT) == 0
    text = (out / "test.json").read_text()
    r = json.loads(text)
    assert sorted(r["per_seed"]) == [str(1_000_000 + s) for s in range(10)] and r["commit"] == COMMIT
    v = r["per_seed"]["1000004"]
    assert v["sd_f"]["RS"] > v["sd_f"]["O"] and v["own_run_share"]["RS"] > v["own_run_share"]["O"]
    assert v["sd_f"]["placebo"] == v["sd_f"]["O"] and v["placebo_matches_swaps"] is True
    assert not re.search(rf"(?<![\d.]){trs.SECRET_MAX}(?![\d.])", text)  # no count
    res = e57.step0_reading(out)
    assert res["test_pair"]["complete"] is True and res["outcome"].startswith("E57a step 0: ")


# ------------------------------------------------------------------------------------------------ the gate
def test_the_gate_holds_o_to_e55s_committed_figure(e57, tmp_path):
    ref = ({"x1": 14.0}, [])
    floor = np.full(10, 50.0)
    ok = e57.gate_check("x1", {harness.FLOOR: floor, O: floor + 14.0 + np.linspace(-0.004, 0.004, 10)}, ref)
    assert ok["pass"] is True and ok["o_off_e55"] < 0.01
    bad = e57.gate_check("x1", {harness.FLOOR: floor, O: floor + 14.02}, ref)
    assert bad["pass"] is False and "0.0200 off E55" in bad["why"][0]
    assert e57.gate_check("x1", {harness.FLOOR: floor}, ref)["pass"] is False
    assert e57.gate_check("x2", {harness.FLOOR: floor, O: floor}, ref)["pass"] is False
    p = tmp_path / "reading.json"
    p.write_text("{}")
    assert e57.e55_o_per_task(p)[0] is None and "sha256" in e57.e55_o_per_task(p)[1][0]


def test_the_committed_e55_reading_gives_o_on_all_nine_tasks(e57):
    ref, why = e57.e55_o_per_task()
    assert why == [] and sorted(ref) == sorted(TASKS)


# ------------------------------------------------------------------------------------------------ the launcher
def test_the_launcher_text_holds_the_steps_the_names_and_no_deletion():
    sh = LAUNCHER.read_text()
    assert "SCRIPT=experiments/coords/e57_rs_read.py" in sh and 'S="$HOME/v/scratch/e57/$RUN"' in sh
    order = re.findall(r"^begin (\w+)$", sh, flags=re.M)
    assert order == ["preflight", "tests", "step0", "gate", "run"] and "begin read" in sh
    assert sh.index('[ "$HEAD" != "$COMMIT" ]') < sh.index("git status --porcelain") < sh.index(
        "sha256sum -c docs/S2_RULES.sha256") < sh.index("python3 -I -B -c") < sh.index("begin tests")
    assert "finish 96" in sh and "finish 98" in sh and "finish 97" in sh and "finish 94" in sh and "finish 93" in sh
    assert "STEP0_OUT=experiments/coords/results/e57a_step0" in sh and "OUT=experiments/coords/results/e57_rs" in sh
    code = [line.split("#")[0] for line in sh.splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate)\s", line) for line in code)
    assert "upload" not in " ".join(code) and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert {"tests/test_run_sorted.py", "tests/test_e57_rs_read.py"} <= set(listed)
    assert all((REPO / f).exists() for f in listed if f != "tests/test_embryo_rs_launcher.py") and (
        REPO / "tests/test_embryo_rs_launcher.py").exists()
    for line in sh.replace("\\\n", " ").splitlines():
        if re.search(r'\$PY "\$SCRIPT"', line):
            assert "> /dev/null 2>&1" in line, line


RULES = "synthetic rules for the launcher tests\n"


def _checkout(tmp_path: Path, rows=(E57_ROW,), rules_ok=True) -> Path:
    """A throwaway repository: the launcher, a synthetic rules file with its sha256 file, a synthetic ledger."""
    repo = tmp_path / "checkout"
    (repo / "experiments/coords").mkdir(parents=True)
    shutil.copyfile(LAUNCHER, repo / REL)
    (repo / "docs").mkdir()
    (repo / "docs/S2_RULES.md").write_text(RULES)
    import hashlib

    sha = hashlib.sha256(RULES.encode()).hexdigest()
    (repo / "docs/S2_RULES.sha256").write_text(f"{sha if rules_ok else '0' * 64}  docs/S2_RULES.md\n")
    tries = repo / "experiments/harness/TRIES.md"
    tries.parent.mkdir(parents=True)
    tries.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in rows))
    git(repo, "init", "-q")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "row")
    return repo


STAND_IN_UV = """#!/bin/sh
# a stand-in uv: every call is logged; each step's code comes from the environment (default 0)
echo "$*" >> "$HOME/uvcalls"
shift 3
case "$1" in
  python)
    shift 3
    case "$1" in
      -m) exit "${FAKE_TESTS:-0}" ;;
      experiments/coords/e57_rs_read.py)
        case "$2" in
          step0) exit "${FAKE_STEP0:-0}" ;;
          step0-read) exit "${FAKE_STEP0_READ:-0}" ;;
          gate) exit "${FAKE_GATE:-0}" ;;
          run) exit "${FAKE_RUN:-0}" ;;
          read) exit "${FAKE_READ:-0}" ;;
        esac ;;
    esac ;;
esac
exit 1
"""


def _home(tmp_path: Path, name: str = "home", uv: str = "#!/bin/sh\necho \"$*\" >> \"$HOME/uvcalls\"\nexit 1\n"):
    home = tmp_path / name
    fake = home / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text(uv)
    (fake / "uv").chmod(0o755)
    return home


def _run(repo: Path, home: Path, *argv: str, trace: bool = False, **env_extra) -> subprocess.CompletedProcess:
    env = {**git_env(), "HOME": str(home), **env_extra}
    return subprocess.run(["bash", *(["-x"] if trace else []), REL, *argv], cwd=repo, env=env, capture_output=True,
                          text=True, timeout=120)


def _scratch(home: Path, run: str) -> Path:
    return home / "v" / "scratch" / "e57" / run


def _exits(s: Path) -> dict:
    return {p.name: p.read_text().strip() for p in s.glob("*.exit")}


def _refused(repo, home, run, code, words):
    out = _run(repo, home, git(repo, "rev-parse", "HEAD"), run)
    assert out.returncode == code and words in out.stderr, (out.returncode, out.stderr)
    s = _scratch(home, run)
    assert _exits(s) == {"preflight.exit": str(code), "done.exit": str(code)}
    assert not (home / "uvcalls").exists()  # no python of the repository ran


def test_a_wrong_commit_is_refused_with_96_and_its_trace_shows_no_python(tmp_path):
    repo = _checkout(tmp_path)
    home = _home(tmp_path)
    out = _run(repo, home, "0" * 40, "t1", trace=True)
    assert out.returncode == 96 and "HEAD is not" in out.stderr and "finish 96" in out.stderr
    assert not re.search(r"^\++ (uv|python3?) ", out.stderr, flags=re.M)
    assert _exits(_scratch(home, "t1")) == {"preflight.exit": "96", "done.exit": "96"}


def test_bad_arguments_and_a_reused_run_folder_are_refused(tmp_path):
    repo = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    assert _run(repo, home, head).returncode == 2 and not (home / "v").exists()
    assert _run(repo, home, head, "a/b").returncode == 90
    assert _run(repo, home, head, ".x").returncode == 90
    assert _run(repo, home, head[:12], "r0").returncode == 90
    assert _run(repo, home, head.upper(), "r0").returncode == 90
    assert _run(repo, home, "0" * 40, "r0").returncode == 96
    again = _run(repo, home, head, "r0")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert _exits(_scratch(home, "r0"))["preflight.exit"] == "96"


def test_an_unclean_tree_or_a_changed_rules_file_is_refused(tmp_path):
    repo = _checkout(tmp_path / "a")
    (repo / "experiments/harness/TRIES.md").write_text(E57_ROW + "\n| E58 (changed) | x |\n")
    _refused(repo, _home(tmp_path, "h1"), "u1", 98, "not clean")
    repo = _checkout(tmp_path / "b", rules_ok=False)
    _refused(repo, _home(tmp_path, "h2"), "u2", 97, "S2_RULES")


@pytest.mark.parametrize("rows", [(), (E57_ROW, E57_ROW), (E57_ROW.replace("E57a", "the step"),),
                                  (E57_ROW.replace("experiments/coords/run_e57.sh", "the launcher"),),
                                  (E57_ROW.replace("1e-4", "1e-9"),),
                                  (E57_ROW.replace("| E57 (", "| E57a ("),)])
def test_without_one_row_e57_naming_its_code_the_preflight_refuses(tmp_path, rows):
    repo = _checkout(tmp_path, rows=rows)
    _refused(repo, _home(tmp_path), "w1", 94, "row E57")


def test_an_existing_reading_is_refused_before_any_step(tmp_path):
    for k, d in enumerate(("experiments/coords/results/e57a_step0", "experiments/coords/results/e57_rs")):
        repo = _checkout(tmp_path / f"c{k}")
        (repo / d).mkdir(parents=True)
        (repo / d / "reading.txt").write_text("E57: PASS\n")  # untracked: the tree stays clean
        _refused(repo, _home(tmp_path, f"h{k}"), f"e{k}", 93, "already holds a reading")


def test_a_passing_preflight_runs_the_tests_and_a_failing_uv_stops_there(tmp_path):
    repo = _checkout(tmp_path)
    home = _home(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    out = _run(repo, home, head, "g1")
    s = _scratch(home, "g1")
    assert out.returncode == 1 and _exits(s) == {"preflight.exit": "0", "tests.exit": "1", "done.exit": "1"}
    calls = (home / "uvcalls").read_text().splitlines()
    assert len(calls) == 1 and "pytest" in calls[0] and "tests/test_e57_rs_read.py" in calls[0]
    assert (s / "head").read_text().strip() == head


def test_the_steps_run_in_order_and_stop_where_the_row_says(tmp_path):
    repo = _checkout(tmp_path)
    head = git(repo, "rev-parse", "HEAD")
    # every step exits 0: tests, ten step-0 pairs, the step-0 read, nine gates, nine runs, the read
    home = _home(tmp_path, "h0", STAND_IN_UV)
    out = _run(repo, home, head, "all")
    s = _scratch(home, "all")
    assert out.returncode == 0, out.stderr
    assert _exits(s) == {f"{k}.exit": "0" for k in ("preflight", "tests", "step0", "gate", "run", "read", "done")}
    calls = (home / "uvcalls").read_text().splitlines()
    kinds = ["tests" if " -m pytest " in c else re.search(r"e57_rs_read\.py (\S+)", c).group(1) for c in calls]
    assert kinds == ["tests", *["step0"] * 10, "step0-read", *["gate"] * 9, *["run"] * 9, "read"]
    assert all(f"--commit {head}" in c for c in calls[1:]) and "--name test " in " ".join(calls)
    # step 0 stops: E57b and E57c never run, and nothing is read
    home = _home(tmp_path, "h1", STAND_IN_UV)
    out = _run(repo, home, head, "stop", FAKE_STEP0_READ="3")
    assert out.returncode == 3
    assert _exits(_scratch(home, "stop")) == {"preflight.exit": "0", "tests.exit": "0", "step0.exit": "3",
                                              "done.exit": "3"}
    assert not any(" gate " in c or " run " in c or " read " in c for c in (home / "uvcalls").read_text().splitlines())
    # the gate fails: the run never starts, the read is called once (VOID) and the call exits with the gate's code
    home = _home(tmp_path, "h2", STAND_IN_UV)
    out = _run(repo, home, head, "gate", FAKE_GATE="3")
    assert out.returncode == 3
    assert _exits(_scratch(home, "gate")) == {"preflight.exit": "0", "tests.exit": "0", "step0.exit": "0",
                                              "gate.exit": "3", "read.exit": "0", "done.exit": "3"}
    calls = (home / "uvcalls").read_text().splitlines()
    assert sum(" run " in c for c in calls) == 0 and sum(" read " in c for c in calls) == 1
    # a failing test step stops everything
    home = _home(tmp_path, "h3", STAND_IN_UV)
    out = _run(repo, home, head, "tests", FAKE_TESTS="1")
    assert out.returncode == 1 and _exits(_scratch(home, "tests"))["tests.exit"] == "1"
    assert len((home / "uvcalls").read_text().splitlines()) == 1
