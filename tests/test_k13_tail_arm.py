"""Tests for tries K13 (experiments/t3_decomp/k13_tail_arm.py), on synthetic data only: no competition file is opened,
nothing is scored, and no prediction is compared with any target."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

T3 = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import eb_factor as eb  # noqa: E402
import k7_wt_revert as k7  # noqa: E402
import k9_residual as k9  # noqa: E402
import k13_tail_arm as k13  # noqa: E402
import test_k9_residual as t9  # noqa: E402  (K9's synthetic world and reference, reused)

CLIP = t9.CLIP
world = t9.world  # the module fixture: synthetic files in task x1's roles


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


# ----------------------------------------------------------------------------------------------- constants and tasks
def test_the_fixed_constants_and_arms():
    assert k13.TAIL_PB_MIN == 0.25
    assert k13.PERM_SEED == k9.PERM_SEED == 20261113 and k13.STREAM == k9.STREAM == harness.K4_STREAM
    assert (k13.MIN_GAIN, k13.MIN_AHEAD, k13.X2_FLOOR) == (0.3, 8, -0.3)
    assert k13.ARMS == ("raw", "raw_tail", "raw_tail_ctrl") and k13.CONTROLS == {"raw_tail": "raw_tail_ctrl"}
    assert k13.TAIL_STREAM == 2 and not hasattr(k13, "permuted_mask") and not hasattr(k13, "TAIL_PERM")
    assert k13.TASKS is k9.TASKS and k13.RESULTS_K9.is_dir()
    assert [c.name for c in k13.candidates()] == [harness.FLOOR, *k13.ARMS]
    assert all(c.applies(t) for c in k13.candidates()[1:] for t in k13.TASKS.values())
    assert not any(c.applies(t) for c in k13.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)


# ----------------------------------------------------------------------------------------------- the filter
def test_the_filter_zeroes_exactly_the_negative_entries_below_the_level():
    rng = np.random.default_rng(1)
    n = 200
    level = rng.gamma(1.0, 0.5, n)
    y = rng.normal(0, 0.7, n)
    inf = k13.tail_infeasible(level)
    assert np.array_equal(inf, level < 0.25) and 0 < inf.sum() < n
    out, hit = k13.tail_filter(y, inf)
    assert np.array_equal(hit, inf & (y < 0)) and hit.any()
    assert np.all(out[hit] == 0.0)
    assert np.array_equal(out[~hit], y[~hit])  # every other entry as it came, bit for bit
    # positive entries on infeasible genes and negative ones on feasible genes are left
    up, low = inf & (y > 0), ~inf & (y < 0)
    assert np.array_equal(out[up], y[up]) and np.array_equal(out[low], y[low])
    den = rng.gamma(2.0, size=n)
    phi_raw = eb.finalize(y, den, CLIP)
    assert np.array_equal(k13.tail_factor(y, den, CLIP, inf, phi_raw), eb.finalize(out, den, CLIP))
    # the level's edge: exactly 0.25 is feasible
    assert list(k13.tail_infeasible(np.array([0.2499, 0.25, 0.3]))) == [True, False, False]


def test_no_infeasible_gene_or_no_negative_entry_gives_raw_at_the_factor():
    rng = np.random.default_rng(2)
    n = 60
    y = rng.normal(0, 0.7, n)
    den = rng.gamma(2.0, size=n)
    phi_raw = np.clip(np.exp(rng.normal(0, 0.6, n)), *CLIP)  # not finalize(y): the idle arm must return it as it is
    none = np.zeros(n, dtype=bool)
    assert np.array_equal(k13.tail_factor(y, den, CLIP, none, phi_raw), phi_raw)
    pos_only = y > 0  # infeasible genes, but none with a negative entry
    assert np.array_equal(k13.tail_factor(y, den, CLIP, pos_only, phi_raw), phi_raw)
    assert k13.tail_check(y, pos_only, pos_only, 0)["negative"] == 0
    ctrl, short = k13.control_mask(y, pos_only)
    assert not ctrl.any() and short == 0  # the arm changes nothing, so the control changes nothing


def test_no_changed_gene_gives_raw_bit_for_bit_at_the_built_cells(bs):
    rng = np.random.default_rng(6)
    ref = t9._reference(bs, rng)
    y = rng.normal(0, 0.7, t9.G)
    den = rng.gamma(2.0, size=t9.G)
    phi_raw = np.clip(np.exp(rng.normal(0, 0.6, t9.G)), *CLIP)
    none = np.zeros(t9.G, dtype=bool)
    prep = {"phi_raw": phi_raw, "phi_tail": k13.tail_factor(y, den, CLIP, none, phi_raw),
            "phi_tail_ctrl": k13.tail_factor(y, den, CLIP, none, phi_raw), "clip": CLIP}
    raw = k13.arm_built(bs, ref, 200, 7, ref.V.p, "raw", prep)
    k9raw = k9.arm_built(bs, ref, 200, 7, ref.V.p, "raw", phi_raw, None, None, CLIP)
    assert np.array_equal(raw.X, k9raw.X) and np.array_equal(raw.row, k9raw.row)
    for arm in ("raw_tail", "raw_tail_ctrl"):
        b = k13.arm_built(bs, ref, 200, 7, ref.V.p, arm, prep)
        assert np.array_equal(b.X, raw.X) and np.array_equal(b.row, raw.row) and np.array_equal(b.lab, raw.lab)


def test_the_arms_draw_the_same_cells_and_the_filter_moves_them(bs):
    rng = np.random.default_rng(8)
    ref = t9._reference(bs, rng)
    y = rng.normal(0, 0.7, t9.G)
    den = rng.gamma(2.0, size=t9.G)
    phi_raw = eb.finalize(y, den, CLIP)
    inf = (y < 0) & (np.arange(t9.G) % 2 == 0)  # half the negative entries infeasible
    ctrl, _ = k13.control_mask(y, inf)
    prep = {"phi_raw": phi_raw, "phi_tail": k13.tail_factor(y, den, CLIP, inf, phi_raw),
            "phi_tail_ctrl": k13.tail_factor(y, den, CLIP, ctrl, phi_raw), "clip": CLIP}
    built = {a: k13.arm_built(bs, ref, 200, 9, ref.V.p, a, prep) for a in k13.ARMS}
    assert all(np.array_equal(b.row, built["raw"].row) for b in built.values())
    assert not np.array_equal(built["raw_tail"].X, built["raw"].X)


def _panel(n: int, frac_infeasible: float, seed: int):
    rng = np.random.default_rng(seed)
    level = np.where(rng.random(n) < frac_infeasible, rng.uniform(0.0, 0.2499, n), rng.uniform(0.25, 3.0, n))
    return level, rng.normal(0, 0.7, n)


@pytest.mark.parametrize("frac", [0.9, 0.1])
def test_the_control_changes_disjoint_feasible_genes_as_many_as_the_arm_or_all_of_them(frac):
    level, y = _panel(400, frac, 11)
    inf = k13.tail_infeasible(level)
    neg = y < 0
    _, hit = k13.tail_filter(y, inf)
    ctrl, short = k13.control_mask(y, inf)
    _, hit_c = k13.tail_filter(y, ctrl)
    assert np.array_equal(hit_c, ctrl)  # the control's mask is exactly the genes it changes
    assert not (hit & hit_c).any()  # disjoint from the arm's genes
    assert not (ctrl & inf).any() and np.all(neg[ctrl])  # feasible genes with a negative entry only
    pool = int((~inf & neg).sum())
    if pool >= hit.sum():  # enough feasible genes with a negative entry: exactly as many as the arm
        assert frac == 0.1 and short == 0 and hit_c.sum() == hit.sum()
    else:  # too few: all of them, and the shortfall recorded
        assert frac == 0.9 and hit_c.sum() == pool and short == hit.sum() - pool > 0
    chk = k13.tail_check(y, inf, ctrl, short)
    assert chk == {"infeasible": int(inf.sum()), "negative": int(hit.sum()), "control_changed": int(hit_c.sum()),
                   "control_short": short}
    # same seed, same control; another seed, another draw where a choice exists
    again, short2 = k13.control_mask(y, inf)
    assert np.array_equal(again, ctrl) and short2 == short
    if frac == 0.1:
        other, _ = k13.control_mask(y, inf, seed=k9.PERM_SEED + 1)
        assert other.sum() == ctrl.sum() and not np.array_equal(other, ctrl)


def test_the_controls_generator_is_a_stream_of_its_own():
    draw = k13.control_rng().random(64)
    assert np.array_equal(draw, np.random.default_rng([k9.PERM_SEED, k13.TAIL_STREAM]).random(64))
    assert not np.array_equal(draw, np.random.default_rng(k9.PERM_SEED).random(64))
    assert not np.array_equal(draw, np.random.default_rng([k9.PERM_SEED, 1]).random(64))
    # and the control's draw itself differs from the same draw on either of those generators
    level, y = _panel(400, 0.1, 12)
    inf = k13.tail_infeasible(level)
    ctrl, _ = k13.control_mask(y, inf)
    pool = np.flatnonzero(~inf & (y < 0))
    take = int((inf & (y < 0)).sum())
    for rng in (np.random.default_rng(k9.PERM_SEED), np.random.default_rng([k9.PERM_SEED, 1])):
        m = np.zeros(len(y), dtype=bool)
        m[rng.choice(pool, size=take, replace=False)] = True
        assert not np.array_equal(m, ctrl)


def test_the_check_line_fields():
    chk = k13.tail_check(np.array([-1.0, 1.0, -0.5, 0.2, -0.1]), np.array([1, 1, 1, 0, 0], bool),
                         np.array([0, 1, 0, 0, 1], bool), 1)
    assert chk == {"infeasible": 3, "negative": 2, "control_changed": 1, "control_short": 1}
    assert k13.check_line("x1", chk) == ("raw_tail check: task=x1 infeasible=3 negative=2 control_changed=1 "
                                         "control_short=1")
    assert "both_changed" not in k13.check_line("x1", chk)


# ----------------------------------------------------------------------------------------------- end to end, synthetic
def test_end_to_end_raw_is_k9s_raw_and_raw_tail_moves_it(world, monkeypatch):
    monkeypatch.setattr(k7, "MIN_CELLS", 20)  # K9's own build of the same fold needs its groups (as k12's test)
    t = k13.TASKS["x1"]
    k13._PREP.clear()
    k9._PREP.clear()
    ctx = t9._ctx(t, world)
    level = harness.t3_mean_expression(world[ctx.task.reference_role])
    # the synthetic world's levels are high; the level is moved to its median here only to exercise the filter
    monkeypatch.setattr(k13, "TAIL_PB_MIN", float(np.median(level)))
    with harness.guarded_reads() as g:
        out = {arm: k13._build(ctx, arm) for arm in k13.ARMS}
        k9raw = k9._build(ctx, k9.BASE)
    assert not g["refused"]
    prep = k13._PREP["fold"]
    assert np.array_equal(prep["phi_raw"], k9._PREP["fold"]["phi_raw"])
    assert np.array_equal(out["raw"].X, k9raw.X) and np.array_equal(out["raw"].coords, k9raw.coords)
    rec = prep["record"]
    assert json.loads(json.dumps(rec)) == rec
    chk = rec["tail_check"]
    assert chk["infeasible"] == int((level < np.median(level)).sum()) and chk["negative"] > 0
    assert chk["control_changed"] + chk["control_short"] == chk["negative"] and "both_changed" not in rec
    assert not (prep["infeasible"] & prep["control"]).any()
    assert not np.array_equal(out["raw_tail"].X, out["raw"].X)
    # the true level: no gene of the synthetic reference below 0.25, or none negative there, gives raw bit for bit
    monkeypatch.setattr(k13, "TAIL_PB_MIN", -1.0)
    k13._PREP.clear()
    with harness.guarded_reads():
        for a in ("raw_tail", "raw_tail_ctrl"):
            off = k13._build(ctx, a)
            assert np.array_equal(off.X, out["raw"].X) and np.array_equal(off.coords, out["raw"].coords), a
    assert k13._PREP["fold"]["record"]["tail_check"]["negative"] == 0
    k13._PREP.clear()
    k9._PREP.clear()


# ----------------------------------------------------------------------------------------------- the reading
def _result(key, scores, clean=True, commit="abc", negative=5, short=0):
    rows = [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5, "severity": 0.6} for a in scores},
             "t3": {a: {"severity": {"state": "scored"}} for a in scores}} for i, s in enumerate(SEEDS)]
    chk = {"infeasible": 9, "negative": negative, "control_changed": max(negative - short, 0), "control_short": short}
    mass = {"wt_mass_share_arm": 0.01, "wt_mass_share_control": 0.2, "finalize_divisor_raw": 1.3,
            "finalize_divisor_arm": 1.31, "finalize_divisor_control": 1.37}
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
            "factor_record": {"tail_check": chk, "mass": mass}, "tail_check": chk}


def _scores(tail=54.0, perm=50.0, raw=52.0):
    n = len(SEEDS)
    return {harness.FLOOR: [50.0] * n, "raw": [raw] * n, "raw_tail": [tail] * n, "raw_tail_ctrl": [perm] * n}


def _k9_dir(d: Path, raw: float = 52.0) -> Path:
    """A stand-in for K9's committed task results: raw scores `raw` on every seed of every task."""
    d.mkdir()
    for name in k13.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k9.TASKS[name].key, {"raw": [raw] * len(SEEDS)})))
    return d


def _read(d: Path, k9_dir=None, negative=None, short=None, **per_task):
    d.mkdir()
    negative, short = negative or {}, short or {}
    for name, sc in per_task.items():
        (d / f"{name}.json").write_text(json.dumps(_result(k13.TASKS[name].key, sc, negative=negative.get(name, 5),
                                                           short=short.get(name, 0))))
    return k13.reading(d, k9_dir if k9_dir is not None else _k9_dir(d / "k9"))


def test_the_reading_supported(tmp_path):
    yes = _read(tmp_path / "yes", x1=_scores(), r12=_scores(), x2=_scores(tail=51.8))
    assert not yes["void"] and yes["supported"] and yes["verdict"].startswith("SUPPORTED")
    assert yes["k9_reproduction"] == "reproduced" and yes["idle_folds"] == []
    txt = k13.text(yes)
    assert "NOT the adoption rule" in txt and "raw_tail check: task=x1 infeasible=9 negative=5" in txt


def test_the_reading_closed_by_x2(tmp_path):
    res = _read(tmp_path / "x2", x1=_scores(), r12=_scores(), x2=_scores(tail=51.6))
    assert not res["void"] and not res["supported"] and res["verdict"].startswith("CLOSED")
    assert not res["checks"]["x2: raw_tail - raw not below -0.3"]


@pytest.mark.parametrize("fold", ["x1", "r12"])
def test_the_reading_closed_by_the_control(tmp_path, fold):
    sc = {n: _scores() for n in k13.TASKS}
    sc[fold] = _scores(perm=53.9)
    res = _read(tmp_path / fold, **sc)
    assert not res["void"] and not res["supported"] and not res["checks"][f"{fold}: raw_tail - raw_tail_ctrl"]
    assert res["checks"][f"{fold}: raw_tail - raw"] and res["verdict"] == "CLOSED: no request"


def test_the_reading_needs_eight_seeds_ahead(tmp_path):
    sc = _scores()
    sc["raw_tail"] = [70.0] * 7 + [40.0] * (len(SEEDS) - 7)
    assert not _read(tmp_path / "a", x1=sc, r12=_scores(), x2=_scores())["supported"]


@pytest.mark.parametrize("fold", ["x1", "r12"])
def test_the_reading_idle_on_a_deciding_fold_is_closed(tmp_path, fold):
    # the scores would pass; the fold's check says the filter changed no gene there
    res = _read(tmp_path / fold, negative={fold: 0}, x1=_scores(), r12=_scores(), x2=_scores())
    assert not res["void"] and not res["supported"] and res["idle_folds"] == [fold]
    assert res["verdict"].startswith("CLOSED") and "idle" in res["verdict"] and fold in res["verdict"]
    assert f"{fold} " in k13.text(res) and "IDLE" in k13.text(res)


@pytest.mark.parametrize("fold", ["x1", "r12"])
def test_the_reading_with_a_short_control_on_a_deciding_fold_is_closed(tmp_path, fold):
    # the scores would pass; the fold's control found fewer feasible genes than the arm changed
    res = _read(tmp_path / fold, short={fold: 2}, x1=_scores(), r12=_scores(), x2=_scores())
    assert not res["void"] and not res["supported"] and res["control_short_folds"] == [fold]
    assert res["verdict"].startswith("CLOSED") and "falls short" in res["verdict"] and fold in res["verdict"]
    assert all(res["checks"].values())  # every score clause holds; the control clause is what cannot be read


def test_the_reading_with_a_short_control_on_x2_alone_decides_nothing_and_the_text_shows_the_mass(tmp_path):
    res = _read(tmp_path / "x2", short={"x2": 3}, x1=_scores(), r12=_scores(), x2=_scores(tail=52.0))
    assert not res["void"] and res["supported"] and res["control_short_folds"] == []
    assert "WT mass share on the changed genes: arm 0.0100, control 0.2000" in k13.text(res)


def test_the_mass_record_shows_what_each_mask_moves():
    y = np.array([-1.0, -0.5, 0.3, -0.2, 0.0, -2.0])
    den = np.array([0.01, 0.02, 0.3, 0.4, 0.07, 0.2])
    infeasible = np.array([True, True, False, False, False, False])
    control = np.array([False, False, False, True, False, True])
    m = k13.mass_record(y, den, infeasible, control)
    assert m["wt_mass_share_arm"] == pytest.approx(0.03) and m["wt_mass_share_control"] == pytest.approx(0.6)
    assert m["finalize_divisor_raw"] == pytest.approx(float(den @ np.exp(y)))  # den sums to 1 here
    ya, yc = y.copy(), y.copy()
    ya[[0, 1]] = 0.0
    yc[[3, 5]] = 0.0
    assert m["finalize_divisor_arm"] == pytest.approx(float(den @ np.exp(ya)))
    assert m["finalize_divisor_control"] == pytest.approx(float(den @ np.exp(yc)))
    raw = m["finalize_divisor_raw"]
    assert m["finalize_divisor_control"] - raw > m["finalize_divisor_arm"] - raw  # the control moves more mass


def test_the_reading_idle_on_x2_alone_is_reported_and_decides_nothing(tmp_path):
    res = _read(tmp_path / "x2", negative={"x2": 0}, x1=_scores(), r12=_scores(), x2=_scores(tail=52.0))
    assert not res["void"] and res["supported"] and res["idle_folds"] == ["x2"]


def test_the_reading_is_void_on_a_raw_mismatch(tmp_path):
    k9_dir = _k9_dir(tmp_path / "k9")
    ok = _read(tmp_path / "ok", k9_dir, x1=_scores(), r12=_scores(), x2=_scores())
    assert not ok["void"] and ok["k9_reproduction"] == "reproduced"
    off = _read(tmp_path / "off", k9_dir, x1=_scores(raw=52.01), r12=_scores(), x2=_scores())
    assert off["void"] and any("does not reproduce" in w for w in off["why"])
    assert k13.text(off).startswith("VOID")


def test_the_reading_is_void_without_k9s_committed_scores(tmp_path):
    res = _read(tmp_path / "r", tmp_path / "absent", x1=_scores(), r12=_scores(), x2=_scores())
    assert res["void"] and res["why"] == [f"K9's committed raw scores ({tmp_path / 'absent'}) are absent from this "
                                          f"checkout"]


def test_the_reading_is_void_on_missing_arm_seed_check_or_dirty_tree(tmp_path):
    k9_dir = _k9_dir(tmp_path / "k9")
    assert _read(tmp_path / "m", x1=_scores(), r12=_scores())["void"]
    sc = _scores()
    del sc["raw_tail_ctrl"]
    assert _read(tmp_path / "arm", x1=sc, r12=_scores(), x2=_scores())["void"]
    for case in ("seed", "dirty", "check"):
        d = tmp_path / case
        d.mkdir()
        for name in k13.TASKS:
            r = _result(k13.TASKS[name].key, _scores(), clean=not (case == "dirty" and name == "x2"))
            if case == "seed" and name == "x1":
                r["rows"] = r["rows"][:-1]
            if case == "check" and name == "r12":
                del r["tail_check"]
            (d / f"{name}.json").write_text(json.dumps(r))
        res = k13.reading(d, k9_dir)
        assert res["void"], case


# ----------------------------------------------------------------------------------------------- the commands
def test_run_writes_the_check_numbers_and_prints_the_line(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(harness, "run_task", lambda *a, **kw: {"rows": []})
    chk = {"infeasible": 7, "negative": 3, "control_changed": 2, "control_short": 1}
    monkeypatch.setitem(k13._PREP, "fold", {"record": {"tail_check": chk}})
    runs = tmp_path / "runs"
    k13.main(["run", "--task", "r12", "--out", str(runs), "--data", str(tmp_path / "nowhere")])
    res = json.loads((runs / "r12.json").read_text())
    assert res["tail_check"] == chk and res["factor_record"]["tail_check"] == chk
    assert ("raw_tail check: task=r12 infeasible=7 negative=3 control_changed=2 control_short=1"
            in capsys.readouterr().out)


def test_run_and_read_refuse_to_overwrite(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(harness, "run_task", lambda *a, **kw: called.append(a))
    runs = tmp_path / "runs"
    runs.mkdir()
    (runs / "x1.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k13.main(["run", "--task", "x1", "--out", str(runs), "--data", str(tmp_path / "nowhere")])
    assert not called and (runs / "x1.json").read_text() == "{}"
    out = tmp_path / "out"
    out.mkdir()
    for f in ("reading.json", "reading.txt"):
        (out / f).write_text("{}")
        with pytest.raises(SystemExit, match="exists"):
            k13.main(["read", "--runs", str(runs), "--out", str(out)])
        (out / f).unlink()
