"""Tests for tries K15a and K15b (experiments/t3_decomp/k15_file_form_half.py), on synthetic data only: no competition
file is opened, nothing is scored, and no prediction is compared with any target."""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS

T3 = Path(__file__).resolve().parents[1] / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import k7_wt_revert as k7  # noqa: E402
import k9_residual as k9  # noqa: E402
import k10_s1_form as k10  # noqa: E402
import k14_form_arms as k14  # noqa: E402
import k15_file_form_half as k15  # noqa: E402
import test_k9_residual as t9  # noqa: E402  (K9's synthetic reference and build context)
import test_k10_s1_form as t10  # noqa: E402  (K10's synthetic replicate world)

G = t9.G
world = t10.world  # the module fixture: synthetic files in r12's and r21's roles
SOURCE = T3 / "k15_file_form_half.py"
REAL_PINS = dict(k15.K10_SHA256)  # try K10's committed files' sha256, before any test re-pins synthetic comparators
PINS: dict[str, str] = {}


@pytest.fixture(autouse=True)
def _synthetic_pins(monkeypatch):
    """Synthetic comparators are pinned by their own sha256 as _k10_dir writes them; the real pins are checked against
    the committed files in test_k10s_committed_files_are_on_the_branch_and_pass_the_comparator_checks."""
    PINS.clear()
    PINS.update(REAL_PINS)
    monkeypatch.setattr(k15, "K10_SHA256", PINS)


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


def _setup(bs, seed=4):
    rng = np.random.default_rng(seed)
    ref = t9._reference(bs, rng)
    phi = np.clip(np.exp(rng.normal(0, 0.6, G)), *t9.CLIP)
    ratio = np.array([0.7, 1.3, 1.0])
    basis = rng.gamma(2.0, size=G)
    return ref, phi, ratio, basis


# ----------------------------------------------------------------------------------------------- constants and tasks
def test_the_fixed_constants_arms_and_folds():
    assert k15.ARMS == ("s1x", "s1x_half", "s1x_half_perm") and k15.BASE == k10.BASE
    assert k15.PERM_SEED == k9.PERM_SEED == 20261113 and k15.STREAM == k10.STREAM == k9.STREAM
    assert (k15.MIN_GAIN, k15.MIN_AHEAD, k15.FLOOR_MEAN) == (0.3, 8, -0.3)
    assert set(k15.TASKS) == {"r12", "r21"} and all(k15.TASKS[n] is k10.TASKS[n] for n in k15.TASKS)
    assert k15.SOURCE is k10.SOURCE and k15.RESULTS_K10 == T3 / "results" / "k10_s1_form" / "tasks"
    assert k15.K10_SOURCE == "experiments/t3_decomp/k10_s1_form.py"
    assert k15.K10_COMMIT == "69057bc64054ee8b6f3efce04062f0b728adb0c1"
    assert k15.K10_ARMS == {"s1x", "s1x_res", "s1x_res_perm"}
    assert [c.name for c in k15.candidates()] == [harness.FLOOR, *k15.ARMS]
    assert all(c.applies(t) for c in k15.candidates()[1:] for t in k15.TASKS.values())
    assert not any(c.applies(t) for c in k15.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)
    assert not any(c.applies(t) for c in k15.candidates()[1:] for t in k9.TASKS.values())


def test_the_source_opens_no_file_beyond_the_tasks_inputs():
    src = SOURCE.read_text()
    for call in ("read_h5ad", "load_stage", "read_stage", "h5py", "np.load", "open(", "boards/", "runs/", "data/"):
        assert call not in src, call
    assert re.findall(r'"[^"]*\.h5ad"', src) == []


def test_k10s_committed_files_are_on_the_branch_and_pass_the_comparator_checks():
    """must_fix 1: try K10's committed task files are on this branch, byte for byte the pinned ones, and the
    comparator checks accept them (only their provenance fields are read here, never a score)."""
    assert k15.RESULTS_K10.is_dir()
    for name in k15.TASKS:
        p = k15.RESULTS_K10 / f"{name}.json"
        assert hashlib.sha256(p.read_bytes()).hexdigest() == REAL_PINS[name]
        r, why, rec = k15.k10_comparator(name, p)
        assert why == [] and r is not None, why
        assert rec == {"path": str(p), "sha256": REAL_PINS[name], "commit": k15.K10_COMMIT}


# ----------------------------------------------------------------------------------------------- the factors
def test_the_half_factor_is_k14bs_form():
    phi = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
    assert np.array_equal(k15.half_factor(phi), [0.625, 0.75, 1.0, 1.5, 2.5])
    assert np.array_equal(k15.half_factor(phi), k14.half_ctrl_factor(phi))
    assert np.array_equal(k15.half_factor(np.ones(7)), np.ones(7))


def test_the_control_permutes_within_k9s_deciles_with_its_seed_and_keeps_the_distribution():
    rng = np.random.default_rng(3)
    n = 57
    half = k15.half_factor(np.exp(rng.normal(0, 0.7, n)))
    basis = rng.gamma(2.0, size=n)
    perm = k15.permuted_half(half, basis)
    assert np.array_equal(perm, harness.t3_decile_permutation(half, basis, k9.PERM_SEED))
    assert np.array_equal(perm, harness.t3_decile_permutation(half, basis, k15.PERM_SEED, harness.K2_PERM_BINS))
    assert not np.array_equal(perm, half)
    assert np.array_equal(np.sort(perm), np.sort(half))  # the factor's distribution, value for value
    for genes in np.array_split(np.argsort(basis, kind="stable"), harness.K2_PERM_BINS):
        assert np.array_equal(np.sort(perm[genes]), np.sort(half[genes]))  # each decile keeps its values
    # the same permutation as of the deviation from one
    dev = harness.t3_decile_permutation(half - 1.0, basis, k9.PERM_SEED)
    assert np.allclose(perm - 1.0, dev, rtol=0, atol=1e-15)
    idx = harness.t3_decile_permutation(np.arange(n, dtype=np.float64), basis, k9.PERM_SEED).astype(int)
    assert np.array_equal(perm, half[idx])
    assert not np.array_equal(k15.permuted_half(half, basis, k9.PERM_SEED + 1), perm)


def test_the_depth_summary_and_the_counts_spy():
    totals = np.array([10.0, 20.0, 30.0, 50.0])
    lab = np.array([0, 0, 2, 2])
    d = k15.depth_summary(totals, lab, ["a", "b", "c"])
    assert d == {"overall": 27.5, "per_label": {"a": 15.0, "c": 40.0}}

    class T:
        @staticmethod
        def lognorm(C):
            return np.asarray(C) * 2.0

    t = T()
    real = t.lognorm
    with k15.counts_seen(t) as seen:
        out = t.lognorm(np.array([[1, 2], [3, 4]]))
    assert np.array_equal(out, [[2.0, 4.0], [6.0, 8.0]]) and t.lognorm is real
    assert len(seen) == 1 and np.array_equal(seen[0], [3.0, 7.0])
    fig = k15.depth_figures({"s1x": {"1": {"overall": 10.0, "per_label": {"a": 8.0, "b": 12.0}},
                                     "2": {"overall": 12.0, "per_label": {"a": 10.0}}},
                             "s1x_half": {"1": {"overall": 9.0, "per_label": {"a": 9.0, "b": 9.0}},
                                          "2": {"overall": 13.0, "per_label": {"a": 9.0}}}})
    assert fig["s1x"]["overall"] == 11.0 and fig["s1x"]["per_label"] == {"a": 9.0, "b": 12.0}
    assert fig["s1x_half"]["overall_vs_s1x"] == 1.0
    assert fig["s1x_half"]["per_label_vs_s1x"] == {"a": 1.0, "b": 0.75}


def test_the_share_beyond_one_and_a_half_fold():
    assert k15.share_beyond(np.array([1.0, 1.5, 1.51, 1 / 1.6, 0.9, 3.0])) == pytest.approx(3 / 6)
    phi = np.array([4.0, 0.25, 1.2, 2.0])
    assert k15.share_beyond(phi) == 0.75 and k15.share_beyond(k15.half_factor(phi)) == 0.5


# ----------------------------------------------------------------------------------------------- the arms
def test_s1x_is_k10s_s1x_and_the_half_is_its_call_with_only_the_factor_changed(bs, monkeypatch):
    ref, phi, ratio, basis = _setup(bs)
    perm = k15.permuted_half(k15.half_factor(phi), basis)
    base = k15.arm_built(bs, ref, 200, 7, ref.V.p, "s1x", phi, ratio, perm)
    want = k10.arm_built(bs, ref, 200, 7, ref.V.p, "s1x", phi, ratio, None, None, t9.CLIP)
    assert np.array_equal(base.X, want.X) and np.array_equal(base.row, want.row)
    calls = []
    real = bs.construct

    def spy(ref_, n, rng, **kw):
        calls.append((ref_, n, rng.bit_generator.state, kw))
        return real(ref_, n, rng, **kw)

    monkeypatch.setattr(bs, "construct", spy)
    built = {a: k15.arm_built(bs, ref, 200, 7, ref.V.p, a, phi, ratio, perm) for a in k15.ARMS}
    assert len(calls) == 3
    (r0, n0, s0, kw0), *rest = calls
    for (r, n, s, kw), fac in zip(rest, (k15.half_factor(phi), perm), strict=True):
        assert r is r0 and n == n0 and s == s0  # the same reference, size and generator state
        assert set(kw) == set(kw0) == {"f", "phi", "ratio"}
        assert kw["f"] is kw0["f"] and kw["ratio"] is kw0["ratio"]  # the depth step unchanged
        assert np.array_equal(kw["phi"], fac)
    assert np.array_equal(kw0["phi"], phi)
    for a in k15.ARMS:
        assert np.array_equal(built[a].row, base.row) and np.array_equal(built[a].lab, base.lab)
    assert not np.array_equal(built["s1x_half"].X, base.X)
    assert not np.array_equal(built["s1x_half"].X, built["s1x_half_perm"].X)


def test_a_factor_of_one_everywhere_gives_s1x_bit_for_bit(bs):
    ref, _, ratio, basis = _setup(bs, 5)
    one = np.ones(G)
    perm = k15.permuted_half(k15.half_factor(one), basis)
    base = k15.arm_built(bs, ref, 200, 9, ref.V.p, "s1x", one, ratio, perm)
    for arm in ("s1x_half", "s1x_half_perm"):
        b = k15.arm_built(bs, ref, 200, 9, ref.V.p, arm, one, ratio, perm)
        assert np.array_equal(b.X, base.X) and np.array_equal(b.row, base.row) and np.array_equal(b.lab, base.lab)


def test_the_arm_refusals(bs):
    ref, phi, ratio, _ = _setup(bs)
    with pytest.raises(ValueError, match="permuted"):
        k15.arm_built(bs, ref, 50, 1, ref.V.p, "s1x_half_perm", phi, ratio, None)
    with pytest.raises(ValueError, match="unknown arm"):
        k15.arm_built(bs, ref, 50, 1, ref.V.p, "s1x_res", phi, ratio, phi)


@pytest.mark.parametrize("name", ["r12", "r21"])
def test_end_to_end_s1x_is_k10s_build_and_the_arms_take_the_folds_own_source(world, monkeypatch, bs, name):
    monkeypatch.setattr(k7, "MIN_CELLS", 20)
    t = k15.TASKS[name]
    ctx = t9._ctx(t, world[name])
    k15._PREP.clear()
    k10._PREP.clear()
    with harness.guarded_reads() as g:
        out = {arm: k15._build(ctx, arm) for arm in k15.ARMS}
        theirs = k10._build(ctx, "s1x")
    assert not g["refused"]
    assert np.array_equal(out["s1x"].X, theirs.X) and np.array_equal(out["s1x"].coords, theirs.coords)
    assert not np.array_equal(out["s1x_half"].X, out["s1x"].X)
    assert not np.array_equal(out["s1x_half"].X, out["s1x_half_perm"].X)
    prep = k15._PREP["fold"]
    assert np.array_equal(prep["phi_src"], k10._PREP["fold"]["phi_src"])
    assert np.array_equal(prep["ratio_src"], k10._PREP["fold"]["ratio_src"])
    pool = world[name]["wt"]
    want = harness.t3_decile_permutation(k15.half_factor(prep["phi_src"]), harness.t3_mean_expression(pool),
                                         k9.PERM_SEED)
    assert np.array_equal(prep["phi_perm"], want)
    rec = prep["record"]
    assert json.loads(json.dumps(rec)) == rec and rec["labels_with_a_depth_ratio"] >= 1
    assert rec["share_beyond_1.5_fold"]["half"] <= rec["share_beyond_1.5_fold"]["full"]
    # each arm's counts per cell after the build, recorded without changing a value (s1x equals K10's above)
    depth = prep["depth"]
    assert set(depth) == set(k15.ARMS) and all(set(v) == {str(ctx.seed)} for v in depth.values())
    assert json.loads(json.dumps(depth)) == depth
    labels = set(prep["ref"].V.labels)
    for v in depth.values():
        one = v[str(ctx.seed)]
        assert one["overall"] > 0 and set(one["per_label"]) <= labels and one["per_label"]
    assert depth["s1x_half"] != depth["s1x"]
    k15._PREP.clear()
    k10._PREP.clear()


def test_the_controls_decile_basis_is_the_references_not_the_source_wt(world, monkeypatch):
    """should_fix: the source WT (wt_train) and the reference role get different objects with different mean
    expression; the control's basis is the reference's."""
    monkeypatch.setattr(k7, "MIN_CELLS", 20)
    t = k15.TASKS["r12"]
    inputs = dict(world["r12"])
    ref_role = harness.TaskView.of(t).reference_role
    assert ref_role != "wt_train"
    other = inputs["wt_train"].copy()
    other.X = other.X[:, ::-1].copy()  # the same genes and cells, another gene profile in .X
    inputs["wt_train"] = other
    pool = inputs[ref_role]
    b_ref, b_src = harness.t3_mean_expression(pool), harness.t3_mean_expression(other)
    assert not np.allclose(b_ref, b_src)
    k15._PREP.clear()
    prep = k15._prepared(t9._ctx(t, inputs))
    half = k15.half_factor(prep["phi_src"])
    assert np.array_equal(prep["phi_perm"], harness.t3_decile_permutation(half, b_ref, k9.PERM_SEED))
    assert not np.array_equal(prep["phi_perm"], harness.t3_decile_permutation(half, b_src, k9.PERM_SEED))
    k15._PREP.clear()


# ----------------------------------------------------------------------------------------------- the reading
def _rows(scores: dict[str, list[float]], gated: tuple = ()):
    return [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5 + (0.01 if a == k15.HALF else 0.0) - (0.02 if a == k15.PERM else 0.0),
                            "mmd_u": 0.4} for a in scores},
             "t3": {a: {"severity": {"state": "r2_gate" if (a, i) in gated else "scored"}} for a in scores}}
            for i, s in enumerate(SEEDS)]


def _result(name, scores, clean=True, commit="abc", source=k15.TASK_SOURCE, gated=()):
    return {"task": {"key": k15.TASKS[name].key}, "rows": _rows(scores, gated), "task_source": source,
            "provenance": {"commit": commit, "tree_clean": clean},
            "factor_record": {"labels_with_a_depth_ratio": 5, "share_beyond_1.5_fold": {"full": 0.2, "half": 0.05}}}


def _scores(half, perm=50.0, base=52.0):
    n = len(SEEDS)
    h = list(half) if isinstance(half, list) else [half] * n
    return {harness.FLOOR: [50.0] * n, k15.BASE: [base] * n, k15.HALF: h, k15.PERM: [perm] * n}


def _k10_file(name: str, base: float = 52.0, **over) -> dict:
    """A synthetic file in try K10's committed form (its task_source, commit, clean tree, key and arms)."""
    r = {"task": {"key": k15.TASKS[name].key}, "task_source": k15.K10_SOURCE,
         "provenance": {"commit": k15.K10_COMMIT, "tree_clean": True},
         "rows": [{"seed": s, "scores": {"floor": 50.0, "s1x": base, "s1x_res": 0.0, "s1x_res_perm": 0.0}}
                  for s in SEEDS]}
    for k, v in over.items():
        r[k] = v
    return r


def _k10_dir(d: Path, base: float = 52.0, files: dict | None = None, pin: bool = True) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    for name in k15.TASKS:
        blob = json.dumps((files or {}).get(name) or _k10_file(name, base)).encode()
        (d / f"{name}.json").write_bytes(blob)
        if pin:
            PINS[name] = hashlib.sha256(blob).hexdigest()
    return d


def _read(d: Path, k10_dir=None, results=None, **per_task):
    d.mkdir()
    for name, sc in per_task.items():
        (d / f"{name}.json").write_text(json.dumps(_result(name, sc) if results is None else results[name]))
    return k15.reading(d, _k10_dir(d.parent / f"{d.name}_k10") if k10_dir is None else k10_dir)


def test_the_reading_on_the_replicate_folds_says_not_readable_on_the_cross_knockout_folds(tmp_path):
    res = _read(tmp_path / "a", r12=_scores(54.0), r21=_scores([53.0] * 9 + [51.0]))
    assert not res["void"] and res["supported"] is None and res["checks"] == {} and res["failing_folds"] == []
    assert res["verdict"].startswith("NOT READABLE ON THE CROSS-KNOCKOUT FOLDS") and "not on record" in res["verdict"]
    assert res["k10_reproduction"] == "reproduced" and res["not_readable_folds"] == []
    r12, r21 = res["tasks"]["r12"], res["tasks"]["r21"]
    assert r12["readable"] and not r12["counted_fold"]
    d = r12["paired"]["s1x_half - s1x"]
    assert (d["mean"], d["ahead"], d["n"], d["min"], d["max"], d["passes"]) == (
        pytest.approx(2.0), 10, 10, pytest.approx(2.0), pytest.approx(2.0), True)
    assert r12["paired"]["s1x_half - s1x_half_perm"]["mean"] == pytest.approx(4.0)
    assert r12["paired"]["s1x_half_perm - s1x"]["mean"] == pytest.approx(-2.0)
    e = r21["paired"]["s1x_half - s1x"]
    assert e["mean"] == pytest.approx(0.8) and e["ahead"] == 9 and e["min"] == pytest.approx(-1.0)
    assert e["max"] == pytest.approx(1.0) and e["passes"]
    assert r12["mean_score"] == pytest.approx({"floor": 50.0, "s1x": 52.0, "s1x_half": 54.0, "s1x_half_perm": 50.0})
    assert r12["skill_diff"]["s1x_half - s1x"] == pytest.approx({"de_score": 0.01, "mmd_u": 0.0})
    assert r12["skill_diff"]["s1x_half - s1x_half_perm"] == pytest.approx({"de_score": 0.03, "mmd_u": 0.0})
    assert r12["severity_states"]["s1x_half"] == {"scored": 10}
    assert r12["record"]["share_beyond_1.5_fold"] == {"full": 0.2, "half": 0.05}
    txt = k15.text(res)
    for phrase in ("NOT READABLE ON THE CROSS-KNOCKOUT FOLDS", "NOT the adoption rule", "changes no file",
                   "captain's own ruling", "without the file's composition term", "two samples on one slide",
                   "no released fold has the test's structure", "which fold the test resembles is not known",
                   "differ in depth against their references", "one factor, ten seeds of draw noise",
                   "per-seed range", "half 0.050, full 0.200", "replicate fold, never counted",
                   "Replicate-only", "prices shrinking the source factor", "identity check",
                   "Neither sign by itself can support or reject the half factor at the test",
                   "whole effect in the file's form",
                   "draw noise only", "mirrors of one slide, not two independent replications",
                   "Run as: DGX ~/v/vec-campaign", "CPU only, 8 BLAS threads, one process per fold",
                   "not below the guard margin -0.3", f"commit {k15.K10_COMMIT}"):
        assert phrase in txt, phrase
    assert r12["half_vs_s1x_within_guard_margin"] and res["run_as"] == k15.RUN_AS
    cmp = res["k10_comparator"]
    assert cmp["commit"] == k15.K10_COMMIT and set(cmp["files"]) == {"r12", "r21"}
    assert all(v["sha256"] == PINS[n] and v["commit"] == k15.K10_COMMIT for n, v in cmp["files"].items())
    for v in cmp["files"].values():
        assert v["sha256"] in txt


def test_a_fold_with_a_gated_seed_is_reported_not_readable(tmp_path):
    results = {"r12": _result("r12", _scores(54.0), gated=((k15.PERM, 3),)), "r21": _result("r21", _scores(54.0))}
    res = _read(tmp_path / "g", results=results, r12=None, r21=None)
    assert not res["void"] and res["not_readable_folds"] == ["r12"]
    assert not res["tasks"]["r12"]["readable"] and res["tasks"]["r21"]["readable"]
    assert "r12 (replicate fold, never counted), NOT READABLE" in k15.text(res)
    # the floor's severity state does not decide readability
    results = {"r12": _result("r12", _scores(54.0), gated=((harness.FLOOR, 0),)), "r21": _result("r21", _scores(54))}
    assert _read(tmp_path / "f", results=results, r12=None, r21=None)["not_readable_folds"] == []


def _fold(half_base: list[float], half_perm: float = 4.0, readable: bool = True) -> dict:
    # half_base is the per-seed s1x_half - s1x; the control sits half_perm below s1x_half's mean
    rows = _rows({harness.FLOOR: [0.0] * len(SEEDS), k15.BASE: [0.0] * len(SEEDS), k15.HALF: list(half_base),
                  k15.PERM: [float(np.mean(half_base)) - half_perm] * len(SEEDS)})
    f = k15.fold_figures(rows, counted=False)
    f["readable"] = readable
    return f


def _per(**over) -> dict:
    good = [1.0] * 10
    per = {"x1": _fold(good), "r12": _fold(good), "x2": _fold([-0.2] * 10), "r21": _fold([-0.3] * 10)}
    return per | over


def test_decide_supported_and_each_failing_fold_named():
    ok, checks, fails = k15.decide(_per())
    assert ok and fails == [] and all(checks.values()) and len(checks) == 6
    assert k15.verdict(_per())["verdict"] == "SUPPORTED IN THE FILE'S OWN FORM"
    cases = {
        "x1 against s1x": (_per(x1=_fold([0.3] * 10)), ["x1"]),
        "r12 against the control": (_per(r12=_fold([1.0] * 10, half_perm=0.2)), ["r12"]),
        "x1 seven seeds ahead": (_per(x1=_fold([2.0] * 7 + [-0.1] * 3)), ["x1"]),
        "x2 below -0.3": (_per(x2=_fold([-0.31] * 10)), ["x2"]),
        "r21 below -0.3": (_per(r21=_fold([-0.5] * 10)), ["r21"]),
        "r12 not readable": (_per(r12=_fold([1.0] * 10, readable=False)), ["r12"]),
        "x2 not readable": (_per(x2=_fold([1.0] * 10, readable=False)), ["x2"]),
        "two folds": (_per(x1=_fold([0.0] * 10), r21=_fold([-1.0] * 10)), ["x1", "r21"]),
    }
    for what, (per, want) in cases.items():
        ok, _, fails = k15.decide(per)
        assert not ok and fails == want, what
        v = k15.verdict(per)
        assert v["supported"] is False and v["verdict"] == "NOT SUPPORTED (failing folds: " + ", ".join(want) + ")"
    # eight seeds ahead with a mean above 0.3 passes
    assert k15.decide(_per(x1=_fold([1.0] * 8 + [-0.1] * 2)))[0]
    # a guarded fold whose paired mean is exactly FLOOR_MEAN in floating point is not below it; the next float down is
    for n in k15.GUARDED:
        per = _per()
        per[n]["paired"]["s1x_half - s1x"]["mean"] = k15.FLOOR_MEAN
        assert k15.FLOOR_MEAN == -0.3 and k15.decide(per)[0], n
        per[n]["paired"]["s1x_half - s1x"]["mean"] = float(np.nextafter(k15.FLOOR_MEAN, -1.0))
        ok, _, fails = k15.decide(per)
        assert not ok and fails == [n]
    # an absent cross-knockout fold: not readable there, never decided
    per = _per()
    del per["x2"]
    assert k15.verdict(per)["supported"] is None


def test_the_reading_is_void_on_missing_dirty_mixed_or_foreign_results(tmp_path):
    good = {"r12": _scores(54.0), "r21": _scores(54.0)}
    assert any("missing task r21" in w for w in _read(tmp_path / "m", r12=good["r12"])["why"])
    for what, results in {
        "dirty": {"r12": _result("r12", good["r12"]), "r21": _result("r21", good["r21"], clean=False)},
        "commits": {"r12": _result("r12", good["r12"]), "r21": _result("r21", good["r21"], commit="def")},
        "another script": {"r12": _result("r12", good["r12"], source="experiments/t3_decomp/k10_s1_form.py"),
                           "r21": _result("r21", good["r21"])},
    }.items():
        res = _read(tmp_path / what.replace(" ", "_"), results=results, r12=None, r21=None)
        assert res["void"], what
    no_arm = _result("r12", good["r12"])
    for q in no_arm["rows"][4:5]:
        del q["scores"][k15.PERM]
    res = _read(tmp_path / "arm", results={"r12": no_arm, "r21": _result("r21", good["r21"])}, r12=None, r21=None)
    assert res["void"] and any("missing arms ['s1x_half_perm']" in w for w in res["why"])
    no_seed = _result("r12", good["r12"])
    no_seed["rows"] = no_seed["rows"][:-1]
    assert _read(tmp_path / "seed", results={"r12": no_seed, "r21": _result("r21", good["r21"])},
                 r12=None, r21=None)["void"]
    other_key = _result("r12", good["r12"])
    other_key["task"]["key"] = k9.TASKS["r12"].key
    assert _read(tmp_path / "key", results={"r12": other_key, "r21": _result("r21", good["r21"])},
                 r12=None, r21=None)["void"]


def test_the_reading_is_void_unless_s1x_reproduces_k10s_committed_scores(tmp_path):
    good = {"r12": _scores(54.0), "r21": _scores(54.0)}
    absent = _read(tmp_path / "a", k10_dir=tmp_path / "nowhere", **good)
    assert absent["void"] and any("absent from this checkout" in w for w in absent["why"])
    off = _read(tmp_path / "o", k10_dir=_k10_dir(tmp_path / "k10_off", base=52.0 + 2e-6), **good)
    assert off["void"] and any("does not reproduce" in w for w in off["why"])
    near = _read(tmp_path / "n", k10_dir=_k10_dir(tmp_path / "k10_near", base=52.0 + 5e-7), **good)
    assert not near["void"]
    partial = _k10_dir(tmp_path / "k10_part")
    (partial / "r21.json").unlink()
    res = _read(tmp_path / "p", k10_dir=partial, **good)
    assert res["void"] and any("r21: try K10's committed result is missing" in w for w in res["why"])


def test_the_comparator_must_be_k10s_own_file(tmp_path):
    """must_fix 2: a self-compare, try K14's or K9's file, another commit, a dirty tree, another key, try K10's arms
    missing, or a file whose bytes are not the pinned ones each void the reading, even when s1x's scores match."""
    good = {"r12": _scores(54.0), "r21": _scores(54.0)}
    runs = tmp_path / "self"
    # self-compare: K15's own task files as the comparator (the same keys, s1x at the same scores)
    runs.mkdir()
    for name, sc in good.items():
        (runs / f"{name}.json").write_text(json.dumps(_result(name, sc)))
    res = k15.reading(runs, runs)
    assert res["void"] and any("not try K10's committed task file" in w for w in res["why"])
    assert res["k10_comparator"]["dir"] == str(runs)
    k14_file = _k10_file("r12", task_source="experiments/t3_decomp/k14_form_arms.py",
                         task={"key": k9.TASKS["r12"].key},
                         rows=[{"seed": s, "scores": {"raw": 52.0, "raw_half": 1.0, "s1x": 52.0}} for s in SEEDS])
    k9_file = _k10_file("r12", task_source="experiments/t3_decomp/k9_residual.py")
    wrong_commit = _k10_file("r12", provenance={"commit": "b45cd928616d642d855ad57b3e1bb159fd937094",
                                                "tree_clean": True})
    dirty = _k10_file("r12", provenance={"commit": k15.K10_COMMIT, "tree_clean": False})
    other_key = _k10_file("r12", task={"key": k15.TASKS["r21"].key})
    no_res = _k10_file("r12", rows=[{"seed": s, "scores": {"s1x": 52.0}} for s in SEEDS])
    for what, f in {"k14": k14_file, "k9": k9_file, "commit": wrong_commit, "dirty": dirty, "key": other_key,
                    "no s1x_res": no_res}.items():
        d = _k10_dir(tmp_path / f"cmp_{what.replace(' ', '_')}", files={"r12": f})
        res = _read(tmp_path / f"run_{what.replace(' ', '_')}", k10_dir=d, **good)
        assert res["void"] and any("r12: the comparator" in w and "not try K10's" in w for w in res["why"]), what
        assert not any(w.startswith("r21") for w in res["why"]), what
    # right fields, but bytes that are not the pinned file's
    d = _k10_dir(tmp_path / "cmp_unpinned", pin=False)
    res = _read(tmp_path / "run_unpinned", k10_dir=d, **good)
    assert res["void"] and any("sha256 is not that of try K10's committed file" in w for w in res["why"])
    # a comparator that is not JSON voids, never raises
    d = _k10_dir(tmp_path / "cmp_garbage")
    (d / "r12.json").write_text("not json")
    PINS["r12"] = hashlib.sha256(b"not json").hexdigest()
    res = _read(tmp_path / "run_garbage", k10_dir=d, **good)
    assert res["void"] and any("r12: the comparator" in w for w in res["why"])


def test_a_row_without_skills_or_severity_voids_and_does_not_raise(tmp_path):
    good = {"r12": _scores(54.0), "r21": _scores(54.0)}
    for what in ("skills", "t3", "empty skills", "severity"):
        r = _result("r12", good["r12"])
        q = r["rows"][3]
        if what == "skills":
            del q["skills"][k15.PERM]
        elif what == "t3":
            del q["t3"][k15.BASE]
        elif what == "empty skills":
            q["skills"][k15.HALF] = {}
        else:
            del q["t3"][k15.HALF]["severity"]
        res = _read(tmp_path / what.replace(" ", "_"), results={"r12": r, "r21": _result("r21", good["r21"])},
                    r12=None, r21=None)
        assert res["void"] and any("r12: a row lacks an arm's skills or T3 severity state" in w
                                   for w in res["why"]), what


def test_cmd_run_records_the_machine_and_each_arms_depth(tmp_path, monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "8")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "8")
    monkeypatch.delenv("MKL_NUM_THREADS", raising=False)
    depth = {"s1x": {"7": {"overall": 1.0, "per_label": {"a": 1.0}}}}

    def fake_run(t, cands, seeds, data_dir):
        k15._PREP["fold"] = {"record": {"labels_with_a_depth_ratio": 1}, "depth": depth}
        return {"task": {"key": t.key}, "rows": []}

    monkeypatch.setattr(k15.harness, "run_task", fake_run)
    k15.main(["run", "--task", "r12", "--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])
    r = json.loads((tmp_path / "r12.json").read_text())
    assert r["depth_record"] == depth and r["task_source"] == k15.TASK_SOURCE
    env = r["run_env"]
    assert set(env) == {"host", *k15.RUN_ENV_VARS} and env["OMP_NUM_THREADS"] == "8" and env["MKL_NUM_THREADS"] is None
    assert env["host"] == os.uname().nodename
    k15._PREP.clear()


def test_the_reading_writes_no_per_gene_value(tmp_path):
    res = _read(tmp_path / "w", r12=_scores(54.0), r21=_scores(53.0))
    blob = json.dumps(res) + k15.text(res)
    assert not any(g in re.findall(r"[A-Za-z0-9_]+", blob) for g in t9.GENES)
    assert "phi_src" not in blob and "phi_perm" not in blob


def test_run_and_read_refuse(tmp_path, monkeypatch):
    (tmp_path / "r12.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k15.main(["run", "--task", "r12", "--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])
    for task in ("x1", "x2"):
        with pytest.raises(SystemExit, match="design not on record"):
            k15.main(["run", "--task", task, "--out", str(tmp_path / "x"), "--data", str(tmp_path / "nowhere")])
        assert not (tmp_path / "x" / f"{task}.json").exists()
    runs = tmp_path / "runs"
    _read(runs, r12=_scores(54.0), r21=_scores(53.0))
    k10_dir = tmp_path / "runs_k10"
    out = tmp_path / "out"
    with pytest.raises(SystemExit):  # no operator path for the comparator
        k15.main(["read", "--runs", str(runs), "--out", str(out), "--k10", str(k10_dir)])
    assert not out.exists()
    monkeypatch.setattr(k15, "RESULTS_K10", k10_dir)  # the CLI reads the module's committed directory only
    k15.main(["read", "--runs", str(runs), "--out", str(out)])
    got = json.loads((out / "reading.json").read_text())
    assert got["verdict"].startswith("NOT READABLE") and got["k10_comparator"]["dir"] == str(k10_dir)
    with pytest.raises(SystemExit, match="exists"):
        k15.main(["read", "--runs", str(runs), "--out", str(out)])
    (out / "reading.json").unlink()
    with pytest.raises(SystemExit, match="exists"):  # reading.txt alone also refuses
        k15.main(["read", "--runs", str(runs), "--out", str(out)])


# ----------------------------------------------------------------------------------------------- review fixes (fix2)
def _with_depth(name, sc, half_ratio=0.95):
    r = _result(name, sc)
    r["depth_record"] = {"s1x": {"1": {"overall": 100.0, "per_label": {"a": 100.0}}},
                         "s1x_half": {"1": {"overall": 100.0 * half_ratio, "per_label": {"a": 100.0 * half_ratio}}},
                         "s1x_half_perm": {"1": {"overall": 99.0, "per_label": {"a": 99.0}}}}
    return r


@pytest.mark.parametrize("half", [54.0, 50.0])
def test_the_replicate_sign_is_measured_never_by_construction(tmp_path, half):
    """must_fix: the replicate folds' s1x_half - s1x sign is a measurement of the source factor's noise against its
    signal; no text says it is fixed by construction, whichever sign comes out, and nothing voids on either sign."""
    res = _read(tmp_path / "s", r12=_scores(half), r21=_scores(half))
    assert not res["void"] and res["supported"] is None
    txt = k15.text(res)
    blob = txt + json.dumps(res) + k15.__doc__
    for gone in ("by construction", "right factor", "fails whatever the test's truth is", "loss price"):
        assert gone not in blob, gone
    for phrase in ("Its sign is measured, not given", "shrinkage wins where the factor's per-gene estimation noise "
                   "exceeds its signal", "Try K14b's factor-alone r12 read was -2.28", "it decides nothing",
                   "changes the advice on the captain's half-factor question", "its sign measured, not given"):
        assert phrase in txt, phrase
    assert res["tasks"]["r12"]["half_vs_s1x_within_guard_margin"] is (half - 52.0 >= k15.FLOOR_MEAN)
    assert "shrink_cost_within_guard" not in json.dumps(res)


def test_the_docstring_says_why_r12_is_not_a_transfer_read_and_keeps_the_rule():
    """design attack: the r12 condition is reworded as the noise-versus-signal hurdle it is; the bars are unchanged."""
    doc = " ".join(k15.__doc__.split())
    assert "the r12 condition does NOT" in doc and "noise-versus-signal hurdle" in doc
    assert "cannot be reached on this commit because x1 and x2 are not built" in doc
    assert "byte copies of main's commit 3fea49d" in doc
    assert k15.DECIDING == ("x1", "r12") and k15.GUARDED == ("x2", "r21")
    assert (k15.MIN_GAIN, k15.MIN_AHEAD, k15.FLOOR_MEAN) == (0.3, 8, -0.3)


def test_a_failed_identity_check_is_reported_never_voided(tmp_path):
    """should_fix: the control beating the half arm on a replicate fold is reported as a failure, not a void."""
    res = _read(tmp_path / "i", r12=_scores(54.0, perm=56.0), r21=_scores(54.0, perm=56.0))
    assert not res["void"]
    d = res["tasks"]["r12"]["paired"]["s1x_half - s1x_half_perm"]
    assert d["mean"] == pytest.approx(-2.0) and not d["passes"]
    txt = k15.text(res)
    assert "a pass is expected where the factor carries gene identity, and a failure is reported as such, never as a " \
           "void" in txt
    assert "identity check: which gene carries which half value, not size; a failure is reported as such" in txt


def test_the_limits_name_the_uncoupled_draws_and_the_missing_size_control(tmp_path):
    res = _read(tmp_path / "l", r12=_scores(54.0), r21=_scores(54.0))
    txt = " ".join(k15.text(res).split())
    assert "the arms share their drawn cells, but the factor's and the depth step's own draws are not coupled across " \
           "arms, so the per-seed range includes that noise" in txt
    assert "no arm prices size apart from depth" in txt


def test_the_depth_shift_is_printed_next_to_the_half_difference(tmp_path):
    """should_fix: s1x_half's counts per cell against s1x's sit on the s1x_half - s1x line."""
    res = _read(tmp_path / "d", results={"r12": _with_depth("r12", _scores(54.0)),
                                         "r21": _with_depth("r21", _scores(54.0), 0.909)}, r12=None, r21=None)
    lines = k15.text(res).splitlines()
    got = [ln for ln in lines if ln.strip().startswith("s1x_half - s1x [")]
    assert len(got) == 2
    assert got[0].endswith("s1x_half's mean counts per cell after the build: 0.9500 of s1x's")
    assert got[1].endswith("s1x_half's mean counts per cell after the build: 0.9090 of s1x's")
    assert not any("counts per cell" in ln for ln in lines if ln.strip().startswith("s1x_half - s1x_half_perm ["))
    res = _read(tmp_path / "n", r12=_scores(54.0), r21=_scores(54.0))  # no depth record: said so, no raise
    assert "counts per cell after the build: not recorded" in k15.text(res)


@pytest.mark.parametrize("change", ["missing", "extra"])
def test_a_k10_comparator_with_other_seeds_voids_and_does_not_raise(tmp_path, change):
    """should_fix: a pinned K10 file whose rows miss a seed or carry an extra one voids with 'does not reproduce'."""
    f = _k10_file("r12")
    if change == "missing":
        f["rows"] = f["rows"][:-1]
    else:
        f["rows"].append({"seed": max(SEEDS) + 1, "scores": {"floor": 50.0, "s1x": 52.0, "s1x_res": 0.0,
                                                               "s1x_res_perm": 0.0}})
    d = _k10_dir(tmp_path / f"cmp_{change}", files={"r12": f})
    res = _read(tmp_path / f"run_{change}", k10_dir=d, r12=_scores(54.0), r21=_scores(54.0))
    assert res["void"] and any(w.startswith("r12:") and "does not reproduce" in w for w in res["why"])
    assert not any(w.startswith("r21") for w in res["why"])
