"""Tests for tries K10a and K10b (experiments/t3_decomp/k10_s1_form.py), on synthetic data only: no competition file is
opened, nothing is scored, and no prediction is compared with any target."""

from __future__ import annotations

import json
import re
import sys
from dataclasses import fields
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from vec import harness
from vec.adoption import SEEDS
from vec.prep import load_stage

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import k7_wt_revert as k7  # noqa: E402
import k9_residual as k9  # noqa: E402
import k10_s1_form as k10  # noqa: E402
from test_k9_residual import E875, GENES, _cells, _ctx, _label, _reference  # noqa: E402

SOURCE = HERE / "k10_s1_form.py"
G = 16
CLIP = (0.25, 4.0)
REP1, REP2 = "T3/E8.75_gata4_ko_1.h5ad", "T3/E8.75_gata4_ko_2.h5ad"


@pytest.fixture(scope="module")
def bs():
    return harness._t3_head()


# ----------------------------------------------------------------------------------------------- constants and tasks
def test_the_fixed_constants():
    assert k10.RES_EXPONENT == k9.RES_EXPONENT == 0.5
    assert k10.PERM_SEED == k9.PERM_SEED == 20261113 and k10.STREAM == k9.STREAM == harness.K4_STREAM
    assert (k10.MIN_GAIN, k10.MIN_AHEAD) == (0.3, 8)
    assert k10.ARMS == ("s1x", "s1x_res", "s1x_res_perm")
    assert k10.SOURCE == harness.K4_FOLDS["T3/Gata4rep1->rep2"]
    assert (k10.SOURCE["ko_roles"], k10.SOURCE["source_preset"], k10.SOURCE["reference_preset"]) == (
        ("ko_train",), "gata4", "gata4")


def test_r12_is_k9s_r12_and_r21_its_mirror():
    r12, r21, k9r = k10.TASKS["r12"], k10.TASKS["r21"], k9.TASKS["r12"]
    assert set(k10.TASKS) == {"r12", "r21"} and r12.key != k9r.key and r21.key != r12.key
    for f in fields(r12):
        if f.name not in ("key", "notes"):
            assert getattr(r12, f.name) == getattr(k9r, f.name), f.name
    assert r12.target == REP2 and r12.inputs["ko_train"] == REP1 and r12.reference == "T2/E8.75.h5ad"
    assert r21.target == REP1 and r21.inputs["ko_train"] == REP2
    assert {r: f for r, f in r21.inputs.items() if r != "ko_train"} == {r: f for r, f in r12.inputs.items()
                                                                         if r != "ko_train"}
    for f in fields(r12):
        if f.name not in ("key", "notes", "target", "inputs", "forbidden", "validation_board", "twin_board"):
            assert getattr(r21, f.name) == getattr(r12, f.name), f.name
    assert r12.validation_board == "T3:gata4" and r21.validation_board is None and r21.twin_board == "T3:gata4"
    swap = {REP1: REP2, REP2: REP1}
    assert set(r21.forbidden) == {swap.get(f, f) for f in r12.forbidden} | {REP1}
    for t in k10.TASKS.values():
        assert harness.check_task(t) is t
        assert t.target not in t.inputs.values() and t.replicate and not t.counted and not t.adoption
        assert t.test_board == "T3:bcatenin" and harness.TaskView.of(t).reference_role == "wt"
        assert {f for r, f in t.inputs.items() if r.startswith(k9.SERIES_PREFIX)} == set(
            k9.series_files(t.reference).values())
    assert REP1 in r21.forbidden
    assert all(c.applies(t) for c in k10.candidates()[1:] for t in k10.TASKS.values())
    assert not any(c.applies(t) for c in k10.candidates()[1:] for ts in harness.load_tasks().values() for t in ts)
    assert not any(c.applies(t) for c in k10.candidates()[1:] for t in k9.TASKS.values())


def test_the_source_opens_no_file_beyond_the_tasks_inputs():
    src = SOURCE.read_text()
    for call in ("read_h5ad", "load_stage", "read_stage", "h5py", "np.load", "open(", "boards/", "runs/", "data/"):
        assert call not in src, call
    assert re.findall(r'"[^"]*\.h5ad"', src) == []


# ----------------------------------------------------------------------------------------------- the arms
def _setup(bs, rng):
    ref = _reference(bs, rng)
    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, ref.pool.C.shape[0])]
    phi = np.clip(np.exp(rng.normal(0, 0.6, G)), *CLIP)
    ratio = np.array([0.7, 1.3, 1.0])
    pool = np.exp(rng.normal(0, 0.3, G))
    per = {"A": pool * np.exp(rng.normal(0, 0.8, G)), "B": pool * np.exp(rng.normal(0, 0.8, G))}
    return ref, group, phi, ratio, pool, per


def _spies(bs, monkeypatch):
    calls = []
    real_f, real_d = bs.t3d.apply_gene_factor, bs.t3d.rescale_depth

    def fac(C, phi, P, rng):
        calls.append(("factor", C.copy(), np.array(phi)))
        return real_f(C, phi, P, rng)

    def dep(C, ratio, P, rng):
        calls.append(("depth", C.copy(), np.array(ratio)))
        return real_d(C, ratio, P, rng)

    monkeypatch.setattr(bs.t3d, "apply_gene_factor", fac)
    monkeypatch.setattr(bs.t3d, "rescale_depth", dep)
    return calls


def test_s1x_is_construct_with_the_factor_and_the_depth_ratio(bs):
    rng = np.random.default_rng(4)
    ref, group, phi, ratio, _, _ = _setup(bs, rng)
    got = k10.arm_built(bs, ref, 200, 7, ref.V.p, "s1x", phi, ratio, None, group, CLIP)
    want = bs.construct(ref, 200, bs.t3d.stable_rng(7, k10.STREAM), f=ref.V.p, phi=phi, ratio=ratio)
    assert np.array_equal(got.X, want.X) and np.array_equal(got.row, want.row) and np.array_equal(got.lab, want.lab)
    no_depth = bs.construct(ref, 200, bs.t3d.stable_rng(7, k10.STREAM), f=ref.V.p, phi=phi)
    assert not np.array_equal(got.X, no_depth.X)  # the depth step moves the cells


def test_every_arm_draws_the_same_cells_and_applies_each_groups_factor_then_the_depth(bs, monkeypatch):
    rng = np.random.default_rng(5)
    ref, group, phi, ratio, pool, per = _setup(bs, rng)
    rho = k9.residual_factors(per, pool, k10.RES_EXPONENT)
    rho_perm = k9.permuted_residuals(rho, rng.gamma(2.0, size=G), k10.PERM_SEED)
    calls = _spies(bs, monkeypatch)
    built = {arm: k10.arm_built(bs, ref, 200, 7, ref.V.p, arm, phi, ratio, r, group, CLIP)
             for arm, r in (("s1x", None), ("s1x_res", rho), ("s1x_res_perm", rho_perm))}
    rows, lab = built["s1x"].row, built["s1x"].lab
    assert all(np.array_equal(b.row, rows) and np.array_equal(b.lab, lab) for b in built.values())
    assert not np.array_equal(built["s1x"].X, built["s1x_res"].X)
    assert not np.array_equal(built["s1x_res"].X, built["s1x_res_perm"].X)
    # s1x: one factor call and one depth call; each residual arm: one factor call per group, then one depth call
    assert [c[0] for c in calls] == ["factor", "depth"] + (["factor"] * 3 + ["depth"]) * 2
    assert np.array_equal(calls[0][1], ref.pool.C[rows]) and np.array_equal(calls[0][2], phi)
    assert all(np.array_equal(calls[k][2], ratio[lab]) for k in (1, 5, 9))
    drawn = group[rows]
    for k, r in ((2, rho), (6, rho_perm)):
        for (_, C, fac), g in zip(calls[k:k + 3], ("A", "B", "Unknown"), strict=True):
            assert np.array_equal(C, ref.pool.C[rows][drawn == g])
            assert np.array_equal(fac, k9.combined_factor(phi, r.get(g), CLIP))
    again = k10.arm_built(bs, ref, 200, 7, ref.V.p, "s1x_res", phi, ratio, rho, group, CLIP)
    assert np.array_equal(again.X, built["s1x_res"].X)


@pytest.mark.parametrize("case", ["exponent_zero", "per_equals_pool"])
def test_without_a_residual_s1x_res_is_s1x_bit_for_bit(bs, monkeypatch, case):
    rng = np.random.default_rng(6)
    ref, group, phi, ratio, pool, per = _setup(bs, rng)
    if case == "exponent_zero":
        monkeypatch.setattr(k10, "RES_EXPONENT", 0.0)
    else:
        per = {"A": pool.copy(), "B": pool.copy()}
    rho = k9.residual_factors(per, pool, k10.RES_EXPONENT)
    rho_perm = k9.permuted_residuals(rho, rng.gamma(2.0, size=G), k10.PERM_SEED)
    base = k10.arm_built(bs, ref, 200, 9, ref.V.p, "s1x", phi, ratio, None, group, CLIP)
    for arm, r in (("s1x_res", rho), ("s1x_res_perm", rho_perm)):
        res = k10.arm_built(bs, ref, 200, 9, ref.V.p, arm, phi, ratio, r, group, CLIP)
        assert np.array_equal(res.X, base.X) and np.array_equal(res.row, base.row) and np.array_equal(res.lab, base.lab)


# ----------------------------------------------------------------------------------------------- end to end, synthetic
@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Synthetic files in r12's and r21's roles (WT E8.75, the two Gata4 replicates shallower than the WT and moved by
    a factor, the wild-type series), written and loaded as the harness loads them."""
    d = tmp_path_factory.mktemp("k10")
    rng = np.random.default_rng(0)
    prof = {k: rng.lognormal(0, 1.0, len(GENES)) for k in sorted({_label(c) for c in E875})}
    fac = rng.lognormal(0, 0.4, len(GENES))
    files = {"wt": _cells(E875, prof, 700, rng, 150, rng.lognormal(0, 0.3, len(GENES)), name="ref"),
             "rep1": _cells(E875, prof, 600, rng, 90, fac, name="rep1"),
             "rep2": _cells(E875, prof, 600, rng, 100, fac, name="rep2")}
    for r in k10.TASKS["r12"].inputs:
        if r.startswith(k9.SERIES_PREFIX):
            files[r] = _cells(E875, prof, 500, rng, 120, rng.lognormal(0, 0.1, len(GENES)), name=r)
    loaded = {}
    for r, a in files.items():
        p = d / f"{r}.h5ad"
        a.write_h5ad(p)
        loaded[r] = load_stage(p, genes=GENES, drop_duplicates=True)
    series = {r: a for r, a in loaded.items() if r.startswith(k9.SERIES_PREFIX)}
    wt_train = load_stage(d / "wt.h5ad", genes=GENES, drop_duplicates=True)
    return {"r12": {"wt": loaded["wt"], "wt_train": wt_train, "ko_train": loaded["rep1"], **series},
            "r21": {"wt": loaded["wt"], "wt_train": wt_train, "ko_train": loaded["rep2"], **series}}


def test_ratio_src_is_the_frozen_estimates_ratio(world, monkeypatch):
    """t3_decompose.estimate called directly on the source side (the frame's cells and labels, the WT in the same
    vocabulary; the composition bracket stubbed, no stage pools) gives the same ratio and the same factor."""
    t3d = harness._t3_head().t3d
    rules = harness._t3_head().load_rules()
    inputs = world["r12"]
    phi_src, ratio_src, labels = k10.source_side(inputs["wt_train"], [inputs["ko_train"]], rules, "gata4")
    frame = harness.t3_source_frame(inputs["wt_train"], [inputs["ko_train"]], rules, "gata4")
    w = harness._shallow(inputs["wt_train"])
    t3d.add_label_column(w, "combined")
    V, wlab = harness._t3_vocab(t3d, w, rules)
    Cw = t3d.dense_counts(w)
    wpool = t3d.Pool(Cw, wlab, None, V.L)
    Xw = t3d.lognorm(Cw).astype(np.float64)
    wt = SimpleNamespace(pool=wpool, X=Xw, M=t3d.type_means(Xw, wlab, V.L), mu=t3d.props_type_means(Cw, wlab, V.L),
                         med=t3d.label_medians(wlab, Cw.sum(1), V.L), count=wpool.count)
    q = V.p[V.model_idx] / V.p[V.model_idx].sum()
    monkeypatch.setattr(t3d, "bracket", lambda calls, model, rows=None: {e: q.copy() for e in t3d.EST})
    ctx = SimpleNamespace(V=V, wt=wt, rules=rules, harm=SimpleNamespace(model=lambda key: None),
                          p_oof={None: {e: q.copy() for e in t3d.EST}}, stages=[])
    T = SimpleNamespace(C=frame.C, lab=frame.lab, calls=None, key=None)
    est = t3d.estimate(ctx, T, np.arange(len(frame.lab)), t3d.stable_rng(0, "full:test"))
    assert labels == list(V.labels)
    assert np.array_equal(ratio_src, est.ratio) and np.array_equal(phi_src, est.phi_nm)
    assert np.array_equal(phi_src, frame.phi)
    assert (ratio_src != 1.0).any() and (ratio_src < 1.0).any()  # the shallower replicate gives ratios below 1
    dmin = rules["components"]["depth"]["min_cells"]
    small = (np.bincount(frame.lab[frame.lab >= 0], minlength=V.L) < dmin) | (wpool.count < dmin)
    assert np.all(ratio_src[small] == 1.0)


@pytest.mark.parametrize("name", ["r12", "r21"])
def test_end_to_end_s1x_is_the_frozen_stand_in_and_the_residual_moves_it(world, monkeypatch, bs, name):
    monkeypatch.setattr(k7, "MIN_CELLS", 20)
    t = k10.TASKS[name]
    k10._PREP.clear()
    ctx = _ctx(t, world[name])
    with harness.guarded_reads() as g:
        out = {arm: k10._build(ctx, arm) for arm in k10.ARMS}
    assert not g["refused"]
    prep = k10._PREP["fold"]
    rec = prep["record"]
    assert rec["groups_with_a_factor"] >= 2 and -1.0 <= rec["cosine_log_phi_raw_mean_log_rho"] <= 1.0
    assert rec["labels_with_a_depth_ratio"] >= 1
    assert json.loads(json.dumps(rec)) == rec
    rules = bs.load_rules()
    phi_src, ratio_src, _ = k10.source_side(world[name]["wt_train"], [world[name]["ko_train"]], rules, "gata4")
    assert np.array_equal(prep["phi_src"], phi_src) and np.array_equal(prep["ratio_src"], ratio_src)
    ref = harness.t3_reference(world[name]["wt"], "gata4", rules)
    n = harness.count_for(ctx.board, world[name]["wt"].n_obs)
    built = bs.construct(ref, n, bs.t3d.stable_rng(ctx.seed, k10.STREAM), f=ref.V.p, phi=phi_src, ratio=ratio_src)
    X, coords, _ = bs.assemble(ref, built, list(ctx.board.genes), bs.t3d.stable_rng(ctx.seed, f"{k10.STREAM}:assemble"))
    want = harness.finish(ctx, world[name]["wt"], X, coords)
    assert np.array_equal(out["s1x"].X, want.X) and np.array_equal(out["s1x"].coords, want.coords)
    assert not np.array_equal(out["s1x_res"].X, out["s1x"].X)
    k10._PREP.clear()
    monkeypatch.setattr(k10, "RES_EXPONENT", 0.0)
    with harness.guarded_reads():
        none = k10._build(ctx, "s1x_res")
    assert np.array_equal(none.X, out["s1x"].X) and np.array_equal(none.coords, out["s1x"].coords)
    k10._PREP.clear()


def test_the_two_folds_take_their_own_source(world, monkeypatch):
    rules = harness._t3_head().load_rules()
    a = k10.source_side(world["r12"]["wt_train"], [world["r12"]["ko_train"]], rules, "gata4")
    b = k10.source_side(world["r21"]["wt_train"], [world["r21"]["ko_train"]], rules, "gata4")
    assert not np.array_equal(a[1], b[1]) and a[2] == b[2]


# ----------------------------------------------------------------------------------------------- the reading
def _result(key: str, scores: dict[str, list[float]], clean: bool = True, commit: str = "abc") -> dict:
    rows = [{"seed": s, "scores": {a: v[i] for a, v in scores.items()},
             "skills": {a: {"de_score": 0.5 + (0.01 if a == k10.MAIN else 0.0), "severity": 0.6} for a in scores},
             "t3": {a: {"severity": {"state": "scored" if a == k10.MAIN else "r2_gate"}} for a in scores}}
            for i, s in enumerate(SEEDS)]
    return {"task": {"key": key}, "rows": rows, "provenance": {"commit": commit, "tree_clean": clean},
            "factor_record": {"groups_with_a_factor": 9, "labels_with_a_depth_ratio": 5,
                              "cosine_log_phi_raw_mean_log_rho": 0.4}}


def _scores(main, control: float = 50.0, base: float = 52.0) -> dict[str, list[float]]:
    n = len(SEEDS)
    m = list(main) if isinstance(main, list) else [main] * n
    return {harness.FLOOR: [50.0] * n, k10.BASE: [base] * n, k10.MAIN: m, k10.CONTROL: [control] * n}


def _read(d: Path, **per_task) -> dict:
    d.mkdir()
    for name, sc in per_task.items():
        (d / f"{name}.json").write_text(json.dumps(_result(k10.TASKS[name].key, sc)))
    return k10.reading(d)


def test_the_reading_needs_both_folds_against_both_the_base_and_the_control(tmp_path):
    yes = _read(tmp_path / "yes", r12=_scores(54.0), r21=_scores(53.0))
    assert not yes["void"] and yes["k10_confirmed"] and yes["verdict"].startswith("CONFIRMED")
    assert yes["tasks"]["r12"]["paired"]["s1x_res - s1x"]["mean"] == pytest.approx(2.0)
    assert yes["tasks"]["r12"]["paired"]["s1x - floor"]["mean"] == pytest.approx(2.0)
    assert yes["tasks"]["r12"]["skill_diff_s1x_res_minus_s1x"] == pytest.approx({"de_score": 0.01, "severity": 0.0})
    assert yes["tasks"]["r21"]["severity_states"][k10.MAIN] == {"scored": len(SEEDS)}
    assert not yes["tasks"]["r12"]["counted_fold"] and not yes["tasks"]["r21"]["counted_fold"]
    txt = k10.text(yes)
    for phrase in ("replicate folds on one slide", "never counted by the adoption rule", "frozen stand-in",
                   "factor and depth from the other replicate", "composition term", "fixed by the frozen rules",
                   "magnitude was set from both replicates", "one priced request"):
        assert phrase in txt, phrase
    # either fold failing against either comparison: not confirmed
    assert not _read(tmp_path / "r21b", r12=_scores(54.0), r21=_scores(52.2))["k10_confirmed"]
    no = _read(tmp_path / "r12c", r12=_scores(54.0, control=53.9), r21=_scores(53.0))
    assert not no["k10_confirmed"] and no["verdict"].startswith("NOT CONFIRMED")
    assert "not carried to the first-choice file's form" in k10.text(no)
    assert not _read(tmp_path / "r21c", r12=_scores(54.0), r21=_scores(53.0, control=53.0))["k10_confirmed"]
    assert not _read(tmp_path / "r12b", r12=_scores(52.0), r21=_scores(60.0))["k10_confirmed"]


def test_the_gain_and_seed_boundaries(tmp_path):
    n = len(SEEDS)
    # a mean gain of exactly 0.3 does not pass; just above does (base 52.0, control far below)
    assert not _read(tmp_path / "g0", r12=_scores(52.3), r21=_scores(54.0))["k10_confirmed"]
    assert _read(tmp_path / "g1", r12=_scores(52.31), r21=_scores(54.0))["k10_confirmed"]
    # ahead in 8 of 10 seeds passes, in 7 does not (means well above 0.3 in both)
    eight = [60.0] * 8 + [51.0] * (n - 8)
    seven = [60.0] * 7 + [51.0] * (n - 7)
    r8 = _read(tmp_path / "s8", r12=_scores(eight), r21=_scores(54.0))
    assert r8["k10_confirmed"] and r8["tasks"]["r12"]["paired"]["s1x_res - s1x"]["ahead"] == 8
    r7 = _read(tmp_path / "s7", r12=_scores(54.0), r21=_scores(seven))
    assert not r7["k10_confirmed"] and r7["tasks"]["r21"]["paired"]["s1x_res - s1x"]["ahead"] == 7


def test_the_reading_is_void_on_missing_dirty_or_mixed_results(tmp_path):
    missing = _read(tmp_path / "m", r12=_scores(54.0))
    assert missing["void"] and any("missing task r21" in w for w in missing["why"])
    d = tmp_path / "dirty"
    d.mkdir()
    for name in k10.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k10.TASKS[name].key, _scores(54.0), clean=name != "r21")))
    assert k10.reading(d)["void"]
    d = tmp_path / "commits"
    d.mkdir()
    for i, name in enumerate(k10.TASKS):
        (d / f"{name}.json").write_text(json.dumps(_result(k10.TASKS[name].key, _scores(54.0), commit=f"c{i}")))
    assert k10.reading(d)["void"]
    d = tmp_path / "k9key"
    d.mkdir()
    for name in k10.TASKS:
        (d / f"{name}.json").write_text(json.dumps(_result(k9.TASKS["r12"].key, _scores(54.0))))
    assert k10.reading(d)["void"]


def test_the_reading_writes_no_per_gene_value_or_cell_number(tmp_path):
    res = _read(tmp_path / "w", r12=_scores(54.0), r21=_scores(53.0))
    blob = json.dumps(res) + k10.text(res)
    assert not any(g in re.findall(r"[A-Za-z0-9_]+", blob) for g in GENES)
    assert "n_cells" not in blob and "phi_src" not in json.dumps(res)


def test_run_and_read_refuse_to_overwrite(tmp_path):
    (tmp_path / "r12.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        k10.main(["run", "--task", "r12", "--out", str(tmp_path), "--data", str(tmp_path / "nowhere")])
    runs = tmp_path / "runs"
    _read(runs, r12=_scores(54.0), r21=_scores(53.0))
    out = tmp_path / "out"
    k10.main(["read", "--runs", str(runs), "--out", str(out)])
    assert json.loads((out / "reading.json").read_text())["k10_confirmed"]
    with pytest.raises(SystemExit, match="exists"):
        k10.main(["read", "--runs", str(runs), "--out", str(out)])
    (out / "reading.json").unlink()
    with pytest.raises(SystemExit, match="exists"):  # reading.txt alone also refuses
        k10.main(["read", "--runs", str(runs), "--out", str(out)])



def test_each_task_gets_a_board_whatever_its_truths_size():
    """r12 is scored on the published validation board, whose counts its files must reproduce; the mirror's truth is
    the other replicate, of another size, so it must take the twin form (a synthetic board with its own counts)."""
    from vec.boards import Board

    n_truth, n_ref = 4000, 3000
    pub = Board(key="T3:gata4", task="T3", genes=["g1", "g2"], needs_coords=True, split="validation",
                target_ko="gata4", **harness.pseudo_counts(n_truth, n_ref))
    boards = {"T3:gata4": pub}
    b12 = harness.task_board(k10.TASKS["r12"], boards, n_truth, n_ref)
    assert b12.key == "T3:gata4" and b12.anchors == {}
    with pytest.raises(ValueError, match="published counts"):
        harness.task_board(k10.TASKS["r12"], boards, n_truth + 600, n_ref)
    b21 = harness.task_board(k10.TASKS["r21"], boards, n_truth + 600, n_ref)  # the other replicate's size
    assert b21.key == "pseudo:T3k10/r21" and b21.split == "pseudo" and b21.genes == pub.genes
    assert b21.truth_cells == harness.pseudo_counts(n_truth + 600, n_ref)["truth_cells"] != pub.truth_cells
