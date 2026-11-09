"""The T1 late-proxy tournament (experiments/harness/TRIES.md T13-T18, N20), on synthetic data only: the bracket
argument of feng_dirs (default today's), the candidates' legality and hand-over, the controls and variants, the
within-type reads, the tournament's decision rule and its refusals."""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import GENES
from scipy import sparse

from vec import harness, t1_scope
from vec.boards import Board
from vec.evaluate import Stage
from vec.harness import CANDIDATES, PseudoTask, load_tasks
from vec.prep import log1p_cp10k

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments" / "t1_m1"))
sys.path.insert(0, str(ROOT / "experiments" / "t1_hygiene"))
import feng_dirs  # noqa: E402
import late_proxy  # noqa: E402

LP = ("lp_proxy", "lp_cflip", "lp_cperm", "lp_hyg", "lp_hyg_n", "lp_eb")
STAGES = ("E9.5", "E14.5", "E18.5")
LABELS = ("ventricular_cm", "atrial_cm", "endocardial_ec", "vascular_ec", "epicardial", "blood")


# ---------------------------------------------------------------------------------------------------------------------
# Synthetic Feng statistics


def _stats(path: Path, genes: list[str], *, stages=STAGES, seed: int = 3) -> Path:
    """cd1_stats.npz and its MANIFEST.json, in the clean reduction's layout: per stage|zone|label|half, a noisy mean
    CP10k profile with a stage trend (so the directions are not all zero)."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(1.0, 5.0, len(genes))
    keys, n, cp = [], [], []
    for label in LABELS:
        for i, stage in enumerate(stages):
            trend = np.exp((-1.0) ** (np.arange(len(genes)) % 2) * 0.25 * i)  # odd and even genes move opposite ways
            for half in "AB":
                k = int(rng.integers(30, 60))
                prof = base * trend * rng.lognormal(0.0, 0.2, len(genes))
                keys.append(f"{stage}|LV|{label}|{half}")
                n.append(k)
                cp.append(k * prof)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = np.array(n, dtype=np.float64)
    np.savez(path, genes=np.array(genes), keys=np.array(keys), n=n, S_cp=np.array(cp),
             S_log=np.zeros((len(keys), len(genes))), S_nz=np.zeros((len(keys), len(genes))),
             depth_sum=n * 1000.0)
    sha = __import__("hashlib").sha256(path.read_bytes()).hexdigest()
    (path.parent / "MANIFEST.json").write_text(json.dumps({"outputs": {path.name: {"sha256": sha}}}))
    return path


@pytest.fixture
def stats(tmp_path, monkeypatch) -> Path:
    monkeypatch.setattr(feng_dirs, "REPO", tmp_path)  # feng_directions records the file relative to the repository
    return _stats(tmp_path / "data" / "external" / "feng2022_gse193346_geo" / "cd1_stats.npz", list(GENES))


def _manual(path: Path, a_stage: str, b_stage: str, labels: tuple[str, ...], eps: float = 0.5) -> np.ndarray:
    z = np.load(path, allow_pickle=True)
    keys = [str(k) for k in z["keys"]]

    def mean(stage):
        m = np.array([k.split("|")[0] == stage and k.split("|")[2] in labels for k in keys])
        return z["S_cp"][m].sum(0) / z["n"][m].sum()

    return np.log((mean(b_stage) + eps) / (mean(a_stage) + eps))


# ---------------------------------------------------------------------------------------------------------------------
# feng_dirs: the bracket argument, default today's bracket


def test_the_default_bracket_is_todays_and_the_late_bracket_is_e14_5_to_e18_5(stats):
    out, rec = feng_dirs.feng_directions(list(GENES), stats_path=stats)
    assert rec["bracket"] == ["E9.5", "E14.5"]
    np.testing.assert_allclose(out["CM"], _manual(stats, "E9.5", "E14.5", ("ventricular_cm", "atrial_cm")))
    # the late bracket: the same formula across E14.5 -> E18.5, and every group of the default bracket
    late, rec_late = feng_dirs.feng_directions(list(GENES), stats_path=stats, bracket=("E14.5", "E18.5"))
    assert rec_late["bracket"] == ["E14.5", "E18.5"] and set(late) == set(out)
    np.testing.assert_allclose(late["CM"], _manual(stats, "E14.5", "E18.5", ("ventricular_cm", "atrial_cm")))
    assert feng_dirs.bracket_fraction(10.5, 12.5) == pytest.approx(0.4)
    assert feng_dirs.bracket_fraction(9.5, 10.5) == pytest.approx(0.2)
    assert feng_dirs.bracket_fraction(9.5, 10.5, (14.5, 18.5)) == pytest.approx(0.25)


def test_directions_from_stats_is_what_feng_directions_computes_with_no_file_reads(stats):
    st = feng_dirs.dl.load_stats(stats, fields=feng_dirs.dl.STAT_FIELDS)
    cfg = feng_dirs.load_config()
    a, cells, idx = feng_dirs.directions_from_stats(st, list(GENES), cfg, ("E14.5", "E18.5"))
    b, rec = feng_dirs.feng_directions(list(GENES), cfg, stats, bracket=("E14.5", "E18.5"))
    assert set(a) == set(b) and cells == rec["groups"] and (idx >= 0).sum() == len(GENES)
    for g in a:
        np.testing.assert_array_equal(a[g], b[g])


def test_bracket_stats_keeps_only_the_two_ends_and_refuses_a_missing_one(stats):
    st = feng_dirs.dl.load_stats(stats, fields=feng_dirs.dl.STAT_FIELDS)
    cut = feng_dirs.bracket_stats(st, ("E14.5", "E18.5"))
    assert {k.split("|")[0] for k in cut["keys"]} == {"E14.5", "E18.5"} and len(cut["genes"]) == len(GENES)
    assert cut["S_cp"].shape == (len(cut["keys"]), len(GENES)) and len(cut["n"]) == len(cut["keys"])
    with pytest.raises(ValueError, match="no stage"):
        feng_dirs.bracket_stats(st, ("E14.5", "E19.5"))


def test_eb_shrink_scales_each_decile_by_its_signal_share_and_keeps_nan():
    d = np.array([2.0, -2.0, 0.1, -0.1, np.nan])
    d_a, d_b = d + 0.0, d + 0.0
    d_a[[2, 3]] = [0.3, 0.1]  # the second decile's halves disagree: (0.3 - (-0.1)) / 2 = 0.2 and 0.1
    d_b[[2, 3]] = [-0.1, -0.3]
    basis = np.array([9.0, 9.0, 1.0, 1.0, np.nan])
    out, factors = feng_dirs.eb_shrink(d, d_a, d_b, basis, bins=2)
    assert np.isnan(out[4])  # no value stays no value
    # the low-expression decile (0.1 and -0.1): noise = mean(((d_a - d_b) / 2)^2) = mean(0.04, 0.04) = 0.04; signal =
    # max(mean(d^2) - noise, 0) = max(0.01 - 0.04, 0) = 0, so its genes shrink to 0
    assert factors[0] == 0.0 and out[2] == 0.0 and out[3] == 0.0
    # the high-expression decile: no disagreement, noise 0, factor 1, unchanged
    assert factors[1] == 1.0 and out[0] == 2.0 and out[1] == -2.0


def test_direction_halves_use_each_half_alone_and_the_mean_of_the_two_ends(stats):
    st = feng_dirs.dl.load_stats(stats, fields=feng_dirs.dl.STAT_FIELDS)
    cfg = feng_dirs.load_config()
    h = feng_dirs.direction_halves(st, list(GENES), cfg, ("E14.5", "E18.5"))["CM"]
    z = np.load(stats, allow_pickle=True)
    keys = [str(k) for k in z["keys"]]

    def mean_cp(stage, halves=("A", "B")):
        m = np.array([k.split("|")[0] == stage and k.split("|")[2] in ("ventricular_cm", "atrial_cm")
                      and k.split("|")[3] in halves for k in keys])
        return z["S_cp"][m].sum(0) / z["n"][m].sum()

    np.testing.assert_allclose(h["d_a"], np.log((mean_cp("E18.5", ("A",)) + 0.5) / (mean_cp("E14.5", ("A",)) + 0.5)))
    np.testing.assert_allclose(h["d_b"], np.log((mean_cp("E18.5", ("B",)) + 0.5) / (mean_cp("E14.5", ("B",)) + 0.5)))
    np.testing.assert_allclose(h["mean_cp"], (mean_cp("E14.5") + mean_cp("E18.5")) / 2.0)


# ---------------------------------------------------------------------------------------------------------------------
# The candidates: legal on the three T1 tasks, what they are handed


def _within_task(key: str) -> PseudoTask:
    return {t.key: t for t in load_tasks(late_proxy.CONFIG)["T1:late_proxy"]}[key]


def _t1_tasks() -> dict[str, PseudoTask]:
    return {t.key: t for t in load_tasks()["T1:test"]}


def test_the_late_proxy_candidates_are_legal_on_the_three_t1_tasks_and_declare_their_bracket():
    tasks = _t1_tasks()
    for name in LP:
        c = CANDIDATES[name]
        assert c.externals == (("X002", 14.5, 18.5),) and c.params["scope"] == "chamber"
        for t in tasks.values():
            assert c.applies(t), (name, t.key)
            harness.check_externals(c, t)  # legal: the bracket lies outside the pseudo-task's version of the window


def test_the_candidates_take_the_held_flat_list_only_where_its_stages_are_inputs(monkeypatch, tmp_path):
    tasks = _t1_tasks()
    c = CANDIDATES["lp_proxy"]
    assert c.derived == ("t1_late_proxy",) and c.derived_optional == ("t1_held_flat",)
    assert "t1_held_flat" in c.derived_for(tasks["T1/E9.5->E10.5"])
    assert "t1_held_flat" not in c.derived_for(tasks["T1/E8.5->E9.5"])  # its pool rule was computed on E8.5 and E9.5
    assert c.applies(tasks["T1/E8.5->E9.5"])  # an optional input does not stop a candidate where it cannot go


def test_the_hand_over_holds_only_the_bracket_stages_and_the_module_lists(monkeypatch, tmp_path):
    p = _stats(tmp_path / "x" / "cd1_stats.npz", list(GENES))
    bundle = harness._late_proxy_load(p)
    assert {str(k).split("|")[0] for k in bundle["stats"]["keys"]} == {"E14.5", "E18.5"}  # no other Feng stage
    assert {"Xist", "Ddx3y"} <= set(bundle["named"]) and set(bundle["modules"]) >= {"D", "A1", "N1", "N2"}


def _t1(n, m, labels, rng):
    """A released-style T1 stage: counts layer and log1p CP10k, the released labels the chamber config knows."""
    types = rng.choice(len(labels), size=n)
    C = rng.poisson(m[types] * 30).astype(np.int32)
    C[np.arange(n), rng.integers(0, len(GENES), n)] += 1
    a = ad.AnnData(X=log1p_cp10k(C), obs=pd.DataFrame({"celltype": np.asarray(labels)[types]},
                                                       index=[f"c{i}" for i in range(n)]),
                   var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    return a


@pytest.fixture
def t1_ctx(tmp_path, monkeypatch):
    """A BuildContext on the E9.5 -> E10.5 geometry, synthetic stages with released labels and Feng's statistics
    handed over as the harness would."""
    stats_path = _stats(tmp_path / "feng" / "cd1_stats.npz", list(GENES))
    monkeypatch.setitem(harness.DERIVED, "t1_late_proxy", replace(harness.DERIVED["t1_late_proxy"], path=stats_path))
    rng = np.random.default_rng(11)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    ref = _t1(1500, base, ("V-CM", "Endocardium", "Blood"), rng)
    start = _t1(1500, base * rng.lognormal(0.0, 0.3, size=base.shape), ("V-CM", "Endothelium", "Blood"), rng)
    task = PseudoTask(key="t1lp", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      twin_board="T1:val")
    board = Board(key="pseudo:t1lp", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600,
                  truth_cells=200, ref_cells=300)
    cboard = harness._candidate_board(board)

    def ctx(candidate):
        return harness.BuildContext(task=harness.TaskView.of(task), board=cboard, inputs={"ref": ref, "start": start},
                                    seed=1_000_003, derived=harness.load_derived(candidate, task))

    return ctx, task, stats_path


def _spy(monkeypatch):
    calls = []
    real = t1_scope.build_submission

    def spy(*a, **kw):
        calls.append(kw)
        return real(*a, **kw)

    monkeypatch.setattr(t1_scope, "build_submission", spy)
    return calls


def test_lp_proxy_moves_m1s_groups_along_the_proxy_at_its_time_fraction_and_leaves_the_rest(t1_ctx, monkeypatch):
    ctx_of, task, stats_path = t1_ctx
    calls = _spy(monkeypatch)
    pred = CANDIDATES["lp_proxy"].build(ctx_of(CANDIDATES["lp_proxy"]))
    base = CANDIDATES["b0_scope"].build(ctx_of(CANDIDATES["b0_scope"]))
    kw = calls[0]
    assert set(kw["directions"]) <= set(harness.LP_GROUPS) and "CM" in kw["directions"]  # Blood never moves
    assert kw["direction_fraction"] == pytest.approx(0.25)
    assert kw["direction_groups"]["V-CM"] == "CM" and kw["direction_groups"]["Endocardium"] == "Endocardium"
    assert kw.get("start") is None and kw.get("damping", 0.0) == 0.0  # no momentum: the cells without a direction stay
    np.testing.assert_allclose(kw["directions"]["CM"], _manual(stats_path, "E14.5", "E18.5",
                                                               ("ventricular_cm", "atrial_cm")))
    assert pred.n == base.n
    assert not np.array_equal(pred.X, base.X)
    # the pre-registered held-flat list stands on this task (its stages are inputs), every name of it held
    want = set(t1_scope.held_flat_genes())
    assert set(kw["held_flat"]) >= want


def test_the_held_flat_list_is_named_genes_only_where_its_pool_rule_would_leak(t1_ctx, monkeypatch):
    ctx_of, task, stats_path = t1_ctx
    t = replace(task, key="t1lpe", target="T1/E9.5_RNA.h5ad", reference="T1/E8.5_RNA.h5ad",
                inputs={"ref": "T1/E8.5_RNA.h5ad"})
    c = CANDIDATES["lp_proxy"]
    calls = _spy(monkeypatch)
    rng = np.random.default_rng(13)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    ref = _t1(1200, base, ("LV-CM", "Endothelium", "Blood"), rng)
    ctx = harness.BuildContext(task=harness.TaskView.of(t), board=ctx_of(c).board, inputs={"ref": ref}, seed=1_000_003,
                               derived=harness.load_derived(c, t))
    c.build(ctx)
    held = set(calls[0]["held_flat"])
    assert {"Xist", "Ddx3y", "Cwc22"} <= held  # the named genes
    assert not {"Lin28a", "Abcb7", "Gm7324"} & held  # pool-rule genes are not handed over on this task
    assert calls[0]["direction_fraction"] == pytest.approx(0.25)  # E8.5 -> E9.5, one day of the four-day bracket


def test_the_controls_flip_and_permute_within_deciles_and_nothing_else(t1_ctx, monkeypatch):
    ctx_of, *_ = t1_ctx
    calls = _spy(monkeypatch)
    for name in ("lp_proxy", "lp_cflip", "lp_cperm"):
        CANDIDATES[name].build(ctx_of(CANDIDATES[name]))
    proxy, flip, perm = (calls[i]["directions"] for i in range(3))
    for g in proxy:
        np.testing.assert_allclose(flip[g], -proxy[g])  # Cflip: the sign and nothing else
        assert np.array_equal(np.sort(perm[g][np.isfinite(perm[g])]), np.sort(proxy[g][np.isfinite(proxy[g])]))
        assert not np.array_equal(perm[g], proxy[g])  # same values on other genes (a fixed permutation)
    again = _spy(monkeypatch)
    CANDIDATES["lp_cperm"].build(ctx_of(CANDIDATES["lp_cperm"]))
    for g in perm:
        np.testing.assert_array_equal(again[0]["directions"][g], perm[g])  # deterministic: one permutation, its seed


def test_the_hygiene_variants_hold_their_modules_flat_and_the_shrinkage_is_eb(t1_ctx, monkeypatch):
    ctx_of, task, stats_path = t1_ctx
    calls = _spy(monkeypatch)
    for name in ("lp_proxy", "lp_hyg", "lp_hyg_n", "lp_eb"):
        CANDIDATES[name].build(ctx_of(CANDIDATES[name]))
    proxy, hyg, hyg_n, eb = calls
    bundle = harness._late_proxy_load(stats_path)
    modules = bundle["modules"]
    assert set(hyg["held_flat"]) == set(proxy["held_flat"]) | set(modules["D"]) | set(modules["A1"])
    assert set(hyg_n["held_flat"]) == set(hyg["held_flat"]) | set(modules["N1"]) | set(modules["N2"])
    # the shrinkage: lp_proxy's direction, each group by feng_dirs.eb_shrink (its deciles of the group's mean CP10k)
    st = bundle["stats"]
    halves = feng_dirs.direction_halves(st, list(GENES), bundle["cfg"], harness.LP_BRACKET)
    for g, d in proxy["directions"].items():
        shrunk, _ = feng_dirs.eb_shrink(d, halves[g]["d_a"], halves[g]["d_b"], halves[g]["mean_cp"], harness.LP_BINS)
        np.testing.assert_allclose(eb["directions"][g], shrunk)


def test_a_late_proxy_candidate_runs_inside_the_read_guard_and_records_its_hand_over(tmp_path, monkeypatch):
    stats_path = _stats(tmp_path / "feng" / "cd1_stats.npz", list(GENES))
    monkeypatch.setitem(harness.DERIVED, "t1_late_proxy", replace(harness.DERIVED["t1_late_proxy"], path=stats_path))
    rng = np.random.default_rng(5)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    later = base.copy()
    later[:, :8] *= 4.0
    later[:, 8:16] *= 0.25
    ref_a = _t1(3000, base, ("V-CM", "Endocardium", "Blood"), rng)
    start_a = _t1(2500, base, ("LV-CM", "Endothelium", "Blood"), rng)
    tgt_a = _t1(4000, later, ("V-CM", "Endocardium", "Blood"), rng)
    task = PseudoTask(key="t1lpi", test_board="T1:test", task="T1", kind="extrap", target="T1/E10.5_RNA.h5ad",
                      reference="T1/E9.5_RNA.h5ad", inputs={"ref": "T1/E9.5_RNA.h5ad", "start": "T1/E8.5_RNA.h5ad"},
                      twin_board="T1:val")
    board = Board(key="pseudo:t1lpi", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=600,
                  truth_cells=200, ref_cells=300)
    inputs = {"T1/E9.5_RNA.h5ad": ref_a, "T1/E8.5_RNA.h5ad": start_a}
    target = Stage(X=tgt_a.X.toarray())
    ref = Stage(X=ref_a.X.toarray())

    def read(path, genes):
        st = target if str(path).endswith(task.target) else ref
        return ad.AnnData(X=np.asarray(st.X), var=pd.DataFrame(index=GENES))

    monkeypatch.setattr(harness, "read_stage", read)
    monkeypatch.setattr(harness, "load_stage", lambda path, **kw: inputs[str(path).split(str(tmp_path) + "/")[1]])
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    monkeypatch.setattr(harness, "task_board", lambda t, boards, n_t, n_r: board)
    # the record's relative path needs a file under the repository, which the synthetic statistics are not
    monkeypatch.setattr(harness, "derived_record", lambda c, t: {n: {"file": n} for n in c.derived_for(t)})
    res = harness.run_task(task, [CANDIDATES["b0_scope"], CANDIDATES["lp_proxy"]], seeds=(0,), data_dir=tmp_path,
                           boards={"T1:val": replace(board, key="T1:val")}, log=None)
    assert res["candidates"] == ["floor", "b0_scope", "lp_proxy"] and res["rows"][0]["scores"]["floor"] == 50.0
    assert not any("results" in p or p.startswith("data") for p in res["rows"][0]["reads"]["lp_proxy"])
    assert set(res["params"]["lp_proxy"]["derived"]) == {"t1_late_proxy", "t1_held_flat"}
    assert [tuple(e) for e in res["params"]["lp_proxy"]["externals"]] == [("X002", 14.5, 18.5)]


# ---------------------------------------------------------------------------------------------------------------------
# The within-type reads


def test_the_within_type_config_is_three_reads_that_pass_the_rules():
    tasks = load_tasks(late_proxy.CONFIG)["T1:late_proxy"]
    assert {t.key for t in tasks} == set(late_proxy.WITHIN)
    assert all(harness.check_task(t) is t and t.target_filter["within"] in {"cm", "endocardium"} and t.counted
               for t in tasks)
    with pytest.raises(ValueError, match="category"):
        harness.check_task(replace(tasks[0], target_filter={"within": "not_a_category"}))


def test_a_within_type_read_cuts_target_reference_and_inputs_to_the_category(tmp_path, monkeypatch):
    """The cut (vec.harness.run_task): the target through filter_target, the reference and every input stage alike;
    the board counts follow the cut files, and no surrounding cell reaches a candidate."""
    t = _within_task("T1/E9.5->E10.5:cm")
    rng = np.random.default_rng(3)
    base = rng.gamma(2.0, 2.0, size=(2, len(GENES)))
    target, ref, start = (_t1(n, base, ("V-CM", "Endocardium"), rng) for n in (800, 600, 500))
    files = {t.target: target, t.reference: ref, t.inputs["start"]: start}

    def read(path, genes):
        return files[str(path).split(str(tmp_path) + "/")[1]].copy()

    monkeypatch.setattr(harness, "read_stage", read)
    monkeypatch.setattr(harness, "load_stage", lambda path, **kw: files[str(path).split(str(tmp_path) + "/")[1]].copy())
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    board = Board(key=f"pseudo:{t.key}", task="T1", genes=list(GENES), needs_coords=False, min_cells=10, max_cells=600,
                  truth_cells=20, ref_cells=30)
    seen = {}

    def task_board(t_, boards, n_target, n_ref):
        seen["counts"] = (n_target, n_ref)
        return board

    monkeypatch.setattr(harness, "task_board", task_board)

    def score_seeds(t_, board_, tgt, ref_, inputs, candidates, **kw):
        seen["target"], seen["ref"], seen["inputs"] = tgt, ref_, inputs
        return []

    monkeypatch.setattr(harness, "score_seeds", score_seeds)
    harness.run_task(t, [CANDIDATES["b0_scope"]], seeds=(0,), data_dir=tmp_path, boards={"T1:val": board}, log=None)
    want = {name: (a.obs["celltype"].to_numpy() == "V-CM").sum() for name, a in (("target", target), ("ref", ref))}
    assert seen["counts"] == (want["target"], want["ref"]) and seen["target"].n == want["target"]
    assert seen["ref"].n == want["ref"]
    assert seen["inputs"]["ref"].n_obs == want["ref"] and set(seen["inputs"]["start"].obs["celltype"]) == {"V-CM"}
    assert set(seen["inputs"]["ref"].obs["celltype"]) == {"V-CM"}


def test_a_within_type_read_scores_the_late_proxy_candidates_end_to_end_inside_the_guard(tmp_path, monkeypatch):
    t = _within_task("T1/E9.5->E10.5:cm")
    stats_path = _stats(tmp_path / "feng" / "cd1_stats.npz", list(GENES))
    monkeypatch.setitem(harness.DERIVED, "t1_late_proxy", replace(harness.DERIVED["t1_late_proxy"], path=stats_path))
    monkeypatch.setattr(harness, "derived_record", lambda c, t_: {n: {"file": n} for n in c.derived_for(t_)})
    rng = np.random.default_rng(21)
    base = rng.gamma(2.0, 2.0, size=(3, len(GENES)))
    later = base.copy()
    later[:, :8] *= 4.0
    later[:, 8:16] *= 0.25
    files = {t.target: _t1(6000, later, ("V-CM", "Endocardium", "Blood"), rng),
             t.reference: _t1(4000, base, ("V-CM", "Endocardium", "Blood"), rng),
             t.inputs["start"]: _t1(3000, base, ("LV-CM", "Endothelium", "Blood"), rng)}
    monkeypatch.setattr(harness, "read_stage", lambda p, g: files[str(p).split(str(tmp_path) + "/")[1]].copy())
    monkeypatch.setattr(harness, "load_stage", lambda p, **kw: files[str(p).split(str(tmp_path) + "/")[1]])
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    board = Board(key=f"pseudo:{t.key}", task="T1", genes=list(GENES), needs_coords=False, min_cells=100, max_cells=800,
                  truth_cells=100, ref_cells=130)
    seen = []
    monkeypatch.setattr(harness, "task_board", lambda t_, b, n_t, n_r: seen.append((n_t, n_r)) or board)
    res = harness.run_task(t, [CANDIDATES["b0_scope"], CANDIDATES["lp_proxy"], CANDIDATES["lp_cflip"]], seeds=(0,),
                           data_dir=tmp_path, boards={"T1:val": replace(board, key="T1:val")}, log=None)
    n_cm = {k: int((a.obs["celltype"] == "V-CM").sum()) for k, a in files.items()}
    assert seen == [(n_cm[t.target], n_cm[t.reference])]  # the board counts are the cut cells'
    assert res["candidates"] == ["floor", "b0_scope", "lp_proxy", "lp_cflip"]
    assert res["rows"][0]["scores"]["floor"] == 50.0
    assert not any("results" in p or p.startswith("data") for p in res["rows"][0]["reads"]["lp_proxy"])
    assert res["target_filter_record"]["within"] == "cm" and res["target_filter_record"]["kept"] == n_cm[t.target]
    assert res["rows"][0]["scores"]["lp_proxy"] != res["rows"][0]["scores"]["b0_scope"]


def test_a_whole_population_task_is_unchanged_by_the_within_option(tmp_path, monkeypatch):
    """A task with no `within` loads exactly what it did: the reference and inputs uncut."""
    t = _t1_tasks()["T1/E9.5->E10.5:chamber"]
    assert t.target_filter.get("within") is None
    rng = np.random.default_rng(4)
    base = rng.gamma(2.0, 2.0, size=(2, len(GENES)))
    target, ref, start = (_t1(n, base, ("V-CM", "Foregut"), rng) for n in (800, 600, 500))
    files = {t.target: target, t.reference: ref, t.inputs["start"]: start}
    monkeypatch.setattr(harness, "read_stage", lambda p, g: files[str(p).split(str(tmp_path) + "/")[1]].copy())
    monkeypatch.setattr(harness, "load_stage", lambda p, **kw: files[str(p).split(str(tmp_path) + "/")[1]])
    monkeypatch.setattr(harness, "_sha256", lambda path: "0" * 64)
    board = Board(key=f"pseudo:{t.key}", task="T1", genes=list(GENES), needs_coords=False, min_cells=10, max_cells=600,
                  truth_cells=20, ref_cells=30)
    seen = {}
    monkeypatch.setattr(harness, "task_board", lambda t_, b, n_t, n_r: seen.setdefault("counts", (n_t, n_r)) and board)

    def score_seeds(t_, board_, tgt, ref_, inputs, candidates, **kw):
        seen["ref"], seen["inputs"] = ref_, inputs
        return []

    monkeypatch.setattr(harness, "score_seeds", score_seeds)
    harness.run_task(t, [CANDIDATES["b0_scope"]], seeds=(0,), data_dir=tmp_path, boards={"T1:val": board}, log=None)
    assert seen["ref"].n == ref.n_obs and seen["inputs"]["ref"].n_obs == ref.n_obs
    assert seen["inputs"]["start"].n_obs == start.n_obs and seen["counts"][1] == ref.n_obs


# ---------------------------------------------------------------------------------------------------------------------
# The decision rule


GAINS_OK = [0.6, 0.5, 0.7, 0.4, 0.6, 0.5, 0.45, 0.55, 0.65, 0.5]


def _res(read: str, scores: dict[str, list[float]], groups: bool = True, commit: str = "c" * 40) -> dict:
    rows = [{"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
             "scores": {c: 50.0 + v[s] for c, v in scores.items()},
             "groups": {c: {"de_recovery": 0.5 + v[s] / 100, "direction": 0.5, "distribution": 0.5, "covariance": 0.5}
                        for c, v in scores.items()} if groups else {}}
            for s in harness.SEEDS]
    return {"task": {"key": read, "test_board": "T1:late_proxy"}, "candidates": list(scores), "rows": rows,
            "protocol": "bake", "build": "main", "anchors": "computed", "frac": 0.1, "seeds": list(harness.SEEDS),
            "provenance": {"commit": commit, "tree_clean": True, "config_sha256": "c" * 64,
                           "files_sha256": {"T1/E9.5_RNA.h5ad": "f"}},
            "late_proxy": {"stage": "gate", "inputs": {"feng_stats_sha256": "a" * 64}}}


def _gate_results(over: dict | None = None, reads=late_proxy.READS, **kw) -> dict:
    """Scores per read where lp_proxy beats b0_scope (0) and the controls by about 0.55 in every seed; `over` replaces
    a candidate's gains on a read, e.g. {("T1/E8.5->E9.5:cm", "lp_proxy"): [...]}."""
    over = over or {}
    out = {}
    for read in reads:
        scores = {"b0_scope": [0.0] * 10, "lp_proxy": GAINS_OK, "lp_cflip": [-0.2] * 10, "lp_cperm": [-0.1] * 10}
        for (r, c), v in over.items():
            if r == read:
                scores[c] = v
        out[read] = _res(read, scores, **kw)
    return out


def test_the_gate_passes_only_with_a_held_out_gain_a_within_type_e10_gain_and_both_controls():
    v = late_proxy.gate_verdict(_gate_results())
    assert v["passed"] and all(v["conditions"].values())
    assert set(v["comparisons"]) == {"lp_proxy vs b0_scope", "lp_proxy vs lp_cflip", "lp_proxy vs lp_cperm",
                                     "lp_cperm vs b0_scope", "lp_cflip vs b0_scope", "lp_proxy vs b0 (reported)",
                                     "lp_proxy vs the floor (reported)"}


@pytest.mark.parametrize("over, condition", [
    ({("T1/E8.5->E9.5:cm", "lp_proxy"): [0.1] * 10}, "held_out"),  # the held-out read must pass
    ({(r, "lp_proxy"): [0.1] * 10 for r in late_proxy.WITHIN_E10}, "within_type_e10"),  # one E10.5 read within type
    ({("T1/E8.5->E9.5:cm", "lp_cflip"): [0.5] * 10}, "beats_cflip"),  # and it must beat the flipped control there
    ({("T1/E9.5->E10.5:cm", "lp_cperm"): [0.5] * 10}, "beats_cperm"),  # and the permuted control
    ({("T1/E8.5->E9.5:cm", "lp_cperm"): [0.6] * 10}, "beats_cperm"),  # (a permuted control that gains as well fails)
    ({("T1/E9.5->E10.5:endocardium", "lp_proxy"): [-0.5] * 10}, "no_veto"),  # a within-type read's mean loss above 0.3
])
def test_the_gate_fails_for_each_broken_condition(over, condition):
    v = late_proxy.gate_verdict(_gate_results(over))
    assert not v["passed"] and not v["conditions"][condition]


def test_a_whole_population_loss_is_reported_and_decides_nothing():
    v = late_proxy.gate_verdict(_gate_results({("T1/E9.5->E10.5:chamber", "lp_proxy"): [-0.5] * 10}))
    assert v["passed"] and v["comparisons"]["lp_proxy vs b0_scope"]["T1/E9.5->E10.5:chamber"]["vetoes"]


@pytest.mark.parametrize("control", ["lp_cperm", "lp_cflip"])
def test_the_gate_fails_when_a_control_itself_passes_the_held_out_read(control):
    v = late_proxy.gate_verdict(_gate_results({("T1/E8.5->E9.5:cm", control): [0.5] * 10,
                                               ("T1/E8.5->E9.5:cm", "lp_proxy"): [1.5] * 10}))
    # the control passes on H (0.5 > 0.3) and lp_proxy beats it by more than 0.3 there; the gate still fails
    assert v["conditions"]["beats_cperm"] and v["conditions"]["beats_cflip"] and not v["passed"]
    assert not v["conditions"]["controls_do_not_pass_held_out"]


def test_a_missing_read_fails_its_condition_and_a_complete_set_decides():
    results = _gate_results()
    del results["T1/E8.5->E9.5:cm"]
    v = late_proxy.gate_verdict(results)
    assert not v["passed"] and not v["conditions"]["held_out"] and not v["conditions"]["no_veto"]


def _variant_results(variant_gains: dict[str, list[float]] | None = None) -> dict:
    """Variant-stage scores per read, each variant 0.6 ahead of its reference in every seed unless overridden: lp_hyg
    and lp_eb on lp_proxy, lp_hyg_n on lp_hyg."""
    variant_gains = variant_gains or {}
    out = {}
    for read in late_proxy.READS:
        scores = {"b0_scope": [0.0] * 10, "lp_proxy": list(GAINS_OK)}  # the gate's lp_proxy, as the rerun reproduces it
        scores["lp_hyg"] = [a + b for a, b in zip(scores["lp_proxy"], variant_gains.get("lp_hyg", [0.6] * 10),
                                                  strict=True)]
        scores["lp_eb"] = [a + b for a, b in zip(scores["lp_proxy"], variant_gains.get("lp_eb", [0.6] * 10),
                                                 strict=True)]
        scores["lp_hyg_n"] = [a + b for a, b in zip(scores["lp_hyg"], variant_gains.get("lp_hyg_n", [0.6] * 10),
                                                    strict=True)]
        out[read] = _res(read, scores)
    return out


def test_a_variant_needs_nine_of_ten_seeds_ahead_on_the_held_out_read_and_a_within_type_e10_read():
    res = _variant_results()
    for variant, reference in late_proxy.VARIANTS.items():
        v = late_proxy.variant_verdict(res, variant, reference)
        assert v["passed"], variant
    # 8 of 10 seeds ahead passes the adoption rule's test but not the stricter bar for the variants
    gains = [0.6] * 8 + [-0.05, -0.05]
    res = _variant_results({"lp_hyg": gains})
    v = late_proxy.variant_verdict(res, "lp_hyg", "lp_proxy")
    assert not v["passed"] and not v["conditions"]["held_out"]
    assert v["comparison"]["T1/E8.5->E9.5:cm"]["passes"] and not v["comparison"]["T1/E8.5->E9.5:cm"]["passes_strict"]


def test_a_variant_is_vetoed_by_a_within_type_loss_and_only_reported_for_a_whole_population_loss():
    res = _variant_results()
    for row in res["T1/E9.5->E10.5:endocardium"]["rows"]:
        row["scores"]["lp_eb"] = row["scores"]["lp_proxy"] - 0.5
    v = late_proxy.variant_verdict(res, "lp_eb", "lp_proxy")
    assert not v["passed"] and not v["conditions"]["no_veto"]
    res = _variant_results()
    for row in res["T1/E9.5->E10.5"]["rows"]:
        row["scores"]["lp_eb"] = row["scores"]["lp_proxy"] - 0.5
    v = late_proxy.variant_verdict(res, "lp_eb", "lp_proxy")
    assert v["passed"] and v["comparison"]["T1/E9.5->E10.5"]["vetoes"]  # reported in the table, deciding nothing


def test_group_margins_sum_to_the_score_margin_scale():
    res = _res("T1/E8.5->E9.5:cm", {"a": [1.0] * 10, "b": [0.0] * 10})
    m = late_proxy.group_margins(res, "a", "b")
    assert set(m) == {"de_recovery", "direction", "distribution", "covariance"}
    assert m["de_recovery"] == pytest.approx(100 * 0.25 * 0.01)


# ---------------------------------------------------------------------------------------------------------------------
# The commands: reads, refusals


def _write_results(folder: Path, results: dict) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for read, r in results.items():
        (folder / f"{late_proxy.safe(read)}.json").write_text(json.dumps(r))


def _ready(monkeypatch, tmp_path, proceed: bool = True) -> str:
    """The rows and N20's gate that `run` checks before anything else."""
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"N20", "T13", "T14", "T15", "T16", "T17", "T18"})
    n20 = tmp_path / "n20.json"
    n20.write_text(json.dumps({"proceed_to_tries": proceed}))
    return str(n20)


def test_run_needs_the_rows_and_n20s_yes(tmp_path, monkeypatch):
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"T13"})
    n20 = tmp_path / "n20.json"
    with pytest.raises(SystemExit, match="not in experiments/harness/TRIES.md"):
        late_proxy.main(["run", "--stage", "gate", "--out", str(tmp_path), "--n20", str(n20)])
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"N20", "T13", "T14", "T15"})
    with pytest.raises(SystemExit, match="N20 runs before the gate"):
        late_proxy.main(["run", "--stage", "gate", "--out", str(tmp_path), "--n20", str(n20)])
    n20.write_text(json.dumps({"proceed_to_tries": False}))
    with pytest.raises(SystemExit, match="said no"):
        late_proxy.main(["run", "--stage", "gate", "--out", str(tmp_path), "--n20", str(n20)])
    # the variants run only after a passed gate
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"T16", "T17", "T18"})
    with pytest.raises(SystemExit, match="only after the gate passed"):
        late_proxy.main(["run", "--stage", "variants", "--out", str(tmp_path), "--n20", str(n20)])
    (tmp_path / "decision_gate.json").write_text(json.dumps({"verdict": {"passed": False}}))
    with pytest.raises(SystemExit, match="only after the gate passed"):
        late_proxy.main(["run", "--stage", "variants", "--out", str(tmp_path), "--n20", str(n20)])


def test_registered_rows_reads_the_ledgers_table_rows(tmp_path, monkeypatch):
    ledger = tmp_path / "TRIES.md"
    ledger.write_text("# ledger\n\n| # | Board |\n|---|---|\n| T13 | T1 |\n| N20 (diagnostic, not a try) | T1 |\n"
                      "| R2 (ruling) | x |\n")
    monkeypatch.setattr(late_proxy, "TRIES", ledger)
    assert late_proxy.registered_rows() == {"T13", "N20", "R2"}


def test_run_writes_each_read_stamped_and_refuses_to_overwrite(tmp_path, monkeypatch):
    calls = []

    def run_task(t, cands, **kw):
        calls.append(t.key)
        return _res(t.key, {"b0_scope": [0.0] * 10, "lp_proxy": [0.5] * 10})

    monkeypatch.setattr(late_proxy.harness, "run_task", run_task)
    monkeypatch.setattr(late_proxy, "input_hashes", lambda: {"feng_stats_sha256": "a" * 64})
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: ("c" * 40, []))
    n20 = _ready(monkeypatch, tmp_path)
    late_proxy.main(["run", "--stage", "gate", "--out", str(tmp_path), "--task", "T1_E8.5_E9.5_cm", "--n20", n20])
    path = tmp_path / "gate" / "T1_E8.5_E9.5_cm.json"
    res = json.loads(path.read_text())
    assert calls == ["T1/E8.5->E9.5:cm"] and res["late_proxy"]["stamp"]["tree_clean"] is True
    with pytest.raises(SystemExit, match="exists"):
        late_proxy.main(["run", "--stage", "gate", "--out", str(tmp_path), "--task", "T1/E8.5->E9.5:cm", "--n20", n20])


def test_a_failed_read_is_recorded_and_passes_nothing(tmp_path, monkeypatch):
    def run_task(t, cands, **kw):
        raise ValueError("inverted anchors")

    monkeypatch.setattr(late_proxy.harness, "run_task", run_task)
    monkeypatch.setattr(late_proxy, "input_hashes", lambda: {"feng_stats_sha256": "a" * 64})
    n20 = _ready(monkeypatch, tmp_path)
    with pytest.raises(SystemExit, match="failed"):
        late_proxy.main(["run", "--stage", "gate", "--out", str(tmp_path), "--task", "T1/E8.5->E9.5:cm", "--n20", n20])
    assert (tmp_path / "gate" / "T1_E8.5_E9.5_cm.failed.txt").read_text().startswith("ValueError: inverted anchors")


def _decide(tmp_path, monkeypatch, stage, results):
    monkeypatch.setattr(late_proxy, "check_consistent", lambda res, tasks: None)
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: ("c" * 40, []))
    _write_results(tmp_path / stage, results)
    late_proxy.main(["decide", "--stage", stage, "--runs", str(tmp_path), "--out", str(tmp_path)])
    return json.loads((tmp_path / f"decision_{stage}.json").read_text())


def test_decide_writes_the_gate_and_refuses_variants_before_a_passed_gate(tmp_path, monkeypatch, capsys):
    out = _decide(tmp_path, monkeypatch, "gate", _gate_results())
    assert out["verdict"]["passed"] and out["constants"]["min_ahead"] == 8 and out["constants"]["variant_ahead"] == 9
    assert "PASSED" in capsys.readouterr().out
    with pytest.raises(SystemExit, match="never overwritten"):
        late_proxy.main(["decide", "--stage", "gate", "--runs", str(tmp_path), "--out", str(tmp_path)])
    # a variants decision without a gate decision is refused; after a failed gate too
    other = tmp_path / "other"
    other.mkdir()
    _write_results(tmp_path / "variants", _variant_results())
    with pytest.raises(SystemExit, match="after the gate"):
        late_proxy.main(["decide", "--stage", "variants", "--runs", str(tmp_path), "--out", str(other)])
    failed = _gate_results({("T1/E8.5->E9.5:cm", "lp_proxy"): [0.0] * 10})
    third = tmp_path / "third"
    gate = _decide(third, monkeypatch, "gate", failed)
    assert not gate["verdict"]["passed"]
    _write_results(third / "variants", _variant_results())
    with pytest.raises(SystemExit, match="did not pass"):
        late_proxy.main(["decide", "--stage", "variants", "--runs", str(third), "--out", str(third)])


def test_variants_after_a_passed_gate_are_read_at_the_gates_commit(tmp_path, monkeypatch, capsys):
    _decide(tmp_path, monkeypatch, "gate", _gate_results())
    out = _decide(tmp_path, monkeypatch, "variants", _variant_results())
    assert {v["variant"]: v["passed"] for v in out["verdict"]["variants"]} == dict.fromkeys(late_proxy.VARIANTS, True)
    capsys.readouterr()
    # a variants run at another commit is refused
    other = tmp_path / "x"
    _decide(other, monkeypatch, "gate", _gate_results())
    _write_results(other / "variants", {k: {**r, "provenance": {**r["provenance"], "commit": "d" * 40}}
                                        for k, r in _variant_results().items()})
    with pytest.raises(SystemExit, match="one commit"):
        late_proxy.main(["decide", "--stage", "variants", "--runs", str(other), "--out", str(other)])


def test_registry_rows_follow_the_registry_header(tmp_path):
    results = _gate_results()
    for r in results.values():
        r["board"] = {"key": "pseudo:x"}
        r["params"] = {}
        r["provenance"]["utc"] = "2026-11-11 00:00:00"
    rows = late_proxy.registry_rows(results, tmp_path / "gate", "gate")
    assert {tuple(sorted(r)) for r in rows} == {tuple(sorted(late_proxy.registry.HEADER))}
    assert len(rows) == 4 * len(late_proxy.READS) and rows[0]["run_id"] == "late_proxy_v1/gate"
    assert late_proxy.registry.append(rows, tmp_path / "registry.csv") == len(rows)


def test_every_read_must_carry_the_consistent_harness_stamp():
    results = _gate_results()
    results["T1/E9.5->E10.5"]["provenance"]["commit"] = "d" * 40
    with pytest.raises(harness.InconsistentResults, match="commits"):
        late_proxy.check_consistent(results, late_proxy.all_tasks())


# ---------------------------------------------------------------------------------------------------------------------
# N20: the proxy's target-free profile and its gate


def _n20_inputs(tmp_path, monkeypatch, stages=STAGES):
    import n20_proxy_profile as n20

    genes = [*GENES, "Hba-syn1", "mt-Syn1"]
    stats_path = _stats(tmp_path / "data" / "external" / "feng2022_gse193346_geo" / "cd1_stats.npz", genes,
                        stages=stages)
    monkeypatch.setattr(n20.feng_dirs, "clean_stats_path", lambda: stats_path)
    monkeypatch.setattr(n20.provenance, "git_state", lambda repo: ("c" * 40, []))
    monkeypatch.setattr(n20.n12, "tries_rows", lambda *a: {"N20"})
    t1_file = tmp_path / "E9.5_RNA.h5ad"
    a = ad.AnnData(X=sparse.csr_matrix(np.zeros((3, len(genes)), np.float32)))
    a.var_names, a.obs_names = genes, ["c0", "c1", "c2"]
    a.write_h5ad(t1_file)
    return n20, ["--t1-file", str(t1_file), "--out", str(tmp_path / "n20")]


def test_n20_profiles_the_proxy_says_proceed_and_refuses_a_rerun(tmp_path, monkeypatch, capsys):
    n20, argv = _n20_inputs(tmp_path, monkeypatch)
    monkeypatch.setattr(n20, "MIN_MAPPED", 40)  # the gate's 20,000-gene floor, scaled to the synthetic panel
    n20.main(argv)
    out = tmp_path / "n20"
    report = json.loads((out / "n20_proxy_profile.json").read_text())
    assert report["proceed_to_tries"] is True and report["gate"]["cm_at_both_ends"] is True
    assert set(report["groups"]) == {"CM", "Endocardium", "Endothelium", "Epicardium"}  # M1's four groups, not Blood
    assert set(report) == {"label", "gate", "proceed_to_tries", "proxy_cells", "real_cells", "groups", "n_t1_genes",
                           "inputs", "stamp"}  # a gene count and Feng's cells, never a T1 cell count
    cm = report["groups"]["CM"]
    assert cm["genes_with_value"] == len(GENES) + 2 and len(cm["eb_factors"]) == harness.LP_BINS
    assert cm["spearman_with_real"] is not None and -1.0 <= cm["spearman_with_real"] <= 1.0
    assert set(cm["modules"]) == {"top_400_up", "top_400_down"}
    assert report["stamp"]["commit"] == "c" * 40 and report["stamp"]["tree_clean_at_end"] is True
    assert (out / "n20_proxy_profile.txt").read_text().startswith("N20: proceed to the tries: YES")
    assert capsys.readouterr().out.startswith("N20: proceed to the tries: YES")
    with pytest.raises(SystemExit, match="already written"):
        n20.main(argv)


def test_n20_needs_its_row_before_it_runs(tmp_path, monkeypatch):
    n20, argv = _n20_inputs(tmp_path, monkeypatch)
    monkeypatch.setattr(n20.n12, "tries_rows", lambda *a: {"T13"})
    with pytest.raises(SystemExit, match="no row"):
        n20.main(argv)


def test_n20_stops_the_tries_when_the_proxy_bracket_is_missing_or_the_genes_too_few(tmp_path, monkeypatch):
    n20, argv = _n20_inputs(tmp_path, monkeypatch, stages=("E9.5", "E14.5"))
    with pytest.raises(SystemExit, match="says no"):
        n20.main(argv)
    assert json.loads((tmp_path / "n20" / "n20_proxy_profile.json").read_text())["gate"]["stages_present"] is False
    other = tmp_path / "more"
    other.mkdir()
    n20, argv = _n20_inputs(other, monkeypatch)  # all stages, but the panel's 62 genes are below the 20,000 floor
    with pytest.raises(SystemExit, match="says no"):
        n20.main(argv)
    assert json.loads((other / "n20" / "n20_proxy_profile.json").read_text())["gate"]["mapped_genes_at_least"] is False
