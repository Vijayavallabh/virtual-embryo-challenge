"""The T1 late-proxy tournament (experiments/harness/TRIES.md T13-T18, N20), on synthetic data only: the bracket
argument of feng_dirs (default today's), the candidates' legality and hand-over, the controls and variants, the
within-type reads, the tournament's decision rule and its refusals."""

from __future__ import annotations

import csv
import json
import subprocess
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
STAGES = ("E9.5", "E14.5", "E18")
LABELS = ("ventricular_cm", "atrial_cm", "endocardial_ec", "vascular_ec", "epicardial", "blood")


# ---------------------------------------------------------------------------------------------------------------------
# Synthetic Feng statistics


def _stats(path: Path, genes: list[str], *, stages=STAGES, seed: int = 3, small: tuple[str, str] | None = None) -> Path:
    """cd1_stats.npz and its MANIFEST.json, in the clean reduction's layout: per stage|zone|label|half, a noisy mean
    CP10k profile with a stage trend (so the directions are not all zero). `small` (label, stage) leaves that label too
    few cells at that stage for a direction."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(1.0, 5.0, len(genes))
    keys, n, cp = [], [], []
    for label in LABELS:
        for i, stage in enumerate(stages):
            trend = np.exp((-1.0) ** (np.arange(len(genes)) % 2) * 0.25 * i)  # odd and even genes move opposite ways
            for half in "AB":
                k = 2 if small == (label, stage) else int(rng.integers(30, 60))
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
    # the late bracket: the same formula across E14.5 -> E18, and every group of the default bracket
    late, rec_late = feng_dirs.feng_directions(list(GENES), stats_path=stats, bracket=("E14.5", "E18"))
    assert rec_late["bracket"] == ["E14.5", "E18"] and set(late) == set(out)
    np.testing.assert_allclose(late["CM"], _manual(stats, "E14.5", "E18", ("ventricular_cm", "atrial_cm")))
    assert feng_dirs.bracket_fraction(10.5, 12.5) == pytest.approx(0.4)
    assert feng_dirs.bracket_fraction(9.5, 10.5) == pytest.approx(0.2)
    assert feng_dirs.bracket_fraction(9.5, 10.5, (14.5, 18.0)) == pytest.approx(2 / 7)


def test_directions_from_stats_is_what_feng_directions_computes_with_no_file_reads(stats):
    st = feng_dirs.dl.load_stats(stats, fields=feng_dirs.dl.STAT_FIELDS)
    cfg = feng_dirs.load_config()
    a, cells, idx = feng_dirs.directions_from_stats(st, list(GENES), cfg, ("E14.5", "E18"))
    b, rec = feng_dirs.feng_directions(list(GENES), cfg, stats, bracket=("E14.5", "E18"))
    assert set(a) == set(b) and cells == rec["groups"] and (idx >= 0).sum() == len(GENES)
    for g in a:
        np.testing.assert_array_equal(a[g], b[g])


def test_bracket_stats_keeps_only_the_two_ends_and_refuses_a_missing_one(stats):
    st = feng_dirs.dl.load_stats(stats, fields=feng_dirs.dl.STAT_FIELDS)
    cut = feng_dirs.bracket_stats(st, ("E14.5", "E18"))
    assert {k.split("|")[0] for k in cut["keys"]} == {"E14.5", "E18"} and len(cut["genes"]) == len(GENES)
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
    h = feng_dirs.direction_halves(st, list(GENES), cfg, ("E14.5", "E18"))["CM"]
    z = np.load(stats, allow_pickle=True)
    keys = [str(k) for k in z["keys"]]

    def mean_cp(stage, halves=("A", "B")):
        m = np.array([k.split("|")[0] == stage and k.split("|")[2] in ("ventricular_cm", "atrial_cm")
                      and k.split("|")[3] in halves for k in keys])
        return z["S_cp"][m].sum(0) / z["n"][m].sum()

    np.testing.assert_allclose(h["d_a"], np.log((mean_cp("E18", ("A",)) + 0.5) / (mean_cp("E14.5", ("A",)) + 0.5)))
    np.testing.assert_allclose(h["d_b"], np.log((mean_cp("E18", ("B",)) + 0.5) / (mean_cp("E14.5", ("B",)) + 0.5)))
    np.testing.assert_allclose(h["mean_cp"], (mean_cp("E14.5") + mean_cp("E18")) / 2.0)


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
        assert c.externals == (("X002", 14.5, 18.0),) and c.params["scope"] == "chamber"
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
    assert {str(k).split("|")[0] for k in bundle["stats"]["keys"]} == {"E14.5", "E18"}  # no other Feng stage
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
    assert kw["direction_fraction"] == pytest.approx(2 / 7)
    assert kw["direction_groups"]["V-CM"] == "CM" and kw["direction_groups"]["Endocardium"] == "Endocardium"
    assert kw.get("start") is None and kw.get("damping", 0.0) == 0.0  # no momentum: the cells without a direction stay
    np.testing.assert_allclose(kw["directions"]["CM"], _manual(stats_path, "E14.5", "E18",
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
    assert calls[0]["direction_fraction"] == pytest.approx(2 / 7)  # E8.5 -> E9.5, one day of the 3.5-day bracket


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
    assert [tuple(e) for e in res["params"]["lp_proxy"]["externals"]] == [("X002", 14.5, 18.0)]


# ---------------------------------------------------------------------------------------------------------------------
# The within-type reads


def test_the_within_type_config_is_three_reads_that_pass_the_rules():
    tasks = load_tasks(late_proxy.CONFIG)["T1:late_proxy"]
    assert {t.key for t in tasks} == set(late_proxy.WITHIN)
    assert all(harness.check_task(t) is t and t.target_filter["within"] in {"cm", "endocardium"} and t.counted
               for t in tasks)
    with pytest.raises(ValueError, match="category"):
        harness.check_task(replace(tasks[0], target_filter={"within": "not_a_category"}))


def test_read_h_is_named_with_its_known_feature_and_is_not_called_held_out():
    """Read H was read with Feng's real direction by the readout's smoke: the config has to say so, in numbers."""
    note = {t.key: t.notes for t in load_tasks(late_proxy.CONFIG)["T1:late_proxy"]}[late_proxy.READ_H]
    assert "not a clean read" in note and "per_type.csv" in note and "0.717" in note and "0.084" in note
    assert "held-out" not in note and "cleanest" not in note
    assert "held-out" not in late_proxy.__doc__ and "not a clean" in late_proxy.__doc__
    assert "sign test" not in late_proxy.__doc__  # seeds check draw noise; multiplicity is the try count's


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


def _res(read: str, scores: dict[str, list[float]], groups: bool = True, commit: str = "c" * 40,
         on: str = "de_recovery") -> dict:
    """A result of `late_proxy.py run` for one read. `on` is the scored group the score's margin sits on."""
    rows = [{"seed": s, "build_seed": harness.BUILD_SEED_OFFSET + s,
             "scores": {c: 50.0 + v[s] for c, v in scores.items()},
             "groups": {c: {g: 0.5 + (v[s] / 100 if g == on else 0.0)
                            for g in ("de_recovery", "direction", "distribution", "covariance")}
                        for c, v in scores.items()} if groups else {}}
            for s in harness.SEEDS]
    return {"task": {"key": read, "test_board": "T1:late_proxy"}, "candidates": list(scores), "rows": rows,
            "protocol": "bake", "build": "main", "anchors": "computed", "frac": 0.1, "seeds": list(harness.SEEDS),
            "board": {"key": f"pseudo:{read}"}, "params": {},
            "provenance": {"commit": commit, "tree_clean": True, "config_sha256": "c" * 64,
                           "utc": "2026-11-11 00:00:00", "files_sha256": {"T1/E9.5_RNA.h5ad": "f"}},
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


def test_the_gate_passes_only_with_a_gain_on_read_h_a_within_type_e10_gain_and_both_controls():
    v = late_proxy.gate_verdict(_gate_results())
    assert v["passed"] and all(v["conditions"].values())
    assert set(v["comparisons"]) == {"lp_proxy vs b0_scope", "lp_proxy vs lp_cflip", "lp_proxy vs lp_cperm",
                                     "lp_cperm vs b0_scope", "lp_cflip vs b0_scope", "lp_proxy vs b0 (reported)",
                                     "lp_proxy vs the floor (reported)"}


@pytest.mark.parametrize("over, condition", [
    ({("T1/E8.5->E9.5:cm", "lp_proxy"): [0.1] * 10}, "read_h"),  # read H must pass
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
def test_the_gate_fails_when_a_control_itself_passes_read_h(control):
    v = late_proxy.gate_verdict(_gate_results({("T1/E8.5->E9.5:cm", control): [0.5] * 10,
                                               ("T1/E8.5->E9.5:cm", "lp_proxy"): [1.5] * 10}))
    # the control passes on H (0.5 > 0.3) and lp_proxy beats it by more than 0.3 there; the gate still fails
    assert v["conditions"]["beats_cperm"] and v["conditions"]["beats_cflip"] and not v["passed"]
    assert not v["conditions"]["controls_do_not_pass_h"]


def test_a_missing_read_fails_its_condition_and_a_complete_set_decides():
    results = _gate_results()
    del results["T1/E8.5->E9.5:cm"]
    v = late_proxy.gate_verdict(results)
    assert not v["passed"] and not v["conditions"]["read_h"] and not v["conditions"]["no_veto"]


@pytest.mark.parametrize("group", ["distribution", "covariance"])
def test_a_gain_that_sits_off_the_rank_metrics_is_a_fit_and_fails_the_gate(group):
    """The row's fit column is a condition of the pass: a gain on mmd_u or the variogram cannot unlock T16-T18."""
    v = late_proxy.gate_verdict(_gate_results(on=group))
    assert not v["passed"] and not v["conditions"]["gain_on_rank_metrics"]
    assert all(ok for name, ok in v["conditions"].items() if name != "gain_on_rank_metrics")  # the fit alone stops it


def test_a_gain_on_direction_counts_as_a_rank_metric_gain():
    assert late_proxy.gate_verdict(_gate_results(on="direction"))["passed"]


@pytest.mark.parametrize("margins, ok", [
    ({"de_recovery": 0.3, "direction": 0.3, "distribution": 0.4, "covariance": 0.0}, True),
    ({"de_recovery": 0.25, "direction": 0.25, "distribution": 0.5, "covariance": 0.0}, True),  # exactly half
    ({"de_recovery": 0.1, "direction": 0.1, "distribution": 0.6, "covariance": 0.2}, False),  # under half
    ({"de_recovery": -0.1, "direction": 0.05, "distribution": -0.9, "covariance": 0.0}, False),  # rank margin <= 0
])
def test_rank_share_needs_a_positive_rank_margin_of_at_least_half_the_total(margins, ok):
    assert late_proxy.rank_share_ok({"complete": True, "group_margins": margins}) is ok


@pytest.mark.parametrize("t", [{"complete": False}, {"complete": True}, {"complete": True, "group_margins": {}}])
def test_rank_share_is_false_for_an_incomplete_read_or_no_margins(t):
    assert late_proxy.rank_share_ok(t) is False


def _variant_results(variant_gains: dict[str, list[float]] | None = None, **kw) -> dict:
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
        out[read] = _res(read, scores, **kw)
    return out


def test_a_variant_needs_nine_of_ten_seeds_ahead_on_read_h_and_both_within_type_e10_reads():
    res = _variant_results()
    for variant, reference in late_proxy.VARIANTS.items():
        v = late_proxy.variant_verdict(res, variant, reference)
        assert v["passed"], variant
    # 8 of 10 seeds ahead passes the adoption rule's test but not the stricter bar for the variants
    gains = [0.6] * 8 + [-0.05, -0.05]
    res = _variant_results({"lp_hyg": gains})
    v = late_proxy.variant_verdict(res, "lp_hyg", "lp_proxy")
    assert not v["passed"] and not v["conditions"]["read_h"]
    assert v["comparison"]["T1/E8.5->E9.5:cm"]["passes"] and not v["comparison"]["T1/E8.5->E9.5:cm"]["passes_strict"]


@pytest.mark.parametrize("read", late_proxy.WITHIN_E10)
def test_a_variant_must_pass_on_both_within_type_e10_reads_not_one(read):
    """The bar binds on something the seeds cannot fake: the gain has to hold on each of the three reads."""
    res = _variant_results()
    for row in res[read]["rows"]:
        row["scores"]["lp_eb"] = row["scores"]["lp_proxy"] + 0.1  # a small gain: no strict pass on this one read
    v = late_proxy.variant_verdict(res, "lp_eb", "lp_proxy")
    assert not v["passed"] and not v["conditions"]["both_within_type_e10"]
    assert v["conditions"]["read_h"] and v["conditions"]["no_veto"] and len(v["passed_within_e10"]) == 1


@pytest.mark.parametrize("group", ["distribution", "covariance"])
def test_a_variants_gain_off_the_rank_metrics_is_a_fit(group):
    v = late_proxy.variant_verdict(_variant_results(on=group), "lp_hyg", "lp_proxy")
    assert not v["passed"] and not v["conditions"]["gain_on_rank_metrics"]
    assert v["conditions"]["read_h"] and v["conditions"]["both_within_type_e10"] and v["conditions"]["no_veto"]


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


def _n20_report(proceed: bool = True, close: bool = True, t18: bool = True, commit: str = "c" * 40) -> dict:
    """What late_proxy.py reads of N20's report."""
    return {"proceed_to_tries": proceed, "proxy_close_to_real": close, "t18_runs": t18,
            "landing_floor": {"group": "CM", "spearman_cross_half": 0.4 if close else 0.05, "floor": 0.2, "met": close},
            "stamp": {"commit": commit, "tree_clean_at_end": True}}


def _ready(monkeypatch, tmp_path, proceed: bool = True, t18: bool = True, commit: str = "c" * 40) -> str:
    """Everything `run` checks before anything else: the rows, X002's amendment, a clean tree, N20's stamped yes and,
    for the gate, a passed preflight at this commit."""
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"N20b", "T13", "T14", "T15", "T16", "T17", "T18"})
    monkeypatch.setattr(feng_dirs, "x002_amended", lambda *a: True)
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: (commit, []))
    n20 = tmp_path / "n20.json"
    n20.write_text(json.dumps(_n20_report(proceed, t18=t18, commit=commit)))
    (tmp_path / "preflight.json").write_text(json.dumps({"passed": True, "reads": {},
                                                         "stamp": {"commit": commit, "tree_clean": True}}))
    return str(n20)


def test_run_needs_the_rows_and_n20s_yes(tmp_path, monkeypatch):
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"T13"})
    monkeypatch.setattr(feng_dirs, "x002_amended", lambda *a: True)
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: ("c" * 40, []))
    n20 = tmp_path / "n20.json"
    argv = ["run", "--stage", "gate", "--out", str(tmp_path), "--n20", str(n20)]
    with pytest.raises(SystemExit, match="not in experiments/harness/TRIES.md"):
        late_proxy.main(argv)
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"N20b", "T13", "T14", "T15"})
    with pytest.raises(SystemExit, match="N20 runs before the gate"):
        late_proxy.main(argv)
    n20.write_text(json.dumps(_n20_report(proceed=False)))
    with pytest.raises(SystemExit, match="said no"):
        late_proxy.main(argv)
    # the variants run only after a passed gate
    n20.write_text(json.dumps(_n20_report()))
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"T16", "T17", "T18"})
    argv = ["run", "--stage", "variants", "--out", str(tmp_path), "--n20", str(n20)]
    with pytest.raises(SystemExit, match="only after the gate passed"):
        late_proxy.main(argv)
    (tmp_path / "decision_gate.json").write_text(json.dumps({"verdict": {"passed": False}}))
    with pytest.raises(SystemExit, match="only after the gate passed"):
        late_proxy.main(argv)


@pytest.mark.parametrize("command", ["preflight", "run"])
def test_a_stage_refuses_without_x002s_amendment_on_a_dirty_tree_or_off_n20s_commit(tmp_path, monkeypatch, command):
    n20 = _ready(monkeypatch, tmp_path)
    argv = [command, *(["--stage", "gate"] if command == "run" else []), "--out", str(tmp_path), "--n20", n20]
    monkeypatch.setattr(feng_dirs, "x002_amended", lambda *a: False)
    with pytest.raises(SystemExit, match="X002 row does not carry the amendment"):
        late_proxy.main(argv)
    monkeypatch.setattr(feng_dirs, "x002_amended", lambda *a: True)
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: ("c" * 40, ["results/registry.csv"]))
    with pytest.raises(SystemExit, match="not clean"):
        late_proxy.main(argv)
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: ("d" * 40, []))
    with pytest.raises(SystemExit, match="N20 must have run clean at this commit"):
        late_proxy.main(argv)
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: ("c" * 40, []))
    report = _n20_report()
    report["stamp"]["tree_clean_at_end"] = False
    Path(n20).write_text(json.dumps(report))
    with pytest.raises(SystemExit, match="N20 must have run clean at this commit"):
        late_proxy.main(argv)


def test_x002_amended_reads_x002s_own_row_only(tmp_path):
    register = tmp_path / "EXTERNAL_DATA.md"
    register.write_text("| X001 | E14.5 -> E18 elsewhere |\n| X002 | Feng, the plain row |\n")
    assert feng_dirs.x002_amended(register) is False
    register.write_text("| X001 | other |\n| X002 | Feng; from T13-T18: Feng's E14.5 -> E18 change is a stand-in |\n")
    assert feng_dirs.x002_amended(register) is True
    register.write_text("| X001 | other |\n")
    assert feng_dirs.x002_amended(register) is False


def test_the_preflight_guards_the_gate_and_a_read_that_cannot_be_scored_stops_it_uncounted(tmp_path, monkeypatch):
    n20 = _ready(monkeypatch, tmp_path)
    (tmp_path / "preflight.json").unlink()
    run = ["run", "--stage", "gate", "--out", str(tmp_path), "--n20", n20]
    with pytest.raises(SystemExit, match="the preflight runs before the gate"):
        late_proxy.main(run)
    seen = []

    def run_task(t, cands, **kw):
        seen.append((t.key, list(cands), kw["seeds"], kw["log"]))
        if t.key.endswith("endocardium"):
            raise ValueError("the floor scores 41.0 (inverted ['mmd_u'], missing anchors [])")
        return {}

    monkeypatch.setattr(late_proxy.harness, "run_task", run_task)
    with pytest.raises(SystemExit, match="not run and do not count"):
        late_proxy.main(["preflight", "--out", str(tmp_path), "--n20", n20])
    assert {k for k, *_ in seen} == set(late_proxy.WITHIN)  # the deciding reads, and only those
    assert all(cands == [] and seeds == tuple(harness.SEEDS) and log is None for _, cands, seeds, log in seen)
    report = json.loads((tmp_path / "preflight.json").read_text())
    assert report["passed"] is False and report["reads"]["T1/E9.5->E10.5:endocardium"].startswith("ValueError")
    assert report["reads"]["T1/E8.5->E9.5:cm"] == "ok" and report["stamp"]["tree_clean"] is True
    with pytest.raises(SystemExit, match="cannot be scored"):
        late_proxy.main(run)
    with pytest.raises(SystemExit, match="never overwritten"):
        late_proxy.main(["preflight", "--out", str(tmp_path), "--n20", n20])


def test_a_preflight_from_another_commit_or_a_dirty_tree_does_not_open_the_gate(tmp_path, monkeypatch):
    n20 = _ready(monkeypatch, tmp_path)
    run = ["run", "--stage", "gate", "--out", str(tmp_path), "--n20", n20]
    for pre_stamp in ({"commit": "d" * 40, "tree_clean": True}, {"commit": "c" * 40, "tree_clean": False}):
        (tmp_path / "preflight.json").write_text(json.dumps({"passed": True, "reads": {}, "stamp": pre_stamp}))
        with pytest.raises(SystemExit, match="must have passed clean at this commit"):
            late_proxy.main(run)


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


def _decide(tmp_path, monkeypatch, stage, results, n20: dict | None = None):
    monkeypatch.setattr(late_proxy, "check_consistent", lambda res, tasks: None)
    monkeypatch.setattr(late_proxy.provenance, "git_state", lambda repo: ("c" * 40, []))
    _write_results(tmp_path / stage, results)
    path = tmp_path / "n20.json"
    if n20 is not None or not path.exists():
        path.write_text(json.dumps(n20 or _n20_report()))
    late_proxy.main(["decide", "--stage", stage, "--runs", str(tmp_path), "--out", str(tmp_path), "--n20", str(path)])
    return json.loads((tmp_path / f"decision_{stage}.json").read_text())


def test_decide_writes_the_gate_and_refuses_variants_before_a_passed_gate(tmp_path, monkeypatch, capsys):
    out = _decide(tmp_path, monkeypatch, "gate", _gate_results())
    assert out["verdict"]["passed"] and out["constants"]["min_ahead"] == 8 and out["constants"]["variant_ahead"] == 9
    assert out["constants"]["min_rank_share"] == 0.5 and out["constants"]["read_h"] == late_proxy.READ_H
    assert out["verdict"]["n20"] == {"proxy_close_to_real": True, "t18_runs": True,
                                     "landing_floor": _n20_report()["landing_floor"]}
    assert "PASSED" in capsys.readouterr().out
    n20 = ["--n20", str(tmp_path / "n20.json")]
    with pytest.raises(SystemExit, match="never overwritten"):
        late_proxy.main(["decide", "--stage", "gate", "--runs", str(tmp_path), "--out", str(tmp_path), *n20])
    # a variants decision without a gate decision is refused; after a failed gate too
    other = tmp_path / "other"
    other.mkdir()
    _write_results(tmp_path / "variants", _variant_results())
    with pytest.raises(SystemExit, match="after the gate"):
        late_proxy.main(["decide", "--stage", "variants", "--runs", str(tmp_path), "--out", str(other), *n20])
    failed = _gate_results({("T1/E8.5->E9.5:cm", "lp_proxy"): [0.0] * 10})
    third = tmp_path / "third"
    gate = _decide(third, monkeypatch, "gate", failed)
    assert not gate["verdict"]["passed"]
    _write_results(third / "variants", _variant_results())
    with pytest.raises(SystemExit, match="did not pass"):
        late_proxy.main(["decide", "--stage", "variants", "--runs", str(third), "--out", str(third), *n20])


def test_decide_refuses_without_n20s_report(tmp_path, monkeypatch):
    monkeypatch.setattr(late_proxy, "check_consistent", lambda res, tasks: None)
    _write_results(tmp_path / "gate", _gate_results())
    with pytest.raises(SystemExit, match="similarity floor"):
        late_proxy.main(["decide", "--stage", "gate", "--runs", str(tmp_path), "--out", str(tmp_path), "--n20",
                         str(tmp_path / "absent.json")])
    assert not (tmp_path / "decision_gate.json").exists()


def test_variants_after_a_passed_gate_are_read_at_the_gates_commit(tmp_path, monkeypatch, capsys):
    _decide(tmp_path, monkeypatch, "gate", _gate_results())
    out = _decide(tmp_path, monkeypatch, "variants", _variant_results())
    assert {v["variant"]: v["passed"] for v in out["verdict"]["variants"]} == dict.fromkeys(late_proxy.VARIANTS, True)
    assert out["verdict"]["not_run"] == {}
    landings = {v["variant"]: v["landing"] for v in out["verdict"]["variants"]}
    assert "first counted support for R2" in landings["lp_hyg"]
    assert "recommended to the captain" in landings["lp_hyg_n"] and "recommended to the captain" in landings["lp_eb"]
    capsys.readouterr()
    # a variants run at another commit is refused
    other = tmp_path / "x"
    _decide(other, monkeypatch, "gate", _gate_results())
    _write_results(other / "variants", {k: {**r, "provenance": {**r["provenance"], "commit": "d" * 40}}
                                        for k, r in _variant_results().items()})
    with pytest.raises(SystemExit, match="one commit"):
        late_proxy.main(["decide", "--stage", "variants", "--runs", str(other), "--out", str(other), "--n20",
                         str(other / "n20.json")])


def test_a_pass_below_n20s_similarity_floor_recommends_nothing_for_m1(tmp_path, monkeypatch, capsys):
    below = _n20_report(close=False)
    _decide(tmp_path, monkeypatch, "gate", _gate_results(), n20=below)
    out = _decide(tmp_path, monkeypatch, "variants", _variant_results())
    assert out["verdict"]["n20"]["proxy_close_to_real"] is False
    for v in out["verdict"]["variants"]:
        assert v["passed"] and "platform evidence only" in v["landing"] and "recommended" not in v["landing"], v
    assert "no recommendation for M1" in capsys.readouterr().out


def test_a_variant_not_passing_lands_nothing():
    assert late_proxy.landing({"variant": "lp_eb", "passed": False}, close=True) == "not passed: nothing lands"


@pytest.mark.parametrize("flag", [False, None])
def test_t18_leaves_the_variants_when_n20s_stop_says_so_and_is_not_counted(tmp_path, monkeypatch, capsys, flag):
    n20 = _ready(monkeypatch, tmp_path, t18=bool(flag))
    if flag is None:  # a report with no flag is refused the same way: nothing lets T18 run by default
        report = _n20_report()
        del report["t18_runs"]
        Path(n20).write_text(json.dumps(report))
    (tmp_path / "decision_gate.json").write_text(json.dumps({"verdict": {"passed": True}}))
    built = []

    def run_task(t, cands, **kw):
        built.append([c.name for c in cands])
        return _variant_results()[t.key]

    monkeypatch.setattr(late_proxy.harness, "run_task", run_task)
    monkeypatch.setattr(late_proxy, "input_hashes", lambda: {"feng_stats_sha256": "a" * 64})
    late_proxy.main(["run", "--stage", "variants", "--out", str(tmp_path), "--task", "T1/E8.5->E9.5:cm", "--n20", n20])
    assert built == [[harness.FLOOR, "b0_scope", "lp_proxy", "lp_hyg", "lp_hyg_n"]]  # lp_eb is not built


def test_t18_not_run_is_recorded_in_the_decision_and_is_not_read(tmp_path, monkeypatch, capsys):
    stop = _n20_report(t18=False)
    _decide(tmp_path, monkeypatch, "gate", _gate_results(), n20=stop)
    out = _decide(tmp_path, monkeypatch, "variants", _variant_results())
    assert [v["variant"] for v in out["verdict"]["variants"]] == ["lp_hyg", "lp_hyg_n"]
    assert set(out["verdict"]["not_run"]) == {"lp_eb"} and out["verdict"]["n20"]["t18_runs"] is False
    assert "lp_eb not run" in capsys.readouterr().out
    assert late_proxy.stage_candidates("gate", stop) == late_proxy.CANDIDATES["gate"]  # its stop leaves the gate alone


def test_registry_rows_follow_the_registry_header(tmp_path):
    rows = late_proxy.registry_rows(_gate_results(), tmp_path / "gate", "gate")
    assert {tuple(sorted(r)) for r in rows} == {tuple(sorted(late_proxy.registry.HEADER))}
    assert len(rows) == 4 * len(late_proxy.READS) and rows[0]["run_id"] == "late_proxy_v1/gate"
    assert late_proxy.registry.append(rows, tmp_path / "registry.csv") == len(rows)


def test_every_read_must_carry_the_consistent_harness_stamp():
    results = _gate_results()
    results["T1/E9.5->E10.5"]["provenance"]["commit"] = "d" * 40
    with pytest.raises(harness.InconsistentResults, match="commits"):
        late_proxy.check_consistent(results, late_proxy.all_tasks())
    results = _gate_results()
    results["T1/E9.5->E10.5"]["provenance"]["tree_clean"] = False
    with pytest.raises(harness.InconsistentResults, match="dirty tree"):
        late_proxy.check_consistent(results, late_proxy.all_tasks())


def test_register_appends_a_decided_stage_once(tmp_path, monkeypatch):
    reg = tmp_path / "registry.csv"
    argv = ["register", "--stage", "gate", "--runs", str(tmp_path), "--out", str(tmp_path), "--registry", str(reg)]
    monkeypatch.setattr(late_proxy, "check_consistent", lambda res, tasks: None)
    _write_results(tmp_path / "gate", _gate_results())
    with pytest.raises(SystemExit, match="registered after its decision"):
        late_proxy.main(argv)
    (tmp_path / "decision_gate.json").write_text("{}")
    late_proxy.main(argv)
    with reg.open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 4 * len(late_proxy.READS) and {r["run_id"] for r in rows} == {"late_proxy_v1/gate"}
    with pytest.raises(SystemExit, match="registered once"):
        late_proxy.main(argv)
    assert len(reg.read_text().splitlines()) == len(rows) + 1  # the refusal appended nothing


# ---------------------------------------------------------------------------------------------------------------------
# The command order, in a real repository (registry.csv tracked, check_consistent not mocked)


def _repo(tmp_path) -> tuple[Path, Path, str]:
    """A git repository with a tracked results/registry.csv (the header alone); returns it, the registry, its HEAD."""
    repo = tmp_path / "repo"
    (repo / "results").mkdir(parents=True)
    reg = repo / "results" / "registry.csv"
    reg.write_text(",".join(late_proxy.registry.HEADER) + "\n")
    for cmd in (["init", "-q"], ["add", "results/registry.csv"],
                ["-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "registry"]):
        subprocess.run(["git", *cmd], cwd=repo, check=True, capture_output=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, check=True, capture_output=True, text=True)
    return repo, reg, head.stdout.strip()


def _in_repo(monkeypatch, tmp_path, repo: Path, head: str) -> tuple[Path, str]:
    """late_proxy.py pointed at the temporary repository: run_task records the state harness._git sees there, exactly
    as vec.harness.run_task does, and every other input of the checks is in place. check_consistent is the real one."""
    def git_there(*args):
        return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=True).stdout.strip()

    def run_task(t, cands, **kw):
        names = [c.name for c in cands]
        if not names:  # the preflight: the floor and anchors only
            return {}
        res = (_variant_results() if "lp_hyg" in names else _gate_results())[t.key]
        res["provenance"].update(commit=harness._git("rev-parse", "HEAD"),
                                 tree_clean=harness._git("status", "--porcelain", "--untracked-files=no") == "")
        return res

    monkeypatch.setattr(harness, "_git", git_there)
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(late_proxy.provenance, "REPO_ROOT", repo)
    monkeypatch.setattr(late_proxy, "registered_rows", lambda: {"N20b", "T13", "T14", "T15", "T16", "T17", "T18"})
    monkeypatch.setattr(feng_dirs, "x002_amended", lambda *a: True)
    monkeypatch.setattr(late_proxy, "input_hashes", lambda: {"feng_stats_sha256": "a" * 64})
    out = repo / "experiments" / "t1_m1" / "results" / "late_proxy_v1"  # untracked and exempt, as on the DGX
    n20 = repo / "experiments" / "t1_m1" / "results" / "n20b_proxy_profile" / "n20_proxy_profile.json"
    n20.parent.mkdir(parents=True)
    n20.write_text(json.dumps(_n20_report(commit=head)))
    return out, str(n20)


def test_the_tournament_runs_through_a_passed_gate_in_a_real_repository_and_registers_last(tmp_path, monkeypatch):
    repo, reg, head = _repo(tmp_path)
    out, n20 = _in_repo(monkeypatch, tmp_path, repo, head)
    common = ["--out", str(out), "--n20", n20]
    late_proxy.main(["preflight", *common])
    late_proxy.main(["run", "--stage", "gate", *common])
    late_proxy.main(["decide", "--stage", "gate", "--runs", str(out), *common])
    assert harness._git("status", "--porcelain", "--untracked-files=no") == ""  # the decision changed no tracked file
    late_proxy.main(["run", "--stage", "variants", *common])  # so the tree is still clean for the variants
    late_proxy.main(["decide", "--stage", "variants", "--runs", str(out), *common])  # the real check_consistent passes
    decision = json.loads((out / "decision_variants.json").read_text())
    assert decision["stamp"]["reads_commit"] == head and all(v["passed"] for v in decision["verdict"]["variants"])
    late_proxy.main(["register", "--stage", "gate", "--runs", str(out), "--out", str(out), "--registry", str(reg)])
    late_proxy.main(["register", "--stage", "variants", "--runs", str(out), "--out", str(out), "--registry", str(reg)])
    with reg.open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert {r["run_id"] for r in rows} == {"late_proxy_v1/gate", "late_proxy_v1/variants"}
    assert {r["commit"] for r in rows} == {head} and {r["tree_clean"] for r in rows} == {"True"}
    assert harness._git("status", "--porcelain", "--untracked-files=no") != ""  # the tracked registry changed last


def test_registering_the_gate_before_the_variants_is_refused_before_any_work(tmp_path, monkeypatch):
    """The trap in the old command order: a tracked registry.csv change before the variants run makes them dirty."""
    repo, reg, head = _repo(tmp_path)
    out, n20 = _in_repo(monkeypatch, tmp_path, repo, head)
    common = ["--out", str(out), "--n20", n20]
    late_proxy.main(["preflight", *common])
    late_proxy.main(["run", "--stage", "gate", *common])
    late_proxy.main(["decide", "--stage", "gate", "--runs", str(out), *common])
    late_proxy.main(["register", "--stage", "gate", "--runs", str(out), "--out", str(out), "--registry", str(reg)])
    with pytest.raises(SystemExit, match="not clean"):
        late_proxy.main(["run", "--stage", "variants", *common])
    assert not (out / "variants").exists()  # refused before a single read ran


# ---------------------------------------------------------------------------------------------------------------------
# N20: the proxy's target-free profile and its gate


def _n20_inputs(tmp_path, monkeypatch, stages=STAGES, small=None):
    import n20_proxy_profile as n20

    genes = [*GENES, "Hba-syn1", "mt-Syn1"]
    stats_path = _stats(tmp_path / "data" / "external" / "feng2022_gse193346_geo" / "cd1_stats.npz", genes,
                        stages=stages, small=small)
    monkeypatch.setattr(n20.feng_dirs, "clean_stats_path", lambda: stats_path)
    monkeypatch.setattr(n20.feng_dirs, "x002_amended", lambda *a: True)
    monkeypatch.setattr(n20.provenance, "git_state", lambda repo: ("c" * 40, []))
    monkeypatch.setattr(n20.n12, "tries_rows", lambda *a: {"N20b"})
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
    assert report["gate"]["every_read_has_a_direction"] is True
    assert report["gate"]["direction_available"] == dict.fromkeys(late_proxy.WITHIN, True)
    assert set(report["groups"]) == {"CM", "Endocardium", "Endothelium", "Epicardium"}  # M1's four groups, not Blood
    assert set(report) == {"label", "gate", "proceed_to_tries", "landing_floor", "proxy_close_to_real", "t18_stop",
                           "t18_runs", "proxy_cells", "real_cells", "groups", "n_t1_genes",
                           "inputs", "stamp"}  # a gene count and Feng's cells, never a T1 cell count
    assert report["landing_floor"]["floor"] == 0.2 and report["landing_floor"]["group"] == "CM"
    # the synthetic directions share the same odd/even trend, so the cross-half correlation is strongly positive
    assert report["landing_floor"]["spearman_cross_half"] > 0.2 and report["proxy_close_to_real"] is True
    assert set(report["t18_stop"]["below_no_op"]) == {"CM", "Endocardium"}
    assert report["t18_runs"] is report["t18_stop"]["runs"]
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


def test_n20_refuses_before_reading_anything_without_x002s_amendment_or_on_a_dirty_tree(tmp_path, monkeypatch):
    n20, argv = _n20_inputs(tmp_path, monkeypatch)
    monkeypatch.setattr(n20.feng_dirs, "x002_amended", lambda *a: False)
    monkeypatch.setattr(n20.feng_dirs, "load_config", lambda *a: pytest.fail("N20 read past its guard"))
    with pytest.raises(SystemExit, match="X002 row does not carry the amendment"):
        n20.main(argv)
    monkeypatch.setattr(n20.feng_dirs, "x002_amended", lambda *a: True)
    for state in (("c" * 40, ["src/x.py"]), ("unknown", None)):
        monkeypatch.setattr(n20.provenance, "git_state", lambda repo, state=state: state)
        with pytest.raises(SystemExit, match="not clean"):
            n20.main(argv)
    assert not (tmp_path / "n20").exists()  # nothing was written, so the one allowed run is not spent


def test_n20s_read_groups_are_the_tournaments_within_type_reads(tmp_path, monkeypatch):
    n20, _ = _n20_inputs(tmp_path, monkeypatch)
    assert set(n20.READ_GROUPS) == set(late_proxy.WITHIN)
    assert set(n20.READ_GROUPS.values()) <= set(harness.LP_GROUPS)


def test_n20_stops_the_tries_when_the_endocardium_read_has_no_direction(tmp_path, monkeypatch):
    n20, argv = _n20_inputs(tmp_path, monkeypatch, small=("endocardial_ec", "E18"))
    monkeypatch.setattr(n20, "MIN_MAPPED", 40)
    with pytest.raises(SystemExit, match="says no"):
        n20.main(argv)
    report = json.loads((tmp_path / "n20" / "n20_proxy_profile.json").read_text())
    assert report["gate"]["cm_at_both_ends"] is True and report["gate"]["every_read_has_a_direction"] is False
    assert report["gate"]["direction_available"]["T1/E9.5->E10.5:endocardium"] is False
    assert report["gate"]["direction_available"]["T1/E8.5->E9.5:cm"] is True
    assert report["proceed_to_tries"] is False and report["t18_runs"] is False  # no Endocardium factors, no T18


def test_n20_landing_floor_is_recorded_not_met_above_the_correlation(tmp_path, monkeypatch):
    n20, argv = _n20_inputs(tmp_path, monkeypatch)
    monkeypatch.setattr(n20, "MIN_MAPPED", 40)
    monkeypatch.setattr(n20, "LANDING_FLOOR", 0.9999)
    n20.main(argv)
    report = json.loads((tmp_path / "n20" / "n20_proxy_profile.json").read_text())
    assert report["proceed_to_tries"] is True  # the floor conditions the landing, not the tries
    assert report["landing_floor"]["met"] is False and report["proxy_close_to_real"] is False
    assert report["landing_floor"]["spearman_cross_half"] is not None


def test_n20_landing_floor_is_the_cross_half_correlation_of_the_two_directions(tmp_path, monkeypatch):
    """The two directions share E14.5 with opposite signs, so the floor correlates disjoint halves of Feng's cells: the
    proxy from half A against the real direction from half B, and the reverse, their mean."""
    n20, _ = _n20_inputs(tmp_path, monkeypatch)
    st = n20.dl.load_stats(n20.feng_dirs.clean_stats_path(), fields=n20.dl.STAT_FIELDS)
    genes, cfg = [str(g) for g in st["genes"]], n20.feng_dirs.load_config()
    bracket = n20.feng_dirs.bracket_stats(st, n20.PROXY)
    p = n20.feng_dirs.direction_halves(bracket, genes, cfg, n20.PROXY, groups={"CM"})["CM"]
    r = n20.feng_dirs.direction_halves(st, genes, cfg, n20.REAL, groups={"CM"})["CM"]
    ok = np.isfinite(p["d_a"]) & np.isfinite(p["d_b"]) & np.isfinite(r["d_a"]) & np.isfinite(r["d_b"])
    want = (n20.n12.spearman(p["d_a"][ok], r["d_b"][ok]) + n20.n12.spearman(p["d_b"][ok], r["d_a"][ok])) / 2.0
    rec = n20.landing_floor(st, bracket, genes, cfg)
    assert rec["group"] == "CM" and rec["floor"] == 0.2
    assert rec["spearman_cross_half"] == pytest.approx(want) and rec["met"] is bool(want >= 0.2)
    assert n20.landing_floor(st, None, genes, cfg) == {"group": "CM", "spearman_cross_half": None, "floor": 0.2,
                                                       "met": False}  # no proxy bracket, nothing to read it from


@pytest.mark.parametrize("factors, runs", [
    ({"CM": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95],  # decile 4 below 0.9, expressed deciles shrink
      "Endocardium": [0.5] * 10}, True),
    ({"CM": [0.1] * 3 + [0.95] * 7, "Endocardium": [0.5] * 10}, False),  # CM barely moves above decile 3: cannot pass
    ({"CM": [0.5] * 10, "Endocardium": [0.1] * 3 + [0.9] * 7}, False),  # the same for Endocardium; 0.9 is no shrinkage
    ({"CM": [0.5] * 10}, False),  # a group without factors stops T18
])
def test_t18s_stop_reads_the_expressed_deciles_of_cm_and_endocardium(tmp_path, monkeypatch, factors, runs):
    n20, _ = _n20_inputs(tmp_path, monkeypatch)
    stop = n20.t18_stop({g: {"eb_factors": f} for g, f in factors.items()})
    assert stop["runs"] is runs and stop["no_op_at"] == 0.9 and stop["deciles"] == "4-10"


def test_n20_stops_the_tries_when_the_proxy_bracket_is_missing_or_the_genes_too_few(tmp_path, monkeypatch):
    n20, argv = _n20_inputs(tmp_path, monkeypatch, stages=("E9.5", "E14.5"))
    with pytest.raises(SystemExit, match="says no"):
        n20.main(argv)
    assert json.loads((tmp_path / "n20" / "n20_proxy_profile.json").read_text())["gate"]["stages_present"] is False
    report = json.loads((tmp_path / "n20" / "n20_proxy_profile.json").read_text())
    assert report["landing_floor"]["met"] is False and report["t18_runs"] is False  # nothing to read them from
    other = tmp_path / "more"
    other.mkdir()
    n20, argv = _n20_inputs(other, monkeypatch)  # all stages, but the panel's 62 genes are below the 20,000 floor
    with pytest.raises(SystemExit, match="says no"):
        n20.main(argv)
    assert json.loads((other / "n20" / "n20_proxy_profile.json").read_text())["gate"]["mapped_genes_at_least"] is False
