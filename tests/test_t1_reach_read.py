"""experiments/t1_m1/reach_read.py (try T60 of experiments/harness/TRIES.md) on synthetic stages and hand-made task
files: the split, the projection, the arms, the control, the guard (with its dry run on committed reproduction files)
and the reading. No released file, statistics file or board file is read; the guard's dry run reads committed practice
results only."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest
from test_t1_scope import GENES, HELD
from test_t1_state_transport import states_stage
from test_t1_step_read import _arms_world

from vec import harness, t1_scope
from vec.boards import Board
from vec.core import counts_of
from vec.evaluate import score_board
from vec.io import dense
from vec.t1_scope import classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "reach_read.py"
SCOPE_RESULTS = ROOT / "experiments" / "t1_scope" / "results"
SIZES = {"V-CM": 300, "Endocardium": 120, "Proepicardium": 60, "BEC": 50, "Blood": 40, "Foregut": 80}
SHARES = {"cm": 0.5, "endocardium": 0.2, "epicardium": 0.1, "vascular_endothelium": 0.1, "blood": 0.1}
FENG = ("CM", "Endocardium", "Endothelium", "Epicardium")


@pytest.fixture
def rr(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("reach_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod.gv, "stamp", lambda started, files: dict(STAMP))
    return mod


STAMP = {"commit": "f" * 40, "tree_clean": True, "utc": "x", "inputs": {}}


def _world(rr, seed=3, n=400):
    """An E9.5-like source, an E10.5-like half B with a real change, a base at a chamber mix, M1's groups, four random
    Feng directions (a few genes without a value) and the shift's keywords as the run passes them."""
    spec = load_spec()
    start = states_stage(11, SIZES, later=True)
    classify(start, source="E9.5_RNA")
    sib = states_stage(12, SIZES, later=True, target=True)
    classify(sib, source="E10.5_RNA")
    within = {"surrounding": list(rr.CONTAMINATION_LABELS)}
    base, _ = reweight(start, SHARES, n, seed, within=within, fallback=spec.fallbacks("E9.5_RNA"))
    groups = rr.feng_dirs.label_groups(start.obs["celltype"].astype(str).unique())
    rng = np.random.default_rng(0)
    directions = {g: rng.normal(0.0, 0.5, len(GENES)) for g in FENG}
    for d in directions.values():
        d[[3, 17]] = np.nan
    kw = {"within": within, "start_fallback": spec.fallbacks("E9.5_RNA"), "ref_fallback": spec.fallbacks("E10.5_RNA"),
          "min_cells": rr.MIN_CELLS, "eps": rr.EPS}
    return start, sib, base, groups, directions, kw, list(spec.categories)


# ---------------------------------------------------------------------------------------------------------------------
# The split


def test_the_split_is_row_t27s_rule_for_any_library(rr):
    names = [f"c{i}_{lib}" for i, lib in enumerate("3121332133123" * 7)]
    for lib in ("1", "2", "3"):
        a, b = rr.split_library(names, lib)
        rows = np.flatnonzero(rr.pt.library_of(names) == lib)
        perm = np.random.default_rng(20261112).permutation(rows)
        assert np.array_equal(a, np.sort(perm[: len(perm) // 2])) and np.array_equal(b, np.sort(perm[len(perm) // 2:]))
        assert not set(a) & set(b) and set(a) | set(b) == set(rows)
        assert len(a) == len(rows) // 2
    assert rr.SPLIT_SEED == 20261112  # kind_gap.SPLIT_SEED, row T27's
    a2, _ = rr.split_library(names, "3")
    assert np.array_equal(a2, rr.split_library(names, "3")[0])


# ---------------------------------------------------------------------------------------------------------------------
# The projection


def test_the_projection_is_the_least_squares_multiple_over_the_moved_genes(rr):
    rng = np.random.default_rng(1)
    d = rng.normal(size=20)
    d[4] = np.nan
    held = np.zeros(20, dtype=bool)
    held[[0, 1]] = True
    moved = rr.moved_genes(d, held)
    assert not moved[[0, 1, 4]].any() and moved.sum() == 17
    logf = 2.5 * np.nan_to_num(d) + 0.0
    logf[[0, 1, 4]] = [3.0, -2.0, 1.0]  # off the moved genes: never reached
    out, rec = rr.project(logf, d, moved)
    assert rec["coefficient"] == pytest.approx(2.5) and rec["explained_share"] == pytest.approx(1.0)
    assert np.all(out[~moved] == 0.0) and np.allclose(out[moved], 2.5 * d[moved])
    assert rec["size2_off_moved"] == pytest.approx(9.0 + 4.0 + 1.0)
    noise = rng.normal(size=20)
    out, rec = rr.project(noise, d, moved)
    dm = d[moved]
    assert rec["coefficient"] == pytest.approx(noise[moved] @ dm / (dm @ dm))
    assert 0.0 <= rec["explained_share"] <= 1.0
    out, rec = rr.project(noise, np.zeros(20), moved)  # a zero direction reaches nothing
    assert rec["coefficient"] == 0.0 and np.all(out == 0.0) and rec["explained_share"] == 0.0


def test_a_type_takes_one_direction_or_none_and_a_mixed_type_is_refused(rr):
    start, _, base, groups, directions, _, _ = _world(rr)
    tg = rr.type_groups(base, directions, groups)
    assert tg == {"blood": None, "cm": "CM", "endocardium": "Endocardium", "epicardium": "Epicardium",
                  "vascular_endothelium": "Endothelium"}
    mixed = base.copy()
    labels = mixed.obs["celltype"].astype(str).to_numpy()
    cm = np.flatnonzero(mixed.obs["category"].astype(str).to_numpy() == "cm")
    labels[cm[0]] = "Blood"
    mixed.obs["celltype"] = labels
    with pytest.raises(ValueError, match="not one direction"):
        rr.type_groups(mixed, directions, groups)


# ---------------------------------------------------------------------------------------------------------------------
# The arms


def test_the_capture_filter_leaves_the_oracle_shift_as_pseudo_tasks_shift_1(rr):
    start, sib, base, _, _, kw, _ = _world(rr)
    plain, rep = t1_scope.momentum_shift(start, sib, base, 1.0, HELD, seed=[3, 1], **kw)
    caught, rep2, logf = rr.capture_shift(start, sib, base, HELD, [3, 1], **kw)
    assert (counts_of(plain) != counts_of(caught)).nnz == 0
    assert set(logf) == {r["type"] for r in rep.rows if r["status"] == "shifted"} and len(logf) == 5
    held = np.isin(GENES, HELD)
    for v in logf.values():
        assert np.all(v[held] == 0.0) and np.abs(v).max() > 0.1


_FACTORS: dict = {}


def _factors(start):
    """g_mom's factors on the synthetic source and an E8.5-like start, fitted once for the module."""
    if "f" not in _FACTORS:
        early = states_stage(13, {"V-CM": 200, "Endothelium": 90, "Blood": 40, "Foregut": 60})
        _FACTORS["f"] = harness.g_mom_pair_factors(start, "E9.5_RNA", early, "E8.5_RNA")[0]
    return _FACTORS["f"]


TYPE_GROUP = {"cm": "CM", "endocardium": "Endocardium", "epicardium": "Epicardium",
              "vascular_endothelium": "Endothelium", "blood": "Blood"}


def test_proj_shift_is_oracle_shift_when_the_direction_is_the_oracle_shift_itself(rr):
    """The guard of the brief, through the run's own arms: with each type's direction replaced by its own oracle log
    shift (blood given one too) and M1's held set the shift's, the projection is the shift, coefficient 1, and
    proj_shift equals oracle_shift count for count."""
    start, sib, base, groups, _, kw, order = _world(rr)
    _, _, logf = rr.capture_shift(start, sib, base, HELD, [3, 1], **kw)
    directions = {TYPE_GROUP[t]: v.copy() for t, v in logf.items()}
    arms, rec = rr.split_arms(base, start, sib, directions, groups, HELD, HELD, _factors(start), 3, order, kw)
    for t in logf:
        assert rec[t]["coefficient"] == 1.0 and rec[t]["explained_share"] == 1.0
    assert (counts_of(arms["proj_shift"]) != counts_of(arms["oracle_shift"])).nnz == 0
    assert np.array_equal(dense(arms["proj_shift"].X), dense(arms["oracle_shift"].X))
    proj, _, _ = rr.projections(logf, rr.type_groups(base, directions, groups), directions, np.isin(GENES, HELD),
                                start, order)
    assert all(np.array_equal(proj[t], logf[t]) for t in logf)


def test_proj_shift_is_oracle_mix_when_the_direction_is_zero(rr):
    start, sib, base, groups, directions, kw, order = _world(rr)
    zero = {g: np.zeros(len(GENES)) for g in directions}
    arms, rec = rr.split_arms(base, start, sib, zero, groups, HELD, [*HELD, "g5"], _factors(start), 3, order, kw)
    for k in ("proj_shift", "proj_perm"):
        assert (counts_of(arms[k]) != counts_of(base)).nnz == 0, k
        assert np.allclose(dense(arms[k].X), dense(base.X), atol=1e-6), k
    assert all(v["coefficient"] in (0.0, None) for v in rec.values())


def test_the_five_arms_share_the_draw_and_move_as_designed(rr):
    start, sib, base, groups, directions, kw, order = _world(rr)
    held_m1 = [*HELD, "g5", "g6"]
    arms, rec = rr.split_arms(base, start, sib, directions, groups, HELD, held_m1, _factors(start), 3, order, kw)
    assert list(arms) == ["oracle_mix", "oracle_shift", "proj_shift", "proj_perm", "entry_step"]
    assert arms["oracle_mix"] is base
    for name, a in arms.items():  # the same drawn cells, in the same order
        assert list(a.obs["source_obs"]) == list(base.obs["source_obs"]), name
    C = {k: counts_of(v) for k, v in arms.items()}
    for k in ("oracle_shift", "proj_shift", "proj_perm", "entry_step"):
        assert (C[k] != C["oracle_mix"]).nnz > 0, k
    assert (C["proj_shift"] != C["proj_perm"]).nnz > 0
    # held genes of M1 (R2) and blood (R1) are unmoved by the projection arms; the oracle shift moves both
    g56 = np.isin(GENES, ["g5", "g6"])
    blood = base.obs["category"].astype(str).to_numpy() == "blood"
    for k in ("proj_shift", "proj_perm"):  # R2's genes exactly as drawn, outside the depth constant too
        assert (C[k][:, g56] != C["oracle_mix"][:, g56]).nnz == 0
        assert (C[k][blood] != C["oracle_mix"][blood]).nnz == 0  # R1: blood has no Feng group and stays as drawn
    assert (C["oracle_shift"][:, g56] != C["oracle_mix"][:, g56]).nnz > 0
    assert rec["blood"]["group"] is None and rec["cm"]["group"] == "CM"
    for t in ("cm", "endocardium", "epicardium", "vascular_endothelium"):
        assert np.isfinite(rec[t]["coefficient"]) and rec[t]["coefficient"] != rec[t]["coefficient_perm"]
    assert (C["oracle_shift"][blood] != C["oracle_mix"][blood]).nnz > 0


def test_the_control_permutes_the_same_values_within_deciles_and_differs_from_its_arm(rr):
    start, _, base, groups, directions, _, order = _world(rr)
    d = directions["CM"]
    held = np.isin(GENES, HELD)
    moved = rr.moved_genes(d, held)
    basis = rr.type_basis(start, "cm")
    dp = rr.permuted_direction(d, basis, moved, rr.type_seed("cm", order))
    assert np.array_equal(dp[~moved], d[~moved], equal_nan=True)  # held and no-value genes keep theirs
    assert np.array_equal(np.sort(dp[moved]), np.sort(d[moved]))  # the same values: matched in size
    assert np.linalg.norm(dp[moved]) == pytest.approx(np.linalg.norm(d[moved]))
    idx = np.flatnonzero(moved)
    for run in np.array_split(idx[np.argsort(basis[idx], kind="stable")], rr.PERM_BINS):  # within each decile
        assert np.array_equal(np.sort(dp[run]), np.sort(d[run]))
    assert not np.array_equal(dp, d, equal_nan=True)
    again = rr.permuted_direction(d, basis, moved, rr.type_seed("cm", order))
    assert np.array_equal(again, dp, equal_nan=True)  # one permutation per type, fixed by its seed
    assert rr.type_seed("cm", order) != rr.type_seed("endocardium", order)
    assert rr.type_seed("cm", order)[0] == rr.PERM_SEED == 20261119


def test_entry_step_is_m1_r_gmom_on_the_same_base(rr):
    """On B0's own draw, entry_arm with g_mom_pair_factors is gmom_combo.m1_r_arms's m1_r_gmom bit for bit: the arm is
    the entry's form, only its base differs in the run."""
    ctx, groups, dirs, held = _arms_world()
    want = rr.gc.m1_r_arms(ctx, dirs, groups, 0.2, held, held)["m1_r_gmom"]
    role = ctx.task.reference_role
    src, start_src = ctx.task.role_stems[role], ctx.task.role_stems["start"]
    spec = load_spec()
    ref = harness._shallow(ctx.inputs[role])
    classify(ref, source=src, spec=spec)
    shares = harness.t1_chamber_shares(ref, src)
    base, _ = reweight(ref, shares, harness.count_for(ctx.board, ref.n_obs), ctx.seed, fallback=spec.fallbacks(src),
                       spec=spec)
    factors = harness.g_mom_pair_factors(ctx.inputs[role], src, ctx.inputs["start"], start_src)[0]
    got = rr.entry_arm(base, dirs, groups, 0.2, held, held, factors, ctx.seed)
    X = dense(got[:, ctx.board.genes].X).astype(np.float32)
    assert np.array_equal(np.clip(X, 0.0, None), want.X)
    assert rr.FRACTION == pytest.approx(0.2)


def test_every_arm_is_scored_against_the_first_arms_anchors(rr):
    """score_arms computes one seed's anchors on its first arm and scores every arm against them, which gives each arm
    the score score_board gives it with anchors computed on the same draws (pseudo_task's scoring)."""
    start, sib, base, *_ = _world(rr)
    board = Board(key="T1:synthetic", task="T1", genes=list(GENES), needs_coords=False, min_cells=10, max_cells=2000)
    arms = {"first": rr.pt.to_stage(start[np.arange(0, 600, 2)]), "second": rr.pt.to_stage(base)}
    truth, ref = rr.pt.to_stage(sib), rr.pt.to_stage(start)
    res = rr.score_arms(board, arms, truth, ref, 2)
    for name, pred in arms.items():
        direct = score_board(board, pred, truth, ref, seed=2)
        assert res["arms"][name]["score"] == pytest.approx(direct.score, abs=1e-9), name
    assert set(res["arms"]["second"]["groups"]) == {"de_recovery", "direction", "distribution", "covariance"}


# ---------------------------------------------------------------------------------------------------------------------
# The guard


def _scores(base=60.0, shift=0.0, seeds=range(10), per_seed=None):
    out = {}
    for i, a in enumerate(("copy_last", "oracle", "momentum_0.5")):
        out[a] = {s: base + i + shift + (per_seed[s] if per_seed else 0.0) for s in seeds}
    return out


def test_the_guard_holds_inside_its_tolerances_and_fails_outside(rr):
    c = _scores()
    assert rr.guard(_scores(), c)["holds"]
    assert rr.guard(_scores(shift=0.099), c)["holds"]
    assert not rr.guard(_scores(shift=0.11), c)["holds"]
    one = [0.0] * 10
    one[4] = 0.49
    assert rr.guard(_scores(per_seed=one), c)["holds"]
    one[4] = 0.51  # mean 0.051, inside; one seed outside
    g = rr.guard(_scores(per_seed=one), c)
    assert not g["holds"] and g["arms"]["oracle"]["largest_seed_difference"] == pytest.approx(0.51)
    assert not rr.guard(_scores(seeds=range(9)), c)["holds"]  # a missing seed fails
    assert not rr.guard({}, c)["holds"] and not rr.guard(c, {})["holds"]
    noise = np.random.default_rng(0).normal(0, 1.0, 10)  # pure noise of a point per seed fails the seed tolerance
    assert not rr.guard(_scores(per_seed=dict(enumerate(noise))), c)["holds"]
    assert rr.GUARD_MEAN_TOL == 0.1 and rr.GUARD_SEED_TOL == 0.5


def test_the_guards_dry_run_on_committed_reproduction_files(rr):
    """Rule 9, J7: the guard's arithmetic on committed reproduction data other than the pair the row reads. O1's
    reproduction gate reran copy_last, oracle and shift_1 of the E8.5 -> E9.5 runs at n_cells 5,000 under new tags
    (experiments/t1_scope/results/o1_repro_lib1, _lib2): the guard holds there with every difference zero, so it could
    not have failed there; on two different tasks (T28's library 1 against its library 2) it fails, so it can fail."""
    arms = ("copy_last", "oracle", "shift_1")
    for lib in ("1", "2"):
        fresh = rr.committed_scores(SCOPE_RESULTS / f"o1_repro_lib{lib}" / "scores.csv", arms)
        old = rr.committed_scores(SCOPE_RESULTS / f"pseudo_E8.5_E9.5_lib{lib}" / "scores.csv", arms)
        g = rr.guard(fresh, old)
        assert g["holds"] and max(v["largest_seed_difference"] for v in g["arms"].values()) < 1e-9
    t28 = {lib: rr.committed_scores(rr.t28_path(lib)) for lib in rr.LIBRARIES}
    for lib, s in t28.items():  # T28's files hold every guard arm on every seed
        assert all(sorted(s[a]) == list(range(10)) for a in rr.GUARD_ARMS), lib
    assert not rr.guard(t28["1"], t28["2"])["holds"]


# ---------------------------------------------------------------------------------------------------------------------
# The reading


def _arm(score):
    return {"score": score, "groups": {"de_recovery": score / 100, "direction": 0.5, "distribution": 0.5,
                                       "covariance": 0.5},
            "skills": {}, "metrics": {}}


def _lib(gains, per_seed=None, seeds=range(10), guard_holds=True, commit="c" * 40, clean=True, perm_equal=False):
    """A task file: mix 60; oracle_shift = mix + gains['os']; proj_shift = mix + gains['ps']; proj_perm = mix +
    gains['pp']; entry_step = mix + gains['es']. per_seed: key -> ten values replacing a gain."""

    def g(k, i):
        return per_seed[k][i] if per_seed and k in per_seed else gains.get(k, 0.0)

    rows = []
    for i, s in enumerate(seeds):
        ps = g("ps", i)
        arms = {"copy_last": _arm(50.0), "oracle_mix": _arm(60.0), "oracle_shift": _arm(60.0 + g("os", i)),
                "proj_shift": _arm(60.0 + ps), "proj_perm": _arm(60.0 + (ps if perm_equal else g("pp", i))),
                "entry_step": _arm(60.0 + g("es", i))}
        rows.append({"seed": s, "arms": arms, "inverted_anchors": []})
    return {"label": "synthetic", "row": "T60", "library": "x", "seeds": list(seeds),
            "guard": {"holds": guard_holds}, "rows": rows,
            "projection": {"cm": {"group": "CM", "coefficient": 0.8, "coefficient_perm": 0.1, "explained_share": 0.3,
                                  "explained_share_perm": 0.01, "size2_on_moved": 4.0, "size2_off_moved": 1.0},
                           "blood": {"group": None, "coefficient": None, "coefficient_perm": None}},
            "truth": {"n_rows": 4321, "n_half_b": 4322},
            "stamp": {"commit": commit, "tree_clean": clean, "utc": "x", "inputs": {}}}


def _write(tmp_path, libs):
    out = tmp_path / "out"
    (out / "tasks").mkdir(parents=True)
    for lib, rec in libs.items():
        (out / "tasks" / f"lib{lib}.json").write_text(json.dumps(rec))
    return out


ROOM = {"os": 4.0, "ps": 1.0, "pp": 0.2, "es": 0.6}  # A 3.0, C 0.8, B 0.4
NONE = {"os": 3.0, "ps": 2.5, "pp": 0.5, "es": 1.0}  # A 0.5, C 2.0


def test_the_statements_follow_their_fixed_bars(rr, tmp_path):
    res = rr.read(_write(tmp_path, {"1": _lib(NONE), "2": _lib(NONE), "3": _lib(ROOM)}))
    r = res["reading"]
    assert r["valid"] and r["statements"] == {rr.STATEMENT_1: True, rr.STATEMENT_2: True}
    c = r["libraries"]["3"]["comparisons"]
    assert c["A"]["mean"] == pytest.approx(3.0) and c["C"]["mean"] == pytest.approx(0.8)
    assert c["B"]["mean"] == pytest.approx(0.4) and c["A"]["ahead"] == 10
    assert r["libraries"]["3"]["reached_share"] == pytest.approx(0.25)
    assert r["libraries"]["1"]["statements"] == {rr.STATEMENT_1: False, rr.STATEMENT_2: True}
    # the reading is library 3's: libraries 1 and 2 carry no statement into the outcome
    res2 = rr.t60_reading({"1": _lib(ROOM), "2": _lib(ROOM), "3": _lib(NONE)})
    assert res2["statements"] == {rr.STATEMENT_1: False, rr.STATEMENT_2: True}


def test_the_bars_are_strict_and_need_eight_seeds_ahead(rr):
    at_bar = rr.t60_reading({"3": _lib({"os": 1.5, "ps": 0.0, "pp": -0.3})})
    assert at_bar["valid"]
    assert at_bar["statements"] == {rr.STATEMENT_1: False, rr.STATEMENT_2: False}  # equal to a bar is not above it
    seven = [3.0] * 7 + [-0.5] * 3  # mean 1.95 above the bar, ahead on 7 of 10
    r = rr.t60_reading({"3": _lib({"ps": 0.0, "pp": -0.1}, per_seed={"os": seven})})
    assert r["libraries"]["3"]["comparisons"]["A"]["mean"] > rr.A_BAR
    assert r["statements"][rr.STATEMENT_1] is False
    eight = [2.5] * 8 + [-0.5] * 2
    r = rr.t60_reading({"3": _lib({"ps": 0.0, "pp": -0.1}, per_seed={"os": eight})})
    assert r["statements"][rr.STATEMENT_1] is True


def _noise_rates(rr, sd, draws, seed=20261119):
    rng = np.random.default_rng(seed)
    made = {rr.STATEMENT_1: 0, rr.STATEMENT_2: 0}
    for _ in range(draws):
        per = {k: rng.normal(0.0, sd, 10) for k in ("os", "ps", "pp", "es")}
        st = rr.t60_reading({"3": _lib({}, per_seed=per)})["statements"]
        for k in made:
            made[k] += st[k]
    return {k: v / draws for k, v in made.items()}


def test_on_pure_noise_the_statements_are_made_at_these_rates(rr):
    """Every arm the mix plus independent noise per seed, no effect anywhere. At 0.25 points per arm (a paired sd of
    0.35, above the per-seed sd of a paired difference on the committed library-3 files: 0.15 to 0.23 in T46, 0.18 in
    T28) statement (1) is never made in 2,000 draws and statement (2) in under 1% (0.5% here). At 1 point per arm (a
    paired sd of 1.4) statement (1) is made once in 2,000 draws, since its bar is 1.5 points; statement (2), whose bar
    is 0.3, is made in about 5% of draws: 8 of 10 seeds ahead happens by chance 5.5% of the time, and a paired sd of
    1.4 points puts the mean above 0.3 in most of those. The row states both."""
    low = _noise_rates(rr, 0.25, 2000)
    assert low[rr.STATEMENT_1] == 0.0 and low[rr.STATEMENT_2] < 0.01
    high = _noise_rates(rr, 1.0, 2000)
    assert high[rr.STATEMENT_1] < 0.002
    assert 0.02 < high[rr.STATEMENT_2] < 0.07


def test_the_void_conditions(rr):
    ok = {"1": _lib(NONE), "2": _lib(NONE), "3": _lib(ROOM)}
    assert rr.t60_reading(ok)["valid"]
    for bad, why in ((dict(guard_holds=False), "guard"), (dict(clean=False), "not clean"),
                     (dict(perm_equal=True), "proj_perm scores exactly as proj_shift"),
                     (dict(seeds=range(9)), "seeds")):
        r = rr.t60_reading({**ok, "3": _lib(ROOM, **bad)})
        assert r["outcome"] == "VOID" and any(why in x for x in r["reasons"]), why
    r = rr.t60_reading({**ok, "3": _lib({"os": 0.0, "ps": 1.0, "pp": 0.2})})
    assert r["outcome"] == "VOID" and any("oracle_shift scores exactly as oracle_mix" in x for x in r["reasons"])
    r = rr.t60_reading({"1": None, "2": _lib(NONE), "3": None})
    assert r["outcome"] == "VOID" and "no task file" in r["reasons"][0]
    # a void portion library is named and reported; the reading on library 3 stands
    r = rr.t60_reading({**ok, "1": _lib(NONE, guard_holds=False)})
    assert r["valid"] and not r["libraries"]["1"]["valid"] and r["libraries"]["2"]["valid"]
    r = rr.t60_reading({**ok, "1": None})
    assert r["valid"] and r["libraries"]["1"]["reasons"] == ["no task file"]
    # the valid libraries at two commits: void
    r = rr.t60_reading({**ok, "2": _lib(NONE, commit="d" * 40)})
    assert r["outcome"] == "VOID" and any("2 commits" in x for x in r["reasons"])
    missing = _lib(ROOM)
    del missing["rows"][3]["arms"]["proj_perm"]
    r = rr.t60_reading({**ok, "3": missing})
    assert r["outcome"] == "VOID" and any("missing" in x for x in r["reasons"])
    inv = _lib(ROOM)
    inv["rows"][0]["inverted_anchors"] = ["mmd_u"]
    assert rr.t60_reading({**ok, "3": inv})["outcome"] == "VOID"


def test_the_text_holds_the_statements_and_limits_and_no_count(rr, tmp_path):
    res = rr.read(_write(tmp_path, {"1": _lib(NONE), "2": _lib(NONE), "3": _lib(ROOM)}))
    txt = (tmp_path / "out" / "t60.txt").read_text()
    assert f"{rr.STATEMENT_1}: true" in txt and f"{rr.STATEMENT_2}: true" in txt
    for line in rr.LIMITS:
        assert line in txt
    body = txt.rsplit("\nCommit ", 1)[0]
    assert not re.search(r"(?<![\d.])(4321|4322)(?![\d.])", body)  # the truth's row counts are never printed
    assert json.loads((tmp_path / "out" / "t60.json").read_text())["reading"]["valid"] == res["reading"]["valid"]
    with pytest.raises(SystemExit):  # it refuses to write over its outputs
        rr.read(tmp_path / "out")


def test_the_run_is_refused_without_its_row(rr, monkeypatch):
    monkeypatch.setattr(rr.n12, "tries_rows", lambda: {"T59"})
    with pytest.raises(SystemExit, match="T60 has no row"):
        rr.main(["run", "--library", "3", "--out", "x", "--data", "/nonexistent"])


def test_the_constants_are_the_rows(rr):
    assert rr.COUNTED_ARMS == ("oracle_shift", "proj_shift", "proj_perm", "entry_step")
    assert rr.REFERENCE_ARMS == ("copy_last", "oracle_mix")
    assert (rr.A_BAR, rr.C_BAR, rr.MIN_AHEAD, rr.READ_LIBRARY) == (1.5, 0.3, 8, "3")
    assert (rr.N_CELLS, rr.DAMPING, rr.EPS, rr.MIN_CELLS, rr.PERM_BINS) == (5000, 1.0, 0.5, 20, 10)
    assert rr.READINGS == {"A": ("oracle_shift", "proj_shift"), "C": ("proj_shift", "proj_perm"),
                           "B": ("proj_shift", "entry_step")}
