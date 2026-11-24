"""Row E55 (as specified by the cold challenger, round_1124_t2): the composed shape form overlay_stretch and its two
harness-only placebos (vec.shape_interp), its builder option (`build_test.py embryo --shape-interp overlay_stretch`),
the read script experiments/coords/e55_overlay_stretch_read.py and its launcher experiments/coords/run_e55.sh. Synthetic
clouds, synthetic stages and hand-made run and gate files only: nothing is read from data/ and nothing is scored.

Pinned: the stretched union's normalised SDs equal the target along the blocks; the size rule's RMS is unchanged; the
pool is overlay's pool under the same affine map with overlay's weights; the anti placebo takes the anti target, the
axes placebo the arm's own f along V @ Q with Q from default_rng([seed, AXES_STREAM]); determinism per seed; the
existing tuples and forms unchanged; the builder's refusals; each clause P1 to P7 and each void condition of the read;
the author's five rulings (the O guard per task against E49, the gate record and its 1e-9 recompute, P7's merge, a
placebo refusal is a STOP, the count reading changes nothing); the row check; the launcher's syntax and refusals."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import harness
from vec import shape_interp as si
from vec import t2_build as tb

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "coords" / "run_e55.sh"
SEEDS = list(range(10))
COMMIT = "a" * 40
ARM, ANTI, AXES = "fgw_shape0x_ov_sr", "fgw_shape0x_ov_sr_anti", "fgw_shape0x_ov_sr_axes"
O = "fgw_shape0x_ov"  # noqa: E741  (the specification names the control O)
C0 = "fgw_shape0x"
METRICS = ("de_score", "de_direction", "mmd_u", "variogram", "d2_shape", "occupancy_dice", "scale_log_ratio",
           "neighborhood_mmd")


@pytest.fixture(scope="module")
def e55():
    spec = importlib.util.spec_from_file_location("e55_overlay_stretch_read",
                                                  REPO / "experiments" / "coords" / "e55_overlay_stretch_read.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _rot(seed: int) -> np.ndarray:
    from scipy.spatial.transform import Rotation

    return Rotation.random(random_state=np.random.default_rng(seed)).as_matrix()


def _ellipsoid(n: int, sds, seed: int, R=None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 3))
    X -= X.mean(0)
    w, V = np.linalg.eigh(X.T @ X / n)
    X = X @ V @ np.diag(w ** -0.5) @ V.T
    X = X * np.asarray(sds, float)
    return X @ (np.eye(3) if R is None else R)


def _frame(C):
    return (C - C.mean(0)) / tb.rms_radius(C)


@pytest.fixture(scope="module")
def pair():
    """Two brackets whose axis 1 is well separated from a round 2-3 plane (blocks 1|2,3), in the earlier frame."""
    Zb = _frame(_ellipsoid(2000, [3.0, 1.0, 0.97], 1, _rot(5)))
    Za = _frame(_ellipsoid(1600, [3.6, 1.0, 0.98], 2, _rot(6)))
    rows_b, rows_a = np.arange(0, 300), np.arange(100, 400)
    perm = np.random.default_rng(3).permutation(600)
    return Zb, Za, rows_b, rows_a, perm, _rot(7)


def _form(pair, form, seed=1_000_001, t=0.5):
    Zb, Za, rows_b, rows_a, perm, R = pair
    return si.apply_form(form, Z=None, Zb=Zb, Za=Za, t=t, R=R, seed=seed, rows_b=rows_b, rows_a=rows_a, perm=perm,
                         t_mix=0.5)


# ------------------------------------------------------------------------------------------------ the form's algebra
def test_the_new_names_live_in_their_own_tuples_and_the_old_ones_are_unchanged():
    assert si.FORMS == ("stretch_ratio", "stretch_registered", "overlay")
    assert si.CONTROLS == ("stretch_ratio_anti", "stretch_registered_anti", "stretch_registered_rot", "overlay_rot")
    assert si.ALL_FORMS == (*si.FORMS, *si.CONTROLS)
    assert si.COMPOSED_FORMS == ("overlay_stretch",)
    assert si.COMPOSED_CONTROLS == ("overlay_stretch_anti", "overlay_stretch_axes")
    assert not set(si.COMPOSED_FORMS + si.COMPOSED_CONTROLS) & set(si.ALL_FORMS + si.GEODESIC_FORMS)
    assert harness.OSR_CANDIDATES == {ARM: "overlay_stretch", ANTI: "overlay_stretch_anti",
                                      AXES: "overlay_stretch_axes"}
    for name, form in harness.OSR_CANDIDATES.items():
        c = harness.CANDIDATES[name]
        assert c.params["shape_form"] == form and c.params["positions"] == tb.POSITIONS[1]
        assert c.requires is harness.CANDIDATES[C0].requires


def test_the_axes_stream_is_used_by_no_build_step():
    """AXES_STREAM appears in no other module's default_rng call or stream constant under src/vec."""
    assert si.AXES_STREAM == 56 and si.AXES_STREAM not in (si.IDENT_STREAM, si.ROT_STREAM)
    for p in (REPO / "src" / "vec").glob("*.py"):
        text = p.read_text()
        if p.name == "shape_interp.py":
            continue
        assert not re.search(r"STREAMS?\s*=\s*\(?[^\n#]*\b56\b", text), p
        assert not re.search(r"default_rng\(\[[^\]]*\b56\]", text), p


def test_the_stretched_union_reaches_the_target_along_the_blocks(pair):
    Zb, Za, rows_b, rows_a, perm, R = pair
    out = _form(pair, "overlay_stretch")
    U = si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.5).Z
    blocks = si.bracket_blocks(Zb, Za)
    assert si.blocks_text(blocks) == "1|2,3" == out.record["blocks"]
    target = si.log_interp_sds(si.normalised_sds(Zb), si.normalised_sds(Za), 0.5)
    s_o, V = si.principal(U)
    f = si.stretch_factors(s_o, target, blocks)
    M = si.stretch_matrix(V, f)
    np.testing.assert_allclose(out.Z, U.mean(0) + (U - U.mean(0)) @ M, atol=1e-13)
    # along the blocks: axis 1 at the target, the 2-3 plane's sum of squares at the target's
    sd = np.sqrt(((out.Z - out.Z.mean(0)) @ V.T).var(0))
    sd = sd / np.sqrt((sd ** 2).sum())
    assert sd[0] == pytest.approx(target[0], rel=1e-10)
    assert sd[1] ** 2 + sd[2] ** 2 == pytest.approx(target[1] ** 2 + target[2] ** 2, rel=1e-10)
    assert sd[1] / sd[2] == pytest.approx(s_o[1] / s_o[2], rel=1e-10)  # the plane keeps its own shape
    rec = out.record
    assert [rec["f1"], rec["f2"], rec["f3"]] == pytest.approx(list(f))
    assert rec["g"] == pytest.approx(np.log(f[0] / f[1]))
    assert [rec[f"target_sd{j}"] for j in (1, 2, 3)] == pytest.approx(list(target))
    assert rec["thin_sd"] == pytest.approx(si.normalised_sds(out.Z)[-1])
    assert rec["form"] == "overlay_stretch" and rec["ident_stream"] == si.IDENT_STREAM


def test_with_three_identifiable_axes_the_normalised_sds_are_the_target():
    Zb = _frame(_ellipsoid(2000, [3.0, 1.8, 0.7], 11, _rot(12)))
    Za = _frame(_ellipsoid(2000, [3.5, 1.6, 0.6], 13, _rot(14)))
    rows_b, rows_a = np.arange(400), np.arange(400)
    perm = np.random.default_rng(0).permutation(800)
    out = si.overlay_stretch(Zb, Za, rows_b, rows_a, perm, np.eye(3), 0.5, t=0.5)
    assert out.record["blocks"] == "1|2|3"
    target = si.log_interp_sds(si.normalised_sds(Zb), si.normalised_sds(Za), 0.5)
    np.testing.assert_allclose(si.normalised_sds(out.Z), target, rtol=1e-10)


def test_the_pool_is_overlays_pool_under_the_same_map_with_overlays_weights(pair):
    Zb, Za, rows_b, rows_a, perm, R = pair
    ov = si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.5)
    for form, _variant in (("overlay_stretch", "arm"), ("overlay_stretch_anti", "anti"),
                          ("overlay_stretch_axes", "axes")):
        out = _form(pair, form)
        U, m = ov.Z, ov.Z.mean(0)
        M = np.linalg.lstsq(U - m, out.Z - m, rcond=None)[0]
        np.testing.assert_allclose(out.pool, m + (ov.pool - m) @ M, atol=1e-10)
        np.testing.assert_array_equal(out.pool_weights, ov.pool_weights)
        assert np.linalg.det(M) > 0 and np.allclose(M, M.T, atol=1e-12)  # symmetric, no mirror image


def test_the_anti_placebo_takes_the_anti_target_and_nothing_else(pair):
    Zb, Za, rows_b, rows_a, perm, R = pair
    arm, anti = _form(pair, "overlay_stretch"), _form(pair, "overlay_stretch_anti")
    U = si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.5).Z
    blocks = si.bracket_blocks(Zb, Za)
    s_o, V = si.principal(U)
    t_anti = si.log_interp_sds(si.normalised_sds(Zb), si.normalised_sds(Za), 0.5, anti=True)
    f = si.stretch_factors(s_o, t_anti, blocks)
    np.testing.assert_allclose(anti.Z, U.mean(0) + (U - U.mean(0)) @ si.stretch_matrix(V, f), atol=1e-13)
    assert [anti.record[f"target_sd{j}"] for j in (1, 2, 3)] == pytest.approx(list(t_anti))
    assert np.sign(anti.record["g"]) == -np.sign(arm.record["g"])  # the step mirrored
    assert "axes_stream" not in anti.record


def test_the_axes_placebo_takes_the_arms_own_f_along_a_random_frame(pair):
    Zb, Za, rows_b, rows_a, perm, R = pair
    seed = 1_000_004
    arm, axes = _form(pair, "overlay_stretch", seed), _form(pair, "overlay_stretch_axes", seed)
    U = si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.5).Z
    _, V = si.principal(U)
    f = np.array([arm.record[k] for k in ("f1", "f2", "f3")])
    assert [axes.record[k] for k in ("f1", "f2", "f3")] == pytest.approx(list(f))
    from scipy.spatial.transform import Rotation

    Q = Rotation.random(random_state=np.random.default_rng([seed, 56])).as_matrix()
    np.testing.assert_allclose(si.axes_rotation(seed), Q)
    M_ax = (V @ Q).T @ np.diag(f) @ (V @ Q)
    np.testing.assert_allclose(axes.Z, U.mean(0) + (U - U.mean(0)) @ M_ax, atol=1e-13)
    assert np.allclose(np.sort(np.linalg.eigvalsh(M_ax)), np.sort(f))  # the same magnitude of distortion
    assert not np.allclose(axes.Z, arm.Z)
    assert axes.record["axes_stream"] == si.AXES_STREAM == 56
    assert not np.allclose(si.axes_rotation(seed), si.random_rotation(seed))  # not ROT_STREAM's Q


def test_each_form_is_deterministic_per_seed(pair):
    for form in (*si.COMPOSED_FORMS, *si.COMPOSED_CONTROLS):
        a, b = _form(pair, form, 1_000_002), _form(pair, form, 1_000_002)
        np.testing.assert_array_equal(a.Z, b.Z)
        np.testing.assert_array_equal(a.pool, b.pool)
        assert a.record == b.record
    assert not np.allclose(_form(pair, "overlay_stretch_axes", 1_000_002).Z,
                           _form(pair, "overlay_stretch_axes", 1_000_003).Z)
    np.testing.assert_array_equal(_form(pair, "overlay_stretch", 1_000_002).Z,
                                  _form(pair, "overlay_stretch", 1_000_003).Z)  # the arm uses no seed


def test_the_existing_forms_are_unchanged_by_the_dispatch(pair):
    """Every existing form through apply_form gives exactly its own function's output (no new branch catches it)."""
    Zb, Za, rows_b, rows_a, perm, R = pair
    Z = si.overlay(Zb, Za, rows_b, rows_a, perm, np.eye(3), 0.5).Z
    seed = 1_000_001
    want = {"stretch_ratio": si.stretch_ratio(Z, Zb, Za, 0.5),
            "stretch_ratio_anti": si.stretch_ratio(Z, Zb, Za, 0.5, anti=True),
            "stretch_registered_anti": si.stretch_registered(Z, Zb, Za, 0.5, R, variant="anti"),
            "stretch_registered_rot": si.stretch_registered(Z, Zb, Za, 0.5, R, variant="rot",
                                                            Q=si.random_rotation(seed)),
            "overlay": si.overlay(Zb, Za, rows_b, rows_a, perm, R, 0.5),
            "overlay_rot": si.overlay(Zb, Za, rows_b, rows_a, perm, si.random_rotation(seed), 0.5)}
    for form, w in want.items():
        got = si.apply_form(form, Z=Z, Zb=Zb, Za=Za, t=0.5, R=R, seed=seed, rows_b=rows_b, rows_a=rows_a, perm=perm,
                            t_mix=0.5)
        np.testing.assert_array_equal(got.Z, w.Z)
        np.testing.assert_array_equal(got.pool, w.pool)
        assert {k: v for k, v in got.record.items() if k not in ("form", "rot_stream", "ident_stream")} == w.record


# ------------------------------------------------------------------------------------------------ the hook
@pytest.fixture(scope="module")
def stages():
    from test_occupancy_frame import _stage

    return _stage(1500, 1, 0.0), _stage(900, 2, 0.3)


def _interp(stages, **kw):
    b, a = stages
    return tb.interpolate_clouds(np.asarray(b.obsm["spatial_3D"], float), np.asarray(a.obsm["spatial_3D"], float),
                                 b.X.toarray(), a.X.toarray(), 0.5, 400, 1_000_001, target_rms=6.0, t_pos=0.0,
                                 positions="floor_exact", n_anchor=150, **kw)


@pytest.mark.parametrize("form", [*si.COMPOSED_FORMS, *si.COMPOSED_CONTROLS])
def test_the_hook_keeps_rows_order_count_and_the_size_rule(stages, form):
    ov, r = _interp(stages, shape_form="overlay"), _interp(stages, shape_form=form)
    np.testing.assert_array_equal(r.rows, ov.rows)
    np.testing.assert_array_equal(r.stage, ov.stage)
    assert tb.rms_radius(r.coords) == pytest.approx(6.0, rel=1e-13)
    assert tb.rms_radius(r.coords) == pytest.approx(tb.rms_radius(ov.coords), rel=1e-13)
    np.testing.assert_allclose(r.coords.mean(0), ov.coords.mean(0), atol=1e-12)
    assert not np.allclose(r.coords, ov.coords)
    rec = r.record
    assert rec["shape_interp"] == form and rec["floor_exact"] == ov.record["floor_exact"]
    for k in ("f1", "f2", "f3", "g", "target_sd1", "target_sd2", "target_sd3", "thin_sd"):
        assert isinstance(rec[f"shape_interp_{k}"], float)
    assert r.pool_weights.sum() == pytest.approx(1.0)
    np.testing.assert_array_equal(r.pool_weights, ov.pool_weights)


# ------------------------------------------------------------------------------------------------ the builder
@pytest.fixture(scope="module")
def bt():
    from test_occupancy_frame import REPO as OF_REPO

    sys.path.insert(0, str(OF_REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def embryo_env(bt, monkeypatch, own_repo, tmp_path):
    from test_occupancy_frame import install

    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path, reads)
    return reads


def test_the_builder_offers_overlay_stretch_off_by_default_and_refuses_the_placebos_by_name(bt):
    from test_occupancy_frame import ARMS

    fe = ARMS["floor_exact"]
    assert bt.build_parser().parse_args(["embryo", *fe, "--out", "x.h5ad"]).shape_interp == "none"
    a = bt.build_parser().parse_args(["embryo", *fe, "--shape-interp", "overlay_stretch", "--out", "x.h5ad"])
    assert a.shape_interp == "overlay_stretch"
    for placebo in si.COMPOSED_CONTROLS:
        with pytest.raises(SystemExit):
            bt.build_parser().parse_args(["embryo", *fe, "--shape-interp", placebo, "--out", "x.h5ad"])


def test_the_builder_refuses_a_placebo_by_name_even_past_the_parser(bt, embryo_env, tmp_path, monkeypatch):
    from test_occupancy_frame import ARMS

    a = bt.build_parser().parse_args(["embryo", *ARMS["floor_exact"], "--shape-interp", "overlay_stretch",
                                      "--shape-interp-test-ok", "--out", str(tmp_path / "p" / "T2embryo.h5ad")])
    for placebo in (*si.COMPOSED_CONTROLS, "overlay_rot"):
        a.shape_interp = placebo
        with pytest.raises(SystemExit, match="never built"):
            a.func(a) if hasattr(a, "func") else bt.build_embryo(a)
    assert not (tmp_path / "p").exists() and embryo_env == []


def test_overlay_stretch_on_a_test_board_needs_the_named_flag(bt, embryo_env, tmp_path, capsys):
    from test_occupancy_frame import ARMS, SECRET_MAX, SECRET_TRUTH, _shows

    out = tmp_path / "g" / "T2embryo.h5ad"
    with pytest.raises(SystemExit, match="needs the captain's yes") as e:
        bt.main(["embryo", *ARMS["floor_exact"], "--shape-interp", "overlay_stretch", "--out", str(out)])
    assert not out.exists() and embryo_env == []
    assert not _shows(str(e.value) + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)
    with pytest.raises(SystemExit, match="floor_exact only"):
        bt.main(["embryo", *ARMS["moved"], "--shape-interp", "overlay_stretch", "--shape-interp-test-ok",
                 "--out", str(out)])


def test_overlay_stretch_builds_positions_only_and_records_its_block(bt, embryo_env, tmp_path, capsys):
    from test_occupancy_frame import ARMS, SECRET_MAX, SECRET_TRUTH, _build, _file, _shows

    fe = ARMS["floor_exact"]
    base, _, _ = _build(bt, tmp_path, "base", *fe, "--shape-interp", "overlay", "--shape-interp-test-ok")
    out, rec, text = _build(bt, tmp_path, "osr", *fe, "--shape-interp", "overlay_stretch", "--shape-interp-test-ok")
    X0, C0_ = _file(base)
    X1, C1 = _file(out)
    np.testing.assert_array_equal(X1, X0)
    assert len(C1) == len(C0_) and not np.allclose(C1, C0_)
    assert tb.rms_radius(C1.astype(float)) == pytest.approx(tb.rms_radius(C0_.astype(float)), rel=1e-6)
    blk = rec["shape_interp"]
    assert blk["form"] == "overlay_stretch" and blk["row"] == "E55 (as specified)" and blk["ruling"] is None
    for k in ("f1", "f2", "f3", "g", "target_sd1", "target_sd2", "target_sd3", "thin_sd"):
        assert k in blk
    assert all(not isinstance(v, (dict, list)) for v in blk.values())
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)


# ------------------------------------------------------------------------------------------------ the read: constants
def test_the_constants_and_the_rule_text(e55):
    assert e55.ARM == ARM and e55.PLACEBOS == (ANTI, AXES) and e55.O == O and e55.C0 == C0
    assert e55.TRY_IDS == {ARM: "E55a", ANTI: "E55b", AXES: "E55c"}
    assert e55.BUILDS == ("main", "prefix") and e55.TEST_LIKE == ("half", "late", "new1", "new2")
    assert (e55.BAR, e55.AHEAD, e55.G_MIN, e55.EXPR_TOL, e55.SCALE_TOL, e55.GATE_TOL, e55.O_TOL, e55.G_TOL) == \
        (0.3, 8, 0.05, 1e-9, 1e-6, 0.01, 0.01, 1e-9)
    assert e55.GATE_COMMIT == "d145d7882e4113b450c643a5c23ca225f7c57ff0"
    assert set(e55.CLAUSES) == {"P1", "P1b", "P2", "P3", "P4", "P5", "P6", "P7"}
    for k, v in e55.CLAUSES.items():
        assert v in e55.PASS_RULE, k  # word for word from the specification's pass rule
    assert "selected on" in e55.ROW_FIRST_LINE and e55.LABEL == "diagnostic, selected on the nine tasks"
    assert "--shape-interp overlay_stretch --shape-interp-test-ok" in e55.CANDIDATE_COMMAND
    assert list(e55.SEEDS) == SEEDS
    assert e55.GATE_SHA256 == "7c0fc6b6398092d54f6e4318bd8c79c1e092c4ec7e1edbbb6b86fcbede7b4d67"
    assert e55.GATE_RECORD == REPO / "experiments" / "coords" / "results" / "osr_gate_g1"


def test_the_row_check(e55, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    for text, want in (("| E54 (one try) | embryo |\n", False), ("| E55a | x |\n", False), ("E55 in prose\n", False),
                       ("| E55 | x |\n", True), ("| E55 (three tries) | embryo |\n", True)):
        tries.write_text(text)
        assert e55.has_row(tries) is want, text
    assert not e55.has_row(tmp_path / "missing.md")
    tries.write_text("| E54 | x |\n")
    monkeypatch.setattr(e55, "TRIES", tries)
    with pytest.raises(SystemExit, match="E55 has no row"):
        e55._refuse_unless_ready(COMMIT)
    with pytest.raises(SystemExit, match="E55 has no row"):
        e55.main(["read", "--commit", COMMIT, "--runs", str(tmp_path), "--out", str(tmp_path / "o"),
                  "--gate-record", str(tmp_path), "--gate-sha256", "0" * 64])


def test_an_unclean_tree_is_refused(e55, tmp_path, monkeypatch):
    tries = tmp_path / "TRIES.md"
    tries.write_text("| E55 | x |\n")
    monkeypatch.setattr(e55, "TRIES", tries)
    monkeypatch.setattr(e55, "repo_state", lambda repo: type("S", (), {"unclean": ("x.py",), "commit": COMMIT})())
    with pytest.raises(SystemExit, match="not clean"):
        e55._refuse_unless_ready(COMMIT)


# ------------------------------------------------------------------------------------------------ the arms' wrappers
def _ctx():
    from types import SimpleNamespace

    return SimpleNamespace(seed=harness.BUILD_SEED_OFFSET + 3)


def _stage(k: float):
    from vec.evaluate import Stage

    return Stage(X=np.full((50, 2), k, np.float32), coords=np.full((50, 3), k), genes=["a", "b"])


def test_a_refused_placebo_is_recorded_and_os_file_stands_in(e55):
    from vec.boards import Board

    board = Board(key="pseudo:T2e/x", task="T2", genes=["a"], needs_coords=True, split="pseudo", min_cells=1,
                  max_cells=50, truth_cells=10)
    o = harness.Candidate(O, frozenset({"T2"}), (), lambda ctx: _stage(1), {})

    def refuse(ctx):
        raise si.ShapeRefused("out of range")

    def guard(ctx):
        raise ValueError("pseudo: cloud guard flags ['coverage']")

    with e55.captured() as seen:
        seen["board"] = board
        ref = e55.observed(o, seen).build(_ctx())
        for name, fn in ((ANTI, refuse), (AXES, guard)):
            out = e55.guarded(harness.Candidate(name, frozenset({"T2"}), (), fn, {}), seen).build(_ctx())
            np.testing.assert_array_equal(out.coords, ref.coords)
        with pytest.raises(ValueError, match="else"):
            e55.guarded(harness.Candidate(ARM, frozenset({"T2"}), (), lambda c: (_ for _ in ()).throw(
                ValueError("something else")), {}), seen).build(_ctx())
    assert seen["refused"][(ANTI, 3)].startswith("shape:") and seen["refused"][(AXES, 3)] == "cloud guard"
    test = Board(**{**board.__dict__, "key": "T2:embryo:test_interp", "split": "test_interp"})
    with e55.captured() as seen:
        seen["board"] = test
        with pytest.raises(RuntimeError, match="test board"):
            e55.guarded(o, seen).build(_ctx())
        with pytest.raises(RuntimeError, match="test board"):
            e55.observed(o, seen).build(_ctx())


# ------------------------------------------------------------------------------------------- the O guard (ruling 2)
def _e49_reading(tmp_path, per_task, commit="d76ccd5" + "0" * 33, void=False):
    p = tmp_path / "e49.json"
    p.write_text(json.dumps({"void": void, "commit": commit,
                             "tasks": {n: {O: {"board": {"mean": v}}} for n, v in per_task.items()}}))
    return p


def test_the_o_guard_is_per_task_against_e49_only(e55, tmp_path):
    names = list(e55.nt.all_tasks())
    ref = dict.fromkeys(names, 2.0)
    o_ref = e55.e49_o_per_task(_e49_reading(tmp_path, ref))
    assert o_ref[1] == [] and o_ref[0] == ref
    key = e55.nt.all_tasks()["e75"].key
    old, why = e55.e47_scores(key)
    assert why == []
    scores = {c: np.asarray(old[c], float).copy() for c in e55.GATE_ARMS}
    scores[O] = scores[C0] + 2.0
    g = e55.gate_check("e75", key, scores, o_ref)
    assert g["pass"] and g["o_off_e49"] == pytest.approx(0.0, abs=1e-12)
    wobble = scores[O].copy()
    wobble[3] += 0.05  # per seed O moves; its ten-seed mean moves by 0.005: no per-seed check of O
    wobble[4] -= 0.0
    g = e55.gate_check("e75", key, {**scores, O: wobble}, o_ref)
    assert g["pass"]
    g = e55.gate_check("e75", key, {**scores, O: scores[O] + 0.011}, o_ref)
    assert not g["pass"] and any("off E49's ten-seed mean" in w for w in g["why"])
    moved = {**scores, C0: scores[C0].copy()}
    moved[C0][4] += 0.011  # C0 per seed against E47, as E49's gate
    moved[O] = moved[C0] + 2.0
    assert not e55.gate_check("e75", key, moved, o_ref)["pass"]
    void = e55.e49_o_per_task(_e49_reading(tmp_path, ref, void=True))
    assert void[1] and not e55.gate_check("e75", key, scores, void)["pass"]


def test_the_committed_e49_reading_serves_as_the_reference(e55):
    ref, why = e55.e49_o_per_task()
    assert why == [] and set(ref) == set(e55.nt.all_tasks())


def test_gate_files_are_required_at_the_run_commit(e55, tmp_path):
    assert len(e55.gate_problems(tmp_path, COMMIT)) == 9
    for name in e55.nt.all_tasks():
        p = e55._gate_path(tmp_path, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"pass": True, "commit": COMMIT}))
    assert e55.gate_problems(tmp_path, COMMIT) == []
    assert len(e55.gate_problems(tmp_path, "c" * 40)) == 9


# ------------------------------------------------------------------------------------------- the gate record (ruling 4)
GATE_MEDIANS = {"test": 0.0383, "half": 0.0313, "late": 0.0740, "new1": 0.0120, "new2": 0.0460, "e75": 0.0078,
                "e725": 0.0144, "x1": 0.0977, "x2": 0.1037, "x3": 0.1618}  # hand-made, the shape of a record


def _gate_folder(e55, folder: Path, commit: str, shift: float = 0.0, verdict="PASS") -> str:
    folder.mkdir(parents=True, exist_ok=True)
    seeds = [str(s) for s in e55.og.SEEDS]
    for name, med in GATE_MEDIANS.items():
        per = {s: {"g": med + 0.001 * (i - 4.5) + shift, "blocks": "1|2,3"} for i, s in enumerate(seeds)}
        (folder / f"{name}.json").write_text(json.dumps({"name": name, "commit": commit, "per_seed": per,
                                                         "median_g": med + shift}))
    gate = {"verdict": verdict, "stopped_by": [] if verdict == "PASS" else ["S1"],
            "median_g": {n: m + shift for n, m in GATE_MEDIANS.items()}, "commits": [commit]}
    (folder / "gate.json").write_text(json.dumps(gate))
    return hashlib.sha256((folder / "gate.json").read_bytes()).hexdigest()


@pytest.fixture
def gate_dirs(e55, tmp_path):
    rec, new = tmp_path / "record", tmp_path / "runs" / "gate0"
    sha = _gate_folder(e55, rec, e55.GATE_COMMIT)
    _gate_folder(e55, new, COMMIT)
    return rec, sha, new


def test_the_gate_record_and_its_recompute_pass(e55, gate_dirs):
    rec, sha, new = gate_dirs
    why, gate = e55.gate0_problems(rec, sha, new, COMMIT)
    assert why == [] and gate["median_g"]["late"] == pytest.approx(0.074) and len(gate["test_pair_g"]) == 10


@pytest.mark.parametrize("breach", ["missing", "sha", "record_commit", "recompute_commit", "g_off", "g_within",
                                    "stop_record", "stop_recompute", "recompute_missing_pair", "seeds"])
def test_the_gate_record_void_conditions(e55, gate_dirs, tmp_path, breach):
    rec, sha, new = gate_dirs
    args = [rec, sha, new, COMMIT]
    if breach == "missing":
        args[0] = tmp_path / "nowhere"
    elif breach == "sha":
        args[1] = "0" * 64
    elif breach == "record_commit":
        r = json.loads((rec / "x2.json").read_text())
        r["commit"] = "b" * 40
        (rec / "x2.json").write_text(json.dumps(r))
    elif breach == "recompute_commit":
        args[3] = "c" * 40
    elif breach in ("g_off", "g_within"):
        r = json.loads((new / "late.json").read_text())
        r["per_seed"]["1000005"]["g"] += 2e-9 if breach == "g_off" else 5e-10
        (new / "late.json").write_text(json.dumps(r))
    elif breach == "stop_record":
        args[1] = _gate_folder(e55, rec, e55.GATE_COMMIT, verdict="STOP")
    elif breach == "stop_recompute":
        _gate_folder(e55, new, COMMIT, verdict="STOP")
    elif breach == "recompute_missing_pair":
        (new / "x3.json").rename(new / "x3.moved")
    else:
        r = json.loads((new / "half.json").read_text())
        r["per_seed"].pop("1000003")
        (new / "half.json").write_text(json.dumps(r))
    why, _ = e55.gate0_problems(*args)
    if breach == "g_within":
        assert why == []
    else:
        assert why


def test_p7s_task_set_and_the_test_pair_statement(e55):
    groups4 = e55.merges()["four"]
    at, below = e55.p7_tasks({k: v for k, v in GATE_MEDIANS.items() if k != "test"}, groups4)
    assert at == ["late", "x1", "x2", "x3"] and below == ["e75", "half", "new1", "new2"]
    assert e55.restricted_groups(groups4, at) == {"T2/E8.0.h5ad": ["late"], "T2/E8.25_late.h5ad": ["x1"],
                                                  "T2/E8.5.h5ad": ["x2", "x3"]}
    st = e55.tp_statement(GATE_MEDIANS, [0.0154, 0.0582, 0.0135, 0.0174, 0.0522, 0.0426, 0.0339, 0.0562, 0.0208,
                                         0.1121])
    assert st["below_p7"] and st["seeds_below"] == 6 and st["nearest"] == ["half", "new2"]
    assert st["small_g_test_like"] == ["half", "new1", "new2"]
    assert "range from zero" in st["text"] and "0.0383" in st["text"]


# ------------------------------------------------------------------------------------------------ decide: P1 to P7
def _good():
    up = np.full(10, 0.5)
    return {"board4": up.copy(), "board5": up.copy(), "sg4": np.full(10, 0.01), "board_tl": up.copy(),
            "sg_tl": np.full(10, 0.01),
            "groups5": {g: np.zeros(10) for g in ("expression_change", "cell_state", "local_organisation")},
            "groups4": {g: np.zeros(10) for g in ("expression_change", "cell_state", "local_organisation")},
            "S_vs": {ANTI: np.full(10, 0.01), AXES: np.full(10, 0.01)}, "board4_g": np.full(10, 0.1),
            "arm_refused": 0, "placebo_refused": {ANTI: 0, AXES: 0}}


def test_decide_passes_when_all_eight_clauses_hold(e55):
    r = e55.decide(_good())
    assert r["passes"] and set(r["clauses"]) == set(e55.CLAUSES) and all(r["clauses"].values())


def _seven_ahead(mean_target=0.5):
    """Mean above the bar but only 7 of 10 seeds above 0."""
    D = np.full(10, 1.5)
    D[:3] = -0.5 - 0.0
    return D * (mean_target / D.mean())


@pytest.mark.parametrize("case, clause", [
    ("P1_bar", "P1"), ("P1_seeds", "P1"), ("P1b_bar", "P1b"), ("P1b_seeds", "P1b"), ("P2", "P2"),
    ("P3_board", "P3"), ("P3_shape", "P3"), ("P4_five", "P4"), ("P4_four", "P4"), ("P5_mean", "P5"),
    ("P5_seeds", "P5"), ("P5_axes", "P5"), ("P5_placebo_refused", "P5"), ("P6", "P6"), ("P7_none", "P7"),
    ("P7_below", "P7")])
def test_each_clause_can_stop_the_arm(e55, case, clause):
    d = _good()
    if case == "P1_bar":
        d["board4"] = np.full(10, 0.29)
    elif case == "P1_seeds":
        d["board4"] = _seven_ahead()
    elif case == "P1b_bar":
        d["board5"] = np.full(10, 0.3)  # "above +0.3": the bar itself fails
    elif case == "P1b_seeds":
        d["board5"] = _seven_ahead()
    elif case == "P2":
        d["sg4"] = np.full(10, -0.001)
    elif case == "P3_board":
        d["board_tl"] = np.zeros(10)
    elif case == "P3_shape":
        d["sg_tl"] = np.full(10, -0.01)
    elif case == "P4_five":
        D = np.full(10, -0.0006)
        D[:2] = 0.0001
        d["groups5"]["local_organisation"] = D  # 8 of 10 below 0, mean below 0: behind, no tolerance
    elif case == "P4_four":
        d["groups4"]["cell_state"] = np.full(10, -1e-6)
    elif case == "P5_mean":
        d["S_vs"][ANTI] = np.full(10, -0.01)
    elif case == "P5_seeds":
        D = np.full(10, 0.02)
        D[:3] = -0.001
        d["S_vs"][ANTI] = D
    elif case == "P5_axes":
        d["S_vs"][AXES] = np.zeros(10)
    elif case == "P5_placebo_refused":
        d["placebo_refused"][AXES] = 1
    elif case == "P6":
        d["arm_refused"] = 1
    elif case == "P7_none":
        d["board4_g"] = None
    else:
        d["board4_g"] = np.full(10, -0.01)
    r = e55.decide(d)
    assert not r["passes"] and not r["clauses"][clause]
    assert all(v for k, v in r["clauses"].items() if k != clause)


def test_p4_is_not_behind_with_seven_seeds_below_and_p7_needs_no_seed_count(e55):
    d = _good()
    D = np.full(10, -0.001)
    D[:3] = 0.0001
    d["groups5"]["expression_change"] = D  # 7 below: not behind by E49's definition
    D7 = np.full(10, -0.1)
    D7[0] = 1.5  # mean above 0 with one seed ahead: P7 holds ("above 0", no seeds-ahead count)
    d["board4_g"] = D7
    r = e55.decide(d)
    assert r["clauses"]["P4"] and r["clauses"]["P7"] and r["passes"]


# ------------------------------------------------------------------------------------------------ the reading
GROUPS = ("shape_growth", "expression_change", "cell_state", "local_organisation")
NAMES = ("floor", O, ARM, ANTI, AXES)


def _row(seed, shift, s_shift, grp_shift):
    scores, skills, groups = {}, {}, {}
    for n in NAMES:
        scores[n] = 50.0 if n == "floor" else 47.0 + 0.01 * seed + shift.get(n, 0.0)
        sk = dict.fromkeys(METRICS, 0.5)
        sk["d2_shape"] = sk["occupancy_dice"] = 0.5 + s_shift.get(n, 0.0)
        skills[n] = sk
        groups[n] = {g: 0.4 + (s_shift.get(n, 0.0) if g == "shape_growth" else grp_shift.get((n, g), 0.0))
                     for g in GROUPS}
    return {"seed": seed, "scores": scores, "skills": skills, "groups": groups}


def _facts(refused=None):
    names = (O, ARM, ANTI, AXES)
    recs = {n: {"n_jittered": 0, "fallback": False} for n in names}
    recs[ARM]["shape_interp_g"] = 0.05
    return {"board_split": "pseudo", "x": dict.fromkeys(names, "x"), "n": dict.fromkeys(names, 100),
            "refused": {n: (refused or {}).get(n) for n in names[1:]}, "records": recs,
            "hands": {n: {"prefix_mode": "proper", "dice_proper": 0.8, "dice_mirror": 0.7, "dice_gap": 0.1}
                      for n in names},
            "files_identical_across_builds": True, "refusals_identical_across_builds": True}


def _recs(e55, shift=None, s_shift=None, grp_shift=None, task_shift=None, refused_on=None):
    out = {}
    for name in e55.nt.all_tasks():
        rs = []
        for s in SEEDS:
            sh = dict(shift or {})
            for n, v in ((task_shift or {}).get(name) or {}).items():
                sh[n] = sh.get(n, 0.0) + v
            row = _row(s, sh, s_shift or {}, grp_shift or {})
            builds = {b: {"provenance": {"commit": COMMIT, "tree_clean": True}, "build": b, "protocol": "bake",
                          "board": {"split": "pseudo"}, "rows": [json.loads(json.dumps(row))]} for b in e55.BUILDS}
            rs.append({"seed": s, "run_commit": COMMIT, "builds": builds,
                       "facts": _facts((refused_on or {}).get((name, s)))})
        out[name] = rs
    return out


GATE = {"median_g": GATE_MEDIANS, "test_pair_g": [0.0154, 0.0582, 0.0135, 0.0174, 0.0522, 0.0426, 0.0339, 0.0562,
                                                   0.0208, 0.1121]}
GOOD = {ARM: 0.8}
GOOD_S = {ARM: 0.02}


def test_a_pass_names_one_labelled_question(e55):
    res = e55.reading(_recs(e55, GOOD, GOOD_S), [], COMMIT, GATE)
    assert not res["void"] and res["passes"] and res["outcome"] == "PASS"
    q = res["question"]
    assert q["label"] == "diagnostic, selected on the nine tasks"
    assert q["four_truth_mean"] == pytest.approx(0.8) and q["test_like_mean"] == pytest.approx(0.8)
    assert len(q["per_truth"]) == 5 and set(q["gains"]) == set(q["per_truth"]) and q["losses"] == {}
    assert res["p7"]["tasks_at_or_above"] == ["late", "x1", "x2", "x3"]
    text = e55.text(res)
    assert "selected on" in text.split("PASS:")[1].splitlines()[0]
    assert "all eight clauses" in text and "range from zero" in text and "0.0383" in text
    assert "--shape-interp overlay_stretch" in text and "VOID" not in text


def test_four_truth_is_primary_a_known_feature_gain_on_e725_alone_does_not_pass(e55):
    res = e55.reading(_recs(e55, task_shift={"e725": {ARM: 3.0}}, s_shift=GOOD_S), [], COMMIT, GATE)
    assert not res["passes"] and res["clauses"]["P1b"] and not res["clauses"]["P1"]
    assert res["question"] is None and "STOP: no pass, no file." in e55.text(res)


def test_p7_reads_the_three_truth_restricted_merge(e55):
    """A gain only on the large-g tasks' truths still passes P7; a loss there fails it even with gains elsewhere."""
    ts = {t: {ARM: -1.0} for t in ("late", "x1", "x2", "x3")}
    ts.update({t: {ARM: 3.0} for t in ("half", "new1", "new2", "e75")})
    res = e55.reading(_recs(e55, task_shift=ts, s_shift=GOOD_S), [], COMMIT, GATE)
    assert not res["clauses"]["P7"] and res["merged"]["board4_g"]["mean"] == pytest.approx(-1.0)
    assert res["p7"]["board4_below"]["mean"] > 0


def test_a_placebo_refusal_is_a_stop_not_a_void(e55):
    res = e55.reading(_recs(e55, GOOD, GOOD_S, refused_on={("x2", 5): {AXES: "cloud guard"}}), [], COMMIT, GATE)
    assert not res["void"] and not res["passes"] and res["outcome"] == "STOP"
    assert not res["clauses"]["P5"] and res["clauses"]["P6"] and res["unreadable_because_refused"] == [AXES]
    assert "unreadable" in e55.text(res) and "STOP, not void" in e55.text(res)
    res = e55.reading(_recs(e55, GOOD, GOOD_S, refused_on={("late", 2): {ARM: "shape: x"}}), [], COMMIT, GATE)
    assert not res["void"] and not res["clauses"]["P6"] and res["outcome"] == "STOP"


@pytest.mark.parametrize("breach", ["seed", "commit", "unclean", "expression", "count", "expr_skill", "scale",
                                    "builds", "refusal_builds", "gate", "test_board", "test_board_build",
                                    "missing_arm", "gate_record"])
def test_the_void_conditions(e55, breach):
    recs = _recs(e55, GOOD, GOOD_S)
    why, gate = [], GATE
    r = recs["late"][3]
    f = r["facts"]
    if breach == "seed":
        recs["late"].pop(3)
    elif breach == "commit":
        r["builds"]["prefix"]["provenance"]["commit"] = "b" * 40
    elif breach == "unclean":
        r["builds"]["main"]["provenance"]["tree_clean"] = False
    elif breach == "expression":
        f["x"][ARM] = "y"
    elif breach == "count":
        f["n"][ANTI] = 99
    elif breach == "expr_skill":
        r["builds"]["main"]["rows"][0]["skills"][AXES]["mmd_u"] += 2e-9
    elif breach == "scale":
        r["builds"]["main"]["rows"][0]["skills"][ARM]["scale_log_ratio"] += 2e-6
    elif breach == "builds":
        f["files_identical_across_builds"] = False
    elif breach == "refusal_builds":
        f["refusals_identical_across_builds"] = False
    elif breach == "gate":
        why = ["the gate did not pass on late at the run commit"]
    elif breach == "test_board":
        f["board_split"] = "test_interp"
    elif breach == "test_board_build":
        r["builds"]["main"]["board"]["split"] = "test_interp"
    elif breach == "missing_arm":
        r["builds"]["main"]["rows"][0]["scores"].pop(AXES)
    else:
        gate = None
    res = e55.reading(recs, why, COMMIT, gate)
    assert res["void"] and res["why"]
    assert e55.text(res).startswith(f"Row E55 (run commit {COMMIT}): VOID")


def test_tolerances_hold_at_their_edges_and_a_jittered_seed_is_exempt(e55):
    recs = _recs(e55, GOOD, GOOD_S)
    r = recs["late"][3]
    r["builds"]["main"]["rows"][0]["skills"][AXES]["mmd_u"] += 5e-10
    r["builds"]["main"]["rows"][0]["skills"][ARM]["scale_log_ratio"] += 5e-7
    assert not e55.reading(recs, [], COMMIT, GATE)["void"]
    r["builds"]["main"]["rows"][0]["skills"][ARM]["scale_log_ratio"] += 1e-5
    r["facts"]["records"][ARM]["n_jittered"] = 2
    assert not e55.reading(recs, [], COMMIT, GATE)["void"]


def test_the_read_refuses_to_overwrite(e55, tmp_path):
    (tmp_path / "reading.txt").write_text("x")
    with pytest.raises(SystemExit, match="refusing to overwrite"):
        e55.main(["read", "--commit", COMMIT, "--runs", str(tmp_path), "--out", str(tmp_path),
                  "--gate-record", str(tmp_path), "--gate-sha256", "0" * 64])


def test_the_merges_are_five_four_and_test_like(e55):
    m = e55.merges()
    assert set(m["five"]) == {"T2/E8.0.h5ad", "T2/E8.25_late.h5ad", "T2/E8.5.h5ad", "T2/E7.5.h5ad", "T2/E7.25.h5ad"}
    assert m["four"] == {k: v for k, v in m["five"].items() if k != "T2/E7.25.h5ad"}
    assert m["five"]["T2/E8.0.h5ad"] == ["half", "late"] and m["five"]["T2/E8.5.h5ad"] == ["new2", "x2", "x3"]
    assert sorted(n for v in m["test_like"].values() for n in v) == ["half", "late", "new1", "new2"]


# ------------------------------------------------------------------------------------------------ the launcher
def _launch(home: Path, *argv: str, cwd: Path = REPO, **extra) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env.update(HOME=str(home), **extra)
    return subprocess.run(["bash", str(LAUNCHER), *argv], cwd=cwd, env=env, capture_output=True, text=True)


def test_the_launcher_syntax_and_contents():
    assert subprocess.run(["bash", "-n", str(LAUNCHER)]).returncode == 0
    sh = LAUNCHER.read_text()
    assert not re.search(r"(^|[\s;|&(])rm\s", sh)
    assert "OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3" in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh and 'PATH="$HOME/.local/bin:$PATH"' in sh
    assert 'TMPDIR="$S/tmp"' in sh and 'S="$HOME/v/scratch/e55/$RUN"' in sh
    tasks = re.search(r"^TASKS=\(([^)]*)\)$", sh, flags=re.M).group(1).split()
    assert sorted(tasks) == sorted(["half", "late", "new1", "new2", "e75", "e725", "x1", "x2", "x3"])
    assert 3 * len(tasks) <= 32 and 3 * 10 <= 32  # nine (and the gate's ten) processes at 3 threads
    for line in sh.splitlines():
        if "$PY " in line and not line.strip().startswith("#"):
            joined = line
            assert "/dev/null" in joined or line.rstrip().endswith("\\"), line
    assert "upload" not in sh.split("\nset -u", 1)[1]
    assert "tests/test_e55_overlay_stretch.py" in sh


def test_the_launcher_usage_and_odd_names(tmp_path):
    assert _launch(tmp_path, "c" * 40).returncode == 2
    assert _launch(tmp_path, "c" * 40, "a/b").returncode == 90 and not (tmp_path / "v").exists()


def test_the_launcher_refuses_an_existing_run_folder(tmp_path):
    s = tmp_path / "v" / "scratch" / "e55" / "r1"
    s.mkdir(parents=True)
    assert _launch(tmp_path, "c" * 40, "r1").returncode == 91 and list(s.iterdir()) == []


def test_a_refusal_leaves_an_exit_record(tmp_path):
    r = _launch(tmp_path, "c" * 40, "r2", E55_GATE_RECORD=str(tmp_path / "no_record"))
    s = tmp_path / "v" / "scratch" / "e55" / "r2"
    assert r.returncode == 93
    assert (s / "preflight.exit").read_text().strip() == "93" and (s / "done.exit").read_text().strip() == "93"
    work = tmp_path / "notarepo"
    work.mkdir()
    (tmp_path / "record").mkdir()
    (tmp_path / "record" / "gate.json").write_text("{}")
    r = _launch(tmp_path, "c" * 40, "r3", cwd=work, E55_GATE_RECORD=str(tmp_path / "record"), E55_GATE_SHA256="0" * 64,
                GIT_CEILING_DIRECTORIES=str(tmp_path))
    s = tmp_path / "v" / "scratch" / "e55" / "r3"
    assert r.returncode == 98 and (s / "preflight.exit").read_text().strip() == "98"
    assert (s / "done.exit").read_text().strip() == "98" and not (s / "tests.exit").exists()


# ------------------------------------------------------------------------------------------------ the second rulings
def test_s1_the_small_g_range_is_one_merge_and_says_no_gain_at_or_below_zero(e55):
    ts = {"half": {ARM: 0.6}, "new1": {ARM: -0.3}, "new2": {ARM: 0.0}}
    res = e55.reading(_recs(e55, GOOD, GOOD_S, task_shift=ts), [], COMMIT, GATE)
    tp = res["test_pair"]
    # equal-weight merge over three truths, each holding one small-g task: (1.4 + 0.5 + 0.8) / 3
    assert tp["gain_on_small_g_test_like"]["mean"] == pytest.approx((1.4 + 0.5 + 0.8) / 3)
    assert not tp["no_gain_at_the_test"] and tp["small_g_per_task"]["new1"] == pytest.approx(0.5)
    assert "upper end of that range" in e55.text(res)
    ts = {"half": {ARM: -0.8}, "new1": {ARM: -0.8}, "new2": {ARM: -0.8}}
    res = e55.reading(_recs(e55, GOOD, GOOD_S, task_shift=ts), [], COMMIT, GATE)
    assert res["test_pair"]["gain_on_small_g_test_like"]["mean"] == pytest.approx(0.0, abs=1e-12)
    assert res["test_pair"]["no_gain_at_the_test"] and "No gain at the test." in e55.text(res)


def test_s2_void_wins_over_a_stop(e55):
    recs = _recs(e55, GOOD, GOOD_S, refused_on={("x2", 5): {AXES: "cloud guard"}})
    recs["late"][3]["facts"]["n"][ARM] = 99
    res = e55.reading(recs, [], COMMIT, GATE)
    assert res["void"] and "clauses" not in res and "outcome" not in res
    assert e55.text(res).startswith(f"Row E55 (run commit {COMMIT}): VOID")


def test_s3_a_step0_failure_makes_the_run_step_refuse_to_start(e55, tmp_path, monkeypatch, gate_dirs):
    rec, sha, new = gate_dirs
    runs = new.parent
    monkeypatch.setattr(e55, "_refuse_unless_ready", lambda commit: None)
    built = []
    monkeypatch.setattr(e55, "one_seed", lambda *a, **k: built.append(a))
    for name in e55.nt.all_tasks():  # the reproduction gate passes
        p = e55._gate_path(runs, name)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"pass": True, "commit": COMMIT}))
    r = json.loads((new / "x1.json").read_text())
    r["per_seed"]["1000002"]["g"] += 1e-6  # a g off the record
    (new / "x1.json").write_text(json.dumps(r))
    with pytest.raises(SystemExit, match="refused"):
        e55.main(["run", "--task", "late", "--commit", COMMIT, "--runs", str(runs), "--data", str(tmp_path),
                  "--gate-record", str(rec), "--gate-sha256", sha])
    with pytest.raises(SystemExit, match="refused"):
        e55.main(["run", "--task", "late", "--commit", COMMIT, "--runs", str(runs), "--data", str(tmp_path),
                  "--gate-record", str(tmp_path / "none"), "--gate-sha256", sha])
    assert built == [] and not (runs / "run").exists()  # nothing started, nothing counted


def test_s4_a_refusal_of_o_or_a_control_is_not_caught(e55):
    from vec.boards import Board

    board = Board(key="pseudo:T2e/x", task="T2", genes=["a"], needs_coords=True, split="pseudo", min_cells=1,
                  max_cells=50, truth_cells=10)

    def refuse(ctx):
        raise si.ShapeRefused("out of range")

    def guard(ctx):
        raise ValueError("pseudo: cloud guard flags ['coverage']")

    for fn in (refuse, guard):
        with e55.captured() as seen:
            seen["board"] = board
            with pytest.raises(ValueError):
                e55.observed(harness.Candidate(O, frozenset({"T2"}), (), fn, {}), seen).build(_ctx())
    cands = e55.candidates({"board": board, "stage": {}, "x": {}, "file": {}, "refused": {}})
    assert cands[0] is harness.CANDIDATES[harness.FLOOR]  # the floor unwrapped: nothing catches its refusal
    assert [c.name for c in cands] == ["floor", O, ARM, ANTI, AXES]
    assert e55.GATE_RUN == (harness.FLOOR, C0, O)  # the gate builds them first, unwrapped
    recs = _recs(e55, GOOD, GOOD_S)
    recs["x1"] = recs["x1"][:4]  # a control's refusal stopped the task: the missing seeds void the read
    assert e55.reading(recs, [], COMMIT, GATE)["void"]


def test_s5_the_g_tolerance_is_fixed(e55, gate_dirs):
    assert e55.G_TOL == 1e-9
    import inspect

    assert "tol" not in inspect.signature(e55.gate0_problems).parameters
    ap_help = subprocess.run([sys.executable, "-I", "-B", str(REPO / "experiments/coords/e55_overlay_stretch_read.py"),
                              "read", "--help"], capture_output=True, text=True).stdout
    assert "--gate-record" in ap_help and "tol" not in ap_help.lower().replace("tolerance", "x")


def test_s6_a_void_reading_says_a_rerun_needs_a_new_row(e55):
    res = e55.reading(_recs(e55, GOOD, GOOD_S), ["the gate did not pass on late at the run commit"], COMMIT, GATE)
    text = e55.text(res)
    assert res["void"] and "A rerun needs a new row placed before it" in text
    assert "No task or seed is re-run selectively" in text and "stay counted" in text
