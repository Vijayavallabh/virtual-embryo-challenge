"""The off-by-default builder option `build_test.py embryo --shape-interp` (row E49 as drafted) and the harness's
E49 candidates, on synthetic E7.5 and E8.0 stages and a synthetic test board whose limits are never shown (the set-up
of tests/test_occupancy_frame.py). Nothing is read from data/ and nothing is scored.

Pinned: the default build writes the bytes it wrote before the option existed, with the option absent and with
--shape-interp none; a form changes positions only (the same expression rows in the same order, the same count, the
size rule's RMS) and writes a shape_interp block only when given; the option refuses without --method fgw_shape0
--positions floor_exact, and on a test board without --shape-interp-test-ok before anything loads; it composes with
--row-order copy_parity (form first, row order last); the harness has the seven candidates."""

from __future__ import annotations

import sys

import numpy as np
import pytest
import test_occupancy_frame as of
from test_occupancy_frame import (
    ARMS,
    PRE_OPTION_SHA256,
    SECRET_MAX,
    SECRET_TRUTH,
    _build,
    _file,
    _shows,
    embryo_test_board,
    install,
    sha,
)

from vec import harness
from vec import shape_interp as si
from vec import t2_build as tb
from vec.boards import Board

FLOOR_EXACT = ARMS["floor_exact"]


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(of.REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def embryo_env(bt, monkeypatch, own_repo, tmp_path):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path, reads)
    return reads


def test_the_default_path_writes_the_old_bytes_with_the_option_absent_or_none(bt, embryo_env, tmp_path, capsys):
    a, rec_a, text_a = _build(bt, tmp_path, "a", *FLOOR_EXACT)
    b, rec_b, text_b = _build(bt, tmp_path, "b", *FLOOR_EXACT, "--shape-interp", "none")
    assert sha(a) == sha(b) == PRE_OPTION_SHA256["floor_exact"]
    assert text_a == text_b and "shape_interp" not in rec_a and "shape_interp" not in text_a
    assert bt.build_parser().parse_args(["embryo", *FLOOR_EXACT, "--out", "x.h5ad"]).shape_interp == "none"
    m, _, _ = _build(bt, tmp_path, "m", *ARMS["moved"])
    assert sha(m) == PRE_OPTION_SHA256["moved"]


@pytest.mark.parametrize("form", si.FORMS)
def test_a_form_changes_positions_only_and_records_its_block(bt, embryo_env, tmp_path, capsys, form):
    base, _, _ = _build(bt, tmp_path, "base", *FLOOR_EXACT)
    out, rec, text = _build(bt, tmp_path, form, *FLOOR_EXACT, "--shape-interp", form, "--shape-interp-test-ok")
    X0, C0 = _file(base)
    X1, C1 = _file(out)
    np.testing.assert_array_equal(X1, X0)  # the same expression rows in the same order: the pairing of rows kept
    assert len(C1) == len(C0) and not np.allclose(C1, C0)
    assert tb.rms_radius(C1.astype(float)) == pytest.approx(tb.rms_radius(C0.astype(float)), rel=1e-6)  # float32
    blk = rec["shape_interp"]
    assert blk["form"] == form and blk["ruling"] is None and "cloud guard" in blk["order"]
    assert all(not isinstance(v, (dict, list)) for v in blk.values())
    assert rec["record"]["shape_interp"] == form
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)


def test_the_option_refuses_without_the_exact_positions(bt, embryo_env, tmp_path):
    out = tmp_path / "r" / "T2embryo.h5ad"
    for argv in (ARMS["moved"], ["--method", "fgw"]):
        with pytest.raises(SystemExit, match="floor_exact only"):
            bt.main(["embryo", *argv, "--shape-interp", "overlay", "--shape-interp-test-ok", "--out", str(out)])
    assert not out.exists() and embryo_env == []
    with pytest.raises(SystemExit):  # names only: a control is not a choice
        bt.build_parser().parse_args(["embryo", *FLOOR_EXACT, "--shape-interp", "overlay_rot", "--out", "x.h5ad"])


def test_the_option_on_a_test_board_needs_the_named_flag_and_refuses_before_anything_loads(bt, embryo_env, tmp_path,
                                                                                           capsys):
    out = tmp_path / "g" / "T2embryo.h5ad"
    with pytest.raises(SystemExit, match="needs the captain's yes") as e:
        bt.main(["embryo", *FLOOR_EXACT, "--shape-interp", "stretch_ratio", "--out", str(out)])
    assert not out.exists() and embryo_env == []
    assert not _shows(str(e.value) + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)
    with pytest.raises(SystemExit, match="with a form only"):
        bt.main(["embryo", *FLOOR_EXACT, "--shape-interp-test-ok", "--out", str(out)])
    assert bt.build_parser().parse_args(["embryo", *FLOOR_EXACT, "--out", "x.h5ad"]).shape_interp_test_ok is False


def test_the_option_on_a_validation_board_needs_no_flag(bt, embryo_env, tmp_path, monkeypatch):
    val = Board(**{**embryo_test_board().__dict__, "key": "T2:embryo:val_interp", "split": "val_interp"})
    monkeypatch.setattr(bt, "load_boards", lambda: {val.key: val})
    out = tmp_path / "v" / "T2embryo_dry.h5ad"
    bt.main(["embryo", *FLOOR_EXACT, "--board", val.key, "--shape-interp", "stretch_registered", "--out", str(out)])
    assert out.exists()


def test_it_composes_with_the_row_order_in_the_fixed_order(bt, embryo_env, tmp_path, capsys):
    shaped, rec_s, _ = _build(bt, tmp_path, "s", *FLOOR_EXACT, "--shape-interp", "stretch_ratio",
                              "--shape-interp-test-ok")
    both, rec_b, text = _build(bt, tmp_path, "b", *FLOOR_EXACT, "--shape-interp", "stretch_ratio",
                               "--shape-interp-test-ok", "--row-order", "copy_parity", "--row-order-test-ok")
    Xs, Cs = _file(shaped)
    Xb, Cb = _file(both)

    def cells(X, C):
        return sorted(map(tuple, np.hstack([C.astype(np.float64), X.astype(np.float64)]).tolist()))

    assert cells(Xb, Cb) == cells(Xs, Cs)  # the row order permutes the shaped file's cells, nothing else
    assert rec_b["record"]["row_order"] == "copy_parity" and rec_b["shape_interp"]["form"] == "stretch_ratio"
    assert not _shows(text + capsys.readouterr().out, SECRET_MAX, SECRET_TRUTH)


def test_the_harness_has_the_seven_candidates():
    assert len(harness.SHAPE_CANDIDATES) == 7
    for name, form in harness.SHAPE_CANDIDATES.items():
        c = harness.CANDIDATES[name]
        assert c.params["shape_form"] == form and c.params["positions"] == tb.POSITIONS[1]
        assert c.requires is harness.CANDIDATES["fgw_shape0x"].requires
