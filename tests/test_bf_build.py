"""Row R71's files T2embryo_s2x_ovbf.h5ad and T2embryo_s2xs_ovbf.h5ad: the builder's off-by-default `--shape-interp
overlay_bayes` (experiments/builds/build_test.py embryo), on synthetic stages and a synthetic test board whose limits
are never shown. Nothing is read from data/ or boards/, a synthetic TRIES.md stands in for the ledger (one test asks, by
code, whether the real ledger holds row R71 as the builder needs it), and nothing is scored.

Pinned: the option is taken only as O's or SO's command with it and --captain-ruling R71, written under its own file's
name, while row R71's first cell is the captain's order naming BF on O and on SO only and the ruling; every refusal
comes before anything is read; the file keeps its base's rows, order, count, expression and size rule and moves
positions only; its prior is the seven released specimens, read after every refusal; two builds give the same bytes;
its record gains the bayes block and shows no count; and every other embryo command gives what the builder of BASE
(main before the option) gives, run in this same process from git's copy of build_test.py: the same h5ad bytes and
record, or the same refusal. Built through experiments/builds/build_test.py."""

from __future__ import annotations

import hashlib
import json
import sys
import types
from functools import partial
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import git
from scipy import sparse
from test_mix_lean_build import (
    EXTRA,
    PANEL,
    R48_ROW,
    R59_ROW,
    R60_RS_ROW,
    SECRET_MAX,
    SECRET_TRUTH,
    _board,
    _file,
    _shows,
    _stage,
)

from vec import shape_interp as si

REPO = Path(__file__).resolve().parents[1]
BASE = "c15e09de74bff8ee27f2e0d2400b436fd4cf4761"  # main as this branch left it, without --shape-interp overlay_bayes
FLOOR = ["--method", "fgw_shape0", "--positions", "floor_exact"]
OV = ["--shape-interp", "overlay", "--shape-interp-test-ok"]
BF = ["--shape-interp", "overlay_bayes", "--shape-interp-test-ok"]
O_ARGS = [*FLOOR, *OV]
SO_ARGS = [*FLOOR, "--expression", "series", *OV]
BFO_ARGS = [*FLOOR, *BF, "--captain-ruling", "R71"]
BFSO_ARGS = [*FLOOR, "--expression", "series", *BF, "--captain-ruling", "R71"]
BFO, BFSO = "T2embryo_s2x_ovbf.h5ad", "T2embryo_s2xs_ovbf.h5ad"
R71_ROW = ("| R71 (the captain's order of 2026-11-25, ruling, no run: \"My named exception on the embryo's overlay "
           "form and the z-extent frozen rule for each BF file, as R59, with the challenger's fixes.\"; the embryo's "
           "consensus-form files BF on O and on SO only; built only after their gate row reads PASS, only with "
           "`--captain-ruling R71`; not a try) | T2 embryo | synthetic |")
OTHER_STAGES = ("E6.75", "E8.25_late", "E8.75", "E9.5")


def _specimen(name: str, seed: int) -> ad.AnnData:
    """A released specimen the prior reads (coordinates only matter): an elongated blob of its own shape."""
    rng = np.random.default_rng(seed)
    n = 260 + 20 * seed % 7
    C = rng.poisson(1.0, size=(n, len(PANEL) + len(EXTRA))) + 1
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(C).astype(np.float32)),
                   obs=pd.DataFrame({"celltype": ["Epi"] * n, "group": ["Epiblast"] * n},
                                    index=[f"{name}_{i}" for i in range(n)]),
                   var=pd.DataFrame(index=PANEL + EXTRA))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.obsm["spatial_3D"] = (rng.normal(size=(n, 3)) * [3.0 + seed % 3, 1.5, 0.6]).astype(np.float32)
    return a


@pytest.fixture(scope="module")
def stages():
    out = {"E7.25.h5ad": _stage("E7.25", 11), "E7.5.h5ad": _stage("E7.5", 12), "E8.0.h5ad": _stage("E8.0", 13)}
    out.update({f"{s}.h5ad": _specimen(s, 20 + k) for k, s in enumerate(OTHER_STAGES)})
    return out


def install(mod, setattr_, setitem_, shared: Path, reads: list, stages: dict, ledger_rows=None) -> None:
    """Point a build_test module (this checkout's or BASE's) at the synthetic stages, boards, size JSON and ledger; the
    interpolation runs with 150 anchors and the exact assignment on 600 draws (the code paths are the builder's)."""
    ledger_rows = (R48_ROW, R59_ROW, R60_RS_ROW, R71_ROW) if ledger_rows is None else ledger_rows
    boards = {b.key: b for b in (_board(), _board("T2:embryo:val_interp", "val_interp"))}

    def load(path, genes=None, lineage=True, drop_duplicates=False):
        name = Path(path).name
        assert Path(path).parent.name == "T2" and drop_duplicates and genes in (None, PANEL), path
        reads.append(name)
        return stages[name].copy()

    def read(path, genes=None):
        reads.append("released:" + Path(path).name)
        a = stages[Path(path).name]
        return (a if genes is None else a[:, list(genes)]).copy()

    shared.mkdir(parents=True, exist_ok=True)
    size = shared / "size_embryo.json"
    size.write_text(json.dumps({"board": "embryo", "valid": True,
                                "preregistered": {"z_mode": "raw", "valid": True, "pred_rms_full": 6.0}}))
    ledger = shared / "TRIES.md"
    ledger.write_text("| ID | what |\n|---|---|\n| E59 (gate; synthetic) | T2 embryo |\n"
                      + "".join(r + "\n" for r in ledger_rows))
    setattr_(mod, "load_stage", load)
    setattr_(mod, "read_stage", read)
    setattr_(mod, "load_boards", lambda: boards)
    setattr_(mod, "_t2_count", lambda b, n: {"decision": "synthetic"})
    setattr_(mod.t2_build, "interpolate_clouds", partial(mod.t2_build.interpolate_clouds, n_anchor=150))
    setattr_(si, "BAYES_ANCHORS", 600)
    setitem_(mod.SIZES, "embryo", size)
    setattr_(mod, "TRIES_FILE", ledger)


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def env(bt, monkeypatch, own_repo, tmp_path, stages):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, monkeypatch.setitem, tmp_path / "shared", reads, stages)
    return reads


def _build(mod, tmp_path: Path, folder: str, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / folder / name
    mod.main(["embryo", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


# ------------------------------------------------------------------------------------------------ the option
def test_both_files_keep_their_bases_rows_and_expression_and_move_positions_only(bt, env, tmp_path, capsys):
    o, rec_o, _ = _build(bt, tmp_path, "o", "T2embryo.h5ad", *O_ARGS)
    so, rec_so, _ = _build(bt, tmp_path, "so", "T2embryo.h5ad", *SO_ARGS)
    reads_before = list(env)
    bfo, rec_bfo, text_o = _build(bt, tmp_path, "bfo", BFO, *BFO_ARGS)
    bfso, rec_bfso, text_so = _build(bt, tmp_path, "bfso", BFSO, *BFSO_ARGS)
    for base, arm in ((o, bfo), (so, bfso)):
        Xb, Cb, nb, vb = _file(base)
        Xa, Ca, na, va = _file(arm)
        np.testing.assert_array_equal(Xa, Xb)  # the base's rows, order, count and expression
        assert na == nb and va == vb and not np.allclose(Ca, Cb)
        rms = float(np.sqrt(((Ca.astype(np.float64) - Ca.mean(0)) ** 2).sum(1).mean()))
        assert rms == pytest.approx(6.0, rel=1e-3)  # the size rule, after the form
    # the prior: the seven released specimens, each read once, after the brackets
    new = env[len(reads_before):]
    assert new[:2] == ["E7.5.h5ad", "E8.0.h5ad"]
    assert new[2:9] == [f"{s}.h5ad" for s in si.BAYES_PRIOR_STAGES]
    for rec, name, base in ((rec_bfo, BFO, "T2embryo_s2x_ov.h5ad (O)"), (rec_bfso, BFSO, "T2embryo_s2xs_ov.h5ad (SO)")):
        blk = rec["bayes"]
        assert blk["file"] == name and blk["base"] == base and blk["captain_ruling"] == "R71"
        assert blk["command"] == bt.EMBRYO_BF_COMMANDS[name] and "E59" in blk["gate"]
        assert blk["prior"] == "|".join(si.BAYES_PRIOR_STAGES)
        assert all(not isinstance(v, (dict, list)) for v in blk.values())
        s = rec["shape_interp"]
        assert s["form"] == "overlay_bayes" and s["ruling"] == "R71" and s["variant"] == "arm"
        assert s["prior"] == "|".join(si.BAYES_PRIOR_STAGES) and s["flip_copies"] == 8
        assert rec["n_cells"] == "max_cells" and not set(bt.COUNT_FIELDS) & set(rec["record"])
    assert rec_bfo["record"]["target_rms"] == rec_o["record"]["target_rms"] == 6.0
    n = len(_file(bfo)[0])
    shown = text_o + text_so + capsys.readouterr().out
    assert not _shows(shown, SECRET_MAX, SECRET_TRUTH, n)
    again, _, _ = _build(bt, tmp_path, "bfo2", BFO, *BFO_ARGS)  # built again: the same bytes
    assert again.read_bytes() == bfo.read_bytes()


def test_the_launchers_spelling_with_expression_mixture_and_the_test_board_is_the_same_command(bt, env, tmp_path):
    a, _, _ = _build(bt, tmp_path, "a", BFO, *BFO_ARGS)
    b, _, _ = _build(bt, tmp_path, "b", BFO, *FLOOR, "--expression", "mixture", *BF, "--captain-ruling", "R71")
    c, _, _ = _build(bt, tmp_path, "c", BFO, *BFO_ARGS, "--board", "T2:embryo:test_interp", "--seed", "1000001")
    assert a.read_bytes() == b.read_bytes() == c.read_bytes()


def _refused(bt, env, tmp_path, argv, match, name=BFO):
    out = tmp_path / "r" / name
    with pytest.raises(SystemExit, match=match) as e:
        bt.main(["embryo", *argv, "--out", str(out)])
    assert not out.exists() and env == []  # refused before any stage was read
    assert not _shows(str(e.value), SECRET_MAX, SECRET_TRUTH)
    return str(e.value)


@pytest.mark.parametrize("argv", [
    [*FLOOR, *BF],  # no ruling
    [*FLOOR, *BF, "--captain-ruling", "R59"],
    [*FLOOR, *BF, "--captain-ruling", "R48"],
    [*FLOOR, *BF, "--captain-ruling", "R70"],
    [*FLOOR, *BF, "--captain-ruling", "R71a"],
    [*FLOOR, *BF, "--captain-ruling", "E59"],
])
def test_the_ruling_must_be_r71(bt, env, tmp_path, argv):
    msg = _refused(bt, env, tmp_path, argv, "overlay_bayes|goes with --mix-lean later only")
    assert "R71" in msg or "R59 goes with --mix-lean later only" in msg  # R59's own check refuses its ruling first


@pytest.mark.parametrize("argv", [
    [*FLOOR, "--expression", "render_hi", *BF, "--captain-ruling", "R71"],  # O_hi: the depth line, not ordered
    [*FLOOR, "--expression", "series_back", *BF, "--captain-ruling", "R71"],  # W: the depth line, not ordered
    [*FLOOR, "--shape-interp", "overlay_bayes", "--captain-ruling", "R71"],  # no yes flag
    ["--method", "fgw_shape0", *BF, "--captain-ruling", "R71"],
    [*BFO_ARGS, "--row-order", "copy_parity", "--row-order-test-ok"],
    [*BFO_ARGS, "--herd"],
    [*BFO_ARGS, "--draw", "midpoint"],
    [*BFO_ARGS, "--move-coupling", "aligned"],
    [*BFO_ARGS, "--copy-rows", "--copy-rows-test-ok"],
    [*BFO_ARGS, "--seed", "1000002"],
    [*BFO_ARGS, "--neighbourhood", "run_sorted"],
])
def test_any_command_but_os_and_sos_is_refused(bt, env, tmp_path, argv):
    _refused(bt, env, tmp_path, argv, "O's or SO's command|refused")


def test_the_mix_lean_with_r71_and_a_dry_run_are_refused(bt, env, tmp_path):
    _refused(bt, env, tmp_path, [*BFO_ARGS, "--mix-lean", "later"], "R59")
    _refused(bt, env, tmp_path, [*BFO_ARGS, "--board", "T2:embryo:val_interp"], "O's or SO's command",
             "T2embryo_s2x_ovbf_dry.h5ad")


def test_each_command_writes_its_own_file_only(bt, env, tmp_path):
    _refused(bt, env, tmp_path, BFO_ARGS, "--out must name it", BFSO)
    _refused(bt, env, tmp_path, BFSO_ARGS, "--out must name it", BFO)
    _refused(bt, env, tmp_path, BFO_ARGS, "--out must name it", "T2embryo.h5ad")


def test_ruling_r71_without_the_option_is_refused(bt, env, tmp_path):
    for argv in ([*O_ARGS, "--captain-ruling", "R71"], [*SO_ARGS, "--captain-ruling", "R71"],
                 ["--method", "fgw", "--captain-ruling", "R71"]):
        _refused(bt, env, tmp_path, argv, "--captain-ruling R71 goes with --shape-interp overlay_bayes only",
                 "T2embryo.h5ad")


@pytest.mark.parametrize("ledger", [
    (R48_ROW,),  # no row R71
    (R71_ROW, R71_ROW),  # two rows R71
    (R71_ROW.replace("BF on O and on SO only", "BF on O, SO, O_hi and W"),),
    (R71_ROW.replace("`--captain-ruling R71`", "--captain-ruling R71"),),
    (R71_ROW.replace("ruling, no run", "ruling"),),
    (R71_ROW.replace("(the captain's order", "(the AI assistant's order"),),
    (R71_ROW.replace("| R71 ", "| R71a "),),
])
def test_the_ledger_must_hold_exactly_one_row_r71_ordering_the_files(bt, env, tmp_path, monkeypatch, ledger):
    p = tmp_path / "ledger" / "TRIES.md"
    p.parent.mkdir()
    p.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in ledger))
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    _refused(bt, env, tmp_path, BFO_ARGS, "exactly one row R71")


def test_the_real_ledger_holds_row_r71_as_the_builder_needs_it(bt):
    assert bt._bayes_ordered() is True  # R71's real first cell, read by the builder's own check


def test_the_option_is_off_by_default_and_named_in_the_parser(bt):
    a = bt.build_parser().parse_args(["embryo", *O_ARGS, "--out", "x.h5ad"])
    assert a.shape_interp == "overlay" and a.captain_ruling is None
    assert bt.BUILDABLE_SHAPES == (*si.FORMS, *si.COMPOSED_FORMS, "overlay_bayes")
    assert bt.EMBRYO_BF_COMMANDS == {BFO: "embryo " + " ".join(BFO_ARGS), BFSO: "embryo " + " ".join(BFSO_ARGS)}
    for control in si.BAYES_CONTROLS:  # the gate's controls are never built
        with pytest.raises(SystemExit):
            bt.build_parser().parse_args(["embryo", *FLOOR, "--shape-interp", control, "--out", "x.h5ad"])


def test_the_launcher_list_names_this_file():
    sh = (REPO / "experiments" / "external" / "n54_step.sh").read_text()
    assert "tests/test_bf_build.py" in sh and "tests/test_bf_launchers.py" in sh


# ------------------------------------------------------------------------------------------------ every other option
def _module(name: str, rel: str) -> types.ModuleType:
    """A module made from BASE's text of `rel`, registered under `name` (dataclasses look their module up there)."""
    mod = types.ModuleType(name)
    mod.__file__ = str(REPO / rel)
    sys.modules[name] = mod
    exec(compile(git(REPO, "show", f"{BASE}:{rel}") + "\n", str(REPO / rel), "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def base_bt(bt):
    """BASE's builder: every vec module is this checkout's in both builders, so a difference below is the option's or
    a change to vec's other forms."""
    mod = _module("_build_test_at_base_r71", "experiments/builds/build_test.py")
    assert not hasattr(mod, "EMBRYO_BF_ROW") and not hasattr(mod, "BUILDABLE_SHAPES")
    return mod


DRY = "T2:embryo:val_interp"
COMMANDS = {  # name: (argv, file name, what the builder of BASE does with it: "built", "refused" or None for either)
    "O": (O_ARGS, "T2embryo.h5ad", "built"),
    "SO": (SO_ARGS, "T2embryo.h5ad", "built"),
    "O_hi": ([*FLOOR, "--expression", "render_hi", *OV], "T2embryo.h5ad", "built"),
    "W": ([*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R48"], "T2embryo.h5ad", "built"),
    "Ot23": ([*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R59"], "T2embryo.h5ad", "built"),
    "s1": (["--method", "fgw"], "T2embryo.h5ad", "built"),
    "s2": (["--method", "fgw_shape0"], "T2embryo.h5ad", "built"),
    "P": (FLOOR, "T2embryo.h5ad", "built"),
    "overlay_stretch": ([*FLOOR, "--shape-interp", "overlay_stretch", "--shape-interp-test-ok"], "T2embryo.h5ad",
                        None),
    "stretch_ratio": ([*FLOOR, "--shape-interp", "stretch_ratio", "--shape-interp-test-ok"], "T2embryo.h5ad", None),
    "dry_O": ([*FLOOR, "--shape-interp", "overlay", "--board", DRY], "T2embryo_dry.h5ad", "built"),
    "W_no_ruling": ([*FLOOR, "--expression", "series_back", *OV], "T2embryo.h5ad", "refused"),
    "ruling_R48_alone": ([*O_ARGS, "--captain-ruling", "R48"], "T2embryo.h5ad", "refused"),
    "form_needs_ok": ([*FLOOR, "--shape-interp", "overlay"], "T2embryo.h5ad", "refused"),
    "seed_low": ([*O_ARGS, "--seed", "5"], "T2embryo.h5ad", "refused"),
}


def _outcome(mod, tmp: Path, stages, which: str, argv: list[str], name: str):
    """What one builder does with one command: the h5ad's sha256 and its record's text, or ("refused", its message)."""
    with pytest.MonkeyPatch.context() as mp:
        install(mod, mp.setattr, mp.setitem, tmp / "shared", [], stages)
        out = tmp / which / name
        try:
            mod.main(["embryo", *argv, "--out", str(out)])
        except SystemExit as e:
            return "refused", str(e.code).replace(str(out), "<out>")
        except Exception as e:  # a step's own refusal: the same exception from both builders
            return "refused", f"{type(e).__name__}: {e}"
        return hashlib.sha256(out.read_bytes()).hexdigest(), out.with_suffix(".build.json").read_text()


@pytest.mark.parametrize("key", list(COMMANDS))
def test_every_other_command_does_what_base_does(bt, base_bt, own_repo, tmp_path, stages, key):
    argv, name, kind = COMMANDS[key]
    want = _outcome(base_bt, tmp_path, stages, "base", argv, name)
    got = _outcome(bt, tmp_path, stages, "new", argv, name)
    base_did = "refused" if want[0] == "refused" else "built"
    assert kind is None or base_did == kind, (base_did, kind)
    assert got[0] == want[0] and got[1] == want[1]  # the same bytes and record, or the same refusal


def test_base_has_no_bayes_form_and_refuses_ruling_r71(base_bt, own_repo, tmp_path, stages):
    got = _outcome(base_bt, tmp_path, stages, "base", [*O_ARGS, "--captain-ruling", "R71"], "T2embryo.h5ad")
    assert got[0] == "refused"
    with pytest.raises(SystemExit):
        base_bt.build_parser().parse_args(["embryo", *BFO_ARGS, "--out", "x.h5ad"])
