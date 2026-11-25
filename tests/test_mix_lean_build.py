"""Row R59's file T2embryo_s2x_ovt23.h5ad (Ot23): the builder's off-by-default `--mix-lean later`
(experiments/builds/build_test.py embryo), on synthetic stages and a synthetic test board whose limits are never shown.
Nothing is read from data/ or boards/, a synthetic TRIES.md stands in for the ledger (one test asks, by code, whether
the real ledger holds row R59 as the builder needs it), and nothing is scored.

Pinned: the option passes t_mix = 2/3 exactly (the nearest double, the same float as min(1, 1/2 + 1/6)) to
vec.t2_build.interpolate_clouds and nothing else changes in the call (t = 1/2, O's seed, streams, form and positions);
the file keeps O's count, genes and size rule, with a third of its rows real E7.5 cells and two thirds real E8.0 cells;
its record gains the mix_lean block and shows no count; it is taken only on O's command with --captain-ruling R59 while
row R59's first cell is the captain's order naming the option and the ruling, and every refusal comes before anything is
read; and every other embryo command, O's and W's (row R48) among them, gives what the builder of BASE (main before the
option) gives, run in this same process from git's copy of build_test.py: the same h5ad bytes and record, or the same
refusal. Built through experiments/builds/build_test.py."""

from __future__ import annotations

import hashlib
import json
import re
import sys
import types
from fractions import Fraction
from functools import partial
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from conftest import git
from scipy import sparse

from vec.boards import Board

REPO = Path(__file__).resolve().parents[1]
BASE = "ffb2d65776f56973f03f1740430ba5de51f7c9dc"  # main as merged into this branch, without --mix-lean (RS in)
PANEL = [f"g{i}" for i in range(12)]
EXTRA = ["x0", "x1", "x2"]
SECRET_MAX, SECRET_TRUTH = 400, 120  # stand for the test board's max_cells and truth_cells: never shown
WORLD = {"EXE-Endoderm": ("Endoderm", 300, 250, 0), "ExM": ("ExtraEmbryonicMeso", 300, 250, 1),
         "Epi": ("Epiblast", 300, 200, 2), "Meso": ("Mesoderm", 300, 200, 2), "Neural": ("Ectoderm", 0, 100, 1)}
BUMP = {"EXE-Endoderm": [0, 4], "ExM": [1, 5], "Epi": [8, 10], "Meso": [2, 9], "Neural": [3, 11]}
RUN = np.r_[np.full(4, 2.0), np.full(4, 0.5), np.ones(4), np.ones(3)]
DEPTH = {"E7.25": 50.0, "E7.5": 60.0, "E8.0": 200.0}
SHIFT = {"E7.25": -0.1, "E7.5": 0.0, "E8.0": 0.3}
FLOOR = ["--method", "fgw_shape0", "--positions", "floor_exact"]
OV = ["--shape-interp", "overlay", "--shape-interp-test-ok"]
O_ARGS = [*FLOOR, *OV]
OT_ARGS = [*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R59"]
OT_LAUNCHER_ARGS = [*FLOOR, "--expression", "mixture", *OV, "--mix-lean", "later", "--captain-ruling", "R59"]
W_ARGS = [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R48"]
R48_ROW = ("| R48 (the captain's order of 2026-11-24, ruling, no run: \"Build T2embryo_s2x_ovsb.h5ad as named\"; the "
           "embryo builder is to take `--expression series_back` on O's command only with `--captain-ruling R48`; not "
           "a try) | T2 embryo | synthetic |")
R59_ROW = ("| R59 (the captain's order of 2026-11-24, ruling, no run: \"Yes to T2embryo_s2x_ovt23.h5ad (Ot23): my "
           "named exception on the embryo's overlay form for this file, bound to its gate G1.\"; the embryo builder is "
           "to take `--mix-lean later` on O's command only with `--captain-ruling R59`; not a try) | T2 embryo | "
           "synthetic |")
R60_RS_ROW = ("| R60 (the captain's ruling of 2026-11-25, ruling, no run: build T2embryo_s2x_ovrs.h5ad; the embryo "
              "builder is to take `--neighbourhood run_sorted` on O's command only with `--captain-ruling R60`) | T2 "
              "embryo | synthetic |")  # as tests/test_run_sorted.py's synthetic ruling row


# ------------------------------------------------------------------------------------------------ the synthetic world
def _stage(name: str, seed: int) -> ad.AnnData:
    """A bracket as a released file holds it after the duplicate filter (as tests/test_series_back_build.py makes it):
    Poisson counts over the panel and three genes outside it, obs celltype and group, three spatial clusters."""
    rng = np.random.default_rng(seed)
    centres = np.array([[-4.0, 0.0, 0.0], [0.0, 3.0, 0.0], [4.0, 0.0, 1.0]]) * (1.0 + SHIFT.get(name, 0.0))
    C, labs, groups, P = [], [], [], []
    for lab, (grp, n_early, n_late, cl) in WORLD.items():
        n = n_late if name == "E8.0" else n_early
        if not n:
            continue
        prof = np.r_[np.ones(8), np.full(4, 0.3), np.ones(3)]
        prof[BUMP.get(lab, [])] += 2.0
        if name == "E8.0":
            prof = prof * RUN
        c = rng.poisson(prof / prof.sum() * DEPTH[name], size=(n, len(prof)))
        c[c[:, :len(PANEL)].sum(1) == 0, 0] += 1
        C.append(c)
        labs += [lab] * n
        groups += [grp] * n
        P.append(centres[cl] + rng.normal(size=(n, 3)) * [1.5, 1.0, 0.5])
    order = rng.permutation(sum(len(c) for c in C))
    C, labs, groups, P = np.vstack(C)[order], np.asarray(labs)[order], np.asarray(groups)[order], np.vstack(P)[order]
    N = C.sum(1)
    a = ad.AnnData(X=sparse.csr_matrix(np.log1p(1e4 * C / N[:, None]).astype(np.float32)),
                   obs=pd.DataFrame({"celltype": labs, "group": groups}, index=[f"{name}_{i}" for i in range(len(C))]),
                   var=pd.DataFrame(index=PANEL + EXTRA))
    a.layers["counts"] = sparse.csr_matrix(C.astype(np.int32))
    a.obsm["spatial_3D"] = np.asarray(P, dtype=np.float32)
    return a


@pytest.fixture(scope="module")
def stages():
    return {"E7.25.h5ad": _stage("E7.25", 11), "E7.5.h5ad": _stage("E7.5", 12), "E8.0.h5ad": _stage("E8.0", 13)}


def _board(key="T2:embryo:test_interp", split="test_interp") -> Board:
    return Board(key=key, task="T2", genes=list(PANEL), needs_coords=True, split=split, setting="embryo",
                 min_cells=50, max_cells=SECRET_MAX, truth_cells=SECRET_TRUTH)


def install(mod, setattr_, setitem_, shared: Path, reads: list, stages: dict,
            ledger_rows=(R48_ROW, R59_ROW, R60_RS_ROW)) -> None:
    """Point a build_test module (this checkout's or BASE's) at the synthetic stages, boards, size JSON and ledger; the
    interpolation runs with 150 anchors (the code path is the builder's)."""
    boards = {b.key: b for b in (_board(), _board("T2:embryo:val_interp", "val_interp"))}

    def load(path, genes=None, drop_duplicates=False):
        name = Path(path).name
        assert Path(path).parent.name == "T2" and genes == PANEL and drop_duplicates, path
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
    ledger.write_text("| ID | what |\n|---|---|\n| E58 (gate G1; synthetic) | T2 embryo |\n"
                      + "".join(r + "\n" for r in ledger_rows))
    setattr_(mod, "load_stage", load)
    setattr_(mod, "read_stage", read)
    setattr_(mod, "load_boards", lambda: boards)
    setattr_(mod, "_t2_count", lambda b, n: {"decision": "synthetic"})
    setattr_(mod.t2_build, "interpolate_clouds", partial(mod.t2_build.interpolate_clouds, n_anchor=150))
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


def _shows(text: str, *values: int) -> bool:
    """Whether any of the integers appears in the text as a number of its own (not as digits inside a float)."""
    return any(re.search(rf"(?<![\d.]){v}(?![\d.])", text) for v in values)


def _build(mod, tmp_path: Path, name: str, *argv: str) -> tuple[Path, dict, str]:
    out = tmp_path / name / "T2embryo.h5ad"
    mod.main(["embryo", *argv, "--out", str(out)])
    text = out.with_suffix(".build.json").read_text()
    return out, json.loads(text), text


def _file(path: Path) -> tuple[np.ndarray, np.ndarray, list, list]:
    a = ad.read_h5ad(path)
    X = np.asarray(a.X.toarray() if sparse.issparse(a.X) else a.X)
    return X, np.asarray(a.obsm["spatial_3D"]), list(a.obs_names), list(a.var_names)


def _capture(bt, monkeypatch) -> list:
    """Record each interpolation the builder makes: its arguments and its result."""
    seen = []
    inner = bt.t2_build.interpolate_clouds

    def spy(*a, **kw):
        res = inner(*a, **kw)
        seen.append((a, kw, res))
        return res

    monkeypatch.setattr(bt.t2_build, "interpolate_clouds", spy)
    return seen


# ------------------------------------------------------------------------------------------------ the option
def test_t_mix_is_two_thirds_exactly_and_the_only_change_to_os_call(bt, env, tmp_path, monkeypatch):
    seen = _capture(bt, monkeypatch)
    _build(bt, tmp_path, "o", *O_ARGS)
    _build(bt, tmp_path, "ot", *OT_ARGS)
    (a_o, kw_o, _), (a_t, kw_t, res_t) = seen
    assert "t_mix" not in kw_o  # O's call as it was: t_mix unset, so t_mix = t
    assert kw_t["t_mix"] == 2 / 3 == float(Fraction(2, 3)) == min(1.0, 0.5 + 1 / 6) == bt.EMBRYO_OT23_T_MIX
    assert type(kw_t["t_mix"]) is float
    assert {k: v for k, v in kw_t.items() if k != "t_mix"} == kw_o  # form, positions, t_pos, seed, ... as O's
    assert a_t[4] == a_o[4] == 0.5 and a_t[5:] == a_o[5:]  # t = 1/2 still sets the moves; O's count and seed
    for x, y in zip(a_t[:4], a_o[:4], strict=True):
        np.testing.assert_array_equal(x, y)
    assert res_t.record["t_mix"] == 2 / 3 and res_t.record["t"] == 0.5 and res_t.record["t_pos"] == 0.0


def test_the_file_keeps_os_count_genes_and_size_with_two_thirds_e80_rows(bt, env, tmp_path, monkeypatch, capsys,
                                                                         stages):
    seen = _capture(bt, monkeypatch)
    o, rec_o, _ = _build(bt, tmp_path, "o", *O_ARGS)
    ot, rec_t, text = _build(bt, tmp_path, "ot", *OT_ARGS)
    Xo, Co, no, vo = _file(o)
    Xt, Ct, nt_, vt = _file(ot)
    res_o, res_t = seen[0][2], seen[1][2]
    n = len(res_o.rows)
    assert len(Xt) == len(Xo) == len(res_t.rows) == n and vt == vo and nt_ == no  # O's count, genes and names
    assert int((res_t.stage == 0).sum()) == round(n / 3) and int((res_t.stage == 1).sum()) == n - round(n / 3)
    assert int((res_o.stage == 0).sum()) == round(n / 2)
    for C in (Co, Ct):  # the size rule, after the form: the RMS radius about the centroid is the size JSON's target
        rms = float(np.sqrt(((C.astype(np.float64) - C.mean(0)) ** 2).sum(1).mean()))
        assert rms == pytest.approx(6.0, rel=1e-3)
    assert not np.array_equal(Ct, Co) and not np.array_equal(Xt, Xo)
    # the overlay's pool weights follow t_mix: a third on the E7.5 cells, two thirds on the E8.0 cells
    w, n75 = res_t.pool_weights, stages["E7.5.h5ad"].n_obs
    assert len(w) == n75 + stages["E8.0.h5ad"].n_obs
    assert w[:n75].sum() == pytest.approx(1 / 3) and w[n75:].sum() == pytest.approx(2 / 3)
    assert res_o.pool_weights[:n75].sum() == pytest.approx(1 / 2)
    # the record: O's but for t_mix and what follows from the other rows (the union's scale and thinnest SD); the
    # coupling and its rotation R are O's; the mix_lean block; no cell number
    r, ro = rec_t["record"], rec_o["record"]
    assert r["t_mix"] == 2 / 3 and ro["t_mix"] == 0.5 and set(r) == set(ro)
    assert {k for k in r if r[k] != ro[k]} <= {"t_mix", "scale_from_canonical", "shape_interp_thin_sd"}
    assert r["rotation_deg"] == ro["rotation_deg"] and r["coupling_nonzero"] == ro["coupling_nonzero"]
    assert r["target_rms"] == ro["target_rms"] == 6.0 and not set(bt.COUNT_FIELDS) & set(r)
    blk = rec_t["mix_lean"]
    assert blk["captain_ruling"] == "R59" and blk["file"] == "T2embryo_s2x_ovt23.h5ad" and blk["t_mix"] == 2 / 3
    assert blk["command"] == bt.EMBRYO_OT23_COMMAND == "embryo " + " ".join(OT_ARGS)
    assert blk["o_command"] == "embryo " + " ".join(O_ARGS) and "E58" in blk["gate"] and blk["t_mix_exact"] == "2/3"
    assert all(not isinstance(v, (dict, list)) for v in blk.values())
    si_t, si_o = (dict(x["shape_interp"]) for x in (rec_t, rec_o))
    assert si_t.pop("thin_sd") != si_o.pop("thin_sd") and si_t == si_o  # the overlay form: its union moved
    assert rec_t["n_cells"] == rec_o["n_cells"] == "max_cells"
    assert "mix_lean" not in json.dumps(rec_o)
    shown = text + capsys.readouterr().out
    # the count and the two brackets' row numbers (n // 2 is not checked: the record's max_iter is 200 here)
    assert not _shows(shown, SECRET_MAX, SECRET_TRUTH, round(n / 3), n - round(n / 3))
    ot2, _, _ = _build(bt, tmp_path, "ot2", *OT_ARGS)  # built again: the same bytes
    assert ot2.read_bytes() == ot.read_bytes()


def test_the_launchers_spelling_with_expression_mixture_is_the_same_command(bt, env, tmp_path):
    ot, _, _ = _build(bt, tmp_path, "a", *OT_ARGS)
    ot2, _, _ = _build(bt, tmp_path, "b", *OT_LAUNCHER_ARGS)
    assert ot.read_bytes() == ot2.read_bytes()
    ot3, _, _ = _build(bt, tmp_path, "c", *OT_ARGS, "--board", "T2:embryo:test_interp", "--seed", "1000001")
    assert ot3.read_bytes() == ot.read_bytes()


def _refused(bt, env, tmp_path, argv, match, name="T2embryo.h5ad"):
    out = tmp_path / "r" / name
    with pytest.raises(SystemExit, match=match) as e:
        bt.main(["embryo", *argv, "--out", str(out)])
    assert not out.exists() and env == []  # refused before any stage was read
    assert not _shows(str(e.value), SECRET_MAX, SECRET_TRUTH)
    return str(e.value)


@pytest.mark.parametrize("argv", [
    [*O_ARGS, "--mix-lean", "later"],  # no ruling
    [*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R48"],
    [*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R58"],
    [*O_ARGS, "--mix-lean", "later", "--captain-ruling", "R59a"],
    [*O_ARGS, "--mix-lean", "later", "--captain-ruling", "E58"],
])
def test_the_ruling_must_be_r59(bt, env, tmp_path, argv):
    assert "R59" in _refused(bt, env, tmp_path, argv, "--mix-lean later")


@pytest.mark.parametrize("argv", [
    ["--method", "fgw", "--mix-lean", "later", "--captain-ruling", "R59"],
    ["--method", "fgw_shape0", *OV, "--mix-lean", "later", "--captain-ruling", "R59"],
    [*FLOOR, "--mix-lean", "later", "--captain-ruling", "R59"],  # no form
    [*FLOOR, "--shape-interp", "overlay", "--mix-lean", "later", "--captain-ruling", "R59"],  # no yes flag
    [*FLOOR, "--shape-interp", "stretch_ratio", "--shape-interp-test-ok", "--mix-lean", "later", "--captain-ruling",
     "R59"],
    [*FLOOR, "--shape-interp", "overlay_stretch", "--shape-interp-test-ok", "--mix-lean", "later", "--captain-ruling",
     "R59"],
    [*FLOOR, "--expression", "series", *OV, "--mix-lean", "later", "--captain-ruling", "R59"],
    [*FLOOR, "--expression", "series_back", *OV, "--mix-lean", "later", "--captain-ruling", "R59"],
    [*FLOOR, "--expression", "render_hi", *OV, "--mix-lean", "later", "--captain-ruling", "R59"],
    [*OT_ARGS, "--row-order", "copy_parity", "--row-order-test-ok"],
    [*OT_ARGS, "--herd"],
    [*OT_ARGS, "--draw", "midpoint"],
    [*OT_ARGS, "--move-coupling", "aligned"],
    [*OT_ARGS, "--copy-rows", "--copy-rows-test-ok"],
    [*OT_ARGS, "--seed", "1000002"],
    [*OT_ARGS, "--seed", "1000000"],
    [*OT_ARGS, "--neighbourhood", "run_sorted"],  # row E57's option is another file's
])
def test_any_command_but_os_is_refused(bt, env, tmp_path, argv):
    _refused(bt, env, tmp_path, argv, "O's command")


def test_no_dry_run_either(bt, env, tmp_path):
    _refused(bt, env, tmp_path, [*OT_ARGS, "--board", "T2:embryo:val_interp"], "O's command", "T2embryo_dry.h5ad")


def test_the_earlier_lean_is_struck_and_refused_by_the_parser(bt, env, tmp_path, capsys):
    out = tmp_path / "r" / "T2embryo.h5ad"
    for value in ("earlier", "Later", "2/3", ""):
        with pytest.raises(SystemExit) as e:
            bt.main(["embryo", *O_ARGS, "--mix-lean", value, "--captain-ruling", "R59", "--out", str(out)])
        assert e.value.code == 2 and "invalid choice" in capsys.readouterr().err
    assert not out.exists() and env == [] and bt.MIX_LEANS == ("later",)


def test_ruling_r59_without_the_option_is_refused(bt, env, tmp_path):
    for argv in ([*O_ARGS, "--captain-ruling", "R59"], ["--method", "fgw", "--captain-ruling", "R59"],
                 [*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R59"]):
        _refused(bt, env, tmp_path, argv, "--captain-ruling R59 goes with --mix-lean later only")


@pytest.mark.parametrize("ledger", [
    (R48_ROW,),  # no row R59
    (R48_ROW, R59_ROW, R59_ROW),  # two rows R59
    (R59_ROW.replace("`--mix-lean later`", "--mix-lean later"),),  # not in backticks
    (R59_ROW.replace("`--captain-ruling R59`", "--captain-ruling R59"),),
    (R59_ROW.replace("ruling, no run", "ruling"),),
    (R59_ROW.replace("(the captain's order", "(the AI assistant's order"),),
    ("| R59 (the captain's order of 2026-11-24, ruling, no run) | `--mix-lean later` `--captain-ruling R59` |",),
    (R59_ROW.replace("| R59 ", "| R59a "),),
    (R59_ROW.replace("`--mix-lean later`", "`--mix-lean earlier`"),),
])
def test_the_ledger_must_hold_exactly_one_row_r59_naming_the_option_and_the_ruling(bt, env, tmp_path, monkeypatch,
                                                                                   ledger):
    p = tmp_path / "ledger" / "TRIES.md"
    p.parent.mkdir()
    p.write_text("| ID | what |\n|---|---|\n" + "".join(r + "\n" for r in ledger))
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    _refused(bt, env, tmp_path, OT_ARGS, "exactly one row R59")


def test_the_real_ledger_holds_row_r59_as_the_builder_needs_it(bt):
    assert bt._mix_lean_ordered() is True  # R59's real first cell, read by the builder's own check


def test_the_option_is_off_by_default_and_named_in_the_parser(bt):
    a = bt.build_parser().parse_args(["embryo", *O_ARGS, "--out", "x.h5ad"])
    assert a.mix_lean is None and a.captain_ruling is None
    assert bt.EMBRYO_O_COMMAND == "embryo " + " ".join(O_ARGS)
    assert bt.EMBRYO_OT23_COMMAND == bt.EMBRYO_O_COMMAND + " --mix-lean later --captain-ruling R59"
    assert bt.EMBRYO_OT23_OPTION == "`--mix-lean later`" and bt.EMBRYO_OT23_ROW == "R59"
    src = (REPO / "experiments" / "builds" / "build_test.py").read_text()
    assert "mix_lean = _mix_lean_refusals(a)" in src and src.count('"t_mix": EMBRYO_OT23_T_MIX') == 2


def test_the_launcher_list_names_this_file():
    sh = (REPO / "experiments" / "external" / "n54_step.sh").read_text()
    assert "tests/test_mix_lean_build.py" in sh


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
    """BASE's builder: the option changes build_test.py only (vec.t2_build already took t_mix), and every module is
    this checkout's in both builders, so a difference below is the option's."""
    mod = _module("_build_test_at_base_r59", "experiments/builds/build_test.py")
    assert "mix_lean" not in mod.__doc__ and not hasattr(mod, "EMBRYO_OT23_T_MIX")
    return mod


DRY = "T2:embryo:val_interp"
# name: (argv, file name, what the builder of BASE does with it: "built", "refused" or None for either)
COMMANDS = {
    "O": (O_ARGS, "T2embryo.h5ad", "built"),
    "O_launcher": ([*FLOOR, "--expression", "mixture", *OV], "T2embryo.h5ad", "built"),
    "W": (W_ARGS, "T2embryo.h5ad", "built"),
    # row E57's file: on these synthetic groups vec.run_sort refuses (its pools); both builders must do the same
    "RS": ([*O_ARGS, "--neighbourhood", "run_sorted", "--captain-ruling", "R60"], "T2embryo.h5ad", None),
    "RS_no_ruling": ([*O_ARGS, "--neighbourhood", "run_sorted"], "T2embryo.h5ad", "refused"),
    "s1": (["--method", "fgw"], "T2embryo.h5ad", "built"),
    "s2": (["--method", "fgw_shape0"], "T2embryo.h5ad", "built"),
    "P": (FLOOR, "T2embryo.h5ad", "built"),
    "SO": ([*FLOOR, "--expression", "series", *OV], "T2embryo.h5ad", "built"),
    "O_hi": ([*FLOOR, "--expression", "render_hi", *OV], "T2embryo.h5ad", "built"),
    "overlay_stretch": ([*FLOOR, "--shape-interp", "overlay_stretch", "--shape-interp-test-ok"], "T2embryo.h5ad",
                        None),
    "series": (["--method", "fgw", "--expression", "series"], "T2embryo.h5ad", "built"),
    "midpoint": (["--method", "fgw", "--draw", "midpoint"], "T2embryo.h5ad", "built"),
    "copy_parity": ([*FLOOR, "--row-order", "copy_parity", "--row-order-test-ok"], "T2embryo.h5ad", None),
    "copy_rows": ([*FLOOR, "--copy-rows", "--copy-rows-test-ok"], "T2embryo.h5ad", None),
    "seed": ([*O_ARGS, "--seed", "1000003"], "T2embryo.h5ad", "built"),
    "dry_O": ([*FLOOR, "--shape-interp", "overlay", "--board", DRY], "T2embryo_dry.h5ad", "built"),
    # the refusals as they were, row R48's among them
    "W_no_ruling": ([*FLOOR, "--expression", "series_back", *OV], "T2embryo.h5ad", "refused"),
    "W_other_ruling": ([*FLOOR, "--expression", "series_back", *OV, "--captain-ruling", "R47"], "T2embryo.h5ad",
                       "refused"),
    "W_seed": ([*W_ARGS, "--seed", "1000002"], "T2embryo.h5ad", "refused"),
    "ruling_R48_alone": ([*O_ARGS, "--captain-ruling", "R48"], "T2embryo.h5ad", "refused"),
    "ruling_R60_alone": ([*O_ARGS, "--captain-ruling", "R60"], "T2embryo.h5ad", "refused"),
    "form_needs_ok": ([*FLOOR, "--shape-interp", "overlay"], "T2embryo.h5ad", "refused"),
    "positions": (["--method", "fgw", "--positions", "floor_exact"], "T2embryo.h5ad", "refused"),
    "seed_low": ([*O_ARGS, "--seed", "5"], "T2embryo.h5ad", "refused"),
    "dry_name": ([*FLOOR, "--board", DRY], "T2embryo.h5ad", "refused"),
    "test_name": (FLOOR, "T2embryo_dry.h5ad", "refused"),
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
        except Exception as e:  # a step's own refusal (vec.run_sort's, say): the same exception from both builders
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


def test_base_has_no_mix_lean_and_refuses_ruling_r59(base_bt, own_repo, tmp_path, stages):
    got = _outcome(base_bt, tmp_path, stages, "base", [*O_ARGS, "--captain-ruling", "R59"], "T2embryo.h5ad")
    assert got[0] == "refused" and "series_back only" in got[1]
    with pytest.raises(SystemExit):
        base_bt.build_parser().parse_args(["embryo", *OT_ARGS, "--out", "x.h5ad"])
