"""Row R68's four hedge files: the builder experiments/builds/hedge_fixed.py and its DGX launcher
experiments/builds/hedge_fixed_line.sh (HEDGE_FIXED_design.md, section 3).

The builder's pinned table is held to the records' committed lists and the design's pins. Synthetic files (small stages
written here, a temporary board folder, a temporary ledger holding row R68's first cell, a clean throwaway repository
for the tree and the stamps) run each pick, each refusal (exit 2, nothing written), the writers' own refusal (exit 1)
and two builds giving one sha256. The launcher runs for real in throwaway repositories under a temporary home with a
stand-in `uv`: every refusal path, a whole run of each name, and each failure after the gate. Nothing here reads the
repository's own ledger or boards, or any built file of a test board."""

import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
from conftest import git_env
from scipy import sparse

REPO = Path(__file__).resolve().parents[1]
BUILDER = REPO / "experiments" / "builds" / "hedge_fixed.py"
LAUNCHER = REPO / "experiments" / "builds" / "hedge_fixed_line.sh"
RESULTS = REPO / "experiments" / "builds" / "results"


def _load():
    spec = importlib.util.spec_from_file_location("hedge_fixed", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["hedge_fixed"] = mod
    spec.loader.exec_module(mod)
    return mod


hf = _load()

# Row R68's first cell, as experiments/harness/TRIES.md holds it (the captain's order of 2026-11-25)
REAL_R68 = ("R68 (the captain's order of 2026-11-25, ruling, no run: \"Yes to the four hedge files of "
            "HEDGE_FIXED_design.md as named.\" (his pick \"All four\" in a dialog); the hedge files with parents fixed "
            "by recipe by a score-blind designer, 50:50, each row whole from one parent: `T3_c2_s2_mxe2.h5ad` (C2 and "
            "xe2), `T3_s2_mx_xe2_ck.h5ad` (xe2 and ck), `T2embryo_s2x_ovmxso.h5ad` (O and SO per cell), "
            "`T2embryo_s2x_ovbkso.h5ad` (O and SO by spatial blocks), each built by "
            "`experiments/builds/hedge_fixed.py` only with `--captain-ruling R68`; decision log 2026-11-25 (qqq); "
            "not a try)")
ROW_R68 = f"| {REAL_R68} | T3, T2 | ruling |\n"
T3, EMB = "T3:bcatenin", "T2:embryo:test_interp"
# name -> (anchor, first, second, board, pick), as the design's section 3 commands them
ROLES = {"T3_c2_s2_mxe2.h5ad": ("C2", "C2", "xe2", T3, "parity"),
         "T3_s2_mx_xe2_ck.h5ad": ("C2", "xe2", "ck", T3, "parity"),
         "T2embryo_s2x_ovmxso.h5ad": ("O", "O", "SO", EMB, "parity"),
         "T2embryo_s2x_ovbkso.h5ad": ("O", "O", "SO", EMB, "blocks")}
NAMES = list(ROLES)
# the design's pins (HEDGE_FIXED_design.md: "47111ac1...", "77ff0b94...", ...) and the records' lines naming them
DESIGN_PINS = {"C2": "47111ac1", "xe2": "77ff0b94", "ck": "e470aaee", "O": "1cfc91be", "SO": "77062082"}
RECORDS = {"C2": ("t3two_fb01487", "c2-s2.h5ad"), "xe2": ("t3k23_49f3fd5", "xall-s2.h5ad"),
           "ck": ("t3ck_20a0e1a", "ck-s2.h5ad"), "O": ("embovsb_20b2381_w1", "O/T2embryo_s2x_ov.h5ad"),
           "SO": ("embovsb_20b2381_w1", "SO/T2embryo_s2xs_ov.h5ad")}
SHA_VARS = {"C2": "C2_SHA", "xe2": "XE2_SHA", "ck": "CK_SHA", "O": "O_SHA", "SO": "SO_SHA"}
GENES = [f"g{i}" for i in range(12)]
N = 64


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


# ------------------------------------------------------------------ the pinned table


def test_the_pins_are_the_records_lines_and_the_designs():
    for key, (record, file) in RECORDS.items():
        sha = getattr(hf, SHA_VARS[key])
        assert _sums(f"{record}.SHA256SUMS")[file] == sha, key
        assert sha.startswith(DESIGN_PINS[key]) and re.fullmatch(r"[0-9a-f]{64}", sha), key
        assert hf.PARENTS[sha] == (key, file, record)
    assert len(hf.PARENTS) == 5


def test_the_table_is_the_designs_four_files():
    sha = {k: getattr(hf, v) for k, v in SHA_VARS.items()}
    assert list(hf.TABLE) == NAMES
    for name, (anchor, first, second, board, pick) in ROLES.items():
        assert hf.TABLE[name] == hf.Hedge(board, (sha[first], sha[second]), pick, sha[anchor]), name
        assert name in re.findall(r"`([^`]*)`", REAL_R68)
    assert hf.RULING_ROW == "R68" and hf.SCRIPT == "experiments/builds/hedge_fixed.py"
    assert {h.board for h in hf.TABLE.values()} == set(hf.WRITERS) == {T3, EMB}


def test_the_builder_puts_its_own_src_first_and_imports_vec_provenance_before_anything_else():
    src = BUILDER.read_text()
    code = src[src.index('"""', src.index('"""') + 3) + 3:]  # after the module docstring
    imports = [ln for ln in code.splitlines() if re.match(r"(import|from) \S", ln)]
    assert imports[:4] == ["from __future__ import annotations", "import sys", "from pathlib import Path",
                           imports[3]] and imports[3].startswith("import vec.provenance  #")
    assert code.index('sys.path.insert(0, str(REPO / "src"))') < code.index("import vec.provenance")
    assert 'REPO = Path(__file__).resolve().parents[2]' in code
    assert "captain_order_names(RULING_ROW, f\"`{name}`\", p)" in code  # the check is build_outside_form's own


# ------------------------------------------------------------------ synthetic files


def _write(path: Path, X, coords, *, names=None, genes=GENES, obs_cols=None, more_obsm=None) -> str:
    """A synthetic file: make_submission's form (names cell_i, spatial_3D only) unless told otherwise; its sha256."""
    a = ad.AnnData(X=X)
    a.var_names = list(genes)
    a.obs_names = list(names) if names is not None else [f"cell_{i}" for i in range(X.shape[0])]
    if coords is not None:
        a.obsm["spatial_3D"] = coords
    for k, v in (obs_cols or {}).items():
        a.obs[k] = v
    for k, v in (more_obsm or {}).items():
        a.obsm[k] = v
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path, compression="gzip")
    return _sha(path)


def _boards(tmp_path: Path, genes=GENES, t3_min: int = 10) -> Path:
    d = tmp_path / "boards"
    d.mkdir(exist_ok=True)
    (d / "genes.txt").write_text("\n".join(genes) + "\n")
    spec = {"genes_file": "genes.txt", "needs_coords": True, "max_cells": None, "obsm_required": ["spatial_3D"]}
    (d / "index.json").write_text(json.dumps({T3: {**spec, "task": "T3", "min_cells": t3_min},
                                              EMB: {**spec, "task": "T2", "min_cells": 10}}))
    return d


def _cloud(rng, n: int = N) -> np.ndarray:
    return (rng.normal(size=(n, 3)) * [3.0, 2.0, 1.0] + [4.0, -1.0, 2.0]).astype(np.float32)


@pytest.fixture
def world(tmp_path, monkeypatch, own_repo):
    """Five synthetic parents (C2, xe2 and ck on one cloud, O and SO on another), the table pinned to their sha256,
    a temporary board folder and a ledger holding row R68; the tree and the stamps are own_repo's (clean)."""
    rng = np.random.default_rng(68)
    clouds = {T3: _cloud(rng), EMB: _cloud(rng)}
    xs, paths, shas = {}, {}, {}
    for key, board in (("C2", T3), ("xe2", T3), ("ck", T3), ("O", EMB), ("SO", EMB)):
        xs[key] = rng.gamma(1.0, 1.0, size=(N, len(GENES))).astype(np.float32)
        xs[key][rng.random(xs[key].shape) < 0.3] = 0.0
        paths[key] = tmp_path / "parents" / f"{key}.h5ad"
        shas[key] = _write(paths[key], xs[key], clouds[board])
    monkeypatch.setattr(hf, "TABLE", {name: hf.Hedge(board, (shas[f], shas[s]), pick, shas[a])
                                      for name, (a, f, s, board, pick) in ROLES.items()})
    monkeypatch.setattr(hf, "PARENTS", {shas[k]: (k, f"{k}.h5ad", f"rec_{k}") for k in shas})
    monkeypatch.setattr(hf, "BOARDS_DIR", _boards(tmp_path))
    ledger = tmp_path / "TRIES.md"
    ledger.write_text(f"# ledger\n{ROW_R68}")
    monkeypatch.setattr(hf, "TRIES_FILE", ledger)
    return SimpleNamespace(tmp=tmp_path, xs=xs, paths=paths, shas=shas, clouds=clouds, ledger=ledger, repo=own_repo,
                           monkeypatch=monkeypatch)


def _argv(w, name: str, out: Path | None = None, ruling: str | None = "R68", **given) -> list[str]:
    a, f, s, _, _ = ROLES.get(name, ROLES[NAMES[0]])
    roles = {"anchor": w.paths[a], "first": w.paths[f], "second": w.paths[s], **given}
    out = w.tmp / "out" / "a" / name if out is None else out
    argv = ["build", "--name", name, *(x for r in ("anchor", "first", "second") for x in (f"--{r}", str(roles[r]))),
            "--out", str(out)]
    return argv + (["--captain-ruling", ruling] if ruling is not None else [])


def _pin(w, name: str, **keys_or_paths) -> None:
    """Re-pin `name`'s anchor, first or second to the sha256 of the given files (a world key or a path)."""
    h = hf.TABLE[name]
    sha = {r: _sha(w.paths[v] if isinstance(v, str) else v) for r, v in keys_or_paths.items()}
    w.monkeypatch.setitem(hf.TABLE, name, h._replace(anchor=sha.get("anchor", h.anchor),
                                                     parents=(sha.get("first", h.parents[0]),
                                                              sha.get("second", h.parents[1]))))


def _blocks_by_svd(coords) -> np.ndarray:
    """The blocks pick computed apart from the builder: SVD axes, each signed by its largest loading, RMS units."""
    C = np.asarray(coords, dtype=np.float64)
    X = C - C.mean(0)
    Vt = np.linalg.svd(X, full_matrices=False)[2]
    for k in range(3):
        if Vt[k, np.argmax(np.abs(Vt[k]))] < 0:
            Vt[k] = -Vt[k]
    Z = X @ Vt.T / np.sqrt((X ** 2).sum(1).mean())
    return np.floor(Z).astype(int).sum(1) % 2 == 1


def _nothing_written(w) -> None:
    out = w.tmp / "out"
    assert not out.exists() or not [p for p in out.rglob("*") if p.is_file()]


# ------------------------------------------------------------------ the picks


@pytest.mark.parametrize("name", NAMES)
def test_each_row_comes_whole_from_its_picked_parent_on_the_anchors_rows(world, name, capsys):
    out = world.tmp / "out" / "a" / name
    assert hf.main(_argv(world, name, out)) == 0
    printed = capsys.readouterr().out.split()
    assert printed == ["ok", _sha(out)]
    a_key, f_key, s_key, board, pick = ROLES[name]
    anchor = ad.read_h5ad(world.paths[a_key])
    second = np.arange(N) % 2 == 1 if pick == "parity" else _blocks_by_svd(anchor.obsm["spatial_3D"])
    assert 0 < second.sum() < N  # both parents give rows
    want = np.where(second[:, None], world.xs[s_key], world.xs[f_key])
    b = ad.read_h5ad(out)
    assert isinstance(b.X, np.ndarray) and b.X.dtype == np.float32
    assert b.X.tobytes() == want.tobytes()  # each row bit for bit from its parent
    for i in range(N):
        src = world.xs[s_key] if second[i] else world.xs[f_key]
        assert b.X[i].tobytes() == src[i].tobytes()
    assert list(b.obs_names) == list(anchor.obs_names) == [f"cell_{i}" for i in range(N)]
    assert b.obs.shape[1] == 0 and list(b.var_names) == GENES and set(b.obsm) == {"spatial_3D"}
    assert b.obsm["spatial_3D"].dtype == np.float32
    assert b.obsm["spatial_3D"].tobytes() == anchor.obsm["spatial_3D"].tobytes()
    assert (out.parent / f"{name}.stamp.json").is_file()  # the board's writer stamps the file


def test_parity_takes_the_first_parent_on_even_rows_and_the_second_on_odd_rows(world):
    h = hf.TABLE[NAMES[0]]
    anchor = ad.read_h5ad(world.paths["C2"])
    got = hf.second_rows(h, anchor)
    assert got.tolist() == [i % 2 == 1 for i in range(N)]


def test_blocks_match_a_separate_computation_and_take_no_count(world):
    h = hf.TABLE["T2embryo_s2x_ovbkso.h5ad"]
    anchor = ad.read_h5ad(world.paths["O"])
    got = hf.second_rows(h, anchor)
    assert got.tolist() == _blocks_by_svd(anchor.obsm["spatial_3D"]).tolist()
    # a row's pick is a function of its own canonical coordinates only: the floor-sum parity there
    Z = hf.canonical(anchor.obsm["spatial_3D"])
    assert got.tolist() == [bool((np.floor(z).astype(int).sum()) % 2) for z in Z]
    assert "n_obs" not in hf.canonical.__code__.co_names and "n_obs" not in hf.PICK_RULES["blocks"]


def test_the_blocks_frame_is_the_principal_axes_signed_by_their_largest_loading_in_rms_units():
    """A full grid has its centroid at 0 and a diagonal covariance with variances x > y > z: after any rotation R
    (and a shift), the canonical coordinates are the grid's, each axis signed by the largest entry of R's column, over
    the grid's RMS radius."""
    g = np.stack(np.meshgrid(np.arange(-3, 4) * 1.3, np.arange(-2, 3) * 0.7, np.arange(-1, 2) * 0.3, indexing="ij"),
                 -1).reshape(-1, 3)
    r = np.sqrt((g ** 2).sum(1).mean())
    for seed in range(4):
        R = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))[0]
        C = g @ R.T + [5.0, -2.0, 1.5]
        s = np.array([np.sign(R[np.argmax(np.abs(R[:, k])), k]) for k in range(3)])
        assert np.allclose(hf.canonical(C), g * s / r, atol=1e-9), seed
        assert np.allclose(hf.canonical(-C), -g * s / r, atol=1e-9), seed  # mirrored: every axis keeps its sign


# ------------------------------------------------------------------ two builds, the sidecar


@pytest.mark.parametrize("name", NAMES)
def test_two_builds_give_one_sha256_and_one_sidecar(world, name, capsys):
    got = []
    for folder in ("a", "b"):
        out = world.tmp / "out" / folder / name
        assert hf.main(_argv(world, name, out)) == 0
        printed = capsys.readouterr().out.split()
        assert printed == ["ok", _sha(out)]
        got.append((printed[1], out.with_suffix(".json").read_text()))
    assert got[0][0] == got[1][0]
    assert got[0][1] == got[1][1]  # the sidecar names no output path, time or count


def _numbers(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _numbers(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _numbers(v)
    elif isinstance(obj, int | float) and not isinstance(obj, bool):
        yield obj


@pytest.mark.parametrize("name", NAMES)
def test_the_sidecar_holds_names_sha256s_the_pick_the_scripts_sha256_and_the_row_never_a_count(world, name):
    out = world.tmp / "out" / "a" / name
    assert hf.main(_argv(world, name, out)) == 0
    side = json.loads(out.with_suffix(".json").read_text())
    a_key, f_key, s_key, board, pick = ROLES[name]
    assert not list(_numbers(side))
    assert side["name"] == name and side["board"] == board and side["pick"] == pick
    assert side["pick_rule"] == hf.PICK_RULES[pick] and side["writer"] == hf.WRITERS[board]
    assert side["script"] == "experiments/builds/hedge_fixed.py" and side["script_sha256"] == _sha(BUILDER)
    for role, key in (("anchor", a_key), ("first", f_key), ("second", s_key)):
        assert side[role]["sha256"] == world.shas[key] and side[role]["parent"] == key
        assert side[role]["path"] == str(world.paths[key])
    assert side["captain_ruling"] == {"row": "R68", "ledger": "experiments/harness/TRIES.md",
                                      "ledger_sha256": _sha(world.ledger)}
    assert side["output"] == {"file": name, "sha256": _sha(out), "format_check": "passed"}
    assert side["tree_at_start"]["clean"] is True
    assert set(side) == {"script", "script_sha256", "design", "name", "board", "pick", "pick_rule", "anchor", "first",
                         "second", "obs_obsm", "writer", "captain_ruling", "tree_at_start", "output"}


# ------------------------------------------------------------------ refusals: exit 2, nothing written


def _refused(w, argv, capsys, says: str) -> None:
    assert hf.main(argv) == 2
    io = capsys.readouterr()
    assert io.out == "" and io.err.startswith("refused, nothing written: ") and says in io.err, io.err
    _nothing_written(w)


@pytest.mark.parametrize("ruling, says", [
    (None, "no --captain-ruling R68"), ("R67", "'R67' is refused"), ("r68", "'r68' is refused"),
    ("R68 ", "'R68 ' is refused"), ("R22", "'R22' is refused")])
def test_no_or_another_ruling_is_refused(world, capsys, ruling, says):
    _refused(world, _argv(world, NAMES[0], ruling=ruling), capsys, says)


@pytest.mark.parametrize("name", ["T3_c2_s2.h5ad", "T2embryo_s2x_ov.h5ad", "t3_c2_s2_mxe2.h5ad",
                                  "T3_c2_s2_mxe2", "T2embryo_s2x_ovmxso.h5ad ", "T2embryo_s2xs_ov.h5ad", ""])
def test_an_unknown_name_is_refused(world, capsys, name):
    out = world.tmp / "out" / "a" / (name.strip() or "x.h5ad")
    _refused(world, _argv(world, name, out), capsys, "is not one of row R68's four files")


@pytest.mark.parametrize("tries", [
    "",
    ROW_R68 + ROW_R68,  # two rows with one ID
    ROW_R68.replace("| R68 (", "| R68b ("),  # another ID
    ROW_R68.replace("R68 (the captain's order", "R68 (a proposal"),
    ROW_R68.replace("R68 (the captain's order", "R68 awaits (the captain's order"),
    ROW_R68.replace("ruling, no run", "a try"),
    ROW_R68.replace("`T3_c2_s2_mxe2.h5ad`", "T3_c2_s2_mxe2.h5ad"),  # not backticked
    ROW_R68.replace("`T3_c2_s2_mxe2.h5ad`", "`T3_c2_s2_mxe2.h5ad x`"),  # inside a longer string
    ROW_R68.replace("`T3_c2_s2_mxe2.h5ad` (C2 and xe2), ", "") + "| R69 | `T3_c2_s2_mxe2.h5ad` |\n",
    "| R68 (the captain's order of X, ruling, no run: `--captain-ruling R68` `experiments/builds/hedge_fixed.py`) "
    "| `T3_c2_s2_mxe2.h5ad` |\n",  # the name outside the first cell only
    ROW_R68.replace("`experiments/builds/hedge_fixed.py`", "`experiments/builds/hedge_mix.py`"),
    ROW_R68.replace("`--captain-ruling R68`", "`--captain-ruling R67`"),
    ROW_R68.replace("`--captain-ruling R68`", "`--captain-ruling R68 now`"),
])
def test_a_ledger_without_row_r68_naming_the_file_is_refused(world, capsys, tries):
    world.ledger.write_text(f"# ledger\n{tries}")
    _refused(world, _argv(world, NAMES[0]), capsys, "holds no single row R68")


def test_a_short_row_in_its_own_words_passes(world):
    world.ledger.write_text("| R68 (the captain's order of X, ruling, no run: `T3_c2_s2_mxe2.h5ad` by "
                            "`experiments/builds/hedge_fixed.py` with `--captain-ruling R68`) | T3 |\n")
    assert hf.main(_argv(world, NAMES[0])) == 0
    world.ledger.write_text("")
    with pytest.raises(hf.Refused):
        hf.check_ruling("R68", NAMES[0])


def test_an_out_that_is_misnamed_or_exists_is_refused_and_left_as_it_was(world, capsys):
    name = NAMES[2]
    _refused(world, _argv(world, name, world.tmp / "out" / "a" / "T2embryo_s2x_ovbkso.h5ad"), capsys, "is not named")
    for taken in (name, name.replace(".h5ad", ".json"), f"{name}.stamp.json"):  # the file, its sidecar, its stamp
        d = world.tmp / f"taken_{taken}"
        d.mkdir()
        q = d / taken
        q.write_text("keep")
        assert hf.main(_argv(world, name, d / name)) == 2
        assert "exists: nothing is replaced" in capsys.readouterr().err
        assert [p.name for p in d.iterdir()] == [q.name] and q.read_text() == "keep"


@pytest.mark.parametrize("dirt", ["untracked", "tracked"])
def test_an_unclean_tree_is_refused(world, capsys, dirt):
    if dirt == "untracked":
        (world.repo / "stray.py").write_text("x = 2\n")
    else:
        (world.repo / "method.py").write_text("x = 2\n")
    _refused(world, _argv(world, NAMES[3]), capsys, "the tree is not clean")


@pytest.mark.parametrize("name, given, says", [
    ("T3_c2_s2_mxe2.h5ad", {"anchor": "xe2"}, "--anchor"),
    ("T3_s2_mx_xe2_ck.h5ad", {"anchor": "ck"}, "--anchor"),
    ("T3_c2_s2_mxe2.h5ad", {"first": "ck"}, "--first"),
    ("T3_s2_mx_xe2_ck.h5ad", {"second": "C2"}, "--second"),
    ("T2embryo_s2x_ovmxso.h5ad", {"second": "O"}, "--second"),
    ("T2embryo_s2x_ovbkso.h5ad", {"anchor": "SO"}, "--anchor"),
    ("T3_c2_s2_mxe2.h5ad", {"first": "xe2", "second": "C2"}, "the parents are in the wrong order"),
    ("T3_s2_mx_xe2_ck.h5ad", {"first": "ck", "second": "xe2"}, "the parents are in the wrong order"),
    ("T2embryo_s2x_ovmxso.h5ad", {"first": "SO", "second": "O"}, "the parents are in the wrong order"),
    ("T2embryo_s2x_ovbkso.h5ad", {"first": "SO", "second": "O"}, "the parents are in the wrong order"),
])
def test_a_file_off_its_pin_or_parents_in_the_wrong_order_are_refused(world, capsys, name, given, says):
    _refused(world, _argv(world, name, **{r: world.paths[k] for r, k in given.items()}), capsys, says)


def test_a_missing_file_is_refused(world, capsys):
    _refused(world, _argv(world, NAMES[1], second=world.tmp / "parents" / "none.h5ad"), capsys, "is not a file")


def _odd_parent(w, name: str, role: str, **kw) -> Path:
    """A parent like the world's (the role's own X, the anchor's cloud) but for `kw`, pinned for `name`."""
    a_key, f_key, s_key, board, _ = ROLES[name]
    key = {"first": f_key, "second": s_key}[role]
    X = kw.pop("X", w.xs[key])
    coords = kw.pop("coords", w.clouds[board])
    p = w.tmp / "odd" / f"{role}.h5ad"
    _write(p, X, coords, **kw)
    _pin(w, name, **{role: p})
    return p


def test_a_parent_of_another_count_is_refused(world, capsys):
    name = NAMES[0]
    p = _odd_parent(world, name, "second", X=world.xs["xe2"][:-1], coords=world.clouds[T3][:-1])
    _refused(world, _argv(world, name, second=p), capsys, "its n_obs is not the anchor's")


@pytest.mark.parametrize("change, says", [
    ({"names": [f"cell_{i}" for i in range(N - 1)] + ["cell_x"]}, "its obs_names are not the anchor's"),
    ({"names": [f"cell_{i}" for i in reversed(range(N))]}, "its obs_names are not the anchor's"),
    ({"genes": GENES[::-1]}, "its var_names are not the anchor's"),
    ({"genes": GENES[:-1] + ["other"]}, "its var_names are not the anchor's"),
    ("nudged", "its obsm['spatial_3D'] is not the anchor's"),
    ("float64", "its obsm['spatial_3D'] is not the anchor's"),
    ({"coords": None}, "its obsm['spatial_3D'] is not the anchor's"),
])
@pytest.mark.parametrize("name, role", [("T3_s2_mx_xe2_ck.h5ad", "first"), ("T2embryo_s2x_ovbkso.h5ad", "second")])
def test_a_parent_not_on_the_anchors_rows_names_genes_or_positions_is_refused(world, capsys, change, says, name,
                                                                                role):
    board = ROLES[name][3]
    if change == "nudged":
        c = world.clouds[board].copy()
        c[5, 1] = np.nextafter(c[5, 1], np.float32(np.inf))
        change = {"coords": c}
    elif change == "float64":
        change = {"coords": world.clouds[board].astype(np.float64)}
    p = _odd_parent(world, name, role, **change)
    _refused(world, _argv(world, name, **{role: p}), capsys, says)


def test_var_names_other_than_the_boards_genes_are_refused(world, capsys, monkeypatch):
    (world.tmp / "b2").mkdir()
    monkeypatch.setattr(hf, "BOARDS_DIR", _boards(world.tmp / "b2", genes=GENES[::-1]))
    _refused(world, _argv(world, NAMES[2]), capsys, "its var_names are not the genes of T2:embryo:test_interp")


@pytest.mark.parametrize("bad, says", [
    ("sparse", "X is not a dense array"), ("float64", "X is float64, not float32"),
    ("nan", "X holds NaN or infinite values"), ("inf", "X holds NaN or infinite values"),
    ("negative", "X holds negative values")])
@pytest.mark.parametrize("name, role", [("T3_c2_s2_mxe2.h5ad", "second"), ("T3_s2_mx_xe2_ck.h5ad", "first"),
                                        ("T2embryo_s2x_ovmxso.h5ad", "second")])
def test_a_parents_x_not_dense_float32_finite_and_non_negative_is_refused(world, capsys, bad, says, name, role):
    key = ROLES[name][{"first": 1, "second": 2}[role]]
    X = world.xs[key].copy()
    if bad == "sparse":
        X = sparse.csr_matrix(X)
    elif bad == "float64":
        X = X.astype(np.float64)
    else:
        X[3, 4] = {"nan": np.nan, "inf": np.inf, "negative": -0.5}[bad]
    p = _odd_parent(world, name, role, X=X)
    _refused(world, _argv(world, name, **{role: p}), capsys, says)


@pytest.mark.parametrize("odd, says", [
    ({"names": [f"c{i}" for i in range(N)]}, "the anchor's obs names are not cell_0"),
    ({"obs_cols": {"celltype": ["a"] * N}}, "the anchor's obs holds columns"),
    ({"more_obsm": {"X_umap": np.zeros((N, 2), dtype=np.float32)}}, "the anchor's obsm holds other keys"),
    ("float64", "is not (n, 3) float32"),
    ("four", "is not (n, 3) float32"),
    ("point", "is not finite with spread"),
])
def test_an_anchor_whose_obs_or_obsm_the_writer_would_not_give_back_is_refused(world, capsys, odd, says):
    """The anchor and both parents share the oddity (so they agree with each other), pinned for an embryo file."""
    name = NAMES[2]
    coords = world.clouds[EMB]
    if odd == "float64":
        odd = {"coords": coords.astype(np.float64)}
    elif odd == "four":
        odd = {"coords": np.hstack([coords, np.zeros((N, 1), dtype=np.float32)])}
    elif odd == "point":
        odd = {"coords": np.zeros((N, 3), dtype=np.float32)}
    paths = {}
    for role, key in (("first", "O"), ("second", "SO")):
        kw = dict(odd)
        paths[role] = world.tmp / "odd" / f"{role}.h5ad"
        _write(paths[role], world.xs[key], kw.pop("coords", coords), **kw)
    _pin(world, name, anchor=paths["first"], first=paths["first"], second=paths["second"])
    _refused(world, _argv(world, name, anchor=paths["first"], **paths), capsys, says)


# ------------------------------------------------------------------ the writers' own refusals: exit 1


def test_the_embryo_writers_refusal_exits_1_and_writes_nothing(world, capsys):
    """t2_build.build_submission refuses repeated positions before writing: exit 1, digits masked."""
    name = NAMES[3]
    coords = world.clouds[EMB].copy()
    coords[7] = coords[3]
    paths = {}
    for role, key in (("first", "O"), ("second", "SO")):
        paths[role] = world.tmp / "odd" / f"{role}.h5ad"
        _write(paths[role], world.xs[key], coords)
    _pin(world, name, anchor=paths["first"], first=paths["first"], second=paths["second"])
    assert hf.main(_argv(world, name, anchor=paths["first"], **paths)) == 1
    io = capsys.readouterr()
    assert io.out == "" and io.err.startswith("failed: the board's writer refused: ")
    assert not re.search(r"\d", io.err)
    _nothing_written(world)


def test_a_failed_t3_format_check_exits_1_without_a_sidecar(world, capsys, monkeypatch):
    (world.tmp / "b3").mkdir()
    monkeypatch.setattr(hf, "BOARDS_DIR", _boards(world.tmp / "b3", t3_min=N + 1))
    out = world.tmp / "out" / "a" / NAMES[0]
    assert hf.main(_argv(world, NAMES[0], out)) == 1
    io = capsys.readouterr()
    assert io.out == "" and io.err.startswith("failed: the format check failed: ") and not re.search(r"\d", io.err)
    assert not out.with_suffix(".json").exists()


# ------------------------------------------------------------------ the launcher: what it runs


def _sh() -> str:
    return LAUNCHER.read_text()


def _code(sh: str) -> str:
    """The launcher without its comment lines and trailing comments ("$#" stays)."""
    return "\n".join(re.sub(r"\s#(\s.*)?$", "", line) for line in sh.splitlines()
                     if not line.lstrip().startswith("#"))


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def _launcher_table(sh: str) -> dict[str, tuple[str, str, str, str]]:
    shas = dict(re.findall(r"^([A-Z0-9]+_SHA)=([0-9a-f]{64})  #", sh, flags=re.M))
    rows = re.findall(r"^  (\S+)\) BOARD=(\S+) ANCHOR_PIN=\$(\w+) FIRST_PIN=\$(\w+) SECOND_PIN=\$(\w+) ;;$", sh,
                      flags=re.M)
    return {n: (b, shas[a], shas[f], shas[s]) for n, b, a, f, s in rows}


def test_the_launchers_table_and_pins_are_the_builders():
    sh = _sh()
    assert _launcher_table(sh) == {n: (h.board, h.anchor, *h.parents) for n, h in hf.TABLE.items()}
    for var in SHA_VARS.values():
        assert f"{var}={getattr(hf, var)}  #" in sh and sh.count(getattr(hf, var)) == 1
    assert '*) echo "refused: NAME is not one of row R68\'s four files" >&2; exit 89 ;;' in sh


def test_the_launcher_runs_the_designs_command_twice_and_checks_each_file():
    sh = _sh()
    assert ('$PY experiments/builds/hedge_fixed.py build --captain-ruling R68 --name "$NAME" --anchor "$ANCHOR" \\\n'
            '    --first "$FIRST" --second "$SECOND" --out "$d/$NAME" > /dev/null 2>&1') in sh
    assert 'uv run --frozen --no-sync vec check "$d/$NAME" --board "$BOARD" > /dev/null 2>&1 || return 13' in sh
    assert 'for STEP in a b; do begin; build "$STEP" || finish $?; step_ok; done' in sh
    done = sh[sh.index("STEP=done"):]
    assert '(cd "$O" && sha256sum "a/$NAME" "b/$NAME") > "$S/SHA256SUMS" || finish 40' in done
    assert 'cmp -s "$O/a/$NAME" "$O/b/$NAME" || finish 41' in done
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh and "O=runs/hedgefix_${head:0:7}_$RUN" in sh
    assert "S=$HOME/v/scratch/hedgefix/$RUN" in sh and 'TMPDIR="$S/tmp"' in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    first_uv = sh.index("uv run")
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < first_uv
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert "upload" not in sh.replace("Nothing here uploads", "")


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = _code(_sh()).splitlines()
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 6
    for msg in usages:  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    for line in sh.splitlines():
        if re.search(r"hedge_fixed\.py build|vec check \"", line) and not line.lstrip().startswith("#"):
            assert "> /dev/null 2>&1" in line or line.rstrip().endswith("\\"), line
    code = _code(sh)
    assert code.count(".json") == 1 and '"$1.stamp.json"' in code  # the stamp only, never the sidecar
    for word in ("du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -", "tail ",
                 "cat ", "curl", "wget", "git fetch", "git pull"):
        assert word not in code, word
    assert 'print(s["commit"], s["clean"])' in sh and '"$1.stamp.json"' in sh


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    assert sh[sh.index("STEP=done\n"):].split("\n", 2)[1] == "begin"
    for code in ("exit 90", "exit 89", "exit 91"):
        assert sh.index(code) < sh.index("STEP=gate"), code
    gate = sh[sh.index("STEP=gate"):sh.index("stamp() {")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")
    assert [int(c) for c in re.findall(r"finish (\d+)", gate)] == [96, 98, 95, 95, 95, 95, 81, 82, 83, 88, 93]


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", _code(sh))}
    documented = {int(c) for c in re.findall(r"(?:^|[\s,#])(\d\d) ", doc, flags=re.M)}
    assert used == documented == {90, 89, 91, 96, 98, 95, 81, 82, 83, 88, 93, 11, 12, 13, 14, 15, 40, 41}


# ------------------------------------------------------------------ the launcher, run for real

STAND_IN = {k: f"parent {k}\n".encode() for k in SHA_VARS}  # the stand-in parents' bytes
FAILING_UV = "#!/bin/sh\nexit 1\n"
# A stand-in `uv` for runs past the gate: the builder writes "built <name>" (plus " other" in b/ with B_DIFFERS) and a
# stamp naming HEAD and clean (CLEAN=false: not clean; NO_STAMP: none), or exits 2 (REFUSES) or 1 (FAILS); `vec check`
# passes unless CHECK_FAILS; every builder and check call is listed in calls.txt; the stamp reader is python itself.
FULL_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) echo "$*" >> calls.txt; [ -n "${CHECK_FAILS:-}" ] && exit 1; exit 0 ;;
  python)
    shift 3
    case "$1" in
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/hedge_fixed.py)
        echo "$*" >> calls.txt
        [ -n "${REFUSES:-}" ] && exit 2
        [ -n "${FAILS:-}" ] && exit 1
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        mkdir -p "$(dirname "$out")"
        body="built $(basename "$out")"
        case "$out" in */b/*) [ -n "${B_DIFFERS:-}" ] && body="$body other" ;; esac
        echo "$body" > "$out"
        [ -n "${NO_STAMP:-}" ] || printf '{"commit": "%s", "clean": %s}' "$(git rev-parse HEAD)" "${CLEAN:-true}" \\
          > "$out.stamp.json"
        exit 0 ;;
    esac ;;
esac
exit 1
"""


def _stand_in_launcher() -> str:
    """The launcher with its pins set to the stand-in parents' sha256."""
    sh = _sh()
    for key, var in SHA_VARS.items():
        sh = sh.replace(f"{var}={getattr(hf, var)}", f"{var}={hashlib.sha256(STAND_IN[key]).hexdigest()}")
    assert sh != _sh()
    return sh


def _env(tmp_path: Path) -> dict:
    env = git_env()
    env.update(HOME=str(tmp_path / "home"), GIT_CEILING_DIRECTORIES=str(tmp_path),
               UV_PROJECT_ENVIRONMENT=str(tmp_path / "elsewhere"))  # the launcher unsets the last
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R68, venv: str | None = "own", uv: str = FAILING_UV,
          launcher: str | None = None) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher (its pins the stand-in parents') and the files its gate reads, the
    stand-in parents under tmp_path/parents, and a home whose ~/.local/bin holds the given `uv`."""
    env = _env(tmp_path)
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True)
    (fake / "uv").write_text(uv)
    (fake / "uv").chmod(0o755)
    repo = tmp_path / "repo"
    (repo / "experiments" / "builds").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "src").mkdir()
    (repo / "src" / "keep.txt").write_text("src\n")
    (repo / "experiments" / "builds" / "hedge_fixed_line.sh").write_text(launcher or _stand_in_launcher())
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / ".gitignore").write_text("runs/\n.venv/\n")
    if venv is not None:
        site = repo / ".venv" / "lib" / "python3.12" / "site-packages"
        site.mkdir(parents=True)
        (site / "vec.pth").write_text(f"{repo / 'src' if venv == 'own' else tmp_path / 'other' / 'src'}\n")
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    (tmp_path / "parents").mkdir()
    for key, data in STAND_IN.items():
        (tmp_path / "parents" / f"{key}.h5ad").write_bytes(data)
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _files(tmp_path: Path, name: str, **swap: str) -> list[str]:
    """NAME and its anchor, first and second stand-in paths, as the table orders them (swap: role -> another key)."""
    a, f, s, _, _ = ROLES[name]
    keys = {"anchor": a, "first": f, "second": s, **swap}
    return [name, *(str(tmp_path / "parents" / f"{keys[r]}.h5ad") for r in ("anchor", "first", "second"))]


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "experiments/builds/hedge_fixed_line.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "hedgefix" / run


def _exits(tmp_path: Path, run: str = "r1") -> dict[str, str]:
    return {p.name[:-5]: p.read_text().strip() for p in _scratch(tmp_path, run).glob("*.exit")}


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert _exits(tmp_path, run) == {"gate": str(code)}
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no head file, no later step's mark
    assert not list((s / "tmp").iterdir())


ONE = "T2embryo_s2x_ovbkso.h5ad"  # the file the gate's refusals are run on


@pytest.mark.parametrize("argv, code", [
    ([], None), (["r1"], None), (["r1", ONE], None), (["r1", ONE, "a", "f"], None), (["r1", ""], None),
    (["a/b", "N", "a", "f", "s"], 90), (["r 1", "N", "a", "f", "s"], 90), (["$(id)", "N", "a", "f", "s"], 90),
    (["r1", ONE, "a", "f", "s", "extra"], 89),
    *[(["r1", bad, "a", "f", "s"], 89) for bad in (
        "T3_c2_s2.h5ad", "T2embryo_s2x_ov.h5ad", "T2embryo_s2xs_ov.h5ad", "T3_c2_s2_mxe2", "t3_c2_s2_mxe2.h5ad",
        "T3_c2_s2_mxe2.h5ad ", "runs/T3_c2_s2_mxe2.h5ad", "$(id)", "*")],
])
def test_odd_arguments_and_unknown_names_are_refused_before_any_folder_is_made(tmp_path, argv, code):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, *argv)
    assert out.returncode != 0 if code is None else out.returncode == code, out.stderr
    assert not (tmp_path / "home" / "v").exists() and not (repo / "runs").exists()


def test_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head, "r1", *_files(tmp_path, ONE)).returncode == 11
    s = _scratch(tmp_path, "r1")
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1", *_files(tmp_path, NAMES[0]))
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    assert _run(repo, env, head, "r2", *_files(tmp_path, ONE)).returncode == 11


def test_a_wrong_commit_and_a_prefix_are_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1", *_files(tmp_path, ONE))
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()
    assert _run(repo, env, head[:12], "r2", *_files(tmp_path, ONE)).returncode == 96
    _stopped_at_the_gate(tmp_path, "r2", 96)


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R68}more\n")
    out = _run(repo, env, head, "r1", *_files(tmp_path, ONE))
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


@pytest.mark.parametrize("tries, says", [
    ("", "exactly one row R68"),
    (ROW_R68 + ROW_R68, "exactly one row R68"),
    (ROW_R68.replace("| R68 (", "| R68b ("), "exactly one row R68"),
    (ROW_R68.replace("R68 (the captain's order", "R68 (a proposal"), "not the captain's ruling"),
    (ROW_R68.replace("R68 (the captain's", "R68 awaits (the captain's"), "not the captain's ruling"),
    (ROW_R68.replace("ruling, no run", "a try"), "does not say ruling, no run"),
    (ROW_R68.replace(f"`{ONE}`", ONE), f"does not name `{ONE}`"),
    (ROW_R68.replace(f"`{ONE}`", f"`{ONE} x`"), f"does not name `{ONE}`"),
    (f"| R68 (the captain's order of X, ruling, no run: `experiments/builds/hedge_fixed.py` `--captain-ruling R68`) "
     f"| `{ONE}` |\n", f"does not name `{ONE}`"),  # named outside the first cell only
    (ROW_R68.replace("`experiments/builds/hedge_fixed.py`", "`experiments/builds/hedge_mix.py`"),
     "does not name `experiments/builds/hedge_fixed.py`"),
    (ROW_R68.replace("`--captain-ruling R68`", "`--captain-ruling R60`"), "does not name `--captain-ruling R68`"),
    (ROW_R68.replace("`--captain-ruling R68`", "`--captain-ruling R68 now`"), "does not name `--captain-ruling R68`"),
])
def test_each_row_refusal(tmp_path, tries, says):
    repo, env, head = _repo(tmp_path, tries=tries)
    out = _run(repo, env, head, "r1", *_files(tmp_path, ONE))
    assert out.returncode == 95, out.stderr
    assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 95)
    assert not (repo / "runs").exists()


@pytest.mark.parametrize("name, swap, code", [
    ("T3_c2_s2_mxe2.h5ad", {"anchor": "xe2"}, 81),
    ("T3_s2_mx_xe2_ck.h5ad", {"anchor": "none"}, 81),  # no such file
    ("T3_c2_s2_mxe2.h5ad", {"first": "ck"}, 82),
    ("T3_s2_mx_xe2_ck.h5ad", {"first": "ck", "second": "xe2"}, 82),  # the parents in the wrong order
    ("T2embryo_s2x_ovmxso.h5ad", {"first": "SO", "second": "O"}, 82),
    ("T3_s2_mx_xe2_ck.h5ad", {"second": "C2"}, 83),
    ("T2embryo_s2x_ovbkso.h5ad", {"second": "O"}, 83),
])
def test_each_pin_refusal(tmp_path, name, swap, code):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1", *_files(tmp_path, name, **swap))
    assert out.returncode == code, out.stderr
    assert "is not a file with its pinned sha256" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs").exists()


def test_the_real_pins_refuse_the_stand_in_parents(tmp_path):
    """The launcher as committed (its real pins): the stand-in parents are other bytes."""
    repo, env, head = _repo(tmp_path, launcher=_sh())
    out = _run(repo, env, head, "r1", *_files(tmp_path, NAMES[0]))
    assert out.returncode == 81
    _stopped_at_the_gate(tmp_path, "r1", 81)


@pytest.mark.parametrize("venv", [None, "other"])
def test_a_venv_of_another_checkout_or_none_is_refused(tmp_path, venv):
    repo, env, head = _repo(tmp_path, venv=venv)
    out = _run(repo, env, head, "r1", *_files(tmp_path, ONE))
    assert out.returncode == 88 and "no .venv of its own" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 88)


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"hedgefix_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1", *_files(tmp_path, ONE))
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]


@pytest.mark.parametrize("name", NAMES)
def test_each_name_passes_the_gate_and_stops_at_a_failing_build_with_nothing_built(tmp_path, name):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1", *_files(tmp_path, name))
    assert out.returncode == 11, out.stderr
    s = _scratch(tmp_path, "r1")
    assert _exits(tmp_path) == {"gate": "0", "a": "11"}
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["a.exit", "gate.exit", "head", "stamps", "tmp"]
    o = repo / "runs" / f"hedgefix_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())


@pytest.mark.parametrize("name", NAMES)
def test_a_whole_stand_in_run_passes_every_step_in_order(tmp_path, name):
    repo, env, head = _repo(tmp_path, uv=FULL_UV)
    files = _files(tmp_path, name)
    out = _run(repo, env, head, "r1", *files)
    assert out.returncode == 0, out.stderr
    assert _exits(tmp_path) == {"gate": "0", "a": "0", "b": "0", "done": "0"}
    s = _scratch(tmp_path, "r1")
    r = f"runs/hedgefix_{head[:7]}_r1"
    o = repo / r
    assert sorted(p.name for p in o.iterdir()) == ["SHA256SUMS", "a", "b"]
    sha = hashlib.sha256(f"built {name}\n".encode()).hexdigest()
    assert (o / "SHA256SUMS").read_text() == (s / "SHA256SUMS").read_text() == f"{sha}  a/{name}\n{sha}  b/{name}\n"
    assert (s / "stamps").read_text() == f"a/{name} {head} True\nb/{name} {head} True\n"
    _, anchor, first, second = files
    board = ROLES[name][3]
    built = [f"experiments/builds/hedge_fixed.py build --captain-ruling R68 --name {name} --anchor {anchor} "
             f"--first {first} --second {second} --out {r}/{d}/{name}" for d in ("a", "b")]
    checks = [f"vec check {r}/{d}/{name} --board {board}" for d in ("a", "b")]
    assert (repo / "calls.txt").read_text().splitlines() == [built[0], checks[0], built[1], checks[1]]


@pytest.mark.parametrize("flag, step, code", [
    ({"REFUSES": "1"}, "a", 12), ({"FAILS": "1"}, "a", 11), ({"CHECK_FAILS": "1"}, "a", 13),
    ({"NO_STAMP": "1"}, "a", 14), ({"CLEAN": "false"}, "a", 15), ({"B_DIFFERS": "1"}, "done", 41),
])
def test_each_failure_after_the_gate_stops_the_run_at_its_step(tmp_path, flag, step, code):
    repo, env, head = _repo(tmp_path, uv=FULL_UV)
    out = _run(repo, {**env, **flag}, head, "r1", *_files(tmp_path, ONE))
    assert out.returncode == code, out.stderr
    steps = ["gate", "a", "b", "done"]
    ran = steps[:steps.index(step) + 1]
    assert _exits(tmp_path) == {**{x: "0" for x in ran[:-1]}, step: str(code)}
