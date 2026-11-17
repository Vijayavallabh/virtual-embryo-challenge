"""experiments/builds/build_test.py t1 --direction-form (campaign branch campaign/t1-w; no ruling): the default path
is main's builder byte for byte, and `common` hands build_submission try T55's form cn
(experiments/t1_m1/part_read.py) and changes nothing else. Synthetic stages only (the fixtures of
tests/test_hidden_counts.py): no released file, board file, Feng statistic or held-flat list is read."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from test_hidden_counts import T1_GENES, T1_HELD, T1_SECRET, UNION, bt, t1_env  # noqa: F401

from vec.boards import REPO_ROOT

# main at the time the option was built: the builder the default path must equal (pinned, so the comparison does not
# turn into a comparison with itself once the option is merged)
BASE_COMMIT = "c0c3971f95df5b15e0b89a328ef45941969ae6e1"
S1 = ["t1", "--g-mom", "--base", "lib3"]  # experiments/builds/build_all.sh s1
S2 = ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"]  # experiments/builds/build_all.sh s2
COMMON = ["--direction-form", "common"]


@pytest.fixture
def old_bt(bt, tmp_path):  # noqa: F811
    """build_test.py as committed at BASE_COMMIT, loaded beside the current one (its experiment modules are the ones
    already imported, so the fixtures' patches on them reach it too)."""
    try:
        src = subprocess.run(["git", "show", f"{BASE_COMMIT}:experiments/builds/build_test.py"], cwd=REPO_ROOT,
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as e:
        pytest.skip(f"commit {BASE_COMMIT} not readable: {e}")
    path = tmp_path / "build_test_base.py"
    path.write_text(src)
    spec = importlib.util.spec_from_file_location("build_test_base", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    yield mod
    sys.modules.pop(spec.name, None)


@pytest.fixture
def world(bt, t1_env, monkeypatch):  # noqa: F811
    """t1_env with what s1 and s2 need: library-numbered cell names, the E9.5 file for --g-mom, R2's union, and
    directions that share a part (so that sum(c d_g) > 0 in every group); `calls` logs each build_submission call."""
    ref = t1_env.ref
    lab = ref.obs["celltype"].astype(str).to_numpy()
    lib = np.where(lab == "V-CM", np.where(np.arange(ref.n_obs) % 4 == 0, 1, 3), np.arange(ref.n_obs) % 3 + 1)
    ref.obs_names = [f"bc{i}-1_{k}" for i, k in enumerate(lib)]
    rng = np.random.default_rng(11)
    shared = rng.normal(0.0, 1.0, len(T1_GENES))
    t1_env.dirs = {g: shared + rng.normal(0.0, 0.6, len(T1_GENES)) for g in sorted(t1_env.dirs)}
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: UNION)

    def load(path, **kw):
        return (t1_env.ref if Path(path).name.startswith("E10.5") else t1_env.start).copy()

    monkeypatch.setattr(bt, "load_stage", load)
    t1_env.calls = []
    real = t1_env.build

    def spy(*a, **k):
        seen = {}
        hook = k.get("post_cells")

        def wrapped(pop):  # the drawn and moved population, with its labels, as g_mom's hook receives it
            seen["labels"] = pop.obs["celltype"].astype(str).to_numpy()
            return hook(pop)

        sub, rep = real(*a, **({**k, "post_cells": wrapped} if hook else k))
        t1_env.calls.append((a, k, sub, seen.get("labels")))
        return sub, rep

    monkeypatch.setattr(bt.t1_scope, "build_submission", spy)
    return t1_env


def _patch_old(old, bt, monkeypatch):  # noqa: F811
    """The module-level names the fixtures patch on the current builder, patched the same way on the old one."""
    for name in ("load_stage", "_t1_momentum", "load_boards"):
        monkeypatch.setattr(old, name, getattr(bt, name))


def _same(x, y, what):
    """Equal by value and type: AnnData by names, obs, X and counts; arrays exactly; callables by qualified name
    (each build makes its own g_mom hook)."""
    assert type(x) is type(y), what
    if isinstance(x, ad.AnnData):
        assert list(x.obs_names) == list(y.obs_names) and list(x.var_names) == list(y.var_names), what
        assert x.obs.equals(y.obs), what
        assert np.array_equal(np.asarray(x.X.todense() if hasattr(x.X, "todense") else x.X),
                              np.asarray(y.X.todense() if hasattr(y.X, "todense") else y.X)), what
    elif isinstance(x, np.ndarray):
        assert x.dtype == y.dtype and np.array_equal(x, y, equal_nan=True), what
    elif isinstance(x, dict):
        assert list(x) == list(y), what
        for k in x:
            _same(x[k], y[k], f"{what}[{k!r}]")
    elif isinstance(x, list | tuple):
        assert len(x) == len(y), what
        for i, (a, b) in enumerate(zip(x, y, strict=True)):
            _same(a, b, f"{what}[{i}]")
    elif callable(x) and not isinstance(x, type):
        assert x.__qualname__ == y.__qualname__, what
    else:
        assert x == y, what


def _cn(directions, genes, held):
    """try T55's cn directions and its record, by part_read's own function."""
    import part_read

    forms, rec = part_read.all_directions(directions, list(genes), held)
    return forms["cn"], rec


def _sha(v) -> str:
    import zone_read

    return zone_read.direction_sha256(v)


# ----------------------------------------------------------------------------------------------- the default path
def test_the_option_is_entry_by_default_and_t1_only(bt):  # noqa: F811
    ap = bt.build_parser()
    assert ap.parse_args([*S1, "--out", "x.h5ad"]).direction_form == "entry"
    assert ap.parse_args([*S2, *COMMON, "--out", "x.h5ad"]).direction_form == "common"
    with pytest.raises(SystemExit):
        ap.parse_args(["heart", "--entry", "s1", *COMMON, "--out", "x.h5ad"])


@pytest.mark.parametrize("argv", [S1, S2], ids=["s1", "s2"])
@pytest.mark.parametrize("explicit", [False, True], ids=["absent", "entry"])
def test_default_path_is_mains_builder_byte_for_byte(bt, old_bt, world, tmp_path, monkeypatch, capsys,  # noqa: F811
                                                    argv, explicit):
    """s1's and s2's commands, with the option absent or `entry`: the call to build_submission, the written file and
    the record are main's, the record with direction_form "entry" added and nothing else; part_read is not imported."""
    _patch_old(old_bt, bt, monkeypatch)
    was_imported = "part_read" in sys.modules
    old_out, new_out = tmp_path / "old" / "T1.h5ad", tmp_path / "new" / "T1.h5ad"
    old_bt.main([*argv, "--out", str(old_out)])
    bt.main([*argv, *(["--direction-form", "entry"] if explicit else []), "--out", str(new_out)])
    assert len(world.calls) == 2
    (a0, k0, _, _), (a1, k1, _, _) = world.calls
    _same(a0, a1, "positional")
    _same(k0, k1, "keywords")
    assert {g: _sha(v) for g, v in k1["directions"].items()} == {g: _sha(v) for g, v in world.dirs.items()}
    assert new_out.read_bytes() == old_out.read_bytes()
    rec_old = json.loads(old_out.with_suffix(".build.json").read_text())
    rec_new = json.loads(new_out.with_suffix(".build.json").read_text())
    assert rec_new.pop("direction_form") == "entry" and "direction_form_common" not in rec_new
    assert rec_new == rec_old
    if not was_imported:
        assert "part_read" not in sys.modules  # the default path imports nothing new
    assert str(T1_SECRET) not in capsys.readouterr().out


# ------------------------------------------------------------------------------------------------------ common
@pytest.mark.parametrize("argv", [S1, S2], ids=["s1", "s2"])
def test_common_hands_over_t55s_cn_directions_and_nothing_else_changes(bt, world, tmp_path, capsys, argv):  # noqa: F811
    bt.main([*argv, "--out", str(tmp_path / "entry.h5ad")])
    bt.main([*argv, *COMMON, "--out", str(tmp_path / "common.h5ad")])
    (a0, k0, _, _), (a1, k1, _, _) = world.calls
    # the directions: try T55's cn made from the same statistics and held list, sha256 per group as part_read gives it
    held = sorted(set(T1_HELD) | set(UNION["genes"]))
    assert k0["held_flat"] == held
    want, rec = _cn(world.dirs, T1_GENES, held)
    got = {g: _sha(v) for g, v in k1["directions"].items()}
    assert got == {g: _sha(v) for g, v in want.items()} == {g: rec["cn"][g]["sha256"] for g in rec["cn"]}
    assert sorted(got) == sorted(world.dirs) and all(got[g] != _sha(world.dirs[g]) for g in got)
    for g in rec["cn"]:  # cn is what it claims: the entry's squared size, an exact multiple of c on M
        assert rec["cn"][g]["size_ratio"] == pytest.approx(1.0, abs=1e-12)
    # every other argument equal by value and type
    _same(a0, a1, "positional")
    _same({k: v for k, v in k0.items() if k != "directions"}, {k: v for k, v in k1.items() if k != "directions"},
          "keywords")
    # the record: direction_form, per group a_n and the sha256 as T55 records them; nothing else differs
    r0 = json.loads((tmp_path / "entry.build.json").read_text())
    r1 = json.loads((tmp_path / "common.build.json").read_text())
    assert r0["direction_form"] == "entry" and r1["direction_form"] == "common"
    c = r1.pop("direction_form_common")
    assert c["a_n"] == {g: rec["entry"][g]["a_n"] for g in sorted(world.dirs)}
    assert c["sha256"] == {g: rec["cn"][g]["sha256"] for g in sorted(world.dirs)}
    r0.pop("direction_form"), r1.pop("direction_form")
    # g_mom's median depth scale is measured on the cells it edits, which the direction moved: an outcome, not a setting
    assert r0["g_mom"].pop("median_depth_scale") > 0 and r1["g_mom"].pop("median_depth_scale") > 0
    assert r0 == r1
    text = (tmp_path / "common.build.json").read_text() + capsys.readouterr().out
    assert str(T1_SECRET) not in text and r1["n_cells"] == "max_cells"


@pytest.mark.parametrize("argv", [S1, S2], ids=["s1", "s2"])
def test_common_end_to_end_same_cells_in_order_differences_only_in_covered_cells(bt, world, tmp_path, argv):  # noqa: F811
    bt.main([*argv, "--out", str(tmp_path / "entry.h5ad")])
    bt.main([*argv, *COMMON, "--out", str(tmp_path / "common.h5ad")])
    (_, k0, sub0, lab0), (_, _, sub1, lab1) = world.calls
    f0, f1 = ad.read_h5ad(tmp_path / "entry.h5ad"), ad.read_h5ad(tmp_path / "common.h5ad")
    assert f0.shape == f1.shape and list(f0.obs_names) == list(f1.obs_names)
    assert list(f0.var_names) == list(f1.var_names) and f0.obs.equals(f1.obs)
    assert list(sub0.obs_names) == list(f0.obs_names) and sub0.obs.equals(sub1.obs)
    # the drawn cells' labels as g_mom's hook (the last edit) receives them: the same cells, in the same order
    assert len(lab0) == f0.n_obs and np.array_equal(lab0, lab1)
    labels = lab0
    groups, dirs = k0["direction_groups"], k0["directions"]
    covered = np.array([groups.get(lab) in dirs for lab in labels])
    assert covered.any() and (~covered).any()
    x0, x1 = np.asarray(f0.X), np.asarray(f1.X)
    assert np.array_equal(x0[~covered], x1[~covered])  # every cell Feng's direction does not cover is unchanged
    assert not np.array_equal(x0[covered], x1[covered])


# ------------------------------------------------------------------------------------------------------ guards
@pytest.mark.parametrize("extra, msg", [
    (["t1", "--base", "lib3"], "refused without --g-mom"),
    (["t1", "--g-mom", "--strong-null"], "refused with --strong-null"),
    (["t1", "--g-mom", "--herd"], "refused with --herd"),
    (["t1", "--g-mom", "--no-hold-unbracketed", "--no-feng-hygiene"], "refused with --no-hold-unbracketed"),
    (["t1", "--g-mom", "--base", "lib3", "--cm-share", "0.6", "--lib3-pair"], "refused with --lib3-pair"),
    (["t1", "--g-mom", "--base", "lib3", "--no-feng-hygiene"], "refused with --no-feng-hygiene"),
])
def test_each_flag_guard_refuses_before_anything_is_read(bt, t1_env, tmp_path, monkeypatch, extra, msg):  # noqa: F811
    def no_read(*a, **k):
        pytest.fail("nothing is read before the refusal")

    for name in ("load_stage", "load_boards", "_t1_momentum"):
        monkeypatch.setattr(bt, name, no_read)
    monkeypatch.setattr(bt.hygiene, "load_union", no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", no_read)
    out = tmp_path / "T1.h5ad"
    with pytest.raises(SystemExit, match=f"--direction-form common is {msg}"):
        bt.main([*extra, *COMMON, "--out", str(out)])
    assert not out.exists() and not out.with_suffix(".build.json").exists()


def test_refused_where_a_groups_moved_genes_are_not_m(bt, world, tmp_path):  # noqa: F811
    g = sorted(world.dirs)[0]
    world.dirs[g] = world.dirs[g].copy()
    world.dirs[g][3] = np.nan  # g3 is moved in the other groups only
    out = tmp_path / "T1.h5ad"
    with pytest.raises(SystemExit, match=r"--direction-form common is refused: T55 is refused: the moved genes of"):
        bt.main([*S1, *COMMON, "--out", str(out)])
    assert world.calls == [] and not out.exists() and not out.with_suffix(".build.json").exists()


def test_refused_where_sum_c_d_is_not_positive_in_a_group(bt, world, tmp_path):  # noqa: F811
    groups = sorted(world.dirs)
    rng = np.random.default_rng(5)
    base = rng.normal(0.0, 1.0, len(T1_GENES))
    world.dirs = {g: base.copy() for g in groups}
    world.dirs[groups[0]] = -0.5 * base  # c = 0.625 base, opposite to this group's d_g only
    out = tmp_path / "T1.h5ad"
    with pytest.raises(SystemExit, match=r"sum\(c d_g\) is not positive in \['" + groups[0]):
        bt.main([*S2, *COMMON, "--out", str(out)])
    assert world.calls == [] and not out.exists() and not out.with_suffix(".build.json").exists()


def test_the_determinism_check_hashes_the_files_and_not_the_records():
    """The record is written after the file, and the freeze-morning check compares the reference list's files only:
    the list names .h5ad files, so an added record key cannot fail it."""
    names = [line.split()[1] for line in
             (REPO_ROOT / "experiments/builds/results/builds_dc7b6ae.SHA256SUMS").read_text().splitlines() if line]
    assert names and all(n.endswith(".h5ad") for n in names)
    assert 'sha256sum -c "$sums"' in (REPO_ROOT / "experiments/builds/determinism_check.sh").read_text()
    text = (REPO_ROOT / "experiments/builds/build_all.sh").read_text()
    assert "--direction-form" not in text  # build_all.sh is untouched: no upload carries the option
