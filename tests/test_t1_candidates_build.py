"""`build_test.py t1 --g-mom-damping` and `--feng-fraction-mult`: the T1 candidates of docs/PLAN.md decision log
2026-11-20 (i), each a single change to the recomputed last file (`t1 --g-mom --base lib3 --feng-groups recomputed`).
A: g_mom's factors finished at damping 1.0 in place of 0.5; B: the outside atlas's direction at 1.5 times the calendar
fraction; C: both. Synthetic stages and directions only (the build_t1 fixture of tests/test_hidden_counts.py and the
recomputed setup of tests/test_feng_relabel.py): no released file, outside-atlas file, board or score is read."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
from test_feng_relabel import GROUPS4, TARGET, _recomputed_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)

from vec import harness  # noqa: E402
from vec import t1_state_transport as ts  # noqa: E402

RC = [*TARGET, "--feng-groups", "recomputed"]  # the recomputed last file's command, after "t1"
# sha256 of the fixture's recomputed file on all four groups, built by main at 018d7bd before these options existed
# (the same synthetic build as `_build(..., [])` below, run twice in separate processes with the same result). It pins
# today's bytes for the default command.
TODAY_RC_SHA = "a5899e2d0a4b97f8fcabca6f8781a613afdd070194c0a3c76e0b14d92249d132"
THREE = ("CM", "Endocardium", "Endothelium")  # Epicardium without a direction: its cells must never change


@pytest.fixture
def bt():
    import build_test

    return build_test


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _dense(x) -> np.ndarray:
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _spy(bt, monkeypatch) -> list[dict]:
    """Around the fixture's build_submission: per call, its direction_fraction and the population g_mom's hook gets
    and returns (obs and counts)."""
    seen = []
    inner = bt.t1_scope.build_submission

    def capture(*a, **k):
        entry = {"fraction": k.get("direction_fraction"), "groups": k.get("direction_groups"),
                 "directions": sorted(k.get("directions") or {})}
        hook = k.get("post_cells")
        if hook is not None:
            def wrapped(pop):
                entry["before"] = pop.copy()
                out = hook(pop)
                entry["after"] = out.copy()
                return out

            k = {**k, "post_cells": wrapped}
        seen.append(entry)
        return inner(*a, **k)

    monkeypatch.setattr(bt.t1_scope, "build_submission", capture)
    return seen


def _build(bt, tmp_path: Path, name: str, extra: list[str]) -> tuple[Path, dict]:
    out = tmp_path / f"{name}.h5ad"
    bt.main(["t1", *RC, *extra, "--out", str(out)])
    return out, json.loads(out.with_suffix(".build.json").read_text())


@pytest.fixture
def rc4(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch)
    return _spy(bt, monkeypatch)


@pytest.fixture
def rc3(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch, groups=THREE)
    return _spy(bt, monkeypatch)


def _same_frame(f1: Path, f2: Path) -> tuple[np.ndarray, np.ndarray]:
    """Two files with the same cells in the same order, the same obs and var: returns their expression."""
    a, b = ad.read_h5ad(f1), ad.read_h5ad(f2)
    assert list(a.obs_names) == list(b.obs_names) and a.obs.equals(b.obs)
    assert list(a.var_names) == list(b.var_names) and a.var.equals(b.var)
    assert a.obsm.keys() == b.obsm.keys() and a.uns.keys() == b.uns.keys()
    return _dense(a.X), _dense(b.X)


def _rows_changed(x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
    return (x1 != x2).any(axis=1)


# ------------------------------------------------------------------ defaults: today's bytes


def test_the_options_are_off_by_default(bt):
    a = bt.build_parser().parse_args(["t1", *RC, "--out", "x.h5ad"])
    assert a.g_mom_damping is None and a.feng_fraction_mult is None
    assert bt.G_MOM_DAMPINGS == (harness.B0_DAMPING, 1.0) == (0.5, 1.0) and bt.FENG_FRACTION_MULTS == (1.0, 1.5)
    assert bt.T1_CANDIDATE_COMMAND == "t1 --g-mom --base lib3 --feng-groups recomputed"


def test_the_recomputed_file_is_todays_with_and_without_the_default_values(bt, rc4, tmp_path):
    """No option: today's bytes and today's record (pinned against main's build of the same fixture). The default
    values given explicitly: the same bytes; the record differs only by the t1_candidate field and g_mom's damping."""
    f0, r0 = _build(bt, tmp_path, "T1_rc", [])
    assert _sha(f0) == TODAY_RC_SHA
    assert "t1_candidate" not in r0 and "damping" not in r0["g_mom"]
    assert r0["fraction"] == bt.feng_dirs.bracket_fraction(10.5, 12.5) == 0.4
    f1, r1 = _build(bt, tmp_path, "T1_rc_defaults", ["--g-mom-damping", "0.5", "--feng-fraction-mult", "1.0"])
    assert f1.read_bytes() == f0.read_bytes()
    assert r1["t1_candidate"]["g_mom_damping"] == 0.5 and r1["t1_candidate"]["feng_fraction_mult"] == 1.0
    assert r1["t1_candidate"]["applied_fraction"] == r1["t1_candidate"]["calendar_fraction"] == r1["fraction"] == 0.4
    assert r1["g_mom"]["damping"] == 0.5 and r1["g_mom"]["rule"] == r0["g_mom"]["rule"]
    strip = {k: v for k, v in r1.items() if k != "t1_candidate"}
    strip["g_mom"] = {k: v for k, v in r1["g_mom"].items() if k != "damping"}
    assert strip == r0
    assert [e["fraction"] for e in rc4] == [0.4, 0.4]


def test_without_the_options_the_candidate_step_is_none(bt):
    """build_t1 reads the options only through _t1_candidate, which returns None when neither is given: every other
    command takes today's path (its identity tests are in the files the hand-back lists)."""
    from types import SimpleNamespace

    assert bt._t1_candidate(SimpleNamespace(g_mom_damping=None, feng_fraction_mult=None)) is None
    assert bt._t1_candidate(SimpleNamespace()) is None


# ------------------------------------------------------------------ each option alone changes expression only


def test_damping_changes_only_cells_g_mom_edits(bt, rc3, tmp_path):
    f0, _ = _build(bt, tmp_path, "T1_rc", [])
    fa, ra = _build(bt, tmp_path, "T1_rcA", ["--g-mom-damping", "1.0"])
    x0, xa = _same_frame(f0, fa)
    changed = _rows_changed(x0, xa)
    assert changed.any()
    pop = rc3[-1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    cats = pop.obs["category"].astype(str).to_numpy()
    held = bt._labels_without_direction(np.unique(labels), list(THREE))
    editable = np.isin(cats, bt.G_MOM_EDIT_CATEGORIES) & ~np.isin(labels, sorted(held))
    assert not (changed & ~editable).any()  # only cells g_mom edits change
    # the direction step is the default's: the hook gets the same population in both builds
    b0, ba = rc3[0]["before"], rc3[-1]["before"]
    assert list(b0.obs_names) == list(ba.obs_names) and b0.obs.equals(ba.obs)
    assert (_dense(b0.layers["counts"]) == _dense(ba.layers["counts"])).all()
    assert rc3[0]["fraction"] == rc3[-1]["fraction"] == 0.4
    assert ra["g_mom"]["damping"] == 1.0 and "damped by 1 in place of B0_DAMPING" in ra["g_mom"]["rule"]
    assert ra["t1_candidate"]["g_mom_damping"] == 1.0 and ra["t1_candidate"]["feng_fraction_mult"] == 1.0
    assert ra["fraction"] == 0.4


def test_fraction_changes_only_cells_with_a_direction(bt, rc3, tmp_path):
    f0, _ = _build(bt, tmp_path, "T1_rc", [])
    fb, rb = _build(bt, tmp_path, "T1_rcB", ["--feng-fraction-mult", "1.5"])
    x0, xb = _same_frame(f0, fb)
    changed = _rows_changed(x0, xb)
    assert changed.any()
    pop = rc3[-1]["before"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    groups = rc3[-1]["groups"]
    with_dir = np.array([groups.get(lab) in THREE for lab in labels])
    assert not with_dir.all() and not (changed & ~with_dir).any()  # Epicardium and the rest keep their counts
    assert rc3[0]["directions"] == rc3[-1]["directions"] == sorted(THREE)
    # the cells, their order and obs are the default's at every step (the draw and the streams are unchanged)
    b0 = rc3[0]["before"]
    assert list(b0.obs_names) == list(pop.obs_names) and b0.obs.equals(pop.obs)
    assert rb["t1_candidate"]["feng_fraction_mult"] == 1.5 and rb["t1_candidate"]["g_mom_damping"] == 0.5
    assert "damping" not in rb["g_mom"] or rb["g_mom"]["damping"] == 0.5


def test_the_multiplier_gives_the_direction_step_0_6_where_the_default_gives_0_4(bt, rc4, tmp_path):
    _, r0 = _build(bt, tmp_path, "T1_rc", [])
    _, rb = _build(bt, tmp_path, "T1_rcB", ["--feng-fraction-mult", "1.5"])
    f0, fb = rc4[0]["fraction"], rc4[-1]["fraction"]
    assert f0 == 0.4 and fb == pytest.approx(0.6, abs=1e-12)
    assert fb == 1.5 * bt.feng_dirs.bracket_fraction(10.5, 12.5)  # MULTIPLIER x fraction, as double_step_read forms it
    assert rb["fraction"] == fb and rb["t1_candidate"]["applied_fraction"] == fb
    assert rb["t1_candidate"]["calendar_fraction"] == 0.4 and r0["fraction"] == 0.4


def test_candidate_c_is_both_changes_on_the_same_cells(bt, rc4, tmp_path):
    f0, _ = _build(bt, tmp_path, "T1_rc", [])
    fc, rc = _build(bt, tmp_path, "T1_rcC", ["--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5"])
    x0, xc = _same_frame(f0, fc)
    assert _rows_changed(x0, xc).any()
    assert rc4[-1]["fraction"] == pytest.approx(0.6, abs=1e-12)
    assert rc["t1_candidate"]["g_mom_damping"] == 1.0 and rc["t1_candidate"]["feng_fraction_mult"] == 1.5
    assert rc["g_mom"]["damping"] == 1.0


# ------------------------------------------------------------------ the practice arms' constructions


def test_the_damping_hook_equals_the_practice_g_mom_at_damping_1(bt, t1_env, rc4, tmp_path):  # noqa: F811
    """The builder's hook at --g-mom-damping 1.0 on the cells it receives equals recomputed_read.rc_hook (the T62/T63
    practice hook under the recomputed groups) on the factors bets_read.gmom_fits finishes at GMOM_STRONG (arm
    gmom_s1), fitted on the same two stages; and gmom_fits' 'file' factors are the default build's."""
    import bets_read as br
    import recomputed_read as rcr

    _build(bt, tmp_path, "T1_rcA", ["--g-mom-damping", "1.0"])
    entry = rc4[-1]
    pop_in, got = entry["before"], entry["after"]
    assert br.GMOM_STRONG == 1.0 and br.GMOM_STREAM == 5
    fits = br.gmom_fits(t1_env.ref.copy(), "E10.5_RNA", t1_env.start.copy(), "E9.5_RNA")
    labels = pop_in.obs["celltype"].astype(str).unique()
    _, covered, with_direction = rcr.direction_kit(labels, {g: None for g in GROUPS4})
    held = bt.t1_scope.held_flat_genes()  # the hook's own list (the fixture's stand-in for the pre-registered one)
    want = rcr.rc_hook(fits["strong"], harness.BUILD_SEED_OFFSET + 1, held, with_direction, covered)(pop_in.copy())
    assert list(want.obs_names) == list(got.obs_names)
    assert (_dense(want.layers["counts"]) == _dense(got.layers["counts"])).all()
    assert (_dense(want.X) == _dense(got.X)).all()
    # and the factors themselves: finish_momentum of the same raw fit at 1.0
    _, raw = harness.g_mom_pair_factors(t1_env.ref.copy(), "E10.5_RNA", t1_env.start.copy(), "E9.5_RNA")
    strong = ts.finish_momentum(raw, 1.0)
    for c in strong.categories:
        assert np.array_equal(strong.categories[c].logf, fits["strong"].categories[c].logf)


# ------------------------------------------------------------------ refusals, before anything is read


def _no_read(*a, **k):
    pytest.fail("nothing is read before the refusal")


@pytest.mark.parametrize("extra, msg", [
    (["--g-mom-damping", "0.75"], "--g-mom-damping takes 0.5"),
    (["--g-mom-damping", "0"], "--g-mom-damping takes 0.5"),
    (["--g-mom-damping", "2"], "--g-mom-damping takes 0.5"),
    (["--feng-fraction-mult", "2.0"], "--feng-fraction-mult takes 1"),
    (["--feng-fraction-mult", "0.5"], "--feng-fraction-mult takes 1"),
])
def test_other_values_are_refused(bt, monkeypatch, tmp_path, extra, msg):
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    with pytest.raises(SystemExit, match=msg):
        bt.main(["t1", *RC, *extra, "--out", str(tmp_path / "T1.h5ad")])


@pytest.mark.parametrize("argv", [
    ["t1"],  # the authors' groups, today's s1 design without g_mom
    ["t1", *TARGET],  # the first upload's command
    ["t1", *TARGET, "--feng-groups", "none"],  # the none file
    ["t1", "--g-mom", "--base", "mid", "--cm-share", "0.60"],  # T1's second upload
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "recomputed", "--seed", "1000002"],  # another seed
    ["t1", "--g-mom", "--base", "lib3", "--feng-groups", "recomputed", "--board", "T1:val"],  # a dry run
])
@pytest.mark.parametrize("opt", [["--g-mom-damping", "1.0"], ["--feng-fraction-mult", "1.5"],
                                 ["--g-mom-damping", "0.5"], ["--feng-fraction-mult", "1.0"]])
def test_the_options_are_refused_on_any_other_command(bt, monkeypatch, tmp_path, argv, opt):
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match="change T1's recomputed last file only"):
        bt.main([*argv, *opt, "--out", str(tmp_path / name)])


@pytest.mark.parametrize("extra", [["--cm-share", "0.6"], ["--herd"], ["--no-feng-hygiene"], ["--strong-null"],
                                   ["--lib3-pair"], ["--blood-share", "0.1"]])
def test_the_recomputed_commands_own_guards_still_hold_with_the_options(bt, monkeypatch, tmp_path, extra):
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    with pytest.raises(SystemExit, match="--feng-groups recomputed builds one file"):
        bt.main(["t1", *RC, *extra, "--g-mom-damping", "1.0", "--feng-fraction-mult", "1.5", "--out",
                 str(tmp_path / "T1.h5ad")])


def test_without_g_mom_the_recomputed_guard_refuses_the_options(bt, monkeypatch, tmp_path):
    monkeypatch.setattr(bt, "load_stage", _no_read)
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _no_read)
    with pytest.raises(SystemExit, match="--feng-groups recomputed builds one file"):
        bt.main(["t1", "--base", "lib3", "--feng-groups", "recomputed", "--g-mom-damping", "1.0", "--out",
                 str(tmp_path / "T1.h5ad")])
