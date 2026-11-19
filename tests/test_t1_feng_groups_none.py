"""`build_test.py t1 --feng-groups none` and `feng_dirs.feng_directions(groups="none")`: the entry with Feng's direction
taken out and nothing put in its place (docs/PLAN.md decision log 2026-11-19 (a) to (d)). Synthetic arrays and
hand-made files only: no Feng file, released file, board or score is read."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/builds"):
    sys.path.insert(0, str(ROOT / _p))

import feng_dirs  # noqa: E402
import filterio as fio  # noqa: E402
from test_g_mom_verify import _g_mom_loads  # noqa: E402
from test_hidden_counts import T1_HELD, _reference_file, t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)

REAL_FENG_DIRECTIONS = feng_dirs.feng_directions  # taken at import, before any fixture replaces it on the module
GROUPS4 = ("CM", "Endocardium", "Endothelium", "Epicardium")
NONE_RECORD = {"file": None, "sha256": None, "bracket": ["E9.5", "E14.5"], "groups": {}, "feng_groups": "none"}


@pytest.fixture
def bt():
    import build_test

    return build_test


def _sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _real_none(bt, env, monkeypatch) -> None:
    """The fixture answers every call for Feng's directions with its synthetic ones. Here a call with groups="none"
    goes to the real function instead (and is logged), so the option's own path is the one under test."""
    synthetic = bt.feng_dirs.feng_directions

    def directions(genes, *a, **k):
        if k.get("groups") == "none":
            env.reads.append(("none", list(genes)))
            return REAL_FENG_DIRECTIONS(genes, *a, **k)
        return synthetic(genes, *a, **k)

    monkeypatch.setattr(bt.feng_dirs, "feng_directions", directions)


def _fail(what: str):
    return lambda *a, **k: pytest.fail(what)


def _kinds(reads) -> list[str]:
    """The fixture's log, by kind: "E10.5" for a load of the reference, "Feng" for a synthetic answer, "none"."""
    return [r if isinstance(r, str) else r[0] for r in reads]


# feng_dirs.feng_directions ---------------------------------------------------------------------------------------


def test_the_groups_are_authors_then_none():
    assert feng_dirs.FENG_GROUPS == ("authors", "none")


def test_none_gives_no_direction_and_reads_nothing(monkeypatch):
    for name in ("load_config", "clean_stats_path", "directions_from_stats"):
        monkeypatch.setattr(feng_dirs, name, _fail(f"groups='none' must not call {name}"))
    monkeypatch.setattr(feng_dirs.fio, "verify_manifest_output", _fail("groups='none' reads no manifest"))
    monkeypatch.setattr(feng_dirs.dl, "load_stats", _fail("groups='none' reads no Feng statistics"))
    d, rec = feng_dirs.feng_directions(["g1", "g2"], groups="none")
    assert d == {} and rec == NONE_RECORD
    d, rec = feng_dirs.feng_directions(["g1"], bracket=("E14.5", "E18"), groups="none")
    assert d == {} and rec["bracket"] == ["E14.5", "E18"] and rec["file"] is None


def test_an_unknown_groups_value_is_refused_before_anything_is_read(monkeypatch):
    monkeypatch.setattr(feng_dirs, "load_config", _fail("refused before the config is read"))
    for bad in ("recomputed", "jax", "", "None"):
        with pytest.raises(ValueError, match="groups must be one of"):
            feng_dirs.feng_directions(["g1"], groups=bad)


def test_the_authors_default_is_the_call_as_before(tmp_path, monkeypatch):
    """No groups argument and groups='authors' are the same call on the same file, and the record carries no new key."""
    monkeypatch.setattr(feng_dirs, "REPO", tmp_path)
    genes = [f"f{i}" for i in range(30)]
    stats = feng_dirs.clean_stats_path()
    rng = np.random.default_rng(3)
    keys, n, cp = [], [], []
    for lab in ("ventricular_cm", "atrial_cm", "endocardial_ec", "vascular_ec", "epicardial"):
        for stage in ("E9.5", "E14.5"):
            for half in "AB":
                keys.append(f"{stage}|LV|{lab}|{half}")
                n.append(20)
                cp.append(20 * rng.gamma(0.6, 20.0, len(genes)))
    stats.parent.mkdir(parents=True, exist_ok=True)
    n = np.array(n, dtype=float)
    np.savez(stats, genes=np.array(genes), keys=np.array(keys), n=n, S_cp=np.array(cp),
             S_log=np.zeros((len(keys), 30)), S_nz=np.zeros((len(keys), 30)), depth_sum=n)
    fio.write_manifest(stats.parent, {"outputs": {stats.name: {"sha256": fio.sha256_file(stats)}}})
    d0, r0 = feng_dirs.feng_directions(genes)
    d1, r1 = feng_dirs.feng_directions(genes, groups="authors")
    assert r0 == r1 and "feng_groups" not in r0 and set(d0) == set(GROUPS4)
    assert all(np.array_equal(d0[g], d1[g], equal_nan=True) for g in d0)


# build_test.py t1 --feng-groups ----------------------------------------------------------------------------------


def test_the_option_is_off_by_default_and_takes_two_values(bt):
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "x.h5ad"]).feng_groups == "authors"
    assert ap.parse_args(["t1", "--feng-groups", "none", "--out", "x.h5ad"]).feng_groups == "none"
    for bad in ("recomputed", "jax"):
        with pytest.raises(SystemExit):
            ap.parse_args(["t1", "--feng-groups", bad, "--out", "x.h5ad"])
    assert bt.FENG_GROUPS == ("authors", "none") and bt.FENG_GROUPS_DEFAULT == "authors"


def test_the_default_makes_the_old_call_with_no_keyword(bt, monkeypatch):
    calls = []
    monkeypatch.setattr(bt.feng_dirs, "feng_directions",
                        lambda genes, *a, **k: calls.append((genes, a, k)) or ({"CM": np.zeros(1)}, {"file": "f"}))
    bt._feng_directions(["g"])
    bt._feng_group_names()
    bt._feng_directions(["g"], "none")
    assert calls == [(["g"], (), {}), ([""], (), {}), (["g"], (), {"groups": "none"})]


def test_the_default_build_and_its_record_are_as_before(bt, t1_env, tmp_path):  # noqa: F811
    """With the option not given, and with it given as authors: the same bytes as the entry's own builder writes, the
    same reads in the same order, and a record with no feng_groups field."""
    today = _reference_file(bt, t1_env, tmp_path / "today.h5ad")
    out, out2 = tmp_path / "T1.h5ad", tmp_path / "T1_authors.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--out", str(out)])
    reads = list(t1_env.reads)
    bt.main(["t1", "--no-feng-hygiene", "--feng-groups", "authors", "--out", str(out2)])
    assert _sha(out) == _sha(out2) == _sha(today)
    assert t1_env.reads == reads + reads and "Feng" in _kinds(reads) and "none" not in _kinds(reads)
    for o in (out, out2):
        assert "feng_groups" not in json.loads(o.with_suffix(".build.json").read_text())


def test_none_is_the_entry_with_no_direction_and_reads_no_outside_file(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """The file equals the one the entry's builder writes when handed no direction at all (R1 on every label), with
    the pre-registered held list alone; R2's union and the Blood guard's call for Feng's group names are never made."""
    _real_none(bt, t1_env, monkeypatch)
    monkeypatch.setattr(bt.hygiene, "load_union", _fail("R2's union is not read under --feng-groups none"))
    monkeypatch.setattr(bt, "_feng_group_names", _fail("Feng's group names are not asked for under none"))
    out = tmp_path / "T1_none.h5ad"
    bt.main(["t1", "--feng-groups", "none", "--out", str(out)])
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["feng_groups"] == {"groups": "none", "file": None, "sha256": None}
    assert rec["feng_hygiene"] is None and rec["hold_unbracketed"] is True and rec["external"] == NONE_RECORD
    assert t1_env.held_calls == [list(T1_HELD)]
    assert _kinds(t1_env.reads) == ["E10.5", "none"]  # one E10.5 load, one call, no synthetic Feng answer
    t1_env.dirs = {}
    assert _sha(out) == _sha(_reference_file(bt, t1_env, tmp_path / "no_direction.h5ad"))


def test_none_draws_the_entrys_cells_and_changes_only_the_cells_feng_moved(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """Against the entry on the same inputs and seed: the same cells in the same order; a cell whose label has no Feng
    group is the same in both files; and every cell of the none file carries the expression of an E10.5 cell as it
    was measured (no cell is moved)."""
    _real_none(bt, t1_env, monkeypatch)
    entry, none = tmp_path / "T1.h5ad", tmp_path / "T1_none.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--out", str(entry)])
    bt.main(["t1", "--feng-groups", "none", "--out", str(none)])
    a, b = ad.read_h5ad(entry), ad.read_h5ad(none)
    assert list(a.obs_names) == list(b.obs_names) and a.n_obs > 0
    xa, xb = _dense(a.X), _dense(b.X)
    groups = bt.feng_dirs.label_groups(t1_env.ref.obs["celltype"].astype(str).unique())
    moved = {lab for lab, g in groups.items() if g in t1_env.dirs}
    label = _labels(b, t1_env)
    held = np.array([lab not in moved for lab in label])
    assert held.any() and (~held).any()
    assert np.array_equal(xa[held], xb[held])  # the labels R1 holds: untouched by the option
    assert not np.array_equal(xa[~held], xb[~held])  # the entry moved these along the (synthetic) direction
    measured = {r.tobytes() for r in np.round(_dense(t1_env.ref.X), 4).astype(np.float32)}
    assert all(r.tobytes() in measured for r in np.round(xb, 4).astype(np.float32))


def _dense(x) -> np.ndarray:
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _labels(sub, env) -> list[str]:
    """Each submitted cell's E10.5 label, by matching its expression row to the reference's (the none file moves no
    cell, so every row is a measured one)."""
    ref = np.round(_dense(env.ref.X), 4).astype(np.float32)
    by_row = {r.tobytes(): lab for r, lab in zip(ref, env.ref.obs["celltype"].astype(str), strict=True)}
    return [by_row[r.tobytes()] for r in np.round(_dense(sub.X), 4).astype(np.float32)]


def test_none_with_g_mom_edits_the_labels_of_r1s_four_groups(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    _real_none(bt, t1_env, monkeypatch)
    _g_mom_loads(bt, t1_env, monkeypatch)
    out = tmp_path / "T1_none_gmom.h5ad"
    bt.main(["t1", "--g-mom", "--feng-groups", "none", "--out", str(out)])
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["g_mom"]["categories_with_factors"] and rec["feng_groups"]["groups"] == "none"
    groups = {"V-CM": "CM", "Endocardium-1": "Endocardium", "BEC": "Endothelium", "Blood": "Blood",
              "Epicardium": "Epicardium", "Macrophage": "Macrophage"}
    assert bt._g_mom_covered(groups, {}, "none") == {"V-CM", "Endocardium-1", "BEC", "Epicardium"}
    assert bt._g_mom_covered(groups, {"CM": 1}) == {"V-CM"}  # the default: the labels Feng's direction moves
    assert bt._g_mom_covered(groups, {"CM": 1}, "authors") == {"V-CM"}


def test_none_is_refused_with_the_options_that_act_on_fengs_direction(bt, monkeypatch, tmp_path):
    monkeypatch.setattr(bt, "load_stage", _fail("nothing is read before the refusal"))
    monkeypatch.setattr(bt.hygiene, "load_union", _fail("nothing is read before the refusal"))
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _fail("nothing is read before the refusal"))
    ap, out = bt.build_parser(), str(tmp_path / "x.h5ad")
    for flag, msg in (("--no-hold-unbracketed", "refused with --no-hold-unbracketed"),
                      ("--feng-hygiene", "refused with --feng-hygiene"),
                      ("--strong-null", "refused with --strong-null")):
        with pytest.raises(SystemExit, match=f"--feng-groups none is {msg}"):
            bt.build_t1(ap.parse_args(["t1", flag, "--feng-groups", "none", "--out", out]))
