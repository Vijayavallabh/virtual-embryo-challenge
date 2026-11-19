"""`build_test.py t1 --feng-groups none` and `feng_dirs.feng_directions(groups="none")`: the first upload's command
(t1 --g-mom --base lib3) with Feng's direction taken out and nothing put in its place (docs/PLAN.md decision log
2026-11-19 (a) to (d)). Synthetic arrays and hand-made files only: no Feng file, released file, board or score is
read."""

from __future__ import annotations

import hashlib
import json
import re
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

from vec import t1_state_transport as ts  # noqa: E402

REAL_FENG_DIRECTIONS = feng_dirs.feng_directions  # taken at import, before any fixture replaces it on the module
GROUPS4 = ("CM", "Endocardium", "Endothelium", "Epicardium")
NONE_RECORD = {"file": None, "sha256": None, "groups": {}, "feng_groups": "none"}
TARGET = ["--g-mom", "--base", "lib3"]  # the first upload's command, the one design --feng-groups none builds
EDITED = {"V-CM", "Endocardium-1", "EndoMT Endocardium"}  # the fixture's cm and endocardium labels but Early VIC
ALLOWED = re.escape("--feng-groups none builds one file: the first upload's command (t1 --g-mom --base lib3)")


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


def _lib3(env) -> None:
    """The fixture's E10.5 cell names with the dissected library's suffix: two cells in three end in _3, the rest in
    _1, so --base lib3 has a pool of every label (the names end in _<library number>, as the released file's)."""
    env.ref.obs_names = [f"{n}_{'3' if i % 3 else '1'}" for i, n in enumerate(env.ref.obs_names)]


def _target_env(bt, env, monkeypatch) -> list:
    """The fixture made ready for the target command: library suffixes, the real none call, the E9.5 load g_mom's fit
    needs, and a spy that keeps, for each build_submission call with a hook, the population handed to the hook
    (with its counts) and the population the hook returns."""
    _lib3(env)
    _real_none(bt, env, monkeypatch)
    env.loads = _g_mom_loads(bt, env, monkeypatch)  # this load logs file names, not env.reads
    seen = []
    spy = bt.t1_scope.build_submission  # the fixture's spy (it logs the held lists), around the real builder

    def capture(*a, **k):
        hook = k.get("post_cells")
        if hook is None:
            return spy(*a, **k)
        entry = {"kwargs": k}

        def wrapped(pop):
            entry["before"] = (pop.copy(), _dense(pop.layers["counts"]))
            out = hook(pop)
            entry["after"] = (out.copy(), _dense(out.layers["counts"]))
            return out

        seen.append(entry)
        return spy(*a, **{**k, "post_cells": wrapped})

    monkeypatch.setattr(bt.t1_scope, "build_submission", capture)
    return seen


def _fail(what: str):
    return lambda *a, **k: pytest.fail(what)


def _kinds(reads) -> list[str]:
    """The fixture's log, by kind: "E10.5" for a load of the reference, "Feng" for a synthetic answer, "none"."""
    return [r if isinstance(r, str) else r[0] for r in reads]


def _dense(x) -> np.ndarray:
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _released(env, pop) -> np.ndarray:
    """The released counts of each cell of `pop`, by its source cell's name in the reference."""
    rows = env.ref.obs_names.get_indexer(pop.obs["source_obs"].astype(str))
    assert (rows >= 0).all()
    return _dense(env.ref.layers["counts"])[rows]


def _edited_labels(seen_entry) -> set[str]:
    pop, before = seen_entry["before"]
    _, after = seen_entry["after"]
    labels = pop.obs["celltype"].astype(str).to_numpy()
    return set(labels[(before != after).any(axis=1)])


# feng_dirs.feng_directions ---------------------------------------------------------------------------------------


def test_the_groups_are_authors_then_none():
    # "recomputed" (the permitted-only relabelling, tests/test_feng_relabel.py) comes after them: authors stays first
    assert feng_dirs.FENG_GROUPS == ("authors", "none", "recomputed")


def test_none_gives_no_direction_and_reads_nothing(monkeypatch):
    for name in ("load_config", "clean_stats_path", "directions_from_stats"):
        monkeypatch.setattr(feng_dirs, name, _fail(f"groups='none' must not call {name}"))
    monkeypatch.setattr(feng_dirs.fio, "verify_manifest_output", _fail("groups='none' reads no manifest"))
    monkeypatch.setattr(feng_dirs.dl, "load_stats", _fail("groups='none' reads no Feng statistics"))
    d, rec = feng_dirs.feng_directions(["g1", "g2"], groups="none")
    assert d == {} and rec == NONE_RECORD
    d, rec = feng_dirs.feng_directions(["g1"], bracket=("E14.5", "E18"), groups="none")
    assert d == {} and rec == NONE_RECORD  # no step, so no bracket is recorded


def test_an_unknown_groups_value_is_refused_before_anything_is_read(monkeypatch):
    monkeypatch.setattr(feng_dirs, "load_config", _fail("refused before the config is read"))
    for bad in ("Recomputed", "jax", "", "None"):  # "recomputed" is a value since the relabeller's port
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
    for bad in ("Recomputed", "jax"):  # "recomputed" is a value since the relabeller's port (test_feng_relabel.py)
        with pytest.raises(SystemExit):
            ap.parse_args(["t1", "--feng-groups", bad, "--out", "x.h5ad"])
    assert bt.FENG_GROUPS == ("authors", "none", "recomputed") and bt.FENG_GROUPS_DEFAULT == "authors"
    assert bt.FENG_NONE_COMMAND == "t1 --g-mom --base lib3"


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
        rec = json.loads(o.with_suffix(".build.json").read_text())
        assert "feng_groups" not in rec and rec["fraction"] == bt.feng_dirs.bracket_fraction(10.5, 12.5)


def test_the_default_g_mom_record_is_as_before(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """The entry's g_mom record keeps its wording: the none wording is written only under none."""
    _lib3(t1_env)
    _g_mom_loads(bt, t1_env, monkeypatch)
    out = tmp_path / "T1.h5ad"
    bt.main(["t1", *TARGET, "--no-feng-hygiene", "--out", str(out)])
    g = json.loads(out.with_suffix(".build.json").read_text())["g_mom"]
    assert "applied last (after Feng's direction)" in g["rule"]
    assert g["held_flat"] == "the candidate's pre-registered list (R2's union is held in Feng's direction only)"


def test_none_is_the_entry_with_no_direction_and_reads_no_outside_file(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """The target command's file equals t1_scope.build_submission called directly with directions={} and g_mom's hook,
    on the entry's pool (library 3), shares, count and seed, with the pre-registered held list alone; R2's union and
    the Blood guard's call for Feng's group names are never made."""
    _target_env(bt, t1_env, monkeypatch)
    monkeypatch.setattr(bt.hygiene, "load_union", _fail("R2's union is not read under --feng-groups none"))
    monkeypatch.setattr(bt, "_feng_group_names", _fail("Feng's group names are not asked for under none"))
    out = tmp_path / "T1_none.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(out)])
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["feng_groups"] == {"groups": "none", "file": None, "sha256": None,
                                  "direction_step": "removed: no cell is moved along any direction; \"entry\": \"M1\" "
                                                    "names the design this file comes from"}
    assert rec["entry"] == "M1" and rec["fraction"] is None and rec["base"] == "lib3" and rec["check_ok"]
    assert rec["feng_hygiene"] is None and rec["hold_unbracketed"] is True and rec["external"] == NONE_RECORD
    assert t1_env.held_calls == [list(T1_HELD)]
    assert _kinds(t1_env.reads) == ["none"]  # one call, no synthetic Feng answer
    assert t1_env.loads == ["E10.5_RNA.h5ad", "E9.5_RNA.h5ad"]  # E10.5, then E9.5 for g_mom's fit, nothing else
    # (b) the same bytes from a direct call
    ref = t1_env.ref.copy()
    fit = [c for c, s in bt.harness.t1_chamber_shares(bt.harness._shallow(ref), "E10.5_RNA").items() if s > 0]
    in3 = np.array([str(x).rsplit("_", 1)[-1] == "3" for x in ref.obs_names])
    pool = ref[np.flatnonzero(in3)].copy()
    shares = bt.harness.t1_chamber_shares(bt.harness._shallow(pool), "E10.5_RNA")
    groups = bt.feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    seed = bt.harness.BUILD_SEED_OFFSET + 1
    hook, _ = bt._t1_g_mom(ref, fit, seed, bt._g_mom_covered(groups, {}, "none"), feng_direction=False)
    sub, _ = t1_env.build(pool, t1_env.board, shares, bt.harness.count_for(t1_env.board, ref.n_obs), seed,
                          held_flat=list(T1_HELD), ref_source="E10.5_RNA", directions={}, direction_fraction=0.0,
                          direction_groups=groups, herd=False, post_cells=hook)
    direct = tmp_path / "direct.h5ad"
    bt.t1_scope.write_checked(sub, direct, t1_env.board)
    assert _sha(out) == _sha(direct)


def test_none_draws_the_entrys_cells_and_changes_only_the_cells_feng_moved(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """(a) Against the entry on the same inputs and seed, both with the target command: the hook gets the same source
    cells in the same order; under none their counts before the hook are the released ones; the hook edits the same
    labels in both (V-CM, Endocardium-1, EndoMT Endocardium) and leaves Early VIC alone; every cell it does not edit
    keeps its released counts; a cell whose label has no Feng group is the same in both files; and the g_mom record
    differs only in its wording."""
    seen = _target_env(bt, t1_env, monkeypatch)
    entry, none = tmp_path / "T1.h5ad", tmp_path / "T1_none.h5ad"
    bt.main(["t1", *TARGET, "--no-feng-hygiene", "--out", str(entry)])
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(none)])
    assert len(seen) == 2 and seen[0]["kwargs"]["directions"] and seen[1]["kwargs"]["directions"] == {}
    (pop_e, before_e), (pop_n, before_n) = seen[0]["before"], seen[1]["before"]
    assert list(pop_e.obs["source_obs"]) == list(pop_n.obs["source_obs"]) and pop_n.n_obs > 0
    assert np.array_equal(before_n, _released(t1_env, pop_n))  # no step before the hook under none
    assert not np.array_equal(before_e, _released(t1_env, pop_e))  # the entry's direction step ran before its hook
    assert _edited_labels(seen[1]) == _edited_labels(seen[0]) == EDITED
    for s in seen:
        pop, _ = s["before"]
        labels = pop.obs["celltype"].astype(str).to_numpy()
        _, after = s["after"]
        vic = labels == "Early VIC"
        assert vic.any() and np.array_equal(after[vic], _released(t1_env, pop)[vic])
    labels = pop_n.obs["celltype"].astype(str).to_numpy()
    unedited = ~np.isin(labels, sorted(EDITED))
    assert unedited.any() and np.array_equal(seen[1]["after"][1][unedited], _released(t1_env, pop_n)[unedited])
    a, b = ad.read_h5ad(entry), ad.read_h5ad(none)
    assert list(a.obs_names) == list(b.obs_names)
    xa, xb = _dense(a.X), _dense(b.X)
    groups = bt.feng_dirs.label_groups(t1_env.ref.obs["celltype"].astype(str).unique())
    moved = {lab for lab, g in groups.items() if g in t1_env.dirs}
    held = np.array([lab not in moved for lab in labels])
    assert held.any() and (~held).any()
    assert np.array_equal(xa[held], xb[held])  # the labels R1 holds: untouched by the option
    assert not np.array_equal(xa[~held], xb[~held])  # the entry moved these along the (synthetic) direction
    ge, gn = (json.loads(o.with_suffix(".build.json").read_text())["g_mom"] for o in (entry, none))
    wording = {"rule": ge["rule"].replace("(after Feng's direction)", "(no Feng direction in this file)")
               .replace("R1's held labels unedited", "labels the readout config's regexes give no group unedited"),
               "held_flat": "the candidate's pre-registered list"}
    assert "mask" not in ge and "is not read" in gn["mask"]
    assert "R1's held labels" in ge["rule"] and "R1's held labels" not in gn["rule"]
    assert {**ge, **wording, "mask": gn["mask"]} == gn


def test_none_with_g_mom_edits_the_entrys_labels(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    _target_env(bt, t1_env, monkeypatch)
    out = tmp_path / "T1_none_gmom.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(out)])
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["g_mom"]["categories_with_factors"] and rec["feng_groups"]["groups"] == "none"
    assert "applied last (no Feng direction in this file)" in rec["g_mom"]["rule"]
    assert "Feng's direction)" not in rec["g_mom"]["rule"] and "R2" not in rec["g_mom"]["held_flat"]
    groups = {"V-CM": "CM", "Endocardium-1": "Endocardium", "BEC": "Endothelium", "Blood": "Blood",
              "Epicardium": "Epicardium", "Macrophage": "Macrophage", "Early VIC": "other:Early VIC"}
    assert bt._g_mom_covered(groups, {}, "none") == set(groups) - {"Early VIC"}  # every label the map gives a group
    assert bt._g_mom_covered(groups, {"CM": 1}) == {"V-CM"}  # the default: the labels Feng's direction moves
    assert bt._g_mom_covered(groups, {"CM": 1}, "authors") == {"V-CM"}


def test_none_does_not_read_r1_feng_groups(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """Under none g_mom's mask is the readout config's own rule. The target command writes the same bytes whether
    R1_FENG_GROUPS (a list fixed by counting Feng cells) is as committed, empty, or every harmonised group of the
    readout config, and R1's mask function is never called."""
    _target_env(bt, t1_env, monkeypatch)
    out = tmp_path / "T1_none.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(out)])
    cfg = feng_dirs.load_config()
    every = tuple(sorted({g for g, _ in cfg["harmonise"]} | {"CM"}))
    assert set(ts.R1_FENG_GROUPS) < set(every)
    for name, groups in (("empty", ()), ("every", every), ("no_endocardium", ("CM", "Endothelium", "Epicardium"))):
        monkeypatch.setattr(ts, "R1_FENG_GROUPS", groups)
        other = tmp_path / f"T1_none_{name}.h5ad"
        bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(other)])
        assert _sha(out) == _sha(other), name
    monkeypatch.setattr(ts, "r1_unbracketed_labels", _fail("R1's mask is not asked for under --feng-groups none"))
    again = tmp_path / "T1_none_again.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(again)])
    assert _sha(out) == _sha(again)
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert "list of groups with a Feng direction is not read" in rec["g_mom"]["mask"]


def test_the_entry_still_uses_r1s_mask_and_its_record_has_no_mask_field(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """With the option not given the hook asks R1's mask function, as before, and never the readout config's own
    rule; and dropping Endocardium from R1_FENG_GROUPS changes the entry's bytes (the list binds there, which is why
    none does not read it)."""
    _target_env(bt, t1_env, monkeypatch)
    calls = []
    real = ts.r1_unbracketed_labels
    monkeypatch.setattr(ts, "r1_unbracketed_labels", lambda labels, *a, **k: calls.append(1) or real(labels, *a, **k))
    monkeypatch.setattr(bt, "_labels_without_group", _fail("the entry's mask is R1's"))
    entry = tmp_path / "T1.h5ad"
    bt.main(["t1", "--no-feng-hygiene", *TARGET, "--out", str(entry)])
    assert calls and "mask" not in json.loads(entry.with_suffix(".build.json").read_text())["g_mom"]
    monkeypatch.setattr(ts, "R1_FENG_GROUPS", ("CM", "Endothelium", "Epicardium"))
    other = tmp_path / "T1_no_endocardium.h5ad"
    bt.main(["t1", "--no-feng-hygiene", *TARGET, "--out", str(other)])
    assert _sha(entry) != _sha(other)


def test_the_two_masks_agree_on_every_label_the_committed_map_can_place_in_cm_or_endocardium(bt):
    """Committed configs only, no data: on every label configs/t1_chamber_classes.yaml can place in a category g_mom
    edits (its labels, per-file overrides and stand-ins), R1's mask and the readout config's own rule hold back the same
    labels, Early VIC alone; and every other such label maps to CM, Endocardium or Endothelium."""
    import yaml

    cfg = yaml.safe_load((ROOT / "configs" / "t1_chamber_classes.yaml").read_text())
    edit = set(bt.G_MOM_EDIT_CATEGORIES)
    labels = {str(lab) for lab, cat in cfg["labels"].items() if cat in edit}
    for per_file in cfg["file_overrides"].values():
        labels |= {str(lab) for lab, cat in per_file.items() if cat in edit}
    for per_file in cfg["file_fallbacks"].values():
        for cat, stand_ins in per_file.items():
            if cat in edit:
                labels |= set(map(str, stand_ins))
    assert {"V-CM", "Endocardium", "Early VIC"} <= labels
    groups = feng_dirs.label_groups(sorted(labels))
    own_rule = bt._labels_without_group(labels)
    assert own_rule == ts.r1_unbracketed_labels(labels) == {"Early VIC"}
    assert {groups[lab] for lab in labels - own_rule} <= {"CM", "Endocardium", "Endothelium"}


def test_none_with_g_mom_reads_no_r2_union_and_asks_for_no_feng_group_names(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """(c) and (d): with --g-mom, R2's union and Feng's group names are never asked for, and an explicit
    --no-feng-hygiene writes the same bytes as the option alone."""
    _target_env(bt, t1_env, monkeypatch)
    monkeypatch.setattr(bt.hygiene, "load_union", _fail("R2's union is not read under --feng-groups none"))
    monkeypatch.setattr(bt, "_feng_group_names", _fail("Feng's group names are not asked for under none"))
    out, off = tmp_path / "T1_none.h5ad", tmp_path / "T1_none_r2off.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(out)])
    bt.main(["t1", *TARGET, "--no-feng-hygiene", "--feng-groups", "none", "--out", str(off)])
    assert _sha(out) == _sha(off)
    assert t1_env.held_calls == [list(T1_HELD), list(T1_HELD)]


def test_the_target_command_builds_on_names_with_the_dissected_librarys_suffix(bt, t1_env, monkeypatch, tmp_path):  # noqa: F811
    """(e) "--g-mom --base lib3 --feng-groups none" builds; its pool is library 3's cells alone."""
    seen = _target_env(bt, t1_env, monkeypatch)
    out = tmp_path / "T1_s2none.h5ad"
    bt.main(["t1", *TARGET, "--feng-groups", "none", "--out", str(out)])
    assert out.exists() and json.loads(out.with_suffix(".build.json").read_text())["check_ok"]
    pop, _ = seen[0]["before"]
    assert all(str(x).endswith("_3") for x in pop.obs["source_obs"])


def _no_reads(bt, monkeypatch) -> None:
    for obj, name in ((bt, "load_stage"), (bt, "load_boards"), (bt, "_feng_group_names"),
                      (bt.hygiene, "load_union"), (bt.feng_dirs, "feng_directions"),
                      (bt.t1_scope, "held_flat_genes"), (bt.strong_null, "check_approved")):
        monkeypatch.setattr(obj, name, _fail("nothing is read before the refusal"))


@pytest.mark.parametrize("argv", [
    ["--base", "lib3"],  # without --g-mom
    [],
    ["--g-mom"],  # the pooled base
    ["--g-mom", "--base", "pooled"],
    ["--g-mom", "--base", "mid", "--cm-share", "0.60"],
    ["--g-mom", "--base", "mid"],
    ["--g-mom", "--base", "lib3", "--cm-share", "0.60"],
    ["--g-mom", "--base", "lib3", "--blood-share", "0.10"],
    ["--g-mom", "--base", "lib3", "--lib3-pair", "--cm-share", "0.60"],
    ["--g-mom", "--base", "lib3", "--lib3-pair"],
    ["--g-mom", "--base", "lib3", "--herd"],
    ["--g-mom", "--herd"],
    ["--g-mom", "--cm-share", "0.60"],
    ["--g-mom", "--blood-share", "0.10"],
])
def test_none_builds_only_the_first_uploads_command(bt, monkeypatch, tmp_path, argv):
    _no_reads(bt, monkeypatch)
    out = tmp_path / "x.h5ad"
    with pytest.raises(SystemExit, match=ALLOWED):
        bt.main(["t1", *argv, "--feng-groups", "none", "--out", str(out)])
    assert not out.exists() and not out.with_suffix(".build.json").exists()


def test_none_is_refused_with_the_options_that_act_on_fengs_direction(bt, monkeypatch, tmp_path):
    _no_reads(bt, monkeypatch)
    ap, out = bt.build_parser(), str(tmp_path / "x.h5ad")
    for flag, msg in (("--no-hold-unbracketed", "refused with --no-hold-unbracketed"),
                      ("--feng-hygiene", "refused with --feng-hygiene"),
                      ("--strong-null", "refused with --strong-null")):
        with pytest.raises(SystemExit, match=f"--feng-groups none is {msg}"):
            bt.build_t1(ap.parse_args(["t1", *TARGET, flag, "--feng-groups", "none", "--out", out]))
        with pytest.raises(SystemExit, match=f"--feng-groups none is {msg}"):
            bt.main(["t1", *TARGET, flag, "--feng-groups", "none", "--out", out])
