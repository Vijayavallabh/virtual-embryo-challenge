"""Try T72 (experiments/harness/TRIES.md), candidate MX2s: experiments/t1_m1/mx2s_read.py and its launcher mx2s.sh on
synthetic stages, synthetic tables, hand-made task files and temporary git repositories only. No released file, outside
table, statistics file, board file or score is read; the committed T63 task folder is checked for its file names only.
It reuses row T64's and T71's fixtures."""

from __future__ import annotations

import hashlib
import importlib.util
import inspect
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from test_g_mom_verify import _small_inputs
from test_t1_scope import GENES, HELD
from test_t1_state_transport import _gmom_ctx

from vec import harness
from vec.core import counts_of
from vec.t1_scope import classify, load_spec, reweight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "t1_m1" / "mx2s_read.py"
LAUNCHER = ROOT / "experiments" / "t1_m1" / "mx2s.sh"
SHARES = {"cm": 0.6, "endocardium": 0.25, "blood": 0.15}
DIRECTED = ("CM", "Endocardium", "Endothelium")
HELD_M1 = sorted({*HELD, GENES[3], GENES[7]})
WITHIN = {"surrounding": ["Foregut", "Epithelium"]}
C40 = "c" * 40
FIGURE = r"(?<![E\d])\d+\.\d"  # a printed figure; a stage name such as E9.5 is not one
YSHA = "y" * 64


@pytest.fixture
def mx(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("mx2s_read_module", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "RESULTS", tmp_path / "results")
    monkeypatch.setattr(mod, "ATLASES_YAML_SHA256", YSHA)
    derived = tmp_path / "data" / "derived" / "t1_ext_t72"
    monkeypatch.setattr(mod, "DERIVED", derived)
    monkeypatch.setattr(mod, "TABLES", derived / "tables")
    monkeypatch.setattr(mod, "TABLES_MANIFEST", derived / "tables" / "MANIFEST.json")
    monkeypatch.setattr(mod, "CONSENSUS", derived / "consensus")
    return mod


def _dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def _classified():
    ref, start, _ = _small_inputs()
    classify(ref, source="E9.5_RNA")
    classify(start, source="E8.5_RNA")
    return ref, start


def _directions(n_genes, groups=DIRECTED, size=0.4):
    rng = np.random.default_rng(11)
    return {g: rng.normal(0.0, size, n_genes) for g in groups}


def _fb():
    spec = load_spec()
    return {"E8.5": spec.fallbacks("E8.5_RNA"), "E9.5": spec.fallbacks("E9.5_RNA")}


def _entries(dirs, cons=None, mask=None, perm=None, groups=DIRECTED):
    """Each counted arm's changed-group entries; by default copies of X002's own directions."""
    same = {g: dirs[g].copy() for g in groups}
    return {"mx_cons": cons or same, "mx_mask": mask or {g: v.copy() for g, v in same.items()},
            "mx_perm": perm or {g: v.copy() for g, v in same.items()}}


def _lib(mx, e95, e85, entries=None, seed=4):
    factors = mx.br.gmom_fits(e95, "E9.5_RNA", e85, "E8.5_RNA")["file"]
    dirs = _directions(e95.n_vars)
    entries = _entries(dirs) if entries is None else entries(dirs)
    arms, rec = mx.library_arms(e95, SHARES, factors, HELD, HELD_M1, seed, WITHIN, _fb(), dirs, 0.2, entries,
                                n_cells=300)
    base0, _ = reweight(e95, SHARES, 300, seed, within=WITHIN, fallback=_fb()["E9.5"])
    return arms, rec, base0, factors, dirs


def _cut(mx, entries=None, seed=1000002):
    ref, start, _ = _small_inputs()
    ctx = _gmom_ctx(ref, start, seed=seed)
    fits = mx.br.gmom_fits(ref, "E9.5_RNA", start, "E8.5_RNA")
    dirs = _directions(ref.n_vars)
    entries = _entries(dirs) if entries is None else entries(dirs)
    arms, rec = mx.cut_arms(ctx, fits["file"], dirs, 0.2, HELD_M1, entries)
    return arms, rec, ctx, fits, dirs


def _group_mask(mx, labels, dirs, groups):
    g, _, _ = mx.rc.direction_kit(np.unique(labels), dirs)
    return np.array([g.get(lab) in groups for lab in labels], dtype=bool)


def _cut_draw_labels(ctx):
    ref = harness._shallow(ctx.inputs[ctx.task.reference_role])
    src = ctx.task.role_stems[ctx.task.reference_role]
    classify(ref, source=src)
    shares = harness.t1_chamber_shares(ref, src)
    pop, _ = reweight(ref, shares, harness.count_for(ctx.board, ref.n_obs), ctx.seed,
                      fallback=load_spec().fallbacks(src))
    return pop.obs["celltype"].astype(str).to_numpy()


# ---------------------------------------------------------------------------------------------------------------------
# The arms: the base byte for byte when the consensus is X002's direction; only the changed groups' cells change


def test_with_the_consensus_equal_to_x002s_direction_every_arm_is_the_base_on_both_routes(mx):
    e95, e85 = _classified()
    arms, rec, _, _, _ = _lib(mx, e95, e85)
    shas = {k: mx.cj.prediction_sha(v.X) for k, v in arms.items()}
    assert shas["mx_cons"] == shas["mx_mask"] == shas["mx_perm"] == shas[mx.BASE]
    assert all(rec[n]["cells_changed"] == 0 and rec[n]["changed_cells"] > 0 for n in mx.COUNTED)
    arms, _, _, _, _ = _cut(mx)
    shas = {k: mx.cj.prediction_sha(v.X) for k, v in arms.items()}
    assert shas["mx_cons"] == shas["mx_mask"] == shas["mx_perm"] == shas[mx.BASE]


def test_a_mask_with_every_bin_count_zero_is_x002s_direction_and_gives_the_base(mx):
    rng = np.random.default_rng(3)
    n = len(GENES)
    free = ~np.isin(np.asarray(GENES), HELD_M1)
    e95, e85 = _classified()

    def entries(dirs):
        mask = {}
        for g in DIRECTED:
            c = dirs[g] + 1.0  # no exact zero: every bin's count is 0
            m, bins = mx.mask_control(dirs[g], c, rng.gamma(1.0, 10.0, n), free, mx.harness.LP_GROUPS.index(g))
            assert all(b["consensus_zeros"] == 0 == b["mask_zeros"] for b in bins)
            out, k = mx.rescale(m, dirs[g], mx.s_genes(dirs[g], free))
            assert k == 1.0 and out.tobytes() == dirs[g].tobytes()
            mask[g] = out
        return _entries(dirs, mask=mask)

    arms, _, _, _, _ = _lib(mx, e95, e85, entries=entries)
    assert mx.cj.prediction_sha(arms["mx_mask"].X) == mx.cj.prediction_sha(arms[mx.BASE].X)


def test_on_the_library_route_only_the_changed_groups_cells_change(mx):
    e95, e85 = _classified()
    rng = np.random.default_rng(5)

    def entries(dirs):  # CM and Endocardium changed in every arm; Endothelium kept by its absence
        return {a: {g: rng.normal(0, 1.0, e95.n_vars) for g in ("CM", "Endocardium")} for a in mx.COUNTED}

    arms, rec, base0, _, dirs = _lib(mx, e95, e85, entries=entries)
    changed = _group_mask(mx, base0.obs["celltype"].astype(str).to_numpy(), dirs, ("CM", "Endocardium"))
    assert changed.any() and not changed.all()
    a = _dense(counts_of(arms[mx.BASE]))
    for name in mx.COUNTED:
        differs = (a != _dense(counts_of(arms[name]))).any(axis=1)
        assert differs.any() and not differs[~changed].any(), name
        assert rec[name]["outside_equal"] is True and rec[name]["cells_changed"] == int(differs.sum())


def test_on_the_cut_route_only_the_changed_groups_cells_change(mx):
    rng = np.random.default_rng(6)

    def entries(dirs):
        return {a: {"CM": rng.normal(0, 1.0, len(GENES))} for a in mx.COUNTED}

    arms, rec, ctx, _, dirs = _cut(mx, entries=entries)
    changed = _group_mask(mx, _cut_draw_labels(ctx), dirs, ("CM",))
    assert changed.any()
    for name in mx.COUNTED:
        differs = (_dense(arms[mx.BASE].X) != _dense(arms[name].X)).any(axis=1)
        assert differs.any() and not differs[~changed].any(), name
        assert rec[name]["outside_equal"] is True


def test_a_cell_outside_the_changed_groups_that_differs_stops_the_run(mx):
    labels = np.array(["V-CM", "Blood", "Endocardium", "Epicardium"])
    groups = {"V-CM": "CM", "Endocardium": "Endocardium", "Epicardium": "Epicardium"}
    X = np.arange(12, dtype=np.float32).reshape(4, 3)
    Y = X.copy()
    Y[0, 1] += 1.0
    Y[2, 0] += 1.0
    r = mx.byte_record(labels, groups, ["CM", "Endocardium"], X, Y, "mx_cons")
    assert r == {"changed_cells": 2, "cells_changed": 2, "outside_equal": True}
    with pytest.raises(SystemExit, match="outside the changed groups"):
        mx.byte_record(labels, groups, ["CM"], X, Y, "mx_cons")  # Endocardium is not changed here
    Y[3, 2] += 1.0
    with pytest.raises(SystemExit, match="T72"):
        mx.byte_record(labels, groups, ["CM", "Endocardium"], X, Y, "mx_cons")


def test_the_base_is_t63s_doubled_step_on_both_routes(mx):
    rng = np.random.default_rng(7)

    def entries(dirs):
        return {a: {"CM": rng.normal(size=len(GENES))} for a in mx.COUNTED}

    arms, _, ctx, fits, dirs = _cut(mx, entries=entries, seed=1000003)
    x2 = mx.rc.cut_pair(ctx, fits["file"], dirs, 0.4, HELD_M1)["m1rc_gmom"]  # T63's x2 is this call at 2/5
    assert arms[mx.BASE].X.tobytes() == x2.X.tobytes()
    e95, e85 = _classified()
    arms, _, _, factors, dirs = _lib(mx, e95, e85, seed=5)
    t62 = mx.rc.library_pair(e95, SHARES, factors, HELD, HELD_M1, 5, WITHIN, _fb(), dirs, 0.4, n_cells=300)
    assert _dense(arms[mx.BASE].X).tobytes() == _dense(t62["m1rc_gmom"].X).tobytes()
    for fn in (mx.cut_arms, mx.library_arms):  # no truth reaches an arm
        assert not {"truth", "target", "e105"} & set(inspect.signature(fn).parameters)


def test_the_arm_dict_changes_the_changed_groups_entries_alone(mx):
    dirs = _directions(30, groups=(*DIRECTED, "Epicardium"))
    got = mx.arm_directions(dirs, {"CM": np.ones(30), "Endothelium": np.zeros(30)})
    assert set(got) == set(dirs) and np.array_equal(got["CM"], np.ones(30))
    assert got["Endocardium"] is dirs["Endocardium"] and got["Epicardium"] is dirs["Epicardium"]
    with pytest.raises(ValueError):
        mx.arm_directions({"Endocardium": np.ones(3)}, {"CM": np.ones(3)})


# ---------------------------------------------------------------------------------------------------------------------
# The sign rule, the rescale and the controls


NAN = float("nan")


@pytest.mark.parametrize("x002,others,want", [
    (1.0, [], 0.0),                       # X002 alone: n = 1
    (1.0, [NAN, NAN], 0.0),               # n = 1
    (1.0, [0.0, 0.0], 0.0),               # exact zeros are not votes
    (NAN, [1.0, 1.0], NAN),               # NaN in X002 stays NaN
    (0.0, [1.0, 3.0], 2.0),               # X002 zero: two members agree, n = 2
    (1.0, [3.0], 2.0),                    # n = 2, same sign: the mean
    (1.0, [-3.0], 0.0),                   # n = 2, opposite signs
    (1.0, [2.0, 3.0], 2.0),               # 3-0: (1 + 2 + 3) / 3
    (1.0, [2.0, -3.0], 1.0),              # 2-1: the two kept, divided by 3
    (-1.0, [2.0, -5.0], -2.0),            # 2-1 negative
    (1.0, [-2.0, NAN], 0.0),              # n = 2 opposite
    (-2.0, [0.0, -4.0], -3.0),            # n = 2 (zero not counted), same sign
])
def test_the_sign_rules_truth_table(mx, x002, others, want):
    c = mx.sign_consensus(np.array([x002]), [np.array([o]) for o in others])
    if np.isnan(want):
        assert np.isnan(c[0])
    else:
        assert c[0] == pytest.approx(want) and (want != 0.0 or c[0] == 0.0)


def test_the_rescale_holds_x002s_l2_norm_over_s_g(mx):
    rng = np.random.default_rng(2)
    n = 500
    d = rng.normal(size=n)
    d[:7] = np.nan
    free = np.ones(n, dtype=bool)
    free[20:40] = False
    c = mx.sign_consensus(d, [d + rng.normal(0, 0.7, n), d + rng.normal(0, 0.7, n)])
    c[20:40] = 50.0  # held genes outside S_g do not enter the factor
    s = mx.s_genes(d, free)
    out, k = mx.rescale(c, d, s)
    assert np.linalg.norm(out[s]) == pytest.approx(np.linalg.norm(d[s]), rel=1e-12)
    assert np.allclose(out, c * k, equal_nan=True) and np.array_equal(np.isnan(out), np.isnan(d))
    assert mx.rescale(np.where(np.isnan(d), np.nan, 0.0), d, s) [1] is None  # zero norm: no factor


def test_the_mask_keeps_each_bins_zero_count_and_its_seed(mx):
    rng = np.random.default_rng(4)
    n = 600
    d = rng.normal(size=n)
    d[:5] = np.nan
    d[5:8] = 0.0
    level = rng.gamma(1.0, 10.0, n)
    level[8:12] = np.nan
    free = np.ones(n, dtype=bool)
    free[50:60] = False
    c = np.where(rng.random(n) < 0.3, 0.0, d)
    m, rec = mx.mask_control(d, c, level, free, 2)
    s = mx.s_genes(d, free) & (d != 0) & np.isfinite(level)
    bins = mx.bins_of(np.flatnonzero(s), level)
    assert len(rec) == 10 and [r["genes"] for r in rec] == [b.size for b in bins]
    for b, r in zip(bins, rec, strict=True):
        assert r["consensus_zeros"] == int(np.sum(c[b] == 0)) == int(np.sum(m[b] == 0)) == r["mask_zeros"]
    assert np.array_equal(m[~s], d[~s], equal_nan=True)  # only S'_g is masked
    again, _ = mx.mask_control(d, c, level, free, 2)
    assert np.array_equal(again, m, equal_nan=True)
    other, _ = mx.mask_control(d, c, level, free, 1)
    assert not np.array_equal(other, m, equal_nan=True)
    assert mx.MASK_SEED == 20261123 and mx.BINS == 10


def test_the_permutation_keeps_bins_and_nans_with_the_rows_seeds(mx):
    rng = np.random.default_rng(8)
    n = 700
    d = rng.normal(size=n)
    d[:9] = np.nan
    level = rng.gamma(1.0, 10.0, n)
    level[9:14] = np.nan
    free = np.ones(n, dtype=bool)
    free[100:120] = False
    p, rec = mx.perm_member(d, level, free, 1, 0)
    g = free & np.isfinite(d) & np.isfinite(level)
    for b in mx.bins_of(np.flatnonzero(g), level):
        assert np.array_equal(np.sort(p[b]), np.sort(d[b]))
    assert rec["bins_keep_values"] and rec["nan_kept"] and rec["outside_untouched"]
    assert np.array_equal(p[~g], d[~g], equal_nan=True) and not np.array_equal(p[g], d[g])
    seeds = {(m, g2): mx.perm_member(d, level, free, mx.MEMBERS.index(m), mx.harness.LP_GROUPS.index(g2))[1]["seed"]
             for m, g2 in (("X019", "CM"), ("X019", "Endocardium"), ("X019", "Endothelium"), ("X020", "CM"))}
    assert seeds == {("X019", "CM"): 20261134, ("X019", "Endocardium"): 20261135,
                     ("X019", "Endothelium"): 20261136, ("X020", "CM"): 20261144}
    assert mx.harness.LP_GROUPS.index("CM") == 0 and mx.MEMBERS == ("X002", "X019", "X020")


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate on synthetic tables


N = 400
GN = [f"g{i}" for i in range(N)]


def _held():
    h = np.zeros(N, dtype=bool)
    h[:5] = True
    return h


def _x002():
    rng = np.random.default_rng(21)
    d = {g: rng.normal(size=N) for g in DIRECTED}
    for v in d.values():
        v[5:8] = np.nan
    cp = {g: {"E9.5": rng.gamma(2.0, 2.0, N), "E14.5": rng.gamma(2.0, 2.0, N)} for g in DIRECTED}
    return d, cp


def _member_d(d, seed, noise=0.8, nan_share=0.25):
    """A member correlated with X002, NaN on a share of genes (symbols the atlas lacks)."""
    rng = np.random.default_rng(seed)
    out = d + rng.normal(0, noise, N)
    out[rng.random(N) < nan_share] = np.nan
    return out


def _tables(d, cells=(600, 900), cm_x020=None, over=None):
    """{(register, group): table} for X019 in three groups and X020 in CM, correlated with X002."""
    out, over = {}, over or {}
    for i, (reg, g) in enumerate((("X019", "CM"), ("X020", "CM"), ("X019", "Endocardium"),
                                  ("X019", "Endothelium"))):
        dd = _member_d(d[g], 100 + i)
        if reg == "X020" and cm_x020 is not None:
            dd = cm_x020
        out[(reg, g)] = {"D": over.get((reg, g), dd), "m_early": np.random.default_rng(200 + i).gamma(1.0, 5.0, N),
                         "n_early": cells[0], "n_late": cells[1], "file": f"{mx_acc(reg)}/{g}.npz",
                         "sha256": "t" * 64}
    return out


def mx_acc(reg):
    return {"X019": "GSE244783", "X020": "GSE230531"}[reg]


def _manifest(mx, **change):
    atl = {}
    for acc, ends in mx.PINNED.items():
        reg = {"GSE244783": "X019", "GSE230531": "X020"}[acc]
        groups = ("CM", "Endocardium", "Endothelium") if reg == "X019" else ("CM",)
        atl[acc] = {"register": reg, "status": "use", "skm_late_share": 0.05,
                    "samples": {g: {"end": e, "pinned": True, "loaded": True} for e, gs in ends.items() for g in gs},
                    "tables": {g: {"file": f"{acc}/{g}.npz", "sha256": "t" * 64, "n_cells_early": 600,
                                   "n_cells_late": 900} for g in groups}}
    man = {"tasks": ["T1"], "row": "T72", "commit": C40, "clean": True, "yaml_sha256": YSHA,
           "author_labels_read": False, "atlases": atl,
           "constants": {"direction_scale": 0.625, "pseudocount": 0.5, "min_cells": 10, "skm_score": 0.3,
                         "cell_calls": json.loads(json.dumps(mx.CELL_CALLS))}}
    man.update(change)
    return man


X002_N = {"CM": {"E9.5": 700, "E14.5": 1767}}


def _pg(mx, man=None, tables=None, x002_n=None, d=None, table_bad=(), genes_ok=True, yaml_bad=()):
    dd, cp = _x002()
    d = dd if d is None else d
    return mx.pregate_reading(GN, d, cp, X002_N if x002_n is None else x002_n, _held(),
                              _manifest(mx) if man is None else man, _tables(d) if tables is None else tables,
                              list(table_bad), genes_ok, C40, list(yaml_bad))


def test_the_pregate_passes_with_three_changed_groups_and_records_every_figure(mx):
    r = _pg(mx)
    assert r["passes"], r["failed"]
    assert r["changed_groups"] == ["CM", "Endocardium", "Endothelium"] and r["drops"] == []
    cm = r["groups"]["CM"]
    assert cm["members"] == ["X002", "X019", "X020"] and sorted(cm["perm"]) == ["X019", "X020"]
    assert set(cm["arms"]) == set(mx.COUNTED) and all(cm["arms"][a]["k_g"] > 0 for a in mx.COUNTED)
    assert r["groups"]["Endocardium"]["members"] == ["X002", "X019"]
    for g, ga in r["groups"].items():
        s = mx.s_genes(_x002()[0][g], ~_held())
        for v in ga["vectors"].values():
            assert np.linalg.norm(v[s]) == pytest.approx(np.linalg.norm(_x002()[0][g][s]))
        assert len({ga["arms"][a]["vector_sha256"] for a in mx.COUNTED}) == 3
    words = mx.pregate_words(r)
    assert words == ["Pre-gate: PASSES (changed groups: CM, Endocardium, Endothelium; members dropped: none)"]
    assert r["x002_m1"] == {"source": "row N54's record (group_cells)", "cm_cells": {"E9.5": 700, "E14.5": 1767},
                            "bar": 200, "failed_ends": []}
    for g, ga in r["groups"].items():  # the seeds the pre-gate records are the row's formulas, as the read checks
        assert mx.seed_problems(g, ga) == []
        assert ga["mask_seed"] == 20261123 + mx.harness.LP_GROUPS.index(g)


@pytest.mark.parametrize("n,fails", [({"E9.5": 200, "E14.5": 200}, []), ({"E9.5": 199, "E14.5": 900}, ["E9.5"]),
                                     ({"E9.5": 900, "E14.5": 199}, ["E14.5"]), ({"E9.5": 900}, ["E14.5"]),
                                     ({}, ["E9.5", "E14.5"])])
def test_x002s_m1_count_is_recorded_and_a_count_below_200_at_either_end_fails_the_pregate(mx, n, fails):
    r = _pg(mx, x002_n={"CM": n} if n else {})
    assert r["x002_m1"]["cm_cells"] == {e: n.get(e) for e in ("E9.5", "E14.5")}
    assert r["x002_m1"]["failed_ends"] == fails and r["passes"] is (not fails)
    assert "X002" not in {d["member"] for d in r["drops"]}  # X002 is never dropped
    for e in fails:
        assert f"(m1) X002 CM at {e}: below 200 labelled cells in N54's record" in r["failed"]
    if fails:
        words = mx.pregate_words(r)
        assert words[0].startswith("Pre-gate: FAILS ((m1) X002 CM at ") and words[1] == mx.FOLLOWS_PREGATE
        assert not re.search(FIGURE, "\n".join(words))


def test_one_member_dropped_leaves_the_others(mx):
    man = _manifest(mx)
    man["atlases"]["GSE230531"]["skm_late_share"] = 0.11  # X020 fails (m2) in CM
    r = _pg(mx, man=man)
    assert r["passes"] and r["changed_groups"] == ["CM", "Endocardium", "Endothelium"]
    assert [(x["member"], x["group"], x["clauses"]) for x in r["drops"]] == [("X020", "CM", ["(m2)"])]
    assert r["groups"]["CM"]["members"] == ["X002", "X019"]
    assert mx.pregate_words(r)[0].endswith("members dropped: X020, CM, (m2))")


def test_a_group_left_with_one_member_keeps_x002_byte_for_byte(mx):
    d, _ = _x002()
    tables = _tables(d, over={("X019", "Endocardium"): -d["Endocardium"]})  # Spearman -1 fails (m3)
    r = _pg(mx, tables=tables)
    assert r["passes"], r["failed"]
    assert r["changed_groups"] == ["CM", "Endothelium"] and "Endocardium" not in r["groups"]
    assert ("X019", "Endocardium", ["(m3)"]) in [(x["member"], x["group"], x["clauses"]) for x in r["drops"]]
    # no arm table for Endocardium: every arm keeps X002's entry, the same object, so its cells are the base's
    entries = {a: {g: r["groups"][g]["vectors"][a] for g in r["changed_groups"]} for a in mx.COUNTED}
    for a in mx.COUNTED:
        got = mx.arm_directions(d, entries[a])
        assert got["Endocardium"] is d["Endocardium"]


def test_both_new_members_dropped_or_screened_out_fails_the_pregate(mx):
    man = _manifest(mx)
    man["atlases"]["GSE244783"]["status"] = "dropped"
    man["atlases"]["GSE230531"]["status"] = "dropped"
    r = _pg(mx, man=man)
    assert not r["passes"] and "(p2) no changed group" in r["failed"]
    assert {(x["member"], x["group"], tuple(x["clauses"])) for x in r["drops"]} == {
        ("X019", g, ("screens",)) for g in DIRECTED} | {("X020", "CM", ("screens",))}
    words = mx.pregate_words(r)
    assert words[0].startswith("Pre-gate: FAILS (") and "(p2)" in words[0] and words[1] == mx.FOLLOWS_PREGATE
    d, _ = _x002()  # every member dropped by its own clauses: the same failure
    r = _pg(mx, tables=_tables(d, cells=(20, 20)))
    assert not r["passes"] and "(p2) no changed group" in r["failed"] and len(r["drops"]) == 4


def test_a_member_with_no_table_is_dropped_under_m1(mx):
    d, _ = _x002()
    tables = _tables(d)
    del tables[("X020", "CM")]
    r = _pg(mx, tables=tables)
    assert r["passes"] and ("X020", "CM", ["(m1)"]) in [(x["member"], x["group"], x["clauses"]) for x in r["drops"]]


@pytest.mark.parametrize("case,part", [
    ("cm_cells", "(m1)"), ("ec_cells", "(m1)"), ("purity", "(m2)"), ("purity_nan", "(m2)"), ("low", "(m3)"),
    ("high", "(m3)"), ("held", "(m4)")])
def test_each_member_clause_drops_its_member(mx, case, part):
    d, _ = _x002()
    man = _manifest(mx)
    tables = _tables(d)
    reg, g = ("X019", "Endocardium") if case == "ec_cells" else ("X020", "CM")
    t = tables[(reg, g)]
    if case == "cm_cells":
        t["n_late"] = 199
    elif case == "ec_cells":
        t["n_early"] = 49
    elif case == "purity":
        man["atlases"]["GSE230531"]["skm_late_share"] = 0.1000001
    elif case == "purity_nan":
        man["atlases"]["GSE230531"]["skm_late_share"] = None
    elif case == "low":
        t["D"] = np.random.default_rng(1).normal(size=N)
    elif case == "high":
        t["D"] = d[g] * 2.0
    elif case == "held":
        t["D"] = np.where(_held(), 1000.0, t["D"])
    r = _pg(mx, man=man, tables=tables)
    assert (reg, g, [part]) in [(x["member"], x["group"], x["clauses"]) for x in r["drops"]], r["drops"]


def test_the_member_bars_that_pass(mx):
    d, _ = _x002()
    man = _manifest(mx)
    man["atlases"]["GSE230531"]["skm_late_share"] = 0.10
    tables = _tables(d)
    tables[("X020", "CM")].update(n_early=200, n_late=200)
    tables[("X019", "Endocardium")].update(n_early=50, n_late=50)
    r = _pg(mx, man=man, tables=tables)
    assert r["passes"] and r["drops"] == [], r["drops"]


@pytest.mark.parametrize("case,part", [
    ("x002_cm", "(m1) X002"), ("zero_low", "(p3) CM"), ("zero_high", "(p3) Endothelium"), ("spearman", "(p4) CM"),
    ("cons_held", "(p5) CM"), ("p1", "(p1)"), ("record_commit", "record"), ("record_tasks", "record"),
    ("record_yaml", "record"), ("record_authors", "record"), ("record_never", "record"),
    ("record_unpinned", "record"), ("record_constants", "record"), ("table_bad", "record"), ("genes", "record"),
    ("yaml_disk", "record")])
def test_each_line_clause_fails_the_pregate(mx, case, part):
    d, _ = _x002()
    man, tables, kw = _manifest(mx), _tables(d), {}
    if case == "x002_cm":
        kw["x002_n"] = {"CM": {"E9.5": 199, "E14.5": 1767}}
    elif case == "zero_low":  # three members finite everywhere: a majority sign always exists, no zero
        for i, k in enumerate((("X019", "CM"), ("X020", "CM"))):
            tables[k]["D"] = _member_d(d["CM"], 300 + i, nan_share=0.0)
    elif case == "zero_high":  # X019 finite on a quarter of the genes: X002 alone gives 0 elsewhere
        tables[("X019", "Endothelium")]["D"] = _member_d(d["Endothelium"], 310, nan_share=0.75)
    elif case == "spearman":  # members close to X002: their mean is closer still
        for i, k in enumerate((("X019", "CM"), ("X020", "CM"))):
            tables[k]["D"] = _member_d(d["CM"], 320 + i, noise=0.3, nan_share=0.1)
    elif case == "cons_held":  # X002 large on the held genes; the members small there, of the same sign
        d = {g: v.copy() for g, v in d.items()}
        d["CM"][:5] = 200.0
        tables = _tables(d)
        for k in (("X019", "CM"), ("X020", "CM")):
            tables[k]["D"][:5] = 1.0
        kw["d"] = d
    elif case == "p1":
        man["atlases"]["GSE230531"]["samples"]["GSM7226274"]["loaded"] = False
    elif case == "record_commit":
        man["commit"] = "d" * 40
    elif case == "record_tasks":
        man["tasks"] = ["T1", "T2"]
    elif case == "record_yaml":
        man["yaml_sha256"] = "0" * 64
    elif case == "record_authors":
        man["author_labels_read"] = True
    elif case == "record_never":
        man["atlases"]["GSE244783"]["samples"]["GSM7829701"] = {"end": "late", "pinned": True, "loaded": True}
    elif case == "record_unpinned":
        man["atlases"]["GSE230531"]["samples"]["GSM9999999"] = {"end": "late", "pinned": True, "loaded": True}
    elif case == "record_constants":
        man["constants"]["pseudocount"] = 1.0
    elif case == "table_bad":
        kw["table_bad"] = ["GSE230531 CM: the table is missing or not the one the manifest names"]
    elif case == "genes":
        kw["genes_ok"] = False
    elif case == "yaml_disk":
        kw["yaml_bad"] = ["atlases.yaml's sha256 is x, not the row's"]
    r = _pg(mx, man=man, tables=tables, **kw)
    assert not r["passes"] and any(f.startswith(part) for f in r["failed"]), (case, r["failed"])
    words = mx.pregate_words(r)
    assert words[0].startswith("Pre-gate: FAILS (") and words[1] == mx.FOLLOWS_PREGATE


def test_a_rescale_factor_that_does_not_exist_fails_p6(mx):
    d, _ = _x002()
    tables = _tables(d)
    tables[("X019", "Endothelium")]["D"] = -d["Endothelium"]  # (m3) would drop it: force the group's arms instead
    ga = mx.group_arms("Endothelium", d["Endothelium"], {"X019": tables[("X019", "Endothelium")]},
                       _x002()[1]["Endothelium"]["E9.5"], ~_held(), _held())
    assert ga["arms"]["mx_cons"]["k_g"] is None and ga["arms"]["mx_cons"]["zero_share"] == 1.0


def test_the_pregate_words_carry_no_figure(mx):
    for r in (_pg(mx), _pg(mx, man=_manifest(mx, commit="d" * 40))):
        text = "\n".join(mx.pregate_words(r))
        assert not re.search(FIGURE, text)


# ---------------------------------------------------------------------------------------------------------------------
# The pre-gate step end to end: exclusive outputs, exit 3, a refusal on its own outputs


def _write_tables(mx, d, man):
    for (reg, g), t in _tables(d).items():
        acc = mx_acc(reg)
        path = mx.TABLES / acc / f"{g}.npz"
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, genes=np.array(GN), m_early=t["m_early"], m_late=t["m_early"], D=t["D"],
                 n_cells_early=np.array(t["n_early"]), n_cells_late=np.array(t["n_late"]))
        man["atlases"][acc]["tables"][g]["sha256"] = mx.sha256_file(path)
    mx.TABLES_MANIFEST.write_text(json.dumps(man))


def _tries(tmp_path, row="T72", copies=1):
    t = tmp_path / "TRIES.md"
    t.write_text("".join(f"| {row} (one row) | T1 | x |\n" for _ in range(copies)))
    return t


def _pregate_ready(mx, tmp_path, monkeypatch, man=None):
    d, cp = _x002()
    monkeypatch.setattr(mx, "TRIES", _tries(tmp_path))
    monkeypatch.setattr(mx.gv, "tree_state", lambda: {"commit": C40, "clean": True})
    monkeypatch.setattr(mx.gv, "stamp", lambda started, fs: {"commit": C40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(mx, "get_board", lambda b: type("B", (), {"genes": GN})())
    monkeypatch.setattr(mx.rc, "recomputed_inputs", lambda g: ({**d, "Epicardium": np.ones(N)}, {}))
    monkeypatch.setattr(mx.rc, "_held_m1", lambda p: ([GN[i] for i in range(5)], {}))
    monkeypatch.setattr(mx, "x002_group_means", lambda genes, g: cp[g])
    monkeypatch.setattr(mx, "x002_counts", lambda: X002_N)
    monkeypatch.setattr(mx, "yaml_problems", lambda *a: [])
    monkeypatch.setattr(mx.feng_dirs, "recomputed_stats_path", lambda *a: tmp_path / "stats.npz")
    _write_tables(mx, d, _manifest(mx) if man is None else man)

    class A:
        out, commit = "t72", C40
    return A


def test_the_pregate_step_writes_its_outputs_exclusively_and_prints_only_words(mx, tmp_path, monkeypatch, capsys):
    a = _pregate_ready(mx, tmp_path, monkeypatch)
    assert mx.cmd_pregate(a) == 0
    printed = capsys.readouterr().out
    assert printed.startswith("Pre-gate: PASSES (changed groups: CM, Endocardium, Endothelium")
    assert not re.search(FIGURE, printed)
    tasks = mx.RESULTS / "t72" / "tasks"
    text = (tasks / "pregate.json").read_text()
    assert (tasks / "pregate.json.sha256").read_text().split()[0] == hashlib.sha256(text.encode()).hexdigest()
    rec = json.loads(text)
    assert rec["reading"]["passes"] and sorted(rec["arm_tables"]) == ["CM", "Endocardium", "Endothelium"]
    assert len(list(mx.CONSENSUS.glob("*.npz"))) == 9 and mx.arm_table_problems(rec) == []
    assert rec["constants"] == json.loads(json.dumps(mx.pregate_constants()))
    for g, arms in rec["arm_tables"].items():
        assert set(arms) == set(mx.COUNTED) and rec["reading"]["groups"][g]["arms"]["mx_cons"]["k_g"] > 0
    with pytest.raises(SystemExit, match="runs once"):
        mx.cmd_pregate(a)
    (tasks / "pregate.json").rename(tasks / "moved.json")  # the consensus folder alone also refuses
    with pytest.raises(SystemExit, match="runs once"):
        mx.cmd_pregate(a)


def test_the_pregate_step_exits_3_on_a_failure(mx, tmp_path, monkeypatch, capsys):
    man = _manifest(mx)
    for acc in man["atlases"]:
        man["atlases"][acc]["status"] = "dropped"
    a = _pregate_ready(mx, tmp_path, monkeypatch, man=man)
    assert mx.cmd_pregate(a) == mx.PREGATE_FAIL == 3
    printed = capsys.readouterr().out.splitlines()
    assert printed == ["Pre-gate: FAILS ((p2) no changed group)", mx.FOLLOWS_PREGATE]
    assert json.loads((mx.RESULTS / "t72" / "tasks" / "pregate.json").read_text())["reading"]["passes"] is False


def test_the_pregate_refuses_before_reading_anything(mx, tmp_path, monkeypatch):
    monkeypatch.setattr(mx, "TRIES", _tries(tmp_path, "T71"))
    monkeypatch.setattr(mx.gv, "tree_state", lambda: pytest.fail("nothing before the row's refusal"))

    class A:
        out, commit = "t72", C40
    with pytest.raises(SystemExit, match="exactly one row"):
        mx.cmd_pregate(A)
    monkeypatch.setattr(mx, "TRIES", _tries(tmp_path, "T72", copies=2))
    with pytest.raises(SystemExit, match="exactly one row"):
        mx.cmd_pregate(A)
    assert not mx.CONSENSUS.exists()


# ---------------------------------------------------------------------------------------------------------------------
# The run steps' refusals and the start record


def test_a_run_refuses_a_second_start_and_reads_the_arm_tables_by_gene_name(mx, tmp_path, monkeypatch):
    a = _pregate_ready(mx, tmp_path, monkeypatch)
    assert mx.cmd_pregate(a) == 0
    monkeypatch.setattr(mx.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    monkeypatch.setattr(mx, "all_source_bad", lambda: [])
    tasks = mx.RESULTS / "t72" / "tasks"
    started, _, pg, sha = mx.refuse(tasks, "lib1.json", C40)
    assert started["clean"] is True and sha == mx.sha256_file(tasks / "pregate.json")
    mx.write_started(tasks, "lib1.json", C40)
    assert json.loads((tasks / "lib1.started").read_text())["commit"] == C40
    with pytest.raises(SystemExit, match="start record exists"):
        mx.refuse(tasks, "lib1.json", C40)
    with pytest.raises(FileExistsError):
        mx.write_started(tasks, "lib1.json", C40)
    with pytest.raises(SystemExit, match="HEAD"):
        mx.refuse(tasks, "lib2.json", "d" * 40)
    rev = GN[::-1]
    got = mx.arm_entries(pg, rev)
    with np.load(mx.CONSENSUS / pg["arm_tables"]["CM"]["mx_cons"]["file"]) as z:
        assert np.array_equal(got["mx_cons"]["CM"], z["direction"][::-1], equal_nan=True)
    assert sorted(got) == sorted(mx.COUNTED) and sorted(got["mx_perm"]) == ["CM", "Endocardium", "Endothelium"]
    with pytest.raises(SystemExit, match="genes"):
        mx.arm_entries(pg, GN[:-1])
    f = mx.CONSENSUS / pg["arm_tables"]["CM"]["mx_mask"]["file"]
    f.chmod(0o644)
    with open(f, "ab") as h:
        h.write(b"x")
    with pytest.raises(SystemExit, match="arm table"):
        mx.refuse(tasks, "lib3.json", C40)


def test_a_run_is_refused_without_a_passed_pregate_or_on_an_existing_task_file(mx, tmp_path, monkeypatch):
    monkeypatch.setattr(mx, "TRIES", _tries(tmp_path))
    monkeypatch.setattr(mx.gv, "tree_state", lambda: {"commit": C40, "clean": True})
    monkeypatch.setattr(mx.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    monkeypatch.setattr(mx, "yaml_problems", lambda *a: [])
    tasks = tmp_path / "tasks"
    tasks.mkdir()
    with pytest.raises(SystemExit, match="pregate.json is missing"):
        mx.refuse(tasks, "x.json", C40)
    text = json.dumps({"reading": {"passes": False}, "stamp": {"commit": C40, "tree_clean": True}})
    (tasks / "pregate.json").write_text(text)
    (tasks / "pregate.json.sha256").write_text(hashlib.sha256(text.encode()).hexdigest() + "  pregate.json\n")
    with pytest.raises(SystemExit, match="has not passed"):
        mx.refuse(tasks, "x.json", C40)
    (tasks / "x.json").write_text("{}")
    with pytest.raises(SystemExit, match="exists"):
        mx.refuse(tasks, "x.json", C40)
    monkeypatch.setattr(mx.gv, "tree_state", lambda: {"commit": C40, "clean": False})
    with pytest.raises(SystemExit, match="not clean"):
        mx.refuse(tasks, "y.json", C40)


def test_yaml_problems_names_an_unfilled_constant_and_a_file_off_it(mx, tmp_path, monkeypatch):
    y = tmp_path / "atlases.yaml"
    y.write_text("status: frozen\n")
    monkeypatch.setattr(mx, "ATLASES_YAML_SHA256", None)
    assert mx.yaml_problems(y) and "placement" in mx.yaml_problems(y)[0]
    monkeypatch.setattr(mx, "ATLASES_YAML_SHA256", mx.sha256_file(y))
    assert mx.yaml_problems(y) == []
    y.write_text("status: frozen\nx: 1\n")
    assert mx.yaml_problems(y)


@pytest.mark.parametrize("const", ["unfilled", "wrong"])
def test_the_pregate_refuses_on_the_yaml_constant_before_it_creates_anything(mx, tmp_path, monkeypatch, capsys, const):
    real = mx.yaml_problems
    y = tmp_path / "atlases.yaml"
    y.write_text("status: frozen\n")
    a = _pregate_ready(mx, tmp_path, monkeypatch, man=_manifest(mx, yaml_sha256=mx.sha256_file(y)))
    monkeypatch.setattr(mx, "ATLASES_YAML", y)
    monkeypatch.setattr(mx, "yaml_problems", real)
    monkeypatch.setattr(mx, "ATLASES_YAML_SHA256", None if const == "unfilled" else "0" * 64)
    with pytest.raises(SystemExit) as e:
        mx.cmd_pregate(a)
    assert e.value.code not in (0, 3, None) and "atlases.yaml" in str(e.value.code)
    tasks = mx.RESULTS / "t72" / "tasks"
    assert not (tasks / "pregate.json").exists() and not (tasks / "pregate.json.sha256").exists()
    assert not mx.CONSENSUS.exists() and not tasks.exists()
    assert "Pre-gate" not in capsys.readouterr().out
    monkeypatch.setattr(mx, "ATLASES_YAML_SHA256", mx.sha256_file(y))  # filled and matching: the step runs
    assert mx.cmd_pregate(a) == 0 and (tasks / "pregate.json").exists()


def test_an_unfilled_yaml_constant_step_by_step(mx, tmp_path, monkeypatch, capsys):
    """pregate: refuses (non-zero, not 3) before it creates anything (the test above); run-cut and run-lib: refuse
    before their start record; check-pregate: refuses; read: follows the row (a void reason), reached only after a
    passed pre-gate and a start record."""
    real = mx.yaml_problems
    a = _pregate_ready(mx, tmp_path, monkeypatch)
    assert mx.cmd_pregate(a) == 0  # a passed pre-gate, made with the constant filled
    capsys.readouterr()
    tasks = mx.RESULTS / "t72" / "tasks"
    y = tmp_path / "atlases.yaml"
    y.write_text("status: frozen\n")
    monkeypatch.setattr(mx, "ATLASES_YAML", y)
    monkeypatch.setattr(mx, "yaml_problems", real)
    monkeypatch.setattr(mx, "ATLASES_YAML_SHA256", None)
    monkeypatch.setattr(mx.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    monkeypatch.setattr(mx, "all_source_bad", lambda: [])
    monkeypatch.setattr(mx, "load_boards", lambda: pytest.fail("nothing loads after a refusal"))
    monkeypatch.setattr(mx, "get_board", lambda b: pytest.fail("nothing loads after a refusal"))

    class Run:
        out, commit, seeds, data = "t72", C40, "0", str(tmp_path)
        cut, library = "chamber", "1"
    for cmd, task in ((mx.cmd_run_cut, mx.CUT_FILES["chamber"]), (mx.cmd_run_lib, "lib1.json")):
        with pytest.raises(SystemExit, match="placement"):
            cmd(Run)
        assert not (tasks / mx.started_name(task)).exists() and not (tasks / task).exists()
    with pytest.raises(SystemExit, match="placement"):
        mx.cmd_check_pregate(a)
    # the read: the row's void list holds "atlases.yaml's sha256 is not the row's"
    (tasks / "lib1.started").write_text("{}")
    monkeypatch.setattr(mx, "committed", lambda f, *x: {s: 51.0 for s in range(10)})
    res = mx.read(mx.RESULTS / "t72", C40)
    assert res["reading"]["outcome"] == "VOID"
    assert any("atlases.yaml's sha256 is not the row's" in x for x in res["reading"]["reasons"])
    assert any("placement" in x for x in res["reading"]["reasons"])
    assert capsys.readouterr().out.startswith("Outcome: VOID (")


# ---------------------------------------------------------------------------------------------------------------------
# The reading on hand-made task files


def _arm(score):
    return {"score": score, "groups": {g: score / 100 for g in ("de_recovery", "direction", "distribution",
                                                               "covariance")}}


def _relabel(mx, **change):
    return {"sha256": mx.rc.STATS_SHA256, "stats_sha256_hashed": mx.rc.STATS_SHA256, "yaml_sha256": mx.rc.YAML_SHA256,
            "route": "recomputed", "record_route": "recomputed", "authors_labels_read": False, **change}


ARM_TABLES = {"CM": {a: {"file": f"CM_{a}.npz", "sha256": a[3] * 64, "vector_sha256": a[3] * 64}
                     for a in ("mx_cons", "mx_mask", "mx_perm")}}
PSHA = "p" * 64


def _pregate(mx, passes=True, commit=C40, clean=True, **change):
    pg = {"reading": {"passes": passes, "changed_groups": ["CM"], "groups": {"CM": {
        "group_index": 0, "members": ["X002", "X019", "X020"], "mask_seed": 20261123,
        "mask_bins": [{"consensus_zeros": 3, "mask_zeros": 3}] * 10,
        "perm": {"X019": {"bins_keep_values": True, "nan_kept": True, "seed": 20261134},
                 "X020": {"bins_keep_values": True, "nan_kept": True, "seed": 20261144}},
        "arms": {a: {"vector_sha256": a[3] * 64} for a in mx.COUNTED}}}},
        "stamp": {"commit": commit, "tree_clean": clean}, "constants": json.loads(json.dumps(mx.pregate_constants())),
        "data_constants": {"direction_scale": 0.625, "pseudocount": 0.5, "min_cells": 10, "skm_score": 0.3,
                           "cell_calls": json.loads(json.dumps(mx.CELL_CALLS))},
        "yaml_sha256": YSHA, "author_labels_read": False, "arm_tables": ARM_TABLES,
        "samples": {acc: sorted(g for e in ends.values() for g in e) for acc, ends in mx.PINNED.items()}}
    pg.update(change)
    return pg


def _task(mx, gain=1.0, mask=0.0, perm=0.0, lib=False, commit=C40, clean=True, seeds=range(10), gains=None, **top):
    rows = []
    for s in seeds:
        b = 51.0 + s / 100
        g = gain if gains is None else gains[s]
        arms = {mx.BASE: _arm(b), "mx_cons": _arm(b + g), "mx_mask": _arm(b + mask), "mx_perm": _arm(b + perm)}
        sha = {a: f"{i}" * 64 for i, a in enumerate(arms)}
        rec = {n: {"outside_equal": True, "changed_cells": 30, "cells_changed": 30} for n in mx.COUNTED}
        if lib:
            arms["copy_last"], sha["copy_last"] = _arm(40.0), "c" * 64
        rows.append({"seed": s, "arms": arms, "sha256": sha, "records": rec})
    pg = _pregate(mx)
    out = {"rows": rows, "stamp": {"commit": commit, "tree_clean": clean}, "feng_fraction": 0.2, "multiplier": 2.0,
           "applied_fraction": 0.4, "constants": json.loads(json.dumps(mx.constants(HELD_M1, pg))),
           "relabel": _relabel(mx), "pregate_sha256": PSHA, "arm_tables": ARM_TABLES, "yaml_sha256": YSHA}
    out.update(top)
    return out


def _gate():
    return {k: {s: 51.0 + s / 100 for s in range(10)} for k in ("chamber", "whole", "1", "2", "3")}


def _loaded(mx, **kw):
    out = {"chamber": _task(mx, **kw), "whole": _task(mx, **kw)}
    out.update({lib: _task(mx, lib=True, **kw) for lib in ("1", "2", "3")})
    return out


OK = {"reasons": []}


def _read(mx, loaded, gate=None, pregate=None, checks=OK, **kw):
    return mx.t72_reading(loaded, _gate() if gate is None else gate, checks,
                          _pregate(mx) if pregate is None else pregate, PSHA, **kw)


def test_ahead_when_all_three_clauses_hold(mx):
    r = _read(mx, _loaded(mx))
    assert r["outcome"] == "AHEAD" and r["failed"] == []
    lines = mx.outcome_lines(r)
    assert lines[0] == ("Outcome: MX2s is AHEAD (clauses (i), (ii) and (iii) hold at the test's fraction 2/5)")
    assert all(b in lines for b in mx.BESIDE) and lines[-1] == mx.FOLLOWS_AHEAD
    assert ("Deciding nothing: mx_mask would not meet clause (i) against the base; mx_perm would not meet clause (i) "
            "against the base") in lines
    assert r["j6_range"] == [0.0, pytest.approx(1.0)]


def test_not_ahead_names_each_failed_clause_and_closes_the_line(mx):
    r = _read(mx, _loaded(mx, gain=0.1))
    assert r["outcome"] == "NOT AHEAD" and r["follows"] == mx.FOLLOWS_END
    line = mx.outcome_lines(r)[0]
    assert line.startswith("Outcome: MX2s is NOT AHEAD ((i) the chamber") and "(ii)" in line and "(iii)" in line
    assert "'multi-atlas sign consensus'" in r["follows"]


def _p(mean, ahead=10):
    return {"mean": mean, "ahead": ahead}


def _good(v=1.0):
    return {k: _p(v) for k in ("chamber", "whole", "1", "2", "3")}


@pytest.mark.parametrize("part,ab,am,ap", [
    ("(i) the chamber", {"chamber": _p(0.2999999999)}, {}, {}),
    ("(i) the chamber", {"chamber": _p(1.0, 7)}, {}, {}),
    ("(i) the whole", {"whole": _p(0.0)}, {}, {}),
    ("(i) the mean over the three libraries", {"1": _p(0.2), "2": _p(0.2), "3": _p(0.4)}, {}, {}),
    ("(i) library 3", {"3": _p(1.0, 6)}, {}, {}),
    ("(i) library 3", {"3": _p(0.0)}, {}, {}),
    ("(i) libraries", {"2": _p(-0.31), "1": _p(1.0), "3": _p(1.0)}, {}, {}),
    ("(ii) the chamber", {}, {"chamber": _p(0.19)}, {}),
    ("(ii) the chamber", {}, {"chamber": _p(1.0, 7)}, {}),
    ("(ii) the library", {}, {"1": _p(0.1), "2": _p(0.1), "3": _p(0.3)}, {}),
    ("(iii) the chamber", {}, {}, {"chamber": _p(0.19)}),
    ("(iii) the chamber", {}, {}, {"chamber": _p(1.0, 7)})])
def test_not_ahead_on_each_clause(mx, part, ab, am, ap):
    ahead, failed = mx.statement({**_good(), **ab}, {**_good(), **am}, {**_good(), **ap})
    assert not ahead and any(f.startswith(part) for f in failed), failed
    assert len(failed) == 1


def test_the_boundaries_that_pass_library_3_needs_7_and_the_chamber_8(mx):
    ab = {"chamber": _p(0.3, 8), "whole": _p(1e-12), "1": _p(0.6), "2": _p(-0.3), "3": _p(0.6, 7)}
    am = {"chamber": _p(0.2, 8), "whole": _p(-5.0), "1": _p(0.2), "2": _p(0.2), "3": _p(0.2, 0)}
    ap = {"chamber": _p(0.2, 8), "whole": _p(-5.0), "1": _p(-9.0), "2": _p(-9.0), "3": _p(-9.0, 0)}
    ahead, failed = mx.statement(ab, am, ap)
    assert ahead, failed
    assert (mx.MIN_AHEAD, mx.LIB3_MIN_AHEAD) == (8, 7)


def test_an_exact_zero_is_not_ahead(mx):
    gains = [0.0, 0.0, 0.0] + [1.0] * 7  # seven of ten ahead on the chamber cut
    r = _read(mx, _loaded(mx, gains=gains))
    assert r["outcome"] == "NOT AHEAD" and r["tasks"]["mx_cons_minus_m1rc_gmom_x2"]["chamber"]["ahead"] == 7
    assert r["tasks"]["mx_cons_minus_m1rc_gmom_x2"]["chamber"]["behind"] == 0


def test_a_control_ahead_of_the_base_is_printed_and_decides_nothing(mx):
    r = _read(mx, _loaded(mx, gain=2.0, mask=1.0, perm=0.5))
    assert r["outcome"] == "AHEAD" and r["controls_meet_clause_i"] == {"mx_mask": True, "mx_perm": True}
    assert "mx_mask would meet clause (i)" in mx.outcome_lines(r)[-2]


def test_an_idle_task_is_named_for_every_counted_arm(mx):
    loaded = _loaded(mx)
    for x in loaded["whole"]["rows"]:
        x["sha256"]["mx_cons"] = x["sha256"][mx.BASE]
    for x in loaded["2"]["rows"]:
        x["sha256"]["mx_mask"] = x["sha256"][mx.BASE]
    for x in loaded["chamber"]["rows"]:
        x["sha256"]["mx_perm"] = x["sha256"][mx.BASE]
    loaded["3"]["rows"][4]["sha256"]["mx_perm"] = loaded["3"]["rows"][4]["sha256"][mx.BASE]  # one seed: not idle
    r = _read(mx, loaded)
    lines = mx.outcome_lines(r)
    idle = [ln for ln in lines if "moved no cell" in ln]
    assert idle == ["mx_cons: the arm moved no cell on whole: nothing is priced there",
                    "mx_mask: the arm moved no cell on 2: nothing is priced there",
                    "mx_perm: the arm moved no cell on chamber: nothing is priced there"]
    txt = mx.t72_text({"reading": r, "label": "x", "tries": 15, "stamp": {"commit": C40, "tree_clean": True,
                                                                          "utc": "x"}})
    assert all(txt.count(ln) == 1 for ln in idle)


VOID_CASES = ["missing_task", "started_no_task", "task_no_start", "missing_seed", "missing_arm", "missing_sha",
              "missing_copy_last", "score_nan", "two_commits", "unclean", "read_commit", "gate", "gate_seed",
              "multiplier", "fraction", "applied", "constant", "pregate_constant", "data_constant", "stats", "settings",
              "route", "authors_relabel", "checks", "source", "yaml_task", "yaml_pregate", "yaml_read", "never",
              "unpinned", "author_labels", "mask_bin", "perm_values", "perm_nan", "mask_equals_cons",
              "perm_equals_cons", "byte", "pregate_missing", "pregate_failed", "pregate_unclean", "pregate_commit",
              "pregate_tables", "pregate_sha", "pregate_bad", "mask_seed", "perm_seed_x019", "perm_seed_x020",
              "perm_members", "group_index", "group_unknown"]


@pytest.mark.parametrize("case", VOID_CASES)
def test_every_void_reason_fires(mx, case):
    loaded, gate, pregate, checks = _loaded(mx), _gate(), _pregate(mx), OK
    kw = {"started": {k: True for k in loaded}, "source_bad": [], "yaml_bad": [], "read_commit": C40,
          "pregate_bad": []}
    psha = PSHA
    assert mx.t72_reading(loaded, gate, checks, pregate, psha, **kw)["outcome"] != "VOID"
    if case == "missing_task":
        loaded["2"] = None
        kw["started"]["2"] = False
    elif case == "started_no_task":
        loaded["2"] = None
    elif case == "task_no_start":
        kw["started"]["chamber"] = False
    elif case == "missing_seed":
        loaded["whole"] = _task(mx, seeds=range(9))
    elif case == "missing_arm":
        del loaded["chamber"]["rows"][2]["arms"]["mx_perm"]
    elif case == "missing_sha":
        del loaded["whole"]["rows"][5]["sha256"]["mx_mask"]
    elif case == "missing_copy_last":
        del loaded["3"]["rows"][0]["arms"]["copy_last"]
    elif case == "score_nan":
        loaded["1"]["rows"][3]["arms"]["mx_cons"]["score"] = float("nan")
    elif case == "two_commits":
        loaded["1"] = _task(mx, lib=True, commit="d" * 40)
    elif case == "unclean":
        loaded["chamber"] = _task(mx, clean=False)
    elif case == "read_commit":
        kw["read_commit"] = "e" * 40
    elif case == "gate":
        loaded["2"]["rows"][4]["arms"][mx.BASE]["score"] += 2e-9
    elif case == "gate_seed":
        del gate["3"][7]
    elif case == "multiplier":
        loaded["whole"]["multiplier"] = 1.5
    elif case == "fraction":
        loaded["chamber"]["feng_fraction"] = 0.4
    elif case == "applied":
        loaded["chamber"]["applied_fraction"] = 0.6
    elif case == "constant":
        loaded["2"]["constants"]["mask_seed"] = 1
    elif case == "pregate_constant":
        pregate["constants"]["bins"] = 5
    elif case == "data_constant":
        pregate["data_constants"]["direction_scale"] = 0.5
    elif case == "stats":
        loaded["chamber"]["relabel"] = _relabel(mx, stats_sha256_hashed="0" * 64)
    elif case == "settings":
        loaded["chamber"]["relabel"] = _relabel(mx, yaml_sha256="0" * 64)
    elif case == "route":
        loaded["1"]["relabel"] = _relabel(mx, route="none")
    elif case == "authors_relabel":
        loaded["1"]["relabel"] = _relabel(mx, authors_labels_read=True)
    elif case == "checks":
        checks = {"reasons": ["the settings file's sha256 is x"]}
    elif case == "source":
        kw["source_bad"] = ["the source names a banned token"]
    elif case == "yaml_task":
        loaded["3"]["yaml_sha256"] = "0" * 64
    elif case == "yaml_pregate":
        pregate["yaml_sha256"] = "0" * 64
    elif case == "yaml_read":
        kw["yaml_bad"] = ["atlases.yaml's sha256 is x, not the row's"]
    elif case == "never":
        pregate["samples"]["GSE244783"].append("GSM7829701")
    elif case == "unpinned":
        pregate["samples"]["GSE230531"].append("GSM9999999")
    elif case == "author_labels":
        pregate["author_labels_read"] = True
    elif case == "mask_bin":
        pregate["reading"]["groups"]["CM"]["mask_bins"] = [{"consensus_zeros": 3, "mask_zeros": 2}]
    elif case == "perm_values":
        pregate["reading"]["groups"]["CM"]["perm"]["X020"]["bins_keep_values"] = False
    elif case == "perm_nan":
        pregate["reading"]["groups"]["CM"]["perm"]["X019"]["nan_kept"] = False
    elif case == "mask_equals_cons":
        pregate["reading"]["groups"]["CM"]["arms"]["mx_mask"]["vector_sha256"] = "c" * 64
    elif case == "perm_equals_cons":
        pregate["reading"]["groups"]["CM"]["arms"]["mx_perm"]["vector_sha256"] = "c" * 64
    elif case == "byte":
        loaded["2"]["rows"][6]["records"]["mx_perm"]["outside_equal"] = False
    elif case == "pregate_missing":
        pregate = None
    elif case == "pregate_failed":
        pregate = _pregate(mx, passes=False)
    elif case == "pregate_unclean":
        pregate = _pregate(mx, clean=False)
    elif case == "pregate_commit":
        pregate = _pregate(mx, commit="e" * 40)
    elif case == "pregate_tables":
        pregate = _pregate(mx, arm_tables={"CM": {**ARM_TABLES["CM"], "mx_cons": {"file": "x", "sha256": "0" * 64}}})
    elif case == "pregate_sha":
        psha = "0" * 64
    elif case == "pregate_bad":
        kw["pregate_bad"] = ["pregate.json's sha256 is not the one kept beside it"]
    elif case == "mask_seed":
        pregate["reading"]["groups"]["CM"]["mask_seed"] = 20261124
    elif case == "perm_seed_x019":
        pregate["reading"]["groups"]["CM"]["perm"]["X019"]["seed"] = 20261135
    elif case == "perm_seed_x020":
        del pregate["reading"]["groups"]["CM"]["perm"]["X020"]["seed"]
    elif case == "perm_members":
        del pregate["reading"]["groups"]["CM"]["perm"]["X020"]
    elif case == "group_index":
        pregate["reading"]["groups"]["CM"]["group_index"] = 1
    elif case == "group_unknown":
        pregate["reading"]["groups"]["Epicardium"] = pregate["reading"]["groups"]["CM"]
    r = mx.t72_reading(loaded, gate, checks, pregate, psha, **kw)
    assert r["outcome"] == "VOID" and r["reasons"] and r["follows"] == mx.FOLLOWS_VOID, case
    lines = mx.outcome_lines(r)
    assert lines[0].startswith("Outcome: VOID (") and lines[-1] == mx.FOLLOWS_VOID
    assert not any(b in lines for b in mx.BESIDE)


def test_the_gate_names_its_largest_difference_and_holds_at_1e_9(mx):
    loaded = _loaded(mx)
    loaded["2"]["rows"][4]["arms"][mx.BASE]["score"] += 1e-9 * 0.5
    assert _read(mx, loaded)["outcome"] != "VOID"
    loaded["2"]["rows"][4]["arms"][mx.BASE]["score"] += 0.25
    r = _read(mx, loaded)
    assert r["outcome"] == "VOID" and any("the gate failed" in x and "0.25" in x for x in r["reasons"])
    assert mx.REPRO_TOL == 1e-9 != mx.br.REPRO_TOL


def test_the_counts_and_constants(mx):
    assert mx.TRIES_COUNT == 15 and mx.COUNTED == ("mx_cons", "mx_mask", "mx_perm") and mx.BASE == "m1rc_gmom_x2"
    assert mx.LIB_ARMS[0] == "copy_last" and mx.ROW == "T72"
    assert (mx.MULTIPLIER, mx.PRACTICE_FRACTION, mx.TEST_FRACTION, mx.FRACTION_TOL) == (2.0, 0.2, 0.4, 1e-12)
    assert (mx.MASK_SEED, mx.PERM_SEED, mx.BINS) == (20261123, 20261124, 10)
    assert mx.MIN_CELLS == {"CM": 200, "Endocardium": 50, "Endothelium": 50}
    assert (mx.SKM_SCORE, mx.MAX_SKM_SHARE, mx.SPEARMAN_LOW, mx.SPEARMAN_HIGH, mx.X002_MIN_CP) == (
        0.3, 0.10, 0.2, 0.95, 1.0)
    assert (mx.ZERO_LOW, mx.ZERO_HIGH, mx.MAX_CONS_SPEARMAN, mx.MAX_HELD_SHARE) == (0.05, 0.60, 0.98, 0.5)
    assert (mx.CHAMBER_BAR, mx.LIB_BAR, mx.LIB_FLOOR, mx.CTRL_CHAMBER_BAR, mx.CTRL_LIB_BAR) == (
        0.3, 0.3, -0.3, 0.2, 0.2)
    assert (mx.DIRECTION_SCALE, mx.PSEUDOCOUNT, mx.MIN_TABLE_CELLS) == (5 / 8, 0.5, 10)
    assert mx.CELL_CALLS == {"GSE230531": {"min_umi": 1000, "min_genes": 500, "max_mito": 0.20},
                             "GSE244783": {"min_umi": 500, "min_genes": 300, "max_mito": 0.05}}
    assert mx.GROUP_MEMBERS == {"CM": ("X002", "X019", "X020"), "Endocardium": ("X002", "X019"),
                                "Endothelium": ("X002", "X019")}
    assert mx.T63_TASKS == SCRIPT.parent / "results" / "t63_double_step" / "tasks"
    assert all((mx.T63_TASKS / f).exists() for f in mx.task_files().values())
    assert set(mx.NEVER_FETCHED) >= {"GSM7829701", "GSE244783_RAW.tar", "GSE230531_RAW.tar", "GSM7226279"}
    assert mx.ATLASES_YAML_SHA256 == YSHA  # the fixture's; the module leaves it to the placement fill


def test_every_output_is_created_exclusively_and_the_source_guard(mx, tmp_path):
    p = tmp_path / "lib1.json"
    mx.write_new(p, "first")
    with pytest.raises(FileExistsError):
        mx.write_new(p, "second")
    text = SCRIPT.read_text()
    assert ".write_text(" not in text and "os.mkdir(CONSENSUS)" in text
    assert mx.source_guard() == [] and mx.rc.source_guard(SCRIPT) == [] and mx.all_source_bad() == []
    for line in ("import double_step_read\n", "import perm_read\n", "import " + "n5" + "4_relabel_listing\n",
                 "import " + "n4" + "8_m1_direction_shrink\n", "import jaxblood_read\n"):
        bad = tmp_path / "bad.py"
        bad.write_text(line)
        assert mx.source_guard(bad), line


def test_no_t2_or_t3_path_or_the_shared_core_imports_this_read_or_its_tables():
    pat = re.compile(r"mx2s_read|t1_ext_t72|mx2s\.sh")
    roots = [ROOT / "src" / "vec", *sorted((ROOT / "experiments").glob("t2*")),
             *sorted((ROOT / "experiments").glob("t3*"))]
    hits = [str(p) for r in roots for p in r.rglob("*.py") if pat.search(p.read_text(errors="ignore"))]
    assert hits == [] and (ROOT / "src" / "vec").is_dir()


def test_the_read_prints_only_words_writes_the_outcome_line_and_never_overwrites(mx, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mx, "TRIES", _tries(tmp_path))
    monkeypatch.setattr(mx.gv, "tree_state", lambda: {"commit": C40, "clean": True})
    monkeypatch.setattr(mx.gv, "stamp", lambda started, fs: {"commit": C40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(mx.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    monkeypatch.setattr(mx, "yaml_problems", lambda *a: [])
    monkeypatch.setattr(mx, "arm_table_problems", lambda pg: [])
    out = mx.RESULTS / "t72"
    tasks = out / "tasks"
    tasks.mkdir(parents=True)
    text = json.dumps(_pregate(mx))
    (tasks / "pregate.json").write_text(text)
    sha = hashlib.sha256(text.encode()).hexdigest()
    (tasks / "pregate.json.sha256").write_text(sha + "  pregate.json\n")
    loaded = _loaded(mx, gain=1.234567)
    for k, f in mx.task_files().items():
        loaded[k]["pregate_sha256"] = sha
        (tasks / f).write_text(json.dumps(loaded[k]))
        (tasks / mx.started_name(f)).write_text("{}")
    gate = _gate()
    files = {v: k for k, v in mx.task_files().items()}
    monkeypatch.setattr(mx, "committed", lambda f, *a: gate[files[f]])
    res = mx.read(out, C40)
    assert res["reading"]["outcome"] == "AHEAD" and res["tries"] == 15
    printed = capsys.readouterr().out
    assert printed.splitlines()[0] == mx.AHEAD_WORDS and not re.search(FIGURE, printed)
    txt = (out / "t72.txt").read_text()
    assert "\nOUTCOME AHEAD\n" in txt and "+1.23" in txt
    assert json.loads((out / "t72.json").read_text())["reading"]["outcome"] == "AHEAD"
    with pytest.raises(SystemExit, match="already holds"):
        mx.read(out, C40)


def _read_ready(mx, tmp_path, monkeypatch, pregate=True, passes=True):
    monkeypatch.setattr(mx, "TRIES", _tries(tmp_path))
    monkeypatch.setattr(mx.gv, "tree_state", lambda: {"commit": C40, "clean": True})
    monkeypatch.setattr(mx.gv, "stamp", lambda started, fs: {"commit": C40, "tree_clean": True, "utc": "x"})
    monkeypatch.setattr(mx.rc, "relabel_checks", lambda *a, **k: {"passes": True, "reasons": []})
    monkeypatch.setattr(mx, "committed", lambda f, *a: {s: 51.0 for s in range(10)})
    out = mx.RESULTS / "t72"
    (out / "tasks").mkdir(parents=True)
    if pregate:
        text = json.dumps(_pregate(mx, passes=passes))
        (out / "tasks" / "pregate.json").write_text(text)
        (out / "tasks" / "pregate.json.sha256").write_text(hashlib.sha256(text.encode()).hexdigest() + "  x\n")
    return out


def test_the_read_says_void_when_a_start_record_has_no_task_file(mx, tmp_path, monkeypatch, capsys):
    out = _read_ready(mx, tmp_path, monkeypatch)
    (out / "tasks" / "lib2.started").write_text("{}")
    res = mx.read(out, C40)
    assert res["reading"]["outcome"] == "VOID"
    assert any("2: no task file (its start record exists)" in x for x in res["reading"]["reasons"])
    assert capsys.readouterr().out.startswith("Outcome: VOID (")
    assert "\nOUTCOME VOID\n" in (out / "t72.txt").read_text()


@pytest.mark.parametrize("case", ["pregate_missing", "pregate_failed", "no_start_record"])
def test_the_read_refuses_and_writes_nothing_without_a_passed_pregate_or_a_start_record(mx, tmp_path, monkeypatch,
                                                                                         capsys, case):
    out = _read_ready(mx, tmp_path, monkeypatch, pregate=case != "pregate_missing", passes=case != "pregate_failed")
    if case != "no_start_record":
        (out / "tasks" / "T1_E9.5_E10.5_chamber.started").write_text("{}")
    with pytest.raises(SystemExit) as e:
        mx.read(out, C40)
    assert e.value.code not in (0, None) and "nothing is written" in str(e.value.code)
    assert not (out / "t72.json").exists() and not (out / "t72.txt").exists()
    assert "the tries count" not in capsys.readouterr().out


# ---------------------------------------------------------------------------------------------------------------------
# The launcher, traced under a temporary HOME in a temporary repository


INPUTS = ("data/derived/feng_relabel_c/cd1_stats_relabel.npz", "data/raw/T1/E8.5_RNA.h5ad",
          "data/raw/T1/E9.5_RNA.h5ad", "data/raw/T1/E10.5_RNA.h5ad")
COMMITTED_INPUTS = ("experiments/t1_ext/atlases.yaml", "experiments/t1_ext/results/extra_competing.json")


def _git_env(home):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(home), GIT_CEILING_DIRECTORIES=str(home.parent), GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
    return env


def _repo(tmp_path, row="T72", copies=1, inputs=True):
    repo, home = tmp_path / "repo", tmp_path / "home"
    (repo / "experiments" / "t1_m1").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    home.mkdir()
    shutil.copy(LAUNCHER, repo / "experiments" / "t1_m1" / "mx2s.sh")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(
        "".join(f"| {row} (one row) | T1 | x |\n" for _ in range(copies)))
    (repo / ".gitignore").write_text("data/\n")
    if inputs:
        for f in COMMITTED_INPUTS:
            (repo / f).parent.mkdir(parents=True, exist_ok=True)
            (repo / f).write_text("x\n")
    env = _git_env(home)
    for cmd in (["init", "-q"], ["add", "-A"], ["commit", "-q", "-m", "row"]):
        subprocess.run(["git", *cmd], cwd=repo, env=env, check=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, check=True, capture_output=True,
                          text=True).stdout.strip()
    _fake_uv(home, code=7)  # nothing real ever runs from a launcher test
    return repo, home, env, head


def _data(repo):
    for f in INPUTS:
        (repo / f).parent.mkdir(parents=True, exist_ok=True)
        (repo / f).write_text("x")


def _run(repo, env, *args, cwd=None):
    r = subprocess.run(["bash", str(repo / "experiments" / "t1_m1" / "mx2s.sh"), *args], cwd=cwd or repo, env=env,
                       capture_output=True, text=True)
    return r.returncode


def _s(home):
    return home / "v" / "scratch" / "t72"


def _refused(home, name):
    return (_s(home) / f"{name}.refused").read_text().strip()


def _no_record(home, name):
    return not (_s(home) / f"{name}.exit").exists() and not (_s(home) / f"{name}.head").exists()


def _passed(home, head, *names, code="0"):
    _s(home).mkdir(parents=True, exist_ok=True)
    for n in names:
        (_s(home) / f"{n}.exit").write_text(code + "\n")
        (_s(home) / f"{n}.head").write_text(head + "\n")


def _fake_uv(home, code=0):
    """A fake uv in $HOME/.local/bin (the launcher puts it first on PATH): it logs its arguments and exits `code`."""
    b = home / ".local" / "bin"
    b.mkdir(parents=True, exist_ok=True)
    (b / "uv").write_text(f'#!/usr/bin/env bash\necho "$*" >> "$HOME/calls.txt"\nexit {code}\n')
    (b / "uv").chmod(0o755)
    return home / "calls.txt"


def test_the_launcher_refusals_in_the_rows_order(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    assert _run(repo, env) == 2 and _refused(home, "usage") == "2"
    assert _run(repo, env, "cut", head) == 2
    assert _run(repo, env, "bogus", head) == 97 and _refused(home, "bogus") == "97"
    assert _run(repo, env, "tests", "f" * 40) == 95 and _no_record(home, "tests")
    assert _run(repo, env, "tests", head[:12]) == 95 and _refused(home, "tests") == "95"
    assert _run(repo, env, "tests", head.upper()) == 95
    assert _run(repo, env, "tests", head) == 93  # the released stages and N54's file are not in place
    _data(repo)
    order = ["tests", "fetch", "load", "tables", "pregate", "cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3",
             "read"]
    args = {"cut_chamber": ("cut", "chamber"), "cut_whole": ("cut", "whole"), "lib_1": ("lib", "1"),
            "lib_2": ("lib", "2"), "lib_3": ("lib", "3")}
    for prev, step in zip(order, order[1:], strict=False):
        assert _run(repo, env, *args.get(step, (step,)), head) == 94 and _no_record(home, step), step
        _passed(home, head, prev, code="3")
        assert _run(repo, env, *args.get(step, (step,)), head) == 94, step  # exit 3 is not 0
        _passed(home, head, prev)
    _passed(home, head, *order)
    for step in order:  # every step that ran once is refused a second call
        assert _run(repo, env, *args.get(step, (step,)), head) == 91 and _refused(home, step) == "91", step
    (_s(home) / "cut_chamber.head").write_text("0" * 40 + "\n")  # exit 0 at another commit is not passed
    assert _run(repo, env, "lib", "2", head) == 91  # lib_1 has passed; lib_2 ran
    (repo / "stray.txt").write_text("x")
    assert _run(repo, env, "read", head) == 98 and _refused(home, "read") == "98"


def test_the_data_steps_and_the_pregate_refuse_on_their_outputs_and_inputs(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    _passed(home, head, "tests")
    (repo / "data" / "external" / "t1_ext_t72").mkdir(parents=True)
    assert _run(repo, env, "fetch", head) == 91 and _no_record(home, "fetch")
    _passed(home, head, "fetch")
    (repo / "data" / "external" / "t1_ext_t72").rename(repo / "data" / "external" / "moved")
    assert _run(repo, env, "load", head) == 93
    _passed(home, head, "load")
    (repo / "data" / "derived" / "t1_ext_t72" / "tables").mkdir(parents=True)
    assert _run(repo, env, "tables", head) == 91
    _passed(home, head, "tables")
    assert _run(repo, env, "pregate", head) == 93  # no tables manifest
    (repo / "data" / "derived" / "t1_ext_t72" / "tables" / "MANIFEST.json").write_text("{}")
    (repo / "data" / "derived" / "t1_ext_t72" / "consensus").mkdir()
    assert _run(repo, env, "pregate", head) == 91 and _no_record(home, "pregate")
    (repo / "experiments" / "t1_ext" / "atlases.yaml").write_text("changed\n")
    assert _run(repo, env, "pregate", head) == 98  # a tracked input changed: the tree is not clean


def test_the_launcher_refuses_without_exactly_one_row(tmp_path):
    repo, home, env, head = _repo(tmp_path, row="T71")
    assert _run(repo, env, "tests", head) == 96 and _refused(home, "tests") == "96"
    repo2, home2, env2, head2 = _repo(tmp_path / "two", copies=2)
    assert _run(repo2, env2, "tests", head2) == 96


def test_a_committed_input_missing_refuses(tmp_path):
    repo, home, env, head = _repo(tmp_path, inputs=False)
    _data(repo)
    assert _run(repo, env, "tests", head) == 93


def test_a_step_that_runs_marks_itself_running_then_writes_its_exit_with_uv_from_local_bin(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    b = home / ".local" / "bin"
    seen = tmp_path / "seen.txt"
    (b / "uv").write_text(f'#!/usr/bin/env bash\ncat "$HOME/v/scratch/t72/tests.exit" "$HOME/v/scratch/t72/'
                          f'tests.head" > "{seen}"\npwd >> "{seen}"\necho "$TMPDIR $NUMPY_MADVISE_HUGEPAGE '
                          f'$OMP_NUM_THREADS" >> "{seen}"\necho "$*" >> "{seen}"\nexit 4\n')
    (b / "uv").chmod(0o755)
    assert _run(repo, env, "tests", head, cwd=repo / "experiments" / "t1_m1") == 4
    lines = seen.read_text().splitlines()
    assert lines[:3] == ["running", head, str(repo)]
    assert lines[3] == f"{_s(home)}/tmp 0 8" and "tests/test_mx2s_read.py" in lines[4]
    assert (_s(home) / "tests.exit").read_text().strip() == "4"
    assert _run(repo, env, "tests", head) == 91  # a crashed or failed step is not run again


def test_a_run_step_refuses_on_its_task_file_or_start_record(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    calls = _fake_uv(home)
    _passed(home, head, "tests", "fetch", "load", "tables", "pregate")
    tasks = repo / "experiments" / "t1_m1" / "results" / "t72_mx2s" / "tasks"
    tasks.mkdir(parents=True)
    (tasks / "T1_E9.5_E10.5_chamber.started").write_text("{}")
    assert _run(repo, env, "cut", "chamber", head) == 91 and _no_record(home, "cut_chamber")
    (tasks / "T1_E9.5_E10.5_chamber.started").rename(tasks / "moved")
    assert _run(repo, env, "cut", "chamber", head) == 0
    assert calls.read_text().split("\n")[0].endswith(f"run-cut --cut chamber --out t72_mx2s --commit {head}")


def _runs_ready(tmp_path, chamber="1"):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    calls = _fake_uv(home)
    _passed(home, head, "tests", "fetch", "load", "tables", "pregate")
    _passed(home, head, "cut_chamber", code=chamber)
    tasks = repo / "experiments" / "t1_m1" / "results" / "t72_mx2s" / "tasks"
    tasks.mkdir(parents=True)
    return repo, home, env, head, calls, tasks


def test_the_runs_entry_restarts_once_after_a_refused_cut_and_runs_the_five_steps_in_order(tmp_path):
    repo, home, env, head, calls, _ = _runs_ready(tmp_path)
    assert _run(repo, env, "runs", head) == 0
    got = [ln.split("mx2s_read.py ")[1].split(" --out")[0] for ln in calls.read_text().splitlines()]
    assert got == ["check-pregate", "run-cut --cut chamber", "run-cut --cut whole", "run-lib --library 1",
                   "run-lib --library 2", "run-lib --library 3"]
    for n in ("cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3"):
        assert (_s(home) / f"runs_{n}.exit").read_text().strip() == "0"
    assert (_s(home) / "runs.exit").read_text().strip() == "0"
    assert _run(repo, env, "runs", head) == 91  # once in all
    assert _run(repo, env, "read", head) == 0  # the read follows the runs entry's lib 3


@pytest.mark.parametrize("case", ["chamber_passed", "chamber_running", "chamber_never", "start_record",
                                  "task_file", "read_ran", "check_fails", "other_start"])
def test_the_runs_entry_refuses_outside_the_rows_conditions(tmp_path, case):
    repo, home, env, head, calls, tasks = _runs_ready(tmp_path, chamber={"chamber_passed": "0",
                                                                         "chamber_running": "running"}.get(case, "1"))
    if case == "chamber_never":
        (_s(home) / "cut_chamber.exit").rename(_s(home) / "moved")
    elif case == "start_record":
        (tasks / "lib1.started").write_text("{}")
    elif case == "other_start":
        (tasks / "anything.started").write_text("{}")
    elif case == "task_file":
        (tasks / "T1_E9.5_E10.5.json").write_text("{}")
    elif case == "read_ran":
        _passed(home, head, "read", code="1")
    elif case == "check_fails":
        _fake_uv(home, code=1)
    assert _run(repo, env, "runs", head) == 92 and _no_record(home, "runs"), case


def test_the_launcher_runs_the_tests_first_sets_path_and_never_prints_a_log_or_removes_a_file():
    sh = LAUNCHER.read_text()
    tests = sh.split("  tests)\n")[1].split(";;")[0]
    assert "$PY -m pytest" in tests
    for f in ("test_mx2s_read.py", "test_cmjax_read.py", "test_jax_cm.py", "test_t62_recomputed_read.py",
              "test_t63_double_step_read.py", "test_core.py", "test_t1_ext_load.py", "test_t1_ext_stage.py",
              "test_t1_ext_tables.py"):
        assert f"tests/{f}" in tests
    assert all((ROOT / "tests" / f).exists() for f in ("test_mx2s_read.py", "test_cmjax_read.py", "test_core.py"))
    assert sh.index('export PATH="$HOME/.local/bin:$PATH"') < sh.index("uv run")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh and "OMP_NUM_THREADS=8" in sh
    runs = [ln for ln in sh.splitlines() if "$PY" in ln and ".py" in ln and "check-pregate" not in ln
            and "runs_check" not in ln] + [tests]
    assert all('> "$S/$NAME.out" 2> "$S/$NAME.err"' in r or '> "$S/$N.out" 2> "$S/$N.err"' in r for r in runs), runs
    assert not re.search(r"\brm\b|rmdir|unlink|:\?", sh)
    assert "/v/vec-campaign" not in sh and "cat \"$S/$NAME.out\"" not in sh
    assert subprocess.run(["bash", "-n", str(LAUNCHER)]).returncode == 0


# ---------------------------------------------------------------------------------------------------------------------
# The read after a stop (blocker 1), the clean-tree exemption and the start records (row v3)


def _tasks_dir(repo):
    t = repo / "experiments" / "t1_m1" / "results" / "t72_mx2s" / "tasks"
    t.mkdir(parents=True, exist_ok=True)
    return t


def test_the_launchers_task_files_are_the_scripts(mx):
    sh = LAUNCHER.read_text()
    names = re.search(r'TASK_FILES="([^"]+)"', sh).group(1).split()
    assert sorted(names) == sorted(mx.task_files().values())
    assert {mx.CUT_FILES[c] for c in ("chamber", "whole")} == {"T1_E9.5_E10.5_chamber.json", "T1_E9.5_E10.5.json"}


def test_a_stopped_second_task_step_with_the_first_steps_start_record_lets_the_read_run_and_it_says_void(
        mx, tmp_path, monkeypatch, capsys):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    calls = _fake_uv(home)
    _passed(home, head, "tests", "fetch", "load", "tables", "pregate", "cut_chamber")
    _passed(home, head, "cut_whole", code="1")  # cut whole stopped after cut chamber loaded and wrote its file
    tasks = _tasks_dir(repo)
    (tasks / "T1_E9.5_E10.5_chamber.started").write_text("{}")
    (tasks / "T1_E9.5_E10.5_chamber.json").write_text("{}")
    for step in (("cut", "whole"), ("lib", "1")):  # the chain does not go on
        assert _run(repo, env, *step, head) in (91, 94)
    assert _run(repo, env, "read", head) == 0
    assert calls.read_text().splitlines()[-1].endswith(f"mx2s_read.py read --out t72_mx2s --commit {head}")
    assert (_s(home) / "read.exit").read_text().strip() == "0"
    assert _run(repo, env, "read", head) == 91  # once
    # the script's read on the same state: the row's one VOID, every missing task named
    out = _read_ready(mx, tmp_path, monkeypatch)
    (out / "tasks" / "T1_E9.5_E10.5_chamber.started").write_text("{}")
    (out / "tasks" / "T1_E9.5_E10.5_chamber.json").write_text(json.dumps(_task(mx)))
    res = mx.read(out, C40)
    assert res["reading"]["outcome"] == "VOID"
    for k in ("whole", "1", "2", "3"):
        assert f"{k}: no task file" in res["reading"]["reasons"]
    printed = capsys.readouterr().out.splitlines()
    assert printed[0].startswith("Outcome: VOID (") and printed[-1] == mx.FOLLOWS_VOID


def test_a_stop_inside_the_runs_entry_lets_the_read_run(tmp_path):
    repo, home, env, head, calls, tasks = _runs_ready(tmp_path)  # cut chamber refused before it loaded (exit 1)
    _passed(home, head, "runs_cut_chamber")
    _passed(home, head, "runs_cut_whole", "runs", code="1")
    (tasks / "T1_E9.5_E10.5_chamber.started").write_text("{}")
    (tasks / "T1_E9.5_E10.5.started").write_text("{}")
    assert _run(repo, env, "read", head) == 0


@pytest.mark.parametrize("case", ["pregate_failed", "pregate_never", "pregate_other_commit", "no_start_record",
                                  "running", "runs_running", "stop_at_other_commit", "no_stop"])
def test_the_read_is_refused_outside_the_chains_end_or_stop(tmp_path, case):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    calls = _fake_uv(home)
    _passed(home, head, "tests", "fetch", "load", "tables")
    _passed(home, head, "pregate", code={"pregate_failed": "3"}.get(case, "0"))
    if case == "pregate_never":
        (_s(home) / "pregate.exit").rename(_s(home) / "moved")
    if case == "pregate_other_commit":
        (_s(home) / "pregate.head").write_text("0" * 40 + "\n")
    _passed(home, head, "cut_chamber")
    _passed(home, head, "cut_whole", code={"running": "running", "no_stop": "0"}.get(case, "1"))
    if case == "stop_at_other_commit":
        (_s(home) / "cut_whole.head").write_text("0" * 40 + "\n")
    if case == "runs_running":
        _passed(home, head, "runs_lib_1", code="running")
    if case != "no_start_record":
        (_tasks_dir(repo) / "T1_E9.5_E10.5_chamber.started").write_text("{}")
    assert _run(repo, env, "read", head) == 94 and _refused(home, "read") == "94", case
    assert _no_record(home, "read") and not calls.exists()


def test_the_read_after_lib_3_still_needs_the_pregate(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    _fake_uv(home)
    _passed(home, head, "tests", "fetch", "load", "tables", "cut_chamber", "cut_whole", "lib_1", "lib_2", "lib_3")
    _passed(home, head, "pregate", code="3")
    assert _run(repo, env, "read", head) == 94
    _passed(home, head, "pregate")
    assert _run(repo, env, "read", head) == 0


@pytest.mark.parametrize("where,code", [
    ("experiments/t1_m1/results/t72_mx2s/tasks/lib1.started", 94),
    ("experiments/t1_m1/results/t72_mx2s/tasks/pregate.json", 94),
    ("experiments/t1_m1/results/t72_mx2s/t72.txt", 94),
    ("experiments/t1_m1/results/other/x.json", 98),
    ("experiments/t1_m1/results/t72_mx2s_other/x.json", 98),
    ("experiments/t1_m1/x.started", 98),
    ("stray.txt", 98)])
def test_the_clean_check_exempts_only_untracked_files_of_the_rows_results_folder(tmp_path, where, code):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    _passed(home, head, "tests")
    (repo / where).parent.mkdir(parents=True, exist_ok=True)
    (repo / where).write_text("x")
    assert _run(repo, env, "load", head) == code  # 94: the tree passed and fetch has not run; 98: unclean
    assert _refused(home, "load") == str(code)


def test_a_changed_tracked_file_in_the_rows_results_folder_still_refuses(tmp_path):
    repo, home, env, head = _repo(tmp_path)
    _data(repo)
    f = _tasks_dir(repo) / "kept.txt"
    f.write_text("x")
    for cmd in (["add", "-A"], ["commit", "-q", "-m", "tracked"]):
        subprocess.run(["git", *cmd], cwd=repo, env=env, check=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, env=env, check=True, capture_output=True,
                          text=True).stdout.strip()
    f.write_text("changed")
    assert _run(repo, env, "tests", head) == 98


def test_nothing_deletes_or_moves_a_start_record():
    for path in (SCRIPT, LAUNCHER):
        text = path.read_text()
        for token in (".unlink(", "os.remove(", "rmtree", ".rmdir(", ".rename(", "os.replace(", "shutil.move",
                      "os.truncate"):
            assert token not in text, (path, token)
    sh = LAUNCHER.read_text()
    assert not re.search(r"\b(rm|mv|rmdir|unlink|truncate)\b", sh)
    assert "> \"$TASKS" not in sh  # the launcher never writes into the tasks folder