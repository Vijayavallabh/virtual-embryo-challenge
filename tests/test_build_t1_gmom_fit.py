"""--g-mom-fit of experiments/builds/build_test.py t1 (off by default, no entry): g_mom's test-time factors fitted on
E10.5's dissected-heart library (names ending in _3) in place of the whole reference. Synthetic stages only: every
input of build_t1 is monkeypatched as in tests/test_hidden_counts.py, and no released file, board file or Feng
statistic is read."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from test_hidden_counts import FENG_GROUPS, T1_SECRET, bt, t1_env  # noqa: F401

from vec import t1_scope
from vec import t1_state_transport as ts
from vec.boards import Board

SEED = 1_000_001
# the g_mom record's keys as build_t1 wrote them before --g-mom-fit (main cfa7e3f), and the fields only lib3 adds
BEFORE_KEYS = {"rule", "edit_categories", "categories_with_factors", "categories_without_factors", "stream",
               "held_flat", "median_depth_scale"}
LIB3_KEYS = {"fit", "fit_reference", "metacells", "reading", "edited_share_pool", "edited_share"}


def _libraries(a, lib_of):
    """A copy of stage `a` whose cell names end in _<library>, lib_of(i, label) giving each cell's library."""
    out = a.copy()
    labels = out.obs["celltype"].astype(str).to_numpy()
    out.obs_names = [f"c{i}_{lib_of(i, lab)}" for i, lab in enumerate(labels)]
    return out


def _alternate(i, lab):  # half of every label in library 3, the rest in libraries 1 and 2
    return 3 if i % 2 == 0 else (1 if i % 4 == 1 else 2)


def _loads(bt, env, monkeypatch):  # noqa: F811
    """E10.5 and E9.5 as build_t1 reads them; `returned` holds each object handed out, by file name."""
    returned = {}

    def load(path, **kw):
        name = Path(path).name
        obj = (env.ref if name.startswith("E10.5") else env.start).copy()
        returned.setdefault(name, []).append(obj)
        return obj

    monkeypatch.setattr(bt, "load_stage", load)
    return returned


def _fit_spy(bt, monkeypatch):  # noqa: F811
    calls = []
    real = bt.harness.g_mom_pair_factors

    def spy(ref, src, start, start_src, **kw):
        calls.append({"ref": ref, "ref_obs": ref.obs.copy(), "src": src, "start": start, "start_src": start_src, **kw})
        out = real(ref, src, start, start_src, **kw)
        calls[-1]["factors"] = out[0]  # the fit as returned, kept beside the call (not a key the call compares)
        return out

    monkeypatch.setattr(bt.harness, "g_mom_pair_factors", spy)
    return calls


def _counts(a):
    return sparse.csr_matrix(a.layers["counts"]).toarray()


def test_the_option_defaults_to_whole_and_needs_g_mom(bt):  # noqa: F811
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "runs/x/T1_s1.h5ad"]).g_mom_fit == "whole"
    assert ap.parse_args(["t1", "--g-mom", "--g-mom-fit", "lib3", "--out", "runs/x/T1_s1.h5ad"]).g_mom_fit == "lib3"
    with pytest.raises(SystemExit):
        ap.parse_args(["t1", "--g-mom-fit", "other", "--out", "runs/x/T1_s1.h5ad"])


def test_lib3_without_g_mom_is_refused_before_anything_is_read(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    monkeypatch.setattr(bt, "load_stage", lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    out = tmp_path / "T1_s1.h5ad"
    with pytest.raises(SystemExit, match="--g-mom-fit lib3 is refused without --g-mom"):
        bt.main(["t1", "--g-mom-fit", "lib3", "--out", str(out)])
    assert not out.exists() and t1_env.reads == []


def test_default_fits_on_the_whole_reference_object_and_writes_the_file_as_before(bt, t1_env, tmp_path,  # noqa: F811
                                                                                  monkeypatch):
    """The default path hands g_mom_pair_factors the very E10.5 object build_t1 loaded, and the file is byte for byte
    the one the build path before the option writes (build_submission with the hook of _t1_g_mom on the whole
    reference, called as before the option); an explicit --g-mom-fit whole writes the same bytes. The record's g_mom
    field holds exactly the keys it held before the option: the default adds no field."""
    t1_env.ref = _libraries(t1_env.ref, _alternate)
    returned = _loads(bt, t1_env, monkeypatch)
    calls = _fit_spy(bt, monkeypatch)
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--out", str(out)])
    assert len(calls) == 1 and calls[0]["ref"] is returned["E10.5_RNA.h5ad"][0]
    rec = json.loads(out.with_suffix(".build.json").read_text())["g_mom"]
    assert set(rec) == BEFORE_KEYS  # the record as on main before the option
    assert "n_cells" not in json.dumps(rec) and str(T1_SECRET) not in json.dumps(rec)

    explicit = tmp_path / "T1_s1_whole.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--g-mom-fit", "whole", "--out", str(explicit)])
    assert explicit.read_bytes() == out.read_bytes()
    assert json.loads(explicit.with_suffix(".build.json").read_text())["g_mom"] == rec

    # the build path as it was before the option, step by step
    ref = t1_env.ref.copy()
    shares = bt.harness.t1_chamber_shares(bt.harness._shallow(ref), "E10.5_RNA")
    cats = [c for c, s in shares.items() if s > 0]
    groups = bt.feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique())
    covered = {lab for lab, g in groups.items() if g in t1_env.dirs}
    hook, _ = bt._t1_g_mom(ref, cats, SEED, covered)
    sub, _ = t1_env.build(ref, t1_env.board, shares, bt.harness.count_for(t1_env.board, ref.n_obs), SEED,
                          held_flat=["held_a", "held_b"], ref_source="E10.5_RNA", directions=t1_env.dirs,
                          direction_fraction=bt.feng_dirs.bracket_fraction(10.5, 12.5), direction_groups=groups,
                          herd=False, post_cells=hook)
    before = tmp_path / "T1_before.h5ad"
    bt.t1_scope.write_checked(sub, before, t1_env.board)
    assert before.read_bytes() == out.read_bytes()


def test_lib3_fit_receives_exactly_the_library3_rows_and_nothing_else_differs(bt, t1_env, tmp_path,  # noqa: F811
                                                                              monkeypatch):
    t1_env.ref = _libraries(t1_env.ref, _alternate)
    returned = _loads(bt, t1_env, monkeypatch)
    calls = _fit_spy(bt, monkeypatch)
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--out", str(tmp_path / "T1_whole.h5ad")])
    out = tmp_path / "T1_lib3.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--g-mom-fit", "lib3", "--out", str(out)])
    whole, lib3 = calls
    loaded = returned["E10.5_RNA.h5ad"][1]
    in3 = np.array([str(n).endswith("_3") for n in loaded.obs_names])
    assert 0 < in3.sum() < loaded.n_obs
    assert list(lib3["ref"].obs_names) == list(loaded.obs_names[in3])  # library 3's rows, in file order, and no other
    assert (_counts(lib3["ref"]) == _counts(loaded)[in3]).all()
    assert list(lib3["ref"].var_names) == list(loaded.var_names)
    assert list(lib3["ref_obs"].columns) == list(whole["ref_obs"].columns)  # the obs as the fit received it
    assert lib3["ref_obs"].equals(whole["ref_obs"][in3])
    for k in ("src", "start_src", "categories"):  # the same start stage, categories rule and stems
        assert lib3[k] == whole[k]
    assert set(lib3) == set(whole)
    assert list(lib3["start"].obs_names) == list(whole["start"].obs_names)
    assert (_counts(lib3["start"]) == _counts(whole["start"])).all()
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["g_mom"]["fit"] == "lib3" and "library-3" in rec["g_mom"]["fit_reference"]
    assert rec["base"] == "pooled" and rec["g_mom"]["stream"] == [SEED, 5]
    assert set(rec["g_mom"]) == (BEFORE_KEYS | LIB3_KEYS) - {"edited_share"}  # edited_share: stripped on a test board
    # the record's metacells are the fit's metacell counts (never cell numbers), each within METACELL_MAX
    factors = lib3["factors"]
    assert rec["g_mom"]["metacells"] == {c: factors.categories[c].n_metacells for c in factors.categories}
    assert all(0 < v <= ts.METACELL_MAX for v in rec["g_mom"]["metacells"].values())
    assert rec["g_mom"]["reading"] == bt.G_MOM_LIB3_READING
    for words in ("size-matched", "zero by construction", "is not the same step without the collection part",
                  "unknown in sign"):
        assert words in rec["g_mom"]["reading"]


def test_on_a_reference_all_of_library_3_the_two_fits_give_the_same_file(bt, t1_env, tmp_path,  # noqa: F811
                                                                         monkeypatch):
    """The idle case: every reference cell is of library 3, so the lib3 fit's reference is the whole one's rows."""
    t1_env.ref = _libraries(t1_env.ref, lambda i, lab: 3)
    _loads(bt, t1_env, monkeypatch)
    a, b = tmp_path / "T1_whole.h5ad", tmp_path / "T1_lib3.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--out", str(a)])
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--g-mom-fit", "lib3", "--out", str(b)])
    assert a.read_bytes() == b.read_bytes()
    ra = json.loads(a.with_suffix(".build.json").read_text())
    rb = json.loads(b.with_suffix(".build.json").read_text())
    assert set(ra["g_mom"]) == BEFORE_KEYS  # the default adds no field
    assert ra["g_mom"] == {k: v for k, v in rb["g_mom"].items() if k not in LIB3_KEYS}
    assert {k: v for k, v in ra.items() if k != "g_mom"} == {k: v for k, v in rb.items() if k != "g_mom"}
    # every pool cell sits in the fit, so the lib3 record's pool share is R1's unheld share of each category
    pool = t1_env.ref.copy()
    t1_scope.classify(pool, source="E10.5_RNA")
    cat = pool.obs["category"].astype(str).to_numpy()
    lab = pool.obs["celltype"].astype(str).to_numpy()
    held = np.isin(lab, sorted(ts.r1_unbracketed_labels(np.unique(lab))))
    for c in bt.G_MOM_EDIT_CATEGORIES:
        assert rb["g_mom"]["edited_share_pool"][c] == pytest.approx((~held[cat == c]).mean())


def test_cells_whose_source_is_outside_the_fit_keep_their_counts_and_enter_the_share(bt, t1_env,  # noqa: F811
                                                                                     monkeypatch):
    """The lib3 fit on a pooled base: a drawn cell from library 1 or 2 sits in no fitted metacell, keeps every count,
    and the record's shares (drawn and pool) count it among R3b's cells that the step does not edit; the pool's "all"
    weights the two categories by the build's shares, not by their numbers in the pool."""
    ref = _libraries(t1_env.ref, _alternate)
    monkeypatch.setattr(bt, "load_stage", lambda *a, **k: t1_env.start.copy())
    shares = bt.harness.t1_chamber_shares(bt.harness._shallow(ref), "E10.5_RNA")
    cats = [c for c, s in shares.items() if s > 0]
    covered = {lab for lab, g in bt.feng_dirs.label_groups(ref.obs["celltype"].astype(str).unique()).items()
               if g in FENG_GROUPS}
    in3 = np.array([str(n).endswith("_3") for n in ref.obs_names])
    build_shares = t1_scope.with_share(shares, "cm", 0.60)  # the build's target shares differ from the pool's mix
    hook, rec = bt._t1_g_mom(ref[np.flatnonzero(in3)].copy(), cats, SEED, covered, fit="lib3", pool=ref,
                             shares=build_shares)
    pool = ref.copy()
    t1_scope.classify(pool, source="E10.5_RNA")
    pop, _ = t1_scope.reweight(pool, {c: 1 / len(cats) for c in cats}, 400, 5,
                               fallback=t1_scope.load_spec().fallbacks("E10.5_RNA"))
    out = hook(pop)
    cat = pop.obs["category"].astype(str).to_numpy()
    labels = pop.obs["celltype"].astype(str).to_numpy()
    held = np.isin(labels, sorted(ts.r1_unbracketed_labels(np.unique(labels))))
    src3 = np.array([str(s).endswith("_3") for s in pop.obs["source_obs"]])
    scope = np.isin(cat, bt.G_MOM_EDIT_CATEGORIES)
    outside = scope & ~src3
    assert outside.any() and (scope & src3 & ~held).any()
    before, after = _counts(pop), _counts(out)
    assert (before[outside] == after[outside]).all()  # outside the fit: every count kept
    changed = (before != after).any(axis=1)
    assert changed[scope & src3 & ~held].mean() > 0.9 and not changed[outside].any()
    expected = (scope & src3 & ~held).sum() / scope.sum()
    assert rec["edited_share"]["all"] == pytest.approx(expected) and rec["edited_share"]["all"] < 1.0
    for c in bt.G_MOM_EDIT_CATEGORIES:
        m = cat == c
        assert rec["edited_share"][c] == pytest.approx((m & src3 & ~held).sum() / m.sum())
    # the pool's share: the base's cells of R3b's categories whose name sits in a fitted metacell, R1's labels apart
    pcat = pool.obs["category"].astype(str).to_numpy()
    plab = pool.obs["celltype"].astype(str).to_numpy()
    pheld = np.isin(plab, sorted(ts.r1_unbracketed_labels(np.unique(plab))))
    per = {c: ((pcat == c) & in3 & ~pheld).sum() / (pcat == c).sum() for c in bt.G_MOM_EDIT_CATEGORIES}
    for c in bt.G_MOM_EDIT_CATEGORIES:
        assert rec["edited_share_pool"][c] == pytest.approx(per[c])
    w = {c: build_shares[c] for c in bt.G_MOM_EDIT_CATEGORIES}
    expected_all = sum(w[c] * per[c] for c in w) / sum(w.values())
    assert rec["edited_share_pool"]["all"] == pytest.approx(expected_all)
    assert 0 < rec["edited_share_pool"]["all"] < 1
    pscope = np.isin(pcat, bt.G_MOM_EDIT_CATEGORIES)
    pool_weighted = (pscope & in3 & ~pheld).sum() / pscope.sum()
    assert abs(expected_all - pool_weighted) > 1e-6  # the case tells the two weightings apart


def test_a_test_board_record_keeps_the_pool_share_and_drops_the_drawn_share(bt, t1_env, tmp_path,  # noqa: F811
                                                                            monkeypatch):
    t1_env.ref = _libraries(t1_env.ref, _alternate)
    _loads(bt, t1_env, monkeypatch)
    out = tmp_path / "T1_s1.h5ad"
    bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--g-mom-fit", "lib3", "--out", str(out)])
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert "edited_share" not in rec["g_mom"] and 0 < rec["g_mom"]["edited_share_pool"]["all"] < 1
    assert str(T1_SECRET) not in out.with_suffix(".build.json").read_text()
    # on a board whose counts are open (a dry run's), the drawn cells' share stays
    open_board = Board(key="T1:val", task="T1", genes=list(t1_env.board.genes), needs_coords=False, min_cells=100,
                       max_cells=T1_SECRET)
    full = {"n_cells": T1_SECRET, "g_mom": {"edited_share": {"all": 0.5}, "edited_share_pool": {"all": 0.5}}}
    assert not open_board.counts_hidden and bt._hide_counts(full, open_board) == full
    hidden = bt._hide_counts(full, t1_env.board)
    assert "edited_share" not in hidden["g_mom"] and hidden["g_mom"]["edited_share_pool"] == {"all": 0.5}
    assert full["g_mom"]["edited_share"] == {"all": 0.5}  # the input record is left alone


def test_lib3_fit_stops_when_library_3_leaves_a_category_r3b_edits_without_factors(bt, t1_env, tmp_path,  # noqa: F811
                                                                                   monkeypatch):
    endo = {"Endocardium-1", "EndoMT Endocardium", "Early VIC"}  # every label of the endocardium category
    t1_env.ref = _libraries(t1_env.ref, lambda i, lab: 1 if lab in endo else _alternate(i, lab))
    _loads(bt, t1_env, monkeypatch)
    out = tmp_path / "T1_s1.h5ad"
    with pytest.raises(ValueError, match="no factors for a category R3b edits"):
        bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--g-mom-fit", "lib3", "--out", str(out)])
    assert not out.exists()


def test_lib3_fit_stops_when_the_names_do_not_name_libraries(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    t1_env.ref = _libraries(t1_env.ref, lambda i, lab: 1 if i % 2 else 2)  # no library 3
    _loads(bt, t1_env, monkeypatch)
    with pytest.raises(SystemExit, match="--g-mom-fit lib3: the reference's cell names do not separate library 3"):
        bt.main(["t1", "--no-feng-hygiene", "--g-mom", "--g-mom-fit", "lib3", "--out", str(tmp_path / "T1_s1.h5ad")])


def _v_cm_mostly_3(i, lab):  # cardiomyocytes mostly library 3, the rest spread over 1-3: the categories' shares differ
    if lab == "V-CM":
        return 1 if i % 4 == 0 else 3
    return i % 3 + 1


def _pool_spy(bt, monkeypatch):  # noqa: F811
    seen = []
    inner = bt.t1_scope.build_submission

    def spy(pool, board, shares, n, seed, **kw):
        seen.append({"pool": pool.copy(), "shares": dict(shares)})
        return inner(pool, board, shares, n, seed, **kw)

    monkeypatch.setattr(bt.t1_scope, "build_submission", spy)
    return seen


def _pool_shares(bt, pool):  # noqa: F811
    """Per R3b category: the share of the pool's cells of library 3 whose label R1 does not hold (the lib3 fit edits
    those), and each category's library-3 share."""
    p = pool.copy()
    t1_scope.classify(p, source="E10.5_RNA")
    cat = p.obs["category"].astype(str).to_numpy()
    lab = p.obs["celltype"].astype(str).to_numpy()
    held = np.isin(lab, sorted(ts.r1_unbracketed_labels(np.unique(lab))))
    in3 = np.array([str(n).endswith("_3") for n in p.obs_names])
    return ({c: ((cat == c) & in3 & ~held).sum() / (cat == c).sum() for c in bt.G_MOM_EDIT_CATEGORIES},
            {c: ((cat == c) & in3).sum() / (cat == c).sum() for c in bt.G_MOM_EDIT_CATEGORIES})


@pytest.mark.parametrize("base", [["--base", "mid", "--cm-share", "0.60"], ["--base", "lib3"]])
def test_on_the_two_t1_bases_the_lib3_fit_gets_library_3s_rows_of_the_loaded_object_not_the_pool(
        bt, t1_env, tmp_path, monkeypatch, base):  # noqa: F811
    """The two command lines the brief names: the fit receives exactly library 3's rows of the E10.5 object build_t1
    loaded (on mid, not the mid pool, which also holds drawn cells of libraries 1 and 2), and edited_share_pool is the
    pool's library-3, R1-unheld share per category, its "all" weighted by the shares the build draws with."""
    t1_env.ref = _libraries(t1_env.ref, _v_cm_mostly_3)
    returned = _loads(bt, t1_env, monkeypatch)
    calls = _fit_spy(bt, monkeypatch)
    seen = _pool_spy(bt, monkeypatch)
    out = tmp_path / "T1_s2.h5ad"
    bt.main(["t1", "--no-feng-hygiene", *base, "--g-mom", "--g-mom-fit", "lib3", "--out", str(out)])
    (call,), (built,) = calls, seen
    loaded = returned["E10.5_RNA.h5ad"][0]
    in3 = np.array([str(n).endswith("_3") for n in loaded.obs_names])
    assert list(call["ref"].obs_names) == list(loaded.obs_names[in3])
    assert (_counts(call["ref"]) == _counts(loaded)[in3]).all()
    pool = built["pool"]
    pool3 = np.array([str(n).endswith("_3") for n in pool.obs_names])
    assert sorted(pool.obs_names[pool3]) == sorted(loaded.obs_names[in3])  # every library-3 cell is in the pool
    rec = json.loads(out.with_suffix(".build.json").read_text())
    assert rec["base"] == base[1] and rec["g_mom"]["fit"] == "lib3"
    edited, lib3_share = _pool_shares(bt, pool)
    for c in bt.G_MOM_EDIT_CATEGORIES:
        assert rec["g_mom"]["edited_share_pool"][c] == pytest.approx(edited[c])
    w = {c: built["shares"][c] for c in bt.G_MOM_EDIT_CATEGORIES}
    assert rec["g_mom"]["edited_share_pool"]["all"] == pytest.approx(sum(w[c] * edited[c] for c in w) / sum(w.values()))
    if base[1] == "mid":
        assert not pool3.all()  # the mid pool is not library 3's rows: a fit on it would differ from the call above
        assert len(call["ref"]) < pool.n_obs
        assert all(lib3_share[c] < 1 for c in bt.G_MOM_EDIT_CATEGORIES)  # the dose cut the record shows
        assert all(rec["g_mom"]["edited_share_pool"][c] < 1 for c in bt.G_MOM_EDIT_CATEGORIES)
    else:
        assert pool3.all() and all(lib3_share[c] == 1 for c in bt.G_MOM_EDIT_CATEGORIES)


def test_edited_shares_weights_all_by_the_given_shares_and_by_the_cells_without_them(bt):  # noqa: F811
    cats = np.array(["cm"] * 8 + ["endocardium"] * 2 + ["other"] * 5)
    n = {"cm": 2, "endocardium": 2}
    drawn = bt._edited_shares(cats, n)
    assert drawn == {"all": 0.4, "cm": 0.25, "endocardium": 1.0}  # over the cells as they stand
    pool = bt._edited_shares(cats, n, weights={"cm": 0.2, "endocardium": 0.6, "other": 0.2})
    assert pool["all"] == pytest.approx((0.2 * 0.25 + 0.6 * 1.0) / 0.8) and pool["cm"] == 0.25
    assert bt._edited_shares(cats, n, weights={})["all"] is None
    assert bt._edited_shares(np.array(["cm"] * 4), {"cm": 1}, weights={"cm": 0.5, "endocardium": 0.5}) \
        == {"all": 0.25, "cm": 0.25, "endocardium": None}


def test_help_text_names_the_other_inputs_the_dose_cut_and_the_unknown_sign(bt):  # noqa: F811
    sub = next(a for a in bt.build_parser()._actions if a.choices and "t1" in a.choices)
    opt = next(a for a in sub.choices["t1"]._actions if "--g-mom-fit" in a.option_strings)
    for words in ("fewer metacells, its own PCA, OT and per-category mean", "size-matched control",
                  "a cut in the step's dose", "(1 + s) / 2", "weighted by the build's shares", "unknown in sign",
                  "With whole the record is as before the option"):
        assert words in " ".join(opt.help.split())
