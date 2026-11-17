"""build_test.py t1 --direction-fraction F (campaign branch campaign/t1-n; no ruling): F replaces the bracket's share of
days in M1's direction step and nothing else changes. Synthetic only: the T1 builder runs end to end on the stages of
tests/test_hidden_counts.py, every input monkeypatched."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from test_hidden_counts import T1_GENES, T1_HELD, T1_SECRET, UNION, bt, t1_env  # noqa: F401

BRACKET = 0.4  # feng_dirs.bracket_fraction(10.5, 12.5)
NEW_KEYS = ("bracket_fraction", "fraction_source")


def _loads(bt, env, monkeypatch):  # noqa: F811
    """E10.5 and E9.5 (the --g-mom fit), the union of R2 (on, the default form), and a spy on build_submission that
    keeps every call's arguments and the labels of the population g_mom's hook receives."""
    monkeypatch.setattr(bt, "load_stage",
                        lambda path, **kw: (env.ref if Path(path).name.startswith("E10.5") else env.start).copy())
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: UNION)
    calls = []

    def spy(*a, **k):
        call = {"args": a, "kw": dict(k)}
        hook = k.get("post_cells")
        if hook is not None:
            def wrapped(pop):
                call["labels"] = pop.obs["celltype"].astype(str).to_numpy()
                return hook(pop)
            k = {**k, "post_cells": wrapped}
        calls.append(call)
        return env.build(*a, **k)

    monkeypatch.setattr(bt.t1_scope, "build_submission", spy)
    return calls


def _same(x, y) -> bool:
    if isinstance(x, ad.AnnData):
        return (list(x.obs_names) == list(y.obs_names) and list(x.var_names) == list(y.var_names)
                and (x.X != y.X).nnz == 0 and x.obs.equals(y.obs))
    if isinstance(x, dict):
        return type(x) is type(y) and list(x) == list(y) and all(_same(x[k], y[k]) for k in x)
    if isinstance(x, np.ndarray):
        return x.dtype == y.dtype and np.array_equal(x, y, equal_nan=True)
    if callable(x):
        return callable(y) and x.__qualname__ == y.__qualname__
    return type(x) is type(y) and x == y


def test_the_option_is_absent_by_default_and_parses_a_float(bt):  # noqa: F811
    ap = bt.build_parser()
    assert ap.parse_args(["t1", "--out", "x.h5ad"]).direction_fraction is None
    assert ap.parse_args(["t1", "--direction-fraction", "0.8", "--out", "x.h5ad"]).direction_fraction == 0.8
    assert bt.feng_dirs.bracket_fraction(10.5, 12.5) == BRACKET


@pytest.mark.parametrize("base", [[], ["--base", "lib3"], ["--base", "mid", "--cm-share", "0.60"]])
def test_default_path_passes_the_bracket_fraction_and_the_option_changes_only_it(bt, t1_env, tmp_path,  # noqa: F811
                                                                                 monkeypatch, capsys, base):
    """Without the option build_submission receives the bracket's share as a float; with it, F, and every other
    argument equals the default call's (same types, same values; g_mom's hook is a new closure per build, the same
    function). The record: fraction = the applied share, bracket_fraction 0.4, fraction_source; nothing else differs.
    Only the cells Feng's direction covers change; the count, the cells drawn and their labels do not."""
    ref = t1_env.ref
    lab = ref.obs["celltype"].astype(str).to_numpy()
    lib = np.where(lab == "V-CM", np.where(np.arange(ref.n_obs) % 4 == 0, 1, 3), np.arange(ref.n_obs) % 3 + 1)
    ref.obs_names = [f"bc{i}-1_{k}" for i, k in enumerate(lib)]
    calls = _loads(bt, t1_env, monkeypatch)
    argv = ["t1", "--g-mom", *base]
    bt.main([*argv, "--out", str(tmp_path / "default.h5ad")])
    bt.main([*argv, "--direction-fraction", "0.8", "--out", str(tmp_path / "option.h5ad")])
    d, o = calls
    assert type(d["kw"]["direction_fraction"]) is float and d["kw"]["direction_fraction"] == BRACKET
    assert type(o["kw"]["direction_fraction"]) is float and o["kw"]["direction_fraction"] == 0.8
    assert len(d["args"]) == len(o["args"]) == 5 and all(_same(x, y) for x, y in zip(d["args"], o["args"], strict=True))
    assert set(d["kw"]) == set(o["kw"])
    for k in d["kw"]:
        if k != "direction_fraction":
            assert _same(d["kw"][k], o["kw"][k]), k
    assert callable(d["kw"]["post_cells"])  # g_mom's hook on both
    # the record
    rd, ro = (json.loads((tmp_path / f"{k}.build.json").read_text()) for k in ("default", "option"))
    assert (rd["fraction"], rd["bracket_fraction"], rd["fraction_source"]) == (BRACKET, BRACKET, "bracket")
    assert (ro["fraction"], ro["bracket_fraction"], ro["fraction_source"]) == (0.8, BRACKET, "option")
    assert list(rd) == list(ro)
    differ = {k for k in rd if rd[k] != ro[k]}
    assert differ <= {"fraction", "fraction_source", "g_mom"}
    if "g_mom" in differ:  # the hook's median depth scale is read off the cells it edits, which the step moved
        assert {k for k in rd["g_mom"] if rd["g_mom"][k] != ro["g_mom"][k]} == {"median_depth_scale"}
    assert str(T1_SECRET) not in capsys.readouterr().out
    # the files: the same cells in the same order; only the cells Feng's direction covers differ
    labels = d["labels"]
    assert np.array_equal(labels, o["labels"])
    groups = d["kw"]["direction_groups"]
    covered = np.array([groups.get(x) in d["kw"]["directions"] for x in labels])
    fd, fo = ad.read_h5ad(tmp_path / "default.h5ad"), ad.read_h5ad(tmp_path / "option.h5ad")
    xd, xo = np.asarray(fd.X.todense() if hasattr(fd.X, "todense") else fd.X), \
        np.asarray(fo.X.todense() if hasattr(fo.X, "todense") else fo.X)
    assert xd.shape == xo.shape and covered.any() and (~covered).any()
    assert np.array_equal(xd[~covered], xo[~covered])
    assert not np.array_equal(xd[covered], xo[covered])


def test_the_option_at_the_bracket_share_rebuilds_the_default_file(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    """F = 0.4 given explicitly writes the default file byte for byte: the option touches nothing but the share."""
    _loads(bt, t1_env, monkeypatch)
    bt.main(["t1", "--g-mom", "--out", str(tmp_path / "default.h5ad")])
    bt.main(["t1", "--g-mom", "--direction-fraction", "0.4", "--out", str(tmp_path / "option.h5ad")])
    assert (tmp_path / "default.h5ad").read_bytes() == (tmp_path / "option.h5ad").read_bytes()
    ro = json.loads((tmp_path / "option.build.json").read_text())
    assert ro["fraction_source"] == "option" and ro["fraction"] == ro["bracket_fraction"] == BRACKET


def test_the_default_record_is_todays_plus_the_two_keys(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    """Without the option the record holds today's keys in today's order, with the two new keys after "fraction"."""
    _loads(bt, t1_env, monkeypatch)
    bt.main(["t1", "--g-mom", "--out", str(tmp_path / "d.h5ad")])
    rec = json.loads((tmp_path / "d.build.json").read_text())
    keys = [k for k in rec if k not in NEW_KEYS]
    assert keys == ["entry", "board", "dry_run", "n_cells", "seed", "fraction", "herd", "hold_unbracketed",
                    "feng_hygiene", "base", "cm_share_of_heart", "blood_share_of_chamber", "shares", "external",
                    "check_ok", "g_mom"]
    assert list(rec).index("bracket_fraction") == list(rec).index("fraction") + 1


@pytest.mark.parametrize("value", ["nan", "inf", "-inf", "0", "0.0", "-0.2", "1.0000001", "2"])
def test_the_option_refuses_a_share_outside_zero_one(bt, t1_env, tmp_path, monkeypatch, value):  # noqa: F811
    monkeypatch.setattr(bt, "load_stage", lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    with pytest.raises(SystemExit, match="must be finite and in \\(0, 1\\]"):
        bt.main(["t1", "--g-mom", f"--direction-fraction={value}", "--out", str(tmp_path / "T1.h5ad")])  # = : "-0.2"
    assert not (tmp_path / "T1.h5ad").exists()


def test_the_option_accepts_the_whole_bracket(bt, t1_env, tmp_path, monkeypatch):  # noqa: F811
    calls = _loads(bt, t1_env, monkeypatch)
    bt.main(["t1", "--g-mom", "--direction-fraction", "1", "--out", str(tmp_path / "T1.h5ad")])
    assert calls[0]["kw"]["direction_fraction"] == 1.0


@pytest.mark.parametrize("extra, msg", [
    ([], "refused without --g-mom"),
    (["--g-mom", "--strong-null"], "refused with --strong-null"),
    (["--g-mom", "--herd"], "refused with --herd"),
    (["--g-mom", "--no-hold-unbracketed", "--no-feng-hygiene"], "refused with --no-hold-unbracketed"),
    (["--g-mom", "--base", "lib3", "--cm-share", "0.6", "--lib3-pair"], "refused with --lib3-pair"),
])
def test_the_option_is_refused_outside_the_form_the_reads_read(bt, t1_env, tmp_path, monkeypatch, extra,  # noqa: F811
                                                                msg):
    monkeypatch.setattr(bt, "load_stage", lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    monkeypatch.setattr(bt.hygiene, "load_union", lambda *a, **k: pytest.fail("nothing is read before the refusal"))
    with pytest.raises(SystemExit, match=msg):
        bt.main(["t1", *extra, "--direction-fraction", "0.8", "--out", str(tmp_path / "T1.h5ad")])
    assert not (tmp_path / "T1.h5ad").exists()


def test_the_freeze_builds_do_not_carry_the_option():
    text = (Path(__file__).resolve().parents[1] / "experiments" / "builds" / "build_all.sh").read_text()
    assert "--direction-fraction" not in text
    assert T1_HELD and T1_GENES  # the shared synthetic inputs are the ones the other builder tests use
