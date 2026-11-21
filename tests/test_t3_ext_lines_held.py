"""Tests for the one field row K24's restated (b4) adds to experiments/t3_ext/combine.py: each bulk table's manifest
records lines_held_by_study, the number of each study's lines holding each panel gene (the share's denominator), beside
line_agreement_by_study and under the same condition. Nothing else changes: the bulk tables' csv bytes are those of
the earlier combine.py (loaded from `git show b50bf57:experiments/t3_ext/combine.py`, the last commit before the
field), and the single-cell family gains nothing.
Synthetic units only; no outside or competition data."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
T3X = REPO / "experiments" / "t3_ext"
if str(T3X) not in sys.path:
    sys.path.insert(0, str(T3X))
import combine as cb  # noqa: E402

GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES")


BASELINE = "b50bf57"   # combine.py as it was before lines_held_by_study (sha256 pinned in the loader)
BASELINE_SHA256 = "bf56610088cd7dfad8f66d4371a7a699e3e9fc05bc1911006098a22f2563c830"

def _ref(n=40, seed=3):
    rng = np.random.default_rng(seed)
    return {"genes": [f"G{i}" for i in range(n)], "groups": ["a", "b"], "shares": {"a": 0.6, "b": 0.4},
            "mu": {"a": rng.random(n) + 0.1, "b": rng.random(n) + 0.1}, "file_sha256": "r" * 64}


def _bulk_units(ref, study="MGP-DMDD", lines=13, seed=4, acc="E-MTAB-7322", tier="unrelated"):
    """Bulk units with holes: line i leaves out the panel genes i % 7 :: 7 (so each gene is held in a varying number
    of lines), genes 0 and 1 are held by no line, and one line (L0) holds two units of one gene."""
    rng = np.random.default_rng(seed)
    n = len(ref["genes"])
    common = rng.normal(0, 0.2, n)
    units = []
    for i in range(lines):
        for rep in range(2 if i == 0 else 1):
            L = common + rng.normal(0, 0.15, n)
            L[i % 7::7] = np.nan
            L[:2] = np.nan
            if rep == 1:
                L[5] = np.nan  # the second unit of L0 misses gene 5; the line still holds it through its first unit
            units.append({"accession": acc, "gene": f"L{i}", "study": study, "tier": tier, "L": {"__whole__": L},
                          "sha256": f"{i:02d}{rep}" * 21 + "u", "modality": "bulk", "k": {"E9.5": 0.64}, "m": {},
                          "libraries": {"ko": ["k1", "k2"], "control": ["c1", "c2"]},
                          "guards": {"E9.5": "passed"}})
    return units


def _sc_units(ref, seed=6):
    rng = np.random.default_rng(seed)
    n = len(ref["genes"])
    units = []
    for st, acc in (("S1", "GSE1"), ("S2", "GSE2")):
        for gene in ("Ko1", "Ko2"):
            L = {g: rng.normal(0, 0.2, n) for g in ref["groups"]}
            L["a"][:3] = np.nan
            units.append({"accession": acc, "gene": gene, "study": st, "tier": "unrelated", "L": L,
                          "sha256": f"{acc}{gene}".ljust(64, "s")})
    return units


def _main_combine(tmp_path: Path):
    """main's combine.py, loaded under another module name from `git show` (its imports resolve to the modules already
    loaded from this checkout's experiments/t3_ext, which this change leaves alone)."""
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS} | {"GIT_OPTIONAL_LOCKS": "0"}
    # the baseline is the last commit of combine.py BEFORE the field was added (b50bf57), not branch main: once
    # the change is merged, main's combine.py is the new one and would compare with itself
    r = subprocess.run(["git", "show", f"{BASELINE}:experiments/t3_ext/combine.py"], cwd=REPO, env=env,
                       capture_output=True)
    if r.returncode != 0:
        pytest.skip(f"commit {BASELINE} is not in this checkout: the byte comparison with the earlier combine.py "
                    "cannot be made")
    assert hashlib.sha256(r.stdout).hexdigest() == BASELINE_SHA256
    d = tmp_path / "main_combine"
    d.mkdir()
    p = d / "combine_main.py"
    p.write_bytes(r.stdout)
    spec = importlib.util.spec_from_file_location("combine_main", p)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    finally:
        while str(d) in sys.path:  # main's combine.py puts its own folder on sys.path; take it off again
            sys.path.remove(str(d))
    assert mod.fe is cb.fe and mod.cm is cb.cm and mod.tb is cb.tb
    return mod


def _manifest(out: Path, name: str) -> dict:
    """The manifest apart from its 'created' time, with the output folder's path written as <out>."""
    m = json.loads((out / name / f"{name}_manifest.json").read_text().replace(str(out), "<out>"))
    m.pop("created", None)
    return m


def test_the_field_is_the_count_of_finite_values_per_gene():
    ref = _ref()
    units = _bulk_units(ref)
    st = cb.bulk_study_tables(units, len(ref["genes"]))["MGP-DMDD"]
    # by hand: per line the mean over its units (a gene a line holds in any unit is held), then count the lines
    per_line = []
    for g in sorted({u["gene"] for u in units}):
        rows = np.vstack([u["L"]["__whole__"] for u in units if u["gene"] == g])
        per_line.append(np.isfinite(rows).any(0))
    want = np.vstack(per_line).sum(0)
    assert st["lines_held"] == [int(v) for v in want]
    assert all(type(v) is int for v in st["lines_held"]) and len(st["lines_held"]) == len(ref["genes"])
    # 13 lines, each leaving out one gene in seven: every other gene is held in 11 or 12 lines
    assert st["lines_held"][0] == st["lines_held"][1] == 0 and set(st["lines_held"][2:]) == {11, 12}
    assert len(set(st["lines_held"][2:])) > 1  # the holes give genes different counts
    # it is the share's denominator: a gene held in no line has no share, every other gene a share in k / n_g
    for n_g, sh in zip(st["lines_held"], st["line_agreement"], strict=True):
        assert (sh is None) == (n_g == 0)
        if sh is not None:
            assert abs(sh * n_g - round(sh * n_g)) < 1e-9


def test_the_manifests_record_it_beside_line_agreement_under_the_same_condition(tmp_path):
    ref = _ref()
    units = _bulk_units(ref) + _bulk_units(ref, study="OTHER", lines=1, seed=9, acc="E-MTAB-0000")
    out = tmp_path / "t"
    cb.bulk_tables(units, ref, {"left_out": [], "not_cleared": [], "inputs": {}}, out)
    studies = cb.bulk_study_tables(units, len(ref["genes"]))
    for name in ("bulk_all", "bulk_study_MGP-DMDD", "bulk_study_OTHER", "bulk_all__loso_OTHER",
                 "bulk_all__loso_MGP-DMDD"):
        m = _manifest(out, name)
        assert set(m["lines_held_by_study"]) == set(m["line_agreement_by_study"]), name
        for s, v in m["lines_held_by_study"].items():
            assert v == studies[s]["lines_held"] and all(type(x) is int for x in v)
    # a one-line study gets neither record
    assert _manifest(out, "bulk_study_OTHER")["lines_held_by_study"] == {}
    assert set(_manifest(out, "bulk_all")["lines_held_by_study"]) == {"MGP-DMDD"}


def test_the_bulk_csv_bytes_are_main_s_and_only_the_field_is_added(tmp_path):
    main_cb = _main_combine(tmp_path)
    ref = _ref()
    units = _bulk_units(ref) + _bulk_units(ref, study="S3", lines=4, seed=11, acc="E-MTAB-1111") \
        + _bulk_units(ref, study="S4", lines=3, seed=12, acc="E-MTAB-2222")
    rec = {"left_out": [], "not_cleared": [], "inputs": {}}
    a, b = tmp_path / "mine", tmp_path / "main"
    wa = cb.bulk_tables(units, ref, rec, a)
    wb = main_cb.bulk_tables(units, ref, rec, b)
    assert set(wa) == set(wb) and len(wa) == 1 + 3 + 3
    for name in wa:
        csvs = sorted(p.relative_to(a) for p in (a / name).rglob("*.csv"))
        assert csvs and csvs == sorted(p.relative_to(b) for p in (b / name).rglob("*.csv"))
        for p in csvs:
            assert (a / p).read_bytes() == (b / p).read_bytes(), p
        ma, mb = _manifest(a, name), _manifest(b, name)
        assert "lines_held_by_study" in ma and "lines_held_by_study" not in mb
        ma.pop("lines_held_by_study")
        assert ma == mb, name
    # the sign rule's bulk_study_tables values are main's too
    sa, sb = cb.bulk_study_tables(units, len(ref["genes"])), main_cb.bulk_study_tables(units, len(ref["genes"]))
    for s in sa:
        assert np.array_equal(sa[s]["L"]["__whole__"], sb[s]["L"]["__whole__"], equal_nan=True)
        assert sa[s]["line_agreement"] == sb[s]["line_agreement"] and "lines_held" not in sb[s]


def test_the_single_cell_tables_gain_nothing(tmp_path):
    main_cb = _main_combine(tmp_path)
    ref = _ref()
    units = _sc_units(ref)
    rec = {"left_out": [], "not_cleared": [], "inputs": {}}
    a, b = tmp_path / "mine", tmp_path / "main"
    wa = cb.single_cell_tables(units, ref, rec, a)
    wb = main_cb.single_cell_tables(units, ref, rec, b)
    assert set(wa) == set(wb) and "all" in wa
    for name in wa:
        ma, mb = _manifest(a, name), _manifest(b, name)
        assert ma == mb and "lines_held_by_study" not in ma, name
        for p in sorted((a / name).rglob("*.csv")):
            assert p.read_bytes() == (b / p.relative_to(a)).read_bytes(), p
