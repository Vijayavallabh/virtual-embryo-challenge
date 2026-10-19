"""Tests for the T3 submission builder (experiments/t3_decomp/build_submission.py).

Offline, synthetic and fast: `uv run pytest experiments/t3_decomp/test_build_submission.py -q`. No competition data: a
small synthetic WT file stands in for WT E8.75 and synthetic factor files for a run's factors/. Decision files are made
by t3_decompose.decide_submissions itself, so every table row a test selects is the row the release-day run would
write. The end-to-end test through t3_decompose.build_context (`factors`, then `build` with stage mixing) is skipped
unless RUN_T3_E2E=1 (run it on the DGX; it takes a minute or two).
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import os
import shutil
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("build_submission", HERE / "build_submission.py")
bs = importlib.util.module_from_spec(_spec)
sys.modules["build_submission"] = bs
_spec.loader.exec_module(bs)
t3d = bs.t3d
RULES = yaml.safe_load((HERE / "rules.yaml").read_text())
SCRIPT_SHA = t3d.sha256(HERE / "t3_decompose.py")
RULES_SHA = t3d.sha256(HERE / "rules.yaml")
HASH_OK = {"mismatch": []}
G = 40
GENES = [f"G{j}" for j in range(G)]
# (celltype, cm_celltype, share): in the combined label mode the vocabulary holds every depleted scenario state
# (OFT/RV-CM, aSHF, pSHF, aPHM, pPHM), a spared one (vCM1), three strict non-Mesp1 labels and two other labels
CELLS = [("V-CM", "OFT/RV-CM", 0.10), ("V-CM", "vCM1", 0.14), ("aSHF", "Unknown", 0.08), ("pSHF", "Unknown", 0.06),
         ("aPHM", "Unknown", 0.08), ("pPHM", "Unknown", 0.08), ("NCC", "Unknown", 0.12), ("SE", "Unknown", 0.12),
         ("Foregut", "Unknown", 0.08), ("Peri", "Unknown", 0.07), ("Endo", "Unknown", 0.07)]
OUTCOMES = ("D1", "D2", "D2u", "D4", "D5", "D6", "D7", "D7m", "D8", "D9", "D10", "D11a", "D11b")
M, M_RAW, M7M, W1 = 0.8, 0.5, 0.7, (0.75, 1.0)      # D7 m, D7 m_raw, D7m m, W1 (chosen, runner-up)
PARAMS = {"depth_ratio_pooled": {"vCM1": 1.5, "NCC": 0.7},
          "stage_alphas_mean": {"E9.5": 0.2, "E8.25_late": 0.1},
          "coverage_shares_mean": {"NCC": 0.2, "SE": 0.05, "Foregut": 0.1}}


def _phi(up=(), down=(), up_by=2.0):
    phi = np.ones(G)
    phi[list(up)] = up_by
    phi[list(down)] = 0.5
    return phi


PHIS = {"gata4_pooled": _phi(range(0, 10), range(10, 20)), "reference_R": _phi(range(10, 20), range(0, 10)),
        "transfer_mab21l2_all": _phi(range(20, 30), range(30, 35)), "mab21l2_pooled": _phi(range(0, 5), up_by=3.0)}


def _write_wt(path: Path, rng, n: int = 900) -> None:
    import anndata as ad
    from scipy import sparse

    from vec.prep import log1p_cp10k

    prof = rng.lognormal(0, 1.0, size=(len(CELLS), G))
    w = np.array([c[2] for c in CELLS])
    lab = rng.choice(len(CELLS), n, p=w / w.sum())
    C = np.vstack([rng.poisson(prof[i] / prof[i].sum() * rng.lognormal(np.log(150), 0.3)) for i in lab])
    C = C.astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [CELLS[i][0] for i in lab], "cm_celltype": [CELLS[i][1] for i in lab]},
                       index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 400, n), rng.uniform(0, 400, n), rng.integers(0, 12, n) * 10.0]
    a.write_h5ad(path)


def _write_factors(run: Path, names, preset: str, wt_sha: str) -> dict:
    fdir = run / "factors"
    fdir.mkdir(parents=True)
    files = {}
    for k in names:
        p = fdir / f"{k}.csv"
        pd.DataFrame({"gene": GENES, "factor": PHIS[k]}).to_csv(p, index=False)
        files[k] = {"file": str(p), "sha256": t3d.sha256(p)}
    rec = {"preset": preset, "files": files, "parameters": PARAMS, "inputs": {"wt": {"path": "wt", "sha256": wt_sha}}}
    (fdir / "factors.json").write_text(json.dumps(rec))
    return files


def _with_stages(ref, rng):
    """The reference plus two synthetic stage pools in the gata4 preset's order (every label has >= 30 cells)."""
    stages = []
    for name in RULES["presets"]["gata4"]["stage_pools"]:
        lab = np.repeat(np.arange(ref.V.L), 60)
        C = rng.poisson(5, size=(len(lab), G)).astype(np.int32)
        C[C.sum(1) == 0, 0] = 1
        pool = t3d.Pool(C, lab, (C / C.sum(1, keepdims=True)).astype(np.float32), ref.V.L)
        stages.append(t3d.StagePool(name, pool, np.zeros((ref.V.L, G)), pool.count, {}))
    return replace(ref, stages=stages)


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    from vec.boards import Board

    d = tmp_path_factory.mktemp("t3build")
    rng = np.random.default_rng(0)
    wt = d / "wt.h5ad"
    _write_wt(wt, rng)
    wt_sha = t3d.sha256(wt)
    run = d / "run"
    files = _write_factors(run, ["gata4_pooled", "reference_R", "transfer_mab21l2_all"], "gata4", wt_sha)
    run_smoke = d / "run_smoke"
    _write_factors(run_smoke, ["mab21l2_pooled", "reference_R", "transfer_mab21l2_all"], "mab21l2", wt_sha)
    ref = bs.load_reference(wt, "gata4", RULES, threads=2)
    board = Board(key="T3:synthetic", task="T3", genes=GENES, needs_coords=True, min_cells=200, max_cells=3000)
    boards_dir = d / "boards"
    boards_dir.mkdir()
    (boards_dir / "T3__synthetic.genes.txt").write_text("\n".join(GENES) + "\n")
    (boards_dir / "index.json").write_text(json.dumps({"T3:synthetic": {
        "task": "T3", "genes_file": "T3__synthetic.genes.txt", "needs_coords": True, "min_cells": 200,
        "max_cells": 3000}}))
    return SimpleNamespace(d=d, wt=wt, run=run, run_smoke=run_smoke, params={**PARAMS, "factor_files": files},
                           ref=ref, ref_stages=_with_stages(ref, rng), board=board, boards_dir=boards_dir)


def rule_block(outcomes=None, m=M, m_raw=M_RAW, m7m=M7M, w1=W1) -> dict:
    o = dict.fromkeys(OUTCOMES, False)
    o.update({k: v for k, v in (outcomes or {}).items() if k in OUTCOMES})
    r = {k: {"outcome": v} for k, v in o.items()}
    r["D7"].update(m=m if o["D7"] else 0.0, m_raw=m_raw, n_counted=int(bool(o["D7"])))
    r["D7m"].update(m=m7m if o["D7m"] else 0.0, n_counted=int(bool(o["D7m"])))
    r["m_raw_positive"] = {"outcome": bool((outcomes or {}).get("m_raw_positive", False))}
    r["D3"] = {"outcome": {"factor": True}}              # a dict outcome, as in real files: not a table name
    r["W1"] = {"grid": [], "chosen": w1[0], "runner_up": w1[1], "applicable": True}
    return r


def make_decision(world, outcomes=None, preset="gata4", table=True, rules_kw=None, **over) -> dict:
    rules = rule_block(outcomes, **(rules_kw or {}))
    dec = {"preset": preset, "rules_sha256": RULES_SHA, "script_sha256": SCRIPT_SHA, "protocol_as_registered": True,
           "seeds": list(range(10)), "rules": rules}
    if table:
        dec["submissions"] = bs.roundtrip(t3d.decide_submissions(RULES, rules, None,
                                                                 {"submission_parameters": world.params}))
    else:
        dec["submissions"] = "not applicable (smoke test on the training knockout; the gata4 run decides)"
    dec.update(over)
    return dec


def write_decision(path: Path, dec: dict) -> Path:
    path.write_text(json.dumps(dec))
    return path


def table_names(dec: dict) -> dict:
    return {k: v["outcome"] for k, v in dec["rules"].items()
            if isinstance(v, dict) and "outcome" in v and not isinstance(v["outcome"], dict)}


EXTRA_KEYS = ("depth", "stage", "coverage")              # rules.yaml submissions.allowed.gata4_extras
GATA4 = ("gata4", "gata4_unvalidated")                   # the factor sources that are the pooled Gata4 factor


def comp_from_text(kind: str, dec: dict) -> dict:
    """A composition written from the rules text: the scenario's depleted states at 2^-(m x weight), m from D7."""
    if kind == "none":
        return {"kind": "none"}
    d7 = dec["rules"]["D7"]
    l2 = (d7["m_raw"] if kind == "mild_raw" else d7["m"]) * {"central": 1.0, "mild": 0.5, "mild_raw": 0.5}[kind]
    scen = RULES["presets"]["gata4"]["scenario"]
    return {"kind": kind, "depleted_states": scen["depleted_states"], "spared_states": scen["spared_states"],
            "log2_magnitude": l2, "depletion_factor": 2.0 ** -l2, "magnitude_source": "Gata4 replicates (D7)"}


def allowed_from_text(dec: dict, s: str) -> dict:
    """The allowed branch's s1 or s2 written out from the text of rules.yaml notes 6, 7 and 8, independently of the
    builder: the stored table, replaced by note 8 when no factor check passed, then note 6 as note 7 reads it."""
    o = table_names(dec)
    off = dict.fromkeys(EXTRA_KEYS, False)
    if not (o["D2"] or o["D2u"] or o["D10"]):                     # note 8: the factor table's last row
        if s == "s1":
            f = {"source": "mab21l2", "exponent": 1.0} if o["D11b"] is True else {"source": "none", "exponent": 0.0}
            return {"factor": f, "composition": comp_from_text("central" if o["D7"] else "none", dec), **off}
        kind = ("mild" if o["D8"] else "central") if o["D7"] else "none"         # composition_s2_other_factor
        return {"factor": {"source": "gata4_unvalidated", "exponent": 0.5}, "composition": comp_from_text(kind, dec),
                **off}
    tab = dec["submissions"]["allowed"]
    s1, both = tab["s1"], {tab["s1"]["factor"]["source"], tab["s2"]["factor"]["source"]} <= set(GATA4)
    if s == "s1" or o["D11b"] is True or not both:
        return tab[s]
    if o["D7"]:                                                   # note 6 with note 7's D8 hedge
        return {"factor": {"source": "none", "exponent": 0.0},
                "composition": comp_from_text("mild" if o["D8"] is True else "central", dec), **off}
    # note 6 without D7: the Gata4 factor at 0.5 and s1's composition, no extras; at s1's own dose, half of it (note 7)
    return {"factor": {"source": s1["factor"]["source"], "exponent": 0.25 if s1["factor"]["exponent"] == 0.5 else 0.5},
            "composition": s1["composition"], **off}


def check_honoured(world, ref, dec: dict, branch: str, spec: dict, kw: dict) -> None:
    """kw (what t3_decompose.build receives) implements the spec: the named factor file and exponent, the scenario
    composition at the decision file's magnitude, and the extras with the run's parameters."""
    src = spec["factor"]["source"]
    if src == "none":
        assert kw["phi"] is None
    else:
        assert np.array_equal(kw["phi"], PHIS[bs.factor_name(src, dec["preset"])])
        assert kw["w"] == spec["factor"]["exponent"]
    comp = spec["composition"]
    f = ref.V.p.copy()
    if comp["kind"] != "none":
        m = dec["rules"]["D7m" if branch != "allowed" else "D7"]["m"]
        base = dec["rules"]["D7"]["m_raw"] if comp["kind"] == "mild_raw" else m
        assert np.isclose(comp["log2_magnitude"], base * RULES["weights"]["composition"][comp["kind"]])
        f[ref.V.index(comp["depleted_states"])] *= 2.0 ** -comp["log2_magnitude"]
        f /= f.sum()
    if spec.get("coverage"):
        fc = f.copy()
        for lab, s in PARAMS["coverage_shares_mean"].items():
            fc[ref.V.index([lab])] = s
        nm = ref.V.nm
        f = np.where(nm, fc, f * (1.0 - fc[nm].sum()) / f[~nm].sum())
        f /= f.sum()
    assert np.allclose(kw["f"], f, rtol=0, atol=1e-12)
    if spec.get("depth"):
        expect = np.ones(ref.V.L)
        for lab, r in PARAMS["depth_ratio_pooled"].items():
            expect[ref.V.index([lab])] = r
        assert np.array_equal(kw["ratio"], expect)
    else:
        assert kw["ratio"] is None
    if spec.get("stage"):
        assert np.array_equal(kw["alphas"], [0.2, 0.1])
    else:
        assert kw["alphas"] is None


# ------------------------------------------------------------------------------------------------- decision tables

def test_every_table_row_is_honoured(world):
    """Every row of every allowed/refused table is reached by some combination of rule outcomes (with D8 also not
    evaluated, a Gata4 extra on and W1 at 0.5), and the builder builds what the rules text gives: the stored table on
    the refused branch, the table with notes 6, 7 and 8 on the allowed branch (allowed_from_text)."""
    sub = RULES["submissions"]
    keys = ("D2", "D2u", "D10", "D11b", "D7", "D7m", "m_raw_positive", "D4")
    hit, note_cases = set(), set()
    for combo, d8, w1 in itertools.product(itertools.product([True, False], repeat=len(keys)), (True, False, None),
                                           (W1, (0.5, 1.0))):
        dec = make_decision(world, {**dict(zip(keys, combo, strict=True)), "D8": d8}, rules_kw={"w1": w1})
        names = table_names(dec)
        for branch in ("allowed", "refused"):
            b = sub[branch]
            row = t3d.first_match(b["factor_table"], names)
            hit.add((branch, "factor_table", b["factor_table"].index(row)))
            c1 = t3d.first_match(b["composition_s1"], names)
            hit.add((branch, "composition_s1", b["composition_s1"].index(c1)))
            tab = "composition_s2_other_factor" if row["s1"] != row["s2"] else "composition_s2_same_factor"
            c2 = t3d.first_match(b[tab], names)
            hit.add((branch, tab, b[tab].index(c2)))
            specs = {}
            for s, c in (("s1", c1), ("s2", c2)):
                spec, info = bs.resolve_spec(dec, RULES, branch, s, HASH_OK)
                specs[s] = spec
                assert info["recomputed_matches"] is True
                assert ("note6" in info) == ("note8" in info) == (branch == "allowed")
                if branch == "allowed":
                    assert spec == allowed_from_text(dec, s)
                    if s == "s2":
                        note_cases.add((info["note8"]["applies"], info["note6"]["applies"], names["D11b"] is True,
                                        bool(names["D7"])))
                else:
                    assert spec == dec["submissions"][branch][s]
                    assert spec["factor"]["source"] == row[s] and spec["composition"]["kind"] == c["value"]
                kw, _ = bs.inputs_for(spec, dec, world.ref_stages, world.run)
                check_honoured(world, world.ref_stages, dec, branch, spec, kw)
            if branch == "allowed" and names["D11b"] is not True and specs["s1"]["factor"]["source"] in GATA4:
                assert specs["s2"]["factor"] != specs["s1"]["factor"]     # never s1's Gata4 factor at s1's dose
    every = {(br, k, i) for br in ("allowed", "refused") for k, tab in sub[br].items() if isinstance(tab, list)
             for i in range(len(tab))}
    assert every <= hit, sorted(every - hit)
    # note 8 and note 6 each apply somewhere, never together; note 6 with and without D7, never when D11b passes
    assert not any(n8 and n6 for n8, n6, _, _ in note_cases) and not any(n6 and d for _, n6, d, _ in note_cases)
    assert {(True, False, False, True), (True, False, False, False), (True, False, True, True)} <= note_cases
    assert {(False, True, False, True), (False, True, False, False), (False, False, True, True)} <= note_cases


def test_rules_not_evaluated_count_as_failed_in_the_notes(world):
    """Every combination of True, False and None (not evaluated) for the rules the notes read: the built specs are the
    ones the notes' text gives (an outcome that could not be evaluated counts as false), and a D7 that was not
    evaluated carries no magnitude, as in real files."""
    keys = ("D2", "D2u", "D10", "D11b", "D7", "D8")
    for combo in itertools.product((True, False, None), repeat=len(keys)):
        dec = make_decision(world, dict(zip(keys, combo, strict=True)))
        if dec["rules"]["D7"]["outcome"] is None:
            dec["rules"]["D7"] = {"outcome": None}
            dec["submissions"] = bs.roundtrip(t3d.decide_submissions(RULES, dec["rules"], None,
                                                                     {"submission_parameters": world.params}))
        for s in ("s1", "s2"):
            assert bs.resolve_spec(dec, RULES, "allowed", s, HASH_OK)[0] == allowed_from_text(dec, s), (combo, s)
    # the three readings a mutant could change without failing the True/False sweep, written out
    dec = make_decision(world, {"D11b": None})                  # note 8, D11b not evaluated: WT cells, not Mab21l2
    assert bs.resolve_spec(dec, RULES, "allowed", "s1", HASH_OK)[0]["factor"] == {"source": "none", "exponent": 0.0}
    dec = make_decision(world, {"D2": None, "D2u": None, "D10": None})   # none evaluated: note 8 applies
    assert bs.resolve_spec(dec, RULES, "allowed", "s1", HASH_OK)[1]["note8"]["applies"]
    dec = make_decision(world, {"D2": True, "D7": None})        # note 6 without D7: the Gata4 factor at 0.5
    s1, s2 = (bs.resolve_spec(dec, RULES, "allowed", s, HASH_OK)[0] for s in ("s1", "s2"))
    assert s2["factor"] == {"source": "gata4", "exponent": 0.5} and s2["composition"] == s1["composition"]


def test_records_describe_what_is_built(world):
    """Once a note changes the pair, the table's note and source move aside (table_note, table_source) and the
    replaced submission names the note."""
    dec = make_decision(world, {"D2": True})                    # note 6 replaces s2; the table's s2 equalled its s1
    (s1, i1), (s2, i2) = (bs.resolve_spec(dec, RULES, "allowed", s, HASH_OK) for s in ("s1", "s2"))
    table_note = dec["submissions"]["allowed"]["note"]
    assert table_note and i1["table_note"] == i2["table_note"] == table_note
    assert i1["note"].startswith("rules.yaml note 6") and i2["note"].startswith("rules.yaml note 6")
    assert i1["source"] == "decision file submissions.allowed.s1" and i1["factor_row"] is not None   # s1 as stored
    assert i2["source"] == "rules.yaml note 6 (replacing decision file submissions.allowed.s2)"
    assert i2["factor_row"] is None and i2["table_factor_row"] == dec["submissions"]["allowed"]["factor_row"]
    dec = make_decision(world, {})                              # note 8 replaces both
    for s in ("s1", "s2"):
        info = bs.resolve_spec(dec, RULES, "allowed", s, HASH_OK)[1]
        assert info["source"] == f"rules.yaml note 8 (replacing decision file submissions.allowed.{s})"
        assert info["table_factor_row"]["s1"] == "gata4_unvalidated" and info["factor_row"] is None
    info = bs.resolve_spec(make_decision(world, {"D2": True, "D10": True}), RULES, "allowed", "s2", HASH_OK)[1]
    assert "table_note" not in info and info["source"] == "decision file submissions.allowed.s2"   # no note applied


def test_note6_gates_the_second_submission(world):
    """rules.yaml note 6: unless D11b passes, the allowed s1 and s2 do not both rely on the Gata4 factor."""
    per_target = {"rep1": {"t3": 57.0, "sev_floor_frac": 0.0}, "rep2": {"t3": 58.5, "sev_floor_frac": 0.1}}

    def resolve(outcomes, **rules_kw):
        dec = make_decision(world, outcomes, rules_kw=rules_kw)
        dec["rules"]["D11b"]["per_target"] = per_target
        return dec, [bs.resolve_spec(dec, RULES, "allowed", s, HASH_OK) for s in ("s1", "s2")]

    dec, [(s1, _), (s2, info)] = resolve({"D2": True, "D11b": True})
    assert s2 == dec["submissions"]["allowed"]["s2"] and not info["note6"]["applies"]
    assert s1["factor"]["source"] == s2["factor"]["source"] == "gata4"         # D11b passed: both may use it
    for d11b in (False, None):                                                  # failed or not evaluated
        dec, [(s1, i1), (s2, info)] = resolve({"D2": True, "D7": True, "D11b": d11b})
        assert s1 == dec["submissions"]["allowed"]["s1"] and s1["composition"]["kind"] == "central"
        assert s2["factor"] == {"source": "none", "exponent": 0.0} and s2["composition"]["kind"] == "central"
        assert np.isclose(s2["composition"]["log2_magnitude"], M) and info["note6"]["applies"]
        assert info["note6"]["D11b_per_target"] == per_target and i1["note6"]["D11b_outcome"] is d11b
    dec, [(s1, _), (s2, info)] = resolve({"D2": True, "D7": True, "D8": True, "D11b": False})   # note 7: D8 hedge kept
    assert s2["factor"] == {"source": "none", "exponent": 0.0} and s2["composition"]["kind"] == "mild"
    assert np.isclose(s2["composition"]["log2_magnitude"], M / 2) and "D8" in info["note6"]["reason"]
    dec, [(s1, _), (s2, info)] = resolve({"D2": True})
    assert s2 == {"factor": {"source": "gata4", "exponent": 0.5}, "composition": s1["composition"],
                  **dict.fromkeys(bs.EXTRAS, False)}
    table_note = dec["submissions"]["allowed"]["note"]                          # the table's s2 equalled its s1
    assert table_note and info["table_note"] == table_note and info["note"].startswith("rules.yaml note 6")
    dec, [(s1, _), (s2, info)] = resolve({"D2": True}, w1=(0.5, 1.0))           # s2 would equal s1
    assert s1["factor"]["exponent"] == 0.5 and s2["factor"]["exponent"] == 0.25 and "half" in info["note6"]["reason"]
    dec, [(s1, _), (s2, info)] = resolve({"D2": True, "D4": True}, w1=(0.5, 1.0))   # s1 has depth: same rule (note 7)
    assert s1["factor"] == {"source": "gata4", "exponent": 0.5} and s1["depth"] is True
    assert s2 == {"factor": {"source": "gata4", "exponent": 0.25}, "composition": {"kind": "none"}, "depth": False,
                  "stage": False, "coverage": False}
    dec, [(s1, _), (s2, info)] = resolve({"D2": True, "m_raw_positive": True})   # the table's s2 had mild_raw
    assert dec["submissions"]["allowed"]["s2"]["composition"]["kind"] == "mild_raw"
    assert s2["composition"] == {"kind": "none"} and info["note6"]["applies"]   # note 7: s1's composition
    dec, [(s1, _), (s2, info)] = resolve({})                  # the last row: note 8 replaces both, note 6 is idle
    assert s2["factor"] == {"source": "gata4_unvalidated", "exponent": 0.5}
    assert info["note8"]["applies"] and not info["note6"]["applies"]
    dec, [(s1, _), (s2, info)] = resolve({"D2": True, "D10": True})            # s2 is R: the table stands
    assert s2 == dec["submissions"]["allowed"]["s2"] and s2["factor"]["source"] == "reference"
    assert not info["note6"]["applies"]
    kw, _ = bs.inputs_for(resolve({"D2": True, "D7": True})[1][1][0], dec, world.ref, world.run)
    assert kw["phi"] is None and kw["ratio"] is None and kw["alphas"] is None
    # refused and barred branches carry no note 6
    dec = make_decision(world, {"D10": True, "D11b": False})
    assert "note6" not in bs.resolve_spec(dec, RULES, "refused", "s2", HASH_OK)[1]
    assert "note6" not in bs.resolve_spec(dec, RULES, "barred", "s2", HASH_OK)[1]


def test_note8_takes_the_bet_off_s1_when_no_check_supports_the_gata4_factor(world):
    """rules.yaml note 8 (the user's decision): when the allowed factor table reaches its last row, s1 is WT cells (or
    the Mab21l2 factor when D11b passes) and s2 the Gata4 factor at exponent 0.5. Literal expectations."""
    scen = RULES["presets"]["gata4"]["scenario"]
    off = dict.fromkeys(EXTRA_KEYS, False)

    def comp(kind, l2):
        return {"kind": kind, "depleted_states": scen["depleted_states"], "spared_states": scen["spared_states"],
                "log2_magnitude": l2, "depletion_factor": 2.0 ** -l2, "magnitude_source": "Gata4 replicates (D7)"}

    def resolve(outcomes):
        dec = make_decision(world, outcomes)
        return dec, [bs.resolve_spec(dec, RULES, "allowed", s, HASH_OK) for s in ("s1", "s2")]

    dec, [(s1, i1), (s2, i2)] = resolve({})                    # nothing passes
    assert dec["submissions"]["allowed"]["s1"]["factor"]["source"] == "gata4_unvalidated"   # the table's choice
    assert s1 == {"factor": {"source": "none", "exponent": 0.0}, "composition": {"kind": "none"}, **off}
    assert s2 == {"factor": {"source": "gata4_unvalidated", "exponent": 0.5}, "composition": {"kind": "none"}, **off}
    assert i1["note8"]["applies"] and i1["note"].startswith("rules.yaml note 8") and "plain WT" in i1["warning"]
    assert not i2["note6"]["applies"] and i2["note"].startswith("rules.yaml note 8")
    _, [(s1, i1), (s2, _)] = resolve({"D7": True})           # composition_s1 and composition_s2_other_factor
    assert s1["composition"] == comp("central", M) and s2["composition"] == comp("central", M)
    assert i1["warning"] is None
    _, [(s1, _), (s2, _)] = resolve({"D7": True, "D8": True})
    assert s1["composition"] == comp("central", M) and s2["composition"] == comp("mild", M / 2)
    _, [(s1, i1), (s2, _)] = resolve({"D11b": True, "D7": True})   # the Mab21l2 factor carries over: s1 uses it
    assert s1["factor"] == {"source": "mab21l2", "exponent": 1.0} and s1["composition"] == comp("central", M)
    assert s2["factor"] == {"source": "gata4_unvalidated", "exponent": 0.5} and i1["warning"] is None
    for outcomes in ({"D2": True}, {"D2u": True}, {"D10": True}):   # any factor check passing: the table stands
        dec, [(s1, i1), (s2, _)] = resolve(outcomes)
        assert not i1["note8"]["applies"] and s1 == dec["submissions"]["allowed"]["s1"]
    dec, [(s1, _), (s2, _)] = resolve({})
    kw1, _ = bs.inputs_for(s1, dec, world.ref, world.run)
    kw2, _ = bs.inputs_for(s2, dec, world.ref, world.run)
    assert kw1["phi"] is None and np.array_equal(kw2["phi"], PHIS["gata4_pooled"]) and kw2["w"] == 0.5
    for branch in ("refused", "barred"):
        assert "note8" not in bs.resolve_spec(dec, RULES, branch, "s1", HASH_OK)[1]


def test_plan_lists_every_branch_and_builds_nothing(world, tmp_path):
    from typer.testing import CliRunner

    dec = make_decision(world, {"D7": True})
    path = write_decision(tmp_path / "gata4_decisions.json", dec)
    out = tmp_path / "plan.json"
    res = CliRunner().invoke(bs.app, ["plan", str(path), "--out", str(out)])
    assert res.exit_code == 0, res.output
    plan = json.loads(out.read_text())
    assert plan["decision_file"]["sha256"] == t3d.sha256(path) and plan["decision_file"]["protocol_as_registered"]
    assert plan["hashes"]["build_submission"] == t3d.sha256(HERE / "build_submission.py") and plan["created"]
    assert set(plan["branches"]) == set(bs.BRANCHES)
    assert all(set(v) == {"s1", "s2"} for v in plan["branches"].values())
    allowed = plan["branches"]["allowed"]
    assert allowed["s1"]["spec"]["factor"] == {"source": "none", "exponent": 0.0}
    assert allowed["s2"]["spec"]["factor"] == {"source": "gata4_unvalidated", "exponent": 0.5}
    assert plan["branches"]["barred"]["s1"]["spec"]["factor"]["source"] == "mab21l2"
    assert "note 8" in res.output
    assert sorted(p.name for p in tmp_path.iterdir()) == ["gata4_decisions.json", "plan.json"]   # nothing built
    bad = write_decision(tmp_path / "other.json", {**dec, "rules_sha256": "0" * 64})
    assert CliRunner().invoke(bs.app, ["plan", str(bad)]).exit_code == 2
    loose = write_decision(tmp_path / "loose.json", {**dec, "protocol_as_registered": False})
    res = CliRunner().invoke(bs.app, ["plan", str(loose)])
    assert res.exit_code == 0 and "registered protocol" in res.output       # shown, with a warning (build refuses)


def test_gata4_extras_follow_their_rules(world):
    dec = make_decision(world, {"D2": True, "D10": True, "D4": True, "D5": True, "D6": True})
    s1, _ = bs.resolve_spec(dec, RULES, "allowed", "s1", HASH_OK)
    s2, _ = bs.resolve_spec(dec, RULES, "allowed", "s2", HASH_OK)
    assert s1["factor"] == {"source": "gata4", "exponent": W1[0]} and all(s1[k] for k in bs.EXTRAS)
    assert s2["factor"]["source"] == "reference" and not any(s2[k] for k in bs.EXTRAS)
    for spec in (s1, s2):
        kw, info = bs.inputs_for(spec, dec, world.ref_stages, world.run)
        check_honoured(world, world.ref_stages, dec, "allowed", spec, kw)
    assert info["factor"]["recorded_in"].startswith("decision file")
    # depth alone
    dec = make_decision(world, {"D2": True, "D4": True})
    spec, _ = bs.resolve_spec(dec, RULES, "allowed", "s1", HASH_OK)
    assert spec["depth"] and not spec["stage"] and not spec["coverage"]


def test_barred_branch_is_rules_note_3(world):
    dec = make_decision(world, {"D2": True, "D10": True, "D7": True, "D7m": True})
    s1, info = bs.resolve_spec(dec, RULES, "barred", "s1", HASH_OK)
    s2, _ = bs.resolve_spec(dec, RULES, "barred", "s2", HASH_OK)
    assert info["source"].startswith("rules.yaml note 3")
    assert s1["factor"] == s2["factor"] == {"source": "mab21l2", "exponent": 1.0}
    assert s1["composition"]["kind"] == "central" and np.isclose(s1["composition"]["log2_magnitude"], M7M)
    assert s2["composition"]["kind"] == "mild" and np.isclose(s2["composition"]["log2_magnitude"], M7M / 2)
    for spec in (s1, s2):
        kw, _ = bs.inputs_for(spec, dec, world.ref, world.run)
        check_honoured(world, world.ref, dec, "barred", spec, kw)
        assert np.array_equal(kw["phi"], PHIS["transfer_mab21l2_all"])
    # D7m without a counted state: m = 0, so s2 would equal s1 and takes half the exponent (identical_fallback)
    dec = make_decision(world, {"D7m": False})
    s1, _ = bs.resolve_spec(dec, RULES, "barred", "s1", HASH_OK)
    s2, info = bs.resolve_spec(dec, RULES, "barred", "s2", HASH_OK)
    assert s1["composition"]["log2_magnitude"] == 0.0 and s2["factor"]["exponent"] == 0.5 and info["note"]


def test_smoke_decision_without_table_is_derived(world):
    dec = make_decision(world, {"D7m": True}, preset="mab21l2", table=False)
    derived = bs.roundtrip(t3d.decide_submissions(RULES, dec["rules"], None, {"submission_parameters": {}}))
    for branch in ("allowed", "refused"):
        for s in ("s1", "s2"):
            spec, info = bs.resolve_spec(dec, RULES, branch, s, HASH_OK)
            assert info.get("table_source", info["source"]).startswith("derived")
            assert spec == (allowed_from_text({**dec, "submissions": derived}, s) if branch == "allowed"
                            else derived[branch][s])
    spec, _ = bs.resolve_spec(dec, RULES, "allowed", "s2", HASH_OK)    # last row: note 8's half-strength bet
    assert spec["factor"] == {"source": "gata4_unvalidated", "exponent": 0.5}
    kw, info = bs.inputs_for(spec, dec, world.ref, world.run_smoke)
    assert np.array_equal(kw["phi"], PHIS["mab21l2_pooled"]) and info["factor"]["recorded_in"].endswith(
        "factors.json")
    # a gata4 decision file must carry its table
    with pytest.raises(bs.Refused, match="without a submissions table"):
        bs.resolve_spec(make_decision(world, table=False), RULES, "allowed", "s1", HASH_OK)


def test_changed_table_refused_unless_hash_mismatch_accepted(world):
    dec = make_decision(world, {"D2": True})
    dec["submissions"]["allowed"]["s1"]["factor"]["source"] = "reference"
    with pytest.raises(bs.Refused, match="changed after the run"):
        bs.resolve_spec(dec, RULES, "allowed", "s1", HASH_OK)
    spec, info = bs.resolve_spec(dec, RULES, "allowed", "s1", {"mismatch": ["rules"]})
    assert spec["factor"]["source"] == "reference" and info["recomputed_matches"] is False


def test_unknown_branch_or_submission_refused(world, tmp_path):
    dec = make_decision(world)
    for branch, s in (("maybe", "s1"), ("Allowed", "s1"), ("allowed", "s3")):
        with pytest.raises(bs.Refused, match="unknown"):
            bs.resolve_spec(dec, RULES, branch, s, HASH_OK)
    path = write_decision(tmp_path / "dec.json", dec)
    with pytest.raises(bs.Refused, match="unknown branch"):
        bs.build_submission(path, "maybe", "s1", world.wt, world.board, 600, tmp_path / "x.h5ad", run_dir=world.run)
    assert not (tmp_path / "x.h5ad").exists()


def test_validate_spec_rejects_unknowns():
    good = {"factor": {"source": "gata4", "exponent": 1.0}, "composition": {"kind": "none"}}
    bs.validate_spec(good, RULES)
    for bad in ({**good, "factor": {"source": "literature", "exponent": 1.0}},
                {**good, "composition": {"kind": "strong"}},
                {**good, "composition": {**bs.composition_spec("central", 1.0, RULES, "x"), "depleted_states": ["X"]}},
                {**good, "composition": {**bs.composition_spec("central", 1.0, RULES, "x"), "log2_magnitude": None}},
                {**good, "factor": {"source": "gata4", "exponent": None}},
                {**good, "depth": "yes"}):
        with pytest.raises(bs.Refused):
            bs.validate_spec(bad, RULES)


# ------------------------------------------------------------------------------------------------- provenance checks

def test_hash_mismatch_refused_unless_allowed_and_recorded(world, tmp_path):
    for field, what in (("rules_sha256", "rules"), ("script_sha256", "script")):
        path = write_decision(tmp_path / f"dec_{what}.json", make_decision(world, {"D2": True}, **{field: "0" * 64}))
        with pytest.raises(bs.Refused, match=what):
            bs.build_submission(path, "allowed", "s1", world.wt, world.board, 600, tmp_path / f"{what}.h5ad",
                                run_dir=world.run, threads=2)
        side = bs.build_submission(path, "allowed", "s1", world.wt, world.board, 600, tmp_path / f"{what}.h5ad",
                                   run_dir=world.run, allow_hash_mismatch=True, threads=2)
        assert side["hashes"]["mismatch"] == [what] and side["hashes"]["allow_hash_mismatch"] is True
        assert side["hashes"][f"decision_{field}"] == "0" * 64
        assert json.loads((tmp_path / f"{what}.json").read_text())["hashes"]["allow_hash_mismatch"] is True


def test_factor_file_must_match_its_record(world, tmp_path):
    run = tmp_path / "run"
    shutil.copytree(world.run, run)
    df = pd.read_csv(run / "factors" / "gata4_pooled.csv")
    df.loc[0, "factor"] = 3.0
    df.to_csv(run / "factors" / "gata4_pooled.csv", index=False)
    dec = make_decision(world, {"D2": True})
    spec, _ = bs.resolve_spec(dec, RULES, "allowed", "s1", HASH_OK)
    with pytest.raises(bs.Refused, match="differs from the one recorded"):
        bs.inputs_for(spec, dec, world.ref, run)
    with pytest.raises(bs.Refused, match="pass --run-dir"):
        bs.inputs_for(spec, dec, world.ref, None)


def test_wt_file_must_match_the_run_record(world, tmp_path):
    run = tmp_path / "run"
    shutil.copytree(world.run, run)
    fj = run / "factors" / "factors.json"
    rec = json.loads(fj.read_text())
    rec["inputs"]["wt"]["sha256"] = "f" * 64
    fj.write_text(json.dumps(rec))
    path = write_decision(tmp_path / "dec.json", make_decision(world, {"D2": True}))
    with pytest.raises(bs.Refused, match="WT file differs"):
        bs.build_submission(path, "allowed", "s1", world.wt, world.board, 600, tmp_path / "x.h5ad", run_dir=run,
                            threads=2)


def test_cell_count_protocol_and_preset_checks(world, tmp_path):
    path = write_decision(tmp_path / "dec.json", make_decision(world, {"D2": True}))
    for n in (100, 5000):
        with pytest.raises(bs.Refused, match="outside the board's limits"):
            bs.preflight(path, "allowed", "s1", world.board, n)
    with pytest.raises(bs.Refused, match="another preset"):
        bs.preflight(path, "allowed", "s1", world.board, 600, preset="mab21l2")
    assert bs.preflight(path, "barred", "s1", world.board, 600, preset="mab21l2")["preset"] == "mab21l2"
    path = write_decision(tmp_path / "dec_u.json", make_decision(world, {"D2": True}, protocol_as_registered=False))
    with pytest.raises(bs.Refused, match="registered protocol"):
        bs.preflight(path, "allowed", "s1", world.board, 600)
    pf = bs.preflight(path, "allowed", "s1", world.board, 600, allow_unregistered_protocol=True)
    assert any("non-registered" in w for w in pf["warnings"])
    t1 = replace(world.board, task="T1", key="T1:x")
    with pytest.raises(bs.Refused, match="T1 board"):
        bs.preflight(write_decision(tmp_path / "dec2.json", make_decision(world)), "allowed", "s1", t1, 600)


# ------------------------------------------------------------------------------------------------- construction

def test_construct_is_t3d_build_and_tracks_source_rows(world):
    ref = world.ref
    phi = PHIS["gata4_pooled"]
    a = bs.construct(ref, 700, t3d.stable_rng(1, "x"), f=ref.V.p, phi=phi, w=0.75)
    b = t3d.build(ref.ctx, 700, t3d.stable_rng(1, "x"), f=ref.V.p, phi=phi, w=0.75)
    assert np.array_equal(a.X, b)
    assert set(a.source) == {"wt"} and np.array_equal(a.lab, ref.pool.lab[a.row])
    # no factor, no depth: every cell is its source WT cell
    c = bs.construct(ref, 700, t3d.stable_rng(2, "x"), f=ref.V.p)
    assert np.array_equal(c.X, t3d.lognorm(ref.pool.C[c.row]))
    # the factor moves the genes it names, in count space
    d = a.X.mean(0) - t3d.lognorm(ref.pool.C[a.row]).mean(0)
    assert d[:10].mean() > 0.1 and d[10:20].mean() < -0.1


def test_composition_and_stage_rows_follow_the_inputs(world):
    ref = world.ref_stages
    spec = {"factor": {"source": "none", "exponent": 0.0},
            "composition": bs.composition_spec("central", 1.0, RULES, "test"), "stage": True}
    kw, info = bs.prediction_inputs(ref, spec, PARAMS, None)
    n = 1000
    built = bs.construct(ref, n, t3d.stable_rng(3, "y"), **kw)
    k = t3d.alloc(n, kw["f"])
    assert np.array_equal(np.bincount(built.lab, minlength=ref.V.L), k)
    assert info["states_missing"] == [] and np.isclose(info["depleted_share"] / info["depleted_share_wt"],
                                                       0.5 / (1 - 0.5 * info["depleted_share_wt"]))
    for t in np.flatnonzero(k):
        e95 = int(np.sum((built.lab == t) & (built.source == "stage:E9.5")))
        e825 = int(np.sum((built.lab == t) & (built.source == "stage:E8.25_late")))
        assert e95 == int(round(k[t] * 0.2)) and e825 == min(int(round(k[t] * 0.1)), int(k[t]) - e95)


def test_assemble_orders_panel_and_moves_duplicate_positions(world):
    from sklearn.neighbors import NearestNeighbors

    ref = world.ref
    built = bs.construct(ref, 2000, t3d.stable_rng(4, "z"), f=ref.V.p)       # > pool size: rows repeat
    panel = GENES[::-1]
    X, coords, info = bs.assemble(ref, built, panel, t3d.stable_rng(4, "z:assemble"))
    assert X.shape == (2000, G) and np.unique(coords, axis=0).shape[0] == 2000
    assert info["duplicate_positions_moved"] > 0 and info["jitter_radius"] > 0
    assert np.allclose(np.sort(X.sum(1)), np.sort(built.X.sum(1)), atol=1e-4)   # same cells, new order
    assert np.isclose(X[:, 0].sum(), built.X[:, G - 1].sum(), rtol=1e-5)        # columns in panel order
    # every cell sits on its source WT cell's section, within the jitter radius of a WT position in xy
    assert np.isin(coords[:, 2], ref.coords[:, 2]).all()
    d, _ = NearestNeighbors(n_neighbors=1).fit(ref.coords).kneighbors(coords.astype(np.float64))
    assert d.max() <= info["jitter_radius"] + 1e-3
    assert info["jitter_radius"] <= bs.JITTER_NN * bs.median_nn(ref.coords) + 1e-12
    with pytest.raises(bs.Refused, match="not in the WT reference"):
        bs.assemble(ref, built, [*GENES, "Extra"], t3d.stable_rng(0, "z"))


# ------------------------------------------------------------------------------------------------- files

def test_writes_a_valid_submission_and_sidecar(world, tmp_path):
    from vec.io import check_submission

    dec = make_decision(world, {"D2": True, "D10": True, "D7": True})
    path = write_decision(tmp_path / "gata4_decisions.json", dec)
    out = tmp_path / "bcat_allowed_s1.h5ad"
    side = bs.build_submission(path, "allowed", "s1", world.wt, world.board, 2500, out, run_dir=world.run,
                               threads=2)
    rep = check_submission(out, world.board)
    assert rep.ok, rep.errors
    assert rep.n_cells == 2500 and not any("repeat the exact coordinates" in w for w in rep.warnings)
    s = json.loads(out.with_suffix(".json").read_text())
    assert s == bs.roundtrip(side)
    assert s["decision_file"]["sha256"] == t3d.sha256(path) and s["branch"] == "allowed" and s["submission"] == "s1"
    assert s["hashes"]["build_submission_sha256"] == t3d.sha256(HERE / "build_submission.py")
    assert s["hashes"]["t3_decompose_on_disk"] == SCRIPT_SHA and s["hashes"]["mismatch"] == []
    assert s["inputs"]["wt"]["sha256"] == t3d.sha256(world.wt) and s["inputs"]["n_cells"] == 2500
    assert s["inputs"]["factor"]["sha256"] == world.params["factor_files"]["gata4_pooled"]["sha256"]
    assert s["inputs"]["board"]["key"] == "T3:synthetic" and s["output"]["check"]["ok"]
    assert s["spec"] == dec["submissions"]["allowed"]["s1"] and s["spec"]["composition"]["kind"] == "central"
    assert s["output"]["sha256"] == t3d.sha256(out)
    with pytest.raises(bs.Refused, match="exists"):
        bs.build_submission(path, "allowed", "s1", world.wt, world.board, 2500, out, run_dir=world.run, threads=2)


def test_cli_exit_codes(world, tmp_path):
    from typer.testing import CliRunner

    from vec.io import check_submission

    path = write_decision(tmp_path / "dec.json", make_decision(world, {"D10": True}))
    base = ["build", str(path), "--submission", "s2", "--wt", str(world.wt), "--board", "T3:synthetic",
            "--boards-dir", str(world.boards_dir), "--n-cells", "1500", "--run-dir", str(world.run), "--threads", "2"]
    runner = CliRunner()
    res = runner.invoke(bs.app, [*base, "--branch", "maybe", "--out", str(tmp_path / "a.h5ad")])
    assert res.exit_code == 2 and "unknown branch" in res.output
    res = runner.invoke(bs.app, [*base, "--branch", "refused", "--out", str(tmp_path / "b.h5ad")])
    assert res.exit_code == 0, res.output
    assert check_submission(tmp_path / "b.h5ad", world.board).ok
    side = json.loads((tmp_path / "b.json").read_text())
    assert side["spec"]["factor"]["source"] == "reference" and side["branch"] == "refused"


def test_compare_with_records(tmp_path):
    rng = np.random.default_rng(5)
    t3 = 80 + rng.normal(0, 0.7, 10)
    pd.DataFrame({"target": "mab21l2", "seed": range(10), "variant": "fac:factor", "t3": t3}).to_csv(
        tmp_path / "rec.csv", index=False)
    c = bs.compare_with_records({s: t3[s] + 0.2 for s in range(10)}, tmp_path / "rec.csv", "fac:factor", None)
    assert np.isclose(c["diff_mean"], 0.2) and c["within_seed_noise"] and c["target"] == "mab21l2"
    c = bs.compare_with_records({s: t3[s] + 5 for s in range(10)}, tmp_path / "rec.csv", "fac:factor", None)
    assert not c["within_seed_noise"]


# ------------------------------------------------------------------------------------------------- end to end (DGX)

def _write_released_like(path: Path, labels: list[str], prof: dict, n: int, rng, depth: float, factor=None,
                         cm: dict | None = None) -> None:
    import anndata as ad
    from scipy import sparse

    from vec.prep import log1p_cp10k

    lab = rng.choice(len(labels), n)
    rows = [rng.poisson(prof[labels[i]] * (1.0 if factor is None else factor) / prof[labels[i]].sum()
                        * rng.lognormal(np.log(depth), 0.3)) for i in lab]
    C = np.asarray(rows, dtype=np.int32)
    C[C.sum(1) == 0, 0] = 1
    obs = pd.DataFrame({"celltype": [labels[i] for i in lab]}, index=[f"{path.stem}_{i}" for i in range(n)])
    if cm is not None:
        obs["cm_celltype"] = [cm.get(labels[i], "Unknown") for i in lab]
    a = ad.AnnData(X=log1p_cp10k(C), obs=obs, var=pd.DataFrame(index=GENES))
    a.layers["counts"] = sparse.csr_matrix(C)
    a.obsm["spatial_3D"] = np.c_[rng.uniform(0, 500, n), rng.uniform(0, 500, n), rng.integers(0, 20, n) * 10.0]
    a.write_h5ad(path)


@pytest.mark.skipif(os.environ.get("RUN_T3_E2E") != "1", reason="set RUN_T3_E2E=1 (run on the DGX)")
def test_end_to_end_factors_and_build(world, tmp_path, monkeypatch):
    """`factors` on synthetic released-like files, then `build` from a smoke decision file (allowed s1: the
    knockout's own factor) and the barred branch on the gata4 preset; the stage-pool reference (build_context) holds
    the same WT pool as the light one."""
    monkeypatch.setattr(bs, "DATA", tmp_path / "no_released_data")     # never default to a real released file
    rng = np.random.default_rng(1)
    wt95 = ["V-CM", "A-CM", "Peri", "SE", "Foregut", "NCC", "ST", "Proepi", "aPHM", "pPHM", "Endo"]
    e875 = ["V-CM", "aPHM", "pPHM", "aSHF", "pSHF", "NCC", "SE", "Foregut", "Peri", "Endo", "Neural Tube"]
    prof = {k: rng.lognormal(0, 1.0, G) for k in sorted(set(wt95) | set(e875))}
    raw = tmp_path / "raw"
    raw.mkdir()
    cm = {"V-CM": "OFT/RV-CM"}
    _write_released_like(raw / "E9.5.h5ad", wt95, prof, 2400, rng, 80)
    _write_released_like(raw / "E8.75.h5ad", e875, prof, 2400, rng, 130, cm=cm)
    _write_released_like(raw / "E8.25_late.h5ad", e875, prof, 2000, rng, 200, cm=cm)
    _write_released_like(raw / "ko.h5ad", wt95, prof, 2400, rng, 95, factor=rng.lognormal(0, 0.4, G))
    rec = bs.write_factors("mab21l2", bs.preset_paths("mab21l2", raw / "E9.5.h5ad", raw / "ko.h5ad", None,
                                                       raw / "E8.75.h5ad", raw / "E8.25_late.h5ad"),
                           tmp_path / "smoke", threads=2)
    assert {"mab21l2_pooled", "reference_R", "transfer_mab21l2_all"} <= set(rec["files"])
    smoke = make_decision(world, {"D7m": True}, preset="mab21l2", table=False)
    path = write_decision(tmp_path / "mab21l2_decisions.json", smoke)
    side = bs.build_submission(path, "allowed", "s2", raw / "E9.5.h5ad", world.board, 2000, tmp_path / "s2.h5ad",
                               run_dir=tmp_path / "smoke", threads=2)        # last row: note 8's s2
    assert side["output"]["check"]["ok"] and side["inputs"]["factor"]["name"] == "mab21l2_pooled"
    assert side["spec"]["factor"]["exponent"] == 0.5 and side["spec_source"]["note8"]["applies"]
    # gata4 preset without any Gata4 file: transfer factor and R on WT E8.75; then stage mixing through build_context
    rec = bs.write_factors("gata4", bs.preset_paths("gata4", raw / "E8.75.h5ad", None, raw / "ko.h5ad", None,
                                                     raw / "E8.25_late.h5ad", raw / "E9.5.h5ad"),
                           tmp_path / "g4", threads=2)
    assert "transfer_mab21l2_vs95" in rec["files"]
    side = bs.build_submission(path, "barred", "s1", raw / "E8.75.h5ad", world.board, 2000, tmp_path / "b1.h5ad",
                               run_dir=tmp_path / "g4", preset="gata4", threads=2)
    assert side["output"]["check"]["ok"] and side["spec"]["factor"]["source"] == "mab21l2"
    ref = bs.load_reference(raw / "E8.75.h5ad", "gata4", RULES, {"E9.5": raw / "E9.5.h5ad",
                                                                 "E8.25_late": raw / "E8.25_late.h5ad"}, threads=2)
    light = bs.load_reference(raw / "E8.75.h5ad", "gata4", RULES, threads=2)
    assert ref.V.labels == light.V.labels and np.array_equal(ref.pool.C, light.pool.C)
    assert np.array_equal(ref.pool.lab, light.pool.lab) and [s.name for s in ref.stages] == ["E9.5", "E8.25_late"]
