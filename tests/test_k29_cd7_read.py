"""Try K29 (experiments/harness/TRIES.md; experiments/t3_decomp/k29_cd7_read.py and run_k29.sh): the practice read that
T3's file cD7 on C2 is bound to. Synthetic vocabularies and pools, hand-made score tables, hand-made task files and
hand-made K18 files whose pins are computed here, a hand-made ledger, a fixed clock and a throwaway git repository for
the launcher. No competition file, no score of a fold, no read of the committed ledger and no pin of bytes built from
released data. The frozen rules.yaml is read for its two state lists."""

from __future__ import annotations

import functools
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from vec import harness
from vec.adoption import SEEDS

REPO = Path(__file__).resolve().parents[1]
T3 = REPO / "experiments" / "t3_decomp"
if str(T3) not in sys.path:
    sys.path.insert(0, str(T3))
import k7_wt_revert as k7  # noqa: E402
import k18_built_files as k18  # noqa: E402
import k21_factor_forms as k21  # noqa: E402
import k29_cd7_read as k29  # noqa: E402

SCEN = ["OFT/RV-CM", "aSHF", "pSHF", "aPHM", "pPHM"]
OWN = ["vCM1", "vCM2", "aCM1", "aCM2", "OFT/RV-CM", "JCF/SV-CM", "FHF/JCF", "V-CM", "IFT-CM"]
LISTS = {"own": OWN, "oth": SCEN}
LABELS = ["vCM1", "vCM2", "aCM1", "aCM2", "OFT/RV-CM", "JCF/SV-CM", "FHF/JCF", "aSHF", "pSHF", "aPHM", "pPHM", "Endo",
          "Blood"]  # seven of the nine own labels: V-CM and IFT-CM absent, as in the gata4 run's vocabulary
SPEC = {"factor": {"source": "mab21l2", "exponent": 1.0},
        "composition": {"kind": "central", "depleted_states": SCEN,
                        "spared_states": ["vCM1", "vCM2", "JCF/SV-CM", "FHF/JCF"], "log2_magnitude": k29.L2_F,
                        "depletion_factor": 2.0 ** -k29.L2_F, "magnitude_source": "Mab21l2 (D7m)"},
        "depth": False, "stage": False, "coverage": False}
ALLOWED = {"kind": "central", "depleted_states": SCEN, "log2_magnitude": k29.L2_D, "depletion_factor": k29.DEPLETION_D}
DEC = {"submissions": {"allowed": {"s1": {"composition": ALLOWED}}}}
RULES = {"presets": {"gata4": {"cm_states": OWN, "scenario": {"depleted_states": SCEN}}},
         "protocol": {"seeds": list(SEEDS)}}
ROW = "| K29 (one row, one script) | T3 | `experiments/t3_decomp/k29_cd7_read.py` with its tests |\n"
G, N_POOL = 12, 80
GENES = [f"g{i}" for i in range(G)]
BASE = 56.0
GOOD = {"own_f": 0.5, "own_d": 1.0}  # G +0.5 and Gd +1.0 in every cell, C 0: nothing strikes and nothing is labelled
BEFORE = datetime(2026, 11, 24, 12, 0, tzinfo=UTC)
RESULTS_AS_SHIPPED = k29.RESULTS  # before the autouse fixture points it at an empty folder


@pytest.fixture(autouse=True)
def clean_state(monkeypatch, tmp_path):
    monkeypatch.setattr(k29, "_now", lambda: BEFORE)  # no test reads the machine's clock
    monkeypatch.setattr(k29, "RESULTS", tmp_path / "k29_results_none")  # nor meets a reading on disk
    yield
    k29._K29.clear()
    k18._FOLD.clear()


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _vocab(labels=LABELS, seed=0):
    t3d = harness._t3_head().t3d
    count = np.random.default_rng(seed).integers(20, 200, len(labels))
    L = len(labels)
    return t3d.Vocab(labels=list(labels), group=np.array(["g"] * L), mesp1=np.array(["partial"] * L),
                     nm=np.zeros(L, dtype=bool), broad=np.ones(L, dtype=bool), model_idx=np.arange(L),
                     p=count / count.sum(), count=count)


def _hand(V, states, l2) -> np.ndarray:
    """The composition by hand: each listed state the vocabulary holds times 2^-m, then renormalised."""
    f = np.asarray(V.p, dtype=np.float64).copy()
    for t in states:
        if t in V.labels:
            f[V.labels.index(t)] *= 2.0 ** -l2
    return f / f.sum()


@functools.cache
def _comps_json() -> str:
    return json.dumps(k29.arm_compositions(harness._t3_head(), _vocab(), SPEC, LISTS)["figures"])


def _comps() -> dict:
    return json.loads(_comps_json())


# ------------------------------------------------------------------------------------------------ constants, pins


def test_arms_constants_outcome_lines_and_the_pins():
    assert k29.ARMS == ("oth_f", "own_h", "own_f", "own_d", "oth_h", "oth_d") and k29.GUARD_ARM == "oth_f"
    assert k29.COUNTED == ("own_h", "own_f", "own_d", "oth_h", "oth_d")
    assert (k29.L2_H, k29.L2_F, k29.L2_D) == (0.358906935493353, 0.717813870986706, 1.3826521560990648)
    assert k29.L2_H == k29.L2_F / 2 and [round(2.0 ** -x, 4) for x in (k29.L2_H, k29.L2_F, k29.L2_D)] == [
        0.7798, 0.608, 0.3835]
    assert k29.DEPLETION_D == 0.3835131214774091 == 2.0 ** -k29.L2_D and k29.GATE_REL == 1e-12
    assert k29.ARM_OF == {"own_h": ("own", k29.L2_H), "own_f": ("own", k29.L2_F), "own_d": ("own", k29.L2_D),
                          "oth_h": ("oth", k29.L2_H), "oth_f": ("oth", k29.L2_F), "oth_d": ("oth", k29.L2_D)}
    assert k29.STEPS == {"G": ("own_f", "own_h"), "Gd": ("own_d", "own_h"), "C": ("oth_d", "oth_h")}
    assert k29.GATED_PAIRS == (("own_d", "own_h"), ("oth_d", "oth_h")) and k29.FLIP_PAIR == ("own_d", "own_h")
    assert k29.OWN_ABSENT == ("V-CM", "IFT-CM") and k29.OWN_FOUND == 7
    assert (k29.MIN_GAIN, k29.MIN_AHEAD, k29.MAX_BEHIND, k29.MIN_CELLS, k29.MAX_EXTRA_GATED) == (0.3, 8, -0.3, 2, 1)
    assert (k29.GUARD_MEAN, k29.GUARD_SEED) == (0.1, 0.5)
    assert k29.CUTOFF == datetime(2026, 11, 25, 6, 0, tzinfo=UTC) and k29.CUTOFF_TEXT == "2026-11-25 06:00 UTC"
    assert (k29.EXIT_REFUSED, k29.EXIT_INCOMPLETE, k29.EXIT_EXISTS, k29.EXIT_GATE, k29.EXIT_GUARD, k29.EXIT_LATE) == (
        2, 3, 4, 5, 6, 7)
    assert k29.STREAM == k18.STREAM == "t3_s2_form" and k29.TASKS == k18.TASKS and k29.WORLDS == k18.WORLDS
    assert k29.CELLS == ("plain g1", "plain g2", "dm g1", "dm g2") and list(SEEDS) == list(range(10))
    assert k29.VOID_LINE.format("x") == "Try K29: VOID: x; no file."
    assert k29.STRUCK_LINE.format("S1, S3") == "Try K29: STRUCK (S1, S3): no file; C2 stays."
    assert k29.NOT_STRUCK_LINE == ("Try K29: NOT STRUCK: the file may be built and uploaded under the captain's bound "
                                   "yes (ruling row R36).")
    assert k29.GATE_VOID.format("x") == "the target-free gate failed (x)" and k29.LATE == "INCOMPLETE past its cutoff"
    # K18's three committed task files: the pins rows K21 to K24 read them at, with the prefixes row K29 quotes
    assert k29.K18_SHA256 == k21.K18_SHA256
    for f, pre in (("plain_g1.json", "5ec09ef58386"), ("plain_g2.json", "cfc5ede7e718"), ("dm.json", "1ef3e26347c6")):
        assert k29.K18_SHA256[f].startswith(pre) and re.fullmatch(r"[0-9a-f]{64}", k29.K18_SHA256[f])
    assert k29.K18_TASKS == T3 / "results" / "k18_built_files" / "tasks"
    assert RESULTS_AS_SHIPPED == T3 / "results" / "k29_cd7_read"
    assert k29.SCRIPT == "experiments/t3_decomp/k29_cd7_read.py"


def test_the_lists_are_the_frozen_rules_and_the_magnitudes_are_checked_against_the_files():
    assert k29.state_lists(harness._t3_head().load_rules()) == LISTS  # nine own labels, the scenario's five
    assert k29.state_lists(RULES) == LISTS
    k29.check_spec(SPEC, DEC, LISTS)
    comp = SPEC["composition"]
    bad = [({**SPEC, "composition": {**comp, "log2_magnitude": 0.7178}}, DEC, LISTS),
           ({**SPEC, "composition": {**comp, "depleted_states": SCEN[:4]}}, DEC, LISTS),
           ({**SPEC, "composition": {"kind": "none"}}, DEC, LISTS),
           (SPEC, {"submissions": {"allowed": {"s1": {"composition": {**ALLOWED, "log2_magnitude": 1.3826}}}}}, LISTS),
           (SPEC, {"submissions": {"allowed": {"s1": {"composition": {**ALLOWED, "depletion_factor": 0.3835}}}}},
            LISTS),
           (SPEC, {"submissions": {"allowed": {"s1": {"composition": {"kind": "none"}}}}}, LISTS),
           (SPEC, {}, LISTS),
           (SPEC, DEC, {"own": OWN[:7], "oth": SCEN})]
    for spec, dec, lists_ in bad:
        with pytest.raises(k29.Refused):
            k29.check_spec(spec, dec, lists_)


# ------------------------------------------------------------------------------------------------ the target-free gate


def test_the_gate_passes_on_the_rows_vocabulary_and_prints_counts_and_ratios_only():
    res = k29.gate_eval(harness._t3_head(), _vocab(), DEC, RULES)
    assert res["passes"] and res["failed"] == [] and res["checks"] == {"a": True, "b": True, "c": True}
    assert res["scenario_found"] == 5 and res["scenario_missing"] == [] and res["uncut_labels"] == len(LABELS) - 5
    assert res["own_found"] == 7 and res["own_listed"] == 9 and res["own_missing"] == ["V-CM", "IFT-CM"]
    assert set(res["ratios"]) == set(SCEN)
    for v in res["ratios"].values():  # (f_i / p_i) / (f_j / p_j) over every uncut label j
        assert abs(v["min"] / k29.DEPLETION_D - 1) <= 1e-12 and abs(v["max"] / k29.DEPLETION_D - 1) <= 1e-12
    lines = k29.gate_lines(res)
    shown = "\n".join(lines)
    assert lines[0] == "K29 target-free gate: passes" and "found 5 of 5" in shown and "7 of 9 labels found" in shown
    for x in re.findall(r"\d+\.\d+(?:e-?\d+)?", shown):  # every decimal printed is a ratio: no share, no cell count
        assert abs(float(x) / k29.DEPLETION_D - 1) < 1e-9, x


@pytest.mark.parametrize("labels,dec,failed", [
    ([x for x in LABELS if x != "pPHM"], DEC, ["a", "b"]),  # a scenario state missing: (a), and (b) with it
    (LABELS + ["V-CM"], DEC, ["c"]),  # the own list finds eight labels
    ([x for x in LABELS if x != "vCM1"], DEC, ["c"]),  # six
    (LABELS, {"submissions": {"allowed": {"s1": {"composition": {**ALLOWED, "log2_magnitude": 1.38}}}}}, ["b"]),
    (LABELS, {"submissions": {"allowed": {"s1": {"composition": {**ALLOWED, "depletion_factor": 0.38}}}}}, ["b"]),
])
def test_each_gate_check_fails(labels, dec, failed):
    res = k29.gate_eval(harness._t3_head(), _vocab(labels), dec, RULES)
    assert not res["passes"] and res["failed"] == [k29.GATE_CHECKS[k] for k in failed]
    assert k29.gate_lines(res)[0].startswith("K29 target-free gate: FAILS (")


@pytest.fixture
def gate_env(monkeypatch, tmp_path):
    """cmd_gate with the preconditions passed, the decision file by hand and the builder's loader returning a synthetic
    vocabulary."""
    st = {"commit": "abc", "tree_clean": True, "machine": "m1", "utc": "x"}
    bs = harness._t3_head()
    env = SimpleNamespace(out=tmp_path / "run", V=_vocab(), loads=[])

    def load_reference(wt, preset, rules, threads=16):
        env.loads.append((Path(wt).name, preset, threads))
        return SimpleNamespace(V=env.V)

    monkeypatch.setattr(k29, "preconditions", lambda commit: ("K29", st))
    monkeypatch.setattr(k18, "barred_spec", lambda decision=None: (SPEC, DEC, RULES))
    monkeypatch.setattr(bs, "load_reference", load_reference)
    monkeypatch.setattr(bs, "sha256", lambda p: "s" * 64)
    return env


def test_the_gate_writes_gate_json_once_and_exits_5_on_a_failed_check(gate_env, capsys):
    argv = ["gate", "--commit", "abc", "--out", str(gate_env.out)]
    k29.main(argv)
    g = json.loads((gate_env.out / "gate.json").read_text())
    assert g["passes"] and g["row"] == "K29" and g["composition"] == ALLOWED and g["provenance"]["commit"] == "abc"
    assert gate_env.loads == [("E8.75.h5ad", "gata4", 8)]  # WT E8.75 through the builder's own loader
    assert capsys.readouterr().out.splitlines()[0] == "K29 target-free gate: passes"
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 4  # made once
    gate_env.V = _vocab([x for x in LABELS if x != "aPHM"])
    out2 = gate_env.out.parent / "run2"
    with pytest.raises(SystemExit) as e:
        k29.main(["gate", "--commit", "abc", "--out", str(out2)])
    g = json.loads((out2 / "gate.json").read_text())
    assert e.value.code == 5 and g["passes"] is False and g["failed"] == [k29.GATE_CHECKS["a"], k29.GATE_CHECKS["b"]]
    assert "FAILS (a: the five scenario states found" in capsys.readouterr().out


def test_a_crash_of_the_gate_is_masked_and_writes_nothing_and_a_refusal_writes_nothing(gate_env, monkeypatch, capsys):
    def boom(*a, **k):
        raise MemoryError("share 0.1234 of 5678 cells")

    monkeypatch.setattr(harness._t3_head(), "load_reference", boom)
    with pytest.raises(SystemExit) as e:
        k29.main(["gate", "--commit", "abc", "--out", str(gate_env.out)])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "MemoryError" in err and "0.1234" not in err and not gate_env.out.exists()
    comp = {**SPEC["composition"], "log2_magnitude": 0.75}
    monkeypatch.setattr(k18, "barred_spec", lambda decision=None: ({**SPEC, "composition": comp}, DEC, RULES))
    with pytest.raises(SystemExit) as e:
        k29.main(["gate", "--commit", "abc", "--out", str(gate_env.out)])
    assert e.value.code == 2 and "not the row's second" in capsys.readouterr().err and not gate_env.out.exists()


# ------------------------------------------------------------------------------------------------ the compositions


def test_each_arm_is_the_barred_spec_with_only_its_composition_block_changed():
    bs = harness._t3_head()
    V = _vocab()
    out = k29.arm_compositions(bs, V, SPEC, LISTS)
    rest = {k: v for k, v in SPEC.items() if k != "composition"}
    for a in k29.ARMS:
        which, l2 = k29.ARM_OF[a]
        s = out["specs"][a]
        assert {k: v for k, v in s.items() if k != "composition"} == rest
        assert (s["composition"]["kind"], s["composition"]["depleted_states"], s["composition"]["log2_magnitude"]) == (
            "central", LISTS[which], l2)
        assert np.allclose(out["f"][a], _hand(V, LISTS[which], l2), rtol=0, atol=1e-15)
        fig = out["figures"][a]
        assert fig["states_missing"] == (["V-CM", "IFT-CM"] if which == "own" else [])
        assert len(fig["states"]) == (7 if which == "own" else 5)
        for v in fig["states"].values():
            assert v["ratio_before_renormalising"] == pytest.approx(2.0 ** -l2, rel=1e-12)
            assert v["share"] == pytest.approx(v["wt_share"] * v["ratio_after"], rel=1e-12)
    f_barred, _ = bs.composition(V, SPEC, {})
    assert np.array_equal(out["f"]["oth_f"], f_barred)  # oth_f is the barred s1's composition value for value
    assert SPEC["composition"]["depleted_states"] == SCEN and SPEC["composition"]["log2_magnitude"] == k29.L2_F
    assert k29.composition_checks(out["figures"]) == []


@pytest.mark.parametrize("drop,add,why", [("pPHM", None, "oth_f: the scenario misses ['pPHM']"),
                                          ("aSHF", None, "oth_d: the scenario misses ['aSHF']"),
                                          (None, "V-CM", "own_h: the own list misses ['IFT-CM']"),
                                          ("vCM2", None, "own_d: the own list misses ['vCM2', 'V-CM', 'IFT-CM']")])
def test_the_arms_composition_checks_raise_before_anything_is_drawn(drop, add, why):
    labels = [x for x in LABELS if x != drop] + ([add] if add else [])
    with pytest.raises(RuntimeError, match=re.escape(why)):
        k29.arm_compositions(harness._t3_head(), _vocab(labels), SPEC, LISTS)


def test_the_arms_checks_read_each_cut_state_before_renormalising():
    figs = _comps()
    assert k29.composition_checks(figs) == []
    figs["oth_d"]["states"]["aSHF"]["ratio_before_renormalising"] = 0.3842
    assert "oth_d: a cut state is not at 2^-m of its wild-type share before renormalising" in \
        k29.composition_checks(figs)
    assert k29.composition_checks({})  # nothing recorded: every arm fails


# ------------------------------------------------------------------------------------------------ the builds


@pytest.fixture
def fold(monkeypatch):
    """A synthetic reference whose vocabulary holds the row's labels, the builder's draw replaced by a recorder (the
    keywords, the stream's state), K18's plain-world preparation stubbed as K18's own tests stub it."""
    bs = harness._t3_head()
    V = _vocab()
    rng = np.random.default_rng(5)
    C = rng.poisson(3.0, size=(N_POOL, G)).astype(np.int32)
    C[C.sum(1) == 0, 0] = 1
    lab = rng.integers(0, V.L, N_POOL)
    P = (C / C.sum(1, keepdims=True)).astype(np.float32)
    ref = bs.Reference(preset="gata4", rules={}, genes=GENES, V=V, pool=bs.t3d.Pool(C, lab, P, V.L),
                       coords=rng.normal(size=(N_POOL, 3)))
    calls = []

    def construct(r, n, g, f=None, phi=None, w=1.0, ratio=None, alphas=None):
        calls.append({"state": json.dumps(g.bit_generator.state, default=int), "f": np.asarray(f).copy(),
                      "phi": None if phi is None else np.asarray(phi).copy(), "w": w, "ratio": ratio,
                      "alphas": alphas, "n": n})
        rows = g.integers(0, N_POOL, n)
        return bs.Built(X=bs.t3d.lognorm(C[rows]), lab=lab[rows], source=np.full(n, "wt", dtype=object), row=rows)

    group = np.array(["A", "B", "Unknown"])[rng.integers(0, 3, N_POOL)]
    per = {"A": np.exp(rng.normal(0, 0.5, G)), "B": np.exp(rng.normal(0, 0.5, G))}
    phi = np.exp(rng.normal(0, 0.5, G))
    file_genes = GENES[::-1] + ["extra"]
    barred_file = np.exp(rng.normal(0, 0.5, G + 1))
    phi_b = k18.in_gene_order(barred_file, file_genes, GENES)
    pool = ad.AnnData(X=bs.t3d.lognorm(C).astype(np.float32), var=pd.DataFrame(index=GENES))
    prep = {"ref": ref, "phi": phi, "per": per, "group": group, "record": {}}
    monkeypatch.setattr(bs, "construct", construct)
    monkeypatch.setattr(bs, "assemble", lambda r, built, panel, g: (built.X, None, {}))
    monkeypatch.setattr(k7, "_prepared", lambda ctx: prep)
    monkeypatch.setattr(harness, "count_for", lambda board, n: 40)
    monkeypatch.setattr(harness, "finish", lambda ctx, p, X, coords: X)
    monkeypatch.setattr(k18, "_BARRED", {"phi": barred_file, "genes": file_genes, "spec": SPEC, "record": {}})
    monkeypatch.setattr(k29, "_K29", {"lists": LISTS})
    F = k29.fold_factors(bs, ref, per, phi, group, np.asarray(C.mean(0), dtype=float), phi_b, SPEC)

    def ctx(task: str, seed: int):
        return SimpleNamespace(task=SimpleNamespace(key=k18.TASKS[task].key, reference_role="wt"),
                               inputs={"wt": pool}, board=SimpleNamespace(genes=GENES, needs_coords=False), seed=seed)

    return SimpleNamespace(bs=bs, ref=ref, V=V, F=F, calls=calls, phi_b=phi_b, ctx=ctx)


def test_every_arm_is_k18s_barred_call_and_oth_f_is_barred_form_value_for_value(fold):
    want = json.dumps(fold.bs.t3d.stable_rng(3, k18.STREAM).bit_generator.state, default=int)
    for a in k29.ARMS:
        k29.build_arm(fold.bs, fold.ref, a, 40, 3, fold.F)
    k18.build_arm(fold.bs, fold.ref, k18.BARRED, 40, 3, fold.F)  # K18's own barred_form on the same fold
    *mine, theirs = fold.calls
    for a, c in zip(k29.ARMS, mine, strict=True):
        which, l2 = k29.ARM_OF[a]
        assert c["state"] == want and c["n"] == 40  # one stream, one seed: the same draw generator in every arm
        assert np.array_equal(c["phi"], fold.phi_b) and c["w"] == 1.0 and c["ratio"] is None and c["alphas"] is None
        assert np.allclose(c["f"], _hand(fold.V, LISTS[which], l2), rtol=0, atol=1e-15)
        if a != "oth_f":
            assert not np.array_equal(c["f"], theirs["f"])
    oth_f = mine[k29.ARMS.index("oth_f")]
    for k, v in oth_f.items():
        assert (np.array_equal(v, theirs[k]) if isinstance(v, np.ndarray) else v == theirs[k]), k
    with pytest.raises(RuntimeError, match="unknown arm"):
        k29.build_arm(fold.bs, fold.ref, k18.BARRED, 40, 3, fold.F)


def test_k18s_plain_build_draws_this_trys_arms_inside_the_context_only(fold):
    with k29.k18_arms(k29.ARMS):
        for a in k29.ARMS:
            k18._build_plain(fold.ctx("g1", 2), a)
        assert "k29" in k18._FOLD["F"] and set(k18._FOLD["F"]["k29"]["figures"]) == set(k29.ARMS)
    i = LABELS.index("aSHF")
    for a, c in zip(k29.ARMS, fold.calls, strict=True):
        assert c["f"][i] == pytest.approx(_hand(fold.V, LISTS[k29.ARM_OF[a][0]], k29.ARM_OF[a][1])[i], rel=1e-12)
    assert k18.ARMS == k29._K18_ORIGINAL["ARMS"] and k18._FOLD == {}


def test_k18s_names_are_replaced_inside_the_context_and_restored_on_exit_and_on_an_exception():
    orig = (k18.ARMS, k18.build_arm, k18.fold_factors)
    with k29.k18_arms(k29.ARMS):
        assert k18.ARMS == k29.ARMS and k18.build_arm is k29.build_arm and k18.fold_factors is k29.fold_factors
    assert (k18.ARMS, k18.build_arm, k18.fold_factors) == orig
    with pytest.raises(ValueError, match="boom"):
        with k29.k18_arms(("oth_f",)):
            assert k18.ARMS == ("oth_f",)
            raise ValueError("boom")
    assert (k18.ARMS, k18.build_arm, k18.fold_factors) == orig
    for bad in (("own_h",), ("oth_f", k18.BARRED)):
        with pytest.raises(RuntimeError), k29.k18_arms(bad):
            pass
    assert k29.part_arms("guard") == ("oth_f",) and k29.part_arms("score") == k29.ARMS
    with pytest.raises(RuntimeError, match="not loaded"):
        k29.fold_factors(None, None, {}, None, None, None, None, SPEC)


# ------------------------------------------------------------------------------------------------ the guard


def _refs() -> dict:
    """K18's committed barred_form rows by hand."""
    return {tuple(c.split()): [{"seed": s, "scores": {k18.BARRED: BASE + 0.1 * s}} for s in SEEDS] for c in k29.CELLS}


def _guard_rows(shift=None, one=None, drop=None) -> dict:
    out = {}
    for c in k29.CELLS:
        rows = [{"seed": s, "scores": {"oth_f": BASE + 0.1 * s + (shift or {}).get(c, 0.0)}} for s in SEEDS]
        if one and c in one:
            rows[3]["scores"]["oth_f"] += one[c]
        out[c] = rows[:9] if drop == c else rows
    return out


def test_the_guard_holds_within_k17s_tolerance_in_all_four_cells():
    refs = _refs()
    ok = k29.guard_eval(_guard_rows(), refs)
    assert ok["passes"] and ok["failed"] == [] and set(ok["cells"]) == set(k29.CELLS)
    assert all(v["distance"]["max_abs_diff"] == 0.0 for v in ok["cells"].values())
    assert k29.guard_eval(_guard_rows(shift={"dm g2": 0.09}, one={"plain g1": 0.45}), refs)["passes"]
    for kw in ({"shift": {"plain g2": 0.11}}, {"shift": {"dm g1": -0.11}}, {"one": {"dm g2": 0.6}},
               {"drop": "plain g1"}):
        res = k29.guard_eval(_guard_rows(**kw), refs)
        assert not res["passes"] and res["failed"] == ["oth_f"], kw
        assert sum(not v["passes"] for v in res["cells"].values()) == 1
    res = k29.guard_eval(_guard_rows(), {k: v for k, v in refs.items() if k != ("dm", "g1")})
    assert not res["passes"] and not res["cells"]["dm g1"]["passes"]  # every cell is needed


def test_the_run_steps_oth_f_is_matched_with_the_guard_steps():
    g = _guard_rows()
    assert k29.match_eval(_guard_rows(shift={"dm g1": 0.09}), g)["passes"]
    for kw in ({"shift": {"plain g1": 0.11}}, {"one": {"dm g2": -0.6}}, {"drop": "dm g1"}):
        res = k29.match_eval(_guard_rows(**kw), g)
        assert not res["passes"] and sum(not v["passes"] for v in res["cells"].values()) == 1, kw


# ------------------------------------------------------------------------------------------------ the rule


def _table(offsets: dict | None = None) -> dict[str, list[float]]:
    """Every arm at BASE + 0.1 seed, plus each named arm's offset (a number, or ten numbers by seed)."""
    t = {a: [BASE + 0.1 * s for s in SEEDS] for a in k29.ARMS}
    for a, v in (offsets or {}).items():
        t[a] = [BASE + 0.1 * s + x for s, x in zip(SEEDS, v if isinstance(v, list) else [v] * len(SEEDS),
                                                    strict=True)]
    return t


def _state(gated: dict | None, arm: str, seed: int) -> str:
    """An arm's severity state on a seed: `gated` maps the arm to seeds (r2_gate) or to {seed: state}."""
    g = (gated or {}).get(arm, ())
    if isinstance(g, dict):
        return g.get(seed, "scored")
    return "r2_gate" if seed in g else "scored"


def _row_set(table: dict[str, list[float]], gated: dict | None = None) -> list[dict]:
    out = []
    for i, s in enumerate(SEEDS):
        sc = {a: float(v[i]) for a, v in table.items()}
        out.append({"seed": s, "scores": sc,
                    "skills": {a: {"de_score": round(sc[a] / 100.0, 6), "mmd_u": 0.5} for a in sc},
                    "t3": {a: {"severity": {"state": _state(gated, a, s), "beta": 0.1, "r2": 0.05}} for a in sc}})
    return out


def _cell_table(c: str, offsets: dict) -> dict:
    t = _table(offsets)
    if c.startswith("plain"):
        t[harness.FLOOR] = [50.0 + 0.1 * s for s in SEEDS]
    return t


def _per(cells: dict | None = None, gated: dict | None = None) -> dict:
    per = {}
    for c in k29.CELLS:
        t = _cell_table(c, (cells or {}).get(c, GOOD))
        per[c] = k29.summarise(_row_set(t, (gated or {}).get(c)), tuple(t))
    return per


PLAIN, DM = ("plain g1", "plain g2"), ("dm g1", "dm g2")


def test_nothing_strikes_and_the_table_and_reported_pairs():
    per = _per()
    dec = k29.decide(per)
    assert dec["clauses"] == [] and not dec["struck"] and dec["outcome"] == k29.NOT_STRUCK_LINE
    assert all(dec["table"]["Gd"][c]["AHEAD"] and dec["table"]["G"][c]["AHEAD"] for c in k29.CELLS)
    assert not any(dec["table"]["C"][c]["AHEAD"] or dec["table"]["C"][c]["BEHIND"] for c in k29.CELLS)
    assert dec["label"]["statement"] is None and dec["c_four_cell_mean"] == pytest.approx(0.0, abs=1e-9)
    assert set(per["plain g1"]["paired"]) == {"G", "Gd", "C", "own_h - oth_f", "own_f - oth_f", "own_d - oth_f",
                                              "oth_h - oth_f", "oth_d - oth_f", "own_d - own_f", "oth_f - floor"}
    assert "oth_f - floor" not in per["dm g1"]["paired"]
    assert per["dm g2"]["paired"]["own_d - own_f"]["mean"] == pytest.approx(0.5)
    assert per["dm g2"]["per_term"]["Gd"]["de_score"] == pytest.approx(0.01)
    assert per["dm g1"]["flips"] == {"seeds": [], "flip_free_seeds": 10, "flip_free_mean_gd": pytest.approx(1.0)}
    lines = k29.table_lines(dec)
    for s in ("G = own_f - own_h: ", "Gd = own_d - own_h: ", "C = oth_d - oth_h: "):
        assert any(x.startswith(f"  {s}") for x in lines), s


@pytest.mark.parametrize("cells,clauses", [
    ({c: {**GOOD, "own_d": 0.0} for c in PLAIN}, ["S1"]),  # Gd AHEAD in no plain cell
    ({c: {**GOOD, "own_d": 0.0} for c in DM}, ["S1"]),  # Gd AHEAD in no dm cell
    ({c: {**GOOD, "own_d": [2.0] * 7 + [-0.5] * 3} for c in PLAIN}, ["S1"]),  # mean 1.25, 7 of 10 positive: not AHEAD
    ({c: {**GOOD, "own_d": 0.29} for c in DM}, ["S1"]),  # a mean not above 0.3: not AHEAD
    ({c: {**GOOD, "own_f": -0.5} for c in ("plain g2", "dm g1")}, ["S2"]),  # G BEHIND in two cells
    ({c: {**GOOD, "own_f": -0.5} for c in k29.CELLS}, ["S2"]),  # G BEHIND in four cells
    ({c: {**GOOD, "oth_d": -0.5} for c in ("plain g1", "dm g2")}, ["S3"]),  # C BEHIND in two cells, mean -0.25
    ({"plain g1": {**GOOD, "oth_d": -1.3}, **{c: {**GOOD, "oth_d": -0.1} for c in ("plain g2", *DM)}}, ["S3"]),
])
def test_each_strike_clause(cells, clauses):
    dec = k29.decide(_per(cells))
    assert dec["clauses"] == clauses and dec["struck"]
    assert dec["outcome"] == f"Try K29: STRUCK ({', '.join(clauses)}): no file; C2 stays."


def test_s3_on_the_four_cell_mean_alone_and_the_edges_of_each_clause():
    cells = {"plain g1": {**GOOD, "oth_d": -1.3}, **{c: {**GOOD, "oth_d": -0.1} for c in ("plain g2", *DM)}}
    dec = k29.decide(_per(cells))
    assert dec["c_four_cell_mean"] == pytest.approx(-0.4) and sum(v["BEHIND"] for v in dec["table"]["C"].values()) == 1
    # AHEAD with exactly 8 of 10 seeds positive and a mean above 0.3
    assert k29.decide(_per({c: {**GOOD, "own_d": [2.0] * 8 + [-0.5] * 2} for c in k29.CELLS}))["clauses"] == []
    # G BEHIND in one cell only, C BEHIND in one cell with a four-cell mean of -0.215: nothing strikes
    dec = k29.decide(_per({"plain g1": {**GOOD, "own_f": -0.5, "oth_d": -0.86}}))
    assert dec["clauses"] == [] and dec["c_four_cell_mean"] == pytest.approx(-0.215)
    # -0.29 is not BEHIND
    assert k29.decide(_per({c: {**GOOD, "own_f": -0.29, "oth_d": -0.29} for c in k29.CELLS}))["clauses"] == []
    # one AHEAD cell in each world is enough for S1
    assert k29.decide(_per({c: {**GOOD, "own_d": 0.0} for c in ("plain g2", "dm g1")}))["clauses"] == []


@pytest.mark.parametrize("gated,struck", [
    ({"dm g1": {"own_d": [0, 1]}}, True),  # 2 > 0 + 1
    ({"plain g2": {"oth_d": [0, 1, 2], "oth_h": [5]}}, True),  # 3 > 1 + 1
    ({"dm g1": {"own_d": [0]}}, False),  # 1 is not more than 0 + 1
    ({"dm g2": {"oth_d": [0, 1], "oth_h": [7]}}, False),  # 2 is not more than 1 + 1
    ({"dm g1": {"own_h": [0, 1, 2]}}, False),  # fewer gated in the _d arm never strikes
    ({"plain g1": {"own_f": [0, 1, 2, 3], "oth_f": [0, 1, 2]}}, False),  # other arms' gates are not S4's
    ({"dm g2": {"own_d": {0: "beta_gate", 1: "no_response"}}}, True),  # every gated state counts
])
def test_s4_counts_gated_seeds_against_the_h_arm_plus_one(gated, struck):
    dec = k29.decide(_per(gated=gated))
    assert dec["strike"]["S4"] is struck and ("S4" in dec["clauses"]) is struck and bool(dec["s4_cells"]) is struck
    assert dec["clauses"] == (["S4"] if struck else [])


def test_several_clauses_are_named_in_order():
    cells = {c: {**GOOD, "own_d": 0.0, "oth_d": -0.5} for c in PLAIN}
    dec = k29.decide(_per(cells, gated={"dm g2": {"own_d": [0, 1, 2]}}))
    assert dec["clauses"] == ["S1", "S3", "S4"]
    assert dec["outcome"] == "Try K29: STRUCK (S1, S3, S4): no file; C2 stays."


FLIPS = [2.0] * 5 + [0.2] * 5  # Gd +2.0 on seeds 0 to 4 (a gate flip there), +0.2 on the flip-free seeds 5 to 9
FLIPPED = {"own_h": [0, 1, 2, 3, 4]}  # own_h gated on seeds 0 to 4, own_d scored: five flips


def test_l1_c_ahead_in_at_least_as_many_cells_as_gd():
    dec = k29.decide(_per({c: {**GOOD, "oth_d": 1.0} for c in k29.CELLS}))  # C 4, Gd 4
    assert dec["label"]["L1"] and dec["label"]["statement"] == k29.FIT_LABEL and dec["clauses"] == []
    assert dec["outcome"] == k29.NOT_STRUCK_LINE  # the label decides nothing
    dec = k29.decide(_per({**{c: {**GOOD, "oth_d": 1.0} for c in ("plain g1", "dm g1")},
                           **{c: {**GOOD, "own_d": 0.0} for c in ("plain g2", "dm g2")}}))  # C 2, Gd 2
    assert dec["label"]["L1"] and dec["label"]["ahead_cells"] == {"Gd": 2, "C": 2}
    assert not k29.decide(_per({"dm g1": {**GOOD, "oth_d": 1.0}}))["label"]["L1"]  # C 1, Gd 4
    assert k29.decide(_per({c: {**GOOD, "own_d": 0.0} for c in k29.CELLS}))["label"]["L1"]  # C 0, Gd 0, as worded


def test_l2_the_gains_sit_in_gate_flips_of_own_d_and_own_h():
    cells = {"dm g1": {**GOOD, "own_d": FLIPS}, "dm g2": {**GOOD, "own_d": 0.0}}  # Gd AHEAD in dm g1 alone
    dec = k29.decide(_per(cells, gated={"dm g1": FLIPPED}))
    lab = dec["label"]
    assert dec["clauses"] == [] and lab["L2"] and lab["statement"] == k29.FIT_LABEL
    assert lab["flip_seeds"]["dm g1"] == [0, 1, 2, 3, 4] and lab["flip_free_mean_gd"]["dm g1"] == pytest.approx(0.2)
    assert not k29.decide(_per(cells))["label"]["L2"]  # no flip: the flip-free mean is the cell's own
    up = {"dm g1": {**GOOD, "own_d": [2.0] * 5 + [0.5] * 5}}  # a flip-free mean of 0.5
    assert not k29.decide(_per(up, gated={"dm g1": FLIPPED}))["label"]["L2"]
    every = {"dm g1": {**GOOD, "own_d": FLIPS}}  # dm g2 AHEAD without a flip: not every AHEAD dm cell sits in flips
    assert not k29.decide(_per(every, gated={"dm g1": FLIPPED}))["label"]["L2"]
    none_free = {"own_h": {s: "beta_gate" for s in SEEDS}, "own_d": {s: "r2_gate" for s in SEEDS}}
    dec = k29.decide(_per(cells, gated={"dm g1": none_free}))  # gated in both arms, in different states: ten flips
    assert dec["label"]["L2"] and dec["label"]["flip_free_mean_gd"]["dm g1"] is None  # no flip-free seed: 0.3 or less
    assert dec["label"]["flip_seeds"]["dm g1"] == list(SEEDS)


def test_l3_the_gains_hold_in_one_world_only():
    cells = {"dm g1": {**GOOD, "own_d": FLIPS}, "dm g2": {**GOOD, "own_d": 0.0}}
    lab = k29.decide(_per(cells, gated={"dm g1": FLIPPED}))["label"]
    assert lab["L3"] and lab["gd_cells_counted"] == list(PLAIN)  # dm g1 is not counted: its flip-free mean is 0.2
    lab = k29.decide(_per({c: {**GOOD, "own_f": 0.0} for c in PLAIN}))["label"]  # G AHEAD in the dm world only
    assert lab["L3"] and lab["g_ahead_cells"] == list(DM) and not lab["L1"] and not lab["L2"]
    assert not k29.decide(_per())["label"]["L3"]  # Gd and G AHEAD in both worlds
    lab = k29.decide(_per({c: {"own_d": 0.0} for c in k29.CELLS}))["label"]  # no AHEAD cell of Gd or of G
    assert not lab["L3"] and lab["gd_cells_counted"] == [] and lab["g_ahead_cells"] == []


# ------------------------------------------------------------------------------------------------ the reading


def _prov(commit="abc", clean=True, machine="m1") -> dict:
    return {"commit": commit, "tree_clean": clean, "machine": machine, "utc": "x"}


def _body(world: str, task: str, row: dict, part: str, arms, prov: dict) -> dict:
    if world == "plain":
        return {"task": {"key": k29.TASKS[task].key}, "rows": [row], "seed": row["seed"], "world": "plain",
                "part": part, "arms": list(arms), "compositions": _comps(), "provenance": prov}
    return {**row, "target": k18.DM_TARGETS[task], "task": task, "world": "dm", "part": part, "arms": list(arms),
            "compositions": _comps(), "provenance": prov}


def _write_run(root: Path, cells=None, gcells=None, commit="abc", clean=True, machine="m1", gate_passes=True,
               guard_passes=True, score=True, bad=None, gated=None) -> SimpleNamespace:
    """A run folder by hand (gate.json, the guard part, guard.json, run_started.json and the scoring part) and K18's
    committed files by hand beside it, their pins computed here. `bad` maps a file's path to a function that spoils its
    body (None: the file is not written)."""
    cells = cells or {}
    root.mkdir(parents=True)
    prov = _prov(commit, clean, machine)
    gate = {"passes": gate_passes, "failed": [] if gate_passes else [k29.GATE_CHECKS["a"], k29.GATE_CHECKS["b"]],
            "scenario_found": 5 if gate_passes else 4, "own_found": 7, "own_listed": 9, "ratios": {},
            "depletion_factor": k29.DEPLETION_D, "provenance": prov}
    k18.write_task_file(root / "gate.json", gate)
    gsha = {}
    for part, prefix, arms in (("guard", "guard/", ("oth_f",)), ("score", "", k29.ARMS)):
        if part == "score" and not score:
            continue
        src = (cells if gcells is None else gcells) if part == "guard" else cells
        for c in k29.CELLS:
            world, task = c.split()
            t = _cell_table(c, src.get(c, GOOD))
            want = ((harness.FLOOR,) if world == "plain" else ()) + tuple(arms)
            for r in _row_set({a: t[a] for a in want}, (gated or {}).get(c) if part == "score" else None):
                rel = prefix + (f"plain_{task}/seed_{r['seed']}.json" if world == "plain" else
                                f"dm/{task}_seed_{r['seed']}.json")
                body = _body(world, task, r, part, arms, prov)
                if bad and rel in bad:
                    body = bad[rel](body)
                if body is None:
                    continue
                p = root / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                k18.write_task_file(p, body)
                if part == "guard":
                    gsha[rel] = _sha(p)
        if part == "guard":
            k18.write_task_file(root / "guard.json", {"passes": guard_passes,
                                                      "failed": [] if guard_passes else ["oth_f"], "cells": {},
                                                      "files_sha256": gsha, "provenance": prov})
        else:
            (root / "run_started.json").write_text(json.dumps({"row": "K29", "commit": commit, "utc": "x"}))
    refs = root.parent / f"k18_{root.name}"
    refs.mkdir()
    rows = _refs()
    for n in k29.TASKS:
        (refs / f"plain_{n}.json").write_text(json.dumps({"rows": rows[("plain", n)]}))
    (refs / "dm.json").write_text(json.dumps({"rows": {n: rows[("dm", n)] for n in k29.TASKS}}))
    return SimpleNamespace(root=root, refs=refs, pins={f: _sha(refs / f) for f in k29.K18_SHA256}, k18rows=rows)


def _read(w, monkeypatch) -> dict:
    monkeypatch.setattr(k29, "K18_SHA256", w.pins)
    return k29.reading(w.root, w.refs)


def _drop(arm: str):
    def spoil(body):
        row = body["rows"][0] if "rows" in body else body
        for k in ("scores", "skills", "t3"):
            row[k].pop(arm)
        return body

    return spoil


def _set_prov(**kw):
    def spoil(body):
        body["provenance"] = {**body["provenance"], **kw}
        return body

    return spoil


def test_the_reading_of_a_run_folder_that_nothing_strikes(tmp_path, monkeypatch):
    w = _write_run(tmp_path / "ok")
    res = _read(w, monkeypatch)
    assert not res["void"], res.get("why")
    assert res["outcome"] == k29.NOT_STRUCK_LINE and res["decision"]["clauses"] == []
    assert res["tries_counted"] == 5 and res["counted_arms"] == list(k29.COUNTED) and res["run_step_started"]
    assert res["guard"]["passes"] and res["guard_match"]["passes"] and res["gate"]["passes"]
    assert res["commit"] == "abc" and res["machine"] == "m1" and res["run_started"]["row"] == "K29"
    assert len(res["files_sha256"]) == 2 + 40 + 40 and res["compositions"]["dm"] == _comps()
    txt = k29.text(res)
    assert txt.splitlines()[0] == k29.NOT_STRUCK_LINE
    for x in ("G = own_f - own_h", "Gd = own_d - own_h", "C = oth_d - oth_h", "own_d - own_f", "own_h - oth_f",
              "oth_f - floor", "per-term skill, Gd", "realised cut-state shares, dm world", "gate flips of own_d",
              "tries counted: 5 (own_h, own_f, own_d, oth_h, oth_d)", "strike: S1 no, S2 no, S3 no, S4 no",
              "labelled, deciding nothing: no label; L1 False", "target-free gate: passes True"):
        assert x in txt, x


def test_a_struck_reading(tmp_path, monkeypatch):
    res = _read(_write_run(tmp_path / "s", cells={c: {**GOOD, "own_d": 0.0} for c in PLAIN}), monkeypatch)
    line = "Try K29: STRUCK (S1): no file; C2 stays."
    assert not res["void"] and res["outcome"] == line and res["tries_counted"] == 5
    assert k29.text(res).splitlines()[0] == line


def test_a_failed_gate_is_void_naming_the_check_and_counts_nothing(tmp_path, monkeypatch):
    res = _read(_write_run(tmp_path / "g", gate_passes=False, score=False), monkeypatch)
    want = (f"Try K29: VOID: the target-free gate failed ({k29.GATE_CHECKS['a']}; {k29.GATE_CHECKS['b']}); no file.")
    assert res["void"] and res["outcome"] == want and res["tries_counted"] == 0 and res["counted_arms"] == []
    assert not res["run_step_started"] and k29.text(res).splitlines()[0] == want


@pytest.mark.parametrize("kw,why", [
    ({"bad": {"plain_g1/seed_3.json": _drop("own_d")}}, "plain_g1/seed_3.json: missing arms ['own_d']"),
    ({"bad": {"dm/g2_seed_7.json": _drop("oth_h")}}, "dm/g2_seed_7.json: missing arms ['oth_h']"),
    ({"bad": {"plain_g2/seed_0.json": _drop("floor")}}, "plain_g2/seed_0.json: missing arms ['floor']"),
    ({"bad": {"guard/dm/g1_seed_5.json": _drop("oth_f")}}, "guard/dm/g1_seed_5.json: missing arms ['oth_f']"),
    ({"bad": {"dm/g1_seed_4.json": _set_prov(commit="def")}}, "results from different commits"),
    ({"bad": {"plain_g1/seed_9.json": _set_prov(machine="m2")}}, "results from different machines"),
    ({"bad": {"dm/g2_seed_1.json": _set_prov(tree_clean=False)}}, "dm/g2_seed_1.json: a dirty tree"),
    ({"clean": False}, "gate.json: a dirty tree"),
    ({"gcells": {"dm g2": {**GOOD, "oth_f": 0.2}}}, "the guard step's oth_f, read again, is off"),
    ({"cells": {"plain g2": {**GOOD, "oth_f": 0.11}}, "gcells": {}}, "the run step's oth_f is off the guard step's in "
                                                                     "plain g2"),
    ({"cells": {"dm g1": {**GOOD, "oth_f": [0.0] * 3 + [0.6] + [0.0] * 6}}, "gcells": {}},
     "the run step's oth_f is off the guard step's in dm g1"),
])
def test_each_void_condition(tmp_path, monkeypatch, kw, why):
    res = _read(_write_run(tmp_path / "v", **kw), monkeypatch)
    assert res["void"] and any(why in x for x in res["why"]), res["why"]
    assert res["outcome"].startswith("Try K29: VOID: ") and res["outcome"].endswith("; no file.")
    assert why in res["outcome"] and "decision" not in res
    assert res["tries_counted"] == 5 and res["counted_arms"] == list(k29.COUNTED)  # a full run counts five, exactly
    assert "VOID (checked first; no strike clause is applied)" in k29.text(res)


def test_the_run_step_is_matched_with_the_guard_step_not_with_k18s_rows(tmp_path, monkeypatch):
    # the guard step 0.08 off K18's rows (inside the guard), the run step 0.15 off them and 0.07 off the guard step's
    w = _write_run(tmp_path / "m", cells={"plain g1": {**GOOD, "oth_f": 0.15}},
                   gcells={"plain g1": {**GOOD, "oth_f": 0.08}})
    res = _read(w, monkeypatch)
    assert not res["void"], res["why"]
    assert res["guard_match"]["cells"]["plain g1"]["distance"]["mean_diff"] == pytest.approx(0.07)
    assert res["tries_counted"] == 5


def test_k18s_committed_files_off_their_pins_void(tmp_path, monkeypatch):
    w = _write_run(tmp_path / "p")
    (w.refs / "dm.json").write_text((w.refs / "dm.json").read_text() + " ")
    res = _read(w, monkeypatch)
    assert res["void"] and "K18's committed dm.json is missing or not at its pinned sha256" in res["outcome"]
    assert res["tries_counted"] == 5 and res["guard"] is None  # the guard is not read again; oth_f did not fail it
    w2 = _write_run(tmp_path / "q")
    (w2.refs / "plain_g1.json").unlink()
    assert any("plain_g1.json is missing" in x for x in _read(w2, monkeypatch)["why"])


@pytest.mark.parametrize("kw", [
    {"cells": {"plain g1": {**GOOD, "oth_f": 0.2}}, "gcells": {}},  # the run step's oth_f off the guard step's
    {"gcells": {"dm g1": {**GOOD, "oth_f": 0.3}}},  # the guard step's, read again, off K18's rows
    {"bad": {"dm/g1_seed_0.json": _drop("own_d")}},  # a malformed file drops a seed of oth_f's match as well
    {"bad": {"guard/plain_g2/seed_4.json": _drop("oth_f")}},
    {"bad": {"dm/g1_seed_0.json": _set_prov(machine="m9")}},
])
def test_a_full_run_counts_five_whatever_its_files_hold(tmp_path, monkeypatch, kw):
    res = _read(_write_run(tmp_path / "a", **kw), monkeypatch)
    assert res["void"] and res["tries_counted"] == 5 and res["counted_arms"] == list(k29.COUNTED)


def test_a_guard_that_failed_in_the_guard_step_is_void_with_one_counted_try(tmp_path, monkeypatch):
    res = _read(_write_run(tmp_path / "g", guard_passes=False, score=False), monkeypatch)
    assert res["void"] and res["tries_counted"] == 1 and res["counted_arms"] == ["oth_f"]
    assert not res["run_step_started"]
    assert res["outcome"] == ("Try K29: VOID: the guard failed in the guard step: oth_f is off K18's committed "
                              "barred_form; no file.")
    assert k29.text(res).splitlines()[0] == res["outcome"]
    res = _read(_write_run(tmp_path / "h", guard_passes=False), monkeypatch)  # scoring files there all the same
    assert res["void"] and res["tries_counted"] == 1 and res["counted_arms"] == ["oth_f"]


@pytest.mark.parametrize("kw", [{"clean": False}, {"bad": {"dm/g1_seed_2.json": _drop("own_h")}},
                                {"bad": {"plain_g1/seed_1.json": _set_prov(commit="later")}},
                                {"guard_passes": False}, {"gate_passes": False}])
def test_void_is_checked_first_and_wins_over_a_strike(tmp_path, monkeypatch, kw):
    strike = {c: {**GOOD, "own_d": 0.0} for c in k29.CELLS}  # S1 strikes on these tables
    assert _read(_write_run(tmp_path / "s", cells=strike), monkeypatch)["outcome"].startswith("Try K29: STRUCK (S1")
    res = _read(_write_run(tmp_path / "v", cells=strike, **kw), monkeypatch)
    assert res["void"] and res["outcome"].startswith("Try K29: VOID: ") and "decision" not in res
    assert "STRUCK" not in k29.text(res) and "NOT STRUCK" not in res["outcome"]


def test_missing_files_are_not_a_reading_and_no_seed_file_is_opened(tmp_path, monkeypatch):
    w = _write_run(tmp_path / "i", bad={"plain_g2/seed_7.json": lambda b: None})
    (w.root / "dm" / "g1_seed_0.json").write_text("{ not JSON: the reading must not open this file")
    assert _read(w, monkeypatch) == {"incomplete": True, "void": False, "missing": ["plain_g2/seed_7.json"],
                                     "outcome": k29.INCOMPLETE}
    w2 = _write_run(tmp_path / "j", bad={"guard/dm/g1_seed_2.json": lambda b: None})
    (w2.root / "guard.json").write_text("{ not JSON: not opened while a guard file is missing")
    assert _read(w2, monkeypatch)["missing"] == ["guard/dm/g1_seed_2.json"]
    (w2.root / "guard.json").unlink()
    assert _read(w2, monkeypatch)["missing"] == ["guard.json", "guard/dm/g1_seed_2.json"]
    (w2.root / "gate.json").unlink()
    assert _read(w2, monkeypatch)["missing"] == ["gate.json"]
    assert k29.text(_read(w2, monkeypatch)).startswith("Try K29: INCOMPLETE: not a reading")


def test_past_the_cutoff_a_missing_file_makes_the_void_reading(tmp_path, monkeypatch):
    monkeypatch.setattr(k29, "_now", lambda: k29.CUTOFF)  # at the cutoff itself
    w = _write_run(tmp_path / "late", bad={"plain_g2/seed_7.json": lambda b: None})
    res = _read(w, monkeypatch)
    assert res["void"] and res["late"] and res["tries_counted"] == 5 and res["run_step_started"]
    assert res["outcome"] == ("Try K29: VOID: INCOMPLETE past its cutoff (2026-11-25 06:00 UTC): 1 file(s) missing, "
                              "e.g. plain_g2/seed_7.json; no file.")
    w2 = _write_run(tmp_path / "early", score=False)  # the guard step done, the run step never started
    (w2.root / "guard.json").unlink()
    res = _read(w2, monkeypatch)
    assert res["void"] and res["tries_counted"] == 0 and res["outcome"].startswith("Try K29: VOID: INCOMPLETE past")
    monkeypatch.setattr(k29, "_now", lambda: datetime(2026, 11, 25, 5, 59, 59, tzinfo=UTC))
    assert _read(w, monkeypatch)["incomplete"]


# ------------------------------------------------------------------------------------------------ the steps


@pytest.fixture
def ready(tmp_path, monkeypatch):
    """The row in a hand-made ledger, a clean stamp at commit abc on m1, a hand-made run folder and K18's files."""
    tries = tmp_path / "TRIES.md"
    tries.write_text("| K28 | T3 | another row |\n" + ROW)
    monkeypatch.setattr(k29, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: _prov())
    w = _write_run(tmp_path / "run")
    monkeypatch.setattr(k29, "K18_SHA256", w.pins)
    monkeypatch.setattr(k29, "reading", lambda runs, f=k29.reading: f(runs, w.refs))
    return w


def test_read_writes_once_prints_the_outcome_table_and_labels_and_never_overwrites(tmp_path, ready, capsys):
    out = tmp_path / "reading"
    argv = ["read", "--commit", "abc", "--runs", str(ready.root), "--out", str(out)]
    k29.main(argv)
    body = json.loads((out / "reading.json").read_text())
    assert body["outcome"] == k29.NOT_STRUCK_LINE and body["row"] == "K29" and body["tries_counted"] == 5
    assert body["reading_code"]["script_sha256"] == _sha(T3 / "k29_cd7_read.py")
    assert body["reading_code"]["commit"] == "abc"
    assert (out / "reading.txt").read_text().splitlines()[0] == k29.NOT_STRUCK_LINE
    shown = capsys.readouterr().out.splitlines()
    assert shown[0] == k29.NOT_STRUCK_LINE and shown[-1].startswith("wrote ")
    for s in ("G = own_f - own_h: ", "Gd = own_d - own_h: ", "C = oth_d - oth_h: "):
        assert any(x.startswith(f"  {s}") for x in shown), s
    assert any(x.startswith("labelled, deciding nothing: no label; L1 False") for x in shown)
    first = {p.name: p.read_bytes() for p in out.iterdir()}
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 4 and {p.name: p.read_bytes() for p in out.iterdir()} == first
    assert "never overwritten" in capsys.readouterr().err
    (out / "reading.json").unlink()  # reading.txt alone is still a reading
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 4 and not (out / "reading.json").exists()


def test_an_incomplete_folder_writes_nothing_and_exits_3_until_the_cutoff(tmp_path, ready, monkeypatch, capsys):
    (ready.root / "dm" / "g2_seed_4.json").unlink()
    out = tmp_path / "reading"
    argv = ["read", "--commit", "abc", "--runs", str(ready.root), "--out", str(out)]
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 3 and not out.exists() and "not a reading" in capsys.readouterr().out
    monkeypatch.setattr(k29, "_now", lambda: datetime(2026, 11, 25, 6, 0, 1, tzinfo=UTC))
    k29.main(argv)
    body = json.loads((out / "reading.json").read_text())
    assert body["void"] and body["outcome"].startswith("Try K29: VOID: INCOMPLETE past its cutoff")
    assert capsys.readouterr().out.splitlines()[0] == body["outcome"]


def test_a_crash_in_read_is_masked_and_writes_nothing(tmp_path, ready, monkeypatch, capsys):
    def boom(runs):
        raise KeyError("own_d mean 57.1234 against own_h 55.4321")

    monkeypatch.setattr(k29, "reading", boom)
    out = tmp_path / "reading"
    with pytest.raises(SystemExit) as e:
        k29.main(["read", "--commit", "abc", "--runs", str(ready.root), "--out", str(out)])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "KeyError" in err and "57.1234" not in err and not out.exists()


def _steps(tmp_path) -> list[list[str]]:
    o = str(tmp_path / "o")
    return [["gate", "--commit", "abc", "--out", o],
            ["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--out", o],
            ["run", "--commit", "abc", "--world", "dm", "--out", o],
            ["guard", "--commit", "abc", "--runs", o],
            ["read", "--commit", "abc", "--runs", o, "--out", str(tmp_path / "reading")]]


def test_row_k29_and_every_step_refuses_without_it(tmp_path, monkeypatch, capsys):
    p = tmp_path / "TRIES.md"
    script = "`experiments/t3_decomp/k29_cd7_read.py`"
    for text_, ok in (("no rows\n", False), ("", False),
                      (f"| K28 (one row) | T3 | {script} |\n", False),
                      ("| K29 (one row) | T3 | `experiments/t3_decomp/k28_other.py` |\n", False),
                      (f"| K29a (one row) | T3 | {script} |\n", False),
                      (ROW * 2, False),
                      (ROW + f"| R37 (a ruling) | cites {script} |\n", True)):
        p.write_text(text_)
        if ok:
            assert k29.row_present(p) == "K29"
        else:
            with pytest.raises(SystemExit) as e:
                k29.row_present(p)
            assert e.value.code == 2, text_
    assert "row K29 of TRIES.md names" in capsys.readouterr().err
    with pytest.raises(SystemExit) as e:
        k29.row_present(tmp_path / "missing.md")
    assert e.value.code == 2
    monkeypatch.setattr(k29, "TRIES_MD", tmp_path / "none.md")
    seen = []
    monkeypatch.setattr(k18, "stamp", lambda: seen.append("stamp") or _prov())
    monkeypatch.setattr(k18, "load_barred_for_plain", lambda *a: seen.append("load"))
    monkeypatch.setattr(k18, "barred_spec", lambda *a: seen.append("spec"))
    for argv in _steps(tmp_path):
        with pytest.raises(SystemExit) as e:
            k29.main(argv)
        assert e.value.code == 2, argv
    assert seen == [] and not (tmp_path / "o").exists() and not (tmp_path / "reading").exists()


def _writing_steps(tmp_path) -> list[list[str]]:
    o = str(tmp_path / "o")
    return [["gate", "--commit", "abc", "--out", o],
            ["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--out", o],
            ["run", "--commit", "abc", "--part", "guard", "--world", "dm", "--out", o],
            ["run", "--commit", "abc", "--world", "plain", "--task", "g2", "--seeds", "7", "--out", o],
            ["run", "--commit", "abc", "--world", "dm", "--out", o],
            ["guard", "--commit", "abc", "--runs", o]]


def test_gate_run_and_guard_refuse_at_and_after_the_cutoff_and_write_nothing(tmp_path, monkeypatch, capsys):
    tries = tmp_path / "TRIES.md"
    tries.write_text(ROW)
    monkeypatch.setattr(k29, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: _prov())
    seen = []
    monkeypatch.setattr(k29.k22, "frozen_ok", lambda: seen.append("frozen"))
    monkeypatch.setattr(k18, "barred_spec", lambda *a: seen.append("spec"))
    for now in (k29.CUTOFF, datetime(2026, 11, 26, 0, 0, tzinfo=UTC)):
        monkeypatch.setattr(k29, "_now", lambda now=now: now)
        for argv in _writing_steps(tmp_path):
            with pytest.raises(SystemExit) as e:
                k29.main(argv)
            assert e.value.code == 7 and "past the cutoff (2026-11-25 06:00 UTC)" in capsys.readouterr().err, argv
    assert seen == [] and not (tmp_path / "o").exists()
    monkeypatch.setattr(k29, "_now", lambda: datetime(2026, 11, 25, 5, 59, 59, tzinfo=UTC))  # a second before it
    with pytest.raises(SystemExit) as e:
        k29.main(_writing_steps(tmp_path)[0])
    assert e.value.code != 7 and seen == ["frozen", "spec"]  # the clock let it on to the next checks


def test_a_folder_left_incomplete_at_the_cutoff_stays_so_and_reads_void(tmp_path, ready, monkeypatch):
    for f in ("seed_7.json", "seed_7.json.prov.json"):
        (ready.root / "plain_g2" / f).unlink()
    monkeypatch.setattr(k29, "_now", lambda: k29.CUTOFF)
    with pytest.raises(SystemExit) as e:  # a rerun of the missing seed, at the cutoff
        k29.main(["run", "--commit", "abc", "--world", "plain", "--task", "g2", "--seeds", "7", "--out",
                  str(ready.root)])
    assert e.value.code == 7 and not (ready.root / "plain_g2" / "seed_7.json").exists()
    out = tmp_path / "reading"
    k29.main(["read", "--commit", "abc", "--runs", str(ready.root), "--out", str(out)])
    body = json.loads((out / "reading.json").read_text())
    assert body["outcome"] == ("Try K29: VOID: INCOMPLETE past its cutoff (2026-11-25 06:00 UTC): 1 file(s) missing, "
                               "e.g. plain_g2/seed_7.json; no file.")
    assert body["tries_counted"] == 5


@pytest.mark.parametrize("stamp,why", [({"tree_clean": False}, "the tree is not clean"),
                                       ({"commit": "def"}, "HEAD is def, not the run commit abc")])
def test_every_step_refuses_an_unclean_tree_or_another_head(tmp_path, monkeypatch, capsys, stamp, why):
    tries = tmp_path / "TRIES.md"
    tries.write_text(ROW)
    monkeypatch.setattr(k29, "TRIES_MD", tries)
    monkeypatch.setattr(k18, "stamp", lambda: {**_prov(), **stamp})
    for argv in _steps(tmp_path):
        with pytest.raises(SystemExit) as e:
            k29.main(argv)
        assert e.value.code == 2 and why in capsys.readouterr().err, argv
    assert not (tmp_path / "o").exists() and not (tmp_path / "reading").exists()


def _write_gate(out: Path, passes=True, **prov) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "gate.json").write_text(json.dumps({"passes": passes, "provenance": {**_prov(), **prov}}))


@pytest.fixture
def checker(tmp_path, monkeypatch):
    """cmd_guard over a hand-made guard part: the preconditions passed, K18's rows by hand."""
    w = _write_run(tmp_path / "gc", score=False)
    for f in ("guard.json", "guard.json.prov.json"):
        (w.root / f).unlink()
    monkeypatch.setattr(k29, "preconditions", lambda commit: ("K29", _prov()))
    monkeypatch.setattr(k29, "guard_refs", lambda refs=None: w.k18rows)
    return w


def test_the_guard_check_writes_guard_json_once(checker, capsys):
    argv = ["guard", "--commit", "abc", "--runs", str(checker.root)]
    k29.main(argv)
    g = json.loads((checker.root / "guard.json").read_text())
    assert g["passes"] and g["failed"] == [] and g["row"] == "K29" and set(g["files_sha256"]) == set(k29.guard_files())
    assert all(v["distance"]["max_abs_diff"] == 0.0 for v in g["cells"].values())
    out = capsys.readouterr().out
    assert out == "guard: passes\n"
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 4  # made once


def test_the_guard_check_exits_6_when_the_guard_fails(tmp_path, monkeypatch, capsys):
    w = _write_run(tmp_path / "gf", gcells={"dm g2": {**GOOD, "oth_f": 0.2}}, score=False)
    for f in ("guard.json", "guard.json.prov.json"):
        (w.root / f).unlink()
    monkeypatch.setattr(k29, "preconditions", lambda commit: ("K29", _prov()))
    monkeypatch.setattr(k29, "guard_refs", lambda refs=None: w.k18rows)
    with pytest.raises(SystemExit) as e:
        k29.main(["guard", "--commit", "abc", "--runs", str(w.root)])
    g = json.loads((w.root / "guard.json").read_text())
    assert e.value.code == 6 and g["passes"] is False and g["failed"] == ["oth_f"] and not g["cells"]["dm g2"]["passes"]
    assert "FAILS" in capsys.readouterr().out


@pytest.mark.parametrize("spoil,code,why", [
    ("missing", 3, ""),
    ("commit", 2, "not this commit's and machine's"),
    ("pins", 2, "not at its pinned sha256"),
    ("gate", 2, "the target-free gate did not pass"),
])
def test_the_guard_check_refuses_or_waits_and_writes_nothing(checker, monkeypatch, capsys, spoil, code, why):
    f = checker.root / "guard" / "dm" / "g1_seed_6.json"
    if spoil == "missing":
        f.unlink()
    elif spoil == "commit":
        body = json.loads(f.read_text())
        body["provenance"]["commit"] = "def"
        f.write_text(json.dumps(body))
    elif spoil == "pins":
        def off(refs=None):
            raise k29.Refused("K18's committed dm.json is missing or not at its pinned sha256: refused")

        monkeypatch.setattr(k29, "guard_refs", off)
    else:
        _write_gate(checker.root, passes=False)
    with pytest.raises(SystemExit) as e:
        k29.main(["guard", "--commit", "abc", "--runs", str(checker.root)])
    assert e.value.code == code and why in capsys.readouterr().err
    assert not (checker.root / "guard.json").exists()


@pytest.fixture
def fake_run(monkeypatch, tmp_path):
    """cmd_run with the preconditions passed, a passing gate.json, the barred spec and factor by hand and the harness
    stubbed: the fake harness logs each seed's line in the harness's own format (every arm's score) and leaves the
    fold's compositions where K18's plain build leaves them."""
    calls = {"names": [], "seeds": []}
    lines = {}
    st = _prov()

    def run_task(t, cands, seeds, data_dir, log=print):
        calls["names"] = [c.name for c in cands]
        calls["seeds"] += list(seeds)
        k18._FOLD["F"] = {"k29": {"figures": _comps()}}
        for sd in seeds:
            lines[sd] = f"{t.key} seed {sd}: floor 50.00, oth_f {60 - sd / 9:.2f}, own_d {61 + sd / 7:.2f}"
            log(lines[sd])
        return {"task": {"key": t.key}, "rows": [{"seed": sd} for sd in seeds],
                "provenance": {"commit": "abc", "tree_clean": True}}

    def load(*a):
        k18._BARRED.clear()
        k18._BARRED.update(phi=np.ones(3), genes=["a", "b", "c"], spec=SPEC, record={})
        return {"file": "f"}

    monkeypatch.setattr(k29, "preconditions", lambda commit: ("K29", st))
    monkeypatch.setattr(k29, "guard_refs", lambda refs=None: {})
    monkeypatch.setattr(k18, "barred_spec", lambda decision=None: (SPEC, DEC, RULES))
    monkeypatch.setattr(k18, "load_barred_for_plain", load)
    monkeypatch.setattr(harness, "run_task", run_task)
    monkeypatch.setattr(k18, "stamp", lambda: st)
    out = tmp_path / "run"
    _write_gate(out)
    yield SimpleNamespace(calls=calls, lines=lines, out=out, st=st)
    k18._BARRED.clear()


def test_the_guard_part_runs_oth_f_only_and_logs_progress_only(fake_run, capsys):
    argv = ["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--seeds", "0-1", "--out",
            str(fake_run.out)]
    k29.main(argv)
    assert fake_run.calls["names"] == [harness.FLOOR, "oth_f"] and fake_run.calls["seeds"] == [0, 1]
    printed = capsys.readouterr()
    for x in printed.out.splitlines():
        assert re.fullmatch(r"K29 guard plain g1 seed \d: scored|wrote .+/guard/plain_g1/seed_\d\.json: commit abc, "
                            r"clean True, machine m1", x), x
    for line in fake_run.lines.values():
        for score in re.findall(r"\d+\.\d\d", line):
            assert score not in printed.out and score not in printed.err
    body = json.loads((fake_run.out / "guard" / "plain_g1" / "seed_0.json").read_text())
    assert body["part"] == "guard" and body["arms"] == ["oth_f"] and body["compositions"] == _comps()
    assert body["provenance"]["machine"] == "m1" and k29._K29["lists"] == LISTS
    assert k18.ARMS == k29._K18_ORIGINAL["ARMS"] and not (fake_run.out / "run_started.json").exists()
    k29.main(argv)  # a rerun skips the seeds written
    assert fake_run.calls["seeds"] == [0, 1] and "present, skipped" in capsys.readouterr().out


def test_the_guard_part_refuses_without_a_passing_gate_of_this_commit_and_machine(fake_run, capsys):
    argv = ["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--seeds", "0", "--out",
            str(fake_run.out)]
    for kw, why in (({"passes": False}, "the target-free gate did not pass"),
                    ({"commit": "def"}, "gate.json is from another commit or machine"),
                    ({"machine": "m2"}, "gate.json is from another commit or machine")):
        _write_gate(fake_run.out, **kw)
        with pytest.raises(SystemExit) as e:
            k29.main(argv)
        assert e.value.code == 2 and why in capsys.readouterr().err
    (fake_run.out / "gate.json").unlink()
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 2 and "the target-free gate runs and passes first" in capsys.readouterr().err
    assert fake_run.calls["seeds"] == []


def test_a_scoring_part_refuses_without_a_passing_guard_and_marks_the_run_steps_start(fake_run, capsys):
    out = fake_run.out
    argv = ["run", "--commit", "abc", "--world", "plain", "--task", "g2", "--seeds", "0", "--out", str(out)]
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 2 and "the guard part runs and passes first" in capsys.readouterr().err
    for body, why in (({"passes": False, "provenance": fake_run.st}, "the guard part did not pass"),
                      ({"passes": True, "provenance": {**fake_run.st, "commit": "def"}}, "another commit or machine"),
                      ({"passes": True, "provenance": {**fake_run.st, "machine": "m2"}}, "another commit or machine"),
                      ("[]", "did not pass"), ("not json", "is not JSON")):
        (out / "guard.json").write_text(body if isinstance(body, str) else json.dumps(body))
        with pytest.raises(SystemExit) as e:
            k29.main(argv)
        assert e.value.code == 2 and why in capsys.readouterr().err
    assert fake_run.calls["seeds"] == [] and not (out / "run_started.json").exists()  # no count while refused
    (out / "guard.json").write_text(json.dumps({"passes": True, "provenance": fake_run.st}))
    k29.main(argv)
    assert fake_run.calls["names"] == [harness.FLOOR, *k29.ARMS] and fake_run.calls["seeds"] == [0]
    body = json.loads((out / "plain_g2" / "seed_0.json").read_text())
    assert body["part"] == "score" and body["arms"] == list(k29.ARMS) and body["task_source"] == k29.SCRIPT
    mark = json.loads((out / "run_started.json").read_text())
    assert mark["counted"] == list(k29.COUNTED) and mark["commit"] == "abc"
    k29.main([*argv[:-4], "--seeds", "1", "--out", str(out)])  # a second part finds the mark there, unchanged
    assert json.loads((out / "run_started.json").read_text()) == mark


def test_a_run_refuses_once_a_reading_is_there_whatever_the_folder(fake_run, monkeypatch, capsys, tmp_path):
    results = tmp_path / "results"
    results.mkdir()
    monkeypatch.setattr(k29, "RESULTS", results)
    (fake_run.out / "guard.json").write_text(json.dumps({"passes": True, "provenance": fake_run.st}))
    argvs = [["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--seeds", "0", "--out",
              str(fake_run.out)],
             ["run", "--commit", "abc", "--world", "plain", "--task", "g2", "--seeds", "0", "--out", str(fake_run.out)],
             ["run", "--commit", "abc", "--world", "dm", "--out", str(tmp_path / "another_run")]]
    for name in ("reading.json", "reading.txt"):
        (results / name).write_text("{}")
        for argv in argvs:
            with pytest.raises(SystemExit) as e:
                k29.main(argv)
            assert e.value.code == 4 and "no arm is scored again" in capsys.readouterr().err, argv
        (results / name).unlink()
    assert fake_run.calls["seeds"] == [] and not (fake_run.out / "run_started.json").exists()
    assert not (fake_run.out / "guard").exists() and not (tmp_path / "another_run").exists()
    k29.main(argvs[1])  # no reading there: the part runs
    assert fake_run.calls["seeds"] == [0]


def test_a_run_refuses_a_decision_file_off_the_row_and_k18s_files_off_their_pins(fake_run, monkeypatch, capsys):
    argv = ["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--seeds", "0", "--out",
            str(fake_run.out)]
    comp = {**SPEC["composition"], "log2_magnitude": 0.75}
    monkeypatch.setattr(k18, "barred_spec", lambda decision=None: ({**SPEC, "composition": comp}, DEC, RULES))
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 2 and "not the row's second" in capsys.readouterr().err

    def off(refs=None):
        raise k29.Refused("K18's committed dm.json is missing or not at its pinned sha256: refused")

    monkeypatch.setattr(k29, "guard_refs", off)
    with pytest.raises(SystemExit) as e:
        k29.main(argv)
    assert e.value.code == 2 and "pinned sha256" in capsys.readouterr().err
    assert fake_run.calls["seeds"] == [] and not (fake_run.out / "guard").exists()


def test_a_crash_is_masked_writes_no_task_file_and_restores_k18(fake_run, monkeypatch, capsys):
    def boom(*a, **k):
        raise RuntimeError("oth_f 61.234 against barred_form 61.2341")

    monkeypatch.setattr(harness, "run_task", boom)
    with pytest.raises(SystemExit) as e:
        k29.main(["run", "--commit", "abc", "--part", "guard", "--world", "plain", "--task", "g1", "--seeds", "0",
                  "--out", str(fake_run.out)])
    err = capsys.readouterr().err
    assert e.value.code == 1 and "61.234" not in err and "<num>" in err and "RuntimeError" in err
    assert (k18.ARMS, k18.build_arm, k18.fold_factors) == (k29._K18_ORIGINAL["ARMS"], k29._K18_ORIGINAL["build_arm"],
                                                           k29._K18_ORIGINAL["fold_factors"])
    assert not (fake_run.out / "guard" / "plain_g1" / "seed_0.json").exists()


# ------------------------------------------------------------------------------------------------ the launcher

SH = T3 / "run_k29.sh"
OUT_REL = "experiments/t3_decomp/results/k29_cd7_read"
SCRIPT_REL = "experiments/t3_decomp/k29_cd7_read.py"
GIT_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_OBJECT_DIRECTORY", "GIT_COMMON_DIR", "GIT_NAMESPACE",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES")


def test_the_launcher_text():
    sh = SH.read_text()
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    assert ":?" not in sh and not re.search(r"\$\{[^}]*:\?[^}]*'[^}]*\}", sh)  # no ${var:?message}, so no apostrophe
    assert not re.search(r"(^|[\s;|&(])rm(\s|$)", sh, flags=re.M)  # nothing deleted
    threads = [int(x) for x in re.findall(r"^\s*threads (\d+)", sh, flags=re.M)]
    assert threads and max(threads) <= 16 and "--workers 8 --threads 2" in sh and "--threads 8" in sh
    assert f"SCRIPT={SCRIPT_REL}" in sh and f"OUT={OUT_REL}" in sh
    assert 'S="$HOME/v/scratch/k29/$RUN"' in sh and 'R="runs/k29_$RUN"' in sh
    assert "RUN_DIR=$HOME/v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z" in sh
    for step, prev in (("tests", "preflight"), ("gate", "tests"), ("guard", "gate"), ("run", "guard"),
                       ("read", "run")):
        assert f"begin {step} {prev}" in sh
    order = [sh.index(f"begin {s} ") for s in ("tests", "gate", "guard", "run", "read")]
    assert order == sorted(order)
    assert "tests/test_k29_cd7_read.py" in sh and "tests/test_k18_built_files.py" in sh
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh and "set -u" in sh
    # past the cutoff: a step's exit 7 makes the reading once and finishes 27; a reading there: no part starts (44)
    assert sh.count("if [ $rc = 7 ]; then void_read 27; fi") == 2 and "if late guard_; then void_read 27; fi" in sh
    assert 'if late ""; then void_read 27; fi' in sh
    guard_step, run_step = sh.split("begin guard gate", 1)[1], sh.split("begin run guard", 1)[1]
    assert guard_step.index("if reading_there; then") < guard_step.index("parts guard_")
    assert run_step.index("if reading_there; then") < run_step.index('parts ""')


def _git_env() -> dict:
    env = {k: v for k, v in os.environ.items() if k not in GIT_VARS}
    return {**env, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1", "GIT_OPTIONAL_LOCKS": "0"}


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                        "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, env=_git_env(),
                       text=True)
    return r.stdout.strip()


UV_STUB = """#!{py} -I
# uv run --frozen --no-sync python -I -B ARGS: the tests exit $HOME/tests_rc; a -c one-liner runs in python -I -B; every
# call of the script is logged in $HOME/script_calls; `gate` and `guard` exit $HOME/gate_rc or guard_rc when it is
# there (writing nothing), else write gate.json or guard.json (failing, exit 5 or 6, when $HOME/gate_fails or
# guard_fails is there); `run` exits $HOME/run_rc and `read` $HOME/read_rc (0 when absent)
import json, os, sys
from pathlib import Path
a = sys.argv[1:]
if a[:6] != ["run", "--frozen", "--no-sync", "python", "-I", "-B"]:
    sys.exit(97)
a = a[6:]
home = Path(os.environ["HOME"])
def code(name):
    p = home / name
    return int(p.read_text()) if p.exists() else 0
if a[0] == "-m":
    sys.exit(code("tests_rc"))
if a[0] == "-c":
    os.execv(sys.executable, [sys.executable, "-I", "-B", *a])
with open(home / "script_calls", "a") as fh:
    fh.write(" ".join(a) + "\\n")
opt = dict(zip(a[2::2], a[3::2]))
if a[1] in ("gate", "guard") and code(a[1] + "_rc"):
    sys.exit(code(a[1] + "_rc"))
if a[1] in ("gate", "guard"):
    ok = not (home / f"{a[1]}_fails").exists()
    where = opt["--out"] if a[1] == "gate" else opt["--runs"]
    Path(where).mkdir(parents=True, exist_ok=True)
    Path(where, f"{a[1]}.json").write_text(json.dumps({"passes": ok, "provenance": {"commit": opt["--commit"]}}))
    sys.exit(0 if ok else (5 if a[1] == "gate" else 6))
sys.exit(code(a[1] + "_rc"))
"""


def _launcher_repo(tmp_path: Path) -> SimpleNamespace:
    repo, home, bindir = tmp_path / "repo", tmp_path / "home", tmp_path / "bin"
    (home / "v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z").mkdir(parents=True)
    files = {"experiments/t3_decomp/run_k29.sh": SH.read_text(), "experiments/harness/TRIES.md": ROW,
             ".gitignore": "runs/\n",
             **{f"experiments/t3_decomp/results/k18_built_files/tasks/{f}": "{}\n"
                for f in ("plain_g1.json", "plain_g2.json", "dm.json")}}
    repo.mkdir()
    _git(repo, "init", "-q")
    for rel, t in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(t)
    _git(repo, "add", ".")
    _git(repo, "commit", "-q", "-m", "init")
    bindir.mkdir()
    (bindir / "uv").write_text(UV_STUB.replace("{py}", sys.executable))
    (bindir / "uv").chmod(0o755)
    return SimpleNamespace(repo=repo, home=home, bin=bindir, head=_git(repo, "rev-parse", "HEAD"))


def _launch(L, *args) -> tuple[int, str]:
    env = {**_git_env(), "HOME": str(L.home), "PATH": f"{L.bin}:{os.environ['PATH']}",
           "GIT_CEILING_DIRECTORIES": str(L.repo.parent)}
    r = subprocess.run(["bash", "experiments/t3_decomp/run_k29.sh", *args], cwd=L.repo, env=env, capture_output=True,
                       text=True, timeout=120)
    return r.returncode, r.stderr


def _exit(L, run: str, step: str) -> str:
    p = L.home / "v" / "scratch" / "k29" / run / f"{step}.exit"
    return p.read_text().strip() if p.exists() else "absent"


def _calls(L) -> list[str]:
    p = L.home / "script_calls"
    return p.read_text().splitlines() if p.exists() else []


def _kinds(L) -> list[str]:
    return [c.split()[1] for c in _calls(L)]


def test_the_launcher_refusals_in_a_throwaway_repository(tmp_path):
    L = _launcher_repo(tmp_path)
    tries = L.repo / "experiments/harness/TRIES.md"
    code, err = _launch(L)
    assert code == 2 and "usage: run_k29.sh COMMIT RUN" in err
    assert _launch(L, L.head)[0] == 2
    for bad in (".hidden", "a/b", "x y"):
        assert _launch(L, L.head, bad)[0] == 90
    code, err = _launch(L, "0" * 40, "r1")
    assert code == 97 and "HEAD is not" in err
    assert _exit(L, "r1", "done") == "97" and _exit(L, "r1", "preflight") == "97"
    code, err = _launch(L, L.head, "r1")  # r1 was made for another commit
    assert code == 91 and "another commit" in err
    (L.home / "tests_rc").write_text("1")
    assert _launch(L, L.head, "r2")[0] == 10 and _exit(L, "r2", "tests") == "10" and _exit(L, "r2", "preflight") == "0"
    (L.home / "tests_rc").write_text("0")
    tries.write_text(ROW + "an edit\n")
    code, err = _launch(L, L.head, "r3")
    assert code == 98 and "not clean" in err
    _git(L.repo, "checkout", "--", "experiments/harness/TRIES.md")
    tries.write_text(ROW.replace("| K29 ", "| K28 "))
    _git(L.repo, "commit", "-q", "-am", "no row K29")
    code, err = _launch(L, _git(L.repo, "rev-parse", "HEAD"), "r4")
    assert code == 95 and "no TRIES.md row K29" in err
    tries.write_text(ROW)
    _git(L.repo, "rm", "-q", "--cached", "experiments/t3_decomp/results/k18_built_files/tasks/dm.json")
    _git(L.repo, "commit", "-q", "-am", "untrack K18's dm.json")
    code, err = _launch(L, _git(L.repo, "rev-parse", "HEAD"), "r5")
    assert code == 96 and "dm.json is not tracked" in err
    _git(L.repo, "add", "experiments/t3_decomp/results/k18_built_files/tasks/dm.json")
    _git(L.repo, "commit", "-q", "-m", "track it again")
    head = _git(L.repo, "rev-parse", "HEAD")
    run_dir = L.home / "v/virtual-embryo-challenge/runs/t3_decomp/gata4_20261028T122425Z"
    run_dir.rename(run_dir.with_name("elsewhere"))
    code, err = _launch(L, head, "r6")
    assert code == 96 and "gata4 run folder" in err
    run_dir.with_name("elsewhere").rename(run_dir)
    assert _calls(L) == []  # no refusal called the script


def test_the_launcher_runs_tests_gate_guard_run_and_read_in_order_with_exit_records(tmp_path):
    L = _launcher_repo(tmp_path)
    code, err = _launch(L, L.head, "k29a")
    assert code == 0, err
    for step in ("preflight", "tests", "gate", "guard", "run", "read", "done"):
        assert _exit(L, "k29a", step) == "0", step
        if step != "done":
            assert (L.home / "v/scratch/k29/k29a" / f"{step}.head").read_text().strip() == L.head
    calls = _calls(L)
    assert _kinds(L) == ["gate"] + ["run"] * 5 + ["guard"] + ["run"] * 5 + ["read"]
    assert all(c.startswith(f"{SCRIPT_REL} ") and f"--commit {L.head}" in c for c in calls)
    assert calls[0] == f"{SCRIPT_REL} gate --commit {L.head} --out runs/k29_k29a --threads 8"
    assert all("--part guard" in c for c in calls[1:6]) and all("--part score" in c for c in calls[7:12])
    assert sorted(c.split("--world ")[1].split()[0] for c in calls[1:6]) == ["dm", "plain", "plain", "plain", "plain"]
    assert sum("--world dm --workers 8 --threads 2" in c for c in calls) == 2
    assert calls[6] == f"{SCRIPT_REL} guard --commit {L.head} --runs runs/k29_k29a"
    assert calls[12] == f"{SCRIPT_REL} read --commit {L.head} --runs runs/k29_k29a --out {OUT_REL}"
    (L.home / "script_calls").unlink()
    code, err = _launch(L, L.head, "k29a")  # the same RUN at the same commit: gate.json and guard.json read again
    assert code == 0, err
    assert _kinds(L) == ["run"] * 5 + ["read"]


def test_a_failed_gate_makes_the_void_reading_once_and_scores_nothing(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.home / "gate_fails").write_text("")
    code, err = _launch(L, L.head, "g")
    assert code == 25 and _exit(L, "g", "done") == "25" and _exit(L, "g", "gate") == "25", err
    assert _exit(L, "g", "guard") == "absent" and _exit(L, "g", "void_read") == "0"
    assert _kinds(L) == ["gate", "read"]
    (L.repo / OUT_REL).mkdir(parents=True)
    (L.repo / OUT_REL / "reading.json").write_text("{}")
    (L.home / "script_calls").unlink()
    assert _launch(L, L.head, "g")[0] == 25 and _calls(L) == []  # the reading is there: 25 again, nothing called
    R = L.repo / "runs" / "k29_g"
    (R / "gate.json").write_text(json.dumps({"passes": True, "provenance": {"commit": "f" * 40}}))
    assert _launch(L, L.head, "g")[0] == 21 and _calls(L) == []
    (R / "gate.json").write_text("not json")
    assert _launch(L, L.head, "g")[0] == 20 and _calls(L) == []


@pytest.mark.parametrize("rc,want", [(1, 23), (2, 22), (4, 20)])
def test_a_crash_or_refusal_of_the_gate_ends_the_call_without_a_reading(tmp_path, rc, want):
    L = _launcher_repo(tmp_path)
    (L.home / "gate_rc").write_text(str(rc))
    assert _launch(L, L.head, "c")[0] == want and _exit(L, "c", "gate") == str(want)
    assert _kinds(L) == ["gate"] and not (L.repo / "runs" / "k29_c" / "gate.json").exists()


def test_a_failed_guard_makes_the_void_reading_once_and_starts_no_scoring_part(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.home / "guard_fails").write_text("")
    code, err = _launch(L, L.head, "g")
    assert code == 26 and _exit(L, "g", "done") == "26" and _exit(L, "g", "guard") == "26", err
    assert _exit(L, "g", "run") == "absent" and _exit(L, "g", "read") == "absent"
    assert _kinds(L) == ["gate"] + ["run"] * 5 + ["guard", "read"]
    assert _calls(L)[-1] == f"{SCRIPT_REL} read --commit {L.head} --runs runs/k29_g --out {OUT_REL}"
    (L.repo / OUT_REL).mkdir(parents=True)
    (L.repo / OUT_REL / "reading.json").write_text("{}")
    (L.home / "script_calls").unlink()
    assert _launch(L, L.head, "g")[0] == 26 and _calls(L) == []  # the reading is there: 26 again, nothing called
    R = L.repo / "runs" / "k29_g"
    (R / "guard.json").write_text(json.dumps({"passes": True, "provenance": {"commit": "f" * 40}}))
    assert _launch(L, L.head, "g")[0] == 21 and _calls(L) == []
    (R / "guard.json").write_text("not json")
    assert _launch(L, L.head, "g")[0] == 20 and _calls(L) == []


def test_a_failed_part_stops_the_call_before_the_next_step(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.home / "run_rc").write_text("1")
    assert _launch(L, L.head, "a")[0] == 30 and _exit(L, "a", "guard") == "30"
    assert _kinds(L) == ["gate"] + ["run"] * 5  # the guard parts only: no check, no scoring part
    (L.home / "script_calls").unlink()
    R = L.repo / "runs" / "k29_b"
    R.mkdir(parents=True)
    for f in ("gate.json", "guard.json"):
        (R / f).write_text(json.dumps({"passes": True, "provenance": {"commit": L.head}}))
    assert _launch(L, L.head, "b")[0] == 30 and _exit(L, "b", "run") == "30" and _exit(L, "b", "read") == "absent"
    assert _kinds(L) == ["run"] * 5


@pytest.mark.parametrize("rc,want", [(3, 43), (4, 44), (2, 42), (1, 41)])
def test_the_read_steps_codes(tmp_path, rc, want):
    L = _launcher_repo(tmp_path)
    (L.home / "read_rc").write_text(str(rc))
    assert _launch(L, L.head, "r")[0] == want and _exit(L, "r", "read") == str(want) and _exit(L, "r", "run") == "0"


def _passing(L, run: str) -> None:
    R = L.repo / "runs" / f"k29_{run}"
    R.mkdir(parents=True)
    for f in ("gate.json", "guard.json"):
        (R / f).write_text(json.dumps({"passes": True, "provenance": {"commit": L.head}}))


@pytest.mark.parametrize("stub,passing,kinds,step", [
    ({"gate_rc": "7"}, False, ["gate", "read"], "gate"),
    ({"run_rc": "7"}, False, ["gate"] + ["run"] * 5 + ["read"], "guard"),  # the guard parts
    ({"guard_rc": "7"}, False, ["gate"] + ["run"] * 5 + ["guard", "read"], "guard"),  # the guard check
    ({"run_rc": "7"}, True, ["run"] * 5 + ["read"], "run"),  # the scoring parts
])
def test_a_refusal_past_the_cutoff_makes_the_reading_once_and_finishes_27(tmp_path, stub, passing, kinds, step):
    L = _launcher_repo(tmp_path)
    for name, value in stub.items():
        (L.home / name).write_text(value)
    if passing:
        _passing(L, "late")
    code, err = _launch(L, L.head, "late")
    assert code == 27 and _exit(L, "late", step) == "27" and _exit(L, "late", "void_read") == "0", err
    assert _kinds(L) == kinds
    assert _calls(L)[-1] == f"{SCRIPT_REL} read --commit {L.head} --runs runs/k29_late --out {OUT_REL}"


def test_a_second_run_after_a_reading_scores_no_arm(tmp_path):
    L = _launcher_repo(tmp_path)
    (L.repo / OUT_REL).mkdir(parents=True)
    (L.repo / OUT_REL / "reading.json").write_text("{}")
    code, err = _launch(L, L.head, "second")  # a fresh RUN: the gate (target-free) runs, no guard part starts
    assert code == 44 and _exit(L, "second", "guard") == "44" and _kinds(L) == ["gate"], err
    (L.home / "script_calls").unlink()
    _passing(L, "third")  # a RUN whose gate and guard passed: no scoring part starts
    code, err = _launch(L, L.head, "third")
    assert code == 44 and _exit(L, "third", "run") == "44" and _calls(L) == [], err
